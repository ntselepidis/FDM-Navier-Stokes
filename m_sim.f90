module m_sim
  use m_bc,     only: apply_boundary_conditions
  use m_fdm
  use m_mg,     only: MGsolve_2DPoisson
  use m_stream, only: compute_velocity
  use m_nvtx,   only: nvtx_push, nvtx_pop
#ifdef NVHPC
  use curandex
#endif
  implicit none

  private
  public :: simulation

  ! ----------------------
  ! simulation type
  ! ----------------------
  type :: simulation
    real :: dt
    real :: total_time
  contains
    procedure, private, nopass :: compute_dt
    procedure, public, pass    :: read_inputs
    procedure, public, nopass  :: initialize
    procedure, public, pass    :: evolve
    procedure, public, nopass  :: export_results
    procedure, public, nopass  :: destroy
  end type simulation

  ! -------------------------
  ! constants
  ! -------------------------
  real, parameter :: pi = 3.141592653589793238462643383279502884 ! pi
  real, parameter :: k = 1.0 ! thermal diffusivity
  integer, parameter :: niters = 50 ! maximum number of iterations of MG solver

  ! -------------------------
  ! ids of input/output files
  ! -------------------------
  integer, parameter :: input_file_id = 1
  integer, parameter :: v_file_id = 2
  integer, parameter :: S_file_id = 1
  integer, parameter :: W_file_id = 1
  integer, parameter :: T_file_id = 1

  ! -------------------------
  ! default parameter file
  ! -------------------------
  character(len=100) :: input_file_name = 'parameters.txt'

  ! ----------------------------------------------
  ! input parameters to be set from parameter file
  ! ----------------------------------------------
  integer :: nx = 257, ny = 65          ! number of grid-points in x- and y-direction
  real :: total_time = 0.1              ! integration time
  real :: Ra = 1.e6                     ! Rayleigh number
  real :: Pr = 0.01                     ! Prandtl number
  real :: beta = 0.5                    ! explicit/implicit/semi-implicit
  real :: err = 1.e-3                   ! prescribed tolerance for MG solver
  real :: a_dif = 0.15, a_adv = 0.4     ! diffusion and advection timestep parameters
  character(len=15) :: Tinit = 'cosine' ! initial condition for temperature

  ! ----------------------------------------------
  ! utility variables
  ! ----------------------------------------------
  real :: h, dt_dif, dt_adv

  ! ----------------------------------------------
  ! field arrays
  ! ----------------------------------------------
  real, allocatable :: S(:, :), vx(:, :), vy(:, :), v(:, :)
  real, allocatable :: T(:, :), dT2(:, :), dTx(:, :), dTy(:, :), T_rhs(:, :)
  real, allocatable :: W(:, :), dW2(:, :), dWx(:, :), dWy(:, :), W_rhs(:, :)
  real, allocatable :: Ra_dTdx(:, :)

contains

  subroutine read_inputs(this)
    implicit none
    ! arguments
    class(simulation), intent(inout) :: this
    ! local variables
    namelist /inputs/ nx, ny, total_time, Ra, Pr, beta, err, a_dif, a_adv, Tinit

    ! if parameter file is passed as command line argument, use this one
    if (command_argument_count() > 0) then
      call get_command_argument(1, input_file_name)
    end if
    print*, 'Reading inputs from ', input_file_name

    ! read inputs from file
    open(input_file_id, file=input_file_name)
    read(input_file_id, inputs)
    close(input_file_id)
    write(*, inputs) ! echo values to stdout

    this%total_time = total_time

  end subroutine

  subroutine initialize()
    implicit none
    ! local variables
    real :: width
    integer :: i, j

    h = 1.0 / ( ny - 1.0 ) ! mesh size
    width = ( nx - 1.0 ) / ( ny - 1.0 ) ! width

    ! compute diffusive timestep
    dt_dif = ( a_dif * min(h, h)**2 ) / max(k, Pr)

    ! allocate field arrays
    allocate(S(nx, ny), vx(nx, ny), vy(nx, ny), v(nx, ny))
    allocate(T(nx, ny), dT2(nx, ny), dTx(nx, ny), dTy(nx, ny), T_rhs(nx, ny))
    allocate(W(nx, ny), dW2(nx, ny), dWx(nx, ny), dWy(nx, ny), W_rhs(nx, ny))
    allocate(Ra_dTdx(nx, ny))

    do concurrent (i=1:nx, j=1:ny)
      S(i, j) = 0.0
      vx(i, j) = 0.0
      vy(i, j) = 0.0
      v(i, j) = 0.0
      T(i, j) = 0.0
      dT2(i, j) = 0.0
      dTx(i, j) = 0.0
      dTy(i, j) = 0.0
      T_rhs(i, j) = 0.0
      dW2(i, j) = 0.0
      dWx(i, j) = 0.0
      dWy(i, j) = 0.0
      W_rhs(i, j) = 0.0
      Ra_dTdx(i, j) = 0.0
    end do

    ! initialize temperature T
    if (Tinit == 'cosine') then
      do concurrent (i=1:nx, j=1:ny)
        T(i, j) = 0.5 * ( 1.0 + cos((3.0*pi*(i-1)*h)/width) )
      end do
    else
      call random_number(T)
    end if

    ! randomly initialize vorticity W
    call random_number(W)

    ! open velocity file
    open(v_file_id, file='v.bin', form='unformatted', access='stream')

  end subroutine initialize

  subroutine evolve(this, time)
    implicit none
    ! arguments
    class(simulation), intent(inout) :: this
    real, intent(in) :: time
    ! local variables
    real :: r_rms, v_max, c
    integer :: i, j

    ! solve for stream function S: D S = W (Dirichlet BCs = 0)
    call nvtx_push('MGsolve_S')
    r_rms = MGsolve_2DPoisson(S, W, h, 0.0, err, niters, .false.)
    call nvtx_pop()

    ! compute velocity field (vx, vy) from stream function S
    call nvtx_push('compute_velocity')
    call compute_velocity(S, h, h, vx, vy)
    call nvtx_pop()

    ! compute velocity magnitude v
    call nvtx_push('v_magnitude')
    v_max = 0.0
    do concurrent (i=2:nx-1, j=2:ny-1) reduce(max:v_max)
      v(i, j) = sqrt(vx(i, j)**2 + vy(i, j)**2)
      v_max = max(v_max, v(i, j))
    end do
    write(v_file_id) time, v_max
    call nvtx_pop()

    ! compute timestep dt
    this%dt = compute_dt(v_max)

    ! apply boundary conditions
    call apply_boundary_conditions(T)

    ! compute Ra * dT / dx
    call nvtx_push('Ra_dTdx')
    do concurrent (i=2:nx-1, j=2:ny-1)
      Ra_dTdx(i, j) = Ra * ( ( T(i+1, j) - T(i-1, j) ) / (2 * h) )
    end do
    call nvtx_pop()

    ! diffusion terms for temperature T and vorticity W
    if (beta /= 1.0) then
      call nvtx_push('diffusion')
      call diffusion2d(T, h, k, dT2)
      call diffusion2d(W, h, Pr, dW2)
      call nvtx_pop()
    end if

    ! advection terms for temperature T and vorticity W
    call nvtx_push('advection')
    call advection2dx(T, h, vx, dTx)
    call advection2dy(T, h, vy, dTy)
    call advection2dx(W, h, vx, dWx)
    call advection2dy(W, h, vy, dWy)
    call nvtx_pop()

    ! Euler step for temperature T and vorticity W
    if (beta > 0.0) then
      ! semi-implicit step for temperature T
      c = 1.0 / (beta * this%dt)
      call nvtx_push('T_rhs')
      do concurrent (i=2:nx-1, j=2:ny-1)
        T_rhs(i, j) = -c * ( T(i,j) + this%dt * ( (1.0 - beta) * dT2(i, j) - dTx(i, j) - dTy(i, j) ) )
      end do
      call nvtx_pop()
      call nvtx_push('MGsolve_T')
      r_rms = MGsolve_2DPoisson(T, T_rhs, h, c, err, niters, .true.)
      call nvtx_pop()
      ! semi-implicit step for vorticity W
      c = c / Pr
      call nvtx_push('W_rhs')
      do concurrent (i=2:nx-1, j=2:ny-1)
        W_rhs(i, j) = -c * ( W(i, j) + this%dt * ( (1.0 - beta) * dW2(i, j) - dWx(i, j) - dWy(i, j) - Pr * Ra_dTdx(i, j) ) )
      end do
      call nvtx_pop()
      call nvtx_push('MGsolve_W')
      r_rms = MGsolve_2DPoisson(W, W_rhs, h, c, err, niters, .false.)
      call nvtx_pop()
    else
      ! explicit step for temperature T and vorticity W
      call nvtx_push('explicit_step')
      do concurrent (i=2:nx-1, j=2:ny-1)
        T(i, j) = T(i, j) + this%dt * ( dT2(i, j) - dTx(i, j) - dTy(i, j) )
        W(i, j) = W(i, j) + this%dt * ( dW2(i, j) - dWx(i, j) - dWy(i, j) - Pr * Ra_dTdx(i, j) )
      end do
      call nvtx_pop()
    end if

  end subroutine evolve

  real function compute_dt(v_max) result(dt_)
    implicit none
    ! arguments
    real, intent(in)  :: v_max
    ! local variables
    real    :: max_absvx, max_absvy
    integer :: i, j

    if (v_max == 0.0) then
      ! set timestep equal to diffusive timestep (avoid division by zero)
      dt_ = dt_dif
    else
      ! compute advective timestep
      max_absvx = 0.0
      max_absvy = 0.0
      do concurrent (i=2:nx-1, j=2:ny-1) reduce(max:max_absvx) reduce(max:max_absvy)
        max_absvx = max(max_absvx, abs(vx(i, j)))
        max_absvy = max(max_absvy, abs(vy(i, j)))
      end do
      dt_adv = a_adv * min(h / max_absvx, h / max_absvy)
      ! compute timestep
      if (beta >= 0.5) then
        dt_ = dt_adv ! allow big diffusive timesteps
      else
        dt_ = min(dt_dif, dt_adv)
      end if
    end if

  end function compute_dt

  subroutine export_results()
    implicit none

    ! print final S to file
    open(S_file_id, file='S.bin', form='unformatted', access='stream')
    write (S_file_id) nx, ny, S
    close(S_file_id)

    ! print final W to file
    open(W_file_id, file='W.bin', form='unformatted', access='stream')
    write (W_file_id) nx, ny, W
    close(W_file_id)

    ! print final T to file
    open(T_file_id, file='T.bin', form='unformatted', access='stream')
    write (T_file_id) nx, ny, T
    close(T_file_id)

  end subroutine export_results

  subroutine destroy()
    implicit none

    ! close velocity file
    close(v_file_id)
    ! deallocate field arrays
    deallocate(S, vx, vy, v)
    deallocate(T, dT2, dTx, dTy, T_rhs)
    deallocate(W, dW2, dWx, dWy, W_rhs)
    deallocate(Ra_dTdx)

  end subroutine destroy

end module m_sim
