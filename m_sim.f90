module m_sim
  use m_bc,     only: apply_boundary_conditions
  use m_fdm
  use m_mg,     only: MGsolve_2DPoisson
  use m_stream, only: compute_velocity
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
    integer :: i

    h = 1.0 / ( ny - 1.0 ) ! mesh size
    width = ( nx - 1.0 ) / ( ny - 1.0 ) ! width

    ! compute diffusive timestep
    dt_dif = ( a_dif * min(h, h)**2 ) / max(k, Pr)

    ! allocate field arrays
    allocate(S(nx, ny), vx(nx, ny), vy(nx, ny), v(nx, ny))
    allocate(T(nx, ny), dT2(nx, ny), dTx(nx, ny), dTy(nx, ny), T_rhs(nx, ny))
    allocate(W(nx, ny), dW2(nx, ny), dWx(nx, ny), dWy(nx, ny), W_rhs(nx, ny))
    allocate(Ra_dTdx(nx, ny))

    S = 0.0; vx = 0.0; vy = 0.0; v = 0.0
    T = 0.0; dT2 = 0.0; dTx = 0.0; dTy = 0.0; T_rhs = 0.0
    W = 0.0; dW2 = 0.0; dWx = 0.0; dWy = 0.0; W_rhs = 0.0
    Ra_dTdx = 0.0

    ! initialize temperature T
    if (Tinit == 'cosine') then
      do concurrent (i=1:nx)
        T(i, :) = 0.5 * ( 1.0 + cos((3.0*pi*(i-1)*h)/width) )
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

    ! solve for stream function S: D S = W (Dirichlet BCs = 0)
    r_rms = MGsolve_2DPoisson(S, W, h, 0.0, err, niters, .false.)

    ! compute velocity field (vx, vy) from stream function S
    call compute_velocity(S, h, h, vx, vy)

    ! compute velocity magnitude v
    v = sqrt(vx**2 + vy**2)
    v_max = maxval(v)
    write(v_file_id) time, v_max

    ! compute timestep dt
    this%dt = compute_dt(v_max)

    ! apply boundary conditions
    call apply_boundary_conditions(T)

    ! compute Ra * dT / dx
    Ra_dTdx(2:nx-1, 2:ny-1) = Ra * ( ( T(3:nx, 2:ny-1) - T(1:nx-2, 2:ny-1) ) / (2 * h) )

    ! diffusion terms for temperature T and vorticity W
    if (beta /= 1.0) then
      dT2 = diffusion2d(T, h, k)
      dW2 = diffusion2d(W, h, Pr)
    end if

    ! advection terms for temperature T and vorticity W
    dTx = advection2dx(T, h, vx)
    dTy = advection2dy(T, h, vy)
    dWx = advection2dx(W, h, vx)
    dWy = advection2dy(W, h, vy)

    ! Euler step for temperature T and vorticity W
    if (beta > 0.0) then
      ! semi-implicit step for temperature T
      c = 1.0 / (beta * this%dt)
      T_rhs = -c * ( T + this%dt * ( (1.0 - beta) * dT2 - dTx - dTy ) )
      r_rms = MGsolve_2DPoisson(T, T_rhs, h, c, err, niters, .true.)
      ! semi-implicit step for vorticity W
      c = c / Pr
      W_rhs = -c * ( W + this%dt * ( (1.0 - beta) * dW2 - dWx - dWy - Pr * Ra_dTdx ) )
      r_rms = MGsolve_2DPoisson(W, W_rhs, h, c, err, niters, .false.)
    else
      ! explicit step for temperature T and vorticity W
      T = T + this%dt * ( dT2 - dTx - dTy )
      W = W + this%dt * ( dW2 - dWx - dWy - Pr * Ra_dTdx )
    end if

  end subroutine evolve

  real function compute_dt(v_max) result(dt_)
    implicit none
    ! arguments
    real, intent(in)  :: v_max

    if (v_max == 0.0) then
      ! set timestep equal to diffusive timestep (avoid division by zero)
      dt_ = dt_dif
    else
      ! compute advective timestep
      dt_adv = a_adv * min(h / maxval(vx), h / maxval(vy))
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
