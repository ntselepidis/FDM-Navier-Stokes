module m_mg
  use m_bc
  use m_nvtx, only: nvtx_push, nvtx_pop
  ! no variables in this module
contains

  ! solves $(\nabla^2 - c) u = f$ using V-cycle multigrid
  function MGsolve_2DPoisson(u, f, h, c, tol, niters, apply_BCs) result(r_rms)
    implicit none
    ! arguments
    real, intent(inout) :: u(:, :)
    real, intent(in)    :: f(:, :)
    real, intent(in)    :: h
    real, intent(in)    :: c
    real, intent(in)    :: tol
    integer, intent(in) :: niters
    logical, intent(in) :: apply_BCs
    ! local variables
    real                :: r_rms, f_rms, tolf
    integer             :: nx, ny, iter, i, j

    ! set nx, ny
    nx = size(u, 1)
    ny = size(u, 2)

    ! compute f_rms
    f_rms = 0.0
    do concurrent (i=2:nx-1, j=2:ny-1) reduce(+:f_rms)
      f_rms = f_rms + f(i,j)**2
    end do
    f_rms = sqrt(f_rms / (nx*ny))

    tolf = tol * f_rms

    r_rms = 0.0
    do iter = 1, niters
      ! apply Dirichlet and Neumann BCs for temperature T
      if (apply_BCs) then
        call apply_boundary_conditions(u)
      end if
      ! execute V-cycle iteration
      call nvtx_push('Vcycle')
      r_rms = Vcycle_2DPoisson(u, f, h, c, apply_BCs, 0)
      call nvtx_pop()
      ! print*, iter, r_rms / f_rms
      if (r_rms < tolf) then
        ! print*, 'V-cycle multigrid converged in', iter, 'iterations.'
        exit
      end if
    end do

    if (r_rms > tolf) then
      print*, 'V-cycle multigrid failed to converge within', niters, 'iterations.'
    end if

  end function MGsolve_2DPoisson

  ! performs one Gauss-Seidel iteration on field u
  ! returns rms residual
  function iteration_2DPoisson(u, f, h, c, alpha, res, compute_norm) result(r_rms)
    implicit none
    ! arguments
    real, intent(inout) :: u(:, :)
    real, intent(in)    :: f(:, :)
    real, intent(in)    :: h
    real, intent(in)    :: c
    real, intent(in)    :: alpha
    real, intent(inout) :: res(:, :) ! local but preallocated
    logical, intent(in) :: compute_norm
    ! local variables
    real                :: r_rms
    integer             :: nx, ny, i, j

    call nvtx_push('iteration_2DPoisson')

    nx = size(u, 1)
    ny = size(u, 2)

    ! compute residual
    call residual_2DPoisson(u, f, h, c, res)

    ! update u
    do concurrent (i=2:nx-1, j=2:ny-1)
        u(i, j) = u(i, j) + alpha * (h**2 / (4.0 + c * h**2)) * res(i, j)
    end do

    if (compute_norm) then
      r_rms = 0.0
      do concurrent (i=2:nx-1, j=2:ny-1) reduce(+:r_rms)
        r_rms = r_rms + res(i,j)**2
      end do
      r_rms = sqrt(r_rms / (nx * ny))
    else
      r_rms = -1.0
    end if

    call nvtx_pop()

  end function iteration_2DPoisson

  ! computes the residual $R = (\nabla^2 - c) u - f$ in array res
  subroutine residual_2DPoisson(u, f, h, c, res)
    implicit none
    ! arguments
    real, intent(in)  :: u(:, :)
    real, intent(in)  :: f(:, :)
    real, intent(in)  :: h
    real, intent(in)  :: c
    real, intent(out) :: res(:, :)
    ! local variables
    integer           :: i, j
    integer           :: nx, ny

    call nvtx_push('residual_2DPoisson')

    nx = size(u, 1)
    ny = size(u, 2)

    do concurrent (i=2:nx-1, j=2:ny-1)
      res(i, j) = ( u(i+1, j) + u(i-1, j) + u(i, j+1) + u(i, j-1) - (4.0 + c * h**2) * u(i, j) ) / h**2 - f(i, j)
    end do

    call nvtx_pop()

  end subroutine

  ! copies every other point in fine into coarse
  subroutine restrict(fine, coarse, apply_BCs)
    implicit none
    ! arguments
    real, intent(in)  :: fine(:, :)
    real, intent(out) :: coarse(:, :)
    logical, intent(in) :: apply_BCs
    ! local variables
    integer           :: nx, ny, nxc, nyc
    integer           :: i, j, ic, jc
    real, parameter   :: a4  = 1.0 / 4.0
    real, parameter   :: a8  = 1.0 / 8.0
    real, parameter   :: a16 = 1.0 / 16.0

    call nvtx_push('restrict')

    nx = size(fine, 1)
    ny = size(fine, 2)

    nxc = size(coarse, 1)
    nyc = size(coarse, 2)

    ! initalize coarse to zero (+ Dirichlet(0) BCs)
    do concurrent (i=1:nxc, j=1:nyc)
      coarse(i, j) = 0.0
    end do

    ! apply restriction stencil
    do concurrent (i=3:nx-2, j=3:ny-2)
      if ((mod(i, 2) == 1) .and. (mod(j, 2) == 1)) then
        ic  = (i-1)/2+1
        jc  = (j-1)/2+1
        coarse(ic, jc) = fine(i, j)
        ! coarse(ic, jc) = a4 * fine(i, j) &
        !                + a8 * fine(i+1, j) &
        !                + a8 * fine(i-1, j) &
        !                + a8 * fine(i, j+1) &
        !                + a8 * fine(i, j-1) &
        !                + a16 * fine(i+1, j+1) &
        !                + a16 * fine(i+1, j-1) &
        !                + a16 * fine(i-1, j+1) &
        !                + a16 * fine(i-1, j-1)
      end if
    end do

    ! apply Neumann BCs for temperature T
    if (apply_BCs) then
      call apply_neumann_boundary_conditions(coarse)
    end if

    call nvtx_pop()

  end subroutine restrict

  ! copies coarse into every other point in fine
  ! fill the other points using linear interpolation
  subroutine prolongate(coarse, fine, apply_BCs)
    implicit none
    ! arguments
    real, intent(in)    :: coarse(:, :)
    real, intent(out)   :: fine(:, :)
    logical, intent(in) :: apply_BCs
    ! local variables
    integer             :: nx, ny
    integer             :: i, j, ic, jc
    real, parameter     :: a2 = 1.0 / 2.0
    real, parameter     :: a4 = 1.0 / 4.0

    call nvtx_push('prolongate')

    nx = size(fine, 1)
    ny = size(fine, 2)

    ! initialize fine to zero (+ Dirichlet(0) BCs)
    do concurrent (i=1:nx, j=1:ny)
      fine(i, j) = 0.0
    end do

    ! Needs atomics !!!
    do concurrent (i=3:nx-2, j=3:ny-2)
      if ((mod(i, 2) == 1) .and. (mod(j, 2) == 1)) then
        ic  = (i-1)/2+1
        jc  = (j-1)/2+1
        fine(i, j) = coarse(ic, jc)
        !$acc atomic update
        fine(i+1, j) = fine(i+1, j) + a2 * coarse(ic, jc)
        !$acc atomic update
        fine(i-1, j) = fine(i-1, j) + a2 * coarse(ic, jc)
        !$acc atomic update
        fine(i, j+1) = fine(i, j+1) + a2 * coarse(ic, jc)
        !$acc atomic update
        fine(i, j-1) = fine(i, j-1) + a2 * coarse(ic, jc)
        !$acc atomic update
        fine(i+1, j+1) = fine(i+1, j+1) + a4 * coarse(ic, jc)
        !$acc atomic update
        fine(i+1, j-1) = fine(i+1, j-1) + a4 * coarse(ic, jc)
        !$acc atomic update
        fine(i-1, j+1) = fine(i-1, j+1) + a4 * coarse(ic, jc)
        !$acc atomic update
        fine(i-1, j-1) = fine(i-1, j-1) + a4 * coarse(ic, jc)
      end if
    end do

    ! apply Neumann BCs for temperature T
    if (apply_BCs) then
      call apply_neumann_boundary_conditions(fine)
    end if

    call nvtx_pop()

  end subroutine prolongate

  ! Vcycle multigrid
  recursive function Vcycle_2DPoisson(u_f, rhs, h, c, apply_BCs, level) result (resV)
    implicit none
    real resV
    ! arguments
    real, intent(inout) :: u_f(:, :)
    real, intent(in)    :: rhs(:, :), h, c
    logical, intent(in) :: apply_BCs
    integer, intent(in) :: level
    ! local variables
    integer             :: nx, ny, nxc, nyc, i, j, iter
    real, allocatable   :: res_c(:, :), corr_c(:, :), res_f(:, :), corr_f(:, :)
    real                :: alpha=4.0/5.0, res_rms
    character(len=16)   :: label

    write(label, '(A,I0)') 'Vcycle_L', level
    call nvtx_push(trim(label))

    nx = size(u_f, 1); ny = size(u_f, 2)  ! must be power of 2 plus 1

    if( nx-1 /= 2*((nx-1)/2) .or. ny-1 /= 2*((ny-1)/2) ) then
      stop 'ERROR:not a power of 2'
    end if

    nxc = 1+(nx-1)/2; nyc = 1+(ny-1)/2  ! coarse grid size

    if (min(nx, ny) > 5) then  ! not the coarsest level

      allocate(res_f(nx, ny), corr_f(nx, ny), &
        corr_c(nxc, nyc), res_c(nxc, nyc))

      !---------- take 2 iterations on the fine grid--------------
      res_rms = iteration_2DPoisson(u_f, rhs, h, c, alpha, res_f, .false.)
      res_rms = iteration_2DPoisson(u_f, rhs, h, c, alpha, res_f, .false.)

      !--------- restrict the residual to the coarse grid --------
      call residual_2DPoisson(u_f, rhs, h, c, res_f)
      call restrict(res_f, res_c, apply_BCs)

      !---------- solve for the coarse grid correction -----------
      do concurrent (i=1:nxc, j=1:nyc)
        corr_c(i, j) = 0.
      end do

      res_rms = Vcycle_2DPoisson(corr_c, res_c, h*2, c, apply_BCs, level+1) ! *RECURSIVE CALL*

      !---- prolongate (interpolate) the correction to the fine grid
      call prolongate(corr_c, corr_f, apply_BCs)

      !---------- correct the fine-grid solution -----------------
      do concurrent (i=1:nx, j=1:ny)
        u_f(i, j) = u_f(i, j) - corr_f(i, j)
      end do

      !---------- two more smoothing iterations on the fine grid---
      res_rms = iteration_2DPoisson(u_f, rhs, h, c, alpha, res_f, .false.)
      res_rms = iteration_2DPoisson(u_f, rhs, h, c, alpha, res_f, level==0)

      deallocate(res_f, corr_f, res_c, corr_c)

    else

      !----- coarsest level (ny=5): iterate to get 'exact' solution
      allocate(res_f(nx, ny))

      !do i = 1, 100
      !  res_rms = iteration_2DPoisson(u_f, rhs, h, c, alpha, res_f)
      !end do

      do concurrent (i=2:nx-1, j=2:ny-1)
        do iter = 1, 100
          res_f(i, j) = ( u_f(i+1, j) + u_f(i-1, j) + u_f(i, j+1) + u_f(i, j-1) - (4.0 + c * h**2) * u_f(i, j) ) / h**2 - rhs(i, j)
          u_f(i, j) = u_f(i, j) + alpha * (h**2 / (4.0 + c * h**2)) * res_f(i, j)
        end do
      end do

      res_rms = -1.0
      !res_rms = 0.0
      !do concurrent (i=2:nx-1, j=2:ny-1) reduce(+:res_rms)
      !  res_rms = res_rms + res_f(i,j)**2
      !end do
      !res_rms = sqrt(res_rms / (nx*ny))

      deallocate(res_f)

    end if

    resV = res_rms   ! returns the rms. residual

    call nvtx_pop()

  end function Vcycle_2DPoisson

end module m_mg
