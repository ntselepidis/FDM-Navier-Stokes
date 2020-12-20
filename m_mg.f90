module m_mg
  use m_bc
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
    do j = 1, ny
      do i = 1, nx
        f_rms = f_rms + f(i, j)**2
      end do
    end do
    f_rms = sqrt(f_rms / (nx * ny))

    tolf = tol * f_rms

    r_rms = 0.0
    do iter = 1, niters
      ! apply Dirichlet and Neumann BCs for temperature T
      if (apply_BCs) then
        call apply_boundary_conditions(u)
      end if
      ! execute V-cycle iteration
      r_rms = Vcycle_2DPoisson(u, f, h, c, apply_BCs)
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
  function iteration_2DPoisson(u, f, h, c, alpha) result(r_rms)
    implicit none
    ! arguments
    real, intent(inout) :: u(:, :)
    real, intent(in)    :: f(:, :)
    real, intent(in)    :: h
    real, intent(in)    :: c
    real, intent(in)    :: alpha
    ! local variables
    real                :: r, r_rms
    integer             :: nx, ny, i, j

    nx = size(u, 1)
    ny = size(u, 2)

    ! Gauss-Seidel iteration
    r_rms = 0.0
    do j = 2, ny - 1
      do i = 2, nx - 1

        ! compute residual
        r = ( u(i+1, j) &
            + u(i-1, j) &
            + u(i, j+1) &
            + u(i, j-1) &
            - (4.0 + c * h**2) * u(i, j) ) / h**2 - f(i, j)

        ! update u
        u(i, j) = u(i, j) + alpha * (h**2 / (4.0 + c * h**2)) * r

        ! update r_rms
        r_rms = r_rms + r**2

      end do
    end do

    r_rms = sqrt(r_rms / (nx * ny))

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
    integer           :: nx, ny

    nx = size(u, 1)
    ny = size(u, 2)

    res(2:nx-1, 2:ny-1) = ( u(3:nx,   2:ny-1) &
                          + u(1:nx-2, 2:ny-1) &
                          + u(2:nx-1, 3:ny  ) &
                          + u(2:nx-1, 1:ny-2) &
                          - (4.0 + c * h**2) * u(2:nx-1, 2:ny-1) ) / h**2 &
                          - f(2:nx-1, 2:ny-1)

  end subroutine

  ! copies every other point in fine into coarse
  subroutine restrict(fine, coarse, apply_BCs)
    implicit none
    ! arguments
    real, intent(in)  :: fine(:, :)
    real, intent(out) :: coarse(:, :)
    logical, intent(in) :: apply_BCs
    ! local variables
    integer           :: nx, ny
    integer           :: i, j, ic, jc
    real, parameter   :: a4  = 1.0 / 4.0
    real, parameter   :: a8  = 1.0 / 8.0
    real, parameter   :: a16 = 1.0 / 16.0

    nx = size(fine, 1)
    ny = size(fine, 2)

    ! initalize coarse to zero (+ Dirichlet(0) BCs)
    coarse = 0.0

    ! apply restriction stencil
    jc = 2
    do j = 3, ny-2, 2
      ic = 2
      do i = 3, nx-2, 2
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
        ic = ic + 1
      end do
      jc = jc + 1
    end do

    ! apply Neumann BCs for temperature T
    if (apply_BCs) then
      call apply_neumann_boundary_conditions(coarse)
    end if

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

    nx = size(fine, 1)
    ny = size(fine, 2)

    ! initialize fine to zero (+ Dirichlet(0) BCs)
    fine = 0.0

    jc = 2
    do j = 3, ny-2, 2
      ic = 2
      do i = 3, nx-2, 2
        fine(i, j) = coarse(ic, jc)
        fine(i+1, j) = fine(i+1, j) + a2 * coarse(ic, jc)
        fine(i-1, j) = fine(i-1, j) + a2 * coarse(ic, jc)
        fine(i, j+1) = fine(i, j+1) + a2 * coarse(ic, jc)
        fine(i, j-1) = fine(i, j-1) + a2 * coarse(ic, jc)
        fine(i+1, j+1) = fine(i+1, j+1) + a4 * coarse(ic, jc)
        fine(i+1, j-1) = fine(i+1, j-1) + a4 * coarse(ic, jc)
        fine(i-1, j+1) = fine(i-1, j+1) + a4 * coarse(ic, jc)
        fine(i-1, j-1) = fine(i-1, j-1) + a4 * coarse(ic, jc)
        ic = ic + 1
      end do
      jc = jc + 1
    end do

    ! apply Neumann BCs for temperature T
    if (apply_BCs) then
      call apply_neumann_boundary_conditions(fine)
    end if

  end subroutine prolongate

  ! Vcycle multigrid
  recursive function Vcycle_2DPoisson(u_f, rhs, h, c, apply_BCs) result (resV)
    implicit none
    real resV
    ! arguments
    real, intent(inout) :: u_f(:, :)
    real, intent(in)    :: rhs(:, :), h, c
    logical, intent(in) :: apply_BCs
    ! local variables
    integer             :: nx, ny, nxc, nyc, i
    real, allocatable   :: res_c(:, :), corr_c(:, :), res_f(:, :), corr_f(:, :)
    real                :: alpha=1.0, res_rms

    nx = size(u_f, 1); ny = size(u_f, 2)  ! must be power of 2 plus 1

    if( nx-1 /= 2*((nx-1)/2) .or. ny-1 /= 2*((ny-1)/2) ) then
      stop 'ERROR:not a power of 2'
    end if

    nxc = 1+(nx-1)/2; nyc = 1+(ny-1)/2  ! coarse grid size

    if (min(nx, ny) > 5) then  ! not the coarsest level

      allocate(res_f(nx, ny), corr_f(nx, ny), &
        corr_c(nxc, nyc), res_c(nxc, nyc))

      !---------- take 2 iterations on the fine grid--------------
      res_rms = iteration_2DPoisson(u_f, rhs, h, c, alpha)
      res_rms = iteration_2DPoisson(u_f, rhs, h, c, alpha)

      !--------- restrict the residual to the coarse grid --------
      call residual_2DPoisson(u_f, rhs, h, c, res_f)
      call restrict(res_f, res_c, apply_BCs)

      !---------- solve for the coarse grid correction -----------
      corr_c = 0.
      res_rms = Vcycle_2DPoisson(corr_c, res_c, h*2, c, apply_BCs) ! *RECURSIVE CALL*

      !---- prolongate (interpolate) the correction to the fine grid
      call prolongate(corr_c, corr_f, apply_BCs)

      !---------- correct the fine-grid solution -----------------
      u_f = u_f - corr_f

      !---------- two more smoothing iterations on the fine grid---
      res_rms = iteration_2DPoisson(u_f, rhs, h, c, alpha)
      res_rms = iteration_2DPoisson(u_f, rhs, h, c, alpha)

      deallocate(res_f, corr_f, res_c, corr_c)

    else

      !----- coarsest level (ny=5): iterate to get 'exact' solution

      do i = 1, 100
        res_rms = iteration_2DPoisson(u_f, rhs, h, c, alpha)
      end do

    end if

    resV = res_rms   ! returns the rms. residual

  end function Vcycle_2DPoisson

end module m_mg
