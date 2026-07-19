module m_fdm
  ! no variables in this module
contains

  subroutine diffusion2d(T, h, k, dT2)
    implicit none
    ! arguments
    real, intent(in) :: h, k
    real, intent(in) :: T(:, :)
    ! result
    real, intent(out) :: dT2(:, :)
    ! local variables
    integer :: i, j
    integer :: nx, ny
    nx = size(T, 1)
    ny = size(T, 2)

    ! apply 2D diffusion
    do concurrent (i=2:nx-1, j=2:ny-1)
      dT2(i, j) = k * ( T(i+1, j) + T(i-1, j) + T(i, j+1) + T(i, j-1) - 4.0 * T(i, j) ) / h**2
    end do

  end subroutine diffusion2d

  subroutine advection2dx(T, hx, vx, dTx)
    implicit none
    ! arguments
    real, intent(in) :: hx
    real, intent(in) :: T(:, :)
    real, intent(in) :: vx(:, :)
    ! result
    real, intent(out) :: dTx(:, :)
    ! local variables
    integer :: i, j
    integer :: nx, ny
    nx = size(T, 1)
    ny = size(T, 2)

    do concurrent (i=2:nx-1, j=2:ny-1)
      if (vx(i, j) > 0) then
        dTx(i, j) = vx(i, j) * (T(i, j) - T(i-1, j)) / hx
      else
        dTx(i, j) = vx(i, j) * (T(i+1, j) - T(i, j)) / hx
      end if
    end do

  end subroutine advection2dx

  subroutine advection2dy(T, hy, vy, dTy)
    implicit none
    ! arguments
    real, intent(in) :: hy
    real, intent(in) :: T(:, :)
    real, intent(in) :: vy(:, :)
    ! result
    real, intent(out) :: dTy(:, :)
    ! local variables
    integer :: i, j
    integer :: nx, ny
    nx = size(T, 1)
    ny = size(T, 2)

    do concurrent (i=2:nx-1, j=2:ny-1)
      if (vy(i, j) > 0) then
        dTy(i, j) = vy(i, j) * (T(i, j) - T(i, j-1)) / hy
      else
        dTy(i, j) = vy(i, j) * (T(i, j+1) - T(i, j)) / hy
      end if
    end do

  end subroutine advection2dy

end module m_fdm
