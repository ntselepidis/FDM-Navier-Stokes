module m_fdm
  ! no variables in this module
contains

  function diffusion2d(T, h, k)
    implicit none
    ! arguments
    real, intent(in) :: h, k
    real, intent(in) :: T(:, :)
    ! result
    real, dimension(size(T, 1), size(T, 2)) :: diffusion2d
    ! local variables
    integer :: nx, ny
    nx = size(T, 1)
    ny = size(T, 2)

    ! apply 2D diffusion
    diffusion2d(2:nx-1, 2:ny-1) = k * ( T(3:nx,   2:ny-1) &
                                      + T(1:nx-2, 2:ny-1) &
                                      + T(2:nx-1, 3:ny  ) &
                                      + T(2:nx-1, 1:ny-2) &
                                      - 4.0 * T(2:nx-1, 2:ny-1) ) / h**2

  end function diffusion2d

  function advection2dx(T, hx, vx)
    implicit none
    ! arguments
    real, intent(in) :: hx
    real, intent(in) :: T(:, :)
    real, intent(in) :: vx(:, :)
    ! result
    real, dimension(size(T, 1), size(T, 2)) :: advection2dx
    ! local variables
    integer :: i, j
    integer :: nx, ny
    nx = size(T, 1)
    ny = size(T, 2)

    do concurrent (i=2:nx-1, j=2:ny-1)
      if (vx(i, j) > 0) then
        advection2dx(i, j) = vx(i, j) * (T(i, j) - T(i-1, j)) / hx
      else
        advection2dx(i, j) = vx(i, j) * (T(i+1, j) - T(i, j)) / hx
      end if
    end do

  end function advection2dx

  function advection2dy(T, hy, vy)
    implicit none
    ! arguments
    real, intent(in) :: hy
    real, intent(in) :: T(:, :)
    real, intent(in) :: vy(:, :)
    ! result
    real, dimension(size(T, 1), size(T, 2)) :: advection2dy
    ! local variables
    integer :: i, j
    integer :: nx, ny
    nx = size(T, 1)
    ny = size(T, 2)

    do concurrent (i=2:nx-1, j=2:ny-1)
      if (vy(i, j) > 0) then
        advection2dy(i, j) = vy(i, j) * (T(i, j) - T(i, j-1)) / hy
      else
        advection2dy(i, j) = vy(i, j) * (T(i, j+1) - T(i, j)) / hy
      end if
    end do

  end function advection2dy

end module m_fdm
