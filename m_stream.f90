module m_stream
  ! no variables in this module
contains

  subroutine compute_stream_function(nx, ny, hx, hy, xmin, xmax, ymin, ymax, B, S)
    implicit none
    ! arguments
    integer, intent(in) :: nx, ny
    real, intent(in) :: hx, hy
    real, intent(in) :: xmin, xmax, ymin, ymax
    real, intent(in) :: B
    real, intent(out) :: S(:, :)
    ! local variables
    real, parameter :: pi = 3.141592653589793238462643383279502884
    real :: x = 0.0, y = 0.0
    integer :: i = 0, j = 0

    do concurrent (i=1:nx, j=1:ny)
      x = xmin + (i-1) * hx
      y = ymin + (j-1) * hy
      S(i, j) = B * sin((pi * x) / xmax) * sin((pi * y) / ymax)
    end do

  end subroutine compute_stream_function

  subroutine compute_velocity(S, hx, hy, vx, vy)
    implicit none
    ! arguments
    real, intent(in) :: hx, hy
    real, intent(in) :: S(:, :)
    real, intent(out) :: vx(:, :)
    real, intent(out) :: vy(:, :)
    ! local variables
    integer :: i, j
    integer :: nx, ny
    nx = size(S, 1)
    ny = size(S, 2)

    ! Initialize vx and vy
    do concurrent (i=1:nx, j=1:ny)
      vx(i, j) = 0.0
      vy(i, j) = 0.0
    end do

    ! Compute velocity at each interior point using centered finite differences
    do concurrent (i=2:nx-1, j=2:ny-1)
      vx(i, j) =  (S(i, j+1) - S(i, j-1)) / (2 * hy)
      vy(i, j) = -(S(i+1, j) - S(i-1, j)) / (2 * hx)
    end do

  end subroutine compute_velocity

end module m_stream
