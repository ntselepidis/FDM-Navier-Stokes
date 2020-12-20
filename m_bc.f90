module m_bc
  ! no variables in this module
contains

  subroutine apply_boundary_conditions(T)
    implicit none
    real, intent(inout) :: T(:, :)

    call apply_dirichlet_boundary_conditions(T)
    call apply_neumann_boundary_conditions(T)

  end subroutine apply_boundary_conditions

  subroutine apply_dirichlet_boundary_conditions(T)
    implicit none
    real, intent(inout) :: T(:, :)
    integer :: ny
    ny = size(T, 2)

    ! apply Dirichlet BCs on bottom and top boundaries
    T(:, 1) = 1.0
    T(:, ny) = 0.0

  end subroutine apply_dirichlet_boundary_conditions

  subroutine apply_neumann_boundary_conditions(T)
    implicit none
    real, intent(inout) :: T(:, :)
    integer :: nx
    nx = size(T, 1)

    ! apply Neumann BCs on left and right boundaries
    T(1, :) = T(2, :)
    T(nx, :) = T(nx-1, :)

  end subroutine apply_neumann_boundary_conditions

end module m_bc
