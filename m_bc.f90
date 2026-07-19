module m_bc
  use m_nvtx, only: nvtx_push, nvtx_pop
  ! no variables in this module
contains

  subroutine apply_boundary_conditions(T)
    implicit none
    real, intent(inout) :: T(:, :)

    call nvtx_push('apply_boundary_conditions')
    call apply_dirichlet_boundary_conditions(T)
    call apply_neumann_boundary_conditions(T)
    call nvtx_pop()

  end subroutine apply_boundary_conditions

  subroutine apply_dirichlet_boundary_conditions(T)
    implicit none
    real, intent(inout) :: T(:, :)
    integer :: nx, ny, i
    nx = size(T, 1)
    ny = size(T, 2)

    ! apply Dirichlet BCs on bottom and top boundaries
    do concurrent (i=1:nx)
      T(i, 1) = 1.0
      T(i, ny) = 0.0
    end do

  end subroutine apply_dirichlet_boundary_conditions

  subroutine apply_neumann_boundary_conditions(T)
    implicit none
    real, intent(inout) :: T(:, :)
    integer :: nx, ny, j
    nx = size(T, 1)
    ny = size(T, 2)

    ! apply Neumann BCs on left and right boundaries
    do concurrent (j=1:ny)
      T(1, j) = T(2, j)
      T(nx, j) = T(nx-1, j)
    end do

  end subroutine apply_neumann_boundary_conditions

end module m_bc
