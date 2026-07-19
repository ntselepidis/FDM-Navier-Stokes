module m_nvtx
#ifdef NVHPC
  use iso_c_binding
#endif
  implicit none

#ifdef NVHPC
  interface
    subroutine nvtxRangePushA(name) bind(C, name='nvtxRangePushA')
      use iso_c_binding
      character(kind=c_char) :: name(*)
    end subroutine nvtxRangePushA

    subroutine nvtxRangePop() bind(C, name='nvtxRangePop')
    end subroutine nvtxRangePop
  end interface
#endif
contains

  subroutine nvtx_push(name)
    character(len=*), intent(in) :: name
#ifdef NVHPC
    call nvtxRangePushA(trim(name) // c_null_char)
#endif
  end subroutine nvtx_push

  subroutine nvtx_pop()
#ifdef NVHPC
    call nvtxRangePop()
#endif
  end subroutine nvtx_pop

end module m_nvtx
