! fftw_mod.f90

! this is copied over from J.B, check over and make consistent with my naming schemes

module fftw3
  use, intrinsic :: iso_c_binding
  include 'fftw3.f03'

contains
   
  subroutine allocate_3d_array(L,M,N, f, fk)
    integer :: L,M,N
    real(C_DOUBLE), pointer :: f(:,:,:)
    complex(C_DOUBLE_COMPLEX), pointer :: fk(:,:,:)
!    type(C_PTR) :: fptr, fkptr
!    integer :: LL

    call allocate_3d_real_array(L,M,N,f)
    call allocate_3d_fourier_array(L,M,N,fk)
  end subroutine allocate_3d_array

  subroutine allocate_3d_real_array(L,M,N,f)
    integer :: L,M,N
    real(C_DOUBLE), pointer :: f(:,:,:)
    type(C_PTR) :: fptr
    
    fptr = fftw_alloc_real(int(L*M*N, C_SIZE_T))
    call c_f_pointer(fptr,f,[L,M,N])
  end subroutine allocate_3d_real_array
  
  subroutine allocate_3d_fourier_array(L,M,N,fk)
    integer :: L,M,N
    complex(C_DOUBLE_COMPLEX), pointer :: fk(:,:,:)
    type(C_PTR) :: fkptr
    integer :: LL
    
    LL = L/2+1
    fkptr = fftw_alloc_complex(int(LL*M*N, C_SIZE_T))
    call c_f_pointer(fkptr, fk, [LL,M,N])
  end subroutine allocate_3d_fourier_array
  
end module fftw3
