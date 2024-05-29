! vars.f90

module vars
#include "macros.h"
  use params
  use fftw3

  implicit none

  real(dl), allocatable :: fldtt(:,:,:,:,:)  ! tensor field
  real(dl), allocatable :: spec_l(:,:,:), spec_r(:,:,:)  ! spectra for L/R polarization/helicity (f1,f2,mode)
  
  type(C_PTR) :: planf, planb
  real(C_DOUBLE), pointer :: f1(:,:,:), f2(:,:,:), f3(:,:,:)!, f4(:,:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: Fk1(:,:,:), Fk2(:,:,:), Fk3(:,:,:), Fk4(:,:,:)

contains

  ! Subroutine to initialize fft plans and arrays
  subroutine init_fft_ptrs()
    call allocate_3d_array(nx, ny, nz, f1, Fk1)
    call allocate_3d_array(nx, ny, nz, f2, Fk2)
    call allocate_3d_array(nx, ny, nz, f3, Fk3)
    call allocate_3d_fourier_array(nx, ny, nz, Fk4)
    planf = fftw_plan_dft_r2c_3d(nz, ny, nx, f1, Fk1, FFTW_MEASURE+FFTW_DESTROY_INPUT) 
    planb = fftw_plan_dft_c2r_3d(nz, ny, nx, Fk1, f1, FFTW_MEASURE+FFTW_DESTROY_INPUT) 
  end subroutine init_fft_ptrs

  ! Subroutine to allocate array for tensor field variables
  subroutine allocate_ar()
    allocate(fldtt(nfld,6,IRANGE))
    allocate(spec_l(nkos,nfld,nfld), spec_r(nkos,nfld,nfld))
  end subroutine allocate_ar

end module vars
