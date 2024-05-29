! filter_mod.f90

module filter_mod
#include "macros.h"
  use params
  use fftw3

  implicit none

  interface
     function filt_template(k_in, param_r) result(f)
       import
       real(dl) :: k_in
       real(dl) :: param_r(:)
       real(dl) :: f
     end function filt_template
  end interface
  
contains
  
  ! Subroutine to filter an input field
  subroutine lat_filter(filt, kfilt, f, fk, dk, planf, planb)
    procedure(filt_template) :: filt                 ! filter function
    real(dl), intent(in) :: kfilt(:)                 ! filter parameters
    real(C_DOUBLE), pointer :: f(:,:,:)
    complex(C_DOUBLE_COMPLEX), pointer :: fk(:,:,:)
    real(dl), intent(in) :: dk
    type(C_PTR) :: planf, planb

    integer :: i,j,k, ii, jj, kk
    real(dl) :: rad
    
    call fftw_execute_dft_r2c(planf, f, fk)
    do k=1,nz; if (k>nnz) then; kk = nz+1-k; else; kk=k-1; endif
       do j=1,ny; if (j>nny) then; jj = ny+1-j; else; jj=j-1; endif 
          do i=1,nnx; ii=i-1
             rad = sqrt(dble(ii**2)+dble(jj**2)+dble(kk**2))
             fk(IRANGE) = fk(IRANGE) * filt(dk*rad, kfilt)
          end do
       end do
    end do
    call fftw_execute_dft_c2r(planb, fk, f)
    f = f/nvol
  end subroutine lat_filter

  ! Function for a top_hat filter
  function top_hat(k_in, param_r) result(f)
    real(dl) :: k_in
    real(dl) :: param_r(:)
    real(dl) :: f

    f = 0.5_dl + sign(0.5_dl, param_r(1)-k_in)
  end function top_hat

  ! Function for a bowler hat filter with a cubic tail
  ! n.b. must have param_r(1) .lt. param_r(2)
  function bowler_hat_cubic(k_in, param_r) result(f)
    real(dl) :: k_in
    real(dl) :: param_r(:)
    real(dl) :: f

    real(dl) :: temp
    
    temp = (2._dl*k_in - (param_r(1) + param_r(2))) / ((param_r(2) - param_r(1)))
    f = 0.5_dl + sign(0.5_dl, param_r(1)-k_in) &
         + (sign(0.5_dl, k_in-param_r(1)) + sign(0.5_dl, param_r(2)-k_in)) &
         * (0.5_dl + 0.25_dl*(temp**3 - 3._dl*temp))
  end function bowler_hat_cubic
  
end module filter_mod
