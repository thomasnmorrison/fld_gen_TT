! fld_gen_tt.f90

! to do: set up spectra to be read from file

program fld_gen_tt
#include "macros.h"
  use omp_lib
  use params
  use vars
  use io
  use lat_init
  use filter_mod
  
  implicit none

  integer :: i,j
  integer :: k, ii, jj, kk
  complex(C_DOUBLE_COMPLEX), dimension(3,3) :: pol_w, pol_sum
  integer :: seed = 1  ! RNG seed
  real(dl), dimension(2) :: kfilt = dk*(/dble(nn-1)/3._dl, dble(nn-1)/2._dl/)
  
  integer :: terror    ! error code
  
  terror = fftw_init_threads()
  print*,"Error code is ",terror
  terror = omp_get_max_threads()
  print*,"Num Threads = ",terror
  call fftw_plan_with_nthreads(omp_get_max_threads())

  call init_fft_ptrs()
  call allocate_ar()
  
!  call init_input()
  call init_output()

  ! initialialize spec_l
  ! call read_spec(spec_l)
  call init_spec_invar(spec_l, 1._dl)
  ! initialialize spec_r
  spec_r = spec_l
  
  call lat_init_tt(spec_l, spec_r, kos, (nn-1)*kos, seed, fldtt, f1, fk1, planb, norm)

  do i=1,nfld
     do j=1,6
        f1(:,:,:) = fldtt(i,j,IRANGE)
        call lat_filter(bowler_hat_cubic, kfilt, f1, fk1, dk, planf, planb)
        fldtt(i,j,IRANGE) = f1(:,:,:)
     end do
  end do

  ! Make output
!  call summary_stat()
  call make_output()

!  pol_sum = (0._dl, 0._dl)
!  do k=1,nz; if (k>nnz) then; kk = k-nz-1; else; kk=k-1; endif
!     do j=1,ny; if (j>nny) then; jj = j-ny-1; else; jj=j-1; endif 
!        do i=1,nx; if (i>nnx) then; ii = i-nx-1; else; ii=i-1; endif   
!           if (sqrt(dble(ii**2 + jj**2 + kk**2)) .gt. nn/2._dl) then
!              cycle
!           end if
!           if (ii .eq. 0 .and. jj .eq. 0 .and. kk .eq. 0) then
!              cycle
!           end if
!           call get_pol_l(pol_w, (/ii,jj,kk/))
!           pol_sum = pol_sum + pol_w * conjg(pol_w)
!        end do
!     end do
!  end do

!  print*, pol_sum/nvol

contains

  ! Subroutine to initialize a scale invariant spectrum
  subroutine init_spec_invar(spec, amp)
    real(dl) :: spec(:,:,:)
    real(dl) :: amp
    integer :: i

    do i=1,nkos
       spec(i,:,:) = amp/(i*dkos)**3
    end do
  end subroutine init_spec_invar
  
end program fld_gen_tt
