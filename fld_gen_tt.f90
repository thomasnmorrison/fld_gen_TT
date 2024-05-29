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
  spec_l(:,:,:) = 1._dl  ! adjust with norm
  ! initialialize spec_r
  spec_r = spec_l
  
  call lat_init_tt(spec_l, spec_r, kos, nkos-1, seed, fldtt, f1, fk1, planb, norm)

  do i=1,nfld
     do j=1,6
        f1(:,:,:) = fldtt(i,j,IRANGE)
        call lat_filter(bowler_hat_cubic, kfilt, f1, fk1, dk, planf, planb)
        fldtt(i,j,IRANGE) = f1(:,:,:)
     end do
  end do

  ! Make output
  call make_output()
  
contains

  
end program fld_gen_tt
