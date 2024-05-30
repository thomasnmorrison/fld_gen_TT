module params
#include "macros.h"
  integer, parameter :: dl = kind(1.d0)
  real(dl), parameter :: twopi = 6.2831853071795864769252867665590

  integer, parameter :: nx=16, ny=16, nz=16
  integer, parameter :: nnx=nx/2+1, nny=ny/2+1, nnz=nz/2+1
  integer, parameter :: nn = min(nnx,nny,nnz)

  real(dl), parameter :: nvol = dble(nx)*dble(ny)*dble(nz)
  real(dl), parameter :: len = 2._dl
  real(dl), parameter :: dx = len/dble(nx)
  real(dl), parameter :: dk = twopi / len
  real(dl), parameter :: mpl = 1.e5

  integer, parameter :: nfld = NFIELD

  integer, parameter :: kos = 2**2            ! oversampling of k modes
  integer, parameter :: nkos = kos*nn         ! number of oversampled modes
  real(dl), parameter :: dkos = dk/dble(kos)  ! spacing of over sampled k modes
  real(dl), parameter :: norm = sqrt(nvol/dx**3/mpl**2)          ! normalization factor for spectra
  
end module params
