! lat_init.f90

! to do: check that complex is actually double complex
! to do: put tensor_rot into a module

module lat_init
#include "macros.h"
  use fftw3
  use params
  use vars
  use grv_mod
  use tensor_mod
  
  implicit none

  complex(C_DOUBLE_COMPLEX), parameter, dimension(3,3) :: pol_plus = reshape((/(1.,0.),(0.,0.),(0.,0.), &
                                                             (0.,0.), (-1.,0.), (0.,0.), &
                                                             (0.,0.), (0.,0.), (0.,0.)/), (/3,3/))
  complex(C_DOUBLE_COMPLEX), parameter, dimension(3,3) :: pol_cross = reshape((/(0.,0.),(1.,0.),(0.,0.), &
                                                              (1.,0.), (0.,0.), (0.,0.), &
                                                              (0.,0.), (0.,0.), (0.,0.)/), (/3,3/))
  complex(C_DOUBLE_COMPLEX), parameter, dimension(3,3) :: pol_l = sqrt(0.5)*(pol_plus - (0.,1.)*pol_cross)
  complex(C_DOUBLE_COMPLEX), parameter, dimension(3,3) :: pol_r = sqrt(0.5)*(pol_plus + (0.,1.)*pol_cross)
  
contains

  ! to do: it appears the polarization weighting is bugged (need to check tensor rotation)
  subroutine lat_init_tt(spec_l, spec_r, kos, kcut, seed, fldtt_init, f, fk, planb, norm)
    real(dl), intent(in) :: spec_l(:,:,:), spec_r(:,:,:)
    integer, intent(in) :: kos, kcut, seed
    real(dl) :: fldtt_init(:,:,:,:,:)
    real(C_DOUBLE), pointer :: f(:,:,:)
    complex(C_DOUBLE_COMPLEX), pointer :: fk(:,:,:)
    type(C_PTR) :: planb  ! fftw c2r plan
    real(dl), intent(in) :: norm

    complex(C_DOUBLE_COMPLEX), dimension(nfld,nfld) :: spec_interp
    complex(C_DOUBLE_COMPLEX), allocatable :: fk_mat_tt(:,:,:,:,:)
    complex(C_DOUBLE_COMPLEX), dimension(nfld) :: grv

    complex(C_DOUBLE_COMPLEX), dimension(3,3) :: pol_w ! polarization weights
    real(dl) :: rad
    integer :: i,j,k,ii,jj,kk,l

    call init_rng(seed)  ! initialize rng
    
    allocate(fk_mat_tt(nfld,6,nnx,ny,nz))  ! allocate a complex array to hold Fourier transforms

    ! loop over wavenumber
    do k=1,nz; if (k>nnz) then; kk = k-nz-1; else; kk=k-1; endif
       do j=1,ny; if (j>nny) then; jj = j-ny-1; else; jj=j-1; endif 
          do i=1,nnx; ii=i-1
             rad = sqrt(dble(ii**2)+dble(jj**2)+dble(kk**2))  ! find wavenumber magnitude
             l = floor(rad*kos)  ! get wavenumber index for supplied spectra
             if (l == 0 .or. l .gt. kcut) then  ! zero the zero/cut mode and cycle the loop
                fk_mat_tt(:,:,LATIND) = (0._dl, 0._dl)! set Fourier transform to zero for this wavenumber  
                cycle
             end if
             ! interpolate spec_l
             spec_interp(:,:) = (1._dl,0._dl)*((1._dl + l - rad*kos)*spec_l(l,:,:) + (rad*kos - l)*spec_l(l+1,:,:))
             ! get L polarization index weighting
             call get_pol_l(pol_w, (/ii,jj,kk/))
             ! random draw for L polarization
             call zpotrf('L',nfld,spec_interp,nfld,l)  ! Choleski decomposition
             grv(:) = get_grv_complex(nfld)
             call ztrmv('L','N','N', nfld, spec_interp, nfld, grv, 1)  ! multiply grv and Choleski
             fk_mat_tt(:,1,LATIND) = pol_w(1,1)*grv(:)
             fk_mat_tt(:,2,LATIND) = pol_w(1,2)*grv(:)
             fk_mat_tt(:,3,LATIND) = pol_w(1,3)*grv(:)
             fk_mat_tt(:,4,LATIND) = pol_w(2,2)*grv(:)
             fk_mat_tt(:,5,LATIND) = pol_w(2,3)*grv(:)
             fk_mat_tt(:,6,LATIND) = pol_w(3,3)*grv(:)
             ! interpolate spec_r
             spec_interp(:,:) = (1._dl,0._dl)*((1._dl + l - rad*kos)*spec_r(l,:,:) + (rad*kos - l)*spec_r(l+1,:,:))
             ! get R polarization index weighting
             call get_pol_r(pol_w, (/ii,jj,kk/))
             ! random draw for R polarization
             call zpotrf('L',nfld,spec_interp,nfld,l)  ! Choleski decomposition
             grv(:) = get_grv_complex(nfld)
             call ztrmv('L','N','N', nfld, spec_interp, nfld, grv, 1)  ! multiply grv and Choleski
             fk_mat_tt(:,1,LATIND) = pol_w(1,1)*grv(:)
             fk_mat_tt(:,2,LATIND) = pol_w(1,2)*grv(:)
             fk_mat_tt(:,3,LATIND) = pol_w(1,3)*grv(:)
             fk_mat_tt(:,4,LATIND) = pol_w(2,2)*grv(:)
             fk_mat_tt(:,5,LATIND) = pol_w(2,3)*grv(:)
             fk_mat_tt(:,6,LATIND) = pol_w(3,3)*grv(:)
          end do
       end do
    end do
             
    ! set ii=0 components to be complex conjugate
    do k=1,nz; if (k.ne.1) then; kk = nz+2-k; else; kk=k; endif
       do j=2,nny; jj = ny+2-j
          fk_mat_tt(:,:,1,jj,kk) = conjg(fk_mat_tt(:,:,1,j,k))
       end do
    end do

    do k=2,nnz; kk = nz+2-k
       fk_mat_tt(:,:,1,1,kk) = conjg(fk_mat_tt(:,:,1,1,k))
    end do

    ! invert FFT and add fluctuations to homogeneous fields
    fk_mat_tt = fk_mat_tt*norm/nvol
    do i=1,nfld
       do j=1,6
          fk(:,:,:) = fk_mat_tt(i,j,:,:,:)
          call fftw_execute_dft_c2r(planb, fk, f)
          fldtt_init(i,j,IRANGE) = f(:,:,:)
       end do
    end do
    
    deallocate(fk_mat_tt)  ! de-allocate Fourier transform array
  end subroutine lat_init_tt

  subroutine get_pol_plus(s, v)
    complex(C_DOUBLE_COMPLEX) :: s(:,:)  ! components of symmetric spatial tensor
    integer, intent(in) :: v(:)

    complex(C_DOUBLE_COMPLEX) :: c1, c2, s1, s2  ! respectively $\cos(\theta)$, $\cos(\phi)$, $\sin(\theta)$, $\sin(\phi)$ with angles as in my notes

    if (v(1) .eq. 0 .and. v(2) .eq. 0) then
       c1 = (1._dl, 0._dl)
       s1 = (0._dl, 0._dl)
    else
       c1 = dble(v(1))/sqrt(dble(v(1)**2 + v(2)**2)) * (1._dl, 0._dl)
       s1 = dble(v(2))/sqrt(dble(v(1)**2 + v(2)**2)) * (1._dl, 0._dl)
    end if
    c2 = dble(v(3))/sqrt(dble(v(1)**2 + v(2)**2 + v(3)**2)) * (1._dl, 0._dl)
    s2 = sqrt(1._dl - c2**2) * (1._dl, 0._dl)
    s = reshape((/c2**2*c1**2 - s1**2, (c2**2 + (1._dl,0._dl))*c1*s1, -c1*c2*s2, &
         (c2**2 + (1._dl,0._dl))*c1*s1, c2**2*s1**2 - c1**2, -c2*s1*s2, &
         -c1*c2*s2, -c2*s1*s2, s2**2/),(/3,3/))
  end subroutine get_pol_plus

  subroutine get_pol_cross(s,v)
    complex(C_DOUBLE_COMPLEX) :: s(:,:)  ! components of symmetric spatial tensor
    integer, intent(in) :: v(:)

    complex(C_DOUBLE_COMPLEX) :: c1, c2, s1, s2  ! respectively $\cos(\theta)$, $\cos(\phi)$, $\sin(\theta)$, $\sin(\phi)$ with angles as in my notes

    if (v(1) .eq. 0 .and. v(2) .eq. 0) then
       c1 = (1._dl, 0._dl)
       s1 = (0._dl, 0._dl)
    else
       c1 = dble(v(1))/sqrt(dble(v(1)**2 + v(2)**2)) * (1._dl, 0._dl)
       s1 = dble(v(2))/sqrt(dble(v(1)**2 + v(2)**2)) * (1._dl, 0._dl)
    end if
    c2 = dble(v(3))/sqrt(dble(v(1)**2 + v(2)**2 + v(3)**2)) * (1._dl, 0._dl)
    s2 = sqrt(1._dl - c2**2) * (1._dl, 0._dl)
    s = reshape((/-2._dl*c1*c2*s1, c2*(c1**2-s1**2), s1*s2, &
         c2*(c1**2 - s1**2), 2._dl*c1*c2*s1, -c1*s2, &
         s1*s2, -c1*s2, (0._dl,0._dl)/),(/3,3/))
  end subroutine get_pol_cross

  subroutine get_pol_l(s,v)
    complex(C_DOUBLE_COMPLEX) :: s(:,:)  ! components of symmetric spatial tensor
    integer, intent(in) :: v(:)

    complex(C_DOUBLE_COMPLEX) :: temp(3,3)

    call get_pol_plus(s,v)
    call get_pol_cross(temp,v)
    s = sqrt(0.5_dl)*(s - (0._dl,1._dl)*temp)
  end subroutine get_pol_l

  subroutine get_pol_r(s,v)
    complex(C_DOUBLE_COMPLEX) :: s(:,:)  ! components of symmetric spatial tensor
    integer, intent(in) :: v(:)

    complex(C_DOUBLE_COMPLEX) :: temp(3,3)
    
    call get_pol_plus(s,v)
    call get_pol_cross(temp,v)
    s = sqrt(0.5_dl)*(s + (0._dl,1._dl)*temp)
  end subroutine get_pol_r
  
end module lat_init
