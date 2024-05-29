! tensor_mod.f90

! Module for dealing with 3x3 tensors. Right now supports $T^i_j$ type tensors (actually need to double check that if metric is not diagonal).

module tensor_mod
#include "macros.h"
  use fftw3
  use params
  
  implicit none

contains

  ! Subroutine to express components of a (^3,_3) spatial tensor in a new coordinate system with $\hat{z}$
  ! rotated onto $\hat{v}$.
  ! to do: document exactly what rotation is being considered
  subroutine tensor_rot(s, v)
    complex(C_DOUBLE_COMPLEX) :: s(:,:)  ! components of symmetric spatial tensor
    integer, intent(in) :: v(:)

    complex(C_DOUBLE_COMPLEX) :: c1, c2, s1, s2  ! respectively $\cos(\theta)$, $\cos(\phi)$, $\sin(\theta)$, $\sin(\phi)$ with angles as in my notes
    complex(C_DOUBLE_COMPLEX), dimension(3,3) :: rot  ! rotation matrix (from x to x'')
    complex(C_DOUBLE_COMPLEX), dimension(3,3) :: temp  ! temp array for LAPACK subroutines

    if (v(1) == 0 .and. v(2) == 0) then
       return  ! no rotation if v is already aligned along z-axis
    else
       c1 = dble(v(1))/sqrt(dble(v(1)**2 + v(2)**2)) * (1._dl, 0._dl)
       s1 = dble(v(2))/sqrt(dble(v(1)**2 + v(2)**2)) * (1._dl, 0._dl)
       c2 = dble(v(3))/sqrt(dble(v(1)**2 + v(2)**2 + v(3)**2)) * (1._dl, 0._dl)
       s2 = sqrt(1._dl - c2**2) * (1._dl, 0._dl)
       rot = reshape((/c2*c1, -s1, s2*c1, c2*s1, c1, s2*s1, -s2, (0._dl,0._dl), c2/),(/3,3/))
       
       ! call matrix multiplication (note use of symmetric matrix multiplication)
       call zsymm('R', 'U', 3, 3, (1._dl,0._dl), s, 3, rot, 3, (0._dl,0._dl), temp, 3)
       s = temp
       ! call matrix multiplication
       call zgemm('N', 'T', 3, 3, 3, (1._dl,0._dl), s, 3, rot, 3, (0._dl,0._dl), temp, 3)
       s = temp
    end if
  end subroutine tensor_rot
  
end module tensor_mod
