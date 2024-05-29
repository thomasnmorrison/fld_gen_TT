! io.f90

! to do: decide on input format for spectra
! to do: write subroutine for reading spectra

module io
#include "macros.h"
  use params
  use vars
  
  implicit none

  character(len=*), parameter :: ident = "output"   ! identifier appended to output fiel names
  character(len=*), parameter :: f_spec = "spec_file"  ! file name for spectra

  ! File numbers
  integer, parameter :: nfile_temp = 99
  integer, parameter :: nfile_spec = 98
  integer, parameter :: nfile_fldtt = 97
  
contains

  ! to do: specify file format
  subroutine init_input()
    open(unit=nfile_spec,file=f_spec,status="old")
  end subroutine init_input

  ! to do: name output file
  subroutine init_output()
    open(unit=nfile_fldtt,file="fld_TT"//ident//".out",form="unformatted",access="stream")
  end subroutine init_output

  subroutine make_output()
    write(nfile_fldtt) fldtt(:,:,IRANGE)
  end subroutine make_output
  
end module io
