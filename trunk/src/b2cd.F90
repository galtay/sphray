!> \file b2cd.F90

!> \brief Handles the impact parameter to column depth table

module b2cd_mod
use myf90_mod
implicit none
private

public :: b2cdfac
public :: read_b2cd_file

integer(i4b) :: Nbins     !< impact parameter bins
integer(i4b) :: Nentries  !< table entries
real(r8b) :: db           !< spacing in impact parameter
real(r8b), allocatable :: b2cd_table(:)  !< stores line integrations


contains

!----------------------------------------------------
!>  read in the impact parameter -> column depth file
subroutine read_b2cd_file(b2cd_file)

  character(clen), intent(in) :: b2cd_file  ! impact parameter to CD table

  integer, parameter :: verb = 1
  character(clen) :: str
  integer(i8b) :: lun, err, i
  logical :: fthere

  character(clen), parameter :: myname = 'read_b2cd_file'
  logical, parameter :: crash = .true.

  write(str,'(A,T27,A)') 'using b2cd file: ', trim(b2cd_file)
  call mywrite(str, verb) 

  inquire(file=b2cd_file, exist=fthere)
  if (.not. fthere) then
     call myerr('cant find b2cd file: ' // trim(b2cd_file), myname, crash)
  end if

  call open_formatted_file_r(b2cd_file,lun)
  read(lun,*) Nentries
  write(str,'(A,I5)') '  b2cd table entries:    ', Nentries
  call mywrite('', verb+1)
  call mywrite(str, verb+1) 

  allocate ( b2cd_table(Nentries), stat=err )
  if(err/=0) call myerr('cant allocate b2cd_table', myname, crash)
     
  str = '  reading in SPH impact parameter -> column depth table'
  call mywrite(str, verb+1)

  do i = 1,Nentries
     read(lun,*) b2cd_table(i)
  end do
  Nbins = Nentries - 1
  db = 1.0d0 / Nbins

  close(lun)
  call sum_b2cd_table()


end subroutine read_b2cd_file

!----------------------------------------------------------------
!> checks normalization of impact parameter to column depth table
subroutine sum_b2cd_table()
use physical_constants_mod, only: pi

  real(r8b) :: sum,b
  integer(i4b) :: i

  integer, parameter :: verb=2
  character(clen) :: str

  sum = 0.0d0
  
  do i = 1,Nbins-1
     b = (i+0.5) * db 
     sum = sum + 2.0d0 * pi * b * b2cd_table(i) * db 
  end do

  write(str,'(A,ES10.4)') &
       "  normalization of b2cd table (should be close to 1.0) = ", sum
  
  call mywrite(str, verb) 
  call mywrite('', verb)
  
end subroutine sum_b2cd_table

!> scales the table entry to the column depth factor using the 
!! smoothing length, the code cgs length, and the impact parameter
!! normalized to the smoothing length
function b2cdfac(b_norm,hsml,cgs_len) result(cdfac)

  real(r8b) :: b_norm  !< normalized impact parameter
  real(r8b) :: hsml    !< smoothing length
  real(r8b) :: cgs_len !< cgs length
  real(r8b) :: cdfac   !< output column depth factor

  integer :: bindx

  bindx = floor(Nbins * b_norm) + 1
  if (bindx.gt.Nentries) bindx = Nentries
  cdfac = b2cd_table(bindx) / ( hsml*hsml*cgs_len*cgs_len )

end function b2cdfac


end module b2cd_mod
