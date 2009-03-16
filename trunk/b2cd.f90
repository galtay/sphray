!> \file b2cd.f90

!> \brief Handles the impact parameter to column depth table

module b2cd_mod
use myf90_mod
implicit none
private

public :: b2cdfac, read_b2cd_file

integer(i4b) :: Nbins     !< impact parameter bins
integer(i4b) :: Nentries  !< table entries
real(r8b) :: db           !< spacing in impact parameter
real(r8b), allocatable :: b2cd_table(:)  !< stores line integrations


contains

!----------------------------------------------------
!>  read in the impact parameter -> column depth file
subroutine read_b2cd_file(verbose,b2cd_file)

  logical, intent(inout) :: verbose !< screen output if true
  character(clen), intent(in) :: b2cd_file  ! impact parameter to CD table

  integer(i8b) :: lun, err, i
  logical :: fthere

    verbose = .false.
    write(*,'(A,T29,A)') "using b2cd file: ", trim(b2cd_file)

    inquire(file=b2cd_file,exist=fthere)
    if (.not. fthere) then
       write(*,*) "cant find b2cd file: ", trim(b2cd_file)
       stop
    end if

    call open_formatted_file_r(b2cd_file,lun)
    read(lun,*) Nentries
    if (verbose) write(*,'(A,I5)') "b2cd table entries: ", Nentries

    allocate ( b2cd_table(Nentries), stat=err )
    if(err/=0) then
       write(*,*) "   b2cd.f90>> cant allocate b2cd_table"
       stop
    end if

    if (verbose) write(*,'(A)',advance="no") &
         "reading in SPH impact parameter -> column depth table ..."
    do i = 1,Nentries
       read(lun,*) b2cd_table(i)
    end do
    Nbins = Nentries - 1
    db = 1.0d0 / Nbins

    if (verbose) write(*,*) "done"
    close(lun)
    call sum_b2cd_table(verbose)
    verbose = .true.

end subroutine read_b2cd_file

!----------------------------------------------------------------
!> checks normalization of impact parameter to column depth table
subroutine sum_b2cd_table(verbose)
use physical_constants_mod, only: pi

  logical, intent(in) :: verbose  !< screen output if true

  real(r8b) :: sum,b
  integer(i4b) :: i

    sum = 0.0d0
  
    do i = 0,Nbins-1
       b = (i+0.5) * db 
       sum = sum + 2.0d0 * pi * b * b2cd_table(i) * db 
    end do

    if (verbose) write(*,'(A,ES10.4)') "normalization of b2cd table (should be close to 1.0) = ", sum

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
