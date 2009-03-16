!> \file spectra.f90

!> \brief spectra handling
!<

module spectra_mod
use myf90_mod
implicit none
private

  public :: rn2freq, read_spectra_file

  integer(i8b), parameter :: maxfreqbins = 1000  !< maximum frequency bins

  integer(i8b) :: Nspectra      !< number of user defined spectra
  integer(i8b) :: Nfreqs        !< number of entries

  real(r8b), allocatable :: nus(:,:)  !< frequency in HI ionizing units (SpecNum,freq)
  real(r8b), allocatable :: cdf(:,:)  !< cumulative distribution function (SpecNum,cdf)

contains


!>   reads all user defined spectra 
!---------------------------------------
  subroutine read_spectra_file(verbose,spectra_file)

    logical, intent(inout) :: verbose !< report to screen? 
    character(clen), intent(in) :: spectra_file  !< file containing spectra

    real(r8b) :: Numin, Numax
    integer(i8b) :: lun, err, i, j
    logical :: fthere

    verbose = .false.
    write(*,'(A,T29,A)') "using spectra file: ", trim(spectra_file)

    inquire(file=spectra_file,exist=fthere)
    if (.not. fthere) then
       write(*,'(2A)') "cant find spectra file: ", trim(spectra_file)
       stop
    end if

    call open_formatted_file_r(spectra_file,lun)
    read(lun,*) Nspectra
    if (Nspectra<=0) Nspectra=1  ! just to make sure spectra is allocated
    if (verbose) write(*,'(A,I5)') "number of user defined spectra = ", Nspectra

    100 format(A,I5,A,I5,A) 

    do i = 1,Nspectra

       read(lun,*) Nfreqs, Numin, Numax

       if (verbose) then
          write(*,100) "user defined spectrum ",i," has ", Nfreqs," entries"
          write(*,*) 
       end if

       if (i == 1) then
          allocate ( nus(Nfreqs,Nspectra), cdf(Nfreqs,Nspectra), stat=err )
          if(err/=0) then
             write(*,*) "   spectra.f90>> can't allocate nus and cdf"
             stop
          end if
       end if

       do j = 1,Nfreqs
          read(lun,*) nus(j,i), cdf(j,i)
       end do

    end do

    close(lun)
    if (verbose) write(*,*) 
    verbose = .true.

  end subroutine read_spectra_file


  !> returns a frequency given a spectral type
  !-------------------------------------------
  function rn2freq(spectype) result(freq)
  use mt19937_mod, only: genrand_real1

    real(r4b) :: spectype  !< spectral user defined type
    real(r4b) :: freq      !< returned frequency

    real(r8b) :: rn
    integer(i8b) :: specbin
    integer(i8b) :: specnum

    if (spectype <= 0.0) then
       freq = abs(spectype)
       return
    else
       specnum = int(spectype)
       rn = genrand_real1()
       do specbin = 1,Nfreqs
          if ( rn <= cdf(specbin,specnum) ) exit
       end do
       freq = nus(specbin,specnum)  
    end if

   end function rn2freq

end module spectra_mod
