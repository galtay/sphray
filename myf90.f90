!> \file myf90.f90

!> \brief my fortran 90 module
!!
!! this module contains preconnected logical unit numbers (lun), 
!! portable precision definitions and subroutines for reading and 
!! parsing files 
!<

module myf90_mod
implicit none

  public

  interface scanfile
     module procedure scanfile_r4b, scanfile_r8b, scanfile_i4b, &
                      scanfile_i8b, scanfile_chr, scanfile_log
  end interface

  integer, parameter :: stderr = 0 !< preconnected standard error lun 
  integer, parameter :: stdin  = 5 !< preconnected standard in lun
  integer, parameter :: stdout = 6 !< preconnected standard output lun

! selected_int_kind(r)     -10^r < n < 10^r
!------------------------------------------------------------------------

  integer, parameter :: i1b = selected_int_kind(2)  !< 1 byte integer type
  integer, parameter :: i2b = selected_int_kind(4)  !< 2 byte integer type
  integer, parameter :: i4b = selected_int_kind(9)  !< 4 byte integer type
  integer, parameter :: i8b = selected_int_kind(18) !< 8 byte integer type

! selected_real_kind(p,r)  p is decimal precision, r is exponent range
!------------------------------------------------------------------------
  integer, parameter :: r4b  = selected_real_kind(6,37)    !< 4 byte real type

  integer, parameter :: r8b  = selected_real_kind(15,307)  !< 8 byte real type
  integer, parameter :: r16b = selected_real_kind(33,4931) !< 16 byte real type

! set default precision
!------------------------
  integer, parameter :: rdp = r4b
  integer, parameter :: idp = i4b

  integer, parameter :: clen  = 200 !< default character variable length


contains


!-----------------------------------------------
!> reports an error message and stops execution
   subroutine myerr(str,routine,crash)
     character(*), intent(in) :: str
     character(*), intent(in) :: routine
     logical, intent(in) :: crash

     write(*,*) "*****************************"
     write(*,*) "error detected"
     write(*,*) "err string: ", trim(str)
     write(*,*) "routine   : ", trim(routine)
     write(*,*) "*****************************"
     if (crash) stop

   end subroutine myerr


!------------------------------------
!> returns a free logical unit number
   function get_free_lun() result(lun)

    integer(i8b) :: lun                        !< free lun
    integer(i8b) :: i                          !< loop counter
    integer(i8b), parameter :: minlun = 110    !< min lun to check
    integer(i8b), parameter :: maxlun = 1000   !< max lun to check
    logical :: badlun                     !< true if lun is already connected
 
    do i = minlun,maxlun
      inquire(unit=i, opened=badlun)
      if (.not. badlun) then
         lun = i
         return
      end if
    end do

    write(*,*) '  get_free_lun> no free logical unit numbers'
    write(*,*) '  checked luns from, ', minlun, ' to ', maxlun
    stop 

  end function get_free_lun


!----------------------------------------
!> returns the current working directory
  function get_current_dir() result(path)
 
    character(clen) :: path         !< the current working dir
    integer(i8b) :: lun                  !< lun for temp file
    character(clen) :: temp_file    !< name of temp file
    character(clen) :: system_call  !< string passed to system

    temp_file = 'get_current_dir_temp_file'
    lun = get_free_lun() 
    
    write(system_call,100) trim(temp_file)
    100 format ("pwd > ", A)
    call system (system_call)

    call open_formatted_file_r(temp_file,lun)
    read(lun,'(A)') path
    close(lun)

    write(system_call,200) trim(temp_file)
    200 format ("rm -f ", A)
    call system(system_call)

  end function get_current_dir


!--------------------------------------------------------
!> opens a formatted file for reading and returns the lun
  subroutine open_formatted_file_r(filename,lun)
  
    character(*), intent(in) :: filename !< name of file to open
    integer(i8b), intent(out) :: lun          !< lun 
    logical :: fthere

    inquire(file=filename,exist=fthere)
    if (.not. fthere) then
       write(*,*) "cant find file: ", trim(filename)
       stop
    end if

    lun = get_free_lun() 
    open(unit=lun, file=filename, action="read")

  end subroutine open_formatted_file_r


!--------------------------------------------------------
!> opens a formatted file for writing and returns the lun
  subroutine open_formatted_file_w(filename,lun)
  
    character(*), intent(in) :: filename !< name of file to open
    integer(i8b), intent(out) :: lun          !< lun 

    lun = get_free_lun()
    open(unit=lun, file=filename, action="write")

  end subroutine open_formatted_file_w

!-------------------------------------------------------------
!> opens an unformatted file for reading and returns the lun
  subroutine open_unformatted_file_r(filename,lun)
  
    character(*), intent(in) :: filename  !< name of file to open
    integer(i8b), intent(out)  :: lun          !< lun
    logical :: fthere

    inquire(file=filename,exist=fthere)
    if (.not. fthere) then
       write(*,*) "cant find file: ", trim(filename)
       stop "file error"
    end if

    lun = get_free_lun()
    open(unit=lun, file=filename, form='unformatted', action="read")

  end subroutine open_unformatted_file_r

!--------------------------------------------------------------------
!> opens an unformatted file for stream reading and returns the lun
  subroutine open_unformatted_file_sr(filename,lun)
  
    character(*), intent(in) :: filename  !< name of file to open
    integer(i8b), intent(out)  :: lun          !< lun
    logical :: fthere

    inquire(file=filename,exist=fthere)
    if (.not. fthere) then
       write(*,*) "cant find file: ", trim(filename)
       stop
    end if

    lun = get_free_lun()
    open(unit=lun, file=filename, form="unformatted", &
         action="read", access="stream", position="rewind")

  end subroutine open_unformatted_file_sr

!------------------------------------------------------------
!> opens an unformatted file for writing and returns the lun
  subroutine open_unformatted_file_w(filename,lun)
  
    character(*), intent(in) :: filename  !< name of file to open
    integer(i8b), intent(out)  :: lun             !< lun

    lun = get_free_lun()
    open(unit=lun, file=filename, form='unformatted', &
         action="write")
   
  end subroutine open_unformatted_file_w

!-------------------------------------------------------------------
!> opens an unformatted file for stream writing and returns the lun
  subroutine open_unformatted_file_sw(filename,lun)
  
    character(*), intent(in) :: filename  !< name of file to open
    integer(i8b), intent(out)  :: lun             !< lun

    lun = get_free_lun()
    open(unit=lun, file=filename, form='unformatted', &
         action="write", access="stream", position="rewind")

  end subroutine open_unformatted_file_sw



  ! scan file routines, one for each type
  !====================================================================

  subroutine scanfile_r4b(filename,keyword,var) 
    character(*), intent(in) :: filename
    character(*), intent(in) :: keyword
    real(r4b), intent(out) :: var
    character(26), parameter :: tmpfile = "file_scanner_temp_file.out"
    character(200) :: cmnd

    integer(i8b) :: lun, Nfound
        
    ! produce awk command and send it to the system
    !-----------------------------------------------
    100 format("awk '{if ($1 ~ /^",A,"$/) print $2}' ",A, " > ",A)
    write(cmnd,100) trim(keyword), trim(filename), tmpfile
    call system(cmnd)
    
    Nfound = scanfile_count_args(tmpfile,filename,keyword)

    ! read variable from temporary file
    !-----------------------------------
    call open_formatted_file_r(tmpfile,lun)
    read(lun,*) var
    close(lun)

    ! remove temporary file
    !-----------------------
    200 format("rm -f ", A)
    write(cmnd,200) tmpfile
    call system(cmnd)

  end subroutine scanfile_r4b


  subroutine scanfile_r8b(filename,keyword,var) 
    character(*), intent(in) :: filename
    character(*), intent(in) :: keyword
    real(r8b), intent(out) :: var
    character(26), parameter :: tmpfile = "file_scanner_temp_file.out"
    character(200) :: cmnd

    integer(i8b) :: lun, Nfound
    
    ! produce awk command and send it to the system
    !-----------------------------------------------
    100 format("awk '{if ($1 ~ /^",A,"$/) print $2}' ",A, " > ",A)
    write(cmnd,100) trim(keyword), trim(filename), tmpfile
    call system(cmnd)

    Nfound = scanfile_count_args(tmpfile,filename,keyword)
    
    ! read variable from temporary file
    !-----------------------------------
    call open_formatted_file_r(tmpfile,lun)
    read(lun,*) var
    close(lun)

    ! remove temporary file
    !-----------------------
    200 format("rm -f ", A)
    write(cmnd,200) tmpfile
    call system(cmnd)

  end subroutine scanfile_r8b


  subroutine scanfile_i4b(filename,keyword,var) 
    character(*), intent(in) :: filename
    character(*), intent(in) :: keyword
    integer(i4b), intent(out) :: var
    character(26), parameter :: tmpfile = "file_scanner_temp_file.out"
    character(200) :: cmnd

    integer(i8b) :: lun, Nfound
        
    ! produce awk command and send it to the system
    !-----------------------------------------------
    100 format("awk '{if ($1 ~ /^",A,"$/) print $2}' ",A, " > ",A)
    write(cmnd,100) trim(keyword), trim(filename), tmpfile
    call system(cmnd)
    
    Nfound = scanfile_count_args(tmpfile,filename,keyword)

    ! read variable from temporary file
    !-----------------------------------
    call open_formatted_file_r(tmpfile,lun)
    read(lun,*) var
    close(lun)

    ! remove temporary file
    !-----------------------
    200 format("rm -f ", A)
    write(cmnd,200) tmpfile
    call system(cmnd)

  end subroutine scanfile_i4b



  subroutine scanfile_i8b(filename,keyword,var) 
    character(*), intent(in) :: filename
    character(*), intent(in) :: keyword
    integer(i8b), intent(out) :: var
    character(26), parameter :: tmpfile = "file_scanner_temp_file.out"
    character(200) :: cmnd
    
    integer(i8b) :: lun, Nfound
    
    
    ! produce awk command and send it to the system
    !-----------------------------------------------
    100 format("awk '{if ($1 ~ /^",A,"$/) print $2}' ",A, " > ",A)
    write(cmnd,100) trim(keyword), trim(filename), tmpfile
    call system(cmnd)
    
    Nfound = scanfile_count_args(tmpfile,filename,keyword)

    ! read variable from temporary file
    !-----------------------------------
    call open_formatted_file_r(tmpfile,lun)
    read(lun,*) var
    close(lun)

    ! remove temporary file
    !-----------------------
    200 format("rm -f ", A)
    write(cmnd,200) tmpfile
    call system(cmnd)

  end subroutine scanfile_i8b



  subroutine scanfile_chr(filename,keyword,var)
    character(*), intent(in) :: filename
    character(*), intent(in) :: keyword
    character(clen), intent(out) :: var
    character(26), parameter :: tmpfile = "file_scanner_temp_file.out"
    character(200) :: cmnd

    integer(i8b) :: lun, Nfound    
    
    ! produce awk command and send it to the system
    !-----------------------------------------------
    100 format("awk '{if ($1 ~ /^",A,"$/) print $2}' ",A, " > ",A)
    write(cmnd,100) trim(keyword), trim(filename), tmpfile
    call system(cmnd)
    
    Nfound = scanfile_count_args(tmpfile,filename,keyword)

    ! read variable from temporary file
    !-----------------------------------
    call open_formatted_file_r(tmpfile,lun)
    read(lun,'(A)') var
    close(lun)

    ! remove temporary file
    !-----------------------
    200 format("rm -f ", A)
    write(cmnd,200) tmpfile
    call system(cmnd)

  end subroutine scanfile_chr



  subroutine scanfile_log(filename,keyword,var)
    character(*), intent(in) :: filename
    character(*), intent(in) :: keyword
    logical, intent(out) :: var
    character(26), parameter :: tmpfile = "file_scanner_temp_file.out"
    character(200) :: cmnd

    integer(i8b) :: lun, Nfound
        
    ! produce awk command and send it to the system
    !-----------------------------------------------
    100 format("awk '{if ($1 ~ /^",A,"$/) print $2}' ",A, " > ",A)
    write(cmnd,100) trim(keyword), trim(filename), tmpfile
    call system(cmnd)
    
    Nfound = scanfile_count_args(tmpfile,filename,keyword)

    ! read variable from temporary file
    !-----------------------------------
    call open_formatted_file_r(tmpfile,lun)
    read(lun,*) var
    close(lun)

    ! remove temporary file
    !-----------------------
    200 format("rm -f ", A)
    write(cmnd,200) tmpfile
    call system(cmnd)

  end subroutine scanfile_log


  ! retunrs the number of matchs awk found
  function scanfile_count_args(awkoutfile,paramfile,keyword) result(N)
    character(*), intent(in) :: awkoutfile
    character(*), intent(in) :: paramfile
    character(*), intent(in) :: keyword
    integer(i8b) :: N
    character(27), parameter :: tmpfile = "file_scanner_temp_file.out1"
    character(200) :: cmnd
    integer(i8b) :: lun
 
    100 format("awk 'END {print NR}' ",A, " > ",A)
    write(cmnd,100) awkoutfile, tmpfile
    call system(cmnd)

    ! read variable from temporary file
    !-----------------------------------
    call open_formatted_file_r(tmpfile,lun)
    read(lun,*) N
    close(lun)

    ! remove temporary file
    !-----------------------
    200 format("rm -f ", A)
    write(cmnd,200) tmpfile
    call system(cmnd)

    if (N==1) return

    write(*,*) " *** error while parsing file ", trim(paramfile)
        
    if (N==0) then
       write(*,*) "keyword ", trim(keyword), " does not appear in file"
    else if (N > 1) then
       write(*,*) "keyword ", trim(keyword), " is defined ", N , " times"
    end if

    ! remove other temporary file
    !-----------------------------
    write(cmnd,200) awkoutfile
    call system(cmnd)

    stop

  end function scanfile_count_args


end module myf90_mod


