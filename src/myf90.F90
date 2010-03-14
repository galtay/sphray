!> \file myf90.F90

!> \brief My Fortran 90 module
!!
!! Contains preconnected logical unit numbers (lun), 
!! portable precision definitions and subroutines for reading and 
!! parsing files 
!<

module myf90_mod
use iso_fortran_env    
implicit none

  public

  interface scanfile
     module procedure scanfile_r4b, scanfile_r8b, scanfile_i4b, &
                      scanfile_i8b, scanfile_chr, scanfile_log
  end interface


  integer, parameter :: stderr = error_unit  !< preconnected std error lun 
  integer, parameter :: stdin  = input_unit  !< preconnected std in lun
  integer, parameter :: stdout = output_unit !< preconnected std out lun

! selected_int_kind(r)     -10^r < n < 10^r
!------------------------------------------------------------------------

  integer, parameter :: i1b = selected_int_kind(2)  !< 1 byte integer type
  integer, parameter :: i2b = selected_int_kind(4)  !< 2 byte integer type
  integer, parameter :: i4b = selected_int_kind(9)  !< 4 byte integer type
  integer, parameter :: i8b = selected_int_kind(18) !< 8 byte integer type

! selected_real_kind(p,r)  p is decimal precision, r is exponent range
!------------------------------------------------------------------------
  integer, parameter :: r4b  = selected_real_kind(p=6,r=37)    !< 4 byte real type
  integer, parameter :: r8b  = selected_real_kind(p=15,r=307)  !< 8 byte real type
  integer, parameter :: r16b = selected_real_kind(p=33,r=4931) !< 16 byte real type


  integer, parameter :: clen  = 500 !< default character variable length

  integer :: myf90_verbosity  !< global verbosity threshold

!> command line type. 
!=========================
  type command_line_type
     integer :: nargs       !< number of arguments excluding executable
     integer :: len         !< number of chars in whole command line 
     character(clen) :: str !< the whole command line  
     character(clen), allocatable :: args(:)  !< individual command line args
     integer, allocatable :: arglens(:)       !< number of chars in each arg
  end type command_line_type


contains


!-----------------------------------------------
!> fills in the command line variable

  function initialize_command_line(verb) result(cmnd)
    type(command_line_type) :: cmnd
    integer :: verb 
    integer :: stat
    integer :: i
    character(clen) :: str

    character(clen), parameter :: myname="initialize_command_line"
    logical, parameter :: crash = .true.


    cmnd%nargs = command_argument_count()
    call get_command( command=cmnd%str, length=cmnd%len, status=stat )
    if (stat /= 0) call myerr("stat /= 0 after get_command", myname, crash)

    call mywrite('', verb) 
    call mywrite('command string: ' // trim(cmnd%str), verb)
    write(str,'(A,I5)') 'command nargs :', cmnd%nargs
    call mywrite( str, verb) 
    
    allocate( cmnd%args( 0:cmnd%nargs ), cmnd%arglens( 0:cmnd%nargs ) )

    do i = 0,cmnd%nargs
       call get_command_argument( i, &
            value=cmnd%args(i), &
            length=cmnd%arglens(i), &
            status=stat)

       if (stat /= 0) then
          call myerr("stat /= 0 after get_command_argument", myname, crash)
       end if

       write(str,'(A,I2,A,A)') 'arg ', i, ': ', trim(cmnd%args(i))
       call mywrite(str, verb)

       call mywrite('', verb)
    end do
    
  end function initialize_command_line


!-----------------------------------------------
!> reports an error message and stops execution
   subroutine myerr(str,routine,crash)
     character(*), intent(in) :: str
     character(*), intent(in) :: routine
     logical, intent(in) :: crash

     write(*,*) "*****************************"
     write(*,*) "error detected"
     write(*,*) "routine: ", trim(routine)
     write(*,*) trim(str)
     write(*,*) "*****************************"
     if (crash) stop

   end subroutine myerr


!-----------------------------------------------
!> verbosity dependent write
   subroutine mywrite( str, verb, lun, fmt, adv )
     character(*), intent(in) :: str
     integer, intent(in) :: verb
     integer, optional, intent(in) :: lun
     character(*), optional, intent(in) :: fmt
     logical, optional, intent(in) :: adv
     character(3) :: sadv

     character(clen), parameter :: myname='mywrite'
     logical, parameter :: crash=.true.

     if (present(adv)) then
        if (adv) then
           sadv = "yes"
        else
           sadv = "no "
        end if
     else
        sadv = "yes"
     end if

     if (present(adv) .and. .not. present(fmt) ) then
        call myerr('advance can only be used w/ formatted output', &
             myname, crash)
     end if

     if ( .not. present(lun) .and. .not. present(fmt) ) then
        if (verb <= myf90_verbosity) write(*,*) trim(str)
     else if ( present(lun) .and. .not. present(fmt) ) then 
        if (verb <= myf90_verbosity) write(lun,*) trim(str)
     else if ( .not. present(lun) .and. present(fmt) ) then 
        if (verb <= myf90_verbosity) write(*,fmt,advance=sadv) trim(str)
     else if ( present(lun) .and. present(fmt) ) then 
        if (verb <= myf90_verbosity) write(lun,fmt,advance=sadv) trim(str)
     end if

   end subroutine mywrite


!------------------------------------
!> returns a free logical unit number
   function get_free_lun() result(lun)

    integer(i4b) :: lun                        !< free lun
    integer(i4b) :: i                          !< loop counter
    integer(i4b), parameter :: minlun = 110    !< min lun to check
    integer(i4b), parameter :: maxlun = 1000   !< max lun to check
    logical :: badlun                          !< true if lun is already connected

    character(clen), parameter :: myname="get_free_lun"
 
    do i = minlun,maxlun
      inquire(unit=i, opened=badlun)
      if (.not. badlun) then
         lun = i
         return
      end if
    end do

    write(*,*) '  checked luns from, ', minlun, ' to ', maxlun

    call myerr( str=" no free logical unit numbers", &
                routine=myname, &
                crash=.true. )

  end function get_free_lun


!----------------------------------------
!> returns the current working directory
  function get_current_dir() result(path)
    character(clen) :: path  !< the current working dir

    call get_environment_variable('PWD',  value=path,  trim_name=.true.)

  end function get_current_dir


!--------------------------------------------------------
!> opens a formatted file for reading and returns the lun
  subroutine open_formatted_file_r(filename,lun)
  
    character(*), intent(in) :: filename !< name of file to open
    integer(i4b), intent(out) :: lun     !< lun 
    logical :: fthere

    character(clen), parameter :: myname="open_formatted_file_r"
    character(clen) :: str

    inquire(file=filename,exist=fthere)
    if (.not. fthere) then
       write(str,"(A,A)") "cant find file: ",  trim(filename)
       call myerr(str, myname, .true.)
    end if

    lun = get_free_lun() 
    open(unit=lun, file=filename, action="read")

  end subroutine open_formatted_file_r


!-------------------------------------------------------------
!> opens an unformatted file for reading and returns the lun
  subroutine open_unformatted_file_r(filename,lun)
  
    character(*), intent(in) :: filename  !< name of file to open
    integer(i4b), intent(out)  :: lun          !< lun
    logical :: fthere

    character(clen), parameter :: myname="open_unformatted_file_r"
    character(clen) :: str

    inquire(file=filename,exist=fthere)
    if (.not. fthere) then
       write(str,"(A,A)") "cant find file: ",  trim(filename)
       call myerr(str, myname, .true.)
    end if

    lun = get_free_lun()
    open(unit=lun, file=filename, form='unformatted', action="read")

  end subroutine open_unformatted_file_r

!--------------------------------------------------------------------
!> opens an unformatted file for stream reading and returns the lun
  subroutine open_unformatted_file_sr(filename,lun)
  
    character(*), intent(in) :: filename  !< name of file to open
    integer(i4b), intent(out)  :: lun          !< lun
    logical :: fthere

    character(clen), parameter :: myname="open_unformatted_file_sr"
    character(clen) :: str

    inquire(file=filename,exist=fthere)
    if (.not. fthere) then
       write(str,"(A,A)") "cant find file: ",  trim(filename)
       call myerr(str, myname, .true.)
    end if

    lun = get_free_lun()
    open(unit=lun, file=filename, form="unformatted", &
         action="read", access="stream", position="rewind")

  end subroutine open_unformatted_file_sr


!--------------------------------------------------------
!> opens a formatted file for writing and returns the lun
  subroutine open_formatted_file_w(filename,lun)
  
    character(*), intent(in) :: filename !< name of file to open
    integer(i4b), intent(out) :: lun     !< lun 

    lun = get_free_lun()
    open(unit=lun, file=filename, action="write")

  end subroutine open_formatted_file_w


!------------------------------------------------------------
!> opens an unformatted file for writing and returns the lun
  subroutine open_unformatted_file_w(filename,lun)
  
    character(*), intent(in) :: filename  !< name of file to open
    integer(i4b), intent(out)  :: lun             !< lun

    lun = get_free_lun()
    open(unit=lun, file=filename, form='unformatted', &
         action="write")
   
  end subroutine open_unformatted_file_w

!-------------------------------------------------------------------
!> opens an unformatted file for stream writing and returns the lun
  subroutine open_unformatted_file_sw(filename,lun)
  
    character(*), intent(in) :: filename  !< name of file to open
    integer(i4b), intent(out)  :: lun             !< lun

    lun = get_free_lun()
    open(unit=lun, file=filename, form='unformatted', &
         action="write", access="stream", position="rewind")

  end subroutine open_unformatted_file_sw



  ! scan file routines, one for each type
  !==================================================================
  subroutine scanfile_r4b(filename,keyword,var) 
    character(*), intent(in) :: filename
    character(*), intent(in) :: keyword
    real(r4b), intent(out) :: var

    integer :: keylen
    integer :: strlen
    integer :: keypos
    character(clen) :: string
    integer :: ierr
    integer :: lun
    integer :: count

    count=0
    keylen = len(adjustl(trim(keyword)))
    call open_formatted_file_r(filename,lun)
    do while( ierr==0 ) 
       read(lun,'(A)',iostat=ierr) string
       strlen = len(adjustl(trim(string)))
       keypos = index( string, adjustl(trim(keyword)) )
       if (keypos == 1) then
          count = count + 1
          read(string(keypos+keylen:strlen),*) var
       endif
    end do
    close(lun)

    if (count/=1) then
       write(*,'(A,A)') ' keyword ', trim(keyword)
       write(*,'(A,A)') ' file:   ', trim(filename)
       if (count==0) then
          write(*,*) ' *** keyword not found in file'
          stop
       else if (count>1) then
          write(*,*) ' *** keyword fount ', count, ' times'
          stop
       endif
    endif

  end subroutine scanfile_r4b


  subroutine scanfile_r8b(filename,keyword,var) 
    character(*), intent(in) :: filename
    character(*), intent(in) :: keyword
    real(r8b), intent(out) :: var

    integer :: keylen
    integer :: strlen
    integer :: keypos
    character(clen) :: string
    integer :: ierr
    integer :: lun
    integer :: count

    count=0
    keylen = len(adjustl(trim(keyword)))
    call open_formatted_file_r(filename,lun)
    do while( ierr==0 ) 
       read(lun,'(A)',iostat=ierr) string
       strlen = len(adjustl(trim(string)))
       keypos = index( string, adjustl(trim(keyword)) )
       if (keypos == 1) then
          count = count + 1
          read(string(keypos+keylen:strlen),*) var
       endif
    end do
    close(lun)

    if (count/=1) then
       write(*,'(A,A)') ' keyword ', trim(keyword)
       write(*,'(A,A)') ' file:   ', trim(filename)
       if (count==0) then
          write(*,*) ' *** keyword not found in file'
          stop
       else if (count>1) then
          write(*,*) ' *** keyword fount ', count, ' times'
          stop
       endif
    endif

  end subroutine scanfile_r8b


  subroutine scanfile_i4b(filename,keyword,var) 
    character(*), intent(in) :: filename
    character(*), intent(in) :: keyword
    integer(i4b), intent(out) :: var

    integer :: keylen
    integer :: strlen
    integer :: keypos
    character(clen) :: string
    integer :: ierr
    integer :: lun
    integer :: count

    count=0
    keylen = len(adjustl(trim(keyword)))
    call open_formatted_file_r(filename,lun)
    do while( ierr==0 ) 
       read(lun,'(A)',iostat=ierr) string
       strlen = len(adjustl(trim(string)))
       keypos = index( string, adjustl(trim(keyword)) )
       if (keypos == 1) then
          count = count + 1
          read(string(keypos+keylen:strlen),*) var
       endif
    end do
    close(lun)

    if (count/=1) then
       write(*,'(A,A)') ' keyword ', trim(keyword)
       write(*,'(A,A)') ' file:   ', trim(filename)
       if (count==0) then
          write(*,*) ' *** keyword not found in file'
          stop
       else if (count>1) then
          write(*,*) ' *** keyword fount ', count, ' times'
          stop
       endif
    endif

  end subroutine scanfile_i4b



  subroutine scanfile_i8b(filename,keyword,var) 
    character(*), intent(in) :: filename
    character(*), intent(in) :: keyword
    integer(i8b), intent(out) :: var

    integer :: keylen
    integer :: strlen
    integer :: keypos
    character(clen) :: string
    integer :: ierr
    integer :: lun
    integer :: count

    count=0
    keylen = len(adjustl(trim(keyword)))
    call open_formatted_file_r(filename,lun)
    do while( ierr==0 ) 
       read(lun,'(A)',iostat=ierr) string
       strlen = len(adjustl(trim(string)))
       keypos = index( string, adjustl(trim(keyword)) )
       if (keypos == 1) then
          count = count + 1
          read(string(keypos+keylen:strlen),*) var
       endif
    end do
    close(lun)

    if (count/=1) then
       write(*,'(A,A)') ' keyword ', trim(keyword)
       write(*,'(A,A)') ' file:   ', trim(filename)
       if (count==0) then
          write(*,*) ' *** keyword not found in file'
          stop
       else if (count>1) then
          write(*,*) ' *** keyword fount ', count, ' times'
          stop
       endif
    endif

  end subroutine scanfile_i8b



  subroutine scanfile_chr(filename,keyword,var)
    character(*), intent(in) :: filename
    character(*), intent(in) :: keyword
    character(clen), intent(out) :: var

    integer :: keylen
    integer :: strlen
    integer :: keypos
    character(clen) :: string
    integer :: ierr
    integer :: lun
    integer :: count

    count=0
    keylen = len(adjustl(trim(keyword)))
    call open_formatted_file_r(filename,lun)
    do while( ierr==0 ) 
       read(lun,'(A)',iostat=ierr) string
       strlen = len(adjustl(trim(string)))
       keypos = index( string, adjustl(trim(keyword)) )
       if (keypos == 1) then
          count = count + 1
          read(string(keypos+keylen:strlen),'(A)') var
          var = adjustl(trim(var))
       endif
    end do
    close(lun)

    if (count/=1) then
       write(*,'(A,A)') ' keyword ', trim(keyword)
       write(*,'(A,A)') ' file:   ', trim(filename)
       if (count==0) then
          write(*,*) ' *** keyword not found in file'
          stop
       else if (count>1) then
          write(*,*) ' *** keyword fount ', count, ' times'
          stop
       endif
    endif
  end subroutine scanfile_chr



  subroutine scanfile_log(filename,keyword,var)
    character(*), intent(in) :: filename
    character(*), intent(in) :: keyword
    logical, intent(out) :: var

    integer :: keylen
    integer :: strlen
    integer :: keypos
    character(clen) :: string
    integer :: ierr
    integer :: lun
    integer :: count

    count=0
    keylen = len(adjustl(trim(keyword)))
    call open_formatted_file_r(filename,lun)
    do while( ierr==0 ) 
       read(lun,'(A)',iostat=ierr) string
       strlen = len(adjustl(trim(string)))
       keypos = index( string, adjustl(trim(keyword)) )
       if (keypos == 1) then
          count = count + 1
          read(string(keypos+keylen:strlen),*) var
       endif
    end do
    close(lun)

    if (count/=1) then
       write(*,'(A,A)') ' keyword ', trim(keyword)
       write(*,'(A,A)') ' file:   ', trim(filename)
       if (count==0) then
          write(*,*) ' *** keyword not found in file'
          stop
       else if (count>1) then
          write(*,*) ' *** keyword fount ', count, ' times'
          stop
       endif
    endif

  end subroutine scanfile_log





end module myf90_mod


