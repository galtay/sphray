module hdf5_wrapper

  use hdf5

  implicit none
  save
  public

  ! Parameters that may be altered from within user programs:
  ! VERBOSITY, 0=Quiet, 1=Talkative
  integer, public :: HDF_VERBOSITY = 0

  ! 0=Don't suppress errors, 1=Suppress Errors
  integer, public :: HDF_ERROR_SUPPRESS = 1

  ! maximum number of files which can be open
  integer, private, parameter :: nfilemax = 1000

  ! table of currently open files
  integer, private :: nopen = 0
  integer(hid_t), dimension(nfilemax), private :: file_id
  logical, dimension(nfilemax), private        :: read_only
  character(len=1000), dimension(nfilemax)     :: fname_table

  ! is this the first call?
  logical, private :: needInit = .true.

  ! Data types for reading integers (there's no KIND=... support in HDF5!)
  integer(hid_t), private :: NATIVE_INTEGER4
  integer(hid_t), private :: NATIVE_INTEGER8

  ! Data types for reading reals (default REAL and DOUBLE PRECISION may be the 
  ! same, so we can't use those if the library is to be portable)
  integer(hid_t), private :: NATIVE_REAL4
  integer(hid_t), private :: NATIVE_REAL8

  ! Information about the currently open dataset
  logical,             private :: open_dset
  integer(hid_t),      private :: open_dset_id
  character(len=1000), private :: open_file_name, open_dset_name

  ! WHether to automatically open/close HDF5
  logical, private :: do_open = .true.

  ! Interface for generic read data routine
  interface hdf5_read_data
     module procedure read_scalar_real
     module procedure read_1d_real_array
     module procedure read_2d_real_array
     module procedure read_3d_real_array
     module procedure read_4d_real_array
     module procedure read_5d_real_array
     module procedure read_6d_real_array
     module procedure read_7d_real_array
     module procedure read_scalar_double
     module procedure read_1d_double_array
     module procedure read_2d_double_array
     module procedure read_3d_double_array
     module procedure read_4d_double_array
     module procedure read_5d_double_array
     module procedure read_6d_double_array
     module procedure read_7d_double_array
     module procedure read_scalar_integer4
     module procedure read_1d_integer4_array
     module procedure read_2d_integer4_array
     module procedure read_3d_integer4_array
     module procedure read_4d_integer4_array
     module procedure read_5d_integer4_array
     module procedure read_6d_integer4_array
     module procedure read_7d_integer4_array
     module procedure read_scalar_integer8
     module procedure read_1d_integer8_array
     module procedure read_2d_integer8_array
     module procedure read_3d_integer8_array
     module procedure read_4d_integer8_array
     module procedure read_5d_integer8_array
     module procedure read_6d_integer8_array
     module procedure read_7d_integer8_array
     module procedure read_scalar_string
     module procedure read_1d_string_array
     module procedure read_2d_string_array
     module procedure read_3d_string_array
     module procedure read_4d_string_array
     module procedure read_5d_string_array
     module procedure read_6d_string_array
     module procedure read_7d_string_array
  end interface

  ! interface for generic write_data routine
  interface hdf5_write_data
     module procedure write_1d_integer4_array
     module procedure write_2d_integer4_array
     module procedure write_3d_integer4_array
     module procedure write_4d_integer4_array
     module procedure write_5d_integer4_array
     module procedure write_6d_integer4_array
     module procedure write_7d_integer4_array
     module procedure write_scalar_integer4
     module procedure write_1d_integer8_array
     module procedure write_2d_integer8_array
     module procedure write_3d_integer8_array
     module procedure write_4d_integer8_array
     module procedure write_5d_integer8_array
     module procedure write_6d_integer8_array
     module procedure write_7d_integer8_array
     module procedure write_scalar_integer8
     module procedure write_1d_real_array
     module procedure write_2d_real_array
     module procedure write_3d_real_array
     module procedure write_4d_real_array
     module procedure write_5d_real_array
     module procedure write_6d_real_array
     module procedure write_7d_real_array
     module procedure write_scalar_real
     module procedure write_1d_double_array
     module procedure write_2d_double_array
     module procedure write_3d_double_array
     module procedure write_4d_double_array
     module procedure write_5d_double_array
     module procedure write_6d_double_array
     module procedure write_7d_double_array
     module procedure write_scalar_double
     module procedure write_1d_string_array
     module procedure write_2d_string_array
     module procedure write_3d_string_array
     module procedure write_4d_string_array
     module procedure write_5d_string_array
     module procedure write_6d_string_array
     module procedure write_7d_string_array
     module procedure write_scalar_string
  end interface

  ! interface for generic read_attribute routine
  interface hdf5_read_attribute
     module procedure read_1d_integer4_array_attribute
     module procedure read_2d_integer4_array_attribute
     module procedure read_3d_integer4_array_attribute
     module procedure read_4d_integer4_array_attribute
     module procedure read_5d_integer4_array_attribute
     module procedure read_6d_integer4_array_attribute
     module procedure read_7d_integer4_array_attribute
     module procedure read_scalar_integer4_attribute
     module procedure read_1d_integer8_array_attribute
     module procedure read_2d_integer8_array_attribute
     module procedure read_3d_integer8_array_attribute
     module procedure read_4d_integer8_array_attribute
     module procedure read_5d_integer8_array_attribute
     module procedure read_6d_integer8_array_attribute
     module procedure read_7d_integer8_array_attribute
     module procedure read_scalar_integer8_attribute
     module procedure read_1d_real_array_attribute
     module procedure read_2d_real_array_attribute
     module procedure read_3d_real_array_attribute
     module procedure read_4d_real_array_attribute
     module procedure read_5d_real_array_attribute
     module procedure read_6d_real_array_attribute
     module procedure read_7d_real_array_attribute
     module procedure read_scalar_real_attribute
     module procedure read_1d_double_array_attribute
     module procedure read_2d_double_array_attribute
     module procedure read_3d_double_array_attribute
     module procedure read_4d_double_array_attribute
     module procedure read_5d_double_array_attribute
     module procedure read_6d_double_array_attribute
     module procedure read_7d_double_array_attribute
     module procedure read_scalar_double_attribute
     module procedure read_1d_string_array_attribute
     module procedure read_2d_string_array_attribute
     module procedure read_3d_string_array_attribute
     module procedure read_4d_string_array_attribute
     module procedure read_5d_string_array_attribute
     module procedure read_6d_string_array_attribute
     module procedure read_7d_string_array_attribute
     module procedure read_scalar_string_attribute
  end interface

  interface hdf5_write_attribute
     module procedure write_1d_integer4_array_attribute
     module procedure write_2d_integer4_array_attribute
     module procedure write_3d_integer4_array_attribute
     module procedure write_4d_integer4_array_attribute
     module procedure write_5d_integer4_array_attribute
     module procedure write_6d_integer4_array_attribute
     module procedure write_7d_integer4_array_attribute
     module procedure write_scalar_integer4_attribute
     module procedure write_1d_integer8_array_attribute
     module procedure write_2d_integer8_array_attribute
     module procedure write_3d_integer8_array_attribute
     module procedure write_4d_integer8_array_attribute
     module procedure write_5d_integer8_array_attribute
     module procedure write_6d_integer8_array_attribute
     module procedure write_7d_integer8_array_attribute
     module procedure write_scalar_integer8_attribute
     module procedure write_1d_real_array_attribute
     module procedure write_2d_real_array_attribute
     module procedure write_3d_real_array_attribute
     module procedure write_4d_real_array_attribute
     module procedure write_5d_real_array_attribute
     module procedure write_6d_real_array_attribute
     module procedure write_7d_real_array_attribute
     module procedure write_scalar_real_attribute
     module procedure write_1d_double_array_attribute
     module procedure write_2d_double_array_attribute
     module procedure write_3d_double_array_attribute
     module procedure write_4d_double_array_attribute
     module procedure write_5d_double_array_attribute
     module procedure write_6d_double_array_attribute
     module procedure write_7d_double_array_attribute
     module procedure write_scalar_double_attribute
     module procedure write_1d_string_array_attribute
     module procedure write_2d_string_array_attribute
     module procedure write_3d_string_array_attribute
     module procedure write_4d_string_array_attribute
     module procedure write_5d_string_array_attribute
     module procedure write_6d_string_array_attribute
     module procedure write_7d_string_array_attribute
     module procedure write_scalar_string_attribute
  end interface

  public :: hdf5_auto_open

contains

  subroutine define_types
    !
    ! Create type IDs for 4 and 8 byte INTEGER and REAL data
    !
    implicit none
    integer :: i, hdf_err

    call byte_order(i)

    if(i.eq.1)then
       ! Little endian machine
       if(HDF_VERBOSITY.eq.1) THEN
          write(*,*)'[define_types] This machine is little endian.'
       ENDIF
       NATIVE_INTEGER8 = H5T_STD_I64LE
       NATIVE_INTEGER4 = H5T_STD_I32LE
       NATIVE_REAL4    = H5T_IEEE_F32LE
       NATIVE_REAL8    = H5T_IEEE_F64LE
    else
       ! Big endian machine
       if(HDF_VERBOSITY.eq.1) THEN
          write(*,*)'[define_types] This machine is big endian.'
       ENDIF
       NATIVE_INTEGER8 = H5T_STD_I64BE
       NATIVE_INTEGER4 = H5T_STD_I32BE
       NATIVE_REAL4    = H5T_IEEE_F32BE
       NATIVE_REAL8    = H5T_IEEE_F64BE
    end if

    return
  end subroutine define_types


  subroutine hdf5_auto_open(hdf5_open)
!
! Set whether to call h5open/h5close
!
    implicit none
    logical :: hdf5_open

    do_open = hdf5_open

    return
  end subroutine hdf5_auto_open


  subroutine hdf5_open_file(handle, fname, readonly)
    !
    ! open an existing file for reading or writing and return an integer
    ! file handle (a default integer, not an integer(hid_t)!)
    !  
    implicit none
    integer                      :: handle
    character(len=*), intent(in) :: fname
    integer(hid_t)               :: hdf_file_id
    integer                      :: hdf_err, hdf_err2
    integer                      :: i
    logical, optional            :: readonly
    logical                      :: need_write_access

    if(needInit) then
       needInit = .false.
       fname_table(1:nfilemax) = '*'
       if(do_open)call h5open_f(hdf_err)
       if (HDF_ERROR_SUPPRESS .eq. 0) then
          call h5eset_auto_f(1, hdf_err2)
       else
          !Verbosity = 0, supress HDF5 error messages
          call h5eset_auto_f(0, hdf_err2)
       endif
       if(hdf_err.lt.0)call hdf5_abort( &
            'Unable to initialize HDF5 to open file',fname=fname)
       file_id=-1
       call define_types()
    end if

    need_write_access = .true.
    if (present(readonly)) then
       if(readonly)then
          need_write_access = .false.
          call h5fopen_f (fname,H5F_ACC_RDONLY_F,hdf_file_id,hdf_err)
       else
          call h5fopen_f (fname,H5F_ACC_RDWR_F,hdf_file_id,hdf_err)
       end if
    else
       call h5fopen_f (fname,H5F_ACC_RDWR_F,hdf_file_id,hdf_err)
    endif

    if(hdf_err.lt.0.and.need_write_access)then
       ! Open in READ/WRITE mode has failed. See if read only works.
       call h5fopen_f(fname,H5F_ACC_RDONLY_F,hdf_file_id,hdf_err)
       ! If it did, the problem is most likely that we don't have write
       ! permission
       if(hdf_err.ge.0)then
          call hdf5_abort("Unable to open HDF5 file for writing"&
               //" in open_file() - check file permissions, or "&
               //"set readonly=.true. if you don't need write access", &
               fname=fname) 
       endif
    endif

    if(hdf_err.lt.0) call hdf5_abort( &
         'Unable to open HDF5 file in open_file()!',fname=fname)

    i=1
    do while(file_id(i).ge.0)
       i=i+1
       if(i.gt.nfilemax) call hdf5_abort( &
            'Too many files open in open_file()!',fname=fname)
    end do
    file_id(i)=hdf_file_id
    nopen=nopen+1

    handle=i

    read_only(handle)=.false.
    if(present(readonly))then
       if(readonly)then
          read_only(handle)=.true.
       endif
    endif

    fname_table(handle) = fname
    
    if(need_write_access)then
       if (HDF_VERBOSITY .ge. 1) then
          write(*,*)'[hdf5_open_file] File ',trim(fname),' Opened writeable'
       endif
    else
       if (HDF_VERBOSITY .ge. 1) then
          write(*,*)'[hdf5_open_file] File ',trim(fname),' Opened Read Only'
       endif
    endif

    return
  end subroutine hdf5_open_file



  subroutine hdf5_create_file(handle, fname)
    !
    ! Create a new file and return an integer
    ! file handle (a default integer, not an integer(hid_t)!)
    !  
    implicit none
    integer :: handle
    character(LEN=*), intent(IN) :: fname
    integer(HID_T) :: hdf_file_id
    integer :: hdf_err, hdf_err2
    integer :: i

    if(needInit)then
       needInit=.false.
       fname_table(1:nfilemax) = '*'
       if(do_open)call h5open_f(hdf_err)
       if (HDF_ERROR_SUPPRESS .gt. 0) then
          call h5eset_auto_f(1, hdf_err2)
       else
          !Verbosity = 0, supress HDF5 error messages
          call h5eset_auto_f(0, hdf_err2)
       endif
       if(hdf_err.lt.0) call hdf5_abort( &
            'Unable to initialize HDF5 in create_file()!',fname=fname)
       file_id=-1
       call define_types
    end if
    call h5fcreate_f(fname, H5F_ACC_TRUNC_F, hdf_file_id, hdf_err, H5P_DEFAULT_F, H5P_DEFAULT_F)

    if(hdf_err.lt.0) call hdf5_abort( &
         'Unable to open HDF5 file in create_file()!',fname=fname)

    if (HDF_VERBOSITY .ge. 1) then
       write(*,*)'[hdf5_create_file] File ',trim(fname),' Created'
    endif

    i=1
    do while(file_id(i).ge.0)
       i=i+1
       if(i.gt.nfilemax) call hdf5_abort( &
            'Too many files open in create_file()!',fname=fname)
    end do
    file_id(i)=hdf_file_id
    nopen=nopen+1
    handle=i

    read_only(handle) = .false.
    fname_table(handle) = fname

    return
  end subroutine hdf5_create_file


  subroutine hdf5_close_file(ifile)
    !
    ! Close a currently open HDF5 file
    !  
    implicit none
    integer, intent(IN) :: ifile
    integer :: hdf_err

    if(ifile.lt.1.or.ifile.gt.nfilemax)call hdf5_abort( &
         'Invalid file handle in close_file')

    if(file_id(ifile).lt.0)call hdf5_abort( &
         'File is not open in close_file()!')

    if (HDF_VERBOSITY .ge. 1) then
       write(*,*)'[hdf5_close_file] Closed file '
    endif

    ! If we left a dataset open, close it
    if(open_dset.and.fname_table(ifile).eq.open_file_name)then
       call h5dclose_f(open_dset_id,   hdf_err)
       open_dset = .false.
    endif

    call h5fclose_f(file_id(ifile), hdf_err) 

    if(hdf_err.lt.0)call hdf5_abort( &
         'Unable to close HDF5 file in close_file()!',fname=fname_table(ifile))

    nopen=nopen-1
    file_id(ifile)=-1

    if(nopen.eq.0) then 
       if(do_open)call h5close_f(hdf_err)
       needInit = .true.
    endif

    fname_table(ifile) = '*'

    return
  end subroutine hdf5_close_file


  subroutine hdf5_get_dimensions(ifile,name,rank,dims)
!
! Routine to find out the rank and dimensions of a dataset
! or attribute.
!
    implicit none
    integer, intent(IN) :: ifile
    character(len=*), intent(IN) :: name
    integer, intent(OUT) :: rank
    integer, dimension(*), intent(out) :: dims

    integer(hid_t) :: dset_id, dspace_id, attr_id, loc_id
    integer :: hdf_err, hdf_err2
    integer(hsize_t), dimension(10) :: hdf_dims,hdf_maxdims
    logical :: isattribute, isgroup
    integer :: nslash
    character(len=200) :: loc_name,attr_name

    if(ifile.lt.1.or.ifile.gt.nfilemax)call hdf5_abort( &
         'Invalid file handle in get_dimensions',name=name)

    if(file_id(ifile).lt.0) call hdf5_abort( &
         'File is not open in get_dimensions()!', &
         name=name)

    ! Try to open 'name' as a dataset then as an attribute
    call h5eset_auto_f(0, hdf_err2)
    call h5dopen_f(file_id(ifile),name,dset_id,hdf_err)
    if(HDF_ERROR_SUPPRESS .eq. 0)call h5eset_auto_f(1, hdf_err2)

    if(hdf_err.lt.0)then
       ! Its not a dataset - try it as an attribute. Need to identify
       ! what it is an attribute of!
       ! Split 'name' into parent object name and attribute name
       nslash=index(name,"/",.true.)
       if(nslash.eq.0)call hdf5_abort( &
            'Invalid dataset/attribute name in get_dimensions()!', &
            fname=fname_table(ifile),name=name)
       loc_name=name(1:nslash-1)
       attr_name=name(nslash+1:len_trim(name))

       ! Try to open loc_name as a group then as a dataset
       isgroup=.false.
       call h5eset_auto_f(0, hdf_err2)
       call h5dopen_f(file_id(ifile),loc_name,loc_id,hdf_err)
       if(HDF_ERROR_SUPPRESS .eq. 0)call h5eset_auto_f(1, hdf_err2)
       if(hdf_err.lt.0)then
          call h5gopen_f(file_id(ifile),loc_name, loc_id, hdf_err)
          if(hdf_err.lt.0)call hdf5_abort( &
               'Unable to open attribute parent object in get_dimensions()!', &
               name=name,fname=fname_table(ifile))
          isgroup=.true.
       end if

       ! Open the attribute and get its dataspace
       call h5aopen_name_f(loc_id,attr_name,attr_id,hdf_err) 
       if(hdf_err.lt.0)call hdf5_abort( &
            'Cannot open datset/attribute in get_dimensions()!', &
            name=name,fname=fname_table(ifile))
       call h5aget_space_f(attr_id, dspace_id, hdf_err)
       isattribute=.true.

    else
       call h5dget_space_f(dset_id, dspace_id, hdf_err)
       isattribute=.false.
    end if

    call h5sget_simple_extent_ndims_f(dspace_id,rank,hdf_err)
    if(rank.gt.0) then
       call h5sget_simple_extent_dims_f(dspace_id,hdf_dims,hdf_maxdims, &
            hdf_err)
    else
       hdf_dims=0
       hdf_maxdims=0
    endif
    call h5sclose_f(dspace_id,hdf_err)

    if(isattribute)then
       call h5aclose_f(attr_id,hdf_err)
       if(isgroup)then
          call h5gclose_f(loc_id,hdf_err)
       else
          call h5dclose_f(loc_id,hdf_err)
       end if
    else
       call h5dclose_f(dset_id,hdf_err)
    end if

    ! Convert hdf_dims from hsize_t to default integer
    dims(1:rank)=INT(hdf_dims(1:rank),KIND(dims))

    return
  end subroutine hdf5_get_dimensions


  subroutine hdf5_get_type(ifile,name,datatype,size)
!
! Routine to find out the type of a dataset
! or attribute.
!
    implicit none
    integer, intent(IN) :: ifile
    character(len=*), intent(IN) :: name
    character(len=*), intent(out) :: datatype
    integer :: class

    integer(hid_t) :: dset_id, dspace_id, attr_id, loc_id, dtype_id
    integer :: hdf_err, hdf_err2
    integer(hsize_t), dimension(10) :: hdf_dims,hdf_maxdims
    logical :: isattribute, isgroup
    integer :: nslash
    character(len=200) :: loc_name,attr_name
    integer, intent(out) :: size
    integer(kind=size_t) :: h5size

    if(ifile.lt.1.or.ifile.gt.nfilemax)call hdf5_abort( &
         'Invalid file handle in get_type',name=name)

    if(file_id(ifile).lt.0) call hdf5_abort( &
         'File is not open in get_type()!', &
         name=name)

    ! Try to open 'name' as a dataset then as an attribute
    call h5eset_auto_f(0, hdf_err2)
    call h5dopen_f(file_id(ifile),name,dset_id,hdf_err)
    if(HDF_ERROR_SUPPRESS .eq. 0)call h5eset_auto_f(1, hdf_err2)

    if(hdf_err.lt.0)then
       ! Its not a dataset - try it as an attribute. Need to identify
       ! what it is an attribute of!
       ! Split 'name' into parent object name and attribute name
       nslash=index(name,"/",.true.)
       if(nslash.eq.0)call hdf5_abort( &
            'Invalid dataset/attribute name in get_type()!', &
            fname=fname_table(ifile),name=name)
       loc_name=name(1:nslash-1)
       attr_name=name(nslash+1:len_trim(name))

       ! Try to open loc_name as a group then as a dataset
       isgroup=.false.
       call h5eset_auto_f(0, hdf_err2)
       call h5dopen_f(file_id(ifile),loc_name,loc_id,hdf_err)
       if(HDF_ERROR_SUPPRESS .eq. 0)call h5eset_auto_f(1, hdf_err2)
       if(hdf_err.lt.0)then
          call h5gopen_f(file_id(ifile),loc_name, loc_id, hdf_err)
          if(hdf_err.lt.0)call hdf5_abort( &
               'Unable to open attribute parent object in get_type()!', &
               name=name,fname=fname_table(ifile))
          isgroup=.true.
       end if

       ! Open the attribute and get its datatype
       call h5aopen_name_f(loc_id,attr_name,attr_id,hdf_err) 
       if(hdf_err.lt.0)call hdf5_abort( &
            'Cannot open datset/attribute in get_type()!', &
            name=name,fname=fname_table(ifile))
       call h5aget_type_f(attr_id, dtype_id, hdf_err) 
       isattribute=.true.

    else
       call h5dget_type_f(dset_id, dtype_id, hdf_err) 
       isattribute=.false.
    end if

    ! Determine type of the data
    call h5tget_class_f(dtype_id, class, hdf_err) 
    if(class.eq.H5T_INTEGER_F) then
       datatype = "INTEGER"
    else if(class.eq.H5T_FLOAT_F) then
       datatype = "REAL"
    else if(class.eq.H5T_STRING_F) then
       datatype = "CHARACTER"
    else
       datatype = "OTHER"
    end if
                
    ! Find the size of this data type
    call h5tget_size_f(dtype_id, h5size, hdf_err) 
    size = h5size

    call h5tclose_f(dtype_id,hdf_err)

    if(isattribute)then
       call h5aclose_f(attr_id,hdf_err)
       if(isgroup)then
          call h5gclose_f(loc_id,hdf_err)
       else
          call h5dclose_f(loc_id,hdf_err)
       end if
    else
       call h5dclose_f(dset_id,hdf_err)
    end if

    return
  end subroutine hdf5_get_type

!
! Stop cleanly and print an error message
!
  subroutine hdf5_abort(message,name,fname)
    implicit none
    character(len=*)           :: message
    character(len=*), optional :: name, fname
    integer                    :: hdf_err
    integer, parameter :: ncol = 80

    if(do_open)call h5close_f(hdf_err)

    write(*,*)' '
    call write_wrapped_text('***Abort HDF5 : '//trim(message), &
         ncol, indent=16)
    write(*,*)' '
    if (present(name)) then
       call write_wrapped_text('   identifier : '//trim(name), &
            ncol)
    endif
    if (present(fname)) then
       call write_wrapped_text('   file name  : '//trim(fname),&
            ncol)
    end if
    write(*,*)' '
    stop
  end subroutine hdf5_abort


!
! Function for creating groups:
!
  subroutine hdf5_create_group(ifile,name)
    implicit none
    integer          :: ifile
    character(len=*) :: name
    character, parameter            :: slash='/'
    integer, parameter              :: MAX_NEST = 10
    integer(hid_t)                  :: id_chain(MAX_NEST)
    integer(hid_t)                  :: group_id, loc_id
    integer, parameter              :: LEN_STR = 256

    character(len=LEN_STR)          :: thestring,currentgroup
    integer          :: i,n_nest, starti, endi, hdf_err

    if(ifile.lt.1.or.ifile.gt.nfilemax)call hdf5_abort( &
         'Invalid file handle in create_group',name=name)

    if(read_only(ifile))call hdf5_abort( &
         "Attempt to create group in read only file!", &
         name=name,fname=fname_table(ifile));

    thestring = name
    i = index(thestring,slash)
    if (i .eq. 1) then
       !strip leading slash
       thestring = thestring(2:LEN_STR)
       i = index(thestring,slash)
    endif

    !Ensure there is a trailing / on the end:
    i = index(thestring,slash,.true.)
    if (i .lt. len_trim(thestring)) then
       thestring(len_trim(thestring)+1:len_trim(thestring)+1) = '/'
    endif
    i = index(thestring,slash)

    loc_id = file_id(ifile)

    starti = 1
    endi   = i
    
    n_nest = 0
    do while(i > 0)
       !Get the name of the current group:
       currentgroup = thestring(1:i-1)
       
       !Cut the current group from the start of the string:
       thestring = thestring(i+1:LEN_STR)
       
       !Update position of the first /
       i = index(thestring,slash)
       
       !Create groups:
       call h5eset_auto_f(0, hdf_err)
       call h5gcreate_f(loc_id, currentgroup, group_id, hdf_err, 0_size_t)
       if (hdf_err .eq. 0 .and. HDF_VERBOSITY .ge. 1) then
          write(*,*)	'[hdf5_create_group] Creating group ',trim(currentgroup)
       endif
       if(HDF_ERROR_SUPPRESS.eq.0)call h5eset_auto_f(1, hdf_err)
       
       call h5gopen_f(loc_id, currentgroup, group_id, hdf_err)
       if(hdf_err.lt.0) then
          call hdf5_abort('Unable to open group in create_group()!', &
               name=name,fname=fname_table(ifile))
       endif
       
       n_nest = n_nest + 1
       id_chain(n_nest) = group_id
       loc_id = group_id
    enddo
    
    ! Close the nest of groups:
    do i=n_nest,1,-1
       call h5gclose_f(id_chain(i),hdf_err)
       if(hdf_err .lt. 0) call hdf5_abort( &
            'Unable to close group in write_data()!', &
            name=name,fname=fname_table(ifile))
    enddo
    
  end subroutine hdf5_create_group

!
! Subroutine to generate a list of the datasets in a group
!
  subroutine hdf5_list_datasets(ifile, name, ndatasets, dataset_names)

    implicit none
    ! Parameters
    integer, intent(in)  :: ifile
    integer, intent(out) :: ndatasets
    character(len=*), intent(in)                          :: name
    character(len=*), dimension(:), intent(out), optional :: dataset_names
    ! Internal
    integer(hid_t) :: group_id
    integer :: hdf_err, nmembers, obj_type
    integer :: i, nd
    character(len=500) :: obj_name

    ! Check file handle is ok
    if(ifile.lt.1.or.ifile.gt.nfilemax)call hdf5_abort( &
         'Invalid file handle in list_datasets',name=name)
    if(file_id(ifile).lt.0) call hdf5_abort( &
         'File is not open in list_datasets()!', &
         name=name)
    
    ! Try to open the group
    call h5gopen_f(file_id(ifile), name, group_id, hdf_err)
    if(hdf_err .lt. 0) call hdf5_abort( &
         'Unable to open group in list_datasets()!', &
         name=name,fname=fname_table(ifile))

    ! Determine number of objects in the group
    call h5gn_members_f(file_id(ifile), name, nmembers, hdf_err)
    if(hdf_err.ne.0) call hdf5_abort( &
         'Unable to retrieve object count in list_datasets()!', &
         name=name,fname=fname_table(ifile))

    ! Find out which objects are datasets
    nd = 0
    do i = 1, nmembers, 1
       call h5gget_obj_info_idx_f(file_id(ifile), name, i-1, & 
            obj_name, obj_type, hdf_err)
       if(hdf_err.ne.0) call hdf5_abort( &
            'Unable to retrieve object info in list_datasets()!', &
            name=name,fname=fname_table(ifile))
       if(obj_type.eq.H5G_DATASET_F)then
          nd = nd + 1
          if(present(dataset_names))then
             dataset_names(nd) = trim(obj_name)
          endif
       endif
    end do
    ndatasets = nd

    ! Close the group
    call h5gclose_f(group_id, hdf_err)
    
    return
  end subroutine hdf5_list_datasets

!
! Subroutine to generate a list of the groups inside another group
!
  subroutine hdf5_list_groups(ifile, name, ngroups, group_names)

    implicit none
    ! Parameters
    integer, intent(in)  :: ifile
    integer, intent(out) :: ngroups
    character(len=*), intent(in)                          :: name
    character(len=*), dimension(:), intent(out), optional :: group_names
    ! Internal
    integer(hid_t) :: group_id
    integer :: hdf_err, nmembers, obj_type
    integer :: i, nd
    character(len=500) :: obj_name

    ! Check file handle is ok
    if(ifile.lt.1.or.ifile.gt.nfilemax)call hdf5_abort( &
         'Invalid file handle in list_groups',name=name)
    if(file_id(ifile).lt.0) call hdf5_abort( &
         'File is not open in list_groups()!', &
         name=name)
    
    ! Try to open the group
    call h5gopen_f(file_id(ifile), name, group_id, hdf_err)
    if(hdf_err .lt. 0) call hdf5_abort( &
         'Unable to open group in list_groups()!', &
         name=name,fname=fname_table(ifile))

    ! Determine number of objects in the group
    call h5gn_members_f(file_id(ifile), name, nmembers, hdf_err)
    if(hdf_err.ne.0) call hdf5_abort( &
         'Unable to retrieve object count in list_groups()!', &
         name=name,fname=fname_table(ifile))

    ! Find out which objects are groups
    nd = 0
    do i = 1, nmembers, 1
       call h5gget_obj_info_idx_f(file_id(ifile), name, i-1, & 
            obj_name, obj_type, hdf_err)
       if(hdf_err.ne.0) call hdf5_abort( &
            'Unable to retrieve object info in list_groups()!', &
            name=name,fname=fname_table(ifile))
       if(obj_type.eq.H5G_GROUP_F)then
          nd = nd + 1
          if(present(group_names))then
             group_names(nd) = trim(obj_name)
          endif
       endif
    end do
    ngroups = nd

    ! Close the group
    call h5gclose_f(group_id, hdf_err)
    
    return
  end subroutine hdf5_list_groups

!
! Subroutine to generate a list of the attributes attached to a group or dataset
!
  subroutine hdf5_list_attributes(ifile, name, nattribs, attrib_names)

    implicit none
    ! Parameters
    integer, intent(in)  :: ifile
    integer, intent(out) :: nattribs
    character(len=*), intent(in)                          :: name
    character(len=*), dimension(:), intent(out), optional :: attrib_names
    ! Internal
    integer(hid_t) :: loc_id, attr_id
    integer :: hdf_err, nmembers, obj_type
    integer :: i, nd
    character(len=500) :: obj_name
    logical :: is_group

    ! Check file handle is ok
    if(ifile.lt.1.or.ifile.gt.nfilemax)call hdf5_abort( &
         'Invalid file handle in list_attributes',name=name)
    if(file_id(ifile).lt.0) call hdf5_abort( &
         'File is not open in list_attributes()!', &
         name=name)
    
    ! Try to open the location as a group
    call h5eset_auto_f(0, hdf_err)
    is_group = .true.
    call h5gopen_f(file_id(ifile), name, loc_id, hdf_err)
    if(hdf_err.ne.0)then
       ! and if that fails, try it as a dataset
       call h5dopen_f(file_id(ifile), name, loc_id, hdf_err)
       is_group = .false.
    endif
    if(HDF_ERROR_SUPPRESS .eq. 0)call h5eset_auto_f(1, hdf_err)
    if(hdf_err .lt. 0) call hdf5_abort( &
         'Unable to open parent object in list_attributes()!', &
         name=name,fname=fname_table(ifile))
    
    ! Determine number of attributes
    call h5aget_num_attrs_f(loc_id, nattribs, hdf_err)
    if(hdf_err.ne.0) call hdf5_abort( &
         'Unable to retrieve attribute count in list_attributes()!', &
         name=name,fname=fname_table(ifile))

    ! Get the names of the attributes
    if(present(attrib_names))then
       do i = 1, nattribs, 1
          call h5aopen_idx_f(loc_id, i-1, attr_id, hdf_err)
          if(hdf_err.ne.0) call hdf5_abort( &
               'Unable to open attribute in list_attributes()!', &
               name=name,fname=fname_table(ifile))
          call h5aget_name_f(attr_id, int(len(attrib_names(i)),kind=size_t),&
               attrib_names(i), hdf_err)
          call h5aclose_f(attr_id, hdf_err)
       end do
    endif

    ! Close the group or dataset
    if(is_group)then
       call h5gclose_f(loc_id, hdf_err)
    else
       call h5dclose_f(loc_id, hdf_err)
    endif

    return
  end subroutine hdf5_list_attributes

!
! This function returns a HDF5 file ID given a wrapper file handle
! (useful if you want to use the wrapper and direct HDF5 calls
! to the same file)
!
  integer(kind=hid_t) function hdf5_get_file_id(ifile)
    
    implicit none
    integer, intent(in) :: ifile
    
    ! Check file handle is ok
    if(ifile.lt.1.or.ifile.gt.nfilemax)call hdf5_abort( &
         'Invalid file handle in get_file_id')
    if(file_id(ifile).lt.0) call hdf5_abort( &
         'File is not open in get_file_id()')
    
    ! Return the file ID
    hdf5_get_file_id = file_id(ifile)
    
    return
  end function hdf5_get_file_id


!
! Reading and writing subroutines
!
#include "read_dataset_preprocessor.inc"

#include "read_attribute_preprocessor.inc"

#include "write_dataset_preprocessor.inc"
 
#include "write_attribute_preprocessor.inc"



  subroutine write_wrapped_text(str, n, iunit, indent)

    ! Writes the specified text with lines wrapped at the
    ! specified maximum length (handy for long error
    ! messages!)

    implicit none
    character(len=*), intent(in) :: str
    integer, intent(in) :: n
    integer, intent(in), optional :: iunit, indent
    integer :: i, j, nchar
    integer :: nspace
    character(len=len_trim(str)+1) :: buf

    nspace = 0
    buf = trim(adjustl(str))

    i = 1
    do while(i.le.len_trim(buf))

       if(present(indent).and.i.gt.1)nspace = indent

       j = min(i+n-2-nspace,len(buf))
       nchar = scan(buf(i:j), " ", back=.true.)

       ! If we can't split the line at a space, just split after
       ! n characters (taking any indentation into account)
       if(nchar.lt.1)nchar = min(n-2-nspace,len(buf)-i)

       if(.not.present(iunit))then
          write(*,*)repeat(" ",nspace)//trim(adjustl(buf(i:i+nchar-1)))
       else
          write(iunit,*)repeat(" ",nspace)//trim(adjustl(buf(i:i+nchar-1)))
       endif
       i = i + nchar
    end do

    return
  end subroutine write_wrapped_text


end module hdf5_wrapper
