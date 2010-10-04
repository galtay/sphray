  subroutine SUB_NAME (ifile,name,arr,start,count)
!
! Read an n-D array dataset from an open file. This assumes some
! preprocessor variables have been set - see hdf_wrapper.F90.
!
! Try to read all data first as a dataset, and failing that as an
! attribute
!
    implicit none
    ! Parameters
    character(LEN=*), intent(IN) :: name 
    integer, dimension(:), intent(IN), optional :: start,count
    integer,intent(IN) :: ifile
#ifdef SCALAR
    ARR_TYPE :: arr
#else
    ARR_TYPE , dimension ARRAY_DIM :: arr
#endif
    ! HDF5 error status
    integer               :: hdf_err
    ! HDF5 identifiers
    integer(kind=hid_t)   :: memtype_id, dset_id, dtype_id
    ! Information about the type of data in the file
    integer(kind=size_t)  :: dtype_size
    integer               :: dtype_class
    ! Size of one element in the data array
    integer(kind=size_t)  :: nativesize
    ! Information about the data in the file
    integer(kind=hid_t)   :: file_dspace_id
    integer               :: file_rank
    integer(kind=hsize_t) :: file_dims(7), file_maxdims(7)
    ! Information about the dimensions of the data to read
    integer :: eff_read_rank
    integer(kind=hsize_t) :: read_dims(7), eff_read_dims(7), read_dim_idx(7)
    ! Information about the dimensions of the array
    integer               :: array_rank, eff_array_rank
    integer(kind=hsize_t) :: array_dims(7),  eff_array_dims(7), array_dim_idx(7)
    ! Dimensions of the memory dataspace
    integer(kind=hsize_t) :: memspace_dims(7), maxdims(7)
    ! Memory dataspace handle
    integer(kind=hid_t)   :: memspace_id
    ! Arrays for select hyperslab
    integer(kind=hsize_t) :: hs_start(7), hs_count(7)

    ! Loop indexes etc
    integer :: i, j

    ! Check inputs
    if(present(start).and.(.not.present(count)))call hdf5_abort( &
         'If parameter start is present count must be present in read_dataset')
    if(present(count).and.(.not.present(start)))call hdf5_abort( &
         'If parameter count is present start must be present in read_dataset')
    if(ifile.lt.1.or.ifile.gt.nfilemax)call hdf5_abort( &
         'Invalid file handle in read_data',name=name)
    if(file_id(ifile).lt.0) &
         call hdf5_abort('File is not open in read_data()!', &
         name=name)

    ! Determine datatype for array in memory
#ifndef STRING
    call h5tcopy_f( NATIVE_TYPE, memtype_id, hdf_err)
    if(hdf_err.lt.0)call hdf5_abort( &
         'Unable to open datatype in read_data()!', &
         name=name,fname=fname_table(ifile))
#else
    call h5tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, hdf_err)
    call h5tset_size_f(memtype_id, INT(len(arr),KIND=SIZE_T), hdf_err) 
#endif

    ! Open datset and get size. Dataset may already be open
    if(trim(adjustl(fname_table(ifile))).ne.trim(adjustl(open_file_name)) &
         .or.trim(adjustl(name)).ne.trim(adjustl(open_dset_name)))then
       if(open_dset)then
          open_dset = .false.
          call h5dclose_f(open_dset_id, hdf_err)
       endif
    endif
    if(open_dset)then
       dset_id = open_dset_id
    else
       call h5dopen_f(file_id(ifile),name, dset_id, hdf_err)
       open_dset_id = dset_id
       open_dset = .true.
       open_dset_name = trim(name)
       open_file_name = trim(fname_table(ifile))
    endif

    if (hdf_err .lt. 0) then
       call hdf5_abort('Unable to open dataset in read_dataset', &
            name=name,fname=fname_table(ifile))
    endif
    call h5dget_type_f(dset_id, dtype_id, hdf_err)
    call h5tget_size_f(dtype_id, dtype_size, hdf_err)
    call h5tget_class_f(dtype_id, dtype_class, hdf_err) 
    call h5tget_size_f(memtype_id, nativesize, hdf_err) 
    call h5dget_space_f(dset_id, file_dspace_id, hdf_err) 
    call h5sget_simple_extent_ndims_f(file_dspace_id, file_rank, hdf_err) 
    call h5sget_simple_extent_dims_f(file_dspace_id, file_dims, file_maxdims, hdf_err) 

    ! Check we can read this type of data into this type of array
    if(dtype_class.ne. DATA_CLASS ) then
       call hdf5_abort('Dataset is incompatible type in read_data()', &
            name=name,fname=fname_table(ifile))
    endif
#ifndef STRING
    if(dtype_size.gt.nativesize) then
       call hdf5_abort &
            ('Precision in file is greater than precision of array in read_data()', &
            name=name,fname=fname_table(ifile))
    endif
#else
    if(dtype_size.gt.nativesize) then
       call hdf5_abort &
            ('String in file doesn''t fit in string in memory in read_data()',&
            name=name,fname=fname_table(ifile))
    endif
#endif

    ! Make sure start and count have enough dimensions if present
    if(present(start))then
       if(size(start).lt.file_rank) call hdf5_abort( &
            'start parameter has incorrect dimensions in read_data', &
            name=name,fname=fname_table(ifile))
       if(size(count).lt.file_rank) call hdf5_abort( &
            'count parameter has incorrect dimensions in read_data', &
            name=name,fname=fname_table(ifile))
    endif

    ! Get the dimensions of the requested data
    if(present(count))then
       ! count parameter gives the size of the data to be read
       read_dims(1:file_rank) = count(1:file_rank)
    else
       ! otherwise we want to read the whole dataset
       read_dims(1:file_rank) = file_dims(1:file_rank)
    endif
    ! Throw away any of these dimensions with size 1
    eff_read_rank = 0
    do i = 1, file_rank, 1
       if(read_dims(i).gt.1)then
          eff_read_rank = eff_read_rank + 1
          eff_read_dims(eff_read_rank) = read_dims(i)
          ! For each effective dimension this gives the index of the 'real'
          ! dimension in the dataset
          read_dim_idx(eff_read_rank) = i
       endif
    end do

    ! Repeat this for the array dimensions
    array_rank          = NDIMS

#ifndef SCALAR
    do i = 1, NDIMS, 1
       array_dims(i) = size(arr,i)
    end do
#endif

    eff_array_rank = 0
    do i = 1, array_rank, 1
       if(array_dims(i).gt.1)then
          eff_array_rank = eff_array_rank + 1
          eff_array_dims(eff_array_rank) = array_dims(i)
          ! For each effective dimension this gives the index of the 'real'
          ! dimension in the array
          array_dim_idx(eff_array_rank) = i
       endif
    end do

    ! All dimensions with size>1 have to be at least as big
    ! in the array as in the file dataset, otherwise the data
    ! isn't going to fit.
    if(eff_array_rank.ne.eff_read_rank)then
       write(0,*)'Array rank (omitting dimensions with size=1) : ',eff_array_rank
       write(0,*)'Data  rank (omitting dimensions with size=1) : ',eff_read_rank
       call hdf5_abort( &
            'Supplied array does not have the same rank as the data to be read', &
            name=name,fname=fname_table(ifile))
    endif
    ! Dimensions of the array need to be big enough in each direction
    do i = 1, eff_array_rank, 1
       if(eff_array_dims(i).lt.eff_read_dims(i))call hdf5_abort( &
            'Data to be read does not fit in the supplied array', &
            name=name,fname=fname_table(ifile))
    end do

    ! Now we need to create a dataspace to describe the array in memory,
    ! then select the region that the data to be read in will occupy.
    ! First, figure out the dimensions of this data space. It has to
    ! have the same rank as the dataset in the file so we can't just use
    ! the dimensions of the array (e.g. may want to read a 2D slice of
    ! a 3D dataset, or read a single value from a 1D dataset into a scalar).
    memspace_dims(1:7) = 1
    do i = 1, eff_array_rank, 1
       j = read_dim_idx(i)
       memspace_dims(j) = eff_array_dims(i)
    end do

    ! Create the memory dataspace
    maxdims(1:file_rank) = memspace_dims(1:file_rank)
    call h5screate_simple_f(file_rank,memspace_dims,memspace_id,hdf_err,maxdims) 
    if(hdf_err.lt.0)call hdf5_abort( &
         'H5Screate_simple() call failed (possible wrapper bug?)',&
         name=name,fname=fname_table(ifile))

    ! Select the part of this dataspace that will hold the data
    if(eff_array_rank.gt.0)then
       hs_start(1:file_rank) = 0
       hs_count(1:file_rank) = read_dims(1:file_rank)
       call h5sselect_hyperslab_f(memspace_id, H5S_SELECT_SET_F, &
            hs_start, hs_count, hdf_err)
       if(hdf_err.lt.0)call hdf5_abort( &
            'H5Sselect_hyperslab call #1 failed (possible wrapper bug?)',&
            name=name,fname=fname_table(ifile))
    endif

    ! Select the part of the file dataspace to be read
    if(present(count))then
       hs_start(1:file_rank) = start(1:file_rank)-1
       hs_count(1:file_rank) = count(1:file_rank)
       call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, &
            hs_start, hs_count, hdf_err)
       if(hdf_err.lt.0)call hdf5_abort( &
            'H5Sselect_hyperslab() call #2 failed (possible wrapper bug?)',&
            name=name,fname=fname_table(ifile))
    endif

    ! Read the dataset
    call read_hdf5_dataset(dset_id, memtype_id , memspace_id, &
         file_dspace_id, arr, hdf_err)

    if(hdf_err.lt.0) then
       call hdf5_abort('Unable to read dataset in read_data()', &
            name=name,fname=fname_table(ifile))
    endif

    if (HDF_VERBOSITY .ge. 1) then
       write(*,*)'[read_dataset] Read dataset ',trim(name) 
    endif

    ! Close everything, except the dataset which will be left
    ! open in case there's another read to the same file
    call h5tclose_f(dtype_id,hdf_err)
    call h5tclose_f(memtype_id,hdf_err)
    call h5sclose_f(file_dspace_id,hdf_err)
    call h5sclose_f(memspace_id,hdf_err)    

    return
  end subroutine SUB_NAME
