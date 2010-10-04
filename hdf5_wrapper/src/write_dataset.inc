subroutine SUB_NAME (ifile,name,arr,gzip,overwrite,szip,extensible,start,count, &
     initial_size)
  !
  ! Write an n-D array dataset to an open file. This assumes some
  ! preprocessor variables have been set - see hdf_wrapper.F90.
  !
  ! Compress data if 'gzip' parameter is set - value is compression level, 0-9,
  ! as used by the unix gzip command.
  !
  implicit none
  ! Parameters
  character(len=*), intent(in)    :: name 
  integer, intent(in)             :: ifile
  integer, intent(in), optional   :: gzip
  logical, intent(in), optional   :: szip
  logical, intent(in), optional   :: extensible
#ifndef SCALAR
  ARR_TYPE , dimension ARRAY_DIM  :: arr
#else
  ARR_TYPE                        :: arr
#endif
  logical, optional               :: overwrite
  integer, dimension(:), optional :: start, count, initial_size
  ! Loop indices etc
  integer :: i, j
  ! Gzip/szip
  logical :: gzip_shuffle, szip_on
  integer :: gzip_level
  logical :: extensible_on
  ! Maximum chunk size to use in each dimension
  integer, parameter                  :: nchunk_max = 1000
  integer(kind=hsize_t), dimension(7) :: nchunk
  ! HDF5 identifiers
  integer(kind=hid_t)                 :: dset_id, dspace_id, dtype_id, memtype_id
  integer(kind=hid_t)                 :: memspace_id, group_id
  ! Dataset opening / creation
  logical :: dset_opened
  integer(hid_t)                      :: prop_id
  ! Size of the input array
  integer                             :: array_rank, select_rank
  integer(kind=hsize_t)               :: array_dims(7), select_dims(7)
  integer(kind=hsize_t)               :: hs_start(7), hs_count(7)
  !
  integer(kind=hsize_t)               :: eff_select_rank, eff_array_rank, pad_dims(7)
  integer(kind=hsize_t)               :: eff_select_dims(7), eff_array_dims(7), eff_select_idx(7)

  ! Dimensions of the file dataspace
  integer                             :: file_rank
  integer(kind=hsize_t)               :: file_dims(7), maxdims(7)
  integer :: hdf_err
  ! Data type sizes and classes
  integer(kind=size_t)  :: dtype_size,  memtype_size
  integer               :: dtype_class, memtype_class

  ! Check input parameters
  ! ----------------------

  if(present(gzip).and.present(szip))then
     if(gzip.gt.0.and.szip)then
        call hdf5_abort( &
             'Can only specify one compression method in write_data!')
     endif
  endif

  ! Add 10 to the gzip level to apply shuffle filter
  gzip_shuffle = .false.
  gzip_level = 0
  szip_on = .false.
  if(present(gzip))then
     gzip_level = gzip
     if(gzip.gt.9)then
        gzip_level = gzip - 10
        gzip_shuffle = .TRUE.
     endif
  endif
  if(present(szip))szip_on = szip
  extensible_on = .false.
  if(present(extensible))extensible_on = .true.

  if(ifile.lt.1.or.ifile.gt.nfilemax)call hdf5_abort( &
       'Invalid file handle in write_data',name=name)

  if(file_id(ifile).lt.0) call hdf5_abort( &
       'File is not open in write_data()!', &
       name=name)

  if(read_only(ifile))call hdf5_abort( &
       "Attempt to write dataset to read only file!", &
       name=name,fname=fname_table(ifile))

  if(present(start).and.(.not.present(count)))call hdf5_abort( &
       'If parameter start is present count must be present in read_dataset')
  if(present(count).and.(.not.present(start)))call hdf5_abort( &
       'If parameter count is present start must be present in read_dataset')

  ! Get the data type for the data in memory
#ifndef STRING
  call h5tcopy_f( NATIVE_TYPE, memtype_id, hdf_err)
  if(hdf_err.lt.0)call hdf5_abort( &
       'Unable to open datatype in write_data()!', &
       name=name,fname=fname_table(ifile))
#else
  call h5tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, hdf_err)
  call h5tset_size_f(memtype_id, INT(LEN(arr),KIND=SIZE_T), hdf_err) 
#endif

  ! Open or create the dataset
  ! --------------------------
  
  call h5eset_auto_f(0, hdf_err)
  call h5dopen_f(file_id(ifile), name, dset_id, hdf_err)
  dset_opened = (hdf_err.eq.0)
  if(HDF_ERROR_SUPPRESS .eq. 0)call h5eset_auto_f(1, hdf_err)
  
  if(dset_opened.and.present(overwrite))then
     if(overwrite)then
        ! If the dataset exists but we have overwrite = .true.,
        ! delete the dataset
        call h5dclose_f(dset_id, hdf_err)
        call h5gunlink_f(file_id(ifile),name,hdf_err)
        dset_opened = .false.
     endif
  endif

  ! Store size of the input array
  array_rank          = NDIMS
#ifndef SCALAR
  do i = 1, NDIMS, 1
     array_dims(i) = size(arr,i)
  end do
#endif

  ! Why doesn't this work?!
  ! array_dims(1:NDIMS) = shape(arr(1:NDIMS))

  ! If we don't have a dataset open at this point, we'll have to
  ! create a new one. First need to determine dataset creation
  ! properties.
  if(.not.dset_opened)then

     ! Create the group that will contain the dataset if necessary
     i = index(name,"/",back=.true.)
     if(i.gt.1)then
        call h5eset_auto_f(0, hdf_err)
        call h5gopen_f(file_id(ifile), name(1:i-1), group_id, hdf_err)
        if(hdf_err.ne.0)then
           call hdf5_create_group(ifile,trim(name(1:i-1)))
        else
           call h5gclose_f(group_id, hdf_err)
        endif
        if(HDF_ERROR_SUPPRESS .eq. 0)call h5eset_auto_f(1, hdf_err)
     endif

     ! Create the data space for the new dataset using the dimensions of
     ! the supplied array, or start and count if present.
     file_rank = array_rank
     
     if(present(initial_size))then
        if(size(initial_size).lt.file_rank)call hdf5_abort( &
             'initial_size parameter does not have enough elements in write_data()',&
             name=name,fname=fname_table(ifile))
        ! Use initial size if specified
        file_dims(1:file_rank) = initial_size(1:file_rank)
     else
        ! Otherwise, just make the dataset big enough to contain the
        ! data
        if(present(start))then
           if(size(count).lt.file_rank.or.size(start).lt.file_rank)call hdf5_abort( &
                'start/count parameters do not have enough elements in write_data()',&
                name=name,fname=fname_table(ifile))
           file_dims(1:file_rank) = start(1:file_rank)+count(1:file_rank)-1
        else
           file_dims(1:file_rank) = array_dims(1:file_rank)
        endif
     endif
        
     if(extensible_on)then
        maxdims(1:file_rank) = H5S_UNLIMITED_F
     else
        maxdims(1:file_rank) = file_dims(1:file_rank)
     endif
     call h5screate_simple_f(file_rank,file_dims,dspace_id,hdf_err,maxdims) 
     if(hdf_err.ne.0)call hdf5_abort( &
          'Unable to create dataspace in write_data()!', &
          name=name,fname=fname_table(ifile))

     
     ! Dataset creation property list
     call h5pcreate_f(H5P_DATASET_CREATE_F, prop_id, hdf_err) 
     ! Decide on chunk size based on the initial size of the dataset
     do i = 1, file_rank, 1
        nchunk(i) = min(file_dims(i),nchunk_max)
     end do
     ! Only use chunking if we need it (saves space for small datasets)
     if(gzip_level.gt.0.or.szip_on.or.extensible_on) &
          call h5pset_chunk_f(prop_id, NDIMS, nchunk, hdf_err)
     ! Enable compression if required
     if(gzip_level.gt.0)then
        if(gzip_shuffle)then
           call h5pset_shuffle_f(prop_id, hdf_err)
           if(hdf_err.lt.0)write(*,*) &
                '[hdf5_write_dataset] Dataset ',trim(name),' will use the shuffle filter'
        endif
        call h5pset_deflate_f(prop_id,gzip_level,hdf_err)
        if(hdf_err.eq.0)then
           if(HDF_VERBOSITY .ge. 1)write(*,*) &
                '[hdf5_write_dataset] Dataset ',trim(name),' will be gzipped'
        else
           if(hdf_err.lt.0)write(*,*) &
                '[hdf5_write_dataset] Unable to gzip data.'//&
                ' Will write uncompressed data.'
        endif
     endif
     ! Szip compression isn't always enabled in the HDF5 library
#ifdef USE_SZIP
     if(szip_on)call H5Pset_szip_f(prop_id, H5_SZIP_NN_OM_F, 8, hdf_err)
#else
     if(szip_on.and.HDF_VERBOSITY.ge.1)write(*,*) &
          "[hdf5_write_dataset] This version of the wrapper was compiled"//&
          " without szip support. Dataset will not be compressed."
#endif
     
     ! Set the datatype for the new dataset
     call h5tcopy_f(memtype_id, dtype_id, hdf_err)

     ! Create the dataset
     call h5dcreate_f(file_id(ifile), name, dtype_id, &
          dspace_id, dset_id, hdf_err,prop_id)
     if(hdf_err.lt.0)call hdf5_abort( &
          'Unable to create dataset in write_data()!', &
          name=name,fname=fname_table(ifile))
     
     ! Finished with the creation property list
     call h5pclose_f(prop_id, hdf_err)

  else

     ! Get the rank and dimensions of the existing dataset
     call h5dget_type_f(dset_id, dtype_id, hdf_err)
     call h5dget_space_f(dset_id, dspace_id, hdf_err) 
     call h5sget_simple_extent_ndims_f(dspace_id, file_rank, hdf_err) 
     call h5sget_simple_extent_dims_f(dspace_id, file_dims, maxdims, hdf_err)

     ! The input array must be of the same class and precision as the
     ! dataset
     call h5tget_size_f(dtype_id,    dtype_size,    hdf_err)
     call h5tget_class_f(dtype_id,   dtype_class,   hdf_err) 
     call h5tget_size_f(memtype_id,  memtype_size,  hdf_err)
     call h5tget_class_f(memtype_id, memtype_class, hdf_err) 
     if(memtype_class.ne.dtype_class)call hdf5_abort( &
          'Unable to write this type of data to existing data set',&
          name=name,fname=fname_table(ifile))
#ifndef STRING
     if(memtype_size.ne.dtype_size)call hdf5_abort( &
          'Precision of supplied array doesn''t match precision of existing data set',&
          name=name,fname=fname_table(ifile))
#else
     if(memtype_size.gt.dtype_size)call hdf5_abort( &
          'Length of string(s) to write is greater than length of string(s) in data set',&
          name=name,fname=fname_table(ifile))    
#endif

     ! Check start/count have enough elements if present
     if(present(start))then
        if(size(start).lt.file_rank.or.size(count).lt.file_rank) call hdf5_abort( &
             'start/count parameters do not have enough elements in write_data', &
             name=name,fname=fname_table(ifile))
     else
        ! Don't allow writing to an old dataset without start/count
        call hdf5_abort( &
             'Must specify start and count parameters to write to an existing dataset',&
             name=name,fname=fname_table(ifile))
     endif

     ! Extend the dataset if necessary
     hs_count(1:file_rank) = &
          max(file_dims(1:file_rank), start(1:file_rank)+count(1:file_rank)-1)
     call h5dextend_f(dset_id, hs_count, hdf_err)
     if(hdf_err.lt.0)call hdf5_abort( &
          "Unable to extend dataset in write_data() - "//&
          "try creating it with extensible=.true.", &
          name=name,fname=fname_table(ifile))

     ! Update the dimensions of the dataset
     call h5sclose_f(dspace_id,hdf_err)
     call h5dget_space_f(dset_id, dspace_id, hdf_err)
     call h5sget_simple_extent_ndims_f(dspace_id, file_rank, hdf_err)
     call h5sget_simple_extent_dims_f(dspace_id, file_dims, maxdims, hdf_err)

  endif

!
! At this point the dataset in the file has been opened/created, and we have:
!
! dset_id    - dataset handle
! dtype_id   - data type handle for data type in file
! dspace_id  - data space handle for data in file
! memtype_id - data type handle for data type in memory
! 
! file_rank / file_dims - dimensions of dataset in the file
!

  ! Get the dimensions of the part of the dataset to write
  select_rank = file_rank
  if(present(start))then
     select_dims(1:select_rank) = count(1:select_rank)
  else
     select_dims(1:select_rank) = file_dims(1:select_rank)
  endif

  ! Get the dimensions of the array in memory. These may need some
  ! padding with 1 element dimensions so that the rank is the same
  ! as the dataset in the file.
  j = 0
  do i = 1, file_rank, 1
     if(select_dims(i).gt.1)then
        j = j + 1
        eff_select_dims(j) = select_dims(i)
        eff_select_idx(j) = i
     endif
  end do
  eff_select_rank = j

  j = 0
  do i = 1, array_rank, 1
     if(array_dims(i).gt.1)then
        j = j + 1
        eff_array_dims(j) = array_dims(i)
     endif
  end do
  eff_array_rank = j

  ! Check dataset and array have the same rank, ignoring dimensions of size 1
  if(eff_array_rank.ne.eff_select_rank)call hdf5_abort( &
       'The supplied data does not have the same number'//&
       ' of dimensions as the dataset to write',&
       name=name,fname=fname_table(ifile))    

  ! Check that the array contains enough data
  do i = 1, eff_array_rank, 1
     if(eff_array_dims(i).lt.eff_select_dims(i))call hdf5_abort( &
          'Size of the data array is inconsistent with '//&
          'the count parameter in write_data',&
          name=name,fname=fname_table(ifile))    
  end do

  ! Generate the padded dimensions array.
  pad_dims(1:7) = 1
  do i = 1, eff_array_rank, 1
     pad_dims(eff_select_idx(i)) = eff_array_dims(i)
  end do

  ! Use this to create the dataspace to describe the array in memory
  maxdims(1:file_rank) = pad_dims(1:file_rank)
  call h5screate_simple_f(file_rank,pad_dims,memspace_id,hdf_err,maxdims) 
  
  ! Select the part of the array in memory which will be written out
  if(present(start))then
     hs_start(1:file_rank) = 0
     hs_count(1:file_rank) = count(1:file_rank)
  else
     hs_start(1:file_rank) = 0
     hs_count(1:file_rank) = file_dims(1:file_rank)
  endif
  if(file_rank.gt.0)then ! Can't use select hyperslab on scalar dataspace
     call h5sselect_hyperslab_f(memspace_id, H5S_SELECT_SET_F, &
          hs_start, hs_count, hdf_err)
     if(hdf_err.ne.0)call hdf5_abort( &
          'Unable to select array elements to write in write_data',&
          name=name,fname=fname_table(ifile))    
  endif

  ! Select the part of the dataset in the file that will be written to
  if(present(start))then
     ! Select region specified by start/count
     hs_start(1:file_rank) = start(1:file_rank) - 1
     hs_count(1:file_rank) = count(1:file_rank)
  else
     ! Select whole dataset
     hs_start(1:file_rank) = 0
     hs_count(1:file_rank) = file_dims(1:file_rank)
  endif
  if(file_rank.gt.0)then
     call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, &
          hs_start, hs_count, hdf_err)
     if(hdf_err.ne.0)call hdf5_abort( &
          'Unable to select part of dataset to write in write_data',&
          name=name,fname=fname_table(ifile)) 
  endif

  ! Write out the data
  ! ------------------

  if (HDF_VERBOSITY .ge. 1) then
     write(*,*)'[hdf5_write_dataset] Writing dataset ',trim(name)
  endif
 
  call write_hdf5_dataset(dset_id, memtype_id, dspace_id, memspace_id, &
       arr, hdf_err)
  if(hdf_err.ne.0)call hdf5_abort( &
       'Unable to write dataset in write_data()',&
       name=name,fname=fname_table(ifile))

  call h5tclose_f(memtype_id,hdf_err)
  call h5tclose_f(dtype_id,hdf_err)
  call h5sclose_f(dspace_id,hdf_err)
  call h5sclose_f(memspace_id,hdf_err)
  call h5dclose_f(dset_id, hdf_err)

  return
end subroutine SUB_NAME
