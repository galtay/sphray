subroutine SUB_NAME (ifile,name,arr)
  !
  ! Read an n-D array dataset from an open file. This assumes some
  ! preprocessor variables have been set - see write_attribute_preprocessor.F90.
  !
      implicit none
      integer, parameter              :: LEN_STR = 256
      integer, parameter              :: MAX_NEST = 10
      character(len=*), intent(in)    :: name 
      integer,intent(in)              :: ifile
      integer                         :: hdf_err
      integer(hid_t)                  :: loc_id, attr_id, dspace_id
      integer(hid_t)                  :: dtype_id, memtype_id
#ifndef SCALAR
      integer(hsize_t)                :: dimensions(NDIMS)
      ARR_TYPE , dimension ARRAY_DIM  :: arr
#else
      ARR_TYPE                        :: arr
#endif
      character(len=LEN_STR)          :: loc_name
      character(len=LEN_STR)          :: attr_name
      integer                         :: nslash
      integer                         :: itype,i
!     
      if(ifile.lt.1.or.ifile.gt.nfilemax)call hdf5_abort( &
           'Invalid file handle in read_attribute',name=name)

      if(file_id(ifile).lt.0) then 
         call hdf5_abort('File is not open in hdf5_read_attribute()!',&
              name=name)
      endif

#ifndef STRING
      call h5tcopy_f( NATIVE_TYPE, memtype_id, hdf_err)
      if(hdf_err.lt.0)call hdf5_abort( &
           'Unable to open datatype in read_data()!', &
           name=name,fname=fname_table(ifile))
#else
      call h5tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, hdf_err)
      call h5tset_size_f(memtype_id, INT(LEN(arr),KIND=SIZE_T), hdf_err) 
#endif

!     
      ! Generate array of dimensions:
#ifndef SCALAR
      do i=1,NDIMS
         dimensions(i) = ubound(arr,i) - lbound(arr,i) + 1
      enddo
#endif
      
      ! Split name into path and attribute name
      nslash=index(name,"/",.true.)
      if(nslash.eq.0) then
         call hdf5_abort('Invalid attribute name in read_attribute!', &
              name=name,fname=fname_table(ifile))
      endif
      
      loc_name=name(1:nslash-1)
      attr_name=name(nslash+1:len_trim(name))
      ! Try to open loc_name as a group and as a dataset
      itype=1
      call h5eset_auto_f(0, hdf_err)      
      call h5dopen_f(file_id(ifile),loc_name,loc_id,hdf_err)
      if(hdf_err.lt.0)then
         call h5gopen_f(file_id(ifile),loc_name, loc_id, hdf_err)
         if(hdf_err.lt.0) then
            call hdf5_abort( &
                 'Unable to open attribute parent object in read_attribute!', &
                 name=name,fname=fname_table(ifile))
         endif 
         itype=2 
      end if 
      if (HDF_ERROR_SUPPRESS .eq. 0) call h5eset_auto_f(1, hdf_err) 
      
      call h5aopen_name_f(loc_id, attr_name, attr_id, hdf_err) 
	   if(hdf_err.lt.0) then  
              call hdf5_abort('Unable to open attribute in read_attribute()!',&
                   name=name,fname=fname_table(ifile))
      endif

      !Call the c version of this routine
      CALL read_hdf5_attribute(attr_id, memtype_id , arr, hdf_err)
      if(hdf_err.lt.0) then 
         call hdf5_abort &
              ('Unable to read attribute in read_attribute()!', &
              name=name,fname=fname_table(ifile))
      endif

      call h5aclose_f(attr_id,hdf_err)
      if(hdf_err.lt.0) then
         call hdf5_abort &
              ('Unable to close attribute in read_attribute()!', &
              name=name,fname=fname_table(ifile))
      endif

      if(itype.eq.1)then
         call h5dclose_f(loc_id,hdf_err)
         if(hdf_err.lt.0) then
            call hdf5_abort( &
                 'Unable to close dataset in read_attribute()!', &
                 name=name,fname=fname_table(ifile))
         endif
      else
         call h5gclose_f(loc_id,hdf_err)
         if(hdf_err.lt.0) then
            call hdf5_abort( &
                 'Unable to close group in read_attribute()!', &
                 name=name,fname=fname_table(ifile))
         endif
      end if
      
end subroutine SUB_NAME 


