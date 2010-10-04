PROGRAM test

  USE hdf5_wrapper
  IMPLICIT NONE
  CHARACTER(LEN=256) :: fname
  INTEGER :: ifile
  INTEGER*4, DIMENSION(5,3)   :: iarr
  integer*4                   :: i4scalarattr
  INTEGER*4, DIMENSION(5,3)   :: i4attr
  INTEGER*8, DIMENSION(3,3,3) :: i8attr
  REAL*4, DIMENSION(5,3)      :: r4attr
  REAL*8, DIMENSION(5,3)      :: r8attr
  INTEGER*8, DIMENSION(2,6)   :: i8arr
  REAL*8, DIMENSION(100000)       :: r8arr
  REAL*4, DIMENSION(2,2,2)    :: r4arr
  INTEGER :: i,j,k
  CHARACTER(LEN=20) :: str(10,10)

  REAL :: rscalar = 3.14159

  i4scalarattr = 10

  DO i=1,5
     DO j=1,3
        iarr(i,j)=1000*i+j
        i4attr(i,j) = i+j
        r4attr(i,j) = 2*i+j
        r8attr(i,j) = 4*i+j
     END DO
  END DO

  DO i=1,2
     DO j=2,6
        i8arr(i,j)=(i+j)*2
     END DO
  END DO

  DO i=1,10,1
     r8arr(i)=1.0d0*i
  END DO

  DO i=1,2
     DO j=1,2
        DO k=1,2
           r4arr(i,j,k)=REAL(i+j+k)
        END DO
     END DO
  END DO

  DO i=1,3
     DO j=1,3
        DO k=1,3
           i8attr(i,j,k)=(i+j+k)
        END DO
     END DO
  END DO

  i=12

  r8arr=12

  DO i=1,2,1
     DO j=1,2,1
        WRITE(str(i,j),'("String element ",i1,",",i1)')i,j
     END DO
  END DO

  !call hdf5_auto_open(.false.)

  fname = "./test.hdf5"
  HDF_VERBOSITY=2
  HDF_ERROR_SUPPRESS = 0
  CALL hdf5_create_file(ifile, fname)
  CALL hdf5_write_data(ifile,"Pi",rscalar)
  CALL hdf5_write_data(ifile,"Dozen",12)
  CALL hdf5_write_data(ifile,"/Arrays/I4_Array",iarr)
  CALL hdf5_write_data(ifile,"/Arrays/I8_Array",i8arr,overwrite=.true.)
  CALL hdf5_write_data(ifile,"/Arrays/Real/R8_Array",r8arr,szip=.false.)
  CALL hdf5_write_data(ifile,"/Arrays/Real/R4_Array",r4arr)
  CALL hdf5_write_data(ifile,"/Arrays/Test-String",str(1:2,1:2))
  CALL hdf5_write_attribute(ifile,"/Arrays/I4_Attribute",i4attr, &
       overwrite=.true.)
  CALL hdf5_write_attribute(ifile,"/Arrays/I8_Attribute",i8attr)
  CALL hdf5_write_attribute(ifile,"/Arrays/R4_Attribute",r4attr)
  CALL hdf5_write_attribute(ifile,"/Arrays/R8_Attribute",r8attr)
  CALL hdf5_write_attribute(ifile,"/Arrays/I4_Scalar",i4scalarattr)
  CALL hdf5_create_group(ifile,"Group1/Group2")
  CALL hdf5_create_group(ifile,"Group1/Group2/GroupA")
  CALL hdf5_create_group(ifile,"Group1/Group2/GroupB")
  CALL hdf5_close_file(ifile)


  ! Re-open and try to overwrite an attribute and a dataset
  CALL hdf5_open_file(ifile, "test.hdf5", readonly=.false.)
  CALL hdf5_write_attribute(ifile,"/Arrays/I4_Scalar",20,overwrite=.true.)
  CALL hdf5_write_data(ifile,"/Arrays/Real/R4_Array",r4arr*0.0+16.0, &
       overwrite=.true.)
  CALL hdf5_close_file(ifile)


END PROGRAM test


