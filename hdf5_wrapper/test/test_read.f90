PROGRAM test

  USE hdf5_wrapper
  IMPLICIT NONE
  CHARACTER(LEN=256)          :: fname
  INTEGER                     :: ifile
  INTEGER*4                   :: dozen,intattr
  integer*4, dimension(5,3)   :: i4attr
  real, dimension(5,3)        :: r4attr
  real*8, dimension(3,100000) :: r8arr
  integer*4, dimension(10,20) :: i4arr
  INTEGER :: i,j,k
  CHARACTER(LEN=20) :: str(5,5)
  integer, dimension(7) :: start,count
  integer, dimension(10) :: i1d
  character(len=50) :: dtype
  integer :: size
  integer :: nsets
  character(len=200), dimension(100) :: set_name

  fname = "./test.hdf5"

  HDF_VERBOSITY=2

  CALL hdf5_open_file(ifile, fname, readonly=.true.)

  call hdf5_read_data(ifile,"Dozen",Dozen)
  write(*,*)'One Dozen is : ',dozen

  call hdf5_read_data(ifile,"/Arrays/Real/R8_Array",r8arr(1:1,:))
  write(*,*)'R8_Array(1:10) = ',r8arr(1:1,1:10)

  call hdf5_read_attribute(ifile,"Arrays/I4_Attribute",i4attr)
  write(*,*)'The 2d integer*4 attribute is : ',i4attr

  call hdf5_read_attribute(ifile,"Arrays/I4_Scalar",intattr)
  write(*,*)'The scalar integer*4 attribute is : ',intattr

  call hdf5_read_attribute(ifile,"Arrays/R4_Attribute",r4attr)
  write(*,*)'The 2d real*4 attribute is : ',r4attr

  call hdf5_read_data(ifile,"Arrays/Test-String",str)
  DO i=1,2
     DO j=1,2
        WRITE(*,*)'String: ',i,j,TRIM(str(i,j))
     END DO
  END DO

  call hdf5_read_data(ifile,"/Arrays/I4_Array",i4arr)
  write(*,*)'I4_Array:'
  write(*,*)i4arr(1,1),i4arr(2,1),i4arr(3,1),i4arr(4,1),i4arr(5,1)
  write(*,*)i4arr(1,2),i4arr(2,2),i4arr(3,2),i4arr(4,2),i4arr(5,2)
  write(*,*)i4arr(1,3),i4arr(2,3),i4arr(3,3),i4arr(4,3),i4arr(5,3)

  i4arr=i4arr*0
  start(1)=3
  start(2)=2
  count(1)=3
  count(2)=2
  call hdf5_read_data(ifile,"/Arrays/I4_Array",i4arr,start=start,count=count)
  write(*,*)'Subsection from I4_Array:'
  write(*,*)i4arr(1,1),i4arr(2,1),i4arr(3,1),i4arr(4,1),i4arr(5,1)
  write(*,*)i4arr(1,2),i4arr(2,2),i4arr(3,2),i4arr(4,2),i4arr(5,2)
  write(*,*)i4arr(1,3),i4arr(2,3),i4arr(3,3),i4arr(4,3),i4arr(5,3)

  start(1)=4
  start(2)=2
  count(1)=1
  count(2)=1
  call hdf5_read_data(ifile,"/Arrays/I4_Array",i,start=start,count=count)
  write(*,*)'One element from I4_Array:'
  write(*,*)i

  start(1)=2
  start(2)=2
  count(1)=3
  count(2)=1
  i4arr=i4arr*0
  call hdf5_read_data(ifile,"/Arrays/I4_Array",i4arr(3,1:3),start=start,count=count)
  write(*,*)'Part of a row from I4_Array read into the third column of a 2D array:'
  write(*,*)i4arr(1,1),i4arr(2,1),i4arr(3,1),i4arr(4,1),i4arr(5,1)
  write(*,*)i4arr(1,2),i4arr(2,2),i4arr(3,2),i4arr(4,2),i4arr(5,2)
  write(*,*)i4arr(1,3),i4arr(2,3),i4arr(3,3),i4arr(4,3),i4arr(5,3)

  start(1)=2
  start(2)=1
  count(1)=1
  count(2)=3
  i4arr=i4arr*0
  call hdf5_read_data(ifile,"/Arrays/I4_Array",i4arr(1:3,3),start=start,count=count)
  write(*,*)'A column from I4_Array read into the third row of a 2D array:'
  write(*,*)i4arr(1,1),i4arr(2,1),i4arr(3,1),i4arr(4,1),i4arr(5,1)
  write(*,*)i4arr(1,2),i4arr(2,2),i4arr(3,2),i4arr(4,2),i4arr(5,2)
  write(*,*)i4arr(1,3),i4arr(2,3),i4arr(3,3),i4arr(4,3),i4arr(5,3)

  CALL hdf5_get_type(ifile,"Dozen",dtype,size)
  write(*,*)'Type and size of "Dozen" are ',trim(dtype),size

  CALL hdf5_get_type(ifile,"/Arrays/Real/R8_Array",dtype,size)
  write(*,*)'Type and size of "R8_Array" are ',trim(dtype),size

  write(*,*)'Datasets in root group:'
  call hdf5_list_datasets(ifile,"/",nsets,set_name)
  do i = 1, nsets, 1
     write(*,*)"  ",trim(set_name(i))
  end do

  write(*,*)'Datasets in group "Arrays":'
  call hdf5_list_datasets(ifile,"/Arrays",nsets,set_name)
  do i = 1, nsets, 1
     write(*,*)"  ",trim(set_name(i))
  end do

  write(*,*)'Groups in root group:'
  call hdf5_list_groups(ifile,"/",nsets,set_name)
  do i = 1, nsets, 1
     write(*,*)"  ",trim(set_name(i))
  end do

  write(*,*)'Attributes in group "Arrays":'
  call hdf5_list_attributes(ifile,"/Arrays",nsets,set_name)
  do i = 1, nsets, 1
     write(*,*)"  ",trim(set_name(i))
  end do

  CALL hdf5_close_file(ifile)



END PROGRAM test
