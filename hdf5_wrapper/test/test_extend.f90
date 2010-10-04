program test_extend
!
! Test random access and extensible dataset functions
! in the hdf5 wrapper
!
  use hdf5_wrapper
  implicit none
  ! File parameters
  character(len=500) :: fname
  integer            :: ifile
  ! Data to write
  integer, parameter       :: nmax = 100
  integer, dimension(nmax) :: idata
  integer, parameter       :: nx = 10, ny = 5
  real,    dimension(nx,ny) :: rdata
  real, dimension(nmax) :: r1d
  integer :: i, j
  character(len=3), dimension(3,5) :: str

  fname = "./test2.hdf5"
  HDF_VERBOSITY=2
  HDF_ERROR_SUPPRESS = 0
  
  str = "abc"

  do i = 1, nmax, 1
     idata(i) = i+(10000*i)
  end do

  do i = 1, nx, 1
     do j = 1, ny, 1
        rdata(i,j) = i+1000*j
     end do
  end do

  call hdf5_create_file(ifile, fname)
  
  ! Write out the data sequentially
  do i = 1, nmax, 10
     call hdf5_write_data(ifile,"idata",idata(i:i+9),extensible=.true., &
          start=(/i/), count = (/10/))
  end do

  ! Try overwriting elements 50-59 with -1
  idata(1:10) = -1
  call hdf5_write_data(ifile,"idata",idata,extensible=.true., &
       start=(/50/), count = (/10/))
  
  ! Write out the complete 2D array
  call hdf5_write_data(ifile,"rdata1",rdata)

  ! Write out the 2D array one row at a time
  do i = 1, nx, 1
     call hdf5_write_data(ifile,"rdata2",rdata(i:i,1:ny),extensible=.true., &
          start=(/i,1/), count = (/1,ny/))
  end do

  ! Overwrite some elements at the end of rdata1
  rdata = -1000
  call hdf5_write_data(ifile,"rdata1",rdata(1:3,1:2),extensible=.true., &
       start=(/8,4/), count = (/3,2/))

  ! Create a dataset with some unused elements
  rdata = 10
  call hdf5_write_data(ifile,"rdata3",rdata,start=(/1,ny+1/),count=(/nx,ny/))

  ! Create a dataset thats larger than the data initially written to it
  call hdf5_write_data(ifile,"strdata", str, initial_size=(/4, 8/), &
       start=(/1,1/), count=(/3,5/))
  str="def"
  ! Write some more data to another part of the dataset
  call hdf5_write_data(ifile,"strdata", str, start=(/2,4/), count=(/3,5/))

  ! Write out a 2D array of 0's
  rdata = 0.0
  call hdf5_write_data(ifile,"Array2D",rdata)
  
  ! Overwrite some 1D subsections of the 2D dataset
  r1d = 1.0
  call hdf5_write_data(ifile,"Array2D",r1d(1:nx),start=(/1,1/),count=(/nx-2,1/))
  r1d = 2.0
  call hdf5_write_data(ifile,"Array2D",r1d(1:ny),start=(/3,1/),count=(/1,ny/))

  ! Overwrite a single element with a scalar
  call hdf5_write_data(ifile,"Array2D",3.0,start=(/4,4/),count=(/1,1/))


  ! Write one element of a 1D integer dataset
  idata(1:5) = (/ 1,2,3,4,5 /)
  call hdf5_write_data(ifile,"itest",idata(1:1),start=(/1/),count=(/1/),&
       initial_size=(/5/))
  call hdf5_close_file(ifile)

end program test_extend
