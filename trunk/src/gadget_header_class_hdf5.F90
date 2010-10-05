!> \file gadget_header_class_hdf5.f90

!> \brief Handles HDF5 GADGET style headers.  
!!
!! Contains the means to read/write headers from a file.  
!< 

module gadget_header_class_hdf5
use myf03_mod
use gadget_public_header_class
use hdf5_io_mod
use hdf5
implicit none
private

public :: gadget_public_header_hdf5_type



type, extends(gadget_public_header_type) :: gadget_public_header_hdf5_type
 contains
   procedure :: read_file_hdf5 => read_header_from_file_hdf5
   procedure :: write_file_hdf5 => write_header_to_file_hdf5
end type gadget_public_header_hdf5_type


interface rw_attr
   module procedure rw_attr_int,  rw_attr_int_arr, &
                    rw_attr_real, rw_attr_real_arr
end interface



contains



!> reads a gadget header from an hdf5 file
!--------------------------------------------------------------
subroutine read_header_from_file_hdf5(this, snapfile)
  class(gadget_public_header_hdf5_type) :: this
  character(*), intent(in) :: snapfile

  integer :: ierr
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: attr_id

  integer(HSIZE_T), dimension(1) :: dims

  dims(1) = 6

  ! initialize fortran interface
  call h5open_f(ierr)

  ! open file
  call h5fopen_f(snapfile, H5F_ACC_RDONLY_F, file_id, ierr)

  ! open group
  call h5gopen_f(file_id, '/Header', grp_id, ierr) 

  ! read header attributes
  call rw_attr(grp_id, 'NumPart_ThisFile', H5T_NATIVE_INTEGER, this%npar_file, 'r')
  call rw_attr(grp_id, 'NumPart_Total', H5T_NATIVE_INTEGER, this%npar_all, 'r')
  call rw_attr(grp_id, 'NumPart_Total_HW', H5T_NATIVE_INTEGER, this%npar_hw, 'r')
  call rw_attr(grp_id, 'MassTable', H5T_NATIVE_DOUBLE, this%mass, 'r')
  call rw_attr(grp_id, 'Time', H5T_NATIVE_DOUBLE, this%a, 'r')
  call rw_attr(grp_id, 'Redshift', H5T_NATIVE_DOUBLE, this%z, 'r')
  call rw_attr(grp_id, 'BoxSize', H5T_NATIVE_DOUBLE, this%boxlen, 'r')
  call rw_attr(grp_id, 'NumFilesPerSnapshot', H5T_NATIVE_INTEGER, this%nfiles, 'r')
  call rw_attr(grp_id, 'Omega0', H5T_NATIVE_DOUBLE, this%OmegaM, 'r')
  call rw_attr(grp_id, 'OmegaLambda', H5T_NATIVE_DOUBLE, this%OmegaL, 'r')
  call rw_attr(grp_id, 'HubbleParam', H5T_NATIVE_DOUBLE, this%h, 'r')
  call rw_attr(grp_id, 'Flag_Sfr', H5T_NATIVE_INTEGER, this%flag_sfr, 'r')
  call rw_attr(grp_id, 'Flag_Cooling', H5T_NATIVE_INTEGER, this%flag_cooling, 'r')
  call rw_attr(grp_id, 'Flag_StellarAge', H5T_NATIVE_INTEGER, this%flag_age, 'r')
  call rw_attr(grp_id, 'Flag_Metals', H5T_NATIVE_INTEGER, this%flag_metals, 'r')
  call rw_attr(grp_id, 'Flag_Feedback', H5T_NATIVE_INTEGER, this%flag_feedback, 'r')
  call rw_attr(grp_id, 'Flag_Entropy_ICs', H5T_NATIVE_INTEGER, this%flag_entr_ics, 'r')


  ! close everything else up
  call h5gclose_f(grp_id, ierr)
  call h5fclose_f(file_id, ierr)
  call h5close_f(ierr)



end subroutine read_header_from_file_hdf5



!> writes a gadget header to an hdf5 file
!--------------------------------------------------------------
subroutine write_header_to_file_hdf5(this, snapfile)
  class(gadget_public_header_hdf5_type) :: this
  character(*), intent(in) :: snapfile

  integer :: ierr
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: attr_id
  integer(HID_T) :: space_id

  integer :: rank
  integer(HSIZE_T), dimension(1) :: dims
  
  rank = 1
  dims(1) = 6

  ! initialize fortran interface
  call h5open_f(ierr)

  ! create a new file (EXCL fails if file exist, TRUNC overwrites)
  call h5fcreate_f(snapfile, H5F_ACC_EXCL_F, file_id, ierr)
  if (ierr < 0) then
     write(*,*) 
     write(*,*) ' ** Error ** '
     write(*,*) 'Routine: write_header_to_file_hdf5'
     write(*,*) 'File already exists: '
     write(*,*) trim(snapfile)
     write(*,*) 
     stop
  endif

  ! create the group
  call h5gcreate_f(file_id, '/Header', grp_id, ierr) 

  ! write header attributes
  call rw_attr(grp_id, 'NumPart_ThisFile', H5T_NATIVE_INTEGER, this%npar_file, 'w')
  call rw_attr(grp_id, 'NumPart_Total', H5T_NATIVE_INTEGER, this%npar_all, 'w')
  call rw_attr(grp_id, 'NumPart_Total_HW', H5T_NATIVE_INTEGER, this%npar_hw, 'w')
  call rw_attr(grp_id, 'MassTable', H5T_NATIVE_DOUBLE, this%mass, 'w')
  call rw_attr(grp_id, 'Time', H5T_NATIVE_DOUBLE, this%a, 'w')
  call rw_attr(grp_id, 'Redshift', H5T_NATIVE_DOUBLE, this%z, 'w')
  call rw_attr(grp_id, 'BoxSize', H5T_NATIVE_DOUBLE, this%boxlen, 'w')
  call rw_attr(grp_id, 'NumFilesPerSnapshot', H5T_NATIVE_INTEGER, this%nfiles, 'w')
  call rw_attr(grp_id, 'Omega0', H5T_NATIVE_DOUBLE, this%OmegaM, 'w')
  call rw_attr(grp_id, 'OmegaLambda', H5T_NATIVE_DOUBLE, this%OmegaL, 'w')
  call rw_attr(grp_id, 'HubbleParam', H5T_NATIVE_DOUBLE, this%h, 'w')
  call rw_attr(grp_id, 'Flag_Sfr', H5T_NATIVE_INTEGER, this%flag_sfr, 'w')
  call rw_attr(grp_id, 'Flag_Cooling', H5T_NATIVE_INTEGER, this%flag_cooling, 'w')
  call rw_attr(grp_id, 'Flag_StellarAge', H5T_NATIVE_INTEGER, this%flag_age, 'w')
  call rw_attr(grp_id, 'Flag_Metals', H5T_NATIVE_INTEGER, this%flag_metals, 'w')
  call rw_attr(grp_id, 'Flag_Feedback', H5T_NATIVE_INTEGER, this%flag_feedback, 'w')
  call rw_attr(grp_id, 'Flag_Entropy_ICs', H5T_NATIVE_INTEGER, this%flag_entr_ics, 'w')

  ! close everything else up
  call h5gclose_f(grp_id, ierr)
  call h5fclose_f(file_id, ierr)
  call h5close_f(ierr)

end subroutine write_header_to_file_hdf5




! read a scalar, integer, group attribute
!------------------------------------------------------------------
subroutine rw_attr_int(grp_id, attr_name, attr_mem_type, buf, rorw)

  integer(HID_T) :: grp_id
  character(*) :: attr_name
  integer(HID_T) :: attr_mem_type
  integer(i4b) :: buf
  character(1) :: rorw
 
  integer(HSIZE_T), dimension(1), parameter :: dims = [1]
  integer(i4b) :: ierr
  integer(HID_T) :: attr_id
  integer(HID_T) :: space_id
  
  if (rorw == 'r') then 
     call h5aopen_f(grp_id, attr_name, attr_id, ierr)
     call h5aread_f(attr_id, attr_mem_type, buf, dims, ierr)
     call h5aclose_f(attr_id, ierr)
  else if (rorw == 'w') then
     call h5screate_f(H5S_SCALAR_F, space_id, ierr) 
     call h5acreate_f(grp_id, attr_name, attr_mem_type, space_id, attr_id, ierr)
     call h5awrite_f(attr_id, attr_mem_type, buf, dims, ierr)
     call h5sclose_f(space_id, ierr)
     call h5aclose_f(attr_id, ierr)
  else
     write(*,*) 
     write(*,*) ' ** Error ** '
     write(*,*) 'Routine: rw_attr_int'
     write(*,*) 'rorw variable not "r" or "w" '
     write(*,*) 
     stop     
  endif

end subroutine rw_attr_int


! read a vector, integer, group attribute
!------------------------------------------------------------------
subroutine rw_attr_int_arr(grp_id, attr_name, attr_mem_type, buf, rorw)

  integer(HID_T) :: grp_id
  character(*) :: attr_name
  integer(HID_T) :: attr_mem_type
  integer(i4b) :: buf(0:5)
  character(1) :: rorw

  integer(HSIZE_T), dimension(1), parameter :: dims = [6]
  integer(i4b), parameter :: rank = 1
  integer(i4b) :: ierr
  integer(HID_T) :: attr_id
  integer(HID_T) :: space_id
  
  if (rorw == 'r') then 
     call h5aopen_f(grp_id, attr_name, attr_id, ierr)
     call h5aread_f(attr_id, attr_mem_type, buf, dims, ierr)
     call h5aclose_f(attr_id, ierr)
  else if (rorw == 'w') then
     call h5screate_simple_f(rank, dims, space_id, ierr)
     call h5acreate_f(grp_id, attr_name, attr_mem_type, space_id, attr_id, ierr)
     call h5awrite_f(attr_id, attr_mem_type, buf, dims, ierr)
     call h5sclose_f(space_id, ierr)
     call h5aclose_f(attr_id, ierr)
  else
     write(*,*) 
     write(*,*) ' ** Error ** '
     write(*,*) 'Routine: rw_attr_int_arr'
     write(*,*) 'rorw variable not "r" or "w" '
     write(*,*) 
     stop     
  endif

end subroutine rw_attr_int_arr



! read a scalar, real, group attribute
!------------------------------------------------------------------
subroutine rw_attr_real(grp_id, attr_name, attr_mem_type, buf, rorw)

  integer(HID_T) :: grp_id
  character(*) :: attr_name
  integer(HID_T) :: attr_mem_type
  real(r8b) :: buf
  character(1) :: rorw
 
  integer(HSIZE_T), dimension(1), parameter :: dims = [1]
  integer(i4b) :: ierr
  integer(HID_T) :: attr_id
  integer(HID_T) :: space_id
  
  if (rorw == 'r') then 
     call h5aopen_f(grp_id, attr_name, attr_id, ierr)
     call h5aread_f(attr_id, attr_mem_type, buf, dims, ierr)
     call h5aclose_f(attr_id, ierr)
  else if (rorw == 'w') then
     call h5screate_f(H5S_SCALAR_F, space_id, ierr) 
     call h5acreate_f(grp_id, attr_name, attr_mem_type, space_id, attr_id, ierr)
     call h5awrite_f(attr_id, attr_mem_type, buf, dims, ierr)
     call h5sclose_f(space_id, ierr)
     call h5aclose_f(attr_id, ierr)
  else
     write(*,*) 
     write(*,*) ' ** Error ** '
     write(*,*) 'Routine: rw_attr_real'
     write(*,*) 'rorw variable not "r" or "w" '
     write(*,*) 
     stop     
  endif

end subroutine rw_attr_real


! read a vector, real, group attribute
!------------------------------------------------------------------
subroutine rw_attr_real_arr(grp_id, attr_name, attr_mem_type, buf, rorw)

  integer(HID_T) :: grp_id
  character(*) :: attr_name
  integer(HID_T) :: attr_mem_type
  real(r8b) :: buf(0:5)
  character(1) :: rorw

  integer(HSIZE_T), dimension(1), parameter :: dims = [6]
  integer(i4b), parameter :: rank = 1
  integer(i4b) :: ierr
  integer(HID_T) :: attr_id
  integer(HID_T) :: space_id
  
  if (rorw == 'r') then 
     call h5aopen_f(grp_id, attr_name, attr_id, ierr)
     call h5aread_f(attr_id, attr_mem_type, buf, dims, ierr)
     call h5aclose_f(attr_id, ierr)
  else if (rorw == 'w') then
     call h5screate_simple_f(rank, dims, space_id, ierr)
     call h5acreate_f(grp_id, attr_name, attr_mem_type, space_id, attr_id, ierr)
     call h5awrite_f(attr_id, attr_mem_type, buf, dims, ierr)
     call h5sclose_f(space_id, ierr)
     call h5aclose_f(attr_id, ierr)
  else
     write(*,*) 
     write(*,*) ' ** Error ** '
     write(*,*) 'Routine: rw_attr_real_arr'
     write(*,*) 'rorw variable not "r" or "w" '
     write(*,*) 
     stop     
  endif

end subroutine rw_attr_real_arr

end module gadget_public_header_class_hdf5
