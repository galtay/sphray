!> \file gadget_public_header_class_hdf5.f90

!> \brief Handles HDF5 GADGET 2.0 Public version style headers.  
!!
!! Contains the means to read/write headers from a file.  
!< 

module gadget_public_header_class_hdf5
use myf03_mod
use gadget_public_header_class
#ifdef useHDF5
use hdf5_wrapper
#endif
implicit none
private
 
public :: gadget_public_header_hdf5_type


type, extends(gadget_public_header_type) :: gadget_public_header_hdf5_type
 contains
   procedure :: read_lun_hdf5 => read_header_from_lun_hdf5
   procedure :: read_file_hdf5 => read_header_from_file_hdf5
   procedure :: read_file1_hdf5 => read_header_from_file1_hdf5
   procedure :: write_lun_hdf5 => write_header_to_lun_hdf5
end type gadget_public_header_hdf5_type




contains


!> reads a gadget header from an hdf5 file
!--------------------------------------------------------------
subroutine read_header_from_lun_hdf5(this, fh)
  class(gadget_public_header_hdf5_type), intent(out) :: this
  integer, intent(in) :: fh

#ifdef useHDF5
  call hdf5_read_attribute(fh,'Header/NumPart_ThisFile',this%npar_file)
  call hdf5_read_attribute(fh,'Header/NumPart_Total',this%npar_all)
  call hdf5_read_attribute(fh,'Header/NumPart_Total_HW',this%npar_hw)
  call hdf5_read_attribute(fh,'Header/MassTable',this%mass)
  call hdf5_read_attribute(fh,'Header/Time',this%a)

  call hdf5_read_attribute(fh,'Header/Redshift',this%z)
  call hdf5_read_attribute(fh,'Header/BoxSize',this%boxlen)
  call hdf5_read_attribute(fh,'Header/NumFilesPerSnapshot',this%nfiles)
  call hdf5_read_attribute(fh,'Header/Omega0',this%OmegaM)
  call hdf5_read_attribute(fh,'Header/OmegaLambda',this%OmegaL)
  call hdf5_read_attribute(fh,'Header/HubbleParam',this%h)

  call hdf5_read_attribute(fh,'Header/Flag_Sfr',this%flag_sfr)
  call hdf5_read_attribute(fh,'Header/Flag_Cooling',this%flag_cooling)
  call hdf5_read_attribute(fh,'Header/Flag_StellarAge',this%flag_age)
  call hdf5_read_attribute(fh,'Header/Flag_Metals',this%flag_metals)
  call hdf5_read_attribute(fh,'Header/Flag_Feedback',this%flag_feedback)
  call hdf5_read_attribute(fh,'Header/Flag_Entropy_ICs',this%flag_entr_ics)
#endif

end subroutine read_header_from_lun_hdf5


!> reads a gadget header from an hdf5 file
!--------------------------------------------------------------
subroutine read_header_from_file_hdf5(this, snapfile)
  class(gadget_public_header_hdf5_type) :: this
  character(*), intent(in) :: snapfile
  integer :: fh

#ifdef useHDF5
  call hdf5_open_file( fh, snapfile, readonly=.true. )
  call read_header_from_lun_hdf5(this, fh)
  call hdf5_close_file( fh )
#endif

end subroutine read_header_from_file_hdf5



!> writes a gadget header to an hdf5 file
!--------------------------------------------------------------
subroutine write_header_to_lun_hdf5(this, fh)
  class(gadget_public_header_hdf5_type), intent(in) :: this
  integer, intent(in) :: fh

#ifdef useHDF5
  call hdf5_write_attribute(fh,'Header/NumPart_ThisFile',this%npar_file)
  call hdf5_write_attribute(fh,'Header/NumPart_Total',this%npar_all)
  call hdf5_write_attribute(fh,'Header/NumPart_Total_HW',this%npar_hw)
  call hdf5_write_attribute(fh,'Header/MassTable',this%mass)
  call hdf5_write_attribute(fh,'Header/Time',this%a)

  call hdf5_write_attribute(fh,'Header/Redshift',this%z)
  call hdf5_write_attribute(fh,'Header/BoxSize',this%boxlen)
  call hdf5_write_attribute(fh,'Header/NumFilesPerSnapshot',this%nfiles)
  call hdf5_write_attribute(fh,'Header/Omega0',this%OmegaM)
  call hdf5_write_attribute(fh,'Header/OmegaLambda',this%OmegaL)
  call hdf5_write_attribute(fh,'Header/HubbleParam',this%h)

  call hdf5_write_attribute(fh,'Header/Flag_Sfr',this%flag_sfr)
  call hdf5_write_attribute(fh,'Header/Flag_Cooling',this%flag_cooling)
  call hdf5_write_attribute(fh,'Header/Flag_StellarAge',this%flag_age)
  call hdf5_write_attribute(fh,'Header/Flag_Metals',this%flag_metals)
  call hdf5_write_attribute(fh,'Header/Flag_Feedback',this%flag_feedback)
  call hdf5_write_attribute(fh,'Header/Flag_Entropy_ICs',this%flag_entr_ics)
#endif

end subroutine write_header_to_lun_hdf5


!> reads the first header in a snapshot. snapbase is everything up until 
!! the snapshot number so that for /data/snapshot_009/snap_009[.0][.hdf5], 
!! snapbase = /data/snapshot_009/snap
!---------------------------------------------------------------------------
subroutine read_header_from_file1_hdf5(this, snapbase)

  character(*), intent(in) :: snapbase
  class(gadget_public_header_hdf5_type) :: this

  character(len(snapbase)+10) :: snapfile
  logical :: fthere

  snapfile = trim(snapbase) // ".hdf5"
  inquire(file=snapfile, exist=fthere)
  if (.not. fthere) then
     snapfile = trim(snapbase) // ".0.hdf5"
  endif

  call this%read_file_hdf5(snapfile)
  call this%print_lun(stdout)

end subroutine read_header_from_file1_hdf5





end module gadget_public_header_class_hdf5
