!> \file ionization_tables.F90

!> \brief the module that handles cloudy ionization tables 

module ionization_tables_mod
use myf90_mod
#ifdef hdf5
use hdf5_wrapper
#endif
implicit none

real(r8b), allocatable :: H1_table(:,:,:)
real(r8b), allocatable :: He1_table(:,:,:)
real(r8b), allocatable :: He2_table(:,:,:)


contains


subroutine read_cloudy_ionization_table(ion)
  character(5), intent(in) :: ion
  character(clen) :: filename
  integer(i4b) :: fh

  character(clen) :: str

  filename = "../data/ionization_tables/" // trim(ion) // ".hdf5"

  call hdf5_open_file(fh,filename)
  call hdf5_read_attribute(fh,'Header/cloudy_version',str)
  call hdf5_close_file(fh)

  stop "in read_cloudy_ionization_table"

end subroutine read_cloudy_ionization_table




end module ionization_tables_mod
