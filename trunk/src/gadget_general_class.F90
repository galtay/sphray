!> \file gadget_general_class.f90

!> \brief Handles universal GADGET stuff.  
!!
!! Provides types to handle units and constants in a 
!! standard way.  Note that where 6 element arrays are needed that
!! correspond to particle types, we have indexed them from 0-5 as 
!! opposed to the default Fortran 1-6. 
!< 

module gadget_general_class
use myf03_mod

#ifdef useMPI
  use mpi
#endif

implicit none
private


public :: ptype_names
public :: gadget_constants_type
public :: gadget_units_type
public :: form_gadget_snapshot_file_name
public :: broadcast_gadget_units



!> Particle type names
!-----------------------------------------
character(5), parameter :: ptype_names(0:5) = (/"gas  ","halo ","disk ",&
                                                "bulge","stars","bndry"/)

!> Constants type ( 23 doubles )
!-----------------------------------------
type gadget_constants_type
   real(r8b) :: PI = 3.1415927d0                 !< [pure]                  
   real(r8b) :: GAMMA = 1.6666667d0              !< [pure]
   real(r8b) :: GRAVITY = 6.6720000d-08          !< [cm^3 g^-1s^-2]    
   real(r8b) :: SOLAR_MASS = 1.9890000d+33       !< [g]                        
   real(r8b) :: SOLAR_LUM = 3.8260000d+33        !< [erg/s]
   real(r8b) :: RAD_CONST = 7.5650000d-15        !< [erg cm^-3 K^-4]  
   real(r8b) :: AVOGADRO = 6.0222000d+23         !< [pure] 
   real(r8b) :: BOLTZMANN = 1.3806000d-16        !< [cm^2 g s^-2 K^-1] 
   real(r8b) :: GAS_CONST = 83142500.d0          !< [erg K^-1 mol^-1]  
   real(r8b) :: C = 2.9979000d+10                !< [cm/s]   
   real(r8b) :: PLANCK = 6.6262000d-27           !< [cm^2 g s^-1]  
   real(r8b) :: CM_PER_MPC = 3.0856780d+24       !< [pure]    
   real(r8b) :: PROTONMASS = 1.6726000d-24       !< [g]
   real(r8b) :: ELECTRONMASS = 9.1095300d-28     !< [g] 
   real(r8b) :: ELECTRONCHARGE = 4.8032000d-10   !< [esu]
   real(r8b) :: HUBBLE = 3.2407789d-18           !< [s^-1 h]
   real(r8b) :: T_CMB0 = 2.7280000d0             !< [K]
   real(r8b) :: SEC_PER_MEGAYEAR = 3.1550000d+13 !< [pure]
   real(r8b) :: SEC_PER_YEAR = 31550000.d0       !< [pure]
   real(r8b) :: STEFAN = 7.5657000d-15     !< = rad. const. [erg cm^-3 K^-4]  
   real(r8b) :: THOMPSON = 6.6524587d-25   !< [cm^2]
   real(r8b) :: EV_TO_ERG = 1.6021765d-12  !< [pure]
   real(r8b) :: Z_SOLAR = 0.012663729d0    !< [Mass Fraction] 
end type gadget_constants_type


!> Units type (default= kpc/h, 1.0e10 Msolar/h, km/s)
!-----------------------------------------------------
type gadget_units_type
   real(r8b) :: cgs_length   = 3.0856780d21  !<  [cm h^-1]
   real(r8b) :: cgs_mass     = 1.989d43      !<  [g h^-1]
   real(r8b) :: cgs_velocity = 1.0d5         !<  [cm s^-1]
   real(r8b) :: cgs_time     = 3.085678d16   !<  [s h^-1]
   real(r8b) :: cgs_density  = 6.769911d-22  !<  [g cm^-3 h^2]
   real(r8b) :: cgs_pressure = 6.769911d-12  !<  [ba = g cm^-1 s^-2 h^2]
   real(r8b) :: cgs_energy   = 1.989d53      !<  [erg = g cm^2 s^-2 h^-1] 
 contains
   procedure :: set_user => set_gadget_units_user !< sets user supplied values
   procedure :: print_lun => print_units_to_lun   !< formatted print to lun
end type gadget_units_type



contains


!> sets user supplied units
!--------------------------------------------------------------
subroutine set_gadget_units_user(this, cgs_length, cgs_mass, cgs_velocity)
  class(gadget_units_type) :: this
  real(r8b) :: cgs_length
  real(r8b) :: cgs_mass
  real(r8b) :: cgs_velocity

  this%cgs_length   = cgs_length
  this%cgs_mass     = cgs_mass
  this%cgs_velocity = cgs_velocity

  this%cgs_density  = this%cgs_mass / this%cgs_length**3
  this%cgs_energy   = this%cgs_mass * this%cgs_velocity**2
  this%cgs_time     = this%cgs_length / this%cgs_velocity
  this%cgs_pressure = this%cgs_mass / &
       (this%cgs_length**3 / this%cgs_velocity**2)

end subroutine set_gadget_units_user


!> formatted print of units to lun (including standard out)
!---------------------------------------------------------------
subroutine print_units_to_lun(this, lun, h)
  class(gadget_units_type), intent(in) :: this 
  integer(i4b), intent(in) :: lun
  real(r8b), intent(in) :: h

  type(gadget_constants_type) :: gconst

  character(clen) :: n1
  character(clen) :: star_fmt
  character(clen) :: unit_fmt
  character(clen) :: line_fmt

  real(r8b), parameter :: cm_per_km = 1.0d5

  ! binding energy of 1.0d8 [Msolar/h] halo @ z=0 
  ! Delta_c = 18 pi^2, units = erg/h
  ! http://arxiv.org/pdf/astro-ph/0010468v3
  real(r8b), parameter :: E8 = 5.45d53 * 1.0d-1

  real(r8b) :: rho_crit_0

  rho_crit_0 = 3.0d0 * gconst%hubble**2 / (8.0d0 * gconst%pi * gconst%gravity)

  star_fmt = "(78('='))"
  line_fmt = "(78('-'))"
  unit_fmt = "(A,T30,A)"

  write(lun,*)
  write(lun,star_fmt)
  write(lun,*) "   Units"
  write(lun,line_fmt)
  write(lun,*) "  h = ", h, " = H0[km/s/Mpc] / 100 "
  write(n1,'(ES20.6)') this%cgs_length
  write(lun,unit_fmt) "length [cm/h]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_mass
  write(lun,unit_fmt) "mass [g/h]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_velocity
  write(lun,unit_fmt) "velocity [cm/s]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_density
  write(lun,unit_fmt) "density [g/cm^3 h^2]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_energy
  write(lun,unit_fmt) "energy [erg/h]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_time
  write(lun,unit_fmt) "time [s/h]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_pressure
  write(lun,unit_fmt) "pressure [ba h^2]:", trim(adjustl(n1)) 

  write(lun,*) 

  write(n1,'(ES20.6)') rho_crit_0
  write(lun,unit_fmt) "rho_crit_0 [g/cm^3 h^2]:", trim(adjustl(n1))
  write(n1,'(ES20.6)') E8
  write(lun,unit_fmt) "E8 [erg/h]:", trim(adjustl(n1))

  write(lun,*)

  write(n1,'(ES20.6)') this%cgs_length / gconst%cm_per_mpc
  write(lun,unit_fmt) "length [Mpc/h]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_mass / gconst%solar_mass
  write(lun,unit_fmt) "mass [Msolar/h]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_velocity / cm_per_km
  write(lun,unit_fmt) "velocity [km/s]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_density / rho_crit_0
  write(lun,unit_fmt) "density [rho_crit_0 h^2]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_energy / E8
  write(lun,unit_fmt) "energy [E8]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_time / gconst%sec_per_megayear
  write(lun,unit_fmt) "time [Myr/h]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_pressure / gconst%boltzmann 
  write(lun,unit_fmt) "pressure/k_b [cm^-3 K h^2]:", trim(adjustl(n1)) 

  write(lun,star_fmt)

end subroutine print_units_to_lun


!> forms a snapshot file name from, path, base, SnapNum, and FileNum.
!! For example, path = "/home/galtay/data/snapshots", base = "snap",
!! SnapNum = 2 and FileNum = 15 would return in SnapFile, 
!! "/home/galtay/data/snapshots/snap_002.15".  Setting the hdf5bool to 
!! true appends ".hdf5" to the file. 
!------------------------------------------------------------------------
subroutine form_gadget_snapshot_file_name(path, base, SnapNum, FileNum, &
     SnapFile, hdf5bool)

  character(*), intent(in)     :: path        !< path to snapshot dir
  character(*), intent(in)     :: base        !< base snapshot name
  integer(i4b), intent(in)     :: SnapNum     !< snapshot number
  integer(i4b), intent(in)     :: FileNum     !< file number of snapshot
  character(clen), intent(out) :: SnapFile    !< snapshot filename
  logical, intent(in)          :: hdf5bool    !< hdf5 file? 

  character(clen) :: SnapFileTmp
  character(10) :: FileNumChar
  character(clen) :: fmt
  logical :: Fthere

  write(FileNumChar,"(I6)") FileNum
  fmt = "(A,'/',A,'_',I3.3)"

  ! first write a file with no extension
  !--------------------------------------
  write(SnapFileTmp,fmt) trim(path), trim(base), SnapNum
  SnapFile = trim(SnapFileTmp)
  if (hdf5bool) SnapFile = trim(SnapFile) // ".hdf5"
  inquire( file=SnapFile, exist=Fthere )


  ! if the file number is 0 and a file with no extension exists then return
  ! otherwise append the FileNum to the file name
  !-------------------------------------------------------------------------
  if (FileNum == 0 .and. Fthere) then
     return
  else
     SnapFile = trim(SnapFileTmp) // "." // trim(adjustl(FileNumChar))
     if (hdf5bool) SnapFile = trim(SnapFile) // ".hdf5"
  end if


end subroutine form_gadget_snapshot_file_name



!> broadcasts units
!-----------------------------------------
subroutine broadcast_gadget_units( units )
  type(gadget_units_type) :: units
  integer :: count
  integer :: root
  integer :: ierr

#ifdef useMPI
  count = 7
  root = 0
  call mpi_bcast( units, count, mpi_double_precision, root, mpi_comm_world, ierr )
  if (ierr /= 0) then
     call mpi_barrier(mpi_comm_world, ierr)
     call mpi_finalize(ierr)
  endif
#endif

end subroutine broadcast_gadget_units






end module gadget_general_class


