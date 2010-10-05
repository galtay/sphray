!> \file gadget_public_header_class.f90

!> \brief Handles GADGET 2.0 Public version style headers.  
!!
!! Contains the means to read/write headers from a file or open lun.  
!! In addition, provides types to handle units and constants in a 
!! standard way.  Note that where 6 element arrays are needed that
!! correspond to particle types, we have indexed them from 0-5 as 
!! opposed to the default Fortran 1-6. 
!< 

module gadget_public_header_class
use myf03_mod
#ifdef useMPI
  use mpi
#endif

implicit none
private

public :: gadget_public_constants_type
public :: gadget_public_units_type
public :: gadget_public_header_type

public :: form_gadget_snapshot_file_name
public :: ptype_names

public :: broadcast_gadget_header
public :: broadcast_gadget_units
public :: broadcast_gadget_constants


!> Particle type names
!-----------------------------------------
character(5), parameter :: ptype_names(0:5) = (/"gas  ","halo ","disk ",&
                                                "bulge","stars","bndry"/)


!> Constants type ( 23 doubles )
!-----------------------------------------
type gadget_public_constants_type
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
end type gadget_public_constants_type



!> Units type (default= kpc/h, 1.0e10 Msolar/h, km/s)
!-----------------------------------------------------
type gadget_public_units_type
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
end type gadget_public_units_type



!> Header type
!-----------------------------------------
type gadget_public_header_type
   integer(i4b) :: npar_file(0:5)   !< number of particles in snapshot file
   real(r8b)    :: mass(0:5)        !< mass of each particle type if constant
   real(r8b)    :: a                !< scale factor or time
   real(r8b)    :: z                !< redshift
   integer(i4b) :: flag_sfr         !< flag for star formation
   integer(i4b) :: flag_feedback    !< flag for feedback
   integer(i4b) :: npar_all(0:5)    !< number of particles in whole snapshot
   integer(i4b) :: flag_cooling     !< flag for radiative cooling
   integer(i4b) :: nfiles           !< number of files in a this snapshot
   real(r8b)    :: boxlen           !< box length
   real(r8b)    :: OmegaM           !< omega matter
   real(r8b)    :: OmegaL           !< omega lambda
   real(r8b)    :: h                !< little hubble
   integer(i4b) :: flag_age         !< flag for stellar age
   integer(i4b) :: flag_metals      !< flag for metallicity
   integer(i4b) :: npar_hw(0:5)     !< 64 bit part of npar 
   integer(i4b) :: flag_entr_ics    !< ICs contain entropy instead of energy
   integer(i4b) :: unused(15)       !< padding to make 256 bytes
 contains
   procedure :: read_file  => read_header_from_file  !< read from unopened file
   procedure :: read_file1 => read_header_from_file1 !< read 1st file in snap
   procedure :: read_lun   => read_header_from_lun   !< read from open file
   procedure :: write_lun  => write_header_to_lun    !< write to open file
   procedure :: print_lun  => print_header_to_lun    !< formatted print to lun
   procedure :: return_gyr     
end type gadget_public_header_type



contains


!> uses the scale factor and OmegaM+OmegaL to 
!! calculate the time in Gyr since the Big Bang
!----------------------------------------------------
function return_gyr(this) result(t_gyr)
  class(gadget_public_header_type) :: this
  type(gadget_public_constants_type) :: const
  real(r8b) :: t_gyr
  real(r8b) :: aeq, pre, arg, H0

  H0 = this%h * 100.0d0 / ( const%CM_PER_MPC / 1.0d5 )
  aeq = (this%OmegaM/this%OmegaL)**(1.0d0/3.0d0)
  pre = 2.0d0 / (3.0d0 * sqrt(this%OmegaL) )
  arg = (this%a/aeq)**(3.0d0/2.0d0) + sqrt(1 + (this%a/aeq)**3 )
  t_gyr = pre * log(arg) / H0  ! in s
  t_gyr = t_gyr / ( const%SEC_PER_MEGAYEAR * 1.0d3 )

end function return_gyr



!> sets user supplied units
!--------------------------------------------------------------
subroutine set_gadget_units_user(this, cgs_length, cgs_mass, cgs_velocity)
  class(gadget_public_units_type) :: this
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
  class(gadget_public_units_type), intent(in) :: this 
  integer(i4b), intent(in) :: lun
  real(r8b), intent(in) :: h

  type(gadget_public_constants_type) :: gconst

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



!> forms a snapshot file name from,
!!  path [e.g. "/home/galtay/data/snapshots"],
!!  base [e.g. "snap"],
!!  snapshot number [e.g. 2], and
!!  file number [e.g. 15].  
!!  The above choices would return,
!!  "/home/galtay/data/snapshots/snap_002.15
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



 


!> reads a gadget header from an open raw binary file 
!--------------------------------------------------------------
subroutine read_header_from_lun(this, lun)
  class(gadget_public_header_type) :: this
  integer(i4b), intent(in) :: lun

  read(lun) this%npar_file(:), this%mass(:), this%a, this%z, &              
       this%flag_sfr, this%flag_feedback, this%npar_all(:), &    
       this%flag_cooling, this%nfiles, this%boxlen, this%OmegaM, &         
       this%OmegaL, this%h, this%flag_age, this%flag_metals, &    
       this%npar_hw(:), this%flag_entr_ics, this%unused(:)      

end subroutine read_header_from_lun



!> reads a gadget header from snapfile and closes it afterwards
!--------------------------------------------------------------
subroutine read_header_from_file(this, snapfile)
  class(gadget_public_header_type) :: this
  character(*), intent(in) :: snapfile
  integer(i4b) :: lun    

  call open_unformatted_file_r( snapfile, lun )
  call this%read_lun(lun)
  close(lun)

end subroutine read_header_from_file


!> writes a gadget header to an open file
!--------------------------------------------------------------
subroutine write_header_to_lun(this, lun)
  class(gadget_public_header_type) :: this
  integer(i4b), intent(in) :: lun

  write(lun) this%npar_file(:), this%mass(:), this%a, this%z, &              
       this%flag_sfr, this%flag_feedback, this%npar_all(:), &    
       this%flag_cooling, this%nfiles, this%boxlen, this%OmegaM, &         
       this%OmegaL, this%h, this%flag_age, this%flag_metals, &    
       this%npar_hw(:), this%flag_entr_ics, this%unused(:)      

end subroutine write_header_to_lun



!> reads the first header in a snapshot. snapbase is everything up until 
!! the snapshot number so that for /data/snapshot_009/snap_009[.0][.hdf5], 
!! snapbase = /data/snapshot_009/snap
!---------------------------------------------------------------------------
subroutine read_header_from_file1(this, snapbase)

  character(*), intent(in) :: snapbase
  class(gadget_public_header_type) :: this

  character(len(snapbase)+10) :: snapfile
  logical :: fthere

  snapfile = trim(snapbase)
  inquire(file=snapfile, exist=fthere)
  if (.not. fthere) then
     snapfile = trim(snapbase) // ".0"
  endif
  
  call this%read_file(snapfile)
  call this%print_lun(stdout)

end subroutine read_header_from_file1




!> formatted print of header to lun (including standard out)
!---------------------------------------------------------------
subroutine print_header_to_lun(this, lun)
  class(gadget_public_header_type), intent(in) :: this 
  integer(i4b), intent(in) :: lun
  integer(i8b) :: i

  character(clen) :: n1,n2,n3
  character(clen) :: clmn_fmt
  character(clen) :: type_fmt
  character(clen) :: totl_fmt
  character(clen) :: scal_fmt
  character(clen) :: flag_fmt
  character(clen) :: star_fmt
  character(clen) :: line_fmt

  clmn_fmt = "(T17,A,T32,A,T47,A)"
  type_fmt = "(A,I1,A,A,A,T17,A,T32,A,T47,A)"
  totl_fmt = "(A,T17,A,T47,A)"
  scal_fmt = "(A,A,T25,A,A,T45,A,A)"
  flag_fmt = "(A,6(A,I1))"
  star_fmt = "(78('='))"
  line_fmt = "(78('-'))"

  write(lun,*)
  write(lun,star_fmt)
  write(lun,clmn_fmt) "n_all", "mass", "n_file"
  write(lun,line_fmt)
  do i = 0,5
     write(n1,'(I20)') this%npar_all(i)
     write(n2,'(E12.5)') this%mass(i)
     write(n3,'(I20)') this%npar_file(i)
     write(lun,type_fmt) "type",i,"(",ptype_names(i),")",&
          trim(adjustl(n1)), trim(adjustl(n2)), trim(adjustl(n3))
  end do
  write(lun,line_fmt)
  write(n1,'(I20)') sum(this%npar_all)
  write(n2,'(I20)') sum(this%npar_file)
  write(lun,totl_fmt) "total:", trim(adjustl(n1)), trim(adjustl(n2)) 
  write(lun,*)

  write(n1,'(F12.5)') this%a
  write(n2,'(F12.5)') this%z
  write(n3,'(F12.5)') this%h

  write(lun,scal_fmt) "a = ", trim(adjustl(n1)), &
       "z = ", trim(adjustl(n2)), &
       "h = ", trim(adjustl(n3))

  write(n1,'(F12.5)') this%OmegaM
  write(n2,'(F12.5)') this%OmegaL
  write(n3,'(F12.5)') this%boxlen

  write(lun,scal_fmt) "OmegaM = ", trim(adjustl(n1)), &
       "OmegaL = ", trim(adjustl(n2)), &
       "BoxSize = ", trim(adjustl(n3))
  write(lun,*)

  write(n1,'(I20)') this%nfiles
  write(lun,'(A,A)') "nfiles = ", trim(adjustl(n1))
  write(lun,flag_fmt) "/flags/ ",&
       "  sfr=",this%flag_sfr, &
       ", feedback=",this%flag_feedback, &
       ", cooling=",this%flag_cooling, &
       ", age=",this%flag_age, &
       ", metals=",this%flag_metals, &
       ", entr_ics=", this%flag_entr_ics
  write(lun,*) 
  write(lun,star_fmt)

end subroutine print_header_to_lun






!> broadcasts header
!-----------------------------------------
subroutine broadcast_gadget_header( head )
  type(gadget_public_header_type) :: head
  integer :: count
  integer :: root
  integer :: ierr

  root = 0

#ifdef useMPI
  count = 6
  call mpi_bcast( head%npar_file, count, mpi_integer, root, mpi_comm_world, ierr )
  call mpi_bcast( head%mass, count, mpi_double_precision, root, mpi_comm_world, ierr )
  call mpi_bcast( head%npar_all, count, mpi_integer, root, mpi_comm_world, ierr )
  call mpi_bcast( head%npar_hw, count, mpi_integer, root, mpi_comm_world, ierr )

  count = 1
  call mpi_bcast( head%a, count, mpi_double_precision, root, mpi_comm_world, ierr )
  call mpi_bcast( head%z, count, mpi_double_precision, root, mpi_comm_world, ierr )
  call mpi_bcast( head%flag_sfr, count, mpi_integer, root, mpi_comm_world, ierr )
  call mpi_bcast( head%flag_feedback, count, mpi_integer, root, mpi_comm_world, ierr )
  call mpi_bcast( head%flag_cooling, count, mpi_integer, root, mpi_comm_world, ierr )
  call mpi_bcast( head%nfiles, count, mpi_integer, root, mpi_comm_world, ierr )
  call mpi_bcast( head%boxlen, count, mpi_double_precision, root, mpi_comm_world, ierr )
  call mpi_bcast( head%OmegaM, count, mpi_double_precision, root, mpi_comm_world, ierr )
  call mpi_bcast( head%OmegaL, count, mpi_double_precision, root, mpi_comm_world, ierr )
  call mpi_bcast( head%h, count, mpi_double_precision, root, mpi_comm_world, ierr )
  call mpi_bcast( head%flag_age, count, mpi_integer, root, mpi_comm_world, ierr )
  call mpi_bcast( head%flag_metals, count, mpi_integer, root, mpi_comm_world, ierr )
  call mpi_bcast( head%flag_entr_ics, count, mpi_integer, root, mpi_comm_world, ierr )
  call mpi_bcast( head%OmegaB, count, mpi_double_precision, root, mpi_comm_world, ierr )
#endif

end subroutine broadcast_gadget_header


!> broadcasts units
!-----------------------------------------
subroutine broadcast_gadget_units( units )
  type(gadget_public_units_type) :: units
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



!> broadcasts constants
!-----------------------------------------
subroutine broadcast_gadget_constants( const )
  type(gadget_public_constants_type) :: const
  integer :: count
  integer :: root
  integer :: ierr

#ifdef useMPI
  count = 23
  root = 0
  call mpi_bcast( const, count, mpi_double_precision, root, mpi_comm_world, ierr )
  if (ierr /= 0) then
     call mpi_barrier(mpi_comm_world, ierr)
     call mpi_finalize(ierr)
  endif
#endif

end subroutine broadcast_gadget_constants

end module gadget_public_header_class
