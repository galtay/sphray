!> \file gadget_header_class.F90

!> \brief Handles GADGET style headers, constants, and units. 
!<

module gadget_header_class
use myf03_mod

#ifdef useHDF5
  use hdf5_wrapper
#endif

#ifdef useMPI
  use mpi
#endif

implicit none
private

public :: gadget_header_type
public :: gadget_units_type
public :: gadget_constants_type

public :: form_gadget_snapshot_file_name

public :: read_gadget_header_file_hdf5
public :: read_gadget_header_lun_hdf5
public :: write_gadget_header_lun_hdf5

public :: read_gadget_header_file
public :: read_gadget_header_lun
public :: write_gadget_header_lun

public :: read_first_header
public :: read_gadget_units
public :: read_gadget_constants
public :: gadget_header_to_file

public :: broadcast_gadget_header
public :: broadcast_gadget_units
public :: broadcast_gadget_constants





!> Gadget particle type names
!-----------------------------------------
character(5), parameter :: ptype_names(6) = (/"gas  ","halo ","disk ",&
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
   procedure :: set_owls_gimic => set_gadget_units_owls_gimic !< sets len = Mpc/h
   procedure :: set_user => set_gadget_units_user !< sets user supplied values
   procedure :: print_lun => print_units_to_lun   !< formatted print to lun
end type gadget_units_type


!> Gadget header type
!-----------------------------------------
type gadget_header_type
   integer(i4b) :: npar_file(6)    !< number of particles in snapshot file
   real(r8b)    :: mass(6)         !< mass of each particle type if constant
   real(r8b)    :: a               !< scale factor or time
   real(r8b)    :: z               !< redshift
   integer(i4b) :: flag_sfr        !< flag for star formation
   integer(i4b) :: flag_feedback   !< flag for feedback
   integer(i4b) :: npar_all(6)     !< number of particles in whole snapshot
   integer(i4b) :: flag_cooling    !< flag for radiative cooling
   integer(i4b) :: nfiles          !< number of files in a this snapshot
   real(r8b)    :: boxlen          !< box length
   real(r8b)    :: OmegaM          !< omega matter
   real(r8b)    :: OmegaL          !< omega lambda
   real(r8b)    :: h               !< little hubble
   integer(i4b) :: flag_age        !< flag for stellar age
   integer(i4b) :: flag_metals     !< flag for metallicity
   integer(i4b) :: npar_hw(6)      !< 64 bit part of npar 
   integer(i4b) :: flag_entr_ics   !< flag for entropic initial conditions
                                  
   real(r8b)    :: OmegaB          !< omega baryon 
   integer(i8b) :: rays_traced     !< number of rays traced
   integer(i4b) :: flag_Hmf        !< output has Hydorgen mass fraction? 
   integer(i4b) :: flag_Hemf       !< output has Helium mass fraction? 
   integer(i4b) :: flag_helium     !< output has xHeI and xHeII?
   integer(i4b) :: flag_gammaHI    !< output has HI photoionization rate? 
   integer(i4b) :: flag_cloudy     !< output has Cloudy xHIeq? 
   integer(i4b) :: flag_eos        !< output has EOS info?
   real(r8b)    :: time_myr        !< time since BB from cosmo variables
   real(i4b)    :: unused(3)       !< spacer
end type gadget_header_type






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


!> sets units used in OWLS/GIMIC runs
!--------------------------------------------------------------
subroutine set_gadget_units_owls_gimic(this)
  class(gadget_units_type) :: this

  this%cgs_length   = 3.0856780d24
  this%cgs_mass     = 1.989d43    
  this%cgs_velocity = 1.0d5

  this%cgs_density  = this%cgs_mass / this%cgs_length**3
  this%cgs_energy   = this%cgs_mass * this%cgs_velocity**2
  this%cgs_time     = this%cgs_length / this%cgs_velocity
  this%cgs_pressure = this%cgs_mass / &
       (this%cgs_length**3 / this%cgs_velocity**2)

end subroutine set_gadget_units_owls_gimic



!> formatted print of units to lun (including standard out)
!---------------------------------------------------------------
subroutine print_units_to_lun(this, lun)
  class(gadget_units_type), intent(in) :: this 
  integer(i4b), intent(in) :: lun

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
subroutine form_gadget_snapshot_file_name(path, base, SnapNum, FileNum, SnapFile, hdf5bool)
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



 
!> reads a gadget header from hdf5 snapfile and closes it afterwards
!---------------------------------------------------------------------
subroutine read_gadget_header_file_hdf5(snapfile, ghead)
  character(*), intent(in) :: snapfile
  type(gadget_header_type) :: ghead
  integer(i4b) :: fh

#ifdef useHDF5
  call hdf5_open_file( fh, snapfile, readonly=.true. )
  call hdf5_read_attribute(fh,'Header/NumPart_ThisFile',ghead%npar_file)
  call hdf5_read_attribute(fh,'Header/MassTable',ghead%mass)
  call hdf5_read_attribute(fh,'Header/ExpansionFactor',ghead%a)
  call hdf5_read_attribute(fh,'Header/Redshift',ghead%z)
  call hdf5_read_attribute(fh,'Header/Flag_Sfr',ghead%flag_sfr)
  call hdf5_read_attribute(fh,'Header/Flag_Feedback',ghead%flag_feedback)
  call hdf5_read_attribute(fh,'Header/NumPart_Total',ghead%npar_all)
  call hdf5_read_attribute(fh,'Header/Flag_Cooling',ghead%flag_cooling)
  call hdf5_read_attribute(fh,'Header/NumFilesPerSnapshot',ghead%nfiles)
  call hdf5_read_attribute(fh,'Header/BoxSize',ghead%boxlen)
  call hdf5_read_attribute(fh,'Header/Omega0',ghead%OmegaM)
  call hdf5_read_attribute(fh,'Header/OmegaLambda',ghead%OmegaL)
  call hdf5_read_attribute(fh,'Header/HubbleParam',ghead%h)
  call hdf5_read_attribute(fh,'Header/Flag_StellarAge',ghead%flag_age)
  call hdf5_read_attribute(fh,'Header/Flag_Metals',ghead%flag_metals)
  call hdf5_read_attribute(fh,'Header/NumPart_Total_HighWord',ghead%npar_hw)
  call hdf5_read_attribute(fh,'Header/OmegaBaryon',ghead%OmegaB)
  call hdf5_close_file( fh )
#else
  write(*,*) 'cant call read_gadget_header_file_hdf5 w/o defining hdf5 macro'
  stop
#endif

end subroutine read_gadget_header_file_hdf5



!> reads a gadget header from an already open file associated with lun 
!----------------------------------------------------------------------
subroutine read_gadget_header_lun_hdf5(lun, ghead)
  integer(i4b), intent(in) :: lun    
  type(gadget_header_type) :: ghead

#ifdef useHDF5
  call hdf5_read_attribute(lun,'Header/NumPart_ThisFile',ghead%npar_file)
  call hdf5_read_attribute(lun,'Header/MassTable',ghead%mass)
  call hdf5_read_attribute(lun,'Header/ExpansionFactor',ghead%a)
  call hdf5_read_attribute(lun,'Header/Redshift',ghead%z)
  call hdf5_read_attribute(lun,'Header/Flag_Sfr',ghead%flag_sfr)
  call hdf5_read_attribute(lun,'Header/Flag_Feedback',ghead%flag_feedback)
  call hdf5_read_attribute(lun,'Header/NumPart_Total',ghead%npar_all)
  call hdf5_read_attribute(lun,'Header/Flag_Cooling',ghead%flag_cooling)
  call hdf5_read_attribute(lun,'Header/NumFilesPerSnapshot',ghead%nfiles)
  call hdf5_read_attribute(lun,'Header/BoxSize',ghead%boxlen)
  call hdf5_read_attribute(lun,'Header/Omega0',ghead%OmegaM)
  call hdf5_read_attribute(lun,'Header/OmegaLambda',ghead%OmegaL)
  call hdf5_read_attribute(lun,'Header/HubbleParam',ghead%h)
  call hdf5_read_attribute(lun,'Header/Flag_StellarAge',ghead%flag_age)
  call hdf5_read_attribute(lun,'Header/Flag_Metals',ghead%flag_metals)
  call hdf5_read_attribute(lun,'Header/NumPart_Total_HighWord',ghead%npar_hw)
  call hdf5_read_attribute(lun,'Header/OmegaBaryon',ghead%OmegaB)
#else
  write(*,*) 'cant call read_gadget_header_lun_hdf5 w/o defining hdf5 macro'
  stop
#endif

end subroutine read_gadget_header_lun_hdf5



!> writes a gadget header from an already open file associated with lun 
!----------------------------------------------------------------------
subroutine write_gadget_header_lun_hdf5(lun, ghead)
  integer(i4b), intent(in) :: lun    
  type(gadget_header_type) :: ghead

#ifdef useHDF5
  call hdf5_write_attribute(lun,'Header/NumPart_ThisFile',ghead%npar_file)
  call hdf5_write_attribute(lun,'Header/MassTable',ghead%mass)
  call hdf5_write_attribute(lun,'Header/ExpansionFactor',ghead%a)
  call hdf5_write_attribute(lun,'Header/Redshift',ghead%z)
  call hdf5_write_attribute(lun,'Header/Flag_Sfr',ghead%flag_sfr)
  call hdf5_write_attribute(lun,'Header/Flag_Feedback',ghead%flag_feedback)
  call hdf5_write_attribute(lun,'Header/NumPart_Total',ghead%npar_all)
  call hdf5_write_attribute(lun,'Header/Flag_Cooling',ghead%flag_cooling)
  call hdf5_write_attribute(lun,'Header/NumFilesPerSnapshot',ghead%nfiles)
  call hdf5_write_attribute(lun,'Header/BoxSize',ghead%boxlen)
  call hdf5_write_attribute(lun,'Header/Omega0',ghead%OmegaM)
  call hdf5_write_attribute(lun,'Header/OmegaLambda',ghead%OmegaL)
  call hdf5_write_attribute(lun,'Header/HubbleParam',ghead%h)
  call hdf5_write_attribute(lun,'Header/Flag_StellarAge',ghead%flag_age)
  call hdf5_write_attribute(lun,'Header/Flag_Metals',ghead%flag_metals)
  call hdf5_write_attribute(lun,'Header/NumPart_Total_HighWord',ghead%npar_hw)
  call hdf5_write_attribute(lun,'Header/OmegaBaryon',ghead%OmegaB)
#else
  write(*,*) 'cant call write_gadget_header_lun_hdf5 w/o defining hdf5 macro'
  stop
#endif

end subroutine write_gadget_header_lun_hdf5



!> reads a gadget header from snapfile and closes it afterwards
!--------------------------------------------------------------
subroutine read_gadget_header_file(snapfile, ghead)
  character(*), intent(in) :: snapfile
  type(gadget_header_type) :: ghead
  integer(i4b) :: lun    

  call open_unformatted_file_r( snapfile, lun )
  read(lun) ghead%npar_file(:), ghead%mass(:), ghead%a, ghead%z, &              
       ghead%flag_sfr, ghead%flag_feedback, ghead%npar_all(:), &    
       ghead%flag_cooling, ghead%nfiles, ghead%boxlen, ghead%OmegaM, &         
       ghead%OmegaL, ghead%h, ghead%flag_age, ghead%flag_metals, &    
       ghead%npar_hw(:), ghead%flag_entr_ics, ghead%OmegaB, &         
       ghead%rays_traced, ghead%flag_Hmf, ghead%flag_Hemf, &
       ghead%flag_helium, ghead%flag_gammaHI, ghead%flag_cloudy, &
       ghead%flag_eos, ghead%unused(:)      
  close(lun)


end subroutine read_gadget_header_file



!> reads a gadget header from an already open file associated with lun 
!----------------------------------------------------------------------
subroutine read_gadget_header_lun(lun, ghead)
  integer(i4b), intent(in) :: lun    
  type(gadget_header_type) :: ghead

  read(lun) ghead%npar_file(:), ghead%mass(:), ghead%a, ghead%z, &              
       ghead%flag_sfr, ghead%flag_feedback, ghead%npar_all(:), &    
       ghead%flag_cooling, ghead%nfiles, ghead%boxlen, ghead%OmegaM, &         
       ghead%OmegaL, ghead%h, ghead%flag_age, ghead%flag_metals, &    
       ghead%npar_hw(:), ghead%flag_entr_ics, ghead%OmegaB, &         
       ghead%rays_traced, ghead%flag_Hmf, ghead%flag_Hemf, &
       ghead%flag_helium, ghead%flag_gammaHI, ghead%flag_cloudy, &
       ghead%flag_eos, ghead%unused(:)      

end subroutine read_gadget_header_lun




!> writes a gadget header from an already open file associated with lun 
!----------------------------------------------------------------------
subroutine write_gadget_header_lun(lun, ghead)
  integer(i4b), intent(in) :: lun    
  type(gadget_header_type) :: ghead

  write(lun) ghead%npar_file(:), ghead%mass(:), ghead%a, ghead%z, &              
       ghead%flag_sfr, ghead%flag_feedback, ghead%npar_all(:), &    
       ghead%flag_cooling, ghead%nfiles, ghead%boxlen, ghead%OmegaM, &         
       ghead%OmegaL, ghead%h, ghead%flag_age, ghead%flag_metals, &    
       ghead%npar_hw(:), ghead%flag_entr_ics, ghead%OmegaB, &         
       ghead%rays_traced, ghead%flag_Hmf, ghead%flag_Hemf, &
       ghead%flag_helium, ghead%flag_gammaHI, ghead%flag_cloudy, &
       ghead%flag_eos, ghead%unused(:)      


end subroutine write_gadget_header_lun












!> reads the first header in a snapshot. snapbase is everything up until the 
!! snapshot number so that for /data/snapshot_009/snap_009[.0][.hdf5], 
!! snapbase = /data/snapshot_009/snap
!---------------------------------------------------------------------------
subroutine read_first_header( snapbase, ghead, gunits, gconst, hdf5bool, simname  )
  character(*), intent(in) :: snapbase
  type(gadget_header_type) :: ghead
  type(gadget_units_type) :: gunits
  type(gadget_constants_type) :: gconst
  logical, intent(in) :: hdf5bool
  character(*), intent(in), optional :: simname

  character(len(snapbase)+10) :: snapfile
  logical :: fthere

  if (hdf5bool) then
     snapfile = trim(snapbase) // ".hdf5"
     inquire(file=snapfile, exist=fthere)
     if (.not. fthere) then
        snapfile = trim(snapbase) // ".0.hdf5"
     endif
     call read_gadget_header_file_hdf5(snapfile, ghead)
     call read_gadget_units(snapfile, gunits, hdf5bool) 
     call read_gadget_constants(snapfile, gconst, hdf5bool)
  else
     snapfile = trim(snapbase)
     inquire(file=snapfile, exist=fthere)
     if (.not. fthere) then
        snapfile = trim(snapbase) // ".0"
     endif
     call read_gadget_header_file(snapfile, ghead)
     call read_gadget_units(snapfile, gunits, hdf5bool, simname ) 
     call read_gadget_constants(snapfile, gconst, hdf5bool, simname )
  endif


  call gadget_header_to_file(ghead,stdout) 
 

end subroutine read_first_header






!> reads gadget units from an HDF5 file
!-----------------------------------------
subroutine read_gadget_units_hdf5(snapfile, gunits)
  character(*), intent(in) :: snapfile
  type(gadget_units_type) :: gunits
  integer(i4b) :: fh
  
  if (hdf5bool) then 
#ifdef useHDF5
     call hdf5_open_file( fh, snapfile, readonly=.true. )
     call hdf5_read_attribute(fh,'Units/UnitLength_in_cm',gunits%len)
     call hdf5_read_attribute(fh,'Units/UnitMass_in_g',gunits%mass)
     call hdf5_read_attribute(fh,'Units/UnitVelocity_in_cm_per_s',gunits%vel)
     call hdf5_read_attribute(fh,'Units/UnitDensity_in_cgs',gunits%rho)
     call hdf5_read_attribute(fh,'Units/UnitEnergy_in_cgs',gunits%energy)
     call hdf5_read_attribute(fh,'Units/UnitPressure_in_cgs',gunits%prs)
     call hdf5_read_attribute(fh,'Units/UnitTime_in_s',gunits%time)
     call hdf5_close_file( fh )
#endif
  else
     if ( .not. present(simname) ) then
        write(*,*) " if not using hdf5 must supply simname [e.g. 'gimic']"
        stop
     else
        select case (trim(simname))
        case("gimic")
           call set_gadget_owls_gimic_units(gunits)
        case("owls")
           call set_gadget_owls_gimic_units(gunits)
        case("default")
           call set_gadget_default_units(gunits)
        case default
           write(*,*) "  sim name: ", trim(simname), " not recognized"
           stop
        end select
     endif
  endif

end subroutine read_gadget_units








!> reads an hdf5 gadget constants
!-----------------------------------------
subroutine read_gadget_constants(snapfile, gconst, hdf5bool, simname)
  character(*), intent(in) :: snapfile
  type(gadget_constants_type) :: gconst
  logical, intent(in) :: hdf5bool
  character(*), intent(in), optional :: simname

  integer(i4b) :: fh
  integer(i8b) :: lun

  if (hdf5bool) then 
#ifdef useHDF5
     call hdf5_open_file( fh, snapfile, readonly=.true. )
     call hdf5_read_attribute(fh,'Constants/PI',gconst%pi)
     call hdf5_read_attribute(fh,'Constants/GAMMA',gconst%gamma)
     call hdf5_read_attribute(fh,'Constants/GRAVITY',gconst%gravity)
     call hdf5_read_attribute(fh,'Constants/SOLAR_MASS',gconst%solar_mass)
     call hdf5_read_attribute(fh,'Constants/SOLAR_LUM',gconst%solar_lum)
     call hdf5_read_attribute(fh,'Constants/RAD_CONST',gconst%rad_const)
     call hdf5_read_attribute(fh,'Constants/AVOGADRO',gconst%avogadro)
     call hdf5_read_attribute(fh,'Constants/BOLTZMANN',gconst%boltzmann)
     call hdf5_read_attribute(fh,'Constants/GAS_CONST',gconst%gas_const)
     call hdf5_read_attribute(fh,'Constants/C',gconst%c)
     call hdf5_read_attribute(fh,'Constants/PLANCK',gconst%planck)
     call hdf5_read_attribute(fh,'Constants/CM_PER_MPC',gconst%cm_per_mpc)
     call hdf5_read_attribute(fh,'Constants/PROTONMASS',gconst%protonmass)
     call hdf5_read_attribute(fh,'Constants/ELECTRONMASS',gconst%electronmass)
     call hdf5_read_attribute(fh,'Constants/ELECTRONCHARGE',gconst%electroncharge)
     call hdf5_read_attribute(fh,'Constants/HUBBLE',gconst%hubble)
     call hdf5_read_attribute(fh,'Constants/T_CMB0',gconst%t_cmb0)
     call hdf5_read_attribute(fh,'Constants/SEC_PER_MEGAYEAR',gconst%sec_per_megayear)
     call hdf5_read_attribute(fh,'Constants/SEC_PER_YEAR',gconst%sec_per_year)
     call hdf5_read_attribute(fh,'Constants/STEFAN',gconst%stefan)
     call hdf5_read_attribute(fh,'Constants/THOMPSON',gconst%thompson)
     call hdf5_read_attribute(fh,'Constants/EV_TO_ERG',gconst%ev_to_erg)
     call hdf5_read_attribute(fh,'Constants/Z_Solar',gconst%z_solar)
     call hdf5_close_file( fh )
#endif
  else
     if ( .not. present(simname) ) then
        write(*,*) " if not using hdf5 must supply simname [e.g. 'gimic']"
        stop
     else
        select case (trim(simname))
        case("gimic")
           call set_gadget_owls_gimic_constants(gconst)
        case("owls")
           call set_gadget_owls_gimic_constants(gconst)
        case("default")
           call set_gadget_owls_gimic_constants(gconst)
        case default
           write(*,*) "  sim name: ", trim(simname), " not recognized"
           stop
        end select
     endif
  endif
  
end subroutine read_gadget_constants







!> writes gadget header data to the screen
!-----------------------------------------
subroutine gadget_header_to_file(ghead,lun)
  type(gadget_header_type), intent(in) :: ghead !< particle header to print
  integer(i4b), intent(in) :: lun
  integer(i8b) :: i


98  format (A,T12,I12,T36,I12)
99  format (T17,A,T30,A,T40,A)
100 format (A,I1,A,A,A,I12,ES12.3,I12)
101 format (A,F8.5,A,F10.5,A,F8.4)  
102 format (A,F8.5,A,F8.5,A,F13.4)  
104 format (A,I3)
103 format (A,6(A,I1))

  write(lun,*)
  write(lun,99) "n_all","mass","n_file"
  do i = 1,6
     write(lun,100) "type",i,"(",ptype_names(i),")",&
          ghead%npar_all(i),ghead%mass(i),ghead%npar_file(i)
  end do
  write(lun,98) "total:", sum(ghead%npar_all), sum(ghead%npar_file)

  write(lun,101) "a   = ", ghead%a, &
       ", z    = ", ghead%z, &
       ", h       = ", ghead%h
  write(lun,102) "OmegaM = ", ghead%OmegaM, &
       ", OmegaL = ", ghead%OmegaL, & 
       ", boxsize = ", ghead%boxlen 
  write(lun,104) "nfiles  = ", ghead%nfiles
  write(lun,103) "/flags/ ",&
       "  sfr=",ghead%flag_sfr, &
       ", feedback=",ghead%flag_feedback, &
       ", cooling=",ghead%flag_cooling, &
       ", age=",ghead%flag_age, &
       ", metals=",ghead%flag_metals, &
       ", entr_ics=", ghead%flag_entr_ics
  write(lun,*) 


end subroutine gadget_header_to_file





!> broadcasts header
!-----------------------------------------
subroutine broadcast_gadget_header( head )
  type(gadget_header_type) :: head
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



!> broadcasts constants
!-----------------------------------------
subroutine broadcast_gadget_constants( const )
  type(gadget_constants_type) :: const
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

end module gadget_header_class
