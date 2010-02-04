!> \file gadget_header_class.F90

!> \brief Handles GADGET style headers, constants, and units. 
!<

module gadget_header_class
use myf90_mod

#ifdef hdf5
  use hdf5_wrapper
#endif

#ifdef usempi
  use mpi
#endif

implicit none
private

public :: gadget_header_type
public :: gadget_units_type
public :: gadget_constants_type

public :: form_gadget_snapshot_file_name
public :: read_gadget_header_file
public :: read_gadget_header_lun

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

!> Gadget header type
!-----------------------------------------
type gadget_header_type
   sequence
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
   real(i4b)    :: unused(5)       !< spacer
end type gadget_header_type


!> gadget units type ( 7 doubles )
!-----------------------------------------
type gadget_units_type
   sequence
   real(r8b) :: len      
   real(r8b) :: mass
   real(r8b) :: vel
   real(r8b) :: rho
   real(r8b) :: energy
   real(r8b) :: prs
   real(r8b) :: time
end type gadget_units_type


!> gadget constants type ( 23 doubles )
!-----------------------------------------
type gadget_constants_type
   sequence
   real(r8b) :: PI              
   real(r8b) :: GAMMA           
   real(r8b) :: GRAVITY         
   real(r8b) :: SOLAR_MASS      
   real(r8b) :: SOLAR_LUM       
   real(r8b) :: RAD_CONST       
   real(r8b) :: AVOGADRO        
   real(r8b) :: BOLTZMANN       
   real(r8b) :: GAS_CONST       
   real(r8b) :: C               
   real(r8b) :: PLANCK          
   real(r8b) :: PROTONMASS      
   real(r8b) :: ELECTRONMASS    
   real(r8b) :: ELECTRONCHARGE  
   real(r8b) :: HUBBLE          
   real(r8b) :: T_CMB0          
   real(r8b) :: SEC_PER_MEGAYEAR                 
   real(r8b) :: SEC_PER_YEAR    
   real(r8b) :: STEFAN          
   real(r8b) :: THOMPSON        
   real(r8b) :: EV_TO_ERG       
   real(r8b) :: Z_SOLAR         
   real(r8b) :: CM_PER_MPC        
end type gadget_constants_type



contains


!> forms a snapshot file name from,
!!  path [e.g. "/home/galtay/data/snapshots"],
!!  base [e.g. "snap"],
!!  snapshot number [e.g. 2], and
!!  file number [e.g. 15].  
!!  The above choices would return,
!!  "/home/galtay/data/snapshots/snap_002.15
!------------------------------------------------------------------------
subroutine form_gadget_snapshot_file_name(path,base,SnapNum,FileNum,SnapFile)
  character(*), intent(in)     :: path        !< path to snapshot dir
  character(*), intent(in)     :: base        !< base snapshot name
  integer(i8b), intent(in)     :: SnapNum     !< snapshot number
  integer(i8b), intent(in)     :: FileNum     !< file number of snapshot
  character(clen), intent(out) :: SnapFile    !< snapshot filename

  character(10) :: FileNumChar
  character(clen) :: fmt
  logical :: Fthere

  write(FileNumChar,"(I6)") FileNum
  fmt = "(A,'/',A,'_',I3.3)"

  ! first write a file with no extension
  !--------------------------------------
  write(SnapFile,fmt) trim(path), trim(base), SnapNum
  inquire( file=SnapFile, exist=Fthere )

  ! if the file number is 0 and a file with no extension exists then return
  ! otherwise append the FileNum to the file name
  !-------------------------------------------------------------------------
  if (FileNum == 0 .and. Fthere) then
     return
  else
     SnapFile = trim(SnapFile) // "." // trim(adjustl(FileNumChar))
  end if

  ! if we are using the hdf5 lib then also append .hdf5
  !------------------------------------------------------
#ifdef hdf5
  SnapFile = trim(SnapFile) // ".hdf5"
#endif


end subroutine form_gadget_snapshot_file_name




!> reads a gadget header from snapfile and closes it afterwards
!--------------------------------------------------------------
subroutine read_gadget_header_file(snapfile, ghead)
  character(*), intent(in) :: snapfile
  type(gadget_header_type) :: ghead
  integer(i8b) :: lun    
  integer(i4b) :: fh

#ifdef hdf5
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
  call open_unformatted_file_r( snapfile, lun )
  read(lun) ghead
  close(lun)
#endif

end subroutine read_gadget_header_file



!> reads a gadget header from an already open file associated with lun 
!----------------------------------------------------------------------
subroutine read_gadget_header_lun(lun, ghead)
  integer(i4b), intent(in) :: lun    
  type(gadget_header_type) :: ghead

#ifdef hdf5
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
  read(lun) ghead
#endif

end subroutine read_gadget_header_lun



!> reads the first header in a snapshot
!-----------------------------------------
subroutine read_first_header( snapbase, ghead, gunits, gconst, simname )
  character(*) :: snapbase
  type(gadget_header_type) :: ghead
  type(gadget_units_type) :: gunits
  type(gadget_constants_type) :: gconst
  character(*), intent(in), optional :: simname
  character(len(snapbase)+10) :: snapfile
  logical :: fthere

#ifdef hdf5
  snapfile = trim(snapbase) // ".hdf5"
  inquire(file=snapfile, exist=fthere)
  if (.not. fthere) then
     snapfile = trim(snapbase) // ".0.hdf5"
  endif
#else
  snapfile = trim(snapbase)
  inquire(file=snapfile, exist=fthere)
  if (.not. fthere) then
     snapfile = trim(snapbase) // ".0"
  endif
#endif

  call read_gadget_header_file(snapfile, ghead)
  if (present(simname)) then
     call read_gadget_units(snapfile, gunits, simname) 
     call read_gadget_constants(snapfile, gconst, simname)
  else
     call read_gadget_units(snapfile, gunits) 
     call read_gadget_constants(snapfile, gconst)
  endif
  call gadget_header_to_file(ghead,stdout) 
 

end subroutine read_first_header






!> reads an hdf5 gadget units
!-----------------------------------------
subroutine read_gadget_units(snapfile, gunits, simname)
  character(*), intent(in) :: snapfile
  type(gadget_units_type) :: gunits
  character(*), intent(in), optional :: simname
  integer(i4b) :: fh
  integer(i8b) :: lun


#ifdef hdf5
  call hdf5_open_file( fh, snapfile, readonly=.true. )
  call hdf5_read_attribute(fh,'Units/UnitLength_in_cm',gunits%len)
  call hdf5_read_attribute(fh,'Units/UnitMass_in_g',gunits%mass)
  call hdf5_read_attribute(fh,'Units/UnitVelocity_in_cm_per_s',gunits%vel)
  call hdf5_read_attribute(fh,'Units/UnitDensity_in_cgs',gunits%rho)
  call hdf5_read_attribute(fh,'Units/UnitEnergy_in_cgs',gunits%energy)
  call hdf5_read_attribute(fh,'Units/UnitPressure_in_cgs',gunits%prs)
  call hdf5_read_attribute(fh,'Units/UnitTime_in_s',gunits%time)
  call hdf5_close_file( fh )
#else
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
#endif

end subroutine read_gadget_units


!> sets the gadget units used in the OWLS and GIMIC runs
!---------------------------------------------------------
subroutine set_gadget_owls_gimic_units( gunits )
  type(gadget_units_type) :: gunits
 
  gunits%len = 3.0856780d+24
  gunits%mass = 1.9890000d+43
  gunits%vel = 1.0d+5
  gunits%rho = 6.7699112d-31
  gunits%energy = 1.9890000d+53
  gunits%prs = 6.7699112d-21 
  gunits%time = 3.0856780d+19

end subroutine set_gadget_owls_gimic_units


!> sets the default gadget units from the users guide 
!---------------------------------------------------------
subroutine set_gadget_default_units( gunits )
  type(gadget_units_type) :: gunits
 
  gunits%len = 3.0856780d+21
  gunits%mass = 1.9890000d+43
  gunits%vel = 1.0d+5 
  gunits%time = gunits%len / gunits%vel
  gunits%rho = gunits%mass / gunits%len**3
  gunits%energy = gunits%mass * gunits%len**2 / gunits%time**2
  gunits%prs = gunits%mass / gunits%len / gunits%time**2

end subroutine set_gadget_default_units



!> reads an hdf5 gadget constants
!-----------------------------------------
subroutine read_gadget_constants(snapfile, gconst, simname)
  character(*), intent(in) :: snapfile
  type(gadget_constants_type) :: gconst
  character(*), intent(in), optional :: simname
  integer(i4b) :: fh
  integer(i8b) :: lun


#ifdef hdf5
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
#else
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
#endif

end subroutine read_gadget_constants


!> sets the gadget constants used in the OWLS and GIMIC runs
!-------------------------------------------------------------
subroutine set_gadget_owls_gimic_constants( gconst )
  type(gadget_constants_type) :: gconst

  gconst%pi = 3.1415927d0
  gconst%gamma = 1.6666667d0
  gconst%gravity = 6.6720000d-08
  gconst%solar_mass = 1.9890000d+33
  gconst%solar_lum = 3.8260000d+33 
  gconst%rad_const = 7.5650000d-15 
  gconst%avogadro = 6.0222000d+23
  gconst%boltzmann = 1.3806000d-16
  gconst%gas_const = 83142500.d0
  gconst%c = 2.9979000d+10
  gconst%planck = 6.6262000d-27
  gconst%cm_per_mpc = 3.0856780d+24
  gconst%protonmass = 1.6726000d-24
  gconst%electronmass = 9.1095300d-28 
  gconst%electroncharge = 4.8032000d-10 
  gconst%hubble = 3.2407789d-18
  gconst%t_cmb0 = 2.7280000d0
  gconst%sec_per_megayear = 3.1550000d+13
  gconst%sec_per_year = 31550000.d0 
  gconst%stefan = 7.5657000d-15
  gconst%thompson = 6.6524587d-25
  gconst%ev_to_erg = 1.6021765d-12
  gconst%z_solar = 0.012663729d0

end subroutine set_gadget_owls_gimic_constants




!> writes gadget header data to the screen
!-----------------------------------------
subroutine gadget_header_to_file(ghead,lun)
  type(gadget_header_type), intent(in) :: ghead !< particle header to print
  integer(i8b), intent(in) :: lun

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

#ifdef usempi
  count = 6
  call mpi_bcast( head%npar_file, count, mpi_integer, root, mpi_comm_world, ierr )
  call mpi_bcast( head%mass, count, mpi_double_precision, root, mpi_comm_world, ierr )

  count = 1
  call mpi_bcast( head%a, count, mpi_double_precision, root, mpi_comm_world, ierr )
  call mpi_bcast( head%z, count, mpi_double_precision, root, mpi_comm_world, ierr )
  call mpi_bcast( head%flag_sfr, count, mpi_integer, root, mpi_comm_world, ierr )
  call mpi_bcast( head%flag_feedback, count, mpi_integer, root, mpi_comm_world, ierr )

  count = 6
  call mpi_bcast( head%npar_all, count, mpi_integer, root, mpi_comm_world, ierr )

  count = 1
  call mpi_bcast( head%flag_cooling, count, mpi_integer, root, mpi_comm_world, ierr )
  call mpi_bcast( head%nfiles, count, mpi_integer, root, mpi_comm_world, ierr )
  call mpi_bcast( head%boxlen, count, mpi_double_precision, root, mpi_comm_world, ierr )
  call mpi_bcast( head%OmegaM, count, mpi_double_precision, root, mpi_comm_world, ierr )
  call mpi_bcast( head%OmegaL, count, mpi_double_precision, root, mpi_comm_world, ierr )
  call mpi_bcast( head%h, count, mpi_double_precision, root, mpi_comm_world, ierr )
  call mpi_bcast( head%flag_age, count, mpi_integer, root, mpi_comm_world, ierr )
  call mpi_bcast( head%flag_metals, count, mpi_integer, root, mpi_comm_world, ierr )

  count = 6
  call mpi_bcast( head%npar_hw, count, mpi_integer, root, mpi_comm_world, ierr )

  count = 1
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

#ifdef usempi
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

#ifdef usempi
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
