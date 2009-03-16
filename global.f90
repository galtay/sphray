!> \file global.f90

!> \brief the module that handles global variables
!<
module global_mod
use myf90_mod 
use particle_system_mod, only: particle_system_type
use oct_tree_mod, only: oct_tree_type
use raylist_mod, only: raylist_type
use physical_constants_mod, only: kpc2cm, L_solar, M_solar, Myr2sec
implicit none


type run_planning_type  
   real(r8b)    :: TimeAtSnap        !< t at snapshot
   real(r8b)    :: ScalefacAtSnap    !< scale factor at snapshot
   real(r8b)    :: TimeBetSnaps      !< TimeAtSnap(n+1) - TimeAtSnap(n)
   real(r8b)    :: SnapRunTimes      !< duration to trace each snapshot
   real(r8b)    :: SnapStartTimes    !< t0 for tracing each snapshot 
   integer(i8b) :: RaysFromSrcHeader !< rays listed in source header
   integer(i8b) :: SrcRaysForSnap    !< rays to trace for each snap
end type run_planning_type


! default code units (these are the Gadget system of units)
!-------------------------------------------------------------
real(r8b), parameter :: dflt_cgs_len  = 3.085678d21   !< 1 kpc/h 
real(r8b), parameter :: dflt_cgs_mass = 1.989d43      !< 1e10 Msun/h 
real(r8b), parameter :: dflt_cgs_vel  = 1.0d5         !< 1.0 km/s 
real(r8b), parameter :: dflt_cgs_time = dflt_cgs_len / dflt_cgs_vel                       !< 9.77814e8 years/h
real(r8b), parameter :: dflt_cgs_rho  = dflt_cgs_mass / dflt_cgs_len**3                   !< 1e10 Msun/kpc^3 h^2
real(r8b), parameter :: dflt_cgs_prs  = dflt_cgs_mass / dflt_cgs_len / dflt_cgs_time**2   !< 6.769624e-12 barye h^2
real(r8b), parameter :: dflt_cgs_enrg = dflt_cgs_mass * (dflt_cgs_vel)**2                 !< 1.98892e53 ergs/h
real(r8b), parameter :: dflt_cgs_lum  = dflt_cgs_enrg / dflt_cgs_time                     !< 1.6847e3 Lsun
real(r8b), parameter :: dflt_cgs_mass_solar = dflt_cgs_mass / M_solar                     !< Code Mass Msun/h
real(r8b), parameter :: dflt_photon_per_sec_unit = 1.0d50                                 !< Code photons/s

! global variables
!=====================
type(particle_system_type) :: psys
type(raylist_type) :: globalraylist         !< ray/particle intersections
type(oct_tree_type) :: tree

type(run_planning_type), allocatable :: PLAN(:) !< run plan
real(r8b), allocatable :: OutputTimes(:)        !< times to make outputs at
 
!> global variables type. 
!=========================
type global_variables_type


  ! these are read in directly from the config file
  !-------------------------------------------------

   logical         :: DoTestScenario      !< [Config File] set true if performing a test problem
   character(clen) :: TestScenario        !< [Config File] one of {iliev_test1, iliev_test2, iliev_test3, iliev_test4}

   logical         :: JustInit            !< [Config File] set true to stop after initialization
   logical         :: Comoving            !< [Config File] set true if values to be read are in comoving coords

   real(r8b)       :: IsoTemp             !< [Config File] if non zero all pars fixed @ IsoTemp
   logical         :: FixSnapTemp         !< [Config File] T = fix temp at snapshot values (ignored if IsoTemp /= 0)
                                          !!               F = evolve temperature and hand off between snapshots 

   integer(i8b)    :: IntSeed             !< [Config File] seed for mersenne twister

   real(r8b)       :: StaticFieldSimTime  !< [Config File] sim time for single snapshot jobs
   character(clen) :: StaticSimTimeUnit   !< [Config File] one of {codetime,myr}


   integer(i8b)    :: InputType           !< [Config File] one of {1:Sphray 2:Gadget Public 3:Gadget Tiziana}
   character(clen) :: SnapPath            !< [Config File] dir where particle snapshots are
   character(clen) :: SourcePath          !< [Config File] dir where source snapshots are


   character(clen) :: SpectraFile         !< [Config File] file containing spectra tables
   character(clen) :: b2cdFile            !< [Config File] file containing b2cd tables
   character(clen) :: AtomicRatesFile     !< [Config File] file containnig atomic rates tables

   character(clen) :: ParFileBase         !< [Config File] particle snapshot file base
   character(clen) :: SourceFileBase      !< [Config File] source snapshot file base

   integer(i8b)    :: StartSnapNum        !< [Config File] snapshot to start with
   integer(i8b)    :: EndSnapNum          !< [Config File] snapshot to end with

   integer(i8b)    :: ParFilesPerSnap     !< [Config File] files per particle snapshot
   integer(i8b)    :: SourceFilesPerSnap  !< [Config File] files per source snapshot


   character(clen) :: RayScheme           !< [Config File] one of {raynum, header}
   real(r8b)       :: ForcedRayNumber     !< [Config File] number of rays to trace if RayScheme = raynum

   logical         :: RayStats            !< [Config File] T = massive output file on ray statistics in raystats.dat
   integer(i8b)    :: BndryCond           !< [Config File] one of {-1:reflecting 0:vacuum 1:periodic}

   real(r8b)       :: RecRayTol           !< [Config File] minimum recombination fraction to make a recomb ray
   real(r8b)       :: RayPhotonTol        !< [Config File] fractional ray depletion to stop ray


   logical         :: OnTheSpotH          !< [Config File] T = on the spot approximation for Hydrogen
   logical         :: OnTheSpotHe         !< [Config File] T = on the spot approximation for Helium

   logical         :: HydrogenCaseA       !< [Config File] T = use case A for Hydrogen OTS (ignored if not using OTS)
   logical         :: HeliumCaseA         !< [Config File] T = use case A for Helium OTS (ignored if not using OTS)

   integer(i8b)    :: IonTempSolver       !< [Config File] one of {1:euler, 2:bdf}

   real(r8b)       :: Tfloor              !< [Config File] minimum allowed temperature
   real(r8b)       :: Tceiling            !< [Config File] maximum allowed temperature

   real(r8b)       :: xfloor              !< [Config File] minimum allowed ionization fraction
   real(r8b)       :: xceiling            !< [Config File] maximum allowed ionization fraction

   real(r8b)       :: NeBackGround        !< [Config File] constant background electron number density from metals

   integer(i8b)    :: UpdateAllRays       !< [Config File] update all particles every UpdateAllRays rays
   integer(i8b)    :: RecRaysPerSrcRay    !< [Config File] number of recomb ray for each source ray

   real(r8b)       :: H_mf                !< [Config File] hydrogen mass fraction
   real(r8b)       :: He_mf               !< [Config File] helium mass fraction

   character(clen) :: OutputDir           !< [Config File] path to output directory
   character(clen) :: OutputFileBase      !< [Config File] output file base

   integer(i8b)    :: OutputType          !< [Config File] one of {1:Sphray 2:Gadget}

   character(clen) :: OutputTiming        !< [Config File] one of {standard, forced}
   integer(i8b)    :: NumStdOuts          !< [Config File] if OutputTiming = "standard", # of outputs (plus 1 initial)

   logical         :: DoInitialOutput     !< [Config File] produces output before any raytracing
   integer(i8b)    :: IonFracOutRays      !< [Config File] do mini output every IonFracOutRays src rays

   character(clen) :: ForcedOutFile       !< [Config File] file with forced output times
   character(clen) :: ForcedUnits         !< [Config File] one of {codetime, myr, mwionfrac}

   integer(i8b)    :: PartPerCell         !< [Config File] minimum particles in a tree leaf


   ! these are all set in get_planning_data in main_input.f90
   !-----------------------------------------------------------
   integer(i8b)              :: Nsnaps               !< snap count (particle=source)

   real(r8b) :: cgs_len  !< code length [cm/h]
   real(r8b) :: cgs_mass !< code mass [g/h]
   real(r8b) :: cgs_vel  !< code velocity [cm/s]

   real(r8b) :: cgs_time !< code time [s/h]
   real(r8b) :: cgs_rho  !< code density [g/cm^3 h^2]
   real(r8b) :: cgs_prs  !< code pressure [dyne/cm^2 h^2]
   real(r8b) :: cgs_enrg !< code energy [ergs/h]
   real(r8b) :: cgs_lum  !< code luminosity [ergs/s]
   real(r8b) :: Lunit    !< code source luminosity unit [photons/s]

   real(r8b) :: OmegaM   !< matter / critical density z=0
   real(r8b) :: OmegaB   !< baryon / critical density z=0
   real(r8b) :: OmegaL   !< lambda / critical density z=0
   real(r8b) :: LittleH  !< Hubble parameter z=0 in units of 100 km/s/Mpc

   ! these are set in do_output_planning and do_ray_planning in intialize.f90
   !--------------------------------------------------------------------------
   integer(i8b) :: OutputIndx                     !< keeps track of outputs   
   real(r8b)    :: TotalSimTime                   !< total time to ray trace
   integer(i8b) :: NumTotOuts                     !< total outputs to do
 
   ! these should be reset each time a new snapshot is read in
   !------------------------------------------------------------
   real(r8b) :: dtray_code        !< time between each ray (code units)
   real(r8b) :: dtray_s           !< time between each ray (seconds)

   real(r8b) :: BoxLowers(3)      !< coords of lower x,y,z corner
   real(r8b) :: BoxUppers(3)      !< coords of upper x,y,z corner

   real(r8b) :: BoxLengths(3)     !< side lengths in code units (xf-xi,yf-yi,zf-zi)
   real(r8b) :: BoxLengths_cm(3)  !< box side lengths in cm (xf-xi,yf-yi,zf-zi)
   real(r8b) :: BoxLengths_kpc(3) !< box side lengths in kpc (xf-xi,yf-yi,zf-zi)

   real(r8b) :: BoxVol            !< box volume in code units (xlen*ylen*zlen)
   real(r8b) :: BoxVol_cm         !< box volume in cm^3 (xlen*ylen*zlen)
   real(r8b) :: BoxVol_kpc        !< box volume in kpc^3 (xlen*ylen*zlen)

   real(r8b) :: total_mass      !< summed mass of all particles in a snapshot
   real(r8b) :: total_lum       !< summed luminosity of sources in a snapshot
   real(r8b) :: total_atoms     !< sum of all atoms in computational volume
   real(r8b) :: total_photons   !< sum of all photons to be released
   real(r8b) :: Tcmb_cur        !< CMB temperature for the current snapshot

   
   ! these are updated continuosly while the code runs 
   ! and most are initialized in initialize.f90
   !----------------------------------------------------
   character(clen) :: ionfrac_file !< file where mini outputs are put 
   integer(i8b) :: ionlun          !< lun for ionfrac_file

   character(clen) :: raystat_file !< file where ray stats are put
   integer(i8b) :: raystatlun      !< lun for ray stat file

   integer(i8b) :: CurSnapNum      !< current snapshot number
   integer(i8b) :: rayn            !< current ray number (src + recomb)
   integer(i8b) :: src_rayn        !< current source ray number

   real(r8b) :: time_code          !< current time in code units
   real(r8b) :: time_s             !< current time in seconds

   real(r8b) :: time_elapsed_code  !< elapsed time in code units
   real(r8b) :: time_elapsed_s     !< elapsed time in seconds

   real(r8b) :: mwionfrac          !< mass weighted ionization fraction

   real(r8b) :: TotalSourceRaysCast    !< total rays traced from user def. sources
   real(r8b) :: TotalDiffuseRaysCast   !< total recombination rays traced
   real(r8b) :: IonizingPhotonsPerSec  !< ionizing photons emitted per second
   real(r8b) :: TotalPhotonsCast       !< total number of photons emitted
   real(r8b) :: TotalPhotonsAbsorbed   !< total number of photons absorbed
   real(r8b) :: PhotonsLeavingBox      !< total photons leaving the box
   real(r8b) :: TotalIonizations       !< total number of photoionizations
   real(r8b) :: TotalRecombinations    !< total number of recombinations
   
   real(r8b) :: PeakUpdates            !< max updates for a single particle
   real(r8b) :: AverageUpdatesPerPar   !< average number of updates per particle
   real(r8b) :: ParticleCrossings      !< number of ray / particle intersections
   real(r8b) :: TotalDerivativeCalls   !< times the solver being used has run
   
end type global_variables_type
 

type(global_variables_type) :: GV           !< global variables          



contains


!> reads the config file 
!==============================
subroutine read_config_file(config_file,GV)

  character(clen), intent(in) :: config_file       !< file to read config vars from
  type(global_variables_type), intent(inout) :: GV !< global variables
  character(clen) :: keyword


    write(*,'(A,T29,A)') "using configuration file: ", trim(config_file)


    keyword = "DoTestScenario:"
    call scanfile(config_file,keyword,GV%DoTestScenario)

    keyword = "TestScenario:"
    call scanfile(config_file,keyword,GV%TestScenario)
!-----------------------

    keyword = "JustInit:"
    call scanfile(config_file,keyword,GV%JustInit)

    keyword = "Comoving:"
    call scanfile(config_file,keyword,GV%Comoving)

    keyword = "IsoTemp:"
    call scanfile(config_file,keyword,GV%IsoTemp)

    keyword = "FixSnapTemp:"
    call scanfile(config_file,keyword,GV%FixSnapTemp)

    keyword = "IntSeed:"
    call scanfile(config_file,keyword,GV%IntSeed)

    keyword = "StaticFieldSimTime:"
    call scanfile(config_file,keyword,GV%StaticFieldSimTime)

    keyword = "StaticSimTimeUnit:"
    call scanfile(config_file,keyword,GV%StaticSimTimeUnit)    

    !   input snapshot information
    !------------------------------
    keyword = "InputType:"
    call scanfile(config_file,keyword,GV%InputType)

    keyword = "SnapPath:"
    call scanfile(config_file,keyword,GV%SnapPath)

    keyword = "SourcePath:"
    call scanfile(config_file,keyword,GV%SourcePath)

    keyword = "SpectraFile:"
    call scanfile(config_file,keyword,GV%SpectraFile)

    keyword = "b2cdFile:"
    call scanfile(config_file,keyword,GV%b2cdFile)

    keyword = "AtomicRatesFile:"
    call scanfile(config_file,keyword,GV%AtomicRatesFile)

    keyword = "ParFileBase:"
    call scanfile(config_file,keyword,GV%ParFileBase)

    keyword = "SourceFileBase:"
    call scanfile(config_file,keyword,GV%SourceFileBase)

    keyword = "StartSnapNum:"
    call scanfile(config_file,keyword,GV%StartSnapNum)

    keyword = "EndSnapNum:"
    call scanfile(config_file,keyword,GV%EndSnapNum)

    keyword = "ParFilesPerSnap:"
    call scanfile(config_file,keyword,GV%ParFilesPerSnap)

    keyword = "SourceFilesPerSnap:"
    call scanfile(config_file,keyword,GV%SourceFilesPerSnap)


    !   ray tracing
    !----------------------------
    keyword = "RayScheme:"
    call scanfile(config_file,keyword,GV%RayScheme)

    keyword = "ForcedRayNumber:"
    call scanfile(config_file,keyword,GV%ForcedRayNumber)
 
    keyword = "RayStats:"
    call scanfile(config_file,keyword,GV%RayStats)

    keyword = "BndryCond:"
    call scanfile(config_file,keyword,GV%BndryCond)

    keyword = "RecRayTol:"
    call scanfile(config_file,keyword,GV%RecRayTol)

    keyword = "RayPhotonTol:"
    call scanfile(config_file,keyword,GV%RayPhotonTol)


    !   ion/temp solving
    !----------------------------
    keyword = "OnTheSpotH:"
    call scanfile(config_file,keyword,GV%OnTheSpotH)

    keyword = "OnTheSpotHe:"
    call scanfile(config_file,keyword,GV%OnTheSpotHe)

    keyword = "HydrogenCaseA:"
    call scanfile(config_file,keyword,GV%HydrogenCaseA)

    keyword = "HeliumCaseA:"
    call scanfile(config_file,keyword,GV%HeliumCaseA)

    keyword = "IonTempSolver:"
    call scanfile(config_file,keyword,GV%IonTempSolver)

    keyword = "Tfloor:"
    call scanfile(config_file,keyword,GV%Tfloor)

    keyword = "Tceiling:"
    call scanfile(config_file,keyword,GV%Tceiling)

    keyword = "xfloor:"
    call scanfile(config_file,keyword,GV%xfloor)

    keyword = "xceiling:"
    call scanfile(config_file,keyword,GV%xceiling)

    keyword = "NeBackGround:"
    call scanfile(config_file,keyword,GV%NeBackGround)

    keyword = "UpdateAllRays:"
    call scanfile(config_file,keyword,GV%UpdateAllRays)

    keyword = "RecRaysPerSrcRay:"
    call scanfile(config_file,keyword,GV%RecRaysPerSrcRay)

    keyword = "H_mf:"
    call scanfile(config_file,keyword,GV%H_mf)

    keyword = "He_mf:"
    call scanfile(config_file,keyword,GV%He_mf)

    !   output
    !-------------
    keyword = "OutputDir:"
    call scanfile(config_file,keyword,GV%OutputDir)

    keyword = "OutputFileBase:"
    call scanfile(config_file,keyword,GV%OutputFileBase)

    keyword = "OutputType:"
    call scanfile(config_file,keyword,GV%OutputType)

    keyword = "OutputTiming:"
    call scanfile(config_file,keyword,GV%OutputTiming)

    keyword = "NumStdOuts:"
    call scanfile(config_file,keyword,GV%NumStdOuts)

    keyword = "DoInitialOutput:"
    call scanfile(config_file,keyword,GV%DoInitialOutput)

    keyword = "IonFracOutRays:"
    call scanfile(config_file,keyword,GV%IonFracOutRays)

    keyword = "ForcedOutFile:"
    call scanfile(config_file,keyword,GV%ForcedOutFile)

    keyword = "ForcedUnits:"
    call scanfile(config_file,keyword,GV%ForcedUnits)
 
!--------------------
    keyword = "PartPerCell:"
    call scanfile(config_file,keyword,GV%PartPerCell)

  
    call dummy_check_config_variables(config_file,GV)
    call config_info_to_file(GV)

end subroutine read_config_file


!> run dummy checks on config variables
!========================================
subroutine dummy_check_config_variables(config_file,GV)

  character(clen), intent(in) :: config_file  !< configuration file
  type(global_variables_type), intent(in) :: GV !< global variables
  character(clen) :: Cwarning, Mwarning
  logical :: config_good
  logical :: charmatch

  Cwarning = "please edit " // trim(config_file)
  Mwarning = "please edit Makefile"
  config_good = .true. 

  if (GV%InputType /= 1 .and. GV%InputType /= 2 .and. GV%InputType /= 3) then
     write(*,*) "Input Type ", GV%InputType, " not recognized"
     write(*,*) "must be 1 (Sphray), 2 (Gadget Public), or 3 (Gadget Tiziana)"
     config_good = .false. 
  end if

  if (GV%OutputType /= 1 .and. GV%OutputType /= 2) then
     write(*,*) "Input Type ", GV%OutputType, " not recognized"
     write(*,*) "must be 1 (Sphray) or 2 (Gadget)"
     config_good = .false. 
  end if


  if (GV%Tfloor < 0.0 .or. GV%Tceiling < 0.0) then
     write(*,*) "Tfloor and Tceiling must be greater than or equal to 0.0"
     config_good = .false. 
  end if

  if (GV%Tfloor > 1.0e9 .or. GV%Tceiling > 1.0e9) then
     write(*,*) "Tfloor and Tceiling must be less than or equal to 1.0e9"
     config_good = .false. 
  end if

  if (GV%Tfloor > GV%Tceiling) then
     write(*,*) "Tceiling must be greater than Tfloor"
     config_good = .false. 
  end if

  if (GV%IsoTemp > 0.0) then
  
     if (GV%IsoTemp < GV%Tfloor) then
        write(*,*) "IsoTemp cannot be set lower than Tfloor"
        config_good = .false. 
     end if
     
     if (GV%IsoTemp > GV%Tceiling) then
        write(*,*) "IsoTemp cannot be set higher than Tceiling"
        config_good = .false. 
     end if

  end if

  if (GV%StartSnapNum < 0 .or. GV%EndSnapNum < 0) then
     write(*,*) "Starting and Ending snapshot numbers must be > 0"
     config_good = .false. 
  end if

  if (GV%StartSnapNum > GV%EndSnapNum) then
     write(*,*) "Starting snapshot number cannot be > Ending snapshot number"
     config_good = .false. 
  end if


  if ( GV%IonTempSolver /= 1 .and. GV%IonTempSolver /= 2) then
     config_good = .false.
     write(*,*) "IonTempSolver: ", GV%IonTempSolver, " not recognized"
     write(*,*) "must be '1'=euler or '2'=bdf "
  end if


  charmatch = .false.
  if ( trim(GV%StaticSimTimeUnit) == "codetime" ) charmatch = .true.
  if ( trim(GV%StaticSimTimeUnit) == "myr"      ) charmatch = .true.
  if (.not. charmatch) then
     write(*,*) "StaticSimTimeUnit: ", trim(GV%StaticSimTimeUnit), " not recognized"
     write(*,*) "must be 'codetime' or 'myr' "
     config_good = .false.
  end if

  charmatch = .false.
  if ( trim(GV%RayScheme) == "raynum" ) charmatch = .true.
  if ( trim(GV%RayScheme) == "header"   ) charmatch = .true.
  if (.not. charmatch) then
     write(*,*) "RayScheme: ", trim(GV%RayScheme), " not recognized"
     write(*,*) "must be 'raynum' or 'header' "
     config_good = .false.
  end if

  charmatch = .false.
  if ( trim(GV%OutputTiming) == "standard" ) charmatch = .true.
  if ( trim(GV%OutputTiming) == "forced"   ) charmatch = .true.
  if (.not. charmatch) then
     write(*,*) "OutputTiming: ", trim(GV%OutputTiming), " not recognized"
     write(*,*) "must be 'standard' or 'forced' "
     config_good = .false.
  end if

  charmatch = .false.
  if ( trim(GV%ForcedUnits) == "codetime" ) charmatch = .true.
  if ( trim(GV%ForcedUnits) == "myr"      ) charmatch = .true.
  if ( trim(GV%ForcedUnits) == "mwionfrac") charmatch = .true.
  if (.not. charmatch) then
     write(*,*) "ForcedUnits: ", trim(GV%ForcedUnits), " not recognized"
     write(*,*) "must be 'codetime', 'myr', or 'mwionfrac' "
     config_good = .false.
  end if

#ifdef incHe
  if (GV%He_mf == 0.0) then
     write(*,*) "You have defined the incHe macro in the Makefile, but"
     write(*,*) "the Helium mass fraction is set to 0.0 in the config file."
     write(*,*) "Please comment out incHe in the Makefile or set the Helium"
     write(*,*) "mass fraction to something greater than zero in the config file."
     config_good = .false.
  end if
#else
  if (GV%He_mf /= 0.0) then
     write(*,*) "You have not defined the incHe macro in the Makefile, but"
     write(*,*) "the Helium mass fraction is non-zero in the config file."
     write(*,*) "Please uncomment the incHe line in the Makefile or set the"
     write(*,*) "Helium mass fraction to zero in the config file."
     config_good = .false.
  end if
#endif




  if (GV%DoTestScenario) then  
     charmatch = .false.

     if ( trim(GV%TestScenario) == "iliev_test1" ) then
        charmatch = .true.
        if (GV%IsoTemp /= 1.0d4) then
           config_good = .false.
           write(*,*) "iliev test 1 must have IsoTemp = 1.0e4"
        end if
     end if

     if ( trim(GV%TestScenario) == "iliev_test1He" ) then
        charmatch = .true.
        if (GV%IsoTemp /= 1.0d4) then
           config_good = .false.
           write(*,*) "iliev test 1 He must have IsoTemp = 1.0e4"
        end if
     end if

     if ( trim(GV%TestScenario) == "iliev_test2" ) then
        charmatch = .true.
        if (GV%IsoTemp /= 0.0) then
           config_good = .false.
           write(*,*) "iliev test 2 must have IsoTemp = 0.0"
        end if
     end if

     if ( trim(GV%TestScenario) == "iliev_test3" ) then
        charmatch = .true.
        if (GV%IsoTemp /= 0.0) then
           config_good = .false.
           write(*,*) "iliev test 3 must have IsoTemp = 0.0"
        end if
     end if

     if ( trim(GV%TestScenario) == "iliev_test4" ) then
        charmatch = .true.
        if (GV%IsoTemp /= 0.0) then
           config_good = .false.
           write(*,*) "iliev test 4 must have IsoTemp = 0.0"
        end if
     end if


     if (.not. charmatch) then
        write(*,*) "TestScenario: ", trim(GV%TestScenario), " not recognized"
        write(*,*) "must be 'iliev_test1(He)', 'iliev_test2', 'iliev_test3' or 'iliev_test4' "
        config_good = .false.
     end if
  end if

  if (.not. config_good) then
     write(*,*) Cwarning
     stop
  end if


end subroutine dummy_check_config_variables




!> writes the configuration file information to the output directory
!====================================================================
subroutine config_info_to_file(GV)

  type(global_variables_type), intent(in) :: GV !< global variables

  character(200) :: config_log_file
  integer(i8b) :: lun

  config_log_file = trim(GV%OutputDir) // "/config_values_used.log"
  
  call open_formatted_file_w(config_log_file,lun)

  105 format(T2,A,I3.3)
  111 format(T2,A,I10)
  120 format(T2,A,ES10.3)

  write(lun,*)"============================================"
  write(lun,*)"The SPHRAY configuration file has been read "
  write(lun,*)"The SPHRAY configuration file variables are "
  write(lun,*)"============================================"
  write(lun,*)"Do a test scenario? " , GV%DoTestScenario
  write(lun,*)
  write(lun,*)"Which test scenario? " , trim(GV%TestScenario)
  write(lun,*)
  write(lun,*)"Just initialize? " , GV%JustInit
  write(lun,*)
  write(lun,*)"Input is in comoving coords? " , GV%Comoving
  write(lun,*)
  write(lun,*)"Iso temperature (0.0 = variable temperature): ", GV%IsoTemp
  write(lun,*) 
  write(lun,*)"Fix temperature at snapshot values?: ", GV%FixSnapTemp
  write(lun,*) 
  write(lun,*)"Integer Seed for RNG ", GV%IntSeed
  write(lun,*) 
  write(lun,*)"Static field simulation time: ", GV%StaticFieldSimTime
  write(lun,*) 
  write(lun,*)"Static field time unit: ", GV%StaticSimTimeUnit
  write(lun,*) 

  write(lun,*)

  write(lun,*)"Input Type (1=Sphray, 2=Gadget Public, 3=Gadget Tiziana)", GV%InputType
  write(lun,*)
  write(lun,*)"Path to Snapshot File(s):"
  write(lun,*)trim(GV%SnapPath)
  write(lun,*)
  write(lun,*)"Path to Source File(s):"
  write(lun,*)trim(GV%SourcePath)
  write(lun,*) 
  write(lun,*)"Path to the Spectra file:"
  write(lun,*)trim(GV%SpectraFile)
  write(lun,*) 
  write(lun,*)"Path to the impact parameter -> column depth file:"
  write(lun,*)trim(GV%b2cdFile)
  write(lun,*) 
  write(lun,*)"Path to the atomic rates file:"
  write(lun,*)trim(GV%AtomicRatesFile)

  write(lun,*)
  write(lun,*)"Particle file base: ", trim(GV%ParFileBase)
  write(lun,*)"Source file base:   ", trim(GV%SourceFileBase)
  write(lun,*)
  write(lun,"(A,I3.3)") "Starting snapshot number:", GV%StartSnapNum
  write(lun,"(A,I3.3)") "Ending snapshot number:  ", GV%EndSnapNum
  write(lun,"(A,I3)") "Par Files per snapshot:    ", GV%ParFilesPerSnap
  write(lun,"(A,I10)") "Source Files per snapshot:", GV%SourceFilesPerSnap

  write(lun,*)
  write(lun,*)
  write(lun,*)  "Ray Scheme : ", trim(GV%RayScheme)
  write(lun,120) "Forced Ray Number : ", real(GV%ForcedRayNumber)

  write(lun,*) "Report ray statistics in raystats.dat?", GV%RayStats

  if(GV%BndryCond==-1) then
  write(lun,*)  "Boundry Conditions : ", "reflecting"
  else if(GV%BndryCond==0) then
  write(lun,*)  "Boundry Conditions : ", "vacuum"
  else if(GV%BndryCond==1) then  
  write(lun,*)  "Boundry Conditions : ", "periodic"
  end if

  write(lun,*)  "Rec Ray Frac Tol   : ", GV%RecRayTol
  write(lun,*)  "Ray Photon Tol     : ", GV%RayPhotonTol

  write(lun,*)  "Use On The Spot for Hydrogen? : ", GV%OnTheSpotH
  write(lun,*)  "Use On The Spot for Helium?   : ", GV%OnTheSpotHe

  write(lun,*)  "Use Case A recombination rates for Hydrogen? :", GV%HydrogenCaseA
  write(lun,*)  "Use Case A recombination rates for Helium?   :", GV%HeliumCaseA

  write(lun,*)  "Ionization and temperature solver :", GV%IonTempSolver

  write(lun,*)  "Temperature floor   : ", GV%Tfloor
  write(lun,*)  "Temperature ceiling : ", GV%Tceiling
  write(lun,*)  "Ionization floor    : ", GV%xfloor
  write(lun,*)  "Ionization ceiling  : ", GV%xceiling


  write(lun,*)  "ne background      : ", GV%NeBackGround
  write(lun,*)  "Rays between all particle update:  ", GV%UpdateAllRays
  write(lun,*)  "Recomb ray / src ray: ", GV%RecRaysPerSrcRay

  write(lun,*) 
  write(lun,*) 
  write(lun,*)  "Hydrogen Mass Fraction: ", GV%H_mf
  write(lun,*)  "Helium Mass Fraction: ", GV%He_mf
  
  write(lun,*)
  write(lun,*)
  write(lun,*)  "Output Dir         : ", trim(GV%OutputDir)
  write(lun,*)  "Output File Base   : ", trim(GV%OutputFileBase)
  write(lun,*)  "Output Type (1=Sphray, 2=Gadget) : ", GV%OutputType

  write(lun,*)  "Output timing plan : ", trim(GV%OutputTiming)

  write(lun,*)  "Number Std Outs    : ", GV%NumStdOuts
  write(lun,*)  "Do Initial Output  : ", GV%DoInitialOutput

  write(lun,*)  "Ion Frac Out Rays  : ", GV%IonFracOutRays

  write(lun,*)  "ForcedOutFile      : ", trim(GV%ForcedOutFile)
  write(lun,*)  "ForcedUnits        : ", trim(GV%ForcedUnits)

  write(lun,*)  "Particles Per Tree Cell: ", GV%PartPerCell


  write(lun,*)"====================================="
  write(lun,*)"   End SPHRAY Configuration Output   "
  write(lun,*)"====================================="
  write(lun,*)
  write(lun,*)

  if (GV%FixSnapTemp) then
     write(lun,*) "***********************************************************"
     write(lun,*) "you are running a constant temperature simulation."
     write(lun,*) "the temperatures are fixed at the readin snapshot values"
     write(lun,*) "***********************************************************"
  else
     if (GV%IsoTemp /= 0.0) then
        write(lun,*) "***********************************************************"
        write(lun,*) "you are running a constant temperature simulation."
        write(lun,*) "the temperature is fixed at T (K) = ", GV%IsoTemp
        write(lun,*) "***********************************************************"
     end if
  end if

 
  close(lun)

end subroutine config_info_to_file






end module global_mod
