!> \file global.F90

!> \brief the module that handles global variables
!<
module global_mod
use myf90_mod 
use gadget_header_class, only: gadget_header_type
use gadget_header_class, only: gadget_units_type
use gadget_header_class, only: gadget_constants_type
use particle_system_mod, only: particle_system_type
use oct_tree_mod, only: oct_tree_type
use raylist_mod, only: raylist_type
use atomic_rates_mod, only: atomic_rates_table_type
use atomic_rates_mod, only: atomic_rates_type
implicit none


!> snapshot information type. 
!==============================
type snap_info_type
   real(r8b)    :: TimeAt            !< t at snapshot
   real(r8b)    :: ScalefacAt        !< scale factor at snapshot
   real(r8b)    :: TimeToNext        !< TimeAt(n+1) - TimeAt(n)
   real(r8b)    :: RunTime           !< duration to trace snapshot
   real(r8b)    :: StartTime         !< t0 for tracing snapshot 
   integer(i8b) :: RaysFromSrcHeader !< rays listed in source header
   integer(i8b) :: SrcRays           !< rays to trace for snapshot
end type snap_info_type


!> run planning type.
!======================
type run_planning_type  
   type(snap_info_type), allocatable :: snap(:) !< snap information
   real(r8b), allocatable :: OutputTimes(:)     !< times to make outputs 
end type run_planning_type


! global variables
!=====================
type(particle_system_type) :: psys        !< particles + sources + box
type(raylist_type) :: globalraylist       !< ray/particle intersections
type(oct_tree_type) :: tree               !< octree

type(atomic_rates_table_type) :: rtable   !< rates read in from file
type(atomic_rates_type) :: isoT_k         !< static rates for iso-temperature run
type(atomic_rates_type) :: cmbT_k         !< static rates for cmb-temperature
type(atomic_rates_type) :: xHII_k         !< static rates for xHII-temperature 

type(run_planning_type) :: PLAN           !< run plan

type(gadget_header_type), allocatable :: saved_gheads(:,:) !< all headers (nsnaps,nfiles)
type(gadget_constants_type) :: gconst                      !< gadget constants

 
!> global variables type. 
!=========================
type global_variables_type


  ! these are read in directly from the config file
  !-------------------------------------------------

   integer(i4b)    :: Verbosity           !< [Config File] 0=silent, 1=whisper, 2=talk, 3=debug

   logical         :: DoTestScenario      !< [Config File] set true if performing a test problem
   character(clen) :: TestScenario        !< [Config File] one of {iliev_test1, iliev_test2, iliev_test3, iliev_test4}

   logical         :: JustInit            !< [Config File] set true to stop after initialization
   logical         :: Comoving            !< [Config File] set true if values to be read are in comoving coords

   real(r8b)       :: IsoTemp             !< [Config File] if non zero all pars fixed @ IsoTemp
   logical         :: FixSnapTemp         !< [Config File] T = fix temp at snapshot values (ignored if IsoTemp /= 0)

   real(r8b)       :: EOStemp             !< [Config File] if non-negative, initialize EOS particles w/ T = EOStemp 
   real(r8b)       :: InitxHI             !< [Config File] if non-negative, all xHI initialized to this value
   logical         :: RayDepletion        !< [Config File] remove photons from rays as they travel?

   integer(i8b)    :: IntSeed             !< [Config File] seed for mersenne twister

   real(r8b)       :: StaticFieldSimTime  !< [Config File] sim time for single snapshot jobs
   character(clen) :: StaticSimTimeUnit   !< [Config File] one of {codetime,myr}


   integer(i4b)    :: InputType           !< [Config File] one of Gadget {1: Public 2: Cooling 3: HDF5 4: Bromm}
   character(clen) :: SnapPath            !< [Config File] dir where particle snapshots are
   character(clen) :: SourcePath          !< [Config File] dir where source snapshots are


   character(clen) :: SpectraFile         !< [Config File] file containing spectra tables
   character(clen) :: b2cdFile            !< [Config File] file containing b2cd tables
   character(clen) :: AtomicRatesFile     !< [Config File] file containnig atomic rates tables

   character(clen) :: ParFileBase         !< [Config File] particle snapshot file base
   character(clen) :: SourceFileBase      !< [Config File] source snapshot file base

   integer(i4b)    :: StartSnapNum        !< [Config File] snapshot to start with
   integer(i4b)    :: EndSnapNum          !< [Config File] snapshot to end with

   integer(i4b)    :: ParFilesPerSnap     !< [Config File] files per particle snapshot
   integer(i4b)    :: SourceFilesPerSnap  !< [Config File] files per source snapshot


   character(clen) :: RayScheme           !< [Config File] one of {raynum, header}
   real(r8b)       :: ForcedRayNumber     !< [Config File] number of rays to trace if RayScheme = raynum

   logical         :: RayStats            !< [Config File] T = massive output file on ray statistics in raystats.dat
   integer(i4b)    :: BndryCond           !< [Config File] one of {-1:reflecting 0:vacuum 1:periodic}

   real(r8b)       :: RecRayTol           !< [Config File] minimum recombination fraction to make a recomb ray
   real(r8b)       :: RayPhotonTol        !< [Config File] fractional ray depletion to stop ray


   logical         :: OnTheSpotH          !< [Config File] T = on the spot approximation for Hydrogen
   logical         :: OnTheSpotHe         !< [Config File] T = on the spot approximation for Helium

   logical         :: HydrogenCaseA       !< [Config File] T = use case A for Hydrogen OTS (ignored if not using OTS)
   logical         :: HeliumCaseA         !< [Config File] T = use case A for Helium OTS (ignored if not using OTS)

   integer(i4b)    :: IonTempSolver       !< [Config File] one of {1:euler, 2:bdf}

   real(r8b)       :: Tfloor              !< [Config File] minimum allowed temperature
   real(r8b)       :: Tceiling            !< [Config File] maximum allowed temperature

   real(r8b)       :: xfloor              !< [Config File] minimum allowed ionization fraction
   real(r8b)       :: xceiling            !< [Config File] maximum allowed ionization fraction

   real(r8b)       :: NeBackGround        !< [Config File] constant background electron number density from metals

   integer(i8b)    :: NraysUpdateNoHits   !< [Config File] update all pars not hit by a ray in last NraysUpdateNoHits
   integer(i8b)    :: RecRaysPerSrcRay    !< [Config File] number of recomb rays for each source ray

   real(r8b)       :: H_mf                !< [Config File] hydrogen mass fraction
   real(r8b)       :: He_mf               !< [Config File] helium mass fraction

   character(clen) :: OutputDir           !< [Config File] path to output directory
   character(clen) :: OutputFileBase      !< [Config File] output file base

   integer(i4b)    :: OutputType          !< [Config File] one of {1:Standard Binary Gadget 2:HDF5 Gadget}

   character(clen) :: OutputTiming        !< [Config File] one of {standard, forced}
   integer(i4b)    :: NumStdOuts          !< [Config File] if OutputTiming = "standard", # of outputs (maybe +1 initial)

   logical         :: DoInitialOutput     !< [Config File] produces output before any raytracing
   integer(i8b)    :: IonFracOutRays      !< [Config File] do mini output every IonFracOutRays src rays

   character(clen) :: ForcedOutFile       !< [Config File] file with forced output times
   character(clen) :: ForcedUnits         !< [Config File] one of {codetime, myr, mwionfrac, vwionfrac}

   integer(i4b)    :: PartPerCell         !< [Config File] minimum particles in a tree leaf

   character(clen) :: config_file         !< name of the config file


   ! particle and source sizes
   !-----------------------------------------------------------
   integer(i4b) :: bytesperpar  !< bytes per particle
   integer(i4b) :: bytespersrc  !< bytes per source
   real(r8b) :: MB              !< tracks memory consumption


   ! these are all set in get_planning_data in main_input.f90
   !-----------------------------------------------------------
   integer(i4b) :: Nsnaps !< snap count (particle=source)

   real(r8b) :: cgs_len   !< code length [cm/h]
   real(r8b) :: cgs_mass  !< code mass [g/h]
   real(r8b) :: cgs_vel   !< code velocity [cm/s] 

   real(r8b) :: cgs_time  !< code time [s/h]
   real(r8b) :: cgs_rho   !< code density [g/cm^3 h^2]
   real(r8b) :: cgs_prs   !< code pressure [dyne/cm^2 h^2]
   real(r8b) :: cgs_enrg  !< code energy [ergs/h]
   real(r8b) :: Lunit     !< code source luminosity unit [photons/s]

   real(r8b) :: LenFac_cm  !< code (input) length -> cm = cgs_len * a / h
   real(r8b) :: MassFac_g  !< code (input) mass -> g = cgs_mass / h
   real(r8b) :: TimeFac_s  !< code (input) time -> s = cgs_time / h

   real(r8b) :: OmegaM    !< matter / critical density z=0
   real(r8b) :: OmegaB    !< baryon / critical density z=0
   real(r8b) :: OmegaL    !< lambda / critical density z=0
   real(r8b) :: LittleH   !< Hubble parameter z=0 in units of 100 km/s/Mpc

   ! these are set in do_output_planning and do_ray_planning in intialize.f90
   !--------------------------------------------------------------------------
   integer(i4b) :: OutputIndx                     !< keeps track of outputs   
   real(r8b)    :: TotalSimTime                   !< total time to ray trace
   integer(i4b) :: NumTotOuts                     !< total outputs to do
 
   ! these should be reset each time a new snapshot is read in
   !------------------------------------------------------------
   real(r8b) :: dtray_code        !< time between each ray (code units)
   real(r8b) :: dtray_s           !< time between each ray (seconds)

   real(r8b) :: BoxLwrsComoh(3)   !< comoving h-mangled input coords of lower x,y,z corner
   real(r8b) :: BoxLwrsComo(3)    !< comoving non h-mangled input coords of lower x,y,z corner
   real(r8b) :: BoxLwrsPhysh(3)   !< physical h-mangled input coords of lower x,y,z corner
   real(r8b) :: BoxLwrsPhys(3)    !< physical non h-mangled input coords of lower x,y,z corner

   real(r8b) :: BoxUprsComoh(3)   !< comoving h-mangled input coords of upper x,y,z corner
   real(r8b) :: BoxUprsComo(3)    !< comoving non h-mangled input coords of upper x,y,z corner
   real(r8b) :: BoxUprsPhysh(3)   !< physical h-mangled input coords of upper x,y,z corner
   real(r8b) :: BoxUprsPhys(3)    !< physical non h-mangled input coords of upper x,y,z corner

   real(r8b) :: BoxLensComoh(3)   !< comoving h-mangled side lengths (xf-xi,yf-yi,zf-zi)
   real(r8b) :: BoxLensComo(3)    !< comoving non h-mangled side lengths (xf-xi,yf-yi,zf-zi)
   real(r8b) :: BoxLensPhysh(3)   !< physical h-mangled side lengths (xf-xi,yf-yi,zf-zi)
   real(r8b) :: BoxLensPhys(3)    !< physical non h-mangled side lengths (xf-xi,yf-yi,zf-zi)

   real(r8b) :: BoxLensComoh_cm(3)  !< comoving h-mangled side lengths [cm/h] (xf-xi,yf-yi,zf-zi)
   real(r8b) :: BoxLensComo_cm(3)   !< comoving non h-mangled side lengths [cm] (xf-xi,yf-yi,zf-zi)
   real(r8b) :: BoxLensPhysh_cm(3)  !< physical h-mangled side lengths [cm/h] (xf-xi,yf-yi,zf-zi)
   real(r8b) :: BoxLensPhys_cm(3)   !< physical non h-mangled side lengths [cm] (xf-xi,yf-yi,zf-zi)

   real(r8b) :: BoxLensComoh_kpc(3)  !< comoving h-mangled side lengths [kpc/h] (xf-xi,yf-yi,zf-zi)
   real(r8b) :: BoxLensComo_kpc(3)   !< comoving non h-mangled side lengths [kpc] (xf-xi,yf-yi,zf-zi)
   real(r8b) :: BoxLensPhysh_kpc(3)  !< physical h-mangled side lengths [kpc/h] (xf-xi,yf-yi,zf-zi)
   real(r8b) :: BoxLensPhys_kpc(3)   !< physical non h-mangled side lengths [kpc] (xf-xi,yf-yi,zf-zi)

   real(r8b) :: BoxVolComoh       !< comoving h-mangled box volume (xlen*ylen*zlen)
   real(r8b) :: BoxVolComo        !< comoving non h-mangled box volume (xlen*ylen*zlen)
   real(r8b) :: BoxVolPhysh       !< physical h-mangled box volume (xlen*ylen*zlen)
   real(r8b) :: BoxVolPhys        !< physical non h-mangled box volume (xlen*ylen*zlen)

   real(r8b) :: BoxVolComoh_cm       !< comoving h-mangled box volume [cm^3/h^3] (xlen*ylen*zlen)
   real(r8b) :: BoxVolComo_cm        !< comoving non h-mangled box volume [cm^3] (xlen*ylen*zlen)
   real(r8b) :: BoxVolPhysh_cm       !< physical h-mangled box volume [cm^3/h^3] (xlen*ylen*zlen)
   real(r8b) :: BoxVolPhys_cm        !< physical non h-mangled box volume [cm^3] (xlen*ylen*zlen)

   real(r8b) :: BoxVolComoh_kpc       !< comoving h-mangled box volume [kpc^3/h^3] (xlen*ylen*zlen)
   real(r8b) :: BoxVolComo_kpc        !< comoving non h-mangled box volume [kpc^3] (xlen*ylen*zlen)
   real(r8b) :: BoxVolPhysh_kpc       !< physical h-mangled box volume [kpc^3/h^3] (xlen*ylen*zlen)
   real(r8b) :: BoxVolPhys_kpc        !< physical non h-mangled box volume [kpc^3] (xlen*ylen*zlen)

   real(r8b) :: total_mass      !< summed mass of all particles in a snapshot
   real(r8b) :: total_lum       !< summed luminosity of sources in a snapshot
   real(r8b) :: total_atoms     !< sum of all atoms in computational volume
   real(r8b) :: total_photons   !< sum of all photons to be released
   real(r8b) :: Tcmb_cur        !< CMB temperature for the current snapshot

   real(r8b) :: sf_gamma_eos    !< index for polytropic equation of state for star forming gas.

   real(r8b) :: UVB_gammaHI_cloudy !< magnitude of UVB from cloudy ionization table
   
   ! these are updated continuosly while the code runs 
   ! and most are initialized in initialize.f90
   !----------------------------------------------------
   character(clen) :: ionfrac_file !< file where mini outputs are put 
   integer(i4b) :: ionlun          !< lun for ionfrac log file

   character(clen) :: raystat_file !< file where ray stats are put
   integer(i4b) :: raystatlun      !< lun for ray stat log file

   character(clen) :: pardata_file !< file with particle data summaries
   integer(i4b) :: pardatalun      !< lun for particle data log file
   
   character(clen) :: srcdata_file !< file with source data summaries
   integer(i4b) :: srcdatalun      !< lun for source data log file

   character(clen) :: rayplan_file  !< file for ray planning
   integer(i4b) :: rayplanlun      !< lun for ray planning file

   character(clen) :: outplan_file  !< file for out planning
   integer(i4b) :: outplanlun      !< lun for out planning file

   integer(i4b) :: CurSnapNum      !< current snapshot number
   integer(i8b) :: rayn            !< current ray number (src + recomb)
   integer(i8b) :: src_rayn        !< current source ray number

   real(r8b) :: time_code          !< current time in code units
   real(r8b) :: time_s             !< current time in seconds

   real(r8b) :: time_elapsed_code  !< elapsed time in code units
   real(r8b) :: time_elapsed_s     !< elapsed time in seconds

   real(r8b) :: nwionfrac          !< number weighted ionization fraction
   real(r8b) :: mwionfrac          !< mass weighted ionization fraction
   real(r8b) :: vwionfrac          !< volume weighted ionization fraction

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










end module global_mod
