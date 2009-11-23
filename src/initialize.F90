!> \file initialize.F90 

!> \brief Initialization module
!! 
!! This Module Contains the subroutines to initialize all of the non  
!! particle system variables.     
!<
module initialize_mod
use myf90_mod
use physical_constants_mod
use config_mod, only: read_config_file
use mt19937_mod, only: init_mersenne_twister
use b2cd_mod, only: read_b2cd_file
use spectra_mod, only: read_spectra_file
use atomic_rates_mod, only: read_atomic_rates_file
use atomic_rates_mod, only: set_xHII_atomic_rates
use atomic_rates_mod, only: set_iso_atomic_rates
use iliev_comparison_project_mod, only: initialize_iliev_tests
use main_input_mod, only: get_planning_data
use global_mod, only: GV, PLAN
use ray_mod, only: curface
implicit none

  contains

!> Main initialization routine
!---------------------------------
subroutine initialize(config_file)

  character(clen), intent(in) :: config_file !< configuration file
  
  ! these subroutines are in other files 
  ! indicated in the use statements above

  call read_config_file(config_file)

  call init_mersenne_twister(GV%IntSeed)
  call read_b2cd_file(GV%b2cdFile)
  call read_spectra_file(GV%SpectraFile)

  call read_atomic_rates_file(GV%AtomicRatesFile, GV%OutputDir)
  call set_xHII_atomic_rates(1.0d4)
  if (GV%IsoTemp /= 0.0) then
     call set_iso_atomic_rates(GV%IsoTemp)
  end if

  call initialize_iliev_tests()

  call get_planning_data()

  call do_output_planning()

  call do_ray_planning()

  call initialize_global_variables()


          
end subroutine initialize

  
!-----------------------------------------------------------------------------
!> Plans output times.  For a single snapshot (static field) this 
!! amounts to setting the start time to the value in the particle header 
!! and the run time to the value in the config file.  
!! For multiple snapshots this amounts to using the times (or scale factors
!! for cosmological runs) in the headers to plan the raytracing such that
!! the time represented by the snapshot will be in the middle of the raytracing
!! block.  Note that this means that the time can start negative if the 
!! initial snapshot has t=0.  

subroutine do_output_planning()

  character(clen), parameter :: myname="do_output_planning"
  logical, parameter :: crash=.true.
  integer, parameter :: verb=1

  integer(i8b) :: i, lun
  integer(i8b) :: Ni, Nf
  logical :: fthere

  integer(i8b) :: loglun
  character(clen) :: logfile
  character(clen) :: str


  call mywrite(" doing output planning",verb) 
  
  Ni = GV%StartSnapNum
  Nf = GV%EndSnapNum

  GV%Nsnaps = Nf - Ni + 1
  
  ! At this point TimeAt, ScalefacAt, and SrcRays should have been readin.
  ! It is important that the times in the particle headers be correct.  
  ! * single snapshot runs  
  !   ** start time - particle header
  !   ** total time - config file 
  ! * multiple snapshot runs 
  !   ** cosmological (Comoving=true in config file)
  !       scale factors are used to calculate start time and sim time
  !   ** newtonian (Comoving=false in config file)
  !       time in particle headers are used to calculate start time and 
  !       sim time
  

  ! convert static field sim time to code time if need be
  !-------------------------------------------------------
  if ( trim(GV%StaticSimTimeUnit) == "myr" ) then
     GV%StaticFieldSimTime = GV%StaticFieldSimTime * Myr2sec     ! [s]
     GV%StaticFieldSimTime = GV%StaticFieldSimTime * GV%LittleH  ! [s/h]
     GV%StaticFieldSimTime = GV%StaticFieldSimTime / GV%cgs_time ! [code]
  end if
  

  ! for a single snapshot
  !-----------------------
  if (GV%Nsnaps==1) then          

     PLAN%snap(Ni)%RunTime   = GV%StaticFieldSimTime    ! from config
     PLAN%snap(Ni)%StartTime = PLAN%snap(Ni)%TimeAt ! from header

     GV%TotalSimTime = GV%StaticFieldSimTime ! from config

  ! for multiple snapshots
  !------------------------
  else if (GV%Nsnaps/=1) then    


     ! calculate code time between snapshots
     !----------------------------------------
     do i = Ni, Nf-1
        PLAN%snap(i)%TimeToNext = PLAN%snap(i+1)%TimeAt - PLAN%snap(i)%TimeAt
     end do

     ! figure out the amount of time each snapshot should be raytraced.
     !-----------------------------------------------------------------
     PLAN%snap(Ni)%RunTime = PLAN%snap(Ni)%TimeToNext   ! for 1st
     PLAN%snap(Nf)%RunTime = PLAN%snap(Nf-1)%TimeToNext ! and last
     do i = Ni+1, Nf-1
        PLAN%snap(i)%RunTime = 0.5 *(PLAN%snap(i-1)%TimeToNext + PLAN%snap(i)%TimeToNext)
     end do
     GV%TotalSimTime = sum(PLAN%snap(:)%RunTime)
     
     ! figure out starting times for each snap (different than TimeAt) 
     ! snap(i)%StartTime + 0.5 * snap(i)%RunTime = snap(i)%TimeAt
     !---------------------------------------------------------------------
     PLAN%snap(Ni)%StartTime = PLAN%snap(Ni)%TimeAt - (0.5 * PLAN%snap(Ni)%RunTime) 
     do i = Ni+1, Nf
        PLAN%snap(i)%StartTime = PLAN%snap(i-1)%StartTime + PLAN%snap(i-1)%RunTime
     end do
     
  end if


  ! output planning ! 
  !-----------------!
  GV%OutputIndx = 1 
  
  ! == for standard output == 
  if (trim(GV%OutputTiming)=="standard") then
     GV%NumTotOuts = GV%NumStdOuts
     allocate( PLAN%OutputTimes(1:GV%NumTotOuts) )
     do i = 1,GV%NumTotOuts
        PLAN%OutputTimes(i) = PLAN%snap(Ni)%StartTime + i * GV%TotalSimTime / GV%NumTotOuts
     end do
     PLAN%OutputTimes(GV%NumTotOuts) = PLAN%snap(Ni)%StartTime + GV%TotalSimTime

     
  ! == for forced output == 
  else if (trim(GV%OutputTiming)=="forced") then
     inquire(file=GV%ForcedOutFile,exist=fthere)
     if (.not. fthere) then
        str = "  cannot find forced output times file: " // trim(GV%ForcedOutFile)
        call myerr(str,myname,crash)
     end if
     call open_formatted_file_r(GV%ForcedOutFile,lun)
     read(lun,*) GV%NumTotOuts
     allocate( PLAN%OutputTimes(GV%NumTotOuts) )
     do i = 1,GV%NumTotOuts
        read(lun,*) PLAN%OutputTimes(i)
     end do
     close(lun)

     ! convert forced output times to code units if need be
     !------------------------------------------------------
     if ( trim(GV%ForcedUnits) == "myr" ) then
        PLAN%OutputTimes = PLAN%OutputTimes * Myr2sec     ! [s]
        PLAN%OutputTimes = PLAN%OutputTimes * GV%LittleH  ! [s/h]
        PLAN%OutputTimes = PLAN%OutputTimes / GV%cgs_time ! [code]
     end if

  end if   


  logfile = trim(GV%OutputDir) // "/" // "output_planning.log"
  call open_formatted_file_w(logfile,loglun)
  write(loglun,'(A,A,A)') "planning output times (OutputTiming = ", trim(GV%OutputTiming), ")"
  
  write(loglun,'(A,I3,A)') "SPHRAY will generate", GV%NumTotOuts, " full outputs"
  if (GV%DoInitialOutput) write(loglun,'(A)') "plus one initial full output"
  write(loglun,*) 

  write(loglun,*) "start time [code]: ", PLAN%snap(Ni)%StartTime
  write(loglun,*) "start time [Myr]:  ", PLAN%snap(Ni)%StartTime * GV%cgs_time / GV%LittleH / Myr2sec
  write(loglun,*)

  105 format (T3,A,I3.3,A,T20,ES12.5,T40,ES12.5)    
  110 format(T20,A,T40,A)
  write(loglun,110) "(code units)", "(Myrs)"
  do i = GV%OutputIndx,GV%NumTotOuts
     write(loglun,105) "output ", i, " : ", PLAN%OutputTimes(i), &
                                            PLAN%OutputTimes(i) * GV%cgs_time / GV%LittleH / Myr2sec
  end do

  write(loglun,*)
  write(loglun,*) "total simulation time [code] = ", GV%TotalSimTime
  write(loglun,*) "total simulation time [Myr]  = ", GV%TotalSimTime * GV%cgs_time / GV%LittleH / Myr2sec

  close(loglun)
    
end subroutine do_output_planning

!-------------------------------------------------------------
!> determine the number of rays to be traced in each snapshot. 
subroutine do_ray_planning()

  integer(i8b) :: i
  real(r8b) :: code2Myr
  integer(i8b) :: loglun
  character(200) :: logfile

  write(*,'(A)') "doing raytracing planning"

  code2Myr = GV%cgs_time / GV%LittleH / Myr2sec

  logfile = trim(GV%OutputDir) // "/" // "raytracing_planning.log"
  call open_formatted_file_w(logfile,loglun)


  if (trim(GV%RayScheme)=="header") then
     ! this works for multiple snapshots.  the number of rays to trace
     ! has already been read from the source files
     write(loglun,*) "RayScheme in cofig file = ", trim(GV%RayScheme)
     write(loglun,*) "Using ray numbers in source file header(s)"
     PLAN%snap(:)%SrcRays = PLAN%snap(:)%RaysFromSrcHeader

  else if(trim(GV%RayScheme)=="raynum") then
     if (GV%StartSnapNum==GV%EndSnapNum) then ! static field simulation
        PLAN%snap(GV%StartSnapNum)%SrcRays = GV%ForcedRayNumber
     else if (GV%StartSnapNum/=GV%EndSnapNum) then ! multiple snapshots
        write(*,*) "ForcedRayNumber ray planning for multiple snapshot runs "
        write(*,*) "not supported yet.  Please specify the number of rays "
        write(*,*) "in the source file headers and set RayScheme = 'header' "
        write(*,*) "in the config file.  Ideally the number of rays should"
        write(*,*) "be proportional to the total luminosity of the snapshot"
        stop
        ! use physical snapshot times and total luminosities to distribute 
        ! the rays to the snapshots so that each ray will have roughly the
        ! same energy
     end if

  else
     write(*,*) "keyword RayScheme not specefied correctly"
     write(*,*) "RayScheme = ", trim(GV%RayScheme)
     write(*,*) "please edit the config file. "
     stop
  end if


  99 format(A,I3,4A,L5,A)
  write(loglun,99) "ray planning: snapshots = ", GV%Nsnaps, &
                   " (RayScheme = ", trim(GV%RayScheme),")", &
                   " (Comoving =", GV%Comoving,")"

  write(loglun,*) 
  write(loglun,*) 



  100 format(A,  T8,A,      T22,A,      T36,A,      T50,A,      T64,A)
  101 format(I3, T8,ES12.6, T22,ES12.6, T36,ES12.6, T50,ES12.6, T64,ES12.6)

  write(loglun,*) "ray plan in code units"
  write(loglun,100) "snap", "t @ snap", "t start", "t end", "dt", "Src rays"
  do i = GV%StartSnapNum,GV%EndSnapNum
     write(loglun,101) i, & 
                       PLAN%snap(i)%TimeAt, &
                       PLAN%snap(i)%StartTime, &
                       PLAN%snap(i)%StartTime + PLAN%snap(i)%RunTime, &
                       PLAN%snap(i)%RunTime, &
                       real(PLAN%snap(i)%SrcRays)
  end do


  write(loglun,*) 
  write(loglun,*) 

  write(loglun,*) "ray plan in Myrs"
  write(loglun,100) "snap", "t @ snap", "t start", "t end", "dt", "Src rays"
  do i = GV%StartSnapNum,GV%EndSnapNum
     write(loglun,101) i, & 
                       PLAN%snap(i)%TimeAt * code2Myr, &
                       PLAN%snap(i)%StartTime * code2Myr, &
                       (PLAN%snap(i)%StartTime + PLAN%snap(i)%RunTime) * code2Myr, &
                       PLAN%snap(i)%RunTime * code2Myr, &
                       real(PLAN%snap(i)%SrcRays)
  end do


  write(loglun,*)
  write(loglun,*) "total simulation time [code] = ", GV%TotalSimTime
  write(loglun,*) "total simulation time [Myr]  = ", GV%TotalSimTime * code2Myr
  close(loglun)


end subroutine do_ray_planning


!=======================================================================
!> when a new particle snapshot is read in, this routine should be called
!! to initialize some global variables.
subroutine initialize_global_variables()
  use particle_system_mod, only: calc_bytes_per_particle
  use ray_mod, only: raystatbuffsize

  logical :: verbose 
  integer(i8b) :: StartSnapNum, EndSnapNum
  integer(i8b) :: bpp !< bytes per particle
  integer(i8b) :: bps !< bytes per source

  write(*,'(A)') "initializing global variables"

  if (GV%IonTempSolver == 1) then
     write(*,'(A)') "using Implicit Euler ionization solver"
  else if (GV%IonTempSolver == 2) then
     write(*,'(A)') "using Backwards Difference ionization solver"
  end if

  verbose = .true.

  StartSnapNum = GV%StartSnapNum
  EndSnapNum = GV%EndSnapNum

  GV%ionfrac_file = trim(GV%OutputDir) // "/ionfrac.log"
  call open_formatted_file_w(GV%ionfrac_file,GV%ionlun)

  GV%raystat_file = trim(GV%OutputDir) // "/raystats.dat"
  if (GV%raystats) then
     call open_unformatted_file_w(GV%raystat_file, GV%raystatlun)
     write(GV%raystatlun) sum(PLAN%snap(:)%SrcRays), raystatbuffsize 
  end if

                 
  GV%CurSnapNum = GV%StartSnapNum
  GV%rayn = 0
  GV%src_rayn = 0
 
  GV%time_code = PLAN%snap(GV%StartSnapNum)%StartTime
  GV%time_s = GV%time_code * GV%cgs_time / GV%LittleH
      
  GV%TotalSourceRaysCast = 0.0
  GV%TotalDiffuseRaysCast = 0.0
  GV%IonizingPhotonsPerSec = 0.0
  GV%TotalPhotonsCast = 0.0
  GV%TotalPhotonsAbsorbed = 0.0
  GV%PhotonsLeavingBox = 0.0
  GV%TotalIonizations = 0.0
  GV%TotalRecombinations = 0.0
  
  GV%PeakUpdates = 0.0
  GV%AverageUpdatesPerPar = 0.0
  GV%ParticleCrossings = 0.0
  GV%TotalDerivativeCalls = 0.0
  
  ! set background radiation temperature
  GV%Tcmb_cur = CMBtempNow / PLAN%snap(StartSnapNum)%ScalefacAt  


  call calc_bytes_per_particle(bpp, bps)
  GV%bytesperpar = bpp
  GV%bytespersrc = bps

  ! only important for runs with a background source (emisprf=-3)
  curface = 0


end subroutine initialize_global_variables


end module initialize_mod
