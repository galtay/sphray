!> \file gadget_input.F90

!> \brief Handles readin of GADGET formatted files
!<

module gadget_input_mod
use myf90_mod
use gadget_header_class
use gadget_input_hdf5_mod
use particle_system_mod, only: particle_system_type
use particle_system_mod, only: set_collisional_ionization_equilibrium, set_ye
use global_mod, only: psys, PLAN, GV
use global_mod, only: saved_gheads, gconst
implicit none
private

public :: get_planning_data_gadget
public :: read_Gpub_particles
public :: read_Gcool_particles
public :: read_Gbromm_particles
public :: update_particles
public :: set_temp_from_u
public :: set_u_from_temp


contains


!>   gets run planning data from Gadget Headers
!===============================================
subroutine get_planning_data_gadget()
  use cosmology_mod, only: tsinceBB

  type(gadget_header_type) :: ghead
  type(gadget_units_type) :: gunits

  integer(i4b) :: iSnap, fSnap    ! initial and final snapshot numbers
  integer(i4b) :: pfiles          ! files/snap for particles    
  integer(i4b) :: i,j             ! counters
  character(clen) :: snapfile     ! snapshot file name

  integer(i4b) :: loglun
  character(clen) :: logfile
  real(r8b) :: kpc2cm
  real(r8b) :: km2cm

  logical :: hdf5bool

  ! set hdf5 boolean
  !======================================================
  hdf5bool = .false.

  ! open up the planning data log file
  !======================================================
  logfile = trim(GV%OutputDir) // "/" // "particle_headers.log"
  call open_formatted_file_w(logfile,loglun)


  ! these global variables are read from the config file
  !======================================================
  iSnap  = GV%StartSnapNum
  fSnap  = GV%EndSnapNum
  pfiles = GV%ParFilesPerSnap

  if ( allocated(saved_gheads) ) deallocate(saved_gheads)
  allocate( saved_gheads(iSnap:fSnap, 0:pfiles-1) )
 
  call read_gadget_constants( "dummysnapname", gconst, hdf5bool, "default" ) 
  call read_gadget_units( "dummysnapname", gunits, hdf5bool, "default" )  

  ! read all particle headers and write to log file
  !===================================================
  write(loglun,'(A)') "reading all GADGET particle header(s) ... "
  do i = iSnap,fSnap
     do j = 0,pfiles-1


        call form_gadget_snapshot_file_name(GV%SnapPath, GV%ParFileBase, i, j, snapfile, hdf5bool)
        write(loglun,'(I3,"  ",A)') i, trim(snapfile)

        call read_gadget_header_file(snapfile,ghead)
        call gadget_header_to_file(ghead,loglun)
        saved_gheads(i,j) = ghead

        ! make sure there is gas in this snapshot
        !-------------------------------------------
        if (.not. ghead%npar_all(1) > 0) then
           write(*,*) "Gadget snapshot does not contain any gas particles"
           write(*,*) "Sphray cannot read dark matter particles directly, "
           write(*,*) "please calculate smoothing lengths for these particles"
           write(*,*) "and write them as gas particles. "
           stop
        end if

        ! IsoTemp comes from config file

        GV%BoxLwrsComoh(:) = 0.0d0
        GV%BoxUprsComoh(:) = ghead%boxlen

        GV%OmegaM = ghead%OmegaM
        GV%OmegaL = ghead%OmegaL

        GV%OmegaB = 0.04  ! this should be moved to the config file
                          ! but its not used in the code now

        GV%LittleH = ghead%h

        GV%cgs_len  = gunits%len
        GV%cgs_mass = gunits%mass
        GV%cgs_vel  = gunits%vel

        GV%cgs_time = gunits%time
        GV%cgs_rho  = gunits%rho
        GV%cgs_prs  = gunits%prs
        GV%cgs_enrg = gunits%energy


        if (GV%Comoving) then
           PLAN%snap(i)%ScalefacAt = ghead%a
           PLAN%snap(i)%TimeAt = tsinceBB(ghead%a, GV%OmegaM, GV%LittleH)       ! in seconds
           PLAN%snap(i)%TimeAt = PLAN%snap(i)%TimeAt * GV%LittleH / GV%cgs_time ! in code units
        else
           PLAN%snap(i)%ScalefacAt = 1.0d0 / (1.0d0 + ghead%z)
           PLAN%snap(i)%TimeAt = ghead%a
        end if

     end do
  end do


  ! close headers log file
  !========================
  close(loglun)


  ! write units to log file
  !===================================================

  logfile = trim(GV%OutputDir) // "/" // "code_units.log"
  call open_formatted_file_w(logfile,loglun)

  kpc2cm = gconst%cm_per_mpc * 1.0d-3
  km2cm = 1.0d5

  106    format(A,ES12.5,A)
  write(loglun,*) 
  write(loglun,'(A)') "setting code units ..."
  write(loglun,106) "  Hubble:     ", GV%LittleH, " = H0[km/s/Mpc] / 100 "
  write(loglun,106) "  length:     ", GV%cgs_len / kpc2cm, " [kpc/h]"
  write(loglun,106) "  mass:       ", GV%cgs_mass / gconst%solar_mass, " [Msun/h]"
  write(loglun,106) "  velocity    ", GV%cgs_vel / km2cm, " [km/s]"
  write(loglun,106) "  time:       ", GV%cgs_time / gconst%sec_per_megayear, " [Myr/h]"
  write(loglun,106) "  density:    ", GV%cgs_rho / (gconst%solar_mass / kpc2cm**3), " [Msun/kpc^3 h^2]"
  write(loglun,106) "  pressure:   ", GV%cgs_prs, " [dyne/cm^2 h^2]"
  write(loglun,106) "  energy:     ", GV%cgs_enrg, " [ergs/h]"
  write(loglun,*) 

  close(loglun)


end subroutine get_planning_data_gadget



!> reads a Gadget-2 (public version) snapshot into a particle array  
!========================================================================
subroutine read_Gpub_particles()

  character(clen), parameter :: myname="read_Gpub_particles"
  logical, parameter :: crash=.true.
  integer, parameter :: verb=2
  character(clen) :: str,fmt

  real(r4b), allocatable :: rblck(:)
  real(r4b), allocatable :: rblck3(:,:)
  integer(i4b), allocatable :: iblck(:)
  integer(i8b) :: ngasread

  type(gadget_header_type) :: ghead
  character(clen) :: snapfile
  integer(i4b) :: lun
  integer(i8b) :: i
  integer(i4b) :: err

  integer(i8b) :: npar, ngas, nmass
  integer(i8b) :: npar1, ngas1, nmass1
  logical :: varmass(6)
  integer(i4b) :: fn

  real(r8b) :: meanweight
  logical :: caseA(2)
  real(r8b) :: MB 

  logical :: hdf5bool

  ! set hdf5 boolean
  !======================================================
  hdf5bool = .false.

  ! set local particle numbers
  !============================
  ghead = saved_gheads( GV%CurSnapNum, 0 )
  varmass = (ghead%npar_all > 0 .and. ghead%mass == 0)
  npar = sum(ghead%npar_all)
  ngas = ghead%npar_all(1)
  nmass = sum(ghead%npar_all, mask=varmass)

  ! do Gadget dummy checks
  !============================
  if (ngas == 0) call myerr("snapshot has no gas particles",myname,crash)

  ! calculate bytes per particle and allocate particle array
  !===========================================================
  MB = GV%bytesperpar * real(ngas) / 2.**20
  GV%MB = GV%MB + MB

  fmt="(A,F10.4,A,I10,A)"
  write(str,fmt) "  allocating ", MB, " MB for ", ngas, " particles"
  call mywrite(str,verb)

  allocate (psys%par(ngas), stat=err)
  if (err /= 0) call myerr("failed to allocate par",myname,crash)

  ! now read all snapshot files
  !==============================          

  ngasread = 0
  files: do fn = 0, ghead%nfiles-1

     ! recall the header info
     !-----------------------------------------------------------!  
     ghead   = saved_gheads( GV%CurSnapNum, fn )
     varmass = (ghead%npar_file > 0 .and. ghead%mass == 0)
     npar1   = sum(ghead%npar_file)
     ngas1   = ghead%npar_file(1)
     nmass1  = sum(ghead%npar_file, mask=varmass)
     if (ngas1 == 0) cycle

     ! begin read
     !-----------------------------------------------------------!  
     call form_gadget_snapshot_file_name(GV%SnapPath,GV%ParFileBase,GV%CurSnapNum,fn,snapfile,hdf5bool)
     call mywrite("   reading public gadget particle snapshot file "//trim(snapfile), verb,fmt="(A)")
     call open_unformatted_file_r( snapfile, lun )
     call read_gadget_header_lun(lun,ghead)


     ! read positions
     !-----------------------------------------------------------!  
     allocate(rblck3(3,ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck3 for pos",myname,crash)
     read(lun, iostat=err) rblck3
     if (err/=0) call myerr("reading rblk3 for pos",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%pos(1) = rblck3(1,i)
     forall(i=1:ngas1) psys%par(ngasread+i)%pos(2) = rblck3(2,i)
     forall(i=1:ngas1) psys%par(ngasread+i)%pos(3) = rblck3(3,i)
     deallocate(rblck3)

     ! read velocities 
     !-----------------------------------------------------------!  
#ifdef incVel
     allocate(rblck3(3,ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck3 for vel",myname,crash)
     read(lun, iostat=err) rblck3
     if (err/=0) call myerr("reading rblk3 for vel",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%vel(1) = rblck3(1,i)
     forall(i=1:ngas1) psys%par(ngasread+i)%vel(2) = rblck3(2,i)
     forall(i=1:ngas1) psys%par(ngasread+i)%vel(3) = rblck3(3,i)
     deallocate(rblck3)
#else
     read(lun, iostat=err)
     if (err/=0) call myerr("dummy reading vel",myname,crash) 
#endif

     ! read id's 
     !-----------------------------------------------------------!  
     allocate(iblck(ngas1), stat=err )
     if(err/=0) call myerr("allocating iblck for ID",myname,crash)
     read(lun, iostat=err) iblck  
     if (err/=0) call myerr("reading iblk for ID",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%id = iblck(i)
     deallocate(iblck)

     ! read masses 
     !-----------------------------------------------------------!  

     ! if there are variable mass particles of any type
     if (nmass1 > 0) then  

        ! if there are variable mass gas particles
        if (varmass(1)) then  
           allocate(rblck(ngas1), stat=err)
           if(err/=0) call myerr("allocating rblck for mass",myname,crash)
           read(lun, iostat=err) rblck 
           if (err/=0) call myerr("reading rblk for mass",myname,crash) 
           forall(i=1:ngas1) psys%par(ngasread+i)%mass = rblck(i)
           deallocate(rblck)

        ! just dummy read non gas particles
        else 
           read(lun)  
           psys%par(ngasread+1:ngasread+ngas1)%mass = ghead%mass(1)
        end if

     ! if none of the particles are variable mass
     else  
        psys%par(ngasread+1:ngasread+ngas1)%mass = ghead%mass(1)
     end if

     ! read temperature (internal energy / unit mass for now)
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for u",myname,crash)
     read(lun, iostat=err) rblck
     if (err/=0) call myerr("reading rblk for u",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%T = rblck(i) 
     deallocate(rblck)


     ! read density 
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for rho",myname,crash)
     read(lun, iostat=err) rblck
     if (err/=0) call myerr("reading rblk for rho",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%rho = rblck(i)
     deallocate(rblck)

     ! read smoothing lengths 
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for hsml",myname,crash)
     read(lun, iostat=err) rblck
     if (err/=0) call myerr("reading rblk for hsml",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%hsml = rblck(i)
     deallocate(rblck)


     ngasread = ngasread + ngas1
     close(lun)

  end do files


  ! There is no cooling in the public version of Gadget-2 so we
  ! will make some simple assumptions to calculate temperature 
  ! from the internal energy / unit mass


  ! assume the particles are neutral to compute the temperature
  !---------------------------------------------------------------
  psys%par(:)%xHI = 1.0d0 - 1.0d-5
  psys%par(:)%xHII = 1.0d-5
  psys%par(:)%ye = psys%par(:)%xHII

  call set_temp_from_u(psys, GV%H_mf, GV%cgs_enrg, GV%cgs_mass)



  ! set caseA true or false for collisional equilibrium
  !-----------------------------------------------------
  caseA = .false.
  if (.not. GV%OnTheSpotH  .or. GV%HydrogenCaseA) caseA(1) = .true.
  if (.not. GV%OnTheSpotHe .or. GV%HeliumCaseA)   caseA(2) = .true.

  call set_collisional_ionization_equilibrium(psys, caseA, GV%IsoTemp, DoHydrogen=.true., fit="hui")
  call set_ye(psys, GV%H_mf, GV%He_mf, GV%NeBackGround)




end subroutine read_Gpub_particles

 

!> reads a Gadget snapshot that includes ye and xHI (i.e. cooling) into a particle array  
!=========================================================================================
subroutine read_Gcool_particles()

  character(clen), parameter :: myname="read_Gcool_particles"
  logical, parameter :: crash=.true.
  integer, parameter :: verb=2
  character(clen) :: str,fmt 

  real(r4b), allocatable :: rblck(:)
  real(r4b), allocatable :: rblck3(:,:)
  integer(i4b), allocatable :: iblck(:)
  integer(i8b) :: ngasread

  type(gadget_header_type) :: ghead
  character(clen) :: snapfile
  integer(i4b) :: lun,i
  integer(i4b) :: err

  integer(i8b) :: npar, ngas, nmass
  integer(i8b) :: npar1, ngas1, nmass1
  logical :: varmass(6)
  integer(i4b) :: fn

  real(r8b) :: meanweight
  logical :: caseA(2)
  real(r8b) :: xvec(5)
  real(r8b) :: Tdum
  real(r8b) :: MB 

  logical :: hdf5bool

  ! set hdf5 boolean
  !======================================================
  hdf5bool = .false.

  ! set local particle numbers
  !============================
  ghead = saved_gheads( GV%CurSnapNum, 0 )
  varmass = (ghead%npar_all > 0 .and. ghead%mass == 0)
  npar = sum(ghead%npar_all)
  ngas = ghead%npar_all(1)
  nmass = sum(ghead%npar_all, mask=varmass)

  ! do Gadget dummy checks
  !============================
  if (ngas .EQ. 0) call myerr("snapshot has no gas particles",myname,crash)

  ! calculate bytes per particle and allocate particle array
  !===========================================================
  MB = GV%bytesperpar * real(ngas) / 2.**20
  GV%MB = GV%MB + MB

  fmt="(A,F10.4,A,I10,A)"
  write(str,fmt) "allocating ", MB, " MB for ", ngas, " particles"
  call mywrite(str,verb)

  allocate (psys%par(ngas), stat=err)
  if (err /= 0) call myerr("failed to allocate par",myname,crash)



  ! now read all snapshot files
  !==============================          
  ngasread = 0
  files: do fn = 0, ghead%nfiles-1

     ! recall the header info
     !-----------------------------------------------------------!  
     ghead   = saved_gheads( GV%CurSnapNum, fn )
     varmass = (ghead%npar_file > 0 .and. ghead%mass == 0)
     npar1   = sum(ghead%npar_file)
     ngas1   = ghead%npar_file(1)
     nmass1  = sum(ghead%npar_file, mask=varmass)
     if (ngas1 == 0) cycle

     ! begin read
     !-----------------------------------------------------------!  
     call form_gadget_snapshot_file_name(GV%SnapPath,GV%ParFileBase,GV%CurSnapNum,fn,snapfile,hdf5bool)
     call mywrite("reading gadget w/ cooling snapshot file "//trim(snapfile), verb)
     call open_unformatted_file_r( snapfile, lun )
     read(lun) ghead


     ! read positions 
     !-----------------------------------------------------------!  
     allocate(rblck3(3,ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck3 for pos",myname,crash)
     read(lun, iostat=err) rblck3
     if (err/=0) call myerr("reading rblk3 for pos",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%pos(1) = rblck3(1,i)
     forall(i=1:ngas1) psys%par(ngasread+i)%pos(2) = rblck3(2,i)
     forall(i=1:ngas1) psys%par(ngasread+i)%pos(3) = rblck3(3,i)
     deallocate(rblck3)

     ! read velocities 
     !-----------------------------------------------------------!  
#ifdef incVel
     allocate(rblck3(3,ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck3 for vel",myname,crash)
     read(lun, iostat=err) rblck3
     if (err/=0) call myerr("reading rblk3 for vel",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%vel(1) = rblck3(1,i)
     forall(i=1:ngas1) psys%par(ngasread+i)%vel(2) = rblck3(2,i)
     forall(i=1:ngas1) psys%par(ngasread+i)%vel(3) = rblck3(3,i)
     deallocate(rblck3)
#else
     read(lun, iostat=err)
     if (err/=0) call myerr("dummy reading vel",myname,crash) 
#endif


     ! read id's 
     !-----------------------------------------------------------!  
     allocate(iblck(ngas1), stat=err )
     if(err/=0) call myerr("allocating iblck for ID",myname,crash)
     read(lun, iostat=err) iblck  
     if (err/=0) call myerr("reading iblk for ID",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%id = iblck(i)
     deallocate(iblck)

     ! read masses 
     !-----------------------------------------------------------!  

     ! if there are variable mass particles of any type
     if (nmass1 > 0) then  

        ! if there are variable mass gas particles
        if (varmass(1)) then  
           allocate(rblck(ngas1), stat=err)
           if(err/=0) call myerr("allocating rblck for mass",myname,crash)
           read(lun, iostat=err) rblck 
           if (err/=0) call myerr("reading rblk for mass",myname,crash) 
           forall(i=1:ngas1) psys%par(ngasread+i)%mass = rblck(i)
           deallocate(rblck)

        ! just dummy read non gas particles
        else 
           read(lun)  
           psys%par(ngasread+1:ngasread+ngas1)%mass = ghead%mass(1)
        end if

     ! if none of the particles are variable mass
     else  
        psys%par(ngasread+1:ngasread+ngas1)%mass = ghead%mass(1)
     end if


     ! read temperature (internal energy / unit mass for now)
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for u",myname,crash)
     read(lun, iostat=err) rblck
     if (err/=0) call myerr("reading rblk for u",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%T = rblck(i) 
     deallocate(rblck)


     ! read density 
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for rho",myname,crash)
     read(lun, iostat=err) rblck
     if (err/=0) call myerr("reading rblk for rho",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%rho = rblck(i)
     deallocate(rblck)

     ! read electron fraction ye = ne/nH
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for ye",myname,crash)
     read(lun, iostat=err) rblck
     if (err/=0) call myerr("reading rblck for ye",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%ye = rblck(i)
     deallocate(rblck)

     ! read neutral hydrogen fraction (xHI)
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for xHI",myname,crash)
     read(lun, iostat=err) rblck
     if (err/=0) call myerr("reading rblk for xHI",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%xHI = rblck(i)
     deallocate(rblck)

     ! read smoothing lengths 
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for hsml",myname,crash)
     read(lun, iostat=err) rblck
     if (err/=0) call myerr("reading rblk for hsml",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%hsml = rblck(i)
     deallocate(rblck)


     ! the following data blocks may be in the snapshot file
     ! but we make no use of them for the time being

     ! star formation rate (ngas)
     ! stellar age (nstar)
     ! metallicity (ngas+nstar)
     ! black hole mass (nbh)
     ! black hole accretion rate (nbh)

     ngasread = ngasread + ngas1
     close(lun)

  end do files


  ! convert the internal energy to temperature using ye from snap
  !----------------------------------------------------------------!  
  call set_temp_from_u(psys, GV%H_mf, GV%cgs_enrg, GV%cgs_mass)


  ! set xHII from xHI 
  !----------------------------------------------------------------!  
  psys%par%xHII = 1.0d0 - psys%par%xHI


  ! set caseA true or false for collisional equilibrium
  !-----------------------------------------------------
  caseA = .false.
  if (.not. GV%OnTheSpotH  .or. GV%HydrogenCaseA) caseA(1) = .true.
  if (.not. GV%OnTheSpotHe .or. GV%HeliumCaseA)   caseA(2) = .true.


  ! if Helium, initialize ionization fractions to collisional equilibrium
  !------------------------------------------------------------------------
#ifdef incHe
  call set_collisional_ionization_equilibrium(psys, caseA, GV%IsoTemp, DoHydrogen=.false., fit="hui")
#endif


  ! set the electron fractions from the ionization fractions
  !----------------------------------------------------------
  !  call set_ye(psys, GV%H_mf, GV%He_mf)



end subroutine read_Gcool_particles




!> reads a Gadget snapshot from Volker Bromm's group
!=========================================================================================
subroutine read_Gbromm_particles()

  character(clen), parameter :: myname="read_Gbromm_particles"
  logical, parameter :: crash=.true.
  integer, parameter :: verb=2
  character(clen) :: str,fmt 

  real(r4b), allocatable :: rblck(:)
  real(r4b), allocatable :: rblck3(:,:)
  integer(i4b), allocatable :: iblck(:)
  integer(i8b) :: ngasread

  type(gadget_header_type) :: ghead
  character(clen) :: snapfile
  integer(i4b) :: lun,i
  integer(i4b) :: err

  integer(i8b) :: npar, ngas, nmass
  integer(i8b) :: npar1, ngas1, nmass1
  logical :: varmass(6)
  integer(i4b) :: fn

  real(r8b) :: meanweight
  real(r8b) :: nH_over_nHe
  logical :: caseA(2)
  real(r8b) :: xvec(5)
  real(r8b) :: Tdum
  real(r8b) :: MB 

  logical :: hdf5bool

  ! set hdf5 boolean
  !======================================================
  hdf5bool = .false.

  ! set local particle numbers
  !============================
  ghead = saved_gheads( GV%CurSnapNum, 0 )
  varmass = (ghead%npar_all > 0 .and. ghead%mass == 0)
  npar = sum(ghead%npar_all)
  ngas = ghead%npar_all(1)
  nmass = sum(ghead%npar_all, mask=varmass)

  ! do Gadget dummy checks
  !============================
  if (ngas .EQ. 0) call myerr("snapshot has no gas particles",myname,crash)

  ! calculate bytes per particle and allocate particle array
  !===========================================================
  MB = GV%bytesperpar * real(ngas) / 2.**20
  GV%MB = GV%MB + MB

  fmt="(A,F10.4,A,I10,A)"
  write(str,fmt) "allocating ", MB, " MB for ", ngas, " particles"
  call mywrite(str,verb)

  allocate (psys%par(ngas), stat=err)
  if (err /= 0) call myerr("failed to allocate par",myname,crash)



  ! now read all snapshot files
  !==============================          
  ngasread = 0
  files: do fn = 0, ghead%nfiles-1

     ! recall the header info
     !-----------------------------------------------------------!  
     ghead   = saved_gheads( GV%CurSnapNum, fn )
     varmass = (ghead%npar_file > 0 .and. ghead%mass == 0)
     npar1   = sum(ghead%npar_file)
     ngas1   = ghead%npar_file(1)
     nmass1  = sum(ghead%npar_file, mask=varmass)
     if (ngas1 == 0) cycle

     ! begin read
     !-----------------------------------------------------------!  
     call form_gadget_snapshot_file_name(GV%SnapPath,GV%ParFileBase,GV%CurSnapNum,fn,snapfile,hdf5bool)
     call mywrite("reading bromm gadget snapshot file "//trim(snapfile), verb)
     call open_unformatted_file_r( snapfile, lun )
     read(lun) ghead


     ! read positions 
     !-----------------------------------------------------------!  
     allocate(rblck3(3,ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck3 for pos",myname,crash)
     read(lun, iostat=err) rblck3
     if (err/=0) call myerr("reading rblk3 for pos",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%pos(1) = rblck3(1,i)
     forall(i=1:ngas1) psys%par(ngasread+i)%pos(2) = rblck3(2,i)
     forall(i=1:ngas1) psys%par(ngasread+i)%pos(3) = rblck3(3,i)
     deallocate(rblck3)

     ! read velocities 
     !-----------------------------------------------------------!  
#ifdef incVel
     allocate(rblck3(3,ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck3 for vel",myname,crash)
     read(lun, iostat=err) rblck3
     if (err/=0) call myerr("reading rblk3 for vel",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%vel(1) = rblck3(1,i)
     forall(i=1:ngas1) psys%par(ngasread+i)%vel(2) = rblck3(2,i)
     forall(i=1:ngas1) psys%par(ngasread+i)%vel(3) = rblck3(3,i)
     deallocate(rblck3)
#else
     read(lun, iostat=err)
     if (err/=0) call myerr("dummy reading vel",myname,crash) 
#endif


     ! read id's 
     !-----------------------------------------------------------!  
     allocate(iblck(ngas1), stat=err )
     if(err/=0) call myerr("allocating iblck for ID",myname,crash)
     read(lun, iostat=err) iblck  
     if (err/=0) call myerr("reading iblk for ID",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%id = iblck(i)
     deallocate(iblck)

     ! read masses 
     !-----------------------------------------------------------!  

     ! if there are variable mass particles of any type
     if (nmass1 > 0) then  

        ! if there are variable mass gas particles
        if (varmass(1)) then  
           allocate(rblck(ngas1), stat=err)
           if(err/=0) call myerr("allocating rblck for mass",myname,crash)
           read(lun, iostat=err) rblck 
           if (err/=0) call myerr("reading rblk for mass",myname,crash) 
           forall(i=1:ngas1) psys%par(ngasread+i)%mass = rblck(i)
           deallocate(rblck)

        ! just dummy read non gas particles
        else 
           read(lun)  
           psys%par(ngasread+1:ngasread+ngas1)%mass = ghead%mass(1)
        end if

     ! if none of the particles are variable mass
     else  
        psys%par(ngasread+1:ngasread+ngas1)%mass = ghead%mass(1)
     end if


     ! read temperature (internal energy / unit mass for now)
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for u",myname,crash)
     read(lun, iostat=err) rblck
     if (err/=0) call myerr("reading rblk for u",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%T = rblck(i) 
     deallocate(rblck)


     ! read density 
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for rho",myname,crash)
     read(lun, iostat=err) rblck
     if (err/=0) call myerr("reading rblk for rho",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%rho = rblck(i)
     deallocate(rblck)


     ! read smoothing lengths 
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for hsml",myname,crash)
     read(lun, iostat=err) rblck
     if (err/=0) call myerr("reading rblk for hsml",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%hsml = rblck(i)
     deallocate(rblck)
     
     ! CAFG: First, big array containing all the abundances.
     allocate(rblck(10*ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for all abundances",myname,crash)
     read(lun, iostat=err) rblck
     if (err/=0) call myerr("reading rblk for all abundances",myname,crash) 
     
     ! CAFG: Here is where we should assign the desired abundances to the
     ! relevant particle data.
     forall(i=1:ngas1) psys%par(ngasread+i)%xHII = rblck(10*(i-1)+2)
     
 
#ifdef incHe
     ! read HeII = nHeII / nH
     !-----------------------------------------------------------!  
     forall(i=1:ngas1) psys%par(ngasread+i)%xHeII = rblck(10*(i-1)+9)


     ! read HeIII = nHeIII / nH
     !-----------------------------------------------------------!  
     forall(i=1:ngas1) psys%par(ngasread+i)%xHeIII = rblck(10*(i-1)+10)

#endif
     
     deallocate(rblck)

     ngasread = ngasread + ngas1
     close(lun)

  end do files



  ! set xHI from xHII 
  !----------------------------------------------------------------!  
  psys%par(:)%xHI = 1.0d0 - psys%par(:)%xHII


  ! use Hydrogen and Helium mass fractions from config file to
  ! set the Helium ionization fractions
  !----------------------------------------------------------------!  
#ifdef incHe
  nH_over_nHe = 4 * GV%H_mf / GV%He_mf
  psys%par(:)%xHeII  = psys%par(:)%xHeII  * nH_over_nHe
  psys%par(:)%xHeIII = psys%par(:)%xHeIII * nH_over_nHe
  psys%par(:)%xHeI   = 1.0d0 - psys%par(:)%xHeII - psys%par(:)%xHeIII
#endif 

  ! set ye from ionization fractions
  !----------------------------------------------------------------!  
  psys%par(:)%ye = psys%par(:)%xHII
#ifdef incHe
  psys%par(:)%ye = psys%par(:)%ye + psys%par(:)%xHeII + 2 * psys%par(:)%xHeIII
#endif

  ! convert the internal energy to temperature using ye 
  !----------------------------------------------------------------!  
  
  call set_temp_from_u_bromm(psys, GV%H_mf, GV%cgs_enrg, GV%cgs_mass)


  ! shrink particles with negative IDs to 
  !----------------------------------------------------------------!  
  where( psys%par(:)%id < 0 ) psys%par(:)%hsml = 0.0




end subroutine read_Gbromm_particles



!> Reads in an update snapshot.  The idea is to keep the ionization fractions
!! and temperature from the already loaded snapshot while updating the 
!! positions, velocities, smoothing lengths ... from the snapshot on disk.
!! The particles are reordered during the raytracing so when a new snapshot
!! is readin they must be matched
!===========================================================================
subroutine update_particles()

  character(clen), parameter :: myname="update_particles" 
  logical, parameter :: crash=.true.

  integer(i4b), allocatable :: idold(:)
  real(r4b), allocatable :: Told(:)
  real(r4b), allocatable :: yeold(:)
  real(r4b), allocatable :: xHIold(:)
  real(r4b), allocatable :: xHIIold(:)

  real(r4b), allocatable :: xHeIold(:)
  real(r4b), allocatable :: xHeIIold(:)
  real(r4b), allocatable :: xHeIIIold(:)

  integer(i8b), allocatable :: lasthitold(:)
  integer(i8b), allocatable :: orderold(:)  

  integer(i8b) :: minIDold, minIDnew
  integer(i8b) :: maxIDold, maxIDnew

  integer(i8b) :: i, idnew, indx, err




  ! store carry over variables
  !==============================

  allocate( idold(size(psys%par)), stat=err)
  if (err /= 0) call myerr("failed to allocate IDold",myname,crash)     
  idold = psys%par%id
  minIDold = minval(idold)
  maxIDold = maxval(idold)

  allocate( Told(size(psys%par)), stat=err)
  if (err /= 0) call myerr("failed to allocate Told",myname,crash)     
  Told = psys%par%T

  allocate( yeold(size(psys%par)), stat=err)
  if (err /= 0) call myerr("failed to allocate yeold",myname,crash)     
  yeold = psys%par%ye

  allocate( xHIold(size(psys%par)), stat=err)
  if (err /= 0) call myerr("failed to allocate xHIold",myname,crash)     
  xHIold = psys%par%xHI

  allocate( xHIIold(size(psys%par)), stat=err)
  if (err /= 0) call myerr("failed to allocate xHIIold",myname,crash)     
  xHIIold = psys%par%xHII

#ifdef incHe
  allocate( xHeIold(size(psys%par)), stat=err)
  if (err /= 0) call myerr("failed to allocate xHeIold",myname,crash)     
  xHeIold = psys%par%xHeI

  allocate( xHeIIold(size(psys%par)), stat=err)
  if (err /= 0) call myerr("failed to allocate xHeIIold",myname,crash)     
  xHeIIold = psys%par%xHeII

  allocate( xHeIIIold(size(psys%par)), stat=err)
  if (err /= 0) call myerr("failed to allocate xHeIIIold",myname,crash)     
  xHeIIIold = psys%par%xHeIII
#else
  allocate( xHeIold(1), xHeIIold(1), xHeIIIold(1) )
#endif



  allocate( lasthitold(size(psys%par)), stat=err)
  if (err /= 0) call myerr("failed to allocate lasthitold",myname,crash)     
  lasthitold = psys%par%lasthit


  ! deallocate the current particle array and read new particle data
  !==================================================================
  deallocate( psys%par )
  if (GV%InputType == 1) then
     call read_Gpub_particles()
  else if (GV%InputType == 2) then
     call read_Gcool_particles()
  else if (GV%InputType == 3) then
     call read_Ghdf5_particles()
  else if (GV%InputType == 4) then
     call read_Gbromm_particles()
  end if
  minIDnew = minval(psys%par%id)
  maxIDnew = maxval(psys%par%id)


  ! match id's from new snapshot data to the old IDs
  !===================================================


  ! create an indexing array for the old particle array ID's
  ! ID's that are no longer present will be indexed with minIDold - 1
  !--------------------------------------------------------------------
  write(*,*) "creating indexing array for current particle system..."
  write(*,*) "min/max ID (old snap) = ", minIDold, maxIDold

  allocate( orderold(minIDold:maxIDold), stat=err)
  if (err/=0) call myerr("allocating orderold array",myname,crash)

  orderold = minIDold - 1
  do i = 1,size(idold)
     orderold(idold(i)) = i
  end do

  ! transfer all carryover variables
  !-----------------------------------
  do i = 1,size(psys%par)

     idnew = psys%par(i)%id   
     indx  = orderold(idnew)

     if ( idnew /= idold(indx) ) then
        write(*,*) "i,idnew,idold", i, idnew, idold(indx)
        if (err/=0) call myerr("reordering id error",myname,crash)
     end if

     if (.not. GV%FixSnapTemp) then
        psys%par(i)%T = Told(indx)
     end if

     psys%par(i)%ye   = yeold(indx) 
     psys%par(i)%xHI  = xHIold(indx)
     psys%par(i)%xHII = xHIIold(indx)

#ifdef incHe
     psys%par(i)%xHeI   = xHeIold(indx)
     psys%par(i)%xHeII  = xHeIIold(indx)
     psys%par(i)%xHeIII = xHeIIIold(indx)
#endif

     psys%par(i)%lasthit = lasthitold(indx)

  end do

  deallocate (idold, Told, yeold, xHIold, xHIIold, lasthitold, orderold)

  deallocate (xHeIold, xHeIIold, xHeIIIold)



end subroutine update_particles


!> converts internal energies / unit mass to temperature K using Hmf and ye = ne/nH
!=======================================================================================
subroutine set_temp_from_u(psys, dfltH_mf, cgs_enrg, cgs_mass)

  type(particle_system_type) :: psys
  real(r8b), intent(in) :: dfltH_mf
  real(r8b), intent(in) :: cgs_enrg
  real(r8b), intent(in) :: cgs_mass
  integer(i8b) :: i
  real(r8b) :: Hmf
  real(r8b) :: mu
  real(r8b) :: Tdum

  do i = 1,size(psys%par)

#ifdef incHmf
     Hmf = psys%par(i)%Hmf
#else     
     Hmf = dfltH_mf    
#endif

     mu = 4.0d0 / (3.0d0 * Hmf + 1.0d0 + 4.0d0 * Hmf * psys%par(i)%ye)
     Tdum = mu * gconst%PROTONMASS / gconst%BOLTZMANN * (gconst%GAMMA - 1.0d0) * psys%par(i)%T
     psys%par(i)%T = Tdum * cgs_enrg / cgs_mass

  end do


end subroutine set_temp_from_u


!> converts internal energies / unit mass to temperature K using Hmf and ye = ne/nH
!=======================================================================================
subroutine set_temp_from_u_bromm(psys, dfltH_mf, cgs_enrg, cgs_mass)

  type(particle_system_type) :: psys
  real(r8b), intent(in) :: dfltH_mf
  real(r8b), intent(in) :: cgs_enrg
  real(r8b), intent(in) :: cgs_mass
  integer(i8b) :: i
  real(r8b) :: Hmf
  real(r8b) :: mu
  real(r8b) :: Tdum

  do i = 1,size(psys%par)

#ifdef incHmf
     Hmf = psys%par(i)%Hmf
#else
     
     Hmf = dfltH_mf    

#endif

     mu = 4.0d0 / (3.0d0 * Hmf + 1.0d0 + 4.0d0 * Hmf * psys%par(i)%ye)

     ! CAFG: artificially 'correct' abnormally large values, which are
     ! producing floating overflow
     if (psys%par(i)%T .GT. 1000000000.0d0) psys%par(i)%T = 1000.0d0

     Tdum = mu * gconst%PROTONMASS / gconst%BOLTZMANN * (gconst%GAMMA - 1.0d0) * psys%par(i)%T
     psys%par(i)%T = Tdum * cgs_enrg / cgs_mass

  end do


end subroutine set_temp_from_u_bromm



!> converts temperature K to internal energies / unit mass 
!=======================================================================================
subroutine set_u_from_temp(psys, dfltH_mf, cgs_enrg, cgs_mass)

  type(particle_system_type) :: psys
  real(r8b), intent(in) :: dfltH_mf
  real(r8b), intent(in) :: cgs_enrg
  real(r8b), intent(in) :: cgs_mass
  integer(i8b) :: i
  real(r8b) :: Hmf
  real(r8b) :: mu
  real(r8b) :: Udum

  do i = 1,size(psys%par)

#ifdef incHmf
     Hmf = psys%par(i)%Hmf
#else
     Hmf = dfltH_mf
#endif

     mu = 4.0d0 / (3.0d0 * Hmf + 1.0d0 + 4.0d0 * Hmf * psys%par(i)%ye)
     Udum = gconst%BOLTZMANN * psys%par(i)%T /( (gconst%GAMMA - 1.0d0) * mu * gconst%PROTONMASS )
     psys%par(i)%T = Udum * cgs_mass / cgs_enrg

  end do


end subroutine set_u_from_temp












end module gadget_input_mod
