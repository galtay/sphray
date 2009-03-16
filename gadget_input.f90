!> \file gadget_input.f90

!> \brief Handles readin of GADGET formatted files
!<

module gadget_input_mod
use myf90_mod
use particle_system_mod
use source_input_mod, only: source_header_type
use global_mod, only: global_variables_type
use atomic_rates_mod, only: calc_colion_eq_table
use atomic_rates_mod, only: calc_colion_eq_fits
use physical_constants_mod
use global_mod, only: psys, PLAN
implicit none

!> gadget particle types
character(5), parameter :: ptype_names(6) = (/"gas  ","halo ","disk ",&
                                              "bulge","stars","bndry"/)
  
real(r8b), parameter :: GLEN = 3.085678d21  !< in cm/h (1 kpc/h)
real(r8b), parameter :: GMASS = 1.989d43    !< in g/h  (10^10 Solar/h)
real(r8b), parameter :: GVEL = 1.0d5        !< in km/s 

real(r8b), parameter :: GTIME = GLEN / GVEL                !< gadget time unit
real(r8b), parameter :: GRHO = GMASS / GLEN**3             !< gadget density unit
real(r8b), parameter :: GPRS = GMASS / GLEN / GTIME**2     !< gadget pressure unit
real(r8b), parameter :: GENRG = GMASS * GLEN**2 / GTIME**2 !< gadget energy unit
real(r8b), parameter :: GLUM = GENRG / GTIME               !< gadget luminosity unit

real(r8b), parameter :: GMASS_solar = GMASS / M_solar

real(r8b), parameter, private :: BOLTZMANN = 1.3806d-16   !< boltzmann's constant
real(r8b), parameter, private :: PROTONMASS = 1.6726d-24  !< proton mass (cgs)
real(r8b), parameter, private :: GAMMA = 5.0d0/3.0d0


!> gadget header type
type gadget_header_type
   integer(i4b) :: npar_file(6)    !< number of particles in snapshot file
   real(r8b) :: mass(6)            !< mass of each particle type if constant
   real(r8b) :: a                  !< scale factor or time
   real(r8b) :: z                  !< redshift
   integer(i4b) :: flag_sfr        !< flag for star formation
   integer(i4b) :: flag_feedback   !< flag for feedback
   integer(i4b) :: npar_all(6)     !< number of particles in whole snapshot
   integer(i4b) :: flag_cooling    !< flag for radiative cooling
   integer(i4b) :: nfiles          !< number of files in a this snapshot
   real(r8b) :: boxlen             !< box length
   real(r8b) :: OmegaM             !< omega matter
   real(r8b) :: OmegaL             !< omega lambda
   real(r8b) :: h                  !< little hubble
   integer(i4b) :: flag_age        !< flag for stellar age
   integer(i4b) :: flag_metals     !< flag for metallicity
   integer(i4b) :: npar_hw(6)      !< 64 bit part of npar 
   integer(i4b) :: flag_entr_ics   !< flag for entropic initial conditions
   
   real(r4b)    :: OmegaB          !< omega baryon
   integer(i8b) :: rays_traced     !< number of rays traced so far
   integer(i4b) :: flag_helium     !< is helium present
   integer(i4b) :: flag_gamma      !< is gammaHI present
   integer(i4b) :: unused(9)       !< spacer
end type gadget_header_type


type(gadget_header_type), allocatable :: saved_gheads(:,:)



contains

!>   reads in a particle header, prints the header to the screen if verbose 
!==========================================================================
subroutine read_gadget_header(snapfile,verbose,closefile,ghead,lun)
  character(*), intent(in) :: snapfile  !< snapshot file name
  logical, intent(in) :: verbose        !< report to screen?
  logical, intent(in) :: closefile      !< close file when done
  type(gadget_header_type), intent(out) :: ghead  !< particle header
  integer(i8b), intent(out) :: lun           !< logical unit number used

  call open_unformatted_file_r(snapfile,lun)
  read(lun) ghead
  if (closefile) close(lun)  
  if (verbose) call gadget_header_to_screen(ghead)

end subroutine read_gadget_header


!>   reads box data from a gadget snapshot 
!==========================================================================
subroutine read_gadget_box(GV,box)
  type(global_variables_type), intent(in) :: GV !< global variables
  type(box_type), intent(out) :: box            !< box 

  integer(i8b) :: fn
  character(200) :: snapfile   
  type(gadget_header_type) :: ghead  
  integer(i8b) :: lun    

  fn=0
  call form_Gsnapshot_file_name(GV%SnapPath,GV%ParFileBase,GV%CurSnapNum,fn,snapfile)
  call open_unformatted_file_r(snapfile,lun)
  read(lun) ghead
  close(lun)  

  box%bot = 0.0
  box%top = ghead%boxlen
  box%bbound = GV%BndryCond
  box%tbound = GV%BndryCond

end subroutine read_gadget_box

!>   gets run planning data from Gadget Headers
!===============================================
subroutine get_planning_data_gadget(GV)
use physical_constants_mod
use cosmology_mod, only: tsinceBB
use source_input_mod, only: read_source_header

  type(global_variables_type) :: GV  !< global variables

  logical :: verbose
  logical :: closefile
  integer(i8b) :: lun
  type(gadget_header_type) :: ghead
  type(source_header_type) :: shead

  integer(i8b) :: iSnap, fSnap    ! initial and final snapshot numbers
  integer(i8b) :: pfiles, sfiles  ! files/snap for particles and sources   
  integer(i8b) :: i,j             ! counters
  character(200) :: snapfile ! snapshot file name

  integer(i8b) :: loglun
  character(200) :: logfile


  logfile = trim(GV%OutputDir) // "/" // "headers.log"
  call open_formatted_file_w(logfile,loglun)


  ! these global variables must have already been set
  !====================================================
  iSnap = GV%StartSnapNum
  fSnap = GV%EndSnapNum
    
  pfiles = GV%ParFilesPerSnap
  sfiles = GV%SourceFilesPerSnap
    
  allocate( saved_gheads(iSnap:fSnap, 0:pfiles-1) )

  ! read in and report to screen
  !==============================
  closefile = .true.
  verbose = .false.

  ! read all particle headers and write to log file
  !===================================================
  write(loglun,'(A)') "reading all GADGET particle header(s) ... "
    do i = iSnap,fSnap
       do j = 0,pfiles-1

          call form_Gsnapshot_file_name(GV%SnapPath,GV%ParFileBase,i,j,snapfile)
          write(loglun,'(I3,"  ",A)') i,trim(snapfile)
          call read_gadget_header(snapfile,verbose,closefile,ghead,lun)
          call gadget_header_to_file(ghead,loglun)

          saved_gheads(i,j) = ghead

          ! make sure there is gas in this snapshot
          if (.not. ghead%npar_all(1) > 0) then
             write(*,*) "Gadget snapshot does not contain any gas particles"
             write(*,*) "Sphray cannot read dark matter particles directly, "
             write(*,*) "please use the utilities in the /gadget directory to"
             write(*,*) "convert collisionless snapshot"
             stop
          end if

          ! IsoTemp comes from config file

          GV%BoxLowers(:) = 0.0d0
          GV%BoxUppers(:) = ghead%boxlen

          GV%OmegaM = ghead%OmegaM
          GV%OmegaL = ghead%OmegaL
          GV%OmegaB = 0.04
          GV%LittleH = ghead%h

          GV%cgs_len  = GLEN
          GV%cgs_mass = GMASS
          GV%cgs_vel  = GVEL

          GV%cgs_time = GTIME
          GV%cgs_rho  = GRHO
          GV%cgs_prs  = GPRS
          GV%cgs_enrg = GENRG
          GV%cgs_lum  = GLUM

          if (GV%Comoving) then
             PLAN(i)%ScalefacAtSnap = ghead%a
             PLAN(i)%TimeAtSnap = tsinceBB(ghead%a, GV%OmegaM, GV%LittleH)
             PLAN(i)%TimeAtSnap = PLAN(i)%TimeAtSnap * GV%LittleH / GV%cgs_time
          else
             PLAN(i)%ScalefacAtSnap = 1.0d0
             PLAN(i)%TimeAtSnap = ghead%a
          end if

       end do
    end do

    ! read all source headers and write to log file
    !===================================================
    write(loglun,*) 
    write(loglun,'(A)') "reading all Gadget source header(s) ... "
    do i = iSnap,fSnap
       do j = 1,sfiles
          call form_snapshot_file_name(GV%SourcePath,GV%SourceFileBase,i,j,snapfile)
          write(loglun,'(I3,"  ",A)') i,trim(snapfile)
          call read_source_header(snapfile,verbose,closefile,shead,lun)

          PLAN(i)%RaysFromSrcHeader = shead%TotalRays
          GV%Lunit = shead%Lunit

       end do
    end do

    close(loglun)


  ! write units to log file
  !===================================================

  logfile = trim(GV%OutputDir) // "/" // "code_units.log"
  call open_formatted_file_w(logfile,loglun)

  106    format(A,ES12.5,A)
  write(loglun,*) 
  write(loglun,'(A)') "setting code units ..."
  write(loglun,106) "  Hubble:     ", GV%LittleH, " = H0[km/s/Mpc] / 100 "
  write(loglun,106) "  length:     ", GV%cgs_len / kpc2cm, " [kpc/h]"
  write(loglun,106) "  mass:       ", GV%cgs_mass / M_solar, " [Msun/h]"
  write(loglun,106) "  velocity    ", GV%cgs_vel / km2cm, " [km/s]"
  write(loglun,106) "  time:       ", GV%cgs_time / Myr2sec, " [Myr/h]"
  write(loglun,106) "  density:    ", GV%cgs_rho / (M_solar / kpc2cm**3), " [Msun/kpc^3 h^2]"
  write(loglun,106) "  pressure:   ", GV%cgs_prs, " [dyne/cm^2 h^2]"
  write(loglun,106) "  energy:     ", GV%cgs_enrg, " [ergs/h]"
  write(loglun,106) "  luminosity: ", GV%cgs_lum / L_solar, " [Lsun]"
  write(loglun,*) 
  
  close(loglun)

end subroutine get_planning_data_gadget



!> reads a Gadget-2 (public version) snapshot into a particle array  
!========================================================================
subroutine read_Gpub_particles(GV,MB)
  type(global_variables_type), intent(in) :: GV !< global variables
  real(r8b), intent(out) :: MB !< MB allocated for particle storage

  character(200) :: myname 
  logical :: crash

  real(r4b), allocatable :: rblck(:)
  real(r4b), allocatable :: rblck3(:,:)
  integer(i4b), allocatable :: iblck(:)
  integer(i8b) :: ngasread

  type(gadget_header_type) :: ghead
  character(200) :: snapfile
  integer(i8b) :: lun,i
  integer(i4b) :: err

  integer(i8b) :: npar, ngas, nmass
  integer(i8b) :: npar1, ngas1, nmass1
  integer(i8b) :: bytesperpar
  logical :: closefile
  logical :: header_verbose
  logical :: varmass(6)
  integer(i8b) :: fn

  real(r8b) :: meanweight
  logical :: caseA(2)
  real(r8b) :: xvec(5)
  real(r8b) :: Tdum

  real(r8b) :: nHe_over_nH

  ! set error handling variables
  !===============================
  myname = "read_Gpub_particles"
  crash = .true.

  ! read header from first file
  !============================
  fn=0
  call form_Gsnapshot_file_name(GV%SnapPath,GV%ParFileBase,GV%CurSnapNum,fn,snapfile)
  closefile=.true.
  header_verbose=.false.
  call read_gadget_header(snapfile,header_verbose,closefile,ghead,lun)

  ! set local particle numbers
  !============================
  varmass = (ghead%npar_all > 0 .and. ghead%mass == 0)
  npar = sum(ghead%npar_all)
  ngas = ghead%npar_all(1)
  nmass = sum(ghead%npar_all, mask=varmass)

  ! do Gadget dummy checks
  !============================
  if (ngas .EQ. 0) call myerr("snapshot has no gas particles",myname,crash)

  ! calculate bytes per particle and allocate particle array
  !===========================================================
  call calc_bytes_per_particle(bytesperpar)
  MB = bytesperpar * real(ngas) / 2.**20

  101  format(A,F10.4,A,I10,A)
  write(*,101) "allocating ", MB, " MB for ", ngas, " particles"
  
  allocate (psys%par(ngas), stat=err)
  if (err /= 0) call myerr("failed to allocate par",myname,crash)

  ! now read all snapshot files
  !==============================          

  ngasread = 0
  closefile=.false.
  header_verbose=.false.
  files: do fn = 0, ghead%nfiles-1

     call form_Gsnapshot_file_name(GV%SnapPath,GV%ParFileBase,GV%CurSnapNum,fn,snapfile)
     write(*,'(A,A)') "reading public gadget particle snapshot file ", trim(snapfile)
     call read_gadget_header(snapfile,header_verbose,closefile,ghead,lun)

     varmass = (ghead%npar_file > 0 .and. ghead%mass == 0)
     npar1 = sum(ghead%npar_file)
     ngas1 = ghead%npar_file(1)
     nmass1 = sum(ghead%npar_file, mask=varmass)

     if (ngas1 == 0) cycle


     ! read positions
     !-----------------------------------------------------------!  
     allocate(rblck3(3,ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck3 for pos",myname,crash)
     read(lun, iostat=err) rblck3
     if (err/=0) call myerr("reading rblk3 for pos",myname,crash) 

     do i = 1,ngas1
        psys%par(ngasread+i)%pos(1) = rblck3(1,i)
        psys%par(ngasread+i)%pos(2) = rblck3(2,i)
        psys%par(ngasread+i)%pos(3) = rblck3(3,i)
     end do
  
     deallocate(rblck3)
     
     ! read velocities 
     !-----------------------------------------------------------!  
     allocate(rblck3(3,ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck3 for vel",myname,crash)
     read(lun, iostat=err) rblck3
     if (err/=0) call myerr("reading rblk3 for vel",myname,crash) 

     do i = 1,ngas1
        psys%par(ngasread+i)%vel(1) = rblck3(1,i)
        psys%par(ngasread+i)%vel(2) = rblck3(2,i)
        psys%par(ngasread+i)%vel(3) = rblck3(3,i)
     end do
  
     deallocate(rblck3)

     ! read id's 
     !-----------------------------------------------------------!  
     allocate(iblck(ngas1), stat=err )
     if(err/=0) call myerr("allocating iblck for ID",myname,crash)
     read(lun, iostat=err) iblck  
     if (err/=0) call myerr("reading iblk for ID",myname,crash) 

     do i = 1,ngas1
        psys%par(ngasread+i)%id = iblck(i)
     end do

     deallocate(iblck)

     ! read masses 
     !-----------------------------------------------------------!  
     if (nmass1 > 0) then  ! if there are massive particles

        if (varmass(1)) then  ! if there are gas particles
           allocate(rblck(ngas1), stat=err)
           if(err/=0) call myerr("allocating rblck for mass",myname,crash)
           read(lun, iostat=err) rblck 
           if (err/=0) call myerr("reading rblk for mass",myname,crash) 

           do i = 1,ngas1
              psys%par(ngasread+i)%mass = rblck(i)
           end do

           deallocate(rblck)

        else 
           ! just dummy read non gas particles
           read(lun)
           psys%par(ngasread+i)%mass = ghead%mass(1)
        end if

     else
        psys%par(ngasread+1:ngasread+ngas1)%mass = ghead%mass(1)

     end if

     ! read temperature (internal energy until we change it)
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for u",myname,crash)
     read(lun, iostat=err) rblck
     if (err/=0) call myerr("reading rblk for u",myname,crash) 

     do i = 1,ngas1
        psys%par(ngasread+i)%T = rblck(i) 
     end do

     deallocate(rblck)


     ! read density 
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for rho",myname,crash)
     read(lun, iostat=err) rblck
     if (err/=0) call myerr("reading rblk for rho",myname,crash) 

     do i = 1, ngas1
        psys%par(ngasread+i)%rho = rblck(i)
     end do
     
     deallocate(rblck)
     
     ! read smoothing lengths 
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for hsml",myname,crash)
     read(lun, iostat=err) rblck
     if (err/=0) call myerr("reading rblk for hsml",myname,crash) 

     do i = 1, ngas1
        psys%par(ngasread+i)%hsml = rblck(i)
     end do
     
     deallocate(rblck)
     


     ngasread = ngasread + ngas1
     close(lun)

  end do files


  ! if not iso-temp convert internal energies / unit mass to temperature K
  !-------------------------------------------------------------------------
  if (GV%IsoTemp == 0.0) then
     do i = 1,ngas
        meanweight = 4.0d0 / (3.0d0 * GV%H_mf + 1.0d0 + 4.0d0)
        Tdum = meanweight * PROTONMASS / BOLTZMANN * (GAMMA - 1.0d0) * psys%par(i)%T
        psys%par(i)%T = Tdum * GENRG / GMASS
     end do
  else
     psys%par(:)%T = GV%IsoTemp
  end if



  ! set caseA true or false for collisional equilibrium
  !-----------------------------------------------------
  caseA = .false.
  if (.not. GV%OnTheSpotH  .or. GV%HydrogenCaseA) caseA(1) = .true.
  if (.not. GV%OnTheSpotHe .or. GV%HeliumCaseA)   caseA(2) = .true.

  
  ! if we have a single temperature
  !------------------------------------
  if (GV%IsoTemp /= 0.0) then
     call calc_colion_eq_fits(GV%IsoTemp, caseA, xvec)
     psys%par(:)%xHI = xvec(1)
     psys%par(:)%xHII = xvec(2)
#ifdef incHe
     psys%par(:)%xHeI = xvec(3)
     psys%par(:)%xHeII = xvec(4)
     psys%par(:)%xHeIII = xvec(5)
#endif

  ! if we have individual temperatures
  !------------------------------------
  else

     do i = 1,ngas
        Tdum = psys%par(i)%T
        call calc_colion_eq_fits(Tdum, caseA, xvec)
        psys%par(i)%xHI  = xvec(1)
        psys%par(i)%xHII = xvec(2)
#ifdef incHe
        psys%par(i)%xHeI   = xvec(3)
        psys%par(i)%xHeII  = xvec(4)
        psys%par(i)%xHeIII = xvec(5)
#endif
     end do

  end if

  psys%par(:)%ye = psys%par(:)%xHII 
#ifdef incHe
  nHe_over_nH = 0.25d0 * GV%He_mf / GV%H_mf
  psys%par(:)%ye = psys%par(:)%ye + ( psys%par(:)%xHeII + 2.0d0 * psys%par(:)%xHeIII ) * nHe_over_nH
#endif


end subroutine read_Gpub_particles



!> reads a Gadget-2 (Tiziana version) snapshot into a particle array  
!========================================================================
subroutine read_Gtiz_particles(GV,MB)
  type(global_variables_type), intent(in) :: GV !< global variables
  real(r8b), intent(out) :: MB         !< MB allocated for particle storage

  character(200) :: myname 
  logical :: crash

  real(r4b), allocatable :: rblck(:)
  real(r4b), allocatable :: rblck3(:,:)
  integer(i4b), allocatable :: iblck(:)
  integer(i8b) :: ngasread

  type(gadget_header_type) :: ghead
  character(200) :: snapfile
  integer(i8b) :: lun,i
  integer(i4b) :: err

  integer(i8b) :: npar, ngas, nmass
  integer(i8b) :: npar1, ngas1, nmass1
  integer(i8b) :: bytesperpar
  logical :: closefile
  logical :: header_verbose
  logical :: varmass(6)
  integer(i8b) :: fn

  real(r8b) :: meanweight
  logical :: caseA(2)
  real(r8b) :: xvec(5)
  real(r8b) :: Tdum


 
  ! set error handling variables
  !===============================
  myname = "read_Gtiz_particles"
  crash = .true.

  ! read header from first file
  !============================
  fn=0
  call form_Gsnapshot_file_name(GV%SnapPath,GV%ParFileBase,GV%CurSnapNum,fn,snapfile)
  closefile=.true.
  header_verbose=.false.
  call read_gadget_header(snapfile,header_verbose,closefile,ghead,lun)

  ! set local particle numbers
  !============================
  varmass = (ghead%npar_all > 0 .and. ghead%mass == 0)
  npar = sum(ghead%npar_all)
  ngas = ghead%npar_all(1)
  nmass = sum(ghead%npar_all, mask=varmass)

  ! do Gadget dummy checks
  !============================
  if (ngas .EQ. 0) call myerr("snapshot has no gas particles",myname,crash)


  ! calculate bytes per particle and allocate particle array
  !===========================================================
  call calc_bytes_per_particle(bytesperpar)
  MB = bytesperpar * real(ngas) / 2.**20

  101  format(A,F10.4,A,I10,A)
  write(*,101) "allocating ", MB, " MB for ", ngas, " particles"

  allocate (psys%par(ngas), stat=err)
  if (err /= 0) call myerr("failed to allocate par",myname,crash)


  ! now read all snapshot files
  !==============================          

  ngasread = 0
  closefile=.false.
  header_verbose=.false.
  files: do fn = 0, ghead%nfiles-1

     call form_Gsnapshot_file_name(GV%SnapPath,GV%ParFileBase,GV%CurSnapNum,fn,snapfile)
     write(*,'(A,A)') "reading tiziana gadget snapshot file ", trim(snapfile)
     call read_gadget_header(snapfile,header_verbose,closefile,ghead,lun)

     varmass = (ghead%npar_file > 0 .and. ghead%mass == 0)
     npar1 = sum(ghead%npar_file)
     ngas1 = ghead%npar_file(1)
     nmass1 = sum(ghead%npar_file, mask=varmass)

     if (ngas1 == 0) cycle

     ! read positions 
     !-----------------------------------------------------------!  
     allocate(rblck3(3,ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck3 for pos",myname,crash)
     read(lun, iostat=err) rblck3
     if (err/=0) call myerr("reading rblk3 for pos",myname,crash) 

     do i = 1,ngas1
        psys%par(ngasread+i)%pos(1) = rblck3(1,i)
        psys%par(ngasread+i)%pos(2) = rblck3(2,i)
        psys%par(ngasread+i)%pos(3) = rblck3(3,i)
     end do
  
     deallocate(rblck3)
     
     ! read velocities 
     !-----------------------------------------------------------!  
     allocate(rblck3(3,ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck3 for vel",myname,crash)
     read(lun, iostat=err) rblck3
     if (err/=0) call myerr("reading rblk3 for vel",myname,crash) 

     do i = 1,ngas1
        psys%par(ngasread+i)%vel(1) = rblck3(1,i)
        psys%par(ngasread+i)%vel(2) = rblck3(2,i)
        psys%par(ngasread+i)%vel(3) = rblck3(3,i)
     end do
  
     deallocate(rblck3)

     ! read id's 
     !-----------------------------------------------------------!  
     allocate(iblck(ngas1), stat=err )
     if(err/=0) call myerr("allocating iblck for ID",myname,crash)
     read(lun, iostat=err) iblck  
     if (err/=0) call myerr("reading iblk for ID",myname,crash) 

     do i = 1,ngas1
        psys%par(ngasread+i)%id = iblck(i)
     end do

     deallocate(iblck)

     ! read masses 
     !-----------------------------------------------------------!  
     if (nmass1 > 0) then  ! if there are massive particles

        if (varmass(1)) then  ! if there are gas particles
           allocate(rblck(ngas1), stat=err)
           if(err/=0) call myerr("allocating rblck for mass",myname,crash)
           read(lun, iostat=err) rblck 
           if (err/=0) call myerr("reading rblk for mass",myname,crash) 

           do i = 1,ngas1
              psys%par(ngasread+i)%mass = rblck(i)
           end do

           deallocate(rblck)

        else 
           ! just dummy read non gas particles
           read(lun)
           psys%par(ngasread+i)%mass = ghead%mass(1)
        end if

     else
        psys%par(ngasread+1:ngasread+ngas1)%mass = ghead%mass(1)

     end if

     ! read internal energy (internal energy until we change it)
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for u",myname,crash)
     read(lun, iostat=err) rblck
     if (err/=0) call myerr("reading rblk for u",myname,crash) 

     do i = 1,ngas1
        psys%par(ngasread+i)%T = rblck(i) 
     end do

     deallocate(rblck)

     ! read density 
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for rho",myname,crash)
     read(lun, iostat=err) rblck
     if (err/=0) call myerr("reading rblk for rho",myname,crash) 

     do i = 1, ngas1
        psys%par(ngasread+i)%rho = rblck(i)
     end do
     
     deallocate(rblck)

        
     ! read electron fraction 
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for ye",myname,crash)
     read(lun, iostat=err) rblck
     if (err/=0) call myerr("reading rblck for ye",myname,crash) 

     do i = 1, ngas1
        psys%par(ngasread+i)%ye = rblck(i)
     end do

     deallocate(rblck)

     ! read neutral hydrogen fraction (xHI)
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for xHI",myname,crash)
     read(lun, iostat=err) rblck
     if (err/=0) call myerr("reading rblk for xHI",myname,crash) 

     do i = 1, ngas1
        psys%par(ngasread+i)%xHI = rblck(i)
     end do
     
     deallocate(rblck)
     

     ! read smoothing lengths 
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for hsml",myname,crash)
     read(lun, iostat=err) rblck
     if (err/=0) call myerr("reading rblk for hsml",myname,crash) 

     do i = 1, ngas1
        psys%par(ngasread+i)%hsml = rblck(i)
     end do
     
     deallocate(rblck)
     

     ! the following data blocks are also in the snapshot file
     ! but we make no use of them for the time being

     ! star formation rate (ngas)
     ! stellar age (nstar)
     ! metallicity (ngas+nstar)
     ! black hole mass (nbh)
     ! black hole accretion rate (nbh)

     ngasread = ngasread + ngas1
     close(lun)

  end do files


  ! initialize xHII
  !-------------------------------------------------------------------------
  psys%par(:)%xHII = 1.0d0 - psys%par(:)%xHI

  ! if not iso-temp convert internal energies / unit mass to temperature K
  !-------------------------------------------------------------------------
  if (GV%IsoTemp == 0.0) then
     do i = 1,ngas
        meanweight = 4.0d0 / (3.0d0 * GV%H_mf + 1.0d0 + 4.0d0 * GV%H_mf * psys%par(i)%ye)
        Tdum = meanweight * PROTONMASS / BOLTZMANN * (GAMMA - 1.0d0) * psys%par(i)%T
        psys%par(i)%T = Tdum * GENRG / GMASS
     end do
  else
     psys%par(:)%T = GV%IsoTemp
  end if

  

  ! if we have Helium, initialize the ionization fractions to collisional equilibrium
  !----------------------------------------------------------------------------------
#ifdef incHe

  ! set caseA true or false for collisional equilibrium
  !-----------------------------------------------------
  caseA = .false.
  if (.not. GV%OnTheSpotH  .or. GV%HydrogenCaseA) caseA(1) = .true.
  if (.not. GV%OnTheSpotHe .or. GV%HeliumCaseA  ) caseA(2) = .true.

  ! if we have a single temperature
  !------------------------------------
  if (GV%IsoTemp /= 0.0) then

     call calc_colion_eq_fits(GV%IsoTemp, caseA, xvec)
     psys%par(:)%xHeI = xvec(3)
     psys%par(:)%xHeII = xvec(4)
     psys%par(:)%xHeIII = xvec(5)

  ! if we have individual temperatures
  !------------------------------------
  else

     do i = 1,ngas
        Tdum = psys%par(i)%T
        call calc_colion_eq_fits(Tdum, caseA, xvec)
        psys%par(i)%xHeI   = xvec(3)
        psys%par(i)%xHeII  = xvec(4)
        psys%par(i)%xHeIII = xvec(5)
     end do

  end if

#endif

 
end subroutine read_Gtiz_particles



!> Reads in an update snapshot.  The idea is to keep the ionization fractions
!! and temperature from the already loaded snapshot while updating the 
!! positions, velocities, smoothing lengths ... from the snapshot on disk.
!! The particles are reordered during the raytracing so when a new snapshot
!! is readin they must be matched
!===========================================================================
subroutine update_particles(GV,MB)
  type(global_variables_type), intent(in) :: GV !< global variables
  real(r8b), intent(out) :: MB   !< MB allocated for particle storage

  character(200) :: myname 
  logical :: crash

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


  ! set error handling variables
  !==============================
  myname = "update_particles"
  crash = .true.


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
#endif
  
  allocate( lasthitold(size(psys%par)), stat=err)
  if (err /= 0) call myerr("failed to allocate lasthitold",myname,crash)     
  lasthitold = psys%par%lasthit
  

  ! deallocate the current particle array and read new particle data
  !==================================================================
  deallocate( psys%par )
  if (GV%InputType == 2) then
     call read_Gpub_particles(GV,MB)
  else if (GV%InputType == 3) then
     call read_Gtiz_particles(GV,MB)
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
#ifdef incHe
  deallocate (xHeIold, xHeIIold, xHeIIIold)
#endif


end subroutine update_particles





!> forms a snapshot name from a path, a file base, a snapshot number
!> and a file number.
!===================================================================
subroutine form_Gsnapshot_file_name(path,base,SnapNum,FileNum,SnapFile)
  character(200), intent(in) :: path        !< path to snapshot dir
  character(200), intent(in) :: base        !< base snapshot name
  integer(i8b), intent(in) :: SnapNum       !< snapshot number
  integer(i8b), intent(in) :: FileNum       !< file number of snapshot
  character(200), intent(out) :: SnapFile   !< snapshot filename
  
  character(10) :: FileNumChar
  logical :: Fthere

  write(FileNumChar,"(I6)") FileNum
  100 format(A,"/",A,"_",I3.3)

  ! first write a file with no extension
  !--------------------------------------
  write(SnapFile,100) trim(path), trim(base), SnapNum
  inquire( file=SnapFile, exist=Fthere )

  ! if the file number is 0 and a file with no extension exists then return
  !-------------------------------------------------------------------------
  if (FileNum == 0 .and. Fthere) then
     return
  else
     SnapFile = trim(SnapFile) // "." // trim(adjustl(FileNumChar))
  end if

end subroutine form_Gsnapshot_file_name


!> forms a snapshot name from a path, a file base, a snapshot number
!! and a file number.
!===================================================================
subroutine form_snapshot_file_name(Path,FileBase,SnapNum,FileNum,SnapFile)
  character(200), intent(in) :: Path       !< path to snapshot dir
  character(200), intent(in) :: FileBase   !< file base names
  integer(i8b), intent(in) :: SnapNum      !< snapshot number
  integer(i8b), intent(in) :: FileNum      !< file number in snapshot
  character(200), intent(out) :: SnapFile  !< file name to return
  
  character(10) :: FileNumChar
  
  write(FileNumChar,"(I6)") FileNum
  
  100 format(A,"/",A,"_",I3.3,".",A)
  write(SnapFile,100) trim(Path), trim(FileBase), SnapNum, &
                      trim(adjustl(FileNumChar))
  
end subroutine form_snapshot_file_name


!> writes gadget header data to the screen
!-----------------------------------------
subroutine gadget_header_to_screen(ghead)
  type(gadget_header_type), intent(in) :: ghead !< particle header to print
  integer(i8b) :: i

98  format (A,T12,I12,T36,I12)
99  format (T17,A,T30,A,T40,A)
100 format (A,I1,A,A,A,I12,ES12.3,I12)
101 format (A,F8.5,A,F10.5,A,F8.4)  
102 format (A,F8.5,A,F8.5,A,F13.4)  
104 format (A,I3)
103 format (A,6(A,I1))

  write(*,*)
  write(*,99) "n_all","mass","n_file"
  do i = 1,6
     write(*,100) "type",i,"(",ptype_names(i),")",&
                   ghead%npar_all(i),ghead%mass(i),ghead%npar_file(i)
  end do
  write(*,98) "total:", sum(ghead%npar_all), sum(ghead%npar_file)

  write(*,101) "a   = ", ghead%a, &
               ", z    = ", ghead%z, &
               ", h       = ", ghead%h
  write(*,102) "OmegaM = ", ghead%OmegaM, &
               ", OmegaL = ", ghead%OmegaL, & 
               ", boxsize = ", ghead%boxlen 
  write(*,104) "nfiles  = ", ghead%nfiles
  write(*,103) "/flags/ ",&
             "  sfr=",ghead%flag_sfr, &
             ", feedback=",ghead%flag_feedback, &
             ", cooling=",ghead%flag_cooling, &
             ", age=",ghead%flag_age, &
             ", metals=",ghead%flag_metals, &
             ", entr_ics=", ghead%flag_entr_ics
  write(*,*) 

end subroutine gadget_header_to_screen


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

end module gadget_input_mod
