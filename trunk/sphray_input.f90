!> \file sphray_input.f90

!> \brief Handles readin of SPHRAY formatted snapshot files
!<

module sphray_input_mod
use myf90_mod 
use particle_system_mod, only: particle_type
use particle_system_mod, only: box_type
use particle_system_mod, only: calc_bytes_per_particle
use source_input_mod, only: source_header_type
use source_input_mod, only: read_source_header
use source_input_mod, only: source_header_to_file
use global_mod, only: global_variables_type
use global_mod, only: psys, PLAN, NgasInFile
implicit none



integer, parameter :: NdataBlocks = 30 !< number of possible data blks
integer, parameter :: TagLen = 24      !< # of characters in data tags

type sphray_header_type
   integer(i8b) :: NparSnap     !< number of particles in snapshot
   integer(i8b) :: NparFile     !< number of particles in file
   integer(i8b) :: Nfiles       !< number of files for this snapshot
   integer(i8b) :: BndryCond    !< boundary conditions employed
   integer(i8b) :: RaysTraced   !< number of rays traced in this snapshot

   real(r8b) :: Time         !< code time or scale factor
   real(r8b) :: BoxLower(3)  !< coordinates of lower box corner
   real(r8b) :: BoxUpper(3)  !< coordinates of upper box corner
   real(r8b) :: CGSlen       !< code length [cm/h] (Gadget Units)
   real(r8b) :: CGSmass      !< code mass [g/h] (Gadget Units)
   real(r8b) :: CGSvel       !< code velocity [cm/s] (Gadget Units)
   real(r8b) :: OmegaM       !< Omega Matter 
   real(r8b) :: OmegaB       !< Omega Baryon
   real(r8b) :: OmegaL       !< Omega Lambda
   real(r8b) :: LittleH      !< Hubble Parameter Now = 100 * LittleH km/s/Mpc

   integer(i4b) :: DataPresent(NdataBlocks)  !< 1 if data block present
end type sphray_header_type


type(sphray_header_type), allocatable :: saved_sheads(:,:)


!> names of the possible snapshot data fields
character(taglen), parameter :: data_tags(NdataBlocks) = &  
  (/ "id                      ", "pos                     ", &  ! 1,2
     "vel                     ", "hsml                    ", &  ! 3,4
     "rho                     ", "mass                    ", &  ! 5,6
     "temperature             ", "xHI                     ", &  ! 7,8
     "xHeI                    ", "xHeII                   ", &  ! 9,10
     "lasthit                 ", "gammaHI                 ", &  ! 11,12
     "undefined               ", "undefined               ", &  ! 13,14
     "undefined               ", "undefined               ", &  ! 15,16
     "undefined               ", "undefined               ", &  ! 17,18
     "undefined               ", "undefined               ", &  ! 19,20
     "undefined               ", "undefined               ", &  ! 21,22
     "undefined               ", "undefined               ", &  ! 23,24
     "undefined               ", "undefined               ", &  ! 25,26
     "undefined               ", "undefined               ", &  ! 27,28
     "undefined               ", "undefined               " /)  ! 29,30


contains




!>  reads in a particle header, prints the header to the screen if verbose 
!==========================================================================
subroutine read_sphray_header(snapfile,verbose,closefile,phead,lun)
  character(*), intent(in) :: snapfile   !< file containing particle header
  logical, intent(in) :: verbose         !< report to screen?
  logical, intent(in) :: closefile       !< close file when complete?
  type(sphray_header_type), intent(out) :: phead  !< particle header to read
  integer(i8b), intent(out) :: lun    !< output lun assigned to snapfile
   
  call open_unformatted_file_r(snapfile,lun)
  read(lun) phead
  if (closefile) close(lun)
  if (verbose) call sphray_header_to_screen(phead)

end subroutine read_sphray_header


!>  reads box data from a sphray snapshot 
!==========================================================================
subroutine read_sphray_box(GV,box)
  type(global_variables_type), intent(in) :: GV !< global variables
  type(box_type), intent(out) :: box    !< box 

  integer(i8b) :: fn
  character(200) :: snapfile   
  type(sphray_header_type) :: phead  
  integer(i8b) :: lun    

  fn=1
  call form_snapshot_file_name(GV%SnapPath,GV%ParFileBase,GV%CurSnapNum,fn,snapfile) 
  call open_unformatted_file_r(snapfile,lun)
  read(lun) phead
  close(lun) 

  box%bot = phead%BoxLower
  box%top = phead%BoxUpper
  box%bbound = GV%BndryCond
  box%tbound = GV%BndryCond

end subroutine read_sphray_box


!> puts code units read in from a particle snapshot header
!! into the global variables
!=========================================================
subroutine get_planning_data_sphray(GV)
use physical_constants_mod
use cosmology_mod, only: tsinceBB

  type(global_variables_type), intent(inout) :: GV !< global variables

  logical :: verbose
  logical :: closefile
  integer(i8b) :: lun
  type(sphray_header_type) :: phead
  type(source_header_type) :: shead

  integer(i8b) :: iSnap, fSnap    ! initial and final snapshot numbers
  integer(i8b) :: pfiles, sfiles  ! files/snap for particles and sources   
  integer(i8b) :: i,j        ! counters
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

  allocate( saved_sheads(iSnap:fSnap, 1:pfiles) )
  allocate( NgasInFile(iSnap:fSnap, 1:pfiles) )

  ! read in and report to screen
  !===============================
  closefile = .true.
  verbose = .false.

  ! read all particle headers and write to log file
  !===================================================
  write(loglun,'(A)') "reading all SPHRAY particle header(s) ... "
  do i = iSnap,fSnap
     do j = 1,pfiles

        call form_snapshot_file_name(GV%SnapPath,GV%ParFileBase,i,j,snapfile)
        write(loglun,'(I3,"  ",A)') i,trim(snapfile)
        call read_sphray_header(snapfile,verbose,closefile,phead,lun)
        call sphray_header_to_file(phead,loglun)
        
        saved_sheads(i,j) = phead
        NgasInFile(i,j) = phead%NparFile

        ! IsoTemp comes from config file

        GV%BoxLowers(1:3) = phead%BoxLower(1:3)
        GV%BoxUppers(1:3) = phead%BoxUpper(1:3)

        GV%OmegaM   = phead%OmegaM
        GV%OmegaB   = phead%OmegaB
        GV%OmegaL   = phead%OmegaL
        GV%LittleH  = phead%LittleH

        GV%cgs_len  = phead%CGSlen
        GV%cgs_mass = phead%CGSmass
        GV%cgs_vel  = phead%CGSvel

        GV%cgs_time = GV%cgs_len  / GV%cgs_vel
        GV%cgs_rho  = GV%cgs_mass / GV%cgs_len**3 
        GV%cgs_prs  = GV%cgs_mass / GV%cgs_len / GV%cgs_time**2
        GV%cgs_enrg = GV%cgs_mass * GV%cgs_vel**2
        GV%cgs_lum  = GV%cgs_enrg / GV%cgs_time

        if (GV%Comoving) then
           PLAN(i)%ScalefacAtSnap = phead%Time
           PLAN(i)%TimeAtSnap = tsinceBB(phead%Time, GV%OmegaM, GV%LittleH) 
           PLAN(i)%TimeAtSnap = PLAN(i)%TimeAtSnap * GV%LittleH / GV%cgs_time
        else
           PLAN(i)%ScalefacAtSnap = 1.0d0
           PLAN(i)%TimeAtSnap = phead%Time
        end if
        
     end do
  end do
  

  ! read all source headers and write to log file
  !===================================================
  write(loglun,*) 
  write(loglun,'(A)') "reading all Sphray source header(s) ... "
  do i = iSnap,fSnap
     do j = 1,sfiles
        call form_snapshot_file_name(GV%SourcePath,GV%SourceFileBase,i,j,snapfile)
        write(loglun,'(I3,"  ",A)') i,trim(snapfile)
        call read_source_header(snapfile,verbose,closefile,shead,lun)
        call source_header_to_file(shead,loglun)

        PLAN(i)%RaysFromSrcHeader = shead%TotalRays
        GV%Lunit = shead%Lunit
     end do
  end do

  close(loglun)


  ! write units to log file
  !===================================================

  logfile = trim(GV%OutputDir) // "/" // "code_units.log"
  call open_formatted_file_w(logfile,loglun)


  106 format(A,ES12.5,A)
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


 end subroutine get_planning_data_sphray


!> Reads in a SPHRAY particle snapshot.  Required data fields are
!! id, pos, vel, hsml, rho, mass, temp, xHI, lasthit.  The rest are optional.  
!===========================================================================
subroutine read_sphray_particles(GV,verbose,MB)
  type(global_variables_type), intent(in) :: GV   !< global variables
  logical :: verbose                    !< report to screen?
  real(r8b), intent(out) :: MB   !< MB allocated for particle storage

  character(200) :: myname 
  logical :: crash

  real(r4b), allocatable :: rblk(:)
  real(r4b), allocatable :: rblk3(:,:)
  integer(i8b), allocatable :: iblk(:)
  integer(i8b) :: readin

  type(sphray_header_type) :: phead
  character(200) :: snapfile
  integer(i8b) :: lun
  integer(i4b) :: err
  integer(i8b) :: Nfile, Nsnap

  integer(i8b) :: bytesperpar
  integer(i8b) :: fn
  integer(i8b) :: p 
  logical :: closefile
  logical :: header_verbose
  integer(i8b) :: ii,ff


  ! set error handling variables
  !===============================  
  myname = "read_sphray_particles"
  crash = .true.

  ! read header from first file
  !============================
  fn=1
  call form_snapshot_file_name(GV%SnapPath,GV%ParFileBase,GV%CurSnapNum,fn,snapfile)
  closefile=.true.
  header_verbose=.false.
  call read_sphray_header(snapfile,header_verbose,closefile,phead,lun)
 
  ! set local particle number and allocate particle array
  !=======================================================
  Nsnap = phead%NparSnap  
  call calc_bytes_per_particle(bytesperpar)
  MB = bytesperpar * real(Nsnap) / 2.**20

  101  format(A,F10.4,A,I10,A) 
  write(*,101) "allocating ", MB, " MB for ", Nsnap, " particles"

  allocate (psys%par(Nsnap),stat=err)
  if (err /= 0) call myerr("failed to allocate par",myname,crash)

  ! now read all snapshot files
  !==============================  

  readin = 0
  closefile=.false.
  header_verbose=.false.
  do fn = 1, phead%Nfiles

     call form_snapshot_file_name(GV%SnapPath,GV%ParFileBase,GV%CurSnapNum,fn,snapfile)

     if (verbose) then 
        write(*,'(A,A)') "reading sphray particle snapshot file: ", trim(snapfile)
     end if

     call read_sphray_header(snapfile,header_verbose,closefile,phead,lun)
     
     Nfile = phead%NparFile  
     ii = readin + 1
     ff = readin + Nfile

     ! read ID's (REQUIRED)
     !-----------------------------------------------------------!  
     if ( phead%DataPresent(1) == 0 ) then
        call myerr("snapshots do not contain ID data",myname,crash)
     else
        allocate( iblk(Nfile), stat=err )
        if (err/=0) call myerr("allocating iblk for IDs",myname,crash)
        read( lun, iostat=err) iblk
        if (err/=0) call myerr("reading iblk for IDs",myname,crash)
        psys%par(ii:ff)%id = iblk
        deallocate( iblk )
     end if

     ! read positions (REQUIRED)
     !-----------------------------------------------------------!  
     if ( phead%DataPresent(2) == 0 ) then
        call myerr("snapshots do not contain pos data",myname,crash)
     else
        allocate( rblk3(3,Nfile), stat=err )
        if (err/=0) call myerr("allocating rblk3 for pos",myname,crash)
        read( lun, iostat=err) rblk3
        if (err/=0) call myerr("reading rblk3 for pos",myname,crash)
        do p = 1,Nfile
           psys%par(readin+p)%pos(1) = rblk3(1,p)
           psys%par(readin+p)%pos(2) = rblk3(2,p)
           psys%par(readin+p)%pos(3) = rblk3(3,p)
        end do
        deallocate( rblk3 )
     end if

     ! read velocities (REQUIRED)
     !-----------------------------------------------------------!  
     if ( phead%DataPresent(3) == 0 ) then
        call myerr("snapshots do not contain vel data",myname,crash)
     else
        allocate( rblk3(3,Nfile), stat=err )
        if (err/=0) call myerr("allocating rblk3 for vel",myname,crash)
        read( lun, iostat=err) rblk3
        if (err/=0) call myerr("reading rblk3 for vel",myname,crash)
        do p = 1,Nfile
           psys%par(readin+p)%vel(1) = rblk3(1,p)
           psys%par(readin+p)%vel(2) = rblk3(2,p)
           psys%par(readin+p)%vel(3) = rblk3(3,p)
        end do
        deallocate( rblk3 )
     end if

     ! read smoothing lengths (REQUIRED)
     !-----------------------------------------------------------!  
     if ( phead%DataPresent(4) == 0 ) then
        call myerr("snapshots do not contain hsml data",myname,crash)
     else
        allocate( rblk(Nfile), stat=err )
        if (err/=0) call myerr("allocating rblk for hsml",myname,crash)
        read( lun, iostat=err) rblk
        if (err/=0) call myerr("reading rblk for hsml",myname,crash)
        psys%par(ii:ff)%hsml = rblk
        deallocate( rblk )
     end if

     ! read density (REQUIRED)
     !-----------------------------------------------------------!  
     if ( phead%DataPresent(5) == 0 ) then
        call myerr("snapshots do not contain rho data",myname,crash)
     else
        allocate( rblk(Nfile), stat=err )
        if (err/=0) call myerr("allocating rblk for rho",myname,crash)
        read( lun, iostat=err) rblk
        if (err/=0) call myerr("reading rblk for rho",myname,crash)
        psys%par(ii:ff)%rho = rblk
        deallocate( rblk )
     end if

     ! read masses (REQUIRED)  
     !-----------------------------------------------------------!  
     if ( phead%DataPresent(6) == 0 ) then
        call myerr("snapshots do not contain mass data",myname,crash)
     else
        allocate( rblk(Nfile), stat=err )
        if (err/=0) call myerr("allocating rblk for mass",myname,crash)
        read( lun, iostat=err) rblk
        if (err/=0) call myerr("reading rblk for mass",myname,crash)
        psys%par(ii:ff)%mass = rblk
        deallocate( rblk )
     end if

     ! read temperature (REQUIRED)
     !-----------------------------------------------------------!  
     if ( phead%DataPresent(7) == 0 ) then
        call myerr("snapshots do not contain temperature data",myname,crash)
     else
        allocate( rblk(Nfile), stat=err )
        if (err/=0) call myerr("allocating rblk for temperature",myname,crash)
        read( lun, iostat=err) rblk
        if (err/=0) call myerr("reading rblk for temperature",myname,crash)
        psys%par(ii:ff)%T = rblk
        deallocate( rblk )
     end if

     ! read xHI (REQUIRED)  
     !-----------------------------------------------------------!  
     if ( phead%DataPresent(8) == 0 ) then
        call myerr("snapshots do not contain xHI data",myname,crash)
     else
        allocate( rblk(Nfile), stat=err )
        if (err/=0) call myerr("allocating rblk for xHI",myname,crash)
        read( lun, iostat=err) rblk
        if (err/=0) call myerr("reading rblk for xHI",myname,crash)
        psys%par(ii:ff)%xHI = rblk
        deallocate( rblk )
     end if

     ! read xHeI (OPTIONAL)
     !-----------------------------------------------------------!  
     if ( phead%DataPresent(9) == 1 .and. GV%datamask(9) == 0 ) then
        read(lun) 
     else if ( phead%DataPresent(9) == 0 .and. GV%datamask(9) == 1 ) then
        call myerr("xHeI is unmasked, but not present",myname,crash)
     else if ( phead%DataPresent(9) == 1 .and. GV%datamask(9) == 1 ) then
        allocate( rblk(Nfile), stat=err )
        if (err/=0) call myerr("allocating rblk for xHeI",myname,crash)
        read( lun, iostat=err) rblk
        if (err/=0) call myerr("reading rblk for xHeI",myname,crash)
#ifdef incHe
        psys%par(ii:ff)%xHeI = rblk
#endif
        deallocate( rblk )
     end if

     ! read xHeII (OPTIONAL)
     !-----------------------------------------------------------!  
     if ( phead%DataPresent(10) == 1 .and. GV%datamask(10) == 0 ) then
        read(lun) 
     else if ( phead%DataPresent(10) == 0 .and. GV%datamask(10) == 1 ) then
        call myerr("xHeII is unmasked, but not present",myname,crash)
     else if ( phead%DataPresent(10) == 1 .and. GV%datamask(10) == 1 ) then
        allocate( rblk(Nfile), stat=err )
        if (err/=0) call myerr("allocating rblk for xHeII",myname,crash)
        read( lun, iostat=err) rblk
        if (err/=0) call myerr("reading rblk for xHeII",myname,crash)
#ifdef incHe
        psys%par(ii:ff)%xHeII = rblk
#endif
        deallocate( rblk )
     end if


     ! read last hit (REQUIRED)
     !-----------------------------------------------------------!  
     if ( phead%DataPresent(11) == 0 ) then
        call myerr("snapshots do not contain lasthit data",myname,crash)
     else
        allocate( iblk(Nfile), stat=err )
        if (err/=0) call myerr("allocating iblk for lasthit",myname,crash)
        read( lun, iostat=err) iblk
        if (err/=0) call myerr("reading iblk for lasthit",myname,crash)
        psys%par(ii:ff)%lasthit = iblk
        deallocate( iblk )
     end if

     readin = readin + Nfile
     close(lun)     

     psys%par(:)%xHII = 1.0d0 - psys%par(:)%xHI

#ifdef incHe
     psys%par(:)%xHeIII = 1.0d0 - psys%par(:)%xHeI - psys%par(:)%xHeII
#endif

  end do


end subroutine read_sphray_particles



!> Reads in an update snapshot.  The idea is to keep the ionization fractions
!! and temperature from the already loaded snapshot while updating the 
!! positions, velocities, smoothing lengths ... from the snapshot on disk.
!! The particles are reordered during the raytracing so when a new snapshot
!! is readin they must be matched
!===========================================================================
subroutine update_sphray_particles(GV,verbose,MB)
use m_mrgrnk, only: mrgrnk

  type(global_variables_type), intent(in) :: GV   !< global variables
  logical :: verbose !< report to screen?
  real(r8b), intent(out) :: MB   !< MB allocated for particle storage

  character(200) :: myname 
  logical :: crash

  integer(i8b), allocatable :: idold(:)
  real(r4b), allocatable :: Told(:)
  real(r4b), allocatable :: xHIold(:)
  real(r4b), allocatable :: xHeIold(:)
  real(r4b), allocatable :: xHeIIold(:)
  integer(i8b), allocatable :: lasthitold(:)
  integer(i8b), allocatable :: orderold(:)  

  integer(i8b) :: minIDold, minIDnew
  integer(i8b) :: maxIDold, maxIDnew

  integer(i8b) :: i, idnew, indx, err
  

  ! set error handling variables
  !==============================
  myname = "update_sphray_particles"
  crash = .true.


  ! store carry over variables
  !==============================

  allocate( idold(size(psys%par)), stat=err )
  if (err /= 0) call myerr("failed to allocate IDold",myname,crash)     
  idold = psys%par%id
  minIDold = minval(idold)
  maxIDold = maxval(idold)
  
  allocate( Told(size(psys%par)), stat=err )
  if (err /= 0) call myerr("failed to allocate Told",myname,crash)     
  Told = psys%par%T
  
  allocate( xHIold(size(psys%par)), stat=err )
  if (err /= 0) call myerr("failed to allocate xHIold",myname,crash)     
  xHIold = psys%par%xHI
  
#ifdef incHe
  allocate( xHeIold(size(psys%par)), stat=err )
  if (err /= 0) call myerr("failed to allocate xHeIold",myname,crash)     
  xHeIold = psys%par%xHeI

  allocate( xHeIIold(size(psys%par)), stat=err )
  if (err /= 0) call myerr("failed to allocate xHeIIold",myname,crash)     
  xHeIIold = psys%par%xHeII
#endif

  allocate( lasthitold(size(psys%par)), stat=err )
  if (err /= 0) call myerr("failed to allocate lasthitold",myname,crash)     
  lasthitold = psys%par%lasthit
  

  ! deallocate the current particle array and read new particle data
  !==================================================================
  deallocate( psys%par )
  call read_sphray_particles(GV,verbose,MB)
  minIDnew = minval(psys%par%id)
  maxIDnew = maxval(psys%par%id)

  ! match id's from new snapshot data to the old IDs
  !===================================================
 
  
  ! create an indexing array for the old particle array ID's
  ! ID's that are no longer present will be indexed with minIDold - 1
  !--------------------------------------------------------------------
  write(*,*) "creating indexing array for current particle system..."
  write(*,*) "min/max ID (old snap) = ", minIDold, maxIDold

  allocate( orderold(minIDold:maxIDold), stat=err )
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


     psys%par(i)%xHI = xHIold(indx) 
#ifdef incHe
     psys%par(i)%xHeI  = xHeIold(indx)
     psys%par(i)%xHeII = xHeIIold(indx)
#endif

     psys%par(i)%lasthit = lasthitold(indx)

  end do

  deallocate (idold, Told, xHIold, lasthitold, orderold)
#ifdef incHe
  deallocate (xHeIold, xHeIIold)
#endif

end subroutine update_sphray_particles




!> outputs particle snapshot header to the screen
!===================================================
subroutine sphray_header_to_screen(phead)
  type(sphray_header_type), intent(in) :: phead !< particle header to print
  integer(i8b) :: i

  97  format(T2,"=",T58,"=")
  98  format(T2,"= ",A,T58,"=")
  99  format(T2,57("="))
  201 format(T2,"=",T4,I3,T8,A,T58,"=")
  205 format(T2,"=",T4,A,T20,I20,T58,"=")
  210 format(T2,"=",T4,A,T25,ES20.8,T58,"=")
  215 format(T2,"=",T4,A,T15,3ES14.4,T58,"=")

  write(*,99)
  write(*,98) "sphray header data"
  write(*,99)
  write(*,97) 
  
  write(*,205) "NparSnap:    ", phead%NparSnap
  write(*,205) "NparFile:    ", phead%NparFile
  write(*,205) "Nfiles:      ", phead%Nfiles
  write(*,205) "BndryCond:   ", phead%BndryCond
  write(*,205) "RaysTraced:  ", phead%RaysTraced
  write(*,97) 
  write(*,210) "Time:        ", phead%Time
  write(*,97) 
  write(*,215) "BoxLower:    ", phead%BoxLower
  write(*,215) "BoxUpper:    ", phead%BoxUpper
  write(*,97) 
  write(*,210) "cgs length[cm/h]:   ", phead%CGSlen
  write(*,210) "cgs mass[g/h]:      ", phead%CGSmass
  write(*,210) "cgs velocity[cm/s]: ", phead%CGSvel
  write(*,210) "Omega Matter:     ", phead%OmegaM
  write(*,210) "Omega Baryon:     ", phead%OmegaB
  write(*,210) "Omega Lambda:     ", phead%OmegaL
  write(*,210) "Hubble Parameter: ", phead%LittleH 
  
  write(*,97) 
  write(*,98) "The snapshot contains the following data blocks"
  do i = 1,NdataBlocks
     if (phead%DataPresent(i)==1) write(*,201) i, data_tags(i)
  end do
  write(*,97)
  
  write(*,99)  
  write(*,*) 
  

end subroutine sphray_header_to_screen


!> outputs particle snapshot header to the screen
!===================================================
subroutine sphray_header_to_file(phead,lun)
  type(sphray_header_type), intent(in) :: phead !< particle header to print
  integer(i8b), intent(in) :: lun

  integer(i8b) :: i

  97  format(T2,"=",T58,"=")
  98  format(T2,"= ",A,T58,"=")
  99  format(T2,57("="))
  201 format(T2,"=",T4,I3,T8,A,T58,"=")
  205 format(T2,"=",T4,A,T20,I20,T58,"=")
  210 format(T2,"=",T4,A,T25,ES20.8,T58,"=")
  215 format(T2,"=",T4,A,T15,3ES14.4,T58,"=")

  write(lun,99)
  write(lun,98) "sphray header data"
  write(lun,99)
  write(lun,97) 
  
  write(lun,205) "NparSnap:    ", phead%NparSnap
  write(lun,205) "NparFile:    ", phead%NparFile
  write(lun,205) "Nfiles:      ", phead%Nfiles
  write(lun,205) "BndryCond:   ", phead%BndryCond
  write(lun,205) "RaysTraced:  ", phead%RaysTraced
  write(lun,97) 
  write(lun,210) "Time:        ", phead%Time
  write(lun,97) 
  write(lun,215) "BoxLower:    ", phead%BoxLower
  write(lun,215) "BoxUpper:    ", phead%BoxUpper
  write(lun,97) 
  write(lun,210) "cgs length[cm/h]:   ", phead%CGSlen
  write(lun,210) "cgs mass[g/h]:      ", phead%CGSmass
  write(lun,210) "cgs velocity[cm/s]: ", phead%CGSvel
  write(lun,210) "Omega Matter:     ", phead%OmegaM
  write(lun,210) "Omega Baryon:     ", phead%OmegaB
  write(lun,210) "Omega Lambda:     ", phead%OmegaL
  write(lun,210) "Hubble Parameter: ", phead%LittleH 
  
  write(lun,97) 
  write(lun,98) "The snapshot contains the following data blocks"
  do i = 1,NdataBlocks
     if (phead%DataPresent(i)==1) write(lun,201) i, data_tags(i)
  end do
  write(lun,97)
  
  write(lun,99)  
  write(lun,*) 
  

end subroutine sphray_header_to_file


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




end module sphray_input_mod
