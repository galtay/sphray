!> \file sphray_input.f90

!> \brief the module that handles readin of SPHRAY formatted files
!<

module sphray_input_mod
use myf90_mod 
use particle_system_mod
use global_mod, only: global_variables_type
implicit none

private

integer, parameter :: NdataBlocks = 30 !< number of possible data blks
integer, parameter :: TagLen = 24      !< # of characters in data tags

type sphray_header_type
   integer(i8b) :: NparSnap     !< number of particles in snapshot
   integer(i8b) :: NparFile     !< number of particles in file
   integer(i8b) :: Nfiles       !< number of files for this snapshot
   integer(i8b) :: NsphNnb      !< number of smoothing neighbors
   integer(i8b) :: BndryCond    !< boundary conditions employed
   integer(i8b) :: OnTheSpot    !< on the spot approximation?
   integer(i8b) :: RaysTraced   !< number of rays traced through this snapshot

   real(r8b) :: Time         !< code time or scale factor
   real(r8b) :: IsoMass      !< if > 0 then unique particle mass
   real(r8b) :: IsoTemp      !< if > 0 then unique particle temperature
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


!> names of the possible snapshot data fields
character(taglen), parameter :: data_tags(NdataBlocks) = &  
  (/ "id                      ", "pos                     ", &  ! 1,2
     "vel                     ", "hsml                    ", &  ! 3,4
     "rho                     ", "mass                    ", &  ! 5,6
     "temperature             ", "xHII                    ", &  ! 7,8
     "xHeII                   ", "xHeIII                  ", &  ! 9,10
     "xHIIrc                  ", "xHeIIrc                 ", &  ! 11,12
     "xHeIIIrc                ", "lasthit                 ", &  ! 13,14
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
subroutine read_particle_header(snapfile,verbose,closefile,phead,lun)
  character(*), intent(in) :: snapfile   !< file containing particle header
  logical, intent(in) :: verbose         !< report to screen?
  logical, intent(in) :: closefile       !< close file when complete?
  type(particle_header_type), intent(out) :: phead  !< particle header to read
  integer, intent(out) :: lun    !< output lun assigned to snapfile

  integer :: i
  
  call open_unformatted_file_r(snapfile,lun)
  read(lun) phead
  if (closefile) close(lun)
  if (verbose) call particle_header_to_screen(phead)

end subroutine read_particle_header




!> routine for reading in an integer data block from the particle snapshot
!=========================================================================
subroutine read_particle_iblock(lun,Nbl,N,iblk)
  integer, intent(in) :: lun              !< logical unit number to read
  integer, intent(in) :: Nbl              !< data block num (1=id, 2=xpos, ...)
  integer, intent(in) :: N                !< size of data block
  integer(i8b), intent(inout) :: iblk(N)  !< integer data block
  
  character(taglen) :: tag
  integer :: ioerr
  
  read(lun,iostat=ioerr) tag,iblk
  if (ioerr/=0) then
     write(*,*) "iostat error reading iblock ", trim(tag)
     stop
  end if
  
  if (tag /= Dtags(Nbl) ) then
     write(*,*) "  *** error - corrupt snapshot ***  "
     write(*,*) "expected data label: ", Dtags(Nbl)
     write(*,*) "read data label:     ", tag
     stop
  end if
  
end subroutine read_particle_iblock

!> routine for reading in a real data block from the particle snapshot
!=====================================================================
subroutine read_particle_rblock(lun,Nbl,N,rblk)
  integer, intent(in) :: lun          !< logical unit number to read
  integer, intent(in) :: Nbl          !< data block num (1=id, 2=xpos, ...)
  integer, intent(in) :: N            !< size of data block
  real(r4b), intent(inout) :: rblk(N) !< real data block
  
  character(taglen) :: tag
  integer :: ioerr
  
  read(lun,iostat=ioerr) tag,rblk
  if (ioerr/=0) then
     write(*,*) "iostat error reading rblock ", trim(tag)
     stop
  end if
  
  if (tag /= Dtags(Nbl) ) then
     write(*,*) "  *** error - corrupt snapshot ***  "
     write(*,*) "expected data label: ", Dtags(Nbl)
     write(*,*) "read data label:     ", tag
     stop
  end if
  
end subroutine read_particle_rblock


!> Makefile wants data not in snapshot - stop
!----------------------------------------------
subroutine make_want_data(str1,str2)
  character(*) :: str1  !< data block snapshots dont contain
  character(*) :: str2  !< Makefile flag leading SPRHAY to think they do
  
  write(stderr,*) "********************** ERROR **********************"
  write(stderr,*) "these snapshots dont include ", str1, " information"
  write(stderr,*) "please recompile SPHRAY without the Makefile flag(s) ", str2
  write(stderr,*) "********************** ERROR **********************"
  stop
  
end subroutine make_want_data

!> snapshot contains data not requested by Makefile - warning
!--------------------------------------------------------------
subroutine snap_has_xtra_data(str1,str2)
  character(*) :: str1  !< Makefile flag that is commented out
  character(*) :: str2  !< data block that could potentially be used

  write(stderr,*) "---------------------------------------------"
  write(stderr,*) "** WARNING ** ", str2, " flag commented out, "
  write(stderr,*) "but snapshots contain ", str1, " information"
  write(stderr,*) "you could use this data if you wanted"
  write(stderr,*) "---------------------------------------------"    

end subroutine snap_has_xtra_data

!> snapshot is missing a required data block - stop
!---------------------------------------------------------------------
subroutine need_datablock(Nbl)
  integer, intent(in) :: Nbl  !< necessary block number that is not present
  
  write(stderr,*) "********************** ERROR **********************"
  write(stderr,*) Dtags(Nbl), " not in snapshot, cannot proceed."
  write(stderr,*) "********************** ERROR **********************"
  stop
  
end subroutine need_datablock


!> reads in the initial particle snapshot.  The required data fields
!! for an initial snapshot are pos, hsml, rho, xHII, lasthit.  The rest are
!! optional.  
!===========================================================================
subroutine read_initial_par_snapshot(path,base,Ns,verbose,psys,MB)
  character(200), intent(in) :: path   !< path to dir containing snapshots
  character(200), intent(in) :: base   !< base name of snapshots
  integer, intent(in) :: Ns            !< number of snapshot to read
  logical :: verbose                   !< report to screen?
  type(particle_system_type), intent(inout) :: psys  !< particle system to load
  real, intent(out) :: MB   !< MB allocated for particle storage

  real(r4b), allocatable :: rblk(:)
  integer(i8b), allocatable :: iblk(:)

  character(200) :: snapfile
  type(particle_header_type) :: phead
  integer :: lun,err
  integer :: N, Nall

  logical :: DataInSnap(NdataBlocks)
  integer :: bytesperpar, Nbl
  logical :: closefile
  logical :: header_verbose
  integer :: fn

  fn=1
  call form_snapshot_file_name(path,base,Ns,fn,snapfile)
  write(*,'(A,A)') "reading initial particle snapshot from ", trim(snapfile)

  closefile=.false.
  header_verbose=.false.
  call read_particle_header(snapfile,header_verbose,closefile,phead,lun)
  DataInSnap = .false.
  where (phead%DataPresent==1) DataInSnap=.true.

  psys%box%bot(1:3) = phead%real(5:7)%var  
  psys%box%top(1:3) = phead%real(8:10)%var 
  ! bndry conditions are set in main_input.f90

  ! do particle numbers
  Nall = phead%int(1)%var  ! NparSnap
  N = phead%int(2)%var     ! NparFile

  psys%npar = Nall
  write(*,*) "psys%npar = ", psys%npar
  call calc_bytes_per_particle(bytesperpar)
  MB = bytesperpar * real(Nall) / 2**20

  100  format(A,I3)
  101  format(A,F10.4,A,I10,A)
  write(*,100) "bytes per particle = ", bytesperpar
  write(*,101) "allocating ", MB, " MB for ", Nall, " particles"

  ! allocations
  allocate (psys%par(Nall),stat=err)
  if (err /= 0) stop "cant allocate mem for particles in sphray_input.f90"

  allocate (rblk(N), stat=err)
  if (err/=0) stop "cant allocate rblock in sphray_input.f90"
  allocate(iblk(N), stat=err)
  if(err/=0) stop "cant allocate iblock in sphray_input.f90"

! read id's (OPTIONAL)
!-----------------------------------------------------------!  
  Nbl=1
  if (DataInSnap(Nbl)) then
     call read_particle_iblock(lun,Nbl,N,iblk)
#ifndef incid 
     call snap_has_xtra_data("id","incid")
#endif
  end if
   
#ifdef incid
  if (.not. DataInSnap(Nbl)) call make_want_data("id","incid")
  psys%par(1:N)%id = iblk(1:N)  
  if (verbose) write(*,*)"read ", Dtags(Nbl)
#endif

! read positions (REQUIRED)
!-----------------------------------------------------------!  
  Nbl=2
  if (.not. DataInSnap(Nbl)) call need_datablock(Nbl)
  call read_particle_rblock(lun,Nbl,N,rblk)
  psys%par(1:N)%pos(1) = rblk(1:N) 
  if (verbose) write(*,*)"read ", Dtags(Nbl)
  
  Nbl=3
  if (.not. DataInSnap(Nbl)) call need_datablock(Nbl)
  call read_particle_rblock(lun,Nbl,N,rblk)
  psys%par(1:N)%pos(2) = rblk(1:N) 
  if (verbose) write(*,*)"read ", Dtags(Nbl)

  Nbl=4
  if (.not. DataInSnap(Nbl)) call need_datablock(Nbl)
  call read_particle_rblock(lun,Nbl,N,rblk)
  psys%par(1:N)%pos(3) = rblk(1:N) 
  if (verbose) write(*,*)"read ", Dtags(Nbl)

! read velocities (OPTIONAL)
!-----------------------------------------------------------!  
  Nbl=5
  if (DataInSnap(Nbl)) then
     call read_particle_rblock(lun,Nbl,N,rblk)
#ifndef incvel 
     call snap_has_xtra_data("vel","incvel")  
#endif
  end if

#ifdef incvel
  if (.not. DataInSnap(Nbl)) call make_want_data("x-vel","incvel")
  psys%par(1:N)%vel(1) = rblk(1:N) 
  if (verbose) write(*,*)"read ", Dtags(Nbl)
#endif        
        
  Nbl=6
  if (DataInSnap(Nbl)) call read_particle_rblock(lun,Nbl,N,rblk)

#ifdef incvel
  if (.not. DataInSnap(Nbl)) call make_want_data("y-vel","incvel")
  psys%par(1:N)%vel(2) = rblk(1:N) 
  if (verbose) write(*,*)"read ", Dtags(Nbl)
#endif        

  Nbl=7
  if (DataInSnap(Nbl)) call read_particle_rblock(lun,Nbl,N,rblk)

#ifdef incvel
  if (.not. DataInSnap(Nbl)) call make_want_data("z-vel","incvel")
  psys%par(1:N)%vel(3) = rblk(1:N) 
  if (verbose) write(*,*)"read ", Dtags(Nbl)
#endif        

! read smoothing lengths (REQUIRED)
!-----------------------------------------------------------!  
  Nbl=8
  if (.not. DataInSnap(Nbl)) call need_datablock(Nbl)
  call read_particle_rblock(lun,Nbl,N,rblk)
  psys%par(1:N)%hsml = rblk(1:N) 
  if (verbose) write(*,*)"read ", Dtags(Nbl)

! read density (REQUIRED)
!-----------------------------------------------------------!  
  Nbl=9
  if (.not. DataInSnap(Nbl)) call need_datablock(Nbl)
  call read_particle_rblock(lun,Nbl,N,rblk)
  psys%par(1:N)%rho = rblk(1:N) 
  if (verbose) write(*,*)"read ", Dtags(Nbl)

! read masses (OPTIONAL)  
!-----------------------------------------------------------!  
  Nbl=10
  if (DataInSnap(Nbl)) then
     call read_particle_rblock(lun,Nbl,N,rblk)
#ifndef incmass 
     call snap_has_xtra_data("mass","incmass")  
#endif
  end if

#ifdef incmass
  if (.not. DataInSnap(Nbl)) call make_want_data("mass","incmass")
  psys%par(1:N)%mass = rblk(1:N) 
  if (verbose) write(*,*)"read ", Dtags(Nbl)
#endif        

! read temperature (OPTIONAL)
!-----------------------------------------------------------!  
  Nbl=11
  if (DataInSnap(Nbl)) then
     call read_particle_rblock(lun,Nbl,N,rblk)
#ifndef incT 
     call snap_has_xtra_data("temp","incT")  
#endif
  end if

#ifdef incT
  if (.not. DataInSnap(Nbl)) call make_want_data("temp","incT")
  psys%par(1:N)%T = rblk(1:N) 
  if (verbose) write(*,*)"read ", Dtags(Nbl)
#endif        

! read xHII (REQUIRED)  
!-----------------------------------------------------------!  
  Nbl=12
  if (.not. DataInSnap(Nbl)) call need_datablock(Nbl)
  call read_particle_rblock(lun,Nbl,N,rblk)
  psys%par(1:N)%xHII = rblk(1:N) 
  if (verbose) write(*,*)"read ", Dtags(Nbl)

! read xHeII (OPTIONAL)
!-----------------------------------------------------------!  
  Nbl=13
  if (DataInSnap(Nbl)) then
     call read_particle_rblock(lun,Nbl,N,rblk)
#ifndef incHe 
     call snap_has_xtra_data("He","incHe")  
#endif
  end if

#ifdef incHe
  if (.not. DataInSnap(Nbl)) call make_want_data("xHeII","incHe")
  psys%par(1:N)%xHeII = rblk(1:N) 
  if (verbose) write(*,*)"read ", Dtags(Nbl)
#endif        

! read xHeIII (OPTIONAL)
!-----------------------------------------------------------!  
  Nbl=14
  if (DataInSnap(Nbl)) call read_particle_rblock(lun,Nbl,N,rblk)

#ifdef incHe
  if (.not. DataInSnap(Nbl)) call make_want_data("xHeIII","incHe")
  psys%par(1:N)%xHeIII = rblk(1:N) 
  if (verbose) write(*,*)"read ", Dtags(Nbl)
#endif        

! read xHIIrc (OPTIONAL) 
!-----------------------------------------------------------!  
  Nbl=15
  if (DataInSnap(Nbl)) then
     call read_particle_rblock(lun,Nbl,N,rblk)
#ifndef increc 
     call snap_has_xtra_data("xHIIrc","increc")  
#endif
  end if

#ifdef increc
!  if (.not. DataInSnap(Nbl)) call make_want_data("xHIIrc","increc")
  psys%par(1:N)%xHIIrc = 0.0 ! rblk(1:N) 
  if (verbose) write(*,*)"set xHIIrc to zero"
#endif        

! read xHeIIrc (OPTIONAL)
!-----------------------------------------------------------!  
  Nbl=16
  if (DataInSnap(Nbl)) then
     call read_particle_rblock(lun,Nbl,N,rblk)
#ifndef incHe 
     call snap_has_xtra_data("He recomb.","incHe")
#endif
#ifndef increc 
     call snap_has_xtra_data("He recomb.","increc")
#endif
  end if

#ifdef incHe
#ifdef increc
  if (.not. DataInSnap(Nbl)) call make_want_data("xHeIIrc","incHe increc")
  psys%par(1:N)%xHeIIrc = rblk(1:N) 
  if (verbose) write(*,*)"read ", Dtags(Nbl)
#endif      
#endif  

! read xHeIIIrc (OPTIONAL)
!-----------------------------------------------------------!  
  Nbl=17
  if (DataInSnap(Nbl)) call read_particle_rblock(lun,Nbl,N,rblk)

#ifdef incHe
#ifdef increc
  if (.not. DataInSnap(Nbl)) call make_want_data("xHeIIIrc","incHe increc")
  psys%par(1:N)%xHeIIIrc = rblk(1:N) 
  if (verbose) write(*,*)"read ", Dtags(Nbl)
#endif      
#endif  

! read last hit (REQUIRED)
!-----------------------------------------------------------!  
  Nbl=18
  if (.not. DataInSnap(Nbl)) call need_datablock(Nbl)
  call read_particle_iblock(lun,Nbl,N,iblk)
  psys%par(1:N)%lasthit = iblk(1:N) 
  if (verbose) write(*,*)"read ", Dtags(Nbl)



  deallocate (rblk)
  deallocate (iblk)

  close(lun)

  if (Verbose) write(*,*)"N particles read in = ", Nall

end subroutine read_initial_par_snapshot



!> reads in an update snapshot with a possible reordering of the 
!! particle data.  calls read_initial_snapshot for the first read in.
!! This routine only reads in the data fields that SPHRAY does not
!! track (IDs, positions, velocities, smoothing lengths, densities
!! and mass)
!===========================================================================
subroutine read_order_par_snapshot(path,base,Ns,verbose,first,psys,MB)
use sort_mod, only: mrgrnk

  character(200), intent(in) :: path  !< path to snapshot dir
  character(200), intent(in) :: base  !< base name of snapshots
  integer, intent(in) :: Ns           !< snapshot number
  logical :: verbose                  !< report to screen?
  logical :: first                    !< initial readin?
  type(particle_system_type), intent(inout) :: psys  !< particle system to load
  real, intent(out) :: MB   !< MBs allocated to particle storage

  real(r4b), allocatable     :: rblk(:)
  integer(i8b), allocatable  :: iblk(:)


  type(particle_system_type) :: psysnew
  type(particle_header_type) :: phead
  logical :: reorder

  character(200) :: snapfile
  integer :: lun,err
  integer :: N, Nall

  logical :: DataInSnap(NdataBlocks)
  integer :: bytesperpar, Nbl
  logical :: closefile
  logical :: header_verbose
  integer :: fn

#ifdef incid
  integer(i8b) :: i, id, indx
  integer(i8b), allocatable  :: onow(:)
#endif


  ! if this is the first readin call the simple snapshot readin routine
  if (first) then
     call read_initial_par_snapshot(path,base,Ns,verbose,psys,MB)
     return
  end if

  fn=1
  call form_snapshot_file_name(path,base,Ns,fn,snapfile)
  write(*,'(A,A)') "updating particle system from ", trim(snapfile)

  closefile=.false.
  header_verbose = .false.
  call read_particle_header(snapfile,header_verbose,closefile,phead,lun)
  DataInSnap = .false.
  where (phead%DataPresent==1) DataInSnap=.true.

  psys%box%bot(1:3) = phead%real(5:7)%var  
  psys%box%top(1:3) = phead%real(8:10)%var 
  ! bndry conditions are set in main_input.f90

  Nall = phead%int(1)%var    ! NparSnap
  N = phead%int(2)%var       ! NparFile

  ! check if we need to reorder particles.  This needs to happen if
  ! the order of the particles in the snapshot has changed OR if the 
  ! number of particles in the snapshot has changed.  If we are using
  ! IDs then we will assume the order is different.  This check cannot
  ! tell if the number of particles is the same but the order has changed
  ! as well.  It is up to the user to include particle IDs if the number
  ! of particles remains the same, but the order in the snapshot changes

  reorder = .false.
#ifdef incid
  reorder = .true.
#endif

  if (psys%npar .NE. Nall) then
     write(*,*) "number of gas particles has changed"
     reorder = .true.
#ifndef incID
     write(*,*) "particle IDs have not been used, cannot proceed"
     stop
#endif
  end if


  if (reorder) then
     call calc_bytes_per_particle(bytesperpar)
     MB = bytesperpar * real(Nall) / 2**20

     100  format(A,I3)
     101  format(A,F10.4,A,I10,A)
     write(*,100) "bytes per particle = ", bytesperpar
     write(*,101) "allocating ", MB, " MB for update of", Nall, " particles"

     if (Nall > psys%npar) then
        write(*,*) "Npar in new snap > Npar in old snap"
        write(*,*) "Nall new = ", Nall
        write(*,*) "Nall old = ", psys%npar
        stop
     end if
     allocate (psysnew%par(Nall))
  end if

  ! allocate arrays for reading in new data blocks 
  allocate (rblk(N), stat=err)
  if (err/=0) stop "cant allocate update rblock in sphray_input.f90"
  allocate(iblk(N), stat=err)
  if (err/=0) stop "cant allocate update iblock in sphray_input.f90"

! read id's (OPTIONAL)
!-----------------------------------------------------------!  
  Nbl=1
  if(reorder .and. .not. DataInSnap(Nbl)) call need_datablock(Nbl)
  if (DataInSnap(Nbl)) call read_particle_iblock(lun,Nbl,N,iblk)
#ifdef incid
     psysnew%par(:)%id = iblk
     if (verbose) write(*,*)"read ", Dtags(Nbl)
#endif
 

! read positions (REQUIRED)
!-----------------------------------------------------------!  
  Nbl=2
  if (.not. DataInSnap(Nbl)) call need_datablock(Nbl)
  call read_particle_rblock(lun,Nbl,N,rblk)
  if(reorder) then
     psysnew%par(:)%pos(1) = rblk
  else
     psys%par(:)%pos(1) = rblk
  endif
  if (verbose) write(*,*)"updated ", Dtags(Nbl), " from snapshot"

  Nbl=3
  if (.not. DataInSnap(Nbl)) call need_datablock(Nbl)
  call read_particle_rblock(lun,Nbl,N,rblk)
  if(reorder) then
     psysnew%par(:)%pos(2) = rblk
  else
     psys%par(:)%pos(2) = rblk
  endif
  if (verbose) write(*,*)"updated ", Dtags(Nbl), " from snapshot"

  Nbl=4
  if (.not. DataInSnap(Nbl)) call need_datablock(Nbl)
  call read_particle_rblock(lun,Nbl,N,rblk)
  if(reorder) then
     psysnew%par(:)%pos(3) = rblk
  else
     psys%par(:)%pos(3) = rblk
  endif 
  if (verbose) write(*,*)"updated ", Dtags(Nbl), " from snapshot"

! read velocities (OPTIONAL)
!-----------------------------------------------------------!  
  Nbl=5
  if (DataInSnap(Nbl)) call read_particle_rblock(lun,Nbl,N,rblk)

#ifdef incvel
  if (.not. DataInSnap(Nbl)) call make_want_data("x-vel","incvel")
  if(reorder) then
     psysnew%par(:)%vel(1) = rblk
  else
     psys%par(:)%vel(1) = rblk
  end if
  if (verbose) write(*,*)"updated ", Dtags(Nbl), " from snapshot"
#endif        
        
  Nbl=6
  if (DataInSnap(Nbl)) call read_particle_rblock(lun,Nbl,N,rblk)

#ifdef incvel
  if (.not. DataInSnap(Nbl)) call make_want_data("y-vel","incvel")
  if(reorder) then
     psysnew%par(:)%vel(2) = rblk
  else
     psys%par(:)%vel(2) = rblk
  end if
  if (verbose) write(*,*)"updated ", Dtags(Nbl), " from snapshot"
#endif        

  Nbl=7
  if (DataInSnap(Nbl)) call read_particle_rblock(lun,Nbl,N,rblk)

#ifdef incvel
  if (.not. DataInSnap(Nbl)) call make_want_data("z-vel","incvel")
  if(reorder) then
     psysnew%par(:)%vel(3) = rblk
  else
     psys%par(:)%vel(3) = rblk
  end if
  if (verbose) write(*,*)"updated ", Dtags(Nbl), " from snapshot"
#endif        

! read smoothing lengths (REQUIRED)
!-----------------------------------------------------------!  
  Nbl=8
  if (.not. DataInSnap(Nbl)) call need_datablock(Nbl)
  call read_particle_rblock(lun,Nbl,N,rblk)
  if(reorder) then
     psysnew%par(:)%hsml = rblk
  else
     psys%par(:)%hsml = rblk
  endif
  if (verbose) write(*,*)"updated ", Dtags(Nbl), " from snapshot"

! read density (REQUIRED)
!-----------------------------------------------------------!  
  Nbl=9
  if (.not. DataInSnap(Nbl)) call need_datablock(Nbl)
  call read_particle_rblock(lun,Nbl,N,rblk)
  if(reorder) then
     psysnew%par(:)%rho = rblk
  else
     psys%par(:)%rho = rblk
  endif
  if (verbose) write(*,*)"updated ", Dtags(Nbl), " from snapshot"

! read masses (OPTIONAL)  
!-----------------------------------------------------------!  
  Nbl=10
  if (DataInSnap(Nbl)) call read_particle_rblock(lun,Nbl,N,rblk)

#ifdef incmass
  if (.not. DataInSnap(Nbl)) call make_want_data("mass","incmass")
  if(reorder) then
     psysnew%par(:)%mass = rblk
  else
     psys%par(:)%mass = rblk
  endif
  if (verbose) write(*,*)"updated ", Dtags(Nbl), " from snapshot"
#endif

  ! deallocate readin arrays
  deallocate( iblk, rblk )

#ifdef incid

     ! create an indexing array for the current particle system id's
     ! id's that are no longer present are indexed with a GV%minid-1
     write(*,'(A)') "creating indexing array for current particle system..."
     write(*,*) "min/max ID (first snap) = ", psys%minIDfirst, psys%maxIDfirst
     allocate(onow(psys%minIDfirst:psys%maxIDfirst), stat=err)
     if(err/=0) stop "sphray_input.f90> cant allocate onow array"
     onow = psys%minIDfirst-1
     do i = 1,psys%npar
        onow(psys%par(i)%id) = i
     end do
     write(*,'(A)') "done"

     ! transfer all carryover variables
     do i = 1,Nall

        id = psysnew%par(i)%id   
        indx = onow(id)
     
        if (id.NE.psys%par(indx)%id) then
           write(*,*) "Problem transfering new particles into psys"
           write(*,*) "idnew,idold", id, psys%par(indx)%id
           write(*,*) "i = ", i
           stop
        end if

#ifdef incT
        psysnew%par(i)%T = psys%par(indx)%T 
#endif

        psysnew%par(i)%xHII = psys%par(indx)%xHII 
#ifdef incHe
        psysnew%par(i)%xHeII = psys%par(indx)%xHeII 
        psysnew%par(i)%xHeIII = psys%par(indx)%xHeIII
#endif

#ifdef increc
        psysnew%par(i)%xHIIrc = psys%par(indx)%xHIIrc
#ifdef incHe
        psysnew%par(i)%xHeIIrc = psys%par(indx)%xHeIIrc
        psysnew%par(i)%xHeIIIrc = psys%par(indx)%xHeIIIrc
#endif
#endif

        psysnew%par(i)%lasthit = psys%par(indx)%lasthit 

     end do

     ! move data into shrunken psys
     deallocate( psys%par )
     allocate( psys%par(Nall) )

     do i = 1,Nall
        psys%par(i) = psysnew%par(i)
     end do

     psys%npar = Nall

     deallocate( psysnew%par )
     deallocate (onow)

#endif

  close(lun)  


end subroutine read_order_par_snapshot

!>   reads a source snapshot into psys (psys%src is allocated in this routine)
!=============================================================================
subroutine read_src_snapshot(path,base,Ns,verbose,psys,rayn,MB)

  character(200), intent(in) :: path  !< path to source snapshot
  character(200), intent(in) :: base  !< base name of snapshot
  integer, intent(in) :: Ns           !< snapshot number
  logical, intent(in) :: verbose      !< report to screen?
  type(particle_system_type), intent(inout) :: psys !< particle system to load
  integer, intent(in) :: rayn   !< ray number at readin
  real, intent(out) :: MB  !< MB's allocated to source storage

  character(200) :: snapfile  
  type(source_header_type) :: shead
  integer :: lun, err
  integer :: i, N, Nall
  integer :: bytespersource
  logical :: closefile
  integer :: fn

  fn=1
  call form_snapshot_file_name(path,base,Ns,fn,snapfile)
  write(*,'(A,A)') "reading in source snapshot from ", trim(snapfile)
      
  bytespersource = 6 * 4 + 8
  closefile=.false.
  call read_source_header(snapfile,verbose,closefile,shead,lun)

  N = shead%NsrcFile
  Nall = shead%NsrcSnap        
  MB = bytespersource * real(Nall) / 2**20
100 format(A,I3)
101 format(A,F10.4,A,I10,A)
  write(*,100) "bytes per source = ", bytespersource
  write(*,101) "allocating ", MB, " MB for ", Nall, " sources"

  if (allocated(psys%src)) deallocate(psys%src)
  allocate (psys%src(Nall),stat=err)
  if (err /= 0) stop "cant allocate sources in particle system"
  psys%nsrc = Nall
     
  do i = 1,N 
     read(lun,*) psys%src(i)%pos(1), & 
                 psys%src(i)%pos(2), & 
                 psys%src(i)%pos(3), &
                 psys%src(i)%L, &
                 psys%src(i)%SpcType, &
                 psys%src(i)%EmisPrf
  end do
  close(lun)
 
  psys%src(1:N)%lastemit = rayn

  if (verbose) call source_info_to_screen(psys)
  
end subroutine read_src_snapshot



!===========================================================================
!===========================================================================
!
! Screen Output for Dummy Checking
!
!===========================================================================
!===========================================================================


!> outputs particle snapshot header to the screen
!-----------------------------------------------
  subroutine particle_header_to_screen(phead)
    type(particle_header_type), intent(in) :: phead !< particle header to print
    integer(i4b) :: i

    97 format(T2,"=",T58,"=")
    98 format(T2,"= ",A,T58,"=")
    99 format(T2,57("="))
    100 format(T2,"=",T4,I3,T8,A,I12,T58,"=")
    200 format(T2,"=",T4,I3,T8,A,ES12.5,T58,"=")
    201 format(T2,"=",T4,I3,T8,A,T58,"=")

    write(*,*)
    write(*,99)
    write(*,98) "begin particle header data"
    write(*,99)
    write(*,97) 

    do i = 1,Niheadp
       write(*,100) i,phead%int(i)%tag,phead%int(i)%var
    end do
    write(*,97) 
    do i = 1,Nrheadp
       write(*,200) i,phead%real(i)%tag,phead%real(i)%var
    end do
    write(*,97) 
    write(*,98) "The snapshot contains the following data blocks"
    do i = 1,ndatablocks
       if (phead%DataPresent(i)==1) write(*,201) i, Dtags(i)
    end do
    write(*,97)

    write(*,99) 
    write(*,98) "end particle header data"
    write(*,99) 
    write(*,*) 


end subroutine particle_header_to_screen


!> outputs source snapshot header to the screen
!-----------------------------------------------
subroutine source_header_to_screen(shead)
  type(source_header_type), intent(in) :: shead !< source header to print

97 format(T2,"=",T58,"=")
98 format(T2,"= ",A,T58,"=")
99 format(T2,57("="))
100 format(T2,"=",T4,A,T30,I6,T58,"=")
200 format(T2,"=",T4,A,T30,ES12.5,T58,"=")
201 format(T2,"=",T4,A,T58,"=")
  
  write(*,99) 
  write(*,98) "source header data"
  write(*,99) 
  write(*,97) 
  write(*,100) "sources in snapshot:", shead%NsrcSnap
  write(*,100) "sources in file:", shead%NsrcFile
  write(*,100) "number of files in snap:", shead%NumFiles
  write(*,200) "rays for this snap:*", shead%TotalRays
  write(*,201) "*(only if RayScheme='header' in config file)"
  write(*,97) 
  write(*,99) 
  
end subroutine source_header_to_screen


!> puts code units read in from a particle snapshot header
!! into the global variables
!=========================================================
subroutine get_planning_data_sphray(verbose,GV)
use physical_constants_mod

  logical, intent(in) :: verbose   !< report to screen?
  type(global_variables_type), intent(inout) :: GV !< global variables

  logical :: closefile
  integer :: lun
  type(particle_header_type) :: phead
  type(source_header_type) :: shead

  integer :: iSnap, fSnap    ! initial and final snapshot numbers
  integer :: pfiles, sfiles  ! files/snap for particles and sources   
  integer :: i,j             ! counters
  character(200) :: snapfile ! snapshot file name


  iSnap = GV%StartSnapNum
  fSnap = GV%EndSnapNum
    
  pfiles = GV%ParFilesPerSnap
  sfiles = GV%SourceFilesPerSnap
    

  ! read in and report to screen
  !==============================
  closefile = .true.
  write(*,*)    
  ! first just report the files that will be read
  write(*,'(A)') "reading all SPHRAY particle header(s) ... "
    do i = iSnap,fSnap
       do j = 1,pfiles
          call form_snapshot_file_name(GV%SnapPath,GV%ParFileBase,i,j,snapfile)
          write(*,'(I3,"  ",A)') i,trim(snapfile)
          call read_particle_header(snapfile,verbose,closefile,phead,lun)


          GV%TimeAtSnap(i)     = phead%real(1)%var      
          ! IsoMass comes from config file
          ! IsoTemp comes from config file
          GV%ScalefacAtSnap(i) = phead%real(4)%var      
          GV%BoxLows(:)  = phead%real(5:7)%var
          GV%BoxHighs(:) = phead%real(8:10)%var
          GV%cgs_len     = phead%real(11)%var
          GV%cgs_lum     = phead%real(12)%var
          GV%cgs_time    = phead%real(13)%var
          GV%OmegaM      = phead%real(14)%var
          GV%OmegaB      = phead%real(15)%var
          GV%OmegaL      = phead%real(16)%var
          GV%LittleH     = phead%real(17)%var

       end do
    end do

    GV%cgs_mass = GV%cgs_lum * GV%cgs_time**3 / GV%cgs_len**2
    GV%cgs_enrg = GV%cgs_lum * GV%cgs_time
    GV%cgs_vel  = GV%cgs_len / GV%cgs_time

    write(*,*) 
    write(*,'(A)') "reading all SPHRAY source header(s) ... "
    do i = iSnap,fSnap
       do j = 1,sfiles
          call form_snapshot_file_name(GV%SourcePath,GV%SourceFileBase,i,j,snapfile)
          write(*,'(I3,"  ",A)') i,trim(snapfile)
          call read_source_header(snapfile,verbose,closefile,shead,lun)

          GV%RaysFromSrcHeader(i) = shead%TotalRays

       end do
    end do

   if (verbose) then
      106   format(A,F12.5,A)
      write(*,'(A)') "setting code units ..."
      write(*,106) "  length:     ", GV%cgs_len * cm2kpc, " kpc"
      write(*,106) "  luminosity: ", GV%cgs_lum / L_solar, " solar"
      write(*,106) "  time:       ", GV%cgs_time * s2Myr, " Myr"
      write(*,106) "  mass:       ", GV%cgs_mass / M_solar, " solar"
      write(*,106) "  velocity    ", GV%cgs_vel / km2cm, " km/s"
      write(*,106) "  energy:     ", GV%cgs_enrg / solarLMyrs, &
                                    " solar luminosity * Myr"
      write(*,*) 
   end if

 end subroutine get_planning_data_sphray


!> forms a snapshot name from a path, a file base, a snapshot number
!! and a file number.
!===================================================================
subroutine form_snapshot_file_name(Path,FileBase,SnapNum,FileNum,SnapFile)
  character(200), intent(in) :: Path       !< path to snapshot dir
  character(200), intent(in) :: FileBase   !< file base names
  integer, intent(in) :: SnapNum           !< snapshot number
  integer, intent(in) :: FileNum           !< file number in snapshot
  character(200), intent(out) :: SnapFile  !< file name to return
  
  character(10) :: FileNumChar
  
  write(FileNumChar,"(I6)") FileNum
  
100 format(A,"/",A,"_",I3.3,".",A)
  write(SnapFile,100) trim(Path), trim(FileBase), SnapNum, &
       trim(adjustl(FileNumChar))
  
end subroutine form_snapshot_file_name




end module sphray_input_mod
