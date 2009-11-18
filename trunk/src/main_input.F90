!> \file main_input.F90

!> \brief The module that calls the specific input routines
!<

module main_input_mod
use myf90_mod
use gadget_input_mod
use source_input_mod
use particle_system_mod, only: particle_type
use particle_system_mod, only: source_type
use particle_system_mod, only: box_type
use particle_system_mod, only: scale_comoving_to_physical
use particle_system_mod, only: scale_physical_to_comoving
use global_mod, only: psys, PLAN, GV
implicit none

contains


!> Read in planning data from the header of all snapshots 
!========================================================
subroutine get_planning_data()

  ! allocate global planning arrays
  GV%Nsnaps = GV%EndSnapNum - GV%StartSnapNum + 1

  allocate( PLAN%snap(GV%StartSnapNum : GV%EndSnapNum) )

  write(*,*) 
  write(*,'(A)') "reading headers "

  call get_planning_data_gadget()

end subroutine get_planning_data


!> read in particle, box, and source data 
!============================================
subroutine readin_snapshot(box,MBalloc)
use particle_system_mod, only: enforce_x_and_T_minmax
use particle_system_mod, only: particle_info_to_screen
use physical_constants_mod, only: M_H, erg2eV, CMBtempNow, cm2kpc3, cm2kpc
use atomic_rates_mod, only: set_cmb_atomic_rates

type(box_type), intent(inout) :: box   !< simulation vol. 
real(r8b), intent(out) :: MBalloc      !< MB allocated

logical :: verbose, first
real(r8b) :: MB
integer(i8b) :: i

character(clen), parameter :: myname="readin_snapshot"
logical, parameter :: crash=.true.
integer, parameter :: verb=1

  call mywrite(' input type = ', verb, fmt='(A)', adv=.false.)
  if (GV%InputType==1) then
   call mywrite(" Gadget-2 public (SnapFormat=1)", verb)
  else if (GV%InputType==2) then
     call mywrite(" Gadget w/ cooling, stars, UVB (SnapFormat=1)", verb)
  else if (GV%InputType==3) then
     call mywrite(" Gadget HDF5", verb)
  end if


  ! read in the particle and box data
  !============================================
  MBalloc = 0.0

  verbose=.true.
  first = .false.
  if (GV%CurSnapNum == GV%StartSnapNum) then
     first = .true.
  else
     GV%CurSnapNum = GV%CurSnapNum + 1
  end if


! public gadget
!---------------------------------------------------------------
  if (GV%InputType == 1) then

     if (first) then
        call read_Gpub_particles(MB)
        psys%par(:)%lasthit = 0
     else
        call update_particles(MB)
     end if

     call read_gadget_box(box)

     MBalloc = MBalloc + MB

! gadget w/ cooling, stars, met
!---------------------------------------------------------------
  else if (GV%InputType == 2) then

     if (first) then
        call read_Gtiz_particles(MB)
        psys%par(:)%lasthit = 0
     else
        call update_particles(MB)
     end if
     
     call read_gadget_box(box)

     MBalloc = MBalloc + MB

! gadget w/ HDF5
!---------------------------------------------------------------
  else if (GV%InputType == 3) then 
     
     write(*,*) "set up HDF5 input" 
     stop

! not recognized
!---------------------------------------------------------------
  else
     write(*,*) "<readin_snapshot> input type, ", GV%InputType, "not recognized"
     stop
  end if


  if (GV%DoTestScenario) then 
     if ( trim(GV%TestScenario) == "iliev_test1He" ) then
        psys%par(:)%xHI = 1.0d0
        psys%par(:)%xHII = 0.0d0
#ifdef incHe
        psys%par(:)%xHeI = 1.0d0
        psys%par(:)%xHeII = 0.0d0
        psys%par(:)%xHeIII = 0.0d0
#endif
     end if
  end if



  if (GV%IsoTemp /= 0.0) psys%par(:)%T = GV%IsoTemp

  ! these quantities track the photoionization rate.  they are 
  ! rezeroed at inputs (because new source files are loaded) and 
  ! outputs (for time dependence)
#ifdef outGamma
  psys%par(:)%gammaHI = 0.0
  psys%par(:)%time = 0.0
#endif


  ! read in the source data
  !============================================
  call read_src_snapshot(GV,verbose,MB)
  MBalloc = MBalloc + MB


  call particle_info_to_screen(psys%par)
   
  ! now that we've got the data we will set the global variables
  !===============================================================
  
  GV%dtray_code = PLAN%snap(GV%CurSnapNum)%RunTime / PLAN%snap(GV%CurSnapNum)%SrcRays
  GV%dtray_s    = GV%dtray_code * GV%cgs_time / GV%LittleH

  GV%BoxLowers = box%bot  
  GV%BoxUppers = box%top

  GV%BoxLengths(1:3) = GV%BoxUppers(1:3) - GV%BoxLowers(1:3)
  GV%BoxLengths_cm   = GV%BoxLengths * GV%cgs_len / GV%LittleH
  GV%BoxLengths_kpc  = GV%BoxLengths_cm * cm2kpc 

  GV%BoxVol = GV%BoxLengths(1) * GV%BoxLengths(2) * GV%BoxLengths(3) 
  GV%BoxVol_cm  = GV%BoxLengths_cm(1) * GV%BoxLengths_cm(2) * GV%BoxLengths_cm(3) 
  GV%BoxVol_kpc = GV%BoxLengths_kpc(1) * GV%BoxLengths_kpc(2) * GV%BoxLengths_kpc(3) 
  
  GV%Tcmb_cur = CMBtempNow / PLAN%snap(GV%CurSnapNum)%ScalefacAt
  call set_cmb_atomic_rates(GV%Tcmb_cur)  

  GV%total_mass = 0.0d0
  do i = 1,size(psys%par)
     GV%total_mass = GV%total_mass + psys%par(i)%mass
  end do
  
  GV%total_lum = 0.0d0
  do i = 1,size(psys%src)
     GV%total_lum = GV%total_lum + psys%src(i)%L
  end do
  
  GV%total_atoms = GV%total_mass * GV%cgs_mass / GV%LittleH * (GV%H_mf / M_H + GV%He_mf / M_He)
  GV%total_photons = GV%TotalSimTime * GV%cgs_time / GV%LittleH * GV%total_lum * GV%Lunit
  
  write(*,*) 
  write(*,'(A,ES12.5)') "dt/ray [code] = ", GV%dtray_code
  write(*,'(A,ES12.5)') "dt/ray [s]    = ", GV%dtray_s
  write(*,'(A,ES12.5)') "dt/ray [Myr]  = ", GV%dtray_s * s2Myr
  write(*,*)
  write(*,'(A,ES12.5)') "total photons = ", GV%total_photons
  write(*,'(A,ES12.5)') "total atoms   = ", GV%total_atoms
  write(*,'(A,ES12.5)') "photons / atoms = ", GV%total_photons / GV%total_atoms
  write(*,*)

  call order_sources_lum(psys%src)
  call enforce_x_and_T_minmax(psys%par, GV%xfloor, GV%xceiling, GV%Tfloor, GV%Tceiling)
  if(GV%Comoving) call scale_comoving_to_physical(PLAN%snap(GV%CurSnapNum)%ScalefacAt, psys%par, psys%src, box )
  if (verbose) call particle_info_to_screen(psys%par)  
 

end subroutine readin_snapshot





end module main_input_mod
