!> \file main_input.F90

!> \brief The module that calls the specific input routines
!<

module main_input_mod
use myf90_mod
use gadget_input_mod
use gadget_input_hdf5_mod
use source_input_mod
use particle_system_mod, only: particle_type
use particle_system_mod, only: source_type
use particle_system_mod, only: box_type
use particle_system_mod, only: scale_comoving_to_physical
use particle_system_mod, only: scale_physical_to_comoving
use particle_system_mod, only: set_ye
use particle_system_mod, only: enforce_x_and_T_minmax
use particle_system_mod, only: particle_info_to_screen
use atomic_rates_mod, only: get_atomic_rates
use global_mod, only: PLAN, GV, rtable, cmbT_k
use global_mod, only: psys, gconst, saved_gheads
implicit none

contains


!> Read in planning data from the header of all snapshots 
!========================================================
subroutine get_planning_data(skewers)
  logical, optional, intent(in) :: skewers
  
  character(clen), parameter :: myname="get_planning_data"
  logical, parameter :: crash=.true.
  integer, parameter :: verb=1

  logical :: skew


  if (present(skewers)) then
     skew = skewers
  else
     skew = .false.
  endif

  call mywrite("getting planning data:", verb)
  call mywrite("",verb) 


  GV%Nsnaps = GV%EndSnapNum - GV%StartSnapNum + 1
  if (allocated(PLAN%snap)) deallocate(PLAN%snap)
  allocate( PLAN%snap(GV%StartSnapNum : GV%EndSnapNum) )

  if (GV%InputType == 1 .or. GV%InputType == 2 .or. GV%InputType == 4) then
     call get_planning_data_gadget()
  else if (GV%InputType == 3) then
     call get_planning_data_gadget_hdf5()
  end if


  if (.not. skew) then
     call get_planning_data_sources()
  endif

end subroutine get_planning_data


!> read in particle, box, and source data 
!============================================
subroutine readin_snapshot(skewers)

  
  logical, optional, intent(in) :: skewers

  character(clen), parameter :: myname="readin_snapshot"
  logical, parameter :: crash=.true.
  integer, parameter :: verb=2
  character(clen) :: str,fmt
  
  logical :: skew
  logical :: first
  real(r8b) :: MB
  integer(i8b) :: i
  real(r8b) :: a     !< scale factor
  real(r8b) :: h     !< Hubble paraemter (little H)
  real(r8b) :: Flux  !< photons / s from planar sources
  character(clen) :: snpbase
  character(clen) :: srcbase
  real :: cm2kpc


  if (present(skewers)) then
     skew = skewers
  else
     skew = .false.
  endif
  
  call mywrite("reading in particle and source snapshots:", verb-1)
  call mywrite("",verb-1) 

  ! set local variables
  !======================
  a = PLAN%snap(GV%CurSnapNum)%ScalefacAt 
  h = GV%LittleH

  ! report readin type
  !======================
  call mywrite('   input type = ', verb, fmt='(A)', adv=.false.)
  if (GV%InputType==1) then
     call mywrite(" Gadget-2 public (SnapFormat=1)", verb)
  else if (GV%InputType==2) then
     call mywrite(" Gadget w/ cooling, stars, UVB (SnapFormat=1)", verb)
  else if (GV%InputType==3) then
     call mywrite(" Gadget HDF5", verb)
  else if (GV%InputType==4) then
     call mywrite(" Gadget Bromm", verb)
  end if



  ! read in the particle data
  !=============================
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
        call read_Gpub_particles()
        psys%par(:)%lasthit = 0
     else
        call update_particles()
     end if

  ! gadget w/ cooling (i.e. ye and xHI)
  !---------------------------------------------------------------
  else if (GV%InputType == 2) then

     if (first) then
        call read_Gcool_particles()
        psys%par(:)%lasthit = 0
     else
        call update_particles()
     end if
     
  ! gadget w/ HDF5
  !---------------------------------------------------------------
  else if (GV%InputType == 3) then 
     
     if (first) then
        call read_Ghdf5_particles()
        psys%par(:)%lasthit = 0
     else
        call update_particles()
     end if
     
  ! gadget w/ ions from Volker Bromm's group 
  !---------------------------------------------------------------
  else if (GV%InputType == 4) then

     if (first) then
        call read_Gbromm_particles()
        psys%par(:)%lasthit = 0
     else
        call update_particles()
     end if

  ! not recognized
  !---------------------------------------------------------------
  else
     write(str,*) "input type, ", GV%InputType, "not recognized" 
     call myerr(str,myname,crash)
  end if


  ! copy over box properties
  !==================================================
  psys%box%top = GV%BoxUprsComoh
  psys%box%bot = GV%BoxLwrsComoh
  psys%box%tbound = GV%BndryCond
  psys%box%bbound = GV%BndryCond

  ! tighten box - add this to config file also

!  do i = 1,3
!     psys%box%top(i) = maxval( psys%par(:)%pos(i) )
!     psys%box%bot(i) = minval( psys%par(:)%pos(i) )
!     GV%BoxUprsComoh(i) = psys%box%top(i)
!     GV%BoxLwrsComoh(i) = psys%box%bot(i)
!  end do


  ! read in the source data
  !============================================
  if (.not. skew) then
     call read_src_snapshot()
     call order_sources_lum(psys%src)
     psys%src%lastemit = GV%rayn
  endif


  ! set par and src file bases for output to logfiles
  !====================================================
  fmt = "(A,'/',A,'_',I3.3)"
  write(snpbase,fmt) trim(GV%SnapPath),   trim(GV%ParFileBase),    GV%CurSnapNum
  write(srcbase,fmt) trim(GV%SourcePath), trim(GV%SourceFileBase), GV%CurSnapNum
 

  ! write fresh reads to the particle_data.log and source_data.log files
  !========================================================================  
  fmt = "(A,A)"
  write(str,fmt) "Fresh read from ", trim(snpbase)
  call particle_info_to_screen(psys,str,GV%pardatalun)
  write(GV%pardatalun,*)
  write(GV%pardatalun,*)

  if (.not. skew) then
#ifdef hdf5
     write(GV%srcdatalun,*) "================================================================="
     write(GV%srcdatalun,*) " HM01 G+C gammaHI for z = ", saved_gheads(GV%CurSnapNum,0)%z, ": ", &
                              GV%UVB_gammaHI_cloudy
     write(GV%srcdatalun,*) "================================================================="
     write(GV%srcdatalun,*) 
#endif
     write(str,fmt) "Fresh read from ", trim(srcbase)
     call source_info_to_screen(psys,str,GV%srcdatalun)
     write(GV%srcdatalun,*)
     write(GV%srcdatalun,*)
  endif

  ! take care of all the box variables
  !===============================================================
  cm2kpc = 1.0d0 / (gconst%cm_per_mpc * 1.0d-3)

  GV%BoxLwrsComo = GV%BoxLwrsComoh / h
  GV%BoxLwrsPhys = GV%BoxLwrsComoh / h * a

  GV%BoxUprsComo = GV%BoxUprsComoh / h
  GV%BoxUprsPhys = GV%BoxUprsComoh / h * a
 
  GV%BoxLensComoh = GV%BoxUprsComoh - GV%BoxLwrsComoh
  GV%BoxLensComo  = GV%BoxUprsComo  - GV%BoxLwrsComo
  GV%BoxLensPhysh = GV%BoxUprsPhysh - GV%BoxLwrsPhysh
  GV%BoxLensPhys  = GV%BoxUprsPhys  - GV%BoxLwrsPhys

  GV%BoxLensComoh_cm = GV%BoxLensComoh * GV%cgs_len
  GV%BoxLensComo_cm  = GV%BoxLensComo  * GV%cgs_len
  GV%BoxLensPhysh_cm = GV%BoxLensPhysh * GV%cgs_len
  GV%BoxLensPhys_cm  = GV%BoxLensPhys  * GV%cgs_len
 
  GV%BoxLensComoh_kpc = GV%BoxLensComoh_cm * cm2kpc
  GV%BoxLensComo_kpc  = GV%BoxLensComo_cm  * cm2kpc
  GV%BoxLensPhysh_kpc = GV%BoxLensPhysh_cm * cm2kpc
  GV%BoxLensPhys_kpc  = GV%BoxLensPhys_cm  * cm2kpc

  GV%BoxVolComoh = product( GV%BoxLensComoh )
  GV%BoxVolComo  = product( GV%BoxLensComo  )
  GV%BoxVolPhysh = product( GV%BoxLensPhysh )
  GV%BoxVolPhys  = product( GV%BoxLensPhys  )

  GV%BoxVolComoh_cm = product( GV%BoxLensComoh_cm )
  GV%BoxVolComo_cm  = product( GV%BoxLensComo_cm  )
  GV%BoxVolPhysh_cm = product( GV%BoxLensPhysh_cm )
  GV%BoxVolPhys_cm  = product( GV%BoxLensPhys_cm  )

  GV%BoxVolComoh_kpc = product( GV%BoxLensComoh_kpc )
  GV%BoxVolComo_kpc  = product( GV%BoxLensComo_kpc  )
  GV%BoxVolPhysh_kpc = product( GV%BoxLensPhysh_kpc )
  GV%BoxVolPhys_kpc  = product( GV%BoxLensPhys_kpc  )


  ! convert number density to flux for planar sources
  !==========================================================
  if (.not. skew) then
     do i = 1,size(psys%src)
        
        if (psys%src(i)%EmisPrf == -3 .or. &
            psys%src(i)%EmisPrf == -2 .or. &
            psys%src(i)%EmisPrf == -1) then 
           ! if this is true the input luminosity is a number density 
           ! [photons/cm^3].  we want the flux that would produce this 
           ! number density in an optically thin volume

           write(*,*) 
           write(*,*) "  converting a photon number density to a flux"
           write(*,*) "  n_photon/cm^3                = ", psys%src(i)%L
           
           Flux = psys%src(i)%L * gconst%c * GV%BoxLensPhys_cm(1)**2 
           Flux = Flux / GV%Lunit
           psys%src(i)%L = Flux

           write(*,*) "  photons/s from walls [1.e50] = ",  psys%src(i)%L
           write(*,*)            

        end if

     end do
  endif

  ! check test conditionals 
  !==========================================================
  if (GV%DoTestScenario) then

     if ( trim(GV%TestScenario) == "iliev_test1") then
        psys%par(:)%xHII = 1.2d-3
        psys%par(:)%xHI = 1.0d0 - psys%par(:)%xHII
        psys%par(:)%T = 1.0d4

     else if ( trim(GV%TestScenario) == "iliev_test2" ) then
        psys%par(:)%xHII = 0.0d0
        psys%par(:)%xHI = 1.0d0 - psys%par(:)%xHII
        psys%par(:)%T = 1.0d2

     else if ( trim(GV%TestScenario) == "iliev_test1He" ) then
        psys%par(:)%xHI = 1.0d0
        psys%par(:)%xHII = 0.0d0
#ifdef incHe
        psys%par(:)%xHeI = 1.0d0
        psys%par(:)%xHeII = 0.0d0
        psys%par(:)%xHeIII = 0.0d0
#endif
        psys%par(:)%T = 1.0d4

     end if
  end if


  ! set EOS particles to EOS temp if you want
  !=======================================================
#ifdef incEOS
  if (first) then
     do i = 1, size(psys%par(:))
        if (GV%EOStemp > 0.0) then
           if (psys%par(i)%eos == 1.0) then
              psys%par(i)%T = GV%EOStemp
           endif
        endif
     enddo
  endif
#endif



  ! set neutral or ionized if we need to
  !=======================================================
if (first) then

   if (GV%InitxHI > 0.0) then
      psys%par(:)%xHI  = GV%InitxHI
      psys%par(:)%xHII = 1.0d0 - GV%InitxHI
      call set_ye(psys, GV%H_mf, GV%He_mf, GV%NeBackGround)
   endif
   
endif



  ! set constant temperature if we have one
  !=======================================================
  if (GV%IsoTemp /= 0.0) psys%par(:)%T = GV%IsoTemp


  ! these quantities track the photoionization rate.  they are 
  ! rezeroed at inputs (because new source files are loaded) and 
  ! outputs (for time dependence)
  !===============================================================
#ifdef outGammaHI
  psys%par(:)%gammaHI = 0.0
  psys%par(:)%time = 0.0
#endif

  ! cap the ionization fractions and temperatures if we have to
  !================================================================
  call enforce_x_and_T_minmax(psys%par, GV%xfloor, GV%xceiling, GV%Tfloor, GV%Tceiling)

  ! write data after above conditionals to the particle and source log files
  !==========================================================================
  fmt = "(A,A)"
  write(str,fmt) "After test conditionals from ", trim(snpbase)
  call particle_info_to_screen(psys,str,GV%pardatalun)
  write(GV%pardatalun,*)
  write(GV%pardatalun,*)
 
  if (.not. skew) then
     write(str,fmt) "After test conditionals from ", trim(srcbase)
     call source_info_to_screen(psys,str,GV%srcdatalun)
     write(GV%srcdatalun,*)
     write(GV%srcdatalun,*)
  endif

  ! scale the data if we need to
  !=====================================================================
  if(GV%Comoving) then
     call scale_comoving_to_physical(a, psys%par, psys%src, psys%box, h )
  endif

  ! write data after rescaling to the particle_data.log file
  !=====================================================================
  fmt = "(A,F5.3,A,F5.3,A,A)"
  write(str,fmt) "After rescaling (a=",a,",h=",h,") from ", trim(snpbase)
  call particle_info_to_screen(psys,str,GV%pardatalun)
  write(GV%pardatalun,*)
  write(GV%pardatalun,*)

  if (.not. skew) then
     write(str,fmt) "After rescaling (a=",a,",h=",h,") from ", trim(srcbase)
     call source_info_to_screen(psys,str,GV%srcdatalun)
     write(GV%srcdatalun,*)
     write(GV%srcdatalun,*)
  endif


   

  

  ! and the rest of the stuff
  !===============================================================
  if (.not. skew) then 
     GV%dtray_code = PLAN%snap(GV%CurSnapNum)%RunTime / PLAN%snap(GV%CurSnapNum)%SrcRays
     GV%dtray_s    = GV%dtray_code * GV%cgs_time / h
  endif

  GV%Tcmb_cur = gconst%t_cmb0 / a
  call get_atomic_rates(GV%Tcmb_cur, rtable, cmbT_k)

  GV%total_mass = 0.0d0
  do i = 1,size(psys%par)
     GV%total_mass = GV%total_mass + psys%par(i)%mass
  end do
  
  GV%total_lum = 0.0d0
  do i = 1,size(psys%src)
     GV%total_lum = GV%total_lum + psys%src(i)%L
  end do
  
  GV%total_atoms = GV%total_mass * GV%cgs_mass * &
       (GV%H_mf  / (gconst%protonmass) + &
        GV%He_mf / (4*gconst%protonmass) )


  GV%total_photons = (GV%TotalSimTime * GV%cgs_time / GV%LittleH) * (GV%total_lum * GV%Lunit)

  
  ! write some final data to the source log file
  !=====================================================================
  if (.not. skew) then
     fmt = "(A,A)"
     100  format(72("-"))
     
     write(GV%srcdatalun,100)
     write(GV%srcdatalun,fmt) "Ray / Luminosity info from ", trim(srcbase)
     write(GV%srcdatalun,*) 
     write(GV%srcdatalun,'(A,ES12.5)') "dt/ray [code] = ", GV%dtray_code
     write(GV%srcdatalun,'(A,ES12.5)') "dt/ray [s]    = ", GV%dtray_s
     write(GV%srcdatalun,'(A,ES12.5)') "dt/ray [Myr]  = ", GV%dtray_s / gconst%sec_per_megayear
     write(GV%srcdatalun,*)
     write(GV%srcdatalun,'(A,ES12.5)') "total photons = ", GV%total_photons
     write(GV%srcdatalun,'(A,ES12.5)') "total atoms   = ", GV%total_atoms
     write(GV%srcdatalun,'(A,ES12.5)') "photons / atoms = ", GV%total_photons / GV%total_atoms
     write(GV%srcdatalun,*)
     write(GV%srcdatalun,100)
  endif

  call mywrite("",verb)
 

end subroutine readin_snapshot





end module main_input_mod
