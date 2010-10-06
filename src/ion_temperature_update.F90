!> \file ion_temperature_update.F90

!> \brief Module that is the starting point for ionization
!! and temperature updates.  It calls a runge kutta or a 
!! backwards difference formula (depending on the Makefile
!! flags)  
!<
module ion_temperature_update
use myf03_mod
use raylist_mod, only: raylist_type
use ray_mod, only: ray_type
use ray_mod, only: make_recomb_ray, make_probe_ray
use particle_system_mod, only: particle_system_type
use particle_system_mod, only: particle_type
use particle_system_mod, only: box_type
use raylist_mod, only: trace_ray
use oct_tree_mod, only: oct_tree_type
use raylist_mod, only: raylist_type
use ionpar_mod
use spectra_mod, only: rn2freq
use euler_mod, only: eulerint, recombeulerint
use bdf_mod, only: bdfint
use atomic_rates_mod, only: get_atomic_rates
use physical_constants_mod
use global_mod, only: GV, saved_gheads, rtable
implicit none

private
public :: update_raylist
public :: update_no_hits
 
contains


!> sets the logicals that can be determined from global variables
!!-----------------------------------------------------------------
subroutine set_bools( He, caseA, isoT, fixT )
  logical, intent(out) :: He
  logical, intent(out) :: caseA(2)
  logical, intent(out) :: isoT
  logical, intent(out) :: fixT

#ifdef incHe
  He = .true.
#else
  He = .false.
#endif

  caseA = .false.
  if (.not. GV%OnTheSpotH  .or. GV%HydrogenCaseA) caseA(1) = .true.
  if (.not. GV%OnTheSpotHe .or. GV%HeliumCaseA)   caseA(2) = .true.

  if (GV%IsoTemp > 0.0) then
     isoT = .true.
  else
     isoT = .false.
  end if

  if (GV%FixSnapTemp) then
     fixT = .true.
  else
     fixT = .false.
  end if

end subroutine set_bools



!> updates the particles intersected by a ray 
!!-----------------------------------------------------------------
subroutine update_raylist(raylist, pars, box, srcray)

  type(raylist_type), intent(inout) :: raylist !< ray/particle intersections
  type(particle_type), intent(inout) :: pars(:)  !< particle system
  type(box_type), intent(in) :: box  !< particle system
  logical, intent(in) :: srcray !< is this update for a source ray?

  type(particle_type) :: par
  type(ionpart_type) :: ipar
  integer(i8b) :: impact  
  integer(i8b) :: scalls  ! number of calls to solver
  logical :: photo
  logical :: He
  logical :: caseA(2)
  logical :: isoT
  logical :: fixT
  integer(i8b) :: index


  ! set booleans
  !-------------------------------------------------------
  call set_bools( He, caseA, isoT, fixT )
  photo = .true.

  
  ! loop through the ray particle intersections
  !-------------------------------------------------------
  impact_loop: do impact = 1,raylist%nnb

     GV%ParticleCrossings = GV%ParticleCrossings + 1     
     index = raylist%intersection(impact)%pindx
     par = pars(index)


     

     ! check we dont have double intersections when we shouldn't
     !-------------------------------------------------------------
     if (srcray) then
        if (box%tbound(1)==1 .and. GV%itime == par%lasthit) then
           ! here we have periodic BCs and a particle has been hit
           ! twice by the same ray so we stop tracing 
           GV%PhotonsLeavingBox = GV%PhotonsLeavingBox + raylist%ray%pcnt
           raylist%lastnnb = impact-1
           exit
        else if (box%tbound(1)==0 .and. GV%itime == par%lasthit) then
           ! here we have transmissive BCs and a particle has been hit
           ! twice by the same ray so something is wrong
           write(*,*) "transmissive BCs and particle hit twice in one ray!"
           write(*,*) "impact  = ", impact
           ipar%strtag = "check boundary conditions"
           call ionpar2screen(ipar)
           stop
        end if
     end if

     call initialize_ionpar(ipar,par,index,srcray,He,raylist,impact)


!     write(*,*) "d,dl:", raylist%intersection(impact)%d, ipar%dl
!     write(*,"(A,4F12.6)") "pos,nH: ", ipar%pos, ipar%nH
!     write(*,*) "inside: ", ipar%inside
!     write(*,*) 

     if (srcray) then
        if (GV%IonTempSolver==1) then
           call eulerint(ipar,scalls,photo,caseA,He,isoT,fixT)
           ipar%strtag = "on_eulerint_output"
        else if (GV%IonTempSolver==2) then
           call bdfint(ipar,scalls,photo,caseA,He,isoT,fixT)
           ipar%strtag = "on_bdfint_output"
        end if
     else
        call recombeulerint(ipar,scalls)
        ipar%strtag = "on_recombeulerint_output"
     end if
     call check_x(ipar)

     raylist%lastnnb = impact

     GV%TotalDerivativeCalls = GV%TotalDerivativeCalls + scalls
     if (scalls .GT. GV%PeakUpdates) GV%PeakUpdates = scalls


     !  put the updated particle data into the particle system
     !===========================================================
     call ionpar2par(ipar,par)
     if (par%T < GV%Tfloor) par%T = GV%Tfloor

#ifdef outGammaHI
     par%gammaHI = pars(ipar%index)%gammaHI + ipar%gammaHI * ipar%dt_s
     par%time    = pars(ipar%index)%time + ipar%dt_s
#endif

     pars(ipar%index) = par

     if (srcray) then
        pars(ipar%index)%lasthit = GV%itime 
     end if


     !  use the solution to set some global variables
     !=================================================
     GV%TotalPhotonsAbsorbed = GV%TotalPhotonsAbsorbed + ipar%pdeps
     GV%TotalIonizations     = GV%TotalIonizations + &
                               (ipar%xHII - ipar%xHII_in) * ipar%Hcnt
#ifdef incHe
     GV%TotalIonizations     = GV%TotalIonizations + &
                               (ipar%xHeII - ipar%xHeII_in) * ipar%Hecnt     
     GV%TotalIonizations     = GV%TotalIonizations + &
                               (ipar%xHeIII - ipar%xHeIII_in) * ipar%Hecnt     
#endif


     ! if the particle satisfies the rec ray tol put it on the recomb list
     !=====================================================================
#ifdef incHrec
     if (srcray) then
        if (.not. pars(ipar%indx)%OnRecList) then
           if (pars(ipar%indx)%xHIIrc > GV%RecRayTol) then
              pars(ipar%indx)%OnRecList = .true.
              GV%recpt = GV%recpt + 1
              reclist(GV%recpt) = ipar%indx
           end if
        end if
     end if
#endif
     
     
     
     !  determine if we move to the next particle and track some ray stats
     !=====================================================================

     if (GV%RayDepletion) then
        raylist%ray%pcnt = raylist%ray%pcnt - ipar%pdeps
     endif
     
     
     ! if photons are exhausted
     !-------------------------------
     if (raylist%ray%pini > 0.0) then
        if (raylist%ray%pcnt / raylist%ray%pini < GV%RayPhotonTol) then
           exit
        end if
     end if
     
     ! if vacuum BCs and exiting box
     !-------------------------------
     if(box%tbound(1)==0) then
        if(impact==raylist%nnb) then
           GV%PhotonsLeavingBox = GV%PhotonsLeavingBox + raylist%ray%pcnt
        end if
     end if

     
  end do impact_loop

!  stop "finished impact loop"

end subroutine update_raylist


!==============================================================================
!> this routine is called to update particles that have not been hit by rays.
!! it reverse ray traces the non-hit particles.  in other words it traces rays
!! from those particles to the walls (the background sources) and updates the 
!! particle using a sum of those fluxes.  might add point sources to this
!! update later. 

subroutine update_no_hits(psys, tree)

  integer, parameter :: MaxSteps = 50000000

  integer, parameter :: ray_dirs(3,6) = reshape( (/  1, 0, 0, &
                                                     0, 1, 0, &
                                                     0, 0, 1, &
                                                    -1, 0, 0, &
                                                     0,-1, 0, &
                                                     0, 0,-1 /), (/3, 6/) )
  
  type(particle_system_type), intent(inout) :: psys !< particle system
  type(oct_tree_type), intent(in) :: tree !< oct-tree to search  

  type(particle_type) :: par
  type(ionpart_type) :: ipar  
  type(ray_type) :: ray(6)
  type(raylist_type) :: raylist(6)
  real(r8b) :: tauHI(6)
  real(r8b) :: parflux(6)
  

  integer(i8b) :: ip              ! particle counter.
  integer(i8b) :: rays_since_hit  ! rays traced since this particle has been hit.
  integer(i8b) :: scalls          ! number of times solver was called
  logical :: photo 
  logical :: srcray
  logical :: He
  logical :: caseA(2)
  logical :: isoT
  logical :: fixT

  integer :: i, j, isrc
  integer :: nsrcs

  real :: Hmf
  real(r8b) :: dir(3)
  real(r8b) :: pos(3)
  real(r8b) :: delta
  real(r8b) :: sigmaHI
  real(r8b) :: nHI
  real(r8b) :: wallflux
  real(r8b) :: freq
  real(r8b) :: enrg

  integer :: pindx_here
  integer :: pindx_last
  integer :: itravel

  type(gadget_constants_type) :: gconst
 

  photo=.true.
  

  ! find background source
  !--------------------------------------------
  do i = 1, size(psys%src)
     if ( psys%src(i)%EmisPrf == -3 ) isrc = i
     exit
  end do  
  freq = rn2freq( psys%src(isrc)%SpcType )
  enrg = freq * HI_th_erg
  sigmaHI = Verner_HI_photo_cs(freq)    


  wallflux = psys%src(isrc)%L * GV%Lunit 
      
 

  ! loop through all the particles
  !--------------------------------------------
  do ip = 1, size(psys%par)

     rays_since_hit = GV%rayn - psys%par(ip)%lasthit


     ! dont update if particle has been hit by a ray recently
     !--------------------------------------------------------
     if (rays_since_hit < GV%NraysUpdateNoHits) then

        cycle

     ! find the flux from each wall at this particle
     !--------------------------------------------------------        
     else

        par = psys%par(ip)
        pos = par%pos

        tauHI(:)   = 0.0d0
        parflux(:) = 0.0d0

        do i = 1,6

           dir = ray_dirs(:,i)
           call make_probe_ray( pos, dir, ray(i) )
           call trace_ray(ray(i), raylist(i), psys, tree) 

           itravel = sum ( maxloc( abs( dir ) ) ) ! x=1, y=2, z=3
           
           ! loop over particles in ray from source particle to wall
           ! first in list will be source particle
           !---------------------------------------------------------
           do j = 2,raylist(i)%nnb  

              pindx_here = raylist(i)%intersection(j)%pindx
              pindx_last = raylist(i)%intersection(j-1)%pindx

#ifdef incHmf
              Hmf = psys%par(pindx_here)%Hmf
#else
              Hmf = GV%H_mf
#endif

              nHI = psys%par(pindx_here)%rho * GV%cgs_rho * &
                   Hmf / gconst%PROTONMASS * psys%par(pindx_here)%xHI

              delta = abs( psys%par(pindx_here)%pos(itravel) - &
                           psys%par(pindx_last)%pos(itravel)   )

              delta = delta * GV%cgs_len
              
              tauHI(i) = tauHI(i) + nHI * delta * sigmaHI

           end do

           parflux(i) = wallflux * exp(-tauHI(i))           

        end do


        ! calculate He, caseA, isoT, and fixT
        !--------------------------------------------
        call set_bools( He, caseA, isoT, fixT )

        srcray = .false.
        call initialize_ionpar(ipar, par, ip, srcray, He)

        ! we need to initialize the variables that are 
        ! usually initialized from the raylist here.
        !--------------------------------------------

        ipar%bnorm = 0.0d0   
        ipar%cdfac = b2cdfac( ipar%bnorm, ipar%hsml, GV%cgs_len )

        ipar%dt_code = GV%dt_code * GV%NraysUpdateNoHits
        ipar%dt_s    = GV%dt_s    * GV%NraysUpdateNoHits

        ipar%pdeps    = 0.0d0
        ipar%pdeps_eq = 0.0d0

        ipar%penrg = enrg

        ipar%sigmaHI = Verner_HI_photo_cs(freq)    
        if (He) then
           ipar%sigmaHeI = Osterbrok_HeI_photo_cs(freq * HI_th_Hz)    
           ipar%sigmaHeII = Osterbrok_HeII_photo_cs(freq * HI_th_Hz)
        else
           ipar%sigmaHeI = 0.0d0
           ipar%sigmaHeII = 0.0d0
        end if

        ipar%pflux = sum(parflux)

        ! override the densities set in initialize_ionpar
        !---------------------------------------------------
        ipar%nH = ipar%rho * GV%cgs_rho * ipar%H_mf / gconst%PROTONMASS
        if (He) then
           ipar%nHe = ipar%rho * GV%cgs_rho * ipar%He_mf / (4*gconst%PROTONMASS)
        endif

        ! call solver
        !---------------------------------------------------
        call eulerint(ipar,scalls,photo,caseA,He,isoT,fixT)
        ipar%strtag = "on_reverse_trace_eulerint_output"


        GV%TotalDerivativeCalls = GV%TotalDerivativeCalls + scalls
        if (scalls .GT. GV%PeakUpdates) GV%PeakUpdates = scalls


        !  put the updated particle data into the particle system
        !===========================================================
        call ionpar2par(ipar,par)
        if (par%T < GV%Tfloor  ) par%T = GV%Tfloor
        if (par%T > GV%Tceiling) par%T = GV%Tceiling

#ifdef outGammaHI
        par%gammaHI = psys%par(ip)%gammaHI + ipar%gammaHI * ipar%dt_s
        par%time    = psys%par(ip)%time    + ipar%dt_s
#endif

        psys%par(ip) = par
        psys%par(ip)%lasthit = GV%itime



        
     end if



     
  end do
  
end subroutine update_no_hits







end module ion_temperature_update
