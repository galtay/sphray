!> \file ion_temperature_update.f90

!> \brief Module that is the starting point for ionization
!! and temperature updates.  It calls a runge kutta or a 
!! backwards difference formula (depending on the Makefile
!! flags)  
!<
module ion_temperature_update
use myf90_mod
use global_mod, only: global_variables_type
use raylist_mod, only: raylist_type
use ray_mod, only: ray_type, make_recomb_ray
use particle_system_mod, only: particle_type, box_type
use raylist_mod, only: trace_ray
use ionpar_mod
use euler_mod, only: eulerint, recombeulerint
use bdf_mod, only: bdfint
use physical_constants_mod
implicit none

private
public :: update_raylist, non_photo_update_all
 
contains

!> updates the particles intersected by a ray using 
!! one of the photo solvers
subroutine update_raylist(GV,raylist,pars,box,srcray)

  type(global_variables_type), intent(inout) :: GV !< global variables
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


#ifdef incHe
  He = .true.
#else
  He = .false.
#endif


  caseA = .false.
  if (.not. GV%OnTheSpotH  .or. GV%HydrogenCaseA) caseA(1) = .true.
  if (.not. GV%OnTheSpotHe .or. GV%HeliumCaseA)   caseA(2) = .true.

  if (GV%IsoTemp /= 0.0) then
     isoT = .true.
  else
     isoT = .false.
  end if

  if (GV%FixSnapTemp) then
     fixT = .true.
  else
     fixT = .false.
  end if


  ! loop through the ray particle intersections
  impact_loop: do impact = 1,raylist%nnb
     
     ipar%indx = raylist%intersection(impact)%pindx
     par = pars(ipar%indx)

     if (srcray) then
        if (box%tbound(1)==1 .and. GV%rayn.EQ.par%lasthit) then
           ! here we have periodic BCs and a particle has been hit
           ! twice by the same ray so we stop tracing 
           GV%PhotonsLeavingBox = GV%PhotonsLeavingBox + raylist%ray%pcnt
           raylist%lastnnb = impact-1
           exit
        else if (box%tbound(1)==0 .and. GV%rayn.EQ.par%lasthit) then
           ! here we have transmissive BCs and a particle has been hit
           ! twice by the same ray so something is wrong
           write(*,*) "transmissive BCs and particle hit twice in one ray!"
           write(*,*) "impact  = ", impact
           ipar%strtag = "check boundary conditions"
           call ionpar2screen(ipar)
           stop
        end if
     end if

     photo = .true.
     GV%ParticleCrossings = GV%ParticleCrossings + 1
     call initialize_ionpar(ipar,par,GV,srcray,He,raylist,impact)

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

     raylist%lastnnb = impact

     GV%TotalDerivativeCalls = GV%TotalDerivativeCalls + scalls
     if (scalls .GT. GV%PeakUpdates) GV%PeakUpdates = scalls

     ipar%strtag = "just exited eulerint"
     call check_x(ipar)


!    put the updated particle data into the particle system
!===========================================================
     call ionpar2par(ipar,par)
     if (par%T < GV%Tfloor) par%T = GV%Tfloor

#ifdef outGamma
     par%gammaHI = pars(ipar%indx)%gammaHI + ipar%gammaHI * ipar%dt_s
     par%time    = pars(ipar%indx)%time + ipar%dt_s
#endif

     pars(ipar%indx) = par

     if (srcray) then
        pars(ipar%indx)%lasthit = GV%rayn 
     end if


!  use the solution to set some global variables
!=================================================
     GV%TotalPhotonsAbsorbed = GV%TotalPhotonsAbsorbed + ipar%pdeps
     GV%TotalIonizations = GV%TotalIonizations + &
                           (ipar%xHII - ipar%xHII_in) * ipar%Hcnt
#ifdef incHe
     GV%TotalIonizations = GV%TotalIonizations + &
                           (ipar%xHeII - ipar%xHeII_in) * ipar%Hecnt     
     GV%TotalIonizations = GV%TotalIonizations + &
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

     raylist%ray%pcnt = raylist%ray%pcnt - ipar%pdeps
 
     ! if photons are exhausted
     if (raylist%ray%pini > 0.0) then
        if (raylist%ray%pcnt/raylist%ray%pini < GV%RayPhotonTol) then
           exit
        end if
     end if

     ! if vacuum BCs and exiting box
     if(box%tbound(1)==0) then
        if(impact==raylist%nnb) then
           GV%PhotonsLeavingBox = GV%PhotonsLeavingBox + raylist%ray%pcnt
        end if
     end if

     
  end do impact_loop


end subroutine update_raylist


!==============================================================================
!> this routine is called to update particles that have not been hit by rays.
!! it loops through every particle and does a non-photo update on all those
!! that havent been hit by a ray in the last GV%NonPhotoUpdateAllRays rays
!! traced.
subroutine non_photo_update_all(GV,pars)
  
  type(global_variables_type), intent(inout) :: GV !< global variables
  type(particle_type), intent(inout) :: pars(:) !< particle system
  
  type(particle_type) :: par
  type(ionpart_type) :: ipar  
  
  integer(i8b) :: p               ! particle counter.
  integer(i8b) :: rays_since_hit  ! rays traced since this particle has been hit.
  integer(i8b) :: scalls          ! number of times solver was called
  logical :: photo 
  logical :: srcray
  logical :: He
  logical :: caseA(2)
  logical :: isoT
  logical :: fixT

#ifdef incHe
  He = .true.
#else
  He = .false.
#endif

  caseA = .false.
  if (.not. GV%OnTheSpotH  .or. GV%HydrogenCaseA) caseA(1) = .true.
  if (.not. GV%OnTheSpotHe .or. GV%HeliumCaseA)   caseA(2) = .true.

  if (GV%IsoTemp /= 0.0) then
     isoT = .true.
  else
     isoT = .false.
  end if

  if (GV%FixSnapTemp) then
     fixT = .true.
  else
     fixT = .false.
  end if



  ipar%dt_code = GV%dtray_code * GV%UpdateAllRays
  ipar%dt_s = GV%dtray_s * GV%UpdateAllRays
  
  ! loop through all the particles
  do p = 1, size(pars)
     rays_since_hit = GV%rayn - pars(p)%lasthit
     
     ! dont update if particle has been hit by a ray recently
     if (rays_since_hit < GV%UpdateAllRays) then
        cycle
        
     ! update ionization and temperature
     else


        !  set the initial conditions of the system
        !  and call the appropriate solver
        !=============================================
        par = pars(p)
        srcray = .false.
        call initialize_ionpar(ipar,par,GV,srcray,He)
        photo = .false.
        call eulerint(ipar,scalls,photo,caseA,He,isoT,fixT)
        call ionpar2par(ipar,par)

        if (par%T < GV%Tfloor) par%T = GV%Tfloor

        
        pars(p) = par
        pars(p)%lasthit = GV%rayn
        
     end if
     
  end do
  
end subroutine non_photo_update_all







end module ion_temperature_update
