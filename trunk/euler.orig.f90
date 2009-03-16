!> \file euler.f90

!> \brief The module for the fancy euler solver.
!!
!<

module euler_mod
use atomic_rates_mod, only: atomic_rates_type
use ionpar_mod, only: ionpart_type
implicit none

  integer, parameter :: MaxSteps = 50000 !< maximum number of steps
  real, parameter :: ION_FAC = 0.001     !< ionization time factor
  real, parameter :: COOL_FAC = 0.001    !< cooling time factor

contains

!> this routine takes in the initial values in ipar 
!! {xHII,xHeII,xHeIII,T,pdep} and evolves the ionization particle
!! for a timestep ip%dt.  The final values are in the variables
!! ip%xHII, ip%xHeII, ip%xHeIII, ip%T, ip%pdeps
!========================================================================
subroutine eulerint(ip,scalls,photo)
use CenAtomicRates, only: Haiman_Bremss_cool
use CenAtomicRates, only: Haiman_Comp_Heol
use atomic_rates_mod, only: get_atomic_rates
use atomic_rates_mod, only: iso_k
use ionpar_mod, only: par2ionpar, ionpar2par, ionpar2screen
use ionpar_mod, only: set_gammas, set_cooling_func, set_dnedt
use physical_constants_mod

  type(ionpart_type), intent(inout) :: ip !< ionization particle
  integer, intent(out) :: scalls !< number of runs through the solver
  logical, intent(in) :: photo   !< consider photo ionization / heating ?
 
  type(atomic_rates_type) :: k
  real :: t_tot
  real :: dt2
  real :: dt_i
  real :: t_min
  integer :: i

  scalls = 0   
  t_tot = 0.0
  dt2 = ip%dt_s/2.


!---------------------------------------------------------
! if constant temperature run, set constant atomic rates
#ifndef incT
  k = iso_k
#endif

!---------------------------------------------------------
! do first round of calculations

  ip%nHII = ip%nH * ip%xHII
  ip%nHI  = ip%nH - ip%nHII
  ip%xHI  = ip%nHI / ip%nH
  ip%ne = ip%nHII + ip%NeBckgnd
#ifdef incHe
  ip%nHeIII = ip%nHe * ip%xHeIII
  ip%nHeII  = ip%nHe * ip%xHeII
  ip%nHeI   = ip%nHe - ip%nHeIII - ip%nHeII
  ip%xHeI  = ip%nHeI / ip%nHe
  ip%ne = ip%nHII + ip%nHeII + 2.0 * ip%nHeIII + ip%NeBckgnd
#endif  

  if (photo) then
     call set_gammas(ip)
  else
     ip%gammaHI = 0.0
     ip%gammaHeI = 0.0
     ip%gammaHeII = 0.0
  end if

!---------------------------------------------------------  
  
  do i = 1,MaxSteps

     ip%iter = i

#ifdef incT
     ip%u = 1.5 * (ip%nH + ip%nHe + ip%ne) * k_erg_K * ip%T 
     call get_atomic_rates(ip%T,k)
     call set_cooling_func(ip,k)
     ip%dudt = - (ip%CIC + ip%CEC + ip%RCCb + ip%BREM + ip%COMP)
     ip%dudt = ip%PH + ip%dudt
     ip%dTdt = 2.0 * ip%dudt / (3.0 * (ip%nH + ip%nHe + ip%ne) * k_erg_K)
#endif

     !----------------------------------------------------------  
     ! calculate the time rate of change of the electron density
     call set_dnedt(ip,k)

     !----------------------------------------------------------  
     ! calculate the time step
     ip%tion  = abs(ip%ne / ip%dnedt)
#ifdef incT
     ip%tcool = abs(ip%u / ip%dudt)
     dt_i = min(ION_FAC * ip%tion, COOL_FAC * ip%tcool)
#else
     dt_i = ION_FAC * ip%tion
#endif
     
     if (dt_i > ip%dt_s - t_tot) dt_i = ip%dt_s - t_tot
     if (dt_i > dt2) dt_i = dt2

     !----------------------------------------------------------  
     ! take the time step
     call take_euler_step(ip,dt_i)


     !----------------------------------------------------------  
     ! check bounds - this method is guaranteed to stay positive so 
     ! if any are greater than 1 use the small part to set x = 1-small

     !  if (ip%xHI < ip%xHII) then
     !     ip%xHII = 1.0 - ip%xHI
     !  else
     !     ip%xHI = 1.0 - ip%xHII
     !  end if


     !----------------------------------------------------------  
     ! track photons, recombinations and update total time

     ! photo ionizations

     if (photo) then
        ip%pdeps = ip%pdeps + ip%gammaHI * ip%HIcnt * dt_i
#ifdef incHe
        ip%pdeps = ip%pdeps + &
             (ip%gammaHeI * ip%HeIcnt + ip%gammaHeII * ip%HeIIcnt) * dt_i
#endif
     end if

     ! recombination fractions

#ifdef increc 
     ip%xHIIrc   = ip%xHIIrc   +  (k%HIIrcA - k%HIIrcB)     * ip%ne * dt_i

     ! the Helium has to be done better because it is not just direct 
     ! recombination to the ground level that is important here
#ifdef incHe   
     ip%xHeIIrc  = ip%xHeIIrc  +  (k%HeIIrcA - k%HeIIrcB)   * ip%ne * dt_i
     ip%xHeIIIrc = ip%xHeIIIrc +  (k%HeIIIrcA - k%HeIIIrcB) * ip%ne * dt_i
#endif
#endif

     ! time step

     t_tot = t_tot + dt_i
     t_min = huge(1.0)
     if (t_tot < t_min) t_min = t_tot
     if (abs(ip%dt_s - t_min) < 0.0001 * ip%dt_s) exit
    
     if (i == MaxSteps) then
        write(*,*) " ************* EULER ERROR ************* "
        write(*,*) "finished ", MaxSteps, " little steps before "
        write(*,*) "completing large timestep.  "
        write(*,*) "t/dt = ", t_tot/ip%dt_s, "rayn = ", ip%rayn
        write(*,*) "ipar = "
        call ionpar2screen(ip)
        write(*,*) " *************  euler.f90  *************"
        stop
     end if

     scalls = scalls + 1
     if(photo) call set_gammas(ip)

  end do

  ip%xHII = ip%nHII   / ip%nH
#ifdef incHe
  ip%xHeII = ip%nHeII  / ip%nHe
  ip%xHeIII = ip%nHeIII / ip%nHe
#endif


  
end subroutine eulerint
 

!> this routine takes in the initial values in ipar 
!! {xHII,xHeII,xHeIII,T,pdep} and simply deposits recombination 
!! photons based on the optical depth.  The final values are in the variables
!! ip%xHII, ip%xHeII, ip%xHeIII, ip%T, ip%pdeps
!========================================================================
subroutine recombeulerint(ip,scalls)
use CenAtomicRates, only: Haiman_Bremss_cool
use CenAtomicRates, only: Haiman_Comp_Heol
use atomic_rates_mod, only: get_atomic_rates
use atomic_rates_mod, only: iso_k
use ionpar_mod, only: par2ionpar, ionpar2par, ionpar2screen
use ionpar_mod, only: set_gammas, set_cooling_func, set_dnedt
use physical_constants_mod

  real, parameter :: TAU_TOL = 1.0e-2

  type(ionpart_type), intent(inout) :: ip !< ionization particle
  integer, intent(out) :: scalls !< number of runs through the solver
 
  scalls = 0   

!---------------------------------------------------------
! set initial number and ionization fractions

  ip%nHII = ip%nH * ip%xHII
  ip%nHI  = ip%nH - ip%nHII
  ip%xHI  = ip%nHI / ip%nH
  ip%ne = ip%nHII + ip%NeBckgnd
#ifdef incHe
  ip%nHeIII = ip%nHe * ip%xHeIII
  ip%nHeII  = ip%nHe * ip%xHeII
  ip%nHeI   = ip%nHe - ip%nHeIII - ip%nHeII
  ip%xHeI  = ip%nHeI / ip%nHe
  ip%ne = ip%nHII + ip%nHeII + 2.0 * ip%nHeIII + ip%NeBckgnd
#endif  

!---------------------------------------------------------
! simlple count of atoms that are absorbing photons

  ip%HIcnt = ip%Hcnt * ( ip%nHI / ip%nH )       ! N HI atoms
  ip%Allcnt = ip%HIcnt
#ifdef incHe
  ip%HeIcnt = ip%Hecnt * ( ip%nHeI / ip%nHe )   ! N HeI atoms
  ip%HeIIcnt = ip%Hecnt * ( ip%nHeII / ip%nHe ) ! N HeII atoms
  if (ip%penrg .GE. HeI_th_erg) then
     ip%Allcnt = ip%HIcnt + ip%HeIcnt
  else if (ip%penrg .GE. HeII_th_erg) then
     ip%Allcnt = ip%HIcnt + ip%HeIcnt + ip%HeIIcnt
  end if
#endif  

!---------------------------------------------------------
! optical depths

  ip%tauHI = ip%cdfac * ip%HIcnt * ip%sigmaHI
  ip%tausum = ip%tauHI
  if (ip%tauHI .LT. TAU_TOL) then
     ip%HItaufac = ip%tauHI
  else
     ip%HItaufac = 1.0 - exp(-ip%tauHI)
  end if
  ip%taufacsum = ip%HItaufac  

#ifdef incHe
  ip%tauHeI = ip%cdfac * ip%HeIcnt * ip%sigmaHeI
  ip%tauHeII = ip%cdfac * ip%HeIIcnt * ip%sigmaHeII
  ip%tausum = ip%tauHI + ip%tauHeI + ip%tauHeII 
  if (ip%tauHeI .LT. TAU_TOL) then
     ip%HeItaufac = ip%tauHeI
  else
     ip%HeItaufac = 1.0 - exp(-ip%tauHeI)
  end if
  if (ip%tauHeII .LT. TAU_TOL) then
     ip%HeIItaufac = ip%tauHeII
  else
     ip%HeIItaufac = 1.0 - exp(-ip%tauHeII)
  end if
  ip%taufacsum = ip%HItaufac + ip%HeItaufac + ip%HeIItaufac
#endif 

!---------------------------------------------------------
! deposit photons

  ip%fracabsorb = 1.0 - exp(-ip%tausum)

  ip%pdeps = ip%pflux * ip%fracabsorb

  ! if more photons going in then absorbers
  if (ip%pdeps >= ip%HIcnt) then
     ip%pdeps = ip%HIcnt
     ip%xHII = 1.0
     ip%xHI = 0.0
  else
     ip%xHII = ip%xHII + ip%pdeps / ip%Hcnt
     ip%xHI = 1.0 - ip%xHII
  end if

  scalls = scalls + 1
     
  
end subroutine recombeulerint


!> this routine evolves the ionization fractions and possibly the  
!! temperature of the ionization particle forward in time by dt
!================================================================
subroutine take_euler_step(ip,dt)
use physical_constants_mod, only: k_erg_K

  type(ionpart_type), intent(inout) :: ip !< ionization particle  
  real, intent(inout) :: dt               !< time step

  real :: det
  real :: RRHII
  real :: GGHI

#ifdef incHe
  real :: HeI,HeII,HeIII  
  real :: RRHeII, RRHeIII
  real :: GGHeI, GGHeII
#endif



#ifdef increc
  RRHII   = ip%RRHIIa
#else
  RRHII   = ip%RRHIIb
#endif


  GGHI   = ip%GGHI
  det = 1 + (GGHI + RRHII) * dt

  ! HI
  ip%xHI = (1 + dt*RRHII) * ip%xHI   + &
             (dt*RRHII)   * ip%xHII
  ip%xHI = ip%xHI/det
  

  ! HII
  ip%xHII = (dt*GGHI)      * ip%xHI  + & 
            (1 + dt*GGHI)  * ip%xHII
  ip%xHII = ip%xHII/det

  ip%nHI  = ip%nH * ip%xHI
  ip%nHII = ip%nH * ip%xHII
  ip%ne = ip%nHII + ip%NeBckgnd


! Helium Code
!-------------
#ifdef incHe

#ifdef increc
  RRHeII  = ip%RRHeIIa
  RRHeIII = ip%RRHeIIIa
#else
  RRHeII  = ip%RRHeIIb
  RRHeIII = ip%RRHeIIIb
#endif

  GGHeI  = ip%GGHeI
  GGHeII = ip%GGHeII
  
  det = 1.0 + (GGHeI + GGHeII + RRHeII + RRHeIII) * dt + &
            + (GGHeI * GGHeII + GGHeI * RRHeIII + RRHeII * RRHeIII) * (dt*dt)

  ! HeI
  HeI   = 1.0 + (GGHeII + RRHeII + RRHeIII + dt * RRHeII * RRHeIII)* dt
  HeII  = dt * RRHeII * (1.0 + dt * RRHeIII)
  HeIII = dt * dt * RRHeII * RRHeIII

  ip%xHeI = HeI * ip%xHeI + HeII * ip%xHeII + HeIII * ip%xHeIII
  ip%xHeI = ip%xHeI/det

  ! HeII
  HeI   = GGHeI * dt * (1.0 + RRHeIII * dt)
  HeII  = (1.0 + GGHeI * dt) * (1.0 + dt * RRHeIII)
  HeIII = dt * (1.0 + GGHeI * dt) * RRHeIII

  ip%xHeII = HeI * ip%xHeI + HeII * ip%xHeII + HeIII * ip%xHeIII
  ip%xHeII = ip%xHeII/det

  ! HeIII
  HeI   = GGHeI * GGHeII * dt * dt
  HeII  = GGHeII * dt * (1.0 + GGHeI * dt)
  HeIII = 1.0 + (GGHeI + GGHeII + dt * GGHeI * GGHeII + RRHeII ) * dt

  ip%xHeIII = HeI * ip%xHeI + HeII * ip%xHeII + HeIII * ip%xHeIII
  ip%xHeIII = ip%xHeIII/det

  ip%nHeI   = ip%nHe * ip%xHeI
  ip%nHeII  = ip%nHe * ip%xHeII
  ip%nHeIII = ip%nHe * ip%xHeIII
  ip%ne = ip%nHII + ip%nHeII + 2.0 * ip%nHeIII + ip%NeBckgnd

#endif
! End Helium code


#ifdef incT
  ip%T = ip%T + ip%dTdt * dt
  ip%u = 1.5d0 * (ip%nH + ip%nHe + ip%ne) * k_erg_K * ip%T
#endif
 
end subroutine take_euler_step


end module euler_mod
