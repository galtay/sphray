!> \file ionpar.f90

!> \brief the module that handles ionization particles
!<

module ionpar_mod
use particle_system_mod, only: particle_type
use global_mod, only: global_variables_type
use raylist_mod, only: raylist_type
use myf90
implicit none

!---------------------------
!> ionization particle type. 
type ionpart_type

   integer(i8b) :: id      !< particle id
   real    :: pos(3)       !< x,y,z coordinates
   real    :: vel(3)       !< x,y,z velocities
   real    :: hsml         !< smoothing length
   real    :: rho          !< density = mass * NsphNnb / hsml^3 
   real    :: mass         !< particle mass
   real    :: T            !< temperature in K       
   real    :: xHII         !< HII ionization fraction
   real    :: xHeII        !< HeII ionization fraction
   real    :: xHeIII       !< HeIII ionization fraction
   real    :: xHIIrc       !< HII recombination fraction
   real    :: xHeIIrc      !< HeII recombination fraction
   real    :: xHeIIIrc     !< HeIII recombination fraction
   integer(i8b) :: lasthit !< last ray to cross this particle

   integer :: rayn   !< ray number
   integer :: iter   !< number of the iteration in the solver
   integer :: impact !< index of the particle in the raylist
   integer :: indx   !< index of the particle in the psys
   real :: d         !< distance along ray (code)
   real :: b         !< impact parameter 

   real :: bnorm     !< impact parameter normalized to smoothing length
   real :: cdfac     !< column depth conversion factor
   real :: gpercm3   !< density in cgs
   real :: cm3       !< volume in cgs

   real :: dt_code  !< time step for this par (code)
   real :: dt_s     !< time step for this par (s)

   real :: pflux       !< total photon flux arriving at the particle
   real :: fracabsorb  !< fraction of all photons absorbed
   real :: pdepr       !< rate that photons are being deposited
   real :: penrg       !< energy of one photon in the ray
   real :: pdeps       !< total photons deposited into particle

!   real :: HIcolions   !< total HI collisional ionizations in the particle
!   real :: HeIcolions  !< total HI collisional ionizations in the particle
!   real :: HeIIcolions !< total HI collisional ionizations in the particle

!   real :: HIIrcmbsB   !< total HII recombinations excluding to n=1 lvl
!   real :: HIIrcmbsA   !< total HII recombinations to all levels
!   real :: HeIIrcmbsB  !< total HeII recombinations excluding to n=1 lvl
!   real :: HeIIrcmbsA  !< total HeII recombinations to all levels
!   real :: HeIIIrcmbsB !< total HeIII recombinations excluding to n=1 lvl
!   real :: HeIIIrcmbsA !< total HeIII recombinations to all levels


   real :: Tcmb     !< background radiation field temperature
   real :: H_mf     !< Hydrogen mass fraction
   real :: He_mf    !< Helium mass fraction

   real :: nH       !< number density of H
   real :: nHe      !< number density of He
   real :: Hcnt     !< number of H nuclei
   real :: Hecnt    !< number of He nuclei

   real :: xHI      !< nHI/nH

   real :: NeBckgnd !< number density of metallic electrons
   real :: nHI      !< number density of HI
   real :: nHII     !< number density of HII
   real :: ne       !< number density of electrons
   real :: dnedt    !< time rate of change of ne

   real :: HIcnt    !< number of HI atoms
   real :: HIIcnt   !< number of HII atoms

   real :: xHeI     !< nHeI/nHe

   real :: nHeI     !< number density of HeI
   real :: nHeII    !< number density of HeII
   real :: nHeIII   !< number density of HeIII

   real :: HeIcnt   !< number of HeI atoms
   real :: HeIIcnt  !< number of HeII atoms
   real :: HeIIIcnt !< number of HeIII atoms
   real :: Allcnt   !< number of all photo absorbing species (HI,HeI,HeII)

   real :: tauHI    !< HI optical depth
   real :: tauHeI   !< HeI optical depth
   real :: tauHeII  !< HeII optical depth
   real :: tausum   !< sum of all optical depths

   real :: HItaufac   !< 1-exp(-tauHI)
   real :: HeItaufac  !< 1-exp(-tauHeI)
   real :: HeIItaufac !< 1-exp(-tauHeII)
   real :: taufacsum  !< sum of all tau factors

   real :: HIfrac   !< fraction of photons absorbed by HI
   real :: HeIfrac  !< fraction of photons absorbed by HeI
   real :: HeIIfrac !< fraction of photons absorbed by HeII

   real :: sigmaHI   !< HI photo cross-section
   real :: sigmaHeI  !< HeI photo cross-section
   real :: sigmaHeII !< HeII photo cross-section

   real :: gammaHI   !< HI photoionization rate
   real :: gammaHeI  !< HeI photoionization rate
   real :: gammaHeII !< HeII photoionization rate
   real :: gammasum  !< sum of all photoionization rates

   real :: GGHI     !< Photoionization rate + CI * ne
   real :: GGHeI    !< Photoionization rate + CI * ne
   real :: GGHeII   !< Photoionization rate + CI * ne
   real :: RRHIIa   !< Recombination rate A * ne
   real :: RRHeIIa  !< Recombination rate A * ne
   real :: RRHeIIIa !< Recombination rate A * ne
   real :: RRHIIb   !< Recombination rate B * ne
   real :: RRHeIIb  !< Recombination rate B * ne
   real :: RRHeIIIb !< Recombination rate B * ne

   real :: CIC  !< collisional ionization cooling rate
   real :: CEC  !< collisional excitation cooling rate
   real :: RCCa !< recombination cooling case A
   real :: RCCb !< recombination cooling case B
   real :: PH   !< photo heating rate
   real :: BREM !< bremsstrahlung cooling rate
   real :: COMP !< compton heating/cooling rate

   real :: tion     !< ionization time (s)
   real :: tcool    !< cooling time (s)
   real :: u        !< energy of particle (ergs)
   real :: dudt     !< time rate of change of energy
   real :: dTdt     !< time rate of change of temperature

   real :: xHI_in      !< initial xHI
   real :: xHII_in     !< initial xHII
   real :: xHeI_in     !< initial xHeI
   real :: xHeII_in    !< initial xHeII
   real :: xHeIII_in   !< initial xHeIII
   real :: T_in        !< initial T
   
   real :: pdeps_eq    !< photons deposited in equilibrium conditions

   character(200) :: deriv_call !< labels call to derivative in case of crash

end type ionpart_type

contains

!> initializes the ionization particle values
subroutine initialize_ionpar(ipar,par,GV,srcray,raylist,impact)
use b2cd_mod, only: b2cdfac
use physical_constants_mod, only: HI_th_nu, M_H, M_He
use CenAtomicRates, only: Osterbrok_HI_photo_cs
use CenAtomicRates, only: Osterbrok_HeI_photo_cs
use CenAtomicRates, only: Osterbrok_HeII_photo_cs

  type(ionpart_type), intent(inout) :: ipar           !< ionization particle
  type(particle_type), intent(in) :: par              !< standard particle
  type(global_variables_type), intent(in) :: GV       !< global variables
  logical, intent(in) :: srcray                       !< source ray update ?
  type(raylist_type), intent(in), optional :: raylist !< optional raylist
  integer, intent(in), optional :: impact             !< optional impact number

  call par2ionpar(GV%IsoMass,GV%IsoTemp,par,ipar)

  ipar%rayn = GV%rayn
  if (present(impact)) ipar%impact = impact
  ipar%NeBckgnd = GV%NeBackground
  ipar%H_mf = GV%H_mf
  ipar%He_mf = GV%He_mf
  ipar%Tcmb = GV%Tcmb_cur

  ipar%gpercm3 = ipar%rho * GV%cgs_mass / &
                (GV%cgs_len * GV%cgs_len * GV%cgs_len)

  ipar%cm3 = ipar%mass * GV%cgs_mass / ipar%gpercm3
  ipar%nH = ipar%gpercm3 * ipar%H_mf  / M_H
  ipar%Hcnt = ipar%mass * ipar%H_mf * GV%cgs_mass  / M_H  

#ifdef incHe
  ipar%nHe = ipar%gpercm3 * ipar%He_mf / M_He
  ipar%Hecnt = ipar%mass * ipar%He_mf * GV%cgs_mass / M_He
#else
  ipar%nHe = 0.0d0
  ipar%Hecnt = 0.0d0
#endif


  if (present(raylist)) then

     ipar%d = raylist%intersection(impact)%d
     ipar%b = raylist%intersection(impact)%b
     ipar%bnorm = ipar%b/ipar%hsml
     ipar%cdfac = b2cdfac(ipar%bnorm,ipar%hsml,GV%cgs_len)


     ipar%dt_code = (GV%rayn - ipar%lasthit) * GV%dtray_code
     ipar%dt_s    = (GV%rayn - ipar%lasthit) * GV%dtray_s
     
     if (srcray) then        
        ipar%pflux = raylist%ray%pcnt / ipar%dt_s 
     else
        ipar%pflux = raylist%ray%pcnt 
     end if

     ipar%penrg = raylist%ray%enrg
     ipar%pdeps = 0.0
     ipar%pdeps_eq = 0.0

     ipar%sigmaHI = Osterbrok_HI_photo_cs(raylist%ray%freq * HI_th_nu)    
#ifdef incHe
     ipar%sigmaHeI = Osterbrok_HeI_photo_cs(raylist%ray%freq * HI_th_nu)    
     ipar%sigmaHeII = Osterbrok_HeII_photo_cs(raylist%ray%freq * HI_th_nu)
#else
     ipar%sigmaHeI = 0.0
     ipar%sigmaHeII = 0.0
#endif
     
  end if
     
end subroutine initialize_ionpar

!-----------------------------------------------------------------------
!> copies the basic particle data into an ionization particle
subroutine par2ionpar(dfltMass,dfltTemp,par,ipar)
 real, intent(in) :: dfltMass               !< default mass for isomass
 real, intent(in) :: dfltTemp               !< default temperature for isotemp
 type(particle_type), intent(in) :: par     !< input particle
 type(ionpart_type), intent(inout) :: ipar  !< output ionization particle

 real :: mass
 real :: temp

 mass = dfltMass
#ifdef incmass
 mass = par%mass
#endif

 temp = dfltTemp
#ifdef incT
 temp = par%T
#endif

 ipar%pos = par%pos
 ipar%hsml = par%hsml
 ipar%rho = par%rho
 ipar%xHII = par%xHII
 ipar%lasthit = par%lasthit

!------------
#ifdef incid
 ipar%id = par%id
#else
 ipar%id = 0
#endif

!------------
#ifdef incvel
 ipar%vel = par%vel
#else
 ipar%vel = 0.0d0
#endif

!------------
 ipar%mass = mass

!------------
 ipar%T = temp 

!------------
#ifdef incHe
 ipar%xHeII  = par%xHeII
 ipar%xHeIII = par%xHeIII
#else
 ipar%xHeII  = 0.0
 ipar%xHeIII = 0.0
#endif

!------------
#ifdef increc
 ipar%xHIIrc = par%xHIIrc
#endif

#ifdef incHe
#ifdef increc
 ipar%xHeIIrc = par%xHeIIrc
 ipar%xHeIIIrc = par%xHeIIIrc
#endif
#endif


 
 ipar%xHI_in = 1.0 - par%xHII
 ipar%xHII_in = par%xHII
#ifdef incHe
 ipar%xHeI_in = 1.0 - par%xHeII - par%xHeIII
 ipar%xHeII_in = par%xHeII
 ipar%xHeIII_in = par%xHeIII
#else
 ipar%xHeI_in = 0.0
 ipar%xHeII_in = 0.0
 ipar%xHeIII_in = 0.0
#endif

#ifdef incT
 ipar%T_in = par%T
#else
 ipar%T_in = dfltTemp
#endif


end subroutine par2ionpar

!-----------------------------------------------------------------------
!> copies the ionization particle data into a basic particle 
subroutine ionpar2par(ipar,par)
 type(ionpart_type), intent(in) :: ipar    !< input ionization particle
 type(particle_type), intent(inout) :: par !< output particle

#ifdef incT
 par%T = ipar%T
#endif

 par%xHII = ipar%xHII

#ifdef incHe
 par%xHeII  = ipar%xHeII
 par%xHeIII = ipar%xHeIII
#endif

#ifdef increc
 par%xHIIrc = ipar%xHIIrc
#endif
#ifdef incHe
#ifdef increc
 par%xHeIIrc = ipar%xHeIIrc
 par%xHeIIIrc = ipar%xHeIIIrc
#endif
#endif

 par%lasthit = ipar%lasthit

end subroutine ionpar2par



!-----------------------------------------------------
!> prints ionization particle information to screen
subroutine ionpar2screen(ipar)

  type(ionpart_type), intent(in) :: ipar  !< ionization particle
  
   95 format(A,T25,I15)
   100 format(A,T25,3F12.4)
   105 format(A,T25,1F12.4)
   110 format(A,T25,3ES18.9)
   115 format(A,T25,2ES18.9)
   write(*,*) 
   write(*,95) "ID", ipar%id
   write(*,100) "pos", ipar%pos
   write(*,100) "vel", ipar%vel
   write(*,105) "hsml", ipar%hsml
   write(*,105) "rho", ipar%rho
   write(*,105) "mass", ipar%mass
   write(*,105) "T", ipar%T
   write(*,95) "lasthit", ipar%lasthit
   write(*,*) 
   write(*,95) "ray num", ipar%rayn
   write(*,95) "psys indx", ipar%indx
   write(*,95) "impact", ipar%impact
   write(*,95) "iteration", ipar%iter
   write(*,105) "raydist", ipar%d
   write(*,*) 
   write(*,*) "deriv call = ", trim(ipar%deriv_call)
   write(*,105) "xHI+xHII", ipar%xHI + ipar%xHII
   write(*,'(A,T25,ES18.9)') "T_in", ipar%T_in
   write(*,110) "b, hsml, bnorm",ipar%b, ipar%hsml, ipar%bnorm
   write(*,*) 
   write(*,110) "xHII_in,HeII_in,HeIII_in", &
                 ipar%xHII_in, ipar%xHeII_in, ipar%xHeIII_in
   write(*,110) "xHI_in,HeI_in,HeII_in", &
                 ipar%xHI_in, ipar%xHeI_in, ipar%xHeII_in
   write(*,*) 
   write(*,110) "xHII,xHeII,xHeIII", ipar%xHII, ipar%xHeII, ipar%xHeIII
   write(*,110) "nHII,nHeII,nHeIII", ipar%nHII, ipar%nHeII, ipar%nHeIII
   write(*,*) 
   write(*,110) "xHI,xHeI,xHeII", ipar%xHI, ipar%xHeI, ipar%xHeII
   write(*,110) "nHI,nHeI,nHeII", ipar%nHI, ipar%nHeI, ipar%nHeII
   write(*,*) 
   write(*,110) "HIIrc,HeIIrc,HeIIIrc", & 
                 ipar%xHIIrc, ipar%xHeIIrc, ipar%xHeIIIrc
   write(*,'(A,T25,ES18.9)') "RCCa", ipar%RCCa
   write(*,*) 
   write(*,110) "tauHI,HeI,HeII", ipar%tauHI, ipar%tauHeI, ipar%tauHeII
   write(*,110) "HI,HeI,HeIItaufac", &
                 ipar%HItaufac, ipar%HeItaufac, ipar%HeIItaufac
   write(*,*) 
   write(*,110) "nH,nHe,ne", ipar%nH, ipar%nHe, ipar%ne
   write(*,110) "Hcnt,Hecnt,dnedt", ipar%Hcnt, ipar%Hecnt, ipar%dnedt
   write(*,*) 
   write(*,110) "HIcnt,HeIcnt,HeIIcnt", ipar%HIcnt, ipar%HeIcnt, ipar%HeIIcnt 
   write(*,110) "HI,HeI,HeIIfrac", ipar%HIfrac,ipar%HeIfrac,ipar%HeIIfrac
   write(*,*)
   write(*,110) "sigmaHI,HeI,HeII", ipar%sigmaHI, ipar%sigmaHeI, ipar%sigmaHeII
   write(*,110) "GammaHI,HeI,HeII", ipar%gammaHI, ipar%gammaHeI, ipar%gammaHeII
   write(*,*)
   write(*,110) "GGHI,HeI,HeII", ipar%GGHI, ipar%GGHeI, ipar%GGHeII
   write(*,110) "RRHII,HeII,HeIII", ipar%RRHIIb, ipar%RRHeIIb, ipar%RRHeIIIb
   write(*,*) 
   write(*,110) "CIC,CEC,RCCb", ipar%CIC, ipar%CEC, ipar%RCCb
   write(*,110) "PH,BREM,COMP", ipar%PH, ipar%BREM, ipar%COMP
   write(*,*) 
   write(*,110) "pdepr,penrg,pdeps", ipar%pdepr, ipar%penrg, ipar%pdeps
!   write(*,110) "colions HI,HeI,HeII", & 
!                 ipar%HIcolions,ipar%HeIcolions,ipar%HeIIcolions
!   write(*,*)
!   write(*,110) "rcmbsA HII,HeII,HeIII", &
!               ipar%HIIrcmbsA,ipar%HeIIrcmbsA,ipar%HeIIIrcmbsA
!   write(*,110) "rcmbsB HII,HeII,HeIII", &
!               ipar%HIIrcmbsB,ipar%HeIIrcmbsB,ipar%HeIIIrcmbsB
   write(*,*)
   write(*,115) "cdfac,d", ipar%cdfac,  ipar%d
   write(*,*) 
   write(*,115) "gpercm3,cm3", ipar%gpercm3, ipar%cm3
   write(*,115) "pflux, fracabsorb", ipar%pflux, ipar%fracabsorb
   write(*,*) 
   write(*,115) "u,dudt", ipar%u, ipar%dudt 
   write(*,115) "t,dTdt", ipar%T, ipar%dTdt
   write(*,*) 
   write(*,115) "Allcnt,GammaSum", ipar%Allcnt,ipar%gammasum
   write(*,115) "tausum,taufacsum",ipar%tausum, ipar%taufacsum
   write(*,*) 
   write(*,115) "tion,tcool", ipar%tion, ipar%tcool
   write(*,115) "dt_code, dt_s", ipar%dt_code, ipar%dt_s

   write(*,'(A,T25,ES18.9)') "pdeps_eq", ipar%pdeps_eq
   

 end subroutine ionpar2screen


!================================================================
!> dummy checks the bounds on HII, HeII, HeIII, and T
subroutine check_x(ip) 

  real, parameter :: TOL = 1.0e-6
  type(ionpart_type), intent(inout) :: ip  !< ionization particle
  logical :: bad

  100 format(A,I2,A)
  101 format(A,I2,A,F15.8)
  102 format(A,I15)
  103 format(A,2F15.8)

  bad = .false.
 
  if ( ip%xHI .LT. 0.0 - TOL ) then
     bad = .true.
     write(*,100) "xHI < 0.0 in deriv"
  end if
     
  if ( ip%xHI .GT. 1.0 + TOL ) then
     bad = .true.
     write(*,100) "xHI > 1.0 in deriv"
  end if

  if ( ip%xHII .LT. 0.0 - TOL ) then
     bad = .true.
     write(*,100) "xHII < 0.0 in deriv"
  end if
     
  if ( ip%xHII .GT. 1.0 + TOL ) then
     bad = .true.
     write(*,100) "xHII > 1.0 in deriv"
  end if
     

  if ( ip%xHeI .LT. 0.0 - TOL ) then
     bad = .true.
     write(*,100) "xHeI < 0.0 in deriv"
  end if
     
  if ( ip%xHeI .GT. 1.0 + TOL ) then
     bad = .true.
     write(*,100) "xHeI > 1.0 in deriv"
  end if

  if ( ip%xHeII .LT. 0.0 - TOL ) then
     bad = .true.
     write(*,100) "xHeII < 0.0 in deriv"
  end if
     
  if ( ip%xHeII .GT. 1.0 + TOL ) then
     bad = .true.
     write(*,100) "xHeII > 1.0 in deriv"
  end if

  if ( ip%xHeIII .LT. 0.0 - TOL ) then
     bad = .true.
     write(*,100) "xHeIII < 0.0 in deriv"
  end if
     
  if ( ip%xHeIII .GT. 1.0 + TOL ) then
     bad = .true.
     write(*,100) "xHeIII > 1.0 in deriv"
  end if

#ifdef incT
  if ( ip%T .LE. 0.0 ) then
     bad = .true.
     write(*,*) "T <= 0.0 in rk_deriv"
  end if
#endif

  if (bad) then
     call ionpar2screen(ip)
     stop
  end if

end subroutine check_x

!> sets the photoionization rate for an ionization particle
!================================================================
subroutine set_gammas(ip)
use physical_constants_mod, only: HI_th_erg, HeI_th_erg, HeII_th_erg

  real, parameter :: TAU_TOL = 1.0e-2
  type(ionpart_type), intent(inout) :: ip !< ionization particle  

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

  ! fraction absorbed and deposition rates
  ip%fracabsorb = 1.0 - exp(-ip%tausum)
  ip%pdepr = ip%pflux * ip%fracabsorb
    
  if( ip%taufacsum > 0.0 ) then
     ip%gammasum  = ip%pdepr / ip%Allcnt
     ip%HIfrac    = ip%HItaufac / ip%taufacsum
     ip%gammaHI   = ip%gammasum * ip%HIfrac 
#ifdef incHe
     ip%HeIfrac   = ip%HeItaufac  / ip%taufacsum
     ip%HeIIfrac  = ip%HeIItaufac / ip%taufacsum
     ip%gammaHeI  = ip%gammasum * ip%HeIfrac 
     ip%gammaHeII = ip%gammasum * ip%HeIIfrac 
#endif      

  else 
     ip%HIfrac = 0.0
     ip%gammaHI = 0.0
#ifdef incHe
     ip%HeIfrac = 0.0
     ip%HeIIfrac = 0.0
     ip%gammaHeI = 0.0
     ip%gammaHeII = 0.0
#endif 
  endif
  
end subroutine set_gammas

!> sets the global cooling rate for an ionization particle
!================================================================
subroutine set_cooling_func(ip,k)
use physical_constants_mod, only: HI_th_erg, HeI_th_erg, HeII_th_erg
use CenAtomicRates, only: Haiman_Bremss_cool
use CenAtomicRates, only: Haiman_Comp_Heol
use atomic_rates_mod, only: atomic_rates_type

  type(ionpart_type), intent(inout) :: ip   !< ionization particle  
  type(atomic_rates_type), intent(in) :: k  !< rates

  ! CIC  = collisional ionization cooling
  ! CEC  = collisional excitation cooling
  ! RCCb = recombination cooling case B
  ! PH   = photo heating

  ip%CIC =  ( k%HIcic   * ip%nHI  ) * ip%ne 
  ip%CEC =  ( k%HIcec   * ip%nHI  ) * ip%ne
  ip%RCCb = ( k%HIIrccB * ip%nHII ) * ip%ne

  ip%PH = ip%pdepr * ip%HIfrac * (ip%penrg - HI_th_erg)

#ifdef incHe
  ip%CIC  = ip%CIC + &
          ( k%HeIcic  * ip%nHeI  + &
            k%HeIIcic * ip%nHeII ) * ip%ne
  ip%CEC  = ip%CEC + &
          ( k%HeIcec  * ip%nHeII * ip%ne + & ! yes there is an extra ne
            k%HeIIcec * ip%nHeII ) * ip%ne
  ip%RCCb = ip%RCCb + &
          ( k%HeIIrccB  * ip%nHeII  + &
            k%HeIIIrccB * ip%nHeIII ) * ip%ne 

  ip%PH = ip%PH + ip%pdepr * ip%HeIfrac  * (ip%penrg - HeI_th_erg ) 
  ip%PH = ip%PH + ip%pdepr * ip%HeIIfrac * (ip%penrg - HeII_th_erg) 
#endif

  ip%PH = ip%PH / ip%cm3
     
  ip%BREM = Haiman_Bremss_cool(ip%T,ip%nHII,ip%nHeII,ip%nHeIII,ip%ne)
  ip%COMP = Haiman_Comp_Heol(ip%T,ip%Tcmb,ip%ne)

end subroutine set_cooling_func

!> sets the time rate of change of the electron number density
!================================================================
subroutine set_dnedt(ip,k)
use atomic_rates_mod, only: atomic_rates_type

  type(ionpart_type), intent(inout) :: ip   !< ionization particle
  type(atomic_rates_type), intent(in) :: k  !< rates

  real :: RRHII

#ifdef incHe
  real :: RRHeII, RRHeIII
#endif

  ip%GGHI   = k%HIci   * ip%ne + ip%gammaHI 
  ip%RRHIIa = k%HIIrcA * ip%ne
  ip%RRHIIb = k%HIIrcB * ip%ne

#ifdef increc
  RRHII   = ip%RRHIIa
#else
  RRHII   = ip%RRHIIb
#endif

  ip%dnedt = ip%GGHI * ip%nHI - RRHII * ip%nHII




#ifdef incHe

  ip%GGHeI  = ip%gammaHeI  + k%HeIci * ip%ne
  ip%GGHeII = ip%gammaHeII + k%HeIci * ip%ne
  ip%RRHeIIa  = k%HeIIrcA  * ip%ne
  ip%RRHeIIIa = k%HeIIIrcA * ip%ne
  ip%RRHeIIb  = k%HeIIrcB  * ip%ne
  ip%RRHeIIIb = k%HeIIIrcB * ip%ne

#ifdef increc
  RRHeII  = ip%RRHeIIa
  RRHeIII = ip%RRHeIIIa
#else
  RRHeII  = ip%RRHeIIb
  RRHeIII = ip%RRHeIIIb
#endif

  ip%dnedt =       ip%dnedt                          &
               +   ip%GGHeI    * ip%nHeI             &
               +   ip%GGHeII   * ip%nHeII            &
               -   RRHeII      * ip%nHeII            &
               -   RRHeIII     * ip%nHeIII
#endif

end subroutine set_dnedt

end module ionpar_mod
