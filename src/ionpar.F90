!> \file ionpar.F90

!> \brief the module that handles ionization particles
!<

module ionpar_mod
use atomic_rates_mod, only: atomic_rates_type
use particle_system_mod, only: particle_type
use global_mod, only: global_variables_type
use raylist_mod, only: raylist_type
use myf90_mod
use physical_constants_mod
use b2cd_mod, only: b2cdfac
use cen_atomic_rates_mod, only: Osterbrok_HI_photo_cs
use cen_atomic_rates_mod, only: Osterbrok_HeI_photo_cs
use cen_atomic_rates_mod, only: Osterbrok_HeII_photo_cs
use cen_atomic_rates_mod, only: Haiman_Bremss_cool
use cen_atomic_rates_mod, only: Haiman_Comp_Heol
implicit none

!---------------------------
!> ionization particle type. 
type ionpart_type

   integer(i8b) :: id      !< particle id
   real(r8b) :: pos(3)       !< x,y,z coordinates
   real(r8b) :: vel(3)       !< x,y,z velocities
   real(r8b) :: hsml         !< smoothing length
   real(r8b) :: rho          !< density = mass * NsphNnb / hsml^3 
   real(r8b) :: mass         !< particle mass
   real(r8b) :: T            !< temperature in K       
   real(r8b) :: xHI     !< HII ionization fraction
   real(r8b) :: xHII    !< nHII/nH
   real(r8b) :: xHeI    !< HeII ionization fraction
   real(r8b) :: xHeII   !< HeIII ionization fraction
   real(r8b) :: xHeIII  !< nHeIII/nHe
   integer(i8b) :: lasthit !< last ray to cross this particle

   real(r8b) :: xHI_in      !< initial xHI
   real(r8b) :: xHII_in     !< initial xHII
   real(r8b) :: xHeI_in     !< initial xHeI
   real(r8b) :: xHeII_in    !< initial xHeII
   real(r8b) :: xHeIII_in   !< initial xHeIII
   real(r8b) :: T_in        !< initial T

   real(r8b) :: DH(2,2)   !< dxH/dt array
   real(r8b) :: DHe(3,3)  !< dxHe/dt array

   integer(i8b) :: rayn   !< ray number
   integer(i8b) :: iter   !< number of the iteration in the solver
   integer(i8b) :: impact !< index of the particle in the raylist
   integer(i8b) :: indx   !< index of the particle in the psys
   real(r8b) :: d         !< distance along ray (code)
   real(r8b) :: b         !< impact parameter 

   real(r8b) :: bnorm     !< impact parameter normalized to smoothing length
   real(r8b) :: cdfac     !< column depth conversion factor
   real(r8b) :: gpercm3   !< density in cgs
   real(r8b) :: cm3       !< volume in cgs

   real(r8b) :: dt_code  !< time step for this par (code)
   real(r8b) :: dt_s     !< time step for this par (s)

   real(r8b) :: pflux       !< photons per second arriving at the particle
   real(r8b) :: fracabsorb  !< fraction of all photons absorbed
   real(r8b) :: pdepr       !< photons per second deposited in particle
   real(r8b) :: penrg       !< energy of one photon in the ray (ergs)
   real(r8b) :: pdeps       !< total photons deposited into particle

!   real(r8b) :: HIcolions   !< total HI collisional ionizations in the particle
!   real(r8b) :: HeIcolions  !< total HI collisional ionizations in the particle
!   real(r8b) :: HeIIcolions !< total HI collisional ionizations in the particle

!   real(r8b) :: HIIrcmbsB   !< total HII recombinations excluding to n=1 lvl
!   real(r8b) :: HIIrcmbsA   !< total HII recombinations to all levels
!   real(r8b) :: HeIIrcmbsB  !< total HeII recombinations excluding to n=1 lvl
!   real(r8b) :: HeIIrcmbsA  !< total HeII recombinations to all levels
!   real(r8b) :: HeIIIrcmbsB !< total HeIII recombinations excluding to n=1 lvl
!   real(r8b) :: HeIIIrcmbsA !< total HeIII recombinations to all levels


   real(r8b) :: Tcmb     !< background radiation field temperature
   real(r8b) :: H_mf     !< Hydrogen mass fraction
   real(r8b) :: He_mf    !< Helium mass fraction

   real(r8b) :: nH       !< number density of H
   real(r8b) :: nHe      !< number density of He
   real(r8b) :: Hcnt     !< number of H nuclei
   real(r8b) :: Hecnt    !< number of He nuclei

   real(r8b) :: NeBckgnd !< number density of metallic electrons
   real(r8b) :: ne       !< number density of electrons
   real(r8b) :: dnedt    !< time rate of change of ne

   real(r8b) :: nHI      !< number density of HI
   real(r8b) :: nHII     !< number density of HII
   real(r8b) :: nHeI     !< number density of HeI
   real(r8b) :: nHeII    !< number density of HeII
   real(r8b) :: nHeIII   !< number density of HeIII

   real(r8b) :: HIcnt    !< number of HI atoms
   real(r8b) :: HIIcnt   !< number of HII atoms
   real(r8b) :: HeIcnt   !< number of HeI atoms
   real(r8b) :: HeIIcnt  !< number of HeII atoms
   real(r8b) :: HeIIIcnt !< number of HeIII atoms
   real(r8b) :: Allcnt   !< number of all photo absorbing species (HI,HeI,HeII)

   real(r8b) :: tauHI    !< HI optical depth
   real(r8b) :: tauHeI   !< HeI optical depth
   real(r8b) :: tauHeII  !< HeII optical depth
   real(r8b) :: tausum   !< sum of all optical depths

   real(r8b) :: HItaufac   !< 1-exp(-tauHI)
   real(r8b) :: HeItaufac  !< 1-exp(-tauHeI)
   real(r8b) :: HeIItaufac !< 1-exp(-tauHeII)
   real(r8b) :: taufacsum  !< sum of all tau factors

   real(r8b) :: HIfrac   !< fraction of photons absorbed by HI
   real(r8b) :: HeIfrac  !< fraction of photons absorbed by HeI
   real(r8b) :: HeIIfrac !< fraction of photons absorbed by HeII

   real(r8b) :: sigmaHI   !< HI photo cross-section
   real(r8b) :: sigmaHeI  !< HeI photo cross-section
   real(r8b) :: sigmaHeII !< HeII photo cross-section

   real(r8b) :: gammaHI   !< HI photoionization rate
   real(r8b) :: gammaHeI  !< HeI photoionization rate
   real(r8b) :: gammaHeII !< HeII photoionization rate
   real(r8b) :: gammasum  !< sum of all photoionization rates

   real(r8b) :: GG     !< HIci * ne
   real(r8b) :: GGp    !< HI photoionization rate + HIci * ne
   real(r8b) :: RRa    !< HII recombination rate case A * ne
   real(r8b) :: RRb    !< HII recombination rate case B * ne

   real(r8b) :: GGI    !< HeIci * ne
   real(r8b) :: GGII   !< HeIIci * ne
   real(r8b) :: GGIp   !< HeI photoionization rate + HeIci * ne
   real(r8b) :: GGIIp  !< HeII photoionization rate + HeIIci * ne
   real(r8b) :: RRIIa  !< HeII recombination rate case A * ne
   real(r8b) :: RRIIIa !< HeIII recombination rate case A * ne
   real(r8b) :: RRIIb  !< HeII recombination rate case B * ne
   real(r8b) :: RRIIIb !< HeIII recombination rate case B * ne

   real(r8b) :: CIC   !< collisional ionization cooling rate
   real(r8b) :: CEC   !< collisional excitation cooling rate
   real(r8b) :: RCC   !< recombination cooling 
   real(r8b) :: PH    !< photo heating rate
   real(r8b) :: BREM  !< bremsstrahlung cooling rate
   real(r8b) :: COMP  !< compton heating/cooling rate
   real(r8b) :: COOL  !< total cooling function w/o photo heating
   real(r8b) :: COOLp !< total cooling function w photo heating

   real(r8b) :: tion     !< ionization time (s)
   real(r8b) :: tcool    !< cooling time (s)
   real(r8b) :: u        !< energy of particle (ergs)
   real(r8b) :: dudt     !< time rate of change of energy
   real(r8b) :: dTdt     !< time rate of change of temperature
   
   real(r8b) :: pdeps_eq    !< photons deposited in equilibrium conditions

   character(200) :: strtag !< labels code landmarks for debugging 

end type ionpart_type


contains

!> initializes the ionization particle values
subroutine initialize_ionpar(ipar,par,GV,srcray,He,raylist,impact)

  type(ionpart_type), intent(inout) :: ipar           !< ionization particle
  type(particle_type), intent(in) :: par              !< standard particle
  type(global_variables_type), intent(in) :: GV       !< global variables
  logical, intent(in) :: srcray                       !< source ray update ?
  logical, intent(in) :: He                           !< update Helium ?
  type(raylist_type), intent(in), optional :: raylist !< optional raylist
  integer(i8b), intent(in), optional :: impact        !< optional impact number
  
  real(r8b) :: mass_cgs
  real(r8b) :: rho_cgs

  call par2ionpar(par,ipar)
  ipar%iter = 0

  ipar%xHI_in = ipar%xHI
  ipar%xHII_in = ipar%xHII
  
  if (He) then
     ipar%xHeI_in = ipar%xHeI
     ipar%xHeII_in = ipar%xHeII
     ipar%xHeIII_in = ipar%xHeIII
  else
     ipar%xHeI_in = 0.0d0
     ipar%xHeII_in = 0.0d0
     ipar%xHeIII_in = 0.0d0
  end if

  ipar%T_in = ipar%T

  ipar%rayn = GV%rayn
  if (present(impact)) ipar%impact = impact
  ipar%NeBckgnd = GV%NeBackground
  ipar%H_mf = GV%H_mf
  ipar%He_mf = GV%He_mf
  ipar%Tcmb = GV%Tcmb_cur

  mass_cgs = ipar%mass * GV%cgs_mass / GV%LittleH
  rho_cgs  = ipar%rho * (GV%cgs_mass / (GV%cgs_len*GV%cgs_len*GV%cgs_len)) * GV%LittleH * GV%LittleH
 
  ipar%gpercm3 = rho_cgs
  ipar%cm3 = mass_cgs / rho_cgs

  ipar%nH = rho_cgs * ipar%H_mf  / M_H
  ipar%Hcnt = mass_cgs * ipar%H_mf / M_H  

  if (He) then
     ipar%nHe = rho_cgs * ipar%He_mf / M_He
     ipar%Hecnt = mass_cgs * ipar%He_mf / M_He
  else
     ipar%nHe = 0.0d0
     ipar%Hecnt = 0.0d0
  end if


  if (present(raylist)) then

     ipar%d = raylist%intersection(impact)%d
     ipar%b = raylist%intersection(impact)%b
     ipar%bnorm = ipar%b/ipar%hsml
     ipar%cdfac = b2cdfac(ipar%bnorm,ipar%hsml,GV%cgs_len/GV%LittleH)

     ipar%dt_code = (GV%rayn - ipar%lasthit) * GV%dtray_code
     ipar%dt_s    = (GV%rayn - ipar%lasthit) * GV%dtray_s
     
     if (srcray) then        
        ipar%pflux = raylist%ray%pcnt / ipar%dt_s 
     else
        ipar%pflux = raylist%ray%pcnt 
     end if

     ipar%penrg = raylist%ray%enrg
     ipar%pdeps = 0.0d0
     ipar%pdeps_eq = 0.0d0

     ipar%sigmaHI = Osterbrok_HI_photo_cs(raylist%ray%freq * HI_th_Hz)    
     if (He) then
        ipar%sigmaHeI = Osterbrok_HeI_photo_cs(raylist%ray%freq * HI_th_Hz)    
        ipar%sigmaHeII = Osterbrok_HeII_photo_cs(raylist%ray%freq * HI_th_Hz)
     else
        ipar%sigmaHeI = 0.0d0
        ipar%sigmaHeII = 0.0d0
     end if
     
  end if


  ipar%strtag = "in initialize_ionpar"
  call check_x(ipar)
     
end subroutine initialize_ionpar

!-----------------------------------------------------------------------
!> copies the basic particle data into an ionization particle
subroutine par2ionpar(par,ipar)
 type(particle_type), intent(in) :: par     !< input particle
 type(ionpart_type), intent(inout) :: ipar  !< output ionization particle


 ipar%id = par%id
 ipar%pos = par%pos

#ifdef incVel
 ipar%vel = par%vel
#else
 ipar%vel = 0.0
#endif

 ipar%hsml = par%hsml
 ipar%rho = par%rho
 ipar%mass = par%mass
 ipar%T = par%T

 ipar%xHI = par%xHI
 ipar%xHII = par%xHII

#ifdef incHe
 ipar%xHeI = par%xHeI
 ipar%xHeII = par%xHeII
 ipar%xHeIII = par%xHeIII
#else
 ipar%xHeI  = 0.0d0
 ipar%xHeII = 0.0d0
 ipar%xHeIII = 0.0d0
#endif

 ipar%lasthit = par%lasthit

end subroutine par2ionpar

!-----------------------------------------------------------------------
!> copies the ionization particle data into a basic particle 
subroutine ionpar2par(ipar,par)
 type(ionpart_type), intent(in) :: ipar    !< input ionization particle
 type(particle_type), intent(inout) :: par !< output particle

 par%T = ipar%T

 par%xHI = ipar%xHI
 par%xHII = ipar%xHII

#ifdef incHe
 par%xHeI = ipar%xHeI
 par%xHeII = ipar%xHeII
 par%xHeIII = ipar%xHeIII
#endif

 par%lasthit = ipar%lasthit

end subroutine ionpar2par



!-----------------------------------------------------
!> prints ionization particle information to screen
subroutine ionpar2screen(ipar)

  type(ionpart_type), intent(inout) :: ipar  !< ionization particle
  
   95 format(A,T25,I15)
   100 format(A,T25,3F12.4)
   102 format(A,T25,2F12.4)
   105 format(A,T25,1F12.4)
   110 format(A,T25,3ES18.9)
   115 format(A,T25,2ES18.9)
   write(*,*) 
   write(*,*) "ID", ipar%id
   write(*,*) "pos", ipar%pos
   write(*,*) "vel", ipar%vel
   write(*,*) "hsml", ipar%hsml
   write(*,*) "rho", ipar%rho
   write(*,*) "mass", ipar%mass
   write(*,*) "T", ipar%T
   write(*,*) "lasthit", ipar%lasthit
   write(*,*) 
   write(*,*) "ray num", ipar%rayn
   write(*,*) "psys indx", ipar%indx
   write(*,*) "impact", ipar%impact
   write(*,*) "iteration", ipar%iter
   write(*,*) "raydist", ipar%d
   write(*,*) 
   write(*,*) "string tag = ", trim(ipar%strtag)
   write(*,*) "xHI+xHII", ipar%xHI + ipar%xHII
   write(*,*) "xHeI+xHeII+xHeIII", ipar%xHeI + ipar%xHeII + ipar%xHeIII
   write(*,'(A,T25,ES18.9)') "T_in", ipar%T_in
   write(*,110) "b, hsml, bnorm",ipar%b, ipar%hsml, ipar%bnorm
   write(*,*) 
   write(*,*) "penrg/(HIth,HeIth,HeIIth)"
   write(*,*) ipar%penrg / HI_th_erg, ipar%penrg / HeI_th_erg, &
              ipar%penrg / HeII_th_erg
   write(*,*) 
   write(*,*) "input(xHI,xHII / xHeI,xHeII,xHeIII)"
   write(*,*) ipar%xHI_in, ipar%xHII_in
   write(*,*) ipar%xHeI_in, ipar%xHeII_in, ipar%xHeIII_in
   write(*,*) 
   write(*,*) "current(xHI,xHII / xHeI,xHeII,xHeIII)"
   write(*,*) ipar%xHI, ipar%xHII
   write(*,*) ipar%xHeI, ipar%xHeII, ipar%xHeIII
   write(*,*) 
   write(*,*) "current(nHI,nHII / nHeI,nHeII,nHeIII)"
   write(*,*) ipar%nHI, ipar%nHII
   write(*,*) ipar%nHeI, ipar%nHeII, ipar%nHeIII
   write(*,*) 
   write(*,*) "optical depths / tau facs (HI,HeI,HeII)"
   write(*,*) ipar%tauHI, ipar%tauHeI, ipar%tauHeII
   write(*,*) ipar%HItaufac, ipar%HeItaufac, ipar%HeIItaufac
   write(*,*) 
   write(*,*) "nH,nHe,ne / Hcnt,Hecnt,dnedt"
   write(*,*)  ipar%nH, ipar%nHe, ipar%ne
   write(*,*)  ipar%Hcnt, ipar%Hecnt, ipar%dnedt
   write(*,*) 
   write(*,*) "HIcnt,HeIcnt,HeIIcnt / HIfrac,HeIfrac,HeIIfrac"
   write(*,*) ipar%HIcnt, ipar%HeIcnt, ipar%HeIIcnt 
   write(*,*) ipar%HIfrac,ipar%HeIfrac,ipar%HeIIfrac
   write(*,*)
   write(*,*) "sigmas (HI,HeI,HeII) / Gammas (HI,HeI,HeII)"
   write(*,*) ipar%sigmaHI, ipar%sigmaHeI, ipar%sigmaHeII
   write(*,*) ipar%gammaHI, ipar%gammaHeI, ipar%gammaHeII
   write(*,*)
   write(*,*) "GGp (HI,HeI,HeII) / RRb (HII,HeII,HeIII)"
   write(*,*) ipar%GGp, ipar%GGIp, ipar%GGIIp
   write(*,*) ipar%RRb, ipar%RRIIb, ipar%RRIIIb
   write(*,*)
   write(*,*) "CIC,CEC,RCC / PH,BREM,COMP / COOLp, COOL"
   write(*,*) ipar%CIC, ipar%CEC, ipar%RCC
   write(*,*) ipar%PH, ipar%BREM, ipar%COMP
   write(*,*) ipar%COOLp, ipar%COOL
   write(*,*) 
   write(*,*) "pflux, fracabsorb, pdepr"
   write(*,*) ipar%pflux, ipar%fracabsorb, ipar%pdepr
   write(*,*) 
   write(*,110) "pdepr,penrg,pdeps", ipar%pdepr, ipar%penrg, ipar%pdeps

   write(*,*)
   write(*,115) "cdfac,d", ipar%cdfac,  ipar%d
   write(*,*) 
   write(*,115) "gpercm3,cm3", ipar%gpercm3, ipar%cm3

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

   ipar%strtag = ""

 end subroutine ionpar2screen


!================================================================
!> dummy checks the bounds on HII, HeII, HeIII, and T
subroutine check_x(ip) 

  real(r8b), parameter :: zero = 0.0d0
  real(r8b), parameter :: one = 1.0d0

  real(r8b), parameter :: TOL = 0.0d+0
  type(ionpart_type), intent(inout) :: ip  !< ionization particle
  logical :: bad

  100 format(A,I2,A)
  101 format(A,I2,A,F15.8)
  102 format(A,I15)
  103 format(A,2F15.8)

  bad = .false.
 
  if ( ip%xHI .LT. zero - TOL ) then
     bad = .true.
     write(*,100) "xHI < zero in check_x"
  end if
     
  if ( ip%xHI > one + TOL ) then
     bad = .true.
     write(*,*) "xHI = ", ip%xHI
     write(*,100) "xHI > one in check_x"
  end if

  if ( ip%xHII .LT. zero - TOL ) then
     bad = .true.
     write(*,100) "xHII < zero in check_x"
  end if
     
  if ( ip%xHII .GT. one + TOL ) then
     bad = .true.
     write(*,100) "xHII > one in check_x"
  end if
     

  if ( ip%xHeI .LT. zero - TOL ) then
     bad = .true.
     write(*,100) "xHeI < zero in check_x"
  end if
     
  if ( ip%xHeI .GT. one + TOL ) then
     bad = .true.
     write(*,100) "xHeI > one in check_x"
  end if

  if ( ip%xHeII .LT. zero - TOL ) then
     bad = .true.
     write(*,100) "xHeII < zero in check_x"
  end if
     
  if ( ip%xHeII .GT. one + TOL ) then
     bad = .true.
     write(*,100) "xHeII > one in check_x"
  end if

  if ( ip%xHeIII .LT. zero - TOL ) then
     bad = .true.
     write(*,100) "xHeIII < zero in check_x"
  end if
     
  if ( ip%xHeIII .GT. one + TOL ) then
     bad = .true.
     write(*,100) "xHeIII > one in check_x"
  end if

  if ( ip%T .LE. zero ) then
     bad = .true.
     write(*,*) "T <= zero in check_x"
  end if


  if (bad) then
     call ionpar2screen(ip)
     stop
  end if

end subroutine check_x


!> sets the photoionization rate for an ionization particle
!================================================================
subroutine set_taus(ip,He)

  real(r8b), parameter :: zero = 0.0d0
  real(r8b), parameter :: one = 1.0d0
  real(r8b), parameter :: TAU_LOW = 1.0d-4
  real(r8b), parameter :: TAU_HIGH = 3.0d1
  type(ionpart_type), intent(inout) :: ip !< ionization particle  
  logical, intent(in) :: He

  logical :: HI,HeI,HeII


  ! check which species will be absorbing
  !---------------------------------------
  HI   = .false.
  HeI  = .false.
  HeII = .false.

  if (ip%xHI > zero) HI = .true.

  if (He) then
     if (ip%xHeI  > zero .and. ip%penrg > HeI_th_erg)  HeI = .true.
     if (ip%xHeII > zero .and. ip%penrg > HeII_th_erg) HeII = .true.
  end if


  ! calculate atom counts
  !---------------------------------------
  ip%HIcnt = ip%Hcnt * ip%xHI       
  if (He) then
     ip%HeIcnt  = ip%Hecnt * ip%xHeI        
     ip%HeIIcnt = ip%Hecnt * ip%xHeII          
  else
     ip%HeIcnt = zero
     ip%HeIIcnt = zero
  end if

  ip%Allcnt = ip%HIcnt + ip%HeIcnt + ip%HeIIcnt


  ! calculate taus
  !---------------------------------------
  if (HI) then
     ip%tauHI = ip%cdfac * ip%HIcnt * ip%sigmaHI
  else
     ip%tauHI = zero
  end if

  !---------------------------------------
  if (HeI) then
     ip%tauHeI = ip%cdfac * ip%HeIcnt * ip%sigmaHeI
  else
     ip%tauHeI = zero
  end if
     
  !---------------------------------------
  if (HeII) then
     ip%tauHeII = ip%cdfac * ip%HeIIcnt * ip%sigmaHeII
  else 
     ip%tauHeII = zero
  end if
  
  ip%tausum = ip%tauHI + ip%tauHeI + ip%tauHeII


  ! calculate absorption ratios
  !--------------------------------
  if (HI) then

     if (ip%tauHI < TAU_LOW) then
        ip%HItaufac = ip%tauHI
     else if (ip%tauHI > TAU_HIGH) then
        ip%HItaufac = one
     else
        ip%HItaufac = one - exp(-ip%tauHI)
     end if

  else
     
     ip%HItaufac = zero

  end if

  !-------------------------------------
  if (HeI) then

     if (ip%tauHeI < TAU_LOW) then
        ip%HeItaufac = ip%tauHeI
     else if (ip%tauHeI > TAU_HIGH) then
        ip%HeItaufac = one
     else
        ip%HeItaufac = one - exp(-ip%tauHeI)
     end if

  else
     
     ip%HeItaufac = zero

  end if

  !-------------------------------------
  if (HeII) then

     if (ip%tauHeII < TAU_LOW) then
        ip%HeIItaufac = ip%tauHeII
     else if (ip%tauHeII > TAU_HIGH) then
        ip%HeIItaufac = one
     else
        ip%HeIItaufac = one - exp(-ip%tauHeII)
     end if

  else
     
     ip%HeIItaufac = zero

  end if

  ip%taufacsum = ip%HItaufac + ip%HeItaufac + ip%HeIItaufac

  if (ip%taufacsum > zero) then
  
     if (.not. He) then
        ip%HIfrac   = one
        ip%HeIfrac  = zero
        ip%HeIIfrac = zero
     else
        ip%HIfrac   = ip%HItaufac   / ip%taufacsum
        ip%HeIfrac  = ip%HeItaufac  / ip%taufacsum
        ip%HeIIfrac = ip%HeIItaufac / ip%taufacsum
     end if

  else

     ip%HIfrac   = zero
     ip%HeIfrac  = zero
     ip%HeIIfrac = zero
     
  end if

end subroutine set_taus


!> sets the photoionization rate for an ionization particle
!================================================================
subroutine set_gammas(ip,He)

  real(r8b), parameter :: zero = 0.0d0
  real(r8b), parameter :: one = 1.0d0
  real(r8b), parameter :: TAU_LOW = 1.0d-4
  real(r8b), parameter :: TAU_HIGH = 3.0d+1
  type(ionpart_type), intent(inout) :: ip !< ionization particle  
  logical, intent(in) :: He

  
  call set_taus(ip,He)

  ! photoionization rates
  !------------------------
  if (ip%tausum > zero) then

     if (ip%tausum < TAU_LOW) then
        ip%fracabsorb = ip%tausum
     else if (ip%tausum > TAU_HIGH) then
        ip%fracabsorb = one
     else
        ip%fracabsorb = one - exp(-ip%tausum)
     end if

     ip%pdepr    = ip%pflux * ip%fracabsorb
     ip%gammasum = ip%pdepr / ip%Allcnt

     ip%gammaHI = ip%gammasum * ip%HIfrac 
     if (He) then
        ip%gammaHeI  = ip%gammasum * ip%HeIfrac 
        ip%gammaHeII = ip%gammasum * ip%HeIIfrac 
     end if

  else

     ip%gammasum = zero
     ip%gammaHI = zero
     if (He) then
        ip%gammaHeI = zero
        ip%gammaHeII = zero
     end if
     return

  end if
  
end subroutine set_gammas




!> sets the time rate of change of the electron number density
!================================================================
subroutine set_dnedt(ip,photo,caseA,He)

  type(ionpart_type), intent(inout) :: ip   !< ionization particle
  logical, intent(in) :: photo
  logical, intent(in) :: caseA(2)
  logical, intent(in) :: He

  real(r8b) :: GG,RR
  real(r8b) :: GGI,GGII,RRII,RRIII
  
  if (photo) then
     GG = ip%GGp
  else
     GG = ip%GG
  end if

  if (caseA(1)) then
     RR = ip%RRa
  else
     RR = ip%RRb
  end if

  ip%dnedt = GG * ip%nHI - RR * ip%nHII

  if (He) then

     if (photo) then
        GGI  = ip%GGIp
        GGII = ip%GGIIp
     else
        GGI  = ip%GGI
        GGII = ip%GGII
     end if

     if (caseA(2)) then
        RRII  = ip%RRIIa
        RRIII = ip%RRIIIa
     else
        RRII  = ip%RRIIb
        RRIII = ip%RRIIIb
     end if

     ip%dnedt = ip%dnedt +  (GGI   * ip%nHeI ) + (GGII  * ip%nHeII )   &
                         -  (RRII  * ip%nHeII) - (RRIII * ip%nHeIII)

  end if

end subroutine set_dnedt
 
!> sets the GG's and RR's from the atomic rates, ne's and gammas
!=================================================================
subroutine set_ionization_func(ip,k,photo,caseA,He)

  type(ionpart_type), intent(inout) :: ip   !< ionization particle  
  type(atomic_rates_type), intent(in) :: k  !< rates
  logical, intent(in) :: photo
  logical, intent(in) :: caseA(2)
  logical, intent(in) :: He
 
  ip%GG  = k%HIci * ip%ne 
  if (photo) ip%GGp = ip%GG + ip%gammaHI
  
  ip%RRb = k%HIIrcB * ip%ne
  if (caseA(1)) ip%RRa = k%HIIrcA * ip%ne

  if (He) then
     ip%GGI    = k%HeIci  * ip%ne 
     ip%GGII   = k%HeIIci * ip%ne 
     if (photo) then
        ip%GGIp  = ip%GGI  + ip%gammaHeI  
        ip%GGIIp = ip%GGII + ip%gammaHeII  
     end if

     ip%RRIIb  = k%HeIIrcB  * ip%ne
     ip%RRIIIb = k%HeIIIrcB * ip%ne

     if (caseA(2)) then
        ip%RRIIa  = k%HeIIrcA  * ip%ne
        ip%RRIIIa = k%HeIIIrcA * ip%ne
     end if
  end if

end subroutine set_ionization_func

 
!> sets the global cooling rate for an ionization particle
!================================================================
subroutine set_cooling_func(ip,k,photo,caseA,He)

  type(ionpart_type), intent(inout) :: ip   !< ionization particle  
  type(atomic_rates_type), intent(in) :: k  !< rates
  logical, intent(in) :: photo
  logical, intent(in) :: caseA(2)
  logical, intent(in) :: He

  ! CIC  = collisional ionization cooling
  ! CEC  = collisional excitation cooling
  ! RCC  = recombination cooling 
  ! BREM = free free cooling
  ! COMP = compton cooling
  ! PH   = photo heating  

  ! cooling from Hydrogen
  !-----------------------
  ip%PH  = 0.0d0
  ip%CIC = ( k%HIcic * ip%nHI  ) * ip%ne 
  ip%CEC = ( k%HIcec * ip%nHI  ) * ip%ne

  if (caseA(1)) then
     ip%RCC = ( k%HIIrccA * ip%nHII ) * ip%ne
  else
     ip%RCC = ( k%HIIrccB * ip%nHII ) * ip%ne
  end if

  if (photo) then
     ip%PH = ip%pdepr * ip%HIfrac * (ip%penrg - HI_th_erg)
  end if

  ! cooling from Helium (yes there is an extra ne in He CEC)
  !----------------------------------------------------------
  if (He) then
     ip%CIC = ip%CIC + ( k%HeIcic  * ip%nHeI  + k%HeIIcic * ip%nHeII ) * ip%ne
     ip%CEC = ip%CEC + ( k%HeIcec  * ip%nHeII * ip%ne + k%HeIIcec * ip%nHeII ) * ip%ne

     if (caseA(2)) then
        ip%RCC = ip%RCC + ( k%HeIIrccA * ip%nHeII + k%HeIIIrccA * ip%nHeIII ) * ip%ne 
     else
        ip%RCC = ip%RCC + ( k%HeIIrccB * ip%nHeII + k%HeIIIrccB * ip%nHeIII ) * ip%ne 
     end if

     if (photo) then
        ip%PH = ip%PH + ip%pdepr * ip%HeIfrac  * (ip%penrg - HeI_th_erg ) 
        ip%PH = ip%PH + ip%pdepr * ip%HeIIfrac * (ip%penrg - HeII_th_erg) 
     end if
  end if

  ip%PH = ip%PH / ip%cm3
  ip%BREM = Haiman_Bremss_cool(ip%T,ip%nHII,ip%nHeII,ip%nHeIII,ip%ne)
  ip%COMP = Haiman_Comp_Heol(ip%T,ip%Tcmb,ip%ne)

  ip%COOL  = -( ip%CIC + ip%CEC + ip%RCC + ip%BREM + ip%COMP )
  ip%COOLp = ip%PH + ip%COOL 
 
  
end subroutine set_cooling_func 

!> sets the Hydrogen derivative matrix from the values in GGRR
!======================================================================
subroutine setDH(ip,photo,caseA)

  type(ionpart_type), intent(inout) :: ip   !< ionization particle  
  logical, intent(in) :: photo
  logical, intent(in) :: caseA
  real(r8b) :: GG,RR
 
  if (photo) then
     GG = ip%GGp
  else
     GG = ip%GG
  end if

  if (caseA) then
     RR = ip%RRa
  else
     RR = ip%RRb
  end if

  ip%DH(1,1) = -GG
  ip%DH(2,1) = GG

  ip%DH(1,2) = RR
  ip%DH(2,2) = -RR

end subroutine setDH


!> sets the Helium derivative matrix from the values in GGRR
!======================================================================
subroutine setDHe(ip,photo,caseA)

  real(r8b), parameter :: zero = 0.0d0
  type(ionpart_type), intent(inout) :: ip   !< ionization particle  
  logical, intent(in) :: photo
  logical, intent(in) :: caseA

  real(r8b) :: GGI,GGII,RRII,RRIII

  if (photo) then
     GGI  = ip%GGIp
     GGII = ip%GGIIp
  else
     GGI  = ip%GGI
     GGII = ip%GGII     
  end if

  if (caseA) then
     RRII  = ip%RRIIa
     RRIII = ip%RRIIIa
  else
     RRII  = ip%RRIIb
     RRIII = ip%RRIIIb
  end if

  ip%DHe(1,1) = -GGI
  ip%DHe(2,1) = GGI
  ip%DHe(3,1) = zero

  ip%DHe(1,2) = RRII
  ip%DHe(2,2) = -(GGII + RRII)
  ip%DHe(3,2) = GGII

  ip%DHe(1,3) = zero
  ip%DHe(2,3) = RRIII
  ip%DHe(3,3) = -RRII

  
end subroutine setDHe



end module ionpar_mod
