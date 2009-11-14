!> \file iliev_comparison_project.f90

!> \brief Module that does special output for the 
!! Iliev Comparison Project (ICP) tests
!<
module iliev_comparison_project_mod
use myf90_mod
use raylist_mod, only: raylist_type
use ray_mod, only: ray_type
use particle_system_mod, only: particle_system_type
use oct_tree_mod, only: oct_tree_type
use global_mod, only: global_variables_type
use physical_constants_mod
implicit none

private

public :: initialize_iliev_tests
public :: iliev_test1_screen_out
public :: iliev_test2_screen_out
public :: iliev_test3_screen_out
public :: iliev_test4_screen_out

! some variables from astro-ph/0603199 the Cosmological RT Comparison 
! Static Density Field Tests
!---------------------------------------------------------------------

    real(r8b), parameter :: CmpPrjRecRate = 2.59d-13                       !< ICP RR rate (B) [cm^3/s]
    real(r8b), parameter :: CmpPrjColRate = 3.16d-16                       !< ICP CI rate [cm^3/s]
    real(r8b), parameter :: CmpPrjTest1xHII = 1.2d-3                       !< xHII(t=0) for T1
    real(r8b), parameter :: CmpPrjTest2xHII = 0.0                          !< xHII(t=0) for T2
    real(r8b), parameter :: CmpPrjTest1Dens = 1.0d-3                       !< nH for T1 [cm^-3]
    real(r8b), parameter :: CmpPrjTest2Dens = 1.0d-3                       !< nH for T2 [cm^-3]
    real(r8b), parameter :: CmpPrjTest3OutDens = 2.0d-4                    !< nH_out for T3 [cm^-3]
    real(r8b), parameter :: CmpPrjTest3InDens  = 200 * CmpPrjTest3OutDens  !< nH_in for T3 [cm^-3] 
    real(r8b), parameter :: CmpPrjTest1Flux = 5.0d48                       !< T1 photons / s
    real(r8b), parameter :: CmpPrjTest2Flux = 5.0d48                       !< T2 photons / s
    real(r8b), parameter :: CmpPrjTest3Flux = 1.0d6                        !< T3 photons / s / cm^2
    real(r8b), parameter :: CmpPrjTest1Temp = 1.0d4                        !< T(t=0) for T1[K]
    real(r8b), parameter :: CmpPrjTest2Temp = 1.0d2                        !< T(t=0) for test 2 [K]
    real(r8b), parameter :: CmpPrjTest3OutTemp = 8.0d3                     !< T_out(t=0) test 3 [K]
    real(r8b), parameter :: CmpPrjTest3InTemp = 4.0d1                      !< T_in(t=0) test 3 [K]

    real(r8b), parameter :: CmpPrjRecTime_s = 1.0 / (CmpPrjRecRate * CmpPrjTest1Dens)      !< T1 rec time [s]
    real(r8b), parameter :: CmpPrjColTime_s = 1.0 / (CmpPrjColRate * CmpPrjTest1Dens)      !< T1 col ion time [s]
    real(r8b), parameter :: CmpPrjRecTime_yrs = CmpPrjRecTime_s / year2sec                 !< T1 rec time [yrs]
    real(r8b), parameter :: CmpPrjColTime_yrs = CmpPrjColTime_s / year2sec                 !< T1 col ion time [yrs]
    real(r8b), parameter :: CmpPrjEqxHII = CmpPrjColRate / (CmpPrjRecRate + CmpPrjColRate) !< T1 col ion eq.

    
    real(r8b), save :: CmpPrjStromRad3      !< T1 Strom. Radius Cubed [cm^3]
    real(r8b), save :: CmpPrjStromRad_cm    !< T1 Strom. Radius [cm]
    real(r8b), save :: CmpPrjStromRad_kpc   !< T1 Strom. Radius [kpc]


contains

!------------------------------------------------
!> initialize iliev comparison project variables
  subroutine initialize_iliev_tests()
   
    CmpPrjStromRad3    = (3 * CmpPrjTest1Flux) / (4 * pi * CmpPrjRecRate * CmpPrjTest1Dens**2)  
    CmpPrjStromRad_cm  = exp( 1./3. * log(CmpPrjStromRad3) ) 
    CmpPrjStromRad_kpc = CmpPrjStromRad_cm / kpc2cm                

  end subroutine initialize_iliev_tests


!-------------------------------
!> iliev T1screen output
  subroutine iliev_test1_screen_out(psys,searchtree,GV)
  use ray_mod, only: make_probe_ray
  use raylist_mod, only: trace_ray
  use raylist_mod, only: kill_raylist, prepare_raysearch
  
    type(particle_system_type), intent(in) :: psys  !< input particle system  
    type(oct_tree_type), intent(in) :: searchtree !< oct-tree to search
    type(global_variables_type), intent(in) :: GV  !< global variables

    type(ray_type) :: ray
    type(raylist_type) :: raylist

    real(r8b) :: NWionfrac, time_ratio
    real(r8b) :: ionRa, ionRn, ionVn
    real(r8b) :: pos(3), dir(3)
    real(r8b) :: d,d1,d2,d3,d4
    integer :: i,pindx
    
!   compute analytic Stromgren radius
    time_ratio = GV%time_s / CmpPrjRecTime_s
    ionRa = CmpPrjStromRad_kpc * (1.0d0 - exp(-time_ratio) )**(1.0/3.0)
 
!   compute numerical Stromgren radius (Volume method)
    NWionfrac = number_weight_ionfrac(psys)
    if (psys%src(1)%EmisPrf == -2) then
       ionVn = ( (NWionfrac - CmpPrjTest1xHII) * 6.6**3 ) * 8  
    else
       ionVn = (NWionfrac - CmpPrjTest1xHII) * GV%BoxVol_kpc
    end if
    ionRn = (3.0 * ionVn / (4 * pi) )**(1.0/3.0)

!   compute numerical Stromgren radius (Ray probe method)
    call prepare_raysearch(psys,raylist)
    pos = (/0.0,0.0,0.0/)

    dir = (/1.0,2.0,1.0/)
    call make_probe_ray( pos,dir,ray )
    call trace_ray(ray,raylist,psys,searchtree) 
    do i = 1,raylist%nnb
       pindx = raylist%intersection(i)%pindx
       d = raylist%intersection(i)%d
       if (psys%par(pindx)%xHII < 0.5) exit
    end do
    d1 = raylist%intersection(i-1)%d + &
        (raylist%intersection(i)%d - raylist%intersection(i-1)%d)/2 
         
    dir = (/2.0,1.0,1.0/)
    call make_probe_ray( pos,dir,ray )
    call trace_ray(ray,raylist,psys,searchtree) 
    do i = 1,raylist%nnb
       pindx = raylist%intersection(i)%pindx
       d = raylist%intersection(i)%d
       if (psys%par(pindx)%xHII < 0.5) exit
    end do
    d2 = raylist%intersection(i-1)%d + &
        (raylist%intersection(i)%d - raylist%intersection(i-1)%d)/2 

    dir = (/1.0,1.0,2.0/)
    call make_probe_ray( pos,dir,ray )
    call trace_ray(ray,raylist,psys,searchtree) 
    do i = 1,raylist%nnb
       pindx = raylist%intersection(i)%pindx
       d = raylist%intersection(i)%d
       if (psys%par(pindx)%xHII < 0.5) exit
    end do
    d3 = raylist%intersection(i-1)%d + &
        (raylist%intersection(i)%d - raylist%intersection(i-1)%d)/2 

    dir = (/1.0,1.0,1.0/)
    call make_probe_ray( pos,dir,ray )
    call trace_ray(ray,raylist,psys,searchtree) 
    do i = 1,raylist%nnb
       pindx = raylist%intersection(i)%pindx
       d = raylist%intersection(i)%d
       if (psys%par(pindx)%xHII < 0.5) exit
    end do
    d4 = raylist%intersection(i-1)%d + &
        (raylist%intersection(i)%d - raylist%intersection(i-1)%d)/2 

    d = (d1+d2+d3+d4)/4

    call  kill_raylist(raylist)

    write(*,*) "=============================================================="
    write(*,*) " Iliev Test 1 Specific Monitoring "
    write(*,*) "=============================================================="
    write(*,*) "Stromgren Radii (kpc)"
    write(*,*) "Analytic                  = ", ionRa
    write(*,*) "Numerical (Volume Method) = ", ionRn
    write(*,*) "Numerical (Ray Method)    = ", d
    write(*,*) "numerical/analytic (Volume) = ", ionRn / ionRa
    write(*,*) "numerical/analytic (Ray)    = ", d / ionRa
    write(*,*) "=============================================================="
   
  end subroutine iliev_test1_screen_out


!-------------------------------
!> iliev test 2 screen output
  subroutine iliev_test2_screen_out(psys,searchtree,GV)
  use ray_mod, only: make_probe_ray
  use raylist_mod, only: trace_ray
  use raylist_mod, only: kill_raylist, prepare_raysearch
  
    type(particle_system_type), intent(in) :: psys !< input particle system  
    type(oct_tree_type), intent(in) :: searchtree !< oct-tree to search
    type(global_variables_type), intent(in) :: GV  !< global variables

    type(ray_type) :: ray
    type(raylist_type) :: raylist

    real(r8b) :: NWionfrac, time_ratio
    real(r8b) :: ionRa, ionRn, ionVn
    real(r8b) :: pos(3), dir(3)
    real(r8b) :: d,d1,d2,d3,d4
    integer :: i,pindx

!   compute analytic Stromgren radius
    time_ratio = GV%time_s / CmpPrjRecTime_s
    ionRa = CmpPrjStromRad_kpc * (1.0d0 - exp(-time_ratio) )**(1.0/3.0)
 
!   compute numerical Stromgren radius
    NWionfrac = number_weight_ionfrac(psys)
    ionVn = (NWionfrac - CmpPrjTest2xHII) * GV%BoxVol_kpc
    ionRn = (3.0 * ionVn / (4 * pi) )**(1.0/3.0)

!   compute numerical Stromgren radius (Ray probe method)
    call prepare_raysearch(psys,raylist)
    pos = (/0.0,0.0,0.0/)

    dir = (/1.0,0.0,0.0/)
    call make_probe_ray( pos,dir,ray )
    call trace_ray(ray,raylist,psys,searchtree) 
    do i = 1,raylist%nnb
       pindx = raylist%intersection(i)%pindx
       d = raylist%intersection(i)%d
       if (psys%par(pindx)%xHII < 0.5) exit
    end do
    d1 = raylist%intersection(i-1)%d + &
        (raylist%intersection(i)%d - raylist%intersection(i-1)%d)/2 
         
    dir = (/0.0,1.0,0.0/)
    call make_probe_ray( pos,dir,ray )
    call trace_ray(ray,raylist,psys,searchtree) 
    do i = 1,raylist%nnb
       pindx = raylist%intersection(i)%pindx
       d = raylist%intersection(i)%d
       if (psys%par(pindx)%xHII < 0.5) exit
    end do
    d2 = raylist%intersection(i-1)%d + &
        (raylist%intersection(i)%d - raylist%intersection(i-1)%d)/2 

    dir = (/0.0,0.0,1.0/)
    call make_probe_ray( pos,dir,ray )
    call trace_ray(ray,raylist,psys,searchtree) 
    do i = 1,raylist%nnb
       pindx = raylist%intersection(i)%pindx
       d = raylist%intersection(i)%d
       if (psys%par(pindx)%xHII < 0.5) exit
    end do
    d3 = raylist%intersection(i-1)%d + &
        (raylist%intersection(i)%d - raylist%intersection(i-1)%d)/2 

    dir = (/1.0,1.0,1.0/)
    call make_probe_ray( pos,dir,ray )
    call trace_ray(ray,raylist,psys,searchtree) 
    do i = 1,raylist%nnb
       pindx = raylist%intersection(i)%pindx
       d = raylist%intersection(i)%d
       if (psys%par(pindx)%xHII < 0.5) exit
    end do
    d4 = raylist%intersection(i-1)%d + &
        (raylist%intersection(i)%d - raylist%intersection(i-1)%d)/2 

    d = (d1+d2+d3+d4)/4

    call  kill_raylist(raylist)

    write(*,*) "=============================================================="
    write(*,*) " Iliev Test 2 Specific Monitoring "
    write(*,*) "=============================================================="
    write(*,*) "Stromgren Radii (kpc)"
    write(*,*) "Analytic                  = ", ionRa
    write(*,*) "Numerical (Volume Method) = ", ionRn
    write(*,*) "Numerical (Ray Method)    = ", d
    write(*,*) "numerical/analytic (Volume) = ", ionRn / ionRa
    write(*,*) "numerical/analytic (Ray)    = ", d / ionRa
    write(*,*) "=============================================================="


  end subroutine iliev_test2_screen_out


!-------------------------------
!> iliev test 3 screen output
  subroutine iliev_test3_screen_out(psys,GV)
 
    type(particle_system_type), intent(in) :: psys !< input particle system  
    type(global_variables_type), intent(in) :: GV  !< global variables

    real(r8b), parameter :: cpos(3) = (/5.0,3.3,3.3/)
    integer :: i, incount, outcount
    real(r8b) :: d, cionsum, oionsum, Tinsum, Toutsum

    ! real(r8b) :: pos(3)

    

    ! ionization fraction of the clump
    cionsum = 0.0
    Tinsum = 0.0
    incount = 0
    oionsum = 0.0
    Toutsum = 0.0
    outcount = 0
    do i = 1,size(psys%par)
       d = sqrt( sum( (cpos - psys%par(i)%pos)**2 ) )
       if (d < 0.8) then
          incount = incount + 1
          cionsum = cionsum + psys%par(i)%xHII
#ifdef incT
            Tinsum = Tinsum + psys%par(i)%T
#endif
       else 
          outcount = outcount + 1
          oionsum = oionsum + psys%par(i)%xHII
#ifdef incT
            Toutsum = Toutsum + psys%par(i)%T
#endif
       end if
    end do
    oionsum = oionsum / outcount 
    Toutsum = Toutsum / outcount
    cionsum = cionsum / incount 
    Tinsum = Tinsum / incount

    100 format(A,T40,2ES12.5)
    write(*,*) "=============================================================="
    write(*,*) GV%TestScenario, " Specific Monitoring "
    write(*,*) "=============================================================="
    write(*,100) "Ionization fraction in/outside clump: ", cionsum, oionsum
    write(*,100) "Temperature in/outisde clump:         ", Tinsum, Toutsum
    write(*,*) "=============================================================="

  end subroutine iliev_test3_screen_out


!-------------------------------
!> iliev test 4 screen output
  subroutine iliev_test4_screen_out(psys,GV)  
    type(particle_system_type), intent(in) :: psys !< input particle system  
    type(global_variables_type), intent(in) :: GV  !< global variables

    real(r8b) :: NWionfrac

    NWionfrac = number_weight_ionfrac(psys)

    write(*,*) "=============================================================="
    write(*,*) GV%TestScenario, " Specific Monitoring "
    write(*,*) "=============================================================="

  end subroutine iliev_test4_screen_out


!-----------------------------------------------------------------------
!> function to calculate the global number weighted ionization fraction
  function number_weight_ionfrac(psys) result(numionfrac)

     type(particle_system_type), intent(in) :: psys !< input particle system  
     real(r8b) :: numionfrac !< number weighted global ionization fraction
     integer :: i

       numionfrac = 0.0d0
       do i = 1,size(psys%par)
          numionfrac = numionfrac + psys%par(i)%xHII
       end do
       numionfrac = numionfrac / size(psys%par)

  end function number_weight_ionfrac



end module iliev_comparison_project_mod




