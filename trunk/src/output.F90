!> \file output.f90

!> \brief The module that handles output
!<

module output_mod
use myf90_mod
use gadget_input_mod, only: gadget_header_type, saved_gheads
use particle_system_mod, only: particle_system_type, particle_type
use oct_tree_mod, only: oct_tree_type
use physical_constants_mod
use global_mod, only: PLAN, GV
implicit none

contains
 
!=======================================================================
!=======================================================================
!
! OUTPUT ROUTINES
!
!=======================================================================
!=======================================================================




!> put the config choices in a sphray header for output.
! ==========================================================
 
  subroutine config_to_ghead(fnum, ghead)

    integer, intent(in) :: fnum                    !< which snapshot file
    type(gadget_header_type), intent(out) :: ghead !< gadget header

    ghead%npar_all = 0
    ghead%npar_all(1) = saved_gheads(GV%CurSnapNum,fnum)%npar_all(1)     

    ghead%npar_file = 0
    ghead%npar_file(1) = saved_gheads(GV%CurSnapNum,fnum)%npar_file(1)     

    ghead%mass = 0.0

    ghead%nfiles    = saved_gheads(GV%CurSnapNum,fnum)%nfiles

    if (GV%Comoving) then
       ghead%a = PLAN%snap(GV%CurSnapNum)%ScalefacAt
       ghead%z = 1.0d0 / ghead%a - 1.0d0
    else
       ghead%a = GV%time_code
       ghead%z = 1.0d0 / ghead%a - 1.0d0
    end if

    ghead%boxlen  = saved_gheads(GV%CurSnapNum,fnum)%boxlen
    ghead%OmegaM  = saved_gheads(GV%CurSnapNum,fnum)%OmegaM
    ghead%OmegaL  = saved_gheads(GV%CurSnapNum,fnum)%OmegaL
    ghead%h       = saved_gheads(GV%CurSnapNum,fnum)%h


    ghead%npar_hw = 0
    ghead%npar_hw(1) = saved_gheads(GV%CurSnapNum,fnum)%npar_hw(1)
       
    ghead%flag_sfr      = saved_gheads(GV%CurSnapNum,fnum)%flag_sfr
    ghead%flag_feedback = saved_gheads(GV%CurSnapNum,fnum)%flag_feedback
    ghead%flag_cooling  = saved_gheads(GV%CurSnapNum,fnum)%flag_cooling
    ghead%flag_age      = saved_gheads(GV%CurSnapNum,fnum)%flag_age
    ghead%flag_metals   = saved_gheads(GV%CurSnapNum,fnum)%flag_metals
    ghead%flag_entr_ics = saved_gheads(GV%CurSnapNum,fnum)%flag_entr_ics

    
    ghead%OmegaB = GV%OmegaB
    ghead%rays_traced = GV%src_rayn

#ifdef incHe
    ghead%flag_helium = 1
#else
    ghead%flag_helium = 0
#endif

#ifdef outGamma
    ghead%flag_gamma = 1
#else
    ghead%flag_gamma = 0
#endif

    ghead%unused = 0

  end subroutine config_to_ghead

!> outputs the whole shebang to a file
!=======================================

  subroutine output_total_snap(pars)
  use particle_system_mod, only: scale_comoving_to_physical  
  use particle_system_mod, only: scale_physical_to_comoving  

     character(clen), parameter :: myname="output_total_snap"
     logical, parameter :: crash = .true.

     type(particle_type), intent(inout) :: pars(:)     !< particles 
     type(gadget_header_type) :: ghead

     character(3) :: label
     character(4) :: ext
     character(clen) :: filename
     integer(i8b) :: lun, Nread, Nfile
     integer(i4b) :: i, ifile

     real(r8b) :: nHe_over_nH
     real(r4b), allocatable :: mu(:), uint(:)
     real(r4b), allocatable :: rblock3(:,:)


     write(*,*) 'output number: ', GV%OutputIndx
     write(*,*) 'output type:   ', GV%OutputType
     write(*,*) "writing total state of system"
     write(*,*) "time (elapsed code) ", GV%time_elapsed_code
     write(*,*) "time (elapsed myr)  ", GV%time_elapsed_s * s2Myr

     if (GV%Comoving) call scale_physical_to_comoving(PLAN%snap(GV%CurSnapNum)%ScalefacAt, pars)

     100 format(I3.3)
     write(label,100) GV%OutputIndx


     ! GADGET formatted output
     !================================================================
     !================================================================
     Nread = 0
        
     do ifile = 0, GV%ParFilesPerSnap-1
           
        call config_to_ghead( ifile, ghead )

        !    form file name
        write(ext,'(I3)') ifile

        filename = trim(GV%OutputDir) // "/" // trim(GV%OutputFileBase) // "_" // label 
        if (GV%ParFilesPerSnap > 1) then
           filename = trim(filename) // "." // trim(adjustl(ext))
        end if
        write(*,*) "writing snapshot state to ", trim(filename)
        write(*,*) 
        call open_unformatted_file_w(filename,lun)
        
        Nfile = ghead%npar_file(1)

        write(lun) ghead     

        allocate( rblock3(3,Nfile) )
        rblock3(1,:) = pars(Nread+1:Nread+Nfile)%pos(1) 
        rblock3(2,:) = pars(Nread+1:Nread+Nfile)%pos(2) 
        rblock3(3,:) = pars(Nread+1:Nread+Nfile)%pos(3) 
        write(lun) rblock3

        rblock3(1,:) = pars(Nread+1:Nread+Nfile)%vel(1) 
        rblock3(2,:) = pars(Nread+1:Nread+Nfile)%vel(2) 
        rblock3(3,:) = pars(Nread+1:Nread+Nfile)%vel(3) 
        write(lun) rblock3
        deallocate( rblock3 )
        
        write(lun) pars(Nread+1:Nread+Nfile)%id 
        write(lun) pars(Nread+1:Nread+Nfile)%mass


        ! calculate some quantities for the next round of outputs
        allocate( mu(Nfile), uint(Nfile) )

        pars(:)%ye = pars(:)%xHII 
#ifdef incHe
        nHe_over_nH = 0.25d0 * GV%He_mf / GV%H_mf
        pars(:)%ye = pars(:)%ye + ( pars(:)%xHeII + 2.0d0 * pars(:)%xHeIII ) * nHe_over_nH
#endif

        mu = 4.0d0 / (3.0d0 * GV%H_mf + 1.0d0 + 4.0d0 * GV%H_mf * pars(Nread+1:Nread+Nfile)%ye)
        uint = ( k_erg_K * pars(Nread+1:Nread+Nfile)%T ) / ( mu * M_p * 2.0d0/3.0d0 * GV%cgs_enrg / GV%cgs_mass) 
        
        write(lun) uint(1:Nfile)
        write(lun) pars(Nread+1:Nread+Nfile)%rho         
        write(lun) pars(Nread+1:Nread+Nfile)%ye
        write(lun) pars(Nread+1:Nread+Nfile)%xHI
        write(lun) pars(Nread+1:Nread+Nfile)%hsml 

        write(lun) pars(Nread+1:Nread+Nfile)%T
#ifdef incHe
        write(lun) pars(Nread+1:Nread+Nfile)%xHeI
        write(lun) pars(Nread+1:Nread+Nfile)%xHeII
#endif
#ifdef outGamma
        do i = Nread+1, Nread+Nfile
           if (pars(i)%time > 0.0) then
              pars(i)%gammaHI = pars(i)%gammaHI / pars(i)%time
           else
              pars(i)%gammaHI = 0.0
           end if
        end do
        write(lun) pars(Nread+1:Nread+Nfile)%gammaHI
        
        pars(Nread+1:Nread+Nfile)%gammaHI = 0.0
        pars(Nread+1:Nread+Nfile)%time = 0.0
#endif

        write(lun) pars(Nread+1:Nread+Nfile)%lasthit


        close(lun)
        deallocate( mu, uint )
        Nread = Nread + Nfile

     end do


     !================================================================
     !================================================================
     

     if (GV%Comoving) call scale_comoving_to_physical(PLAN%snap(GV%CurSnapNum)%ScalefacAt,pars)

  end subroutine output_total_snap



  
!> writes the global ionization state of the system to a file
!=============================================================
  subroutine ion_frac_out(psys,tree)
  use iliev_comparison_project_mod

     type(particle_system_type), intent(in) :: psys    !< particle system
     type(oct_tree_type), intent(in) :: tree

     real :: Nionfrac, Mionfrac, Vionfrac
     real :: CallsPerCross

     100 format(A,T24,A,T44,A,T64,A)
     101 format(A,T44,A,T64,A)
     105 format(T21,ES12.3,T43,ES12.3,T63,ES12.3)
     106 format(T43,ES12.3,T63,ES12.3)
     110 format(T21,I15,T43,I15,T63,ES11.5)

     write(*,100) "time:", "code units", "Myrs", "seconds"
     write(*,105) GV%time_code, &
                  GV%time_s * s2Myr, &
                  GV%time_s
     write(*,*) 

     write(*,100) "time elspsed:", "code units", "Myrs", "seconds"
     write(*,105) GV%time_elapsed_code, &
                  GV%time_elapsed_s * s2Myr, &
                  GV%time_elapsed_s
     write(*,*) 


!    calculate the ionized number, volume, and mass fractions
     Nionfrac = number_weight_ionfrac(psys%par)
     Mionfrac = mass_weight_ionfrac(psys%par) 
     Vionfrac = volume_weight_ionfrac(psys%par)

     GV%nwionfrac = Nionfrac
     GV%mwionfrac = Mionfrac
     GV%vwionfrac = Vionfrac


     write(*,100) "neutral fraction:", "mass weighted"
     write(*,105) 1.0d0 - Mionfrac
     write(*,*) 

     150 format (4ES15.5)
     write(GV%ionlun,150) GV%time_elapsed_code, &
                          GV%time_elapsed_s * s2Myr, &
                          Mionfrac

     write(*,100) "rays cast:", "source", "diffuse", "diffuse/source"
     write(*,105) GV%TotalSourceRaysCast, GV%TotalDiffuseRaysCast, &
                  GV%TotalDiffuseRaysCast/GV%TotalSourceRaysCast 
     write(*,*) 
     
     if (GV%TotalPhotonsCast > 0.) then
        write(*,100) "photons:", "total cast", "per second", "% leaving box"
        write(*,105) GV%TotalPhotonsCast, GV%IonizingPhotonsPerSec, &
                     100.0 * GV%PhotonsLeavingBox / GV%TotalPhotonsCast
     else
        write(*,*) "zero luminosity source test"
     end if
     write(*,*)

     if (GV%ParticleCrossings .GT. 0) then
        CallsPerCross = GV%TotalDerivativeCalls / GV%ParticleCrossings
     else
        CallsPerCross = 0.0
     end if
     write(*,101) "solver(since last screen output):", &
                  "par crossings", "evals/crossings"
     write(*,106) GV%ParticleCrossings, CallsPerCross
     GV%TotalDerivativeCalls = 0.0
     GV%ParticleCrossings = 0.0

                  
     write(*,*)
     write(*,*) "Peak Updates = ", GV%PeakUpdates
     write(*,*) 



     160 format(T1,A,T20,ES12.5,T34,ES12.5)
     161 format(T1,A,T20,ES12.5)
     write(*,160) "Min/Max xHI     = ", minval(psys%par(:)%xHI), &
                                        maxval(psys%par(:)%xHI)

     write(*,160) "Min/Max xHII    = ", minval(psys%par(:)%xHII), &
                                        maxval(psys%par(:)%xHII)

#ifdef incHe
     write(*,160) "Min/Max xHeI    = ", minval(psys%par(:)%xHeI), &
                                        maxval(psys%par(:)%xHeI)
     write(*,160) "Min/Max xHeII   = ", minval(psys%par(:)%xHeII), &
                                        maxval(psys%par(:)%xHeII)
     write(*,160) "Min/Max xHeIII  = ", minval(psys%par(:)%xHeIII), &
                                        maxval(psys%par(:)%xHeIII)
#endif

     write(*,160) "Min/Max T       = ", minval(psys%par(:)%T), &
                                        maxval(psys%par(:)%T)




! do test specific outputs
     if (GV%DoTestScenario) then
        if (GV%TestScenario=="iliev_test1") then
           call iliev_test1_screen_out(psys,tree,GV)
        else if (GV%TestScenario=="iliev_test2") then
           call iliev_test2_screen_out(psys,tree,GV)
        else if (GV%TestScenario=="iliev_test3") then
           call iliev_test3_screen_out(psys,GV)
        else if (GV%TestScenario=="iliev_test4") then
           call iliev_test4_screen_out(psys,GV)
        end if
     end if

     write(*,*) 

  end subroutine ion_frac_out

!> calculates the number weighted ionization fraction
!========================================================

  function number_weight_ionfrac(pars) result(numionfrac)
     type(particle_type), intent(in) :: pars(:)   !< particle system
     real(r8b) :: numionfrac                      !< ion fraction n weighted

     integer :: i

!      because the ionization fractions are stored as single
!      precision variables, the intrinsic sum function returns 
!      a single precision answer.  I think it is better to use 
!      the loop here.

       numionfrac = 0.0d0
       do i = 1,size(pars)
          numionfrac = numionfrac + pars(i)%xHII 
       end do
       numionfrac = numionfrac / size(pars)

  end function number_weight_ionfrac

!> calculates the mass weighted ionization fraction
!========================================================

  function mass_weight_ionfrac(pars) result(massionfrac)
     type(particle_type), intent(in) :: pars(:) !< particle system
     real(r8b) :: massionfrac                   !< ion fraction m weighted
     real(r8b) :: masstot                       !< total volume
     integer :: i

     massionfrac = 0.0d0
     masstot = 0.0d0
     do i = 1,size(pars)
        massionfrac = massionfrac + pars(i)%mass * pars(i)%xHII
        masstot = masstot + pars(i)%mass
     end do
     massionfrac = massionfrac / masstot


  end function mass_weight_ionfrac


!> calculates the volume weighted ionization fraction
!========================================================

  function volume_weight_ionfrac(pars) result(volionfrac)
     type(particle_type), intent(in) :: pars(:) !< particle system
     real(r8b) :: volionfrac                    !< ion fraction v weighted
     real(r8b) :: voltot                        !< total volume
     real(r8b) :: h3                            !< hsml^3
     integer :: i

     volionfrac = 0.0d0
     voltot = 0.0d0
     do i = 1,size(pars)
        h3 = pars(i)%hsml * pars(i)%hsml * pars(i)%hsml
        volionfrac = volionfrac + h3 * pars(i)%xHII
        voltot = voltot + h3
     end do
     volionfrac = volionfrac / voltot

   end function volume_weight_ionfrac


end module output_mod
