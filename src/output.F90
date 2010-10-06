!> \file output.F90

!> \brief The module that handles output
!<

module output_mod
use myf03_mod
use gadget_general_class
use gadget_public_header_class
use gadget_public_input_hdf5_mod
use gadget_sphray_header_class
use particle_system_mod, only: particle_system_type
use particle_system_mod, only: particle_type
use particle_system_mod, only: set_ye_pars
use particle_system_mod, only: number_weight_ionfrac
use particle_system_mod, only: mass_weight_ionfrac
use particle_system_mod, only: volume_weight_ionfrac
use oct_tree_mod, only: oct_tree_type
use physical_constants_mod
use global_mod, only: PLAN, GV
use global_mod, only: saved_gheads
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
    type(gadget_sphray_header_type), intent(out) :: ghead !< gadget header

    ghead%npar_all = 0
    ghead%npar_all(0) = saved_gheads(GV%CurSnapNum,fnum)%npar_all(0)     

    ghead%npar_file = 0
    ghead%npar_file(0) = saved_gheads(GV%CurSnapNum,fnum)%npar_file(0)     

    ghead%mass = 0.0

    ghead%nfiles    = saved_gheads(GV%CurSnapNum,fnum)%nfiles

    if (GV%Comoving) then
       ghead%a = PLAN%snap(GV%CurSnapNum)%ScalefacAt
       ghead%z = 1.0d0 / ghead%a - 1.0d0
    else
       ghead%a = GV%time_elapsed_myr
       ghead%z = saved_gheads(GV%CurSnapNum,fnum)%z
    end if

    ghead%boxlen   = saved_gheads(GV%CurSnapNum,fnum)%boxlen
    ghead%OmegaM   = saved_gheads(GV%CurSnapNum,fnum)%OmegaM
    ghead%OmegaL   = saved_gheads(GV%CurSnapNum,fnum)%OmegaL
    ghead%OmegaB   = saved_gheads(GV%CurSnapNum,fnum)%OmegaB
    ghead%time_gyr = saved_gheads(GV%CurSnapNum,fnum)%time_gyr
    ghead%h        = saved_gheads(GV%CurSnapNum,fnum)%h

    ghead%npar_hw = 0
    ghead%npar_hw(1) = saved_gheads(GV%CurSnapNum,fnum)%npar_hw(1)
       
    ghead%flag_sfr      = saved_gheads(GV%CurSnapNum,fnum)%flag_sfr
    ghead%flag_feedback = saved_gheads(GV%CurSnapNum,fnum)%flag_feedback
    ghead%flag_cooling  = saved_gheads(GV%CurSnapNum,fnum)%flag_cooling
    ghead%flag_age      = saved_gheads(GV%CurSnapNum,fnum)%flag_age
    ghead%flag_metals   = saved_gheads(GV%CurSnapNum,fnum)%flag_metals
    ghead%flag_entr_ics = saved_gheads(GV%CurSnapNum,fnum)%flag_entr_ics
    

    ghead%rays_traced = GV%src_rayn

#ifdef incHmf
    ghead%flag_Hmf = 1
#else
    ghead%flag_Hmf = 0
#endif
!-------------------------
#ifdef incHemf
    ghead%flag_Hemf = 1
#else
    ghead%flag_Hemf = 0
#endif
!-------------------------
#ifdef incHe
    ghead%flag_helium = 1
#else
    ghead%flag_helium = 0
#endif
!-------------------------
#ifdef outGammaHI
    ghead%flag_gammaHI = 1
#else
    ghead%flag_gammaHI = 0
#endif
!-------------------------
#ifdef cloudy
    ghead%flag_cloudy = 1
#else
    ghead%flag_cloudy = 0
#endif
!-------------------------
#ifdef incEOS
    ghead%flag_eos = 1
#else
    ghead%flag_eos = 0
#endif
!-------------------------
#ifdef incSFR
    ghead%flag_incsfr = 1
#else
    ghead%flag_incsfr = 0
#endif
!-------------------------

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
     type(gadget_sphray_header_type) :: ghead
     type(gadget_constants_type) :: gconst

     character(3) :: label
     character(4) :: ext
     character(clen) :: filename
     integer(i4b) :: lun
     integer(i8b) :: Nread, Nfile
     integer(i4b) :: ipar, ifile

     integer(i8b) :: i,j

     real(r8b) :: scale
     real(r8b) :: hub
     real(r8b) :: nHe_over_nH
     real(r8b) :: Hmf, Hemf
     real(r4b) :: mu
     real(r4b), allocatable :: uint(:)
     real(r4b), allocatable :: rblock3(:,:)

     ! set defaults
     !==============
     nHe_over_nH = 0.0
     ipar=0

     scale = PLAN%snap(GV%CurSnapNum)%ScalefacAt
     hub   = GV%LittleH

     write(*,*) 'output number: ', GV%OutputIndx
     write(*,*) 'output type:   ', GV%OutputType
     write(*,*) "writing total state of system"
     write(*,*) "time (elapsed code) ", GV%time_elapsed_code
     write(*,*) "time (elapsed myr)  ", GV%time_elapsed_myr

     if (GV%Comoving) call scale_physical_to_comoving(scale, hub, pars)

     100 format(I3.3)
     write(label,100) GV%OutputIndx


     call set_ye_pars( pars, GV%H_mf, GV%He_mf, GV%NeBackground )


     !================================================================
     ! GADGET standard formatted output
     !================================================================
     !================================================================
     if (GV%OutputType == 1) then
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

           call ghead%write_Gsphray_header_lun(lun)

           allocate( rblock3(3,Nfile) )
           rblock3(1,:) = pars(Nread+1:Nread+Nfile)%pos(1) 
           rblock3(2,:) = pars(Nread+1:Nread+Nfile)%pos(2) 
           rblock3(3,:) = pars(Nread+1:Nread+Nfile)%pos(3) 
           write(lun) rblock3

#ifdef incVel
           rblock3(1,:) = pars(Nread+1:Nread+Nfile)%vel(1) 
           rblock3(2,:) = pars(Nread+1:Nread+Nfile)%vel(2) 
           rblock3(3,:) = pars(Nread+1:Nread+Nfile)%vel(3) 
#else
           rblock3 = 0.0e0
#endif

           write(lun) rblock3
           deallocate( rblock3 )

           write(lun) pars(Nread+1:Nread+Nfile)%id 
           write(lun) pars(Nread+1:Nread+Nfile)%mass


           ! do some conversions
           !------------------------
           allocate( uint(Nfile) )

           j=1
           do i = Nread+1, Nread+Nfile
!------------------------------
#ifdef incHmf
              Hmf=pars(i)%Hmf
#else
              Hmf=GV%H_mf
#endif     
!-------------------------------
              mu      = 4.0d0 / ( 3.0d0 * Hmf + 1.0d0 + 4.0d0 * Hmf * pars(i)%ye )
              uint(j) = ( gconst%BOLTZMANN * pars(i)%T ) / &
                        ( (gconst%GAMMA - 1.0d0) * mu * gconst%PROTONMASS  )
              uint(j) = uint(j) * GV%cgs_mass / GV%cgs_enrg

              j = j+1
           end do


           write(lun) uint(1:Nfile)
           deallocate( uint )

           write(lun) pars(Nread+1:Nread+Nfile)%rho         
           write(lun) pars(Nread+1:Nread+Nfile)%ye
           write(lun) pars(Nread+1:Nread+Nfile)%xHI
           write(lun) pars(Nread+1:Nread+Nfile)%hsml 

           write(lun) pars(Nread+1:Nread+Nfile)%T


#ifdef incHmf
           write(lun) pars(Nread+1:Nread+Nfile)%Hmf
#endif


#ifdef incHemf
           write(lun) pars(Nread+1:Nread+Nfile)%Hemf
#endif


#ifdef incHe
           write(lun) pars(Nread+1:Nread+Nfile)%xHeI
           write(lun) pars(Nread+1:Nread+Nfile)%xHeII
#endif


#ifdef outGammaHI
           do ipar = Nread+1, Nread+Nfile
              if (pars(ipar)%time > 0.0) then
                 pars(ipar)%gammaHI = pars(ipar)%gammaHI / pars(ipar)%time
              else
                 pars(ipar)%gammaHI = 0.0
              end if
           end do
           write(lun) pars(Nread+1:Nread+Nfile)%gammaHI

           pars(Nread+1:Nread+Nfile)%gammaHI = 0.0
           pars(Nread+1:Nread+Nfile)%time = 0.0
#endif

 
#ifdef cloudy
           write(lun) pars(Nread+1:Nread+Nfile)%xHI_cloudy
#endif

#ifdef incEOS
           write(lun) pars(Nread+1:Nread+Nfile)%eos
#endif

#ifdef incSFR
           write(lun) pars(Nread+1:Nread+Nfile)%sfr
#endif

           write(lun) pars(Nread+1:Nread+Nfile)%lasthit


           close(lun)
           Nread = Nread + Nfile

        end do

     end if
 

     !================================================================     
     ! GADGET HDF5 formatted output
     !================================================================
     !================================================================
     if (GV%OutputType == 2) then
        


     end if
 
     !================================================================
     !================================================================


     if (GV%Comoving) call scale_comoving_to_physical(scale, hub, pars)

  end subroutine output_total_snap



  
!> writes the global ionization state of the system to a file
!=============================================================
  subroutine ion_frac_out(psys,tree)
  use iliev_comparison_project_mod

     type(particle_system_type), intent(in) :: psys    !< particle system
     type(oct_tree_type), intent(in) :: tree

     real :: Nionfrac, Mionfrac, Vionfrac
     real :: CallsPerCross
     real :: minGHI, maxGHI, GHI
     integer :: i
     integer :: nothit

     100 format(A,T24,A,T44,A,T64,A)
     101 format(A,T44,A,T64,A)
     105 format(T21,ES12.3,T43,ES12.3,T63,ES12.3)
     106 format(T43,ES12.3,T63,ES12.3)
     110 format(T21,I15,T43,I15,T63,ES11.5)

     write(*,100) "time:", "code units", "Myrs", "seconds"
     write(*,105) GV%start_time_code + GV%itime * GV%dt_code, &
                  GV%start_time_myr  + GV%itime * GV%dt_myr, &
                  GV%start_time_s    + GV%itime * GV%dt_s
     write(*,*) 

     write(*,100) "time elapsed:", "code units", "Myrs", "seconds"
     write(*,105) GV%time_elapsed_code, &
                  GV%time_elapsed_myr, &
                  GV%time_elapsed_s
     write(*,*) 


!    calculate the ionized number, volume, and mass fractions
     Nionfrac = number_weight_ionfrac(psys)
     Mionfrac = mass_weight_ionfrac(psys) 
     Vionfrac = volume_weight_ionfrac(psys)

     GV%nwionfrac = Nionfrac
     GV%mwionfrac = Mionfrac
     GV%vwionfrac = Vionfrac


     write(*,100) "neutral fraction:", "number weighted", "mass weighted", "volume weighted"
     write(*,105) 1.0d0-Nionfrac, 1.d0-Mionfrac, 1.0d0-Vionfrac
     write(*,*) 

     150 format (6ES15.5)
     write(GV%ionlun,150) GV%start_time_myr, &
                          GV%time_elapsed_myr, &
                          1.0d0-Nionfrac, 1.0d0-Mionfrac, 1.0d0-Vionfrac

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
     162 format(T1,A,T20,I12,T34,I12)
     write(*,160) "Min/Max xHI      = ", minval(psys%par(:)%xHI), &
                                        maxval(psys%par(:)%xHI)

     write(*,160) "Min/Max xHII     = ", minval(psys%par(:)%xHII), &
                                        maxval(psys%par(:)%xHII)

#ifdef incHe
     write(*,160) "Min/Max xHeI     = ", minval(psys%par(:)%xHeI), &
                                        maxval(psys%par(:)%xHeI)
     write(*,160) "Min/Max xHeII    = ", minval(psys%par(:)%xHeII), &
                                        maxval(psys%par(:)%xHeII)
     write(*,160) "Min/Max xHeIII   = ", minval(psys%par(:)%xHeIII), &
                                        maxval(psys%par(:)%xHeIII)
#endif

     write(*,160) "Min/Max T        = ", minval(psys%par(:)%T), &
                                        maxval(psys%par(:)%T)

#ifdef outGammaHI
     minGHI = huge(1.0)
     maxGHI = tiny(1.0)
     do i = 1, size(psys%par(:))
        if (psys%par(i)%time > 0.0) then
           GHI = psys%par(i)%gammaHI / psys%par(i)%time
           if (GHI < minGHI) minGHI = GHI
           if (GHI > maxGHI) maxGHI = GHI
        endif
     enddo
     write(*,160) "Min/Max GHI      = ", minGHI, maxGHI
#endif

     write(*,162) "Min/Max LastHit  = ", minval(psys%par(:)%lasthit), &
                                        maxval(psys%par(:)%lasthit)


     nothit = 0
     do i = 1, size(psys%par)
        if (psys%par(i)%lasthit == 0) then
           nothit = nothit + 1
        endif
     end do
     write(*,161) "Fraction Not Hit = ", real(nothit) / size(psys%par)


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




end module output_mod
