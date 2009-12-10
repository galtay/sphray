!> \file skewers.F90

!> \brief Analysis tool for line of sight quantities
!! 
!<

module skewers_mod
use myf90_mod
use main_input_mod, only: readin_snapshot
use oct_tree_mod, only: buildtree, setparticleorder
use raylist_mod, only: prepare_raysearch, kill_raylist, trace_ray
use ray_mod, only: make_skewer_ray
use ray_mod, only: ray_type, curface
use b2cd_mod, only: b2cdfac
use mt19937_mod, only: genrand_real1
use output_mod, only: output_total_snap, ion_frac_out
use physical_constants_mod

! variables
use global_mod, only: psys
use global_mod, only: globalraylist
use global_mod, only: tree
use global_mod, only: GV, PLAN
implicit none

contains


subroutine do_skewers()
  character(clen), parameter :: myname="do_skewers"
  logical, parameter :: crash=.true.
  integer, parameter :: verb=1
  character(clen) :: str,fmt

  !  local counters 
  !-----------------  
  integer(i8b) :: snapn !< snapshot counter
  integer(i8b) :: rayn  !< ray counter
  integer(i8b) :: srcn  !< source counter
  integer(i8b) :: i,j,k
    
  real(r8b) :: rn       !< random number

  integer(i8b) :: pindx
  real(r8b) :: d,b,bnorm,cdfac
  real(r8b) :: hsml,mass,rho,xHI
  real(r8b) :: mass_cgs, rho_cgs, gpercm3, cm3
  real(r8b) :: Hcnt, HIcnt, nH, NHI, NHItot
  
    
  ! work variables
  !-----------------
  type(ray_type) :: ray
  real(r8b) :: MB


  ! loop over the snapshots 
  !=========================
  snaps: do snapn = GV%StartSnapNum, GV%EndSnapNum
     

     open(file="NHI.txt",unit=200)


     call readin_snapshot(skewers=.true.)
     
     ! create custom src for skewering
     !---------------------------------
     if (allocated(psys%src) ) deallocate(psys%src)
     allocate(psys%src(1))
     
     psys%src(1)%pos = (/0.0,0.0,0.0/)
     psys%src(1)%EmisPrf = -3

     !  build oct tree.  only need to do this once per snap (for now)
     !----------------------------------------------------------------
     call buildtree(psys,tree,MB,GV%PartPerCell)
     GV%MB = GV%MB + MB
     call setparticleorder(psys, tree)             
     call prepare_raysearch(psys, globalraylist)


     GV%ForcedRayNumber = 60000

     ! begin ray tracing 
     !------------------------- 
     src_rays: do rayn = 1, GV%ForcedRayNumber

        if (mod(rayn,10000)==0) write(*,*) "rayn = ", rayn
          
          GV%rayn = GV%rayn + 1
          GV%src_rayn = GV%src_rayn + 1
          GV%TotalSourceRaysCast = GV%TotalSourceRaysCast + 1      
          
          curface=0
       
          !  create a source ray and calc the impacts
          call make_skewer_ray(psys%box, ray)


          globalraylist%ray = ray
          call trace_ray(globalraylist%ray, globalraylist, psys, tree) 

          NHItot = 0.0d0
          
          100 format(A,I4,6F8.4, 2E14.5)
          do i = 1,globalraylist%nnb

             d = globalraylist%intersection(i)%d
             b = globalraylist%intersection(i)%b
             pindx = globalraylist%intersection(i)%pindx

             hsml = psys%par(pindx)%hsml
             mass = psys%par(pindx)%mass
             rho  = psys%par(pindx)%rho
             xHI  = psys%par(pindx)%xHI

             bnorm = b / hsml
             cdfac = b2cdfac(bnorm,hsml,GV%cgs_len)

             mass_cgs = mass * GV%cgs_mass 
             rho_cgs  = rho *  GV%cgs_rho 
             
             gpercm3 = rho_cgs
             cm3 = mass_cgs / rho_cgs
             
             nH = rho_cgs * GV%H_mf  / M_H
             Hcnt = mass_cgs * GV%H_mf / M_H
             HIcnt = Hcnt * xHI
             NHI = cdfac * HIcnt

             NHItot = NHItot + NHI

          end do

 
          write(200,"(I,4ES14.4)") curface, NHItot, ray%start

          


       end do src_rays

       close(200)


       stop

  end do snaps
  
  
end subroutine do_skewers


end module skewers_mod






program skewers
  use myf90_mod
  use initialize_mod, only: initialize_skewers
  use skewers_mod
  implicit none

  character(clen) :: config_file    !< configuration file
  type(command_line_type) :: cmnd   !< command line variable

  integer, parameter :: verb = 3    !< verbosity before config file is read
  logical, parameter :: crash = .true.             !< stop on error
  character(clen), parameter :: myname = "skewers"  !< routine name



  ! set initial myf90_verbosity (will be changed once config file is read)
  ! 3 lets any statement with verbosity 3 or less through
  !========================================================================
  myf90_verbosity = verb

  ! initialize command line variable
  !========================================================================
  cmnd = initialize_command_line(verb)

  ! check command line arguments
  !========================================================================
  if (cmnd%nargs /= 1 ) then 
     call myerr( 'usage: ./skewers config_file', myname, crash )
  end if

  ! clear terminal and print splash screen
  !========================================================================
  call system("clear")

  write(*,*) "**************************************************************"
  write(*,*) "**************************************************************"
  write(*,*) "         SKEWERS Line of Sight Analysis                       "
  write(*,*)
  write(*,*) "  developed by,                                               "
  write(*,*)
  write(*,*) "  Gabriel Altay (gabriel.altay@gmail.com)                     "
  write(*,*) "**************************************************************"
  write(*,*) "**************************************************************"
  write(*,*) 
  write(*,*) 

  ! initialize and call mainloop()
  !========================================================================
  config_file = cmnd%args(1)
  call initialize_skewers(config_file)

  call do_skewers()


end program skewers
