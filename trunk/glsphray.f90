!> \file glsphray.f90

!> \brief SPHRAY basic operation program + OpenGL viewer
!! 
!! this program does the same thing as sphray.f90 but contains an 
!! extra subroutine to launch and kill the OpenGL viewer.  
!<
program glsphray
use initialize_mod, only: initialize
use mainloop_mod, only: mainloop
implicit none

   character(200) :: config_file !< configuration file

!  get the config file from the runtime argument
   call getarg(1,config_file)
   call system("clear")

   write(*,*) "**************************************************************"
   write(*,*) "**************************************************************"
   write(*,*) "         Smoothed Particle Hydrodynamics Ray Tracer           "
   write(*,*) "                         SPHRAY                               "
   write(*,*) "                  (with OpenGL Viewer)                        "
   write(*,*) "**************************************************************"
   write(*,*) "**************************************************************"
   write(*,*) 
   write(*,*) 
   write(*,*) 
   
   call initialize(config_file)
   
   call viewer()

   call mainloop()

   print*,' *** GAME OVER ***'
   print*,'you can close the viewer safely now..'

   call viewer()


end program glsphray
