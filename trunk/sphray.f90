!> \file sphray.f90

!> \brief SPHRAY basic operation program
!! 
!! This program calls mainloop without initializing the OpenGL viewer.
!<
program sphray
use initialize_mod, only: initialize
use mainloop_mod, only: mainloop
use physical_constants_mod
implicit none

   character(200) :: config_file !< configuration file

   call getarg(1,config_file)
   call system("clear")

   write(*,*) "**************************************************************"
   write(*,*) "**************************************************************"
   write(*,*) "         Smoothed Particle Hydrodynamics Ray Tracer           "
   write(*,*) "                        SPHRAY 1.0                            "
   write(*,*)
   write(*,*) "  developed by,                                               "
   write(*,*)
   write(*,*) "  Gabriel Altay (gabriel.altay@gmail.com) and                 "
   write(*,*) "  Inti Pelupessy                                              "
   write(*,*) "**************************************************************"
   write(*,*) "**************************************************************"
   write(*,*) 
   write(*,*) 

   call initialize(config_file)

   call mainloop()

end program sphray
