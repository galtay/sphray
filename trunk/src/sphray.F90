!> \file sphray.F90

!> \brief SPHRAY basic operation program
!! 
!! This program calls mainloop without initializing the OpenGL viewer.
!<
program sphray
use myf90_mod
use initialize_mod, only: initialize
use mainloop_mod, only: mainloop
implicit none

   character(clen) :: config_file    !< configuration file
   type(command_line_type) :: cmnd   !< command line variable

   integer, parameter :: verb = 3    !< verbosity before config file is read
   logical, parameter :: crash = .true.             !< stop on error
   character(clen), parameter :: myname = "sphray"  !< routine name



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
      call myerr( 'usage: ./sphray config_file', myname, crash )
   end if

   ! clear terminal and print splash screen
   !========================================================================
   write(*,*) 
   write(*,*) "**************************************************************"
   write(*,*) "**************************************************************"
   write(*,*) "         Smoothed Particle Hydrodynamics Ray Tracer           "
   write(*,*) "                        SPHRAY 0.9                            "
   write(*,*)
   write(*,*) "  developed by,                                               "
   write(*,*)
   write(*,*) "  Gabriel Altay (gabriel.altay@gmail.com) and                 "
   write(*,*) "  Inti Pelupessy                                              "
   write(*,*) "**************************************************************"
   write(*,*) "**************************************************************"
   write(*,*) 
   write(*,*) 

   ! initialize and call mainloop()
   !========================================================================
   config_file = cmnd%args(1)
   call initialize(config_file)

   call mainloop()

end program sphray
