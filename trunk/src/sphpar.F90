!> \file sphpar.F90

!> \brief module for determining smoothing lengths of particles
!<

module sphpar_mod
use myf90_mod
use particle_system_mod, only: particle_type
use particle_system_mod, only: transformation_type

!---------------------------------
!> local particle type. 
  type, extends(particle_type) :: sphpar_type

#ifndef incVel
     real(r4b) :: vel
#endif

#ifndef incHe
     real(r4b) :: xHeI
     real(r4b) :: xHeII
     real(r4b) :: xHeIII
#endif

     integer(i4b) :: nnb     !< number of SPH neighbors
     real(r4b) :: gradrho(3) !< vector gradient of the density
     real(r4b) :: drhodh     !< derivative of rho wrt to h
     real(r4b) :: fi         !< used for finding smoothing length
     real(r4b) :: dfi        !< used for finding smoothing length
 end type sphpar_type

contains

!-----------------------------------------------------------------------
!> creates a transormed sph particle from an input particle
  subroutine par2sphpar(sphpar,par,pt)
    type(sphpar_type) :: sphpar  !< output sphpar
    type(particle_type) :: par   !< input particle
    type(transformation_type),optional :: pt  !< possible transformation

    sphpar%pos = par%pos

#ifdef incVel
    sphpar%vel = par%vel
#else
    sphpar%vel = 0.0
#endif

    sphpar%id = par%id
    sphpar%hsml = par%hsml
    sphpar%rho = par%rho
    sphpar%mass = par%mass
    sphpar%T = par%T 
    
    sphpar%xHI = par%xHI
    sphpar%xHII = par%xHII

#ifdef incHe
    sphpar%xHeI = par%xHeI
    sphpar%xHeII = par%xHeII
    sphpar%xHeIII = par%xHeIII
#else
    sphpar%xHeI = 0.0
    sphpar%xHeII = 0.0
    sphpar%xHeIII = 0.0
#endif

    sphpar%lasthit = par%lasthit

    if(present(pt)) sphpar%pos=pt%fac*(par%pos-pt%shift)

  end subroutine par2sphpar


end module sphpar_mod
