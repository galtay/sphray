!> \file sphpar.f90

!> \brief module for determining smoothing lengths of particles
!<

module sphpar_mod
use myf90_mod
use particle_system_mod, only: particle_type
use particle_system_mod, only: transformation_type

!---------------------------------
!> local particle type. 
  type sphpar_type
     integer(i8b) :: id !< particle id
     real(r4b) :: pos(3) !< x,y,z coordinates
     real(r4b) :: vel(3) !< x,y,z velocities
     real(r4b) :: hsmooth !< smoothing length
     real(r4b) :: rho !< density = mass * NsphNnb / hsml^3 
     real(r4b) :: mass !< particle mass 
     real(r4b) :: T !< particle temperature in K
     real(r4b) :: xHI !< HI ionization fraction
     real(r4b) :: xHII !< HII ionization fraction
     real(r4b) :: xHeI !< HeI ionization fraction
     real(r4b) :: xHeII !< HeII ionization fraction
     real(r4b) :: xHeIII !< HeIII ionization fraction
     integer(i8b) :: lasthit !< last ray to cross this particle
     integer(i4b) :: nnb !< number of SPH neighbors
     real(r4b) :: gradrho(3) !< vector gradient of the density
     real(r4b) :: drhodh !< derivative of rho wrt to h
     real(r4b) :: fi  !< used for finding smoothing length
     real(r4b) :: dfi !< used for finding smoothing length
 end type sphpar_type

contains

!-----------------------------------------------------------------------
!> creates a transormed sph particle from an input particle
  subroutine par2sphpar(sphpar,par,pm)
    type(sphpar_type) :: sphpar  !< output sphpar
    type(particle_type) :: par   !< input particle
    type(transformation_type),optional :: pm  !< possible transformation

    sphpar%id = par%id
    sphpar%pos = par%pos
    sphpar%vel = par%vel
    sphpar%hsmooth = par%hsml
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

    if(present(pm)) sphpar%pos=pm%fac*(par%pos-pm%shift)

  end subroutine par2sphpar


end module sphpar_mod
