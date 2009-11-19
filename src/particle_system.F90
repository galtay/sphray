!> \file particle_system.F90

!> \brief particles, sources, and box, types and subroutines
!! 
!<

!> particle, source, and box descriptions
module particle_system_mod
use myf90_mod
implicit none
private

public :: particle_type
public :: source_type
public :: box_type
public :: particle_system_type
public :: transformation_type
public :: orderparticles
public :: orderpsys
public :: calc_bytes_per_particle
public :: scale_comoving_to_physical
public :: scale_physical_to_comoving
public :: enforce_x_and_T_minmax
public :: particle_info_to_screen
public :: adjustbox

!> basic particle type. 
!=========================
type particle_type

   real(r4b)    :: pos(3)    !< x,y,z coordinates

#ifdef incVel
   real(r4b)    :: vel(3)    !< x,y,z velocities
#endif

   integer(i4b) :: id        !< particle id
   real(r4b)    :: mass      !< particle mass
   real(r4b)    :: T         !< temperature in K       
   real(r4b)    :: rho       !< density 
   real(r4b)    :: ye        !< electron fraction
   real(r4b)    :: xHI       !< nHI/nH 
   real(r4b)    :: xHII      !< * nHII/nH
   real(r4b)    :: hsml      !< smoothing length

#ifdef incHmf   
   real(r4b)    :: Hmf       !< Hydrogen mass fraction
#endif

#ifdef incHe
   real(r4b)    :: xHeI      !< * nHeI/nHe 
   real(r4b)    :: xHeII     !< * nHeII/nHe
   real(r4b)    :: xHeIII    !< * nHeIII/nHe
#endif

#ifdef incHemf
   real(r4b)    :: Hemf      !< Helium mass fraction
#endif

#ifdef outGamma
   real(r4b)    :: gammaHI   !< * time averaged HI photoionization rate
   real(r4b)    :: time      !< * elapsed time in seconds - reset at outputs
#endif

   integer(i8b) :: lasthit   !< * indx of last ray to cross this particle
  
end type particle_type



!> source type
!================
type source_type
   real(r4b)    :: pos(3)    !< x,y,z coordinates

#ifdef incVel
   real(r4b)    :: vel(3)    !< x,y,z velocities
#endif

   real(r4b)    :: L         !< luminosity
   real(r4b)    :: SpcType   !< spectral type
   integer(i4b) :: EmisPrf   !< emission profile
   real(r4b)    :: Lcdf      !< relates this luminosity to the other sources
   integer(i8b) :: lastemit  !< last ray emitted from this source
end type source_type


!> simulation box and boundary conditions
!==========================================
type box_type
   real(r8b)    :: top(1:3)     !< upper x,y,z coordinates
   real(r8b)    :: bot(1:3)     !< lower x,y,z coordinates
   integer(i8b) :: bbound(1:3)  !< BCs for upper faces (0:vac 1:per -1: ref) 
   integer(i8b) :: tbound(1:3)  !< BCs for lower faces (0:vac 1:per -1: ref) 
end type box_type


!----------------------------
!> particles, sources, and box
  type particle_system_type
     type(box_type) :: box                           !< the simulation box     
     type(particle_type), allocatable :: par(:)      !< all particles
     type(source_type), allocatable :: src(:)        !< all sources
  end type particle_system_type


!> transformation type
!=======================
type transformation_type
   integer(i8b) :: fac(1:3)  !< newpos = fac * (oldpos - shift)
   real(i8b) :: shift(1:3)   !< newpos = fac * (oldpos - shift)
end type transformation_type



contains


! this routine rearranges the particles in pars so that they are stored
! in the sequence given by the array order.  ie order = [3,1,2] takes the 
! third particle to the first position, the first particle to the second
! position, and the second particle to the third position.  the array 
! order is not preserved during the routine
!> reorders the particles according to the array order
!===========================================================================
subroutine orderparticles(pars,order)
  type(particle_type), intent(inout) :: pars(:) !< input particles
  integer(i8b), intent(inout) :: order(:)       !< desired order

  type(particle_type) :: par
  integer(i8b) :: i
  integer(i8b) :: goal
  integer(i8b) :: npar
  
  if (size(pars) /= size(order)) stop "size(pars) /= size(order)"
  npar = size(pars)
  
  do i=1,npar 
     par=pars(i)
     goal=order(i)
     do while(goal < i)
        goal=order(goal)
        order(i)=goal
     enddo
     pars(i)=pars(goal)
     pars(goal)=par 
  enddo
  do i=1,npar
     order(i)=i
  enddo

end subroutine orderparticles

!-------------------------------------------------------------------------- 
! this routine rearranges the particles in a particle system  according to 
! the array order as in the routine above
!> reorders the particles in a particle system
    subroutine orderpsys(psys,order)
      type(particle_system_type), intent(inout) :: psys  !< inout particle sys
      integer(i8b), intent(inout) :: order(1:size(psys%par))  !< desired order

      type(particle_type) :: par
      integer(i8b) :: i
      integer(i8b) :: goal

        do i=1,size(psys%par)
           par=psys%par(i)
           goal=order(i)
           do while(goal < i)
              goal=order(goal)
              order(i)=goal
           enddo
           psys%par(i)=psys%par(goal)
           psys%par(goal)=par 
        enddo
        do i=1,size(psys%par)
           order(i)=i
        enddo

      end subroutine orderpsys

!> creates a transformed particle from an input particle
!============================================================
subroutine copypart(parclone,par,transform)
  type(particle_type), intent(out) :: parclone  !< particle copy
  type(particle_type), intent(in)  :: par       !< input particle
  type(transformation_type), intent(in), optional :: transform !< optional transformation
  parclone=par
  if(present(transform)) parclone%pos = transform%fac * (par%pos - transform%shift)
end subroutine copypart


!> resets the box limits where the BCs are vacuum
!==================================================================
subroutine adjustbox(box,bot,top)
  type(box_type), intent(inout) :: box !< input box
  real(r4b), intent(in) :: bot(3)           !< new bottoms
  real(r4b), intent(in) :: top(3)           !< new tops

  where (box%bbound==0) box%bot = bot
  where (box%tbound==0) box%top = top

end subroutine adjustbox


!> scales particles, sources, and the box from comoving to physical values
!> (the velocity is taken from comoving to peculiar)
!==========================================================================
subroutine scale_comoving_to_physical(a,par,src,box)

  real(r8b), intent(in) :: a                           !< scale factor
  type(particle_type), intent(inout) :: par(:)         !< particles
  type(source_type), optional, intent(inout) :: src(:) !< sources
  type(box_type), optional, intent(inout) :: box       !< box 

  write(*,'(A,F12.4)') "scaling comoving to physical coordinates, a = ", a

  par%pos(1) = par%pos(1) * a
  par%pos(2) = par%pos(2) * a
  par%pos(3) = par%pos(3) * a

#ifdef incVel
  par%vel(1) = par%vel(1) * a
  par%vel(2) = par%vel(2) * a
  par%vel(3) = par%vel(3) * a
#endif

  par%hsml = par%hsml * a
  par%rho  = par%rho / (a*a*a)

  if (present(src)) then  
     src%pos(1) = src%pos(1) * a
     src%pos(2) = src%pos(2) * a
     src%pos(3) = src%pos(3) * a

#ifdef incVel
     src%vel(1) = src%vel(1) * a
     src%vel(2) = src%vel(2) * a
     src%vel(3) = src%vel(3) * a
#endif

  end if

  if (present(box)) then
     box%top = box%top * a
     box%bot = box%bot * a
  end if

end subroutine scale_comoving_to_physical


!> scales particles, sources, and the box from physical to comoving  values
!> (the velocity is taken from peculiar to comoving)
!==========================================================================
subroutine scale_physical_to_comoving(a,par,src,box)

  real(r8b), intent(in) :: a                                 !< scale factor 
  type(particle_type), intent(inout) :: par(:)          !< particles
  type(source_type), optional, intent(inout) :: src(:)  !< sources
  type(box_type), optional, intent(inout) :: box        !< box 

  write(*,'(A,F8.4)') "scaling physical to comoving coordinates, a = ", a
  
  par%pos(1) = par%pos(1) / a
  par%pos(2) = par%pos(2) / a
  par%pos(3) = par%pos(3) / a

#ifdef incVel
  par%vel(1) = par%vel(1) / a
  par%vel(2) = par%vel(2) / a
  par%vel(3) = par%vel(3) / a
#endif

  par%hsml = par%hsml / a
  par%rho = par%rho * (a*a*a)

  if (present(src)) then
     src%pos(1) = src%pos(1) / a
     src%pos(2) = src%pos(2) / a
     src%pos(3) = src%pos(3) / a

#ifdef incVel
     src%vel(1) = src%vel(1) / a
     src%vel(2) = src%vel(2) / a
     src%vel(3) = src%vel(3) / a
#endif

  end if

  if (present(box)) then
     box%top = box%top / a
     box%bot = box%bot / a
  end if

end subroutine scale_physical_to_comoving





!>   enforces a minimum and maximum value on the ionization fractions and temperatures
!=======================================================================================
subroutine enforce_x_and_T_minmax(par,xmin,xmax,tmin,tmax)

  type(particle_type), intent(inout) :: par(:)  !< inout particle system
  real(r8b), intent(in) :: xmin, xmax, tmin, tmax
  integer(i8b) :: i

  do i = 1,size(par)
     if (par(i)%xHI < xmin) par(i)%xHI = xmin
     if (par(i)%xHI > xmax) par(i)%xHI = xmax
#ifdef incHe
     if (par(i)%xHeI  < xmin) par(i)%xHeI  = xmin
     if (par(i)%xHeII < xmin) par(i)%xHeII = xmin
     if (par(i)%xHeI  > xmax) par(i)%xHeI  = xmax
     if (par(i)%xHeII > xmax) par(i)%xHeII = xmax
#endif

     if (par(i)%T < tmin) par(i)%T = tmin
     if (par(i)%T > tmax) par(i)%T = tmax

  end do

end subroutine enforce_x_and_T_minmax


!> figures out how many bytes of RAM are needed per particle
!=================================================================
subroutine calc_bytes_per_particle(bpp, bps)
  integer(i8b), intent(out) :: bpp !< bytes per particle
  integer(i8b), intent(out) :: bps !< bytes per source
  integer(i8b) :: bytes  
  integer, parameter :: verb=2

  character(clen) :: str

  bytes = 12         ! positions

#ifdef incVel
  bytes = bytes + 12 ! velocities  
#endif

  bytes = bytes + 4  ! ID
  bytes = bytes + 4  ! mass
  bytes = bytes + 4  ! temperature
  bytes = bytes + 4  ! rho
  bytes = bytes + 4  ! ye
  bytes = bytes + 8  ! H ionization fractions
  bytes = bytes + 4  ! hsml

#ifdef incHe
  bytes = bytes + 12 ! He ionization fractions
#endif

#ifdef outGamma
  bytes = bytes + 8  ! GammaHI tracking
#endif

  bytes = bytes + 8  ! last hit index


  write(str,'(A,I4)') "bytes per particle = ", bytes
  call mywrite(str, verb)
  bpp = bytes
  

 
  bps = 10 * 4 + 8
  write(*,'(A,I3)') "bytes per src = ", bps


 end subroutine calc_bytes_per_particle


!> outputs currently loaded particle data to the screen
!=================================================================
subroutine particle_info_to_screen(par)

  type(particle_type), intent(in) :: par(:)  !< particle array

99 format(72("="))  
100 format(A,T10,3ES15.5)
101 format(A,T10,2I15,ES15.5)
102 format(A,T10,2I15)
  write(*,99) 
  write(*,*) " particle data  "
  write(*,99) 
  write(*,*) 


     write(*,100) "xpos", minval(par%pos(1)), maxval(par%pos(1)), &
                          meanval_real(par%pos(1))

     write(*,100) "ypos", minval(par%pos(2)), maxval(par%pos(2)), &
                          meanval_real(par%pos(2))

     write(*,100) "zpos", minval(par%pos(3)), maxval(par%pos(3)), &
                          meanval_real(par%pos(3))   

#ifdef incVel
     write(*,100) "xvel", minval(par%vel(1)), maxval(par%vel(1)), &
                          meanval_real(par%vel(1))

     write(*,100) "yvel", minval(par%vel(2)), maxval(par%vel(2)), &
                          meanval_real(par%vel(2))

     write(*,100) "zvel", minval(par%vel(3)), maxval(par%vel(3)), &
                          meanval_real(par%vel(3))
#endif

     write(*,102) "id",   minval(par%id), maxval(par%id)

     write(*,100) "mass", minval(par%mass), maxval(par%mass), &
                          meanval_real(par%mass)

     write(*,100) "T",    minval(par%T), maxval(par%T), &
                          meanval_real(par%T)

     write(*,100) "rho",  minval(par%rho), maxval(par%rho), &
                          meanval_real(par%rho)

     write(*,100) "ye",    minval(par%ye), maxval(par%ye), &
                            meanval_real(par%ye)

     write(*,100) "xHI",    minval(par%xHI), maxval(par%xHI), &
                            meanval_real(par%xHI)

     write(*,100) "xHII",   minval(par%xHII), maxval(par%xHII), &
                            meanval_real(par%xHII)

     write(*,100) "hsml", minval(par%hsml), maxval(par%hsml), &
                          meanval_real(par%hsml)

#ifdef incHe
     write(*,100) "xHeI",   minval(par%xHeI), maxval(par%xHeI), &
                            meanval_real(par%xHeI)

     write(*,100) "xHeII",  minval(par%xHeII), maxval(par%xHeII), &
                            meanval_real(par%xHeII)

     write(*,100) "xHeIII", minval(par%xHeIII), maxval(par%xHeIII), &
                            meanval_real(par%xHeIII)
#endif

#ifdef outGamma
     write(*,100) "gammaHI",   minval(par%gammaHI), maxval(par%gammaHI), &
                            meanval_real(par%gammaHI)

     write(*,100) "time(s)",   minval(par%time), maxval(par%time), &
                            meanval_real(par%time)
#endif

     write(*,101) "lasthit", minval(par%lasthit), maxval(par%lasthit)
  
  write(*,99) 

end subroutine particle_info_to_screen

!> calculates the mean value w/o using the intrinsics
!=================================================================
function meanval_real(arr) result (mean)
  real(r4b), dimension(:), intent(in) :: arr  !< array to average
  real(r8b) :: mean                           !< mean value to return
  integer(i8b) :: i
  
  mean = 0.0d0
  do i = 1,size(arr)
     mean = mean + arr(i)
  end do
  mean = mean / size(arr)

end function meanval_real





end module particle_system_mod
