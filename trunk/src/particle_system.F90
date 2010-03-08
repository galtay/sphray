!> \file particle_system.F90

!> \brief Particles, sources, and box, types and subroutines.
!! 
!<

!> particle, source, and box descriptions
module particle_system_mod
use myf90_mod
use atomic_rates_mod, only: calc_colion_eq_fits

implicit none
private

public :: particle_type
public :: source_type
public :: box_type
public :: particle_system_type
public :: transformation_type
public :: orderparticles
public :: orderpsys
public :: calc_bytes_per_particle_and_source
public :: scale_comoving_to_physical
public :: scale_physical_to_comoving
public :: set_ye
public :: set_ye_pars
public :: set_collisional_ionization_equilibrium
public :: enforce_x_and_T_minmax
public :: particle_info_to_screen
public :: adjustbox

!> Particle type. 
!=========================
type particle_type

   real(r4b)    :: pos(3)     !< x,y,z coordinates

#ifdef incVel
   real(r4b)    :: vel(3)     !< x,y,z velocities
#endif

   integer(i4b) :: id         !< particle id
   real(r4b)    :: mass       !< particle mass
   real(r4b)    :: T          !< temperature in K       
   real(r4b)    :: rho        !< density 
   real(r4b)    :: ye         !< electron fraction
   real(r4b)    :: xHI        !< nHI/nH 
   real(r4b)    :: xHII       !< * nHII/nH
   real(r4b)    :: hsml       !< smoothing length

#ifdef cloudy
   real(r4b)    :: xHI_cloudy !< cloudy eq solutions
#endif

#ifdef incHmf   
   real(r4b)    :: Hmf        !< Hydrogen mass fraction
#endif

#ifdef incHe
   real(r4b)    :: xHeI       !< * nHeI/nHe 
   real(r4b)    :: xHeII      !< * nHeII/nHe
   real(r4b)    :: xHeIII     !< * nHeIII/nHe
#endif

#ifdef incHemf
   real(r4b)    :: Hemf       !< Helium mass fraction
#endif

#ifdef outGammaHI
   real(r4b)    :: gammaHI    !< * time averaged HI photoionization rate
   real(r4b)    :: time       !< * elapsed time in seconds - reset at outputs
#endif

#ifdef incEOS
   real(r4b)    :: eos        !< equation of state variable
#endif

   integer(i8b) :: lasthit    !< * indx of last ray to cross this particle
  
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
  real(r4b), intent(in) :: bot(3)      !< new bottoms
  real(r4b), intent(in) :: top(3)      !< new tops

  where (box%bbound==0) box%bot = bot
  where (box%tbound==0) box%top = top

end subroutine adjustbox


!> scales particles, sources, and the box from comoving to physical values
!==========================================================================
subroutine scale_comoving_to_physical(a,par,src,box,hub)

  character(clen), parameter :: myname="scale_comoving_to_physical"
  integer, parameter :: verb=2
  character(clen) :: str,fmt

  real(r8b), intent(in) :: a                           !< scale factor
  type(particle_type), intent(inout) :: par(:)         !< particles
  type(source_type), optional, intent(inout) :: src(:) !< sources
  type(box_type), optional, intent(inout) :: box       !< box 
  real(r8b), optional, intent(in) :: hub               !< hubble parameter (little h)

  real(r8b) :: h

  if (present(hub)) then 
     h = hub
  else
     h = 1.0d0
  end if

  call mywrite("  scaling comoving to physical coordinates", verb)
  fmt = "(A,F12.5,T22,A,T25,F12.5)"
  write(str,fmt) "  a = ", a, "h = ", h
  call mywrite(str,verb)

  
  par%pos(1) = par%pos(1) * a / h
  par%pos(2) = par%pos(2) * a / h
  par%pos(3) = par%pos(3) * a / h

#ifdef incVel
  par%vel(1) = par%vel(1) * sqrt(a) 
  par%vel(2) = par%vel(2) * sqrt(a) 
  par%vel(3) = par%vel(3) * sqrt(a) 
#endif

  par%mass = par%mass / h
  par%hsml = par%hsml * a / h
  par%rho  = ( par%rho / (a*a*a) ) * (h*h)

  if (present(src)) then  
     src%pos(1) = src%pos(1) * a / h
     src%pos(2) = src%pos(2) * a / h
     src%pos(3) = src%pos(3) * a / h

#ifdef incVel
     src%vel(1) = src%vel(1) * sqrt(a) 
     src%vel(2) = src%vel(2) * sqrt(a)
     src%vel(3) = src%vel(3) * sqrt(a)
#endif

  end if

  if (present(box)) then
     box%top = box%top * a / h
     box%bot = box%bot * a / h
  end if

end subroutine scale_comoving_to_physical


!> scales particles, sources, and the box from physical to comoving  values
!==========================================================================
subroutine scale_physical_to_comoving(a,par,src,box,hub)

  character(clen), parameter :: myname="scale_physical_to_comoving"
  integer, parameter :: verb=2
  character(clen) :: str,fmt

  real(r8b), intent(in) :: a                            !< scale factor 
  type(particle_type), intent(inout) :: par(:)          !< particles
  type(source_type), optional, intent(inout) :: src(:)  !< sources
  type(box_type), optional, intent(inout) :: box        !< box 
  real(r8b), optional, intent(in) :: hub                !< hubble parameter (little h)

  real(r8b) :: h

  if (present(hub)) then 
     h = hub
  else
     h = 1.0d0
  end if

  call mywrite("  scaling physical to comoving coordinates", verb)
  fmt = "(A,F12.5,T22,A,T25,F12.5)"
  write(str,fmt) "  a = ", a, "h = ", h
  call mywrite(str,verb)

  
  par%pos(1) = par%pos(1) / a * h
  par%pos(2) = par%pos(2) / a * h
  par%pos(3) = par%pos(3) / a * h

#ifdef incVel
  par%vel(1) = par%vel(1) / sqrt(a)
  par%vel(2) = par%vel(2) / sqrt(a)
  par%vel(3) = par%vel(3) / sqrt(a)
#endif

  par%mass = par%mass * h
  par%hsml = par%hsml / a * h
  par%rho = par%rho * (a*a*a) / (h*h)

  if (present(src)) then
     src%pos(1) = src%pos(1) / a * h
     src%pos(2) = src%pos(2) / a * h
     src%pos(3) = src%pos(3) / a * h

#ifdef incVel
     src%vel(1) = src%vel(1) / sqrt(a)
     src%vel(2) = src%vel(2) / sqrt(a)
     src%vel(3) = src%vel(3) / sqrt(a)
#endif

  end if

  if (present(box)) then
     box%top = box%top / a * h
     box%bot = box%bot / a * h
  end if

end subroutine scale_physical_to_comoving



!> set electron fraction from ionization fractions
!=======================================================================================
subroutine set_ye(psys, dfltH_mf, dfltHe_mf, ne_bckgnd)

  type(particle_system_type) :: psys
  real(r8b), intent(in) :: dfltH_mf
  real(r8b), intent(in) :: dfltHe_mf
  real(r8b), intent(in) :: ne_bckgnd
  integer(i8b) :: i
  real(r8b) :: Hmf
  real(r8b) :: Hemf
  real(r8b) :: nHe_over_nH


  psys%par(:)%ye = psys%par(:)%xHII + ne_bckgnd


#ifdef incHe

  do i = 1,size(psys%par)
     
!--------------------------
#ifdef incHmf
     Hmf = psys%par(i)%Hmf
#else
     Hmf = dfltH_mf
#endif
!--------------------------
#ifdef incHemf
     Hemf = psys%par(i)%Hemf
#else
     Hemf = dfltHe_mf
#endif
!--------------------------     

     nHe_over_nH = 0.25d0 * Hemf / Hmf
     psys%par(i)%ye = psys%par(i)%ye + &
          ( psys%par(i)%xHeII + 2.0d0 * psys%par(i)%xHeIII ) * nHe_over_nH
     
  end do

#endif


end subroutine set_ye


subroutine set_ye_pars(pars, dfltH_mf, dfltHe_mf, ne_bckgnd)

  type(particle_type) :: pars(:)
  real(r8b), intent(in) :: dfltH_mf
  real(r8b), intent(in) :: dfltHe_mf
  real(r8b), intent(in) :: ne_bckgnd
  integer(i8b) :: i
  real(r8b) :: Hmf
  real(r8b) :: Hemf
  real(r8b) :: nHe_over_nH


  pars(:)%ye = pars(:)%xHII + ne_bckgnd


#ifdef incHe

  do i = 1,size(pars)
     
!--------------------------
#ifdef incHmf
     Hmf = pars(i)%Hmf
#else
     Hmf = dfltH_mf
#endif
!--------------------------
#ifdef incHemf
     Hemf = pars(i)%Hemf
#else
     Hemf = dfltHe_mf
#endif
!--------------------------     

     nHe_over_nH = 0.25d0 * Hemf / Hmf
     pars(i)%ye = pars(i)%ye + &
          ( pars(i)%xHeII + 2.0d0 * pars(i)%xHeIII ) * nHe_over_nH
     
  end do

#endif

end subroutine set_ye_pars




!> sets ionization fractions to their collisional equilibrium values
!=======================================================================================
subroutine set_collisional_ionization_equilibrium(psys, caseA, IsoTemp, DoHydrogen, fit)

  type(particle_system_type) :: psys
  logical, intent(in) :: caseA(2)   
  real(r8b), intent(in) :: IsoTemp
  logical, intent(in) :: DoHydrogen
  character(*), intent(in) :: fit
  real(r8b) :: xvec(5)
  real(r8b) :: Tdum
  integer(i8b) :: i


  ! if we have a single temperature
  !------------------------------------
  if (IsoTemp /= 0.0) then
     call calc_colion_eq_fits(fit, IsoTemp, caseA, xvec)
     if (DoHydrogen) then
        psys%par(:)%xHI    = xvec(1)
        psys%par(:)%xHII   = xvec(2)
     endif
#ifdef incHe
     psys%par(:)%xHeI   = xvec(3)
     psys%par(:)%xHeII  = xvec(4)
     psys%par(:)%xHeIII = xvec(5)
#endif

  ! if we have individual temperatures
  !------------------------------------
  else
     do i = 1,size(psys%par)
        Tdum = psys%par(i)%T
        call calc_colion_eq_fits(fit, Tdum, caseA, xvec)
        if (DoHydrogen) then
           psys%par(i)%xHI    = xvec(1)
           psys%par(i)%xHII   = xvec(2)
        endif
#ifdef incHe
        psys%par(i)%xHeI   = xvec(3)
        psys%par(i)%xHeII  = xvec(4)
        psys%par(i)%xHeIII = xvec(5)
#endif
     end do

  end if

end subroutine set_collisional_ionization_equilibrium


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
subroutine calc_bytes_per_particle_and_source(bpp, bps)

  character(clen), parameter :: myname="calc_bytes_per_particle_and_source"
  logical, parameter :: crash=.true.
  integer, parameter :: verb=2

  integer(i4b), intent(out) :: bpp !< bytes per particle
  integer(i4b), intent(out) :: bps !< bytes per source
  character(clen) :: str

  bpp = 12       ! positions

#ifdef incVel
  bpp = bpp + 12 ! velocities  
#endif

  bpp = bpp + 4  ! ID
  bpp = bpp + 4  ! mass
  bpp = bpp + 4  ! temperature
  bpp = bpp + 4  ! rho
  bpp = bpp + 4  ! ye
  bpp = bpp + 8  ! H ionization fractions
  bpp = bpp + 4  ! hsml

#ifdef cloudy
  bpp = bpp + 4  ! cloudy table xHI
#endif

#ifdef incHmf
  bpp = bpp + 4  ! H mass fraction
#endif

#ifdef incHe
  bpp = bpp + 12 ! He ionization fractions
#endif

#ifdef incHemf
  bpp = bpp + 4  ! He mass fraction
#endif

#ifdef outGammaHI
  bpp = bpp + 8  ! GammaHI tracking
  bpp = bpp + 8  ! time var
#endif

#ifdef incEOS
  bpp = bpp + 4  ! Equation of State variable
#endif

  bpp = bpp + 8  ! last hit index


  write(str,'(A,I4)') "  bytes per particle = ", bpp
  call mywrite(str, verb)

  
  bps = 10 * 4 + 8
  write(str,'(A,I4)') "  bytes per source = ", bps
  call mywrite(str, verb)


end subroutine calc_bytes_per_particle_and_source


!> outputs currently loaded particle data to the screen
!=================================================================
subroutine particle_info_to_screen(psys,str,lun)

  type(particle_system_type), intent(in) :: psys     !< particle system
  character(*), optional, intent(in) :: str          !< arbitrary string
  integer(i4b), optional, intent(in) :: lun          !< if present goes to file
  integer(i4b) :: outlun
  
  
  outlun=stdout
  if (present(lun)) outlun=lun


  99  format(72("-"))  
  100 format(A,T10,3ES15.5)
  101 format(A,T10,2I15,ES15.5)
  102 format(A,T10,2I15)
  103 format(A,T10,3I15)


  write(outlun,99) 
  if (present(str)) write(outlun,"(A)") trim(str)
  write(outlun,"(A,I15,A)") "particle data for ", size(psys%par), "  particles"
!  write(outlun,99) 
  write(outlun,*) 


  write(outlun,100) "xpos", minval(psys%par%pos(1)), maxval(psys%par%pos(1)), &
       meanval_real(psys%par%pos(1))

  write(outlun,100) "ypos", minval(psys%par%pos(2)), maxval(psys%par%pos(2)), &
       meanval_real(psys%par%pos(2))
  
  write(outlun,100) "zpos", minval(psys%par%pos(3)), maxval(psys%par%pos(3)), &
       meanval_real(psys%par%pos(3))   
  
#ifdef incVel
  write(outlun,100) "xvel", minval(psys%par%vel(1)), maxval(psys%par%vel(1)), &
       meanval_real(psys%par%vel(1))
  
  write(outlun,100) "yvel", minval(psys%par%vel(2)), maxval(psys%par%vel(2)), &
       meanval_real(psys%par%vel(2))
  
  write(outlun,100) "zvel", minval(psys%par%vel(3)), maxval(psys%par%vel(3)), &
       meanval_real(psys%par%vel(3))
#endif
  
  write(outlun,102) "id",   minval(psys%par%id), maxval(psys%par%id)
  
  write(outlun,100) "mass", minval(psys%par%mass), maxval(psys%par%mass), &
       meanval_real(psys%par%mass)
  
  write(outlun,100) "T",    minval(psys%par(:)%T), maxval(psys%par(:)%T), &
       meanval_real(psys%par%T)
  
  write(outlun,100) "rho",  minval(psys%par%rho), maxval(psys%par%rho), &
       meanval_real(psys%par%rho)
  
  write(outlun,100) "ye",    minval(psys%par%ye), maxval(psys%par%ye), &
       meanval_real(psys%par%ye)
  
  write(outlun,100) "xHI",    minval(psys%par%xHI), maxval(psys%par%xHI), &
       meanval_real(psys%par%xHI)

#ifdef cloudy
  write(outlun,100) "xHI_cld",    minval(psys%par%xHI_cloudy), maxval(psys%par%xHI_cloudy), &
       meanval_real(psys%par%xHI_cloudy)
#endif
  
  write(outlun,100) "xHII",   minval(psys%par%xHII), maxval(psys%par%xHII), &
       meanval_real(psys%par%xHII)

#ifdef incHmf  
  write(outlun,100) "Hmf",    minval(psys%par%Hmf), maxval(psys%par%Hmf), &
       meanval_real(psys%par%Hmf)
#endif

  write(outlun,100) "hsml", minval(psys%par%hsml), maxval(psys%par%hsml), &
       meanval_real(psys%par%hsml)
  
#ifdef incHe
  write(outlun,100) "xHeI",   minval(psys%par%xHeI), maxval(psys%par%xHeI), &
       meanval_real(psys%par%xHeI)
  
  write(outlun,100) "xHeII",  minval(psys%par%xHeII), maxval(psys%par%xHeII), &
       meanval_real(psys%par%xHeII)
  
  write(outlun,100) "xHeIII", minval(psys%par%xHeIII), maxval(psys%par%xHeIII), &
       meanval_real(psys%par%xHeIII)
#endif


#ifdef incHemf  
  write(outlun,100) "Hemf",    minval(psys%par%Hemf), maxval(psys%par%Hemf), &
       meanval_real(psys%par%Hemf)
#endif


#ifdef outGammaHI
  write(outlun,100) "gammaHI",   minval(psys%par%gammaHI), maxval(psys%par%gammaHI), &
       meanval_real(psys%par%gammaHI)
  
  write(outlun,100) "time(s)",   minval(psys%par%time), maxval(psys%par%time), &
       meanval_real(psys%par%time)
#endif

#ifdef incEOS
  write(outlun,100) "eos",   minval(psys%par%eos), maxval(psys%par%eos), &
       meanval_real(psys%par%eos)  
#endif


  write(outlun,101) "lasthit", minval(psys%par%lasthit), maxval(psys%par%lasthit)

  write(outlun,*)
  
  write(outlun,100) "Box Uppers = ", psys%box%top
  write(outlun,100) "Box Lowers = ", psys%box%bot
  write(outlun,103) "Upr BCs    = ", psys%box%tbound
  write(outlun,103) "Lwr BCs    = ", psys%box%bbound
  
  write(outlun,99) 
  
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
