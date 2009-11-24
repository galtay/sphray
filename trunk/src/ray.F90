!> \file ray.F90

!> \brief the ray module 
!!
!<

module ray_mod
use myf90_mod
use particle_system_mod
use oct_tree_mod
use mt19937_mod, only: genrand_real1
implicit none

integer, parameter :: raystatbuffsize = 5000

!> a ray to be traced through the density field while depositing photons
!-----------------------------------------------------------------------
  type ray_type
     integer(i8b) :: class !< determines octant
     real(r8b) :: start(3) !< starting position
     real(r8b) :: dir(3)   !< direction
     real(r8b) :: length   !< length (determines when to stop tracing)
     real(r8b) :: freq     !< freq in HI ionizing units
     real(r8b) :: enrg     !< enrg of a single photon in ergs
     real(r8b) :: pcnt     !< photon count (changes as the ray is depleted)
     real(r8b) :: pini     !< initial photons
  end type ray_type

  type raystat_type
     integer(i4b) :: srcn
     real(r4b) :: dir(3)
     real(r4b) :: freq
     real(r4b) :: temit
     real(r8b) :: pini
     real(r8b) :: pcnt
     integer(i8b) :: pindx
     real(r4b) :: b
     real(r4b) :: d
  end type raystat_type

  private
  public :: raystatbuffsize
  public :: curface
  public :: ray_type
  public :: raystat_type
  public :: set_ray
  public :: cell_intersection
  public :: part_intersection
  public :: transform_ray
  public :: dist2ray
  public :: make_source_ray
  public :: make_probe_ray
  public :: make_recomb_ray


  integer :: curface  !< used to track rays emitted from box faces

contains

!> creates a source ray (as opposed to recombination)
!-----------------------------------------------------------------------  
  subroutine make_source_ray(src,rayn,dtray,LumFac,box,ray)
  use particle_system_mod, only: source_type, box_type
  use spectra_mod, only: rn2freq
  use physical_constants_mod, only: HI_th_erg

    type(source_type), intent(inout) :: src !< the source
    integer(i8b), intent(in) :: rayn        !< the ray indx
    real(r8b), intent(in) :: dtray          !< the time between rays
    real(r8b), intent(in) :: LumFac         !< pt. src src%lum -> photons/s
    type(box_type), intent(in) :: box       !< simulation box
    type(ray_type), intent(out) :: ray      !< output ray
  
    real(r8b) :: xx,yy,zz,r
    real(r8b) :: rn, rn1, rn2
    real(r8b) :: prate
    integer :: i


  
!  set the direction of the ray from the emmission profile (src%EmisPrf)
!     0  = isotropic
!    -1  = z=0 plane towards +z
!    -2  = z=blxlen plane towards -z
!    -3  = all planes into box
!
!  note that for all sources the luminosity is interpreted as a Flux 
!  [photons/s].  For point sources (src%EmisPrf >= 0) this flux is just 
!  src%L * LumFac.  For extended sources, (src%EmisPrf < 0) the luminosity 
!  that is read in is expected to be a photon number density, n 
!  [photons/cm^3].  In this case this input value is converted in 
!  read_src_snapshot to a Flux using the area of the emitting planes such that 
!  the flux F (photons/s) emitted from the planes would produce the number 
!  density, n,  in an optically thin volume.


    select case (src%EmisPrf)


       ! this cycles through all planes and casts into the box
       ! the position coordinate is used to track which face is next
       case(-3)

          curface = curface + 1
          if (curface > 6) curface = 1

          rn1 = genrand_real1()
          rn2 = genrand_real1()

          ray%dir = 0.0

          select case (curface)
             
             case(1)
                ray%dir(1) = 1.0            
                ray%start(1) = box%bot(1)
                ray%start(2) = box%bot(2) + rn1 * (box%top(2)-box%bot(2))
                ray%start(3) = box%bot(3) + rn2 * (box%top(3)-box%bot(3))
             case(2)
                ray%dir(2) = 1.0             
                ray%start(2) = box%bot(2)
                ray%start(3) = box%bot(3) + rn1 * (box%top(3)-box%bot(3))
                ray%start(1) = box%bot(1) + rn2 * (box%top(1)-box%bot(1))
             case(3)
                ray%dir(3) = 1.0             
                ray%start(3) = box%bot(3)
                ray%start(1) = box%bot(1) + rn1 * (box%top(1)-box%bot(1))
                ray%start(2) = box%bot(2) + rn2 * (box%top(2)-box%bot(2))     
             case(4)
                ray%dir(1) = -1.0             
                ray%start(1) = box%top(1)
                ray%start(2) = box%bot(2) + rn1 * (box%top(2)-box%bot(2))
                ray%start(3) = box%bot(3) + rn2 * (box%top(3)-box%bot(3))
             case(5)
                ray%dir(2) = -1.0
                ray%start(2) = box%top(2)
                ray%start(3) = box%bot(3) + rn1 * (box%top(3)-box%bot(3))
                ray%start(1) = box%bot(1) + rn2 * (box%top(1)-box%bot(1))
             case(6)
                ray%dir(3) = -1.0
                ray%start(3) = box%top(3)
                ray%start(1) = box%bot(1) + rn1 * (box%top(1)-box%bot(1))
                ray%start(2) = box%bot(2) + rn2 * (box%top(2)-box%bot(2))
             case default 
                stop "curface out of bounds"
          end select


       ! this makes the z=boxlen plane a source rays go in -z direction  
       case(-2)

          ray%dir(1) = 0.0
          ray%dir(2) = 0.0
          ray%dir(3) = -1.0
          
          ray%start(3) = box%top(3)
          
          rn = genrand_real1()
          ray%start(2) = box%bot(2) + rn * (box%top(2)-box%bot(2))
          rn = genrand_real1()
          ray%start(1) = box%bot(1) + rn * (box%top(1)-box%bot(1))


       ! this makes the bottom z=0 plane a source rays go in +z direction 
       case(-1)
       
          ray%dir(1) = 0.0
          ray%dir(2) = 0.0
          ray%dir(3) = 1.0
          
          ray%start(3) = box%bot(3)
          
          rn = genrand_real1()
          ray%start(2) = box%bot(2) + rn * (box%top(2)-box%bot(2))
          rn = genrand_real1()
          ray%start(1) = box%bot(1) + rn * (box%top(1)-box%bot(1))


       ! random direction on the unit sphere
       case(0)

          ray%start = src%pos
          r=2. 
          do while ( r .GT. 1.0 .and. r .NE. 0.0 )
             xx=(2*genrand_real1()-1)   
             yy=(2*genrand_real1()-1)   
             zz=(2*genrand_real1()-1)   
             r=xx*xx+yy*yy+zz*zz
          enddo
          r = sqrt(r)
          ray%dir(1) = xx/r  ! it is important that ray%dir be a unit vector
          ray%dir(2) = yy/r
          ray%dir(3) = zz/r          
          
       case default

          write(*,*) "emission profile not recognized"
          write(*,*) "profile = ", src%EmisPrf
          stop          

    end select




    ray%length=sqrt(sum(ray%dir**2))

  
!   set the class of the ray (what octant is it going into)
    ray%class=0
    do i=1,3
       if(ray%dir(i).GE.0) ray%class=ray%class+2**(i-1)
    enddo

!   set the frequency and energy / photon of the ray
    ray%freq = rn2freq(src%SpcType)
    ray%enrg = ray%freq * HI_th_erg  


!   set the number of photons in the ray
    if (ray%enrg > 0.) then
       prate = src%L * LumFac ! photons per sec
    else 
       prate = 0.
    end if
    if (rayn > src%lastemit) then
       ray%pini = prate * dtray * (rayn - src%lastemit)
    else
       write(*,*) "make_source_ray> rayn .LT. src%lastemit in ray.f90"
       stop
    end if
    ray%pcnt = ray%pini
    src%lastemit = rayn


  end subroutine make_source_ray

!> creates a recombination ray (as opposed to source)
!-----------------------------------------------------------------------  
  subroutine make_recomb_ray(par, DfltMass, Munit, H_mf, He_mf, ray)
  use physical_constants_mod, only: HI_th_erg, M_H, M_He

    type(particle_type), intent(in) :: par  !< recombining particle
    real(r8b), intent(in) :: DfltMass            !< default mass
    real(r8b), intent(in) :: Munit               !< cgs mass unit
    real(r8b), intent(in) :: H_mf                !< Hydrogen mass fraction
    real(r8b), intent(in) :: He_mf               !< Helium mass fraction
    type(ray_type), intent(out) :: ray      !< output ray

    real(r8b) :: mass, H_nuclei, He_nuclei
    real(r8b) :: r,xx,yy,zz
    integer :: i

    mass = DfltMass
#ifdef incmass
    mass = par%mass
#endif

    H_nuclei  = mass * Munit * H_mf  / M_H
    He_nuclei = mass * Munit * He_mf / M_He


    !  set the direction of the ray (random direction on unit sphere)  
    ray%start = par%pos
    r=2. 
    do while ( r .GT. 1.0 .and. r .NE. 0.0 )
       xx=(2*genrand_real1()-1)   
       yy=(2*genrand_real1()-1)   
       zz=(2*genrand_real1()-1)   
       r=xx*xx+yy*yy+zz*zz
    enddo
    r = sqrt(r)
    ray%dir(1) = xx/r  ! it is important that ray%dir be a unit vector
    ray%dir(2) = yy/r
    ray%dir(3) = zz/r
    
    ray%length=sqrt(sum(ray%dir**2))

  
!   set the class of the ray (what octant is it going into)
    ray%class=0
    do i=1,3
       if(ray%dir(i).GE.0) ray%class=ray%class+2**(i-1)
    enddo

!   set the frequency and energy / photon of the ray
    ray%freq = 1.0
    ray%enrg = ray%freq * HI_th_erg  

!   set the number of photons in the ray.  need the macro for compiler issues
#ifdef increc
    ray%pini = H_nuclei * par%xHIIrc
#endif
    ray%pcnt = ray%pini

  end subroutine make_recomb_ray


!> creates a probe ray w/o photons (useful for probing the state of particles)
!------------------------------------------------------------------------------
  subroutine make_probe_ray(pos,dir,ray)

    real(r8b), intent(in) :: pos(3)         !< starting position
    real(r8b), intent(in) :: dir(3)         !< direction
    type(ray_type), intent(out) :: ray !< output ray
  
    integer :: i
    
    ray%start = pos
    ray%dir = dir
    ray%length = sqrt(sum(dir*dir))
    ray%dir = ray%dir / ray%length
    ray%length = 1.0


!   set the class of the ray (what quadrant is it going into)
    ray%class=0
    do i=1,3
       if(ray%dir(i).GE.0) ray%class=ray%class+2**(i-1)
    enddo

    ray%freq = 0.0
    ray%enrg = 0.0
    ray%pcnt = 0.0
    ray%pini = 0.0

  end subroutine make_probe_ray

!> properly sets the starting point, direction, class, and length of a ray
!--------------------------------------------------------------------------
  subroutine set_ray(ray,start,dir,len)
    type(ray_type) :: ray   !< inout ray
    real(r8b) :: start(3)        !< starting position
    real(r8b) :: dir(3)          !< direction 
    real(r8b),optional :: len    !< length
    integer :: i
    
    ray%start=start
    ray%length=SQRT(sum(dir**2))
    ray%dir=dir/ray%length  
    if(present(len)) ray%length=len  
    ray%class=0
    do i=1,3
       if(dir(i).GE.0) ray%class=ray%class+2**(i-1)
    enddo
  end subroutine set_ray
  
!> returns a transformed ray while preserving the initial ray
!--------------------------------------------------------------
  subroutine transform_ray(ray,inray,pm)
    type(ray_type) :: ray           !< output transformed ray
    type(ray_type) :: inray         !< input ray
    type(transformation_type) :: pm !< transformation
    integer :: i
    ray%start=inray%start*pm%fac+pm%shift
    ray%dir=inray%dir*pm%fac
    ray%length=inray%length  
    ray%class=0
    do i=1,3
       if(ray%dir(i).GE.0) ray%class=ray%class+2**(i-1)
    enddo
  end subroutine transform_ray


!> returns the distance^2 between a particle and a ray
!--------------------------------------------------------------  
  function pdist2ray(ray,part,dd) result(dist2)
    type(ray_type) :: ray        !< input ray
    type(particle_type) :: part  !< input particle
    real(r8b) :: dist2                !< distance^2
    real(r8b), optional :: dd         !< GT 0 when particle is in front of ray start
    real(r8b) :: d

    d=SUM((part%pos-ray%start)*ray%dir)           ! sp dot dir
    dist2=SUM((part%pos-d*ray%dir-ray%start)**2)  ! dist^2
    if(present(dd)) dd=d
  end function pdist2ray

!> returns the distance^2 between a point and a ray
!--------------------------------------------------------------  
  function dist2ray(ray,pos,dd) result(dist2)
    type(ray_type) :: ray  !< input ray
    real(r4b) :: pos(3)    !< input position
    real(r8b) :: dist2          !< distance^2
    real(r8b), optional :: dd   !< GT 0 when point is in front of ray start
    real(r8b) :: dotp

    dotp=SUM((pos-ray%start)*ray%dir)
    dist2=SUM((pos-dotp*ray%dir-ray%start)**2)
    if(present(dd)) dd=dotp
  end function dist2ray
  
!> tests for ray / particle intersection. 
!----------------------------------------------  
  function part_intersection(ray,part) result(x)
    logical :: x                 !< true or false result
    type(ray_type) :: ray        !< ray
    type(particle_type) :: part  !< particle
    real(r8b) dist2,dotp
    dist2=dist2ray(ray,part%pos,dotp)
    x=dist2.LT.part%hsml**2.and.dotp.GE.0
  end function part_intersection
  
!> tests for ray / cell intersection. 
!----------------------------------------------  
  function cell_intersection(ray,cell) result(x)
    logical :: x             !< true or false result
    type(ray_type) :: ray    !< ray
    type(cell_type) :: cell  !< cell
    real(r8b) bot(3),top(3)
    bot=cell%botrange-ray%start
    top=cell%toprange-ray%start
    x=pluecker(ray,bot,top)
  end function cell_intersection
  

!> pluecker test for line segment / cell intersection
!-----------------------------------------------------    
  function pluecker(ray,bot,top) result(x)
    logical :: x                !< true or false result
    type(ray_type) :: ray       !< input ray
    real(r8b) bot(3)                 !< lower cell corner
    real(r8b) top(3)                 !< upper cell corner

    real(r8b) dir(3)

    dir=ray%dir  
    x=.FALSE.
    select case(ray%class)
    case(0)
       if(bot(1).GT.0.OR.bot(2).GT.0.OR.bot(3).GT.0) return
       if(dir(1)*bot(2)-dir(2)*top(1).LT.0.OR. &
            dir(1)*top(2)-dir(2)*bot(1).GT.0.OR. &
            dir(1)*top(3)-dir(3)*bot(1).GT.0.OR. &
            dir(1)*bot(3)-dir(3)*top(1).LT.0.OR. &
            dir(2)*bot(3)-dir(3)*top(2).LT.0.OR. &
            dir(2)*top(3)-dir(3)*bot(2).GT.0.) return
    case(1)
       if(top(1).LT.0.OR.bot(2).GT.0.OR.bot(3).GT.0) return
       if(dir(1)*top(2)-dir(2)*top(1).LT.0.OR. &
            dir(1)*bot(2)-dir(2)*bot(1).GT.0.OR. &
            dir(1)*bot(3)-dir(3)*bot(1).GT.0.OR. &
            dir(1)*top(3)-dir(3)*top(1).LT.0.OR. &
            dir(2)*bot(3)-dir(3)*top(2).LT.0.OR. &
            dir(2)*top(3)-dir(3)*bot(2).GT.0.) return
    case(2)
       if(bot(1).GT.0.OR.top(2).LT.0.OR.bot(3).GT.0) return
       if(dir(1)*bot(2)-dir(2)*bot(1).LT.0.OR. &
            dir(1)*top(2)-dir(2)*top(1).GT.0.OR. &
            dir(1)*top(3)-dir(3)*bot(1).GT.0.OR. &
            dir(1)*bot(3)-dir(3)*top(1).LT.0.OR. &
            dir(2)*top(3)-dir(3)*top(2).LT.0.OR. &
            dir(2)*bot(3)-dir(3)*bot(2).GT.0.) return
    case(3)
       if(top(1).LT.0.OR.top(2).LT.0.OR.bot(3).GT.0) return
       if(dir(1)*top(2)-dir(2)*bot(1).LT.0.OR. &
            dir(1)*bot(2)-dir(2)*top(1).GT.0.OR. &
            dir(1)*bot(3)-dir(3)*bot(1).GT.0.OR. &
            dir(1)*top(3)-dir(3)*top(1).LT.0.OR. &
            dir(2)*top(3)-dir(3)*top(2).LT.0.OR. &
            dir(2)*bot(3)-dir(3)*bot(2).GT.0.) return
    case(4)
       if(bot(1).GT.0.OR.bot(2).GT.0.OR.top(3).LT.0) return
       if(dir(1)*bot(2)-dir(2)*top(1).LT.0.OR. &
            dir(1)*top(2)-dir(2)*bot(1).GT.0.OR. &
            dir(1)*top(3)-dir(3)*top(1).GT.0.OR. &
            dir(1)*bot(3)-dir(3)*bot(1).LT.0.OR. &
            dir(2)*bot(3)-dir(3)*bot(2).LT.0.OR. &
            dir(2)*top(3)-dir(3)*top(2).GT.0.) return
    case(5)
       if(top(1).LT.0.OR.bot(2).GT.0.OR.top(3).LT.0) return
       if(dir(1)*top(2)-dir(2)*top(1).LT.0.OR. &
            dir(1)*bot(2)-dir(2)*bot(1).GT.0.OR. &
            dir(1)*bot(3)-dir(3)*top(1).GT.0.OR. &
            dir(1)*top(3)-dir(3)*bot(1).LT.0.OR. &
            dir(2)*bot(3)-dir(3)*bot(2).LT.0.OR. &
            dir(2)*top(3)-dir(3)*top(2).GT.0.) return
    case(6)
       if(bot(1).GT.0.OR.top(2).LT.0.OR.top(3).LT.0) return
       if(dir(1)*bot(2)-dir(2)*bot(1).LT.0.OR. &
            dir(1)*top(2)-dir(2)*top(1).GT.0.OR. &
            dir(1)*top(3)-dir(3)*top(1).GT.0.OR. &
            dir(1)*bot(3)-dir(3)*bot(1).LT.0.OR. &
            dir(2)*top(3)-dir(3)*bot(2).LT.0.OR. &
            dir(2)*bot(3)-dir(3)*top(2).GT.0.) return
    case(7)
       if(top(1).LT.0.OR.top(2).LT.0.OR.top(3).LT.0) return
       if(dir(1)*top(2)-dir(2)*bot(1).LT.0.OR. &
            dir(1)*bot(2)-dir(2)*top(1).GT.0.OR. &
            dir(1)*bot(3)-dir(3)*top(1).GT.0.OR. &
            dir(1)*top(3)-dir(3)*bot(1).LT.0.OR. &
            dir(2)*top(3)-dir(3)*bot(2).LT.0.OR. &
            dir(2)*bot(3)-dir(3)*top(2).GT.0.) return
    case default
       call rayError('ray class.')
    end select
    x=.TRUE.
  end function pluecker

!> error handling
!-----------------------------      
  subroutine rayError(string,i)
    character(*) :: string  !< error string
    integer, optional :: i  !< error number
    
    print*,' Error detected:'
    
    if(present(i)) then
       print*,string,i
    else
       print*,string
    endif
    
    stop
  end subroutine rayError
  
end module ray_mod
