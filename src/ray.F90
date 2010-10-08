!> \file ray.F90

!> \brief the ray module 
!!
!<

module ray_mod
use myf03_mod
use particle_system_mod
use oct_tree_mod
use mt19937_mod, only: genrand_real1
use spectra_mod, only: rn2freq
use physical_constants_mod, only: HI_th_erg, M_H, M_He
implicit none

  private
  public :: raystatbuffsize
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
  public :: make_skewer_ray


  integer(i8b) :: raystatbuffsize 
  real(r8b), parameter :: one = 1.0d0
  real(r8b), parameter :: zero = 0.0d0
  

!> a ray to be traced through the density field while depositing photons
!-----------------------------------------------------------------------
  type ray_type
     integer(i8b) :: class  !< determines octant
     real(r8b) :: start(3)  !< starting position
     real(r8b) :: dir(3)    !< direction
     real(r8b) :: length    !< length (determines when to stop tracing)
     real(r8b) :: freq      !< freq in HI ionizing units
     real(r8b) :: enrg      !< enrg of a single photon in ergs
     real(r8b) :: pcnt      !< photon count (changes as the ray is depleted)
     real(r8b) :: pini      !< initial photons
     integer(i8b) :: weight !< weight relative to deepest level = 1
     real(r8b) :: dt_s      !< time step associated with ray [s]
  end type ray_type


!> ray stats that can be output for each ray
!-----------------------------------------------------------------------
  type raystat_type
     integer(i4b) :: srcn
     real(r4b) :: start(3)
     real(r4b) :: ryd
  end type raystat_type




contains






!> creates a source ray (as opposed to recombination)
!-----------------------------------------------------------------------  
  subroutine make_source_ray(src, rayn, dtray_s, LumFac, box, ray, length)

    type(source_type), intent(inout) :: src   !< the source
    integer(i8b), intent(in) :: rayn          !< the ray indx
    real(r8b), intent(in) :: dtray_s          !< the time between rays [s]   
    real(r8b), intent(in) :: LumFac           !< pt. src src%lum -> photons/s
    type(box_type), intent(in) :: box         !< simulation box
    type(ray_type), intent(out) :: ray        !< output ray
    real(r8b), intent(in), optional :: length !< optional length (default=huge)

    real(r8b) :: start(3)
    real(r8b) :: xx,yy,zz,r
    real(r8b) :: rn1, rn2
    real(r8b) :: prate
    real(r8b) :: weight
    integer :: i

    integer(i8b) :: sobol_dim
  
!  set the direction of the ray from the emmission profile (src%EmisPrf)
!     0  = isotropic
!    -1  = towards +z
!    -2  = towards -z
!
!  note that for all point sources the luminosity is interpreted as a Flux 
!  [photons/s].  



    select case (src%EmisPrf)

       
    ! this makes rays go in -z direction  
    !-----------------------------------------------------------------------  
    case(-2)

       ray%dir(1) = 0.0
       ray%dir(2) = 0.0
       ray%dir(3) = -1.0       
       weight = 1.0d0
       
    ! this makes rays go in +z direction 
    !-----------------------------------------------------------------------  
    case(-1)

       ray%dir(1) = 0.0
       ray%dir(2) = 0.0
       ray%dir(3) = 1.0
       weight = 1.0d0
                 
    ! random direction on the unit sphere
    !-----------------------------------------------------------------------  
    case(0)

       weight = 1.0d0
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


    ! set the ray starting position
    start = src%pos
    ray%start = start   
    ray%weight = weight

    if ( present(length) ) then
       ray%length = length
    else
       ray%length = huge(1.0d0)
    endif

  
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
       ray%dt_s = dtray_s * (rayn - src%lastemit) * weight
       ray%pini = prate * ray%dt_s
    else
       write(*,*) "make_source_ray> rayn .LT. src%lastemit in ray.f90"
       stop
    end if
    ray%pcnt = ray%pini
    src%lastemit = rayn


  end subroutine make_source_ray

!> creates a recombination ray (as opposed to source)
!-----------------------------------------------------------------------  
  subroutine make_recomb_ray(par, Dflt_Mass, Munit, Dflt_Hmf, Dflt_Hemf, ray, length)

    type(particle_type), intent(in) :: par    !< recombining particle
    real(r8b), intent(in) :: Dflt_Mass        !< default mass
    real(r8b), intent(in) :: Munit            !< cgs mass unit
    real(r8b), intent(in) :: Dflt_Hmf         !< Hydrogen mass fraction
    real(r8b), intent(in) :: Dflt_Hemf        !< Helium mass fraction
    type(ray_type), intent(out) :: ray        !< output ray
    real(r8b), intent(in), optional :: length !< ray length (default=huge)

    real(r8b) :: mass, H_mf, He_mf, H_nuclei, He_nuclei
    real(r8b) :: r,xx,yy,zz
    integer :: i

    mass = Dflt_Mass
#ifdef incmass
    mass = par%mass
#endif

    H_mf = Dflt_Hmf
#ifdef incHmf
    H_mf = par%Hmf
#endif
    
    He_mf = Dflt_Hemf
#ifdef incHemf
    He_mf = par%Hemf
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
    
    if (present(length)) then
       ray%length = length
    else
       ray%length = huge(1.0d0)
    endif
  
!   set the class of the ray (what octant is it going into)
    ray%class=0
    do i=1,3
       if(ray%dir(i).GE.0) ray%class=ray%class+2**(i-1)
    enddo

!   set the frequency and energy / photon of the ray
    ray%freq = 1.0
    ray%enrg = ray%freq * HI_th_erg  

!   set the number of photons in the ray.  
#ifdef incHrec
    ray%pini = H_nuclei * par%xHIIrc
#endif
    ray%pcnt = ray%pini

  end subroutine make_recomb_ray


!> creates a probe ray w/o photons (useful for probing the state of particles)
!------------------------------------------------------------------------------
  subroutine make_probe_ray(pos,dir,ray,length)

    real(r8b), intent(in) :: pos(3)           !< starting position
    real(r8b), intent(in) :: dir(3)           !< direction
    type(ray_type), intent(out) :: ray        !< output ray
    real(r8b), intent(in), optional :: length !< ray length (default=huge)

    real(r8b) :: unit(3)
    integer :: i
    
    unit = dir / sqrt( dot_product(dir,dir) ) 

    ray%start = pos
    ray%dir   = unit

    if (present(length)) then
       ray%length = length
    else
       ray%length = huge(1.0d0)
    endif


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


!> creates a skewer ray
!-----------------------------------------------------------------------  
  subroutine make_skewer_ray(box,ray,length)

    type(box_type), intent(in) :: box       !< simulation box
    type(ray_type), intent(out) :: ray      !< output ray
    real(r8b), intent(in), optional :: length !< ray length (default=huge)  

    real(r8b) :: xx,yy,zz,r
    real(r8b) :: rn, rn1, rn2
    integer :: i
    integer, save :: curface = 1


    curface = curface + 1
    if (curface > 6) curface = 1

    rn1 = genrand_real1()
    rn2 = genrand_real1()

    ray%dir = 0.0

    select case (curface)
       
    case(1)
       ray%dir(1) = 1.0            
       ray%start(1) = box%bots(1)
       ray%start(2) = box%bots(2) + rn1 * (box%tops(2)-box%bots(2))
       ray%start(3) = box%bots(3) + rn2 * (box%tops(3)-box%bots(3))
    case(2)
       ray%dir(2) = 1.0             
       ray%start(2) = box%bots(2)
       ray%start(3) = box%bots(3) + rn1 * (box%tops(3)-box%bots(3))
       ray%start(1) = box%bots(1) + rn2 * (box%tops(1)-box%bots(1))
    case(3)
       ray%dir(3) = 1.0             
       ray%start(3) = box%bots(3)
       ray%start(1) = box%bots(1) + rn1 * (box%tops(1)-box%bots(1))
       ray%start(2) = box%bots(2) + rn2 * (box%tops(2)-box%bots(2))     
    case(4)
       ray%dir(1) = -1.0             
       ray%start(1) = box%tops(1)
       ray%start(2) = box%bots(2) + rn1 * (box%tops(2)-box%bots(2))
       ray%start(3) = box%bots(3) + rn2 * (box%tops(3)-box%bots(3))
    case(5)
       ray%dir(2) = -1.0
       ray%start(2) = box%tops(2)
       ray%start(3) = box%bots(3) + rn1 * (box%tops(3)-box%bots(3))
       ray%start(1) = box%bots(1) + rn2 * (box%tops(1)-box%bots(1))
    case(6)
       ray%dir(3) = -1.0
       ray%start(3) = box%tops(3)
       ray%start(1) = box%bots(1) + rn1 * (box%tops(1)-box%bots(1))
       ray%start(2) = box%bots(2) + rn2 * (box%tops(2)-box%bots(2))
    case default 
       stop "curface out of bounds"
    end select


    if (present(length)) then
       ray%length = length
    else
       ray%length = huge(1.0d0)
    endif


  
!   set the class of the ray (what octant is it going into)
    ray%class=0
    do i=1,3
       if(ray%dir(i).GE.0) ray%class=ray%class+2**(i-1)
    enddo

!   set rest to zero 
    ray%freq = 0.0
    ray%enrg = 0.0
    ray%pini = 0.0
    ray%pcnt = 0.0


  end subroutine make_skewer_ray





!> properly sets the starting point, direction, class, and length of a ray
!--------------------------------------------------------------------------
  subroutine set_ray(ray,start,dir,len)
    type(ray_type) :: ray     !< inout ray
    real(r8b) :: start(3)     !< starting position
    real(r8b) :: dir(3)       !< direction 
    real(r8b),optional :: len !< length
    integer :: i
    
    real(r8b) :: unit(3)

    unit = dir / sqrt( dot_product(dir,dir) )

    ray%start = start
    ray%dir   = unit
    ray%class = 0

    do i=1,3
       if(dir(i).GE.0) ray%class=ray%class+2**(i-1)
    enddo

    if (present(len)) then
       ray%length = len
    else
       ray%length = huge(1.0d0)
    endif


  end subroutine set_ray
  
!> returns a transformed ray while preserving the initial ray
!--------------------------------------------------------------
  subroutine transform_ray(ray, inray, pm)
    type(ray_type) :: ray           !< output transformed ray
    type(ray_type) :: inray         !< input ray
    type(transformation_type) :: pm !< transformation
    integer :: i
    ray%start  = inray%start * pm%fac + pm%shift
    ray%dir    = inray%dir   * pm%fac
    ray%length = inray%length  
    ray%class  = 0
    do i=1,3
       if (ray%dir(i) >= 0) ray%class=ray%class+2**(i-1)
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

    d     = sum( (part%pos - ray%start) * ray%dir  )    ! sp dot dir
    dist2 = sum( (part%pos-d*ray%dir-ray%start)**2 )    ! dist^2

    if (present(dd)) dd = d

  end function pdist2ray


!> returns the distance^2 between a point and a ray
!--------------------------------------------------------------  
  function dist2ray(ray, pos, dd) result(dist2)
    type(ray_type) :: ray       !< input ray
    real(r4b) :: pos(3)         !< input position
    real(r8b), optional :: dd   !< GT 0 when point is in front of ray start
    real(r8b) :: dist2          !< distance^2

    real(r8b) :: dotp
    real(r8b) :: diff(3)

    real(r8b) :: vec1(3)
    real(r8b) :: vec2(3)

    vec1 = pos - ray%start
    vec2 = ray%dir

    dotp  = dot_product( vec1, vec2 )
    diff  = vec1 - dotp * ray%dir 
    dist2 = dot_product( diff, diff )

    if( present(dd) ) dd = dotp

  end function dist2ray
  
!> tests for ray / particle intersection. 
!----------------------------------------------  
  function part_intersection(ray, part) result(hit)
    logical :: hit               !< true or false result
    type(ray_type) :: ray        !< ray
    type(particle_type) :: part  !< particle
    real(r8b) :: dist2           !< distance to particle squared
    real(r8b) :: dotp            !< projected distance along ray

    dist2 = dist2ray(ray, part%pos, dotp)
    hit = dist2 < part%hsml * part%hsml .and. dotp >= zero .and. dotp <= ray%length

  end function part_intersection
  
!> tests for ray / cell intersection. 
!----------------------------------------------  
  function cell_intersection(ray,cell) result(hit)
    logical :: hit           !< true or false result
    type(ray_type) :: ray    !< ray
    type(cell_type) :: cell  !< cell
    real(r8b) bot(3),top(3)
    bot = cell%botrange - ray%start
    top = cell%toprange - ray%start
    hit = pluecker(ray, bot, top)
  end function cell_intersection
  

!> pluecker test for line segment / cell intersection
!-----------------------------------------------------    
  function pluecker(ray,bot,top) result(hit)
    logical :: hit            !< true or false result
    type(ray_type) :: ray     !< input ray
    real(r8b) :: bot(3)       !< lower cell corner
    real(r8b) :: top(3)       !< upper cell corner

    real(r8b) :: dir(3)

    dir=ray%dir  
    hit=.FALSE.
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
    hit=.TRUE.
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
