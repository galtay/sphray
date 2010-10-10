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
  public :: base_ray_type
  public :: src_ray_type
  public :: raystat_type



  integer(i8b) :: raystatbuffsize 
  real(r8b), parameter :: zero = 0.0d0
  real(r8b), parameter :: one = 1.0d0

  
!> bare bones ray class
!---------------------------------------------------------------------
type base_ray_type
   real(r8b) :: start(3)  !< starting position
   real(r8b) :: dir(3)    !< unit vector direction
   real(r8b) :: length    !< length (determines when to stop tracing)
   integer(i4b) :: class  !< based on direction signs (MMM, PMM, ...)
 contains
   procedure :: class_from_dir => ray_class_from_dir !< computes class from dir
   procedure :: dist2pt        !< perpendicular distance to point
   procedure :: pluecker       !< pluecker AABB intersection test
   procedure :: part_intersection !< intersects a particle? 
   procedure :: cell_intersection !< intersects an AABB?
end type base_ray_type


!> adds source and photon properties to base ray class
!-----------------------------------------------------------------------
  type, extends(base_ray_type) :: src_ray_type
     real(r8b) :: freq      !< freq in HI ionizing units
     real(r8b) :: enrg      !< enrg of a single photon in ergs
     real(r8b) :: pcnt      !< photon count (changes as the ray is depleted)
     real(r8b) :: pini      !< initial photons
     real(r8b) :: dt_s      !< time step associated with ray [s]
   contains
     procedure :: make_src_ray => make_src_ray            !< creates ray
     procedure :: transform => transform_src_ray  !< transform for BCs
  end type src_ray_type


!> ray stats that can be output for each ray
!-----------------------------------------------------------------------
  type raystat_type
     integer(i4b) :: srcn
     real(r4b) :: start(3)
     real(r4b) :: ryd
  end type raystat_type




contains




!> creates a source ray 
!-----------------------------------------------------------------------  
  subroutine make_src_ray(ray, src, rayn, dtray_s, Lunit, box, length)

    class(src_ray_type), intent(out) :: ray   !< ray to make
    type(source_type), intent(inout) :: src   !< source
    integer(i8b), intent(in) :: rayn          !< ray indx
    real(r8b), intent(in) :: dtray_s          !< time between rays [s]   
    real(r8b), intent(in) :: Lunit            !< converts src%lum -> photons/s
    type(box_type), intent(in) :: box         !< simulation box

    real(r8b), intent(in), optional :: length !< optional length (default=huge)

    real(r8b) :: xx,yy,zz,r
    real(r8b) :: rn1, rn2
    real(r8b) :: prate
    integer :: i

  
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
       
    ! this makes rays go in +z direction 
    !-----------------------------------------------------------------------  
    case(-1)

       ray%dir(1) = 0.0
       ray%dir(2) = 0.0
       ray%dir(3) = 1.0
                 
    ! random direction on the unit sphere
    !-----------------------------------------------------------------------  
    case(0)

       r=2.0d0 
       do while ( r .GT. 1.0d0 .and. r .NE. 0.0d0 )
          xx=(2.0d0 * genrand_real1()-1.0d0)   
          yy=(2.0d0 * genrand_real1()-1.0d0)   
          zz=(2.0d0 * genrand_real1()-1.0d0)   
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
    !-----------------------------------------------------------------------  

    ray%start = src%pos

    if ( present(length) ) then
       ray%length = length
    else
       ray%length = huge(1.0d0)
    endif

  
!   set the class of the ray (what octant is it going into)
    call ray%class_from_dir()

!   set the frequency and energy / photon of the ray
    ray%freq = rn2freq(src%SpcType)
    ray%enrg = ray%freq * HI_th_erg  


!   set the number of photons in the ray
    if (ray%enrg > 0.) then
       prate = src%L * Lunit ! photons per sec
    else 
       prate = 0.
    end if
    if (rayn > src%lastemit) then
       ray%dt_s = dtray_s * (rayn - src%lastemit) 
       ray%pini = prate * ray%dt_s
    else
       write(*,*) "make_source_ray> rayn .LT. src%lastemit in ray.f90"
       stop
    end if
    ray%pcnt = ray%pini
    src%lastemit = rayn

  end subroutine make_src_ray


! pre computes the class of the ray for the Pluecker test
! ray label    class
!   MMM          0
!   PMM          1
!   MPM          2
!   PPM          3
!   MMP          4
!   PMP          5
!   MPP          6
!   PPP          7
!-----------------------------------------------------------
subroutine ray_class_from_dir( ray ) 
  class(base_ray_type) :: ray
  integer(i4b) :: i

  ray%class = 0
  do i = 1, 3
     if ( ray%dir(i) >= zero ) ray%class = ray%class + 2**(i-1)
  end do
  
end subroutine ray_class_from_dir




!> returns a transformed ray while preserving the initial ray
!--------------------------------------------------------------
  subroutine transform_src_ray(inray, outray, trans)
    class(src_ray_type) :: inray       !< input ray
    type(src_ray_type) :: outray       !< output ray
    type(transformation_type) :: trans !< transformation
    
    outray%start  = inray%start * trans%fac + trans%shift
    outray%dir    = inray%dir   
    outray%length = inray%length  
    outray%class  = inray%class

    outray%freq = inray%freq
    outray%enrg = inray%enrg
    outray%pcnt = inray%pcnt
    outray%pini = inray%pini
    outray%dt_s = inray%dt_s

  end subroutine transform_src_ray


!> returns the perpendicular distance between a point and a ray
!--------------------------------------------------------------  
function dist2pt(ray, pos, proj) result(dist)
  class(base_ray_type), intent(in) :: ray  !< calling ray
  real(r4b), intent(in) :: pos(3)    !< point
  real(r8b), optional :: proj        !< distance along ray 
  real(r8b) :: dist                  !< distance perp. to ray
  
  real(r8b) :: dotp
  real(r8b) :: diff(3)
  real(r8b) :: dist2
  real(r8b) :: vec1(3)
  real(r8b) :: vec2(3)
  
  vec1 = pos - ray%start
  vec2 = ray%dir
  
  dotp  = dot_product( vec1, vec2 )
  diff  = vec1 - dotp * ray%dir 
  dist2 = dot_product( diff, diff )
  dist = sqrt(dist2)
  
  if( present(proj) ) proj = dotp
  
end function dist2pt



!> tests for ray / particle intersection.  
!----------------------------------------------  
  function part_intersection(ray, part) result(hit)
    class(base_ray_type) :: ray    !< ray
    type(particle_type) :: part    !< particle
    logical :: hit                 !< true or false result

    real(r8b) :: start2cen       !< ray start to particle position
    real(r8b) :: dist            !< perpendicular distance to particle 
    real(r8b) :: proj            !< projected distance along ray

    start2cen = sqrt( sum( (part%pos - ray%start)*(part%pos - ray%start) ) )
    if (start2cen < part%hsml) then
       hit = .true.
       return
    endif

    dist = ray%dist2pt(part%pos, proj)
    if (dist >= part%hsml) then
       hit = .false.
       return
    endif

    hit = dist < part%hsml .and. proj >= zero 

  end function part_intersection


!> tests for ray - AABB intersection. 
!----------------------------------------------  
  function cell_intersection(ray,cell) result(hit)
    class(base_ray_type) :: ray    !< ray
    type(cell_type) :: cell        !< cell
    logical :: hit                 !< true or false result
    real(r8b) :: bot(3)
    real(r8b) :: top(3)
    bot = cell%botrange - ray%start
    top = cell%toprange - ray%start
    hit = ray%pluecker(bot, top)
  end function cell_intersection


!> pluecker test for line segment / cell intersection
!-----------------------------------------------------    
function pluecker(ray, s2b, s2t) result( hit )

  logical :: hit              !< true or false result
  class(base_ray_type) :: ray !< input ray
  real(r8b) :: s2b(3)         !< vector from ray start to lower cell corner
  real(r8b) :: s2t(3)         !< vector from ray start to upper cell corner
  
  real(r8b) :: dir(3)
  real(r8b) :: dist

  real(r8b) :: e2b(3)       !< vector from ray end to lower cell corner
  real(r8b) :: e2t(3)       !< vector from ray end to upper cell corner

  dir  = ray%dir  
  dist = ray%length

  e2b = s2b - dir * dist
  e2t = s2t - dir * dist

  hit = .false.

  ! branch on ray direction
  !---------------------------
  select case( ray%class )

     ! MMM
     !-----------
  case(0)

     if(s2b(1) > zero .or. s2b(2) > zero .or. s2b(3) > zero) return ! on negative part of ray 
     if(e2t(1) < zero .or. e2t(2) < zero .or. e2t(3) < zero) return ! past length of ray      

     if ( dir(1)*s2b(2) - dir(2)*s2t(1) < zero .or.  &
          dir(1)*s2t(2) - dir(2)*s2b(1) > zero .or.  &
          dir(1)*s2t(3) - dir(3)*s2b(1) > zero .or.  &
          dir(1)*s2b(3) - dir(3)*s2t(1) < zero .or.  &
          dir(2)*s2b(3) - dir(3)*s2t(2) < zero .or.  &
          dir(2)*s2t(3) - dir(3)*s2b(2) > zero       ) return
     
     ! PMM
     !-----------
  case(1)
     
     if(s2t(1) < zero .or. s2b(2) > zero .or. s2b(3) > zero) return ! on negative part of ray 
     if(e2b(1) > zero .or. e2t(2) < zero .or. e2t(3) < zero) return ! past length of ray      
     
     if ( dir(1)*s2t(2) - dir(2)*s2t(1) < zero .or.  &
          dir(1)*s2b(2) - dir(2)*s2b(1) > zero .or.  &
          dir(1)*s2b(3) - dir(3)*s2b(1) > zero .or.  &
          dir(1)*s2t(3) - dir(3)*s2t(1) < zero .or.  &
          dir(2)*s2b(3) - dir(3)*s2t(2) < zero .or.  &
          dir(2)*s2t(3) - dir(3)*s2b(2) > zero       ) return
     
     ! MPM
     !-----------
  case(2)
     
     if(s2b(1) > zero .or. s2t(2) < zero .or. s2b(3) > zero) return ! on negative part of ray 
     if(e2t(1) < zero .or. e2b(2) > zero .or. e2t(3) < zero) return ! past length of ray      
     
     if ( dir(1)*s2b(2) - dir(2)*s2b(1) < zero .or.  &
          dir(1)*s2t(2) - dir(2)*s2t(1) > zero .or.  &
          dir(1)*s2t(3) - dir(3)*s2b(1) > zero .or.  &
          dir(1)*s2b(3) - dir(3)*s2t(1) < zero .or.  &
          dir(2)*s2t(3) - dir(3)*s2t(2) < zero .or.  &
          dir(2)*s2b(3) - dir(3)*s2b(2) > zero       ) return
     
     ! PPM
     !-----------
  case(3)
     
     if(s2t(1) < zero .or. s2t(2) < zero .or. s2b(3) > zero) return ! on negative part of ray 
     if(e2b(1) > zero .or. e2b(2) > zero .or. e2t(3) < zero) return ! past length of ray      
     
     if ( dir(1)*s2t(2) - dir(2)*s2b(1) < zero .or.  &
          dir(1)*s2b(2) - dir(2)*s2t(1) > zero .or.  &
          dir(1)*s2b(3) - dir(3)*s2b(1) > zero .or.  &
          dir(1)*s2t(3) - dir(3)*s2t(1) < zero .or.  &
          dir(2)*s2t(3) - dir(3)*s2t(2) < zero .or.  &
          dir(2)*s2b(3) - dir(3)*s2b(2) > zero       ) return
     
     ! MMP
     !-----------
  case(4)
     
     if(s2b(1) > zero .or. s2b(2) > zero .or. s2t(3) < zero) return ! on negative part of ray 
     if(e2t(1) < zero .or. e2t(2) < zero .or. e2b(3) > zero) return ! past length of ray      
     
     if ( dir(1)*s2b(2) - dir(2)*s2t(1) < zero .or.  &
          dir(1)*s2t(2) - dir(2)*s2b(1) > zero .or.  &
          dir(1)*s2t(3) - dir(3)*s2t(1) > zero .or.  &
          dir(1)*s2b(3) - dir(3)*s2b(1) < zero .or.  &
          dir(2)*s2b(3) - dir(3)*s2b(2) < zero .or.  &
          dir(2)*s2t(3) - dir(3)*s2t(2) > zero       ) return
     
     
     ! PMP
     !-----------
  case(5)
     
     if(s2t(1) < zero .or. s2b(2) > zero .or. s2t(3) < zero) return ! on negative part of ray 
     if(e2b(1) > zero .or. e2t(2) < zero .or. e2b(3) > zero) return ! past length of ray      
     
     if ( dir(1)*s2t(2) - dir(2)*s2t(1) < zero .or.  &
          dir(1)*s2b(2) - dir(2)*s2b(1) > zero .or.  &
          dir(1)*s2b(3) - dir(3)*s2t(1) > zero .or.  &
          dir(1)*s2t(3) - dir(3)*s2b(1) < zero .or.  &
          dir(2)*s2b(3) - dir(3)*s2b(2) < zero .or.  &
          dir(2)*s2t(3) - dir(3)*s2t(2) > zero       ) return
     
     
     ! MPP
     !-----------
  case(6)
     
     if(s2b(1) > zero .or. s2t(2) < zero .or. s2t(3) < zero) return ! on negative part of ray 
     if(e2t(1) < zero .or. e2b(2) > zero .or. e2b(3) > zero) return ! past length of ray      
     
     if ( dir(1)*s2b(2) - dir(2)*s2b(1) < zero .or.  &
          dir(1)*s2t(2) - dir(2)*s2t(1) > zero .or.  &
          dir(1)*s2t(3) - dir(3)*s2t(1) > zero .or.  &
          dir(1)*s2b(3) - dir(3)*s2b(1) < zero .or.  &
          dir(2)*s2t(3) - dir(3)*s2b(2) < zero .or.  &
          dir(2)*s2b(3) - dir(3)*s2t(2) > zero       ) return
     
     ! PPP
     !-----------
  case(7)
     
     if(s2t(1) < zero .or. s2t(2) < zero .or. s2t(3) < zero) return ! on negative part of ray 
     if(e2b(1) > zero .or. e2b(2) > zero .or. e2b(3) > zero) return ! past length of ray      
     
     if ( dir(1)*s2t(2) - dir(2)*s2b(1) < zero .or.  &
          dir(1)*s2b(2) - dir(2)*s2t(1) > zero .or.  &
          dir(1)*s2b(3) - dir(3)*s2t(1) > zero .or.  &
          dir(1)*s2t(3) - dir(3)*s2b(1) < zero .or.  &
          dir(2)*s2t(3) - dir(3)*s2b(2) < zero .or.  &
          dir(2)*s2b(3) - dir(3)*s2t(2) > zero       ) return
     
  case default
     call rayError('ray class.')
     
  end select
  
  hit=.true.
  
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
