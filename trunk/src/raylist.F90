!> \file raylist.f90

!> \brief The raylist module 
!!
!<

module raylist_mod
use myf90_mod
use ray_mod, only: ray_type
use particle_system_mod, only: particle_system_type
use particle_system_mod, only: particle_type
use particle_system_mod, only: transformation_type
use oct_tree_mod, only: oct_tree_type
implicit none

 private
 public :: intersection_type
 public :: trace_ray
 public :: raylist_type
 public :: prepare_raysearch
 public :: kill_raylist


 integer,parameter :: MAX_RAYLIST_LENGTH=1000000  !< default maximum
 integer,parameter :: ndim=3                      !< number of dimensions



!> holds a particle index, an impact parameter, and a distance along a ray
!! usefull for the ionization and temperature updating
! -----------------------------------------------------------------------
   type intersection_type
      integer(i8b) :: pindx   !< particle index
      real :: b          !< impact parameter
      real :: d          !< distance along ray
   end type intersection_type
 
!> grouping of all things ray + impacts
!--------------------------------------- 
   type raylist_type
      integer :: nnb              !< number of intersections
      integer :: maxnnb           !< maximum number of intersections
      integer :: lastnnb          !< intersection where we stopped
      integer(i8b) :: searchcell  !< index of cell being searched
      logical :: reuseable        !< is this ray reusable?
      integer :: searchimage      !< which transformation of the particles?
      integer :: nsearchimages    !< how many images to search
      type(ray_type) :: ray       !< ray 
      type(transformation_type) :: trafo(3**ndim)  !< transformations  
      type(intersection_type), allocatable :: intersection(:) !< ray/par 
   end type raylist_type



contains


 
!> set intersection values
!-----------------------------------------
 function set_intersection(ray,part,i) result(t)
 use ray_mod, only: dist2ray

   type(intersection_type) :: t  !< the intersection
   type(ray_type) :: ray         !< the ray
   type(particle_type) :: part   !< the particle
   integer(i8b)  :: i            !< the particle index
   real(r8b) :: d
   
     t%pindx=i
     t%b=sqrt(dist2ray(ray,part%pos,d))
     t%d=d

 end function set_intersection
 
!> initialize raylist variables and set search images
!------------------------------------------------------
 subroutine prepare_raysearch(psys,raylist)

   type(particle_system_type) psys !< particle system
   type(raylist_type) raylist      !< ray list
  
     call make_raylist(MAX_RAYLIST_LENGTH,raylist)
     call setsearchimages(psys,raylist)

 end subroutine

!> check all of the search images for intersection
!---------------------------------------------------
 subroutine fullsearch(psys,searchtree,raylist)

   type(particle_system_type) :: psys  !< particle system
   type(oct_tree_type) :: searchtree   !< oct-tree to search
   type(raylist_type) raylist          !< raylist
 
     if(raylist%searchimage.EQ.0) call raylistError('raylist init.')
     do while(raylist%searchimage.LE.raylist%nsearchimages)   
        call raysearch(psys,searchtree,raylist)
        if(raylist%searchcell.NE.0) return 
        raylist%searchimage=raylist%searchimage+1
        raylist%searchcell=1
     enddo
     raylist%searchcell=0

 end subroutine fullsearch

!> checks a single search image for intersection
!-----------------------------------------------
 subroutine raysearch(psys,searchtree,raylist)
 use ray_mod, only: part_intersection, cell_intersection, transform_ray

   type(particle_system_type) :: psys          !< particle system
   type(oct_tree_type), target :: searchtree   !< oct-tree to search
   type(raylist_type) :: raylist               !< raylist

   type(oct_tree_type),pointer :: tree
   type(ray_type) :: curay
   integer(i8b) :: this, daughter, next
   integer(i8b) :: i, si, orderindx

!   integer :: k  

     tree=>searchtree
     si=raylist%searchimage
     call transform_ray(curay,raylist%ray,raylist%trafo(si))
     next=raylist%searchcell

     do while(next.ne.0)
        this=next
        daughter=tree%cell(this)%daughter
        next=tree%cell(this)%next

        if(daughter.EQ.0) then
           if(raylist%nnb+tree%cell(next)%start-tree%cell(this)%start .GT. &
           raylist%maxnnb) then 

              raylist%reuseable=.FALSE.
              raylist%searchcell=this
              return

           endif

           do i = tree%cell(this)%start, tree%cell(next)%start-1
              orderindx = tree%partorder(i)
              if(part_intersection(curay,psys%par(orderindx))) then
                 raylist%nnb = raylist%nnb+1
                 raylist%intersection(raylist%nnb)= &
                      set_intersection(curay,psys%par(orderindx),orderindx)
              endif
           enddo

        else
           if(cell_intersection(curay,tree%cell(this))) next=daughter  
        endif
     enddo

     raylist%searchcell=next

 end subroutine raysearch

!> initializes values in the raylist
!-------------------------------------------
 subroutine make_raylist(maxnnb,raylist,ray)
 use ray_mod, only: set_ray

   integer            :: maxnnb       !< maximum number of intersections
   type(raylist_type)      :: raylist !< the raylist
   type(ray_type),optional :: ray     !< the ray
   real(r8b) :: start(ndim)
   real(r8b) :: dir(ndim)
  
     raylist%nnb=0
     raylist%maxnnb=maxnnb
     raylist%searchcell=1
     raylist%reuseable=.FALSE.
     raylist%searchimage=1
     raylist%nsearchimages=1
     raylist%trafo(1)%fac=1
     raylist%trafo(1)%shift=0 

     start = (/0.,0.,0./)
     dir   = (/1.,0.,0./) 
  
     call set_ray(raylist%ray, start, dir)
     allocate(raylist%intersection(maxnnb))
     if(present(ray)) raylist%ray=ray

 end subroutine make_raylist

!> reset an already initialized raylist
!---------------------------------------
 subroutine reset_raylist(raylist,ray)

   type(raylist_type) :: raylist   !< the raylist
   type(ray_type), optional :: ray !< the ray
  
     raylist%nnb=0
     raylist%searchcell=1
     raylist%searchimage=1
     raylist%reuseable=.FALSE.
     if(present(ray)) raylist%ray=ray

 end subroutine reset_raylist

!> kill a raylist
!---------------------------------
 subroutine kill_raylist(raylist)
 use ray_mod, only: set_ray

   type(raylist_type) :: raylist !< the raylist to kill
   real(r8b) :: start(ndim)
   real(r8b) :: dir(ndim) 
  
     raylist%nnb=0
     raylist%maxnnb=MAX_RAYLIST_LENGTH
     raylist%searchcell=0
     start = (/0.,0.,0./)
     dir   = (/1.,0.,0./) 
     call set_ray(raylist%ray, start, dir)
     deallocate(raylist%intersection)

 end subroutine kill_raylist

!> create the transformations for the searchimages
!-------------------------------------------------
 subroutine setsearchimages(psys,raylist)

   type(raylist_type)     :: raylist
   type(particle_system_type) :: psys
   integer :: i,j,k(ndim)
   real :: top(ndim),bot(ndim)
   integer :: bbound(ndim),tbound(ndim),l,nsearchimages
   real :: lx(ndim),hx(ndim)
   logical :: boxexists

     top=psys%box%top
     bot=psys%box%bot
     bbound=psys%box%bbound
     tbound=psys%box%tbound
 

     nsearchimages=0

     do i=0,3**ndim-1
        l=i

        do j=1,ndim
           k(j)=mod(l,3)-1 
           l=l/3
        enddo

        boxexists=.TRUE.


        do j=1,ndim
           lx(j)=bot(j)+k(j)*(top(j)-bot(j))
           hx(j)=top(j)+k(j)*(top(j)-bot(j))
           if(k(j).eq.-1.AND.bbound(j).EQ.0) boxexists=.FALSE.  
           if(k(j).eq.1.AND.tbound(j).EQ.0) boxexists=.FALSE.  
        enddo

        if(boxexists) then
           nsearchimages=nsearchimages+1
  
           do j=1,ndim
              raylist%trafo(nsearchimages)%fac(j)=1
              raylist%trafo(nsearchimages)%shift(j)=0.
              if(k(j).eq.-1) then
                 raylist%trafo(nsearchimages)%fac(j)=bbound(j)
                 raylist%trafo(nsearchimages)%shift(j) = &
                   (-3*raylist%trafo(nsearchimages)%fac(j)+1)*bot(j)/2 + &
                    (1+raylist%trafo(nsearchimages)%fac(j))*top(j)/2
              endif
              if(k(j).eq.1) then
                 raylist%trafo(nsearchimages)%fac(j)=tbound(j)
                 raylist%trafo(nsearchimages)%shift(j) = &
                   (-3*raylist%trafo(nsearchimages)%fac(j)+1)*top(j)/2 + &
                    (1+raylist%trafo(nsearchimages)%fac(j))*bot(j)/2
              endif

           enddo

        endif ! boxexists


     enddo

     raylist%nsearchimages=nsearchimages


 end subroutine setsearchimages


!> error handling
!---------------------------------
 subroutine raylistError(string,i)
   character(*) :: string  !< error message
   integer, optional :: i  !< error number

     print*,' Error detected:'
  
     if(present(i)) then
        print*,string,i
     else
        print*,string
     endif
  
     stop

 end subroutine raylistError

!> sorts a raylist according to the distance along the ray with the particles 
!! closest to the origin of the ray first in the list.   The corresponding 
!! changes are made to pindx and b
!---------------------------------------------------------------------------
 subroutine sort3_raylist(raylist)
 use m_mrgrnk, only: mrgrnk
 implicit none

   type(raylist_type) :: raylist     !< raylist to sort

   integer(i8b) :: N 
   integer(i8b) :: indexx(raylist%nnb)
   real(r8b) :: darr(raylist%nnb)

   N = raylist%nnb
   darr(1:N)=raylist%intersection(1:N)%d
   call mrgrnk(darr,indexx)

   raylist%intersection(1:N)%d=raylist%intersection(indexx(1:N))%d  
   raylist%intersection(1:N)%b=raylist%intersection(indexx(1:N))%b 
   raylist%intersection(1:N)%pindx=raylist%intersection(indexx(1:N))%pindx  

 end subroutine sort3_raylist

!> given a ray creates a raylist with intersections
!------------------------------------------------------
 subroutine trace_ray(ray,raylist,psys,searchtree) 
   type(ray_type), intent(in) :: ray               !< the ray to trace
   type(raylist_type), intent(inout) :: raylist    !< the returned raylist
   type(particle_system_type), intent(in) :: psys  !< the particle system
   type(oct_tree_type), intent(in) :: searchtree   !< the oct-tree to search
   
   if (raylist%maxnnb.ne.MAX_RAYLIST_LENGTH) then      
      call prepare_raysearch(psys,raylist)
   end if
   call reset_raylist(raylist,ray)
   call fullsearch(psys,searchtree,raylist) 
   call sort3_raylist(raylist)
 end subroutine trace_ray




end module raylist_mod
