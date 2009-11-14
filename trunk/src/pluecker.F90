!> \file pluecker.f90

!> \brief the pluecker test module 
!!
!<

module pluecker_mod
use myf90_mod
implicit none
private
  
public :: pluecker


contains

!> pluecker test for line segment / cell intersection
!-----------------------------------------------------    
 function pluecker(dir,class,bot,top) result(x)
   real(r4b) :: dir(3)      !< direction vector of ray
   integer(i4b) :: class    !< class of ray (determines octant)
   real(r8b) :: bot(3)      !< lower cell corner
   real(r8b) :: top(3)      !< upper cell corner
   logical :: x             !< true or false result
   
   x = .FALSE.

   select case(class)
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
      stop "ray class in pluecker test out of bounds" 
   end select
   
   x = .TRUE.
   
 end function pluecker
 

 
end module pluecker_mod
