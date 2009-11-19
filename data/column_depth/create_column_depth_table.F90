!> \file create_column_depth_table.F90

!> \brief Stored modules and program to make column depth look up tables 
!!
!<


!> Module for raytracing smoothing kernels
!===========================================================================
module kernel_functions_mod
use myf90_mod
implicit none
private

public :: set_hsmooth, set_bimpact
public :: monlat, monlatvol, monlatline
public :: tophat, tophatvol, tophatline
public :: pi, zero, one, two, three, four, six, eight, half


real(r8b), parameter :: pi = 3.141592653589793238462643383279502884197_r8b  
real(r8b), parameter :: zero = 0.0_r8b
real(r8b), parameter :: one = 1.0_r8b
real(r8b), parameter :: two = 2.0_r8b
real(r8b), parameter :: three = 3.0_r8b
real(r8b), parameter :: four = 4.0_r8b
real(r8b), parameter :: six = 6.0_r8b
real(r8b), parameter :: eight = 8.0_r8b
real(r8b), parameter :: half = 0.5_r8b

real(r8b), parameter :: coef1 = eight / pi
real(r8b), parameter :: coef2 = eight / pi * six
real(r8b), parameter :: coef3 = eight / pi * two


real(r8b), private :: hsmooth
real(r8b), private :: hsmooth3
real(r8b), private :: bimpact


contains



!> Call to set the smoothing length of the current kernel
!---------------------------------------------------------------------------
subroutine set_hsmooth(h)
real(r8b), intent(in) :: h
hsmooth = h
hsmooth3 = h*h*h
end subroutine set_hsmooth


!> Call to set the impact parameter of the current integral
!---------------------------------------------------------------------------
subroutine set_bimpact(b)
real(r8b), intent(in) :: b
bimpact = b
end subroutine set_bimpact


!> Monaghan and Lattanzio spline (vectorized)
!> returns W(r,h) eq. 4 from the Gadget-2 Paper
!---------------------------------------------------------------------------
function monlat(r) result(W)
real(r8b), dimension(:), intent(in) :: r
real(r8b), dimension(size(r)) :: W
real(r8b), dimension(size(r)) :: x 
real(r8b), dimension(size(r)) :: x2

integer :: i

x = r / hsmooth

if ( any( x < zero .or. x > one ) ) then
   write(*,*) "one of the x entries in monlatline is out of range"
   do i = 1,size(x)
      write(*,*) "i,x(i) = ", i, x(i)
   end do
   stop
end if

x2 = x * x

where ( x >= zero .and. x <= half ) 
   W = coef1 - coef2 * x2 + coef2 * x2 * x
elsewhere 
   W = coef3 * (one - x) * (one - x) * (one - x) 
end where

W = W / hsmooth3

end function monlat



!> This function returns W(s,h) where s is the distance along the third
!> side of the right triangle formed from the impact parameter (b) and the 
!> smoothing length (h) measured from the right angle.  It can be used to 
!> determine the line integral thru the kernel from -smax to smax where
!> smax = sqrt( h^2 - b^2 )
!---------------------------------------------------------------------------
function monlatline(s) result(W)
real(r8b), dimension(:), intent(in) :: s
real(r8b), dimension(size(s)) :: W
real(r8b), dimension(size(s)) :: r
real(r8b), dimension(size(s)) :: r2
real(r8b), dimension(size(s)) :: x 
real(r8b), dimension(size(s)) :: x2

integer :: i

r2 = bimpact * bimpact + s * s
r = sqrt( r2 )
x = r / hsmooth
 
if ( any( x < zero .or. x > one ) ) then
   write(*,*) "one of the x entries in monlatline is out of range"
   do i = 1,size(x)
      write(*,*) "i,x(i) = ", i, x(i)
   end do
   stop
end if

x2 = x * x

where ( x >= zero .and. x <= half ) 
   W = coef1 - coef2 * x2 + coef2 * x2 * x
elsewhere 
   W = coef3 * (one - x) * (one - x) * (one - x) 
end where

W = W / hsmooth3


end function monlatline


!> function that can be integrated to do the volume integral over monaghan 
!> and lattanzio spline
!---------------------------------------------------------------------------
function monlatvol(r) result(W)
real(r8b), dimension(:), intent(in) :: r
real(r8b), dimension(size(r)) :: W

   W = 4 * pi * r * r * monlat(r) 

end function monlatvol


!> tophat kernel
!---------------------------------------------------------------------------
function tophat(r) result(W)
real(r8b), dimension(:), intent(in) :: r
real(r8b), dimension(size(r)) :: W
real(r8b), dimension(size(r)) :: x 

integer :: i

x = r / hsmooth

if ( any( x < zero .or. x > one ) ) then
   write(*,*) "one of the x entries in tophat is out of range"
   do i = 1,size(x)
      write(*,*) "i,x(i) = ", i, x(i)
   end do
   stop
end if

W = three / (four * pi * hsmooth3)

end function tophat


!> line integral through tophat
!---------------------------------------------------------------------------
function tophatline(l) result(W)
real(r8b), dimension(:), intent(in) :: l
real(r8b), dimension(size(l)) :: W
real(r8b), dimension(size(l)) :: x 

integer :: i

x = sqrt( bimpact * bimpact + l * l )
x = x / hsmooth

if ( any( x < zero .or. x > one ) ) then
   write(*,*) "one of the x entries in tophat is out of range"
   do i = 1,size(x)
      write(*,*) "i,x(i) = ", i, x(i)
   end do
   stop
end if

W = three / (four * pi * hsmooth3)

end function tophatline


!>  function that can be integrated to do the tophat volume integral
!---------------------------------------------------------------------------
function tophatvol(r) result(W)
real(r8b), dimension(:), intent(in) :: r
real(r8b), dimension(size(r)) :: W

W = 4 * pi * r * r * tophat(r)

end function tophatvol



end module kernel_functions_mod



!> Program to calculate the column depth look up tables
!========================================================
program create_column_depth_table
use myf90_mod
use kernel_functions_mod
use nr, only: qromb
implicit none



real(r8b), parameter :: bmin = zero
real(r8b), parameter :: bmax = one
integer(i4b), parameter :: bevals = 51
real(r8b), parameter :: db = (bmax - bmin) / (bevals-1)

real(r8b) :: hsml
real(r8b) :: bimp
real(r8b) :: kint
real(r8b) :: smax

integer(i4b) :: i

real(r8b) :: res(1)


hsml = one
call set_hsmooth(hsml)

bimp = zero
call set_bimpact(bimp)


write(*,*)
write(*,*) "==============================================================="
write(*,*) "     Calculating Impact Paramaeter -> Column Depth Table       "
write(*,*) "==============================================================="
write(*,*) 
write(*,*) " checking simple evaluations of kernel "
write(*,*) "---------------------------------------"


res = monlat( (/zero/) ) 
write(*,*) " monlat(0.0) = ", res

res = monlat( (/one/) ) 
write(*,*) " monlat(1.0) = ", res
write(*,*) 

res = monlatline( (/zero/) ) 
write(*,*) " monlatline(0.0) = ", res

res = monlatline( (/one/) ) 
write(*,*) " monlatline(1.0) = ", res
write(*,*) 


res = tophat( (/zero/) ) 
write(*,*) " tophat(0.0) = ", res

res = tophat( (/one/) ) 
write(*,*) " tophat(1.0) = ", res
write(*,*) 


res = tophatline( (/zero/) ) 
write(*,*) " tophatline(0.0) = ", res

res = tophatline( (/one/) ) 
write(*,*) " tophatline(1.0) = ", res
write(*,*) 



write(*,*)  
write(*,*) " checking kernel normalizations (volume integrals) "
write(*,*) "---------------------------------------------------"

kint = qromb(monlatvol, zero, hsml)
write(*,*) " monlat volume integral = ", kint 

kint = qromb(tophatvol, zero, hsml)
write(*,*) " tophat volume integral = ", kint 

write(*,*) 
write(*,*) 
write(*,*) " calculating line integrals at b=0 (as reference) "

kint = qromb(monlat, zero, hsml)
write(*,*) " monlat line integral through center = ", 2.0 * kint 

kint = qromb(tophat, zero, hsml)
write(*,*) " tophat line integral through center = ", 2.0 * kint 
write(*,*) 


bimp = zero
call set_bimpact(bimp)

kint = qromb(monlatline, zero, hsml)
write(*,*) " monlat line check = ", 2.0 * kint 

kint = qromb(tophatline, zero, hsml)
write(*,*) " tophat line check = ", 2.0 * kint 
write(*,*) 


!----------------------------------------------
hsml = one
call set_hsmooth(hsml)

open(unit=10,file="latmon_b2cd_table.txt")
open(unit=11,file="tophat_b2cd_table.txt")

write(10,*) bevals
write(11,*) bevals
do i = 0,bevals-1

   bimp = bmin + i * db

   call set_bimpact(bimp)
   smax = sqrt( hsml*hsml - bimp*bimp )

   kint = qromb(monlatline, zero, smax)
   write(10,*) 2.0 * kint 

   kint = qromb(tophatline, zero, smax)
   write(11,*) 2.0 * kint 

end do
close(10)
close(11)

 
end program create_column_depth_table
