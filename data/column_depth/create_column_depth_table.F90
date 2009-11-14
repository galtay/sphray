! this file can be used to generate column depth look up tables
!===============================================================

module kernel_functions
use myf90_mod
private

public :: set_hsmooth, set_bimpact
public :: monlatspline, monlatvol, monlatline
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

real(r8b) :: hsmooth
real(r8b) :: bimpact


contains

! so that the kernel functions can be integrated using the 
! NR routines the smoothing length and impact parameter need 
! to be set through calls to these subroutines
!---------------------------------------------------------------
subroutine set_hsmooth(h)
real(r8b), intent(in) :: h
hsmooth = h
end subroutine set_hsmooth

subroutine set_bimpact(b)
real(r8b), intent(in) :: b
bimpact = b
end subroutine set_bimpact


! monaghan and lattanzio spline (vectorized)
! this returns W(r,h) eq. 4 from the Gadget-2 Paper
!----------------------------------------------------------
function monlatspline(r) result(W)
real(r8b), dimension(:), intent(in) :: r
real(r8b), dimension(size(r)) :: W
real(r8b), dimension(size(r)) :: x 

x = r / hsmooth

do i = 1,size(x)
   if (x(i) < zero .or. x(i) > one) then
      write(*,*) "one of the x entries in monlatspline is out of range"
      write(*,*) "i,x(i) = ", i, x(i)
      stop
   end if
end do

where ( x >= zero .and. x <= half ) 
   W = one - six * x * x + six * x * x * x
elsewhere 
   W = two * (one - x) * (one - x) * (one - x) 
end where

W = W * eight / (pi * hsmooth * hsmooth * hsmooth)

end function monlatspline

! this function is the one that can be used to calculate a line
! integral through the monaghan and lattanzio spline kernel.
! The input is s, which is the length of the third side of the 
! right triangle formed along with the impact parameter and radius. 
! Once an impact parameter and smoothing length is set the line integral
! extends from -smax to smax where smax = sqrt( h^2 - b^2 )
!-------------------------------------------------------------------

function monlatline(s) result(W)
real(r8b), dimension(:), intent(in) :: s
real(r8b), dimension(size(s)) :: r
real(r8b), dimension(size(s)) :: x 
real(r8b), dimension(size(s)) :: W

r = sqrt( bimpact * bimpact + s * s )
x = r / hsmooth

do i = 1,size(x)
   if (x(i) < zero .or. x(i) > one) then
      write(*,*) "one of the x entries in monlatspline is out of range"
      write(*,*) "i,x(i) = ", i, x(i)
      stop
   end if
end do

where ( x >= zero .and. x <= half ) 
   W = one - six * x * x + six * x * x * x
elsewhere 
   W = two * (one - x) * (one - x) * (one - x) 
end where

W = W * eight / (pi * hsmooth * hsmooth * hsmooth)

end function monlatline


! volume integral over monaghan and lattanzio spline
!----------------------------------------------------
function monlatvol(r) result(W)
real(r8b), dimension(:), intent(in) :: r
real(r8b), dimension(size(r)) :: W

   W = 4 * pi * r * r * monlatspline(r) 

end function monlatvol


! tophat kernel
!----------------------------------------------------
function tophat(r) result(W)
real(r8b), dimension(:), intent(in) :: r
real(r8b), dimension(size(r)) :: W
real(r8b), dimension(size(r)) :: x 

x = r / hsmooth

do i = 1,size(x)
   if (x(i) < zero .or. x(i) > one) then
      write(*,*) "one of the x entries in monlatspline is out of range"
      write(*,*) "i,x(i) = ", i, x(i)
      stop
   end if
end do

W = one

W = W * three / (four * pi * hsmooth * hsmooth * hsmooth)

end function tophat

! line integral through tophat
!----------------------------------------------------
function tophatline(l) result(W)
real(r8b), dimension(:), intent(in) :: l
real(r8b), dimension(size(l)) :: W
real(r8b), dimension(size(l)) :: x 

x = sqrt( bimpact * bimpact + l * l )
x = x / hsmooth

do i = 1,size(x)
   if (x(i) < zero .or. x(i) > one) then
      write(*,*) "one of the x entries in tophat is out of range"
      write(*,*) "i,x(i) = ", i, x(i)
      stop
   end if
end do

W = one
W = W * three / (four * pi * hsmooth * hsmooth * hsmooth)

end function tophatline

! tophat volume kernel
!----------------------------------------------------
function tophatvol(r) result(W)
real(r8b), dimension(:), intent(in) :: r
real(r8b), dimension(size(r)) :: W

W = 4 * pi * r * r * tophat(r)

end function tophatvol



end module kernel_functions


program create_column_depth_table
use myf90_mod
use kernel_functions
use nr, only: qromb
implicit none



real(r8b), parameter :: bmin = zero
real(r8b), parameter :: bmax = one
integer(i4b), parameter :: bevals = 1001
real(r8b), parameter :: db = (bmax - bmin) / (bevals-1)

real(r8b) :: hsml
real(r8b) :: bimp
real(r8b) :: kint
real(r8b) :: lmax

integer(i4b) :: i


write(*,*)
write(*,*) "==============================================================="
write(*,*) "     Calculating Impact Paramaeter -> Column Depth Table       "
write(*,*) "==============================================================="
write(*,*)  
write(*,*) " checking kernel normalizations (volume integrals)"

hsml = one
call set_hsmooth(hsml)

kint = qromb(monlatvol, zero, hsml)
write(*,*) " monlat volume integral = ", kint 

kint = qromb(tophatvol, zero, hsml)
write(*,*) " tophat volume integral = ", kint 

write(*,*) 
write(*,*) 
write(*,*) " calculating line integrals at b=0 (as reference) "

kint = qromb(monlatspline, zero, hsml)
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
   lmax = sqrt( hsml*hsml - bimp*bimp )

   kint = qromb(monlatline, zero, lmax)
   write(10,*) 2.0 * kint 

   kint = qromb(tophatline, zero, lmax)
   write(11,*) 2.0 * kint 

end do
close(10)
close(11)

 
end program create_column_depth_table
