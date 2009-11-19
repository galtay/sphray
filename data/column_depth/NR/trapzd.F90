module trapzd_mod
use nrtype
use nrutil
implicit none


interface trapzd
   module procedure trapzd_sp, trapzd_dp
end interface


contains

  subroutine trapzd_sp(func,a,b,s,n)
    real(SP), intent(IN) :: a,b
    real(SP), intent(INOUT) :: s
    integer(I4B), intent(IN) :: n

    interface
       function func(x)
         use nrtype
         real(SP), dimension(:), intent(IN) :: x
         real(SP), dimension(size(x)) :: func
       end function func
    end interface

    real(SP) :: del,fsum
    integer(I4B) :: it

    if (n == 1) then
       s=0.5_sp*(b-a)*sum(func( (/ a,b /) ))
    else
       it=2**(n-2)
       del=(b-a)/it
       fsum=sum(func(arth(a+0.5_sp*del,del,it)))
       s=0.5_sp*(s+del*fsum)
    end if

  end subroutine trapzd_sp
  
  subroutine trapzd_dp(func,a,b,s,n)
    real(DP), intent(IN) :: a,b
    real(DP), intent(INOUT) :: s
    integer(I4B), intent(IN) :: n

    interface
       function func(x)
         use nrtype
         real(DP), dimension(:), intent(IN) :: x
         real(DP), dimension(size(x)) :: func
       end function func
    end interface

    real(DP) :: del,fsum
    integer(I4B) :: it

    if (n == 1) then
       s=0.5_dp*(b-a)*sum(func( (/ a,b /) ))
    else
       it=2**(n-2)
       del=(b-a)/it
       fsum=sum(func(arth(a+0.5_dp*del,del,it)))
       s=0.5_dp*(s+del*fsum)
    end if
  end subroutine trapzd_dp

  
end module trapzd_mod
