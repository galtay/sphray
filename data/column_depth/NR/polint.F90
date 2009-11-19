module polint_mod
use nrtype
use nrutil
implicit none


interface polint
   module procedure polint_sp, polint_dp
end interface

contains

  subroutine polint_sp(xa,ya,x,y,dy)
    real(SP), dimension(:), intent(IN) :: xa,ya
    real(SP), intent(IN) :: x
    real(SP), intent(OUT) :: y,dy
    integer(I4B) :: m,n,ns
    real(SP), dimension(size(xa)) :: c,d,den,ho
    n=assert_eq(size(xa),size(ya),'polint_sp')
    c=ya
    d=ya
    ho=xa-x
    ns=iminloc(abs(x-xa))
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
       den(1:n-m)=ho(1:n-m)-ho(1+m:n)
       if (any(den(1:n-m) == 0.0)) &
            call nrerror('polint_sp: calculation failure')
       den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
       d(1:n-m)=ho(1+m:n)*den(1:n-m)
       c(1:n-m)=ho(1:n-m)*den(1:n-m)
       if (2*ns < n-m) then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       end if
       y=y+dy
    end do
  end subroutine polint_sp


  subroutine polint_dp(xa,ya,x,y,dy)
    real(DP), dimension(:), intent(IN) :: xa,ya
    real(DP), intent(IN) :: x
    real(DP), intent(OUT) :: y,dy
    integer(I4B) :: m,n,ns
    real(DP), dimension(size(xa)) :: c,d,den,ho
    n=assert_eq(size(xa),size(ya),'polint_dp')
    c=ya
    d=ya
    ho=xa-x
    ns=iminloc(abs(x-xa))
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
       den(1:n-m)=ho(1:n-m)-ho(1+m:n)
       if (any(den(1:n-m) == 0.0_dp)) &
            call nrerror('polint_dp: calculation failure')
       den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
       d(1:n-m)=ho(1+m:n)*den(1:n-m)
       c(1:n-m)=ho(1:n-m)*den(1:n-m)
       if (2*ns < n-m) then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       end if
       y=y+dy
    end do
  end subroutine polint_dp


end module polint_mod
