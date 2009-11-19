module qromb_mod
use nrtype
use nrutil
use trapzd_mod
use polint_mod
implicit none


interface qromb
   module procedure qromb_sp, qromb_dp
end interface


contains
 
  function qromb_sp(func,a,b)
    real(SP), intent(IN) :: a,b
    real(SP) :: qromb_sp

    interface
       function func(x)
         use nrtype
         real(SP), dimension(:), intent(IN) :: x
         real(SP), dimension(size(x)) :: func
       end function func
    end interface

    integer(I4B), parameter :: JMAX=25,JMAXP=JMAX+1,K=5,KM=K-1
    real(SP), parameter :: EPS=1.0e-6_sp
    real(SP), dimension(JMAXP) :: h,s
    real(SP) :: dqromb
    integer(I4B) :: j

    h(1)=1.0_sp
    do j=1,JMAX
       call trapzd(func,a,b,s(j),j)
       if (j >= K) then
          call polint(h(j-KM:j),s(j-KM:j),0.0_sp,qromb_sp,dqromb)
          if (abs(dqromb) <= EPS*abs(qromb_sp)) return
       end if
       s(j+1)=s(j)
       h(j+1)=0.25_sp*h(j)
    end do
    call nrerror('qromb_sp: too many steps')
  end function qromb_sp


  function qromb_dp(func,a,b)
    real(DP), intent(IN) :: a,b
    real(DP) :: qromb_dp

    interface
       function func(x)
         use nrtype
         real(DP), dimension(:), intent(IN) :: x
         real(DP), dimension(size(x)) :: func
       end function func
    end interface

    integer(I4B), parameter :: JMAX=25,JMAXP=JMAX+1,K=5,KM=K-1
    real(DP), parameter :: EPS=1.0e-6_dp
    real(DP), dimension(JMAXP) :: h,s
    real(DP) :: dqromb
    integer(I4B) :: j
    h(1)=1.0_dp
    do j=1,JMAX
       call trapzd(func,a,b,s(j),j)
       if (j >= K) then
          call polint(h(j-KM:j),s(j-KM:j),0.0_dp,qromb_dp,dqromb)
          if (abs(dqromb) <= EPS*abs(qromb_dp)) return
       end if
       s(j+1)=s(j)
       h(j+1)=0.25_dp*h(j)
    end do
    call nrerror('qromb_dp: too many steps')
  end function qromb_dp


end module qromb_mod
