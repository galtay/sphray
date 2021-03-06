	FUNCTION qromb_sp(func,a,b)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE nr, ONLY : polint,trapzd
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP) :: qromb_sp
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: JMAX=25,JMAXP=JMAX+1,K=5,KM=K-1
	REAL(SP), PARAMETER :: EPS=1.0e-6_sp
	REAL(SP), DIMENSION(JMAXP) :: h,s
	REAL(SP) :: dqromb
	INTEGER(I4B) :: j
	h(1)=1.0
	do j=1,JMAX
		call trapzd(func,a,b,s(j),j)
		if (j >= K) then
			call polint(h(j-KM:j),s(j-KM:j),0.0_sp,qromb_sp,dqromb)
			if (abs(dqromb) <= EPS*abs(qromb_sp)) RETURN
		end if
		s(j+1)=s(j)
		h(j+1)=0.25_sp*h(j)
	end do
	call nrerror('qromb_sp: too many steps')
        END FUNCTION qromb_sp

	FUNCTION qromb_dp(func,a,b)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE nr, ONLY : polint,trapzd
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: a,b
	REAL(DP) :: qromb_dp
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		REAL(DP), DIMENSION(:), INTENT(IN) :: x
		REAL(DP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: JMAX=25,JMAXP=JMAX+1,K=5,KM=K-1
	REAL(DP), PARAMETER :: EPS=1.0e-12_dp
	REAL(DP), DIMENSION(JMAXP) :: h,s
	REAL(DP) :: dqromb
	INTEGER(I4B) :: j

	h(1)=1.0
	do j=1,JMAX
		call trapzd(func,a,b,s(j),j)
		if (j >= K) then
			call polint(h(j-KM:j),s(j-KM:j),0.0_dp,qromb_dp,dqromb)
			if (abs(dqromb) <= EPS*abs(qromb_dp)) RETURN
		end if
		s(j+1)=s(j)
		h(j+1)=0.25_dp*h(j)
	end do
	call nrerror('qromb_dp: too many steps')
	END FUNCTION qromb_dp
