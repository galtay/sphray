program test_rates_vs_iontables
use atomic_rates_mod
use ion_table_class
implicit none

character(clen) :: rates_file
character(clen) :: ion_file

real :: z, logT, logD
real :: dlogT, logTmin, logTmax
real :: x_i, x_h, x_c
real*8 :: T, D

real :: RR, QQ, PP, det
real :: GHI, CI, RC, y
integer :: i, nsteps


type(ion_table_type) :: itable            !< CLOUDY ionization table

type(atomic_rates_table_type) :: rtable   !< rates read in from file
type(atomic_rates_type) :: k





rates_file = "atomic_rates_Cen.txt"
ion_file = "../ionization_tables/h1.hdf5"


call read_atomic_rates_file(rtable, rates_file)
call read_ion_table_file(ion_file, itable)


z = 2.013
logD = 1.0e-3
D = 10**(logD)

y=0.00
GHI = itable%ihead%ispec%gammaHI(24)

write(*,*) "z(24) = ", itable%z(24)
write(*,*) "G(24) = ", GHI


logTmin=3.0
logTmax=6.0
nsteps=100
dlogT= (logTmax - logTmin)/(nsteps-1)


open(unit=10, file="rate_comp.txt")
do i = 0,nsteps-1
   logT = logTmin + dlogT * i
   T = 10**(logT)

   x_i = interpolate_ion_table( itable, z, logT, logD )


   call get_atomic_rates( T, rtable, k )
   CI = Hui_HI_col_ion(T)

   RC = Hui_HII_recombA(T)
   RR = -(CI + RC) * D
   QQ = CI * D - GHI - (CI + RC) * D * y
   PP = GHI + CI * D * y
   det = QQ**2 - 4 * RR * PP
   x_h = 1.0 - ( -QQ - sqrt(det) ) / (2*RR)

      
   RC = Cen_HII_recombA(T)
   RR = -(CI + RC) * D
   QQ = CI * D - GHI - (CI + RC) * D * y
   PP = GHI + CI * D * y
   det = QQ**2 - 4 * RR * PP
   x_c = 1.0 - ( -QQ - sqrt(det) ) / (2*RR)


   write(10,*) logT, x_h,x_c,x_i

end do
close(10)












end program test_rates_vs_iontables
