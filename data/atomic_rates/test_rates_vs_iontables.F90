program test_rates_vs_iontables
use atomic_rates_mod
use ion_table_class
implicit none

character(clen) :: rates_file
character(clen) :: ion_file

real :: z, logT, logD(2)
real :: dlogT, logTmin, logTmax
real*8 :: x_i, x_hA, x_cA, x_hB, x_cB
real*8 :: T, D(2)

real*8 :: RR, QQ, PP, det
real*8 :: GHI, CI, RC, y
integer :: i, j, nsteps

character(clen) :: files(2)

type(ion_table_type) :: itable            !< CLOUDY ionization table

type(atomic_rates_table_type) :: rtable   !< rates read in from file
type(atomic_rates_type) :: k





rates_file = "atomic_rates_Cen.txt"
ion_file = "../ionization_tables/h1.hdf5"


call read_atomic_rates_file(rtable, rates_file)
call read_ion_table_file(ion_file, itable)

files(1) = "rate_comp_D1n5.txt"
files(2) = "rate_comp_D1p2.txt"

y=0.00
GHI = itable%ihead%ispec%gammaHI(24)
z = 2.013

write(*,*) "z(24) = ", itable%z(24)
write(*,*) "G(24) = ", GHI


logD = (/ -5.0, 2.0 /)
D = 10**(logD)


logTmin=3.0
logTmax=6.0
nsteps=100
dlogT= (logTmax - logTmin)/(nsteps-1)




do j = 1,2

   open(unit=10, file=files(j))

   do i = 0,nsteps-1
      logT = logTmin + dlogT * i
      T = 10**(logT)
      
      x_i = interpolate_ion_table( itable, z, logT, logD(j) )


      call get_atomic_rates( T, rtable, k )
      CI = Hui_HI_col_ion(T)

      RC = Hui_HII_recombA(T)
      RR = -(CI + RC) * D(j)
      QQ = CI * D(j) - GHI - (CI + RC) * D(j) * y
      PP = GHI + CI * D(j) * y
      det = QQ**2 - 4 * RR * PP
      x_hA = 1.0 - ( -QQ - sqrt(det) ) / (2*RR)
      
      
      RC = Cen_HII_recombA(T)
      RR = -(CI + RC) * D(j)
      QQ = CI * D(j) - GHI - (CI + RC) * D(j) * y
      PP = GHI + CI * D(j) * y
      det = QQ**2 - 4 * RR * PP
      x_cA = 1.0 - ( -QQ - sqrt(det) ) / (2*RR)
      

      RC = Hui_HII_recombB(T)
      RR = -(CI + RC) * D(j)
      QQ = CI * D(j) - GHI - (CI + RC) * D(j) * y
      PP = GHI + CI * D(j) * y
      det = QQ**2 - 4 * RR * PP
      x_hB = 1.0 - ( -QQ - sqrt(det) ) / (2*RR)
      
      
      write(10,"(5ES12.5)") logT, x_hA, x_cA, x_i, x_hB

   end do
   close(10)
   
end do













end program test_rates_vs_iontables
