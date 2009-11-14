program write_rates
use myf90_mod
use HuiGnedinAtomicRates
use CenAtomicRates
implicit none

real(r8b), parameter :: T1 = 1.0d0
real(r8b), parameter :: T2 = 1.0d9
real(r8b), parameter :: logT1 = log10(T1)
real(r8b), parameter :: logT2 = log10(T2)
integer(i8b), parameter :: Tbins = 2251
real(r8b), parameter :: dTlog = (logT2 - logT1) / real(Tbins-1)
integer(i8b) :: i
real(r8b) :: T
real(r8b) :: logT

character(200) :: ratesfile

character(14), parameter :: names(25) = (/ "temperature" , &
                                          "HIci",    "HeIci",    "HeIIci", &
                                          "HIIrcA",  "HeIIrcA",  "HeIIIrcA", &
                                          "HIIrcB",  "HeIIrcB",  "HeIIIrcB", &
                                          "HeDrc", &
                                          "HIIrccA", "HeIIrccA", "HeIIIrccA", &
                                          "HIIrccB", "HeIIrccB", "HeIIIrccB", &
                                          "HeDrcc", &
                                          "HIcic",   "HeIcic",   "HeIIcic", &
                                          "He23cic", &
                                          "HIcec",   "HeIcec",   "HeIIcec"  /)

!--------------
!> atomic rates
type atomic_rates_type
   real(r4b) :: logT      !< log10 temperature

   real(r4b) :: HIci      !< HI   collisional ionization rate 
   real(r4b) :: HeIci     !< HeI  collisional ionization rate 
   real(r4b) :: HeIIci    !< HeII collisional ionization rate 

   real(r4b) :: HIIrcA    !< HII   recombination rate (case A)
   real(r4b) :: HeIIrcA   !< HeII  recombination rate (case A)
   real(r4b) :: HeIIIrcA  !< HeIII recombination rate (case A)

   real(r4b) :: HIIrcB    !< HII   recombination rate (case B)
   real(r4b) :: HeIIrcB   !< HeII  recombination rate (case B)
   real(r4b) :: HeIIIrcB  !< HeIII recombination rate (case B)
   real(r4b) :: HeDrc     !< dielectronic He recombination rate   

   real(r4b) :: HIIrccA   !< HII   recombination cooling rate (case A)
   real(r4b) :: HeIIrccA  !< HeII  recombination cooling rate (case A)
   real(r4b) :: HeIIIrccA !< HeIII recombination cooling rate (case A)

   real(r4b) :: HIIrccB   !< HII   recombination cooling rate (case B)
   real(r4b) :: HeIIrccB  !< HeII  recombination cooling rate (case B)
   real(r4b) :: HeIIIrccB !< HeIII recombination cooling rate (case B)
   real(r4b) :: HeDrcc    !< dielectronic He recombination cooling rate
   
   real(r4b) :: HIcic     !< HI   collisional ionization cooling rate
   real(r4b) :: HeIcic    !< HeI  collisional ionization cooling rate
   real(r4b) :: HeIIcic   !< HeII collisional ionization cooling rate
   real(r4b) :: He23cic   !< He23 collisional ionization cooling rate
   
   real(r4b) :: HIcec     !< HI   collisional excitation cooling rate
   real(r4b) :: HeIcec    !< HeI  collisional excitation cooling rate
   real(r4b) :: HeIIcec   !< HeII collisional excitation cooling rate
      
end type atomic_rates_type



type(atomic_rates_type) :: k

write(*,*) "dTlog = ", dTlog
write(*,*) "logT2 - logT2 = ", logT2 - logT1


ratesfile = "atomic_rates.txt"
open(unit=10,file=ratesfile)

150 format(2ES12.5,I8,ES12.5)
write(10,150) logT1, logT2, Tbins, dTlog
200 format(1X,25(A))
write(10,200) names


do i = 0,Tbins-1

   logT = logT1 + i*dTlog
   T = 10**logT

   k%logT = logT

   k%HIci   = Cen_HI_col_ion(T)   
   k%HeIci  = Cen_HeI_col_ion(T)  
   k%HeIIci = Cen_HeII_col_ion(T) 

   k%HIIrcA   = Hui_HII_recombA(T)
   k%HeIIrcA  = Hui_HeII_recombA(T) 
   k%HeIIIrcA = Hui_HeIII_recombA(T) 

   k%HIIrcB   = Hui_HII_recombB(T)
   k%HeIIrcB  = Hui_HeII_recombB(T) 
   k%HeIIIrcB = Hui_HeIII_recombB(T)      
   k%HeDrc    = Hui_He_dielec_recomb(T)

   k%HIIrccA   = Hui_HII_rec_coolA(T)
   k%HeIIrccA  = Hui_HeII_rec_coolA(T) 
   k%HeIIIrccA = Hui_HeIII_rec_coolA(T) 

   k%HIIrccB   = Hui_HII_rec_coolB(T)
   k%HeIIrccB  = Hui_HeII_rec_coolB(T) 
   k%HeIIIrccB = Hui_HeIII_rec_coolB(T)      
   k%HeDrcc    = Hui_He_dielec_recomb_cool(T)

   k%HIcic     = Hui_HI_col_ion_cool(T)
   k%HeIcic    = Hui_HeI_col_ion_cool(T)
   k%HeIIcic   = Hui_HeII_col_ion_cool(T)
   k%He23cic   = Cen_He23s_col_ion_cool(T)

   k%HIcec     =  Cen_HI_col_ext_cool(T)
   k%HeIcec    =  Cen_HeI_col_ext_cool(T)
   k%HeIIcec   =  Cen_HeII_col_ext_cool(T)

100 format(25(ES12.5,2X))
    write(10,100) k

end do

close(10)



end program write_rates
