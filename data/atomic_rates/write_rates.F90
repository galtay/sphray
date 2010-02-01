program write_rates
use myf90_mod
use hui_gnedin_atomic_rates_mod
use cen_atomic_rates_mod
use atomic_rates_mod
implicit none

real(r8b), parameter :: T1 = 1.0d0
real(r8b), parameter :: T2 = 1.0d9
real(r8b), parameter :: logT1 = log10(T1)
real(r8b), parameter :: logT2 = log10(T2)
integer(i8b), parameter :: Tbins = 2251
real(r8b), parameter :: dlogT = (logT2 - logT1) / real(Tbins-1)
integer(i8b) :: i, j
real(r8b) :: T
real(r8b) :: logT

character(200) :: ratesfile


type(atomic_rates_table_type) :: table

write(*,*) "dlogT = ", dlogT
write(*,*) "logT2 - logT2 = ", logT2 - logT1

table%logT1 = logT1
table%logT2 = logT2
table%Tbins = Tbins
table%dlogT = dlogT

call allocate_atomic_rates_table( table, table%Tbins )

table%HIci%source      = "hui97" 
table%HeIci%source     = "hui97" 
table%HeIIci%source    = "hui97" 

table%HIIrcA%source    = "hui97" 
table%HeIIrcA%source   = "hui97" 
table%HeIIIrcA%source  = "hui97" 

table%HIIrcB%source    = "hui97"
table%HeIIrcB%source   = "hui97"
table%HeIIIrcB%source  = "hui97" 
table%HeDrc%source     = "hui97" 

table%HIIrccA%source   = "hui97" 
table%HeIIrccA%source  = "hui97" 
table%HeIIIrccA%source = "hui97" 

table%HIIrccB%source   = "hui97"
table%HeIIrccB%source  = "hui97"
table%HeIIIrccB%source = "hui97"   
table%HeDrcc%source    = "hui97" 

table%HIcic%source     = "hui97" 
table%HeIcic%source    = "hui97" 
table%HeIIcic%source   = "hui97" 
table%He23cic%source   = "cen92" 

table%HIcec%source     = "cen92" 
table%HeIcec%source    = "cen92" 
table%HeIIcec%source   = "cen92" 





do i = 0,Tbins-1

   logT = logT1 + i*dlogT
   T = 10**logT

   j = i+1
   table%logT(j) = logT


   table%HIci%rate(j)   = Hui_HI_col_ion(T)   
   table%HeIci%rate(j)  = Hui_HeI_col_ion(T)  
   table%HeIIci%rate(j) = Hui_HeII_col_ion(T) 

   table%HIIrcA%rate(j)   = Hui_HII_recombA(T)
   table%HeIIrcA%rate(j)  = Hui_HeII_recombA(T) 
   table%HeIIIrcA%rate(j) = Hui_HeIII_recombA(T) 

   table%HIIrcB%rate(j)   = Hui_HII_recombB(T)
   table%HeIIrcB%rate(j)  = Hui_HeII_recombB(T) 
   table%HeIIIrcB%rate(j) = Hui_HeIII_recombB(T)      
   table%HeDrc%rate(j)    = Hui_He_dielec_recomb(T)

   table%HIIrccA%rate(j)   = Hui_HII_rec_coolA(T)
   table%HeIIrccA%rate(j)  = Hui_HeII_rec_coolA(T) 
   table%HeIIIrccA%rate(j) = Hui_HeIII_rec_coolA(T) 

   table%HIIrccB%rate(j)   = Hui_HII_rec_coolB(T)
   table%HeIIrccB%rate(j)  = Hui_HeII_rec_coolB(T) 
   table%HeIIIrccB%rate(j) = Hui_HeIII_rec_coolB(T)      
   table%HeDrcc%rate(j)    = Hui_He_dielec_recomb_cool(T)

   table%HIcic%rate(j)     = Hui_HI_col_ion_cool(T)
   table%HeIcic%rate(j)    = Hui_HeI_col_ion_cool(T)
   table%HeIIcic%rate(j)   = Hui_HeII_col_ion_cool(T)
   table%He23cic%rate(j)   = Cen_He23s_col_ion_cool(T)

   table%HIcec%rate(j)     =  Cen_HI_col_ext_cool(T)
   table%HeIcec%rate(j)    =  Cen_HeI_col_ext_cool(T)
   table%HeIIcec%rate(j)   =  Cen_HeII_col_ext_cool(T)


end do


call write_atomic_rates_table_to_file( table, "atomic_rates.txt" )




end program write_rates
