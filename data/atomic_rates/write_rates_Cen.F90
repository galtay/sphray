module atomic_rates
use myf90_mod
implicit none

character(14), parameter :: names(25) = (/ "Log T" , &
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

!> single atomic rate 
!---------------------------
type atomic_rate_type
   character(14) :: source
   real(r4b), allocatable :: rate(:)   !< allocate to number of temperature bins
end type atomic_rate_type



!> atomic rates table
!---------------------------
type atomic_rates_table

   real(r4b) :: logT1
   real(r4b) :: logT2
   integer   :: Tbins
   real(r4b) :: dTlog

   real(r4b), allocatable :: logT(:)   !< log10 temperature

   type(atomic_rate_type) :: HIci      !< HI   collisional ionization rate 
   type(atomic_rate_type) :: HeIci     !< HeI  collisional ionization rate 
   type(atomic_rate_type) :: HeIIci    !< HeII collisional ionization rate 

   type(atomic_rate_type) :: HIIrcA    !< HII   recombination rate (case A)
   type(atomic_rate_type) :: HeIIrcA   !< HeII  recombination rate (case A)
   type(atomic_rate_type) :: HeIIIrcA  !< HeIII recombination rate (case A)

   type(atomic_rate_type) :: HIIrcB    !< HII   recombination rate (case B)
   type(atomic_rate_type) :: HeIIrcB   !< HeII  recombination rate (case B)
   type(atomic_rate_type) :: HeIIIrcB  !< HeIII recombination rate (case B)
   type(atomic_rate_type) :: HeDrc     !< dielectronic He recombination rate   

   type(atomic_rate_type) :: HIIrccA   !< HII   recombination cooling rate (case A)
   type(atomic_rate_type) :: HeIIrccA  !< HeII  recombination cooling rate (case A)
   type(atomic_rate_type) :: HeIIIrccA !< HeIII recombination cooling rate (case A)

   type(atomic_rate_type) :: HIIrccB   !< HII   recombination cooling rate (case B)
   type(atomic_rate_type) :: HeIIrccB  !< HeII  recombination cooling rate (case B)
   type(atomic_rate_type) :: HeIIIrccB !< HeIII recombination cooling rate (case B)
   type(atomic_rate_type) :: HeDrcc    !< dielectronic He recombination cooling rate
   
   type(atomic_rate_type) :: HIcic     !< HI   collisional ionization cooling rate
   type(atomic_rate_type) :: HeIcic    !< HeI  collisional ionization cooling rate
   type(atomic_rate_type) :: HeIIcic   !< HeII collisional ionization cooling rate
   type(atomic_rate_type) :: He23cic   !< He23 collisional ionization cooling rate
   
   type(atomic_rate_type) :: HIcec     !< HI   collisional excitation cooling rate
   type(atomic_rate_type) :: HeIcec    !< HeI  collisional excitation cooling rate
   type(atomic_rate_type) :: HeIIcec   !< HeII collisional excitation cooling rate
      
end type atomic_rates_table


contains

  subroutine allocate_atomic_rates_table( table, n_bins )
    type(atomic_rates_table) :: table
    integer :: n_bins

    allocate( table%logT(n_bins) )

    allocate( table%HIci%rate(n_bins) )
    allocate( table%HeIci%rate(n_bins) )
    allocate( table%HeIIci%rate(n_bins) )
    
    allocate( table%HIIrcA%rate(n_bins) )
    allocate( table%HeIIrcA%rate(n_bins) )
    allocate( table%HeIIIrcA%rate(n_bins) )
    
    allocate( table%HIIrcB%rate(n_bins) )
    allocate( table%HeIIrcB%rate(n_bins) )
    allocate( table%HeIIIrcB%rate(n_bins) )
    allocate( table%HeDrc%rate(n_bins) )
    
    allocate( table%HIIrccA%rate(n_bins) )
    allocate( table%HeIIrccA%rate(n_bins) )
    allocate( table%HeIIIrccA%rate(n_bins) )
    
    allocate( table%HIIrccB%rate(n_bins) )
    allocate( table%HeIIrccB%rate(n_bins) )
    allocate( table%HeIIIrccB%rate(n_bins) )
    allocate( table%HeDrcc%rate(n_bins) )
    
    allocate( table%HIcic%rate(n_bins) )
    allocate( table%HeIcic%rate(n_bins) )
    allocate( table%HeIIcic%rate(n_bins) )
    allocate( table%He23cic%rate(n_bins) )
    
    allocate( table%HIcec%rate(n_bins) )
    allocate( table%HeIcec%rate(n_bins) )
    allocate( table%HeIIcec%rate(n_bins) )
    
  end subroutine allocate_atomic_rates_table






subroutine write_atomic_rates_table_to_file( table, file )
  type(atomic_rates_table) :: table
  character(*) :: file

  character(clen) :: fmt
  integer :: j

  
  open(unit=10,file=file)

  fmt = "(2ES12.5,I8,ES12.5)"
  write(10,fmt) table%logT1, table%logT2, table%Tbins, table%dTlog

  fmt = "(1X,25(A))"
  write(10,fmt) names

  write(10,fmt) &
       "              ", & 

       table%HIci%source, &       
       table%HeIci%source, &      
       table%HeIIci%source, &     
       
       table%HIIrcA%source, &     
       table%HeIIrcA%source, &    
       table%HeIIIrcA%source, &   
       
       table%HIIrcB%source, &     
       table%HeIIrcB%source, &    
       table%HeIIIrcB%source, &     
       table%HeDrc%source, &      
       
       table%HIIrccA%source, &    
       table%HeIIrccA%source, &    
       table%HeIIIrccA%source, &  
       
       table%HIIrccB%source, &    
       table%HeIIrccB%source, &   
       table%HeIIIrccB%source, &      
       table%HeDrcc%source, & 
       
       table%HIcic%source, &      
       table%HeIcic%source, &     
       table%HeIIcic%source, &    
       table%He23cic%source, &    
       
       table%HIcec%source, &      
       table%HeIcec%source, &     
       table%HeIIcec%source   


  fmt = "(25(ES12.5,2X))"

  do j = 1, table%Tbins
     write(10,fmt) &

          table%logT(j), & 
          
          table%HIci%rate(j), &       
          table%HeIci%rate(j), &      
          table%HeIIci%rate(j), &     
          
          table%HIIrcA%rate(j), &     
          table%HeIIrcA%rate(j), &    
          table%HeIIIrcA%rate(j), &   
          
          table%HIIrcB%rate(j), &     
          table%HeIIrcB%rate(j), &    
          table%HeIIIrcB%rate(j), &     
          table%HeDrc%rate(j), &      
          
          table%HIIrccA%rate(j), &    
          table%HeIIrccA%rate(j), &    
          table%HeIIIrccA%rate(j), &  
          
          table%HIIrccB%rate(j), &    
          table%HeIIrccB%rate(j), &   
          table%HeIIIrccB%rate(j), &      
          table%HeDrcc%rate(j), & 
          
          table%HIcic%rate(j), &      
          table%HeIcic%rate(j), &     
          table%HeIIcic%rate(j), &    
          table%He23cic%rate(j), &    
          
          table%HIcec%rate(j), &      
          table%HeIcec%rate(j), &     
          table%HeIIcec%rate(j)   
     
  end do
     
end subroutine write_atomic_rates_table_to_file







  
end module atomic_rates



program write_rates_Cen
use myf90_mod
use hui_gnedin_atomic_rates_mod
use cen_atomic_rates_mod
use atomic_rates
implicit none

real(r8b), parameter :: T1 = 1.0d0
real(r8b), parameter :: T2 = 1.0d9
real(r8b), parameter :: logT1 = log10(T1)
real(r8b), parameter :: logT2 = log10(T2)
integer(i8b), parameter :: Tbins = 2251
real(r8b), parameter :: dTlog = (logT2 - logT1) / real(Tbins-1)
integer(i8b) :: i, j
real(r8b) :: T
real(r8b) :: logT

character(200) :: ratesfile




type(atomic_rates_table) :: k

write(*,*) "dTlog = ", dTlog
write(*,*) "logT2 - logT2 = ", logT2 - logT1

k%logT1 = logT1
k%logT2 = logT2
k%Tbins = Tbins
k%dTlog = dTlog

call allocate_atomic_rates_table( k, k%Tbins )


k%HIci%source      = "cen92" 
k%HeIci%source     = "cen92" 
k%HeIIci%source    = "cen92" 

k%HIIrcA%source    = "cen92" 
k%HeIIrcA%source   = "cen92" 
k%HeIIIrcA%source  = "cen92" 

k%HIIrcB%source    = "hui97"
k%HeIIrcB%source   = "hui97"
k%HeIIIrcB%source  = "hui97" 
k%HeDrc%source     = "cen92" 

k%HIIrccA%source   = "cen92" 
k%HeIIrccA%source  = "cen92" 
k%HeIIIrccA%source = "cen92" 

k%HIIrccB%source   = "hui97"
k%HeIIrccB%source  = "hui97"
k%HeIIIrccB%source = "hui97"   
k%HeDrcc%source    = "cen92" 

k%HIcic%source     = "cen92" 
k%HeIcic%source    = "cen92" 
k%HeIIcic%source   = "cen92" 
k%He23cic%source   = "cen92" 

k%HIcec%source     = "cen92" 
k%HeIcec%source    = "cen92" 
k%HeIIcec%source   = "cen92" 



do i = 0,Tbins-1

   logT = logT1 + i*dTlog
   T = 10**logT

   j = i+1
   k%logT(j) = logT

   k%HIci%rate(j)      = Cen_HI_col_ion(T)   
   k%HeIci%rate(j)     = Cen_HeI_col_ion(T)  
   k%HeIIci%rate(j)    = Cen_HeII_col_ion(T) 

   k%HIIrcA%rate(j)    = Cen_HII_recombA(T)
   k%HeIIrcA%rate(j)   = Cen_HeII_recombA(T) 
   k%HeIIIrcA%rate(j)  = Cen_HeIII_recombA(T) 

   k%HIIrcB%rate(j)    = Hui_HII_recombB(T)
   k%HeIIrcB%rate(j)   = Hui_HeII_recombB(T) 
   k%HeIIIrcB%rate(j)  = Hui_HeIII_recombB(T)      
   k%HeDrc%rate(j)     = Cen_He_dielec_recomb(T)

   k%HIIrccA%rate(j)   = Cen_HII_recomb_coolA(T)
   k%HeIIrccA%rate(j)  = Cen_HeII_recomb_coolA(T) 
   k%HeIIIrccA%rate(j) = Cen_HeIII_recomb_coolA(T) 

   k%HIIrccB%rate(j)   = Hui_HII_rec_coolB(T)
   k%HeIIrccB%rate(j)  = Hui_HeII_rec_coolB(T) 
   k%HeIIIrccB%rate(j) = Hui_HeIII_rec_coolB(T)      
   k%HeDrcc%rate(j)    = Cen_He_dielec_recomb_cool(T)

   k%HIcic%rate(j)     = Cen_HI_col_ion_cool(T)
   k%HeIcic%rate(j)    = Cen_HeI_col_ion_cool(T)
   k%HeIIcic%rate(j)   = Cen_HeII_col_ion_cool(T)
   k%He23cic%rate(j)   = Cen_He23s_col_ion_cool(T)

   k%HIcec%rate(j)     = Cen_HI_col_ext_cool(T)
   k%HeIcec%rate(j)    = Cen_HeI_col_ext_cool(T)
   k%HeIIcec%rate(j)   = Cen_HeII_col_ext_cool(T)

end do

call write_atomic_rates_table_to_file( k, "atomic_rates_Cen.txt" )

close(10)



end program write_rates_Cen
