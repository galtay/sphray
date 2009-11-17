!> \file atomic_rates.f90

!> \brief the module that handles atomic rates

module atomic_rates_mod
use myf90_mod
use hui_gnedin_atomic_rates_mod
use cen_atomic_rates_mod
implicit none

  character(14), parameter :: ratenames(25) = &
       (/ "temperature   ",                                     &
          "HIci          ", "HeIci         ", "HeIIci        ", &
          "HIIrcA        ", "HeIIrcA       ", "HeIIIrcA      ", &
          "HIIrcB        ", "HeIIrcB       ", "HeIIIrcB      ", &
          "HeDrc         ",                                     &
          "HIIrccA       ", "HeIIrccA      ", "HeIIIrccA     ", &
          "HIIrccB       ", "HeIIrccB      ", "HeIIIrccB     ", &
          "HeDrcc        ",                                     &
          "HIcic         ", "HeIcic        ", "HeIIcic       ", &
          "He23cic       ",                                     &
          "HIcec         ", "HeIcec        ", "HeIIcec       "  /) !< rate names

  integer(i8b), private, parameter :: RateHeaderLines = 0 !< lines in header
  real(r8b) :: logTfloor   !< log10 of the min temperature in the rates table
  real(r8b) :: logTceiling !< log10 of the max temperature in the rates table
  real(r8b) :: RateLogTempSpacing !< entry spacing in table in log10 T units
  integer(i8b) :: RateFileEntries !< number of entries in the rate file


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




  type(atomic_rates_type), allocatable :: rtable(:)  !< atomic rates
  type(atomic_rates_type) :: iso_k  !< static rates for iso-temperature run
  type(atomic_rates_type) :: cmb_k  !< static rates for cmb-temperature
  type(atomic_rates_type) :: xHII_k !< static rates for xHII-temperature 
  
contains


  !---------------------------------------------
  !> open and read in the atomic rate file
  subroutine read_atomic_rates_file(AtomicRatesFile, OutputDir)
    character(clen) :: AtomicRatesFile
    character(clen) :: OutputDir

    character(clen) :: logfile
    integer(i8b) :: lun,err,i,loglun
    real(r8b) :: Tdum
    logical :: fthere
    
    real(r8b) :: GGHI,RRHII
    real(r8b) :: GGHeI, RRHeII
    real(r8b) :: GGHeII, RRHeIII
    real(r8b) :: xHI,xHII,xHeI,xHeII,xHeIII,den

    character(clen) :: str
    integer, parameter :: verb=1
    character(clen), parameter :: myname = "read_atomic_rates_file"
    logical, parameter :: crash = .true.

    type(atomic_rates_type) :: k

    logfile = trim(OutputDir) // '/' // 'atomic_rates.log'
    call open_formatted_file_w(logfile,loglun)

    write(str,'(A,T27,A)') 'using rates file: ', trim(AtomicRatesFile)
    call mywrite(str, verb) 
    call mywrite('', verb)

    inquire(file=AtomicRatesFile, exist=fthere)
    if (.not. fthere) call myerr('cant find atomic rates file',myname,crash)

    call open_formatted_file_r(AtomicRatesFile,lun)
    do i = 1,RateHeaderLines
       read(lun,*)
    end do

    write(loglun,'(A)') 'reading in the atomic rates header'
    read(lun,*) logTfloor, logTceiling, RateFileEntries, RateLogTempSpacing
    read(lun,*) ! column names

    write(loglun,'(A,ES10.4)') "log10 of temperature floor   = ", logTfloor
    write(loglun,'(A,ES10.4)') "log10 of temperature ceiling = ", logTceiling
    write(loglun,'(A,I5)')     "rate file entries            = ", RateFileEntries
    write(loglun,'(A,ES10.4)') "log10 spacing of temperature = ", RateLogTempSpacing
    write(loglun,*)    

    allocate(rtable(RateFileEntries), stat=err)
    if(err .ne. 0) call myerr("cant allocate rates table",myname,crash)
    do i = 1,RateFileEntries
       read(lun,*) rtable(i)
    end do

    ! write out characteristic rates @ T = 10,000 K
    Tdum = 1.0E4
    call get_atomic_rates(Tdum,k)

    100 format(A,ES12.5,A)
    101 format(A,":",T32,ES12.5,A)
    102 format(A,T15,ES12.3)

    write(loglun,*) 
    write(loglun,*) "-----------------------------------------------------------"
    write(loglun,100) "Atomic Rates @ ", Tdum," K from rate tables"
    write(loglun,*) 
    write(loglun,101) "HI Col Ion", k%HIci,   " cm^3/s"

    write(loglun,101) "HeI Col Ion", k%HeIci,  " cm^3/s"
    write(loglun,101) "HeII Col Ion", k%HeIIci, " cm^3/s"

    write(loglun,*) 
    write(loglun,101) "HII Rec (A)", k%HIIrcA,   " cm^3/s"
    write(loglun,101) "HII Rec (B)", k%HIIrcB,   " cm^3/s"
    write(loglun,*) 
    write(loglun,101) "HeII Rec (A)", k%HeIIrcA,  " cm^3/s"
    write(loglun,101) "HeII Rec (B)", k%HeIIrcB,  " cm^3/s"
    write(loglun,*) 
    write(loglun,101) "HeIII Rec (A)", k%HeIIIrcA, " cm^3/s"
    write(loglun,101) "HeIII Rec (B)", k%HeIIIrcB, " cm^3/s"
    write(loglun,*) 
    write(loglun,101) "He D Rec", k%HeDrc, "cm^3/s"
    write(loglun,*) 

    write(loglun,101) "HII Rec Cool (A)", k%HIIrccA, " ergs cm^3/s"
    write(loglun,101) "HII Rec Cool (B)", k%HIIrccB, " ergs cm^3/s"
    write(loglun,*) 
    write(loglun,101) "HeII Rec Cool (A)", k%HeIIrccA,  " ergs cm^3/s"
    write(loglun,101) "HeII Rec Cool (B)", k%HeIIrccB,  " ergs cm^3/s"
    write(loglun,*) 
    write(loglun,101) "HeIII Rec Cool (A)", k%HeIIIrccA, " ergs cm^3/s"
    write(loglun,101) "HeIII Rec Cool (B)", k%HeIIIrccB, " ergs cm^3/s"
    write(loglun,*) 
    write(loglun,101) "He D Rec Cool", k%HeDrc, "cm^3/s"
    write(loglun,*) 

    write(loglun,101) "HI Col Ion Cool", k%HIcic, " ergs cm^3/s"
    write(loglun,101) "HeI Col Ion Cool", k%HeIcic, " ergs cm^3/s"
    write(loglun,101) "HeII Col Ion Cool", k%HeIIcic, " ergs cm^3/s"
    write(loglun,101) "He 23 Col Ion Cool", k%He23cic, " ergs cm^3/s"
    write(loglun,*) 


    write(loglun,101) "HI Col Ext Cool", k%HIcec, " ergs cm^3/s"
    write(loglun,101) "HeI Col Ext Cool", k%HeIcec, " ergs cm^6/s"
    write(loglun,101) "HeII Col Ext Cool", k%HeIIcec, " ergs cm^3/s"

    write(loglun,*) 

      
    write(loglun,'(A)') "Collisional Equilibrium"
    write(loglun,*) 
      
    GGHI   = k%HIci 
    GGHeI  = k%HeIci
    GGHeII = k%HeIIci
    
    RRHII   = k%HIIrcA
    RRHeII  = k%HeIIrcA
    RRHeIII = k%HeIIIrcA
    
    den = (RRHII + GGHI)
    xHI = RRHII / den
    xHII = GGHI / den
    
    den = (RRHeII * RRHeIII + RRHeIII * GGHeI + GGHeI * GGHeII)
    xHeI   = RRHeII * RRHeIII / den
    xHeII  = RRHeIII * GGHeI  / den 
    xHeIII = GGHeI * GGHeII   / den
    
    write(loglun,*)   " CASE A"
    write(loglun,102) "xHIeq    = ", xHI
    write(loglun,102) "xHIIeq   = ", xHII
    write(loglun,102) "xHeIeq   = ", xHeI
    write(loglun,102) "xHeIIeq  = ", xHeII
    write(loglun,102) "xHeIIIeq = ", xHeIII
    write(loglun,*) 
    
    RRHII   = k%HIIrcB
    RRHeII  = k%HeIIrcB
    RRHeIII = k%HeIIIrcB
    
    den = (RRHII + GGHI)
    xHI = RRHII / den
    xHII = GGHI / den
    
    den = (RRHeII * RRHeIII + RRHeIII * GGHeI + GGHeI * GGHeII)
    xHeI   = RRHeII * RRHeIII / den
    xHeII  = RRHeIII * GGHeI  / den 
    xHeIII = GGHeI * GGHeII   / den
    
    write(loglun,*)   " CASE B"
    write(loglun,102) "xHIeq    = ", xHI
    write(loglun,102) "xHIIeq   = ", xHII
    write(loglun,102) "xHeIeq   = ", xHeI
    write(loglun,102) "xHeIIeq  = ", xHeII
    write(loglun,102) "xHeIIIeq = ", xHeIII
    write(loglun,*) 
    
  
  end subroutine read_atomic_rates_file


  !-------------------------------------------------------------------------------
  !> Calculates collisional ionization equilibrium for all species 
  !! at given temperature using the rates table
  subroutine calc_colion_eq_table(Tin,caseA,xvec)
    real(r8b), intent(in) :: Tin
    logical, intent(in) :: caseA(2)
    real(r4b), intent(out) :: xvec(5)
    
    type(atomic_rates_type) :: k
    real(r8b) :: GGHI, GGHeI, GGHeII
    real(r8b) :: RRHII, RRHeII, RRHeIII
    real(r8b) :: den
    
    call get_atomic_rates(Tin,k)
    
    GGHI   = k%HIci 
    GGHeI  = k%HeIci
    GGHeII = k%HeIIci
    
    if (caseA(1)) then
       RRHII   = k%HIIrcA 
    else
       RRHII   = k%HIIrcB 
    end if
    
    if (caseA(2)) then
       RRHeII  = k%HeIIrcA 
       RRHeIII = k%HeIIIrcA 
    else
       RRHeII  = k%HeIIrcB 
       RRHeIII = k%HeIIIrcB 
    end if
    
    den = (RRHII + GGHI)
    xvec(1) = RRHII / den
    xvec(2) = GGHI / den
    
    den = (RRHeII * RRHeIII + RRHeIII * GGHeI + GGHeI * GGHeII)
    xvec(3) = RRHeII * RRHeIII / den
    xvec(4) = RRHeIII * GGHeI  / den 
    xvec(5) = GGHeI * GGHeII   / den
    
  end subroutine calc_colion_eq_table
  

  !-------------------------------------------------------------------------------
  !> Calculates collisional ionization equilibrium for all species 
  !! at given temperature using the fits directly
  subroutine calc_colion_eq_fits(Tin,caseA,xvec)
    real(r8b), intent(in) :: Tin
    logical, intent(in) :: caseA(2)
    real(r8b), intent(out) :: xvec(5)
    
    real(r8b) :: Tfloor, Tceil, T
    real(r8b) :: GGHI, GGHeI, GGHeII
    real(r8b) :: RRHII, RRHeII, RRHeIII
    real(r8b) :: den
    
    Tfloor = 10**(logTfloor)
    Tceil = 10**(logTceiling)
    T = Tin
    
    if (T > Tceil ) T = Tceil
    if (T < Tfloor) T = Tfloor
    
    GGHI   = Hui_HI_col_ion(T)
    GGHeI  = Hui_HeI_col_ion(T)
    GGHeII = Hui_HeII_col_ion(T)
    
    if (caseA(1)) then
       RRHII   = Hui_HII_recombA(T)
    else
       RRHII   = Hui_HII_recombB(T)
    end if
    
    if (caseA(2)) then
       RRHeII  = Hui_HeII_recombA(T)
       RRHeIII = Hui_HeIII_recombA(T)
    else
       RRHeII  = Hui_HeII_recombB(T)
       RRHeIII = Hui_HeIII_recombB(T)
    end if
    
    den = (RRHII + GGHI)
    xvec(1) = RRHII / den
    xvec(2) = GGHI / den
    
    den = (RRHeII * RRHeIII + RRHeIII * GGHeI + GGHeI * GGHeII)
    xvec(3) = RRHeII * RRHeIII / den
    xvec(4) = RRHeIII * GGHeI  / den 
    xvec(5) = GGHeI * GGHeII   / den
    
  end subroutine calc_colion_eq_fits



  !--------------------------------------------------------------------
  !> takes a temperature and calculates the atomic rate table index and
  !! how far the true temperature is past the table temperature
  subroutine get_Tindx_and_Tremain(T,Tindx,Tremain)
  
    real(r8b), intent(in) :: T         !< input temperature
    integer(i8b), intent(out) :: Tindx !< atomic rate table index
    real(r8b), intent(out) :: Tremain  !< rates%(T) + Tremain = T
    real(r8b) :: Tnum
    
    if (T < 0.0) stop "T < 0.0 in get_Tindx_and_Tremain"
    
    Tnum = ( log10(T) - logTfloor ) / RateLogTempSpacing
    Tindx = ceiling(Tnum)
    
    if (Tindx == 0) then
       Tindx = 1
       Tremain = 0.0d0
    end if
    Tremain = Tnum - (Tindx-1)
    
    if(Tindx < 1)                 stop "Tindx < 1 "
    if(Tindx > RateFileEntries-1) stop "Tindx > RateFileEntries - 1 "
    
  end subroutine get_Tindx_and_Tremain



  !> sets static rates for isotemperature run
  !================================================================
  subroutine set_iso_atomic_rates(Tiso)
    real(r8b), intent(in) :: Tiso !< iso temperature 

    call get_atomic_rates(Tiso, iso_k)

  end subroutine set_iso_atomic_rates
  
  
  !> sets cmb temperature atomic rates 
  !================================================================
  subroutine set_cmb_atomic_rates(Tcmb)
    real(r8b), intent(in) :: Tcmb !< cmb temperature 

    call get_atomic_rates(Tcmb, cmb_k)

  end subroutine set_cmb_atomic_rates

  !> sets ionized region atomic rates
  !================================================================
  subroutine set_xHII_atomic_rates(TxHII)
    real(r8b), intent(in) :: TxHII !< xHII temperature 

    call get_atomic_rates(TxHII, xHII_k)

  end subroutine set_xHII_atomic_rates


  !=================================
  !> get atomic rates from table 
  subroutine get_atomic_rates(T,k)

    real(r8b), intent(in) :: T                          !< input temperature
    type(atomic_rates_type), intent(out) :: k  !< returns necessary rates
  
    integer(i8b) :: Tindx
    real(r8b) :: Tremain, rdif

    call get_Tindx_and_Tremain(T,Tindx,Tremain)


    ! Collisional Ionization

    rdif = (rtable(Tindx+1)%HIci - rtable(Tindx)%HIci) * Tremain
    k%HIci = (rtable(Tindx)%HIci + rdif) 
    
    rdif = (rtable(Tindx+1)%HeIci - rtable(Tindx)%HeIci) * Tremain
    k%HeIci = (rtable(Tindx)%HeIci + rdif) 
    
    rdif = (rtable(Tindx+1)%HeIIci - rtable(Tindx)%HeIIci) * Tremain
    k%HeIIci = (rtable(Tindx)%HeIIci + rdif) 
    
    ! Recombination
    
    rdif = (rtable(Tindx+1)%HIIrcA - rtable(Tindx)%HIIrcA) * Tremain
    k%HIIrcA = (rtable(Tindx)%HIIrcA + rdif) 
    
    rdif = (rtable(Tindx+1)%HIIrcB - rtable(Tindx)%HIIrcB) * Tremain
    k%HIIrcB = (rtable(Tindx)%HIIrcB + rdif) 
    
    rdif = (rtable(Tindx+1)%HeIIrcA - rtable(Tindx)%HeIIrcA) * Tremain
    k%HeIIrcA = (rtable(Tindx)%HeIIrcA + rdif) 
    
    rdif = (rtable(Tindx+1)%HeIIrcB - rtable(Tindx)%HeIIrcB) * Tremain
    k%HeIIrcB = (rtable(Tindx)%HeIIrcB + rdif) 
    
    rdif = (rtable(Tindx+1)%HeIIIrcA - rtable(Tindx)%HeIIIrcA) * Tremain
    k%HeIIIrcA = (rtable(Tindx)%HeIIIrcA + rdif) 
    
    rdif = (rtable(Tindx+1)%HeIIIrcB - rtable(Tindx)%HeIIIrcB) * Tremain
    k%HeIIIrcB = (rtable(Tindx)%HeIIIrcB + rdif) 
    
    rdif = (rtable(Tindx+1)%HeDrc - rtable(Tindx)%HeDrc) * Tremain
    k%HeDrc = (rtable(Tindx)%HeDrc + rdif) 
    
    ! Recombination Cooling
    
    rdif = (rtable(Tindx+1)%HIIrccA - rtable(Tindx)%HIIrccA) * Tremain
    k%HIIrccA = (rtable(Tindx)%HIIrccA + rdif) 
    
    rdif = (rtable(Tindx+1)%HIIrccB - rtable(Tindx)%HIIrccB) * Tremain
    k%HIIrccB = (rtable(Tindx)%HIIrccB + rdif) 
    
    rdif = (rtable(Tindx+1)%HeIIrccA - rtable(Tindx)%HeIIrccA) * Tremain
    k%HeIIrccA = (rtable(Tindx)%HeIIrccA + rdif) 
    
    rdif = (rtable(Tindx+1)%HeIIrccB - rtable(Tindx)%HeIIrccB) * Tremain
    k%HeIIrccB = (rtable(Tindx)%HeIIrccB + rdif) 
    
    rdif = (rtable(Tindx+1)%HeIIIrccA - rtable(Tindx)%HeIIIrccA) * Tremain
    k%HeIIIrccA = (rtable(Tindx)%HeIIIrccA + rdif) 
    
    rdif = (rtable(Tindx+1)%HeIIIrccB - rtable(Tindx)%HeIIIrccB) * Tremain
    k%HeIIIrccB = (rtable(Tindx)%HeIIIrccB + rdif) 
    
    rdif = (rtable(Tindx+1)%HeDrcc - rtable(Tindx)%HeDrcc) * Tremain
    k%HeDrcc = (rtable(Tindx)%HeDrcc + rdif) 
    
    ! Collisional Ionization Cooling
    
    rdif = (rtable(Tindx+1)%HIcic - rtable(Tindx)%HIcic) * Tremain
    k%HIcic = (rtable(Tindx)%HIcic + rdif) 
    
    rdif = (rtable(Tindx+1)%HeIcic - rtable(Tindx)%HeIcic) * Tremain
    k%HeIcic = (rtable(Tindx)%HeIcic + rdif) 
    
    rdif = (rtable(Tindx+1)%HeIIcic - rtable(Tindx)%HeIIcic) * Tremain
    k%HeIIcic = (rtable(Tindx)%HeIIcic + rdif) 
    
    rdif = (rtable(Tindx+1)%He23cic - rtable(Tindx)%He23cic) * Tremain
    k%He23cic = (rtable(Tindx)%He23cic + rdif) 
    
    ! Collisional Excitation Cooling
    
    rdif = (rtable(Tindx+1)%HIcec - rtable(Tindx)%HIcec) * Tremain
    k%HIcec = (rtable(Tindx)%HIcec + rdif) 
    
    rdif = (rtable(Tindx+1)%HeIcec - rtable(Tindx)%HeIcec) * Tremain
    k%HeIcec = (rtable(Tindx)%HeIcec + rdif) 
    
    rdif = (rtable(Tindx+1)%HeIIcec - rtable(Tindx)%HeIIcec) * Tremain
    k%HeIIcec = (rtable(Tindx)%HeIIcec + rdif) 
    
    
  end subroutine get_atomic_rates



end module atomic_rates_mod
