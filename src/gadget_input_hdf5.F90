!> \file gadget_input_hdf5.F90

!> \brief Handles readin of GADGET HDF5 formatted files
!<

module gadget_input_hdf5_mod
use myf90_mod
use gadget_header_class
use ion_table_class
use particle_system_mod, only: set_ye
use global_mod, only: psys, PLAN, GV
use global_mod, only: saved_gheads, gconst


#ifdef hdf5
use hdf5_wrapper
#endif
implicit none
private


public :: get_planning_data_gadget_hdf5
public :: read_Ghdf5_particles
public :: gadget_output_hdf5

contains

#ifndef hdf5

! these are dummy subroutines so that the calls outside this file
! dont have to be wrapped with pre processor macros.

subroutine read_gadget_header_hdf5()
  call myerr("this routine shuold not have been called","hdf5dummy",crash=.true.)
end subroutine read_gadget_header_hdf5

subroutine get_planning_data_gadget_hdf5()
  call myerr("this routine shuold not have been called","hdf5dummy",crash=.true.)
end subroutine get_planning_data_gadget_hdf5

subroutine read_Ghdf5_particles()
  call myerr("this routine shuold not have been called","hdf5dummy",crash=.true.)
end subroutine read_Ghdf5_particles

subroutine gadget_output_hdf5()
  call myerr("this routine shuold not have been called","hdf5dummy",crash=.true.)
end subroutine gadget_output_hdf5


#else



!>   gets run planning data from Gadget Headers
!===============================================
subroutine get_planning_data_gadget_hdf5()

  type(gadget_header_type) :: ghead
  type(gadget_units_type) :: gunits

  integer(i4b) :: fh      !< hdf5 file handle
  integer(i8b) :: iSnap   !< first snap number
  integer(i8b) :: fSnap   !< last snap number
  integer(i8b) :: pfiles  !< nfiles for particle snapshots
  integer(i8b) :: i,j     !< counters

  character(clen) :: snapfile

  real(r8b) :: Time_GYR
  integer(i8b) :: lun
  integer(i8b) :: loglun
  character(clen) :: logfile
  real(r8b) :: kpc2cm
  real(r8b) :: km2cm
 
  ! open up the planning data log file
  !======================================================
  logfile = trim(GV%OutputDir) // "/" // "particle_headers.log"
  call open_formatted_file_w(logfile,loglun)

  
  ! these global variables are read from the config file
  !======================================================
  iSnap  = GV%StartSnapNum
  fSnap  = GV%EndSnapNum    
  pfiles = GV%ParFilesPerSnap

  if ( allocated(saved_gheads) ) deallocate(saved_gheads)
  allocate( saved_gheads(iSnap:fSnap, 0:pfiles-1) )


  ! read all particle headers and write to log file
  !===================================================
  write(loglun,'(A)') "reading all GADGET particle header(s) ... "
  do i = iSnap,fSnap
     do j = 0,pfiles-1


        call form_gadget_snapshot_file_name(GV%SnapPath,GV%ParFileBase,i,j,snapfile)
        write(loglun,'(I3,"  ",A)') i,trim(snapfile)
        write(*,*) 'snapfile = ', trim(snapfile) 

        call read_gadget_header_file(snapfile, ghead)
        call gadget_header_to_file(ghead,loglun)
        saved_gheads(i,j) = ghead

        ! make sure there is gas in this snapshot
        if (.not. ghead%npar_all(1) > 0) then
           write(*,*) "Gadget snapshot does not contain any gas particles"
           write(*,*) "Sphray cannot read dark matter particles directly, "
           write(*,*) "please calculate smoothing lengths for these particles"
           write(*,*) "and write them as gas particles. "
           stop
        end if
        
        ! IsoTemp comes from config file

        GV%BoxLwrsComoh(:) = 0.0d0
        GV%BoxUprsComoh(:) = ghead%boxlen
                
        GV%OmegaM = ghead%OmegaM
        GV%OmegaL = ghead%OmegaL
        GV%OmegaB = ghead%OmegaB
        GV%LittleH = ghead%h

        call read_gadget_units( snapfile, gunits )
        call read_gadget_constants( snapfile, gconst )

        GV%cgs_len  = gunits%len
        GV%cgs_mass = gunits%mass
        GV%cgs_vel  = gunits%vel

        GV%cgs_time = gunits%time
        GV%cgs_rho  = gunits%rho
        GV%cgs_prs  = gunits%prs
        GV%cgs_enrg = gunits%energy


        call hdf5_open_file( fh, snapfile, readonly=.true. )
        call hdf5_read_attribute(fh,'Header/Time_GYR', Time_GYR)
        call hdf5_close_file( fh )

                
        if (GV%Comoving) then
           PLAN%snap(i)%ScalefacAt = ghead%a
           PLAN%snap(i)%TimeAt = Time_GYR * gconst%sec_per_megayear * 1.0d3      ! in seconds
           PLAN%snap(i)%TimeAt = PLAN%snap(i)%TimeAt * GV%LittleH / GV%cgs_time  ! in code units
        else
           PLAN%snap(i)%ScalefacAt = 1.0d0 / (1.0d0 + ghead%z)
           PLAN%snap(i)%TimeAt = Time_GYR * gconst%sec_per_megayear * 1.0d3      ! in seconds
           PLAN%snap(i)%TimeAt = PLAN%snap(i)%TimeAt * GV%LittleH / GV%cgs_time  ! in code units
        end if
        
     end do
  end do



  ! close headers log file
  !========================
  close(loglun)





  ! write units to log file
  !===================================================

  logfile = trim(GV%OutputDir) // "/" // "code_units.log"
  call open_formatted_file_w(logfile,loglun)

  kpc2cm = gconst%cm_per_mpc * 1.0d-3
  km2cm = 1.0d5

  106    format(A,ES12.5,A)
  write(loglun,*) 
  write(loglun,'(A)') "setting code units ..."
  write(loglun,106) "  Hubble:     ", GV%LittleH, " = H0[km/s/Mpc] / 100 "
  write(loglun,106) "  length:     ", GV%cgs_len / kpc2cm, " [kpc/h]"
  write(loglun,106) "  mass:       ", GV%cgs_mass / gconst%solar_mass, " [Msun/h]"
  write(loglun,106) "  velocity    ", GV%cgs_vel / km2cm, " [km/s]"
  write(loglun,106) "  time:       ", GV%cgs_time / gconst%sec_per_megayear, " [Myr/h]"
  write(loglun,106) "  density:    ", GV%cgs_rho / (gconst%solar_mass / kpc2cm**3), " [Msun/kpc^3 h^2]"
  write(loglun,106) "  pressure:   ", GV%cgs_prs, " [dyne/cm^2 h^2]"
  write(loglun,106) "  energy:     ", GV%cgs_enrg, " [ergs/h]"
  write(loglun,*) 
  
  close(loglun)

end subroutine get_planning_data_gadget_hdf5







!> reads a Gadget HDF5 snapshot into a particle array  
!========================================================================
subroutine read_Ghdf5_particles()

  character(clen), parameter :: myname="read_Ghdf5_particles" 
  logical, parameter :: crash=.true.
  integer, parameter :: verb=2
  character(clen) :: str,fmt
  
  real(r4b), allocatable :: rblck(:)
  real(r4b), allocatable :: rblck3(:,:)
  integer(i4b), allocatable :: iblck(:)
  integer(i8b) :: ngasread

  type(gadget_header_type) :: ghead
  character(clen) :: snapfile, VarName, GroupName
  integer(i8b) :: i
  integer(i4b) :: err,fh

  integer(i8b) :: npar, ngas, nmass
  integer(i8b) :: npar1, ngas1, nmass1
  logical :: varmass(6)
  integer(i8b) :: fn

  real(r8b) :: meanweight
  logical :: caseA(2)
  real(r8b) :: Tdum
  real(r8b) :: MB
  real(r8b) :: Hmf
  real(r8b) :: nH8
  real(r4b) :: nH4
  real(r4b) :: temp

  type(ion_table_type) :: itab
  real :: redshift



  ! set local particle numbers
  !============================
  ghead = saved_gheads( GV%CurSnapNum, 0 )
  varmass = (ghead%npar_all > 0 .and. ghead%mass == 0)
  npar = sum(ghead%npar_all)
  ngas = ghead%npar_all(1)
  nmass = sum(ghead%npar_all, mask=varmass)

  ! do Gadget dummy checks
  !============================
  if (ngas == 0) call myerr("snapshot has no gas particles",myname,crash)

  ! calculate bytes per particle and allocate particle array
  !===========================================================
  MB = GV%bytesperpar * real(ngas) / 2.**20
  GV%MB = GV%MB + MB

  fmt="(A,F10.4,A,I10,A)"
  write(str,fmt) "  allocating ", MB, " MB for ", ngas, " particles"
  call mywrite(str,verb) 

  allocate (psys%par(ngas), stat=err)
  if (err /= 0) call myerr("failed to allocate par",myname,crash)


  ! now read all snapshot files
  !==============================          
  ngasread = 0
  GroupName = 'PartType0/'
  files: do fn = 0, ghead%nfiles-1


     ! recall the header info
     !-----------------------------------------------------------!  
     ghead   = saved_gheads( GV%CurSnapNum, fn )
     varmass = (ghead%npar_file > 0 .and. ghead%mass == 0)
     npar1   = sum(ghead%npar_file)
     ngas1   = ghead%npar_file(1)
     nmass1  = sum(ghead%npar_file, mask=varmass)
     if (ngas1 == 0) cycle


     ! begin read
     !-----------------------------------------------------------!  
     call form_gadget_snapshot_file_name(GV%SnapPath,GV%ParFileBase,GV%CurSnapNum,fn,snapfile)
     call mywrite("   reading particle snapshot file: "//trim(snapfile), verb, fmt="(A)")
     call hdf5_open_file(fh, snapfile, readonly=.true.)


     ! read positions 
     !-----------------------------------------------------------!  
     allocate(rblck3(3,ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck3 for pos",myname,crash)
     VarName = 'Coordinates'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck3)
     forall(i=1:ngas1) psys%par(ngasread+i)%pos(1) = rblck3(1,i)
     forall(i=1:ngas1) psys%par(ngasread+i)%pos(2) = rblck3(2,i)
     forall(i=1:ngas1) psys%par(ngasread+i)%pos(3) = rblck3(3,i)
     deallocate(rblck3)

     ! read velocities 
     !-----------------------------------------------------------!  
#ifdef incVel
     allocate(rblck3(3,ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck3 for vel",myname,crash)
     VarName = 'Velocity'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck3)
     forall(i=1:ngas1) psys%par(ngasread+i)%vel(1) = rblck3(1,i)
     forall(i=1:ngas1) psys%par(ngasread+i)%vel(2) = rblck3(2,i)
     forall(i=1:ngas1) psys%par(ngasread+i)%vel(3) = rblck3(3,i)
     deallocate(rblck3)
#endif

     ! read id's 
     !-----------------------------------------------------------!  
     allocate(iblck(ngas1), stat=err )
     if(err/=0) call myerr("allocating iblck for ID",myname,crash)
     VarName = 'ParticleIDs'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),iblck)
     forall(i=1:ngas1) psys%par(ngasread+i)%id = iblck(i)
     deallocate(iblck)

     ! read masses 
     !-----------------------------------------------------------!  

     ! if gas particles are variable mass
     if (varmass(1)) then  
        allocate(rblck(ngas1), stat=err)
        if(err/=0) call myerr("allocating rblck for mass",myname,crash)
        VarName = 'Mass'
        call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
        forall(i=1:ngas1) psys%par(ngasread+i)%mass = rblck(i)
        deallocate(rblck)

     ! if gas particles are isomass
     else
        psys%par(ngasread+1:ngasread+ngas1)%mass = ghead%mass(1)       
     end if


     ! read temperature
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for T",myname,crash)
     VarName = 'Temperature'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
     forall(i=1:ngas1) psys%par(ngasread+i)%T = rblck(i)
     deallocate(rblck)
 

     ! read EOS
     !-----------------------------------------------------------!  
#ifdef incEOS
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for EOS",myname,crash)
     VarName = 'OnEquationOfState'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
     forall(i=1:ngas1) psys%par(ngasread+i)%eos = rblck(i)
     deallocate(rblck)
#endif

     ! read density 
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for rho",myname,crash)
     VarName = 'Density'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
     forall(i=1:ngas1) psys%par(ngasread+i)%rho = rblck(i)
     deallocate(rblck)


     ! read smoothing lengths 
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for hsml",myname,crash)
     VarName = 'SmoothingLength'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
     forall(i=1:ngas1) psys%par(ngasread+i)%hsml = rblck(i)
     deallocate(rblck)


     ! read Hydrogen mass fractions
     !-----------------------------------------------------------!  
#ifdef incHmf
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for Hmf",myname,crash)
     VarName = 'ElementAbundance/Hydrogen'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
     forall(i=1:ngas1) psys%par(ngasread+i)%Hmf = rblck(i)
     deallocate(rblck)
#endif


     ! read Helium mass fractions
     !-----------------------------------------------------------!  
#ifdef incHemf
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for Hemf",myname,crash)
     VarName = 'ElementAbundance/Helium'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
     forall(i=1:ngas1) psys%par(ngasread+i)%Hemf = rblck(i)
     deallocate(rblck)
#endif


     ngasread = ngasread + ngas1
     call hdf5_close_file(fh)


  end do files



  ! calculate xHI from CLOUDY iontables
  !-----------------------------------------------------------!  
  write(*,*) 
  write(*,*) " calculating xHI from CLOUDY tables"

  call read_ion_table_file( "../data/ionization_tables/h1.hdf5", itab )
  redshift = ghead%z

  ! first get the gammaHI from the uniform UVB at this redshift
  GV%UVB_gammaHI_cloudy = return_gammaHI_at_z( itab, redshift )

 
  do i = 1,ngas

#ifdef incHmf
     Hmf = psys%par(i)%Hmf
#else
     Hmf = GV%H_mf
#endif

     nH8 = psys%par(i)%rho * GV%cgs_rho * ghead%h**2 / ghead%a**3 * &
           Hmf / gconst%PROTONMASS
     nH4 = nH8

     temp = psys%par(i)%T

#ifdef incEOS
     if (GV%EOStemp > 0.0) then
        if (psys%par(i)%eos == 1.0) then
           temp = GV%EOStemp
        endif
     endif
#endif

     psys%par(i)%xHI = &
          interpolate_ion_table( itab, redshift, log10(temp), log10(nH4) )
  end do



#ifdef cloudy
  psys%par(:)%xHI_cloudy = psys%par(:)%xHI
  write(*,*) " min/max xHI_cloudy = ", minval( psys%par%xHI_cloudy ), maxval( psys%par%xHI_cloudy )
#endif

  write(*,*) " min/max xHI        = ", minval( psys%par%xHI ), maxval( psys%par%xHI )

  

  ! set xHII from xHI 
  !-----------------------------------------------------------!  
  psys%par%xHII = 1.0d0 - psys%par%xHI


  ! set caseA true or false for collisional equilibrium
  !-----------------------------------------------------
  caseA = .false.
  if (.not. GV%OnTheSpotH  .or. GV%HydrogenCaseA) caseA(1) = .true.
  if (.not. GV%OnTheSpotHe .or. GV%HeliumCaseA)   caseA(2) = .true.


  ! if Helium, initialize ionization fractions to collisional equilibrium
  !------------------------------------------------------------------------
#ifdef incHe
  call set_collisional_ionization_equilibrium(psys, caseA, GV%IsoTemp, DoHydrogen=.false., fit="hui")
#endif


  ! set the electron fractions from the ionization fractions
  !----------------------------------------------------------
  call set_ye(psys, GV%H_mf, GV%He_mf, GV%NeBackGround)



end subroutine read_Ghdf5_particles




subroutine gadget_output_hdf5()

  stop "gadget_output_hdf5 not implemented"

end subroutine gadget_output_hdf5



#endif


end module gadget_input_hdf5_mod
