!> \file gadget_input_hdf5.F90

!> \brief Handles readin of GADGET HDF5 formatted files
!<

module gadget_input_hdf5_mod
use global_mod, only: GV
use gadget_input_mod
#ifdef hdf5
use hdf5_wrapper
#endif
implicit none



contains

#ifndef hdf5

! these are dummy subroutines so that the calls outside this file
! dont have to be wrapped with pre processor macros.

subroutine get_planning_data_gadget_hdf5()
  call myerr("this routine shuold not have been called","hdf5dummy",crash=.true.)
end subroutine get_planning_data_gadget_hdf5

subroutine read_Ghdf5_particles()
  call myerr("this routine shuold not have been called","hdf5dummy",crash=.true.)
end subroutine read_Ghdf5_particles

subroutine update_hdf5_particles()
  call myerr("this routine shuold not have been called","hdf5dummy",crash=.true.)
end subroutine update_hdf5_particles

subroutine gadget_output_hdf5()
  call myerr("this routine shuold not have been called","hdf5dummy",crash=.true.)
end subroutine gadget_output_hdf5


#else

!>   reads in an HDF5 particle header
!==========================================================================
subroutine read_gadget_header_hdf5(snapfile, ghead, fh, closefile)
  character(*), intent(in) :: snapfile
  type(gadget_header_type), intent(out) :: ghead
  integer(i4b), intent(out) :: fh
  logical, optional :: closefile

  logical :: closef

  if (.not. present(closefile) ) then
     closef = .true.
  else
     closef = closefile
  end if

  call hdf5_open_file(fh,snapfile)
  call hdf5_read_attribute(fh,'Header/NumPart_ThisFile',ghead%npar_file)
  call hdf5_read_attribute(fh,'Header/MassTable',ghead%mass)
  call hdf5_read_attribute(fh,'Header/ExpansionFactor',ghead%a)
  call hdf5_read_attribute(fh,'Header/Redshift',ghead%z)
  call hdf5_read_attribute(fh,'Header/Flag_Sfr',ghead%flag_sfr)
  call hdf5_read_attribute(fh,'Header/Flag_Feedback',ghead%flag_feedback)
  call hdf5_read_attribute(fh,'Header/NumPart_Total',ghead%npar_all)
  call hdf5_read_attribute(fh,'Header/Flag_Cooling',ghead%flag_cooling)
  call hdf5_read_attribute(fh,'Header/NumFilesPerSnapshot',ghead%nfiles)
  call hdf5_read_attribute(fh,'Header/BoxSize',ghead%boxlen)
  call hdf5_read_attribute(fh,'Header/Omega0',ghead%OmegaM)
  call hdf5_read_attribute(fh,'Header/OmegaLambda',ghead%OmegaL)
  call hdf5_read_attribute(fh,'Header/HubbleParam',ghead%h)
  call hdf5_read_attribute(fh,'Header/Flag_StellarAge',ghead%flag_age)
  call hdf5_read_attribute(fh,'Header/Flag_Metals',ghead%flag_metals)
  call hdf5_read_attribute(fh,'Header/NumPart_Total_HighWord',ghead%npar_hw)
  call hdf5_read_attribute(fh,'Header/OmegaBaryon',ghead%OmegaB)
  if (closef) call hdf5_close_file(fh)



end subroutine read_gadget_header_hdf5


!>   gets run planning data from Gadget Headers
!===============================================
subroutine get_planning_data_gadget_hdf5()
use cosmology_mod, only: tsinceBB
use source_input_mod, only: read_source_header

  integer(i4b) :: fh      !< hdf5 file handle
  integer(i8b) :: iSnap   !< first snap number
  integer(i8b) :: fSnap   !< last snap number
  integer(i8b) :: pfiles  !< nfiles for particle snapshots
  integer(i8b) :: sfiles  !< nfiles for source snapshosts
  integer(i8b) :: i,j     !< counters

  type(gadget_header_type) :: ghead
  type(source_header_type) :: shead

  character(clen) :: snapfile

  real(r8b) :: Time_GYR
  integer(i8b) :: lun
  integer(i8b) :: loglun
  character(clen) :: logfile
 
  ! open up the planning data log file
  !======================================================
  logfile = trim(GV%OutputDir) // "/" // "headers.log"
  call open_formatted_file_w(logfile,loglun)

  
  ! these global variables are read from the config file
  !======================================================
  iSnap = GV%StartSnapNum
  fSnap = GV%EndSnapNum
    
  pfiles = GV%ParFilesPerSnap
  sfiles = GV%SourceFilesPerSnap

  allocate( saved_gheads(iSnap:fSnap, 0:pfiles-1) )


  ! read all particle headers and write to log file
  !===================================================
  write(loglun,'(A)') "reading all GADGET particle header(s) ... "
  do i = iSnap,fSnap
     do j = 0,pfiles-1

        call form_Gsnapshot_file_name(GV%SnapPath,GV%ParFileBase,i,j,snapfile)
        write(snapfile,"(A,A)") trim(snapfile), ".hdf5"
        write(loglun,'(I3,"  ",A)') i,trim(snapfile)

        call read_gadget_header_hdf5(snapfile, ghead, fh, closefile=.false.)

        call hdf5_read_attribute(fh,'Parameters/StellarEvolutionParameters/SF_EOSGammaEffective',GV%sf_gamma_eos)

        call hdf5_read_attribute(fh,'Units/UnitLength_in_cm',GV%cgs_len)
        call hdf5_read_attribute(fh,'Units/UnitMass_in_g',GV%cgs_mass)
        call hdf5_read_attribute(fh,'Units/UnitVelocity_in_cm_per_s',GV%cgs_vel)
        call hdf5_read_attribute(fh,'Units/UnitTime_in_s',GV%cgs_time)        
        call hdf5_read_attribute(fh,'Units/UnitDensity_in_cgs',GV%cgs_rho)
        call hdf5_read_attribute(fh,'Units/UnitEnergy_in_cgs',GV%cgs_enrg)
        call hdf5_read_attribute(fh,'Units/UnitPressure_in_cgs',GV%cgs_prs)
        GV%cgs_lum  = GV%cgs_enrg / GV%cgs_time

        call hdf5_read_attribute(fh,'Header/Time_GYR', Time_GYR)

        call hdf5_close_file(fh)
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
                
        if (GV%Comoving) then
           PLAN%snap(i)%ScalefacAt = ghead%a
           PLAN%snap(i)%TimeAt = Time_GYR * Gyr2sec                              ! in seconds
           PLAN%snap(i)%TimeAt = PLAN%snap(i)%TimeAt * GV%LittleH / GV%cgs_time  ! in code units
        else
           PLAN%snap(i)%ScalefacAt = 1.0d0 / (1.0d0 + ghead%z)
           PLAN%snap(i)%TimeAt = Time_GYR * Gyr2sec                              ! in seconds
           PLAN%snap(i)%TimeAt = PLAN%snap(i)%TimeAt * GV%LittleH / GV%cgs_time  ! in code units
        end if
        
     end do
  end do



  ! read all source headers and write to log file
  !===================================================
  write(loglun,*) 
  write(loglun,'(A)') "reading all Gadget source header(s) ... "
  do i = iSnap,fSnap
     do j = 1,sfiles
        call form_snapshot_file_name(GV%SourcePath,GV%SourceFileBase,i,j,snapfile)
        write(loglun,'(I3,"  ",A)') i,trim(snapfile)
        call read_source_header(snapfile,shead,lun,closefile=.true.)
        
        PLAN%snap(i)%RaysFromSrcHeader = shead%TotalRays
        GV%Lunit = shead%Lunit
        
     end do
  end do


    ! close headers log file
  !========================
  close(loglun)


  ! write units to log file
  !===================================================

  logfile = trim(GV%OutputDir) // "/" // "code_units.log"
  call open_formatted_file_w(logfile,loglun)

  106    format(A,ES12.5,A)
  write(loglun,*) 
  write(loglun,'(A)') "setting code units ..."
  write(loglun,106) "  Hubble:     ", GV%LittleH, " = H0[km/s/Mpc] / 100 "
  write(loglun,106) "  length:     ", GV%cgs_len / kpc2cm, " [kpc/h]"
  write(loglun,106) "  mass:       ", GV%cgs_mass / M_solar, " [Msun/h]"
  write(loglun,106) "  velocity    ", GV%cgs_vel / km2cm, " [km/s]"
  write(loglun,106) "  time:       ", GV%cgs_time / Myr2sec, " [Myr/h]"
  write(loglun,106) "  density:    ", GV%cgs_rho / (M_solar / kpc2cm**3), " [Msun/kpc^3 h^2]"
  write(loglun,106) "  pressure:   ", GV%cgs_prs, " [dyne/cm^2 h^2]"
  write(loglun,106) "  energy:     ", GV%cgs_enrg, " [ergs/h]"
  write(loglun,106) "  luminosity: ", GV%cgs_lum / L_solar, " [Lsun]"
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
  real(r4b), allocatable :: eosblck(:)
  real(r4b), allocatable :: rblck3(:,:)
  integer(i4b), allocatable :: iblck(:)
  integer(i8b) :: ngasread

  type(gadget_header_type) :: ghead
  character(clen) :: snapfile, VarName, GroupName
  integer(i8b) :: i
  integer(i4b) :: err,fh

  integer(i8b) :: npar, ngas, nmass
  integer(i8b) :: npar1, ngas1, nmass1
  logical :: closefile
  logical :: header_verbose
  logical :: varmass(6)
  integer(i8b) :: fn

  real(r8b) :: meanweight
  logical :: caseA(2)
  real(r8b) :: xvec(5)
  real(r8b) :: Tdum
  real(r8b) :: MB
  real(r8b) :: Hmf

  ! set default values
  !===============================
  caseA=.false.
  xvec=0.0

  ! set local particle numbers
  !============================
  ghead = saved_gheads( GV%CurSnapNum, 0 )
  varmass = (ghead%npar_all > 0 .and. ghead%mass == 0)
  npar = sum(ghead%npar_all)
  ngas = ghead%npar_all(1)
  nmass = sum(ghead%npar_all, mask=varmass)

  ! do Gadget dummy checks
  !============================
  if (ngas .EQ. 0) call myerr("snapshot has no gas particles",myname,crash)

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

     call form_Gsnapshot_file_name(GV%SnapPath,GV%ParFileBase,GV%CurSnapNum,fn,snapfile)
     write(snapfile,"(A,A)") trim(snapfile), ".hdf5"
     call mywrite("   reading particle snapshot file: "//trim(snapfile), verb, fmt="(A)")
     call hdf5_open_file(fh,snapfile)

     ghead = saved_gheads( GV%CurSnapNum, fn )

     varmass = (ghead%npar_file > 0 .and. ghead%mass == 0)
     npar1 = sum(ghead%npar_file)
     ngas1 = ghead%npar_file(1)
     nmass1 = sum(ghead%npar_file, mask=varmass)

     if (ngas1 == 0) cycle

     ! read positions 
     !-----------------------------------------------------------!  
     allocate(rblck3(3,ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck3 for pos",myname,crash)

     VarName = 'Coordinates'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck3)

     do i = 1,ngas1
        psys%par(ngasread+i)%pos(1) = rblck3(1,i)
        psys%par(ngasread+i)%pos(2) = rblck3(2,i)
        psys%par(ngasread+i)%pos(3) = rblck3(3,i)
     end do
  
     deallocate(rblck3)

     ! read velocities 
     !-----------------------------------------------------------!  
#ifdef incVel
     allocate(rblck3(3,ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck3 for vel",myname,crash)

     VarName = 'Velocity'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck3)


     do i = 1,ngas1
        psys%par(ngasread+i)%vel(1) = rblck3(1,i)
        psys%par(ngasread+i)%vel(2) = rblck3(2,i)
        psys%par(ngasread+i)%vel(3) = rblck3(3,i)
     end do
  
     deallocate(rblck3)
#endif

     ! read id's 
     !-----------------------------------------------------------!  
     allocate(iblck(ngas1), stat=err )
     if(err/=0) call myerr("allocating iblck for ID",myname,crash)

     VarName = 'ParticleIDs'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),iblck)

     do i = 1,ngas1
        psys%par(ngasread+i)%id = iblck(i)
     end do

     deallocate(iblck)


     ! read masses 
     !-----------------------------------------------------------!  
     if (varmass(1)) then  ! if gas particles are variable mass
        allocate(rblck(ngas1), stat=err)
        if(err/=0) call myerr("allocating rblck for mass",myname,crash)

        VarName = 'Mass'
        call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)

        do i = 1,ngas1
           psys%par(ngasread+i)%mass = rblck(i)
        end do

        deallocate(rblck)
 
     else

        psys%par(ngasread+1:ngasread+ngas1)%mass = ghead%mass(1)
        
     end if


     ! read temperature
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for T",myname,crash)

     VarName = 'Temperature'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)

#ifdef OWLS
     allocate(eosblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating eosblck for EOS",myname,crash)

     VarName = 'OnEquationOfState'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),eosblck)
#endif

     do i = 1,ngas1
        psys%par(ngasread+i)%T = rblck(i) 
#ifdef OWLS
        if ( eosblck(i) == 1.0 ) psys%par(ngasread+i)%T = 1.0e4
#endif
     end do

     deallocate(rblck)
     if (allocated(eosblck)) deallocate(eosblck)


     ! read density 
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for rho",myname,crash)

     VarName = 'Density'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)

     do i = 1, ngas1
        psys%par(ngasread+i)%rho = rblck(i)
     end do
     
     deallocate(rblck)


     ! read smoothing lengths 
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for hsml",myname,crash)

     VarName = 'SmoothingLength'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)

     do i = 1, ngas1
        psys%par(ngasread+i)%hsml = rblck(i)
     end do
     
     deallocate(rblck)


     ! read Hydrogen mass fractions
     !-----------------------------------------------------------!  
#ifdef incHmf
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for Hmf",myname,crash)

     VarName = 'ElementAbundance/Hydrogen'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)

     do i = 1, ngas1
        psys%par(ngasread+i)%Hmf = rblck(i)
     end do
     
     deallocate(rblck)
#endif


     ! read Helium mass fractions
     !-----------------------------------------------------------!  
#ifdef incHemf
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for Hemf",myname,crash)

     VarName = 'ElementAbundance/Helium'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)

     do i = 1, ngas1
        psys%par(ngasread+i)%Hemf = rblck(i)
     end do
     
     deallocate(rblck)
#endif


     ! read xHI (really HI mass and convert) 
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for xHI",myname,crash)

     VarName = 'IonMass/h1'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)

     do i = 1, ngas1
#ifdef incHmf
        Hmf = psys%par(ngasread+i)%Hmf
#else
        Hmf = GV%H_mf
#endif
        psys%par(ngasread+i)%xHI = rblck(i) / ( psys%par(ngasread+i)%mass * Hmf )
     end do
     
     deallocate(rblck)





     ngasread = ngasread + ngas1
     call hdf5_close_file(fh)


  end do files

 
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
  call set_x_eq(psys, caseA, GV%IsoTemp, JustHe=.true.)


  ! set the electron fractions from the ionization fractions
  !----------------------------------------------------------
  call set_ye(psys, GV%H_mf, GV%He_mf)



end subroutine read_Ghdf5_particles





subroutine update_hdf5_particles()
end subroutine update_hdf5_particles



subroutine gadget_output_hdf5()

  stop "gadget_output_hdf5 not implemented"

end subroutine gadget_output_hdf5



#endif


end module gadget_input_hdf5_mod
