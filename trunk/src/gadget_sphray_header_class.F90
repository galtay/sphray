!> \file gadget_sphray_header_class.f90

!> \brief Handles Gadget style headers with extra SPHRAY info.  
!!
!!
!< 

module gadget_sphray_header_class
use myf03_mod
use gadget_public_header_class
use gadget_owls_header_class_hdf5
implicit none
private

public :: gadget_sphray_header_type

!> Gadget header type
!-----------------------------------------
type gadget_sphray_header_type
   integer(i4b) :: npar_file(0:5)  !< number of particles in snapshot file
   real(r8b)    :: mass(0:5)       !< mass of each particle type if constant
   real(r8b)    :: a               !< scale factor or time
   real(r8b)    :: z               !< redshift
   integer(i4b) :: flag_sfr        !< flag for star formation
   integer(i4b) :: flag_feedback   !< flag for feedback
   integer(i4b) :: npar_all(0:5)   !< number of particles in whole snapshot
   integer(i4b) :: flag_cooling    !< flag for radiative cooling
   integer(i4b) :: nfiles          !< number of files in a this snapshot
   real(r8b)    :: boxlen          !< box length
   real(r8b)    :: OmegaM          !< omega matter
   real(r8b)    :: OmegaL          !< omega lambda
   real(r8b)    :: h               !< little hubble
   integer(i4b) :: flag_age        !< flag for stellar age
   integer(i4b) :: flag_metals     !< flag for metallicity
   integer(i4b) :: npar_hw(0:5)    !< 64 bit part of npar 
   integer(i4b) :: flag_entr_ics   !< ICs contain entropy instead of energy

   real(r8b)    :: OmegaB          !< omega baryon 
   integer(i8b) :: rays_traced     !< number of rays traced
   integer(i4b) :: flag_Hmf        !< output has Hydorgen mass fraction? 
   integer(i4b) :: flag_Hemf       !< output has Helium mass fraction? 
   integer(i4b) :: flag_helium     !< output has xHeI and xHeII?
   integer(i4b) :: flag_gammaHI    !< output has HI photoionization rate? 
   integer(i4b) :: flag_cloudy     !< output has Cloudy xHIeq? 
   integer(i4b) :: flag_eos        !< output has EOS info?
   integer(i4b) :: flag_incsfr     !< output has SFR info?
   real(r8b)    :: time_gyr        !< time since BB from cosmo variables [Gyr]
   integer(i4b) :: unused(2)       !< spacer
 contains
   procedure :: copy_ghead_public => ghead_to_gshead
   procedure :: copy_ghead_owls => ohead_to_gshead
   procedure :: write_lun => write_header_to_lun
end type gadget_sphray_header_type


contains


!> writes a gadget header to an open file
!--------------------------------------------------------------
subroutine write_header_to_lun(this, lun)
  class(gadget_sphray_header_type) :: this
  integer(i4b), intent(in) :: lun

  write(lun) this%npar_file(:), this%mass(:), this%a, this%z, &              
       this%flag_sfr, this%flag_feedback, this%npar_all(:), &    
       this%flag_cooling, this%nfiles, this%boxlen, this%OmegaM, &         
       this%OmegaL, this%h, this%flag_age, this%flag_metals, &    
       this%npar_hw(:), this%flag_entr_ics, this%OmegaB, &
       this%rays_traced, this%flag_Hmf, this%flag_Hemf, this%flag_helium, &
       this%flag_gammaHI, this%flag_cloudy, this%flag_eos, this%flag_sfr, &
       this%time_gyr, this%unused(:)      

end subroutine write_header_to_lun



subroutine ohead_to_gshead( this, ohead )
  class(gadget_sphray_header_type) :: this
  class(gadget_owls_header_type) :: ohead

  this%npar_file(0:5)   = ohead%npar_file(0:5)   
  this%mass(0:5)        = ohead%mass(0:5)        
  this%a                = ohead%a                
  this%z                = ohead%z                
  this%flag_sfr         = ohead%flag_sfr         
  this%flag_feedback    = ohead%flag_feedback    
  this%npar_all(0:5)    = ohead%npar_all(0:5)    
  this%flag_cooling     = ohead%flag_cooling     
  this%nfiles           = ohead%nfiles           
  this%boxlen           = ohead%boxlen           
  this%OmegaM           = ohead%OmegaM           
  this%OmegaL           = ohead%OmegaL           
  this%h                = ohead%h                
  this%flag_age         = ohead%flag_age         
  this%flag_metals      = ohead%flag_metals      
  this%npar_hw(0:5)     = ohead%npar_hw(0:5)     

  this%OmegaB           = ohead%OmegaB
  this%time_gyr         = ohead%time_gyr

end subroutine ohead_to_gshead



subroutine ghead_to_gshead( this, ghead )
  class(gadget_sphray_header_type) :: this
  class(gadget_public_header_type) :: ghead

  this%npar_file(0:5)   = ghead%npar_file(0:5)   
  this%mass(0:5)        = ghead%mass(0:5)        
  this%a                = ghead%a                
  this%z                = ghead%z                
  this%flag_sfr         = ghead%flag_sfr         
  this%flag_feedback    = ghead%flag_feedback    
  this%npar_all(0:5)    = ghead%npar_all(0:5)    
  this%flag_cooling     = ghead%flag_cooling     
  this%nfiles           = ghead%nfiles           
  this%boxlen           = ghead%boxlen           
  this%OmegaM           = ghead%OmegaM           
  this%OmegaL           = ghead%OmegaL           
  this%h                = ghead%h                
  this%flag_age         = ghead%flag_age         
  this%flag_metals      = ghead%flag_metals      
  this%npar_hw(0:5)     = ghead%npar_hw(0:5)     
  this%flag_entr_ics    = ghead%flag_entr_ics    

end subroutine ghead_to_gshead


end module gadget_sphray_header_class
