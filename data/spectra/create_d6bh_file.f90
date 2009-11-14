! this module contains functions to make powerlaw and thermal spectra
! the files that Sphray reads are actually cumulative distribution 
! functions and so are calculated by integrating the spectra.
! The user can choose (by setting the parameters in the program
! section below) how many spectra to put in a file, their type
! (thermal or powerlaw), the frequency limits, and the name of the
! file that stores them.  Two files will be produced.  One labeled
! 'spec' is for plotting purposes and contains the frequency vs 
! normalized intensity.  The other is labeled 'cdf' and is the one
! that should be indicated in the sphray config file.
!======================================================================

module specfuncs
use myf90_mod
use physical_constants_mod
implicit none

private

public :: set_slope
public :: powerlaw
public :: set_temperature
public :: blackbody

public :: pi, zero, one, two, three, four, six, eight, ten, half


real(r8b), parameter :: zero = 0.0_r8b
real(r8b), parameter :: one = 1.0_r8b
real(r8b), parameter :: two = 2.0_r8b
real(r8b), parameter :: three = 3.0_r8b
real(r8b), parameter :: four = 4.0_r8b
real(r8b), parameter :: six = 6.0_r8b
real(r8b), parameter :: eight = 8.0_r8b
real(r8b), parameter :: ten = 10.0_r8b
real(r8b), parameter :: half = 0.5_r8b

real(r8b) :: alpha
real(r8b) :: temperature
real(r8b) :: kT_eV
real(r8b) :: kT_erg
real(r8b) :: nuwien

contains

! power law spectrum. input is frequency in units of 
! the hydrogen ionizing frequency. 
!----------------------------------------------------
function powerlaw(nu) result(P)
  real(r8b), dimension(:), intent(in) :: nu
  real(r8b), dimension(size(nu)) :: P
  real(r8b), dimension(size(nu)) :: nuHz
  real(r8b) :: nufac

  nuHz = nu * HI_th_Hz
  P = nu**(alpha)
 
end function powerlaw

! sets the slope for the power law spectra
!-------------------------------------------
subroutine set_slope(a)
  real(r8b), intent(in) :: a
  alpha = a
end subroutine set_slope

! black body spectrum. input is frequency in units of 
! the hydrogen ionizing frequency. 
!----------------------------------------------------
function blackbody(nu) result(B)
  real(r8b), dimension(:), intent(in) :: nu
  real(r8b), dimension(size(nu)) :: B
  real(r8b), dimension(size(nu)) :: nuHz
  real(r8b) :: nufac

  nuHz = nu * HI_th_Hz
  nufac = HI_th_Hz * h_eV_s / kT_eV
  B = two * h_eV_s * nuHz * nuHz * nuHz / ( (c * c) * (exp(nu * nufac) - one) )
 
end function blackbody

! sets the temperature for the blackbody spectrum
!--------------------------------------------------
subroutine set_temperature(T)
  real(r8b), intent(in) :: T
  temperature = T
  kT_eV  = k_eV_K  * T
  kT_erg = k_erg_K * T
  nuwien = 5.879e10_r8b * T
  write(*,*) " nuwien = ", nuwien / HI_th_Hz
end subroutine set_temperature


end module specfuncs


program create_spectra_file
use myf90_mod
use physical_constants_mod
use specfuncs
use nr, only: qromb
implicit none

! these parameters are meant to be adjusted by the user
! Nspectra =  number of spectra to put in the file
! NuMax = an array indicating the maximum freq to integrate
!         up to in units of HI_th
! SpecType = an array, t = thermal, p = power law
! SpecParam = an array, for SpecType=t this is temperature
!                       for SpecType=p this is slope
! filebase = the base of the file the spectra will be written to
!            <filebase>.spec and <filebase>.cdf will be produced
!=================================================================

integer(i4b), parameter   :: Nspectra = 2
real(r8b), parameter      :: NuMax(Nspectra) = (/ 36.0d0, 36.0d0 /)
character(1), parameter   :: SpecType(Nspectra) = (/ "p", "p" /)
real(r8b), parameter      :: SpecParam(Nspectra) = (/ -2.0, -1.76  /)
character(200), parameter :: filebase = "d6bh"


real(r8b), parameter :: numin = one
integer(i4b), parameter :: nuevals = 1001
real(r8b), parameter :: dnu(Nspectra) = (NuMax - numin) / (nuevals-1)


real(r8b) :: nuimp
real(r8b) :: kint
real(r8b) :: slope
real(r8b) :: temp
real(r8b) :: nu
real(r8b) :: intensity
real(r8b) :: norm
real(r8b) :: cdf

integer(i4b) :: i, j

character(200) :: specfile 
character(200) :: cdffile 

specfile = trim(filebase) // ".spec"
cdffile = trim(filebase) // ".cdf"

write(*,*)
write(*,*) "==============================================================="
write(*,*) "     Creating Spectra File                                     "
write(*,*) "==============================================================="
write(*,*)  
write(*,*) " number of spectra:     ", Nspectra
write(*,*) " number of evaluations: ", nuevals
write(*,*) " min nu (Hth units):    ", numin
write(*,*) " max nu (Hth units):    ", numax
write(*,*) 
write(*,*) " spectra file (for plotting) : ", trim(specfile)
write(*,*) " cdf file (read in by Sphray): ", trim(cdffile)
write(*,*) 


100 format(A,I3,ES12.5,2X,A,2X,ES12.5)
write(*,*) " spectrum descriptions"
do i = 1,Nspectra
   write(*,100) "  i,numax,type,slope/temp:", i, NuMax(i),SpecType(i),SpecParam(i)
end do
write(*,*) 



open(unit=10,file=specfile)
open(unit=11,file=cdffile)

write(11,*) Nspectra

do i = 1,Nspectra


   write(11,*) nuevals, numin, NuMax(i)   

   if (SpecType(i)=="p") then
      slope = SpecParam(i) 
      call set_slope(slope)

      norm = qromb(powerlaw, numin, NuMax(i))
      write(*,*) " i, norm = ", i, norm

      do j = 0,nuevals-1
         nu = numin + j * dnu(i)
         intensity = sum( powerlaw( (/nu/) ) )
         write(10,*) nu, intensity/norm
         if (j==nuevals) then
            cdf = norm
         else
            cdf = qromb(powerlaw, numin, nu)
         end if
         write(11,*) nu, cdf/norm
      end do

   else if (SpecType(i)=="t") then
      temp = SpecParam(i)
      call set_temperature(temp)

      norm = qromb(blackbody, numin, NuMax(i))
      write(*,*) " i, norm = ", i, norm

      do j = 0,nuevals-1
         nu = numin + j * dnu(i)
         intensity = sum( blackbody( (/nu/) ) )
         write(10,*) nu, intensity/norm
         if (j==nuevals) then
            cdf = norm
         else
            cdf = qromb(blackbody, numin, nu)
         end if
         write(11,*) nu, cdf/norm
      end do

   else
      write(*,*) "SpecType(i)=", SpecType(i), " not recognized"
      stop
      
   end if

end do

close(10)
close(11)



end program create_spectra_file
