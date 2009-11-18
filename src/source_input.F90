!> \file source_input.F90

!> \brief Handles readin of sourcefiles
!<

module source_input_mod
use myf90_mod
use particle_system_mod, only: source_type
use global_mod, only: global_variables_type
use global_mod, only: psys
implicit none

!> source header type 
!-------------------------
type source_header_type 
   integer(i8b) :: NsrcSnap  !< number of sources in snapshot
   integer(i8b) :: NsrcFile  !< number of sources in file
   integer(i8b) :: NumFiles  !< number of files per snapshot
   real(r8b) :: TotalRays    !< sum of all rays to be cast this snap
   real(r8b) :: Lunit        !< Luminosities are photons/s * Lunit
end type source_header_type

contains

!> reads in a source header, prints the header to the screen if verbose
!========================================================================
subroutine read_source_header(snapfile, verbose, closefile, shead, lun)
  character(*), intent(in) :: snapfile  !< file containing source header
  logical, intent(in) :: verbose        !< report to screen?
  logical, intent(in) :: closefile      !< close file when done?
  type(source_header_type), intent(inout) :: shead  !< source header to read
  integer(i8b), intent(out) :: lun  !< output lun assigned to snapfile

  call open_formatted_file_r(snapfile,lun)

  read(lun,*) shead%NsrcSnap
  read(lun,*) shead%NsrcFile
  read(lun,*) shead%NumFiles
  read(lun,*) shead%TotalRays
  read(lun,*) shead%Lunit
  read(lun,*)
  read(lun,*)
  if (closefile) close(lun)
  
  if (verbose) call source_header_to_screen(shead)
  
end subroutine read_source_header


!> reads a source snapshot into arc array (src is allocated in this routine)
!=============================================================================
subroutine read_src_snapshot(GV,verbose,MB)
  type(global_variables_type), intent(in) :: GV
  logical, intent(in) :: verbose      !< report to screen?
  real(r8b), intent(out) :: MB  !< MB's allocated to source storage

  character(200) :: snapfile  
  type(source_header_type) :: shead
  integer(i8b) :: lun, err
  integer(i8b) :: i, N, Nall
  integer(i8b) :: bytespersource
  logical :: closefile
  integer(i8b) :: fn


  fn=1
  call form_snapshot_file_name(GV%SourcePath,GV%SourceFileBase,GV%CurSnapNum,fn,snapfile)
      
  bytespersource = 10 * 4 + 8
  closefile = .false.
  call read_source_header(snapfile,verbose,closefile,shead,lun)

  Nall = shead%NsrcSnap        
  N = shead%NsrcFile
  MB = bytespersource * real(Nall) / 2**20

  101 format(A,F10.4,A,I10,A)
  write(*,101) "allocating ", MB, " MB for ", Nall, " sources"

  if (allocated(psys%src)) deallocate(psys%src)
  allocate (psys%src(Nall),stat=err)
  if (err /= 0) stop "cant allocate sources in particle system"

  write(*,'(A,A)') "reading source snapshot file: ", trim(snapfile)
     
  do i = 1,N 
     read(lun,*) psys%src(i)%pos(1), psys%src(i)%pos(2), psys%src(i)%pos(3), &
                 psys%src(i)%vel(1), psys%src(i)%vel(2), psys%src(i)%vel(3), &
                 psys%src(i)%L, &
                 psys%src(i)%SpcType, &
                 psys%src(i)%EmisPrf
  end do

  close(lun)
 
  psys%src(1:N)%lastemit = GV%rayn

  
end subroutine read_src_snapshot


!> Reorders the sources according to the array order. 
!! The array order is not preserved.
!========================================================
subroutine order_sources(srcs,order) 
  type(source_type), intent(inout) :: srcs(:)  !< sources to order
  integer(i8b), intent(inout) :: order(:)         !< desired order
  
  type(source_type) :: dumsrc
  integer(i8b) :: i,goal,nsrc
  
  if (size(srcs) /= size(order)) stop "size(srcs) /= size(order)"
  nsrc = size(srcs)

  do i = 1,nsrc
     dumsrc = srcs(i)
     goal = order(i)
     do while(goal < i)
        goal=order(goal)
        order(i)=goal
     end do
     srcs(i) = srcs(goal)
     srcs(goal)=dumsrc
  end do
  
end subroutine order_sources



!> Orders the sources from brightest to dimmest.  Also
!! sets the cumulative luminosity function entry for each source
!================================================================
subroutine order_sources_lum(src) 
  use m_mrgrnk, only: mrgrnk
  
  type(source_type), intent(inout) :: src(:) !< source array
  
  type(source_type), allocatable :: srclist(:)
  real(r8b), allocatable :: lumlist(:)
  integer(i8b), allocatable :: order(:)
  integer(i8b) :: i,N
  integer(i8b) :: err
  real(r8b) :: Ltot
  
  !      write(*,'(A)') "sorting the sources from brightest to dimmest ... "
  
  Ltot = 0.0
  N = size(src)

  allocate( srclist(N), stat=err )
  if(err/=0) stop "  order_sources_lum> cant allocate srclist"
      
  allocate( lumlist(N) , stat=err)
  if(err/=0) stop "  order_sources_lum> cant allocate lumlist"
    
  allocate( order(N) , stat=err)
  if(err/=0) stop "  order_sources_lum> cant allocate order"

  srclist = src(1:N)
  lumlist = src(1:N)%L
  call mrgrnk(-lumlist(1:N), order)
  
  ! sorting the negative of the luminosities puts the brightest first

  Ltot = 0.0d0
  do i = 1,N
     src(i) = srclist(order(i))
     Ltot = Ltot + src(i)%L
  end do
 
  lumlist = lumlist / Ltot
  
  do i = 1,N
     src(i)%Lcdf = lumlist(order(i))
  end do
  
  do i = 2,N
     src(i)%Lcdf = src(i-1)%Lcdf + src(i)%Lcdf 
  end do
  
  deallocate( srclist, lumlist, order )

end subroutine order_sources_lum


!> forms a snapshot name from a path, a file base, a snapshot number
!! and a file number.
!===================================================================
subroutine form_snapshot_file_name(Path,FileBase,SnapNum,FileNum,SnapFile)
  character(200), intent(in) :: Path       !< path to snapshot dir
  character(200), intent(in) :: FileBase   !< file base names
  integer(i8b), intent(in) :: SnapNum      !< snapshot number
  integer(i8b), intent(in) :: FileNum      !< file number in snapshot
  character(200), intent(out) :: SnapFile  !< file name to return
  
  character(10) :: FileNumChar
  
  write(FileNumChar,"(I6)") FileNum
  
  100 format(A,"/",A,"_",I3.3,".",A)
  write(SnapFile,100) trim(Path), trim(FileBase), SnapNum, &
                      trim(adjustl(FileNumChar))
  
end subroutine form_snapshot_file_name


!> outputs source snapshot header to the screen
!-----------------------------------------------
subroutine source_header_to_screen(shead)
  type(source_header_type), intent(in) :: shead !< source header to print

97 format(T2,"=",T58,"=")
98 format(T2,"= ",A,T58,"=")
99 format(T2,57("="))
100 format(T2,"=",T4,A,T30,I6,T58,"=")
200 format(T2,"=",T4,A,T30,ES12.5,T58,"=")
201 format(T2,"=",T4,A,T58,"=")
  
  write(*,99) 
  write(*,98) "source header data"
  write(*,99) 
  write(*,97) 
  write(*,100) "sources in snapshot:", shead%NsrcSnap
  write(*,100) "sources in file:", shead%NsrcFile
  write(*,100) "number of files in snap:", shead%NumFiles
  write(*,200) "rays for this snap:*", shead%TotalRays
  write(*,201) "*(only if RayScheme='header' in config file)"
  write(*,200) "luminosity unit (photons/s):", shead%Lunit
  write(*,97) 
  write(*,99) 
  
end subroutine source_header_to_screen


!> outputs source snapshot header to file
!-----------------------------------------------
subroutine source_header_to_file(shead,lun)
  type(source_header_type), intent(in) :: shead !< source header to print
  integer(i8b), intent(in) :: lun

97 format(T2,"=",T58,"=")
98 format(T2,"= ",A,T58,"=")
99 format(T2,57("="))
100 format(T2,"=",T4,A,T30,I6,T58,"=")
200 format(T2,"=",T4,A,T30,ES12.5,T58,"=")
201 format(T2,"=",T4,A,T58,"=")
  
  write(lun,99) 
  write(lun,98) "source header data"
  write(lun,99) 
  write(lun,97) 
  write(lun,100) "sources in snapshot:", shead%NsrcSnap
  write(lun,100) "sources in file:", shead%NsrcFile
  write(lun,100) "number of files in snap:", shead%NumFiles
  write(lun,200) "rays for this snap:*", shead%TotalRays
  write(lun,201) "*(only if RayScheme='header' in config file)"
  write(lun,200) "luminosity unit (photons/s):", shead%Lunit
  write(lun,97) 
  write(lun,99) 
  
end subroutine source_header_to_file

!> outputs the source data currently loaded in the particle system to screen
!! -------------------------------------------------------------------------
subroutine source_info_to_screen(src)

  integer(i8b), parameter :: srclimit = 20     !< max number of sources to screen
  type(source_type), intent(in) :: src(:) !< get sources from here
  integer(i8b) :: i
  
  100 format(72("="))
  110 format(T2,A,3ES10.3)
  111 format(T2,A,ES10.3)
  
  write(*,100)
  write(*,*) " source data  "
  write(*,100)
  write(*,*) 
  write(*,'(A,I3,A)') "writing first ", srclimit, &
                      " sources to screen (code units) ..."
  
  write(*,*) 
  120 format(T1,A, T15,A,     T37,A,     T51,A,     T65,A)
  121 format(T1,I2,T5,3ES10.1,T37,ES10.3,T51,ES10.3,T65,I3)
  write(*,120) "Src","Position","Luminosity","Spectrum","Emis Prf"
  write(*,100)
  
  do i = 1, size(src)
     write(*,121) i, src(i)%pos, src(i)%L, src(i)%SpcType, &
          src(i)%EmisPrf
     if (i==srclimit) exit
  end do
  
  write(*,100) 
  write(*,*)
  
  
end subroutine source_info_to_screen

end module source_input_mod
