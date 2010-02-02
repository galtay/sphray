module ion_table_class
use myf90_mod

#ifdef hdf5
use hdf5_wrapper
#endif

#ifdef usempi
use mpi
#endif

implicit none
private

 
public :: ion_table_type
public :: ion_header_type
public :: ion_spectrum_type
public :: read_ion_table_file
public :: interpolate_ion_table
public :: broadcast_ion_table


type ion_spectrum_type
   sequence
   integer :: s_gammaHI
   integer :: s_z
   integer :: s_logryd
   integer :: s_logflux
   real, allocatable :: gammaHI(:)
   real, allocatable :: z(:)
   real, allocatable :: logryd(:)    ! [Rydbergs]
   real, allocatable :: logflux(:,:) ! [ergs/s/Hz/cm^2/sr]
   character(clen) :: model_name
end type ion_spectrum_type

type ion_header_type
   sequence
   character(clen) :: cloudy_version
   type(ion_spectrum_type) :: ispec
end type ion_header_type

type ion_table_type
   sequence
   type(ion_header_type) :: ihead
   integer :: s_ibal
   integer :: s_z
   integer :: s_logt
   integer :: s_logd
   real :: logd_max
   real :: logt_max
   real, allocatable :: ibal(:,:,:)  ! ibal( iz, it, id )
   real, allocatable :: z(:)
   real, allocatable :: logt(:)
   real, allocatable :: logd(:)
end type ion_table_type





contains

  subroutine broadcast_ion_table( MyPE, itab )
    integer, intent(in) :: MyPE
    type(ion_table_type) :: itab
    integer :: count
    integer :: root
    integer :: ierr
    
    root = 0

#ifdef usempi

    count = 1
    call mpi_bcast( itab%logd_max, count, mpi_real, root, mpi_comm_world, ierr )
    call mpi_bcast( itab%logt_max, count, mpi_real, root, mpi_comm_world, ierr )

    call mpi_bcast( itab%s_logd,    count, mpi_integer, root, mpi_comm_world, ierr )
    call mpi_bcast( itab%s_logt,    count, mpi_integer, root, mpi_comm_world, ierr )
    call mpi_bcast( itab%s_z,       count, mpi_integer, root, mpi_comm_world, ierr )
    call mpi_bcast( itab%s_ibal,    count, mpi_integer, root, mpi_comm_world, ierr )

    call mpi_bcast( itab%ihead%ispec%s_logflux, count, mpi_integer, root, mpi_comm_world, ierr )
    call mpi_bcast( itab%ihead%ispec%s_logryd,  count, mpi_integer, root, mpi_comm_world, ierr )
    call mpi_bcast( itab%ihead%ispec%s_z,       count, mpi_integer, root, mpi_comm_world, ierr )
    call mpi_bcast( itab%ihead%ispec%s_gammaHI, count, mpi_integer, root, mpi_comm_world, ierr )


    call mpi_barrier(mpi_comm_world, ierr)
    if (MyPE /= root) then
       allocate( itab%logd( itab%s_logd ) )
       allocate( itab%logt( itab%s_logt ) )
       allocate( itab%z   ( itab%s_z    ) )
       allocate( itab%ibal( itab%s_z, itab%s_logt, itab%s_logd ) )

       allocate( itab%ihead%ispec%logflux ( itab%ihead%ispec%s_z, itab%ihead%ispec%s_logryd  ) )
       allocate( itab%ihead%ispec%logryd  ( itab%ihead%ispec%s_logryd  ) )
       allocate( itab%ihead%ispec%z       ( itab%ihead%ispec%s_z       ) )
       allocate( itab%ihead%ispec%gammaHI ( itab%ihead%ispec%s_gammaHI ) )
    endif
    call mpi_barrier(mpi_comm_world, ierr)


    count = itab%s_logd  
    call mpi_bcast( itab%logd, count, mpi_real, root, mpi_comm_world, ierr )

    count = itab%s_logt  
    call mpi_bcast( itab%logt, count, mpi_real, root, mpi_comm_world, ierr )

    count = itab%s_z 
    call mpi_bcast( itab%z, count, mpi_real, root, mpi_comm_world, ierr )

    count = itab%s_ibal 
    call mpi_bcast( itab%ibal, count, mpi_real, root, mpi_comm_world, ierr )


    count = clen
    call mpi_bcast( itab%ihead%cloudy_version,   count, mpi_character, root, mpi_comm_world, ierr )
    call mpi_bcast( itab%ihead%ispec%model_name, count, mpi_character, root, mpi_comm_world, ierr )

    count = itab%ihead%ispec%s_logflux
    call mpi_bcast( itab%ihead%ispec%logflux(1,1), count, mpi_real, root, mpi_comm_world, ierr )

    count = itab%ihead%ispec%s_logryd 
    call mpi_bcast( itab%ihead%ispec%logryd, count, mpi_real, root, mpi_comm_world, ierr )

    count = itab%ihead%ispec%s_z 
    call mpi_bcast( itab%ihead%ispec%z, count, mpi_real, root, mpi_comm_world, ierr )

    count = itab%ihead%ispec%s_gammaHI 
    call mpi_bcast( itab%ihead%ispec%gammaHI, count, mpi_real, root, mpi_comm_world, ierr )



#endif


  end subroutine broadcast_ion_table



  subroutine read_ion_table_file( file, itab )
    character(*) :: file
    type(ion_table_type) :: itab
    integer :: fh

#ifdef hdf5    
    integer :: rank
    integer :: dims(3)
    character(clen) :: fmt

    integer :: nz, nt, nd

    fmt = "(' ',A, 4I8)"
    write(*,*)

    call hdf5_open_file( fh, file, readonly=.true. )

    call hdf5_read_attribute(fh, 'header/cloudy_version', itab%ihead%cloudy_version)
    call hdf5_read_attribute(fh, 'header/spectrum/model_name', itab%ihead%ispec%model_name)

    write(*,*) 
    write(*,*) "cloudy version = ", trim(itab%ihead%cloudy_version)
    write(*,*) "model name     = ", trim(itab%ihead%ispec%model_name)
    write(*,*) 

    dims=0
    call hdf5_get_dimensions(fh, 'header/spectrum/gammahi', rank, dims)
    allocate( itab%ihead%ispec%gammaHI( dims(1) ) )
    call hdf5_read_data(fh, 'header/spectrum/gammahi', itab%ihead%ispec%gammaHI )
    write(*,fmt) "rank and dims of gammaHI:       ", rank, dims
    itab%ihead%ispec%s_gammaHI = dims(1) 

    dims=0
    call hdf5_get_dimensions(fh, 'header/spectrum/logenergy_ryd', rank, dims)
    allocate( itab%ihead%ispec%logryd( dims(1) ) )
    call hdf5_read_data(fh, 'header/spectrum/logenergy_ryd', itab%ihead%ispec%logryd )
    write(*,fmt) "rank and dims of logenergy_ryd: ", rank, dims
    itab%ihead%ispec%s_logryd = dims(1)

    dims=0
    call hdf5_get_dimensions(fh, 'header/spectrum/logflux', rank, dims)
    allocate( itab%ihead%ispec%logflux( dims(1), dims(2) ) )
    call hdf5_read_data(fh, 'header/spectrum/logflux', itab%ihead%ispec%logflux )
    write(*,fmt) "rank and dims of logflux:       ", rank, dims
    itab%ihead%ispec%s_logflux = product( dims(1:2) )

    dims=0
    call hdf5_get_dimensions(fh, 'header/spectrum/redshift', rank, dims)
    allocate( itab%ihead%ispec%z( dims(1) ) )
    call hdf5_read_data(fh, 'header/spectrum/redshift', itab%ihead%ispec%z )
    write(*,fmt) "rank and dims of redshift:      ", rank, dims
    itab%ihead%ispec%s_z = dims(1)

    dims=0
    call hdf5_get_dimensions(fh, 'ionbal', rank, dims)
    allocate( itab%ibal( dims(1), dims(2), dims(3) ) )
    call hdf5_read_data(fh, 'ionbal', itab%ibal )
    write(*,fmt) "rank and dims of ionbal:        ", rank, dims
    itab%s_ibal = product( dims(1:3) )

    dims=0
    call hdf5_get_dimensions(fh, 'logd', rank, dims)
    allocate( itab%logd( dims(1) ) )
    call hdf5_read_data(fh, 'logd', itab%logd )
    write(*,fmt) "rank and dims of logd:          ", rank, dims
    nd = dims(1)
    itab%s_logd = dims(1) 

    dims=0
    call hdf5_get_dimensions(fh, 'logt', rank, dims)
    allocate( itab%logt( dims(1) ) )
    call hdf5_read_data(fh, 'logt', itab%logt )
    write(*,fmt) "rank and dims of logt:          ", rank, dims
    nt = dims(1)
    itab%s_logt = dims(1) 

    dims=0
    call hdf5_get_dimensions(fh, 'redshift', rank, dims)
    allocate( itab%z( dims(1) ) )
    call hdf5_read_data(fh, 'redshift', itab%z )
    write(*,fmt) "rank and dims of z:             ", rank, dims
    nz = dims(1)
    itab%s_z = dims(1) 

    call hdf5_close_file( fh )

    itab%ibal = log10( itab%ibal )

    itab%logd_max = itab%logd(nd)
    itab%logt_max = itab%logt(nt)


    write(*,*) 
    write(*,*) "spectrum:"
    write(*,*) "min/max G_HI    = ", minval(itab%ihead%ispec%gammaHI), maxval(itab%ihead%ispec%gammaHI)
    write(*,*) "min/max logryd  = ", minval(itab%ihead%ispec%logryd), maxval(itab%ihead%ispec%logryd)
    write(*,*) "min/max logflux = ", minval(itab%ihead%ispec%logflux), maxval(itab%ihead%ispec%logflux)
    write(*,*) "min/max z       = ", minval(itab%ihead%ispec%z), maxval(itab%ihead%ispec%z)
    write(*,*) 
    write(*,*) "ionization balance:"
    write(*,*) "min/max log10(ionbal)  = ", minval(itab%ibal), maxval(itab%ibal)
    write(*,*) "min/max logd           = ", minval(itab%logd), maxval(itab%logd)
    write(*,*) "min/max logt           = ", minval(itab%logt), maxval(itab%logt)
    write(*,*) "min/max z              = ", minval(itab%z), maxval(itab%z)
    write(*,*) 

#else

    stop "trying to read hdf5 ion table, but haven't defined hdf5 macro in Makefile"

#endif

  end subroutine read_ion_table_file




  function interpolate_ion_table( itab, z, logt, logd ) result (ioneq)
    type(ion_table_type) :: itab
    real :: z      !< input redshift
    real :: logt   !< input log temperature
    real :: logd   !< input log density
    real :: ioneq  !< returned ionization equilibrium

    integer :: nz, nt, nd !< size of interpolation grid
    integer :: iz, it, id !< lower bracketing index
    integer :: fz, ft, fd !< upper bracketing index
    real    :: rz, rt, rd !< real index (between bracketing index)

    real :: wiz, wfz, wit, wft, wid, wfd  !< fractional distances between bracketing indices
    real :: w111, w211, w121, w221, w112, w212, w122, w222  !< weights
    real :: dz, dt, dd                    !< difference between bracketing indices

    integer :: i !< general loop index

    ! get sizes of interpolation grid
    !----------------------------------
    nz = size(itab%z)
    nt = size(itab%logt)
    nd = size(itab%logd)

    ! do some checks
    !-------------------    
    if ( z < itab%z(1) )  then
       write(*,*) "warning z below limit in table"
       stop
    endif

    if ( z > itab%z(nz) ) then
       write(*,*) "warning z above limit in table"
       stop
    endif

    if ( logt < itab%logt(1) )  then
!       write(*,*) "warning logt below limit in table"
       logt = itab%logt(1)
    endif

    if ( logt > itab%logt(nt) ) then
!       write(*,*) "warning logt above limit in table"
       logt = itab%logt(nt)
    endif

    if ( logd < itab%logd(1) )  then
!       write(*,*) "warning logd below limit in table"
       logd = itab%logd(1)
    endif

    if ( logd > itab%logd(nd) ) then
!       write(*,*) "warning logd above limit in table"
       logd = itab%logd(nd)
    endif

    
    

    ! find bracketing indices
    !--------------------------
    do i = 1, nz
       if (itab%z(i) > z) then
          fz = i
          iz = fz-1
          exit
       end if
    end do

    do i = 1, nt
       if (itab%logt(i) > logt) then
          ft = i
          it = ft-1
          exit
       end if
    end do

    do i = 1, nd
       if (itab%logd(i) > logd) then
          fd = i
          id = fd-1
          exit
       end if
    end do


    if (z==itab%z(1)) then
       iz = 1
       fz = 2
    else if (z==itab%z(nz)) then
       fz = nz
       iz = nz-1
    endif

    if (logt==itab%logt(1)) then
       it = 1
       ft = 2
    else if (logt==itab%logt(nt)) then
       ft = nt
       it = nt-1
    endif

    if (logd==itab%logd(1)) then
       id = 1
       fd = 2
    else if (logd==itab%logd(nd)) then
       fd = nd
       id = nd-1
    endif


    dz = itab%z(fz) - itab%z(iz)
    wiz = ( z - itab%z(iz) ) / dz
    wfz = ( itab%z(fz) - z ) / dz

    dt = itab%logt(ft) - itab%logt(it)
    wit = ( logt - itab%logt(it) ) / dt
    wft = ( itab%logt(ft) - logt ) / dt

    dd = itab%logd(fd) - itab%logd(id)
    wid = ( logd - itab%logd(id) ) / dd
    wfd = ( itab%logd(fd) - logd ) / dd

    rz = iz + wiz
    rt = it + wit
    rd = id + wid

 
    w111 = wfz * wft * wfd
    w211 = wiz * wft * wfd
    w121 = wfz * wit * wfd
    w221 = wiz * wit * wfd
    w112 = wfz * wft * wid
    w212 = wiz * wft * wid
    w122 = wfz * wit * wid
    w222 = wiz * wit * wid


    ioneq = &
         w111 * itab%ibal( iz, it, id ) + &
         w211 * itab%ibal( fz, it, id ) + &
         w121 * itab%ibal( iz, ft, id ) + &
         w221 * itab%ibal( fz, ft, id ) + &
         w112 * itab%ibal( iz, it, fd ) + &
         w212 * itab%ibal( fz, it, fd ) + &
         w122 * itab%ibal( iz, ft, fd ) + &
         w222 * itab%ibal( fz, ft, fd ) 

    
    ioneq = 10**ioneq



   
  end function interpolate_ion_table




end module ion_table_class
