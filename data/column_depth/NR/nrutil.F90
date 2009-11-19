module nrutil
  use nrtype
  implicit none
  integer(I4B), parameter :: NPAR_ARTH=16,NPAR2_ARTH=8
  integer(I4B), parameter :: NPAR_GEOP=4,NPAR2_GEOP=2
  integer(I4B), parameter :: NPAR_CUMSUM=16
  integer(I4B), parameter :: NPAR_CUMPROD=8
  integer(I4B), parameter :: NPAR_POLY=8
  integer(I4B), parameter :: NPAR_POLYTERM=8

  interface array_copy
     module procedure array_copy_r, array_copy_d, array_copy_i
  end interface
  interface swap
     module procedure swap_i,swap_r,swap_rv,swap_c, &
          swap_cv,swap_cm,swap_z,swap_zv,swap_zm, &
          masked_swap_rs,masked_swap_rv,masked_swap_rm
  end interface
  interface reallocate
     module procedure reallocate_rv,reallocate_rm,&
          reallocate_iv,reallocate_im,reallocate_hv
  end interface
  interface imaxloc
     module procedure imaxloc_r,imaxloc_i
  end interface
  interface iminloc
     module procedure iminloc_sp,iminloc_dp
  end interface
  interface assert
     module procedure assert1,assert2,assert3,assert4,assert_v
  end interface
  interface assert_eq
     module procedure assert_eq2,assert_eq3,assert_eq4,assert_eqn
  end interface
  interface arth
     module procedure arth_r, arth_d, arth_i
  end interface
  interface geop
     module procedure geop_r, geop_d, geop_i, geop_c, geop_dv
  end interface
  interface cumsum
     module procedure cumsum_r,cumsum_i
  end interface
  interface poly
     module procedure poly_rr,poly_rrv,poly_dd,poly_ddv,&
          poly_rc,poly_cc,poly_msk_rrv,poly_msk_ddv
  end interface
  interface poly_term
     module procedure poly_term_rr,poly_term_cc
  end interface
  interface outerprod
     module procedure outerprod_r,outerprod_d
  end interface
  interface outerdiff
     module procedure outerdiff_r,outerdiff_d,outerdiff_i
  end interface
  interface scatter_add
     module procedure scatter_add_r,scatter_add_d
  end interface
  interface scatter_max
     module procedure scatter_max_r,scatter_max_d
  end interface
  interface diagadd
     module procedure diagadd_rv,diagadd_r
  end interface
  interface diagmult
     module procedure diagmult_rv,diagmult_r
  end interface
  interface get_diag
     module procedure get_diag_rv, get_diag_dv
  end interface
  interface put_diag
     module procedure put_diag_rv, put_diag_r
  end interface
contains
  !BL
  subroutine array_copy_r(src,dest,n_copied,n_not_copied)
    real(SP), dimension(:), intent(IN) :: src
    real(SP), dimension(:), intent(OUT) :: dest
    integer(I4B), intent(OUT) :: n_copied, n_not_copied
    n_copied=min(size(src),size(dest))
    n_not_copied=size(src)-n_copied
    dest(1:n_copied)=src(1:n_copied)
  end subroutine array_copy_r
  !BL
  subroutine array_copy_d(src,dest,n_copied,n_not_copied)
    real(DP), dimension(:), intent(IN) :: src
    real(DP), dimension(:), intent(OUT) :: dest
    integer(I4B), intent(OUT) :: n_copied, n_not_copied
    n_copied=min(size(src),size(dest))
    n_not_copied=size(src)-n_copied
    dest(1:n_copied)=src(1:n_copied)
  end subroutine array_copy_d
  !BL
  subroutine array_copy_i(src,dest,n_copied,n_not_copied)
    integer(I4B), dimension(:), intent(IN) :: src
    integer(I4B), dimension(:), intent(OUT) :: dest
    integer(I4B), intent(OUT) :: n_copied, n_not_copied
    n_copied=min(size(src),size(dest))
    n_not_copied=size(src)-n_copied
    dest(1:n_copied)=src(1:n_copied)
  end subroutine array_copy_i
  !BL
  !BL
  subroutine swap_i(a,b)
    integer(I4B), intent(INOUT) :: a,b
    integer(I4B) :: dum
    dum=a
    a=b
    b=dum
  end subroutine swap_i
  !BL
  subroutine swap_r(a,b)
    real(SP), intent(INOUT) :: a,b
    real(SP) :: dum
    dum=a
    a=b
    b=dum
  end subroutine swap_r
  !BL
  subroutine swap_rv(a,b)
    real(SP), dimension(:), intent(INOUT) :: a,b
    real(SP), dimension(size(a)) :: dum
    dum=a
    a=b
    b=dum
  end subroutine swap_rv
  !BL
  subroutine swap_c(a,b)
    complex(SPC), intent(INOUT) :: a,b
    complex(SPC) :: dum
    dum=a
    a=b
    b=dum
  end subroutine swap_c
  !BL
  subroutine swap_cv(a,b)
    complex(SPC), dimension(:), intent(INOUT) :: a,b
    complex(SPC), dimension(size(a)) :: dum
    dum=a
    a=b
    b=dum
  end subroutine swap_cv
  !BL
  subroutine swap_cm(a,b)
    complex(SPC), dimension(:,:), intent(INOUT) :: a,b
    complex(SPC), dimension(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
  end subroutine swap_cm
  !BL
  subroutine swap_z(a,b)
    complex(DPC), intent(INOUT) :: a,b
    complex(DPC) :: dum
    dum=a
    a=b
    b=dum
  end subroutine swap_z
  !BL
  subroutine swap_zv(a,b)
    complex(DPC), dimension(:), intent(INOUT) :: a,b
    complex(DPC), dimension(size(a)) :: dum
    dum=a
    a=b
    b=dum
  end subroutine swap_zv
  !BL
  subroutine swap_zm(a,b)
    complex(DPC), dimension(:,:), intent(INOUT) :: a,b
    complex(DPC), dimension(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
  end subroutine swap_zm
  !BL
  subroutine masked_swap_rs(a,b,mask)
    real(SP), intent(INOUT) :: a,b
    logical(LGT), intent(IN) :: mask
    real(SP) :: swp
    if (mask) then
       swp=a
       a=b
       b=swp
    end if
  end subroutine masked_swap_rs
  !BL
  subroutine masked_swap_rv(a,b,mask)
    real(SP), dimension(:), intent(INOUT) :: a,b
    logical(LGT), dimension(:), intent(IN) :: mask
    real(SP), dimension(size(a)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  end subroutine masked_swap_rv
  !BL
  subroutine masked_swap_rm(a,b,mask)
    real(SP), dimension(:,:), intent(INOUT) :: a,b
    logical(LGT), dimension(:,:), intent(IN) :: mask
    real(SP), dimension(size(a,1),size(a,2)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  end subroutine masked_swap_rm
  !BL
  !BL
  function reallocate_rv(p,n)
    real(SP), dimension(:), pointer :: p, reallocate_rv
    integer(I4B), intent(IN) :: n
    integer(I4B) :: nold,ierr
    allocate(reallocate_rv(n),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_rv: problem in attempt to allocate memory')
    if (.not. associated(p)) return
    nold=size(p)
    reallocate_rv(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  end function reallocate_rv
  !BL
  function reallocate_iv(p,n)
    integer(I4B), dimension(:), pointer :: p, reallocate_iv
    integer(I4B), intent(IN) :: n
    integer(I4B) :: nold,ierr
    allocate(reallocate_iv(n),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_iv: problem in attempt to allocate memory')
    if (.not. associated(p)) return
    nold=size(p)
    reallocate_iv(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  end function reallocate_iv
  !BL
  function reallocate_hv(p,n)
    character(1), dimension(:), pointer :: p, reallocate_hv
    integer(I4B), intent(IN) :: n
    integer(I4B) :: nold,ierr
    allocate(reallocate_hv(n),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_hv: problem in attempt to allocate memory')
    if (.not. associated(p)) return
    nold=size(p)
    reallocate_hv(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  end function reallocate_hv
  !BL
  function reallocate_rm(p,n,m)
    real(SP), dimension(:,:), pointer :: p, reallocate_rm
    integer(I4B), intent(IN) :: n,m
    integer(I4B) :: nold,mold,ierr
    allocate(reallocate_rm(n,m),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_rm: problem in attempt to allocate memory')
    if (.not. associated(p)) return
    nold=size(p,1)
    mold=size(p,2)
    reallocate_rm(1:min(nold,n),1:min(mold,m))=&
         p(1:min(nold,n),1:min(mold,m))
    deallocate(p)
  end function reallocate_rm
  !BL
  function reallocate_im(p,n,m)
    integer(I4B), dimension(:,:), pointer :: p, reallocate_im
    integer(I4B), intent(IN) :: n,m
    integer(I4B) :: nold,mold,ierr
    allocate(reallocate_im(n,m),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_im: problem in attempt to allocate memory')
    if (.not. associated(p)) return
    nold=size(p,1)
    mold=size(p,2)
    reallocate_im(1:min(nold,n),1:min(mold,m))=&
         p(1:min(nold,n),1:min(mold,m))
    deallocate(p)
  end function reallocate_im
  !BL
  function ifirstloc(mask)
    logical(LGT), dimension(:), intent(IN) :: mask
    integer(I4B) :: ifirstloc
    integer(I4B), dimension(1) :: loc
    loc=maxloc(merge(1,0,mask))
    ifirstloc=loc(1)
    if (.not. mask(ifirstloc)) ifirstloc=size(mask)+1
  end function ifirstloc
  !BL
  function imaxloc_r(arr)
    real(SP), dimension(:), intent(IN) :: arr
    integer(I4B) :: imaxloc_r
    integer(I4B), dimension(1) :: imax
    imax=maxloc(arr(:))
    imaxloc_r=imax(1)
  end function imaxloc_r
  !BL
  function imaxloc_i(iarr)
    integer(I4B), dimension(:), intent(IN) :: iarr
    integer(I4B), dimension(1) :: imax
    integer(I4B) :: imaxloc_i
    imax=maxloc(iarr(:))
    imaxloc_i=imax(1)
  end function imaxloc_i
  !BL
  function iminloc_sp(arr)
    real(SP), dimension(:), intent(IN) :: arr
    integer(I4B), dimension(1) :: imin
    integer(I4B) :: iminloc_sp
    imin=minloc(arr(:))
    iminloc_sp=imin(1)
  end function iminloc_sp
  !BL
  function iminloc_dp(arr)
    real(DP), dimension(:), intent(IN) :: arr
    integer(I4B), dimension(1) :: imin
    integer(I4B) :: iminloc_dp
    imin=minloc(arr(:))
    iminloc_dp=imin(1)
  end function iminloc_dp
  !BL
  subroutine assert1(n1,string)
    character(LEN=*), intent(IN) :: string
    logical, intent(IN) :: n1
    if (.not. n1) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       stop 'program terminated by assert1'
    end if
  end subroutine assert1
  !BL
  subroutine assert2(n1,n2,string)
    character(LEN=*), intent(IN) :: string
    logical, intent(IN) :: n1,n2
    if (.not. (n1 .and. n2)) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       stop 'program terminated by assert2'
    end if
  end subroutine assert2
  !BL
  subroutine assert3(n1,n2,n3,string)
    character(LEN=*), intent(IN) :: string
    logical, intent(IN) :: n1,n2,n3
    if (.not. (n1 .and. n2 .and. n3)) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       stop 'program terminated by assert3'
    end if
  end subroutine assert3
  !BL
  subroutine assert4(n1,n2,n3,n4,string)
    character(LEN=*), intent(IN) :: string
    logical, intent(IN) :: n1,n2,n3,n4
    if (.not. (n1 .and. n2 .and. n3 .and. n4)) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       stop 'program terminated by assert4'
    end if
  end subroutine assert4
  !BL
  subroutine assert_v(n,string)
    character(LEN=*), intent(IN) :: string
    logical, dimension(:), intent(IN) :: n
    if (.not. all(n)) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       stop 'program terminated by assert_v'
    end if
  end subroutine assert_v
  !BL
  function assert_eq2(n1,n2,string)
    character(LEN=*), intent(IN) :: string
    integer, intent(IN) :: n1,n2
    integer :: assert_eq2
    if (n1 == n2) then
       assert_eq2=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       stop 'program terminated by assert_eq2'
    end if
  end function assert_eq2
  !BL
  function assert_eq3(n1,n2,n3,string)
    character(LEN=*), intent(IN) :: string
    integer, intent(IN) :: n1,n2,n3
    integer :: assert_eq3
    if (n1 == n2 .and. n2 == n3) then
       assert_eq3=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       stop 'program terminated by assert_eq3'
    end if
  end function assert_eq3
  !BL
  function assert_eq4(n1,n2,n3,n4,string)
    character(LEN=*), intent(IN) :: string
    integer, intent(IN) :: n1,n2,n3,n4
    integer :: assert_eq4
    if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
       assert_eq4=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       stop 'program terminated by assert_eq4'
    end if
  end function assert_eq4
  !BL
  function assert_eqn(nn,string)
    character(LEN=*), intent(IN) :: string
    integer, dimension(:), intent(IN) :: nn
    integer :: assert_eqn
    if (all(nn(2:) == nn(1))) then
       assert_eqn=nn(1)
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       stop 'program terminated by assert_eqn'
    end if
  end function assert_eqn
  !BL
  subroutine nrerror(string)
    character(LEN=*), intent(IN) :: string
    write (*,*) 'nrerror: ',string
    stop 'program terminated by nrerror'
  end subroutine nrerror
  !BL
  function arth_r(first,increment,n)
    real(SP), intent(IN) :: first,increment
    integer(I4B), intent(IN) :: n
    real(SP), dimension(n) :: arth_r
    integer(I4B) :: k,k2
    real(SP) :: temp
    if (n > 0) arth_r(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_r(k)=arth_r(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_r(k)=arth_r(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  end function arth_r
  !BL
  function arth_d(first,increment,n)
    real(DP), intent(IN) :: first,increment
    integer(I4B), intent(IN) :: n
    real(DP), dimension(n) :: arth_d
    integer(I4B) :: k,k2
    real(DP) :: temp
    if (n > 0) arth_d(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_d(k)=arth_d(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_d(k)=arth_d(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  end function arth_d
  !BL
  function arth_i(first,increment,n)
    integer(I4B), intent(IN) :: first,increment,n
    integer(I4B), dimension(n) :: arth_i
    integer(I4B) :: k,k2,temp
    if (n > 0) arth_i(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_i(k)=arth_i(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_i(k)=arth_i(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  end function arth_i
  !BL
  !BL
  function geop_r(first,factor,n)
    real(SP), intent(IN) :: first,factor
    integer(I4B), intent(IN) :: n
    real(SP), dimension(n) :: geop_r
    integer(I4B) :: k,k2
    real(SP) :: temp
    if (n > 0) geop_r(1)=first
    if (n <= NPAR_GEOP) then
       do k=2,n
          geop_r(k)=geop_r(k-1)*factor
       end do
    else
       do k=2,NPAR2_GEOP
          geop_r(k)=geop_r(k-1)*factor
       end do
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       do
          if (k >= n) exit
          k2=k+k
          geop_r(k+1:min(k2,n))=temp*geop_r(1:min(k,n-k))
          temp=temp*temp
          k=k2
       end do
    end if
  end function geop_r
  !BL
  function geop_d(first,factor,n)
    real(DP), intent(IN) :: first,factor
    integer(I4B), intent(IN) :: n
    real(DP), dimension(n) :: geop_d
    integer(I4B) :: k,k2
    real(DP) :: temp
    if (n > 0) geop_d(1)=first
    if (n <= NPAR_GEOP) then
       do k=2,n
          geop_d(k)=geop_d(k-1)*factor
       end do
    else
       do k=2,NPAR2_GEOP
          geop_d(k)=geop_d(k-1)*factor
       end do
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       do
          if (k >= n) exit
          k2=k+k
          geop_d(k+1:min(k2,n))=temp*geop_d(1:min(k,n-k))
          temp=temp*temp
          k=k2
       end do
    end if
  end function geop_d
  !BL
  function geop_i(first,factor,n)
    integer(I4B), intent(IN) :: first,factor,n
    integer(I4B), dimension(n) :: geop_i
    integer(I4B) :: k,k2,temp
    if (n > 0) geop_i(1)=first
    if (n <= NPAR_GEOP) then
       do k=2,n
          geop_i(k)=geop_i(k-1)*factor
       end do
    else
       do k=2,NPAR2_GEOP
          geop_i(k)=geop_i(k-1)*factor
       end do
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       do
          if (k >= n) exit
          k2=k+k
          geop_i(k+1:min(k2,n))=temp*geop_i(1:min(k,n-k))
          temp=temp*temp
          k=k2
       end do
    end if
  end function geop_i
  !BL
  function geop_c(first,factor,n)
    complex(SP), intent(IN) :: first,factor
    integer(I4B), intent(IN) :: n
    complex(SP), dimension(n) :: geop_c
    integer(I4B) :: k,k2
    complex(SP) :: temp
    if (n > 0) geop_c(1)=first
    if (n <= NPAR_GEOP) then
       do k=2,n
          geop_c(k)=geop_c(k-1)*factor
       end do
    else
       do k=2,NPAR2_GEOP
          geop_c(k)=geop_c(k-1)*factor
       end do
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       do
          if (k >= n) exit
          k2=k+k
          geop_c(k+1:min(k2,n))=temp*geop_c(1:min(k,n-k))
          temp=temp*temp
          k=k2
       end do
    end if
  end function geop_c
  !BL
  function geop_dv(first,factor,n)
    real(DP), dimension(:), intent(IN) :: first,factor
    integer(I4B), intent(IN) :: n
    real(DP), dimension(size(first),n) :: geop_dv
    integer(I4B) :: k,k2
    real(DP), dimension(size(first)) :: temp
    if (n > 0) geop_dv(:,1)=first(:)
    if (n <= NPAR_GEOP) then
       do k=2,n
          geop_dv(:,k)=geop_dv(:,k-1)*factor(:)
       end do
    else
       do k=2,NPAR2_GEOP
          geop_dv(:,k)=geop_dv(:,k-1)*factor(:)
       end do
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       do
          if (k >= n) exit
          k2=k+k
          geop_dv(:,k+1:min(k2,n))=geop_dv(:,1:min(k,n-k))*&
               spread(temp,2,size(geop_dv(:,1:min(k,n-k)),2))
          temp=temp*temp
          k=k2
       end do
    end if
  end function geop_dv
  !BL
  !BL
  recursive function cumsum_r(arr,seed) result(ans)
    real(SP), dimension(:), intent(IN) :: arr
    real(SP), optional, intent(IN) :: seed
    real(SP), dimension(size(arr)) :: ans
    integer(I4B) :: n,j
    real(SP) :: sd
    n=size(arr)
    if (n == 0_i4b) return
    sd=0.0_sp
    if (present(seed)) sd=seed
    ans(1)=arr(1)+sd
    if (n < NPAR_CUMSUM) then
       do j=2,n
          ans(j)=ans(j-1)+arr(j)
       end do
    else
       ans(2:n:2)=cumsum_r(arr(2:n:2)+arr(1:n-1:2),sd)
       ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
    end if
  end function cumsum_r
  !BL
  recursive function cumsum_i(arr,seed) result(ans)
    integer(I4B), dimension(:), intent(IN) :: arr
    integer(I4B), optional, intent(IN) :: seed
    integer(I4B), dimension(size(arr)) :: ans
    integer(I4B) :: n,j,sd
    n=size(arr)
    if (n == 0_i4b) return
    sd=0_i4b
    if (present(seed)) sd=seed
    ans(1)=arr(1)+sd
    if (n < NPAR_CUMSUM) then
       do j=2,n
          ans(j)=ans(j-1)+arr(j)
       end do
    else
       ans(2:n:2)=cumsum_i(arr(2:n:2)+arr(1:n-1:2),sd)
       ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
    end if
  end function cumsum_i
  !BL
  !BL
  recursive function cumprod(arr,seed) result(ans)
    real(SP), dimension(:), intent(IN) :: arr
    real(SP), optional, intent(IN) :: seed
    real(SP), dimension(size(arr)) :: ans
    integer(I4B) :: n,j
    real(SP) :: sd
    n=size(arr)
    if (n == 0_i4b) return
    sd=1.0_sp
    if (present(seed)) sd=seed
    ans(1)=arr(1)*sd
    if (n < NPAR_CUMPROD) then
       do j=2,n
          ans(j)=ans(j-1)*arr(j)
       end do
    else
       ans(2:n:2)=cumprod(arr(2:n:2)*arr(1:n-1:2),sd)
       ans(3:n:2)=ans(2:n-1:2)*arr(3:n:2)
    end if
  end function cumprod
  !BL
  !BL
  function poly_rr(x,coeffs)
    real(SP), intent(IN) :: x
    real(SP), dimension(:), intent(IN) :: coeffs
    real(SP) :: poly_rr
    real(SP) :: pow
    real(SP), dimension(:), allocatable :: vec
    integer(I4B) :: i,n,nn
    n=size(coeffs)
    if (n <= 0) then
       poly_rr=0.0_sp
    else if (n < NPAR_POLY) then
       poly_rr=coeffs(n)
       do i=n-1,1,-1
          poly_rr=x*poly_rr+coeffs(i)
       end do
    else
       allocate(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       do
          vec(n+1)=0.0_sp
          nn=ishft(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          if (nn == 1) exit
          pow=pow*pow
          n=nn
       end do
       poly_rr=vec(1)
       deallocate(vec)
    end if
  end function poly_rr
  !BL
  function poly_dd(x,coeffs)
    real(DP), intent(IN) :: x
    real(DP), dimension(:), intent(IN) :: coeffs
    real(DP) :: poly_dd
    real(DP) :: pow
    real(DP), dimension(:), allocatable :: vec
    integer(I4B) :: i,n,nn
    n=size(coeffs)
    if (n <= 0) then
       poly_dd=0.0_dp
    else if (n < NPAR_POLY) then
       poly_dd=coeffs(n)
       do i=n-1,1,-1
          poly_dd=x*poly_dd+coeffs(i)
       end do
    else
       allocate(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       do
          vec(n+1)=0.0_dp
          nn=ishft(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          if (nn == 1) exit
          pow=pow*pow
          n=nn
       end do
       poly_dd=vec(1)
       deallocate(vec)
    end if
  end function poly_dd
  !BL
  function poly_rc(x,coeffs)
    complex(SPC), intent(IN) :: x
    real(SP), dimension(:), intent(IN) :: coeffs
    complex(SPC) :: poly_rc
    complex(SPC) :: pow
    complex(SPC), dimension(:), allocatable :: vec
    integer(I4B) :: i,n,nn
    n=size(coeffs)
    if (n <= 0) then
       poly_rc=0.0_sp
    else if (n < NPAR_POLY) then
       poly_rc=coeffs(n)
       do i=n-1,1,-1
          poly_rc=x*poly_rc+coeffs(i)
       end do
    else
       allocate(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       do
          vec(n+1)=0.0_sp
          nn=ishft(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          if (nn == 1) exit
          pow=pow*pow
          n=nn
       end do
       poly_rc=vec(1)
       deallocate(vec)
    end if
  end function poly_rc
  !BL
  function poly_cc(x,coeffs)
    complex(SPC), intent(IN) :: x
    complex(SPC), dimension(:), intent(IN) :: coeffs
    complex(SPC) :: poly_cc
    complex(SPC) :: pow
    complex(SPC), dimension(:), allocatable :: vec
    integer(I4B) :: i,n,nn
    n=size(coeffs)
    if (n <= 0) then
       poly_cc=0.0_sp
    else if (n < NPAR_POLY) then
       poly_cc=coeffs(n)
       do i=n-1,1,-1
          poly_cc=x*poly_cc+coeffs(i)
       end do
    else
       allocate(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       do
          vec(n+1)=0.0_sp
          nn=ishft(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          if (nn == 1) exit
          pow=pow*pow
          n=nn
       end do
       poly_cc=vec(1)
       deallocate(vec)
    end if
  end function poly_cc
  !BL
  function poly_rrv(x,coeffs)
    real(SP), dimension(:), intent(IN) :: coeffs,x
    real(SP), dimension(size(x)) :: poly_rrv
    integer(I4B) :: i,n,m
    m=size(coeffs)
    n=size(x)
    if (m <= 0) then
       poly_rrv=0.0_sp
    else if (m < n .or. m < NPAR_POLY) then
       poly_rrv=coeffs(m)
       do i=m-1,1,-1
          poly_rrv=x*poly_rrv+coeffs(i)
       end do
    else
       do i=1,n
          poly_rrv(i)=poly_rr(x(i),coeffs)
       end do
    end if
  end function poly_rrv
  !BL
  function poly_ddv(x,coeffs)
    real(DP), dimension(:), intent(IN) :: coeffs,x
    real(DP), dimension(size(x)) :: poly_ddv
    integer(I4B) :: i,n,m
    m=size(coeffs)
    n=size(x)
    if (m <= 0) then
       poly_ddv=0.0_dp
    else if (m < n .or. m < NPAR_POLY) then
       poly_ddv=coeffs(m)
       do i=m-1,1,-1
          poly_ddv=x*poly_ddv+coeffs(i)
       end do
    else
       do i=1,n
          poly_ddv(i)=poly_dd(x(i),coeffs)
       end do
    end if
  end function poly_ddv
  !BL
  function poly_msk_rrv(x,coeffs,mask)
    real(SP), dimension(:), intent(IN) :: coeffs,x
    logical(LGT), dimension(:), intent(IN) :: mask
    real(SP), dimension(size(x)) :: poly_msk_rrv
    poly_msk_rrv=unpack(poly_rrv(pack(x,mask),coeffs),mask,0.0_sp)
  end function poly_msk_rrv
  !BL
  function poly_msk_ddv(x,coeffs,mask)
    real(DP), dimension(:), intent(IN) :: coeffs,x
    logical(LGT), dimension(:), intent(IN) :: mask
    real(DP), dimension(size(x)) :: poly_msk_ddv
    poly_msk_ddv=unpack(poly_ddv(pack(x,mask),coeffs),mask,0.0_dp)
  end function poly_msk_ddv
  !BL
  !BL
  recursive function poly_term_rr(a,b) result(u)
    real(SP), dimension(:), intent(IN) :: a
    real(SP), intent(IN) :: b
    real(SP), dimension(size(a)) :: u
    integer(I4B) :: n,j
    n=size(a)
    if (n <= 0) return
    u(1)=a(1)
    if (n < NPAR_POLYTERM) then
       do j=2,n
          u(j)=a(j)+b*u(j-1)
       end do
    else
       u(2:n:2)=poly_term_rr(a(2:n:2)+a(1:n-1:2)*b,b*b)
       u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
    end if
  end function poly_term_rr
  !BL
  recursive function poly_term_cc(a,b) result(u)
    complex(SPC), dimension(:), intent(IN) :: a
    complex(SPC), intent(IN) :: b
    complex(SPC), dimension(size(a)) :: u
    integer(I4B) :: n,j
    n=size(a)
    if (n <= 0) return
    u(1)=a(1)
    if (n < NPAR_POLYTERM) then
       do j=2,n
          u(j)=a(j)+b*u(j-1)
       end do
    else
       u(2:n:2)=poly_term_cc(a(2:n:2)+a(1:n-1:2)*b,b*b)
       u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
    end if
  end function poly_term_cc
  !BL
  !BL
  function zroots_unity(n,nn)
    integer(I4B), intent(IN) :: n,nn
    complex(SPC), dimension(nn) :: zroots_unity
    integer(I4B) :: k
    real(SP) :: theta
    zroots_unity(1)=1.0
    theta=TWOPI/n
    k=1
    do
       if (k >= nn) exit
       zroots_unity(k+1)=cmplx(cos(k*theta),sin(k*theta),SPC)
       zroots_unity(k+2:min(2*k,nn))=zroots_unity(k+1)*&
            zroots_unity(2:min(k,nn-k))
       k=2*k
    end do
  end function zroots_unity
  !BL
  function outerprod_r(a,b)
    real(SP), dimension(:), intent(IN) :: a,b
    real(SP), dimension(size(a),size(b)) :: outerprod_r
    outerprod_r = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  end function outerprod_r
  !BL
  function outerprod_d(a,b)
    real(DP), dimension(:), intent(IN) :: a,b
    real(DP), dimension(size(a),size(b)) :: outerprod_d
    outerprod_d = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  end function outerprod_d
  !BL
  function outerdiv(a,b)
    real(SP), dimension(:), intent(IN) :: a,b
    real(SP), dimension(size(a),size(b)) :: outerdiv
    outerdiv = spread(a,dim=2,ncopies=size(b)) / &
         spread(b,dim=1,ncopies=size(a))
  end function outerdiv
  !BL
  function outersum(a,b)
    real(SP), dimension(:), intent(IN) :: a,b
    real(SP), dimension(size(a),size(b)) :: outersum
    outersum = spread(a,dim=2,ncopies=size(b)) + &
         spread(b,dim=1,ncopies=size(a))
  end function outersum
  !BL
  function outerdiff_r(a,b)
    real(SP), dimension(:), intent(IN) :: a,b
    real(SP), dimension(size(a),size(b)) :: outerdiff_r
    outerdiff_r = spread(a,dim=2,ncopies=size(b)) - &
         spread(b,dim=1,ncopies=size(a))
  end function outerdiff_r
  !BL
  function outerdiff_d(a,b)
    real(DP), dimension(:), intent(IN) :: a,b
    real(DP), dimension(size(a),size(b)) :: outerdiff_d
    outerdiff_d = spread(a,dim=2,ncopies=size(b)) - &
         spread(b,dim=1,ncopies=size(a))
  end function outerdiff_d
  !BL
  function outerdiff_i(a,b)
    integer(I4B), dimension(:), intent(IN) :: a,b
    integer(I4B), dimension(size(a),size(b)) :: outerdiff_i
    outerdiff_i = spread(a,dim=2,ncopies=size(b)) - &
         spread(b,dim=1,ncopies=size(a))
  end function outerdiff_i
  !BL
  function outerand(a,b)
    logical(LGT), dimension(:), intent(IN) :: a,b
    logical(LGT), dimension(size(a),size(b)) :: outerand
    outerand = spread(a,dim=2,ncopies=size(b)) .and. &
         spread(b,dim=1,ncopies=size(a))
  end function outerand
  !BL
  subroutine scatter_add_r(dest,source,dest_index)
    real(SP), dimension(:), intent(OUT) :: dest
    real(SP), dimension(:), intent(IN) :: source
    integer(I4B), dimension(:), intent(IN) :: dest_index
    integer(I4B) :: m,n,j,i
    n=assert_eq2(size(source),size(dest_index),'scatter_add_r')
    m=size(dest)
    do j=1,n
       i=dest_index(j)
       if (i > 0 .and. i <= m) dest(i)=dest(i)+source(j)
    end do
  end subroutine scatter_add_r
  subroutine scatter_add_d(dest,source,dest_index)
    real(DP), dimension(:), intent(OUT) :: dest
    real(DP), dimension(:), intent(IN) :: source
    integer(I4B), dimension(:), intent(IN) :: dest_index
    integer(I4B) :: m,n,j,i
    n=assert_eq2(size(source),size(dest_index),'scatter_add_d')
    m=size(dest)
    do j=1,n
       i=dest_index(j)
       if (i > 0 .and. i <= m) dest(i)=dest(i)+source(j)
    end do
  end subroutine scatter_add_d
  subroutine scatter_max_r(dest,source,dest_index)
    real(SP), dimension(:), intent(OUT) :: dest
    real(SP), dimension(:), intent(IN) :: source
    integer(I4B), dimension(:), intent(IN) :: dest_index
    integer(I4B) :: m,n,j,i
    n=assert_eq2(size(source),size(dest_index),'scatter_max_r')
    m=size(dest)
    do j=1,n
       i=dest_index(j)
       if (i > 0 .and. i <= m) dest(i)=max(dest(i),source(j))
    end do
  end subroutine scatter_max_r
  subroutine scatter_max_d(dest,source,dest_index)
    real(DP), dimension(:), intent(OUT) :: dest
    real(DP), dimension(:), intent(IN) :: source
    integer(I4B), dimension(:), intent(IN) :: dest_index
    integer(I4B) :: m,n,j,i
    n=assert_eq2(size(source),size(dest_index),'scatter_max_d')
    m=size(dest)
    do j=1,n
       i=dest_index(j)
       if (i > 0 .and. i <= m) dest(i)=max(dest(i),source(j))
    end do
  end subroutine scatter_max_d
  !BL
  subroutine diagadd_rv(mat,diag)
    real(SP), dimension(:,:), intent(INOUT) :: mat
    real(SP), dimension(:), intent(IN) :: diag
    integer(I4B) :: j,n
    n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagadd_rv')
    do j=1,n
       mat(j,j)=mat(j,j)+diag(j)
    end do
  end subroutine diagadd_rv
  !BL
  subroutine diagadd_r(mat,diag)
    real(SP), dimension(:,:), intent(INOUT) :: mat
    real(SP), intent(IN) :: diag
    integer(I4B) :: j,n
    n = min(size(mat,1),size(mat,2))
    do j=1,n
       mat(j,j)=mat(j,j)+diag
    end do
  end subroutine diagadd_r
  !BL
  subroutine diagmult_rv(mat,diag)
    real(SP), dimension(:,:), intent(INOUT) :: mat
    real(SP), dimension(:), intent(IN) :: diag
    integer(I4B) :: j,n
    n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagmult_rv')
    do j=1,n
       mat(j,j)=mat(j,j)*diag(j)
    end do
  end subroutine diagmult_rv
  !BL
  subroutine diagmult_r(mat,diag)
    real(SP), dimension(:,:), intent(INOUT) :: mat
    real(SP), intent(IN) :: diag
    integer(I4B) :: j,n
    n = min(size(mat,1),size(mat,2))
    do j=1,n
       mat(j,j)=mat(j,j)*diag
    end do
  end subroutine diagmult_r
  !BL
  function get_diag_rv(mat)
    real(SP), dimension(:,:), intent(IN) :: mat
    real(SP), dimension(size(mat,1)) :: get_diag_rv
    integer(I4B) :: j
    j=assert_eq2(size(mat,1),size(mat,2),'get_diag_rv')
    do j=1,size(mat,1)
       get_diag_rv(j)=mat(j,j)
    end do
  end function get_diag_rv
  !BL
  function get_diag_dv(mat)
    real(DP), dimension(:,:), intent(IN) :: mat
    real(DP), dimension(size(mat,1)) :: get_diag_dv
    integer(I4B) :: j
    j=assert_eq2(size(mat,1),size(mat,2),'get_diag_dv')
    do j=1,size(mat,1)
       get_diag_dv(j)=mat(j,j)
    end do
  end function get_diag_dv
  !BL
  subroutine put_diag_rv(diagv,mat)
    real(SP), dimension(:), intent(IN) :: diagv
    real(SP), dimension(:,:), intent(INOUT) :: mat
    integer(I4B) :: j,n
    n=assert_eq2(size(diagv),min(size(mat,1),size(mat,2)),'put_diag_rv')
    do j=1,n
       mat(j,j)=diagv(j)
    end do
  end subroutine put_diag_rv
  !BL
  subroutine put_diag_r(scal,mat)
    real(SP), intent(IN) :: scal
    real(SP), dimension(:,:), intent(INOUT) :: mat
    integer(I4B) :: j,n
    n = min(size(mat,1),size(mat,2))
    do j=1,n
       mat(j,j)=scal
    end do
  end subroutine put_diag_r
  !BL
  subroutine unit_matrix(mat)
    real(SP), dimension(:,:), intent(OUT) :: mat
    integer(I4B) :: i,n
    n=min(size(mat,1),size(mat,2))
    mat(:,:)=0.0_sp
    do i=1,n
       mat(i,i)=1.0_sp
    end do
  end subroutine unit_matrix
  !BL
  function upper_triangle(j,k,extra)
    integer(I4B), intent(IN) :: j,k
    integer(I4B), optional, intent(IN) :: extra
    logical(LGT), dimension(j,k) :: upper_triangle
    integer(I4B) :: n
    n=0
    if (present(extra)) n=extra
    upper_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) < n)
  end function upper_triangle
  !BL
  function lower_triangle(j,k,extra)
    integer(I4B), intent(IN) :: j,k
    integer(I4B), optional, intent(IN) :: extra
    logical(LGT), dimension(j,k) :: lower_triangle
    integer(I4B) :: n
    n=0
    if (present(extra)) n=extra
    lower_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) > -n)
  end function lower_triangle
  !BL
  function vabs(v)
    real(SP), dimension(:), intent(IN) :: v
    real(SP) :: vabs
    vabs=sqrt(dot_product(v,v))
  end function vabs
  !BL
end module nrutil
