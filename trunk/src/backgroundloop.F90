!> \file backgroundloop.F90

!> \brief the background program loop
!! 
!! Loops through the snapshots calling the raytracing and output routines 
!<

module backgroundloop_mod
  use myf90_mod

  ! routines
  use main_input_mod, only: readin_snapshot
  use oct_tree_mod, only: buildtree
  use oct_tree_mod, only: setparticleorder
  use raylist_mod, only: prepare_raysearch
  use output_mod, only: output_total_snap
  use ray_mod, only: make_healpix_ray
  use particle_system_mod, only: particle_info_to_screen

  use ion_table_class

  use raylist_mod, only: kill_raylist, trace_ray
  use ray_mod, only: make_source_ray, make_recomb_ray
  use ray_mod, only: ray_type, raystat_type, raystatbuffsize
  use ion_temperature_update, only: update_raylist, update_no_hits
  use mt19937_mod, only: genrand_real1
  use output_mod, only: ion_frac_out
  use global_mod, only: set_dt_from_dtcode, set_time_elapsed_from_itime

  
  ! variables and types
  use global_mod, only: psys
  use global_mod, only: globalraylist
  use global_mod, only: tree
  use global_mod, only: GV
  use global_mod, only: PLAN
  use global_mod, only: gconst
  use global_mod, only: saved_gheads
  
  use raylist_mod, only: raylist_type
  use ionpar_mod

  implicit none
  
  real(r8b), parameter :: zero = 0.0d0
  real(r8b), parameter :: four_pi = 12.5663706d0
  integer(i4b), parameter :: nside = 2
  integer(i4b), parameter :: npix = 12 * nside**2
  real(r8b), parameter :: sky_frac = 1.0d0 / npix

  

  
contains
  




  !> this is the main driver of SPHRAY
  !======================================
  subroutine backgroundloop()
    implicit none
    
    character(clen), parameter :: myname="backgroundloop"
    logical, parameter :: crash=.true.
    integer, parameter :: verb=1
    character(clen) :: str,fmt

    type(raystat_type) :: raystats(raystatbuffsize)
    integer(i8b) :: raystatcnt

    !  local counters 
    !-----------------
    
    integer(i8b) :: snapn    !< snapshot counter
    integer(i8b) :: rayn     !< ray counter
    integer(i8b) :: srcn     !< source counter
    integer(i8b) :: slv_indx ! index of particle being solved
    integer(i8b) :: hit_indx ! index of a ray particle
    integer(i8b) :: ihit     !< intersection counter
    
    real(r8b) :: rn       !< random number
    real(r8b) :: MB       !< MBs for memory consumption tracking
    
    ! work variables
    !-----------------
    logical :: srcray
    type(ray_type) :: ray
    real(r8b) :: outmark
    
    integer(i8b) :: pindx
    integer(i8b) :: lasthitcheck


    type(raylist_type) :: raylist_group(0:npix-1)   
    integer(i4b) :: ipix
    type(ionpart_type) :: hit_ionpar
    type(ionpart_type) :: slv_ionpar
    real(r8b) :: tausum(0:npix-1)

    character(clen) :: ion_table_file
    type( ion_table_type ) :: itab
    type( mini_spec_type ) :: minispec
    real(r8b) :: redshift
    real(r8b) :: GHIint_thin
    real(r8b) :: GHIint_shield
    real(r8b) :: GHItable
    real(r8b) :: GHItotal_thin
    real(r8b) :: GHItotal_shield

    real(r8b) :: min_logryd
    real(r8b) :: max_logryd

    integer(i4b) :: iter
    
    ion_table_file = '../data/ionization_tables/h1.hdf5'

    
    ! loop over the snapshots 
    !=========================
    snaps: do snapn = GV%StartSnapNum, GV%EndSnapNum

       !  read in particle and source snapshot
       !---------------------------------------------------------------- 
       call readin_snapshot()
              

       !  get background flux from Cloudy table
       !---------------------------------------------------------------- 
       call read_ion_table_file( ion_table_file, itab )
       redshift = saved_gheads(snapn,0)%z

       min_logryd = 0.0d0
       max_logryd = 1.0d0
       call set_mini_spectrum( min_logryd, max_logryd, redshift, itab, minispec)
       GHIint_thin = gammaHI_from_mini_spec_thin( minispec )       
       GHItable = minispec%gammaHI

       write(*,*) 'GHI int', GHIint_thin
       write(*,*) 'GHI tab', GHItable
       write(*,*) 'GHIi/GHIt = ', GHIint_thin/GHItable


       !  build oct tree.  only need to do this once per snap (for now)
       !----------------------------------------------------------------
       call buildtree(psys,tree,MB,GV%PartPerCell)
       GV%MB = GV%MB + MB
       call setparticleorder(psys, tree)             
       call prepare_raysearch(psys, globalraylist)
       

          
       if(GV%JustInit) then
          write(str,"(A,F10.2)") "total memory allocated [MB] = ", GV%MB
          call mywrite(str,verb)
          call mywrite("just initializing", verb)
          call mywrite("",verb)
          stop
       end if

       
       if (GV%DoInitialOutput) then
          GV%OutputIndx = 0
          call output_total_snap(psys%par)      
          GV%OutputIndx = 1
       end if
       call particle_info_to_screen(psys)


       over_iterations: do iter = 1, 5
       
       
          ! begin ray tracing 
          !------------------------- 
          over_pars: do slv_indx = 1, size(psys%par)
             
             GV%itime = GV%itime + 1
             call set_time_elapsed_from_itime( GV )
             GV%IonizingPhotonsPerSec = GV%TotalPhotonsCast / GV%time_elapsed_s
             
             call initialize_ionpar( slv_ionpar, psys%par(slv_indx), slv_indx, &
                  srcray= .false., He= .false. )
                          
             tausum = zero
             GHItotal_thin = zero
             GHItotal_shield = zero
             over_pixels: do ipix = 0, npix-1
                GV%src_rayn = GV%src_rayn + 1
                call make_healpix_ray( psys%par(slv_indx), nside, ipix, psys%box, ray )
                raylist_group(ipix)%ray = ray
                call trace_ray(ray, raylist_group(ipix), psys, tree) 
                
                over_hits: do ihit = raylist_group(ipix)%nnb, 2, -1
                   hit_indx = raylist_group(ipix)%intersection(ihit)%pindx
                   call initialize_ionpar( hit_ionpar, psys%par(hit_indx), hit_indx, &
                        srcray= .false., He= .false., &
                        raylist= raylist_group(ipix), impact= ihit)
                   call set_taus(hit_ionpar, He=.false.)
                   tausum(ipix) = tausum(ipix) + hit_ionpar%tauHI
                   hit_ionpar%nHI = hit_ionpar%nH * hit_ionpar%xHI
                                                        
                end do over_hits
                
                GHIint_shield = gammaHI_from_mini_spec_shield( minispec, tausum(ipix) )
                
                GHItotal_thin   = GHItotal_thin   + sky_frac * GHIint_thin
                GHItotal_shield = GHItotal_shield + sky_frac * GHIint_shield 
        
             end do over_pixels
                         
             slv_ionpar%gammaHI = GHItotal_shield
             
             call set_ionpar_xH_eq( slv_ionpar )
             call ionpar2par( slv_ionpar, psys%par(slv_indx) )
             psys%par(slv_indx)%gammaHI = GHItotal_shield
             psys%par(slv_indx)%time = 1.0d0
             
          end do over_pars
          

          call particle_info_to_screen(psys)
          call output_total_snap(psys%par)      
          GV%OutputIndx = GV%OutputIndx + 1



       end do over_iterations










       stop



       src_rays: do rayn = 1, PLAN%snap(snapn)%SrcRays

          
          GV%rayn                = GV%rayn + 1
          GV%src_rayn            = GV%src_rayn + 1
          GV%TotalSourceRaysCast = GV%TotalSourceRaysCast + 1                
          
          !  select a source randomly (weighted by their luminosity)
          rn = genrand_real1() * psys%src(size(psys%src))%Lcdf
          srcn=1
          do while(psys%src(srcn)%Lcdf.LT.rn)
             srcn=srcn+1
             if(srcn.GT.size(psys%src)) then
                write(*,*) srcn, rn, psys%src(size(psys%src))%Lcdf, size(psys%src)
                stop "src num > number of sources in backgroundloop.f90"
             endif
          enddo
                    
          !  create a source ray and calc the impacts
          call make_source_ray(psys%src(srcn), GV%rayn, GV%WallSampling, GV%dt_s, GV%Lunit, psys%box, ray)
          
        
          GV%itime = GV%itime + ray%weight
          call set_time_elapsed_from_itime( GV )
          GV%IonizingPhotonsPerSec = GV%TotalPhotonsCast / GV%time_elapsed_s



          globalraylist%ray = ray
          call trace_ray(globalraylist%ray, globalraylist, psys, tree) 
          
          GV%TotalPhotonsCast = GV%TotalPhotonsCast + globalraylist%ray%pini

          
          srcray = .true.
          call update_raylist(globalraylist,psys%par,psys%box,srcray)
          
          
          if (GV%raystats) then
             
             raystatcnt = raystatcnt + 1
             
             raystats(raystatcnt)%srcn  = srcn
             raystats(raystatcnt)%start = globalraylist%ray%start  
             raystats(raystatcnt)%ryd   = globalraylist%ray%freq
             
             if (raystatcnt == raystatbuffsize) then
                write(GV%raystatlun) raystats
                flush(GV%raystatlun)
                raystatcnt = 0
             end if
                          
          end if
          
          
          ! recombination ray loop
#ifdef incHrec
          if ( mod(GV%rayn,GV%NrecombRaysPerSourceRay)==0 ) then
             
             do rindx = 1,GV%recpt
                pindx = reclist(rindx)
                
                ! reset this
                psys%par(pindx)%OnRecList = .false.
                
                GV%TotalDiffuseRaysCast = GV%TotalDiffuseRaysCast + 1
                
                !               xHIIrc = psys%par(pindx)%xHIIrc
                
                
                call make_recomb_ray( psys%par(pindx), GV%IsoMass, GV%cgs_mass, &
                     GV%H_mf, GV%He_mf, ray)
                
                call trace_ray(ray,globalraylist,psys,tree) 
                
                
                psys%par(pindx)%xHIIrc = 0.0
#ifdef incHe
                psys%par(pindx)%xHeIIrc = 0.0
                psys%par(pindx)%xHeIIIrc = 0.0
#endif
                
                lasthitcheck = psys%par(pindx)%lasthit 
                
                srcray = .false.
                call update_raylist(globalraylist,psys,srcray)
                
                if (.not. psys%par(pindx)%lasthit == lasthitcheck) then
                   write(*,*) "lasthit check failed"
                   stop
                end if
                
             end do
             
             ! reset reclist
             reclist(1:GV%recpt) = 0
             GV%recpt = 0
             
          end if
          
#endif
          

          if (GV%NraysUpdateNoHits == 0) GV%NraysUpdateNoHits = huge(1_i8b)
          if ( mod(GV%itime, GV%NraysUpdateNoHits)==0 ) then
             write(*,*) "updating particles not hit by rays"
             call update_no_hits(psys, tree)
          end if
          
          

          
          
          !        output routines
          !------------------------
          
          !        check if this time step requires a small ionization frac output
          !        these are small outputs done often to monitor certain quantites
          if( mod(GV%rayn,GV%IonFracOutRays)==0 ) then
             call ion_frac_out(psys, tree )
          end if
          
          
          ! check if this time step requires a full output
          if ( GV%OutputIndx <= GV%NumTotOuts ) then
             
             ! set correct time marker unit
             if ( trim(GV%OutputTiming) == "standard" ) then
                
                outmark = GV%start_time_code + GV%time_elapsed_code
                
             else if ( trim(GV%OutputTiming) == "forced" ) then
                
                if (trim(GV%ForcedUnits) == "mwionfrac") then
                   outmark = GV%mwionfrac
                else 
                   outmark = GV%time_elapsed_code
                end if
                
             else
                
                write(*,*) "output type ", trim(GV%OutputTiming), "not recognized"
                stop 
                
             end if


             
             ! check outmark against tabulated output "times"
             if ( outmark >= PLAN%OutputTimes(GV%OutputIndx) ) then
                call output_total_snap(psys%par)
                GV%OutputIndx = GV%OutputIndx + 1 
             end if
             
          end if
          
          
          ! if we are on the last ray and we havent gotten to the last 
          ! output, do a full output.
          if ( snapn == GV%EndSnapNum ) then
             if ( GV%OutputIndx <= GV%NumTotOuts ) then
                if ( GV%src_rayn==PLAN%snap(snapn)%SrcRays ) then
                   write(*,*) "doing an output on the last ray"
                   call output_total_snap(psys%par)
                end if
             end if
          end if
          
          
          
          
       end do src_rays


       
       ! free up the memory from the globalraylist.
       call kill_raylist(globalraylist)
       
    end do snaps
    
    close(GV%ionlun)
    
    
    if (GV%raystats) then
       close(GV%raystatlun)
    end if
    
    
  end subroutine backgroundloop
  
end module backgroundloop_mod
