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
  use particle_system_mod, only: create_particle_density_access_list
  use particle_system_mod, only: create_particle_random_access_list

  use ion_table_class

  use raylist_mod, only: trace_ray
  use ray_mod, only: ray_type
  use output_mod, only: ion_frac_out
  use global_mod, only: set_dt_from_dtcode
  use global_mod, only: set_time_elapsed_from_itime

  
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
  integer(i4b), parameter :: nside = 1
  integer(i4b), parameter :: npix = 12 * nside**2
  real(r8b), parameter :: sky_frac = 1.0d0 / npix
  real(r8b), parameter :: TAU_MAX = 30.0d0


  
contains
  


!> this is the main driver of SPHRAY
!======================================
subroutine backgroundloop()    
  character(clen), parameter :: myname="backgroundloop"
  logical, parameter :: crash=.true.
  integer, parameter :: verb=1
  character(clen) :: str,fmt

  !  local counters 
  !-----------------  
  integer(i8b) :: snapn     !< snapshot counter
  integer(i8b) :: rayn      !< ray counter
  integer(i8b) :: list_indx !< index in particle access list
  integer(i8b) :: par_indx  !< index in particle system of particle being solved for
  integer(i8b) :: hit_indx  !< index of a ray particle
  integer(i8b) :: ihit      !< intersection counter
        
  real(r8b) :: MB

  ! work variables
  !-----------------
  type(ray_type) :: ray
  type(raylist_type) :: raylist_group(0:npix-1)   
  integer(i4b) :: ipix
  type(tau_particle_type) :: tau_par
  type(bckgnd_particle_type) :: bck_par
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
  integer(i4b) :: err

  integer(i4b) :: i
  integer(i4b) :: indx
  character(clen) :: rayfile
  
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
       

          ! create list to access particles from east to most dense
          !-----------------------------------------------------------
          if (iter==1) then
             call create_particle_density_access_list( psys )
          else
             call create_particle_random_access_list( psys )
          endif

       
          ! begin ray tracing 
          !------------------------- 
          over_pars: do list_indx = 1, size(psys%par)

             if (mod(list_indx,size(psys%par)/100)==0) then
                write(*,*) 'done ', list_indx, 'of ', size(psys%par)
             endif

             par_indx = psys%acc_list(list_indx)
             
             GV%itime = GV%itime + 1
             call set_time_elapsed_from_itime( GV )
             GV%IonizingPhotonsPerSec = GV%TotalPhotonsCast / GV%time_elapsed_s
                                       
             tausum          = zero
             GHItotal_thin   = zero
             GHItotal_shield = zero
             over_pixels: do ipix = 0, npix-1
                
                GV%src_rayn = GV%src_rayn + 1
                if (GV%MaxRayDist > zero) then
                   call make_healpix_ray( psys%par(par_indx), nside, ipix, psys%box, ray, length=GV%MaxRayDist )
                else
                   call make_healpix_ray( psys%par(par_indx), nside, ipix, psys%box, ray)
                endif
                raylist_group(ipix)%ray = ray

                call trace_ray(ray, raylist_group(ipix), psys, tree, dosort=.false.) 
                                
                over_hits: do ihit = 1, raylist_group(ipix)%nnb

                   hit_indx = raylist_group(ipix)%intersection(ihit)%pindx
                   if (hit_indx == par_indx) cycle
                   tau_par = initialize_tau_particle( psys%par(hit_indx), raylist_group(ipix)%intersection(ihit) )
                   tausum(ipix) = tausum(ipix) + tau_par%tauHI_th
                   if (tausum(ipix) > TAU_MAX) exit over_hits
                                                        
                end do over_hits

                if (tausum(ipix) > TAU_MAX) then
                   GHIint_shield = zero
                else
                   GHIint_shield = gammaHI_from_mini_spec_shield( minispec, tausum(ipix) )
                endif

                GHItotal_thin   = GHItotal_thin   + sky_frac * GHIint_thin
                GHItotal_shield = GHItotal_shield + sky_frac * GHIint_shield 
        
             end do over_pixels
                             
             bck_par = initialize_bckgnd_particle( psys%par(par_indx), GHItotal_shield )        
             err = set_bckgnd_particle_xH_eq(bck_par)

             psys%par(par_indx)%xHI     = bck_par%xHI
             psys%par(par_indx)%xHII    = bck_par%xHII
             psys%par(par_indx)%gammaHI = GHItotal_shield
             psys%par(par_indx)%time    = 1.0d0
             
          end do over_pars
          

          call particle_info_to_screen(psys)
          call output_total_snap(psys%par)      
          GV%OutputIndx = GV%OutputIndx + 1


       end do over_iterations
    
    end do snaps
    
       
  end subroutine backgroundloop
  
end module backgroundloop_mod
