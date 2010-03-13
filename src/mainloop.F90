!> \file mainloop.F90

!> \brief the main program loop
!! 
!! Loops through the snapshots calling the raytracing and output routines 
!<

module mainloop_mod
  use myf90_mod

  ! routines
  use main_input_mod, only: readin_snapshot
  use oct_tree_mod, only: buildtree, setparticleorder
  use raylist_mod, only: prepare_raysearch, kill_raylist, trace_ray
  use ray_mod, only: make_source_ray, make_recomb_ray
  use ray_mod, only: ray_type, raystat_type, raystatbuffsize
  use ion_temperature_update, only: update_raylist, update_no_hits
  use mt19937_mod, only: genrand_real1
  use output_mod, only: output_total_snap, ion_frac_out
  
  ! variables
  use global_mod, only: psys
  use global_mod, only: globalraylist
  use global_mod, only: tree
  use global_mod, only: GV
  use global_mod, only: PLAN
  use global_mod, only: gconst
  
  implicit none
  
  integer(i8b), parameter :: one = 1
  
contains
  
  !> this is the main driver of SPHRAY
  !======================================
  subroutine mainloop()
    implicit none
    
    character(clen), parameter :: myname="mainloop"
    logical, parameter :: crash=.true.
    integer, parameter :: verb=1
    character(clen) :: str,fmt

    type(raystat_type) :: raystats(raystatbuffsize)
    integer(i8b) :: raystatcnt

    !  local counters 
    !-----------------
    
    integer(i8b) :: snapn !< snapshot counter
    integer(i8b) :: rayn  !< ray counter
    integer(i8b) :: srcn  !< source counter
    
    real(r8b) :: rn       !< random number
    real(r8b) :: MB       !< MBs for memory consumption tracking
    
    ! work variables
    !-----------------
    logical :: srcray
    type(ray_type) :: ray
    real(r8b) :: outmark
    
#ifdef incHrec
    integer(i8b) :: rindx
    integer(i8b) :: pindx
    integer(i8b) :: lasthitcheck
#endif
    
    
    raystatcnt = 0
    
    ! loop over the snapshots 
    !=========================
    snaps: do snapn = GV%StartSnapNum, GV%EndSnapNum
              
       call readin_snapshot()
       
#ifdef incHrec
       if (snapn == GV%StartSnapNum) then
          allocate( reclist(psys%npar) )
       end if
       reclist = 0
       GV%recpt = 0
       psys%par(:)%OnRecList = .false.
#endif
       
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
       
       
       ! begin ray tracing 
       !------------------------- 
       src_rays: do rayn = one, PLAN%snap(snapn)%SrcRays

          
          GV%rayn     = GV%rayn + 1
          GV%src_rayn = GV%src_rayn + 1
          
          
          !  select a source randomly (weighted by their luminosity)
          rn = genrand_real1() * psys%src(size(psys%src))%Lcdf
          srcn=1
          do while(psys%src(srcn)%Lcdf.LT.rn)
             srcn=srcn+1
             if(srcn.GT.size(psys%src)) then
                write(*,*) srcn, rn, psys%src(size(psys%src))%Lcdf, size(psys%src)
                stop "src num > number of sources in mainloop.f90"
             endif
          enddo
          
          GV%time_code = GV%time_code + GV%dtray_code
          GV%time_s    = GV%time_s    + GV%dtray_s
          GV%time_elapsed_code = GV%time_elapsed_code + GV%dtray_code
          GV%time_elapsed_s = GV%time_elapsed_s + GV%dtray_s
          GV%TotalSourceRaysCast = GV%TotalSourceRaysCast + 1      
          
          
          !  create a source ray and calc the impacts
          call make_source_ray(psys%src(srcn), GV%rayn, GV%dtray_s, GV%Lunit, psys%box, ray)


!          ray%dir = (/1.,0.,0./)
!          write(*,"(A,3F6.3)") "ray start: ", ray%start


          globalraylist%ray = ray
          call trace_ray(globalraylist%ray, globalraylist, psys, tree) 
          
          GV%TotalPhotonsCast = GV%TotalPhotonsCast + globalraylist%ray%pini
          GV%IonizingPhotonsPerSec = GV%TotalPhotonsCast / GV%time_elapsed_s
          
          srcray = .true.
          call update_raylist(globalraylist,psys%par,psys%box,srcray)
          
          
          if (GV%raystats) then
             
             raystatcnt = raystatcnt + 1
             
             raystats(raystatcnt)%srcn  = srcn
             raystats(raystatcnt)%dir   = globalraylist%ray%dir                                 
             raystats(raystatcnt)%freq  = globalraylist%ray%freq
             raystats(raystatcnt)%temit = GV%time_elapsed_s * gconst%sec_per_megayear
             raystats(raystatcnt)%pini  = globalraylist%ray%pini
             raystats(raystatcnt)%pcnt  = globalraylist%ray%pcnt
             raystats(raystatcnt)%pindx = globalraylist%intersection(globalraylist%lastnnb)%pindx
             raystats(raystatcnt)%b     = globalraylist%intersection(globalraylist%lastnnb)%b
             raystats(raystatcnt)%d     = globalraylist%intersection(globalraylist%lastnnb)%d
             
             if (raystatcnt == raystatbuffsize) then
                write(GV%raystatlun) raystats
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
          if ( mod(GV%rayn, GV%NraysUpdateNoHits)==0 ) then
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
             
             ! set proper time marker unit
             if ( trim(GV%OutputTiming) == "standard" ) then
                
                outmark = GV%time_code
                
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
    
    
  end subroutine mainloop
  
end module mainloop_mod
