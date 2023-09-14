!> \file GFS_stochastics.f90
!! This file contains code previously in GFS_stochastics_driver.

!>\defgroup gfs_stoch GFS Stochastics Physics Module
!! This module
    module GFS_stochastics

      contains

!> \section arg_table_GFS_stochastics_init Argument Table
!! \htmlinclude GFS_stochastics_init.html
!!
!>\section gfs_stochy_general GFS_stochastics_init General Algorithm
!! This is the GFS stochastic physics initialization.
!! -# define vertical tapering for CA global
      subroutine GFS_stochastics_init (si,vfact_ca,km,do_ca,ca_global, errmsg, errflg)

      use machine,               only: kind_phys

      implicit none
      real(kind_phys), dimension(:),         intent(in)    :: si
      real(kind_phys), dimension(:),         intent(inout) :: vfact_ca
      integer,                               intent(in)    :: km
      logical,                               intent(in)    :: do_ca
      logical,                               intent(in)    :: ca_global
      character(len=*),                      intent(out)   :: errmsg
      integer,                               intent(out)   :: errflg
      integer :: k,nz

      errmsg = ''
      errflg = 0
      if (do_ca .and. ca_global) then
         nz=min(km,size(vfact_ca))
         vfact_ca(:)=0.0
         do k=1,nz
            if (si(k) .lt. 0.1 .and. si(k) .gt. 0.025) then
               vfact_ca(k) = (si(k)-0.025)/(0.1-0.025)
            else if (si(k) .lt. 0.025) then
               vfact_ca(k) = 0.0
            else
               vfact_ca(k) = 1.0
            endif
         enddo
         vfact_ca(2)=vfact_ca(3)*0.5
         vfact_ca(1)=0.0
      endif
      end subroutine GFS_stochastics_init

      subroutine GFS_stochastics_finalize()
      end subroutine GFS_stochastics_finalize


!> \section arg_table_GFS_stochastics_run Argument Table
!! \htmlinclude GFS_stochastics_run.html
!!
!>\section gfs_stochy_general GFS_stochastics_run General Algorithm
!! This is the GFS stochastic physics driver.
!! Routines are called prior to radiation and physics steps to handle:
!! -# sets up various time/date variables
!! -# sets up various triggers
!! -# defines random seed indices for radiation (in a reproducible way)
!! -# interpolates coefficients for prognostic ozone calculation
!! -# performs surface data cycling via the GFS gcycle routine
      subroutine GFS_stochastics_run (im, km, kdt, delt, do_sppt, pert_mp, use_zmtnblck, &
                                      do_shum ,do_skeb, do_mlp,do_mlp_cnv,do_mlp_mp,do_mlp_pbl,do_mlp_shalcnv,&
                                      do_ca,ca_global,ca1,vfact_ca,    &
                                      zmtnblck, sppt_wts, skebu_wts, skebv_wts, shum_wts,&
                                      diss_est, ugrs, vgrs, tgrs, qgrs_wv,               &
                                      qgrs_cw, qgrs_rw, qgrs_sw, qgrs_iw, qgrs_gl,       &
                                      gu0, gv0, gt0, gq0_wv, dtdtnp,                     &
!mlp
                                      mlp_pert_tcnv,mlp_pert_qcnv,mlp_pert_ucnv,mlp_pert_vcnv,&
                                      mlp_pert_tmp,mlp_pert_qmp, &
                                      mlp_pert_tpbl,mlp_pert_qpbl,mlp_pert_upbl,mlp_pert_vpbl,&
                                      mlp_pert_tshalcnv,mlp_pert_qshalcnv,mlp_pert_ushalcnv,mlp_pert_vshalcnv,&
                                      gq0_cw, gq0_rw, gq0_sw, gq0_iw, gq0_gl,            &
                                      rain, rainc, tprcp, totprcp, cnvprcp,              &
                                      totprcpb, cnvprcpb, cplflx,                        &
                                      rain_cpl, snow_cpl, drain_cpl, dsnow_cpl,          &
                                      ntcw,ntrw,ntsw,ntiw,ntgl,                          &
                                      errmsg, errflg)

         use machine,               only: kind_phys

         implicit none

         integer,                               intent(in)    :: im
         integer,                               intent(in)    :: km
         integer,                               intent(in)    :: kdt
         real(kind_phys),                       intent(in)    :: delt     
         logical,                               intent(in)    :: do_sppt
         logical,                               intent(in)    :: pert_mp
         logical,                               intent(in)    :: do_ca
         logical,                               intent(in)    :: ca_global
         logical,                               intent(in)    :: use_zmtnblck
         logical,                               intent(in)    :: do_shum
         logical,                               intent(in)    :: do_skeb
!jwb/sam mlp
         logical,                               intent(in)    :: do_mlp
         logical,                               intent(in)    :: do_mlp_cnv
         logical,                               intent(in)    :: do_mlp_mp
         logical,                               intent(in)    :: do_mlp_shalcnv
         logical,                               intent(in)    :: do_mlp_pbl
!jwb/sam mlp
         real(kind_phys), dimension(:),         intent(in)    :: zmtnblck
         ! sppt_wts only allocated if do_sppt == .true.
         real(kind_phys), dimension(:,:),       intent(inout) :: sppt_wts
         ! skebu_wts, skebv_wts only allocated if do_skeb == .true.
         real(kind_phys), dimension(:,:),       intent(in)    :: skebu_wts
         real(kind_phys), dimension(:,:),       intent(in)    :: skebv_wts
         ! shum_wts only allocated if do_shum == .true.
         real(kind_phys), dimension(:,:),       intent(in)    :: shum_wts
         real(kind_phys), dimension(:,:),       intent(in)    :: diss_est
         real(kind_phys), dimension(:,:),       intent(in)    :: ugrs
         real(kind_phys), dimension(:,:),       intent(in)    :: vgrs
         real(kind_phys), dimension(:,:),       intent(in)    :: tgrs
         real(kind_phys), dimension(:,:),       intent(in)    :: qgrs_wv
         real(kind_phys), dimension(:,:),       intent(in)    :: qgrs_cw
         real(kind_phys), dimension(:,:),       intent(in)    :: qgrs_rw
         real(kind_phys), dimension(:,:),       intent(in)    :: qgrs_sw
         real(kind_phys), dimension(:,:),       intent(in)    :: qgrs_iw
         real(kind_phys), dimension(:,:),       intent(in)    :: qgrs_gl
!jwb/mlp
         real(kind_phys), dimension(1:im,1:km), intent(inout)    :: mlp_pert_ucnv
         real(kind_phys), dimension(1:im,1:km), intent(inout)    :: mlp_pert_vcnv
         real(kind_phys), dimension(1:im,1:km), intent(inout)    :: mlp_pert_tcnv
         real(kind_phys), dimension(1:im,1:km), intent(inout)    :: mlp_pert_qcnv
         real(kind_phys), dimension(1:im,1:km), intent(inout)    :: mlp_pert_tmp
         real(kind_phys), dimension(1:im,1:km), intent(inout)    :: mlp_pert_qmp
         real(kind_phys), dimension(1:im,1:km), intent(inout)    :: mlp_pert_ushalcnv
         real(kind_phys), dimension(1:im,1:km), intent(inout)    :: mlp_pert_vshalcnv
         real(kind_phys), dimension(1:im,1:km), intent(inout)    :: mlp_pert_tshalcnv
         real(kind_phys), dimension(1:im,1:km), intent(inout)    :: mlp_pert_qshalcnv
         real(kind_phys), dimension(1:im,1:km), intent(inout)    :: mlp_pert_upbl
         real(kind_phys), dimension(1:im,1:km), intent(inout)    :: mlp_pert_vpbl
         real(kind_phys), dimension(1:im,1:km), intent(inout)    :: mlp_pert_tpbl
         real(kind_phys), dimension(1:im,1:km), intent(inout)    :: mlp_pert_qpbl
!mlp
         real(kind_phys), dimension(:,:),       intent(inout) :: gu0
         real(kind_phys), dimension(:,:),       intent(inout) :: gv0
         real(kind_phys), dimension(:,:),       intent(inout) :: gt0
         real(kind_phys), dimension(:,:),       intent(inout) :: gq0_wv
         real(kind_phys), dimension(:,:),       intent(inout) :: gq0_cw
         real(kind_phys), dimension(:,:),       intent(inout) :: gq0_rw
         real(kind_phys), dimension(:,:),       intent(inout) :: gq0_sw
         real(kind_phys), dimension(:,:),       intent(inout) :: gq0_iw
         real(kind_phys), dimension(:,:),       intent(inout) :: gq0_gl
         integer, intent(in) ::      ntcw
         integer, intent(in) ::      ntrw
         integer, intent(in) ::      ntsw
         integer, intent(in) ::      ntiw
         integer, intent(in) ::      ntgl
         real(kind_phys), dimension(:,:),       intent(inout) :: dtdtnp
         real(kind_phys), dimension(:),         intent(in)    :: rain
         real(kind_phys), dimension(:),         intent(in)    :: rainc
         real(kind_phys), dimension(:),         intent(inout) :: tprcp
         real(kind_phys), dimension(:),         intent(inout) :: totprcp
         real(kind_phys), dimension(:),         intent(inout) :: cnvprcp
         real(kind_phys), dimension(:),         intent(inout) :: totprcpb
         real(kind_phys), dimension(:),         intent(inout) :: cnvprcpb
         logical,                               intent(in)    :: cplflx
         ! rain_cpl, snow_cpl only allocated if cplflx == .true. or cplchm == .true.
         real(kind_phys), dimension(:),         intent(inout) :: rain_cpl
         real(kind_phys), dimension(:),         intent(inout) :: snow_cpl
         ! drain_cpl, dsnow_cpl only allocated if cplflx == .true. or cplchm == .true.
         real(kind_phys), dimension(:),         intent(in)    :: drain_cpl
         real(kind_phys), dimension(:),         intent(in)    :: dsnow_cpl
         real(kind_phys), dimension(:),         intent(in)    :: vfact_ca
         real(kind_phys), dimension(:),         intent(in)    :: ca1
         character(len=*),                      intent(out)   :: errmsg
         integer,                               intent(out)   :: errflg

         !--- local variables
         integer :: k, i
         real(kind=kind_phys) :: upert, vpert, tpert, qpert, qnew, sppt_vwt
!jwb/sam mlp 
         real(kind=kind_phys) :: upert_cnv, vpert_cnv, tpert_cnv, qpert_cnv
         real(kind=kind_phys) :: upert_mp, vpert_mp, tpert_mp, qpert_mp
         real(kind=kind_phys) :: upert_shalcnv, vpert_shalcnv, tpert_shalcnv, qpert_shalcnv
         real(kind=kind_phys) :: upert_pbl, vpert_pbl, tpert_pbl, qpert_pbl
!mlp
         real(kind=kind_phys), dimension(1:im,1:km) :: ca

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (do_sppt) then
           do k=1,km
             do i=1,im
               sppt_vwt=1.0
               if (zmtnblck(i).EQ.0.0) then
                  sppt_vwt=1.0
               else 
                  if (k.GT.zmtnblck(i)+2) then
                     sppt_vwt=1.0
                  endif
                  if (k.LE.zmtnblck(i)) then
                     sppt_vwt=0.0
                  endif
                  if (k.EQ.zmtnblck(i)+1) then
                     sppt_vwt=0.333333
                  endif
                  if (k.EQ.zmtnblck(i)+2) then
                     sppt_vwt=0.666667
                  endif
               endif

               if (use_zmtnblck)then
                  sppt_wts(i,k)=(sppt_wts(i,k)-1)*sppt_vwt+1.0
               endif

         if (do_mlp) then
          !print*,'doing pert gfs_stocastic',do_mlp_cnv,do_mlp_mp
! for total
            if(do_mlp_cnv)then
              upert_cnv = mlp_pert_ucnv(i,k)
              vpert_cnv = mlp_pert_vcnv(i,k)
              tpert_cnv = mlp_pert_tcnv(i,k)
              qpert_cnv = mlp_pert_qcnv(i,k)
      
            else
              upert_cnv = 0.
              vpert_cnv = 0.
              tpert_cnv = 0.
              qpert_cnv = 0.
            endif
            if(do_mlp_mp)then
              upert_mp = 0.
              vpert_mp = 0.
              tpert_mp = mlp_pert_tmp(i,k)
              qpert_mp = mlp_pert_qmp(i,k)
            else
              upert_mp = 0.
              vpert_mp = 0.
              tpert_mp = 0.
              qpert_mp = 0.
            endif
!
            if(do_mlp_shalcnv)then
              upert_shalcnv = mlp_pert_ushalcnv(i,k)
              vpert_shalcnv = mlp_pert_vshalcnv(i,k)
              tpert_shalcnv = mlp_pert_tshalcnv(i,k)
              qpert_shalcnv = mlp_pert_qshalcnv(i,k)
            else
              upert_shalcnv = 0.
              vpert_shalcnv = 0.
              tpert_shalcnv = 0.
              qpert_shalcnv = 0.
            endif

! make mlppblfac?
           if(do_mlp_pbl)then
              upert_pbl = 0.6*mlp_pert_upbl(i,k)
              vpert_pbl = 0.6*mlp_pert_vpbl(i,k)
              tpert_pbl = 0.6*mlp_pert_tpbl(i,k)
              qpert_pbl = 0.6*mlp_pert_qpbl(i,k)
           ! if(mlp_pert_qpbl(i,k).ne.0.)then
          !print*,'pbl pert t mlp_pert_tpbl in gfs_stocastic',mlp_pert_tpbl(i,k)
          !  endif
            else
              upert_pbl = 0.
              vpert_pbl = 0.
              tpert_pbl = 0.
              qpert_pbl = 0.
            endif

           upert=upert_cnv + upert_mp + upert_shalcnv + upert_pbl
           vpert=vpert_cnv + vpert_mp + vpert_shalcnv + vpert_pbl
           tpert=tpert_cnv + tpert_mp + tpert_shalcnv + tpert_pbl
           qpert=qpert_cnv + qpert_mp + qpert_shalcnv + qpert_pbl


           gu0(i,k)  = gu0(i,k)+upert
           gv0(i,k)  = gv0(i,k)+vpert
         !negative humidity check

           qnew = gq0_wv(i,k)+qpert
           if (qnew >= 1.0e-10) then
             gq0_wv(i,k) = qnew
             gt0(i,k)   = gt0(i,k) + tpert
           endif
         ! if(qpert.gt.0)then
         ! print*,' pert q not zero mlp_pert_q in gfs_stocastic',qpert,qpert_cnv,qpert_mp,qpert_shalcnv,qpert_pbl,gq0_wv(i,k)
         ! endif
! sam mlp
         else
               !print*,"not doing mlp perts "

               upert = (gu0(i,k) - ugrs(i,k))   * sppt_wts(i,k)
               vpert = (gv0(i,k) - vgrs(i,k))   * sppt_wts(i,k)
               tpert = (gt0(i,k) - tgrs(i,k) - (delt*dtdtnp(i,k))) * sppt_wts(i,k)
               qpert = (gq0_wv(i,k) - qgrs_wv(i,k)) * sppt_wts(i,k)

               gu0(i,k)  = ugrs(i,k)+upert
               gv0(i,k)  = vgrs(i,k)+vpert

               !negative humidity check
               qnew = qgrs_wv(i,k)+qpert
               if (qnew >= 1.0e-10) then
                  gq0_wv(i,k) = qnew
                  gt0(i,k) = tgrs(i,k) + tpert + (delt*dtdtnp(i,k))
               endif
               if (pert_mp) then
                  if (ntcw>0) then
                     qpert = (gq0_cw(i,k) - qgrs_cw(i,k)) * sppt_wts(i,k)
                     qnew = qgrs_cw(i,k)+qpert
                     gq0_cw(i,k) = qnew
                     if (qnew < 0.0) then
                        gq0_cw(i,k) = 0.0
                     endif
                  endif
                  if (ntrw>0) then
                     qpert = (gq0_rw(i,k) - qgrs_rw(i,k)) * sppt_wts(i,k)
                     qnew = qgrs_rw(i,k)+qpert
                     gq0_rw(i,k) = qnew
                     if (qnew < 0.0) then
                        gq0_rw(i,k) = 0.0
                     endif
                  endif
                  if (ntsw>0) then
                     qpert = (gq0_sw(i,k) - qgrs_sw(i,k)) * sppt_wts(i,k)
                     qnew = qgrs_sw(i,k)+qpert
                     gq0_sw(i,k) = qnew
                     if (qnew < 0.0) then
                        gq0_sw(i,k) = 0.0
                     endif
                  endif
                  if (ntiw>0) then
                     qpert = (gq0_iw(i,k) - qgrs_iw(i,k)) * sppt_wts(i,k)
                     qnew = qgrs_iw(i,k)+qpert
                     gq0_iw(i,k) = qnew
                     if (qnew < 0.0) then
                        gq0_iw(i,k) = 0.0
                     endif
                  endif
                  if (ntgl>0) then
                     qpert = (gq0_gl(i,k) - qgrs_gl(i,k)) * sppt_wts(i,k)
                     qnew = qgrs_gl(i,k)+qpert
                     gq0_gl(i,k) = qnew
                     if (qnew < 0.0) then
                        gq0_gl(i,k) = 0.0
                     endif
                  endif
               endif
! mlp
             endif
             enddo
           enddo

           ! instantaneous precip rate going into land model at the next time step
           tprcp(:) = sppt_wts(:,15)*tprcp(:)
           totprcp(:) = totprcp(:) + (sppt_wts(:,15) - 1 )*rain(:)
           ! acccumulated total and convective preciptiation
           cnvprcp(:) = cnvprcp(:) + (sppt_wts(:,15) - 1 )*rainc(:)
           ! bucket precipitation adjustment due to sppt
           totprcpb(:) = totprcpb(:) + (sppt_wts(:,15) - 1 )*rain(:)
           cnvprcpb(:) = cnvprcpb(:) + (sppt_wts(:,15) - 1 )*rainc(:)

           if (cplflx) then
               rain_cpl(:) = rain_cpl(:) + (sppt_wts(:,15) - 1.0)*drain_cpl(:)
               snow_cpl(:) = snow_cpl(:) + (sppt_wts(:,15) - 1.0)*dsnow_cpl(:)
           endif
           !zero out radiative heating tendency for next physics step
           dtdtnp(:,:)=0.0

         endif

         if (do_ca .and. ca_global) then

          !if(kdt == 1)then
          !endif
   
            do k = 1,km
               do i = 1,im
                  sppt_vwt=1.0
                  if (zmtnblck(i).EQ.0.0) then
                     sppt_vwt=1.0
                  else
                     if (k.GT.zmtnblck(i)+2) then
                        sppt_vwt=1.0
                     endif
                     if (k.LE.zmtnblck(i)) then
                        sppt_vwt=0.0
                     endif
                     if (k.EQ.zmtnblck(i)+1) then
                        sppt_vwt=0.333333
                     endif
                     if (k.EQ.zmtnblck(i)+2) then
                        sppt_vwt=0.666667
                     endif
                  endif

                  ca(i,k)=((ca1(i)-1.)*sppt_vwt*vfact_ca(k))+1.0

                  upert = (gu0(i,k)   - ugrs(i,k))   * ca(i,k)
                  vpert = (gv0(i,k)   - vgrs(i,k))   * ca(i,k)
                  tpert = (gt0(i,k)   - tgrs(i,k) - (delt*dtdtnp(i,k))) * ca(i,k)
                  qpert = (gq0_wv(i,k)   - qgrs_wv(i,k)) * ca(i,k)
                  gu0(i,k)  = ugrs(i,k)+upert
                  gv0(i,k)  = vgrs(i,k)+vpert
                  !negative humidity check                                                                                                                                                                                                                     
                  qnew = qgrs_wv(i,k)+qpert
                  if (qnew >= 1.0e-10) then
                     gq0_wv(i,k) = qnew
                     gt0(i,k)   = tgrs(i,k) + tpert + (delt*dtdtnp(i,k))
                  endif
                  if (pert_mp) then
                     if (ntcw>0) then
                        qpert = (gq0_cw(i,k) - qgrs_cw(i,k)) * ca(i,k)
                        qnew = qgrs_cw(i,k)+qpert
                        gq0_cw(i,k) = qnew
                        if (qnew < 0.0) then
                           gq0_cw(i,k) = 0.0
                        endif
                     endif
                     if (ntrw>0) then
                        qpert = (gq0_rw(i,k) - qgrs_rw(i,k)) * ca(i,k)
                        qnew = qgrs_rw(i,k)+qpert
                        gq0_rw(i,k) = qnew
                        if (qnew < 0.0) then
                           gq0_rw(i,k) = 0.0
                        endif
                     endif
                     if (ntsw>0) then
                        qpert = (gq0_sw(i,k) - qgrs_sw(i,k)) * ca(i,k)
                        qnew = qgrs_sw(i,k)+qpert
                        gq0_sw(i,k) = qnew
                        if (qnew < 0.0) then
                           gq0_sw(i,k) = 0.0
                        endif
                     endif
                     if (ntiw>0) then
                        qpert = (gq0_iw(i,k) - qgrs_iw(i,k)) * ca(i,k)
                        qnew = qgrs_iw(i,k)+qpert
                        gq0_iw(i,k) = qnew
                        if (qnew < 0.0) then
                           gq0_iw(i,k) = 0.0
                        endif
                     endif
                     if (ntgl>0) then
                        qpert = (gq0_gl(i,k) - qgrs_gl(i,k)) * ca(i,k)
                        qnew = qgrs_gl(i,k)+qpert
                        gq0_gl(i,k) = qnew
                        if (qnew < 0.0) then
                           gq0_gl(i,k) = 0.0
                        endif
                     endif
                  endif
               enddo
            enddo
       
            ! instantaneous precip rate going into land model at the next time step                                                                                                                                                                         
            tprcp(:) = ca(:,15)*tprcp(:)
            totprcp(:) = totprcp(:) + (ca(:,15) - 1 )*rain(:)
            ! acccumulated total and convective preciptiation                                                                                                                                                                                               
            cnvprcp(:) = cnvprcp(:)      + (ca(:,15) - 1 )*rainc(:)
            ! bucket precipitation adjustment due to sppt                                                                                                                                                                                                   
            totprcpb(:)      = totprcpb(:)      + (ca(:,15) - 1 )*rain(:)
            cnvprcpb(:)      = cnvprcpb(:)      + (ca(:,15) - 1 )*rainc(:)
            
            if (cplflx) then
               rain_cpl(:) = rain_cpl(:) + (ca(:,15) - 1.0)*drain_cpl(:)
               snow_cpl(:) = snow_cpl(:) + (ca(:,15) - 1.0)*dsnow_cpl(:)
            endif
            !zero out radiative heating tendency for next physics step
            dtdtnp(:,:)=0.0


         endif

         if (do_shum) then
           do k=1,km
             gq0_wv(:,k) = gq0_wv(:,k)*(1.0 + shum_wts(:,k))
           end do
         endif
         
         if (do_skeb) then
           do k=1,km
             gu0(:,k) = gu0(:,k)+skebu_wts(:,k)*(diss_est(:,k))
             gv0(:,k) = gv0(:,k)+skebv_wts(:,k)*(diss_est(:,k))
           enddo
         endif

      end subroutine GFS_stochastics_run

    end module GFS_stochastics
