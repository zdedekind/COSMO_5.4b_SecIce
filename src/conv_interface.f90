!+ Module conv_interface for organizing calls to blocked convection schemes
!------------------------------------------------------------------------------

MODULE conv_interface

!------------------------------------------------------------------------------
!
! Description:
!
! This module initializes the copy-to-block mechanism and organizes the calls to
! the blocked versions of the convection schemes:
!   - Tiedtke convection (used for coarser resolutions)
!   - Tiedtke-Bechtold   (new IFS scheme also used in ICON)
!   - shallow convection (used for convection permitting resolutions)
!
!
! Current Code Owner: MeteoSwiss, Xavier Lapillonne
!  phone: +41 58 460 9237
!  fax  : +41 58 460 9278
!  email:  xavier.lapillonne@meteoswiss.ch
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V5_4b        2016-07-12 Katherine Osterried, Xavier Lapillonne, 
!                         Jochen Foerstner, Ulrich Schaettler
!  Initial release for blocked version of convection in COSMO
! @VERSION@    @DATE@     Ulrich Schaettler
!  Removed test print outs
!  Fixed some syntax errors in the GPU part
!  Variable zcutke_b needs to be treated only in the Tiedtke part
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:
USE kind_parameters,    ONLY :    &
  wp           ! KIND-type parameter for real variables

USE conv_data,          ONLY :    &
#ifdef _OPENACC
  conv_wkarr_alloc, conv_wkarr_dealloc,                 &
#endif
  ctmelt , cc2    , c5hlccp, c5hlscp, chlcdcp, chlsdcp, &
  cu_frr,       & ! fraction of grid area covered with precipitation
  cu_evap,      & ! factor for evaporation of rain
  entr_sc,      & ! mean entrainment rate for shallow convection
  thick_sc,     & ! limit for convective clouds to be "shallow" (in Pa)
  ntrcr_con, itrcr_con, trcr_tend_b, trcr_tran_b,       &
  ! for Tiedtke-Bechtold scheme
  phy_params

USE data_block_fields,  ONLY:     &
  t_b, qv_b, p0_b, ps_b, pp_b, p0hl_b, hhl_b,  &
  dqvdt_conv_b,clc_con_b,clw_con_b,bas_con_b,  &
  top_con_b,tt_conv_b,qvt_conv_b,tket_conv_b,  &
  t_g_b, shfl_s_b, lhfl_s_b,  mflx_con_b,      &
  cape_con_b, qcvg_con_b, fif_b, llandmask_b,  &
  fih_b, pf_b, ph_b,ttdiab_conv_b, prne_con_b, &
  prr_con_b, prs_con_b, qct_conv_b, qit_conv_b,&
  qhfl_b, tke_con_b, u_m_b, v_m_b, ut_conv_b,  &
  vt_conv_b, vgust_con_b, w_conv_b, tke_b,     &
  tt_lheat_new_b, llakemask_b, ttens_conv_b,   &
  cloud_num_b, mtn_mask_b, trop_mask_b,        &
  conv_k850_b, conv_k950_b, qc_b, qi_b, w_b,   &
  rho_b


USE data_constants,     ONLY:     &
  r_earth,      & ! radius of the earth
  pi,           & ! circle constant
  g,            & ! gravity constant
  r_d,          & ! gas constant for dry air
  rdv,          & ! r_d / r_v
  rvd_m_o,      & ! r_v/r_d - 1
  cp_d,         & ! specific heat of dry air at constant pressure
  lh_v,         & ! latent heat of vapourization
  lh_f,         & ! latent heat of fusion
  lh_s,         & ! latent heat of sublimation
  g,            & ! acceleration due to gravity
  b1,           & ! variables for computing the saturation vapour pressure
  b2w,          & ! over water (w) and ice (i)
  b2i,          & !               -- " --
  b3,           & !               -- " --
  b4w,          & !               -- " --
  b4i             !               -- " --

USE data_fields,        ONLY: dqvdt, dqvdt_conv, rmy, w, qvsflx, w_conv, qhfl, &
                              ttens_conv, ttens, sohr, thhr, tt_sso, p0, rlat, &
                              crlat, hhl, llandmask

#ifdef NUDGING
USE data_lheat_nudge,   ONLY: llhn   ! on/off switch for latent heat nudging (lhn)
#endif

USE data_modelconfig,   ONLY: dt, dt2, istart, iend, jstart, jend, ke1, ke,    &
                              edadlat, dlon, dlat, ie, je, je_tot, raddeg

USE data_parallel,      ONLY: my_cart_id

USE data_runcontrol ,   ONLY: nnow, nold, nnew, ntke, loutput_diab, lconf_avg, &
                              lconv_inst, l2tls, nblock, nproma, nlastproma,   &
                              itype_conv, nincconv, y_conv_closure, lctke

USE src_block_fields,   ONLY: CopylistStruct, copyToBlockF, copyFromBlockF,    &
                              register_copy, init_copy_list,                   &
                              mind_ilon, mind_jlat

USE src_tracer,         ONLY: trcr_get, trcr_get_ntrcr, trcr_meta_get,         &
                              trcr_errorstr

USE data_tracer,        ONLY: T_CONV_ID, T_CONV_ON

USE environment,        ONLY: model_abort

USE vgrid_refatm_utils, ONLY: vcoord

USE conv_shallow,       ONLY: cu_shallow
USE conv_tiedtke,       ONLY: cu_tiedtke

USE conv_cuparameters,  ONLY: sucst, su_yoethf, sucumf,       &
                              suphli, suvdf, suvdfs, sucldp
USE conv_cumaster,      ONLY: cumastrn

!==============================================================================

IMPLICIT NONE

!==============================================================================

PRIVATE

! Declare public entities
PUBLIC ::  &
  conv_init,      &
  conv_init_copy, &
  conv_prepare,   &
  conv_organize,  &
  conv_finalize,  &
#ifdef _OPENACC
  conv_in_wkarr_alloc, conv_in_wkarr_dealloc, &
#endif
  convCopyList

!==============================================================================

! Module variables
TYPE(CopylistStruct) :: convCopyList

! Local variables for the interface
!----------------------------------

#ifdef _OPENACC
! Declarations of working  arrays for the convection interface
! (for Tiedtke and shallow)
REAL (KIND=wp), ALLOCATABLE :: &
  zcumflx_b     (:)  ,     & ! cu_base massflux
  zcucape_b     (:)  ,     & ! convective available energy
  zcuconv_b     (:)          ! moisture convergence to be used in cu_base massflux

REAL (KIND=wp), ALLOCATABLE :: &
  zmfu_b        (:,:),     & !
  zdt_con_b     (:,:),     & ! convective tendency of temperature
  zdqv_con_b    (:,:),     & ! convective tendency of specific humidity
  zdtke_con_b   (:,:),     & ! convective tendency of cloud water
  zco_clw_b     (:,:),     & ! convective cloud water
  zdtdiab_con_b (:,:)        ! conv T-tend. Due to LH Exchanges

INTEGER,        ALLOCATABLE :: &
  mbas_con_b    (:)  ,     & ! cloud base level index
  mtop_con_b    (:)  ,     & ! cloud top  level index
  ktype_b       (:)          !

LOGICAL,        ALLOCATABLE :: &
  locum_b       (:),       & ! indicator for convection at gridpoints
  ldshcv_b      (:)

! and only for Tiedtke
REAL (KIND=wp), ALLOCATABLE :: &
  zcutke_b      (:)  ,     & ! convective turbulent energy
  zvddraf_b     (:)  ,     & ! maximum possible convective gust
  zpr_con_b     (:)  ,     & ! convective precipitation rate of rain
  zps_con_b     (:)  ,     & ! convective precipitation rate of snow
  zpre_con_b    (:)          ! convective precipitation rate without evaporation

REAL (KIND=wp), ALLOCATABLE :: &
  zqct_con_b    (:,:),     & ! convective tendency of specific humidity
  zqit_con_b    (:,:),     & ! convective tendency of cloud water
  zdu_con_b     (:,:),     & ! convective tendency of u
  zdv_con_b     (:,:),     & ! convective tendency of v
  ztke_b        (:,:),     & ! to store tke_b for lctke
  zdtlhn_con_b  (:,:)        ! conv T-tend. due to LH exchanges/no rain evap.

! and for Tiedtke-Bechtold
REAL (KIND=wp), ALLOCATABLE :: &
  zplitot_b     (:,:),     &
  zomega_b      (:,:),     &
  zqhfl_b       (:,:),     &
  zshfl_b       (:,:),     &
  zdpres_b      (:,:),     &
  zdgeopot_b    (:,:),     &
  ztt_conv_b    (:,:),     &
  zqvt_conv_b   (:,:),     &
  zprr_con_b    (:,:),     &
  zprs_con_b    (:,:),     &
  zcon_udd_b    (:,:,:)

#endif

!==============================================================================
! Module Procedures in conv_shallow_interface
!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure "conv_shallow_init" to initialize calls to shallow conv.
!------------------------------------------------------------------------------

SUBROUTINE conv_init

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine handles the initialization of the convection schemes.
!
! Method:
!
!------------------------------------------------------------------------------

INTEGER :: i,j,k, izerror, iztrcr, ib, iv, ivend, ii, jj

 REAL(KIND=wp) :: hag               ! height above ground
 REAL(KIND=wp) :: zlat              ! absolute value of latitude in degrees
 REAL(KIND=wp) :: z_rsltn           ! mean approximate resolution in meters
 REAL(KIND=wp) :: zn_ij, zn_avg, zn_sq, zn_min, zn_rms

 ! (for calculation of) reference pressure
 REAL(KIND=wp) :: zpref(ke)

! height of 850 and 950hPa surface for US standard atmosphere in m
! For derivation, see documentation of US standard atmosphere
REAL(KIND=wp), PARAMETER ::          &
  h850_standard = 1457.235199_wp,  &
  h950_standard =  540.3130233_wp

CHARACTER (LEN=255) :: yzerrmsg 
CHARACTER (LEN=  9) :: yzroutine='conv_init'

!------------------------------------------------------------------------------
! Begin Subroutine conv_init
!------------------------------------------------------------------------------ 

  izerror = 0

  ! Set utilitiy constants for the convection schemes
  ctmelt  = 273.16_wp                ! tripel point
  cc2     = b1 * rdv                 ! definition of utitility constants
  c5hlccp = b2w*(b3-b4w)*lh_v/cp_d   ! for saturation humidity
  c5hlscp = b2i*(b3-b4i)*lh_s/cp_d   !
  chlcdcp = lh_v/cp_d                !
  chlsdcp = lh_s/cp_d                !

  ! Set fraction of grid area covered with precipitation
  cu_frr  = 2.0_wp / SQRT ( 0.001_wp/edadlat )
  cu_frr  = MIN ( cu_frr, 1.0_wp )

  ! factor for evaporation of rain
  ALLOCATE (cu_evap(ke))
  cu_evap(:) = 0.0_wp
  DO  k = 1, ke
    cu_evap(k)   = 1.93E-6_wp * 261.0_wp * SQRT( 1000.0_wp/   &
          (38.3_wp*0.293_wp)* SQRT(vcoord%sigm_coord(k)) )*0.5_wp / g
  ENDDO

  ! Check which tracers should be transported by (Tiedtke) convection
  ! US: taken from organize_conv_tiedtke

  IF (itype_conv == 0) THEN
#ifndef MESSY
    ! Get the value (ON/OFF) of the convection metadata for all tracers
    ALLOCATE (itrcr_con(trcr_get_ntrcr()), STAT=izerror)
    itrcr_con = 0

    CALL trcr_meta_get(izerror, T_CONV_ID, itrcr_con)
    IF (izerror /= 0) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF

    ! Find the tracers which should be transported by convection
    ntrcr_con = 0
    DO iztrcr = 1, trcr_get_ntrcr() ! loop over tracers
      IF (itrcr_con(iztrcr) == T_CONV_ON) THEN
        ntrcr_con = ntrcr_con + 1
        itrcr_con(ntrcr_con) = iztrcr
      ENDIF
    ENDDO

    ! Allocate local structure containing the convective tendency of the tracers
    ALLOCATE( trcr_tend_b(nproma, ke, ntrcr_con), STAT=izerror)
    trcr_tend_b = 0.0_wp

    ALLOCATE( trcr_tran_b    (nproma, ke, ntrcr_con), STAT=izerror)
    trcr_tran_b = 0.0_wp
#else
    ntrcr_con = 1
#endif
  ENDIF

#ifdef _OPENACC
  CALL conv_wkarr_alloc (nproma, ke, itype_conv, izerror)
  CALL conv_in_wkarr_alloc (nproma, ke)
#endif

  IF (itype_conv == 2) THEN
    ! Initializations for Tiedtke-Bechtold convection scheme
#ifdef NUDGING
    IF (my_cart_id == 0) THEN
      IF (llhn) THEN
        PRINT *, '*** WARNING: Tiedtke-Bechtold IFS-/ICON- convection scheme may  ***'
        PRINT *, '***          conflict with latent heat nudging!                 ***'
      ENDIF
    ENDIF
#endif

    ! Call the initialization routines of the Tiedtke-Bechtold
    ! IFS-/ICON-convection parametrization

    CALL sucst(54,0,0,0)    ! 54: logical unit for output
                            !     date in the form AAAAMMDD
                            !     number of seconds in the day
                            !     printing level
    CALL su_yoethf


    ! mean pressure of each model level pref(:) (used for NJKT* in sucumf)
    ! US: this should be done in src_input or whereever we have p0 as a full field????
!US DO k = 1, ke
!US   DO j = 1, je
!US     z_p0_sum (j,k) = SUM( p0(:,j,k) ) / REAL(ie,wp)
!US   ENDDO
!US   zpref(k) = SUM( z_p0_sum(:,k) ) / REAL(je,wp)
!US ENDDO

    DO k = 1, ke
      zpref(k) = SUM( p0(:,:,k) ) / REAL(ie*je,wp)
    ENDDO

    ! Please take care for scale-dependent initializations!
    ! Calculate mean approximate model resolution (m)
    ! (WARNING: assumes a small variance of resolution over the domain):
!US z_rsltn = 2.0_wp*pi * r_earth  &
!US       * SQRT( (dlon/360.0_wp)*(dlat/360.0_wp) * crlat(INT(je/2.0_wp),1) )
!US ! make sure z_rsltn (and thereby RTAU) is identical for all processors:
!US IF (num_compute > 1) THEN
!US   z_rsltn = z_rsltn / num_compute
!US   CALL global_values(z_rsltn,1,'SUM',imp_reals,icomm_cart,-1,yzerrmsg,izerror)
!US ENDIF
    z_rsltn = 2.0_wp*pi * r_earth  &
          * SQRT( (dlon/360.0_wp)*(dlat/360.0_wp) * crlat(INT(je_tot/2.0_wp),1) )

    ! sucumf sets (among other things) the namelist variables into phy_params
    ! the last parameter indicates, whether only shallow convection should run
    ! to choose this is not yet an option in the COSMO-Model!
    CALL sucumf (z_rsltn, ke, zpref, phy_params, lshallow_only=.FALSE.)
    CALL suphli
    CALL suvdf
    CALL suvdfs
    CALL sucldp

    ! Initialize fields conv_k850 and conv_k950, which are required for
    ! computing the convective contribution to wind gusts
    conv_k850_b(:,:) = ke
    conv_k950_b(:,:) = ke

    DO ib = 1, nblock
      IF (ib==nblock) THEN
        ivend=nlastproma    !last block
      ELSE
        ivend=nproma        !other blocks
      END IF

      DO iv = 1, ivend
        i = mind_ilon(iv,ib)
        j = mind_jlat(iv,ib)

        DO k = ke, 1, -1
          ! height above ground
          hag = 0.5_wp*( hhl(i,j,k) + hhl(i,j,k+1) ) - hhl(i,j,ke1)

          IF (hag < h850_standard) THEN
            conv_k850_b(iv,ib) = k
          ENDIF
          IF (hag < h950_standard) THEN
            conv_k950_b(iv,ib) = k
          ENDIF
        ENDDO
        ! security measure
        conv_k950_b(iv,ib) = MAX(conv_k950_b(iv,ib),2)
        conv_k850_b(iv,ib) = MAX(conv_k850_b(iv,ib),2)
      ENDDO
    ENDDO

    ! Initialize mask field to distinguish between tropics and extratropics
    ! (needed for some tuning measures)
    DO ib = 1, nblock
      IF (ib==nblock) THEN
        ivend=nlastproma    !last block
      ELSE
        ivend=nproma        !other blocks
      END IF

      DO iv = 1, ivend
        i = mind_ilon(iv,ib)
        j = mind_jlat(iv,ib)

        zlat = ABS(rlat(i,j)*raddeg)
        IF (zlat < 25.0_wp) THEN
          trop_mask_b(iv,ib) = 1.0_wp
        ELSE IF (zlat > 30.0_wp) THEN
          trop_mask_b(iv,ib) = 0.0_wp
        ELSE
          trop_mask_b(iv,ib) = (30.0_wp-zlat)/5.0_wp
        ENDIF
      ENDDO
    ENDDO

    ! compute mask field for mountain or upper slope grid points
    DO ib = 1, nblock
      IF (ib==nblock) THEN
        ivend=nlastproma    !last block
      ELSE
        ivend=nproma        !other blocks
      END IF

      DO iv = 1, ivend
        i = mind_ilon(iv,ib)
        j = mind_jlat(iv,ib)

        IF (i==1 .OR. i==ie .OR. j==1 .OR. j==je) THEN
            mtn_mask_b(iv,ib) = 0.0_wp
        ELSE
          zn_ij  =   hhl(i,j,ke1)
          zn_avg = - zn_ij
          zn_sq  = - zn_ij*zn_ij
          zn_min =   zn_ij

          DO ii = i-1, i+1
            DO jj = j-1, j+1
              zn_avg = zn_avg + hhl(ii,jj,ke1)
              zn_sq  = zn_sq  + hhl(ii,jj,ke1)*hhl(ii,jj,ke1)
              zn_min = MIN(zn_min, hhl(ii,jj,ke1))
            ENDDO
          ENDDO

          zn_avg = zn_avg / 8.0_wp
          zn_sq  = zn_sq  / 8.0_wp
          zn_rms = SQRT(MAX(0.0_wp, zn_sq - zn_avg*zn_avg))
          IF ( zn_ij-zn_min > 200.0_wp .AND. zn_ij > zn_avg ) THEN
            mtn_mask_b(iv,ib) = MIN(1.0_wp, 1.0E5_wp * MAX(0.0_wp, MAX(zn_rms,zn_ij-zn_avg)-200.0_wp) &
                                                            / (z_rsltn*z_rsltn))
          ELSE
            mtn_mask_b(iv,ib) = 0.0_wp
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    ! initialize conv_cloud_num
    DO ib = 1, nblock
      IF (ib==nblock) THEN
        ivend=nlastproma    !last block
      ELSE
        ivend=nproma        !other blocks
      END IF

      DO iv = 1, ivend
        i = mind_ilon(iv,ib)
        j = mind_jlat(iv,ib)

        IF (llandmask(i,j)) THEN
          cloud_num_b(:,:) = 200.0E+06_wp
        ELSE
          cloud_num_b(:,:) =  50.0E+06_wp
        ENDIF
      ENDDO
    ENDDO

  ENDIF  ! itype_conv=2: Tiedtke-Bechtold initializations

!------------------------------------------------------------------------------
! End of module procedure conv_init
!------------------------------------------------------------------------------

END SUBROUTINE conv_init

!==============================================================================
!==============================================================================
!+ Module procedure "conv_init_copy" to register blocks for convection
!------------------------------------------------------------------------------

SUBROUTINE conv_init_copy

!------------------------------------------------------------------------------
!
! Description:
!
! This subroutine handles the registration of the block fields required for
! the convection schemes.
!
! Method: Calls to register_copy
!
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine conv_init_copy
!------------------------------------------------------------------------------

  !Initialize the list of fields to be copied to/from the block structure
  CALL init_copy_list(convCopyList)

  !Register the fields

  ! Variables with intent IN for all convection schemes (Tiedtke / Tiedtke-Bechtold and shallow)
  CALL register_copy(hhl_b          ,convCopyList,copyToBlockF)
  CALL register_copy(t_b            ,convCopyList,copyToBlockF)
  CALL register_copy(qv_b           ,convCopyList,copyToBlockF)
  CALL register_copy(p0_b           ,convCopyList,copyToBlockF)
  CALL register_copy(pp_b           ,convCopyList,copyToBlockF)
  CALL register_copy(p0hl_b         ,convCopyList,copyToBlockF)
  CALL register_copy(ps_b           ,convCopyList,copyToBlockF)
  CALL register_copy(dqvdt_conv_b   ,convCopyList,copyToBlockF)
  CALL register_copy(t_g_b          ,convCopyList,copyToBlockF)
  CALL register_copy(shfl_s_b       ,convCopyList,copyToBlockF)
  CALL register_copy(lhfl_s_b       ,convCopyList,copyToBlockF)

  IF (itype_conv == 0 .OR. itype_conv == 2) THEN    ! (only for Tiedtke / Tiedtke-Bechtold)
  CALL register_copy(u_m_b          ,convCopyList,copyToBlockF)
  CALL register_copy(v_m_b          ,convCopyList,copyToBlockF)
  CALL register_copy(llandmask_b    ,convCopyList,copyToBlockF)
  ENDIF

  IF (itype_conv == 0) THEN    ! (only for Tiedtke)
  CALL register_copy(w_conv_b       ,convCopyList,copyToBlockF)
  CALL register_copy(qhfl_b         ,convCopyList,copyToBlockF)
  ENDIF

  IF (itype_conv == 2) THEN    ! (only for Tiedtke-Bechtold)
  CALL register_copy(qc_b           ,convCopyList,copyToBlockF)
  CALL register_copy(qi_b           ,convCopyList,copyToBlockF)
  CALL register_copy(llakemask_b    ,convCopyList,copyToBlockF)
  CALL register_copy(rho_b          ,convCopyList,copyToBlockF)
  CALL register_copy(w_b            ,convCopyList,copyToBlockF)
  CALL register_copy(ttens_conv_b   ,convCopyList,copyToBlockF)
  ENDIF

  !Variables with intent INOUT (for all schemes)
  CALL register_copy(bas_con_b      ,convCopyList,copyToBlockF)
  CALL register_copy(bas_con_b      ,convCopyList,copyFromBlockF)
  CALL register_copy(top_con_b      ,convCopyList,copyToBlockF)
  CALL register_copy(top_con_b      ,convCopyList,copyFromBlockF)
  CALL register_copy(tket_conv_b    ,convCopyList,copyToBlockF)
  CALL register_copy(tket_conv_b    ,convCopyList,copyFromBlockF)
  IF (itype_conv == 0 .OR. itype_conv == 2) THEN    ! (only for Tiedtke / Tiedtke-Bechtold)
  CALL register_copy(vgust_con_b    ,convCopyList,copyToBlockF)
  CALL register_copy(vgust_con_b    ,convCopyList,copyFromBlockF)
  ENDIF

  ! latent heat nudging only to be used in Tiedtke scheme
  IF (itype_conv == 0 .AND. llhn) THEN
  CALL register_copy(tt_lheat_new_b ,convCopyList,copyToBlockF)
  CALL register_copy(tt_lheat_new_b ,convCopyList,copyFromBlockF)
  ENDIF

  !Variables with intent OUT (all schemes)
  CALL register_copy(clc_con_b      ,convCopyList,copyFromBlockF)
  CALL register_copy(clw_con_b      ,convCopyList,copyFromBlockF)
  CALL register_copy(tt_conv_b      ,convCopyList,copyFromBlockF)
  CALL register_copy(ttdiab_conv_b  ,convCopyList,copyFromBlockF)
  CALL register_copy(qvt_conv_b     ,convCopyList,copyFromBlockF)
  CALL register_copy(mflx_con_b     ,convCopyList,copyFromBlockF)
  CALL register_copy(cape_con_b     ,convCopyList,copyFromBlockF)
  CALL register_copy(qcvg_con_b     ,convCopyList,copyFromBlockF)

  IF (itype_conv == 0 .OR. itype_conv == 2) THEN    ! (only for Tiedtke / Tiedtke-Bechtold)
  CALL register_copy(prr_con_b      ,convCopyList,copyFromBlockF)
  CALL register_copy(prs_con_b      ,convCopyList,copyFromBlockF)
  CALL register_copy(prne_con_b     ,convCopyList,copyFromBlockF)
  CALL register_copy(qct_conv_b     ,convCopyList,copyFromBlockF)
  CALL register_copy(qit_conv_b     ,convCopyList,copyFromBlockF)
  CALL register_copy( ut_conv_b     ,convCopyList,copyFromBlockF)
  CALL register_copy( vt_conv_b     ,convCopyList,copyFromBlockF)
  CALL register_copy(tke_con_b      ,convCopyList,copyFromBlockF)
  ENDIF

!------------------------------------------------------------------------------
! End of module procedure conv_shallow_init_copy
!------------------------------------------------------------------------------

END SUBROUTINE conv_init_copy

!==============================================================================
!==============================================================================
!+ Module procedure "conv_shallow_prepare" to compute a transformed dqvdt
!------------------------------------------------------------------------------

SUBROUTINE conv_prepare

!------------------------------------------------------------------------------
!
! Description:
!
! This subroutine computes the transformed dqvdt needed as input for the
! shallow convection.
!
! Method:
!
!
!------------------------------------------------------------------------------

 IMPLICIT NONE

! Subroutine arguments:
! --------------------


! Local parameters:
! ----------------
REAL    (KIND=wp   ) ::  &  ! Weights for horizontal averaging
  zcent = 0.2500_wp, & ! centre weight in a nine point stencil
  zside = 0.1250_wp, & ! weight for side points
  zedge = 0.0625_wp    ! weight for edge points

INTEGER :: i,j,k, nx

!------------------------------------------------------------------------------
! Begin Subroutine conv_prepare
!------------------------------------------------------------------------------

  !GPU_TMP_UPDATE : this update should be remove once the physics run fully on GPU
#ifdef _OPENACC
  print*, "GPU-info : update device before conv_shallow_prepare"
#endif
  !!!!US  !$acc update device(dqvdt_conv, dqvdt, rmy, w_conv, w, qhfl, qvsflx)
  !!!!US  why dqvdt_conv ?? it is computed here
  !$acc update device(rmy, dqvdt, w, qvsflx, ttens_conv, ttens, sohr, thhr, tt_sso)

! Prepare input arrays for the convection scheme
  !$acc data &
  !$acc present (rmy, dqvdt_conv, dqvdt, w_conv, w, qhfl, qvsflx)        &
  !$acc present (ttens_conv, ttens, sohr, thhr, tt_sso)

  ! Select timelevel and timestep of the computation
  IF ( l2tls ) THEN
    nx  = nnow
  ELSE
    nx  = nold
  ENDIF

  IF (lconf_avg) THEN ! Replace convective forcings by area average

    !$acc parallel
    DO  k = 1, ke
      !$acc loop gang
      DO j = jstart, jend
        !$acc loop vector
        DO i = istart, iend
          dqvdt_conv(i,j,k)  = (1.0_wp-rmy(i,j,1)) * ( zcent*dqvdt(i,j,k)     &
                          + zside*( dqvdt(i-1,j,k  ) + dqvdt(i+1,j,k  )       &
                                  + dqvdt(i,j-1,k  ) + dqvdt(i,j+1,k  )  )    &
                          + zedge*( dqvdt(i-1,j-1,k) + dqvdt(i+1,j-1,k)       &
                                  + dqvdt(i-1,j+1,k) + dqvdt(i+1,j+1,k)  ) )
        ENDDO
      ENDDO
    ENDDO
    !$acc end parallel

    IF (itype_conv == 0) THEN
      ! additional preparations necessary for Tiedtke convection
      !$acc parallel
      DO  k = 1, ke
        !$acc loop gang
        DO j = jstart, jend
          !$acc loop vector
          DO i = istart, iend
            w_conv(i,j,k) = 0.5_wp * &
                           (  zcent*( w(i,j,k,nx) + w(i,j,k+1,nx)  )           &
                            + zside*( w(i-1,j  ,k  ,nx) + w(i+1,j  ,k  ,nx)    &
                                    + w(i,j-1  ,k  ,nx) + w(i  ,j+1,k  ,nx)  ) &
                            + zedge*( w(i-1,j-1,k  ,nx) + w(i+1,j-1,k  ,nx)    &
                                    + w(i-1,j+1,k  ,nx) + w(i+1,j+1,k  ,nx)  ) &
                            + zside*( w(i-1,j  ,k+1,nx) + w(i+1,j  ,k+1,nx)    &
                                    + w(i  ,j-1,k+1,nx) + w(i  ,j+1,k+1,nx)  ) &
                            + zedge*( w(i-1,j-1,k+1,nx) + w(i+1,j-1,k+1,nx)    &
                                    + w(i-1,j+1,k+1,nx) + w(i+1,j+1,k+1,nx)  ) )
          ENDDO
        ENDDO
      ENDDO
      !$acc end parallel

      !$acc parallel
      DO j = jstart, jend
        !$acc loop vector
        DO i = istart, iend
          qhfl (i,j)      = (1.0_wp-rmy(i,j,1)) * ( zcent*qvsflx(i,j)    &
                          + zside*( qvsflx(i-1,j  ) + qvsflx(i+1,j  )    &
                                  + qvsflx(i,j-1  ) + qvsflx(i,j+1  )  ) &
                          + zedge*( qvsflx(i-1,j-1) + qvsflx(i+1,j-1)    &
                                  + qvsflx(i-1,j+1) + qvsflx(i+1,j+1)  ) )
        ENDDO
      ENDDO
      !$acc end parallel
    ENDIF  ! itype_conv==0

  ELSE     ! not lconf_avg

    !$acc parallel
    DO  k = 1, ke
      !$acc loop gang
      DO j = jstart, jend
        !$acc loop vector
        DO i = istart, iend
          dqvdt_conv(i,j,k)  = dqvdt(i,j,k)*(1.0_wp-rmy(i,j,1))
        ENDDO
      ENDDO
    ENDDO
    !$acc end parallel

    IF (itype_conv == 0) THEN
      !$acc parallel
      DO  k = 1, ke
        !$acc loop gang
        DO j = jstart, jend
          !$acc loop vector
          DO i = istart, iend
            w_conv(i,j,k) = 0.5_wp * (w(i,j,k,nx) + w(i,j,k+1,nx))
          ENDDO
        ENDDO
      ENDDO
      !$acc end parallel

      !$acc parallel
      DO j = jstart, jend
        !$acc loop vector
        DO i = istart, iend
          qhfl (i,j)      = (1.0_wp-rmy(i,j,1)) * qvsflx(i,j)
        ENDDO
      ENDDO
      !$acc end parallel
    ENDIF  ! itype_conv==0

  ENDIF    ! lconf_avg

  IF (itype_conv == 2) THEN
    ! sum of temperature tendencies
    !$acc parallel
    DO  k = 1, ke
      !$acc loop gang
      DO j = jstart, jend
        !$acc loop vector
        DO i = istart, iend
          ttens_conv(i,j,k) = ttens_conv(i,j,k)  &
                      + ttens(i,j,k) + sohr(i,j,k) + thhr(i,j,k) + tt_sso(i,j,k)
        ENDDO
      ENDDO
    ENDDO
    !$acc end parallel
  ENDIF
  !$acc end data

  !GPU_TMP_UPDATE : this update should be remove once the physics run fully on GPU
#ifdef _OPENACC
  print*, "GPU-info : update host after conv_shallow_prepare"
#endif
  !!!!US why dqvdt, rmy???  !$acc update host(dqvdt_conv, dqvdt, rmy)
  !$acc update host(dqvdt_conv, w_conv, qhfl, ttens_conv)

!------------------------------------------------------------------------------
! End of module procedure conv_shallow_prepare
!------------------------------------------------------------------------------

END SUBROUTINE conv_prepare

!==============================================================================
!==============================================================================
!+ Module procedure "conv_shallow_organize" to organize calls to shallow conv.
!------------------------------------------------------------------------------

SUBROUTINE conv_organize (ipend, iblock, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!
! This subroutine organizes the calls to the convection schemes.
!
!
! Method:
!
!
!------------------------------------------------------------------------------

 IMPLICIT NONE

! Subroutine arguments:
! --------------------
INTEGER, INTENT(IN)  ::  &
  ipend,                 & ! compute domain in horizontal direction
  iblock                   ! number of block to work on

INTEGER, INTENT(OUT) ::  &
  ierror

CHARACTER(LEN=*),        INTENT(OUT) ::  &
  yerrmsg

! Local parameters:
! ----------------

INTEGER :: iz, ip,k, km1,nx, mtop, mbas, zk850, zk950
REAL(KIND=wp) :: zdt, zbas, ztop, zdttop, zdtbot, zdepth_max, zdepth_now,     &
                 ticeini, lfocpd, u850, u950,v850, v950, wfac

CHARACTER (LEN= 25)     :: &
  yzroutine

#ifndef _OPENACC
! Declarations of working  arrays for the convection interface
REAL (KIND=wp) :: &
  zcumflx_b     (nproma)   ,    & ! cu_base massflux
  zcucape_b     (nproma)   ,    & ! convective available energy
  zcuconv_b     (nproma)          ! moisture convergence to be used in cu_base massflux

REAL (KIND=wp) :: &
  zmfu_b        (nproma,ke),    & !
  zdt_con_b     (nproma,ke),    & ! convective tendency of temperature
  zdqv_con_b    (nproma,ke),    & ! convective tendency of specific humidity
  zdtke_con_b   (nproma,ke),    & ! convective tendency of cloud water
  zco_clw_b     (nproma,ke),    & ! convective cloud water
  zdtdiab_con_b (nproma,ke)       ! conv T-tend. Due to LH Exchanges

INTEGER        :: &
  mbas_con_b    (nproma)   ,    & ! cloud base level index
  mtop_con_b    (nproma)   ,    & ! cloud top  level index
  ktype_b       (nproma)          ! cloud top  level index

LOGICAL        :: &
  locum_b       (nproma)   ,    & ! indicator for convection at gridpoints
  ldshcv_b      (nproma)

! and only for Tiedtke
REAL (KIND=wp) :: &
  zcutke_b      (nproma)  ,     & ! convective turbulent energy
  zvddraf_b     (nproma)  ,     & ! maximum possible convective gust
  zpr_con_b     (nproma)  ,     & ! convective precipitation rate of rain
  zps_con_b     (nproma)  ,     & ! convective precipitation rate of snow
  zpre_con_b    (nproma)          ! convective precipitation rate without evaporation

REAL (KIND=wp) :: &
  zqct_con_b    (nproma,ke),    & ! convective tendency of specific humidity
  zqit_con_b    (nproma,ke),    & ! convective tendency of cloud water
  zdu_con_b     (nproma,ke),    & ! convective tendency of u
  zdv_con_b     (nproma,ke),    & ! convective tendency of v
  ztke_b        (nproma,ke),    & ! to store tke_b for lctke
  zdtlhn_con_b  (nproma,ke)       ! conv T-tend. due to LH exchanges/no rain evap.

! and for Tiedtke-Bechtold
REAL (KIND=wp) :: &
  zplitot_b     (nproma,ke) ,   &
  zomega_b      (nproma,ke) ,   &
  zqhfl_b       (nproma,ke1),   &
  zshfl_b       (nproma,ke1),   &
  zdpres_b      (nproma,ke) ,   &
  zdgeopot_b    (nproma,ke) ,   &
  ztt_conv_b    (nproma,ke) ,   &
  zqvt_conv_b   (nproma,ke) ,   &
  zprr_con_b    (nproma,ke1),   &
  zprs_con_b    (nproma,ke1),   &
  zcon_udd_b    (nproma,ke,6)
#endif

!------------------------------------------------------------------------------
! Begin Subroutine conv_organize
!------------------------------------------------------------------------------

  ! XA: Update device lines for Shallow Convection CALL
  !GPU_TMP_UPDATE : this update should be remove once the physics run fully on GPU
#ifdef _OPENACC
  print*, "GPU-info : update device before conv_shallow_organize"
#endif
  !$acc update device(fif_b, fih_b, pf_b, ph_b, hhl_b, p0_b, pp_b, p0hl_b, pp_b)
  !$acc update device(ps_b, t_b, qv_b, dqvdt_conv_b, clc_con_b, clw_con_b)
  !$acc update device(bas_con_b, top_con_b, tt_conv_b, qvt_conv_b, tket_conv_b)
  !$acc update device(mflx_con_b, cape_con_b, qcvg_con_b)
  !$acc update device(t_g_b, shfl_s_b, lhfl_s_b)
! for Tiedtke
  !$acc update device(tke_b, qct_conv_b, qit_conv_b, ut_conv_b, vt_conv_b)
  !$acc update device(prr_con_b, prs_con_b, prne_con_b, tt_lheat_new_b)
  !$acc update device(tke_con_b, vgust_con_b, llandmask_b)
! for Tiedtke-Bechtold
  !$acc update device(trop_mask_b, mtn_mask_b, cloud_num_b, conv_k850_b)
  !$acc update device(conv_k950_b, llakemask_b, qc_b, qi_b, ttens_conv_b)

  !$acc data                                                                  &
  !$acc present (fif_b, fih_b, pf_b, ph_b, hhl_b, p0_b, pp_b, p0hl_b, pp_b)   &
  !$acc present (ps_b, t_b, qv_b, dqvdt_conv_b, clc_con_b, clw_con_b)         &
  !$acc present (bas_con_b, top_con_b, tt_conv_b, qvt_conv_b, tket_conv_b)    &
  !$acc present (mflx_con_b, cape_con_b, qcvg_con_b)                          &
  !$acc present (t_g_b, shfl_s_b, lhfl_s_b, ttdiab_conv_b)                    &
  !$acc present (zcumflx_b , zcuconv_b    , zco_clw_b , zmfu_b)               &
  !$acc present (zdt_con_b , zdtdiab_con_b, zdqv_con_b, zdtke_con_b)          &
  !$acc present (mbas_con_b, mtop_con_b   , locum_b)                          &
! for Tiedtke
  !$acc present (tke_b, qct_conv_b, qit_conv_b, ut_conv_b, vt_conv_b)         &
  !$acc present (prr_con_b, prs_con_b, prne_con_b, tt_lheat_new_b)            &
  !$acc present (tke_con_b, vgust_con_b, llandmask_b)                         &
! for Tiedtke-Bechtold
  !$acc present (trop_mask_b, mtn_mask_b, cloud_num_b, conv_k850_b)           &
  !$acc present (conv_k950_b, llakemask_b, qc_b, qi_b, ttens_conv_b)          &
  !$acc present (rho_b, w_b, u_m_b, v_m_b)

  ierror   = 0
  yerrmsg  = ''
  yzroutine= 'conv_organize'

  ! Select timelevel and timestep of the computation
  IF ( l2tls ) THEN
    nx  = nnow
    zdt = dt
  ELSE
    nx  = nold
    zdt = dt2
  ENDIF

  ! Reset the  convective cloud cover each time when the convection scheme is called
  !$acc parallel
  DO k = 1, ke
    !$acc loop gang vector
    DO ip = 1, ipend
      clc_con_b (ip,k) = 0.0_wp
    ENDDO
  ENDDO
  !$acc end parallel

  IF ( lconv_inst ) THEN
    !$acc parallel
    !$acc loop gang vector
    DO ip = 1, ipend
      bas_con_b (ip) = 0.0_wp
      top_con_b (ip) = 0.0_wp
    ENDDO
    !$acc end parallel
  ENDIF

  ! Preset the output arrays of the convection scheme
  !$acc parallel
  DO k = 1, ke
    !$acc loop gang vector
    DO ip = 1, ipend
      zdt_con_b    (ip,k) = 0.0_wp
      zdqv_con_b   (ip,k) = 0.0_wp
      zdtke_con_b  (ip,k) = 0.0_wp
      zco_clw_b    (ip,k) = 0.0_wp
      zdtdiab_con_b(ip,k) = 0.0_wp
    ENDDO
  ENDDO
  !$acc end parallel

  IF (itype_conv == 0) THEN
    ! for Tiedtke also preset ztke for lctke
    IF (lctke) THEN
      !$acc parallel
      DO k = 1, ke
        !$acc loop gang vector
        DO ip = 1, ipend
          ztke_b (ip,k) = tke_b(ip,k,ntke)
        ENDDO
      ENDDO
      !$acc end parallel
    ENDIF

    !$acc parallel
    !$acc loop gang vector
    DO ip = 1, ipend
      zcutke_b   (ip) = 0.0_wp
    ENDDO
    !$acc end parallel

!   ! for Tiedtke also preset zdtlhn_con
!   !$acc parallel
!   DO k = 1, ke
!     !$acc loop gang vector
!     DO ip = 1, ipend
!       zdtlhn_con_b (ip,k) = 0.0_wp
!     ENDDO
!   ENDDO
!   !$acc end parallel
  ENDIF

  !$acc parallel
  !$acc loop gang vector
  DO ip = 1, ipend
    mbas_con_b (ip)   = 0
    mtop_con_b (ip)   = 0
    locum_b    (ip)   = .FALSE.
    ! Preset the additional output fields for convection
    zcumflx_b  (ip) = 0.0_wp
    zcucape_b  (ip) = 0.0_wp
    zcuconv_b  (ip) = 0.0_wp
  ENDDO
  !$acc end parallel

  !$acc parallel
  DO k = 1, ke
    km1 = MAX ( 1, k-1 )
    !$acc loop gang vector
    DO ip = 1, ipend
      fif_b  (ip,k)  = 0.5_wp*g*( hhl_b(ip,k) + hhl_b(ip,k+1) )
      fih_b  (ip,k)  = g* hhl_b(ip,k)
      pf_b   (ip,k)  = p0_b(ip,k) + pp_b(ip,k)
      ph_b   (ip,k)  = p0hl_b(ip,k) + 0.5_wp*(pp_b(ip,k) + pp_b(ip,km1))
    ENDDO
  ENDDO
  !$acc end parallel

  !$acc parallel
  !$acc loop gang vector
  DO ip = 1, ipend
    fih_b  (ip,ke1)  = g* hhl_b(ip,ke1)
    ph_b   (ip,ke1)  = ps_b(ip)
  ENDDO
  !$acc end parallel

  ! More initializations for Tiedtke-Bechtold
  IF (itype_conv == 2) THEN

    lfocpd  = lh_f/cp_d
    ticeini = 261.15_wp

    !$acc parallel
    !$acc loop gang vector
    DO ip = 1, ipend
      ldshcv_b(ip) = .TRUE.
    ENDDO
    !$acc end parallel

    ! cloud water + cloud ice for entrainment computation
    !$acc parallel
    DO k = 1, ke
      !$acc loop gang vector
      DO ip = 1, ipend
        zplitot_b(ip,k)  = qc_b(ip,k) + qi_b(ip,k)
      ENDDO
    ENDDO
    !$acc end parallel

    ! compute omega (vertical velocity p-dot in pressure coordinate system)
    !$acc parallel
    DO k = 1, ke
      !$acc loop gang vector
      DO ip = 1, ipend
        zomega_b(ip,k) = - rho_b(ip,k) * g * 0.5_wp * ( w_b(ip,k) + w_b(ip,k+1) )
      ENDDO
    ENDDO
    !$acc end parallel


    ! compute latent and sensible heat flux forcing (surface fluxes only)
    !$acc parallel
    DO k = 1, ke
      km1 = MAX ( 1, k-1 )
      !$acc loop gang vector
      DO ip = 1, ipend
        zqhfl_b(ip,k)  = 0.0_wp
        zshfl_b(ip,k)  = 0.0_wp
      ENDDO
    ENDDO
    !$acc end parallel

    !$acc parallel
    !$acc loop gang vector
    DO ip = 1, ipend
      zqhfl_b(ip,ke1)  = lhfl_s_b(ip) / lh_v  ! convention in ECMWF IFS: w*q*
      zshfl_b(ip,ke1)  = shfl_s_b(ip)         ! convention in ECMWF IFS: c_p w*T*
    ENDDO
    !$acc end parallel

    ! layer thickness with respect to pressure and geopotential
    !$acc parallel
    DO k = 1, ke
      !$acc loop gang vector
      DO ip = 1, ipend
        zdpres_b(ip,k)   = ph_b(ip,k+1) - ph_b(ip,k)
        zdgeopot_b(ip,k) = g * ( hhl_b(ip,k) - hhl_b(ip,k+1) )
      ENDDO
    ENDDO
    !$acc end parallel

    ! save the tendencies of temperature and moisture
    !$acc parallel
    DO k = 1, ke
      !$acc loop gang vector
      DO ip = 1, ipend
        ztt_conv_b (ip,k) = ttens_conv_b(ip,k)
        zqvt_conv_b(ip,k) = dqvdt_conv_b(ip,k)
      ENDDO
    ENDDO
    !$acc end parallel

    ! The following input fields must be reset to zero because the convective
    ! tendencies are added to them
    !$acc parallel
    DO k = 1, ke
      !$acc loop gang vector
      DO ip = 1, ipend
        ut_conv_b   (ip,k) = 0.0_wp
        vt_conv_b   (ip,k) = 0.0_wp
        qct_conv_b  (ip,k) = 0.0_wp
        qit_conv_b  (ip,k) = 0.0_wp
        tket_conv_b (ip,k) = 0.0_wp   ! neccessary?
        zprr_con_b  (ip,k) = 0.0_wp
        zprs_con_b  (ip,k) = 0.0_wp
      ENDDO
    ENDDO
    !$acc end parallel

  ENDIF  ! itype_conv == 2

  SELECT CASE (itype_conv)

  CASE (0)
    CALL cu_tiedtke (                                          &
         ! input
         t_b         , qv_b        , trcr_tran_b , u_m_b       , v_m_b       , &
         w_conv_b    , ztke_b      , fif_b       , fih_b       , pf_b        , &
         ph_b        , t_g_b       , dqvdt_conv_b, qhfl_b      , llandmask_b , &
         zdt         , nproma      , ke          , 1           , ipend       , &
         ntrcr_con   , iblock      ,                                           &

         ! output
         zcumflx_b   , zcucape_b   , zcuconv_b   , zcutke_b    ,               &
         zvddraf_b   , zco_clw_b   ,                                           &
         zdt_con_b   , zdtlhn_con_b,zdtdiab_con_b, zdqv_con_b  , trcr_tend_b , &
         zdu_con_b   , zdv_con_b   , zqct_con_b  , zqit_con_b  , zdtke_con_b , &
         zpr_con_b   , zps_con_b   , mbas_con_b  , mtop_con_b  , zpre_con_b  , &
         locum_b)

  CASE (2)
    CALL cumastrn &
      ( kidia = 1, kfdia = ipend, klon = nproma,     & !! IN
        ktdia = 1, klev = ke,                        & !! IN
        ldland      = llandmask_b (:),               & !! IN
        ptsphy      = zdt*nincconv,                  & !! IN
        ldlake      = llakemask_b (:),               & !! IN
        k950        = conv_k950_b (:,iblock),        & !! IN
        phy_params  = phy_params,                    & !! IN
        capdcfac    = trop_mask_b (:,iblock),        & !! IN
        mtnmask     = mtn_mask_b  (:,iblock),        & !! IN
        pten        = t_b         (:,:),             & !! IN
        pqen        = qv_b        (:,:),             & !! IN
        puen        = u_m_b       (:,:),             & !! IN
        pven        = v_m_b       (:,:),             & !! IN
        plitot      = zplitot_b   (:,:),             & !! IN
        pvervel     = zomega_b    (:,:),             & !! IN
        pqhfl       = zqhfl_b     (:,:),             & !! IN
        pahfs       = zshfl_b     (:,:),             & !! IN
        pap         = pf_b        (:,:),             & !! IN
        paph        = ph_b        (:,:),             & !! IN
        pgeo        = fif_b       (:,:),             & !! IN
        pgeoh       = fih_b       (:,:),             & !! IN
        zdph        = zdpres_b    (:,:),             & !! IN
        zdgeoh      = zdgeopot_b  (:,:),             & !! IN
        pcloudnum   = cloud_num_b (:,iblock),        & !! IN
        ptent       = ztt_conv_b  (:,:),             & !! INOUT
        ptenu       = ut_conv_b   (:,:),             & !! OUT
        ptenv       = vt_conv_b   (:,:),             & !! OUT
        ptenq       = zqvt_conv_b (:,:),             & !! INOUT
        ptenl       = qct_conv_b  (:,:),             & !! OUT
        pteni       = qit_conv_b  (:,:),             & !! OUT
        ldcum       = locum_b     (:),               & !! OUT
        ktype       = ktype_b     (:),               & !! OUT
        kcbot       = mbas_con_b  (:),               & !! OUT
        kctop       = mtop_con_b  (:),               & !! OUT
        ldshcv      = ldshcv_b    (:),               & !! IN
        pmfu        = zcon_udd_b  (:,:,1),           & !! OUT
        pmfd        = zcon_udd_b  (:,:,2),           & !! OUT
        pmfude_rate = zcon_udd_b  (:,:,3),           & !! OUT
        pmfdde_rate = zcon_udd_b  (:,:,4),           & !! OUT
        ptu         = zcon_udd_b  (:,:,5),           & !! OUT
        pqu         = zcon_udd_b  (:,:,6),           & !! OUT
        plu         = clw_con_b   (:,:),             & !! OUT
        pmflxr      = zprr_con_b  (:,:),             & !! OUT
        pmflxs      = zprs_con_b  (:,:),             & !! OUT
        prain       = prne_con_b  (:),               & !! OUT
        pdtke_con   = tket_conv_b (:,:),             & !! OUT
        pcape       = cape_con_b  (:),               & !! OUT
        ktrac       = 0 )                              !! IN 


  CASE (3)
    CALL cu_shallow (                                          &
         t_b       , qv_b         ,                            &
         fif_b     , fih_b        , pf_b      , ph_b        ,  &
         t_g_b     , shfl_s_b     , lhfl_s_b  , dqvdt_conv_b,  &
         y_conv_closure           , thick_sc  , entr_sc     ,  &
         nproma    , ke           , 1         , ipend       ,  &
         zcumflx_b , zcuconv_b    , zco_clw_b , zmfu_b      ,  &
         zdt_con_b , zdtdiab_con_b, zdqv_con_b, zdtke_con_b ,  &
         mbas_con_b, mtop_con_b   , locum_b)

  CASE DEFAULT
    PRINT *, '*** ERROR: type of scheme not implemented:  ', itype_conv
    ierror = 10
  END SELECT

  ! Store the output from the convection scheme on the corresponding global arrays

  ! for Tiedtke and for shallow convection
  ! --------------------------------------

  IF (itype_conv == 0 .OR. itype_conv == 3) THEN
    !$acc parallel
    DO  k = 1, ke
      !$acc loop gang vector
      DO  ip = 1, ipend
        IF( locum_b(ip) ) THEN
          tt_conv_b  (ip,k) = zdt_con_b (ip,k)
          qvt_conv_b (ip,k) = zdqv_con_b(ip,k)
          tket_conv_b(ip,k) = zdtke_con_b(ip,k)
          clw_con_b  (ip,k) = zco_clw_b (ip,k)
        ELSE
          tt_conv_b  (ip,k) = 0.0_wp
          qvt_conv_b (ip,k) = 0.0_wp
          tket_conv_b(ip,k) = 0.0_wp
          clw_con_b  (ip,k) = 0.0_wp
        ENDIF
      ENDDO
    ENDDO
    !$acc end parallel

    ! for Tiedtke and for shallow convection
    IF (loutput_diab) THEN
      !$acc parallel
      !$acc loop gang
      DO  k = 1, ke
        !$acc loop vector
        DO ip = 1, ipend
         IF (locum_b(ip)) THEN
           ttdiab_conv_b(ip,k) = zdtdiab_con_b (ip,k)
         ELSE
           ttdiab_conv_b(ip,k) = 0.0_wp
         ENDIF
        ENDDO
      ENDDO
    !$acc end parallel
    ENDIF

    ! for Tiedtke and for shallow convection
    !$acc parallel
    !$acc loop gang vector
    DO  ip = 1, ipend
      IF( locum_b(ip) .AND. mtop_con_b(ip) > 0 ) THEN
        mtop_con_b (ip) = MAX   ( mtop_con_b(ip)-2, 2 )
        zdepth_max      = bas_con_b(ip) - top_con_b(ip)
        zdepth_now      = REAL ( mbas_con_b(ip) - mtop_con_b(ip), wp )
        IF ( zdepth_max < zdepth_now ) THEN
          top_con_b(ip) = REAL ( mtop_con_b(ip), wp )
          bas_con_b(ip) = REAL ( mbas_con_b(ip), wp )
        ENDIF
      ENDIF
    ENDDO
    !$acc end parallel

    !$acc parallel
    !$acc loop gang vector
    DO ip = 1, ipend
      IF (locum_b(ip)) THEN
        mflx_con_b(ip) = zcumflx_b(ip)
        cape_con_b(ip) = zcucape_b(ip)
        qcvg_con_b(ip) = zcuconv_b(ip)
      ELSE
        mflx_con_b(ip) = 0.0_wp
        cape_con_b(ip) = 0.0_wp
        qcvg_con_b(ip) = 0.0_wp
      END IF
    END DO
    !$acc end parallel
  ENDIF

  ! only for Tiedtke convection
  ! ---------------------------

  IF (itype_conv == 0) THEN
    !$acc parallel
    DO  k = 1, ke
      !$acc loop gang vector
      DO  ip = 1, ipend
        IF( locum_b(ip) ) THEN
          qct_conv_b (ip,k) = zqct_con_b(ip,k)
          qit_conv_b (ip,k) = zqit_con_b(ip,k)
           ut_conv_b (ip,k) =  zdu_con_b(ip,k)
           vt_conv_b (ip,k) =  zdv_con_b(ip,k)
        ELSE
          qct_conv_b (ip,k) = 0.0_wp
          qit_conv_b (ip,k) = 0.0_wp
           ut_conv_b (ip,k) = 0.0_wp
           vt_conv_b (ip,k) = 0.0_wp
        ENDIF
      ENDDO
    ENDDO
    !$acc end parallel

    ! only for Tiedtke convection
    !  US: in principle we could use prr_con_b, prs_con_b, prne_con_b directly as
    !      arguments for cu_tiedtke, but their setting in the Tiedtke convection is
    !      not related to locum. So this could change results
    !$acc parallel
    !$acc loop gang vector
    DO  ip = 1, ipend
      IF( locum_b(ip) .AND. mtop_con_b(ip) > 0 ) THEN
        prr_con_b (ip) = zpr_con_b(ip)
        prs_con_b (ip) = zps_con_b(ip)
        prne_con_b(ip) = zpre_con_b(ip)
      ELSE
        prr_con_b (ip) = 0.0_wp
        prs_con_b (ip) = 0.0_wp
        prne_con_b(ip) = 0.0_wp
      ENDIF
    ENDDO
    !$acc end parallel

#ifdef NUDGING
    IF (llhn) THEN
      !$acc parallel
      !$acc loop gang
      DO  k = 1, ke
        !$acc loop vector
        DO ip = 1, ipend
         IF (locum_b(ip)) THEN
           tt_lheat_new_b(ip,k) = tt_lheat_new_b(ip,k) + zdt*zdtlhn_con_b(ip,k)
         ENDIF
        ENDDO
      ENDDO
      !$acc end parallel
    ENDIF
#endif

    !$acc parallel
    !$acc loop gang vector
    DO ip = 1, ipend
      IF (locum_b(ip)) THEN
        tke_con_b (ip) = MAX (0.0_wp, zcutke_b(ip))
        ! correction of zvddraf_b (max. possible convective gust)
        IF (3600.0_wp*(zpr_con_b(ip)+zps_con_b(ip)) <= 0.015_wp) zvddraf_b(ip)=0.0_wp
        vgust_con_b(ip) = MAX(vgust_con_b(ip),zvddraf_b(ip))
      ELSE
        tke_con_b (ip) = 0.0_wp
      END IF
    END DO
    !$acc end parallel
  ENDIF   ! itype_conv == 0

  ! only for Tiedtke-Bechtold convection
  ! ------------------------------------

  IF (itype_conv == 2) THEN
    ! temperature tendency due to convection
    !$acc parallel
    !$acc loop gang
    DO  k = 1, ke
      !$acc loop vector
      DO ip = 1, ipend
        tt_conv_b(ip,k)  =  ztt_conv_b(ip,k) - ttens_conv_b(ip,k)

        ! moisture convection due to convection
        qvt_conv_b(ip,k) = zqvt_conv_b(ip,k) - dqvdt_conv_b(ip,k)
      ENDDO
    ENDDO
    !$acc end parallel

    ! Convert detrained cloud ice into cloud water if the temperature is only slightly below freezing
    ! and convective cloud top is not cold enough for substantial ice initiation
    !$acc parallel
    !$acc loop gang
    DO  k = 1, ke
      !$acc loop vector
      DO ip = 1, ipend
        IF (qit_conv_b(ip,k) > 0._wp .AND. t_b(ip,k) > ticeini) THEN
          wfac = MAX( 0._wp, MIN(1._wp,0.25_wp*(t_b(ip,k)-ticeini)) + &
                      0.25_wp*MIN(0._wp,t_b(ip,mtop_con_b(ip))-ticeini) )
          qct_conv_b(ip,k) = qct_conv_b(ip,k) + wfac        * qit_conv_b(ip,k)
          tt_conv_b (ip,k) = tt_conv_b (ip,k) - lfocpd*wfac * qit_conv_b(ip,k)
          qit_conv_b(ip,k) = (1._wp-wfac)                   * qit_conv_b(ip,k)
        ENDIF
      ENDDO
    ENDDO
    !$acc end parallel

    ! precipitation fluxes / rates  3d -> 2d
    !$acc parallel
    !$acc loop gang vector
    DO  ip = 1, ipend
      prr_con_b(ip) = zprr_con_b(ip,ke1)
      prs_con_b(ip) = zprs_con_b(ip,ke1)
    ENDDO
    !$acc end parallel

    ! other diagnostics:
    !$acc parallel
    !$acc loop gang vector
    DO ip = 1, ipend
      IF (locum_b(ip)) THEN
        top_con_b (ip) = REAL(mtop_con_b(ip),wp)
        bas_con_b (ip) = REAL(mbas_con_b(ip),wp)
        mflx_con_b(ip) = zcon_udd_b(ip,mbas_con_b(ip),1) ! updraft massflux at cloud base
      ELSE
        top_con_b (ip) = 0.0_wp
        bas_con_b (ip) = 0.0_wp
        mflx_con_b(ip) = 0.0_wp
      ENDIF
    ENDDO
    !$acc end parallel

    ! convective contribution to wind gust 
    ! (based on simple parameterization by Peter Bechthold)
    !$acc parallel
    !$acc loop gang vector
    DO ip = 1, ipend
      IF ( ktype_b(ip) == 1 )  THEN   ! penetrative convection
        zk850 = conv_k850_b(ip,iblock)
        zk950 = conv_k950_b(ip,iblock)
        ! We take the arithmetic mean of u_m_b(i,zk850) and u_m_b(i,zk850-1)
        ! as well as v_m_b(i,zk850) and v_m_b(i,zk850-1), since the levels 
        ! zk850 and zk950 are located just below the respective threshold
        ! heights.
        u850 = 0.5_wp * (u_m_b(ip,zk850) + u_m_b(ip,zk850-1))
        u950 = 0.5_wp * (u_m_b(ip,zk950) + u_m_b(ip,zk950-1))
        v850 = 0.5_wp * (v_m_b(ip,zk850) + v_m_b(ip,zk850-1))
        v950 = 0.5_wp * (v_m_b(ip,zk950) + v_m_b(ip,zk950-1))
        ! convective gust
        vgust_con_b(ip) = nwp_con_gust( u850, u950, v850, v950 )
      ELSE
        vgust_con_b(ip) = 0.0_wp
      ENDIF
    ENDDO
    !$acc end parallel

  ENDIF   ! itype_conv == 2

#ifndef MESSY
  ! loop over tracers (which undergo convective transport)
  !   US: this has to be added later, when we really have tracer support 
  !       in the copy-in/copy-out mechanism
#endif

  SELECT CASE (y_conv_closure)
    ! for Tiedtke and Tiedtke-Bechtold convection only the "standard" is implemented
    ! for shallow convection also the Boeing closure is possible

  CASE ("standard")
    ! Calculate the convective cloud cover by a simple empirical
    ! relation (following B.Ritter, FE14). Anvils are assumed for
    ! a temperature increase at top level

    !$acc parallel
    !$acc loop gang vector
    DO  ip = 1, ipend
      IF ( locum_b(ip) .AND. mtop_con_b(ip) > 0 ) THEN
!Tiedtke-Bechtold     IF ( locum_b(ip) .AND. mtop_con_b(ip) > 0 .AND. mtop_con(ip) > (ke-1) ) THEN
        mtop = mtop_con_b(ip)
        mbas = mbas_con_b(ip)
        zbas = 0.5_wp*( hhl_b(ip,mbas) + hhl_b(ip,mbas+1) )
        ztop = 0.5_wp*( hhl_b(ip,mtop) + hhl_b(ip,mtop+1) )
        !$acc loop seq
        DO  k = mtop, mbas-1
          clc_con_b(ip,k) = 0.35_wp*(ztop-zbas)/5000.0_wp
          IF ( k == mtop ) THEN
            zdtbot = t_b(ip,k+1) - t_b(ip,k  )
            zdttop = t_b(ip,k  ) - t_b(ip,k-1)
            IF ( zdtbot > 0.0_wp .AND. zdttop <= 0.0_wp ) THEN
              clc_con_b(ip,k) = 2.0_wp*clc_con_b(ip,k)
            ENDIF
          ENDIF
          clc_con_b(ip,k) = MIN ( 1.0_wp, MAX(0.05_wp, clc_con_b(ip,k)) )
        ENDDO
      ENDIF
    ENDDO
    !$acc end parallel

  CASE ("Boeing")
    !$acc parallel
    !$acc loop gang vector
    DO  ip = 1, ipend
      IF( locum_b(ip) .AND. mtop_con_b(ip) > 0 ) THEN
        mtop = mtop_con_b(ip)
        mbas = mbas_con_b(ip)
        zbas = 0.5_wp*( hhl_b(ip,mbas) + hhl_b(ip,mbas+1) )
        ztop = 0.5_wp*( hhl_b(ip,mtop) + hhl_b(ip,mtop+1) )
        !$acc loop seq
        DO  k = mtop, mbas-1
          IF (zmfu_b(ip,mbas) > 1.0e-6_wp) THEN
            ! SB14 Assume decrease in mass flux corresponds to decrease in cloud cover
            ! With little change in vertical velocity (as found by e.g. Neggers 2009)
            ! A Dual Mass Flux Framework for Boundary Layer Convection. Part II: Clouds
            ! Prefactor from Ouwersloot et al 2013
            ! Quantifying the transport of subcloud layer reactants by shallow
            ! cumulus clouds over the Amazon
            clc_con_b(ip,k) = (zmfu_b(ip,k)/zmfu_b(ip,mbas))*0.03_wp/0.84_wp
            zdtbot = t_b(ip,k+1) - t_b(ip,k  )
            zdttop = t_b(ip,k  ) - t_b(ip,k-1)
            IF ( zdtbot > 0.0_wp .AND. zdttop <= 0.0_wp .or. k==mbas-1) THEN
              ! SB14 Assume total cloud fraction (for radiation etc)
              ! is 3.3 the cloud core fraction
              ! Near an inversion (anvils), but also cloud base!
              ! 1.7 is a typical value for BOMEX (see Siebesma 2003)
              ! 2.45 is more typical for ARM (Ouwersloot et al 2013)
              ! Use a higher factor here to account for the fact that cloud cover
              ! Gets underestimated
              ! See also Ouwersloot et al 2013 for an example of the relation
              ! between cloud cover and cloud core/area
              clc_con_b(ip,k) = 3.33_wp*clc_con_b(ip,k)
              ! Correct (ad hoc) for the fact that the passive parts
              ! contain less liquid water, but do not overdo this
              clw_con_b(ip,k) = 0.6_wp*clw_con_b(ip,k)
            ELSE
              ! SB14 Assume total cloud fraction (for radiation etc)
              ! is twice the cloud core fraction
              clc_con_b(ip,k) = 2.0_wp*clc_con_b(ip,k)
              ! Correct (ad hoc) for the fact that the passive parts
              ! contain less liquid water, but do not overdo this
              clw_con_b(ip,k) = 0.7_wp*clw_con_b(ip,k)
            ENDIF
            ! SB14 Removed lower limiter on cloud fraction
            clc_con_b(ip,k) = MIN ( 1.0_wp, MAX(0.00_wp, clc_con_b(ip,k)) )
          ENDIF
        ENDDO
      ENDIF
    ENDDO
    !$acc end parallel

  END SELECT

  !$acc end data


  !GPU_TMP_UPDATE : this update should be remove once the physics run fully on GPU
#ifdef _OPENACC
  print*, "GPU-info : update host after conv_shallow_organize"
#endif
  !$acc update host(fif_b, fih_b, pf_b, ph_b, hhl_b, p0_b, pp_b, p0hl_b, pp_b)
  !$acc update host(ps_b, t_b, qv_b, dqvdt_conv_b, clc_con_b, clw_con_b)
  !$acc update host(bas_con_b, top_con_b, tt_conv_b, qvt_conv_b, tket_conv_b)
  !$acc update host(mflx_con_b, cape_con_b, qcvg_con_b)
  !$acc update host(t_g_b, shfl_s_b, lhfl_s_b)
! for Tiedtke
  !$acc update host(tke_b, qct_conv_b, qit_conv_b, ut_conv_b, vt_conv_b)
  !$acc update host(prr_con_b, prs_con_b, prne_con_b, tt_lheat_new_b)
  !$acc update host(tke_con_b, vgust_con_b)
! for Tiedtke-Bechtold
  !$acc update host(ttens_conv_b)

!------------------------------------------------------------------------------
! End of module procedure conv_shallow_organize
!------------------------------------------------------------------------------

END SUBROUTINE conv_organize

!==============================================================================
!==============================================================================
!+ Module procedure "conv_finalize" to finalize calls to shallow conv.
!------------------------------------------------------------------------------

SUBROUTINE conv_finalize

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine handles the finalization of the shallow convection module.
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

! Local parameters:
! ----------------

INTEGER :: izerror

!------------------------------------------------------------------------------
! Begin Subroutine conv_shallow_finalize
!------------------------------------------------------------------------------

  izerror = 0
  ! Deallocate the working arrays

#ifdef _OPENACC
  CALL conv_in_wkarr_dealloc
  CALL conv_wkarr_dealloc (itype_conv, izerror)

  IF (izerror /= 0) THEN
    CALL model_abort(my_cart_id, izerror, 'conv_wkarr_dealloc', 'conv_finalize')
  ENDIF
#endif

!------------------------------------------------------------------------------
! End of module procedure conv_finalize
!------------------------------------------------------------------------------

END SUBROUTINE conv_finalize

!==============================================================================

!-------------------------------------------------------------------------
!!
!! Calculate convective contribution to the wind gusts
!!     gust_conv = \alpha MAX(0,U_850 - U_950)
!! where \alpha=0.6 is a tunable constant and U_850-U_950 is the difference between
!! the 850 hPa and 950 hPa wind speeds, which represents the low-level wind shear.
!!
!! @par Literature
!! Bechthold, P. and J. Bidlot (2009): Parameterization of convective gusts. 
!! ECMWF Newsletter No. 119
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2014-03-25)
!!
ELEMENTAL FUNCTION nwp_con_gust( u_850, u_950, v_850, v_950) RESULT(vgust_con)

  REAL(wp), INTENT(IN) :: u_850, &    ! zonal wind component at 850 hPa [m/s]
    &                     u_950, &    ! zonal wind component at 950 hPa [m/s]
    &                     v_850, &    ! meridional wind component at 850 hPa [m/s]
    &                     v_950       ! meridional wind component at 950 hPa [m/s]

  REAL(wp) :: vgust_con               ! convective contribution to the wind gusts [m/s]

  REAL(wp), PARAMETER :: alpha = 0.6_wp ! convective mixing parameter

  vgust_con = alpha * MAX(0._wp, SQRT((u_850**2 + v_850**2)) - SQRT((u_950**2 + v_950**2)))

END FUNCTION nwp_con_gust

!==============================================================================

#ifdef _OPENACC
!==============================================================================
!+ Module procedure for allocating working arrays for the convection interface
!------------------------------------------------------------------------------

SUBROUTINE conv_in_wkarr_alloc (nproma, ke)

!------------------------------------------------------------------------------
!
! Purpose:   Allocation of the working arrays required in the convection
!            interface
!
!------------------------------------------------------------------------------

! Declarations:
INTEGER, INTENT(IN) :: &
  nproma, ke

INTEGER :: izl, ist

! End of header
!==============================================================================

  izl = 0
  ist = 0

  ALLOCATE ( zcumflx_b     (nproma)     ,STAT=izl ); zcumflx_b     = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( zcucape_b     (nproma)     ,STAT=izl ); zcucape_b     = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( zcuconv_b     (nproma)     ,STAT=izl ); zcuconv_b     = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( zco_clw_b     (nproma, ke) ,STAT=izl ); zco_clw_b     = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( zdt_con_b     (nproma, ke) ,STAT=izl ); zdt_con_b     = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( zdtdiab_con_b (nproma, ke) ,STAT=izl ); zdtdiab_con_b = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( zdqv_con_b    (nproma, ke) ,STAT=izl ); zdqv_con_b    = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( zdtke_con_b   (nproma, ke) ,STAT=izl ); zdtke_con_b   = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( mbas_con_b    (nproma)     ,STAT=izl ); mbas_con_b    = 0       ; ist=ist+izl
  ALLOCATE ( mtop_con_b    (nproma)     ,STAT=izl ); mtop_con_b    = 0       ; ist=ist+izl
  ALLOCATE ( locum_b       (nproma)     ,STAT=izl ); locum_b       = .FALSE. ; ist=ist+izl

  IF     (itype_conv == 0) THEN
    ! only for Tiedtke
    ALLOCATE ( zcutke_b    (nproma)     ,STAT=izl ); zcutke_b      = 0.0_wp  ; ist=ist+izl
    ALLOCATE ( zvddraf_b   (nproma)     ,STAT=izl ); zvddraf_b     = 0.0_wp  ; ist=ist+izl
    ALLOCATE ( zdtlhn_con_b(nproma, ke) ,STAT=izl ); zdtlhn_con_b  = 0.0_wp  ; ist=ist+izl
    ALLOCATE ( zdu_con_b   (nproma, ke) ,STAT=izl ); zdu_con_b     = 0.0_wp  ; ist=ist+izl
    ALLOCATE ( zdv_con_b   (nproma, ke) ,STAT=izl ); zdv_con_b     = 0.0_wp  ; ist=ist+izl
    ALLOCATE ( zqct_con_b  (nproma, ke) ,STAT=izl ); zqct_con_b    = 0.0_wp  ; ist=ist+izl
    ALLOCATE ( zqit_con_b  (nproma, ke) ,STAT=izl ); zqit_con_b    = 0.0_wp  ; ist=ist+izl
    ALLOCATE ( ztke_b      (nproma, ke) ,STAT=izl ); ztke_b        = 0.0_wp  ; ist=ist+izl
    ALLOCATE ( zpr_con_b   (nproma)     ,STAT=izl ); zpr_con_b     = 0.0_wp  ; ist=ist+izl
    ALLOCATE ( zps_con_b   (nproma)     ,STAT=izl ); zps_con_b     = 0.0_wp  ; ist=ist+izl
    ALLOCATE ( zpre_con_b  (nproma)     ,STAT=izl ); zpre_con_b    = 0.0_wp  ; ist=ist+izl
  ELSEIF (itype_conv == 3) THEN
    ! only for shallow
    ALLOCATE ( zmfu_b      (nproma, ke) ,STAT=izl ); zmfu_b        = 0.0_wp  ; ist=ist+izl
  ENDIF

  ! and for Tiedtke-Bechtold
  IF     (itype_conv == 2) THEN
    ALLOCATE ( ldshcv_b    (nproma)     ,STAT=izl ); ldshcv_b      = .FALSE. ; ist=ist+izl
    ALLOCATE ( zplitot_b   (nproma,ke)  ,STAT=izl ); zplitot_b     = 0.0_wp  ; ist=ist+izl
    ALLOCATE ( zomega_b    (nproma,ke)  ,STAT=izl ); zomega_b      = 0.0_wp  ; ist=ist+izl
    ALLOCATE ( zqhfl_b     (nproma,ke1) ,STAT=izl ); zqhfl_b       = 0.0_wp  ; ist=ist+izl
    ALLOCATE ( zshfl_b     (nproma,ke1) ,STAT=izl ); zshfl_b       = 0.0_wp  ; ist=ist+izl
    ALLOCATE ( zdpres_b    (nproma,ke)  ,STAT=izl ); zdpres_b      = 0.0_wp  ; ist=ist+izl
    ALLOCATE ( zdgeopot_b  (nproma,ke)  ,STAT=izl ); zdgeopot_b    = 0.0_wp  ; ist=ist+izl
    ALLOCATE ( ztt_conv_b  (nproma,ke)  ,STAT=izl ); ztt_conv_b    = 0.0_wp  ; ist=ist+izl
    ALLOCATE ( zqvt_conv_b (nproma,ke)  ,STAT=izl ); zqvt_conv_b   = 0.0_wp  ; ist=ist+izl
    ALLOCATE ( zprr_con_b  (nproma,ke1) ,STAT=izl ); zprr_con_b    = 0.0_wp  ; ist=ist+izl
    ALLOCATE ( zprs_con_b  (nproma,ke1) ,STAT=izl ); zprs_con_b    = 0.0_wp  ; ist=ist+izl
    ALLOCATE ( zcon_udd_b  (nproma,ke,6),STAT=izl ); zcon_udd_b    = 0.0_wp  ; ist=ist+izl
  ENDIF

  IF (ist /= 0) THEN
    CALL model_abort(my_cart_id, ist, 'Allocation of the interface working arrays failed', &
                     'conv_in_wkarr_alloc')
  ENDIF

  !$acc enter data                                                         &
  !$acc create (zcumflx_b, zcucape_b, zcuconv_b, zdt_con_b, zdtdiab_con_b) &
  !$acc create (zdqv_con_b, zdtke_con_b, zco_clw_b)                        &
  !$acc create (mbas_con_b, mtop_con_b, zmfu_b, locum_b)                   &
  !$acc create (zcutke_b, zvddraf_b, zdtlhn_con_b, zdu_con_b, zdv_con_b)   &
  !$acc create (zqct_con_b, zqit_con_b, zpr_con_b, zps_con_b, zpre_con_b)

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE conv_in_wkarr_alloc

!==============================================================================
!==============================================================================
!+ Module procedure for deallocating convection interface working arrays
!------------------------------------------------------------------------------

SUBROUTINE conv_in_wkarr_dealloc

!------------------------------------------------------------------------------
!
! Purpose:   Deallocation of the working arrays required in the convection
!            interface
!
!------------------------------------------------------------------------------

INTEGER :: izl, ist

! End of header
!==============================================================================

  izl = 0
  ist = 0

  !$acc exit data                                                          &
  !$acc delete (zcumflx_b, zcucape_b, zcuconv_b, zdt_con_b, zdtdiab_con_b) &
  !$acc delete (zdqv_con_b, zdtke_con_b, zco_clw_b)                        &
  !$acc delete (mbas_con_b, mtop_con_b, zmfu_b, locum_b)                   &
  !$acc delete (zcutke_b, zvddraf_b, zdtlhn_con_b, zdu_con_b, zdv_con_b)   &
  !$acc delete (zqct_con_b, zqit_con_b, zpr_con_b, zps_con_b, zpre_con_b)

  DEALLOCATE ( zcumflx_b      , STAT=izl ); ist=ist+izl
  DEALLOCATE ( zcucape_b      , STAT=izl ); ist=ist+izl
  DEALLOCATE ( zcuconv_b      , STAT=izl ); ist=ist+izl
  DEALLOCATE ( zco_clw_b      , STAT=izl ); ist=ist+izl
  DEALLOCATE ( zdt_con_b      , STAT=izl ); ist=ist+izl
  DEALLOCATE ( zdtdiab_con_b  , STAT=izl ); ist=ist+izl
  DEALLOCATE ( zdqv_con_b     , STAT=izl ); ist=ist+izl
  DEALLOCATE ( zdtke_con_b    , STAT=izl ); ist=ist+izl
  DEALLOCATE ( mbas_con_b     , STAT=izl ); ist=ist+izl
  DEALLOCATE ( mtop_con_b     , STAT=izl ); ist=ist+izl
  DEALLOCATE ( locum_b        , STAT=izl ); ist=ist+izl

  IF     (itype_conv == 0) THEN
    ! only for Tiedtke
    DEALLOCATE ( zcutke_b       , STAT=izl ); ist=ist+izl
    DEALLOCATE ( zvddraf_b      , STAT=izl ); ist=ist+izl
    DEALLOCATE ( zdtlhn_con_b   , STAT=izl ); ist=ist+izl
    DEALLOCATE ( zdu_con_b      , STAT=izl ); ist=ist+izl
    DEALLOCATE ( zdv_con_b      , STAT=izl ); ist=ist+izl
    DEALLOCATE ( zqct_con_b     , STAT=izl ); ist=ist+izl
    DEALLOCATE ( zqit_con_b     , STAT=izl ); ist=ist+izl
    DEALLOCATE ( ztke_b         , STAT=izl ); ist=ist+izl
    DEALLOCATE ( zpr_con_b      , STAT=izl ); ist=ist+izl
    DEALLOCATE ( zps_con_b      , STAT=izl ); ist=ist+izl
    DEALLOCATE ( zpre_con_b     , STAT=izl ); ist=ist+izl
  ELSEIF (itype_conv == 3) THEN
    ! only for shallow
    DEALLOCATE ( zmfu_b         , STAT=izl ); ist=ist+izl
  ENDIF

  ! and for Tiedtke-Bechtold
  IF     (itype_conv == 2) THEN
    DEALLOCATE ( ldshcv_b       , STAT=izl ); ist=ist+izl
    DEALLOCATE ( zplitot_b      , STAT=izl ); ist=ist+izl
    DEALLOCATE ( zomega_b       , STAT=izl ); ist=ist+izl
    DEALLOCATE ( zqhfl_b        , STAT=izl ); ist=ist+izl
    DEALLOCATE ( zshfl_b        , STAT=izl ); ist=ist+izl
    DEALLOCATE ( zdpres_b       , STAT=izl ); ist=ist+izl
    DEALLOCATE ( zdgeopot_b     , STAT=izl ); ist=ist+izl
    DEALLOCATE ( ztt_conv_b     , STAT=izl ); ist=ist+izl
    DEALLOCATE ( zqvt_conv_b    , STAT=izl ); ist=ist+izl
    DEALLOCATE ( zprr_con_b     , STAT=izl ); ist=ist+izl
    DEALLOCATE ( zprs_con_b     , STAT=izl ); ist=ist+izl
    DEALLOCATE ( zcon_udd_b     , STAT=izl ); ist=ist+izl
  ENDIF

  IF (ist /= 0) THEN
    CALL model_abort(my_cart_id, ist, 'Deallocation of the interface working arrays failed', &
                     'conv_in_wkarr_dealloc')
  ENDIF

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE conv_in_wkarr_dealloc

!==============================================================================
#endif
!==============================================================================

END MODULE conv_interface
