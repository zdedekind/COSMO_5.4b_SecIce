!+ Module src_block_fields_org for organizing block fields
!------------------------------------------------------------------------------

MODULE src_block_fields_org

!------------------------------------------------------------------------------
!
! Description:
!
! This module performs allocation/deallocation of block fields used in the
! Physics as well as set the association between fields and block fields
! using routine register_block_fields_all
!
!
! Current Code Owner: MeteoSwiss, Xavier Lapillonne
!  phone: +41 58 460 9237
!  fax  : +41 58 460 9278
!  email:  xavier.lapillonne@meteoswiss.ch
!
! History:
! Version      Date       Name
! ----------   ---------- ----
! V5_1         2014-11-28 Xavier Lapillonne
!  Initial release
! V5_2         2015-05-21 Ulrich Schaettler
!  Introduced some ifdef NUDGING directives to be able to compile without Nudging
!  Also allocate rho_b in SR block_fields_deallocate
! V5_3         2015-10-09 Ulrich Schaettler, Xavier Lapillonne
!  Adaptations to MCH version
! V5_3a        2015-11-24 Ulrich Schaettler
!  Only allocate variables necessary at the moment and extend the list when 
!  adding more parameterizations of the COSMO-ICON physics
!  Updated code owner information
! V5_4         2016-03-10 Xavier Lapillonne
!  Restructured module as a preparation for future optimization
!  (multi-copy)
! V5_4a        2016-05-10 Ulrich Schaettler
!  Extensions for turbulence scheme in blocked form
! V5_4b        2016-07-12 Ulrich Schaettler, Ulrich Blahak
!                         Xavier Lapillonne, Jochen Foerstner
!  Corrections for running with OPENACC (US)
!  Added h0noise for lartif_data (UB)
!  Added fields for surface schemes (US)
!  Added fields for convection schemes (XL,US,JF)
! @VERSION@    @DATE@     Ulrich Schaettler
!  Fixed some syntax errors in the GPU part
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Declarations:
!
! Modules used:
!
USE kind_parameters,  ONLY : wp
USE data_constants,   ONLY : rprecision

USE data_runcontrol,  ONLY : nproma,nlastproma,nblock,nnow,nnew,nold,l2tls, &
                             lmulti_snow,lmulti_layer,nhori, lgsp_first,    &
                             itype_vdif, lartif_data, llake, itype_conv

USE data_modelconfig, ONLY : ke,idt_qv,idt_qc,idt_qi,idt_qr,idt_qs,idt_qg,  &
                             ke_snow, ke_soil, kcm,ke1, lalloc_tke

USE data_parallel,    ONLY : my_cart_id

USE environment,      ONLY : model_abort

!------------------------------------------------------------------------------

USE data_fields,      ONLY : &
  ! constant fields from reference atmosphere 
  p0, p0hl, rho0, dp0, hhl,                                          &

  ! external parameters, needed in several parameterizations
  hsurf, fr_land, isoiltyp, plcov, rootdp, sai, tai, eai, l_pat,     &
  rsmin2d, sso_stdh, depth_lk, fc, llandmask, llakemask,             &

  ! dynamical variables
  u_m, v_m, w, t, pp, tke, rho, utens, vtens, ttens, tketens,        &

  ! fields for surface values
  ps, t_s, t_g, qv_s,                                                &

  ! from turbulence scheme
  l_hori, gz0, tcm, tch, tfm, tfh, tfv, tvm, tvh, tkr, tkred_sfc,    &
  tkvm, tkvh, rcld, tkhm, tkhh, hdef2, hdiv, dwdx, dwdy, edr,        &
  tket_sso, tket_hshr, tket_adv, dqvdt, ut_sso, vt_sso,              &
  t_2m, qv_2m, td_2m, rh_2m, u_10m, v_10m, ut_turb, vt_turb,         &
  shfl_s, lhfl_s, qvsflx, h0noise,                                   &

  ! from microphysics
  prr_gsp, prs_gsp, prg_gsp, tinc_lh, qrs,                           &

  ! from the surface schemes
  ! already above: t_s, t_g, ps, qv_s, u_10m, v_10m, 
  !                prr_gsp, prs_gsp, prg_gsp, tch, tcm, tfv, 
  !                shfl_s, lhfl_s
  t_snow, w_snow, h_snow, rho_snow, freshsnow, fr_snow, w_i, w_p,    &
  w_s, t_so, w_so, w_so_ice, t_snow_mult, rho_snow_mult, wliq_snow,  &
  w_snow_mult, dzh_snow_mult, runoff_s, runoff_g, lhfl_bs, lhfl_pl,  &
  rstom, h_ice, t_ice,                                               &
  ! not yet: meltrate, conv_frac, h_snow_gp:     only in interface
  !          shfl_snow, lhfl_snow:               only in interface
  !          shfl_sfc, lhfl_sfc, qhfl_sfc: Optional
  ! from FLake
  qmomflux, fr_lake, depth_lk, fetch_lk, dp_bs_lk, t_bs_lk, gamso_lk,&
  t_mnw_lk, t_wml_lk, t_bot_lk, t_b1_lk, c_t_lk, h_ml_lk, h_b1_lk,   &

  ! from convection schemes
  dqvdt_conv, clc_con, clw_con, prr_con, prs_con, prne_con,          &
  bas_con, top_con, tt_conv, ttdiab_conv, qvt_conv, qct_conv,        &
  qit_conv, ut_conv, vt_conv, tket_conv, mflx_con, cape_con,         &
  qcvg_con, tke_con, vgust_con, qhfl, w_conv, ttens_conv,            &

  ! special treatment for radiation variables
  sobs, thbs, pabs

USE src_stoch_physics, ONLY : pertstoph

#ifdef NUDGING
USE data_lheat_nudge,  ONLY : tt_lheat, qrsflux
#endif

USE data_block_fields       ! all fields are used

!------------------------------------------------------------------------------

USE src_tracer, ONLY       : trcr_get_block,trcr_errorstr,trcr_get

USE src_block_fields, ONLY : register_block_field, finalize_block_fields, &
                             !for testing:
                             register_copy, CopylistStruct,               &
                             copyToBlockF, copyFromBlockF,                &
                             mind_ilon, mind_jlat, request_copy,          &
                             copy_to_block, copy_from_block,              &
                             finalize_copy, init_copy_list

!------------------------------------------------------------------------------

IMPLICIT NONE

!------------------------------------------------------------------------------

PRIVATE

!==============================================================================

!------------------------------------------------------------------------------
! Global (i.e. public) Declarations:

PUBLIC :: block_fields_allocate, block_fields_deallocate,                 &
          block_fields_register_all, block_fields_cleanup,                &
          block_fields_copytoblock_tke, block_fields_copyfromblock_tke

! public methods for testing module correctness
PUBLIC :: block_fields_test_copy

!==============================================================================
! ACC variables : list of variables which should be allocated and deallocated
!                 on the accelerator. A macro is used in order to avoid
!                 duplication in multiple subroutines
!------------------------------------------------------------------------------

!Note:  add new variables if necessary (try to keep a grouping of variables)
!       if a new ACC_LIST is introduced it should be added in all routines
!       in this file which uses ACC_LIST

! constant fields and external parameters
#define ACC_LIST01    p0_b       , p0hl_b     , rho0_b     , dp0_b      , dz_b       , hhl_b
#define ACC_LIST02    hsurf_b    , fr_land_b  , isoiltyp_b , plcov_b    , rootdp_b
#define ACC_LIST03    sai_b      , tai_b      , eai_b      , l_pat_b    , rsmin2d_b  , depth_lk_b
#define ACC_LIST04    fr_lake_b  , sso_stdh_b , fc_b       , llandmask_b, llakemask_b

! dynamical variables
#define ACC_LIST05    ptot_b        , pp_b          , t_b           , u_m_b          , v_m_b
#define ACC_LIST06    w_b           , rho_b         , tke_b         , t_new_b        , pp_new_b
#define ACC_LIST07    qv_b          , qc_b          , qi_b          , qr_b           , qs_b          , qg_b
#define ACC_LIST08    qv_new_b      , qc_new_b      , qi_new_b      , qr_new_b       , qs_new_b      , qg_new_b
#define ACC_LIST09    utens_b       , vtens_b       , ttens_b       , qvtens_b       , qctens_b
#define ACC_LIST10    qitens_b      , tketens_b

! surface values
#define ACC_LIST11    ps_b          , t_s_b         , t_s_new_b     , t_g_b          , t_g_new_b
#define ACC_LIST12    qv_s_b        , qv_s_new_b

! from microphysics
#define ACC_LIST13    prr_gsp_b     , prs_gsp_b     , prg_gsp_b     , tinc_lh_b      , qrs_b
#define ACC_LIST14    tt_lheat_b    , tt_lheat_new_b, qrsflux_b

! from turbulence scheme
#define ACC_LIST15    l_hori_b      , gz0_b         , tcm_b         , tch_b          , tfm_b
#define ACC_LIST16    tfh_b         , tfv_b         , tvm_b         , tvh_b          , tkr_b
#define ACC_LIST17    tkred_sfc_b   , tkvm_b        , tkvh_b        , rcld_b         , tkhm_b
#define ACC_LIST18    tkhh_b        , hdef2_b       , hdiv_b        , dwdx_b         , dwdy_b
#define ACC_LIST19    edr_b         , tket_sso_b    , tket_hshr_b   , tket_adv_b
#define ACC_LIST20    dqvdt_b       , ut_turb_b     , vt_turb_b     , ut_sso_b       , vt_sso_b
#define ACC_LIST21    t_2m_b        , qv_2m_b       , td_2m_b       , rh_2m_b        , u_10m_b
#define ACC_LIST22    v_10m_b       , shfl_s_b      , lhfl_s_b      , qvsflx_b

! from surface schemes
#define ACC_LIST23    t_snow_b      , t_snow_new_b  , w_snow_b      , w_snow_new_b   , h_snow_b
#define ACC_LIST24    h_snow_new_b  , rho_snow_b    , rho_snow_new_b, freshsnow_b    , fr_snow_b
#define ACC_LIST25    w_i_b         , w_i_new_b     , w_p_b         , w_p_new_b      , w_s_b
#define ACC_LIST26    w_s_new_b     , t_so_b        , t_so_new_b    , w_so_b         , w_so_new_b
#define ACC_LIST27    w_so_ice_b    , w_so_ice_new_b
#define ACC_LIST28    runoff_s_b    , runoff_g_b    , lhfl_bs_b     , lhfl_pl_b      , rstom_b
#define ACC_LIST29    t_snow_mult_b , w_snow_mult_b , wliq_snow_b   , rho_snow_mult_b,dzh_snow_mult_b
#define ACC_LIST30    h_ice_b       , h_ice_new_b   , t_ice_b       , t_ice_new_b
#define ACC_LIST31      t_snow_mult_new_b   ,   w_snow_mult_new_b   , wliq_snow_new_b
#define ACC_LIST32    rho_snow_mult_new_b   , dzh_snow_mult_new_b

! from flake or seaice scheme
#define ACC_LIST33    fetch_lk_b    , dp_bs_lk_b    , t_bs_lk_b     , gamso_lk_b     , t_mnw_lk_b
#define ACC_LIST34    t_mnw_lk_new_b, t_wml_lk_b    , t_wml_lk_new_b, t_bot_lk_b     , t_bot_lk_new_b
#define ACC_LIST35    t_b1_lk_b     , t_b1_lk_new_b , c_t_lk_b      , c_t_lk_new_b   , h_ml_lk_b
#define ACC_LIST36    h_ml_lk_new_b , h_b1_lk_b     , h_b1_lk_new_b
#define ACC_LIST37    qmomflux_b

! from convection schemes
#define ACC_LIST38    fif_b         , fih_b         , pf_b          , ph_b
#define ACC_LIST39    w_conv_b      , dqvdt_conv_b  , qhfl_b        , ttens_conv_b
#define ACC_LIST40    clc_con_b     , clw_con_b     , prr_con_b     , prs_con_b      , prne_con_b
#define ACC_LIST41    bas_con_b     , top_con_b     , tt_conv_b     , ttdiab_conv_b  , qvt_conv_b
#define ACC_LIST42    qct_conv_b    , qit_conv_b    , ut_conv_b     , vt_conv_b      , tket_conv_b
#define ACC_LIST43    mflx_con_b    , cape_con_b    , qcvg_con_b    , tke_con_b      , vgust_con_b

! miscellaneous fields
#define ACC_LIST44    pertstoph_b   , h0noise_b

!       if a new ACC_LIST is introduced it should be added in all routines 
!       in this file which uses ACC_LIST

!==============================================================================
! Module procedures in src_block_fields_org
!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure "block_fields_allocate" in "src_block_fields_org"
!------------------------------------------------------------------------------

SUBROUTINE block_fields_allocate( ist )

!------------------------------------------------------------------------------
!
! Description:
!   Allocate block fields used in the physics
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------
  INTEGER, INTENT(OUT) :: ist

! Locals:
! -------
  INTEGER :: izl
  REAL(KIND=wp), PARAMETER :: &
!   init_val=TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_wp)
                              ! Produce a NaN  with  most compiler    
    r_init_val=-99999999_wp   ! It is an error to use a block field
                              ! without assigning it a value before.
                              ! All fields are therefore initialized to
                              ! non realistic value
  INTEGER      , PARAMETER :: &
    i_init_val=-99999999      ! non-realistic integer values

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine block_fields_allocate
!------------------------------------------------------------------------------

  ist=0
  izl=0

  !PHYS_EDIT: add new variables if necessary

  ! constant fields for the reference atmosphere
  ALLOCATE(p0_b           (nproma,ke)  , STAT=izl); p0_b           = r_init_val; ist=ist+izl
  ALLOCATE(p0hl_b         (nproma,ke+1), STAT=izl); p0hl_b         = r_init_val; ist=ist+izl
  ALLOCATE(rho0_b         (nproma,ke)  , STAT=izl); rho0_b         = r_init_val; ist=ist+izl
  ALLOCATE(dp0_b          (nproma,ke)  , STAT=izl); dp0_b          = r_init_val; ist=ist+izl
  ALLOCATE(hhl_b          (nproma,ke+1), STAT=izl); hhl_b          = r_init_val; ist=ist+izl
  ALLOCATE(dz_b           (nproma,ke)  , STAT=izl); dz_b           = r_init_val; ist=ist+izl

  ! external parameters, needed in several parameterizations
  ALLOCATE(hsurf_b        (nproma)     , STAT=izl); hsurf_b        = r_init_val; ist=ist+izl
  ALLOCATE(fr_land_b      (nproma)     , STAT=izl); fr_land_b      = r_init_val; ist=ist+izl
  ALLOCATE(llandmask_b    (nproma)     , STAT=izl); llandmask_b    = .FALSE.   ; ist=ist+izl
  ALLOCATE(llakemask_b    (nproma)     , STAT=izl); llakemask_b    = .FALSE.   ; ist=ist+izl
  ALLOCATE(isoiltyp_b     (nproma)     , STAT=izl); isoiltyp_b     = i_init_val; ist=ist+izl
  ALLOCATE(plcov_b        (nproma)     , STAT=izl); plcov_b        = r_init_val; ist=ist+izl
  ALLOCATE(rootdp_b       (nproma)     , STAT=izl); rootdp_b       = r_init_val; ist=ist+izl
  IF (itype_vdif > -2) THEN
  ALLOCATE(sso_stdh_b     (nproma)     , STAT=izl); sso_stdh_b     = r_init_val; ist=ist+izl
  ENDIF
  ALLOCATE(sai_b          (nproma)     , STAT=izl); sai_b          = r_init_val; ist=ist+izl
  ALLOCATE(tai_b          (nproma)     , STAT=izl); tai_b          = r_init_val; ist=ist+izl
  ALLOCATE(eai_b          (nproma)     , STAT=izl); eai_b          = r_init_val; ist=ist+izl
  ALLOCATE(l_pat_b        (nproma)     , STAT=izl); l_pat_b        = r_init_val; ist=ist+izl
  ALLOCATE(rsmin2d_b      (nproma)     , STAT=izl); rsmin2d_b      = r_init_val; ist=ist+izl
  ALLOCATE(depth_lk_b     (nproma)     , STAT=izl); depth_lk_b     = r_init_val; ist=ist+izl
  ALLOCATE(fr_lake_b      (nproma)     , STAT=izl); fr_lake_b      = r_init_val; ist=ist+izl
  ALLOCATE(fc_b           (nproma)     , STAT=izl); fc_b           = r_init_val; ist=ist+izl

  ! dynamical variables
  ALLOCATE(ptot_b         (nproma,ke)  , STAT=izl); ptot_b         = r_init_val; ist=ist+izl
  ALLOCATE(pp_b           (nproma,ke)  , STAT=izl); pp_b           = r_init_val; ist=ist+izl
  ALLOCATE(t_b            (nproma,ke)  , STAT=izl); t_b            = r_init_val; ist=ist+izl
  ALLOCATE(u_m_b          (nproma,ke)  , STAT=izl); u_m_b          = r_init_val; ist=ist+izl
  ALLOCATE(v_m_b          (nproma,ke)  , STAT=izl); v_m_b          = r_init_val; ist=ist+izl
  ALLOCATE(w_b            (nproma,ke+1), STAT=izl); w_b            = r_init_val; ist=ist+izl
  ALLOCATE(rho_b          (nproma,ke)  , STAT=izl); rho_b          = r_init_val; ist=ist+izl
  ALLOCATE(qv_b           (nproma,ke)  , STAT=izl); qv_b           = r_init_val; ist=ist+izl
  ALLOCATE(qc_b           (nproma,ke)  , STAT=izl); qc_b           = r_init_val; ist=ist+izl
  ALLOCATE(qi_b           (nproma,ke)  , STAT=izl); qi_b           = r_init_val; ist=ist+izl
  ALLOCATE(qr_b           (nproma,ke)  , STAT=izl); qr_b           = r_init_val; ist=ist+izl
  ALLOCATE(qs_b           (nproma,ke)  , STAT=izl); qs_b           = r_init_val; ist=ist+izl
  ALLOCATE(qg_b           (nproma,ke)  , STAT=izl); qg_b           = r_init_val; ist=ist+izl
  IF (ALLOCATED(tke)) THEN
    ALLOCATE(tke_b(nproma,ke+1,SIZE(tke,4))); tke_b = r_init_val ; ist=ist+izl
  ELSE !this is required for the exit data if lphys=.TRUE. and ltur=.FALSE.
    ALLOCATE(tke_b(nproma,ke+1,1)); tke_b = r_init_val ; ist=ist+izl
  END IF
  ALLOCATE(t_new_b        (nproma,ke)  , STAT=izl); t_new_b        = r_init_val; ist=ist+izl
  ALLOCATE(pp_new_b       (nproma,ke)  , STAT=izl); pp_new_b       = r_init_val; ist=ist+izl
  ALLOCATE(qv_new_b       (nproma,ke)  , STAT=izl); qv_new_b       = r_init_val; ist=ist+izl
  ALLOCATE(qc_new_b       (nproma,ke)  , STAT=izl); qc_new_b       = r_init_val; ist=ist+izl
  ALLOCATE(qi_new_b       (nproma,ke)  , STAT=izl); qi_new_b       = r_init_val; ist=ist+izl
  ALLOCATE(qr_new_b       (nproma,ke)  , STAT=izl); qr_new_b       = r_init_val; ist=ist+izl
  ALLOCATE(qs_new_b       (nproma,ke)  , STAT=izl); qs_new_b       = r_init_val; ist=ist+izl
  ALLOCATE(qg_new_b       (nproma,ke)  , STAT=izl); qg_new_b       = r_init_val; ist=ist+izl

  ALLOCATE(utens_b        (nproma,ke)  , STAT=izl); utens_b        = r_init_val; ist=ist+izl
  ALLOCATE(vtens_b        (nproma,ke)  , STAT=izl); vtens_b        = r_init_val; ist=ist+izl
  ALLOCATE(ttens_b        (nproma,ke)  , STAT=izl); ttens_b        = r_init_val; ist=ist+izl
  ALLOCATE(tketens_b      (nproma,ke+1), STAT=izl); tketens_b      = r_init_val; ist=ist+izl
  ALLOCATE(qvtens_b       (nproma,ke)  , STAT=izl); qvtens_b       = r_init_val; ist=ist+izl
  ALLOCATE(qctens_b       (nproma,ke)  , STAT=izl); qctens_b       = r_init_val; ist=ist+izl
  ALLOCATE(qitens_b       (nproma,ke)  , STAT=izl); qitens_b       = r_init_val; ist=ist+izl

  ! fields for surface values
  ALLOCATE(ps_b           (nproma)     , STAT=izl); ps_b           = r_init_val; ist=ist+izl
  ALLOCATE(t_s_b          (nproma)     , STAT=izl); t_s_b          = r_init_val; ist=ist+izl
  ALLOCATE(t_s_new_b      (nproma)     , STAT=izl); t_s_new_b      = r_init_val; ist=ist+izl
  ALLOCATE(t_g_b          (nproma)     , STAT=izl); t_g_b          = r_init_val; ist=ist+izl
  ALLOCATE(t_g_new_b      (nproma)     , STAT=izl); t_g_new_b      = r_init_val; ist=ist+izl
  ALLOCATE(qv_s_b         (nproma)     , STAT=izl); qv_s_b         = r_init_val; ist=ist+izl
  ALLOCATE(qv_s_new_b     (nproma)     , STAT=izl); qv_s_new_b     = r_init_val; ist=ist+izl

  ! from microphysics
  ALLOCATE(tinc_lh_b      (nproma,ke)  , STAT=izl); tinc_lh_b      = r_init_val; ist=ist+izl
  ALLOCATE(qrs_b          (nproma,ke)  , STAT=izl); qrs_b          = r_init_val; ist=ist+izl
  ALLOCATE(prr_gsp_b      (nproma)     , STAT=izl); prr_gsp_b      = r_init_val; ist=ist+izl
  ALLOCATE(prs_gsp_b      (nproma)     , STAT=izl); prs_gsp_b      = r_init_val; ist=ist+izl
  ALLOCATE(prg_gsp_b      (nproma)     , STAT=izl); prg_gsp_b      = r_init_val; ist=ist+izl
#ifdef NUDGING
  ALLOCATE(tt_lheat_b     (nproma,ke)  , STAT=izl); tt_lheat_b     = r_init_val; ist=ist+izl
  ALLOCATE(tt_lheat_new_b (nproma,ke)  , STAT=izl); tt_lheat_new_b = r_init_val; ist=ist+izl
  ALLOCATE(qrsflux_b      (nproma,ke)  , STAT=izl); qrsflux_b      = r_init_val; ist=ist+izl
#endif

! from radiation scheme (needed in surface TERRA)
  ALLOCATE(sobs_b         (nproma)     , STAT=izl); sobs_b         = r_init_val; ist=ist+izl
  ALLOCATE(thbs_b         (nproma)     , STAT=izl); thbs_b         = r_init_val; ist=ist+izl
  ALLOCATE(pabs_b         (nproma)     , STAT=izl); pabs_b         = r_init_val; ist=ist+izl

  ! from turbulence scheme
  ALLOCATE(l_hori_b       (nproma)     , STAT=izl); l_hori_b       = r_init_val; ist=ist+izl
  ALLOCATE(gz0_b          (nproma)     , STAT=izl); gz0_b          = r_init_val; ist=ist+izl
  ALLOCATE(tcm_b          (nproma)     , STAT=izl); tcm_b          = r_init_val; ist=ist+izl
  ALLOCATE(tch_b          (nproma)     , STAT=izl); tch_b          = r_init_val; ist=ist+izl
  ALLOCATE(tfm_b          (nproma)     , STAT=izl); tfm_b          = r_init_val; ist=ist+izl
  ALLOCATE(tfh_b          (nproma)     , STAT=izl); tfh_b          = r_init_val; ist=ist+izl
  ALLOCATE(tfv_b          (nproma)     , STAT=izl); tfv_b          = r_init_val; ist=ist+izl
  ALLOCATE(tvm_b          (nproma)     , STAT=izl); tvm_b          = r_init_val; ist=ist+izl
  ALLOCATE(tvh_b          (nproma)     , STAT=izl); tvh_b          = r_init_val; ist=ist+izl
  ALLOCATE(tkr_b          (nproma)     , STAT=izl); tkr_b          = r_init_val; ist=ist+izl
  ALLOCATE(tkred_sfc_b    (nproma)     , STAT=izl); tkred_sfc_b    = r_init_val; ist=ist+izl
  ALLOCATE(rcld_b         (nproma,ke+1), STAT=izl); rcld_b         = r_init_val; ist=ist+izl
  ALLOCATE(tkvm_b         (nproma,ke+1), STAT=izl); tkvm_b         = r_init_val; ist=ist+izl
  ALLOCATE(tkvh_b         (nproma,ke+1), STAT=izl); tkvh_b         = r_init_val; ist=ist+izl
  ALLOCATE(tkhm_b         (nproma,ke+1), STAT=izl); tkhm_b         = r_init_val; ist=ist+izl
  ALLOCATE(tkhh_b         (nproma,ke+1), STAT=izl); tkhh_b         = r_init_val; ist=ist+izl
  ALLOCATE(hdef2_b        (nproma,ke+1), STAT=izl); hdef2_b        = r_init_val; ist=ist+izl
  ALLOCATE(hdiv_b         (nproma,ke+1), STAT=izl); hdiv_b         = r_init_val; ist=ist+izl
  ALLOCATE(dwdx_b         (nproma,ke+1), STAT=izl); dwdx_b         = r_init_val; ist=ist+izl
  ALLOCATE(dwdy_b         (nproma,ke+1), STAT=izl); dwdy_b         = r_init_val; ist=ist+izl
  ALLOCATE(edr_b          (nproma,ke+1), STAT=izl); edr_b          = r_init_val; ist=ist+izl
  ALLOCATE(tket_sso_b     (nproma,ke+1), STAT=izl); tket_sso_b     = r_init_val; ist=ist+izl
  ALLOCATE(tket_hshr_b    (nproma,ke+1), STAT=izl); tket_hshr_b    = r_init_val; ist=ist+izl
  ALLOCATE(tket_adv_b     (nproma,ke+1), STAT=izl); tket_adv_b     = r_init_val; ist=ist+izl
  ALLOCATE(dqvdt_b        (nproma,ke)  , STAT=izl); dqvdt_b        = r_init_val; ist=ist+izl
  ALLOCATE(ut_turb_b      (nproma,ke)  , STAT=izl); ut_turb_b      = r_init_val; ist=ist+izl
  ALLOCATE(vt_turb_b      (nproma,ke)  , STAT=izl); vt_turb_b      = r_init_val; ist=ist+izl
  ALLOCATE(ut_sso_b       (nproma,ke)  , STAT=izl); ut_sso_b       = r_init_val; ist=ist+izl
  ALLOCATE(vt_sso_b       (nproma,ke)  , STAT=izl); vt_sso_b       = r_init_val; ist=ist+izl

  ALLOCATE(t_2m_b         (nproma)     , STAT=izl); t_2m_b         = r_init_val; ist=ist+izl
  ALLOCATE(qv_2m_b        (nproma)     , STAT=izl); qv_2m_b        = r_init_val; ist=ist+izl
  ALLOCATE(td_2m_b        (nproma)     , STAT=izl); td_2m_b        = r_init_val; ist=ist+izl
  ALLOCATE(rh_2m_b        (nproma)     , STAT=izl); rh_2m_b        = r_init_val; ist=ist+izl
  ALLOCATE(u_10m_b        (nproma)     , STAT=izl); u_10m_b        = r_init_val; ist=ist+izl
  ALLOCATE(v_10m_b        (nproma)     , STAT=izl); v_10m_b        = r_init_val; ist=ist+izl
  ALLOCATE(shfl_s_b       (nproma)     , STAT=izl); shfl_s_b       = r_init_val; ist=ist+izl
  ALLOCATE(lhfl_s_b       (nproma)     , STAT=izl); lhfl_s_b       = r_init_val; ist=ist+izl
  ALLOCATE(qvsflx_b       (nproma)     , STAT=izl); qvsflx_b       = r_init_val; ist=ist+izl

  ! from the surface schemes
  ALLOCATE(t_snow_b       (nproma)     , STAT=izl); t_snow_b       = r_init_val; ist=ist+izl
  ALLOCATE(t_snow_new_b   (nproma)     , STAT=izl); t_snow_new_b   = r_init_val; ist=ist+izl
  ALLOCATE(w_snow_b       (nproma)     , STAT=izl); w_snow_b       = r_init_val; ist=ist+izl
  ALLOCATE(w_snow_new_b   (nproma)     , STAT=izl); w_snow_new_b   = r_init_val; ist=ist+izl
  ALLOCATE(h_snow_b       (nproma)     , STAT=izl); h_snow_b       = r_init_val; ist=ist+izl
  ALLOCATE(h_snow_new_b   (nproma)     , STAT=izl); h_snow_new_b   = r_init_val; ist=ist+izl
  ALLOCATE(rho_snow_b     (nproma)     , STAT=izl); rho_snow_b     = r_init_val; ist=ist+izl
  ALLOCATE(rho_snow_new_b (nproma)     , STAT=izl); rho_snow_new_b = r_init_val; ist=ist+izl
  ALLOCATE(freshsnow_b    (nproma)     , STAT=izl); freshsnow_b    = r_init_val; ist=ist+izl
  ALLOCATE(fr_snow_b      (nproma)     , STAT=izl); fr_snow_b      = r_init_val; ist=ist+izl
  ALLOCATE(w_i_b          (nproma)     , STAT=izl); w_i_b          = r_init_val; ist=ist+izl
  ALLOCATE(w_i_new_b      (nproma)     , STAT=izl); w_i_new_b      = r_init_val; ist=ist+izl
  ALLOCATE(w_p_b          (nproma)     , STAT=izl); w_p_b          = r_init_val; ist=ist+izl
  ALLOCATE(w_p_new_b      (nproma)     , STAT=izl); w_p_new_b      = r_init_val; ist=ist+izl
  ALLOCATE(w_s_b          (nproma)     , STAT=izl); w_s_b          = r_init_val; ist=ist+izl
  ALLOCATE(w_s_new_b      (nproma)     , STAT=izl); w_s_new_b      = r_init_val; ist=ist+izl

  ALLOCATE(t_so_b         (nproma,0:ke_soil+1), STAT=izl); t_so_b         = r_init_val; ist=ist+izl
  ALLOCATE(t_so_new_b     (nproma,0:ke_soil+1), STAT=izl); t_so_new_b     = r_init_val; ist=ist+izl
  ALLOCATE(w_so_b         (nproma,  ke_soil+1), STAT=izl); w_so_b         = r_init_val; ist=ist+izl
  ALLOCATE(w_so_new_b     (nproma,  ke_soil+1), STAT=izl); w_so_new_b     = r_init_val; ist=ist+izl
  ALLOCATE(w_so_ice_b     (nproma,  ke_soil+1), STAT=izl); w_so_ice_b     = r_init_val; ist=ist+izl
  ALLOCATE(w_so_ice_new_b (nproma,  ke_soil+1), STAT=izl); w_so_ice_new_b = r_init_val; ist=ist+izl

  ALLOCATE(t_snow_mult_b      (nproma,0:ke_snow), STAT=izl); t_snow_mult_b      = r_init_val; ist=ist+izl
  ALLOCATE(t_snow_mult_new_b  (nproma,0:ke_snow), STAT=izl); t_snow_mult_new_b  = r_init_val; ist=ist+izl
  ALLOCATE(w_snow_mult_b      (nproma,  ke_snow), STAT=izl); w_snow_mult_b      = r_init_val; ist=ist+izl
  ALLOCATE(w_snow_mult_new_b  (nproma,  ke_snow), STAT=izl); w_snow_mult_new_b  = r_init_val; ist=ist+izl
  ALLOCATE(wliq_snow_b        (nproma,  ke_snow), STAT=izl); wliq_snow_b        = r_init_val; ist=ist+izl
  ALLOCATE(wliq_snow_new_b    (nproma,  ke_snow), STAT=izl); wliq_snow_new_b    = r_init_val; ist=ist+izl
  ALLOCATE(rho_snow_mult_b    (nproma,  ke_snow), STAT=izl); rho_snow_mult_b    = r_init_val; ist=ist+izl
  ALLOCATE(rho_snow_mult_new_b(nproma,  ke_snow), STAT=izl); rho_snow_mult_new_b= r_init_val; ist=ist+izl
  ALLOCATE(dzh_snow_mult_b    (nproma,  ke_snow), STAT=izl); dzh_snow_mult_b    = r_init_val; ist=ist+izl
  ALLOCATE(dzh_snow_mult_new_b(nproma,  ke_snow), STAT=izl); dzh_snow_mult_new_b= r_init_val; ist=ist+izl

  ALLOCATE(runoff_s_b     (nproma)     , STAT=izl); runoff_s_b     = r_init_val; ist=ist+izl
  ALLOCATE(runoff_g_b     (nproma)     , STAT=izl); runoff_g_b     = r_init_val; ist=ist+izl
  ALLOCATE(lhfl_bs_b      (nproma)     , STAT=izl); lhfl_bs_b      = r_init_val; ist=ist+izl
  ALLOCATE(lhfl_pl_b   (nproma,ke_soil), STAT=izl); lhfl_pl_b      = r_init_val; ist=ist+izl
  ALLOCATE(rstom_b        (nproma)     , STAT=izl); rstom_b        = r_init_val; ist=ist+izl

  ! from flake or seaice scheme
  IF (llake) THEN
  ALLOCATE(qmomflux_b     (nproma)     , STAT=izl); qmomflux_b     = r_init_val; ist=ist+izl
  ALLOCATE(fetch_lk_b     (nproma)     , STAT=izl); fetch_lk_b     = r_init_val; ist=ist+izl
  ALLOCATE(dp_bs_lk_b     (nproma)     , STAT=izl); dp_bs_lk_b     = r_init_val; ist=ist+izl
  ALLOCATE(t_bs_lk_b      (nproma)     , STAT=izl); t_bs_lk_b      = r_init_val; ist=ist+izl
  ALLOCATE(gamso_lk_b     (nproma)     , STAT=izl); gamso_lk_b     = r_init_val; ist=ist+izl
  ALLOCATE(t_mnw_lk_b     (nproma)     , STAT=izl); t_mnw_lk_b     = r_init_val; ist=ist+izl
  ALLOCATE(t_mnw_lk_new_b (nproma)     , STAT=izl); t_mnw_lk_new_b = r_init_val; ist=ist+izl
  ALLOCATE(t_wml_lk_b     (nproma)     , STAT=izl); t_wml_lk_b     = r_init_val; ist=ist+izl
  ALLOCATE(t_wml_lk_new_b (nproma)     , STAT=izl); t_wml_lk_new_b = r_init_val; ist=ist+izl
  ALLOCATE(t_bot_lk_b     (nproma)     , STAT=izl); t_bot_lk_b     = r_init_val; ist=ist+izl
  ALLOCATE(t_bot_lk_new_b (nproma)     , STAT=izl); t_bot_lk_new_b = r_init_val; ist=ist+izl
  ALLOCATE(t_b1_lk_b      (nproma)     , STAT=izl); t_b1_lk_b      = r_init_val; ist=ist+izl
  ALLOCATE(t_b1_lk_new_b  (nproma)     , STAT=izl); t_b1_lk_new_b  = r_init_val; ist=ist+izl
  ALLOCATE(c_t_lk_b       (nproma)     , STAT=izl); c_t_lk_b       = r_init_val; ist=ist+izl
  ALLOCATE(c_t_lk_new_b   (nproma)     , STAT=izl); c_t_lk_new_b   = r_init_val; ist=ist+izl
  ALLOCATE(h_ml_lk_b      (nproma)     , STAT=izl); h_ml_lk_b      = r_init_val; ist=ist+izl
  ALLOCATE(h_ml_lk_new_b  (nproma)     , STAT=izl); h_ml_lk_new_b  = r_init_val; ist=ist+izl
  ALLOCATE(h_b1_lk_b      (nproma)     , STAT=izl); h_b1_lk_b      = r_init_val; ist=ist+izl
  ALLOCATE(h_b1_lk_new_b  (nproma)     , STAT=izl); h_b1_lk_new_b  = r_init_val; ist=ist+izl
  ENDIF

  ALLOCATE(h_ice_b        (nproma)     , STAT=izl); h_ice_b        = r_init_val; ist=ist+izl
  ALLOCATE(h_ice_new_b    (nproma)     , STAT=izl); h_ice_new_b    = r_init_val; ist=ist+izl
  ALLOCATE(t_ice_b        (nproma)     , STAT=izl); t_ice_b        = r_init_val; ist=ist+izl
  ALLOCATE(t_ice_new_b    (nproma)     , STAT=izl); t_ice_new_b    = r_init_val; ist=ist+izl

  ! from the convection schemes
  ALLOCATE(fif_b          (nproma,ke)  , STAT=izl); fif_b          = r_init_val; ist=ist+izl
  ALLOCATE(fih_b          (nproma,ke+1), STAT=izl); fih_b          = r_init_val; ist=ist+izl
  ALLOCATE(pf_b           (nproma,ke)  , STAT=izl); pf_b           = r_init_val; ist=ist+izl
  ALLOCATE(ph_b           (nproma,ke+1), STAT=izl); ph_b           = r_init_val; ist=ist+izl
  ALLOCATE(dqvdt_conv_b   (nproma,ke)  , STAT=izl); dqvdt_conv_b   = r_init_val; ist=ist+izl
  IF (itype_conv == 0) THEN
    ALLOCATE(w_conv_b     (nproma,ke)  , STAT=izl); w_conv_b       = r_init_val; ist=ist+izl
    ALLOCATE(qhfl_b       (nproma)     , STAT=izl); qhfl_b         = r_init_val; ist=ist+izl
  ENDIF

  IF (itype_conv == 2) THEN
    ! from Tiedtke-Bechtold convection scheme: these fields have to be allocated 
    ! with nblock also, because they are constant through time and we do not want
    ! to compute (or copy) them all the time
    ALLOCATE(cloud_num_b  (nproma,nblock), STAT=izl); cloud_num_b    = r_init_val; ist=ist+izl
    ALLOCATE(trop_mask_b  (nproma,nblock), STAT=izl); trop_mask_b    = r_init_val; ist=ist+izl
    ALLOCATE(mtn_mask_b   (nproma,nblock), STAT=izl); mtn_mask_b     = r_init_val; ist=ist+izl
    ALLOCATE(conv_k850_b  (nproma,nblock), STAT=izl); conv_k850_b    = i_init_val; ist=ist+izl
    ALLOCATE(conv_k950_b  (nproma,nblock), STAT=izl); conv_k950_b    = i_init_val; ist=ist+izl
  ENDIF

  ALLOCATE(clc_con_b      (nproma,ke)  , STAT=izl); clc_con_b      = r_init_val; ist=ist+izl
  ALLOCATE(clw_con_b      (nproma,ke)  , STAT=izl); clw_con_b      = r_init_val; ist=ist+izl
  ALLOCATE(bas_con_b      (nproma)     , STAT=izl); bas_con_b      = r_init_val; ist=ist+izl
  ALLOCATE(top_con_b      (nproma)     , STAT=izl); top_con_b      = r_init_val; ist=ist+izl
  ALLOCATE(tt_conv_b      (nproma,ke)  , STAT=izl); tt_conv_b      = r_init_val; ist=ist+izl
  ALLOCATE(ttens_conv_b   (nproma,ke)  , STAT=izl); ttens_conv_b   = r_init_val; ist=ist+izl
  ALLOCATE(ttdiab_conv_b  (nproma,ke)  , STAT=izl); ttdiab_conv_b  = r_init_val; ist=ist+izl
  ALLOCATE(qvt_conv_b     (nproma,ke)  , STAT=izl); qvt_conv_b     = r_init_val; ist=ist+izl
  ALLOCATE(mflx_con_b     (nproma)     , STAT=izl); mflx_con_b     = r_init_val; ist=ist+izl
  ALLOCATE(cape_con_b     (nproma)     , STAT=izl); cape_con_b     = r_init_val; ist=ist+izl
  ALLOCATE(qcvg_con_b     (nproma)     , STAT=izl); qcvg_con_b     = r_init_val; ist=ist+izl

  ! this fields are used in Tiedtke (not shallow) but also in turbulence scheme
  ALLOCATE( ut_conv_b     (nproma,ke)  , STAT=izl);  ut_conv_b     = r_init_val; ist=ist+izl
  ALLOCATE( vt_conv_b     (nproma,ke)  , STAT=izl);  vt_conv_b     = r_init_val; ist=ist+izl
  ALLOCATE(tket_conv_b    (nproma,ke+1), STAT=izl); tket_conv_b    = r_init_val; ist=ist+izl

  IF (itype_conv == 0 .OR. itype_conv == 2) THEN
    ALLOCATE(prr_con_b    (nproma)     , STAT=izl); prr_con_b      = r_init_val; ist=ist+izl
    ALLOCATE(prs_con_b    (nproma)     , STAT=izl); prs_con_b      = r_init_val; ist=ist+izl
    ALLOCATE(prne_con_b   (nproma)     , STAT=izl); prne_con_b     = r_init_val; ist=ist+izl
    ALLOCATE(qct_conv_b   (nproma,ke)  , STAT=izl); qct_conv_b     = r_init_val; ist=ist+izl
    ALLOCATE(qit_conv_b   (nproma,ke)  , STAT=izl); qit_conv_b     = r_init_val; ist=ist+izl
    ALLOCATE(tke_con_b    (nproma)     , STAT=izl); tke_con_b      = r_init_val; ist=ist+izl
    ALLOCATE(vgust_con_b  (nproma)     , STAT=izl); vgust_con_b    = r_init_val; ist=ist+izl
  ENDIF

  ! miscellaneous fields
  ALLOCATE(pertstoph_b    (nproma,ke),   STAT=izl); pertstoph_b    = r_init_val; ist=ist+izl
  IF (lartif_data) THEN
  ALLOCATE(h0noise_b      (nproma)     , STAT=izl); h0noise_b      = r_init_val; ist=ist+izl
  ENDIF

  !Set microphysics pointer variable depending on calling order
  IF (lgsp_first) THEN
    t_nx_b  => t_b
    pp_nx_b => pp_b
    qv_nx_b => qv_b
    qc_nx_b => qc_b
    qi_nx_b => qi_b
    qr_nx_b => qr_b
    qs_nx_b => qs_b
    qg_nx_b => qg_b
    tt_lheat_nx_b => tt_lheat_b
  ELSE
    t_nx_b  => t_new_b
    pp_nx_b => pp_new_b
    qv_nx_b => qv_new_b
    qc_nx_b => qc_new_b
    qi_nx_b => qi_new_b
    qr_nx_b => qr_new_b
    qs_nx_b => qs_new_b
    qg_nx_b => qg_new_b
    tt_lheat_nx_b => tt_lheat_new_b
  ENDIF


!ACC_LISTX are defined in the module header
!$acc enter data           &
!$acc copyin (ACC_LIST01)  &
!$acc copyin (ACC_LIST02)  &
!$acc copyin (ACC_LIST03)  &
!$acc copyin (ACC_LIST04)  &
!$acc copyin (ACC_LIST05)  &
!$acc copyin (ACC_LIST06)  &
!$acc copyin (ACC_LIST07)  &
!$acc copyin (ACC_LIST08)  &
!$acc copyin (ACC_LIST09)  &
!$acc copyin (ACC_LIST10)  &
!$acc copyin (ACC_LIST11)  &
!$acc copyin (ACC_LIST12)  &
!$acc copyin (ACC_LIST13)  &
!$acc copyin (ACC_LIST14)  &
!$acc copyin (ACC_LIST15)  &
!$acc copyin (ACC_LIST16)  &
!$acc copyin (ACC_LIST17)  &
!$acc copyin (ACC_LIST18)  &
!$acc copyin (ACC_LIST19)  &
!$acc copyin (ACC_LIST20)  &
!$acc copyin (ACC_LIST21)  &
!$acc copyin (ACC_LIST22)  &
!$acc copyin (ACC_LIST23)  &
!$acc copyin (ACC_LIST24)  &
!$acc copyin (ACC_LIST25)  &
!$acc copyin (ACC_LIST26)  &
!$acc copyin (ACC_LIST27)  &
!$acc copyin (ACC_LIST28)  &
!$acc copyin (ACC_LIST29)  &
!$acc copyin (ACC_LIST30)  &
!$acc copyin (ACC_LIST31)  &
!$acc copyin (ACC_LIST32)  &
!$acc copyin (ACC_LIST33)  &
!$acc copyin (ACC_LIST34)  &
!$acc copyin (ACC_LIST35)  &
!$acc copyin (ACC_LIST36)  &
!$acc copyin (ACC_LIST37)  &
!$acc copyin (ACC_LIST38)  &
!$acc copyin (ACC_LIST39)  &
!$acc copyin (ACC_LIST40)  &
!$acc copyin (ACC_LIST41)  &
!$acc copyin (ACC_LIST42)  &
!$acc copyin (ACC_LIST43)  &
!$acc copyin (ACC_LIST44)

!------------------------------------------------------------------------------
! End of module procedure block_fields_allocate
!------------------------------------------------------------------------------

END SUBROUTINE block_fields_allocate

!==============================================================================
!+ Module procedure "block_fields_deallocate" in "src_block_fields_org"
!------------------------------------------------------------------------------

SUBROUTINE block_fields_deallocate( ist )

!------------------------------------------------------------------------------
!
! Description:
!   Deallocation block fields
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  INTEGER, INTENT(OUT) :: ist

! Locals:
! --------------------
  INTEGER :: izl

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine block_fields_deallocate
!------------------------------------------------------------------------------

  ist=0
  izl=0

  !ACC_LISTX are define in the module header
  !$acc exit data           &
  !$acc delete (ACC_LIST01)  &
  !$acc delete (ACC_LIST02)  &
  !$acc delete (ACC_LIST03)  &
  !$acc delete (ACC_LIST04)  &
  !$acc delete (ACC_LIST05)  &
  !$acc delete (ACC_LIST06)  &
  !$acc delete (ACC_LIST07)  &
  !$acc delete (ACC_LIST08)  &
  !$acc delete (ACC_LIST09)  &
  !$acc delete (ACC_LIST10)  &
  !$acc delete (ACC_LIST11)  &
  !$acc delete (ACC_LIST12)  &
  !$acc delete (ACC_LIST13)  &
  !$acc delete (ACC_LIST14)  &
  !$acc delete (ACC_LIST15)  &
  !$acc delete (ACC_LIST16)  &
  !$acc delete (ACC_LIST17)  &
  !$acc delete (ACC_LIST18)  &
  !$acc delete (ACC_LIST19)  &
  !$acc delete (ACC_LIST20)  &
  !$acc delete (ACC_LIST21)  &
  !$acc delete (ACC_LIST22)  &
  !$acc delete (ACC_LIST23)  &
  !$acc delete (ACC_LIST24)  &
  !$acc delete (ACC_LIST25)  &
  !$acc delete (ACC_LIST26)  &
  !$acc delete (ACC_LIST27)  &
  !$acc delete (ACC_LIST28)  &
  !$acc delete (ACC_LIST29)  &
  !$acc delete (ACC_LIST30)  &
  !$acc delete (ACC_LIST31)  &
  !$acc delete (ACC_LIST32)  &
  !$acc delete (ACC_LIST33)  &
  !$acc delete (ACC_LIST34)  &
  !$acc delete (ACC_LIST35)  &
  !$acc delete (ACC_LIST36)  &
  !$acc delete (ACC_LIST37)  &
  !$acc delete (ACC_LIST38)  &
  !$acc delete (ACC_LIST39)  &
  !$acc delete (ACC_LIST40)  &
  !$acc delete (ACC_LIST41)  &
  !$acc delete (ACC_LIST42)  &
  !$acc delete (ACC_LIST43)  &
  !$acc delete (ACC_LIST44)
!PHYS_EDIT: add new variables if necessary

  ! constant fields for the reference atmosphere
  DEALLOCATE(p0_b           , STAT=izl);     ist=ist+izl
  DEALLOCATE(p0hl_b         , STAT=izl);     ist=ist+izl
  DEALLOCATE(rho0_b         , STAT=izl);     ist=ist+izl
  DEALLOCATE(dp0_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(hhl_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(dz_b           , STAT=izl);     ist=ist+izl

  ! external parameters, needed in several parameterizations
  DEALLOCATE(hsurf_b        , STAT=izl);     ist=ist+izl
  DEALLOCATE(fr_land_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(llandmask_b    , STAT=izl);     ist=ist+izl
  DEALLOCATE(llakemask_b    , STAT=izl);     ist=ist+izl
  DEALLOCATE(isoiltyp_b     , STAT=izl);     ist=ist+izl
  DEALLOCATE(plcov_b        , STAT=izl);     ist=ist+izl
  DEALLOCATE(rootdp_b       , STAT=izl);     ist=ist+izl
  IF (itype_vdif > -2) THEN
  DEALLOCATE(sso_stdh_b     , STAT=izl);     ist=ist+izl
  ENDIF
  DEALLOCATE(sai_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(tai_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(eai_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(l_pat_b        , STAT=izl);     ist=ist+izl
  DEALLOCATE(rsmin2d_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(depth_lk_b     , STAT=izl);     ist=ist+izl
  DEALLOCATE(fr_lake_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(fc_b           , STAT=izl);     ist=ist+izl

  ! dynamical variables
  DEALLOCATE(ptot_b         , STAT=izl);     ist=ist+izl
  DEALLOCATE(pp_b           , STAT=izl);     ist=ist+izl
  DEALLOCATE(t_b            , STAT=izl);     ist=ist+izl
  DEALLOCATE(u_m_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(v_m_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(w_b            , STAT=izl);     ist=ist+izl
  DEALLOCATE(rho_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(qv_b           , STAT=izl);     ist=ist+izl
  DEALLOCATE(qc_b           , STAT=izl);     ist=ist+izl
  DEALLOCATE(qi_b           , STAT=izl);     ist=ist+izl
  DEALLOCATE(qr_b           , STAT=izl);     ist=ist+izl
  DEALLOCATE(qs_b           , STAT=izl);     ist=ist+izl
  DEALLOCATE(qg_b           , STAT=izl);     ist=ist+izl
  DEALLOCATE(tke_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(t_new_b        , STAT=izl);     ist=ist+izl
  DEALLOCATE(pp_new_b       , STAT=izl);     ist=ist+izl
  DEALLOCATE(qv_new_b       , STAT=izl);     ist=ist+izl
  DEALLOCATE(qc_new_b       , STAT=izl);     ist=ist+izl
  DEALLOCATE(qi_new_b       , STAT=izl);     ist=ist+izl
  DEALLOCATE(qr_new_b       , STAT=izl);     ist=ist+izl
  DEALLOCATE(qs_new_b       , STAT=izl);     ist=ist+izl
  DEALLOCATE(qg_new_b       , STAT=izl);     ist=ist+izl

  DEALLOCATE(utens_b        , STAT=izl);     ist=ist+izl
  DEALLOCATE(vtens_b        , STAT=izl);     ist=ist+izl
  DEALLOCATE(ttens_b        , STAT=izl);     ist=ist+izl
  DEALLOCATE(qvtens_b       , STAT=izl);     ist=ist+izl
  DEALLOCATE(qctens_b       , STAT=izl);     ist=ist+izl
  DEALLOCATE(qitens_b       , STAT=izl);     ist=ist+izl
  DEALLOCATE(tketens_b      , STAT=izl);     ist=ist+izl

  ! fields for surface values
  DEALLOCATE(ps_b           , STAT=izl);     ist=ist+izl
  DEALLOCATE(t_s_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(t_s_new_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(t_g_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(t_g_new_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(qv_s_b         , STAT=izl);     ist=ist+izl
  DEALLOCATE(qv_s_new_b     , STAT=izl);     ist=ist+izl

  ! from microphysics
  DEALLOCATE(tinc_lh_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(qrs_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(prr_gsp_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(prs_gsp_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(prg_gsp_b      , STAT=izl);     ist=ist+izl
#ifdef NUDGING
  DEALLOCATE(tt_lheat_b     , STAT=izl);     ist=ist+izl
  DEALLOCATE(tt_lheat_new_b , STAT=izl);     ist=ist+izl
  DEALLOCATE(qrsflux_b      , STAT=izl);     ist=ist+izl
#endif

  ! from turbulence scheme
  DEALLOCATE(l_hori_b       , STAT=izl);     ist=ist+izl
  DEALLOCATE(gz0_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(tcm_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(tch_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(tfm_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(tfh_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(tfv_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(tvm_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(tvh_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(tkr_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(tkred_sfc_b    , STAT=izl);     ist=ist+izl
  DEALLOCATE(rcld_b         , STAT=izl);     ist=ist+izl
  DEALLOCATE(tkvm_b         , STAT=izl);     ist=ist+izl
  DEALLOCATE(tkvh_b         , STAT=izl);     ist=ist+izl
  DEALLOCATE(tkhm_b         , STAT=izl);     ist=ist+izl
  DEALLOCATE(tkhh_b         , STAT=izl);     ist=ist+izl
  DEALLOCATE(hdef2_b        , STAT=izl);     ist=ist+izl
  DEALLOCATE(hdiv_b         , STAT=izl);     ist=ist+izl
  DEALLOCATE(dwdx_b         , STAT=izl);     ist=ist+izl
  DEALLOCATE(dwdy_b         , STAT=izl);     ist=ist+izl
  DEALLOCATE(edr_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(tket_sso_b     , STAT=izl);     ist=ist+izl
  DEALLOCATE(tket_hshr_b    , STAT=izl);     ist=ist+izl
  DEALLOCATE(tket_adv_b     , STAT=izl);     ist=ist+izl
  DEALLOCATE(dqvdt_b        , STAT=izl);     ist=ist+izl
  DEALLOCATE(ut_turb_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(vt_turb_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(ut_sso_b       , STAT=izl);     ist=ist+izl
  DEALLOCATE(vt_sso_b       , STAT=izl);     ist=ist+izl

  DEALLOCATE(t_2m_b         , STAT=izl);     ist=ist+izl
  DEALLOCATE(qv_2m_b        , STAT=izl);     ist=ist+izl
  DEALLOCATE(td_2m_b        , STAT=izl);     ist=ist+izl
  DEALLOCATE(rh_2m_b        , STAT=izl);     ist=ist+izl
  DEALLOCATE(u_10m_b        , STAT=izl);     ist=ist+izl
  DEALLOCATE(v_10m_b        , STAT=izl);     ist=ist+izl
  DEALLOCATE(shfl_s_b       , STAT=izl);     ist=ist+izl
  DEALLOCATE(lhfl_s_b       , STAT=izl);     ist=ist+izl
  DEALLOCATE(qvsflx_b       , STAT=izl);     ist=ist+izl

  ! from the surface schemes
  DEALLOCATE(t_snow_b       , STAT=izl);     ist=ist+izl
  DEALLOCATE(t_snow_new_b   , STAT=izl);     ist=ist+izl
  DEALLOCATE(w_snow_b       , STAT=izl);     ist=ist+izl
  DEALLOCATE(w_snow_new_b   , STAT=izl);     ist=ist+izl
  DEALLOCATE(h_snow_b       , STAT=izl);     ist=ist+izl
  DEALLOCATE(h_snow_new_b   , STAT=izl);     ist=ist+izl
  DEALLOCATE(rho_snow_b     , STAT=izl);     ist=ist+izl
  DEALLOCATE(rho_snow_new_b , STAT=izl);     ist=ist+izl
  DEALLOCATE(freshsnow_b    , STAT=izl);     ist=ist+izl
  DEALLOCATE(fr_snow_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(w_i_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(w_i_new_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(w_p_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(w_p_new_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(w_s_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(w_s_new_b      , STAT=izl);     ist=ist+izl

  DEALLOCATE(t_so_b         , STAT=izl);     ist=ist+izl
  DEALLOCATE(t_so_new_b     , STAT=izl);     ist=ist+izl
  DEALLOCATE(w_so_b         , STAT=izl);     ist=ist+izl
  DEALLOCATE(w_so_new_b     , STAT=izl);     ist=ist+izl
  DEALLOCATE(w_so_ice_b     , STAT=izl);     ist=ist+izl
  DEALLOCATE(w_so_ice_new_b , STAT=izl);     ist=ist+izl

  DEALLOCATE(t_snow_mult_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(t_snow_mult_new_b  , STAT=izl);     ist=ist+izl
  DEALLOCATE(w_snow_mult_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(w_snow_mult_new_b  , STAT=izl);     ist=ist+izl
  DEALLOCATE(wliq_snow_b        , STAT=izl);     ist=ist+izl
  DEALLOCATE(wliq_snow_new_b    , STAT=izl);     ist=ist+izl
  DEALLOCATE(rho_snow_mult_b    , STAT=izl);     ist=ist+izl
  DEALLOCATE(rho_snow_mult_new_b, STAT=izl);     ist=ist+izl
  DEALLOCATE(dzh_snow_mult_b    , STAT=izl);     ist=ist+izl
  DEALLOCATE(dzh_snow_mult_new_b, STAT=izl);     ist=ist+izl

  DEALLOCATE(runoff_s_b     , STAT=izl);     ist=ist+izl
  DEALLOCATE(runoff_g_b     , STAT=izl);     ist=ist+izl
  DEALLOCATE(lhfl_bs_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(lhfl_pl_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(rstom_b        , STAT=izl);     ist=ist+izl

  ! from flake or seaice scheme
  IF (llake) THEN
  DEALLOCATE(qmomflux_b     , STAT=izl);     ist=ist+izl
  DEALLOCATE(fetch_lk_b     , STAT=izl);     ist=ist+izl
  DEALLOCATE(dp_bs_lk_b     , STAT=izl);     ist=ist+izl
  DEALLOCATE(t_bs_lk_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(gamso_lk_b     , STAT=izl);     ist=ist+izl
  DEALLOCATE(t_mnw_lk_b     , STAT=izl);     ist=ist+izl
  DEALLOCATE(t_mnw_lk_new_b , STAT=izl);     ist=ist+izl
  DEALLOCATE(t_wml_lk_b     , STAT=izl);     ist=ist+izl
  DEALLOCATE(t_wml_lk_new_b , STAT=izl);     ist=ist+izl
  DEALLOCATE(t_bot_lk_b     , STAT=izl);     ist=ist+izl
  DEALLOCATE(t_bot_lk_new_b , STAT=izl);     ist=ist+izl
  DEALLOCATE(t_b1_lk_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(t_b1_lk_new_b  , STAT=izl);     ist=ist+izl
  DEALLOCATE(c_t_lk_b       , STAT=izl);     ist=ist+izl
  DEALLOCATE(c_t_lk_new_b   , STAT=izl);     ist=ist+izl
  DEALLOCATE(h_ml_lk_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(h_ml_lk_new_b  , STAT=izl);     ist=ist+izl
  DEALLOCATE(h_b1_lk_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(h_b1_lk_new_b  , STAT=izl);     ist=ist+izl
  ENDIF

  DEALLOCATE(h_ice_b        , STAT=izl);     ist=ist+izl
  DEALLOCATE(h_ice_new_b    , STAT=izl);     ist=ist+izl
  DEALLOCATE(t_ice_b        , STAT=izl);     ist=ist+izl
  DEALLOCATE(t_ice_new_b    , STAT=izl);     ist=ist+izl

  ! from the convection schemes
  DEALLOCATE(fif_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(fih_b          , STAT=izl);     ist=ist+izl
  DEALLOCATE(pf_b           , STAT=izl);     ist=ist+izl
  DEALLOCATE(ph_b           , STAT=izl);     ist=ist+izl
  DEALLOCATE(dqvdt_conv_b   , STAT=izl);     ist=ist+izl
  IF (itype_conv == 0) THEN
    DEALLOCATE(w_conv_b     , STAT=izl);     ist=ist+izl
    DEALLOCATE(qhfl_b       , STAT=izl);     ist=ist+izl
  ENDIF

  IF (itype_conv == 2) THEN
    ! from Tiedtke-Bechtold convection scheme
    DEALLOCATE(cloud_num_b  , STAT=izl);     ist=ist+izl
    DEALLOCATE(trop_mask_b  , STAT=izl);     ist=ist+izl
    DEALLOCATE(mtn_mask_b   , STAT=izl);     ist=ist+izl
    DEALLOCATE(conv_k850_b  , STAT=izl);     ist=ist+izl
    DEALLOCATE(conv_k950_b  , STAT=izl);     ist=ist+izl
  ENDIF

  DEALLOCATE(clc_con_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(clw_con_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(bas_con_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(top_con_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(tt_conv_b      , STAT=izl);     ist=ist+izl
  DEALLOCATE(ttens_conv_b   , STAT=izl);     ist=ist+izl
  DEALLOCATE(ttdiab_conv_b  , STAT=izl);     ist=ist+izl
  DEALLOCATE(qvt_conv_b     , STAT=izl);     ist=ist+izl
  DEALLOCATE(mflx_con_b     , STAT=izl);     ist=ist+izl
  DEALLOCATE(cape_con_b     , STAT=izl);     ist=ist+izl
  DEALLOCATE(qcvg_con_b     , STAT=izl);     ist=ist+izl

  ! this fields are used in Tiedtke (not shallow) but also in turbulence scheme
  DEALLOCATE( ut_conv_b     , STAT=izl);     ist=ist+izl
  DEALLOCATE( vt_conv_b     , STAT=izl);     ist=ist+izl
  DEALLOCATE(tket_conv_b    , STAT=izl);     ist=ist+izl

  IF (itype_conv == 0) THEN
    DEALLOCATE(prr_con_b    , STAT=izl);     ist=ist+izl
    DEALLOCATE(prs_con_b    , STAT=izl);     ist=ist+izl
    DEALLOCATE(prne_con_b   , STAT=izl);     ist=ist+izl
    DEALLOCATE(qct_conv_b   , STAT=izl);     ist=ist+izl
    DEALLOCATE(qit_conv_b   , STAT=izl);     ist=ist+izl
    DEALLOCATE(tke_con_b    , STAT=izl);     ist=ist+izl
    DEALLOCATE(vgust_con_b  , STAT=izl);     ist=ist+izl
  ENDIF

  ! miscellaneous fields
  DEALLOCATE(pertstoph_b    , STAT=izl);     ist=ist+izl
  IF (lartif_data) THEN
  DEALLOCATE(h0noise_b      , STAT=izl);     ist=ist+izl
  ENDIF

!------------------------------------------------------------------------------
! End of module procedure block_fields_deallocate
!------------------------------------------------------------------------------

END SUBROUTINE block_fields_deallocate

!==============================================================================
!+ Module procedure "register_block_fields_all" in "src_block_fields_org"
!------------------------------------------------------------------------------

SUBROUTINE block_fields_register_all( ierror )

!------------------------------------------------------------------------------
!
! Description:
!   Register all block fields with corresponding i,j,k field
!   name convention:  "t_b" corresponds to time level nnow (or nold for l-f)
!                     "t_new_b" to time level nnew
!   tracer need to be registered using their tracer
!   index
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------
  INTEGER, INTENT(OUT) :: ierror
INTEGER, POINTER :: pnx

!----------- End of header ----------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine block_fields_register_all
!------------------------------------------------------------------------------

  ierror = 0     

  ! select time level according to the integration scheme used
  IF ( l2tls ) THEN
    pnx  => nnow
  ELSE
    pnx  => nold
  ENDIF

  !PHYS_EDIT: add new variables if necessary

  ! constant fields for the reference atmosphere
  CALL register_block_field ("p0"            , p0             , p0_b           )
  CALL register_block_field ("p0hl"          , p0hl           , p0hl_b         )
  CALL register_block_field ("rho0"          , rho0           , rho0_b         )
  CALL register_block_field ("dp0"           , dp0            , dp0_b          )
  CALL register_block_field ("hhl"           , hhl            , hhl_b          )
  ! about dz_b:  is computed in gscp_interface on the fly

  ! external parameters, needed in several parameterizations
  CALL register_block_field ("hsurf"         , hhl            , hhl_b          )
  CALL register_block_field ("fr_land"       , fr_land        , fr_land_b      )
  CALL register_block_field ("llandmask"     , llandmask      , llandmask_b    )
  CALL register_block_field ("llakemask"     , llakemask      , llakemask_b    )
  CALL register_block_field ("isoiltyp"      , isoiltyp       , isoiltyp_b     )
  CALL register_block_field ("plcov"         , plcov          , plcov_b        )
  CALL register_block_field ("rootdp"        , rootdp         , rootdp_b       )
  IF (itype_vdif > -2) THEN
  CALL register_block_field ("sso_stdh"      , sso_stdh       , sso_stdh_b     )
  ENDIF
  CALL register_block_field ("sai"           , sai            , sai_b          )
  CALL register_block_field ("tai"           , tai            , tai_b          )
  CALL register_block_field ("eai"           , eai            , eai_b          )
  CALL register_block_field ("l_pat"         , l_pat          , l_pat_b        )
  CALL register_block_field ("rsmin2d"       , rsmin2d        , rsmin2d_b      )
  CALL register_block_field ("depth_lk"      , depth_lk       , depth_lk_b     )
  CALL register_block_field ("fr_lake"       , fr_lake        , fr_lake_b      )
  CALL register_block_field ("fc"            , fc             , fc_b           )

  ! dynamical variables (time nnow)
  ! about ptot: is computed in gscp_interface, turb_interface on the fly
  CALL register_block_field ("pp"            , pp             , pp_b           , pnx)
  CALL register_block_field ("t"             , t              , t_b            , pnx)
  CALL register_block_field ("u_m"           , u_m            , u_m_b          )
  CALL register_block_field ("v_m"           , v_m            , v_m_b          )
  CALL register_block_field ("w"             , w              , w_b            , pnx)
  CALL register_block_field ("rho"           , rho            , rho_b          )
  CALL register_block_field_tracer ("qv"     , idt_qv         , qv_b           , pnx)
  CALL register_block_field_tracer ("qc"     , idt_qc         , qc_b           , pnx)
  IF(idt_qi>=0) CALL register_block_field_tracer ("qi", idt_qi, qi_b           , pnx)
  IF(idt_qr>=0) CALL register_block_field_tracer ("qr", idt_qr, qr_b           , pnx)
  IF(idt_qs>=0) CALL register_block_field_tracer ("qs", idt_qs, qs_b           , pnx)
  IF(idt_qg>=0) CALL register_block_field_tracer ("qg", idt_qg, qg_b           , pnx)
!XL_TODO  CALL register_block_field("tke",tke,tke_b,?) tke is not supported : timelevel

  !Time new
  CALL register_block_field ("t_new"         , t              , t_new_b        , nnew)
  CALL register_block_field ("pp_new"        , pp             , pp_new_b       , nnew)
  CALL register_block_field_tracer ("qv_new" , idt_qv         , qv_new_b       , nnew)
  CALL register_block_field_tracer ("qc_new" , idt_qc         , qc_new_b       , nnew)
  IF(idt_qi>=0) CALL register_block_field_tracer ("qi_new",idt_qi,qi_new_b     , nnew)
  IF(idt_qr>=0) CALL register_block_field_tracer ("qr_new",idt_qr,qr_new_b     , nnew)
  IF(idt_qs>=0) CALL register_block_field_tracer ("qs_new",idt_qs,qs_new_b     , nnew)
  IF(idt_qg>=0) CALL register_block_field_tracer ("qg_new",idt_qg,qg_new_b     , nnew)

  CALL register_block_field ("utens"         , utens          , utens_b        )
  CALL register_block_field ("vtens"         , vtens          , vtens_b        )
  CALL register_block_field ("tttens"        , ttens          , ttens_b        )
  IF (lalloc_tke) CALL register_block_field ("tketens"       , tketens        , tketens_b      )
  CALL register_block_field_tracer_tens ("qvtens"  , idt_qv   , qvtens_b       )
  CALL register_block_field_tracer_tens ("qctens"  , idt_qc   , qctens_b       )
  CALL register_block_field_tracer_tens ("qitens"  , idt_qi   , qitens_b       )

  ! fields for surface values
  CALL register_block_field ("ps"            , ps             , ps_b           , pnx )
  CALL register_block_field ("t_s"           , t_s            , t_s_b          , pnx )
  CALL register_block_field ("t_s_new"       , t_s            , t_s_new_b      , nnew)
  CALL register_block_field ("t_g"           , t_g            , t_g_b          , pnx )
  CALL register_block_field ("t_g_new"       , t_g            , t_g_new_b      , nnew)
  CALL register_block_field ("qv_s"          , qv_s           , qv_s_b         , pnx )
  CALL register_block_field ("qv_s_new"      , qv_s           , qv_s_new_b     , nnew)

  ! from microphysics
  IF (ALLOCATED(tinc_lh))   CALL register_block_field ("tinc_lh"  , tinc_lh  , tinc_lh_b)
  CALL register_block_field ("qrs"           , qrs            , qrs_b          )
  CALL register_block_field ("prr_gsp"       , prr_gsp        , prr_gsp_b      )
  IF (idt_qs>=0)            CALL register_block_field ("prs_gsp"  , prs_gsp  , prs_gsp_b)
  IF (ALLOCATED(prg_gsp))   CALL register_block_field ("prg_gsp"  , prg_gsp  , prg_gsp_b)
#ifdef NUDGING
  IF (ALLOCATED(tt_lheat))  CALL register_block_field ("tt_lheat"    , tt_lheat, tt_lheat_b    , pnx)
  IF (ALLOCATED(tt_lheat))  CALL register_block_field ("tt_lheat_new", tt_lheat, tt_lheat_new_b, nnew)
  IF (ALLOCATED(qrsflux))   CALL register_block_field ("qrsflux"     , qrsflux , qrsflux_b     )
#endif

  ! from turbulence scheme
  CALL register_block_field ("l_hori"        , l_hori         , l_hori_b       )
  CALL register_block_field ("gz0"           , gz0            , gz0_b          )
  CALL register_block_field ("tcm"           , tcm            , tcm_b          )
  CALL register_block_field ("tch"           , tch            , tch_b          )
  CALL register_block_field ("tfm"           , tfm            , tfm_b          )
  CALL register_block_field ("tfh"           , tfh            , tfh_b          )
  CALL register_block_field ("tfv"           , tfv            , tfv_b          )
  CALL register_block_field ("tvm"           , tvm            , tvm_b          )
  CALL register_block_field ("tvh"           , tvh            , tvh_b          )
  CALL register_block_field ("tkr"           , tkr            , tkr_b          )
  CALL register_block_field ("tkred_sfc"     , tkred_sfc      , tkred_sfc_b    )
  IF (lalloc_tke) CALL register_block_field ("rcld"          , rcld           , rcld_b         )
  CALL register_block_field ("tkvm"          , tkvm           , tkvm_b         )
  CALL register_block_field ("tkvh"          , tkvh           , tkvh_b         )
  CALL register_block_field ("tkhm"          , tkhm           , tkhm_b         )
  CALL register_block_field ("tkhh"          , tkhh           , tkhh_b         )
  CALL register_block_field ("hdef2"         , hdef2          , hdef2_b        )
  CALL register_block_field ("hdiv"          , hdiv           , hdiv_b         )
  CALL register_block_field ("dwdx"          , dwdx           , dwdx_b         )
  CALL register_block_field ("dwdy"          , dwdy           , dwdy_b         )
  IF (lalloc_tke) CALL register_block_field ("edr"           , edr            , edr_b          )
  CALL register_block_field ("tket_sso"      , tket_sso       , tket_sso_b     )
  CALL register_block_field ("tket_hshr"     , tket_hshr      , tket_hshr_b    )
  IF (lalloc_tke) CALL register_block_field ("tket_adv"      , tket_adv       , tket_adv_b     )
  CALL register_block_field ("dqvdt"         , dqvdt          , dqvdt_b        )
  CALL register_block_field ("ut_turb"       , ut_turb        , ut_turb_b      )
  CALL register_block_field ("vt_turb"       , vt_turb        , vt_turb_b      )
  CALL register_block_field ("ut_sso"        , ut_sso         , ut_sso_b       )
  CALL register_block_field ("vt_sso"        , vt_sso         , vt_sso_b       )

  CALL register_block_field ("t_2m"          , t_2m           , t_2m_b         )
  CALL register_block_field ("qv_2m"         , qv_2m          , qv_2m_b        )
  CALL register_block_field ("td_2m"         , td_2m          , td_2m_b        )
  CALL register_block_field ("rh_2m"         , rh_2m          , rh_2m_b        )
  CALL register_block_field ("u_10m"         , u_10m          , u_10m_b        )
  CALL register_block_field ("v_10m"         , v_10m          , v_10m_b        )
  CALL register_block_field ("shfl_s"        , shfl_s         , shfl_s_b       )
  CALL register_block_field ("lhfl_s"        , lhfl_s         , lhfl_s_b       )
  CALL register_block_field ("qvsflx"        , qvsflx         , qvsflx_b       )

  ! from the surface schemes
  CALL register_block_field ("t_snow"        , t_snow         , t_snow_b       , pnx )
  CALL register_block_field ("t_snow_new"    , t_snow         , t_snow_new_b   , nnew)
  CALL register_block_field ("w_snow"        , w_snow         , w_snow_b       , pnx )
  CALL register_block_field ("w_snow_new"    , w_snow         , w_snow_new_b   , nnew)
  CALL register_block_field ("h_snow"        , h_snow         , h_snow_b       , pnx )
  CALL register_block_field ("h_snow_new"    , h_snow         , h_snow_new_b   , nnew)
  CALL register_block_field ("rho_snow"      , rho_snow       , rho_snow_b     , pnx )
  CALL register_block_field ("rho_snow_new"  , rho_snow       , rho_snow_new_b , nnew)
  CALL register_block_field ("freshsnow"     , freshsnow      , freshsnow_b    )
  CALL register_block_field ("fr_snow"       , fr_snow        , fr_snow_b      )
  CALL register_block_field ("w_i"           , w_i            , w_i_b          , pnx )
  CALL register_block_field ("w_i_new"       , w_i            , w_i_new_b      , nnew)
  CALL register_block_field ("w_p"           , w_p            , w_p_b          , pnx )
  CALL register_block_field ("w_p_new"       , w_p            , w_p_new_b      , nnew)
  CALL register_block_field ("w_s"           , w_s            , w_s_b          , pnx )
  CALL register_block_field ("w_s_new"       , w_s            , w_s_new_b      , nnew)

  CALL register_block_field ("t_so"          , t_so           , t_so_b         , pnx )
  CALL register_block_field ("t_so_new"      , t_so           , t_so_new_b     , nnew)
  CALL register_block_field ("w_so"          , w_so           , w_so_b         , pnx )
  CALL register_block_field ("w_so_new"      , w_so           , w_so_new_b     , nnew)
  CALL register_block_field ("w_so_ice"      , w_so_ice       , w_so_ice_b     , pnx )
  CALL register_block_field ("w_so_ice_new"  , w_so_ice       , w_so_ice_new_b , nnew)

  CALL register_block_field ("t_snow_mult"      , t_snow_mult    , t_snow_mult_b       , pnx )
  CALL register_block_field ("t_snow_mult_new"  , t_snow_mult    , t_snow_mult_new_b   , nnew)
  CALL register_block_field ("w_snow_mult"      , w_snow_mult    , w_snow_mult_b       , pnx )
  CALL register_block_field ("w_snow_mult_new"  , w_snow_mult    , w_snow_mult_new_b   , nnew)
  CALL register_block_field ("wliq_snow"        , wliq_snow      , wliq_snow_b         , pnx )
  CALL register_block_field ("wliq_snow_new"    , wliq_snow      , wliq_snow_new_b     , nnew)
  CALL register_block_field ("rho_snow_mult"    , rho_snow_mult  , rho_snow_mult_b     , pnx )
  CALL register_block_field ("rho_snow_mult_new", rho_snow_mult  , rho_snow_mult_new_b , nnew)
  CALL register_block_field ("dzh_snow_mult"    , dzh_snow_mult  , dzh_snow_mult_b     , pnx )
  CALL register_block_field ("dzh_snow_mult_new", dzh_snow_mult  , dzh_snow_mult_new_b , nnew)

  CALL register_block_field ("runoff_s"      , runoff_s       , runoff_s_b     )
  CALL register_block_field ("runoff_g"      , runoff_g       , runoff_g_b     )
  CALL register_block_field ("lhfl_bs"       , lhfl_bs        , lhfl_bs_b      )
  CALL register_block_field ("lhfl_pl"       , lhfl_pl        , lhfl_pl_b      )
  CALL register_block_field ("rstom"         , rstom          , rstom_b        )

  ! needed from radiation
  CALL register_block_field ("sobs"          , sobs           , sobs_b         )
  CALL register_block_field ("thbs"          , thbs           , thbs_b         )
  CALL register_block_field ("pabs"          , pabs           , pabs_b         )

  ! from flake and seaice scheme
  IF (llake) THEN
  CALL register_block_field ("qmomflux"      , qmomflux       , qmomflux_b     )
  CALL register_block_field ("fetch_lk"      , fetch_lk       , fetch_lk_b     )
  CALL register_block_field ("dp_bs_lk"      , dp_bs_lk       , dp_bs_lk_b     )
  CALL register_block_field ("t_bs_lk"       , t_bs_lk        , t_bs_lk_b      )
  CALL register_block_field ("gamso_lk"      , gamso_lk       , gamso_lk_b     )

  CALL register_block_field ("t_mnw_lk"      , t_mnw_lk       , t_mnw_lk_b     , pnx )
  CALL register_block_field ("t_mnw_lk_new"  , t_mnw_lk       , t_mnw_lk_new_b , nnew)
  CALL register_block_field ("t_wml_lk"      , t_wml_lk       , t_wml_lk_b     , pnx )
  CALL register_block_field ("t_wml_lk_new"  , t_wml_lk       , t_wml_lk_new_b , nnew)
  CALL register_block_field ("t_bot_lk"      , t_bot_lk       , t_bot_lk_b     , pnx )
  CALL register_block_field ("t_bot_lk_new"  , t_bot_lk       , t_bot_lk_new_b , nnew)
  CALL register_block_field ("t_b1_lk"       , t_b1_lk        , t_b1_lk_b      , pnx )
  CALL register_block_field ("t_b1_lk_new"   , t_b1_lk        , t_b1_lk_new_b  , nnew)
  CALL register_block_field ("c_t_lk"        , c_t_lk         , c_t_lk_b       , pnx )
  CALL register_block_field ("c_t_lk_new"    , c_t_lk         , c_t_lk_new_b   , nnew)
  CALL register_block_field ("h_ml_lk"       , h_ml_lk        , h_ml_lk_b      , pnx )
  CALL register_block_field ("h_ml_lk_new"   , h_ml_lk        , h_ml_lk_new_b  , nnew)
  CALL register_block_field ("h_b1_lk"       , h_b1_lk        , h_b1_lk_b      , pnx )
  CALL register_block_field ("h_b1_lk_new"   , h_b1_lk        , h_b1_lk_new_b  , nnew)
  ENDIF

  IF (ALLOCATED( t_ice ))   CALL register_block_field("t_ice"    , t_ice , t_ice_b     , pnx )
  IF (ALLOCATED( t_ice ))   CALL register_block_field("t_ice_new", t_ice , t_ice_new_b , nnew)
  IF (ALLOCATED( h_ice ))   CALL register_block_field("h_ice"    , h_ice , h_ice_b     , pnx )
  IF (ALLOCATED( h_ice ))   CALL register_block_field("h_ice_new", h_ice , h_ice_new_b , nnew)

  ! from the convection schemes
  CALL register_block_field ("dqvdt_conv"    , dqvdt_conv     , dqvdt_conv_b   )
  IF (itype_conv == 0) THEN
  CALL register_block_field ("w_conv"        , w_conv         , w_conv_b       )
  CALL register_block_field ("qhfl"          , qhfl           , qhfl_b         )
  ENDIF
  CALL register_block_field ("bas_con"       , bas_con        , bas_con_b      )
  CALL register_block_field ("top_con"       , top_con        , top_con_b      )
  CALL register_block_field ("clc_con"       , clc_con        , clc_con_b      )
  CALL register_block_field ("clw_con"       , clw_con        , clw_con_b      )
  CALL register_block_field ("tt_conv"       , tt_conv        , tt_conv_b      )
  CALL register_block_field ("ttens_conv"    , ttens_conv     , ttens_conv_b   )
  CALL register_block_field ("ttdiab_conv"   , ttdiab_conv    , ttdiab_conv_b  )
  CALL register_block_field ("qvt_conv"      , qvt_conv       , qvt_conv_b     )
  CALL register_block_field ("mflx_con"      , mflx_con       , mflx_con_b     )
  CALL register_block_field ("cape_con"      , cape_con       , cape_con_b     )
  CALL register_block_field ("qcvg_con"      , qcvg_con       , qcvg_con_b     )

  ! this fields are used in Tiedtke (not shallow) but also in turbulence scheme
  CALL register_block_field ("ut_conv"       , ut_conv        , ut_conv_b      )
  CALL register_block_field ("vt_conv"       , vt_conv        , vt_conv_b      )
  CALL register_block_field ("tket_conv"     , tket_conv      , tket_conv_b    )

  IF (itype_conv == 0 .OR. itype_conv == 2) THEN
  CALL register_block_field ("prr_con"       , prr_con        , prr_con_b      )
  CALL register_block_field ("prs_con"       , prs_con        , prs_con_b      )
  CALL register_block_field ("prne_con"      , prne_con       , prne_con_b     )
  CALL register_block_field ("qct_conv"      , qct_conv       , qct_conv_b     )
  CALL register_block_field ("qit_conv"      , qit_conv       , qit_conv_b     )
  CALL register_block_field ("tke_con"       , tke_con        , tke_con_b      )
  CALL register_block_field ("vgust_con"     , vgust_con      , vgust_con_b    )
  ENDIF

  ! miscellaneous fields
  IF (ALLOCATED(pertstoph)) CALL register_block_field ("pertstoph"   , pertstoph  , pertstoph_b   )
  IF (lartif_data)          CALL register_block_field ("h0noise"     , h0noise    , h0noise_b      )

CONTAINS 

!==============================================================================
! wrapper routine to register tracer field using trcr_get_block
! as for register_block_field, this routine may abort
!------------------------------------------------------------------------------

  SUBROUTINE register_block_field_tracer(yname,idt_trcr,field_b, ntime)

    ! Subroutine arguments:
    ! --------------------

      CHARACTER (LEN=*), INTENT(IN)  :: yname
      INTEGER, INTENT(IN)   ::    idt_trcr  
      REAL(KIND=wp), TARGET, INTENT(IN)   :: field_b(:,:)
      INTEGER, TARGET, INTENT(IN) :: ntime 

    ! Locals :
    ! --------

       REAL(KIND=wp), POINTER :: pdata_trcr(:,:,:,:,:)
       INTEGER :: izerror
       CHARACTER (LEN=250)       :: yzerrmsg

    !----------- End of header ------------------------------------------------

    CALL trcr_get_block(izerror, idt_trcr, idt_trcr,  ptr_wtlev=pdata_trcr)

    IF ( izerror /= 0 ) THEN
      yzerrmsg = trcr_errorstr( izerror )
      CALL model_abort( my_cart_id, izerror, yzerrmsg, 'trcr_init' )
    ENDIF

    !the trcr_get_block routine returns a 5d pointer. Since we call it with
    !argument idt_trcr, idt_trcr the 4th index (tracer) is from 1 to 1. We
    !therefore need to pass pdata_trcr(:,:,:,1,:)
    CALL register_block_field(yname,pdata_trcr(:,:,:,1,:),field_b,ntime)

  END SUBROUTINE register_block_field_tracer

!==============================================================================
! wrapper routine to register tendency tracer field using trcr_get_block
! as for register_block_field, this routine may abort
!------------------------------------------------------------------------------

  SUBROUTINE register_block_field_tracer_tens(name,idt_trcr,field_b)

    ! Subroutine arguments:
    ! --------------------

      CHARACTER (LEN=*), INTENT(IN)    :: name
      INTEGER,           INTENT(IN)    :: idt_trcr
      REAL(KIND=wp), TARGET, INTENT(IN):: field_b(:,:)

    ! Locals :
    ! --------

      REAL(KIND=wp), POINTER :: p_tens(:,:,:,:)
      INTEGER :: izerror
      CHARACTER (LEN=250)       :: yzerrmsg

    !----------- End of header ------------------------------------------------

    CALL trcr_get_block(izerror, idt_trcr, idt_trcr,  ptr_tens=p_tens)

    IF ( izerror /= 0 ) THEN
      yzerrmsg = trcr_errorstr( izerror )
      CALL model_abort( my_cart_id, izerror, yzerrmsg, 'register_block_field_tracer' )
    ENDIF

    !the trcr_get_block routine returns a 4d pointer. Since we call it with
    !argument idt_trcr, idt_trcr the 3rd index (tracer) is from 1 to 1. We
    !therefore need to pass pdata_trcr(:,:,:,1)
    CALL register_block_field(name,p_tens(:,:,:,1),field_b)

  END SUBROUTINE register_block_field_tracer_tens

!==============================================================================

!------------------------------------------------------------------------------
! End of module procedure block_fields_register_all
!------------------------------------------------------------------------------

END SUBROUTINE block_fields_register_all

!==============================================================================
!==============================================================================
!+ Module procedure "register_block_fields_all" in "src_block_fields_org" 
!------------------------------------------------------------------------------

SUBROUTINE block_fields_cleanup( ierror )

!------------------------------------------------------------------------------
!
! Description:
!   Deallocation block fields and cleanup data structure
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------
  INTEGER, INTENT(OUT) :: ierror

!----------- End of header ----------------------------------------------------

  ierror=0

  CALL finalize_block_fields(ierror)

  CALL block_fields_deallocate(ierror)

END SUBROUTINE block_fields_cleanup

!==============================================================================
!==============================================================================
!+ Module procedure "block_fields_copytoblock_tke" in "src_block_fields_org"
!------------------------------------------------------------------------------

SUBROUTINE block_fields_copytoblock_tke(ipend,ib)

  INTEGER, INTENT(IN) :: ipend,ib
  INTEGER :: i,j,k,ip, itim, ntim

  ntim=SIZE(tke,4)

  ! not before turbulence is running on GPUs
  !xacc data present(tke, tke_b, mind_ilon, mind_jlat)
  !xacc parallel
  DO itim = 1,ntim
    DO k = 1,ke+1
      !$acc loop gang vector
      DO ip = 1, ipend
        i = mind_ilon(ip,ib)
        j = mind_jlat(ip,ib)
        tke_b(ip,k,itim) = tke(i,j,k,itim)
      END DO
    END DO
  END DO
  !xacc end parallel
  !xacc end data

END SUBROUTINE block_fields_copytoblock_tke

!==============================================================================
!==============================================================================
!+ Module procedure "block_fields_copyfromblock_tke" in "src_block_fields_org"
!------------------------------------------------------------------------------

SUBROUTINE block_fields_copyfromblock_tke(ipend,ib)

  INTEGER, INTENT(IN) :: ipend,ib
  INTEGER :: i,j,k,ip, itim, ntim

  ntim=SIZE(tke,4)

  ! not before turbulence is running on GPUs
  !xacc data present(tke, tke_b, mind_ilon, mind_jlat)
  !xacc parallel
  DO itim = 1,ntim
    DO k = 1, ke+1
      !$acc loop gang vector
      DO ip = 1, ipend
        i = mind_ilon(ip,ib)
        j = mind_jlat(ip,ib)
        tke(i,j,k,itim)= tke_b(ip,k,itim)
      END DO
    END DO
  END DO
  !xacc end parallel
  !xacc end data

END SUBROUTINE block_fields_copyfromblock_tke

!==============================================================================
!==============================================================================
!+ Module procedure "block_fields_test_copy" in "src_block_fields_org" 
!------------------------------------------------------------------------------

SUBROUTINE block_fields_test_copy

!------------------------------------------------------------------------------
!
! Description:
!   Build in test to check block field interface. Copy various field
!   back and forth to block. Could be called from lmorg after block 
!   initialization
!
!------------------------------------------------------------------------------

! Locals :
! --------

  TYPE(CopylistStruct)   :: testcopyList
  INTEGER                :: ib,ipend,ip,i,j
  REAL(KIND=wp)          :: t_ref, qv_ref, t_res, qv_res
  INTEGER                :: ierror
  CHARACTER(LEN=64)      :: yerrmsg

  INTEGER, POINTER       :: ptlev
  REAL(KIND=wp), POINTER :: qv(:,:,:) 

!----------- End of header ----------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine   block_fields_test_copy
!------------------------------------------------------------------------------

  IF (my_cart_id==0) PRINT*, ' TESTING COPY TO BLOCK'
  !init
  CALL init_copy_list(testcopyList)
  ptlev=>nnew

  !get handle qv
  CALL trcr_get(ierror, idt_qv,  ptr_tlev = ptlev, ptr = qv)
  IF (ierror /= 0) THEN
    yerrmsg = trcr_errorstr(ierror)
    CALL model_abort(my_cart_id, ierror, yerrmsg, 'block_fields_test_copy')
  ENDIF

  !save checksum
  t_ref=sum(t)/size(t)
  qv_ref=sum(qv)/size(qv)

  ! register some copies
  !IN
  CALL register_copy(t_new_b  ,testcopyList, copyToBlockF)
  CALL register_copy(qv_new_b ,testcopyList, copyToBlockF)
  !OUT
  CALL register_copy(t_new_b  ,testcopyList, copyFromBlockF)
  CALL register_copy(qv_new_b ,testcopyList, copyFromBlockF)
    

  ! do block loop
  DO ib=1,nblock
    ipend=nproma
    if (ib==nblock) ipend=nlastproma

    !copy set flag and request list
    CALL request_copy(testcopyList,ierror,yerrmsg)
    IF (ierror>0) CALL model_abort(my_cart_id,ierror,yerrmsg,&
            'block_fields_test_copy:request_copy')

    !Apply copy to block
    CALL copy_to_block(testcopyList,ipend,ib,ierror,yerrmsg)
    IF (ierror>0) CALL model_abort(my_cart_id,ierror,yerrmsg,&
         'testcopy:copy_to_block')

    !Erase copyied field to make sure the copy is working          
    DO ip=1,ipend
      i = mind_ilon(ip,ib)
      j = mind_jlat(ip,ib)  
      t(i,j,:,nnew)=0._wp
      qv(i,j,:)=0._wp
    END DO

    !Copy back
    CALL copy_from_block(testcopyList,ipend,ib,ierror,yerrmsg)
    IF (ierror>0) CALL model_abort(my_cart_id,ierror,yerrmsg,&
         'testcopy:copy_to_block')

  END DO

  CALL finalize_copy(ierror, yerrmsg)
  IF (ierror>0) CALL model_abort(my_cart_id,ierror,yerrmsg,&
          'testcopy:finalize_copy_step')

  !check results
  t_res=sum(t)/size(t)
  qv_res=sum(qv)/size(qv)

  IF(abs(t_ref-t_res)+&
     abs(qv_ref-qv_res) > 10_wp*rprecision) THEN
    WRITE(*,"(I4,A)") my_cart_id, ' TEST COPY FAIL'
    WRITE(*,"(I4,ES16.4)") my_cart_id, abs(t_ref-t_res)
    WRITE(*,"(I4,ES16.4)") my_cart_id, abs(qv_ref-qv_res)
  ELSE
    IF (my_cart_id==0) PRINT*, 'TEST COPY OK'
  END IF

!------------------------------------------------------------------------------
! End of module procedure block_fields_test_copy
!------------------------------------------------------------------------------

END SUBROUTINE  block_fields_test_copy

!==============================================================================

END MODULE src_block_fields_org
