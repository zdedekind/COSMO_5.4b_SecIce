!+ Module acc_global_data for organizing accelerator memory 
!------------------------------------------------------------------------------

MODULE acc_global_data

!------------------------------------------------------------------------------
!
! Description:
!
! This module performs allocation/deallocation of global fields on the 
! accelerator using OpenACC directives. Fields are used from
!   - data_fields.f90
!   - data_lheat_nudge.f90
!   - data_soil.f90
!   - radiation_data.f90
!   - vgrid_refatm_utils
!   - grid_metrics_utilities
!   - src_stoch_physics
! 
! Tools to synchronize host and cpu memory are also provided:
!   - acc_update_device_global_data
!   - acc_update_host_global_data
!
! Method:
!
! Fortran preprocessor variables are used to define several lists 
! of variables once, which are then used in the routines in this module.
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
! V5_4         2016-03-10 Xavier Lapillonne
!  Initial release
! V5_4b        2016-07-12 Xavier Lapillonne, Ulrich Schaettler
!  Added new fields for convection schemes (XL)
!  Added new fields for COSMO-ICON versions of surface and turbulence schemes (US)
! @VERSION@    @DATE@     Ulrich Schaettler
!  Corrected some number in the ACC_LISTs
!  Eliminated a code block in SR acc_update_device_global_data, which was written
!    twice erroneously
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

USE data_fields              !all fields
USE data_parallel,           ONLY: my_cart_id
USE data_runcontrol,         ONLY: idbg_level
USE data_soil,               ONLY: csalb,cporv,cadp

#ifdef NUDGING
USE data_lheat_nudge,        ONLY:  tt_lheat, qrsflux
#endif

USE radiation_data,          ONLY:    &
    zaea, zaes, zaeg, zaef,           &
    pvdaes, pvdael, pvdaeu, pvdaed,   &
    zlwg, zlww, zlwe, zlwemn, zlwemx, &
    ziwg, ziww, ziwe, ziwemn, ziwemx, &
    coai, cobi, coali, cobti,         &
    rad_csalbw, zrsc, solant,         &
    mind_ilon_rad, mind_jlat_rad,     &
    gp_radcoarse_loc

USE vgrid_refatm_utils,      ONLY :  &
    vcoord

USE grid_metrics_utilities,  ONLY:   &
    sqrtg_r_s, sqrtg_r_u, sqrtg_r_v, sqrtg_r_w

USE acc_utilities,           ONLY: acc_copyin_1d

USE src_stoch_physics,       ONLY: pertstoph

USE src_tracer,              ONLY:   &
    trcr_acc_update_device, trcr_acc_update_host

!------------------------------------------------------------------------------

IMPLICIT NONE

!------------------------------------------------------------------------------

PRIVATE

!==============================================================================

!------------------------------------------------------------------------------
! Global (i.e. public) Declarations:
  
PUBLIC :: acc_allocate_data_fields, acc_deallocate_data_fields, &
          acc_update_device_global_data, acc_update_host_global_data

!==============================================================================
! ACC variables : list of variables which should be allocated and deallocated
!                 on the accelerator. A preprocessor #define is used in order 
!                 to avoid duplication in multiple subroutines
!------------------------------------------------------------------------------

! NOTE: Add new gpu variable here (keep grouping of variables as in data_fields)
!       if a new ACC_LIST is introduced it should be added in all routines 
!       in this file which uses ACC_LIST

! 1. constant fields for the reference atmosphere
! -----------------------------------------------

#define ACC_LIST0101  hhl         , p0         , p0hl       , dp0        , rho0



! 2. external parameter fields
! ----------------------------

#define ACC_LIST0201  hsurf       , fr_land    , plcov      , rootdp     , sso_stdh
#define ACC_LIST0202  sai         , tai        , eai        , l_pat      , rsmin2d
#define ACC_LIST0203  fr_lake     , depth_lk   , fc         , llandmask  , isoiltyp
#define ACC_LIST0204  vio3        , hmo3       , rlat       , rlon       , rmy
#define ACC_LIST0205  alb_dry     , alb_sat    , alb_dif    , for_e      , for_d
#define ACC_LIST0206  emis_rad    , aerlan     , aerurb     , aerdes     , aersea
#define ACC_LIST0207  aer_su      , aer_du     , aer_or     , aer_bc     , aer_ss
#define ACC_LIST0208  llakemask



! 3. prognostic variables
! -----------------------

#define ACC_LIST0301  pp          , t          , u          , v          , w
#define ACC_LIST0302  u_m         , v_m



! 4. tendency fields for the prognostic variables
! -----------------------------------------------

#define ACC_LIST0401  ttens       , utens      , vtens      , tketens



! 5. fields for surface values
! ----------------------------

#define ACC_LIST0501  ps          , t_s        , t_g        , qv_s



! 6. fields that are computed in the parametrization
! --------------------------------------------------

#define ACC_LIST0601  ttens_diab  , rho

! from microphysics
#define ACC_LIST0602  qrs         , prr_gsp    , prs_gsp    , prg_gsp    , tt_lheat
#define ACC_LIST0603  qrsflux


! from radiation
#define ACC_LIST0611  sohr        , sotr       , sotr_par   , thhr       , sobs
#define ACC_LIST0612  thbs        , pabs       , sobt       , thbt       , qc_rad
#define ACC_LIST0613  qi_rad      , alb_rad    , alb_rad_rc , clc_sgs    , clch
#define ACC_LIST0614  clcm        , clcl       , clct       , sodwddm    , swdir_s
#define ACC_LIST0615  swdifd_s    , swdifu_s   , swtrdir_s  , swtrdifd_s , swtrdifu_s
#define ACC_LIST0616  lwd_s       , lwu_s      , swdir_cor  , sod_t      , asod_t
#define ACC_LIST0617  sun_el      , sun_azi    , skyview    , slo_asp    , slo_ang
#define ACC_LIST0618  horizon

! from radiation_data
#define ACC_LIST0621  coai        , cobi       , coali      , cobti      , zrsc
#define ACC_LIST0622  solant      , zaea       , zaes       , zaeg       , zaef
#define ACC_LIST0623  pvdaes      , pvdael     , pvdaeu     , pvdaed     , rad_csalbw
#define ACC_LIST0624  zlwe        , zlww       , zlwg       , zlwemn     , zlwemx
#define ACC_LIST0625  ziwe        , ziww       , ziwg       , ziwemn     , ziwemx
#define ACC_LIST0626  csalb       , cporv      , cadp
#define ACC_LIST0627  mind_ilon_rad, mind_jlat_rad, gp_radcoarse_loc

!   fields from the sub-grid scale orography
#define ACC_LIST0631  tt_sso      , vt_sso     , ut_sso     , tket_sso


! from turbulence
#define ACC_LIST0641  l_hori      , gz0        , tcm        , tcm        , tfm
#define ACC_LIST0642  tfh         , tfv        , tvm        , tvh        , tkr
#define ACC_LIST0643  tkred_sfc   , rcld       , tkvm       , tkvh       , tkhm
#define ACC_LIST0644  tkhh        , hdef2      , hdiv       , dwdx       , dwdy
#define ACC_LIST0645  edr         , tket_hshr  , tket_adv   , ut_turb    , vt_turb



! from the surface schemes (TERRA, FLake, SeaIce)
#define ACC_LIST0651  t_snow      , w_snow     , h_snow     , rho_snow   , freshsnow
#define ACC_LIST0652  w_i         , w_p        , w_s        , t_so       , w_so      
#define ACC_LIST0653  w_so_ice    , t_snow_mult, w_snow_mult, wliq_snow  , rho_snow_mult
#define ACC_LIST0654 dzh_snow_mult, runoff_s   , runoff_g   , lhfl_bs    , lhfl_pl
#define ACC_LIST0655  rstom
#define ACC_LIST0656  qmomflux    , fetch_lk   , dp_bs_lk   , t_bs_lk    , gamso_lk
#define ACC_LIST0657  t_mnw_lk    , t_wml_lk   , t_bot_lk   , t_b1_lk    , c_t_lk
#define ACC_LIST0658  h_ml_lk     , h_b1_lk    , t_ice      , h_ice


! from convection
#define ACC_LIST0661  clc_con     , clw_con    , prr_con    , prs_con    , prne_con
#define ACC_LIST0662  bas_con     , top_con    , tt_conv    , ttens_conv , ttdiab_conv
#define ACC_LIST0663  qvt_conv    , qct_conv   , qit_conv   , ut_conv    , vt_conv
#define ACC_LIST0664  tket_conv   , mflx_con   , cape_con   , qcvg_con   , tke_con
#define ACC_LIST0665  vgust_con   , w_conv     , dqvdt_conv , qhfl

! fields that are computed in the dynamics and / or physics
#define ACC_LIST0671  lhfl_s      , shfl_s     , dqvdt_conv , dqvdt      , qvsflx

! 7. fields for model output and diagnostics
! ------------------------------------------

#define ACC_LIST0701  t_2m        , qv_2m      , td_2m      , rh_2m      , u_10m
#define ACC_LIST0702  v_10m

! 8. fields for the boundary values
! ---------------------------------

#define ACC_LIST0801  u_bd        , v_bd       , w_bd       , t_bd       , pp_bd
#define ACC_LIST0802  t_s_bd      , qv_s_bd    , t_snow_bd  , w_snow_bd  , hmo3_bd
#define ACC_LIST0803  vio3_bd     , w_cl_bd    , t_cl_bd    , lai_bd     , rootdp_bd
#define ACC_LIST0804  plcov_bd    , aer_du_bd  , aer_su_bd  , aer_ss_bd  , aer_bc_bd
#define ACC_LIST0805  aer_or_bd


!==============================================================================
! Module procedures in src_block_fields_org
!==============================================================================

CONTAINS 

!==============================================================================
!+ Module procedure "acc_allocate_data_fields" in "acc_global_data" 
!------------------------------------------------------------------------------

SUBROUTINE acc_allocate_data_fields()

!------------------------------------------------------------------------------
!
! Description:
!   Allocates and intializes data fields on the accelerator
!   Allocation status of each variable is the same as on the host, fields which 
!   are not allocated at this point on the host will not be allocated on the
!   accelerator.
!
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine acc_allocate_data_fields
!------------------------------------------------------------------------------

#ifdef _OPENACC
  IF ( (my_cart_id == 0 ) .AND. (idbg_level > 5) ) THEN
    WRITE(*,*) ' GPUINFO: ALLOCATING fields on accelerator' 
  ENDIF
#endif

  !Allocates and initializes variables to host values
  !ACC_LISTX are defined in the module header
  !$acc enter data           &
  !$acc copyin(ACC_LIST0101) &
  !$acc copyin(ACC_LIST0201) &

  !$acc copyin(ACC_LIST0202) &
  !$acc copyin(ACC_LIST0203) &
  !$acc copyin(ACC_LIST0204) &
  !$acc copyin(ACC_LIST0205) &
  !$acc copyin(ACC_LIST0206) &
  !$acc copyin(ACC_LIST0207) &
  !$acc copyin(ACC_LIST0208) &

  !$acc copyin(ACC_LIST0301) &
  !$acc copyin(ACC_LIST0302) &

  !$acc copyin(ACC_LIST0401) &

  !$acc copyin(ACC_LIST0501) &

  !$acc copyin(ACC_LIST0601) &
  !$acc copyin(ACC_LIST0602) &
  !$acc copyin(ACC_LIST0603) &

  !$acc copyin(ACC_LIST0611) &
  !$acc copyin(ACC_LIST0612) &
  !$acc copyin(ACC_LIST0613) &
  !$acc copyin(ACC_LIST0614) &
  !$acc copyin(ACC_LIST0615) &
  !$acc copyin(ACC_LIST0616) &
  !$acc copyin(ACC_LIST0617) &
  !$acc copyin(ACC_LIST0618) &

  !$acc copyin(ACC_LIST0621) &
  !$acc copyin(ACC_LIST0622) &
  !$acc copyin(ACC_LIST0623) &
  !$acc copyin(ACC_LIST0624) &
  !$acc copyin(ACC_LIST0625) &
  !$acc copyin(ACC_LIST0626) &
  !$acc copyin(ACC_LIST0627) &

  !$acc copyin(ACC_LIST0631) &

  !$acc copyin(ACC_LIST0641) &
  !$acc copyin(ACC_LIST0642) &
  !$acc copyin(ACC_LIST0643) &
  !$acc copyin(ACC_LIST0644) &
  !$acc copyin(ACC_LIST0645) &

  !$acc copyin(ACC_LIST0651) &
  !$acc copyin(ACC_LIST0652) &
  !$acc copyin(ACC_LIST0653) &
  !$acc copyin(ACC_LIST0654) &
  !$acc copyin(ACC_LIST0655) &
  !$acc copyin(ACC_LIST0656) &
  !$acc copyin(ACC_LIST0657) &
  !$acc copyin(ACC_LIST0658) &

  !$acc copyin(ACC_LIST0661) &
  !$acc copyin(ACC_LIST0662) &
  !$acc copyin(ACC_LIST0663) &
  !$acc copyin(ACC_LIST0664) &
  !$acc copyin(ACC_LIST0665) &

  !$acc copyin(ACC_LIST0671) &

  !$acc copyin(ACC_LIST0701) &
  !$acc copyin(ACC_LIST0702) &

  !$acc copyin(ACC_LIST0801) &
  !$acc copyin(ACC_LIST0802) &
  !$acc copyin(ACC_LIST0803) &
  !$acc copyin(ACC_LIST0804) &
  !$acc copyin(ACC_LIST0805)
   
  !XL_TODO : add in update / delete ...
  CALL acc_copyin_1d(vcoord%vert_coord)

!------------------------------------------------------------------------------
! End of module procedure acc_allocate_data_fields
!------------------------------------------------------------------------------

END SUBROUTINE acc_allocate_data_fields

!==============================================================================
!+ Module procedure "acc_deallocate_data_fields" in "acc_global_data" 
!------------------------------------------------------------------------------

SUBROUTINE acc_deallocate_data_fields()

!------------------------------------------------------------------------------
!
! Description:
!   Deallocates data fields on the accelerator. The routine should be called 
!   before deallocation of fields on the host.
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine acc_deallocate_data_fields
!------------------------------------------------------------------------------

#ifdef _OPENACC
   IF ( (my_cart_id == 0 ) .AND. (idbg_level > 5) ) THEN
     WRITE(*,*) ' GPUINFO: DEALLOCATING fields on accelerator' 
   ENDIF
#endif

!XL_TODO: currently (24.07.2014) this produces an error with CCE.
!         The GPU deallocation is for the moment deactivated
!         Waiting for cray answer
  !ACC_LISTX are defined in the module header
  !NOTUSEDacc exit data            &
  !NOTUSEDacc delete (ACC_LIST0101) &
  !NOTUSEDacc delete (ACC_LIST0201) &

  !NOTUSEDacc delete (ACC_LIST0202) &
  !NOTUSEDacc delete (ACC_LIST0203) &
  !NOTUSEDacc delete (ACC_LIST0204) &
  !NOTUSEDacc delete (ACC_LIST0205) &
  !NOTUSEDacc delete (ACC_LIST0206) &
  !NOTUSEDacc delete (ACC_LIST0207) &
  !NOTUSEDacc delete (ACC_LIST0208) &

  !NOTUSEDacc delete (ACC_LIST0301) &
  !NOTUSEDacc delete (ACC_LIST0302) &

  !NOTUSEDacc delete (ACC_LIST0401) &

  !NOTUSEDacc delete (ACC_LIST0501) &

  !NOTUSEDacc delete (ACC_LIST0601) &
  !NOTUSEDacc delete (ACC_LIST0602) &
  !NOTUSEDacc delete (ACC_LIST0603) &

  !NOTUSEDacc delete (ACC_LIST0611) &
  !NOTUSEDacc delete (ACC_LIST0612) &
  !NOTUSEDacc delete (ACC_LIST0613) &
  !NOTUSEDacc delete (ACC_LIST0614) &
  !NOTUSEDacc delete (ACC_LIST0615) &
  !NOTUSEDacc delete (ACC_LIST0616) &
  !NOTUSEDacc delete (ACC_LIST0617) &
  !NOTUSEDacc delete (ACC_LIST0618) &

  !NOTUSEDacc delete (ACC_LIST0621) &
  !NOTUSEDacc delete (ACC_LIST0622) &
  !NOTUSEDacc delete (ACC_LIST0623) &
  !NOTUSEDacc delete (ACC_LIST0624) &
  !NOTUSEDacc delete (ACC_LIST0625) &
  !NOTUSEDacc delete (ACC_LIST0626) &
  !NOTUSEDacc delete (ACC_LIST0627) &

  !NOTUSEDacc delete (ACC_LIST0631) &

  !NOTUSEDacc delete (ACC_LIST0641) &
  !NOTUSEDacc delete (ACC_LIST0642) &
  !NOTUSEDacc delete (ACC_LIST0643) &
  !NOTUSEDacc delete (ACC_LIST0644) &
  !NOTUSEDacc delete (ACC_LIST0645) &

  !NOTUSEDacc delete (ACC_LIST0651) &
  !NOTUSEDacc delete (ACC_LIST0652) &
  !NOTUSEDacc delete (ACC_LIST0653) &
  !NOTUSEDacc delete (ACC_LIST0654) &
  !NOTUSEDacc delete (ACC_LIST0655) &
  !NOTUSEDacc delete (ACC_LIST0656) &
  !NOTUSEDacc delete (ACC_LIST0657) &
  !NOTUSEDacc delete (ACC_LIST0658) &

  !NOTUSEDacc delete (ACC_LIST0661) &
  !NOTUSEDacc delete (ACC_LIST0662) &
  !NOTUSEDacc delete (ACC_LIST0663) &
  !NOTUSEDacc delete (ACC_LIST0664) &
  !NOTUSEDacc delete (ACC_LIST0665) &

  !NOTUSEDacc delete (ACC_LIST0671) &

  !NOTUSEDacc delete (ACC_LIST0701) &
  !NOTUSEDacc delete (ACC_LIST0702) &

  !NOTUSEDacc delete (ACC_LIST0801) &
  !NOTUSEDacc delete (ACC_LIST0802) &
  !NOTUSEDacc delete (ACC_LIST0803) &
  !NOTUSEDacc delete (ACC_LIST0804) &
  !NOTUSEDacc delete (ACC_LIST0805)

!------------------------------------------------------------------------------
! End of module procedure acc_deallocate_data_fields
!------------------------------------------------------------------------------

END SUBROUTINE acc_deallocate_data_fields

!==============================================================================
!+ Module procedure "acc_update_device_global_data" in "acc_global_data" 
!------------------------------------------------------------------------------

SUBROUTINE acc_update_device_global_data(mtag)

!------------------------------------------------------------------------------
!
! Description:
!   Copy all global_data and tracer variables from host to device   
!
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

!- End of header
!==============================================================================

  CHARACTER (LEN=*),    INTENT(IN)    ::   mtag      ! call-site information

!------------------------------------------------------------------------------
! Begin Subroutine acc_update_device_global_data
!------------------------------------------------------------------------------

#ifdef _OPENACC
  IF ( ( my_cart_id == 0 ) .AND. ( idbg_level >= 0 ) ) THEN
    WRITE(*,*) ' GPUINFO: CPU-GPU copy, all fields, ', TRIM(mtag)
  ENDIF
#endif

  !ACC_LISTX are defined in the module header
  !$acc update device(ACC_LIST0101)
  !$acc update device(ACC_LIST0201)

  !$acc update device(ACC_LIST0202)
  !$acc update device(ACC_LIST0203)
  !$acc update device(ACC_LIST0204)
  !$acc update device(ACC_LIST0205)
  !$acc update device(ACC_LIST0206)
  !$acc update device(ACC_LIST0207)
  !$acc update device(ACC_LIST0208)

  !$acc update device(ACC_LIST0301)
  !$acc update device(ACC_LIST0302)

  !$acc update device(ACC_LIST0401)

  !$acc update device(ACC_LIST0501)

  !$acc update device(ACC_LIST0601)
  !$acc update device(ACC_LIST0602)
  !$acc update device(ACC_LIST0603)

  !$acc update device(ACC_LIST0611)
  !$acc update device(ACC_LIST0612)
  !$acc update device(ACC_LIST0613)
  !$acc update device(ACC_LIST0614)
  !$acc update device(ACC_LIST0615)
  !$acc update device(ACC_LIST0616)
  !$acc update device(ACC_LIST0617)
  !$acc update device(ACC_LIST0618)

  !$acc update device(ACC_LIST0621)
  !$acc update device(ACC_LIST0622)
  !$acc update device(ACC_LIST0623)
  !$acc update device(ACC_LIST0624)
  !$acc update device(ACC_LIST0625)
  !$acc update device(ACC_LIST0626)
  !$acc update device(ACC_LIST0627)

  !$acc update device(ACC_LIST0631)

  !$acc update device(ACC_LIST0641)
  !$acc update device(ACC_LIST0642)
  !$acc update device(ACC_LIST0643)
  !$acc update device(ACC_LIST0644)
  !$acc update device(ACC_LIST0645)

  !$acc update device(ACC_LIST0651)
  !$acc update device(ACC_LIST0652)
  !$acc update device(ACC_LIST0653)
  !$acc update device(ACC_LIST0654)
  !$acc update device(ACC_LIST0655)
  !$acc update device(ACC_LIST0656)
  !$acc update device(ACC_LIST0657)
  !$acc update device(ACC_LIST0658)

  !$acc update device(ACC_LIST0661)
  !$acc update device(ACC_LIST0662)
  !$acc update device(ACC_LIST0663)
  !$acc update device(ACC_LIST0664)
  !$acc update device(ACC_LIST0665)

  !$acc update device(ACC_LIST0671)

  !$acc update device(ACC_LIST0701)
  !$acc update device(ACC_LIST0702)

  !$acc update device(ACC_LIST0801)
  !$acc update device(ACC_LIST0802)
  !$acc update device(ACC_LIST0803)
  !$acc update device(ACC_LIST0804)
  !$acc update device(ACC_LIST0805)

  CALL trcr_acc_update_device

!------------------------------------------------------------------------------
! End of module procedure acc_update_device_global_data
!------------------------------------------------------------------------------

END SUBROUTINE acc_update_device_global_data

!==============================================================================


!==============================================================================
!+ Module procedure "acc_update_host_global_data" in "acc_global_data" 
!------------------------------------------------------------------------------

SUBROUTINE acc_update_host_global_data(mtag)

!------------------------------------------------------------------------------
!
! Description:
!   Copy all global_data GPU variables from device to host   
!
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

!- End of header
!==============================================================================

  CHARACTER (LEN=*),    INTENT(IN)    ::   mtag      ! call-site information

!------------------------------------------------------------------------------
! Begin Subroutine acc_update_host_global_data
!------------------------------------------------------------------------------

#ifdef _OPENACC
  IF ( ( my_cart_id == 0 ) .AND. ( idbg_level >= 0 ) ) THEN
    WRITE(*,*) ' GPUINFO: GPU-CPU copy, all fields, ', TRIM(mtag)
  ENDIF
#endif

  !ACC_LISTX are defined in the module header
  !$acc update host(ACC_LIST0101)
  !$acc update host(ACC_LIST0201)

  !$acc update host(ACC_LIST0202)
  !$acc update host(ACC_LIST0203)
  !$acc update host(ACC_LIST0204)
  !$acc update host(ACC_LIST0205)
  !$acc update host(ACC_LIST0206)
  !$acc update host(ACC_LIST0207)
  !$acc update host(ACC_LIST0208)

  !$acc update host(ACC_LIST0301)
  !$acc update host(ACC_LIST0302)

  !$acc update host(ACC_LIST0401)

  !$acc update host(ACC_LIST0501)

  !$acc update host(ACC_LIST0601)
  !$acc update host(ACC_LIST0602)
  !$acc update host(ACC_LIST0603)

  !$acc update host(ACC_LIST0611)
  !$acc update host(ACC_LIST0612)
  !$acc update host(ACC_LIST0613)
  !$acc update host(ACC_LIST0614)
  !$acc update host(ACC_LIST0615)
  !$acc update host(ACC_LIST0616)
  !$acc update host(ACC_LIST0617)
  !$acc update host(ACC_LIST0618)

  !$acc update host(ACC_LIST0621)
  !$acc update host(ACC_LIST0622)
  !$acc update host(ACC_LIST0623)
  !$acc update host(ACC_LIST0624)
  !$acc update host(ACC_LIST0625)
  !$acc update host(ACC_LIST0626)
  !$acc update host(ACC_LIST0627)

  !$acc update host(ACC_LIST0631)

  !$acc update host(ACC_LIST0641)
  !$acc update host(ACC_LIST0642)
  !$acc update host(ACC_LIST0643)
  !$acc update host(ACC_LIST0644)
  !$acc update host(ACC_LIST0645)

  !$acc update host(ACC_LIST0651)
  !$acc update host(ACC_LIST0652)
  !$acc update host(ACC_LIST0653)
  !$acc update host(ACC_LIST0654)
  !$acc update host(ACC_LIST0655)
  !$acc update host(ACC_LIST0656)
  !$acc update host(ACC_LIST0657)
  !$acc update host(ACC_LIST0658)

  !$acc update host(ACC_LIST0661)
  !$acc update host(ACC_LIST0662)
  !$acc update host(ACC_LIST0663)
  !$acc update host(ACC_LIST0664)
  !$acc update host(ACC_LIST0665)

  !$acc update host(ACC_LIST0671)

  !$acc update host(ACC_LIST0701)
  !$acc update host(ACC_LIST0702)

  !$acc update host(ACC_LIST0801)
  !$acc update host(ACC_LIST0802)
  !$acc update host(ACC_LIST0803)
  !$acc update host(ACC_LIST0804)
  !$acc update host(ACC_LIST0805)

  CALL trcr_acc_update_host

!------------------------------------------------------------------------------
! End of module procedure acc_update_host_global_data
!------------------------------------------------------------------------------

END SUBROUTINE acc_update_host_global_data

!==============================================================================


END MODULE acc_global_data
