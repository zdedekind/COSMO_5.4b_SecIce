!>
!! cloud microphysics
!!
!! !------------------------------------------------------------------------------
!!
!! @par Description of *gscp_cloudice*:
!!   This module procedure calculates the rates of change of temperature, cloud
!!   water, cloud ice, water vapor, rain and snow due to cloud microphysical
!!   processes related to the formation of grid scale precipitation. This
!!   includes the sedimentation of rain and snow. The precipitation fluxes at
!!   the surface are also calculated here.
!!
!! Method:
!!   Prognostic one-moment bulk microphysical parameterization.
!!   The sedimentation of rain and snow is computed implicitly.
!!
!! Current Code Owner: DWD, A. Seifert
!!    phone: +49-69-8062-2729,  fax:   +49-69-8062-3721
!!    email: Axel.Seifert@dwd.de
!!
!! @author Guenther Doms
!! @author Axel Seifert
!!
! History:
! Version    Date       Name
! ---------- ---------- ----
! V5_1         2014-11-28 Felix Rieper
!  Modifications from Felix Rieper for modifications to improve supercooled 
!  liquid water (SLW) prediction in cloudice
!      - reduced number of ice crystal Ni(T), now according to Cooper (1986)
!      - reduced rain freezing rate srfrzr, now according to Bigg (1953)
!      - reduced depositional growth of ice and snow (zsidep, zssdep) 
!        at cloud top, according to R. Forbes (2013) -> IFS model!  
! 
!! @par Copyright
!! 2002-2009 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
!! @par reference   This is an adaption of subroutine cloudice in file src_gscp.f90
!!  of the COSMO-Model. Equation numbers refer to
!!  Doms, Foerstner, Heise, Herzog, Raschendorfer, Schrodin, Reinhardt, Vogel
!!    (September 2005): "A Description of the Nonhydrostatic Regional Model LM",
!!
!------------------------------------------------------------------------------

MODULE gscp_cloudice

!------------------------------------------------------------------------------
!>
!! Description:
!!
!!   The subroutine in this module calculates the rates of change of
!!   temperature, cloud condensate and water vapor due to cloud microphysical
!!   processes related to the formation of grid scale clouds and precipitation.
!!
!==============================================================================
!
! Declarations:
!
! Modules used:

!------------------------------------------------------------------------------
! Microphysical constants and variables
!------------------------------------------------------------------------------

#ifdef __COSMO__
USE kind_parameters, ONLY :   &
  wp ,       &    !! KIND-type parameter for real variables
  i4              !! KIND-type parameter for standard integer vaiables

USE data_constants  , ONLY :   &
!   pi,           & !!
!! 2. physical constants and related variables
!! -------------------------------------------
    t0=>t0_melt,  & !! melting temperature of ice
!   r_d,          & !! gas constant for dry air
    r_v,          & !! gas constant for water vapour
    rdv,          & !! r_d / r_v
    o_m_rdv,      & !! 1 - r_d/r_v
    rvd_m_o,      & !! r_v/r_d - 1
    cp_d,         & !! specific heat of dry air
    cpdr,         & !! 1 / cp_d
    lh_v,         & !! latent heat of vapourization
!   lh_f,         & !! latent heat of fusion
    lh_s,         & !! latent heat of sublimation
!   g,            & !! acceleration due to gravity
!   rho_w,        & !! density of liquid water (kg/m^3)
! 3. constants for parametrizations
!! ---------------------------------
    b1,           & !! variables for computing the saturation vapour pressure
    b2w,          & !! over water (w) and ice (i)
    b2i,          & !!               -- " --
    b3,           & !!               -- " --
    b4w,          & !!               -- " --
    b234w,        & !!               -- " --
    b4i,          & !!               -- " --
!> 4. tuning constants for radiation, cloud physics, turbulence
!! ------------------------------------------------------------
    qi0,          & !! cloud ice threshold for autoconversion
    qc0,          & !! cloud water threshold for autoconversion
!> 5. Precision-dependent security parameters (epsilons)
!! ------------------------------------------------------
    repsilon        !! precision of 1.0 in current floating point format

!! end of data_constants

USE data_runcontrol , ONLY :   &
    ldiabf_lh,      & ! include diabatic forcing due to latent heat in RK-scheme
    lsuper_coolw,   & ! switch for supercooled liquid water
    lsppt,          & ! switch, if .true., perturb the physical tendencies
    itype_qxpert_rn,& ! define which hum variables tend. are perturbed
    itype_qxlim_rn    ! type of reduction/removal of the perturbation 
                      ! in case of negative (qv, qc, qi) or 
                      ! supersaturated (qv) values

!------------------------------------------------------------------------------
!! COSMO environment modules
!------------------------------------------------------------------------------

USE utilities,                ONLY : message
USE meteo_utilities,          ONLY : satad       !! saturation adjustment
USE src_stoch_physics,        ONLY : apply_tqx_tend_adj
#endif
! of ifdef __COSMO__

!------------------------------------------------------------------------------

#ifdef NUDGING
USE data_lheat_nudge,           ONLY :  &
    llhn,         & ! main switch for latent heat nudging
    llhnverif,    & ! main switch for latent heat nudging
    lhn_qrs,      & ! use integrated precipitaion flux as reference
    tt_lheat,     & ! latent heat release of the model
    qrsflux         ! total precipitation flux
#endif

!------------------------------------------------------------------------------

#ifdef __ICON__
USE mo_kind,               ONLY: wp         , &
                                 i4
!USE mo_math_constants    , ONLY: pi
USE mo_physical_constants, ONLY: r_v   => rv    , & !> gas constant for water vapour
                                 r_d   => rd    , & !! gas constant for dry air
                                 rvd_m_o=>vtmpc1, & !! r_v/r_d - 1
                                 o_m_rdv        , & !! 1 - r_d/r_v
                                 rdv            , & !! r_d / r_v
                                 lh_v  => alv   , & !! latent heat of vapourization
                                 lh_s  => als   , & !! latent heat of sublimation
!                                lh_f  => alf   , & !! latent heat of fusion
                                 cp_d  => cpd   , & !! specific heat of dry air at constant press
                                 cpdr  => rcpd  , & !! (spec. heat of dry air at constant press)^-1
                                 cvdr  => rcvd  , & !! (spec. heat of dry air at const vol)^-1
                                 b3    => tmelt , & !! melting temperature of ice/snow
!                                rho_w => rhoh2o, & !! density of liquid water (kg/m^3)
                                 g     => grav  , & !! acceleration due to gravity
                                 t0    => tmelt     !! melting temperature of ice/snow

USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config

USE mo_convect_tables,     ONLY: b1    => c1es  , & !! constants for computing the sat. vapour
                                 b2w   => c3les , & !! pressure over water (l) and ice (i)
                                 b2i   => c3ies , & !!               -- " --
                                 b4w   => c4les , & !!               -- " --
                                 b4i   => c4ies , & !!               -- " --
                                 b234w => c5les     !!               -- " --
USE mo_satad,              ONLY: satad_v_3d,     &  !! new saturation adjustment
                                 sat_pres_water, &  !! saturation vapor pressure w.r.t. water
                                 sat_pres_ice!,   &  !! saturation vapor pressure w.r.t. ice
USE mo_exception,          ONLY: message, message_text
#endif

!------------------------------------------------------------------------------

! this can be used by ICON and COSMO
USE gscp_data, ONLY: &          ! all variables are used here

    ccsrim,    ccsagg,    ccsdep,    ccsvel,    ccsvxp,    ccslam,       &
    ccslxp,    ccsaxp,    ccsdxp,    ccshi1,    ccdvtp,    ccidep,       &
    ccswxp,    zconst,    zcev,      zbev,      zcevxp,    zbevxp,       &
    zvzxp,     zvz0r,                                                    &

    v0snow,    mu_rain,   rain_n0_factor,       cloud_num,               &

    x13o8,     x1o2,      x27o16,    x3o4,      x7o4,      x7o8,         &
    zbms,      zbvi,      zcac,      zccau,     zciau,     zcicri,       &
    zcrcri,    zcrfrz,    zcrfrz1,   zcrfrz2,   zeps,      zkcac,        &
    zkphi1,    zkphi2,    zkphi3,    zmi0,      zmimax,    zmsmin,       &
    zn0s0,     zn0s1,     zn0s2,     znimax_thom,          zqmin,        &
    zrho0,     zthet,     zthn,      ztmix,     ztrfrz,    zv1s,         &
    zvz0i,     x13o12,    x2o3,      x5o24,     zams,      zasmel,       &
    zbsmel,    zcsmel,                                                   &
    iautocon,  isnow_n0temp, dist_cldtop_ref,   reduce_dep_ref,          &
    tmin_iceautoconv,     zceff_fac, zceff_min,                          &
    mma, mmb

#ifdef __ICON__
! this is (at the moment) an ICON part
USE gscp_data, ONLY: &          ! all variables are used here
    vtxexp,    & !  kc_c1,     & !
    kc_c2,     & !
    kc_alpha,  & !
    kc_beta,   & !
    kc_gamma,  & !
    kc_sigma,  & !
    do_i,      & !
    co_i
#endif

!==============================================================================

IMPLICIT NONE
PRIVATE

!------------------------------------------------------------------------------
!! Public subroutines
!------------------------------------------------------------------------------

PUBLIC :: cloudice

LOGICAL, PARAMETER :: &
#ifdef __COSMO__
  lorig_icon   = .FALSE. , &  ! switch for original ICON setup (only for cloudice)
                              ! XL : should be false for COSMO ?
#else
  lorig_icon   = .TRUE.  , &  ! switch for original ICON setup (only for cloudice)
#endif

  lsedi_ice    = .TRUE.  , &  ! switch for sedimentation of cloud ice (Heymsfield & Donner 1990 *1/3)
  lstickeff    = .TRUE.       ! switch for sticking coeff. (work from Guenther Zaengl)

!------------------------------------------------------------------------------
!> Parameters and variables which are global in this module
!------------------------------------------------------------------------------

#ifdef __COSMO__
CHARACTER(132) :: message_text = ''
#endif

!==============================================================================

CONTAINS

!==============================================================================
!> Module procedure "cloudice" in "gscp_cloudice" for computing effects of
!!  grid scale precipitation including cloud water, cloud ice, rain and snow
!------------------------------------------------------------------------------

SUBROUTINE cloudice (             &
  nvec,ke,                           & !> array dimensions
  ivstart,ivend, kstart,             & !! optional start/end indicies
  idbg,                              & !! optional debug level
  zdt, dz,                           & !! numerics parameters
  t,p,rho,qv,qc,qi,qr,qs,            & !! prognostic variables
#ifdef __ICON__
  !xxx: this should become a module variable, e.g. in a new module mo_gscp_data.f90
  qi0,qc0,                           & !! cloud ice/water threshold for autoconversion
#endif
  prr_gsp,prs_gsp,                   & !! surface precipitation rates
#ifdef __COSMO__
  tinc_lh,                           & !  t-increment due to latent heat 
  pstoph,                            & !  stochastic multiplier of physics tendencies
#endif
#ifdef NUDGING
  tt_lheat,                          & !  t-increments due to latent heating (nud) 
  qrsflux,                           & !  total precipitation flux
#endif
  l_cv,                              &
  ddt_tend_t     , ddt_tend_qv     , &
  ddt_tend_qc    , ddt_tend_qi     , & !> ddt_tend_xx are tendencies
  ddt_tend_qr    , ddt_tend_qs     , & !!    necessary for dynamics
  ddt_diag_au    , ddt_diag_ac     , & !!
  ddt_diag_ev    , ddt_diag_nuc    , & !! ddt_diag_xxx are optional
  ddt_diag_idep  , ddt_diag_sdep   , & !!   diagnostic tendencies of all
  ddt_diag_agg   , ddt_diag_rim    , & !!   individual microphysical processes
  ddt_diag_rcri  , ddt_diag_icri   , & !!
  ddt_diag_dau   , ddt_diag_iau    , & !!
  ddt_diag_imelt , ddt_diag_smelt  , & !!
  ddt_diag_cfrz  , ddt_diag_rfrz   , & !!
  ddt_diag_shed                      ) !!

!------------------------------------------------------------------------------
!>
!! Description:
!!   This module procedure calculates the rates of change of temperature, cloud
!!   water, cloud ice, water vapor, rain and snow due to cloud microphysical
!!   processes related to the formation of grid scale precipitation. The
!!   variables are updated in this subroutine. Rain and snow are prognostic
!!   variables. The precipitation fluxes at the surface are stored on the
!!   corresponding global fields.
!!
!! Method:
!!
!------------------------------------------------------------------------------
!! Declarations:
!!
!------------------------------------------------------------------------------
!! Modules used: These are declared in the module declaration section
!! -------------

!! Subroutine arguments:
!! --------------------

  INTEGER, INTENT(IN) ::  &
    nvec          ,    & !> number of horizontal points
    ke                     !! number of grid points in vertical direction

  INTEGER, INTENT(IN), OPTIONAL ::  &
    ivstart   ,    & !> optional start index for horizontal direction
    ivend     ,    & !! optional end index   for horizontal direction
    kstart    ,    & !! optional start index for the vertical index
    idbg             !! optional debug level

  REAL(KIND=wp), INTENT(IN) :: &
    zdt                    !> time step for integration of microphysics     (  s  )

#ifdef __ICON__
  REAL(KIND=wp), INTENT(IN) :: &
    qi0,qc0          !> cloud ice/water threshold for autoconversion
#endif

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  ! note: that these are actually intent(in)
  !       declared as intent(inout) to avoid copying
  REAL(KIND=wp), DIMENSION(nvec,ke), INTENT(IN) ::      &
#else
  REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) ::      &   ! (ie,ke)
#endif
    dz              ,    & !> layer thickness of full levels                (  m  )
    rho             ,    & !! density of moist air                          (kg/m3)
    p                      !! pressure                                      ( Pa  )

  LOGICAL, INTENT(IN), OPTIONAL :: &
    l_cv                   !! if true, cv is used instead of cp

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  REAL(KIND=wp), DIMENSION(nvec,ke), INTENT(INOUT) ::   &
#else
  REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT) ::   &   ! dim (ie,ke)
#endif
    t               ,    & !> temperature                                   (  K  )
    qv              ,    & !! specific water vapor content                  (kg/kg)
    qc              ,    & !! specific cloud water content                  (kg/kg)
    qi              ,    & !! specific cloud ice   content                  (kg/kg)
    qr              ,    & !! specific rain content                         (kg/kg)
    qs                     !! specific snow content                         (kg/kg)

#ifdef __COSMO__
  REAL(KIND=wp), INTENT(INOUT) :: &
       tinc_lh(:,:)    ! temperature increments due to heating              ( K/s )

  REAL(KIND=wp), INTENT(IN)    :: &
       pstoph(:,:)     ! stochastic multiplier of physics tendencies
#endif  

#ifdef NUDGING
  REAL(KIND=wp), INTENT(INOUT) :: &
       tt_lheat(:,:)  ,  & !  t-increments due to latent heating (nudg) ( K/s )
       qrsflux(:,:)       ! total precipitation flux (nudg)
#endif

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  REAL(KIND=wp), DIMENSION(nvec), INTENT(INOUT) ::   &
#else
  REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) ::   &   ! dim (ie)
#endif
    prr_gsp,             & !> precipitation rate of rain, grid-scale        (kg/(m2*s))
    prs_gsp                !! precipitation rate of snow, grid-scale        (kg/(m2*s))

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  ! note: that these are actually intent(out)
  !       declared as intent(inout) to avoid copying
  REAL(KIND=wp), DIMENSION(nvec,ke), INTENT(OUT), OPTIONAL ::   &
#else
  REAL(KIND=wp), DIMENSION(:,:), INTENT(OUT), OPTIONAL ::   &     ! dim (ie,ke)
#endif
    ddt_tend_t      , & !> tendency T                                       ( 1/s )
    ddt_tend_qv     , & !! tendency qv                                      ( 1/s )
    ddt_tend_qc     , & !! tendency qc                                      ( 1/s )
    ddt_tend_qi     , & !! tendency qi                                      ( 1/s )
    ddt_tend_qr     , & !! tendency qr                                      ( 1/s )
    ddt_tend_qs         !! tendency qs                                      ( 1/s )

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  ! note: that these are actually intent(out)
  !       declared as intent(inout) to avoid copying
  REAL(KIND=wp), DIMENSION(nvec,ke), INTENT(OUT), OPTIONAL ::   &
#else
  REAL(KIND=wp), DIMENSION(:,:), INTENT(OUT), OPTIONAL ::   &   ! dim (ie,ke)
#endif
    ddt_diag_au     , & !> optional output autoconversion rate cloud to rain           ( 1/s )
    ddt_diag_ac     , & !! optional output accretion rate cloud to rain                ( 1/s )
    ddt_diag_ev     , & !! optional output evaporation of rain                         ( 1/s )
    ddt_diag_nuc    , & !! optional output mass nucleation of cloud ice                ( 1/s )
    ddt_diag_idep   , & !! optional output depositional growth of cloud ice            ( 1/s )
    ddt_diag_sdep   , & !! optional output depositional growth of snow                 ( 1/s )
    ddt_diag_agg    , & !! optional output aggregation snow collecting cloud ice       ( 1/s )
    ddt_diag_rim    , & !! optional output riming of snow by cloud water               ( 1/s )
    ddt_diag_rcri   , & !! optional output cloud ice + rain -> snow (rcri is sink qr)  ( 1/s )
    ddt_diag_icri   , & !! optional output cloud ice + rain -> snow (icri is sink qi)  ( 1/s )
    ddt_diag_dau    , & !! optional output depositional cloud ice autoconversion       ( 1/s )
    ddt_diag_iau    , & !! optional output aggregational cloud ice autoconversion      ( 1/s )
    ddt_diag_imelt  , & !! optional output melting of cloud ice                        ( 1/s )
    ddt_diag_smelt  , & !! optional output melting of snow                             ( 1/s )
    ddt_diag_cfrz   , & !! optional output freezing of cloud water                     ( 1/s )
    ddt_diag_rfrz   , & !! optional output rainwater freezing                          ( 1/s )
    ddt_diag_shed       !! optional output shedding                                    ( 1/s )


  !! Local parameters: None, parameters are in module header, gscp_data or data_constants
  !! ----------------
  
  !> Local scalars:
  !! -------------
  
  INTEGER (KIND=i4) :: &
    iv, k             !> loop indices

  REAL    (KIND=wp   ) :: nnr

  REAL    (KIND=wp   ) :: z_heat_cap_r !! reciprocal of cpdr or cvdr (depending on l_cv)

  INTEGER ::  &
    iv_start     ,    & !> start index for horizontal direction
    iv_end       ,    & !! end index for horizontal direction
    k_start      ,    & !! model level where computations start
    izdebug             !! debug level

  REAL    (KIND=wp   ) ::  &
    fpvsw, fpvsi, fqvs,& ! name of statement functions
    fxna ,             & ! statement function for ice crystal number
    fxna_cooper ,      & ! statement function for ice crystal number, Cooper(1986) 
    ztx  , zpv  , zpx ,& ! dummy arguments for statement functions
    znimax,            & ! maximum number of cloud ice crystals
    znimix,            & ! number of ice crystals at ztmix -> threshold temp for mixed-phase clouds
    zpvsw0,            & ! sat.vap. pressure at melting temperature
    zqvsw0,            & ! sat.specific humidity at melting temperature
    zdtr ,             & ! reciprocal of timestep for integration
    zscau , zscac  , zscrim , zscshe,         & ! local values of the
    zsiau , zsagg  , zsidep , zsicri, zsrcri, & ! transfer rates
    zsdau , zssdep , zssmelt, & ! defined below
    zscsum, zscmax, zcorr,  & ! terms for limiting  total cloud water depletion
    zsrsum, zsrmax,    & ! terms for limiting  total rain water depletion
    zssmax,            & ! term for limiting snow depletion
    znin,              & ! number of cloud ice crystals at nucleation
    fnuc,              & !FR: coefficient needed for Forbes (2012) SLW layer parameterization 
    znid,              & ! number of cloud ice crystals for deposition
    zmi ,              & ! mass of a cloud ice crystal
    zsvidep, zsvisub,  & ! deposition, sublimation of cloud ice
    zsimax , zsisum , zsvmax,   & ! terms for limiting total
    zqvsw,             & ! sat. specitic humidity at ice and water saturation
    zztau, zxfac, zx1, zx2, ztt, &   ! some help variables
    ztau, zphi, zhi, zdvtp, ztc, zlog_10

  REAL    (KIND=wp   ) ::  &
    zqct   ,& ! layer tendency of cloud water
    zqvt   ,& ! layer tendency of water vapour
    zqit   ,& ! layer tendency of cloud ice
    zqrt   ,& ! layer tendency of rain
    zqst      ! layer tendency of snow

  REAL    (KIND=wp   ) ::  &
    zlnqrk,zlnqsk,     & !
    zlnlogmi,qcgk_1,               & !
    qcg,tg,qvg,qrg, qsg,qig,rhog,ppg, alf,bet,m2s,m3s,hlp,maxevap,temp_c,stickeff

  LOGICAL :: &
    llqr,llqs,llqc,llqi  !   switch for existence of qr, qs, qc, qi

  REAL(KIND=wp), DIMENSION(nvec,ke) ::   &
    t_in               ,    & !> temperature                                   (  K  )
    qv_in              ,    & !! specific water vapor content                  (kg/kg)
    qc_in              ,    & !! specific cloud water content                  (kg/kg)
    qi_in              ,    & !! specific cloud ice   content                  (kg/kg)
    qr_in              ,    & !! specific rain content                         (kg/kg)
    qs_in                     !! specific snow content                         (kg/kg)


!! Local (automatic) arrays:
!! -------------------------

  REAL    (KIND=wp   ) ::  &
    zqvsi             ,     & !> sat. specitic humidity at ice and water saturation
    zvzr        (nvec),     & !
    zvzs        (nvec),     & !
    zvzi        (nvec),     & ! terminal fall velocity of ice
    zpkr        (nvec),     & !
    zpks        (nvec),     & !
    zpki        (nvec),     & ! precipitation flux of ice
    zprvr       (nvec),     & !
    zprvs       (nvec),     & !
    zprvi       (nvec),     & !
#ifdef __COSMO__
    zdummy      (nvec,8),   & !
#endif
    zcsdep            ,     & !
    zcidep            ,     & !
    zvz0s             ,     & !
    zcrim             ,     & !
    zcagg             ,     & !
    zbsdep            ,     & !
    zcslam            ,     & !
    zn0s              ,     & !
    zimr              ,     & !
    zims              ,     & !
    zimi              ,     & !
    zzar              ,     & !
    zzas              ,     & !
    zzai              ,     & !
    zqrk              ,     & !
    zqsk              ,     & !
    zqik              ,     & !
    zdtdh             ,     & !
    z1orhog           ,     & ! 1/rhog
    zrho1o2           ,     & ! (rho0/rhog)**1/2
    zeln7o8qrk        ,     & !
    zeln7o4qrk        ,     & ! FR new     
    zeln27o16qrk      ,     & !
    zeln13o8qrk       ,     & !
!    zeln3o16qrk       ,     & !
    zeln13o12qsk      ,     & !
    zeln5o24qsk       ,     & !
    zeln2o3qsk                !

  REAL    (KIND=wp   ) ::  &
    scau         , & ! transfer rate due to autoconversion of cloud water
    scac         , & ! transfer rate due to accretion of cloud water
    snuc         , & ! transfer rate due nucleation of cloud ice
    scfrz        , & ! transfer rate due homogeneous freezing of cloud water
    simelt       , & ! transfer rate due melting of cloud ice
    sidep        , & ! transfer rate due depositional growth of cloud ice
    ssdep        , & ! transfer rate due depositional growth of snow
    sdau         , & ! transfer rate due depositional cloud ice autoconversion
    srim         , & ! transfer rate due riming of snow
    sshed        , & ! transfer rate due shedding
    sicri        , & ! transfer rate due cloud ice collection by rain (sink qi)
    srcri        , & ! transfer rate due cloud ice collection by rain (sink qr)
    sagg         , & ! transfer rate due aggregation of snow and cloud ice
    siau         , & ! transfer rate due autoconversion of cloud ice
    ssmelt       , & ! transfer rate due melting of snow
    sev          , & ! transfer rate due evaporation of rain
    srfrz        , & ! transfer rate due to rainwater freezing
    reduce_dep   , & ! FR: coefficient: reduce deposition at cloud top (Forbes 2012)
    dist_cldtop(nvec)! FR: distance from cloud top layer

  LOGICAL :: ldum

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
!! Begin Subroutine cloudice
!------------------------------------------------------------------------------

!> Statement functions
! -------------------

! saturation vapour pressure over water (fpvsw), over ice (fpvsi)
! and specific humidity at vapour saturation (fqvs)
  fpvsw(ztx)     = b1*EXP( b2w*(ztx-b3)/(ztx-b4w) )
  fpvsi(ztx)     = b1*EXP( b2i*(ztx-b3)/(ztx-b4i) )
  fqvs (zpv,zpx) = rdv*zpv/( zpx - o_m_rdv*zpv )
  
! Number of activate ice crystals;  ztx is temperature
  fxna(ztx)   = 1.0E2_wp * EXP(0.2_wp * (t0 - ztx))
  fxna_cooper(ztx) = 5.0E+0_wp * EXP(0.304_wp * (t0 - ztx))   ! FR: Cooper (1986) used by Greg Thompson(2008)

! Define reciprocal of heat capacity of dry air (at constant pressure vs at constant volume)

#ifdef __COSMO__
  z_heat_cap_r = cpdr
#endif

#ifdef __ICON__
  IF (PRESENT(l_cv)) THEN
    IF (l_cv) THEN
      z_heat_cap_r = cvdr
    ELSE
      z_heat_cap_r = cpdr
    ENDIF
  ELSE
    z_heat_cap_r = cpdr
  ENDIF
#endif

!------------------------------------------------------------------------------
!  Section 1: Initial setting of local and global variables
!------------------------------------------------------------------------------

! Some constant coefficients
  IF( lsuper_coolw) THEN
    znimax = znimax_Thom         !znimax_Thom = 250.E+3_wp,
    znimix = fxna_cooper(ztmix) ! number of ice crystals at temp threshold for mixed-phase clouds
  ELSEIF(lorig_icon) THEN
    znimax = 150.E+3_wp     ! from previous ICON code 
    znimix = fxna_cooper(ztmix) ! number of ice crystals at temp threshold for mixed-phase clouds
  ELSE
    znimax = fxna(zthn) ! Maximum number of cloud ice crystals
    znimix = fxna(ztmix) ! number of ice crystals at temp threshold for mixed-phase clouds
  END IF

  zpvsw0 = fpvsw(t0)  ! sat. vap. pressure for t = t0
  zlog_10 = LOG(10._wp) ! logarithm of 10

! Optional arguments

  IF (PRESENT(ddt_tend_t)) THEN
    ! save input arrays for final tendency calculation
    t_in  = t
    qv_in = qv
    qc_in = qc
    qi_in = qi
    qr_in = qr
    qs_in = qs
  END IF
  IF (PRESENT(ivstart)) THEN
    iv_start = ivstart
  ELSE
    iv_start = 1
  END IF
  IF (PRESENT(ivend)) THEN
    iv_end = ivend
  ELSE
    iv_end = nvec
  END IF
  IF (PRESENT(kstart)) THEN
    k_start = kstart
  ELSE
    k_start = 1
  END IF
  IF (PRESENT(idbg)) THEN
    izdebug = idbg
  ELSE
    izdebug = 0
  END IF

! timestep for calculations
  zdtr  = 1.0_wp / zdt

! output for various debug levels
  IF (izdebug > 15) CALL message('gscp_cloudice: ','Start of cloudice')
  IF (izdebug > 20) THEN
    WRITE (message_text,*) '   nvec = ',nvec       ; CALL message('',message_text)
    WRITE (message_text,*) '   ke = ',ke           ; CALL message('',message_text)
    WRITE (message_text,*) '   ivstart = ',ivstart ; CALL message('',message_text)
    WRITE (message_text,*) '   ivend   = ',ivend   ; CALL message('',message_text)
  END IF
  IF (izdebug > 50) THEN
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN dz  = ',MAXVAL(dz),MINVAL(dz)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN T   = ',MAXVAL(t),MINVAL(t)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN p   = ',MAXVAL(p),MINVAL(p)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN rho = ',MAXVAL(rho),MINVAL(rho)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qv  = ',MAXVAL(qv),MINVAL(qv)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qc  = ',MAXVAL(qc),MINVAL(qc)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qr  = ',MAXVAL(qr),MINVAL(qr)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qi  = ',MAXVAL(qi),MINVAL(qi)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qs  = ',MAXVAL(qs),MINVAL(qs)
    CALL message('',message_text)
  ENDIF

  DO iv = iv_start, iv_end

    ! Delete precipitation fluxes from previous timestep
    prr_gsp (iv) = 0.0_wp
    prs_gsp (iv) = 0.0_wp
    zpkr(iv)     = 0.0_wp
    zpks(iv)     = 0.0_wp
    zpki(iv)     = 0.0_wp
    zprvr(iv)    = 0.0_wp
    zprvs(iv)    = 0.0_wp
    zprvi(iv)    = 0.0_wp
    zvzr(iv)     = 0.0_wp
    zvzs(iv)     = 0.0_wp
    zvzi(iv)     = 0.0_wp
    dist_cldtop(iv) = 0.0_wp

  END DO

! *********************************************************************
! Loop from the top of the model domain to the surface to calculate the
! transfer rates  and sedimentation terms
! *********************************************************************

  loop_over_levels: DO  k = k_start, ke

    DO iv = iv_start, iv_end !loop over horizontal domain

#ifdef NUDGING
      ! add part of latent heating calculated in subroutine cloudice to model latent
      ! heating field: subtract temperature from model latent heating field
      IF (llhn .OR. llhnverif) THEN
        IF (lhn_qrs) THEN
          qrsflux(iv,k) = 0.0_wp
        ENDIF
        tt_lheat(iv,k) = tt_lheat(iv,k) - t(iv,k)
        ! replaces: CALL get_gs_lheating ('add',1,ke)
      ENDIF
#endif

#ifdef __COSMO__
      IF ( ldiabf_lh ) THEN
        ! initialize temperature increment due to latent heat
        tinc_lh(iv,k) = tinc_lh(iv,k) - t(iv,k)
      ENDIF
#endif

      !----------------------------------------------------------------------------
      ! Section 2: Check for existence of rain and snow
      !            Initialize microphysics and sedimentation scheme
      !----------------------------------------------------------------------------

      zcrim  = 0.0_wp
      zcagg  = 0.0_wp
      zbsdep = 0.0_wp
      zvz0s  = 0.0_wp
      zn0s   = zn0s0
      reduce_dep = 1.0_wp  !FR: Reduction coeff. for dep. growth of rain and ice  

      qvg  = qv(iv,k)
      qcg  = qc(iv,k)
      qrg  = qr(iv,k)
      qsg  = qs(iv,k)
      qig  = qi(iv,k)
      tg   = t(iv,k)
      ppg  = p(iv,k)
      rhog = rho(iv,k)

      !..for density correction of fall speeds
      z1orhog = 1.0_wp/rhog
      zrho1o2 = EXP(LOG(zrho0*z1orhog)*x1o2)

      zqrk = qrg * rhog
      zqsk = qsg * rhog
      zqik = qig * rhog

      llqr = zqrk > zqmin
      llqs = zqsk > zqmin
      llqi = zqik > zqmin

      zdtdh = 0.5_wp * zdt / dz(iv,k)

      zzar   = zqrk/zdtdh + zprvr(iv) + zpkr(iv)
      zzas   = zqsk/zdtdh + zprvs(iv) + zpks(iv)
      zzai   = zqik/zdtdh + zprvi(iv) + zpki(iv)

      ! qs_prepare:
      !------------
      IF (llqs) THEN
        IF (isnow_n0temp == 1) THEN
          ! Calculate n0s using the temperature-dependent
          ! formula of Field et al. (2005)
          ztc = tg - t0
          ztc = MAX(MIN(ztc,0.0_wp),-40.0_wp)
          zn0s = zn0s1*EXP(zn0s2*ztc)
          zn0s = MIN(zn0s,1.0E9_wp)
          zn0s = MAX(zn0s,1.0E6_wp)
        ELSEIF (isnow_n0temp == 2) THEN
          ! Calculate n0s using the temperature-dependent moment
          ! relations of Field et al. (2005)
          ztc = tg - t0
          ztc = MAX(MIN(ztc,0.0_wp),-40.0_wp)

          nnr  = 3.0_wp
          hlp = mma(1) + mma(2)*ztc + mma(3)*nnr + mma(4)*ztc*nnr &
              + mma(5)*ztc**2 + mma(6)*nnr**2 + mma(7)*ztc**2*nnr &
              + mma(8)*ztc*nnr**2 + mma(9)*ztc**3 + mma(10)*nnr**3
          alf = EXP(hlp*zlog_10) ! 10.0_wp**hlp
          bet = mmb(1) + mmb(2)*ztc + mmb(3)*nnr + mmb(4)*ztc*nnr &
              + mmb(5)*ztc**2 + mmb(6)*nnr**2 + mmb(7)*ztc**2*nnr &
              + mmb(8)*ztc*nnr**2 + mmb(9)*ztc**3 + mmb(10)*nnr**3

          ! Here is the exponent bms=2.0 hardwired! not ideal! (Uli Blahak)
          m2s = qsg * rhog / zams   ! UB rho added as bugfix
          m3s = alf*EXP(bet*LOG(m2s))

          hlp  = zn0s1*EXP(zn0s2*ztc)
          zn0s = 13.50_wp * m2s * (m2s / m3s) **3
          zn0s = MAX(zn0s,0.5_wp*hlp)
          zn0s = MIN(zn0s,1.0E2_wp*hlp)
          zn0s = MIN(zn0s,1.0E9_wp)
          zn0s = MAX(zn0s,1.0E6_wp)
        ELSE
          ! Old constant n0s
          zn0s = 8.0E5_wp
        ENDIF
        zcrim  = ccsrim*zn0s
        zcagg  = ccsagg*zn0s
        zbsdep = ccsdep*SQRT(v0snow)
        zvz0s  = ccsvel*EXP(ccsvxp * LOG(zn0s))

        !      IF (zvzs(iv) == 0.0_wp) THEN
        zvzs(iv) = zvz0s * EXP (ccswxp * LOG (0.5_wp*zqsk)) * zrho1o2
        !      ENDIF
      ENDIF

      ! qr_sedi:
      !---------
      IF (llqr) THEN
        !      IF (zvzr(iv) == 0.0_wp) THEN
        !replaced: zvzr(iv) = zvz0r * EXP (x1o8  * LOG (0.5_wp*zqrk)) * zrho1o2
        zvzr(iv) = zvz0r * EXP (zvzxp  * LOG (0.5_wp*zqrk)) * zrho1o2
        !      ENDIF
      ENDIF

      ! qi_sedi:
      !---------
      IF (llqi) THEN
        !      IF (zvzi(iv) == 0.0_wp .AND. lnew_ice_sedi) THEN
        !! density correction not needed zrho1o2
        zvzi(iv) = zvz0i * EXP (zbvi  * LOG (0.5_wp*zqik))
        !      ENDIF
      ENDIF

      !----------------------------------------------------------------------------
      ! Section 3:
      !----------------------------------------------------------------------------

!     zeln7o8qrk    = 0.0_wp
!     zeln7o4qrk    = 0.0_wp
!     zeln27o16qrk  = 0.0_wp
!     zeln13o8qrk   = 0.0_wp
!     zeln13o12qsk  = 0.0_wp
!     zeln5o24qsk   = 0.0_wp
!     zeln2o3qsk    = 0.0_wp

!     zcsdep        = 3.367E-2_wp    
!     zcidep        = 1.3E-5_wp
!     zcslam        = 1.0E10_wp

!     scau          = 0.0_wp
!     scac          = 0.0_wp
!     snuc          = 0.0_wp
!     scfrz         = 0.0_wp
!     simelt        = 0.0_wp
!     sidep         = 0.0_wp
!     ssdep         = 0.0_wp
!     sdau          = 0.0_wp
!     srim          = 0.0_wp
!     sshed         = 0.0_wp
!     sicri         = 0.0_wp
!     srcri         = 0.0_wp
!     sagg          = 0.0_wp
!     siau          = 0.0_wp
!     ssmelt        = 0.0_wp
!     sev           = 0.0_wp
!     srfrz         = 0.0_wp

#ifdef __COSMO__
      zqvsi = fqvs( fpvsi(tg), ppg )
#endif

#ifdef __ICON__
      zqvsi = sat_pres_ice(tg)/(rhog * r_v * tg)
#endif

      llqr = zqrk > zqmin
      llqs = zqsk > zqmin
      llqi = zqik > zqmin

      IF (llqr) THEN
        !US reported by Thorsten Reinhardt: Multiplication with zrho1o2 was missing
        zpkr(iv) = zqrk * zvz0r * EXP (zvzxp * LOG (zqrk)) * zrho1o2
      ELSE
        zpkr(iv) = 0.0_wp
      ENDIF

      IF (llqs) THEN
        !US reported by Thorsten Reinhardt: Multiplication with zrho1o2 was missing
        zpks(iv) = zqsk * zvz0s * EXP (ccswxp * LOG (zqsk)) * zrho1o2
      ELSE
        zpks(iv) = 0.0_wp
      ENDIF

      IF (llqi) THEN
        zpki(iv) = zqik * zvz0i * EXP (zbvi * LOG (zqik)) * zrho1o2
      ELSE
        zpki(iv) = 0.0_wp
      ENDIF

      zpkr(iv)   = MIN( zpkr(iv) , zzar )
      zpks(iv)   = MIN( zpks(iv) , zzas )
      zpki(iv)   = MIN( zpki(iv) , zzai )

      zzar = zdtdh * (zzar-zpkr(iv))
      zzas = zdtdh * (zzas-zpks(iv))
      zzai = zdtdh * (zzai-zpki(iv))

      zimr = 1.0_wp / (1.0_wp + zvzr(iv) * zdtdh)
      zims = 1.0_wp / (1.0_wp + zvzs(iv) * zdtdh)
      zimi = 1.0_wp / (1.0_wp + zvzi(iv) * zdtdh)

      zqrk = zzar*zimr
      zqsk = zzas*zims
      zqik = zzai*zimi

      llqr = zqrk > zqmin
      llqs = zqsk > zqmin
      llqi = zqik > zqmin
      llqc =  qcg > zqmin

      ! qr:
      !----
      IF (llqr) THEN
        llqi =  qig > zqmin

        zlnqrk       = LOG (zqrk)
        IF ( qig+qcg > zqmin ) THEN
          zeln7o8qrk   = EXP (x7o8   * zlnqrk)
        ELSE
          zeln7o8qrk   = 0.0_wp
        ENDIF
        IF ( tg < ztrfrz ) THEN
          zeln7o4qrk   = EXP (x7o4   * zlnqrk)
          zeln27o16qrk = EXP (x27o16 * zlnqrk)
        ELSE
          zeln7o4qrk   = 0.0_wp
          zeln27o16qrk = 0.0_wp
        ENDIF
        IF (llqi) THEN
          zeln13o8qrk  = EXP (x13o8  * zlnqrk)
        ELSE
          zeln13o8qrk  = 0.0_wp
        ENDIF
        !      IF (qcg <= 0.0_wp ) THEN
        !        zeln3o16qrk  = EXP (x3o16  * zlnqrk)
        !      ELSE
        !        zeln3o16qrk  = 0.0_wp
        !      ENDIF
      ELSE
        zeln7o8qrk   = 0.0_wp
        zeln7o4qrk   = 0.0_wp
        zeln27o16qrk = 0.0_wp
        zeln13o8qrk  = 0.0_wp
        ! zeln3o16qrk  = 0.0_ireals
      ENDIF

      ! qs:
      !----
      IF (llqs) THEN
        zlnqsk       = LOG (zqsk)
        IF (qig+qcg > zqmin) THEN
          zeln13o12qsk = EXP (x13o12 * zlnqsk)
        ELSE
          zeln13o12qsk = 0.0_wp
        ENDIF
        zeln5o24qsk  = EXP (x5o24  * zlnqsk)
        zeln2o3qsk   = EXP (x2o3   * zlnqsk)
      ELSE
        zeln13o12qsk = 0.0_wp
        zeln5o24qsk  = 0.0_wp
        zeln2o3qsk   = 0.0_wp
      ENDIF

      ! qi_qs:
      !-------
      IF (llqi .OR. llqs) THEN
        zdvtp  = ccdvtp * EXP(1.94_wp * LOG(tg)) / ppg          
        zhi    = ccshi1*zdvtp*rhog*zqvsi/(tg*tg)
        hlp    = zdvtp / (1.0_wp + zhi)
        zcidep = ccidep * hlp
        IF (llqs) THEN
          zcslam = EXP(ccslxp * LOG(ccslam * zn0s / zqsk ))
          zcslam = MIN(zcslam,1.0E15_wp)
          zcsdep = 4.0_wp * zn0s * hlp
        ELSE
          zcslam        = 1.0E10_wp
          zcsdep        = 3.367E-2_wp
        ENDIF
      ELSE
        zcidep        = 1.3E-5_wp
        zcslam        = 1.0E10_wp
        zcsdep        = 3.367E-2_wp
      ENDIF

      !--------------------------------------------------------------------------
      ! Section 4: Initialize the conversion rates with zeros in every layer
      !            Deposition nucleation for low temperatures below a threshold
      !--------------------------------------------------------------------------

      scau          = 0.0_wp
      scac          = 0.0_wp
      snuc          = 0.0_wp
      scfrz         = 0.0_wp
      simelt        = 0.0_wp
      sidep         = 0.0_wp
      ssdep         = 0.0_wp
      sdau          = 0.0_wp
      srim          = 0.0_wp
      sshed         = 0.0_wp
      sicri         = 0.0_wp
      srcri         = 0.0_wp
      sagg          = 0.0_wp
      siau          = 0.0_wp
      ssmelt        = 0.0_wp
      sev           = 0.0_wp
      srfrz         = 0.0_wp

      ! icenucleation:
      !---------------
      IF ( tg < zthet .AND. qvg >  8.E-6_wp &    !XL_CHECK <- 8.E-6 maybe an issue with single precsion
                      .AND. qig <= 0.0_wp)   THEN
        IF( qvg > zqvsi ) THEN
          IF( lsuper_coolw .OR. lorig_icon) THEN
            znin  = MIN( fxna_cooper(tg), znimax )
          ELSE
            znin  = MIN( fxna(tg), znimax )
          ENDIF
          snuc = zmi0 * z1orhog * znin * zdtr
        ENDIF
      ENDIF

      !--------------------------------------------------------------------------
      ! Section 5: Search for cloudy grid points with cloud water and
      !            calculation of the conversion rates involving qc
      !--------------------------------------------------------------------------

      ! qc:
      !----
      IF ( llqc ) THEN
        IF (iautocon == 0) THEN
          ! Kessler (1969) autoconversion rate
          zscau  = zccau * MAX( qcg - qc0, 0.0_wp )
          zscac  = zcac  * qcg * zeln7o8qrk
        ELSEIF (iautocon == 1) THEN
          ! Seifert and Beheng (2001) autoconversion rate
          ! with constant cloud droplet number concentration cloud_num
          IF (qcg > 1.0E-6_wp) THEN
            ztau   = MIN(1.0_wp-qcg/(qcg+qrg),0.9_wp)
            ztau   = MAX(ztau,1.0E-30_wp)
            hlp    = EXP(zkphi2*LOG(ztau))
            zphi   = zkphi1 * hlp * (1.0_wp - hlp)**3
            zscau  = zconst * qcg*qcg*qcg*qcg &
                   * (1.0_wp + zphi/(1.0_wp - ztau)**2)
            zphi   = (ztau/(ztau+zkphi3))**4
            zscac  = zkcac * qcg * qrg * zphi !* zrho1o2
          ELSE
            zscau  = 0.0_wp
            zscac  = 0.0_wp
          ENDIF
        ENDIF

        IF (llqs) THEN
          zscrim = zcrim * EXP(ccsaxp * LOG(zcslam)) * qcg !* zrho1o2
        ELSE
          zscrim = 0.0_wp
        ENDIF

        zscshe = 0.0_wp
        IF( tg >= t0 ) THEN
          zscshe = zscrim
          zscrim = 0.0_wp
        ENDIF
        ! Check for maximum depletion of cloud water and adjust the
        ! transfer rates accordingly
        zscmax = qcg*zdtr 
        zscsum = zscau + zscac + zscrim + zscshe 
        zcorr  = zscmax / MAX( zscmax, zscsum )
        IF( tg <= zthn ) THEN
          scfrz = zscmax
        ELSE
          scau  = zcorr*zscau
          scac  = zcorr*zscac
          srim  = zcorr*zscrim
          sshed = zcorr*zscshe 
        ENDIF

        ! Calculation of heterogeneous nucleation of cloud ice.
        ! This is done in this section, because we require water saturation
        ! for this process (i.e. the existence of cloud water) to exist.
        ! Hetrogeneous nucleation is assumed to occur only when no
        ! cloud ice is present and the temperature is below a nucleation
        ! threshold.
        IF( (tg <= 267.15_wp) .AND. (qig <= 0.0_wp)) THEN
          IF  (lsuper_coolw .OR. lorig_icon) THEN
            znin = MIN( fxna_cooper(tg), znimax )
            snuc = zmi0 * z1orhog * znin * zdtr
          ELSE 
            znin = MIN( fxna(tg), znimax )
            snuc = zmi0 * z1orhog * znin * zdtr
          END IF
        ENDIF
        ! Calculation of in-cloud rainwater freezing
        IF (tg < ztrfrz)  THEN
          IF (lsuper_coolw) THEN
            srfrz = zcrfrz1*(EXP(zcrfrz2*(ztrfrz-tg))-1.0_wp ) * zeln7o4qrk
          ELSE 
            srfrz = zcrfrz*SQRT( (ztrfrz-tg)**3 )* zeln27o16qrk
          ENDIF
        ENDIF

        ! Calculation of reduction of depositional growth at cloud top (Forbes 2012)
        IF( k>1 .AND. k<ke .AND. lsuper_coolw ) THEN
          znin = MIN(fxna_cooper(tg), znimax )
          fnuc = MIN(znin/znimix, 1.0_wp)

          qcgk_1 = qc(iv,k-1)

          !! distance from cloud top
          IF( qcgk_1 < zqmin ) THEN     ! upper cloud layer
            dist_cldtop(iv) = 0.0_wp    ! reset distance to upper cloud layer
          ELSE
            dist_cldtop(iv) = dist_cldtop(iv) + dz(iv,k)
          END IF

          ! with asymptotic behaviour dz -> 0 (xxx)
          !        reduce_dep = MIN(fnuc + (1.0_wp-fnuc)*(reduce_dep_ref + &
          !                             dist_cldtop(iv)/dist_cldtop_ref + &
          !                             (1.0_wp-reduce_dep_ref)*(zdh/dist_cldtop_ref)**4), 1.0_wp)

          ! without asymptotic behaviour dz -> 0
          reduce_dep = MIN(fnuc + (1.0_wp-fnuc)*(reduce_dep_ref + &
               dist_cldtop(iv)/dist_cldtop_ref), 1.0_wp)

        ENDIF ! Reduction of dep. growth of snow/ice 
      ENDIF

      !------------------------------------------------------------------------
      ! Section 6: Search for cold grid points with cloud ice and/or snow and
      !            calculation of the conversion rates involving qi and ps
      !------------------------------------------------------------------------

      ! qi_qs:
      !-------
      IF ((zqik > zqmin) .OR. llqs) THEN
        ! careful: lqi is qig>zqmin, and not zqik > zqmin here: therefore cannot replace
        ! and set llqi again (above it was only set for llqr)
        llqi =  qig > zqmin

        IF (tg<=t0) THEN
          IF( lsuper_coolw .OR. lorig_icon) THEN
            znin   = MIN( fxna_cooper(tg), znimax )
          ELSE
            znin   = MIN( fxna(tg), znimax )
          END IF
          IF (lstickeff .OR. lorig_icon) THEN
            stickeff = MIN(EXP(0.09_wp*(tg-t0)),1.0_wp)
            stickeff = MAX(stickeff, zceff_min, zceff_fac*(tg-tmin_iceautoconv))
          ELSE !original sticking efficiency of cloud ice
            stickeff = MIN(EXP(0.09_wp*(tg-t0)),1.0_wp)
            stickeff = MAX(stickeff,0.2_wp)
          END IF
          zmi      = MIN( rhog*qig/znin, zmimax )
          zmi      = MAX( zmi0, zmi )
          zsvmax   = (qvg - zqvsi) * zdtr
          zsagg    = zcagg * EXP(ccsaxp*LOG(zcslam)) * qig
          zsagg    = MAX( zsagg, 0.0_wp ) * stickeff
          znid     = rhog * qig/zmi
          IF (llqi) THEN
            zlnlogmi = LOG (zmi)
            zsidep   = zcidep * znid * EXP(0.33_wp * zlnlogmi) * (qvg - zqvsi)
          ELSE
            zsidep = 0.0_wp
          ENDIF
          zsvidep   = 0.0_wp
          zsvisub   = 0.0_wp
          zsimax    = qig*zdtr 
          IF( zsidep > 0.0_wp ) THEN
            IF (lsuper_coolw ) THEN
              zsidep = zsidep * reduce_dep  !FR new: SLW reduction
            END IF
            zsvidep = MIN( zsidep, zsvmax )
          ELSEIF (zsidep < 0.0_wp ) THEN
            zsvisub = - MAX(-zsimax, zsvmax )
          ENDIF
          zsiau = zciau * MAX( qig - qi0, 0.0_wp ) * stickeff
          IF (llqi) THEN
            zlnlogmi = LOG(zmsmin/zmi)
            zztau    = 1.5_wp*( EXP(0.66_wp*zlnlogmi) - 1.0_wp)
            zsdau    = zsvidep/MAX(zztau,zeps)
          ELSE
            zsdau    =  0.0_wp
          ENDIF
          zsicri    = zcicri * qig * zeln7o8qrk
          zsrcri    = zcrcri * (qig/zmi) * zeln13o8qrk
          zxfac     = 1.0_wp + zbsdep * EXP(ccsdxp*LOG(zcslam))
          zssdep    = zcsdep * zxfac * ( qvg - zqvsi ) / (zcslam+zeps)**2

          ! Check for maximal depletion of vapor by sdep
          IF (zssdep > 0.0_wp) THEN
            IF (lsuper_coolw) THEN
              zssdep = zssdep*reduce_dep  !FR new: SLW reduction
            END IF
            zssdep = MIN(zssdep, zsvmax-zsvidep)
          END IF

          ! Check for maximal depletion of snow by sdep
          IF (zssdep < 0.0_wp) zssdep = MAX(zssdep, -qsg*zdtr)

          zsisum = zsiau + zsdau + zsagg + zsicri + zsvisub
          zcorr  = 0.0_wp
          IF( zsimax > 0.0_wp ) zcorr  = zsimax / MAX( zsimax, zsisum )
          sidep  = zsvidep - zcorr*zsvisub
          sdau   = zcorr*zsdau
          siau   = zcorr*zsiau
          sagg   = zcorr*zsagg
          ssdep  = zssdep
          srcri  = zsrcri
          sicri  = zcorr*zsicri

        ELSE ! tg > 0

          !------------------------------------------------------------------------
          ! Section 7: Search for warm grid points with cloud ice and/or snow and
          !            calculation of the melting rates of qi and qs
          !------------------------------------------------------------------------

          simelt = qig*zdtr
          zqvsw0 = zpvsw0 / (rhog * r_v * tg)
          zx1    = (tg - t0) + zasmel*(qvg - zqvsw0)
          zx2    = 1.0_wp + zbsmel * zeln5o24qsk
          zssmelt= zcsmel * zx1 * zx2 * zeln2o3qsk
          ssmelt = MAX( zssmelt, 0.0_wp )
        ENDIF ! tg

      ENDIF ! qi_qs

      !--------------------------------------------------------------------------
      ! Section 8: Search for grid points with rain in subsaturated areas
      !            and calculation of the evaporation rate of rain
      !--------------------------------------------------------------------------

      ! qr_nocloud:
      !------------
      IF (llqr .AND. qcg <= 0.0_wp) THEN

#ifdef __COSMO__
        zqvsw    = fqvs( fpvsw(tg), ppg )
#endif

#ifdef __ICON__
        zqvsw    = sat_pres_water(tg)/(rhog * r_v *tg)
#endif

        zlnqrk      = LOG (zqrk)
        zx1         = 1.0_wp + zbev * EXP (zbevxp  * zlnqrk)
        ! Limit evaporation rate in order to avoid overshoots towards supersaturation
        ! the pre-factor approximates (esat(T_wb)-e)/(esat(T)-e) at temperatures between 0 degC and 30 degC
        temp_c  = tg - t0
        maxevap = (0.61_wp-0.0163_wp*temp_c+1.111e-4_wp*temp_c**2)*(zqvsw-qvg)/zdt
        sev     = MIN(zcev*zx1*(zqvsw - qvg) * EXP (zcevxp  * zlnqrk), maxevap)
        !      sev    = zcev*zx1*(zqvsw - qvg) * EXP (zcevxp  * zlnqrk)

        !      zqvsw  = fqvs( fpvsw(tg), ppg )
        !      zx1    = 1.0_wp + zbev* zeln3o16qrk
        !      zsev   = zcev*zx1*(zqvsw - qvg)*SQRT(zqrk)
        !      sev    = MAX( zsev, 0.0_wp )

        ! Calculation of below-cloud rainwater freezing
        IF ( tg < ztrfrz ) THEN
          IF ( lsuper_coolw ) THEN
            !FR new: reduced rain freezing rate
            srfrz = zcrfrz1*(EXP(zcrfrz2*(ztrfrz-tg))-1.0_wp ) * zeln7o4qrk
          ELSE
            srfrz = zcrfrz*SQRT( (ztrfrz-tg)**3 ) * zeln27o16qrk
          ENDIF
        ENDIF

      ENDIF ! qr_nocloud

      !--------------------------------------------------------------------------
      ! Section 9: Calculate the total tendencies of the prognostic variables.
      !            Update the prognostic variables in the interior domain.
      !--------------------------------------------------------------------------

      zsrmax = zzar*z1orhog*zdtr
      zssmax = zzas*z1orhog*zdtr
      zsrsum = sev + srfrz + srcri
      zcorr  = 1.0_wp
      IF(zsrsum > 0) THEN
        zcorr  = zsrmax / MAX( zsrmax, zsrsum )
      ENDIF
      sev    = zcorr*sev
      srfrz  = zcorr*srfrz
      srcri  = zcorr*srcri
      ssmelt = MIN(ssmelt, zssmax)
      IF (ssdep < 0.0_wp ) THEN
        ssdep = MAX(ssdep, - zssmax)
      ENDIF
      zqvt = sev   - sidep - ssdep  - snuc
      zqct = simelt- scau  - scfrz  - scac   - sshed - srim 
      zqit = snuc  + scfrz - simelt - sicri  + sidep - sdau  - sagg  - siau
      zqrt = scau  + sshed + scac   + ssmelt - sev   - srcri - srfrz 
      zqst = siau  + sdau  + sagg   - ssmelt + sicri + srcri + srim  + ssdep + srfrz
      ztt = z_heat_cap_r*( lh_v*(zqct+zqrt) + lh_s*(zqit+zqst) )

#ifdef __COSMO__
      IF (lsppt) THEN
!US how to check?          IF(ntstep>0.AND.i>=istart.and.i<=iend.and.j>=jstart.and.j<=jend) THEN
          CALL apply_tqx_tend_adj(itype_qxpert_rn, itype_qxlim_rn, ppg,    &
                                  tg,  qvg,  qcg,  qig,  qrg,  qsg,        &
                  pstoph(iv,k),  ztt, zqvt, zqct, zqit, zqrt, zqst, ldum)
!US                        ENDIF
      ENDIF
#endif

      ! Update variables and add qi to qrs for water loading 
      IF (lsedi_ice .OR. lorig_icon) THEN
        qig = MAX ( 0.0_wp, (zzai*z1orhog + zqit*zdt)*zimi)
      ELSE
        qig = MAX ( 0.0_wp, qig + zqit*zdt)
      END IF
      qrg = MAX ( 0.0_wp, (zzar*z1orhog + zqrt*zdt)*zimr)
      qsg = MAX ( 0.0_wp, (zzas*z1orhog + zqst*zdt)*zims)

      !----------------------------------------------------------------------
      ! Section 10: Complete time step
      !----------------------------------------------------------------------

      IF ( k /= ke) THEN
        ! Store precipitation fluxes and sedimentation velocities 
        ! for the next level
        zprvr(iv) = qrg*rhog*zvzr(iv)
        zprvs(iv) = qsg*rhog*zvzs(iv)
        zprvi(iv) = qig*rhog*zvzi(iv)

        IF (zprvr(iv) <= zqmin) zprvr(iv)=0.0_wp
        IF (zprvs(iv) <= zqmin) zprvs(iv)=0.0_wp
        IF (zprvi(iv) <= zqmin) zprvi(iv)=0.0_wp

#ifdef NUDGING
        ! for the latent heat nudging
        IF ((llhn .OR. llhnverif) .AND. lhn_qrs) THEN
          IF (lsedi_ice .OR. lorig_icon) THEN
            qrsflux(iv,k) = zprvr(iv)+zprvs(iv)+zprvi(iv)
            qrsflux(iv,k) = 0.5_wp*(qrsflux(iv,k)+zpkr(iv)+zpks(iv)+zpki(iv)) 
          ELSE 
            qrsflux(iv,k) = zprvr(iv)+zprvs(iv)
            qrsflux(iv,k) = 0.5_wp*(qrsflux(iv,k)+zpkr(iv)+zpks(iv))
          ENDIF
        ENDIF
#endif

        IF (qrg+qr(iv,k+1) <= zqmin) THEN
          zvzr(iv)= 0.0_wp
        ELSE
          zvzr(iv)= zvz0r * EXP(zvzxp*LOG((qrg+qr(iv,k+1))*0.5_wp*rhog)) * zrho1o2
        ENDIF
        IF (qsg+qs(iv,k+1) <= zqmin) THEN
          zvzs(iv)= 0.0_wp
        ELSE
          zvzs(iv)= zvz0s * EXP(zv1s/(zbms+1.0_wp)*LOG((qsg+qs(iv,k+1))*0.5_wp*rhog)) * zrho1o2
        ENDIF
        IF (qig+qi(iv,k+1) <= zqmin ) THEN
          zvzi(iv)= 0.0_wp
        ELSE
          !! density correction not needed
          zvzi(iv)= zvz0i * EXP(zbvi*LOG((qig+qi(iv,k+1))*0.5_wp*rhog))
        ENDIF

      ELSE
        ! Precipitation fluxes at the ground
        prr_gsp(iv) = 0.5_wp * (qrg*rhog*zvzr(iv) + zpkr(iv))
        IF (lsedi_ice .OR. lorig_icon) THEN
          prs_gsp(iv) = 0.5_wp * (rhog*(qsg*zvzs(iv)+qig*zvzi(iv)) + zpks(iv)+zpki(iv))
        ELSE
          prs_gsp(iv) = 0.5_wp * (qsg*rhog*zvzs(iv) + zpks(iv))
        END IF

#ifdef NUDGING
        ! for the latent heat nudging
        IF ((llhn .OR. llhnverif) .AND. lhn_qrs) THEN
          qrsflux(iv,k) = prr_gsp(iv)+prs_gsp(iv)
        ENDIF
#endif
      ENDIF

      ! Update of prognostic variables or tendencies
      qr (iv,k) = qrg
      qs (iv,k) = MAX ( 0.0_wp, qsg )
      qi (iv,k) = qig
      !        qrs(iv,k   ) = qrg+qsg+qig       !qrs is now computed outside
      t  (iv,k) = t (iv,k) + ztt*zdt 
      qv (iv,k) = MAX ( 0.0_wp, qv(iv,k) + zqvt*zdt )
      qc (iv,k) = MAX ( 0.0_wp, qc(iv,k) + zqct*zdt )

      IF (izdebug > 15) THEN
        ! Check for negative values
        IF (qr(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(A)') ' WARNING: cloudice, negative value in qr'
          CALL message('',message_text)
        ENDIF
        IF (qc(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(A)') ' WARNING: cloudice, negative value in qc'
          CALL message('',message_text)
        ENDIF
        IF (qi(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(A)') ' WARNING: cloudice, negative value in qi'
          CALL message('',message_text)
        ENDIF
        IF (qs(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(A)') ' WARNING: cloudice, negative value in qs'
          CALL message('',message_text)
        ENDIF
        IF (qv(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(A)') ' WARNING: cloudice, negative value in qv'
          CALL message('',message_text)
        ENDIF
      ENDIF

    ENDDO !loop over iv
  ENDDO loop_over_levels

#if defined (__COSMO__)
! Do a final saturation adjustment for new values of t, qv and qc
DO  k = k_start, ke
    CALL satad ( 1, t(:,k), qv(:,k),                   &
               qc(:,k), t(:,k), p(:,k),                &
               zdummy(:,1),zdummy(:,2),zdummy(:,3),    &
               zdummy(:,4),zdummy(:,5),zdummy(:,6),    &
               zdummy(:,7),zdummy(:,8),                &
               b1, b2w, b3, b4w, b234w, rdv, o_m_rdv,  &
               rvd_m_o, lh_v, z_heat_cap_r, cp_d,      &
               nvec, 1, iv_start, iv_end, 1 , 1 )
ENDDO !loop over k


DO  k = k_start, ke
  DO iv=iv_start,iv_end

    IF ( ldiabf_lh ) THEN
      ! compute temperature increment due to latent heat
      tinc_lh(iv,k) = tinc_lh(iv,k) + t(iv,k)
    ENDIF
#ifdef NUDGING
    ! add part of latent heating calculated in subroutine hydci to model latent
    ! heating field: add temperature to model latent heating field
    IF (llhn .OR. llhnverif) THEN
      ! CALL get_gs_lheating ('inc',1,ke)  !XL :this should be called from within the block
      tt_lheat(iv,k) = tt_lheat(iv,k) + t(iv,k)
    ENDIF
#endif

  END DO
END DO
#endif

#ifdef __ICON__

CALL satad_v_3d (                             &
               & maxiter  = 10_i4    ,& !> IN
               & tol      = 1.e-3_wp ,& !> IN
               & te       = t        ,&
               & qve      = qv       ,&
               & qce      = qc       ,&
               & rhotot   = rho      ,&
               & idim     = nvec     ,&
               & kdim     = ke       ,&
               & ilo      = iv_start ,&
               & iup      = iv_end   ,&
               & klo      = k_start  ,&
               & kup      = ke        &
               )

#endif

!------------------------------------------------------------------------------
! final tendency calculation for ICON
!
! Note: as soon as we have a new satad subroutine in ICON, this tendency
! calculation will be done in the k-loop and the original 3D variables wont
! be used to store the new values. Then we wont need the _in variables anymore.
!------------------------------------------------------------------------------

  IF (PRESENT(ddt_tend_t)) THEN

    DO k=k_start,ke
       DO iv=iv_start,iv_end

          ! calculated pseudo-tendencies
          ddt_tend_t (iv,k) = (t (iv,k) - t_in (iv,k))*zdtr
          ddt_tend_qv(iv,k) = MAX(-qv_in(iv,k)*zdtr,(qv(iv,k) - qv_in(iv,k))*zdtr)
          ddt_tend_qc(iv,k) = MAX(-qc_in(iv,k)*zdtr,(qc(iv,k) - qc_in(iv,k))*zdtr)
          ddt_tend_qr(iv,k) = MAX(-qr_in(iv,k)*zdtr,(qr(iv,k) - qr_in(iv,k))*zdtr)
          ddt_tend_qs(iv,k) = MAX(-qs_in(iv,k)*zdtr,(qs(iv,k) - qs_in(iv,k))*zdtr)
          ddt_tend_qi(iv,k) = MAX(-qi_in(iv,k)*zdtr,(qi(iv,k) - qi_in(iv,k))*zdtr)

      END DO
    END DO

  END IF

  IF (izdebug > 15) THEN
    CALL message('gscp_cloudice', 'UPDATED VARIABLES')
   WRITE(message_text,'(A,2E20.9)') 'cloudice T= ',&
    MAXVAL( t(:,:)), MINVAL(t(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(A,2E20.9)') 'cloudice qv= ',&
    MAXVAL( qv(:,:)), MINVAL(qv(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(A,2E20.9)') 'cloudice qc= ',&
    MAXVAL( qc(:,:)), MINVAL(qc(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(A,2E20.9)') 'cloudice qi= ',&
    MAXVAL( qi(:,:)), MINVAL(qi(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(A,2E20.9)') 'cloudice qr= ',&
    MAXVAL( qr(:,:)), MINVAL(qr(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(A,2E20.9)') 'cloudice qs= ',&
    MAXVAL( qs(:,:)), MINVAL(qs(:,:) )
    CALL message('', TRIM(message_text))
  ENDIF

!------------------------------------------------------------------------------
! End of subroutine cloudice
!------------------------------------------------------------------------------

END SUBROUTINE cloudice

!==============================================================================

END MODULE gscp_cloudice
