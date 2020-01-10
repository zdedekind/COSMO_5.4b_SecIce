!+ Data module for variables of the convection parameterizations
!------------------------------------------------------------------------------

MODULE conv_data

!------------------------------------------------------------------------------
!
! Description:
!  This module contains variables that are used in the various convection
!  parameterizations.
!
! Current Code Owner: DWD, Dmitrii Mironov
!  phone:  +49  69  8062 2705
!  fax:    +49  69  8062 3721
!  email:  Dmitrii.Mironov@dwd.de
!
! History from former module data_convection:
! -------------------------------------------
! V4_5         2008/09/10 Ulrich Schaettler
!  Initial release
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_18        2011/05/26 Ulrich Schaettler
!  New variable thick_sc (by Martin Koehler)
!  Changed the code owner
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V5_4b        2016-07-12 Ulrich Schaettler
!  Initial release for blocked version of convection in COSMO
!
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

USE kind_parameters, ONLY :   &
    wp           ! KIND-type parameter for real variables

!==============================================================================

IMPLICIT NONE

!==============================================================================

  ! The following parameters are tunable constants for the cumulus convection
  ! scheme. 

  REAL (KIND=wp)     ::           &
    ! mean entrainment rate for shallow convection
    entr_sc  = 0.00030_wp,        &

    ! limit for convective clouds to be "shallow" (in Pa)
    ! Shallow convection parameterization becomes active only if cloud
    ! thickness from cloud base to cloud top exceeds a threshold.  To evaluate
    ! this condition a parcel is launched.  This threshold is typically set to
    ! values between 200hPa and 300hPa with a COSMO DE default of 250hPa.
    thick_sc = 25000.0_wp,        &

    ! COSMO-DE default (by Guenther Doms) : 250 hPa
    ! IFS      default (by Peter Bechtold): 200 hPa
    ! reasonable values:  between 100 hPa and 450 hPa

    ! Threshold of cloud thickness for precipitating deep convection mode in 
    ! Tiedtke-Bechtold deep convection parameterization
    thick_dc = 20000.0_wp

  LOGICAL                     :: &
    lmfmid        = .TRUE.    ,  & ! switch for mid-level convection
    lmfdd         = .TRUE.    ,  & ! switch for inclusion of downdrafts
    lmfdudv       = .TRUE.         ! switch for cumulus effects on momentum

  REAL (KIND=wp), PARAMETER   :: &
    entrpen       = 0.00010_wp,  & ! mean entrainment rate for deep convection
    entrmid       = 0.00010_wp,  & ! mean entrainment rate for mid-level convection
    entrdd        = 0.00020_wp,  & ! mean entrainment rate for downdrafts
    entr_pb       = 1.E-04_wp ,  & ! constant in extended formulation (DM)
                                 ! turbulent entrainment/detrainment rate

    cmfcmax       = 1.0_wp    ,  & ! maximum mass flux
    cmfcmin       = 1.E-10_wp ,  & ! minimum mass flux  (security)
    cmfctop       = 0.33_wp   ,  & ! relative mass flux above level of non-buoyancy
    cmfdeps       = 0.3_wp    ,  & ! relative mass flux at level of free-sinking

    capeconstant  = 0.05_wp   ,  & ! time-constant cu_base massflux from CAPE
    capemin       = 1.0_wp    ,  & ! minimum CAPE for cu_base massflux determination
    cmbtke        = 1.0_wp    ,  & ! time-constant cu_base massflux from ctke
    ctkemin       = 0.05_wp   ,  & ! mimimum ctke for cu_base massflux determination

    cprcon        = 0.0002_wp      ! conversion rate cloud water to precipitation

  ! Set Tmpmin, Tmpmax and exp_mp as in the ECMWF IFS (CY31r1).
  REAL (KIND=wp), PARAMETER   :: &
    Tmpmin        = 250.16_wp ,  & ! Minium temperature of the mixed-phase temperature range [K]
                                 ! US: Tmpmin = ctmelt!!!!!!
    Tmpmax        = 273.16_wp ,  & ! Maximum temperature of the mixed-phase temperature range [K]
    exp_mp        = 2.0_wp         ! Exponent in the interpolation formula
                                 ! for the mixed-phase water fraction [-]

  REAL (KIND=wp)              :: &
    cu_frr                         ! fraction of grid area covered with precipitation

  REAL (KIND=wp), ALLOCATABLE :: &
    cu_evap (:)                    ! factor for evaporation of rain

!------------------------------------------------------------------------------

  ! Other parameters for Tiedtke and shallow convection
  REAL (KIND = wp)        :: &
    ctmelt               , & ! tripel point
    cc2                  , & ! definition of utitility constants
    c5hlccp              , & ! for saturation humidity
    c5hlscp              , & !
    chlcdcp              , & !
    chlsdcp                  !

!------------------------------------------------------------------------------

  ! Other parameters for Tiedtke-Bechtold convection
  LOGICAL, PARAMETER ::                                 &
    ltuning_kessler = .TRUE.

  INTEGER ::                                            &
    ! type of coupling between aerosols and convection scheme
    icpl_aero_conv

  REAL (KIND = wp)        :: &
    ! Fraction of CAPE diurnal cycle correction applied in the extratropics 
    ! (relevant only if icapdcycl = 3)
    tune_capdcfac_et = 0.0_wp      , &

    ! relative humidity threshold for onset of evaporation below cloud base ...
    tune_rhebc_land  = 0.75_wp     , & ! ... over land
    tune_rhebc_ocean = 0.85_wp     , & ! ... over sea

    ! Excess value for temperature used in test parcel ascent
    tune_texc        = 0.125_wp    , & !

    ! Excess fraction of grid-scale QV used in test parcel ascent
    tune_qexc        = 1.25e-2_wp  , & !

    ! Convective area fraction
    tune_rcucov      = 0.05_wp     , & !

    ! Entrainment parameter for deep convection valid at dx=20 km 
    tune_entrorg     = 1.825e-3_wp , & !

    ! The following switches allow separate tuning for evaporation 
    ! below cloud base in the tropics

    ! relative humidity threshold for onset of evaporation below cloud base ...
    tune_rhebc_land_trop  = 0.70_wp, & ! ... over tropical land 
                                       ! (relevant only if smaller than rhebc_land)
    tune_rhebc_ocean_trop = 0.80_wp, & ! ... over tropical sea  
                                       ! (relevant only if smaller than rhebc_ocean)

    ! Convective area fraction in the tropics
    tune_rcucov_trop      = 0.05_wp    ! (relevant only if smaller than rcucov)



  ! Definition of a type for the interface of the Tiedtke-Bechtold scheme
  ! This type holds the above tuning constants
  TYPE t_phy_params
    INTEGER  :: kcon1, kcon2           ! Level parameters for convection scheme
    REAL(wp) :: tau, mfcfl, tau0       ! resolution-dependent parameters for convection scheme

    REAL(wp) :: rhebc_land, rhebc_ocean, rhebc_land_trop, rhebc_ocean_trop
                                       ! relative humidity below which sub-cloud evaporation of rain starts
                                       ! over land and water, respectively

    REAL(wp) :: texc, qexc             ! 'excess values' of temperature and QV used for convection
                                       ! triggering (test parcel ascent)

    REAL(wp) :: rcucov, rcucov_trop    ! fractional area covered by convective precipitation

    REAL(wp) :: entrorg                ! tuning coefficient for organized entrainment of deep convection

    LOGICAL  :: lmfscv, lmfmid, lmfpen ! switches for activation of shallow, midlevel and deep convection

  END TYPE t_phy_params

  TYPE (t_phy_params) :: phy_params

!------------------------------------------------------------------------------

  ! organization of tracers which undergo passive transport by convection
  ! this is for the Tiedtke-scheme
  INTEGER :: &
    ntrcr_con      ! number of tracers 

  INTEGER, ALLOCATABLE :: &
    itrcr_con(:)  ! Index list of tracers

  ! this should perhaps be put to data_block_fields
  REAL (KIND=wp),     ALLOCATABLE ::  &
    trcr_tran_b (:,:,:), & ! tracers for convective transport at full levels
    trcr_tend_b (:,:,:)    ! convective tendency of that tracers


  ! And this is for the Tiedtke-Bechtold scheme: TODO: Have to unify these approaches
!_JF: TODO: implement convective transport of ART tracers via IFS-/ICON-convection
!_JF:       temporary dummy definition of tracer structure used in ICON
  ! Type to pass pointer arrays to convection (and turbulent diffusion) subroutines
  TYPE t_ptr_tracer
    REAL(wp), POINTER :: ptr(:,:)
    INTEGER           :: idx_tracer
  END TYPE t_ptr_tracer

!=======================================================================

#ifdef _OPENACC
  ! Declarations of working  arrays for the shallow convection
  INTEGER, ALLOCATABLE :: &
    mlab   (:,:)    ,   & !
    mtype  (:)      ,   & !
    ! And fields that are also necessary for the Tiedtke convection
    mclab  (:,:)    ,   & !
    mdtop  (:)      ,   & !
    mctop0 (:)      ,   & !
    mlwmin (:)      ,   & !
    mtmelt (:)      ,   & ! melting level
    kcptop (:)            ! highest unstable layer in CAPE-determination

  REAL (KIND=wp), ALLOCATABLE :: &
    ztu    (:,:)    ,   & !
    zqu    (:,:)    ,   & !
    zlu    (:,:)    ,   & !
    zlude  (:,:)    ,   & !
    zmfus  (:,:)    ,   & !
    zmfuq  (:,:)    ,   & !
    zmful  (:,:)    ,   & !
    zqsen  (:,:)    ,   & ! saturation specific humidiy at full levels
    ztenh  (:,:)    ,   & !
    zqenh  (:,:)    ,   & !
    zqsenh (:,:)    ,   & !
    zrhgdz (:,:)    ,   & ! rho*g*dz
    z1dp   (:)      ,   & !
    zrho   (:)      ,   & !
    zqold  (:)      ,   & !
    zentr  (:)      ,   & !
    zmfub  (:)      ,   & !
    zdqpbl (:)      ,   & !
    zdmfen (:)      ,   & !
    zdmfde (:)            !

  LOGICAL, ALLOCATABLE  :: &
    loflag (:)            !

  ! And fields that are also necessary for the Tiedtke convection
  REAL (KIND=wp), ALLOCATABLE :: &
    zuu    (:,:)    ,   & !
    zvu    (:,:)    ,   & !
    zmfu   (:,:)    ,   & !
    zdmfup (:,:)    ,   & !
    ztd    (:,:)    ,   & !
    zqd    (:,:)    ,   & !
    zud    (:,:)    ,   & !
    zvd    (:,:)    ,   & !
    zmfd   (:,:)    ,   & !
    zmfds  (:,:)    ,   & !
    zmfdq  (:,:)    ,   & !
    zdmfdp (:,:)    ,   & !
    zcptu  (:,:)    ,   & ! updraft temperature in CAPE-determination
    zcpqu  (:,:)    ,   & ! updraft humidity in CAPE-determination
    zmfuu  (:,:)    ,   & ! mass flux up: u
    zmfuv  (:,:)    ,   & ! mass flux up: v
    zmfdu  (:,:)    ,   & ! mass flux down: u
    zmfdv  (:,:)    ,   & ! mass flux down: v
    zqenwb (:,:)    ,   & !
    ztenwb (:,:)    ,   & !
    zcape  (:)      ,   & ! convective available energy CAPE
    zcond  (:)      ,   & !
    zrfl   (:)      ,   & !
    zrneva (:)      ,   & !
    zhcbase(:)      ,   & !
    zmfub1 (:)      ,   & !
    zdqcv  (:)      ,   & !
    zwmax  (:)            !

  ! Fields necessary for tracer support in Tiedtke
  REAL (KIND=wp), ALLOCATABLE :: &
    ztru   (:,:,:)  ,   & ! tracer updraft
    ztrd   (:,:,:)  ,   & ! tracer downdraft
    zmfutr (:,:,:)  ,   & ! mass flux up: tracers
    zmfdtr (:,:,:)        ! mass flux down: tracers
#endif

!==============================================================================
! Module procedures in conv_data
!==============================================================================

CONTAINS

!==============================================================================
#ifdef _OPENACC
!==============================================================================
!+ Module procedure for allocation of working arrays
!------------------------------------------------------------------------------

SUBROUTINE conv_wkarr_alloc ( nproma, ke , itype_conv, ierror)

!------------------------------------------------------------------------------
!
! Purpose:   Allocation of the working arrays required in the subroutine
!            cu_shallow
!
! Method:
!
!------------------------------------------------------------------------------

! Declarations:
INTEGER, INTENT(IN)  :: nproma, ke, itype_conv
INTEGER, INTENT(OUT) :: ierror

INTEGER :: izl, ist


! End of header
!==============================================================================

  ierror = 0
  izl = 0
  ist = 0

  ALLOCATE ( mlab      (nproma, ke) ,STAT=izl ); mlab    = 0     ; ist=ist+izl
  ALLOCATE ( mtype     (nproma)     ,STAT=izl ); mtype   = 0     ; ist=ist+izl

  IF (itype_conv == 0) THEN
    ALLOCATE ( mclab   (nproma, ke) ,STAT=izl ); mclab   = 0     ; ist=ist+izl
    ALLOCATE ( mdtop   (nproma)     ,STAT=izl ); mdtop   = 0     ; ist=ist+izl
    ALLOCATE ( mctop0  (nproma)     ,STAT=izl ); mctop0  = 0     ; ist=ist+izl
    ALLOCATE ( mlwmin  (nproma)     ,STAT=izl ); mlwmin  = 0     ; ist=ist+izl
    ALLOCATE ( mtmelt  (nproma)     ,STAT=izl ); mtmelt  = 0     ; ist=ist+izl
    ALLOCATE ( kcptop  (nproma)     ,STAT=izl ); kcptop  = 0     ; ist=ist+izl
  ENDIF

  ALLOCATE ( ztu       (nproma, ke) ,STAT=izl ); ztu     = 0.0_wp; ist=ist+izl
  ALLOCATE ( zqu       (nproma, ke) ,STAT=izl ); zqu     = 0.0_wp; ist=ist+izl
  ALLOCATE ( zlu       (nproma, ke) ,STAT=izl ); zlu     = 0.0_wp; ist=ist+izl
  ALLOCATE ( zlude     (nproma, ke) ,STAT=izl ); zlude   = 0.0_wp; ist=ist+izl
  ALLOCATE ( zmfus     (nproma, ke) ,STAT=izl ); zmfus   = 0.0_wp; ist=ist+izl
  ALLOCATE ( zmfuq     (nproma, ke) ,STAT=izl ); zmfuq   = 0.0_wp; ist=ist+izl
  ALLOCATE ( zmful     (nproma, ke) ,STAT=izl ); zmful   = 0.0_wp; ist=ist+izl
  ALLOCATE ( zqsen     (nproma, ke) ,STAT=izl ); zqsen   = 0.0_wp; ist=ist+izl
  ALLOCATE ( ztenh     (nproma, ke) ,STAT=izl ); ztenh   = 0.0_wp; ist=ist+izl
  ALLOCATE ( zqenh     (nproma, ke) ,STAT=izl ); zqenh   = 0.0_wp; ist=ist+izl
  ALLOCATE ( zqsenh    (nproma, ke) ,STAT=izl ); zqsenh  = 0.0_wp; ist=ist+izl
  ALLOCATE ( zrhgdz    (nproma, ke) ,STAT=izl ); zrhgdz  = 0.0_wp; ist=ist+izl
  ALLOCATE ( z1dp      (nproma)     ,STAT=izl ); z1dp    = 0.0_wp; ist=ist+izl
  ALLOCATE ( zrho      (nproma)     ,STAT=izl ); zrho    = 0.0_wp; ist=ist+izl
  ALLOCATE ( zqold     (nproma)     ,STAT=izl ); zqold   = 0.0_wp; ist=ist+izl
  ALLOCATE ( zentr     (nproma)     ,STAT=izl ); zentr   = 0.0_wp; ist=ist+izl
  ALLOCATE ( zmfub     (nproma)     ,STAT=izl ); zmfub   = 0.0_wp; ist=ist+izl
  ALLOCATE ( zdqpbl    (nproma)     ,STAT=izl ); zdqpbl  = 0.0_wp; ist=ist+izl
  ALLOCATE ( zdmfen    (nproma)     ,STAT=izl ); zdmfen  = 0.0_wp; ist=ist+izl
  ALLOCATE ( zdmfde    (nproma)     ,STAT=izl ); zdmfde  = 0.0_wp; ist=ist+izl
  ALLOCATE ( loflag    (nproma)     ,STAT=izl ); loflag  =.FALSE.; ist=ist+izl

  IF (itype_conv == 0) THEN
    ALLOCATE ( zuu     (nproma,ke)  ,STAT=izl ); zuu     = 0.0_wp; ist=ist+izl
    ALLOCATE ( zvu     (nproma,ke)  ,STAT=izl ); zvu     = 0.0_wp; ist=ist+izl
    ALLOCATE ( zmfu    (nproma,ke)  ,STAT=izl ); zmfu    = 0.0_wp; ist=ist+izl
    ALLOCATE ( zdmfup  (nproma,ke)  ,STAT=izl ); zdmfup  = 0.0_wp; ist=ist+izl
    ALLOCATE ( ztd     (nproma,ke)  ,STAT=izl ); ztd     = 0.0_wp; ist=ist+izl
    ALLOCATE ( zqd     (nproma,ke)  ,STAT=izl ); zqd     = 0.0_wp; ist=ist+izl
    ALLOCATE ( zud     (nproma,ke)  ,STAT=izl ); zud     = 0.0_wp; ist=ist+izl
    ALLOCATE ( zvd     (nproma,ke)  ,STAT=izl ); zvd     = 0.0_wp; ist=ist+izl
    ALLOCATE ( zmfd    (nproma,ke)  ,STAT=izl ); zmfd    = 0.0_wp; ist=ist+izl
    ALLOCATE ( zmfds   (nproma,ke)  ,STAT=izl ); zmfds   = 0.0_wp; ist=ist+izl
    ALLOCATE ( zmfdq   (nproma,ke)  ,STAT=izl ); zmfdq   = 0.0_wp; ist=ist+izl
    ALLOCATE ( zdmfdp  (nproma,ke)  ,STAT=izl ); zdmfdp  = 0.0_wp; ist=ist+izl
    ALLOCATE ( zcptu   (nproma,ke)  ,STAT=izl ); zcptu   = 0.0_wp; ist=ist+izl
    ALLOCATE ( zcpqu   (nproma,ke)  ,STAT=izl ); zcpqu   = 0.0_wp; ist=ist+izl
    ALLOCATE ( zmfuu   (nproma,ke)  ,STAT=izl ); zmfuu   = 0.0_wp; ist=ist+izl
    ALLOCATE ( zmfuv   (nproma,ke)  ,STAT=izl ); zmfuv   = 0.0_wp; ist=ist+izl
    ALLOCATE ( zmfdu   (nproma,ke)  ,STAT=izl ); zmfdu   = 0.0_wp; ist=ist+izl
    ALLOCATE ( zmfdv   (nproma,ke)  ,STAT=izl ); zmfdv   = 0.0_wp; ist=ist+izl
    ALLOCATE ( zqenwb  (nproma,ke)  ,STAT=izl ); zqenwb  = 0.0_wp; ist=ist+izl
    ALLOCATE ( ztenwb  (nproma,ke)  ,STAT=izl ); ztenwb  = 0.0_wp; ist=ist+izl
    ALLOCATE ( zcape   (nproma)     ,STAT=izl ); zcape   = 0.0_wp; ist=ist+izl
    ALLOCATE ( zcond   (nproma)     ,STAT=izl ); zcond   = 0.0_wp; ist=ist+izl
    ALLOCATE ( zrfl    (nproma)     ,STAT=izl ); zrfl    = 0.0_wp; ist=ist+izl
    ALLOCATE ( zrneva  (nproma)     ,STAT=izl ); zrneva  = 0.0_wp; ist=ist+izl
    ALLOCATE ( zhcbase (nproma)     ,STAT=izl ); zhcbase = 0.0_wp; ist=ist+izl
    ALLOCATE ( zmfub1  (nproma)     ,STAT=izl ); zmfub1  = 0.0_wp; ist=ist+izl
    ALLOCATE ( zdqcv   (nproma)     ,STAT=izl ); zdqcv   = 0.0_wp; ist=ist+izl
    ALLOCATE ( zwmax   (nproma)     ,STAT=izl ); zwmax   = 0.0_wp; ist=ist+izl


    ALLOCATE ( ztru    (nproma,ke,ntrcr_con), STAT=izl ); ztru    = 0.0_wp; ist=ist+izl
    ALLOCATE ( ztrd    (nproma,ke,ntrcr_con), STAT=izl ); ztrd    = 0.0_wp; ist=ist+izl
    ALLOCATE ( zmfutr  (nproma,ke,ntrcr_con), STAT=izl ); zmfutr  = 0.0_wp; ist=ist+izl
    ALLOCATE ( zmfdtr  (nproma,ke,ntrcr_con), STAT=izl ); zmfdtr  = 0.0_wp; ist=ist+izl
  ENDIF

  IF (ist /= 0) THEN
    PRINT *, '*** ERROR: Allocation of working arrays for convection failed ***'
    ierror = 1
  ENDIF

  !$acc enter data                                                              &
  !$acc create (mlab, mtype, ztu, zqu, zlu, zlude, zmfus, zmfuq)                &
  !$acc create (zmful, zqsen, ztenh, zqenh, zqsenh, zrhgdz, z1dp, zrho)         &
  !$acc create (zqold, zentr, zmfub, zdqpbl, zdmfen, zdmfde, loflag)

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE conv_wkarr_alloc

!==============================================================================
!+ Module procedure for deallocation of working arrays
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------

SUBROUTINE conv_wkarr_dealloc (itype_conv, ierror)

!------------------------------------------------------------------------------
!
! Purpose:   Deallocation of the working arrays required in the subroutine
!            cu_shallow
!
! Method:
!
!------------------------------------------------------------------------------

! Declarations:

INTEGER, INTENT(IN)  :: itype_conv
INTEGER, INTENT(OUT) :: ierror

INTEGER :: izl, ist


! End of header
!==============================================================================

  ierror = 0
  izl = 0
  ist = 0

  !$acc exit data                                                             &
  !$acc delete (mlab, mtype, ztu, zqu, zlu, zlude, zmfus, zmfuq)              &
  !$acc delete (zmful, zqsen, ztenh, zqenh, zqsenh, zrhgdz, z1dp, zrho)       &
  !$acc delete (zqold, zentr, zmfub, zdqpbl, zdmfen, zdmfde, loflag)

  DEALLOCATE ( mlab     , STAT=izl );  ist=ist+izl
  DEALLOCATE ( mtype    , STAT=izl );  ist=ist+izl

  IF (itype_conv == 0) THEN
  DEALLOCATE ( mclab    , STAT=izl );  ist=ist+izl
  DEALLOCATE ( mdtop    , STAT=izl );  ist=ist+izl
  DEALLOCATE ( mctop0   , STAT=izl );  ist=ist+izl
  DEALLOCATE ( mlwmin   , STAT=izl );  ist=ist+izl
  DEALLOCATE ( mtmelt   , STAT=izl );  ist=ist+izl
  DEALLOCATE ( kcptop   , STAT=izl );  ist=ist+izl
  ENDIF
 
  DEALLOCATE ( ztu      , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zqu      , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zlu      , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zlude    , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zmfus    , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zmfuq    , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zmful    , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zqsen    , STAT=izl );  ist=ist+izl
  DEALLOCATE ( ztenh    , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zqenh    , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zqsenh   , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zrhgdz   , STAT=izl );  ist=ist+izl
  DEALLOCATE ( z1dp     , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zrho     , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zqold    , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zentr    , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zmfub    , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zdqpbl   , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zdmfen   , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zdmfde   , STAT=izl );  ist=ist+izl
  DEALLOCATE ( loflag   , STAT=izl );  ist=ist+izl

  IF (itype_conv == 0) THEN
  DEALLOCATE ( zuu      , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zvu      , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zmfu     , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zdmfup   , STAT=izl );  ist=ist+izl
  DEALLOCATE ( ztd      , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zqd      , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zud      , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zvd      , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zmfd     , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zmfds    , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zmfdq    , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zdmfdp   , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zcptu    , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zcpqu    , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zmfuu    , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zmfuv    , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zmfdu    , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zmfdv    , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zqenwb   , STAT=izl );  ist=ist+izl
  DEALLOCATE ( ztenwb   , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zcape    , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zcond    , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zrfl     , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zrneva   , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zhcbase  , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zmfub1   , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zdqcv    , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zwmax    , STAT=izl );  ist=ist+izl


  DEALLOCATE ( ztru     , STAT=izl );  ist=ist+izl
  DEALLOCATE ( ztrd     , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zmfutr   , STAT=izl );  ist=ist+izl
  DEALLOCATE ( zmfdtr   , STAT=izl );  ist=ist+izl
  ENDIF

  IF (ist /= 0) THEN
    PRINT *, '*** ERROR: Deallocation of working arrays for convection failed ***'
    ierror = 1
  ENDIF

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE conv_wkarr_dealloc

!==============================================================================
#endif
!==============================================================================

END MODULE conv_data
