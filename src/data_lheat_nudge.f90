!+ Data module for global fields needed for latent heat nudging
!-------------------------------------------------------------------------------

MODULE data_lheat_nudge

!-------------------------------------------------------------------------------
!
! Description:
!  This module declares all namelist variables and fields needed for the
!  latent heat nudging which are also needed by other modules or have to
!  be held in the long term storage.
!  Variables and fields only needed temporarily for the lhn are declared in the
!  latent heat nudging module.
!  Remark:
!  This module is only meant as a temporary solution for more flexible tests
!  of lhn and should be included later into : data_nudge_all.f90  to have one
!  common namelist for all parameters needed for tasks related to 
!  data assimilation (change corresponding USE statements !!).
!    - NAMELIST variables
!    - 3-dim. array for the storage of model latent heating profiles
!    - 2-dim. arrays for the storage of observations
!
!  The fields are declared as allocatable arrays. They are allocated in the
!  lhn module at its first call.
!
! Current Code Owner: DWD, Klaus Stephan
!  phone:  +49  69  8062 2689
!  fax:    +49  69  8062 3721
!  email:  klaus.stephan@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 3.13       2004/12/03 Christina Koepken
!  Initial version
! 3.18       2006/03/03 Klaus Stephan
!  Namelist parameter moved from data_nudge_all to avoid to many dependencies
!  New varibles for integrated precipitation flux and radar blacklist
! 3.21       2006/12/04 Klaus Stephan
!  New Namelist variable lhn_wweight to apply weighting with respect to mean
!  horizontal wind.
!  New variables pr_mod_sum, pr_ref_sum for computations on cumulated precipitation
! V3_24        2007/04/26 Klaus Stephan
!  Eliminated Namelist variables lhn_diag, rlhn_scale_dp
! V3_28        2007/07/18 Klaus Stephan
!  New Namelist variables for
!  - data input: noobs_date,n_noobs
!  - reference precipitaion: rqrsgflux
!  - vertical restriction: ktop_temp
! V4_1         2007/12/04 Klaus Stephan
!  New Namelist variables and variables for LHN:
!  - use radar beam height: lhn_height, height_file, dxheight(:,:)
!  - bright band detection algorithm: lhn_bright, brightband(:,:)
! V4_12        2010/05/11 Daniel Leuenberger
!  - removed namelist switch lhn_radar_dx
!  - added namelist switch lhn_spqual       (default value: false)
!  - added namelist lhn_dt_obs              (default vaule: 5min)
!  - added namelist nradar                  (default value: 32)
!  - replaced field dist by field spqual
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  Initialize lhn_qrs with FALSE. It is only set in organize_assimilation and is
!     not set, if luseobs=.FALSE. or if model is compiled without NUDGING
!  MESSy interface added
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Eliminated MESSy interface, because MESSy cannot be used with data assimilation
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!===============================================================================
!
! Declarations:
!
! Modules used:
USE data_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    iintegers, & ! KIND-type parameter for standard integer variables
    irealgrib, & ! KIND-type parameter for real variables in the grib library
    intgribf     ! KIND-type parameter for fortran files in the grib library

!==============================================================================

IMPLICIT NONE

!===============================================================================

! 1. Namelist variables controlling the latent heat nudging
! ---------------------------------------------------------

  LOGICAL                          ::          &
    llhn         ,& ! on/off switch for latent heat nudging (lhn)
    llhnverif    ,& ! on/off switch for verification against radar
    lhn_search   ,& ! search for appropriate nearby model heating profile
    lhn_filt     ,& ! vertical filtering of lhn t-increments
    lhn_relax    ,& ! horizontal filtering of lhn t-increments
    lhn_limit    ,& ! impose an absolute limit on lhn t-increments (abs_lhn_lim)
    lhn_hum_adj  ,& ! apply a humidity adjustment along with t-increments
    lhn_incloud  ,& ! apply temperature increment only in cloudy levels
    lhn_spqual   ,& ! switch for the use of a spatial quality function
    lhn_black    ,& ! use blacklist for radar data
    lhn_height   ,& ! use height infos for radar data
    lhn_bright   ,& ! apply bright band detection
    lhn_logscale ,& ! apply logarithmic scaling factors
    lhn_wweight  ,& ! apply a weighting with respect to the mean horizontal wind
    lhn_diag        ! produce more detailed diagnostic output during lhn

  LOGICAL                          ::          &
    lhn_qrs=.FALSE. ! calculate the integrated precipitation flux
                    ! (value has to be initialized when src_gscp is called)

  INTEGER (KIND=iintegers)         ::          &
    nlhn_start      ,& ! start of latent heat nudging period in timesteps
    nlhn_end        ,& ! end of latent heat nudging period in timesteps
    nlhnverif_start ,& ! start of latent heat nudging period in timesteps
    nlhnverif_end   ,& ! end of latent heat nudging period in timesteps
    rlhn_search     ,& ! radius (gridpoints) for profiles search (if lhn_search)
    nlhn_relax      ,& ! number of iteration of the horizontal filter
    ktop_lhn        ,& ! index for uppest model layer for which lhn is performed
    kbot_lhn        ,& ! index for lowest model layer for which lhn is performed
    nradar             ! max. number of radar stations within input data

  REAL (KIND=wp)                   ::           &
    lhn_coef        ,& ! factor for reduction of lhn t-increments
    abs_lhn_lim     ,& ! absolute limit for lhn t-increments (imposed if lhn_limit)
    fac_lhn_search  ,& ! factor for search nearby profiles
    fac_lhn_up      ,& ! limiting factor for upscaling of model heating profile
    fac_lhn_down    ,& ! limiting factor for downscaling of model heating profile
    rad_wobs_lhn    ,& ! max. distance to radar for full observation weight
    thres_lhn       ,& ! threshold of rain rates to be consinderd within lhn approach
    rqrsgmax        ,& ! ratio of maximum of qrsgflux, needed for reference precipitation
    ktop_temp          ! temperature of uppest model layer for which lhn is performed

  REAL (KIND=wp), SAVE             ::           &
    lhn_dt_obs         ! time step of input data in minutes

  CHARACTER (LEN=100)              ::           &
    radar_in        ,& ! directory for reading radar-files and as well for blacklist_file
    blacklist_file  ,& ! name of blacklist file
    height_file        ! name of dxheight file

  INTEGER (KIND=iintegers)         ::           &
    nulhn              ! unit of lhn output file

  CHARACTER (LEN=7)       , PARAMETER  :: &
    yulhn = 'YULHN'     ! name of lhn output file

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    n_noobs = 36       ! max number of missing observations (12 per hour)
                       ! 36 for hlhn_end equal 3 hours

  CHARACTER (LEN=12)               ::           &
    noobs_date(n_noobs)

!-------------------------------------------------------------------------------

! Global (i.e. public) Declarations:

! 2. fields  and related variables                                     (unit)
! --------------------------------

  REAL  (KIND=wp), ALLOCATABLE ::               &
    tt_lheat(:,:,:,:) ,& ! profile of t-increments due to latent heating ( K/s )
                         ! (stored for current and previous timestep)
    tinc_lhn(:,:,:)   ,& ! temperature increments due to lhn             ( K/s )
    qrsflux  (:,:,:)     ! total precipitation flux

!! for testing output only !!
  REAL  (KIND=wp), ALLOCATABLE ::               &
    ttm_lheat(:,:,:)  ,& ! array for cumulated latent heating (grid scale + conv) ( K )
    ttm_cv(:,:,:)     ,& ! array for test output of diverse 2D fields
    tminc_lhn(:,:,:)  ,& ! cumulated temperature increments due to lhn   ( K )
    pr_obs_sum(:,:)   ,& ! cumulated precipitation (hourly)
    pr_mod_sum(:,:)   ,& ! cumulated precipitation (hourly)
    pr_ref_sum(:,:)      ! cumulated precipitation (hourly)
!! for testing output

  REAL  (KIND=wp), ALLOCATABLE ::               &
    obs(:,:,:)    ,& ! observations on model grid at six observation time levels
    spqual(:,:,:)   ,& ! spatial quality function on model grid at two observation time levels
    blacklist(:,:),& ! blacklist for DX radar data
    brightband(:,:),&! bright band mask field
    dxheight(:,:,:)  ! DX radar heights

! Arrays that are allocated during program execution
  INTEGER (KIND=intgribf), ALLOCATABLE :: &
    iblock_rad(:),      & ! array for gribed data
    ibmap_rad (:)         ! array for bit map

  REAL   (KIND=irealgrib), ALLOCATABLE :: &
    dsup_rad (:),       & ! array for special data
    ds_rad   (:)          ! array for unpacked data

!===============================================================================
END MODULE data_lheat_nudge
