!+ Interface module for organizing the radiation
!------------------------------------------------------------------------------

MODULE radiation_interface

!------------------------------------------------------------------------------
!
! Description:
!  The module radiation_interface contains routines to initialize and execute
!  radiation schemes when running the COSMO-model. Radiation schemes are in
!  the blocked data format (dimensions: nproma, ke, nblock) and the interface
!  takes care that all necessary input to the radiation schemes is transformed
!  from the COSMO-Model data format (ie,je,ke) to the blocked data format.
!
!  In contrast to the other parameterizations the "copy-to-block" scheme is 
!  not used, because we also support the usage of a coarse radiation grid,
!  which does not fit into this scheme. Instead, all routines and code parts
!  are written in a way, that they take the COSMO-Model grid as input and 
!  provide the output in the blocked data structure. Because the radiation
!  schemes need different input variables than the other parameterizations,
!  there is no double work necessary.
!
!  The only implemented radiation scheme at the moment is the Ritter-Geleyn
!  radiation. All routines of this scheme are grouped in the module 
!  radiation_rg.f90 (Subroutine fesft and routines called from fesft).
!
!  Routines included in radiation_interface:
!   - radiation_init:
!     initializations for the RG-scheme
!
!   - radiation_radcoarse_init
!     Additional initializations for the coarse grid (if nradcoarse > 1)
!
!   - radiation_organize
!     Computes input variables and calls the (RG) radiation scheme
!     This subroutine contains one other subroutine:
!   - average_to_coarse_grid
!     to average all input variables to the coarse radiation grid
!
!   - radiation_average (lzradstep)
!     Applies a discrete filter to several output variables, if lradf_avg is set
!
!   - radiation_wkarr_alloc (nproma_rad, ke, nradcoarse)
!     Allocation of work arrays
!
!   - radiation_wkarr_dealloc
!     Deallocation of work arrays
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  Ulrich.Schaettler@dwd.de
!
! History:
! Version      Date       Name
! ----------   ---------- ----
! V5_3         2015-10-09 Ulrich Schaettler
!  Initial release
! V5_3a        2015-11-24 Ulrich Schaettler
!  Adaptations in ACC code to run radiation on the GPUs
!  Minor adaptations to COSMO-ART
! V5_4         2016-03-10 Xavier Lapillonne, Heike Vogel (KIT)
!  Adaptations for running radiation on GPUs
!  Adaptations for COSMO-ART
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
    wp,        & ! KIND-type parameter for real variables
    sp,        & ! KIND-type parameter for real variables (single precision)
    dp           ! KIND-type parameter for real variables (double precision)

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------
    ie,           & ! number of grid points in zonal direction
    ie_tot,       & ! the same for the total field
    je_tot,       & ! the same for the total field
    je,           & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction
    ke1,          & ! ke+1
    czmls,        & ! depth of the soil layers in meters

! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from
!    the other ones because of the use of the staggered Arakawa-C-grid.

!   zonal direction
    istart,       & ! start index for the forecast of w, t, qd, qw and pp
    iend,         & ! end index for the forecast of w, t, qd, qw and pp
    istartpar,    & ! start index for computations in the parallel program
    iendpar,      & ! end index for computations in the parallel program

!   meridional direction
    jstart,       & ! start index for the forecast of w, t, qd, qw and pp
    jend,         & ! end index for the forecast of w, t, qd, qw and pp
    jstartpar,    & ! start index for computations in the parallel program
    jendpar,      & ! end index for computations in the parallel program

! 4. variables for the time discretization and related variables
! --------------------------------------------------------------
    dt,           & ! long time-step
    degrad,       & ! factor for transforming degree to rad

! 5. constants for the horizontal rotated grid and related variables
! ------------------------------------------------------------------

    pollon,       & ! longitude of the rotated north pole (in degrees, E>0)
    pollat,       & ! latitude of the rotated north pole (in degrees, N>0)

! 8. Organizational variables to handle the COSMO humidity tracers
! ----------------------------------------------------------------
    idt_qv,  idt_qc,  idt_qi

! end of data_modelconfig
!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &
    pi,           & ! circle constant
    cp_d,         & ! specific heat of dry air at constant pressure
    g,            & ! acceleration due to gravity
    sigma,        & ! Boltzmann-constant
    solc,         & ! solar constant
    rprecision      !

! end of data_constants
!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &
    p0          ,   & ! reference pressure at full levels             ( pa  )
    p0hl        ,   & ! reference pressure at half levels             ( Pa  )
    dp0         ,   & ! pressure thickness
    rlat        ,   & ! geographical latitude                         ( rad )
    rlon        ,   & ! geographical longitude                        ( rad )
    aerlan      ,   & ! aerosol-distribution for rural areas            --
    aerurb      ,   & ! aerosol-distribution for urban areas            --
    aerdes      ,   & ! aerosol-distribution for desert areas           --
    aersea      ,   & ! aerosol-distribution for sea                    --
    t           ,   & ! temperature                                   (  k  )
    pp          ,   & ! deviation from the reference pressure         ( pa  )
    ps          ,   & ! surface pressure                              ( pa  )
    t_g         ,   & ! weighted surface temperature                  (  k  )
    sohr        ,   & ! rate of solar heating                         ( k/s )
    sotr        ,   & ! solar transmissivity
    sotr_par    ,   & ! solar transmissivity, photosynthetic active radiation
!US these fields belong to a development which did never make it into the official
!   version. But if this development will be activated again (see also subroutine
!   radiation_average in radiation_interface) they are needed as global fields
!   flpar_s_dir ,   & ! direct component (aka parallel component)     ( W/m2)
!   flpar_s_difd,   & ! diffuse downward component                    ( W/m2)
!   flpar_s_difu,   & ! diffuse upward component                      ( W/m2)
    thhr        ,   & ! rate of thermal heating                       ( k/s )
    alb_rad     ,   & ! albedo of the ground                            --
    alb_rad_rc  ,   & ! albedo of the ground (on coarse radiation grid) --
    sobs        ,   & ! solar radiation at the ground                 ( w/m2)
    thbs        ,   & ! thermal radiation at the ground               ( w/m2)
    pabs        ,   & ! photosynthetic active radiation at the ground ( w/m2)
    sobt        ,   & ! solar radiation at the upper boundary         ( w/m2)
                      ! of the atmosphere
    thbt        ,   & ! thermal radiation at the upper boundary       ( w/m2)
                      ! of the atmosphere
    sun_el      ,   & ! sun elevation angle                           (deg  )
    sun_azi     ,   & ! sun azimuth  angle                            (deg  )

    ! and for the Climate-LM Version
    sodwddm     ,   & ! downward direct solar radiative flux / smu0   ( W/m2)

!   fields for the radiation correction scheme
    ! these are actual values
    swdir_s     ,   & ! direct comp. of solar radiative flux at surface ( W/m2)
    swdifd_s    ,   & ! diffuse downward comp. of short wave rad. flux  ( W/m2)
    swdifu_s    ,   & ! diffuse upward   comp. of short wave rad. flux  ( W/m2)
    swtrdir_s   ,   & ! direct comp. of solar radiative transmiss. at surface
    swtrdifd_s  ,   & ! diffuse downward comp. of short wave rad. transmiss.
    swtrdifu_s  ,   & ! diffuse upward   comp. of short wave rad. transmiss.
    lwd_s,lwu_s ,   & ! downward/upward comp. of long  wave rad. flux  ( W/m2)

    ! this is the essential correction factor
    swdir_cor   ,   & ! direct short wave radiation correction factor actual value

    ! these are topographic parameters
    skyview     ,   & ! sky view
    slo_asp     ,   & ! slope aspect
    slo_ang     ,   & ! slope angle
    horizon     ,   & ! horizon

! 7. fields for model output and diagnostics                          (unit )
! ------------------------------------------
     sod_t      ,   & ! solar downward radiation at top of atmosphere (     )
    asod_t            ! averaged solar downward radiation at top      (     )

! end of data_fields

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
    nold  ,       & ! corresponds to ntstep - 1
    nnow  ,       & ! corresponds to ntstep 
    nnew  ,       & ! corresponds to ntstep + 1

! 3. controlling the physics
! --------------------------
    itype_aerosol,& ! type of aerosol map internal/external
    l_cosmo_art,  & ! if .TRUE., run the COSMO_ART
    nincrad,      & ! time step increment for running the radiation
    hincrad,      & ! increment for running the radiation in hours
    nextrad,      & ! next step for running the radiation 
    hnextrad,     & ! next step for running the radiation in hours
    nradcoarse,   & ! radiation coarse-grid number of gpts per hor. direction !T.R.
    lradf_avg,    & ! switch for filtering of radiative increments !T.R.
    icldm_rad,    & ! mode of cloud representation in radiation  parametr.
    lradtopo,     & ! if .TRUE., calculate topographic correction of radiation
    nhori,        & ! number of sectors for the horizont array by the topographic
                    ! correction of the radiation
#ifdef TWOMOM_SB
    iradpar_cloud,& ! type of parameterization of radiative transfer parameters (ext., sing. alb., asym.)  
                    ! for grid- and subgrid scale clouds (cloud water, ice water)
                    !   1 = old method
                    !   2 = based on eff. radius (code from Elias Zubler, ETHZ)
#endif

    ! and for blocked physics
    nproma_rad, nlastproma_rad, nblock_rad, &
    nproma, nlastproma, nblock

USE data_runcontrol , ONLY :   &

! 5. additional control variables
! -------------------------------
    ltime,        & ! 
    l2tls     ,   & ! forecast with 1-TL integration scheme        
    lperi_x,      & ! if lartif_data=.TRUE.: periodic boundary conditions in x-dir.
                    !               .FALSE.:   or with Davies conditions
    lperi_y,      & ! if lartif_data=.TRUE.: periodic boundary conditions in y-dir.
                    !               .FALSE.:   or with Davies conditions
    l2dim,        & ! 2 dimensional runs
    lscm,         & ! SCM-run
    lprog_qi,     & ! if .TRUE., running with cloud ice

! 9. Variables for Ascii file handling, time measuring, ...
! ---------------------------------------------------------
    itype_calendar,&! for specifying the calendar used

! 12. controlling verbosity of debug output
! -----------------------------------------
    idbg_level,   & ! to control the verbosity of debug output
    ldebug_rad,   & ! if .TRUE., debug output for radiation
    lprintdeb_all   ! .TRUE.:  all tasks print debug output
                    ! .FALSE.: only task 0 prints debug output

! end of data_runcontrol 

!------------------------------------------------------------------------------

USE data_soil , ONLY :   &
    csalbw    , & !
    ctalb

! end of data_soil       

!------------------------------------------------------------------------------

USE data_parallel,      ONLY :  &
  nprocx,isubpos, &
  num_compute,   & ! number of compute PEs
  nboundlines,   & ! number of boundary lines of the domain for which
                   ! no forecast is computed = overlapping boundary
                   ! lines of the subdomains
  ldatatypes,    & ! if .TRUE.: use MPI-Datatypes for some communications
  ltime_barrier, & ! if .TRUE.: use additional barriers for determining the
                   ! load-imbalance
  ncomm_type,    & ! type of communication
  my_cart_id,    & ! rank of this subdomain in the cartesian communicator
  my_cart_neigh, & ! neighbors of this subdomain in the cartesian grid
  icomm_cart,    & ! communicator for the virtual cartesian topology
                   ! that can be used by MPI_WAIT to identify the send
  imp_reals,     & ! determines the correct REAL type used in the model
                   ! for MPI
  imp_integers,  & ! determines the correct INTEGER type used in the model
                   ! for MPI
  nexch_tag,     & ! tag to be used for MPI boundary exchange
                   !  (in calls to exchg_boundaries)
  sendbuf,       & ! sending buffer for boundary exchange:
                   ! 1-4 are used for sending, 5-8 are used for receiving
  isendbuflen      ! length of one column of sendbuf

!------------------------------------------------------------------------------

USE environment,              ONLY :  &
  exchg_boundaries,        & ! performs the boundary exchange between
                             ! neighboring processors
  model_abort
!------------------------------------------------------------------------------

USE radiation_data, ONLY : &
    aerdis,          & ! subroutine to compute vertical distribution of aerosols

!   idim_rad,        & ! ie-dimension of the coarser grid
!   istartrad,       & ! start- and end-indices for computing the radiation
!   iendrad,         & !   (when running on a coarser grid, the input values for
!   jstartrad,       & !    fesft are computed on all grid points, to compute an
!   jendrad,         & !    average input over several grid points)
!   iendparrad,      & ! end-index just for preparations
!   jendparrad,      & ! end-index just for preparations
!   istartradheat,   & !
!   iendradheat,     & !
!   jstartradheat,   & !
!   jendradheat,     & !

! Imported array data variables with intent (in) for radiation_init:
  jpgas,                   & ! Number of gases
  coali,                   & !
  cobti,                   & !
  pgas,                    & !
  tgas,                    & !
  jpsol   , & ! Number of solar spectral intervals
  jpther  , & ! Number of thermal spectral intervals
  jpspec  , & ! =jpsol+jpther (Total number of spectral intervals)

  solant  , & ! Fraction of solar energy at TOA contained in individual
              ! solar spectral intervals

  planck  , & ! coefficients for the description of the fraction of the total
              ! black body radiation contained in an thermal spectral interval
              ! as a function (2.order polynomial) of temperature
  zketypr , & ! e-type continuum-coefficient for all spectral intervals 
              ! (PA (H2O)**-2) at 296 K
  zketypa , & ! (r) following ROBERTS ET AL. 1976,(a) implicitly derived 
              ! from the AFGL spectral data
  ztetypr , & ! constant for the temperature dependancy of e-type absorption 
              ! for all intervals
  ztetypa , & ! (r) following ROBERTS ET AL. 1976,(a) implicitly derived 
              ! from the AFGL spectral data
  zteref  , & ! reference temperature

  grenze  , & ! Limits of spectral intervals (for information only)

  ! absorption properties of atmospheric gases
  ncgas   , & ! number of coefficients for each interval and gas (maximum=7)
  nfast   , & ! control variable for choice between ESFT/FESFT method in 
              ! each interval
  coai    , & ! weigthing coefficients
  cobi    , & ! absorption coefficients

! array data variables for opt_th and opt_so with intent (in):
  zaea, zaes, zaeg, zaef,         & !
  zlwe, zlwemn, zlwemx,           & !
  zlww, zlwg,                     & !
  ziwe, ziwemn, ziwemx,           & !
  ziww, ziwg,                     & !
  zrsc,                           & !

! for albedo
  rad_csalbw,                     & !

! 6a: parameters for the vertical distribution of aerosols
! --------------------------------------------------------

    pvdaes        , & ! normalized vertical distribution (sea)
    pvdael        , & ! normalized vertical distribution (land)
    pvdaeu        , & ! normalized vertical distribution (urban)
    pvdaed        , & ! normalized vertical distrubution (desert)
    ptrbga        , & ! b. optical depths divided by pressure (tropospheric)
    pvobga        , & ! b. optical depths divided by pressure (volcanic)
    pstbga        , & ! b. optical depths divided by pressure (stratospheric)
    paeops        , & ! total optical depths for vertical varying aerosols (sea)
    paeopl        , & ! total optical depths for vertical varying aerosols (land)
    paeopu        , & ! total optical depths for vertical varying aerosols (urban)
    paeopd        , & ! total optical depths for vertical varying aerosols (desert)
    ptrpt         , & ! temperature exponent for the stratosperic definition
    paeadk        , & ! constants for definition of the quantity of water 
    paeadm        , & ! vapour that will be adsorbed to the dry aerosols to
                      ! form moist aerosols

! 10. Organizational variables / types for a coarse radiation grid
! ----------------------------------------------------------------

  gp_radcoarse, gp_radcoarse_tot, gp_radcoarse_loc, ietot_rc, jetot_rc,          &
  ie_rc, je_rc, mind_ilon_rad, mind_jlat_rad

!------------------------------------------------------------------------------

USE src_tracer,         ONLY :  trcr_get, trcr_errorstr
USE data_tracer,        ONLY :  T_ERR_NOTFOUND

!------------------------------------------------------------------------------

USE src_block_fields,   ONLY :  mind_ilon, mind_jlat
USE radiation_utilities,ONLY :  co2_scenarios_and_absorbers, cloudiness_and_humidity, &
                                surface_albedo, compute_sunshine_conditions,          &
                                calc_rad_corrections
USE radiation_rg,       ONLY :  fesft, radiation_rg_wkarr_dealloc
USE utilities,          ONLY :  get_utc_date, dosleep
USE parallel_utilities, ONLY :  i_global, j_global
USE vgrid_refatm_utils, ONLY :  vcoord

!------------------------------------------------------------------------------

#ifdef COSMOART
USE art_papa,           ONLY :  calcjval

USE data_block_fields_art, ONLY:                &
                                Eup_b,          &
                                Edown_b,        &
                                Edir_b
USE data_cosmo_art,     ONLY:                   &
                                Eup           , & ! upward flux (3. band of GRAALS)            (W m-2)
                                Edown         , & ! downward flux (3. band of GRAALS)          (W m-2)
                                Edir          , & ! direkt intensity (3. band of GRAALS)       (W m-2)
                                mmy           , & ! cosinus of sun-zenith-angle in radiation   (1)
                                lgas

#ifdef TWOMOM_SB
USE src_cloud_opt_reff, ONLY :   &
   calc_cloud_opt,               &
   calc_optthick_cloud,          &
   zlwoda_so_prefac,             &
   zlwods_so_prefac,             &
   zlwb0_so,                     &
   zlwb_so_prefac,               &
   ziwoda_so_prefac,             &
   ziwods_so_prefac,             &
   ziwb0_so,                     &
   ziwb_so_prefac,               &
   zlwoda_th_prefac,             &
   zlwods_th_prefac,             &
   zlwb0_th,                     &
   zlwb_th,                      &
   ziwoda_th_prefac,             &
   ziwods_th_prefac,             &
   ziwb0_th
#endif
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================

! this is now for the blocked data structure
REAL (KIND=wp)       , ALLOCATABLE :: &
! temperature and pressure
  zti    (:,:,:)    , & ! Temperature at layer boundaries
  ztg_rc (:,:)      , & ! ground temperature on the coarse radiation grid
                        ! (needed to compute surface outgoing thermal radiation
                        !  with values from the coarse grid)
  zdpr   (:,:,:)    , & ! Pressure thickness of layers
  zapre  (:,:)      , & ! Surface pressure
! cloudiness_and_humidity
  zclc   (:,:,:)    , & ! Cloud cover in each layer
  zwv    (:,:,:)    , & ! Water vapour mixing ratio
  zsw    (:,:,:)    , & ! Saturation water vapour mixing ratio over water
  zclwc  (:,:,:)    , & ! liquid water mixing ratio
  zciwc  (:,:,:)    , & ! ice mixing ratio
! co2_scenarios_and_absorbers
  zduco2f(:,:,:)    , & ! CO2 content in layer
  zduo3f (:,:,:)    , & ! O3  content in layer
  zaeq1  (:,:,:)    , & ! Type1-Aerosole optical depth at 0.55 micrometer
  zaeq2  (:,:,:)    , & ! Type2  "
  zaeq3  (:,:,:)    , & ! Type3  "
  zaeq4  (:,:,:)    , & ! Type4  "
  zaeq5  (:,:,:)    , & ! Type5  "
! sunshine condition
  zsmu0_fg (:,:)    , & ! zenith angle of sun (retain values on full grid 
                        !                      for computations after fesft)
  zsmu0    (:)      , & ! zenith angle of sun (variable given to fesft)
  zskyview (:,:)    , & ! skyview: used as argument for SR fesft
! surface_albedo
  zalso    (:,:)    , & ! Solar surface albedo
  zalth    (:,:)        ! Thermal surface albedo 

! output of fesft for blocked data structure
REAL (KIND=wp)     , ALLOCATABLE :: &

  zflt   (:,:)       , & ! Thermal radiative flux at layer boundary
  zfls   (:,:)       , & ! Solar radiative flux at layer boundary
  zflsdir(:,:)       , & ! solar direct downward radiative flux at layer boundary
                         ! T.R: a 1d field would be sufficent since only (:,ke1) is used

  ! surface flux of photosynthetic active radiation and components
  zflpar_s     (:)   , & ! surface flux of photosynthetic active radiation
  zflpar_s_dir (:)   , & ! direct ("p"arallel ) component
  zflpar_s_difd(:)   , & ! diffuse downward component
  zflpar_s_difu(:)   , & ! diffuse upward   component

  ! corrected solar and thermal fluxes at layer boundary and components
  zfls_s   (:)       , & ! Corrected solar
  zflt_s   (:)       , & !         thermal
  zflsp    (:)       , & ! direct component of solar radiative flux
  zflsd    (:)       , & ! diffuse downward component of solar flux
  zflsu    (:)       , & ! diffuse upward   component of solar flux
  zfltd    (:)       , & ! diffuse downward component of thermal flux
  zfltu    (:)       , & ! diffuse upward   component of thermal flux

  ! and some work arrays (in (nproma_rad*nradcoarse,(ke),nradcoarse)
  zwork3d1 (:,:,:), &
  zwork3d2 (:,:,:), &
  zwork3d3 (:,:,:), &
  zwork2d1 (:,:),   &
  zwork2d2 (:,:),   &
  zwork2d3 (:,:),   &
  zwork2d4 (:,:),   &
  zwork2d5 (:,:),   &
  zwork2d6 (:,:),   &
  zwork2d7 (:,:),   &
  zwork2d8 (:,:),   &
  zwork3d  (:,:,:)   ! needed for averaging (in (ie,je,ke) if nradcoarse > 1

INTEGER, ALLOCATABLE :: &
  nzwork   (:,:)

! Tracer pointers
!----------------
REAL (KIND=wp),     POINTER ::        &
  qv_locrad  (:,:,:)=> NULL(),           & ! QV at ntl
  qc_locrad  (:,:,:)=> NULL(),           & ! QC at ntl
  qi_locrad  (:,:,:)=> NULL()              ! QI at ntl

REAL   (KIND=wp),   SAVE ::           &
      zsct_save, zdtzgl_save,         &
      zdeksin_save,zdekcos_save, zeit0

REAL    (KIND=wp)       , PARAMETER ::  &
   z_1d7200 = 1._wp/7200._wp ,&
   zcent  = 0.2500_wp, & ! centre weight in a nine point stencil !T.R.
   zside  = 0.1250_wp, & ! weight for side points !T.R.
   zedge  = 0.0625_wp, & ! weight for edge points !T.R.
   zepclc = MAX(1.0E-8_wp,rprecision), &
                         ! avoids cloud cover =1.0 and = 0.0 
   zeph2o = 1.0E-9_wp, & ! minimum value for specific humidity
   zepemu = 1.0E-9_wp, & ! avoids cosine of zenith angle = 0.0
   zclwcm = 1.0E-9_wp, & ! avoids cloud water content = 0.0
   rtod  = 57.2957795_wp ! conversion from radians to degrees

INTEGER , SAVE :: itaja_zsct_previous = 0

!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure in "Radiation" for initializing necessary variables
!------------------------------------------------------------------------------

SUBROUTINE radiation_init

!------------------------------------------------------------------------------
!
! Description:
! The module procedure radiation_init initializes the data necessary for the
! the radiation scheme. It provides the constant aerosol arrays 
! (aersea, aerlan, aerurb and aerdes) and processes the data variables 
! provided by module radiation_data for the calculation of the absorption 
! properties of atmospheric gases.
! The routine is called once at the beginning of a model run in order to
! convert the raw coffeicients from the data variables to those entities used
! in the radiation code.
!
! Method:
! - Computation of the inverse transformations (Legendre and Fourier) of a T10
!   representation of the four aerosol fields
! - scaling of absorption coefficients from the individual reference
!   temperature and pressure to unified conditions, i.e.
!   reference temperature = 273.15 K
!   reference pressure    = 1013.25 hPa
!
!------------------------------------------------------------------------------

! Subroutine arguments: None
! --------------------

! Local arrays and scalars:
! -------------------------
  INTEGER                   ::  &
    i, j, k  , jg  , js , jc ,  & ! loop indices
    jzj, jzm1, jzm2, jzm, jzn,  & ! indices for Legendre coefficients
    imn, imnc, imns, jmm, jnn,  & !
    ist, ip, iplocend, jp, ib

  REAL (KIND=wp)            ::  &
                               ! arrays for the T10 distrubution of
    zaesc(66) , zaess (55) , & ! sea    type aerosols                     
    zaelc(66) , zaels (55) , & ! land   type aerosols
    zaeuc(66) , zaeus (55) , & ! urban  type aerosols
    zaedc(66) , zaeds (55) , & ! desert type aerosols
    zfaes(21) , zfael (21) , & ! coefficients for spectral
    zfaeu(21) , zfaed (21) , & ! expansion
    zalp (66) ,              & !
    zsign(ke1),              & !
    zsinphi   , zcosphi    , & !
    zm, z2m, zre1, ze1, ze2, & !
    zf1m, zf2m, zn, zn2,     & !
    zsin1, zsin2, zsin3, zsin4, zsin5,  & ! 
    zsin6, zsin7, zsin8, zsin9, zsin10, & !
    zcos1, zcos2, zcos3, zcos4, zcos5,  & ! 
    zcos6, zcos7, zcos8, zcos9, zcos10    !

  REAL (KIND=wp)            ::  &
    zdzwb
 
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 0: Data for the Fourier coefficients of the four aerosol types          
!------------------------------------------------------------------------------

  DATA zaesc/                                                                                 &
     +.6688E+00_wp,-.1172E+00_wp,-.1013E+00_wp,+.1636E-01_wp,-.3699E-01_wp,+.1775E-01_wp,     &
                   -.9635E-02_wp,+.1290E-02_wp,+.4681E-04_wp,-.9106E-04_wp,+.9355E-04_wp,     &
     -.7076E-01_wp,-.1782E-01_wp,+.1856E-01_wp,+.1372E-01_wp,+.8210E-04_wp,+.2149E-02_wp,     &
                                 +.4856E-03_wp,+.2231E-03_wp,+.1824E-03_wp,+.1960E-05_wp,     &
     +.2057E-01_wp,+.2703E-01_wp,+.2424E-01_wp,+.9716E-02_wp,+.1312E-02_wp,-.8846E-03_wp,     &
                                               -.3347E-03_wp,+.6231E-04_wp,+.6397E-04_wp,     &
     -.3341E-02_wp,-.1295E-01_wp,-.4598E-02_wp,+.3242E-03_wp,+.8122E-03_wp,-.2975E-03_wp,     &
                                                             -.7757E-04_wp,+.7793E-04_wp,     &
     +.4455E-02_wp,-.1584E-01_wp,-.2551E-02_wp,+.1174E-02_wp,+.1335E-04_wp,+.5112E-04_wp,     &
                                                                           +.5605E-04_wp,     &
     +.7412E-04_wp,+.1857E-02_wp,-.1917E-03_wp,+.4460E-03_wp,+.1767E-04_wp,-.5281E-04_wp,     &
     -.5043E-03_wp,+.2467E-03_wp,-.2497E-03_wp,-.2377E-04_wp,-.3954E-04_wp,                   &
     +.2666E-03_wp,-.8186E-03_wp,-.1441E-03_wp,-.1904E-04_wp,                                 &
     +.3337E-03_wp,-.1696E-03_wp,-.2503E-04_wp,                                               &
     +.1239E-03_wp,-.9983E-04_wp,                                                             &
     -.5283E-04_wp                                                                            /

  DATA zaess/                                                                                 &
     -.3374E-01_wp,-.3247E-01_wp,-.1012E-01_wp,+.6002E-02_wp,+.5190E-02_wp,+.7784E-03_wp,     &
                                 -.1090E-02_wp,+.3294E-03_wp,+.1719E-03_wp,-.5866E-05_wp,     &
     -.4124E-03_wp,-.3742E-01_wp,-.5054E-02_wp,+.3430E-02_wp,+.5513E-03_wp,-.6235E-03_wp,     &
                                               +.2892E-03_wp,-.9730E-04_wp,+.7078E-04_wp,     &
     -.3300E-01_wp,+.5104E-03_wp,-.2156E-02_wp,-.3194E-02_wp,-.5079E-03_wp,-.5517E-03_wp,     &
                                                             +.4632E-04_wp,+.5369E-04_wp,     &
     -.2731E-01_wp,+.5126E-02_wp,+.2241E-02_wp,-.5789E-03_wp,-.3048E-03_wp,-.1774E-03_wp,     &
                                                                           +.1946E-05_wp,     &
     -.8247E-02_wp,+.2338E-02_wp,+.1021E-02_wp,+.1575E-04_wp,+.2612E-05_wp,+.1995E-04_wp,     &
     -.1319E-02_wp,+.1384E-02_wp,-.4159E-03_wp,-.2337E-03_wp,+.5764E-04_wp,                   &
     +.1495E-02_wp,-.3727E-03_wp,+.6075E-04_wp,-.4642E-04_wp,                                 &
     +.5368E-03_wp,-.7619E-04_wp,+.3774E-04_wp,                                               &
     +.1206E-03_wp,-.4104E-06_wp,                                                             &
     +.2158E-04_wp                                                                            /

  DATA zaelc/                                                                                 &
     +.1542E+00_wp,+.8245E-01_wp,-.1879E-03_wp,+.4864E-02_wp,-.5527E-02_wp,-.7966E-02_wp,     &
                   -.2683E-02_wp,-.2011E-02_wp,-.8889E-03_wp,-.1058E-03_wp,-.1614E-04_wp,     &
     +.4206E-01_wp,+.1912E-01_wp,-.9476E-02_wp,-.6780E-02_wp,+.1767E-03_wp,-.5422E-03_wp,     &
                                 -.7753E-03_wp,-.2106E-03_wp,-.9870E-04_wp,-.1721E-04_wp,     &
     -.9536E-02_wp,-.9580E-02_wp,-.1050E-01_wp,-.5747E-02_wp,-.1282E-02_wp,+.2248E-03_wp,     &
                                               +.1694E-03_wp,-.4782E-04_wp,-.2441E-04_wp,     &
     +.5781E-03_wp,+.6212E-02_wp,+.1921E-02_wp,-.1102E-02_wp,-.8145E-03_wp,+.2497E-03_wp,     &
                                                             +.1539E-03_wp,-.2538E-04_wp,     &
     -.3993E-02_wp,+.9777E-02_wp,+.4837E-03_wp,-.1304E-02_wp,+.2417E-04_wp,-.1370E-04_wp,     &
                                                                           -.3731E-05_wp,     &
     +.1922E-02_wp,-.5167E-03_wp,+.4295E-03_wp,-.1888E-03_wp,+.2427E-04_wp,+.4012E-04_wp,     &
     +.1529E-02_wp,-.2120E-03_wp,+.8166E-04_wp,+.2579E-04_wp,+.3488E-04_wp,                   &
     +.2140E-03_wp,+.2274E-03_wp,-.3447E-05_wp,-.1075E-04_wp,                                 &
     -.1018E-03_wp,+.2864E-04_wp,+.3442E-04_wp,                                               &
     -.1002E-03_wp,+.7117E-04_wp,                                                             &
     +.2045E-04_wp                                                                            /

  DATA zaels/                                                                                 &
     +.1637E-01_wp,+.1935E-01_wp,+.1080E-01_wp,+.2784E-02_wp,+.1606E-03_wp,+.1860E-02_wp,     &
                                 +.1263E-02_wp,-.2707E-03_wp,-.2290E-03_wp,-.9761E-05_wp,     &
     -.7317E-02_wp,+.2465E-01_wp,+.6799E-02_wp,-.1913E-02_wp,+.1382E-02_wp,+.6691E-03_wp,     &
                                               +.1414E-03_wp,+.3527E-04_wp,-.5210E-04_wp,     &
     +.1873E-01_wp,+.2977E-02_wp,+.4650E-02_wp,+.2509E-02_wp,+.3680E-03_wp,+.1481E-03_wp,     &
                                                             -.6594E-04_wp,-.5634E-04_wp,     &
     +.1592E-01_wp,-.1875E-02_wp,-.1093E-02_wp,+.3022E-03_wp,+.2625E-03_wp,+.3252E-04_wp,     &
                                                                           -.3803E-04_wp,     &
     +.4218E-02_wp,-.1843E-02_wp,-.1351E-02_wp,-.2952E-03_wp,-.8171E-05_wp,-.1473E-04_wp,     &
     +.9076E-03_wp,-.1057E-02_wp,+.2676E-03_wp,+.1307E-03_wp,-.3628E-04_wp,                   &
     -.9158E-03_wp,+.4335E-03_wp,+.2927E-04_wp,+.6602E-04_wp,                                 &
     -.3570E-03_wp,+.5760E-04_wp,-.3465E-04_wp,                                               &
     -.8535E-04_wp,-.2011E-04_wp,                                                             &
     +.6612E-06_wp                                                                            /

  DATA zaeuc/                                                                                 &
     +.8005E-01_wp,+.7095E-01_wp,+.2014E-01_wp,-.1412E-01_wp,-.2425E-01_wp,-.1332E-01_wp,     &
                   -.2904E-02_wp,+.5068E-03_wp,+.9369E-03_wp,+.4114E-03_wp,+.7549E-04_wp,     &
     +.1922E-01_wp,+.2534E-01_wp,+.2088E-01_wp,+.1064E-01_wp,+.1063E-02_wp,-.2526E-02_wp,     &
                                 -.2091E-02_wp,-.9660E-03_wp,-.2030E-03_wp,+.3865E-04_wp,     &
     -.9900E-02_wp,-.5964E-02_wp,+.2223E-02_wp,+.4941E-02_wp,+.3277E-02_wp,+.1038E-02_wp,     &
                                               -.1480E-03_wp,-.2844E-03_wp,-.1208E-03_wp,     &
     +.3999E-02_wp,+.6282E-02_wp,+.2813E-02_wp,+.1475E-02_wp,+.4571E-03_wp,-.1349E-03_wp,     &
                                                             -.9011E-04_wp,-.1936E-04_wp,     &
     +.1994E-02_wp,+.3540E-02_wp,+.8837E-03_wp,+.1992E-03_wp,+.3092E-04_wp,-.7979E-04_wp,     &
                                                                           -.2664E-04_wp,     &
     -.5006E-04_wp,+.6447E-03_wp,+.5550E-03_wp,+.1197E-03_wp,+.6657E-04_wp,+.1488E-04_wp,     &
     -.9141E-04_wp,-.2896E-03_wp,-.1561E-03_wp,-.6524E-04_wp,-.1559E-04_wp,                   &
     -.1082E-03_wp,-.4126E-03_wp,-.1732E-03_wp,-.8286E-04_wp,                                 &
     -.1993E-04_wp,+.3850E-04_wp,+.2870E-04_wp,                                               &
     +.4493E-04_wp,+.4721E-04_wp,                                                             &
     +.1338E-04_wp                                                                            /

  DATA zaeus/                                                                                 &
     +.6646E-02_wp,+.8373E-02_wp,+.5463E-02_wp,+.4554E-02_wp,+.3301E-02_wp,+.5725E-03_wp,     &
                                 -.7482E-03_wp,-.6222E-03_wp,-.2603E-03_wp,-.5127E-04_wp,     &
     -.3849E-04_wp,+.9741E-02_wp,+.8190E-02_wp,+.5712E-02_wp,+.3039E-02_wp,+.5290E-03_wp,     &
                                               -.2044E-03_wp,-.2309E-03_wp,-.1160E-03_wp,     &
     +.9160E-02_wp,+.1286E-01_wp,+.1170E-01_wp,+.5491E-02_wp,+.1393E-02_wp,-.6288E-04_wp,     &
                                                             -.2715E-03_wp,-.1047E-03_wp,     &
     +.4873E-02_wp,+.3545E-02_wp,+.3069E-02_wp,+.1819E-02_wp,+.6947E-03_wp,+.1416E-03_wp,     &
                                                                           -.1538E-04_wp,     &
     -.4351E-03_wp,-.1907E-02_wp,-.5774E-03_wp,-.2247E-03_wp,+.5345E-04_wp,+.9052E-04_wp,     &
     -.3972E-04_wp,-.9665E-04_wp,+.7912E-04_wp,-.1094E-04_wp,-.6776E-05_wp,                   &
     +.2724E-03_wp,+.1973E-03_wp,+.6837E-04_wp,+.4313E-04_wp,                                 &
     -.7174E-05_wp,+.8527E-05_wp,-.2160E-05_wp,                                               &
     -.7852E-04_wp,+.3453E-06_wp,                                                             &
     -.2402E-05_wp                                                                            /

  DATA zaedc/                                                                                 &
     +.2840E-01_wp,+.1775E-01_wp,-.1069E-01_wp,-.1553E-01_wp,-.3299E-02_wp,+.3583E-02_wp,     &
                   +.2274E-02_wp,+.5767E-04_wp,-.3678E-03_wp,-.1050E-03_wp,+.2133E-04_wp,     &
     +.2326E-01_wp,+.1566E-01_wp,-.3130E-02_wp,-.8253E-02_wp,-.2615E-02_wp,+.1247E-02_wp,     &
                                 +.1059E-02_wp,+.1196E-03_wp,-.1303E-03_wp,-.5094E-04_wp,     &
     +.1185E-01_wp,+.7238E-02_wp,-.1562E-02_wp,-.3665E-02_wp,-.1182E-02_wp,+.4678E-03_wp,     &
                                               +.4448E-03_wp,+.8307E-04_wp,-.3468E-04_wp,     &
     +.5273E-02_wp,+.3037E-02_wp,-.4014E-03_wp,-.1202E-02_wp,-.4647E-03_wp,+.5148E-04_wp,     &
                                                             +.1014E-03_wp,+.2996E-04_wp,     &
     +.2505E-02_wp,+.1495E-02_wp,+.2438E-03_wp,-.1223E-03_wp,-.7669E-04_wp,-.1638E-04_wp,     &
                                                                           +.1869E-05_wp,     &
     +.1094E-02_wp,+.6131E-03_wp,+.1508E-03_wp,+.1765E-04_wp,+.1360E-05_wp,-.7998E-06_wp,     &
     +.4475E-03_wp,+.2737E-03_wp,+.6430E-04_wp,-.6759E-05_wp,-.6761E-05_wp,                   &
     +.1992E-03_wp,+.1531E-03_wp,+.4828E-04_wp,+.5103E-06_wp,                                 &
     +.7454E-04_wp,+.5917E-04_wp,+.2152E-04_wp,                                               &
     +.9300E-05_wp,+.9790E-05_wp,                                                             &
     -.8853E-05_wp                                                                            /

  DATA zaeds/                                                                                 &
     +.9815E-02_wp,+.8436E-02_wp,+.1087E-02_wp,-.2717E-02_wp,-.1755E-02_wp,-.1559E-03_wp,     &
                                 +.2367E-03_wp,+.8808E-04_wp,+.2001E-05_wp,-.1244E-05_wp,     &
     +.1041E-01_wp,+.8039E-02_wp,+.1005E-02_wp,-.1981E-02_wp,-.1090E-02_wp,+.1595E-05_wp,     &
                                               +.1787E-03_wp,+.4644E-04_wp,-.1052E-04_wp,     &
     +.6593E-02_wp,+.3983E-02_wp,-.1527E-03_wp,-.1235E-02_wp,-.5078E-03_wp,+.3649E-04_wp,     &
                                                             +.1005E-03_wp,+.3182E-04_wp,     &
     +.3225E-02_wp,+.1672E-02_wp,-.7752E-04_wp,-.4312E-03_wp,-.1872E-03_wp,-.1666E-04_wp,     &
                                                                           +.1872E-04_wp,     &
     +.1133E-02_wp,+.5643E-03_wp,+.7747E-04_wp,-.2980E-04_wp,-.2092E-04_wp,-.8590E-05_wp,     &
     +.2988E-03_wp,+.6714E-04_wp,-.6249E-05_wp,+.1052E-04_wp,+.8790E-05_wp,                   &
     +.1569E-03_wp,-.1175E-04_wp,-.3033E-04_wp,-.9777E-06_wp,                                 &
     +.1101E-03_wp,+.6827E-05_wp,-.1023E-04_wp,                                               &
     +.4231E-04_wp,+.4905E-05_wp,                                                             &
     +.6229E-05_wp                                                                            /


!------------------------------------------------------------------------------
! Begin Subroutine radiation_init             
!------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------
! Section 1: Compute the start- and end-indices, Initializations
!------------------------------------------------------------------------------

IF (nradcoarse == 1) THEN
  ! use the full grid
  nblock_rad     = nblock
  nproma_rad     = nproma
  nlastproma_rad = nlastproma

  ! Allocate and initialize mind_ilon_rad, mind_jlat_rad: the correspondance 
  ! between the blocked (coarse) radiation grid and the COSMO grid
  ALLOCATE (mind_ilon_rad(nproma_rad, nradcoarse, nblock_rad),      &
            mind_jlat_rad(nproma_rad, nradcoarse, nblock_rad),      &
                                                         STAT=ist)
  IF (ist /= 0) THEN
    CALL model_abort(my_cart_id, 2815, 'Error allocating mind_ilon_rad/jlat_rad_rad',  &
                                       'radiation_init')
  ENDIF

  ! The values of mind_ilon/jlat_rad are in this case the same as mind_ilon/jlat

  DO ib = 1, nblock_rad
    jp = 1
    IF (ib == nblock_rad) THEN
      iplocend = nlastproma_rad
    ELSE
      iplocend = nproma_rad
    ENDIF

    DO ip = 1, iplocend
      mind_ilon_rad (ip,jp,ib) = mind_ilon(ip,ib)
      mind_jlat_rad (ip,jp,ib) = mind_jlat(ip,ib)
    ENDDO
  ENDDO

  ! The ELSE-case (nradcoarse > 1) is considered in radiation_radcoarse_init
  ! and the above values are overwritten then
ENDIF

!------------------------------------------------------------------------------
! Section 2: Calculation of the inverse Legendre and Fourier transformation  
!------------------------------------------------------------------------------

! loops in i and j over model domain

IF (itype_aerosol == 1) THEN

  ! these variables are in double precision in any case
  zaea=RESHAPE((/0.0477_dp, 0.0875_dp,  0.1198_dp, 0.0458_dp, &
                 0.0387_dp, 0.0439_dp,  0.0599_dp, 0.0396_dp, &
                 0.0381_dp, 0.0129_dp,  0.0130_dp, 0.1304_dp, &
                 0.1757_dp, 0.0949_dp,  0.0653_dp, 0.0795_dp, &
                 0.0962_dp, 0.2046_dp,  0.4116_dp, 0.0169_dp, &
                 0.0204_dp, 0.0263_dp,  0.0348_dp, 0.0361_dp, &
                 0.0030_dp, 0.0271_dp,  0.0613_dp, 0.0118_dp, &
                 0.0160_dp, 0.0231_dp,  0.0287_dp, 0.0127_dp, &
                 0.0103_dp, 0.000016_dp,0.0000_dp, 0.0087_dp, &
                 0.0238_dp, 0.0511_dp,  0.0734_dp, 0.0809_dp/),(/8,5/))

  zaes=RESHAPE((/0.1407_dp, 0.4256_dp,  1.0066_dp, 0.0279_dp, &
                 0.0391_dp, 0.0445_dp,  0.0485_dp, 0.0362_dp, &
                 0.6746_dp, 0.8761_dp,  1.0139_dp, 0.0443_dp, &
                 0.0624_dp, 0.0921_dp,  0.1491_dp, 0.2327_dp, &
                 0.0605_dp, 0.2761_dp,  0.7449_dp, 0.0023_dp, &
                 0.0034_dp, 0.0051_dp,  0.0065_dp, 0.0045_dp, &
                 0.0284_dp, 0.5524_dp,  0.9683_dp, 0.0001_dp, &
                 0.0004_dp, 0.0024_dp,  0.0049_dp, 0.0030_dp, &
                 0.0467_dp, 0.3854_dp,  1.1008_dp, 0.0000_dp, &
                 0.00005_dp,0.0004_dp,  0.0006_dp, 0.0006_dp/),(/8,5/))

  zaeg=RESHAPE((/0.6989_dp, 0.6329_dp,  0.6418_dp, 0.6243_dp, &
                 0.7299_dp, 0.7430_dp,  0.7086_dp, 0.8569_dp, &
                 0.7833_dp, 0.7575_dp,  0.7456_dp, 0.4997_dp, &
                 0.6130_dp, 0.7440_dp,  0.7426_dp, 0.7590_dp, &
                 0.5753_dp, 0.5867_dp,  0.5957_dp, 0.6027_dp, &
                 0.6766_dp, 0.6117_dp,  0.5439_dp, 0.6905_dp, &
                 0.5170_dp, 0.6674_dp,  0.7004_dp, 0.0340_dp, &
                 0.0570_dp, 0.1289_dp,  0.1597_dp, 0.1906_dp, &
                 0.3751_dp, 0.6353_dp,  0.7259_dp, 0.0037_dp, &
                 0.0083_dp, 0.0177_dp,  0.0201_dp, 0.0332_dp/),(/8,5/))

  zaef(:,:)= 0.0_dp

  DO j = 1, je
    DO i = 1, ie
 
     ! Calculation of the values zalp for the sine of latitude (zsinphi) of the
     ! normalized Legendre associated functions. The limit wave number is 10.

     IF ((.NOT. lscm) .AND. lperi_x) THEN
       ! In case of periodic BCs, set a constant reference
       ! point for the aerosol distribution to make it equal everywhere
       ! in the domain. This reference point is chosen to be the reference point
       ! of the model domain as determined by pollon, pollat:
       zsinphi  = SIN (degrad*(90.0_wp-ABS(pollat)))
     ELSE
       zsinphi  = SIN(rlat(i,j) )
     ENDIF
     zcosphi  = SQRT(1.0_wp-zsinphi**2)
     jzj      = 2
     zf1m     = SQRT(3.0_wp)
     zalp (1) = 1.0_wp
     zalp (2) = zf1m*zsinphi
     wave_number_loop : DO jzm1 = 1, 11 
       jzm  = jzm1-1
       zm   = REAL (jzm, wp)
       z2m  = zm + zm
       zre1 = SQRT(z2m+3.0_wp)
       ze1  = 1.0_wp/zre1
       IF (jzm.NE.0) THEN     
          zf2m      = zf1m*zcosphi/SQRT(z2m)
          zf1m      = zf2m*zre1
          jzj       = jzj + 1
          zalp(jzj) = zf2m
          IF(jzm ==10) CYCLE wave_number_loop
          jzj       = jzj + 1
          zalp(jzj) = zf1m*zsinphi
          IF(jzm1==10) CYCLE wave_number_loop
       ENDIF  
       jzm2 = jzm+2
       DO jzn = jzm2, 10
          zn        = REAL (jzn, wp)
          zn2       = zn**2
          ze2       = SQRT( (4.0_wp*zn2-1.0_wp)/(zn2-zm**2) )
          jzj       = jzj+1
          zalp(jzj) = ze2*(zsinphi*zalp(jzj-1)-ze1*zalp(jzj-2))
          ze1       = 1.0_wp/ze2
       ENDDO  
     ENDDO wave_number_loop

     ! Legendre transform of aerosols
     zfaes(:) = 0.0_wp
     zfael(:) = 0.0_wp
     zfaeu(:) = 0.0_wp
     zfaed(:) = 0.0_wp

     imn  = 0
     imnc = 0
     imns = 0

     DO jmm = 1, 11
        imn  = imn  + 1
        DO jnn = jmm, 11
           imnc       = imnc + 1
           zfaes(imn) = zfaes(imn)+zalp(imnc)*zaesc(imnc)
           zfael(imn) = zfael(imn)+zalp(imnc)*zaelc(imnc)
           zfaeu(imn) = zfaeu(imn)+zalp(imnc)*zaeuc(imnc)
           zfaed(imn) = zfaed(imn)+zalp(imnc)*zaedc(imnc)
        ENDDO    
        IF(jmm.NE.1) THEN
           imn  = imn+1
           DO jnn = jmm, 11
              imns       = imns + 1
              zfaes(imn) = zfaes(imn)+zalp(imns+11)*zaess(imns)
              zfael(imn) = zfael(imn)+zalp(imns+11)*zaels(imns)
              zfaeu(imn) = zfaeu(imn)+zalp(imns+11)*zaeus(imns)
              zfaed(imn) = zfaed(imn)+zalp(imns+11)*zaeds(imns)
           ENDDO  
        ENDIF
     ENDDO   
 
     ! Inverse Fourier transformation   
     IF ((.NOT. lscm) .AND. (lperi_y .OR. l2dim)) THEN
       ! In case of periodic BCs, set a constant reference
       ! point for the aerosol distribution to make it equal everywhere
       ! in the domain. This reference point is chosen to be the reference point
       ! of the model domain as determined by pollon, pollat:
       zcos1   = COS(degrad*(pollon-SIGN(1.0_wp,pollon)*180.0_wp))
       zsin1   = SIN(degrad*(pollon-SIGN(1.0_wp,pollon)*180.0_wp))
     ELSE
       zcos1   = COS(rlon(i,j) )
       zsin1   = SIN(rlon(i,j) )
     ENDIF
     zcos2   = zcos1*zcos1 - zsin1*zsin1
     zsin2   = zsin1*zcos1 + zcos1*zsin1
     zcos3   = zcos2*zcos1 - zsin2*zsin1
     zsin3   = zsin2*zcos1 + zcos2*zsin1
     zcos4   = zcos3*zcos1 - zsin3*zsin1
     zsin4   = zsin3*zcos1 + zcos3*zsin1
     zcos5   = zcos4*zcos1 - zsin4*zsin1
     zsin5   = zsin4*zcos1 + zcos4*zsin1
     zcos6   = zcos5*zcos1 - zsin5*zsin1
     zsin6   = zsin5*zcos1 + zcos5*zsin1
     zcos7   = zcos6*zcos1 - zsin6*zsin1
     zsin7   = zsin6*zcos1 + zcos6*zsin1
     zcos8   = zcos7*zcos1 - zsin7*zsin1
     zsin8   = zsin7*zcos1 + zcos7*zsin1
     zcos9   = zcos8*zcos1 - zsin8*zsin1
     zsin9   = zsin8*zcos1 + zcos8*zsin1
     zcos10  = zcos9*zcos1 - zsin9*zsin1
     zsin10  = zsin9*zcos1 + zcos9*zsin1
 
     aersea(i,j) = zfaes(1) + 2.0_wp* ( zfaes(2 ) * zcos1 + zfaes(3 ) * zsin1    &
                                  + zfaes(4 ) * zcos2 + zfaes(5 ) * zsin2    &
                                  + zfaes(6 ) * zcos3 + zfaes(7 ) * zsin3    &
                                  + zfaes(8 ) * zcos4 + zfaes(9 ) * zsin4    &
                                  + zfaes(10) * zcos5 + zfaes(11) * zsin5    &
                                  + zfaes(12) * zcos6 + zfaes(13) * zsin6    &
                                  + zfaes(14) * zcos7 + zfaes(15) * zsin7    &
                                  + zfaes(16) * zcos8 + zfaes(17) * zsin8    &
                                  + zfaes(18) * zcos9 + zfaes(19) * zsin9    &
                                  + zfaes(20) * zcos10+ zfaes(21) * zsin10 )

     aerlan(i,j) = zfael(1) + 2.0_wp* ( zfael(2 ) * zcos1 + zfael(3 ) * zsin1    &
                                  + zfael(4 ) * zcos2 + zfael(5 ) * zsin2    &
                                  + zfael(6 ) * zcos3 + zfael(7 ) * zsin3    &
                                  + zfael(8 ) * zcos4 + zfael(9 ) * zsin4    &
                                  + zfael(10) * zcos5 + zfael(11) * zsin5    &
                                  + zfael(12) * zcos6 + zfael(13) * zsin6    &
                                  + zfael(14) * zcos7 + zfael(15) * zsin7    &
                                  + zfael(16) * zcos8 + zfael(17) * zsin8    &
                                  + zfael(18) * zcos9 + zfael(19) * zsin9    &
                                  + zfael(20) * zcos10+ zfael(21) * zsin10 )
     aerurb(i,j) = zfaeu(1) + 2.0_wp* ( zfaeu(2 ) * zcos1 + zfaeu(3 ) * zsin1    &
                                  + zfaeu(4 ) * zcos2 + zfaeu(5 ) * zsin2    &
                                  + zfaeu(6 ) * zcos3 + zfaeu(7 ) * zsin3    &
                                  + zfaeu(8 ) * zcos4 + zfaeu(9 ) * zsin4    &
                                  + zfaeu(10) * zcos5 + zfaeu(11) * zsin5    &
                                  + zfaeu(12) * zcos6 + zfaeu(13) * zsin6    &
                                  + zfaeu(14) * zcos7 + zfaeu(15) * zsin7    &
                                  + zfaeu(16) * zcos8 + zfaeu(17) * zsin8    &
                                  + zfaeu(18) * zcos9 + zfaeu(19) * zsin9    &
                                  + zfaeu(20) * zcos10+ zfaeu(21) * zsin10 )
     aerdes(i,j) = zfaed(1) + 2.0_wp* ( zfaed(2 ) * zcos1 + zfaed(3 ) * zsin1    &
                                  + zfaed(4 ) * zcos2 + zfaed(5 ) * zsin2    &
                                  + zfaed(6 ) * zcos3 + zfaed(7 ) * zsin3    &
                                  + zfaed(8 ) * zcos4 + zfaed(9 ) * zsin4    &
                                  + zfaed(10) * zcos5 + zfaed(11) * zsin5    &
                                  + zfaed(12) * zcos6 + zfaed(13) * zsin6    &
                                  + zfaed(14) * zcos7 + zfaed(15) * zsin7    &
                                  + zfaed(16) * zcos8 + zfaed(17) * zsin8    &
                                  + zfaed(18) * zcos9 + zfaed(19) * zsin9    &
                                  + zfaed(20) * zcos10+ zfaed(21) * zsin10 )

     aersea(i,j) = MAX( 0.0_wp, MIN( 1.0_wp, aersea(i,j) ) )
     aerlan(i,j) = MAX( 0.0_wp, MIN( 1.0_wp, aerlan(i,j) ) )
     aerurb(i,j) = MAX( 0.0_wp, MIN( 1.0_wp, aerurb(i,j) ) )
     aerdes(i,j) = MAX( 0.0_wp, MIN( 1.0_wp, aerdes(i,j) ) )
 
    ! end of loops over model domain      
    ENDDO 
  ENDDO       
ENDIF ! itype_aerosol = 1
 
IF (itype_aerosol == 2) THEN
  ! these variables are in double precision in any case
  zaea=RESHAPE((/0.0345_dp, 0.0511_dp,  0.0847_dp, 0.0336_dp, &
                 0.0499_dp, 0.0364_dp,  0.0382_dp, 0.0260_dp, &
                 0.0457_dp, 0.0018_dp,  0.0015_dp, 0.1361_dp, &
                 0.2346_dp, 0.1177_dp,  0.0684_dp, 0.0808_dp, &
                 0.0707_dp, 0.0689_dp,  0.1557_dp, 0.1258_dp, &
                 0.1588_dp, 0.1973_dp,  0.2766_dp, 0.1134_dp, &
                 0.0597_dp, 0.1077_dp,  0.2095_dp, 0.0299_dp, &
                 0.0456_dp, 0.0358_dp,  0.0377_dp, 0.0304_dp, &
                 0.0103_dp, 0.000016_dp,0.0000_dp, 0.0087_dp, &
                 0.0238_dp, 0.0511_dp,  0.0734_dp, 0.0809_dp/),(/8,5/))

  zaes=RESHAPE((/0.1030_dp, 0.3977_dp,  1.0680_dp, 0.0084_dp, &
                 0.0142_dp, 0.0191_dp,  0.0234_dp, 0.0140_dp, &
                 0.7894_dp, 0.9734_dp,  1.0110_dp, 0.0307_dp, &
                 0.0531_dp, 0.0546_dp,  0.0839_dp, 0.2142_dp, &
                 0.7157_dp, 0.8698_dp,  0.8604_dp, 0.0645_dp, &
                 0.0781_dp, 0.1256_dp,  0.2317_dp, 0.1409_dp, &
                 0.0859_dp, 0.3442_dp,  0.9496_dp, 0.0067_dp, &
                 0.0113_dp, 0.0153_dp,  0.0187_dp, 0.0113_dp, &
                 0.0467_dp, 0.3854_dp,  1.1008_dp, 0.0000_dp, &
                 0.00005_dp,0.0004_dp,  0.0006_dp, 0.0006_dp/),(/8,5/))

  zaeg=RESHAPE((/0.6562_dp, 0.6614_dp,  0.7109_dp, 0.5043_dp, &
                 0.6486_dp, 0.6814_dp,  0.6489_dp, 0.7799_dp, &
                 0.8105_dp, 0.7906_dp,  0.7947_dp, 0.4374_dp, &
                 0.5203_dp, 0.7076_dp,  0.7246_dp, 0.7535_dp, &
                 0.6932_dp, 0.6962_dp,  0.7402_dp, 0.4029_dp, &
                 0.5587_dp, 0.5618_dp,  0.4520_dp, 0.7120_dp, &
                 0.6462_dp, 0.6510_dp,  0.6955_dp, 0.5041_dp, &
                 0.6482_dp, 0.6805_dp,  0.6477_dp, 0.7753_dp, &
                 0.3751_dp, 0.6353_dp,  0.7259_dp, 0.0037_dp, &
                 0.0083_dp, 0.0177_dp,  0.0201_dp, 0.0332_dp/),(/8,5/))

  zaef(:,:) = 0.0_dp
ENDIF ! itype_aerosol = 2

!------------------------------------------------------------------------------
! Section 3: Data for radiative transfer calculations
!------------------------------------------------------------------------------

IF (my_cart_id == 0) THEN
  PRINT *,'*****************************************************'
  PRINT *,'*   Radiative transfer calculations employ data     *'
  PRINT *,'*        provided in routine rad_aibi               *'
  PRINT *,'*****************************************************'
ENDIF

! Include reference pressure and temperature in *cobi*
! calculations in dp!
  DO jg = 1, jpgas
    DO js = 1, jpspec
      DO jc = 1, ncgas(js,jg)
         cobi(jc,js,jg) = cobi(jc,js,jg) * (1.0_dp/pgas(js,jg))**coali(jc,js,jg) &
                                         * (       tgas(js,jg))**cobti(jc,js,jg)
      ENDDO
    ENDDO
  ENDDO

! security settings, if a gas shall not be considered in an intervall
! where the esft will be used.

  DO jg = 1, jpgas
    DO js = 1, jpspec
      IF ( nfast(js) == 0 .AND. ncgas(js,jg) == 0 ) THEN
         ncgas     (js,jg) = 1
         coai    (1,js,jg) = 1.00_dp
         cobi    (1,js,jg) = 0.00_dp
         coali   (1,js,jg) = 1.00_dp
         cobti   (1,js,jg) = 1.00_dp
      END IF
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section 4: Precalculation of albedo of soil type as function of soil water
!            content and depth of upper soil layer
!------------------------------------------------------------------------------

  ! Albedo of soil type (without vegetation) as function of soil water content
  ! (in mH2O) and depth of the upper soil layer      

  zdzwb = 0.0_wp 
  zdzwb = 1.0_wp / (2.0_wp * czmls(1))
  DO ist = 1, 10
    rad_csalbw(ist) = csalbw(ist) * zdzwb
  ENDDO

!------------------------------------------------------------------------------
! Section 5: Initialize background aerosol (aerdis)
!------------------------------------------------------------------------------

  IF (my_cart_id == 0) THEN
    PRINT *, '           initialize background aerosol (aerdis)'
  ENDIF

  ! Allocate fields
  ist = 0
  ALLOCATE ( pvdaes(ke1), pvdael(ke1), pvdaeu(ke1), pvdaed(ke1), STAT=ist)
  IF (ist /= 0) THEN
    CALL model_abort(my_cart_id, 2777, 'Error setting up background aerosol', 'radiation_init')
  ENDIF

  zsign(1) = 0._wp
  DO k = 2, ke1
    zsign(k) = vcoord%sigm_coord(k)
  ENDDO

  CALL aerdis ( zsign, pvdaes, pvdael, pvdaeu, pvdaed, ke1,     &
                ptrbga, pvobga, pstbga, paeops, paeopl, paeopu, &
                paeopd, ptrpt , paeadk, paeadm)

!------------------------------------------------------------------------------
! End of the subroutine 
!------------------------------------------------------------------------------

END SUBROUTINE radiation_init 

!==============================================================================
!==============================================================================
!+ To initialize the organizational variables for the coarse radiation grid
!------------------------------------------------------------------------------

SUBROUTINE radiation_radcoarse_init

!------------------------------------------------------------------------------
!
! Description:
!  This subroutine collects grid points from the full COSMO grid into groups 
!  for the coarse grid. nradcoarse points in x- and in y-direction are taken 
!  for one group. 
!
!  First all groups are determined on the full grid without taking 
!  parallelization into account. 
!   - There are ietot_rc groups in x- and jetot_rc groups in y-direction. 
!   -  Necessary information for each group is stored in the variable 
!      gp_radcoarse_tot (of type gp_radcoarse)
!        = start- and end-indices of the group related to the total grid
!        = number of points in this group
!        = factor for averaging (1 / number of points)
! 
!  Then it is determined, which groups are computed by which subdomain.
!  A group is computed by subdomain A, if all grid points of the group
!  belong to subdomain A, including the boundary lines! Some groups are
!  computed in different subdomains then, but this avoides an additional
!  boundary exchange later on.
!  For every subdomain holds:
!   - There are ie_rc groups in x- and je_rc groups in y-direction
!   -  Necessary information for each group is stored in the variable 
!      gp_radcoarse_loc (of type gp_radcoarse)
!        = start- and end-indices of the group related to the subdomain grid
!        = number of points in this group
!        = factor for averaging (1 / number of points)
!   - The number of blocks (nproma_rad) is computed, taking the same value
!     for nproma as all other parameterizations do: nproma_rad = nproma.
!   - The number of groups in the last block (nlastproma_rad) is computed
!     (which differs from nlastproma for other parameterizations)
!
!  Last, the correspondance between the COSMO full grid and the blocked data
!  structure is computed.
!
!------------------------------------------------------------------------------

! Subroutine arguments: None
! --------------------

! Local arrays and scalars:
! -------------------------
  INTEGER                   ::                                           &
    i, j, ig, jg, ic, jc, iploc, ibloc, i_low, j_low, i_hig, j_hig,      &
    nij, iplocend, izstat, isub, nrb, n, iiii, ib, ip, jp, ijp, jploc

!- End of header
!------------------------------------------------------------------------------

  izstat = 0

!------------------------------------------------------------------------------
! Section 1: Determine number of coarse radiation grid groups 
!------------------------------------------------------------------------------

  ! Number of coarse radiation grid points in west-east direction: ietot_rc
  IF (MOD(ie_tot, nradcoarse) == 0) THEN
    ietot_rc = ie_tot / nradcoarse
  ELSE
    ietot_rc = ie_tot / nradcoarse + 1
  ENDIF

  ! Number of coarse radiation grid points in south-north direction: jetot_rc
  IF (MOD(je_tot, nradcoarse) == 0) THEN
    jetot_rc = je_tot / nradcoarse
  ELSE
    jetot_rc = je_tot / nradcoarse + 1
  ENDIF

  ! There are ietot_rc*jetot_rc groups of coarse radiation grid points:
  ! allocate gp_radcoarse_tot accordingly:

  ALLOCATE (gp_radcoarse_tot (ietot_rc,jetot_rc),           STAT=izstat)
  IF (izstat /= 0) THEN
    CALL model_abort(my_cart_id, 2810, 'Error allocating coarse grid org variables',  &
                                       'radiation_coarse_init')
  ENDIF


  ! Determine start- and end-indices for every group in grid point numbering of the
  ! total domain
  ie_rc = 0
  je_rc = 0

  DO jc = 1, jetot_rc
    DO ic = 1, ietot_rc
      gp_radcoarse_tot(ic,jc)%ilow  = (ic-1)*nradcoarse + 1
      gp_radcoarse_tot(ic,jc)%ihigh = MIN (ic*nradcoarse, ie_tot)
      gp_radcoarse_tot(ic,jc)%jlow  = (jc-1)*nradcoarse + 1
      gp_radcoarse_tot(ic,jc)%jhigh = MIN (jc*nradcoarse, je_tot)

      gp_radcoarse_tot(ic,jc)%ntgp  =                                         &
         (gp_radcoarse_tot(ic,jc)%ihigh - gp_radcoarse_tot(ic,jc)%ilow + 1) * &
         (gp_radcoarse_tot(ic,jc)%jhigh - gp_radcoarse_tot(ic,jc)%jlow + 1)

      gp_radcoarse_tot(ic,jc)%rfacgp = 1.0_wp / REAL (gp_radcoarse_tot(ic,jc)%ntgp, wp)

      ! At the same time determine the number of groups belonging to this 
      ! subdomain my_cart_id
      !   ie_rc: number of groups in i-direction
      !   je_rc: number of groups in j-direction
      IF (jc == 1) THEN
        ! number of groups in i-direction
        IF ( (isubpos(my_cart_id,1)-nboundlines <= gp_radcoarse_tot(ic,jc)%ilow) .AND.  &
             (gp_radcoarse_tot(ic,jc)%ihigh <= isubpos(my_cart_id,3)+nboundlines) ) THEN
          ie_rc = ie_rc + 1
        ENDIF
      ENDIF
        
    ENDDO

    ! number of groups in j-direction
    IF ( (isubpos(my_cart_id,2)-nboundlines <= gp_radcoarse_tot(1 ,jc)%jlow) .AND.      &
         (gp_radcoarse_tot(1 ,jc)%jhigh <= isubpos(my_cart_id,4)+nboundlines) ) THEN
      je_rc = je_rc + 1
    ENDIF
  ENDDO

!------------------------------------------------------------------------------
! Section 2: Determine number of blocks: nblock_rad and nlastproma
!            (different from nblock, nlastproma)
!            and local information for groups computed in this subdomain
!------------------------------------------------------------------------------

  ! Determine number of blocks necessary for the coarse radiation grid: nblock_rad
  ! The same nproma is assumed as in the other full grid parameterizations
  nproma_rad = nproma
  nij        = ie_rc * je_rc
  IF (MOD(nij, nproma_rad) == 0) THEN
    nblock_rad      = INT(nij / nproma_rad)
    nlastproma_rad  = nproma_rad
  ELSE
    nblock_rad      = INT(nij / nproma_rad) + 1
    nlastproma_rad  = MOD(nij,  nproma_rad)
  ENDIF

  ! The structure for storing the local informations about the grid point groups
  ! for the coarse radiation grid is already stored in block structure 
  ! (nproma_rad, nblock_rad)
  ALLOCATE (gp_radcoarse_loc (nproma_rad,nblock_rad),    STAT=izstat)
  IF (izstat /= 0) THEN
    CALL model_abort(my_cart_id, 2812, 'Error allocating coarse grid org variables',    &
                                       'radiation_coarse_init')
  ENDIF

!if (my_cart_id == 0) then
!  print *, ' Information about coarse radiation grid: '
!  print *, '   Number of groups: ietot_rc, jetot_rc = ', ietot_rc, jetot_rc
!  print *, '   List of groups: '
!  nrb = 0
!  do jc = 1, jetot_rc
!    do ic = 1, ietot_rc
!      nrb = nrb + 1
!      WRITE (*,'(A,3I5,A,4I5,I7,F10.2)') '     ', nrb, ic, jc,  '          ',           &
!                               gp_radcoarse_tot(ic,jc)%ilow,                            &
!                               gp_radcoarse_tot(ic,jc)%ihigh,                           &
!                               gp_radcoarse_tot(ic,jc)%jlow,                            &
!                               gp_radcoarse_tot(ic,jc)%jhigh,                           &
!                               gp_radcoarse_tot(ic,jc)%ntgp,                            &
!                               gp_radcoarse_tot(ic,jc)%rfacgp
!    enddo
!  enddo
!endif

  ! Now browse again through the full domain and determine all groups of 
  ! coarse radiation points belonging to this subdomain

  ibloc = 1
  iploc = 0
  DO jc = 1, jetot_rc
    IF ( (isubpos(my_cart_id,2)-nboundlines <= gp_radcoarse_tot(1 ,jc)%jlow) .AND.      &
         (gp_radcoarse_tot(1 ,jc)%jhigh <= isubpos(my_cart_id,4)+nboundlines) ) THEN

      DO ic = 1, ietot_rc
        IF ( (isubpos(my_cart_id,1)-nboundlines <= gp_radcoarse_tot(ic,jc)%ilow) .AND.  &
             (gp_radcoarse_tot(ic,jc)%ihigh <= isubpos(my_cart_id,3)+nboundlines) ) THEN

          iploc = iploc + 1
          IF (iploc > nproma_rad) THEN
            iploc = 1
            ibloc = ibloc + 1
          ENDIF

          ! determine local low- and high-indices for i-direction
          DO i = 1, ie
            ig = i_global(i)
            IF (ig == gp_radcoarse_tot(ic,jc)%ilow ) gp_radcoarse_loc(iploc,ibloc)%ilow  = i
            IF (ig == gp_radcoarse_tot(ic,jc)%ihigh) gp_radcoarse_loc(iploc,ibloc)%ihigh = i
          ENDDO

          ! determine local low- and high-indices for j-direction
          DO j = 1, je
            jg = j_global(j)
            IF (jg == gp_radcoarse_tot(ic,jc)%jlow ) gp_radcoarse_loc(iploc,ibloc)%jlow  = j
            IF (jg == gp_radcoarse_tot(ic,jc)%jhigh) gp_radcoarse_loc(iploc,ibloc)%jhigh = j
          ENDDO

          gp_radcoarse_loc(iploc,ibloc)%ntgp   = gp_radcoarse_tot(ic,jc)%ntgp
          gp_radcoarse_loc(iploc,ibloc)%rfacgp = gp_radcoarse_tot(ic,jc)%rfacgp

        ENDIF
      ENDDO   ! ic-loop
    ENDIF
  ENDDO       ! jc-loop


!do n = 0, num_compute-1
!  if (my_cart_id == n) then
!    print *, ' Information about groups in subdomain:  ', n, ':  ', ie_rc, je_rc
!    nrb = 0
!    do ibloc = 1, nblock_rad
!      if (ibloc == nblock_rad) then
!        iplocend = nlastproma_rad
!      else
!        iplocend = nproma_rad
!      endif
!      do iploc = 1, iplocend
!        nrb = nrb + 1
!        WRITE (*,'(A,3I5,A,4I5,I7,F10.2)') '     ', nrb, iploc, ibloc,':   ', &
!                                 gp_radcoarse_loc(iploc,ibloc)%ilow,          &
!                                 gp_radcoarse_loc(iploc,ibloc)%ihigh,         &
!                                 gp_radcoarse_loc(iploc,ibloc)%jlow,          &
!                                 gp_radcoarse_loc(iploc,ibloc)%jhigh,         &
!                                 gp_radcoarse_loc(iploc,ibloc)%ntgp,          &
!                                 gp_radcoarse_loc(iploc,ibloc)%rfacgp
!      enddo
!    enddo
!  endif
!  iiii  = dosleep(5)
!enddo

!------------------------------------------------------------------------------
! Section 3: Determine correspondance between full COSMO grid and coarse
!            radiation grid
!------------------------------------------------------------------------------

  ! Allocate and initialize mind_ilon_rad, mind_jlat_rad: the correspondance 
  ! between the blocked (coarse) radiation grid and the COSMO grid
  ALLOCATE (mind_ilon_rad(nproma_rad*nradcoarse, nradcoarse, nblock_rad),      &
            mind_jlat_rad(nproma_rad*nradcoarse, nradcoarse, nblock_rad),      &
                                                         STAT=izstat)
  IF (izstat /= 0) THEN
    CALL model_abort(my_cart_id, 2815, 'Error allocating mind_ilon_rad/jlat_rad_rad',  &
                                       'radiation_coarse_init')
  ENDIF


  ! To determine the correspondence we scan through all groups of coarse radiation
  ! grid points and then go over all grid points, starting from ilow and jlow, resp.
  DO ib = 1, nblock_rad
    DO jp = 1, nradcoarse
      IF (ib == nblock_rad) THEN
        iplocend = nlastproma_rad
      ELSE
        iplocend = nproma_rad
      ENDIF
      DO ip = 1, iplocend
        i = gp_radcoarse_loc(ip,ib)%ilow - 1

        j = gp_radcoarse_loc(ip,ib)%jlow - 1 + jp
        IF (j > je) THEN
          ! due to the construction, this can only happen at the northern
          ! border of the total domain
          j = je
        ENDIF

        DO ijp = 1, nradcoarse
          i = i+1
          IF (i > ie) THEN
            ! due to the construction, this can only happen at the eastern
            ! border of the total domain
            i = ie
          ENDIF

          mind_ilon_rad ((ip-1)*nradcoarse+ijp, jp, ib) = i
          mind_jlat_rad ((ip-1)*nradcoarse+ijp, jp, ib) = j
        ENDDO
      ENDDO
    ENDDO
  ENDDO

!do n = 0, num_compute-1
!  if (my_cart_id == n) then
!    print *, ' Information about mind-vars in subdomain:  ', n, ':  '
!    nrb = 0
!
!    do ibloc = 1, nblock_rad
!      if (ibloc == nblock_rad) then
!        iplocend = nlastproma_rad
!      else
!        iplocend = nproma_rad
!      endif
!      do jploc = 1, nradcoarse
!        do iploc = 1, iplocend*nradcoarse
!          WRITE (*,'(A,3I5,A,2I7)') '    Grid Correspondance:  ', iploc, jploc, ibloc, ':   ', &
!                         mind_ilon_rad(iploc,jploc,ibloc), mind_jlat_rad(iploc,jploc,ibloc)
!        enddo
!      enddo
!    enddo
!  endif
!  iiii = dosleep(5)
!enddo

!------------------------------------------------------------------------------
! End of the subroutine 
!------------------------------------------------------------------------------

END SUBROUTINE radiation_radcoarse_init 

!==============================================================================
!==============================================================================
!+ To clean up at the end of the program
!------------------------------------------------------------------------------    

SUBROUTINE radiation_finalize

!------------------------------------------------------------------------------
!
! Description:
!
! Method:
!
!------------------------------------------------------------------------------

  INTEGER :: izerror
  CHARACTER(LEN=80) :: yerrmsg

!----------- End of header ----------------------------------------------------

  izerror = 0

  ! Deallocate mind_ilon/jlat_rad
  DEALLOCATE (mind_ilon_rad, mind_jlat_rad)

  ! Deallocate the working arrays

#ifdef _OPENACC
  CALL radiation_in_wkarr_dealloc(.TRUE., izerror)
  IF (izerror /= 0) THEN
      yerrmsg = 'Error in radiation_in_wkarr_dealloc'
      CALL model_abort(my_cart_id, 666, yerrmsg, 'radiation_init')
  ENDIF

  CALL radiation_rg_wkarr_dealloc(izerror)
  IF (izerror /= 0) THEN
      yerrmsg = 'Error in radiation_rg_wkarr_dealloc'
      CALL model_abort(my_cart_id, 666, yerrmsg, 'radiation_init')
  ENDIF
#endif

!------------------------------------------------------------------------------
! End of the subroutine 
!------------------------------------------------------------------------------

END SUBROUTINE radiation_finalize

!==============================================================================
!==============================================================================
!+ Module procedure in "Radiation" 
!------------------------------------------------------------------------------

SUBROUTINE radiation_organize (ydate_ini, ib, ipend, lzradstep, ntl)

!------------------------------------------------------------------------------
!
! Description:
!
!  The module procedure radiation_organize forms the interface between
!  the model and the radiation schemes (up to now only RG scheme)
!
! Method:
!
!  All variables that are required for the radiation code (i.e. input arrays
!  and scalars) are provided or calculated from the model variables.
!  The results are stored as solar and thermal heating rates on the 
!  corresponding global arrays (e.g. sohr and thhr).
!
!  Thermal variables are only computed in a full radiation step, while
!  solar variables are updated in every step, depending on the zenith angle
!  of the sun.
!
!  For the RG scheme, the following kind of variables are computed:
!  (subroutines contained in parallel_utilities.f90). Similar methods could most
!  probably also be used for other radiation schemes
!
!  For a full radiation step:
!   - Sunshine conditions (zenith angle zsmu0_fg for the time interval that is 
!     considered in the full radiation step)
!   - Surface albedo (solar and thermal)
!   - Cloudiness and humidity conditions: cloud cover is stored on global arrays
!   - Pressure and temperature on input for RG scheme
!   - Choose a CO2-scenario and compute amount of absorbers
!   - Correction factors for radiation in complex topography
!   - Call the radiation scheme
!   - Compute heating rates and radiation budgets (set COSMO variables)
!   - Compute the actual zenith angle for this time step
!   - Calculate solar fluxes and heating rates based on this actual zenith angle
!
!  For the update of solar variables, only the last two steps above are computed:
!   - Compute the actual zenith angle for this time step
!   - Calculate solar fluxes and heating rates based on this actual zenith angle
!
!------------------------------------------------------------------------------
!
! Method for computing the radiation on a coarser grid:
!
! With nradcoarse > 1 it is possible to define a coarser grid for computing 
! the radiation fluxes. Computing for less grid points saves a considerable
! amount of computer time.
!
! Several grid points from the full COSMO grid (nradcoarse*nradcoarse) are grouped
! together and all input fields for the radiation schemes are averaged from all
! grid points of a group. The averaged values are put to a blocked data structure
! to call the routines of the radiation schemes. The results of one group are 
! stored to all COSMO grid points of this group.
!
! The surface outgoing thermal radiation (which is computed using the averaged
! coarse grid ground temperature, stored in ztg_rc) and the reflected solar radiation
! (calculated using the averaged coarse grid albedo, stored in alb_rad_rc) are
! substracted from thbs and sobs, resp., and are added again, but calculated with
! ground temperature and albedo from the full COSMO grid variables.
!
! In between, the radiative increments are filtered with a discrete filter 
! (analogous to lconf_avg), if lradf_avg == .TRUE.
! These computations are done in another subroutine radiation_average.
!
!------------------------------------------------------------------------------
!==============================================================================

! Parameterlist
! -------------

CHARACTER (LEN=14), INTENT(IN)     ::   &
     ydate_ini   ! start of the forecast yyyymmddhh (year, month, day, hour, min, sec)

LOGICAL, INTENT(IN) :: lzradstep    ! true if full radiation timestep, otherwise
                                    ! only recalculation of sohr, sobs, and sobt,
                                    ! based on actual zenith angle
INTEGER, INTENT(IN) :: ib, ipend, & ! number and length of actual block
                       ntl          ! time level to be considered


! Local parameters: 
! ----------------

  LOGICAL, PARAMETER ::          &
    lcrf             = .FALSE. , &
    ldebug           = .FALSE.

! Local scalars:
! -------------

! Input for the radiation routine fesft
! -------------------------------------
  INTEGER                   ::  &
     ki1sd,       & ! start index for first  array dimension
     ki1ed,       & ! end   index for first  array dimension

   ! and the same for the computations
     ki1sc,       & ! start index for first  array computation
     ki1ec          ! end   index for first  array computation

  LOGICAL                  ::  &
     lrady,       & ! for radiation on coarse-grid !T.R.
     lsolar         ! control switch for solar calculations

  REAL    (KIND=wp)        ::  &
     zsct           ! solar constant (at time of year)

  INTEGER                  ::  &
    kzdims(24),            &
    j_rn,nradcoarse_y,     & ! for radiation on coarse grid! T.R.
    ii, ista, iend,        & ! for radiation on coarse grid! T.R.
    i, j, js, k, i_ld,     & ! loop indices over spatial dimensions
        ip, jp,            & ! for blocked data structure
              i_std, i_etd,& ! time level of prognostic variables
    nzrad,                 & !
    jj   , itaja,          & ! output from routine get_utc_dat
    ist  ,                 & ! loop index for soil type  
    ipdim,                 & ! ipend * nradcoarse: i-dimension
    izstata,               & ! error status at allocation !T.R.
    izstatd,               & ! error status at deallocation !T.R.
    izerror, izdebug         ! for error status
 
  REAL    (KIND=wp)        ::  &
    zstunde,                   & ! output from routine get_utc_dat
    zpnf   ,zphf   , zphfo  ,  & !
    zfac

  CHARACTER (LEN=14)   yrad1         ! output from routine get_utc_dat
  CHARACTER (LEN=28)   yrad2         ! output from routine get_utc_dat
  CHARACTER (LEN=255)  yzerrmsg      ! for error message
  CHARACTER (LEN=25)   yzroutine
 
!- End of header
!==============================================================================
 
!------------------------------------------------------------------------------
! Start GPU data region
!------------------------------------------------------------------------------

  !$acc data                                                            &
  !$acc present ( mind_ilon_rad, mind_jlat_rad                        ) &
  !$acc present ( rlat, rlon, sun_azi, sun_el, p0, dp0, p0hl, t, t_g  ) &
  !$acc present ( alb_rad, swdir_cor                                  ) &
  !$acc present ( zalso, zalth, zsmu0_fg, zsmu0, zskyview, zti, zdpr  ) &
  !$acc present ( ztg_rc, zapre, zsw, zwv, zclwc, zciwc, zclc         ) &
  !$acc present ( zduco2f, zduo3f, zaeq1, zaeq2, zaeq3, zaeq4, zaeq5  ) &
  !$acc present ( zwork3d1, zwork3d2, zwork3d3, zwork3d               ) &
  !$acc present ( zwork2d1, zwork2d2, zwork2d3, zwork2d4, zwork2d5    ) &
  !$acc present ( zwork2d6, zwork2d7, zwork2d8, nzwork                ) &
  !$acc present ( zflt,  zfls,  zflt_s, zfls_s, zflsdir               ) &
  !$acc present ( zfltd, zfltu, zflsd,  zflsu,  zflsp                 ) &
  !$acc present ( zflpar_s, zflpar_s_difu, zflpar_s_difd, zflpar_s_dir) &
  !$acc present ( sohr, sotr, sotr_par, thhr, sobs, thbs, pabs, sobt  ) &
  !$acc present ( thbt, alb_rad_rc, sodwddm, swdir_s, swdifd_s        ) &
  !$acc present ( swdifu_s, swtrdir_s, swtrdifd_s, swtrdifu_s         ) &
  !$acc present ( lwd_s, lwu_s, sod_t, asod_t, skyview                )

!------------------------------------------------------------------------------
! Section 0: Basic Initializations
!------------------------------------------------------------------------------

  yzroutine = 'radiation_organize'
  yzerrmsg  = ' '
  izerror   = 0
  izstata   = 0
  ipdim     = ipend * nradcoarse

  IF (ldebug_rad) THEN
    IF (lprintdeb_all) THEN
      izdebug = idbg_level
    ELSE
      IF (my_cart_id == 0) THEN
        izdebug = idbg_level
      ELSE
        izdebug = 0
      ENDIF
    ENDIF
  ELSE
    izdebug = 0
  ENDIF

  IF (lzradstep) THEN

    IF (izdebug > 10) THEN
      PRINT *, '        START radiation_organize for a full radiation step for block: ', ib
    ENDIF

    !----------------------------------------------------------------------------
    ! Section 1: Determine sunshine-condition
    !----------------------------------------------------------------------------

    IF (izdebug > 10) THEN
      PRINT *, '           determine sunshine conditions'
    ENDIF

    ! Compute the sunshine conditions for the computation of the radiative
    ! balance in fesft. The mean conditions over all timesteps until the next 
    ! full radiation timestep is computed.

    CALL compute_sunshine_conditions(                       &
        ydate_ini, izdebug,                                 &
        1        , ipdim     , 1        , nradcoarse, ib,   &
        ntstep, nextrad-1, rlat, rlon, itaja_zsct_previous, &
        zsct_save, zdekcos_save, zdeksin_save, zdtzgl_save, &
        sun_azi, sun_el, zsct, zeit0, zsmu0_fg,             &
        nzwork, zwork2d5, zwork2d6, zwork2d7, zwork2d8, lsolar)

    !------------------------------------------------------------------------------
    ! Section 2:  Calculation of surface albedo taking soil type,              
    !             vegetation and snow/ice conditions into account
    !------------------------------------------------------------------------------

    IF (izdebug > 10) THEN
      PRINT *, '           compute surface albedo'
    ENDIF

    CALL surface_albedo (1, ipdim, 1, nradcoarse, ib, ntl, zalso, zalth)

    ! Save the surface albedo for output: in alb_rad
    DO jp = 1, nradcoarse
      !$acc parallel
      !$acc loop gang vector
      DO ip = 1, ipdim
        ! get i/j indices for 2D COSMO data structure
        i = mind_ilon_rad(ip,jp,ib)
        j = mind_jlat_rad(ip,jp,ib)

        alb_rad (i,j) = zalso (ip,jp) !T.R. fuer Albedokorrektur
      ENDDO
      !$acc end parallel
    ENDDO

    !------------------------------------------------------------------------------
    ! Section 3:  Set cloudiness and humidity on input for fesft; 
    !             Store cloud cover on corresponding global arrays
    !------------------------------------------------------------------------------

    IF (izdebug > 10) THEN
      PRINT *, '           cloudiness and humidity'
    ENDIF

    CALL cloudiness_and_humidity (qv_locrad, qc_locrad, qi_locrad, 1, ipdim,  &
               1  , nradcoarse, ib,   ke, ntl, lprog_qi,                      &
               zsw, zwv, zclwc, zciwc, zclc,                                  &
               zwork3d1, zwork3d2, zwork3d3, zwork2d1)

    !------------------------------------------------------------------------------
    ! Section 4:  Set pressure and temperature on input for fesft; 
    !------------------------------------------------------------------------------

    IF (izdebug > 10) THEN
      PRINT *, '           pressure and temperature'
    ENDIF

    ! Surface pressure, half level pressure and pressure thickness. 
    ! At present, pressure is replaced by model reference pressure 
    ! for radiation calculations

    ! Temperatures at layer boundaries
    !$acc parallel
    !$acc loop gang vector collapse(3)
    DO  k = 2, ke
      DO jp = 1, nradcoarse
        DO ip = 1, ipdim
          ! get i/j indices for 2D COSMO data structure
          i = mind_ilon_rad(ip,jp,ib)
          j = mind_jlat_rad(ip,jp,ib)

          zpnf     = p0hl(i,j,k  )
          zphf     = p0  (i,j,k  )
          zphfo    = p0  (i,j,k-1)
          zti(ip,k,jp) = (t(i,j,k-1,ntl)*zphfo*(zphf - zpnf )      &
            + t(i,j,k  ,ntl)*zphf *(zpnf - zphfo) )    &
            * (1.0_wp/(zpnf *(zphf - zphfo)))

          ! pressure thickness
          zdpr(ip,k,jp) = dp0 (i,j,k)
        ENDDO
      ENDDO
    ENDDO
    !$acc end parallel

    !$acc parallel
    !$acc loop gang vector collapse(2)
    DO jp = 1, nradcoarse
      DO ip = 1, ipdim
        ! get i/j indices for 2D COSMO data structure
        i = mind_ilon_rad(ip,jp,ib)
        j = mind_jlat_rad(ip,jp,ib)

        zpnf    = p0hl  (i,j,2)
        zphf    = p0    (i,j,1)
        zti(ip,ke1,jp) = t_g(i,j,ntl)
        zti(ip,  1,jp) = t  (i,j,1,ntl) - zphf*(t(i,j,1,ntl)-zti(ip,2,jp))/(zphf - zpnf)

        ! pressure at lowest boundary
        zapre(ip,jp) = p0hl(i,j,ke+1)

        ! pressure thickness
        zdpr(ip,1,jp) = dp0 (i,j,1)
      ENDDO
    ENDDO
    !$acc end parallel

    !----------------------------------------------------------------------------
    ! Section 5: Choose CO2 scenario and compute amount of absorbers
    !----------------------------------------------------------------------------

    IF (izdebug > 10) THEN
      PRINT *, '           CO2 scenario and absorbers'
    ENDIF

    IF ((ntstep >= 1) .OR. (nincrad == 1)) THEN
      nzrad   = ntstep + nincrad/2
    ELSEIF (ntstep==0) THEN
      nzrad = 0
    ENDIF

    CALL get_utc_date ( nzrad, ydate_ini, dt, itype_calendar, yrad1, yrad2,  &
                        itaja, zstunde )
    READ (yrad1(1:4),'(I4)') jj

    CALL co2_scenarios_and_absorbers (itaja, jj, 1, ipdim, 1, nradcoarse,    &
              ib , ke, izdebug, zti,                                         &
              zduco2f, zduo3f, zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,            &
              zwork2d1, zwork2d2, zwork2d3, zwork2d4, zwork2d5, zwork2d6,    &
              zwork2d7, zwork2d8 )

    !------------------------------------------------------------------------------
    ! Section 6:  Correction factors for radiation in complex topography
    !------------------------------------------------------------------------------

    IF (lradtopo) THEN
      IF (izdebug > 10) THEN
        PRINT *,'        radiation_organize with lradtopo = ', lradtopo
      ENDIF

      !US:  Careful: this does NOT work with nradcoarse > 1!!!
      CALL calc_rad_corrections (slo_ang, slo_asp, horizon, zsmu0_fg,   &
        rlat, rlon, zdeksin_save, zdekcos_save, zeit0, swdir_cor,       &
        zwork2d1, zwork2d2, zwork2d3,                                   &
        ie, je, nhori, ib, 1    , ipdim  , 1        , nradcoarse,       &
        izdebug, izerror, yzerrmsg)
      IF (izerror /= 0) THEN
        CALL model_abort(my_cart_id, izerror, TRIM(yzerrmsg), 'radiation_interface')
      ENDIF

      ! and set zskyview (1-dimensional slice)
      !$acc parallel
      !$acc loop gang vector collapse(2)
      DO jp = 1, nradcoarse
        DO  ip = 1, ipdim
          ! get i/j indices for 2D COSMO data structure
          i = mind_ilon_rad(ip,jp,ib)
          j = mind_jlat_rad(ip,jp,ib)

          zskyview(ip,jp) = skyview(i,j)
        ENDDO
      ENDDO
      !$acc end parallel
    ELSE
      ! Set default value for skyview
      !$acc parallel
      !$acc loop gang vector collapse(2)
      DO jp = 1, nradcoarse
        DO  ip = 1, ipdim
          zskyview(ip,jp) = 1.0_wp
        ENDDO
      ENDDO
      !$acc end parallel
    ENDIF

    !------------------------------------------------------------------------------
    ! Section 7:  Set start- and end-indices; in case of coarse grid average values
    !             and call RG scheme
    !------------------------------------------------------------------------------

    IF (izdebug > 10) THEN
      PRINT *, '           memory allocation for fesft-interface fields'
    ENDIF

    IF (nradcoarse > 1) THEN

      IF (izdebug > 10) THEN
        PRINT *, '           calculations for radiation averaging  ', nradcoarse
      ENDIF

      ! Setting of array boundaries for routine fesft
      ki1sd=1
      ki1ed=nproma_rad * nradcoarse
      ki1sc=1
      ki1ec=ipend

      CALL average_to_coarse_grid

    ELSE ! no coarse grid

      IF (izdebug > 10) THEN
        PRINT *, '           settings for no radiation averaging'
      ENDIF

      ! Setting of array boundaries for routine fesft
      ki1sd = 1
      ki1ed = nproma
      ki1sc = 1
      ki1ec = ipend

      ! set zsmu0: this is the only field we need twice, because we have to
      ! retain the value on the full grid
      !$acc parallel
      !$acc loop gang vector
      DO ip = 1, ipdim
        zsmu0(ip) = zsmu0_fg(ip,1)
      ENDDO
      !$acc end parallel

#if defined TWOMOM_SB && defined COSMO_ART
      ! Calculate cloud optical properties based on cloud effective radii
      ! ------------------------------------------------------------------
      IF(iradpar_cloud > 1) THEN
        CALL calc_cloud_opt(js,iradpar_cloud)
      ENDIF
#endif
    ENDIF

    CALL fesft                                                                        &
      (zti  (:,:, 1), zdpr (:,:, 1), zclc   (:,:, 1), zwv   (:,:, 1), zsw(:,:, 1),    &
       zclwc(:,:, 1), zciwc(:,:, 1), zduco2f(:,:, 1), zduo3f(:,:, 1),                 &
       zaeq1(:,:, 1), zaeq2(:,:, 1), zaeq3  (:,:, 1), zaeq4 (:,:, 1), zaeq5(:,:, 1),  &
       zapre(:, 1)  , zsmu0(:)     , zalso  (:, 1)  , zalth (:, 1)  , zskyview(:, 1), &
       swdir_cor(:, 1),   sigma,       zsct,                                          &
       ki1sd,     ki1ed,                               1    ,    ke   ,               &
       ki1sc,     ki1ec,                               1    ,    ke   ,               &
       lsolar,    lcrf,             lradtopo,         izdebug,         1,             &
       zflt,      zfls,             zflt_s,           zfls_s,         zflsdir,        &
       zfltd,     zfltu,            zflsd,            zflsu,          zflsp,          &
       zflpar_s,  zflpar_s_difu,    zflpar_s_difd,    zflpar_s_dir,                   &
       izerror, yzerrmsg)

    IF (izerror /= 0) THEN
      CALL model_abort(my_cart_id, izerror, TRIM(yzerrmsg), 'radiation_interface')
    ENDIF

    !------------------------------------------------------------------------------
    ! Section 8:  Heating rates and radiation budget at surface
    !------------------------------------------------------------------------------

    !$acc parallel
    !$acc loop gang vector collapse(2)
    DO  jp = 1, nradcoarse
      DO  ip = 1, ipend
        ! Note: the result fields from fesft are on the coarse radiation grid (if used)
        !       and the ip-loop only is over the coarse blocks (ipend and not ipend*nradcoarse)
        !       To get all COSMO grid points, another loop is necessary now
        !       (when not using a coarse grid, this loop has only one cycle

        IF (nradcoarse > 1) THEN
          ista = (ip-1) * nradcoarse + 1
          iend =  ip    * nradcoarse
        ELSE
          ista = ip
          iend = ip
        ENDIF

        DO ii = ista, iend

          ! get i/j indices for 2D COSMO data structure
          i = mind_ilon_rad(ii,jp,ib)
          j = mind_jlat_rad(ii,jp,ib)

          swtrdir_s (i,j) = 0.0_wp
          swtrdifd_s(i,j) = 0.0_wp
          swtrdifu_s(i,j) = 0.0_wp
          sotr(i,j,ke1)   = 0.0_wp
          sotr_par(i,j)   = 0.0_wp
          sod_t     (i,j) = 0.0_wp
          sodwddm(i,j)    = 0.0_wp
          IF (zsmu0_fg(ii,jp) > zepemu) THEN
            sod_t     (i,j)     = zsmu0_fg (ii,jp)*zsct
            sotr_par  (i,j)     = zflpar_s(ip) / sod_t(i,j) !Calculate transmiss. from flux
            sotr      (i,j,ke1) = zfls_s  (ip) / sod_t(i,j) !Calculate transmiss. from flux
            swtrdir_s (i,j)     = zflsp   (ip) / sod_t(i,j) !Calculate transmiss. from flux
            swtrdifd_s(i,j)     = zflsd   (ip) / sod_t(i,j) !Calculate transmiss. from flux
            swtrdifu_s(i,j)     = zflsu   (ip) / sod_t(i,j) !Calculate transmiss. from flux
            ! for the Climate-LM Version: solar direct parallel radiation at the surface
            sodwddm   (i,j)     = zflsdir (ip,ke1) / zsmu0_fg(ii,jp)
          END IF
          thbs    (i,j) = zflt_s(ip)
          thbt    (i,j) = zflt  (ip,1)
!US       tg_ra   (i,jz1) = tg_rn (i)
          lwd_s   (i,j) = zfltd (ip)
          lwu_s   (i,j) = zfltu (ip)

!US these fields belong to a development which did never make it into the official
!   version. But if this development will be activated again (see also subroutine
!   radiation_average in radiation_interface) they are needed as global fields
!   (also activate in data_fields and src_allocation)
!         flpar_s_dir (i,j) = zflpar_s_dir (i)
!         flpar_s_difd(i,j) = zflpar_s_difd(i)
!         flpar_s_difu(i,j) = zflpar_s_difu(i)

        ENDDO
      ENDDO
    ENDDO
    !$acc end parallel

    IF (nradcoarse > 1) THEN
      ! store back the averaged ground temperature (to compute the surface 
      ! outgoing thermal radiation) and the averaged solar albedo (to compute 
      ! reflected solar radiation) with coarse grid values
      !$acc parallel
      !$acc loop gang vector collapse(2)
      DO  jp = 1, nradcoarse
        DO  ip = 1, ipend
          IF (nradcoarse > 1) THEN
            ista = (ip-1) * nradcoarse + 1
            iend =  ip    * nradcoarse
          ELSE
            ista = ip
            iend = ip
          ENDIF

          DO ii = ista, iend
            ! get i/j indices for 2D COSMO data structure
            i = mind_ilon_rad(ii,jp,ib)
            j = mind_jlat_rad(ii,jp,ib)
            alb_rad_rc(i,j) = zalso(ip,1)
            ztg_rc    (i,j) = zti  (ip,ke1,1)
          ENDDO
        ENDDO
      ENDDO
      !$acc end parallel
    ENDIF

    !$acc parallel
    !$acc loop vector collapse(3)
    DO  k = 1, ke 
      DO  jp = 1, nradcoarse
        DO  ip = 1, ipend
          IF (nradcoarse > 1) THEN
            ista = (ip-1) * nradcoarse + 1
            iend =  ip    * nradcoarse
          ELSE
            ista = ip
            iend = ip
          ENDIF

          DO ii = ista, iend

            ! get i/j indices for 2D COSMO data structure
            i = mind_ilon_rad(ii,jp,ib)
            j = mind_jlat_rad(ii,jp,ib)

            zfac = g / ( cp_d * dp0(i,j,k) )
            sotr(i,j,k) = 0.0_wp      
            IF (zsmu0_fg(ii,jp) > zepemu) THEN
              sotr(i,j,k) = zfls(ip,k) / sod_t(i,j) !Calculate transmiss. from flux
            ENDIF
            thhr(i,j,k)   = zfac * (zflt(ip,k)-zflt(ip,k+1))

          ENDDO

        ENDDO
      ENDDO
    ENDDO
    !$acc end parallel

!US this will not work in the new data structure
!US #if defined TWOMOM_SB && defined COSMO_ART
!US     IF(iradpar_cloud > 1) THEN
!US !very old     CALL calc_optthick_cloud(js,zclc(:,js,:),zclwc(:,js,:),zciwc(:,js,:),lsolar)
!US !fromKIT      CALL calc_optthick_cloud(js,dp0(:,js,:),zclwc(:,js,:),zciwc(:,js,:),lsolar)
!US !possible     CALL calc_optthick_cloud(zdpr(:,:,1),zclwc(:,:,1),zciwc(:,:,1),lsolar)
!US     ENDIF
!US #endif

  ENDIF ! full radiation step:  lzradstep

  !------------------------------------------------------------------------------
  ! Section 10: Calculate actual zenith angle
  !------------------------------------------------------------------------------

  CALL compute_sunshine_conditions(                           &
      ydate_ini, izdebug,                                     &
      1,       ipdim,  1,  nradcoarse,  ib,                   &
      ntstep, ntstep, rlat, rlon, itaja_zsct_previous,        &
      zsct_save, zdekcos_save, zdeksin_save, zdtzgl_save,     &
      sun_azi, sun_el, zsct, zeit0, zsmu0_fg,                 &
      nzwork, zwork2d5, zwork2d6, zwork2d7, zwork2d8, lsolar)

  !------------------------------------------------------------------------------
  ! Section 11: Calculate solar fluxes and heating rates based on actual zenith angle
  !------------------------------------------------------------------------------
  
  !$acc parallel
  !$acc loop gang vector collapse(2)
  DO jp = 1, nradcoarse
    DO  ip = 1, ipend*nradcoarse
      ! Note: for zenith angle update the variable zsmu0_fg is not averaged to
      !       the coarse radiation grid, so we have to go through the full 
      !       field here
      ! get i/j indices for 2D COSMO data structure
      i = mind_ilon_rad(ip,jp,ib)
      j = mind_jlat_rad(ip,jp,ib)

      sod_t   (i,j) = 0.0_wp
      swdir_s (i,j) = 0.0_wp
      swdifd_s(i,j) = 0.0_wp
      swdifu_s(i,j) = 0.0_wp
      sobs    (i,j) = 0.0_wp
      sobt    (i,j) = 0.0_wp
      pabs    (i,j) = 0.0_wp
      IF (zsmu0_fg(ip,jp) > zepemu) THEN
        sod_t   (i,j) = zsmu0_fg  (ip,jp)    * zsct_save
        swdir_s (i,j) = swtrdir_s (i,j)     * sod_t(i,j)
        swdifd_s(i,j) = swtrdifd_s(i,j)     * sod_t(i,j)
        swdifu_s(i,j) = swtrdifu_s(i,j)     * sod_t(i,j)
        sobs    (i,j) = sotr      (i,j,ke1) * sod_t(i,j)
        sobt    (i,j) = sotr      (i,j,1)   * sod_t(i,j)
        pabs    (i,j) = sotr_par  (i,j)     * sod_t(i,j)
      ENDIF
    ENDDO
  ENDDO
  !$acc end parallel

  !$acc parallel
  !$acc loop gang vector collapse(3)
  DO k=1,ke
    DO jp = 1, nradcoarse
      DO  ip = 1, ipend*nradcoarse
        ! get i/j indices for 2D COSMO data structure
        i = mind_ilon_rad(ip,jp,ib)
        j = mind_jlat_rad(ip,jp,ib)

        sohr(i,j,k) = 0.0_wp
        IF (zsmu0_fg(ip,jp) > zepemu) THEN
          zfac        = g / ( cp_d * dp0(i,j,k) )        
          sohr(i,j,k) = (sotr(i,j,k)-sotr(i,j,k+1)) * zfac * sod_t(i,j)
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  !$acc end parallel
  
  !T.R.
  ! Next section (lradtopo.AND.lradave) commented out because
  ! this combination is currently considered not to work anyway.
!!$    IF (lradtopo.AND.lradave) THEN !Careful: this combination does not work!
!!$      ! Storage of individual flux components for thermal radiative surface flux
!!$
!!$      IF (.NOT. lemiss) THEN
!!$        zemissivity =         1.0_wp - ctalb ! surface emissivity
!!$        zemissfac   = sigma * (1.0_wp - ctalb)
!!$      ENDIF
!!$
!!$      DO js = jstartpar,jendpar
!!$        DO i = istartpar,iendpar
!!$
!!$          ! Recompute surface thermal flux components based on 
!!$          ! lower boundary condition
!!$          IF (lemiss) THEN
!!$            lwd_s(i,js) = (thbs(i,js) + sigma*emis_rad(i,js)*t_g(i,js,ntl)**4) / emis_rad(i,js)
!!$          ELSE
!!$            lwd_s(i,js) = (thbs(i,js) + zemissfac*t_g(i,js,ntl)**4) / zemissivity
!!$          ENDIF
!!$
!!$          ! correction for thermal fluxes
!!$          lwu_s(i,js) = lwd_s(i,js) - thbs(i,js)
!!$          lwd_s(i,js) = lwd_s(i,js) *               skyview(i,js) +        &
!!$            lwu_s(i,js) * (1.0_wp - skyview(i,js))
!!$
!!$          thbs(i,js)  = lwd_s(i,js) - lwu_s(i,js)
!!$
!!$          IF (smu0(i,js) > zepemu) THEN
!!$            ! correction for solar fluxes
!!$            ! direct down corrected
!!$            swdir_s(i,js) = swdir_cor(i,js) * swdir_s(i,js)
!!$
!!$            ! diffuse down corrected
!!$            swdifd_s(i,js) = swdifd_s(i,js) *             skyview(i,js)    &
!!$              + swdifu_s(i,js) * (1.0_wp-skyview(i,js))
!!$
!!$            ! diffuse up adapted to new other components
!!$            swdifu_s(i,js) = (swdir_s(i,js) + swdifd_s(i,js)) * alb_rad(i,js)
!!$            sobs    (i,js) =  swdir_s(i,js) + swdifd_s(i,js)  - swdifu_s(i,js)
!!$
!!$            ! correction for solar fluxes of photosynthetic active radiation
!!$            IF ( (zzflpar_s_dir(i,js) > 0.0_wp) .OR.                      &
!!$              (zzflpar_s_difd(i,js) > 0.0_wp)             ) THEN
!!$              zalbradtopo =  zzflpar_s_difu(i,js) /                            &
!!$                (zzflpar_s_dir(i,js)+zzflpar_s_difd(i,js))
!!$
!!$              ! direct down corrected
!!$              zzflpar_s_dir(i,js) = swdir_cor(i,js) * zzflpar_s_dir(i,js)
!!$
!!$              ! diffuse down corrected
!!$              zzflpar_s_difd(i,js) = zzflpar_s_difd(i,js) *             skyview(i,js) &
!!$                + zzflpar_s_difu(i,js) * (1.0_wp-skyview(i,js))
!!$
!!$              ! diffuse up adapted to new other components
!!$              zzflpar_s_difu(i,js) = (zzflpar_s_dir(i,js) + zzflpar_s_difd(i,js))      &
!!$                * zalbradtopo
!!$              pabs(i,js)       = zzflpar_s_dir(i,js) + zzflpar_s_difd(i,js)        &
!!$                - zzflpar_s_difu(i,js)
!!$            ELSE
!!$              pabs(i,js)       = 0.0_wp
!!$            ENDIF
!!$          ELSE
!!$            sobs(i,js) = 0.0_wp
!!$            pabs(i,js) = 0.0_wp
!!$          ENDIF
!!$        ENDDO
!!$      ENDDO
!!$
!!$    ENDIF !lradtopo.AND.lradave

! this must be put to block structure also????????
#ifdef COSMOART    
  IF (lzradstep) THEN
    ! T.R.
    ! It might be considered to call calcjval also every timestep.
    ! Then some adaptations would be necessary. 
    IF(l_cosmo_art .AND. lgas) THEN
     DO k = 1, ke1
         DO ip = 1, ipend
           ! get i/j indices for 2D COSMO data structure
           i = mind_ilon(ip,ib)
           j = mind_jlat(ip,ib)
           mmy(i,j)     = zsmu0_fg(ip,1)
           Edir(i,j,k)  = Edir_b(ip,k)
           Eup(i,j,k)   = Eup_b(ip,k)
           Edown(i,j,k) = Edown_b(ip,k)
          ENDDO
      ENDDO
      CALL calcjval
    ENDIF
  ENDIF !lzradstep
#endif

  !$acc end data

!==============================================================================

CONTAINS

!==============================================================================

SUBROUTINE average_to_coarse_grid

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine averages to a coarser grid:
!   The values of all grid points belonging to a group are summed up and 
!   at the end, the average is build. The results are written back to the
!   original variable zxx(n,k,1), where "n" denotes the first grid point of
!   that group.
!
!------------------------------------------------------------------------------

  INTEGER :: k, ip, izstart, izend, jzend, izsave, jz1, iz1

  REAL (KIND=wp) :: zti_loc     , zdpr_loc    , zclc_loc    , zwv_loc     ,  &
                    zsw_loc     , zclwc_loc   , zciwc_loc   , zduco2f_loc ,  &
                    zduo3f_loc  , zaeq1_loc   , zaeq2_loc   , zaeq3_loc   ,  &
                    zaeq4_loc   , zaeq5_loc   , ztike1_loc  , zalso_loc   ,  &
                    zalth_loc   , zapre_loc   , zsmu0_loc

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Start GPU data region
!------------------------------------------------------------------------------

  !$acc data                                                            &
  !$acc present ( gp_radcoarse_loc                                    ) &
  !$acc present ( zti, zdpr, zclc, zwv, zsw, zclwc, zciwc             ) &
  !$acc present ( zduco2f, zduo3f, zaeq1, zaeq2, zaeq3, zaeq4, zaeq5  ) &
  !$acc present ( zalso, zalth, zsmu0_fg, zsmu0, zapre                )

! Atmosphere                 !XL_COMMENT: this is probably not correct on the GPU because  
  !$acc parallel             !            the loop are not purely nested
                             !            It would require a special treatment for izsave,izstart
  !$acc loop gang
  DO k = 1, ke

    izsave = 0

    izstart = 0
    !acc loop vector
    DO ip = 1, ipend

      izsave   = izsave + 1

      izend = gp_radcoarse_loc(ip,ib)%ihigh - gp_radcoarse_loc(ip,ib)%ilow + 1
      jzend = gp_radcoarse_loc(ip,ib)%jhigh - gp_radcoarse_loc(ip,ib)%jlow + 1

      zti_loc     = 0.0_wp
      zdpr_loc    = 0.0_wp
      zclc_loc    = 0.0_wp
      zwv_loc     = 0.0_wp
      zsw_loc     = 0.0_wp
      zclwc_loc   = 0.0_wp
      zciwc_loc   = 0.0_wp
      zduco2f_loc = 0.0_wp
      zduo3f_loc  = 0.0_wp
      zaeq1_loc   = 0.0_wp
      zaeq2_loc   = 0.0_wp
      zaeq3_loc   = 0.0_wp
      zaeq4_loc   = 0.0_wp
      zaeq5_loc   = 0.0_wp

      DO jz1 = 1, jzend
        DO iz1 = izstart+1, izstart+izend
          zti_loc     = zti_loc     + zti    (iz1,k,jz1)
          zdpr_loc    = zdpr_loc    + zdpr   (iz1,k,jz1)
          zclc_loc    = zclc_loc    + zclc   (iz1,k,jz1)
          zwv_loc     = zwv_loc     + zwv    (iz1,k,jz1)
          zsw_loc     = zsw_loc     + zsw    (iz1,k,jz1)
          zclwc_loc   = zclwc_loc   + zclwc  (iz1,k,jz1)
          zciwc_loc   = zciwc_loc   + zciwc  (iz1,k,jz1)
          zduco2f_loc = zduco2f_loc + zduco2f(iz1,k,jz1)
          zduo3f_loc  = zduo3f_loc  + zduo3f (iz1,k,jz1)
          zaeq1_loc   = zaeq1_loc   + zaeq1  (iz1,k,jz1)
          zaeq2_loc   = zaeq2_loc   + zaeq2  (iz1,k,jz1)
          zaeq3_loc   = zaeq3_loc   + zaeq3  (iz1,k,jz1)
          zaeq4_loc   = zaeq4_loc   + zaeq4  (iz1,k,jz1)
          zaeq5_loc   = zaeq5_loc   + zaeq5  (iz1,k,jz1)
        ENDDO
      ENDDO
      izstart = izstart + nradcoarse  ! end of the last block

      zti    (izsave,k,1) = zti_loc     * gp_radcoarse_loc(ip,ib)%rfacgp
      zdpr   (izsave,k,1) = zdpr_loc    * gp_radcoarse_loc(ip,ib)%rfacgp
      zclc   (izsave,k,1) = zclc_loc    * gp_radcoarse_loc(ip,ib)%rfacgp
      zwv    (izsave,k,1) = zwv_loc     * gp_radcoarse_loc(ip,ib)%rfacgp
      zsw    (izsave,k,1) = zsw_loc     * gp_radcoarse_loc(ip,ib)%rfacgp
      zclwc  (izsave,k,1) = zclwc_loc   * gp_radcoarse_loc(ip,ib)%rfacgp
      zciwc  (izsave,k,1) = zciwc_loc   * gp_radcoarse_loc(ip,ib)%rfacgp
      zduco2f(izsave,k,1) = zduco2f_loc * gp_radcoarse_loc(ip,ib)%rfacgp
      zduo3f (izsave,k,1) = zduo3f_loc  * gp_radcoarse_loc(ip,ib)%rfacgp
      zaeq1  (izsave,k,1) = zaeq1_loc   * gp_radcoarse_loc(ip,ib)%rfacgp
      zaeq2  (izsave,k,1) = zaeq2_loc   * gp_radcoarse_loc(ip,ib)%rfacgp
      zaeq3  (izsave,k,1) = zaeq3_loc   * gp_radcoarse_loc(ip,ib)%rfacgp
      zaeq4  (izsave,k,1) = zaeq4_loc   * gp_radcoarse_loc(ip,ib)%rfacgp
      zaeq5  (izsave,k,1) = zaeq5_loc   * gp_radcoarse_loc(ip,ib)%rfacgp
    ENDDO

  ENDDO
  !$acc end parallel

! Surface
  izsave  = 0
  izstart = 0
  !$acc parallel
  !$acc loop gang vector
  DO ip = 1, ipend

    izsave   = izsave + 1

    izend = gp_radcoarse_loc(ip,ib)%ihigh - gp_radcoarse_loc(ip,ib)%ilow + 1
    jzend = gp_radcoarse_loc(ip,ib)%jhigh - gp_radcoarse_loc(ip,ib)%jlow + 1

    ztike1_loc  = 0.0_wp
    zalso_loc   = 0.0_wp
    zalth_loc   = 0.0_wp
    zapre_loc   = 0.0_wp
    zsmu0_loc   = 0.0_wp

    DO jz1 = 1, jzend
      DO iz1 = izstart+1, izstart+izend
        ztike1_loc  = ztike1_loc  + zti     (iz1,ke1,jz1)
        zalso_loc   = zalso_loc   + zalso   (iz1,    jz1)
        zalth_loc   = zalth_loc   + zalth   (iz1,    jz1)
        zapre_loc   = zapre_loc   + zapre   (iz1,    jz1)
        zsmu0_loc   = zsmu0_loc   + zsmu0_fg(iz1,    jz1)
      ENDDO
    ENDDO
    izstart = izstart + nradcoarse

    zti    (izsave,ke1,1) = ztike1_loc  * gp_radcoarse_loc(ip,ib)%rfacgp
    zalso  (izsave,    1) = zalso_loc   * gp_radcoarse_loc(ip,ib)%rfacgp
    zalth  (izsave,    1) = zalth_loc   * gp_radcoarse_loc(ip,ib)%rfacgp
    zapre  (izsave,    1) = zapre_loc   * gp_radcoarse_loc(ip,ib)%rfacgp
    zsmu0  (ip)           = zsmu0_loc   * gp_radcoarse_loc(ip,ib)%rfacgp
  ENDDO
  !$acc end parallel

!------------------------------------------------------------------------------

  !$acc end data

END SUBROUTINE average_to_coarse_grid

!==============================================================================

!-------------------------------------------------------------------------------
! End of the subroutine
!-------------------------------------------------------------------------------

END SUBROUTINE radiation_organize

!==============================================================================
!==============================================================================
!+ Module procedure in "Radiation" 
!------------------------------------------------------------------------------

SUBROUTINE radiation_average (lzradstep, ntl)

!------------------------------------------------------------------------------
!
! Description:
!
! When running the radiation on a coarser grid, several variables need a
! postprocessing:
!   - thbs (thermal radiation at the ground)
!     the surface outgoing thermal radiation is computed with an average value 
!     of the ground temperature (stored in ztg_rc) and subtracted from thbs.
!     It is added again, but computed with the values from the full COSMO grid.
!     If lradf_avg is set, the values of thbs are filtered before adding the
!     surface outgoing thermal radiation again.
!
!   - sobs (solar radiation at the ground) and
!     pabs (photosynthetic active radiation at the ground)
!     the reflected solar radiation is computed with an average value of the 
!     solar surface albedo (stored in alb_rad_rc) and subtracted from pabs and
!     sobs. It is added again, but computed with the values from the full
!     COSMO grid.
!     If lradf_avg is set, the values of sobs and pabs are filtered before 
!     adding the reflected solar radiation again.
!
!  lradf_avg=.TRUE.: 
!    if this switch is set to .TRUE., the radiative increments are filtered 
!    using a discrete filter (analogous to lconf_avg)
!     
! These computations must be done on the full COSMO grid!
!
!------------------------------------------------------------------------------

LOGICAL, INTENT(IN) :: lzradstep    ! true if full radiation timestep, otherwise
                                    ! only recalculation of sohr, sobs, and sobt,
                                    ! based on actual zenith angle

INTEGER, INTENT(IN) :: ntl          ! time level to be considered


! Local variables

INTEGER        ::  &
  i, j, k, is1, ie1, js1, je1

REAL(KIND=wp)  ::  &
  zalbfak            ! albedo correction factor !T.R.

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 0: Initialization of indices
!------------------------------------------------------------------------------

#ifdef _OPENACC
!XL: not supported, this routine needs to be ported to GPU
CALL model_abort(my_cart_id, 667, "Error: radiation_average not ported to GPU", 'radiation_average')
#endif

! Computations have to be done for the boundary lines at the border of the
! total domain, and at least for one boundary line at the border to neighboring
! domains (for averaging later on)

is1 = MIN (istartpar, istart-1)
ie1 = MAX (iendpar  , iend  +1)
js1 = MIN (jstartpar, jstart-1)
je1 = MAX (jendpar  , jend  +1)

!------------------------------------------------------------------------------
! Section 1: Averaging for thermal variables in a full radiation step
!------------------------------------------------------------------------------

IF (lzradstep) THEN
  DO j = js1, je1
    DO i =  is1, ie1
      thbs(i,j)  = thbs(i,j) + sigma*(1.0_wp - ctalb)*(ztg_rc(i,j)**4)
!                                                     EXP( 4.0_wp * LOG(ztg_rc(i,j)))
    ENDDO
  ENDDO
  
  IF (lradf_avg) THEN

    zwork3d(:,:,:) = thhr(:,:,:)

    DO k = 1, ke
      DO j = jstart, jend
        DO i = istart, iend
          thhr(i,j,k) =                                                                                     &
            ( zcent *  zwork3d(i  ,j  ,k)                                                                   &
            + zside * (zwork3d(i-1,j  ,k) + zwork3d(i+1,j  ,k) + zwork3d(i  ,j-1,k) + zwork3d(i  ,j+1,k) )  &
            + zedge * (zwork3d(i-1,j-1,k) + zwork3d(i+1,j-1,k) + zwork3d(i-1,j+1,k) + zwork3d(i+1,j+1,k) ) )
        ENDDO
      ENDDO
    ENDDO !k

    zwork3d(:,:,1) = thbs        (:,:)
    zwork3d(:,:,2) = thbt        (:,:)
    zwork3d(:,:,3) = sodwddm     (:,:)
    zwork3d(:,:,4) = lwd_s       (:,:)
    zwork3d(:,:,5) = lwu_s       (:,:)

!US these fields belong to a development which did never make it into the official
!   version. But if this development will be activated again (see also subroutine
!   radiation_average in radiation_interface) they are needed as global fields
!   zwork3d(:,:,6) = flpar_s_dir (:,:)
!   zwork3d(:,:,7) = flpar_s_difd(:,:)
!   zwork3d(:,:,8) = flpar_s_difu(:,:)


    DO j=jstart,jend
      DO i = istart,iend

        thbs(i,j)         =                                                                               &
          ( zcent *  zwork3d(i  ,j  ,1)                                                                   &
          + zside * (zwork3d(i-1,j  ,1) + zwork3d(i+1,j  ,1) + zwork3d(i  ,j-1,1) + zwork3d(i  ,j+1,1) )  &
          + zedge * (zwork3d(i-1,j-1,1) + zwork3d(i+1,j-1,1) + zwork3d(i-1,j+1,1) + zwork3d(i+1,j+1,1) ) )
        thbt(i,j)         =                                                                               &
          ( zcent *  zwork3d(i  ,j  ,2)                                                                   &
          + zside * (zwork3d(i-1,j  ,2) + zwork3d(i+1,j  ,2) + zwork3d(i  ,j-1,2) + zwork3d(i  ,j+1,2) )  &
          + zedge * (zwork3d(i-1,j-1,2) + zwork3d(i+1,j-1,2) + zwork3d(i-1,j+1,2) + zwork3d(i+1,j+1,2) ) )
        sodwddm(i,j)      =                                                                               &
          ( zcent *  zwork3d(i  ,j  ,3)                                                                   &
          + zside * (zwork3d(i-1,j  ,3) + zwork3d(i+1,j  ,3) + zwork3d(i  ,j-1,3) + zwork3d(i  ,j+1,3) )  &
          + zedge * (zwork3d(i-1,j-1,3) + zwork3d(i+1,j-1,3) + zwork3d(i-1,j+1,3) + zwork3d(i+1,j+1,3) ) )
        lwd_s(i,j)        =                                                                               &
          ( zcent *  zwork3d(i  ,j  ,4)                                                                   &
          + zside * (zwork3d(i-1,j  ,4) + zwork3d(i+1,j  ,4) + zwork3d(i  ,j-1,4) + zwork3d(i  ,j+1,4) )  &
          + zedge * (zwork3d(i-1,j-1,4) + zwork3d(i+1,j-1,4) + zwork3d(i-1,j+1,4) + zwork3d(i+1,j+1,4) ) )
        lwu_s(i,j)        =                                                                               &
          ( zcent *  zwork3d(i  ,j  ,5)                                                                   &
          + zside * (zwork3d(i-1,j  ,5) + zwork3d(i+1,j  ,5) + zwork3d(i  ,j-1,5) + zwork3d(i  ,j+1,5) )  &
          + zedge * (zwork3d(i-1,j-1,5) + zwork3d(i+1,j-1,5) + zwork3d(i-1,j+1,5) + zwork3d(i+1,j+1,5) ) )
!       flpar_s_dir(i,j)  =                                                                               &
!         ( zcent *  zwork3d(i  ,j  ,6)                                                                   &
!         + zside * (zwork3d(i-1,j  ,6) + zwork3d(i+1,j  ,6) + zwork3d(i  ,j-1,6) + zwork3d(i  ,j+1,6) )  &
!         + zedge * (zwork3d(i-1,j-1,6) + zwork3d(i+1,j-1,6) + zwork3d(i-1,j+1,6) + zwork3d(i+1,j+1,6) ) )
!       flpar_s_difd(i,j) =                                                                               &
!         ( zcent *  zwork3d(i  ,j  ,7)                                                                   &
!         + zside * (zwork3d(i-1,j  ,7) + zwork3d(i+1,j  ,7) + zwork3d(i  ,j-1,7) + zwork3d(i  ,j+1,7) )  &
!         + zedge * (zwork3d(i-1,j-1,7) + zwork3d(i+1,j-1,7) + zwork3d(i-1,j+1,7) + zwork3d(i+1,j+1,7) ) )
!       flpar_s_difu(i,j) =                                                                               &
!         ( zcent *  zwork3d(i  ,j  ,8)                                                                   &
!         + zside * (zwork3d(i-1,j  ,8) + zwork3d(i+1,j  ,8) + zwork3d(i  ,j-1,8) + zwork3d(i  ,j+1,8) )  &
!         + zedge * (zwork3d(i-1,j-1,8) + zwork3d(i+1,j-1,8) + zwork3d(i-1,j+1,8) + zwork3d(i+1,j+1,8) ) )
      ENDDO
    ENDDO

  ENDIF !lradf_avg

! DO j = js1, je1
!   DO i = is1, ie1
  DO j = jstartpar, jendpar
    DO i = istartpar, iendpar
      thbs (i,j) = thbs(i,j) - sigma*(1.0_wp - ctalb)*(t_g(i,j,ntl)**4)
!                                                      EXP( 4.0_wp * LOG(t_g(i,j)))
    ENDDO
  ENDDO

ENDIF  ! lzradstep

!------------------------------------------------------------------------------
! Section 2: Averaging for solar variables in all steps
!------------------------------------------------------------------------------

DO j = js1, je1
  DO i =  is1, ie1
    ! this was eliminated in testsuite 3.4
    ! but keep it for the moment to reproduce results with Version 3.22
    zalbfak      = 1._wp/(1._wp-alb_rad_rc(i,j))
    sobs (i,j)  = sobs(i,j) * zalbfak
    pabs (i,j)  = pabs(i,j) * zalbfak
  ENDDO
ENDDO

IF (lradf_avg) THEN

  zwork3d (:,:,:)   = sohr (:,:,:)

  DO k = 1, ke
    DO j = jstart, jend
      DO i = istart, iend
        sohr(i,j,k) =                                                                                     &
          ( zcent *  zwork3d(i  ,j  ,k)                                                                   &
          + zside * (zwork3d(i-1,j  ,k) + zwork3d(i+1,j  ,k) + zwork3d(i  ,j-1,k) + zwork3d(i  ,j+1,k) )  &
          + zedge * (zwork3d(i-1,j-1,k) + zwork3d(i+1,j-1,k) + zwork3d(i-1,j+1,k) + zwork3d(i+1,j+1,k) ) )
      ENDDO
    ENDDO
  ENDDO

  zwork3d (:,:,1)   = sobs       (:,:)
  zwork3d (:,:,2)   = sobt       (:,:)
  zwork3d (:,:,3)   = pabs       (:,:)
  zwork3d (:,:,4)   = swdir_s    (:,:)
  zwork3d (:,:,5)   = swdifd_s   (:,:)
  zwork3d (:,:,6)   = swdifu_s   (:,:)

  DO j = jstart, jend
    DO i = istart, iend
      sobs(i,j)         =                                                                               &
        ( zcent *  zwork3d(i  ,j  ,1)                                                                   &
        + zside * (zwork3d(i-1,j  ,1) + zwork3d(i+1,j  ,1) + zwork3d(i  ,j-1,1) + zwork3d(i  ,j+1,1) )  &
        + zedge * (zwork3d(i-1,j-1,1) + zwork3d(i+1,j-1,1) + zwork3d(i-1,j+1,1) + zwork3d(i+1,j+1,1) ) )
      sobt(i,j)         =                                                                               &
        ( zcent *  zwork3d(i  ,j  ,2)                                                                   &
        + zside * (zwork3d(i-1,j  ,2) + zwork3d(i+1,j  ,2) + zwork3d(i  ,j-1,2) + zwork3d(i  ,j+1,2) )  &
        + zedge * (zwork3d(i-1,j-1,2) + zwork3d(i+1,j-1,2) + zwork3d(i-1,j+1,2) + zwork3d(i+1,j+1,2) ) )
      pabs(i,j)         =                                                                               &
        ( zcent *  zwork3d(i  ,j  ,3)                                                                   &
        + zside * (zwork3d(i-1,j  ,3) + zwork3d(i+1,j  ,3) + zwork3d(i  ,j-1,3) + zwork3d(i  ,j+1,3) )  &
        + zedge * (zwork3d(i-1,j-1,3) + zwork3d(i+1,j-1,3) + zwork3d(i-1,j+1,3) + zwork3d(i+1,j+1,3) ) )
      swdir_s(i,j)      =                                                                               &
        ( zcent *  zwork3d(i  ,j  ,4)                                                                   &
        + zside * (zwork3d(i-1,j  ,4) + zwork3d(i+1,j  ,4) + zwork3d(i  ,j-1,4) + zwork3d(i  ,j+1,4) )  &
        + zedge * (zwork3d(i-1,j-1,4) + zwork3d(i+1,j-1,4) + zwork3d(i-1,j+1,4) + zwork3d(i+1,j+1,4) ) )
      swdifd_s(i,j)     =                                                                               &
        ( zcent *  zwork3d(i  ,j  ,5)                                                                   &
        + zside * (zwork3d(i-1,j  ,5) + zwork3d(i+1,j  ,5) + zwork3d(i  ,j-1,5) + zwork3d(i  ,j+1,5) )  &
        + zedge * (zwork3d(i-1,j-1,5) + zwork3d(i+1,j-1,5) + zwork3d(i-1,j+1,5) + zwork3d(i+1,j+1,5) ) )
      swdifu_s(i,j)     =                                                                               &
        ( zcent *  zwork3d(i  ,j  ,6)                                                                   &
        + zside * (zwork3d(i-1,j  ,6) + zwork3d(i+1,j  ,6) + zwork3d(i  ,j-1,6) + zwork3d(i  ,j+1,6) )  &
        + zedge * (zwork3d(i-1,j-1,6) + zwork3d(i+1,j-1,6) + zwork3d(i-1,j+1,6) + zwork3d(i+1,j+1,6) ) )
    ENDDO
  ENDDO

ENDIF  ! lradf_avg

DO j = js1, je1
  DO i =  is1, ie1
    ! such it has been tested in testsuite 3.4
    ! but keep it for the moment to reproduce results with Version 3.22
    zalbfak = (1.0_wp-alb_rad(i,j))
    sobs (i,j) = sobs(i,j) * zalbfak
    pabs (i,j) = pabs(i,j) * zalbfak
    ! And this seems to be the better version
    ! swdifu_s(i,j) = ( swdir_s(i,j)     + swdifd_s(i,j) ) * alb_rad(i,j)
    ! sobs(i,j)     =   swdir_s(i,j)     + swdifd_s(i,j)     - swdifu_s(i,j)
    ! pabs(i,j)     =   flpar_s_dir(i,j) + flpar_s_difd(i,j) - flpar_s_difu(i,j)
  ENDDO
ENDDO

!-------------------------------------------------------------------------------
! End of the subroutine
!-------------------------------------------------------------------------------

END SUBROUTINE radiation_average

!==============================================================================
!==============================================================================

SUBROUTINE radiation_in_wkarr_alloc (lradstep, nproma_rad, ke, nradcoarse, ierrloc)

  LOGICAL, INTENT(IN)  :: &
    lradstep                     ! if .TRUE., allocate all variables for full radiation step

  INTEGER, INTENT(IN)  :: &
    nproma_rad, ke, nradcoarse

  INTEGER, INTENT(OUT) :: &
    ierrloc

  INTEGER :: izstat, izdim

  ierrloc = 0
  izstat  = 0
  izdim   = nproma_rad * nradcoarse

  IF (lradstep) THEN
    ! Input of fesft (output of routines from radiation_utilities)
    ALLOCATE (zduco2f(izdim, ke, nradcoarse),                         &
              zduo3f (izdim, ke, nradcoarse),                         &
              zaeq1  (izdim, ke, nradcoarse),                         &
              zaeq2  (izdim, ke, nradcoarse),                         &
              zaeq3  (izdim, ke, nradcoarse),                         &
              zaeq4  (izdim, ke, nradcoarse),                         &
              zaeq5  (izdim, ke, nradcoarse),          STAT=izstat)
    IF (izstat /= 0) THEN
      PRINT *, '*** Error allocating local block variables ***'
      ierrloc = 2771
    ENDIF

    ALLOCATE (zsw    (izdim, ke, nradcoarse),                         &
              zwv    (izdim, ke, nradcoarse),                         &
              zclwc  (izdim, ke, nradcoarse),                         &
              zciwc  (izdim, ke, nradcoarse),                         &
              zclc   (izdim, ke, nradcoarse),          STAT=izstat)
    IF (izstat /= 0) THEN
      PRINT *, '*** Error allocating local block variables ***'
      ierrloc = 2772
    ENDIF

    ALLOCATE (zalso  (izdim, nradcoarse),                             &
              zalth  (izdim, nradcoarse),              STAT=izstat)
    IF (izstat /= 0) THEN
      PRINT *, '*** Error allocating local block variables ***'
      ierrloc = 2773
    ENDIF

    ALLOCATE (zti     (izdim, ke1, nradcoarse),                       &
              zdpr    (izdim, ke,  nradcoarse),                       &
              zapre   (izdim, nradcoarse),                            &
              zskyview(izdim, nradcoarse),                            &
              zsmu0   (izdim),                                        &
!US GPU complains  zsmu0   (nproma_rad),
                                                        STAT=izstat)
    IF (nradcoarse > 1) THEN
      ALLOCATE (ztg_rc  (ie,je),                        STAT=izstat)
      IF (izstat /= 0) THEN
        PRINT *, '*** Error allocating local work  variables ***'
        ierrloc = 2774
      ENDIF
    ENDIF

    ! allocate all the above fesft input arrays for the device
    !$acc enter data                                                    &
    !$acc create ( zduco2f, zduo3f, zaeq1, zaeq2, zaeq3, zaeq4, zaeq5 ) &
    !$acc create ( zsw, zwv, zclwc, zciwc, zclc, zalso, zalth )         &
    !$acc create ( zti, zdpr, zapre, zskyview, zsmu0, ztg_rc )

    ! Output of fesft
    ALLOCATE (zflt         (izdim, ke1),                              &
              zfls         (izdim, ke1),                              &
              zflsdir      (izdim, ke1),                              &
              zflt_s       (izdim),                                   &
              zfls_s       (izdim),                                   &
              zflsp        (izdim),                                   &
              zflsd        (izdim),                                   &
              zflsu        (izdim),                                   &
              zfltd        (izdim),                                   &
              zfltu        (izdim),                                   &
              zflpar_s     (izdim),                                   &
              zflpar_s_dir (izdim),                                   &
              zflpar_s_difd(izdim),                                   &
              zflpar_s_difu(izdim),                                   &
                                                        STAT=izstat)
    IF (izstat /= 0) THEN
      PRINT *, '*** Error allocating local block variables ***'
      ierrloc = 2775
    ENDIF

    ! work arrays for radiation_utilities
    ALLOCATE (zwork3d1(izdim, ke,  nradcoarse),                       &
              zwork3d2(izdim, ke,  nradcoarse),                       &
              zwork3d3(izdim, ke,  nradcoarse),                       &
              zwork2d1(izdim,      nradcoarse),                       &
              zwork2d2(izdim,      nradcoarse),                       &
              zwork2d3(izdim,      nradcoarse),                       &
              zwork2d4(izdim,      nradcoarse),                       &
                                                        STAT=izstat)
    IF (izstat /= 0) THEN
      PRINT *, '*** Error allocating local block variables ***'
      ierrloc = 2776
    ENDIF

    ! allocate all the above fesft output and work arrays for the device
    !$acc enter data                                                  &
    !$acc create ( zflt, zfls, zflsdir, zflt_s, zfls_s, zflsp, zflsd) &
    !$acc create ( zflsu, zfltd, zfltu, zflpar_s, zflpar_s_dir      ) &
    !$acc create ( zflpar_s_difd, zflpar_s_difu                     ) &
    !$acc create ( zwork3d1, zwork3d2, zwork3d3                     ) &
    !$acc create ( zwork2d1, zwork2d2, zwork2d3, zwork2d4           )

  ENDIF

  ! the following variables are also needed during zenith angle update:
  ALLOCATE (zsmu0_fg(izdim, nradcoarse),                       &
            zwork2d5(izdim, nradcoarse),                       &
            zwork2d6(izdim, nradcoarse),                       &
            zwork2d7(izdim, nradcoarse),                       &
            zwork2d8(izdim, nradcoarse),                       &
            nzwork  (izdim, nradcoarse),                       &
                                                        STAT=izstat)
  IF (nradcoarse > 1) THEN
    ALLOCATE (zwork3d(ie,je,ke),                        STAT=izstat)
  ENDIF
  IF (izstat /= 0) THEN
    PRINT *, '*** Error allocating local solar variables ***'
    ierrloc = 2777
  ENDIF

  ! allocate the arrays for the device
  !$acc enter data                                                  &
  !$acc create ( zsmu0_fg, zwork2d5, zwork2d6, zwork2d7, zwork2d8 ) &
  !$acc create ( nzwork, zwork3d )

  !XL: qx_locrad are not associated, if radiation_in_wkarr_alloc is called 
  !    in the initialization phase
  !NOacc create ( nzwork, zwork3d, qc_locrad, qv_locrad, qi_locrad )

END SUBROUTINE radiation_in_wkarr_alloc

!==============================================================================
!==============================================================================

SUBROUTINE radiation_in_wkarr_dealloc(lradstep, ierrloc)

  LOGICAL, INTENT(IN)  :: &
    lradstep                     ! if .TRUE., allocate all variables for full radiation step

  INTEGER, INTENT(OUT) :: ierrloc

  INTEGER :: izstat

  ierrloc = 0
  izstat  = 0

  IF (lradstep) THEN
    DEALLOCATE (zduco2f, zduo3f , zaeq1  , zaeq2  , zaeq3  , zaeq4  , zaeq5, STAT=izstat)
    IF (izstat /= 0) THEN
      PRINT *, '*** Error deallocating local block variables ***'
      ierrloc = 2771
    ENDIF

    DEALLOCATE (zsw, zwv, zclwc, zciwc, zclc,                                STAT=izstat)
    IF (izstat /= 0) THEN
      PRINT *, '*** Error deallocating local block variables ***'
      ierrloc = 2772
    ENDIF

    DEALLOCATE (zalso, zalth,                                                STAT=izstat)
    IF (izstat /= 0) THEN
      PRINT *, '*** Error deallocating local block variables ***'
      ierrloc = 2773
    ENDIF

    DEALLOCATE (zti, zdpr, zapre, zsmu0, zskyview,                           STAT=izstat)
    IF (nradcoarse > 1) THEN
      DEALLOCATE (ztg_rc,                                                    STAT=izstat)
    ENDIF
    IF (izstat /= 0) THEN
      PRINT *, '*** Error deallocating local work  variables ***'
      ierrloc = 2774
    ENDIF

    ! deallocate all the above fesft input arrays for the device
    !$acc exit data                                                     &
    !$acc delete ( zduco2f, zduo3f, zaeq1, zaeq2, zaeq3, zaeq4, zaeq5 ) &
    !$acc delete ( zsw, zwv, zclwc, zciwc, zclc, zalso, zalth )         &
    !$acc delete ( zti, zdpr, zapre, zskyview, zsmu0, ztg_rc )

    DEALLOCATE (zflt, zfls, zflsdir, zfls_s, zflt_s, zflsp, zflsd, zflsu,            &
                zfltd, zfltu, zflpar_s, zflpar_s_dir, zflpar_s_difd,                 &
                zflpar_s_difu,                                               STAT=izstat)
    IF (izstat /= 0) THEN
      PRINT *, '*** Error deallocating local block variables ***'
      ierrloc = 2775
    ENDIF

    DEALLOCATE (zwork3d1, zwork3d2, zwork3d3, zwork2d1, zwork2d2, zwork2d3,          &
                zwork2d4,                                                    STAT=izstat)
    IF (izstat /= 0) THEN
      PRINT *, '*** Error deallocating local block variables ***'
      ierrloc = 2776
    ENDIF
    ! allocate all the above fesft output and work arrays for the device
    !$acc exit data                                                   &
    !$acc delete ( zwork3d1, zwork3d2, zwork3d3                     ) &
    !$acc delete ( zwork2d1, zwork2d2, zwork2d3, zwork2d4           )

  ENDIF

  ! the following variables are also needed during zenith angle update:
  DEALLOCATE (zsmu0_fg, zwork2d5, zwork2d6, zwork2d7, zwork2d8, nzwork, STAT=izstat)
  IF (nradcoarse > 1) THEN
    DEALLOCATE (zwork3d, STAT=izstat)
  ENDIF
  IF (izstat /= 0) THEN
    PRINT *, '*** Error deallocating local solar variables ***'
    ierrloc = 2777
  ENDIF
  
  ! allocate the arrays for the device
  !$acc exit data                                                   &
  !$acc delete ( zsmu0_fg, zwork2d5, zwork2d6, zwork2d7, zwork2d8, nzwork, zwork3d )

END SUBROUTINE radiation_in_wkarr_dealloc

!==============================================================================

END MODULE radiation_interface

!==============================================================================
