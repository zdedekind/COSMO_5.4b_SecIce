!+ Source module for computing local info. on conventional single-level obs
!-------------------------------------------------------------------------------

MODULE src_sing_local

!-------------------------------------------------------------------------------
!
! Description:
!   This module 'src_sing_local' computes all the 'local' information (i.e.
!   observation increments and further parameters on the observations and their 
!   location) on the local single-level data (upper-air single-level, surface-
!   level, and surface pressure data), which is required for the spreading of
!   the observational information later on in the nudging.
!   Also, lists of all reports influencing the current timestep and other
!   preliminaries are done for all types of conventional data (including
!   multi-level data).
!   Specific tasks include:
!    - organizing the application of the (forward and partly inverse) 
!      observation operator (e.g. by vertical interpolation of model values, or
!      observations to the lowest model level (e.g. extrapolation of surface-
!      level obs to the lowest model level using an old surface transfer scheme;
!      or extrapolation of station pressure obs to the model orography)
!    - organizing the threshold quality control,
!    - performing a spatial consistency check of surface pressure obs.
!   Major procedures related to the observation operators and quality control are
!   out-sourced into modules 'src_obs_operator_conv' and 'src_obs_qc_conv' after
!   V5_0.
!
! Method:    
!   This module contains the following procedures:
!    - local_info_aux    : auxiliaries for computing local information on obs
!    - local_sort_reports: produce sorted lists of currently used reports
!    - ps_local_info     : organize obs increments + QC for surface pressure obs
!    - ps_spatial_check  : spatial consistency check of surface pressure obs
!    - surf_local_info   : organize obs increments and QC for surface-level obs
!    - surf_obs          : extrapolate screen-level obs to lowest model level
!    - upair_local_info  : organize obs inc. + QC for upper-air single-level obs
!   All routines are called by 'organize_nudging' (module 'src_obs_use_org.f90')
!   except for 'surf_obs' (called by 'surf_local_info').
!
!   This module also contains elemental functions, formerly statement functions:
!    - fpvsw       : Magnus formula for water: saturation vapour pressure from T
!    - fq2pv       : vapour pressure from specific humidity and air pressure
!    - rmod        : MOD function for positive reals
!    - ibit1       : returns 1 bit at given bit position of given integer word
!
!   Note: This module does not contain any communication between PE's,
!         except for 'model_abort' at retrieving tracer model fields.
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! History:
! Version    Date       Name
! -------    ----       ----
! 1.13       1998/10/22 Christoph Schraff
!  Initial release
! 1.15       1998/11/02 Christoph Schraff
!  Global declaration of allocatable arrays moved to module 'data_nudge_local'.
! 1.19       1998/12/11 Christoph Schraff
!  Setting quality control flags, also for passive data (verification mode).
!  Extrapolation of 2m-humidity data to lowest model level replaced by old
!  scheme (assuming constant relative humidity). ANSI violations removed.
! 1.20       1999/01/07 Guenther Doms
!  Renaming of some global variables
! 1.27       1999/03/29 Christoph Schraff
!  Threshold quality control for surface-level and SYNOP pressure data performed
!  at the station level instead of the lowest model level. Revised ODR format.
! 1.31       1999/07/01 Christoph Schraff
!  Bug correction by initialization of a variable (hint by J.-M. Bettems).
!  Quantities related to MPI communicator 'icomm_world' removed.
! 1.36       2000/02/24 Christoph Schraff
!  Optional storage of observation increments at observation points for VOF.
!  Bug correction in threshold quality control of 2m humidity observations.
! 1.37       2000/03/24 Guenther Doms
!  Minimum limitations for 'zchdcm' and 'zsqcm' introduced.
! 1.38       2000/04/06 Christoph Schraff
!  Surface-level obs increments optionally derived by use of new surface layer
!  parameterization. Inclusion of observed pressure tendency into pressure QC
!  thresholds. Initialisation of height error correlation matrix and QC flags. 
!  Optional storage of obs increments for VOF.
! 1.40       2000/05/23 Christoph Schraff
!  Only active height data used for derivation of TEMP 'surface' pressure.
!  Reduced (by 0.5) weight to pressure tendency for pressure QC thresholds.
! 2.4        2001/01/29 Christoph Schraff
!  Addtion of spatial consistency check for surface pressure observations.
!  Rejection of aircraft temperature, if wind is rejected, and vice versa.
!  Additional factors to QC thresholds for aircrafts and DRIBUs following ECMWF.
! 2.5        2001/06/01 Christoph Schraff
!  Correction of lethal bug at call for computation of surface-level increments
!  (and for quality control), if new analysis increments are computed at each 
!  timestep. Memory deallocation of ODR. Savety test at array allocations.
! 2.7        2001/06/26 Michael Buchhold
!  Bug corrrection. Memory deallocation of ODR only when no surface analysis
!  is carried out.
! 2.13       2002/01/18 Christoph Schraff
!  Correction: Array elements replaced by reals for intent-in subr arguments.
! 2.18       2002/07/16 Ulrich Schaettler
!  lowered index of the lowest horizontal main model level to kk (kk-1 before)
!  (by Guy deMorsier)
! 2.19       2002/10/24 Christoph Schraff + Michael Buchhold
!  - Options for assimilation of 2m humidity data:
!    additional stability-dependent vertical weights (also for 10m wind),
!    quality factor dependent on T-2m observation increments, obs. increments
!    as differences of specific humidity.      (Christoph Schraff)
!  - Shorten temporal range of influence of frequent observations as wind-
!    profiler/rass or VAD-winds                (Michael Buchhold)
! 3.3        2003/04/22 Maria Tomassini + Christoph Schraff
!  Administating and sorting reports for assimilation of GPS-derived IWV.
! 3.6        2003/12/11 Christoph Schraff
!  Introduction of SATOB reports. Minor bug corrections: no T-2m at timestep 0;
!  no integer of infinite ('rmdi') station height.
! 3.12       2004/09/15 Christoph Schraff
!  Extension to (prepare to) include assimilation of satellite retrievals.
!  Correction to allow for verification of all passsive reports.
!  Optional reduction for temporal enhancement of quality control thresholds.
! 3.15       2005/03/03 Ulrich Schaettler
!  Replaced FLOAT by REAL
! 3.16       2005/07/22 Guy DeMorsier
!  Bug correction in SR surf_obs.
! 3.18       2006/03/03 Christoph Schraff
!  Further reduction for temporal enhancement of humidity QC thresholds.
!  Preparation for use of real-data 1DVar satellite retrievals (MSG, NOAA15-18).
!  Adaptations for spatial consistency check of IWV. Flushing of output files.
!  To allow for reproducibility, reset of xoisu, xoiua if report is not taken.
! 3.21       2006/12/04 Ulrich Schaettler
!  Put statement functions which are not used to a comment
! V3_23        2007/03/30 Ulrich Schaettler
!  Introduced ldump_ascii for flushing the ASCII files
!  Moved 'akt' to MODULE data_turbulence (M. Raschendorfer)
! V3_24        2007/04/26 Christoph Schraff
!  Quality control adapted for aircraft humidity observations.
! V4_5         2008/09/10 Christoph Schraff
!  QC enforced the first time an obs is used even if 'lobprcs' is not used (when
!  obs are read from NetCDF files), and element 'nhqcfw' set then from -1 to -2.
!  'o??vip' arrays replaced by 'nbsvi?' or 'nhtvi?' elements of ODR.
!  QC (-flag) lines now looks similar for ps-, single-level and multi-level obs.
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_22        2012/01/31 Christoph Schraff
!  - Observations with the same station id as the quality controlled observation
!    are not used any more in the spatial consistency check of surface pressure.
!  - Bug fix: Surface pressure obs with rejection flag set by pre-processing are
!             also threshold QC flagged now.
!  - Bug fix: to avoid array bound violations: maxmlo replaced by maxmll as size
!             of 'imladm' etc., and checks introduced related to 'nexce??'.
!  - Introduction of non-zero QC threshold for upper-air single-level humidity.
!  - Some variables moved from 'data_nudge_all' to 'data_obs_lib_cosmo', use of
!    indices of SOR (simulated obs record). Variable 'kobtysu' introduced.
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer
!  Replaced qx-variables by using them from the tracer module
! V4_27        2013/03/19 Astrid Kerkweg
!  MESSy interface introduced
! V4_28        2013/07/12 Christoph Schraff, Ulrich Schaettler
!  - Introduction of quality control (individual QC and spatial consistency
!    check 'SCC') of surface pressure data against lateral boundary fields.
!  - Bug fixes for near-surface pressure increments from multi-level reports:
!              QC-rejected increments written to 'dmlhed';  increments for SCC
!              computed also for obs rejected by individual QC;  in SCC, modulus
!              of element 'nhtvip' deployed to determine increment in QC check.
!  - Major procedures related purely to the observation operators and quality 
!    control (except for spatial consistency check of surface pressure) are
!    out-sourced into modules 'src_obs_operator_conv' and 'src_obs_qc_conv'.
!    New, slim interfaces:  Apart from constants, all input and all output
!    information for these out-sourced procedures is via subroutine arguments.
!  - Statement functions replaced by elemental or intrinsic functions.
!  - Improved comments.
!  Use parameters for vertical grid from module vgrid_refatm_utils (US)
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler, Christoph Schraff
!  Unification of MESSy interfaces and COSMO Tracer structure
!  For the COSMO-Model only use vcoord from vgrid_refatm_utils
!  Bug fix: Allocation of 'ssgbdy' etc made unconditional. (CS)
! V5_1         2014-11-28 Christoph Schraff, Oliver Fuhrer
!  Bug fix: If condition introduced to avoid unnecessary caution messages (CS)
!  Processing of Mode-S aircraft obs introduced (CS)
!  Replaced ireals by wp (working precision) (OF)
! V5_4         2016-03-10 Christoph Schraff
!  Variables related to the AOF interface removed.
! V5_4a        2016-05-10 Matthias Raschendorfer
!  Use new module turb_data instead of data_turbulence
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
!
!-------------------------------------------------------------------------------

USE data_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!-------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------

    ie_tot,       & ! number of grid points in zonal direction
    je_tot,       & ! number of grid points in meridional direction
!   ie,           & ! number of grid points in zonal direction
!   je,           & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction

! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the 
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from 
!    the other ones because of the use of the staggered Arakawa-B-grid.
!    
!   zonal direction
    istart,       & ! start index for the forecast of w, t, qv, qc and pp
    iend,         & ! end index for the forecast of w, t, qv, qc and pp

!   meridional direction
    jstart,       & ! start index for the forecast of w, t, qv, qc and pp
    jend,         & ! end index for the forecast of w, t, qv, qc and pp

! 5. variables for the time discretization and related variables
! --------------------------------------------------------------

    dt,           & ! long time-step
    dtdeh,        & ! dt / 3600 seconds

! 8. Organizational variables to handle the COSMO humidity tracers
! ----------------------------------------------------------------
    idt_qv,  idt_qc

! end of data_modelconfig
!-------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 2. physical constants and related variables
! -------------------------------------------

    r_d,          & ! gas constant for dry air
    r_v,          & ! gas constant for water vapor
    rdv,          & ! r_d / r_v
    o_m_rdv,      & ! 1 - r_d/r_v
    rvd_m_o,      & ! r_v/r_d - 1
    cp_d,         & ! specific heat of dry air at constant pressure
    cpdr,         & ! 1 / cp_d
    rdocp,        & ! r_d / cp_d
    g,            & ! acceleration due to gravity

! 3. constants for parametrizations
! ---------------------------------

    b1,           & ! variables for computing the saturation vapour pressure
    b2w,          & ! over water (w) and ice (i)
    b3,           & !               -- " --
    b4w             !               -- " --

! end of data_constants

!-------------------------------------------------------------------------------

USE turb_data       , ONLY :   &

    akt             ! von Karman-constant

! end of turb_data
!-------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------

    rho0       ,    & ! reference density at the full model levels    (kg/m3)
    dp0        ,    & ! reference pressure thickness of layers        ( Pa)
    p0         ,    & ! reference pressure at main levels             ( Pa)
    hhl        ,    & ! geometical height of half model levels        ( m )

! 2. external parameter fields                                        (unit)
! ----------------------------

    gz0        ,    & ! roughness length * g                          (m2/s2)

! 3. prognostic variables                                             (unit)
! -----------------------

    u          ,    & ! zonal wind speed                              ( m/s )
    v          ,    & ! meridional wind speed                         ( m/s )
    t          ,    & ! temperature                                   (  k  )
    pp         ,    & ! deviation from the reference pressure         ( pa  )

! 5. fields for surface values and soil model variables               (unit )
! -----------------------------------------------------

    ps         ,    & ! surface pressure                              ( pa  )
!   ts         ,    & ! temperature of the snow-surface               (  k  )
!   tb         ,    & ! temperature of the ground surface             (  k  )
    t_g        ,    & ! weighted surface temperature                  (  k  )
    qv_s       ,    & ! specific water vapor content on the surface   (kg/kg)

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------

!   turbulent coefficients at the surface
    tcm        ,    & ! turbulent diffusion coefficients for momentum   --
    tch        ,    & ! turbulent diffusion coefficients for heat       --
                      ! and moisture

!   fields from the radiation scheme
    clch       ,    & ! cloud cover with high clouds                    -- 
    clcm       ,    & ! cloud cover with medium clouds                  --
    clcl       ,    & ! cloud cover with low clouds                     --
    clct       ,    & ! total cloud cover                               --
    clc_sgs    ,    & ! subgrid-scale stratiform cloud cover            --
    clc_con    ,    & ! cloud cover due to convection                   --     

!   fields of the precipitation
    qrs        ,    & ! precipitation water (water loading)           (kg/kg)

! 7. fields for model output and diagnostics                          (unit )
! ------------------------------------------

    t_2m       ,    & ! temperature in 2m                             (  K  )
    qv_2m      ,    & ! specific water vapor content in 2m            (kg/kg)
    td_2m      ,    & ! dew-point in 2m                               (  K  )
    u_10m      ,    & ! zonal wind in 10m                             ( m/s )
    v_10m      ,    & ! meridional wind in 10m                        ( m/s )

! 8. fields for the boundary values                                   (unit )
! ---------------------------------

    pp_bd             ! boundary field for pp                         (  pa )

! end of data_fields

!-------------------------------------------------------------------------------

USE src_tracer       , ONLY : trcr_get, trcr_errorstr

!-------------------------------------------------------------------------------

USE environment      , ONLY : model_abort

!-------------------------------------------------------------------------------

USE vgrid_refatm_utils, ONLY :   vcoord

!-------------------------------------------------------------------------------

USE data_parallel    , ONLY : my_cart_id

!-------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------

    nstop,        & ! last time step of the forecast
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
    nnew,         & ! corresponds to ntstep + 1

! 2. boundary definition and update
! ---------------------------------
    nlastbound,   & ! time step of the last boundary update
    nincbound,    & ! time step increment of boundary update
    nbd1,         & ! indices for permutation of the
    nbd2,         & ! two boundary time levels

! 3. controlling the physics
! --------------------------

    lphys,        & ! forecast with physical parametrizations
    lrad,         & ! forecast with radiation
    itype_tran,   & ! type of surface-atmosphere transfer
    ltur,         & ! forecast with vertical diffusion

! 7. additional control variables
! -------------------------------

    ldump_ascii     ! for flushing (close and re-open) the ASCII files

! end of data_runcontrol 

!-------------------------------------------------------------------------------

USE data_nudge_all , ONLY :   &

! 1. Parameters and related variables
! -----------------------------------

    lwonl        ,& ! .TRUE if grid pt. (ionl ,jonl ) lies in the local domain
    lsurfa       ,& !  if surface fields to be analysed
    lfirst       ,& ! .TRUE if 'organize_nudging' is called the first time
    lvofoi       ,& ! .TRUE if observation increments are also written to VOF
    acthr        ,& ! actual forecast hour
    aiwthr       ,& ! model time [hrs] for which temporal weights are valid
!   maxmll       ,& !        : max. number of multi-level reports in the ODR
!   maxsgl       ,& !        : max. number of (surface-level and upper-air)
!                   !                         single-level reports in the ODR
!   maxgpl       ,& !        : max. number of GPS reports in the ODR
!   maxtvl       ,& !        : max. number of satellite retrievals in the ODR

! 2. Namelist variables controlling the data assimilation
! -------------------------------------------------------

    lverif       ,& ! .f. : on - off switch for verification
    mruntyp      ,& ! -1  : type of current model run used for increments in VOF
    nudgend      ,& ! 0   : end of nudging period in timesteps
    nversta      ,& ! 0   : start of verification period in timesteps
    nverend      ,& ! 0   : end of verification period in timesteps
    hversta      ,& ! 0   : start of verification period in 'model integ. hours'
    hverend      ,& ! 0   : end of verification period in 'model integr. hours'
    tconbox      ,& ! 6*dt: timestep [s] for computing analysis increments
                    !       (i.e. time box of constant analysis increments)
    ntpscor      ,& ! 1   : switch for temperature correction when the 'surface'
                    !       pressure (i.e. p. at lowest model level) is nudged
    ptpstop      ,& ! 400.: pressure [hPa] at upper boundary of the temperature
                    !       correction for 'surface' pressure nudging
    gnudg        ,& ! 6,12,6,6*10^-4: nudging coefficients for TEMP / PILOT data
    gnudgar      ,& ! 6, 0,6,0*10^-4: nudging coeffic. for AIRCRAFT data [1/s]
    gnudgms      ,& ! 6, 0,6,0*10^-4: nudging coeffic. for Mode-S aircraft [1/s]
    gnudgsu      ,& ! 6,12,0,6*10^-4: nudging coef. for surface-level data [1/s]
    gnudggp      ,& !        0*10^-4: nudging coef. for GPS-derived IWV [1/s]
    ltipol       ,& ! .t. : .t. ==> linear interpolation in time of upper-air
                    !               data which are less than 'tipolmx' hrs apart
    ltipsu       ,& ! .t. : .t. ==> linear interpolation in time of surface-lev.
                    !               data which are less than 'tipmxsu' hrs apart
    tipolmx      ,& ! 1.0 : max. time span (hrs) allowed for linear interpolat.
                    !       for upper-air data (set tipolmx = 0, if .NOT ltipol)
    tipmxsu      ,& ! 1.0 : max. time span (hrs) allowed for linear interpolat.
                    !       of surface-level data  (tipmxsu = 0, if .NOT ltipsu)
    wtukrsa      ,& ! 3.0 : for saw-tooth shaped temporal weights:
                    !       temporal radius of influence towards the past
                    !       relative to the obs. time .   for TEMP / PILOT
    wtukrse      ,& ! 1.0 : for saw-tooth shaped temporal weights:
                    !       temporal radius of influence towards the future
                    !       relative to the obs. time .   for TEMP / PILOT
    wtukara      ,& ! 1.5 : for saw-tooth-shaped temporal weights:
                    !       temporal radius of influence towards the past
                    !       relative to the obs. time .   for AIRCRAFT data
    wtukare      ,& ! 0.5 : for saw-tooth-shaped temporal weights:
                    !       temporal radius of influence towards the future
                    !       relative to the obs. time .   for AIRCRAFT data
    wtuksua      ,& ! 1.5 : for saw-tooth shaped temporal weights:
                    !       temporal radius of influence towards the past
                    !       relative to the obs. time.   for surface-level data
    wtuksue         ! 0.5 : for saw-tooth shaped temporal weights:
                    !       temporal radius of influence towards the future
                    !       relative to the obs. time.   for surface-level data

USE data_nudge_all , ONLY :   &

    msprpar      ,& ! 1   : switch specifying the surfaces along which observat.
                    !       increments of upper-air data are (primarily) spread 
    msprpsu      ,& ! 0   : switch specifying the surface along which surface-  
                    !       level data increments are primarily spreaded
    vcorls       ,& ! 2*.333,: square of the vertical correlation scale,
                    ! 2*.04    i.e. of the Gaussian vertical influence 'radius'
                    !          in log( pressure ) if (msprpar <= 1), or
                    !          in potential temperature if (msprpar == 2)
    vcutof       ,& ! 2* .75 : cut-off of the vertical correlation.       Units:
                    ! 2*1.     value of correlation at cut-off is [exp(-vcutof)]
    vcorlsu      ,& ! 2*.013,: square of the vertical correlation scale,
                    ! 2*.002   i.e. of the Gaussian vertical influence 'radius'
                    !          in log( pressure ) if (msprpsu == 1) or
                    !          in potential temperature if (msprpsu == 2)
    vcutosu      ,& ! 2* .75 : cut-off of the vertical correlation.       Units:
                    ! 2*4.     value of correlation at cut-off : [exp(-vcutosu)]
    vpblsu       ,& ! 2*99. ,: Gaussian vertical influence 'radius' of potential
                    ! 2*99.    temperature differences between obs. level and
                    !          model level, to render vertical weights depend on
                    !          the stability (in the PBL) even if (msprpsu <= 1)
    rhinfl       ,& ! 0.,70.,: constant part of the 'correlation scale of the
                    ! 0.,0.    autoregressive horiz. correlation function'
                    !          (='COSAC') [km]
    rhvfac       ,& ! 1., 0.,: multiplication factor to the vertically varying
                    ! 2* .83   part of the 'COSAC' (as def. in 'data_nudge_all')
    rhtfac       ,& ! 1.3,   : scaling factor of the total 'COSAC' for the
                    ! 1.43,    beginning and end of the nudging period for 1 obs
                    ! 1.3,     relative to the 'COSAC' at the obs. time as given
                    ! 1.3      by 'rhinfl', 'rhvfac')
    cutofr       ,& ! 4* 3.5 : cut-off in 'COSAC' units of the horizontal
                    !          correlation function
    vcsni        ,& ! 4*2500.: square of Gaussian vertical influence 'radius'
                    !          in potential temperature (if msprpar <= 1) or
                    !          log( pressure ) (if msprpar == 2) on surfaces
                    !          along which obs. increments are spread laterally
    rhiflsu      ,& ! 2 *70.,: constant part of the 'correlation scale of the
                    !   100.,  autoregressive horiz. correlation function'
                    !    70.   (='COSAC') [km]  (at the obs. time)
    rhtfsu       ,& ! 1.,    : scaling factor of the total 'COSAC' for the
                    ! 1.43.,   beginning and end of the nudging period for 1 obs
                    ! 1.,      relative to the 'COSAC' at the obs. time as given
                    ! 1.       by 'rhiflsu')
    cutofsu      ,& ! 2., 3.5: cut-off in 'COSAC' units of the horizontal
                    ! 2., 2.   correlation function
    vcsnisu      ,& ! 2*2500.: square of Gaussian vertical influence 'radius'
                    ! 2*   9.  in potential temperature (if msprpsu <= 1) or
                    !          log( pressure ) (if msprpsu == 2) on surfaces
                    !          along which obs. increments are spread laterally
    loiqv2m      ,& ! .f. : .t. ==> 2-m humidity observation increments as
                    !               differences of specific humidity
    lqfqv2m      ,& ! .f. : .t. ==> quality factor for 2-m humidity observations
                    !               dependent on T-2m differences
    dtqc         ,& ! 720.: timestep (in [s]) for the threshold quality control
    qcc          ,& !  0.,500. constant parts of the quality control thresholds
                    !  0., .7  (='QCT'). (u,v):[m/s], ps: [Pa], T: [k], RH: [ ]
    qcvf         ,& !  5., 1.: multiplication factor to the vertically varying
                    ! 10., 0.  part of the QCT (as def. in 'data_obs_qc_limits')
    qccsu        ,& ! 12.,500. constant parts of the quality control thresholds
                    ! 12., .7  (='QCT'). (u,v):[m/s], ps: [Pa], T: [k], RH: [ ]
    qcflbcp      ,& ! 1.2 : enhancement factor to the threshold used for the
                    !       check of ps against lateral boundary fields
    qcfpst       ,& ! 1.5 : maximum enhancement of the weight for surface
                    !       pressure observations due to the pressure tendency
    doromx       ,& !  100., : cut-off and gaussian radius of height differences
                    !  150.,   between model orography and station height for a
                    !  150.,   factor contributing to the quality weight factor
                    !  150.,   as part of the nudging weights
    maxmlo       ,& !        : max. number of multi-level reports in the ODR
    maxsgo       ,& !        : max. number of (surface-level and upper-air)
                    !                         single-level reports in the ODR
    maxuso       ,& !        : max. number of upper-air single-level rep. in ODR
    maxgpo       ,& !        : max. number of GPS reports within total domain
    maxmlv       ,& !  100   : max. number of obs levels in multi-level ODR
    ionl         ,& !  167   : / grid point coordinates
    jonl            !  103   : \ for standard output on nudging

USE data_nudge_all , ONLY :   &

! 5. Miscellany 
! -------------

    a_u          ,& ! zonal wind speed            on Arakawa A grid ( m/s )
    a_v          ,& ! meridional wind speed       on Arakawa A grid ( m/s )
    a_p          ,& ! pressure (full value)       on main levels    ( Pa  )
    a_z          ,& ! geometrical height          of main levels    (  m  )             
    tmaxbox      ,& ! maximum interval [s] for computing analysis increments
    cqcwbox      ,& ! (maximum) half interval [s] within which observations are
                    ! quality controlled and written on the feedobs file
    qctfpr       ,& ! for VOF output only: time factor for QC thresholds
    liwvssc      ,& ! .t. : spatial consistency check of IWV performed
    gnudgtv      ,& ! 0*10^-4: nudging coefficient for satellite retrieval [1/s]
    maxtvo          !    1   : max. number of sat retrievals within total domain

! end of data_nudge_all

!-------------------------------------------------------------------------------

! USE data_1dvar , ONLY :   &

!   ssat         ,& ! satellite specific information
!   kidsat          ! pointer index to 'ssat' for given Sat ID

! end of data_1dvar

!-------------------------------------------------------------------------------

USE data_obs_lib_cosmo , ONLY :   &

! 1. General parameters
! ---------------------

    c0         ,& ! standard real constant 0.0
    c1         ,& ! standard real constant 1.0
    c2         ,& ! standard real constant 2.0
    c05        ,& ! standard real constant 0.5
    c3600      ,& ! standard real constant 3600.0
    rmdi       ,& ! =-1.E31_wp : commonly used missing data indicator
    rmdich     ,& ! =-1.E30_wp : commonly used check value for miss data
    epsy       ,& ! = 1.E-8_wp : commonly used very small value > 0
    i0         ,& ! standard integer constant 0
    i1         ,& ! standard integer constant 1

! 2. Variables and parameters obtained by calling 'obs_cdf_interface'
! -------------------------------------------------------------------

    maxmll     ,& ! size (report dimension) of the  multi-level (m-l)  ODR
    maxsgl     ,& ! size (report dimension) of the single-level (s-l)  ODR
    maxgpl     ,& ! size (report dimension) of the (ground-based) GPS  ODR
    maxtvl     ,& ! size (report dimension) of the satellite retrieval ODR

! 3. pressure dependent scales and geometry of horizontal correlations
! --------------------------------------------------------------------

    ncolev     ,& ! number of levels in the correlation scale tables
    tabcop     ,& ! levels of the correlation scale tables
    tabcolp    ,& ! LOG( tabcop )
    rhvsond    ,& ! upper-air wind horizontal correlation scales
                  ! (pressure dependent part)
    rhtsond    ,& ! upper-air temperature horiz. correlation scales
    rhqsond    ,& ! upper-air humidity horiz. correlation scales

! 4. I/O device numbers for obs processing / nudging and file names
! -----------------------------------------------------------------

    yuprint    ,& ! file name for all the remaining information
    nupr       ,& ! unit number of file for all the remaining information

! 5. CMA observation type and code type numbers
! ---------------------------------------------

    nsynop     ,& ! SYNOP reports
    nairep     ,& ! AIREP reports (all aircraft reports)
    nsatob     ,& ! SATOB reports
    ndribu     ,& ! DRIBU reports
    ntemp      ,& ! TEMP  reports
    npilot     ,& ! PILOT reports
    ngps       ,& ! GPS   reports
    nmodes     ,& !   mode-s report
    nwp_eu     ,& !   European wind profiler report
    nra_eu     ,& !   European SODAR/RASS report
    nravad     ,& !   Radar VAD wind report
    npr_us     ,& !   US Wind Profiler/RASS report
    nsmsg1     ,& ! MSG_1  satellite retrieval
    nsmsg2     ,& ! MSG_2  satellite retrieval
    nnoa15     ,& ! NOAA15 satellite retrieval  (nnoa15 == 200 + idnoaa15)
    nnoa16     ,& ! NOAA16 satellite retrieval  (nnoa16 == 200 + idnoaa16)
    nnoa17     ,& ! NOAA17 satellite retrieval  (nnoa17 == 200 + idnoaa17)
    nnoa18        ! NOAA18 satellite retrieval  (nnoa18 == 200 + idnoaa18)

! end of data_obs_lib_cosmo

!-------------------------------------------------------------------------------

USE data_obs_record , ONLY :   &

! 1.1 ODR header format
! ---------------------

!   mxrhed       ,& ! header length of multi-level reports
!   mxshed       ,& ! header length of single-level reports
!   mxghed       ,& ! header length of GPS reports
!   mxthed       ,& ! header length of satellite retrieval reports
    nhilon       ,& ! longitude of observing station
    nhjlat       ,& ! latitude  of observing station
    nhalt        ,& ! station altitude [m]
    nhtime       ,& ! time of observat. in forecast hours
    nhsurf       ,& ! height of model grid pt. to which obs. is assigned
    nhzio        ,& ! longitude of obs. station (or lowest datum) in g.pt units
    nhzjo        ,& ! latitude  of obs. station in grid pt. units
    nhvcbu       ,& ! correction factor to vertical correlation scale for wind
                    ! at base of report
    nhvcbt       ,& ! as 'nhvcbu', but for temperature
    nhvcbq       ,& ! as 'nhvcbu', but for humidity
    nhvctu       ,& ! correction factor to vertical correlation scale for wind
                    ! at top of report
    nhvctt       ,& ! as 'nhvctu', but for temperature
    nhvctq       ,& ! as 'nhvctu', but for humidity
    nh1wta       ,& ! 1dvar retriev.: influence radius of temporal weight: past
    nh1wte       ,& ! 1dvar retriev.: influence radius of temporal w.: future
    nh1tip       ,& ! 1dvar retriev.: temporal interval for re-doing minimizat.

!   mxrhdf       ,& ! header length of multi-level reports
!   mxshdf       ,& ! header length of single-level reports
!   mxghdf       ,& ! header length of GPS reports
!   mxthdf       ,& ! header length of satellite retrieval reports
    nhio         ,& ! (local) x-coord. of grid pt. assigned to obs
    nhjo         ,& ! (local) y-coord. of grid pt. assigned to obs
    nhitot       ,& ! global x-coord. of grid pt. assigned to obs
    nhjtot       ,& ! global y-coord. of grid pt. assigned to obs
    nhobtp       ,& ! observation type
    nhcode       ,& ! code type
    nhschr       ,& ! station characteristics                      (see 1.1.4)
    nhpass       ,& ! flag for report being set to 'passive'       (see 1.1.4)
    nhqcfw       ,& ! status of QC and of writing to feedobs files, and
                    !   QC flags for surface pressure increments from TEMPs
!   nhstid       ,& ! station identity number / satellite ID (WMO Table C-5)
    nhnlev       ,& ! number of obs. levels (for multi-level reports)
    nhvqcf       ,& ! for satellite retrieval: threshold quality control flags
    nhaexi       ,& ! for conventional report: flag for exist. of wind or temp.
    nhuexi       ,& ! flag for existence of wind data      in multi-level report
    nhtexi       ,& ! flag for existence of temperat. data in multi-level report
    nhqexi       ,& ! flag for existence of humidity data  in multi-level report
    nhqcps       ,& ! temporary entry for radiosonde only (QC flag for pressure)
    nhtvip       ,& ! observed multi-lev pressure interpol to lowest model level
    nhtviz          ! vertical distance to nearest observation

USE data_obs_record , ONLY :   &

! 1.3 ODR body format
! -------------------

!   maxrsl       ,& ! max. number of levels in multi-level ODR
    maxrtv       ,& ! max. number of levels in satellite retrievals
    mxrbdy       ,& ! body length of multi-level reports
    nbtp         ,& ! pressure [Pa]
    nbtz         ,& ! height [m]
!   nbtzer       ,& ! error of observed height

    mxrbdf       ,& ! body length of multi-level reports
    nbtlid       ,& ! level identity          (bit pattern, see below: 'nb?lid')
    nbterr       ,& ! pre-processing status flags  (bit p., see below: 'nb?err')
    nbtqcf       ,& ! threshold quality control flags      (see below: 'nb?qcf')

    mxsbdy       ,& ! body length of single-level reports
    nbsu         ,& ! u wind component [m/s]
    nbsv         ,& ! v wind component [m/s]
    nbst         ,& ! temperature [K]
    nbsrh        ,& ! relative humidity [/]
    nbsp         ,& ! pressure [Pa]
    nbsz         ,& ! height [m]
    nbsuer       ,& ! error of observed wind component
    nbster       ,& ! error of observed temperature
    nbsqer       ,& ! error of observed relative humidity
!   nbszer       ,& ! error of observed height
    nbscl        ,& ! low cloud cover
    nbsct        ,& ! total cloud cover
    nbspst       ,& ! (3-hourly) pressure tendency [Pa/3h]
    nbsviu       ,& ! 10-m u-wind comp. obs extrapolated to lowest model level
    nbsviv       ,& ! 10-m v-wind comp. obs extrapolated to lowest model level
    nbsvit       ,& !  2-m temperature obs. extrapolated to lowest model level
    nbsviq       ,& !  2-m humidity observ. extrapolated to lowest model level
    nbsvip       ,& ! surface pressure obs. extrapolated to lowest model level
                    !   (if < 0 then obs is not used for 'bias' in SCC)
    nbsviz       ,& ! scaled extrapolation distance for surf. pressure obs
    mxsbdf       ,& ! body length of single-level reports
    nbserr       ,& ! pre-processing status flags  (bit p., see below: 'nb?err')
    nbsqcf       ,& ! threshold quality control flags (as in VOF)
    nbsflg       ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbscwg       ,& ! combined cloud and weather group (as in VOF)
    nbgerr       ,& ! pre-processing status flags  (bit p., see below: 'nb?err')
!   nbgtze       ,& ! error in total zenith delay [mm]

! 1.4 Bit patterns for packed information in ODR body
! ---------------------------------------------------

    nvfzbp       ,& ! bit pos. for main flag on pressure / geopot.         "
    nvfbps       ,& ! bit pattern for main flags:                          "
    nvlidp       ,& ! level id. bit pattern                             nb?lid
    nvlido       ,& ! no. bits occ. by each indicator in level id.        "
    nvnbp        ,& !         "        n   (total cloud cover)
    nvnoc           !           "             n     [Code table 2700]

USE data_obs_record , ONLY :   &

!-------------------------------------------------------------------------------
!       1.4.2  Other bit patt. for packed info in ODR (VOF) body, general words
!              ----------------------------------------------------------------
    nvru         ,& ! bit pos. for status/QC flags for horiz. wind nb?err/nb?qcf
    nvrt         ,& ! bit pos. for status/QC flags for temperature       "
    nvrq         ,& ! bit pos. for status/QC flags for humidity          "
    nvrz         ,& ! bit pos. for status/QC flags for pressure/height   "
    nvrw         ,& ! bit pos. for status/QC flags for vertical wind     "
    nvriwv       ,& ! bit pos. for status/QC flags for IWV               "
    nvrzbc       ,& ! bit pos. for temporary flag: QC ag. LBC pressure   "

! 1.5 Further quantities related to ODR
! -------------------------------------

    imdi         ,& ! missing data indicator for ODR integers (2^31-1)
    ntotml       ,& ! tot. number of stored multi-level reports
    ntotsg       ,& ! tot. number of stored single-level reports
    ntotgp       ,& ! tot. number of stored GPS reports
    ntottv       ,& ! tot. number of stored satellite retrievals
    ilstid       ,& ! character length of the station identity
    ilstidp      ,& ! char. length used for printing the station ID
                    ! Note: (ilstid >= ilstidg >= ilstidp) cf. data_nudge_gather
    ystid        ,& ! obs. station identity to be printed
    fdoro        ,& ! scaling factor to vertical distances betw. model orography
                    ! and station height for (z-obs < z-mod)

! 2. Observation data records (ODR)
! ---------------------------------

    omlbdy       ,& ! body   of multi-level ODR
    omlhed       ,& ! header of multi-level ODR
    osgbdy       ,& ! body   of single-level ODR
    osghed       ,& ! header of single-level ODR
    ogpbdy       ,& ! body   of GPS ODR
    ogphed       ,& ! header of GPS ODR
    otvbdy       ,& ! body   of satellite retrieval ODR
    otvhed       ,& ! header of satellite retrieval ODR
    momlbd       ,& ! body   of multi-level ODR
    momlhd       ,& ! header of multi-level ODR
    mosgbd       ,& ! body   of single-level ODR
    mosghd       ,& ! header of single-level ODR
    mogpbd       ,& ! body   of GPS ODR
    mogphd       ,& ! header of GPS ODR
    motvbd       ,& ! body   of satellite retrieval ODR
    motvhd       ,& ! header of satellite retrieval ODR
    yomlhd       ,& ! header of multi-level ODR
    yosghd       ,& ! header of single-level ODR
    yogphd       ,& ! header of GPS ODR
    yotvhd          ! header of satellite retrieval ODR

USE data_obs_record , ONLY :   &

! 3. Simulated Observation Record (SOR) 
! ------------------------------------- 

    mxsoml       ,& ! SOR body length for multi-level reports
    mxsosg       ,& ! SOR body length for single-level reports
    mxsogp       ,& ! SOR body length for GPS reports
    mxsotv       ,& ! SOR body length for satellite retrieval reports
    mxsops       ,& ! SOR header length for multi-level reports
    nso_u        ,& ! u wind component                               [m/s]
    nso_v        ,& ! v wind component                               [m/s]
    nso_t        ,& ! temperature                                    [K]
    nso_rh       ,& ! relative humidity                              [%] 
    nso_p        ,& ! pressure (sfc.) | geopotential (upper-air)  [Pa] | [m2/s2]
    nso_ct       ,& ! total cloud cover                              [octas]
    nso_cl       ,& ! low cloud cover                                [octas]
    nso_iq       ,& ! integrated water vapour (increment)            [mm]
    nso_ps       ,& ! pressure                                       [Pa]
    smlbdy       ,& ! body of multi-level SOR
    ssgbdy       ,& ! body of single-level SOR
    sgpbdy       ,& ! body of GPS (IWV) SOR
    stvbdy       ,& ! body of satellite retrieval SOR
    dmlhed          ! single-level part of multi-level SOR

! end of data_obs_record

!-------------------------------------------------------------------------------

USE data_nudge_gather , ONLY :   &

! 1. Parameters and general variables
! -----------------------------------

    ilstidg      ,& ! char. length used for gathering the station ID
    zalllow      ,& ! smallest allowed height within model domain
    zallhig      ,& ! largest  allowed height within model domain
    p0r          ,& ! reference pressure for potential temperature
    qst          ,& ! quotient of the mean vertical potential temperature grad.
                    ! in the stratosphere to the analog. tropospheric gradient
    mxispr       ,& !     number of equidistant vertical points
    xezispr      ,& !   / parameters used to
    dezispr      ,& !  /  define the vertically
    sezispr      ,& ! <   equidistant points
    xthispr      ,& !  \  (for non-isotropic
    dthispr      ,& !   \ horizontal correlations)
    vcutnit      ,& ! horiz. correlations are non-isotropic if the
    vcutnip      ,& ! \  vertical scales  'rdsprni' < 'zcutnit,p'
    ktp          ,& ! lowermost purely horizontal model main level
    ktth            ! top model level with spreading along isentropic surfaces

USE data_nudge_gather , ONLY :   &

! 2. The required local information on the observations and their location
! ------------------------------------------------------------------------

! local information on multi-level reports
    maxmloi_tot  ,& ! length of arrays with local info on multi-level reports
    oiml         ,& ! observation increment record OIR for multi-level reports
    zwtml        ,& ! temporal weights
    zspobml      ,& ! spreading parameter at the base / top of the profile
    fcorlml      ,& ! reduction to vertical correl. at the profile's base /top
    zvcutml      ,& ! vertical cut-off at the base / top of the profile
    zrtdgml      ,& ! convertor from height to pressure units for vertical
                    ! correlations at the base / top of the profile
    zpobml       ,& ! pressure (p) at the base / top of the profile
    zemkpml      ,& ! EXP( -R/cp *p ) at the base / top of the profile
    zriflml      ,& ! upper estimate to horizontal radius of influence
    zstpml       ,& ! spreading parameter at the tropopause & obs. location
    zsprml       ,& ! spreading parameter at model levels & obs. location
    zlopml       ,& ! log( pressure ) at model levels & obs. location
    zrtgkml      ,& ! convertor from height to pressure units for vertical
                    ! correlations at model levels & obs. location
    znisml       ,& ! parameter used for non-isotropic horizontal correlations
                    ! at model levels & obs. location
    znismq       ,& ! parameter used for non-isotropic horizontal correlations
                    ! at vertically equidistant points & obs. location
    ioml         ,& ! longitudinal index of station location on local domain
    joml         ,& ! latitudinal  index of station location on local domain
    ioml_tot     ,& ! longitudinal index of station location on global domain
    joml_tot     ,& ! latitudinal  index of station location on global domain
    kviflml      ,& ! vertical range of possibly influenced model levels
    kkoiml       ,& ! level index in the OIR (observation increment record)
    ksprml       ,& ! lowest level with spreading along model levels
    kobtyml      ,& ! observation type
    kcdtyml      ,& ! CMA observation code type
    mszlev       ,& ! number of vertical levels with obs. incr. in the OIR
    mbotlv       ,& ! number of model levels above the surface without
                    ! observation increment
    mtoplv       ,& ! number of model levels at the top of the model without
                    ! observation increment
    ltiml        ,& ! .TRUE if temporal linear interpol. for multi-level rep.
    ystidml         ! station identity of multi-level station

USE data_nudge_gather , ONLY :   &

! local information on upper-air single-level reports
    xoiua        ,& ! observation increments
    zwtua        ,& ! temporal weights
    zsprob       ,& ! spreading parameter at the obs. point
    fcorlua      ,& ! reduction to vertical correl. below / above the obs.
    zvcutua      ,& ! vertical cut-off below / above the obs.
    zrtdgua      ,& ! convertor height to pressure units for vertic. correlat.
    zpobua       ,& ! pressure at the obs. point
    zzobua       ,& ! height at the obs. point
    zriflua      ,& ! upper estimate to horizontal radius of influence
    zvidua       ,& ! vertical interpolation distance of adjacent model levels
    zqualua      ,& ! quality weight factor
    zsprtp       ,& ! spreading parameter at the tropopause & obs. location
    zsprua       ,& ! spreading parameter at model levels & obs. location
    zlopua       ,& ! log( pressure ) at model levels & obs. location
    znisua       ,& ! parameter used for non-isotropic horizontal correlations
                    ! at model levels & obs. location
    znisuq       ,& ! parameter used for non-isotropic horizontal correlations
                    ! at vertically equidistant points & obs. location
    ioua         ,& ! longitudinal index of station location on local domain
    joua         ,& ! latitudinal  index of station location on local domain
    ioua_tot     ,& ! longitudinal index of station location on global domain
    joua_tot     ,& ! latitudinal  index of station location on global domain
    kviflua      ,& ! vertical range of possibly influenced model levels
    kobtyua      ,& ! CMA observation type
    kcdtyua      ,& ! CMA observation code type
    ltiua        ,& ! .TRUE if temporal linear interpol. for upper-air sing.l.
    ystidua         ! station identity of upper-air single-level station

USE data_nudge_gather , ONLY :   &

! local information on surface-level reports
    xoisu        ,& ! observation increments
    zwtsu        ,& ! temporal weights
    zsposu       ,& ! spreading parameter at the obs. point
    zvcutsu      ,& ! vertical cut-off above the obs.
    zrtdgsu      ,& ! convertor height to pressure units for vertic. correlat.
    zpobsu       ,& ! pressure at the obs. point
    zzobsu       ,& ! height at the obs. point
    zriflsu      ,& ! upper estimate to horizontal radius of influence
    zqualsu      ,& ! quality weight factor
    zsprsu       ,& ! spreading parameter at model levels & obs. location
    zpblsu       ,& ! potential temperature difference rel. to the obs. level
    znissu       ,& ! parameter used for non-isotropic horizontal correlations
                    ! at model levels & obs. location
    znissq       ,& ! parameter used for non-isotropic horizontal correlations
                    ! at vertically equidistant points & obs. location
    iosu         ,& ! longitudinal index of station location on local domain
    josu         ,& ! latitudinal  index of station location on local domain
    iosu_tot     ,& ! longitudinal index of station location on global domain
    josu_tot     ,& ! latitudinal  index of station location on global domain
    kviflsu      ,& ! vertical range of possibly influenced model levels
    kobtysu      ,& ! CMA observation type
    kcdtysu      ,& ! CMA observation code type
    ltisu        ,& ! .TRUE if temporal linear interpol. for surface-level rep
    ystidsu         ! station identity of surface-level station

USE data_nudge_gather , ONLY :   &

! local information on 'surface' (i.e. on lowest model level) pressure data
    zoips        ,& ! observation increments
    omykps       ,& ! local part of nudging weight, i.e. exclud. horiz. weight
    r1ifps       ,& ! horizontal correlation scale
    zmassps      ,& ! total mass affected by the 'temperature correction'
    ztdpps       ,& ! 'temperature / pressure' at obs. point
    qcimps       ,& ! model pressure interpolated to obs. point
    qctmps       ,& ! time of all observations
    qcqfps       ,& ! weight factor related to obs. quality (representiveness)
    zoips_b      ,& ! increments: observations minus lateral boundary fields
    iops         ,& ! longitudinal index of station location on local domain
    jops         ,& ! latitudinal  index of station location on local domain
    iops_tot     ,& ! longitudinal index of station location on global domain
    jops_tot     ,& ! latitudinal  index of station location on global domain
    iqclps       ,& ! indicator for current obs. to undergo spatial check now
    iqcfps       ,& ! indicator for obs. to have passed latest threshold QC
    iqcnps       ,& ! index of local administrator of (local) ODR
    ltips        ,& ! .TRUE if temporal linear interpol. for surface pressure
    lmlps        ,& ! .TRUE if datum is derived from multi-level report
    ystidps         ! station identity of 'surface' pressure station

USE data_nudge_gather , ONLY :   &

! local information on integrated water vapour increments (for humidity QC)
    zoiciv       ,& ! observation increments
    zqcfiv       ,& ! local part of nudging weight, i.e. exclud. horiz. weight
    zmodiv       ,& ! model integrated water vapour IWV
    zsativ       ,& ! IWV of saturated model temperature profile
    zactiv       ,& ! observation time
    ioiv         ,& ! longitudinal index of station location on local domain
    joiv         ,& ! latitudinal  index of station location on local domain
    ioiv_tot     ,& ! longitudinal index of station location on global domain
    joiv_tot     ,& ! latitudinal  index of station location on global domain
    iqcliv       ,& ! indicator for current obs. to undergo spatial check now
    iqcfiv       ,& ! indicator for obs. to have passed latest threshold QC
    kiobiv       ,& ! index of local administrator of (local) ODR
    kioiiv       ,& ! index of local information array for multi-level data
    ktypiv       ,& ! observation type
    ystidiv         ! station identity of integrated water vapour station

USE data_nudge_gather , ONLY :   &

! 3. Observation increment record for multi-level reports 'oiml'
! --------------------------------------------------------------

    maxnoi       ,& ! length of observation increment record

! 4. (Local) information gathered by 1 (or 2) nodes for printing for control
! --------------------------------------------------------------------------

    oyqc         ,& ! on data rejected by the threshold quality control
    myqc         ,& ! on data rejected by the threshold quality control
    yyqc         ,& ! on data rejected by the threshold quality control
    oysu         ,& ! on interpol. of surface-level data to lowest model level
    yysu         ,& ! on interpol. of surface-level data to lowest model level

! 5. Variables defining the size of the arrays containing the local information
! -----------------------------------------------------------------------------

    maxoil       ,& ! max.  number of vertical levels in the OIR 'oiml'
    maxpso       ,& ! max.  number of active 'surface' pressure data
    maxivq       ,& ! max.  number of IWV reports used for spatial check
    maxqcp       ,& ! max.  number of rejected data to be printed per timestep
    maxysu       ,& ! max.  number of extrapolated printed surface-level reports
    nuatot       ,& ! total number of active upper-air single-level stations
    nsutot       ,& ! total number of active surface-level stations
    npstot       ,& ! total number of active surface-pressure stations
    nivtot       ,& ! total number of IWV reports used for spatial check
    ntotqc       ,& ! total number of rejected data to be printed per timestep
    ntotys       ,& ! total number of printed interpol. surface-level reports
    nuaex        ,& ! number of local upper-air obs i. not used due to ODR size 
    ktopsu       ,& ! uppermost model level possibly influenced by surface obs
    lnissu       ,& ! non-isotrophic correlations for surface-level data
    lnisua       ,& ! non-isotrophic correlat. for upper-air single-lev. data

! 6. Geometrics and variables used for spatial consistency check of pressure obs
! ------------------------------------------------------------------------------

    pyyd         ,& ! latitudinal (merid.) distance on tangent cone projection
    pxxd2        ,& ! square of zonal distance on tangent cone projection
    pxsalpa      ,& ! factor used for distances on tangent cone projection
    isrtpqc      ,& ! (sorted) list of stations with 'surface' pressure data
    isrtvqc         ! (sorted) list of stations with IWV

! end of data_nudge_gather

!-------------------------------------------------------------------------------

USE data_obs_qc_limits , ONLY :   &

! 2. Table for pressure dependent quality control thresholds
! ----------------------------------------------------------

    nqclev       ,& ! number of levels in the correlation scale tables
    tabqcp       ,& ! levels of the correlation scale tables
    tabqclp      ,& ! LOG( tabqcp(15) )
    qcvairp      ,& ! (root of) AIREP wind error variance
    qctairp      ,& ! (root of) AIREP temperature error variance
    qczcorl      ,& ! radiosonde height error correlation matrix
    nzcorl       ,& ! radiosonde height error correlation matrix [*1000]
    qctf         ,& ! temporal factor for thresholds of upper-air obs.
    qctfsu       ,& ! temporal factor for thresholds of surface-level obs.
    qcvairf      ,& ! additional factor used for airep wind thresholds
    qctairf      ,& ! additional factor used for airep temperature thresholds
    qcvdrbf      ,& ! additional factor used for dribu wind thresholds
    qczdrbf      ,& ! additional factor used for dribu pressure thresholds
    qcvairl      ,& ! background wind speed rejection limit for zero AIREP winds
    qcqfge1      ,& ! first guess humidity error at latitude > qcqlatl
    qcqfgef      ,& ! factors for first guess humidity errors
    qcftend      ,& ! fraction of observed pressure tendency added to threshold
    qcfctnd      ,& ! reduction to namelist threshold, if tendency is observed
    qcfgood      ,& ! reduction of surface pressure threshold for 'good' obs.
    qcfpbad      ,& ! increase of pressure threshold for 'probably bad' obs.
    qcfspat      ,& ! threshold reduction for the spatial consistency check
    qcfbias      ,& ! fraction of bias added to spatial consistency threshold
    dtchkps         ! spatial check: radius of infl. for linear temporal weights

! end of data_obs_qc_limits

!-------------------------------------------------------------------------------

USE data_nudge_local , ONLY :   &

! 3. General parameters and other variables
! -----------------------------------------

    ptropop      ,& ! pressure at tropopause (standard atmosphere)
    thdzt        ,& ! mean vert. pot. temp. gradient in troposphere
    fsvcut       ,& ! fraction with which the adjusted scale S_adj, rather than
                    ! the scale S_nl as specified in the namelist, is used for
                    ! determin. the vertical cut-off in log(p) or theta(e)-units
    wtml         ,& ! temporal nudging weights for multi-level reports
    wtsg         ,& ! temporal nudging weights for single-level reports
    wtgp         ,& ! temporal nudging weights for gps reports
    wttv         ,& ! temporal nudging weights for satellite retrievals
    tmladm       ,& ! administrator of multi-level ODR
    tsgadm       ,& ! administrator of single-level ODR
    tgpadm       ,& ! administrator of gps ODR
    ttvadm       ,& ! administrator of satellite retrievals
    zsdni        ,& ! equidistant pts. used for non-isotropic horiz. correlation
    imladm       ,& ! administrator of multi-level ODR
    isgadm       ,& ! administrator of single-level ODR
    igpadm       ,& ! administrator of gps ODR
    itvadm       ,& ! administrator of satellite retrievals
    nmlsta       ,& ! total number of sorted multi-level stations
    ngpsta       ,& ! total number of sorted gps stations
    ntvsta       ,& ! total number of sorted satellite retrieval 'stations'
    kml300       ,& ! index of full model level corresponding to about 300 hPa
    lqcall          ! .TRUE if threshold quality control for all observations
                    !          at current timestep

! end of data_nudge_local

!-------------------------------------------------------------------------------

USE src_obs_operator_conv, ONLY :   &

    sing_obs_operator   ,& ! forward observation operator for single-level obs
    surf_obs_operator   ,& ! forward observation operator for surface-level obs
    cloud_obs_operator  ,& ! forward observation operator for cloud obs
    ps_obs_operator_sing,& ! obs operator for surface pressure from sfc. station
    ps_obs_operator_mult   ! obs operator for surface pressure from radiosondes

! end of src_obs_operator_conv

!-------------------------------------------------------------------------------

USE src_obs_qc_conv      , ONLY :   &

    sing_quality_cntl   ,& ! quality control of individual single-level observ.
    ps_quality_cntl        ! quality control of individual surface pressure obs

! end of src_obs_qc_conv

!-------------------------------------------------------------------------------
! Local declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE

! 1. Variables
! ------------

  INTEGER (KIND=iintegers)                , PRIVATE :: &
    io     ,jo     ,& ! local  indices of location of observation
    io_tot ,jo_tot ,& ! global indices of location of observation
    ista           ,& ! index of observing station
    itim           ,& ! (time) index over obs. at one observing station
    nuasta         ,& ! total number of sorted upper-air single-level stations
    nsusta         ,& ! total number of sorted surface-level stations
    kml850         ,& ! index of full model level corresponding to about 850 hPa
    kml700            ! index of full model level corresponding to about 700 hPa

  REAL (KIND=wp),     POINTER :: &
    qv  (:,:,:)  => NULL(), &  ! QV at nnew
    qc  (:,:,:)  => NULL()     ! QC at nnew

!===============================================================================

!-------------------------------------------------------------------------------
! Public and Private Subroutines
!-------------------------------------------------------------------------------

CONTAINS

! INCLUDE "local_info_aux.incf"
! INCLUDE "ps_local_info.incf"
! INCLUDE "local_sort_reports.incf"
! INCLUDE "surf_local_info.incf"
! INCLUDE "surf_obs.incf"
! INCLUDE "upair_local_info.incf"

!-------------------------------------------------------------------------------
!+ Module procedure for auxiliaries for computing local information on obs
!-------------------------------------------------------------------------------

SUBROUTINE local_info_aux ( nactio )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure performs various actions and computes auxiliary
!   quantities used for the computation of the local information (on the
!   observation increments and their location). Namely:
!   - computation of time boxes and level indices
!   - allocation of arrays containing the local information
!   - checking if sorting or quality control need to be done at current timestep
!   - equidistant points for non-isotropic 1-dimensional horizontal correlations
!   - getting record with interpolated surface-level reports for printout
!
! Method:
!   Straightforward actions and computations.
!
! Written by        :  Christoph Schraff, DWD  (original version: 14.10.97)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

! Subroutine arguments:
! --------------------

  INTEGER (KIND=iintegers), INTENT (IN)         ::       &
    nactio              ! determines, which quantities are to be computed (when)

! Local parameters:
! ----------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    ktopth  = 0           ! lowest model level where lateral spreading is
                          ! always along model surfaces

  LOGICAL                  , PARAMETER  :: &
    ltopth  = .TRUE.      ! .T. ==> lateral spreading is horizontal at horizon-
                          ! tal model levels even if (msprpar == 2)

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    mv     , kk      ,& ! loop indices
    itic             ,& ! time index of report
    icsu   , icps    ,& ! indices of obs. station in data record
    icua             ,& ! index   of obs. station in data record
    nsgob  , nmlob   ,& ! indices of reports in the ODR (obs data record)
    nvys   , ndys    ,& ! loop index of array for printing surface-level reports
    npspr  , mxprspr ,& ! loop index and limit for printout
    nspe   , nspa    ,& ! interval of indices for printout
    nzerr               ! status of memory (de-)allocation

  REAL    (KIND=wp   )     ::  &
    zcosign             ! sign of the gradient of the vertical coordinate with
                        ! increasing height

  LOGICAL                  ::  &
    lrepsu           ,& ! record element has surface-level data  \ with non-zero
    lrepps              ! record element has surface presssure   / weight

  REAL    (KIND = wp)      , SAVE ::  &
    tminbox             ! min. timestep for computing analysis increments

! Local (automatic) arrays: None
! -------------------------
!
!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine local_info_aux
!-------------------------------------------------------------------------------
  
!-------------------------------------------------------------------------------
!  Section 1: Actions to be done only once at the first timestep
!-------------------------------------------------------------------------------

! Note that these quantities must be saved 'globally', i.e. they must be known
! even after quitting modul 'nudging'

  IF ((nactio == 1) .AND. (lfirst)) THEN

! time boxes
! ----------
    tmaxbox = REAL(INT( tconbox/dt   + c1-epsy ),wp) * dt
    tminbox = REAL(INT( tconbox/dt       +epsy ),wp) * dt
    IF (rmod( tconbox,dt,epsy ) < c2*epsy) THEN
      cqcwbox = REAL(INT( dtqc  /tconbox + c1-epsy ),wp) * tconbox * c05
    ELSE
      cqcwbox = REAL(INT( (dtqc +tmaxbox)/dt - epsy),wp) * dt * c05
    ENDIF
    nversta = INT( hversta /dtdeh - cqcwbox /dt +c1 -epsy )
    nverend = INT( hverend /dtdeh + cqcwbox /dt     +epsy )

! index of the lowest horizontal main model level (ktp == 0 is ok.)
! -----------------------------------------------
    zcosign = SIGN( c1 , vcoord%vert_coord(1) - vcoord%vert_coord(2) )
    ktp = 0
    DO kk = 1 , ke+1
      IF (zcosign*vcoord%vert_coord(kk) >= zcosign*vcoord%vcflat-epsy) ktp = kk
    ENDDO

! uppermost main model level with spreading along isentropic surfaces
! (if spreading along isentropic surfaces is selected)
! -------------------------------------------------------------------
    IF (.NOT. ltopth) ktth = ktopth
    IF (      ltopth) ktth = ktp + 1

! LOG( pressure ) of the levels used for tables for
! horizontal correlations and quality control thresholds
! ------------------------------------------------------
    DO kk = 1 , 11
      tabcolp (kk) = LOG( tabcop(kk) )
    ENDDO
    DO kk = 1 , nqclev
      tabqclp (kk) = LOG( tabqcp(kk) )
    ENDDO
    DO kk = 1 , nqclev
      DO mv = 1 , nqclev
        qczcorl(kk,mv) = nzcorl(kk,mv) * 0.001_wp
      ENDDO
    ENDDO
    IF (qcvf(4) > epsy) THEN
      qctf   (1) = qctf  (1) * 0.50_wp
      qctf   (3) = qctf  (3) * 0.50_wp
      qctf   (4) = qctf  (4) * 0.17_wp
      qctfsu (1) = qctfsu(1) * 0.50_wp
      qctfsu (3) = qctfsu(3) * 0.50_wp
      qctfsu (4) = qctfsu(4) * 0.17_wp
    ENDIF
    DO kk = 1 , 4
      qctfpr (kk  ) = qctf  (kk)
      qctfpr (kk+4) = qctfsu(kk)
    ENDDO

! indices of (full) model levels corresponding to about 850 and 700 hPa
! ---------------------------------------------------------------------
    kml850 = ke
    kml700 = 1
    kml300 = 1
    DO kk = 1 , ke
! for reference pressure   p0sl = 10E5_wp
      IF (      (vcoord%sigm_coord(kk)  <= 0.85_wp)                 &
          .AND. (vcoord%sigm_coord(kk+1) > 0.85_wp)) kml850 = kk
      IF (      (vcoord%sigm_coord(kk)  <= 0.70_wp)                 &
          .AND. (vcoord%sigm_coord(kk+1) > 0.70_wp)) kml700 = kk
      IF (      (vcoord%sigm_coord(kk)  +                               &
                 vcoord%sigm_coord(MAX(kk-i1,i1)) <= 0.601_wp)      &
          .AND. (vcoord%sigm_coord(kk+1) +                              &
                 vcoord%sigm_coord(kk  )           > 0.601_wp)) kml300 = kk
    ENDDO
    kml850 = MAX( kml850 , INT( kml700 + 1 ,iintegers) )

! decide if increments are written to VOF
! ---------------------------------------
    lvofoi  =  (mruntyp >= 0)

! safety allocations (if these array are not allocated here, there is a strange
! ------------------  (lethal) problem in 'local_sort_reports' even though
!                     ntottv is equal to zero and the arrays should not be used)
    ALLOCATE ( otvbdy( 1,1,1 ), STAT=nzerr)
    ALLOCATE ( otvhed( 1,  1 ), STAT=nzerr)
    ALLOCATE ( motvbd( 1,1,1 ), STAT=nzerr)
    ALLOCATE ( motvhd( 1,  1 ), STAT=nzerr)
    ALLOCATE ( yotvhd( 1     ), STAT=nzerr)
    otvbdy (:,:,:) = rmdi
    otvhed (:,:)   = rmdi
    motvbd (:,:,:) = imdi
    motvhd (:,:)   = imdi
    yotvhd (:)     = '        '

  ENDIF

!-------------------------------------------------------------------------------
!  Section 2: Actions to be taken at every timestep before getting
!             the local information on the observations and their location
!-------------------------------------------------------------------------------

  IF (nactio == 2) THEN

! allocation of adiministrator arrays, set in subr. 'local_sort_reports'
! ----------------------------------------------------------------------

    IF (ALLOCATED( imladm )) THEN
      PRINT '("CAUTION in src_sing_local: imladm is already allocated "        &
            &,"at time ",I6)', ntstep
      DEALLOCATE ( imladm , isgadm , igpadm , itvadm , STAT=nzerr )
      DEALLOCATE ( tmladm , tsgadm , tgpadm , ttvadm , STAT=nzerr )
      DEALLOCATE ( wtml   , wtsg   , wtgp   , wttv   , STAT=nzerr )
    ENDIF
    ALLOCATE ( imladm  (MAX( maxmll,i1 ),7) , STAT = nzerr )
    ALLOCATE ( isgadm  (MAX( maxsgl,i1 ),8) , STAT = nzerr )
    ALLOCATE ( igpadm  (MAX( maxgpl,i1 ),6) , STAT = nzerr )
    ALLOCATE ( itvadm  (MAX( maxtvl,i1 ),6) , STAT = nzerr )
    ALLOCATE ( tmladm  (MAX( maxmll,i1 ),3) , STAT = nzerr )
    ALLOCATE ( tsgadm  (MAX( maxsgl,i1 ),3) , STAT = nzerr )
    ALLOCATE ( tgpadm  (MAX( maxgpl,i1 ),3) , STAT = nzerr )
    ALLOCATE ( ttvadm  (MAX( maxtvl,i1 ),3) , STAT = nzerr )

    ALLOCATE ( wtml    (MAX( maxmll,i1 ),4) , STAT = nzerr )
    ALLOCATE ( wtsg    (MAX( maxsgl,i1 ),4) , STAT = nzerr )
    ALLOCATE ( wtgp    (MAX( maxgpl,i1 ),4) , STAT = nzerr )
    ALLOCATE ( wttv    (MAX( maxtvl,i1 ),4) , STAT = nzerr )

! allocation of local info arrays,
! incl. arrays for spatial consistency checks and for QC printout
! ---------------------------------------------------------------

    maxoil      =  maxmlv  +  ke
    maxmloi_tot =  maxmlo  +  maxgpo/2  +  maxtvo
    maxpso      =  maxmlo  +  maxsgo
    maxivq      =  maxmlo  +  maxgpo
    maxqcp      = (maxsgo  +  maxmlo *maxoil  +  maxgpo/2) /20  +  1

    IF (msprpar <= 1) lnisua = (MIN( vcsni(1), vcsni(3), vcsni(4) ) <= vcutnit)
    IF (msprpar == 2) lnisua = (MIN( vcsni(1), vcsni(3), vcsni(4) ) <= vcutnip)
    IF (msprpsu <= 1) THEN
      lnissu  = (MIN( vcsnisu(1) , vcsnisu(3) , vcsnisu(4) ) <= vcutnit)
    ELSE
      lnissu  = (MIN( vcsnisu(1) , vcsnisu(3) , vcsnisu(4) ) <= vcutnip)
    ENDIF

    IF (ALLOCATED( zwtml )) THEN
      PRINT '("WARNING in src_sing_local: zwtml is already allocated "         &
            &,"at time ",I6)', ntstep
      DEALLOCATE ( zwtml  , STAT=nzerr )
    ENDIF

!   ALLOCATE ( oiml    (maxmlo,maxoil,maxnoi) , STAT = nzerr )
    ALLOCATE ( zwtml   (maxmloi_tot,4,3) , STAT = nzerr )
    ALLOCATE ( zspobml (maxmloi_tot,4,3) , STAT = nzerr )
    ALLOCATE ( fcorlml (maxmloi_tot,4,3) , STAT = nzerr )
    ALLOCATE ( zvcutml (maxmloi_tot,4,3) , STAT = nzerr )
    ALLOCATE ( zrtdgml (maxmloi_tot,4,3) , STAT = nzerr )
    ALLOCATE ( zpobml  (maxmloi_tot,4,3) , STAT = nzerr )
    ALLOCATE ( zemkpml (maxmloi_tot,4,3) , STAT = nzerr )
    ALLOCATE ( kviflml (maxmloi_tot,4,3) , STAT = nzerr )
    ALLOCATE ( zriflml (maxmloi_tot,2,3) , STAT = nzerr )
    ALLOCATE ( mbotlv  (maxmloi_tot,2,3) , STAT = nzerr )
    ALLOCATE ( mtoplv  (maxmloi_tot,2,3) , STAT = nzerr )
    ALLOCATE ( ltiml   (maxmloi_tot  ,3) , STAT = nzerr )
    ALLOCATE ( mszlev  (maxmloi_tot,2  ) , STAT = nzerr )
    ALLOCATE ( kkoiml  (maxmloi_tot,2  ) , STAT = nzerr )
    ALLOCATE ( ksprml  (maxmloi_tot    ) , STAT = nzerr )
    ALLOCATE ( kobtyml (maxmloi_tot    ) , STAT = nzerr )
    ALLOCATE ( kcdtyml (maxmloi_tot    ) , STAT = nzerr )
    ALLOCATE ( zstpml  (maxmloi_tot    ) , STAT = nzerr )
    ALLOCATE ( ystidml (maxmloi_tot    ) , STAT = nzerr )
    ALLOCATE ( ioml    (maxmloi_tot    ) , STAT = nzerr )
    ALLOCATE ( joml    (maxmloi_tot    ) , STAT = nzerr )
    ALLOCATE ( ioml_tot(maxmloi_tot    ) , STAT = nzerr )
    ALLOCATE ( joml_tot(maxmloi_tot    ) , STAT = nzerr )
    ALLOCATE ( zsprml  (maxmloi_tot,ke ) , STAT = nzerr )
    ALLOCATE ( zlopml  (maxmloi_tot,ke ) , STAT = nzerr )
    ALLOCATE ( zrtgkml (maxmloi_tot,ke ) , STAT = nzerr )
    IF (lnisua) THEN
      IF (ALLOCATED( znisml )) THEN
        PRINT '("WARNING in src_sing_local: znisml is already allocated "      &
              &,"at time ",I6)', ntstep
        DEALLOCATE ( znisml  , STAT=nzerr )
      ENDIF
      ALLOCATE ( znisml  (maxmloi_tot,ke ) , STAT = nzerr )
    ENDIF
    IF ((lnisua) .AND. (msprpar >= 1))                                         &
      ALLOCATE ( znismq  (maxmloi_tot,mxispr) , STAT = nzerr )

    ALLOCATE ( xoiua   (maxuso,2,4) , STAT = nzerr )
    ALLOCATE ( zwtua   (maxuso,2,3) , STAT = nzerr )
    ALLOCATE ( fcorlua (maxuso,2,3) , STAT = nzerr )
    ALLOCATE ( zvcutua (maxuso,2,3) , STAT = nzerr )
    ALLOCATE ( zriflua (maxuso,2,3) , STAT = nzerr )
    ALLOCATE ( zqualua (maxuso,2,3) , STAT = nzerr )
    ALLOCATE ( kviflua (maxuso,2,3) , STAT = nzerr )
    ALLOCATE ( ltiua   (maxuso  ,3) , STAT = nzerr )
    ALLOCATE ( zsprob  (maxuso    ) , STAT = nzerr )
    ALLOCATE ( zsprtp  (maxuso    ) , STAT = nzerr )
    ALLOCATE ( zrtdgua (maxuso    ) , STAT = nzerr )
    ALLOCATE ( zpobua  (maxuso    ) , STAT = nzerr )
    ALLOCATE ( zzobua  (maxuso    ) , STAT = nzerr )
    ALLOCATE ( zvidua  (maxuso    ) , STAT = nzerr )
    ALLOCATE ( kobtyua (maxuso    ) , STAT = nzerr )
    ALLOCATE ( kcdtyua (maxuso    ) , STAT = nzerr )
    ALLOCATE ( ystidua (maxuso    ) , STAT = nzerr )
    ALLOCATE ( ioua    (maxuso    ) , STAT = nzerr )
    ALLOCATE ( joua    (maxuso    ) , STAT = nzerr )
    ALLOCATE ( ioua_tot(maxuso    ) , STAT = nzerr )
    ALLOCATE ( joua_tot(maxuso    ) , STAT = nzerr )
    ALLOCATE ( zsprua  (maxuso,ke ) , STAT = nzerr )
    ALLOCATE ( zlopua  (maxuso,ke ) , STAT = nzerr )
    IF (lnisua) ALLOCATE ( znisua  (maxuso,ke ) , STAT = nzerr )
    IF ((lnisua) .AND. (msprpar >= 1))                                         &
      ALLOCATE ( znisuq  (maxuso,mxispr) , STAT = nzerr )

    ALLOCATE ( xoisu   (maxsgo,2,4) , STAT = nzerr )
    ALLOCATE ( zwtsu   (maxsgo,2,3) , STAT = nzerr )
    ALLOCATE ( zvcutsu (maxsgo  ,3) , STAT = nzerr )
    ALLOCATE ( zriflsu (maxsgo  ,3) , STAT = nzerr )
    ALLOCATE ( zqualsu (maxsgo,2,3) , STAT = nzerr )
    ALLOCATE ( kviflsu (maxsgo  ,3) , STAT = nzerr )
    ALLOCATE ( ltisu   (maxsgo  ,3) , STAT = nzerr )
    ALLOCATE ( zsposu  (maxsgo    ) , STAT = nzerr )
    ALLOCATE ( zrtdgsu (maxsgo    ) , STAT = nzerr )
    ALLOCATE ( zpobsu  (maxsgo    ) , STAT = nzerr )
    ALLOCATE ( zzobsu  (maxsgo    ) , STAT = nzerr )
    ALLOCATE ( kobtysu (maxsgo    ) , STAT = nzerr )
    ALLOCATE ( kcdtysu (maxsgo    ) , STAT = nzerr )
    ALLOCATE ( ystidsu (maxsgo    ) , STAT = nzerr )
    ALLOCATE ( iosu    (maxsgo    ) , STAT = nzerr )
    ALLOCATE ( josu    (maxsgo    ) , STAT = nzerr )
    ALLOCATE ( iosu_tot(maxsgo    ) , STAT = nzerr )
    ALLOCATE ( josu_tot(maxsgo    ) , STAT = nzerr )
    ALLOCATE ( zsprsu  (maxsgo,ke ) , STAT = nzerr )
    ALLOCATE ( zpblsu  (maxsgo,ke ) , STAT = nzerr )
    IF (lnissu) ALLOCATE ( znissu  (maxsgo,ke ) , STAT = nzerr )
    IF ((lnissu) .AND. (msprpsu >= 1))                                         &
      ALLOCATE ( znissq  (maxsgo,mxispr) , STAT = nzerr )

    ALLOCATE ( zoips   (maxpso,2  ) , STAT = nzerr )
    ALLOCATE ( omykps  (maxpso,2  ) , STAT = nzerr )
    ALLOCATE ( r1ifps  (maxpso    ) , STAT = nzerr )
    ALLOCATE ( zmassps (maxpso    ) , STAT = nzerr )
    ALLOCATE ( ztdpps  (maxpso    ) , STAT = nzerr )
    ALLOCATE ( qcimps  (maxpso,2  ) , STAT = nzerr )
    ALLOCATE ( qctmps  (maxpso,2  ) , STAT = nzerr )
    ALLOCATE ( qcqfps  (maxpso,2  ) , STAT = nzerr )
    ALLOCATE ( zoips_b (maxpso,2  ) , STAT = nzerr )
    ALLOCATE ( ltips   (maxpso    ) , STAT = nzerr )
    ALLOCATE ( lmlps   (maxpso    ) , STAT = nzerr )
    ALLOCATE ( ystidps (maxpso    ) , STAT = nzerr )
    ALLOCATE ( iops    (maxpso    ) , STAT = nzerr )
    ALLOCATE ( jops    (maxpso    ) , STAT = nzerr )
    ALLOCATE ( iops_tot(maxpso    ) , STAT = nzerr )
    ALLOCATE ( jops_tot(maxpso    ) , STAT = nzerr )
    ALLOCATE ( iqclps  (maxpso    ) , STAT = nzerr )
    ALLOCATE ( iqcfps  (maxpso    ) , STAT = nzerr )
    ALLOCATE ( iqcnps  (maxpso    ) , STAT = nzerr )

    IF (ALLOCATED( zoiciv )) THEN
      PRINT '("WARNING in src_sing_local: zoiciv is already allocated "        &
            &,"at time ",I6)', ntstep
      DEALLOCATE ( zoiciv  , STAT=nzerr )
    ENDIF

    ALLOCATE ( zoiciv  (maxivq,2  ) , STAT = nzerr )
    ALLOCATE ( zqcfiv  (maxivq,2  ) , STAT = nzerr )
    ALLOCATE ( zmodiv  (maxivq    ) , STAT = nzerr )
    ALLOCATE ( zsativ  (maxivq    ) , STAT = nzerr )
    ALLOCATE ( zactiv  (maxivq,2  ) , STAT = nzerr )
    ALLOCATE ( ystidiv (maxivq    ) , STAT = nzerr )
    ALLOCATE ( ioiv    (maxivq    ) , STAT = nzerr )
    ALLOCATE ( joiv    (maxivq    ) , STAT = nzerr )
    ALLOCATE ( ioiv_tot(maxivq    ) , STAT = nzerr )
    ALLOCATE ( joiv_tot(maxivq    ) , STAT = nzerr )
    ALLOCATE ( iqcliv  (maxivq    ) , STAT = nzerr )
    ALLOCATE ( iqcfiv  (maxivq    ) , STAT = nzerr )
    ALLOCATE ( kiobiv  (maxivq    ) , STAT = nzerr )
    ALLOCATE ( kioiiv  (maxivq    ) , STAT = nzerr )
    ALLOCATE ( ktypiv  (maxivq    ) , STAT = nzerr )

    IF (ALLOCATED( oyqc )) THEN
      PRINT '("CAUTION in src_sing_local: oyqc is already allocated "          &
            &,"at time ",I6)', ntstep
      DEALLOCATE ( oyqc , myqc , yyqc , STAT=nzerr )
    ENDIF
    ALLOCATE ( oyqc    (maxqcp, 11) , STAT = nzerr )
    ALLOCATE ( myqc    (maxqcp,  2) , STAT = nzerr )
    ALLOCATE ( yyqc    (maxqcp    ) , STAT = nzerr )

    IF (ALLOCATED( zsdni ))  DEALLOCATE ( zsdni , STAT=nzerr )
    ALLOCATE ( zsdni   (mxispr,2) , STAT = nzerr )

! allocation of obs increment arrays used for writing to VOF / feedobs file
! -------------------------------------------------------------------------

    IF (ALLOCATED( ssgbdy )) THEN
      IF (lverif) THEN
        PRINT '("CAUTION in src_sing_local: ssgbdy is already allocated "      &
              &,"at time ",I6)', ntstep
      ENDIF
      DEALLOCATE ( ssgbdy , sgpbdy , stvbdy , smlbdy , dmlhed , STAT=nzerr )
    ENDIF
    ALLOCATE ( ssgbdy (mxsosg        , maxsgl) , STAT = nzerr )
    ALLOCATE ( sgpbdy (mxsogp        , maxgpl) , STAT = nzerr )
    ALLOCATE ( stvbdy (mxsotv, maxrtv, maxtvl) , STAT = nzerr )
    ALLOCATE ( smlbdy (mxsoml, maxmlv, maxmll) , STAT = nzerr )
    ALLOCATE ( dmlhed (mxsops        , maxmll) , STAT = nzerr )
    ssgbdy (:,:)   = rmdi
    sgpbdy (:,:)   = rmdi
    stvbdy (:,:,:) = rmdi
    smlbdy (:,:,:) = rmdi
    dmlhed (:,:)   = rmdi

! other actions
! -------------

    ntotqc = 0

! check if sorting of reports and quality control to be done at current timestep

    IF (lfirst) lqcall = .FALSE.
    IF (.NOT. lqcall) THEN
      lqcall = (lfirst) .OR. (rmod( ntstep*dt,dtqc,epsy ) < tmaxbox)
    ELSE
      lqcall = (rmod( ntstep*dt,dtqc,epsy ) < tminbox)
    ENDIF
    IF (dtqc <= tconbox+epsy) lqcall = .TRUE.

! for non-isotropic 1-dimensional horizontal correlations

    DO mv = 1 , mxispr
      zsdni (mv,1) = -sezispr * LOG( xezispr + (mv-1) *dezispr )
      zsdni (mv,2) =                 xthispr + (mv-1) *dthispr
    ENDDO
    IF ((lwonl) .AND. (lfirst)) THEN
      mxprspr = (mxispr + 6) / 7
      DO npspr = 1 , mxprspr
        nspe = MIN( INT( 7*npspr ,iintegers) , mxispr )
        nspa = nspe - 6
        WRITE( nupr,'(''sprp-lv['',F7.0,'','',F7.0,'']:'',7F7.1)' )            &
               zsdni(nspa,1), zsdni(nspe,1), (zsdni(mv,1), mv=nspa,nspe)
      ENDDO
    ENDIF

! flush YUPRINT file
!   IF ((lwonl) .AND. (lfirst) .AND. (ldump_ascii)) THEN
!     CLOSE (nupr)
!     OPEN  (nupr,FILE=yuprint,FORM='FORMATTED',STATUS='OLD'                   &
!                             ,POSITION='APPEND',IOSTAT=nstat)
!     IF (nstat /= 0) PRINT '("OPENING OF FILE yuprint FAILED, local_info_aux")'
!   ENDIF

  ENDIF

!-------------------------------------------------------------------------------
!  Section 3: Actions to be taken after computing the local information
!             on the observations and their location
!             (but prior to communication related to spatial consistency check)
!-------------------------------------------------------------------------------

  IF (nactio == 3) THEN

! get record with interpolated surface-level (and other single-level) reports
! for later printout for control

    IF ((lfirst) .OR. (rmod( ntstep*dt,c3600,epsy ) < tconbox)) THEN
      ista = 5
      nzerr = 0
      maxysu  =  maxsgo + maxmlo

      IF (ALLOCATED( oysu )) THEN
        PRINT '("CAUTION in src_sing_local: oysu is already allocated "        &
              &,"at time ",I6)', ntstep
        DEALLOCATE ( oysu , yysu , STAT=nzerr )
      ENDIF
      ALLOCATE ( oysu    (maxysu, 19) , STAT = nzerr )
      ALLOCATE ( yysu    (maxysu    ) , STAT = nzerr )
      ntotys = 0

      ndys   = 100
      IF (lfirst) ndys = 1
      DO ista = 1 , nsusta , ndys
        icps = isgadm(ista,7)
        icsu = isgadm(ista,8)
        DO itic = 1, 2
          lrepsu = (icsu > 0)
          IF (lrepsu) lrepsu = (MAX( zwtsu(icsu,itic,1) , zwtsu(icsu,itic,2)   &
                                   , zwtsu(icsu,itic,3) ) > epsy)
          lrepps = (icps > 0)
          IF (lrepps) lrepps = (omykps(icps,itic) > epsy)
          IF (((lrepsu) .OR. (lrepps)) .AND. (ntotys < maxysu)) THEN
            nsgob  = isgadm(ista,itic)
            ntotys = ntotys + 1
            DO nvys = 3 , 16
              oysu (ntotys,nvys) = rmdi
            ENDDO
            IF (lrepsu) THEN
              yysu (ntotys   )  =   ystidsu(icsu) (1:ilstidp)
              oysu (ntotys, 1)  =   osghed(nsgob,nhalt) - osghed(nsgob,nhsurf)
              oysu (ntotys, 7)  =   osgbdy(nsgob,nbsviu)
              oysu (ntotys, 8)  =   osgbdy(nsgob,nbsviv)
              oysu (ntotys, 9)  =   xoisu (icsu,itic,1)
              oysu (ntotys,10)  =   xoisu (icsu,itic,2)
              oysu (ntotys,12)  =   osgbdy(nsgob,nbsvit)
              oysu (ntotys,13)  =   xoisu (icsu,itic,3)
              oysu (ntotys,15)  =   osgbdy(nsgob,nbsviq)
              oysu (ntotys,16)  =   xoisu (icsu,itic,4)
            ENDIF
            IF (lrepps) THEN
              yysu (ntotys   )  =   ystidps(icps) (1:ilstidp)
              oysu (ntotys, 1)  =   osgbdy(nsgob,nbsviz)
              oysu (ntotys, 3)  = - osgbdy(nsgob,nbsvip) / 100.0_wp
              oysu (ntotys, 4)  =   zoips (icps,itic)    / 100.0_wp
            ENDIF
            oysu (ntotys, 2)  =   osgbdy(nsgob,nbsp ) / 100.0_wp
            oysu (ntotys, 5)  =   osgbdy(nsgob,nbsu )
            oysu (ntotys, 6)  =   osgbdy(nsgob,nbsv )
            oysu (ntotys,11)  =   osgbdy(nsgob,nbst )
            oysu (ntotys,14)  =   osgbdy(nsgob,nbsrh)
            oysu (ntotys,17)  =   osghed(nsgob,nhtime)
            oysu (ntotys,18)  =   REAL ( mosghd(nsgob,nhitot), wp )
            oysu (ntotys,19)  =   REAL ( mosghd(nsgob,nhjtot), wp )
          ENDIF
        ENDDO
      ENDDO
      DO ista = 1 , nmlsta
        icps = imladm(ista,7)
        DO itic = 1, 2
          lrepps = (icps > 0)
          IF (lrepps) lrepps = (omykps(icps,itic) > epsy)
          IF ((lrepps) .AND. (ntotys < maxysu)) THEN
            nmlob  = imladm(ista,itic)
            ntotys = ntotys + 1
            yysu (ntotys   )  =   ystidps(icps) (1:ilstidp)
            oysu (ntotys, 1)  =   omlhed(nmlob,nhtviz)
            oysu (ntotys, 2)  =   omlbdy(nmlob,1,nbsp) / 100.0_wp
            oysu (ntotys, 3)  = - omlhed(nmlob,nhtvip) / 100.0_wp
            oysu (ntotys, 4)  =   zoips (icps,itic)    / 100.0_wp
            DO nvys = 5 , 16
              oysu (ntotys,nvys) = rmdi
            ENDDO
            oysu (ntotys,17)  =   omlhed(nmlob,nhtime)
            oysu (ntotys,18)  =   REAL ( momlhd(nmlob,nhitot), wp )
            oysu (ntotys,19)  =   REAL ( momlhd(nmlob,nhjtot), wp )
          ENDIF
        ENDDO
      ENDDO
      DO ista = 1 , nuasta , ndys
        icua = isgadm(nsusta+ista,8)
        DO itic = 1, 2
          lrepsu = (icua > 0)
          IF (lrepsu) lrepsu = (MAX( zwtua(icua,itic,1) , zwtua(icua,itic,2)   &
                                   , zwtua(icua,itic,3) ) > epsy)
          IF ((lrepsu) .AND. (ntotys < maxysu)) THEN
            nsgob  = isgadm(nsusta+ista,itic)
            ntotys = ntotys + 1
            yysu (ntotys   ) = ystidua(icua) (1:ilstidp)
            oysu (ntotys, 1) = 0.0_wp
            oysu (ntotys, 2) = osgbdy(nsgob,nbsp ) / 100.0_wp
            oysu (ntotys, 3) = rmdi
            oysu (ntotys, 4) = rmdi
            oysu (ntotys, 5) = osgbdy(nsgob,nbsu  )
            oysu (ntotys, 6) = osgbdy(nsgob,nbsv  )
            oysu (ntotys, 7) = osgbdy(nsgob,nbsviu)
            oysu (ntotys, 8) = osgbdy(nsgob,nbsviv)
            oysu (ntotys, 9) = xoiua (icua,itic,1)
            oysu (ntotys,10) = xoiua (icua,itic,2)
            oysu (ntotys,11) = osgbdy(nsgob,nbst  )
            oysu (ntotys,12) = osgbdy(nsgob,nbsvit)
            oysu (ntotys,13) = xoiua (icua,itic,3)
            oysu (ntotys,14) = osgbdy(nsgob,nbsrh )
            oysu (ntotys,15) = osgbdy(nsgob,nbsviq)
            oysu (ntotys,16) = xoiua (icua,itic,4)
            oysu (ntotys,17) = osghed(nsgob,nhtime)
            oysu (ntotys,18) = REAL ( mosghd(nsgob,nhitot), wp )
            oysu (ntotys,19) = REAL ( mosghd(nsgob,nhjtot), wp )
          ENDIF
        ENDDO
      ENDDO
    ENDIF

  ENDIF

!-------------------------------------------------------------------------------
!  Section 4: Actions to be taken after the spatial consistency check
!-------------------------------------------------------------------------------

  IF (nactio == 4) THEN

! de-allocation of memory

    DEALLOCATE ( zoiciv  , STAT = nzerr )
    DEALLOCATE ( zqcfiv  , STAT = nzerr )
    DEALLOCATE ( zmodiv  , STAT = nzerr )
    DEALLOCATE ( zsativ  , STAT = nzerr )
    DEALLOCATE ( zactiv  , STAT = nzerr )
    DEALLOCATE ( ystidiv , STAT = nzerr )
    DEALLOCATE ( ioiv    , STAT = nzerr )
    DEALLOCATE ( joiv    , STAT = nzerr )
    DEALLOCATE ( ioiv_tot, STAT = nzerr )
    DEALLOCATE ( joiv_tot, STAT = nzerr )
    DEALLOCATE ( iqcliv  , STAT = nzerr )
    DEALLOCATE ( iqcfiv  , STAT = nzerr )
    DEALLOCATE ( kiobiv  , STAT = nzerr )
    DEALLOCATE ( kioiiv  , STAT = nzerr )
    DEALLOCATE ( ktypiv  , STAT = nzerr )

    DEALLOCATE ( imladm  , STAT = nzerr )
    DEALLOCATE ( isgadm  , STAT = nzerr )
    DEALLOCATE ( igpadm  , STAT = nzerr )
    DEALLOCATE ( itvadm  , STAT = nzerr )
    DEALLOCATE ( tmladm  , STAT = nzerr )
    DEALLOCATE ( tsgadm  , STAT = nzerr )
    DEALLOCATE ( tgpadm  , STAT = nzerr )
    DEALLOCATE ( ttvadm  , STAT = nzerr )

    DEALLOCATE ( wtml    , STAT = nzerr )
    DEALLOCATE ( wtsg    , STAT = nzerr )
    DEALLOCATE ( wtgp    , STAT = nzerr )
    DEALLOCATE ( wttv    , STAT = nzerr )
    DEALLOCATE ( zsdni   , STAT = nzerr )

  ENDIF

!-------------------------------------------------------------------------------
!  Section 5: Actions to be taken after computing the local information,
!             even if analysis increments are not computed at current timestep
!-------------------------------------------------------------------------------

  IF (nactio == 5) THEN

! de-allocation of memory for observation data record

    IF (ntstep == MIN( INT( MAX( nverend , nudgend ),iintegers) , nstop )) THEN

      DEALLOCATE   ( omlbdy , omlhed , momlbd , momlhd , yomlhd , STAT = nzerr )
      IF (.NOT. lsurfa) THEN
        DEALLOCATE ( osgbdy , osghed , mosgbd , mosghd , yosghd , STAT = nzerr )
      ENDIF
      DEALLOCATE   ( ogpbdy , ogphed , mogpbd , mogphd , yogphd , STAT = nzerr )
      DEALLOCATE   ( otvbdy , otvhed , motvbd , motvhd , yotvhd , STAT = nzerr )

    ENDIF

  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure local_info_aux
!-------------------------------------------------------------------------------

END SUBROUTINE local_info_aux


!-------------------------------------------------------------------------------
!+ Module procedure for production of sorted lists of currently used reports
!-------------------------------------------------------------------------------

SUBROUTINE local_sort_reports

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure produces sorted lists of multi-level, single-level,
!   and GPS stations with currently used reports. This allows to use temporal
!   linear interpolation as temporal nudging weights, which are also computed
!   here for the reports.
!
! Method:
!   Multi-level and single-level reports, and aircraft reports are treated
!   separately. 2 reports may be listed to the same station, if they are
!   assigned to the same horizontal grid point and have the same station
!   identity. 
!   If the selected temporal weighting is linear interpolation, then no more
!   than 2 reports per station will be used (per timestep).
!   A condition for 2 reports having the same station index in the list is
!   that the horizontal correlations (scales) must be the same for the 2
!   reports for each variable (which applies if temporal weighting is linear
!   interpolation or the selected correlation scale does not vary in time.
!   (Reason:) Only if this applies, the spreading procedures will be able to
!   process the 2 reports simultaneously).
!   If the selected temporal weighting is linear interpolation, weights accord-
!   ing to the saw-tooth-shaped weight function are additionally computed for
!   the reports. This allows to deal with incomplete vertical profiles, missing
!   data, or data rejected by the threshold quality control.
!
! Written by        :  Christoph Schraff, DWD  (original version: 20.10.97)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments: None
! --------------------

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    istafnd          ,& ! station index of the element in list '???adm' with the
                        ! same station as the current report
    nmlob  , nsgob   ,& ! (loop) indices of reports in the ODR (obs data record)
    nmlobe , nsgobe  ,& ! max. report indices in the ODR for loop over reports
    ngpob  , ngpobe  ,& ! loop index and max. report index in the ODR for GPS
    ntvob  , ntvobe  ,& ! loop index and max. report index in ODR for sat retrv.
    nmlstab, nsustab ,& ! prelimiary number of station elements in list '???adm'
    ngpstab, ntvstab ,& ! prelimiary number of station elements in list '???adm'
    kcdtyp           ,& ! observation code type
!   ksatid           ,& ! satellite identity
    i                ,& ! loop index for printout
    nsgpra , nsgprd  ,& ! first / stepping of sta. index in list for printout
    nstat               ! loop index for printout

  REAL    (KIND=wp   )     ::  &
    wtboxe           ,& ! largest model time [hrs] for which the temporal
                        ! weights should be valid within the current sorting box
    timdif           ,& ! time distance between the 2 reports at same station
    tabdif           ,& ! time distance between the obs. time and the model time
    wtuka  , wtuke   ,& ! influence radii for saw-tooth-shaped temporal weights
    tipmx            ,& ! modified max. time span (hrs) allowed for linear
                        ! interpolat. 0.5 h for frequent obs.
    ztobs1 , ztobs2  ,& ! report times (for printout)
    zpobpr              ! pressure [hPa] (for printout)

  LOGICAL                  ::  &
    lqcflag          ,& ! set threshold quality control flags for current report
    lactive          ,& ! current report is active
    lnotyet          ,& ! station is not yet in the list '???adm'
    lrhtvar          ,& ! correlation scale varies in time
    lprsort             ! print sorted lists to standard output

! Local (automatic) arrays: None
! -------------------------
!
!------------ End of header ----------------------------------------------------

 
!-------------------------------------------------------------------------------
! Begin Subroutine local_sort_reports
!-------------------------------------------------------------------------------

 
!-------------------------------------------------------------------------------
!  Section 0: Preset station counters and loop ranges
!-------------------------------------------------------------------------------

  nmlsta = 0
  nmlobe = ntotml
  nsusta = 0
  nuasta = 0
  nsgobe = ntotsg
  ngpsta = 0
  ngpobe = ntotgp
  ntvsta = 0
  ntvobe = ntottv
  wtboxe = aiwthr

!-------------------------------------------------------------------------------
!  Section 1: Sort multi-level reports from fixed stations (TEMP / PILOT)
!-------------------------------------------------------------------------------

  lrhtvar  =  (MAX( ABS( rhtfac(1) -c1 ) , ABS( rhtfac(2) -c1 )                &
                  , ABS( rhtfac(3) -c1 ) , ABS( rhtfac(4) -c1 ) ) > epsy)

loop_over_multi_level_reports:  DO nmlob = 1 , nmlobe

! modify temporal range of influence of frequent Profiler/RASS/VAD reports.
! influence towards the past: 0.5 hours, towards the future: 0.2 hours.
! (note: this must be the same as in Section 3 of routine 'obs_org_cdf_proc'
!   in module 'rc_obs_proc_cdf.f90' where the time window is also hard-coded
!   for these observation code types !)

  kcdtyp = momlhd(nmlob,nhcode)
  IF (     (kcdtyp == nwp_eu) .OR. (kcdtyp == nra_eu)                          &
      .OR. (kcdtyp == npr_us) .OR. (kcdtyp == nravad)) THEN
    wtuka = MIN( 0.5_wp, wtukrsa )
    wtuke = MIN( 0.2_wp, wtukrse )
    tipmx = 0.5_wp
  ELSE
    wtuka = wtukrsa
    wtuke = wtukrse
    tipmx = tipolmx
  ENDIF
  lactive = .FALSE.
  lqcflag = .FALSE.
  IF (     (momlhd(nmlob,nhobtp) == ntemp)                                     &
      .OR. (momlhd(nmlob,nhobtp) == npilot)) THEN
    lqcflag =       (lverif) .AND. (lqcall)                                    &
              .AND. (momlhd(nmlob,nhqcfw) <= -1)                               &
              .AND. (momlhd(nmlob,nhqcfw) > -99)                               &
              .AND. (ABS( c3600*( omlhed(nmlob,nhtime)-acthr)+epsy ) <= cqcwbox)
    lactive =       (omlhed(nmlob,nhtime) +MAX(wtuke ,tipmx ) > aiwthr)        &
              .AND. (omlhed(nmlob,nhtime) -MAX(wtuka ,tipmx ) < wtboxe)        &
              .AND. (omlhed(nmlob,nhtime) +wtuke  <= nudgend*dtdeh +epsy)      &
              .AND. (momlhd(nmlob,nhpass) == 0)                                &
              .AND. (MAX( momlhd(nmlob,nhuexi) *gnudg(1)                       &
                        , REAL(1- ABS( momlhd(nmlob,nhobtp)-ntemp ),wp) *gnudg(2)&
                        , momlhd(nmlob,nhtexi) *gnudg(3)                       &
                        , momlhd(nmlob,nhqexi) *gnudg(4) ) > epsy)
  ENDIF

  IF ((lactive) .OR. (lqcflag)) THEN
    lnotyet = .TRUE.
    io = momlhd(nmlob,nhio)
    jo = momlhd(nmlob,nhjo)
    ! avoid 'mixing' or 'replacement' between 1 active and 1 passive report
    ! only linear interpolation and 'replacement' between 2 active reports
    IF (.NOT. lactive) io = -io

!   IF ((ltipol) .AND. (lactive)) THEN
    ! check if current report can be added to an existing station in 'imladm':
    !     1.  the current report and the existing station must be active
    !        (to allow for simultaneous processing in the spreading procedures),
    ! and 2a. linear temporal interpolation is done for this 'obs type' (in this
    !         case, only 1 index can be used in 'imladm' for this station)
    ! or  2b. the temporal weight function does not vary in time and the
    !         appropriate slot at index 'ista' is still empty (in this case,
    !         more than 1 index could be used in 'imladm' for this station)
    IF ((lactive) .AND. ((ltipol) .OR. (.NOT. lrhtvar))) THEN
      DO ista = 1 , nmlsta
        IF (      (io == imladm(ista,3)) .AND. (jo == imladm(ista,4))          &
            .AND. (yomlhd(nmlob) == yomlhd(imladm(ista,5)))                    &
            .AND. (kcdtyp == momlhd(imladm(ista,5),nhcode))) THEN
          IF (     (ltipol)                                                    &
              .OR. (      (omlhed(nmlob,nhtime) > aiwthr+epsy)                 &
                    .AND. (imladm(ista,2) == 0))                               &
              .OR. (      (omlhed(nmlob,nhtime) <= aiwthr+epsy)                &
                    .AND. (imladm(ista,1) == 0))) THEN
            lnotyet = .FALSE.
            istafnd = ista
          ENDIF
        ENDIF
      ENDDO
    ENDIF

    IF (lnotyet) THEN
! a new station / report is added to the list 'imladm'
      nmlsta = nmlsta + 1
      imladm (nmlsta,3) = io
      imladm (nmlsta,4) = jo
      imladm (nmlsta,5) = nmlob
      imladm (nmlsta,6) = 0
      IF (omlhed(nmlob,nhtime) > aiwthr+epsy) THEN
        imladm (nmlsta,1) = 0
        imladm (nmlsta,2) = nmlob
        tmladm (nmlsta,1) = c0
        tmladm (nmlsta,2) = wtuka
        tmladm (nmlsta,3) = wtuka
      ELSE
        imladm (nmlsta,1) = nmlob
        imladm (nmlsta,2) = 0
        tmladm (nmlsta,1) = wtuke
        tmladm (nmlsta,2) = c0
        tmladm (nmlsta,3) = c0
      ENDIF

    ELSEIF (omlhed(nmlob,nhtime) > aiwthr+epsy) THEN
! the 'future' report is possibly added to a station existing in list 'imladm'
      IF (     (imladm(istafnd,2) == 0)                                        &
          .OR. (   omlhed(MAX(imladm(istafnd,2),i1),nhtime)                    &
                >= omlhed(nmlob,nhtime))) THEN
        imladm (istafnd,2) = nmlob
        tmladm (istafnd,2) = wtuka
        tmladm (istafnd,3) = wtuka
        IF ((ltipol) .AND. (imladm(istafnd,1) > 0)) THEN
          timdif = omlhed(nmlob,nhtime) - omlhed(imladm(istafnd,1),nhtime)
          IF (timdif <= tipmx +epsy) THEN
            imladm (istafnd,6) = 1
            tmladm (istafnd,1) = timdif
            tmladm (istafnd,2) = timdif
            tmladm (istafnd,3) = timdif
          ELSE
            tmladm (istafnd,3) = MIN( timdif , wtuka )
          ENDIF
        ENDIF
      ENDIF
    ELSE
! the 'past' report is possibly added to a station existing in list 'imladm'
      IF (     (imladm(istafnd,1) == 0)                                        &
          .OR. (   omlhed(MAX(imladm(istafnd,1),i1),nhtime)                    &
                <= omlhed(nmlob,nhtime))) THEN
        imladm (istafnd,1) = nmlob
        tmladm (istafnd,1) = wtuke
        IF ((ltipol) .AND. (imladm(istafnd,2) > 0)) THEN
          timdif = - omlhed(nmlob,nhtime) + omlhed(imladm(istafnd,2),nhtime)
          IF (timdif <= tipmx +epsy) THEN
            imladm (istafnd,6) = 1
            tmladm (istafnd,1) = timdif
            tmladm (istafnd,2) = timdif
            tmladm (istafnd,3) = timdif
          ELSE
            tmladm (istafnd,3) = MIN( timdif , wtuka )
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDIF

ENDDO loop_over_multi_level_reports

! IF (lverif) THEN
!   DO ista = 1 , nmlsta
!     imladm (ista,3) = ABS( imladm(ista,3) )
!   ENDDO
! ENDIF

! If - 2 reports have the same station index ('ista') in the list 'imladm'
!      because linear temporal interpolation is switched on, however
!    - the time distance betw. the 2 reports is too large for linear interpol.
! &  - the correlation scale for ps, (u,v), T, or RH varies in time
! then give the 2 reports 2 different station indices 'ista' (since the
! 2 reports cannot be processed simultaneously in the spreading procedures)
! ----------------------------------------------------------------------------

  IF ((ltipol) .AND. (lrhtvar)) THEN
    nmlstab = nmlsta
    DO ista = 1 , nmlstab
      IF ((imladm(ista,1) > 0) .AND. (imladm(ista,2) > 0)) THEN
        ! modify temporal range of influence of Profiler/RASS/VAD reports.
        kcdtyp = momlhd(imladm(ista,5),nhcode)
        IF (     (kcdtyp == nwp_eu) .OR. (kcdtyp == nra_eu)                    &
            .OR. (kcdtyp == npr_us) .OR. (kcdtyp == nravad)) THEN
          wtuka = MIN( 0.5_wp, wtukrsa )
          wtuke = MIN( 0.2_wp, wtukrse )
          tipmx = 0.5_wp
        ELSE
          wtuka = wtukrsa
          wtuke = wtukrse
          tipmx = tipolmx
        ENDIF
        timdif = omlhed(imladm(ista,2),nhtime) - omlhed(imladm(ista,1),nhtime)
        IF (timdif > tipmx +epsy) THEN
          imladm (ista,6) = 0
          IF (omlhed(imladm(ista,1),nhtime)+wtuke  <= aiwthr) THEN
            imladm (ista,1) = 0
            tmladm (ista,1) = c0
          ELSEIF (omlhed(imladm(ista,2),nhtime)-wtuka  >= wtboxe) THEN
            imladm (ista,2) = 0
            tmladm (ista,2) = c0
            tmladm (ista,3) = c0
          ELSE
            nmlsta = nmlsta + 1
            imladm (nmlsta,1) = 0
            imladm (nmlsta,2) = imladm(ista,2)
            imladm (nmlsta,3) = imladm(ista,3)
            imladm (nmlsta,4) = imladm(ista,4)
            imladm (nmlsta,5) = imladm(ista,5)
            imladm (nmlsta,6) = 0
            tmladm (nmlsta,2) = tmladm(ista,2)
            tmladm (nmlsta,3) = tmladm(ista,3)
            imladm (ista  ,2) = 0
            tmladm (ista  ,2) = c0
            tmladm (ista  ,3) = c0
          ENDIF
        ENDIF
      ENDIF
    ENDDO
  ENDIF


!-------------------------------------------------------------------------------
!  Section 2: Sort multi-level reports from moving stations (AIRCRAFT)
!-------------------------------------------------------------------------------

loop_over_mult_aircraft_reports:  DO nmlob = 1 , nmlobe

  lactive = .FALSE.
  lqcflag = .FALSE.
  IF (momlhd(nmlob,nhobtp) == nairep) THEN
    lqcflag =       (lverif) .AND. (lqcall)                                    &
              .AND. (momlhd(nmlob,nhqcfw) <= -1)                               &
              .AND. (momlhd(nmlob,nhqcfw) > -99)                               &
              .AND. (ABS( c3600*( omlhed(nmlob,nhtime)-acthr)+epsy ) <= cqcwbox)
    lactive =       (omlhed(nmlob,nhtime) +wtukare > aiwthr)                   &
              .AND. (omlhed(nmlob,nhtime) -wtukara < wtboxe)                   &
              .AND. (omlhed(nmlob,nhtime) +wtukare <= nudgend*dtdeh +epsy)     &
              .AND. (momlhd(nmlob,nhpass) == 0)                                &
              .AND. (    (     (momlhd(nmlob,nhcode) /= nmodes)                &
                          .AND.(MAX( momlhd(nmlob,nhuexi) *gnudgar(1)          &
                                   , momlhd(nmlob,nhtexi) *gnudgar(3)          &
                                   , momlhd(nmlob,nhqexi) *gnudgar(4)) > epsy))&
                     .OR.(     (momlhd(nmlob,nhcode) == nmodes)                &
                          .AND.(MAX( momlhd(nmlob,nhuexi) *gnudgms(1)          &
                                   , momlhd(nmlob,nhtexi) *gnudgms(3)          &
                                   , momlhd(nmlob,nhqexi) *gnudgms(4)) > epsy)))
  ENDIF
  IF ((lactive) .OR. (lqcflag)) THEN
    io = momlhd(nmlob,nhio)
    jo = momlhd(nmlob,nhjo)
    IF (.NOT. lactive) io = -io
    nmlsta = nmlsta + 1
    imladm (nmlsta,3) = io
    imladm (nmlsta,4) = jo
    imladm (nmlsta,5) = nmlob
    IF (omlhed(nmlob,nhtime) > aiwthr+epsy) THEN
      imladm (nmlsta,1) = 0
      imladm (nmlsta,2) = nmlob
      tmladm (nmlsta,1) = c0
      tmladm (nmlsta,2) = wtukara
      tmladm (nmlsta,3) = wtukara
    ELSE
      imladm (nmlsta,1) = nmlob
      imladm (nmlsta,2) = 0
      tmladm (nmlsta,1) = wtukare
      tmladm (nmlsta,2) = c0
      tmladm (nmlsta,3) = c0
    ENDIF
  ENDIF

ENDDO loop_over_mult_aircraft_reports


!-------------------------------------------------------------------------------
!  Section 3: Sort surface-level reports from 'fixed' stations (SYNOP/DRIBU/TEMP
!-------------------------------------------------------------------------------

  lrhtvar  =  (MAX( ABS( rhtfsu(1) -c1 ) , ABS( rhtfsu(2) -c1 )                &
                  , ABS( rhtfsu(3) -c1 ) , ABS( rhtfsu(4) -c1 ) ) > epsy)

loop_over_single_level_reports:  DO nsgob = 1 , nsgobe

  lactive = .FALSE.
  lqcflag = .FALSE.
  IF (      (mosghd(nsgob,nhobtp) /= nairep)                                   &
      .AND. (mosghd(nsgob,nhobtp) /= nsatob)) THEN
    lqcflag =       (lverif) .AND. (lqcall)                                    &
              .AND. (mosghd(nsgob,nhqcfw) <= -1)                               &
              .AND. (mosghd(nsgob,nhqcfw) > -99)                               &
              .AND. (ABS( c3600*( osghed(nsgob,nhtime)-acthr)+epsy ) <= cqcwbox)
    lactive =       (osghed(nsgob,nhtime) +MAX(wtuksue,tipmxsu) > aiwthr)      &
              .AND. (osghed(nsgob,nhtime) -MAX(wtuksua,tipmxsu) < wtboxe)      &
              .AND. (osghed(nsgob,nhtime) +wtuksue <= nudgend*dtdeh +epsy)     &
              .AND. (mosghd(nsgob,nhpass) == 0)                                &
              .AND. (     (      (BTEST( mosgbd(nsgob,nbserr),nvru ))          &
                           .AND. (gnudgsu(1) > epsy))                          &
                     .OR. (      (BTEST( mosgbd(nsgob,nbserr),nvrz ))          &
                           .AND. (gnudgsu(2) > epsy))                          &
                     .OR. (      (BTEST( mosgbd(nsgob,nbserr),nvrt ))          &
                           .AND. (gnudgsu(3) > epsy))                          &
                     .OR. (      (BTEST( mosgbd(nsgob,nbserr),nvrq ))          &
                           .AND. (gnudgsu(4) > epsy)))
  ENDIF

  IF ((lactive) .OR. (lqcflag)) THEN
    lnotyet = .TRUE.
    io = mosghd(nsgob,nhio)
    jo = mosghd(nsgob,nhjo)
    ! avoid 'mixing' or 'replacement' between 1 active and 1 passive report
    ! only linear interpolation and 'replacement' between 2 active reports
    IF (.NOT. lactive) io = -io

!   IF ((ltipsu) .AND. (lactive)) THEN
    ! check if current report can be added to an existing station in 'isgadm':
    !     1.  the current report and the existing station must be active
    !        (to allow for simultaneous processing in the spreading procedures),
    ! and 2a. linear temporal interpolation is done for this 'obs type' (in this
    !         case, only 1 index can be used in 'isgadm' for this station)
    ! or  2b. the temporal weight function does not vary in time and the
    !         appropriate slot at index 'ista' is still empty (in this case,
    !         more than 1 index could be used in 'isgadm' for this station)
    IF ((lactive) .AND. ((ltipsu) .OR. (.NOT. lrhtvar))) THEN
      DO ista = 1 , nsusta
        IF (io == isgadm(ista,3)) THEN
          IF (      (jo == isgadm(ista,4))                                     &
              .AND. (yosghd(nsgob) == yosghd(isgadm(ista,5)))                  &
              .AND. (mosghd(nsgob,nhcode) == mosghd(isgadm(ista,5),nhcode))    &
              .AND. (NINT( osghed(nsgob,nhalt)                                 &
                          -osghed(isgadm(ista,5),nhalt) ) == 0)) THEN
                    ! (Note: accuracy of < 1m is consistent with
                    !        the 10Pa check in 'obs_redundancy')
            IF (     (ltipsu)                                                  &
                .OR. (      (osghed(nsgob,nhtime) > aiwthr+epsy)               &
                      .AND. (isgadm(ista,2) == 0))                             &
                .OR. (      (osghed(nsgob,nhtime) <= aiwthr+epsy)              &
                      .AND. (isgadm(ista,1) == 0))) THEN
              lnotyet = .FALSE.
              istafnd = ista
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDIF

    IF (lnotyet) THEN
! a new station / report is added to the list 'isgadm'
      nsusta = nsusta + 1
      isgadm (nsusta,3) = io
      isgadm (nsusta,4) = jo
      isgadm (nsusta,5) = nsgob
      isgadm (nsusta,6) = 0
      IF (osghed(nsgob,nhtime) > aiwthr+epsy) THEN
        isgadm (nsusta,1) = 0
        isgadm (nsusta,2) = nsgob
        tsgadm (nsusta,1) = c0
        tsgadm (nsusta,2) = wtuksua
        tsgadm (nsusta,3) = wtuksua
      ELSE
        isgadm (nsusta,1) = nsgob
        isgadm (nsusta,2) = 0
        tsgadm (nsusta,1) = wtuksue
        tsgadm (nsusta,2) = c0
        tsgadm (nsusta,3) = c0
      ENDIF

    ELSEIF (osghed(nsgob,nhtime) > aiwthr+epsy) THEN
! the 'future' report is possibly added to a station existing in list 'isgadm'
      IF (     (isgadm(istafnd,2) == 0)                                        &
          .OR. (   osghed(MAX(isgadm(istafnd,2),i1),nhtime)                    &
                >= osghed(nsgob,nhtime))) THEN
        isgadm (istafnd,2) = nsgob
        tsgadm (istafnd,2) = wtuksua
        tsgadm (istafnd,3) = wtuksua
        IF ((ltipsu) .AND. (isgadm(istafnd,1) > 0)) THEN
          timdif = osghed(nsgob,nhtime) - osghed(isgadm(istafnd,1),nhtime)
          IF (timdif <= tipmxsu+epsy) THEN
            isgadm (istafnd,6) = 1
            tsgadm (istafnd,1) = timdif
            tsgadm (istafnd,2) = timdif
            tsgadm (istafnd,3) = timdif
          ELSE
            tsgadm (istafnd,3) = MIN( timdif , wtuksua )
          ENDIF
        ENDIF
      ENDIF
    ELSE
! the 'past' report is possibly added to a station existing in list 'isgadm'
      IF (     (isgadm(istafnd,1) == 0)                                        &
          .OR. (   osghed(MAX(isgadm(istafnd,1),i1),nhtime)                    &
                <= osghed(nsgob,nhtime))) THEN
        isgadm (istafnd,1) = nsgob
        tsgadm (istafnd,1) = wtuksue
        IF ((ltipsu) .AND. (isgadm(istafnd,2) > 0)) THEN
          timdif = - osghed(nsgob,nhtime) + osghed(isgadm(istafnd,2),nhtime)
          IF (timdif <= tipmxsu+epsy) THEN
            isgadm (istafnd,6) = 1
            tsgadm (istafnd,1) = timdif
            tsgadm (istafnd,2) = timdif
            tsgadm (istafnd,3) = timdif
          ELSE
            tsgadm (istafnd,3) = MIN( timdif , wtuksua )
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDIF

ENDDO loop_over_single_level_reports

! IF (lverif) THEN
!   DO ista = 1 , nsusta
!     isgadm (ista,3) = ABS( isgadm(ista,3) )
!   ENDDO
! ENDIF

! If - 2 reports have the same station index ('ista') in the list 'isgadm'
!      because linear temporal interpolation is switched on, however
!    - the time distance betw. the 2 reports is too large for linear interpol.
! &  - the correlation scale for ps, (u,v), T, or RH varies in time
! then give the 2 reports 2 different station indices 'ista' (since the
! 2 reports cannot be processed simultaneously in the spreading procedures)
! ----------------------------------------------------------------------------

  IF ((ltipsu) .AND. (lrhtvar)) THEN
    nsustab = nsusta
    DO ista = 1 , nsustab
      IF ((isgadm(ista,1) > 0) .AND. (isgadm(ista,2) > 0)) THEN
        timdif = osghed(isgadm(ista,2),nhtime) - osghed(isgadm(ista,1),nhtime)
        IF (timdif > tipmxsu+epsy) THEN
          isgadm (ista,6) = 0
          IF (osghed(isgadm(ista,1),nhtime)+wtuksue <= aiwthr) THEN
            isgadm (ista,1) = 0
            tsgadm (ista,1) = c0
          ELSEIF (osghed(isgadm(ista,2),nhtime)-wtuksua >= wtboxe) THEN
            isgadm (ista,2) = 0
            tsgadm (ista,2) = c0
            tsgadm (ista,3) = c0
          ELSE
            nsusta = nsusta + 1
            isgadm (nsusta,1) = 0
            isgadm (nsusta,2) = isgadm(ista,2)
            isgadm (nsusta,3) = isgadm(ista,3)
            isgadm (nsusta,4) = isgadm(ista,4)
            isgadm (nsusta,5) = isgadm(ista,5)
            isgadm (nsusta,6) = 0
            tsgadm (nsusta,2) = tsgadm(ista,2)
            tsgadm (nsusta,3) = tsgadm(ista,3)
            isgadm (ista  ,2) = 0
            tsgadm (ista  ,2) = c0
            tsgadm (ista  ,3) = c0
          ENDIF
        ENDIF
      ENDIF
    ENDDO
  ENDIF


!-------------------------------------------------------------------------------
!  Section 4: Sort GPS reports
!-------------------------------------------------------------------------------

  lrhtvar  =  (ABS( rhtfac(4) -c1 ) > epsy)
  wtuka = MAX( wtukrsa, tipolmx, wtuksua, tipmxsu )
  wtuke = MAX( wtukrse,          wtuksue          )

loop_over_gps_reports:  DO ngpob = 1 , ngpobe

  ! correction of known bug when GPS reports have been read from COST ASCII file
  IF (      (.NOT. BTEST( mogpbd(ngpob,nbgerr),nvriwv ))                       &
      .AND. (      BTEST( mogpbd(ngpob,nbgerr),nvru   )))                      &
    mogpbd (ngpob,nbgerr) = IBSET ( mogpbd(ngpob,nbgerr), nvriwv )

  lactive = .FALSE.
  lqcflag = .FALSE.
  lqcflag =       (lverif) .AND. (lqcall)                                      &
            .AND. (mogphd(ngpob,nhqcfw) <= -1)                                 &
            .AND. (mogphd(ngpob,nhqcfw) > -99)                                 &
            .AND. (ABS( c3600*( ogphed(ngpob,nhtime)-acthr)+epsy ) <= cqcwbox)
! lqcflag =       (lverif) .AND. (lqcall)                                      &
!           .AND. (ABS( c3600*( ogphed(ngpob,nhtime)-acthr)+epsy ) <= cqcwbox)
! ! GPS data must go into IWV spatial consistency check whenever lqcall =.TRUE.
! ! (time window does reasonably but not completely cover the time range
! !  for radiosonde plus the time window in the IWV check for the GPS data)
! lqcflag =       ((lqcflag) .OR. ((lqcall) .AND. (liwvssc)))                  &
!           .AND. (ogphed(ngpob,nhtime) +wtuke+c05*c05*dtchkps > aiwthr)       &
!           .AND. (ogphed(ngpob,nhtime) -wtuka-c05*c05*dtchkps < wtboxe)
  ! at least 1 GPS obs from all stations should go into IWV spatial
  ! consistency check whenever lqcall =.TRUE.
  ! --> for data frequency 15 min, time window +/- 12 min is enough)
  lqcflag = ((lqcflag) .OR. (      (lqcall) .AND. (liwvssc)))                  &
                             .AND. (ogphed(ngpob,nhtime) +.2_wp > aiwthr)  &
                             .AND. (ogphed(ngpob,nhtime) -.2_wp < wtboxe)
  lactive =       (ogphed(ngpob,nhtime) +MAX(wtuksue,tipmxsu) > aiwthr)        &
            .AND. (ogphed(ngpob,nhtime) -MAX(wtuksua,tipmxsu) < wtboxe)        &
            .AND. (ogphed(ngpob,nhtime) +wtuksue <= nudgend*dtdeh +epsy)       &
            .AND. (mogphd(ngpob,nhpass) == 0) .AND. (gnudggp > epsy)           &
            .AND. (BTEST( mogpbd(ngpob,nbgerr),nvriwv ))
!           .AND. ((ogpbdy(ngpob,nbgtze) > rmdich) .AND. (gnudggp > epsy))

  IF ((lactive) .OR. (lqcflag)) THEN
    lnotyet = .TRUE.
    io = mogphd(ngpob,nhio)
    jo = mogphd(ngpob,nhjo)
    ! avoid 'mixing' or 'replacement' between 1 active and 1 passive report
    ! only linear interpolation and 'replacement' between 2 active reports
    IF (.NOT. lactive) io = -io

!   IF ((ltipsu) .AND. (lactive)) THEN
    ! check if current report can be added to an existing station in 'igpadm':
    !     1.  the current report and the existing station must be active
    !        (to allow for simultaneous processing in the spreading procedures),
    ! and 2a. linear temporal interpolation is done for this 'obs type' (in this
    !         case, only 1 index can be used in 'igpadm' for this station)
    ! or  2b. the temporal weight function does not vary in time and the
    !         appropriate slot at index 'ista' is still empty (in this case,
    !         more than 1 index could be used in 'igpadm' for this station)
    IF ((lactive) .AND. ((ltipsu) .OR. (.NOT. lrhtvar))) THEN
      DO ista = 1 , ngpsta
        IF (io == igpadm(ista,3)) THEN
          IF (     (jo == igpadm(ista,4))                                      &
             .AND. (yogphd(ngpob) == yogphd(igpadm(ista,5)))                   &
             .AND. (mogphd(ngpob,nhcode) == mogphd(igpadm(ista,5),nhcode))) THEN
                   ! equal code type also means equal GPS processing centre
            IF (     (ltipsu)                                                  &
                .OR. (      (ogphed(ngpob,nhtime) > aiwthr+epsy)               &
                      .AND. (igpadm(ista,2) == 0))                             &
                .OR. (      (ogphed(ngpob,nhtime) <= aiwthr+epsy)              &
                      .AND. (igpadm(ista,1) == 0))) THEN
              lnotyet = .FALSE.
              istafnd = ista
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDIF

    IF (lnotyet) THEN
! a new station / report is added to the list 'igpadm'
      ngpsta = ngpsta + 1
      igpadm (ngpsta,3) = io
      igpadm (ngpsta,4) = jo
      igpadm (ngpsta,5) = ngpob
      igpadm (ngpsta,6) = 0
      IF (ogphed(ngpob,nhtime) > aiwthr+epsy) THEN
        igpadm (ngpsta,1) = 0
        igpadm (ngpsta,2) = ngpob
        tgpadm (ngpsta,1) = c0
        tgpadm (ngpsta,2) = wtuksua
        tgpadm (ngpsta,3) = wtuksua
      ELSE
        igpadm (ngpsta,1) = ngpob
        igpadm (ngpsta,2) = 0
        tgpadm (ngpsta,1) = wtuksue
        tgpadm (ngpsta,2) = c0
        tgpadm (ngpsta,3) = c0
      ENDIF

    ELSEIF (ogphed(ngpob,nhtime) > aiwthr+epsy) THEN
! the 'future' report is possibly added to a station existing in list 'igpadm'
      IF (     (igpadm(istafnd,2) == 0)                                        &
          .OR. (   ogphed(MAX(igpadm(istafnd,2),i1),nhtime)                    &
                >= ogphed(ngpob,nhtime))) THEN
        igpadm (istafnd,2) = ngpob
        tgpadm (istafnd,2) = wtuksua
        tgpadm (istafnd,3) = wtuksua
        IF ((ltipsu) .AND. (igpadm(istafnd,1) > 0)) THEN
          timdif = ogphed(ngpob,nhtime) - ogphed(igpadm(istafnd,1),nhtime)
          IF (timdif <= tipmxsu+epsy) THEN
            igpadm (istafnd,6) = 1
            tgpadm (istafnd,1) = timdif
            tgpadm (istafnd,2) = timdif
            tgpadm (istafnd,3) = timdif
          ELSE
            tgpadm (istafnd,3) = MIN( timdif , wtuksua )
          ENDIF
        ENDIF
      ENDIF
    ELSE
! the 'past' report is possibly added to a station existing in list 'igpadm'
      IF (     (igpadm(istafnd,1) == 0)                                        &
          .OR. (   ogphed(MAX(igpadm(istafnd,1),i1),nhtime)                    &
                <= ogphed(ngpob,nhtime))) THEN
        igpadm (istafnd,1) = ngpob
        tgpadm (istafnd,1) = wtuksue
        IF ((ltipsu) .AND. (igpadm(istafnd,2) > 0)) THEN
          timdif = - ogphed(ngpob,nhtime) + ogphed(igpadm(istafnd,2),nhtime)
          IF (timdif <= tipmxsu+epsy) THEN
            igpadm (istafnd,6) = 1
            tgpadm (istafnd,1) = timdif
            tgpadm (istafnd,2) = timdif
            tgpadm (istafnd,3) = timdif
          ELSE
            tgpadm (istafnd,3) = MIN( timdif , wtuksua )
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDIF

ENDDO loop_over_gps_reports

! if (lqcflag) and not (lactive) then get rid of sign of io = igpadm(ista,3)
! IF (((lverif) .OR. (liwvssc)) .AND. (lqcall)) THEN
!   DO ista = 1 , ngpsta
!     igpadm (ista,3) = ABS( igpadm(ista,3) )
!   ENDDO
! ENDIF

! If - 2 reports have the same station index ('ista') in the list 'igpadm'
!      because linear temporal interpolation is switched on, however
!    - the time distance betw. the 2 reports is too large for linear interpol.
! &  - the correlation scale for ps, (u,v), T, or RH varies in time
! then give the 2 reports 2 different station indices 'ista' (since the
! 2 reports cannot be processed simultaneously in the spreading procedures)
! ----------------------------------------------------------------------------

  IF ((ltipsu) .AND. (lrhtvar)) THEN
    ngpstab = ngpsta
    DO ista = 1 , ngpstab
      IF ((igpadm(ista,1) > 0) .AND. (igpadm(ista,2) > 0)) THEN
        timdif = ogphed(igpadm(ista,2),nhtime) - ogphed(igpadm(ista,1),nhtime)
        IF (timdif > tipmxsu+epsy) THEN
          igpadm (ista,6) = 0
          IF (ogphed(igpadm(ista,1),nhtime)+wtuksue <= aiwthr) THEN
            igpadm (ista,1) = 0
            tgpadm (ista,1) = c0
          ELSEIF (ogphed(igpadm(ista,2),nhtime)-wtuksua >= wtboxe) THEN
            igpadm (ista,2) = 0
            tgpadm (ista,2) = c0
            tgpadm (ista,3) = c0
          ELSE
            ngpsta = ngpsta + 1
            igpadm (ngpsta,1) = 0
            igpadm (ngpsta,2) = igpadm(ista,2)
            igpadm (ngpsta,3) = igpadm(ista,3)
            igpadm (ngpsta,4) = igpadm(ista,4)
            igpadm (ngpsta,5) = igpadm(ista,5)
            igpadm (ngpsta,6) = 0
            tgpadm (ngpsta,2) = tgpadm(ista,2)
            tgpadm (ngpsta,3) = tgpadm(ista,3)
            igpadm (ista  ,2) = 0
            tgpadm (ista  ,2) = c0
            tgpadm (ista  ,3) = c0
          ENDIF
        ENDIF
      ENDIF
    ENDDO
  ENDIF


!-------------------------------------------------------------------------------
!  Section 5: Sort satellite retrievals
!-------------------------------------------------------------------------------

  lrhtvar  =  (MAX( ABS( rhtfac(1) -c1 ) , ABS( rhtfac(2) -c1 )                &
                  , ABS( rhtfac(3) -c1 ) , ABS( rhtfac(4) -c1 ) ) > epsy)

loop_over_sat_retrievals:  DO ntvob = 1 , ntvobe

  kcdtyp  =  motvhd(ntvob,nhcode)
! ksatid  =  motvhd(ntvob,nhstid)
! wtuka   =  ssat(kidsat(ksatid))% wtuka
! wtuke   =  ssat(kidsat(ksatid))% wtuke
! tipmx   =  ssat(kidsat(ksatid))% tipmx
  wtuka   =  otvhed(ntvob,nh1wta)
  wtuke   =  otvhed(ntvob,nh1wte)
  tipmx   =  otvhed(ntvob,nh1tip)
  lactive = .FALSE.
  lqcflag = .FALSE.
  lqcflag =       (lverif) .AND. (lqcall)                                      &
            .AND. (ABS( c3600*( otvhed(ntvob,nhtime)-acthr)+epsy ) <= cqcwbox)
  lactive =       (otvhed(ntvob,nhtime) +MAX(wtuke,tipmx) > aiwthr)            &
            .AND. (otvhed(ntvob,nhtime) -MAX(wtuka,tipmx) < wtboxe)            &
            .AND. (otvhed(ntvob,nhtime) +wtuke <= nudgend*dtdeh +epsy)         &
            .AND. (motvhd(ntvob,nhpass) == 0)                                  &
            .AND. (MAX( motvhd(ntvob,nhuexi) *gnudgtv(1)                       &
                      , motvhd(ntvob,nhtexi) *gnudgtv(3)                       &
                      , motvhd(ntvob,nhqexi) *gnudgtv(4) ) > epsy)

  IF ((lactive) .OR. (lqcflag)) THEN
    lnotyet = .TRUE.
    io = motvhd(ntvob,nhio)
    jo = motvhd(ntvob,nhjo)
    ! avoid 'mixing' or 'replacement' between 1 active and 1 passive report
    ! only linear interpolation and 'replacement' between 2 active reports
    IF (.NOT. lactive) io = -io

!   IF ((ltipol) .AND. (lactive)) THEN
    ! check if current report can be added to an existing station in 'itvadm':
    !     1.  the current report and the existing station must be active
    !        (to allow for simultaneous processing in the spreading procedures),
    ! and 2a. linear temporal interpolation is done for this 'obs type' (in this
    !         case, only 1 index can be used in 'itvadm' for this station)
    ! or  2b. the temporal weight function does not vary in time and the
    !         appropriate slot at index 'ista' is still empty (in this case,
    !         more than 1 index could be used in 'itvadm' for this station)
    IF ((lactive) .AND. ((tipmx > epsy) .OR. (.NOT. lrhtvar))) THEN
      DO ista = 1 , ntvsta
        IF (io == itvadm(ista,3)) THEN
          IF (      (jo == itvadm(ista,4))                                     &
              .AND. (yotvhd(ntvob) == yotvhd(itvadm(ista,5)))                  &
              .AND. (kcdtyp == motvhd(itvadm(ista,5),nhcode))) THEN
            IF (     (tipmx > epsy)                                            &
                .OR. (      (otvhed(ntvob,nhtime) > aiwthr+epsy)               &
                      .AND. (itvadm(ista,2) == 0))                             &
                .OR. (      (otvhed(ntvob,nhtime) <= aiwthr+epsy)              &
                      .AND. (itvadm(ista,1) == 0))) THEN
              lnotyet = .FALSE.
              istafnd = ista
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDIF

    IF (lnotyet) THEN
! a new station / report is added to the list 'itvadm'
      ntvsta = ntvsta + 1
      itvadm (ntvsta,3) = io
      itvadm (ntvsta,4) = jo
      itvadm (ntvsta,5) = ntvob
      itvadm (ntvsta,6) = 0
      IF (otvhed(ntvob,nhtime) > aiwthr+epsy) THEN
        itvadm (ntvsta,1) = 0
        itvadm (ntvsta,2) = ntvob
        ttvadm (ntvsta,1) = c0
        ttvadm (ntvsta,2) = wtuka
        ttvadm (ntvsta,3) = wtuka
      ELSE
        itvadm (ntvsta,1) = ntvob
        itvadm (ntvsta,2) = 0
        ttvadm (ntvsta,1) = wtuke
        ttvadm (ntvsta,2) = c0
        ttvadm (ntvsta,3) = c0
      ENDIF

    ELSEIF (otvhed(ntvob,nhtime) > aiwthr+epsy) THEN
! the 'future' report is possibly added to a station existing in list 'itvadm'
      IF (     (itvadm(istafnd,2) == 0)                                        &
          .OR. (   otvhed(MAX(itvadm(istafnd,2),i1),nhtime)                    &
                >= otvhed(ntvob,nhtime))) THEN
        itvadm (istafnd,2) = ntvob
        ttvadm (istafnd,2) = wtuka
        ttvadm (istafnd,3) = wtuka
        IF ((tipmx > epsy) .AND. (itvadm(istafnd,1) > 0)) THEN
          timdif = otvhed(ntvob,nhtime) - otvhed(itvadm(istafnd,1),nhtime)
          IF (timdif <= tipmx+epsy) THEN
            itvadm (istafnd,6) = 1
            ttvadm (istafnd,1) = timdif
            ttvadm (istafnd,2) = timdif
            ttvadm (istafnd,3) = timdif
          ELSE
            ttvadm (istafnd,3) = MIN( timdif , wtuka )
          ENDIF
        ENDIF
      ENDIF
    ELSE
! the 'past' report is possibly added to a station existing in list 'itvadm'
      IF (     (itvadm(istafnd,1) == 0)                                        &
          .OR. (   otvhed(MAX(itvadm(istafnd,1),i1),nhtime)                    &
                <= otvhed(ntvob,nhtime))) THEN
        itvadm (istafnd,1) = ntvob
        ttvadm (istafnd,1) = wtuke
        IF ((tipmx > epsy) .AND. (itvadm(istafnd,2) > 0)) THEN
          timdif = - otvhed(ntvob,nhtime) + otvhed(itvadm(istafnd,2),nhtime)
          IF (timdif <= tipmx+epsy) THEN
            itvadm (istafnd,6) = 1
            ttvadm (istafnd,1) = timdif
            ttvadm (istafnd,2) = timdif
            ttvadm (istafnd,3) = timdif
          ELSE
            ttvadm (istafnd,3) = MIN( timdif , wtuka )
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDIF

ENDDO loop_over_sat_retrievals

! IF (lverif) THEN
!   DO ista = 1 , ntvsta
!     itvadm (ista,3) = ABS( itvadm(ista,3) )
!   ENDDO
! ENDIF


! If - 2 reports have the same station index ('ista') in the list 'itvadm'
!      because linear temporal interpolation is switched on, however
!    - the time distance betw. the 2 reports is too large for linear interpol.
! &  - the correlation scale for ps, (u,v), T, or RH varies in time
! then give the 2 reports 2 different station indices 'ista' (since the
! 2 reports cannot be processed simultaneously in the spreading procedures)
! ----------------------------------------------------------------------------

! IF ((ltipol) .AND. (lrhtvar)) THEN
! IF ((tipmx > epsy) .AND. (lrhtvar)) THEN
  IF (lrhtvar) THEN
    ntvstab = ntvsta
    DO ista = 1 , ntvstab
      IF ((itvadm(ista,1) > 0) .AND. (itvadm(ista,2) > 0)) THEN
!       ksatid =  motvhd(itvadm(ista,5),nhstid)
!       wtuka  =  ssat(kidsat(ksatid))% wtuka
!       wtuke  =  ssat(kidsat(ksatid))% wtuke
!       tipmx  =  ssat(kidsat(ksatid))% tipmx
        wtuka  =  otvhed(itvadm(ista,5),nh1wta)
        wtuke  =  otvhed(itvadm(ista,5),nh1wte)
        tipmx  =  otvhed(itvadm(ista,5),nh1tip)
        timdif = otvhed(itvadm(ista,2),nhtime) - otvhed(itvadm(ista,1),nhtime)
!       IF (timdif > tipmx+epsy) THEN
        IF ((tipmx > epsy) .AND. (timdif > tipmx+epsy)) THEN
          itvadm (ista,6) = 0
          IF (otvhed(itvadm(ista,1),nhtime)+wtuke <= aiwthr) THEN
            itvadm (ista,1) = 0
            ttvadm (ista,1) = c0
          ELSEIF (otvhed(itvadm(ista,2),nhtime)-wtuka >= wtboxe) THEN
            itvadm (ista,2) = 0
            ttvadm (ista,2) = c0
            ttvadm (ista,3) = c0
          ELSE
            ntvsta = ntvsta + 1
            itvadm (ntvsta,1) = 0
            itvadm (ntvsta,2) = itvadm(ista,2)
            itvadm (ntvsta,3) = itvadm(ista,3)
            itvadm (ntvsta,4) = itvadm(ista,4)
            itvadm (ntvsta,5) = itvadm(ista,5)
            itvadm (ntvsta,6) = 0
            ttvadm (ntvsta,2) = ttvadm(ista,2)
            ttvadm (ntvsta,3) = ttvadm(ista,3)
            itvadm (ista  ,2) = 0
            ttvadm (ista  ,2) = c0
            ttvadm (ista  ,3) = c0
          ENDIF
        ENDIF
      ENDIF
    ENDDO
  ENDIF


!-------------------------------------------------------------------------------
!  Section 6: Sort upper-air single-level reports from moving stations 
!             (AIRCRAFT/SATOB)
!-------------------------------------------------------------------------------

  nuaex = 0

loop_over_sing_aircraft_reports:  DO nsgob = 1 , nsgobe

  lactive = .FALSE.
  lqcflag = .FALSE.
  IF (     (mosghd(nsgob,nhobtp) == nairep)                                    &
      .OR. (mosghd(nsgob,nhobtp) == nsatob)) THEN
    ystid = yosghd(nsgob) (1:ilstidp)
    IF ((lwonl) .AND. (lfirst))                                                &
     WRITE( nupr,'(''airep '',2I5,I4,I5,F6.3,2X,A,4F6.1,F5.2,F4.1,F7.4,F8.0)') &
            nuasta, nsgob, mosghd(nsgob,nhpass), nudgend, dtdeh                &
          , ystid, osghed(nsgob,nhtime), wtukare, wtukara, aiwthr, wtboxe      &
          , osgbdy(nsgob,nbster), gnudgar(3), osgbdy(nsgob,nbsp)
    lqcflag =       (lverif) .AND. (lqcall)                                    &
              .AND. (mosghd(nsgob,nhqcfw) <= -1)                               &
              .AND. (mosghd(nsgob,nhqcfw) > -99)                               &
              .AND. (ABS( c3600*( osghed(nsgob,nhtime)-acthr)+epsy ) <= cqcwbox)
    lactive =       (osghed(nsgob,nhtime) +wtukare > aiwthr)                   &
              .AND. (osghed(nsgob,nhtime) -wtukara < wtboxe)                   &
              .AND. (osghed(nsgob,nhtime) +wtukare <= nudgend*dtdeh +epsy)     &
              .AND. (mosghd(nsgob,nhpass) == 0)
    IF ((lactive) .AND. (momlhd(nmlob,nhcode) /= nmodes)) THEN
      lactive =       (      (BTEST( mosgbd(nsgob,nbserr),nvru ))              &
                       .AND. (gnudgar(1) > epsy))                              &
                 .OR. (      (BTEST( mosgbd(nsgob,nbserr),nvrt ))              &
                       .AND. (gnudgar(3) > epsy))                              &
                 .OR. (      (BTEST( mosgbd(nsgob,nbserr),nvrq ))              &
                       .AND. (gnudgar(4) > epsy))
    ELSEIF ((lactive) .AND. (momlhd(nmlob,nhcode) == nmodes)) THEN
      lactive =       (      (BTEST( mosgbd(nsgob,nbserr),nvru ))              &
                       .AND. (gnudgms(1) > epsy))                              &
                 .OR. (      (BTEST( mosgbd(nsgob,nbserr),nvrt ))              &
                       .AND. (gnudgms(3) > epsy))                              &
                 .OR. (      (BTEST( mosgbd(nsgob,nbserr),nvrq ))              &
                       .AND. (gnudgms(4) > epsy))
    ENDIF
  ENDIF

  IF ((lactive) .OR. (lqcflag)) THEN
!   IF (.NOT. BTEST( mosghd(nsgob,nhschr) , nvpsbp ) THEN
    IF (lwonl .AND. lfirst) WRITE( nupr,'(''airep taken'',2I5)') nuasta,nsusta
    IF (nuasta < maxuso) THEN
      nuasta = nuasta + 1
      ista   = nuasta + nsusta
      io = mosghd(nsgob,nhio)
      jo = mosghd(nsgob,nhjo)
      IF (.NOT. lactive) io = -io
      isgadm (ista,3) = io
      isgadm (ista,4) = jo
      isgadm (ista,5) = nsgob
      isgadm (ista,6) = 0
      IF (osghed(nsgob,nhtime) > aiwthr+epsy) THEN
        isgadm (ista,1) = 0
        isgadm (ista,2) = nsgob
        tsgadm (ista,1) = c0
        tsgadm (ista,2) = wtukara
        tsgadm (ista,3) = wtukara
      ELSE
        isgadm (ista,1) = nsgob
        isgadm (ista,2) = 0
        tsgadm (ista,1) = wtukare
        tsgadm (ista,2) = c0
        tsgadm (ista,3) = c0
      ENDIF
    ELSE
      nuaex = nuaex + 1
      IF ((lfirst) .OR. (rmod( ntstep*dt,c3600,epsy ) < tconbox)) THEN
        ystid = yosghd(isgadm(ista,5)) (1:ilstidp)
        zpobpr = osghed(nsgob,nhtime) / 100.0_wp
        IF (lwonl)                                                             &
          WRITE( nupr,'("CAUTION: airep ",A ," at",F5.1,",",F6.0,"hPa"         &
                      &," cannot be used: maxuso=",I4," is too small")' )      &
                 ystid, osghed(nsgob,nhtime), zpobpr, maxuso
!       PRINT         '("CAUTION: airep ",A ," at",F5.1,",",F6.0,"hPa"         &
!                     &," cannot be used: maxuso=",I4," is too small")' ,      &
!                  ystid, osghed(nsgob,nhtime), zpobpr, maxuso
      ENDIF
    ENDIF
  ENDIF

ENDDO loop_over_sing_aircraft_reports

 
!-------------------------------------------------------------------------------
!  Section 7: Compute temporal weights
!-------------------------------------------------------------------------------

! multi-level reports
! -------------------

loop_over_multi_level_stations:  DO ista = 1 , nmlsta

  nmlob = MAX( imladm(ista,1) , imladm(ista,2) )
  IF (momlhd(nmlob,nhobtp) == nairep) THEN
    wtuka = wtukara
    wtuke = wtukare
  ELSEIF (      (momlhd(nmlob,nhcode) >= nwp_eu)                               &
          .AND. (momlhd(nmlob,nhcode) <= nravad)) THEN
    wtuka = MIN( 0.5_wp , wtukrsa )
    wtuke = MIN( 0.2_wp , wtukrse )
  ELSE
    wtuka = wtukrsa
    wtuke = wtukrse
  ENDIF
  IF ((imladm(ista,1) > 0) .AND. (imladm(ista,3) > 0)) THEN
    tabdif = aiwthr - omlhed(imladm(ista,1),nhtime)
    wtml (ista,1) = MAX( c0 , c1 - tabdif / tmladm(ista,1) )
    wtml (ista,3) = MAX( c0 , c1 - tabdif / wtuke )
  ELSE
    wtml (ista,1) = c0
    wtml (ista,3) = c0
  ENDIF
  IF ((imladm(ista,2) > 0).AND. (imladm(ista,3) > 0))  THEN
    tabdif = omlhed(imladm(ista,2),nhtime) - aiwthr
    wtml (ista,2) = MAX( c0 , c1 - tabdif / tmladm(ista,2) )
    wtml (ista,4) = MAX( c0 , c1 - tabdif / wtuka )
  ELSE
    wtml (ista,2) = c0
    wtml (ista,4) = c0
  ENDIF
  imladm (ista,3) = ABS( imladm(ista,3) )

ENDDO loop_over_multi_level_stations

! single-level reports
! --------------------

loop_over_single_level_stations:  DO ista = 1 , nsusta + nuasta

  nsgob = MAX( isgadm(ista,1) , isgadm(ista,2) )
  wtuka = wtuksua
  wtuke = wtuksue
  IF (     (mosghd(nsgob,nhobtp) == nairep)                                    &
      .OR. (mosghd(nsgob,nhobtp) == nsatob)) THEN
    wtuka = wtukara
    wtuke = wtukare
  ENDIF
  IF ((isgadm(ista,1) > 0) .AND. (isgadm(ista,3) > 0)) THEN
    tabdif = aiwthr - osghed(isgadm(ista,1),nhtime)
    wtsg (ista,1) = MAX( c0 , c1 - tabdif / tsgadm(ista,1) )
    wtsg (ista,3) = MAX( c0 , c1 - tabdif / wtuke )
  ELSE
    wtsg (ista,1) = c0
    wtsg (ista,3) = c0
  ENDIF
  IF ((isgadm(ista,2) > 0) .AND. (isgadm(ista,3) > 0)) THEN
    tabdif = osghed(isgadm(ista,2),nhtime) - aiwthr
    wtsg (ista,2) = MAX( c0 , c1 - tabdif / tsgadm(ista,2) )
    wtsg (ista,4) = MAX( c0 , c1 - tabdif / wtuka )
  ELSE
    wtsg (ista,2) = c0
    wtsg (ista,4) = c0
  ENDIF
  isgadm (ista,3) = ABS( isgadm(ista,3) )

ENDDO loop_over_single_level_stations

! GPS (IWV) reports
! -----------------

loop_over_gps_stations:  DO ista = 1 , ngpsta

  ngpob = MAX( igpadm(ista,1) , igpadm(ista,2) )
  wtuka = wtuksua
  wtuke = wtuksue
  IF ((igpadm(ista,1) > 0) .AND. (igpadm(ista,3) > 0)) THEN
    tabdif = aiwthr - ogphed(igpadm(ista,1),nhtime)
    wtgp (ista,1) = MAX( c0 , c1 - tabdif / tgpadm(ista,1) )
    wtgp (ista,3) = MAX( c0 , c1 - tabdif / wtuke )
  ELSE
    wtgp (ista,1) = c0
    wtgp (ista,3) = c0
  ENDIF
  IF ((igpadm(ista,2) > 0) .AND. (igpadm(ista,3) > 0)) THEN
    tabdif = ogphed(igpadm(ista,2),nhtime) - aiwthr
    wtgp (ista,2) = MAX( c0 , c1 - tabdif / tgpadm(ista,2) )
    wtgp (ista,4) = MAX( c0 , c1 - tabdif / wtuka )
  ELSE
    wtgp (ista,2) = c0
    wtgp (ista,4) = c0
  ENDIF
  igpadm (ista,3) = ABS( igpadm(ista,3) )

ENDDO loop_over_gps_stations
 
! satellite retrievals
! --------------------

loop_over_sat_retriev_stations:  DO ista = 1 , ntvsta

! ksatid =  motvhd(itvadm(ista,5),nhstid)
! wtuka  =  ssat(kidsat(ksatid))% wtuka
! wtuke  =  ssat(kidsat(ksatid))% wtuke
  wtuka  =  otvhed(itvadm(ista,5),nh1wta)
  wtuke  =  otvhed(itvadm(ista,5),nh1wte)
  IF ((itvadm(ista,1) > 0) .AND. (itvadm(ista,3) > 0)) THEN
    tabdif = aiwthr - otvhed(itvadm(ista,1),nhtime)
    wttv (ista,1) = MAX( c0 , c1 - tabdif / ttvadm(ista,1) )
    wttv (ista,3) = MAX( c0 , c1 - tabdif / wtuke )
  ELSE
    wttv (ista,1) = c0
    wttv (ista,3) = c0
  ENDIF
  IF ((itvadm(ista,2) > 0) .AND. (itvadm(ista,3) > 0)) THEN
    tabdif = otvhed(itvadm(ista,2),nhtime) - aiwthr
    wttv (ista,2) = MAX( c0 , c1 - tabdif / ttvadm(ista,2) )
    wttv (ista,4) = MAX( c0 , c1 - tabdif / wtuka )
  ELSE
    wttv (ista,2) = c0
    wttv (ista,4) = c0
  ENDIF
  itvadm (ista,3) = ABS( itvadm(ista,3) )

ENDDO loop_over_sat_retriev_stations


!-------------------------------------------------------------------------------
!  Section 8: Flag when observation increments need to be computed for
!             writing to VOF or NetCDF feedobs file at current timestep
!-------------------------------------------------------------------------------

IF ((lverif) .AND. (lqcall)) THEN

! multi-level reports
! -------------------

  DO ista = 1 , nmlsta
    nmlob = imladm(ista,1)
    IF (nmlob > 0) THEN
      IF (      (ABS( c3600*( omlhed(nmlob,nhtime)-acthr)+epsy ) <= cqcwbox)   &
          .AND. (momlhd(nmlob,nhqcfw) <= -1)                                   &
          .AND. (momlhd(nmlob,nhqcfw) > -99))                                  &
        momlhd (nmlob,nhqcfw) = 0
    ENDIF
    nmlob = imladm(ista,2)
    IF (nmlob > 0) THEN
      IF (      (ABS( c3600*( omlhed(nmlob,nhtime)-acthr)+epsy ) <= cqcwbox)   &
          .AND. (momlhd(nmlob,nhqcfw) <= -1)                                   &
          .AND. (momlhd(nmlob,nhqcfw) > -99))                                  &
        momlhd (nmlob,nhqcfw) = 0
    ENDIF

  ENDDO

! single-level reports
! --------------------

  DO ista = 1 , nsusta + nuasta
    nsgob = isgadm(ista,1)
    IF (nsgob > 0) THEN
      IF (      (ABS( c3600*( osghed(nsgob,nhtime)-acthr)+epsy ) <= cqcwbox)   &
          .AND. (mosghd(nsgob,nhqcfw) <= -1)                                   &
          .AND. (mosghd(nsgob,nhqcfw) > -99))                                  &
        mosghd (nsgob,nhqcfw) = 0
    ENDIF
    nsgob = isgadm(ista,2)
    IF (nsgob > 0) THEN
      IF (      (ABS( c3600*( osghed(nsgob,nhtime)-acthr)+epsy ) <= cqcwbox)   &
          .AND. (mosghd(nsgob,nhqcfw) <= -1)                                   &
          .AND. (mosghd(nsgob,nhqcfw) > -99))                                  &
        mosghd (nsgob,nhqcfw) = 0
    ENDIF
  ENDDO

! GPS (IWV) reports
! -----------------

  DO ista = 1 , ngpsta
    ngpob = igpadm(ista,1)
    IF (ngpob > 0) THEN
      IF (      (ABS( c3600*( ogphed(ngpob,nhtime)-acthr)+epsy ) <= cqcwbox)   &
          .AND. (mogphd(ngpob,nhqcfw) <= -1)                                   &
          .AND. (mogphd(ngpob,nhqcfw) > -99))                                  &
        mogphd (ngpob,nhqcfw) = 0
    ENDIF
    ngpob = igpadm(ista,2)
    IF (ngpob > 0) THEN
      IF (      (ABS( c3600*( ogphed(ngpob,nhtime)-acthr)+epsy ) <= cqcwbox)   &
          .AND. (mogphd(ngpob,nhqcfw) <= -1)                                   &
          .AND. (mogphd(ngpob,nhqcfw) > -99))                                  &
        mogphd (ngpob,nhqcfw) = 0
    ENDIF
  ENDDO
 
! satellite retrievals
! --------------------

  DO ista = 1 , ntvsta
    ntvob = itvadm(ista,1)
    IF (ntvob > 0) THEN
      IF (      (ABS( c3600*( otvhed(ntvob,nhtime)-acthr)+epsy ) <= cqcwbox)   &
          .AND. (motvhd(ntvob,nhqcfw) <= -1)                                   &
          .AND. (motvhd(ntvob,nhqcfw) > -99))                                  &
        motvhd (ntvob,nhqcfw) = 0
    ENDIF
    ntvob = itvadm(ista,2)
    IF (ntvob > 0) THEN
      IF (      (ABS( c3600*( otvhed(ntvob,nhtime)-acthr)+epsy ) <= cqcwbox)   &
          .AND. (motvhd(ntvob,nhqcfw) <= -1)                                   &
          .AND. (motvhd(ntvob,nhqcfw) > -99))                                  &
        motvhd (ntvob,nhqcfw) = 0
    ENDIF
  ENDDO

ENDIF


!-------------------------------------------------------------------------------
!  Section 9: Printout if required
!-------------------------------------------------------------------------------

  lprsort = (ntstep == -1)
  IF (      ((lfirst) .OR. (rmod( ntstep*dt,c3600,epsy ) < tconbox))           &
      .AND. (aiwthr <= nudgend*dtdeh)) THEN
! for multi-level reports
    IF (lwonl)   WRITE( nupr,'("mladm:",I4," multi-level reports")' ) nmlsta
    IF (lprsort) PRINT       '("mladm:",I4," multi-level reports")' , nmlsta
    DO ista = 1 , nmlsta
      ztobs1 = -999.0_wp
      ztobs2 = -999.0_wp
      IF (imladm(ista,1) > 0) ztobs1 = omlhed(imladm(ista,1),nhtime)
      IF (imladm(ista,2) > 0) ztobs2 = omlhed(imladm(ista,2),nhtime)
      ystid = yomlhd(imladm(ista,5)) (1:ilstidp)
      IF (lwonl)   WRITE( nupr,'("mladm ", 5I5, 2X,A , 2F7.2,9F6.2)' )         &
                          ista, (imladm(ista,i),i=1,4), ystid, ztobs1, ztobs2  &
                              , (tmladm(ista,i),i=1,3), (wtml(ista,i),i=1,4)
      IF (lprsort) PRINT       '("mladm ", 5I5, 2X,A , 2F7.2,9F6.2)' ,         &
                          ista, (imladm(ista,i),i=1,4), ystid, ztobs1, ztobs2  &
                              , (tmladm(ista,i),i=1,3), (wtml(ista,i),i=1,4)
    ENDDO
! for single-level reports
    IF (     ((lfirst) .OR. (rmod( ntstep*dt,21600.0_wp,epsy ) < tconbox)) &
        .OR. (aiwthr >= REAL(MIN( nudgend,nstop ),wp)*dtdeh -1.1_wp)) THEN
      nsgpra = 1
      nsgprd = 1
    ELSE
      nsgpra = MAX( nsusta , i1 )
      nsgprd = MAX( nuasta , i1 )
    ENDIF
    IF (lwonl)   WRITE( nupr,'("sgadm:",I5," surface-level and",I4," upper-"   &
                             &,"air single-level reports")' )  nsusta , nuasta
    IF (lprsort) PRINT       '("sgadm:",I5," surface-level and",I4," upper-"   &
                             &,"air single-level reports")' ,  nsusta , nuasta
    DO ista = nsgpra , nsusta + nuasta , nsgprd
      ztobs1 = -999.0_wp
      ztobs2 = -999.0_wp
      IF (isgadm(ista,1) > 0) ztobs1 = osghed(isgadm(ista,1),nhtime)
      IF (isgadm(ista,2) > 0) ztobs2 = osghed(isgadm(ista,2),nhtime)
      ystid = yosghd(isgadm(ista,5)) (1:ilstidp)
      IF (lwonl)   WRITE( nupr,'("sgadm ", 5I5, 2X,A , 2F7.2,9F6.2)' )         &
                          ista, (isgadm(ista,i),i=1,4), ystid, ztobs1, ztobs2  &
                              , (tsgadm(ista,i),i=1,3), (wtsg(ista,i),i=1,4)
      IF (lprsort) PRINT       '("sgadm ", 5I5, 2X,A , 2F7.2,9F6.2)' ,         &
                          ista, (isgadm(ista,i),i=1,4), ystid, ztobs1, ztobs2  &
                              , (tsgadm(ista,i),i=1,3), (wtsg(ista,i),i=1,4)
    ENDDO
! for GPS reports
    IF (lwonl)   WRITE( nupr,'("gpadm:",I5," GPS reports")' )  ngpsta
    IF (lprsort) PRINT       '("gpadm:",I5," GPS reports")' ,  ngpsta
    DO ista = 1 , ngpsta
      ztobs1 = -999.0_wp
      ztobs2 = -999.0_wp
      IF (igpadm(ista,1) > 0) ztobs1 = ogphed(igpadm(ista,1),nhtime)
      IF (igpadm(ista,2) > 0) ztobs2 = ogphed(igpadm(ista,2),nhtime)
      ystid = yogphd(igpadm(ista,5)) (1:ilstidp)
      IF (lwonl)   WRITE( nupr,'("gpadm ", 5I5, 2X,A , 2F7.2,9F6.2)' )         &
                          ista, (igpadm(ista,i),i=1,4), ystid, ztobs1, ztobs2  &
                              , (tgpadm(ista,i),i=1,3), (wtgp(ista,i),i=1,4)
      IF (lprsort) PRINT       '("gpadm ", 5I5, 2X,A , 2F7.2,9F6.2)' ,         &
                          ista, (igpadm(ista,i),i=1,4), ystid, ztobs1, ztobs2  &
                              , (tgpadm(ista,i),i=1,3), (wtgp(ista,i),i=1,4)
    ENDDO
! for satellite retrievals
    IF (lwonl)   WRITE( nupr,'("tvadm:",I5," satellite retrievals")' )  ntvsta
    IF (lprsort) PRINT       '("tvadm:",I5," satellite retrievals")' ,  ntvsta
    DO ista = 1 , ntvsta
      ztobs1 = -999.0_wp
      ztobs2 = -999.0_wp
      IF (itvadm(ista,1) > 0) ztobs1 = otvhed(itvadm(ista,1),nhtime)
      IF (itvadm(ista,2) > 0) ztobs2 = otvhed(itvadm(ista,2),nhtime)
      ystid = yotvhd(itvadm(ista,5)) (1:ilstidp)
      IF (lwonl)   WRITE( nupr,'("tvadm ", 5I5, 2X,A , 2F7.2,9F6.2)' )         &
                          ista, (itvadm(ista,i),i=1,4), ystid, ztobs1, ztobs2  &
                              , (ttvadm(ista,i),i=1,3), (wttv(ista,i),i=1,4)
      IF (lprsort) PRINT       '("tvadm ", 5I5, 2X,A , 2F7.2,9F6.2)' ,         &
                          ista, (itvadm(ista,i),i=1,4), ystid, ztobs1, ztobs2  &
                              , (ttvadm(ista,i),i=1,3), (wttv(ista,i),i=1,4)
    ENDDO
  ENDIF


  IF (ldump_ascii) THEN
    ! flush YUPRINT file
    IF ((lwonl) .AND. (     (lfirst)                                           &
                       .OR. (rmod( ntstep*dt,c3600,epsy ) < tconbox))) THEN
      CLOSE (nupr)
      OPEN  (nupr,FILE=yuprint,FORM='FORMATTED',STATUS='OLD'                   &
                              ,POSITION='APPEND',IOSTAT=nstat)
      IF (nstat /= 0) PRINT '("OPENING OF FILE yuprint FAIL, local_sort_reports")'
    ENDIF
  ENDIF


!-------------------------------------------------------------------------------
! End of module procedure local_sort_reports
!-------------------------------------------------------------------------------

END SUBROUTINE local_sort_reports


!-------------------------------------------------------------------------------
!+ Module procedure to organize obs increments and QC for surface pressure obs
!-------------------------------------------------------------------------------

SUBROUTINE ps_local_info ( nexceps )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure prepares the lateral spreading of observation
!   increments of 'surface' pressure (i.e. pressure at the lowest model level)
!   by computing all the required 'local' information (i.e. information on the
!   observation, and further parameters from its location) which will later be
!   broadcast to other processors.
!
! Method:
!   Observation increments by call of other procedures. 
!   Computation of temporal weights and other quantities. Threshold quality
!   control.
!
! Written by        :  Christoph Schraff, DWD  (original version: 30.06.97)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER (KIND=iintegers) , INTENT (OUT)    ::  &
    nexceps             ! number of obs increment reports exceeding array size

! Local parameters:
! ----------------

  REAL (KIND = wp)         , PARAMETER  :: &
    c300    =  300.0_wp  ,& !
    c2500   = 2500.0_wp     !

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    io, jo           ,& ! indices of location of observation
    icps             ,& ! index over all observing stations with ps data
    itic             ,& ! index over active obs. at one observing station
    nmlob  , nsgob   ,& ! index of observation in the ODR (obs data record)
    nmlob2 , nsgob2  ,& ! ODR index of 2nd obs. at same obs. station
    nmloba , nsgoba  ,& ! non-zero index  'MAX( n??ob , n??ob2 )'
    itima  , itime   ,& ! range of loop over existing reports at current station
    kk               ,& ! loop index over vertical model levels
    ilvvip  (2)      ,& ! 1: index of lowest obs level with p-z obs
                        ! 2: index of (lower) obs level used for interpolation
    kobtyp           ,& ! observation type
    nlev             ,& ! number of vertical levels in report
    ilqc             ,& ! indicator for current obs. to undergo spatial check
    iqcflg           ,& ! indicator for obs. to have passed latest threshold QC
    nobsta           ,& ! number of active pressure reports from surface sta's
    nstat            ,& ! error flag
    mflgqc           ,& ! QC bit flag (first guess + LBC checks)
    np                  ! number of pressure model fields
                        !   (e.g. 2: original field, boundary field)

  REAL    (KIND=wp   )     ::  &
    col_p    (ke,2)  ,& ! pressure (full value) on main levels    ( Pa  )
    col_z    (ke)    ,& ! geometrical height    of main levels    (  m  )
    zobps    (2)     ,& ! observed pressure
    zpsvob   (2)     ,& ! observed pressure interpolated to lowest model level
!   zpsvob_bc        ,& ! as 'zpsvob', but using lateral boundary fields
!   zzviob           ,& ! scaled height distance betw. obs sta. and level 'ke'
    zpsmod           ,& ! model surface pressure at obs. location
    zpsmod_bc        ,& ! surface pressure of lateral BC field at obs. location
    zpsvim   (2)     ,& ! model pressure interpolated to station height
!   zpsvim_bc        ,& ! pressure of LBC field interpolated to station height
    zpsdz            ,& ! scaled height distance level 'ke' and nearest obs.
    zpsoi    (2)     ,& ! surface pressure obs. increments at current station
    zpsoi_bc (2)     ,& ! surface pressure obs. incr. vs. lateral boundary field
    wt1    , wt2     ,& ! temporal nudging weight at current obs. station
    omyk1  , omyk2   ,& ! local part of nudging weights at current station
    omykf  , omykf2  ,& ! local part of the vertical nudging weight
    zqct     (2)     ,& ! time of current observations
    omyqlt   (2)     ,& ! factor to enhance quality weight dep. on obs. tendency
    r1iflp           ,& ! horizontal correlation scale
    timdif           ,& ! time distance to the observation
    tabdif           ,& ! time distance between the 2 obs at same obs. station
    zyqc     (11)    ,& ! information for QC control messages
!   qctps  , qctps_bc,& ! actual thresholds for quality control (f.g., LBC)
    fqclbcp          ,& ! factor to threshold used for QC checks against LBC
    zqgnudg          ,& ! rel. nudging weight factor for surface station obs.
    zdtird              ! factor for temporal linear interpolation of LBC fields

  LOGICAL                  ::  &
    lqc              ,& ! threshold quality control to be applied
    lqcflag          ,& ! set threshold quality control flags in the ODR
    lqcfrst          ,& ! pressure of report used for the first time
    lveridat         ,& ! fill simulated obs body (for writing to feedback file)
    lwrqc            ,& ! data with obs time close to current time, which are
                        !   rejected by the threshold quality control, are
                        !   printed to a separate file at the current timestep
    lti2             ,& ! linear temporal interpolation as time window
    lpassiv             ! passive data are also used

  CHARACTER (LEN= 1)  yeq  ! for diagnostic print-out

  INTEGER (KIND=iintegers) :: izerror
  CHARACTER (LEN=255)      :: yzerrmsg
  CHARACTER (LEN=25)       :: yzroutine = 'ps_local_info'

! Local (automatic) arrays: None
! -------------------------
!
!------------ End of header ----------------------------------------------------
 
!-------------------------------------------------------------------------------
! Begin Subroutine ps_local_info
!-------------------------------------------------------------------------------

! retreive required microphysics tracer fields
! --------------------------------------------

  ! Retrieve the required microphysics tracers (at nnew)
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnew, ptr = qv)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nnew, ptr = qc)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

!-------------------------------------------------------------------------------
!  Section 1: 'Surface' pressure from multi-level reports
!-------------------------------------------------------------------------------
 
  icps    = 0
  nexceps = 0
  zdtird  = REAL( ntstep + 1 - nlastbound, wp) / REAL( nincbound, wp)
  fqclbcp = qcflbcp
  IF (qcflbcp <= epsy)  fqclbcp = c1

  nobsta = nmlsta

loop_over_sounding_stations:  DO ista = 1, nobsta
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  io = imladm(ista,3)
  jo = imladm(ista,4)
  io_tot = momlhd(imladm(ista,5),nhitot)
  jo_tot = momlhd(imladm(ista,5),nhjtot)

!-------------------------------------------------------------------------------
!  Section 1a: This part is required for nudging and LETKF / feedback files:
!              compute the simulated observations H(x) by applying the forward
!              observation operator, and perform the threshold quality control
!              for individual observations (first guess check)
!-------------------------------------------------------------------------------

  kobtyp = momlhd(imladm(ista,5),nhobtp)
  ystid  = yomlhd(imladm(ista,5)) (1:ilstidp)

  ! no surface pressure obs derived from multi-level aircraft reports
  ! (caution: if one wishes to derive ps from m-l aircraft then beware that
  !           for new multi-level aircraft reports made out of old single-level
  !           reports, omlhed(.,nhtvip), omlhed(.,nhtviz) are undefined and
  !           at the same time it is possible that (momlhd(.,nhqcfw) <= -2) !)
  IF (kobtyp == nairep) THEN
    imladm (ista,7) = 0
                                               CYCLE loop_over_sounding_stations
  ENDIF

  zpsmod     =  p0(io,jo,ke)  +  pp(io,jo,ke,nnew)
  zpsmod_bc  =  p0(io,jo,ke)  +  (c1 - zdtird)* pp_bd(io,jo,ke,nbd1)           &
                              +        zdtird * pp_bd(io,jo,ke,nbd2)
  IF (qcflbcp <= epsy)  zpsmod_bc = zpsmod

  ilqc         =  0
  zpsoi    (1) = c0
  zpsoi    (2) = c0
  zpsoi_bc (1) = c0
  zpsoi_bc (2) = c0

  itima = 2 - MIN( 1 , imladm(ista,1) )
  itime = 1 + MIN( 1 , imladm(ista,2) )

  loop_over_time_index:  DO itim = itima , itime

    !   nmlob > 0
    nmlob  = imladm(ista,itim)

    ! check - if the report enters this routine for the first time and the hence
    !            the QC status is not defined yet and QC needs to be done, and 
    !       - if the threshold quality control flags shall be written to the ODR
    ! --------------------------------------------------------------------------
    !   ( if lqcfrst then (MOD(*nhqcfw,4) < 2) or (*nhqcfw == 0) ;
    !     if (MOD(*nhqcfw,4) < 2) then lqcfrst )
    lqcflag  =                 (momlhd(nmlob,nhqcfw) >=  0)
    lqcfrst  =                 (momlhd(nmlob,nhqcfw) <= -1)                    &
                .AND.          (momlhd(nmlob,nhqcfw) > -99)                    &
                .AND. (MOD( ABS(momlhd(nmlob,nhqcfw)), 4 ) < 2)
    IF (lqcfrst)  momlhd (nmlob,nhqcfw)  =       momlhd(nmlob,nhqcfw) - 2
!   IF (lqcflag)  momlbd (nob,.,nbtqcf)  =  MAX( momlbd(nob,.,nbtqcf), 0 )
    lqc  =       (lqcall) .OR. (lqcfrst) .OR. (lqcflag)
    lveridat  =  (lqcflag) .AND. (lvofoi)

    ! forward observation operator to obtain simulated observations (H(x))
    ! inverse observation operator for nudging
    ! (applied whenever the threshold quality control is applied)
    ! --------------------------------------------------------------------

    IF ((lqc) .OR. (lqcflag)) THEN
      timdif    = omlhed(nmlob,nhtime) - acthr
      tabdif   =  ABS( timdif )

      nlev = momlhd(nmlob,nhnlev)
      col_p (:,1)  =  a_p(io,jo,:)
      np  =  1
      !   if required: 'col_p' from lateral boundary fields (for LBC check)
      IF (qcflbcp > epsy) THEN
        col_p (:,2)  =  p0(io,jo,:)  +  (c1 - zdtird)* pp_bd(io,jo,:,nbd1)     &
                                     +        zdtird * pp_bd(io,jo,:,nbd2)
        np  =  2
      ENDIF
      DO kk = 1 , ke
        col_z (kk)  =  c05 * (hhl(io,jo,kk+1) + hhl(io,jo,kk))
      ENDDO
      lpassiv = .FALSE.

      CALL ps_obs_operator_mult ( ke, np, col_p(:,1:np),col_z, t(io,jo,:,nnew) &
                                ,     qv(io,jo,:), qc(io,jo,:), qrs(io,jo,:)   &
                                ,     r_d, g, rdv, doromx(2)                   &
                                , nlev, omlbdy(nmlob,1:nlev,:)                 &
                                , momlbd(nmlob,1:nlev,:), lpassiv, lveridat    &
                                , smlbdy(:,1:nlev,nmlob)                       &
                                , zpsvob(1:np), zpsdz, zpsvim(1:np), ilvvip )
!     =========================

      IF ((ntstep <= 3) .AND. (lwonl)) THEN
        WRITE( nupr,'("nmlob ",I3,", level with lowest z obs",I3               &
                    &,", height diff ",F5.0,", oro",F6.0,", ps",F8.0)' )       &
               nmlob, ilvvip(1), zpsdz, col_z(ke), col_p(ke,1)
        IF (zpsdz < -epsy)                                                     &
          WRITE( nupr,'("hs-obs < hs-mod : ",A ,",lower level of obs.",I3      &
                      &,",dh-min",F5.0,",ipol. ps-obs",F8.0)' )                &
                 ystid, ilvvip(2), zpsdz, zpsvob(1)
        IF (zpsdz >  epsy)                                                     &
          WRITE( nupr,'("hs-obs > hs-mod : ",A ,",dh",F5.0                     &
                      &,",p-obs orig.",F7.0,",ps: mod/obs",2F8.0)' )           &
                 ystid, zpsdz, omlbdy(nmlob,ilvvip(1),nbtp), col_p(ke,1)       &
                      , zpsvob(1)
      ENDIF

      ! threshold quality control (first guess check)
      ! ---------------------------------------------

      mflgqc = 0
      IF (zpsvob(1) > rmdich) THEN
        lwrqc     = (kobtyp == ntemp)
        zyqc (:)  = c0
        zobps (1) = omlbdy(nmlob,ilvvip(1),nbtp)
        IF (np == 2)  zobps (2) = zobps(1)

        CALL ps_quality_cntl ( np, zobps(1:np), zpsvim(1:np), .FALSE., kobtyp  &
                             , tabdif, qcc(2), fqclbcp, rmdi, lwrqc, mxrbdf    &
                             , momlbd(nmlob,ilvvip(1),1:mxrbdf), mflgqc, zyqc )
!       ====================

        ! fill record for later printing for control of QC
        ! ------------------------------------------------

        IF ((zyqc(1) >= epsy) .AND. (ntotqc < maxqcp)) THEN
          ntotqc = ntotqc + 1
          yyqc (ntotqc   ) = ystid
          myqc (ntotqc, 1) = momlhd(nmlob,nhcode)
          myqc (ntotqc, 2) = NINT( zyqc(1) )
          oyqc (ntotqc, 1) = omlhed(nmlob,nhtime)
          oyqc (ntotqc, 2) = zyqc (2)
          oyqc (ntotqc, 3) = omlhed(nmlob,nhjlat)
          oyqc (ntotqc, 4) = omlhed(nmlob,nhilon)
          oyqc (ntotqc, 5) = zyqc (5)
          oyqc (ntotqc, 6) = zyqc (6)
          oyqc (ntotqc, 7) = zyqc (7)
          lwrqc    =  (lqcall) .AND. (ABS( c3600*timdif+epsy ) <= cqcwbox)
          IF (.NOT. lwrqc)  myqc (ntotqc,1) = 0
        ENDIF
      ENDIF

!-------------------------------------------------------------------------------
!  Section 1b: This is used for nudging only:
!              Get the observation increments and temporal weights,
!              as well as additional local information for later spreading
!-------------------------------------------------------------------------------

      !   fill 'omlhed(:,nhtvip)';  determine which data will be used in the SCC
      omlhed (nmlob,nhtvip)  =  zpsvob(1)
      omlhed (nmlob,nhtviz)  =  zpsdz
      IF (kobtyp == ntemp)  momlhd (nmlob,nhqcps)  =  mflgqc

      IF (omlhed (nmlob,nhtvip) > rmdich) THEN
        !   note: these 'IF' result in a dependency betw. f.g.-SCC and LBC-SCC
        IF (mflgqc >= 1)  ilqc = ilqc + itim
        IF (mflgqc >= 2)  omlhed (nmlob,nhtvip) = - ABS( omlhed(nmlob,nhtvip) )
!       IF (mflgqc <= 1)  omlhed (nmlob,nhtvip) = + ABS( omlhed(nmlob,nhtvip) )
      ENDIF
    ENDIF

    ! set QC flag in 'dmlhed' and in ODR header element 'nhqcfw'
    ! for nudging and for VOF file only (not for NetCDF feedobs file)
    ! ---------------------------------------------------------------

    IF (lqcflag) THEN
      IF (omlhed(nmlob,nhtvip) > rmdich) THEN
        !   active data exist within allowed height difference range
        IF (mflgqc >= 2)  momlhd (nmlob,nhqcfw) = 1
        IF (lvofoi) dmlhed (nso_ps,nmlob) = ABS( omlhed(nmlob,nhtvip) ) - zpsmod
      ELSE
        !   no active data within allowed height difference range
        !  --> use passive obs for QC (but not in spatial consistency check SCC)
        lpassiv = .TRUE.
        ilvvip (1) = 0

        IF (kobtyp == ntemp)                                                   &

          CALL ps_obs_operator_mult ( ke, np, col_p(:,1:np), col_z             &
                                    ,     t(io,jo,:,nnew), qv(io,jo,:)         &
                                    ,     qc(io,jo,:), qrs(io,jo,:)            &
                                    ,     r_d, g, rdv, doromx(2)               &
                                    , nlev, omlbdy(nmlob,1:nlev,:)             &
                                    , momlbd(nmlob,1:nlev,:), lpassiv,lveridat &
                                    , smlbdy(:,1:nlev,nmlob)                   &
                                    , zpsvob(1:np), zpsdz, zpsvim(1:np), ilvvip)
!         =========================

        IF ((zpsvob(1) > rmdich) .AND. (ilvvip(1) > 0)) THEN
        ! lwrqc     = (kobtyp == ntemp)
          zobps (1) = omlbdy(nmlob,ilvvip(1),nbtp)
          IF (np == 2)  zobps (2) = zobps(1)

          CALL ps_quality_cntl (np, zobps(1:np), zpsvim(1:np), .FALSE., kobtyp &
                               ,tabdif, qcc(2), fqclbcp, rmdi, .TRUE., mxrbdf  &
                               ,momlbd(nmlob,ilvvip(1),1:mxrbdf), mflgqc, zyqc )
!         ====================

          ! fill record for later printing for control of QC
          ! ------------------------------------------------

          IF ((zyqc(1) >= epsy) .AND. (ntotqc < maxqcp)) THEN
            ntotqc = ntotqc + 1
            yyqc (ntotqc   ) = ystid
            myqc (ntotqc, 1) = momlhd(nmlob,nhcode)
            myqc (ntotqc, 2) = NINT( zyqc(1) )
            oyqc (ntotqc, 1) = omlhed(nmlob,nhtime)
            oyqc (ntotqc, 2) = zyqc (2)
            oyqc (ntotqc, 3) = omlhed(nmlob,nhjlat)
            oyqc (ntotqc, 4) = omlhed(nmlob,nhilon)
            oyqc (ntotqc, 5) = zyqc (5)
            oyqc (ntotqc, 6) = zyqc (6)
            oyqc (ntotqc, 7) = zyqc (7)
            lwrqc    =  (lqcall) .AND. (ABS( c3600*timdif+epsy ) <= cqcwbox)
            IF (.NOT. lwrqc)  myqc (ntotqc,1) = 0
          ENDIF

          !   passive data being good (no usable active data)
          IF (mflgqc <= 1)  momlhd (nmlob,nhqcfw) = 2
          !   passive data not passing quality control (no usable active data)
          IF (mflgqc >= 2)  momlhd (nmlob,nhqcfw) = 3

          IF (lvofoi) dmlhed (nso_ps,nmlob) = zpsvob(1) - zpsmod
        ELSE
          !   no (active or passive) data within allowed height difference range
          momlhd (nmlob,nhqcfw) = 4
        ENDIF
      ENDIF
    ENDIF

    ! get the observation increments
    ! ------------------------------

    !   bug fixed: compute increments also if obs rejected by individual QC
!   IF (omlhed(nmlob,nhtvip) > c0)) THEN
    IF (omlhed(nmlob,nhtvip) > rmdich) THEN
      zpsoi    (itim)  =  ABS( omlhed(nmlob,nhtvip) )  -  zpsmod
      zpsoi_bc (itim)  =  ABS( omlhed(nmlob,nhtvip) )  -  zpsmod_bc
    ENDIF

  ENDDO loop_over_time_index

  ! get local weights, taking into account the threshold quality control
  ! --------------------------------------------------------------------

  nmloba  = 0
  iqcflg  = 0

  DO itic = 1, 2
    nmlob  = imladm(ista,itic)
    nmlob2 = imladm(ista,3-itic)
    IF (nmlob == 0) THEN
      wt2    = c0
      omykf2 = c0
      zqct (itic) = rmdi
    ELSE
      IF ((gnudg(2) <= epsy) .OR. (omlhed(nmlob,nhtvip) <= c0)                 &
                             .OR. (momlhd(nmlob,nhpass) /=  0)) THEN
        wt2 = c0
      ELSEIF ((nmlob2 == 0) .OR. (omlhed(MAX(nmlob2,i1),nhtvip) <= c0)) THEN
        wt2 = wtml(ista,itic+2)
      ELSE
        wt2 = wtml(ista,itic)
      ENDIF
      nmloba = nmlob
      zqct (itic) = omlhed(nmlob,nhtime)
      ! obs are used in SCC, if active after pre-processing
      !                  and unless 'probably bad' in f.g. check or LBC check
      IF (            (omlhed(nmlob,nhtvip) > c0)                              &
          .OR. (      (omlhed(nmlob,nhtvip) > rmdich)                          &
                .AND. (momlhd(nmlob,nhqcps) <= 2) .AND. (kobtyp == ntemp)))    &
        iqcflg = iqcflg + itic + 4*itic
      omykf2 =  EXP( - (omlhed(nmlob,nhtviz) / doromx(2)) **2 )
    ENDIF
    IF (itic == 1) wt1   = wt2
    IF (itic == 1) omykf = omykf2
  ENDDO
  omyk1 = omykf  * wt1
  omyk2 = omykf2 * wt2
  !   assign 'iqcflg' a negative sign for passive observations
  IF ((nmloba > 0) .AND. (momlhd(MAX(nmloba,i1),nhpass) /= 0)) iqcflg = - iqcflg

  ! determine if linear temporal interpolation possible
  ! ---------------------------------------------------

  lti2 = (imladm(ista,6) == 1) .AND. (omyk1 >= epsy) .AND. (omyk2 >= epsy)
  ! lti2 = (ltipol) .AND. (omyk1 >= epsy) .AND. (omyk2 >= epsy)
  ! IF (lti2) THEN
  !   tabdif = ABS( omlhed(nmlob,nhtime) - omlhed(nmlob2,nhtime) )
  !   lti2   = (ABS(tmladm(ista,2) - tabdif) < epsy) .AND. (tabdif <= tipolmx)
  ! ENDIF

  ! determine the horizontal correlation scale
  ! ------------------------------------------

  r1iflp = MAX( (c1 +(rhtfac(2)-c1) *(c1-MAX(wt1,wt2))) * rhinfl(2) , epsy )

  ! store required local information for further processing
  ! -------------------------------------------------------

! IF ((omyk1 >= epsy) .OR. (omyk2 >= epsy)) THEN
  IF (((zqct(1) > rmdich) .OR. (zqct(2) > rmdich)) .AND. (icps < maxpso)) THEN
    icps  =  icps + 1
    iops     (icps) = io
    jops     (icps) = jo
    iops_tot (icps) = io_tot
    jops_tot (icps) = jo_tot
    zoips  (icps,1) = zpsoi   (1)
    zoips  (icps,2) = zpsoi   (2)
    zoips_b(icps,1) = zpsoi_bc(1)
    zoips_b(icps,2) = zpsoi_bc(2)
    omykps (icps,1) = omyk1
    omykps (icps,2) = omyk2
    r1ifps (icps)   = r1iflp
    zmassps(icps)   = zpsmod - ptpstop
    ztdpps (icps)   = t(io,jo,ke,nnew) / zpsmod
    ltips  (icps)   = lti2
    lmlps  (icps)   = .TRUE.
    ystidps(icps)   = yomlhd(imladm(ista,5)) (1:ilstidg)
    imladm (ista,7) = icps
    iqclps (icps)   = ilqc
    iqcfps (icps)   = iqcflg
    iqcnps (icps)   = ista
    qcimps (icps,1) = zpsmod
    qcimps (icps,2) = zpsmod_bc
    qctmps (icps,1) = zqct(1)
    qctmps (icps,2) = zqct(2)
    qcqfps (icps,1) = omykf
    qcqfps (icps,2) = omykf2
    IF ((omykps(icps,1) > c0) .AND. (ABS(zoips(icps,1)) > 999._wp))        &
      PRINT *,'DPML1 ', ystidps(icps), zoips(icps,1), omykps(icps,1)
    IF ((omykps(icps,2) > c0) .AND. (ABS(zoips(icps,2)) > 999._wp))        &
      PRINT *,'DPML2 ', ystidps(icps), zoips(icps,2), omykps(icps,2)
!   IF ((omykps(icps,1) > c0) .OR. (omykps(icps,2) > c0)) THEN
!     PRINT '("DPM1 ",A,4I4,2F8.0,2F6.2,F6.0,F8.0,F7.4,2L2,I4)', ystidps(icps) &
!            , iops(icps), jops(icps), iops_tot(icps), jops_tot(icps)          &
!            , zoips(icps,1), zoips(icps,2), omykps(icps,1), omykps(icps,2)    &
!            , r1ifps(icps), zmassps(icps), ztdpps(icps)                       &
!            , ltips(icps), lmlps(icps), imladm(ista,7)
!     PRINT '("DPM2 ",A,3I4,F8.0,4F8.4)' , ystidps(icps)                       &
!            , iqclps(icps), iqcfps(icps), iqcnps(icps), qcimps(icps,1)        &
!            , qctmps(icps,1), qctmps(icps,2), qcqfps(icps,1), qcqfps(icps,2)
!   ENDIF
  ELSE
    imladm (ista,7) = 0
    IF (((zqct(1) > rmdich) .OR. (zqct(2) > rmdich)) .AND. (icps >= maxpso))   &
      nexceps = nexceps + 1
  ENDIF

ENDDO loop_over_sounding_stations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
!-------------------------------------------------------------------------------
!  Section 2: 'Surface' pressure from surface reports
!-------------------------------------------------------------------------------

  zqgnudg = c1
  IF (gnudg(2) >= epsy) zqgnudg = gnudgsu(2) / gnudg(2)
 
  nobsta = nsusta

loop_over_surface_stations:  DO ista = 1, nobsta
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  io = isgadm(ista,3)
  jo = isgadm(ista,4)
  io_tot = mosghd(isgadm(ista,5),nhitot)
  jo_tot = mosghd(isgadm(ista,5),nhjtot)

!-------------------------------------------------------------------------------
!  Section 2a: This part is required for nudging and LETKF / feedback files:
!              compute the simulated observations H(x) by applying the forward
!              observation operator, and perform the threshold quality control
!              for individual observations (first guess check)
!-------------------------------------------------------------------------------

  kobtyp = mosghd(isgadm(ista,5),nhobtp)
  ystid  = yosghd(isgadm(ista,5)) (1:ilstidp)

  zpsmod     =  p0(io,jo,ke)  +  pp(io,jo,ke,nnew)
  zpsmod_bc  =  p0(io,jo,ke)  +  (c1 - zdtird)* pp_bd(io,jo,ke,nbd1)           &
                              +        zdtird * pp_bd(io,jo,ke,nbd2)
  IF (qcflbcp <= epsy)  zpsmod_bc = zpsmod
  zpsvim (1)  =  zpsmod
  zpsvim (2)  =  zpsmod_bc

  ilqc   = 0
  zpsoi    (1) = c0
  zpsoi    (2) = c0
  zpsoi_bc (1) = c0
  zpsoi_bc (2) = c0

  itima = 2 - MIN( 1 , isgadm(ista,1) )
  itime = 1 + MIN( 1 , isgadm(ista,2) )

  loop_over_time_index_2:  DO itim = itima , itime

    !   nsgob > 0
    nsgob  = isgadm(ista,itim)

    lqc      =  .FALSE.
    lqcflag  =  .FALSE.
    IF (      (osgbdy(nsgob,nbsz) > rmdich)                                    &
        .AND. (osgbdy(nsgob,nbsp) > rmdich)) THEN

      ! check - if the report enters this routine for the first time and hence
      !            the QC status is not defined yet and QC needs to be done, and
      !       - if the threshold quality control flags shall be written to ODR
      ! ------------------------------------------------------------------------
      !   ( if lqcfrst then (MOD(*nhqcfw,4) < 2) or (*nhqcfw == 0) ;
      !     if (MOD(*nhqcfw,4) < 2) then lqcfrst )
      lqcflag  =        (mosghd(nsgob,nhqcfw) >=  0)
      lqcfrst  =        (mosghd(nsgob,nhqcfw) <= -1)                           &
                  .AND. (mosghd(nsgob,nhqcfw) > -99)                           &
                  .AND. (MOD(ABS(mosghd(nsgob,nhqcfw)),4) < 2)
      IF (lqcfrst)  mosghd (nsgob,nhqcfw) =      mosghd(nsgob,nhqcfw) - 2
!     IF (lqcflag)  mosgbd (nsgob,nbsqcf) = 0
      lqc      =  (lqcall) .OR. (lqcfrst) .OR. (lqcflag)
    ENDIF
    lveridat  =  (lqcflag) .AND. (lvofoi)

    ! forward observation operator to obtain simulated observations (H(x))
    ! inverse observation operator for nudging
    ! (applied whenever the threshold quality control is applied)
    ! --------------------------------------------------------------------

    IF (lqc) THEN
      timdif   =  osghed(nsgob,nhtime) - acthr
      tabdif   =  ABS( timdif )
      lwrqc    =  (lqcall) .AND. (ABS( c3600*timdif+epsy ) <= cqcwbox)

      col_p (:,1)  =  a_p(io,jo,:)
      np  =  1
      !   if required: 'col_p' from lateral boundary fields (for LBC check)
      IF (qcflbcp > epsy) THEN
        col_p (:,2)  =  p0(io,jo,:)  +  (c1 - zdtird)* pp_bd(io,jo,:,nbd1)     &
                                     +        zdtird * pp_bd(io,jo,:,nbd2)
        np  =  2
      ENDIF
      DO kk = 1 , ke
        col_z (kk)  =  c05 * (hhl(io,jo,kk+1) + hhl(io,jo,kk))
      ENDDO

      CALL ps_obs_operator_sing ( ke, np, col_p(:,1:np),col_z, t(io,jo,:,nnew) &
                                ,     qv(io,jo,:), qc(io,jo,:), qrs(io,jo,:)   &
                                ,     r_d, g, rdv                              &
                                , osgbdy(nsgob,:), lveridat , ssgbdy(:,nsgob)  &
                                , zpsvob(1:np), zpsdz, zpsvim(1:np) )
!     =========================

      IF ((ntstep <= 1) .AND. (lwonl)) THEN
        yeq = '='
        IF (zpsdz < -epsy)  yeq = '<'
        IF (zpsdz >  epsy)  yeq = '>'
        WRITE( nupr,'("surf: hs-obs ",A1," hs-mod:",A ,", height diff obs-mod" &
                    &,F6.0,", ps mod/obs",2F8.0)' )                            &
               yeq, ystid, zpsdz, col_p(ke,1), zpsvob(1)
      ENDIF

      ! threshold quality control (first guess check)
      ! ---------------------------------------------

      zyqc (:) = c0
      zobps (1) = osgbdy(nsgob,nbsp)
      IF (np == 2)  zobps (2) = zobps(1)

      CALL ps_quality_cntl ( np, zobps(1:np), zpsvim(1:np), .TRUE., kobtyp     &
                           , tabdif, qccsu(2), fqclbcp, osgbdy(nsgob,nbspst)   &
                           , lwrqc, mxsbdf                                     &
                           , mosgbd(nsgob,1:mxsbdf), mflgqc, zyqc )
!     ====================

      ! fill record for later printing for control of QC
      ! ------------------------------------------------

      IF ((zyqc(1) >= epsy) .AND. (ntotqc < maxqcp)) THEN
        ntotqc = ntotqc + 1
        yyqc (ntotqc   ) = ystid
        myqc (ntotqc, 1) = mosghd(nsgob,nhcode)
        myqc (ntotqc, 2) = NINT( zyqc(1) )
        oyqc (ntotqc, 1) = osghed(nsgob,nhtime)
        oyqc (ntotqc, 2) = zyqc (2)
        oyqc (ntotqc, 3) = osghed(nsgob,nhjlat)
        oyqc (ntotqc, 4) = osghed(nsgob,nhilon)
        oyqc (ntotqc, 5) = zyqc (5)
        oyqc (ntotqc, 6) = zyqc (6)
        oyqc (ntotqc, 7) = zyqc (7)
        IF (.NOT. lwrqc)  myqc (ntotqc,1) = 0
      ENDIF

!-------------------------------------------------------------------------------
!  Section 2b: This is used for nudging only:
!              Get the observation increments and temporal weights,
!              as well as additional local information for later spreading
!-------------------------------------------------------------------------------

      osgbdy (nsgob,nbsvip)  =  zpsvob(1)
!     osgbdy (nsgob,nbsviz)  =  zzviob
      IF (BTEST( mosgbd(nsgob,nbsflg), nvfzbp+nvfbps(4) ))                     &
        osgbdy (nsgob,nbsvip) = - ABS( osgbdy(nsgob,nbsvip) )

      !   note: these 'IF' result in a dependency betw. f.g.-SCC and LBC-SCC
      IF (     (      (.NOT. BTEST( mosgbd(nsgob,nbsflg), nvfzbp+nvfbps(4) ))  &
                .AND. (mflgqc >= 2))                                           &
          .OR. (      (      BTEST( mosgbd(nsgob,nbserr), nvrz ))              &
                .AND. (mflgqc >= 1)))  ilqc = ilqc + itim
      IF (      (.NOT. BTEST( mosgbd(nsgob,nbsflg), nvfzbp+nvfbps(4) ))        &
          .AND. (mflgqc >= 3))                                                 &
        osgbdy (nsgob,nbsvip) = - ABS( osgbdy(nsgob,nbsvip) )
    ENDIF

    ! get the observation increments
    ! ------------------------------

    IF (      (      BTEST( mosgbd(nsgob,nbserr), nvrz ))                      &
        .AND. (.NOT. BTEST( mosgbd(nsgob,nbsflg), nvfzbp+nvfbps(4) ))) THEN
      zpsoi    (itim)  =  ABS( osgbdy(nsgob,nbsvip) )  -  zpsmod
      zpsoi_bc (itim)  =  ABS( osgbdy(nsgob,nbsvip) )  -  zpsmod_bc
    ENDIF

  ENDDO loop_over_time_index_2

  ! get local weights, taking into account the threshold quality control
  ! --------------------------------------------------------------------

  nsgoba  = 0
  iqcflg  = 0

  DO itic = 1, 2
    nsgob  = isgadm(ista,itic)
    nsgob2 = isgadm(ista,3-itic)
    omyqlt (itic) = c1
    IF (nsgob == 0) THEN
      wt2 = c0
      zqct (itic) = rmdi
    ELSE
      IF (     (gnudgsu(2) <= epsy) .OR. (mosghd(nsgob,nhpass) /= 0)           &
          .OR. (.NOT. BTEST( mosgbd(nsgob,nbserr), nvrz   ))                   &
          .OR. (      BTEST( mosgbd(nsgob,nbsqcf), nvrz   ))                   &
          .OR. (      BTEST( mosgbd(nsgob,nbsqcf), nvrzbc ))) THEN
        wt2 = c0
      ELSEIF (     (nsgob2 == 0)                                               &
              .OR. (.NOT. BTEST( mosgbd(MAX(nsgob2,i1),nbserr), nvrz   ))      &
              .OR. (      BTEST( mosgbd(MAX(nsgob2,i1),nbsqcf), nvrz   ))      &
              .OR. (      BTEST( mosgbd(MAX(nsgob2,i1),nbsqcf), nvrzbc ))) THEN
        wt2 = wtsg(ista,itic+2)
      ELSE
        wt2 = wtsg(ista,itic)
      ENDIF
      IF (BTEST( mosgbd(nsgob,nbserr), nvrz ))  nsgoba = nsgob
      zqct (itic) = osghed(nsgob,nhtime)
      ! obs are used in SCC, if active after pre-processing
      !                  and unless 'probably bad' in f.g. check or LBC check
      ! (note: this 'IF' results in a dependency betw. f.g.-SCC and LBC-SCC)
      IF (      (BTEST( mosgbd(nsgob,nbserr), nvrz ))                          &
          .AND. (osgbdy(nsgob,nbsvip) > c0))   iqcflg = iqcflg + itic + 4*itic
      ! (lin) factor to quality weight : 1.0    for ABS(tendency) <=  3 hPa/(3h)
      !                                  qcfpst for ABS(tendency) >= 25 hPa/(3h)
      IF  (osgbdy(nsgob,nbspst) > rmdich)                                      &
        omyqlt (itic) = c1 + (MIN( MAX( ABS( osgbdy(nsgob,nbspst) ) , c300 )   &
                                 , c2500 ) - c300) / 2200.0_wp * (qcfpst - c1)
    ENDIF
    IF (itic == 1) wt1 = wt2
  ENDDO
  omykf = c0
  IF ((nsgoba > 0) .AND. (osgbdy(MAX(nsgoba,i1),nbsviz) > rmdich))             &
    omykf = EXP( - (osgbdy(nsgoba,nbsviz) / doromx(2)) **2 ) * zqgnudg
  omyk1 = omykf * wt1 * omyqlt(1)
  omyk2 = omykf * wt2 * omyqlt(2)
  !   assign 'iqcflg' a negative sign for passive observations
  IF ((nsgoba > 0) .AND. (mosghd(MAX(nsgoba,i1),nhpass) /= 0)) iqcflg = - iqcflg

  ! determine if linear temporal interpolation possible
  ! ---------------------------------------------------

  lti2 = (isgadm(ista,6) == 1) .AND. (omyk1 >= epsy) .AND. (omyk2 >= epsy)
  ! lti2 = (ltipsu) .AND. (omyk1 >= epsy) .AND. (omyk2 >= epsy)
  ! IF (lti2) THEN
  !   tabdif = ABS( osghed(nsgob,nhtime) - osghed(nsgob2,nhtime) )
  !   lti2   = (ABS(tsgadm(ista,2) - tabdif) < epsy) .AND. (tabdif <= tipmxsu)
  ! ENDIF

  ! determine the horizontal correlation scale
  ! ------------------------------------------

  r1iflp = MAX( (c1 +(rhtfsu(2)-c1) *(c1-MAX(wt1,wt2))) * rhiflsu(2) , epsy )

  ! store required local information for further processing
  ! -------------------------------------------------------

! IF ((omyk1 >= epsy) .OR. (omyk2 >= epsy)) THEN
  IF (((zqct(1) > rmdich) .OR. (zqct(2) > rmdich)) .AND. (icps < maxpso)) THEN
    icps  =  icps + 1
    iops     (icps) = io
    jops     (icps) = jo
    iops_tot (icps) = io_tot
    jops_tot (icps) = jo_tot
    zoips  (icps,1) = zpsoi   (1)
    zoips  (icps,2) = zpsoi   (2)
    zoips_b(icps,1) = zpsoi_bc(1)
    zoips_b(icps,2) = zpsoi_bc(2)
    omykps (icps,1) = omyk1
    omykps (icps,2) = omyk2
    r1ifps (icps)   = r1iflp
    zmassps(icps)   = zpsmod - ptpstop
    ztdpps (icps)   = t(io,jo,ke,nnew) / zpsmod
    ltips  (icps)   = lti2
    lmlps  (icps)   = .FALSE.
    ystidps(icps)   = yosghd(isgadm(ista,5)) (1:ilstidg)
    isgadm (ista,7) = icps
    iqclps (icps)   = ilqc
    iqcfps (icps)   = iqcflg
    iqcnps (icps)   = ista
    qcimps (icps,1) = zpsvim(1)
    qcimps (icps,2) = zpsvim(2)
    qctmps (icps,1) = zqct(1)
    qctmps (icps,2) = zqct(2)
    qcqfps (icps,1) = omykf * omyqlt(1)
    qcqfps (icps,2) = omykf * omyqlt(2)
    IF ((omykps(icps,1) > c0) .AND. (ABS(zoips(icps,1)) > 999._wp))        &
      PRINT *,'DPSL1 ', ystidps(icps), zoips(icps,1), omykps(icps,1)
    IF ((omykps(icps,2) > c0) .AND. (ABS(zoips(icps,2)) > 999._wp))        &
      PRINT *,'DPSL2 ', ystidps(icps), zoips(icps,2), omykps(icps,2)
!   IF ((omykps(icps,1) > c0) .OR. (omykps(icps,2) > c0)) THEN
!     PRINT '("DPS1 ",A,4I4,2F8.0,2F6.2,F6.0,F8.0,F7.4,2L2,I4)', ystidps(icps) &
!            , iops(icps), jops(icps), iops_tot(icps), jops_tot(icps)          &
!            , zoips(icps,1), zoips(icps,2), omykps(icps,1), omykps(icps,2)    &
!            , r1ifps(icps), zmassps(icps), ztdpps(icps)                       &
!            , ltips(icps), lmlps(icps), isgadm(ista,7)
!     PRINT '("DPS2 ",A,3I4,F8.0,4F8.4)' , ystidps(icps)                       &
!            , iqclps(icps), iqcfps(icps), iqcnps(icps), qcimps(icps,1)        &
!            , qctmps(icps,1), qctmps(icps,2), qcqfps(icps,1), qcqfps(icps,2)
!   ENDIF
  ELSE
    isgadm (ista,7) = 0
    IF (((zqct(1) > rmdich) .OR. (zqct(2) > rmdich)) .AND. (icps >= maxpso))   &
      nexceps = nexceps + 1
  ENDIF

ENDDO loop_over_surface_stations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  npstot = icps

! flush YUPRINT file
  IF ((lwonl) .AND. (lfirst) .AND. (ldump_ascii)) THEN
    CLOSE (nupr)
    OPEN  (nupr,FILE=yuprint,FORM='FORMATTED',STATUS='OLD'                     &
                            ,POSITION='APPEND',IOSTAT=nstat)
    IF (nstat /= 0) PRINT '("OPENING OF FILE yuprint FAILED in ps_local_info")'
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure ps_local_info
!-------------------------------------------------------------------------------

END SUBROUTINE ps_local_info



!-------------------------------------------------------------------------------
!+ Module procedure for spatial consistency check of surface pressure obs
!-------------------------------------------------------------------------------

SUBROUTINE ps_spatial_check

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure performs an additional quality control step for
!   surface pressure observations by checking their spatial consistency.
!
! Method:
!   At the location of each surface pressure observation not passing the
!   threshold quality control with a reduced threshold (i.e. any observation
!   labelled 'probably good', 'probably bad', or 'bad'), a weighted sum of
!   spreaded observation increments from the neighbouring pressure observations
!   (with label 'good', 'probably good' or 'probably bad') is computed. The
!   difference between the original observation increment of the current
!   observation and this weighted sum is then subjected to the threshold
!   quality control as final quality check.
!   (In this sense, the above-noted weighted sum is used as a bias to the
!    observation to be quality controled.)
!   The threshold for this check is reduced depending on the weights used to
!   compute the bias (the reduction is small if the observations have only
!   small weights at the location of the checked station), and then enhanced
!   by a fraction of the bias.
!   (==> The smallest thresholds will occur in data-dense areas with small
!    observation increments at the neighbouring stations.)
!   Formal remarks:
!   The following quantities are modified in this procedure:
!   - the observation error in the (local) ODR, used for the next times
!     when analysis increments are computed (before doing threshold quality
!     control again)
!   - the quality control flag for the (local) VOF, used for verification only
!   - the local weight 'omykps' with global index, used for the spreading.
!   Note, that these output quantities are not used as input for any computa-
!   tions in this procedure, thus yielding results independent from the domain
!   decomposition. 'omykps' is set to zero for all stations lying outside the
!   local model domain to facilitate the broadcasting of the modifications to
!   'omykps' made on the local stations.
!
! Written by        :  Christoph Schraff, DWD  (original version: 25.06.97)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments: None
! --------------------

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    istaps , istaot  ,& ! loop indices of observing station
    istaqc           ,& ! sorted index of observing station
    i      , j       ,& ! horizontal indices of obs. being checked
    io     , jo      ,& ! horizontal indices of obs. used for checking
    idix   , jdix    ,& ! indices depending on the distance from the obs.
    itic             ,& ! loop index over time
    itica  , itice   ,& ! loop range over time for station being checked
    nsgob  , nmlob   ,& ! indices of reports in the ODR (obs data record)
    nsgob2 , nmlob2  ,& ! indices of other reports at same station
!   nzerr            ,& ! status of memory (de-)allocation
    mflgqc           ,& ! QC bit flag (from spatial consistency check)
    mflgqc_bc        ,& ! QC bit flag (from SCC based on lateral BC field)
    ixnvrz              ! =0: p-obs good ;  =1: p-obs rejected by individual QC

  REAL    (KIND=wp   )     ::  &
      ! suffix '_bc' relates to pressure of lateral boundary field
    zpswi  , zpswi_bc,& ! sum of weighted observation increments
    om2ps  , om2ps_bc,& ! sum of square of 'spreading' weights
    zdist  , zdist2  ,& ! (square of) distance from the checked obs. location
    omysc            ,& ! horizontal nudging weight (for 1 obs)
    zsprcor          ,& ! spreading correction factor to ps- obs. increment
    r1iflp           ,& ! horizontal correlation scale
    zcutop , zcutop2 ,& ! (square of) cut-off radius for area of influence
    zcutmlf          ,& ! correction to 'zcutop2' for ps-obs from TEMPs
    zhcut2           ,& ! square of cut-off radius for area of influence
    zreftim          ,& ! time of obs. being checked
    zpke             ,& ! model pressure at lowest model level
    cfitop           ,& ! weighted conveyor of geopot. correl. to p- correlation
    zqmass           ,& ! mass affected by 'temp. corr.' devided by 'zmassps'
    zqmass3          ,& ! zqmass  **3
    omysc2           ,& ! omysc   **2
    zddsp            ,& ! 'zdist' scaled by horizontal correlation scale
    omyk2 (5)        ,& ! (square of) (local) temporal weight of checking obs.
    tabdif           ,& ! time distance between the checking and checked obs.
    zbias            ,& ! bias used for the spatial consistency quality control
    zbias_bc         ,& ! bias used for the SC check against lateral BC field
    zpsmod           ,& ! model surface pressure at obs. location
    zpsmod_bc        ,& ! surface pressure of lateral BC field at obs. location
    zpsvim           ,& ! model pressure interpolated to station height
    zpsvim_bc        ,& ! pressure of LBC field interpolated to station height
    timdif           ,& ! time distance to checked obs.
    qctps  , qctps_bc,& ! actual thresholds for quality control (f.g., LBC)
    fqclbcp          ,& ! factor to threshold used for QC checks against LBC
    wt1    , wt2        ! temporal nudging weight at current obs. station

  LOGICAL                  ::  &
    local            ,& ! obs. station lies on local domain
    lredowt             ! re-compute temporal weights

! Local (automatic) arrays: None
! -------------------------
!
!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine ps_spatial_check
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1: Selection of observations, for which the 
!             spatial consistency check needs to be done
!-------------------------------------------------------------------------------

  r1iflp  = MAX( (c1 +(rhtfsu(2)-c1) *c05) * rhiflsu(2) , epsy )
  zcutop  = cutofsu(2) * r1iflp
  zcutop2 = zcutop  **2
  zcutmlf = (cutofr(2) / cutofsu(2)) **2

  fqclbcp = qcflbcp
  IF (qcflbcp <= epsy)  fqclbcp = c1

loop_over_all_stations:  DO istaps = 1 , npstot
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  istaqc = isrtpqc(istaps)

! check for the criteria for doing the spatial consistency check:
! - the station lies on the local domain
! - the observation is to be quality controlled at the current timestep
! - the obs. has failed to pass the stringent threshold quality control
! ---------------------------------------------------------------------

  i  =  iops (istaqc)
  j  =  jops (istaqc)
  local =       (j <= jend) .AND. (j >= jstart)                                &
          .AND. (i <= iend) .AND. (i >= istart)

  IF ((local) .AND. (iqclps(istaqc) > 0)) THEN

!   io_tot = iops_tot (istaqc)
!   jo_tot = jops_tot (istaqc)
!   r1iflp = r1ifps   (istaqc)
    ystid  = ystidps  (istaqc) (1:ilstidp)
    lredowt = .FALSE.

    itica = 2 - MOD( iqclps(istaqc) , 2 )
    itice = 1 +      iqclps(istaqc) / 2
 
    loop_over_time_index:  DO itic = itica , itice

      om2ps    = c0
      om2ps_bc = c0
      zpswi    = c0
      zpswi_bc = c0
      zreftim  = qctmps(istaqc,itic)

!-------------------------------------------------------------------------------
!  Section 2: Computation of the bias used for the spatial consistency check,
!             which is set to a weighted sum of spreaded observation increments
!             from the neighbouring pressure data.
!-------------------------------------------------------------------------------

      loop_over_other_stations:  DO istaot = 1 , npstot

        ista   = isrtpqc(istaot)

        io     = iops   (ista)
        jo     = jops   (ista)

! check if obs. 'ista' is within (horizontal) area of influence of obs. 'istaqc'
! (and if obs. 'ista' is not already rejected by the threshold quality control)
        zhcut2 = zcutop2
        IF (lmlps(ista)) zhcut2 = zhcut2 * zcutmlf
        idix   = i - io + ie_tot
        jdix   = j - jo + je_tot
        zdist2 = pxxd2(idix,j) + pyyd(jdix) *pyyd(jdix)                        &
                               - pyyd(jdix) *pxsalpa(idix,j)
        IF (      (zdist2 <= zcutop2) .AND. (ista /= istaqc)                   &
            .AND. (ystidps(ista) (1:ilstidp) /= ystid (1:ilstidp))             &
            .AND. (     (iqcfps(ista) > 0)                                     &
! (use passive reports only for checking passive reports) 
                   .OR. ((iqcfps(ista) < 0) .AND. (iqcfps(istaqc) < 0)))) THEN

! get the horizontal weight, the spreading correction factor to obs. increment
! (designed to avoid orographic footprints), and the obs. quality factor

          zpke    = p0(i,j,ke) + pp(i,j,ke,nnew)
          IF (ntpscor >= 1) THEN
            cfitop  = ztdpps(ista)* zpke / t(i,j,ke,nnew)
            zqmass  = (zpke - ptpstop) / zmassps(ista)
            IF (MOD(ntpscor,2) == 1) THEN
              zqmass3 = zqmass *zqmass *zqmass
              zsprcor =   cfitop * zqmass *zqmass                              &
                        * EXP(  (c1-zqmass3) *0.25_wp / MAX(c1 , zqmass3) )
            ELSE
              zsprcor =   cfitop * zqmass * MAX(c1 , zqmass)                   &
                        * (c05*(c1+zqmass)) **NINT(SIGN( c1,c1-zqmass))
            ENDIF
          ELSE
            zsprcor = zpke / (zmassps(ista) + ptpstop)
          ENDIF
          zdist  =  SQRT( zdist2 )
          zddsp  =  zdist / r1iflp
          omysc  =  (c1 + zddsp) *EXP( -zddsp )
          omysc2 =  omysc * omysc

! get the temporal weights (different from the temporal nudging weights)
! and quality weights

          omyk2 (:)  =  c0
!         IF ((qctmps(ista,1) > rmdich) .AND. (MOD(ABS(iqcfps(ista)),2) == 1)) &
          IF (      (qctmps(ista,1) > rmdich)                                  &
              .AND. (     (BTEST( ABS( iqcfps(ista) ), 0 ))                    &
                     .OR. (BTEST( ABS( iqcfps(ista) ), 2 )))) THEN
            tabdif  =  ABS( zreftim - qctmps(ista,1) )
            omyk2 (5)  =  MAX( c0 , c1 - tabdif / dtchkps )  *  qcqfps(ista,1)
            omyk2 (5)  =  omyk2(5) * omyk2(5)
            IF (BTEST( ABS( iqcfps(ista) ), 0 ))  omyk2 (1)  =  omyk2(5)
            IF (BTEST( ABS( iqcfps(ista) ), 2 ))  omyk2 (3)  =  omyk2(5)
          ENDIF
!         IF ((qctmps(ista,2) > rmdich) .AND. (ABS( iqcfps(ista) ) >= 2)) THEN
          IF (      (qctmps(ista,2) > rmdich)                                  &
              .AND. (     (BTEST( ABS( iqcfps(ista) ), 1 ))                    &
                     .OR. (BTEST( ABS( iqcfps(ista) ), 3 )))) THEN
            tabdif  =  ABS( zreftim - qctmps(ista,2) )
            omyk2 (5)  =  MAX( c0 , c1 - tabdif / dtchkps )  *  qcqfps(ista,2)
            omyk2 (5)  =  omyk2(5) * omyk2(5)
            IF (BTEST( ABS( iqcfps(ista) ), 1 ))  omyk2 (2)  =  omyk2(5)
            IF (BTEST( ABS( iqcfps(ista) ), 3 ))  omyk2 (4)  =  omyk2(5)
          ENDIF

! update the sum of weights and weighted obs. increments

!         om2ps     =  om2ps    + omysc2 *(omyk12 + omyk22)
!         zpswi     =  zpswi    + omysc2 *(  omyk12 *zoips(ista,1)             &
!                                          + omyk22 *zoips(ista,2)) *zsprcor
          om2ps     =  om2ps    + omysc2 *(omyk2(1) + omyk2(2))
          om2ps_bc  =  om2ps_bc + omysc2 *(omyk2(3) + omyk2(4))
          zpswi     =  zpswi    + omysc2 *(  omyk2(1) *zoips  (ista,1)         &
                                           + omyk2(2) *zoips  (ista,2)) *zsprcor
          zpswi_bc  =  zpswi_bc + omysc2 *(  omyk2(3) *zoips_b(ista,1)         &
                                           + omyk2(4) *zoips_b(ista,2)) *zsprcor
        ENDIF

      ENDDO loop_over_other_stations

! compute the bias (i.e weighted obs. increment) used for the spatial
! consistency quality control check
! (without close neighbouring obs., the bias will be small due to 'MAX(.,1)')

      zbias     =  zpswi    / MAX( om2ps    , c1 )
      zbias_bc  =  zpswi_bc / MAX( om2ps_bc , c1 )

!-------------------------------------------------------------------------------
!  Section 3: Setting of flags by the spatial consistency quality control check
!-------------------------------------------------------------------------------

! multi-level observations
! ------------------------

      IF (lmlps(istaqc)) THEN

! compute the adjusted 'model value' against which the obs. are to be checked
        zpsmod     =  qcimps(istaqc,1) + zbias
        zpsmod_bc  =  qcimps(istaqc,2) + zbias_bc

! perform the refined quality control check
        nmlob     =  imladm(iqcnps(istaqc),itic)
        timdif    =  omlhed(nmlob,nhtime) - acthr
        qctps     =  (c1 + qctf(2) *ABS( timdif )) * qcc(2)
        qctps_bc  =  (c1 - (c1 - qcfspat) *MIN( om2ps_bc /c2 , c1 )) * qctps   &
                                                                    * fqclbcp
        qctps     =  (c1 - (c1 - qcfspat) *MIN( om2ps    /c2 , c1 )) * qctps
        qctps_bc  =  qctps_bc  +  qcfbias *ABS( zbias_bc )
        qctps     =  qctps     +  qcfbias *ABS( zbias )
        mflgqc    = NINT( c05 +SIGN( c05, ABS( ABS( omlhed(nmlob,nhtvip) )     &
                                              -zpsmod )    - qctps    ) )
        mflgqc_bc = NINT( c05 +SIGN( c05, ABS( ABS( omlhed(nmlob,nhtvip) )     &
                                              -zpsmod_bc ) - qctps_bc ) )
        ixnvrz    = NINT( c05 - SIGN( c05, omlhed(nmlob,nhtvip) ) )
        IF (ixnvrz /= MAX( mflgqc, mflgqc_bc )) THEN

          IF (ntotqc < maxqcp) THEN
            ntotqc = ntotqc + 1
            yyqc (ntotqc   ) = ystid
            myqc (ntotqc, 1) = momlhd(nmlob,nhcode)
            oyqc (ntotqc, 1) = omlhed(nmlob,nhtime)
            oyqc (ntotqc, 2) = ABS( omlhed(nmlob,nhtvip) ) / 100.0_wp
            oyqc (ntotqc, 3) = omlhed(nmlob,nhjlat)
            oyqc (ntotqc, 4) = omlhed(nmlob,nhilon)
            oyqc (ntotqc, 6) = omlhed(nmlob,nhtvip) / 100.0_wp
            ! keep it without 'ABS': neg. value indicates rejected by indiv. QC
!           oyqc (ntotqc, 6) = ABS( omlhed(nmlob,nhtvip) ) / 100.0_wp
            IF (ixnvrz /= mflgqc) THEN
              myqc (ntotqc, 2) = 15
              oyqc (ntotqc, 5) = qctps              / 100.0_wp
              oyqc (ntotqc, 7) = zpsmod             / 100.0_wp
              oyqc (ntotqc, 8) = zbias              / 100.0_wp
              oyqc (ntotqc, 9) = om2ps
            ELSE
              myqc (ntotqc, 2) = 20
              oyqc (ntotqc, 5) = qctps_bc           / 100.0_wp
              oyqc (ntotqc, 7) = zpsmod_bc          / 100.0_wp
              oyqc (ntotqc, 8) = zbias_bc           / 100.0_wp
              oyqc (ntotqc, 9) = om2ps_bc
            ENDIF
            IF (    (.NOT. ((lqcall).AND.(ABS(c3600*timdif+epsy) <= cqcwbox))) &
               .OR. (omlhed(nmlob,nhtvip) < rmdich))                           &
              myqc(ntotqc,1) = 0
          ENDIF
          omlhed (nmlob,nhtvip) = - omlhed(nmlob,nhtvip)

! if required, adjust the quality control flag: replace 0 by 1 , or  1 by 0
          IF ((momlhd(nmlob,nhqcfw) == 0) .OR. (momlhd(nmlob,nhqcfw) == 1))    &
            momlhd (nmlob,nhqcfw) = 1 - momlhd(nmlob,nhqcfw)

! set flag for need for adjustment of temporal nudging weights
          lredowt = .TRUE.
        ENDIF

! single-level observations
! -------------------------

      ELSE

! compute the adjusted vertically interpolated 'model value' used for the check
        zpsvim     =  qcimps(istaqc,1) + zbias
        zpsvim_bc  =  qcimps(istaqc,2) + zbias_bc

! perform the refined quality control check
        nsgob    =  isgadm(iqcnps(istaqc),itic)
        timdif   =  osghed(nsgob,nhtime) - acthr
        IF ((osgbdy(nsgob,nbspst) > rmdich) .AND. (qccsu(2) > epsy)) THEN
!         qctps  =  qcfctnd *qccsu(2)  +  qcftend *ABS( osgbdy(nsgob,nbspst) )
          qctps     =  qcfctnd *qccsu(2)
          qctps_bc  =  fqclbcp *qctps  +  qcftend *ABS( osgbdy(nsgob,nbspst) )
          qctps     =           qctps  +  qcftend *ABS( osgbdy(nsgob,nbspst) )
        ELSE
          qctps     =  (c1 + qctfsu(2) *ABS( timdif )) * ABS( qccsu(2) )
          IF (mosghd(nsgob,nhobtp) == ndribu)  qctps  =  qctps * qczdrbf
          qctps_bc  =  fqclbcp *qctps
        ENDIF
        qctps     =  (c1 - (c1 - qcfspat) *MIN( om2ps    /c2 , c1 )) * qctps
        qctps_bc  =  (c1 - (c1 - qcfspat) *MIN( om2ps_bc /c2 , c1 )) * qctps_bc
        qctps     =  qctps     +  qcfbias *ABS( zbias    )
        qctps_bc  =  qctps_bc  +  qcfbias *ABS( zbias_bc )
        mflgqc    = NINT( c05 +SIGN(c05, ABS(osgbdy(nsgob,nbsp)-zpsvim) -qctps))
        mflgqc_bc = NINT( c05 +SIGN(c05, ABS(osgbdy(nsgob,nbsp)-zpsvim_bc)     &
                                        -qctps_bc) )
        ! note: if (qcflbcp <= epsy) then flag bit 'nvrzbc' is never set to 1
        IF (qcflbcp <= epsy)  mflgqc_bc = 0
        IF (     (ibit1( mosgbd(nsgob,nbsqcf),nvrz  ) /= mflgqc   )            &
            .OR. (ibit1( mosgbd(nsgob,nbsqcf),nvrzbc) /= mflgqc_bc)) THEN

          !   write message only if the SCC alters the total QC status 
          IF (      (   MAX( ibit1( mosgbd(nsgob,nbsqcf),nvrz  )               &
                           , ibit1( mosgbd(nsgob,nbsqcf),nvrzbc) )             &
                     /= MAX( mflgqc , mflgqc_bc ))                             &
              .AND. (ntotqc < maxqcp)) THEN
            ntotqc = ntotqc + 1
            yyqc (ntotqc   ) = ystid
            myqc (ntotqc, 1) = mosghd(nsgob,nhcode)
            oyqc (ntotqc, 1) = osghed(nsgob,nhtime)
            oyqc (ntotqc, 2) = osgbdy(nsgob,nbsp)   / 100.0_wp
            oyqc (ntotqc, 3) = osghed(nsgob,nhjlat)
            oyqc (ntotqc, 4) = osghed(nsgob,nhilon)
            oyqc (ntotqc, 6) = osgbdy(nsgob,nbsp)   / 100.0_wp
            ! negative value of 'oyqc(.,6)' shall indicate rejected by indiv. QC
            IF (     (BTEST( mosgbd(nsgob,nbsqcf), nvrz   ))                   &
                .OR. (BTEST( mosgbd(nsgob,nbsqcf), nvrzbc )))                  &
              oyqc (ntotqc, 6) = - osgbdy(nsgob,nbsp) / 100.0_wp
            IF (ibit1( mosgbd(nsgob,nbsqcf),nvrz ) /= mflgqc) THEN
              myqc (ntotqc, 2) = 15
              oyqc (ntotqc, 5) = qctps              / 100.0_wp
              oyqc (ntotqc, 7) = zpsvim             / 100.0_wp
              oyqc (ntotqc, 8) = zbias              / 100.0_wp
              oyqc (ntotqc, 9) = om2ps
            ELSE
              myqc (ntotqc, 2) = 20
              oyqc (ntotqc, 5) = qctps_bc           / 100.0_wp
              oyqc (ntotqc, 7) = zpsvim_bc          / 100.0_wp
              oyqc (ntotqc, 8) = zbias_bc           / 100.0_wp
              oyqc (ntotqc, 9) = om2ps_bc
            ENDIF
            IF (.NOT. ((lqcall) .AND. (ABS( c3600*timdif+epsy ) <= cqcwbox)))  &
              myqc(ntotqc,1) = 0
          ENDIF
!         osgbdy (nsgob,nbszer) = - osgbdy(nsgob,nbszer)

! update the quality control flag
!         IF (mosghd(nsgob,nhqcfw) >= 0)                                       &
            CALL MVBITS ( mflgqc   , 0, 1, mosgbd(nsgob,nbsqcf), nvrz   )
          IF (qcflbcp >  epsy)                                                 &
            CALL MVBITS ( mflgqc_bc, 0, 1, mosgbd(nsgob,nbsqcf), nvrzbc )

! set flag for need for adjustment of temporal nudging weights
          lredowt = .TRUE.
        ENDIF
      ENDIF

    ENDDO loop_over_time_index

!-------------------------------------------------------------------------------
!  Section 4: Adjustment of the local (temporal and quality) nudging weights
!-------------------------------------------------------------------------------

    IF (lredowt) THEN
      IF (lmlps(istaqc)) THEN
        DO itic = 1 , 2
          nmlob  = imladm(iqcnps(istaqc),itic)
          nmlob2 = imladm(iqcnps(istaqc),3-itic)
          IF ((nmlob == 0) .OR. (gnudg(2) <= epsy)                             &
                           .OR. (omlhed(MAX(nmlob,i1),nhtvip) <= c0)           &
                           .OR. (momlhd(MAX(nmlob,i1),nhpass) /= 0)) THEN
            wt2 = c0
          ELSEIF ((nmlob2 == 0) .OR. (omlhed(MAX(nmlob2,i1),nhtvip) <= c0)) THEN
            wt2 = wtml(iqcnps(istaqc),itic+2)
          ELSE
            wt2 = wtml(iqcnps(istaqc),itic)
          ENDIF
          IF (itic == 1) wt1 = wt2
        ENDDO
      ELSE
        DO itic = 1 , 2
          nsgob  = isgadm(iqcnps(istaqc),itic)
          nsgob2 = isgadm(iqcnps(istaqc),3-itic)
          IF (     (nsgob == 0) .OR. (gnudgsu(2) <= epsy)                      &
              .OR. (mosghd(MAX(nsgob,i1),nhpass) /= 0)                         &
              .OR. (.NOT. BTEST( mosgbd(MAX(nsgob,i1),nbserr), nvrz   ))       &
              .OR. (      BTEST( mosgbd(MAX(nsgob,i1),nbsqcf), nvrz   ))       &
              .OR. (      BTEST( mosgbd(MAX(nsgob,i1),nbsqcf), nvrzbc ))) THEN
!             .OR. (osgbdy(MAX(nsgob,i1),nbszer) <= c0)) THEN
            wt2 = c0
!         ELSEIF ((nsgob2 == 0) .OR. (osgbdy(MAX(nsgob2,i1),nbszer) <= c0)) THEN
          ELSEIF (    (nsgob2 == 0)                                            &
                 .OR. (.NOT. BTEST( mosgbd(MAX(nsgob2,i1),nbserr),nvrz  ))     &
                 .OR. (      BTEST( mosgbd(MAX(nsgob2,i1),nbsqcf),nvrz  ))     &
                 .OR. (      BTEST( mosgbd(MAX(nsgob2,i1),nbsqcf),nvrzbc))) THEN
            wt2 = wtsg(iqcnps(istaqc),itic+2)
          ELSE
            wt2 = wtsg(iqcnps(istaqc),itic)
          ENDIF
          IF (itic == 1) wt1 = wt2
        ENDDO
      ENDIF
!     PRINT '("space cons.:",A ,2I4,":wt:old1,2/new1,2/omykf",6F6.2)'          &
!            , ystid, iops_tot(istaqc), jops_tot(istaqc)                       &
!            ,(omykps(istaqc,jo),jo=1,2), wt1, wt2, (qcqfps(istaqc,jo),jo=1,2)

      omykps (istaqc,1) = wt1 * qcqfps(istaqc,1)
      omykps (istaqc,2) = wt2 * qcqfps(istaqc,2)
    ENDIF

  ELSEIF (.NOT. local) THEN
    omykps (istaqc,1) = c0
    omykps (istaqc,2) = c0
  ENDIF

ENDDO loop_over_all_stations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~


!-------------------------------------------------------------------------------
! End of module procedure ps_spatial_check
!-------------------------------------------------------------------------------

END SUBROUTINE ps_spatial_check



!-------------------------------------------------------------------------------
!+ Module procedure to organize obs increments and QC for surface-level obs
!-------------------------------------------------------------------------------

SUBROUTINE surf_local_info ( nexcesu )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure prepares the lateral spreading of observation
!   increments from surface-level data by computing all the required 'local'
!   information (i.e. information on the observation, and further parameters
!   from its location) which will later be broadcast to other processors.
!
! Method:
!   Straightforward computation of observation increments, temporal weights,
!   and other quantities. Threshold quality control. Vertical interpolation
!   of surface-level observations to the lowest model level by call of other
!   subroutines.
!
! Written by        :  Christoph Schraff, DWD  (original version: 14.07.97)
!                                        (QC done at obs level,
!                                         instead of model level: 04.02.99)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER (KIND=iintegers) , INTENT (OUT)    ::  &
    nexcesu             ! number of obs increment reports exceeding array size

! Local parameters:
! ----------------

  INTEGER (KIND=iintegers) , PARAMETER  ::  &
    ndqc  =  3          ! max. number of QC control messages per report
                        !   ( = first dimension of 'zyqc')

  REAL    (KIND=wp   )      , PARAMETER  ::  &
    qqdts  =  5.0_wp    ! Gaussian radius of T-2m increments used for
                            ! quality weight for RH-2m increments

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    itic             ,& ! time index of report
    itima  , itime   ,& ! range of loop over existing reports at current station
    icsu             ,& ! index of obs. sta. in data record to be used further
    ivrs             ,& ! index of observed quantity: 1=(u,v); 2=T; 3=RH
    ivar             ,& ! index of observed quantity: 1=(u,v); 3=T; 4=RH
    nsgob  , nsgob2  ,& ! indices of reports in the ODR (obs data record)
    nsgoba           ,& ! non-zero index of report in the ODR (obs data record)
    kobtyp , kcdtyp  ,& ! CMA observation type, observation code type
    kal    , kbl     ,& ! model levels adjacent to equidistant points at obs loc
    nbsx             ,& ! ODR indices for present variable
!   nbsxer           ,& ! ODR index   for mean obs. error of present variable
    nvrx             ,& ! bit pos. for status 'active' for variable
    nbsvix           ,& ! ODR index   for vertically extrapolated obs.
    nlqc             ,& ! number of new QC messages
    ntqc             ,& ! number of total QC messages
    nqc              ,& ! loop index for QC messages 
    nwhere           ,& ! indicator how a temporal nudging weight is computed
    kthr             ,& ! index of uppermost possibly influenced model level
    kthrmin          ,& ! min. 'kthr' for one station
    mv     , kk      ,& ! loop indices
    nstat               ! error flag

  REAL    (KIND=wp   )     ::  &
    wtv    , wtv1    ,& ! temporal nudging weight
    wtvp             ,& ! min. non-zero temporal weight for 1 sta. & variable
    timdif           ,& ! time difference to the observation
    tabdif           ,& ! time distance between the 2 obs at same obs. station
    zobop            ,& ! observation operation for 10 m horizontal wind
    zlop0            ,& ! log of reference pressure p0r
    zflop            ,& ! weight factor for vertical interpol. to equidist. pts.
    xmod             ,& ! model values (vector) at the obs. point
    rhmod_ke         ,& ! model relative humidity at level 'ke' and obs location
!   xvexob           ,& ! surface-level obs. extrapolated to lowest model level
!   yvexob           ,& ! 2nd component of obs. extrapol. to lowest model level
    zthref           ,& ! reference potential temperature (at 2m or ke-level)
    zth2mo           ,& ! model 2-m potential temperature
!   zqvvip           ,& ! model specific humidity extrapolated to obs. level
!   zqvobs           ,& ! observed specific humidity at obs. level
    zcut             ,& ! spreading parameter at vertical cut-off
    zptoz2           ,& ! (square of) factor to convert pressure to height
    zcdist2          ,& ! (square of) enhanced cut-off distance
    zdist2a, zdist2b ,& ! (square of) 'effective cut-off distances'
    zfdist           ,& ! factor used to compute 'effective cut-off distance'
    zupper           ,& ! upper height limits where horiz. radii of infl. used
    zrtinfl             ! actual horizontal correl. scale

  LOGICAL                  ::  &
    lqc              ,& ! threshold quality control to be applied
!   lact             ,& ! current datum is active
!   lti2             ,& ! linear temporal interpolation as time window
    ltakeob          ,& ! use data of present station for further processing
    lqcflag (2)      ,& ! set threshold quality control flags in the ODR
    lqcfrst (2)      ,& ! QC needs to be done since QC status is not yet defined
!   lrejx   (2,3)    ,& ! reject obs(itim,ivrs) by threshold quality control
!   lvarexi          ,& ! current variable exists in current observation report
    lveridat         ,& ! fill simulated obs body (for writing to feedback file)
    lwrqc               ! data with obs time close to current time, which are
                        !   rejected by the threshold quality control, are
                        !   printed to a separate file at the current timestep

  INTEGER (KIND=iintegers) :: izerror
  CHARACTER (LEN=255)      :: yzerrmsg
  CHARACTER (LEN=25)       :: yzroutine = 'surf_local_info'

! Local (automatic) arrays:
! -------------------------

  REAL    (KIND=wp   )     ::  &
    svcutof (3)      ,& ! vertical cut-off radii
    zzol    (ke)     ,& ! height on model levels at the obs. location
    zthol   (ke)     ,& ! potential temperature on model levels at the obs. loc.
!   zhhl    (ke+1)   ,& ! height at half levels
    clc     (ke)     ,& ! (grid-scale + convective) cloud cover
    zqualtq (2)      ,& ! quality factor for RH-2m incr. dep. on T-2m increment
    vip     (2,4)    ,& ! simulated obs = vertically interpolated model values
    zpr     (3)      ,& ! auxilliary quantities for later diagnostic print-out
    zyqc    (ndqc,11)   ! information for QC control messages

!
!------------ End of header ----------------------------------------------------
 
!-------------------------------------------------------------------------------
! Begin Subroutine surf_local_info
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1: Initialisations
!-------------------------------------------------------------------------------

! retreive required microphysics tracer fields
! --------------------------------------------

  ! Retrieve the required microphysics tracers (at nnew)
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnew, ptr = qv)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nnew, ptr = qc)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

! initialisations related to observations
! ---------------------------------------

  !   vertical distance between observation point and edge of area of influence
  IF (msprpsu <= 1) THEN
    svcutof (1) = SQRT( vcorlsu(1) * vcutosu(1) )
    svcutof (2) = SQRT( vcorlsu(3) * vcutosu(3) )
    svcutof (3) = SQRT( vcorlsu(4) * vcutosu(4) )
  ELSEIF (msprpsu == 2) THEN
    svcutof (1) = SQRT( vcorlsu(1) * vcutosu(1) / (c1 + vcorlsu(1)             &
                                                       /(vpblsu(1)*vpblsu(1))) )
    svcutof (2) = SQRT( vcorlsu(3) * vcutosu(3) / (c1 + vcorlsu(3)             &
                                                       /(vpblsu(3)*vpblsu(3))) )
    svcutof (3) = SQRT( vcorlsu(4) * vcutosu(4) / (c1 + vcorlsu(4)             &
                                                       /(vpblsu(4)*vpblsu(4))) )
  ENDIF

  DO ivrs   = 1 , 3
    DO ista = 1 , maxsgo
      zwtsu   (ista,1,ivrs) = c0
      zwtsu   (ista,2,ivrs) = c0
      xoisu   (ista,1,ivrs) = c0
      xoisu   (ista,2,ivrs) = c0
      zvcutsu (ista,  ivrs) = c0
      zriflsu (ista,  ivrs) = c0
      zqualsu (ista,1,ivrs) = c0
      zqualsu (ista,2,ivrs) = c0
      kviflsu (ista,  ivrs) = ke + 1
      ltisu   (ista,  ivrs) = .FALSE.
    ENDDO
  ENDDO
  DO ista = 1 , maxsgo
    xoisu   (ista,1,4) = c0
    xoisu   (ista,2,4) = c0
    zpobsu  (ista    ) = c1
    zzobsu  (ista    ) = c0
    zsposu  (ista    ) = c0
    zrtdgsu (ista    ) = c0
    kobtysu (ista    ) = 0
    kcdtysu (ista    ) = 0
  ENDDO

  zlop0  = LOG( p0r )

  ! initialise loop over surface-level 'stations'
  ! ---------------------------------------------

  icsu    = 0
  ktopsu  = ke
  nexcesu = 0
 
loop_over_surface_stations:  DO ista = 1 , nsusta
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  icsu    = icsu + 1
  kthrmin = ke
  ltakeob = .FALSE.

  io = isgadm(ista,3)
  jo = isgadm(ista,4)
  nsgoba = MAX( isgadm(ista,1) , isgadm(ista,2) )
  io_tot = mosghd(nsgoba,nhitot)
  jo_tot = mosghd(nsgoba,nhjtot)

!-------------------------------------------------------------------------------
!  Section 2: This part is required for nudging and LETKF / feedback files:
!             compute the simulated observations H(x) by applying the forward
!             observation operator, and perform the threshold quality control
!             for individual observations (first guess check)
!-------------------------------------------------------------------------------

  kobtyp = mosghd(nsgoba,nhobtp)
  kcdtyp = mosghd(nsgoba,nhcode)
  ystid  = yosghd(isgadm(ista,5)) (1:ilstidg)

  itima = 2 - MIN( 1 , isgadm(ista,1) )
  itime = 1 + MIN( 1 , isgadm(ista,2) )
    
  loop_over_time_index_1:  DO itim = itima , itime

    !   nsgob > 0
    nsgob  =  isgadm(ista,itim)

    ! check - if the report enters this routine for the first time and the hence
    !            the QC status is not defined yet and QC needs to be done, and
    !       - if the threshold quality control flags shall be written to the ODR
    ! --------------------------------------------------------------------------
    !   ( if lqcfrst then (MOD(*nhqcfw,8) < 4) or (*nhqcfw == 0) ;
    !     if (MOD(*nhqcfw,8) < 4) then lqcfrst )
    lqcflag (itim)  =                 (mosghd(nsgob,nhqcfw) >=    0)
    lqcfrst (itim)  =                 (mosghd(nsgob,nhqcfw) <=   -1)           &
                       .AND.          (mosghd(nsgob,nhqcfw) >   -99)           &
                       .AND. (MOD( ABS(mosghd(nsgob,nhqcfw)), 8 ) < 4)
    IF (lqcfrst(itim))  mosghd (nsgob,nhqcfw)  =  mosghd(nsgob,nhqcfw) - 4

    lveridat =  lqcflag(itim) .AND. (lvofoi)
    timdif   =  osghed(nsgob,nhtime) - acthr
    tabdif   =  ABS( timdif )
    lqc      =  (lqcall) .OR. (lqcfrst(itim)) .OR. (lqcflag(itim))
    lwrqc    =  (lqcall) .AND. (ABS( c3600*timdif+epsy ) <= cqcwbox)

    ! forward observation operator to obtain simulated observations (H(x))
    ! --------------------------------------------------------------------

    IF (itype_tran /= 1) THEN

      CALL surf_obs_operator ( ke, u_10m(io,jo), v_10m(io,jo), t_2m(io,jo)     &
                             ,     td_2m(io,jo), t(io,jo,:,nnew), a_z(io,jo,:) &
                             ,     ps(io,jo,nnew), osghed(nsgob,nhsurf)        &
                             , osgbdy(nsgob,:), osghed(nsgob,nhalt)            &
                             , kml700, kml850, lveridat, rdv, b1, b2w, b3, b4w &
                             , ssgbdy(:,nsgob) , vip(itim,:) , loiqv2m , zpr )
    ! ======================

      IF ((zpr(1) > rmdich) .AND. (lwonl) .AND. (MOD( ntstep,90 )) <= 1)       &
        WRITE( nupr,'(A ,'' q2m-qvip'',I2,3F6.3,4F7.2,2F8.5)' )                &
               ystid, itim, osgbdy(nsgob,nbsrh), zpr(1), vip(itim,4)           &
             , osgbdy(nsgob,nbst), t_2m(io,jo), t(io,jo,ke,nnew)               &
             , vip(itim,3), zpr(2), zpr(3)
        !  (zpr(1) =   fpvsw ( td_2m, b1, b2w, b3, b4w )                       &
        !            / fpvsw ( t_2m , b1, b2w, b3, b4w )                       &
        !   zpr(2) = zqvvip   ;   zpr(3) = zqvobs)

    ELSEIF (itype_tran == 1) THEN
      IF (itim == itima)                                                       &
        rhmod_ke  =  fq2pv ( qv(io,jo,ke) + qc(io,jo,ke), a_p(io,jo,ke), rdv ) &
                   / fpvsw ( t(io,jo,ke,nnew), b1, b2w, b3, b4w )

      ! recompute the observed value from the surface-level obs
      ! (by the applying inverse obs operator) when QC is (re-)applied
      ! --------------------------------------------------------------
      IF (lqc) THEN
        IF (osgbdy(nsgob,nbsu) > rmdich) THEN

          CALL surf_obs ( osgbdy(nsgob,nbsu), osgbdy(nsgob,nbsv), 1            &
                        , osgbdy(nsgob,nbsviu), osgbdy(nsgob,nbsviv), zobop )
!         =============
          vip (itim,1)  =  a_u(io,jo,ke) * zobop
          vip (itim,2)  =  a_v(io,jo,ke) * zobop
        ENDIF
        IF (osgbdy(nsgob,nbst) > rmdich)                                       &

          CALL surf_obs ( osgbdy(nsgob,nbst), t(io,jo,ke,nnew), 2              &
                        , osgbdy(nsgob,nbsvit), vip(itim,3) )
!         =============

        IF (osgbdy(nsgob,nbsrh) > rmdich)                                      &

          CALL surf_obs ( osgbdy(nsgob,nbsrh), rhmod_ke, 3                     &
                        , osgbdy(nsgob,nbsviq), vip(itim,4) )
!         =============

        IF (lveridat) THEN
          IF (osgbdy(nsgob,nbsu ) > rmdich)  ssgbdy (nso_u ,nsgob) = vip(itim,1)
          IF (osgbdy(nsgob,nbsv ) > rmdich)  ssgbdy (nso_v ,nsgob) = vip(itim,2)
          IF (osgbdy(nsgob,nbst ) > rmdich)  ssgbdy (nso_t ,nsgob) = vip(itim,3)
          IF (osgbdy(nsgob,nbsrh) > rmdich)  ssgbdy (nso_rh,nsgob) = vip(itim,4)
        ENDIF
      ENDIF
    ENDIF

    !   model cloud cover written to SOR for feedobs / feedback files
    IF ((lrad) .AND. (lveridat)) THEN
!     IF (osgbdy(nsgob,nbscl) > rmdich)                                        &
!       ssgbdy (nso_cl ,nsgob)  =  clcl(io,jo)
!!    IF (IBITS(mosgbd(nsgob,nbscwg),nvnbp,nvnoc) <= 9)                        &
!     IF (osgbdy(nsgob,nbsct) > rmdich)                                        &
!       ssgbdy (nso_ct ,nsgob)  =  clct(io,jo)

!     CALL hhl_col ( ke, a_z(io,jo,:), osghed(nsgob,nhsurf), zhhl )
!     ============         

      clc (:)  =  clc_sgs(io,jo,:) + clc_con(io,jo,:) *(c1 - clc_sgs(io,jo,:))

      CALL cloud_obs_operator ( ke, clc, clct(io,jo), clcl(io,jo), clcm(io,jo) &
                                       , clch(io,jo), hhl(io,jo,:)             &
                              , osgbdy(nsgob,:) , ssgbdy(:,nsgob) )
!     =======================

    ENDIF

    IF (lqc) THEN
    
      ! threshold quality control (first guess check)   ('ps', 'qcvf' can have
      ! ---------------------------------------------    arbitrary values here)
    
      CALL sing_quality_cntl ( osgbdy(nsgob,:), vip(itim,:), ps(io,jo,nnew)    &
                             , -1, kobtyp, tabdif, qcvf, qccsu, ndqc, lwrqc    &
                             , mosgbd(nsgob,:), nlqc, zyqc )
    ! ======================
    
      ! Fill record for later printing for control of QC 
      ! ------------------------------------------------
      !   adjust 'nlqc' if space in arrays 'oyqc' is insufficient
      nlqc  =  MIN( maxqcp - ntotqc , nlqc )
      DO nqc = 1 , nlqc
        ntqc  =  ntotqc + nqc
        yyqc (ntqc   ) = ystid
        myqc (ntqc, 1) = mosghd(nsgob,nhcode) 
        myqc (ntqc, 2) = NINT( zyqc(nqc,1) ) 
        IF (.NOT. lwrqc) myqc (ntqc,1) = 0
        oyqc (ntqc, 1) = osghed(nsgob,nhtime)
        oyqc (ntqc, 2) = zyqc (nqc, 2)
        oyqc (ntqc, 3) = osghed(nsgob,nhjlat)
        oyqc (ntqc, 4) = osghed(nsgob,nhilon)
        oyqc (ntqc, 5) = zyqc (nqc, 5)
        oyqc (ntqc, 6) = zyqc (nqc, 6)
        oyqc (ntqc, 7) = zyqc (nqc, 7)
        oyqc (ntqc, 8) = zyqc (nqc, 8)
        oyqc (ntqc, 9) = zyqc (nqc, 9)
        oyqc (ntqc,10) = zyqc (nqc,10)
!       oyqc (ntqc,11) = zyqc (nqc,11)
        IF ((MOD( myqc(ntqc,2), 4 ) == 1) .AND. (itype_tran == 1))             &
          oyqc (ntqc,11) = c1 / zobop
      ENDDO
      ntotqc  =  ntotqc  +  nlqc
    ENDIF

  ENDDO loop_over_time_index_1


!-------------------------------------------------------------------------------
!  Section 3: This is used for nudging only:
!             Get the observation increments and temporal weights,
!             as well as additional local information for later spreading
!-------------------------------------------------------------------------------

  ! get:
  !   iosu,josu: location of the obs. station
  !   ystidsu  : station id.
  !   zsposu   : values of spreading parameter at the obs. point
  !   zzobsu   : height   at the obs. point
  !   zpobsu   : pressure at the obs. point
  !   zrtdgsu  : (Rd / g) * Tv (conversion factor ln(p) --> z) at the obs. point
  !   zsprsu   : values of spreading parameter at obs. location on model levels
  !   znissq   : vertical profile of the parameter which models the non-isotropy
  !   --------------------------------------------------------------------------

  !   get 'z' at model level of the obs, and profile of 'z' at model levels
  DO kk = 1 , ke
    zzol        (kk) = a_z(io,jo,kk)
    zsprsu (icsu,kk) = zzol(kk)
  ENDDO
  zsposu  (icsu) = zzol(ke)
  zzobsu  (icsu) = zsposu(icsu)
  zpobsu  (icsu) = a_p(io,jo,ke)
  zrtdgsu (icsu) = r_d/g* t(io,jo,ke,nnew)*( c1 + rvd_m_o* qv(io,jo,ke)        &
                                            -qc(io,jo,ke) -qrs(io,jo,ke))
  iosu    (icsu) = io
  josu    (icsu) = jo
  iosu_tot(icsu) = io_tot
  josu_tot(icsu) = jo_tot
  kobtysu (icsu) = kobtyp
  kcdtysu (icsu) = kcdtyp
  ystidsu (icsu) = yosghd(isgadm(ista,5)) (1:ilstidg)
  ystid          = ystidsu(icsu) (1:ilstidp)

  !   if necessary: get profile of 'theta' at model levels
  DO kk = 1 , ke
    zthol (kk) = t(io,jo,kk,nnew) * EXP( rdocp *(zlop0 - LOG(a_p(io,jo,kk))) )
  ENDDO
  IF ((lnissu) .OR. (msprpsu == 2)) THEN
    DO kk = 1 , ke
      IF                 (msprpsu == 2)  zsprsu (icsu,kk) = zthol(kk)
      IF ((lnissu) .AND. (msprpsu <= 1)) znissu (icsu,kk) = zthol(kk)
      IF ((lnissu) .AND. (msprpsu == 2)) znissu (icsu,kk) = zzol (kk)
    ENDDO
  ENDIF
  !   get spreading parameter at obs. level
  IF (msprpsu == 2) zsposu (icsu) = zthol(ke)

  !   if necessary: profiles of 'theta' + 'z' at 'exp(-z/H) -equidistant' levels
  IF ((lnissu) .AND. (msprpsu >= 1)) THEN
    DO mv = 1 , mxispr
      kk = ke
      DO WHILE (  (   (msprpsu <  2) .AND. (zzol (MAX(kk,i1)) <= zsdni(mv,1))  &
                 .OR. (msprpsu == 2) .AND. (zthol(MAX(kk,i1)) <= zsdni(mv,2))) &
            .AND. (kk > 0))
        kk = kk - 1
      ENDDO
      kbl = MIN( INT( kk+1 ,iintegers) , ke )
      kal = MAX(      kk               , i1 )
      zflop = c1
      IF (msprpsu < 2) THEN
        IF (kal /= kbl) zflop =   (zzol(kbl) - zsdni(mv,1))                    &
                                / (zzol(kbl) - zzol(kal))
        znissq (icsu,mv) = (c1-zflop) *zthol(kbl)  + zflop *zthol(kal)
      ELSE
        IF (kal /= kbl) zflop =   (zthol(kbl) - zsdni(mv,2))                   &
                                / (zthol(kbl) - zthol(kal))
        znissq (icsu,mv) = (c1-zflop) *zzol(kbl) + zflop *zzol(kal)
        IF (kal == ke)                                                         &
          znissq (icsu,mv) = zzol(ke) - (zthol(ke) -zsdni(mv,2)) /thdzt
      ENDIF
    ENDDO
  ENDIF

  !   get potential temperature difference to obs. point
  zth2mo = t_2m(io,jo) * EXP( rdocp *(zlop0 - LOG(ps(io,jo,nnew))) )
  zthref = MIN( zth2mo , zthol(ke) )
  !   (T_2m does not exist at first timestep (set to 0.) if (itype_tran == 1) !)
  IF (zth2mo < c1)  zthref = zthol(ke)
  DO kk = 1 , ke
    zpblsu (icsu,kk) = MAX( zthol(kk) - zthref , c0 )
  ENDDO
      
  ! initialise loop over variables
  ! ------------------------------
    
  loop_over_variables:  DO ivrs = 1, 3

    IF (ivrs == 1) THEN
      ivar   = 1
      nbsx   = nbsu
      nvrx   = nvru
      nbsvix = nbsviu
      xmod   = a_u(io,jo,ke)
    ELSEIF (ivrs == 2) THEN
      ivar   = 3
      nbsx   = nbst
      nvrx   = nvrt
      nbsvix = nbsvit
      xmod   = t(io,jo,ke,nnew)
    ELSEIF (ivrs == 3) THEN
      ivar   = 4
      nbsx   = nbsrh
      nvrx   = nvrq
      nbsvix = nbsviq
      IF (      (itype_tran == 1) .AND. (gnudgsu(ivar) > epsy)                 &
          .AND. (itima <= itime))   xmod  = rhmod_ke
    ENDIF

    ! get the observation increments
    ! ------------------------------

    IF (gnudgsu(ivar) > epsy) THEN
      DO itim = itima , itime
        nsgob  = isgadm(ista,itim)
        IF (      (      BTEST( mosgbd(nsgob,nbserr), nvrx ))                  &
            .AND. (.NOT. BTEST( mosgbd(nsgob,nbsqcf), nvrx ))) THEN
          IF (itype_tran == 1) THEN
            xoisu (icsu,itim,ivar)  =  osgbdy(nsgob,nbsvix) - xmod
            IF (ivar == 1)                                                     &
              xoisu (icsu,itim,2)   =  osgbdy(nsgob,nbsviv) - a_v(io,jo,ke)
          ELSE
            xoisu (icsu,itim,ivar)  =  osgbdy(nsgob,nbsx) - vip(itim,ivar)
            IF (ivar == 1)                                                     &
              xoisu (icsu,itim,2)   =  osgbdy(nsgob,nbsv) - vip(itim,2)
          ENDIF
        ENDIF
      ENDDO
    ENDIF

    ! get temporal weights, taking into account the threshold quality control
    ! -----------------------------------------------------------------------

    DO itic = 1, 2
      nsgob  = isgadm(ista,itic)
      nsgob2 = isgadm(ista,3-itic)
!     IF ((nsgob == 0) .OR. (gnudgsu(ivar) <= epsy)                            &
!                      .OR. (osgbdy(MAX(nsgob,i1),nbsxer) <= rmdich)           &
!                      .OR. (mosghd(MAX(nsgob,i1),nhpass) /= 0)) THEN
      IF (     (nsgob == 0) .OR. (gnudgsu(ivar) <= epsy)                       &
          .OR. (mosghd(MAX(nsgob,i1),nhpass) /= 0)                             &
          .OR. (.NOT. BTEST( mosgbd(MAX(nsgob,i1),nbserr), nvrx ))) THEN
        wtv = c0
      ELSE
        nwhere = 2
        wtv = wtsg(ista,itic)
!       IF (osgbdy(nsgob,nbsxer) <= c0) THEN
        IF (BTEST( mosgbd(nsgob,nbsqcf), nvrx )) THEN
          wtv = c0
          nwhere = 3
        ELSEIF (     (nsgob2 == 0)                                             &
                .OR. (.NOT. BTEST( mosgbd(MAX(nsgob2,i1),nbserr), nvrx ))      &
                .OR. (      BTEST( mosgbd(MAX(nsgob2,i1),nbsqcf), nvrx ))) THEN
          wtv = wtsg(ista,itic+2)
          nwhere = 4
        ENDIF
        IF (wtv < epsy)  wtv = c0
        IF ((ntstep <= 1) .AND. (lwonl)) THEN
          IF ((ivar == 1) .AND. (     ((nwhere == 2) .AND. (itic == 1))        &
                                 .OR. (nwhere >= 4)))                          &
            WRITE( nupr,'("sfc: wtu", 3I4,F7.4,2F7.2)' )                       &
                   io_tot, jo_tot, nwhere, wtv, vip(itic,1), osgbdy(nsgob,nbsx)
        ENDIF
      ENDIF
      IF (itic == 1)    wtv1 = wtv
      zwtsu   (icsu,itic,ivrs) = wtv
      IF (wtv >= epsy) THEN
        zqualsu (icsu,itic,ivrs) = c1
!       IF (ivrs <= 2)  zqualsu (icsu,itic,ivrs) = c1
!       IF (ivrs == 3)  zqualsu (icsu,itic,ivrs) = zqualtq(itic)
        IF ((ivrs == 3) .AND. (lqfqv2m) .AND.(osgbdy(nsgob,nbst) > rmdich)) THEN
          !   quality factor for RH-2m depending on temperature obs increment
          zqualtq (itic)  =  (osgbdy(nsgob,nbst) - vip(itic,3)) / qqdts
          zqualsu (icsu,itic,ivrs)  =  exp( - zqualtq(itic) *zqualtq(itic) )
        ENDIF
      ENDIF
    ENDDO

!-------------------------------------------------------------------------------
!  Section 4: This is used for nudging only:
!             Get the vertical upper cut-off for a report,
!             an upper estimate for the range of model levels influenced by a
!             report, and upper estimates for horizontal correlation scales
!             for target grid pts. at the vertical cut-off
!-------------------------------------------------------------------------------

    ! get uppermost model levels within the vertical cut-off at the obs location
    ! --------------------------------------------------------------------------

    IF ((zwtsu(icsu,1,ivrs) >= epsy) .OR. (zwtsu(icsu,2,ivrs) >= epsy)) THEN
      wtvp = c1
      IF (wtv1 >= epsy) wtvp = MIN( wtvp , wtv1 )
      IF (wtv  >= epsy) wtvp = MIN( wtvp , wtv  )
      IF (msprpsu == 2) zcut = zsposu(icsu) + svcutof(ivrs)
!     IF (msprpsu <= 1) zcut = zsposu(icsu) + svcutof(ivrs) *zrtdgsu(icsu)
      IF (msprpsu <= 1) THEN
        !   normal cut-off height (without account for pot. temperature differ.)
        zcut     =  zsposu(icsu) + svcutof(ivrs) *zrtdgsu(icsu)
        !   adapt cut-off height due to PBL (pot. temp. diff) weight correction:
        !   --> obtain (square of) cut-off distance (small cut-off parameter
        !                                            values are enhanced)
        zptoz2   =  zrtdgsu(icsu) *zrtdgsu(icsu)
        zcdist2  =  vcorlsu(ivar) *MAX( vcutosu(ivar), c2*c2 ) *zptoz2
        !   compute 'effective distances' which are defined as normal distances
        !            -------------------  enhanced by scaled potential
        !                                 temperature differences
        !   and set the new cut-off height equal to the height where the
        !   'effective distance' is equal to the cut-off distance 'zcdist2'
        !     (note: zsposu = zzol(ke))
        zfdist   =  vcorlsu(ivar) *zptoz2 /(vpblsu(ivar)*vpblsu(ivar))
        zdist2a  =  zfdist *zpblsu(icsu,ke) *zpblsu(icsu,ke)
        kthr = ke
        DO WHILE ((kthr >  1) .AND. (zdist2a < zcdist2))
          kthr    = kthr - 1
          zdist2b = zdist2a
          zdist2a =  (zzol(kthr)-zsposu(icsu)) *(zzol(kthr)-zsposu(icsu))      &
                   + zfdist *zpblsu(icsu,kthr) *zpblsu(icsu,kthr)
        ENDDO
        zfdist = zcut
        IF (kthr == ke) THEN
          zcut = zzol(ke) + epsy
        ELSEIF (zdist2a >= zcdist2) THEN
          zcut = MIN( zcut , zzol(kthr+1) + (zcdist2 -zdist2b)                 &
                                           /(zdist2a -zdist2b)*( zzol(kthr)    &
                                                                -zzol(kthr+1)) )
        ENDIF
        IF ((lwonl) .AND. (MOD( ntstep,90 ) <= 1) .AND. (ivrs == 3)) THEN
          kbl = MIN( kthr+1 , ke )
          WRITE( nupr,'(A ,'' q2m-zcut'',3F6.0,2F6.1,I3,2F6.1,2F6.0,2F6.1)' )  &
                 ystid, zfdist, zcut, zzol(ke), zpblsu(icsu,ke), zthol(ke)     &
               , kthr, zthol(kbl), zthol(kthr)                                 &
                     , zzol (kbl), zzol (kthr), t(io,jo,ke,nnew), t_2m(io,jo)
        ENDIF
      ENDIF

      kthr = ke - 1
      DO WHILE ((kthr >  1) .AND. (zsprsu(icsu,kthr) < zcut))
        kthr = kthr - 1
      ENDDO
      kthr = kthr + 1

      !   If this range ends at a model level, where the obs increments are
      !   spread along isentropes, then the range of levels influenced by
      !   the obs may be anywhere between the top and lowest level where
      !   spreading is along isentropes.
      IF (msprpsu == 2) kthr = MIN( kthr , ktth )

      !   if (msprpar == 1) : Wanted:
      !   the (upper limit for the) uppermost model level, which intersects
      !   height 'zcut' within the horizontal area of influence.
      !   ==>  the height difference between this uppermost level and the obs
      !        level must not be larger than the orography at the obs. location
      !        ==> 'zupper' = 'zcut' + orography at the obs location

      IF ((msprpsu == 1) .AND. (kthr > ktp)) THEN
        zupper  = zcut + hhl(io,jo,ke+1)
        DO WHILE ((kthr >  1) .AND. (zsprsu(icsu,kthr) < zupper))
          kthr = kthr - 1
        ENDDO
        kthr = - MIN( kthr + 1 , ktp + 1 )
      ENDIF
      zrtinfl = MAX( epsy , (c1 + (rhtfsu(ivar) -c1) *(c1- wtvp))              &
                           *rhiflsu(ivar) )
      zvcutsu (icsu,ivrs) = zcut
      zriflsu (icsu,ivrs) = zrtinfl
      kviflsu (icsu,ivrs) = kthr
      kthrmin             = MIN( kthrmin , ABS( kthr ) )
    ELSE
      zvcutsu (icsu,ivrs) = zalllow
      zriflsu (icsu,ivrs) = - c1
      kviflsu (icsu,ivrs) = ke + 1
    ENDIF

    ! determine if linear temporal interpolation possible
    ! ---------------------------------------------------

!   lti2 = (ltipsu) .AND. (wtv >= epsy) .AND. (zwtsu(icsu,1,ivrs) >= epsy)
!   IF (lti2) THEN
!     tabdif = ABS( osghed(nsgob,nhtime) - osghed(nsgob2,nhtime) )
!     lti2   = (ABS(tsgadm(ista,2) - tabdif) < epsy) .AND. (tabdif <= tipmxsu)
!   ENDIF
!   ltisu (icsu,ivrs) = lti2
    ltisu (icsu,ivrs) =       (isgadm(ista,6) == 1)                            &
                        .AND. (wtv >= epsy) .AND. (zwtsu(icsu,1,ivrs) >= epsy)

    ! determine if data from present station will be used
    ! ---------------------------------------------------

    IF (zwtsu(icsu,1,ivrs) > epsy) ltakeob = .TRUE.
    IF (zwtsu(icsu,2,ivrs) > epsy) ltakeob = .TRUE.

  ENDDO loop_over_variables

  !   use only maxsgo-1 obs incr reports to obtain correct counting of nexcesu;
  !   note: usually local number 'icsu' remains << than global maximum 'maxsgo'
  IF ((ltakeob) .AND. (icsu < maxsgo)) THEN
    ktopsu          = MIN( ktopsu , kthrmin )
    isgadm (ista,8) = icsu
  ELSE
    IF ((ltakeob) .AND. (icsu >= maxsgo))  nexcesu = nexcesu + 1
    DO ivar   = 1 , 4
      xoisu (icsu,1,ivar) = c0
      xoisu (icsu,2,ivar) = c0
    ENDDO
    icsu            = icsu - 1
    isgadm (ista,8) = 0
  ENDIF

ENDDO loop_over_surface_stations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  nsutot = icsu

  IF (ldump_ascii) THEN
    ! flush YUPRINT file
    IF ((lwonl) .AND. (     (lfirst)                                           &
                       .OR. (rmod( ntstep*dt,c3600,epsy ) < tconbox))) THEN
      CLOSE (nupr)
      OPEN  (nupr,FILE=yuprint,FORM='FORMATTED',STATUS='OLD'                   &
                              ,POSITION='APPEND',IOSTAT=nstat)
      IF (nstat /= 0) PRINT '("OPENING OF FILE yuprint FAILED, surf_local_info")'
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure surf_local_info
!-------------------------------------------------------------------------------

END SUBROUTINE surf_local_info


!-------------------------------------------------------------------------------
!+ Module procedure to extrapolate screen-level obs to the lowest model level
!-------------------------------------------------------------------------------

SUBROUTINE surf_obs ( xobs , yobs , ivrs , xke , yke , zobop )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure extrapolates a surface-level observation
!   (of horizontal wind, temperature, or relative humidity) from the observation
!   level (10m or 2m above the ground at station height) to the lowest model
!   level (i.e. an approximation to the inverse of the observation operator
!   is applied).
!
! Method:
!   1.: Height correction, to account for differences between station height
!       and model orography:
!       - for wind and relative humidity, no height correction is applied
!       - temperature is corrected using the lapse rate between about 1500m
!         and 3000m above the ground (i.e. the lapse rate assumed representa-
!         tive for the conditions above the atmospheric boundary layer ABL)
!   2.: Extrapolation of the value at surface-level ('xsl') to the lowest model
!       level (value 'xke'):
!       - for wind                , the non-dimensional turbulence transfer
!         coefficient 'Cj / ABS(vke)' (= model variable 'tcj', j = m or h) is
!         in 1st order independent from the wind speed 'ABS(vke)' at the lowest
!         model level (the dependence is particularly weak near neutral condi-
!         tions). Hence for the computation of 'Cj / ABS(vke)' , model values
!         of wind are used (by using 'tcj') instead of (unknown) 'observed'
!         winds. With this approximation, the (inverse) observation operator,
!         i.e. the relationship between 'xsl' and 'xke', becomes linear.
!       - for temperature             , the model value of dry static energy
!         is interpolated to 2m as in 'near_surface', and converted to a
!         temperature value at 2m. Then, the observation increment of the
!         temperature at 2m is determined, and added to the model value at the
!         lowest model level.
!    \  - for humidity, the model value of specific water content is interpola-
!    \    ted to 2m in the same way as the dry static energy, and converted to a
!    \    relative humidity value at 2m. Then, the observation increment of the
!    \    relative humidity at 2m is determined, and added to the model value 
!    \    at the lowest model level.
!       - for relative humidity, the values at the 2 levels are assumed to be
!         equal (as is done in proc. 'near_surface)
!       Note, that for temperature and relative humidity, the observation incre-
!       ments at 2m and at the lowest model level are assumed to be equal!
!
! Written by        :  Christoph Schraff, DWD  (original version: 21.10.97)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  REAL    (KIND=wp   ),     INTENT (IN)            ::     &
    xobs             ,& ! (observed) quantity at the observation level
                        ! (10m or 2m above the ground at station height)
    yobs                ! for vector quantity: second component at obs. level
                        ! for scalar quantity: model value at lowest model level

  INTEGER (KIND=iintegers), INTENT (IN)            ::     &
    ivrs                ! index of observed quantity: 1=(u,v); 2=T; 3=RH

  REAL    (KIND=wp   ),     INTENT (OUT)           ::     &
    xke              ,& ! quantity 'xobs' extrapolated to lowest model level
    yke                 ! for vector quantity: 'yobs' extrapol. to lowest level
                        ! for scalar quantity: model value extrapolated to
                        !                      obs. station at 2 m

  REAL    (KIND=wp   ),     INTENT (OUT), OPTIONAL ::     &
    zobop               ! observation operator for 10 m horizontal wind

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    nsgob               ! index of reports in the ODR (obs data record)

  REAL    (KIND=wp   )     ::  &
    xsl              ,& ! quantity 'xobs' extrapolated from the surface level
                        ! at the station height to that at the model orography
    ysl              ,& ! second component of a vector quantity at the surface
                        ! level of the model orography
    zz0o             ,& ! roughness length [m] at the obs. location
    zcm              ,& ! non-dim. turbulence transfer coeff. for momentum
    zch              ,& ! non-dim. turbulence transfer coeff. for static energy
    zsqcm            ,& ! SQRT( zcm )
    zchdcm           ,& ! zch / zsqcm / 'von Karman constant'
    z10              ,& ! 10m
    zhmlke           ,& ! height at the lowest model level
    zzke             ,& ! height differ. betw. lowest model level and orography
    zh               ,& ! dry static energy at the lowest model level
    zhs              ,& ! dry static energy at the model orography
    zlapse              ! lapse rate used for the height correction of temperat.
!   zq2mod              ! model value of specific water content at 2m
!   zpsl                ! model pressure at 2m

  REAL    (KIND=wp   )     , SAVE ::  &
    zfext            ,& ! factor, by the difference of dry static energy between
                        ! surface level (2m) and ground level is multiplied to
                        ! render the difference in dry static energy between
                        ! the lowest model level and the ground level
    zt2mod              ! model temperature at 2m
!   zr2mod              ! model value of relative humidity at 2m


! Local (automatic) arrays: None
! -------------------------

  INTEGER (KIND=iintegers) , SAVE  :: &
    istextt  =  0    ,& ! index of last station for which T-2m has been computed
    istextq  =  0    ,& ! index of last station for which RH-2m has been comput.
    ntext    =  0       ! last time step when 'istextq' has been computed

!
!------------ End of header ----------------------------------------------------


 
!-------------------------------------------------------------------------------
! Begin Subroutine surf_obs
!-------------------------------------------------------------------------------

  IF (ntstep /= ntext) istextt = 0
  IF (ntstep /= ntext) istextq = 0

!-------------------------------------------------------------------------------
!  Section 1: Horizontal wind
!-------------------------------------------------------------------------------

IF (ivrs == 1) THEN

! (no) height correction
! ----------------------

  xsl    =  xobs
  ysl    =  yobs

! extrapolation to lowest model level (by inversion of the observation operator)
! -----------------------------------

! check tcm and gz0
  IF ((lphys) .AND. (ltur)) THEN
    zz0o  =  gz0(io,jo) / g
    zcm   =  MAX( tcm(io,jo) , 5.0E-4_wp )
  ELSE
    zz0o  =  MAX( gz0(io,jo) , 1.0E-3_wp ) / g
    zcm   =  1.0E-4_wp
  ENDIF
  z10 = 10.0_wp

! other variables
  zhmlke = c05 * (hhl(io,jo,ke) + hhl(io,jo,ke+1))
  zh     = cp_d * t  (io,jo,ke,nnew) + g* zhmlke
  zhs    = cp_d * t_g(io,jo   ,nnew) + g* hhl(io,jo,ke+1)
  zsqcm  = SQRT( zcm )
  zzke   = zhmlke - hhl(io,jo,ke+1)

  IF (zhs <= zh) THEN
! stable case
    zobop  =              LOG( (z10  +zz0o) /zz0o )                            &
            - z10 /zzke * LOG( (zzke +zz0o) /zz0o )
    zobop  =  z10 /zzke  +  zsqcm /akt * zobop
  ELSE
! instable case
    zsqcm  =  MAX ( zsqcm  , 0.01_wp )
    zobop  =  c1 - zsqcm/akt * LOG( c1 + (EXP( akt/zsqcm ) -c1)                &
                          *MAX(0.0_wp,zzke-z10) *zz0o /(zzke *(z10+zz0o)))
  ENDIF
  xke    =  xsl / zobop
  yke    =  ysl / zobop


!-------------------------------------------------------------------------------
!  Section 2: Temperature
!-------------------------------------------------------------------------------

ELSEIF (ivrs == 2) THEN

! height correction
! -----------------

  zlapse =       (  t  (io,jo,kml700,nnew) - t  (io,jo,kml850,nnew))           &
          * c2 / (  hhl(io,jo,kml700)      - hhl(io,jo,kml850)                 &
                  + hhl(io,jo,kml700+1)    - hhl(io,jo,kml850+1))
! (alternative lapse rates : constant : 0.0065 (e.g. for height diff > 500m , or
!  (Damrath): noon: 0.0110 ; midnight : 0.0045  --> with temporal interpolation)
!      ( for  Tmax: 0.0130 ; for Tmin : 0.0020 )
  nsgob  = isgadm(ista,itim)
  xsl    =  xobs + zlapse *(hhl(io,jo,ke+1) - osghed(nsgob,nhalt))
  ysl    =  yobs - zlapse *(hhl(io,jo,ke+1) - osghed(nsgob,nhalt))

! extrapolation to lowest model level (by constant extrapolation of the obser-
! -----------------------------------  vation increment at 2m above the ground)

! if the model value of T-2m and the factor 'zfext' are not yet available
  IF (ista /= istextq) THEN

! check tcm, tch and gz0
    IF ((lphys) .AND. (ltur)) THEN
      zz0o  =  gz0(io,jo) / g
      zch   =  MAX( tch(io,jo) , 4.0E-5_wp )
      zcm   =  MAX( tcm(io,jo) , 5.0E-4_wp )
    ELSE
      zz0o  =  MAX( gz0(io,jo) , 1.0E-3_wp ) / g
      zch   =  1.0E-4_wp
      zcm   =  1.0E-4_wp
    ENDIF

! other variables
    zhmlke = c05 * (hhl(io,jo,ke) + hhl(io,jo,ke+1))
    zh     = cp_d * t  (io,jo,ke,nnew) + g* zhmlke
    zhs    = cp_d * t_g(io,jo   ,nnew) + g* hhl(io,jo,ke+1)
    zsqcm  = SQRT( zcm )
    zchdcm = zch / (akt * zsqcm)
    zzke   = zhmlke - hhl(io,jo,ke+1)

    IF (zhs <= zh) THEN
! stable case
      zfext = 4.0_wp / (c1 + 4.5_wp /(zzke -c05))
    ELSE
! instable case
      zchdcm    = MAX ( zchdcm , 0.01_wp )
      zfext = c1/ (c1 - zchdcm *LOG( c1 + (EXP( c1/zchdcm ) -c1)             &
                                         *(zzke -c2) *zz0o /(zzke *(c2+zz0o)) ))
    ENDIF

    zt2mod  = cpdr * ((zhs + (zh - zhs) /zfext) - (c2+hhl(io,jo,ke+1)) *g)

  ENDIF

! xke    =  t_g(io,jo,nnew)  +  zfext * (xsl -t_g(io,jo,nnew) + c2 *g *cpdr)  &
!                           -  zzke *g *cpdr
!   (xsl -t_g(io,jo,nnew) + c2*g*cpdr) : static energy diff. betw. 2m and ground
!   zfext *(  dito  )  : static energy diff. betw. lowest model level and ground

  xke      =  yobs  +  xsl - zt2mod
  yke      =  zt2mod  +  ysl - yobs

  istextt  =  ista
  ntext    =  ntstep


  IF (      (lwonl) .AND. (ntstep <= NINT(c2*tconbox/dt))                      &
      .AND. ((ystid(1:3) == '066') .OR. (ystid(1:3) == '075')))                &
    WRITE( nupr,'(A ,''T2m to Tke'',2F6.0,F7.3,F5.2,7F6.1)' )                  &
           ystid, hhl(io,jo,ke+1), osghed(nsgob,nhalt), zlapse, zfext          &
         , xobs, xsl, t(io,jo,ke,nnew), t_g(io,jo,nnew), zt2mod, yobs, xke


!-------------------------------------------------------------------------------
!  Section 3: Relative humidity
!-------------------------------------------------------------------------------

ELSEIF (ivrs == 3) THEN

! (no) height correction
! ----------------------

  xsl    =  xobs
  ysl    =  yobs

! extrapolation to lowest model level (by constant extrapolation of observed
! -----------------------------------  relative humidity at 2m above the ground)

  xke    =  xsl
  yke    =  ysl

! ELSEIF (ivrs == 3) THEN

! (no) height correction
! ----------------------

! xsl    =  xobs

! extrapolation to lowest model level (by constant extrapolation of the obser-
! -----------------------------------  vation increment at 2m above the ground)

!! xke    =  xsl

! if the model value of T-2m and the factor 'zfext' are not yet available
! IF ((ista /= istextq) .AND. (ista /= istextt)) THEN

! check tcm, tch and gz0
!   IF ((lphys) .AND. (ltur)) THEN
!     zz0o  =  gz0(io,jo) / g
!     zch   =  MAX( tch(io,jo) , 4.0E-5_wp )
!     zcm   =  MAX( tcm(io,jo) , 5.0E-4_wp )
!   ELSE
!     zz0o  =  MAX( gz0(io,jo) , 1.0E-3_wp ) / g
!     zch   =  1.0E-4_wp
!     zcm   =  1.0E-4_wp
!   ENDIF

! other variables
!   zhmlke = c05 * (hhl(io,jo,ke) + hhl(io,jo,ke+1))
!   zh     = cp_d * t  (io,jo,ke,nnew) + g* zhmlke
!   zhs    = cp_d * t_g(io,jo   ,nnew) + g* hhl(io,jo,ke+1)
!   zsqcm  = SQRT( zcm )
!   zchdcm = zch / (akt * zsqcm)
!   zzke   = zhmlke - hhl(io,jo,ke+1)

!   IF (zhs <= zh) THEN
! stable case
!     zfext = 4.0_wp / (c1 + 4.5_wp /(zzke -c05))
!   ELSE
! instable case
!     zfext = c1/ (c1 - zchdcm *LOG( c1 + (EXP( c1/zchdcm ) -c1)               &
!                                        *(zzke -c2) *zz0o /(zzke *(c2+zz0o)) ))
!   ENDIF
!    ! model temperature at 2m
!   zt2mod  = cpdr * ((zhs + (zh - zhs) /zfext) - (c2+hhl(io,jo,ke+1)) *g)
! ENDIF

! if the model value of RH-2m is not yet available
! IF (ista /= istextq) THEN
!    ! model value of specific water content at 2m
!   zq2mod = qv_s(io,jo) + (qv(io,jo,ke) - qv_s(io,jo,nnew)) / zfext
!    ! pressure at 2m computed as 'ps' in 'calps', but EXP(x) replaced by (1+x)
!   zpsl = (p0(io,jo,ke) + pp(io,jo,ke,nnew))                                  &
!         *(c1 +  (c05 - c2/(hhl(io,jo,ke) -hhl(io,jo,ke+1))) * dp0(io,jo,ke)  &
!               / (  t(io,jo,ke,nnew) * (c1 + rvd_m_o *qv(io,jo,ke)            &
!                                           - qc(io,jo,ke) - qrs(io,jo,ke))    &
!                  * r * rho0(io,jo,ke)) )
!    ! model value of relative humidity at 2m
!   zr2mod = fq2pv( zq2mod , zpsl , rdv ) / fpvsw( zt2mod, b1, b2w, b3, b4w )
! ENDIF

! xke      =  yobs  +  xsl - zr2mod

! xke      =  MAX( xke , 0.01_wp )

! istextq  =  ista
! ntext    =  ntstep


! IF (      (lwonl) .AND. (ntstep <= NINT(c2*tconbox/dt))                      &
!     .AND. ((ystid(1:3) == '066') .OR. (ystid(1:3) == '075')))                &
!   WRITE( nupr,'(A ,''RH2m to RHke'',3I5,6X, 2F5.2,3F7.4,2F6.2)' )            &
!          ystid, ista, istextt, istextq                                       &
!        , xobs, zr2mod, qv(io,jo,ke), qv_s(io,jo,nnew), zq2mod, yobs, xke

ENDIF


!-------------------------------------------------------------------------------
! End of module procedure surf_obs
!-------------------------------------------------------------------------------

END SUBROUTINE surf_obs


!-------------------------------------------------------------------------------
!+ Module procedure to organize obs increments and QC for upper-air single-level
!-------------------------------------------------------------------------------

SUBROUTINE upair_local_info ( nexceua )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure prepares the lateral spreading of observation
!   increments from upper-air single-level data by computing all the required
!   'local' information (i.e. information on the observation, and further
!   parameters from its location) which will later be broadcast to other
!   processors.
!
! Method:
!   Straightforward computation of observation increments, temporal weights,
!   and other quantities. Threshold quality control. Vertical interpolation
!   linear in ln(p).
!
! Written by        :  Christoph Schraff, DWD  (original version: 14.07.97)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER (KIND=iintegers) , INTENT (OUT)    ::  &
    nexceua             ! number of obs increment reports exceeding array size

! Local parameters:
! ----------------

  INTEGER (KIND=iintegers) , PARAMETER  ::  &
    ndqc  =  3          ! max. number of QC control messages per report
                        !   ( = first dimension of 'zyqc')

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    itic             ,& ! time index of report
    itima  , itime   ,& ! range of loop over existing reports at current station
    icua             ,& ! index of obs. sta. in data record to be used further
    ivrs             ,& ! index of observed quantity: 1=(u,v); 2=T; 3=RH
    ivar             ,& ! index of observed quantity: 1=(u,v); 3=T; 4=RH
    nsgob  , nsgob2  ,& ! indices of reports in the ODR (obs data record)
    nsgoba           ,& ! non-zero index of report in the ODR (obs data record)
    kobtyp , kcdtyp  ,& ! CMA observation type, observation code type
    ka     , kb      ,& ! model levels adjacent to the obs.
    kat    , kbt     ,& ! model levels adjacent to the standard atm. tropopause
    kal    , kbl     ,& ! model levels adjacent to equidistant points at obs loc
    nbsx             ,& ! ODR indices for present variable
    nvrx             ,& ! bit pos. for status 'active' for variable
    nhvcbx           ,& ! /  correction factor to ver-   \  -  at base of report
    nhvctx           ,& ! \  tical correlation scale     /  -  at top  of report
    nlqc             ,& ! number of new QC messages
    ntqc             ,& ! number of total QC messages
    nqc              ,& ! loop index for QC messages
    nwhere           ,& ! indicator how a temporal nudging weight is computed
    kthr             ,& ! level index at upper limits where use horiz infl radii
    ilva   , ilvb    ,& ! indices of correlation scale table levels that are
                        ! adjacent to the specified model level
    mv     , kk      ,& ! loop indices
    ni2    , nbotop  ,& ! loop indices
    nstat               ! error flag

  REAL    (KIND=wp   )     ::  &
    wtv    , wtv1    ,& ! temporal nudging weight
    wtvp             ,& ! min. non-zero temporal weight for 1 sta. & variable
    timdif           ,& ! time difference to the observation
    tabdif           ,& ! time distance between the 2 obs at same obs. station
    ztvob            ,& ! virtual temperature at obs. point
    zlop0            ,& ! log of reference pressure p0r
    zloptp           ,& ! log of pressure at tropopause of standard atmosphere
    zlopf            ,& ! weight factor for vertical interpol. to the obs. point
    zlptf            ,& ! weight factor for vertical interpol. to the tropopause
    zflop            ,& ! weight factor for vertical interpol. to equidist. pts.
    ztp              ,& ! value of spreading parameter at standard tropopause
    zdvia  , zdvib   ,& ! vertical distance of model levels from obs. level
    psign            ,& ! sign, depending if lower or upper cut-off is computed
    svcutos          ,& ! adjusted vertical cut-off radius
    zcut             ,& ! spreading parameter at vertical cut-off
    zupper           ,& ! upper height limits where horiz. radii of infl. used
    zlopco           ,& ! log( pressure ) at level relevant for area of influen.
    zf               ,& ! weight factor for vertical interpol. to that level
    rvinfl           ,& ! vertically dependent part of horizontal correl. scale
    zrtinfl             ! actual horizontal correl. scale

  LOGICAL                  ::  &
    lqc              ,& ! threshold quality control to be applied
!   lact             ,& ! current report is (or can become) active
!   lti2             ,& ! linear temporal interpolation as time window
    ltakeob          ,& ! use data of present station for further processing
    lqcflag (2)      ,& ! set threshold quality control flags in the ODR
    lqcfrst (2)      ,& ! QC needs to be done since QC status is not yet defined
    lveridat         ,& ! fill simulated obs body (for writing to feedback file)
    lwrqc               ! data with obs time close to current time, which are
                        !   rejected by the threshold quality control, are
                        !   printed to a separate file at the current timestep

  INTEGER (KIND=iintegers) :: izerror
  CHARACTER (LEN=255)      :: yzerrmsg
  CHARACTER (LEN=25)       :: yzroutine = 'upair_local_info'

! Local (automatic) arrays:
! -------------------------

  REAL    (KIND=wp   )     ::  &
    svcutof (3)      ,& ! vertical cut-off radii
    zzol    (ke)     ,& ! height on model levels at the obs. location
    zthol   (ke)     ,& ! potential temperature on model levels at the obs loc.
    vip     (2,4)    ,& ! simulated obs = vertically interpolated model values
    zyqc    (ndqc,11)   ! information for QC control messages
!
!------------ End of header ----------------------------------------------------

 
!-------------------------------------------------------------------------------
! Begin Subroutine upair_local_info
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1: Initialisations
!-------------------------------------------------------------------------------

! retreive required microphysics tracer fields
! --------------------------------------------

  ! Retrieve the required microphysics tracers (at nnew)
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnew, ptr = qv)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nnew, ptr = qc)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

! initialisations related to observations
! ---------------------------------------

  svcutof (1) = SQRT( vcorls(1) * vcutof(1) )
  svcutof (2) = SQRT( vcorls(3) * vcutof(3) )
  svcutof (3) = SQRT( vcorls(4) * vcutof(4) )

  DO ni2      = 1 , 2
    DO ivrs   = 1 , 3
      DO ista = 1 , maxuso
        zwtua   (ista,ni2,ivrs) = c0
        xoiua   (ista,ni2,ivrs) = c0
        fcorlua (ista,ni2,ivrs) = c1
        zvcutua (ista,ni2,ivrs) = c0
        zriflua (ista,ni2,ivrs) = c0
        zqualua (ista,ni2,ivrs) = c0
        kviflua (ista,ni2,ivrs) = 0
      ENDDO
    ENDDO
  ENDDO
  DO ista = 1 , maxuso
    xoiua   (ista,1,4) = c0
    xoiua   (ista,2,4) = c0
    ltiua   (ista  ,1) = .FALSE.
    ltiua   (ista  ,2) = .FALSE.
    ltiua   (ista  ,3) = .FALSE.
    zpobua  (ista    ) = c1
    zzobua  (ista    ) = c0
    zsprob  (ista    ) = c0
    zsprtp  (ista    ) = c0
    zrtdgua (ista    ) = c0
    zvidua  (ista    ) = c0
    kobtyua (ista    ) = 0
    kcdtyua (ista    ) = 0
  ENDDO
  DO kk = 1 , ke
    DO ista = 1 , maxuso
      zsprua  (ista,kk) = c0
      zlopua  (ista,kk) = c0
      IF (lnisua) znisua  (ista,kk) = c0
    ENDDO
  ENDDO

  zlop0  = LOG( p0r )
  zloptp = LOG( ptropop )

  ! initialise loop over upper-air single-level 'stations'
  ! ------------------------------------------------------

  icua = 0
  nexceua = 0
 
loop_over_upper_air_stations:  DO ista = nsusta + 1 , nsusta + nuasta
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  icua = icua + 1
  ltakeob = .FALSE.

  io = isgadm(ista,3)
  jo = isgadm(ista,4)
  nsgoba = MAX( isgadm(ista,1) , isgadm(ista,2) )
  io_tot = mosghd(nsgoba,nhitot)
  jo_tot = mosghd(nsgoba,nhjtot)

!-------------------------------------------------------------------------------
!  Section 2: This part is required for nudging and LETKF / feedback files:
!             compute the simulated observations H(x) by applying the forward
!             observation operator, and perform the threshold quality control
!             for individual observations (first guess check)
!-------------------------------------------------------------------------------

  kobtyp = mosghd(nsgoba,nhobtp)
  kcdtyp = mosghd(nsgoba,nhcode)
  ystid  = yosghd(isgadm(ista,5)) (1:ilstidg)

! check - if the report enters this routine for the first time and the hence
!            the QC status is not defined yet and QC needs to be done, and
!       - if the threshold quality control flags are to be written to the ODR
! ---------------------------------------------------------------------------

! DO itic = 1, 2
!   nsgob  =  isgadm(ista,itic)
!
!!  ( if lqcfrst then (MOD(*nhqcfw,8) < 4) or (*nhqcfw == 0) ;
!!    if (MOD(*nhqcfw,8) < 4) then lqcfrst )
!   lqcflag (itic) = (nsgob > 0) .AND. (mosghd(MAX(nsgob,i1),nhqcfw) >=  0)
!   lqcfrst (itic) = (nsgob > 0) .AND. (mosghd(MAX(nsgob,i1),nhqcfw) <= -1)    &
!                                .AND. (mosghd(MAX(nsgob,i1),nhqcfw) > -99)    &
!                        .AND. (MOD(ABS(mosghd(MAX(nsgob,i1),nhqcfw)),8) < 4)
!   IF (lqcfrst(itic))  mosghd (nsgob,nhqcfw) =      mosghd(nsgob,nhqcfw) - 4
! ENDDO

  kb = 0

  itima = 2 - MIN( 1 , isgadm(ista,1) )
  itime = 1 + MIN( 1 , isgadm(ista,2) )

  loop_over_time_index:  DO itim = itima , itime

    !   nsgob > 0
    nsgob  =  isgadm(ista,itim)

    ! check - if the report enters this routine for the first time and the hence
    !            the QC status is not defined yet and QC needs to be done, and
    !       - if the threshold quality control flags shall be written to the ODR
    ! --------------------------------------------------------------------------
    !   ( if lqcfrst then (MOD(*nhqcfw,8) < 4) or (*nhqcfw == 0) ;
    !     if (MOD(*nhqcfw,8) < 4) then lqcfrst )
    lqcflag (itim)  =                 (mosghd(nsgob,nhqcfw) >=    0)
    lqcfrst (itim)  =                 (mosghd(nsgob,nhqcfw) <=   -1)           &
                       .AND.          (mosghd(nsgob,nhqcfw) >   -99)           &
                       .AND. (MOD( ABS(mosghd(nsgob,nhqcfw)), 8 ) < 4)
    IF (lqcfrst(itim))  mosghd (nsgob,nhqcfw)  =  mosghd(nsgob,nhqcfw) - 4

    lveridat =  lqcflag(itim) .AND. (lvofoi)
    timdif   =  osghed(nsgob,nhtime) - acthr
    tabdif   =  ABS( timdif )
    lqc      =  (lqcall) .OR. (lqcfrst(itim)) .OR. (lqcflag(itim))
    lwrqc    =  (lqcall) .AND. (ABS( c3600*timdif+epsy ) <= cqcwbox)

    ! forward observation operator to obtain simulated observations (H(x))
    ! --------------------------------------------------------------------

    CALL sing_obs_operator ( ke, a_u(io,jo,:), a_v(io,jo,:), t(io,jo,:,nnew)   &
                           ,     qv(io,jo,:), qc(io,jo,:), a_p(io,jo,:)        &
                           , osgbdy(nsgob,:), lveridat, rdv, b1, b2w, b3, b4w  &
                           , ssgbdy(:,nsgob) , vip(itim,:), kb, zlopf )
  ! ======================

    IF (lqc) THEN

      ! threshold quality control (first guess check)
      ! ---------------------------------------------

      CALL sing_quality_cntl ( osgbdy(nsgob,:), vip(itim,:), ps(io,jo,nnew)    &
                             , kb, kobtyp, tabdif, qcvf, qcc, ndqc, lwrqc      &
                             , mosgbd(nsgob,:), nlqc, zyqc )
    ! ======================

      ! Fill record for later printing for control of QC
      ! ------------------------------------------------
      !   adjust 'nlqc' if space in arrays 'oyqc' is insufficient
      nlqc  =  MIN( maxqcp - ntotqc , nlqc )
      DO nqc = 1 , nlqc
        ntqc  =  ntotqc + nqc
        yyqc (ntqc   ) = ystid
        myqc (ntqc, 1) = mosghd(nsgob,nhcode)
        myqc (ntqc, 2) = NINT( zyqc(nqc,1) )
        IF (.NOT. lwrqc) myqc (ntqc,1) = 0
        oyqc (ntqc, 1) = osghed(nsgob,nhtime)
        oyqc (ntqc, 2) = zyqc (nqc, 2)
        oyqc (ntqc, 3) = osghed(nsgob,nhjlat)
        oyqc (ntqc, 4) = osghed(nsgob,nhilon)
        oyqc (ntqc, 5) = zyqc (nqc, 5)
        oyqc (ntqc, 6) = zyqc (nqc, 6)
        oyqc (ntqc, 7) = zyqc (nqc, 7)
        oyqc (ntqc, 8) = zyqc (nqc, 8)
        oyqc (ntqc, 9) = zyqc (nqc, 9)
        oyqc (ntqc,10) = zyqc (nqc,10)
!       oyqc (ntqc,11) = zyqc (nqc,11)
      ENDDO
      ntotqc  =  ntotqc  +  nlqc
    ENDIF

  ENDDO loop_over_time_index

!-------------------------------------------------------------------------------
!  Section 3: This is used for nudging only:
!             Get the observation increments and temporal weights,
!             as well as additional local information for later spreading
!-------------------------------------------------------------------------------

  ! get:
  !   ioua,joua: location of the obs. station
  !   ystidua  : station id.
  !   zsprob   : values of spreading parameter at the obs. point
  !   zsprtp   : values of spreading parameter at the tropopause
  !   zzobua   : height   at the obs. point
  !   zpobua   : pressure at the obs. point
  !   zrtdgua  : (Rd / g) * Tv (conversion factor ln(p) --> z) at the obs. point
  !   zsprua   : values of spreading parameter at obs. location on model levels
  !   zlopua   :      log( pressure )          at obs. location on model levels
  !   znisuq   : vertical profile of the parameter which models the non-isotropy
  !   zvidua   : effective vertical distance for the vertical interpolation
  !   --------------------------------------------------------------------------

  IF (kb > 0) THEN
    ka = kb - 1
    !   get 'z' at obs level, and profile of 'z' at model levels
    DO kk = 1 , ke
      zzol        (kk) = c05 * (hhl(io,jo,kk) + hhl(io,jo,kk+1))
      zsprua (icua,kk) = zzol(kk)
      zlopua (icua,kk) = LOG( a_p(io,jo,kk) )
    ENDDO
    zsprob (icua) = (c1-zlopf) *zzol(kb) + zlopf *zzol(ka)
    zzobua (icua) = zsprob(icua)
    zpobua (icua) = osgbdy(nsgoba,nbsp)
    ztvob =  (c1-zlopf) *t(io,jo,kb,nnew) *( c1 + rvd_m_o *qv(io,jo,kb)        &
                                            -qc(io,jo,kb) -qrs(io,jo,kb))      &
            +    zlopf  *t(io,jo,ka,nnew) *( c1 + rvd_m_o *qv(io,jo,ka)        &
                                            -qc(io,jo,ka) -qrs(io,jo,ka))
    zrtdgua (icua) = r_d * ztvob / g
    !   get 'z' at 'tropopause'
    kbt = 2
    DO WHILE ((kbt < ke) .AND. (a_p(io,jo,kbt) <= ptropop))
      kbt = kbt + 1
    ENDDO
    kat = kbt - 1
    zlptf = (zlopua(icua,kbt) - zloptp) / (zlopua(icua,kbt) - zlopua(icua,kat))
    zsprtp  (icua) = (c1-zlptf) *zzol(kbt) + zlptf *zzol(kat)
    ioua    (icua) = io
    joua    (icua) = jo
    ioua_tot(icua) = io_tot
    joua_tot(icua) = jo_tot
    kobtyua (icua) = kobtyp
    kcdtyua (icua) = kcdtyp
    ystidua (icua) = yosghd(isgadm(ista,5)) (1:ilstidg)
    ystid          = ystidua(icua) (1:ilstidp)
  ENDIF

  !   if necessary: get profile of 'theta' at model levels
  IF ((kb > 0) .AND. ((lnisua) .OR. (msprpar == 2))) THEN
    DO kk = 1 , ke
      zthol (kk) = t(io,jo,kk,nnew) * EXP( rdocp *(zlop0 - zlopua(icua,kk)) )
      IF                 (msprpar == 2)  zsprua (icua,kk) = zthol(kk)
      IF ((lnisua) .AND. (msprpar <= 1)) znisua (icua,kk) = zthol(kk)
      IF ((lnisua) .AND. (msprpar == 2)) znisua (icua,kk) = zzol (kk)
    ENDDO
    !   get 'theta' at obs. level and at tropopause
    IF (msprpar == 2) THEN
      zsprob (icua) = (c1-zlopf) *zthol(kb)  + zlopf *zthol(ka)
      zsprtp (icua) = (c1-zlptf) *zthol(kbt) + zlptf *zthol(kat)
    ENDIF
  ENDIF

  !   if necessary: profiles of 'theta' and 'z' at 'exp(-z) -equidistant' levels
  IF ((kb > 0) .AND. (lnisua) .AND. (msprpar >= 1)) THEN
    DO mv = 1 , mxispr
      kk = ke
      DO WHILE (  (   (msprpar <  2) .AND. (zzol (MAX(kk,i1)) <= zsdni(mv,1))  &
                 .OR. (msprpar == 2) .AND. (zthol(MAX(kk,i1)) <= zsdni(mv,2))) &
              .AND. (kk > 0))
        kk = kk - 1
      ENDDO
      kbl = MIN( INT( kk+1 ,iintegers) , ke )
      kal = MAX(      kk               , i1 )
      zflop = c1
      IF (msprpar < 2) THEN
        IF (kal /= kbl) zflop =   (zzol(kbl) - zsdni(mv,1))                    &
                                / (zzol(kbl) - zzol(kal))
        znisuq (icua,mv) = (c1-zflop) *zthol(kbl)  + zflop *zthol(kal)
      ELSE
        IF (kal /= kbl) zflop =   (zthol(kbl) - zsdni(mv,2))                   &
                                / (zthol(kbl) - zthol(kal))
        znisuq (icua,mv) = (c1-zflop) *zzol(kbl) + zflop *zzol(kal)
        IF (kal == ke)                                                         &
          znisuq (icua,mv) = zzol(ke) - (zthol(ke) -zsdni(mv,2)) /thdzt
      ENDIF
    ENDDO
  ENDIF

  !   effective vertical distance for the vertical interpolation
  IF (kb > 0) THEN
    IF (msprpar <= 1) THEN
      zdvia = (  zsprua(icua,ka) - zsprob(icua)) / zrtdgua(icua)
      zdvib = (- zsprua(icua,kb) + zsprob(icua)) / zrtdgua(icua)
    ELSE
      ztp   = zsprtp(icua)
      zdvia =  MAX( MIN( zsprua(icua,ka),ztp ) - zsprob(icua) , c0 )           &
             + MAX( zsprua(icua,ka) - MAX( zsprob(icua),ztp ) , c0 ) /qst
      zdvib =  MAX( MIN( zsprob(icua),ztp ) - zsprua(icua,kb) , c0 )           &
             + MAX( zsprob(icua) - MAX( zsprua(icua,kb),ztp ) , c0 ) /qst
    ENDIF
    zvidua (icua) = ABS( zdvia * zdvib / MAX( zdvia + zdvib , epsy ) )
  ENDIF

  ! initialise loop over variables
  ! ------------------------------

  loop_over_variables:  DO ivrs = 1, 3

    IF (ivrs == 1) THEN
      ivar   = 1
      nbsx   = nbsu
      nvrx   = nvru
      nhvcbx = nhvcbu
      nhvctx = nhvctu
    ELSEIF (ivrs == 2) THEN
      ivar   = 3
      nbsx   = nbst
      nvrx   = nvrt
      nhvcbx = nhvcbt
      nhvctx = nhvctt
    ELSEIF (ivrs == 3) THEN
      ivar   = 4
      nbsx   = nbsrh
      nvrx   = nvrq
      nhvcbx = nhvcbq
      nhvctx = nhvctq
    ENDIF

    ! get the observation increments
    ! ------------------------------

    IF (      (     ((kcdtyp /= nmodes) .AND. (gnudgar(ivar) > epsy))          &
               .OR. ((kcdtyp == nmodes) .AND. (gnudgms(ivar) > epsy)))         &
        .AND. (kb > 0)) THEN
      DO itim = itima , itime
        nsgob  = isgadm(ista,itim)
        IF (      (      BTEST( mosgbd(nsgob,nbserr), nvrx ))                  &
            .AND. (.NOT. BTEST( mosgbd(nsgob,nbsqcf), nvrx ))) THEN
          xoiua (icua,itim,ivar)  =  osgbdy(nsgob,nbsx) - vip(itim,ivar)
          IF (ivar == 1)  xoiua (icua,itim,2) = osgbdy(nsgob,nbsv) - vip(itim,2)
        ENDIF
      ENDDO
    ENDIF

    ! get temporal weights, taking into account the threshold quality control
    ! -----------------------------------------------------------------------

    DO itic = 1, 2
      nsgob  = isgadm(ista,itic)
      nsgob2 = isgadm(ista,3-itic)
!     IF ((kb == 0) .OR. (gnudgar(ivar) <= epsy) .OR. (nsgob == 0)             &
!                   .OR. (osgbdy(MAX(nsgob,i1),nbsxer) <= rmdich)              &
!                   .OR. (mosghd(MAX(nsgob,i1),nhpass) /= 0)) THEN
      IF (     (nsgob == 0) .OR. (kb == 0)                                     &
          .OR. ((gnudgar(ivar) <= epsy) .AND. (kcdtyp /= nmodes))              &
          .OR. ((gnudgms(ivar) <= epsy) .AND. (kcdtyp == nmodes))              &
          .OR. (mosghd(MAX(nsgob,i1),nhpass) /= 0)                             &
          .OR. (.NOT. BTEST( mosgbd(MAX(nsgob,i1),nbserr), nvrx ))) THEN
        wtv = c0
      ELSE
        nwhere = 2
        wtv = wtsg(ista,itic)
        IF (BTEST( mosgbd(nsgob,nbsqcf), nvrx )) THEN
          wtv = c0
          nwhere = 3
        ELSEIF (     (nsgob2 == 0)                                             &
                .OR. (.NOT. BTEST( mosgbd(MAX(nsgob2,i1),nbserr), nvrx ))      &
                .OR. (      BTEST( mosgbd(MAX(nsgob2,i1),nbsqcf), nvrx ))) THEN
          wtv = wtsg(ista,itic+2)
          nwhere = 4
        ENDIF
        IF (wtv < epsy)  wtv = c0
        IF ((ntstep <= 1) .AND. (lwonl) .AND. (ivar == 1))                     &
          WRITE( nupr,'(''air: wtu'', 3I4,F7.4,F7.2)' )                        &
                 io_tot, jo_tot, nwhere, wtv, xoiua(icua,1,ivar)
!US itim is a loop index above and no longer set here: io_tot, jo_tot, nwhere, wtv, xoiua(icua,itim,ivar)
      ENDIF
      IF (itic == 1)    wtv1 = wtv
      zwtua   (icua,itic,ivrs) = wtv
      IF (wtv >= epsy) zqualua (icua,itic,ivrs) = c1
    ENDDO

!-------------------------------------------------------------------------------
!  Section 4: This is used for nudging only:
!             Get the vertical upper & lower cut-offs for a report,
!             an upper estimate for the range of model levels influenced by a
!             report (except for the lower limit with horizontal spreading), &
!             upper estimates for horizontal correlation scales for target grid
!             pts. at the vertical cut-offs
!-------------------------------------------------------------------------------

    ! get range of model levels within the vertical cut-off at the obs. location
    ! --------------------------------------------------------------------------

    IF ((zwtua(icua,1,ivrs) >= epsy) .OR. (zwtua(icua,2,ivrs) >= epsy)) THEN
      wtvp = c1
      IF (wtv1 >= epsy) wtvp = MIN( wtvp , wtv1 )
      IF (wtv  >= epsy) wtvp = MIN( wtvp , wtv  )
      !   (nbotop == 1) for lower cut-off, (nbotop == 2) for upper cut-off
      DO nbotop = 1 , 2
        psign = REAL ( 2 *nbotop - 3, wp )
        IF (nbotop == 1) fcorlua (icua,1,ivrs) = osghed(nsgoba,nhvcbx)
        IF (nbotop == 2) fcorlua (icua,2,ivrs) = osghed(nsgoba,nhvctx)
        svcutos = svcutof(ivrs) * ( (c1-fsvcut)                                &
                                   +    fsvcut * fcorlua(icua,nbotop,ivrs) )
        IF (msprpar <= 1) THEN
          zcut = zsprob(icua) + psign *zrtdgua(icua) *svcutos
        ELSE
          IF (zsprob(icua) <= zsprtp(icua)) THEN
            zcut  =  zsprob(icua) +   psign *       svcutos
            IF (zcut > zsprtp(icua)+epsy)                                      &
              zcut = zsprtp(icua) + SQRT( (         svcutos          )**2      &
                                         -(zsprtp(icua) -zsprob(icua))**2) *qst
          ELSE
            zcut  =  zsprob(icua) +   psign * qst * svcutos
            IF (zcut < zsprtp(icua)-epsy)                                      &
              zcut = zsprtp(icua) - SQRT( (   qst * svcutos          )**2      &
                                         -(zsprtp(icua) -zsprob(icua))**2) /qst
          ENDIF
        ENDIF
        kthr = ke + (nbotop-1) *(kb-ke)
        DO WHILE ((kthr >  1) .AND. (zsprua(icua,kthr) < zcut))
          kthr = kthr - 1
        ENDDO
        IF ((nbotop == 2) .AND. (zsprua(icua,1) > zcut)) kthr = kthr + 1
        IF ((nbotop == 1) .AND. (kthr == ke))   kthr = ke + 1
        IF ((nbotop == 2) .AND. (msprpar >= 1)) kthr = MIN( kthr , ke )

        !   If this range ends at a model level, where the obs increments are
        !   spread along isentropes, then the range of levels influenced by
        !   the obs may be anywhere between the top and lowest level where
        !   spreading is along isentropes.
        IF ((msprpar == 2) .AND. (kthr >= ktth) .AND. (kthr <= ke)) THEN
          IF (nbotop == 1) kthr = ke
          IF (nbotop == 2) kthr = ktth
        ENDIF

        ! if (msprpar == 1) : Wanted:
        ! upper limits to the 2 horizontal radii of influence used to determine
        ! (an upper limit of) the range of model levels influenced by the obs.
        ! ---------------------------------------------------------------------
        !     Needed: the (lower limit for the) pressure at the obs. location of
        !             the uppermost model level, which intersects height 'zcut'
        !             within the horizontal area of influence.
        !     ==> 1.: the height difference between this uppermost level and the
        !             obs level must not be larger than the orography at the obs
        !             location
        !             ==> 'zupper' = 'zcut' + orography at the obs. location
        !         2.: get pressure at the obs. location of the model level which
        !             is immediately above / below (nbotop == 1 / 2) 'zupper'
        !       IF ('zcut' < orography) THEN lower level limit is known: 'ke+1'.

        IF ((msprpar == 1) .AND. (kthr > ktp) .AND. (kthr <= ke)) THEN
          zupper  = zcut + hhl(io,jo,ke+1)
          DO WHILE ((kthr >  1) .AND. (zsprua(icua,kthr) < zupper))
            kthr = kthr - 1
          ENDDO
          IF ((nbotop == 2) .AND. (zsprua(icua,1) > zupper)) kthr = kthr + 1
          kthr = - MIN( kthr , INT( ktp + 1 ,iintegers) )
        ENDIF
        IF ((kthr < 0) .OR. (nbotop == 2)) THEN
          zlopco = LOG( MAX( a_p(io,jo,MIN(ABS(kthr),ke)) , 1000.0_wp ) )
          ilva = 1
          IF ((zlopco > tabcolp(1)) .OR. (zlopco <= tabcolp(ncolev))) THEN
            IF (zlopco <= tabcolp(ncolev)) ilva = ncolev
            IF (ivrs == 1) rvinfl = rhvsond(ilva)
            IF (ivrs == 2) rvinfl = rhtsond(ilva)
            IF (ivrs == 3) rvinfl = rhqsond(ilva)
          ELSE
            DO WHILE (zlopco <= tabcolp(ilva))
              ilva = ilva + 1
            ENDDO
            ilvb = ilva - 1
            zf   = (tabcolp(ilvb) - zlopco) / (tabcolp(ilvb) - tabcolp(ilva))
            IF (ivrs == 1) rvinfl = (c1-zf) *rhvsond(ilvb) + zf *rhvsond(ilva)
            IF (ivrs == 2) rvinfl = (c1-zf) *rhtsond(ilvb) + zf *rhtsond(ilva)
            IF (ivrs == 3) rvinfl = (c1-zf) *rhqsond(ilvb) + zf *rhqsond(ilva)
          ENDIF

          zrtinfl = MAX(  (c1 + (rhtfac(ivar)-c1) *(c1-wtvp))                  &
                        * (rhinfl(ivar) + rhvfac(ivar) *rvinfl) , epsy )
        ELSE
          zrtinfl = -c1
        ENDIF
        zvcutua (icua,nbotop,ivrs) = zcut
        zriflua (icua,nbotop,ivrs) = zrtinfl
        kviflua (icua,nbotop,ivrs) = kthr
      ENDDO
    ELSE
      zvcutua (icua,1,ivrs) = zallhig
      zvcutua (icua,2,ivrs) = zalllow
      zriflua (icua,1,ivrs) = - c1
      zriflua (icua,2,ivrs) = - c1
      kviflua (icua,1,ivrs) = 0
      kviflua (icua,2,ivrs) = ke + 1
    ENDIF

    ! determine if linear temporal interpolation possible
    ! ---------------------------------------------------

!   lti2 = (ltipol) .AND. (wtv >= epsy) .AND. (zwtua(icua,1,ivrs) >= epsy)
!   IF (lti2) THEN
!     tabdif = ABS( osghed(nsgob,nhtime) - osghed(nsgob2,nhtime) )
!     lti2   = (ABS(tsgadm(ista,2) - tabdif) < epsy) .AND. (tabdif <= tipolmx)
!   ENDIF
!   ltiua (icua,ivrs) = lti2
    ltiua (icua,ivrs) =       (isgadm(ista,6) == 1)                            &
                        .AND. (wtv >= epsy) .AND. (zwtua(icua,1,ivrs) >= epsy)

    ! determine if data from present station will be used
    ! ---------------------------------------------------

    IF (MAX( zwtua(icua,1,ivrs), zwtua(icua,2,ivrs) ) > epsy) ltakeob = .TRUE.

  ENDDO loop_over_variables

  !   use only maxuso-1 obs incr reports to obtain correct counting of nexceua;
  !   note: usually local number 'icua' remains << than global maximum 'maxuso'
  IF ((ltakeob) .AND. (icua < maxuso)) THEN
    isgadm (ista,8) = icua
  ELSE
    IF ((ltakeob) .AND. (icua >= maxuso))  nexceua = nexceua + 1
    DO ivar   = 1 , 4
      xoiua (icua,1,ivar) = c0
      xoiua (icua,2,ivar) = c0
    ENDDO
    icua            = icua - 1
    isgadm (ista,8) = 0
  ENDIF

  IF ((lwonl) .AND. (ntstep <= 1)) THEN
    ystid = yosghd(isgadm(ista,5)) (1:ilstidp)
    nuatot = MAX( icua , i1 )
    WRITE( nupr,'(A ,'': UPAIR_LOCAL'',2I4,L2,I4,F8.3,2F5.2,2I3,2F5.2,2I3 )' ) &
           ystid, icua, isgadm(ista,8), ltakeob, kb, zvidua(nuatot)            &
         , (zwtua(nuatot,ni2,1),ni2=1,2), (kviflua(nuatot,ni2,1),ni2=1,2)      &
         , (zwtua(nuatot,ni2,2),ni2=1,2), (kviflua(nuatot,ni2,2),ni2=1,2)
  ENDIF


ENDDO loop_over_upper_air_stations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  nuatot = icua

  IF (ldump_ascii) THEN
    ! flush YUPRINT file
    IF ((lwonl) .AND. (lfirst)) THEN
      CLOSE (nupr)
      OPEN  (nupr,FILE=yuprint,FORM='FORMATTED',STATUS='OLD'                   &
                              ,POSITION='APPEND',IOSTAT=nstat)
      IF (nstat /= 0) PRINT '("OPENING OF FILE yuprint FAILED, upair_local_info")'
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure upair_local_info
!-------------------------------------------------------------------------------

END SUBROUTINE upair_local_info

!===============================================================================

ELEMENTAL REAL (KIND=wp) FUNCTION fpvsw  ( zt, b1, b2w, b3, b4w )
  !----------------------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  zt, b1, b2w, b3, b4w
  !----------------------------------------------------------------------------
  ! Magnus formula for water:  input  'zt'  : temperature
  !                            output 'fpvsw': saturation water vapour pressure
  ! Magnus formula for ice  :  if constants 'b2i', 'b4i' for ice are used for
  !                            'b2w', 'b4w' in the call
  !----------------------------------------------------------------------------
  !
  fpvsw  =  b1 * EXP( b2w*(zt-b3)/(zt-b4w) )
  !
END FUNCTION fpvsw

!-------------------------------------------------------------------------------

ELEMENTAL REAL (KIND=wp) FUNCTION fq2pv  ( zqv, zp, rdv )
  !--------------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  zqv, zp, rdv
  !---------------------------------------------------------------------------
  ! water vapour pressure (at T = saturation vapour press. at Td)
  ! from specific humidity and air pressure
  !---------------------------------------------------------------------------
  !
  fq2pv  =  MAX( epsy , zqv ) * zp / (rdv + zqv*(c1-rdv))
  ! 
END FUNCTION fq2pv

!-------------------------------------------------------------------------------

ELEMENTAL REAL (KIND=wp) FUNCTION rmod  ( ztti, ziv, epsy )
  !----------------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  ztti, ziv, epsy
  !---------------------------------------------------------------------------
  ! MOD function for positive reals
  !---------------------------------------------------------------------------
  !
  rmod  =  ztti  -  ziv * REAL(INT( ztti/ziv + epsy ),wp) + epsy
  !
END FUNCTION rmod

!-------------------------------------------------------------------------------

ELEMENTAL INTEGER FUNCTION ibit1  ( invar, ibp )
  !---------------------------------------------
  INTEGER (KIND=iintegers)  , INTENT (IN)  ::  invar, ibp
  !---------------------------------------------------------------------------
  ! returns 1 bit at bit position 'ibp' of integer word 'invar'
  !---------------------------------------------------------------------------
  !
  ibit1 = IAND( ISHFT( invar,-ibp ), 1 )
  !
END FUNCTION ibit1

!-------------------------------------------------------------------------------

END MODULE src_sing_local
