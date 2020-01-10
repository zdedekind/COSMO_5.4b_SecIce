!+ Source module for computing local information on conventional multi-level obs
!-------------------------------------------------------------------------------

MODULE src_mult_local

!-------------------------------------------------------------------------------
!
! Description:
!   This module 'src_mult_local' computes all the 'local' information (i.e.
!   observation increments and further parameters on the observations and their
!   location) on the local multi-level data which is required for the spreading
!   of the observational information later on in the nudging.
!   Specific tasks include:
!    - organizing the application of the (forward and partly inverse)
!      observation operator (e.g. by vertical interpolation of model values and
!      observations),
!    - organizing the threshold quality control,
!    - performing a spatial consistency check of derived IWV data, and
!    - building vertical profiles of observation increments for multi-level data,
!      as well as profiles of humidity increments for IWV data.
!   Major procedures related to the observation operators and quality control are
!   out-sourced into modules 'src_obs_operator_conv' and 'src_obs_qc_conv' after
!   V5_0.
!
!   This module contains the following procedures:
!   - mult_obs_increment : build vertical profiles of observation increments
!                          from multi-level data
!   - mult_obs_gps       : build profiles of humidity increm. from GPS IWV data
!   - mult_org_localinfo : organizes all tasks, incl.  obs operators and QC
!   - mult_iwv_horicheck : spatial consistency check for IWV (radiosonde + GPS)
!   - mult_vprof_2_iwv   : compute IWV from vertical humidity profiles
!   Driving routine 'mult_org_localinfo' and routine 'mult_iwv_horicheck' are
!   called by 'organize_nudging' of module 'src_obs_use_org.f90'.
!   The other routines are called by 'mult_org_localinfo'.
!
!   This module also contains elemental functions, formerly statement functions:
!   - fpvsw        : Magnus formula for water: saturation vapour pressure from T
!   - ftd          : inverse of Magnus formula: dewpoint T from vapour pressure
!   - fpv2q        : specific humidity from vapour pressure and air pressure
!   - fq2pv        : vapour pressure from specific humidity and air pressure
!   - ibit1        : returns 1 bit at given bit position of given integer word
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
!  Bug corrections, if no observation increments (at model levels) existed.
!  ANSI violations removed.
! 1.20       1999/01/07 Guenther Doms
!  Renaming of some global variables
! 1.27       1999/03/29 Christoph Schraff
!  Revised ODR format.
! 1.31       1999/07/01 Christoph Schraff
!  Quantities related to MPI communicator 'icomm_world' removed. Definition of
!  variable used as length of allocable arrays moved from module to procedure.
! 1.36       2000/02/24 Christoph Schraff
!  Bug correction: Initialization of some allocated fields.
! 1.38       2000/04/06 Christoph Schraff
!  Threshold checks of hydrostatic height and thickness derived from observed
!  temperature profiles as additional quality control for temperature data.
!  Optional storage of obs increments for VOF. 
! 1.40       2000/05/23 Christoph Schraff
!  Correction to the vertical scale adjustment of multi-level observations.
!  Setting of quality control flags: 1. of temperature due to the height/ thick-
!  ness check, 2. of height in the absence of active temperature data. Marking
!  of height observation error for height data rejected by the quality control.
! 1.42       2000/06/19 Christoph Schraff
!  Inclusion of a multi-level check as additional quality control.
!  Bug correction in the threshold check of hydrostatic height and thickness.
! 2.4        2001/01/29 Christoph Schraff
!  Rejection of aircraft temperature, if wind is rejected, and vice versa.
!  Additional factors to QC thresholds for aircrafts and PILOTs following ECMWF.
!  Bug correction for reports with top obs. level below lowest model level.
!  Bug correction of indices for 'zrtgkml' at quality control of derived height.
! 2.5        2001/06/01 Christoph Schraff
!  Savety test at array allocations.
! 2.18       2002/07/16 Christoph Schraff
!  Bug correction for model set-up without flat model main levels.
! 3.3        2003/04/22 Maria Tomassini + Christoph Schraff
!  Extension for assimilation of GPS-derived IWV.
! 3.6        2003/12/11 Christoph Schraff
!  Horizontal correlation scale for GPS data given by namelist variable: rhfgps.
!  For GPS stations below orography, constant extrapolation of model values at
!  lowest model level in order to determine model-derived values of IWV. 
!  Current model humidity profile instead of model saturation profile used for 
!  ice-to-water correction of observed IWV. Time-dependent QC threshold for IWV.
!  Zero IWV observation increments also used. QC flags also used if .not. lqc.
!  Correction to decrease excess in frquency of zero qv-increments from IWV obs.
! 3.12       2004/09/15 Christoph Schraff
!  Extension to (prepare to) include assimilation of satellite retrievals.
!  New QC option (qcvf(4) > 0) with reduced, stability-dependent quality control
!  thresholds for humidity and a new formulation of multi-level check (like OI).
! 3.15       2005/03/03 Ulrich Schaettler
!  Replaced FLOAT by REAL
! 3.18       2006/03/03 Christoph Schraff
!  Revision of stability-dependent quality control thresholds for humidity.
!  Increase of basic humidity quality control thresholds ( = new GME-OI values).
!  Introduction of spatial consistency check of IWV from TEMP humidity and GPS.
!  Removal of variable 'ltvprcs' (satellite BT may be processed at any timestep)
!  Flushing of output files. To allow for reproducibility, reset of some local
!  information (zwtml, ltiml, mbotlv, mtoplv) if report is not taken.
! 3.21       2006/12/04 Ulrich Schaettler
!  Put statement functions which are not used to a comment
! V3_23        2007/03/30 Ulrich Schaettler
!  Introduced ldump_ascii for flushing the ASCII files
! V3_24        2007/04/26 Christoph Schraff
!  QC VOF flags corrected: flag zero aircraft winds, flag RH if T is flagged.
! V4_4         2008/07/16 Christoph Schraff
!  Modifications to avoid illegal (NEC) operations
! V4_5         2008/09/10 Christoph Schraff
!  QC enforced the first time an obs is used even if 'lobprcs' is not used (when
!  obs are read from NetCDF files), and element 'nhqcfw' set then from -1 to -2.
!  'o??vip' arrays replaced by 'nbsvi?' or 'nhtvi?' elements of ODR.
!  QC (-flag) lines now looks similar for ps-, single-level and multi-level obs.
! V4_7         2008/12/12 Christoph Schraff
!  Small bug corrections on managing IWV check and on control output for GPS.
! V4_10        2009/09/11 Davide Cesari
!  Add characters after a backslash in comments; g95 interprets that as
!  continuation line
! V4_12        2010/05/11 Ulrich Schaettler
!  Renamed t0 to t0_melt because of conflicting names
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_22        2012/01/31 Christoph Schraff
!  - Observed 'temperature variable' (e.g. virtual temperature) instead of
!    always temperature is stored in ODR and now converted into temperature in
!    routine 'mult_vertic_intpol'.
!  - Observations with the same station id as the quality controlled observation
!    are not used any more in the spatial consistency check for IWV.
!  - Refined option for computing net increments and weights separately for
!    different pre-specified sets of observing systems.
!  - Bug fix: avoid array bound checks introduced related to 'nexce??'.
!  - Some variables moved from 'data_nudge_all' to 'data_obs_lib_cosmo', use of
!    indices of SOR (simulated observation record) introduced.
!  - Further preparations to reduce the required data modules in the routines
!    containing the observation operators and the quality control. This includes
!    the introduction of more extended subroutine argument lists.
!  - GPS: Always prefer ZPD-derived IWV over reported IWV. Add QC-flags for ZPD.
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer 
!  Replaced qx-variables by using them from the tracer module
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  Introduced MESSy interface
! V4_28        2013/07/12 Christoph Schraff
!  - Major procedures related purely to the observation operators and quality
!    control (except for spatial consistency check of IWV) are out-sourced into
!    modules 'src_obs_operator_conv' and 'src_obs_qc_conv'.
!    New, slim interfaces:  Apart from constants, all input and all output
!    information for these out-sourced procedures is via subroutine arguments,
!    which are now prepared in 'mult_org_localinfo'.
!  - Surface pressure increment in height and thickness check taken into account
!    now even if surface pressure has been rejected by threshold QC.
!  - Statement functions replaced by elemental or intrinsic functions.
!  - Improved comments.
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
! V5_1         2014-11-28 Christoph Schraff, Oliver Fuhrer
!  Bug fix to avoid array bound violation by 'nivtot' = 0 as array index. (CS)
!  Processing of Mode-S aircraft obs introduced (CS)
!  Replaced ireals by wp (working precision) (OF)
! V5_2         2015-05-21 Ulrich Schaettler
!  Initialize character variable ystidml in any case, because it could be used
!   uninitialized otherwise
! V5_4         2016-03-10 Christoph Schraff
!  Variables related to the AOF interface removed.
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

    istart,       & ! start index for the forecast of w, t, qd, qw and pp
    iend,         & ! end index for the forecast of w, t, qd, qw and pp
    jstart,       & ! start index for the forecast of w, t, qd, qw and pp
    jend,         & ! end index for the forecast of w, t, qd, qw and pp

! 4. constants for the horizontal rotated grid and related variables
! ------------------------------------------------------------------

    degrad,       & ! factor for transforming degree to rad

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

    t0_melt,      & ! melting temperature of ice
    r_d,          & ! gas constant for dry air
    r_v,          & ! gas constant for water vapor
    rdv,          & ! r_d / r_v
    o_m_rdv,      & ! 1 - r_d/r_v
    rvd_m_o,      & ! r_v/r_d - 1
    rdocp,        & ! r_d / cp_d
    g,            & ! acceleration due to gravity

! 3. constants for parametrizations
! ---------------------------------

    b1,           & ! variables for computing the saturation vapour pressure
    b2w, b2i,     & ! over water (w) and ice (i)
    b3,           & !               -- " --
    b4w, b4i        !               -- " --

! end of data_constants

!-------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------

    rho0       ,    & ! reference density at the full model levels    (kg/m3)
    dp0        ,    & ! reference pressure thickness of layers        ( Pa)
    p0         ,    & ! reference pressure at main levels             ( Pa)
    hhl        ,    & ! geometical height of half model levels        ( m )

! 3. prognostic variables                                             (unit)
! -----------------------

    u          ,    & ! zonal wind speed                              ( m/s )
    v          ,    & ! meridional wind speed                         ( m/s )
    t          ,    & ! temperature                                   (  k  )
    pp         ,    & ! deviation from the reference pressure         ( pa  )

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------

    qrs               ! precipitation water (water loading)           (kg/kg)

! end of data_fields

!-------------------------------------------------------------------------------

USE data_parallel    , ONLY : my_cart_id 
    
!-------------------------------------------------------------------------------

USE environment      , ONLY : model_abort

!-------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------

    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
!   nnow,         & ! corresponds to ntstep
    nnew,         & ! corresponds to ntstep + 1

! 7. additional control variables
! -------------------------------

    ldump_ascii     ! for flushing (close and re-open) the ASCII files

! end of data_runcontrol 

!-------------------------------------------------------------------------------

USE data_nudge_all , ONLY :   &

! 1. Parameters and related variables
! -----------------------------------

    lwonl        ,& ! .TRUE if grid pt. (ionl ,jonl ) lies in the local domain
    lfirst       ,& ! .TRUE if 'organize_nudging' is called the first time
    lvofoi       ,& ! .TRUE if observation increments are also written to VOF
    acthr        ,& ! actual forecast hour
    aiwthr          ! model time [hrs] for which temporal weights are valid

USE data_nudge_all , ONLY :   &

! 2. Namelist variables controlling the data assimilation
! -------------------------------------------------------

    lverif       ,& ! .f. : on - off switch for verification
    mruntyp      ,& ! -1  : type of current model run used for increments in VOF
    tconbox      ,& ! 6*dt: timestep [s] for computing analysis increments
                    !       (i.e. time box of constant analysis increments)
    gnudg        ,& ! 6,12,6,6*10^-4: nudging coefficients for TEMP / PILOT data
    gnudgar      ,& ! 6,0,6,0 *10^-4: nudging coeffic. for AIRCRAFT data [1/s]
    gnudgms      ,& ! 6, 0,6,0*10^-4: nudging coeffic. for Mode-S aircraft [1/s]
    gnudggp      ,& !        0*10^-4: nudging coef. for GPS-derived IWV [1/s]
!   ltipol       ,& ! .t. : .t. ==> linear interpolation in time of upper-air
!                   !               data which are less than 'tipolmx' hrs apart
!   ltipsu       ,& ! .f. : .t. ==> linear interpolation in time of surface-lev.
!                   !               data which are less than 'tipmxsu' hrs apart
!   tipolmx      ,& ! 1.0 : max. time span (hrs) allowed for linear interpolat.
!                   !       for upper-air data (set tipolmx = 0, if .NOT ltipol)
!   tipmxsu      ,& ! 1.0 : max. time span (hrs) allowed for linear interpolat.
!                   !       of surface-level data  (tipmxsu = 0, if .NOT ltipsu)
    msprpar      ,& ! 1   : switch specifying the surfaces along which observat.
                    !       increments of upper-air data are (primarily) spread 
    vcorls       ,& ! 2*.333,: square of the vertical correlation scale,
                    ! 2*.04    i.e. of the Gaussian vertical influence 'radius'
                    !          in log( pressure ) if (msprpar <= 1), or
                    !          in potential temperature if (msprpar == 2)
    vcutof       ,& ! 2* .75 : cut-off of the vertical correlation.       Units:
                    ! 2*1.     value of correlation at cut-off is [exp(-vcutof)]
    rhinfl       ,& ! 0.,70.,: constant part of the 'correlation scale of the
                    ! 0.,0.    autoregressive horiz. correlation function'
                    !          (='COSAC') [km]
    rhvfac       ,& ! 1., 0.,: multiplication factor to the vertically varying
                    ! 2* .83   part of the 'COSAC' (as def. in 'data_nudge_all')
    rhtfac       ,& ! 1.3,   : scaling factor of the total 'COSAC' for the
                    ! 1.43,    beginning and end of the nudging period for 1 obs
                    ! 1.3,     relative to the 'COSAC' at the obs. time as given
                    ! 1.3      by 'rhinfl', 'rhvfac')
    rhfgps       ,& ! 0.45   : scaling (reduction) factor of the total 'COSAC'
                    !          for humidity derived from GPS IWV
    cutofr       ,& ! 4* 3.5 : cut-off in 'COSAC' units of the horizontal
                    !          correlation function
    vcsni        ,& ! 4*2500.: square of Gaussian vertical influence 'radius'
                    !          in potential temperature (if msprpar <= 1) or
                    !          log( pressure ) (if msprpar == 2) on surfaces
                    !          along which obs. increments are spread laterally
    lscadj       ,& ! .T.,   : .F. ==> linear vertical interpolation (in log(p))
                    ! .T.,     instead of vertical scale adjustment (by vertical
                    ! .T.,     averaging over the model layer) for conveying the
                    ! .F.      observational information to the model levels
    topobs       ,& !  849., : threshold [hPa]: above this level (p < topobs),
                    ! 1099.,   only obs. increments at model levels are used,
                    !  799.,699.  i.e. obs. incr. at obs. levels are not used
    botmod       ,& ! 3*1099.: threshold [hPa]: below this level (p > botmod),
                    !    899.  only obs. increments at obs. levels are used,
                    !          i.e. obs. incr. at model levels are not computed
!   dtqc         ,& ! 720.: timestep (in [s]) for the threshold quality control
    qcc          ,& !  0.,500. constant parts of the quality control thresholds
                    !  0., .7. (='QCT'). (u,v):[m/s], ps: [Pa], T: [k], RH: [ ]
    qcvf         ,& !  5., 1.: multiplication factor to the vertically varying
                    ! 10., 0.  part of the QCT (as def. in 'data_obs_qc_limits')
    qcciq        ,& !  1.   constant part of QC threshold for IWV
    qcsiq        ,& !  .15  IWV QC threshold, as a fraction of IWV of saturated
                    !       model temperature profile
    doromx       ,& !  100., : cut-off and gaussian radius of height differences
                    !  150.,   between model orography and station height for a
                    !  150.,   factor contributing to the quality weight factor
                    !  150.,   as part of the nudging weights
    maxmlo       ,& !        : max. number of multi-level reports in the ODR
    maxgpo       ,& !        : max. number of GPS reports in ODR on total domain
    maxmlv       ,& !  100   : max. number of obs levels in multi-level ODR
    lgps         ,& ! .f.    : .t. if GPS   data is used
    ionl         ,& !  167   : / grid point coordinates
    jonl            !  103   : \ for standard output on nudging

USE data_nudge_all , ONLY :   &

! 5. Miscellany
! -------------

    a_u          ,& ! zonal wind speed            on Arakawa A grid ( m/s )
    a_v          ,& ! meridional wind speed       on Arakawa A grid ( m/s )
    a_p          ,& ! pressure (full value)       on main levels    ( Pa  )
    a_z          ,& ! geometrical height          of main levels    (  m  )
    tmaxbox      ,& ! max. timestep [s] for computing analysis increments
    cqcwbox      ,& ! (maximum) half interval [s] within which observations are
                    ! quality controlled and written on the feedobs file
    liwvssc      ,& ! .t. : spatial consistency check of IWV performed
    gnudgtv      ,& ! 0*10^-4: nudging coefficient for satellite retrieval [1/s]
    rhfrtv       ,& ! 3* 1., : scaling (reduction) factor of the total 'COSAC'
                    !    0.5   for temperature / humidity satellite retrievals
    topobtv      ,& ! 4*1099.: threshold [hPa]: above this level (p < topobs),
                    !          only obs. increments at model levels are used,
                    !             i.e. obs. incr. at obs. levels are not used
    botmotv         ! 4*1099.: threshold [hPa]: below this level (p > botmod),
                    !          only obs. increments at obs. levels are used,
                    !          i.e. obs. incr. at model levels are not computed

! end of data_nudge_all

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

! 2. scalar variables related to the model environment
! ----------------------------------------------------
    madj_hum   ,& ! if = 1 : adjust observed humidity (by ratio of saturation
                  !          vapour pressure over water to the one over ice,
                  !          to be applied if cloud ice is not a state variable)

! 3. pressure dependent scales and geometry of horizontal correlations
! --------------------------------------------------------------------

    ncolev     ,& ! number of levels in the correlation scale tables
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

    nairep     ,& ! AIREP reports (all aircraft reports)
    ntemp      ,& ! TEMP  reports
    npilot     ,& ! PILOT reports
    nsattv     ,& ! SATEM reports
    ngps       ,& ! GPS   reports
    nmodes     ,& !   mode-s report
    nwp_eu     ,& !   European wind profiler report
    nra_eu     ,& !   European SODAR/RASS report
    nravad     ,& !   Radar VAD wind report
    npr_us        !   US Wind Profiler/RASS report

! end of data_obs_lib_cosmo

!-------------------------------------------------------------------------------

USE data_obs_record , ONLY :   &

! 1.1 ODR report header format
! ----------------------------

    mxrhed       ,& ! header length of multi-level reports
    mxghed       ,& ! header length of GPS reports
    mxthed       ,& ! header length of satellite retrieval reports
    nhilon       ,& ! longitude of observing station
    nhjlat       ,& ! latitude  of observing station
    nhalt        ,& ! station altitude [m]
    nhsurf       ,& ! height of model grid pt. to which obs. is assigned
    nhtime       ,& ! time of observat. in forecast hours
    nhzio        ,& ! latitude  of obs. station (or lowest datum) in g.pt units
    nhzjo        ,& ! longitude of obs. station in grid pt. units
    nhvcbu       ,& ! correction factor to vertical correlation scale for wind
                    ! at base of report
    nhvcbt       ,& ! as 'nhvcbu', but for temperature
    nhvcbq       ,& ! as 'nhvcbu', but for humidity
    nhvctu       ,& ! correction factor to vertical correlation scale for wind
                    ! at top of report
    nhvctt       ,& ! as 'nhvctu', but for temperature
    nhvctq       ,& ! as 'nhvctu', but for humidity
    nhtvip       ,& ! observed multi-lev pressure interpol to lowest model level

!   mxrhdf       ,& ! header length of multi-level reports
    mxghdf       ,& ! header length of GPS reports
!   mxthdf       ,& ! header length of satellite retrieval reports
    nhio         ,& ! (local) x-coord. of grid pt. assigned to obs
    nhjo         ,& ! (local) y-coord. of grid pt. assigned to obs
    nhitot       ,& ! global x-coord. of grid pt. assigned to obs
    nhjtot       ,& ! global y-coord. of grid pt. assigned to obs
    nhobtp       ,& ! observation type
    nhcode       ,& ! code type
    nhschr       ,& ! station characteristics (packed as in VOF)
!   nhflag       ,& ! report flags (obs type, surf., alt., sta ID)
    nhpass       ,& ! flag for report being set to 'passive' (as in VOF)
    nhqcfw       ,& ! status of QC and of writing to feedobs files (+ p-QC flag)
    nhnlev       ,& ! number of obs. levels (for multi-level reports)
    nhvqcf       ,& ! for satellite retrieval: threshold quality control flags
    nhaexi       ,& ! for conventional report: flag for exist. of wind or temp.
    nhuexi       ,& ! flag for existence of wind data      in multi-level report
    nhtexi       ,& ! flag for existence of temperat. data in multi-level report
    nhqexi          ! flag for existence of humidity data  in multi-level report

USE data_obs_record , ONLY :   &

! 1.3 ODR body format
! -------------------

!   maxrsl       ,& ! max. number of levels in multi-level ODR
    maxarl       ,& ! max. number of levels for multi-level aircraft reports
    maxrtv       ,& ! max. number of levels in satellite retrieval  reports
    mxrbdy       ,& ! body length of multi-level reports
    nbtu         ,& ! u wind component [m/s]
    nbtv         ,& ! v wind component [m/s]
    nbtt         ,& ! temperature [K]
    nbtrh        ,& ! relative humidity [/]
    nbtp         ,& ! pressure [Pa]
    nbtz         ,& ! height [m]
    nbtuer       ,& ! error of observed wind component
    nbtter       ,& ! error of observed temperature
    nbtqer       ,& ! error of observed relative humidity
    nbtzer       ,& ! error of observed height
    nbtzio       ,& ! longitude in grid pt. units
    nbtzjo       ,& ! latitude  in grid pt. units
    nbtlop       ,& ! LOG( pressure )
    mxrbdf       ,& ! body length of multi-level reports
    nbtflg       ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbterr       ,& ! pre-processing status flags  (bit p., see below: 'nb?err')
    nbtqcf       ,& ! threshold quality control flags      (see below: 'nb?qcf')
    nbtlsg       ,& ! level id (bit pattern, as in NetCDF statistics file)
    nbtlid       ,& ! level identity          (bit pattern, see below: 'nb?lid')
    mxgbdy       ,& ! body length of GPS reports
    nbgtze       ,& ! error in total zenith delay [mm]
    nbgzpd       ,& ! zenith path delay (total zenith delay)             [mm]
    nbgzwd       ,& ! zenith wet delay [mm]
    nbgiwv       ,& ! integrated water vapour [mm]
    nbgbia       ,& ! bias correction to integrated water vapour [mm]
    nbgiwa       ,& ! adjusted (bias corrected) integrated water vapour [mm]
    nbgp         ,& ! pressure [Pa]
    nbgt         ,& ! temperature [K]
    nbgrh        ,& ! relative humidity [/]
    mxgbdf       ,& ! body length of GPS reports
    nbgflg       ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbgerr       ,& ! pre-processing status flags  (bit p., see below: 'nb?err')
    nbgqcf       ,& ! threshold quality control flags
    nbglid          ! level identity          (bit pattern, see below: 'nb?lid')

USE data_obs_record , ONLY :   &

    mxtbdy       ,& ! body length of multi-level reports
    nbvt         ,& ! temperature [K]
    nbvrh        ,& ! relative humidity [/]
    nbvp         ,& ! pressure [Pa]
    mxtbdf       ,& ! body length of sat retrieval reports
    nbvflg       ,& ! main flag word (bit pattern as in VOF body, word 'nvbmfw')
    nbverr       ,& ! pre-processing status flags  (bit p., see below: 'nb?err')

! 1.4 Bit patterns for packed information in ODR body
! ---------------------------------------------------

    nvfgbp       ,& ! bit pos. for main flag on IWV / ZPD               "
    nvfbps       ,& ! bit pos. for main flags: 4: above 300 hPa         "
    nvfboc       ,& ! no. of bits occ. for main flags                   "

!       1.4.2  Other bit patt. for packed info in ODR (VOF) body, general words
!              ----------------------------------------------------------------
    nvru         ,& ! bit pos. for status/QC flags for horiz. wind nb?err/nb?qcf
    nvrt         ,& ! bit pos. for status/QC flags for temperature       "
    nvrq         ,& ! bit pos. for status/QC flags for humidity          "
    nvrz         ,& ! bit pos. for status/QC flags for pressure/height   "
    nvriwv       ,& ! bit pos. for status/QC flags for IWV               "
    nvrzpd       ,& ! bit pos. for status/QC flags for ZPD               "

! 1.5 Further quantities related to ODR
! -------------------------------------

    imdi         ,& ! missing data indicator for ODR integers (2^31 -1)
    ilstid       ,& ! character length of the station identity
    ilstidp      ,& ! char. length used for printing the station ID
                    ! Note: (ilstid >= ilstidg >= ilstidp) cf. data_nudge_gather
    ystid        ,& ! obs. station identity to be printed
    fdoro           ! scaling factor to vertical distances betw. model orography
                    ! and station height for (z-obs < z-mod)

USE data_obs_record , ONLY :   &

! 2. Observation data records (ODR)
! ---------------------------------

    omlbdy       ,& ! body   of multi-level ODR
    omlhed       ,& ! header of multi-level ODR
    ogpbdy       ,& ! body   of GPS ODR
    ogphed       ,& ! header of GPS ODR
    otvbdy       ,& ! body   of satellite retrieval ODR
    otvhed       ,& ! header of satellite retrieval ODR
    momlbd       ,& ! body   of multi-level ODR
    momlhd       ,& ! header of multi-level ODR
    mogpbd       ,& ! body   of GPS ODR
    mogphd       ,& ! header of GPS ODR
    motvbd       ,& ! body   of satellite retrieval ODR
    motvhd       ,& ! header of satellite retrieval ODR
    yomlhd       ,& ! header of multi-level ODR
    yogphd       ,& ! header of GPS ODR
    yotvhd       ,& ! header of satellite retrieval ODR

! 3. Simulated Observation Record (SOR)
! -------------------------------------

    mxsoml       ,& ! SOR body length for multi-level reports
    nso_iq       ,& ! integrated water vapour (increment)            [mm]
    nso_zpd      ,& ! zenith total path delay                        [mm]
    smlbdy       ,& ! body of multi-level SOR
    sgpbdy       ,& ! body of GPS (IWV) SOR
    stvbdy          ! body of satellite retrieval SOR

! end of data_obs_record

!-------------------------------------------------------------------------------

USE data_obs_qc_limits , ONLY :   &

! Table for pressure dependent quality control thresholds
! -------------------------------------------------------

    qcvsond      ,& ! (root of) radiosonde (TEMP, PILOT) wind error variance
    qctsond      ,& ! (root of) radiosonde temperature error variance
    qctf         ,& ! temporal factor for thresholds of upper-air obs. 
    qcftsiv      ,& ! reduction factor to the threshold used for IWV-SSC
    qcfbias      ,& ! fraction of bias added to spatial consistency threshold
    dtchkps         ! spatial check: radius of infl. for linear temporal weights

! end of data_obs_qc_limits

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
    kobtyml      ,& ! CMA observation type
    kcdtyml      ,& ! CMA observation code type
    mszlev       ,& ! number of vertical levels with obs. incr. in the OIR
    mbotlv       ,& ! number of model levels above the surface without
                    ! observation increment
    mtoplv       ,& ! number of model levels at the top of the model without
                    ! observation increment
    ltiml        ,& ! .TRUE if temporal linear interpol. for multi-level rep.
    ystidml         ! station identity of multi-level station

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
    noiu         ,& ! zonal wind observation increment
    noiv         ,& ! meridional wind observation increment
    noit         ,& ! temperature observation increment
    noiqd        ,& ! specific humidity observation increment
    noirh        ,& ! relative humidity observation increment
    noiz         ,& ! height (of obs. increment level)
    noith        ,& ! potential temperature
    noilp        ,& ! log( pressure )
    noiuqc       ,& ! quality weight to wind obs. increment
    noitqc       ,& ! quality weight to temperature obs. increment
    noiqqc       ,& ! quality weight to humidity obs. increment
    noivcr       ,& ! normalization factor to vertical interpolation distances
    noiulr       ,& ! vertical interpolat. distance to obs. incr. level for wind
    noitlr       ,& ! vertical int'pol distance to obs. incr. level for temper.
    noiqlr       ,& ! vertical int'pol distance to obs. incr. level for humidity

! 4. (Local) information gathered by 1 (or 2) nodes for printing for control
! --------------------------------------------------------------------------

    oyqc         ,& ! info. on data rejected by the threshold quality control
    myqc         ,& ! info. on data rejected by the threshold quality control
    yyqc         ,& ! info. on data rejected by the threshold quality control

! 5. Variables defining the size of the arrays containing the local information
! -----------------------------------------------------------------------------

    maxivq       ,& ! max.  number of IWV reports used for spatial check
    maxqcp       ,& ! max.  number of rejected data to be printed per timestep
    nmloit       ,& ! total number of active profiles of obs. incr. in the OIR
    nmltot       ,& ! total number of active multi-level stations
    nivtot       ,& ! total number of IWV reports used for spatial check
    ntotqc       ,& ! total number of rejected data to be printed per timestep
    lnisua       ,& ! non-isotrophic correlat. for upper-air single-lev. data

! 6. Geometrics and variables used for spatial consistency check of pressure obs
! ------------------------------------------------------------------------------

    pyyd         ,& ! latitudinal (merid.) distance on tangent cone projection
    pxxd2        ,& ! square of zonal distance on tangent cone projection
    pxsalpa      ,& ! factor used for distances on tangent cone projection
    isrtvqc         ! (sorted) list of stations with IWV

! end of data_nudge_gather

!-------------------------------------------------------------------------------

USE data_nudge_local , ONLY :   &

! General parameters and other variables
! --------------------------------------

    ptropop      ,& ! pressure at tropopause (standard atmosphere)
    thdzt        ,& ! mean vertical potential temperat. gradient in troposphere
    fsvcut       ,& ! fraction with which the adjusted scale S_adj, rather than
                    ! the scale S_nl as specified in the namelist, is used for
                    ! determin. the vertical cut-off in log(p) or theta(e)-units
    wtml         ,& ! temporal nudging weights for multi-level reports
    wtgp         ,& ! temporal nudging weights for GPS reports
    wttv         ,& ! temporal nudging weights for satellite retrievals
    tmladm       ,& ! administrator of multi-level ODR
!   tgpadm       ,& ! administrator of GPS ODR
!   ttvadm       ,& ! administrator of satellite retrievals
    zsdni        ,& ! equidistant pts. used for non-isotropic horiz. correlation
    imladm       ,& ! administrator of multi-level ODR
    igpadm       ,& ! administrator of GPS ODR
    itvadm       ,& ! administrator of satellite retrievals
    nmlsta       ,& ! total number of sorted multi-level stations
!   nmasta       ,& ! total number of sorted multi-level aircraft reports
    ngpsta       ,& ! total number of sorted GPS stations
    ntvsta       ,& ! total number of sorted satellite retrieval 'stations'
    kml300       ,& ! index of full model level corresponding to about 300 hPa
    lqcall          ! .TRUE if threshold quality control for all observations
                    !          at current timestep

! end of data_nudge_local

!-------------------------------------------------------------------------------

USE mo_fdbk_tables, ONLY :   &

    ot_airep     ,& ! observation type for aircraft obs
    ot_temp      ,& ! observation type for TEMP radiosonde
    ot_pilot     ,& ! observation type for PILOT (incl. wind profiler, RASS,...)
    ot_satem        ! observation type for satellite retrievals

! end of mo_fdbk_tables

!-------------------------------------------------------------------------------

USE src_obs_operator_conv, ONLY :   &

    q2rh_col            ,& ! computes relative humidity from specific humidity
    tvirt_col           ,& ! computes model columns of virtual temperature
    hhl_col             ,& ! computes height on model half levels
    lhyphl_col          ,& ! computes LOG of hydrostatic pressure on half levels
    prep_vi_mo2ob       ,& ! prepares interpol. from model levels to obs levels
    mult_obs_operator   ,& ! forward observation operator for multi-level obs
    mult_obs_2_modlev   ,& ! inverse observation operator for multi-level obs
    mult_obs_operator_z ,& ! supplementary obs operator for multi-level T, z
    zpd_iwv_obs_operator   ! observation operator for GPS ZPD (and IWV) obs

! end of src_obs_operator_conv

!-------------------------------------------------------------------------------

USE src_obs_qc_conv      , ONLY :   &

    mult_obs_qc_fg      ,& ! quality control for individual observation
    mult_obs_qc_dz      ,& ! height / thickness QC for multi-level temperature
    iwv_quality_cntl       ! quality control of individual IWV observations

! end of src_obs_qc_conv

!-------------------------------------------------------------------------------

USE src_tracer       , ONLY : trcr_get, trcr_errorstr

!-------------------------------------------------------------------------------
! Local declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE

! 2. Parameters
! -------------

  INTEGER (KIND=iintegers) , PARAMETER ::  &
    mxvimo  = 5  ,& ! number of variables in 'vimtoob'
    mxviom  = 4     ! number of variables in 'viobtom'

! 2. Variables
! ------------

  INTEGER (KIND=iintegers)                , PRIVATE :: &
    io     ,jo     ,& ! local  indices of location of observation
    io_tot ,jo_tot ,& ! global indices of location of observation
    ista           ,& ! index of observing station
    itim           ,& ! (time) index over obs. at one observing station
    nlev           ,& ! number of observation levels in the multi-level report
    kobtyp            ! observation type

  LOGICAL                                 , PRIVATE :: &
    levdiff           ! obs. increments at obs. levels of multi-level reports
                      ! are used for spreading

  REAL (KIND=wp),     POINTER :: &
    qv  (:,:,:) => NULL() , &    ! QV at tlev=nnew
    qc  (:,:,:) => NULL()        ! QC at tlev=nnew

!===============================================================================

!-------------------------------------------------------------------------------
! Public and Private Subroutines
!-------------------------------------------------------------------------------

CONTAINS


!-------------------------------------------------------------------------------
!+ Module procedure for building vertical profiles of observation increments
!-------------------------------------------------------------------------------

SUBROUTINE mult_obs_increment ( icml, ktopoi, kcdtyp, ke, col_t, col_qv, col_p &
                              , col_z, col_qc, col_qrs, col_th                 &
                              , nlev, zobbdy, mzobbd, vim2ob, viob2m, lobinc   &
                              , numlev )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure builds vertical profiles of observation increments
!   for multi-level reports after all the vertical interpolation has already
!   been done. It also supplies some additional quantities at the top and the
!   base of the profile.
!
! Method:
!   Merging of the difference between the vertically interpolated observation
!   profile and the model value profile with the difference between the obser-
!   vation profile and the vertically interpolated model value profile into
!   one vertical profile of observation increments (differences).
!   For the computation of observation increments of 'specific humidity'
!   (used for horizontal spreading or spreading along isentropic surfaces), 
!   'observed relative humidity times model saturation vapour pressure' is
!   taken as observed vapour pressure instead of the (more) directly observed
!   vapour pressure. In this way, it is achieved that the model's relative
!   humidity is nudged (at the obs. location) towards the observed relative
!   humidity (irrespective of the temperature difference between model and
!   observation). In this sense, the resulting specific humidity increment is
!   equivalent to the relative humidity increment.
!
! Written by        :  Christoph Schraff, DWD  (original version: 02.10.97)
! Current Code Owner:  Christoph Schraff, DWD
!
! Modification comment:
!   Date:   13.08.98   Bug corrected: if no obs. increments at model levels
!                      existed, then the increments at observation levels would
!                      have been omitted erroneously.    Christoph Schraff
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

  INTEGER (KIND=iintegers), INTENT (IN)         ::       &
    icml             ,& ! index of obs. sta. in data record to be used further
    ktopoi           ,& ! number vertical levels used in the observation incre-
                        ! ment array 'oiml' prior to the current report
    kcdtyp           ,& ! CMA observation code type
    nlev             ,& ! number of vertical levels in observation report
    ke                  ! number of vertical (main) levels in model column

  REAL    (KIND=wp   ),     INTENT (IN)         ::       &
                        ! model column (at the observation location) of:
    col_t     (ke)   ,& ! temperature                                   (  K  )
    col_th    (ke)   ,& ! potential temperature                         (  K  )
    col_qv    (ke)   ,& ! specific water vapor content                  (kg/kg)
    col_qc    (ke)   ,& ! specific cloud water content                  (kg/kg)
    col_qrs   (ke)   ,& ! spec. cont. of hydrometeors excl. cloud water (kg/kg)
    col_p     (ke)   ,& ! pressure (full value)       on main levels    ( Pa  )
    col_z     (ke)      ! geometrical height          of main levels    (  m  )

  REAL    (KIND=wp   ),     INTENT (IN)         ::       &
    zobbdy (nlev,mxrbdy) ,& ! multi-level obs. body (format: see 'omlbdy')
    vim2ob (nlev,mxvimo) ,& ! model values interpolated to obs. levels
    viob2m (ke  ,mxviom)    ! obs. (incr.) interpolated to model levels
                            !   (wind: increments; T, RH: full obs values)

  INTEGER (KIND=iintegers), INTENT (IN)         ::       &
    mzobbd (nlev,mxrbdf)   ! multi-level obs. body (format: see 'momlbd')

  LOGICAL                 , INTENT (IN)         ::       &
    lobinc    (3)       ! the quantities to be interpolated are observation
                        ! increments rather than observed values

  INTEGER (KIND=iintegers), INTENT (OUT)        ::       &
    numlev              ! number of levels in the obs. increment profile

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    modespr          ,& ! mode of spreading (horizontal, along model levels
                        !                   ,along isentropes)
    kbotlev          ,& ! number of obs levels below lowest    model level
    ktoplev          ,& ! number of obs levels above uppermost model level
    ivrs             ,& ! index of observed quantity: 1=(u,v); 2=T; 3=RH
    ivar             ,& ! index of observed quantity: 1=(u,v); 3=T; 4=RH
!   nmlob            ,& ! index of report in the ODR (obs data record)
!   nmloi            ,& ! index of report in the obs. increment record
    km               ,& ! vertical (loop) index of model level
    klev             ,& ! vertical (loop) index of obs. increment level
    kobs             ,& ! vertical (loop) index of observation level
    klvi             ,& ! vertical index in 'oiml'
    nobs             ,& ! index of top obs. level with interpolated model values
    mdabot           ,& ! index of lowest model level with interpolated obs.
    mdabotl          ,& ! index of lowest model level possibly with obs. below
    mdatop           ,& ! index of uppermost model level with interpolated obs.
    mdatopo          ,& ! index of uppermost model level with obs. just below
    klevmax          ,& ! max. possible level number in current obs incr.profile
    nbotop , nba     ,& ! loop indices
    kboti            ,& ! first guess for lower index limit of unstable region
    ksbot  , kstop   ,& ! index range in obs. incr. record of unstable region
    kll              ,& ! vertical loop index over obs. increment levels
    noix             ,& ! index of variable in obs. increment record
    nois             ,& ! index of spreading parameter in obs. increment record
    nvrx             ,& ! bit pos. for status 'active' for variable
    ivio             ,& ! index of interpolated quantity: 3=T; 4=RH
    ilvip (ke,3,2)   ,& ! indices of obs levels adjacent to target model level
    ka     , kb      ,& ! indices of model levels adjacent to target obs level
    ilva   , ilvb    ,& ! indices of obs levels adjacent to target model level
    ilvsta              ! (lower) start index for loop over obs levels

  REAL    (KIND=wp   )     ::  &
    zftv2t_m (ke)    ,& ! inverse obs operator for virtual temp. (model levels)
    zftv2t_o         ,& ! inverse obs operator for virtual temp. (obs. levels)
    zpbot            ,& ! pressure of level below which obs increments at model
                        !   levels are not used for nudging
    zstarg           ,& ! value of spreading parameter at target level
    zdvia  , zdvib   ,& ! vertical distance of interpol.levels from target level
    ztp              ,& ! value of spreading parameter at standard tropopause
    zqts             ,& ! quotient of the tropospheric to the stratospheric mean
                        ! vertical gradient of the spreading parameter
    zpvmod           ,& ! model      water vapour pressure
    zpvsat           ,& ! saturation water vapour pressure
    zrhobs           ,& ! observed relative humidity
    zrhmod           ,& ! model    relative humidity
    zqvobs           ,& ! observed specific humidity
    zqvmov           ,& ! model    specific humidity
    qcupr  , qctpr   ,& ! thresholds for obs. increment beyond which they may be
             qcqpr   ,& !  \ printed for control
    zeklop0          ,& ! scaling factor for converting temperature to pot. temp
    sprmax           ,& ! max. pot. temperature within unstable region
    sprmin           ,& ! min. pot. temperature within unstable region
    sprmxi           ,& ! max. pot. temperature within & above unstable region
    zloprev          ,& ! log( pressure ) of previous obs increment level
    zlopf               ! weight factor for vert. interpol. to obs. incr. level

  LOGICAL                  ::  &
    lwrite           ,& ! printout for control at present station / node / time
    lwruv  , lwrt    ,& ! obs. increment of wind , temperature , or humidity may
             lwrq    ,& !  \ be printed for control
    lextu  , lextt   ,& ! obs. increment of wind , temperature , or humidity is
             lextq   ,& !  \ to be determined on current level if possible
    lexisto          ,& ! at least 1 obs. increment exists on current level
    lmono            ,& ! vertical profile stable (in terms of pot. temperature)
    ltop                ! top of (model or obs.) profile reached


! Local (automatic) arrays:
! -------------------------

  INTEGER (KIND=iintegers) ::  &
    kbz     (nlev)          ! level indices for vertical interpol. to obs levels

  REAL    (KIND=wp   )     ::  &
    zcol_lnp  (ke)       ,& ! log( pressure )
    zlpf    (nlev)       ,& ! weight factors for vert. interpolat. to obs levels
    zinv_ob (nlev,mxrbdy),& ! inverse obs operator applied to observed values
    z_spr_o (nlev)       ,& ! spreading parameter at obs level
    zuz2p_o (nlev)       ,& ! convertor from height to pressure units at obs lev
    zuz2por (nlev)       ,& ! convertor from pressure to height units at obs lev
    zwvip_o (nlev)       ,& ! inverse combined interpol. distance to obs level
    zwvip_m (ke  ,3)     ,& ! inverse combined interpol. distance to model level
    zpresoi (maxmlv+ke)  ,& ! pressure                  on obs. increment levels
    zemkpoi (maxmlv+ke)  ,& ! part. factor for convert. temperature to pot. temp
    zrtdgoi (maxmlv+ke)  ,& ! R/g * virtual temperature on obs. increment levels
    zpothoi (maxmlv+ke)  ,& ! pot. temperature (uncorrected in unstable regions)
    vimtoob (nlev,mxvimo),& ! model values interpolated to obs. levels
    viobtom (ke  ,mxviom)   ! obs. (incr.) interpolated to model levels
                            !   (wind: increments; T, RH: full obs values)
!
!------------ End of header ----------------------------------------------------

 
!-------------------------------------------------------------------------------
! Begin Subroutine mult_obs_increment
!-------------------------------------------------------------------------------
  
!-------------------------------------------------------------------------------
!  Section 1: Preliminaries
!-------------------------------------------------------------------------------

! nmloi = knoiml(icml,itim)
! nmlob = imladm(ista,itim)

  vimtoob (:,:)  =  vim2ob(:,:)
  viobtom (:,:)  =  viob2m(:,:)

  modespr = msprpar

  zeklop0 = EXP( rdocp * LOG( p0r ) )

! thresholds for control output
! (qcvsond(4 --> 500hPa) = 3.4m/s ;   qctsond(3 --> 700hPa) = 0.7K)

  qcupr  = 0.85_wp * (qcc(1) + qcvf(1) *qcvsond(4))
  qctpr  = 0.85_wp * (qcc(3) + qcvf(3) *qctsond(3))
  qcqpr  = 0.85_wp *  qcc(4)

  lwrite = (MOD( ntstep,MAX( 1, 6/NINT(tconbox/dt) ) ) == 0)
  lwruv  = (lwrite) .AND. (zwtml(icml,itim,1) > 0.8_wp) .AND. (lwonl)
  lwrt   = (lwrite) .AND. (zwtml(icml,itim,2) > 0.8_wp) .AND. (lwonl)
  lwrq   = (lwrite) .AND. (zwtml(icml,itim,3) > 0.8_wp) .AND. (lwonl)

  lwrite = (lfirst) .AND. (io == ionl) .AND. (jo == jonl) .AND. (lwonl)

  DO km = 1 , ke
    zcol_lnp (km) = LOG( col_p(km) )
  ENDDO

  ! prepare vertical interpolation from model level to obs level
  ! by computing 'kbotlev', 'ktoplev', 'kbz', and 'zlpf'

  CALL prep_vi_mo2ob ( ke, col_p, zcol_lnp, nlev, zobbdy                       &
                     , kbotlev, ktoplev, kbz, zlpf )
! ==================

! adjust range of levels where observation increments are used for nudging

  DO ivrs = 1 , 3
    ivar = ivrs + MIN( ivrs-1 , 1 )
    ! pressure level below which obs increments at model levels are not used
    IF (ksprml(icml) == ke) THEN
      zpbot  =  109900.0_wp
    ELSE
      IF (kobtyp /= nsattv) zpbot = MAX( botmod (ivar) , c1 )
      IF (kobtyp == nsattv) zpbot = MAX( botmotv(ivar) , c1 )
      ! ksprml: level at and above which spreading is along model levels
      zpbot  =  MAX( zpbot , col_p(ksprml(icml)) )
    ENDIF
    ! adjust 'mbotlv'
!   mbot = mbotlv(icml,itim,ivrs)
    DO km = ke , 1 , -1
      IF (col_p(km) > zpbot+epsy)                                              &
!       mbot = MAX( mbot , ke-km+1 )
        mbotlv (icml,itim,ivrs) = MAX( mbotlv(icml,itim,ivrs) , ke-km+1 )
!     IF ((col_p(km) > zpbot+epsy) .AND. (mbotlv(icml,itim,ivrs) < ke-km+1))   &
!       WRITE(0,*) 'MBOT ',ystid, ntstep, km, ivrs                             &
!                         ,mbotlv(icml,itim,ivrs), mtoplv(icml,itim,ivrs)      &
!                         ,col_p(km), zpbot
    ENDDO
  ENDDO
! IF (ntstep <  9) WRITE(0,*) 'MBOT7 ',ystid, ntstep, mbotlv(icml,itim,:)


!-------------------------------------------------------------------------------
!  Section 2: Apply the inverse observation operator
!             - to the original and the interpolated observed values
!             - to the interpolated model value to which the forward observation
!               operator has previously been applied for quality control.
!
!             Currently, the only inverse observation operator applied
!             is the conversion of RASS virtual temperature into temperature.
!-------------------------------------------------------------------------------

  DO km = 1 , ke
    IF (     (kcdtyp == nra_eu) .OR. (kcdtyp == nwp_eu)                        &
        .OR. (kcdtyp == nravad) .OR. (kcdtyp == npr_us)) THEN
      zftv2t_m (km) = c1 / (c1 + rvd_m_o *col_qv(km) - col_qc(km) - col_qrs(km))
      IF (viobtom(km,3) > rmdich)                                              &
      ! apply inverse observation operator to interpolated observation value
        viobtom (km,3)  =  viobtom(km,3) * zftv2t_m(km)
    ENDIF
  ENDDO

  DO kobs = 1 , nlev
    zinv_ob (kobs,nbtu ) = zobbdy (kobs,nbtu )
    zinv_ob (kobs,nbtv ) = zobbdy (kobs,nbtv )
    zinv_ob (kobs,nbtt ) = zobbdy (kobs,nbtt )
    zinv_ob (kobs,nbtrh) = zobbdy (kobs,nbtrh)
    IF (     (kcdtyp == nra_eu) .OR. (kcdtyp == nwp_eu)                        &
        .OR. (kcdtyp == nravad) .OR. (kcdtyp == npr_us)) THEN
      zftv2t_o =  (c1 -zlpf(kobs)) * zftv2t_m(kbz(kobs)  )                     &
                 +     zlpf(kobs)  * zftv2t_m(kbz(kobs)-1)
      IF (zinv_ob (kobs,nbtt) > rmdich)                                        &
      ! apply inverse observation operator to original observation
        zinv_ob (kobs,nbtt)  =  zinv_ob(kobs,nbtt) * zftv2t_o
      IF (vimtoob(kobs,3)     > rmdich)                                        &
      ! apply inverse observation operator to interpolated model value
      ! to which the forward observation operator has been previously applied
      ! for quality control
        vimtoob (kobs,3)     =  vimtoob(kobs,3)    * zftv2t_o
    ENDIF
  ENDDO


!-------------------------------------------------------------------------------
!  Section 3: For quality weighting depending on interpolation range:
!             Measure the 'combined distance' of the interpolation points
!             to the target level
!            (The 'combined distance' is computed from the 2 single distances in
!             the same way as the total resistance from 2 parallel resistances.)
!-------------------------------------------------------------------------------

  ztp = zstpml(icml)
  IF (modespr <= 1) zqts = c1
  IF (modespr == 2) zqts = c1 / qst
  zwvip_o  =  rmdi
  zwvip_m  =  rmdi

! on observation levels: measure the 'combined distance'
!                        of the adjacent model levels to the obs level

  CALL prep_vi_mo2ob ( ke, col_p, zcol_lnp, nlev, zobbdy(1:nlev,:)             &
                     , kbotlev, ktoplev, kbz, zlpf )
! ==================

  DO kobs = 1 , nlev
    kb  =  kbz(kobs)
    ka  =  kbz(kobs) - 1
    z_spr_o (kobs)   =   (c1-zlpf(kobs)) *zsprml (icml,kb)                     &
                        +    zlpf(kobs)  *zsprml (icml,ka)
    zuz2p_o (kobs)   =   (c1-zlpf(kobs)) *zrtgkml(icml,kb)                     &
                        +    zlpf(kobs)  *zrtgkml(icml,ka)
    IF (modespr <= 1)  zuz2por (kobs) =  c1 / zuz2p_o(kobs)
  ENDDO

  DO kobs = 1 + kbotlev , nlev - ktoplev
    kb  =  kbz(kobs)
    ka  =  kbz(kobs) - 1
    IF (modespr <= 1) THEN
! take  ABS( zdvi? )  since  (cdpext > 0)    (==> approximation)
      zdvia = ABS(   zsprml(icml,ka) - z_spr_o(kobs) ) / zrtgkml(icml,ka)
      zdvib = ABS( - zsprml(icml,kb) + z_spr_o(kobs) ) / zrtgkml(icml,kb)
    ELSE
      zdvia =  MAX( MIN( zsprml(icml,ka),ztp ) - z_spr_o(kobs) , c0 )          &
             + MAX( zsprml(icml,ka) - MAX( z_spr_o(kobs),ztp ) , c0 ) *zqts
      zdvib =  MAX( MIN( z_spr_o(kobs),ztp ) - zsprml(icml,kb) , c0 )          &
             + MAX( z_spr_o(kobs) - MAX( zsprml(icml,kb),ztp ) , c0 ) *zqts
    ENDIF
    zwvip_o (kobs)   = zdvia * zdvib / MAX( zdvia + zdvib , epsy )
  ENDDO

! on model levels: measure the 'combined distance'
!                  of the adjacent obs levels to the model level

  DO ivrs = 1 , 3
    ivio = ivrs + MIN( ivrs-1 , 1 )
    IF (ivrs == 1) nvrx = nvru
    IF (ivrs == 2) nvrx = nvrt
    IF (ivrs == 3) nvrx = nvrq

    ! Get indices of obs levels adjacent to target model level
!   ilvsta = 1 + kbotlev
    ilvsta = 1
    DO klev = 1+ mbotlv(icml,itim,ivrs) , ke- mtoplv(icml,itim,ivrs)
      km   = ke + 1 - klev
      ilvip (km,ivrs,1) = 0
      ilvip (km,ivrs,2) = 0
      zstarg  =  zsprml(icml,km)
      ilvb = ilvsta
      loop_over_obs_levels: DO ilva = ilvsta+1, nlev
        IF (      (      BTEST( mzobbd(ilva,nbterr),nvrx ))                    &
            .AND. (.NOT. BTEST( mzobbd(ilva,nbtqcf),nvrx ))) THEN
          IF (      (z_spr_o(ilva) >= zstarg)                                  &
              .AND. (z_spr_o(ilvb) <= zstarg)) THEN
            ilvip (km,ivrs,1) = ilvb
            ilvip (km,ivrs,2) = ilva
                                                       EXIT loop_over_obs_levels
          ENDIF
          ilvb = ilva
        ENDIF
      ENDDO loop_over_obs_levels
      ilvsta = ilvb
!     IF (ilva > nlev-ktoplev)  ilva = 0
!     IF ((     (ilvip2(km,ivrs,1) /= ilvip(km,ivrs,1))                         &
!          .OR. (ilvip2(km,ivrs,2) /= ilvip(km,ivrs,2))) .AND. (ntstep < 99))   &
!       WRITE(0,*) 'ZLEV ',ystid, ntstep, km,ivrs, ilvip2(km,ivrs,1), ilvip2(km,ivrs,2) &
!                                        , ilvip (km,ivrs,1), ilvip (km,ivrs,2) &
!                         ,mbotlv(icml,itim,ivrs), mtoplv(icml,itim,ivrs)
    ENDDO

    DO klev = 1+ mbotlv(icml,itim,ivrs) , ke- mtoplv(icml,itim,ivrs)
      km   = ke + 1 - klev
      ilvb = ilvip(km,ivrs,1) 
      ilva = ilvip(km,ivrs,2) 
      zstarg  =  zsprml(icml,km)
      IF ((ilva > 0) .AND. (ilvb > 0) .AND. (viobtom(km,ivio) > rmdich)) THEN
        IF (modespr <= 1) THEN
          zdvia = (  z_spr_o(ilva) - zstarg) * zuz2por(ilva)
          zdvib = (- z_spr_o(ilvb) + zstarg) * zuz2por(ilvb)
        ELSE
          zdvia =  MAX( MIN( z_spr_o(ilva),ztp ) - zstarg , c0 )               &
                 + MAX( z_spr_o(ilva) - MAX( zstarg,ztp ) , c0 ) *zqts
          zdvib =  MAX( MIN( zstarg,ztp ) - z_spr_o(ilvb) , c0 )               &
                 + MAX( zstarg - MAX( z_spr_o(ilvb),ztp ) , c0 ) *zqts
        ENDIF
        IF (lobinc(ivrs))  zdvia = zdvia + zwvip_o(ilva)
        IF (lobinc(ivrs))  zdvib = zdvib + zwvip_o(ilvb)
        zwvip_m (km,ivrs) = zdvia * zdvib / MAX( zdvia + zdvib , epsy )
      ELSEIF (viobtom(km,ivio) < rmdich) THEN
        zwvip_m (km,ivrs) = c0
      ENDIF
    ENDDO
  ENDDO


!-------------------------------------------------------------------------------
!  Section 4: Profile with observation increments only on main model levels
!-------------------------------------------------------------------------------

! check if obs. increments at obs. levels will be used for spreading
  levdiff = .TRUE.
! levdiff =       (modespr > 0) .AND. (.NOT. lretv)                            &
!           .AND. (MIN( topobs(1) , topobs(3) , topobs(4) ) < 1049.9_wp)   &
!      .OR.       (modespr > 0) .AND. (      lretv)                            &
!           .AND. (MIN( topobtv(1), topobtv(3), topobtv(4)) < 1049.9_wp)
! IF (lwrite) WRITE( nupr,'(''mult_obs_increment: levdiff '',L1)' ) levdiff

  IF (.NOT. levdiff) THEN

! preset 'missing values' and determine pressure, height, and potential temper.

    DO klev = 1 , ke
      klvi = ktopoi + klev
      km   = ke + 1 - klev
      oiml (noiu  ,klvi) = rmdi
      oiml (noiv  ,klvi) = rmdi
      oiml (noit  ,klvi) = rmdi
      oiml (noiqd ,klvi) = rmdi
      oiml (noirh ,klvi) = rmdi
      oiml (noiz  ,klvi) = col_z (km)
      oiml (noith ,klvi) = col_th(km)
      oiml (noiuqc,klvi) = c1
      oiml (noitqc,klvi) = c1
      oiml (noiqqc,klvi) = c1
      oiml (noilp ,klvi) = zlopml(icml,km)
      oiml (noiulr,klvi) = c0
      oiml (noitlr,klvi) = c0
      oiml (noiqlr,klvi) = c0
      oiml (noivcr,klvi) = c1 / zrtgkml(icml,km)
      IF (modespr == 2) oiml (noivcr,klvi) = c1
      zpresoi    (klev)        = col_p(km)
    ENDDO
    IF (lwrite) WRITE( nupr,'(''mult_obs_increment: before u,v'')' )

! get obs. increments of horizontal wind

    DO klev = 1+ mbotlv(icml,itim,1) , ke- mtoplv(icml,itim,1)
      klvi = ktopoi + klev
      km   = ke + 1 - klev
      oiml (noiu  ,klvi) = viobtom(km,1)
      oiml (noiv  ,klvi) = viobtom(km,2)
      oiml (noiulr,klvi) = zwvip_m(km,1)
      IF ((lwruv) .AND. (     (oiml(noiu,klvi) > qcupr)                        &
                         .OR. (oiml(noiv,klvi) > qcupr)))                      &
        WRITE( nupr,'(''Du,Dv '',A ,I5,2F8.2,F8.0)' ) ystid, klev              &
             , oiml(noiu,klvi), oiml(noiv,klvi), col_p(km)
    ENDDO
    IF (lwrite) WRITE( nupr,'(''mult_obs_increment: before T'')' )

! get obs. increments of temperature

    DO klev = 1+ mbotlv(icml,itim,2) , ke- mtoplv(icml,itim,2)
      klvi = ktopoi + klev
      km   = ke + 1 - klev
      oiml (noit  ,klvi) = viobtom(km,3) - col_t(km)
      oiml (noitlr,klvi) = zwvip_m(km,2)
      IF ((lwrt ) .AND. (oiml(noit,klvi) > qctpr))                             &
        WRITE( nupr,'(''DT    '',A ,I5, F8.2,F8.0)' ) ystid, klev              &
             , oiml(noit,klvi), col_p(km)
    ENDDO

! get obs. increments of specific & relative humidity

    DO klev = 1+ mbotlv(icml,itim,3) , ke- mtoplv(icml,itim,3)
      klvi = ktopoi + klev
      km   = ke + 1 - klev
      zpvmod = fq2pv ( col_qv(km) , col_p(km) , rdv )
      zpvsat = fpvsw ( col_t (km) , b1, b2w, b3, b4w )
      oiml (noirh ,klvi) = viobtom(km,4) - MIN( zpvmod /zpvsat , c1 )
      oiml (noiqd ,klvi) =   fpv2q( viobtom(km,4)*zpvsat , col_p(km) , rdv )   &
                           - col_qv(km)
      oiml (noiqlr,klvi) = zwvip_m(km,3)
      IF ((lwrq ) .AND. (oiml(noirh,klvi) > qcqpr))                            &
        WRITE( nupr,'(''DRH   '',A ,I5, F8.4,F8.0)' ) ystid, klev              &
             , oiml(noirh,klvi), col_p(km)
    ENDDO

! get number of levels in the obs. increment profile

    numlev = ke

! get additional local information for later spreading

    DO nbotop = 1 , 2
      nba = 2*(itim-1) + nbotop
      DO ivrs = 1 , 3
        IF (nbotop == 1) km = ke - mbotlv(icml,itim,ivrs)
        IF (nbotop == 2) km = 1  + mtoplv(icml,itim,ivrs)
        zspobml(icml,nba,ivrs) = zsprml(icml,km)
        zpobml (icml,nba,ivrs) = col_p(km)
        zemkpml(icml,nba,ivrs) = EXP( - rdocp * zlopml(icml,km) )
        zrtdgml(icml,nba,ivrs) = zrtgkml(icml,km)
      ENDDO
    ENDDO

! printout for control

    IF (lwrite) THEN
      DO klev = 1 , numlev
        klvi = ktopoi + klev
        km   = ke + 1 - klev
        WRITE( nupr,'(3F8.3, F10.5,F8.3, 2F8.0,F7.2)' )                        &
               oiml(noiu ,klvi), oiml(noiv ,klvi)                              &
             , oiml(noit ,klvi), oiml(noiqd,klvi)                              &
             , oiml(noirh,klvi), col_p(km)                                     &
             , oiml(noiz ,klvi), oiml(noith,klvi)
      ENDDO
    ENDIF

  ELSE


!-------------------------------------------------------------------------------
!  Section 3: Profile with observat. increments on model and observation levels
!-------------------------------------------------------------------------------

! initialize indices

    kobs    = 1 + kbotlev
    nobs    = nlev - ktoplev
    ltop    = .FALSE.
    mdabot  = ke - MIN( mbotlv(icml,itim,1) , mbotlv(icml,itim,2)              &
                      , mbotlv(icml,itim,3) )
    mdatop  = 1  + MIN( mtoplv(icml,itim,1) , mtoplv(icml,itim,2)              &
                      , mtoplv(icml,itim,3) )
    mdatopo = MAX( mdatop - 1 ,  1 )
    mdabotl = MAX( mdabot     , i1 )
!   IF (mdatop <= ktp) mdatopo = mdatop

! preset 'missing values'

    klevmax = nobs - kbotlev + mdabot - mdatop + 1
    DO klev = 1 , klevmax
      klvi = ktopoi + klev
      oiml (noiu  ,klvi) = rmdi
      oiml (noiv  ,klvi) = rmdi
      oiml (noit  ,klvi) = rmdi
      oiml (noiqd ,klvi) = rmdi
      oiml (noirh ,klvi) = rmdi
      oiml (noiuqc,klvi) = c1
      oiml (noitqc,klvi) = c1
      oiml (noiqqc,klvi) = c1
      oiml (noivcr,klvi) = c1
      oiml (noiulr,klvi) = c0
      oiml (noitlr,klvi) = c0
      oiml (noiqlr,klvi) = c0
    ENDDO
    klev    = 0
    klvi    = 0
    lexisto = .FALSE.

! Loop over model levels
! ----------------------

    loop_over_model_levels:  DO km = mdabotl , mdatopo , -1


!  Section 3a: Add obs. increments at observat. levels below current model level
!-------------------------------------------------------------------------------

! decide which obs. increment will be required

      lextu =      (km == ke-mbotlv(icml,itim,1))                              &
              .OR. (km ==    mtoplv(icml,itim,1))                              &
              .OR. (zobbdy(kobs,nbtp) >= topobs (1)) .AND. (kobtyp /= nsattv)  &
              .OR. (zobbdy(kobs,nbtp) >= topobtv(1)) .AND. (kobtyp == nsattv)
      lextt =      (km == ke-mbotlv(icml,itim,2))                              &
              .OR. (km ==    mtoplv(icml,itim,2))                              &
              .OR. (zobbdy(kobs,nbtp) >= topobs (3)) .AND. (kobtyp /= nsattv)  &
              .OR. (zobbdy(kobs,nbtp) >= topobtv(3)) .AND. (kobtyp == nsattv)
      lextq =      (km == ke-mbotlv(icml,itim,3))                              &
              .OR. (km ==    mtoplv(icml,itim,3))                              &
              .OR. (zobbdy(kobs,nbtp) >= topobs (4)) .AND. (kobtyp /= nsattv)  &
              .OR. (zobbdy(kobs,nbtp) >= topobtv(4)) .AND. (kobtyp == nsattv)

! WHILE-loop over remaining observation levels below current model level

      DO WHILE ((zobbdy(kobs,nbtlop) >= zlopml(icml,km)) .AND. (.NOT. ltop))

        zloprev = rmdi
        IF (klvi > 0)  zloprev = oiml(noilp,klvi)
        
        IF (      ((lextu) .OR. (lextt) .OR. (lextq))                          &
            .AND. (     (zobbdy(kobs,nbtlop)+0.0001_wp < zloprev)          &
                   .OR. (kobtyp /= nsattv) .OR. (klev == 0))) THEN

! increase number of obs. increment levels

          klev    = klev + 1
          klvi    = ktopoi + klev
          lexisto = .FALSE.

! get obs. increments of horizontal wind

          IF ((lextu) .AND. (      BTEST( mzobbd(kobs,nbterr),nvru ))          &
                      .AND. (.NOT. BTEST( mzobbd(kobs,nbtqcf),nvru ))          &
                      .AND. (vimtoob(kobs,1) > rmdich)) THEN
            lexisto = .true.
            oiml (noiu  ,klvi) = zinv_ob(kobs,nbtu) - vimtoob(kobs,1)
            oiml (noiv  ,klvi) = zinv_ob(kobs,nbtv) - vimtoob(kobs,2)
            oiml (noiulr,klvi) = zwvip_o(kobs)
            IF ((lwruv) .AND. (     (oiml(noiu,klvi) > qcupr)                  &
                               .OR. (oiml(noiv,klvi) > qcupr)))                &
              WRITE( nupr,'(''Du,Dv '',A ,I5,2F8.2,F8.0)' ) ystid, klev        &
                   , oiml(noiu,klvi), oiml(noiv,klvi),col_p(km)
          ENDIF

! get obs. increments of temperature

          IF ((lextt) .AND. (      BTEST( mzobbd(kobs,nbterr),nvrt ))          &
                      .AND. (.NOT. BTEST( mzobbd(kobs,nbtqcf),nvrt ))          &
                      .AND. (vimtoob(kobs,3) > rmdich)) THEN
            lexisto = .true.
            oiml (noit  ,klvi) = zinv_ob(kobs,nbtt) - vimtoob(kobs,3)
            oiml (noitlr,klvi) = zwvip_o(kobs)
            IF ((lwrt ) .AND. (oiml(noit,klvi) > qctpr))                       &
              WRITE( nupr,'(''DT    '',A ,I5, F8.2,F8.0)' ) ystid, klev        &
                   , oiml(noit,klvi), col_p(km)
          ENDIF

! get obs. increments of specific & relative humidity

          IF ((lextq) .AND. (      BTEST( mzobbd(kobs,nbterr),nvrq ))          &
                      .AND. (.NOT. BTEST( mzobbd(kobs,nbtqcf),nvrq ))          &
                      .AND. (vimtoob(kobs,4) > rmdich)) THEN
            lexisto = .true.
            zrhobs  = MIN( c1 , zinv_ob(kobs,nbtrh) )
            zrhmod  = MIN( c1 , vimtoob(kobs,4) )
            oiml (noirh ,klvi) = zrhobs - zrhmod
            zpvsat  = fpvsw ( vimtoob(kobs,3) , b1, b2w, b3, b4w )
            zqvobs  = fpv2q ( zrhobs *zpvsat , zobbdy(kobs,nbtp) , rdv )
            zqvmov  = fpv2q ( zrhmod *zpvsat , zobbdy(kobs,nbtp) , rdv )
            oiml (noiqd ,klvi) = zqvobs - zqvmov
            oiml (noiqlr,klvi) = zwvip_o(kobs)
            IF ((lwrq ) .AND. (oiml(noirh,klvi) > qcqpr))                      &
              WRITE( nupr,'(''DRH   '',A ,I5, F8.4,F8.0)' ) ystid, klev        &
                   , oiml(noirh,klvi), col_p(km)
          ENDIF

! get additional information if at least 1 obs. increment
! exists at the current level, else omit this level

          IF (lexisto) THEN
            zemkpoi (klev) = EXP( - rdocp * zobbdy(kobs,nbtlop) )
            zrtdgoi (klev) = zuz2p_o(kobs)
            zpresoi (klev) = zobbdy(kobs,nbtp)
            oiml (noiz  ,klvi) = vimtoob(kobs,5)
            oiml (noith ,klvi) = vimtoob(kobs,3) *zeklop0 *zemkpoi(klev)
!           oiml (noiuqc,klvi) = c1
!           oiml (noitqc,klvi) = c1
!           oiml (noiqqc,klvi) = c1
            oiml (noilp ,klvi) = zobbdy(kobs,nbtlop)
            IF (modespr <= 1) oiml (noivcr,klvi) = c1 / zuz2p_o(kobs)
          ELSE
            klev = klev - 1
          ENDIF
        ENDIF

! finalize WHILE-loop: increase index of obs. level

        ltop = (kobs == nobs)
        kobs = MIN( INT( kobs + 1 ,iintegers) , nobs )
      ENDDO

!  Section 3b: Add observation increments at current model level
!---------------------------------------------------------------

! if obs. increment exists: increase number of obs. increment levels

      IF ((km >= mdatop) .AND. (km <= mdabot)) THEN
        lextu = ((viobtom(km,1) > rmdich) .AND. (ke-km >= mbotlv(icml,itim,1)) &
                                          .AND. (km    >  mtoplv(icml,itim,1)))
        lextt = ((viobtom(km,3) > rmdich) .AND. (ke-km >= mbotlv(icml,itim,2)) &
                                          .AND. (km    >  mtoplv(icml,itim,2)))
        lextq = ((viobtom(km,4) > rmdich) .AND. (ke-km >= mbotlv(icml,itim,3)) &
                                          .AND. (km    >  mtoplv(icml,itim,3)))
        lexisto = (lextu) .OR. (lextt) .OR. (lextq)
        ! current model level overwrites previous obs (sat retrieval) level,
        !   if it lies less than 0.0001*p (ca 0.1 hPa) above it  
        IF ((lexisto) .AND. (klev >= 1) .AND. (kobtyp == nsattv)) THEN
          klvi    = ktopoi + klev
          IF (oiml(noilp,klvi) < zlopml(icml,km) +0.0001_wp)               &
            klev = klev - 1
        ENDIF   
        klev    = klev + 1
        klvi    = ktopoi + klev
        IF (lwrite) WRITE( nupr,'(''mult_obs_incr: klev, km'',2I5) ') klev, km

! get obs. increments of horizontal wind

        IF (lextu) THEN
          oiml (noiu  ,klvi) = viobtom(km,1)
          oiml (noiv  ,klvi) = viobtom(km,2)
          oiml (noiulr,klvi) = zwvip_m(km,1)
          IF ((lwruv) .AND. (     (oiml(noiu,klvi) > qcupr)                    &
                             .OR. (oiml(noiv,klvi) > qcupr)))                  &
            WRITE( nupr,'(''Du,Dv '',A ,I5,2F8.2,F8.0)' ) ystid, klev          &
                 , oiml(noiu,klvi), oiml(noiv,klvi), col_p(km)
        ENDIF

! get obs. increments of temperature

        IF (lextt) THEN
          oiml (noit  ,klvi) = viobtom(km,3) - col_t(km)
          oiml (noitlr,klvi) = zwvip_m(km,2)
          IF ((lwrt ) .AND. (oiml(noit,klvi) > qctpr))                         &
            WRITE( nupr,'(''DT    '',A ,I5, F8.2,F8.0)' ) ystid, klev          &
                 , oiml(noit,klvi), col_p(km)
        ENDIF

! get obs. increments of specific & relative humidity

        IF (lextq) THEN
          zpvmod = fq2pv ( col_qv(km) , col_p(km) , rdv )
          zpvsat = fpvsw ( col_t (km) , b1, b2w, b3, b4w )
          oiml (noirh ,klvi) = viobtom(km,4) - MIN( zpvmod /zpvsat , c1 )
!         oiml (noirh ,klvi) = viobtom(km,4) - MIN( col_rh(km) , c1 )
          oiml (noiqd ,klvi) =   fpv2q( viobtom(km,4) *zpvsat, col_p(km), rdv) &
                               - col_qv(km)
          oiml (noiqlr,klvi) = zwvip_m(km,3)
          IF ((lwrq ) .AND. (oiml(noirh,klvi) > qcqpr))                        &
            WRITE( nupr,'(''DRH   '',A ,I5, F8.4,F8.0)' ) ystid, klev          &
                 , oiml(noirh,klvi), col_p(km)
        ENDIF

! get additional information if at least 1 obs. increment
! exists at the current level, else omit this level

        IF (lexisto) THEN
          zemkpoi (klev) = EXP( - rdocp * zlopml(icml,km) )
          zrtdgoi (klev) = zrtgkml(icml,km)
          zpresoi (klev) = col_p (km)
          oiml (noiz  ,klvi) = col_z (km)
          oiml (noith ,klvi) = col_th(km)
          oiml (noiuqc,klvi) = c1
          oiml (noitqc,klvi) = c1
          oiml (noiqqc,klvi) = c1
          oiml (noilp ,klvi) = zlopml(icml,km)
          IF (modespr <= 1) oiml (noivcr,klvi) = c1 / zrtgkml(icml,km)
        ELSE
          klev = klev - 1
        ENDIF
      ENDIF

! end of loop over model levels

    ENDDO loop_over_model_levels

! get number of levels in the obs. increment profile

    numlev = klev

! printout for control

    IF ((lfirst) .AND. (lwonl)) THEN
      WRITE( nupr,'(''  u-incr  v-incr  T-incr  Qd-incr RH-incr pressure''     &
                  &,'' height pot. T. '',A )' ) ystid
      DO klev = 1 , numlev
        klvi  = ktopoi + klev
!US  valgrind complained about too small values
!US     IF (ANY(oiml(:,klvi) < 1E-10_wp)) THEN
!US       WRITE( nupr,'(A)') ' some values of oiml not defined'
!US     ELSE
        WRITE( nupr,'(3F8.3, F10.5,F8.3, 2F8.0,F7.2)' )                        &
               oiml(noiu ,klvi), oiml(noiv ,klvi)                              &
             , oiml(noit ,klvi), oiml(noiqd,klvi)                              &
             , oiml(noirh,klvi), zpresoi   (klev)                              &
             , oiml(noiz ,klvi), oiml(noith,klvi)
!US     ENDIF
      ENDDO
    ENDIF

  ENDIF


!-------------------------------------------------------------------------------
!  Section 4: If required, modification of the potential temperature profile to
!             yield the spreading parameter a monotonous function of pressure
!-------------------------------------------------------------------------------

! check if the profile is stable, and needs to be modified

  DO klev = 2 , numlev
    klvi  = ktopoi + klev
    oiml (noiz,klvi) = MAX( oiml(noiz,klvi) , oiml(noiz,klvi-1) +0.1_wp )
  ENDDO
  lmono = .TRUE.
  mono: DO klev = 2 , numlev
    klvi  = ktopoi + klev
    lmono = (oiml(noith,klvi) >= oiml(noith,klvi-1))
    IF (.NOT. lmono) EXIT mono
  ENDDO mono
  klev = MIN( klev , numlev )
  IF (lwrite) WRITE( nupr,'(''mult_obs_incr: lmono '',L1,'', msprpar '',I2     &
                          &,'', instable below p='',F8.0)' )                   &
                     lmono, modespr, zpresoi(klev)
  IF ((.NOT. lmono) .AND. (modespr == 2)) THEN
    klvi = ktopoi + 2

! WHILE-loop over obs. incr. levels for checking for instabilities

    DO WHILE (klvi <= ktopoi+numlev)
      IF (oiml(noith,klvi) >= oiml(noith,klvi-1)) THEN
        klvi = klvi + 1
      ELSE

! If the layer is instable: construct a stable profile
! First get the minimum & maximum value of pot. temperature within the unstable
! part of the profile. Then get the first level further above / below with pot.
! temperature greater / smaller than that maximum / minimum

        kboti  = klvi - 1
        sprmax = oiml(noith,klvi-1)
        sprmxi = sprmax
        sprmin = oiml(noith,klvi  )
        kstop  = 0
        DO kll = klvi+1 , ktopoi+numlev
          sprmxi = MAX( sprmxi , oiml(noith,kll-1) )
          sprmin = MIN( sprmin , oiml(noith,kll  ) )
          IF (      (oiml(noith,kll-1) <  sprmax)                              &
              .AND. (oiml(noith,kll  ) >= sprmax)) THEN
            sprmax = sprmxi
            kstop  = kll
          ENDIF
        ENDDO
        IF (kstop == 0) THEN
          kstop  = ktopoi + numlev + 1
          sprmax = sprmxi
        ENDIF
        ksbot  = ktopoi
        DO kll = ktopoi+1 , kboti
          IF (oiml(noith,kll) < sprmin-epsy) ksbot = kll
        ENDDO
        IF ((ntstep <= 2) .AND. (lwonl))                                       &
          WRITE( nupr,'(''mult_obs_incr: instability: '',A                     &
                      &,'', ksbot, kstop'',3I4)') ystid, ksbot, kstop, ktopoi
! construct a stable (pot. temp. lin. in ln(p)) profile betw. the 2 levels found
        DO kll = ksbot+1 , kstop-1
          zlopf =  (oiml(noilp,ksbot+1) - oiml(noilp,kll    ))                 &
                 / (oiml(noilp,ksbot+1) - oiml(noilp,kstop-1))
          zpothoi   (kll-ktopoi)  =  oiml(noith,kll)
          oiml (noith,kll)  =  (c1-zlopf) *sprmin  +  zlopf *sprmax
        ENDDO
! printout for control
        IF ((lfirst) .AND. (lwonl)) THEN
          ksbot = MAX( ksbot , INT( ktopoi + 1      ,iintegers) )
          kstop = MIN( kstop , INT( ktopoi + numlev ,iintegers) )
          DO kll = ksbot , kstop
            klev = kll - ktopoi
            WRITE( nupr,'(3F8.3, F8.3, 2F8.0,2F7.2)' )                         &
                   oiml (noiu,kll), oiml(noiv ,kll)                            &
                 , oiml (noit,kll), oiml(noirh,kll)                            &
                 , zpresoi  (klev), oiml(noiz ,kll)                            &
                 , zpothoi  (klev), oiml(noith,kll)
          ENDDO
        ENDIF
        klvi = kstop
      ENDIF
    ENDDO
  ENDIF


!-------------------------------------------------------------------------------
!  Section 5: If profile does also have obs. increments on observation levels:
!             Get additional local information for later spreading
!-------------------------------------------------------------------------------

  IF (levdiff) THEN
!   ivrsa = 2 - MAX( momlhd(nmlob,nhuexi) , i0 )
!   ivrse = 2 + MAX( momlhd(nmlob,nhqexi) , i0 )
!   ivrsd = 1
!   IF (ivrsa == 2) ivrsa = 3 - MAX( momlhd(nmlob,nhtexi) , i0 )
!   IF (ivrse == 2) ivrse = 1 + MAX( momlhd(nmlob,nhtexi) , i0 )
!   IF (ivrse-ivrsa == 2) ivrsd = 2 - MAX( momlhd(nmlob,nhtexi) , i0 )
!   DO ivrs = ivrsa , ivrse , ivrsd

    IF (modespr <= 1) nois = noiz
    IF (modespr == 2) nois = noith
    DO ivrs = 1 , 3
      IF (ivrs == 1) noix = noiu
      IF (ivrs == 2) noix = noit
      IF (ivrs == 3) noix = noirh
      lexisto = .FALSE.
      DO klev = 1 , numlev
        klvi  = ktopoi + klev
        IF ((.NOT. lexisto) .AND. (oiml(noix,klvi) > rmdich)) THEN
          lexisto = .TRUE.
          nba = 2*itim - 1
          zspobml(icml,nba,ivrs) = oiml(nois,klvi)
          zpobml (icml,nba,ivrs) = zpresoi(klev)
          zemkpml(icml,nba,ivrs) = zemkpoi(klev)
          zrtdgml(icml,nba,ivrs) = zrtdgoi(klev)
        ENDIF
        IF (oiml(noix,klvi) > rmdich) THEN
          nba = 2*itim
          zspobml(icml,nba,ivrs) = oiml(nois,klvi)
          zpobml (icml,nba,ivrs) = zpresoi(klev)
          zemkpml(icml,nba,ivrs) = zemkpoi(klev)
          zrtdgml(icml,nba,ivrs) = zrtdgoi(klev)
        ENDIF
      ENDDO
    ENDDO
  ENDIF


!-------------------------------------------------------------------------------
! End of module procedure mult_obs_increment
!-------------------------------------------------------------------------------

END SUBROUTINE mult_obs_increment


!-------------------------------------------------------------------------------
!+ Module procedure for building profiles of humidity increments for GPS reports
!-------------------------------------------------------------------------------

SUBROUTINE mult_obs_gps ( icml , lqc , lqconly , ktopoi , ke , col_t, col_qv   &
                        , col_p, col_z, col_qc, col_qrs, col_th, col_hhl       &
                        , nexceiv , numlev )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure builds vertical profiles of humidity observation
!   increments for Global Positioning System (GPS) reports. A threshold
!   quality control is included.
!
! Method:
!   The basic idea is to build "pseudo" GPS observations of specific humidity
!   from the observed integrated water vapour (IWV) or zenith path delay (ZPD).
!   - The model humidity field is used as first guess profile. This profile is
!     iteratively adjusted until its integration on all model levels matches the
!     observed GPS IWV value. When the adjusted specific humidity exceeds the 
!     saturation threshold, the value is set to the model saturation specific 
!     humidity. The procedure is terminated after 20 iterations.
!   - Note that the model humidity is integrated from the top down to the height
!     of the GPS antenna for the determination of a model IWV value.
!   - A threshold quality control is applied to the IWV value.
!   - Height dependent quality weights are assigned to the humidity observation
!     increments. At each level k the quality weight is w(k)= dz * qsat / w max.
!   - The observation increments are set to missing above about 300 hPa.
!   Also included: - a seasonal daytime-dependent bias correction
!                  - adjustment due to difference between water / ice saturation
!                    to make observed IWV compatible to model without cloud ice
!
! Written by        :  Maria Tomassini  , DWD  (original version: 01.01.01)
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

  INTEGER (KIND=iintegers), INTENT (IN)         ::       &
    icml             ,& ! index of obs. sta. in data record to be used further
    ktopoi           ,& ! number vertical levels used in the observation incre-
                        ! ment array 'oiml' prior to the current report
    ke                  ! number of vertical (main) levels in model column

  LOGICAL                 , INTENT (IN)         ::       &
    lqc              ,& ! threshold quality control to be applied
    lqconly             ! threshold quality control without other processing

  REAL    (KIND=wp   ),     INTENT (IN)         ::       &
                        ! model column (at the observation location) of:
    col_t     (ke)   ,& ! temperature                                   (  K  )
    col_th    (ke)   ,& ! potential temperature                         (  K  )
    col_qv    (ke)   ,& ! specific water vapor content                  (kg/kg)
    col_qc    (ke)   ,& ! specific cloud water content                  (kg/kg)
    col_qrs   (ke)   ,& ! spec. cont. of hydrometeors excl. cloud water (kg/kg)
    col_p     (ke)   ,& ! pressure (full value)       on main levels    ( Pa  )
    col_z     (ke)   ,& ! geometrical height          of main levels    (  m  )
    col_hhl   (ke+1)    ! geometrical height          of half levels    (  m  )

  INTEGER (KIND=iintegers), INTENT (INOUT)      ::       &
    nexceiv             ! number of IWV obs incr. reports exceeding array size

  INTEGER (KIND=iintegers), INTENT (OUT)        ::       &
    numlev              ! number of levels in the obs. increment profile


! Local parameters: 
! ----------------
  INTEGER (KIND=iintegers)   ,  PARAMETER ::  &
    itmax  = 20_iintegers      ! maximum number of iterations 
  
  REAL    (KIND=wp   )       ,  PARAMETER ::  &
    c1000r =  0.001_wp  ,& ! difference limit to continue iteration  
                               ! (if this limit is too large (e.g. 0.01) then
                               !  observation increments are set to 0 too often)
    thriq  =  2.0_wp    ,& ! minimum observed IWV value for assimilation 
    zfrmin = 0.60_wp       ! minimum value of zfracb to consider
                               ! closest bottom level for GPS increments 
    
     
! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    modespr          ,& ! mode of spreading (horizontal, along model levels
                        !                   ,along isentropes)
    ngpob            ,& ! index of report in the ODR (obs data record)
    km               ,& ! vertical (loop) index of model level (== obs)
    klev             ,& ! vertical (loop) index of obs. increment level
    klvi             ,& ! vertical index in 'oiml'
    ivrs             ,& ! index of observed quantity: 1=(u,v); 2=T; 3=RH
    mdabot           ,& ! index of lowest model level with GPS obs. increments
    mdabotl          ,& ! index of lowest model level with GPS pseudo obs. 
    mdatop           ,& ! index of uppermost model level with GPS obs increments
    mdatopo          ,& ! index of uppermost model level with GPS pseudo obs 
    nbotop , nba     ,& ! loop indices
    kkk              ,& ! counter 
    iter             ,& ! iteration counter 
    isat             ,& ! number of levels with obs at saturation  
    irflg               ! rejection flag 0=not rej
                        !                1=rej with qctic
                        !                2=rej with thriq

  REAL    (KIND=wp   )     ::  &
    ziqobs           ,& ! obs.  integrated mass of water vapour q: obs.  IWV
    ziqgps           ,& ! IWV reported in the GPS input file
    ziqint           ,& ! IWV computed from ZPD using interpolated press + temp
    ziqmod           ,& ! model integrated mass of water vapour q: model IWV
    ziqsat           ,& ! IWV from model at saturation q-sat (from GPS height)
    ztqsat           ,& ! as ziqsat, but from model orography
    ziqmic           ,& ! IWV from ice-to-water corrected model q-ice
    ziqite           ,& ! IWV from pseudo-obs. of q at iteration step 'iter'
    zqiq             ,& ! quotient 'ziqobs / ziqite' used to derive q-pseudo-obs
    zpv              ,& ! water vapour pressure
    omykf            ,& ! local part of the vertical nudging weight
    fisd             ,& ! height diff betw. station altitude and model orography
    timdif           ,& ! difference of the observation time to the model time
    tabdif           ,& ! distance   of the observation time to the model time
    zfracb           ,& ! weight for the contribution to model IWV  
                        ! from the lowest layer  
    zmxwei           ,& ! maximum weight to normalize quality weights
    qctiq            ,& ! time-dependent quality control threshold for IWV
    pgps             ,& ! pressure value from GPS report [hPa]
    tgps             ,& ! temperature value from GPS report [K]
    pint             ,& ! model pressure interpolated at GPS height [Pa]
    tint                ! model temperature interpolated at GPS height [K]

  LOGICAL                  ::  &
    lveridat         ,& ! fill simulated obs body (for writing to feedback file)
    lexisto          ,& ! obs increments exist
    lexchek          ,& ! obs incr. exist, which may become active in spatial
                        ! consistency check
    lwrite           ,& ! printout for control at present station / node / time
    lqciwv           ,& ! use obs for IWV spatial consistency checking
    lwrqc               ! data which are close the current time, and which are
                        ! rejected by the threshold quality control, are printed
                        ! to a separate file at the current timestep

! Local (automatic) arrays:
! -------------------------

  REAL    (KIND=wp   )     ::  &
    zqvobs   (ke)    ,& ! observed specific humidity 
                        !(model    specific humidity before iteration)
    zqvsat   (ke)    ,& ! model    specific humidity at saturation
    zrhobs   (ke)    ,& ! observed relative humidity  (not generalized)
    zrhmod   (ke)    ,& ! model    relative humidity  (not generalized)
    ztd      (ke)    ,& ! model dewpoint (temperature)
    zrhodz   (ke)    ,& ! density of moisty air * thickness of model layer
    zgpwei   (ke)    ,& ! quality weights for the GPS observation increments
    zprbuf   (ke)    ,& ! buffer for printout only
    zyqc     (11)       ! information for QC control messages
!
!------------ End of header ----------------------------------------------------

 
!-------------------------------------------------------------------------------
! Begin Subroutine mult_obs_gps
!-------------------------------------------------------------------------------

  ngpob = igpadm(ista,itim)

!-------------------------------------------------------------------------------
!  Section 1: Application of the 'observation operator':
!             - Determination of 'observed' GPS IWV either by computing it from
!               observed ZPD and interpolated model pressure and temperature
!               (i.e. don't use the reported IWV value !);
!             - Determination of simulated zenith total path delay (ZPD) derived
!               from model profiles
!             - Computation of the model profiles and vertically integrated
!               model quantities that are required for the quality control,
!               spatial consistency check, and nudging later on
!             - (Also determination of the levels used to compute the model IWV
!                (from closest level below GPS to top model level))
!-------------------------------------------------------------------------------

  lveridat  = ((mogphd(ngpob,nhqcfw) >= 0) .AND. (lvofoi))
  timdif    =   ogphed(ngpob,nhtime) - acthr
  tabdif    =   ABS( timdif )
  lwrqc     = (lqcall) .AND. (ABS( c3600*timdif+epsy ) <= cqcwbox)

  CALL zpd_iwv_obs_operator ( ke, col_p, col_hhl, col_t, col_qv, col_qc        &
                            ,     col_qrs, g, r_d, r_v, t0_melt, degrad        &
                            ,     b1, b2w, b3, b4w, b2i, b4i, madj_hum         &
                            , ogpbdy(ngpob,nbgzpd), ogphed(ngpob,nhalt)        &
                            , ogphed(ngpob,nhjlat), lveridat , sgpbdy(:,ngpob) &
                            , ziqint, ziqobs, ziqmod, ziqmic, ziqsat, ztqsat   &
                            , zqvsat, zrhmod, zrhodz, mdabotl, zfracb, pint )
! =========================

  !   apply bias correction to GPS Integrated water vapour
! IF (lgpsbias)       ziqobs  =  ziqobs + ogpbdy(ngpob,nbgbia)

  !   set negative values to zero  (values < 2 mm will be not nudged anyway!)
! IF (ziqobs <= c0)   ziqobs  =  c1000r

! Set adjusted (bias-corrected) IWV value in ODR !
! ----------------------------------------------

  ogpbdy (ngpob,nbgiwa) = ziqobs

! Current version: set bias correction to 'ziqobs' minus reported IWV 'ziqpgs'
! ----------------------------------------------------------------------------

  ziqgps = ogpbdy (ngpob,nbgiwv)
  IF (ziqgps > rmdich)  ogpbdy (ngpob,nbgbia) = ziqobs - ziqgps

!-------------------------------------------------------------------------------
!  Section 2:  Threshold quality control of IWV (including minimum limit check)
!-------------------------------------------------------------------------------

  IF (lqc) THEN

    CALL iwv_quality_cntl ( ziqobs, ziqmod, ziqsat, qcsiq, qcciq, qctf(4)      &
                          , tabdif, lwrqc , mogpbd(ngpob,1:mxgbdf)             &
                          ,                 mogphd(ngpob,1:mxghdf) , zyqc )
  ! =====================

    ! fill record for later printing for control of QC
    ! ------------------------------------------------

    IF ((zyqc(1) >= epsy) .AND. (ntotqc < maxqcp)) THEN
      ntotqc = ntotqc + 1
      yyqc (ntotqc   ) = ystid
      myqc (ntotqc, 1) = mogphd(ngpob,nhcode)
      myqc (ntotqc, 2) = NINT( zyqc(1) )
      oyqc (ntotqc, 1) = ogphed(ngpob,nhtime)
      oyqc (ntotqc, 2) = ogpbdy(ngpob,nbgp)   * 0.01_wp
      oyqc (ntotqc, 3) = ogphed(ngpob,nhjlat)
      oyqc (ntotqc, 4) = ogphed(ngpob,nhilon)
      oyqc (ntotqc, 5) = zyqc (5)
      oyqc (ntotqc, 6) = zyqc (6)
      oyqc (ntotqc, 7) = zyqc (7)
      IF (.NOT. lwrqc)  myqc (ntotqc,1) = 0
    ENDIF
  ENDIF

!CS: this may be replaced by a check on wet delay in the obs pre-processing
! IF (ziqobs <= thriq) THEN
!   mogpbd (ngpob,nbgerr) = IBCLR ( mogpbd(ngpob,nbgerr), nvriwv )
!   mogpbd (ngpob,nbgerr) = IBCLR ( mogpbd(ngpob,nbgerr), nvrzpd )
!   mogpbd (ngpob,nbgflg) = IBSET ( mogpbd(ngpob,nbgflg), nvfgbp+nvfbps(3) )
! ENDIF

!-------------------------------------------------------------------------------
!              Subsequent sections are used for nudging only
!-------------------------------------------------------------------------------
!  Section 3:  Preparation of spatial consistency check for IWV.
!              To be compatible with IWV-QC of radiosonde humidity, the
!              IWV obs increment with respect to model orography is stored here.
!              (Also determination of a quality weight depending on the height
!               difference between GPS antenna and model orography.)
!-------------------------------------------------------------------------------

  lqciwv = (lqc) .AND. (liwvssc)

  IF (lqciwv) THEN
    kkk = kiobiv(MAX(nivtot,1))
    IF (nivtot <= 0) THEN
      nivtot = 1
      iqcliv (nivtot)  =  0
      iqcfiv (nivtot)  =  0
    ELSEIF (     (kiobiv(MAX(nivtot,1)) /= ista)                               &
            .OR. (ktypiv(MAX(nivtot,1)) /= kobtyp)) THEN
      IF (nivtot < maxivq) THEN
        nivtot = nivtot + 1
        iqcliv (nivtot)  =  0
        iqcfiv (nivtot)  =  0
      ELSE
        nexceiv = nexceiv + 1
      ENDIF
    ENDIF
    !   quality weight dep. on station height difference
    fisd  =  ogphed(ngpob,nhalt) - col_hhl(ke+1)
    fisd  = ((fdoro(2)-c1)/c2 + SIGN( (fdoro(2)+c1)/c2 , fisd) ) * fisd
    omykf =  EXP( - (fisd / doromx(2)) **2 )
    !   QC threshold, should be the same as in 'iwv_quality_cntl'
    !                 except that 'ziqsat' should be replaced by 'ztqsat'
    qctiq = (c1 + qctf(4) *ABS(timdif)) * (qcsiq *ziqsat + qcciq)
!   qctiq = (c1 + qctf(4) *ABS(timdif)) * (qcsiq *ztqsat + qcciq)
    ystidiv  (nivtot)     =  ystid
    ioiv     (nivtot)     =  io
    joiv     (nivtot)     =  jo
    ioiv_tot (nivtot)     =  io_tot
    joiv_tot (nivtot)     =  jo_tot
    ktypiv   (nivtot)     =  ngps
    kioiiv   (nivtot)     =  icml
    kiobiv   (nivtot)     =  ista
    zactiv (nivtot,itim)  =  ogphed(ngpob,nhtime)
    zoiciv (nivtot,itim)  =  (ziqobs - ziqmod) * ztqsat / ziqsat
    zqcfiv (nivtot,itim)  =  omykf
    zsativ   (nivtot)     =  ztqsat
    zmodiv   (nivtot)     =  ziqmod * ztqsat / ziqsat
    iqcliv   (nivtot)     =  iqcliv(nivtot) + itim
    IF ((ABS( ziqmod -ziqobs) <= qctiq) .AND. (ziqobs > thriq)) THEN
      IF ((mogphd(ngpob,nhpass) > 0) .AND. (lgps) .AND. (gnudggp > epsy)) THEN
! assign 'iqcflg' a negative sign for passive observations (only if lgps,
!                 and if the other obs with index 'nivtot' is not active)
        IF (iqcfiv(nivtot) <= 0)  iqcfiv (nivtot)  =  iqcfiv(nivtot) - itim
      ELSE
        iqcfiv (nivtot)   =  MAX( iqcfiv(nivtot),0 ) + itim
      ENDIF
    ENDIF
    IF ((MOD(NINT(ntstep*dtdeh*3600), 1800) <= NINT(tconbox)-1) .AND. (lwonl)) &
      WRITE( nupr,'("IWV oi:",A ,8I4,F6.2,3F6.2,F6.3,F6.2)' )  ystid, kobtyp   &
                   , ista, icml, nivtot, kkk, iqcliv(nivtot), iqcfiv(nivtot)   &
           , ntstep, ogphed(ngpob,nhtime), zoiciv(nivtot,itim), zmodiv(nivtot) &
                   , ztqsat, omykf, ziqmod
  ENDIF

!-------------------------------------------------------------------------------
!  Section 6:  Computation of "pseudo-observed" GPS humidity profile.
!              If the observed GPS IWV is greater than model IWV at saturation,
!              the observed profile is set equal to the model humidity profile
!              at saturation.
!-------------------------------------------------------------------------------

  lwrite = (lfirst) .AND. (io == ionl) .AND. (jo == jonl) .AND. (lwonl)
  mdatopo = 1
! Station pressure in hPa
  pgps   = ogpbdy (ngpob,nbgp) / 100._wp
! Station temperature (K)
  tgps   = ogpbdy (ngpob,nbgt)

  DO km = 1 , ke
    !--   Model profile as as first guess of observed profile
    zqvobs (km)  =  col_qv(km)
  ENDDO
  zmxwei = c1000r
  DO km = ke , 1 , -1
    ! Quality weight to be assigned to the observation increments
    zgpwei(km) = (col_hhl(km) - col_hhl(km+1)) * zqvsat(km)
    IF (zmxwei < zgpwei(km)) zmxwei = zgpwei(km)
    zpv    = fq2pv( col_qv(km) , col_p(km) , rdv )
    IF (lwrite) ztd(km) = ftd( zpv, b1, b2w, b3, b4w )
  ENDDO
  IF (lwrite) THEN
    DO km = mdabotl , mdatopo , -1
      WRITE( nupr,'(3F9.1,F9.2,I4)' )                                          &
             col_t(km), ztd(km) , col_p(km), zrhmod(km), km
    ENDDO
  ENDIF

  iter    = 0
  irflg   = 0
  lexisto = .TRUE.
  lexchek = .FALSE.
  IF ((lqc) .AND. (ziqobs > thriq))  lexchek = .TRUE.

! IF ((ogpbdy(ngpob,nbgtze) < c0) .OR. (lqconly)) THEN
  IF (     (.NOT. BTEST( mogpbd(ngpob,nbgerr),nvriwv ))                        &
      .OR. (      BTEST( mogpbd(ngpob,nbgqcf),nvriwv )) .OR. (lqconly)) THEN
!   IF (ogpbdy(ngpob,nbgtze) < c0)  irflg  = 3
    IF (     (.NOT. BTEST( mogpbd(ngpob,nbgerr),nvriwv ))                      &
        .OR. (      BTEST( mogpbd(ngpob,nbgqcf),nvriwv )))  irflg  = 3
    lexisto =.FALSE.
    IF (.NOT. lexchek)  iter = itmax + 2
  ENDIF
 
  zprbuf = rmdi
  isat   = i0
  ziqite = ziqmod

! Observed IWV greater than model IWV at saturation 
  IF (((lexisto) .OR. (lexchek)) .AND. (ziqobs >= ziqsat-epsy)) THEN
    DO km = mdabotl, mdatopo, -1
      zqvobs(km) =  zqvsat(km)
      zrhobs(km) =  c1
    ENDDO
    ziqite = ziqsat
    iter   = itmax + 1
  ENDIF

  comp_gps_qv : DO WHILE ( iter < itmax )

! Quotient to scale guess of q- pseudo-obs profile from previous iteration
    zqiq  =  ziqobs / ziqite

! Stop iteration if IWV from pseudo-obs q is very close to observed IWV 
    IF( ABS( zqiq - c1 ) < c1000r ) THEN
      EXIT comp_gps_qv
    ELSE

! Compute new guess for pseudo-obs q-profile and corresponding IWV
      ziqite = c0
      isat   = i0
      DO km = mdabotl , mdatopo , -1
        zqvobs(km) = zqvobs(km) * zqiq
        IF ( zqvobs(km) > zqvsat(km) ) THEN
          zqvobs(km) = zqvsat(km)
          isat = isat + 1
        ENDIF
        ziqite = ziqite + zrhodz(km) * zqvobs(km)
      ENDDO
      iter = iter + 1
    ENDIF
    tint  =  col_t(ke)
    IF ((lwonl) .AND. (ntstep <= 1)) THEN
      WRITE( nupr,'(''MULT_GPS '',A5,I1,F4.1,I4,2I3,7F6.2,F5.2,2F8.1,2F6.1)' ) &
             ystid,irflg,ogphed(ngpob,nhtime),ntstep,iter,isat,ziqgps,ziqint   &
            ,ziqobs,ziqite,ziqmod,ziqmic,ziqsat,zfracb,pint,pgps,tint,tgps
    ENDIF

  END DO comp_gps_qv 

  IF ( iter <= itmax ) THEN
    DO km = mdabotl , mdatopo , -1
      zrhobs(km) =   fq2pv ( zqvobs(km) , col_p(km) , rdv )                    &
                   / fpvsw ( col_t(km) , b1, b2w, b3, b4w )
    ENDDO
  ENDIF

! Print info every half hour at GPS time
! IF ((ntstep <= 1) .OR.(ntstep==48) .OR.(ntstep==138) .OR.(ntstep==228))      &
!   PRINT '(4X,A5,I1,F4.1,2I3,7F6.2,F4.1,2F7.1)' , ystid                       &
!           ,irflg,ogphed(ngpob,nhtime),iter,isat,ziqgps,ziqint                &
!           ,ziqobs,ziqite,ziqmod,ziqmic,ziqsat,zfracb,pint,pgps
! Prepare for printing info every half hour at GPS time
  IF ((ntstep <= 1) .OR.(ntstep==48) .OR.(ntstep==138) .OR.(ntstep==228)) THEN
    IF (ke >= 13) THEN
      zprbuf (1) = ogphed(ngpob,nhtime)
      zprbuf (2) = REAL ( iter, wp )
      zprbuf (3) = REAL ( isat, wp )
      zprbuf (4) = ziqgps
      zprbuf (5) = ziqint
      zprbuf (6) = ziqobs
      zprbuf (7) = ziqite
      zprbuf (8) = ziqmod
      zprbuf (9) = ziqmic
      zprbuf(10) = ziqsat
      zprbuf(11) = zfracb
      zprbuf(12) = pint
      zprbuf(13) = pgps
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
!  Section 7: Profile with observation increments on main model levels
!             between model level close to GPS and approx. 300 hPa
!  Section 2: Determination of the levels with observation increments:
!             - lowest model level is - closest level below GPS if zfracb > 0.6
!                                     - closest level above GPS if zfracb <=0.6
!             - uppermost model level is at about 300 hPa
!-------------------------------------------------------------------------------

  modespr = msprpar

  mdabot = mdabotl

! If weight for lowest layer too small decrease index of lowest level
! with obs increments
  IF ( zfracb < zfrmin ) mdabot = mdabot - 1

! Number of model levels above surface without obs increments
  mbotlv(icml,itim,3) = ke - mdabot

! Number of model levels at the top without obs increments
  mdatop = kml300
  mtoplv(icml,itim,3) = mdatop - 1

! Default number of levels in the obs increments profile
  numlev = mdabot - mdatop + 1

  IF ((lexisto) .OR. (lexchek)) THEN

! Get GPS obs increments of specific and relative humidity, set all other
! increments as missing values, determine pressure, height and potential temp.,
! get quality weight of the GPS observation increment

    klev = 0
    DO km = mdabot , mdatop , -1
      klev = klev + 1
      klvi = ktopoi + klev
      oiml (noiu  ,klvi) = rmdi
      oiml (noiv  ,klvi) = rmdi
      oiml (noit  ,klvi) = rmdi
      oiml (noiqd ,klvi) = zqvobs(km) - col_qv(km)
      oiml (noirh ,klvi) = MIN(c1,zrhobs(km)) - MIN(c1,zrhmod(km))
      oiml (noiz  ,klvi) = col_z (km)
      oiml (noith ,klvi) = col_th(km)
      oiml (noiuqc,klvi) = c1
      oiml (noitqc,klvi) = c1
      oiml (noiqqc,klvi) = zgpwei(km) / zmxwei
      oiml (noilp ,klvi) = zlopml(icml,km)
      oiml (noiulr,klvi) = c0
      oiml (noitlr,klvi) = c0
      oiml (noiqlr,klvi) = c0
      oiml (noivcr,klvi) = c1 / zrtgkml(icml,km)
      IF (modespr == 2) oiml (noivcr,klvi) = c1
! misuse 'noiv' for printout temporarily
      oiml (noiv  ,klvi) = zprbuf(klev)
    ENDDO

! get additional local information for later spreading

    DO nbotop = 1 , 2
      nba = 2*(itim-1) + nbotop
      DO ivrs = 1 , 3
        IF (nbotop == 1) km = ke - mbotlv(icml,itim,ivrs)
        IF (nbotop == 2) km = 1  + mtoplv(icml,itim,ivrs)
        zspobml(icml,nba,ivrs) = zsprml(icml,km)
        zpobml (icml,nba,ivrs) = col_p(km)
        zemkpml(icml,nba,ivrs) = EXP( - rdocp * zlopml(icml,km) )
        zrtdgml(icml,nba,ivrs) = zrtgkml(icml,km)
      ENDDO
    ENDDO

! printout for control of vertical weights

!   IF (ntstep <= 1) THEN
!     PRINT '(''GPS_OIML '', A5,''noiqqc,noiqd,noirh'')',  ystid
!     klev = 0
!     DO km = mdabot , mdatop , -1
!       klev = klev + 1
!       klvi = ktopoi + klev
!       PRINT '(2F10.5,F8.3,I4)',                                              &
!              oiml(noiqqc,klvi), oiml(noiqd,klvi), oiml(noirh,klvi), km
!     ENDDO
!   ENDIF

    IF (.NOT. lexisto) THEN
      ivrs   = 3
      zwtml (icml,itim  ,ivrs) = c0
      zwtml (icml,itim+2,ivrs) = c0
    ENDIF

  ELSE
    numlev = 0
    ivrs   = 3
    zwtml (icml,itim  ,ivrs) = c0
    zwtml (icml,itim+2,ivrs) = c0
  ENDIF

! printout for control

  IF ((lwrite) .AND. ((lexisto) .OR. (lexchek))) THEN
    klev = 0
    DO km = mdabot , mdatop , -1
      klev = klev + 1
      klvi = ktopoi + klev
      WRITE( nupr,'(3F8.3, F10.5,F8.3, 2F8.0,F7.2)' )                          &
             oiml(noiu ,klvi), oiml(noiv ,klvi)                                &
           , oiml(noit ,klvi), oiml(noiqd,klvi)                                &
           , oiml(noirh,klvi), col_p(km)                                       &
           , oiml(noiz ,klvi), oiml(noith,klvi)
    ENDDO
  ENDIF


!-------------------------------------------------------------------------------
! End of module procedure mult_obs_gps
!-------------------------------------------------------------------------------

END SUBROUTINE mult_obs_gps


!-------------------------------------------------------------------------------
!+ Module procedure for organizing obs operators and QC for multi-level reports
!-------------------------------------------------------------------------------

SUBROUTINE mult_org_localinfo ( nexceml , nexceiv )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure prepares the lateral spreading of observation
!   increments from multi-level reports by computing (or organizing) all the
!   required 'local' information (i.e. observation increments, and further
!   parameters on the observations and their location) which will later
!   be broadcast to other processors.
!
! Method:
!   First computation of temporal weights, and determination whether a station
!   is processed further (condition: the report must contain at least 1 obser-
!   vation with non-zero temporal weight (prior to the threshold quality con-
!   trol)).
!   If so, observation increments (incl. additional information on the profile)
!   are obtained by call of procedures, and other quantities depending only on
!   the horizontal station location as well as some quantities on the vertical
!   extent of the influence of the report are computed straightforwardly.
!
! Written by        :  Christoph Schraff, DWD  (original version: 05.09.97)
! Current Code Owner:  Christoph Schraff, DWD
!
! Modification comment:
!   Date:   13.08.98   Quality control also for passive data, for setting flags.
!                      Bug correction for setting temporal weights to zero if no
!                      obs. increment levels exist.      Christoph Schraff
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
    nexceml          ,& ! number of obs increment reports exceeding array size
    nexceiv             ! number of IWV obs incr. reports exceeding array size

! Local parameters: None
! ----------------

! INTEGER (KIND=iintegers) , PARAMETER ::  &
!   mxoutrs = 20    ! max. number of reports where spreading is as by 'msprpar'

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    mxoiml           ,& ! max. number of levels in OIR (obs. increment record)
    istaml           ,& ! loop index of observing station
    itic             ,& ! time index of report
    icml             ,& ! index of obs. sta. in data record to be used further
    ivrs             ,& ! index of observed quantity: 1=(u,v); 2=T; 3=RH
    ivar             ,& ! index of observed quantity: 1=(u,v); 3=T; 4=RH
    nmlob  , nmlo2   ,& ! indices of reports in the ODR (obs data record)
    nmloba           ,& ! non-zero index of report in the ODR (obs data record)
    ngpob  , ngpo2   ,& ! indices of reports in the ODR (obs data record)
    ngpoba           ,& ! non-zero index of report in the ODR (obs data record)
    ntvob  , ntvo2   ,& ! indices of reports in the ODR (obs data record)
    ntvoba           ,& ! non-zero index of report in the ODR (obs data record)
!   nmloi            ,& ! index of report in the obs. increment record
    mckvifl          ,& ! initial values for lowest/uppermost influenced levels
    numlev           ,& ! number of levels in the obs. increment profile
    ktopoi           ,& ! number vertical levels used in the observation incre-
                        ! ment array 'oiml' prior to the current report
    kcdtyp           ,& ! CMA observation code type
    kbotlev          ,& ! number of obs levels below lowest    model level
    ktoplev          ,& ! number of obs levels above uppermost model level
    kat    , kbt     ,& ! model levels adjacent to the standard atm. tropopause
    kal    , kbl     ,& ! model levels adjacent to equidistant points at obs loc
    iv1    , iv1b    ,& ! start indices related to processed variables
    mexi   (3)       ,& ! =1: active ; =0: no ; =-1: passive obs exist
                        !  (index 1: wind ; 2: temperature ; 3: humidity)
    mbotlt , mtoplt  ,& ! lowest / top model level with T-obs increment
    nlqc             ,& ! number of new QC messages
    ntqc             ,& ! number of total QC messages
    nqc              ,& ! loop index for QC messages
    nhxexi           ,& ! ODR index for existance of 'good' obs of present var.
    nhvcbx           ,& ! /  correction factor to ver-   \  -  at base of report
    nhvctx           ,& ! \  tical correlation scale     /  -  at top  of report
    nuprin              ! file unit number for control output (no output if < 0)

  INTEGER (KIND=iintegers) ::  &
    kthr             ,& ! level index at upper limits where use horiz infl radii
    ilva   , ilvb    ,& ! indices of correlation scale table levels that are
                        ! adjacent to the specified model level
    kk     , mv      ,& ! loop indices
    ni2    , ni4     ,& ! loop indices
    nba    , nbotop  ,& ! loop indices
             ilev    ,& ! loop indices
    itica  , itice   ,& ! loop limits
    ivrsa  , ivrse   ,& ! loop limits
    npspr  , mxprspr ,& ! loop index and limit for printout
    nspe   , nspa    ,& ! interval of indices for printout
    istat  , nstat      ! status of memory (de-)allocation

  REAL    (KIND=wp   )     ::  &
    zrwt             ,& ! temporal weight relevant for horiz. correlation scale
    timdif           ,& ! time difference between the 2 obs at same obs. station
    tabdif           ,& ! time distance between the 2 obs at same obs. station
    tdiffw           ,& ! time differ. to the begin. of the nudging time window
    zoblat           ,& ! latitude of observation
    zobdps           ,& ! pressure obs increment valid at lowest model level
    rdg              ,& ! R / g
    zlop0            ,& ! log of reference pressure p0r
    zloptp           ,& ! log of pressure at tropopause of standard atmosphere
    zlptf            ,& ! weight factor for vertical interpol. to the tropopause
    zflop            ,& ! weight factor for vertical interpol. to equidist. pts.
    psign            ,& ! sign, depending if lower or upper cut-off is computed
    svcutof (3)      ,& ! vertical cut-off radii
    svcutos          ,& ! adjusted vertical cut-off radius
    zcut             ,& ! spreading parameter at vertical cut-off
    zsob             ,& ! spreading parameter at observation level
    zupper           ,& ! upper height limits where horiz. radii of infl. used
    zlopco           ,& ! log( pressure ) at level relevant for area of influen.
    zf               ,& ! weight factor for vertical interpol. to that level
    rvinfl           ,& ! vertically dependent part of horizontal correl. scale
    zrtinfl          ,& ! actual horizontal correl. scale
    zczvcut          ,& ! initial values for vertical cut-offs
    zispra , zisprb  ,& ! values of spreading parameter at vertic. equidist. pts
    ztexi  , zqexi   ,& ! values indicating existence of data
    zvcs  (2)        ,& ! vertical correlation scale below / above profile
    zvcut (2)        ,& ! lower / upper height limit influenced by T-increments
    zfcut               ! factor applied to scaled vertical cut-off radius

  LOGICAL                  ::  &
    lqc              ,& ! threshold quality control to be applied
    lqconly          ,& ! threshold quality control without other processing
    lqcflag (2)      ,& ! set threshold quality control flags in the ODR
    lqcfrst (2)      ,& ! pressure of report used for the first time
    ltakeob (2)      ,& ! use data of present station for further processing
    lqcpass (2)      ,& ! do quality control of passive data at obs. time
    lraso            ,& ! observation type is TEMP or PILOT
    lobcnv           ,& ! observation type is conventional: radiosonde, aircraft
    lobgps           ,& ! observation type is GPS
    lobrtv           ,& ! observation type is satellite retrieval
    lvirt            ,& ! virtual temperature observed, instead of temperature
    lreset           ,& ! reset icml and variables depending on icml
    lqcdz            ,& ! thickness quality control to perform / performed
    lredo            ,& ! redo section 3c of *mult_vertic_intpol*
    lrejtot          ,& ! reject total (temperature) profile due to (d)z-QC
    lobinc  (3)      ,& ! interpolated quantities are obs incr., not obs values
    lveridat         ,& ! fill simulated obs body (for writing to feedback file)
    lwrqc            ,& ! data with obs time close to current time, which are
                        ! rejected by the threshold quality control, are printed
                        ! to a separate file at the current timestep
    lprqcm  (2)      ,& ! printout 1 for control at present station/ node/ time
    lprfrs  (2)         ! printout 2 for control at present station/ node/ time

  CHARACTER (LEN=25)       ::  yzroutine
  CHARACTER (LEN=255)      ::  yzerrmsg
  INTEGER (KIND=iintegers) ::  izerror

! possibility of different spreading for different multi-level reports is not
! available
! LOGICAL                  ::  &
!   lthall          ,& ! .T. ==> spreading of all upper-air data not along model
!                      !         levels, but along levels according to 'msprpar'
! INTEGER (KIND=iintegers) ::  &
!   ithspr(mxoutrs) ,& ! grid pts of multi-level reports where spreading is not
!   jthspr(mxoutrs) ,& ! along model levels, but along levels according to
!                      ! 'msprpar' (active only if .NOT.lthall)
!   nthspr          ,& ! number of multi-level reports where spreading is
!                      ! along levels as specified by 'msprpar'
!                      ! (active only if .NOT. lthall)
! DATA  lthall/.TRUE./ , ithspr/ mxoutrs*0/ , jthspr/ mxoutrs*0/ , nthspr/ 0/


! Local (automatic) arrays: None
! -------------------------

  REAL    (KIND=wp   )     ::  &
                                ! --- model columns at the obs location --- :
    col_rh   (ke)            ,& ! model relative humidity                [     ]
    col_tv   (ke)            ,& ! virtual temperature                    [  K  ]
    col_t_og (ke)            ,& ! observed variable related to temperature
                                !   (temperature or virtual temperature) [  K  ]
                                !   (used only for stability-dep. humidity
                                !    QC thresholds)
    zzml     (ke)            ,& ! height
    zthml    (ke)            ,& ! potential temperature
    zlhyphl  (ke+1)          ,& ! LOG( hydrostatic pressure at half levels )
    zhhl     (ke+1)          ,& ! height at half levels
                                !
    zobbdy   (maxmlv,mxrbdy) ,& ! multi-level obs. body (format: see 'omlbdy')
    zdzob    (maxmlv)        ,& ! profile of height observation increments
    zt_o     (maxmlv)        ,& ! observation-derived (dry bulb) temperature
    zyqc     (maxmlv*4,11)   ,& ! information for QC control messages
    vimtoob  (maxmlv,mxvimo) ,& ! model values interpolated to obs. levels
    viobtom  (ke    ,mxviom)    ! obs. (incr.) interpolated to model levels
                                !   (wind: increments; T, RH: full obs values)

  INTEGER (KIND=iintegers) ::  &
    mzobbd   (maxmlv,mxrbdf) ,& ! multi-level obs. body (format: see 'momlbd')
    ilvpr    (ke)               ! obs level used for vertical interpolation
!
!------------ End of header ----------------------------------------------------


 
!-------------------------------------------------------------------------------
! Begin Subroutine mult_org_localinfo
!-------------------------------------------------------------------------------
  
!-------------------------------------------------------------------------------
!  Section 1: Preliminaries
!-------------------------------------------------------------------------------

! retreive required microphysics tracer fields
! --------------------------------------------

  izerror   = 0_iintegers
  nstat     = 0_iintegers  !US added because of warning in valgrind
  yzerrmsg  = ''
  yzroutine = 'mult_org_localinfo'

  ! retrieve the required microphysics tracers (at nnew)
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

! initialize quantities that are related to observations
! ------------------------------------------------------

  mxoiml = 0
  DO ista = 1 , nmlsta
    DO itic = 1, 2
      nmlob = imladm(ista,itic)
      IF (nmlob > 0) THEN
        IF (momlhd(nmlob,nhobtp) == nairep) THEN
          mxoiml  =  mxoiml + maxarl + ke
        ELSE
          mxoiml  =  mxoiml + maxmlv + ke
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  ! observation increment of GPS (max number as model levels)
  DO ista = 1 , ngpsta
    DO itic = 1, 2
      ngpob = igpadm(ista,itic)
      IF (ngpob > 0) THEN
          mxoiml  =  mxoiml + ke
      ENDIF
    ENDDO
  ENDDO

  ! satellite retrievals
  DO ista = 1 , ntvsta
    DO itic = 1, 2
      ntvob = itvadm(ista,itic)
      IF (ntvob > 0) THEN
          mxoiml  =  mxoiml + maxrtv + ke
      ENDIF
    ENDDO
  ENDDO

  IF (ALLOCATED( oiml )) THEN
    PRINT '("CAUTION in src_mult_local: oiml is already allocated "            &
          &,"at time ",I6)', ntstep
    DEALLOCATE ( oiml   , STAT=istat )
  ENDIF
  ALLOCATE ( oiml   (maxnoi,mxoiml) , STAT=istat )

  oiml    (:,:) = rmdi
  viobtom (:,:) = rmdi
  vimtoob (:,:) = rmdi
  zobbdy  (:,:) = rmdi
  mzobbd  (:,:) = 0
  lobinc  (1:3) = .TRUE.

! IF (msprpar <= 1) lnisua  = (MIN( vcsni(1) , vcsni(3) , vcsni(4) ) <= vcutnit)
! IF (msprpar == 2) lnisua  = (MIN( vcsni(1) , vcsni(3) , vcsni(4) ) <= vcutnip)

  svcutof (1) = SQRT( vcorls(1) * vcutof(1) )
  svcutof (2) = SQRT( vcorls(3) * vcutof(3) )
  svcutof (3) = SQRT( vcorls(4) * vcutof(4) )

  DO ni4      = 1 , 4
    ni2       = MOD( ni4-1 , 2 ) + 1
    mckvifl   = (ni2-1) * (ke+1)
    zczvcut   = (2-ni2) * zallhig  +  (ni2-1) * zalllow
    DO ivrs   = 1 , 3
      DO ista = 1 , maxmloi_tot
        zwtml   (ista,ni4,ivrs) = c0
        fcorlml (ista,ni4,ivrs) = c1
        zvcutml (ista,ni4,ivrs) = zczvcut
        zrtdgml (ista,ni4,ivrs) = c0
        kviflml (ista,ni4,ivrs) = mckvifl
        zspobml (ista,ni4,ivrs) = zczvcut
        zpobml  (ista,ni4,ivrs) = c1
        zemkpml (ista,ni4,ivrs) = c1
      ENDDO
    ENDDO
  ENDDO
  DO ivrs   = 1 , 3
    DO ista = 1 , maxmloi_tot
      zriflml (ista,1,ivrs) = c0
      zriflml (ista,2,ivrs) = c0
      ltiml   (ista  ,ivrs) = .FALSE.
      mbotlv  (ista,1,ivrs) = 0
      mbotlv  (ista,2,ivrs) = 0
      mtoplv  (ista,1,ivrs) = 0
      mtoplv  (ista,2,ivrs) = 0
    ENDDO
  ENDDO
  DO ista = 1 , maxmloi_tot
!   zsprob  (ista    ) = c0
    zstpml  (ista    ) = c0
    ksprml  (ista    ) = 0
    kkoiml  (ista,1  ) = 0
    kkoiml  (ista,2  ) = 0
    kobtyml (ista    ) = 0
    kcdtyml (ista    ) = 0
    mszlev  (ista,1  ) = 0
    mszlev  (ista,2  ) = 0
  ENDDO
  DO kk = 1 , ke
    DO ista = 1 , maxmloi_tot
      zrtgkml (ista,kk) = c0
      zsprml  (ista,kk) = c0
      zlopml  (ista,kk) = c0
      IF (lnisua) znisml  (ista,kk) = c0
    ENDDO
  ENDDO
  !   initialize also these variables as they are first used in an if statement
  DO ista = 1 , maxivq
    kiobiv (ista) = 0
    ktypiv (ista) = 0
  ENDDO

  rdg   = r_d / g
  zlop0 = LOG( p0r )
  zloptp= LOG( ptropop )

!-------------------------------------------------------------------------------
!  Section 2: Initialize the loop over multi-level data (incl. GPS) stations,
!             compute the temporal weights, and
!             determine whether a station is processed further (i.e. if it has
!             at least 1 obs. with non-zero temporal weight (incl. quality cntl)
!-------------------------------------------------------------------------------

  icml  = 0
! nmloi = 0
  ktopoi = 0
  nivtot = 0
  nexceml = 0
  nexceiv = 0

  IF (lwonl) WRITE( nupr,'(''ante loop_over_multi_level_stations'',3I4)' )     &
                    nmlsta, ngpsta, ntvsta
! IF (lwonl) PRINT       '(''ante loop_over_multi_level_stations'',3I4)' ,     &
!                   nmlsta, ngpsta, ntvsta
! IF ( ((ntstep<= 1) .OR.(ntstep== 48) .OR.(ntstep== 138) .OR.(ntstep== 228))  &
!     .AND. (lwonl)) THEN
!   PRINT '("          --> IWV: repor- re- adjusted total model model model"   &
!         &," extrapol         ")' 
!   PRINT '("    GPS-STA  hr      ted trieved obs. qv-obs   |   w.ice satur"   &
!         &," fact p-int  p-rep")' 
! ENDIF

loop_over_multi_level_stations:  DO istaml = 1 , nmlsta + ngpsta + ntvsta
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  icml = icml + 1

  lobcnv  =                                   (istaml <= nmlsta)
  lobgps  =  (istaml >  nmlsta)         .AND. (istaml <= nmlsta +ngpsta)
  lobrtv  =  (istaml >  nmlsta +ngpsta) .AND. (istaml <= nmlsta +ngpsta +ntvsta)
  IF (lobcnv)  ista  =  istaml
  IF (lobgps)  ista  =  istaml - nmlsta
  IF (lobrtv)  ista  =  istaml - nmlsta - ngpsta

  lreset         = .FALSE.

! First get temporal weights. The reason is that if (nmlob > 0), the report must
! contain data of type ivar with a non-zero nudging coefficient, and (unless all
! of these data are rejected by the threshold quality control,) it follows that
! there are data to be nudged <==> the temporal weight for these data is > zero.
! Determine also if linear temporal interpolation is applied as temporal weight,
! and if a report contains data that may have non-zero temporal weights (which
! is the condition that this report will be processed further).
! ------------------------------------------------------------------------------

! WRITE( nupr,'(3I4,A)' ) ista, imladm(ista,5), ilstidp, yomlhd(imladm(ista,5))
! ystid = yomlhd(imladm(ista,5)) (1:ilstidp)
! IF ((lwonl) .AND. (ntstep <= 1)) THEN
!   nmltot = MAX( icml , 1 )
!   WRITE( nupr,'(A ,'': MULT_LOCAL'',2I4,2L2,2I4,I3,4F5.2,4I3,2F5.2,4I3 )' )  &
!          ystid, icml, nmloi, ltakeob(1), ltakeob(2)                          &
!        , mszlev(nmltot,1), mszlev(nmltot,2), ksprml(nmltot)                  &
!        , (zwtml(nmltot,ni4,1),ni4=1,4), (kviflml(nmltot,ni4,1),ni4=1,4)      &
!        , (zwtml(nmltot,ni4,2),ni4=1,2), (kviflml(nmltot,ni4,2),ni4=1,4)
! ENDIF

! conventional multi-level data
! -----------------------------

  IF (lobcnv) THEN
    DO itic = 1, 2
      ltakeob (itic) = .FALSE.
      lqcpass (itic) = .FALSE.
      lqcfrst (itic) = .FALSE.
      lqcflag (itic) = .FALSE.
      nmlob = imladm(ista,itic)
      IF (nmlob > 0) THEN
        kobtyp = momlhd(nmlob,nhobtp)
        kcdtyp = momlhd(nmlob,nhcode)
        lraso  = (kobtyp == ntemp) .OR. (kobtyp == npilot)
        DO ivrs = 1, 3
          ivar = ivrs + MIN( ivrs-1 , 1 )
          IF (ivrs == 1) nhxexi = nhuexi
          IF (ivrs == 2) nhxexi = nhtexi
          IF (ivrs == 3) nhxexi = nhqexi
          zwtml (icml,itic  ,ivrs) = c0
          zwtml (icml,itic+2,ivrs) = c0
          IF (      (     (    (lraso           ).AND. (gnudg  (ivar) > epsy)) &
                     .OR. (    (kobtyp == nairep)                              &
                          .AND.(kcdtyp /= nmodes).AND. (gnudgar(ivar) > epsy)) &
                     .OR. (    (kobtyp == nairep)                              &
                          .AND.(kcdtyp == nmodes).AND. (gnudgms(ivar) > epsy)))&
              .AND. (momlhd(nmlob,nhxexi) > 0)                                 &
              .AND. (momlhd(nmlob,nhpass) == 0)) THEN
            zwtml (icml,itic+2,ivrs) = wtml(ista,itic+2)
            nmlo2 = imladm(ista,3-itic)
            IF ((nmlo2 == 0) .OR. (momlhd(MAX(nmlo2,i1),nhxexi) <= 0)) THEN
              zwtml (icml,itic,ivrs) = wtml(ista,itic+2)
            ELSE
              zwtml (icml,itic,ivrs) = wtml(ista,itic  )
              IF (itic == 2)  ltiml (icml,ivrs) = (imladm(ista,6) == 1)
!             IF ((itic == 2) .AND. (ltipol)) THEN
!               tabdif = ABS( omlhed(nmlob,nhtime) - omlhed(nmlo2,nhtime) )
!               ltiml (icml,ivrs) =      (ABS(tmladm(ista,2) - tabdif) < epsy) &
!                                  .AND. (tabdif <= tipolmx)
!             ENDIF
            ENDIF
            IF (MAX( zwtml(icml,itic  ,ivrs)                                   &
                   , zwtml(icml,itic+2,ivrs) ) > epsy) ltakeob (itic) = .TRUE.
          ENDIF
        ENDDO
!       ( if lqcfrst then (MOD(*nhqcfw,8) < 4) or (*nhqcfw == 0) ;
!         if (MOD(*nhqcfw,8) < 4) then lqcfrst )
        lqcflag (itic) =        (momlhd(nmlob,nhqcfw) >=  0)
        lqcfrst (itic) =        (momlhd(nmlob,nhqcfw) <= -1)                   &
                          .AND. (momlhd(nmlob,nhqcfw) > -99)                   &
                          .AND. (MOD(ABS(momlhd(nmlob,nhqcfw)),8) < 4)
        lqcpass (itic) = (.NOT. ltakeob(itic)) .AND. (lqcflag(itic))
        IF ((lqcfrst(itic)) .AND. ((ltakeob(itic)) .OR. (lqcpass(itic))))      &
          momlhd (nmlob,nhqcfw)  =  momlhd(nmlob,nhqcfw) - 4
      ENDIF
    ENDDO
  ENDIF

! GPS-derived IWV data
! --------------------
!CSC: to be adapted

  IF (lobgps) THEN
    DO itic = 1, 2
      ltakeob (itic) = .FALSE.
      lqcpass (itic) = .FALSE.
      lqcfrst (itic) = .FALSE.
      lqcflag (itic) = .FALSE.
      ngpob = igpadm(ista,itic)
      IF (ngpob > 0) THEN
        kobtyp = mogphd(ngpob,nhobtp)
        kcdtyp = mogphd(ngpob,nhcode)
        DO ivrs = 1, 3
          zwtml (icml,itic  ,ivrs) = c0
          zwtml (icml,itic+2,ivrs) = c0
          IF ((ivrs == 3) .AND. (gnudggp > epsy)                               &
!                         .AND. (ogpbdy(ngpob,nbgtze) > rmdich)                &
                          .AND. (BTEST( mogpbd(ngpob,nbgerr),nvriwv ))         &
                          .AND. (mogphd(ngpob,nhpass) == 0)) THEN
            zwtml (icml,itic+2,ivrs) = wtgp(ista,itic+2)
            ngpo2 = igpadm(ista,3-itic)
            IF ((ngpo2 == 0) .OR. (mogphd(MAX(ngpo2,i1),nhpass) > 0)) THEN
              zwtml (icml,itic,ivrs) = wtgp(ista,itic+2)
            ELSE
              zwtml (icml,itic,ivrs) = wtgp(ista,itic  )
              IF (itic == 2)  ltiml (icml,ivrs) = (igpadm(ista,6) == 1)
!             IF ((itic == 2) .AND. (ltipsu)) THEN
!               tabdif = ABS( ogphed(ngpob,nhtime) - ogphed(ngpo2,nhtime) )
!               ltiml (icml,ivrs) =      (ABS(tgpadm(ista,2) - tabdif) < epsy) &
!                                  .AND. (tabdif <= tipmxsu)
!             ENDIF
            ENDIF
            IF (MAX( zwtml(icml,itic  ,ivrs)                                   &
                   , zwtml(icml,itic+2,ivrs) ) > epsy) THEN
              ltakeob (itic) = .TRUE.
            ENDIF
          ENDIF
        ENDDO
!       ( if lqcfrst then (MOD(*nhqcfw,8) < 4) or (*nhqcfw == 0) ;
!         if (MOD(*nhqcfw,8) < 4) then lqcfrst )
        lqcflag (itic) =        (mogphd(ngpob,nhqcfw) >=  0)
        lqcfrst (itic) =        (mogphd(ngpob,nhqcfw) <= -1)                   &
                          .AND. (mogphd(ngpob,nhqcfw) > -99)                   &
                          .AND. (MOD(ABS(mogphd(ngpob,nhqcfw)),8) < 4)
!   all good GPS obs shall be used in consistency check for radiosonde RH
!   --> set 'lqcpass' not only when lqcflag, but also when lqcfrst or lqcall
        lqcpass (itic) =        (.NOT. ltakeob(itic))                          &
                          .AND. (     (lqcflag(itic)) .OR. (lqcall)            &
                                 .OR. (lqcfrst(itic)))
        IF ((lqcfrst(itic)) .AND. ((ltakeob(itic)) .OR. (lqcpass(itic))))      &
          mogphd (ngpob,nhqcfw)  =  mogphd(ngpob,nhqcfw) - 4
      ENDIF
    ENDDO
  ENDIF

! satellite retrievals
! --------------------

  IF (lobrtv) THEN
    DO itic = 1, 2
      ltakeob (itic) = .FALSE.
      lqcpass (itic) = .FALSE.
      lqcfrst (itic) = .FALSE.
      lqcflag (itic) = .FALSE.
      ntvob = itvadm(ista,itic)
      IF (ntvob > 0) THEN
        kobtyp = motvhd(ntvob,nhobtp)
        kcdtyp = motvhd(ntvob,nhcode)
        lraso  = (kobtyp == ntemp) .OR. (kobtyp == npilot)
        DO ivrs = 1, 3
          ivar = ivrs + MIN( ivrs-1 , 1 )
          IF (ivrs == 1) nhxexi = nhuexi
          IF (ivrs == 2) nhxexi = nhtexi
          IF (ivrs == 3) nhxexi = nhqexi
          zwtml (icml,itic  ,ivrs) = c0
          zwtml (icml,itic+2,ivrs) = c0
          IF ((gnudgtv(ivar) > epsy) .AND. (motvhd(ntvob,nhxexi) > 0)          &
                                     .AND. (motvhd(ntvob,nhpass) == 0)) THEN
            zwtml (icml,itic+2,ivrs) = wttv(ista,itic+2)
            ntvo2 = itvadm(ista,3-itic)
            IF ((ntvo2 == 0) .OR. (motvhd(MAX(ntvo2,i1),nhxexi) <= 0)) THEN
              zwtml (icml,itic,ivrs) = wttv(ista,itic+2)
            ELSE
              zwtml (icml,itic,ivrs) = wttv(ista,itic  )
              IF (itic == 2)  ltiml (icml,ivrs) = (itvadm(ista,6) == 1)
!             IF ((itic == 2) .AND. (ltipol)) THEN
!               tabdif = ABS( otvhed(ntvob,nhtime) - otvhed(ntvo2,nhtime) )
!               ltiml (icml,ivrs) =      (ABS(ttvadm(ista,2) - tabdif) < epsy) &
!                                  .AND. (tabdif <= tipolmx)
!             ENDIF
            ENDIF
            IF (MAX( zwtml(icml,itic  ,ivrs)                                   &
                   , zwtml(icml,itic+2,ivrs) ) > epsy) ltakeob (itic) = .TRUE.
          ENDIF
! quality control for satellite data has already been done within obs processing
!         IF (      (.NOT. ltakeob(itic)) .AND. (lverif)                       &
!             .AND. (motvhd(ntvob,nhxexi) /= 0) .AND. (lqcall)) THEN
!           IF (ABS( c3600*(otvhed(ntvob,nhtime) - acthr)+epsy ) <= cqcwbox)   &
!             lqcpass (itic) = .TRUE.
!         ENDIF
        ENDDO
!       ( if lqcfrst then (MOD(*nhqcfw,8) < 4) or (*nhqcfw == 0) ;
!         if (MOD(*nhqcfw,8) < 4) then lqcfrst )
        lqcflag (itic) =        (motvhd(ntvob,nhqcfw) >=  0)
        lqcfrst (itic) =        (motvhd(ntvob,nhqcfw) <= -1)                   &
                          .AND. (motvhd(ntvob,nhqcfw) > -99)                   &
                          .AND. (MOD(ABS(motvhd(ntvob,nhqcfw)),8) < 4)
        lqcpass (itic) = (.NOT. ltakeob(itic)) .AND. (lqcflag(itic))
        IF ((lqcfrst(itic)) .AND. ((ltakeob(itic)) .OR. (lqcpass(itic))))      &
          motvhd (ntvob,nhqcfw)  =  motvhd(ntvob,nhqcfw) - 4
      ENDIF
    ENDDO
  ENDIF


! Note:  ('zwtml' > epsy) is possible only if 'ltakeob' is TRUE.
! IF (.NOT. ((ltakeob(1)) .OR. (ltakeob(2)))) icml = icml - 1

!-------------------------------------------------------------------------------
!  Section 3: Compute some quantities which depend
!             only on the station and its location
!-------------------------------------------------------------------------------

  IF ((ltakeob(1)) .OR. (ltakeob(2)) .OR. (lqcpass(1)) .OR. (lqcpass(2))) THEN

    IF (lobcnv) THEN
      io = imladm(ista,3)
      jo = imladm(ista,4)
      nmloba = MAX( imladm(ista,1) , imladm(ista,2) )
      io_tot = momlhd(nmloba,nhitot)
      jo_tot = momlhd(nmloba,nhjtot)
      zoblat = omlhed(nmloba,nhjlat)
      ystidml (icml) = yomlhd(imladm(ista,5)) (1:ilstidg)
    ELSEIF (lobgps) THEN
      io = igpadm(ista,3)
      jo = igpadm(ista,4)
      ngpoba = MAX( igpadm(ista,1) , igpadm(ista,2) )
      io_tot = mogphd(ngpoba,nhitot)
      jo_tot = mogphd(ngpoba,nhjtot)
      zoblat = ogphed(ngpoba,nhjlat)
      ystidml (icml) = yogphd(igpadm(ista,5)) (1:ilstidg)
    ELSEIF (lobrtv) THEN
      io = itvadm(ista,3)
      jo = itvadm(ista,4)
      ntvoba = MAX( itvadm(ista,1) , itvadm(ista,2) )
      io_tot = motvhd(ntvoba,nhitot)
      jo_tot = motvhd(ntvoba,nhjtot)
      zoblat = otvhed(ntvoba,nhjlat)
      ystidml (icml) = yotvhd(itvadm(ista,5)) (1:ilstidg)
    ENDIF

! get: ioml,joml: location of the obs. station
!      ystidml: station id.
! \    zsprob : values of spreading parameter at the obs. point
!      zstpml : values of spreading parameter at the tropopause
! \    zzobml : height   at the obs. point
! \    zpobml : pressure at the obs. point
! \    zrtdgml: (Rd / g) * Tv (conversion factor ln(p) --> z) at the obs. point
!      zrtgkml: (Rd / g) * Tv (factor ln(p) --> z) at obs loc. on model levels
!      zsprml : values of spreading parameter at obs. location on model levels
!      zlopml :      log( pressure )          at obs. location on model levels
!      znisml : parameter modeling the non-isotropy at obs loc on model levels
!      znismq : parameter modeling the non-isotropy at vertically equidist. pts.
! -----------------------------------------------------------------------------

    ioml    (icml) = io
    joml    (icml) = jo
    ioml_tot(icml) = io_tot
    joml_tot(icml) = jo_tot
    kobtyml (icml) = kobtyp
    kcdtyml (icml) = kcdtyp
    ystid          = ystidml(icml) (1:ilstidp)
    DO kk = 1 , ke
      zzml        (kk) = c05 * (hhl(io,jo,kk) + hhl(io,jo,kk+1))
!     zpmoml      (kk) = a_p(io,jo,kk)
      zsprml (icml,kk) = zzml(kk)
      zrtgkml(icml,kk) = rdg*t(io,jo,kk,nnew)*( c1 + rvd_m_o*qv(io,jo,kk)      &
                                               -qc(io,jo,kk) -qrs(io,jo,kk))
! pressure on half levels: base state + vert. averaging operator applied to 'pp'
!     zphl        (kk) =  p0(io,jo,kk) - c05 *dp0(io,jo,kk)
!                       + ( dp0(io,jo,kk  ) *pp(io,jo,kk-1,nnew)
!                          +dp0(io,jo,kk-1) *pp(io,jo,kk  ,nnew))
!                        /(dp0(io,jo,kk-1) + dp0(io,jo,kk))
! LOG of (non-hydrostatic) model pressure on main levels
      zlopml (icml,kk) = LOG( a_p(io,jo,kk) )
      zthml    (kk) = t(io,jo,kk,nnew) * EXP( rdocp *(zlop0 - zlopml(icml,kk)) )
    ENDDO
!   zphl   (ke+1) = ps(io,jo,nnew)
! get 'z' at 'tropopause'
    kbt = 2
    DO WHILE ((kbt < ke) .AND. (a_p(io,jo,kbt) <= ptropop))
      kbt = kbt + 1
    ENDDO
    kat = kbt - 1
    zlptf = (zlopml(icml,kbt) - zloptp) / (zlopml(icml,kbt) - zlopml(icml,kat))
    zstpml  (icml) = (c1-zlptf) *zzml(kbt) + zlptf *zzml(kat)
! if necessary: get profile of 'theta' at model levels
    IF ((lnisua) .OR. (msprpar == 2)) THEN
      DO kk = 1 , ke
        IF                 (msprpar == 2)  zsprml (icml,kk) = zthml(kk)
        IF ((lnisua) .AND. (msprpar <= 1)) znisml (icml,kk) = zthml(kk)
        IF ((lnisua) .AND. (msprpar == 2)) znisml (icml,kk) = zzml (kk)
      ENDDO
!               get 'theta' at obs. level and at tropopause
      IF (msprpar == 2) THEN
        zstpml (icml) = (c1-zlptf) *zthml(kbt) + zlptf *zthml(kat)
      ENDIF
    ENDIF
! if necessary: get profiles of 'theta' and 'z' at 'exp(-z) -equidistant' levels
    IF ((lnisua) .AND. (msprpar >= 1)) THEN
      DO mv = 1 , mxispr
        kk = ke
        DO WHILE ( (  (msprpar <  2) .AND. (zzml (MAX(kk,i1)) <= zsdni(mv,1))  &
                 .OR. (msprpar == 2) .AND. (zthml(MAX(kk,i1)) <= zsdni(mv,2))) &
             .AND. (kk > 0))
          kk = kk - 1
        ENDDO
        kbl = MIN( INT( kk+1 ,iintegers) , ke )
        kal = MAX(      kk               , i1 )
        zflop = c1
        IF (msprpar < 2) THEN
          IF (kal /= kbl) zflop =   (zzml(kbl) - zsdni(mv,1))                  &
                                  / (zzml(kbl) - zzml(kal))
          znismq (icml,mv) = (c1-zflop) *zthml(kbl)  + zflop *zthml(kal)
        ELSE
          IF (kal /= kbl) zflop =   (zthml(kbl) - zsdni(mv,2))                 &
                                  / (zthml(kbl) - zthml(kal))
          znismq (icml,mv) = (c1-zflop) *zzml(kbl) + zflop *zzml(kal)
          IF (kal == ke)                                                       &
            znismq (icml,mv) = zzml(ke) - (zthml(ke) -zsdni(mv,2)) /thdzt
!           znismq (icml,mv) = zzml(ke) - (zthml(ke) -zsdni(mv,2)) /thdlnpt
        ENDIF
      ENDDO
      IF (lfirst) THEN
        IF ((lwonl) .AND. (io == ionl) .AND. (jo == jonl)) THEN
          mxprspr = (mxispr + 6) / 7
          DO npspr = 1, mxprspr
            nspe = MIN( INT( 7*npspr ,iintegers) , mxispr )
            nspa = nspe - 6
            IF (msprpar <= 1) THEN
              zisprb = xezispr + (nspa-1) * dezispr
              zispra = xezispr + (nspe-1) * dezispr
            ELSE
              zisprb = xthispr + (nspa-1) * dthispr
              zispra = xthispr + (nspe-1) * dthispr
            ENDIF
            WRITE( nupr,'(''spr-par['',F7.1,''..],'',A ,'' :'',7F8.1)' )       &
                   zisprb, ystid ,  (znismq(icml,mv), mv=nspa,nspe)
          ENDDO
        ENDIF
      ENDIF
    ENDIF

! Determination of lowest model level containing target grid pts.
! for which lateral spreading is along model levels
! ---------------------------------------------------------------

!   lthij  =  (msprpar > 0) .AND. (lthall)
!   IF ((msprpar > 0) .AND. (.NOT. lthall)) THEN
!     DO mthij = 1 , MIN( nthspr , maxmlo )
!       IF ((io == ithspr(mthij)) .AND. (jo == jthspr(mthij))) lthij = .TRUE.
!     ENDDO
!   ENDIF
    ksprml (icml) = ke
!   IF (lthij) THEN
      IF (msprpar == 1) ksprml (icml) = ktp
      IF (msprpar == 2) ksprml (icml) = ktth - 1
!   ENDIF

  ELSE
    ! should initialize ystidml somehow, because it is used later
    ystidml(icml) = 'not def: '
  ENDIF    ! IF ltakeob(1/2) .OR. lqcpass(1/2)

! flush YUPRINT file
  IF ((lwonl) .AND. (lfirst) .AND. (istaml == 1) .AND. (ldump_ascii)) THEN
    WRITE( nupr,'(A ,'': MULT_LOCL '',I4,4X,2L2,11X,4F5.2,4I3,2F5.2,4I3 )' )   &
           ystid, icml,        ltakeob(1), ltakeob(2)                          &
         , (zwtml(icml,ni4,1),ni4=1,4), (kviflml(icml,ni4,1),ni4=1,4)          &
         , (zwtml(icml,ni4,2),ni4=1,2), (kviflml(icml,ni4,2),ni4=1,4)
    CLOSE (nupr)
    OPEN  (nupr,FILE=yuprint,FORM='FORMATTED',STATUS='OLD'                     &
                            ,POSITION='APPEND',IOSTAT=nstat)
    IF (nstat /= 0) PRINT '("OPENING OF FILE yuprint FAILS mult_org_localinfo")'
  ENDIF


!-------------------------------------------------------------------------------
!  Section 4: Compute vertical profiles of observation increments
!             including the threshold quality control
!             for conventional multi-level data
!-------------------------------------------------------------------------------
!  Section 4.1: This part is required for nudging and LETKF / feedback files:
!               - apply the forward observation operator and do threshold
!                 quality control for individual observation and multi-level
!                 checks
!               - apply the inverse observation operator to:
!                  - wind, temperature, and humidity for nudging only
!                  - temperature for the hydrostatic height / thickness QC check
!                  - humidity for the spatial consistency check of IWV
!               - perform the hydrostatic height / thickness QC check for
!                 multi-level temperature
!-------------------------------------------------------------------------------

  IF (lobcnv) THEN

    loop_over_time_index:  DO itim = 1, 2

      IF (      (.NOT. ltakeob(itim))                                          &
          .AND. (.NOT. lqcpass(itim)))                CYCLE loop_over_time_index

      IF (ltakeob(itim)) THEN

!       nmloi = nmloi + 1
!       knoiml (icml,itim) = nmloi
        kkoiml (icml,itim) = ktopoi + 1

        !   check if threshold quality control is applied for current timestep
        !   and report
        lqc     = (lqcall) .OR. (lqcfrst(itim)) .OR. (lqcflag(itim))
        lqconly = .FALSE.

      ELSEIF (lqcpass(itim)) THEN

        lqc     = .TRUE.
        lqconly = .TRUE.
      ENDIF

      nmlob  = imladm(ista,itim)
      nlev   = momlhd(nmlob,nhnlev)
!     kobtyp = momlhd(nmlob,nhobtp)
      kcdtyp = momlhd(nmlob,nhcode)

!     lprever  =  (momlhd(nmlob,nhqcfw) >= 0)
      lveridat = (momlhd(nmlob,nhqcfw) >= 0) .AND. (lvofoi)
      timdif   =  omlhed(nmlob,nhtime) - acthr
      tabdif   =  ABS( timdif )
      lwrqc    =  (lqcall) .AND. (ABS( c3600*timdif+epsy ) <= cqcwbox)
      tdiffw   =  omlhed(nmlob,nhtime) - aiwthr - (itim-1) *tmladm(ista,2)
      !   lprint, lprmlc
      lprqcm (1) = (lwrqc) .OR. (ABS( c3600*tdiffw +epsy) < tmaxbox)
      lprqcm (1) = (lqc) .AND. (lprqcm(1)) .AND. (qcvf(4) > epsy) .AND. (lwonl)
      lprqcm (2) = (ystid(1:5) =='11520') .OR. ((io == ionl) .AND. (jo == jonl))
      lprfrs (1) = (lfirst) .AND. (lwonl)
      lprfrs (2) = (lprfrs(1)) .AND. (io == ionl) .AND. (jo == jonl)
      mexi (1)  =  momlhd(nmlob,nhuexi)
      mexi (2)  =  momlhd(nmlob,nhtexi)
      mexi (3)  =  momlhd(nmlob,nhqexi)
      zobdps    =  c0
      IF (omlhed(nmlob,nhtvip) > rmdich)                                       &
        zobdps  =  ABS( omlhed(nmlob,nhtvip) ) - a_p(io,jo,ke)
      lvirt  =  (     (kcdtyp == nra_eu) .OR. (kcdtyp == nwp_eu)               &
                 .OR. (kcdtyp == nravad) .OR. (kcdtyp == npr_us))
      nuprin = -1
      IF (lprfrs(1))  nuprin = nupr

      CALL q2rh_col ( ke, t(io,jo,:,nnew), qv(io,jo,:), qc(io,jo,:)            &
                    , a_p(io,jo,:), rdv, b1, b2w, b3, b4w, nuprin , col_rh(:) )
!     =============

!     IF ((lvirt) .OR. (lqc)) THEN

        CALL tvirt_col ( ke, t(io,jo,:,nnew), qv(io,jo,:), qc(io,jo,:)         &
                       , qrs(io,jo,:), rdv , col_tv(:) )
!       ==============

!     ENDIF
      IF (lvirt) THEN
        col_t_og (:)  =  col_tv(:)
      ELSE
        col_t_og (:)  =  t(io,jo,:,nnew)
      ENDIF
      nuprin = -1
      IF (lprfrs(2))  nuprin = nupr

! Vertical interpolation of observations to model levels and of model values to
! observation levels where required. The threshold quality control is included.
! -----------------------------------------------------------------------------

      CALL mult_obs_operator ( ke, a_u(io,jo,:), a_v(io,jo,:), col_t_og(:)     &
                             ,     col_rh(:), a_p(io,jo,:), a_z(io,jo,:)       &
                             , nlev, omlbdy(nmlob,1:nlev,:)                    &
                             ,       momlbd(nmlob,1:nlev,:)                    &
                             , kobtyp, lveridat, lqc, nuprin                   &
                             , smlbdy(:,1:nlev,nmlob) , vimtoob(1:nlev,1:5)    &
                             , lvirt, t(io,jo,:,nnew) , zt_o(1:nlev)           &
                             , kbotlev, ktoplev )
!     ======================

      IF (lqc) THEN

        CALL mult_obs_qc_fg ( nlev, omlbdy(nmlob,1:nlev,:), zt_o(1:nlev)       &
                            , kbotlev, ktoplev, kobtyp, tabdif, mexi, zoblat   &
                            , qcvf, r_d, g, lqc, lwrqc, lprqcm, nupr, ystid    &
                            , momlbd(nmlob,1:nlev,:), vimtoob(1:nlev,1:5)      &
                            , nlqc, zyqc(1:4*nlev,:), qcc )
!       ===================

        IF (lprfrs(1)) WRITE( nupr,'("Sta. ",A ,3X,I4,": kbotlev, ktoplev "    &
                                   &,2I4)' )  ystid, kbotlev, ktoplev
      ENDIF

      CALL hhl_col ( ke, a_z(io,jo,:), omlhed(nmlob,nhsurf), zhhl )
!     ============

      CALL lhyphl_col ( ke, a_p(io,jo,:), col_tv(:), hhl(io,jo,:), r_d, g      &
                      , zlhyphl )
!     ===============

      iv1 = 1
      IF (lqconly)  iv1 = 2
      iv1b = 2*iv1 - 1

      CALL mult_obs_2_modlev ( ke, a_u(io,jo,:), a_v(io,jo,:), col_t_og(:)     &
                             ,     col_rh(:), a_p(io,jo,:), zlhyphl            &
                             , nlev, omlbdy(nmlob,1:nlev,:)                    &
                             , momlbd(nmlob,1:nlev,:), vimtoob(1:nlev,1:5)     &
                             , kobtyp, nuprin, iv1, 3, lscadj(iv1b:4)          &
                             , viobtom(1:ke,iv1b:4), mbotlv(icml,itim,iv1:3)   &
                             , mtoplv(icml,itim,iv1:3), lobinc(iv1:3), ilvpr )
!     ======================

      !   control output
      IF (lprfrs(1)) THEN
        WRITE( nupr,'("sta. ",A,2X,": mbotlv, lobinc",6I4,3(2X,L1))' )         &
               ystid, (mbotlv(icml,itim,ivrs),ivrs=1,3)                        &
                    , (mtoplv(icml,itim,ivrs),ivrs=1,3), (lobinc(ivrs),ivrs=1,3)
        mv = ke - MINVAL( mbotlv )
        DO kk = mv , 1 , -1
!US  valgrind complained about too small values
!US       IF (ANY(viobtom(kk,:) < 1E-10_wp)) THEN
!US         WRITE( nupr,'(A         ,F8.0,I4)' )                                 &
!US               ' values less than 1E-10_wp    ', a_p(io,jo,kk), ilvpr(kk)
!US       ELSE
          WRITE( nupr,'(3F8.1,F6.2,F8.0,I4)' )                                 &
                (viobtom(kk,ivar), ivar=1,4), a_p(io,jo,kk), ilvpr(kk)
!US       ENDIF
        ENDDO
      ENDIF

      IF ((lqc) .OR. (lveridat)) THEN

        !   vertical correlation scales: optional input for 'mult_obs_operator_z'
        mbotlt = mbotlv(icml,itim,2)
        mtoplt = mtoplv(icml,itim,2)
        zvcut  = - c1
        zvcs   = - c1
          !   below the observed profile
        IF ((mbotlt > 0) .AND. (mbotlt+mtoplt < ke)) THEN
          IF (msprpar <= 1) zvcs (1) = SQRT( vcorls(3)) *zrtgkml(icml,ke-mbotlt)
          IF (msprpar == 2) zvcs (1) = SQRT( vcorls(3)) /thdzt
          zfcut     = (c1-fsvcut) + fsvcut *omlhed(nmlob,nhvcbt)
          zvcut (1) =  a_z(io,jo,ke-mbotlt)  - zvcs(1) *SQRT( vcutof(3) ) *zfcut
          zvcs  (1) =  zvcs(1) *omlhed(nmlob,nhvcbt)
        ENDIF
          !   above the observed profile
        IF ((mtoplt > 0) .AND. (mbotlt+mtoplt < ke)) THEN
          IF (msprpar <= 1) zvcs (2) = SQRT( vcorls(3)) *zrtgkml(icml,1 +mtoplt)
          IF (msprpar == 2) zvcs (2) = SQRT( vcorls(3)) /thdzt
          zfcut     = (c1-fsvcut) + fsvcut *omlhed(nmlob,nhvctt)
          zvcut (2) =  a_z(io,jo,1 +mtoplt)  + zvcs(2) *SQRT( vcutof(3) ) *zfcut
          zvcs  (2) =  zvcs(2) *omlhed(nmlob,nhvctt)
        ENDIF

        CALL mult_obs_operator_z ( ke, t(io,jo,:,nnew), a_p(io,jo,:)           &
                                 ,     a_z(io,jo,:), col_tv(:), hhl(io,jo,:)   &
                                 , nlev, omlbdy(nmlob,1:nlev,:)                &
                                 ,       momlbd(nmlob,1:nlev,:), zobdps        &
                                 , viobtom(1:ke,3:3), mbotlt, mtoplt           &
                                 , kobtyp, r_d, g, lveridat                    &
                                 , smlbdy(:,1:nlev,nmlob), zdzob(1:nlev)       &
                                 , lqcdz , zvcs, zvcut )
      ! ========================

        lrejtot = .FALSE.
        lredo   = .FALSE.
        IF (lqcdz)                                                             &

          CALL mult_obs_qc_dz ( nlev, omlbdy(nmlob,1:nlev,:), zdzob(1:nlev)    &
                              , tabdif, qcvf(2), r_d, g                        &
                              , momlbd(nmlob,1:nlev,:), zyqc(nlqc+1:nlqc+1,:)  &
                              , .FALSE., lrejtot, lredo )
        ! ===================

        IF (zyqc(nlqc+1,1) >= epsy)  nlqc = nlqc + 1

! Fill record for later printing for control of QC
! ------------------------------------------------
        !   adjust 'nlqc' if space in arrays 'oyqc' is insufficient
        nlqc    = MIN( maxqcp - ntotqc , nlqc )
        DO nqc = 1 , nlqc
          ntqc  =  ntotqc + nqc
          yyqc (ntqc   ) = ystid
          myqc (ntqc, 1) = momlhd(nmlob,nhcode)
          myqc (ntqc, 2) = NINT( zyqc(nqc,1) )
          IF (.NOT. lwrqc) myqc (ntqc,1) = 0
          oyqc (ntqc, 1) = omlhed(nmlob,nhtime)
          oyqc (ntqc, 2) = zyqc (nqc, 2)
          oyqc (ntqc, 3) = omlhed(nmlob,nhjlat)
          oyqc (ntqc, 4) = omlhed(nmlob,nhilon)
          oyqc (ntqc, 5) = zyqc (nqc, 5)
          oyqc (ntqc, 6) = zyqc (nqc, 6)
          oyqc (ntqc, 7) = zyqc (nqc, 7)
          oyqc (ntqc, 8) = zyqc (nqc, 8)
          oyqc (ntqc, 9) = zyqc (nqc, 9)
          oyqc (ntqc,10) = zyqc (nqc,10)
          oyqc (ntqc,11) = zyqc (nqc,11)
        ENDDO
        ntotqc  =  ntotqc  +  nlqc

!-------------------------------------------------------------------------------
!  Section 4.2: This part is required for nudging (and for QC in nudging) only:
!               - re-apply the inverse obs operator for temperature + humidity
!               - prepare the spatial consistency check of IWV
!               - fill observation increments and auxilliary quantities in
!                 arrays 'oiml', 'zwtml' etc. for the spreading of obs info
!-------------------------------------------------------------------------------

        !   epilogue of 'mult_obs_qc_dz':
        !   update 'mbotlv', 'mtoplv', 'viobtom' if total profile rejected
        IF (lrejtot) THEN
          mbotlv (icml,itim,2) = ke - mtoplv(icml,itim,2)
          mbotlv (icml,itim,3) = ke - mtoplv(icml,itim,3)
          viobtom (:,3)        = rmdi
          viobtom (:,4)        = rmdi
        ENDIF
      ELSE
        lredo = .FALSE.
      ENDIF

      IF ((lredo) .AND. (.NOT. lqconly))                                       &

        CALL mult_obs_2_modlev ( ke, a_u(io,jo,:), a_v(io,jo,:), col_t_og(:)   &
                               ,     col_rh(:), a_p(io,jo,:), zlhyphl          &
                               , nlev, omlbdy(nmlob,1:nlev,:)                  &
                               , momlbd(nmlob,1:nlev,:), vimtoob(1:nlev,1:5)   &
                               , kobtyp, -1, 2, 3, lscadj(3:4)                 &
                               , viobtom(1:ke,3:4), mbotlv(icml,itim,2:3)      &
                               , mtoplv(icml,itim,2:3), lobinc(2:3) )
      ! ======================

      IF ((lqc) .AND. (liwvssc))                                               &

        CALL mult_vprof_2_iwv ( icml, lobrtv, lqc, lqconly, lprqcm, ke         &
                              , t(io,jo,:,nnew), col_rh(:), qv(io,jo,:)        &
                              , qc(io,jo,:), qrs(io,jo,:), a_p(io,jo,:)        &
                              , a_z(io,jo,:), hhl(io,jo,:), viobtom(1:ke,1:4)  &
                              , nexceiv )
      ! =====================

! Redo interpolation of temperature and humidity observations to model levels,
! if part of profile has been rejected by height / thickness quality control

      IF (ltakeob(itim)) THEN

! Merging of data at observation and model levels into one profile of
! observation increments (full observed and model values only for moisture)
! -------------------------------------------------------------------------

        CALL mult_obs_increment ( icml, ktopoi, kcdtyp, ke, t(io,jo,:,nnew)    &
                                , qv(io,jo,:), a_p(io,jo,:), a_z(io,jo,:)      &
                                , qc(io,jo,:), qrs(io,jo,:), zthml(:)          &
                                , nlev, omlbdy(nmlob,1:nlev,:)                 &
                                ,       momlbd(nmlob,1:nlev,:)                 &
                                ,       vimtoob(1:nlev,:), viobtom, lobinc     &
                                , numlev )
      ! =======================

! If there are no active levels of type 'ivrs', set temporal weights to zero
        DO ivrs = 1 , 3
          IF (    ((.NOT.levdiff) .AND. ( mbotlv(icml,itim,ivrs)               &
                                         +mtoplv(icml,itim,ivrs) >= ke))       &
              .OR.((     levdiff) .AND. (zpobml(icml,2*itim-1,ivrs) < c2))) THEN
            zwtml (icml,itim  ,ivrs) = c0
            zwtml (icml,itim+2,ivrs) = c0
          ENDIF
        ENDDO

! Number (amount) of obs. increment levels
        mszlev (icml,itim) = numlev
        ktopoi = ktopoi + numlev

! Without any active levels, the current report is not processed any further
        IF (numlev == 0) ltakeob (itim) = .FALSE.

      ENDIF

    ENDDO loop_over_time_index

! Note:  ('zwtml' > epsy) is possible only if 'ltakeob' is TRUE.

    IF (.NOT. ((ltakeob(1)) .OR. (ltakeob(2))))  lreset = .TRUE.

  ENDIF


!-------------------------------------------------------------------------------
!  Section 5: Compute vertical profiles of GPS observation increments
!             including the threshold quality control
!             for GPS data
!-------------------------------------------------------------------------------

  IF (lobgps) THEN

    loop_over_time_gps:  DO itim = 1, 2

      IF (ltakeob(itim)) THEN

!       nmloi = nmloi + 1
!       knoiml (icml,itim) = nmloi
        kkoiml (icml,itim) = ktopoi + 1

        ngpob = igpadm(ista,itim)

! check if threshold quality control is applied for current timestep and report

        lqc     = (lqcall) .OR. (lqcfrst(itim)) .OR. (lqcflag(itim))
        lqconly = .FALSE.

! Computation of humidity observation increment profile from GPS integrated
! water vapour on model levels. The threshold quality control is included.
! -------------------------------------------------------------------------

        CALL mult_obs_gps ( icml , lqc , lqconly , ktopoi , ke                 &
                          , t(io,jo,:,nnew), qv(io,jo,:), a_p(io,jo,:)         &
                          , a_z(io,jo,:), qc(io,jo,:), qrs(io,jo,:)            &
                          , zthml(:), hhl(io,jo,:) , nexceiv , numlev )
!       =================

! Number (amount) of obs. increment levels
        mszlev (icml,itim) = numlev
        ktopoi = ktopoi + numlev

! Without any active levels, the current report is not processed any further
        IF (numlev == 0) ltakeob (itim) = .FALSE.

      ELSEIF (lqcpass(itim)) THEN

        lqc     = .TRUE.
        lqconly = .TRUE.

! Computation of humidity observation increment profile from GPS integrated
! water vapour on model levels. The threshold quality control is included.
! -------------------------------------------------------------------------

        CALL mult_obs_gps ( icml , lqc , lqconly , ktopoi , ke                 &
                          , t(io,jo,:,nnew), qv(io,jo,:), a_p(io,jo,:)         &
                          , a_z(io,jo,:), qc(io,jo,:), qrs(io,jo,:)            &
                          , zthml(:), hhl(io,jo,:) , nexceiv , numlev )
!       =================

      ENDIF

    ENDDO loop_over_time_gps

! Note:  ('zwtml' > epsy) is possible only if 'ltakeob' is TRUE.

    IF (.NOT. ((ltakeob(1)) .OR. (ltakeob(2))))  lreset = .TRUE.

!   PRINT * , 'GPSpost ', ntstep, lreset, ltakeob, lqcpass, lqcfrst, lqcflag   &
!                       , ystidml(icml), igpadm(ista,1), igpadm(ista,2)        &
!                 , NINT( ogphed(MAX(1,igpadm(ista,1)),nhtime) *100._wp )  &
!                 , NINT( ogphed(MAX(1,igpadm(ista,2)),nhtime) *100._wp )  &
!                       , mogpbd(MAX(1,igpadm(ista,1)),nbgqcf)                 &
!                       , mogpbd(MAX(1,igpadm(ista,2)),nbgqcf)                 &
!                       , mogpbd(MAX(1,igpadm(ista,1)),nbgerr)                 &
!                       , mogpbd(MAX(1,igpadm(ista,2)),nbgerr)  

  ENDIF


!-------------------------------------------------------------------------------
!  Section 6: Compute vertical profiles of observation increments
!             for satellite retrievals
!             without (!) the threshold quality control (TQC) since the TQC
!             has already been accomplished within the observation processing
!-------------------------------------------------------------------------------

  IF (lobrtv) THEN

    loop_over_time_retv:  DO itim = 1, 2

      IF (      (.NOT. ltakeob(itim))                                          &
          .AND. (.NOT. lqcpass(itim)))                 CYCLE loop_over_time_retv

      IF (ltakeob(itim)) THEN

!       nmloi = nmloi + 1
!       knoiml (icml,itim) = nmloi
        kkoiml (icml,itim) = ktopoi + 1
        ntvob = itvadm(ista,itim)

        !   check if threshold quality control is applied for current timestep
        !   and report
        lqc     = (lqcall) .OR. (lqcfrst(itim)) .OR. (lqcflag(itim))
        lqconly = .FALSE.

      ELSEIF (lqcpass(itim)) THEN

        lqc     = .TRUE.
        lqconly = .TRUE.
      ENDIF

!     IF ((ltakeob(itim)) .OR. (lqcpass(itim))) THEN
        ntvob  = itvadm(ista,itim)
        nlev   = motvhd(ntvob,nhnlev)
        kcdtyp = motvhd(ntvob,nhcode)

        lveridat = (motvhd(ntvob,nhqcfw) >= 0) .AND. (lvofoi)
        timdif   =  otvhed(ntvob,nhtime) - acthr
        tabdif   =  ABS( timdif )
        lwrqc    =  (lqcall) .AND. (ABS( c3600*timdif+epsy ) <= cqcwbox)
        lprqcm (:) = .FALSE.
        lprfrs (1) = (lfirst) .AND. (lwonl)
        lprfrs (2) = (lprfrs(1)) .AND. (io == ionl) .AND. (jo == jonl)
        zobdps     =  c0
        lvirt      =  .FALSE.

        ztexi  = rmdi
        zqexi  = rmdi
        IF (motvhd(ntvob,nhtexi) >= 1)  ztexi = c1
        IF (motvhd(ntvob,nhqexi) >= 1)  zqexi = c1
        IF ((MOD( motvhd(ntvob,nhvqcf),4 ) >= 2) .AND. (ztexi > c0)) ztexi = -c1
        IF ((MOD( motvhd(ntvob,nhvqcf),8 ) >= 4) .AND. (zqexi > c0)) zqexi = -c1
        DO ilev = 1 , nlev
          zobbdy (ilev,nbtt  )  =  otvbdy(ntvob,ilev,nbvt  )
          zobbdy (ilev,nbtrh )  =  otvbdy(ntvob,ilev,nbvrh )
          zobbdy (ilev,nbtp  )  =  otvbdy(ntvob,ilev,nbvp  )
          zobbdy (ilev,nbtlop)  =  LOG( otvbdy(ntvob,ilev,nbvp) )
          zobbdy (ilev,nbtter)  =  ztexi
          zobbdy (ilev,nbtqer)  =  zqexi
          zobbdy (ilev,nbtzio)  =  otvhed(ntvob,nhzio)
          zobbdy (ilev,nbtzjo)  =  otvhed(ntvob,nhzjo)
          zobbdy (ilev,nbtu  )  =  rmdi
          zobbdy (ilev,nbtv  )  =  rmdi
          zobbdy (ilev,nbtz  )  =  rmdi
          zobbdy (ilev,nbtuer)  =  rmdi
          zobbdy (ilev,nbtzer)  =  rmdi
          mzobbd (ilev,nbtflg)  =  motvbd(ntvob,ilev,nbvflg)
          mzobbd (ilev,nbterr)  =  motvbd(ntvob,ilev,nbverr)
          mzobbd (ilev,nbtqcf)  =  0
          mzobbd (ilev,nbtlsg)  =  0
          mzobbd (ilev,nbtlid)  =  0
        ENDDO
!     ENDIF
      nuprin = -1
      IF (lprfrs(1))  nuprin = nupr

      CALL q2rh_col ( ke, t(io,jo,:,nnew), qv(io,jo,:), qc(io,jo,:)              &
                    , a_p(io,jo,:), rdv, b1, b2w, b3, b4w, nuprin , col_rh(:) )
!     =============

      col_t_og (:)  =  t(io,jo,:,nnew)
      nuprin = -1
      IF (lprfrs(2))  nuprin = nupr

! Vertical interpolation of observations to model levels and of model values to
! observation levels where required. The threshold quality control is included.
! -----------------------------------------------------------------------------

      CALL mult_obs_operator ( ke, a_u(io,jo,:), a_v(io,jo,:), col_t_og(:)     &
                             ,     col_rh(:), a_p(io,jo,:), a_z(io,jo,:)       &
                             , nlev, zobbdy(1:nlev,:), mzobbd(1:nlev,:)        &
                             , kobtyp, lveridat, .FALSE., nuprin               &
                             , stvbdy(:,1:nlev,ntvob) , vimtoob(1:nlev,1:5) )
!     ======================

      CALL tvirt_col ( ke, t(io,jo,:,nnew), qv(io,jo,:), qc(io,jo,:)           &
                     , qrs(io,jo,:), rdv , col_tv(:) )
!     ==============

      CALL lhyphl_col ( ke, a_p(io,jo,:), col_tv(:), hhl(io,jo,:), r_d, g      &
                      , zlhyphl )
!     ===============

      CALL mult_obs_2_modlev ( ke, a_u(io,jo,:), a_v(io,jo,:), col_t_og(:)     &
                             ,     col_rh(:), a_p(io,jo,:), zlhyphl            &
                             , nlev, zobbdy(1:nlev,:)                          &
                             ,       mzobbd(1:nlev,:), vimtoob(1:nlev,1:5)     &
                             , kobtyp, nuprin, 2, 3, lscadj(3:4)               &
                             , viobtom(1:ke,3:4), mbotlv(icml,itim,2:3)        &
                             , mtoplv(icml,itim,2:3), lobinc(2:3) )
!     ======================

      IF (lqc) THEN

        !   vertical correlation scales:  optional input for 'mult_obs_qc_dz'
        mbotlt = mbotlv(icml,itim,2)
        mtoplt = mtoplv(icml,itim,2)
        zvcut  = - c1
        zvcs   = - c1
        IF ((mbotlt > 0) .AND. (mbotlt+mtoplt < ke)) THEN
          IF (msprpar <= 1) zvcs (1) = SQRT( vcorls(3)) *zrtgkml(icml,ke-mbotlt)
          IF (msprpar == 2) zvcs (1) = SQRT( vcorls(3)) /thdzt
          zvcut (1)  =  a_z(io,jo,ke-mbotlt) -  zvcs(1) *SQRT( vcutof(3) )
        ENDIF
        IF ((mtoplt > 0) .AND. (mbotlt+mtoplt < ke)) THEN
          IF (msprpar <= 1) zvcs (2) = SQRT( vcorls(3)) *zrtgkml(icml,1 +mtoplt)
          IF (msprpar == 2) zvcs (2) = SQRT( vcorls(3)) /thdzt
          zvcut (2)  =  a_z(io,jo,1+ mtoplt)  +  zvcs(2) *SQRT( vcutof(3) )
        ENDIF

        CALL mult_obs_operator_z ( ke, t(io,jo,:,nnew), a_p(io,jo,:)           &
                                 ,     a_z(io,jo,:), col_tv(:), hhl(io,jo,:)   &
                                 , nlev, zobbdy(1:nlev,:), mzobbd(1:nlev,:)    &
                                 , zobdps, viobtom(1:ke,3:3), mbotlt, mtoplt   &
                                 , kobtyp, r_d, g, .FALSE.                     &
                                 , stvbdy(:,1:nlev,ntvob), zdzob(1:nlev)       &
                                 , lqcdz , zvcs, zvcut )
      ! ========================

        lrejtot = .FALSE.
        IF (lqcdz)                                                             &

          CALL mult_obs_qc_dz ( nlev, zobbdy(1:nlev,:), zdzob(1:nlev)          &
                              , tabdif, qcvf(2), r_d, g                        &
                              , mzobbd(1:nlev,:), zyqc(nlqc+1:nlqc+1,:)        &
                              , .TRUE., lrejtot )
        ! ===================

        IF (zyqc(nlqc+1,1) >= epsy)  nlqc = nlqc + 1

! Fill record for later printing for control of QC
! ------------------------------------------------
        !   adjust 'nlqc' if space in arrays 'oyqc' is insufficient
        nlqc    = MIN( maxqcp - ntotqc , nlqc )
        DO nqc = 1 , nlqc
          ntqc  =  ntotqc + nqc
          yyqc (ntqc   ) = ystid
          myqc (ntqc, 1) = motvhd(ntvob,nhcode)
          myqc (ntqc, 2) = NINT( zyqc(nqc,1) )
          IF (.NOT. lwrqc) myqc (ntqc,1) = 0
          oyqc (ntqc, 1) = otvhed(ntvob,nhtime)
          oyqc (ntqc, 2) = zyqc (nqc, 2)
          oyqc (ntqc, 3) = otvhed(ntvob,nhjlat)
          oyqc (ntqc, 4) = otvhed(ntvob,nhilon)
          oyqc (ntqc, 5) = zyqc (nqc, 5)
          oyqc (ntqc, 6) = zyqc (nqc, 6)
          oyqc (ntqc, 7) = zyqc (nqc, 7)
          oyqc (ntqc, 8) = zyqc (nqc, 8)
          oyqc (ntqc, 9) = zyqc (nqc, 9)
          oyqc (ntqc,10) = zyqc (nqc,10)
          oyqc (ntqc,11) = zyqc (nqc,11)
        ENDDO
        ntotqc  =  ntotqc  +  nlqc

!-------------------------------------------------------------------------------
!  Section 6.2: This part is required for nudging (and for QC in nudging) only:
!               fill observation increments and auxilliary quantities in
!               arrays 'oiml', 'zwtml' etc. for the spreading of obs info
!-------------------------------------------------------------------------------

        !   epilogue of 'mult_obs_qc_dz':
        !   update 'mbotlv', 'mtoplv', 'viobtom' if total profile rejected
        IF (lrejtot) THEN
          mbotlv (icml,itim,2) = ke - mtoplv(icml,itim,2)
          mbotlv (icml,itim,3) = ke - mtoplv(icml,itim,3)
          viobtom (:,3)        = rmdi
          viobtom (:,4)        = rmdi
        ENDIF
      ENDIF

      IF (ltakeob(itim)) THEN

! Merging of data at observation and model levels into one profile of
! observation increments (full observed and model values only for moisture)
! -------------------------------------------------------------------------

        CALL mult_obs_increment ( icml, ktopoi, kcdtyp, ke, t(io,jo,:,nnew)    &
                                , qv(io,jo,:), a_p(io,jo,:), a_z(io,jo,:)      &
                                , qc(io,jo,:), qrs(io,jo,:), zthml(:)          &
                                , nlev, zobbdy(1:nlev,:), mzobbd(1:nlev,:)     &
                                , vimtoob(1:nlev,:), viobtom, lobinc , numlev )
      ! =======================

! If there are no active levels of type 'ivrs', set temporal weights to zero
        DO ivrs = 1 , 3
          IF (    ((.NOT.levdiff) .AND. ( mbotlv(icml,itim,ivrs)               &
                                         +mtoplv(icml,itim,ivrs) >= ke))       &
              .OR.((     levdiff) .AND. (zpobml(icml,2*itim-1,ivrs) < c2))) THEN
            zwtml (icml,itim  ,ivrs) = c0
            zwtml (icml,itim+2,ivrs) = c0
          ENDIF
        ENDDO

! Number (amount) of obs. increment levels
        mszlev (icml,itim) = numlev
        ktopoi = ktopoi + numlev

! Without any active levels, the current report is not processed any further
        IF (numlev == 0) ltakeob (itim) = .FALSE.

      ENDIF

!CS: maybe to be modified again
      IF ((ltakeob(itim)) .OR. (lqcpass(itim))) THEN
        DO ilev = 1 , nlev
          motvhd (ntvob,nhvqcf) = MAX( motvhd(ntvob,nhvqcf)                    &
                                     , mzobbd(ilev,nbtqcf) )
!         motvbd (ntvob,ilev,nbtqcf) = mzobbd(ilev,nbtqcf)
        ENDDO
! If there are no active levels of type 'ivrs', set temporal weights to zero
        DO ivrs = 1 , 3
          IF (motvhd(ntvob,nhvqcf) >= 1) THEN
            zwtml (icml,itim  ,ivrs) = c0
            zwtml (icml,itim+2,ivrs) = c0
          ENDIF
        ENDDO
      ENDIF

    ENDDO loop_over_time_retv

! Note:  ('zwtml' > epsy) is possible only if 'ltakeob' is TRUE.

    IF (.NOT. ((ltakeob(1)) .OR. (ltakeob(2))))  lreset = .TRUE.

  ENDIF


!-------------------------------------------------------------------------------
!  Section 7: Reset local information array if station is not taken
!-------------------------------------------------------------------------------

  ! check array bound
  IF ((icml >= maxmloi_tot) .AND. (.NOT. lreset)) THEN
    ! use only maxmloi_tot-1 obs incr reports for correct counting of nexceml;
    ! note: usually local number 'icml' is << than global maximum 'maxmloi_tot'
    nexceml = nexceml + 1
    lreset = .FALSE.
  ENDIF

  IF (lreset) THEN
    ltakeob (1) = .FALSE.
    ltakeob (2) = .FALSE.
    DO ivrs   = 1 , 3
      zwtml   (icml,1,ivrs) = c0
      zwtml   (icml,2,ivrs) = c0
      zwtml   (icml,3,ivrs) = c0
      zwtml   (icml,4,ivrs) = c0
      ltiml   (icml  ,ivrs) = .FALSE.
      mbotlv  (icml,1,ivrs) = 0
      mbotlv  (icml,2,ivrs) = 0
      mtoplv  (icml,1,ivrs) = 0
      mtoplv  (icml,2,ivrs) = 0
    ENDDO
    kkoiml  (icml,1  ) = 0
    kkoiml  (icml,2  ) = 0

    IF (      ((lobgps) .OR. (lobcnv)) .AND. (liwvssc)                         &
        .AND. (nivtot >= 1) .AND. (kioiiv(MAX(nivtot,1)) == icml)) THEN
      kioiiv (nivtot) = 0
    ENDIF

    icml  =  icml - 1
  ENDIF


!-------------------------------------------------------------------------------
!  Section 8: Get the vertical upper & lower cut-offs for a report,
!             an upper estimate for the range of model levels influenced by a
!             report (except for the lower limit with horizontal spreading), & 
!             upper estimates for horizontal correlation scales for target grid
!             pts. at the vertical cut-offs
!-------------------------------------------------------------------------------

! get range of model levels within the vertical cut-off at the obs. location
! --------------------------------------------------------------------------

  itica = 2
  itice = 1
  IF (ltakeob(1)) itica = 1
  IF (ltakeob(2)) itice = 2
  ivrsa = 1
  ivrse = 3
  IF (lobgps) ivrsa = 3
  IF (lobrtv) ivrsa = 2

! start loop over variables

  loop_over_variables:  DO ivrs = ivrsa , ivrse

    IF ((ltakeob(1)) .OR. (ltakeob(2))) THEN

      ivar = ivrs + MIN( ivrs-1 , 1 )
      IF (lobcnv) THEN
        IF (ivrs == 1)  nhvcbx = nhvcbu
        IF (ivrs == 2)  nhvcbx = nhvcbt
        IF (ivrs == 3)  nhvcbx = nhvcbq
        IF (ivrs == 1)  nhvctx = nhvctu
        IF (ivrs == 2)  nhvctx = nhvctt
        IF (ivrs == 3)  nhvctx = nhvctq
      ENDIF

! get vertical cut-offs (in spreading parameter units)

! (nbotop == 1) for lower cut-off, (nbotop == 2) for upper cut-off
      DO nbotop = 1 , 2
        DO itic = itica , itice
          nba = nbotop + 2*(itic-1)
          psign = REAL ( 2 *nbotop - 3, wp )
          IF (lobcnv) THEN
            nmloba = MAX( imladm(ista,1) , imladm(ista,2) )
            IF (nbotop == 1) fcorlml (icml,nba,ivrs) = omlhed(nmloba,nhvcbx)
            IF (nbotop == 2) fcorlml (icml,nba,ivrs) = omlhed(nmloba,nhvctx)
          ELSE
            IF (nbotop == 1) fcorlml (icml,nba,ivrs) = c1
            IF (nbotop == 2) fcorlml (icml,nba,ivrs) = c1
          ENDIF
          svcutos = svcutof(ivrs) * ( (c1-fsvcut)                              &
                                     +    fsvcut * fcorlml(icml,nba,ivrs) )
          IF (msprpar <= 1) THEN
            zcut = zspobml(icml,nba,ivrs) +psign*zrtdgml(icml,nba,ivrs) *svcutos
          ELSE
            zsob = zspobml(icml,nba,ivrs)
            IF (zsob <= zstpml(icml)) THEN
              zcut  =  zsob         +   psign *       svcutos
              IF (zcut > zstpml(icml)+epsy)                                    &
                zcut = zstpml(icml) + SQRT( (         svcutos   )**2           &
                                           -(zstpml(icml) - zsob)**2) *qst
            ELSE
              zcut  =  zsob         +   psign * qst * svcutos
              IF (zcut < zstpml(icml)-epsy)                                    &
                zcut = zstpml(icml) - SQRT( (   qst * svcutos   )**2           &
                                           -(zstpml(icml) - zsob)**2) /qst
            ENDIF
          ENDIF
          zvcutml (icml,nba,ivrs) = zcut

! get the range of influenced vertical levels
!         kthr = ke - mtoplv(icml,...)
          kthr = ke
          DO WHILE ((kthr > 1) .AND. (zsprml(icml,kthr) < zcut))
            kthr = kthr - 1
          ENDDO
          IF ((nbotop == 2) .AND. (zsprml(icml,1) > zcut)) kthr = kthr + 1
          IF ((nbotop == 1) .AND. (kthr == ke))   kthr = ke + 1
          IF ((nbotop == 2) .AND. (msprpar >= 1)) kthr = MIN( kthr , ke )
! If this range ends at a model level, where obs. increments are spreaded along
! isentropes, the range of levels influenced by the obs. may be anywhere between
! the uppermost and lowermost level where spreading is along isentropes.
          IF ((msprpar == 2) .AND. (kthr >= ktth) .AND. (kthr <= ke)) THEN
            IF (nbotop == 1) kthr = ke
            IF (nbotop == 2) kthr = ktth
          ENDIF
! if (msprpar == 1) :
! Wanted: upper limits to the 2 horizontal radii of influence used to determine
!         (a lower limit of) the range of model levels influenced by the obs.
! -----------------------------------------------------------------------------
!     Needed: the (lower limit for the) pressure at the obs. location of the
!             uppermost model level, which intersects height 'zcut' within
!             the horizontal area of influence.
!     ==> 1.: the height difference between this uppermost level and 'zcut'
!             must not be larger than the orography at the obs. location
!             ==> 'zupper' = 'zcut' + orography at the obs. location
!    !        (assumption: no orographic depression < - 400 m )
!         2.: determ. of pressure at the obs. location of the model level which
!             is immediately above (nbotop == 1) or below (nbotop == 2) 'zupper'
!       IF ('zcut' < orography) THEN lower level limit is known & set to 'ke+1'.
          IF ((msprpar == 1) .AND. (kthr > ktp) .AND. (kthr <= ke)) THEN
            zupper  = zcut + hhl(io,jo,ke+1) + 400.0_wp
            DO WHILE ((kthr >  1) .AND. (zsprml(icml,kthr) < zupper))
              kthr = kthr - 1
            ENDDO
            IF ((nbotop == 2) .AND. (zsprml(icml,1) > zupper)) kthr = kthr + 1
            kk = MIN( kthr , ke )
            IF ((lfirst) .AND. (lwonl))                                        &
              WRITE( nupr,'(''KTHR: '', A , 3I3, 3F7.0, I4)' ) ystid, icml,nba &
                   , ivrs, hhl(io,jo,ke+1), zcut, zsprml(icml,kk), kthr
            kthr = - MIN( kthr , INT( ktp + 1 ,iintegers) )
          ENDIF
          kviflml (icml,nba,ivrs) = kthr

! upper limits for horiz. correl. scales for grid pts. at the vertical cut-offs

          zrwt   = MAX( zwtml(icml,1,ivrs) , zwtml(icml,2,ivrs)                &
                      , zwtml(icml,3,ivrs) , zwtml(icml,4,ivrs) )
          ! note: the largest correlation scale is for smallest zwtml > epsy
          IF (zwtml(icml,1,ivrs) > epsy)  zrwt = MIN( zrwt, zwtml(icml,1,ivrs) )
          IF (zwtml(icml,2,ivrs) > epsy)  zrwt = MIN( zrwt, zwtml(icml,2,ivrs) )
          IF (zwtml(icml,3,ivrs) > epsy)  zrwt = MIN( zrwt, zwtml(icml,3,ivrs) )
          IF (zwtml(icml,4,ivrs) > epsy)  zrwt = MIN( zrwt, zwtml(icml,4,ivrs) )
!         IF (      (.NOT. ltiml(icml,ivrs))                                   &
!             .AND. (MIN( zwtml(icml,1,ivrs) , zwtml(icml,2,ivrs) ) > epsy))   &
!           zrwt = MIN( zwtml(icml,1,ivrs) , zwtml(icml,2,ivrs) )
          IF ((kthr < 0) .OR. (nbotop == 2)) THEN
!           zlopco = LOG( MAX( a_p(io,jo,MIN(ABS(kthr),ke)) , 1000.0_wp ) )
            zlopco = MAX( LOG( a_p(io,jo,MIN(ABS(kthr),ke)) )                  &
                        , tabcolp(ncolev) )
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
            zrtinfl = MAX(  (c1 + (rhtfac(ivar)-c1) *(c1-zrwt))                &
                          * (rhinfl(ivar) + rhvfac(ivar) *rvinfl) , epsy )
            IF (lobgps)  zrtinfl = zrtinfl *rhfgps
            IF (lobrtv)  zrtinfl = zrtinfl *rhfrtv(ivar)

            zriflml (icml,nbotop,ivrs) = MAX( zriflml(icml,nbotop,ivrs)        &
                                            , zrtinfl )
          ENDIF
        ENDDO
      ENDDO

    ENDIF

  ENDDO loop_over_variables

  IF ((lwonl) .AND. (ntstep <= 1)) THEN
    nmltot = MAX( icml , i1 )
    ystid = ystidml(nmltot) (1:ilstidp)
!   IF (lobcnv)                                                                &
    IF ((lobcnv) .OR. (lobrtv))                                                &
      WRITE( nupr,'(A ,'': MULT_LOCAL'',I4,I6,2L2,3I3,4F5.2,4I3,2F5.2,4I3 )' ) &
             ystid, icml, ktopoi, ltakeob(1), ltakeob(2)                       &
           , mszlev(nmltot,1), mszlev(nmltot,2), ksprml(nmltot)                &
           , (zwtml(nmltot,ni4,1),ni4=1,4), (kviflml(nmltot,ni4,1),ni4=1,4)    &
           , (zwtml(nmltot,ni4,2),ni4=1,2), (kviflml(nmltot,ni4,2),ni4=1,4)
    IF (lobgps)                                                                &
      WRITE( nupr,'(A ,'': MULT_LOCAL'',I4,I6,2L2,3I3,4F5.2,4I3,'' GPS '' )' ) &
             ystid, icml, ktopoi, ltakeob(1), ltakeob(2)                       &
           , mszlev(nmltot,1), mszlev(nmltot,2), ksprml(nmltot)                &
           , (zwtml(nmltot,ni4,3),ni4=1,4), (kviflml(nmltot,ni4,3),ni4=1,4)
  ENDIF


ENDDO loop_over_multi_level_stations
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  nmltot = icml
! nmloit = nmloi
  nmloit = ktopoi

! flush YUPRINT file
  IF (lwonl .AND. ldump_ascii) THEN
    CLOSE (nupr)
    OPEN  (nupr,FILE=yuprint,FORM='FORMATTED',STATUS='OLD'                     &
                            ,POSITION='APPEND',IOSTAT=nstat)
    IF (nstat /= 0) PRINT '("OPENING OF FILE yuprint FAIL, mult_org_localinfo")'
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure mult_org_localinfo
!-------------------------------------------------------------------------------

END SUBROUTINE mult_org_localinfo


!-------------------------------------------------------------------------------
!+ Module procedure for spatial consistency check for IWV (radiosonde and GPS)
!-------------------------------------------------------------------------------

SUBROUTINE mult_iwv_horicheck

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure performs an additional quality control step for
!   radiosonde humidity and GPS IWV observations by checking the spatial
!   consistency of IWV derived from these observations.
!
! Method:
!   At the location of each IWV 'observation', a weighted sum of spreaded 
!   observation increments from the neighbouring IWV observations is computed.
!   This sum is used as an 'analysis increment' in the sense that it is added 
!   (only within this check) to the model value to obtain a improved estimate
!   of truth, which is then used in a (first guess type) threshold check of IWV.
!   The threshold for this check is reduced depending on the weights used to
!   compute the 'analysis increment' (the reduction is small if the observations
!   have only small weights at the location of the checked station), and then 
!   enhanced by a fraction of the analysis increment.
!   (==> The smallest thresholds will occur in data-dense areas with small
!        observation increments at the neighbouring stations.)
!   Formal remarks:
!   The following quantities are modified in this procedure:
!   - the observation error in the (local) ODR, used for the next times when
!     real analysis increments are computed (before doing threshold quality
!     control again)
!   - the quality control flag for the (local) VOF, used for verification only
!   - the humidity-related local information (e.g. in 'oiml','zwtml','kviflml').
!   Note, that these output quantities are not used as input for any computa-
!   tions in this procedure, thus yielding results independent from the domain
!   decomposition.
!
! Known shortcoming:
!   Result may (theoretically) depend on the switch 'lverif'.
!
! Written by        :  Christoph Schraff, DWD  (original version: 11.05.05)
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

  INTEGER (KIND=iintegers) , PARAMETER  ::  &
    mwtypqc  =  2       ! type of multiple weighting of observation systems
                        ! for the ad-hoc analysis increment in the spatial
                        ! consistency check for IWV:
                        ! = 0 : net increments from all obs types at once
                        ! = 1 : separate net increments for radiosondes
                        ! = 2 : separate net increments for raso and GPS each

  REAL    (KIND=wp   )     , PARAMETER  ::  &
    fmultw   =  c0      ! 0. or 1., depending on multiple observation weighting

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    istaiv , istaot  ,& ! loop indices of observing station
    istaqc           ,& ! sorted index of observing station
    icml             ,& ! index of obs. sta. in data record to be used further
    ii     , jj      ,& ! horizontal indices of obs. being checked
    io     , jo      ,& ! horizontal indices of obs. used for checking
    idix   , jdix    ,& ! indices depending on the distance from the obs.
    itic             ,& ! loop index over time
    itica  , itice   ,& ! loop range over time for station being checked
    ityp   , ntyp    ,& ! index for observation type
    ilev             ,& ! index for observation level
    ilbot  , iltop   ,& ! range of observation levels with humidity obs
    klvi             ,& ! vertical index in 'oiml'
    ngpob  , nmlob   ,& ! indices of reports in the ODR (obs data record)
    ntvob            ,& ! indices of reports in the ODR (obs data record)
    ngpob2 , nmlob2  ,& ! indices of other reports at same station
    ntvob2           ,& ! indices of other reports at same station
    nqexi  (2)       ,& ! if > 0: humidity obs exists and needs to be checked
    nredowt          ,& ! if > 0: re-compute temporal weights
    mflgqc              ! QC bit flag

  REAL    (KIND=wp   )     ::  &
    ziwvwi (3)       ,& ! sum of weighted observat. increments
    omyiv  (3)       ,& ! sum of 'spreading' weights
    om2iv  (3)       ,& ! sum of square of 'spreading' weights
    zwsum  (4)       ,& ! sum of weights assigned to the different obs. types
    omytp  (3)       ,& ! individual weight assigned to the obs type 'ityp'
    zdist  , zdist2  ,& ! (square of) distance from the checked obs. location
    omysc  , omysc2  ,& ! horizontal nudging weight (for 1 obs) ; omysc **2
    zsprcor          ,& ! spreading correction factor to ps- obs. increment
    r1iflp , r1iflpk ,& ! horizontal correlation scale
    zcutop , zcutop2 ,& ! (square of) cut-off radius for area of influence
    dtchkiq          ,& ! temporal radius of influence for spatial check
!   zcutmlf          ,& ! correction to 'zcutop2' for ps-obs from TEMPs
    zhcut2           ,& ! square of cut-off radius for area of influence
    zreftim          ,& ! time of obs. being checked
    zddsp            ,& ! 'zdist' scaled by horizontal correlation scale
    zincr1 , zincr2  ,& ! obs. incr. of pot. temperature  (for 2 obs at 1 sta.)
    zincrt , zinct2  ,& ! sum of the 2 weighted (**2) obs. incr. of temperature
    omyk1  , omyk2   ,& ! local part of the nudging weights
    tabdif           ,& ! time distance between the checking and checked obs.
    zbias            ,& ! bias used for the spatial consistency quality control
    zadjiwv          ,& ! bias corrected obs increment
    qctiq            ,& ! actual threshold for quality control
    wt1    , wt2        ! temporal nudging weight at current obs. station

  LOGICAL                  ::  &
    local            ,& ! obs. station lies on local domain
    lwrqc               ! data which are close the current time, and which are
                        ! rejected by the threshold quality control, are printed
                        ! to a separate file at the current timestep

! Local (automatic) arrays: None
! -------------------------
!
!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine mult_iwv_horicheck
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1: Selection of observations, for which the 
!             spatial consistency check needs to be done
!-------------------------------------------------------------------------------

! the correlation scale for radiosonde humidity at 850 hPa
! is used here in the spatial check for all types of IWV data
  r1iflp  = MAX(  (c1 +(rhtfac(4)-c1) *c05)                                    &
                * (rhinfl(4) + rhvfac(4) *rhqsond(2)) , epsy )
! zcutop  = cutofr(4) * r1iflp
! zcutop2 = zcutop  **2
! zcutmlf = (cutofr(2) / cutofsu(2)) **2
  ntyp    = 1
  IF (mwtypqc == 1) ntyp = 2
  IF (mwtypqc == 2) ntyp = 3
! fmultw  = REAL ( kwtyp(1) - 1, wp )
! IF (ntpy >= 2)  fmultw  = REAL ( kwtyp(nwtyp+1) - 1, wp )


loop_over_all_stations:  DO istaiv = 1 , nivtot

  istaqc = isrtvqc(istaiv)

! check for the criteria for doing the spatial consistency check:
! - the station lies on the local domain
! - the observation is to be quality controlled at the current timestep
! - the obs. has failed to pass the stringent threshold quality control
! ---------------------------------------------------------------------

  ii  =  ioiv (istaqc)
  jj  =  joiv (istaqc)
  local =       (jj <= jend) .AND. (jj >= jstart)                              &
          .AND. (ii <= iend) .AND. (ii >= istart)

  IF ((local) .AND. (iqcliv(istaqc) > 0)) THEN

!   IF (lwonl) PRINT *,'ZGQ1 ',ystidiv(istaqc), kioiiv(istaqc), kiobiv(istaqc) &
!                             ,ioiv(istaqc), joiv(istaqc), istaqc              &
!                             ,kviflml(MAX(kioiiv(istaqc),1),:,3)

!   io_tot = ioiv_tot (istaqc)
!   jo_tot = joiv_tot (istaqc)
!   r1iflp = r1ifps   (istaqc)
    ystid  = ystidiv  (istaqc) (1:ilstidp)
    kobtyp = ktypiv   (istaqc)
    nredowt = 0
    nqexi (1) = MOD( iqcliv(istaqc) , 2 )
    nqexi (2) =      iqcliv(istaqc) / 2

    itica = 2 - MOD( iqcliv(istaqc) , 2 )
    itice = 1 +      iqcliv(istaqc) / 2

!   IF (ystid == '10739')  PRINT '("IQssc ",A ,4I2)' , ystid, itica,itice, nqexi

    loop_over_time_index:  DO itic = itica , itice

      DO ityp = 1 , ntyp
        omyiv  (ityp) = c0
        om2iv  (ityp) = c0
        omytp  (ityp) = c0
        ziwvwi (ityp) = c0
      ENDDO
      zreftim    = zactiv(istaqc,itic)

!-------------------------------------------------------------------------------
!  Section 2: Computation of the 'analysis increment' or 'bias' used for the 
!             spatial consistency check, which is set to a weighted sum of 
!             spreaded observation increments from the neighbouring IWV data.
!-------------------------------------------------------------------------------

      loop_over_other_stations:  DO istaot = 1 , nivtot

        ista   = isrtvqc(istaot)

        io     = ioiv   (ista)
        jo     = joiv   (ista)
        ityp   = 1
        r1iflpk  =  r1iflp
        dtchkiq  =  dtchkps
        IF (ktypiv(ista) == ngps) THEN
          IF (mwtypqc >= 2)  ityp = 3
          r1iflpk  =  r1iflp * rhfgps
          dtchkiq  =  dtchkps * c05
        ELSEIF (ktypiv(ista) == nsattv) THEN
          r1iflpk  =  r1iflp * rhfrtv(4)
          dtchkiq  =  dtchkps * c05 * c05
        ELSE   !   radiosonde
          IF (mwtypqc >= 1)  ityp = 2
        ENDIF
        zcutop  = cutofr(4) * r1iflpk
        zcutop2 = zcutop  **2

! check if obs. 'ista' is within (horizontal) area of influence of obs. 'istaqc'
! (and if obs. 'ista' is not already rejected by the threshold quality control)
        zhcut2 = zcutop2
!       IF (lmlps(ista)) zhcut2 = zhcut2 * zcutmlf
        idix   = ii - io + ie_tot
        jdix   = jj - jo + je_tot

        zdist2 = pxxd2(idix,jj) + pyyd(jdix) *pyyd(jdix)                       &
                                - pyyd(jdix) *pxsalpa(idix,jj)
        IF (      (zdist2 <= zcutop2) .AND. (iqcfiv(ista) > 0)                 &
            .AND. (ystidiv(ista) (1:ilstidp) /= ystid (1:ilstidp))             &
                  ! do not use buddy GPS obs from same station, even if it is
                  ! processed by a different processing centre
                  ! (GPS station ID: digits 1-4: station; 6-9: process. centre)
            .AND. (     (ktypiv(ista) /= ngps)                                 &
                   .OR. (ystidiv(ista) (1:4) /= ystid (1:4)))) THEN
!           .AND. (ista /= istaqc)) THEN

! get the horizontal weight, the spreading correction factor to obs. increment
! (designed to avoid orographic footprints), and the obs. quality factor

          zsprcor =  zsativ(istaqc) / zsativ(ista)
          zdist   =  SQRT( zdist2 )
          zddsp   =  zdist / r1iflpk
          omysc   =  (c1 + zddsp) *EXP( -zddsp )
          omysc2  =  omysc * omysc

! get the temporal weights (different from the temporal nudging weights,
!   without temporal linear interpolation) and quality weights

          omyk1   =  c0
          omyk2   =  c0
          zincr1  =  c0
          zincr2  =  c0
!         IF ((zactiv(ista,1) > rmdich) .AND. (MOD(iqcfiv(ista),2) == 1)) THEN
          IF (MOD(iqcfiv(ista),2) == 1) THEN
            tabdif  =  ABS( zreftim - zactiv(ista,1) )
            omyk1   =  MAX( c0 , c1 - tabdif / dtchkiq )  *  zqcfiv(ista,1)
            zincr1  =  omyk1 * zoiciv(ista,1) * zsprcor
          ENDIF
!         IF ((zactiv(ista,2) > rmdich) .AND. (    iqcfiv(ista)    >= 2)) THEN
          IF (    iqcfiv(ista)    >= 2) THEN
            tabdif  =  ABS( zreftim - zactiv(ista,2) )
            omyk2   =  MAX( c0 , c1 - tabdif / dtchkiq )  *  zqcfiv(ista,2)
            zincr2  =  omyk2 * zoiciv(ista,2) * zsprcor
          ENDIF
          zincrt  =  zincr1 + zincr2
          zinct2  =  omyk1 *zincr1 + omyk2 *zincr2

! update the sum of weights and weighted obs. increments

          omyiv (ityp)  =  omyiv (ityp) + omysc  *(omyk1         + omyk2)
          om2iv (ityp)  =  om2iv (ityp) + omysc2 *(omyk1 *omyk1  + omyk2 *omyk2)
          ziwvwi(ityp)  =  ziwvwi(ityp) + omysc2 *zinct2 + omysc *zincrt *fmultw

        ENDIF

      ENDDO loop_over_other_stations

! compute the bias (i.e. analysis increment) used for the spatial
! consistency quality control check
! (without close neighbouring obs., the bias will be small due to 'MAX(.,1)')

      zwsum  =  c0
      DO ityp = 1 , ntyp
        IF ((om2iv(ityp) > epsy) .AND. (omyiv(ityp) > epsy)) THEN
! averaged increment for observation type 'ityp'
          ziwvwi (ityp)  =   ziwvwi(ityp) / (om2iv(ityp) + fmultw *omyiv(ityp))
! individual weight 'w(ityp)' assigned to the observation type 'ityp'
          omytp  (ityp)  =  (om2iv(ityp) + fmultw *omyiv(ityp))                &
                          / (omyiv(ityp) + fmultw)
! sum of the individual weights 'w(ityp)' assigned to the different obs. types
          zwsum  (1)     =   omytp(ityp)  +  zwsum(1)
!   (sum of squared individual weight used for scaling QC threshold only)
          zwsum  (4)     =   om2iv(ityp)  +  zwsum(4)
        ELSE
          ziwvwi (ityp)  =   c0
          omytp  (ityp)  =   c0
        ENDIF
      ENDDO
      IF (zwsum(1) > epsy) THEN
        DO ityp = 1 , ntyp
! total multiple weight 'W(ityp)' assigned to observation type 'ityp'
          omyiv  (ityp)  =   omytp(ityp) * (omytp(ityp) + fmultw)              &
                                         / (zwsum(1)    + fmultw)
! weighted averaged increment for observation type 'ityp'
          ziwvwi (ityp)  =   ziwvwi(ityp)  * omyiv(ityp)
! sum of total multiple weights 'W(ityp)' assigned to the different obs. types
          zwsum  (2)     =   omyiv (ityp)  +  zwsum(2)
! sum of the weighted averaged increments for the different observation types
          zwsum  (3)     =   ziwvwi(ityp)  +  zwsum(3)
        ENDDO
      ENDIF

! double-averaged total observation increment
      zbias  =  c0
      IF (zwsum(2) > epsy)  zbias  =  zwsum(3) / zwsum(2)
!     zbias  =  zwsum(3) / MAX( zwsum(2) , c1 )

! if the sum of weights of the observations contributing to the total obs
! increment 'zbias' is small, then the new estimate of the true value should
! not be composed of the model first guess plus the full obs increment 'zbias'.
! Instead, the model first guess should be corrected only by a small fraction
! of 'zbias' (because the value for 'zbias' stands on shaky ground).
! Generally, 'zbias' should be scaled by a factor dep on the sum of weights.

      zbias  =  zbias * zwsum(4) / (c1 + zwsum(4))

!-------------------------------------------------------------------------------
!  Section 3: Setting of flags by the spatial consistency quality control check
!-------------------------------------------------------------------------------

! adjusted increment = obs - 'truth' = obs - (model + zbias) = old-incr - zbias
      zadjiwv  =  zoiciv(istaqc,itic)  -  zbias

      lwrqc    =  (lqcall) .AND. (ABS(c3600*(zreftim - acthr)+epsy) <= cqcwbox)

! QC threshold: 0.15 * IWV sat   (independently of data density 'zwsum(2)')
      qctiq    = qcsiq * zsativ(istaqc)  +  qcciq
      qctiq    = (c1 + qctf(4) *ABS(zreftim - acthr))  * qctiq
      qctiq    = (c1 - (c1 - qcftsiv) *MIN( zwsum(4) /5.0_wp , c1 )) * qctiq
      qctiq    =  qctiq  +  qcfbias *ABS( zbias )

! radiosonde humidity: only set additional rejection flags
! -------------------

      IF ((ABS( zadjiwv ) > qctiq) .AND. (kobtyp /= ngps)                      &
                                   .AND. (kobtyp /= nsattv)) THEN

        nmlob  =  imladm(kiobiv(istaqc),itic)
        nlev   =  momlhd(nmlob,nhnlev)
        ilbot  =  0
        iltop  =  0

! set flag for 'not containing active humidity data'
        nqexi (itic) = 0

        DO ilev = 1 , nlev
! set the quality control flag for humidity
          IF (omlbdy(nmlob,ilev,nbtrh ) > rmdich)                              &
            momlbd (nmlob,ilev,nbtqcf) = IBSET( momlbd(nmlob,ilev,nbtqcf), nvrq)
! determine range of humidity levels for diagnostic output 'oyqc'
          IF (omlbdy(nmlob,ilev,nbtrh ) > rmdich) THEN
            IF (ilbot == 0) ilbot = ilev
            iltop = ilev
          ENDIF
        ENDDO

! set humidity-related entries in oiml to missing
!           (if .not.lreset in mult_org_localinfo) !
        icml   =  kioiiv(istaqc)
        IF (icml > 0) THEN
          nlev   =  kkoiml(icml,itic) + mszlev(icml,itic) - 1
          DO klvi = kkoiml(icml,itic) , nlev
            oiml (noiqd ,klvi) = rmdi
            oiml (noirh ,klvi) = rmdi
            oiml (noiqqc,klvi) = c1
            oiml (noiqlr,klvi) = c0
          ENDDO
          kviflml (icml,2*itic-1,3) = 0
          kviflml (icml,2*itic  ,3) = ke + 1
        ENDIF

        IF (ntotqc < maxqcp) THEN
          ntotqc = ntotqc + 1
          myqc (ntotqc, 2) = 17
          yyqc (ntotqc   ) = ystid
          myqc (ntotqc, 1) = momlhd(nmlob,nhcode)
          IF (.NOT. lwrqc) myqc(ntotqc,1) = 0
          oyqc (ntotqc, 1) = omlhed(nmlob,nhtime)
          oyqc (ntotqc, 2) = omlbdy(nmlob,iltop,nbtp)  / 100.0_wp
          oyqc (ntotqc, 3) = omlhed(nmlob,nhjlat)
          oyqc (ntotqc, 4) = omlhed(nmlob,nhilon)
          oyqc (ntotqc, 5) = qctiq
          oyqc (ntotqc, 6) = zoiciv(istaqc,itic) + zmodiv(istaqc)
          oyqc (ntotqc, 7) = zmodiv(istaqc)
          oyqc (ntotqc, 8) = zbias
          oyqc (ntotqc, 9) = zwsum(4)
        ENDIF

! set flag for need for adjustment of temporal nudging weights
        nredowt = nredowt + itic

! GPS observations: switch flag either to rejection or acceptance where required
! ----------------

      ELSEIF (kobtyp == ngps) THEN
        ngpob  =  igpadm(kiobiv(istaqc),itic)
        mflgqc = NINT( c05 + SIGN(c05, ABS( zadjiwv ) - qctiq) )

        IF (      (ibit1( mogpbd(ngpob,nbgerr),nvriwv ) == 1)                  &
            .AND. (ibit1( mogpbd(ngpob,nbgqcf),nvriwv ) /= mflgqc)) THEN
          IF (ntotqc < maxqcp) THEN
            ntotqc = ntotqc + 1
            myqc (ntotqc, 2) = 17
            yyqc (ntotqc   ) = ystid
            myqc (ntotqc, 1) = mogphd(ngpob,nhcode)
!           myqc (ntotqc, 1) = mogphd(ngpob,nhobtp)
            IF (.NOT. lwrqc) myqc(ntotqc,1) = 0
            oyqc (ntotqc, 1) = ogphed(ngpob,nhtime)
            oyqc (ntotqc, 2) = ogpbdy(ngpob,nbgp)   / 100.0_wp
            oyqc (ntotqc, 3) = ogphed(ngpob,nhjlat)
            oyqc (ntotqc, 4) = ogphed(ngpob,nhilon)
            oyqc (ntotqc, 5) = qctiq
            oyqc (ntotqc, 6) = zoiciv(istaqc,itic) + zmodiv(istaqc)
            oyqc (ntotqc, 7) = zmodiv(istaqc)
            oyqc (ntotqc, 8) = zbias
            oyqc (ntotqc, 9) = zwsum(4)
          ENDIF

! change sign of the observation error according to the check
!         ogpbdy (ngpob,nbgtze) = - ogpbdy(ngpob,nbgtze)

! adjust the quality control flag: replace 0 by 1 , 1 by 0 , 2 by 3 , or 3 by 2
!         IF (mogpbd(ngpob,nbgqcf) >= 0)                                       &
!           mogpbd (ngpob,nbgqcf) =   mogpbd(ngpob,nbgqcf)                     &
!                                   + 1 - 2* MOD( mogpbd(ngpob,nbgqcf) , 2 )

! set flag for need for adjustment of temporal nudging weights
!         nredowt = nredowt + itic
        ENDIF
        CALL MVBITS ( mflgqc, 0, 1, mogpbd(ngpob,nbgqcf), nvriwv )
        CALL MVBITS ( mflgqc, 0, 1, mogpbd(ngpob,nbgqcf), nvrzpd )

! if IWV rejected by spatial consistency check
!    (irrespective of QC in 'mult_obs_gps')
! set humidity-related entries in oiml to missing 
!           (if .not.lreset in mult_org_localinfo) !
        IF ((ABS( zadjiwv) > qctiq) .AND. (kioiiv(istaqc) > 0)) THEN
          icml   =  kioiiv(istaqc)
          nlev   =  kkoiml(icml,itic) + mszlev(icml,itic) - 1
          DO klvi = kkoiml(icml,itic) , nlev
            oiml (noiqd ,klvi) = rmdi
            oiml (noirh ,klvi) = rmdi
            oiml (noiqqc,klvi) = c1
            oiml (noiqlr,klvi) = c0
          ENDDO
          kviflml (icml,2*itic-1,3) = 0
          kviflml (icml,2*itic  ,3) = ke + 1
!       ELSEIF (ABS( zadjiwv) <= qctiq) THEN
! 'oiml' was already set in 'mult_obs_gps' for those GPS obs which 
! could possibly be re-set to active in this spatial consistency check
        ENDIF

! Satellite retrievals
! --------------------

      ELSEIF (kobtyp == nsattv) THEN
        ntvob  =  itvadm(kiobiv(istaqc),itic)
        mflgqc = NINT( c05 + SIGN(c05, ABS( zadjiwv ) - qctiq) )

        IF (ibit1( motvhd(ntvob,nhvqcf),nvrq ) /= mflgqc) THEN
!         IF ((ntotqc < maxqcp) .AND. (motvhd(ntvob,nhqexi) >= 1)) THEN
          IF  (ntotqc < maxqcp) THEN
            ntotqc = ntotqc + 1
            myqc (ntotqc, 2) = 17 
            yyqc (ntotqc   ) = ystid
            myqc (ntotqc, 1) = motvhd(ntvob,nhcode)
!           myqc (ntotqc, 1) = motvhd(ntvob,nhobtp)
            IF (.NOT. lwrqc) myqc(ntotqc,1) = 0
            oyqc (ntotqc, 1) = otvhed(ntvob,nhtime)
            oyqc (ntotqc, 2) = otvbdy(ntvob,1,nbvp)   / 100.0_wp
            oyqc (ntotqc, 3) = otvhed(ntvob,nhjlat)
            oyqc (ntotqc, 4) = otvhed(ntvob,nhilon)
            oyqc (ntotqc, 5) = qctiq
            oyqc (ntotqc, 6) = zoiciv(istaqc,itic) + zmodiv(istaqc)
            oyqc (ntotqc, 7) = zmodiv(istaqc)
            oyqc (ntotqc, 8) = zbias
            oyqc (ntotqc, 9) = zwsum(4)
          ENDIF

! change sign of the observation error according to the check, and
! adjust the quality control flag (for sat retrievals, this is the same)
!         motvhd (ntvob,nhvqcf) =       motvhd(ntvob,nhvqcf)                   &
!                                 - 8 *(motvhd(ntvob,nhvqcf) /4) + 4

! set flag for need for adjustment of temporal nudging weights
          nredowt = nredowt + itic
        ENDIF
! adjust QC flag (this is the only place where element 'nhvqcf' is modified)
        CALL MVBITS ( mflgqc, 0, 1, motvhd(ntvob,nhvqcf), nvrq )

! if IWV rejected by spatial consistency check
! set humidity-related entries in oiml to missing
!           (if .not.lreset in mult_org_localinfo) !
        IF ((ABS( zadjiwv) > qctiq) .AND. (kioiiv(istaqc) > 0)) THEN
          icml   =  kioiiv(istaqc)
          nlev   =  kkoiml(icml,itic) + mszlev(icml,itic) - 1
          DO klvi = kkoiml(icml,itic) , nlev
            oiml (noiqd ,klvi) = rmdi
            oiml (noirh ,klvi) = rmdi
            oiml (noiqqc,klvi) = c1
            oiml (noiqlr,klvi) = c0
          ENDDO
          kviflml (icml,2*itic-1,3) = 0
          kviflml (icml,2*itic  ,3) = ke + 1
        ENDIF
      ENDIF

      IF (lwonl) WRITE( nupr,'("IWV ssc:",A,I2," inc;bias;thr;W",3F6.2,F7.3    &
                             &,I2,3I4,3I2)')  ystid, kobtyp, zadjiwv, zbias    &
                             , qctiq, zwsum(4), nredowt, istaiv, istaqc        &
                             , kioiiv(istaqc), itic, nqexi(itic), iqcliv(istaqc)

    ENDDO loop_over_time_index

!-------------------------------------------------------------------------------
!  Section 4: Adjustment of the local (temporal and quality) nudging weights
!-------------------------------------------------------------------------------

    IF ((nredowt > 0) .AND. (kioiiv(istaqc) > 0)) THEN
      icml  =  kioiiv(istaqc)
! radiosonde: only case 1: wt_old > 0 --> wt_new = 0
!     IF (kobtyp /= ngps) THEN
!       itica = 2 - MOD( nredowt , 2 )
!       itice = 1 +      nredowt / 2
!       DO itic = itica , itice
!       DO itic = 1 , 2
!         nmlob  = imladm(kiobiv(istaqc),itic)
!         nmlob2 = imladm(kiobiv(istaqc),3-itic)
!         zwtml (icml,itic,3) = c0
!         IF (zwtml(icml,3-itic,3) > epsy)                                     &
!         zwtml (icml,3-itic,3) = wtml(kiobiv(istaqc),itic+2)
!       ENDDO
      IF ((kobtyp /= ngps) .AND. (kobtyp /= nsattv)) THEN
        DO itic = 1 , 2
          nmlob  =  imladm(kiobiv(istaqc),itic)
          nmlob2 =  imladm(kiobiv(istaqc),3-itic)
          !  (not separately checking for Mode-S ok as long as it has no RH obs)
          IF (     (      (kobtyp /= nairep) .AND. (gnudg  (4) < epsy))        &
              .OR. (      (kobtyp == nairep) .AND. (gnudgar(4) < epsy))        &
!             .OR. (      (kobtyp == nairep)                                   &
!                   .AND. (kcdtyp /= nmodes) .AND. (gnudgar(4) < epsy))        &
!             .OR. (      (kobtyp == nairep)                                   &
!                   .AND. (kcdtyp == nmodes) .AND. (gnudgms(4) < epsy))        &
              .OR. (momlhd(MAX(nmlob,1),nhpass) /= 0)                          &
              .OR. (nqexi(itic) == 0)) THEN
            wt2 = c0
          ELSEIF ((nmlob2 == 0) .OR. (nqexi(3-itic) == 0)) THEN
            wt2 = wtml(kiobiv(istaqc),itic+2)
          ELSE
            wt2 = wtml(kiobiv(istaqc),itic)
          ENDIF
          IF (itic == 1) wt1 = wt2
        ENDDO
! GPS: case 1: wt_old > 0 --> wt_new = 0 ;  case 2: wt_old = 0 --> wt_new > 0
      ELSEIF (kobtyp == ngps) THEN
        DO itic = 1 , 2
          ngpob  =  igpadm(kiobiv(istaqc),itic)
          ngpob2 =  igpadm(kiobiv(istaqc),3-itic)
!         IF ((gnudggp <= epsy) .OR. (ogpbdy(MAX(ngpob,1),nbgtze) <= c0)       &
!                               .OR. (mogphd(MAX(ngpob,1),nhpass) /= 0)) THEN
          IF (     (gnudggp <= epsy) .OR. (mogphd(MAX(ngpob,1),nhpass) /= 0)   &
              .OR. (.NOT. BTEST( mogpbd(MAX(ngpob,i1),nbgerr),nvriwv ))        &
              .OR. (      BTEST( mogpbd(MAX(ngpob,i1),nbgqcf),nvriwv ))) THEN
            wt2 = c0
!         ELSEIF ((ngpob2 == 0) .OR. (ogpbdy(MAX(ngpob2,i1),nbgtze) <= c0)) THEN
          ELSEIF (    (ngpob2 == 0)                                            &
                 .OR. (.NOT. BTEST( mogpbd(MAX(ngpob2,i1),nbgerr),nvriwv))     &
                 .OR. (      BTEST( mogpbd(MAX(ngpob2,i1),nbgqcf),nvriwv))) THEN
            wt2 = wtgp(kiobiv(istaqc),itic+2)
          ELSE
            wt2 = wtgp(kiobiv(istaqc),itic)
          ENDIF
          IF (itic == 1) wt1 = wt2
        ENDDO
! Satrtv: case 1: wt_old > 0 --> wt_new = 0 ;  case 2: wt_old = 0 --> wt_new > 0
      ELSEIF (kobtyp == nsattv) THEN 
        DO itic = 1 , 2
          ntvob  =  itvadm(kiobiv(istaqc),itic) 
          ntvob2 =  itvadm(kiobiv(istaqc),3-itic)
          IF (    (gnudgtv(4) <= epsy) .OR. (motvhd(MAX(ntvob,1),nhqexi) <= 0) &
                                       .OR. (motvhd(MAX(ntvob,1),nhpass) /= 0) &
             .OR. (BTEST( motvhd(MAX(ntvob,1),nhvqcf),nvrq ))) THEN
            wt2 = c0
          ELSEIF (     (ntvob2 == 0) .OR. (motvhd(MAX(ntvob2,i1),nhqexi) <= 0) &
                  .OR. (BTEST( motvhd(MAX(ntvob2,1),nhvqcf),nvrq ))) THEN
            wt2 = wttv(kiobiv(istaqc),itic+2)
          ELSE
            wt2 = wttv(kiobiv(istaqc),itic)
          ENDIF
          IF (itic == 1) wt1 = wt2
        ENDDO
      ENDIF
      IF (lwonl) WRITE( nupr,'("TQV ssc:",A ,2I4,":wt:old1,2/new1,2",4F6.2)' ) &
                                ystid, ioiv_tot(istaqc), joiv_tot(istaqc)      &
                              ,(zwtml(icml,jo,3),jo=1,2), wt1, wt2
      zwtml (icml,1,3) = wt1
      zwtml (icml,2,3) = wt2
    ENDIF

!   icml = kioiiv(istaqc)
!   PRINT *,'ZQssa ', ystidiv(istaqc), icml, zwtml(MAX(1,icml),1,3)            &
!                                          , zwtml(MAX(2,icml),1,3)            &
!                                        , kviflml(MAX(1,icml),:,3)

  ENDIF

ENDDO loop_over_all_stations


!-------------------------------------------------------------------------------
! End of module procedure mult_iwv_horicheck
!-------------------------------------------------------------------------------

END SUBROUTINE mult_iwv_horicheck


!-------------------------------------------------------------------------------
!+ Module procedure to compute IWV from vertical humidity profiles
!-------------------------------------------------------------------------------

SUBROUTINE mult_vprof_2_iwv ( icml, lretv, lqc, lqconly, lprqcm                &
                            , ke, col_t, col_rh, col_qv, col_qc, col_qrs       &
                            ,     col_p, col_z, col_hhl, viobtom , nexceiv )

!-------------------------------------------------------------------------------
!
! Description:
!   This routine computes observation increments of integrated water vapour
!   (IWV) from vertical humidity profiles in order to prepare a spatial
!   consistency check of IWV.
!
! Method:
!   1. Compute on model levels a vertical profile of relative humidity
!      observation increments (including increments given by the vertical
!      weight function below and above the profile) 
!   2. From this and a temperature increment profile, compute a IWV pseudo
!      observation valid at model orography.
!   3. Also compute the model IWV and the IWV value related to the saturated
!      temperature model profile.
!
! Written by        :  Christoph Schraff, DWD
! Current Code Owner:  Christoph Schraff
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

  INTEGER (KIND=iintegers), INTENT (IN)         ::       &
    icml             ,& ! index of obs. sta. in data record to be used further
    ke                  ! number of vertical (main) levels in model column

  LOGICAL                 , INTENT (IN)         ::       &
    lretv            ,& ! multi-level report is satellite retrieval
    lqc              ,& ! threshold quality control to be applied
    lqconly          ,& ! threshold quality control without other processing
    lprqcm (2)          ! printout for control at present station / node / time

  REAL    (KIND=wp   ),     INTENT (IN)         ::       &
                        ! model column (at the observation location) of:
    col_t     (ke)   ,& ! temperature                                    [  K  ]
    col_rh    (ke)   ,& ! model relative humidity                        [     ]
    col_qv    (ke)   ,& ! specific water vapor content                   [kg/kg]
    col_qc    (ke)   ,& ! specific cloud water content                   [kg/kg]
    col_qrs   (ke)   ,& ! spec. cont. of hydrometeors excl. cloud water  [kg/kg]
    col_p     (ke)   ,& ! pressure (full value)        on main levels    [ Pa  ]
    col_z     (ke)   ,& ! geometrical height           of main levels    [  m  ]
    col_hhl   (ke+1) ,& ! geometrical height           of half levels    [  m  ]
!   col_dp0   (ke)   ,& ! reference pressure thickness of model layers   [ Pa  ]
    viobtom   (ke,4)    ! obs. (incr.) interpolated to model levels

  INTEGER (KIND=iintegers), INTENT (INOUT)      ::       &
    nexceiv             ! number of IWV obs incr. reports exceeding array size

! Local parameters:
! ----------------
! REAL (KIND = wp)         , PARAMETER  :: &
!   c300hpa = 30000.0_wp      ! 300 hPa  in [Pa]

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    modespr          ,& ! mode of spreading (horizontal, along model levels
                        !                   ,along isentropes)
    nmlob  , ntvob   ,& ! index of report in the ODR (obs data record)
    mlev             ,& ! vertical loop index over model levels
    mdabot , mdatop     ! index of lowest / top model level with interpol. obs.

  REAL    (KIND=wp   )     ::  &
    ztimob           ,& ! time of report in hrs rel. to initial model time
    rdsprs           ,& ! vertical correlation scale
    zcutbot, zcuttop ,& ! height range influenced by temperature increments
    ztqobs           ,& ! IWV derived from observed T- and q- profile
    ztvio            ,& ! observed temperature interpolated to model level
    zqvsat           ,& ! model specific humidity at saturation
    zpvsat           ,& ! model saturation vapour pressure
    zqvvio           ,& ! observed specific humidity at model levels
    zrho             ,& ! density of moist air
    zrhodz           ,& ! density of moisty air * thickness of model layer
    ztqmod           ,& ! IWV derived from model T- and q- profile
    ztqsat           ,& ! model IWV when assuming saturation
    ztqsto , ztqsrg  ,& ! model saturation IWV up to 300 hPa / over q-obs range
    omykf            ,& ! local part of the vertical nudging weight
    zvcotq           ,& ! at top  of profile:  /  correction factor to vertical
    zvcobq              ! at base of profile:  \  correlation scale for humidity

  LOGICAL                  ::  &
    lprmlc           ,& ! printout for control of multi-level check
    lqciwv              ! use obs for IWV spatial consistency checking

! Local (automatic) arrays:
! -------------------------

  REAL    (KIND=wp   )     ::  &
    ztoi    (ke)     ,& ! temperature obs. increment profile at model layers
    zrhoi   (ke)        ! rel. humidity obs. increment profile at model layers

!------------ End of header ----------------------------------------------------


 
!-------------------------------------------------------------------------------
! Begin Subroutine mult_vprof_2_iwv
!-------------------------------------------------------------------------------
  
!-------------------------------------------------------------------------------
!  Section 1: Preliminaries
!-------------------------------------------------------------------------------

! nmloi = knoiml(icml,itim)

! initialisations and preliminaries
! ---------------------------------

                           lprmlc = .FALSE.
  IF (kobtyp /= ot_satem)  lprmlc =  (lprqcm(1)) .AND. (lprqcm(2))

  IF (.NOT. lretv) THEN
    nmlob  = imladm(ista,itim)
    ztimob =  omlhed(nmlob,nhtime)
    zvcobq =  omlhed(nmlob,nhvcbq)
    zvcotq =  omlhed(nmlob,nhvctq)
  ELSEIF (lretv) THEN
    ntvob  = itvadm(ista,itim)
    nmlob  = ntvob
    ztimob =  otvhed(ntvob,nhtime)
    zvcobq =  c1
    zvcotq =  c1
  ENDIF

  modespr = msprpar

! auxiliary quantities

  DO mlev = ke , 1 , -1
    ztoi  (mlev) = c0
    zrhoi (mlev) = c0
  ENDDO

!-------------------------------------------------------------------------------
!  Section 5:  Computation of observation increment of IWV (for TEMP only)
!-------------------------------------------------------------------------------
!CS: Note: this section is done only once at the very end of the single or
!          (due to 'lredo') pair of calls of mult_vertic_intpol.

! (not done, if there are no humidity data)
  lqciwv =  (lqc) .AND. (mtoplv(icml,itim,3)+mbotlv(icml,itim,3) < ke)         &
                  .AND. (liwvssc)                     .AND. (.NOT. lretv)
!                 .AND. (liwvssc) .AND. (.NOT. lredo) .AND. (.NOT. lretv)

! IF (ystid == '10739') PRINT '("IQoi ",A,3I5,2L3)',ystid, mtoplv(icml,itim,3) &
!                              ,mbotlv(icml,itim,3), ntstep, lqc , lqciwv

  IF (lqciwv) THEN

    mdatop  =  1 + mtoplv(icml,itim,3)
    mdabot  = ke - mbotlv(icml,itim,3)
    IF (lprmlc)  WRITE( nupr,'("zrhtb ",A,2I4)' ) ystid, mtoplv(icml,itim,3)   &
                                                       , mbotlv(icml,itim,3)

! a quality weight for the IWV increment dep. on vertical range of humidity obs
!   omykf   =  (col_p(mdabot) - MAX( c300hpa + col_dp0(kml300)                 &
!                                  , col_p(mdatop) )) &
!            / (col_p(ke)     -     (c300hpa + col_dp0(kml300)))

! construction of a relative humidity observation increment profile 'zrhoi'
! (in an analogous way as 'ztoi' in section 4)
! (here, temperature increments (e.g. due to vertical spreading)
!  below and above the observed humidity profile are neglected)
! -------------------------------------------------------------------------

    DO mlev = mdatop , mdabot
      ztoi  (mlev) = viobtom(mlev,3) - col_t (mlev)
      zrhoi (mlev) = viobtom(mlev,4) - col_rh(mlev)
    ENDDO

    !   below the observed profile
    IF (mdabot < ke) THEN
      rdsprs  = SQRT( vcorls(4) )
      IF (modespr <= 1) rdsprs = rdsprs * zrtgkml(icml,mdabot)
      IF (modespr == 2) rdsprs = rdsprs / thdzt
      zcutbot = col_z(mdabot) -  rdsprs *SQRT( vcutof(4) )                     &
                              *((c1-fsvcut) + fsvcut *zvcobq)
      rdsprs  = rdsprs * zvcobq
      DO mlev = mdabot + 1 , ke
        IF (col_z(mlev) >= zcutbot)                                            &
          zrhoi (mlev) = zrhoi(mdabot) *EXP( -(  (col_z(mdabot) - col_z(mlev)) &
                                               / rdsprs) **2 )
      ENDDO
    ENDIF

    !   above the observed profile
!   IF (mdatop > 1) THEN         ! is always true for humidity
      rdsprs  = SQRT( vcorls(4) )
      IF (modespr <= 1) rdsprs = rdsprs * zrtgkml(icml,mdatop)
      IF (modespr == 2) rdsprs = rdsprs / thdzt
      zcuttop = col_z(mdatop) +  rdsprs *SQRT( vcutof(4) )                     &
                              *((c1-fsvcut) + fsvcut *zvcotq)
      rdsprs  = rdsprs * zvcotq
      DO mlev = mdatop - 1 , 1 , - 1
        IF (col_z(mlev) <= zcuttop)                                            &
          zrhoi (mlev) = zrhoi(mdatop) *EXP( -(  (col_z(mlev) - col_z(mdatop)) &
                                               / rdsprs) **2 )
      ENDDO
!   ENDIF

! determine observed IWV (with respect to model orography)
! --------------------------------------------------------

    ztqobs = c0
    DO mlev = ke , 1 , -1
! saturation vapour pressure for observed temperature
      ztvio   =  col_t(mlev) + ztoi(mlev)
      zpvsat  =  fpvsw ( ztvio , b1, b2w, b3, b4w )
! observed specific humidity at model levels
      zqvvio  =  fpv2q( (col_rh(mlev) + zrhoi(mlev)) *zpvsat , col_p(mlev), rdv)
! density of moist air
      zrho    =  col_p(mlev) /(r_d * ztvio * (  c1 + rvd_m_o*zqvvio            &
                                              - col_qc(mlev) - col_qrs(mlev)) )
! IWV derived from observed temperature and humidity profile
      zrhodz  =  zrho *(col_hhl(mlev) - col_hhl(mlev+1))
      ztqobs  =  ztqobs  +  zrhodz * zqvvio
    ENDDO

! determine model IWV and IWV of saturated model temperature profile
! ------------------------------------------------------------------

    ztqmod = c0
    ztqsat = c0
    ztqsto = c0
    ztqsrg = c0
    DO mlev = ke , 1 , -1
! model saturation vapour pressure
      zpvsat  =  fpvsw ( col_t(mlev) , b1, b2w, b3, b4w )
! model saturation specific humidity
      zqvsat  =  fpv2q ( zpvsat , col_p(mlev) , rdv )
! density of moist air
      zrho    =  col_p(mlev) /( r_d * col_t(mlev)                              &
                                    * (  c1 + rvd_m_o*col_qv(mlev)             &
                                       - col_qc(mlev) - col_qrs(mlev)) )
! IWV derived from model temperature and humidity profile
      zrhodz  =  zrho *(col_hhl(mlev) - col_hhl(mlev+1))
      ztqmod  =  ztqmod  +  zrhodz * col_qv(mlev)
! model IWV when assuming saturation
      ztqsat  =  ztqsat  +  zrhodz * zqvsat
      IF (mlev >= kml300) THEN
        ztqsto  =  ztqsto  +  zrhodz * zqvsat
        IF ((mlev >= col_p(mdabot)) .AND. (mlev <= col_p(mdatop)))             &
          ztqsrg  =  ztqsrg  +  zrhodz * zqvsat
      ENDIF
    ENDDO

! a quality weight for the IWV increment dep. on vertical range of humidity obs
    omykf   =  ztqsrg / ztqsto

! store IWV obs increment with respect to model orography
! -------------------------------------------------------

!   IF (     (kioiiv(MAX(nivtot,1)) /= icml)                                   &
!       .OR. (ktypiv(MAX(nivtot,1)) /= kobtyp))  nivtot = nivtot + 1
    mlev = kiobiv(MAX(nivtot,1))
    IF (nivtot <= 0) THEN
      nivtot = 1
      iqcliv (nivtot)  =  0
      iqcfiv (nivtot)  =  0
    ELSEIF (     (kiobiv(MAX(nivtot,1)) /= ista)                               &
            .OR. (ktypiv(MAX(nivtot,1)) /= kobtyp)) THEN
      IF (nivtot < maxivq) THEN
        nivtot = nivtot + 1
        iqcliv (nivtot)  =  0
        iqcfiv (nivtot)  =  0
      ELSE
        nexceiv = nexceiv + 1
      ENDIF
    ENDIF

    ystidiv  (nivtot)     =  ystid
    ioiv     (nivtot)     =  io
    joiv     (nivtot)     =  jo
    ioiv_tot (nivtot)     =  io_tot
    joiv_tot (nivtot)     =  jo_tot
    ktypiv   (nivtot)     =  kobtyp
    kioiiv   (nivtot)     =  icml
    kiobiv   (nivtot)     =  ista
    zactiv (nivtot,itim)  =  ztimob
    zoiciv (nivtot,itim)  =  ztqobs - ztqmod
    zqcfiv (nivtot,itim)  =  omykf
    zsativ   (nivtot)     =  ztqsat
    zmodiv   (nivtot)     =  ztqmod
    IF (itim == 1) THEN
      iqcliv (nivtot)     =  2 * (iqcliv(nivtot) / 2)  + 1
      iqcfiv (nivtot)     =  2 * (iqcfiv(nivtot) / 2)  + 1
    ELSE
      iqcliv (nivtot)     =  MOD( iqcliv(nivtot) , 2 ) + 2
      iqcfiv (nivtot)     =  MOD( iqcfiv(nivtot) , 2 ) + 2
    ENDIF

!   PRINT *,'ztqi  ', ystid, ztqobs, ztqmod, ztqsat
!   IF (MOD(ntstep , 90) <= 1)                                                 &
    IF (lprmlc)  WRITE( nupr,'("IWV oi:",A ,7I4,F6.2,3F6.2,F6.3)' )            &
                        ystid, kobtyp, ista, icml, nivtot, mlev                &
                             , iqcliv(nivtot), iqcfiv(nivtot), ztimob          &
                             , zoiciv(nivtot,itim), ztqmod, ztqsat, omykf

  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure mult_vprof_2_iwv
!-------------------------------------------------------------------------------

END SUBROUTINE mult_vprof_2_iwv

!===============================================================================

ELEMENTAL REAL (KIND=wp) FUNCTION fpvsw  ( zt, b1, b2w, b3, b4w )
  !----------------------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  zt, b1, b2w, b3, b4w
  !----------------------------------------------------------------------------
  ! Magnus formula for water:  input  'zt'   : temperature
  !                            output 'fpvsw': saturation water vapour pressure
  ! Magnus formula for ice  :  if constants 'b2i', 'b4i' for ice are used for
  !                            'b2w', 'b4w' in the call
  !- --------------------------------------------------------------------------
  !
  fpvsw  =  b1 * EXP( b2w*(zt-b3)/(zt-b4w) )
  !
END FUNCTION fpvsw

!-------------------------------------------------------------------------------

ELEMENTAL REAL (KIND=wp) FUNCTION ftd  ( zpv, b1, b2w, b3, b4w )
  !---------------------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  zpv, b1, b2w, b3, b4w
  REAL    (KIND=wp)                        ::  zlogpv
  !---------------------------------------------------------------------------
  ! inverse of Magnus formula:  input  'zpv' :  water vapour pressure
  !                             output 'ftd' :  dewpoint temperature
  !---------------------------------------------------------------------------
  !
  zlogpv  =  LOG( zpv / b1 )
  ftd     =  (b3*b2w - zlogpv*b4w) / (b2w - zlogpv)
  !
END FUNCTION ftd

!-------------------------------------------------------------------------------

ELEMENTAL REAL (KIND=wp) FUNCTION fpv2q  ( zpv, zp, rdv )
  !--------------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  zpv, zp, rdv
  !---------------------------------------------------------------------------
  ! specific humidity from water vapour pressure and air pressure
  !---------------------------------------------------------------------------
  !
  fpv2q  =  rdv * zpv / (zp - (c1-rdv)*zpv)
  !
END FUNCTION fpv2q

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

END MODULE src_mult_local
