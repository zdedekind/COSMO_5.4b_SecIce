!+ Data module for the variables only used to compute the local information
!-------------------------------------------------------------------------------

MODULE data_nudge_all

!-------------------------------------------------------------------------------
!
! Description:
!   This module contains those variables (scalars and arrays) which must be
!   generally available in the data assimilation (nudging) part of the code
!   (i.e. which do not meet the restrictions imposed on the other data modules
!    used for the data assimilation). These are
!    - all NAMELIST variables controlling the data assimilation
!    - some parameters (constants) and associated variables
!    - some I/O device numbers and file names for nudging
!    - surface analysis limits and allocatable arrays
!    - some other variables, e.g. the switch for the observation processing
!
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.10       1998/09/29 Christoph Schraff
!  Initial release
! 1.11       1998/10/13 Christoph Schraff
!  Removal of BOZ literal constant.
! 1.13       1998/10/22 Christoph Schraff
!  Introduction of switch for observation processing (section 6)
! 1.19       1998/12/11 Christoph Schraff
!  Introduction of switches for the nudging and for the verification at the
!  current timestep, and of variables defining the verification period.
!  Namelist description improved.
! 1.31       1999/07/01 Christoph Schraff
!  Removal of quanitites related to MPI commnicator 'icomm_world'.
!  Introduction of on/off switch for verification also of passive reports.
! 1.36       2000/02/24 Christoph Schraff
!  Namelist variable 'mruntyp' for writing obs. increments to the VOF.
! 1.38       2000/04/06 Christoph Schraff
!  'qcvf(2)' used now for checking height and thickness.
! 2.4        2001/01/29 Christoph Schraff
!  Namelist variables 'khumbal' and 'qcfpst' added. 'zconbas' + 'zcontop' added.
!  Default values and description of namelist variables updated. 'qctf' and 
!  'qctfsu' removed from the list of namelist parameters and shifted to module
!  'data_nudge_local'.
! 2.5        2001/06/01 Christoph Schraff
!  Introduction of selectable pathname for AOF, 'yaofpath', and of 'aiwthr'.
! 2.13       2002/01/18 Christoph Schraff
!  Introduction of Wind Profiler / RASS reports. New code types for obstype 5.
!  For IBM compiler option 'hot': Data statements replaced by direct assignment.
! 2.19       2002/10/24 Christoph Schraff + Michael Buchhold
!  - New namelist variables: 'qgeotop', 'vpblsu', 'loiqv2m', 'lqfqv2m'.
!  - Introduction of Radar VAD wind profiles. New code type for obstype 5.
!  - Introduction of 10m wind speed analysis
! 3.3        2003/04/22 Christoph Schraff
!  New namelist variables 'gnudggp', 'maxgpo', 'lgps', 'lcd096', 'lgpsbias'
!  and new file unit for assimilation of ground-based GPS-derived integrated
!  water vapour data. Code type numbers moved from 'data_obs_process'.
! 3.6        2003/12/11 Christoph Schraff
!  New namelist variable 'rhfgps' for horizontal correlation scale for GPS data.
! 3.7        2004/02/18 Ulrich Schaettler
!  Eliminated initialization of the unit numbers in Section 5.3
! 3.12       2004/09/15 Christoph Schraff
!  Extension to (prepare to) include assimilation of satellite retrievals:
!  Introduce variables that may become namelist parameters in a later version.
! 3.13       2004/12/03 Klaus Stephan / Stefan Klink
!  Introduced new variables for latent heat nudging
! 3.14       2005/01/25 Klaus Stephan
!  Introduced new Namelist variable rlhn_scale_dp: scaling factor for 
!  diagnostic reference precipitation
! 3.16       2005/07/22 Ulrich Schaettler
!  Bug correction: rlhn_scale_dp was erroneously defined as Integer
! 3.18       2006/03/03 Christoph Schraff
!  New namelist variables 'kmultw' for multiple observations, 'lseparw' for
!  multiple observing systems, and 'qcciq' and 'qcsiq' for QC thresholds of IWV,
!  and 'rhfpsdd' for scaling horizontal correlation scale for surface prssure
!  depending on data density.
!  Preparation for use of real-data 1DVar satellite retrievals (MSG, NOAA15-18).
!  LHN variables moved to data_lheat_nudge to avoid to many dependencies
! V4_5         2008/09/10 Christoph Schraff
!  To allow for reading from NetCDF observation input files instead of an AOF
!  file: new namelist variables 'itype_obfile' and 'ycdfdir'.
! V4_9         2009/07/16 Kathleen Helmert
!  Introduction of a threshold for small VAD-winds (epsy_vad)
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_15        2010/11/19 Klaus Stephan
!  New namelist variables:
!  - 'mqcorr92': switch for bias correction of Vaisala RS92 radiosonde humidity
! V4_22        2012/01/31 Christoph Schraff
!  New namelist variables:
!  - 'mveripr' : switches (code) for writing VOF and/or feedobs files
!  - 'mpsgcor' : switch for applying geostrophic pressure increments
!  - 'qgeops'  : reduction factor to the geostrophic pressure increments
!  - 'lscatt'  : for active use of observation type 'scatterometer'
!  - 'maxmlv'  : max. number of observation levels in multi-level report
!  - 'mxfrep'  : max. number of reports in NetCDF feedobs file
!  - 'mxfobs'  : max. number of observations in NetCDF feedobs file
!  - for use of new observation code types: 'lcd037' (TEMP MOBILE),
!    'lcd038' (PILOT MOBILE), 'lcd140' (METAR, not thoroughly tested),
!    'lcd123' (scatterometer ASCAT), and 'lcd122' (QuickScat) 
!  - 'igpscen' : GPS processing centres assigned as active obs code types,
!                ordered according to their preference in the redundancy check.
!                Default means that no GPS processing centre is used actively.
!  - for new weighting of multiple observations and observation systems:
!    'nwtyp'   : if > 1 then compute net obs. increments for 'nwtyp' different
!                sets of observing systems separately
!    'niwtyp'  : number of observing systems for each set of obs systems
!    'iwtyp'   : defines these sets of observation systems
!    'kwtyp'   : defines the multiple weighting for each of these sets
!  Other changes to the namelist:
!  - Namelist parameters 'lseparw' and 'kmultw' are removed.
!  - The precise meaning and default of 'doromx' are modified (due to a modified
!    meaning of variable 'fdoro').
!  Other changes:
!  - CMA observation type and code type numbers, tables for pressure dependent
!    scales and geometry of horiz. correlations, some I/O device numbers and
!    file names, some general parameters, ODR size, etc are moved away to module
!    'data_obs_lib_cosmo'.
!  - Surface analysis limits and allocatable arrays are moved here from module
!    'data_obs_record'.
!  - Preparation for 1DVAR and the processing of satellite data and production
!    of retrievals. Several new variables are defined which will become namelist
!    variables only when 1DVAR will be included in the official version:
!    - 'l1dvar' : general on/off switch for 1DVar to derive satellite retrievals
!    - 'maxtvo' : ODR size for satellite retrievals
!    - 'gnudgtv' : nudging coefficients for satellite retrievals
!    - 'rhfrtv'  : factor for horizont. weight function for satellite retrievals
!    - 'mcdmsg1', 'mcdmsg2', 'mcdno15', 'mcdno16, 'mcdno17', 'mcdno18' : for use
!                  of new obs types MSG-1 or -2 resp. NOAA-15, -16, -17, or -18.
!    Note: An experimental version for use of satellite data, based on V4_18,
!          is available from christoph.schraff_at_dwd.de .
! V4_25        2012/09/28 Ulrich Schaettler
!  New organizational variables to control output of sfcana-fields
! V4_26        2012/12/06 Ulrich Schaettler
!  Increased icdfdirlen to 250 (length of directory name ycdfdir)
! V4_28        2013/07/12 Christoph Schraff, Ulrich Schaettler
!  New namelist variables:  (CS)
!  - 'qcflbcp'  : factor to threshold used for QC of ps-obs against LBC pressure
!  - 'irun_osse': if feedback file 'fof' present: model run used as obs for OSSE
!  - 'losse_fg' : first guess check flag from 'fof' converted to 'dataset flag'
!  - 'fperturb' : factor to perturb (simulated or original) obs from 'fof' file
!  - 'iseed'    : external seed for random number generator.
!  New namelist variables for near surface analysis (US)
!  - 'ydir_lansfc' : to specify directory where to write the files
!  - 'yform_lansfc': to specify format of the files ('grb1','api1','api2')
! V5_1         2014-11-28 Oliver Fuhrer
!  'gnudgms' and 'lcd146' for Mode-S aircraft data introduced.
!  Replaced ireals by wp (working precision) (OF)
! V5_4         2016-03-10 Christoph Schraff
!  - Removal of AOF interface, therefore removal of namelist variables
!    'itype_obfile', 'yaofpath', 'lpraof', 'dinlat','dislat','diwlon','dielon',
!    'noctrq', 'lgpsbias' and variables 'mrhyes', 'mrhno', 'muvyes', 'muvmo',
!    'yuoafex', 'nuaof', 'nuaofex', 'lobprcs', 'tobnext'. 
!  - New namelist variable: 'yfofdir': directory for output NetCDF 'fof' file.
!  - Initialisation of variables related to the nudging of 1DVAR satellite
!    retrievals (of which the main part has never been included in the official
!    code) moved from 'organize_assimilation' to here.
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!===============================================================================
!
! Modules used:
!
!-------------------------------------------------------------------------------

USE data_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!-------------------------------------------------------------------------------

IMPLICIT NONE

!===============================================================================

! Local Declarations:

!-------------------------------------------------------------------------------
! Global (i.e. public) Declarations:
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Section 1:  General parameters and related variables
!-------------------------------------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    icdfdirlen =250_iintegers ,& ! max. length of name of directory where
                                 !   NetCDF observation input files reside
!   mxda       = 20_iintegers ,& ! max. # data pts 'always' / never to be nudged
    mxgpc      = 20_iintegers ,& ! max. number of GPS processing centres used
    mxwt       = 20_iintegers ,& ! max. number of sets of obs systems for which
                                 ! net obs. increments are computed separately
    mxtyw      = 50_iintegers    ! max. number of obs and code types defining
                                 ! the sets of obs systems

  LOGICAL                    ::       &
    lwonl          ,& ! .TRUE if grid pt. (ionl ,jonl ) lies in the local domain
    lwonl2         ,& ! .TRUE if grid pt. (ionl2,jonl2) lies in the local domain
    lcroot         ,& ! .TRUE if (my_cart_id  == 0) (for print on std. output)
    lvofoi         ,& ! .TRUE if observation increments are also written to VOF
    lfirst = .TRUE.   ! .TRUE if 'organize_nudging' is called the first time

  INTEGER (KIND=iintegers)   ::       &
    onl_cart_id       ! 'my_cart_id'  of node with area containing (ionl,jonl)

  REAL    (KIND=wp)          ::       &
    acthr          ,& ! actual forecast hour
    aiwthr            ! model time [hrs] for which the current temporal weights 
                      ! are valid.(This equals the average time of all timesteps
                      ! in the current 'analysis increment time box' ('AIBT'), 
                      ! i.e. the time box in which the same analysis increments
                      ! are used.)

  LOGICAL                    ::       &
    ltnudge        ,& ! nudging to be done at current timestep
    ltverif           ! verification to be done at current timestep

  REAL    (KIND=wp)          ::       &
    wablua   (4) ,& ! 4*  1. : factor to weights within the ABL
    wablsu   (4)    ! 4*  1. : factor to weights above the ABL
                    !          (atmospheric boundary layer, i.e. mixed layer)

!-------------------------------------------------------------------------------
! Section 2:  Namelist variables controlling the data assimilation
!-------------------------------------------------------------------------------

! For all arrays of length 4: 1: horiz. wind, 2: 'surface' pressure,
!                             3: temperature, 4: (relative) humidity
! (Surface) pressure from TEMP's uses the parameters for upper-air TEMP data;
! T2m, RH2m & UV-10m from TEMP's use  the paremeters for surface-level data.

! ( qdifh    !  1. : reduction factor to horizontal diffusion coefficient )


!      0.    General steering switches
!      -------------------------------

  LOGICAL                    ::       &
    lnudge       ,& ! .f. : on - off switch for nudging
    lverif       ,& ! .f. : on - off switch for verification (writing data to
                    !       the verification observation file VOF as data base)
    lverpas         ! .t. : on - off switch for verif. also of passive reports


!      1.    General variables controlling the nudging
!      -----------------------------------------------

! LOGICAL                    ::       &
!   lobdens         ! .f. : if TRUE then compute obs. density at obs. locations

  INTEGER (KIND=iintegers)   ::       &
    nudgsta      ,& ! 0   : start of nudging period in timesteps
    nudgend      ,& ! 0   : end of nudging period in timesteps
    nversta      ,& ! 0   : start of verification period in timesteps
    nverend      ,& ! 0   : end of verification period in timesteps
    mruntyp      ,& ! -1  : type of current model run used for increments in VOF
                    !       (if -1, then no increments written to VOF or ncfeed)
    mveripr      ,& ! 3   : type of verification/observation file(s) written,
                    !        0 : no file written, equivalent to 'lverif=.false.'
                    !        1 : NetCDF (feedobs/feedbk) file for EnKF or verif.
                    !        2 : ASCII file VOF (YUVERIF)
                    !        3 : both NetCDF and ASCII VOF
    nwtyp        ,& ! 1   : if > 1 then compute net obs. increments for 'nwtyp'
                    !       different sets of observing systems separately
    niwtyp (mxwt),& ! 1,0,0,..: for each of these sets observing systems:
                    !           number of observation or code types which belong
                    !           to that set of observing systems
                    !          (nwtyp = number of non-zero elements in niwtyp+1)
    iwtyp (mxtyw),& ! 0,0,0,..: observation types (for values > 0)
                    !           or     code types (for values < 0, see module
                    !           data_obs_lib_cosmo.f90) belonging to a set of
                    !           obs. systems, specified successively for each
                    !           of these sets in array 'iwtyp';
                    !           '0' denotes all observation and code types
                    !           that are not specified to belong to another
                    !           set of observing systems
                    !           (number of non-zero elements in 'iwtyp' + 1 =
                    !            sum over values in array 'niwtyp')
    kwtyp (mxwt+2)  ! 1   : mode of weights W for multiple observations,
                    !       specified for each set of obseration systems:
                    !       1 : W=  w**2 /sum(w) ;  w= weight of 1 single obs
                    !       2 : W= (w**2 +w) /(1+sum(w))
                    !       (number of non-zero elements in 'kwtyp'
                    !        = 1       if (nwtyp == 1)
                    !        = nwtyp+2 if (nwtyp >= 2) --> in this case, the
                    !          - the first 'nwtyp' entries are used for the
                    !            'nwtyp' specified sets of observing systems,
                    !          specified sets of observing systems, and the
                    !          - the 'nwtyp+1'th entry is used combine the net
                    !            obs. increments from these sets for the final
                    !            analysis increment, and
                    !          - the 'nwtyp+2'th entry is for surface pressure)
                    
  REAL    (KIND=wp)          ::       &
    hversta      ,& ! 0   : start of verification period in 'model integ. hours'
    hverend      ,& ! 0   : end of verification period in 'model integr. hours'
    tconbox         ! 6*dt: timestep [s] for computing analysis increments
                    !       (i.e. time box of constant analysis increments)

!      2.    Corrections to balance the analysis increments
!      ----------------------------------------------------

  LOGICAL                    ::       &
    luvgcor         ! .t. : .t. ==> geostrophic wind correction applied
                    !       (balances near-surface pressure analysis incements)

  INTEGER (KIND=iintegers)   ::       &
    mpsgcor      ,& ! 1   : mode to apply geostrophic pressure correction
                    !       (balances near-surface wind analysis incements
                    !        obtained by isotropic lateral spreading):
                    !       = 0 : no pressure correction
                    !       = 1 : corr. balacing wind from scatterometer only
                    !       = 2 : corr. balacing scatt. + in-situ 10-m wind obs
    ntpscor      ,& ! 1   : switch for temperature correction when the 'surface'
                    !       pressure (pressure at lowest model level) is nudged,
                    !       so that the final geopotential change due to surface
                    !       pressure nudging is zero above the level 'ptpstop'
                    !       = 0 : no temperature correction
                    !       = 1 : cpp(k) = ck(k)^2 * EXP( (1.-ck(k)^3) /8.)
                    !       = 2 : cpp(k) = .5 *ck(k) *(1.+ck(k)) * p(k) / p(ke)
                    !       where  ck(k) =  (p(k) - ptpstop)
                    !                     / (ps - ptpstop)   for p > ptpstop
                    !              ck(k) = 0.                for p < ptpstop
                    !              cpp : pressure correlation betw. levels k, ke
    khumbal         ! 100 : radius (in grid pts. units) of the area around a
                    !       convectively precipitating grid point, at which
                    !       specified instead of relative humidity is preserved
                    !       when nudging temperature.
                    !       If 'khumbal' <= -1 or >= 99, then relative resp. 
                    !       specific humidity is preserved everywhere.
                    !       At the temperature correction (ps-nudging), relative
                    !       humidity is preserved except if 'khumbal' >= 100

  REAL    (KIND=wp)          ::       &
    ptpstop      ,& ! 400.: pressure [hPa] at upper boundary of the temperature
                    !       correction for 'surface' pressure nudging
    qgeo         ,& ! .3  : factor to the geostrophic wind increments at 1000hPa
    qgeotop      ,& ! .5  : factor to the geostrophic wind increments at ptpstop
                    !       (and parabolic interpolation betw. 1000hPa, ptpstop)
    qgeops       ,& ! 1.  : factor to the geostrophic pressure increments

!      3.    Nudging coefficients
!      --------------------------

    gnudg    (4) ,& ! 6,12,6,6*10^-4: nudging coefficients for TEMP / PILOT data
    gnudgar  (4) ,& ! 6, 0,6,0*10^-4: nudging coeffic. for AIRCRAFT data [1/s]
    gnudgms  (4) ,& ! 6, 0,6,0*10^-4: nudging coeffic. for Mode-S aircraft [1/s]
    gnudgsu  (4) ,& ! 6,12,0,6*10^-4: nudging coef. for surface-level data [1/s]
    gnudggp         !        0*10^-4: nudging coefficient for GPS data [1/s]

!      4     Temporal weights
!      ----------------------

  LOGICAL                    ::       &
    ltipol       ,& ! .t. : .t. ==> linear interpolation in time of upper-air
                    !               data which are less than 'tipolmx' hrs apart
                    !           ==> at most 2 reports of same type per station
                    !               used at each timestep
    ltipsu          ! .t. : .t. ==> linear interpolation in time of surface-lev.
                    !               data which are less than 'tipmxsu' hrs apart

  REAL    (KIND=wp)          ::       &
    tipolmx      ,& ! 1.0 : max. time span (hrs) allowed for linear interpolat.
                    !       for upper-air data (set tipolmx = 0, if .NOT ltipol)
    tipmxsu      ,& ! 1.0 : max. time span (hrs) allowed for linear interpolat.
                    !       of surface + GPS data  (tipmxsu = 0, if .NOT ltipsu)
    wtukrsa      ,& ! 3.0 : for saw-tooth shaped temporal weights:
                    !       temporal radius of influence towards the past
                    !       relative to the obs. time     for TEMP / PILOT
    wtukrse      ,& ! 1.0 : for saw-tooth shaped temporal weights:
                    !       temporal radius of influence towards the future
                    !       relative to the obs. time     for TEMP / PILOT
    wtukara      ,& ! 1.5 : for saw-tooth-shaped temporal weights:
                    !       temporal radius of influence towards the past
                    !       relative to the obs. time     for AIRCRAFT data
    wtukare      ,& ! 0.5 : for saw-tooth-shaped temporal weights:
                    !       temporal radius of influence towards the future
                    !       relative to the obs. time     for AIRCRAFT data
    wtuksua      ,& ! 1.5 : for saw-tooth shaped temporal weights:
                    !       temporal radius of influence towards the past
                    !       rel. to the obs. time   for surface-level + GPS data
    wtuksue         ! 0.5 : for saw-tooth shaped temporal weights:
                    !       temporal radius of influence towards the future
                    !       rel. to the obs. time   for surface-level + GPS data

!      5.    Spatial weights : Spreading of observational information
!      --------------------------------------------------------------

!      5.1   Mode of spreading
!      -----------------------

  INTEGER (KIND=iintegers)   ::       &
    msprpar      ,& ! 1   : switch specifying the surfaces along which observat.
                    !       increments of upper-air data are (primarily) spread
                    !       0: spreading along model levels, vertical weights
                    !          depend approx. on differences in log( pressure )
                    !          ('ln(p)', i.e. exactly: scaled height) between
                    !          the point assigned to the obs. increment ('POI')
                    !          and the target level at the observation location
                    !       1: spreading along horizontal surfaces, vertical
                    !          weights depend approx. on 'ln(p)' differences
                    !          between the 'POI' and the target grid point
                    !       2: spreading along isentropic surfaces, vertical
                    !          weights depend potential temperature differences
                    !          between the 'POI' and the target grid point
    msprpsu         ! 0   : switch specifying the surface along which surface-
                    !       level data increments are primarily spreaded
                    !       0: spreading along model levels        \  defined
                    !       1: spreading along horizontal surfaces  > as for
                    !       2: spreading along isentropic surfaces /  'msprpar'
                    !       (here, the 'POI' is always a grid point at the
                    !        lowest full model level)

!      5.2   Vertical weights
!      ----------------------

  REAL    (KIND=wp)          ::       &
!      for upper-air data
    vcorls   (4) ,& ! 2*.333,: square of the vertical correlation scale,
                    ! 2*.04    i.e. of the Gaussian vertical influence 'radius'
                    !          in log( pressure ) if (msprpar <= 1), or
                    !          in potential temperature if (msprpar == 2)
                    !          (reasonable values in latter case: 2*275., 2*33.)
    vcutof   (4) ,& ! 2* .75,: cut-off of the vertical correlation.       Units:
                    ! 2*1.     value of correlation at cut-off is [exp(-vcutof)]
                    !          (c.f UKMO: cut-off at: |ln(p/pobs))| = 0.5
                    !                     correl. fn:  exp(-3 *(ln(p/pobs))**2 )
                    !               ==> value of correl. at cut-off: exp(-.75) )
                    !          (i.e. for (msprpar <= 1) and default values for
                    !                'vcorls', 'vcutof', the cut-off 'radius'
                    !                for wind is 0.5 (density) scale heights.)
!   wablua   (4) ,& ! 4*  1. : factor to weights within the ABL
                    !          (atmospheric boundary layer, i.e. mixed layer)
!      for surface-level data
    vcorlsu  (4) ,& ! 2*.013,: square of the vertical correlation scale,
                    ! 2*.002   i.e. of the Gaussian vertical influence 'radius'
                    !          in log( pressure ) if (msprpsu == 1) or
                    !          in potential temperature if (msprpsu == 2)
                    !          (e-folding decay height of ca. 300m for T, RH, if
                    !          default, or if 2*11.1, 2*1.33 for (msprpsu == 2))
    vcutosu  (4) ,& ! 2* .75,: cut-off of the vertical correlation.       Units:
                    ! 2*4.     value of correlation at cut-off : [exp(-vcutosu)]
    vpblsu   (4)    ! 2*99. ,: Gaussian vertical influence 'radius' of potential
                    ! 2*99.    temperature differences between obs. level and
                    !          model level, to render vertical weights depend on
                    !          the stability (in the PBL) even if (msprpsu <= 1)
!   wablsu   (4)    ! 4*  1. : factor to weights above the ABL
                    !          (atmospheric boundary layer, i.e. mixed layer)

  LOGICAL                    ::       &
    lsvcorl         ! .t. : .t. ==> diminishing of vertical correlation scales
                    !               in the presence of close observations

!      5.3   Horizontal weights
!      ------------------------

  REAL    (KIND=wp)          ::       &
!      for upper-air data
    rhinfl   (4) ,& ! 0.,70.,: constant part of the 'correlation scale of the
                    ! 0., 0.   autoregressive horiz. correlation function'
                    !          (='COSAC') [km]
    rhvfac   (4) ,& ! 1., 0.,: multiplication factor to the vertically varying
                    ! 2* .83   part of the 'COSAC', given on following p-levels:
                    !            1000,850,700,500,400,300,250,200,150,100, 50hPa
                    !          uv: 70, 80, 90,100,100,110,115,120,125,125,125 km
                    !          Tq: 70, 80, 90,100,100,100,100,110,120,120,120 km
    rhtfac   (4) ,& ! 1.3,   : scaling factor of the total 'COSAC' for the
                    ! 1.43,    beginning and end of the nudging period for 1 obs
                    ! 1.3,     relative to the 'COSAC' at the obs. time as given
                    ! 1.3      by 'rhinfl', 'rhvfac')
    rhfgps       ,& ! 0.45   : scaling (reduction) factor of the total 'COSAC'
                    !          for humidity derived from GPS IWV
    rhfpsdd      ,& ! 1.0    : minimum scaling (reduction) factor of the total
                    !          'COSAC' for surface pressure dep. on data density
    cutofr   (4) ,& ! 4* 3.5 : cut-off in 'COSAC' units of the horizontal
                    !          correlation function
    vcsni    (4) ,& ! 4*2500.: square of Gaussian vertical influence 'radius'
                    !          in potential temperature (if msprpar <= 1) or
                    !          log( pressure ) (if msprpar == 2) on surfaces
                    !          along which the obs. increments are spread
                    !          laterally. This determines the non-isotropy of
                    !          the (1-d) horizontal correlation functions.
!      for surface-level data
    rhiflsu  (4) ,& ! 2 *70.,: constant part of the 'correlation scale of the
                    !   100.,  autoregressive horiz. correlation function'
                    !    70.   (='COSAC') [km]  (at the obs. time)
    rhtfsu   (4) ,& ! 1.,    : scaling factor of the total 'COSAC' for the
                    ! 1.43,    beginning and end of the nudging period for 1 obs
                    ! 1.,      relative to the 'COSAC' at the obs. time as given
                    ! 1.       by 'rhiflsu')
    cutofsu  (4) ,& ! 2., 3.5: cut-off in 'COSAC' units of the horizontal
                    ! 2., 2.   correlation function
    vcsnisu  (4) ,& ! 2*2500.: square of Gaussian vertical influence 'radius'
                    ! 2*   9.  in potential temperature (if msprpsu <= 1) or
                    !          log( pressure ) (if msprpsu == 2) on surfaces
                    !          along which the obs. increments are spread
                    !          laterally. This determines the non-isotropy of
                    !          the (1-d) horizontal correlation functions.
!      degree of non-divergence of 2-dim. horizontal wind correlation functions
    cnondiv      ,& ! .1  : constant part of the 'factor to the non-divergence
                    !       correction in the 2-dim wind correlations' (='FNDC')
    fnondiv      ,& ! .8  : multiplication factor to the vertically varying part
                    !       of the 'FNDC', given on following pressure levels:
                    !       hPa : 1000,850,700,500,400,300,250,200,150,100, 50
                    !       V-FNDC: .4, .5, .5, .5, .5, .6,.65, .7,.75,.75,.75
    tnondiv         ! 1.1 : temporal factor to the 'FNDC' for the beginning and
                    !       end of the nudging period for 1 obs relative to that
                    !       given by 'cnondiv', 'fnondiv' for the obs. time

!      6.    Computation of observation increments
!      -------------------------------------------------

!      6.1.  Vertical profiles of observation increments
!      -------------------------------------------------

  LOGICAL                    ::       &
    lscadj   (4)    ! .T.,   : .F. ==> linear vertical interpolation (in log(p))
                    ! .T.,     instead of vertical scale adjustment (by vertical
                    ! .T.,     averaging over the model layer) for conveying the
                    ! .F.      observational information to the model levels
                    !          (for computing obs. increments at model levels)

  REAL    (KIND=wp)          ::       &
    topobs   (4) ,& !  849., : threshold [hPa]: above this level (p < topobs),
                    ! 1099.,   only obs. increments at model levels are used,
                    !  799.,   i.e. obs. increments at obs. levels are not used
                    !  699.    ('topobs' is fixed at 1099. if (msprpar == 0))
    botmod   (4)    ! 1099., : threshold [hPa]: below this level (p > botmod),
                    ! 1099.,   only obs. increments at obs. levels are used, i.e
                    ! 1099.,   obs. increments at model levels are not computed
                    !  899.,   ((botmod >= topobs), and
                    !           'botmod' is fixed at 1099. if (msprpar == 0))

!      6.2   Surface-level observation increments
!      ------------------------------------------

  LOGICAL                    ::       &
    loiqv2m      ,& ! .f. : .t. ==> 2-m humidity observation increments as
                    !               differences of specific humidity instead of
                    !               relative humidity
    lqfqv2m         ! .f. : .t. ==> quality factor for 2-m humidity observations
                    !               dependent on T-2m differences

!      7.    Quality control and quality weights
!      -----------------------------------------

!      7.1   (Threshold) Quality control
!      ---------------------------------

  REAL    (KIND=wp)          ::       &
    dtqc         ,& ! 720.: timestep (in [s]) for the threshold quality control
                    !       (the quality ctrl is always applied to observations
                    !        at the first time when they are used)
                    !       Note: to ensure that all obs are written to VOF or
                    !             NetCDF FOF as required, 'dtqc' should be a
                    !             divisor of min( hverend, hstop )  !
!      for upper-air data
    qcc      (4) ,& !  0.,500: constant part of the 'quality control thresholds'
                    !  0., .7  (='QCT'). (u,v):[m/s], ps: [Pa], T: [k], RH: [ ]
    qcvf     (4) ,& !  5., 1.: multiplication factor to the vertically varying
                    ! 10., 0.  part of the QCT (for height/ thickness instead of
                    !          pressure ps, and not available for humidity RH),
                    !          given on following pressure levels:
                    !          level          ,hPa:1000, 850, 700, 500, 400, 300
                    !               ,250, 200, 150, 100,  70,  50,  30,  20,  10
                    !          TEMP  wind     ,m/s: 2.3, 2.3, 2.5, 3.0, 3.5, 3.7
                    !               ,3.5, 3.5, 3.4, 3.3, 3.2, 3.2, 3.3, 3.6, 4.5
                    !          AIREP wind     ,m/s: 2.5, 2.5, 3.0, 3.5, 4.0, 4.0
                    !               ,4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0
                    !          TEMP  temperature,K: 1.2, 1.0, 0,7, 0.4, 0.4, 0.5
                    !               ,0.5, 0.6, 0,7, 0.8, 0.8, 0,9, 0.9, 1.0, 1.2
                    !          AIREP temperature,K: 1.2 ,1.0, 0.7, 0.5, 0.5 ,0.6
                    !               ,0.6 ,0.7, 0.8, 0.9, 1.0, 1.1, 1.1, 1.2, 1,4
                    !          (cf. 'data_obs_qc_limits' for tables, thickness
                    !           QC and other factors used for thresholds, e.g.
                    !           an additional factor of 0.8 used for AIREP wind)
!      for surface-level data
    qccsu    (4) ,& ! 12.,500: constant parts of the quality control thresholds
                    ! 12., .7  (='QCT'). (u,v):[m/s], ps: [Pa], T: [k], RH: [ ]
!      for integrated water vapour (IWV) derived from GPS or radiosonde data
    qcciq        ,& ! 1.  : constant part of QC threshold for IWV
    qcsiq        ,& ! .15 : IWV QC threshold, as a fraction of IWV of saturated
                    !       model temperature profile
    qcflbcp      ,& ! 1.4 : enhancement factor to the threshold used for the
                    !       check of ps against lateral boundary fields
                    !       (if (qcflbcp <= 0) then no LBC check for ps)

!      7.2   Quality weights
!      ---------------------

    qcfpst          ! 1.5 : maximum enhancement of the (quality) weight for
                    !       surface pressure observations due to the pressure
                    !       tendency (the max. enhancement is used for observa-
                    !                 tions with (an absolute value of) observed
                    !                 3-hourly tendency greater than 25 hPa; the
                    !                 enhancement is always zero for tendencies
                    !                 less than 3 hPa, and is linear in between)
                    ! (for 'doromx', see section 8.1)

!      8.    Observation processing
!      ----------------------------

!      8.0   Reading of observation reports
!      ------------------------------------

  INTEGER (KIND=iintegers)   ::       &
!   itype_obfile ,& !   1    : type of observation input file(s):
                    !          1: AOF ,  2: NetCDF
    irun_osse       !   0    : for data from input feedback file yfofin='fof':
                    !          == 0 : obs from yfofin are used as obs
                    !          >= 1 : simulated obs from yfofin with model run
                    !                 index = irun_osse are used as obs (thus,
                    !                 run 'irun_osse' is nature run for OSSE)

  LOGICAL                    ::       &
    losse_fg        !  .f.   : if true then first guess check flag from 'fof' is
                    !               converted into 'dataset' pre-processing flag
                    !          if false then first guess check flag is discarded
                    !               and obs may be used actively

  CHARACTER (LEN=icdfdirlen) ::       &
    ycdfdir      ,& ! './'   : directory where NetCDF observation input files
                    !            and the blacklist file reside
    yfofdir         ! './'   : directory where NetCDF 'fof' feedobs (feedback)
                    !            file(s) are written

!      8.1   Use of stations / reports
!      -------------------------------

  REAL    (KIND=wp)          ::       &
    obnlat       ,& !   90.  : northern boundary of observation area
    obslat       ,& !  -90.  : southern boundary of observation area
    obwlon       ,& ! -180.  : western boundary of observation area
    obelon       ,& !  180.  : eastern boundary of observation area
    exnlat       ,& !   90.  : northern boundary for exclusion area
    exslat       ,& !  -90.  : southern boundary for exclusion area
    exwlon       ,& ! -180.  : western boundary for exclusion area
    exelon       ,& !  180.  : eastern boundary for exclusion area
    doromx   (4) ,& !  100., : vertical extrapolation cut-off and gaussian
                    !  150.,   radius of height differences between model
                    !  150.,   orography and surface station height for a factor
                    !  150.,   contributing to the quality weight factor as part
                    !          of the nudging weights
                    !          (height diff of SYNOP/ GPS with (z-obs > z-model,
                    !           --> interpolation instead of extrapolation) are
                    !           divided by 4 (fdoro) for surf pressure/ IWV obs)
    altopsu  (4) ,& !  100., : SYNOP obs. above height 'altopsu' are not assimi-
                    ! 3*5000.  lated. If (altopsu == 0.) then SYNOP / surf. TEMP
                    !          assigned to land grid pts. are not assimilated
    thairh       ,& !   20.  : maximum horizontal distance [km] between the
                    !          lowest report and any single level report that
                    !          is added to a multi-level AIRCRAFT report
    fperturb        !    0.  : factor to the obs error variances to define the
                    !          size of random perturbations added to the obs
                    !          (only for data from feedback file yfofin='fof')

  INTEGER (KIND=iintegers)   ::       &
                    !max. # reports currently active (i.e within the time window
                    !given by 'wtuk??a','wtuk??e','tip*mx*') in the total domain
    maxmlo       ,& !  300   : max. number of multi-level reports 
    maxsgo       ,& ! 3000   : max. number of (surface-level and upper-air)
                    !                         single-level reports 
    maxuso       ,& !  900   : max. number of upper-air single-level reports
    maxgpo       ,& ! 3000   : max. number of GPS reports within total domain
    maxmlv       ,& !  100   : max. number of observation levels in multi-level
                    !                         reports
    mxfrep       ,& !   -1   : max. number of reports in NetCDF feedobs file
    mxfobs       ,& !   -1   : max. number of observations in feedobs file
                    !          (if mxfrep, mxfobs <= 0 then reasonable values
                    !           are computed from maxmlo, maxsgo, hverend etc.)
    nolbc        ,& !    5   : number of grid rows at lateral boundaries
                    !          where obs are neglected
    mqcorr92     ,& !    0   : switch for bias correction for Vaisala RS92
                    !            radiosonde humidity
                    !          = 0 : no correction for humidity
                    !          = 1 : correct only solar radiation bias
                    !          = 2 : correct total bias (incl. nighttime bias)
    iseed           !    0   : external seed for random number generator
                    !          (this seed is combined with another seed that
                    !           depends on the initial date and time of the
                    !           model run;  note that for LETKF OSSE, 'iseed'
                    !           should be identical for all ensemble members
                    !           so that the perturbed obs will be identical)

  INTEGER (KIND=iintegers)   , TARGET   ::       &
    igpscen (mxgpc) ! X* -1  : array of processing centres of GPS reports used
                    !          actively (order of centres determines preference
                    !          in redundancy check; '-1' means no active centre)

!-------------------------------------------------------------------------------
!      8.2   Use of observation and code types
!      ---------------------------------------

  LOGICAL                    ::       &
    lsynop       ,& ! .t.    : .t. if SYNOP data is used
    laircf       ,& ! .t.    : .t. if AIREP data is used (aircraft)
    lsatob       ,& ! .false.: .t. if SATOB data is used
    ldribu       ,& ! .t.    : .t. if BUOY  data is used (drifting buoy)
    ltemp        ,& ! .t.    : .t. if TEMP  data is used
    lpilot       ,& ! .t.    : .t. if PILOT data is used
    lsatem       ,& ! .false.: .t. if SATEM data is used
    lgps         ,& ! .false.: .t. if GPS   data is used
    lscatt          ! .t.    : .t. if SCATT data is used (scatterometer)

  LOGICAL                    ::       &
    lcd011       ,& ! .t.    : .t. if synop code  11 data is used (land synop)
    lcd014       ,& ! .t.    : .t. if synop code  14 data is used (automatic)
    lcd021       ,& ! .t.    : .t. if synop code  21 data is used (ship)
    lcd022       ,& ! .t.    : .t. if synop code  22 data is used (ship abbrev.)
    lcd023       ,& ! .t.    : .t. if synop code  23 data is used (shred)
    lcd024       ,& ! .t.    : .t. if synop code  24 data is used (autom. ship)
    lcd140       ,& ! .t.    : .t. if synop code 140 data is used (metar)
    lcd041       ,& ! .t.    : .t. if airep code  41 data is used (codar)
    lcd141       ,& ! .t.    : .t. if airep code 141 data is used (airep)
    lcd241       ,& ! .t.    : .t. if airep code 241 data is used (colba)
    lcd144       ,& ! .t.    : .t. if airep code 144 data is used (amdar)
    lcd146       ,& ! .t.    : .t. if airep code 146 data is used (mode-s)
    lcd244       ,& ! .t.    : .t. if airep code 244 data is used (acars)
    lcd088       ,& ! .t.    : .t. if satob code  88 data is used (satob)
    lcd188       ,& ! .f.    : .t. if satob code 188 data is used (sst)
    lcd063       ,& ! .t.    : .t. if dribu code  63 data is used (bathy)
    lcd064       ,& ! .t.    : .t. if dribu code  64 data is used (tesac)
    lcd165       ,& ! .t.    : .t. if dribu code 165 data is used (drift. buoy)
    lcd035       ,& ! .t.    : .t. if temp  code  35 data is used (land temp)
    lcd036       ,& ! .t.    : .t. if temp  code  36 data is used (temp ship)
    lcd037       ,& ! .t.    : .t. if temp  code  37 data is used (mobile)
    lcd135       ,& ! .t.    : .t. if temp  code 135 data is used (dropsonde)
    lcd039       ,& ! .t.    : .t. if temp  code  39 data is used (rocob)
    lcd040       ,& ! .t.    : .t. if temp  code  40 data is used (rocob ship)
    lcd032       ,& ! .t.    : .t. if pilot code  32 data is used (land pilot)
    lcd033       ,& ! .t.    : .t. if pilot code  33 data is used (pilot ship)
    lcd038       ,& ! .t.    : .t. if pilot code  38 data is used (mobile)
    lcd132       ,& ! .t.    : .t. if pilot code 132 data is used (win-prof eu)
    lcd133       ,& ! .t.    : .t. if pilot code 133 data is used (sod/rass eu)
    lcd136       ,& ! .t.    : .t. if pilot code 136 data is used (pro/rass us)
    lcd137       ,& ! .t.    : .t. if pilot code 137 data is used (Radar VAD)
    lcd086       ,& ! .t.    : .t. if satem code  86 data is used (satem)
    lcd186       ,& ! .t.    : .t. if atovs code 186 data is used (hi-res ATOVS)
    lcd122       ,& ! .t.    : .t. if scatt code 122 data is used (QuickScat)
    lcd123       ,& ! .t.    : .t. if scatt code 123 data is used (ASCAT)
    lcd096          ! .t.    : .t. if gps data from COST ASCII file is used

!      9.    2-D analyses
!      ------------------

  LOGICAL                    ::       &
    lsurfa       ,& ! .f.    : .t. if surface fields are analysed
    lt2m         ,& ! .f.    : .t. if 2m temperat. field is analysed
    lrh2m        ,& ! .f.    : .t. if 2m rel. hum. field is analysed
    lprecp       ,& ! .f.    : .t. if precipitation is analysed
    lff10m          ! .f.    : .t. if 10m wind speed is analysed

  REAL    (KIND=wp)          ::       &
    ht2a         ,& ! 999.   : time of 1. T2m-ana in hours since model start
    ht2i         ,& ! 999.   : time increment to next T2m analysis
    ht2next      ,& ! 999.   : next hour, when T2m-ana shall be done
    hh2a         ,& ! 999.   : time of 1. RH2m-ana in hours since model start
    hh2i         ,& ! 999.   : time increment to next RH2m analysis
    hh2next      ,& ! 999.   : next hour, when RH2m-ana shall be done
    hffa         ,& ! 999.   : time of 1. 10m wind-ana in hours since model start
    hffi         ,& ! 999.   : time increment to next wind analysis
    hffnext      ,& ! 999.   : next hour, when wind-ana shall be done
    hprc         ,& ! 999.   : time of prec-ana in hours since model start
    hprcnext     ,& ! 999.   : time of prec-ana in hours since model start
    raintp          ! 12.    : time period of precipitation analysis

INTEGER (KIND=iintegers)     ::       &
    nt2next      ,& !        : next time step, when T2m-ana shall be done
    nh2next      ,& !        : next time step, when RH2m-ana shall be done
    nffnext      ,& !        : next time step, when wind-ana shall be done
    nprcnext        !        : next time step, when prec-ana shall be done

  CHARACTER (LEN=250)        ::       &
    ydir_lansfc     ! './'   : directory where to write the 2-D analyses

  CHARACTER (LEN=  4)        ::       &
    yform_lansfc    ! 'grb1' : format for the 2-D analyses files

!      10.   Diagnostic output
!      -----------------------

  LOGICAL                    ::       &
    lprodr       ,& ! .t.    : .t. for diagnostic print of obs data records ODR
    ldiasa          ! .f.    : .t. for diagnostics of surface analysis

  INTEGER (KIND=iintegers)   ::       &
    ionl         ,& ! 167    : / grid point coordinates
    jonl         ,& ! 103    : \ for standard output on nudging
    ionl2        ,& ! 167    : / 2nd grid pt coordinates
    jonl2           ! 103    : \ for other standard output on nudging


!-------------------------------------------------------------------------------
! Section 2b:  Variables controlling the 1DVAR for nudging of satellite
!              radiances. They should have become namelist variables as soon as
!              the experimental and private 1DVAR code would be included
!              in the official version, but this will never be the case.
!-------------------------------------------------------------------------------

  LOGICAL                    ::       &
    l1dvar = .FALSE.! .f. : on - off switch for 1DVar
                    !       (if switch is off, then 1DVar is not performed
                    !        irrespective of the values of other parameters,
                    !        and no additional input files are expected
                    !        (e.g. on B-matrix, stratospheric profiles, channel
                    !         selection, bias correction, RTTOV coefficients))
                    ! (Note: no 1dvar is implemented in this official verion !)

  REAL    (KIND=wp)          ::       &
    gnudgtv (4) = (/0.0000_wp, 0.0000_wp  ,&  ! nudging coefficients for 1DVAR
                    0.0006_wp, 0.0006_wp/),&  ! satellite retrievals 
    rhfrtv  (4) = (/1.0_wp, 1.0_wp,&  ! scaling factor of the total 'COSAC' for
                    1.0_wp, 0.5_wp/)  ! temperature / humidity sat. retrievals

  INTEGER (KIND=iintegers)   ::       &
                    !max. # reports currently active (i.e within the time window
                    !given by 'wtuk??a','wtuk??e','tip*mx*') in the total domain
    maxtvo = 1      !    1   : max. number of sat retrievals within total domain
                    !                (however, sat retrievals cannot be produced
                    !                 with current official version !)

  INTEGER (KIND=iintegers)   ::       &
    ! Note: Satellite data cannot yet be produced with current official
    !       version, i.e. the following variables have no effect.
    !       An experimental version for use of satellite data, based on 
    !       V4_18, is available from christoph.schraff_at_dwd.de .
                    ! mcdxxxx = 0 : xxxx is not processed
                    ! mcdxxxx = 1 : xxxx is processed, but not used actively
                    ! mcdxxxx = 2 : xxxx is used actively (and processed)
                    ! --> if (mcdxxxx >= 1) then for satellite xxxx, 3 input
                    !     files are stringently expected (for bias correction,
                    !     channel selection, and RTTOV coefficients)
    mcdmsg1 = 0  ,& !  0     : processing / use of MSG1   code  71 data
    mcdmsg2 = 0  ,& !  0     : processing / use of MSG2   code  72 data
    mcdno15 = 0  ,& !  0     : processing / use of NOAA15 code 206 data
    mcdno16 = 0  ,& !  0     : processing / use of NOAA16 code 207 data
    mcdno17 = 0  ,& !  0     : processing / use of NOAA17 code 208 data
    mcdno18 = 0     !  0     : processing / use of NOAA18 code 209 data

!-------------------------------------------------------------------------------
! Section 3:  I/O device numbers for nudging and file names
!-------------------------------------------------------------------------------

!      3.1  File names
!           ----------  

  CHARACTER (LEN=7)        , PARAMETER  :: &
!   yuaof                   ! see 'yaofpath':  AOF (input)
    yugps   = 'gps    '     ! input GPS observation file
!   yuaofex = 'YUAOFEX'     ! (output:) expanded AOF
!   yucautn = 'YUCAUTN'  ,& ! caution messages if too many obs for ODR size
!   yuquctl = 'YUQUCTL'  ,& ! data rejected by threshold QC at obs time
!   yurejct = 'YUREJCT'  ,& ! direct reporting of rejected obs. reports
!   yustats = 'YUSTATS'  ,& ! statistics of processed reports
!   yuobsdr = 'YUOBSDR'  ,& ! observations stored in the observation data record
!   yuverif = 'YUVERIF'  ,& ! VOF (output): verification observation file (obs.
                            !      incl. quality control flag for verification)
!   yuprint = 'YUPRINT'     ! all the remaining information

!      3.2  Device numbers
!           --------------

  INTEGER (KIND=iintegers) :: &
    nsfc    , & !
    nugps       ! 75 ! GPS observation file unit
!   nuaof   , & ! 65 ! AOF (input): analysis observation file
!   nusatin , & ! 76 ! satellite observations / parameter input files
!   nucautn , & ! 77 ! caution messages if too many obs for current ODR size
!   nuaofex     ! 66 ! (output:) expanded AOF
!   nuqc    , & ! 67 ! data rejected by threshold quality control at obs time
!   nurej   , & ! 68 ! direct reporting of rejected obs. reports
!   nustat  , & ! 69 ! statistics of processed reports
!   nuodr   , & ! 70 ! observations stored in the observation data record
!   nuverif , & ! 71 ! VOF (output): verification observation file (observations
                     !   incl. quality control flag for verification)
!   nupr        ! 64 ! all the remaining information

!-------------------------------------------------------------------------------
!  Section 4 : Surface analysis limits and allocatable arrays
!-------------------------------------------------------------------------------

!       4.1    Surface analysis limits
!              -----------------------

  REAL (KIND=wp)           , PARAMETER  :: &
    rraint(2) = (/  1.0_wp,  24.0_wp/) ,& ! rain period limits
    rrainl    =   100.0_wp             ,& ! amount of rain limit
    rmaxdp(4) = (/ 14.0_wp,   0.8_wp   ,& ! quality control thresholds
                   13.0_wp, 350.0_wp/)    ! for observation departures of
                                          !   t_2m, rh_2m, ff_10m, prec

!       4.2    Allocable arrays used for surface analysis
!              ------------------------------------------

  REAL    (KIND = wp)       , ALLOCATABLE :: &
    alat     (:)   ,& ! model latitudes
    alon     (:)   ,& ! model longitudes
    anal     (:,:) ,& ! global analysis field
    firstg   (:,:) ,& ! global first guess field
    rsurin   (:,:) ,& ! increment field
    deprow   (:,:) ,& ! observations or observed increments for each model row
    rpalto   (:,:) ,& ! original station latitude
    rpalno   (:,:) ,& ! original station longitude
    rpalat   (:,:) ,& ! station latitude (assigned model g.p.)
    rpalon   (:,:) ,& ! station longitude (assigned model g.p.)
    rcslon   (:,:) ,& ! cosine of station longitude
    rsnlon   (:,:) ,& ! sine of station longitude
    rcslat   (:,:) ,& ! cosine of station latitude
    rsnlat   (:,:) ,& ! sine of station latitude
    rpahgt   (:,:) ,& ! station height
    rpaobe   (:,:) ,& ! observation error
    wa       (:)   ,& ! weight * observation increment
    wwa      (:)      ! weights

  INTEGER (KIND = iintegers), ALLOCATABLE :: &
    bufall   (:,:) ,& ! buffer array containing observations from all PEs
    noatab   (:)   ,& ! index of the first obs in a model row
    indexv   (:,:) ,& ! index describing the state of the analysis at a g.p.
    nnodep   (:)   ,& ! number of observations per row
    npagpt   (:,:) ,& ! sea/land type of assigned grid point
    nuprow   (:)   ,& ! northern row of local data selection area
    ndnrow   (:)      ! southern row of local data selection area

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    ilstid_sfc  =  9  ! character length of the station identity

  CHARACTER (LEN = ilstid_sfc)  , ALLOCATABLE :: &
    yidall   (:)   ,& ! buffer of station ids of all observations from all PEs
    ypasid   (:,:)    ! station id


!-------------------------------------------------------------------------------
!  Section 5:  Miscellany
!-------------------------------------------------------------------------------

!      5.0  Full model fields on Arakawa A grid
!           -----------------------------------

  REAL    (KIND = wp)       , ALLOCATABLE :: &
    a_u (:,:,:),    & ! zonal wind speed            on Arakawa A grid ( m/s )
    a_v (:,:,:),    & ! meridional wind speed       on Arakawa A grid ( m/s )
    a_p (:,:,:),    & ! pressure (full value)       on main levels    ( Pa  )
    a_z (:,:,:)       ! geometrical height          of main levels    (  m  )

!      5.1  Switch for observation processing
!           ---------------------------------

  LOGICAL                    ::       &
!   lobprcs    ,& ! .TRUE if observation processing at current timestep
!   ltvprcs    ,& ! .TRUE if 1dvar       processing at current timestep
    lexigps       ! .TRUE if GPS input file exists

  REAL    (KIND=wp)          ::       &
    tmaxbox    ,& ! maximum interval [s] for computing analysis increments
    cqcwbox       ! (maximum) half interval [s] within which observations are
                  ! quality controlled and written on the feedobs file

!      5.2  Areal arrays
!           ------------

  REAL    (KIND=wp)        , ALLOCATABLE :: &
    zconbas  (:,:) ,& ! height of base of convectively instable region
    zcontop  (:,:)    ! height of top  of convectively instable region

!      5.3  Variables
!           ---------

  REAL    (KIND=wp)          ::       &
    qctfpr    (8) ! for VOF output only: time factor for quality cntl thresholds

  INTEGER (KIND=iintegers)   ::       &
    isetyp0    ,& ! index in 'niwtyp' which points to the first index with
                  ! observation type '0' in 'iwtyp' (denoting the set of obs
                  ! systems with contains all ('remaining') obs types that are
                  ! not specified explicitly in 'iwtyp')
    isetyp(mxtyw) ! defines for each observation or code type in 'iwtyp', which
                  !   set of observing systems (index of 'niwtyp') it belongs to

!      5.4  Variables that may become namelist parameters in a later version
!           ----------------------------------------------------------------

  LOGICAL                    ::       &
    liwvssc=.true.  ! .t. : spatial consistency check of IWV performed

!      5.4.1 Variables related to 1DVar
!            --------------------------

  REAL    (KIND=wp)        , PARAMETER  :: &
    cbot    = 109900.0_wp        !

  REAL    (KIND=wp)          ::       &
    topobtv  (4)  & ! 1099., : threshold [hPa]: above this level (p < topobs),
      = (/ cbot,  & ! 1099.,   only obs. increments at model levels are used,
      cbot,cbot,  & ! 1099.,   i.e. obs. increments at obs. levels are not used
      cbot /)    ,& ! 1099.    ('topobs' is fixed at 1099. if (msprpar == 0))
    botmotv  (4)  & ! 1099., : threshold [hPa]: below this level (p > botmod),
      = (/ cbot,  & ! 1099.,   only obs. increments at obs. levels are used, i.e
      cbot,cbot,  & ! 1099.,   obs. increments at model levels are not computed
      cbot /)       ! 1099.,   ((botmod >= topobs), and
                    !           'botmod' is fixed at 1099. if (msprpar == 0))

! LOGICAL                    ::       &
!   lrhBcorr = .FALSE.  ! .f.: .true. if RH correction to B matrix applied 

!-------------------------------------------------------------------------------

END MODULE data_nudge_all
