!+ Source module for the observation processing in the data assimilation mode
!-------------------------------------------------------------------------------

MODULE src_obs_cdfin_mult

!-------------------------------------------------------------------------------
! Description:
!   This module performs the observation pre-processing of multi-level reports
!   read from NetCDF observation input files. The current types of multi-level
!   reports comprise of radiosonde (TEMP and PILOT) and ground-based remote-
!   sensing profiler (Wind Profiler, RASS, Radar VAD) reports.
!   Special tasks are:
!    - reading reports from NetCDF files, selecting, assigning them to model
!      grid points and distributing them to the nodes (precessors) according
!      to the sub-domain which they lie on
!    - pre-processing, including some gross error checking, blacklist checking,
!      assigning observation errors, static bias correction, etc.
!    - storing the reports in the internal ODR arrays
!
!   Note: This module is part of the 'COSMO data assimilation library 1'
!         for reading data from NetCDF observation input files.
!         It is used commonly by the COSMO model and the 3DVAR program package !
!
! Method:
!   This module contains the following module procedures:
!    - obs_cdf_read_temp_pilot: (called by obs_cdf_read_org)
!    - obs_cdf_read_profiler  : (called by obs_cdf_read_org)
!    - obs_cdf_store_multilev :  called by obs_cdf_read_temp_pilot
!                                      and obs_cdf_read_profiler
!
!   This module also contains elemental functions, formerly statement functions:
!   - insert       : inserts bits set in one integer into another integer word
!                    at given bit position
!   - ibit1        : returns 1 bit at given bit position of given integer word
!
!   It uses from:
!    - src_obs_cdfin_comhead: - obs_cdf_read_comhead
!                             - obs_cdf_buffer_comhead
!                             - obs_cdf_store_comhead
!    - src_obs_cdfin_blk:     - obs_cdf_whitelist_local
!                             - obs_cdf_blacklist_local
!    - src_obs_cdfin_util:    - obs_assign_sort_node
!                             - obs_cdf_distrib_reports
!                             - std_atmosphere
!                             - obs_td2rh
!                             - obs_find_level
!                             - f_z2p
!    - utilities:             - phi2phirot
!                             - rla2rlarot
!                             - uv2uvrot_vec
!    - parallel_utilities:    - distribute_values
!    - environment:           - model_abort
!
!   Data modules used:
!    - data_parameters
!    - data_obs_lib_cosmo
!    - data_obs_cdfin 
!    - data_obs_record
!    - mo_fdbk_tables
!
! Current Code Owner (for COSMO and for DWD 3DVAR):
!  DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V4_22        2012/01/31 Christoph Schraff
!  Initial release, extracted from module 'src_obs_proc_cdf' and adapted (e.g.
!    modified routine interfaces and diagnostic arrays, 'obs_pointrs'-->'i_cma',
!    modified ODR (cloud group words, observation status word)).
!  - Call of external routine 'atmhp' (libmisc) replaced by 'std_atmosphere'.
!  - To make observation input more flexible, some variables in the NetCDF
!    files are changed from mandatory to optional.
!  - Observed 'temperature variable' (e.g. virtual temperature) not converted
!    into temperature any more.
!  - Bug correction at selection of surface observations (usage of 'fdoro').
!  - Bug correction allowing for active use of pressure obs from TEMP surface
!    level (in multi-level report only).
!  - Pressure derived from height flagged in main flag word.
!  - Screen-level obs within multi-level reports set to passive (and flagged).
! V4_28        2013/07/12 Christoph Schraff
!  All height obs flagged except lowest active TEMP z-obs (usually surf. level).
!  Surface pressure obs error converted from height to pressure units.
!  Statement functions replaced by elemental or intrinsic functions.
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
! V5_3         2015-10-09 Christoph Schraff
!  More flexible variable names for replication factor and height in wind
!  profilers.
! V5_4         2016-03-10 Christoph Schraff
!  Dimension of 'neventr' and 'neventd' reduced from 3 to 2.
!  Variables related to AOF interface removed.
!
! CAUTION: This module is used commonly by the 3DVAR and COSMO main programs.!!!
!!!        Therefore, anybody wanting to introduce a modification to this    !!!
!!!        module in the context of either of these programs must consult    !!!
!!!        the 'current code owner' of this module for the other program,    !!!
!!!        in order to allow for checking that the modification will comply  !!!
!!!        with both program packages. This must be done before the          !!!
!!!        modification is put into the Version Control System (VCS).        !!!
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

USE netcdf       ! NetCDF F90 interface

USE data_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!-------------------------------------------------------------------------------

USE data_obs_lib_cosmo, ONLY :  &

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
    luse_mlz   ,& ! if false then use multi-level T, not z, and set z-obs to
                  !          passive except for the lowest z-p-obs

! 2. Scalar variables originally defined in other data modules,
!    obtained by calling 'obs_cdf_interface'
! -------------------------------------------------------------

    ! horizontal and vertical sizes of the model fields
    ke             ,& ! number of grid pts in vertical direction (--> 'f_z2p')

    ! variables related to parallelisation / domain decomposition
    num_compute    ,& ! number of compute PEs
    nboundlines    ,& ! number of overlapping boundary lines of the subdomains
    my_cart_id     ,& ! rank of this subdomain in the cartesian communicator
    icomm_cart     ,& ! communicator for the virtual cartesian topology
    imp_reals      ,& ! REAL      type used for MPI
    imp_integers   ,& ! INTEGER   type used for MPI

    ! other variables related to namelist parameters
    maxmlv         ,& ! size (level  dimension) of the  multi-level (m-l)  ODR
    madj_hum          ! = 1 : adjust observed humidity (by ratio of saturation
                      !       vapour pressure over water to the one over ice,
                      !       to be applied if cloud ice is not a state variable

USE data_obs_lib_cosmo, ONLY :  &

    ! constants for the horizontal rotated grid and related variables
    r_pollon       ,& ! longitude of the rotated north pole (in degrees, E>0)
    r_pollat       ,& ! latitude of the rotated north pole (in degrees, N>0) 
    r_polgam       ,& ! angle between the north poles of the systems
    r_dlon         ,& ! grid point distance in zonal direction (in degrees)
    r_dlat         ,& ! grid point distance in meridional direction (in degrees)
    r_startlat_tot ,& ! transformed latitude of the lower left grid point  
                      ! of the total domain (in degrees, N>0)
    r_startlon_tot ,& ! transformed longitude of the lower left grid point 
                      ! of the total domain (in degrees, E>0)
    r_degrad       ,& ! factor for transforming degree to rad

    ! variables related to namelist parameters, and other variables
    doromx         ,& ! SYNOP obs. with height differences betw. model orography
                      !  and station height larger than 'doromx' are set passive
    altopsu        ,& ! SYNOP obs. above height 'altopsu' are set passive

    ! physical constants
    r_g            ,& ! acceleration due to gravity
    tmelt          ,& ! melting temperature of ice
    r_d            ,& ! gas constant for dry air
    b1             ,& ! variables for computing the saturation vapour pressure
    b2w            ,& ! over water (w) and ice (i)
    b2i            ,& !               -- " --
    b3             ,& !               -- " --
    b4w            ,& !               -- " --
    b4i            ,& !               -- " --

    ! switches related to namelist parameters and other
    lverpas        ,& ! write also passive reports on feedback files

! 2b. Pointers for arrays originally defined in other data modules
! ----------------------------------------------------------------

    ! array related to parallelisation / domain decomposition
    i_subpos       ,& ! positions of the subdomains in the total domain
                      ! (i-, j-indices of the lower left + upper right grid pt.
                      !  in the order: i_ll, j_ll, i_ur, j_ur; only the domain
                      !  interior is considered, not the boundary lines)

    ! model fields
    r_p            ,& ! pressure (at main levels)
    r_hhl          ,& ! height   (at half levels)
    r_t_ll         ,& ! temperature at lowest model (main) level
    r_ps              ! surface pressure

USE data_obs_lib_cosmo, ONLY :  &

! 4. I/O device numbers for obs processing / nudging and file names
! -----------------------------------------------------------------
!   nurej      ,& ! direct reporting of rejected obs. reports

! 5. CMA observation type and code type numbers
! ---------------------------------------------

!   nairep     ,& ! AIREP reports (all aircraft reports)
    ntemp      ,& ! TEMP  reports
    npilot     ,& ! PILOT reports
!   nmotcd     ,& !   temp mobile report
!   ntdrop     ,& !   temp drop   report
!   nmopcd     ,& !   pilot mobile report
!   nwp_eu     ,& !   European wind profiler report
!   nra_eu     ,& !   European SODAR/RASS report
!   nravad     ,& !   Radar VAD wind report
!   npr_us     ,& !   US Wind Profiler/RASS report

! 6. Data type with rules for CMA obs and code types
! --------------------------------------------------

!   n_cma      ,& ! number of CMA obs and code types
!   t_cmatyp   ,& ! data type for information on CMA observation and code types
    cma        ,& ! array of meta data on CMA observation and code types

! 7. Functions
! ------------

    i_cma         ! function to determine the index of 'cma'
                  ! referring to a given CMA observation and code type

! end of data_obs_lib_cosmo

!-------------------------------------------------------------------------------

USE data_obs_cdfin, ONLY :  &

! Section 1 : NetCDF Observation Input File formats
!             (this section is used ONLY IF obs data are read from NetCDF files)
!-------------------------------------------------------------------------------

!         1.1   internal attributes of the different NetCDF input files
!               -------------------------------------------------------

    icdfinlen      ,& ! maximum length of NetCDF observation input file name
    iannexlen      ,& ! maximum length of annex of NetCDF obs input file name
    ncdf_temp      ,& ! indicator for processing of NetCDF TEMP       input
    ncdf_tempship  ,& ! indicator for processing of NetCDF TEMPSHIP   input
    ncdf_tempdrop  ,& ! indicator for proc. NetCDF TEMP Dropsonde     input
    ncdf_pilot     ,& ! indicator for proc. NetCDF PILOT (z-levels)   input
    ncdf_pilot_p   ,& ! indicator for proc. NetCDF PILOT (p-levels)   input
    ncdf_amdar_ml  ,& ! indicator for proc. NetCDF AMDAR multi-level  input
    ncdf_amdar_vp  ,& ! indicator for proc. NetCDF AMDAR vert.profile input
    ncdf_amdar     ,& ! indicator for proc. NetCDF AMDAR single-level input
    ncdf_wprof     ,& ! indicator for proc. NetCDF wind profiler      input
    ncdf_rass      ,& ! indicator for proc. NetCDF RASS profiler      input
    ncdf_radar_vad ,& ! indicator for proc. NetCDF radar wind prof.   input
    ycdfin         ,& ! file names of NetCDF observation input files
    icdfin         ,& ! obs file type of NetCDF observation input files
    ncinid         ,& ! unit numbers of NetCDF observation input files
    yncannex       ,& ! annex of NetCDF observation input file names    
    dimids            ! dimension IDs in NetCDF files

USE data_obs_cdfin, ONLY :  &

!         1.2   variables used to read from the data section of the NetCDF files
!               ----------------------------------------------------------------

!         1.2.1 'common' NetCDF header entries and derived variables
!               ----------------------------------------------------

    ilstidn        ,& ! character length of station identity from NetCDF files
! common header entries in NetCDF file
    nc_tisi        ,& ! MTISI  : time significance (BUFR Table 008021)    (dito)
! derived common header variables stored to ODR
    iobstot        ,& ! longitudinal index of grid pt. to which obs is assigned
    jobstot        ,& ! latitudinal  index of grid pt. to which obs is assigned
    iobsloc        ,& ! longitudinal index of grid pt. in local sub-domain
    jobsloc        ,& ! latitudinal  index of grid pt. in local sub-domain
! auxilliary variable, only temporarily available in reader routine
    irproc         ,& ! indices of reports to be processed now

!         1.2.2 other NetCDF header entries
!               ---------------------------

    nc_rstyp       ,& ! NRARA , WMO Common Table C2 : radiosonde type/system
    nc_rad         ,& ! NSR   , BUFR Table B 002013 : solar + IR radiation corr.
    nc_track       ,& ! NSASA , WMO Common Table C7 : tracking technique, status
    nc_na4         ,& ! NA4   , BUFR Table B 002003 : type of measur. equipment
    nc_nix            ! NIX   , BUFR Table B 002001 : station type (man,auto,..)

USE data_obs_cdfin, ONLY :  &

!         1.3   NetCDF body entries
!               -------------------

!         1.3.1  frequent entries
!                ----------------
    nc_nlev        ,& ! MEDRE / MDREP: delayed descriptor replication factor
                      !                (e.g. number of vertical levels)
    nc_dt          ,& ! NLTPD : time [sec] since launch time
    nc_lvtyp       ,& ! MEVSS , BUFR Tab 008042 : extended vertical sounding
                      !                           significance (level identity)
    nc_z           ,& ! NHHHN, NHHH : geopotential height [gpm]    (upper-air)
    nc_dd          ,& ! NDNDN : wind direction      [degree]
    rc_p           ,& ! MPN   , MPPP : pressure
    rc_dlat        ,& ! MLADH : latitude  displacement since launch site
    rc_dlon        ,& ! MLODH : longitude displacement since launch site
    rc_t           ,& ! MTDBT : temperature / dry-bulb temperature
    rc_td          ,& ! MTDNH : dew-point temperature 
    rc_ff          ,& ! NFNFN : wind speed 

!         1.3.3  additional profiler elements
!                ----------------------------

    nc_sinor       ,& ! MSINOR, BUFR Tab 002064 : signal to noise ratio
    nc_qci         ,& ! MQINZ , BUFR Tab 033002 : quality information
    nc_qci2        ,& ! NWPQ  , BUFR Tab 025034 : NOAA QC results   (temporary)
    nc_wce         ,& ! MWCE  , BUFR Tab 025021 : wind computation enhancement
    nc_dto         ,& ! MSETP , time period of measurement [sec]
    nc_dto2        ,& ! NGGTP , time period of measurement [min]    (temporary)
    rc_w           ,& ! MWMPS : vertical velocity (w-component of wind)   [m/s]
    rc_tv          ,& ! MTVIR : virtual temperature                         [K]
    rc_stdff       ,& ! NSTDFF: standard deviation wind speed             [m/s]

!         1.3.4  additional synoptic elements
!                ----------------------------
    nc_clsig       ,& ! MVTSU , BUFR Tab 008002 : vertical significance
    nc_clclm       ,& ! MNH   , BUFR Tab 020011 :(low or mid-level) cloud amount
    nc_ccl         ,& ! MCC   , BUFR Tab 020012 : cloud type (low clouds)
    nc_ccm         ,& ! MCC0  , BUFR Tab 020012 : cloud type (middle clouds)
    nc_cch         ,& ! MCC1  , BUFR Tab 020012 : cloud type (high clouds)
    rc_cbase       ,& ! NH    : height of base of cloud

!         1.4.1  Bit positions for level significance (MEVSS, BUFR Tab 008042)
!                ------------------------------------

    ilv_sfc        ,& ! surface                 level bit
    ilv_std        ,& ! standard                level bit
    ilv_tropo      ,& ! tropopause              level bit
    ilv_max        ,& ! maximum wind            level bit
    ilv_sigt       ,& ! significant temperature level bit
    ilv_sigq       ,& ! significant humidity    level bit
    ilv_sigv       ,& ! significant wind        level bit
    ilv_miss       ,& ! missing value indicator       bit

!         1.5   auxilliary buffer arrays and variables 
!               -------------------------------------- 

    imiss          ,& ! missing value for integers in current NetCDF input file
    rmiss          ,& ! missing value for reals    in current NetCDF input file
    rmisschk          ! value smaller than 'rmiss', to check for missing value

USE data_obs_cdfin, ONLY :  &

! Section 2 : Blacklist and Whitelist
!-------------------------------------------------------------------------------

    ilstid_blk     ,& ! assume 8-character station-IDs in Black-/Whitelist
!   maxintv        ,& ! max. number of vertical blacklist intervals per 1 stat.
    blk_loc           ! blacklists for local reports

USE data_obs_cdfin, ONLY :  &

! Section 3 : Data event counter arrays and diagnostics arrays
!-------------------------------------------------------------------------------

!         3.1     Format of event counters
!                 ------------------------
!         3.1.1   Report event counter array format
!                 ---------------------------------
    nenoda     ,& ! no accepted data in report

!         3.1.2   Data event counter array format
!                 -------------------------------
    mxdeve     ,& ! length of data event counter array
    nelodr     ,& ! level rejected: number of levels exceeding ODR size
    nelmis     ,& ! level rejected: pressure (PILOT: pressure + height) missing
    nelflg     ,& ! level rejected: pressure (PILOT: height) flagged
    nelsfc     ,& ! level rejected: too many surface levels
    nelnop     ,& ! level rejected: PILOT height level not in range of model lev
    nelext     ,& ! level rejected: pressure < 9hPa or level below station alt.
    nelsig     ,& ! level rejected: significant level above a specified limit
    nelrdn     ,& ! level rejected: redundant level in report  (not active yet)
    nepmis     ,& ! pressure (TEMP: height): missing
    nepflg     ,& ! pressure (TEMP: height): flagged
    neprac     ,& ! pressure: bad reporting practice
    nepalt     ,& ! pressure: sta height, or height dist. to orography too large
    nepsdt     ,& ! pressure tendency: flagged, or absolute value too large
    netmis     ,& ! temperature missing
    netflg     ,& ! temperature flagged
    netext     ,& ! temperature too low or too high
    netalt     ,& ! height (diff.) too large for 2m-temp.
    netlps     ,& ! lapse rate of multi-level temperature too large
    neqmis     ,& ! humidity missing
    neqflg     ,& ! humidity flagged
    neqlow     ,& ! humidity too low
    neq300     ,& ! humidity above 300 hpa
    neqbig     ,& ! humidity over allowed value (120%)
    neqsap     ,& ! humidity forced to be saturated (t>0)
    neqsam     ,& ! humidity forced to be saturated (t<0)
    neqclp     ,& ! humidity forced to be <=100% (t>0)
    neqclm     ,& ! humidity forced to be <=100% (t<0)
    neqalt     ,& ! height (diff.) too large for 2m-humid
    nedmis     ,& ! wind direction missing
    nefmis     ,& ! wind speed missing
    nedflg     ,& ! wind direction flagged
    nefflg     ,& ! wind speed flagged
    nefneg     ,& ! wind speed too small  ( < 0 ; DRIBU <= 0 ; VAD < 3m/s )
    nevalt     ,& ! height (diff.) too large for 10m-wind
!   nefshr     ,& ! wind speed shear too large
!   nedshr     ,& ! directional wind shear too large
    nerlim        ! prec.amount exceeds threshold limit

USE data_obs_cdfin, ONLY :  &

!         3.2    Event counter arrays
!                --------------------
    neventr    ,& ! counter array of report events
    neventd    ,& ! counter array of data events

!         3.4    Variables' expectation table
!                ----------------------------
    nzex       ,& ! expect geopotential
    nuex       ,& ! expect horiz. wind
    ntex       ,& ! expect temperature
    ntdex         ! expect humidity

USE data_obs_cdfin, ONLY :  &

!         4.1  Observation error levels
!              ------------------------
    nerlev     ,& ! number of standard error levels
    rolnlv     ,& ! ln(rlevel(15))

!
!         4.2  Observation error constants
!              ---------------------------
    oevsond    ,& ! (root of) radiosonde (TEMP, PILOT) wind error variance
    oezsond    ,& ! (root of) radiosonde (TEMP, PILOT) height error variance
    oetsond    ,& ! (root of) radiosonde temperature error variance
    oeairep    ,& ! (root of) AIREP wind error variance
    oetairp    ,& ! (root of) AIREP temperature error variance
    oevsynp    ,& ! (root of) SYNOP wind error variance
    oezsynp    ,& ! (root of) SYNOP height error variance (land)
    oesatob    ,& ! (root of) SATOB wind error variance
    oevdrib    ,& ! (root of) DRIBU wind   height error variance
    oezdrib    ,& ! (root of) DRIBU height height error variance
    oezship    ,& ! (root of) SHIP (sea SYNOP) height error variance
    rherr1     ,& ! (root of) fixed    / normal conditions
    rherr2     ,& ! relative humidity <  if temperature below 233K
    rherr3        ! error variances    \ if rel. humidity below 20%

USE data_obs_cdfin, ONLY :  &

!         5.2    Temperature / humidity / pressure / height / fog limits
!                -------------------------------------------------------
    rttlim     ,& ! temperature limit below which rel. hum. error is increased
    rrhlim     ,& ! rel. hum. limit below which rel. hum. error is increased
    rhtsat     ,& ! rel. humidity threshold for saturation with real obs
    rtshlm     ,& ! gross error upper limit for relative humidity
!   rerrpp     ,& ! msl pressure above which observed pressure is
!                 ! reckoned to be erroneous
    pminsigt   ,& ! significant-level TEMP / PILOT data are neglected unless
    pminsigv   ,& !   p-obs >= pminsigt (10000 pa) and temperature obs exists or
                  !   p-obs >= pminsigv (20000 pa) and wind obs exists
    pqmin      ,& ! pressure [pa] of level above which moisture obs are not used
    rpplim     ,& ! pressure level obove which observations are not used
!   vfoglim    ,& ! visibility threshold [m] below which the existence of low
!                 !      cloud (fog) is assumed in the presence of precipitation
    fflim_vad  ,& ! lower limit for accepting VAD wind speed

!         7.0    For data rejection messages: Output buffer, size and formats
!                ------------------------------------------------------------
    outbuf     ,& ! buffer containing output for a single node
    nacout     ,& ! actual number of records stored in the output buffer
    nmxoln     ,& ! maximum length of output buffer
    istrej     ,& ! length of strings (station id) in output buffer
!   nfmt3      ,& ! no accepted data
    nfmt4      ,& ! excess of levels
    nfmt5         ! several surface levels

! end of data_obs_cdfin

!-------------------------------------------------------------------------------

USE data_obs_record, ONLY :   &

!       1.     Header formats of ODR reports
!       ------------------------------------

!       1.1.1  Header formats of ODR reports: 'omlhed' and 'osghed'
!              ----------------------------------------------------
!   mxrhed     ,& ! header length of multi-level reports
!   mxshed     ,& ! header length of single-level reports
    nhilon     ,& ! longitude of observing station
    nhjlat     ,& ! latitude  of observing station
    nhalt      ,& ! station altitude [m]
    nhtime     ,& ! (exact) time of observation in forecast hours
    nhsurf     ,& ! height of model grid pt. to which obs. is assigned
!   nhzio      ,& ! longitude of obs. station (or lowest datum) in grid pt. unit
!   nhzjo      ,& ! latitude  of obs. station in grid pt. units
!   nhsynt     ,& ! nominal (synoptic) time of observation in forecast hours
!   nhvcbu     ,& ! correction factor to vertical correlation scale for wind
                  ! at base of report
!   nhvcbt     ,& ! as 'nhvcbu', but for temperature
!   nhvcbq     ,& ! as 'nhvcbu', but for humidity
!   nhvctu     ,& ! correction factor to vertical correlation scale for wind
                  ! at top of report
!   nhvctt     ,& ! as 'nhvctu', but for temperature
!   nhvctq     ,& ! as 'nhvctu', but for humidity

!       1.1.2  Header formats of ODR reports: 'momlhd' and 'mosghd'
!              ----------------------------------------------------
!   mxrhdf     ,& ! header length of multi-level reports
!   mxshdf     ,& ! header length of single-level reports
    nhio       ,& ! (local) x-coord. of grid pt. assigned to obs
    nhjo       ,& ! (local) y-coord. of grid pt. assigned to obs
!   nhitot     ,& ! global x-coord. of grid pt. assigned to obs
!   nhjtot     ,& ! global y-coord. of grid pt. assigned to obs
    nhobtp     ,& ! observation type
    nhcode     ,& ! code type
    nhschr     ,& ! station characteristics                      (see 1.1.4)
    nhpass     ,& ! flag for report being set to 'passive'       (see 1.1.4)
    nhqcfw     ,& ! status of QC and of writing to feedobs files (and p-QC flag)
    nhflag     ,& ! report flags (obs type, surface, altitude, station ID)
!   nhcorr     ,& ! update sequence number (station correction indicator)
!   nhcat      ,& ! data     category (from BUFR Section 1)
!   nhcats     ,& ! data sub-category (from BUFR Section 1)
    nhkz       ,& ! DWD internal classification number (observation type)
!   nhcent     ,& ! originating centre
!   nhstid     ,& ! station identity number
!   nhdate     ,& ! absolute exact observation date [yyyymmdd]
!   nhhrmn     ,& ! absolute exact observation time [hhmm]
!   nhsyhr     ,& ! absolute nominal (synoptic) observation time [yymmddhh]
    nhnlev     ,& ! number of obs. levels (for multi-level reports)
!   nhvqcf     ,& ! for satellite retrieval: threshold quality control flags
    nhaexi     ,& ! flag for exist. of wind or temperature in multi-level report
    nhuexi     ,& ! flag for existence of wind data        in multi-level report
    nhtexi     ,& ! flag for existence of temperature data in multi-level report
    nhqexi     ,& ! flag for existence of humidity data    in multi-level report
    nhrtyp     ,& ! radiosonde type    (NRARA, see WMO common code Table C2)
    nhtrac     ,& ! tracking technique (NSASA, see WMO common code Table C7)
    nhrad      ,& ! solar and IR radiation correction (NSR, BUFR Table 002013)
    nhna4      ,& ! instrument type                   (NA4, BUFR Table 002003)
    nhwce      ,& ! wind comput. enhancement (w-prof, MWCE, BUFR Table 025021)
    nhdt       ,& ! time period of measurement (e.g. w-prof)               [s]
!   nhstyp     ,& ! surface obs: station type (buoy: MQOBL, BUFR Table 002149,
                  !                            else: NIX  , BUFR Table 002001)

!       1.1.3  Header formats of ODR reports: 'yomlhd' and 'yosghd'
!              ----------------------------------------------------

    ilstid     ,& ! character length of the station identity
!   ilstidp    ,& ! char. length used for printing the station ID
                  ! Note: (ilstid >= ilstidg >= ilstidp), cf. data_nudge_gather

!       1.2    Bit patterns for packed information in ODR (and VOF) header
!              -----------------------------------------------------------
    nvsebp     ,& ! bit pos. for report located at sea grid pt.     nhschr
    nvapbp     ,& ! bit pos. for phase of flight (aircraft)           "
    nvapoc     ,& ! no. of bits occ. by phase of flight               "
    nvaabp     ,& ! bit pos. for aircraft roll angle (code)           "
    nvaaoc        ! no. of bits occ. by aircraft roll angle           "

USE data_obs_record, ONLY :   &

!       1.3    ODR body format
!              ---------------

!       1.3.1  Body format of ODR of multi-level reports: 'omlbdy'
!              ---------------------------------------------------
    mxrbdy     ,& ! body length of multi-level reports
    nbtu       ,& ! u wind component [m/s]
    nbtv       ,& ! v wind component [m/s]
    nbtt       ,& ! temperature [K]
    nbtrh      ,& ! relative humidity [/]
    nbtp       ,& ! pressure [Pa]
    nbtz       ,& ! height [m]
    nbtuer     ,& ! error of observed wind component
    nbtter     ,& ! error of observed temperature
    nbtqer     ,& ! error of observed rel. humidity
    nbtzer     ,& ! error of observed height
                  !  Note: nbt?er are set to the negative rms errors, if the
                  !  observations have not passed the threshold quality control
    nbtzio     ,& ! longitude in grid pt. units
    nbtzjo     ,& ! latitude  in grid pt. units
    nbttim     ,& ! observation time relative to report (header) time
    nbtlop     ,& ! LOG( pressure )
    nbtdrh     ,& ! bias correction for relative humidity [/]
    nbtw       ,& ! vertical velocity [m/s]
    nbtsnr     ,& ! signal to noise ratio 
    nbtuac     ,& ! accuracy (std dev from data provider) of horiz. wind [m/s]


!       1.3.2  Body format of ODR of multi-level report flags: 'momlbd'
!              --------------------------------------------------------
    mxrbdf     ,& ! body length of multi-level reports
    nbtflg     ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbterr     ,& ! status flag word        (bit pattern, see below: 'nb?err')
    nbtqcf     ,& ! threshold quality control flags      (see below: 'nb?qcf')
    nbtlsg     ,& ! level id (bit pattern, as in NetCDF statistics file)
    nbtlid        ! level identity          (bit pattern, see below: 'nb?lid')

USE data_obs_record, ONLY :   &

!       1.3.3  Body format of ODR of surface reports: 'osgbdy'
!              -----------------------------------------------
    mxsbdy     ,& ! body length of single-level reports
    nbsu       ,& ! u wind component                                   [m/s]
    nbsv       ,& ! v wind component                                   [m/s]
    nbst       ,& ! temperature                                        [K]
    nbsrh      ,& ! relative humidity                                  [/]
    nbsp       ,& ! pressure                                           [Pa]
    nbsz       ,& ! height                                             [m]
    nbsuer     ,& ! error of observed wind component
    nbster     ,& ! error of observed temperature
    nbsqer     ,& ! error of observed relative humidity
    nbszer     ,& ! error of observed height
    nbspst     ,& ! (3-hourly) pressure tendency                       [Pa/3h]
    nbscbs     ,& ! (lowest) cloud base height                         [m]
    nbscl      ,& ! low       cloud cover        (BUFR Table 020011)   [octas]
    nbscm      ,& ! mid-level cloud cover        (BUFR Table 020011)   [octas]
    nbsch      ,& ! high      cloud cover        (BUFR Table 020011)   [octas]
    nbsct      ,& ! total     cloud cover        (BUFR Table 020011)   [octas]
!   nbsvis     ,& ! (horizontal) visibility                            [m]
!   nbsrr1     ,& ! precipitation amount over 1 hour                   [mm]
!   nbsrr6     ,& ! precipitation amount over 6 hours                  [mm]
!   nbsr12     ,& ! precipitation amount over 12 hours                 [mm]
!   nbsr24     ,& ! precipitation amount over 24 hours                 [mm]
!   nbsfgv     ,& ! max. derived equivalent vertical gust (aircraft)   [m/s]
!   nbsfg1     ,& ! max. wind speed of gusts over 1 hour               [m/s]
!   nbsfg6     ,& ! max. wind speed of gusts over 6 hours              [m/s]
!   nbstn      ,& ! minimum temperature (at 2m during past 12 hrs)     [K]
!   nbstx      ,& ! maximum temperature (at 2m during past 12 hrs)     [K]
!   nbshsw     ,& ! total snow depth                                   [m]
    nbsdrh     ,& ! bias correction for relative humidity [/]

!       1.3.4  Body format of ODR of surface report flags: 'mosgbd'
!              ----------------------------------------------------
    mxsbdf     ,& ! body length of single-level reports
    nbsflg     ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbserr     ,& ! status flag word        (bit pattern, see below: 'nb?err')
    nbsqcf     ,& ! threshold quality control flags      (see below: 'nb?qcf')
    nbslid     ,& ! SYNOP: pressure code (SYNOP)   (code, see below: 'nbslid')
                  ! else : level identity   (bit pattern, see below: 'nb?lid')
    nbscwg     ,& ! combined cloud and weather group (set of classes, s below)
    nbswwe     ,& ! NetCDF read, SYNOP: weather and ground group word  (below)
    nbstur     ,& ! NetCDF read, Aircraft: degree of turbulence WMO Tab 011031
                  !   (not contained in merged multi-level aircraft reports !)
    nbsclg        ! general           cloud       group (code)
!   nbscl1     ,& ! first  individual cloud layer group (code)
!   nbscl2     ,& ! second individual cloud layer group (code)
!   nbscl3     ,& ! third  individual cloud layer group (code)
!   nbscl4        ! forth  individual cloud layer group (code)

USE data_obs_record, ONLY :   &

!       1.4    Bit patterns for packed information in ODR (and VOF) body
!       ------------------------------------------------------------------------
!       1.4.2  Other bit patt. for packed info in ODR (VOF) body, general words
!              ----------------------------------------------------------------
    nvru       ,& ! bit pos. for status/QC flags for horiz. wind   nb?err/nb?qcf
    nvrt       ,& ! bit pos. for status/QC flags for temperature         "
    nvrq       ,& ! bit pos. for status/QC flags for humidity            "
    nvrz       ,& ! bit pos. for status/QC flags for pressure/height     "
    nvrw       ,& ! bit pos. for status/QC flags for vertical wind       "
    nvfubp     ,& ! bit pos. for main flag on wind                     nb?flg
    nvftbp     ,& ! bit pos. for main flag on temperature                "
    nvfqbp     ,& ! bit pos. for main flag on humidity                   "
    nvfzbp     ,& ! bit pos. for main flag on pressure / geopot.         "
    nvfaoc     ,& ! no. of bits occ. by each main flag                   "
    nvfbps     ,& ! bit pattern for main flags:                          "
    nvfboc     ,& ! no. of bits occ. for main flags                      "
    nvflbp     ,& ! bit pos. for level flag: level below surface         "
!   nvfloc     ,& ! no. of bits occ. by level flag                       "
    nvlidp        ! level id. bit pattern                              nb?lid
!   nvlido        ! no. bits occ. by each indicator in level id.         "

USE data_obs_record, ONLY :   & 

!       1.4.3  Bit patterns for 'optional groups' in ODR body 'mosgbd' (and VOF)
!              -----------------------------------------------------------------
                  ! combined cloud and weather group                   nbscwg
    nvchbp     ,& ! bit position for ch  (type of high cloud)
    nvcmbp     ,& !         "        cm  (type of middle cloud)
    nvclbp     ,& !         "        cl  (type of low cloud)
    nvnhbp     ,& !         "        nh  (cover of low, else of middle cloud)
    nvhbp      ,& !         "        h   (cloud base height)
    nvnbp      ,& !         "        n   (total cloud cover)
!   nvwwbp     ,& !         "        ww  (present weather)
                  !                      (see VUB WMO Code tables:)
    nvchoc     ,& ! no. of bits occupied by ch    [Code table 0509]
    nvcmoc     ,& !           "             cm    [Code table 0515]
    nvcloc     ,& !           "             cl    [Code table 0513]
    nvnhoc     ,& !           "             nh    [Code table 2700]
    nvhoc      ,& !           "             h     [Code table 1600]
    nvnoc      ,& !           "             n     [Code table 2700]
!   nvwwoc     ,& !           "             ww    [Code table 4677]
                  ! --> weather and ground group word                  nbswwe
    nvcqbp     ,& ! bit position for refined quality flag on ccl
    nvcqoc     ,& ! no. bits occupied by refined quality flag on ccl
                  !
                  ! --> general    cloud       group word              nbsclg
    nctlbp     ,& ! bit position for low    cloud type code (WMO Table 020012)
    nctmbp     ,& ! bit position for middle cloud type code (WMO Table 020012)
    ncthbp     ,& ! bit position for high   cloud type code (WMO Table 020012)
                  !
                  ! --> individual cloud layer group words             nbscl?
    nxsgbp     ,& ! bit position for vertic. signific. code (WMO Table 008002)
    nxclbp     ,& ! bit position for cloud amount      code (WMO Table 020011)
    nxsgoc     ,& ! no. bits occupied for vert. signf. code (WMO Table 008002)
    nxcloc     ,& ! no. bits occupied for cloud amount code (WMO Table 020011)
    nxctoc        ! no. bits occupied for cloud type   code (WMO Table 020012)

USE data_obs_record, ONLY :   &

!       1.5    Further quantities related to ODR
!              ---------------------------------
    imdi       ,& ! missing data indicator for ODR integers (2^31-1)
    ntotml     ,& ! tot. number of stored multi-level reports
    ntotsg     ,& ! tot. number of stored single-level reports
    fdoro      ,& ! scaling factor to vertical distances betw. model

!       2.     Observation data records (ODR)
!       -------------------------------------

    omlbdy     ,& ! body   of multi-level ODR
    omlhed     ,& ! header of multi-level ODR
    osgbdy     ,& ! body   of single-level ODR
    osghed     ,& ! header of single-level ODR
    momlbd     ,& ! body   of multi-level ODR
    momlhd     ,& ! header of multi-level ODR
    mosgbd     ,& ! body   of single-level ODR
    mosghd     ,& ! header of single-level ODR
    yomlhd     ,& ! header of multi-level ODR
    yosghd     ,& ! header of single-level ODR

!       3.     Masking constants
!       ------------------------

    nibits        ! masking constants

! end of data_obs_record

!-------------------------------------------------------------------------------

  USE mo_fdbk_tables,          ONLY :  &

    FL_NO_OBS     ,& ! no (active) observations in report
    LS_SURFACE    ,& ! surface
    LS_STANDARD   ,& ! standard level
    LS_TROPO      ,& ! tropopause level
    LS_MAX        ,& ! maximum wind level
    LS_SIGN          ! significant level

! end of mo_fdbk_tables

!-------------------------------------------------------------------------------

 USE environment,              ONLY :  &
    model_abort        ! aborts the program in case of errors

!-------------------------------------------------------------------------------

 USE parallel_utilities,       ONLY :  &
!   global_values,   & ! computes global values by operating on local arrays
    distribute_values  ! distributes a set of values from one node to all others

!-------------------------------------------------------------------------------

 USE utilities,                ONLY :  &
    phi2phirot,      & ! converts phi from the real to the rotated system
    rla2rlarot,      & ! converts lambda from the real to the rotated system
    uv2uvrot_vec       ! converts wind components from normal to rotated system

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_util,       ONLY :  &
    obs_assign_sort_node   ,& ! assign node to reports and sort them accordingly
    obs_cdf_distrib_reports,& ! distribute reports to appropriate sub-domains
    std_atmosphere         ,& ! convert variables accord. to standard atmosphere
    obs_td2rh              ,& ! convert dewpoint temperature to relat. humidity
!   obs_qx2rh              ,& ! convert mixing ratio to model-compatible rel hum
!   obs_rhw2rh             ,& ! make (observed) relat. humidity model compatible
    obs_find_level         ,& ! get interpolation levels and factors
    f_z2p                     ! get (model) pressure for a specified height

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_comhead,       ONLY :  &
    obs_cdf_read_comhead   ,& ! reads and evaluates common header information,
                              !   including grid point assignment of reports
    obs_cdf_buffer_comhead ,& ! puts header info in buffer arrays for distribut.
    obs_cdf_store_comhead     ! stores common header info in ODR

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_blk,           ONLY :  &
    obs_cdf_whitelist_local,& ! indic. for local reports if missing on whitelist
    obs_cdf_blacklist_local   ! produces a list of blacklisted vertical
                              !   intervals for each of a set of local reports

!-------------------------------------------------------------------------------
! Local declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE

!===============================================================================

!-------------------------------------------------------------------------------
! Public and Private Subroutines
!-------------------------------------------------------------------------------

CONTAINS


!===============================================================================
!+ Module procedure in "src_obs_cdfin_mult" for reading and distributing reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_read_temp_pilot ( min_sta , min_end , ilcdf                 &
                                   , nodrnew , nexceed)

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_mult" organizes the reading,
!   pre-selection, distribution and storage in ODR of radiosonde (TEMP and
!   PILOT) reports from a given observation time period.
!
! Method:
!   The reports are read by 1 node (processor) and assigned each to its most
!   appropriate grid point of the total model domain. The reports then have
!   to be distributed to those sub-domains which contain the grid point they
!   are assigned to. To do this, the reports are temporarily stored in three
!   long buffer arrays (one each for integer, real, and character elements)
!   in the order according to these sub-domains, and then the appropriate
!   sections of these arrays are distributed to the nodes.
!   Pre-processing steps related to the observation body elements and storing
!   the data in long-term arrays (ODR) are then performed locally, i.e. in
!   parallel mode.
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER (KIND=iintegers) , INTENT (IN)  ::  &
    min_sta       ,& ! start  \  of time interval to be read now
    min_end       ,& ! end    /  (in [minutes] of forecast time)
    ilcdf            ! index (number) of NetCDF observation input file

  INTEGER (KIND=iintegers) , INTENT (INOUT)  ::  &
    nodrnew (2)   ,& ! number of new reports
    nexceed (2)      ! number of reports in excess of array size

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    nodrepn (num_compute) ,& ! number of            reports  \   to be
    nodleni (num_compute) ,& ! number of integer   elements   \  distributed
    nodlenr (num_compute) ,& ! number of real      elements   /  to the
    nodleny (num_compute) ,& ! number of character elements  /   different nodes
    nrepl                 ,& ! number of            reports  \  ( to be
    nlenli                ,& ! number of integer   elements   \   received )
    nlenlr                ,& ! number of real      elements   /  at the
    nlenly                ,& ! number of character elements  /   local node
    ibuflen  (3)          ,& ! length of (i/r/c) buffer arrays for distrib. rep.
    iloclen  (3)             ! length of (i/r/c) buffer arrays received locally

  INTEGER (KIND=iintegers) ::  &
    kcdftyp        ,& ! type of NetCDF input file format (observation type):
                      !  = ncdf_temp      : TEMP
                      !  = ncdf_tempship  : TEMPSHIP
                      !  = ncdf_pilot     : PILOT  (height levels)
                      !  = ncdf_pilot_p   : PILOT  (pressure levels)
    ncid           ,& ! file unit number of current NetCDF obs input file
    mrepsta        ,& ! index (in record dimension) of 1st report to be read
    nrep           ,& ! index interval (in record dimension) which contains
                      !   those reports that should be read now
    nreproc        ,& ! number of reports to be read and processed now
    noffscal (3,2) ,& ! number int/real/char - elements in common header part
                      !                      - scalar elements in total report
    iioffs         ,& ! offset of current report in long integer buffer
    iroffs         ,& ! offset of current report in long real    buffer
    niscal         ,& ! number of scalar integer   elements
    nrscal         ,& ! number of scalar real      elements
    nyscal         ,& ! number of scalar character elements
    nrealml        ,& ! number of multi-level real elements
    inode          ,& ! node to which the current report is assigned
    ipe            ,& ! loop index over nodes
    irps           ,& ! loop index over reports to be processed
    irep           ,& ! record index over reports to be processed
    ilev           ,& ! loop index over vertical levels
    nlev           ,& ! number of vertical levels
    maxlev         ,& ! max number of vertical levels in all processd reports
    mxlev          ,& ! length of vertical level dimension in NetCDF file
    ndims          ,& ! number of dimensions (for 'edition_number' in NetCDF)
    dimid_mxlv     ,& ! dimension ID for vertical levels in NetCDF file
    status         ,& ! NetCDF status variable
    nsta2 (2)      ,& ! start indices       \  of the reports to be read
    ncnt2 (2)      ,& ! number of elements  /  in the external NetCDF
    ktimsig        ,& ! significance of report time
    istat  , ierr     ! error indicators

  INTEGER (KIND=iintegers) ::  &          ! variable ID's in NetCDF file for:
    varid_nlev, varid_lvtyp, varid_z   ,& ! number of levels, level type, height
    varid_dlat, varid_dlon , varid_dt  ,& ! latitude / longitude / time shifts
    varid_p   , varid_t    , varid_td  ,& ! pressure, temperature, dew point
    varid_dd  , varid_ff               ,& ! wind direction / speed
                varid_clsig            ,& ! vertical significance for cloud
    varid_mnh , varid_nh               ,& ! low cloud amount, cloud base height
    varid_mccl, varid_mccm , varid_mcch,& ! type of low / mid-level / high cloud
    varid_na4 , varid_rstyp            ,& ! measurem. equipment, radiosonde type
    varid_rad , varid_track            ,& ! radiation correction, tracking tech.
                varid_mtisi               ! time significance

  CHARACTER (LEN=25)       :: &
    yroutine          ! name of this subroutine
  CHARACTER (LEN=150)      :: &
    yerrmsl           ! error message
  CHARACTER (LEN=80)       :: &
    ymsg              ! control message
  CHARACTER (LEN=80)       :: &
    ytxtob            ! observation type in text messages
  CHARACTER (LEN= icdfinlen + iannexlen)  ::  &
    yfn               ! file name of NetCDF observation input file, with annex

! Local arrays:
! ------------

  INTEGER (KIND=iintegers), ALLOCATABLE :: &
    irnode     (:) ,& ! nodes to which the reports will be distributed
    irsort     (:) ,& ! report indices sorted according to 'irnode'
    iirlen     (:) ,& ! number of integer   elements  \   in each
    irrlen     (:) ,& ! number of real      elements   >  individual
    iyrlen     (:)    ! number of character elements  /   report

  INTEGER (KIND=iintegers), ALLOCATABLE :: &
    ibufsrt    (:) ,& ! integer buffer array with sorted reports read from file
    ibufloc    (:)    ! integer buffer array with local reports (at sub-domain)

  REAL    (KIND=wp)       , ALLOCATABLE :: &
    rbufsrt    (:) ,& ! real    buffer array with sorted reports read from file
    rbufloc    (:)    ! real    buffer array with local reports (at sub-domain)

  CHARACTER (LEN=ilstidn) , ALLOCATABLE :: & 
    ybufsrt    (:) ,& ! character array with sorted reports read from file
    ybufloc    (:)    ! character array with local reports (at sub-domain)
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_read_temp_pilot
!-------------------------------------------------------------------------------

  yroutine = 'obs_cdf_read_temp_pilot'

  ncid     =  ncinid (ilcdf)
  kcdftyp  =  icdfin (ilcdf)
  yfn      =  ycdfin (kcdftyp) (1:LEN_TRIM( ycdfin (kcdftyp) )) //             &
              yncannex (ilcdf) (1:LEN_TRIM( yncannex (ilcdf) ))

! PRINT *, yroutine, kcdftyp, min_sta, min_end

  IF (kcdftyp == ncdf_temp)      ytxtob = 'TEMP '
  IF (kcdftyp == ncdf_tempship)  ytxtob = 'TEMPSHIP'
  IF (kcdftyp == ncdf_tempdrop)  ytxtob = 'TEMPDROP'
  IF (kcdftyp == ncdf_pilot)     ytxtob = 'PILOT'
  IF (kcdftyp == ncdf_pilot_p)   ytxtob = 'PILOT (p-levels)'

  DO ipe = 1 , num_compute
    nodrepn (ipe) = 0
    nodleni (ipe) = 0
    nodlenr (ipe) = 0
    nodleny (ipe) = 0
  ENDDO

! define number of int/real/char elements from common header part 
! to be put into 'Xbufsrt'  (equal for all NetCDF observation input files)
! ------------------------------------------------------------------------
  noffscal (1,1)  =  3 + 18           ! 18 integer  elements + 3 length counters
  noffscal (2,1)  =  9                !  9 real     elements
  noffscal (3,1)  =  1                !  1 character element

! define total number of scalar int/real/char elements 
! to be put into 'Xbufsrt'  (specific for TEMPs here)
! ----------------------------------------------------
  IF ((kcdftyp == ncdf_temp ) .OR. (kcdftyp == ncdf_tempship)) THEN
    noffscal (1,2)  =  noffscal(1,1) + 11
    noffscal (2,2)  =  noffscal(2,1) + 1
    noffscal (3,2)  =  noffscal(3,1) + 0
  ELSEIF (kcdftyp == ncdf_tempdrop) THEN
    noffscal (1,2)  =  noffscal(1,1) + 6
    noffscal (2,2)  =  noffscal(2,1) + 0
    noffscal (3,2)  =  noffscal(3,1) + 0
  ELSEIF ((kcdftyp == ncdf_pilot) .OR. (kcdftyp == ncdf_pilot_p)) THEN
    noffscal (1,2)  =  noffscal(1,1) + 5
    noffscal (2,2)  =  noffscal(2,1) + 0
    noffscal (3,2)  =  noffscal(3,1) + 0
  ENDIF

!-------------------------------------------------------------------------------
! Section 1: Read and pre-process header information common to all obs types:
!            time / lat / lon / station altitude / obs type / station ID /
!            center of origin ...;
!            assign observations to grid points and to nodes (sub-domains)
!-------------------------------------------------------------------------------

  nrep = 0

  IF (my_cart_id == 0) THEN

    CALL obs_cdf_read_comhead ( min_sta , min_end , ilcdf                      &
                              , mrepsta , nrep , nreproc )
!   =========================
  ENDIF

! IF (num_compute > 1)  CALL global_values ( nrep, 1,'MAX',imp_integers        &
!                                                , icomm_cart, -1, yerr,ierr )
!                       ------------------
  IF (num_compute > 1) THEN
    CALL distribute_values (nrep    ,1,0,imp_integers,icomm_cart,ierr)
    CALL distribute_values (imiss   ,1,0,imp_integers,icomm_cart,ierr)
    CALL distribute_values (rmisschk,1,0,imp_reals   ,icomm_cart,ierr)
  ENDIF

  IF ((my_cart_id == 0) .AND. (nrep >= 1)) THEN

! assign the observations to the nodes (sub-domains), determine
! the local grid point indices, and produce a list of indices 
! pointing to the reports sorted according to the nodes
! -------------------------------------------------------------

    ALLOCATE ( irsort (nreproc+1) , STAT=istat )
    ALLOCATE ( irnode (nreproc+1) , STAT=istat )

    CALL obs_assign_sort_node ( nrep, nreproc, irproc, iobstot, jobstot        &
                              , num_compute, i_subpos, nboundlines, my_cart_id &
                              , irnode, irsort, iobsloc, jobsloc )
!   =========================

! 'irnode' has been allocated in 'obs_cdf_read_comhead' and is not used any more
    DEALLOCATE ( irproc , STAT=istat )

!-------------------------------------------------------------------------------
! Section 2: Get entries specific to observation type (TEMP)
!-------------------------------------------------------------------------------

    ALLOCATE ( nc_rstyp (nrep) , STAT=istat )
    ALLOCATE ( nc_track (nrep) , STAT=istat )
    ALLOCATE ( nc_na4   (nrep) , STAT=istat )
    ALLOCATE ( nc_tisi  (nrep) , STAT=istat )
    ALLOCATE ( nc_nlev  (nrep) , STAT=istat )
    DO irep = 1 , nrep
      nc_rstyp(irep)  =  imiss
      nc_track(irep)  =  imiss
      nc_na4  (irep)  =  imiss
      nc_tisi (irep)  =  imiss
!     nc_nlev (irep)  =  imiss
    ENDDO
    IF ((kcdftyp == ncdf_temp ) .OR. (kcdftyp == ncdf_tempship)) THEN
      ALLOCATE ( nc_rad   (nrep) , STAT=istat )
      ALLOCATE ( nc_clsig (nrep) , STAT=istat )
      ALLOCATE ( nc_clclm (nrep) , STAT=istat )
      ALLOCATE ( rc_cbase (nrep) , STAT=istat )
      ALLOCATE ( nc_ccl   (nrep) , STAT=istat )
      ALLOCATE ( nc_ccm   (nrep) , STAT=istat )
      ALLOCATE ( nc_cch   (nrep) , STAT=istat )
      DO irep = 1 , nrep
        nc_rad   (irep)  =  imiss
        nc_clsig (irep)  =  imiss
        nc_clclm (irep)  =  imiss
        rc_cbase (irep)  =  rmiss
        nc_ccl   (irep)  =  imiss
        nc_ccm   (irep)  =  imiss
        nc_cch   (irep)  =  imiss
      ENDDO
    ELSEIF (kcdftyp == ncdf_tempdrop) THEN
      ALLOCATE ( nc_rad   (nrep) , STAT=istat )
    ENDIF

! get header info on type of radiosonde / tracking / radiation correction ...
! ---------------------------------------------------------------------------
! Time significance: 2: time averaged; 3: accumulated; 18: radiosonde launch t.;
!   (MTISI, 008021)  23: monitoring period; 25: nominal reporting time;
!                    26: time of last known position; 31: missing value
!                    (for other values, see WMO BUFR code Table 008021)
! Radiosonde type  : 17: Graw G. (D); 26: Meteolabor Basora (CH);
!   (NRARA)          61: Vaisala RS80 Loran, Digicora I,II or Marwin;
!                    71: Vaisala RS90 Digicora I,II or Marwin; 34: Vinohrady(CZ)
!                    79: Vaisala RS90 Digicora I,II or Marwin;
!                    80: Vaisala RS92 Digicora III; 81: Vaisala RS92 Autosonde
!                    255: missing value
!                    (for other values, see WMO common code Table C2)
! Tracking technique:0: no windfinding; 3: automatic with auxilliary ranging;
!   (NSASA)          2: automatic with auxilliary radio detection finding;
! (Common Table C7)  6: automatic cross chain Loran-C; 8: automatic satellite
!                       navigation; 19: tracking technique not specified;
!                    127: missing value;
!                    (for other values, see WMO common code Table C7)
! Type of measuring: 0: pressure instr. associat. with wind measuring equipment;
!    equipment used: 1: optical theodolite; 2: radio theodolite; 3: radar;
!   (NA4, 002003)    4: VLF-Omega; 5: Loran-C; 6: wind profiler; 
!                    7: satellite navigation; 8: RASS; 9: Sodar; 15: missing
! Solar & infrared : 0: no correction; 1: CIMO solar + CIMO infrared corrected;
! radiation correct: 2: CIMO solar + infrared corr.; 3: CIMO solar corr. only;
!   (NSR, 002013)    4: solar + infrared corr. automatically by radiosonde syst;
!                    5: solar corr. autom. by radios.; 6: solar + infrared corr.
!                       as specified by country; 7: solar corr. by country;
!                    15: missing

! treat radiosonde type as mandatory variables
                              status = nf90_inq_varid (ncid,'NRARA',varid_rstyp)
    IF (status /= nf90_noerr) THEN
      yerrmsl = 'RADIOSONDE TYPE NRARA DOES NOT EXIST IN ' // yfn
      CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
    ENDIF
    status = nf90_get_var (ncid, varid_rstyp, nc_rstyp, (/mrepsta/), (/nrep/))

! treat other variables as optional
    status = nf90_inq_varid (ncid,'NSASA',varid_track)
    IF (status == nf90_noerr)                                                  &
      status = nf90_get_var (ncid, varid_track, nc_track, (/mrepsta/), (/nrep/))
    status = nf90_inq_varid (ncid,'NA4'  ,varid_na4  )
    IF (status == nf90_noerr)                                                  &
      status = nf90_get_var (ncid, varid_na4  , nc_na4  , (/mrepsta/), (/nrep/))
    status = nf90_inq_varid (ncid,'MTISI',varid_mtisi)
    IF (status == nf90_noerr)                                                  &
      status = nf90_get_var (ncid, varid_mtisi, nc_tisi , (/mrepsta/), (/nrep/))
    IF ((kcdftyp == ncdf_temp ) .OR. (kcdftyp == ncdf_tempship)                &
                                .OR. (kcdftyp == ncdf_tempdrop)) THEN
      status = nf90_inq_varid (ncid,'NSR'  ,varid_rad  )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_rad, nc_rad  , (/mrepsta/), (/nrep/))
    ENDIF

! control message if some / all TEMPs report nominal time instead of launch time
    ktimsig = 0
    DO irps = 1 , nreproc
      irep  =  irsort(irps)
      IF (nc_tisi(irep) == 18) ktimsig = 2* (ktimsig /2) + 1
      IF (nc_tisi(irep) == 25) ktimsig = MOD(ktimsig, 2) + 2
      IF (ktimsig == 2) THEN
        WRITE( ymsg,'("NOTE: all ",A ,"s report NOMINAL TIME only (forecast"   &
                    &," minutes:",I5," - ",I5,")")' )                          &
               ytxtob(1:LEN_TRIM(ytxtob)), min_sta, min_end
      ELSEIF (ktimsig == 3) THEN
        WRITE( ymsg,'("NOTE: some ",A ,"s report NOMINAL TIME only (forecast"  &
                    &," minutes:",I5," - ",I5,")")' )                          &
               ytxtob(1:LEN_TRIM(ytxtob)), min_sta, min_end
      ENDIF
      IF (ktimsig >= 2)  PRINT         '(A)' , ymsg
!     IF (ktimsig >= 2)  WRITE( nurej ,'(A)' ) ymsg
    ENDDO

! get cloud information from TEMP: treat as optional
! -------------------------------

    IF ((kcdftyp == ncdf_temp ) .OR. (kcdftyp == ncdf_tempship)) THEN
      status = nf90_inq_varid (ncid,'MVTSU',varid_clsig)
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_clsig, nc_clsig, (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MNH'  ,varid_mnh  )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_mnh  , nc_clclm, (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'NH'   ,varid_nh   )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_nh   , rc_cbase, (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MCC'  ,varid_mccl )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_mccl , nc_ccl  , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MCC0' ,varid_mccm )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_mccm , nc_ccm  , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MCC1' ,varid_mcch )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_mcch , nc_cch  , (/mrepsta/),(/nrep/))
    ENDIF

! get number of vertical levels
! -----------------------------

                              status = nf90_inq_varid (ncid,'MEDRE' ,varid_nlev)
    IF (status /= nf90_noerr) THEN
      yerrmsl = 'NUMBER OF VERTICAL LEVELS DOES NOT EXIST IN ' // yfn
      CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
    ENDIF
    status = nf90_get_var (ncid, varid_nlev, nc_nlev, (/mrepsta/), (/nrep/))
! get maximum number of vertical levels (according to entry 'MEDRE')
    maxlev = 0
    DO irps = 1 , nreproc
      irep  =  irsort(irps)
      maxlev = MAX( nc_nlev(irep) , maxlev )
    ENDDO

! check maximum number of vertical levels
! ---------------------------------------

    IF (maxlev >= 1) THEN
                               status =nf90_inq_varid (ncid,'MEVSS',varid_lvtyp)
      IF (status /= nf90_noerr) THEN
        yerrmsl = 'VERT. SOUND. SIGNIF. MEVSS DOES NOT EXIST IN ' // yfn
        CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
      ENDIF
! get length of vertical dimension in NetCDF file
      status = nf90_Inquire_Variable ( ncid, varid_lvtyp, ndims=ndims          &
                                     , dimids=dimids)
      dimid_mxlv = dimids(1)
      status = nf90_Inquire_Dimension (ncid, dimid_mxlv, len=mxlev)
      IF ((maxlev > mxlev) .OR. (ndims /= 2)) THEN
        PRINT *, 'ERROR with dimensions in file ', ndims, mxlev, maxlev
        yerrmsl = 'NUMBER OF LEVELS EXCEEDS DIMENSION IN ' // yfn
        CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
      ENDIF
    ENDIF

! get profile observations
! ------------------------

!   PRINT *,'varid3 ', ncid, maxlev, ncdf_temp, kcdftyp, nrep
    IF (maxlev >= 1) THEN
      ALLOCATE ( nc_lvtyp (maxlev,nrep) , STAT=istat )
      ALLOCATE ( nc_dt    (maxlev,nrep) , STAT=istat )
      ALLOCATE ( nc_z     (maxlev,nrep) , STAT=istat )
      ALLOCATE ( nc_dd    (maxlev,nrep) , STAT=istat )
      ALLOCATE ( rc_dlat  (maxlev,nrep) , STAT=istat )
      ALLOCATE ( rc_dlon  (maxlev,nrep) , STAT=istat )
      ALLOCATE ( rc_ff    (maxlev,nrep) , STAT=istat )
      ALLOCATE ( rc_p     (maxlev,nrep) , STAT=istat )
!     nc_z    (:,:) = imiss            !  need not exist for (ncdf_pilot_p)
      IF ((kcdftyp == ncdf_temp ) .OR. (kcdftyp == ncdf_tempship)              &
                                  .OR. (kcdftyp == ncdf_tempdrop)) THEN
        ALLOCATE ( rc_t     (maxlev,nrep) , STAT=istat )
        ALLOCATE ( rc_td    (maxlev,nrep) , STAT=istat )
      ENDIF

! 1. variables treated as mandatory
                               status =nf90_inq_varid (ncid,'MEVSS',varid_lvtyp)
      IF(status == nf90_noerr) status =nf90_inq_varid (ncid,'NDNDN',varid_dd   )
      IF(status == nf90_noerr) status =nf90_inq_varid (ncid,'NFNFN',varid_ff   )
      IF ((kcdftyp == ncdf_temp ) .OR. (kcdftyp == ncdf_tempship)              &
                                  .OR. (kcdftyp == ncdf_tempdrop)) THEN
        IF(status == nf90_noerr) status =nf90_inq_varid (ncid,'NHHHN',varid_z  )
        IF(status == nf90_noerr) status =nf90_inq_varid (ncid,'MPN'  ,varid_p  )
        IF(status == nf90_noerr) status =nf90_inq_varid (ncid,'MTDBT',varid_t  )
        IF(status == nf90_noerr) status =nf90_inq_varid (ncid,'MTDNH',varid_td )
      ELSEIF (kcdftyp == ncdf_pilot  ) THEN
        IF(status == nf90_noerr) status =nf90_inq_varid (ncid,'NHHH' ,varid_z  )
      ELSEIF (kcdftyp == ncdf_pilot_p) THEN
        IF(status == nf90_noerr) status =nf90_inq_varid (ncid,'MPN'  ,varid_p  )
      ENDIF
      IF (status /= nf90_noerr) THEN
        yerrmsl = 'STANDARD '// ytxtob // ' PROFILE DATA DO NOT EXIST IN '// yfn
        CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
      ENDIF

      nsta2 (1) = 1
      ncnt2 (1) = maxlev
      nsta2 (2) = mrepsta
      ncnt2 (2) = nrep
      status = nf90_get_var (ncid, varid_lvtyp, nc_lvtyp, nsta2, ncnt2)
      status = nf90_get_var (ncid, varid_dd   , nc_dd   , nsta2, ncnt2)
      status = nf90_get_var (ncid, varid_ff   , rc_ff   , nsta2, ncnt2)
      IF ((kcdftyp == ncdf_temp ) .OR. (kcdftyp == ncdf_tempship)              &
                                  .OR. (kcdftyp == ncdf_tempdrop)) THEN
        status = nf90_noerr
        status = nf90_get_var (ncid, varid_z  , nc_z    , nsta2, ncnt2)
        status = nf90_get_var (ncid, varid_p  , rc_p    , nsta2, ncnt2)
        status = nf90_get_var (ncid, varid_t  , rc_t    , nsta2, ncnt2)
        status = nf90_get_var (ncid, varid_td , rc_td   , nsta2, ncnt2)
!       IF (status /= nf90_noerr) PRINT *,'ppm3 ', TRIM(nf90_strerror(status))
      ELSEIF (kcdftyp == ncdf_pilot  ) THEN
        rc_p  (:,:)  =  rmiss
        status = nf90_get_var (ncid, varid_z  , nc_z    , nsta2, ncnt2)
      ELSEIF (kcdftyp == ncdf_pilot_p) THEN
        nc_z  (:,:)  =  imiss
        status = nf90_get_var (ncid, varid_p  , rc_p    , nsta2, ncnt2)
      ENDIF
!     PRINT *,'press0 ', varid_p, nsta2, ncnt2, maxlev, nrep
!     PRINT '("press1 ",10F8.0)' , (rc_p(ilev,1),ilev=1,5)                     &
!                                , (rc_p(1,irep),irep=1,5) 

! 1. variables treated as optional
      nc_dt   (:,:) = imiss
      rc_dlat (:,:) = rmiss
      rc_dlon (:,:) = rmiss
      status = nf90_inq_varid (ncid,'NLTPD',varid_dt   )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_dt   , nc_dt   , nsta2, ncnt2)
      status = nf90_inq_varid (ncid,'MLADH',varid_dlat )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_dlat , rc_dlat , nsta2, ncnt2)
      status = nf90_inq_varid (ncid,'MLODH',varid_dlon )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_dlon , rc_dlon , nsta2, ncnt2)
      IF     (kcdftyp == ncdf_pilot  ) THEN
        status = nf90_inq_varid (ncid,'MPN'  ,varid_p  )
        IF (status == nf90_noerr)                                              &
          status = nf90_get_var (ncid, varid_p   , rc_p   , nsta2, ncnt2)
      ELSEIF (kcdftyp == ncdf_pilot_p) THEN
        status = nf90_inq_varid (ncid,'NHHH' ,varid_z  )
        IF (status == nf90_noerr)                                              &
          status = nf90_get_var (ncid, varid_z   , nc_z   , nsta2, ncnt2)
      ENDIF
    ENDIF

!-------------------------------------------------------------------------------
! Section 3: Produce long buffer arrays 'Xbufsrt' in which the read reports are
!            sorted according to the nodes to which they will be distributed
!-------------------------------------------------------------------------------

    ALLOCATE ( iirlen (nreproc+1) , STAT=istat )
    ALLOCATE ( irrlen (nreproc+1) , STAT=istat )
    ALLOCATE ( iyrlen (nreproc+1) , STAT=istat )

! compute length (number of elements) for each individual TEMP report
! -------------------------------------------------------------------

! total number of scalar elements
    niscal = noffscal(1,2)
    nrscal = noffscal(2,2)
    nyscal = noffscal(3,2)

    IF ((kcdftyp == ncdf_temp ) .OR. (kcdftyp == ncdf_tempship)                &
                                .OR. (kcdftyp == ncdf_tempdrop))  nrealml = 6
    IF ((kcdftyp == ncdf_pilot) .OR. (kcdftyp == ncdf_pilot_p ))  nrealml = 4
    DO irps = 1 , nreproc
      irep  =  irsort(irps)
! determine length of integer / real buffer for each report
      iirlen (irps)  =  niscal +      4 *nc_nlev(irep)
      irrlen (irps)  =  nrscal + nrealml*nc_nlev(irep)
      iyrlen (irps)  =  nyscal
    ENDDO

! compute length of buffer arrays (i.e length of all reports together + 1)
! -------------------------------

    ibuflen (1)  =  1
    ibuflen (2)  =  1
    ibuflen (3)  =  1
    DO irps = 1 , nreproc
      ibuflen (1)  =  ibuflen(1) + iirlen (irps)
      ibuflen (2)  =  ibuflen(2) + irrlen (irps)
      ibuflen (3)  =  ibuflen(3) + iyrlen (irps)
    ENDDO

! allocate buffer arrays (only for my_cart_id == 0 here)
! ----------------------

    ALLOCATE ( ibufsrt (ibuflen(1)) , STAT=istat )
    ALLOCATE ( rbufsrt (ibuflen(2)) , STAT=istat )
    ALLOCATE ( ybufsrt (ibuflen(3)) , STAT=istat )

! fill the header part which is common to all obs types 
! into the long buffer arrays 'Xbufsrt'
! -----------------------------------------------------

    CALL obs_cdf_buffer_comhead ( nreproc, irsort, iirlen, irrlen, iyrlen      &
                                , ibuflen, ibufsrt, rbufsrt, ybufsrt )
!   ===========================

! fill the remaining scalar elements into the long buffer arrays 'Xbufsrt'
! ------------------------------------------------------------------------

    iioffs = 0
    iroffs = 0

    DO irps = 1 , nreproc
! the following part is common to all observation types
      irep  =  irsort(irps)
      ibufsrt (iioffs+noffscal(1,1)+ 1)  =  nc_nlev (irep)
      ibufsrt (iioffs+noffscal(1,1)+ 2)  =  nc_rstyp(irep)
      ibufsrt (iioffs+noffscal(1,1)+ 3)  =  nc_track(irep)
      ibufsrt (iioffs+noffscal(1,1)+ 4)  =  nc_na4  (irep)
      ibufsrt (iioffs+noffscal(1,1)+ 5)  =  nc_tisi (irep)
      IF ((kcdftyp == ncdf_temp ) .OR. (kcdftyp == ncdf_tempship)) THEN
        ibufsrt (iioffs+noffscal(1,1)+ 6)  =  nc_rad  (irep)
        ibufsrt (iioffs+noffscal(1,1)+ 7)  =  nc_clsig(irep)
        ibufsrt (iioffs+noffscal(1,1)+ 8)  =  nc_clclm(irep)
        ibufsrt (iioffs+noffscal(1,1)+ 9)  =  nc_ccl  (irep)
        ibufsrt (iioffs+noffscal(1,1)+10)  =  nc_ccm  (irep)
        ibufsrt (iioffs+noffscal(1,1)+11)  =  nc_cch  (irep)
        rbufsrt (iroffs+noffscal(2,1)+ 1)  =  rc_cbase(irep)
      ELSEIF (kcdftyp == ncdf_tempdrop) THEN
        ibufsrt (iioffs+noffscal(1,1)+ 6)  =  nc_rad  (irep)
      ENDIF

      iioffs  =  iioffs + iirlen(irps)
      iroffs  =  iroffs + irrlen(irps)
    ENDDO

    DEALLOCATE ( nc_rstyp , STAT=istat )
    DEALLOCATE ( nc_track , STAT=istat )
    DEALLOCATE ( nc_na4   , STAT=istat )
    DEALLOCATE ( nc_tisi  , STAT=istat )
    IF ((kcdftyp == ncdf_temp ) .OR. (kcdftyp == ncdf_tempship)) THEN
      DEALLOCATE ( nc_rad   , STAT=istat )
      DEALLOCATE ( nc_clsig , STAT=istat )
      DEALLOCATE ( nc_clclm , STAT=istat )
      DEALLOCATE ( rc_cbase , STAT=istat )
      DEALLOCATE ( nc_ccl   , STAT=istat )
      DEALLOCATE ( nc_ccm   , STAT=istat )
      DEALLOCATE ( nc_cch   , STAT=istat )
    ELSEIF (kcdftyp == ncdf_tempdrop) THEN
      DEALLOCATE ( nc_rad   , STAT=istat )
    ENDIF

! fill in the multi-level info into the long buffer arrays 'Xbufsrt'
! ------------------------------------------------------------------

    IF (maxlev >= 1) THEN
      iioffs = 0
      iroffs = 0
      DO irps = 1 , nreproc
        irep  =  irsort(irps)
        nlev = nc_nlev(irep)
        DO ilev = 1 , nlev
          ibufsrt (iioffs+niscal       +ilev)  =  nc_lvtyp(ilev,irep)
          ibufsrt (iioffs+niscal+  nlev+ilev)  =  nc_dd   (ilev,irep)
          ibufsrt (iioffs+niscal+2*nlev+ilev)  =  nc_z    (ilev,irep)
          ibufsrt (iioffs+niscal+3*nlev+ilev)  =  nc_dt   (ilev,irep)
          rbufsrt (iroffs+nrscal       +ilev)  =  rc_dlat (ilev,irep)
          rbufsrt (iroffs+nrscal+  nlev+ilev)  =  rc_dlon (ilev,irep)
          rbufsrt (iroffs+nrscal+2*nlev+ilev)  =  rc_ff   (ilev,irep)
          rbufsrt (iroffs+nrscal+3*nlev+ilev)  =  rc_p    (ilev,irep)
          IF ((kcdftyp == ncdf_temp ) .OR. (kcdftyp == ncdf_tempship)          &
                                      .OR. (kcdftyp == ncdf_tempdrop)) THEN
            rbufsrt (iroffs+nrscal+4*nlev+ilev)  =  rc_t    (ilev,irep)
            rbufsrt (iroffs+nrscal+5*nlev+ilev)  =  rc_td   (ilev,irep)
          ENDIF
        ENDDO
!       PRINT '("bufa0 ",4I7   )' , irps, irep, nlev, iroffs
!       PRINT '("bufaa ",8F10.2)' , (rbufsrt(iroffs+istat),istat=1,8)
!       PRINT '("bufab ",8F10.2)' , (rbufsrt(iroffs+istat),istat=9,16)
!       PRINT '("bufac ",8F10.2)' , (rbufsrt(iroffs+nrscal+3*nlev+istat),istat=9,16)
        iioffs  =  iioffs + iirlen(irps)
        iroffs  =  iroffs + irrlen(irps)
      ENDDO

      DEALLOCATE ( nc_lvtyp , STAT=istat )
      DEALLOCATE ( nc_dt    , STAT=istat )
      DEALLOCATE ( nc_z     , STAT=istat )
      DEALLOCATE ( nc_dd    , STAT=istat )
      DEALLOCATE ( rc_dlat  , STAT=istat )
      DEALLOCATE ( rc_dlon  , STAT=istat )
      DEALLOCATE ( rc_ff    , STAT=istat )
      DEALLOCATE ( rc_p     , STAT=istat )
      IF ((kcdftyp == ncdf_temp ) .OR. (kcdftyp == ncdf_tempship)              &
                                  .OR. (kcdftyp == ncdf_tempdrop)) THEN
        DEALLOCATE ( rc_t     , STAT=istat )
        DEALLOCATE ( rc_td    , STAT=istat )
      ENDIF
    ENDIF
    DEALLOCATE ( nc_nlev  , STAT=istat )
!   PRINT '("buf0a ",8F10.2)' , (rbufsrt(istat),istat=1,8)
!   PRINT '("buf0b ",8F10.2)' , (rbufsrt(istat),istat=9,16)

!-------------------------------------------------------------------------------
! Section 4: Distribute the reports
!-------------------------------------------------------------------------------

! determine the number of int / real / char elements 
! to be distributed to the different nodes

    DO irps = 1 , nreproc
      inode =  irnode(irps)
      nodrepn (inode+1) = nodrepn(inode+1) + 1
      nodleni (inode+1) = nodleni(inode+1) + iirlen(irps)
      nodlenr (inode+1) = nodlenr(inode+1) + irrlen(irps)
      nodleny (inode+1) = nodleny(inode+1) + iyrlen(irps)
!     PRINT *,'distr2 ', irps, irnode(irps), iirlen(irps), irrlen(irps), iyrlen(irps)
    ENDDO
!   PRINT *,'distr3a ', nodrepn
!   PRINT *,'distr3b ', nodleni
!   PRINT *,'distr3c ', nodlenr
!   PRINT *,'distr3d ', nodleny

    DEALLOCATE ( iirlen  , STAT=istat )
    DEALLOCATE ( irrlen  , STAT=istat )
    DEALLOCATE ( iyrlen  , STAT=istat )
    DEALLOCATE ( irsort  , STAT=istat )
    DEALLOCATE ( irnode  , STAT=istat )

  ENDIF    ! (my_cart_id == 0)

  IF ((nrep >= 1) .AND. (num_compute > 1)) THEN

    iloclen(1)  =  MAXVAL( nodleni ) + 1
    iloclen(2)  =  MAXVAL( nodlenr ) + 1
    iloclen(3)  =  MAXVAL( nodleny ) + 1

    CALL distribute_values (ibuflen, 3, 0,imp_integers,icomm_cart,ierr)
    CALL distribute_values (iloclen, 3, 0,imp_integers,icomm_cart,ierr)

    IF (my_cart_id /= 0)  ALLOCATE ( ibufsrt (ibuflen(1)) , STAT=istat )
    IF (my_cart_id /= 0)  ALLOCATE ( rbufsrt (ibuflen(2)) , STAT=istat )
    IF (my_cart_id /= 0)  ALLOCATE ( ybufsrt (ibuflen(3)) , STAT=istat )
                          ALLOCATE ( ibufloc (iloclen(1)) , STAT=istat )
                          ALLOCATE ( rbufloc (iloclen(2)) , STAT=istat )
                          ALLOCATE ( ybufloc (iloclen(3)) , STAT=istat )

    CALL obs_cdf_distrib_reports ( ibuflen, ibufsrt, rbufsrt, ybufsrt, ilstidn &
                                 , nodrepn, nodleni, nodlenr, nodleny          &
                                 , nrepl  , nlenli , nlenlr , nlenly           &
                                 , iloclen, ibufloc, rbufloc, ybufloc          &
                                 , num_compute , my_cart_id , icomm_cart       &
                                 , imp_integers, imp_reals )
!   ============================

    DEALLOCATE ( ibufsrt , STAT=istat )
    DEALLOCATE ( rbufsrt , STAT=istat )
    DEALLOCATE ( ybufsrt , STAT=istat )

!-------------------------------------------------------------------------------
! Section 5: Store the local reports in the ODR
!-------------------------------------------------------------------------------

!   PRINT '("buf3a ",8F10.2)' , (rbufloc(istat),istat=1,8)
!   PRINT '("buf3b ",8F10.2)' , (rbufloc(istat),istat=9,16)

    CALL obs_cdf_store_multilev ( kcdftyp , nrepl  , nlenli , nlenlr , nlenly  &
                                , noffscal,          ibufloc, rbufloc, ybufloc &
                                , nodrnew , nexceed )
!   ===========================

    DEALLOCATE ( ibufloc , STAT=istat )
    DEALLOCATE ( rbufloc , STAT=istat )
    DEALLOCATE ( ybufloc , STAT=istat )

  ELSEIF (nrep >= 1) THEN
    nrepl   =  nodrepn(1)
    nlenli  =  nodleni(1)
    nlenlr  =  nodlenr(1)
    nlenly  =  nodleny(1)

    CALL obs_cdf_store_multilev ( kcdftyp , nrepl  , nlenli , nlenlr , nlenly  &
                                , noffscal,          ibufsrt, rbufsrt, ybufsrt &
                                , nodrnew , nexceed )
!   ===========================

    DEALLOCATE ( ibufsrt , STAT=istat )
    DEALLOCATE ( rbufsrt , STAT=istat )
    DEALLOCATE ( ybufsrt , STAT=istat )

  ENDIF

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_read_temp_pilot
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_read_temp_pilot


!===============================================================================
!+ Module procedure in "src_obs_cdfin_mult" for reading and distributing reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_read_profiler ( min_sta , min_end , ilcdf                   &
                                 , nodrnew , nexceed)

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_mult" organizes the reading,
!   pre-selection, distribution and storage in ODR of remote-sensing profiler
!   (wind profiler, RASS radio acoustic sounding system virtual temperature
!   profiler, radar VAD wind profiler) reports from a given observation time
!   period.
!
! Method:
!   The reports are read by 1 node (processor) and assigned each to its most
!   appropriate grid point of the total model domain. The reports then have
!   to be distributed to those sub-domains which contain the grid point they
!   are assigned to. To do this, the reports are temporarily stored in three
!   long buffer arrays (one each for integer, real, and character elements)
!   in the order according to these sub-domains, and then the appropriate
!   sections of these arrays are distributed to the nodes.
!   Pre-processing steps related to the observation body elements and storing
!   the data in long-term arrays (ODR) are then performed locally, i.e. in
!   parallel mode.
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER (KIND=iintegers) , INTENT (IN)  ::  &
    min_sta       ,& ! start  \  of time interval to be read now
    min_end       ,& ! end    /  (in [minutes] of forecast time)
    ilcdf            ! index (number) of NetCDF observation input file

  INTEGER (KIND=iintegers) , INTENT (INOUT)  ::  &
    nodrnew (2)   ,& ! number of new reports
    nexceed (2)      ! number of reports in excess of array size

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    nodrepn (num_compute) ,& ! number of            reports  \   to be
    nodleni (num_compute) ,& ! number of integer   elements   \  distributed
    nodlenr (num_compute) ,& ! number of real      elements   /  to the
    nodleny (num_compute) ,& ! number of character elements  /   different nodes
    nrepl                 ,& ! number of            reports  \  ( to be
    nlenli                ,& ! number of integer   elements   \   received )
    nlenlr                ,& ! number of real      elements   /  at the
    nlenly                ,& ! number of character elements  /   local node
    ibuflen  (3)          ,& ! length of (i/r/c) buffer arrays for distrib. rep.
    iloclen  (3)             ! length of (i/r/c) buffer arrays received locally

  INTEGER (KIND=iintegers) ::  &
    kcdftyp        ,& ! type of NetCDF input file format (observation type):
                      !  = ncdf_wprof     : Wind Profiler
                      !  = ncdf_radar_vad : Radar VAD
                      !  = ncdf_rass      : RASS
    ncid           ,& ! file unit number of current NetCDF obs input file
    mrepsta        ,& ! index (in record dimension) of 1st report to be read
    nrep           ,& ! index interval (in record dimension) which contains
                      !   those reports that should be read now
    nreproc        ,& ! number of reports to be read and processed now
    noffscal (3,2) ,& ! number int/real/char - elements in common header part
                      !                      - scalar elements in total report
    iioffs         ,& ! offset of current report in long integer buffer
    iroffs         ,& ! offset of current report in long real    buffer
    niscal         ,& ! number of scalar integer   elements
    nrscal         ,& ! number of scalar real      elements
    nyscal         ,& ! number of scalar character elements
    nrealml        ,& ! number of multi-level real elements
    nintml         ,& ! number of multi-level int  elements
    inode          ,& ! node to which the current report is assigned
    ipe            ,& ! loop index over nodes
    irps           ,& ! loop index over reports to be processed
    irep           ,& ! record index over reports to be processed
    ilev           ,& ! loop index over vertical levels
    nlev           ,& ! number of vertical levels
    maxlev         ,& ! max number of vertical levels in all processd reports
    mxlev          ,& ! length of vertical level dimension in NetCDF file
    ndims          ,& ! number of dimensions (for 'edition_number' in NetCDF)
    dimid_mxlv     ,& ! dimension ID for vertical levels in NetCDF file
    status         ,& ! NetCDF status variable
    nsta2 (2)      ,& ! start indices       \  of the reports to be read
    ncnt2 (2)      ,& ! number of elements  /  in the external NetCDF
    istat  , ierr     ! error indicators

  INTEGER (KIND=iintegers) ::  & ! variable ID's in NetCDF file for:
    varid_nlev , varid_z      ,& ! number of levels, height (MH)
    varid_dd   , varid_ff     ,& ! wind direction / speed
    varid_w    , varid_tv     ,& ! vertical velocity / virtual temperature
    varid_sinor               ,& ! signal to noise ratio (NSINOR)
    varid_stdff               ,& ! standard deviation wind speed (NSTDFF)
    varid_qci  , varid_qci2   ,& ! quality indices (MQINZ / NWPQ)
    varid_na4                 ,& ! measurement. equipment
    varid_wce                 ,& ! wind computation enhancement (MWCE)
    varid_mtisi               ,& ! time significance
    varid_dto  , varid_dto2      ! time period of obs. [sec] / [min]

  CHARACTER (LEN=25)       :: &
    yroutine          ! name of this subroutine
  CHARACTER (LEN=70)       :: &
    yerrmsl           ! error message
  CHARACTER (LEN= icdfinlen + iannexlen)  ::  &
    yfn               ! file name of NetCDF observation input file, with annex

! Local arrays:
! ------------

  INTEGER (KIND=iintegers), ALLOCATABLE :: &
    irnode     (:) ,& ! nodes to which the reports will be distributed
    irsort     (:) ,& ! report indices sorted according to 'irnode'
    iirlen     (:) ,& ! number of integer   elements  \   in each
    irrlen     (:) ,& ! number of real      elements   >  individual
    iyrlen     (:)    ! number of character elements  /   report

  INTEGER (KIND=iintegers), ALLOCATABLE :: &
    ibufsrt    (:) ,& ! integer buffer array with sorted reports read from file
    ibufloc    (:)    ! integer buffer array with local reports (at sub-domain)

  REAL    (KIND=wp)       , ALLOCATABLE :: &
    rbufsrt    (:) ,& ! real    buffer array with sorted reports read from file
    rbufloc    (:)    ! real    buffer array with local reports (at sub-domain)

  CHARACTER (LEN=ilstidn) , ALLOCATABLE :: & 
    ybufsrt    (:) ,& ! character array with sorted reports read from file
    ybufloc    (:)    ! character array with local reports (at sub-domain)
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_read_profiler
!-------------------------------------------------------------------------------

  yroutine = 'obs_cdf_read_profiler'

  ncid     =  ncinid (ilcdf)
  kcdftyp  =  icdfin (ilcdf)
  yfn      =  ycdfin (kcdftyp) (1:LEN_TRIM( ycdfin (kcdftyp) )) //             &
              yncannex (ilcdf) (1:LEN_TRIM( yncannex (ilcdf) ))

! PRINT *, yroutine, kcdftyp, min_sta, min_end

  DO ipe = 1 , num_compute
    nodrepn (ipe) = 0
    nodleni (ipe) = 0
    nodlenr (ipe) = 0
    nodleny (ipe) = 0
  ENDDO

! define number of int/real/char elements from common header part 
! to be put into 'Xbufsrt'  (equal for all NetCDF observation input files)
! ------------------------------------------------------------------------
  noffscal (1,1)  =  3 + 18           ! 18 integer  elements + 3 length counters
  noffscal (2,1)  =  9                !  9 real     elements
  noffscal (3,1)  =  1                !  1 character element

! define total number of scalar int/real/char elements 
! to be put into 'Xbufsrt'  (specific for Profilers here)
! -------------------------------------------------------
  noffscal (1,2)  =  noffscal(1,1) + 5
  noffscal (2,2)  =  noffscal(2,1) + 0
  noffscal (3,2)  =  noffscal(3,1) + 0

!-------------------------------------------------------------------------------
! Section 1: Read and pre-process header information common to all obs types:
!            time / lat / lon / station altitude / obs type / station ID /
!            center of origin ...;
!            assign observations to grid points and to nodes (sub-domains)
!-------------------------------------------------------------------------------

  nrep = 0

  IF (my_cart_id == 0) THEN

    CALL obs_cdf_read_comhead ( min_sta , min_end , ilcdf                      &
                              , mrepsta , nrep , nreproc )
!   =========================
  ENDIF

! IF (num_compute > 1)  CALL global_values ( nrep, 1,'MAX',imp_integers        &
!                                                , icomm_cart, -1, yerr,ierr )
!                       ------------------
  IF (num_compute > 1) THEN
    CALL distribute_values (nrep    ,1,0,imp_integers,icomm_cart,ierr)
    CALL distribute_values (imiss   ,1,0,imp_integers,icomm_cart,ierr)
    CALL distribute_values (rmisschk,1,0,imp_reals   ,icomm_cart,ierr)
  ENDIF

  IF ((my_cart_id == 0) .AND. (nrep >= 1)) THEN

! assign the observations to the nodes (sub-domains), determine
! the local grid point indices, and produce a list of indices 
! pointing to the reports sorted according to the nodes
! -------------------------------------------------------------

    ALLOCATE ( irsort (nreproc+1) , STAT=istat )
    ALLOCATE ( irnode (nreproc+1) , STAT=istat )

    CALL obs_assign_sort_node ( nrep, nreproc, irproc, iobstot, jobstot        &
                              , num_compute, i_subpos, nboundlines, my_cart_id &
                              , irnode, irsort, iobsloc, jobsloc )
!   =========================

! 'irnode' has been allocated in 'obs_cdf_read_comhead' and is not used any more
    DEALLOCATE ( irproc , STAT=istat )

!-------------------------------------------------------------------------------
! Section 2: Get header entries specific to remote-sensing profilers
!-------------------------------------------------------------------------------

    ALLOCATE ( nc_na4   (nrep) , STAT=istat )
    ALLOCATE ( nc_tisi  (nrep) , STAT=istat )
    ALLOCATE ( nc_dto   (nrep) , STAT=istat )
    ALLOCATE ( nc_dto2  (nrep) , STAT=istat )
    ALLOCATE ( nc_wce   (nrep) , STAT=istat )
    DO irep = 1 , nrep
      nc_na4  (irep)  =  imiss
      nc_tisi (irep)  =  imiss
      nc_dto  (irep)  =  imiss
      nc_dto2 (irep)  =  imiss
      nc_wce  (irep)  =  imiss
    ENDDO

! get header info on type of radiosonde / tracking / radiation correction ...
! ---------------------------------------------------------------------------
! Time significance: 2: time averaged; 31: missing value
!   (MTISI, 008021)  (for other values, see WMO BUFR code Table 008021)
! Type of measuring: 0: pressure instr. associat. with wind measuring equipment;
!    equipment used: 1: optical theodolite; 2: radio theodolite; 3: radar;
!   (NA4, 002003)    4: VLF-Omega; 5: Loran-C; 6: wind profiler; 
!                    7: satellite navigation; 8: RASS; 9: Sodar; 15: missing

    status = nf90_inq_varid (ncid,'NA4'  ,varid_na4  )
    IF (status == nf90_noerr)                                                  &
      status = nf90_get_var (ncid, varid_na4  , nc_na4  , (/mrepsta/), (/nrep/))
    status = nf90_inq_varid (ncid,'MTISI',varid_mtisi)
    IF (status == nf90_noerr)                                                  &
      status = nf90_get_var (ncid, varid_mtisi, nc_tisi , (/mrepsta/), (/nrep/))
    status = nf90_inq_varid (ncid,'MSETP',varid_dto  )
    IF (status == nf90_noerr)                                                  &
      status = nf90_get_var (ncid, varid_dto  , nc_dto  , (/mrepsta/), (/nrep/))
    status = nf90_inq_varid (ncid,'NGGTP',varid_dto2 )
    IF (status == nf90_noerr)                                                  &
      status = nf90_get_var (ncid, varid_dto2 , nc_dto2 , (/mrepsta/), (/nrep/))
    status = nf90_inq_varid (ncid,'MWCE' ,varid_wce  )
    IF (status == nf90_noerr)                                                  &
      status = nf90_get_var (ncid, varid_wce  , nc_wce  , (/mrepsta/), (/nrep/))

! if time significance is not 'time averaged' then unset 'time period of obs'
    DO irps = 1 , nreproc
      irep  =  irsort(irps)
      IF ((nc_dto(irep) == imiss) .AND. (nc_dto2(irep) /= imiss))              &
        nc_dto (irep) = nc_dto2(irep) * 60
      IF (      (nc_tisi(irep) /= imiss) .AND. (nc_tisi(irep) /= 31)           &
          .AND. (nc_tisi(irep) /= 2))   nc_dto (irep) = imiss
    ENDDO

! get number of vertical levels ('MDREP' or 'MEDRE' is assumed to exist always)
! -----------------------------

    ALLOCATE ( nc_nlev  (nrep) , STAT=istat )
                              status = nf90_inq_varid (ncid,'MDREP' ,varid_nlev)
    IF (status /= nf90_noerr) status = nf90_inq_varid (ncid,'MEDRE' ,varid_nlev)
    IF (status /= nf90_noerr) THEN
      yerrmsl = 'NUMBER OF VERTICAL LEVELS DOES NOT EXIST IN ' // yfn
      CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
    ENDIF
    status = nf90_get_var (ncid, varid_nlev, nc_nlev, (/mrepsta/), (/nrep/))
! get maximum number of vertical levels (according to entry 'MDREP')
    maxlev = 0
    DO irps = 1 , nreproc
      irep  =  irsort(irps)
      maxlev = MAX( nc_nlev(irep) , maxlev )
    ENDDO

! check maximum number of vertical levels ('MH' or 'MHVEH' is assumed to exist)
! ---------------------------------------

    IF (maxlev >= 1) THEN
                                status = nf90_inq_varid (ncid,'MH'    ,varid_z)
      IF (status /= nf90_noerr) status = nf90_inq_varid (ncid,'MHVEH' ,varid_z)
      IF (status /= nf90_noerr) THEN
        yerrmsl = 'LEVEL HEIGHT DOES NOT EXIST IN ' // yfn
        CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
      ENDIF
! get length of vertical dimension in NetCDF file
      status = nf90_Inquire_Variable ( ncid, varid_z    , ndims=ndims          &
                                     , dimids=dimids)
      dimid_mxlv = dimids(1)
      status = nf90_Inquire_Dimension (ncid, dimid_mxlv, len=mxlev)
      IF (maxlev > mxlev) THEN
        yerrmsl = 'NUMBER OF LEVELS EXCEEDS DIMENSION IN ' // yfn
        CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
      ENDIF
    ENDIF

! get profile observations
! ------------------------

!   PRINT *,'varid3 ', ncid, maxlev, kcdftyp, nrep
    IF (maxlev >= 1) THEN
      ALLOCATE ( nc_z     (maxlev,nrep) , STAT=istat )
      ALLOCATE ( rc_w     (maxlev,nrep) , STAT=istat )
      ALLOCATE ( nc_sinor (maxlev,nrep) , STAT=istat )
      ALLOCATE ( nc_qci   (maxlev,nrep) , STAT=istat )
      IF ((kcdftyp == ncdf_wprof) .OR. (kcdftyp == ncdf_radar_vad)) THEN
        ALLOCATE ( nc_dd    (maxlev,nrep) , STAT=istat )
        ALLOCATE ( rc_ff    (maxlev,nrep) , STAT=istat )
        ALLOCATE ( rc_stdff (maxlev,nrep) , STAT=istat )
      ELSEIF (kcdftyp == ncdf_rass) THEN
        ALLOCATE ( rc_tv    (maxlev,nrep) , STAT=istat )
      ENDIF

      DO irps = 1 , nreproc
        irep  =  irsort(irps)
        DO ilev = 1 , maxlev
!         nc_z    (ilev,irep) = imiss
          rc_w    (ilev,irep) = rmiss
          nc_qci  (ilev,irep) = imiss
          nc_sinor(ilev,irep) = imiss
          IF ((kcdftyp == ncdf_wprof) .OR. (kcdftyp == ncdf_radar_vad)) THEN
            nc_dd   (ilev,irep) = imiss
            rc_ff   (ilev,irep) = rmiss
            rc_stdff(ilev,irep) = rmiss
          ELSEIF (kcdftyp == ncdf_rass) THEN
            rc_tv   (ilev,irep) = imiss  ! should that be rmiss, because rc_tv is real???
          ENDIF
        ENDDO
      ENDDO
      nsta2 (1) = 1
      ncnt2 (1) = maxlev
      nsta2 (2) = mrepsta
      ncnt2 (2) = nrep
      status = nf90_get_var   (ncid, varid_z    , nc_z    , nsta2, ncnt2)
      status = nf90_inq_varid (ncid,'MWMPS' ,varid_w    )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_w    , rc_w    , nsta2, ncnt2)
      status = nf90_inq_varid (ncid,'MQINZ' ,varid_qci  )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_qci  , nc_qci  , nsta2, ncnt2)
      status = nf90_inq_varid (ncid,'NSINOR',varid_sinor)
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_sinor, nc_sinor, nsta2, ncnt2)
      IF ((kcdftyp == ncdf_wprof) .OR. (kcdftyp == ncdf_radar_vad)) THEN
        status = nf90_inq_varid (ncid,'NDNDN' ,varid_dd   )
        IF (status == nf90_noerr)                                              &
          status = nf90_get_var (ncid, varid_dd   , nc_dd   , nsta2, ncnt2)
        status = nf90_inq_varid (ncid,'NFNFN' ,varid_ff   )
        IF (status == nf90_noerr)                                              &
          status = nf90_get_var (ncid, varid_ff   , rc_ff   , nsta2, ncnt2)
        status = nf90_inq_varid (ncid,'NSTDFF',varid_stdff)
        IF (status == nf90_noerr)                                              &
          status = nf90_get_var (ncid, varid_stdff, rc_stdff, nsta2, ncnt2)
      ELSEIF (kcdftyp == ncdf_rass) THEN
        status = nf90_inq_varid (ncid,'MTVIR' ,varid_tv   )
        IF (status == nf90_noerr)                                              &
          status = nf90_get_var (ncid, varid_tv   , rc_tv   , nsta2, ncnt2)
      ENDIF

! preliminary cheap evaluation steps (to reduce number of elements
! ----------------------------------  to be sent to other nodes)

! exploit 'NOAA WIND PROFILER, QUAL. CONTR. RESULTS' bits (BUFR Table 025034)
      status = nf90_inq_varid (ncid,'NWPQ'  ,varid_qci2 )
      IF (status == nf90_noerr) THEN
        ALLOCATE ( nc_qci2  (maxlev,nrep) , STAT=istat )
        status = nf90_get_var (ncid, varid_qci2 , nc_qci2 , nsta2, ncnt2)
        DO irps = 1 , nreproc
          irep  =  irsort(irps)
          DO ilev = 1 , nc_nlev(irep)
            IF ((nc_qci2(ilev,irep) >= 1) .AND. (nc_qci2(ilev,irep) <= 3))     &
              nc_qci (ilev,irep) = 1
          ENDDO
        ENDDO
        DEALLOCATE ( nc_qci2  , STAT=istat )
      ENDIF

!     PRINT *,'press0 ', varid_p, nsta2, ncnt2, maxlev, nrep
!     PRINT '("press1 ",10F8.0)' , (rc_p(ilev,1),ilev=1,5)                     &
!                                , (rc_p(1,irep),irep=1,5) 
    ENDIF

!-------------------------------------------------------------------------------
! Section 3: Produce long buffer arrays 'Xbufsrt' in which the read reports are
!            sorted according to the nodes to which they will be distributed
!-------------------------------------------------------------------------------

    ALLOCATE ( iirlen (nreproc+1) , STAT=istat )
    ALLOCATE ( irrlen (nreproc+1) , STAT=istat )
    ALLOCATE ( iyrlen (nreproc+1) , STAT=istat )

! compute length (number of elements) for each individual TEMP report
! -------------------------------------------------------------------

! total number of scalar elements
    niscal = noffscal(1,2)
    nrscal = noffscal(2,2)
    nyscal = noffscal(3,2)

    IF ((kcdftyp == ncdf_wprof) .OR. (kcdftyp == ncdf_radar_vad)) nintml  = 4
    IF ((kcdftyp == ncdf_wprof) .OR. (kcdftyp == ncdf_radar_vad)) nrealml = 3
    IF  (kcdftyp == ncdf_rass )                                   nintml  = 3
    IF  (kcdftyp == ncdf_rass )                                   nrealml = 2
    DO irps = 1 , nreproc
      irep  =  irsort(irps)
! determine length of integer / real buffer for each report
      iirlen (irps)  =  niscal + nintml *nc_nlev(irep)
      irrlen (irps)  =  nrscal + nrealml*nc_nlev(irep)
      iyrlen (irps)  =  nyscal
    ENDDO

! compute length of buffer arrays (i.e length of all reports together + 1)
! -------------------------------

    ibuflen (1)  =  1
    ibuflen (2)  =  1
    ibuflen (3)  =  1
    DO irps = 1 , nreproc
      ibuflen (1)  =  ibuflen(1) + iirlen (irps)
      ibuflen (2)  =  ibuflen(2) + irrlen (irps)
      ibuflen (3)  =  ibuflen(3) + iyrlen (irps)
    ENDDO

! allocate buffer arrays (only for my_cart_id == 0 here)
! ----------------------

    ALLOCATE ( ibufsrt (ibuflen(1)) , STAT=istat )
    ALLOCATE ( rbufsrt (ibuflen(2)) , STAT=istat )
    ALLOCATE ( ybufsrt (ibuflen(3)) , STAT=istat )

! fill the header part which is common to all obs types 
! into the long buffer arrays 'Xbufsrt'
! -----------------------------------------------------

    CALL obs_cdf_buffer_comhead ( nreproc, irsort, iirlen, irrlen, iyrlen      &
                                , ibuflen, ibufsrt, rbufsrt, ybufsrt )
!   ===========================

! fill the remaining scalar elements into the long buffer arrays 'Xbufsrt'
! ------------------------------------------------------------------------

    iioffs = 0
    iroffs = 0

    DO irps = 1 , nreproc
! the following part is common to all observation types
      irep  =  irsort(irps)
      ibufsrt (iioffs+noffscal(1,1)+ 1)  =  nc_nlev (irep)
      ibufsrt (iioffs+noffscal(1,1)+ 2)  =  nc_dto  (irep)
      ibufsrt (iioffs+noffscal(1,1)+ 3)  =  nc_wce  (irep)
      ibufsrt (iioffs+noffscal(1,1)+ 4)  =  nc_na4  (irep)
      ibufsrt (iioffs+noffscal(1,1)+ 5)  =  nc_tisi (irep)

      iioffs  =  iioffs + iirlen(irps)
      iroffs  =  iroffs + irrlen(irps)
    ENDDO

    DEALLOCATE ( nc_na4   , STAT=istat )
    DEALLOCATE ( nc_tisi  , STAT=istat )
    DEALLOCATE ( nc_dto   , STAT=istat )
    DEALLOCATE ( nc_dto2  , STAT=istat )
    DEALLOCATE ( nc_wce   , STAT=istat )

! fill in the multi-level info into the long buffer arrays 'Xbufsrt'
! ------------------------------------------------------------------

    IF (maxlev >= 1) THEN
      iioffs = 0
      iroffs = 0
      DO irps = 1 , nreproc
        irep  =  irsort(irps)
        nlev = nc_nlev(irep)
        DO ilev = 1 , nlev
          ibufsrt (iioffs+niscal       +ilev)  =  nc_z    (ilev,irep)
          ibufsrt (iioffs+niscal+  nlev+ilev)  =  nc_sinor(ilev,irep)
          ibufsrt (iioffs+niscal+2*nlev+ilev)  =  nc_qci  (ilev,irep)
          rbufsrt (iroffs+nrscal       +ilev)  =  rc_w    (ilev,irep)
          IF ((kcdftyp == ncdf_wprof) .OR. (kcdftyp == ncdf_radar_vad)) THEN
            ibufsrt (iioffs+niscal+3*nlev+ilev)  =  nc_dd   (ilev,irep)
            rbufsrt (iroffs+nrscal+  nlev+ilev)  =  rc_ff   (ilev,irep)
            rbufsrt (iroffs+nrscal+2*nlev+ilev)  =  rc_stdff(ilev,irep)
          ELSEIF (kcdftyp == ncdf_rass) THEN
            rbufsrt (iroffs+nrscal+  nlev+ilev)  =  rc_tv   (ilev,irep)
          ENDIF
        ENDDO
!       PRINT '("bufa0 ",4I7   )' , irps, irep, nlev, iroffs
!       PRINT '("bufaa ",8F10.2)' , (rbufsrt(iroffs+istat),istat=1,8)
!       PRINT '("bufab ",8F10.2)' , (rbufsrt(iroffs+istat),istat=9,16)
!       PRINT '("bufac ",8F10.2)' , (rbufsrt(iroffs+nrscal+3*nlev+istat),istat=9,16)
        iioffs  =  iioffs + iirlen(irps)
        iroffs  =  iroffs + irrlen(irps)
      ENDDO

      DEALLOCATE ( nc_z      , STAT=istat )
      DEALLOCATE ( rc_w      , STAT=istat )
      DEALLOCATE ( nc_sinor  , STAT=istat )
      DEALLOCATE ( nc_qci    , STAT=istat )
      IF ((kcdftyp == ncdf_wprof) .OR. (kcdftyp == ncdf_radar_vad)) THEN
        DEALLOCATE ( nc_dd     , STAT=istat )
        DEALLOCATE ( rc_ff     , STAT=istat )
        DEALLOCATE ( rc_stdff  , STAT=istat )
      ELSEIF (kcdftyp == ncdf_rass) THEN
        DEALLOCATE ( rc_tv     , STAT=istat )
      ENDIF
    ENDIF
    DEALLOCATE ( nc_nlev  , STAT=istat )
!   PRINT '("buf0a ",8F10.2)' , (rbufsrt(istat),istat=1,8)
!   PRINT '("buf0b ",8F10.2)' , (rbufsrt(istat),istat=9,16)

!-------------------------------------------------------------------------------
! Section 4: Distribute the reports
!-------------------------------------------------------------------------------

! determine the number of int / real / char elements 
! to be distributed to the different nodes

    DO irps = 1 , nreproc
      inode =  irnode(irps)
      nodrepn (inode+1) = nodrepn(inode+1) + 1
      nodleni (inode+1) = nodleni(inode+1) + iirlen(irps)
      nodlenr (inode+1) = nodlenr(inode+1) + irrlen(irps)
      nodleny (inode+1) = nodleny(inode+1) + iyrlen(irps)
!     PRINT *,'distr2 ', irps, irnode(irps), iirlen(irps), irrlen(irps), iyrlen(irps)
    ENDDO
!   PRINT *,'distr3a ', nodrepn
!   PRINT *,'distr3b ', nodleni
!   PRINT *,'distr3c ', nodlenr
!   PRINT *,'distr3d ', nodleny

    DEALLOCATE ( iirlen  , STAT=istat )
    DEALLOCATE ( irrlen  , STAT=istat )
    DEALLOCATE ( iyrlen  , STAT=istat )
    DEALLOCATE ( irsort  , STAT=istat )
    DEALLOCATE ( irnode  , STAT=istat )

  ENDIF    ! (my_cart_id == 0)

  IF ((nrep >= 1) .AND. (num_compute > 1)) THEN

    iloclen(1)  =  MAXVAL( nodleni ) + 1
    iloclen(2)  =  MAXVAL( nodlenr ) + 1
    iloclen(3)  =  MAXVAL( nodleny ) + 1

    CALL distribute_values (ibuflen, 3, 0,imp_integers,icomm_cart,ierr)
    CALL distribute_values (iloclen, 3, 0,imp_integers,icomm_cart,ierr)

    IF (my_cart_id /= 0)  ALLOCATE ( ibufsrt (ibuflen(1)) , STAT=istat )
    IF (my_cart_id /= 0)  ALLOCATE ( rbufsrt (ibuflen(2)) , STAT=istat )
    IF (my_cart_id /= 0)  ALLOCATE ( ybufsrt (ibuflen(3)) , STAT=istat )
                          ALLOCATE ( ibufloc (iloclen(1)) , STAT=istat )
                          ALLOCATE ( rbufloc (iloclen(2)) , STAT=istat )
                          ALLOCATE ( ybufloc (iloclen(3)) , STAT=istat )

    CALL obs_cdf_distrib_reports ( ibuflen, ibufsrt, rbufsrt, ybufsrt, ilstidn &
                                 , nodrepn, nodleni, nodlenr, nodleny          &
                                 , nrepl  , nlenli , nlenlr , nlenly           &
                                 , iloclen, ibufloc, rbufloc, ybufloc          &
                                 , num_compute , my_cart_id , icomm_cart       &
                                 , imp_integers, imp_reals )
!   ============================

    DEALLOCATE ( ibufsrt , STAT=istat )
    DEALLOCATE ( rbufsrt , STAT=istat )
    DEALLOCATE ( ybufsrt , STAT=istat )

!-------------------------------------------------------------------------------
! Section 5: Store the local reports in the ODR
!-------------------------------------------------------------------------------

!   PRINT '("buf3a ",8F10.2)' , (rbufloc(istat),istat=1,8)
!   PRINT '("buf3b ",8F10.2)' , (rbufloc(istat),istat=9,16)

    CALL obs_cdf_store_multilev ( kcdftyp , nrepl  , nlenli , nlenlr , nlenly  &
                                , noffscal,          ibufloc, rbufloc, ybufloc &
                                , nodrnew , nexceed )
!   ===========================

    DEALLOCATE ( ibufloc , STAT=istat )
    DEALLOCATE ( rbufloc , STAT=istat )
    DEALLOCATE ( ybufloc , STAT=istat )

  ELSEIF (nrep >= 1) THEN
    nrepl   =  nodrepn(1)
    nlenli  =  nodleni(1)
    nlenlr  =  nodlenr(1)
    nlenly  =  nodleny(1)

    CALL obs_cdf_store_multilev ( kcdftyp , nrepl  , nlenli , nlenlr , nlenly  &
                                , noffscal,          ibufsrt, rbufsrt, ybufsrt &
                                , nodrnew , nexceed )
!   ===========================

    DEALLOCATE ( ibufsrt , STAT=istat )
    DEALLOCATE ( rbufsrt , STAT=istat )
    DEALLOCATE ( ybufsrt , STAT=istat )

  ENDIF


!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_read_profiler
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_read_profiler


!===============================================================================
!+ Module procedure in "src_obs_cdfin_mult" for storing multi-level repo. in ODR
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_store_multilev ( kcdftyp, nrepl, nlenli, nlenlr, nlenly     &
                                  , noffscal      , ibuf  , rbuf  , ybuf       &
                                  , nodrnew, nexceed )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_mult" stores the multi-level
!   reports in the internal ODR arrays.
!
! Method:
!   First, the header part common to all observation types is stored by calling
!   'obs_cdf_store_comhead', and then the rest of the header is added which
!   depends on the observation type.
!   Next, after storing also the cloud information from TEMPs in a surface
!   report, the multi-level observations are copied from the observation type
!   dependent buffers into standardised arrays. Using these, it is decided
!   which vertical levels are discarded depending on pressure (or height),
!   and a vertically sorted list of those levels is compiled which are processed
!   further. Subsequently, the flag patterns are built, and the elements are
!   stored in the multi-level ODR.  
!   To evaluate the blacklist, first, a list of blacklisted vertical intervals
!   for each of the current reports is produced. Then, within the loop over
!   current reports, a blacklist flag 'kblk' is prepared for each level and
!   variable.
!   In addition, a separate surface report is created to accommodate the
!   surface level. (In the multi-level report, 10-m wind, 2-m temperature and
!   humidity are set passive.)
!   Finally, reports without observations are deleted in the ODR.
!
!
! Initial release: Christoph Schraff, 20.12.07
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER (KIND=iintegers) , INTENT (IN)  ::  &
    kcdftyp          ,& ! type of NetCDF input file format (observation type):
                        !  = ncdf_temp      : TEMP
                        !  = ncdf_tempship  : TEMPSHIP
                        !  = ncdf_tempdrop  : TEMP Dropsonde
                        !  = ncdf_pilot     : PILOT (z-levels)
                        !  = ncdf_pilot_p   : PILOT (p-levels)
                        !  = ncdf_wprof     : WIND PROFILER
                        !  = ncdf_radar_vad : RADAR VAD (wind profiles)
                        !  = ncdf_rass      : RASS (virt. temperature profiler)
    nrepl            ,& ! number of            reports  \   at
    nlenli           ,& ! number of integer   elements   \  the
    nlenlr           ,& ! number of real      elements   /  local
    nlenly           ,& ! number of character elements  /   node
    noffscal (3,2)      ! number int/real/char - elements in common header part
                        !                      - scalar elements in total rep.

  INTEGER (KIND=iintegers) , INTENT (IN)  ::  &
    ibuf (nlenli+1)     ! local buffer array containing all integer elements

  REAL    (KIND=wp)        , INTENT (IN)  ::  &
    rbuf (nlenlr+1)     ! local buffer array containing all real elements

  CHARACTER (LEN=10)       , INTENT (IN)  ::  &
    ybuf (nlenly+1)     ! local buffer array containing all character elements

  INTEGER (KIND=iintegers) , INTENT (INOUT)  ::  &
    nodrnew (2)      ,& ! number of new reports
    nexceed (2)         ! number of reports in excess of array size

! Local parameters:
! ----------------

  REAL (KIND=wp)           , PARAMETER  :: &
    fbogus =  0.0_wp    ! if /= 0 then bogus data are produced by applying
                        ! factor 'fbogus' to pre-defined obs 'corrections'

  LOGICAL                  , PARAMETER ::  &
    lmult_shft =.FALSE. ! upper-air observations are not horizontally shifted

  REAL    (KIND=wp)        , PARAMETER ::  &
    c100r = 0.01_wp     ! maximum number of elements to be received

! Local scalars:
! -------------

  LOGICAL                  ::  &
    lblkw  (nrepl+1) ,& ! report blacklisted because sta-ID missing on whitelist
    lchange          ,& ! indicates changes in latest loop of bubblesorting
    lseaobs          ,& ! sea report (ship, buoy, not land report)
    lseaonly         ,& ! use actively only observations from sea platforms
    lamdarml         ,& ! multi-level AMDAR data are processed
    lbogus              ! .true. if bogus data are produced

  INTEGER (KIND=iintegers) ::  &
    nsgob  , nmlob   ,& ! target ODR report indices (single-/ multi-level reps.)
    iioffs (nrepl+1) ,& ! offset of the integer   elements  \  for each
    iroffs (nrepl+1) ,& ! offset of the real      elements   > local
    iyoffs (nrepl+1) ,& ! offset of the character elements  /  report
    jobtyp (nrepl+1) ,& ! observation type (input array for local blacklist)
    jcdtyp (nrepl+1) ,& ! obs  code   type (input array for local whitelist)
    icma             ,& ! pointer index for 'cma' / diagnostic arrays
    niscal , nrscal  ,& ! number of scalar integer / real  elements
    nrepml           ,& ! number of reports which can be stored in the ODR now
    nexceml, nexcesg ,& ! no. of local mult./sfc reports in excess of array size
    irpl             ,& ! loop index over local reports
    ilev             ,& ! loop index over vertical levels (in input reports)
    klev             ,& ! loop index over vertical levels (in ODR)
    intv             ,& ! loop index over blacklisted vertical intervals
    icl              ,& ! loop index over single characters
    maxlev           ,& ! max. number of vertical levels in all local reports
    nlev             ,& ! number of vertical levels in current report
    nlvp             ,& ! number of vertical levels to be processed (for ODR)
    nrejlv           ,& ! number of empty    levels to be rejected
    kobtyp , kcdtyp  ,& ! observation type / observation code type
    mphase           ,& ! phase of flight                     (WMO Table 008004)
    nccl             ,& ! (low or mid-level) cloud cover    (WMO Table 020011)
    nhkey            ,& ! class of cloud base height (VUB WMO Code table 1600)
    ncsig            ,& ! vertical significance of cloud    (WMO Table 008002)
    nclcl  , nclcm   ,& ! derived low / middle cloud amount            (octas)
    nclqf            ,& ! quality flag for derived low cloud amount
    ncc              ,& ! low/middle/high cloud type (VUB T. 0513, 0515, 0509)
    nclwgn , nclwgr  ,& ! combined cloud and weather group (new/old, for verif)
    iob    , job     ,& ! indices of grid point assigned to observation
    iactr               ! stepping for event counters

  INTEGER (KIND=iintegers) ::  &
    ilevsfc          ,& ! level index of surface level (in input array)
    klevsfc          ,& ! level index of surface level (in ODR)
    nlidvof          ,& ! level identity (format of VOF)
    nlidcdf          ,& ! level identity (format of feedback NetCDF)
    nlidin           ,& ! level identity (format of input NetCDF)
    idsurf , idstd   ,& ! indicators for: surface level    / standard level
    idtropo, iduvmax ,& ! indicators for: tropopause level / maximum wind level
    idsigv , idsigq  ,& ! indicators for significant level of: wind / humidity
    idsigt           ,& ! indicators for significant level of: temperature
    nflag            ,& ! total observation flag in ODR
    klvz             ,& ! level with active height obs
    nzpp             ,& ! pressure [pa]
    ilen             ,& ! length of control message
    ilverr           ,& ! nearest standard error level below observation level
    iactx (5)        ,& ! indicates presence of active  data (for variables)
    ipasx (5)        ,& ! indicates presence of passive data (for variables)
    nzaexi           ,& ! -1: only passive data ; 0: no data ; 1: active data
    ilstid_chk       ,& ! length of (part of) station ID that can be checked
    istat  , irm     ,& ! error indicators
    itmp                ! auxilliary buffer
 
  REAL (KIND=wp)           ::  &
    zobhr  (nrepl+1) ,& ! obs  time (input array for local whitelist)
    zclz             ,& ! cloud base height
    zstalt , zzaltsy ,& ! station altitude
    roblat , roblon  ,& ! latitude and longitude in the geograhical system
    rlat   , rlon    ,& ! latitude and longitude in the rotated system
    roberr           ,& ! observation error
    fiperr           ,& ! interpolation weight factor for error level 'ilverr'
    fisd             ,& ! height diff. betw. model orography and sta. altitude
    fisduv           ,& ! modified height difference for 10-m wind data
    fisdtt           ,& ! modified height difference for 2-m temperature data
    fisdrh           ,& ! modified height difference for 2-m humidity data
    c3600r              ! 1 / 3600.0_wp 

! CHARACTER (LEN=25)       :: &
!   yroutine            ! name of this subroutine
  CHARACTER (LEN=ilstid_blk)  :: &
    ystidl  (nrepl+1)   ! observation type (input array for local blacklist)

! Local allocatable arrays: (they include the dimension length 'maxlev+1')
! ------------------------

  LOGICAL                  , ALLOCATABLE ::  &
    lneedt       (:) ,& ! indic. temperature must be available at current level
    lneedq       (:) ,& ! indicates humidity must be available at current level
    lneedv       (:)    ! indicates   wind   must be available at current level

  INTEGER (KIND=iintegers) , ALLOCATABLE ::  &
    ksurfob      (:) ,& ! surface report indicator and report index
    kproclv      (:) ,& ! processing of level: -9: discarded; -1,-2: below surf;
                        !   1: level given by height; 2: normal pressure level
    isrtlvp      (:) ,& ! sorted list of indices of good vertical levels
    izdt         (:) ,& ! time since launch
    izlv         (:) ,& ! level identity (vertical sounding significance)
    iroll        (:) ,& ! roll angle quality flag  (aircraft, BUFR Table 002064)
    iqci         (:) ,& ! quality index            (profiler, BUFR Table 033002,
                        !                                    with use of 025034)
    kblk       (:,:) ,& ! local blacklist (0: ok , 1: blacklisted, for each pro-
                        !                  cessed level and variable separately)
    nflgx        (:) ,& ! total flag for an observation of a variable
    ioberr       (:)    ! active status                  (processed levels only)

  REAL (KIND=wp)           , ALLOCATABLE ::  &
    zpp          (:) ,& ! pressure
    zdlat        (:) ,& ! latitude  displacement
    zdlon        (:) ,& ! longitude displacement
    ztt          (:) ,& ! temperature
    ztd          (:) ,& ! dew-point temperature
    zff          (:) ,& ! wind speed
    zzz          (:) ,& ! height
    zdd          (:) ,& ! wind direction
    zww          (:) ,& ! vertical velocity
    zsinor       (:) ,& ! signal to noise ratio
    zstdff       (:) ,& ! standard deviation wind speed
    ztv          (:) ,& ! virtual temperature            (processed levels only)
    zttl         (:) ,& ! temperature                    (processed levels only)
    ztdl         (:) ,& ! dew-point temperature          (processed levels only)
    zrhw         (:) ,& ! relative humidity over water   (processed levels only)
    zrh          (:) ,& ! model compatible rel humidity  (processed levels only)
    zuu        (:,:) ,& ! zonal      wind                (processed levels only)
    zvv        (:,:) ,& ! meridional wind                (processed levels only)
    zrlat      (:,:) ,& ! latitude  of (processed) observation level
    zrlon      (:,:) ,& ! longitude of (processed) observation level
    zoberr       (:)    ! observation error              (processed levels only)
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_store_multilev
!-------------------------------------------------------------------------------

! yroutine = 'obs_cdf_store_multilev'
  c3600r = c1 / c3600
  lbogus = (ABS( fbogus ) > epsy)

! lrs_temp   =     (kcdftyp == ncdf_temp )    .OR. (kcdftyp == ncdf_tempship)  &
!                                             .OR. (kcdftyp == ncdf_tempdrop)
! lrs_pilot  =     (kcdftyp == ncdf_pilot)    .OR. (kcdftyp == ncdf_pilot_p )
  lamdarml   =     (kcdftyp == ncdf_amdar_ml) .OR. (kcdftyp == ncdf_amdar_vp)
! lprofiler  =     (kcdftyp == ncdf_wprof   ) .OR. (kcdftyp == ncdf_radar_vad) &
!             .OR. (kcdftyp == ncdf_rass    )

! determine the offsets of the different reports in the long (local) array
! ------------------------------------------------------------------------

  iioffs (1) = 0
  iroffs (1) = 0
  iyoffs (1) = 0
  DO irpl = 2 , nrepl
    iioffs (irpl)  =  iioffs(irpl-1) + ibuf(iioffs(irpl-1)+1)
    iroffs (irpl)  =  iroffs(irpl-1) + ibuf(iioffs(irpl-1)+2)
    iyoffs (irpl)  =  iyoffs(irpl-1) + ibuf(iioffs(irpl-1)+3)
!   PRINT *,'bufoff ', irpl, iioffs(irpl), iroffs(irpl), iyoffs(irpl)
  ENDDO
! PRINT '("buf4a ",8F10.2)' , (rbuf(istat),istat=1,8)
! PRINT '("buf4b ",8F10.2)' , (rbuf(istat),istat=9,16)

! total number of scalar elements
! -------------------------------
  niscal = noffscal(1,2)
  nrscal = noffscal(2,2)

!-------------------------------------------------------------------------------
! Section 1: Store report header
!-------------------------------------------------------------------------------

  ALLOCATE ( ksurfob (nrepl+1) , STAT=istat )
  ksurfob (:)  =  0

! store header part, which is common to all observation types, in ODR
! -------------------------------------------------------------------

  CALL obs_cdf_store_comhead ( 1, ntotml, nrepl, nlenli, nlenlr, nlenly        &
                             , ibuf, rbuf, ybuf, nodrnew(1), nexceml, ksurfob )
! ==========================

! store the same header part for single-level surface reports derived from TEMP
!   (here, a surface report header is created for all obs types,
!    but is prepared for cancellation in Section 5 if (ilevsfc == 0);
!    this always applies if obs/code type is not TEMP or PILOT balloon)
! -----------------------------------------------------------------------------

  CALL obs_cdf_store_comhead ( 4, ntotsg, nrepl, nlenli, nlenlr, nlenly        &
                             , ibuf, rbuf, ybuf, nodrnew(2), nexcesg, ksurfob )
! ==========================

! update number of reports which can be stored in the ODR now
! (the following line is equivalent to:  nrepml = nodrnew(1))
  nrepml = nrepl - nexceml

! update counters for insufficient ODR size, for caution messages
  nexceed (1) = nexceed(1) + nexceml
  nexceed (2) = nexceed(2) + nexcesg

! store header part, which is specific to NetCDF
! observation input file type, in multi-level ODR
! -----------------------------------------------

! radiosondes (TEMPs and PILOTs)
  IF (     (kcdftyp == ncdf_temp ) .OR. (kcdftyp == ncdf_tempship)             &
                                   .OR. (kcdftyp == ncdf_tempdrop)             &
      .OR. (kcdftyp == ncdf_pilot) .OR. (kcdftyp == ncdf_pilot_p )) THEN
    DO irpl = 1 , nrepml
      nmlob = ntotml + irpl
      IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 1) /= imiss)                        &
        momlhd (nmlob,nhnlev) = ibuf (iioffs(irpl)+noffscal(1,1)+ 1)
      IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 2) /= imiss)                        &
        momlhd (nmlob,nhrtyp) = ibuf (iioffs(irpl)+noffscal(1,1)+ 2)
      IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 3) /= imiss)                        &
        momlhd (nmlob,nhtrac) = ibuf (iioffs(irpl)+noffscal(1,1)+ 3)
      IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 4) /= imiss)                        &
        momlhd (nmlob,nhna4 ) = ibuf (iioffs(irpl)+noffscal(1,1)+ 4)
!   time significance is not used / stored here currently
!     IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 5) /= imiss)                        &
!       mtimsig (irpl)        = ibuf (iioffs(irpl)+noffscal(1,1)+ 5)
    ENDDO
    IF (     (kcdftyp == ncdf_temp ) .OR. (kcdftyp == ncdf_tempship)           &
                                     .OR. (kcdftyp == ncdf_tempdrop)) THEN
      DO irpl = 1 , nrepml
        nmlob = ntotml + irpl
        IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 6) /= imiss)                      &
          momlhd (nmlob,nhrad ) = ibuf (iioffs(irpl)+noffscal(1,1)+ 6)
      ENDDO
    ENDIF
  ENDIF
! aircrafts (phase of flight)
  IF (lamdarml) THEN
    DO irpl = 1 , nrepml
      nmlob = ntotml + irpl
      IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 1) /= imiss)                        &
        momlhd (nmlob,nhnlev) = ibuf (iioffs(irpl)+noffscal(1,1)+ 1)
      mphase = nibits(nvapoc)
      IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 2) /= imiss)                        &
        mphase                = ibuf (iioffs(irpl)+noffscal(1,1)+ 2)
      CALL MVBITS( mphase , 0 , nvapoc , momlhd(nmlob,nhschr) , nvapbp )
    ENDDO
  ENDIF
! ground-based remote-sensing profilers (wind profiler, radar VAD, RASS)
  IF (     (kcdftyp == ncdf_wprof) .OR. (kcdftyp == ncdf_radar_vad)            &
      .OR. (kcdftyp == ncdf_rass )) THEN
    DO irpl = 1 , nrepml
      nmlob = ntotml + irpl
      IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 1) /= imiss)                        &
        momlhd (nmlob,nhnlev) = ibuf (iioffs(irpl)+noffscal(1,1)+ 1)
      IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 2) /= imiss)                        &
        momlhd (nmlob,nhdt  ) = ibuf (iioffs(irpl)+noffscal(1,1)+ 2)
      IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 3) /= imiss)                        &
        momlhd (nmlob,nhwce ) = ibuf (iioffs(irpl)+noffscal(1,1)+ 3)
      IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 4) /= imiss)                        &
        momlhd (nmlob,nhna4 ) = ibuf (iioffs(irpl)+noffscal(1,1)+ 4)
!   time significance is not used / stored here currently
!     IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 5) /= imiss)                        &
!       mtimsig (irpl)        = ibuf (iioffs(irpl)+noffscal(1,1)+ 5)
    ENDDO
  ENDIF

! get list of blacklisted vertical intervals for each report
! ----------------------------------------------------------

  ilstid_chk  =  MIN( ilstid, ilstid_blk )
  DO irpl = 1 , nrepml
    nmlob = ntotml + irpl
    jobtyp (irpl)  =  momlhd(nmlob,nhobtp)
    jcdtyp (irpl)  =  momlhd(nmlob,nhcode)
    zobhr  (irpl)  =  omlhed(nmlob,nhtime)
    ystidl (irpl)  =  ' '
    ystidl (irpl)  =  yomlhd(nmlob) (1:ilstid_chk)
  ENDDO
  IF (nrepml >= 1) THEN
    ALLOCATE ( blk_loc (nrepml) , STAT=istat )

    CALL obs_cdf_blacklist_local ( nrepml , jobtyp , ilstid_chk , ystidl )
!   ============================

! determine which reports are missing on the whitelists
! -----------------------------------------------------

    CALL obs_cdf_whitelist_local ( nrepml , jobtyp , jcdtyp , zobhr            &
                                 , ilstid_chk , ystidl , lblkw )
!   ============================

  ENDIF

!-------------------------------------------------------------------------------
! Section 2: Report body: Put the observations from the observation type
!            dependent buffers into standardised arrays
!-------------------------------------------------------------------------------

! store cloud part, which is specific to TEMP, in single-level ODR
! ----------------------------------------------------------------

! IF ((kcdftyp == ncdf_temp ) .OR. (kcdftyp == ncdf_tempship)) THEN

!   DO irpl = 1 , nrepml
!     IF (ksurfob(irpl) >= 1) THEN
!       nsgob = ksurfob(irpl)
! restrict to standard observing rules for base of lowest cloud and cloud types
! (yet to be done if required: low / middle / high cloud type)
!       IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 7) == 0) THEN
!         zclz = rbuf(iroffs(irpl)+noffscal(2,1)+ 1)
!         nccl = ibuf(iioffs(irpl)+noffscal(1,1)+ 8)
!         IF ((nccl /= imiss) .AND. (nccl >= 0) .AND. (nccl <= 9)) THEN
!           nhkey = nibits(nvhoc)
!           IF (zclz < rmisschk)          nhkey = 9
!           IF (zclz < 2500._wp-epsy) nhkey = 8
!           IF (zclz < 2000._wp-epsy) nhkey = 7
!           IF (zclz < 1500._wp-epsy) nhkey = 6
!           IF (zclz < 1000._wp-epsy) nhkey = 5
!           IF (zclz <  600._wp-epsy) nhkey = 4
!           IF (zclz <  300._wp-epsy) nhkey = 3
!           IF (zclz <  200._wp-epsy) nhkey = 2
!           IF (zclz <  100._wp-epsy) nhkey = 1
!           IF (zclz <   50._wp-epsy) nhkey = 0
!           IF (nhkey <= 7) THEN
!             osgbdy (nsgob,nbscl ) = REAL( MIN( nccl, 8 ) , wp )
!           ELSEIF (nhkey <= 9) THEN
!             osgbdy (nsgob,nbscl ) = c0
!           ENDIF
!           nclwgr = imdi
!           CALL MVBITS( nhkey , 0 , nvhoc  , nclwgr , nvhbp  )
!           CALL MVBITS( nccl  , 0 , nvnhoc , nclwgr , nvnhbp )
! write combined cloud and weather group to ODR
!           mosgbd (nsgob,nbscwg) = nclwgr
!         ENDIF
!       ENDIF
!     ENDIF
!   ENDDO
! ENDIF

! store cloud part, which is specific to TEMP, in single-level ODR
! ----------------------------------------------------------------

  IF ((kcdftyp == ncdf_temp ) .OR. (kcdftyp == ncdf_tempship)) THEN

    DO irpl = 1 , nrepml
      IF (ksurfob(irpl) >= 1) THEN
        nsgob = ksurfob(irpl)
        ncsig  = ibuf(iioffs(irpl)+noffscal(1,1)+ 7)
        nccl   = ibuf(iioffs(irpl)+noffscal(1,1)+ 8)
        zclz   = rbuf(iroffs(irpl)+noffscal(2,1)+ 1)
! restrict to vertical significance = standard observing rule or clear sky or
!                                     low, middle, or high cloud
!         and cloud cover between 0 and 9
        IF (      (nccl /= imiss) .AND. (nccl >= 0) .AND. (nccl <= 9)          &
            .AND. (     (ncsig == 0) .OR. (ncsig == 7) .OR. (ncsig == 8)       &
                   .OR. (ncsig == 9) .OR. (ncsig == 62))) THEN
          nclwgn = imdi
          nclwgr = imdi
          CALL MVBITS( ncsig , 0 , nxsgoc , nclwgn , nxsgbp )
! cloud base height: convert 'm' into key code
!   (lower limits of bins: 0, 50, 100, 200, 300, 600, 1000, 1500, 2000, 2500 m)
          IF ((zclz < rmdich) .OR. (zclz < c0) .OR. (zclz > 16380._wp)) THEN
            zclz = rmdi
            nhkey = nibits(nvhoc)
          ELSE
                                      nhkey = 9
            IF (zclz < 2500._wp-epsy) nhkey = 8
            IF (zclz < 2000._wp-epsy) nhkey = 7
            IF (zclz < 1500._wp-epsy) nhkey = 6
            IF (zclz < 1000._wp-epsy) nhkey = 5
            IF (zclz <  600._wp-epsy) nhkey = 4
            IF (zclz <  300._wp-epsy) nhkey = 3
            IF (zclz <  200._wp-epsy) nhkey = 2
            IF (zclz <  100._wp-epsy) nhkey = 1
            IF (zclz <   50._wp-epsy) nhkey = 0
            osgbdy (nsgob,nbscbs) = zclz
          ENDIF
          CALL MVBITS( nhkey , 0 , nvhoc , nclwgr , nvhbp )
! low / middle cloud cover
          IF ((nccl == imdi) .OR. (nccl > nibits(nxcloc)) .OR. (nccl < 0))     &
            nccl  =  nibits(nxcloc)
          CALL MVBITS( nccl , 0 , nxcloc , nclwgn , nxclbp )
          CALL MVBITS( nccl , 0 , nvnhoc , nclwgr , nvnhbp )
! derive low and middle cloud cover
          nclcl = nccl
          nclcm = nccl
          nclqf = 0
!     if clear sky (vert. signif. = value not applicable)
          IF ((ncsig == 62) .AND. (nccl == nibits(nvnoc))) THEN
            nclcl = 0
            nclcm = 0
!     if 'MNH'=9 (sky invisible)
          ELSEIF (nclcl == 9) THEN
            nclcl = 8
            nclqf = 1
            nclcm = nibits(nxcloc)
          ELSE
            IF (nhkey <= 7) THEN           ! cloud base height below 2000 m
              nclcm = nibits(nxcloc)
            ELSEIF (nhkey <= 9) THEN       ! cloud base height above 2000 m
              nclcl = 0
            ENDIF
                                           ! vertical significance :
            IF (ncsig == 7) THEN           !   low    cloud  (WMO Table 008002)
              nclcm = nibits(nxcloc)
            ELSEIF (ncsig == 8) THEN       !   middle cloud
              nclcl = 0
            ELSEIF (ncsig == 9) THEN       !   high   cloud
              nclcl = 0
!     (inconsistency between vertical significance and 'MNH')
              IF (nclcm >= 1) nclcm = nibits(nxcloc)
            ENDIF
          ENDIF
!   if cloud base height and vertical significance not consistent:
!   put missing value for low and middle cloud cover
          IF (     ((ncsig == 9) .AND. (nhkey <= 8))                           &
              .OR. ((ncsig == 8) .AND. (nhkey <= 7))                           &
              .OR. ((ncsig <= 7) .AND. (nhkey == 9))) THEN
            nclcl = nibits(nxcloc)
            nclcm = nibits(nxcloc)
          ENDIF
!   broken, scattered, or few cloud
          IF ((nclcl >= 11) .AND. (nclcl <= 13))  nclqf = 1
          IF (nclcl == 12)  nclcl = 6
          IF (nclcl == 11)  nclcl = 2
          IF (nclcl == 13)  nclcl = 1
          IF (nclcm == 12)  nclcm = 6
          IF (nclcm == 11)  nclcm = 2
          IF (nclcm == 13)  nclcm = 1
!   obscured or undefined
          IF (nclcl >   8)  nclcl = nibits(nxcloc)
          IF (nclcm >   8)  nclcm = nibits(nxcloc)
!     (at this point 0 <= nclcl, nclcm <= 8 , or == nibits(nxcloc))
!   write low and middle cloud cover to ODR, 
!   and include low cloud flag in combined flag word
          IF (nclcl < nibits(nxcloc)) osgbdy (nsgob,nbscl) = REAL(nclcl, wp)
          IF (nclcm < nibits(nxcloc)) osgbdy (nsgob,nbscm) = REAL(nclcm, wp)
          CALL MVBITS( nclqf , 0 , nvcqoc , mosgbd(nsgob,nbswwe) , nvcqbp )
! low cloud type
          ncc = ibuf(iioffs(irpl)+noffscal(1,1)+ 9)
          CALL MVBITS( ncc , 0 , nxctoc , nclwgn , nctlbp )
          IF ((ncc >= 30) .AND. (ncc <= 39)) THEN
            ncc = ncc - 30
          ELSEIF ((ncc == 59) .OR. (ncc == 62)) THEN
            ncc = 10
          ELSE
            ncc = nibits(nvcloc)
          ENDIF
          CALL MVBITS( ncc , 0 , nvcloc , nclwgr , nvclbp )
! middle cloud type
          ncc = ibuf(iioffs(irpl)+noffscal(1,1)+10)
          CALL MVBITS( ncc , 0 , nxctoc , nclwgn , nctmbp )
          IF ((ncc >= 20) .AND. (ncc <= 29)) THEN
            ncc = ncc - 20
          ELSEIF ((ncc == 59) .OR. (ncc == 61)) THEN
            ncc = 10
          ELSE
            ncc = nibits(nvcmoc)
          ENDIF
          CALL MVBITS( ncc , 0 , nvcmoc , nclwgr , nvcmbp )
! high cloud type
          ncc = ibuf(iioffs(irpl)+noffscal(1,1)+11)
          CALL MVBITS( ncc , 0 , nxctoc , nclwgn , ncthbp )
          IF ((ncc >= 10) .AND. (ncc <= 19)) THEN
            ncc = ncc - 10
          ELSEIF ((ncc == 59) .OR. (ncc == 60)) THEN
            ncc = 10
          ELSE
            ncc = nibits(nvchoc)
          ENDIF
          CALL MVBITS( ncc , 0 , nvchoc , nclwgr , nvchbp )
! write combined cloud and weather group to ODR
          mosgbd (nsgob,nbsclg) = nclwgn
          mosgbd (nsgob,nbscwg) = nclwgr
        ENDIF
      ENDIF
    ENDDO
  ENDIF

! allocate standardised arrays
! ----------------------------

! get maximum number of vertical levels
  maxlev = 0
  DO irpl = 1 , nrepml
    nmlob = ntotml + irpl
    maxlev = MAX( maxlev, momlhd(nmlob,nhnlev) )
  ENDDO

  ALLOCATE ( zpp     (maxlev+1)   , STAT=istat )
  ALLOCATE ( zdlat   (maxlev+1)   , STAT=istat )
  ALLOCATE ( zdlon   (maxlev+1)   , STAT=istat )
  ALLOCATE ( ztt     (maxlev+1)   , STAT=istat )
  ALLOCATE ( ztd     (maxlev+1)   , STAT=istat )
  ALLOCATE ( zff     (maxlev+1)   , STAT=istat )
  ALLOCATE ( zzz     (maxlev+1)   , STAT=istat )
  ALLOCATE ( zdd     (maxlev+1)   , STAT=istat )
  ALLOCATE ( zww     (maxlev+1)   , STAT=istat )
  ALLOCATE ( zsinor  (maxlev+1)   , STAT=istat )
  ALLOCATE ( zstdff  (maxlev+1)   , STAT=istat )
  ALLOCATE ( izdt    (maxlev+1)   , STAT=istat )
  ALLOCATE ( izlv    (maxlev+1)   , STAT=istat )
  ALLOCATE ( iroll   (maxlev+1)   , STAT=istat )
  ALLOCATE ( iqci    (maxlev+1)   , STAT=istat )
  ALLOCATE ( kproclv (maxlev+1)   , STAT=istat )
  ALLOCATE ( isrtlvp (maxlev+1)   , STAT=istat )

! the following variables are used in section 4 only, and could alternatively be
! re-allocated there for each report with length 'nlvp' (within the report loop)
  ALLOCATE ( kblk    (maxlev+1,4) , STAT=istat )
  ALLOCATE ( nflgx   (maxlev+1)   , STAT=istat )
  ALLOCATE ( ioberr  (maxlev+1)   , STAT=istat )
  ALLOCATE ( ztv     (maxlev+1)   , STAT=istat )
  ALLOCATE ( zttl    (maxlev+1)   , STAT=istat )
  ALLOCATE ( ztdl    (maxlev+1)   , STAT=istat )
  ALLOCATE ( zrhw    (maxlev+1)   , STAT=istat )
  ALLOCATE ( zrh     (maxlev+1)   , STAT=istat )
  ALLOCATE ( zuu     (maxlev+1,1) , STAT=istat )
  ALLOCATE ( zvv     (maxlev+1,1) , STAT=istat )
  ALLOCATE ( zrlat   (maxlev+1,1) , STAT=istat )
  ALLOCATE ( zrlon   (maxlev+1,1) , STAT=istat )
  ALLOCATE ( zoberr  (maxlev+1)   , STAT=istat )
  ALLOCATE ( lneedt  (maxlev+1)   , STAT=istat )
  ALLOCATE ( lneedq  (maxlev+1)   , STAT=istat )
  ALLOCATE ( lneedv  (maxlev+1)   , STAT=istat )
! IF (my_cart_id == 0) PRINT *,'rmisschk', rmisschk, imiss, nrepml

! huge loop over reports
! ----------------------

  DO irpl = 1 , nrepml    
! ~~~~~~~~~~~~~~~~~~~~

    nmlob = ntotml + irpl
    nlev   = momlhd(nmlob,nhnlev)
    kobtyp = momlhd(nmlob,nhobtp)
    kcdtyp = momlhd(nmlob,nhcode)
    iob    = momlhd(nmlob,nhio)
    job    = momlhd(nmlob,nhjo)
    icma   = i_cma ( kobtyp , kcdtyp )
!            =====
    iactr  = 0
    IF (momlhd(nmlob,nhpass) == 0)  iactr  = 1
    zstalt = omlhed(nmlob,nhalt)

! store NetCDF file type dependent observations in standarised temporary arrays
! -----------------------------------------------------------------------------
!   PRINT *,'segmentfault ', irpl, nmlob, maxlev, nlev, kobtyp, kcdtyp, iob, job
    DO ilev = 1 , nlev
      zpp   (ilev) = rmdi
      zdlat (ilev) = rmdi
      zdlon (ilev) = rmdi
      ztt   (ilev) = rmdi
      ztd   (ilev) = rmdi
      zff   (ilev) = rmdi
      zzz   (ilev) = rmdi
      zdd   (ilev) = rmdi
      zww   (ilev) = rmdi
      zsinor(ilev) = rmdi
      zstdff(ilev) = rmdi
      izdt  (ilev) = 0
      izlv  (ilev) = imdi
      iroll (ilev) = imdi
      iqci  (ilev) = imdi
    ENDDO
    IF (     (kcdftyp == ncdf_temp ) .OR. (kcdftyp == ncdf_tempship)           &
                                     .OR. (kcdftyp == ncdf_tempdrop)           &
        .OR. (kcdftyp == ncdf_pilot) .OR. (kcdftyp == ncdf_pilot_p )) THEN
      DO ilev = 1 , nlev
        IF (ABS( rbuf(iroffs(irpl)+nrscal       +ilev) ) < rmisschk) THEN
          zdlat (ilev) =       rbuf (iroffs(irpl)+nrscal       +ilev)
        ENDIF
        IF (ABS( rbuf(iroffs(irpl)+nrscal+  nlev+ilev) ) < rmisschk) THEN
          zdlon (ilev) =       rbuf (iroffs(irpl)+nrscal+  nlev+ilev)
        ENDIF
        IF (ABS( rbuf(iroffs(irpl)+nrscal+2*nlev+ilev) ) < rmisschk) THEN
          zff   (ilev) =       rbuf (iroffs(irpl)+nrscal+2*nlev+ilev)
        ENDIF
        IF (ABS( rbuf(iroffs(irpl)+nrscal+3*nlev+ilev) ) < rmisschk) THEN
          zpp   (ilev) =       rbuf (iroffs(irpl)+nrscal+3*nlev+ilev)
        ENDIF
        IF (     ibuf(iioffs(irpl)+niscal+  nlev+ilev)  /= imiss   ) THEN
          zdd   (ilev) = REAL( ibuf (iioffs(irpl)+niscal+  nlev+ilev) , wp )
        ENDIF
        IF (     ibuf(iioffs(irpl)+niscal+2*nlev+ilev)  /= imiss   ) THEN
          zzz   (ilev) = REAL( ibuf (iioffs(irpl)+niscal+2*nlev+ilev) , wp )
        ENDIF
        IF (     ibuf(iioffs(irpl)+niscal+3*nlev+ilev)  /= imiss   ) THEN
          izdt  (ilev) =       ibuf (iioffs(irpl)+niscal+3*nlev+ilev)
        ENDIF
        IF (                  ( ibuf(iioffs(irpl)+niscal+ilev) /= imiss  )     &
            .AND. (.NOT. BTEST( ibuf(iioffs(irpl)+niscal+ilev),ilv_miss ))) THEN
          izlv  (ilev) =       ibuf (iioffs(irpl)+niscal       +ilev)
        ENDIF
      ENDDO
      IF ((kcdftyp == ncdf_temp ) .OR. (kcdftyp == ncdf_tempship)              &
                                  .OR. (kcdftyp == ncdf_tempdrop)) THEN
        DO ilev = 1 , nlev
          IF (ABS( rbuf(iroffs(irpl)+nrscal+4*nlev+ilev) ) < rmisschk) THEN
            ztt   (ilev) =       rbuf (iroffs(irpl)+nrscal+4*nlev+ilev)
          ENDIF
          IF (ABS( rbuf(iroffs(irpl)+nrscal+5*nlev+ilev) ) < rmisschk) THEN
            ztd   (ilev) =       rbuf (iroffs(irpl)+nrscal+5*nlev+ilev)
          ENDIF
        ENDDO
      ENDIF
    ENDIF
    IF (lamdarml) THEN
      DO ilev = 1 , nlev
        IF (ABS( rbuf(iroffs(irpl)+nrscal       +ilev) ) < rmisschk) THEN
          zdlat (ilev) =       rbuf (iroffs(irpl)+nrscal       +ilev)
        ENDIF
        IF (ABS( rbuf(iroffs(irpl)+nrscal+  nlev+ilev) ) < rmisschk) THEN
          zdlon (ilev) =       rbuf (iroffs(irpl)+nrscal+  nlev+ilev)
        ENDIF
        IF (ABS( rbuf(iroffs(irpl)+nrscal+2*nlev+ilev) ) < rmisschk) THEN
          zff   (ilev) =       rbuf (iroffs(irpl)+nrscal+2*nlev+ilev)
        ENDIF
        IF (ABS( rbuf(iroffs(irpl)+nrscal+3*nlev+ilev) ) < rmisschk) THEN
          ztt   (ilev) =       rbuf (iroffs(irpl)+nrscal+3*nlev+ilev)
        ENDIF
        IF (ABS( rbuf(iroffs(irpl)+nrscal+4*nlev+ilev) ) < rmisschk) THEN
          ztd   (ilev) =       rbuf (iroffs(irpl)+nrscal+4*nlev+ilev)
        ENDIF
        IF (     ibuf(iioffs(irpl)+niscal+  nlev+ilev)  /= imiss   ) THEN
          zdd   (ilev) = REAL( ibuf (iioffs(irpl)+niscal+  nlev+ilev) , wp )
        ENDIF
        IF (     ibuf(iioffs(irpl)+niscal+2*nlev+ilev)  /= imiss   ) THEN
          zzz   (ilev) = REAL( ibuf (iioffs(irpl)+niscal+2*nlev+ilev) , wp )
        ENDIF
        IF (     ibuf(iioffs(irpl)+niscal       +ilev)  /= imiss   ) THEN
          iroll (ilev) =       ibuf (iioffs(irpl)+niscal       +ilev)
        ENDIF
      ENDDO
    ENDIF
    IF (     (kcdftyp == ncdf_wprof) .OR. (kcdftyp == ncdf_radar_vad)          &
        .OR. (kcdftyp == ncdf_rass )) THEN
      DO ilev = 1 , nlev
        IF (ABS( rbuf(iroffs(irpl)+nrscal       +ilev) ) < rmisschk) THEN
          zww   (ilev) =       rbuf (iroffs(irpl)+nrscal       +ilev)
        ENDIF
        IF (     ibuf(iioffs(irpl)+niscal       +ilev)  /= imiss   ) THEN
          zzz   (ilev) = REAL( ibuf (iioffs(irpl)+niscal       +ilev) , wp )
        ENDIF
        IF (     ibuf(iioffs(irpl)+niscal+  nlev+ilev)  /= imiss   ) THEN
          zsinor(ilev) = REAL( ibuf (iioffs(irpl)+niscal+  nlev+ilev) , wp )
        ENDIF
        IF (     ibuf(iioffs(irpl)+niscal+2*nlev+ilev)  /= imiss   ) THEN
          iqci  (ilev) =       ibuf (iioffs(irpl)+niscal+2*nlev+ilev)
        ENDIF
      ENDDO
      IF ((kcdftyp == ncdf_wprof) .OR. (kcdftyp == ncdf_radar_vad)) THEN
        DO ilev = 1 , nlev
          IF (ABS( rbuf(iroffs(irpl)+nrscal+  nlev+ilev) ) < rmisschk) THEN
            zff   (ilev) =       rbuf (iroffs(irpl)+nrscal+  nlev+ilev)
          ENDIF
          IF (ABS( rbuf(iroffs(irpl)+nrscal+2*nlev+ilev) ) < rmisschk) THEN
            zstdff(ilev) =       rbuf (iroffs(irpl)+nrscal+2*nlev+ilev)
          ENDIF
          IF (     ibuf(iioffs(irpl)+niscal+3*nlev+ilev)  /= imiss   ) THEN
            zdd   (ilev) = REAL( ibuf (iioffs(irpl)+niscal+3*nlev+ilev) ,wp)
          ENDIF
        ENDDO
      ELSEIF (kcdftyp == ncdf_rass) THEN
        DO ilev = 1 , nlev
          IF (ABS( rbuf(iroffs(irpl)+nrscal+  nlev+ilev) ) < rmisschk) THEN
            ztt   (ilev) =       rbuf (iroffs(irpl)+nrscal+  nlev+ilev)
          ENDIF
        ENDDO
      ENDIF
    ENDIF

!-------------------------------------------------------------------------------
! Section 3: Decide which vertical levels shall be discarded, and
!            compile a sorted list of levels which shall be processed further
!-------------------------------------------------------------------------------

! check vertical coordinate and level (pressure or height)
! --------------------------------------------------------

    DO ilev = 1 , nlev
      kproclv (ilev) =  2

! if height is given instead of pressure as a vertical coordinate:
      IF (      (zpp(ilev) <= rmdich)                                          &
          .AND. (zzz(ilev) >  rmdich) .AND. (kobtyp /= ntemp)) THEN
! AMDAR: calculate pressure from height according to standard atmosphere
        IF (lamdarml) THEN

          CALL std_atmosphere ( 'z', zzz(ilev), 'p', zpp(ilev), r_g, r_d, irm )
!         ===================

          IF (irm /= 0)  zpp (ilev) = rmdi
        ELSE
! others: use model fields to determine pressure (level) at given height level
          zpp (ilev) = f_z2p ( zzz(ilev), ke, r_hhl(iob,job,:), r_p(iob,job,:) &
                             , r_t_ll(iob,job), r_g, r_d, rmdi )
!                      =====
        ENDIF
        kproclv (ilev) = 1
! discard z-levels below model orography or above top model level
        IF (zpp(ilev) < rmdich) THEN
          kproclv (ilev) = -9
          neventd (nelnop,icma) = neventd(nelnop,icma) + iactr
        ENDIF

      ELSEIF (zpp(ilev) <= rmdich) THEN
! vertical coordinate does not exist
        kproclv (ilev) = -9
        neventd (nelmis,icma) = neventd(nelmis,icma) + iactr
      ENDIF

! discard any levels above 'rpplim' or at more than 40 m below station height
      IF (    ((zpp(ilev) > rmdich) .AND. (zpp(ilev) < rpplim))                &
         .OR. ((zzz(ilev) > rmdich) .AND. (zzz(ilev)+40._wp < zstalt)          &
                                    .AND. (zstalt > rmdich))) THEN
        kproclv (ilev) = -9
        neventd (nelext,icma) = neventd(nelext,icma) + iactr

! discard non-mandatory levels except for levels
!   - below 'pminsigt' (i.e. p-obs >= pminsigt) which contain temperature, or
!   - below 'pminsigv' (i.e. p-obs >= pminsigv) which contain wind obs
      ELSEIF (    (zpp(ilev) > rmdich) .AND. (     (zpp(ilev) < pminsigv-c05)  &
                                              .OR. (zpp(ilev) < pminsigt-c05)) &
             .AND.(    (kcdftyp == ncdf_temp ) .OR. (kcdftyp == ncdf_tempship) &
                                               .OR. (kcdftyp == ncdf_tempdrop) &
                  .OR. (kcdftyp == ncdf_pilot) .OR. (kcdftyp == ncdf_pilot_p)) &
             ) THEN
        nzpp = NINT( zpp(ilev) )
        IF (      (nzpp-85000 /= 0) .AND. (nzpp-70000 /= 0)                    &
            .AND. (nzpp-50000 /= 0) .AND. (nzpp-40000 /= 0)                    &
            .AND. (nzpp-30000 /= 0) .AND. (nzpp-25000 /= 0)                    &
            .AND. (nzpp-20000 /= 0) .AND. (nzpp-15000 /= 0)                    &
            .AND. (nzpp-10000 /= 0) .AND. (nzpp- 7000 /= 0)                    &
            .AND. (nzpp- 5000 /= 0) .AND. (nzpp- 3000 /= 0)                    &
            .AND. (nzpp- 2500 /= 0) .AND. (nzpp- 2000 /= 0)                    &
            .AND. (nzpp- 1000 /= 0)) THEN
          IF (      ((ztt(ilev) < rmdich) .OR. (zpp(ilev) < pminsigt-c05))     &
              .AND. ((zff(ilev) < rmdich) .OR. (zpp(ilev) < pminsigv-c05))) THEN
            kproclv (ilev) = -9
            neventd (nelsig,icma) = neventd(nelsig,icma) + iactr
          ENDIF
        ENDIF
! discard aircraft data below 'model surface pressure + 1hPa'
      ELSEIF ((lamdarml) .AND. (zpp(ilev) > rmdich)                            &
                         .AND. (zpp(ilev) > r_ps(iob,job)+c1)) THEN
        kproclv (ilev) = -9
        neventd (nelsig,icma) = neventd(nelsig,icma) + iactr
      ENDIF
    ENDDO

! check for several surface levels
! --------------------------------

    ilevsfc = 0
    IF (     (kcdftyp == ncdf_temp ) .OR. (kcdftyp == ncdf_tempship)           &
                                     .OR. (kcdftyp == ncdf_tempdrop)           &
        .OR. (kcdftyp == ncdf_pilot) .OR. (kcdftyp == ncdf_pilot_p )) THEN
      DO ilev = 1 , nlev
! if surface level
        IF (     (BTEST( izlv(ilev),ilv_sfc )) .AND. (izlv(ilev) /= imdi)      &
           .AND. (kproclv(ilev) /= -9)) THEN
          IF (ilevsfc == 0) THEN
            ilevsfc = ilev
          ELSE
! several surface levels: take the one which is closer to the station altitude
            IF (     (zzz(ilev) <= rmdich) .OR. (zstalt <= rmdich)             &
                .OR. (ABS(zzz(ilev)-zstalt) > ABS(zzz(ilevsfc)-zstalt))) THEN
              kproclv (ilev)    = -9
            ELSE
              kproclv (ilevsfc) = -9
              ilevsfc           = ilev
            ENDIF
! for statistics and control message only
            neventd(nelsfc,icma) = neventd(nelsfc,icma) + iactr
            ilen = 2 + istrej
            IF (nacout+ilen <= nmxoln) THEN
              outbuf(nacout+1) = ilen
              outbuf(nacout+2) = nfmt5
              DO icl = 1 , istrej
                outbuf(nacout+2+icl) = ICHAR( yomlhd(nmlob) (icl:icl) )
              ENDDO
              nacout  = nacout + ilen
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      IF ((ilevsfc >= 1) .AND.(MIN( zzz(MAX(ilevsfc,1)),zstalt ) > rmdich)) THEN
        IF (ABS( zzz(ilevsfc) -zstalt ) > 99.9_wp)  ilevsfc = 0
      ENDIF
    ENDIF

! set height of surface level equal to station height (for TEMP only)
    IF (      (ilevsfc >= 1)    .AND. (zzz(MAX(ilevsfc,1)) < rmdich)           &
        .AND. (kobtyp == ntemp) .AND. (zstalt > rmdich))  zzz (ilevsfc) = zstalt

! discard any levels below the surface level
    IF (ilevsfc >= 1) THEN
      DO ilev = 1 , nlev
        IF ((zpp(ilev) > zpp(ilevsfc)+c05) .AND. (kproclv(ilev) /= -9)) THEN
          kproclv (ilev) = -9
          neventd (nelext,icma) = neventd(nelext,icma) + iactr
        ENDIF
      ENDDO
    ENDIF

! compile a list of indices of those vertical levels which will be processed
! --------------------------------------------------------------------------

    nlvp = 0
    DO ilev = 1 , nlev
      IF (kproclv(ilev) > -9) THEN
        nlvp  =  nlvp + 1
        isrtlvp (nlvp) = ilev
      ENDIF
    ENDDO

! sort this list of indices in the vertical
! (iterative sorting with bubblesort algorithm)
! ---------------------------------------------
! (this replaces 'obs_sort_levels' in 'src_obs_proc_aof.f90')

    lchange = .TRUE.
    DO WHILE (lchange)
      lchange = .FALSE.
      DO klev = 1 , nlvp-1
        IF (zpp(isrtlvp(klev)) < zpp(isrtlvp(klev+1))) THEN
          itmp             = isrtlvp (klev+1)
          isrtlvp (klev+1) = isrtlvp (klev)
          isrtlvp (klev)   = itmp
          lchange = .TRUE.
        ENDIF
      ENDDO
    ENDDO

! get surface level
    klevsfc = 0
    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)
      IF (ilev == ilevsfc)  klevsfc = klev
!     PRINT *,'press6 ', yomlhd(nmlob), nlvp, klev, ilev, zpp(ilev), ztt(ilev)
    ENDDO
!   PRINT *,'klevsfc ',yomlhd(nmlob),klevsfc,ilevsfc, zpp(MAX(ilevsfc,1))      &
!                                                   , ztt(MAX(ilevsfc,1))

! check if number of levels exceeds ODR size, and if so, reset number of levels
! -----------------------------------------------------------------------------

    IF (nlvp > maxmlv) THEN
! insufficient ODR size is the only report and data events, which is updated
! even for reports set to passive previously
      neventd (nelodr,icma) = neventd(nelodr,icma) + nlvp - maxmlv
      ilen = 3 + istrej
      IF (nacout+ilen <= nmxoln) THEN
        outbuf(nacout+1) = ilen
        outbuf(nacout+2) = nfmt4
        outbuf(nacout+3) = nlvp
        DO icl = 1 , istrej
          outbuf(nacout+3+icl) = ICHAR( yomlhd(nmlob) (icl:icl) )
        ENDDO
        nacout  = nacout + ilen
      ENDIF
    ENDIF
    nlvp  =  MIN( nlvp , maxmlv )

! in order to accommodate TEMP (or PILOT) Parts A, B, C, D in 1 report
! (in obs_cdf_redundancy) limit each part (mainly part B !) to 'maxmlv-15'
! (note: this makes use of the DWD-specific element 'nhkz')
! ------------------------------------------------------------------------

    IF ((momlhd(nmlob,nhkz) > 500) .AND. (momlhd(nmlob,nhkz) < 800))           &
      nlvp  =  MIN( nlvp , maxmlv - 15 )

! adjust the number of vertical levels in the ODR header
! to those levels which are filled in the ODR body
! ------------------------------------------------------

    momlhd (nmlob,nhnlev) = nlvp

! blacklist
! ---------
! TYPE blacklist_loc
!   INTEGER :: ndim            ! number of blacklisted intervals in report
!   INTEGER :: kvar (maxintv)  ! variable type of blacklisted interval
!                              !   1: wind, 2: geopot, 3: temperat., 4: humidity
!   REAL    :: plow (maxintv)  ! lower boundary (pressure) of blacklist interval
!   REAL    :: pup  (maxintv)  ! upper boundary (pressure) of blacklist interval
! END TYPE blacklist_loc
! TYPE (blacklist_loc) , DIMENSION(:), ALLOCATABLE  :: blk_loc(nrepml)
! ALLOCATE ( blk_loc (nrepml) , STAT=istat )

!-------------------------------------------------------------------------------
! Section 4: Determine and store the elements of the multi-level ODR
!            (for the 'good' vertical levels only)
!-------------------------------------------------------------------------------

! ----------------------------------------------------
! auxilliary: observation blacklist for current report
! ----------------------------------------------------
! for each level and variable, determine whether it is in a blacklisted interval

    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)
      kblk (klev,1) = 0
      kblk (klev,2) = 0
      kblk (klev,3) = 0
      kblk (klev,4) = 0
      DO intv = 1 , blk_loc(irpl)% ndim 
        IF (      (zpp(ilev) <= blk_loc(irpl)%plow(intv) +epsy)                &
            .AND. (zpp(ilev) >= blk_loc(irpl)%pup (intv) -epsy)) THEN
          kblk (klev,blk_loc(irpl)%kvar(intv)) = 1
        ENDIF
      ENDDO
    ENDDO
    IF (lblkw(irpl)) THEN
      DO klev = 1 , nlvp
        kblk (klev,1) = 1
        kblk (klev,2) = 1
        kblk (klev,3) = 1
        kblk (klev,4) = 1
      ENDDO
    ENDIF

! -----------------------------
! level ID / level significance
! -----------------------------
! convert level identity from BUFR Table 0 08 042 into:
! - level identity as defined in VOF
! - level identity as defined for NetCDF feedback file format
! -----------------------------------------------------------

    DO klev = 1 , nlvp
      nlidvof  =  0
      nlidcdf  =  0
      lneedt (klev)  =  .FALSE.
      lneedq (klev)  =  .FALSE.
      lneedv (klev)  =  .FALSE.
      IF (izlv(isrtlvp(klev)) /= imdi) THEN
        !   get    level significance
        nlidin  =  izlv(isrtlvp(klev))
        idsurf   =  ibit1 ( nlidin, ilv_sfc  )
        idstd    =  ibit1 ( nlidin, ilv_std  )
        idtropo  =  ibit1 ( nlidin, ilv_tropo)
        iduvmax  =  ibit1 ( nlidin, ilv_max  )
        idsigt   =  ibit1 ( nlidin, ilv_sigt )
        idsigq   =  ibit1 ( nlidin, ilv_sigq )
        idsigv   =  ibit1 ( nlidin, ilv_sigv )
        !   NetCDF level significance
        nlidcdf  =  insert( 0      , idsurf , LS_SURFACE )
        nlidcdf  =  insert( nlidcdf, idstd  , LS_STANDARD)
        nlidcdf  =  insert( nlidcdf, idtropo, LS_TROPO   )
        nlidcdf  =  insert( nlidcdf, iduvmax, LS_MAX     )
        nlidcdf  =  insert( nlidcdf, IOR( IOR(idsigt,idsigq), idsigv ), LS_SIGN)
        !   VOF    level significance
        nlidvof  =  insert( 0      , idsurf , nvlidp(7) )
        IF (zpp(isrtlvp(klev)) >= 10000.0_wp-epsy) THEN
          nlidvof = insert( nlidvof, idstd  , nvlidp(6) )
        ELSE
          nlidvof = insert( nlidvof, idstd  , nvlidp(4) )
        ENDIF
        nlidvof  =  insert( nlidvof, idtropo, nvlidp(2) )
        nlidvof  =  insert( nlidvof, iduvmax, nvlidp(1) )
        nlidvof  =  insert( nlidvof, idsigv , nvlidp(8) )
        nlidvof  =  insert( nlidvof, IOR( idsigt,idsigq ), nvlidp(9) )
        ! decide whether a certain variable must be available at a certain level
        IF (idsurf == 1)  idstd = 1
        lneedt (klev)  =  (idsigt == 1) .OR. (idtropo == 1) .OR. (idstd == 1)
        lneedq (klev)  =  (idsigq == 1)                     .OR. (idstd == 1)
        lneedv (klev)  =  (idsigv == 1) .OR. (iduvmax == 1) .OR. (idstd == 1)
      ENDIF

      ! fill ODR : level ID and level significance
      ! --------
      momlbd (nmlob,klev,nbtlid) = nlidvof
      momlbd (nmlob,klev,nbtlsg) = nlidcdf

      ! fill ODR : initialize status and quality control flag words
      ! --------
      momlbd (nmlob,klev,nbterr) =  0
      momlbd (nmlob,klev,nbtqcf) =  0
    ENDDO

! --------------------------------------------
! latitude / longitude [grid pt. units] / time
! --------------------------------------------
! Note: Consideration of changes in the horizontal observation location
! ----                   as a function of the vertical level
!       - is always done here (i.e. for 'nbtzio', 'nbtzjo')
!       - is not done when pressure is computed from height levels (call'f_z2p')
!       - depends on 'lmult_shft' when computing the rotation of the wind
!         components (call 'uv2uvrot_vec', see below) !

    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)
      roblat  =  omlhed(nmlob,nhjlat)
      roblon  =  omlhed(nmlob,nhilon)
      IF ((zdlat(ilev) > rmdich) .AND. (zdlon(ilev) > rmdich)) THEN
        !   (AMDAR: zdlat, zdlon are full lat./long., i.e. not a displacement)
        IF (lamdarml)  roblat  =  c0
        IF (lamdarml)  roblon  =  c0
        roblat  =  roblat + zdlat(ilev)
        roblon  =  roblon + zdlon(ilev)
      ENDIF

      rlon = rla2rlarot (roblat, roblon, r_pollat, r_pollon, r_polgam)
      rlat = phi2phirot (roblat, roblon, r_pollat, r_pollon)
!            ==========

      ! fill ODR     (observation position in grid point units (total area !)
      ! --------      and (level-dependent) threshold Q.C. flag )
      omlbdy (nmlob,klev,nbtzio) = c1 + (rlon - r_startlon_tot) /r_dlon
      omlbdy (nmlob,klev,nbtzjo) = c1 + (rlat - r_startlat_tot) /r_dlat
      omlbdy (nmlob,klev,nbttim) = omlhed(nmlob,nhtime) + izdt(ilev) *c3600r
    ENDDO

! -------------------
! height and pressure
! -------------------

    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)
      nflag   =  0
      nflgx  (klev)  =   0
      zoberr (klev)  =  c0
      ioberr (klev)  =   1

      IF (     (zzz(ilev) <= rmdich) .OR. (kproclv(ilev) <= -9)                &
          .OR. (zpp(ilev) <= rmdich)) THEN
!         .OR. (zpp(ilev) <= rmdich) .OR. (zstalt <= rmdich)) THEN
        IF ((nzex(kobtyp) >= 1) .AND. (ioberr(klev) == 1))                     &
          neventd (nepmis,icma) = neventd(nepmis,icma) + iactr
        ioberr (klev)  =  0
      ELSE
        !---   if level below surface
        IF (zstalt > rmdich) THEN
          IF (zzz(ilev) <= zstalt-c2) THEN
            IF (ilev /= ilevsfc) THEN
              !   level below surface, means that also T,q,uv shall
              !   not be used actively !  ( --> set kproclv=-kproclv )
              IF (ioberr(klev) == 1)  neventd (nelext,icma) =                  &
                                      neventd (nelext,icma) + iactr
              nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(4) )
              nflag          =  IBSET( nflag      , nvflbp )
              kproclv (ilev) =  - kproclv(ilev)
            ENDIF
            ioberr (klev)  =  0
          ENDIF
        ENDIF
        !---   if pressure is derived from height (or vice versa)
        IF (ABS( kproclv(ilev) ) == 1) THEN
          IF (lamdarml)                                                        &
            nflgx (klev) =  IBSET( nflgx(klev), nvfbps(5) )
          IF ((ABS( kproclv(ilev) ) == 1) .AND. (.NOT. lamdarml))              &
            nflgx (klev) =  IBSET( nflgx(klev), nvfbps(6) )
          ioberr (klev)  =  0
        ENDIF
        IF (ABS( kproclv(ilev) ) == 2) THEN
          !---   geopotential without temperature is set passive
          IF ((ztt(ilev) < rmdich) .OR. (ztt(ilev) < c1)) THEN
            IF (ioberr(klev) == 1)  neventd (nepflg,icma) =                    &
                                    neventd (nepflg,icma) + iactr
            nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(5) )
            ioberr (klev)  =  0
          !---   if geopotential blacklisted
          ELSEIF (kblk(klev,2) == 1) THEN
            IF (ioberr(klev) == 1)  neventd (nepflg,icma) =                    &
                                    neventd (nepflg,icma) + iactr
            nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(2) )
            ioberr (klev)  =  0
          ENDIF
        ENDIF
      ENDIF

      ! fill ODR (except for flags)
      ! --------
      omlbdy (nmlob,klev,nbtp  ) = zpp   (ilev)
      omlbdy (nmlob,klev,nbtz  ) = zzz   (ilev)
      omlbdy (nmlob,klev,nbtzer) = zoberr(klev)
      momlbd (nmlob,klev,nbtflg) = nflag
      !   bogus observations for semi-idealised tests
      IF (lbogus)                                                              &
        omlbdy (nmlob,klev,nbtp) = omlbdy(nmlob,klev,nbtp) + 100._wp*fbogus
      !   auxilliary ODR element
      omlbdy (nmlob,klev,nbtlop) = LOG ( omlbdy(nmlob,klev,nbtp) )

!     IF (zpp(ilev) > 90000.0_wp)                                              &
!       PRINT *,'zzz4 ', yomlhd(nmlob), ilev, omlbdy(nmlob,klev,nbtzer)        &
!                                           , omlbdy(nmlob,klev,nbtz  )
    ENDDO

    !   if (.not. luse_mlz) then
    !   set 'height'/passive flag for all obs except for TEMP surface level,
    !        or (if no surface level exists) for lowest active TEMP z-obs
    klvz  =  klevsfc
    IF ((.NOT. luse_mlz) .AND. (klevsfc == 0) .AND. (kobtyp == ntemp)) THEN
      !   if no surface obs level, get lowest obs level with active p-z obs
      klvz = 1
      DO WHILE ((klvz < nlvp) .AND.(.NOT.BTEST(momlbd(nmlob,klvz,nbterr),nvrz)))
        klvz  =  klvz + 1
      ENDDO
      IF ((klvz == nlvp) .AND. (.NOT.BTEST( momlbd(nmlob,klvz,nbterr),nvrz ))) &
        klvz  =  0
    ENDIF
    IF (kobtyp /= ntemp)  klvz = 0
    DO klev = 1 , nlvp
      IF ((.NOT. luse_mlz) .AND. (klev /= klvz)) THEN
        ilev = isrtlvp(klev)
        IF ((zzz(ilev) > rmdich) .AND. (ABS( kproclv(ilev) ) == 2)) THEN
          IF (ioberr(klev) == 1)  neventd (nepflg,icma) =                      &
                                  neventd (nepflg,icma) + iactr
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(4) )
          ioberr (klev)  =  0
        ENDIF
      ENDIF

      ! fill ODR: flags only
      ! --------
      momlbd (nmlob,klev,nbterr) = insert( momlbd(nmlob,klev,nbterr)           &
                                         , ioberr(klev), nvrz )
      momlbd (nmlob,klev,nbtflg) = insert( momlbd(nmlob,klev,nbtflg)           &
                                         , nflgx (klev), nvfzbp )
    ENDDO

! -----------
! temperature
! -----------

    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)
      nflgx  (klev)  =   0
      zoberr (klev)  =  c0
      ioberr (klev)  =   1

! if temperature missing or below absolute minimum (no flag in 'nflgx' set)
      IF ((ztt(ilev) < rmdich) .OR. (ztt(ilev) < c1)) THEN
        IF ((ntex(kobtyp) >= 1) .AND. (lneedt(klev)))                          &
          neventd (netmis,icma) = neventd(netmis,icma) + iactr
        ztt    (ilev)  =  rmdi
        ioberr (klev)  =  0
      ELSE
! if level below surface
        IF ((kproclv(ilev) == -2) .OR. (kproclv(ilev) == -1)) THEN
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(4) )
          ioberr (klev)  =  0
        ENDIF
! if temperature blacklisted
        IF (kblk(klev,3) == 1) THEN
          IF (ioberr(klev) == 1)  neventd (netflg,icma) =                      &
                                  neventd (netflg,icma) + iactr
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(2) )
          ioberr (klev)  =  0
        ENDIF
! if temperature gross error
        IF (      (ztt(ilev) < tmelt-90._wp)                                   &
            .OR.  (ztt(ilev) > tmelt+60._wp)                                   &
            .OR. ((ztt(ilev) > tmelt+20._wp).AND. (zpp(ilev) < 70000._wp))     &
            .OR. ((ztt(ilev) > tmelt+ 5._wp).AND. (zpp(ilev) < 50000._wp))     &
            .OR. ((ztt(ilev) > tmelt- 5._wp).AND. (zpp(ilev) < 40000._wp))) THEN
          IF (ioberr(klev) == 1)  neventd (netext,icma) =                      &
                                  neventd (netext,icma) + iactr
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(3) )
          ioberr (klev)  =  0
        ENDIF
      ENDIF
      IF (      (kcdftyp == ncdf_rass )                                        &
          .AND. ((ztt(ilev) > rmdich) .OR. (zww(ilev) > rmdich))) THEN
! RASS profiler quality information (also valid for vertical wind)
        IF ((iqci (ilev) < 0) .OR. (iqci(ilev) > nibits(nvfboc(1))))           &
          iqci (ilev) = nibits(nvfboc(1))
        CALL MVBITS( iqci(ilev) , 0 , nvfboc(1) , nflgx(klev) , nvfbps(1) )
        IF (iqci (ilev) == 1) THEN
          IF (ioberr(klev) == 1)  neventd (netflg,icma) =                      &
                                  neventd (netflg,icma) + iactr
          ioberr (klev)  =  0
        ENDIF
! RASS profiler signal to noise ratio (also valid for vertical wind)
        IF ((zsinor(ilev) > rmdich) .AND. (zsinor(ilev) < -20.0_wp)) THEN
          IF (ioberr(klev) == 1)  neventd (netflg,icma) =                      &
                                  neventd (netflg,icma) + iactr
          CALL MVBITS( 1 , 0 , nvfboc(1) , nflgx(klev) , nvfbps(1) )
          ioberr (klev)  =  0
        ENDIF
      ENDIF
    ENDDO

! Convert RASS / SODAR vitual temperatures into normal air temperatures
! ---------------------------------------------------------------------
!TODO ( --> 'tv2t' is not yet vectorised)
!   IF (kcdftyp == ncdf_rass) THEN
!     DO klev = 1 , nlvp
!       ilev = isrtlvp(klev)
!       IF ((ztt(ilev) > rmdich) .AND. (ztt(ilev) > tmelt-90._wp)              &
!                                .AND. (ztt(ilev) < tmelt+60._wp)) THEN
!         ztv (ilev)  =  ztt(ilev)

!         ztt (ilev)  =  tv2t ( iob, job, ztv(ilev), zzz(ilev) )
!                        ====
!       ENDIF
!     ENDDO
!   ENDIF

    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)

! bogus observations for semi-idealised tests
      IF (      (lbogus)             .AND. (zpp(ilev) >= 66999._wp)            &
          .AND. (ztt(ilev) > rmdich) .AND. (zpp(ilev) <= 86000._wp))           &
        ztt (ilev)  =  MAX( 100._wp, ztt(ilev) + 1._wp*fbogus )

! fill ODR
! --------
      omlbdy (nmlob,klev,nbtt  ) = ztt   (ilev)
      omlbdy (nmlob,klev,nbtter) = zoberr(klev)
      momlbd (nmlob,klev,nbterr) = insert( momlbd(nmlob,klev,nbterr)           &
                                         , ioberr(klev), nvrt )
! note: no data base flag set, override flag is set in obs header, not here;
!       therefore only data base flag and level below surface flag are set
      momlbd (nmlob,klev,nbtflg) = insert( momlbd(nmlob,klev,nbtflg)           &
                                         , nflgx(klev), nvftbp )
    ENDDO

! --------
! humidity     (dew point as input, relative humidity is written to ODR)
! --------                          ( --> temperature is also required )

    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)
      nflgx  (klev)  =   0
      zoberr (klev)  =  c0
      ioberr (klev)  =   1

! set dew point to missing if temperature is missing
      IF ((ztt(ilev) <= rmdich) .AND. (ztd(ilev) > rmdich)) THEN
        neventd (neqflg,icma) = neventd (neqflg,icma) + iactr
        ztd    (ilev)  =  rmdi
        ioberr (klev)  =  0
! if dew point temperature missing or below absolute minimum
! (no flag in 'nflgx' set)
      ELSEIF ((ztd(ilev) <= rmdich) .OR. (ztd(ilev) <= b4w+c1)) THEN
        IF ((ntdex(kobtyp) >= 1) .AND. (lneedq(klev))                          &
                                 .AND. (zpp(ilev) >= pqmin))                   &
          neventd (netmis,icma) = neventd(netmis,icma) + iactr
        ztd    (ilev)  =  rmdi
        ioberr (klev)  =  0
      ELSE
! if level below surface
        IF ((kproclv(ilev) == -2) .OR. (kproclv(ilev) == -1)) THEN
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(4) )
          ioberr (klev)  =  0
        ENDIF
! if humidity blacklisted
        IF (kblk(klev,4) == 1) THEN
          IF (ioberr(klev) == 1)  neventd (neqflg,icma) =                      &
                                  neventd (neqflg,icma) + iactr
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(2) )
          ioberr (klev)  =  0
        ENDIF
! if humidity not in valid height range
        IF (zpp(ilev) < pqmin) THEN
          IF (ioberr(klev) == 1)  neventd (neq300,icma) =                      &
                                  neventd (neq300,icma) + iactr
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(4) )
          ioberr (klev)  =  0
        ENDIF
! if dew point gross error
        IF (      (ztd(ilev) < tmelt-150._wp)                                  &
            .OR.  (ztd(ilev) > tmelt+ 40._wp)                                  &
            .OR. ((ztd(ilev) < tmelt- 90._wp) .AND. (ilev == ilevsfc))) THEN
          IF (ioberr(klev) == 1)  neventd (neqlow,icma) =                      &
                                  neventd (neqlow,icma) + iactr
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(3) )
!         ztd    (ilev)  =  rmdi
          ioberr (klev)  =  0
        ENDIF
! set dew point to passive if temperature is not active
        IF (.NOT. BTEST( momlbd(nmlob,klev,nbterr),nvrt )) THEN
          IF (ioberr(klev) == 1)  neventd (neqflg,icma) =                      &
                                  neventd (neqflg,icma) + iactr
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(5) )
          ioberr (klev)  =  0
        ENDIF
      ENDIF
      zttl (klev)  =  ztt(ilev)
      ztdl (klev)  =  ztd(ilev)
    ENDDO

! compute model compatible relative humidity
! ------------------------------------------

    IF (nlvp > 0) THEN

      CALL obs_td2rh ( madj_hum, nlvp, zttl , ztdl                             &
                     , rmdi, b1, b2w, b2i, b3, b4w, b4i, tmelt , zrhw, zrh )
!     ==============
    ENDIF

    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)
! gross error if relative humidity > 102 % (previous version: if > 120 %)
      IF (zrhw(klev) >  rtshlm) THEN
        IF (ioberr(klev) == 1)  neventd (neqbig,icma) =                        &
                                neventd (neqbig,icma) + iactr
        nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(3) )
        ioberr (klev)  =  0

! bias correction of relative humidity (near saturation only)
! ------------------------------------
      ELSEIF ((zrh(klev) >= rhtsat) .AND. (zrh(klev) < c1-epsy)) THEN
        IF ((zttl(klev) <  tmelt) .AND. (ioberr(klev) == 1))                   &
          neventd (neqsam,icma) = neventd(neqsam,icma) + iactr
        IF ((zttl(klev) >= tmelt) .AND. (ioberr(klev) == 1))                   &
          neventd (neqsap,icma) = neventd(neqsap,icma) + iactr
!       nflgx (klev)  =  IBSET( nflgx(klev), nvfbps(5) )
        nflgx (klev)  =  IBSET( nflgx(klev), nvfbps(6) )
      ELSEIF (zrh(klev) > c1+epsy) THEN
        IF ((zttl(klev) <  tmelt) .AND. (ioberr(klev) == 1))                   &
          neventd (neqclm,icma) = neventd(neqclm,icma) + iactr
        IF ((zttl(klev) >= tmelt) .AND. (ioberr(klev) == 1))                   &
          neventd (neqclp,icma) = neventd(neqclp,icma) + iactr
!       nflgx (klev)  =  IBSET( nflgx(klev), nvfbps(6) )
      ENDIF
      IF (zrh(klev) >= rhtsat) THEN
        zrh   (klev)  =  c1
        ztdl  (klev)  =  zttl(klev)
      ENDIF

! bogus observations for semi-idealised tests
      IF (      (lbogus)             .AND. (zpp(ilev) >= 66999._wp)            &
          .AND. (zrh(klev) > rmdich) .AND. (zpp(ilev) <= 86000._wp))           &
        zrh(klev) = MAX( 0.01_wp , MIN( c1, zrh(klev) +0.1_wp*fbogus ) )

! fill ODR
! --------
      omlbdy (nmlob,klev,nbtrh ) = zrh   (klev)
      omlbdy (nmlob,klev,nbtqer) = zoberr(klev)
      momlbd (nmlob,klev,nbterr) = insert( momlbd(nmlob,klev,nbterr)           &
                                         , ioberr(klev), nvrq )
      IF ((zrh(klev) > rmdich) .AND. (zrhw(klev) > rmdich)) THEN
        omlbdy (nmlob,klev,nbtdrh) = zrh (klev)  -  zrhw(klev)
      ELSEIF (zrh(klev) > rmdich) THEN
        omlbdy (nmlob,klev,nbtdrh) = c0
      ENDIF
! note: no data base flag set, override flag is set in obs header, not here;
!       therefore only data base flag and level below surface flag are set
      momlbd (nmlob,klev,nbtflg) = insert( momlbd(nmlob,klev,nbtflg)           &
                                         , nflgx(klev), nvfqbp )
    ENDDO

! ---------------
! horizontal wind
! ---------------

    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)
      nflgx  (klev)  =   0
      zoberr (klev)  =  c0
      ioberr (klev)  =   1
      nflag          =   0

! vertical velocity
      IF (      ((kcdftyp == ncdf_wprof) .OR. (kcdftyp == ncdf_radar_vad))     &
          .AND. (zww(ilev) > rmdich)) THEN
!   (bad reporting practice: zero wind probably means missing value)
        IF (ABS( zww(ilev) ) < epsy)  zww (ilev) = rmdi
      ENDIF

! if wind speed or direction missing (no flag in 'nflgx' set)
      IF ((zff(ilev) <= rmdich) .OR. (zdd(ilev) <= rmdich)) THEN
        IF ((nuex(kobtyp) >= 1) .AND. (lneedv(klev))) THEN
          IF (zff(ilev) <= rmdich) neventd (nefmis,icma) =                     &
                                   neventd (nefmis,icma) + iactr
          IF (zdd(ilev) <= rmdich) neventd (nedmis,icma) =                     &
                                   neventd (nedmis,icma) + iactr
        ENDIF
        zff    (ilev)  =  rmdi
        zdd    (ilev)  =  rmdi
        ioberr (klev)  =  0
      ELSE
! if level below surface
        IF ((kproclv(ilev) == -2) .OR. (kproclv(ilev) == -1)) THEN
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(4) )
          ioberr (klev)  =  0
        ENDIF
! if wind blacklisted
        IF (kblk(klev,1) == 1) THEN
          IF (ioberr(klev) == 1)  neventd (nedflg,icma) =                      &
                                  neventd (nedflg,icma) + iactr
          IF (ioberr(klev) == 1)  neventd (nefflg,icma) =                      &
                                  neventd (nefflg,icma) + iactr
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(2) )
          ioberr (klev)  =  0
        ENDIF
! if wind speed / direction gross error
        IF (      (ABS( zdd(ilev) ) > 360._wp)                                 &
            .OR.  (zff(ilev) > 150._wp)                                        &
            .OR. ((zff(ilev) >  90._wp) .AND. (zpp(ilev) < 70000._wp))         &
            .OR.  (zff(ilev) < -epsy)) THEN
          IF (ioberr(klev) == 1)  neventd (nefneg,icma) =                      &
                                  neventd (nefneg,icma) + iactr
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(3) )
          ioberr (klev)  =  0
          zff    (ilev)  =  MAX( zff(ilev) , c0 )
        ENDIF
        IF (zff(ilev) < epsy)  zdd (ilev)  =  c0
! if VAD wind speed bad reporting practice
        IF ((zff(ilev) < fflim_vad) .AND. (kcdftyp == ncdf_radar_vad)) THEN
          IF (ioberr(klev) == 1)  neventd (nefneg,icma) =                      &
                                  neventd (nefneg,icma) + iactr
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(5) )
          ioberr (klev)  =  0
        ENDIF 
! aircraft roll angle as revised quality flag for wind
        IF (lamdarml) THEN
!         IF (iroll(ilev) == 1) THEN
!           nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(5) )
!           CALL MVBITS( iroll(ilev) , 0 , nvfboc(5) , nflgx(klev) , nvfbps(5) )
!           zoberr (klev)  =  rmdi
!         ENDIF
          IF ((iroll(ilev) < 0) .OR. (iroll(ilev) > nibits(nvfboc(1))))        &
            iroll (ilev) = nibits(nvfboc(1))
          CALL MVBITS( iroll(ilev) , 0 , nvfboc(1) , nflgx(klev) , nvfbps(1) )
          IF (iroll(ilev) == 1)  ioberr (klev)  =  0
! aircraft roll angle in station characteristics header word ('nvaaoc' bits)
          IF (iroll(ilev) == nibits(nvfboc(1)))  iroll (ilev) = nibits(nvaaoc)
          CALL MVBITS( iroll(ilev), 0 , nvaaoc , momlhd(nmlob,nhschr), nvaabp )
        ENDIF
        IF (zstdff(ilev) > 0.1_wp) THEN
! wind profiler / radar VAD standard deviation wind speed defined
          IF (zstdff(ilev) > 10._wp) THEN
            IF (ioberr(klev) == 1)  neventd (nefflg,icma) =                    &
                                    neventd (nefflg,icma) + iactr
            nflag          =  1
            ioberr (klev)  =  0
          ENDIF
!   adjust observation error
          zoberr (klev)  =  zstdff(ilev)
        ENDIF
      ENDIF
      IF (      ((kcdftyp == ncdf_wprof) .OR. (kcdftyp == ncdf_radar_vad))     &
          .AND. ((zff(ilev) > rmdich) .OR. (zww(ilev) > rmdich))) THEN
! wind profiler / radar VAD quality information (also valid for vertical wind)
        IF ((iqci (ilev) < 0) .OR. (iqci(ilev) > nibits(nvfboc(1))))           &
          iqci (ilev) = nibits(nvfboc(1))
        IF (nflag /= 1)  nflag  =  iqci(ilev)
        IF (iqci (ilev) == 1) THEN
          IF (ioberr(klev) == 1)  neventd (nefflg,icma) =                      &
                                  neventd (nefflg,icma) + iactr
          ioberr (klev)  =  0
        ENDIF
! wind profiler (radar VAD) signal to noise ratio (also valid for vertical wind)
        IF ((zsinor(ilev) > rmdich) .AND. (zsinor(ilev) < -20.0_wp)) THEN
          IF (ioberr(klev) == 1)  neventd (nefflg,icma) =                      &
                                  neventd (nefflg,icma) + iactr
          nflag          =  1
          ioberr (klev)  =  0
        ENDIF
      ENDIF
      CALL MVBITS( nflag , 0 , nvfboc(1) , nflgx(klev) , nvfbps(1) )
    ENDDO

! Transformation of wind speed and direction to wind components
! in the rotated model coordinate system
! -------------------------------------------------------------
    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)
      IF (zff(ilev) > rmdich) THEN
        zuu (klev,1)  =  zff(ilev) * (-SIN( zdd(ilev)*r_degrad ))
        zvv (klev,1)  =  zff(ilev) * (-COS( zdd(ilev)*r_degrad ))
      ELSE 
        zuu (klev,1)  =  c0
        zvv (klev,1)  =  c0
      ENDIF
      zrlat (klev,1)  =  omlhed(nmlob,nhjlat)
      zrlon (klev,1)  =  omlhed(nmlob,nhilon)
      IF ((zdlat(ilev) > rmdich) .AND. (lmult_shft))                           &
        zrlat (klev,1)  =  zrlat(klev,1) + zdlat(ilev)
      IF ((zdlon(ilev) > rmdich) .AND. (lmult_shft))                           &
        zrlon (klev,1)  =  zrlon(klev,1) + zdlon(ilev)
    ENDDO
 
    IF (nlvp > 0) THEN

      CALL uv2uvrot_vec ( zuu, zvv, zrlat, zrlon, r_pollat, r_pollon, nlvp, 1 )
!     =================
    ENDIF

    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)

! bogus observations for semi-idealised tests
      IF (lbogus) THEN
        IF ((zpp(ilev) > 55999._wp) .AND. (zpp(ilev) < 86001._wp)) THEN
          zuu (klev,1)  =  zuu (klev,1) + 1._wp*fbogus
          zvv (klev,1)  =  zvv (klev,1) + 1._wp*fbogus
        ELSEIF (zpp(ilev) > 25999._wp .AND. zpp(ilev) < 54001._wp) THEN
          zuu (klev,1)  =  zuu (klev,1) + 1._wp*fbogus
        ELSEIF (zpp(ilev) < 24001._wp) THEN
          zuu (klev,1)  =  zuu (klev,1) - 1._wp*fbogus
          zvv (klev,1)  =  zvv (klev,1) - 1._wp*fbogus
        ENDIF
      ENDIF

! fill ODR
! --------
      IF (zff(ilev) <= rmdich)  zuu (klev,1)  =  rmdi
      IF (zff(ilev) <= rmdich)  zvv (klev,1)  =  rmdi
      omlbdy (nmlob,klev,nbtu  ) = zuu   (klev,1)
      omlbdy (nmlob,klev,nbtv  ) = zvv   (klev,1)
      omlbdy (nmlob,klev,nbtuer) = zoberr(klev)
      momlbd (nmlob,klev,nbterr) = insert( momlbd(nmlob,klev,nbterr)           &
                                         , ioberr(klev), nvru )
! note: no data base flag set, override flag is set in obs header, not here;
!       therefore only data base flag and level below surface flag are set
      momlbd (nmlob,klev,nbtflg) = insert( momlbd(nmlob,klev,nbtflg)           &
                                         , nflgx(klev), nvfubp )
      omlbdy (nmlob,klev,nbtw  ) = zww   (ilev)
      omlbdy (nmlob,klev,nbtsnr) = zsinor(ilev)
      IF (zstdff(ilev) > 0.1_wp)                                               &
        omlbdy (nmlob,klev,nbtuac) = zstdff(ilev)
    ENDDO

! ------------------
! observation errors
! ------------------
! get the observation errors only at the end for those observations 
! for which the errors need to be specified
! (i.e. for which the errors have been set to zero previously
! -----------------------------------------------------------------

    DO klev = 1 , nlvp

      CALL obs_find_level ( nerlev, rolnlv, omlbdy(nmlob,klev,nbtlop)          &
                          , ilverr, fiperr )
!     ===================

      IF (omlbdy(nmlob,klev,nbtz ) > rmdich) THEN
        omlbdy(nmlob,klev,nbtzer) =       fiperr * oezsond(ilverr)             &
                                    + (c1-fiperr)* oezsond(ilverr+1)
        IF (lamdarml)                                                          &
          momlbd (nmlob,klev,nbterr) = IBCLR( momlbd(nmlob,klev,nbterr), nvrz )
      ENDIF
      IF (omlbdy(nmlob,klev,nbtu ) > rmdich) THEN
        IF (lamdarml) THEN
          omlbdy(nmlob,klev,nbtuer) =       fiperr * oeairep(ilverr)           &
                                      + (c1-fiperr)* oeairep(ilverr+1)
        ELSEIF (omlbdy(nmlob,klev,nbtuer) < +epsy) THEN
          omlbdy(nmlob,klev,nbtuer) =       fiperr * oevsond(ilverr)           &
                                      + (c1-fiperr)* oevsond(ilverr+1)
        ELSE
!    standard deviation wind speed + 0.5 * radiosonde wind error
          omlbdy(nmlob,klev,nbtuer) =   omlbdy(nmlob,klev,nbtuer)              &
                                      +     fiperr * oevsond(ilverr)   *c05    &
                                      + (c1-fiperr)* oevsond(ilverr+1) *c05
        ENDIF
      ENDIF
      IF (omlbdy(nmlob,klev,nbtt ) > rmdich) THEN
        IF (lamdarml) THEN
          omlbdy(nmlob,klev,nbtter) =       fiperr * oetairp(ilverr)           &
                                      + (c1-fiperr)* oetairp(ilverr+1)
        ELSE
          omlbdy(nmlob,klev,nbtter) =       fiperr * oetsond(ilverr)           &
                                      + (c1-fiperr)* oetsond(ilverr+1)
        ENDIF
      ENDIF
      IF (omlbdy(nmlob,klev,nbtrh) > rmdich) THEN
                                                roberr  =  rherr1 * c100r
        IF (omlbdy(nmlob,klev,nbtrh) < rrhlim)  roberr  =  rherr2 * c100r
        IF (omlbdy(nmlob,klev,nbtt ) < rttlim)  roberr  =  rherr3 * c100r
        IF (klev == klevsfc)                    roberr  =  rherr1 * c100r *4._wp
        omlbdy(nmlob,klev,nbtqer) = roberr
      ENDIF
!     PRINT '("nbtqer ",A,2I3,2I5,F8.0,2F8.2,I5)'                              &
!           , yomlhd(nmlob), ilevsfc, klevsfc, klev, momlbd(nmlob,klev,nbtlid) &
!           , omlbdy(nmlob,klev,nbtp) , omlbdy(nmlob,klev,nbtqer)              &
!           , omlbdy(nmlob,klev,nbtrh), momlbd(nmlob,klev,nbterr)
    ENDDO

!-------------------------------------------------------------------------------
! Section 5: Create surface report from surface-level TEMP observations
!-------------------------------------------------------------------------------
!   PRINT *,'zqq1 ', yomlhd(nmlob), klevsfc, ilevsfc, ksurfob(irpl)            &
!                  , momlbd(nmlob,klevsfc,nbterr)

    IF ((ksurfob(irpl) >= 1) .AND. (ilevsfc > 0)) THEN
      nsgob = ksurfob(irpl)

! copy surface-level data from multi-level report into single-level ODR
! ---------------------------------------------------------------------

      osgbdy (nsgob,nbsp  ) = omlbdy(nmlob,klevsfc,nbtp  )
      osgbdy (nsgob,nbsu  ) = omlbdy(nmlob,klevsfc,nbtu  )
      osgbdy (nsgob,nbsv  ) = omlbdy(nmlob,klevsfc,nbtv  )
      osgbdy (nsgob,nbst  ) = omlbdy(nmlob,klevsfc,nbtt  )
      osgbdy (nsgob,nbsrh ) = omlbdy(nmlob,klevsfc,nbtrh )
      osgbdy (nsgob,nbsz  ) = omlbdy(nmlob,klevsfc,nbtz  )
      osgbdy (nsgob,nbsuer) = omlbdy(nmlob,klevsfc,nbtuer)
      osgbdy (nsgob,nbster) = omlbdy(nmlob,klevsfc,nbtter)
      osgbdy (nsgob,nbsqer) = omlbdy(nmlob,klevsfc,nbtqer)
      osgbdy (nsgob,nbszer) = omlbdy(nmlob,klevsfc,nbtzer)
      osgbdy (nsgob,nbsdrh) = omlbdy(nmlob,klevsfc,nbtdrh)
!     mosgbd (nsgob,nbslsg) = momlbd(nmlob,klevsfc,nbtlsg)
      mosgbd (nsgob,nbserr) = momlbd(nmlob,klevsfc,nbterr)
      mosgbd (nsgob,nbslid) = momlbd(nmlob,klevsfc,nbtlid)
      mosgbd (nsgob,nbsflg) = momlbd(nmlob,klevsfc,nbtflg)  ! ???
      mosgbd (nsgob,nbsqcf) = 0
! cancel (surface) pressure obs if converted from reported height
      IF (ABS( kproclv(ilev) ) /= 1) THEN
        osgbdy (nsgob,nbsp  ) = rmdi
        osgbdy (nsgob,nbsz  ) = rmdi
        osgbdy (nsgob,nbszer) = rmdi
        mosgbd (nsgob,nbserr) = IBCLR( mosgbd(nsgob,nbserr) , nvrz )
        CALL MVBITS( 0 , 0 , nvfaoc , mosgbd(nsgob,nbsflg) , nvfzbp )
      ENDIF
! height is generally not used from the derived surface-level report 
! because a surface pressure increment is derived from the multi-level report
      mosgbd (nsgob,nbserr) = IBCLR( mosgbd(nsgob,nbserr) , nvrz )
      IF (osgbdy(nsgob,nbsp  ) > rmdich) THEN
        mosgbd (nsgob,nbsflg) = IBSET( mosgbd(nsgob,nbsflg) , nvfzbp+nvfbps(5) )
        IF (osgbdy(nsgob,nbszer) < epsy)                                       &
        osgbdy (nsgob,nbszer) =   oezsond(1)                                   &
                                * osgbdy(nsgob,nbsp) *r_g /(r_d * tmelt)
      ENDIF
!     PRINT *,'klevsfc2 ', yomlhd(nmlob), klevsfc, omlbdy(nmlob,klevsfc,nbtp)

! set observation to passive to corresponding station selection criteria fail
! ---------------------------------------------------------------------------

! get observation height and difference to model orography
      IF (zzz(ilevsfc) > rmdich) THEN
        zzaltsy = zzz(ilevsfc)
      ELSE
        zzaltsy = osghed(nsgob,nhalt)
      ENDIF
!     fisd    = osghed(nsgob,nhsurf) - osghed(nsgob,nhalt)
      fisd    = osghed(nsgob,nhalt) - osghed(nsgob,nhsurf)
      fisdtt  = ((fdoro(3)-c1)/c2 + SIGN( (fdoro(3)+c1)/c2 , fisd )) * fisd
      fisdrh  = ((fdoro(4)-c1)/c2 + SIGN( (fdoro(4)+c1)/c2 , fisd )) * fisd
      fisduv  = ((fdoro(1)-c1)/c2 + SIGN( (fdoro(1)+c1)/c2 , fisd )) * fisd
      lseaobs = BTEST( mosghd(nsgob,nhschr),nvsebp )

! 2-m temperature
      IF (osgbdy(nsgob,nbst) > rmdich) THEN
        IF (lbogus) osgbdy(nsgob,nbst) = MAX( 100._wp,  osgbdy(nsgob,nbst)     &
                                                      - 1._wp*fbogus )
        lseaonly = (ABS(altopsu(3)) < epsy)
        IF (     ((.NOT. lseaonly) .AND. (zzaltsy > altopsu(3)))               &
            .OR. ((      lseaonly) .AND. (.NOT. lseaobs))                      &
            .OR. (fisdtt > doromx(3))) THEN
          IF (BTEST( mosgbd(nsgob,nbserr),nvrt ))                              &
            neventd(netalt,icma) = neventd(netalt,icma) + iactr
          mosgbd (nsgob,nbserr) = IBCLR( mosgbd(nsgob,nbserr) , nvrt )
          mosgbd (nsgob,nbsflg) = IBSET( mosgbd(nsgob,nbsflg), nvftbp+nvfbps(4))
        ENDIF
      ENDIF

! 2-m relative humidity
      IF (osgbdy(nsgob,nbsrh) > rmdich) THEN
        IF (lbogus) osgbdy(nsgob,nbsrh) = MAX( MIN( c1 ,  osgbdy(nsgob,nbsrh)  &
                                                        + 0.1_wp*fbogus )      &
                                             , 0.01_wp )
        lseaonly = (ABS(altopsu(4)) < epsy)
        IF (     ((.NOT. lseaonly) .AND. (zzaltsy > altopsu(4)))               &
            .OR. ((      lseaonly) .AND. (.NOT. lseaobs))                      &
            .OR. (fisdrh > doromx(4))) THEN
          IF (BTEST( mosgbd(nsgob,nbserr),nvrq ))                              &
            neventd(neqalt,icma) = neventd(neqalt,icma) + iactr
          mosgbd (nsgob,nbserr) = IBCLR( mosgbd(nsgob,nbserr) , nvrq )
          mosgbd (nsgob,nbsflg) = IBSET( mosgbd(nsgob,nbsflg), nvfqbp+nvfbps(4))
        ENDIF
      ENDIF

! 10-m horizontal wind
      IF (osgbdy(nsgob,nbsu) > rmdich) THEN
        IF (lbogus) osgbdy(nsgob,nbsu) = osgbdy(nsgob,nbsu) - 1._wp*fbogus
        IF (lbogus) osgbdy(nsgob,nbsv) = osgbdy(nsgob,nbsv) + 1._wp*fbogus
        lseaonly = (ABS(altopsu(1)) < epsy)
        IF (     ((.NOT. lseaonly) .AND. (zzaltsy > altopsu(1)))               &
            .OR. ((      lseaonly) .AND. (.NOT. lseaobs))                      &
            .OR. (fisduv > doromx(1))) THEN
          IF (BTEST( mosgbd(nsgob,nbserr),nvru ))                              &
            neventd(nevalt,icma) = neventd(nevalt,icma) + iactr
          mosgbd (nsgob,nbserr) = IBCLR( mosgbd(nsgob,nbserr) , nvru )
          mosgbd (nsgob,nbsflg) = IBSET( mosgbd(nsgob,nbsflg), nvfubp+nvfbps(4))
        ENDIF
      ENDIF

! set screen-level observations in multi-level report to passive
! (note: the obs operator for multi-level reports in COSMO
!        inter-/extra-polates between model levels and never
!        uses T2m etc.)
! --------------------------------------------------------------

      IF (omlbdy(nmlob,klevsfc,nbtt) > rmdich) THEN
        momlbd (nmlob,klevsfc,nbterr) = IBCLR( momlbd(nmlob,klevsfc,nbterr)    &
                                             , nvrt )
        momlbd (nmlob,klevsfc,nbtflg) = IBSET( momlbd(nmlob,klevsfc,nbtflg)    &
                                             , nvftbp+nvfbps(4) )
!       CALL MVBITS( 1 , 0 , nvfboc(4) , momlbd(nmlob,klevsfc,nbtflg)          &
!                                      , nvftbp+nvfbps(4) )
      ENDIF
      IF (omlbdy(nmlob,klevsfc,nbtrh) > rmdich) THEN
        momlbd (nmlob,klevsfc,nbterr) = IBCLR( momlbd(nmlob,klevsfc,nbterr)    &
                                             , nvrq )
        momlbd (nmlob,klevsfc,nbtflg) = IBSET( momlbd(nmlob,klevsfc,nbtflg)    &
                                             , nvfqbp+nvfbps(4) )
      ENDIF
      IF (omlbdy(nmlob,klevsfc,nbtu) > rmdich) THEN
        momlbd (nmlob,klevsfc,nbterr) = IBCLR( momlbd(nmlob,klevsfc,nbterr)    &
                                             , nvru )
        momlbd (nmlob,klevsfc,nbtflg) = IBSET( momlbd(nmlob,klevsfc,nbtflg)    &
                                             , nvfubp+nvfbps(4) )
      ENDIF

! avoid active status flags for empty reports
! -------------------------------------------

    ELSEIF (ksurfob(irpl) >= 1) THEN
      nsgob = ksurfob(irpl)
      mosghd (nsgob,nhpass) = -1
      mosgbd (nsgob,nbserr) =  0
    ENDIF

!-------------------------------------------------------------------------------
! Section 6: Discard levels which do not contain any data
!            except for pressure derived from height
!-------------------------------------------------------------------------------
! note that 'isrtlvp' must not be used any more after the following loop

    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)
      IF (      (omlbdy(nmlob,klev,nbtu  ) <= rmdich)                          &
          .AND. (omlbdy(nmlob,klev,nbtt  ) <= rmdich)                          &
          .AND. (omlbdy(nmlob,klev,nbtrh ) <= rmdich)                          &
          .AND. (omlbdy(nmlob,klev,nbtw  ) <= rmdich)                          &
          .AND. (     (omlbdy(nmlob,klev,nbtz  ) <= rmdich)                    &
                 .OR. (kproclv(ilev) <= 1)))   kproclv (ilev) = -9
    ENDDO

    nrejlv = 0
    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)
      IF (kproclv(ilev) <= -9) THEN
        nrejlv = nrejlv + 1
      ELSEIF (nrejlv >= 1) THEN
        omlbdy (nmlob,klev-nrejlv,:) = omlbdy(nmlob,klev,:)
        momlbd (nmlob,klev-nrejlv,:) = momlbd(nmlob,klev,:)
      ENDIF
    ENDDO
    nlvp = nlvp - nrejlv
    momlhd (nmlob,nhnlev) = nlvp

!-------------------------------------------------------------------------------
! Section 7: Check existence of data in multi-level report and surface report
!-------------------------------------------------------------------------------

! multi-level report
! ------------------

    iactx = -1
    ipasx = -1
    DO klev = 1 , nlvp
!     ilev = isrtlvp(klev)
      IF (BTEST( momlbd(nmlob,klev,nbterr),nvru )) iactx (1) = 1
      IF (BTEST( momlbd(nmlob,klev,nbterr),nvrz )) iactx (2) = 1
      IF (BTEST( momlbd(nmlob,klev,nbterr),nvrt )) iactx (3) = 1
      IF (BTEST( momlbd(nmlob,klev,nbterr),nvrq )) iactx (4) = 1
      IF (omlbdy(nmlob,klev,nbtu  ) > rmdich) ipasx (1) = 1
      IF (omlbdy(nmlob,klev,nbtz  ) > rmdich) ipasx (2) = 1
      IF (omlbdy(nmlob,klev,nbtt  ) > rmdich) ipasx (3) = 1
      IF (omlbdy(nmlob,klev,nbtrh ) > rmdich) ipasx (4) = 1
    ENDDO
    iactx (5)             = MAX( iactx(1), iactx(3), iactx(4) )
    ipasx (5)             = MAX( ipasx(1), ipasx(3), ipasx(4) )
    iactx (2)             = MAX( iactx(2), iactx(5) )
    ipasx (2)             = MAX( ipasx(2), ipasx(5) )
    momlhd (nmlob,nhuexi) = MAX( MIN( 0 , -ipasx(1) ) , iactx(1) )
    momlhd (nmlob,nhtexi) = MAX( MIN( 0 , -ipasx(3) ) , iactx(3) )
    momlhd (nmlob,nhqexi) = MAX( MIN( 0 , -ipasx(4) ) , iactx(4) )
    momlhd (nmlob,nhaexi) = MAX( MIN( 0 , -ipasx(5) ) , iactx(5) )
    nzaexi                = MAX( MIN( 0 , -ipasx(2) ) , iactx(2) )
! (nzaexi ==  1): active data present
! (nzaexi == -1): only passive data present
! (nzaexi ==  0): no data at all in total report
    IF (nzaexi <= 0)                                                           &
      neventr (nenoda,icma) = neventr(nenoda,icma) + iactr
    IF ((nzaexi == -1) .AND. (iactr == 1) .AND. (lverpas)) THEN
      momlhd (nmlob,nhpass)  =  2
      momlhd (nmlob,nhflag)  =  IBSET( momlhd(nmlob,nhflag) , FL_NO_OBS )
      cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
      cma(icma)%cnt_ps = cma(icma)%cnt_ps + 1

! flag reports which are to be discarded completely
    ELSEIF ((nzaexi == 0) .OR. ((nzaexi == -1) .AND. (.NOT. lverpas))) THEN
      momlhd (nmlob,nhpass)  = -1
      IF (iactr == 1) cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
      IF (iactr == 0) cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
      cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
    ENDIF
!   PRINT *,'noctrj_m1 ', icma, irpl, cma(icma)%cnt_ps, cma(icma)%cnt_rj       &
!                       , iactr, nzaexi, momlhd(nmlob,nhpass), nlvp

! surface report
! --------------

    IF (ksurfob(irpl) >= 1) THEN
      nsgob = ksurfob(irpl)
      iactx = -1
      ipasx = -1
      IF (BTEST( mosgbd(nsgob,nbserr),nvru )) iactx (1) = 1
      IF (BTEST( mosgbd(nsgob,nbserr),nvrt )) iactx (3) = 1
      IF (BTEST( mosgbd(nsgob,nbserr),nvrq )) iactx (4) = 1
      IF (osgbdy(nsgob,nbsu  ) > rmdich) ipasx (1) = 1
      IF (osgbdy(nsgob,nbst  ) > rmdich) ipasx (3) = 1
      IF (osgbdy(nsgob,nbsrh ) > rmdich) ipasx (4) = 1
      iactx (5)             = MAX( iactx(1), iactx(3), iactx(4) )
      ipasx (5)             = MAX( ipasx(1), ipasx(3), ipasx(4) )
      nzaexi                = MAX( MIN( 0 , -ipasx(5) ) , iactx(5) )
      IF (      (nzaexi == -1) .AND. (iactr == 1) .AND. (lverpas)              &
          .AND. (ilevsfc > 0)) THEN
        mosghd (nsgob,nhpass)  =  2
        mosghd (nsgob,nhflag)  =  IBSET( mosghd(nsgob,nhflag) , FL_NO_OBS )
      ENDIF

! flag reports which are to be discarded completely
      IF (     (nzaexi == 0) .OR. ((nzaexi == -1) .AND. (.NOT. lverpas))       &
          .OR. (ilevsfc == 0))   mosghd (nsgob,nhpass) = -1
    ENDIF

  ENDDO  ! huge loop over report
! ~~~~~

! this is done after obs_cdf_redundancy
! TODO : delete reports with (moxxhd(nxxob,nhpass) == -1)
! CALL obs_del_old_reports
! ========================

! clean up
! --------

  IF (nrepml >= 1) DEALLOCATE ( blk_loc , STAT=istat )

  DEALLOCATE ( kblk    , STAT=istat )
  DEALLOCATE ( nflgx   , STAT=istat )
  DEALLOCATE ( ioberr  , STAT=istat )
  DEALLOCATE ( ztv     , STAT=istat )
  DEALLOCATE ( zttl    , STAT=istat )
  DEALLOCATE ( ztdl    , STAT=istat )
  DEALLOCATE ( zrhw    , STAT=istat )
  DEALLOCATE ( zrh     , STAT=istat )
  DEALLOCATE ( zuu     , STAT=istat )
  DEALLOCATE ( zvv     , STAT=istat )
  DEALLOCATE ( zrlat   , STAT=istat )
  DEALLOCATE ( zrlon   , STAT=istat )
  DEALLOCATE ( zoberr  , STAT=istat )
  DEALLOCATE ( lneedt  , STAT=istat )
  DEALLOCATE ( lneedq  , STAT=istat )
  DEALLOCATE ( lneedv  , STAT=istat )

  DEALLOCATE ( zpp     , STAT=istat )
  DEALLOCATE ( zdlat   , STAT=istat )
  DEALLOCATE ( zdlon   , STAT=istat )
  DEALLOCATE ( ztt     , STAT=istat )
  DEALLOCATE ( ztd     , STAT=istat )
  DEALLOCATE ( zff     , STAT=istat )
  DEALLOCATE ( zzz     , STAT=istat )
  DEALLOCATE ( zdd     , STAT=istat )
  DEALLOCATE ( izdt    , STAT=istat )
  DEALLOCATE ( izlv    , STAT=istat )
  DEALLOCATE ( iroll   , STAT=istat )
  DEALLOCATE ( kproclv , STAT=istat )
  DEALLOCATE ( isrtlvp , STAT=istat )

  DEALLOCATE ( ksurfob , STAT=istat )


!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_store_multilev
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_store_multilev

!===============================================================================

ELEMENTAL INTEGER FUNCTION insert  ( invar, inval, ibp )
  !-----------------------------------------------------
  INTEGER (KIND=iintegers)  , INTENT (IN)  ::  invar, inval, ibp
  !----------------------------------------------------------------------------
  ! inserts bits set in 'inval' into integer word 'invar' at bit position 'ibp'
  !----------------------------------------------------------------------------
  !
  insert = IOR( invar , ISHFT( inval, ibp ) )
  !
END FUNCTION insert

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

END MODULE src_obs_cdfin_mult
