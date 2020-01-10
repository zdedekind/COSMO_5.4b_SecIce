!+ Data module for the variables only used to compute the local information
!-------------------------------------------------------------------------------

MODULE data_obs_record

!-------------------------------------------------------------------------------
!
! Description:
!   This module contains the related information on the Observation Data Record
!   (ODR) which is the long-term storage of the observations obtained after
!   reading and pre-processing of the observation reports. As such, the ODR is
!   the interface to those modules which process the observational information
!   further, in particular in order to apply forward observation operators to
!   derive simulated observations resp. observation increments for data
!   assimiation (e.g. nudging).
!   This module also contains the simulated observations in the form of the
!   Simulated Observation Record (SOR). Within COSMO, the SOR is purely for
!   writing to (NetCDF and ASCII) feedobs files (for verification and LETKF).
!
!   Specifically, this module contains:
!    - the format of the observation data records (ODR) and related variables
!    - the ODR themselves
!    - the SOR and its format
!    - the masking constants for packed information
!
!   Note: This module is part of the 'COSMO data assimilation libraries'
!         (library 1 for reading data from NetCDF observation input files,
!          library 2 for observation operators including quality control).
!         It is used commonly by the COSMO model and the 3DVAR program package !
!
! Current Code Owner (for COSMO and for DWD 3DVAR):
!  DWD, Christoph Schraff
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
!  Removal of any integers exceeding 16-bit integer representation.
! 1.13       1998/10/22 Christoph Schraff
!  Switch for observation processing moved to 'data_nudge_all'.
! 1.19       1998/12/11 Christoph Schraff
!  Revised format of the ODR: collective indices and 2. part of flags removed,
!  threshold quality control flags (including verification status) introduced.
! 1.27       1999/03/29 Christoph Schraff
!  VOF format and surface analysis limits and allocatable arrays introduced.
!  Optional groups introduced in ODR format, ODR flag formats modified, and
!  6 bit hollerith station identity removed.
! 1.28       1999/04/19 Christoph Schraff
!  'rraint(1)' set to 1 (hourly precipitation observations can be verified).
! 1.31       1999/07/01 Christoph Schraff
!  Buffer for station id's and modified string length for surface analysis.
! 1.36       2000/02/24 Christoph Schraff
!  Observation increment format of VOF included.
! 1.38       2000/04/06 Christoph Schraff
!  Index for pressure tendency added to ODR and VOF. In-line docu corrected.
!  Station characteristics extended (e.g. by redundant report flag).
! 1.40       2000/05/23 Christoph Schraff
!  'nvbqcf' extended by threshold quality control flag bit for geopotential.
! 2.5        2001/06/01 Christoph Schraff
!  Introduction of the flight track check and flight thinning flags.
!  Special value of instrument specification indicator for ship observations.
! 2.13       2002/01/18 Christoph Schraff
!  Adaptation for 32-bit machines: Re-organisation of station characterisics
!  and main flag word to reduce the number of used bits from 32 to less than 32.
!  ODR combined cloud and weather group word changed from real to integer.
!  Data statements replaced by direct assignment, usually linked with adding
!  the parameter attribute to the corresponding variables.
! 2.19       2002/10/24 Michael Buchhold
!  Change quality control threshold values for surface analysis
! 3.3        2003/04/22 Maria Tomassini + Christoph Schraff
!  Extension of ODR for GPS reports.
!  Flag for Vaisala RS80 humidity bias correction in station characteristics.
! 3.12       2004/09/15 Christoph Schraff
!  Extension to (prepare to) include assimilation of satellite retrievals.
!  ODR header element 'nhadif' (height diff) replaced by 'nhsurf' (orography).
! V4_5         2008/09/10 Christoph Schraff
!  Adaptions to observation from NetCDF files (extension / modification of
!  ODR formats, deletion of o??vip formats, small modfications of VOF format
!  ('nvhflg', 'nvhpas', 'nvbmfw')).
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_22        2012/01/31 Christoph Schraff
!  - Length of station ID increased from 8 to 9.
!  - ODR body for GPS reports complemented.
!  - High cloud cover and general cloud group word introduced in ODR.
!  - ODR header extended by 'nhsolz' for solar zenith angle, and by 'nbtddb'
!    for data base decoding time.
!  - Bit pattern for 'nhflag' adjusted to match revised feedback file table.
!  - Number of bits for aircraft flight phase and roll angle reduced.
!  - ODR body extended by 'nbtuac' for accuracy of wind obs.
!  - General clean-up of the ODR description, information on packed information.
!    moved from VOF to ODR description. VOF description moved to module
!    'src_obs_print_vof', surface analysis limits and allocatable arrays for
!    surface analysis moved to module 'data_nudge_all'.
! V4_28        2013/07/12 Christoph Schraff
!  - Inclusion of the Simulated Observation Record (SOR) from 'data_nudge_local;
!    increments replaced by simulated observations in SOR (except 'dmlhed').
!  - Additional elements 'nso_cm', 'nso_ch', 'nso_cbs', 'nso_zpd' in SOR.
!  - Additional element 'nhqcps' in ODR body. Additional status flag 'nvrcbs'.
!  - Temporary flags 'nvrzbc' and 'nvrqbc' for LBC QC checks introduced.
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
! V5_3         2015-10-09 Christoph Schraff
!  Additional entries added in single-level ODR and SOR for verification
!  (model equivalent calculator) purposes: wind speed + direction, dewpoint
!  temperature, radiation, 3-hourly precip and wind gusts. Max. and min. T-2m
!  over arbitrary periods instead of fixed 12 hours.
! V5_4         2016-03-10 Christoph Schraff
!  Removal of variables related to the AOF interface: 'nhqofl' and 'nbscfw' with
!  related bit patterns, 'nbspr', 'nbsper', 'nbscer', 'nbstg', 'nbsfgu', 'nbsfme'.
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
! Modules used:
!
!-------------------------------------------------------------------------------

USE data_parameters , ONLY :   &
    wp,         & ! kind-type parameter for "normal" real    variables
    iintegers     ! kind-type parameter for "normal" integer variables

!-------------------------------------------------------------------------------

IMPLICIT NONE

!===============================================================================

! Local Declarations:

!-------------------------------------------------------------------------------
! Global (i.e. public) Declarations:
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!  Section 1 : Format of the observation data records (ODR)
!-------------------------------------------------------------------------------


!       ------------------------------------------------------------------------
!       1.1    ODR header format (words)
!       ------------------------------------------------------------------------

!       1.1.1  Header formats of ODR reports:'omlhed','osghed','ogphed','otvhed'
!              -----------------------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    mxrhed = 18  ,& ! header length of multi-level reports
    mxshed = 16  ,& ! header length of single-level reports
    mxghed = 10  ,& ! header length of GPS reports
    mxthed = 13  ,& ! header length of satellite retrieval reports
    nhilon =  1  ,& ! longitude of observing station
    nhjlat =  2  ,& ! latitude  of observing station
    nhalt  =  3  ,& ! station altitude [m]
    nhtime =  4  ,& ! (exact) time of observation in forecast hours
    nhsurf =  5  ,& ! height of model grid pt. to which obs. is assigned
    nhzio  =  6  ,& ! longitude of obs. station (or lowest obs.) in grid pt unit
    nhzjo  =  7  ,& ! latitude  of obs. station in grid pt. units
    nhsynt =  8  ,& ! nominal (synoptic) time of observation in forecast hours
    nhtddb =  9  ,& ! data base decoding time in forecast hours
    nhsolz = 10  ,& ! solar zenith angle [deg]
    nhvcbu = 11  ,& ! correction factor to vertical correlation scale for wind
                    ! at base of report
    nhvcbt = 12  ,& ! as 'nhvcbu', but for temperature
    nhvcbq = 13  ,& ! as 'nhvcbu', but for humidity
    nhvctu = 14  ,& ! correction factor to vertical correlation scale for wind
                    ! at top of report
    nhvctt = 15  ,& ! as 'nhvctu', but for temperature
    nhvctq = 16  ,& ! as 'nhvctu', but for humidity
                    !   the following elements are (re-)determined in the
                    !   observation increment part (in src_sing_local.f90)
                    !   whenever the quality control is re-done (rather than
                    !   once when reading in the observation processing -
                    !   however this info is used over several timesteps
                    !   and therefore part of the ODR array, where it is
                    !   put in the header because it is single-level info)
    nhtvip = 17  ,& ! observed multi-level pressure interpolated to the lowest
                    ! model level
    nhtviz = 18  ,& ! vertical distance to nearest observation
    nh1wta = 11  ,& ! 1dvar retriev.: influence radius of temporal weight: past
    nh1wte = 12  ,& ! 1dvar retriev.: influence radius of temporal w.: future
    nh1tip = 13     ! 1dvar retriev.: temporal interval for re-doing minimizat.

!       1.1.2  Header formats of ODR reports:'momlhd','mosghd','mopghd','motvhd'
!              -----------------------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    mxrhdf = 30  ,& ! header length of multi-level reports
    mxshdf = 20  ,& ! header length of single-level reports
    mxghdf = 19  ,& ! header length of GPS reports
    mxthdf = 25  ,& ! header length of satellite retrieval reports
    nhio   =  1  ,& ! (local) x-coord. of grid pt. assigned to obs
    nhjo   =  2  ,& ! (local) y-coord. of grid pt. assigned to obs
    nhitot =  3  ,& ! global x-coord. of grid pt. assigned to obs
    nhjtot =  4  ,& ! global y-coord. of grid pt. assigned to obs
    nhobtp =  5  ,& ! observation type
    nhcode =  6  ,& ! code type
    nhschr =  7  ,& ! station characteristics                          (see 1.2)
    nhflag =  8  ,& ! report flags (obs type, surf., alt., sta ID)
                    !   (bit pattern, as in NetCDF feedobs/feedback file
                    !    see table 'flags' in 'mo_fdbk_tables.f90')
                    !   (only used for NetCDF- (or feedobs file) reading;
                    !    note that the order of setting flags 'surf',
                    !    'area', 'height', and 'blacklist' in COSMO
                    !    is not exactly as indicated in 'flags')
!   nhqofl =  8  ,& ! report flags (rds) on lat/long/date/time/alt     (see 1.2)
                    !   (only used for AOF-reading)
    nhpass =  9  ,& ! flag for report being set to 'passive'           (see 1.2)
    nhqcfw = 10  ,& ! status of threshold quality control (QC) and of writing to
                    ! feedobs files (VOF or NetCDF), and QC flags for surface
                    ! pressure increments from multi-level reports:
                    !   >= 0 : QC flags and increments are computed and written
                    !          to feedobs files at current timestep    (see 1.2)
                    !   = -1 : QC flags not set yet (no QC done yet)
                    !   =x-2 : QC flags also set for pressure, \  flags not yet
                    !   =x-4 : QC flags also set for wind/T/q,  > written to
                    !   (x=-1 or -5 resp. -3)                  /  feedobs files
                    !   =-99 : report incl. QC flags already written to f.-files
    nhcorr = 11  ,& ! update sequence number (station correction indicator)
    nhcat  = 12  ,& ! data     category (from BUFR Section 1)
    nhcats = 13  ,& ! data sub-category (from BUFR Section 1)
    nhkz   = 14  ,& ! DWD internal classification number (observation type)
    nhcent = 15  ,& ! originating centre  +  (1000* sub-centre)
    nhstid = 16  ,& ! station identity number / satellite ID (WMO Table C-5)
    nhdate = 17  ,& ! absolute exact observation date [yyyymmdd]
    nhhrmn = 18  ,& ! absolute exact observation time [hhmm]
    nhsyhr = 19  ,& ! absolute nominal (synoptic) observation time [yymmddhh]
    nhstyp = 20  ,& ! sing-lv obs: station type (buoy: MQOBL, BUFR Table 002149,
                    !                            else: NIX  , BUFR Table 002001)
    nhrtyp = 20  ,& ! radiosonde type    (NRARA, see WMO common code Table C2)
    nhnlev = 21  ,& ! number of obs. levels (for multi-level reports)
    nhvqcf = 22  ,& ! for satellite retrieval: threshold quality control flags
                    !                            (as in ODR body, word 'nb?qcf')
    nhaexi = 22  ,& ! for conventional multi-level report:
                    !   flag for existence of wind or temperature
                    ! flags 'nhuexi', 'nhtexi', 'nhqexi', 'nhaexi':
                    !   =  1 : active  observations exist
                    !   =  0 : no      observations exist
                    !   = -1 : passive observations only
    nhuexi = 23  ,& ! flag for existence of wind data      in multi-level report
    nhtexi = 24  ,& ! flag for existence of temperature    in multi-level report
    nhqexi = 25  ,& ! flag for existence of humidity data  in multi-level report
    nhtrac = 26  ,& ! tracking technique (NSASA, see WMO common code Table C7)
    nhrad  = 27  ,& ! solar and IR radiation correction (NSR, BUFR Table 002013)
    nhna4  = 28  ,& ! instrument type                   (NA4, BUFR Table 002003)
    nhwce  = 29  ,& ! wind comput. enhancement (w-prof, MWCE, BUFR Table 025021)
    nhdt   = 30  ,& ! time period of measurement (e.g. w-prof)               [s]
    nhqcps = 30     ! temporary entry for radiosonde only (QC flag for pressure)

!       1.1.3  Header formats of ODR reports:'yomlhd','yosghd','yopghd','yotvhd'
!              -----------------------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    ilstid  =  9 ,& ! character length of the station identity
    ilstidp =  9    ! character length used for printing the station ID
                    ! Note: (ilstid >= ilstidg >= ilstidp), cf.data_nudge_gather


!       ------------------------------------------------------------------------
!       1.2    Bit patterns for packed information in ODR (and VOF) header
!       ------------------------------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
!   variable          meaning                                           word no.
!   --------          -------                                           --------
!   --------------- ! bit pos. for flags on   (VOF: nvhsch) ----------- nhschr
                    ! note: all entries in 'nhschr' are used only for
                    !       output in VOF except: 'nvsebp', 'nvscbp',
                    !                             'nvapbp'.
                    ! note: 'nvpabp','nvpsbp' relate to the report status
                    !       (see word 'nvhpas'):
    nvpabp =  0  ,& ! bit pos. for report set passive since it is          "
                    !              used in a multi-level pseudo report
    nvpsbp =  1  ,& ! bit pos. for report set passive since at least       "
                    !              1 of the following 5 flags or
                    !              the flight track flag applies
    nvobbp =  2  ,& ! bit pos. for flag: 'station location out of          "
                    !                     user-specified area'
    nvalbp =  3  ,& ! bit pos. for flag: 'distance model orography -       "
                    !                     station altitude too large'
    nvbkbp =  4  ,& ! bit pos. for flag: 'blacklisted station (ship/AMDAR)'"
    nvexbp =  5  ,& ! bit pos. for flag: 'observation or code type         "
                    !                     excluded at station location'
    nvrdbp =  6  ,& ! bit pos. for flag: 'redundant report'                "
    nvsebp =  7  ,& ! bit pos. for report located at sea grid pt.          "
    nvscbp =  8  ,& ! bit pos. for station correction indicator            "
    nvssbp =  9  ,& ! bit pos. for station suspicion indicator             "
!   nvsibp = 10  ,& ! bit pos. for important station indicator             "
    nvinbp = 13  ,& ! bit pos. for instrument specification indicator      "
                    !   - code for reading obs from AOF:
                    !        0 : missing
                    !        1 : other reports
                    !        2 : surface wind speed from anemometer
                    !        4 : AIREP = ASDAR obtained by ESA
                    !        5 : AMDAR
                    !        6 : SHIP, wind not from anemometer
                    !        7 : SHIP, wind speed from anemometer
                    !       18 : TEMP Vaisala RS80, bias-corrected RH
                    !   - code for reading GPS obs from COST-716 ASCII:
                    !       36 : IWV reported
                    !       38 : IWV not reported, only ZTD reported
                    !   - for reading obs from NetCDF files: not used
!   nvinoc =  7  ,& !   no. of bits occ. by instrument specification       "
    nvhtbp = 20  ,& ! bit pos. for flight track error flag                 "
    nvhtoc =  1  ,& ! no. of bits occ. by flight track error flag          "
    nvhhbp = 21  ,& ! bit pos. for flight thinning flag                    "
    nvhhoc =  1  ,& ! no. of bits occ. by flight thinning flag             "
    nvapbp = 22  ,& ! bit pos. for phase of flight (WMO Table 008004)      "
    nvapoc =  4  ,& ! no. of bits occ. by extended phase of flight:        "
                    !   bits 0-3 : WMO Table 008004 (phase of flight)
                    !   bit  4   : origin of phase of flight
                    !             = 0 : original BUFR report,
                    !             = 1 : flight tracking in assim. code
    nvaabp = 26  ,& ! bit pos. for aircraft roll angle (WMO 002064)        "
    nvaaoc =  2     ! no. of bits occ. by aircraft roll angle              "
!   --------------- ! bit pos. for flags on   (VOF: nvhqfl) ----------- nhqofl
                    !   latitude  longitude  date    time  altitude        "
!   nvhbps(5) = (/          0    ,    6    ,  12   ,  18  ,   24   /)   ,      &
!   nvhboc    = 6 & ! no. of bits occ. by each flag                        "
!                   ! inner bit structure (pos./ no.) for each flag        "
                    !   hMsubstit QCsubstit override flag QC/hM-flag       "
!   nvhibp(5) = (/          0    ,    1    ,   2   ,   3   ,   5   /)   ,      &
!   nvhioc(5) = (/          1    ,    1    ,   1   ,   2   ,   1   /)
!   --------------- ! report status code:     (VOF: nvhpas) ----------- nhpass
                    !   0 - active   report (used by nudging)
                    !   1 - merged   report (i.e. set passive e.g. as
                    !                        a single-level aircraft
                    !                        report, but used in a
                    !                        merged multi-level rep.)
                    !   2 - passive or rejected report
                    !   --> not yet implemented:
                    !       2 - passive  report (not used by nudging, due
                    !                            to obs / code type
                    !                            or thinning)
                    !       3 - rejected report (not used by nudging, due
                    !                            to insufficient quality)
!   --------------- ! code table on           (VOF: nvhqcp) --------- nhqcfw
                    ! QC flag for surface pressure increment from multi-lev rep:
                    !   0 - active data used, value ok
                    !   1 - active data used, value not ok
                    !   2 - no active data usable, passive data used, value ok
                    !   3 - no active data usable, passive data used, not ok
                    !   4 - no data at all usable for extrapolation


!       ------------------------------------------------------------------------
!       1.3    ODR body format (words)
!       ------------------------------------------------------------------------

!       1.3.0  Number of levels in multi-level ODR 'omlbdy', 'momlbd', 'otvbdy'
!              ----------------------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    maxarl = 40     ! max. number of levels in multi-level aircraft reports

  INTEGER (KIND=iintegers)              :: &
!   maxrsl must be >= maxarl, maxrtv !
!   maxrsl is set equal to maxmlv (see data_nudge_all.f90)
    maxrsl = 100 ,& ! max. number of levels allowed for multi-level reports
    maxrtv = 50     ! max. number of levels in satellite retrieval  reports
                    !   ('maxrtv' must be >= 'ke' , ==> reset, not a parameter)

!       1.3.1  Body format of ODR of multi-level reports: 'omlbdy'
!              ---------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    mxrbdy = 18  ,& ! body length of multi-level reports
    nbtu   =  1  ,& ! u wind component [m/s]
    nbtv   =  2  ,& ! v wind component [m/s]
    nbtt   =  3  ,& ! temperature [K]
    nbtrh  =  4  ,& ! relative humidity [/]   (bias corrected)
    nbtp   =  5  ,& ! pressure [Pa]
    nbtz   =  6  ,& ! height [m]
    nbtuer =  7  ,& ! error of observed wind component
    nbtter =  8  ,& ! error of observed temperature
    nbtqer =  9  ,& ! error of observed rel. humidity
    nbtzer = 10  ,& ! error of observed height
                    ! Note: nbt?er are set to the negative rms errors, if the
                    ! observations have not passed the threshold quality control
    nbtzio = 11  ,& ! longitude in grid pt. units
    nbtzjo = 12  ,& ! latitude  in grid pt. units
    nbttim = 13  ,& ! observation time relative to report (header) time
    nbtlop = 14  ,& ! LOG( pressure )
    nbtdrh = 15  ,& ! bias correction for relative humidity [/]
    nbtw   = 16  ,& ! vertical velocity [m/s]
    nbtsnr = 17  ,& ! signal to noise ratio
    nbtuac = 18     ! accuracy (std dev from data provider) of horiz. wind [m/s]

!       1.3.2  Body format of ODR of multi-level report flags: 'momlbd'
!              --------------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    mxrbdf =  5  ,& ! body length of multi-level reports
    nbtflg =  1  ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbterr =  2  ,& ! pre-processing status flags  (bit p., see below: 'nb?err')
    nbtqcf =  3  ,& ! threshold quality control flags      (see below: 'nb?qcf')
    nbtlsg =  4  ,& ! level id (bit pattern, as in NetCDF feedobs/feedback file
                    !           see table 'level_sig' in 'mo_fdbk_tables.f90')
    nbtlid =  5     ! level identity          (bit pattern, see below: 'nb?lid')

!       1.3.3  Body format of ODR of surface reports: 'osgbdy'
!              -----------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    mxsbdy = 41  ,& ! body length of single-level reports
    nbsu   =  1  ,& ! u wind component                                   [m/s]
    nbsv   =  2  ,& ! v wind component                                   [m/s]
    nbst   =  3  ,& ! temperature                                        [K]
    nbsrh  =  4  ,& ! relative humidity  (bias corrected)                [/]
    nbsp   =  5  ,& ! pressure (station pressure, if present)            [Pa]
    nbsz   =  6  ,& ! height                                             [m]
    nbsuer =  7  ,& ! error of observed wind component
    nbster =  8  ,& ! error of observed temperature
    nbsqer =  9  ,& ! error of observed relative humidity
    nbszer = 10  ,& ! error of observed height
    nbspst = 11  ,& ! (3-hourly) pressure tendency                       [Pa/3h]
    nbscbs = 12  ,& ! (lowest) cloud base height                         [m]
    nbscl  = 13  ,& ! low       cloud cover        (BUFR Table 020011)   [octas]
    nbscm  = 14  ,& ! mid-level cloud cover        (BUFR Table 020011)   [octas]
    nbsch  = 15  ,& ! high      cloud cover        (BUFR Table 020011)   [octas]
    nbsct  = 16  ,& ! total     cloud cover        (BUFR Table 020011)   [octas]
    nbsvis = 17  ,& ! (horizontal) visibility                            [m]
    nbsff  = 18  ,& ! wind speed                                         [m/s]
    nbsdd  = 19  ,& ! wind direction                                     [deg]
    nbstd  = 20  ,& ! dewpoint temperature                               [K]
    nbsrr1 = 21  ,& ! precipitation amount over 1  hour                  [mm]
    nbsrr3 = 22  ,& ! precipitation amount over 3  hours                 [mm]
    nbsrr6 = 23  ,& ! precipitation amount over 6  hours                 [mm]
    nbsr12 = 24  ,& ! precipitation amount over 12 hours                 [mm]
    nbsr24 = 25  ,& ! precipitation amount over 24 hours                 [mm]
    nbsfgv = 26  ,& ! max. derived equivalent vertical gust (aircraft)   [m/s]
    nbsfg1 = 26  ,& ! max. wind speed of gusts over 1 hour               [m/s]
    nbsfg3 = 27  ,& ! max. wind speed of gusts over 3 hours              [m/s]
    nbsfg6 = 28  ,& ! max. wind speed of gusts over 6 hours              [m/s]
    nbstn  = 29  ,& ! minimum temperature (at 2m, in period 'nbsttr')    [K]
    nbstx  = 30  ,& ! maximum temperature (at 2m, in period 'nbsttr')    [K]
    nbsrad = 31  ,& ! global    solar    radiation, sum over 1 hour      [J/m2]
    nbsrdd = 32  ,& ! diffuse   solar    radiation, sum over 1 hour      [J/m2]
    nbsrdt = 33  ,& ! long-wave downward radiation, sum over 1 hour      [J/m2]
    nbshsw = 34  ,& ! total snow depth                                   [m]
    nbsdrh = 35     ! bias correction for relative humidity [/]
                    !   the following elements are (re-)determined in the
                    !   observation increment part (in src_sing_local.f90)
                    !   whenever the quality control is re-done (rather than
                    !   once when reading in the observation processing -
                    !   however this info is used over several timesteps
                    !   and therefore part of the ODR array)

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    nbsviu = 36  ,& ! 10-m u-wind comp. obs extrapolated to lowest model level
    nbsviv = 37  ,& ! 10-m v-wind comp. obs extrapolated to lowest model level
    nbsvit = 38  ,& !  2-m temperature obs. extrapolated to lowest model level
    nbsviq = 39  ,& !  2-m humidity observ. extrapolated to lowest model level
    nbsvip = 40  ,& ! surface pressure obs. extrapolated to lowest model level
    nbsviz = 41     ! scaled extrapolat. distance for surf. pressure obs [m]
!     these elements have following meaning only when the obs are read from AOF
!     (note that these elements are not used actively in the nudging):
!   nbspr  = 24  ,& ! precipitation amount (over 12 or 24 hrs, not used) [mm]
!   nbsper = 14  ,& !   error of precip amount
!   nbscer = 15  ,& !   error of low cloud cover
!   nbstg  = 34  ,& !   ground temperature (min. T-5cm over 12 hrs)       [K]
!   nbsfgu = 26  ,& ! ( = nbsfg1 ) max. gust speed                        [m/s]
!   nbsfme = 28     ! ( = nbsfg6 ) max. wind speed of 10 minute mean wind [m/s]

!       1.3.4  Body format of ODR of surface report flags: 'mosgbd'
!              ----------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    mxsbdf = 12  ,& ! body length of single-level reports
    nbsflg =  1  ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbserr =  2  ,& ! pre-processing status flags  (bit p., see below: 'nb?err')
    nbsqcf =  3  ,& ! threshold quality control flags      (see below: 'nb?qcf')
    nbslid =  4  ,& ! SYNOP: pressure code (SYNOP)   (code, see below: 'nbslid')
                    ! else : level identity   (bit pattern, see below: 'nb?lid')
    nbscwg =  5  ,& ! combined cloud and weather group (set of classes, s below)
!   nbscfw =  6  ,& ! AOF read: flags for cloud, weather, precip, extreme temp.
                    !     (never set except for bit pos. 20, 28)     (see below) 
    nbswwe =  6  ,& ! NetCDF read, SYNOP: weather and ground group word  (below)
    nbstur =  6  ,& ! NetCDF read, Aircraft: degree of turbulence WMO Tab 011031
                    !   (not contained in merged multi-level aircraft reports !)
    nbsclg =  7  ,& ! general           cloud       group (code)
    nbscl1 =  8  ,& ! first  individual cloud layer group (code)
    nbscl2 =  9  ,& ! second individual cloud layer group (code)
    nbscl3 = 10  ,& ! third  individual cloud layer group (code)
    nbscl4 = 11  ,& ! forth  individual cloud layer group (code)
    nbsttr = 12     ! time periods    (bit pattern of values, see below: nbsttr)

!       1.3.5  Body format of ODR of GPS reports: 'ogpbdy'
!              -------------------------------------------
 
  INTEGER (KIND=iintegers) , PARAMETER  :: &
    mxgbdy = 15  ,& ! body length of GPS reports
    nbgtze =  1  ,& ! error in total zenith delay                        [mm]
    nbgzpd =  2  ,& ! zenith path delay (total zenith delay)             [mm]
    nbgzwd =  3  ,& ! zenith wet delay                                   [mm]
    nbgiwv =  4  ,& ! integrated water vapour                            [mm]
    nbgp   =  5  ,& ! pressure                                           [Pa]
    nbgt   =  6  ,& ! temperature                                        [K]
    nbgrh  =  7  ,& ! relative humidity                                  [/]
    nbgbia =  8  ,& ! bias correction to integrated water vapour         [mm]
    nbgiwa =  9  ,& ! adjusted (bias corrected) integrated water vapour  [mm]
    nbgz   = 10  ,& ! height 
    nbgzer = 11  ,& ! error of observed height
    nbgqer = 12  ,& ! error of observed relative humidity
    nbgdrh = 13  ,& ! bias correction for observed relative humidity
    nbgter = 14  ,& ! error of observed temperature
    nbgviz = 15

!       1.3.6  Body format of ODR of GPS report flags: 'mogpbd'
!              ------------------------------------------------
 
  INTEGER (KIND=iintegers)  , PARAMETER   :: &
    mxgbdf =  4  ,& ! body length of GPS reports
    nbgflg =  1  ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbgerr =  2  ,& ! pre-processing status flags  (bit p., see below: 'nb?err')
    nbgqcf =  3  ,& ! threshold quality control flags      (see below: 'nb?qcf')
    nbglid =  4     ! level identity          (bit pattern, see below: 'nb?lid')

!       1.3.7  Body format of ODR of sat retrieval reports: 'otvbdy'
!              -----------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    mxtbdy =  3  ,& ! body length of multi-level reports
    nbvt   =  1  ,& ! temperature [K]
    nbvrh  =  2  ,& ! relative humidity [/]
    nbvp   =  3     ! pressure [Pa]  (mandatory)

!       1.3.8  Body format of ODR of sat retrieval report flags: 'motvbd'
!              ----------------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    mxtbdf =  2  ,& ! body length of satellite retrieval reports
    nbvflg =  1  ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbverr =  2     ! pre-processing status flags  (bit p., see below: 'nb?err')
!   nbvqcf =  3  ,& ! threshold quality control flags      (see below: 'nb?qcf')
!   nbvlid =  4     ! level identity          (bit pattern, see below: 'nb?lid')


!       ------------------------------------------------------------------------
!       1.4    Bit patterns for packed information in ODR (and VOF) body
!       ------------------------------------------------------------------------

!       1.4.1  Bit patterns for packed info in ODR which should match those
!              used in NetCDF Observation OUTPUT (feedback) File format
!              --> this is replaced by means using 'mo_fdbk_tables.f90' itself

!   variable               meaning                                    word no.
!   --------               -------                                    --------
!   -------------------- ! bit pos. for flags on            --------- nbtlsg
!   ls_surface   = 0  ,& ! surface
!   ls_standard  = 1  ,& ! standard level
!   ls_tropo     = 2  ,& ! tropopause level
!   ls_max       = 3  ,& ! maximum wind level
!   ls_sign      = 4     ! significant level

!       1.4.2  Other bit patt. for packed info in ODR (VOF) body, general words
!              ----------------------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
!   variable          meaning                                         word no.
!   --------          -------                                         --------
!   --------------- ! pre-processing status:  (VOF: nvberf) --------- nb?err
!                   !   bit = 1: obs is active  after pre-processing    !!
!                   !   bit = 0: obs is passive after pre-processing    \/
!           ------- ! threshold quality ctrl: (VOF: nvbqcf) --------- nb?qcf
!                   !   bit = 1: obs is rejected by QC                  !!
!                   !   bit = 0: obs is accepted by QC                  !!
!                   !     upper-air single-level pressure bit flag:     !!
!                   !     bit = 1: obs level is below model orography   \/
!                   ! bit pos. for bit patterns in nb?qcf and nb?err:    
    nvru   =  0  ,& ! bit pos. for status/QC flags for horiz. wind      "
    nvrt   =  1  ,& ! bit pos. for status/QC flags for temperature      "
    nvrq   =  2  ,& ! bit pos. for status/QC flags for humidity         "
    nvrz   =  3  ,& ! bit pos. for status/QC flags for pressure/height  "
    nvrw   =  4  ,& ! bit pos. for status/QC flags for vertical wind    "
    nvriwv =  5  ,& ! bit pos. for status/QC flags for IWV              "
    nvrzpd =  6  ,& ! bit pos. for status/QC flags for zenith path delay"
    nvrspd =  7  ,& ! bit pos. for status/QC flags for slant path delay "
    nvrct  =  8  ,& ! bit pos. for status/QC flags for (total) cloud    "
    nvrcl  =  9  ,& ! bit pos. for status/QC flags for low     cloud    "
    nvrcm  = 10  ,& ! bit pos. for status/QC flags for middle  cloud    "
    nvrch  = 11  ,& ! bit pos. for status/QC flags for high    cloud    "
    nvrcbs = 12  ,& ! bit pos. for status/QC flags for cloud base height"
                    ! temporary flags: not written to VOF file:
    nvrtmp = 13  ,& ! lowest bit pos. with temporary flags              "
    nvrzbc = 13  ,& ! bit pos. for temporary flag: QC ag. LBC pressure  "
    nvrqbc = 14  ,& ! bit pos. for temporary flag: QC against LBC IWV   "
!   --------------- ! bit pos. for flags on   (VOF: nvbmfw) --------- nb?flg
    nvfubp =  0  ,& ! bit pos. for main flag on wind                    "
    nvftbp =  7  ,& ! bit pos. for main flag on temperature             "
    nvfqbp = 14  ,& ! bit pos. for main flag on humidity                "
    nvfzbp = 21  ,& ! bit pos. for main flag on pressure / geopot.      "
    nvfgbp =  0  ,& ! bit pos. for main flag on integr. water vapour    "
    nvfaoc =  7  ,& ! no. of bits occ. by each main flag                "
                    ! inner bit structure (pos./ no.) for each flag     "
                    !   (note that in the source code, it is often
                    !    assumed that 'nvfboc(2:6) = 1'; hence this
                    !    bit length of the flags must not be changed
                    !    here without changing the source modules)
    nvfbps(6) = (/ 0, 2, 3, 4, 5, 6 /) ,&
    nvfboc(6) = (/ 2, 1, 1, 1, 1, 1 /) ,&
                    ! meaning of bit pattern for each flag:             "
                    ! 0,1 - AIRCRAFT wind: roll angle flag
                    !       otherwise: dataset quality flag,
                    !         e.g. for aircraft temperature and humidity
                    !              or wind profiler or RASS temperature
                    !              (for w-prof + RASS, dataset flag also
                    !               relates to vertical velocity)
                    !       0:good; 1:bad; 2:reserved; 3:missing (also roll an.)
                    !   2 - flag raised by blacklist check
                    !   3 - gross error
                    !   4 - not in valid height range:
                    !         surface-level: height or height distance
                    !                        to orography too large;
                    !         upper-air humidity: above 300 hPa level
                    !         upper-air height  : not lowest/surface height obs
                    !   5 - bad reporting practice:
                    !         SYNOP pressure (p): bad reporting practice
                    !         (upper-air height : not measured, derived from p
                    !            --> not used any more)
                    !         aircraft  height  : not measured, derived from p
                    !                             using standard atmosphere
                    !         upper-air height  : without temperature obs
                    !         buoy      wind: zero  wind speed
                    !         radar VAD wind: small wind speed
                    !         aircraft  wind: bad roll angle quality
                    !         dew point     : temperature not active
                    !         mixing ratio  : temperature or pressure not active
                    !         rel. humidity : temperature not active (if needed)
                    !         general: sensor not at appropriate height
                    !         general: measurement duration not appropriate
                    !   6 - gross error for upper-air temperature + wind:
                    !         upper-air temperature: lapse rate check
                    !         upper-air wind: wind speed shear or
                    !                         directional shear check
                    !     - humidity: > 96% and bias-corrected to saturation,
                    !                 or bias corrected (Vaisala RS92)
                    !                 (note: this does not lead to rejection
                    !                  of obs, hence flag not written to fdbk)
                    !     - upper-air pressure: derived from reported height
                    !                           in obs pre-processing
                    !         (--> this is not set for aircraft where p is
                    !              (bi)uniquely derived by using std atmosphere)
    nvflbp = 28  ,& ! bit pos. for level flag: level below surface      "
    nvfloc =  1  ,& ! no. of bits occ. by level flag                    "
!   --------------- ! bit pos. for periods on               --------- nbsttr
    ntxbp  =  0  ,& ! bit pos. for time period for T-max (1,..,24 hr)   "
    ntnbp  =  5  ,& ! bit pos. for time period for T-min (1,..,24 hr)   "
    ntxoc  =  5  ,& ! no. of bits occ. by each time period              "
!   --------------- ! bit pos. for flags on   (VOF: nvblid) --------- nb?lid
                    ! if not SYNOP: level id. bit pattern:
    nvlidp(9) = (/ 0, 1, 2, 3, 4, 5, 6, 7, 8 /)  ,&
    nvlido    = 1   ! no. bits occ. by each indicator in level id.      "
                    !   nvlidp(1) = maximum wind level
                    !   nvlidp(2) = tropopause
                    !   nvlidp(3) = D part
                    !   nvlidp(4) = C part
                    !   nvlidp(5) = B part
                    !   nvlidp(6) = A part
                    !   nvlidp(7) = surface level
                    !   nvlidp(8) = significant wind level
                    !   nvlidp(9) = significant temperature level
                    ! if     SYNOP: pressure code:                    nbslid
                    !   0 - sea level
                    !   1 - station level pressure
                    !   2 - 850mb level geopotential
                    !   3 - 700mb level geopotential
                    !   4 - 500gpm level pressure
                    !   5 - 1000gpm level pressure
                    !   6 - 2000gpm level pressure
                    !   7 - 3000gpm level pressure
                    !   8 - 8000gpm level pressure
                    !   9 - 900mb level geopotential
                    !   10- 1000mb level geopotential
                    !   11- 500mb level geopotential
                    !   12- 925mb level geopotential

!       1.4.3  Bit patterns for 'optional groups' in ODR body 'mosgbd' (and VOF)
!              -----------------------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
!   variable          meaning                                         word no.
!   --------          -------                                         --------
!   --------------- ! bit positions           (VOF: nvbcwg) --------- nbscwg
                    ! combined cloud and weather group
    nvchbp =  0  ,& ! bit position for ch  (type of high cloud)
    nvcmbp =  4  ,& !         "        cm  (type of middle cloud)
    nvclbp = 12  ,& !         "        cl  (type of low cloud)
    nvnhbp = 16  ,& !         "        nh  (cover of low, else of middle cloud)
    nvhbp  =  8  ,& !         "        h   (cloud base height)
    nvnbp  = 20  ,& !         "        n   (total cloud cover)
    nvwwbp = 24  ,& !         "        ww  (present weather)
                    !                      (see VUB WMO Code tables:)
    nvchoc =  4  ,& ! no. of bits occupied by ch    [Code table 0509]
    nvcmoc =  4  ,& !           "             cm    [Code table 0515]
    nvcloc =  4  ,& !           "             cl    [Code table 0513]
    nvnhoc =  4  ,& !           "             nh    [Code table 2700]
    nvhoc  =  4  ,& !           "             h     [Code table 1600]
    nvnoc  =  4  ,& !           "             n     [Code table 2700]
    nvwwoc =  7  ,& !           "             ww    [Code table 4677]
!   --------------- ! bit positions           (VOF: nvbcfw) --------- nbscfw
                    ! combined optional group flag(!) word
                    !   (only for output on VOF, if AOF read)
                    !   bp: var  bp: variable        bit pos: variable
                    !    0: ch   14: v   (visibility)     25: fxgu (max. gusts)
                    !    2: cm   16: rr  (total precip)   26: fxme (max. 10'
                    !    4: h    18: ccl (low cloud cov)            mean wind)
                    !    6: cl   20: refined flag on ccl  27: ffff (global rad)
                    !    8: nh   22: tgtg (ground T)      28: precip measurem.
                    !   10: n    23: tntn (min. T)            duration (VUB Code
                    !   12: ww   24: txtn (max. T)            4019, keys 0-7)

!   --------------- ! bit positions           ----------------------- nbswwe
                    ! --> weather and ground group word
    nvw0bp =  0  ,& ! bit position for ww (present wea.) code (WMO Table 020003)
    nvw1bp =  9  ,& ! bit position for w  (past weather) code (WMO Table 020004)
    nvwtbp = 14  ,& ! bit position for time period of w                  [h]
    nvcqbp = 20  ,& ! bit position for refined quality flag on ccl
    nveebp = 22  ,& ! bit position for state of ground        (WMO Table 020062)
    nrtrbp = 28  ,& ! bit position for code of precipitation
                    !     measurement duration       [Code table 4019, keys 0-7]
    nvw0oc =  9  ,& ! no. bits occupied for ww code           (WMO Table 020003)
    nvw1oc =  5  ,& ! no. bits occupied for w  code           (WMO Table 020004)
    nvwtoc =  5  ,& ! no. bits occupied for time period of w             [h]
    nvcqoc =  2  ,& ! no. bits occupied by refined quality flag on ccl
    nveeoc =  5  ,& ! no. bits occupied for state of ground   (WMO Table 020062)
    nrtroc =  3  ,& ! no. bits occupied by precip obs. duration code
                    !
!   --------------- ! bit positions           ----------------------- nbsclg
                    ! --> general    cloud       group word
!   nxclbp =  0  ,& ! bit position for low/middle cloud amount(WMO Table 020011)
    nctlbp =  4  ,& ! bit position for low    cloud type code (WMO Table 020012)
    nctmbp = 10  ,& ! bit position for middle cloud type code (WMO Table 020012)
    ncthbp = 16  ,& ! bit position for high   cloud type code (WMO Table 020012)
    nxsgbp = 22  ,& ! bit position for vertic. signific. code (WMO Table 008002)
!   nxsgoc =  6  ,& ! no. bits occupied for vert. signf. code (WMO Table 008002)
!   nxcloc =  4  ,& ! no. bits occupied for low/mid cld amount(WMO Table 020011)
!   nxctoc =  6  ,& ! no. bits occupied for   cloud type code (WMO Table 020012)
                    !
!   --------------- ! bit positions           ----------------------- nbscl?
                    ! --> individual cloud layer group words
    nxclbp =  0  ,& ! bit position for cloud amount      code (WMO Table 020011)
    nxctbp =  4  ,& ! bit position for cloud type        code (WMO Table 020012)
    nxbsbp = 10  ,& ! bit position for cloud base height                 [m]
    nxsibp = 24  ,& ! bit position for vertic. signific. code (WMO Table 008002)
    nxcloc =  4  ,& ! no. bits occupied for cloud amount code (WMO Table 020011)
    nxctoc =  6  ,& ! no. bits occupied for cloud type   code (WMO Table 020012)
    nxbsoc = 14  ,& ! no. bits occupied for cloud base height            [m]
    nxsgoc =  6     ! no. bits occupied for vert. signf. code (WMO Table 008002)


!       ------------------------------------------------------------------------
!       1.5    Further quantities related to ODR
!       ------------------------------------------------------------------------

!       1.5.0  Missing data indicators in ODR's
!              --------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    imdi   =  2147483647    ! missing data indicator for ODR integers (2^31 -1)
!   imdi   =  HUGE( 0_iintegers )   ! missing data indicator for ODR integers (2^31 -1)

!       1.5.1  Total number of stored reports in ODR's
!              ---------------------------------------

  INTEGER (KIND=iintegers)  ::    &
    ntotml =  0  ,& ! tot. number of stored multi-level reports
    ntotsg =  0  ,& ! tot. number of stored single-level reports
    ntotgp =  0  ,& ! tot. number of stored GPS reports
    ntottv =  0     ! tot. number of stored satellite retrievals

!       1.5.2  Format of printed station id
!              ----------------------------

  CHARACTER (LEN=ilstidp) ystid      ! obs. station identity to be printed


!       1.5.3  Use of surface-level data with orography differences
!              ----------------------------------------------------

  REAL (KIND=wp)           , PARAMETER  :: &
    fdoro (5) = (/ 1.00_wp ,& ! scaling factor to vertical distances
                   0.25_wp ,& ! between station height and model orography
                   1.00_wp ,& ! for (z-obs > oro-mod):
                   1.00_wp ,& ! 1 - 4: for nudging weights (see 'doromx)
                   0.50_wp /) !   5  : when assigning a grid pt. to the obs.


!-------------------------------------------------------------------------------
!  Section 2 : Observation Data Records (ODR)
!-------------------------------------------------------------------------------

  REAL    (KIND = wp)       , ALLOCATABLE :: &
    omlbdy (:,:,:) ,& ! body   of multi-level ODR
    omlhed   (:,:) ,& ! header of multi-level ODR
    osgbdy   (:,:) ,& ! body   of single-level ODR
    osghed   (:,:) ,& ! header of single-level ODR
    ogpbdy   (:,:) ,& ! body   of GPS ODR
    ogphed   (:,:) ,& ! header of GPS ODR
    otvbdy (:,:,:) ,& ! body   of satellite retrieval ODR
    otvhed   (:,:)    ! header of satellite retrieval ODR

  INTEGER (KIND = iintegers), ALLOCATABLE :: &
    momlbd (:,:,:) ,& ! body   of multi-level ODR
    momlhd   (:,:) ,& ! header of multi-level ODR
    mosgbd   (:,:) ,& ! body   of single-level ODR
    mosghd   (:,:) ,& ! header of single-level ODR
    mogpbd   (:,:) ,& ! body   of GPS ODR
    mogphd   (:,:) ,& ! header of GPS ODR
    motvbd (:,:,:) ,& ! body   of satellite retrieval ODR
    motvhd   (:,:)    ! header of satellite retrieval ODR

  CHARACTER (LEN = ilstid)  , ALLOCATABLE :: &
    yomlhd     (:) ,& ! header of multi-level ODR
    yosghd     (:) ,& ! header of single-level ODR
    yogphd     (:) ,& ! header of GPS ODR
    yotvhd     (:)    ! header of satellite retrieval ODR

!-------------------------------------------------------------------------------
! Section 3:  Simulated Observation Record (SOR)
!-------------------------------------------------------------------------------

!      3.1  Format of SOR (for model values or observation increments)
!           -------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    mxsoml =  5  ,& ! SOR body length for multi-level reports
    mxsosg = 13  ,& ! SOR body length for single-level reports
    mxsogp =  2  ,& ! SOR body length for GPS reports
    mxsotv =  4  ,& ! SOR body length for satellite retrieval reports
    mxsops =  1  ,& ! SOR header length for multi-level reports
    nso_u  =  1  ,& ! u wind component                               [m/s]
    nso_v  =  2  ,& ! v wind component                               [m/s]
    nso_t  =  3  ,& ! temperature                                    [K]
    nso_rh =  4  ,& ! relative humidity                              [%]
    nso_p  =  5  ,& ! pressure (sfc.) | geopotential (upper-air)  [Pa] | [m2/s2]
    nso_ct =  6  ,& ! total cloud cover                              [octas]
    nso_cl =  7  ,& ! low cloud cover                                [octas]
    nso_cm =  8  ,& ! mid-level cloud cover                          [octas]
    nso_ch =  9  ,& ! high cloud cover                               [octas]
    nso_cbs= 10  ,& ! cloud base height (above surface)              [m]
    nso_ff = 11  ,& ! wind speed                                     [m/s]
    nso_dd = 12  ,& ! wind direction                                 [deg]
    nso_td = 13  ,& ! dewpoint temperature                           [K]
    nso_iq =  1  ,& ! integrated water vapour                        [mm]
    nso_zpd=  2  ,& ! zenith total path delay                        [mm]
    nso_ps =  1     ! pressure                                       [Pa]

!      3.2  SOR arrays (for model values or observation increments)
!           ----------
!           (length of first dimension is given by 'SOR body length' above,
!            length of last dimension by ODR report size ('maxmll', etc),
!            length of intermediate dim. by number of vertical levels in ODR)

  REAL    (KIND = wp)       , ALLOCATABLE :: &
    smlbdy (:,:,:) ,& ! body of multi-level SOR
    ssgbdy   (:,:) ,& ! body of single-level SOR
    sgpbdy   (:,:) ,& ! body of GPS (IWV) SOR
    stvbdy (:,:,:) ,& ! body of satellite retrieval SOR
    dmlhed   (:,:)    ! single-level part of multi-level SOR

!-------------------------------------------------------------------------------
!  Section 4 : Masking constants
!-------------------------------------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    nibits(31)    & ! masking constants   (2^j -1 , j=1,..31)
               = (/           1,           3,           7,          15,        &
                             31,          63,         127,         255,        &
                            511,        1023,        2047,        4095,        &
                           8191,       16383,       32767,       65535,        &
                         131071,      262143,      524287,     1048575,        &
                        2097151,     4194303,     8388607,    16777215,        &
                       33554431,    67108863,   134217727,   268435455,        &
                      536870911,  1073741823,  2147483647             /)
!   nibits(32)    & ! masking constants
!              = (/                      O"1",                   O"3",         &
!                                        O"7",                  O"17",         &
!                                       O"37",                  O"77",         &
!                                      O"177",                 O"377",         &
!                                      O"777",                O"1777",         &
!                                     O"3777",                O"7777",         &
!                                    O"17777",               O"37777",         &
!                                    O"77777",              O"177777",         &
!                                   O"377777",              O"777777",         &
!                                  O"1777777",             O"3777777",         &
!                                  O"7777777",            O"17777777",         &
!                                 O"37777777",            O"77777777",         &
!                                O"177777777",           O"377777777",         &
!                                O"777777777",          O"1777777777",         &
!                               O"3777777777",          O"7777777777",         &
!                              O"17777777777",         O"37777777777"/)


!-------------------------------------------------------------------------------

END MODULE data_obs_record
