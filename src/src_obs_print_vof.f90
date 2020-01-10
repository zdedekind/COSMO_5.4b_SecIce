!+ Source module for the observation processing in the data assimilation mode
!-------------------------------------------------------------------------------

MODULE src_obs_print_vof

!-------------------------------------------------------------------------------
! Description:
!   This module writes observation reports and increments (or model values
!   projected onto the observation space) to an ASCII file for verification
!   purposes.
!
!   This module contains the following module procedures:
!    - obs_print_vof    : write ASCII formatted file 'VOF'
!
!   This module also contains an elemental function, formerly statement funct.:
!    - rmod             : MOD function for positive reals
!
!   It uses from:
!    - src_obs_cdfin_util: - obs_gather_buffer
!    - utilities:          - get_utc_date
!                          - sortrx
!    - parallel_utilities: - global_values
!    - environment:        - model_abort
!
!   Data modules used: - data_parameters
!                      - data_obs_lib_cosmo
!                      - data_obs_record
!                      - mo_fdbk_tables
!
!   Note: This module contains MPI-based communication between PE's.
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! V4_22        2012/01/31 Christoph Schraff
!  Initial release, extracted from module 'src_obs_processing' and adapted,
!    e.g. subroutine argument list introduced to make routine independent from
!    'data_nudge_all', 'data_modelconfig', 'data_parallel', and 'data_io'.
!  - Increased efficiency e.g. by use of quicksort instead of bubblesort.
!  - Bug fix in sorting reports.
!  - 'merged report' bit flag set in 'nhflag' word.
!  - In VOF Version 2, flag words are written in octal instead of integer
!    format. Obs code type, latitude, observation status flag and quality
!    control flag, as well as wind, temperature and pressure increments
!    are written with an increased number of digits.
!  - In VOF Version 3, some entries have an additional digit to enhance
!    readability of VOF (and facilitate parsing).
! V4_24        2012/06/22 Hendrik Reich
!  Adapted length of strings for date variables
! V4_25        2012/09/28 Christoph Schraff
!  For the lmstat-postprocessing write only 10 digits of the date strings
! V4_27        2013/03/19 Christoph Schraff
!  VOF Version 4: write date strings with 14 digits. Distinguish between 'NUDGE'
!  and other model run types e.g. for 'initial date and hour of the run'.
!  Allow for run times unequal to a full hour.
! V4_28        2013/07/12 Christoph Schraff
!  - Temporary flags for LBC QC checks processed.
!  - Full simulated obs values instead of increments are now stored in 'smlbdy'.
!    Increments for the VOF need to be computed now, using height obs values.
!    The condititon to write a radiosonde height increment is now the existence
!    of the height obs instead of the temperature obs, i.e. height increments
!    are more rare now.
!  - Also write simulated ZPD to VOF. (For this purpose, no new VOF version is
!    defined - lmstat will not read and write simulated ZPD.)
!  - For wind profiler and RASS, unset dataset flag, if it refers to vertical
!    velocity only, and horizontal wind resp. temperature is missing.
!  - Statement functions replaced by elemental or intrinsic functions.
! V4_29        2013/10/04 Christoph Schraff
!  Bug fix: condition 'lvofoi' removed for de-allocation of 'ssgbdy' etc.
! V5_1         2014-11-28 Christoph Schraff, Oliver Fuhrer
!  'lrefend' only if (mruntyp == 2); run type 'forecast' added.
!  Replaced ireals by wp (working precision) (OF)
! V5_4         2016-03-10 Christoph Schraff
!  Removal of VOF variables related to the AOF interface: 'nvhqfl' and 'nvbcfw'
!  with related bit patterns, 'nvbpr', 'nvbtg', 'nvbfgu', 'nvbfme'.
!  (All of them have already been replaced in previous versions by the following
!   variables if the obs are read from NetCDF files: 'nvhflg', 'nvbwwe',
!   'nvbr12', 'nvbhsw', 'nvbfg1', 'nvbfg6'.) 
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

USE data_obs_lib_cosmo, ONLY :  &

! 1. General parameters
! ---------------------

    c0         ,& ! standard real constant 0.0
    c1         ,& ! standard real constant 1.0
    c05        ,& ! standard real constant 0.5
    rmdi       ,& ! =-1.E31_wp : commonly used missing data indicator
    rmdich     ,& ! =-1.E30_wp : commonly used check value for miss data
    epsy       ,& ! = 1.E-8_wp : commonly used very small value > 0

! 4. I/O device numbers for obs processing / nudging and file names
! -----------------------------------------------------------------

    yuverif    ,& ! VOF (output): verification observation file (obs.
                  !      incl. quality control flag for verification)
    nuverif    ,& ! VOF (output): verification observation file (observations
                  !   incl. quality control flag for verification)

! 5. CMA observation type and code type numbers
! ---------------------------------------------

!   nsynop     ,& ! SYNOP reports
    nairep     ,& ! AIREP reports (all aircraft reports)
    nsatob     ,& ! SATOB reports
!   ndribu     ,& ! DRIBU reports
!   ntemp      ,& ! TEMP  reports
!   npilot     ,& ! PILOT reports
!   nsatem     ,& ! SATEM reports
!   nsattv     ,& ! SATEM reports
    ngps       ,& ! GPS   reports
    nacar      ,& !   acar  report
    nacar_vof     !   acar  report (only in VOF)

! end of data_obs_lib_cosmo

!-------------------------------------------------------------------------------

USE data_obs_record, ONLY :   &
 
!       1.1    ODR header format
!       ------------------------

!       1.1.1  Header formats of ODR reports:'omlhed','osghed','ogphed','otvhed'
!              -----------------------------------------------------------------
!   mxrhed,     & ! header length of multi-level reports
!   mxshed,     & ! header length of single-level reports
!   mxghed,     & ! header length of GPS reports
!   mxthed,     & ! header length of satellite retrieval reports
    nhilon,     & ! longitude of observing station
    nhjlat,     & ! latitude  of observing station
    nhalt ,     & ! station altitude [m]
    nhtime,     & ! time of observat. in forecast hours
    nhsurf,     & ! height of model grid pt. to which obs. is assigned

!       1.1.2  Header formats of ODR reports:'momlhd','mosghd','mopghd','motvhd'
!              -----------------------------------------------------------------
!   mxrhdf,     & ! header length of multi-level reports
!   mxshdf,     & ! header length of single-level reports
!   mxghdf,     & ! header length of GPS reports
!   mxthdf,     & ! header length of satellite retrieval reports
    nhitot,     & ! global x-coord. of grid pt. assigned to obs
    nhjtot,     & ! global y-coord. of grid pt. assigned to obs
    nhobtp,     & ! observation type
    nhcode,     & ! code type
    nhschr,     & ! station characteristics                      (see 1.1.4)
    nhflag,     & ! report flags (obs type, surf., alt., sta ID) (see 1.2.1)
    nhpass,     & ! flag for report being set to 'passive'       (see 1.1.4)
    nhqcfw,     & ! status of QC and of writing to feedobs files, and
                  !   QC flags for surface pressure increments from TEMPs
!   nhcorr,     & ! update sequence number (station correction indicator)
!   nhcat ,     & ! data     category (from BUFR Section 1)
!   nhcats,     & ! data sub-category (from BUFR Section 1)
!   nhcent,     & ! originating centre  +  (1000* sub-centre)
!   nhstid,     & ! station identity number
!   nhhrmn,     & ! absolute exact observation time [hhmm]
!   nhsyhr,     & ! absolute nominal (synoptic) observation time [yymmddhh]
!   nhstyp,     & ! sing-lv obs: station type (buoy: MQOBL, BUFR Table 002149,
                  !                            else: NIX  , BUFR Table 002001)
!   nhrtyp,     & ! radiosonde type    (NRARA, see WMO common code Table C2)
    nhnlev,     & ! number of obs. levels (for multi-level reports)
    nhvqcf,     & ! for satellite retrieval: threshold quality control flags
!   nhtrac,     & ! tracking technique (NSASA, see WMO common code Table C7)
!   nhrad ,     & ! solar and IR radiation correction (NSR, BUFR Table 002013)
!   nhna4 ,     & ! instrument type                   (NA4, BUFR Table 002003)

!       1.1.3  Header formats of ODR reports:'yomlhd','yosghd','yopghd','yotvhd'
!              -----------------------------------------------------------------
    ilstid,     & ! character length of the station identity
    ilstidp       ! char. length used for printing the station ID
                  ! Note: (ilstid >= ilstidg >= ilstidp), cf. data_nudge_gather

USE data_obs_record, ONLY :   &

!       1.2    Bit patterns for packed information in ODR (and VOF) header
!       ------------------------------------------------------------------
    nvpabp        ! bit pos. for report set passive since it is     nvschr
                  !              used in a multi-level pseudo report
!   nvpsbp,     & ! bit pos. for report set passive since at least    "
                  !              1 of the following 5 flags or
                  !              the flight track flag applies
!   nvobbp,     & ! bit pos. for flag: 'station location out of       "
                  !                     user-specified area'
!   nvalbp,     & ! bit pos. for flag: 'distance model orography -    "
                  !                     station altitude too large'
!   nvbkbp,     & ! bit pos. for flag: 'blacklisted station (ship)'   "
!   nvexbp,     & ! bit pos. for flag: 'observation or code type      "
                  !                     excluded at station location'
!   nvrdbp,     & ! bit pos. for flag: 'redundant report'             "
!   nvsebp,     & ! bit pos. for report located at sea grid pt.       "
!   nvscbp,     & ! bit pos. for station correction indicator         "
!   nvssbp,     & ! bit pos. for station suspicion indicator          "
!   nvinbp,     & ! bit pos. for instrument specification indicator   "
!   nvinoc,     & ! no. of bits occ. by instrument specif. ind.       "
!   nvhtbp,     & ! bit pos. for flight track error flag              "
!   nvhtoc,     & ! no. of bits occ. by flight track error flag       "
!   nvhhbp,     & ! bit pos. for flight thinning flag                 "
!   nvhhoc,     & ! no. of bits occ. by flight thinning flag          "
!   nvapbp,     & ! bit pos. for phase of flight (aircraft)           "
!   nvapoc,     & ! no. of bits occ. by phase of flight               "
!   nvaabp,     & ! bit pos. for aircraft roll angle (code)           "
!   nvaaoc        ! no. of bits occ. by aircraft roll angle           "

USE data_obs_record, ONLY :   &

!       1.3    ODR body format
!              ---------------

!       1.3.0  Number of levels in multi-level ODR 'omlbdy', 'momlbd'
!              ------------------------------------------------------
!   maxrsl,     & ! max. number of levels in multi-level ODR


!       1.3.1  Body format of ODR of multi-level reports: 'omlbdy'
!              ---------------------------------------------------
!   mxrbdy,     & ! body length of multi-level reports
    nbtu  ,     & ! u wind component [m/s]
    nbtv  ,     & ! v wind component [m/s]
    nbtt  ,     & ! temperature [K]
    nbtrh ,     & ! relative humidity [/]
    nbtp  ,     & ! pressure [Pa]
    nbtz  ,     & ! height [m]
!   nbtuer,     & ! error of observed wind component
!   nbtter,     & ! error of observed temperature
!   nbtqer,     & ! error of observed rel. humidity
!   nbtzer,     & ! error of observed height
                  !  Note: nbt?er are set to the negative rms errors, if the
                  !  observations have not passed the threshold quality control
!   nbtzio,     & ! longitude in grid pt. units
!   nbtzjo,     & ! latitude  in grid pt. units
!   nbtlop,     & ! LOG( pressure )
!   nbtdrh,     & ! bias correction for relative humidity [/]
    nbtw  ,     & ! vertical velocity [m/s]
!   nbtsnr,     & ! signal to noise ratio


!       1.3.2  Body format of ODR of multi-level report flags: 'momlbd'
!              --------------------------------------------------------
!   mxrbdf,     & ! body length of multi-level reports
    nbtflg,     & ! main flag word          (bit pattern, see below: 'nb?flg')
    nbterr,     & ! pre-processing status flags  (bit p., see below: 'nb?err')
    nbtqcf,     & ! threshold quality control flags      (see below: 'nb?qcf')
!   nbtlsg,     & ! level id (bit pattern, as in NetCDF statistics file)
    nbtlid        ! level identity          (bit pattern, see below: 'nb?lid')

USE data_obs_record, ONLY :   &

!       1.3.3  Body format of ODR of surface reports: 'osgbdy'
!              -----------------------------------------------
!   mxsbdy,     & ! body length of single-level reports
    nbsu  ,     & ! u wind component                                   [m/s]
    nbsv  ,     & ! v wind component                                   [m/s]
    nbst  ,     & ! temperature                                        [K]
    nbsrh ,     & ! relative humidity                                  [/]
    nbsp  ,     & ! pressure                                           [Pa]
    nbsz  ,     & ! height                                             [m]
!   nbsuer,     & ! error of observed wind component
!   nbster,     & ! error of observed temperature
!   nbsqer,     & ! error of observed relative humidity
!   nbszer,     & ! error of observed height
    nbspst,     & ! (3-hourly) pressure tendency                       [Pa/3h]
!   nbscbs,     & ! (lowest) cloud base height                         [m]
    nbscl ,     & ! low       cloud cover        (BUFR Table 020011)   [octas]
!   nbscm ,     & ! mid-level cloud cover        (BUFR Table 020011)   [octas]
!   nbsch ,     & ! high      cloud cover        (BUFR Table 020011)   [octas]
!   nbsct ,     & ! total     cloud cover        (BUFR Table 020011)   [octas]
    nbsvis,     & ! (horizontal) visibility                            [m]
!   nbsrr1,     & ! precipitation amount over 1  hour                  [mm]
!   nbsrr6,     & ! precipitation amount over 6  hours                 [mm]
    nbsr12,     & ! precipitation amount over 12 hours                 [mm]
!   nbsr24,     & ! precipitation amount over 24 hours                 [mm]
!   nbsfgv,     & ! max. derived equivalent vertical gust (aircraft)   [m/s]
    nbsfg1,     & ! max. wind speed of gusts over 1 hour               [m/s]
    nbsfg6,     & ! max. wind speed of gusts over 6 hours              [m/s]
    nbstn ,     & ! minimum temperature (at 2m during past 12 hrs)     [K]
    nbstx ,     & ! maximum temperature (at 2m during past 12 hrs)     [K]
    nbshsw,     & ! total snow depth                                   [m]
!   nbsdrh,     & ! bias correction for relative humidity [/]
    nbsrad        ! global radiation, sum over 1 hour                  [J/m2]

USE data_obs_record, ONLY :   &

!       1.3.4  Body format of ODR of surface report flags: 'mosgbd'
!              ----------------------------------------------------
!   mxsbdf,     & ! body length of single-level reports
    nbsflg,     & ! main flag word          (bit pattern, see below: 'nb?flg')
    nbserr,     & ! pre-processing status flags  (bit p., see below: 'nb?err')
    nbsqcf,     & ! threshold quality control flags      (see below: 'nb?qcf')
    nbslid,     & ! SYNOP: pressure code (SYNOP)   (code, see below: 'nbslid')
                  ! else : level identity   (bit pattern, see below: 'nb?lid')
    nbscwg,     & ! combined cloud and weather group (set of classes, s below)
    nbswwe,     & ! NetCDF read, SYNOP: weather and ground group word  (below)
!   nbstur,     & ! NetCDF read, Aircraft: degree of turbulence WMO Tab 011031
                  !   (not contained in merged multi-level aircraft reports !)
!   nbsclg,     & ! general           cloud       group (code)
!   nbscl1,     & ! first  individual cloud layer group (code)
!   nbscl2,     & ! second individual cloud layer group (code)
!   nbscl3,     & ! third  individual cloud layer group (code)
!   nbscl4,     & ! forth  individual cloud layer group (code)

!       1.3.5  Body format of ODR of GPS reports: 'ogpbdy'
!              -------------------------------------------
!   mxgbdy,     & ! body length of GPS reports
!   nbgtze,     & ! error in total zenith delay [mm]
    nbgzpd,     & ! zenith path delay (total zenith delay)             [mm]
    nbgzwd,     & ! zenith wet delay [mm]
    nbgiwv,     & ! integrated water vapour [mm]
    nbgp  ,     & ! pressure [Pa]
    nbgt  ,     & ! temperature [K]
    nbgrh ,     & ! relative humidity [/]
!   nbgbia,     & ! bias correction to integrated water vapour [mm]
    nbgiwa,     & ! adjusted (bias corrected) integrated water vapour [mm]

!       1.3.6  Body format of ODR of GPS report flags: 'mogpbd'
!              -------------------------------------------------
!   mxgbdf,     & ! body length of GPS reports
    nbgflg,     & ! main flag word          (bit pattern, see below: 'nb?flg')
    nbgerr,     & ! pre-processing status flags  (bit p., see below: 'nb?err')
    nbgqcf,     & ! threshold quality control flags      (see below: 'nb?qcf')
    nbglid,     & ! level identity          (bit pattern, see below: 'nb?lid')

!       1.3.7  Body format of ODR of sat retrieval reports: 'otvbdy'
!              -----------------------------------------------------
!   mxtbdy,     & ! body length of multi-level reports
    nbvt  ,     & ! temperature [K]
    nbvrh ,     & ! relative humidity [/]
    nbvp  ,     & ! pressure [Pa]  (mandatory)

!       1.3.8  Body format of ODR of sat retrieval report flags: 'motvbd'
!              ----------------------------------------------------------
!   mxtbdf,     & ! body length of sat retrieval reports
    nbvflg,     & ! main flag word          (bit pattern, see below: 'nb?flg')
    nbverr        ! pre-processing status flags  (bit p., see below: 'nb?err')

USE data_obs_record, ONLY :   &

!       1.4    Bit patterns for packed information in ODR (and VOF) body
!              ---------------------------------------------------------

!       1.4.2  Other bit patt. for packed info in ODR (VOF) body, general words
!              ----------------------------------------------------------------
!   nvru  ,     & ! bit pos. for status/QC flags for horiz. wind  nb?err/nb?qcf
!   nvrt  ,     & ! bit pos. for status/QC flags for temperature        "
    nvrq  ,     & ! bit pos. for status/QC flags for humidity           "
    nvrz  ,     & ! bit pos. for status/QC flags for pressure/height    "
    nvrw  ,     & ! bit pos. for status/QC flags for vertical wind      "
    nvriwv,     & ! bit pos. for status/QC flags for IWV                "
    nvrzpd,     & ! bit pos. for status/QC flags for zenith path delay  "
!   nvrspd,     & ! bit pos. for status/QC flags for slant path delay   "
!   nvrct ,     & ! bit pos. for status/QC flags for (total) cloud      "
!   nvrcl ,     & ! bit pos. for status/QC flags for low     cloud      "
!   nvrcm ,     & ! bit pos. for status/QC flags for middle  cloud      "
!   nvrch ,     & ! bit pos. for status/QC flags for high    cloud      "
    nvrtmp,     & ! lowest bit pos. with temporary flags              "
    nvrzbc,     & ! bit pos. for temporary flag: QC ag. LBC pressure  "
    nvrqbc,     & ! bit pos. for temporary flag: QC against LBC IWV   "
    nvfubp,     & ! bit pos. for main flag on wind                    nb?flg
    nvftbp,     & ! bit pos. for main flag on temperature               "
!   nvfqbp,     & ! bit pos. for main flag on humidity                  "
    nvfzbp,     & ! bit pos. for main flag on pressure / geopot.        "
!   nvfaoc,     & ! no. of bits occ. by each main flag                  "
    nvfbps,     & ! bit pattern for main flags:                         "
    nvfboc,     & ! no. of bits occ. for main flags                     "
!   nvflbp,     & ! bit pos. for level flag: level below surface        "
!   nvfloc,     & ! no. of bits occ. by level flag                      "
    nvlidp,     & ! level id. bit pattern                             nb?lid
    nvlido,     & ! no. bits occ. by each indicator in level id.        "

!       1.5    Further quantities related to ODR
!              ---------------------------------
    imdi  ,     & ! missing data indicator for ODR integers (2^31-1)
    ntotml,     & ! tot. number of stored multi-level reports
    ntotsg,     & ! tot. number of stored single-level reports
    ntotgp,     & ! tot. number of stored GPS reports
    ntottv        ! tot. number of stored satellite retrievals

USE data_obs_record, ONLY :   &

!       2.     Observation data records (ODR)
!       -------------------------------------
    omlbdy,     & ! body of multi-level ODR
    momlbd,     & ! body of multi-level ODR
    omlhed,     & ! header of multi-level ODR
    momlhd,     & ! header of multi-level ODR
    osgbdy,     & ! body of single-level ODR
    mosgbd,     & ! body of single-level ODR
    osghed,     & ! header of single-level ODR
    mosghd,     & ! header of single-level ODR
    ogpbdy,     & ! body of GPS ODR
    mogpbd,     & ! body of GPS ODR
    ogphed,     & ! header of GPS ODR
    mogphd,     & ! header of GPS ODR
    otvbdy,     & ! body of satellite retrieval ODR
    motvbd,     & ! body of satellite retrieval ODR
    otvhed,     & ! header of satellite retrieval ODR
    motvhd,     & ! header of satellite retrieval ODR
    yomlhd,     & ! header of multi-level ODR
    yosghd,     & ! header of single-level ODR
    yogphd,     & ! header of GPS ODR
    yotvhd,     & ! header of satellite retrieval ODR

!       3.1  Format of SOR (for model values or observation increments)
!            -------------
    nso_u      ,& ! u wind component                               [m/s]
    nso_v      ,& ! v wind component                               [m/s]
    nso_t      ,& ! temperature                                    [K]
    nso_rh     ,& ! relative humidity                              [%] 
    nso_p      ,& ! pressure (sfc.) | geopotential (upper-air)  [Pa] | [m2/s2]
    nso_ct     ,& ! total cloud cover                              [octas]
    nso_cl     ,& ! low cloud cover                                [octas]
    nso_iq     ,& ! integrated water vapour (increment)            [mm]
    nso_zpd    ,& ! zenith total path delay                        [mm]
    nso_ps     ,& ! pressure                                       [Pa]

!       3.2  SOR arrays (for model values or observation increments)
!            ----------
    smlbdy     ,& ! body of multi-level SOR
    ssgbdy     ,& ! body of single-level SOR
    sgpbdy     ,& ! body of GPS (IWV) SOR
    stvbdy     ,& ! body of satellite retrieval SOR
    dmlhed        ! single-level part of multi-level SOR

! end of data_obs_record

!-------------------------------------------------------------------------------

  USE mo_fdbk_tables,          ONLY :  &

    FL_MERGE         ! (report used for) merged report (only, e.g. TEMP ABCD
                     !   or single-level aircraft used for multi-level report)

! end of mo_fdbk_tables

!-------------------------------------------------------------------------------

 USE environment,              ONLY :  &
    model_abort        ! aborts the program in case of errors

!-------------------------------------------------------------------------------

 USE parallel_utilities,       ONLY :  &
    global_values      ! computes global values by operating on local arrays

!-------------------------------------------------------------------------------

 USE utilities,                ONLY :  &
    get_utc_date    ,& ! actual date of the forecast in different forms
    sortrx             ! permutation vector for (quick-)sorting a given array

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_util,       ONLY :  &
    obs_gather_buffer  ! gathers elements from buffer arrays from all nodes

!-------------------------------------------------------------------------------
! Local declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE

!===============================================================================

!-------------------------------------------------------------------------------
! Public and Private Subroutines
!-------------------------------------------------------------------------------

CONTAINS


!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_print_vof" for writing to ASCII verif file VOF
!-------------------------------------------------------------------------------

SUBROUTINE obs_print_vof ( lastpr, lrefend, ydate_ref, hversta, hverend        &
                         , mruntyp, mveripr                                    &
                         , num_compute, my_cart_id, icomm_cart, imp_integers   &
                         , dtqc, qcc, qcvf, qccsu, qctfpr                      &
                         , ie_tot, je_tot, ke_tot, pollat, pollon, dlat, dlon  &
                         , startlat_tot, startlon_tot )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "src_obs_print_vof" writes out the
!   observations including the result of the threshold quality control to
!   a separate file in ASCII, when the observation time coincides with the
!   current model time.
!
! Method:
!   If the program is run in parallel mode, the data to be printed are
!   gathered from all PEs at the root node.
!   24.08.98
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! S-story:
! Version    Date       Name
! ---------- ---------- ----
! 1.19       1998/12/11 Christoph Schraff
!  Initial release
! 1.21       1999/01/25 Christoph Schraff
!  Sorting reports temporally first. Height difference replaced by model height.
! 1.27       1999/03/29 Christoph Schraff
!  Revised VOF format, optional groups included. 6-bit holleriths removed.
! 1.28       1999/04/19 Christoph Schraff
!  Revised title line of VOF (to facilitate the search for this title line).
! 1.29       1999/05/11 Ulrich Schaettler
!  Adapted interfaces to utility-modules and prepared use of MPE_IO
! 1.31       1999/07/01 Christoph Schraff
!  MPI calls related to fixed-length arrays replaced by parallel_utility calls;
!  input to 'global_values' adjusted.
! 1.36       2000/02/24 Christoph Schraff
!  Ensuring equivalence of flags with same meaning occurring in different words.
!  Optional writing to VOF of observation increments.
! 1.38       2000/04/06 Christoph Schraff
!  Resetting upper-air quality control flag to missing after printing. Printing
!  pressure tendency, geopotential instead of height deviations, and extended
!  station characteristics to VOF.
! 1.39       2000/05/03 Ulrich Schaettler
!  Changed names for variables concerned to latitude or longitude.
! 1.40       2000/05/23 Christoph Schraff
!  Gravity acceleration for VOF as defined for BUFR.
! 2.4        2001/01/29 Christoph Schraff
!  Formal change of QC variables for output only.
! 2.13       2002/01/18 Christoph Schraff
!  Bug correction: 'dmlhed' used only if defined (i.e. if (lvofoi)).
!  Cloud and weather group word changed from real to integer.
! 2.17       2002/05/08 Ulrich Schaettler
!  Implemented a check for the range of longitude variables
! 3.3        2003/04/22 Christoph Schraff
!  Introduction of GPS reports.
! 3.6        2003/12/11 Christoph Schraff
!  Introduction of SATOB reports. Error flag for GPS TZD revised.
!  Ensuring error message at model_abort.
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!===============================================================================
!
! Declarations:
!
!===============================================================================

IMPLICIT NONE

!===============================================================================

! Subroutine arguments:
! --------------------

  ! these variables are important for the correct functioning of this routine

  LOGICAL                  , INTENT (IN)  :: &
    lastpr      ,& ! last time, something is written to VOF in this model run
    lrefend        ! reference time of model run = end of model run (type NUDGE)

  INTEGER (KIND=iintegers) , INTENT (IN)  :: &
    mruntyp     ,& ! = -1: print only observations ; =2: also print increments
                   ! >= 0: print also increments (only if >= 0: arrays 'smlbdy'
                   !             etc. have to be allocated and filled previously
                   !             and are used and de-allocated here)
                   !       == 2: nudge
                   !       == 1: first guess
                   !       == 0 (or >= 3): forecast
    mveripr     ,& ! 3  : type of verification/observation file(s) written
                   ! (1: ncfeed only ; 2: VOF only ; 3: ncfeed and VOF)
    num_compute ,& ! number of compute PEs
    my_cart_id  ,& ! rank of this subdomain in the cartesian communicator
    icomm_cart  ,& ! communicator for the virtual cartesian topology
    imp_integers   ! INTEGER type used in the model for MPI

  CHARACTER (LEN=14)       , INTENT (IN)  :: &
    ydate_ref      ! reference date (e.g. start of the model run)
                   !                [yyyymmddhhmmss]
                   !                (year, month, day, hour, min., sec.)

  REAL    (KIND=wp)        , INTENT (IN)  :: &
    hversta     ,& ! start of verification period in 'model integration hours'
    hverend        ! end of verification period in 'model integration hours'

  ! following variables are only used for supplementary info in VOF file header

  INTEGER (KIND=iintegers) , INTENT (IN)  :: &
    ie_tot      ,& ! number of grid points in zonal direction
    je_tot      ,& ! number of grid points in meridional direction
    ke_tot         ! number of grid points in vertical direction

  REAL    (KIND=wp)        , INTENT (IN)  :: &
    dtqc        ,& ! timestep (in [s]) for the threshold quality control
    qcc    (4)  ,& ! constant parts of the quality control thresholds (='QCT')
    qcvf   (4)  ,& ! multiplication factor to the vertically varying part of QCT
    qccsu  (4)  ,& ! constant parts of the quality control thresholds (='QCT')
    qctfpr (8)  ,& ! for VOF output only: time factor for QC thresholds
    pollat      ,& ! latitude  of the rotated north pole (in degrees, N>0)
    pollon      ,& ! longitude of the rotated north pole (in degrees, E>0)
    dlat        ,& ! grid point distance in meridional direction (in degrees)
    dlon        ,& ! grid point distance in zonal      direction (in degrees)
    startlat_tot,& ! rotated latitude   \  of the lower left grid point of the
    startlon_tot   ! rotated longitude  /  total domain (in degrees, N,E>0)

! Local parameters:
! ----------------

  ! ------------------------------------------
  ! VOF (Verification Observation File) format
  ! ------------------------------------------
  ! VOF (Verification Observation File) header format (for packed info, see ODR)
  ! -------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    mxvhed = 14  ,& ! header length of all reports
                    !
!   nvhlev =  1  ,& ! number of vertical levels (single-level: 0)
    nvhlon =  2  ,& ! longitude of observing station                 [1/100 deg]
    nvhlat =  3  ,& ! latitude  of observing station                 [1/100 deg]
    nvhtim =  4  ,& ! time of relative to reference hour             [min]
    nvhalt =  5  ,& ! station altitude                               [m]
    nvhsfc =  6  ,& ! height of model orography at station           [m]
    nvhobt =  7  ,& ! observation type
    nvhcdt =  8  ,& ! code type
    nvhsch =  9  ,& ! station characteristics                (see ODR: 'nhschr')
    nvhflg = 10  ,& ! report flags (obs type, surface, altitude, station ID)
                    !   (only if reading from NetCDF)        (see ODR: 'nhflag')
    nvhpas = 11  ,& ! flag for report being set to 'passive' (see ODR: 'nhpass')
    nvhqcp = 12  ,& ! flag of threshold quality control for pressure increments:
                    ! - derived from multi-level report (TEMP):
                    !   0 - active data used, value ok
                    !   1 - active data used, value not ok
                    !   2 - no active data usable, passive data used, value ok
                    !   3 - no active data usable, passive data used, not ok
                    !   4 - no data at all usable for extrapolation
                    ! - value extrapolated from surface station:
                    !   0 - value ok
                    !   1 - value not ok (not passing quality control)
                    !   2 - value not ok (orography difference too large)
    nvhi_t = 13  ,& ! global x-coord.  \  of the grid point of the reference run
    nvhj_t = 14  ,& ! global y-coord.  /  assigned to the obs
                    !
    ilstidv = 9     ! character length of the station identity on the VOF

  ! VOF (Verification Observation File) body format (for packed info, see ODR)
  ! -----------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    mxvbdm = 10  ,& ! body length of multi-level reports
    mxvbdg = 11  ,& ! body length of GPS IWV reports
    mxvbdu = 11  ,& ! body length of (upper-air), non-SYNOP single-level reports
    mxvbdt = 16  ,& ! body length of SYNOP reports with cloud and precipitation
    mxvbds = 22  ,& ! body length of complete (synoptic) SYNOP reports
    nvbu   =  1  ,& ! u wind component                                [1/10 m/s]
    nvbv   =  2  ,& ! v wind component                                [1/10 m/s]
    nvbt   =  3  ,& ! temperature                                     [1/10 K]
    nvbrh  =  4  ,& ! relative humidity                               [1/10 %]
    nvbp   =  5  ,& ! pressure                                        [Pa]
    nvbz   =  6  ,& ! height                                          [m]
    nvberf =  7  ,& ! bit pattern, bits set for finite obs. errors:
                    !   bit 0: (+ 1): wind
                    !   bit 1: (+ 2): temperature
                    !   bit 2: (+ 4): humidity
                    !   bit 3: (+ 8): pressure / geopotential
                    !   bit 5: (+32): IWV (for GPS reports)
                    !   not used any more: bits 4/5: precip / low cloud cover
    nvbqcf =  8  ,& ! flags of threshold quality ctrl. for wind, temp., humidity
                    !                                        (see ODR: 'nb?qcf')
    nvbmfw =  9  ,& ! main flag word for wind, temp., humidity, pressure/geopot
                    !                                        (see ODR: 'nb?flg')
    nvblid = 10  ,& ! level ID (bit pattern), or pressure code (SYNOP)
                    !                                        (see ODR: 'nb?lid')
    nvbpst = 11  ,& ! (3-hourly) pressure tendency              [10 Pa/3h + 500]
    nvbcl  = 12  ,& ! low cloud cover                                 [octas]
    nvbvis = 13  ,& ! (horizontal) visibility                         [10 m]
    nvbcwg = 14  ,& ! combined cloud and weather group (set of classes)
                    !                                        (see ODR: 'nbscwg')
    nvbwwe = 15  ,& ! weather and ground group               (see ODR: 'nbswwe')
    nvbr12 = 16  ,& ! precipitation amount over 12 hours              [1/10 mm]
                    !   (AOF read: nbspr: precip over 12 or 24 hrs)   [1/10 mm]
    nvbtn  = 17  ,& ! minimum temperature (at 2m during past 12 hrs)  [1/10 K]
    nvbtx  = 18  ,& ! maximum temperature (at 2m during past 12 hrs)  [1/10 K]
    nvbhsw = 19  ,& ! total snow depth                                [m]
                    !   (AOF read: nbstg: min. T-5cm in past 12 hrs)  [1/10 K]
    nvbfg1 = 20  ,& ! max. wind speed of gusts over 1 hour            [m/s]
                    !   (AOF read: nbsfgu: max. wind speed of gusts   [m/s]
    nvbfg6 = 21  ,& ! max. wind speed of gusts over 6 hours           [m/s]
                    !   (AOF read: nbsfme: max. 10-min mean speed)    [m/s]
    nvbrad = 22  ,& ! global radiation, sum over 1 hour               [10 kJ/m2]
    nvbiwa =  1  ,& ! for GPS: adjusted (bias corrected) IWV          [1/100 mm]
    nvbiwv =  2  ,& ! for GPS: integrated water vapour (IWV)          [1/100 mm]
    nvbtzd =  6  ,& ! for GPS: total zenith delay                     [mm]
    nvbzwd = 11     ! for GPS: zenith wet delay                       [mm]

  ! VOF (Verification Observation File) increments (deviations) format
  ! ------------------------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    mxvddm =  5  ,& ! deviation body length of multi-level reports
    mxvdds = 14  ,& ! deviation body length of complete (synoptic) SYNOP reports
    mxvddo =  7  ,& ! deviation body length of SYNOP reports in ODR
    mxvddb =  5  ,& ! deviation body length of non-SYNOP surface-level reports
    mxvddu =  4  ,& ! deviation body length of upper-air single-level reports
    mxvddt =  9  ,& ! deviation body length length of SYNOP reports with cloud
    mxvddg =  2  ,& ! deviation body length of GPS reports
    mxvddh =  1  ,& ! deviation header length of multi-level reports
    nvdu   =  1  ,& ! u wind component                               [1/100 m/s]
    nvdv   =  2  ,& ! v wind component                               [1/100 m/s]
    nvdt   =  3  ,& ! temperature                                    [1/100 K]
    nvdrh  =  4  ,& ! relative humidity                              [1/10 %]
    nvdp   =  5  ,& ! pressure (sfc.) | geopot. (upper-air)       [Pa] | [m2/s2]
    nvmct  =  6  ,& ! total cloud cover              (model value)   [octas]
    nvmcl  =  7  ,& ! low cloud cover                (model value)   [octas]
    nvmvis =  8  ,& ! (horizontal) visibility        (model value)   [10 m]
    nvmpr  =  9  ,& ! precipitation amount           (model value)   [1/10 mm]
    nvdtn  = 10  ,& ! minimum temperature (at 2m during past 12 hrs) [1/10 K]
    nvdtx  = 11  ,& ! maximum temperature (at 2m during past 12 hrs) [1/10 K]
    nvdtg  = 12  ,& ! ground temperature (min. T-5cm in past 12 hrs) [1/10 K]
    nvdfgu = 13  ,& ! max. wind speed of gusts                       [m/s]
    nvdrad = 14  ,& ! global radiation, sum over 1 hour              [10 kJ/m2]
    nvdiq  =  1  ,& ! for GPS: integrated water vapour (increment)   [1/100 mm]
    nvdzpd =  2  ,& ! for GPS: zenith total path delay (increment)   [1/100 mm]
    nvdps  =  1     ! pressure                                       [Pa]

  ! other parameters
  ! ----------------

  REAL    (KIND = wp)        , PARAMETER  :: &
    g_bufr  =     9.807_wp   ,& ! gravity acceleration as used in BUFR / VOF
    c1000r  =     0.001_wp   ,& !
    c100r   =     0.01_wp    ,& !
    c10r    =     0.1_wp     ,& !
    c10     =    10.0_wp     ,& !
    c60     =    60.0_wp     ,& !
    c100    =   100.0_wp     ,& !
    c1000   =  1000.0_wp     ,& !
    c3600   =  3600.0_wp     ,& !
    csmall  = 0.00001_wp        ! less than half a second [hrs]

  INTEGER (KIND = iintegers) , PARAMETER  :: &
    nozz    =  -999    ,& ! missing value for heights in the VOF
    nouv    =  9999    ,& ! missing value for winds in the VOF
    noval   =    -9    ,& ! missing value = -9 for non-negative quantities
    nodev   =  9999       ! missing value for deviations in the VOF

  INTEGER (KIND = iintegers) , PARAMETER  :: &
    itype_calendar = 0    ! for specifying the calendar used in 'get_utc_date'

! Local scalars:
! -------------

  INTEGER (KIND = iintegers) , SAVE :: &
    mboxold = -9999       ! time [min] of last (latest) report which had been
                          ! written during the previous calls of the routine

  LOGICAL                    , SAVE :: &
    lfirst = .TRUE.       ! .TRUE., if routine is called for the first time

  CHARACTER (LEN=28)         ::  &
    yakdat2  ! actual date in the form   wd   dd.mm.yy  hh mm ss UTC

  LOGICAL                    :: &
    lvofoi                ! .TRUE if obs increments are also written to VOF

  INTEGER (KIND = iintegers) ::  &
    nodrln             ,& ! dynamic length of the local ODR output buffer
    ireplen            ,& ! length of the current report
    iob      ,irep     ,& ! (loop) indices for reports
    klev     ,kl  ,kk  ,& ! (loop) indices for levels (or variables)
    icl                ,& ! loop index for characters in a string
    ii                 ,& ! scaled observation increment as integer
    ircv               ,& ! rank of the receiving PE
!   iodrerr            ,& ! ODR obs. error flag code
    ipos     ,indst    ,& ! address  within a buffer
    ioutlen            ,& ! number of elements in local buffer
    ilentot            ,& ! number of elements in collected buffer
    nxrep    ,nrep     ,& ! number of reports to be printed
    nlev               ,& ! number of levels in report
    itmp               ,& ! temporary report index
    nfcast             ,& ! number of forecast runs to compare with observations
    nversa             ,& ! start  of verification period in seconds
    ntrun              ,& ! verification length [hrs] (rounded up)
    mxvdd              ,& ! deviation length of report
    nzday              ,& ! Julian day
    mreptim            ,& ! time [min] of current report
    mboxnow            ,& ! time [min] of last (latest) report written currently
    nzerr    ,nstat    ,& ! error status variable
    mpierr                ! for MPI error code

  REAL    (KIND = wp)        :: &
    tveria             ,& ! start  of verification period in hours ('formatted')
    tverie             ,& ! end    of verification period in hours ('formatted')
!   hverint            ,& ! length of verification period in hours
    endlat_tot         ,& ! / upper right corner
    endlon_tot         ,& ! \ of model domain
    zdphir             ,& ! resolution in reciproke degrees
    trun               ,& ! 'ntrun' minus verification length
    zhr                   ! hour

  LOGICAL                    :: &
    lchange               ! .TRUE., if the order of the report list is changed

  CHARACTER (LEN=ilstidp)    :: &
    ystid                 ! station identity

  CHARACTER (LEN=14)       ::  &
    ygrbt                 ! relevant GRIB date ([yyyymmddhhmmss]:
                          ! year,month,day,hr, min., sec.)
  CHARACTER (LEN=11)       ::  &
    ydescr                ! description of model run type
  CHARACTER (LEN=37)       :: &
    yuvofmt               ! format specifier
  CHARACTER (LEN=20)         :: &
    yroutine              ! name of this subroutine
  CHARACTER (LEN=30)         :: &
    yerrmsg               ! error message


! Local (automatic) arrays:
! -------------------------

  INTEGER (KIND = iintegers) , ALLOCATABLE  :: &
    odrbufp        (:) ,& ! buffer containing ODR data of a single PE
    ioffs          (:) ,& ! list of offsets of the reports in odrbufa
    isortrv        (:) ,& ! sorted index list of reports to be printed
    icriterion     (:)    ! criterion accord. to which the index list is sorted

  INTEGER (KIND = iintegers) , POINTER      :: &
    odrbufa        (:)    ! buffer containing ODR data of a all PEs
!
!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine obs_print_vof
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 0: Preliminaries
!-------------------------------------------------------------------------------

  yroutine   = 'obs_print_vof'
  ircv       = 0
  lvofoi     = (mruntyp >= 0)
  zhr        = REAL (INT( hversta + epsy ), wp)

  nfcast     = 0
  IF (lvofoi)  nfcast = 1

!-------------------------------------------------------------------------------
!  Section 1:  Write local ODR data to buffer
!-------------------------------------------------------------------------------

! dynamic allocation of memory for the local integer buffer
! ---------------------------------------------------------
  nodrln = 0
  DO iob = 1, ntotsg
!   IF (mosghd(iob,nhqcfw) >= 0) nodrln = nodrln + 29
    IF (mosghd(iob,nhqcfw) >= 0) nodrln = nodrln  + mxvhed + ilstidv + mxvbds  &
                                        + nfcast *mxvdds
  ENDDO
  DO iob = 1, ntotml
!   IF (momlhd(iob,nhqcfw) >= 0) nodrln = nodrln  + 15 + 13* momlhd(iob,nhnlev)
    IF (momlhd(iob,nhqcfw) >= 0) nodrln = nodrln  + mxvhed + ilstidv           &
                                                  + mxvbdm *momlhd(iob,nhnlev) &
                                        + nfcast *( mxvddm *momlhd(iob,nhnlev) &
                                                   +mxvddh)
  ENDDO
  DO iob = 1, ntotgp
    IF (mogphd(iob,nhqcfw) >= 0) nodrln = nodrln  + mxvhed + ilstidv + mxvbdg  &
                                        + nfcast *mxvddg
  ENDDO
  DO iob = 1, ntottv
    IF (motvhd(iob,nhqcfw) >= 0) nodrln = nodrln  + mxvhed + ilstidv           &
                                                  + mxvbdm *motvhd(iob,nhnlev) &
                                        + nfcast *( mxvddm *motvhd(iob,nhnlev) &
                                                   +mxvddh)
  ENDDO
  ALLOCATE ( odrbufp (nodrln+1)    , STAT = nzerr )
  odrbufp (:) = noval

! store ODR data to be printed in an intermediate integer buffer
! --------------------------------------------------------------
  ipos    =  0
  nxrep   =  0
  ioutlen =  0

single_level_report: DO iob = 1 , ntotsg
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IF ((mosghd(iob,nhqcfw) >= 0) .AND. (osghed(iob,nhtime) >= hversta)          &
                                .AND. (osghed(iob,nhtime) <= hverend)) THEN
    nxrep    =  nxrep + 1
    ireplen  =  mxvhed + ilstidv + mxvbds + nfcast *mxvddu
    IF ((lvofoi) .AND. (mosghd(iob,nhobtp) /= nairep)                          &
                 .AND. (mosghd(iob,nhobtp) /= nsatob))                         &
      ireplen  =  ireplen + mxvdds - mxvddu
    ioutlen  =  ioutlen + ireplen

!   header of single-level reports
!   ------------------------------
!   odrbufp(ipos+     1) = 29    ! (header length = 15, body length = 14)
    odrbufp(ipos+     1) = ireplen
    odrbufp(ipos+nvhlon) = NINT( osghed(iob,nhilon) * c100 )
    odrbufp(ipos+nvhlat) = NINT( osghed(iob,nhjlat) * c100 )
    odrbufp(ipos+nvhtim) = NINT( (osghed(iob,nhtime) - zhr) * c60 )
    IF (osghed(iob,nhalt) > rmdich) THEN
      odrbufp(ipos+nvhalt) = NINT( osghed(iob,nhalt ) )
    ELSE
      odrbufp(ipos+nvhalt) = nozz
    ENDIF
    odrbufp(ipos+nvhsfc) = NINT( osghed(iob,nhsurf) )
    odrbufp(ipos+nvhobt) = mosghd(iob,nhobtp)
    odrbufp(ipos+nvhcdt) = mosghd(iob,nhcode)
    !   in VOF, temporarily set ACARS still to old value 244 instead of 145
    IF (odrbufp(ipos+nvhcdt) == nacar)  odrbufp(ipos+nvhcdt) = nacar_vof
    !   this (correctly) assumes that in entry 'nhschr', the 2 bits starting
    !   from 'nvpabp' are not set
    odrbufp(ipos+nvhsch) = IOR( mosghd(iob,nhschr), ISHFT( mosghd(iob,nhpass)  &
                                                         , nvpabp ) )
    IF (      (mosghd(iob,nhpass) == 1)                                        &
        .AND. (.NOT. BTEST( mosghd(iob,nhflag), FL_MERGE )))                   &
      mosghd (iob,nhflag) = IBSET ( mosghd(iob,nhflag), FL_MERGE )
    odrbufp(ipos+nvhflg) = mosghd(iob,nhflag)
    odrbufp(ipos+nvhpas) = mosghd(iob,nhpass)
    odrbufp(ipos+nvhqcp) = mosghd(iob,nhqcfw)
    IF (BTEST( mosgbd(iob,nbsqcf), nvrz ))                                     &
      odrbufp(ipos+nvhqcp) =      odrbufp(ipos+nvhqcp)/ 2   * 2 + 1
    IF (BTEST( mosgbd(iob,nbsflg), nvfzbp+nvfbps(4) ))                         &
      odrbufp(ipos+nvhqcp) = MOD( odrbufp(ipos+nvhqcp), 2 ) + 2
    odrbufp(ipos+nvhi_t) = mosghd(iob,nhitot)
    odrbufp(ipos+nvhj_t) = mosghd(iob,nhjtot)
    ipos  =  ipos + mxvhed

    DO icl = 1 , ilstidv
      odrbufp(ipos+icl) = ICHAR( yosghd(iob) (icl:icl) )
    ENDDO
    ipos  =  ipos + ilstidv

!   body of single-level reports
!   ----------------------------
    IF ((osgbdy(iob,nbsu) > rmdich) .AND. (osgbdy(iob,nbsv) > rmdich)) THEN
!     roblat = osghed (iob,nhjlat)
!     roblon = osghed (iob,nhilon)
!     zurot  = osgbdy (iob,nbsu  )
!     zvrot  = osgbdy (iob,nbsv  )
!     CALL uvrot2uv ( zurot, zvrot, roblat, rloblon, pollat, pollon, zuu, zvv )
!     =============
!     CALL uv2df ( zuu, zvv, zdd, zff )
!     ==========
!     odrbufp(ipos+16) = NINT( zff * c100 )
!     odrbufp(ipos+17) = NINT( zdd * c10  )
      odrbufp(ipos+nvbu  ) = NINT( osgbdy(iob,nbsu  ) * c10 )
      odrbufp(ipos+nvbv  ) = NINT( osgbdy(iob,nbsv  ) * c10 )
    ELSE
      odrbufp(ipos+nvbu  ) = nouv
      odrbufp(ipos+nvbv  ) = nouv
    ENDIF
    IF (osgbdy(iob,nbst  ) > rmdich)                                           &
      odrbufp(ipos+nvbt  ) = NINT( osgbdy(iob,nbst  ) * c10 )
    IF (osgbdy(iob,nbsrh ) > rmdich)                                           &
      odrbufp(ipos+nvbrh ) = NINT( osgbdy(iob,nbsrh ) * c1000 )
    IF (osgbdy(iob,nbsp  ) > rmdich)                                           &
      odrbufp(ipos+nvbp  ) = NINT( osgbdy(iob,nbsp  ) )
    IF (osgbdy(iob,nbsz  ) > rmdich) THEN
      odrbufp(ipos+nvbz  ) = NINT( osgbdy(iob,nbsz  ) )
    ELSE
      odrbufp(ipos+nvbz  ) = nozz
    ENDIF
    IF ((osgbdy(iob,nbsr12) > epsy) .AND. (osgbdy(iob,nbsr12) < 0.09_wp)) THEN
      odrbufp(ipos+nvbr12) = - 1
    ELSEIF (osgbdy(iob,nbsr12) > rmdich) THEN
      odrbufp(ipos+nvbr12) = NINT( osgbdy(iob,nbsr12) * c10 )
    ENDIF
    IF (osgbdy(iob,nbscl ) > rmdich)                                           &
      odrbufp(ipos+nvbcl ) = NINT( osgbdy(iob,nbscl ) )
    IF (osgbdy(iob,nbsvis) > rmdich)                                           &
      odrbufp(ipos+nvbvis) = NINT( osgbdy(iob,nbsvis) * c10r )
    IF (osgbdy(iob,nbstn ) > rmdich)                                           &
      odrbufp(ipos+nvbtn ) = NINT( osgbdy(iob,nbstn ) * c10 )
    IF (osgbdy(iob,nbstx ) > rmdich)                                           &
      odrbufp(ipos+nvbtx ) = NINT( osgbdy(iob,nbstx ) * c10 )
    IF (osgbdy(iob,nbshsw) > rmdich)                                           &
      odrbufp(ipos+nvbhsw) = NINT( osgbdy(iob,nbshsw) * c10 )
    IF (osgbdy(iob,nbsfg1) > rmdich)                                           &
      odrbufp(ipos+nvbfg1) = NINT( osgbdy(iob,nbsfg1) )
    IF (osgbdy(iob,nbsfg6) > rmdich)                                           &
      odrbufp(ipos+nvbfg6) = NINT( osgbdy(iob,nbsfg6) )
    IF (osgbdy(iob,nbsrad) > rmdich)                                           &
      odrbufp(ipos+nvbrad) = NINT( osgbdy(iob,nbsrad) * 0.0001_wp )
    IF (osgbdy(iob,nbspst) > rmdich)                                           &
      odrbufp(ipos+nvbpst) = NINT( osgbdy(iob,nbspst) * c10r ) + 500
!   iodrerr = 0
!   IF (osgbdy(iob,nbsuer) > rmdich) iodrerr = iodrerr + 1
!   IF (osgbdy(iob,nbster) > rmdich) iodrerr = iodrerr + 2
!   IF (osgbdy(iob,nbsqer) > rmdich) iodrerr = iodrerr + 4
!   IF (osgbdy(iob,nbszer) > rmdich) iodrerr = iodrerr + 8
!   odrbufp(ipos+nvberf) = iodrerr
    odrbufp(ipos+nvberf) = IBITS( mosgbd(iob,nbserr), 0, nvrz+1 )
    ! omit the temporary flags with bit position >= nvrtmp
    odrbufp(ipos+nvbqcf) = IBITS( mosgbd(iob,nbsqcf), 0, nvrtmp )
    odrbufp(ipos+nvbmfw) = mosgbd(iob,nbsflg)
    odrbufp(ipos+nvblid) = mosgbd(iob,nbslid)
    odrbufp(ipos+nvbcwg) = mosgbd(iob,nbscwg)
    odrbufp(ipos+nvbwwe) = mosgbd(iob,nbswwe)

    IF (BTEST( mosgbd(iob,nbsqcf), nvrzbc ))                                   &
      odrbufp(ipos+nvbqcf) = IBSET ( odrbufp(ipos+nvbqcf), nvrz )
    IF (mosghd(iob,nhqcfw)/2 == 1)                                             &
      odrbufp(ipos+nvbmfw) = IBSET ( odrbufp(ipos+nvbmfw), nvfzbp+nvfbps(4) )
    ipos  =  ipos + mxvbds

!   observation increments for single-level reports
!   -----------------------------------------------
    IF (lvofoi) THEN
      odrbufp(ipos+nvdu  ) = nodev
      odrbufp(ipos+nvdv  ) = nodev
      odrbufp(ipos+nvdt  ) = nodev
      odrbufp(ipos+nvdrh ) = nodev
      IF (      (ssgbdy(nso_u ,iob) > rmdich)                                  &
          .AND. (ssgbdy(nso_v ,iob) > rmdich)) THEN
        ii = NINT( (osgbdy(iob,nbsu ) - ssgbdy(nso_u ,iob)) * c100 )
        IF (ABS( ii ) < ABS( nodev ))  odrbufp(ipos+nvdu ) = ii
        ii = NINT( (osgbdy(iob,nbsv ) - ssgbdy(nso_v ,iob)) * c100 )
        IF (ABS( ii ) < ABS( nodev ))  odrbufp(ipos+nvdv ) = ii
      ENDIF
      IF (ssgbdy(nso_t ,iob) > rmdich) THEN
        ii = NINT( (osgbdy(iob,nbst ) - ssgbdy(nso_t ,iob)) * c100 )
        IF (ABS( ii ) < ABS( nodev ))  odrbufp(ipos+nvdt ) = ii
      ENDIF
      IF (ssgbdy(nso_rh,iob) > rmdich) THEN
        ii = NINT( (osgbdy(iob,nbsrh) - ssgbdy(nso_rh,iob)) * c1000 )
        IF (ABS( ii ) < ABS( nodev ))  odrbufp(ipos+nvdrh) = ii
      ENDIF
      IF (      (mosghd(iob,nhobtp) /= nairep)                                 &
          .AND. (mosghd(iob,nhobtp) /= nsatob)) THEN
        odrbufp(ipos+nvdp  ) = nodev
        odrbufp(ipos+nvdtn ) = nodev
        odrbufp(ipos+nvdtx ) = nodev
        odrbufp(ipos+nvdtg ) = nodev
        odrbufp(ipos+nvdfgu) = nodev
        odrbufp(ipos+nvdrad) = nodev
        IF (ssgbdy(nso_p ,iob) > rmdich) THEN
          ii = NINT( osgbdy(iob,nbsp ) - ssgbdy(nso_p ,iob) )
          IF (ABS( ii ) < ABS( nodev ))  odrbufp(ipos+nvdp ) = ii
        ENDIF
        IF (ssgbdy(nso_ct ,iob) > rmdich)                                      &
          odrbufp(ipos+nvmct ) = NINT( ssgbdy(nso_ct ,iob) * c100 )
        IF (ssgbdy(nso_cl ,iob) > rmdich)                                      &
          odrbufp(ipos+nvmcl ) = NINT( ssgbdy(nso_cl ,iob) * c100 )
        ipos  =  ipos + mxvdds
      ELSE
        ipos  =  ipos + mxvddu
      ENDIF
    ENDIF

  ENDIF

! set the flag for having the printing of the current report done
! ---------------------------------------------------------------

  IF ((MOD(mveripr,2) == 0) .AND. (mosghd(iob,nhqcfw) >= 0))                   &
    mosghd (iob,nhqcfw) = - 99

ENDDO single_level_report
!~~~~~~~~~~~~~~~~~~~~~~~~


multi_level_report: DO iob = 1 , ntotml
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IF ((momlhd(iob,nhqcfw) >= 0) .AND. (omlhed(iob,nhtime) >= hversta)          &
                                .AND. (omlhed(iob,nhtime) <= hverend)) THEN
    nxrep    =  nxrep + 1
    ireplen  =  mxvhed + ilstidv + (          mxvbdm                           &
                                    + nfcast *mxvddm)* momlhd(iob,nhnlev)      &
                                    + nfcast *mxvddh
    ioutlen  =  ioutlen + ireplen

!   header of multi-level report (length = 15)
!   ------------------------------------------
    indst  =  ipos +1
!   odrbufp(ipos+ 1) = 15 + 13* momlhd(iob,nhnlev)
    odrbufp(ipos+     1) = mxvhed + ilstidv + mxvbdm* momlhd(iob,nhnlev)
    odrbufp(ipos+nvhlon) = NINT( omlhed(iob,nhilon) * c100 )
    odrbufp(ipos+nvhlat) = NINT( omlhed(iob,nhjlat) * c100 )
    odrbufp(ipos+nvhtim) = NINT( (omlhed(iob,nhtime) - zhr) * c60 )
    IF (omlhed(iob,nhalt) > rmdich) THEN
      odrbufp(ipos+nvhalt) = NINT( omlhed(iob,nhalt ) )
    ELSE
      odrbufp(ipos+nvhalt) = nozz
    ENDIF
    odrbufp(ipos+nvhsfc) = NINT( omlhed(iob,nhsurf) )
    odrbufp(ipos+nvhobt) = - momlhd(iob,nhobtp)
    odrbufp(ipos+nvhcdt) = momlhd(iob,nhcode)
    !   in VOF, temporarily set ACARS still to old value 244 instead of 145
    IF (odrbufp(ipos+nvhcdt) == nacar)  odrbufp(ipos+nvhcdt) = nacar_vof
    odrbufp(ipos+nvhsch) = IOR( momlhd(iob,nhschr), ISHFT( momlhd(iob,nhpass)  &
                                                         , nvpabp ) )
    odrbufp(ipos+nvhflg) = momlhd(iob,nhflag)
    odrbufp(ipos+nvhpas) = momlhd(iob,nhnlev) * (1 - momlhd(iob,nhpass))
    odrbufp(ipos+nvhqcp) = momlhd(iob,nhqcfw)
    IF (lvofoi) THEN
      IF (ABS( dmlhed(nso_ps,iob) ) > (REAL(ABS(nodev), wp)-c05)) odrbufp(ipos+nvhqcp) = 4
    ENDIF
    odrbufp(ipos+nvhi_t) = momlhd(iob,nhitot)
    odrbufp(ipos+nvhj_t) = momlhd(iob,nhjtot)
    ipos  =  ipos + mxvhed

    DO icl = 1 , ilstidv
      odrbufp(ipos+icl) = ICHAR( yomlhd(iob) (icl:icl) )
    ENDDO
    ipos  =  ipos + ilstidv

!   body of multi-level report (length = 13* nlevel)
!   ------------------------------------------------

    Level: DO klev = 1 , momlhd(iob,nhnlev)

      IF (      (omlbdy(iob,klev,nbtu) > rmdich)                               &
          .AND. (omlbdy(iob,klev,nbtv) > rmdich)) THEN
        odrbufp(ipos+nvbu  ) = NINT( omlbdy(iob,klev,nbtu  ) * c10 )
        odrbufp(ipos+nvbv  ) = NINT( omlbdy(iob,klev,nbtv  ) * c10 )
      ELSE
        odrbufp(ipos+nvbu  ) = nouv
        odrbufp(ipos+nvbv  ) = nouv
      ENDIF
      IF (omlbdy(iob,klev,nbtt  ) > rmdich)                                    &
        odrbufp(ipos+nvbt  ) = NINT( omlbdy(iob,klev,nbtt  ) * c10 )
      IF (omlbdy(iob,klev,nbtrh ) > rmdich)                                    &
        odrbufp(ipos+nvbrh ) = NINT( omlbdy(iob,klev,nbtrh ) * c1000 )
      IF (omlbdy(iob,klev,nbtp  ) > rmdich)                                    &
        odrbufp(ipos+nvbp  ) = NINT( omlbdy(iob,klev,nbtp  ) )
      IF (omlbdy(iob,klev,nbtz  ) > rmdich) THEN
        odrbufp(ipos+nvbz  ) = NINT( omlbdy(iob,klev,nbtz  ) )
      ELSE
        odrbufp(ipos+nvbz  ) = nozz
      ENDIF
      odrbufp(ipos+nvberf) = IBITS( momlbd(iob,klev,nbterr), 0, nvrw+1 )
      !   omit the temporary flags with bit position >= nvrtmp
      odrbufp(ipos+nvbqcf) = IBITS( momlbd(iob,klev,nbtqcf), 0, nvrtmp )
      odrbufp(ipos+nvbqcf) = momlbd(iob,klev,nbtqcf)
      odrbufp(ipos+nvbmfw) = momlbd(iob,klev,nbtflg)
      odrbufp(ipos+nvblid) = momlbd(iob,klev,nbtlid)
      !   w-prof, RASS: unset dataset flag, if it refers to 'w', not 'u,v', 'T'
      IF ((omlbdy(iob,klev,nbtw) > rmdich) .AND. (odrbufp(ipos+nvbu) == nouv)) &
        CALL MVBITS( 0, 0, nvfboc(1), odrbufp(ipos+nvbmfw), nvfubp + nvfbps(1) )
      IF ((omlbdy(iob,klev,nbtw) > rmdich) .AND. (odrbufp(ipos+nvbt) ==noval)) &
        CALL MVBITS( 0, 0, nvfboc(1), odrbufp(ipos+nvbmfw), nvftbp + nvfbps(1) )
      IF (BTEST( momlbd(iob,klev,nbtqcf), nvrqbc ))                            &
        odrbufp(ipos+nvbqcf) = IBSET ( odrbufp(ipos+nvbqcf), nvrq )
!     IF (MOD(mveripr,2) == 0)  momlbd(iob,klev,nbtqcf) = -2
      ipos  =  ipos + mxvbdm
    ENDDO Level

!   observation increments for multi-level reports
!   ----------------------------------------------
    IF (lvofoi) THEN
      DO klev = 1 , momlhd(iob,nhnlev)
        odrbufp(ipos+nvdu  ) = nodev
        odrbufp(ipos+nvdv  ) = nodev
        odrbufp(ipos+nvdt  ) = nodev
        odrbufp(ipos+nvdrh ) = nodev
        odrbufp(ipos+nvdp  ) = nodev
        IF (      (smlbdy(nso_u ,klev,iob) > rmdich)                           &
            .AND. (smlbdy(nso_v ,klev,iob) > rmdich)) THEN
          ii = NINT( (omlbdy(iob,klev,nbtu ) - smlbdy(nso_u ,klev,iob)) * c100 )
          IF (ABS( ii ) < ABS( nodev ))  odrbufp(ipos+nvdu ) = ii
          ii = NINT( (omlbdy(iob,klev,nbtv ) - smlbdy(nso_v ,klev,iob)) * c100 )
          IF (ABS( ii ) < ABS( nodev ))  odrbufp(ipos+nvdv ) = ii
        ENDIF
        IF (smlbdy(nso_t ,klev,iob) > rmdich) THEN
          ii = NINT( (omlbdy(iob,klev,nbtt ) - smlbdy(nso_t ,klev,iob)) * c100 )
          IF (ABS( ii ) < ABS( nodev ))  odrbufp(ipos+nvdt ) = ii
        ENDIF
        IF (smlbdy(nso_rh,klev,iob) > rmdich) THEN
          ii = NINT( (omlbdy(iob,klev,nbtrh) - smlbdy(nso_rh,klev,iob))* c1000 )
          IF (ABS( ii ) < ABS( nodev ))  odrbufp(ipos+nvdrh) = ii
        ENDIF
        IF (smlbdy(nso_p ,klev,iob) > rmdich) THEN
          ii = NINT( (omlbdy(iob,klev,nbtz ) - smlbdy(nso_p ,klev,iob))* g_bufr)
          IF (ABS( ii ) < ABS( nodev ))  odrbufp(ipos+nvdp ) = ii
        ENDIF
        ipos  =  ipos + mxvddm
      ENDDO
      odrbufp(ipos+nvdps ) = nodev
      IF (ABS( dmlhed(nso_ps,iob) ) < (REAL(ABS( nodev ),wp)-c05))                      &
        odrbufp(ipos+nvdps ) = NINT( dmlhed(nso_ps,iob) )
      ipos  =  ipos + mxvddh
    ENDIF

! write report length to header
    odrbufp(indst) = ipos-indst+1

  ENDIF

! set the flag for having the printing of the current report done
! ---------------------------------------------------------------

  IF ((MOD(mveripr,2) == 0) .AND. (momlhd(iob,nhqcfw) >= 0))                   &
    momlhd (iob,nhqcfw) = - 99

ENDDO multi_level_report
!~~~~~~~~~~~~~~~~~~~~~~~


sat_retrievals: DO iob = 1 , ntottv
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IF ((motvhd(iob,nhqcfw) >= 0) .AND. (otvhed(iob,nhtime) >= hversta)          &
                                .AND. (otvhed(iob,nhtime) <= hverend)) THEN
    nxrep    =  nxrep + 1
    ireplen  =  mxvhed + ilstidv + (          mxvbdm                           &
                                    + nfcast *mxvddm)* motvhd(iob,nhnlev)      &
                                    + nfcast *mxvddh
    ioutlen  =  ioutlen + ireplen

!   header of sat retrievals, like multi-level report (length = 15)
!   ---------------------------------------------------------------
    indst  =  ipos +1
!   odrbufp(ipos+ 1) = 15 + 13* motvhd(iob,nhnlev)
    odrbufp(ipos+     1) = mxvhed + ilstidv + mxvbdm* motvhd(iob,nhnlev)
    odrbufp(ipos+nvhlon) = NINT( otvhed(iob,nhilon) * c100 )
    odrbufp(ipos+nvhlat) = NINT( otvhed(iob,nhjlat) * c100 )
    odrbufp(ipos+nvhtim) = NINT( (otvhed(iob,nhtime) - zhr) * c60 )
    IF (otvhed(iob,nhalt) > rmdich) THEN
      odrbufp(ipos+nvhalt) = NINT( otvhed(iob,nhalt ) )
    ELSE
      odrbufp(ipos+nvhalt) = nozz
    ENDIF
    odrbufp(ipos+nvhsfc) = NINT( otvhed(iob,nhsurf) )
    odrbufp(ipos+nvhobt) = - motvhd(iob,nhobtp)
    odrbufp(ipos+nvhcdt) = motvhd(iob,nhcode)
    odrbufp(ipos+nvhsch) = IOR( motvhd(iob,nhschr), ISHFT( motvhd(iob,nhpass)  &
                                                         , nvpabp ) )
    odrbufp(ipos+nvhflg) = motvhd(iob,nhflag)
    odrbufp(ipos+nvhpas) = motvhd(iob,nhnlev) * (1 - motvhd(iob,nhpass))
    odrbufp(ipos+nvhqcp) = motvhd(iob,nhqcfw)
    odrbufp(ipos+nvhi_t) = motvhd(iob,nhitot)
    odrbufp(ipos+nvhj_t) = motvhd(iob,nhjtot)
    ipos  =  ipos + mxvhed

    DO icl = 1 , ilstidv
      odrbufp(ipos+icl) = ICHAR( yotvhd(iob) (icl:icl) )
    ENDDO
    ipos  =  ipos + ilstidv

!   body of sat retrieval, like multi-level report (length = 13* nlevel)
!   --------------------------------------------------------------------

    Sat_level: DO klev = 1 , motvhd(iob,nhnlev)

      odrbufp(ipos+nvbu  ) = nouv
      odrbufp(ipos+nvbv  ) = nouv
      IF (otvbdy(iob,klev,nbvt  ) > rmdich)                                    &
        odrbufp(ipos+nvbt  ) = NINT( otvbdy(iob,klev,nbvt  ) * c10 )
      IF (otvbdy(iob,klev,nbvrh ) > rmdich)                                    &
        odrbufp(ipos+nvbrh ) = NINT( otvbdy(iob,klev,nbvrh ) * c1000 )
      IF (otvbdy(iob,klev,nbvp  ) > rmdich)                                    &
        odrbufp(ipos+nvbp  ) = NINT( otvbdy(iob,klev,nbvp  ) )
      odrbufp(ipos+nvbz  ) = nozz
!     iodrerr = 0
!!    IF (      (otvbdy(iob,klev,nbvt  ) > rmdich)                             &
!!        .AND. (MOD( motvhd(iob,nhvqcf),4 ) < 2)) iodrerr = iodrerr + 2
!     IF (otvbdy(iob,klev,nbvt  ) > rmdich) iodrerr = iodrerr + 2
!     IF (otvbdy(iob,klev,nbvrh ) > rmdich) iodrerr = iodrerr + 4
!     odrbufp(ipos+nvberf) = iodrerr
      odrbufp(ipos+nvberf) = motvbd(iob,klev,nbverr)
      odrbufp(ipos+nvbqcf) = motvhd(iob,     nhvqcf)
      odrbufp(ipos+nvbmfw) = motvbd(iob,klev,nbvflg)
      odrbufp(ipos+nvblid) = 0
      ipos  =  ipos + mxvbdm
    ENDDO Sat_level

!   observation increments for sat retrievals, like multi-level reports
!   -------------------------------------------------------------------
    IF (lvofoi) THEN
      DO klev = 1 , motvhd(iob,nhnlev)
        odrbufp(ipos+nvdu  ) = nodev
        odrbufp(ipos+nvdv  ) = nodev
        odrbufp(ipos+nvdt  ) = nodev
        odrbufp(ipos+nvdrh ) = nodev
        odrbufp(ipos+nvdp  ) = nodev
        IF (stvbdy(nso_t ,klev,iob) > rmdich) THEN
          ii = NINT( (otvbdy(iob,klev,nbtt ) - stvbdy(nso_t ,klev,iob)) * c100 )
          IF (ABS( ii ) < ABS( nodev ))  odrbufp(ipos+nvdt ) = ii
        ENDIF
        IF (stvbdy(nso_rh,klev,iob) > rmdich) THEN
          ii = NINT( (otvbdy(iob,klev,nbtrh) - stvbdy(nso_rh,klev,iob))* c1000 )
          IF (ABS( ii ) < ABS( nodev ))  odrbufp(ipos+nvdrh) = ii
        ENDIF
        ipos  =  ipos + mxvddm
      ENDDO
      odrbufp(ipos+nvdps ) = nodev
      ipos  =  ipos + mxvddh
    ENDIF

! write report length to header
    odrbufp(indst) = ipos-indst+1

  ENDIF

! set the flag for having the printing of the current report done
! ---------------------------------------------------------------

  IF ((MOD(mveripr,2) == 0) .AND. (motvhd(iob,nhqcfw) >= 0))                   &
    motvhd (iob,nhqcfw) = - 99

ENDDO sat_retrievals
!~~~~~~~~~~~~~~~~~~~


gps_report: DO iob = 1 , ntotgp
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IF ((mogphd(iob,nhqcfw) >= 0) .AND. (ogphed(iob,nhtime) >= hversta)          &
                                .AND. (ogphed(iob,nhtime) <= hverend)) THEN
    nxrep    =  nxrep + 1
    ireplen  =  mxvhed + ilstidv + mxvbdg + nfcast *mxvddg
    ioutlen  =  ioutlen + ireplen

!   header of GPS reports
!   ---------------------
!   odrbufp(ipos+     1) = 29    ! (header length = 15, body length = 14)
    odrbufp(ipos+     1) = ireplen
    odrbufp(ipos+nvhlon) = NINT( ogphed(iob,nhilon) * c100 )
    odrbufp(ipos+nvhlat) = NINT( ogphed(iob,nhjlat) * c100 )
    odrbufp(ipos+nvhtim) = NINT( (ogphed(iob,nhtime) - zhr) * c60 )
    IF (ogphed(iob,nhalt) > rmdich) THEN
      odrbufp(ipos+nvhalt) = NINT( ogphed(iob,nhalt ) )
    ELSE
      odrbufp(ipos+nvhalt) = nozz
    ENDIF
    odrbufp(ipos+nvhsfc) = NINT( ogphed(iob,nhsurf) )
    odrbufp(ipos+nvhobt) = mogphd(iob,nhobtp)
    odrbufp(ipos+nvhcdt) = mogphd(iob,nhcode)
    odrbufp(ipos+nvhsch) = IOR( mogphd(iob,nhschr), ISHFT( mogphd(iob,nhpass)  &
                                                         , nvpabp ) )
    odrbufp(ipos+nvhflg) = mogphd(iob,nhflag)
    odrbufp(ipos+nvhpas) = mogphd(iob,nhpass)
    odrbufp(ipos+nvhqcp) = mogphd(iob,nhqcfw)
    IF (BTEST( mogpbd(iob,nbgflg), nvfzbp+nvfbps(4) ))                         &
      odrbufp(ipos+nvhqcp) = MOD( odrbufp(ipos+nvhqcp), 2 ) + 2
    odrbufp(ipos+nvhi_t) = mogphd(iob,nhitot)
    odrbufp(ipos+nvhj_t) = mogphd(iob,nhjtot)
    ipos  =  ipos + mxvhed

    DO icl = 1 , ilstidv
      odrbufp(ipos+icl) = ICHAR( yogphd(iob) (icl:icl) )
    ENDDO
    ipos  =  ipos + ilstidv

!   body of GPS reports
!   -------------------
    IF (ogpbdy(iob,nbgiwa) > rmdich) THEN
      odrbufp(ipos+nvbiwa) = NINT( ogpbdy(iob,nbgiwa) * c100 )
    ELSE
      odrbufp(ipos+nvbiwa) = nouv
    ENDIF
    IF (ogpbdy(iob,nbgiwv) > rmdich) THEN
      odrbufp(ipos+nvbiwv) = NINT( ogpbdy(iob,nbgiwv) * c100 )
    ELSE
      odrbufp(ipos+nvbiwv) = nouv
    ENDIF
    IF (ogpbdy(iob,nbgt  ) > rmdich)                                           &
      odrbufp(ipos+nvbt  ) = NINT( ogpbdy(iob,nbgt  ) * c10 )
    IF (ogpbdy(iob,nbgrh ) > rmdich)                                           &
      odrbufp(ipos+nvbrh ) = NINT( ogpbdy(iob,nbgrh ) * c1000 )
    IF (ogpbdy(iob,nbgp  ) > rmdich)                                           &
      odrbufp(ipos+nvbp  ) = NINT( ogpbdy(iob,nbgp  ) )
    IF (ogpbdy(iob,nbgzpd) > rmdich) THEN
      odrbufp(ipos+nvbtzd) = NINT( ogpbdy(iob,nbgzpd)       )
    ELSE
      odrbufp(ipos+nvbtzd) = nozz
    ENDIF
    !   also set 'ZPD' flag to active
    odrbufp(ipos+nvberf) = IBITS( mogpbd(iob,nbgerr), 0, nvrzpd+1 )
    ! omit the temporary flags with bit position >= nvrtmp
    odrbufp(ipos+nvbqcf) = IBITS( mogpbd(iob,nbgqcf), 0, nvrtmp )
    odrbufp(ipos+nvbmfw) = mogpbd(iob,nbgflg)
    odrbufp(ipos+nvblid) = mogpbd(iob,nbglid)
    IF (ogpbdy(iob,nbgzwd) > rmdich)                                           &
      odrbufp(ipos+nvbzwd) = NINT( ogpbdy(iob,nbgzwd)       )
    IF (BTEST( mogpbd(iob,nbgqcf), nvrqbc ))                                   &
      odrbufp(ipos+nvbqcf) = IBSET ( odrbufp(ipos+nvbqcf), nvriwv )
    ipos  =  ipos + mxvbdg

!   observation increments for GPS reports
!   --------------------------------------
    IF (lvofoi) THEN
      odrbufp(ipos+nvdiq ) = nodev
      odrbufp(ipos+nvdzpd) = nodev
!     IF (ABS( dgpbdy(nso_iq,iob) ) < (ABS( nodev )-c05)*c100r)                &
!       odrbufp(ipos+nvdiq ) = NINT( dgpbdy(nso_iq ,iob) * c100 )
      IF (sgpbdy(nso_iq,iob) > rmdich) THEN
        ii = NINT( (ogpbdy(iob,nbgiwa) - sgpbdy(nso_iq,iob)) * c100 )
        IF (ABS( ii ) < ABS( nodev ))  odrbufp(ipos+nvdiq) = ii
      ENDIF
      IF (sgpbdy(nso_zpd,iob) > rmdich) THEN
        ii = NINT( (ogpbdy(iob,nbgzpd) - sgpbdy(nso_zpd,iob)) * c100 )
        IF (ABS( ii ) < ABS( nodev ))  odrbufp(ipos+nvdzpd) = ii
      ENDIF
      ipos  =  ipos + mxvddg
    ENDIF

  ENDIF

! set the flag for having the printing of the current report done
! ---------------------------------------------------------------

  IF ((MOD(mveripr,2) == 0) .AND. (mogphd(iob,nhqcfw) >= 0))                   &
    mogphd (iob,nhqcfw) = - 99

ENDDO gps_report
!~~~~~~~~~~~~~~~


!-------------------------------------------------------------------------------
!  Section 2:  Communication:  Gather VOF data at root node
!-------------------------------------------------------------------------------

  ! note that 'odrbufa' is allocated in obs_gather_buffer

  CALL obs_gather_buffer ( ioutlen, nodrln+1, odrbufp, ilentot, odrbufa, ircv  &
                         , num_compute, my_cart_id, icomm_cart, imp_integers )
! ======================

  DEALLOCATE ( odrbufp   , STAT = nzerr )

  ! gather the number of reports from each node

  IF (num_compute > 1)  CALL global_values ( nxrep, 1, 'SUM', imp_integers     &
                                           , icomm_cart, ircv, yerrmsg, mpierr )
!                       ==================

!-------------------------------------------------------------------------------
!  Section 3:  Make a sorted list of the reports
!              (to make the output independent from the domain decomposition)
!-------------------------------------------------------------------------------

  IF (my_cart_id == ircv) THEN
    ALLOCATE ( isortrv (nxrep + 1) , STAT = nzerr )
    ALLOCATE ( ioffs   (nxrep + 1) , STAT = nzerr )

! Make an (unsorted) list of the offsets of the reports in 'odrbufa',
! initialise the sorting index list, and determine the number of reports
    nrep = 0
    ipos = 0
    Report_position: DO irep = 1 , nxrep + 1
      IF (odrbufa(ipos+1) > 0) THEN
        nrep = nrep + 1
        isortrv (nrep) = nrep
        ioffs   (nrep) = ipos
        ipos = ipos + odrbufa(ipos+1)
      ELSE
        EXIT Report_position
      ENDIF
    ENDDO Report_position
    IF (nrep /= nxrep) PRINT '(''WARNING print_obs_verif: nrep =='',I5         &
                             &,'' /= nxrep =='',I5)' ,  nrep, nxrep

! sort index list
    IF (nrep >= 1) THEN
      ALLOCATE ( icriterion (nxrep + 1) , STAT = nzerr )
      DO irep = 1 , nrep
! sorting criterion: first observation time (min), then longitude, then latitude
!                    then multi-level aircraft report yes/no
        icriterion(irep) = odrbufa(ioffs(irep)+nvhtim) *2*(je_tot+1)*(ie_tot+1)&
                         + odrbufa(ioffs(irep)+nvhi_t) *2*(je_tot+1)           &
                         + odrbufa(ioffs(irep)+nvhj_t) *2                      &
                         + MAX( 1-ABS( odrbufa(ioffs(irep)+nvhobt)+nairep ) ,0 )
      ENDDO
! quick-sort reports according to criterion

      CALL sortrx ( nrep, icriterion, isortrv )
!     ===========

! bubble-sort co-located multi-level aircraft reports
      lchange = .TRUE.
      DO WHILE (lchange)
        lchange = .FALSE.
        DO irep = 1 , nrep-1
! find co-located multi-level aircraft reports
          IF (      (icriterion(isortrv(irep)) == icriterion(isortrv(irep+1))) &
              .AND. (odrbufa(ioffs(isortrv(irep  ))+nvhobt) == - nairep)       &
              .AND. (odrbufa(ioffs(isortrv(irep+1))+nvhobt) == - nairep)) THEN
! sort according to report length
            IF (  odrbufa (ioffs(isortrv(irep  ))+ 1)                          &
                < odrbufa (ioffs(isortrv(irep+1))+ 1)) THEN
              itmp             = isortrv (irep+1)
              isortrv (irep+1) = isortrv (irep)
              isortrv (irep)   = itmp
              lchange = .TRUE.
            ELSEIF (   odrbufa (ioffs(isortrv(irep  ))+1)                      &
                    == odrbufa (ioffs(isortrv(irep+1))+1)) THEN
! for equal report length, sort according to succeeding variables
              Sort_variables: DO kl = 2 , odrbufa(ioffs(isortrv(irep))+1)
                IF (  odrbufa(ioffs(isortrv(irep  ))+kl)                       &
                    < odrbufa(ioffs(isortrv(irep+1))+kl)) THEN
                  itmp             = isortrv (irep+1)
                  isortrv (irep+1) = isortrv (irep)
                  isortrv (irep)   = itmp
                  lchange = .TRUE.
                  EXIT Sort_variables
                ENDIF
                IF (  odrbufa(ioffs(isortrv(irep  ))+kl)                       &
                    > odrbufa(ioffs(isortrv(irep+1))+kl))                      &
                  EXIT Sort_variables
              ENDDO Sort_variables
            ENDIF 
          ENDIF
        ENDDO
      ENDDO
      DEALLOCATE ( icriterion , STAT = nzerr )
    ENDIF

!-------------------------------------------------------------------------------
!  Section 4:  Print out gathered ODR data  
!-------------------------------------------------------------------------------

! open VOF
! --------

    IF ((lfirst) .OR. (lastpr) .OR. (nrep >= 1)) THEN
      OPEN (nuverif, FILE=yuverif ,FORM='FORMATTED',STATUS='UNKNOWN'           &
                                  ,POSITION='APPEND',IOSTAT=nstat )
      IF (nstat /= 0) THEN
        yerrmsg = 'OPENING OF FILE yuverif FAILED'
        CALL model_abort (my_cart_id, 1009, yerrmsg, yroutine)
      ENDIF
    ENDIF

! print VOF file header
! ---------------------

    IF (lfirst) THEN
!     nveria  = NINT( hversta / dtdeh )
      nversa  = NINT( hversta * c3600 )

!     CALL get_utc_date (nversa, ydate_ref, c1, itype_calendar, yakdat1,       &
!                        yakdat2, nzday, zhr)
!     =================

      tveria  = rmod( hversta , c1 ) - epsy
      tverie  = hverend - (hversta - tveria)
      tveria  = tveria - epsy + 0.00005_wp
      tverie  = tverie + epsy - 0.00005_wp
      endlat_tot  =  startlat_tot + (je_tot-1) *dlat
      endlon_tot  =  startlon_tot + (ie_tot-1) *dlon
      ! The longitude values have to be limited to the range (-180.0,+180.0)
      IF (endlon_tot > 180.0_wp) THEN
        endlon_tot = endlon_tot - 360.0_wp
      ENDIF

!     hverint = (nverend - nversta) * dtdeh
!     lverini = (rmod( nveria*dtdeh , c1 ) > c05*dtdeh)
!  In VOF-2, flag words are written in octal instead of integer format.
      WRITE( nuverif,'(''VOF: Verification Observation File: Version 4'')' )
      WRITE( nuverif,'(''Verification period: initial date and hour '',A14)' ) &
             ydate_ref
!            yakdat1(1:10)
      WRITE( nuverif,'(21X,''start:'',F7.4,'' , end:'',F8.4)' )  tveria, tverie
!     WRITE( nuverif,'(''Verification period: initial date and time ''A10      &
!                    &,'' , length '',F5.2,'' [hrs]'')' )                      &
!            yakdat1, hverint
!     WRITE( nuverif,'(''                     data used at exact initial time''&
!                    &,1X,L1)' )      lverini
      WRITE( nuverif,'(''Set-up of the reference model run, used for thresh''  &
                     &,''old quality control (QC):'')' )
      WRITE( nuverif,'(''- LM-grid: pole:'',F7.2,F8.2,'' , lower left corner:''&
                     &,F9.3,F9.3)' )  pollat, pollon, startlat_tot, startlon_tot
      WRITE( nuverif,'(''           resolution:'',2F8.5,'' , upper right:''    &
                     &,F8.3,F9.3)' )  dlat  , dlon  , endlat_tot, endlon_tot
      WRITE( nuverif,'(''           domain size:'',2I5,I4)' )                  &
                                      ie_tot, je_tot, ke_tot
!     WRITE( nuverif,'(''           resolution:'',2F9.6,'' , domain size:''    &
!                    &,2I5,I4)' )     dlat, dlon, ie_tot, je_tot, ke_tot
      WRITE( nuverif,'(''- Initial date and time: '',A14)' )         ydate_ref
      WRITE( nuverif,'(''- QC time step : '',F6.0,'' [s]'')' )       dtqc
      WRITE( nuverif,'(''- QC thresholds: ''                                   &
                     &,    ''upper-air, vertical table:'',4F8.2)' )  qcvf
      WRITE( nuverif,'(17X,''upper-air, constant part :'',4F8.2)' )  qcc
      WRITE( nuverif,'(17X,''upper-air, time factor   :'',4F8.2)' )            &
                                                         (qctfpr(kl), kl=1,4)
      WRITE( nuverif,'(17X,''surface  , constant part :'',4F8.2)' )  qccsu
      WRITE( nuverif,'(17X,''surface  , time factor   :'',4F8.2)' )            &
                                                         (qctfpr(kl), kl=5,8)
      WRITE( nuverif,'('' '')' )
      WRITE( nuverif,'(''Number of model runs to compare with observations:''  &
                     &,I2)' )   nfcast
      IF (lvofoi) THEN
        WRITE( nuverif,'(''Domain used for verification (rotated pole: see ''  &
                       &,''above):'')' )
        WRITE( nuverif,'( 5X,''lower left corner:'',2F8.3,'' , upper right:''  &
                       &,2F8.3)' )  startlat_tot, startlon_tot                 &
                                  ,   endlat_tot,   endlon_tot
        WRITE( nuverif,'( 5X,''domain size:'',2I5)' )     ie_tot, je_tot
        WRITE( nuverif,'(''Types assigned to "model runs" in table below:'')' )
        WRITE( nuverif,'( 5X,''= 0: the "model run" is one straightforward ''  &
                       &,''model integration'')' )
        WRITE( nuverif,'( 5X,''> 0: the "model run" comprises of x-hourly''    &
                       &,'' periods from a series of cycled'')' )
        WRITE( nuverif,'(10X,"integrations starting at x-hour intervals "      &
                       &,"(x = 6, 12 or 24), and")')
        WRITE( nuverif,'(10X,"the initial date and time relates to the latest" &
                       &," of these integrations")' )
        WRITE( nuverif,'( 5X,''= 2: the "model run" is a series of analyses'')')
        WRITE( nuverif,'(''        |             |forecast time | horiz. |''   &
                       &,'' number |'')' )
        WRITE( nuverif,'('' model  |   initial   |at the end of | mesh   |''   &
                       &,''   of   |'')' )
        WRITE( nuverif,'(''  run   |date and hour| verification | width  |''   &
                       &,''vertical| description of'')' )
        WRITE( nuverif,'(''no.|type| of the run  | period [hrs] |[1/deg.]|''   &
                       &,'' levels | the model run'')' )
        zdphir = c1 / dlat
        IF (lrefend) THEN
          !    in this case, it is assumed that (mruntyp == 2) is set correctly
          ydescr   =  'analysis'
          trun     =  c0
          ntrun    = NINT( hverend * c3600 )
!         ntrun  = INT( hverend +1000.99999_wp ) - 1000
!         trun   = hverend - ntrun 
!         ntrun  = NINT( ntrun * c3600 )

          CALL get_utc_date ( ntrun, ydate_ref, c1, itype_calendar, ygrbt      &
                            , yakdat2, nzday, zhr )
!         =================
        ELSE
          !    in this case, it is assumed that (mruntyp /= 2)
          IF (mruntyp == 1) THEN
            ydescr   =  'first guess'
          ELSE
            ydescr   =  'forecast'
          ENDIF
          ygrbt    =  ydate_ref
          trun     =  hverend
        ENDIF
        WRITE( nuverif,'('' 1'',3X,I2.2,2X,A14,2X,F9.4,4X,F7.2,4X,I3,5X,A11)') &
               mruntyp, ygrbt, trun, zdphir, ke_tot, ydescr
      ENDIF
      WRITE( nuverif,'('' '')' )
    ENDIF
    mboxnow = -9999

! print VOF body
! --------------

    DO irep = 1 , nrep
      ipos = ioffs(isortrv(irep))

      mreptim = odrbufa(ipos+nvhtim)
      mboxnow = MAX( mboxnow , mreptim )

      IF ((odrbufa(ipos+nvhobt) > 0) .AND. (mreptim > mboxold)) THEN

  ! print single-level report (incl. GPS report)
  ! --------------------------------------------

        nlev   = 0
        IF (odrbufa(ipos+nvhobt) == ngps) THEN
          nlev   = -4
        ELSEIF (     (odrbufa(ipos+nvhobt) == nairep)                          &
                .OR. (odrbufa(ipos+nvhobt) == nsatob)) THEN
          nlev   = -3
        ELSE
          indst  = ipos + mxvhed + ilstidv
          IF (      (odrbufa(indst+nvbtn ) <   0)                              &
              .AND. (odrbufa(indst+nvbtx ) <   0)                              &
              .AND. (odrbufa(indst+nvbhsw) <   0)                              &
              .AND. (odrbufa(indst+nvbfg1) <   0)                              &
              .AND. (odrbufa(indst+nvbfg6) <   0)                              &
              .AND. (odrbufa(indst+nvbrad) <   0))  nlev  =  - 1
          IF (      (odrbufa(indst+nvbcl ) <   0) .AND. (nlev == -1)           &
              .AND. (odrbufa(indst+nvbvis) <   0)                              &
              .AND. (odrbufa(indst+nvbcwg) <   0)                              &
!             .AND. (odrbufa(indst+nvbwwe) <   1)                              &
              .AND. (odrbufa(indst+nvbr12) <  -1))  nlev  =  - 2
        ENDIF
        DO icl = 1 , ilstidv
          ystid (icl:icl) = CHAR( odrbufa(ipos+mxvhed+icl) )
        ENDDO

! header of single-level report
        WRITE( nuverif,'(I4,1X,A9,1X, 5I6, I3,I4, o11,o11,I2,I3, 2I5)' )       &
               nlev, ystid, (odrbufa(ipos+kl), kl=2,mxvhed)

        ipos  =  ipos + mxvhed + ilstidv
! body of single-level report
        IF (nlev == -4) THEN
          WRITE( nuverif,'(2I5,I5,I5,I7,I6, I5,I5,o11,I4,I4)' )                &
                (odrbufa(ipos+kl), kl=1,mxvbdg)
        ELSEIF (nlev <= -2) THEN
          WRITE( nuverif,'(2I5,I5,I5,I7,I6, I5,I5,o11,I4,I4)' )                &
                (odrbufa(ipos+kl), kl=1,mxvbdu)
        ELSE
          WRITE( nuverif,'(2I5,I5,I5,I7,I6, I5,I5,o11,I4,I4, I3,I5,2o12,I5)' ) &
                (odrbufa(ipos+kl), kl=1,mxvbdt)
          itmp  =  mxvbds - mxvbdt
          IF (nlev == 0) WRITE( nuverif,'(3I5,2I4,I5)' )                       &
                               (odrbufa(ipos+mxvbdt+kl), kl=1,itmp)
        ENDIF
!       IF ((nlev == -1) .AND. (odrbufa(ipos+nvbr12) /= noval))                &
!         PRINT '(I4,1X,A8,2X, I6,I5,3I6, I2,I4, I8,I11,I2,I3, 2I5, I5)' ,     &
!                nlev, ystid, (odrbufa(ipos-mxvhed-ilstidv+kl), kl=2,mxvhed)   &
!                           ,  odrbufa(ipos+nvbr12)
        IF (nlev >= -3) ipos  =  ipos + mxvbds
        IF (nlev == -4) ipos  =  ipos + mxvbdg

! observation increments
        IF (lvofoi) THEN
          IF (nlev == -4) mxvdd = mxvddg
          IF (nlev == -3) mxvdd = mxvddu
          IF (nlev == -2) mxvdd = mxvddb
          IF (nlev == -1) mxvdd = mxvddt
          IF (nlev ==  0) mxvdd = mxvdds
          IF (nlev == -4) yuvofmt = '(I6,I6                              )'
          IF (nlev == -3) yuvofmt = '(I6,I6,I6,I5                        )'
          IF (nlev == -2) yuvofmt = '(I6,I6,I6,I5,I6                     )'
          IF (nlev == -1) yuvofmt = '(I6,I6,I6,I5,I6,2I4,I5,I5           )'
          IF (nlev ==  0) yuvofmt = '(I6,I6,I6,I5,I6,2I4,I5,I5,2I5,2I5,I5)'
          WRITE( nuverif, yuvofmt, IOSTAT=nstat )                              &
                (odrbufa(kl), kl=ipos+1,ipos+mxvdd)
          IF (nlev >= -2) ipos = ipos + mxvdds
          IF (nlev == -3) ipos = ipos + mxvddu
          IF (nlev == -4) ipos = ipos + mxvddg
        ENDIF

      ELSEIF (mreptim > mboxold) THEN
          ! this condition avoids re-printing of multi-level aircraft reports
          ! just after having re-created them again

  ! print multi-level report
  ! ------------------------

        nlev   = ABS( odrbufa(ipos+nvhpas) )
        odrbufa(ipos+nvhpas) = 1 - SIGN( 1 , odrbufa(ipos+nvhpas) )
        odrbufa(ipos+nvhobt) = - odrbufa(ipos+nvhobt)
        DO icl = 1 , ilstidv
          ystid (icl:icl) = CHAR( odrbufa(ipos+mxvhed+icl) )
        ENDDO

! header of multi-level report
        WRITE( nuverif,'(I4,1X,A9,1X, 5I6, I3,I4, o11,o11,I2,I3, 2I5)' )       &
               nlev, ystid, (odrbufa(ipos+kl), kl=2,mxvhed)

        ipos   = ipos + mxvhed + ilstidv
! body of multi-level report (length of one level = mxvbdm)
        DO klev = 1 , nlev
          WRITE( nuverif,'(2I5,I5,I5,I7,I6, I5,I5,o11,I4)' )                   &
                (odrbufa(ipos+kl), kl=1,mxvbdm)
          ipos   = ipos + mxvbdm
        ENDDO

! observation increments
        IF (lvofoi) THEN
          DO klev = 1 , nlev
            WRITE( nuverif, '(I6,I6,I6,I5,I6)', IOSTAT=nstat )                 &
                  (odrbufa(ipos+kl), kl=1,mxvddm)
            ipos   = ipos + mxvddm
          ENDDO
! surface pressure for multi-level reports
          WRITE( nuverif, '(I7)', IOSTAT=nstat )  odrbufa(ipos+1)
          ipos   = ipos + 1
        ENDIF

      ENDIF
    ENDDO

    mboxold = MAX( mboxnow , mboxold )

    DEALLOCATE ( ioffs     , STAT = nzerr )
    DEALLOCATE ( isortrv   , STAT = nzerr )

! print VOF end line
! ------------------

    IF (lastpr) THEN
      kk = 0
      WRITE( nuverif,'(''  -9 xxxxx     '', 5I6, I3,I4, I11,I11,I2,I3          &
                     &, 2I5)' )   kk,kk,kk,kk,kk,kk,kk,kk,kk,kk,kk,kk,kk
    ENDIF

! Close YUVERIF file
! ------------------

    IF ((lfirst) .OR. (lastpr) .OR. (nrep >= 1))  CLOSE (nuverif)

  ENDIF  !  (my_cart_id == ircv)

  DEALLOCATE ( odrbufa   , STAT = nzerr )

! De-allocate smlbdy etc
! ----------------------

  IF (MOD(mveripr,2) == 0) THEN
    DEALLOCATE ( ssgbdy  , STAT = nzerr )
    DEALLOCATE ( sgpbdy  , STAT = nzerr )
    DEALLOCATE ( stvbdy  , STAT = nzerr )
    DEALLOCATE ( smlbdy  , STAT = nzerr )
    DEALLOCATE ( dmlhed  , STAT = nzerr )
  ENDIF

  lfirst = .FALSE.

!-------------------------------------------------------------------------------
! End of module procedure obs_print_vof
!-------------------------------------------------------------------------------
END SUBROUTINE obs_print_vof


!===============================================================================

ELEMENTAL REAL (KIND=wp) FUNCTION rmod  ( ztti, ziv )
  !----------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  ztti, ziv
  !  uses parameter:                           epsy
  !---------------------------------------------------------------------------
  ! MOD function for positive REALS
  !---------------------------------------------------------------------------
  !
  rmod  =  ztti  -  ziv * REAL(INT( ztti/ziv + epsy ), wp) + epsy
  !
END FUNCTION rmod

!-------------------------------------------------------------------------------

END MODULE src_obs_print_vof
