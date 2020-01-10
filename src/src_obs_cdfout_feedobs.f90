!+ Source module for the observation processing in the data assimilation mode
!-------------------------------------------------------------------------------

MODULE src_obs_cdfout_feedobs

!-------------------------------------------------------------------------------
! Description:
!   This module writes observation reports and model equivalents (model values
!   projected onto the observation space) to a NetCDF file which is used either
!   as input for the LETKF (3DVAR package) or for verification purposes.
!
!   This module contains the following module procedures:
!    - obs_org_cdfout_feedobs : organises the writing to NetCDF feed-obs file
!    - obs_neff_fill_buffer   : stores reports in buffer arrays
!    - obs_neff_sort_buffer   : sorts list of reports in buffer arrays
!    - obs_write_cdf_feedobs  : masters the writing to NetCDF feed-obs file
!
!   This module also contains elemental functions, formerly statement functions:
!    - ireplace1       : paste 1 bit from one word into another (integer) word
!    - ibit1           : return 1 bit at given position from an (integer) word
!
!   It uses from:
!    - src_obs_cdfin_util   : - obs_gather_buffer
!                             - ncfeedobs_status
!    - src_obs_operator_conv: - frac2octas
!    - utilities            : - get_utc_date
!                             - sortrx
!    - parallel_utilities   : - global_values
!    - environment          : - model_abort
!    - mo_fdbk              : - setup_fdbk      - close_fdbk
!                             - create_fdbk     - cleanup_fdbk
!                             - add_verification
!                             - (t_fdbk)
!    - mo_netcdf_param
!    - mo_t_table           : - bit1_pos_nr
!    - mo_t_netcdf_file     : - (PLEN)
!    - mo_fdbk_cosmo        : - write_report
!                             - (t_account)
!
!   Data modules used: - data_parameters
!                      - data_obs_lib_cosmo
!                      - data_obs_record
!                      - mo_fdbk_tables
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
!  Initial release, based on code written by Marec Lazanowicz, re-structured
!  and updated.
! V4_23        2012/05/10 Ulrich Schaettler
!  Switched off writing of latex-files
! V4_24        2012/06/22 Hendrik Reich
!  Adapted length of strings for date variables
! V4_28        2013/07/12 Christoph Schraff
!  - GPS station height filled into height of observation body in feedobs file,
!    to avoid problems in LETKF.
!  - Also simulated GPS ZPD written to feedobs file.
!  - Processing of new flags for QC against lateral boundary fields.
!  - 'veri_data' in 't_acc_body' with fixed size 1 instead of pointer, i.e.
!    simulated obs from at most 1 model run can be stored. This avoids an
!    additional non-vectorizable inner loop for memory (de-)allocation of
!    'report'.
!  - Functions 'bit1_pos_nr' , 'ncfeedobs_status' moved to modules 'mo_t_table'
!    resp. 'src_obs_cdfin_util' to make them available to the 3dvar package.
!  - Bug fix: write also flag of lapse rate check in feedobs file.
!  - Bug fix: correct order of tmin and tmax in table 'ibs_varno'.
!  - Flight phase is set in the report header only, not any more in the body.
!  - Dataset flag set for vertical velocity.
!  - Processing of SOR arrays with simulated observations instead of increments.
!  - Statement functions replaced by elemental or intrinsic functions.
! V4_29        2013/10/04 Christoph Schraff
!  Bug fix: condition 'lvofoi' removed for de-allocation of 'ssgbdy' etc.
! V5_1         2014-11-28 Christoph Schraff, Oliver Fuhrer
!  Bug fix for 'iveri_run_type'; set 'lrefend' only if (mruntyp == 2).
!  Bug fix to restore modifications of V4.24 erroneously omitted in V4.28,
!  required to deal with model runs not starting / ending at full hours.
!  Replaced ireals by wp (working precision). But take care for NetCDF output, 
!    which is in NF_FLOAT (sp), so make all types used for NetCDF output to use
!    single precision floats. (OF)
! V5_3         2015-10-09 Christoph Schraff
!  - Writing of additional elements from surface reports for verification
!    (model equivalent calculator: MEC) purposes: wind speed + direction,
!    dewpoint temperature, radiation, 3-hourly precip and wind gusts, max. and
!    min. T-2m over arbitrary periods instead of fixed 12 hours.
!  - Obs for which the obs operator is not implemented in COSMO but exists in
!    MEC (i.e. temporally non-local obs), are set to 'st_passive' (instead of
!    'st_obs_only', even though they have missing value in 'veri-data'!) and
!    flagged as 'fl_obstype'.
!  - Feedback file body entries 'dlat', 'dlon' added (according to 3DVAR V1_42).
!  - 'veri_data' defined as pointer to allow for writing data from more than one
!    model run (as required for satellite radiances, see src_obs_rad.f90).
! V5_4         2016-03-10 Christoph Schraff
!  Only removal of comments related to AOF interface.
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
    sp,        & ! KIND-type parameter for single precision variables
    iintegers    ! KIND-type parameter for standard integer variables

!-------------------------------------------------------------------------------

USE data_obs_lib_cosmo, ONLY :  &

! 1. General parameters
! ---------------------

    c0         ,& ! standard real constant 0.0
!   c1         ,& ! standard real constant 1.0
    c05        ,& ! standard real constant 0.5
    rmdi       ,& ! =-1.E31_wp : commonly used missing data indicator
    rmdich     ,& ! =-1.E30_wp : commonly used check value for miss data
    epsy       ,& ! = 1.E-8_wp : commonly used very small value > 0

! 5. CMA observation type and code type numbers
! ---------------------------------------------

    nsynop     ,& ! SYNOP reports
    nairep     ,& ! AIREP reports (all aircraft reports)
    nsatob     ,& ! SATOB reports
!   ndribu     ,& ! DRIBU reports
    ntemp      ,& ! TEMP  reports
    npilot     ,& ! PILOT reports
!   nsatem     ,& ! SATEM reports
!   nsattv     ,& ! SATEM reports
!   ngps       ,& ! GPS   reports
    nra_eu        !   sodar/rass report (European)

! end of data_obs_lib_cosmo

!-------------------------------------------------------------------------------

USE data_obs_record, ONLY :   &

!       1.1    ODR header format
!       ------------------------

!       1.1.1  Header formats of ODR reports:'omlhed','osghed','ogphed','otvhed'
!              -----------------------------------------------------------------
!   mxrhed     ,& ! header length of multi-level reports
!   mxshed     ,& ! header length of single-level reports
!   mxghed     ,& ! header length of GPS reports
!   mxthed     ,& ! header length of satellite retrieval reports
    nhilon     ,& ! longitude of observing station
    nhjlat     ,& ! latitude  of observing station
    nhalt      ,& ! station altitude [m]
    nhtime     ,& ! time of observat. in forecast hours
    nhsurf     ,& ! height of model grid pt. to which obs. is assigned
    nhtddb     ,& ! data base decoding time in forecast hours
    nhsolz     ,& ! solar zenith angle [deg]

!       1.1.2  Header formats of ODR reports:'momlhd','mosghd','mopghd','motvhd'
!              -----------------------------------------------------------------
!   mxrhdf     ,& ! header length of multi-level reports
!   mxshdf     ,& ! header length of single-level reports
!   mxghdf     ,& ! header length of GPS reports
!   mxthdf     ,& ! header length of satellite retrieval reports
    nhitot     ,& ! global x-coord. of grid pt. assigned to obs
    nhjtot     ,& ! global y-coord. of grid pt. assigned to obs
    nhobtp     ,& ! observation type
    nhcode     ,& ! code type
    nhschr     ,& ! station characteristics                      (see 1.1.4)
    nhflag     ,& ! report flags (obs type, surf., alt., sta ID) (see 1.2.1)
    nhpass     ,& ! flag for report being set to 'passive'       (see 1.1.4)
    nhqcfw     ,& ! status of QC and of writing to feedobs files, and
                  !   QC flags for surface pressure increments from TEMPs
    nhcorr     ,& ! update sequence number (station correction indicator)
    nhcat      ,& ! data     category (from BUFR Section 1)
    nhcats     ,& ! data sub-category (from BUFR Section 1)
    nhcent     ,& ! originating centre  +  (1000* sub-centre)
    nhstid     ,& ! station identity number
    nhhrmn     ,& ! absolute exact observation time [hhmm]
    nhsyhr     ,& ! absolute nominal (synoptic) observation time [yymmddhh]
    nhstyp     ,& ! sing-lv obs: station type (buoy: MQOBL, BUFR Table 002149,
                  !                            else: NIX  , BUFR Table 002001)
    nhrtyp     ,& ! radiosonde type    (NRARA, see WMO common code Table C2)
    nhnlev     ,& ! number of obs. levels (for multi-level reports)
    nhvqcf     ,& ! for satellite retrieval: threshold quality control flags
    nhtrac     ,& ! tracking technique (NSASA, see WMO common code Table C7)
    nhrad      ,& ! solar and IR radiation correction (NSR, BUFR Table 002013)
    nhna4      ,& ! instrument type                   (NA4, BUFR Table 002003)

!       1.1.3  Header formats of ODR reports:'yomlhd','yosghd','yopghd','yotvhd'
!              -----------------------------------------------------------------
    ilstid     ,& ! character length of the station identity
    ilstidp       ! char. length used for printing the station ID
                  ! Note: (ilstid >= ilstidg >= ilstidp), cf. data_nudge_gather

USE data_obs_record, ONLY :   &

!       1.2    Bit patterns for packed information in ODR (and VOF) header
!       ------------------------------------------------------------------
!   nvpabp     ,& ! bit pos. for report set passive since it is     nvschr
                  !              used in a multi-level pseudo report
!   nvpsbp     ,& ! bit pos. for report set passive since 1 of next   "
                  !              5 flags or flight track flag applies
!   nvobbp     ,& ! bit pos. for flag: 'station location out of       "
                  !                     user-specified area'
!   nvalbp     ,& ! bit pos. for flag: 'distance model orography -    "
                  !                     station altitude too large'
!   nvbkbp     ,& ! bit pos. for flag: 'blacklisted station (ship)'   "
!   nvexbp     ,& ! bit pos. for flag: 'observation or code type      "
                  !                     excluded at station location'
!   nvrdbp     ,& ! bit pos. for flag: 'redundant report'             "
    nvsebp     ,& ! bit pos. for report located at sea grid pt.       "
    nvscbp     ,& ! bit pos. for station correction indicator         "
    nvapbp     ,& ! bit pos. for phase of flight (aircraft)           "
    nvapoc     ,& ! no. of bits occ. by (extended) phase of flight    "
    nvaabp     ,& ! bit pos. for aircraft roll angle (code)           "
    nvaaoc        ! no. of bits occ. by aircraft roll angle           "

USE data_obs_record, ONLY :   &

!       1.3    ODR body format
!              ---------------

!       1.3.0  Number of levels in multi-level ODR 'omlbdy', 'momlbd'
!              ------------------------------------------------------
!   maxrsl     ,& ! max. number of levels in multi-level ODR


!       1.3.1  Body format of ODR of multi-level reports: 'omlbdy'
!              ---------------------------------------------------
!   mxrbdy     ,& ! body length of multi-level reports
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
!   nbtzio     ,& ! longitude in grid pt. units
!   nbtzjo     ,& ! latitude  in grid pt. units
!   nbtlop     ,& ! LOG( pressure )
    nbtdrh     ,& ! bias correction for relative humidity [/]
    nbtw       ,& ! vertical velocity [m/s]
    nbtsnr     ,& ! signal to noise ratio
    nbtuac     ,& ! accuracy (std dev from data provider) of horiz. wind [m/s]


!       1.3.2  Body format of ODR of multi-level report flags: 'momlbd'
!              --------------------------------------------------------
!   mxrbdf     ,& ! body length of multi-level reports
    nbtflg     ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbterr     ,& ! pre-processing status flags  (bit p., see below: 'nb?err')
    nbtqcf     ,& ! threshold quality control flags      (see below: 'nb?qcf')
    nbtlsg        ! level id (bit pattern, as in NetCDF statistics file)
!   nbtlid        ! level identity          (bit pattern, see below: 'nb?lid')

USE data_obs_record, ONLY :   &

!       1.3.3  Body format of ODR of surface reports: 'osgbdy'
!              -----------------------------------------------
!   mxsbdy     ,& ! body length of single-level reports
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
    nbsvis     ,& ! (horizontal) visibility                            [m]
    nbsff      ,& ! wind speed                                         [m/s]
    nbsdd      ,& ! wind direction                                     [deg]
    nbstd         ! dewpoint temperature                               [K]

USE data_obs_record, ONLY :   &

    nbsrr1     ,& ! precipitation amount over 1  hour                  [mm]
    nbsrr3     ,& ! precipitation amount over 3  hours                 [mm]
    nbsrr6     ,& ! precipitation amount over 6  hours                 [mm]
    nbsr12     ,& ! precipitation amount over 12 hours                 [mm]
    nbsr24     ,& ! precipitation amount over 24 hours                 [mm]
    nbsfgv     ,& ! max. derived equivalent vertical gust (aircraft)   [m/s]
    nbsfg1     ,& ! max. wind speed of gusts over 1 hour               [m/s]
    nbsfg3     ,& ! max. wind speed of gusts over 3 hours              [m/s]
    nbsfg6     ,& ! max. wind speed of gusts over 6 hours              [m/s]
    nbstn      ,& ! minimum temperature (at 2m, in period 'nbsttr')    [K]
    nbstx      ,& ! maximum temperature (at 2m, in period 'nbsttr')    [K]
    nbsrad     ,& ! global    solar    radiation, sum over 1 hour      [J/m2]
    nbsrdd     ,& ! diffuse   solar    radiation, sum over 1 hour      [J/m2]
    nbsrdt     ,& ! long-wave downward radiation, sum over 1 hour      [J/m2]
    nbshsw     ,& ! total snow depth                                   [m]
    nbsdrh        ! bias correction for relative humidity [/]

USE data_obs_record, ONLY :   &

!       1.3.4  Body format of ODR of surface report flags: 'mosgbd'
!              ----------------------------------------------------
!   mxsbdf     ,& ! body length of single-level reports
    nbsflg     ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbserr     ,& ! pre-processing status flags  (bit p., see below: 'nb?err')
    nbsqcf     ,& ! threshold quality control flags      (see below: 'nb?qcf')
    nbslid     ,& ! SYNOP: pressure code (SYNOP)   (code, see below: 'nbslid')
                  ! else : level identity   (bit pattern, see below: 'nb?lid')
!   nbscwg     ,& ! combined cloud and weather group (set of classes, s below)
    nbswwe     ,& ! NetCDF read, SYNOP: weather and ground group word  (below)
    nbstur     ,& ! NetCDF read, Aircraft: degree of turbulence WMO Tab 011031
                  !   (not contained in merged multi-level aircraft reports !)
    nbsclg     ,& ! general           cloud       group (code)
    nbscl1     ,& ! first  individual cloud layer group (code)
    nbscl2     ,& ! second individual cloud layer group (code)
    nbscl3     ,& ! third  individual cloud layer group (code)
    nbscl4     ,& ! forth  individual cloud layer group (code)
    nbsttr        ! time periods    (bit pattern of values, see below: nbsttr)

USE data_obs_record, ONLY :   &

!       1.3.5  Body format of ODR of GPS reports: 'ogpbdy'
!              -------------------------------------------
!   mxgbdy     ,& ! body length of GPS reports
    nbgtze     ,& ! error in total zenith delay [mm]
    nbgzpd     ,& ! zenith path delay (total zenith delay)             [mm]
    nbgzwd     ,& ! zenith wet delay [mm]
!   nbgiwv     ,& ! integrated water vapour [mm]
    nbgp       ,& ! pressure [Pa]
    nbgt       ,& ! temperature [K]
    nbgrh      ,& ! relative humidity [/]
    nbgz       ,& ! height [m]
    nbgbia     ,& ! bias correction to integrated water vapour [mm]
    nbgiwa     ,& ! adjusted (bias corrected) integrated water vapour [mm]

!       1.3.6  Body format of ODR of GPS report flags: 'mogpbd'
!              -------------------------------------------------
!   mxgbdf     ,& ! body length of GPS reports
    nbgflg     ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbgerr     ,& ! pre-processing status flags  (bit p., see below: 'nb?err')
    nbgqcf     ,& ! threshold quality control flags      (see below: 'nb?qcf')
!   nbglid     ,& ! level identity          (bit pattern, see below: 'nb?lid')

!       1.3.7  Body format of ODR of sat retrieval reports: 'otvbdy'
!              -----------------------------------------------------
!   mxtbdy     ,& ! body length of multi-level reports
    nbvt       ,& ! temperature [K]
    nbvrh      ,& ! relative humidity [/]
    nbvp       ,& ! pressure [Pa]  (mandatory)

!       1.3.8  Body format of ODR of sat retrieval report flags: 'motvbd'
!              ----------------------------------------------------------
!   mxtbdf     ,& ! body length of sat retrieval reports
    nbvflg     ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbverr        ! pre-processing status flags  (bit p., see below: 'nb?err')

USE data_obs_record, ONLY :   &

!       1.4    Bit patterns for packed information in ODR (and VOF) body
!              ---------------------------------------------------------
!       1.4.2  Other bit patt. for packed info in ODR (VOF) body, general words
!              ----------------------------------------------------------------
    nvru       ,& ! bit pos. for status/QC flags for horiz. wind  nb?err/nb?qcf
    nvrt       ,& ! bit pos. for status/QC flags for temperature        "
    nvrq       ,& ! bit pos. for status/QC flags for humidity           "
    nvrz       ,& ! bit pos. for status/QC flags for pressure/height    "
    nvrw       ,& ! bit pos. for status/QC flags for vertical wind      "
    nvriwv     ,& ! bit pos. for status/QC flags for IWV                "
    nvrzpd     ,& ! bit pos. for status/QC flags for zenith path delay  "
!   nvrspd     ,& ! bit pos. for status/QC flags for slant path delay   "
    nvrct      ,& ! bit pos. for status/QC flags for (total) cloud      "
    nvrcl      ,& ! bit pos. for status/QC flags for low     cloud      "
    nvrcm      ,& ! bit pos. for status/QC flags for middle  cloud      "
    nvrch      ,& ! bit pos. for status/QC flags for high    cloud      "
    nvrcbs     ,& ! bit pos. for status/QC flags for cloud base height  "
    nvrzbc     ,& ! bit pos. for temporary flag: QC ag. LBC pressure    "
    nvrqbc     ,& ! bit pos. for temporary flag: QC against LBC IWV     "
    nvfubp     ,& ! bit pos. for main flag on wind                    nb?flg
    nvftbp     ,& ! bit pos. for main flag on temperature               "
    nvfqbp     ,& ! bit pos. for main flag on humidity                  "
    nvfzbp     ,& ! bit pos. for main flag on pressure / geopot.        "
    nvfgbp     ,& ! bit pos. for main flag on integr. water vapour      "
    nvfaoc     ,& ! no. of bits occ. by each main flag                  "
    nvfbps     ,& ! bit pattern for main flags:                         "
    nvfboc     ,& ! no. of bits occ. for main flags                     "
    nvflbp     ,& ! bit pos. for level flag: level below surface        "
!   nvfloc     ,& ! no. of bits occ. by level flag                      "
    ntxbp      ,& ! bit pos. for time period for T-max (1,..,24 hr)   "
    ntnbp      ,& ! bit pos. for time period for T-min (1,..,24 hr)   "
    ntxoc      ,& ! no. of bits occ. by each time period              "
    nvlidp     ,& ! level id. bit pattern                             nb?lid
    nvlido        ! no. bits occ. by each indicator in level id.        "

USE data_obs_record, ONLY :   &

!       1.4.3  Bit patterns for 'optional groups' in ODR body 'mosgbd' (and VOF)
!              -----------------------------------------------------------------
!   ------------- ! bit positions           ----------------------- nbswwe
                  ! --> weather and ground group word
    nvw0bp     ,& ! bit position for ww (present wea.) code (WMO Table 020003)
!   nvw1bp     ,& ! bit position for w  (past weather) code (WMO Table 020004)
!   nvwtbp     ,& ! bit position for time period of w                  [h]
!   nvcqbp     ,& ! bit position for refined quality flag on ccl
!   nveebp     ,& ! bit position for state of ground        (WMO Table 020062)
!   nrtrbp     ,& ! bit position for code of precipitation
                  !     measurement duration       [Code table 4019, keys 0-7]
    nvw0oc     ,& ! no. bits occupied for ww code           (WMO Table 020003)
!   nvw1oc     ,& ! no. bits occupied for w  code           (WMO Table 020004)
!   nvwtoc     ,& ! no. bits occupied for time period of w             [h]
!   nvcqoc     ,& ! no. bits occupied by refined quality flag on ccl
!   nveeoc     ,& ! no. bits occupied for state of ground   (WMO Table 020062)
!   nrtroc     ,& ! no. bits occupied by precip obs. duration code
                  !
!   ------------- ! bit positions           ----------------------- nbsclg
                  ! --> general    cloud       group word
!   nxclbp     ,& ! bit position for low/middle cloud amount(WMO Table 020011)
!   nctlbp     ,& ! bit position for low    cloud type code (WMO Table 020012)
!   nctmbp     ,& ! bit position for middle cloud type code (WMO Table 020012)
!   ncthbp     ,& ! bit position for high   cloud type code (WMO Table 020012)
    nxsgbp     ,& ! bit position for vertic. signific. code (WMO Table 008002)
!!  nxsgoc     ,& ! no. bits occupied for vert. signf. code (WMO Table 008002)
!!  nxcloc     ,& ! no. bits occupied for low/mid cld amount(WMO Table 020011)
!!  nxctoc     ,& ! no. bits occupied for   cloud type code (WMO Table 020012)
                  !
!   ------------- ! bit positions           ----------------------- nbscl?
                  ! --> individual cloud layer group words
    nxclbp     ,& ! bit position for cloud amount      code (WMO Table 020011)
    nxctbp     ,& ! bit position for cloud type        code (WMO Table 020012)
    nxbsbp     ,& ! bit position for cloud base height                 [m]
    nxsibp     ,& ! bit position for vertic. signific. code (WMO Table 008002)
    nxcloc     ,& ! no. bits occupied for cloud amount code (WMO Table 020011)
    nxctoc     ,& ! no. bits occupied for cloud type   code (WMO Table 020012)
    nxbsoc     ,& ! no. bits occupied for cloud base height            [m]
    nxsgoc     ,& ! no. bits occupied for vert. signf. code (WMO Table 008002)

!       1.5    Further quantities related to ODR
!              ---------------------------------
    imdi       ,& ! missing data indicator for ODR integers (2^31-1)
    ntotml     ,& ! tot. number of stored multi-level reports
    ntotsg     ,& ! tot. number of stored single-level reports
    ntotgp     ,& ! tot. number of stored GPS reports
    ntottv        ! tot. number of stored satellite retrievals

USE data_obs_record, ONLY :   &

!       2.     Observation data records (ODR)
!       -------------------------------------
    omlbdy     ,& ! body of multi-level ODR
    momlbd     ,& ! body of multi-level ODR
    omlhed     ,& ! header of multi-level ODR
    momlhd     ,& ! header of multi-level ODR
    osgbdy     ,& ! body of single-level ODR
    mosgbd     ,& ! body of single-level ODR
    osghed     ,& ! header of single-level ODR
    mosghd     ,& ! header of single-level ODR
    ogpbdy     ,& ! body of GPS ODR
    mogpbd     ,& ! body of GPS ODR
    ogphed     ,& ! header of GPS ODR
    mogphd     ,& ! header of GPS ODR
    otvbdy     ,& ! body of satellite retrieval ODR
    motvbd     ,& ! body of satellite retrieval ODR
    otvhed     ,& ! header of satellite retrieval ODR
    motvhd     ,& ! header of satellite retrieval ODR
    yomlhd     ,& ! header of multi-level ODR
    yosghd     ,& ! header of single-level ODR
    yogphd     ,& ! header of GPS ODR
    yotvhd        ! header of satellite retrieval ODR

USE data_obs_record, ONLY :   &

!       3.1  Format of SOR (for model values or observation increments)
!            -------------
    nso_u      ,& ! u wind component                               [m/s]
    nso_v      ,& ! v wind component                               [m/s]
    nso_t      ,& ! temperature                                    [K]
    nso_rh     ,& ! relative humidity                              [%] 
    nso_p      ,& ! pressure (sfc.) | geopotential (upper-air)  [Pa] | [m2/s2]
    nso_ct     ,& ! total cloud cover                              [octas]
    nso_cl     ,& ! low cloud cover                                [octas]
    nso_cm     ,& ! mid-level cloud cover                          [octas]
    nso_ch     ,& ! high cloud cover                               [octas]
    nso_cbs    ,& ! cloud base height (above surface)              [m]
    nso_ff     ,& ! wind speed                                     [m/s]
    nso_dd     ,& ! wind direction                                 [deg]
    nso_td     ,& ! dewpoint temperature                           [K]
    nso_iq     ,& ! integrated water vapour (increment)            [mm]
    nso_zpd    ,& ! zenith total path delay                        [mm]
!   nso_ps     ,& ! pressure                                       [Pa]

!       3.2  SOR arrays (for model values or observation increments)
!            ----------
    smlbdy     ,& ! body of multi-level SOR
    ssgbdy     ,& ! body of single-level SOR
    sgpbdy     ,& ! body of GPS (IWV) SOR
    stvbdy     ,& ! body of satellite retrieval SOR
    dmlhed     ,& ! single-level part of multi-level SOR

!       4.     Masking constants
!       ------------------------
    nibits        ! masking constants

! end of data_obs_record

!-------------------------------------------------------------------------------

  USE mo_fdbk_tables,          ONLY :  &

    VT_FIRSTGUESS     ,& ! run type: first guess
    VT_ANALYSIS       ,& ! run type: analysis
    VT_FORECAST       ,& ! run type: forecast
    RC_ASS            ,& ! run class: assimilation cycle
    ST_ACTIVE         ,& ! status: used in the assimilation
    ST_MERGED         ,& ! status: not used, merged into multi-level report
    ST_PASSIVE        ,& ! status: not used, only monitored
    ST_REJECTED       ,& ! status: not used due to suspicious quality
    ST_PAS_REJ        ,& ! status: passive and rejected
    ST_OBS_ONLY       ,& ! status: obs. only, no model equivalent available
    FL_OBSTYPE        ,& ! passive report type (at obs. location)
    FL_BLACKLIST      ,& ! blacklist (or not whitelist)  (in particular:
                         !   bad (hard blacklisted) aircraft station identity)
    FL_HEIGHT         ,& ! location not in valid height range
    FL_PRACTICE       ,& ! bad reporting practice / insufficient data
    FL_DATASET        ,& ! dataset qality flags
    FL_MERGE          ,& ! (report used for) merged report (only, e.g. TEMP ABCD
                         !   or single-level aircraft used for multi-level rep.)
!   FL_THIN           ,& ! thinning
    FL_GROSS          ,& ! gross error flag
    FL_FG             ,& ! observation minus first guess check
    FL_FG_LBC         ,& ! observation minus lateral boundary field check
    FL_NONE           ,& ! no flag set
    LS_SURFACE        ,& ! surface
    LS_SIGN           ,& ! significant level
    flags                ! flag table (of type t_table)

  USE mo_fdbk_tables,          ONLY :  &
                         !   variable numbers:
    VN_U              ,& ! u-component of wind
    VN_V              ,& ! v-component of wind
    VN_PWC            ,& ! precipitable water content kg/m**2
    VN_RH             ,& ! relative humidity
    VN_RH2M           ,& ! 2 metre relative humidity
    VN_T              ,& ! upper air temperature
    VN_T2M            ,& ! 2 metre temperature
    VN_TD2M           ,& ! 2 metre dew point
    VN_PTEND          ,& ! pressure tendency
    VN_WW             ,& ! present weather
    VN_VV             ,& ! visibility
    VN_NH             ,& ! cloud base height
    VN_N_L            ,& ! low cloud amount
    VN_SDEPTH         ,& ! snow depth
    VN_TRTR           ,& ! time period of information  (h)
    VN_RR             ,& ! precipitation amount        (kg/m^2)
    VN_JJ             ,& ! maximum temperature         (K)
    VN_GCLG           ,& ! general cloud group         Table
    VN_N              ,& ! total cloud amount          (20011)
    VN_PS             ,& ! surface pressure            (Pa)
    VN_DD             ,& ! wind direction              (degree)
    VN_FF             ,& ! wind force                  (m/s)
!   VN_RAWBT          ,& ! brightness temperature      (K)
    VN_U10M           ,& ! 10m u-component of wind     (m/s)
    VN_V10M           ,& ! 10m v-component of wind     (m/s)
    VN_W              ,& ! vertical speed              (m/s)
!   VN_VT             ,& ! virtual temperature         (K)
    VN_HEIGHT         ,& ! height                      (m)
    VN_FLEV              ! nominal flight level        (m)

  USE mo_fdbk_tables,          ONLY :  &
                         !   variable numbers:
!   VN_RREFL          ,& ! radar reflectivity          (Db)
!   VN_PDELAY         ,& ! atmospheric path delay      (m)
    VN_ICLG           ,& ! individual cloud layer group Table
    VN_N_M            ,& ! middle cloud amount          WMO table 020011
    VN_N_H            ,& ! high   cloud amount          WMO table 020011
    VN_RAD_GL         ,& ! 1-h global solar radiation  (J/m2) 
    VN_RAD_DF         ,& ! 1-h diffuse solar radiat.   (J/m2) 
    VN_RAD_LW         ,& ! 1-h long-wave radiation     (J/m2) 
    VN_ZPD            ,& ! zenith path delay
    VN_ZWD            ,& ! zenith wet delay
!   VN_SPD            ,& ! slant path delay
    VN_GUST           ,& ! wind gust                   (m/s)
    VN_P              ,& ! pressure                    (Pa)
    VN_TMIN           ,& ! minimum Temperature         (K)
    VN_NUM               ! ordinal (channel) number    (  )

! end of mo_fdbk_tables

!-------------------------------------------------------------------------------

  USE mo_netcdf_param

!-------------------------------------------------------------------------------

  USE mo_t_table,              ONLY :  &
    bit1_pos_nr          ! first bit = 1 in a word searched in a given order

!-------------------------------------------------------------------------------

  USE mo_t_netcdf_file,        ONLY :  &
!   t_netcdf_file        !
    PLEN                 ! max. length of NetCDF file path name

!-------------------------------------------------------------------------------

  USE mo_fdbk ,                ONLY :  &
    t_fdbk,            & !
    setup_fdbk,        & !
    create_fdbk,       & !
    close_fdbk,        & !
    cleanup_fdbk,      & !
    add_verification     !

!-------------------------------------------------------------------------------

  USE mo_fdbk_cosmo,           ONLY :  &
    write_report      ,& !
    t_account            !
!   t_acc_header      ,&
!   t_acc_body

!-------------------------------------------------------------------------------

 USE environment,              ONLY :  &
    model_abort          ! aborts the program in case of errors

!-------------------------------------------------------------------------------

 USE parallel_utilities,       ONLY :  &
    global_values        ! computes global values by operating on local arrays

!-------------------------------------------------------------------------------

 USE utilities,                ONLY :  &
    get_utc_date      ,& ! actual date of the forecast in different forms
    sortrx               ! permutation vector for (quick-)sorting a given array

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_util,       ONLY :  &
    obs_gather_buffer ,& ! gathers elements from buffer arrays from all nodes
    ncfeedobs_status     ! determines report status in feedback file format

!-------------------------------------------------------------------------------

USE src_obs_operator_conv, ONLY :   &

    frac2octas           ! conversion of (cloud) fraction to octas

! end of src_obs_operator_conv

!-------------------------------------------------------------------------------
! Local declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE

  ! header format in buffer array for NetCDF feedobs file (for packed info,
  ! -----------------------------------------------------  see ODR)

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    mxchiu = 26  ,& ! int. header length of upper-air (incl multi-level) reports
    mxchis = 22  ,& ! int. header length of surface rep (SYNOP, SHIP, BUOY, GPS)
    mxchra =  3  ,& ! real header length
                    !
    ! INTEGER elements
    nchlen =  1  ,& ! length of current report in integer buffer array
    nchler =  2  ,& ! length of current report in real    buffer array
    nchoff =  3  ,& ! offset of current report in real    buffer array
    nchfmt =  4  ,& ! report format in buffer array: 1: upper-air;
                    !                                2: surface (non-GPS);
                    !                                3: ground-based GPS
    nchtim =  5  ,& ! exact obs time of relative to reference time         [min]
    nchalt =  6  ,& ! station altitude                                       [m]
    nchsfc =  7  ,& ! height of model orography at station                   [m]
    nchobt =  8  ,& ! observation type
    nchcdt =  9  ,& ! code type
    nchsch = 10  ,& ! station characteristics                (see ODR: 'nhschr')
    nchflg = 11  ,& ! report flags (obs type, surface, altitude, station ID)
                    !   (only if reading from NetCDF)        (see ODR: 'nhflag')
    nchpas = 12  ,& ! flag for report being set to 'passive' (see ODR: 'nhpass')
    nchi_t = 13  ,& ! global x-coord.  \  of the grid point of the reference run
    nchj_t = 14  ,& ! global y-coord.  /  assigned to the obs
    nchcor = 15  ,& ! update sequence number (station correction indicator)
    nchcat = 16  ,& ! data     category (from BUFR Section 1)
    nchcas = 17  ,& ! data sub-category (from BUFR Section 1)
    nchcen = 18  ,& ! originating centre  +  (1000* sub-centre)
    nchide = 19  ,& ! station identity number
    nchsyn = 20  ,& ! nominal (synoptic) obs time rel. to reference time   [min]
    nchtdb = 21  ,& ! data base decoding time     rel. to reference time   [min]
    nchtyp = 22  ,& ! station type (radiosonde: NRARA, WMO common code Table C2
                    !              ,buoy      : MQOBL, WMO BUFR Table 002149
                    !              ,else      : NIX  , WMO BUFR Table 002001)
    ! INTEGER elements for upper-air reports only
    nchnlv = 23  ,& ! number of obs. levels (for multi-level reports)
!   nchtvq = 24  ,& ! for satellite retrieval: threshold quality control flags
    nchtrc = 24  ,& ! tracking technique (NSASA, see WMO common code Table C7)
    nchrad = 25  ,& ! solar and IR radiation correction (NSR, BUFR Table 002013)
    nchna4 = 26  ,& ! instrument type                   (NA4, BUFR Table 002003)
    ! REAL elements
    nchlon =  1  ,& ! longitude of observing station                       [deg]
    nchlat =  2  ,& ! latitude  of observing station                       [deg]
    nchsol =  3  ,& ! solar zenith angle                                   [deg]
                    !
    ilstidc = 9     ! character length of the station identity on the VOF

  ! body format in buffer array for NetCDF feedobs file (for packed info,
  ! ---------------------------------------------------  see ODR)

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    mxcbiu =  6  ,& ! int. body length of upper-air (incl. multi-level) reports
    mxcbig =  4  ,& ! int. body length of ground-based GPS IWV reports
    mxcbis = 12  ,& ! int. body length of surface-lv reports (SYNOP, SHIP, BUOY)
    mxcbru = 14  ,& ! real body length of upper-air (incl. multi-level) reports
    mxcbrg = 16  ,& ! real body length of ground-based GPS IWV reports
    mxcbrs = 35  ,& ! real body length of surface-lv reports (SYNOP, SHIP, BUOY)
    ! integer elements common to all report types
    ncberf =  1  ,& ! bit pattern: bit=1 --> status active   (see ODR: 'nb?err')
                    !    (currently only bits 0-3,5 (for v,T,q,p|z,IWV) are set)
    ncbqcf =  2  ,& ! bit pattern: bit=1 --> rejected by first guess QC,
                    !     bit pattern is as for ncberf, but currently
                    !     only bits 0,1,2,3,5 (for v,T,q,p|z,IWV) may be set.
    ncbmfw =  3  ,& ! main flag word for wind, temp., humidity, pressure/geopot
                    !                                        (see ODR: 'nb?flg')
    ncblsg =  4  ,& ! level significance (bit pattern), or pressure code (SYNOP)
                    !                                        (see ODR: 'nb?lsg')
    ! integer elements for upper-air observations
    ncbsnr =  5  ,& ! signal to noise ratio
    ncbtur =  6  ,& ! aircraft only: degree of turbulence WMO Tab 011031
                    !   (not contained in merged multi-level aircraft reports !)
    ! integer elements for surface-level observations (SYNOP, SHIP, BUOY):
      ! combined values (code tables)  --> no model values from model run
    ncbwwe =  5  ,& ! NetCDF read, SYNOP: weather and ground group word  (below)
      ! for ncbcl?, see code for 'nbsclg' and 'nbscl?' in data_obs_record
      !             except that vertical significance is excluded)
    ncbclg =  6  ,& ! general           cloud       group (code)
    ncbcl1 =  7  ,& ! first  individual cloud layer group (code)
    ncbcl2 =  8  ,& ! second individual cloud layer group (code)
    ncbcl3 =  9  ,& ! third  individual cloud layer group (code)
    ncbcl4 = 10  ,& ! forth  individual cloud layer group (code)
    ncbcsg = 11  ,& ! (vertical) level significance for all cloud groups
    ncbttr = 12     ! time periods    (bit pattern of values, see below: nbsttr)

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    ! real elements common to all report types
    ncbu   =  1  ,& ! u wind component                                     [m/s]
    ncbv   =  2  ,& ! v wind component                                     [m/s]
    ncbt   =  3  ,& ! temperature                                            [K]
    ncbrh  =  4  ,& ! relative humidity                                      [/]
    ncbp   =  5  ,& ! pressure                                              [Pa]
    ncbz   =  6  ,& ! height                                                 [m]
    ncbuer =  7  ,& ! error of observed wind component                     [m/s]
    ncbter =  8  ,& ! error of observed temperature                          [K]
    ncbqer =  9  ,& ! error of observed rel. humidity                        [/]
    ncbzer = 10  ,& ! error of obs height (upair) | pressure (surf)   [m] | [Pa]
    ncbdrh = 11  ,& ! bias correction for relative humidity                  [/]
    ! real elements for upper-air observations
    ncbw   = 12  ,& ! vertical velocity                                    [m/s]
    ncbfgv = 13  ,& ! max. derived equivalent vertical gust (aircraft)     [m/s]
    ncbuac = 14  ,& ! accuracy (std dev from data provider) of horiz. wind [m/s]
    ! real elements for ground-based GPS observations
    ncbiwa = 12  ,& ! for GPS: adjusted (bias corrected) IWV                [mm]
    ncbiwc = 13  ,& ! for GPS: bias correction for IWV                      [mm]
    ncbzwd = 14  ,& ! for GPS: zenith wet delay                             [mm]
    ncbzpd = 15  ,& ! for GPS: zenith path delay                            [mm]
    ncbzpe = 16  ,& ! for GPS: error in zenith path delay                   [mm]
    ! real elements for surface-level observations (SYNOP, SHIP, BUOY):
      ! instantaneous values
    ncbcbs = 12  ,& ! (lowest) cloud base height                             [m]
    ncbcl  = 13  ,& ! low       cloud cover        (BUFR Table 020011)   [octas]
    ncbcm  = 14  ,& ! mid-level cloud cover        (BUFR Table 020011)   [octas]
    ncbch  = 15  ,& ! high      cloud cover        (BUFR Table 020011)   [octas]
    ncbct  = 16  ,& ! total     cloud cover        (BUFR Table 020011)   [octas]
    ncbvis = 17  ,& ! (horizontal) visibility                                [m]
    ncbff  = 18  ,& ! wind speed                                           [m/s]
    ncbdd  = 19  ,& ! wind direction                                       [deg]
    ncbtd  = 20  ,& ! dewpoint temperature                                   [K]
    ncbhsw = 21     ! total snow depth                                       [m]

  INTEGER (KIND=iintegers) , PARAMETER  :: &
      ! values related to a period     --> no model values from model run
    ncbpst = 22  ,& ! (3-hourly) pressure tendency                       [Pa/3h]
    ncbrr1 = 23  ,& ! precipitation amount over 1  hour                     [mm]
    ncbrr3 = 24  ,& ! precipitation amount over 3  hours                    [mm]
    ncbrr6 = 25  ,& ! precipitation amount over 6  hours                    [mm]
    ncbr12 = 26  ,& ! precipitation amount over 12 hours                    [mm]
    ncbr24 = 27  ,& ! precipitation amount over 24 hours                    [mm]
    ncbfg1 = 28  ,& ! max. wind speed of gusts over 1 hour                 [m/s]
    ncbfg3 = 29  ,& ! max. wind speed of gusts over 3 hours                [m/s]
    ncbfg6 = 30  ,& ! max. wind speed of gusts over 6 hours                [m/s]
    ncbtn  = 31  ,& ! minimum temperature (at 2m during past 12 hrs)         [K]
    ncbtx  = 32  ,& ! maximum temperature (at 2m during past 12 hrs)         [K]
    ncbrad = 33  ,& ! global  solar radiation, sum over 1 hour            [J/m2]
    ncbrdd = 34  ,& ! diffuse solar radiation, sum over 1 hour            [J/m2]
    ncbrdt = 35     ! long-wave     radiation, sum over 1 hour            [J/m2]

  ! format of verification part (veri-part) in buffer array for NetCDF feedobs f.
  ! -----------------------------------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    mxcmru =  6  ,& ! veri_data length of upper-air (incl. multi-level) reports
    mxcmrg =  7  ,& ! veri_data length of ground-based GPS IWV reports
    mxcmrs = 15  ,& ! veri_data length of surface-lv reports (SYNOP, SHIP, BUOY)
    ncmu   =  1  ,& ! u wind component                                     [m/s]
    ncmv   =  2  ,& ! v wind component                                     [m/s]
    ncmt   =  3  ,& ! temperature                                            [K]
    ncmrh  =  4  ,& ! relative humidity                                      [/]
    ncmp   =  5  ,& ! pressure (sfc.) | geopot. (upper-air)           [Pa] | [m]
  ! elements for upper-air observations
    ncmw   =  6  ,& ! vertical velocity                                    [m/s]
  ! elements for ground-based GPS observations
    ncmiq  =  6  ,& ! for GPS: integrated water vapour (increment)          [mm]
    ncmzpd =  7  ,& ! for GPS: integrated water vapour (increment)          [mm]
  ! elements for surface-level observations (SYNOP, SHIP, BUOY)
  !    (temporally non-local elements are not computed inside COSMO)
    ncmcbs =  6  ,& ! (lowest) cloud base height                             [m]
    ncmcl  =  7  ,& ! low       cloud cover        (BUFR Table 020011)   [octas]
    ncmcm  =  8  ,& ! mid-level cloud cover        (BUFR Table 020011)   [octas]
    ncmch  =  9  ,& ! high      cloud cover        (BUFR Table 020011)   [octas]
    ncmct  = 10  ,& ! total     cloud cover        (BUFR Table 020011)   [octas]
    ncmvis = 11  ,& ! (horizontal) visibility                                [m]
    ncmff  = 12  ,& ! wind speed                                           [m/s]
    ncmdd  = 13  ,& ! wind direction                                       [deg]
    ncmtd  = 14  ,& ! dewpoint temperature                                   [K]
    ncmhsw = 15     ! total snow depth                                       [m]

  ! NetCDF feedobs file: bit patterns for cloud group words in body format
  ! ----------------------------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    ncsgbp(0:4) = (/ 0, 6, 12, 18, 24 /) ,&
                    ! bit position for vertic. signific. code (WMO Table 008002)
                    !   index   0  : for general cloud group (ncbclg)
                    !   indices 1-4: for individual cloud groups 1-4 (ncbcl?)
    ncsgoc      = 6 ! no. bits occupied for vert. signf. code (WMO Table 008002)

  ! relation between body and veri-data parts in buffer array
  ! ---------------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
                    ! for each element in the real part of the body, 'ixobs_brx'
                    ! points to the index of the corresponding observed quantity
                    ! in the veri-data part;
                    ! ixobs_brx > 0 : observed quantity exists in veri-data,
                    !                 value = corresponding index in veri-data
                    !           = 0 : observed quantity is to be written to the
                    !                 feedobs file; obs operator does not exist
                    !                 in COSMO but may exist in MEC 
                    !                 (e.g. temporally non-local obs)
                    !                 --> in COSMO: missing value in veri-data,
                    !                 --> set status 'ST_PASSIVE', 'FL_OBSTYPE'.
                    !           =-1 : observed quantity is to be written to the
                    !                 feedobs file; obs operator does not (yet)
                    !                 exist neither in COSMO nor in MEC 
                    !                 --> always missing value in veri-data,
                    !                 --> set status 'ST_OBS_ONLY'.
                    !           =-7 : observed quantity is not to be written to
                    !                 the NetCDF feedobs file
                    !           =-9 : quantity is not an observation (and hence
                    !                 not written to the feedobs file as obs.)
                    !   (note: upper-air height 'obs' at height or flight level
                    !          (i.e. height level converted into pressure level)
                    !           is not written to feedobs file as observation)
                         ! ------> upper-air obs, real elements
                         !    u     v     t     rh    p   z     obs errors b-c
    ixobs_bru (mxcbru) = (/ ncmu, ncmv, ncmt, ncmrh, -9, ncmp, -9,-9,-9,-9,-9 ,&
                            ncmw, -7, -9 /)                                   ,&
                         !   w   v-gust
                         ! ------> surface-level obs, real elements
    ixobs_brs (mxcbrs) = (/ ncmu, ncmv, ncmt, ncmrh, ncmp, -9, -9,-9,-9,-9,-9 ,&
                            ncmcbs, ncmcl, ncmcm, ncmch, ncmct, -1            ,&
                            ncmff, ncmdd, ncmtd, -1   , -1 , 0 , 0            ,&
                            0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0  /)     ,&
                         ! r06 r12 r24 fg1 fg3 fg6 tn  tx  rad rdd rdt 
                         ! ------> GPS obs, real elements
!   ixobs_brg (mxcbrg) = (/ ncmu, ncmv, ncmt, ncmrh, ncmp, -9, -9,-9,-9,-9,-9 ,
    ixobs_brg (mxcbrg) = (/  -1 ,  -1 ,  -1 ,  -1  ,  -1 , -9, -9,-9,-9,-9,-9 ,&
                            ncmiq,-9 , -1 , ncmzpd,-9 /)                      ,&
                         !   iwv  b-c  zwd    zpd  err
                         ! ------> surface-level obs, integer elements
                         !  ---flags--- ww clg cloud layers sign period
    ixobs_bis (mxcbis) = (/ -9,-9,-9,-9, -1, -1 , -1, -1, -1, -1 , -9 , -9 /)

!                        !   w   v-gust
!                        ! ------> surface-level obs, real elements
!   ixobs_brs (mxcbrs) = (/ ncmu, ncmv, ncmt, ncmrh, ncmp, -9, -9,-9,-9,-9,-9 ,&
!                           ncmcbs, ncmcl, ncmcm, ncmch, ncmct, ncmvis, ncmhsw,&
!                           0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 /)              ,&
!                        ! pst rr1 rr6 r12 r24 fg1 fg6 tn  tx

! INTEGER (KIND=iintegers) , PARAMETER  :: &
    ! masks which indicate which elements are observations:
    !   = 0 : auxilliary information (e.g. obs error, bias corr, or level info)
    !   = 1 : observation, instantaneous, model equivalent exists (in veri_data)
    !   = 2 : observation, instantaneous, model equivalent does not exist
    !   = 3 : observation, non-local in time (sum or max), model equiv. exists
    !   = 4 : observation, non-local in time, model equiv. does not exist
      !                     u v t q p z  errors b-c w v-gust
!   ixobs_bru (mxcbru) = (/ 1,1,1,1,0,1, 0,0,0,0,0, 1 , 2 , 0 /)              ,&
!   ixobs_brs (mxcbrs) = (/ 1,1,1,1,1,0, 0,0,0,0,0                            ,&
!                           1 , 1, 1, 1, 1, 1, 2, 4, 4,4,4,4, 4, 4, 4, 4 /)   ,&
      !                    cbh cl cm ch ct vis e pst  4* rr  gu gu tn tx
      !                     u v t q p z  errors b-c iwv b-c zwd zpd err
!   ixobs_brg (mxcbrg) = (/ 1,1,1,1,1,0, 0,0,0,0,0,  1 , 0 , 2 , 2 , 0 /)

  ! relation between body in buffer array and observation variable in FOF
  ! ---------------------------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    VN_VGUST = 253 ,&
    vmiss    =  -1 ,&
    ibu_varno (mxcbru) = (/ VN_U, VN_V, VN_T, VN_RH, VN_P, VN_HEIGHT, vmiss   ,&
                            vmiss, vmiss, vmiss, vmiss, VN_W, VN_VGUST        ,&
                            vmiss /)                                          ,&
    ibs_varno (mxcbrs) = (/ VN_U10M, VN_V10M, VN_T2M, VN_RH2M, VN_PS          ,&
                            VN_HEIGHT, vmiss, vmiss, vmiss, vmiss, vmiss      ,&
                            VN_NH, VN_N_L, VN_N_M, VN_N_H, VN_N, VN_VV        ,&
                            VN_FF, VN_DD, VN_TD2M, VN_SDEPTH                  ,&
                            VN_PTEND, VN_RR, VN_RR, VN_RR, VN_RR, VN_RR       ,&
                            VN_GUST, VN_GUST, VN_GUST, VN_TMIN, VN_JJ         ,&
                            VN_RAD_GL, VN_RAD_DF, VN_RAD_LW /)                ,&
    ibg_varno (mxcbrg) = (/ VN_U10M, VN_V10M, VN_T2M, VN_RH2M, VN_PS          ,&
                            VN_HEIGHT, vmiss, vmiss, vmiss, vmiss, vmiss      ,&
                            VN_PWC, vmiss, VN_ZWD, VN_ZPD, vmiss /)           ,&
    iis_varno (mxcbis) = (/ vmiss, vmiss, vmiss, vmiss, VN_WW, VN_GCLG        ,&
                            VN_ICLG, VN_ICLG, VN_ICLG, VN_ICLG, vmiss, vmiss /)

!===============================================================================

!-------------------------------------------------------------------------------
! Public and Private Subroutines
!-------------------------------------------------------------------------------

CONTAINS


!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_cdfout_feedobs" for organis. feedobs file write
!-------------------------------------------------------------------------------

SUBROUTINE obs_org_cdfout_feedobs ( lastpr, lrefend, ydate_ini                 &
                                  , hversta, hverend                           &
                                  , mruntyp, max_rep, max_body                 &
                                  , num_compute, my_cart_id, icomm_cart        &
                                  , imp_integers, imp_reals                    &
                                  , ie_tot, je_tot, ke_tot, pollat, pollon     &
                                  , dlat, dlon, startlat_tot, startlon_tot     &
                                  , iveri_ens_member, nvers                    &
                                  , yglatt_institution, yglatt_source          &
                                  , yfdbkdir                                   &
                                  , nexce_rep, nexce_bdy )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "src_obs_cdfout_feedobs" organizes the
!   writing of observational reports including the model equivalents onto
!   a NetCDF feedobs file (FOF) in a format according to the 'feedback file
!   definition'.
!
! Method:
!   If the program is run in parallel mode, the data to be printed are
!   gathered from all PEs at the root node. Mainly for this purpose, the
!   required information from the observational reports in the ODR is copied
!   in 1-dimensional arrays. After gathering these arrays at the root node,
!   the data are written onto the feedobs file from these arrays with the
!   help of routines also used in the 3DVar package.
!   03.08.10
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
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

  LOGICAL                  , INTENT (IN)  :: &
    lastpr      ,& ! last time, sth is written to ncfeedobs in this model run
    lrefend        ! .t.: verification reference time set to end of model run 
                   !      (for NUDGE runs only; for FORECAST (LETKF) or
                   !       NUDGECAST, veri_ref_time = start of model run)

  INTEGER (KIND=iintegers) , INTENT (IN)  :: &
    mruntyp     ,& ! =-1: print only observations ; =2: also print increments
                   ! (only if >0: arrays 'smlbdy' etc. have to be allocated and
                   !  filled previously and are used and de-allocated here)
    max_rep     ,& ! max. number of reports      in NetCDF feedobs file (FOF)
    max_body    ,& ! max. number of observations in NetCDF feedobs file (FOF)
    num_compute ,& ! number of compute PEs
    my_cart_id  ,& ! rank of this subdomain in the cartesian communicator
    icomm_cart  ,& ! communicator for the virtual cartesian topology
    imp_integers,& ! INTEGER type used in the model for MPI
    imp_reals      ! REAL    type used in the model for MPI

  CHARACTER (LEN=14)       , INTENT (IN)  :: &
    ydate_ini      ! reference date (e.g. start of the model run)
                   !                [yyyymmddhhmmss]
                   !                (year, month, day, hour, min., sec.)

  REAL    (KIND=wp)        , INTENT (IN)  :: &
    hversta     ,& ! start of verification period in 'model integration hours'
    hverend        ! end of verification period in 'model integration hours'

  ! following variables are only used for supplementary info in nc file header

  INTEGER (KIND=iintegers) , INTENT (IN)  :: &
    ie_tot      ,& ! number of grid points in zonal direction
    je_tot      ,& ! number of grid points in meridional direction
    ke_tot         ! number of grid points in vertical direction

  REAL    (KIND=wp)        , INTENT (IN)  :: &
    pollat      ,& ! latitude  of the rotated north pole (in degrees, N>0)
    pollon      ,& ! longitude of the rotated north pole (in degrees, E>0)
    dlat        ,& ! grid point distance in meridional direction (in degrees)
    dlon        ,& ! grid point distance in zonal      direction (in degrees)
    startlat_tot,& ! rotated latitude   \  of the lower left grid point of the
    startlon_tot   ! rotated longitude  /  total domain (in degrees, N,E>0)

  INTEGER (KIND=iintegers) , INTENT (IN)  :: &
    iveri_ens_member  ,& ! ensemble member ( -1: deterministic)
    nvers                ! exp. id + 16384*(class of model run)

  CHARACTER (LEN= * )      , INTENT (IN)  :: &
    ! NetCDF global attributes:
    yglatt_institution,& ! originating center name
    yglatt_source        ! program name and version
!   yglatt_title      ,& ! title string for the output
!   yglatt_references    ! URL, report etc.

  CHARACTER (LEN= * )      , INTENT (IN)  :: &
    yfdbkdir             ! directory name for feedobs file

  INTEGER (KIND=iintegers) , INTENT (INOUT)  :: &
    nexce_rep   ,& ! number of reports exceeding NetCDF feedobs file size
    nexce_bdy      ! number of obs     exceeding NetCDF feedobs file size

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND = iintegers) , SAVE :: &
    mboxold    = -9999    ! time [min] of last (latest) report which had been
                          ! written during the previous calls of the routine

  INTEGER (KIND = iintegers) ::  &
    nodrlni  ,nodrlnr  ,& ! dynamic length of the local ODR output buffer
    iob                ,& ! (loop) index for reports
    ircv               ,& ! rank of the receiving PE
    ioutleni ,ioutlenr ,& ! number of elements in local buffer
    ialeni   ,ialenr   ,& ! number of elements in global buffer (dim = ialenx+1)
    nxrep    ,n_rep    ,& ! number of reports to be printed
    nfcast             ,& ! number of forecast runs to compare with observations
    nzerr              ,& ! error status variable
    mpierr                ! for MPI error code

  LOGICAL                   :: &
    lvofoi                ! .TRUE  if obs increments are also written to VOF

  CHARACTER (LEN=22)       :: yroutine       ! name of this subroutine
  CHARACTER (LEN=80)       :: yerrmsg        ! error message

! Local arrays:
! ------------

  INTEGER (KIND = iintegers) , ALLOCATABLE  :: &
    iodrbufp  (:) ,& ! buffer containing integer ODR data of a single PE
    ioffsrt   (:)    ! sorted list of offsets of the reports in iodrbufa

  REAL (KIND = wp)           , ALLOCATABLE  :: &
    rodrbufp  (:)    ! buffer containing real ODR data of a single PE
    
  INTEGER (KIND = iintegers) , POINTER      :: &
    iodrbufa  (:)    ! buffer containing integer ODR data of a all PEs

  REAL (KIND = wp)           , POINTER      :: &
    rodrbufa  (:)    ! buffer containing real ODR data of a all PEs
!
!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine obs_org_cdfout_feedobs
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 0: Initialisations
!-------------------------------------------------------------------------------

  yroutine   = 'obs_org_cdfout_feedobs'
  ircv       = 0
  lvofoi     = (mruntyp >= 0)
  nfcast     = MAX( 0 , MIN( 1 , mruntyp+1 ) )

!-------------------------------------------------------------------------------
!  Section 1:  Write local ODR data to local 1-dim. (int / real) buffer arrays
!-------------------------------------------------------------------------------

  ! dynamic allocation of memory for the local real buffer
  ! ------------------------------------------------------
  nodrlni = 1
  nodrlnr = 1
  DO iob = 1, ntotsg
    IF (mosghd(iob,nhqcfw) >= 0) THEN
      IF (      (mosghd(iob,nhobtp) /= nairep)                                 &
          .AND. (mosghd(iob,nhobtp) /= nsatob)) THEN
        nodrlni = nodrlni  + mxchis + ilstidc + mxcbis
        nodrlnr = nodrlnr  + mxchra           + mxcbrs + nfcast *mxcmrs
      ELSE
        nodrlni = nodrlni  + mxchiu + ilstidc + mxcbiu
        nodrlnr = nodrlnr  + mxchra           + mxcbru + nfcast *mxcmru
      ENDIF
    ENDIF
  ENDDO
  DO iob = 1, ntotml
    IF (momlhd(iob,nhqcfw) >= 0) THEN
      nodrlni = nodrlni + mxchiu + ilstidc +         mxcbiu * momlhd(iob,nhnlev)
      nodrlnr = nodrlnr + mxchra           +(        mxcbru                    &
                                             +nfcast*mxcmru)* momlhd(iob,nhnlev)
    ENDIF
  ENDDO
  DO iob = 1, ntottv
    IF (motvhd(iob,nhqcfw) >= 0) THEN
      nodrlni = nodrlni + mxchiu + ilstidc +         mxcbiu * motvhd(iob,nhnlev)
      nodrlnr = nodrlnr + mxchra           +(        mxcbru                    &
                                             +nfcast*mxcmru)* motvhd(iob,nhnlev)
    ENDIF
  ENDDO
  DO iob = 1, ntotgp
    IF (mogphd(iob,nhqcfw) >= 0) THEN
      nodrlni = nodrlni  + mxchis + ilstidc + mxcbig
      nodrlnr = nodrlnr  + mxchra           + mxcbrg + nfcast *mxcmrg
    ENDIF
  ENDDO

  nzerr = 0
                   ALLOCATE ( iodrbufp (nodrlni)    , STAT = nzerr )
  IF (nzerr == 0)  ALLOCATE ( rodrbufp (nodrlnr)    , STAT = nzerr )
  IF (nzerr /= 0) THEN
    yerrmsg = 'ERROR in memory allocation for rodrbufp'
    CALL model_abort (my_cart_id, 19001, yerrmsg, yroutine)
  ENDIF    

  ! store local ODR data in 1-dim. (int / real) buffer arrays
  ! ---------------------------------------------------------

  CALL obs_neff_fill_buffer ( lrefend, mruntyp, hversta, hverend               &
                            , nodrlni, nodrlnr , mboxold, nxrep                &
                            , ioutleni, ioutlenr, iodrbufp, rodrbufp )
! =========================

  ! de-allocate smlbdy etc
  ! ----------------------
  DEALLOCATE ( ssgbdy  , STAT = nzerr )
  DEALLOCATE ( sgpbdy  , STAT = nzerr )
  DEALLOCATE ( stvbdy  , STAT = nzerr )
  DEALLOCATE ( smlbdy  , STAT = nzerr )
  DEALLOCATE ( dmlhed  , STAT = nzerr )

!-------------------------------------------------------------------------------
!  Section 2:  Communication:  Gather NetCDF feedobs file data at root node
!-------------------------------------------------------------------------------

  ! note : 'iodrbufa', 'rodrbufa' are allocated in obs_gather_buffer

  CALL obs_gather_buffer ( ioutleni, nodrlni, iodrbufp, ialeni, iodrbufa, ircv &
                         , num_compute, my_cart_id, icomm_cart, imp_integers )
! ======================

  CALL obs_gather_buffer ( ioutlenr, nodrlnr, rodrbufp, ialenr, rodrbufa, ircv &
                         , num_compute, my_cart_id, icomm_cart, imp_reals )
! ======================

  DEALLOCATE ( iodrbufp   , STAT = nzerr )
  DEALLOCATE ( rodrbufp   , STAT = nzerr )

  ! gather the number of reports from each node

  IF (num_compute > 1)  CALL global_values ( nxrep, 1, 'SUM', imp_integers     &
                                           , icomm_cart, ircv, yerrmsg, mpierr )
!                       ==================

  ! update 'mboxold'

  IF (num_compute > 1)  CALL global_values ( mboxold, 1, 'MAX', imp_integers   &
                                           , icomm_cart, -1, yerrmsg, mpierr )
!                       ==================

!-------------------------------------------------------------------------------
!  Section 3:  Make a sorted list of the reports
!              (independent from the domain decomposition)
!-------------------------------------------------------------------------------

  IF (my_cart_id == ircv) THEN

    ALLOCATE ( ioffsrt (nxrep + 1) , STAT = nzerr )
    IF (nzerr /= 0) THEN
      yerrmsg = 'ERROR in memory allocation for ioffsrt'
      CALL model_abort (my_cart_id, 19001, yerrmsg, yroutine)
    ENDIF    

    CALL obs_neff_sort_buffer ( nxrep, ialeni, iodrbufa, ie_tot, je_tot        &
                              , n_rep, ioffsrt )
!   =========================

!-------------------------------------------------------------------------------
!  Section 4:  Create the NetCDF feedback file 
!-------------------------------------------------------------------------------

    CALL obs_write_cdf_feedobs ( n_rep, ialeni, iodrbufa, ialenr, rodrbufa     &
                               , ioffsrt, lastpr, lrefend                      &
                               , ydate_ini, hversta, hverend                   &
                               , mruntyp, max_rep, max_body                    &
                               , ie_tot, je_tot, ke_tot, pollat, pollon        &
                               , dlat, dlon, startlat_tot, startlon_tot        &
                               , iveri_ens_member, nvers                       &
                               , yglatt_institution, yglatt_source, yfdbkdir   &
                               , nexce_rep, nexce_bdy, nzerr, yerrmsg )
!   ========================

    IF (nzerr /= 0)  CALL model_abort (my_cart_id, nzerr, yerrmsg, yroutine)
!                    ----------------
    DEALLOCATE ( ioffsrt    , STAT = nzerr ) 

  ENDIF  !  (my_cart_id == ircv)

  DEALLOCATE ( iodrbufa   , STAT = nzerr )
  DEALLOCATE ( rodrbufa   , STAT = nzerr )

!-------------------------------------------------------------------------------
! End of module procedure obs_org_cdfout_feedobs
!-------------------------------------------------------------------------------
END SUBROUTINE obs_org_cdfout_feedobs



!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_cdfout_feedobs" for storing reports in buffer
!-------------------------------------------------------------------------------

SUBROUTINE obs_neff_fill_buffer ( lrefend, mruntyp, hversta, hverend           &
                                , nodrlni, nodrlnr , mboxold, nxrep            &
                                , ioutleni, ioutlenr, iodrbufp, rodrbufp )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "src_obs_cdfout_feedobs" stores the reports,
!   which shall be written to the NetCDF feedobs file at the current timestep,
!   to 1-dimensional integer / real buffer arrays.
!
! Method:
!   The data are copied from the ODR for those reports which match the
!   time criteria and for which element 'nhqcfw' >= 0 .
!   The use of 'mboxold', 'mboxnow' prevents (re-produced) multi-level
!   aircraft reports from being written to the feedobs file more than once.
!   Initial version by C. Schraff, based on code written by Marek Lazanowicz
!   on 31.05.09 .
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
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

  LOGICAL                  , INTENT (IN)    :: &
    lrefend        ! .t.: verification reference time set to end of model run 
                   !      (for NUDGE runs only; for FORECAST (LETKF) or
                   !       NUDGECAST, veri_ref_time = start of model run)

  INTEGER (KIND=iintegers) , INTENT (IN)    :: &
    mruntyp     ,& ! =-1: print only observations ; =2: also print increments
    nodrlni     ,& ! length of local integer buffer array 'iodrbufp'
    nodrlnr        ! length of local real    buffer array 'rodrbufp'

  REAL    (KIND=wp)        , INTENT (IN)    :: &
    hversta     ,& ! start of verification period in 'model integration hours'
    hverend        ! end of verification period in 'model integration hours'

  INTEGER (KIND=iintegers) , INTENT (INOUT) :: &
    mboxold        ! time [min] of last (latest) report which had been
                   ! written during the previous calls of the routine

  INTEGER (KIND=iintegers) , INTENT (OUT)   :: &
    nxrep       ,& ! number of reports to be printed (i.e. in 'xodrbufp')
    ioutleni    ,& ! number of integer elements in local buffer array 'iodrbufp'
    ioutlenr    ,& ! number of real    elements in local buffer array 'iodrbufp'
    iodrbufp (nodrlni)   ! buffer containing integer ODR data from a single PE

  REAL    (KIND=wp)        , INTENT (OUT)   :: &
    rodrbufp (nodrlnr)   ! buffer containing real    ODR data from a single PE

! Local parameters:
! ----------------

  REAL    (KIND = wp)        , PARAMETER  :: &
    c60     =    60.0_wp        !

! Local scalars:
! -------------

  INTEGER (KIND = iintegers) ::  &
    irepleni ,ireplenr ,& ! length of the current report (int / real)
    iob                ,& ! (loop) indices for reports
    klev               ,& ! (loop) indices for levels (or variables)
    icl                ,& ! loop index for characters in a string
    iposi    ,iposr    ,& ! address  within a buffer
    nfcast             ,& ! number of forecast runs to compare with observations
    nrefmin            ,& ! diff. betw. verif. ref. time and model initial time
    mindiff            ,& ! diff. betw. synoptic and exact observation time
    mainflag ,mqcflag  ,& ! main flag word  /  first guess QC flag word
    mreptim            ,& ! time [min] of current report
    mboxnow            ,& ! time [min] of last (latest) report written currently
    icsg                  ! (vertical) level significance for all cloud groups

  REAL    (KIND = wp)        :: &
    hour_shift         ,& ! diff. betw. model initial time and verif. ref. time
    zhr                   ! hour

  LOGICAL                   :: &
    lvofoi             ,& ! .TRUE  if obs increments are also written to VOF
    lsurf                 ! .TRUE  if surface-level report

! Local arrays:
! ------------
!
!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine obs_neff_fill_buffer
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 0: Initialisations
!-------------------------------------------------------------------------------

  ! initialise parameters
  lvofoi     = (mruntyp >= 0)
  nfcast     = MAX( 0 , MIN( 1 , mruntyp+1 ) )
  zhr        = REAL (INT( hversta + epsy ), wp)
  mboxnow    = -9999

  ! difference between model initial time and verification reference time
  IF (      lrefend)  nrefmin  =  NINT( hverend *c60 )
  IF (.NOT. lrefend)  nrefmin  =  0
  hour_shift  =  - REAL( nrefmin,wp ) /c60

  ! initialise buffer arrays and counters
  iodrbufp (:) = imdi
  rodrbufp (:) = rmdi
  iposi    =  0
  iposr    =  0
  nxrep    =  0
  ioutleni =  0
  ioutlenr =  0

!-------------------------------------------------------------------------------
!  Section 1:  Single-level reports
!-------------------------------------------------------------------------------
  
!_______________________________________
single_level_report: DO iob = 1 , ntotsg
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  mreptim = NINT( (osghed(iob,nhtime) - zhr) * c60 )

  IF ((mosghd(iob,nhqcfw) >= 0) .AND. (osghed(iob,nhtime) >= hversta)          &
                                .AND. (osghed(iob,nhtime) <= hverend)          &
                                .AND. (mreptim > mboxold)) THEN
    nxrep    =  nxrep + 1
    lsurf  =        (mosghd(iob,nhobtp) /= nairep)                             &
              .AND. (mosghd(iob,nhobtp) /= nsatob)
    IF (lsurf) THEN
      irepleni  =  mxchis + ilstidc + mxcbis
      ireplenr  =  mxchra           + mxcbrs + nfcast *mxcmrs
    ELSE
      irepleni  =  mxchiu + ilstidc + mxcbiu
      ireplenr  =  mxchra           + mxcbru + nfcast *mxcmru
    ENDIF
    ioutleni  =  ioutleni + irepleni
    ioutlenr  =  ioutlenr + ireplenr
    mboxnow   =  MAX( mboxnow , mreptim )

!   header of single-level reports
!   ------------------------------
    iodrbufp(iposi+nchlen) = irepleni
    iodrbufp(iposi+nchler) = ireplenr
    IF (      lsurf) iodrbufp(iposi+nchfmt) = 2
    IF (.NOT. lsurf) iodrbufp(iposi+nchfmt) = 1

    iodrbufp(iposi+nchtim) = NINT( (osghed(iob,nhtime) + hour_shift) *c60 )
    ! assume that exact and nominal (synoptic) obs times differ by < 12 hrs
    mindiff               =  MOD( mosghd(iob,nhsyhr) ,100 )                    &
                            -     mosghd(iob,nhhrmn) /100
    IF (mindiff < -12)  mindiff = mindiff + 24
    IF (mindiff >  12)  mindiff = mindiff - 24
    mindiff               = mindiff *60  -  MOD( mosghd(iob,nhhrmn) ,100 )
    iodrbufp(iposi+nchsyn) = iodrbufp(iposi+nchtim) + mindiff
    IF (osghed(iob,nhtddb) > rmdich)                                           &
      iodrbufp(iposi+nchtdb) = NINT( (osghed(iob,nhtddb) + hour_shift) *c60 )
    IF (osghed(iob,nhalt ) > rmdich)                                           &
      iodrbufp(iposi+nchalt) = NINT( osghed(iob,nhalt ) )
    IF (osghed(iob,nhsurf) > rmdich)                                           &
      iodrbufp(iposi+nchsfc) = NINT( osghed(iob,nhsurf) )
    iodrbufp(iposi+nchobt) = mosghd(iob,nhobtp)
    iodrbufp(iposi+nchcdt) = mosghd(iob,nhcode)
    iodrbufp(iposi+nchsch) = mosghd(iob,nhschr)
    !   not required to include entry 'nhpass' into entry 'nchsch'
!   iodrbufp(iposi+nchsch) = insert( mosghd(iob,nhschr)                        &
!                                  , mosghd(iob,nhpass), nvpabp )
    IF (      (mosghd(iob,nhpass) == 1)                                        &
        .AND. (.NOT. BTEST( mosghd(iob,nhflag), FL_MERGE )))                   &
      mosghd (iob,nhflag) = IBSET ( mosghd(iob,nhflag), FL_MERGE )
    iodrbufp(iposi+nchflg) = mosghd(iob,nhflag)
!   SELECT CASE ( mosghd(iob,nhpass) )
!     CASE ( 0 ) ; iodrbufp(iposi+nchpas) = ST_ACTIVE
!     CASE ( 1 ) ; iodrbufp(iposi+nchpas) = ST_MERGED
!     CASE ( 2 ) ; iodrbufp(iposi+nchpas) = ST_PASSIVE
!     CASE ( 3 ) ; iodrbufp(iposi+nchpas) = ST_REJECTED
!   END SELECT
    iodrbufp(iposi+nchpas) = ncfeedobs_status ( mosghd(iob,nhpass)             &
                                              , mosghd(iob,nhflag) )
    iodrbufp(iposi+nchi_t) = mosghd(iob,nhitot)
    iodrbufp(iposi+nchj_t) = mosghd(iob,nhjtot)
    iodrbufp(iposi+nchcor) = mosghd(iob,nhcorr)
    iodrbufp(iposi+nchcat) = mosghd(iob,nhcat )
    iodrbufp(iposi+nchcas) = mosghd(iob,nhcats)
    iodrbufp(iposi+nchcen) = mosghd(iob,nhcent)
    iodrbufp(iposi+nchide) = mosghd(iob,nhstid)
    iodrbufp(iposi+nchtyp) = mosghd(iob,nhstyp)
    IF (.NOT. lsurf) THEN
      iodrbufp(iposi+nchnlv) = 1
      IF (mosghd(iob,nhobtp) == nairep) THEN
        IF (     (IBITS( mosghd(iob,nhschr),nvapbp,nvapoc ) /= nibits(nvapoc)) &
           .AND. (BTEST( mosghd(iob,nhschr),nvapbp+nvapoc-1)))                 &
          iodrbufp(iposi+nchtrc) = 18
      ENDIF
!     iodrbufp(iposi+nchrad) = mosghd(iob,nhrad )  !  <-- does not exist
!     iodrbufp(iposi+nchna4) = mosghd(iob,nhna4 )  !  <-- does not exist
      iposi  =  iposi + mxchiu
    ELSE
      iposi  =  iposi + mxchis
    ENDIF

    rodrbufp(iposr+nchlon) = osghed(iob,nhilon)
    rodrbufp(iposr+nchlat) = osghed(iob,nhjlat)
    rodrbufp(iposr+nchsol) = osghed(iob,nhsolz)
    iposr  =  iposr + mxchra

    DO icl = 1 , ilstidc
      iodrbufp(iposi+icl) = ICHAR( yosghd(iob) (icl:icl) )
    ENDDO
    iposi  =  iposi + ilstidc

!   body of single-level reports
!   ----------------------------
    iodrbufp(iposi+ncberf) = mosgbd(iob,nbserr)
    mainflag = mosgbd(iob,nbsflg)
    ! upper-air obs pressure level below model surface
    IF ((.NOT. lsurf) .AND. (mosgbd(iob,nbsqcf)/ 8 == 1)) THEN
      mainflag = IBSET ( mainflag, nvfubp+nvfbps(4) )
      mainflag = IBSET ( mainflag, nvftbp+nvfbps(4) )
      mainflag = IBSET ( mainflag, nvfqbp+nvfbps(4) )
      mainflag = IBSET ( mainflag, nvfzbp+nvfbps(4) )
      iodrbufp(iposi+ncberf) = 0
    ENDIF
    ! height differ. betw. surface pressure obs and lowest model level too large
    IF (mosghd(iob,nhqcfw)/ 2 == 1)                                            &
      mainflag = IBSET ( mainflag, nvfzbp+nvfbps(4) )
    ! surface pressure first guess check flag (bit) from ODR header is .NOT.
    ! added to the f.g. check flag word, because in the ODR, this bit IS already
    ! set in the body for single-level reports
!   iodrbufp(iposi+ncbqcf) = ireplace1( mosgbd(iob,nbsqcf), nvrz               &
!                                     , mosghd(iob,nhqcfw), 0 )
    iodrbufp(iposi+ncbqcf) = mosgbd(iob,nbsqcf)
    iodrbufp(iposi+ncbmfw) = mainflag
    IF (.NOT. lsurf) THEN
!     iodrbufp(iposi+ncblsg) = mosgbd(iob,nbslsg)
      iodrbufp(iposi+ncblsg) = IBSET ( 0, LS_SIGN )
      iodrbufp(iposi+ncbsnr) = imdi
      iodrbufp(iposi+ncbtur) = mosgbd(iob,nbstur)
      iposi  =  iposi + mxcbiu
    ELSE
!     iodrbufp(iposi+ncblsg) = IBSET ( 0, LS_SURFACE )
      IF (mosghd(iob,nhobtp) == nsynop) THEN
        ! SYNOP: pressure code
        iodrbufp(iposi+ncblsg) = mosgbd(iob,nbslid)
      ELSE
        ! DRIBU, SCATT, surface level of TEMP or PILOT:
        ! surface-level bit is set in 'nbslid' in subr. 'obs_cdf_store_singlev'
        iodrbufp(iposi+ncblsg) = imdi
        IF (BTEST( mosgbd(iob,nbslid), nvlidp(7) ))                            &
          iodrbufp(iposi+ncblsg) = IBSET ( 0, LS_SURFACE )
      ENDIF
      iodrbufp(iposi+ncbwwe) = mosgbd(iob,nbswwe)
      IF (mosgbd(iob,nbsclg) /= imdi)                                          &
        iodrbufp(iposi+ncbclg) = ISHFT( IBITS( mosgbd(iob,nbsclg), nxclbp      &
                                             , nxcloc+3*nxctoc ) , nxclbp )
      IF (mosgbd(iob,nbscl1) /= imdi)                                          &
        iodrbufp(iposi+ncbcl1) = ISHFT( IBITS( mosgbd(iob,nbscl1), nxclbp      &
                                             , nxcloc+nxctoc+nxbsoc ) , nxclbp )
      IF (mosgbd(iob,nbscl2) /= imdi)                                          &
        iodrbufp(iposi+ncbcl2) = ISHFT( IBITS( mosgbd(iob,nbscl2), nxclbp      &
                                             , nxcloc+nxctoc+nxbsoc ) , nxclbp )
      IF (mosgbd(iob,nbscl3) /= imdi)                                          &
        iodrbufp(iposi+ncbcl3) = ISHFT( IBITS( mosgbd(iob,nbscl3), nxclbp      &
                                             , nxcloc+nxctoc+nxbsoc ) , nxclbp )
      IF (mosgbd(iob,nbscl4) /= imdi)                                          &
        iodrbufp(iposi+ncbcl4) = ISHFT( IBITS( mosgbd(iob,nbscl4), nxclbp      &
                                             , nxcloc+nxctoc+nxbsoc ) , nxclbp )
      icsg                   = imdi
      CALL MVBITS ( mosgbd(iob,nbsclg), nxsgbp, nxsgoc, icsg, ncsgbp(0) )
      CALL MVBITS ( mosgbd(iob,nbscl1), nxsibp, nxsgoc, icsg, ncsgbp(1) )
      CALL MVBITS ( mosgbd(iob,nbscl2), nxsibp, nxsgoc, icsg, ncsgbp(2) )
      CALL MVBITS ( mosgbd(iob,nbscl3), nxsibp, nxsgoc, icsg, ncsgbp(3) )
      CALL MVBITS ( mosgbd(iob,nbscl4), nxsibp, nxsgoc, icsg, ncsgbp(4) )
      iodrbufp(iposi+ncbcsg) = icsg
      iodrbufp(iposi+ncbttr) = mosgbd(iob,nbsttr)
      iposi  =  iposi + mxcbis
    ENDIF

    IF ((osgbdy(iob,nbsu) > rmdich) .AND. (osgbdy(iob,nbsv) > rmdich)) THEN
      rodrbufp(iposr+ncbu  ) = osgbdy(iob,nbsu  )
      rodrbufp(iposr+ncbv  ) = osgbdy(iob,nbsv  )
    ENDIF
    rodrbufp(iposr+ncbt  ) = osgbdy(iob,nbst  )
    rodrbufp(iposr+ncbrh ) = osgbdy(iob,nbsrh )
    rodrbufp(iposr+ncbp  ) = osgbdy(iob,nbsp  )
    rodrbufp(iposr+ncbz  ) = osgbdy(iob,nbsz  )
    rodrbufp(iposr+ncbuer) = osgbdy(iob,nbsuer)
    rodrbufp(iposr+ncbter) = osgbdy(iob,nbster)
    rodrbufp(iposr+ncbqer) = osgbdy(iob,nbsqer)
    rodrbufp(iposr+ncbzer) = osgbdy(iob,nbszer)
    rodrbufp(iposr+ncbdrh) = osgbdy(iob,nbsdrh)
    IF (.NOT. lsurf) THEN
!     rodrbufp(iposr+ncbw  ) = rmdi
      rodrbufp(iposr+ncbfgv) = osgbdy(iob,nbsfgv)
!     rodrbufp(iposr+ncbuac) = rmdi
      iposr  =  iposr + mxcbru
    ELSE
      rodrbufp(iposr+ncbcbs) = osgbdy(iob,nbscbs)
      rodrbufp(iposr+ncbcl ) = osgbdy(iob,nbscl )
      rodrbufp(iposr+ncbcm ) = osgbdy(iob,nbscm )
      rodrbufp(iposr+ncbch ) = osgbdy(iob,nbsch )
      rodrbufp(iposr+ncbct ) = osgbdy(iob,nbsct )
      rodrbufp(iposr+ncbvis) = osgbdy(iob,nbsvis)
      rodrbufp(iposr+ncbff ) = osgbdy(iob,nbsff )
      rodrbufp(iposr+ncbdd ) = osgbdy(iob,nbsdd )
      rodrbufp(iposr+ncbtd ) = osgbdy(iob,nbstd )
      rodrbufp(iposr+ncbhsw) = osgbdy(iob,nbshsw)
      rodrbufp(iposr+ncbpst) = osgbdy(iob,nbspst)
      rodrbufp(iposr+ncbrr1) = osgbdy(iob,nbsrr1)
      rodrbufp(iposr+ncbrr3) = osgbdy(iob,nbsrr3)
      rodrbufp(iposr+ncbrr6) = osgbdy(iob,nbsrr6)
      rodrbufp(iposr+ncbr12) = osgbdy(iob,nbsr12)
      rodrbufp(iposr+ncbr24) = osgbdy(iob,nbsr24)
      rodrbufp(iposr+ncbfg1) = osgbdy(iob,nbsfg1)
      rodrbufp(iposr+ncbfg3) = osgbdy(iob,nbsfg3)
      rodrbufp(iposr+ncbfg6) = osgbdy(iob,nbsfg6)
      rodrbufp(iposr+ncbtn ) = osgbdy(iob,nbstn )
      rodrbufp(iposr+ncbtx ) = osgbdy(iob,nbstx )
      rodrbufp(iposr+ncbrad) = osgbdy(iob,nbsrad)
      rodrbufp(iposr+ncbrdd) = osgbdy(iob,nbsrdd)
      rodrbufp(iposr+ncbrdt) = osgbdy(iob,nbsrdt)
      iposr  =  iposr + mxcbrs
    ENDIF

!   model-derived simulated observations for single-level reports
!   -------------------------------------------------------------
    IF (lvofoi) THEN
!     IF ((osgbdy(iob,nbsu ) > rmdich) .AND. (dsgbdy(nso_u ,iob) > rmdich))    &
      rodrbufp(iposr+ncmu  ) = ssgbdy(nso_u ,iob)
      rodrbufp(iposr+ncmv  ) = ssgbdy(nso_v ,iob)
      rodrbufp(iposr+ncmt  ) = ssgbdy(nso_t ,iob)
      rodrbufp(iposr+ncmrh ) = ssgbdy(nso_rh,iob)
      IF (lsurf) THEN
        rodrbufp(iposr+ncmp  ) = ssgbdy(nso_p  ,iob)
        rodrbufp(iposr+ncmcbs) = ssgbdy(nso_cbs,iob)
        rodrbufp(iposr+ncmcl ) = frac2octas( ssgbdy(nso_cl ,iob) )
        rodrbufp(iposr+ncmcm ) = frac2octas( ssgbdy(nso_cm ,iob) )
        rodrbufp(iposr+ncmch ) = frac2octas( ssgbdy(nso_ch ,iob) )
        rodrbufp(iposr+ncmct ) = frac2octas( ssgbdy(nso_ct ,iob) )
!       rodrbufp(iposr+ncmvis) = ssgbdy(nso_vis,iob)
        rodrbufp(iposr+ncmff ) = ssgbdy(nso_ff ,iob)
        rodrbufp(iposr+ncmdd ) = ssgbdy(nso_dd ,iob)
        rodrbufp(iposr+ncmtd ) = ssgbdy(nso_td ,iob)
!       rodrbufp(iposr+ncmhsw) = ssgbdy(nso_hsw,iob)
        iposr  =  iposr + mxcmrs
      ELSE
        iposr  =  iposr + mxcmru
      ENDIF
    ENDIF

  ENDIF
  
! set the flag for having the printing of the current report done
! ---------------------------------------------------------------

  IF (mosghd(iob,nhqcfw) >= 0)  mosghd (iob,nhqcfw) = - 99

ENDDO single_level_report
!~~~~~~~~~~~~~~~~~~~~~~~~

!-------------------------------------------------------------------------------
!  Section 2:  (Conventional) multi-level reports
!-------------------------------------------------------------------------------

!______________________________________
multi_level_report: DO iob = 1 , ntotml
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  mreptim = NINT( (omlhed(iob,nhtime) - zhr) * c60 )

  IF ((momlhd(iob,nhqcfw) >= 0) .AND. (omlhed(iob,nhtime) >= hversta)          &
                                .AND. (omlhed(iob,nhtime) <= hverend)          &
                                .AND. (mreptim > mboxold)) THEN
      ! the 'mboxold' condition avoids re-printing of multi-level
      ! aircraft reports just after having re-created them again
    nxrep    =  nxrep + 1
    irepleni  =  mxchiu + ilstidc +            mxcbiu * momlhd(iob,nhnlev)
    ireplenr  =  mxchra           + (          mxcbru                          &
                                     + nfcast *mxcmru)* momlhd(iob,nhnlev)
    ioutleni  =  ioutleni + irepleni
    ioutlenr  =  ioutlenr + ireplenr
    mboxnow   =  MAX( mboxnow , mreptim )

!   header of multi-level reports
!   -----------------------------
    iodrbufp(iposi+nchlen) = irepleni
    iodrbufp(iposi+nchler) = ireplenr
    iodrbufp(iposi+nchfmt) = 1

    iodrbufp(iposi+nchtim) = NINT( (omlhed(iob,nhtime) + hour_shift) *c60 )
    ! assume that exact and nominal (synoptic) obs times differ by < 12 hrs
    mindiff               =  MOD( momlhd(iob,nhsyhr) ,100 )                    &
                            -     momlhd(iob,nhhrmn) /100
    IF (mindiff < -12)  mindiff = mindiff + 24
    IF (mindiff >  12)  mindiff = mindiff - 24
    mindiff               = mindiff *60  -  MOD( momlhd(iob,nhhrmn) ,100 )
    iodrbufp(iposi+nchsyn) = iodrbufp(iposi+nchtim) + mindiff
    IF (omlhed(iob,nhtddb) > rmdich)                                           &
      iodrbufp(iposi+nchtdb) = NINT( (omlhed(iob,nhtddb) + hour_shift) *c60 )
    IF (omlhed(iob,nhalt ) > rmdich)                                           &
      iodrbufp(iposi+nchalt) = NINT( omlhed(iob,nhalt ) )
    IF (omlhed(iob,nhsurf) > rmdich)                                           &
      iodrbufp(iposi+nchsfc) = NINT( omlhed(iob,nhsurf) )
    iodrbufp(iposi+nchobt) = momlhd(iob,nhobtp)
    iodrbufp(iposi+nchcdt) = momlhd(iob,nhcode)
    iodrbufp(iposi+nchsch) = momlhd(iob,nhschr)
    !   not required to include entry 'nhpass' into entry 'nchsch'
!   iodrbufp(iposi+nchsch) = insert( momlhd(iob,nhschr)                        &
!                                  , momlhd(iob,nhpass), nvpabp )
    iodrbufp(iposi+nchflg) = momlhd(iob,nhflag)
    iodrbufp(iposi+nchpas) = ncfeedobs_status ( momlhd(iob,nhpass)             &
                                              , momlhd(iob,nhflag) )
    iodrbufp(iposi+nchi_t) = momlhd(iob,nhitot)
    iodrbufp(iposi+nchj_t) = momlhd(iob,nhjtot)
    iodrbufp(iposi+nchcor) = momlhd(iob,nhcorr)
    iodrbufp(iposi+nchcat) = momlhd(iob,nhcat )
    iodrbufp(iposi+nchcas) = momlhd(iob,nhcats)
    iodrbufp(iposi+nchcen) = momlhd(iob,nhcent)
    iodrbufp(iposi+nchide) = momlhd(iob,nhstid)
    iodrbufp(iposi+nchtyp) = momlhd(iob,nhrtyp)
    iodrbufp(iposi+nchnlv) = momlhd(iob,nhnlev)
    iodrbufp(iposi+nchtrc) = momlhd(iob,nhtrac)
    iodrbufp(iposi+nchrad) = momlhd(iob,nhrad )
    iodrbufp(iposi+nchna4) = momlhd(iob,nhna4 )
    IF (momlhd(iob,nhobtp) == nairep) THEN
      IF (      (IBITS( momlhd(iob,nhschr),nvapbp,nvapoc ) /= nibits(nvapoc))  &
          .AND. (BTEST( momlhd(iob,nhschr),nvapbp+nvapoc-1)))                  &
        iodrbufp(iposi+nchtrc) = 18
    ENDIF
    iposi  =  iposi + mxchiu

    rodrbufp(iposr+nchlon) = omlhed(iob,nhilon)
    rodrbufp(iposr+nchlat) = omlhed(iob,nhjlat)
    rodrbufp(iposr+nchsol) = omlhed(iob,nhsolz)

    iposr  =  iposr + mxchra

    DO icl = 1 , ilstidc
      iodrbufp(iposi+icl) = ICHAR( yomlhd(iob) (icl:icl) )
    ENDDO
    iposi  =  iposi + ilstidc

!   body of multi-level reports
!   ---------------------------

    Level: DO klev = 1 , momlhd(iob,nhnlev)

      iodrbufp(iposi+ncberf) = momlbd(iob,klev,nbterr)
      iodrbufp(iposi+ncbqcf) = momlbd(iob,klev,nbtqcf)
      iodrbufp(iposi+ncbmfw) = momlbd(iob,klev,nbtflg)
      iodrbufp(iposi+ncblsg) = momlbd(iob,klev,nbtlsg)
      IF (omlbdy(iob,klev,nbtsnr) > rmdich)                                    &
        iodrbufp(iposi+ncbsnr) = NINT( omlbdy(iob,klev,nbtsnr) )
      iodrbufp(iposi+ncbtur) = imdi
      iposi  =  iposi + mxcbiu

      IF (      (omlbdy(iob,klev,nbtu) > rmdich)                               &
          .AND. (omlbdy(iob,klev,nbtv) > rmdich)) THEN
        rodrbufp(iposr+ncbu  ) = omlbdy(iob,klev,nbtu )
        rodrbufp(iposr+ncbv  ) = omlbdy(iob,klev,nbtv )
      ENDIF
      rodrbufp(iposr+ncbt  ) = omlbdy(iob,klev,nbtt  )
      rodrbufp(iposr+ncbrh ) = omlbdy(iob,klev,nbtrh )
      rodrbufp(iposr+ncbp  ) = omlbdy(iob,klev,nbtp  )
      rodrbufp(iposr+ncbz  ) = omlbdy(iob,klev,nbtz  )
      rodrbufp(iposr+ncbuer) = omlbdy(iob,klev,nbtuer)
      rodrbufp(iposr+ncbter) = omlbdy(iob,klev,nbtter)
      rodrbufp(iposr+ncbqer) = omlbdy(iob,klev,nbtqer)
      rodrbufp(iposr+ncbzer) = omlbdy(iob,klev,nbtzer)
      rodrbufp(iposr+ncbdrh) = omlbdy(iob,klev,nbtdrh)
      rodrbufp(iposr+ncbw  ) = omlbdy(iob,klev,nbtw  )
!     rodrbufp(iposr+ncbfgv) = rmdi
      rodrbufp(iposr+ncbuac) = omlbdy(iob,klev,nbtuac)
      iposr  =  iposr + mxcbru

!     momlbd(iob,klev,nbtqcf) = -2
    ENDDO Level

!   model-derived simulated observations for multi-level reports
!   ------------------------------------------------------------
    IF (lvofoi) THEN
      DO klev = 1 , momlhd(iob,nhnlev)
        IF (smlbdy(nso_u ,klev,iob) > rmdich)                                  &
          rodrbufp(iposr+ncmu )  =  smlbdy(nso_u ,klev,iob)
        IF (smlbdy(nso_v ,klev,iob) > rmdich)                                  &
          rodrbufp(iposr+ncmv )  =  smlbdy(nso_v ,klev,iob)
        IF (smlbdy(nso_t ,klev,iob) > rmdich)                                  &
          rodrbufp(iposr+ncmt )  =  smlbdy(nso_t ,klev,iob)
        IF (smlbdy(nso_rh,klev,iob) > rmdich)                                  &
          rodrbufp(iposr+ncmrh)  =  smlbdy(nso_rh,klev,iob)
        ! multi-level reports: element 'nso_p' contains height increment
        IF (smlbdy(nso_p ,klev,iob) > rmdich)                                  &
          rodrbufp(iposr+ncmp )  =  smlbdy(nso_p ,klev,iob)
        iposr  =  iposr + mxcmru
      ENDDO
    ENDIF

  ENDIF

! set the flag for having the printing of the current report done
! ---------------------------------------------------------------

  IF (momlhd(iob,nhqcfw) >= 0)  momlhd (iob,nhqcfw) = - 99

ENDDO multi_level_report
!~~~~~~~~~~~~~~~~~~~~~~~

!-------------------------------------------------------------------------------
!  Section 3:  Satellite retrievals
!-------------------------------------------------------------------------------

!__________________________________
sat_retrievals: DO iob = 1 , ntottv
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  mreptim = NINT( (otvhed(iob,nhtime) - zhr) * c60 )

  IF ((motvhd(iob,nhqcfw) >= 0) .AND. (otvhed(iob,nhtime) >= hversta)          &
                                .AND. (otvhed(iob,nhtime) <= hverend)          &
                                .AND. (mreptim > mboxold)) THEN

    nxrep    =  nxrep + 1
    irepleni  =  mxchiu + ilstidc +            mxcbiu * motvhd(iob,nhnlev)
    ireplenr  =  mxchra           + (          mxcbru                          &
                                     + nfcast *mxcmru)* motvhd(iob,nhnlev)
    ioutleni  =  ioutleni + irepleni
    ioutlenr  =  ioutlenr + ireplenr
    mboxnow   =  MAX( mboxnow , mreptim )

!   header of satellite retrievals
!   ------------------------------
    iodrbufp(iposi+nchlen) = irepleni
    iodrbufp(iposi+nchler) = ireplenr
    iodrbufp(iposi+nchfmt) = 1

    iodrbufp(iposi+nchtim) = NINT( (otvhed(iob,nhtime) + hour_shift) *c60 )
    ! assume that exact and nominal (synoptic) obs times differ by < 12 hrs
    mindiff               =  MOD( motvhd(iob,nhsyhr) ,100 )                    &
                            -     motvhd(iob,nhhrmn) /100
    IF (mindiff < -12)  mindiff = mindiff + 24
    IF (mindiff >  12)  mindiff = mindiff - 24
    mindiff               = mindiff *60  -  MOD( motvhd(iob,nhhrmn) ,100 )
    iodrbufp(iposi+nchsyn) = iodrbufp(iposi+nchtim) + mindiff
    IF (otvhed(iob,nhtddb) > rmdich)                                           &
      iodrbufp(iposi+nchtdb) = NINT( (otvhed(iob,nhtddb) + hour_shift) *c60 )
    IF (otvhed(iob,nhalt ) > rmdich)                                           &
      iodrbufp(iposi+nchalt) = NINT( otvhed(iob,nhalt ) )
    IF (otvhed(iob,nhsurf) > rmdich)                                           &
      iodrbufp(iposi+nchsfc) = NINT( otvhed(iob,nhsurf) )
    iodrbufp(iposi+nchobt) = motvhd(iob,nhobtp)
    iodrbufp(iposi+nchcdt) = motvhd(iob,nhcode)
    iodrbufp(iposi+nchsch) = motvhd(iob,nhschr)
    iodrbufp(iposi+nchflg) = motvhd(iob,nhflag)
    iodrbufp(iposi+nchpas) = ncfeedobs_status ( motvhd(iob,nhpass)             &
                                              , motvhd(iob,nhflag) )
    iodrbufp(iposi+nchi_t) = motvhd(iob,nhitot)
    iodrbufp(iposi+nchj_t) = motvhd(iob,nhjtot)
    iodrbufp(iposi+nchcor) = motvhd(iob,nhcorr)
    iodrbufp(iposi+nchcat) = motvhd(iob,nhcat )
    iodrbufp(iposi+nchcas) = motvhd(iob,nhcats)
    iodrbufp(iposi+nchcen) = motvhd(iob,nhcent)
    iodrbufp(iposi+nchide) = motvhd(iob,nhstid)
    iodrbufp(iposi+nchtyp) = motvhd(iob,nhrtyp)
    iodrbufp(iposi+nchnlv) = motvhd(iob,nhnlev)
!   iodrbufp(iposi+nchtrc) = motvhd(iob,nhtrac)  !  <-- does not exist
!   iodrbufp(iposi+nchrad) = motvhd(iob,nhrad )  !  <-- does not exist
!   iodrbufp(iposi+nchna4) = motvhd(iob,nhna4 )  !  <-- does not exist
    iposi  =  iposi + mxchiu

    rodrbufp(iposr+nchlon) = otvhed(iob,nhilon)
    rodrbufp(iposr+nchlat) = otvhed(iob,nhjlat)
    rodrbufp(iposr+nchsol) = otvhed(iob,nhsolz)
    iposr  =  iposr + mxchra

    DO icl = 1 , ilstidc
      iodrbufp(iposi+icl) = ICHAR( yotvhd(iob) (icl:icl) )
    ENDDO
    iposi  =  iposi + ilstidc

!   body of satellite retrievals
!   ----------------------------

    Sat_level: DO klev = 1 , motvhd(iob,nhnlev)

!     iodrerr = 0
!     IF (otvbdy(iob,klev,nbvt  ) > rmdich) iodrerr = iodrerr + 2
!     IF (otvbdy(iob,klev,nbvrh ) > rmdich) iodrerr = iodrerr + 4
!     iodrbufp(iposi+ncberf) = iodrerr
      iodrbufp(iposi+ncberf) = motvbd(iob,klev,nbverr)
      iodrbufp(iposi+ncbqcf) = motvhd(iob     ,nhvqcf)
      iodrbufp(iposi+ncbmfw) = motvbd(iob,klev,nbvflg)
      iodrbufp(iposi+ncblsg) = 0
      iodrbufp(iposi+ncbsnr) = imdi
      iodrbufp(iposi+ncbtur) = imdi
      iposi  =  iposi + mxcbiu

!     rodrbufp(iposr+ncbu  ) = rmdi
!     rodrbufp(iposr+ncbv  ) = rmdi
      rodrbufp(iposr+ncbt  ) = otvbdy(iob,klev,nbvt  )
      rodrbufp(iposr+ncbrh ) = otvbdy(iob,klev,nbvrh )
      rodrbufp(iposr+ncbp  ) = otvbdy(iob,klev,nbvp  )
!     rodrbufp(iposr+ncbz  ) = rmdi
!     rodrbufp(iposr+ncbuer) = rmdi
      rodrbufp(iposr+ncbter) = rmdi
      rodrbufp(iposr+ncbqer) = rmdi
!     rodrbufp(iposr+ncbzer) = rmdi
      rodrbufp(iposr+ncbdrh) = c0
!     rodrbufp(iposr+ncbw  ) = rmdi
!     rodrbufp(iposr+ncbfgv) = rmdi
!     rodrbufp(iposr+ncbuac) = rmdi
      iposr  =  iposr + mxcbru

    ENDDO Sat_level

!   model-derived simulated observations for satellite retrievals
!   -------------------------------------------------------------
    IF (lvofoi) THEN
      DO klev = 1 , motvhd(iob,nhnlev)
!       rodrbufp(iposr+ncmu ) = rmdi
!       rodrbufp(iposr+ncmv ) = rmdi
        IF (stvbdy(nso_t ,klev,iob) > rmdich)                                  &
          rodrbufp(iposr+ncmt )  =  stvbdy(nso_t ,klev,iob)
        IF (stvbdy(nso_rh,klev,iob) > rmdich)                                  &
          rodrbufp(iposr+ncmrh)  =  stvbdy(nso_rh,klev,iob)
!       rodrbufp(iposr+ncmp ) = rmdi
        iposr  =  iposr + mxcmru
      ENDDO
    ENDIF

  ENDIF

! set the flag for having the printing of the current report done
! ---------------------------------------------------------------

  IF (motvhd(iob,nhqcfw) >= 0)  motvhd (iob,nhqcfw) = - 99

ENDDO sat_retrievals
!~~~~~~~~~~~~~~~~~~~

!-------------------------------------------------------------------------------
!  Section 4:  GNSS zenith path delay (IWV) reports
!-------------------------------------------------------------------------------

!______________________________
gps_report: DO iob = 1 , ntotgp
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  mreptim = NINT( (ogphed(iob,nhtime) - zhr) * c60 )

  IF ((mogphd(iob,nhqcfw) >= 0) .AND. (ogphed(iob,nhtime) >= hversta)          &
                                .AND. (ogphed(iob,nhtime) <= hverend)          &
                                .AND. (mreptim > mboxold)) THEN
    nxrep    =  nxrep + 1
    irepleni  =  mxchis + ilstidc + mxcbig
    ireplenr  =  mxchra           + mxcbrg + nfcast *mxcmrg
    ioutleni  =  ioutleni + irepleni
    ioutlenr  =  ioutlenr + ireplenr
    mboxnow   =  MAX( mboxnow , mreptim )

!   header of ground-based GPS reports
!   ----------------------------------
    iodrbufp(iposi+nchlen) = irepleni
    iodrbufp(iposi+nchler) = ireplenr
    iodrbufp(iposi+nchfmt) = 3

    iodrbufp(iposi+nchtim) = NINT( (ogphed(iob,nhtime) + hour_shift) *c60 )
    ! assume that exact and nominal (synoptic) obs times differ by < 12 hrs
    mindiff               =  MOD( mogphd(iob,nhsyhr) ,100 )                    &
                            -     mogphd(iob,nhhrmn) /100
    IF (mindiff < -12)  mindiff = mindiff + 24
    IF (mindiff >  12)  mindiff = mindiff - 24
    mindiff               = mindiff *60  -  MOD( mogphd(iob,nhhrmn) ,100 )
    iodrbufp(iposi+nchsyn) = iodrbufp(iposi+nchtim) + mindiff
    IF (ogphed(iob,nhtddb) > rmdich)                                           &
      iodrbufp(iposi+nchtdb) = NINT( (ogphed(iob,nhtddb) + hour_shift) *c60 )
    IF (ogphed(iob,nhalt ) > rmdich)                                           &
      iodrbufp(iposi+nchalt) = NINT( ogphed(iob,nhalt ) )
    IF (ogphed(iob,nhsurf) > rmdich)                                           &
      iodrbufp(iposi+nchsfc) = NINT( ogphed(iob,nhsurf) )
    iodrbufp(iposi+nchobt) = mogphd(iob,nhobtp)
    iodrbufp(iposi+nchcdt) = mogphd(iob,nhcode)
    iodrbufp(iposi+nchsch) = mogphd(iob,nhschr)
    iodrbufp(iposi+nchflg) = mogphd(iob,nhflag)
    iodrbufp(iposi+nchpas) = ncfeedobs_status ( mogphd(iob,nhpass)             &
                                              , mogphd(iob,nhflag) )
    iodrbufp(iposi+nchi_t) = mogphd(iob,nhitot)
    iodrbufp(iposi+nchj_t) = mogphd(iob,nhjtot)
    iodrbufp(iposi+nchcor) = mogphd(iob,nhcorr)
    iodrbufp(iposi+nchcat) = mogphd(iob,nhcat )
    iodrbufp(iposi+nchcas) = mogphd(iob,nhcats)
    iodrbufp(iposi+nchcen) = mogphd(iob,nhcent)
    iodrbufp(iposi+nchide) = mogphd(iob,nhstid)
    iodrbufp(iposi+nchtyp) = imdi
    ! in the future, WMO 002020 (MSACL: satellite classification) could be used
!   iodrbufp(iposi+nchtyp) = mogphd(iob,nhstyp)
    iposi  =  iposi + mxchis

    rodrbufp(iposr+nchlon) = ogphed(iob,nhilon)
    rodrbufp(iposr+nchlat) = ogphed(iob,nhjlat)
    rodrbufp(iposr+nchsol) = ogphed(iob,nhsolz)
    iposr  =  iposr + mxchra

    DO icl = 1 , ilstidc
      iodrbufp(iposi+icl) = ICHAR( yogphd(iob) (icl:icl) )
    ENDDO
    iposi  =  iposi + ilstidc

!   body of ground-based GPS reports
!   --------------------------------
    mainflag = mogpbd(iob,nbgflg)
    ! height differ. betw. surface pressure obs and lowest model level too large
    IF (mogphd(iob,nhqcfw)/ 2 == 1)                                            &
      mainflag = IBSET ( mainflag, nvfzbp+nvfbps(4) )
    iodrbufp(iposi+ncberf) = mogpbd(iob,nbgerr)
    mqcflag  = mogpbd(iob,nbgqcf)
!   ! IWV first guess check flag moved from bit position nvru=0 to nvriwv
!   mqcflag  = ireplace1( mqcflag, nvriwv, mqcflag, nvru )
!   mqcflag  = IBCLR    ( mqcflag, nvru  )
    ! surface pressure first guess check flag added to f.g. check flag word
    iodrbufp(iposi+ncbqcf) = ireplace1( mqcflag, nvrz, mogphd(iob,nhqcfw), 0 )
    iodrbufp(iposi+ncbmfw) = mainflag
    iodrbufp(iposi+ncblsg) = IBSET ( 0, LS_SURFACE )
    iposi  =  iposi + mxcbig

!   rodrbufp(iposr+ncbu  ) = rmdi
!   rodrbufp(iposr+ncbv  ) = rmdi
    rodrbufp(iposr+ncbt  ) = ogpbdy(iob,nbgt  )
    rodrbufp(iposr+ncbrh ) = ogpbdy(iob,nbgrh )
    rodrbufp(iposr+ncbp  ) = ogpbdy(iob,nbgp  )
    rodrbufp(iposr+ncbz  ) = ogpbdy(iob,nbgz  )
!   rodrbufp(iposr+ncbz  ) = rmdi
!   rodrbufp(iposr+ncbuer) = rmdi
    rodrbufp(iposr+ncbter) = rmdi
    rodrbufp(iposr+ncbqer) = rmdi
    rodrbufp(iposr+ncbzer) = rmdi
    rodrbufp(iposr+ncbdrh) = rmdi
    rodrbufp(iposr+ncbiwa) = ogpbdy(iob,nbgiwa)
    rodrbufp(iposr+ncbiwc) = ogpbdy(iob,nbgbia)
    rodrbufp(iposr+ncbzwd) = ogpbdy(iob,nbgzwd)
    rodrbufp(iposr+ncbzpd) = ogpbdy(iob,nbgzpd)
    rodrbufp(iposr+ncbzpe) = ogpbdy(iob,nbgtze)
    iposr  =  iposr + mxcbrg

!   model-derived simulated observations for ground-based GPS reports
!   -----------------------------------------------------------------
    IF (lvofoi) THEN
!     rodrbufp(iposr+ncmu  ) = rmdi
!     rodrbufp(iposr+ncmv  ) = rmdi
      rodrbufp(iposr+ncmt  ) = rmdi
      rodrbufp(iposr+ncmrh ) = rmdi
      rodrbufp(iposr+ncmp  ) = rmdi
      IF ((ogpbdy(iob,nbgiwa) > rmdich) .AND. (sgpbdy(nso_iq,iob) > rmdich))    &
        rodrbufp(iposr+ncmiq ) = sgpbdy(nso_iq,iob)
!       rodrbufp(iposr+ncmiq ) = ogpbdy(iob,nbgiwa) - dgpbdy(nso_iq,iob)
      IF ((ogpbdy(iob,nbgzpd) > rmdich) .AND. (sgpbdy(nso_zpd,iob) > rmdich))   &
        rodrbufp(iposr+ncmzpd) = sgpbdy(nso_zpd,iob)
      iposr  =  iposr + mxcmrg
    ENDIF

  ENDIF

! set the flag for having the printing of the current report done
! ---------------------------------------------------------------

  IF (mogphd(iob,nhqcfw) >= 0)  mogphd (iob,nhqcfw) = - 99

ENDDO gps_report
!~~~~~~~~~~~~~~~

  mboxold = MAX( mboxnow , mboxold )

!-------------------------------------------------------------------------------
! End of module procedure obs_neff_fill_buffer
!-------------------------------------------------------------------------------
END SUBROUTINE obs_neff_fill_buffer


!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_cdfout_feedobs" for sorting report list
!-------------------------------------------------------------------------------

SUBROUTINE obs_neff_sort_buffer ( nxrep, ialeni, iodrbufa, ie_tot, je_tot      &
                                , nrep , ioffsrt )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "src_obs_cdfout_feedobs" makes a sorted list
!   of the offsets of the reports in the 1-dim. buffer arrays. The sorting must
!   be independent from the domain decomposition in order to obtain this
!   independency (reproducibility) also for the subsequent writing onto the
!   NetCDF feedobs file.
!
! Method:
!   It is assumed, that the order is unique for observation reports assigned
!   to the same model grid point and having the same observation time,
!   except for multi-level aircraft reports (which can be (re-)produced
!   several times).
!   In order to get a unique (i.e. independent from the domain decomposition
!   and hence fully reproducible) order, the reports are first sorted
!   according to observation time, then longitudinal grid point, then
!   latitudinal grid point, then multi-level aircraft report no/yes.
!   This is done in one step with combined QuickSort / insertion sort
!   routine (using variable 'icriterion').
!   In a second step, multi-level aircraft reports are sorted according to
!   report length and then content by means of bubble-sort.
!   Finally, the sorted offsets are obtained from the sorted reports and
!   corresponding offsets.
!   03.08.10
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
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

  INTEGER (KIND=iintegers) , INTENT (IN)    :: &
    ie_tot             ,& ! number of grid points in zonal direction
    je_tot             ,& ! number of grid points in meridional direction
    nxrep              ,& ! number of reports in global buffer
    ialeni                ! number of elements in global buffer (dim = ialenx+1)

  INTEGER (KIND=iintegers) , INTENT (INOUT) :: &
    iodrbufa (ialeni+1)   ! buffer containing integer ODR data of a all PEs
                          !   (element 'nchoff' is set in this routine)

  INTEGER (KIND=iintegers) , INTENT (OUT)   :: &
    nrep                  ! number of sorted reports

  INTEGER (KIND=iintegers) , INTENT (OUT)   :: &
    ioffsrt  (nxrep+1)    ! sorted list of offsets of the reports in iodrbufa

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND = iintegers) ::  &
    irep               ,& ! (loop) index for reports
    kl                 ,& ! (loop) index for levels (or variables)
    iposi              ,& ! address  within a buffer
    itmp               ,& ! temporary report index
    ioffset            ,& ! offset of current report in integer buffer array
    irleni   ,irlenr   ,& ! length of current report (in integer / real array)
    irleni1  ,irlenr1  ,& ! length of next    report (in integer / real array)
    ikl      ,iklp1       ! index of current element in current / next report

  LOGICAL                   :: &
    lchange               ! .TRUE. if the order of the report list is changed

  CHARACTER (LEN=20)       :: yroutine       ! name of this subroutine
! CHARACTER (LEN=40)       :: yerrmsg        ! error message

! Local arrays:
! ------------

  INTEGER (KIND = iintegers) ::  &
    isortrv  (nxrep+1) ,& ! sorted list of reports to be printed (dim = nxrep+1)
    ioffs    (nxrep+1) ,& ! list of offsets of the reports in iodrbufa
    icriterion (nxrep+1)  ! criterion accord. to which the index list is sorted
!
!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine obs_neff_sort_buffer
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 0: Initialisations
!-------------------------------------------------------------------------------

  yroutine   = 'obs_neff_sort_buffer'

! IF (my_cart_id == ircv) THEN

! Make an unsorted list of (the initial positions minus 1, i.e. the offsets of)
! the reports in '?odrbufa',
! determine the indices at which the reports start in the real elements array,
! and determine the number of reports
    nrep    = 0
    iposi   = 0
    ioffset = 0
    Report_position: DO irep = 1 , nxrep + 1
      IF (iodrbufa(iposi+nchlen) > 0) THEN
        nrep = nrep + 1
        isortrv  (nrep)         = nrep
        iodrbufa (iposi+nchoff) = ioffset
        ioffset                 = ioffset + iodrbufa(iposi+nchler)
        ioffs    (nrep)         = iposi
        iposi                   = iposi   + iodrbufa(iposi+nchlen)
      ELSE
        EXIT Report_position
      ENDIF
    ENDDO Report_position

    IF (nrep /= nxrep) PRINT '(''WARNING ",A ,": nrep =='',I5                  &
                             &,'' /= nxrep =='',I5)' ,  yroutine, nrep, nxrep
    nrep  =  MIN( nrep , nxrep )

! sort index list
    IF (nrep >= 2) THEN
      DO irep = 1 , nrep
! sorting criterion: first observation time (min), then longitude, then latitude
!                    then multi-level aircraft report yes/no
        icriterion(irep) = 2*(   iodrbufa(ioffs(irep)+nchtim) *je_tot* ie_tot  &
                               + iodrbufa(ioffs(irep)+nchi_t) *je_tot          &
                               + iodrbufa(ioffs(irep)+nchj_t) )
        IF (      (iodrbufa(ioffs(irep) + nchobt) == nairep)                   &
            .AND. (iodrbufa(ioffs(irep) + nchnlv) >= 2))                       &
          icriterion(irep) = icriterion(irep) + 1
      ENDDO

! QuickSort reports according to criterion
!   (at the beginning of 'sortrx', it is set: isortrv(i)=i ;
!    therefore, 'isortrv' rather than 'ioffsrt' has to be sorted in 'sortrx')

      CALL sortrx ( nrep, icriterion, isortrv )
!     ===========

! bubble-sort co-located multi-level aircraft reports
! according to report length and then to integer variables
      lchange = .TRUE.
      DO WHILE (lchange)
        lchange = .FALSE.
        DO irep = 1 , nrep-1
! find co-located multi-level aircraft reports
          IF (      (icriterion(isortrv(irep)) == icriterion(isortrv(irep+1))) &
              .AND. (iodrbufa(ioffs(isortrv(irep  )) + nchobt) == nairep)      &
              .AND. (iodrbufa(ioffs(isortrv(irep+1)) + nchobt) == nairep)      &
              .AND. (iodrbufa(ioffs(isortrv(irep  )) + nchnlv) >= 2)           &
              .AND. (iodrbufa(ioffs(isortrv(irep+1)) + nchnlv) >= 2)) THEN
! sort according to report length
            irleni  = iodrbufa (ioffs(isortrv(irep  ))+ nchlen)
            irleni1 = iodrbufa (ioffs(isortrv(irep+1))+ nchlen)
            irlenr  = iodrbufa (ioffs(isortrv(irep  ))+ nchler)
            irlenr1 = iodrbufa (ioffs(isortrv(irep+1))+ nchler)
            IF (     ( irleni <  irleni1 )                                     &
                .OR. ((irleni == irleni1) .AND. (irlenr <  irlenr1))) THEN
              itmp             = isortrv (irep+1)
              isortrv (irep+1) = isortrv (irep)
              isortrv (irep)   = itmp
              lchange = .TRUE.
            ELSEIF ((irleni == irleni1) .AND. (irlenr == irlenr1)) THEN
! for equal report length, sort according to succeeding variables
              Sort_variables: DO kl = 1 , irleni + irlenr
                IF (kl <= irleni)  ikl   = ioffs(isortrv(irep  )) + kl
                IF (kl <= irleni)  iklp1 = ioffs(isortrv(irep+1)) + kl
                IF (kl >  irleni)  ikl   = ioffs(isortrv(irep  )) + kl - irleni
                IF (kl >  irleni)  iklp1 = ioffs(isortrv(irep+1)) + kl - irleni
                IF (iodrbufa(ikl) < iodrbufa(iklp1)) THEN
                  itmp             = isortrv (irep+1)
                  isortrv (irep+1) = isortrv (irep)
                  isortrv (irep)   = itmp
                  lchange = .TRUE.
                                                             EXIT Sort_variables
                ENDIF
                IF (iodrbufa(ikl) > iodrbufa(iklp1))         EXIT Sort_variables
              ENDDO Sort_variables
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDIF

! convert sorted list of report indices into sorted list of offsets
    ioffsrt (nxrep+1) = 0
    IF (nrep < nxrep)  ioffsrt (nrep+1:nxrep) = 0
    DO irep = 1 , nrep
      ioffsrt (irep)  =  ioffs(isortrv(irep))
    ENDDO
 
! ENDIF  !  (my_cart_id == ircv)

!-------------------------------------------------------------------------------
! End of module procedure obs_neff_sort_buffer
!-------------------------------------------------------------------------------
END SUBROUTINE obs_neff_sort_buffer


!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_cdfout_feedobs" for writ. to NetCDF feedobs file
!-------------------------------------------------------------------------------

SUBROUTINE obs_write_cdf_feedobs ( n_rep, ialeni, iodrbufa, ialenr, rodrbufa   &
                                 , ioffsrt, lastpr, lrefend                    &
                                 , ydate_ini, hversta, hverend                 &
                                 , mruntyp, max_rep, max_body                  &
                                 , ie_tot, je_tot, ke_tot, pollat, pollon      &
                                 , dlat, dlon, startlat_tot, startlon_tot      &
                                 , iveri_ens_member, nvers                     &
                                 , yglatt_institution, yglatt_source, yfdbkdir &
                                 , nexce_rep, nexce_bdy, jerr, yerr )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "src_obs_cdfout_feedobs" writes out the
!   observations including the result of the threshold quality control,
!   as well as the values obtained by applying the observation operators
!   to the model values, to a separate file in ASCII, when the observation
!   time coincides with the current model time.
!
! Method:
!   If the program is run in parallel mode, the data to be printed are
!   gathered from all PEs at the root node.
!   Initial version by Marek Lazanowicz, 31.05.09
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
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

  INTEGER (KIND=iintegers) , INTENT (IN)  :: &
    n_rep               ,& ! number of reports to be printed now
    ialeni              ,& ! number of elements in integer buffer 'iodrbufa'
    ialenr              ,& ! number of elements in real    buffer 'rodrbufa'
    iodrbufa (ialeni+1) ,& ! buffer containing integer ODR data from all PEs
    ioffsrt  (n_rep +1)    ! sorted list of offsets of the reports in iodrbufa

  REAL    (KIND=wp)        , INTENT (IN)  :: &
    rodrbufa (ialenr+1)    ! buffer containing real ODR data from all PEs

  LOGICAL                  , INTENT (IN)  :: &
    lastpr      ,& ! last time, sth is written to ncfeedobs in this model run
    lrefend        ! .t.: verification reference time set to end of model run 
                   !      (for NUDGE runs only; for FORECAST (LETKF) or
                   !       NUDGECAST, veri_ref_time = start of model run)

  INTEGER (KIND=iintegers) , INTENT (IN)  :: &
    max_rep     ,& ! max. number of reports in NetCDF feedobs file
    max_body    ,& ! max. number of obs in NetCDF feedobs file
    mruntyp        ! =-1: print only observations ; =2: also print increments
                   ! (only if >0: arrays 'smlbdy' etc. have to be allocated and
                   !  filled previously and are used and de-allocated here)

  CHARACTER (LEN=14)       , INTENT (IN)  :: &
    ydate_ini      ! reference date (e.g. start of the model run)
                   !                [yyyymmddhhmmss]
                   !                (year, month, day, hour, min., sec.)

  REAL    (KIND=wp)        , INTENT (IN)  :: &
    hversta     ,& ! start of verification period in 'model integration hours'
    hverend        ! end of verification period in 'model integration hours'

  ! following variables are only used for supplementary info in nc file header

  INTEGER (KIND=iintegers) , INTENT (IN)  :: &
    ie_tot      ,& ! number of grid points in zonal direction
    je_tot      ,& ! number of grid points in meridional direction
    ke_tot         ! number of grid points in vertical direction

  REAL    (KIND=wp)        , INTENT (IN)  :: &
    pollat      ,& ! latitude  of the rotated north pole (in degrees, N>0)
    pollon      ,& ! longitude of the rotated north pole (in degrees, E>0)
    dlat        ,& ! grid point distance in meridional direction (in degrees)
    dlon        ,& ! grid point distance in zonal      direction (in degrees)
    startlat_tot,& ! rotated latitude   \  of the lower left grid point of the
    startlon_tot   ! rotated longitude  /  total domain (in degrees, N,E>0)

  INTEGER (KIND=iintegers) , INTENT (IN)  :: &
    iveri_ens_member  ,& ! ensemble member ( -1: deterministic)
    nvers                ! exp. id + 16384*(class of model run)

  CHARACTER (LEN= * )      , INTENT (IN)  :: &
    ! NetCDF global attributes:
    yglatt_institution,& ! originating center name
    yglatt_source        ! program name and version
!   yglatt_title      ,& ! title string for the output
!   yglatt_references    ! URL, report etc.

  CHARACTER (LEN= * )      , INTENT (IN)  :: &
    yfdbkdir             ! directory name for feedobs file

  INTEGER (KIND=iintegers) , INTENT (INOUT)  :: &
    nexce_rep   ,& ! number of reports exceeding FOF size
    nexce_bdy   ,& ! number of obs     exceeding FOF size
    jerr           ! error status variable

    ! it is assumed here that LEN >= 72 !
  CHARACTER (LEN= * )      , INTENT (INOUT)  :: &
    yerr           ! error message

! Local parameters:
! ----------------

  INTEGER (KIND = iintegers) , PARAMETER  :: &
    itype_calendar = 0    ! for specifying the calendar used in 'get_utc_date'

  REAL    (KIND = wp)        , PARAMETER  :: &
    c1000r  =     0.001_wp   ,& !
    c60     =    60.0_wp        !

! Type declarations:
! -----------------

  ! t_entry: like type 't_acc_body', but containing arrays instead of scalars
  TYPE t_entry
    INTEGER                 :: len_body               ! length of report body for
                                                      !   in netcdf feedobs file
    INTEGER       , POINTER :: varno    (:) => NULL() ! type of observed quantity
    REAL(KIND=sp) , POINTER :: obs      (:) => NULL() ! bias corrected observation
    REAL(KIND=sp) , POINTER :: bcor     (:) => NULL() ! bias correction
    REAL(KIND=sp) , POINTER :: e_o      (:) => NULL() ! observation error
    REAL(KIND=sp) , POINTER :: level    (:) => NULL() ! level of observation
    REAL(KIND=sp) , POINTER :: plevel   (:) => NULL() ! pressure level of obs.
    REAL(KIND=sp) , POINTER :: accuracy (:) => NULL() ! accuracy of observation
    REAL(KIND=sp) , POINTER :: dlat     (:) => NULL() ! latitude shift
    REAL(KIND=sp) , POINTER :: dlon     (:) => NULL() ! longitude shift
    INTEGER       , POINTER :: level_typ(:) => NULL() ! type of level information
    INTEGER       , POINTER :: level_sig(:) => NULL() ! level significance
    INTEGER       , POINTER :: state    (:) => NULL() ! status of observation
    INTEGER       , POINTER :: flags    (:) => NULL() ! observation quality check flags
    INTEGER       , POINTER :: check    (:) => NULL() ! check which which raised the
                                                      !   observation status flag value
    INTEGER       , POINTER :: qual     (:) => NULL() ! observation confidence
    REAL(KIND=sp) , POINTER :: veri_data(:,:) => NULL() ! model equivalent
  END TYPE t_entry

  TYPE (t_entry)  , POINTER  :: nc_body_buf(:) => NULL()  ! buffer for body part
                                                          ! feedobs file

  TYPE (t_fdbk)   , SAVE        :: fb         ! feedback meta file      
  TYPE (t_account), ALLOCATABLE :: report(:)  ! reports written to feedback file

! Local scalars:
! -------------

  LOGICAL                    ::  &
    lwrobs         ,& ! .TRUE., if observation type is to be written to FOF
    laddobs           ! .TRUE., if observation      is to be written to FOF

  INTEGER (KIND = iintegers) ::  &
    irep   , icl   ,& ! loop indices (reports, characters)
    iobs   , ilev  ,& ! loop indices (over observations / levels)
    kk             ,& ! counter or loop index over observations for feedobs file
    jv             ,& ! loop index over verification runs for feedobs file
    nfcast         ,& ! number of forecast runs to compare with observations
    nrefmin        ,& ! diff. betw. verif. ref. time and model initial time
    nzday          ,& ! Julian day
    ilenpath       ,& ! length of path of feedobs file
    kerr              ! error status variable

  INTEGER (KIND = iintegers) ::  &
    mxchi (3)      ,& ! integer header length for current report type
    mxcbi (3)      ,& ! integer body   length for current report type
    mxcbr (3)      ,& ! real    body   length for current report type
    mxcmr (3)      ,& ! veri-data      length for current report type
    iposi          ,& ! index offset for report in 'iodrbufa'
    iposr          ,& ! index offset for report in 'rodrbufa'
    iposi_b        ,& ! index offset for integer body
    iposr_b        ,& ! index offset for real body
    iposr_i        ,& ! index offset for increment
    iposi_bdy      ,& ! index offset for level in integer body
    iposr_bdy      ,& ! index offset for level in real body
    iposr_ver      ,& ! index offset for level in increment section
    iposi_obs      ,& ! index address (position) for elements in integer body
    iposr_obs      ,& ! index address (position) for elements in real body
    iposr_inc      ,& ! index address (position) for increment elements
    len_body       ,& ! number of body elements in report used in feedback file
    nw_rep         ,& ! number of reports to be put into FOF
    n_body         ,& ! number of bodys to be put into FOF
    ileni_h        ,& ! 'iodrbufa': length of integer report header
    ileni_b1       ,& ! 'iodrbufa': length of one level in integer body
    ilenr_b1       ,& ! 'iodrbufa': length of one level in real body 
    ilenr_i1       ,& ! 'iodrbufa': length of one level report incements
    kfmt           ,& ! 'iodrbufa': type of report format
    mflgqcf        ,& ! 'iodrbufa': threshold quality control flag word
    mflgerf        ,& ! 'iodrbufa': status flag word
    mflgmfw        ,& ! 'iodrbufa': main flag word
    nvfxbp         ,& ! 'iodrbufa': index for main flag
    nvrx           ,& ! 'iodrbufa': index for status / QC flags
    nvrx2          ,& ! 'iodrbufa': index for temporary QC flags (LBC checks)
!   kschr          ,& ! 'iodrbufa': station characteristics
    kphase         ,& ! 'iodrbufa': flight phase
    ksig           ,& ! 'iodrbufa': level significance
    kobv           ,& ! 'iodrbufa': observation value
    kflg           ,& ! 'nc_body_buf': (quality control) flag word
    kfbit          ,& ! 'nc_body_buf': one (quality control) flag
    kcl            ,& ! 'nc_body_buf': cloud level
    istatus           ! 'nc_body_buf': status of observation

  REAL    (KIND = wp)        ::  &
    zhr               ! hour

  INTEGER                    ::  &
    iveri_exp_id       ,& ! experiment identity number of the model verif. run
    iveri_run_class    ,& ! class (main / ass.) of the model verification run
    iveri_run_type     ,& ! run type (fcst./f.g./ana.) of the verification run
    iveri_ref_date     ,& ! reference date
    iveri_ref_time     ,& ! reference time
    iveri_start        ,& ! start of verification rel. to reference time
    iveri_end          ,& ! end   of verification rel. to reference time
    iveri_forecast_time,& ! forecast time (hhhmm)
    iiveri_ens_member  ,& ! ensemble member ( -1: deterministic)
    idomain (3)        ,& ! model domain size, input for 'create_fdbk'
    imax_rep           ,& ! max. # reports in NetCDF feedobs file, type integer
    imax_body          ,& ! max. # obs in NetCDF feedobs file, type integer
    varid              ,& ! netcdf variable identifier, for 'add_verification'
    inw_rep            ,& ! number of reports to be put into FOF
    in_body            ,& ! number of bodys to be put into FOF
    ibdy_pos           ,& ! time box initial position in body of FOF
    ihdr_pos           ,& ! time box initial position in header of FOF
    imiss                 ! missing data indicator for integers (2^31-1)

  REAL(KIND=sp)         &
    resolution  (2)    ,& ! model resolution (lat, lon)
    pole        (2)    ,& ! pole (lat, lon) of rotated grid
    lower_left  (2)    ,& ! lower left  corner of model domain (lat, lon)
    upper_right (2)    ,& ! upper right corner of model domain (lat, lon)
    rmich                 ! check value for missing real data   (-1.E30)

  CHARACTER (LEN=12)  :: yveri_ref_datim      ! yyyymmddhhmm reference time
  CHARACTER (LEN=12)  :: yveri_initial_date   ! yyyymmddhhmm initial time
  CHARACTER (LEN=45)  :: yveri_descript       ! description of model verif. run
  CHARACTER (LEN=20)  :: yens                 ! description of ensemble member
  CHARACTER (LEN=64)  :: yversion             ! model version number
! CHARACTER (LEN=PLEN-10)     :: yfdbkdir
  CHARACTER (LEN=PLEN) , SAVE :: yfdbkfile
  CHARACTER (LEN=14)  :: yakdat1  ! actual date in the form   yyyymmddhhmmss
  CHARACTER (LEN=28)  :: yakdat2  ! date, form:   wd   dd.mm.yy  hh mm ss UTC

  INTEGER (KIND = iintegers) , SAVE   :: &
    n_tot_verrep  = 0   ,& ! total number of reports
    n_tot_verbody = 0   ,& ! total number of body  reports
    n_veri              ,& ! number of verifications, always 1 in COSMO model
    nfl_                ,& ! number of flags in 'flags' or 'r_flags' entries
    kfl_ (64)     = 0      ! number of flags in 'flags' or 'r_flags' entries

  LOGICAL                    , SAVE   :: &
    lfdbk_crea = .TRUE. ,& ! .TRUE. if feedback file needs to be created
                           !        (i.e. when first obs has to be written)
    lfdbk_open = .FALSE.   ! .TRUE.  if feedback file is open

! Local arrays:
! ------------

  INTEGER (KIND = iintegers) , ALLOCATABLE  :: &
    n_levl    (:) ,& ! number of levels in report
    koffset   (:)    ! buffer containing body offset for feedback file
!
!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine obs_write_cdf_feedobs
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 0: Initialisations
!-------------------------------------------------------------------------------

  jerr       = 0
  yerr       = '  '
  kerr       = 0
  zhr        = REAL (INT( hversta + epsy ), wp)
  nfcast     = MAX( 0 , MIN( 1 , mruntyp+1 ) )

  ! difference between model initial time and verification reference time
  IF (      lrefend)  nrefmin  =  NINT( hverend *c60 )
  IF (.NOT. lrefend)  nrefmin  =  0

! IF ((my_cart_id == ircv) .AND. (n_rep > 0)) THEN

!-------------------------------------------------------------------------------
!  Section 1:  Create the netCDF feedback file 
!-------------------------------------------------------------------------------

  IF ( (lfdbk_crea) .AND. (n_rep > 0) ) THEN

! determine NetCDF global attributes
! ----------------------------------

    CALL get_utc_date ( nrefmin, ydate_ini, c60, itype_calendar, yakdat1       &
                      , yakdat2, nzday, zhr)
!   =================

    READ( yakdat1(1:8) ,'(I8)' )  iveri_ref_date
    READ( yakdat1(9:12),'(I4)' )  iveri_ref_time
    WRITE( yveri_ref_datim,'(I8.8,I4.4)' )  iveri_ref_date, iveri_ref_time
    iveri_end    =  NINT( hverend *c60 )  -  nrefmin
    iveri_start  =  NINT( hversta *c60 )  -  nrefmin
    IF (hversta > epsy)  iveri_start  =  MAX( iveri_start , 1-nrefmin )

    resolution(1)  = REAL (dlat)
    resolution(2)  = REAL (dlon)
    idomain(1)     = ie_tot
    idomain(2)     = je_tot
    idomain(3)     = ke_tot
    pole(1)        = REAL (pollat)
    pole(2)        = REAL (pollon)
    lower_left(1)  = REAL (startlat_tot)                     ! / lower left corner
    lower_left(2)  = REAL (startlon_tot)                     ! \ of model domain
    upper_right(1) = REAL (startlat_tot + (je_tot-1) * dlat) ! / upper right corner
    upper_right(2) = REAL (startlon_tot + (ie_tot-1) * dlon) ! \ of model domain
    yversion = yglatt_source
    IF (yversion(1:5) == 'COSMO') yversion = yversion(7:LEN_TRIM(yversion)-6)
    imax_rep       = max_rep
    imax_body      = max_body

! determine NetCDF variables related to the verification data
! (i.e. to the current model run)
! -----------------------------------------------------------

    yveri_descript    =   '                            '
    IF (nfcast >= 1) THEN
      IF ((nvers >= 0) .AND. (nvers < 1048576)) THEN
        iveri_exp_id     =  MOD( nvers, 16384 )
        iveri_run_class  =       nvers/ 16384
      ELSE
        iveri_exp_id     =  imdi
        iveri_run_class  =  imdi
      ENDIF
      yveri_initial_date  =   yveri_ref_datim

      IF (lrefend) THEN
        !    in this case, it is assumed that (mruntyp == 2) is set correctly
        iveri_forecast_time =   0_iintegers
        yveri_descript      =   'nudging assimilation        '
        iveri_run_type      =   VT_ANALYSIS
      ELSE
        !    in this case, it is assumed that (mruntyp /= 2)
        iveri_forecast_time =   INT( hverend +c1000r ) * 100                   &
                              + MOD( NINT( hverend *c60 ) , 60 )
        IF (iveri_ens_member >= 1)                                             &
          WRITE( yens, '(" ensemble member ",I3)' )  iveri_ens_member
        IF (mruntyp == 1) THEN
          iveri_run_type    =   VT_FIRSTGUESS
          yveri_descript    =   'first guess                 '
          IF (iveri_ens_member >= 1)  yveri_descript  =  'first guess' // yens
        ELSE
          iveri_run_type    =   VT_FORECAST
          yveri_descript    =   'deterministic forecast      '
          IF (iveri_ens_member >= 1)  yveri_descript  =  'forecast' // yens
        ENDIF
      ENDIF
      iiveri_ens_member = iveri_ens_member
    ENDIF

! define the file path name
! -------------------------

    PRINT    '("Creation of NetCDF feedobs file fof_*")'
!   yfdbkfile = 'test_COSMOfeedback.nc'
!CS: wish of Hanisch: do not use ensemble member number as suffix to file name
!   IF ((iveri_ens_member >= 0) .AND. (iveri_ens_member <= 999)) THEN
!     WRITE( yfdbkfile,'("fof_",A12,"00.",I3.3,".nc")' )  yveri_ref_datim      &
!                                                       , iveri_ens_member
!   ELSE
      WRITE( yfdbkfile,'("fof_",A12,"00"      ,".nc")' )  yveri_ref_datim
!   ENDIF
!   yfdbkdir  = ''
    ilenpath   = LEN_TRIM(yfdbkdir)
    IF (ilenpath > PLEN-25-1) THEN
      WRITE( yerr,'("directory name for feedobs/ grib output too long: >",I3)' &
           )  PLEN-26
      jerr = ilenpath
                                                                          RETURN
    ENDIF
    IF (ilenpath >= 1) THEN
      IF (yfdbkdir(ilenpath:ilenpath) == '/')  ilenpath = ilenpath - 1
      yfdbkfile= yfdbkdir(1:ilenpath) // '/' // yfdbkfile(1:LEN_TRIM(yfdbkfile))
    ENDIF

! create the feedback file (and initialise fdbk tables)
! -------------------------

    CALL setup_fdbk   (fb% nc, latex = .false.)
!   ===============

    nfl_  =  flags% n 
    DO kk = 1 , nfl_
      kfl_ (kk)  =  flags% e(kk)% value
    ENDDO

!CS: 'opt' still needs to be complemented

    CALL create_fdbk (fb                      ,&! feedback file meta data
                      yfdbkfile               ,&! path
                      'COSMO'                 ,&! model
                      yversion                ,&! model version
                      yglatt_institution      ,&! institution
                      imax_rep                ,&! d_hdr
                      imax_body               ,&! d_body
                      iveri_ref_date          ,&! reference date
                      iveri_ref_time          ,&! reference time
                      iveri_start             ,&! start of verification
                      iveri_end               ,&! end   of verification
                      resolution              ,&! resolution
                      idomain                 ,&! domain size
                      yveri_descript          ,&! comment (for history)
                      yveri_ref_datim         ,&! time    (for history)
                      pole       =pole        ,&! pole of rotated grid
                      lower_left =lower_left  ,&! l.l. corner of domain
                      upper_right=upper_right ,&! u.r. corner of domain
                      opt='COSMO TEMP PILOT AIREP') ! flag optional variables
!   ================

    IF (nfcast >= 1)                                                           &

      CALL add_verification (fb                 ,&! feedback file meta data
                            'COSMO'             ,&! model 
                            iveri_run_type      ,&! run type
                            iveri_run_class     ,&! run class
                            yveri_initial_date  ,&! initial date 
                            iveri_forecast_time ,&! forecast time (hhhmm)
                            resolution          ,&! resolution
                            idomain             ,&! domain size
                            yveri_descript      ,&! description
                            iiveri_ens_member   ,&! ensemble member 
                            iveri_exp_id        ,&! experiment id
                            varid )               ! variable id
!     ====================

    n_veri          = nfcast  
    lfdbk_open      = .TRUE.
    lfdbk_crea      = .FALSE.
    PRINT    '("OPENED: ",A)' ,  yfdbkfile(1:MIN(72, LEN_TRIM(yfdbkfile)))
  ENDIF   !   (lfdbk_crea) .and (nrep > 0)

! DO irep = 1, n_rep
!   iposi = ioffsrt(irep)
!   iposr = iodrbufa(iposi+nchoff)
!   IF (irep < 10)                                                             &
!     WRITE(0,'("LAT/LON",5I8,2F8.2)') iposi, iposr, iodrbufa(iposi+nchide)    &
!            , iodrbufa(iposi + nchcdt), iodrbufa(iposi+ nchpas)               &
!            , rodrbufa(iposr + nchlat), rodrbufa(iposr+ nchlon) 
! ENDDO

!-------------------------------------------------------------------------------
!  Section 2:  Fill 'nc_body_buf' which contains all the values of the body
!              section of the NetCDF Feedobs File FOF in the correct form;
!              in section 3, 'nc_body_buf' is copied onto 'report% body'
!-------------------------------------------------------------------------------

  IF (n_rep > 0) THEN
    IF (jerr == 0)  ALLOCATE ( n_levl      (n_rep) , STAT = jerr )
    IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 1    
    IF (jerr == 0)  ALLOCATE ( koffset     (n_rep) , STAT = jerr )
    IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 2    
    IF (jerr == 0)  ALLOCATE ( nc_body_buf (n_rep) , STAT = jerr )
    IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 3    
    IF (jerr /= 0) THEN
      IF (kerr==1) WRITE( yerr,'("ERROR in alloc: n_levl (n_rep)  ",I8)' ) n_rep
      IF (kerr==2) WRITE( yerr,'("ERROR in alloc: koffset (n_rep) ",I8)' ) n_rep
      IF (kerr==3) WRITE( yerr,'("ERROR in alloc:nc_body_buf(n_rep)",I8)') n_rep
                                                                          RETURN
    ENDIF
    n_levl  (:) = -999
    koffset (:) =  0
  ENDIF   !   (nrep > 0)

! get number of observation levels for each report
! ------------------------------------------------
  DO irep = 1 , n_rep
    iposi = ioffsrt(irep)
    IF ((iodrbufa(iposi+ nchfmt) == 2) .OR. (iodrbufa(iposi+ nchfmt) == 3)) THEN
    ! surface-level or GPS report
      n_levl (irep) = 1
    ELSE
    ! upper-air (single-level or multi-level) report
      n_levl (irep) = iodrbufa(iposi + nchnlv)
    ENDIF
  ENDDO

! number of elements for upper-air / surface / GPS report in 'xodrbufa'
! ---------------------------------------------------------------------
  mxchi (1) = mxchiu
  mxchi (2) = mxchis
  mxchi (3) = mxchis
!             mxchra
  mxcbi (1) = mxcbiu
  mxcbi (2) = mxcbis
  mxcbi (3) = mxcbig
  mxcbr (1) = mxcbru
  mxcbr (2) = mxcbrs
  mxcbr (3) = mxcbrg
  mxcmr (1) = mxcmru
  mxcmr (2) = mxcmrs
  mxcmr (3) = mxcmrg

! determine body length for each report
! -------------------------------------
  DO irep = 1, n_rep
    iposi = ioffsrt(irep)
    kfmt  = iodrbufa(iposi + nchfmt)
    ilenr_b1 = mxcbr(kfmt)
    iposr_b  = iodrbufa(iposi + nchoff) + mxchra                  
    len_body = 0
    DO  ilev = 1, n_levl(irep)
      DO iobs = 1, ilenr_b1
        iposr_obs = iposr_b + (ilev-1)* ilenr_b1 + iobs
        IF ( rodrbufa(iposr_obs) > rmdich ) THEN
          IF (kfmt == 1) THEN
            IF (ixobs_bru(iobs) >= -1)  len_body = len_body + 1
          ELSEIF (kfmt == 2) THEN
            IF (ixobs_brs(iobs) >= -1)  len_body = len_body + 1
          ELSEIF (kfmt == 3) THEN
            IF (ixobs_brg(iobs) >= -1)  len_body = len_body + 1
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    ! surface-level reports: integer obs (WW, cloud layers)
    IF (kfmt == 2) THEN
      DO iobs = 1, mxcbis
        iposi_obs = iposi + mxchis + ilstidc + iobs
        IF ((iodrbufa(iposi_obs) /= imdi) .AND. (ixobs_bis(iobs) >= 0)) THEN
          IF (iobs == ncbwwe) THEN
            kobv  =  IBITS( iodrbufa(iposi_obs), nvw0bp, nvw0oc )
            IF (kobv /= nibits( nvw0oc )) THEN
              len_body = len_body + 1
            ENDIF
          ELSE
            len_body = len_body + 1
          ENDIF
        ENDIF
      ENDDO
    ENDIF
    nc_body_buf(irep)% len_body = len_body
  ENDDO

! initialize 'nc_body_buf'
! ------------------------
  DO irep = 1, n_rep
    len_body  =  nc_body_buf(irep)% len_body

    kerr = 0
    IF (jerr==0) ALLOCATE( nc_body_buf(irep)% varno  (len_body), STAT=jerr)
    IF ((jerr/=0) .AND. (kerr == 0)) kerr = 11
    IF (jerr==0) ALLOCATE( nc_body_buf(irep)% obs    (len_body), STAT=jerr)
    IF ((jerr/=0) .AND. (kerr == 0)) kerr = 12
    IF (jerr==0) ALLOCATE( nc_body_buf(irep)% bcor   (len_body), STAT=jerr)
    IF ((jerr/=0) .AND. (kerr == 0)) kerr = 13
    IF (jerr==0) ALLOCATE( nc_body_buf(irep)% e_o    (len_body), STAT=jerr)
    IF ((jerr/=0) .AND. (kerr == 0)) kerr = 14
    IF (jerr==0) ALLOCATE( nc_body_buf(irep)% level  (len_body), STAT=jerr)
    IF ((jerr/=0) .AND. (kerr == 0)) kerr = 15
    IF (jerr==0) ALLOCATE( nc_body_buf(irep)% plevel (len_body), STAT=jerr)
    IF ((jerr/=0) .AND. (kerr == 0)) kerr = 16
    IF (jerr==0) ALLOCATE( nc_body_buf(irep)% state  (len_body), STAT=jerr)
    IF ((jerr/=0) .AND. (kerr == 0)) kerr = 17
    IF (jerr==0) ALLOCATE( nc_body_buf(irep)% flags  (len_body), STAT=jerr)
    IF ((jerr/=0) .AND. (kerr == 0)) kerr = 18
    IF (jerr==0) ALLOCATE( nc_body_buf(irep)% check  (len_body), STAT=jerr)
    IF ((jerr/=0) .AND. (kerr == 0)) kerr = 19
    IF (jerr==0) ALLOCATE( nc_body_buf(irep)% qual   (len_body), STAT=jerr)
    IF ((jerr/=0) .AND. (kerr == 0)) kerr = 20
    IF (jerr==0) ALLOCATE( nc_body_buf(irep)%level_typ(len_body),STAT=jerr)
    IF ((jerr/=0) .AND. (kerr == 0)) kerr = 21
    IF (jerr==0) ALLOCATE( nc_body_buf(irep)%level_sig(len_body),STAT=jerr)
    IF ((jerr/=0) .AND. (kerr == 0)) kerr = 22
    IF (jerr==0) ALLOCATE( nc_body_buf(irep)% accuracy(len_body),STAT=jerr)
    IF ((jerr/=0) .AND. (kerr == 0)) kerr = 23
    IF (jerr==0) ALLOCATE( nc_body_buf(irep)% dlat   (len_body), STAT=jerr)
    IF ((jerr/=0) .AND. (kerr == 0)) kerr = 24
    IF (jerr==0) ALLOCATE( nc_body_buf(irep)% dlon   (len_body), STAT=jerr)
    IF ((jerr/=0) .AND. (kerr == 0)) kerr = 25
    IF (jerr==0) ALLOCATE( nc_body_buf(irep)%veri_data(len_body,1),STAT=jerr)
    IF ((jerr/=0) .AND. (kerr == 0)) kerr = 26
    IF (jerr/=0) THEN
      icl = len_body
      IF (kerr==11) WRITE( yerr,'("ERROR in alloc: nc_body_buf%varno ",I8)') icl
      IF (kerr==12) WRITE( yerr,'("ERROR in alloc: nc_body_buf%obs   ",I8)') icl
      IF (kerr==13) WRITE( yerr,'("ERROR in alloc: nc_body_buf%bcor  ",I8)') icl
      IF (kerr==14) WRITE( yerr,'("ERROR in alloc: nc_body_buf%e_o   ",I8)') icl
      IF (kerr==15) WRITE( yerr,'("ERROR in alloc: nc_body_buf%level ",I8)') icl
      IF (kerr==16) WRITE( yerr,'("ERROR in alloc: nc_body_buf%plevel",I8)') icl
      IF (kerr==17) WRITE( yerr,'("ERROR in alloc: nc_body_buf%state ",I8)') icl
      IF (kerr==18) WRITE( yerr,'("ERROR in alloc: nc_body_buf%flags ",I8)') icl
      IF (kerr==19) WRITE( yerr,'("ERROR in alloc: nc_body_buf%check ",I8)') icl
      IF (kerr==20) WRITE( yerr,'("ERROR in alloc: nc_body_buf%qual  ",I8)') icl
      IF (kerr==21) WRITE( yerr,'("ERROR in alloc: nc_body*%level_typ",I8)') icl
      IF (kerr==22) WRITE( yerr,'("ERROR in alloc: nc_body*%level_sig",I8)') icl
      IF (kerr==23) WRITE( yerr,'("ERROR in alloc: nc_body*%accuracy ",I8)') icl
      IF (kerr==24) WRITE( yerr,'("ERROR in alloc: nc_body_buf*dlat  ",I8)') icl
      IF (kerr==25) WRITE( yerr,'("ERROR in alloc: nc_body_buf*dlon  ",I8)') icl
      IF (kerr==26) WRITE( yerr,'("ERROR in alloc: nc_body*%veri_data",I8)') icl
                                                                          RETURN
    ENDIF
    nc_body_buf(irep)% varno     (:) = imdi
    nc_body_buf(irep)% obs       (:) = rmdi
    nc_body_buf(irep)% bcor      (:) = rmdi
    nc_body_buf(irep)% e_o       (:) = rmdi
    nc_body_buf(irep)% level     (:) = rmdi
    nc_body_buf(irep)% plevel    (:) = rmdi
    nc_body_buf(irep)% accuracy  (:) = rmdi
    nc_body_buf(irep)% dlat      (:) = rmdi
    nc_body_buf(irep)% dlon      (:) = rmdi
    nc_body_buf(irep)% level_typ (:) = imdi
    nc_body_buf(irep)% level_sig (:) = imdi
    nc_body_buf(irep)% state     (:) = imdi
    nc_body_buf(irep)% flags     (:) = imdi
    nc_body_buf(irep)% check     (:) = imdi
    nc_body_buf(irep)% qual      (:) = imdi
    nc_body_buf(irep)% veri_data (:,:) = rmdi

! determine offset, lengths etc. related to input buffer arrays 'xodrbufa'
! ------------------------------------------------------------------------
    iposi = ioffsrt(irep)
    kfmt  = iodrbufa(iposi + nchfmt)
!   kschr = iodrbufa(iposi + nchsch)

    ileni_h  =  mxchi(kfmt) + ilstidc
    ileni_b1 =  mxcbi(kfmt)
    ilenr_b1 =  mxcbr(kfmt)
    ilenr_i1 =  mxcmr(kfmt) *nfcast
!   len_body =  nc_body_buf(irep)% len_body

    ! offset for real body
    iposi_b = iposi + ileni_h
    ! offset for real body
    iposr_b = iodrbufa(iposi + nchoff) + mxchra
    ! offset for incements
    iposr_i = iodrbufa(iposi + nchoff) + mxchra + ilenr_b1 *n_levl(irep)

    kk = 0

! fill observations of type real (from 'rodrbufa' into 'nc_body_buf')
! -------------------------------------------------------------------
    DO  ilev = 1, n_levl(irep)
      iposi_bdy = iposi_b + (ilev-1)* ileni_b1
      iposr_bdy = iposr_b + (ilev-1)* ilenr_b1
      iposr_ver = iposr_i + (ilev-1)* ilenr_i1
      iposr_inc = iposr_ver
      mflgqcf   = iodrbufa(iposi_bdy + ncbqcf)
      mflgerf   = iodrbufa(iposi_bdy + ncberf)
      mflgmfw   = iodrbufa(iposi_bdy + ncbmfw)

      DO iobs = 1, ilenr_b1
        lwrobs = .FALSE.
        IF (kfmt == 1) THEN
          lwrobs  =  (ixobs_bru(iobs) >= -1)
        ELSEIF (kfmt == 2) THEN
          lwrobs  =  (ixobs_brs(iobs) >= -1)
        ELSEIF (kfmt == 3) THEN
          lwrobs  =  (ixobs_brg(iobs) >= -1)
        ENDIF
        iposr_obs = iposr_bdy + iobs          ! address for real observation
        IF ((lwrobs) .AND. (rodrbufa(iposr_obs) > rmdich)) THEN
          kk = kk + 1  
          IF (kfmt == 1) THEN
            nc_body_buf(irep)% varno     (kk) = ibu_varno(iobs)
            nc_body_buf(irep)% level     (kk) = REAL(rodrbufa(iposr_bdy + ncbp),sp)
            nc_body_buf(irep)% level_typ (kk) = VN_P
            nc_body_buf(irep)% level_sig (kk) = iodrbufa(iposi_bdy + ncblsg)
            nc_body_buf(irep)% qual      (kk) = iodrbufa(iposi_bdy + ncbsnr)
            IF (iodrbufa(iposi + nchobt) == nairep) THEN
              !   phase will be set in the header word 'phase' instead of here
!             nc_body_buf(irep)% level_sig (kk) = IBITS( kschr, nvapbp         &
!                                                      , nvapoc-1 )
              IF (      (rodrbufa(iposr_bdy + ncbz) > rmdich)                  &
                  .AND. (BTEST( mflgmfw, nvfzbp+nvfbps(5) ))) THEN
                ! aircraft reports: pressure level is converted from reported
                ! flight level using ICAO standard atmosphere
                nc_body_buf(irep)% level_typ (kk) = VN_FLEV
                nc_body_buf(irep)% level     (kk) = REAL(rodrbufa(iposr_bdy + ncbz),sp)
                nc_body_buf(irep)% plevel    (kk) = REAL(rodrbufa(iposr_bdy + ncbp),sp)
              ENDIF
            ELSEIF (BTEST( mflgmfw, nvfzbp+nvfbps(6) )) THEN
              ! pressure level is converted from reported height using model
              ! atmosphere, therefore height should be used as level
              ! because level should be independent from model (ensemble member)
              nc_body_buf(irep)% level_typ (kk) = VN_HEIGHT
              nc_body_buf(irep)% level     (kk) = REAL(rodrbufa(iposr_bdy + ncbz),sp)
              nc_body_buf(irep)% plevel    (kk) = REAL(rodrbufa(iposr_bdy + ncbp),sp)
            ENDIF
            iposr_inc  =  iposr_ver  +  ixobs_bru(iobs)
!           IF (nfcast >= 1)  iposr_inc = iposr_ver + ixobs_bru(iobs)
          ELSEIF (kfmt == 2) THEN
            nc_body_buf(irep)% varno     (kk) = ibs_varno(iobs)
!           nc_body_buf(irep)% level     (kk) = REAL(rodrbufa(iposr_bdy + ncbz),sp)
            nc_body_buf(irep)% level     (kk) = iodrbufa(iposi     + nchalt)
            nc_body_buf(irep)% level_typ (kk) = VN_HEIGHT
            nc_body_buf(irep)% level_sig (kk) = imdi
            IF (      (     (iodrbufa(iposi + nchobt) == ntemp )               &
                       .OR. (iodrbufa(iposi + nchobt) == npilot))              &
                .AND. (     (iobs == ncbu) .OR. (iobs == ncbv)                 &
                       .OR. (iobs == ncbt) .OR. (iobs == ncbrh)))              &
              nc_body_buf(irep)% level_sig (kk) = iodrbufa(iposi_bdy + ncblsg)
            IF (iobs == ncbp)                                                  &
              nc_body_buf(irep)% level_sig (kk) = iodrbufa(iposi_bdy + ncblsg)
            nc_body_buf(irep)% qual      (kk) = imdi
            iposr_inc  =  iposr_ver  +  ixobs_brs(iobs)
          ELSEIF (kfmt == 3) THEN
            nc_body_buf(irep)% varno     (kk) = ibg_varno(iobs)
!           nc_body_buf(irep)% level     (kk) = REAL(rodrbufa(iposr_bdy + ncbz),sp)
            nc_body_buf(irep)% level     (kk) = iodrbufa(iposi     + nchalt)
            nc_body_buf(irep)% level_typ (kk) = VN_HEIGHT
            nc_body_buf(irep)% level_sig (kk) = imdi
            nc_body_buf(irep)% qual      (kk) = imdi
            iposr_inc  =  iposr_ver  +  ixobs_brg(iobs)
          ENDIF
          nc_body_buf(irep)% obs       (kk) = REAL(rodrbufa(iposr_obs),sp)
          nc_body_buf(irep)% bcor      (kk) = c0
          nc_body_buf(irep)% e_o       (kk) = REAL(rmdi,sp)
          nc_body_buf(irep)% accuracy  (kk) = REAL(rmdi,sp)
          nc_body_buf(irep)% dlat      (kk) = REAL(rmdi,sp)
          nc_body_buf(irep)% dlon      (kk) = REAL(rmdi,sp)

          ! set variables that depend on observed quantity
          nvrx   = -1
          nvrx2  = -1
          nvfxbp = -1
              ! horizontal wind
          IF ((iobs == ncbu) .OR. (iobs == ncbv)) THEN
            nc_body_buf(irep)% e_o     (kk) = REAL(rodrbufa(iposr_bdy + ncbuer),sp)
            IF (kfmt == 1)                                                     &
              nc_body_buf(irep)% accuracy(kk) = REAL(rodrbufa(iposr_bdy + ncbuac),sp)
            nvfxbp = nvfubp
            nvrx   = nvru
            IF (iodrbufa(iposi + nchobt) == nairep)                            &
              nc_body_buf(irep)% qual  (kk) = IBITS( mflgmfw                   &
                                                   , nvfubp+nvfbps(1), 2 )
              ! temperature
          ELSEIF (iobs == ncbt) THEN
            nc_body_buf(irep)% e_o     (kk) = REAL(rodrbufa(iposr_bdy + ncbter),sp)
            nvfxbp = nvftbp
            nvrx   = nvrt
              ! humidity
          ELSEIF (iobs == ncbrh) THEN
            nc_body_buf(irep)% e_o     (kk) = REAL(rodrbufa(iposr_bdy + ncbqer),sp)
            nc_body_buf(irep)% bcor    (kk) = REAL(rodrbufa(iposr_bdy + ncbdrh),sp)
            nvfxbp = nvfqbp
            nvrx   = nvrq
            nvrx2  = nvrqbc
              ! pressure / height
          ELSEIF ((iobs == ncbp) .OR. (iobs == ncbz)) THEN
            IF ((kfmt == 2) .OR. (kfmt == 3))                                  &
              nc_body_buf(irep)% level (kk) = REAL(rodrbufa(iposr_bdy + ncbz),sp)
            nc_body_buf(irep)% e_o     (kk) = REAL(rodrbufa(iposr_bdy + ncbzer),sp)
            nvfxbp = nvfzbp
            nvrx   = nvrz
            nvrx2  = nvrzbc
          ENDIF
          ! variables specific to upper-air reports:
                ! vertical velocity
          IF ((kfmt == 1) .AND. (iobs == ncbw  )) THEN
            nvrx   = nvrw
          ! variables specific to ground-based GPS reports
          ELSEIF (kfmt == 3) THEN
                ! IWV
            IF (iobs == ncbiwa) THEN
              IF (rodrbufa(iposr_bdy + ncbiwc) > rmdich)                       &
                nc_body_buf(irep)% bcor  (kk) = REAL(rodrbufa(iposr_bdy + ncbiwc),sp)
            ! compute observation error for IWV !
              IF (rodrbufa(iposr_bdy + ncbzwd) > c05) THEN
                ! if ZWD exists, not close to 0 : IWV-err = ZWD-err * IWV / ZWD
                nc_body_buf(irep)% e_o   (kk) = REAL(rodrbufa(iposr_bdy + ncbzpe)   &
                                                   * rodrbufa(iposr_bdy + ncbiwa)   &
                                                   / rodrbufa(iposr_bdy + ncbzwd),sp)
              ELSE
                ! 0.17 = approx. ratio IWV[mm] / ZWD[mm]
                nc_body_buf(irep)% e_o   (kk) = REAL(rodrbufa(iposr_bdy + ncbzpe)   &
                                              * 0.17_wp,sp)
              ENDIF
              nvfxbp = nvfgbp
              nvrx   = nvriwv
              nvrx2  = nvrqbc
                ! ZWD
            ELSEIF (iobs == ncbzwd) THEN
              nc_body_buf(irep)% e_o     (kk) = REAL(rodrbufa(iposr_bdy + ncbzpe),sp)
              nvfxbp = nvfgbp
              nvrx   = nvrzpd
              nvrx2  = nvrqbc
                ! ZPD
            ELSEIF (iobs == ncbzpd) THEN
              nc_body_buf(irep)% e_o     (kk) = REAL(rodrbufa(iposr_bdy + ncbzpe),sp)
              nvfxbp = nvfgbp
              nvrx   = nvrzpd
              nvrx2  = nvrqbc
            ENDIF
          ! variables specific to surface reports
          ELSEIF (kfmt == 2) THEN
                ! total cloud
            IF (iobs == ncbct ) THEN
              nc_body_buf(irep)% level_sig (kk) = imdi
              icl  =  IBITS( iodrbufa(iposi_bdy + ncbclg), nxsgbp, nxsgoc )
              IF (icl /= 63)  nc_body_buf(irep)% level_sig (kk) = icl
              nvrx   = nvrct
                ! low cloud
            ELSEIF (iobs == ncbcl ) THEN
              nvrx   = nvrcl
              nc_body_buf(irep)% level_sig (kk) = 7
                ! mid-level cloud
            ELSEIF (iobs == ncbcm ) THEN
              nvrx   = nvrcm
              nc_body_buf(irep)% level_sig (kk) = 8
                ! high cloud
            ELSEIF (iobs == ncbch ) THEN
              nvrx   = nvrch
              nc_body_buf(irep)% level_sig (kk) = 9
                ! temporally non-local observations: pressure tendency, precip,
                !                                    gusts, min./max temperature
            ELSEIF ((iobs == ncbff) .OR. (iobs == ncbdd)) THEN
              nvfxbp = nvfubp
              nvrx   = nvru
            ELSEIF (iobs == ncbtd) THEN
              nvfxbp = nvfqbp
              nvrx   = nvrq
              nvrx2  = nvrqbc
            ELSEIF (iobs == ncbcbs) THEN
              nvrx   = nvrcbs
!           ELSEIF (iobs == ncbvis) THEN
!             nvrx2  = -9
!           ELSEIF (iobs == ncbhsw) THEN
!             nvrx2  = -9
            ELSEIF (iobs == ncbpst) THEN
              nc_body_buf(irep)% level     (kk) =  3._sp
              nc_body_buf(irep)% level_typ (kk) = VN_TRTR
            ELSEIF ((iobs == ncbrr1) .OR. (iobs == ncbfg1)) THEN
              nc_body_buf(irep)% level     (kk) =  1._sp
              nc_body_buf(irep)% level_typ (kk) = VN_TRTR
            ELSEIF ((iobs == ncbrr3) .OR. (iobs == ncbfg3)) THEN
              nc_body_buf(irep)% level     (kk) =  3._sp
              nc_body_buf(irep)% level_typ (kk) = VN_TRTR
            ELSEIF ((iobs == ncbrr6) .OR. (iobs == ncbfg6)) THEN
              nc_body_buf(irep)% level     (kk) =  6._sp
              nc_body_buf(irep)% level_typ (kk) = VN_TRTR
            ELSEIF  (iobs == ncbr12) THEN
              nc_body_buf(irep)% level     (kk) = 12._sp
              nc_body_buf(irep)% level_typ (kk) = VN_TRTR
            ELSEIF  (iobs == ncbr24) THEN
              nc_body_buf(irep)% level     (kk) = 24._sp
              nc_body_buf(irep)% level_typ (kk) = VN_TRTR
            ELSEIF  (iobs == ncbtx ) THEN
              icl  =  IBITS( iodrbufa(iposi_bdy + ncbttr), ntxbp, ntxoc )
              nc_body_buf(irep)% level     (kk) = icl
              nc_body_buf(irep)% level_typ (kk) = VN_TRTR
            ELSEIF  (iobs == ncbtn ) THEN
              icl  =  IBITS( iodrbufa(iposi_bdy + ncbttr), ntnbp, ntxoc )
              nc_body_buf(irep)% level     (kk) = icl
              nc_body_buf(irep)% level_typ (kk) = VN_TRTR
            ELSEIF ((iobs == ncbrad) .OR. (iobs == ncbrdd)                     &
                                     .OR. (iobs == ncbrdt)) THEN
              nc_body_buf(irep)% level     (kk) =  1._sp
              nc_body_buf(irep)% level_typ (kk) = VN_TRTR
            ENDIF
!           IF(nc_body_buf(irep)% level_typ (kk) == VN_TRTR)  nvrx2  = -9
          ENDIF

          ! default status is 'obs_only', if obs operator does not exist,
          !   i.e. for which (iposr_inc == iposr_ver-1),
          !   i.e. for which (ixobs_br* == -1) :
          !        snow depth, pressure tendency, ZPD, STD, vertical wind
          ! default status is 'passive' and 'fl_obstype' set, if obs operator is
          ! not implemented in COSMO but exists in MEC: temporally non-local obs
          !   i.e. for which (iposr_inc == iposr_ver),
          !   i.e. for which (ixobs_br* ==  0) :
          !        precip, gusts, min/max temperature, radiation
          ! --> is status is 'passive' then the 'fl_obstype' flag should be set
          istatus = ST_OBS_ONLY
          IF (nvrx   >= 0) THEN
            istatus = ST_ACTIVE
            IF (.NOT. BTEST( mflgerf, nvrx )) THEN
              ! default is 'passive' for obs for which 'nvrx' does not exist:
              !   vertical wind, cloud cover, ZPD, STD
              istatus = ST_PASSIVE
              IF (nvfxbp >= 0)  istatus = ST_REJECTED
            ENDIF
            ! threshold quality control
            IF (BTEST( mflgqcf, nvrx )) THEN
              IF (istatus == ST_ACTIVE )  istatus = ST_REJECTED
              IF (istatus == ST_PASSIVE)  istatus = ST_PAS_REJ
            ENDIF
          ENDIF
          IF (iposr_inc == iposr_ver-1)  istatus = ST_OBS_ONLY
          IF (iposr_inc == iposr_ver  )  istatus = ST_PASSIVE

          ! quality control flags insert
          kflg = 0
          IF (nvfxbp >= 0) THEN
            ! set dataset flag if 2-bit ODR flag = 1
            IF (IBITS( mflgmfw, nvfxbp + nvfbps(1), nvfboc(1) ) == 1)          &
              kflg = IBSET ( kflg, FL_DATASET )
            kflg = ireplace1( kflg, FL_BLACKLIST, mflgmfw, nvfxbp + nvfbps(2) )
            kflg = ireplace1( kflg, FL_PRACTICE , mflgmfw, nvfxbp + nvfbps(5) )
            kfbit =      ibit1( mflgmfw, nvfxbp + nvfbps(3) )
            IF ((kfmt == 1) .AND. ((nvfxbp == nvfubp) .OR.(nvfxbp == nvftbp))) &
              kfbit = MAX( kfbit , ibit1( mflgmfw, nvfxbp + nvfbps(6) ) )
            kflg = ireplace1( kflg, FL_GROSS    , kfbit  , 0 )
            kfbit = MAX( ibit1( mflgmfw, nvfxbp + nvfbps(4) )                  &
                       , ibit1( mflgmfw, nvflbp             ) )
            kflg = ireplace1( kflg, FL_HEIGHT   , kfbit  , 0 )
          ELSEIF (nvrx == nvrw) THEN
            !   vertical velocity (RASS, w-prof): dataset flag from 'T', 'u,v'
            IF (iodrbufa(iposi + nchcdt) == nra_eu) THEN
              IF (IBITS( mflgmfw, nvftbp + nvfbps(1), nvfboc(1) ) == 1)        &
                kflg = IBSET ( kflg, FL_DATASET )
            ELSE
              IF (IBITS( mflgmfw, nvfubp + nvfbps(1), nvfboc(1) ) == 1)        &
                kflg = IBSET ( kflg, FL_DATASET )
            ENDIF
          ENDIF
          IF (nvrx  >= 0)  kflg = ireplace1( kflg, FL_FG    , mflgqcf, nvrx  )
          IF (nvrx2 >= 0)  kflg = ireplace1( kflg, FL_FG_LBC, mflgqcf, nvrx2 )
!         IF (iposr_inc == iposr_ver)  kflg = IBSET ( kflg, FL_OBSTYPE )
          IF ((istatus == ST_PASSIVE) .OR. (istatus == ST_PAS_REJ))            &
            kflg = IBSET ( kflg, FL_OBSTYPE )

          nc_body_buf(irep)% state     (kk) = istatus
          nc_body_buf(irep)% flags     (kk) = kflg
          nc_body_buf(irep)% check     (kk) = bit1_pos_nr( kflg, nfl_          &
                                                         , kfl_(1:nfl_) )
          IF ((nvfxbp <= -1) .AND. (nvrx <= -1) .AND.(iposr_inc /= iposr_ver)) &
            nc_body_buf(irep)% check   (kk) = FL_NONE

          ! veri-data (model equivalent)
          nc_body_buf(irep)% veri_data (kk,1) = REAL(rmdi,sp)
          ! if variable for the model-equivalent exists, i.e. if (ixobs_br* > 0)

          IF ((nfcast >= 1) .AND. (iposr_inc > iposr_ver))                     &
            nc_body_buf(irep)% veri_data (kk,1) = rodrbufa(iposr_inc)

          ! cancel upper-air height obs if pressure was converted from height
          IF (      (nc_body_buf(irep)% varno     (kk) == VN_HEIGHT)           &
              .AND. (     (nc_body_buf(irep)% level_typ (kk) == VN_HEIGHT)     &
                     .OR. (nc_body_buf(irep)% level_typ (kk) == VN_FLEV ))) THEN
            nc_body_buf(irep)% plevel (kk) = REAL(rmdi,sp)
            kk = kk - 1
            nc_body_buf(irep)% len_body    = nc_body_buf(irep)% len_body  -  1

!CS: may be needed to be introduced for programs that process 'fof'-file
!         ELSEIF (nfcast >= 1) THEN
!           IF (      (nc_body_buf(irep)% veri_data (kk,1) == REAL(rmdi,sp))   &
!               .AND. (nc_body_buf(irep)% state     (kk)   /= ST_OBS_ONLY))    &
!             nc_body_buf(irep)% state (kk) = ST_OBS_ONLY
          ENDIF
        ENDIF
      ENDDO  !  iobs
    ENDDO  !  n_levl

! fill observations of type integer (from 'iodrbufa')
! ---------------------------------------------------
    IF (kfmt == 2) THEN
      ksig = iodrbufa(iposi_b + ncbcsg)
      DO iobs = 1, ileni_b1
        kobv = iodrbufa(iposi_b + iobs)
        IF ((kobv /= imdi) .AND. (ixobs_bis(iobs) >= 0)) THEN
          laddobs = .FALSE.
          IF (iobs == ncbwwe) THEN
            kobv  =  IBITS( kobv, nvw0bp, nvw0oc )
            IF (kobv /= nibits( nvw0oc )) THEN
              kk = kk + 1
              laddobs = .TRUE.
              nc_body_buf(irep)% obs       (kk) = REAL(kobv, sp)
              nc_body_buf(irep)% level     (kk) = REAL(c0, sp)
              nc_body_buf(irep)% level_typ (kk) = VN_NUM
              nc_body_buf(irep)% level_sig (kk) = imdi
            ENDIF
          ELSE
            kk = kk + 1
            laddobs = .TRUE.
            IF (iobs == ncbclg) kcl = 0
            IF (iobs == ncbcl1) kcl = 1
            IF (iobs == ncbcl2) kcl = 2
            IF (iobs == ncbcl3) kcl = 3
            IF (iobs == ncbcl4) kcl = 4
            nc_body_buf(irep)% obs       (kk) = REAL(kobv, sp)
            nc_body_buf(irep)% level     (kk) = REAL(kcl, sp)
            nc_body_buf(irep)% level_typ (kk) = VN_NUM
            nc_body_buf(irep)% level_sig (kk) = IBITS(ksig, ncsgbp(kcl), nxsgoc)
          ENDIF
          IF (laddobs) THEN
            nc_body_buf(irep)% varno     (kk) = iis_varno(iobs)
            nc_body_buf(irep)% bcor      (kk) = REAL(c0, sp)
            nc_body_buf(irep)% e_o       (kk) = REAL(rmdi, sp)
            nc_body_buf(irep)% plevel    (kk) = REAL(rmdi, sp)
            nc_body_buf(irep)% accuracy  (kk) = REAL(rmdi, sp)
            nc_body_buf(irep)% dlat      (kk) = REAL(rmdi, sp)
            nc_body_buf(irep)% dlon      (kk) = REAL(rmdi, sp)
            nc_body_buf(irep)% qual      (kk) = imdi
            nc_body_buf(irep)% state     (kk) = ST_OBS_ONLY
            nc_body_buf(irep)% flags     (kk) = 0
            nc_body_buf(irep)% check     (kk) = FL_NONE
            nc_body_buf(irep)% veri_data (kk,1) = REAL(rmdi, sp)
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDDO ! irep;  n_rep

!-------------------------------------------------------------------------------
!  Section 3:  Fill 'report' which is the input for 'write_report'
!-------------------------------------------------------------------------------

! determine 'nw_rep', 'n_body' (number of new reports / observations)
! -------------------------------------------------------------------
  nw_rep  =  n_rep
  n_body  =  0
  IF (n_rep >= 1)  n_body  =  SUM( nc_body_buf(:)% len_body )

  ! check if there is enough space in the NetCDF file (header and body);
  ! if not then write alerting CAUTION messages and adjust 'nw_rep', 'n_body'
  nexce_rep  =  nexce_rep  +  MAX( 0 , n_tot_verrep  + nw_rep - max_rep )
  nexce_bdy  =  nexce_bdy  +  MAX( 0 , n_tot_verbody + n_body - max_body )
  PRINT  '(6X,A,I6,A,I6,A,I7,A,I8)' ,                                          &
         'feedobs file: # newrep', nw_rep, ', oldrep', n_tot_verrep            &
                      ,', newobs', n_body, ', oldobs', n_tot_verbody
  IF (     (n_tot_verrep  + nw_rep > max_rep)                                  &
      .OR. (n_tot_verbody + n_body > max_body)) THEN
    IF (n_tot_verrep + nw_rep > max_rep) THEN
      icl  = n_tot_verrep + nw_rep
      PRINT    '("CAUTION !!!!! total number of reports ",I6                   &
               &," > FOF size max_rep =",I6)' ,  icl, max_rep 
    ENDIF
    IF (n_tot_verbody + n_body > max_body) THEN
      icl = n_tot_verbody + n_body
      PRINT    '("CAUTION !!!!! total number of obs ",I10                      &
               &," > FOF size max_body =",I10)' ,  icl, max_body
    ENDIF

    ! adjust 'nw_rep' and 'n_body'
    nw_rep  =  MIN( nw_rep , max_rep - n_tot_verrep )
    n_body  =  0
    iobs = nw_rep
    addrep: DO irep = 1 , iobs
      IF (n_tot_verbody + n_body + nc_body_buf(irep)% len_body <= max_body) THEN
        n_body  =  n_body  +  nc_body_buf(irep)% len_body
      ELSE
        nw_rep  =  irep - 1
                                                                     EXIT addrep
      ENDIF
    ENDDO addrep
  ENDIF

! set offsets in FOF body
! ------------------------
  IF (nw_rep >= 1)   koffset (1) = n_tot_verbody
  DO irep = 2 , nw_rep
    koffset (irep) = koffset(irep-1) + nc_body_buf(irep-1)% len_body
  ENDDO

! update 'n_tot_verrep' and 'n_tot_verbody' (# total reports / obs in FOF)
! -------------------------------------------------------------------------
  n_tot_verbody  =  n_tot_verbody  +  n_body
  n_tot_verrep   =  n_tot_verrep   +  nw_rep

! initialize filling of 'report'
! ------------------------------
  IF (nw_rep > 0) THEN
    ALLOCATE( report( nw_rep ), STAT=jerr )
    IF (jerr/=0) THEN
      WRITE( yerr,'("ERROR in memory alloc: report(nw_rep) ",I8) ') nw_rep
                                                                          RETURN
    ENDIF
  ENDIF

  DO irep = 1, nw_rep
    iposi = ioffsrt(irep)
    iposr = iodrbufa(iposi + nchoff)
    kfmt  = iodrbufa(iposi + nchfmt)

    report(irep)% len    = nc_body_buf(irep)% len_body
    report(irep)% offset = koffset( irep )

! fill report(irep)% header
! -------------------------

!   ALLOCATE( report(irep)% header , STAT=jerr )
!   IF (jerr/=0) THEN
!     WRITE( yerr,'("ERROR in memory alloc: report(irep)% header")')
!                                                                         RETURN
!   ENDIF

    report(irep)% header% i_body        = koffset( irep )
    report(irep)% header% l_body        = nc_body_buf(irep)% len_body 
    report(irep)% header% n_level       = n_levl(irep)
    report(irep)% header% data_category = iodrbufa( iposi + nchcat )
    report(irep)% header% sub_category  = iodrbufa( iposi + nchcas )
    report(irep)% header% center        = MOD( iodrbufa( iposi+ nchcen ), 1000 )
    report(irep)% header% sub_center    =      iodrbufa( iposi+ nchcen )/ 1000
    report(irep)% header% obstype       = iodrbufa( iposi + nchobt )
    report(irep)% header% codetype      = iodrbufa( iposi + nchcdt )
    report(irep)% header% ident         = iodrbufa( iposi + nchide )
    report(irep)% header% statid(:)     = " "
    DO icl = 1, MIN( ilstidc , 10 )
      report(irep)% header% statid(icl:icl) = CHAR( iodrbufa(iposi+mxchi(kfmt) &
                                                                  +icl) )        
    ENDDO
    report(irep)% header% lat          = REAL(rodrbufa( iposr + nchlat ), sp)
    report(irep)% header% lon          = REAL(rodrbufa( iposr + nchlon ), sp)
    report(irep)% header% sun_zenit    = REAL(rodrbufa( iposr + nchsol ), sp)
    report(irep)% header% time         = iodrbufa( iposi + nchtim )
    report(irep)% header% time_nomi    = iodrbufa( iposi + nchsyn )
    report(irep)% header% time_dbase   = iodrbufa( iposi + nchtdb )
    report(irep)% header% z_station    = iodrbufa( iposi + nchalt )
    report(irep)% header% z_modsurf    = iodrbufa( iposi + nchsfc )
    report(irep)% header% r_state      = iodrbufa( iposi + nchpas )
    kflg                      = MAX( 0 , iodrbufa( iposi + nchflg ) )
    IF (iodrbufa(iposi + nchflg) == imdi)  kflg = 0
    report(irep)% header% r_flags      = kflg
    report(irep)% header% r_check      = bit1_pos_nr( kflg, nfl_, kfl_(1:nfl_) )
    IF (iodrbufa(iposi + nchflg) == imdi)                                      &
      report(irep)% header% r_check    = FL_NONE
    report(irep)% header% sta_corr     = ibit1( iodrbufa(iposi+nchsch), nvscbp )
    report(irep)% header% mdlsfc       = ibit1( iodrbufa(iposi+nchsch), nvsebp )
    report(irep)% header% index_x      = iodrbufa( iposi + nchi_t )
    report(irep)% header% index_y      = iodrbufa( iposi + nchj_t )
    report(irep)% header% instype      = iodrbufa( iposi + nchtyp )
    kphase  =  IBITS( iodrbufa( iposi+nchsch ), nvapbp, nvapoc-1 )
    IF (      (iodrbufa(iposi+nchobt) /= nairep)                               &
        .AND. ((kphase == nibits(nvapoc-1)) .OR. (kphase == 0)))  kphase = imdi
    report(irep)% header% phase        = kphase
    IF (kfmt == 1) THEN
      report(irep)% header% tracking   = iodrbufa( iposi + nchtrc )
      report(irep)% header% meas_type  = iodrbufa( iposi + nchna4 )
      report(irep)% header% rad_corr   = iodrbufa( iposi + nchrad )
    ELSE
      report(irep)% header% tracking   = imdi
      report(irep)% header% meas_type  = imdi
      report(irep)% header% rad_corr   = imdi
    ENDIF
!   WRITE (0,*) 'ZZFOF ', irep, nc_body_buf(irep)% len_body, n_levl(irep)      &
!              , report(irep)% header% obstype                                 &
!              , report(irep)% header% codetype                                &
!              , report(irep)% header% statid

! fill report(irep)% body
! -----------------------

    ALLOCATE( report(irep)% body(nc_body_buf(irep)% len_body) , STAT=jerr )
    IF (jerr/=0) THEN
      WRITE( yerr,'("ERROR in memory alloc: report(irep)% body(len_body)"      &
                  &,2I8)' ) irep, len_body
                                                                          RETURN
    ENDIF
    IF (report(irep)% offset + report(irep)% len > max_body) THEN
      ! this should never happen as 'max_body' is already checked further above
      WRITE( yerr,'("ERROR: MAX_BODY for feedobs file too small",3I8,I6)' )    &
             report(irep)% offset, report(irep)% len, max_body, irep
      jerr = report(irep)% offset + report(irep)% len - max_body
                                                                          RETURN
    ENDIF

    DO kk = 1, nc_body_buf(irep)% len_body
      report(irep)% body(kk)% varno     = nc_body_buf(irep)% varno    (kk)
      report(irep)% body(kk)% obs       = nc_body_buf(irep)% obs      (kk)
      report(irep)% body(kk)% bcor      = nc_body_buf(irep)% bcor     (kk)
      report(irep)% body(kk)% e_o       = nc_body_buf(irep)% e_o      (kk)
      report(irep)% body(kk)% level     = nc_body_buf(irep)% level    (kk)
      report(irep)% body(kk)% plevel    = nc_body_buf(irep)% plevel   (kk)
      report(irep)% body(kk)% accuracy  = nc_body_buf(irep)% accuracy (kk)
      report(irep)% body(kk)% dlat      = nc_body_buf(irep)% dlat     (kk)
      report(irep)% body(kk)% dlon      = nc_body_buf(irep)% dlon     (kk)
      report(irep)% body(kk)% level_typ = nc_body_buf(irep)% level_typ(kk)
      report(irep)% body(kk)% level_sig = nc_body_buf(irep)% level_sig(kk)
      report(irep)% body(kk)% state     = nc_body_buf(irep)% state    (kk)
      report(irep)% body(kk)% flags     = nc_body_buf(irep)% flags    (kk)
      report(irep)% body(kk)% check     = nc_body_buf(irep)% check    (kk)
      report(irep)% body(kk)% qual      = nc_body_buf(irep)% qual     (kk)
    ENDDO

! fill report(irep)% body(kk)% veri_data  (if n_veri >= 1)
! --------------------------------------

    IF (n_veri >= 1) THEN
      jv = 1
      DO kk = 1, nc_body_buf(irep)% len_body
!!      ALLOCATE ( report(irep)% body(kk)% veri_data(n_veri) , STAT=jerr )
!!      IF (jerr/=0) THEN
!!        WRITE( yerr,'("ERROR in alloc: report(irep)%body(kk)%"               &
!!                    &,"veri_data(n_veri)",I6,I8,I3)' ) irep, kk, n_veri
!!                                                                        RETURN
!!      ENDIF
!!!     DO jv = 1,  n_veri
!!!       report(irep)% body(kk)% veri_data(jv) = nc_body_buf(irep)% veri_data(kk)
!!!     ENDDO
        report(irep)% body(kk)% veri_data  =>  nc_body_buf(irep)% veri_data(kk,:)
!       report(irep)% body(kk)% veri_data(1)  =  nc_body_buf(irep)% veri_data(kk)
!!      report(irep)% body(kk)% veri_data  =  nc_body_buf(irep)% veri_data(kk)
      ENDDO
    ENDIF
  ENDDO  !  nw_rep

!-------------------------------------------------------------------------------
!  Section 4:  Write all obs reports (in 'report') into NetCDF feedobs file
!-------------------------------------------------------------------------------

  IF (nw_rep > 0) THEN
    ihdr_pos  =  n_tot_verrep  - nw_rep + 1
    ibdy_pos  =  n_tot_verbody - n_body + 1
    inw_rep   =  nw_rep
    in_body   =  n_body
    imiss     =  imdi
    rmich     =  REAL(rmdich, sp)

    CALL write_report ( fb, report, inw_rep, in_body, ihdr_pos, ibdy_pos       &
                      , imiss, rmich, jerr, yerr )
!   =================

    IF (jerr/=0)                                                          RETURN
  ENDIF  !  nw_rep

  ! de-allocate 'nc_body_buf': loop over 'n_rep', not 'nw_rep' !
  DO irep = 1, n_rep
    IF (jerr==0) DEALLOCATE( nc_body_buf(irep)% varno     , STAT=jerr )
    IF (jerr==0) DEALLOCATE( nc_body_buf(irep)% obs       , STAT=jerr )
    IF (jerr==0) DEALLOCATE( nc_body_buf(irep)% bcor      , STAT=jerr )
    IF (jerr==0) DEALLOCATE( nc_body_buf(irep)% e_o       , STAT=jerr )
    IF (jerr==0) DEALLOCATE( nc_body_buf(irep)% level     , STAT=jerr )
    IF (jerr==0) DEALLOCATE( nc_body_buf(irep)% plevel    , STAT=jerr )
    IF (jerr==0) DEALLOCATE( nc_body_buf(irep)% accuracy  , STAT=jerr )
    IF (jerr==0) DEALLOCATE( nc_body_buf(irep)% dlat      , STAT=jerr )
    IF (jerr==0) DEALLOCATE( nc_body_buf(irep)% dlon      , STAT=jerr )
    IF (jerr==0) DEALLOCATE( nc_body_buf(irep)% state     , STAT=jerr )
    IF (jerr==0) DEALLOCATE( nc_body_buf(irep)% flags     , STAT=jerr )
    IF (jerr==0) DEALLOCATE( nc_body_buf(irep)% check     , STAT=jerr )
    IF (jerr==0) DEALLOCATE( nc_body_buf(irep)% qual      , STAT=jerr )
    IF (jerr==0) DEALLOCATE( nc_body_buf(irep)% level_typ , STAT=jerr )
    IF (jerr==0) DEALLOCATE( nc_body_buf(irep)% level_sig , STAT=jerr )
    IF (jerr==0) DEALLOCATE( nc_body_buf(irep)% veri_data , STAT=jerr )
    IF (jerr/=0) THEN
      WRITE( yerr,'("ERROR in memory de-alloc: nc_body_buf")')
                                                                          RETURN
    ENDIF
  ENDDO  !  n_rep

!-------------------------------------------------------------------------------
!  Section 5:  Clean up
!-------------------------------------------------------------------------------

! de-allocate arrays and derived types
! ------------------------------------

  ! de-allocate 'report'
  DO irep = 1, nw_rep
!   IF (n_veri >= 1) THEN
!     DO kk = 1, report(irep)% len
!       IF (jerr==0) DEALLOCATE( report(irep)% body(kk)% veri_data , STAT=jerr )
!     ENDDO
!   ENDIF
    IF ((jerr/=0) .AND. (kerr == 0)) kerr = 41
    IF (jerr==0) DEALLOCATE( report(irep)% body   , STAT=jerr )
    IF ((jerr/=0) .AND. (kerr == 0)) kerr = 42
!   IF (jerr==0) DEALLOCATE( report(irep)% header , STAT=jerr )
!   IF ((jerr/=0) .AND. (kerr == 0)) kerr = 43
  ENDDO
  IF (nw_rep > 0) THEN
    IF (jerr==0) DEALLOCATE( report , STAT=jerr )
    IF ((jerr/=0) .AND. (kerr == 0)) kerr = 44
  ENDIF
  ! de-allocate 'nc_body_buf' and 'n_levl'
  IF (n_rep > 0) THEN
    IF (jerr==0) DEALLOCATE( nc_body_buf , STAT=jerr )
    IF ((jerr/=0) .AND. (kerr == 0)) kerr = 45
    IF (jerr==0) DEALLOCATE( koffset     , STAT=jerr )
    IF ((jerr/=0) .AND. (kerr == 0)) kerr = 46
    IF (jerr==0) DEALLOCATE( n_levl      , STAT=jerr )
    IF ((jerr/=0) .AND. (kerr == 0)) kerr = 47
  ENDIF
  ! error messages
  IF (jerr/=0) THEN
    IF (kerr==41) WRITE( yerr,'("ERROR in memory de-alloc: *%body% veri_data")')
    IF (kerr==42) WRITE( yerr,'("ERROR in memory de-alloc: report% body     ")')
    IF (kerr==43) WRITE( yerr,'("ERROR in memory de-alloc: report% header   ")')
    IF (kerr==44) WRITE( yerr,'("ERROR in memory de-alloc: report           ")')
    IF (kerr==45) WRITE( yerr,'("ERROR in memory de-alloc: nc_body_buf      ")')
    IF (kerr==46) WRITE( yerr,'("ERROR in memory de-alloc: koffset          ")')
    IF (kerr==47) WRITE( yerr,'("ERROR in memory de-alloc: n_levl           ")')
                                                                          RETURN
  ENDIF

! print statistics and alerting CAUTION messages (if size of FOF is too small)
! ----------------------------------------------
  IF ((lastpr) .AND. (lfdbk_open)) THEN
    IF (LEN_TRIM(yfdbkfile) <= 47) THEN
      PRINT   '(A,A,A)' , "    NetCDF FEEDOBS FILE: "                          &
                        , yfdbkfile(1:LEN_TRIM(yfdbkfile)) , " created"
    ELSEIF (LEN_TRIM(yfdbkfile) <= 120) THEN
      PRINT   '(A)' , "    NetCDF FEEDOBS FILE created: "
      PRINT   '(A)' ,  yfdbkfile(1:MIN(120, PLEN, LEN_TRIM(yfdbkfile)))
    ELSEIF (PLEN >= 121) THEN
      PRINT   '(A)' , "    NetCDF FEEDOBS FILE created: "
      PRINT   '(A)' ,  yfdbkfile(1:120)
      PRINT   '(A)' ,  yfdbkfile(121:MIN(240, PLEN, LEN_TRIM(yfdbkfile)))
    ENDIF
    PRINT   '(A,I7,A,I8)' , "    Total number of reports in fdbk file:"        &
                          , n_tot_verrep,  " ; max allowed : ", max_rep
    PRINT   '(A,I8,A,I8)' , "    Total number of bodys in fdbk file: "         &
                          , n_tot_verbody, " ; max allowed : ", max_body

! set some global attributes and close FOF
! -----------------------------------------
    jerr          = nf_redef ( fb% nc% ncid )
    IF ((jerr/=0) .AND. (kerr == 0)) kerr = 51
    fb% n_hdr     = n_tot_verrep
    fb% nc% error = nf_put_att_int ( fb%nc% ncid, NF_GLOBAL, 'n_hdr',NF_INT, 1 &
                                   , fb% n_hdr )
    fb% n_body    = n_tot_verbody
    fb% nc% error = nf_put_att_int ( fb%nc% ncid, NF_GLOBAL,'n_body',NF_INT, 1 &
                                   , fb% n_body )
    jerr          = nf_enddef ( fb% nc% ncid )
    IF ((jerr/=0) .AND. (kerr == 0)) kerr = 52
    IF (jerr/=0) THEN
      IF (kerr == 51)                                                          &
        WRITE( yerr,'("ERROR in nf_redef  in obs_write_cdf_feedobs")' )
      IF (kerr == 52)                                                          &
        WRITE( yerr,'("ERROR in nf_enddef in obs_write_cdf_feedobs")' )
                                                                          RETURN
    ENDIF

    CALL close_fdbk   (fb)
!   ===============
    CALL cleanup_fdbk (fb)
!   =================

    lfdbk_open = .FALSE.
    PRINT   '("feedback file closed")'
  ENDIF !  lastpr .AND. lfdbk_open

  IF ( (lastpr) .AND. (lfdbk_crea) )                                           &
    PRINT   '("NOTE !!! feedobs file not created as there was no observation")' 

!ENDIF  !  (my_cart_id == ircv)


!-------------------------------------------------------------------------------
! End of module procedure obs_write_cdf_feedobs
!-------------------------------------------------------------------------------
END SUBROUTINE obs_write_cdf_feedobs


!===============================================================================

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

ELEMENTAL INTEGER FUNCTION ireplace1   ( invar, ipos, irepl, ipsr )
  !----------------------------------------------------------------
  INTEGER (KIND=iintegers)  , INTENT (IN)  ::  invar, ipos, irepl, ipsr
  !-----------------------------------------------------------------------------
  ! replaces the bit at position 'ipos' of integer word 'invar'
  ! by the bit at position 'ipsr' from integer word 'irepl'
  !-----------------------------------------------------------------------------
  !
  ireplace1 = IOR( IAND( invar, NOT( ISHFT( 1, ipos ) ) )                      &
                 , ISHFT( IAND( ISHFT( irepl,-ipsr ), 1 ), ipos ) )
  !
END FUNCTION ireplace1

!-------------------------------------------------------------------------------

END MODULE src_obs_cdfout_feedobs
