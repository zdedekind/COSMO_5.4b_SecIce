!+ Source module for the observation processing in the data assimilation mode
!-------------------------------------------------------------------------------

MODULE src_obs_fdbk_in

!-------------------------------------------------------------------------------
! Description:
!   This module performs the reading of observations from a NetCDF feedobs
!   (feedback) file, as well as their pre-selection (according to observation
!   time, location, type, and code type), distribution and storage in COSMO's
!   internal data structre ODR.
!   As an option, the simulated observations from the feedobs file are written
!   to the ODR as observed values, so that they can serve in an OSSE setup as
!   synthetic observations derived from a nature run. As an additional option,
!   these observed values are randomly perturbed.
!
! Method:
!   This module contains the following module procedures:
!    - obs_fof_read           : called by obs_cdf_read_org (src_obs_cdfin_org.f90)
!    - obs_fof_store_1rep     : called by obs_fof_head
!    - obs_fof_mask_reports   : called by obs_fof_head
!    - obs_fof_perturb        : called by obs_fof_head
!    - construct_seed         : called by obs_fof_head
!
!   This module also contains elemental functions, formerly statement functions:
!   - insert       : inserts bits set in one integer into another integer word
!                    at given bit position
!
!   It uses from:
!    - src_obs_cdfin_util: obs_assign_gridpt : assigns obs horiz'ly. to grid pt.
!    - mo_fdbk_io:       - read_fdbk_head    : reads feedback file header
!                        - read_fdbk_body    : reads feedback file body
!                        - read_fdbk_veri    : reads verification data
!                        - pack_fdbk         : packs header + body acc. to masks
!                      \ - associate_hb      : associates related header + body
!                        - scatter_fdbk_data : scatter type t_fdbk_data to PE's
!    - mo_random:        - construct         : new random generator from seed
!                        - random_gauss      : gets Gaussian random numbers
!    - utilities:        - phi2phirot        : converts lat to rotated coord.
!                        - rla2rlarot        : converts lon to rotated coord.
!                        - get_utc_date      : derived date in different forms
!                        - diff_minutes      : difference between 2 dates/times
!    - parallel_utilities: distribute_values : distributes values to all PEs
!    - environment:      - model_abort       : aborts program in case of errors
!    - netcdf
!
!   Data modules used:
!    - data_parameters
!    - data_obs_lib_cosmo
!    - data_obs_cdfin 
!    - data_obs_record
!    - mo_fdbk_tables
!    - mo_netcdf_param
!    - mo_fdbk_io
!    - kind_parameters
!    - mo_random
!
!   Note: This module is part of the 'COSMO data assimilation library 1' 
!         for reading data from NetCDF observation input files.
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V4_28        2013/07/12 Christoph Schraff
!  Initial release.
! V5_1         2014-11-28 Oliver Fuhrer, Christoph Schraff
!  Replaced ireals by wp (working precision) (OF), Ulrich Schaettler
!  Replaced mo_kind by kind_parameters (US)
!  Old random number seed construction (used for observation perturbation)
!  replaced by the more general routine used for SPPT (stochastic physics).
! V5_3         2015-10-09 Christoph Schraff
!  Processing of additional elements from surface reports for verification
!  (model equivalent calculator) purposes: wind speed + direction, dewpoint
!  temperature, radiation, 3-hourly precip and wind gusts. Max. and min. T-2m
!  over arbitrary periods instead of fixed 12 hours.
! V5_4         2016-03-10 Christoph Schraff
!  Dimension of 'neventr' and 'neventd' reduced from 3 to 2.
!  Variable related to AOF interface removed.
! V5_4a        2016-05-10 Ulrich Schaettler
!  Bug fix: wp was used from data_parameters and kind_parameters
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
    iintegers    ! KIND-type parameter for standard integer variables

!-------------------------------------------------------------------------------

USE data_obs_lib_cosmo, ONLY :  &

! 1. General parameters
! ---------------------

!   c0         ,& ! standard real constant 0.0
    c1         ,& ! standard real constant 1.0
!   c2         ,& ! standard real constant 2.0
    c05        ,& ! standard real constant 0.5
!   c3600      ,& ! standard real constant 3600.0
    rmdi       ,& ! =-1.E31_wp : commonly used missing data indicator
    rmdich     ,& ! =-1.E30_wp : commonly used check value for miss data
    epsy       ,& ! = 1.E-8_wp : commonly used very small value > 0

! 2. Scalar variables originally defined in other data modules,
!    obtained by calling 'obs_cdf_interface'
! -------------------------------------------------------------

    ! horizontal and vertical sizes of the model fields
    ke             ,& ! number of grid pts in vertical direction (--> 'f_z2p')
    ie_tot         ,& ! number of grid pts in zonal direction (in total domain)
    je_tot         ,& ! number of grid pts in meridional dir. (in total domain)

    ! variables related to parallelisation / domain decomposition
    num_compute    ,& ! number of compute PEs
    nboundlines    ,& ! number of overlapping boundary lines of the subdomains
    my_cart_id     ,& ! rank of this subdomain in the cartesian communicator
    icomm_cart     ,& ! communicator for the virtual cartesian topology
    imp_reals      ,& ! REAL      type used for MPI
    imp_integers   ,& ! INTEGER   type used for MPI

    ! other variables related to namelist parameters
    maxmlv         ,& ! size (level  dimension) of the  multi-level (m-l)  ODR
    maxmll         ,& ! size (report dimension) of the  multi-level (m-l)  ODR
    maxsgl         ,& ! size (report dimension) of the single-level (s-l)  ODR
    maxgpl         ,& ! size (report dimension) of the (ground-based) GPS  ODR
    nolbc          ,& ! number of grid rows at lateral boundaries
                      !   where obs are neglected
    irun_osse      ,& ! model run to derive obs values from file yfofin='fof':
                      !   == 0 : obs from yfofin are used as obs
                      !   >= 1 : simulated obs from file yfofin with model run
                      !          index = irun_osse are used as obs  (thus, run
                      !          'irun_osse' is the nature run for the OSSE)
    losse_fg       ,& ! if true  then first guess check flag from 'fof' is
                      !          converted into 'dataset' pre-processing flag
                      ! if false then first guess check flag is discarded
                      !          and obs may be used actively
    iseed          ,& ! external seed for random number generator

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
    acthr          ,& ! actual model time [hours] with respect to 'yydate_ref'
    fperturb       ,& ! factor to obs error variances to define size of random 
                      ! perturbations added to the obs (only from yfofin='fof')

    ! switches related to namelist parameters and other
    lverpas        ,& ! write also passive reports on feedback files
    ydate_ref      ,& ! reference date (e.g. start of the forecast)
                      ! yyyymmddhhmmss (year, month, day, hour, min., sec.)

! 2b. Pointers for arrays originally defined in other data modules
! ----------------------------------------------------------------

    ! array related to parallelisation / domain decomposition
    i_subpos       ,& ! positions of the subdomains in the total domain
                      ! (i-, j-indices of the lower left + upper right grid pt.
                      !  in the order: i_ll, j_ll, i_ur, j_ur; only the domain
                      !  interior is considered, not the boundary lines)

    ! model fields
    r_hhl          ,& ! height   (at half levels)
    r_ps              ! surface pressure

USE data_obs_lib_cosmo, ONLY :  &

! 4. I/O device numbers for obs processing / nudging and file names
! -----------------------------------------------------------------
    yucautn    ,& ! caution messages if too many obs for ODR size
    yurejct    ,& ! direct reporting of rejected obs. reports
    nurej      ,& ! direct reporting of rejected obs. reports
    nucautn    ,& ! caution messages if too many obs for current ODR size
    lopen_rej  ,& ! .true. if yurejct is open

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

    ! 1.1  internal attributes of the different NetCDF input files
    ! ------------------------------------------------------------

    icdfinlen  ,& ! maximum length of NetCDF observation input file name
    iannexlen  ,& ! maximum length of annex of NetCDF obs input file name
    yfofin     ,& ! file name of feedobs (feedback) input file(s) (='fof')
    yfofannex  ,& ! annex of NetCDF feedobs (feedback) input file names
    ichar10    ,& ! character length of station identity from feedback files

    ! 2    Blacklist and Whitelist
    ! ----------------------------

!   ilstid_blk ,& ! assume 8-character station-IDs in Black-/Whitelist
!   blk_loc    ,& ! blacklists for local reports

    ! 3.1  Report event counter array format
    ! --------------------------------------

    nenoal     ,& ! no station altitude
    nedbfl     ,& ! data base flag on loc/tim/alt high
    neloca     ,& ! station location out of range
    nezdif     ,& ! distance 'model orography - station altitude' too large
    neblak     ,& ! blacklisted ship
    neobct     ,& ! observation or code type excluded on area with sta. location
    nesodr     ,& ! report number exceeding size of ODR (==> adjust Namelist)
    nenoda     ,& ! no accepted data in report
    neredn     ,& ! redundancy between 2 multi-level, or 2 single-level reports
!   neredx     ,& ! redundancy between 1 multi- and 1 single-level report
!   nenoml     ,& ! multi-levl report not made due to ODR array size
    netrac     ,& ! (flight) track suspicious
    nethin     ,& ! thinning of aircraft (flight) track
    nelodr     ,& ! level rejected: number of levels exceeding ODR size

    ! 3.2  Event counter arrays
    ! -------------------------

    neventr    ,& ! counter array of report events
    neventd    ,& ! counter array of data events

    ! 7    For data rejection messages: Output buffer, size and formats
    ! -----------------------------------------------------------------

    outbuf     ,& ! buffer containing output for a single node
    nacout     ,& ! actual number of records stored in the output buffer
    nmxoln     ,& ! maximum length of output buffer
    istrej     ,& ! length of strings (station id) in output buffer
    nfmt4      ,& ! excess of levels

    ! 8    Temporary global model fields
    ! ----------------------------------

    hsurf_tot  ,& ! total array of model surface height
    fland_tot     ! total array of fraction of land in each grid element

! end of data_obs_cdfin

!-------------------------------------------------------------------------------

USE data_obs_record, ONLY :   &

!       1.     Header formats of ODR reports
!       ------------------------------------

!       1.1.1  Header formats of ODR reports: 'omlhed' and 'osghed'
!              ----------------------------------------------------
!   mxrhed     ,& ! header length of multi-level reports
!   mxshed     ,& ! header length of single-level reports
!   mxghed     ,& ! header length of GPS reports
    nhilon     ,& ! longitude of observing station
    nhjlat     ,& ! latitude  of observing station
    nhalt      ,& ! station altitude [m]
    nhtime     ,& ! (exact) time of observation in forecast hours
    nhsurf     ,& ! height of model grid pt. to which obs. is assigned
    nhzio      ,& ! longitude of obs. station (or lowest datum) in grid pt. unit
    nhzjo      ,& ! latitude  of obs. station in grid pt. units
    nhsynt     ,& ! nominal (synoptic) time of observation in forecast hours
    nhtddb     ,& ! data base decoding time in forecast hours
    nhsolz     ,& ! solar zenith angle [deg]
    nhvcbu     ,& ! correction factor to vertical correlation scale for wind
                  ! at base of report
    nhvcbt     ,& ! as 'nhvcbu', but for temperature
    nhvcbq     ,& ! as 'nhvcbu', but for humidity
    nhvctu     ,& ! correction factor to vertical correlation scale for wind
                  ! at top of report
    nhvctt     ,& ! as 'nhvctu', but for temperature
    nhvctq        ! as 'nhvctu', but for humidity
!   nhtvip     ,& ! obs multi-level pressure interpol. to lowest model level
!   nhtviz        ! vertical distance to nearest observation

USE data_obs_record, ONLY :   &

!       1.1.2  Header formats of ODR reports: 'momlhd' and 'mosghd'
!              ----------------------------------------------------
!   mxrhdf     ,& ! header length of multi-level reports
!   mxshdf     ,& ! header length of single-level reports
!   mxghdf     ,& ! header length of GPS reports
    nhio       ,& ! (local) x-coord. of grid pt. assigned to obs
    nhjo       ,& ! (local) y-coord. of grid pt. assigned to obs
    nhitot     ,& ! global x-coord. of grid pt. assigned to obs
    nhjtot     ,& ! global y-coord. of grid pt. assigned to obs
    nhobtp     ,& ! observation type
    nhcode     ,& ! code type
    nhschr     ,& ! station characteristics                      (see 1.1.4)
    nhpass     ,& ! flag for report being set to 'passive'       (see 1.1.4)
    nhqcfw     ,& ! status of QC and of writing to feedobs files (and p-QC flag)
    nhflag     ,& ! report flags (obs type, surface, altitude, station ID)
    nhcorr     ,& ! update sequence number (station correction indicator)
    nhcat      ,& ! data     category (from BUFR Section 1)
    nhcats     ,& ! data sub-category (from BUFR Section 1)
    nhkz       ,& ! DWD internal classification number (observation type)
    nhcent     ,& ! originating centre
    nhstid     ,& ! station identity number
    nhdate     ,& ! absolute exact observation date [yyyymmdd]
    nhhrmn     ,& ! absolute exact observation time [hhmm]
    nhsyhr     ,& ! absolute nominal (synoptic) observation time [yymmddhh]
    nhnlev     ,& ! number of obs. levels (for multi-level reports)
    nhvqcf     ,& ! for satellite retrieval: threshold quality control flags
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
    nhstyp     ,& ! surface obs: station type (buoy: MQOBL, BUFR Table 002149,
                  !                            else: NIX  , BUFR Table 002001)

!       1.1.3  Header formats of ODR reports: 'yomlhd' and 'yosghd'
!              ----------------------------------------------------

    ilstid        ! character length of the station identity

USE data_obs_record, ONLY :   &

!       1.2    Bit patterns for packed information in ODR (and VOF) header
!              -----------------------------------------------------------
!   nvpabp     ,& ! bit pos. for report passive due to being merged   nhschr
!   nvpsbp     ,& ! bit pos. for report set passive due to next 5 flags "
!   nvobbp     ,& ! bit pos. for flag: 'station loc. out of spec. area' "
    nvalbp     ,& ! bit pos. for flag: 'height distance too large'      "
!   nvbkbp     ,& ! bit pos. for flag: 'blacklisted station (ship)'     "
!   nvexbp     ,& ! bit pos. for flag: 'obs type excluded in area'      "
!   nvrdbp     ,& ! bit pos. for flag: 'redundant report'               "
    nvsebp     ,& ! bit pos. for report located at sea grid pt.         "
    nvscbp     ,& ! bit pos. for station correction indicator           "
    nvapbp     ,& ! bit pos. for phase of flight (aircraft)             "
    nvapoc     ,& ! no. of bits occ. by phase of flight                 "
    nvaabp     ,& ! bit pos. for aircraft roll angle (code)             "
!   nvaaoc     ,& ! no. of bits occ. by aircraft roll angle             "

!       1.3    ODR body format
!              ---------------

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
    nbtzio     ,& ! longitude in grid pt. units
    nbtzjo     ,& ! latitude  in grid pt. units
    nbttim     ,& ! observation time relative to report (header) time
    nbtlop     ,& ! LOG( pressure )
    nbtdrh     ,& ! bias correction for relative humidity [/]
    nbtw       ,& ! vertical velocity [m/s]
    nbtsnr     ,& ! signal to noise ratio 
!   nbtuac     ,& ! accuracy (std dev from data provider) of horiz. wind [m/s]


!       1.3.2  Body format of ODR of multi-level report flags: 'momlbd'
!              --------------------------------------------------------
!   mxrbdf     ,& ! body length of multi-level reports
    nbtflg     ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbterr     ,& ! status flag word        (bit pattern, see below: 'nb?err')
    nbtqcf     ,& ! threshold quality control flags      (see below: 'nb?qcf')
    nbtlsg     ,& ! level id (bit pattern, as in NetCDF statistics file)
    nbtlid        ! level identity          (bit pattern, see below: 'nb?lid')

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

    nbsrr1     ,& ! precipitation amount over 1 hour                   [mm]
    nbsrr3     ,& ! precipitation amount over 3 hours                  [mm]
    nbsrr6     ,& ! precipitation amount over 6 hours                  [mm]
    nbsr12     ,& ! precipitation amount over 12 hours                 [mm]
    nbsr24     ,& ! precipitation amount over 24 hours                 [mm]
!   nbsfgv     ,& ! max. derived equivalent vertical gust (aircraft)   [m/s]
    nbsfg1     ,& ! max. wind speed of gusts over 1 hour               [m/s]
    nbsfg3     ,& ! max. wind speed of gusts over 3 hours              [m/s]
    nbsfg6     ,& ! max. wind speed of gusts over 6 hours              [m/s]
    nbstn      ,& ! minimum temperature (at 2m, in period 'nbsttr')    [K]
    nbstx      ,& ! maximum temperature (at 2m, in period 'nbsttr')    [K]
    nbsrad     ,& ! global    solar    radiation, sum over 1 hour      [J/m2]
    nbsrdd     ,& ! diffuse   solar    radiation, sum over 1 hour      [J/m2]
    nbsrdt     ,& ! long-wave downward radiation, sum over 1 hour      [J/m2]
    nbshsw     ,& ! total snow depth                                   [m]
    nbsdrh     ,& ! bias correction for relative humidity [/]
    nbsviz        ! scaled extrapolat. distance for surf. pressure obs [m]

USE data_obs_record, ONLY :   &

!       1.3.4  Body format of ODR of surface report flags: 'mosgbd'
!              ----------------------------------------------------
!   mxsbdf     ,& ! body length of single-level reports
    nbsflg     ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbserr     ,& ! status flag word        (bit pattern, see below: 'nb?err')
    nbsqcf     ,& ! threshold quality control flags      (see below: 'nb?qcf')
    nbslid     ,& ! SYNOP: pressure code (SYNOP)   (code, see below: 'nbslid')
                  ! else : level identity   (bit pattern, see below: 'nb?lid')
!   nbscwg     ,& ! combined cloud and weather group (set of classes, s below)
    nbswwe     ,& ! NetCDF read, SYNOP: weather and ground group word  (below)
!   nbstur     ,& ! NetCDF read, Aircraft: degree of turbulence WMO Tab 011031
                  !   (not contained in merged multi-level aircraft reports !)
    nbsclg     ,& ! general           cloud       group (code)
    nbscl1     ,& ! first  individual cloud layer group (code)
    nbscl2     ,& ! second individual cloud layer group (code)
    nbscl3     ,& ! third  individual cloud layer group (code)
    nbscl4     ,& ! forth  individual cloud layer group (code)
    nbsttr     ,& ! time periods    (bit pattern of values, see below: nbsttr)

!       1.3.5  Body format of ODR of GPS reports: 'ogpbdy'
!              -------------------------------------------
!   mxgbdy     ,& ! body length of GPS reports
    nbgtze     ,& ! error in total zenith delay                        [mm]
    nbgzpd     ,& ! zenith path delay (total zenith delay)             [mm]
    nbgzwd     ,& ! zenith wet delay                                   [mm]
    nbgiwv     ,& ! integrated water vapour                            [mm]
    nbgp       ,& ! pressure                                           [Pa]
    nbgt       ,& ! temperature                                        [K]
    nbgrh      ,& ! relative humidity                                  [/]
    nbgbia     ,& ! bias correction to integrated water vapour         [mm]
    nbgiwa     ,& ! adjusted (bias corrected) integrated water vapour  [mm]
!   nbgz       ,& ! height
    nbgzer     ,& ! error of observed height
    nbgqer     ,& ! error of observed relative humidity
    nbgdrh     ,& ! bias correction for observed relative humidity
    nbgter     ,& ! error of observed temperature

!       1.3.6  Body format of ODR of GPS report flags: 'mogpbd'
!              ------------------------------------------------
!   mxgbdf     ,& ! body length of GPS reports
    nbgflg     ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbgerr     ,& ! pre-processing status flags  (bit p., see below: 'nb?err')
    nbgqcf     ,& ! threshold quality control flags      (see below: 'nb?qcf')
    nbglid        ! level identity          (bit pattern, see below: 'nb?lid')

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
    nvriwv     ,& ! bit pos. for status/QC flags for IWV                 "
    nvrzpd     ,& ! bit pos. for status/QC flags for zenith path delay   "
    nvrspd     ,& ! bit pos. for status/QC flags for slant path delay    "
    nvrct      ,& ! bit pos. for status/QC flags for (total) cloud       "
    nvrcl      ,& ! bit pos. for status/QC flags for low     cloud       "
    nvrcm      ,& ! bit pos. for status/QC flags for middle  cloud       "
    nvrch      ,& ! bit pos. for status/QC flags for high    cloud       "
    nvrcbs     ,& ! bit pos. for status/QC flags for cloud base height   "
    nvfubp     ,& ! bit pos. for main flag on wind                     nb?flg
    nvftbp     ,& ! bit pos. for main flag on temperature                "
    nvfqbp     ,& ! bit pos. for main flag on humidity                   "
    nvfzbp     ,& ! bit pos. for main flag on pressure / geopot.         "
    nvfgbp     ,& ! bit pos. for main flag on integr. water vapour       "
    nvfaoc     ,& ! no. of bits occ. by each main flag                   "
    nvfbps     ,& ! bit pattern for main flags:                          "
    nvfboc     ,& ! no. of bits occ. for main flags                      "
    nvflbp     ,& ! bit pos. for level flag: level below surface         "
!   nvfloc     ,& ! no. of bits occ. by level flag                       "
    ntxbp      ,& ! bit pos. for time period for T-max (1,..,24 hr)   "
    ntnbp      ,& ! bit pos. for time period for T-min (1,..,24 hr)   "
    ntxoc      ,& ! no. of bits occ. by each time period              "
    nvlidp     ,& ! level id. bit pattern                              nb?lid
!   nvlido     ,& ! no. bits occ. by each indicator in level id.         "

!       1.4.3  Bit patterns for 'optional groups' in ODR body 'mosgbd' (and VOF)
!              -----------------------------------------------------------------
                  ! --> weather and ground group word                  nbswwe
    nvw0bp     ,& ! bit position for ww (present wea.) code (WMO Table 020003)
                  !
                  ! --> general    cloud       group word              nbsclg
    nxsgbp     ,& ! bit position for vertic. signific. code (WMO Table 008002)
                  !
                  ! --> individual cloud layer group words             nbscl?
    nxsibp     ,& ! bit position for vertic. signific. code (WMO Table 008002)
    nxsgoc        ! no. bits occupied for vert. signf. code (WMO Table 008002)

USE data_obs_record, ONLY :   &

!       1.5    Further quantities related to ODR
!              ---------------------------------
    imdi       ,& ! missing data indicator for ODR integers (2^31-1)
    ntotml     ,& ! total number of stored multi-level reports
    ntotsg     ,& ! total number of stored single-level reports
    ntotgp     ,& ! total number of stored GPS ZPD/IWV reports
    fdoro      ,& ! scaling factor to vertical distances betw. model

!       2.     Observation data records (ODR)
!       -------------------------------------

    omlbdy     ,& ! body   of multi-level ODR
    omlhed     ,& ! header of multi-level ODR
    osgbdy     ,& ! body   of single-level ODR
    osghed     ,& ! header of single-level ODR
    ogpbdy     ,& ! body   of GPS ODR
    ogphed     ,& ! header of GPS ODR
    momlbd     ,& ! body   of multi-level ODR
    momlhd     ,& ! header of multi-level ODR
    mosgbd     ,& ! body   of single-level ODR
    mosghd     ,& ! header of single-level ODR
    mogpbd     ,& ! body   of GPS ODR
    mogphd     ,& ! header of GPS ODR
    yomlhd     ,& ! header of multi-level ODR
    yosghd     ,& ! header of single-level ODR
    yogphd     ,& ! header of GPS ODR

!       3.     Masking constants
!       ------------------------

    nibits        ! masking constants

! end of data_obs_record

!-------------------------------------------------------------------------------

  USE mo_fdbk_tables,          ONLY :  &

    OT_SYNOP   ,& ! SYNOP report
    OT_AIREP   ,& ! AIREP report
    OT_SATOB   ,& ! SATOB report (AMV)
    OT_DRIBU   ,& ! DRIBU report
    OT_TEMP    ,& ! TEMP  report
    OT_PILOT   ,& ! PILOT report
    OT_SATEM   ,& ! SATEM report
    OT_PAOB    ,& ! PAOB  report
    OT_SCATT   ,& ! Scatterometer report
    OT_RAD     ,& ! Radiances
    OT_GPSRO   ,& ! GPS Radio occultations
    OT_GPSGB   ,& ! GPS ground based observations
    OT_RADAR   ,& ! RADAR (volume data)
    OC_ATSCD   ,& ! automatic synop surface report
    OC_AHSCD   ,& ! ship synop report
!!! OC_ABSCD   ,& ! ship synop abbreviated report
!!! OC_SHRED   ,& ! shred reportKp synop report
    OC_ATSHS   ,& ! automatic ship synop report
    OC_WP_EU   ,& ! European  wind profiler
    OC_RA_EU   ,& ! European sodar/rass report
    OC_WP_JP   ,& ! Japanese  wind profiler
    OC_PR_US   ,& ! wind/profiler/rass report (USA)
    OC_RAVAD      ! radar VAD wind profile report

  USE mo_fdbk_tables,          ONLY :  &

    VT_FIRSTGUESS  ,& ! run type: first guess
!   VT_ANALYSIS    ,& ! run type: analysis
!   VT_FORECAST    ,& ! run type: forecast
    VE_DETERM      ,& ! ensemble member: derterministic model run
    ST_ACTIVE      ,& ! status: used in the assimilation
    ST_MERGED      ,& ! status: not used, merged into multi-level report
    ST_PASSIVE     ,& ! status: not used, only monitored
!   ST_REJECTED    ,& ! status: not used due to suspicious quality
!   ST_PAS_REJ     ,& ! status: passive and rejected
!   ST_OBS_ONLY    ,& ! status: obs. only, no model equivalent available
    FL_OBSTYPE     ,& ! passive report type (at obs. location)
    FL_BLACKLIST   ,& ! blacklist (or not whitelist)  (in particular:
                      !   bad (hard blacklisted) aircraft station identity)
    FL_SUSP_LOCT   ,& ! suspicious location (3D)
    FL_AREA        ,& ! location not in valid area (i.e. outside user-def. area)
    FL_HEIGHT      ,& ! location not in valid height range
    FL_SURF        ,& ! incorrect surface (land,ice,etc), e.g. sea obs over land
    FL_PRACTICE    ,& ! bad reporting practice / insufficient data
    FL_DATASET     ,& ! dataset qality flags
    FL_REDUNDANT   ,& ! redundant report
    FL_FLIGHTTRACK ,& ! flight track error flag
    FL_MERGE       ,& ! (report used for) merged report (only, e.g. TEMP ABCD
                      !   or single-level aircraft used for multi-level rep.)
    FL_THIN        ,& ! thinning
    FL_GROSS       ,& ! gross error flag
    FL_NO_OBS      ,& ! no (active) observations in report
    FL_FG          ,& ! observation minus first guess check
    FL_FG_LBC      ,& ! observation minus lateral boundary field check
    FL_NONE        ,& ! no flag set
    LS_SURFACE     ,& ! surface
    LS_STANDARD    ,& ! standard level
    LS_TROPO       ,& ! tropopause level
    LS_MAX         ,& ! maximum wind level
    LS_SIGN           ! significant level

  USE mo_fdbk_tables,          ONLY :  &
                      !   variable numbers:
    VN_U           ,& ! u-component of wind
    VN_V           ,& ! v-component of wind
    VN_PWC         ,& ! precipitable water content kg/m**2
    VN_RH          ,& ! relative humidity
    VN_RH2M        ,& ! 2 metre relative humidity
    VN_T           ,& ! upper air temperature
    VN_T2M         ,& ! 2 metre temperature
    VN_TD2M        ,& ! 2 metre dew point
    VN_PTEND       ,& ! pressure tendency
    VN_WW          ,& ! present weather
    VN_VV          ,& ! visibility
    VN_NH          ,& ! cloud base height
    VN_N_L         ,& ! low cloud amount
    VN_SDEPTH      ,& ! snow depth
    VN_TRTR        ,& ! time period of information  (h)
    VN_RR          ,& ! precipitation amount        (kg/m^2)
    VN_JJ          ,& ! maximum temperature         (K)
    VN_GCLG        ,& ! general cloud group         Table
    VN_N           ,& ! total cloud amount          (20011)
    VN_PS          ,& ! surface pressure            (Pa)
    VN_DD          ,& ! wind direction              (degree)
    VN_FF          ,& ! wind force                  (m/s)
!   VN_RAWBT       ,& ! brightness temperature      (K)
    VN_U10M        ,& ! 10m u-component of wind     (m/s)
    VN_V10M        ,& ! 10m v-component of wind     (m/s)
    VN_W           ,& ! vertical speed              (m/s)
!   VN_VT          ,& ! virtual temperature         (K)
    VN_HEIGHT      ,& ! height                      (m)
    VN_FLEV           ! nominal flight level        (m)

  USE mo_fdbk_tables,          ONLY :  &
!   VN_RREFL       ,& ! radar reflectivity          (Db)
!   VN_PDELAY      ,& ! atmospheric path delay      (m)
    VN_ICLG        ,& ! individual cloud layer group Table
    VN_N_M         ,& ! middle cloud amount          WMO table 020011
    VN_N_H         ,& ! high   cloud amount          WMO table 020011
    VN_RAD_GL      ,& ! 1-h global solar radiation  (J/m2) 
    VN_RAD_DF      ,& ! 1-h diffuse solar radiat.   (J/m2) 
    VN_RAD_LW      ,& ! 1-h long-wave radiation     (J/m2) 
    VN_ZPD         ,& ! zenith path delay
    VN_ZWD         ,& ! zenith wet delay
!   VN_SPD         ,& ! slant path delay
    VN_GUST        ,& ! wind gust                   (m/s)
    VN_P           ,& ! pressure                    (Pa)
    VN_TMIN        ,& ! minimum Temperature         (K)
    VN_NUM            ! ordinal (channel) number    (  )

! end of mo_fdbk_tables

!-------------------------------------------------------------------------------

  USE mo_netcdf_param

!-------------------------------------------------------------------------------

  USE mo_fdbk_io,              ONLY :  &
    t_fdbk_data     ,& ! derived type to hold feedback file content
    t_fdbk_head     ,& ! component of 't_fdbk_data' (header table entry)
    t_fdbk_body     ,& ! component of 't_fdbk_data' (body   table entry)
    ifill, rfill    ,& ! fillvalues used generally in this module
    read_fdbk_head  ,& ! read feedback file header
    read_fdbk_body  ,& ! read feedback file body
    read_fdbk_veri  ,& ! read verification data
    pack_fdbk       ,& ! pack header and body according to masks
!   associate_hb    ,& ! associate related header and body entries
    scatter_fdbk_data  ! scatter type 't_fdbk_data' to PE's

!-------------------------------------------------------------------------------

  USE kind_parameters,         ONLY :  &
    wp              ,& ! working precision real kind
    dp                 ! double precision (required for random_gauss)

!-------------------------------------------------------------------------------

  USE mo_random,               ONLY :  &
    random_state_t  ,& ! derived type: random generator state
    construct       ,& ! constructs new random generator (state) from seed
!   destruct        ,& ! destructs random generator
    random_gauss       ! generator for Gaussian random number distribution

!-------------------------------------------------------------------------------

  USE environment,              ONLY :  &
    model_abort        ! aborts the program in case of errors

!-------------------------------------------------------------------------------

  USE utilities,                ONLY :  &
    get_utc_date    ,& ! actual date of the forecast in different forms
    phi2phirot      ,& ! converts phi from the real to the rotated system
    rla2rlarot      ,& ! converts lambda from the real to the rotated system
!   uv2uvrot_vec    ,& ! converts wind components from normal to rotated system
    diff_minutes       ! compute difference in minutes between 2 dates/times

!-------------------------------------------------------------------------------

 USE parallel_utilities,       ONLY :  &
    distribute_values  ! distributes a set of values from one node to all others

!-------------------------------------------------------------------------------

  USE src_obs_cdfin_util,       ONLY :  &
!   obs_assign_sort_node   ,& ! assign node to reports and sort them accordingly
    obs_assign_gridpt         ! horiz. assignment of an obs report to a grid pt.

!-------------------------------------------------------------------------------

  USE src_stoch_physics,        ONLY :  &
    set_seed_rand_numb        ! sets a seed for a random number stream as a
                              ! function of (reference) date / time and an
                              ! external seed number

!-------------------------------------------------------------------------------

! USE src_obs_cdfin_blk,           ONLY :  &
!   obs_cdf_whitelist_local,& ! indic. for local reports if missing on whitelist
!   obs_cdf_blacklist_local   ! produces a list of blacklisted vertical
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
!+ Module procedure in "src_obs_fdbk_in" for reading and distributing reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_fof_read ( min_sta, min_end, ilfof, ycdfdir, icdfdirlen , nodrnew, nexceed )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_fdbk_in" organizes the reading
!   from a feedobs (feedback) file, pre-selection, distribution and storage
!   in the ODR of observations from a given observation time period.
!   As an option, the simulated observations from the feedobs file are written
!   to the ODR as observed values, so that they can serve in an OSSE setup as
!   synthetic observations derived from a nature run. As an additional option,
!   these observed values may be randomly perturbed.
!
! Method:
!   The reports are read by 1 node (processor) by call of feedback file reading
!   routines from the 3DVAR environment (module 'mo_fdbk_io.f90'). Observation
!   time, location, type, and code type are checked to be in the required range.
!   If required, the grid point assignment of the observations to the model
!   grid points is applied.
!   If required, the simulated observations are re-labelled as observations,
!   and Gaussian random perturbations are added using a random error generator
!   from the 3DVAR environment (Donald E. Knuth's portable pseudo-random number
!   generator, module 'mo_random').
!   By use of another routine of the 3DVAR environment, the observation reports
!   are then distributed to those sub-domains which contain the grid point they
!   are assigned to.
!   Finally, the observational reports in feedback file format are converted
!   into the COSMO data structure ODR in parallel mode (section 2) for long-term
!   internal storage.
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

  INTEGER (KIND=iintegers) , INTENT (IN)     ::  &
    min_sta       ,& ! start  \  of time interval to be read now
    min_end       ,& ! end    /  (in [minutes] of forecast time)
    icdfdirlen    ,& ! max. length of name of directory where
                     !   NetCDF observation input files reside
    ilfof            ! index (number) of NetCDF feedobs (feedback) input file

  CHARACTER (LEN= *)       , INTENT (IN)     ::  &
    ycdfdir          ! directory where NetCDF obs input + feedobs files reside

  INTEGER (KIND=iintegers) , INTENT (INOUT)  ::  &
    nodrnew (4)   ,& ! number of new reports
    nexceed (3)      ! number of reports in excess of array size

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  TYPE (t_fdbk_data)  :: fdbk    ! feedobs (feedback) file and data structure;
                                 !   dimensions etc. are set in 'read_fdbk_head'

  INTEGER (KIND=iintegers) ::  &
    ianyy   ,imoyy   ,& ! year     / of the date and time
    ianmm   ,imomm   ,& ! month   /  of the current
    iandd   ,imodd   ,& ! day    <   observation
    ianhh   ,imohh   ,& ! hour    \  resp. of the initial
    ianmin  ,imomin  ,& ! minute   \ model date and time
    imindif          ,& ! time difference [min] between these two times
    ierrf               ! error status for time difference routine

  INTEGER (KIND=iintegers) ::  &
    ilenpath       ,& ! length of path of feedobs file
    mt_sta         ,& ! start of required time range [min rel. to ydate_ref]
    mt_end         ,& ! end   of required time range [min rel. to ydate_ref]
    kveri          ,& ! veri_data index for the input nature run in an OSSE
    irep           ,& ! loop index over reports
    iobt           ,& ! index of first observation of report in 'fdbk' body
    kobtyp         ,& ! CMA observation type
    kcdtyp         ,& ! CMA observation code type
    nlev           ,& ! number of vertical levels
    n_bdy          ,& ! number of observations (body elements) in 'fdbk'
    nt             ,& ! total number of reports (of type kfmt)
    icma           ,& ! pointer index for 'cma' / diagnostic arrays
!   inode          ,& ! node to which the current report is assigned
    ipe            ,& ! loop index over nodes
    ix     (1)     ,& ! model run index in feedobs file
    ndistr (2)     ,& ! distributed values
    ierr              ! error indicators

  LOGICAL                  ::  &
    lassign           ! compute assignment of obs report to model grid point

  CHARACTER (LEN=MIN( icdfdirlen+icdfinlen+iannexlen, 255 ))  ::  &
    yfofpath         ! file name of feedobs (feedback) input file, with annex
  CHARACTER (LEN=25)       :: &
    yroutine          ! name of this subroutine

! Local arrays:
! ------------

  LOGICAL                 , ALLOCATABLE :: &
    lt_mask    (:)    ! if true then report is in required time range

  INTEGER (KIND=iintegers), ALLOCATABLE :: &
    kfmt       (:) ,& ! type of ODR report format:
                      !   = 0 : upper-air single-level ;   = 0 : multi-level
                      !   = 2 : surface-level          ;   = 3 : GPS IWV
    nexce      (:)    ! if = 1 then report is in excess of ODR array size

  REAL    (KIND=wp)       , ALLOCATABLE :: &
    zverdum    (:)    ! dummy for simulated obs from nature run in 'fdbk'

!   INTEGER :: itype_osse = 2  ! == 0 : obs of fdbk are used as obs
                               ! == 1 : simulated obs of 'fdbk' are used as obs;
                               !        obs without simulated obs are discarded
                               ! == 2 : simulated obs of 'fdbk' are used as obs;
                               !        obs without simulated obs used as obs
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_fof_read
!-------------------------------------------------------------------------------

  yroutine = 'obs_fof_read'

! PRINT *, yroutine, kcdftyp, min_sta, min_end

!-------------------------------------------------------------------------------
! Section 1: Read and pre-process header information common to all obs types:
!            time / lat / lon / station altitude / obs type / station ID /
!            center of origin ...;
!            assign observations to grid points and to nodes (sub-domains)
!-------------------------------------------------------------------------------

  IF (my_cart_id == 0) THEN

    ! get full path of input feedobs (feedback) file
    yfofpath = TRIM( yfofin ) // TRIM( yfofannex(ilfof) )
    ilenpath = LEN_TRIM(ycdfdir)
    IF (ilenpath >= 1) THEN
      IF (ycdfdir(ilenpath:ilenpath) == '/')  ilenpath = ilenpath - 1
      yfofpath = ycdfdir(1:ilenpath) // '/' // TRIM( yfofpath )
    ENDIF

    ! read meta data and all header elements from feedobs file
    ! --------------------------------------------------------

    CALL read_fdbk_head ( fdbk, yfofpath )
  ! ===================

    ! determine which reports have observation time within required range
    ! -------------------------------------------------------------------

    !   initial (reference) date and time of the forecast run
    READ (ydate_ref,'(I4,4I2)') imoyy, imomm, imodd, imohh, imomin
    !   reference date and time of the feedobs file
!   READ ( fdbk% f% refdate,'(I4,2I2  )' )  ianyy, ianmm, iandd
!   READ ( fdbk% f% reftime,'(   2I2.2)' )  ianhh, ianmin
    ianyy   =        fdbk% f% refdate / 10000
    ianmm   =  MOD ( fdbk% f% refdate / 100   , 100 )
    iandd   =  MOD ( fdbk% f% refdate         , 100 )
    ianhh   =        fdbk% f% reftime / 100
    ianmin  =  MOD ( fdbk% f% reftime         , 100 )

    !   time difference in minutes between feedobs ref time and model ref time
    CALL diff_minutes ( imoyy, imomm, imodd, imohh, imomin,                    &
                        ianyy, ianmm, iandd, ianhh, ianmin, 0, imindif, ierrf )
  ! =================
    !   required obs time interval in units of feedobs file variable 'time'
    mt_sta  =  min_sta - imindif
    mt_end  =  min_end - imindif

    !   allocate and set 'lt_mask'
    ALLOCATE ( lt_mask ( fdbk% f% n_hdr ) )
    lt_mask = .FALSE.
    WHERE (      (fdbk% h% time >= mt_sta)                                     &
           .AND. (fdbk% h% time <= mt_end))  lt_mask = .TRUE.

    ! update 'lt_mask' to discard report for following reasons
    !   - observation / code type unknown to COSMO
    !   - station altitude of surface report missing
    !   - station location missing or outside model domain
    ! also do horizontal grid point assignment and update
    ! fdbk% h% index_x , index_y, z_modsurf, mdlsfc , if required
    ! -----------------------------------------------------------

    CALL obs_fof_mask_reports ( fdbk% f% n_hdr, yfofpath , fdbk, imindif       &
                              , lt_mask , lassign )
  ! =========================

    ! pack 'fdbk' such that it only contains reports in required time range
    ! (note that after 'pack_fdbk( mask_h )', use 'ib' instead of 'i_body';
    !  also use 'SIZE( fdbk% h )' instead of 'fdbk% f% n_hdr')
    ! ---------------------------------------------------------------------

    CALL pack_fdbk ( fdbk , mask_h = lt_mask )
  ! ==============
    DEALLOCATE ( lt_mask )

    ! read all body elements of the reports within
    ! the target time range from the feedobs file
    ! --------------------------------------------

    CALL read_fdbk_body ( fdbk )
  ! ===================

  ! CALL associate_hb ( fdbk )
  ! =================

    !  after this, use 'SIZE( fdbk% b )' instead of 'fdbk% f% n_body' !

!-------------------------------------------------------------------------------
! Section 2: OSSE mode: label (non-missing) simulated obs as obs
!                       and perturb them if required
!-------------------------------------------------------------------------------

    IF (irun_osse > 0) THEN

      ! read verification data, i.e. simulated obs
      ! ------------------------------------------
      ! (other option: e.g. CALL read_fdbk_veri (fdbk, run_type = VT_FIRSTGUESS
      !                                              , ens_member = VE_DETERM) )
      ix (1) = irun_osse

      CALL read_fdbk_veri ( fdbk, fill=rfill,  ix = ix )
    ! ===================

      ! for OSSE, exchange obs with simulated obs from nature run in 'fdbk'
      ! -------------------------------------------------------------------
      !   note: only 1 run is read into 'fdbk% veri_data', hence it resides
      !         at run index = 1 (i.e. kveri = 1)
      kveri = 1
      ALLOCATE ( zverdum ( SIZE( fdbk% b ) ) )
      zverdum (:)  =  fdbk% veri_data(:,kveri)
!     IF (itype_osse == 1) THEN
!       !    new 'fdbk' shall only contain the entries (obs variables)
!       !    with non-missing simulated obs values from nature run
!       fdbk% veri_data (:,kveri)  =  fdbk% b(:)% obs  
!       fdbk% b(:)% obs            =  zverdum(:)

!       CALL pack_fdbk ( fdbk , mask_b = (zverdum /= rfill) )
!     ! ==============

        !  (note that after 'pack_fdbk( mask_b )', use 'nb' instead of 'l_body')

!     ELSEIF (itype_osse == 2) THEN
!       !    exchange the entries with non-missing simulated obs of nature run,
!       !    and leave the original obs in the other entries in the new 'fdbk'
        WHERE ( zverdum(:) /= rfill )
          fdbk% veri_data (:,kveri)  =  fdbk% b(:)% obs  
          fdbk% b(:)% obs            =  zverdum(:)
        END WHERE
!     ENDIF
      DEALLOCATE ( zverdum )
    ENDIF

    ! perturb observations with obs errors, if required (e.g. for OSSE)
    ! -------------------------------------------------

    IF ((fperturb > epsy) .AND. (SIZE( fdbk% b ) > 0)) THEN
      n_bdy = SIZE( fdbk% b )

      CALL obs_fof_perturb  ( SIZE( fdbk% b ), fperturb, iseed , fdbk% b )
    ! ====================

    ENDIF

!-------------------------------------------------------------------------------
! Section 3: Scatter reports to PE's of sub-domains which they reside in
!-------------------------------------------------------------------------------

    ! determine target node for each report
    ! -------------------------------------
    !   (dim. 1 of pointer 'i_subpos' really is: (0:num_compute-1)
    IF (SIZE( fdbk% h ) > 0) THEN
      DO ipe = 0 , num_compute-1
        DO irep = 1 , SIZE( fdbk% h )
          IF (      (i_subpos(ipe,1) <= fdbk% h(irep)% index_x)                &
              .AND. (i_subpos(ipe,3) >= fdbk% h(irep)% index_x)                &
              .AND. (i_subpos(ipe,2) <= fdbk% h(irep)% index_y)                &
              .AND. (i_subpos(ipe,4) >= fdbk% h(irep)% index_y))               &
            fdbk% h(irep)% pe  =  ipe
        ENDDO
      ENDDO
    ENDIF

    yfofpath = TRIM( yfofin ) // TRIM( yfofannex(ilfof) )
    irep     = SIZE( fdbk% h )
    PRINT'(" processing",I6," reports from",I5," - ",I4," [min] from file ",A)'&
          ,irep, min_sta, min_end, yfofpath(1:LEN_TRIM(yfofpath))
  ENDIF

  ! scatter reports
  ! ---------------

  IF (num_compute > 1) THEN

    CALL scatter_fdbk_data ( fdbk , 0, icomm_cart, num_compute, my_cart_id )
  ! ======================

  ! scatter auxiliary information
  ! -----------------------------

    ndistr (:) = 0
    IF  (my_cart_id == 0)                   ndistr (1) = imindif
    IF ((my_cart_id == 0) .AND. (lassign))  ndistr (2) = 1

    CALL distribute_values ( ndistr , 2, 0, imp_integers, icomm_cart, ierr )
  ! ======================

    imindif  =  ndistr(1)
    lassign  =  (ndistr(2) == 1)
  ENDIF

!-------------------------------------------------------------------------------
! Section 4: Fill ODR
!-------------------------------------------------------------------------------

  ALLOCATE ( kfmt ( SIZE( fdbk% h ) ) )

  ! determine ODR report format 'kfmt': upper-air single-level; multi-level;
  ! ----------------------------------  surface; GPS IWV

  DO irep = 1 , SIZE( fdbk% h )
    !   note that it has already been checked that obs / code types are known
    kobtyp = fdbk% h(irep)% obstype
    kcdtyp = fdbk% h(irep)% codetype
    nlev   = fdbk% h(irep)% n_level
    !   determine type of report
    IF ((nlev >= 2) .AND. (     (kobtyp == OT_TEMP) .OR. (kobtyp == OT_PILOT)  &
                           .OR. (kobtyp == OT_AIREP))) THEN
      !   multi-level obs report
      kfmt (irep)  =  1
    ELSEIF ((kobtyp == OT_AIREP) .OR. (kobtyp == OT_SATOB)) THEN
      !   upper-air single-level obs report
      kfmt (irep)  =  0
    ELSEIF  (kobtyp == OT_GPSGB) THEN
      !   ground-based GPS (IWV) obs report
      kfmt (irep)  =  3
    ELSEIF ((kobtyp == OT_SYNOP) .OR. (kobtyp == OT_DRIBU)                     &
                                 .OR. (kobtyp == OT_SCATT)) THEN
      !   surface-level obs report
      kfmt (irep)  =  2
    ELSEIF (     (kcdtyp == OC_WP_EU) .OR. (kcdtyp == OC_RA_EU)                &
            .OR. (kcdtyp == OC_RAVAD) .OR. (kcdtyp == OC_PR_US)                &
            .OR. (kcdtyp == OC_WP_JP)) THEN
      !   profiler reports are always treated as multi-level (even if nlev = 1)
!     kfmt (irep)  =  0
      kfmt (irep)  =  1
    ELSEIF ((kobtyp == OT_TEMP ) .OR. (kobtyp == OT_PILOT)) THEN
      !   TEMP and PILOT balloon reports:
      !   surface reports are derived in the COSMO obs reading routines and have
      !   to be distinguished from regular upper-air reports with only 1 level
      !   -->  check if first observation is an upper-air variable
      !        (surface variables are VN_U10M, VN_T2M, VN_PS, VN_N_L, etc.)
      iobt = fdbk% h(irep)% ib + 1
      if (     (fdbk% b(iobt)% varno == VN_U     )                             &
          .OR. (fdbk% b(iobt)% varno == VN_V     )                             &
          .OR. (fdbk% b(iobt)% varno == VN_T     )                             &
          .OR. (fdbk% b(iobt)% varno == VN_RH    )                             &
          .OR. (fdbk% b(iobt)% varno == VN_HEIGHT)                             &
          .OR. (fdbk% b(iobt)% varno == VN_W     )) THEN
!         .OR. (fdbk% b(iobt)% varno == VN_VGUST )) THEN
!       kfmt (irep)  =  0
        kfmt (irep)  =  1
      ELSE
        kfmt (irep)  =  2
      ENDIF
    ELSE
      kfmt (irep)  =  -1
    ENDIF
  ENDDO
  DO irep = 1 , SIZE( fdbk% h )
    IF (kfmt(irep) == -1) THEN
      !   obs type unknown or not allowed for current obs operators
                                                                          RETURN
    ENDIF 
  ENDDO

  ! fill ODR, if array size large enough
  ! ------------------------------------

  ALLOCATE ( nexce ( SIZE( fdbk% h ) ) )
  nexce (:) = 0
  DO irep = 1 , SIZE( fdbk% h )

    !   multi-level reports
    IF (kfmt(irep) == 1) THEN
      IF (ntotml+nodrnew(1) < maxmll) THEN
        nodrnew (1)  =  nodrnew(1) + 1
        nt           =  ntotml  +  nodrnew(1)
        nlev         =  MIN( fdbk% h(irep)% n_level , maxmlv )

        CALL obs_fof_store_1rep( irep, fdbk, kfmt(irep), nlev, imindif,lassign &
                               , yomlhd(nt)                                    &
                               , omlhed(nt       ,:), momlhd(nt       ,:)      &
                               , omlbdy(nt,1:nlev,:), momlbd(nt,1:nlev,:))
      ! =======================

        !   this is only to indicate whether need to call 'obs_air_org_mult'
        IF (fdbk% h(irep)% obstype == OT_AIREP)  nodrnew (4)  =  nodrnew(4) + 1
      ELSE
        nexce (irep) = 1
      ENDIF

    !   (upper-air or surface) single-level reports
    ELSEIF ((kfmt(irep) == 0) .OR. (kfmt(irep) == 2)) THEN
      IF (ntotsg+nodrnew(2) < maxsgl) THEN
        nodrnew (2)  =  nodrnew(2) + 1
        nt           =  ntotsg  +  nodrnew(2)

        CALL obs_fof_store_1rep ( irep, fdbk, kfmt(irep), 1, imindif, lassign  &
                                , yosghd(nt)                                   &
                                , osghed(nt   ,:), mosghd(nt   ,:)             &
                                , osgbdy(nt:nt,:), mosgbd(nt:nt,:) )
      ! =======================

        IF (fdbk% h(irep)% obstype == OT_AIREP)  nodrnew (4)  =  nodrnew(4) + 1
      ELSE
        nexce (irep) = 1
      ENDIF

    !   ground-based GPS ZPD reports
    ELSEIF (kfmt(irep) == 3) THEN
      IF (ntotgp+nodrnew(3) < maxgpl) THEN
        nodrnew (3)  =  nodrnew(3) + 1
        nt           =  ntotgp  +  nodrnew(3)

        CALL obs_fof_store_1rep ( irep, fdbk, kfmt(irep), 1, imindif, lassign  &
                                , yogphd(nt)                                   &
                                , ogphed(nt   ,:), mogphd(nt   ,:)             &
                                , ogpbdy(nt:nt,:), mogpbd(nt:nt,:) )
      ! =======================

      ELSE
        nexce (irep) = 1
      ENDIF
    ENDIF
  ENDDO

  ! if ODR array size exceeded: update counters,
  !                             also to prepare control messages
  ! ------------------------------------------------------------

  IF (MAXVAL( nexce ) > 0) THEN
    DO irep = 1 , SIZE( fdbk% h )
      IF (nexce(irep) >= 1) THEN
        IF (kfmt(irep) >= 1)  nexceed (kfmt(irep)) = nexceed(kfmt(irep)) + 1
        IF (kfmt(irep) == 0)  nexceed (2)          = nexceed(2)          + 1
        kobtyp = fdbk% h(irep)% obstype
        kcdtyp = fdbk% h(irep)% codetype
        !   get diagnostic array position
        icma  =  i_cma ( kobtyp , kcdtyp )
               ! =====
        !   event 'insufficient ODR size' is updated even for passive reports !!
        neventr (nesodr,icma) = neventr(nesodr,icma) + 1
        IF (fdbk% h(irep)% r_state <  ST_PASSIVE) THEN
          cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
        ELSE
          cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
        ENDIF
        cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
      ENDIF
    ENDDO
  ENDIF

  DEALLOCATE ( kfmt  )
  DEALLOCATE ( nexce )

!-------------------------------------------------------------------------------
! End Subroutine obs_fof_read
!-------------------------------------------------------------------------------
END SUBROUTINE obs_fof_read


!===============================================================================
!+ Module procedure in "src_obs_fdbk_in" for reading and distributing reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_fof_store_1rep ( irep, fdbk, kfmt, nlev, imindif, lassign       &
                              , yzobhd, zobhed, mzobhd, zobbdy, mzobbd )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_fdbk_in" writes 1 report from the
!   feedobs file structure 'fdbk' into the ODR.
!   Limitations:
!   - All ODR elements are filled except for:
!     - most flags in the station characteristics word
!     - level identity word (nbtlid) for multi-level reports, except surface bit
!     - accuracy (std. dev. from data provider) of horizontal wind
!   - Longitude, latitude, and time variations within multi-level reports are
!     zeroed.
!   - Lapse rate and wind shear errors are indicated only as gross error.
!   - Level below surface flag 'nvflbp' is set only as invalid height flag for
!     all obs variables.
!
! Method:
!   4 types of reports are distinguished (according to 'kfml'):
!   upper-air single-level, multi-level, surface, and GPS IWV reports.
!
! Written by        : DWD, Christoph Schraff, 04.01.2013
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER (KIND=iintegers) , INTENT (IN)     ::  &
    irep          ,& ! index of report in 'fdbk' to be stored in ODR
    kfmt          ,& ! type of ODR report format:
                     !   = 0 : upper-air single-level ;   = 0 : multi-level
                     !   = 2 : surface-level          ;   = 3 : GPS IWV
    nlev          ,& ! number of vertical levels in report
    imindif          ! difference in minutes between fdbk  reference time
                     !                           and model reference time

  TYPE    (t_fdbk_data)    , INTENT (IN)     ::  &
    fdbk             ! feedback file data structure

  LOGICAL                  , INTENT (IN)     ::  &
    lassign          ! assignment of obs report to model grid point computed

  CHARACTER (LEN= *)       , INTENT (INOUT)  ::  &
    yzobhd           ! COSMO ODR header, character part (for station id)

  INTEGER (KIND=iintegers) , INTENT (INOUT)  ::  &
                     !          1st dimension:  (1:nlev)
                     !          2nd dimension:  according to array format 
                     !         array format for kfmt =   1    /  0, 2  /   3
    mzobhd   (:)  ,& ! COSMO ODR header, integer part (momlhd / mosghd / mogphd)
    mzobbd (:,:)     ! COSMO ODR body  , integer part (momlbd / mosgbd / mogpbd)

  REAL    (KIND=wp   )     , INTENT (INOUT)  ::  &
                     !      array format for kfmt =   1    /  0, 2  /   3
    zobhed   (:)  ,& ! COSMO ODR header, real part (omlhed / osghed / ogphed)
    zobbdy (:,:)     ! COSMO ODR body  , real part (omlbdy / osgbdy / ogpbdy)

! Local parameters:
! ----------------

  INTEGER (KIND = iintegers) , PARAMETER  :: &
    itype_calendar = 0    ! for specifying the calendar used in 'get_utc_date'

  REAL    (KIND = wp)        , PARAMETER  :: &
    c60   =   60.0_wp

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    kobtyp      ,& ! observation type
    kcdtyp      ,& ! observation code type
    mintob      ,& ! observation time [min., rel. to model ref. time]
    kphase      ,& ! flight phase
    nbxflg      ,& ! ODR index for main flag word
    nbxerr      ,& ! ODR index for pre-processing status flag word
    nbxqcf      ,& ! ODR index for threshold quality control flag word
    nbxlid      ,& ! ODR index for level identity word
    nbxx        ,& ! ODR index for observed variable
    nbxxer      ,& ! ODR index for observation error
    nbxbia      ,& ! ODR index for bias correction (rel. humidity or IWV)
    nvfxbp      ,& ! ODR index for main flag
    nvrx        ,& ! ODR index for status / QC flags
    ilev        ,& ! vertical level index
    iobs        ,& ! loop index over observations (within 1 report)
    icl         ,& ! loop index over single characters
    iov         ,& ! observed value, integer quantity
    kfl         ,& ! fdbk observation flags
    kflfg       ,& ! fdbk observation flags without first guess check flag
    ilevsig     ,& ! fdbk observation level significance
    mexi   (3)  ,& ! active (+1) or passive (-1) obs exist in report
    icma        ,& ! pointer index for 'cma' / diagnostic arrays
    iactr       ,& ! = 1 if report is active, = 0 otherwise
    nempty      ,& ! number of vertical levels without pressure value
    ivar        ,& ! variable index for namelist parameters
    nflgx       ,& ! total flag for an observation of a variable
    ilen        ,& ! length of control message
    nzday          ! Julian day

  LOGICAL                  ::  &
    lraso       ,& ! true, if balloon sounding (TEMP or PILOT)
    lev_new     ,& ! true, if new vertical level (defined only for kfmt <= 1)
    lexact      ,& ! true, if active data exist in report
    lseaonly    ,& ! use actively only observations from sea platforms
    landsy         ! land synoptic report

  REAL (KIND=wp)           ::  &
    zpp    , zzz    ,& ! pressure / height level
    roblat , roblon ,& ! geographical latitude and longitude
    rlat   , rlon   ,& ! latitude and longitude in the rotated system
    zhr             ,& ! hour
    fisd            ,& ! height diff. betw. model orography and station altitude
    fisdzz          ,& ! scaled extrapolation distance for surface press. value
    fisdnu          ,& ! scaled extrapolation distance for full values and
                       !   increments of surface pressure (for nudging only)
    zzkeml          ,& ! height of the lowest main model level
    zporo              ! (approx.) model surface pressure


! CHARACTER (LEN=25)       :: &
!   yroutine            ! name of this subroutine
  CHARACTER (LEN=14)  :: yakdat1  ! actual date in the form   yyyymmddhh
  CHARACTER (LEN=28)  :: yakdat2  ! actual date, form:   wd   dd.mm.yy  hh UTC

  TYPE (t_fdbk_head)  , POINTER ::  fh        ! pointer to fdbk body data
  TYPE (t_fdbk_body)  , POINTER ::  bd (:)    ! pointer to fdbk body data
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_fof_store_1rep
!-------------------------------------------------------------------------------

  kobtyp = fdbk% h(irep)% obstype
  kcdtyp = fdbk% h(irep)% codetype

  lraso  =       ((kobtyp == OT_TEMP ) .OR. (kobtyp == OT_PILOT))              &
           .AND. (kcdtyp /= OC_WP_EU) .AND. (kcdtyp /= OC_RA_EU)               &
           .AND. (kcdtyp /= OC_RAVAD) .AND. (kcdtyp /= OC_PR_US)               &
           .AND. (kcdtyp /= OC_WP_JP)

!-------------------------------------------------------------------------------
! Section 1: Observation header: fill ODR arrays 'zobhed', 'mzobhd', 'yzobhd'
!-------------------------------------------------------------------------------

  !   for convenience: set a local pointer to the fdbk body of current report

  fh  =>  fdbk% h (irep)

  zobhed   (:) = rmdi
  mzobhd   (:) = imdi

  !   (note: - 'time', 'lat', 'lon', 'obstype', 'codetype' have been checked
  !            before to have non-missing values within reasonable limits;
  !          - for 'index_x', 'index_y', 'z_modsurf', 'mdlsfc' reasonable
  !            values have been produced previously (if required);
  !          - it is assumed that 'r_flags' contains reasonable values !)

                               yzobhd   =  fh% statid (1:MIN(ichar10,ilstid))
                               zobhed (nhtime) = (fh% time       + imindif) /c60
  IF (fh% time_nomi  /= ifill) zobhed (nhsynt) = (fh% time_nomi  + imindif) /c60
  IF (fh% time_dbase /= ifill) zobhed (nhtddb) = (fh% time_dbase + imindif) /c60
                               zobhed (nhjlat) = fh% lat
                               zobhed (nhilon) = fh% lon
  IF (fh% z_station  /= ifill) zobhed (nhalt ) = REAL( fh% z_station , wp )
                               zobhed (nhsurf) = REAL( fh% z_modsurf , wp )
  IF (fh% sun_zenit  /= rfill) zobhed (nhsolz) = fh% sun_zenit
                               mzobhd (nhitot) = ABS ( fh% index_x )
                               mzobhd (nhjtot) = ABS ( fh% index_y )
                               mzobhd (nhobtp) = fh% obstype
                               mzobhd (nhcode) = fh% codetype
                               mzobhd (nhflag) = fh% r_flags
                               mzobhd (nhqcfw) = -1
                               mzobhd (nhcorr) = 0
  IF (fh% sta_corr   /= ifill) mzobhd (nhcorr) = fh% sta_corr
  IF (fh%data_category/=ifill) mzobhd (nhcat ) = fh% data_category
  IF (fh% sub_category/=ifill) mzobhd (nhcats) = fh% sub_category
                               mzobhd (nhcent) = 255
  IF (fh% center     /= ifill) mzobhd (nhcent) = MOD( fh% center  , 1000 )
  IF (fh% sub_center /= ifill) mzobhd (nhcent) =   mzobhd(nhcent)              &
                                                 + fh% sub_center * 1000
  IF (fh% dbkz       /= -1   ) mzobhd (nhkz  ) = fh% dbkz
  IF (fh% ident      /= ifill) mzobhd (nhstid) = fh% ident
  IF (fh% r_state == ST_MERGED) THEN
                               mzobhd (nhpass) = 1
  ELSEIF ((fh% r_state <= ST_ACTIVE) .AND. (fh% r_state >= 0)) THEN
                               mzobhd (nhpass) = 0
  ELSE
                               mzobhd (nhpass) = 2
  ENDIF
  mzobhd (nhio  ) = mzobhd(nhitot) - i_subpos(my_cart_id,1) + 1 + nboundlines
  mzobhd (nhjo  ) = mzobhd(nhjtot) - i_subpos(my_cart_id,2) + 1 + nboundlines

  !   convert into rotated coordinates
  roblat = zobhed(nhjlat)
  roblon = zobhed(nhilon)
  rlon = rla2rlarot (roblat, roblon, r_pollat, r_pollon, r_polgam)
  rlat = phi2phirot (roblat, roblon, r_pollat, r_pollon)
  !   observation position in grid point units (total area !)
  zobhed (nhzio ) = 1._wp + (rlon - r_startlon_tot) /r_dlon
  zobhed (nhzjo ) = 1._wp + (rlat - r_startlat_tot) /r_dlat
  mintob = fh% time       + imindif

  CALL get_utc_date ( mintob, ydate_ref, c60, itype_calendar, yakdat1          &
                    , yakdat2, nzday, zhr)
! =================

  READ( yakdat1(1:8) ,'(I8)' )  mzobhd (nhdate)
! READ( yakdat1(9:12),'(I4)' )  mzobhd (nhhrmn)
  yakdat2 (9:14) = yakdat1(9:10) // '0000'
  READ( yakdat2(9:12),'(I4)' )  mzobhd (nhhrmn)
  mintob = fh% time_nomi  + imindif

  CALL get_utc_date ( mintob, ydate_ref, c60, itype_calendar, yakdat1          &
                    , yakdat2, nzday, zhr)
! =================

  READ( yakdat1(3:10),'(I8)' )  mzobhd (nhsyhr)

  !   the station characteristics element 'nhschr' is not used in the
  !   obs operators, but flags 'nvscbp', 'nvsebp', and 'nvapbp' are
  !   used in module 'src_obs_cdfout_feedobs' for writing feedobs files;
  !   therefore fill these 3 flags
  mzobhd (nhschr) = 0
  IF ((kobtyp == OT_AIREP) .AND. (fh% phase /= -1)) THEN
    kphase = MAX( 0 , MIN( nibits( nvapoc-1 ) , fh% phase ) )
    mzobhd (nhschr) = ISHFT( kphase , nvapbp )
!   mzobhd (nhschr) = insert( mzobhd(nhschr), kphase, nvapbp )
    !   phase determined from flight tracking in COSMO obs processing
!CS: to do
    IF (fh% tracking == 18)                                                    &
      mzobhd (nhschr) = IBSET( mzobhd(nhschr), nvapbp+nvapoc-1 )
  ENDIF
  IF (fh% sta_corr == 1)  mzobhd (nhschr) = IBSET( mzobhd(nhschr), nvscbp )
  IF (fh% mdlsfc   == 1)  mzobhd (nhschr) = IBSET( mzobhd(nhschr), nvsebp )

  IF (kfmt <= 2) THEN
    zobhed (nhvcbu) = c1
    zobhed (nhvcbt) = c1
    zobhed (nhvcbq) = c1
    zobhed (nhvctu) = c1
    zobhed (nhvctt) = c1
    zobhed (nhvctq) = c1
  ENDIF
  IF ((kfmt == 2) .OR. (kfmt == 0)) THEN
    IF (fh% instype    /= ifill) mzobhd (nhstyp) = fh% instype
  ELSEIF (kfmt == 1) THEN
    IF (fh% instype    /= ifill) THEN
                                 mzobhd (nhrtyp) = fh% instype
      !!!   negative value: temporary indicator that bias correction
      !!!                   must not be (re-)done
      IF (kobtyp == OT_TEMP)     mzobhd (nhrtyp) = - fh% instype
    ENDIF
                                 mzobhd (nhnlev) = nlev
    IF (fh% tracking   /= -1   ) mzobhd (nhtrac) = fh% tracking
    IF (fh% rad_corr   /= -1   ) mzobhd (nhrad ) = fh% rad_corr
    IF (fh% meas_type  /= -1   ) mzobhd (nhna4 ) = fh% meas_type
                                 mzobhd (nhaexi) = imdi
                                 mzobhd (nhuexi) = imdi
                                 mzobhd (nhtexi) = imdi
                                 mzobhd (nhqexi) = imdi
                                 mzobhd (nhwce ) = imdi
                                 mzobhd (nhdt  ) = imdi
!                                zobhed (nhtvip) = rmdi
!                                zobhed (nhtviz) = rmdi
  ENDIF

  IF ((kfmt == 2) .AND. (lassign) .AND. (zobhed(nhalt) > rmdich))              &
    landsy  =        ((kobtyp == OT_SYNOP) .OR. (kobtyp == OT_TEMP))           &
               .AND. (.NOT. BTEST( mzobhd(nhschr),nvsebp ))

!-------------------------------------------------------------------------------
! Section 2: Observation body: fill ODR arrays 'zobbdy', 'mzobbd'
!-------------------------------------------------------------------------------

  !   for convenience: set a local pointer to the fdbk body of current report
  bd  =>  fdbk% b ( fdbk% h(irep)% ib                                          &
                  : fdbk% h(irep)% ib + fdbk% h(irep)% nb - 1 )

  !   set indices and initialise arrays related to (COSMO) ODR body
  IF (kfmt == 1) THEN
    nbxqcf = nbtqcf
    nbxerr = nbterr
    nbxflg = nbtflg
    nbxlid = nbtlid
  ELSEIF ((kfmt == 0) .OR. (kfmt == 2)) THEN
    nbxqcf = nbsqcf
    nbxerr = nbserr
    nbxflg = nbsflg
    nbxlid = nbslid
  ELSEIF (kfmt == 3) THEN
    nbxqcf = nbgqcf
    nbxerr = nbgerr
    nbxflg = nbgflg
    nbxlid = nbglid
  ENDIF
! ALLOCATE ( zsobdy (mxsoxx,nlev) )
  zobbdy  (:,:     ) = rmdi
  mzobbd  (:,:     ) = imdi
  mzobbd  (:,nbxqcf) = 0
  mzobbd  (:,nbxerr) = 0
  mzobbd  (:,nbxflg) = 0
  mzobbd  (:,nbxlid) = 0
  IF (kfmt == 2)  mzobbd  (:,nbsttr) = 0
! zsobdy  (:,:) = rmdi

  IF (kfmt <= 1)  lev_new = .true.
  IF (kfmt == 1)  ilev = 0
  IF (kfmt /= 1)  ilev = 1

! loop_over_obs:  DO iobs = 1 , fdbk% h(irep)% l_body
  loop_over_obs:  DO iobs = 1 , fdbk% h(irep)% nb

    ! get vertical level; determine whether current obs is at a new obs level
    ! -----------------------------------------------------------------------

    !   upper-air obs: determine - pressure level
    !                            - whether current obs is at a new obs level
    IF (kfmt <= 1) THEN
      zzz  =  rmdi
      IF     (bd(iobs)% level_typ == VN_P     ) THEN
        zpp  =  bd(iobs)% level
      ELSEIF (bd(iobs)% level_typ == VN_FLEV  ) THEN
        zpp  =  bd(iobs)% plevel
        zzz  =  bd(iobs)% level
      ELSEIF (bd(iobs)% level_typ == VN_HEIGHT) THEN
        zpp  =  bd(iobs)% plevel
        zzz  =  bd(iobs)% level
      ENDIF
      IF (kfmt == 1) THEN
        lev_new = (ilev == 0)
        IF (ilev >= 1) THEN
          IF (      (bd(iobs)% level_typ == VN_HEIGHT)                         &
              .AND. (zobbdy(ilev,nbtz) > rmdich)) THEN
            lev_new  =  (zzz > zobbdy(ilev,nbtz) + 0.005_wp)
          ELSE
            lev_new  =  (zpp < zobbdy(ilev,nbtp) - 0.5_wp)
          ENDIF
        ENDIF
        IF ((lev_new) .AND. (ilev < maxmlv)) THEN
          ilev = ilev + 1
        ELSEIF (lev_new) THEN
          !   get diagnostic array position
          icma   =  i_cma ( kobtyp , kcdtyp )
                  ! =====
          !   insufficient ODR size (event updated even for passive reports)
          neventd (nelodr,icma) =   neventd(nelodr,icma)                       &
                                  + fdbk% h(irep)% n_level - maxmlv
          ilen = 3 + istrej
          IF (nacout+ilen <= nmxoln) THEN
            outbuf(nacout+1) = ilen
            outbuf(nacout+2) = nfmt4
            outbuf(nacout+3) = fdbk% h(irep)% n_level
            DO icl = 1 , istrej
              outbuf(nacout+3+icl) = ICHAR( yzobhd (icl:icl) )
            ENDDO
            nacout  = nacout + ilen
          ENDIF
  EXIT loop_over_obs
        ENDIF
      ENDIF
      IF (lev_new) THEN
        IF (kfmt == 1)  zobbdy (ilev,nbtp) = zpp
        IF (kfmt == 0)  zobbdy (ilev,nbsp) = zpp
      ENDIF
      IF ((lev_new) .AND. (zzz > rmdich)) THEN
        IF (kfmt == 1)  zobbdy (ilev,nbtz) = zzz
        IF (kfmt == 0)  zobbdy (ilev,nbsz) = zzz
        IF (kfmt == 1)  nbxflg = nbtflg
        IF (kfmt == 0)  nbxflg = nbsflg
        IF (bd(iobs)% level_typ == VN_FLEV  ) THEN
          mzobbd (ilev,nbxflg) = IBSET( mzobbd(ilev,nbxflg), nvfzbp +nvfbps(5) )
        ELSEIF (bd(iobs)% level_typ == VN_HEIGHT) THEN
          mzobbd (ilev,nbxflg) = IBSET( mzobbd(ilev,nbxflg), nvfzbp +nvfbps(6) )
        ENDIF
        !   for upper-air obs, 'zobbdy(:,nbtz)' is filled only further below,
        !   when 'varno'=VN_HEIGHT;  for 'varno'=VN_HEIGHT,
        !   'level_typ' must be VN_P and cannot be VN_HEIGHT or VN_FLEV
      ENDIF
    ENDIF
            
    !   surface-level obs: get height level for 'surface' pressure obs
    !   (for other 'surface' obs, height level is 'station altitude')
    IF ((kfmt >= 2) .AND. (bd(iobs)% varno == VN_PS   ))                       &
      zobbdy (ilev,nbsz) = bd(iobs)% level

    ! set ODR indices for observation variable, obs error, flag bit numbers
    ! ---------------------------------------------------------------------

    nbxx   = 0
    nbxxer = 0
    nbxbia = 0
    nvrx   = -1
    nvfxbp = -1 
    IF     (     (bd(iobs)% varno == VN_U     )                                &
            .OR. (bd(iobs)% varno == VN_U10M  )) THEN
      nbxx   = nbsu
      nbxxer = nbsuer
      nvfxbp = nvfubp
      nvrx   = nvru
!     nso_x  = nso_u
      IF (kfmt == 1)  nbxx   = nbtu
      IF (kfmt == 1)  nbxxer = nbtuer
    ELSEIF (     (bd(iobs)% varno == VN_V     )                                &
            .OR. (bd(iobs)% varno == VN_V10M  )) THEN
      nbxx   = nbsv
      nbxxer = nbsuer
      nvfxbp = nvfubp
      nvrx   = nvru
!     nso_x  = nso_v
      IF (kfmt == 1)  nbxx   = nbtv
      IF (kfmt == 1)  nbxxer = nbtuer
    ELSEIF (     (bd(iobs)% varno == VN_T     )                                &
            .OR. (bd(iobs)% varno == VN_T2M   )) THEN
      nbxx   = nbst
      nbxxer = nbster
      nvfxbp = nvftbp
      nvrx   = nvrt
!     nso_x  = nso_t
      IF (kfmt == 1)  nbxx   = nbtt
      IF (kfmt == 1)  nbxxer = nbtter
      IF (kfmt == 3)  nbxx   = nbgt
      IF (kfmt == 3)  nbxxer = nbgter
    ELSEIF (     (bd(iobs)% varno == VN_RH    )                                &
            .OR. (bd(iobs)% varno == VN_RH2M  )) THEN
      nbxx   = nbsrh
      nbxxer = nbsqer
      nbxbia = nbsdrh
      nvfxbp = nvfqbp
      nvrx   = nvrq
!     nso_x  = nso_rh
      IF (kfmt == 1)  nbxx   = nbtrh
      IF (kfmt == 1)  nbxxer = nbtqer
      IF (kfmt == 1)  nbxbia = nbtdrh
      IF (kfmt == 3)  nbxx   = nbgrh
      IF (kfmt == 3)  nbxxer = nbgqer
      IF (kfmt == 3)  nbxbia = nbgdrh
    ELSEIF (     (bd(iobs)% varno == VN_HEIGHT)                                &
            .OR. (bd(iobs)% varno == VN_PS    )) THEN
      nbxx   = nbsz
      nbxxer = nbszer
      nvfxbp = nvfzbp
      nvrx   = nvrz
!     nso_x  = nso_p
      IF (kfmt == 2)  nbxx   = nbsp
      IF (kfmt == 1)  nbxx   = nbtz
      IF (kfmt == 1)  nbxxer = nbtzer
      IF (kfmt == 3)  nbxx   = nbgp
      IF (kfmt == 3)  nbxxer = nbgzer
    ELSEIF (kfmt == 1) THEN
      IF (bd(iobs)% varno == VN_W     )  nbxx   = nbtw
      nvrx   = nvrw
!   ELSEIF (kfmt == 0) THEN
!     !   not (yet) implemented:  vertical gust (aircraft)
!     IF (bd(iobs)% varno == VN_VGUST )  nbxx   = nbsfgv
    ELSEIF (kfmt == 2) THEN
      IF (bd(iobs)% varno == VN_PTEND )  nbxx   = nbspst
      IF (bd(iobs)% varno == VN_NH    )  nbxx   = nbscbs
      IF (bd(iobs)% varno == VN_N_L   )  nbxx   = nbscl
      IF (bd(iobs)% varno == VN_N_M   )  nbxx   = nbscm
      IF (bd(iobs)% varno == VN_N_H   )  nbxx   = nbsch
      IF (bd(iobs)% varno == VN_N     )  nbxx   = nbsct
      IF (bd(iobs)% varno == VN_NH    )  nvrx   = nvrcbs
      IF (bd(iobs)% varno == VN_N_L   )  nvrx   = nvrcl
      IF (bd(iobs)% varno == VN_N_M   )  nvrx   = nvrcm
      IF (bd(iobs)% varno == VN_N_H   )  nvrx   = nvrch
      IF (bd(iobs)% varno == VN_N     )  nvrx   = nvrct
!     IF (bd(iobs)% varno == VN_NH    )  nso_x  = nso_cbs
!     IF (bd(iobs)% varno == VN_N_L   )  nso_x  = nso_cl
!     IF (bd(iobs)% varno == VN_N_M   )  nso_x  = nso_cm
!     IF (bd(iobs)% varno == VN_N_H   )  nso_x  = nso_ch
!     IF (bd(iobs)% varno == VN_N     )  nso_x  = nso_ct
      IF (bd(iobs)% varno == VN_VV    )  nbxx   = nbsvis
      IF (bd(iobs)% varno == VN_FF    )  nbxx   = nbsff
      IF (bd(iobs)% varno == VN_DD    )  nbxx   = nbsdd
      IF (bd(iobs)% varno == VN_TD2M  )  nbxx   = nbstd
!     IF (bd(iobs)% varno == VN_FF    )  nso_x  = nso_ff
!     IF (bd(iobs)% varno == VN_DD    )  nso_x  = nso_dd
!     IF (bd(iobs)% varno == VN_TD2M  )  nso_x  = nso_td
      IF (bd(iobs)% varno == VN_SDEPTH)  nbxx   = nbshsw
      IF (bd(iobs)% varno == VN_JJ    )  nbxx   = nbstx
      IF (bd(iobs)% varno == VN_TMIN  )  nbxx   = nbstn
      IF (bd(iobs)% varno == VN_GUST  ) THEN
        IF (NINT( bd(iobs)% level ) == 1 )  nbxx = nbsfg1
        IF (NINT( bd(iobs)% level ) == 3 )  nbxx = nbsfg3
        IF (NINT( bd(iobs)% level ) == 6 )  nbxx = nbsfg6
      ENDIF
      IF (bd(iobs)% varno == VN_RR    ) THEN
        IF (NINT( bd(iobs)% level ) == 1 )  nbxx = nbsrr1
        IF (NINT( bd(iobs)% level ) == 3 )  nbxx = nbsrr3
        IF (NINT( bd(iobs)% level ) == 6 )  nbxx = nbsrr6
        IF (NINT( bd(iobs)% level ) == 12)  nbxx = nbsr12
        IF (NINT( bd(iobs)% level ) == 24)  nbxx = nbsr24
      ENDIF
      IF (bd(iobs)% varno == VN_WW    )  nbxx   = - nbswwe
      IF (bd(iobs)% varno == VN_GCLG  )  nbxx   = - nbsclg
      IF (bd(iobs)% varno == VN_ICLG  ) THEN
        IF (NINT( bd(iobs)% level ) == 1 )  nbxx = - nbscl1
        IF (NINT( bd(iobs)% level ) == 2 )  nbxx = - nbscl2
        IF (NINT( bd(iobs)% level ) == 3 )  nbxx = - nbscl3
        IF (NINT( bd(iobs)% level ) == 4 )  nbxx = - nbscl4
      ENDIF
      IF (bd(iobs)% varno == VN_RAD_GL)  nbxx   = nbsrad
      IF (bd(iobs)% varno == VN_RAD_DF)  nbxx   = nbsrdd
      IF (bd(iobs)% varno == VN_RAD_LW)  nbxx   = nbsrdt
    ELSEIF (kfmt == 3) THEN
      IF (bd(iobs)% varno == VN_PWC   ) THEN
        nbxx   = nbgiwa
        nbxbia = nbgbia
        nvfxbp = nvfgbp
        nvrx   = nvriwv
!       nso_x  = nso_iq
      ENDIF
      IF (bd(iobs)% varno == VN_ZPD   ) THEN
        nbxx   = nbgzpd
        nbxxer = nbgtze
        nvfxbp = nvfgbp
        nvrx   = nvrzpd
      ENDIF
      IF (bd(iobs)% varno == VN_ZWD   )  nbxx   = nbgzwd
    ENDIF

    ! fill ODR: observation values, errors, bias correction
    ! -----------------------------------------------------

    IF  (nbxx  == 0)                                                           &
  CYCLE loop_over_obs

    !   real elements
    IF  (nbxx   > 0)         zobbdy (ilev,nbxx  )  =  bd(iobs)% obs
    IF ((nbxxer > 0) .AND. (bd(iobs)% e_o  /= rfill))                          &
                             zobbdy (ilev,nbxxer)  =  bd(iobs)% e_o
    IF ((nbxbia > 0) .AND. (bd(iobs)% bcor /= rfill)) THEN
                             zobbdy (ilev,nbxbia)  =  bd(iobs)% bcor
      IF (nbxbia == nbgbia)  zobbdy (ilev,nbgiwv)  =    zobbdy(ilev,nbgiwa)    &
                                                      - zobbdy(ilev,nbgbia)
    ENDIF

    !   integer elements: weather and cloud groups
    IF  (nbxx   < 0) THEN
      iov = NINT( bd(iobs)% obs )
      IF (nbxx == -nbswwe) THEN
        iov = ISHFT( iov , nvw0bp )
      ELSEIF (nbxx == -nbsclg) THEN
        CALL MVBITS( bd(iobs)% level_sig , 0 , nxsgoc , iov , nxsgbp )
      ELSE
        CALL MVBITS( bd(iobs)% level_sig , 0 , nxsgoc , iov , nxsibp )
      ENDIF
      mzobbd (ilev,-nbxx) = iov
    ENDIF

    ! fill ODR: flags and level significance
    ! --------------------------------------

    IF ((nvrx >= 0) .OR. (nvfxbp >= 0))  kfl = bd(iobs)% flags

    !   QC flag:     do not update, because QC is re-done in current run
    !      IF ((nvrx >= 0).and.(BTEST( kfl,FL_FG ).or.BTEST( kfl,FL_FG_LBC ))) &
    !        mzobbd (ilev,nbxqcf) = IBSET( mzobbd(ilev,nbxqcf) , nvrx )

    !   convert first guess check flags into dataset flag, if required
    IF ((nvrx >= 0) .OR. (nvfxbp >= 0)) THEN
      IF ((losse_fg) .AND. (     (BTEST( kfl, FL_FG     ))                     &
                            .OR. (BTEST( kfl, FL_FG_LBC ))))                   &
         kfl = IBSET( kfl , FL_DATASET )
    ENDIF

    !   pre-processing status flag
    IF (nvrx >= 0) THEN
      !   active if  - status flag in 'fof' is 'active', or
      !              - first guess check was only reason for status 'not active'
      kflfg  =  IBCLR( IBCLR( kfl , FL_FG ) , FL_FG_LBC )
      IF (     ((bd(iobs)% state <= ST_ACTIVE) .AND. (bd(iobs)% state >= 0))   &
          .OR. ((kflfg == 0) .AND. (kfl > 0)))                                 &
        mzobbd (ilev,nbxerr) = IBSET( mzobbd(ilev,nbxerr) , nvrx )
    ENDIF

    !   main flag word
    IF (nvfxbp >= 0) THEN
      IF (BTEST( kfl , FL_DATASET  )) THEN
        CALL MVBITS( 1, 0, nvfboc(1), mzobbd(ilev,nbxflg) , nvfxbp + nvfbps(1) )
!CS>: just for testing, maybe to be cancelled; in feedobs file, it's not possible
!     to distinguish between 'dataset flag not set' and 'dataset flag missing'
      ELSEIF (     ((kobtyp == OT_AIREP) .AND. (     (nvfxbp == nvftbp)        &
                                                .OR. (nvfxbp == nvfqbp)))      &
              .OR. ((kobtyp == OT_PILOT) .AND. (.NOT. lraso)                   &
                                         .AND. (kcdtyp /= OC_RAVAD)            &
                                         .AND. (     (nvfxbp == nvfubp)        &
                                                .OR. (nvfxbp == nvftbp)))) THEN
        CALL MVBITS( 3, 0, nvfboc(1), mzobbd(ilev,nbxflg) , nvfxbp + nvfbps(1) )
      ENDIF
!CS<
      IF (BTEST( kfl , FL_BLACKLIST))                                          &
        mzobbd (ilev,nbxflg) = IBSET( mzobbd(ilev,nbxflg) , nvfxbp + nvfbps(2) )
      IF (BTEST( kfl , FL_PRACTICE ))                                          &
        mzobbd (ilev,nbxflg) = IBSET( mzobbd(ilev,nbxflg) , nvfxbp + nvfbps(5) )
      !   lapse rate and wind shear errors are indicated here as gross error
      IF (BTEST( kfl , FL_GROSS    ))                                          &
        mzobbd (ilev,nbxflg) = IBSET( mzobbd(ilev,nbxflg) , nvfxbp + nvfbps(3) )
      !   (in ODR, if 'nvflbp' is set, then also 'nvfbps(4)' is set for v,T,q,z;
      !    'nvflbp' is not used hereafter and is therefore not set)
      IF (BTEST( kfl , FL_HEIGHT   ))                                          &
        mzobbd (ilev,nbxflg) = IBSET( mzobbd(ilev,nbxflg) , nvfxbp + nvfbps(4) )
      IF (      (BTEST( kfl , FL_NO_OBS )) .AND. (nvfxbp == nvfzbp)            &
          .AND. (kfmt == 1) .AND. (kobtyp /= OT_AIREP))                        &
        mzobbd (ilev,nbxflg) = IBSET( mzobbd(ilev,nbxflg) , nvfxbp + nvfbps(6) )
      IF ((nvfxbp == nvfqbp) .AND. (ABS( bd(iobs)% bcor ) > epsy))             &
        mzobbd (ilev,nbxflg) = IBSET( mzobbd(ilev,nbxflg) , nvfxbp + nvfbps(6) )
      IF (      (nvfxbp == nvfubp) .AND. (kobtyp == OT_AIREP)                  &
          .AND. (bd(iobs)% qual >= 0) .AND. (bd(iobs)% qual <= 3)) THEN
        CALL MVBITS( bd(iobs)% qual, 0, nvfboc(1)                              &
                   , mzobbd(ilev,nbxflg) , nvfxbp + nvfbps(1) )
        CALL MVBITS( bd(iobs)% qual, 0, nvfboc(1), mzobhd(nhschr) , nvaabp )
      ENDIF
    ENDIF

    !   level significance    (and ODR level identity)
    ilevsig = bd(iobs)% level_sig
    IF ((ilevsig == -99) .OR. (ilevsig == ifill))  ilevsig = 0
!CS !     to cope with old convention where 'phase' has been put in 'level_sig':
    IF (kobtyp == OT_AIREP)  ilevsig = 0
    IF (kfmt == 1)  mzobbd (ilev,nbtlsg) = ilevsig
    !
    !   ODR level identity  (really used: SYNOP pressure code, surface bit)
    IF (kobtyp == OT_SYNOP) THEN
      IF (nbxx == nbsp)  mzobbd (ilev,nbslid) = ilevsig
    ELSEIF (kfmt == 3) THEN
      mzobbd (ilev,nbglid) = ISHFT( 1 , nvlidp(7) )
    ELSEIF (     (nvrx == nvru) .OR. (nvrx == nvrt)                            &
            .OR. (nvrx == nvrq) .OR. (nvrx == nvrz)) THEN
      IF (mzobbd(ilev,nbxlid) == imdi)  mzobbd (ilev,nbxlid) = 0
      IF (     (kobtyp == OT_DRIBU) .OR. (kobtyp == OT_SCATT)                  &
          .OR. (BTEST( ilevsig , LS_SURFACE ))) THEN
        mzobbd (ilev,nbxlid) = IBSET( mzobbd(ilev,nbxlid) , nvlidp(7) )
      ELSEIF ((BTEST( ilevsig , LS_SIGN )) .AND. (nvrx == nvru)) THEN
        mzobbd (ilev,nbxlid) = IBSET( mzobbd(ilev,nbxlid) , nvlidp(8) )
      ELSEIF ((BTEST( ilevsig , LS_SIGN )) .AND. (nvrx == nvrt)) THEN
        mzobbd (ilev,nbxlid) = IBSET( mzobbd(ilev,nbxlid) , nvlidp(9) )
      ENDIF
      IF ((kfmt == 1) .AND. (lraso)) THEN
        IF ( BTEST( ilevsig , LS_STANDARD )) THEN
          IF (zobbdy(ilev,nbtp) > 9999._wp) THEN
            mzobbd (ilev,nbxlid) = IBSET( mzobbd(ilev,nbxlid) , nvlidp(6) )
          ELSE
            mzobbd (ilev,nbxlid) = IBSET( mzobbd(ilev,nbxlid) , nvlidp(4) )
          ENDIF
        ENDIF
        IF ( BTEST( ilevsig , LS_TROPO    ))                                   &
          mzobbd (ilev,nbxlid) = IBSET( mzobbd(ilev,nbxlid) , nvlidp(2) )
        IF ( BTEST( ilevsig , LS_MAX      ))                                   &
          mzobbd (ilev,nbxlid) = IBSET( mzobbd(ilev,nbxlid) , nvlidp(1) )
      ENDIF
    ENDIF

    ! fill ODR: other elements in multi-level report
    ! ----------------------------------------------
    IF (kfmt == 1) THEN
      IF (bd(iobs)% qual /= -1)                                                &
        zobbdy (ilev,nbtsnr)  =  REAL( bd(iobs)% qual, wp )
!CS: to do: fill next 4 elements correctly as soon as info is available in fdbk
      zobbdy (ilev,nbtzio)  =  zobhed(nhzio)
      zobbdy (ilev,nbtzjo)  =  zobhed(nhzjo)
      zobbdy (ilev,nbttim)  =  0._wp
!     IF (kobtyp == OT_PILOT) THEN
!       IF (bd(iobs)% accuracy /= rfill)                                       &
!         zobbdy (ilev,nbtuac)  =  bd(iobs)% accuracy
!     ENDIF
    ENDIF

    ! fill ODR: time period of T-max and T-min (surface reports)
    ! ----------------------------------------
    IF (kfmt == 2) THEN
      IF (nbxx == nbstx) mzobbd (ilev,nbsttr) = insert( mzobbd(ilev,nbsttr)    &
                                                      , NINT( bd(iobs)% level) &
                                                      , ntxbp )
      IF (nbxx == nbstn) mzobbd (ilev,nbsttr) = insert( mzobbd(ilev,nbsttr)    &
                                                      , NINT( bd(iobs)% level) &
                                                      , ntnbp )
    ENDIF

    ! re-determine height and status flags for screen-level obs,
    ! if grid point assignment has been (re-)done
    ! ----------------------------------------------------------
    IF ((kfmt == 2) .AND. (lassign) .AND. (zobhed(nhalt) > rmdich)) THEN
      IF (      ((nbxx == nbsu) .OR. (nbxx == nbst) .OR. (nbxx == nbsrh))      &
          .AND. (zobbdy(ilev,nbxx) > rmdich)) THEN
        IF (nbxx == nbsu )  ivar = 1
        IF (nbxx == nbst )  ivar = 3
        IF (nbxx == nbsrh)  ivar = 4
        fisd =  zobhed(nhalt) - zobhed(nhsurf)
        IF (fisd >  epsy)  fisd  =   fisd * fdoro(ivar)
        IF (fisd < -epsy)  fisd  = - fisd
        kfl = 0
        lseaonly = (ABS( altopsu(ivar) ) < epsy)
        IF (     ((.NOT. lseaonly) .AND. (zobhed(nhalt) > altopsu(ivar)))      &
            .OR. ((      lseaonly) .AND. (landsy))                             &
            .OR. (fisd > doromx(ivar)))  kfl = 1
        nflgx = IBITS( mzobbd(ilev,nbsflg), nvfxbp, nvfaoc )
        IF ((.NOT. BTEST( nflgx, nvfbps(4) )) .AND. (kfl == 1)) THEN
          !   if the height flag was not set in the feedobs file and
          !   has to be set now, then also set the passive status flag
          mzobbd (ilev,nbsflg) = IBSET( mzobbd(ilev,nbsflg), nvfxbp +nvfbps(4) )
          mzobbd (ilev,nbserr) = IBSET( mzobbd(ilev,nbserr), nvrx )
        ELSEIF ((BTEST( nflgx, nvfbps(4) )) .AND. (kfl == 0)) THEN
          !   if the height flag set in the feedobs file has to be unset now,
          !   then set the obs status to active only if both the 'dataset flag',
          !   the blacklist, gross error, and bad reporting practice flags),
          !   and the temperature or wind shear checks are all -not- set
          mzobbd (ilev,nbsflg) = IBCLR( mzobbd(ilev,nbsflg), nvfxbp +nvfbps(4) )
          IF (       (IBITS( nflgx, nvfbps(1), nvfboc(1) ) /= 1)               &
              .AND.  (.NOT. BTEST( nflgx, nvfbps(2) ))                         &
              .AND.  (.NOT. BTEST( nflgx, nvfbps(3) ))                         &
              .AND.  (.NOT. BTEST( nflgx, nvfbps(5) ))                         &
              .AND. ((.NOT. BTEST( nflgx, nvfbps(6) )) .OR. (nbxx == nbsrh)))  &
            mzobbd (ilev,nbserr) = IBCLR( mzobbd(ilev,nbserr), nvrx )
        ENDIF
      ENDIF
    ENDIF

    IF (kfmt <= 1)  lev_new = .false.
  ENDDO  loop_over_obs

  ! set aircraft rep. to passive if less than ( 25m or)  3hPa  \  above (model)
  ! set  SATOB   rep. to passive if less than (200m or) 20hPa  /   'orography'
  ! ---------------------------------------------------------------------------
  !   --> do not use AMDAR height because it is ficticious and can be well below
  !                                  the real airport height within anticyclones
  IF ((kfmt == 0) .AND. (lassign)) THEN
    zporo = r_ps(mzobhd(nhio),mzobhd(nhjo)) 
    IF (    ((kobtyp == OT_AIREP) .AND. (zobbdy(ilev,nbsp) > zporo- 300._wp))     &
       .OR. ((kobtyp == OT_SATOB) .AND. (zobbdy(ilev,nbsp) > zporo-2000._wp))) THEN
      mzobhd (nhschr)  =  IBSET ( mzobhd(nhschr), nvalbp )
      mzobhd (nhflag)  =  IBSET ( mzobhd(nhflag), FL_HEIGHT )
      !   get diagnostic array position
      icma   =  i_cma ( kobtyp , kcdtyp )
              ! =====
      IF (mzobhd(nhpass) <= 1) THEN
        neventr (nezdif,icma) = neventr(nezdif,icma) + 1
        cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
        cma(icma)%cnt_ps = cma(icma)%cnt_ps + 1
        mzobhd (nhpass)  =  2
      ENDIF
    ENDIF
  ENDIF

  ! required for nudging purposes only, and related to
  !       (surface) pressure from surface reports only:
  !       -  compute ODR element 'nbsviz' (for weights in nudging)
  !       -  and set invalid height flag  (if grid pt. assignment done)
  ! -------------------------------------------------------------------
  IF (kfmt == 2) THEN
    IF (      (zobbdy(ilev,nbsp) > rmdich) .AND. (zobhed(nhalt) > rmdich)      &
        .AND. (zobbdy(ilev,nbsz) > rmdich)) THEN
      !   compute ODR element 'nbsviz' (as in routine 'obs_cdf_store_singlev'
      !                                 of module 'src_obs_cdfin_sing.f90')
    ! zzkeml = zobhed(nhsurf)         ! approximation, to avoid use of 'r_hhl'
      zzkeml = c05* (  r_hhl(mzobhd(nhio),mzobhd(nhjo),ke+1)                   &
                     + r_hhl(mzobhd(nhio),mzobhd(nhjo),ke  ))
      IF (zobbdy(ilev,nbsz) >= zzkeml-epsy) THEN
        fisdzz = ABS( zobhed(nhalt) - zobbdy(ilev,nbsz) )
        fisdnu = fisdzz + fdoro(2) *(zobbdy(ilev,nbsz) - zzkeml)
      ELSEIF (zobbdy(ilev,nbsz) >= zobhed(nhalt)-epsy) THEN
        fisdzz = zzkeml - zobhed(nhalt)
        fisdnu = fisdzz
      ELSE
        fisdzz = MAX( zobhed(nhalt)  - zobbdy(ilev,nbsz)                       &
                    , zzkeml         - zobbdy(ilev,nbsz) )
        fisdnu = fisdzz
      ENDIF
      zobbdy (ilev,nbsviz)  =  fisdnu
      IF (lassign) THEN
        !   only re-set height flag (to reject it in the nudging);
        !   do not set passive status flag (to allow for using it in the LETKF)
        lseaonly = (ABS( altopsu(2) ) < epsy)
        IF (     ((.NOT. lseaonly) .AND. (zobhed(nhalt) > altopsu(2)))         &
            .OR. ((      lseaonly) .AND. (landsy))                             &
            .OR. (fisdzz > doromx(2))) THEN
          mzobbd (ilev,nbsflg) = IBSET( mzobbd(ilev,nbsflg), nvfzbp +nvfbps(4) )
        ELSE
          mzobbd (ilev,nbsflg) = IBCLR( mzobbd(ilev,nbsflg), nvfzbp +nvfbps(4) )
        ENDIF
      ENDIF
    ENDIF
  ENDIF

  ! discard vertical levels without pressure, and other issues
  ! ----------------------------------------------------------

  IF (kfmt == 1) THEN
    nempty = 0
    !   multi-level body elements
    DO ilev = 1 , nlev
      IF (zobbdy(ilev,nbtp) <= rmdich) THEN
        nempty = nempty + 1
      ELSEIF (nempty > 0) THEN
        zobbdy (ilev-nempty,:) = zobbdy(ilev,:)
        mzobbd (ilev-nempty,:) = mzobbd(ilev,:)
      ENDIF
    ENDDO
    mzobhd (nhnlev)  =  mzobhd(nhnlev) - nempty

    !    fill auxilliary element 'LOG( pressure )'
    DO ilev = 1 , mzobhd(nhnlev)
      zobbdy (ilev,nbtlop)  =  LOG( zobbdy(ilev,nbtp) )
    ENDDO
  ENDIF

  !   scatterometer :   hard code zero observation height
  IF ((kfmt == 2) .AND. (kobtyp == OT_SCATT))                                  &
    zobbdy (ilev,nbsz)  =  0.0_wp

!-------------------------------------------------------------------------------
! Section 3: Check existence of active data in report
!-------------------------------------------------------------------------------

  ! multi-level report: fill ODR header elements for existence of data
  ! ------------------------------------------------------------------

  IF (kfmt == 1) THEN
    mexi = 0
    DO ilev = 1 , mzobhd(nhnlev)
      IF (zobbdy(ilev,nbtu ) > rmdich)  mexi (1) = -1
      IF (zobbdy(ilev,nbtt ) > rmdich)  mexi (2) = -1
      IF (zobbdy(ilev,nbtrh) > rmdich)  mexi (3) = -1
    ENDDO
    DO ilev = 1 , mzobhd(nhnlev)
      IF (BTEST( mzobbd(ilev,nbxerr), nvru ))  mexi (1) = 1
      IF (BTEST( mzobbd(ilev,nbxerr), nvrt ))  mexi (2) = 1
      IF (BTEST( mzobbd(ilev,nbxerr), nvrq ))  mexi (3) = 1
    ENDDO
    mzobhd (nhuexi)  =  mexi(1)
    mzobhd (nhtexi)  =  mexi(2)
    mzobhd (nhqexi)  =  mexi(3)
    mzobhd (nhaexi)  =  MAXVAL( mexi )
    IF (mzobhd(nhaexi) == 0)  mzobhd (nhaexi)  =  MINVAL( mexi )
  ENDIF

  ! check existence of active data and update event counters
  ! --------------------------------------------------------

  IF (kfmt == 1)  lexact  =  (mzobhd(  nhaexi) == 1)
  IF (kfmt /= 1)  lexact  =  (mzobbd(1,nbxerr) >= 1)

  IF (.NOT. lexact) THEN
    !   get diagnostic array position
    icma   =  i_cma ( kobtyp , kcdtyp )
            ! =====
    !   (merged reports are counted as active)
    iactr  =  1  -  mzobhd(nhpass) /2
    neventr (nenoda,icma) = neventr(nenoda,icma) + iactr

    IF ((iactr == 1) .AND. (lverpas)) THEN
      !   set report to passive in report header
      mzobhd (nhpass)  =  2
      mzobhd (nhflag)  =  IBSET( mzobhd(nhflag) , FL_NO_OBS )
      cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
      cma(icma)%cnt_ps = cma(icma)%cnt_ps + 1

    ELSEIF (.NOT. lverpas) THEN
      !   flag reports which are to be discarded completely
      mzobhd (nhpass)  = -1
      IF (iactr == 1) cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
      IF (iactr == 0) cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
      cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! End Subroutine obs_fof_store_1rep
!-------------------------------------------------------------------------------
END SUBROUTINE obs_fof_store_1rep


!===============================================================================
!+ Module procedure in "src_obs_fdbk_in" for reading and distributing reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_fof_mask_reports ( n_hdr, yfofpath , fdbk, imindif              &
                                , l_mask , lassign )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_fdbk_in" updates the mask to
!   (either further process or) discard reports in 'fdbk', and determines the
!   report flags. 
!   Reports can be discarded here for the following reasons:
!     - observation / code type unknown to COSMO
!     - station altitude of surface report missing
!     - station location missing or outside model domain.
!   Reports can be flagged to be set passive (or be rejected) for the following
!   reasons:
!     - ship / buoy report assigned to land grid point
!     - report outside specified observation area
!     - station height not in valid range
!     - bad aircraft id
!     - report in exclusion area for observation type
!   Other report flags in 'fdbk' are also prepared to be handed to the ODR.
!   Furthermore, the horizontal assignment of obs reports to a model grid point
!   is done, if the grid point indices in the fdbk report header are missing
!   or if the horizontal model configuration in the fdbk file does not coincide
!   with the one of the current COSMO run. In this case, these indices are
!   updated in 'fdbk'.
!
! Method:
!   Straightforward checks.
!   Processing and report event counters are also updated, and caution messages
!   and messages on discarded reports are written if required.
!
! Written by        : DWD, Christoph Schraff, 19.12.2012
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER (KIND=iintegers) , INTENT (IN)     ::  &
    n_hdr            ,& ! number of reports (header elements) in 'fdbk'
                        !   (n_hdr == fdbk% f% n_hdr  must be satisfied !)
    imindif             ! difference in minutes between fdbk  reference time
                        !                           and model reference time

! CHARACTER (LEN= icdfinlen + iannexlen)     ::  &
  CHARACTER (LEN= *)       , INTENT (IN)     ::  &
    yfofpath            ! full path of input feedobs (feedback) file

  TYPE    (t_fdbk_data)    , INTENT (INOUT)  ::  &
    fdbk                ! feedback file data structure

  LOGICAL                  , INTENT (INOUT)  ::  &
    l_mask (n_hdr)      ! mask: if true  then report is processed further
                        !       if false then report is discarded 

  LOGICAL                  , INTENT (OUT)    ::  &
    lassign             ! compute assignment of obs report to model grid point

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    kproc  (n_hdr)   ,& ! flags which indicate no further processing
    kflag  (n_hdr)   ,& ! flags which indicate rejection or setting to passive
    icmaa  (n_hdr)   ,& ! pointer index for 'cma' / diagnostic arrays
    icma             ,& ! pointer index for 'cma' / diagnostic arrays
    irep             ,& ! record index over reports
    kobtyp , kcdtyp  ,& ! observation type , observation code type
    nerrobt          ,& ! number of reports with unknown data (sub-)category
    ireperr          ,& ! index of last report with unknown data (sub-)category
    nelevt           ,& ! station altitude
    iob_tot, job_tot ,& ! grid point assigned to obs (total model area)
    ffl              ,& ! fdbk report flags 'r_flags'
    istat               ! error indicators

  LOGICAL                  ::  &
    lsfcob           ,& ! orographic surface report conditions are met
    lupair           ,& ! upper-air report  
    lseaobs          ,& ! sea observation type
    lseaob2             ! sea observation type

  REAL (KIND=wp)           ::  &
    zsurf            ,& ! height of model grid pt. to which obs. is assigned
    zio_tot, zjo_tot ,& ! observation location in grid point units (total area)
    roblat , roblon  ,& ! geographical latitude and longitude
    rlat   , rlon    ,& ! latitude and longitude in the rotated system
    zrephr              ! obs time [hours rel. to model reference time]

  LOGICAL                  ::  &
    lnewxy (n_hdr)      ! grid point assignment is (re-)done

  CHARACTER (LEN=ichar10)  :: &
    ystid               ! station identity (char length as in type t_fdbk_head)
  CHARACTER (LEN=25)       :: &
    yroutine            ! name of this subroutine
  CHARACTER (LEN=30)       :: &
    yerr                ! error message
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_fof_mask_reports
!-------------------------------------------------------------------------------

  yroutine = 'obs_fof_mask_reports'

  kproc  (:) = 0
  kflag  (:) = 0
  lnewxy (:) = .FALSE.

  ! check if observation type is known to COSMO
  ! -------------------------------------------
  icmaa (:) = 0
  nerrobt   = 0
  DO irep = 1 , fdbk% f% n_hdr
    IF (l_mask(irep)) THEN
      kobtyp = fdbk% h(irep)% obstype
      kcdtyp = fdbk% h(irep)% codetype
      !   get diagnostic array position
      icmaa (irep)  =  i_cma ( kobtyp , kcdtyp )
                     ! =====
      IF (icmaa(irep) == 0) THEN
        kproc (irep) = kproc(irep) + 2
        nerrobt      = nerrobt + 1 
        IF (nerrobt == 1)  ireperr = irep
      ENDIF
    ENDIF
  ENDDO
  !    write alerting CAUTION messages onto file YUCAUTN if obs types unknown
  IF (nerrobt > 0) THEN
    OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'            &
                               , POSITION='APPEND', IOSTAT=istat)
    IF (istat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
    IF (istat /= 0) CALL model_abort (my_cart_id, 1509, yerr, yroutine)
    WRITE( nucautn,'("CAUTION !!!!! t=",F6.3,":",I5," TIMES WRONG OBS"         &
                   &," (CODE) TYPE IN FILE:",A)' )  acthr, nerrobt, yfofpath
    WRITE( nucautn,'("CAUTION !!!!!   ",6X,"  E.G. OBS TYPE:",I11              &
                   &,", CODE TYPE:",I11)')  fdbk% h(ireperr)% obstype          &
                                          , fdbk% h(ireperr)% codetype
    WRITE( nucautn,'("        ==>  CHECK   NetCDF - FILE ",A)' ) yfofpath
    CLOSE (nucautn)
  ENDIF

  ! check if station altitude is missing where required
  ! ---------------------------------------------------
  DO irep = 1 , fdbk% f% n_hdr
    IF ((l_mask(irep)) .AND. (kproc(irep) == 0)) THEN
      IF (fdbk% h(irep)% z_station == ifill) THEN
        kobtyp = fdbk% h(irep)% obstype
        IF (     (kobtyp == OT_SYNOP) .OR. (kobtyp == OT_GPSGB)                &
            .OR. (kobtyp == OT_DRIBU) .OR. (kobtyp == OT_SCATT))               &
          kproc (irep)  =  kproc(irep) + 32
      ENDIF
    ENDIF
  ENDDO

  !   check if horiz. model geometry in feedobs file coincides with current one
  lassign =      (ABS( fdbk% f% lower_left(1) - r_startlat_tot ) > epsy)       &
            .OR. (ABS( fdbk% f% lower_left(2) - r_startlon_tot ) > epsy)       &
            .OR. (ABS( fdbk% f% pole      (1) - r_pollat       ) > epsy)       &
            .OR. (ABS( fdbk% f% pole      (2) - r_pollon       ) > epsy)       &
            .OR. (ABS( fdbk% f% resolution(1) - r_dlat         ) > epsy)       &
            .OR. (ABS( fdbk% f% resolution(2) - r_dlon         ) > epsy)       &
            .OR. (ABS( fdbk% f% domain    (1) - ie_tot         ) > 0   )       &
            .OR. (ABS( fdbk% f% domain    (2) - je_tot         ) > 0   )

  ! check if station location is within model domain, and perform
  ! horizontal assignment of reports to model grid points, if required
  ! ------------------------------------------------------------------
  DO irep = 1 , fdbk% f% n_hdr
    IF ((l_mask(irep)) .AND. (kproc(irep) == 0)) THEN
      roblat  =  fdbk% h(irep)% lat
      roblon  =  fdbk% h(irep)% lon
      IF (      (ABS( roblat ) <=   90._wp)                                    &
          .AND. (     roblon   <=  360._wp +epsy)                              &
          .AND. (     roblon   >= -180._wp -epsy)) THEN
        !   check if horizontal assignment has to be computed
        IF ((lassign) .OR. (fdbk% h(irep)% index_x   <  1)                     &
                      .OR. (fdbk% h(irep)% index_y   <  1)                     &
                      .OR. (fdbk% h(irep)% index_x   >  ie_tot)                &
                      .OR. (fdbk% h(irep)% index_y   >  je_tot)                &
                      .OR. (fdbk% h(irep)% mdlsfc    == ifill)                 &
                      .OR. (fdbk% h(irep)% z_modsurf == ifill)) THEN
          !   convert into rotated coordinates
          rlon = rla2rlarot (roblat, roblon, r_pollat, r_pollon, r_polgam)
          rlat = phi2phirot (roblat, roblon, r_pollat, r_pollon)
          !   observation position in grid point units (total area)
          zio_tot    = 1._wp + (rlon - r_startlon_tot) /r_dlon
          zjo_tot    = 1._wp + (rlat - r_startlat_tot) /r_dlat
          !   assign observation to grid point
          nelevt = imdi
          IF (fdbk% h(irep)% z_station /= ifill)                               &
            nelevt = fdbk% h(irep)% z_station
          lupair  =       (kobtyp /= OT_SYNOP) .AND. (kobtyp /= OT_GPSGB)      &
                    .AND. (kobtyp /= OT_DRIBU) .AND. (kobtyp /= OT_SCATT)
          lseaobs =      (kobtyp == OT_DRIBU) .OR. (kobtyp == OT_SCATT)        &
                    .OR. (kcdtyp == OC_AHSCD) .OR. (kcdtyp == OC_ATSHS)
!                   .OR. (kcdtyp == OC_ABSCD) .OR. (kcdtyp == OC_SHRED)
          lseaob2 = lseaobs

          CALL obs_assign_gridpt ( zio_tot, zjo_tot, nelevt, lupair, .FALSE.   &
                                 , ie_tot, je_tot, hsurf_tot, fland_tot        &
                                 , r_dlat, r_dlon, r_startlat_tot, r_degrad    &
                                 , doromx, fdoro, imdi, nboundlines , nolbc    &
                                 , lseaobs , iob_tot, job_tot, zsurf, lsfcob )
        ! ======================
          !   (lost info: 'lseaobs /= lseaob2', 'iob_tot < 0', lsfcob)
!                           isurfob (irep) = 0
!         IF (lsfcob)       isurfob (irep) = 1
!         IF (lseaobs)      mrepflg (irep) = IBSET( mrepflg(irep) , nvsebp )
          !   reject report if observation location outside model domain
          IF   (iob_tot == 0) THEN
            kproc (irep)  =  kproc(irep) + 16
          ELSE
            lnewxy (irep) = .TRUE.
            fdbk% h(irep)% index_x    =  ABS( iob_tot )
            fdbk% h(irep)% index_y    =  ABS( job_tot )
            fdbk% h(irep)% z_modsurf  =  NINT( zsurf )
            fdbk% h(irep)% mdlsfc     =  0
            IF (lseaobs)  fdbk% h(irep)% mdlsfc  =  1
            !   assigned surface type (land / sea) not appropriate
            IF ((lseaob2) .AND. (.NOT. lseaobs))                               &
              kflag (irep) = IBSET( kflag(irep) , FL_SURF )
            !   distance "model orography - station altitude"
            !   too large for surface obs
            IF (iob_tot < 0)  kflag (irep) = IBSET( kflag(irep) , FL_HEIGHT )
          ENDIF
        ENDIF
      ELSE
        kproc (irep)  =  kproc(irep) + 8
        OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'        &
                                   , POSITION='APPEND', IOSTAT=istat)
        IF (istat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
        IF (istat /= 0) CALL model_abort (my_cart_id, 1510, yerr, yroutine)
        zrephr = (fdbk% h(irep)% time + imindif) /60._wp
        WRITE( nucautn,'("CAUTION !!!!! t=",F6.3,", FILE:",A,": OBS LOCATION " &
                       &,"UNDEFINED, STA:",A,F6.1," HRS")' )                   &
               acthr, yfofpath, fdbk% h(irep)% statid, zrephr
        CLOSE (nucautn)
      ENDIF
    ENDIF
  ENDDO

  ! second, vectorisable loop over reports to be processed
  ! ------------------------------------------------------
  DO irep = 1 , fdbk% f% n_hdr
    IF ((l_mask(irep)) .AND. (kproc(irep) == 0)) THEN
      roblat = fdbk% h(irep)% lat
      roblon = fdbk% h(irep)% lon
      icma   = icmaa (irep)

      ! observation location outside user-defined area
      ! ----------------------------------------------
      IF (     (roblat > cma(icma)%obnlat) .OR. (roblat < cma(icma)%obslat)    &
          .OR. (      (cma(icma)%obwlon >= cma(icma)%obelon)                   &
                .AND. (roblon < cma(icma)%obwlon)                              &
                .AND. (roblon > cma(icma)%obelon))                             &
          .OR. (      (cma(icma)%obwlon <  cma(icma)%obelon)                   &
                .AND. (     (roblon < cma(icma)%obwlon)                        &
                       .OR. (roblon > cma(icma)%obelon))))                     & 
        kflag (irep) = IBSET( kflag(irep) , FL_AREA )

      ! check for bad aircraft station id 
      ! ---------------------------------
      IF (fdbk% h(irep)% obstype == OT_AIREP) THEN
        IF (     (fdbk% h(irep)% statid (1:3) == "XXX")                        &
            .OR. (fdbk% h(irep)% statid (1:3) == "???")                        &
            .OR. (fdbk% h(irep)% statid (1:3) == "///")                        &
            .OR. (fdbk% h(irep)% statid (1:3) == "***")                        &
            .OR. (fdbk% h(irep)% statid (1:3) == "   "))                       &
          kflag (irep) = IBSET( kflag(irep) , FL_BLACKLIST )
      ENDIF

      ! check whether in exclusion area for current observation or code type
      ! --------------------------------------------------------------------
      IF (      (roblat <= cma(icma)%exnlat+epsy)                              &
          .AND. (roblat >= cma(icma)%exslat-epsy)                              &
          .AND. (roblon >= cma(icma)%exwlon-epsy)                              &
          .AND. (roblon <= cma(icma)%exelon+epsy))                             &
        kflag   (irep) = IBSET( kflag  (irep) , FL_OBSTYPE )
    ENDIF
  ENDDO

  ! update report event counters and write control messages
  ! for reports to be discarded
  ! -------------------------------------------------------
  DO irep = 1 , fdbk% f% n_hdr
    IF ((l_mask(irep)) .AND. (kproc(irep) > 0)) THEN
      IF (.NOT. lopen_rej) THEN
        OPEN (nurej ,FILE=yurejct ,FORM='FORMATTED',STATUS='UNKNOWN'           &
                                  ,POSITION='APPEND',IOSTAT=istat)
        lopen_rej = .TRUE.
        IF (istat /= 0)  yerr = 'OPENING OF FILE yurejct FAILED'
        IF (istat /= 0)  CALL model_abort (my_cart_id, 7004, yerr, yroutine)
      ENDIF
      zrephr = (fdbk% h(irep)% time + imindif) /60._wp
      IF     (MOD(kproc(irep),  4) >=  2) THEN
        WRITE( nurej,'(" STATION ",A ," : UKNOWN OBSERVATION / CODE TYPES",    &
                      &2I11,",",F6.1," HRS")' )   fdbk% h(irep)% statid        &
             , fdbk% h(irep)% obstype, fdbk% h(irep)% codetype, zrephr
      ELSEIF (MOD(kproc(irep), 64) >= 32) THEN
        WRITE( nurej,'(" STATION ",A ," : SYNOP STA. HEIGHT MISSING",          &
                      &F6.1," HRS")' )   fdbk% h(irep)% statid , zrephr
        IF (fdbk% h(irep)% r_state <= ST_ACTIVE)                               &
          neventr (nenoal,icmaa(irep)) = neventr(nenoal,icmaa(irep)) + 1
      ELSEIF (MOD(kproc(irep), 32) >= 16) THEN
        WRITE( nurej,'(" STATION ",A ," : OBS. LOCATION OUT OF DOMAIN "        &
                     &,2F7.1," , ",F6.1," HRS")' )   fdbk% h(irep)% statid     &
             , fdbk% h(irep)% lon, fdbk% h(irep)% lat, zrephr
        IF (fdbk% h(irep)% r_state <= ST_ACTIVE)                               &
          neventr (neloca,icmaa(irep)) = neventr(neloca,icmaa(irep)) + 1
      ELSEIF (MOD(kproc(irep), 16) >=  8) THEN
        IF (fdbk% h(irep)% r_state <= ST_ACTIVE)                               &
          neventr (neloca,icmaa(irep)) = neventr(neloca,icmaa(irep)) + 1
      ENDIF
    ENDIF
  ENDDO

  ! write control messages for reports to be rejected
  ! -------------------------------------------------
  DO irep = 1 , fdbk% f% n_hdr
    IF ((l_mask(irep)) .AND. (kproc(irep) == 0) .AND. (kflag(irep) > 0)) THEN
      IF (.NOT. lopen_rej) THEN
        OPEN (nurej ,FILE=yurejct ,FORM='FORMATTED',STATUS='UNKNOWN'           &
                                  ,POSITION='APPEND',IOSTAT=istat)
        lopen_rej = .TRUE.
        IF (istat /= 0)  yerr = 'OPENING OF FILE yurejct FAILED'
        IF (istat /= 0)  CALL model_abort (my_cart_id, 7004, yerr, yroutine)
      ENDIF
      zrephr = (fdbk% h(irep)% time + imindif) /60._wp
      ystid  =  fdbk% h(irep)% statid
      IF     (BTEST( kflag(irep),FL_SURF      )) THEN
        WRITE( nurej,'(" STATION ",A ," : SEA OBS. IS LOCATED OVER (MODEL) LA" &
                     &,"ND AT ",2F7.1," , ",F6.1," HRS")' )                    &
               ystid, fdbk% h(irep)% lon, fdbk% h(irep)% lat, zrephr
      ELSEIF (BTEST( kflag(irep),FL_AREA      )) THEN
        WRITE( nurej,'(" STATION ",A ," : OBS. LOCATION OUT OF USER-SPECIFIED" &
                     &," AREA ",2F7.1," , ",F6.1," HRS")' )                    &
               ystid, fdbk% h(irep)% lon, fdbk% h(irep)% lat, zrephr
      ELSEIF (BTEST( kflag(irep),FL_HEIGHT    )) THEN
        WRITE( nurej,'(" STATION ",A ," : HEIGHT ", I6 ," DIFF. TO MODEL "     &
                     &, I6 ," TOO LARGE, ",F7.1," HRS")')                      &
               ystid, fdbk% h(irep)% z_station, fdbk% h(irep)% z_modsurf, zrephr
      ELSEIF (BTEST( kflag(irep),FL_BLACKLIST )) THEN
        WRITE( nurej,'(" STATION ",A ," : BAD AIRCRAFT ID")')  ystid
      ENDIF
    ENDIF
  ENDDO

  ! update fdbk report flags and simplified report status
  !   (note: report 'r_check' is not needed later on)
  ! -----------------------------------------------------
  DO irep = 1 , fdbk% f% n_hdr
    IF ((l_mask(irep)) .AND. (kproc(irep) == 0)) THEN
      ffl = fdbk% h(irep)% r_flags
      IF (lnewxy(irep)) CALL MVBITS( kflag(irep), FL_SURF  , 1, ffl, FL_SURF   )
      IF (lnewxy(irep)) CALL MVBITS( kflag(irep), FL_HEIGHT, 1, ffl, FL_HEIGHT )
      CALL MVBITS ( kflag(irep) , FL_AREA      , 1 , ffl , FL_AREA      )
      CALL MVBITS ( kflag(irep) , FL_BLACKLIST , 1 , ffl , FL_BLACKLIST )
      CALL MVBITS ( kflag(irep) , FL_OBSTYPE   , 1 , ffl , FL_OBSTYPE   )
      fdbk% h(irep)% r_flags = ffl

      !   simplified report status: 'pas_rej', 'rejected' are set to 'passive'
                                       fdbk% h(irep)% r_state = ST_ACTIVE
      IF (BTEST( ffl, FL_MERGE )    )  fdbk% h(irep)% r_state = ST_MERGED
      IF (IBCLR( ffl, FL_MERGE ) > 0)  fdbk% h(irep)% r_state = ST_PASSIVE
    ENDIF
  ENDDO

  ! update flags and control messages for reports to be rejected or set passive
  ! ---------------------------------------------------------------------------
  DO irep = 1 , fdbk% f% n_hdr
    IF ((l_mask(irep)) .AND. (kproc(irep) == 0)) THEN
      !   note that the order of setting flags and hence of report events
      !   is not exactly as indicated in 'flags' in 'mo_fdbk_tables.f90'
      !   for flags 'FL_AREA', 'FL_BLACKLIST', 'FL_HEIGHT', and 'FL_SURF'.
      IF     (BTEST( fdbk% h(irep)% r_flags, FL_SUSP_LOCT  )) THEN
        neventr (nedbfl,icmaa(irep)) = neventr(nedbfl,icmaa(irep)) + 1
      ELSEIF (BTEST( fdbk% h(irep)% r_flags, FL_SURF       )) THEN
        neventr (nezdif,icmaa(irep)) = neventr(nezdif,icmaa(irep)) + 1
      ELSEIF (BTEST( fdbk% h(irep)% r_flags, FL_AREA       )) THEN
        neventr (neloca,icmaa(irep)) = neventr(neloca,icmaa(irep)) + 1
      ELSEIF (BTEST( fdbk% h(irep)% r_flags, FL_HEIGHT     )) THEN
        neventr (nezdif,icmaa(irep)) = neventr(nezdif,icmaa(irep)) + 1
      ELSEIF (BTEST( fdbk% h(irep)% r_flags, FL_BLACKLIST  )) THEN
        neventr (neblak,icmaa(irep)) = neventr(neblak,icmaa(irep)) + 1
      ELSEIF (BTEST( fdbk% h(irep)% r_flags, FL_OBSTYPE    )) THEN
        neventr (neobct,icmaa(irep)) = neventr(neobct,icmaa(irep)) + 1
!     ELSEIF (.NOT. lreproc) THEN
!       IF     (BTEST( fdbk% h(irep)% r_flags, FL_REDUNDANT  )) THEN
        ELSEIF (BTEST( fdbk% h(irep)% r_flags, FL_REDUNDANT  )) THEN
          neventr (neredn,icmaa(irep)) = neventr(neredn,icmaa(irep)) + 1
        ELSEIF (BTEST( fdbk% h(irep)% r_flags, FL_FLIGHTTRACK)) THEN
          neventr (netrac,icmaa(irep)) = neventr(netrac,icmaa(irep)) + 1
        ELSEIF (BTEST( fdbk% h(irep)% r_flags, FL_THIN       )) THEN
          neventr (nethin,icmaa(irep)) = neventr(nethin,icmaa(irep)) + 1
!       ENDIF
      ENDIF
    ENDIF
  ENDDO

  ! update counters for processed/ active/ passive (rejected)/ discarded reports
  ! ----------------------------------------------------------------------------
  DO irep = 1 , fdbk% f% n_hdr
    IF (l_mask(irep)) THEN
      icma   = icmaa (irep)
      cma(icma)%cnt_pr = cma(icma)%cnt_pr + 1
      IF (kproc(irep) /= 0) THEN 
        !   means: report outside model domain, or synop without station height
        !   ==> report is always thrown away (not only set passive)
        cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
      ELSEIF (fdbk% h(irep)% r_state >= ST_PASSIVE) THEN
        IF (      lverpas)  cma(icma)%cnt_ps = cma(icma)%cnt_ps + 1
        IF (.NOT. lverpas)  cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
      ELSE
        !   (merged reports are counted as active)
        cma(icma)%cnt_ac = cma(icma)%cnt_ac + 1
      ENDIF
    ENDIF
  ENDDO

  ! update mask for further processing
  ! ----------------------------------
  DO irep = 1 , fdbk% f% n_hdr
    IF (     (kproc(irep) >= 1)                                                &
        .OR. ((fdbk% h(irep)% r_state >= ST_PASSIVE) .AND. (.NOT. lverpas)))   &
      l_mask (irep)  =  .FALSE.
  ENDDO

!-------------------------------------------------------------------------------
! End Subroutine obs_fof_mask_reports
!-------------------------------------------------------------------------------
END SUBROUTINE obs_fof_mask_reports


!===============================================================================
!+ Module procedure in "src_obs_fdbk_in" for perturbing observation values
!-------------------------------------------------------------------------------

SUBROUTINE obs_fof_perturb  ( n_bdy, fperturb, iseed_ex , fb )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure adds perturbations to the observation values.
!   Each perturbation is made up of the observation error as specified in the
!   feedback file observation body, a scaling factor 'fperturb', and a Gaussian
!   random number (with variance = 1).
!   For observation quantities with upper and/or lower limits, the perturbations
!   are reduced to prevent the perturbed value from exceeding this limits. This
!   reduction is done in a way to avoid the introduction of a bias.
!
! Method:
!   In the first call, the Gaussian random number generator is initialized.
!   The Gaussian random numbers are then generated by using module 'mo_random',
!   i.e. using Donald E. Knuth's portable pseudo-random number generator,
!   from Seminumerical Algorithms by D.E. Knuth, 3rd edition (1997), including
!   the modifications made in the 9th printing (2002) (Double precision version
!   only).
!   Note: The regular random number generator with flat distribution is strictly
!         machine-independent. However, in the conversion from flat to Gaussian
!         distribution, machine-dependent rounding errors in a check may rarely
!         lead to a different result, so that long sequences of Gaussian random
!         numbers are often identical in the beginning but may differ (strongly)
!         after some point on different machines. 
!
!   To avoid the introduction of a bias, the reduction of a perturbation for an
!   observation quantity with an upper and/or lower limit takes into account
!   the modulus of the original perturbation. Thus, for an observation near an
!   upper limit, both a large positive or a large negative perturbation would be
!   reduced to the same size.
!
!   Known shortcoming: IWV (PWC) and ZPD are perturbed independently from each
!                      other, even though the obs values are related to each
!                      other.
!                      In COSMO, the reported bias corrected IWV obs value (ODR
!                      entry 'nbgiwa') is replaced later on in the obs operator
!                      by the value derived from observed ZPD. As a result, the
!                      perturbations for bias-corrected IWV and ZPD become
!                      consistent, then. However, the uncorrected IWV obs value
!                      (ODR entry 'nbgiwv') is not replaced later on and remains
!                      therefore inconsistent with 'nbgiwa' and ZPD entries.
!                      While this entry is written only to the VOF, the bias
!                      correction itself (entry 'nbgbia'), which is updated in
!                      the obs operator and therefore becomes meaningless in
!                      this case, is written also to the NetCDF feedobs file!
!
! Written by        : DWD, Christoph Schraff, 20.03.2013
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER (KIND=iintegers) , INTENT (IN)     ::  &
    n_bdy     ,& ! number of observations (body elements) in 'fdbk'
                 !   (n_bdy == SIZE( fdbk% b ) must be satisfied !)
    iseed_ex     ! external seed number to initialize random generator

  REAL    (KIND=wp)        , INTENT (IN)     ::  &
    fperturb     ! factor to obs error variance with which the obs are perturbed

  TYPE    (t_fdbk_body)    , INTENT (INOUT)  ::  &
    fb (n_bdy)   ! processed part of feedback file body

! Local parameters: None
! ----------------

  REAL    (KIND=wp)        , PARAMETER       ::  &
    clarge = 1.E9_wp ,& ! large value which is (much) larger than the modulus
                        !   of any observation value of any observation variable
    clchk  = .9E9_wp    ! check value

! Local scalars:
! -------------

  TYPE (random_state_t) , SAVE ::   rg_state  ! state of random number generator

  LOGICAL               , SAVE ::  & 
    lfirst = .TRUE.     ! true only if routine is called for the first time

  !   for random number generator
  REAL    (KIND=dp)            ::  &
    dp_random_gauss (n_bdy)  ! random numbers

  REAL    (KIND=wp)        ::  &
    zpert          ,& ! regular perturbation value
    zapert         ,& ! modulus of perturbation
    zlim_upp       ,& ! upper limt for observed quantity
    zlim_low          ! lower limt for observed quantity

  INTEGER (KIND=iintegers) ::  &
    iobs           ,& ! loop index over observations
    iseed             ! seed number to initialize random generator

  LOGICAL                  ::  &
    lpert             ! obs value still needs to be perturbed

! CHARACTER (LEN=ichar10)  :: &
!   ystid             ! station identity (char length as in type t_fdbk_head)
  CHARACTER (LEN=25)       :: &
    yroutine          ! name of this subroutine
  CHARACTER (LEN=40)       :: &
    yerr              ! error message
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_fof_perturb
!-------------------------------------------------------------------------------

  yroutine = 'obs_fof_perturb'

  !  initialize Gaussian random number generator
  !  -------------------------------------------

  iseed  =  iseed_ex

! IF (lfirst)  CALL construct_seed ( iseed )              ! construct total seed
             ! ===================
  IF (lfirst)  CALL set_seed_rand_numb ( ydate_ref, 0 , iseed )
             ! =======================

  IF (lfirst)  CALL construct ( rg_state, seed=iseed )    ! initialize sequence
             ! ==============

  lfirst = .FALSE.

  !  produce required array of Gaussian random numbers
  !  -------------------------------------------------

  IF (n_bdy /= SIZE( fb ))  yerr = 'fdbk body size does not match n_bdy'
  IF (n_bdy /= SIZE( fb ))  CALL model_abort (my_cart_id, 1511, yerr, yroutine)

  CALL random_gauss ( dp_random_gauss , rg_state )
! =================

  DO iobs = 1 , n_bdy
    IF (      (fb(iobs)% obs /= rfill)                                         &
        .AND. (fb(iobs)% e_o /= rfill)) THEN

      !  set upper and lower limits for observed quantities
      !  --------------------------------------------------

      zlim_upp =   clarge
      zlim_low = - clarge
      IF (     (fb(iobs)% varno == VN_RH    )                                  &
          .OR. (fb(iobs)% varno == VN_RH2M  )) THEN
        zlim_upp = 1._wp  
        zlim_low = .005_wp
      ELSEIF (     (fb(iobs)% varno == VN_N_L )                                &
              .OR. (fb(iobs)% varno == VN_N_M )                                &
              .OR. (fb(iobs)% varno == VN_N_H )                                &
              .OR. (fb(iobs)% varno == VN_N   )) THEN
        zlim_upp = 1._wp  
        zlim_low = 0._wp
      ELSEIF (     (fb(iobs)% varno == VN_NH    )                              &
              .OR. (fb(iobs)% varno == VN_VV    )                              &
              .OR. (fb(iobs)% varno == VN_SDEPTH)                              &
              .OR. (fb(iobs)% varno == VN_GUST  )                              &
              .OR. (fb(iobs)% varno == VN_RR    )                              &
              .OR. (fb(iobs)% varno == VN_PWC   )                              &
              .OR. (fb(iobs)% varno == VN_ZWD   )) THEN
        zlim_low = 0._wp
      ENDIF

      !  compute regular perturbation value
      !  ----------------------------------

      lpert   =  .TRUE.
      !   regular perturbation value
      zpert   =  fb(iobs)% e_o  * REAL(dp_random_gauss(iobs),wp)  *fperturb
      !   modulus of perturbation
      zapert  =  ABS( zpert )
!     zapert  =  fb(iobs)% e_o *ABS( REAL(dp_random_gauss(iobs),wp)) *fperturb

      !  if observed quantities have an upper and/or lower limit, and if
      !  the obs value plus or minus the perturbation exceeds that limit,
      !  then reduce the perturbation accordingly without introducing a bias
      !  -------------------------------------------------------------------

      IF ((zlim_upp < clchk) .OR. (-zlim_low < clchk)) THEN
        !   no perturbation if the obs value is already at or beyond the limit
        IF (     (fb(iobs)% obs >= zlim_upp)                                   &
            .OR. (fb(iobs)% obs <= zlim_low)) THEN
          lpert = .FALSE.
        !   if adding the modulus of the perturbation to the obs value exceeds
        !   the upper limit, then limit the modulus of the perturbation such
        !   that adding the modulus of the limited perturbation to the obs value
        !   is equal to the upper limit
        !   --> by limiting the modulus of both positive and negative
        !       perturbations, the introduction of a bias is avoided
        ELSEIF (fb(iobs)% obs + zapert > zlim_upp) THEN
          fb(iobs)% obs  =  fb(iobs)% obs  +  SIGN( zlim_upp - fb(iobs)% obs       &
                                                  , REAL(dp_random_gauss(iobs),wp) )
          lpert = .FALSE.
        !   equivalent procedure for the lower limit
        ELSEIF (fb(iobs)% obs - zapert < zlim_low) THEN
          fb(iobs)% obs  =  fb(iobs)% obs  -  SIGN( zlim_low - fb(iobs)% obs       &
                                                  , REAL(dp_random_gauss(iobs),wp) )
          lpert = .FALSE.
        ENDIF
      ENDIF

      !  otherwise: add regular perturbation value
      !  -----------------------------------------

      IF (lpert)  fb(iobs)% obs  =  fb(iobs)% obs  +  zpert
    ENDIF
  ENDDO

!-------------------------------------------------------------------------------
! End Subroutine obs_fof_perturb
!-------------------------------------------------------------------------------
END SUBROUTINE obs_fof_perturb


!===============================================================================
!+ Module procedure in "src_obs_fdbk_in" for perturbing observation values
!-------------------------------------------------------------------------------

SUBROUTINE construct_seed  ( iseed )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure constructs a total seed number for a random number
!   generator based on the input seed and on date and time of the reference
!   time 'ydate_ref'. From the input seed, only the first 'ibts_ex' = 16 bits
!   (from the right) are used.
!
! Method:
!   Based on 'ydate_ref', a different number is constructed for every 5 minutes.
!   From that, only the first 'ibts_t' bits (from the right) are used.
!   For ibts_t = 14 , the same time seed number is produced about every
!   56.9 days.
!
! Written by        : DWD, Christoph Schraff, 29.04.2013
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER (KIND=iintegers) , INTENT (INOUT)  ::  &
    iseed        ! seed number for random number generator:
                 !   input : external seed number (only 'ibts_ex' bits used)
                 !   output: total seed number, which is a function of input
                 !           seed as well as date and time

! Local parameters:
! ----------------

!US  REAL    (KIND=wp)       ,   PARAMETER       ::
!US  these variables are input to the intrinsic procedure IBITS and MUST be
!    integer values!!!
  INTEGER (KIND=iintegers),   PARAMETER       ::  &
    ibts_t   =  14   ,& ! number of bits used from date and time
    ibts_ex  =  16      ! number of bits used from external seed number

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    imoyy          ,& ! year     / of the date
    imomm          ,& ! month   /  and time
    imodd          ,& ! day    <   of the reference
    imohh          ,& ! hour    \  (i.e. initial)
    imomin         ,& ! minute   \ model date and time
    iseed_t           ! 'ibts_t'-bit seed related to date and time
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine construct_seed
!-------------------------------------------------------------------------------

  READ (ydate_ref,'(2X,5I2)') imoyy, imomm, imodd, imohh, imomin
  !   unique value every 5 minutes within century, 0 <= iseed_time < 2^23
  iseed_t = (((imoyy *12 + imomm-1) *31 + imodd-1) *24 + imohh) *12  + imomin/5

  !   reduce value of iseed_t to less than 2^ibts_t : MOD( iseed_t , 2^ibts_t )
  !     (ibts_t = 14  -->  value repeats itself after about 56.9 days)
  !   from input value 'iseed', only first 'ibts_ex' bits (from right) are used
  iseed   =   ISHFT( IBITS( iseed   , 0 , ibts_ex ) , ibts_t )                 &
            +        IBITS( iseed_t , 0 , ibts_t  )

!-------------------------------------------------------------------------------
! End Subroutine construct_seed
!-------------------------------------------------------------------------------
END SUBROUTINE construct_seed

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

END MODULE src_obs_fdbk_in
