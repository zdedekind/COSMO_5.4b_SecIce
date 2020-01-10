!+ Data module for all data, that are used by the data assimilation
!-------------------------------------------------------------------------------

MODULE data_obs_cdfin

!-------------------------------------------------------------------------------
! Description:
!   This module declares and initializes all parametric data and data arrays
!   that are required to read and pre-process observation reports from
!   NetCDF Observation Input Files (except for those variables that reside
!   in module data_obs_lib_cosmo).
!   It contains the following sections:
!
!    1. NetCDF Observation Input File formats, and attributes of input
!       'fof' feedobs (feedback) files
!    2. Blacklist and Whitelist
!    3. Event counters: format, arrays, and description
!    4. Observation errors
!    5. Different kinds of limits
!    6. Variables used for the production of aircraft multi-level reports
!    7. For reporting rejection of data: Output buffer, size and formats
!    8. Temporary global model fields
!
!   Note: This module is part of the 'COSMO data assimilation library 1'
!         for reading data from NetCDF observation input files.
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
! V4_22        2012/01/31 Christoph Schraff
!  Initial version, extracted from module 'data_obs_process' and adapted.
! V4_28        2013/07/12 Christoph Schraff
!  Added variables used for (option.) reading obs from feedobs (feedback) files.
! V5_1         2014-11-28 Christoph Schraff, Oliver Fuhrer
!  Extensions for reading / processing Mode-S aircraft observations, from NetCDF
!  files with 2 different templates: (i) converted from BUFR as received from
!  Mode-S processing centre KNMI, (ii) converted into DWD-ACARS template. (CS)
!  Replaced ireals by wp (working precision) (OF)
! V5_1a        2015-03-30 Christoph Schraff
!  Specific horizontal and vertical check limits ('rhzlshp', 'rvtlshp')
!  introduced for redundancy check of ship and buoy reports.
! V5_3         2015-10-09 Christoph Schraff
!  Variables for solar radiation added to the entries of the NetCDF input file.
!  Specific temporal and vertical check limits ('rtmlsy', 'rvtlsy')
!  introduced for redundancy check of SYNOP reports.
! V5_4         2016-03-10 Christoph Schraff
!  Dimension of 'neventr' and 'neventd' reduced from 3 to 2.
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

!-------------------------------------------------------------------------------

USE netcdf
 
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
! Section 1 : NetCDF Observation Input File formats, and feedobs file attributes
!-------------------------------------------------------------------------------

!         1.1   internal attributes of the different NetCDF input files
!               -------------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    ntype_cdfin = 25  ,& ! number of types of NetCDF observation input files
    mxcdfin     = 80  ,& ! maximum number of NetCDF observation input files
    mxfofin     = 15  ,& ! max. number of NetCDF input feedobs (feedback) files
    icdfinlen   = 16  ,& ! maximum length of NetCDF observation input file name
    iannexlen   =  6  ,& ! maximum length of annex of NetCDF obs input file name
    mxtotin     = mxcdfin + mxfofin ! total maximum number of NetCDF input files

  INTEGER (KIND=iintegers) , PARAMETER  :: &
                        ! (the values of the following parameters have to
                        !  correspond to the entries in 'ycdfin')
    ncdf_temp      =  1 ,& ! indicator for processing of NetCDF TEMP       input
    ncdf_tempship  =  2 ,& ! indicator for processing of NetCDF TEMPSHIP   input
    ncdf_tempdrop  =  3 ,& ! indicator for proc. NetCDF TEMP Dropsonde     input
    ncdf_pilot     =  4 ,& ! indicator for proc. NetCDF PILOT (z-levels)   input
    ncdf_pilot_p   =  5 ,& ! indicator for proc. NetCDF PILOT (p-levels)   input
    ncdf_amdar_ml  =  6 ,& ! indicator for proc. NetCDF AMDAR multi-level  input
    ncdf_amdar_vp  =  7 ,& ! indicator for proc. NetCDF AMDAR vert.profile input
    ncdf_amdar     =  8 ,& ! indicator for proc. NetCDF AMDAR single-level input
    ncdf_acars     =  9 ,& ! indicator for proc. NetCDF ACARS single-level input
    ncdf_acars_uk  = 10 ,& ! indicator for proc. NetCDF ACARS UK + Canada  input
    ncdf_acars_us  = 11 ,& ! indicator for proc. NetCDF ACARS US w. humid. input
    ncdf_modes     = 12 ,& ! indicator for proc. NetCDF MODE-S KNMI format input
    ncdf_modes_acr = 13 ,& ! indicator for proc. NetCDF MODE-S ACARS fmt.  input
    ncdf_wprof     = 14 ,& ! indicator for proc. NetCDF wind profiler      input
    ncdf_rass      = 15 ,& ! indicator for proc. NetCDF RASS profiler      input
    ncdf_radar_vad = 16 ,& ! indicator for proc. NetCDF radar wind prof.   input
    ncdf_synop     = 17 ,& ! indicator for proc. NetCDF SYNOP              input
    ncdf_synop_mob = 18 ,& ! indicator for proc. NetCDF SYNOP mobile       input
    ncdf_ship      = 19 ,& ! indicator for proc. NetCDF SHIP               input
    ncdf_buoy      = 20 ,& ! indicator for proc. NetCDF BUOY               input
    ncdf_metar     = 21 ,& ! indicator for proc. NetCDF METAR sfc aviation input
    ncdf_gps_zenith= 22 ,& ! indicator for proc. NetCDF GPS (ZPD / IWV)    input
    ncdf_ascat     = 23 ,& ! indicator for proc. NetCDF ASCAT scatterom.   input
    ncdf_qscat     = 24 ,& ! indicator for proc. NetCDF QuickScat scatter. input
    ncdf_satob     = 25 ,& ! indicator for proc. NetCDF SATOB wind         input
    ncdf_fof       = -1    ! indicator for proc. 'fof' feedobs (fdbk) file input

                        ! file names of NetCDF observation input files
                        !   (the order of the entries has to correspond to
                        !    the values of the above parameters 'ncdf_*')
  CHARACTER  (LEN=icdfinlen) , PARAMETER  :: &
    ycdfin (ntype_cdfin) = (/'cdfin_temp      '                               ,&
                             'cdfin_tempship  '                               ,&
                             'cdfin_tempdrop  '                               ,&
                             'cdfin_pilot     '                               ,&
                             'cdfin_pilot_p   '                               ,&
                             'cdfin_amdar_ml  '                               ,&
                             'cdfin_amdar_vp  '                               ,&
                             'cdfin_amdar     '                               ,&
                             'cdfin_acars     '                               ,&
                             'cdfin_acars_uk  '                               ,&
                             'cdfin_acars_us  '                               ,&
                             'cdfin_modes     '                               ,&
                             'cdfin_modes_acr '                               ,&
                             'cdfin_wprof     '                               ,&
                             'cdfin_rass      '                               ,&
                             'cdfin_radar_vad '                               ,&
                             'cdfin_synop     '                               ,&
                             'cdfin_synop_mob '                               ,&
                             'cdfin_ship      '                               ,&
                             'cdfin_buoy      '                               ,&
                             'cdfin_metar     '                               ,&
                             'cdfin_gps_zenith'                               ,&
                             'cdfin_ascat     '                               ,&
                             'cdfin_qscat     '                               ,&
                             'cdfin_satob     '                               /)

  CHARACTER  (LEN=icdfinlen) , PARAMETER  :: &
    yfofin    = 'fof'   ! file name of feedobs (feedback) input file(s)
                        !   note: this has no date/time in the file name
                        !         in contrast to the feedobs output file

  INTEGER (KIND=iintegers)  :: &
    n_cdfin              ,& ! number of existing NetCDF observation input files
    n_fofin              ,& ! number of existing input feedobs (feedback) files
    icdfin (mxcdfin) = 0 ,& ! obs file type of NetCDF observation input files
    ncinid (mxcdfin) = 0 ,& ! unit numbers of NetCDF observation input files
    ncid                    ! file unit number of current NetCDF file

  CHARACTER (LEN=iannexlen) :: &
    yncannex  (mxcdfin) ,& ! annex of NetCDF observation input file names
    yfofannex (mxfofin)    ! annex of NetCDF input feedobs (feedback) file names

  INTEGER (KIND=iintegers)  :: &
    dimids (nf90_max_dims) ! dimension IDs in NetCDF files

!         1.2   Variables used to read report header entries from NetCDF files
!               --------------------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    ilstidn  = 10  ,& ! character length of station identity from NetCDF files
    ichar10  = 10     ! character length of station identity from feedback files

!         1.2.1 'common' NetCDF header entries and derived variables
!               ----------------------------------------------------
!               (i.e. entries and variables related to internal (ODR)
!                header elements which are common to all obs types)

  INTEGER (KIND=iintegers), ALLOCATABLE :: &
    nc_dim_len (:) ,& ! length of a dimension          (only temporary variable)
! BUFR section 1 and section 2 entries  (common to all obs types, stored in ODR
!                                             except if only temporary variable)
    kcat       (:) ,& ! data category     (WMO Common Code Table C13)
    kcatsub    (:) ,& ! data sub-category (WMO Common Code Table C13)
    kcentre    (:) ,& ! data centre       (WMO Common Code Table C11)
    kcensub    (:) ,& ! data sub-centre   (WMO C12)    (only temporary variable)
    kupdate    (:) ,& ! update sequence number (indicates station correction)
    kz_dwd     (:) ,& ! DWD-internal classification number (Kennzahl KZ)
    nsynhr     (:) ,& ! nominal (synoptic) hour [yymmddhh] (from s1_date,_time)
! common header entries in NetCDF file
    nc_mjjj    (:) ,& ! MJJJ   : year                  (only temporary variable)
    nc_mmm     (:) ,& ! MMM    : month                 (only temporary variable)
    nc_myy     (:) ,& ! MYY    : day                   (only temporary variable)
    nc_mgg     (:) ,& ! MGG    : hour                  (only temporary variable)
    nc_ngg     (:) ,& ! NGG    : minute                (only temporary variable)
    nc_msec    (:) ,& ! MSEC   : second                (only temporary variable)
    nc_mii     (:) ,& ! MII    : WMO block number      (only temporary variable)
    nc_niii    (:) ,& ! NIII   : WMO station number    (only temporary variable)
    nc_alt     (:) ,& ! MHP, MH: station altitude      (only temporary variable)
    nc_tisi    (:) ,& ! MTISI  : time significance (BUFR Table 008021)    (dito)
    nc_altq    (:) ,& ! MSEQM  : altitude quality (mobil, BUFR 033024)    (dito)
    nc_locq    (:) ,& ! MQOBL  : location quality (buoy , BUFR 033023)    (dito)
! derived common header variables stored to ODR
    iobtyp     (:) ,& ! observation type
    icdtyp     (:) ,& ! observation code type
    istidn     (:) ,& ! station number
    nr_date    (:) ,& ! absolute date [YYYYMMDD]
    nr_time    (:) ,& ! absolute time [HHMM]
    isurfob    (:) ,& ! indic. if orographic surface report conditions are met 
    ksurfob    (:) ,& ! surface report indicator and report index
    mrepflg    (:) ,& ! report flags
    kflag      (:) ,& ! processing flag    (bit pattern as for 'nhflag' in ODR)
    iobstot    (:) ,& ! longitudinal index of grid pt. to which obs is assigned
    jobstot    (:) ,& ! latitudinal  index of grid pt. to which obs is assigned
    iobsloc    (:) ,& ! longitudinal index of grid pt. in local sub-domain
    jobsloc    (:) ,& ! latitudinal  index of grid pt. in local sub-domain
! auxilliary variable, only temporarily available in reader routine
    irproc     (:)    ! indices of reports to be processed now
! derived common header variables, only temporarily available in reader routine
!   kproc      (:) ,& ! flags which indicate no further processing
!   nt_date    (:) ,& ! absolute date [YYYYMMDD]
!   nt_time    (:) ,& ! absolute time [HHMM]
!   min_ob     (:)    ! time relative to model initial time [min]

  REAL    (KIND=wp)       , ALLOCATABLE :: &
! common header entries in NetCDF file     (stored to ODR or temporary variable)
    rc_lat     (:) ,& ! MLAH, MLALA : latitude
    rc_lon     (:) ,& ! MLOH, MLOLO : longitude
    rc_altp    (:) ,& ! MHOBNN      : barometer altitude    (temporary variable)
    rc_alt     (:) ,& ! MHOSNN      : station   altitude    (temporary variable)
! derived common header variables stored to ODR
    zr_hour    (:) ,& ! observation time relative to model initial time [hour]
    zsynhr     (:) ,& ! nominal report time rel.  to model initial time [hour]
    zaltob     (:) ,& ! obs station altitude
    zaltmo     (:) ,& ! model orography (at the grid point assigned to the obs)
    rio_tot    (:) ,& ! longitude in model grid point units
    rjo_tot    (:) ,& ! latitude  in model grid point units
    ztdecdb    (:)    ! data base decoding time rel to model initial time [hour]

  CHARACTER (LEN= nf90_max_name) , ALLOCATABLE :: &
    yc_dim_name(:)    ! name of dimension in NetCDF file

  CHARACTER (LEN= 1)             , ALLOCATABLE :: &
    ync_stid (:,:)    ! YDDDD, YSSOSN, YAIRN, YCCC8, YSOSN: station identity
                      !                          (temporary auxilliary variable)

  CHARACTER (LEN= ilstidn)       , ALLOCATABLE :: &
! derived common header variable stored to ODR
    ystidn     (:)    ! station identity

!         1.2.2 other NetCDF header entries
!               ---------------------------

  INTEGER (KIND=iintegers), ALLOCATABLE :: &
    nc_rstyp   (:) ,& ! NRARA , WMO Common Table C2 : radiosonde type/system
    nc_rad     (:) ,& ! NSR   , BUFR Table B 002013 : solar + IR radiation corr.
    nc_track   (:) ,& ! NSASA , WMO Common Table C7 : tracking technique, status
    nc_na4     (:) ,& ! NA4   , BUFR Table B 002003 : type of measur. equipment
    nc_nix     (:) ,& ! NIX   , BUFR Table B 002001 : station type (man,auto,..)
    nc_tbuoy   (:) ,& ! MTODB , BUFR Table B 002149 : type of data buoy
    nc_phase   (:) ,& ! MPHAI , BUFR Table B 008004 : phase of aircraft flight
    nc_phasd   (:)    ! NDEPF , BUFR Table B 008009 : detailed phase of flight

!         1.3   NetCDF body entries
!               -------------------

!         1.3.1  freqent entries
!                ---------------

  INTEGER (KIND=iintegers), ALLOCATABLE :: &
    nc_nlev    (:) ,& ! MEDRE : delayed descriptor replication factor
                      !         (e.g. number of vertical levels)
    nc_dt    (:,:) ,& ! NLTPD : time [sec] since launch time
    nc_lvtyp (:,:) ,& ! MEVSS , BUFR Tab 008042 : extended vertical sounding
                      !                           significance (level identity)
                      !         (bits in reverse order !!)
    nc_z     (:,:) ,& ! NHHHN : geopotential height [gpm]    (upper-air)
                      ! NFLEV , MIAA : flight level / indicated altitude
    nc_dd    (:,:)    ! NDNDN : wind direction      [degree]
                      !
  REAL    (KIND=wp)       , ALLOCATABLE :: &
    rc_p     (:,:) ,& ! MPN   , MPPP : pressure
    rc_dlat  (:,:) ,& ! MLADH : latitude  displacement since launch site
    rc_dlon  (:,:) ,& ! ML0DH : longitude displacement since launch site
    rc_t     (:,:) ,& ! MTDBT : temperature / dry-bulb temperature
    rc_td    (:,:) ,& ! MTDNH : dew-point temperature
    rc_ff    (:,:) ,& ! NFNFN : wind speed  (or NFF for scatterometer)
    rc_dd    (:,:) ,& ! NDD   : wind direction (for scatterometer only)
    rc_lat0  (:,:) ,& ! MLAH0 : latitude
    rc_lon0  (:,:)    ! ML0H0 : longitude

!         1.3.2  additional aircraft elements
!                ----------------------------

  INTEGER (KIND=iintegers), ALLOCATABLE :: &
    nc_rolla (:,:) ,& ! MQARA , BUFR Tab 002064 : aircraft roll angle quality
    nc_turb    (:) ,& ! MB    , BUFR Tab 011031 : degree of turbulence
    nc_rh      (:) ,& ! MUUU  : relative humidity                           [%]
    nc_qxq     (:)    ! MMRQ  , BUFR Tab 033026 : mixing ratio quality
                      !
  REAL    (KIND=wp)       , ALLOCATABLE :: &
    rc_z       (:) ,& ! MHHH  : height or altitude (vertical location)
    rc_qx      (:) ,& ! MMIXR : mixing ratio                            [kg/kg]
    rc_tq      (:) ,& ! MPOTO : precision of temperature observation        [K]
    rc_vgust   (:)    ! NMDEWX: maximum derived equivalent vertical gust  [m/s]

!         1.3.3  additional profiler elements
!                ----------------------------

  INTEGER (KIND=iintegers), ALLOCATABLE :: &
    nc_sinor (:,:) ,& ! MSINOR, BUFR Tab 002064 : signal to noise ratio
    nc_qci   (:,:) ,& ! MQINZ , BUFR Tab 033002 : quality information
    nc_qci2  (:,:) ,& ! NWPQ  , BUFR Tab 025034 : NOAA QC results   (temporary)
    nc_wce     (:) ,& ! MWCE  , BUFR Tab 025021 : wind computation enhancement
    nc_dto     (:) ,& ! MSETP , time period of measurement [sec]
    nc_dto2    (:)    ! NGGTP , time period of measurement [min]    (temporary)
                      !
  REAL    (KIND=wp)       , ALLOCATABLE :: &
    rc_w     (:,:) ,& ! MWMPS : vertical velocity (w-component of wind)   [m/s]
    rc_tv    (:,:) ,& ! MTVIR : virtual temperature                         [K]
    rc_stdff (:,:)    ! NSTDFF: standard deviation wind speed             [m/s]

!         1.3.4  additional synoptic elements
!                ----------------------------

  INTEGER (KIND=iintegers), ALLOCATABLE :: &
    nc_uuu     (:) ,& ! MUUU  : relative humidity
    nc_fi      (:) ,& ! NHHHN : geopotential height [gpm] of standard level
    nc_vtisi   (:) ,& ! MTISI*: BUFR Tab 008021 : time significance (wind obs)
    nc_vdt     (:) ,& ! NGGTP : time period of wind measurement           [min]
    nc_gudt    (:) ,& ! NGGTP0: time period of (1st) max. wind gust speed [min]
    nc_gudt2   (:) ,& ! NGGTP0: time period of  2nd  max. wind gust speed [min]
    nc_tmxt1   (:) ,& ! MGGTP*: time displacement: start of period of T-max [h]
    nc_tmxt2   (:) ,& ! MGGTP*: time displacement: end   of period of T-max [h]
    nc_tmnt1   (:) ,& ! MGGTP*: time displacement: start of period of T-min [h]
    nc_tmnt2   (:) ,& ! MGGTP*: time displacement: end   of period of T-min [h]
    nc_rrdt    (:) ,& ! MGGTP*: time period of (first) precipitation obs    [h]
    nc_rrdt2   (:) ,& ! MGGTP*: time period of  second precipitation obs    [h]
    nc_clct    (:) ,& ! MN                      : total cloud amount [%]
    nc_clsig   (:) ,& ! MVTSU , BUFR Tab 008002 : vertical significance
    nc_clxsg (:,:) ,& ! MVTSU*, BUFR Tab 008002 :    dito    (additional levels)
    nc_clclm   (:) ,& ! MNH   , BUFR Tab 020011 :(low or mid-level) cloud amount
    nc_clxcl (:,:) ,& ! MNH*  , BUFR Tab 020011 :    dito    (additional levels)
    nc_ccl     (:) ,& ! MCC   , BUFR Tab 020012 : cloud type (low clouds)
    nc_clxct (:,:) ,& ! MCC*  , BUFR Tab 020012 :    dito    (additional levels)
    nc_ccm     (:) ,& ! MCC0  , BUFR Tab 020012 : cloud type (middle clouds)
    nc_cch     (:) ,& ! MCC1  , BUFR Tab 020012 : cloud type (high clouds)
    nc_clev    (:) ,& ! MDREP*, replication factor (for additional cloud levels)
    nc_wdt     (:) ,& ! MGGTP : time period [h] for past weather
    nc_ww      (:) ,& ! WW    , BUFR Tab 020003 : present weather
    nc_w1      (:) ,& ! W1    , BUFR Tab 020004 : past    weather (1)
    nc_w2      (:) ,& ! W2    , BUFR Tab 020005 : past    weather (2)
    nc_e       (:)    ! ME    , BUFR Tab 020062 : state of ground (w/wo snow)

  REAL    (KIND=wp)       , ALLOCATABLE :: &
    rc_pfi     (:) ,& ! MPN   : pressure of standard level  (surface report)
    rc_pmsl    (:) ,& ! MPPPP : pressure reduced to mean sea level
    rc_dpdt    (:) ,& ! NPPP  : 3-hour pressure change
    rc_gust    (:) ,& ! NFXGU : maximum wind speed of gusts
    rc_gust2   (:) ,& ! NFXGU : maximum wind speed of gusts, second period
    rc_tmax    (:) ,& ! MTXTXH: maximum temperature \ height, period specified,
    rc_tmin    (:) ,& ! MTNTNH: minimum temperature / processed if: 2-m , 12-h
    rc_zt      (:) ,& ! MHOSEN*:height of T-sensor above local ground [m]
    rc_vis     (:) ,& ! MVV(VV):horizontal visibility
    rc_rr24    (:) ,& ! MRR24 : total precipitation over past 24 hours
    rc_rr      (:) ,& ! MRRR  : total precipitation (first period)
    rc_rr2     (:) ,& ! MRRR  : total precipitation, second period
    rc_cbase   (:) ,& ! NH    : height of base of cloud
    rc_clxbs (:,:) ,& ! NH*   : height of base of cloud (additional levels)
    rc_radgl   (:) ,& ! MGLSR : global  solar radiation                  [J/m2]
    rc_raddf   (:) ,& ! MDSRH : diffuse solar radiation                  [J/m2]
    rc_radlw   (:) ,& ! MLWR  : long-wave     radiation                  [J/m2]
    rc_hsnow   (:)    ! NSSS  : total snow depth [m]

!         1.3.5  additional gps elements      
!                -----------------------

  INTEGER (KIND=iintegers), ALLOCATABLE :: &
    nc_flg     (:)    ! NQFGD : quality flag

  REAL    (KIND=wp)       , ALLOCATABLE :: &
    rc_el      (:) ,& ! MDE   : elevation                                   [°]
    rc_iwv   (:,:) ,& ! NWLN  : integrated water vapour                 [kg/m²]
    rc_zwd   (:,:) ,& ! NCZWV : zenith wet delay                            [m]
    rc_zpder   (:) ,& ! NEERR : estimated error in zpd                      [m]
    rc_zpd     (:)    ! NADES : zenith path delay                           [m]

  INTEGER (KIND=iintegers) , PARAMETER  :: &

!         1.4.1  Bit positions for level significance (MEVSS, BUFR Tab 008042)
!                ------------------------------------ (bits in reverse order!)

    ilv_sfc   = 17 ,& ! surface                 level bit (i.e. 2^17)
    ilv_std   = 16 ,& ! standard                level bit
    ilv_tropo = 15 ,& ! tropopause              level bit
    ilv_max   = 14 ,& ! maximum wind            level bit
    ilv_sigt  = 13 ,& ! significant temperature level bit
    ilv_sigq  = 12 ,& ! significant humidity    level bit
    ilv_sigv  = 11 ,& ! significant wind        level bit
    ilv_miss  = 18    ! missing value indicator       bit

!         1.5   missing value specification
!               ---------------------------

  INTEGER (KIND=iintegers)  :: &
    imiss             ! missing value for integers in current NetCDF input file

  REAL    (KIND=wp)         :: &
    rmiss          ,& ! missing value for reals    in current NetCDF input file
    rmisschk          ! value smaller than 'rmiss', to check for missing value


!-------------------------------------------------------------------------------
! Section 2 : Blacklist and Whitelist
!-------------------------------------------------------------------------------

!         2.1   global blacklist and whitelists as read from file
!               -------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    ilstid_blk  =  8    ,& ! assume 8-character station-IDs in Black-/Whitelist
    mxot_blk    = 12    ,& ! (max.) number of observation types in blacklist
                           !   ( should be ' = mxobtp ' in module
                           !     'data_obs_lib_cosmo.f90' )
    mx_wits     = 20       ! max. number of whitelists (1 per obs code type)

  CHARACTER (LEN=icdfinlen), PARAMETER  :: &
    yblk_in = 'blklsttmp'  ! file name of blacklist file (which must reside in
                           ! the same directory as the NetCDF obs input files !)

  INTEGER (KIND=iintegers)  :: &
    nblack              ,& ! (total) number of blacklisted reports
    nwhite              ,& !  total  number of whitelisted stations
    n_wits              ,& ! actual  number of whitelisted stations
    kwit_obt  (mx_wits) ,& ! observation type of whitelist
    kwit_cdt  (mx_wits) ,& ! obs  code   type of whitelist
    kwit_frst (mx_wits) ,& ! index of first station of current whitelist
    kwit_last (mx_wits) ,& ! index of last  station of current whitelist
    iblk_frst (mxot_blk),& ! index of first report of current obs type
    iblk_last (mxot_blk)   ! index of last  report of current obs type

  CHARACTER (LEN=ilstid_blk) , ALLOCATABLE :: &
    yblk_id    (:) ,& ! station identity of blacklisted stations
    ywit_id    (:)    ! station identity of whitelisted stations

  INTEGER (KIND=iintegers), ALLOCATABLE :: &
    iblk_len   (:) ,& ! length of 'yblk_id'
    iwit_len   (:) ,& ! length of 'ywit_id'
    iblk_pts   (:) ,& ! number of points (wild cards) in 'yblk_id'
    iwit_pts   (:) ,& ! number of points (wild cards) in 'ywit_id'
    iblk_obtyp (:)    ! observation type of blacklisted report

  REAL    (KIND=wp)       , ALLOCATABLE :: &
    rblk_pzlow (:) ,& ! lower  \  boundary (pressure) of blacklisted interval
    rblk_pzup  (:) ,& ! upper  /                      for geopotential
    rblk_pvlow (:) ,& ! lower  \  boundary (pressure) of blacklisted interval
    rblk_pvup  (:) ,& ! upper  /                      for wind
    rblk_ptlow (:) ,& ! lower  \  boundary (pressure) of blacklisted interval
    rblk_ptup  (:) ,& ! upper  /                      for temperature
    rblk_pqlow (:) ,& ! lower  \  boundary (pressure) of blacklisted interval
    rblk_pqup  (:)    ! upper  /                      for humidity

!         2.2   array of separate blacklists which are each for 1 local report
!               --------------------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    maxintv    = 16     ! max. number of vertical blacklist intervals per 1
                        !   station ID (for each variable, a separate interval
                        !   must be defined)

  TYPE blacklist_loc
    INTEGER       :: ndim            ! number of blacklisted intervals in report
    INTEGER       :: kvar (maxintv)  ! variable type of blacklisted interval
                                     !   1: wind, 2: geopot, 3: temperat., 4: humidity
    REAL(KIND=wp) :: plow (maxintv)  ! lower boundary (pressure) of blacklist interval
    REAL(KIND=wp) :: pup  (maxintv)  ! upper boundary (pressure) of blacklist interval
  END TYPE blacklist_loc

  TYPE (blacklist_loc) , DIMENSION(:) , ALLOCATABLE  ::                        &
    blk_loc(:)      ! blacklists for local reports
 
!-------------------------------------------------------------------------------
! Section 3 : Event counters: format, arrays, and description
!-------------------------------------------------------------------------------
 
!         3.1    Format of event counters
!                ------------------------
  
!         3.1.1  Report event counter array format
!                ---------------------------------
 
  INTEGER (KIND=iintegers) , PARAMETER  :: &
    mxreve = 17  ,& ! length of report event counter array
                    !
    nedbfl =  1  ,& ! data base flag on loc/tim/alt high
    netime =  2  ,& ! time out of range (too old)
    nenoal =  3  ,& ! no station altitude
    neloca =  4  ,& ! station location out of domain or out of user-specif. area
    nezdif =  5  ,& ! distance 'model orography - station altitude' too large
    neblak =  6  ,& ! blacklisted ship
    neobct =  7  ,& ! obs. or code type excluded on area with sta. location
    nesodr =  8  ,& ! report number exceeding size of ODR (==> adjust namelist)
    nenoda =  9  ,& ! no accepted data in report
    nenops = 10  ,& ! pressure too small (< 20hPa), or missing with aircraft rep
    neredn = 11  ,& ! redundancy between 2 multi-level or 2 single-level reports
    netrac = 12  ,& ! (flight) track suspicious
    nethin = 13  ,& ! thinning of aircraft (flight) track 
    neredx = 14  ,& ! redundancy between 1 multi- and 1 single-level report
    neslml = 15  ,& ! one multi-level report made from other reports (either
                    !   from single-level aircraft, or from TEMP parts A,B,C,D)
    neslps = 16  ,& ! report (either single-level aircraft, or TEMP part A,B,C,
                    !   or D) put in multi-level report and set passive
    nenoml = 17     ! multi-levl report not made due to ODR array size

!         3.1.2  Data event counter array format
!                -------------------------------
 
  INTEGER (KIND=iintegers) , PARAMETER  :: &
    mxdeve = 38  ,& ! length of data event counter array
                    !
    nelodr =  1  ,& ! level rejected: number of levels exceeding ODR size
    nelmis =  2  ,& ! level rejected: pressure (PILOT: pressure +height) missing
    nelflg =  3  ,& ! level rejected: pressure (PILOT: height) flagged
    nelsfc =  4  ,& ! level rejected: too many surface levels
    nelnop =  5  ,& ! level rejected: PILOT height level out of model lev. range
    nelext =  6  ,& ! level rejected: pressure < 9hPa or level below sta. height
    nelsig =  7  ,& ! level rejected: significant level above a specified limit
    nelrdn =  8  ,& ! level rejected: redundant level in report (not active yet)
    nepmis =  9  ,& ! pressure (TEMP: height): missing
    nepflg = 10  ,& ! pressure (TEMP: height): flagged
    neprac = 11  ,& ! pressure: bad reporting practice
    nepalt = 12  ,& ! pressure: stat. height, or distance to orography too large
    nepsdt = 13  ,& ! pressure tendency: flagged, or absolute value too large
    netmis = 14  ,& ! temperature missing
    netflg = 15  ,& ! temperature flagged
    netext = 16  ,& ! temperature too low or too high
    netalt = 17  ,& ! height (diff.) too large for 2m-temp.
    netlps = 18  ,& ! lapse rate of multi-level temperature too large
    neqmis = 19  ,& ! humidity missing
    neqflg = 20  ,& ! humidity flagged
    neqlow = 21  ,& ! humidity too low
    neq300 = 22  ,& ! humidity above 300 hpa
    neqbig = 23  ,& ! humidity over allowed value (120%)
    neqsap = 24  ,& ! humidity forced to be saturated (t>0)
    neqsam = 25  ,& ! humidity forced to be saturated (t<0)
    neqclp = 26  ,& ! humidity forced to be <=100% (t>0)
    neqclm = 27  ,& ! humidity forced to be <=100% (t<0)
    neqalt = 28  ,& ! height (diff.) too large for 2m-humid
    nedmis = 29  ,& ! wind direction missing
    nefmis = 30  ,& ! wind speed missing
    nedflg = 31  ,& ! wind direction flagged
    nefflg = 32  ,& ! wind speed flagged
    nefneg = 33  ,& ! wind speed too small  ( < 0 ; DRIBU <= 0 ; VAD < 3m/s )
    nevalt = 34  ,& ! height (diff.) too large for 10m-wind
    nefshr = 35  ,& ! wind speed shear too large
    nedshr = 36  ,& ! directional wind shear too large
    nerlim = 37  ,& ! precipitation amount exceeds threshold limit
    negmis = 38     ! ZPD missing or too small 

!         3.2    Event counter arrays
!                --------------------

  INTEGER (KIND = iintegers), ALLOCATABLE :: &
    neventr  (:,:) ,& ! counter of report events
    neventd  (:,:)    ! counter of data events
 
!         3.3    Character descriptions of events and flags
!                ------------------------------------------

!         3.3.1  Report events
!                -------------

  CHARACTER (LEN=72)       , PARAMETER  :: &
                    ! Description of report events
    crepev (mxreve) =                                                          &
                    !
  (/' 1 = DATA BASE FLAG ON LOCATION / TIME / ALTITUDE HIGH                  ',&
    ' 2 = OBSERVATION TIME OUT OF RANGE (TOO OLD) (OR TIME MISSING)          ',&
    ' 3 = STATION ALTITUDE MISSING                                           ',&
    ' 4 = STATION LOCATION OUT OF DOMAIN OR OUT OF USER-SPECIFIED AREA       ',&
    ' 5 = DISTANCE "MODEL OROGRAPHY - STATION ALTITUDE" TOO LARGE            ',&
    ' 6 = BLACKLISTED SHIP                                                   ',&
    ' 7 = OBSERVATION OR CODE TYPE EXCLUDED IN AREA AROUND STATION LOCATION  ',&
    ' 8 = REPORT NUMBER EXCEEDING SIZE OF ODR (==> ADJUST NAMELIST !!!)      ',&
    ' 9 = NO ACCEPTED DATA IN REPORT                                         ',&
    '10 = PRESSURE TOO SMALL ( < 9 HPA), OR MISSING WITH AIRCRAFT REPORT     ',&
    '11 = REDUNDANCY BETWEEN 2 MULTI-LEVEL, OR 2 SINGLE-LEVEL REPORTS        ',&
    '12 = FLIGHT TRACK SUSPICIOUS, OR EXAGGERATED COLOCATION                 ',&
    '13 = THINNING OF DENSE AIRCRAFT FLIGHT TRACK                            ',&
    '14 = SURFACE LEVEL FROM MULTI-LEVEL REP. REDUNDANT AGAINST OTHER REPORT ',&
    '15 = ONE MULTI-LEVEL REPORT MADE FROM SINGLE-LEVEL REPORTS              ',&
    '16 = SINGLE-LEVEL REPORT OR TEMP PART PUT INTO MULTI-LEVEL REPORT       ',&
    '17 = MULTI-LEVEL REPORT NOT BUILT DUE TO ODR SIZE LIMIT: ADJUST NAMELIST'/)
 

!         3.3.2  Data events
!                -----------

  CHARACTER (LEN=72)       , PARAMETER  :: &
                    ! Description of data events
    cdatev (mxdeve) =                                                          &
                    !
  (/' 1 = LEVEL REJECTED: NUMBER OF LEVELS EXCEEDING ODR SIZE                ',&
    ' 2 = LEVEL REJECTED: PRESSURE (PILOT: PRESSURE AND HEIGHT) MISSING      ',&
    ' 3 = LEVEL REJECTED: PRESSURE (PILOT: HEIGHT) FLAGGED                   ',&
    ' 4 = LEVEL REJECTED: TOO MANY SURFACE LEVELS                            ',&
    ' 5 = LEVEL REJECTED: PILOT HEIGHT LEVEL OUTSIDE RANGE OF MODEL LEVELS   ',&
    ' 6 = LEVEL REJECTED: PRESSURE < 9 HPA, OR LEVEL BELOW STATION HEIGHT    ',&
    ' 7 = LEVEL REJECTED: SIGNIFICANT LEVEL ABOVE A SPECIFIED LIMIT          ',&
    ' 8 = LEVEL REJECTED: REDUNDANT LEVEL IN REPORT  (NOT ACTIVE YET)        ',&
    ' 9 = PRESSURE (TEMP: HEIGHT): MISSING                                   ',&
    '10 = PRESSURE (TEMP: HEIGHT): FLAGGED                                   ',&
    '11 = PRESSURE: BAD REPORTING PRACTICE                                   ',&
    '12 = PRESSURE: HEIGHT DISTANCE TO OROGRAPHY OR STATION HEIGHT TOO LARGE ',&
    '13 = PRESSURE TENDENCY: FLAGGED, OR ABSOLUTE VALUE > 40 HPA/3H          ',&
    '14 = TEMPERATURE: MISSING (TEMP: AT SIGNIFICANT TEMPERATURE LEVELS ONLY)',&
    '15 = TEMPERATURE: FLAGGED                                               ',&
    '16 = TEMPERATURE: < -90 C,  OR  > +60 C  (P < 700HPA: > +20 C , ETC)    ',&
    '17 = TEMPERATURE AT 2M: HEIGHT OR HEIGHT DISTANCE TO OROGRAPHY TOO LARGE',&
    '18 = TEMPERATURE (TEMP ONLY): LAPSE RATE TOO LARGE                      ',&
    '19 = HUMIDITY: MISSING (TEMP: AT SIGNIFICANT LEVELS BELOW 300 HPA LEVEL)',&
    '20 = HUMIDITY: FLAGGED                                                  ',&
    '21 = HUMIDITY: DEWPOINT < -150 C (SURFACE-LEV OBS: < -90 C), OR > +40 C ',&
    '22 = HUMIDITY: ABOVE 300 HPA LEVEL                                      ',&
    '23 = HUMIDITY: EXCEEDING ALLOWED VALUE (120%)                           ',&
    '24 = HUMIDITY: FORCED TO BE SATURATED (T>O)                             ',&
    '25 = HUMIDITY: FORCED TO BE SATURATED (T<O)                             ',&
    '26 = HUMIDITY: FORCED TO BE <= 100% (T>0)                               ',&
    '27 = HUMIDITY: FORCED TO BE <= 100% (T<0)                               ',&
    '28 = HUMIDITY AT 2M: HEIGHT OR HEIGHT DISTANCE TO OROGRAPHY TOO LARGE   ',&
    '29 = WIND DIRECTION: MISSING                                            ',&
    '30 = WIND SPEED: MISSING                                                ',&
    '31 = WIND DIRECTION: FLAGGED , OR ABSOLUTE VALUE > 360 DEGREES          ',&
    '32 = WIND SPEED: FLAGGED                                                ',&
    '33 = WIND SPEED: < 0 (DRIBU: <= 0) , OR > 150 M/S (P > 700HPA: > 90 M/S)',&
    '34 = WIND AT 10M: HEIGHT OR HEIGHT DISTANCE TO OROGRAPHY TOO LARGE      ',&
    '35 = WIND SPEED: SHEAR TOO LARGE                                        ',&
    '36 = WIND DIRECTION: SHEAR TOO LARGE                                    ',&
    '37 = PRECIPITATION: AMOUNT EXCEEDING THRESHOLD LIMIT                    ',& 
    '38 = ZENITH PATH DELAY MISSING OR TOO SMALL                             '/)

!         3.4    Variables' expectation table
!                ----------------------------
!                (event counter for missing value is increased
!                 only if value is expected)

  INTEGER (KIND=iintegers) , PARAMETER  :: &                      ! expect:
    nzex  (mxot_blk) = (/ 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1/) ,& ! geopotential
    nuex  (mxot_blk) = (/ 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0/) ,& ! horiz. wind
    ntex  (mxot_blk) = (/ 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1/) ,& ! temperature
    ntdex (mxot_blk) = (/ 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1/)    ! humidity


!-------------------------------------------------------------------------------
! Section 4 : Observation errors
!-------------------------------------------------------------------------------

!         4.1  Observation error levels
!              ------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    nerlev  =  15    ! number of standard error levels

  REAL (KIND=wp)            ::    &
    rolnlv (nerlev)    ! ln(rlevel(15))

  REAL (KIND=wp)           , PARAMETER  :: &
    rlevel (nerlev)  & ! levels of error tables
              = (/ 100000._wp, 85000._wp, 70000._wp, 50000._wp,&
                    40000._wp, 30000._wp, 25000._wp, 20000._wp,&
                    15000._wp, 10000._wp,  7000._wp,  5000._wp,&
                     3000._wp,  2000._wp,  1000._wp /)


!         4.2  Observation error constants
!              ---------------------------

  REAL (KIND=wp)           , PARAMETER  :: &
                     !
                     ! (root of) radiosonde (TEMP, PILOT) wind error variance
                     !           -----------------------------
    oevsond (nerlev) = (/ 2.0_wp,   2.4_wp,   2.5_wp,   3.4_wp,&
                          3.5_wp,   3.7_wp,   3.5_wp,   3.5_wp,&
                          3.4_wp,   3.3_wp,   3.2_wp,   3.2_wp,&
                          3.3_wp,   3.6_wp,   4.5_wp /)           ,&
                     ! (values in O.I.: 2.0, 2.4, 2.5, 3.4, 3.6, 3.8, 3.2, 3.2
                     !                  2.4, 2.2, 2.0, 2.0, 2.5, 3.0, 4.0 )
                     ! (values up to  : 2.0, 2.4, 2.5, 3.4, 3.6, 3.8, 3.2, 3.2
                     !       LM 2.12    2.4, 2.2, 2.0, 2.0, 2.0, 2.5, 3.0 )
                     ! ( LETKF Desroziers statistics for 1 week in June 2012,
                     !   for TEMP, aircraft, wind profiler + radar VAD together:
                     !                  1.95 2.00 1.88 1.91 2.04 2.48 2.48 2.48
                     !                  2.48 2.48 2.48 2.48 2.48 2.48 2.48 )
                     !
                     ! (root of) radiosonde (TEMP, PILOT) height error variance
                     !           -------------------------------
    oezsond (nerlev) = (/ 4.3_wp,   4.4_wp,   5.2_wp,   8.4_wp,&
                          9.8_wp,  10.7_wp,  11.8_wp,  13.2_wp,&
                         15.2_wp,  18.1_wp,  19.5_wp,  22.5_wp,&
                         25.0_wp,  32.0_wp,  40.0_wp /)           ,&
                     ! ( values for European Radiosondes in O.I.: 1.5*table )
                     ! ( values up to LM 2.12                   : 1.5*table )
                     !
                     ! (root of) radiosonde temperature error variance
                     !           ----------------------
    oetsond (nerlev) = (/ 1.2_wp,   1.0_wp,   0.7_wp,   0.4_wp,&
                          0.4_wp,   0.5_wp,   0.5_wp,   0.6_wp,&
                          0.7_wp,   0.8_wp,   0.8_wp,   0.9_wp,&
                          0.9_wp,   1.0_wp,   1.2_wp /)           ,&
                     ! ( ECMWF: 1.7, 1.5, 1.3, 1.2, 1.2, 1.4, 1.5, 1.5,
                     !          1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.5 )
                     ! ( LETKF Desroziers statistics for 1 week in June 2012,
                     !   for TEMP, aircraft and RASS together:
                     !          1.10 0.85 0.67 0.55 0.54 0.56 0.56 0.56
                     !          0.56 0.56 0.56 0.56 0.56 0.56 0.56 )
                     !
                     ! (root of) radiosonde relative humidity error variance
                     !           ----------------------------
                     ! ( LETKF Desroziers statistics for 1 week in June 2012:
                     !          0.09 0.13 0.12 0.13 0.13 0.14 0.14 0.14
                     !          0.14 0.14 0.14 0.14 0.14 0.14 0.14 )
                     !
                     ! (root of) AIREP wind error variance
                     !           ----------
    oeairep (nerlev) = (/ 2.5_wp,   2.5_wp,   3.0_wp,   3.5_wp,&
                          4.0_wp,   4.0_wp,   4.0_wp,   4.0_wp,&
                          4.0_wp,   4.0_wp,   4.0_wp,   4.0_wp,&
                          4.0_wp,   4.0_wp,   4.0_wp /)           ,&
                     ! (val. up to LM 2.12: 3.0 below 500 hPa, 4.0 above 300hPa)
                     !
                     ! (root of) AIREP temperature error variance
                     !           -----------------
    oetairp (nerlev) = (/ 1.2_wp,   1.0_wp,   0.7_wp,   0.5_wp,&
                          0.5_wp,   0.6_wp,   0.6_wp,   0.7_wp,&
                          0.8_wp,   0.9_wp,   1.0_wp,   1.1_wp,&
                          1.1_wp,   1.2_wp,   1.4_wp /)
                     ! ( ECMWF        : 1.4, 1.3, 1.2, 1.2, 1.2, 1.3, 1.3, 1.4,
                     !                  1.4, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2 )
                     ! ( up to LM 2.12: 1.4, 1.2, 0.9, 0.8, 0.9, 0.9, 1.0, 1.1,
                     !                  1.2, 1.3, 1.3, 1.4, 1.4, 1.5, 1.7 )

  REAL (KIND=wp)           , PARAMETER  :: &
                     !
                     ! (root of) SYNOP wind error variance
                     !           ----------
    oevsynp (nerlev) = (/ 3.6_wp,   3.6_wp,   5.8_wp,   6.8_wp,&
                          7.8_wp,   9.8_wp,  11.0_wp,  11.8_wp,&
                         11.8_wp,  11.8_wp,  11.8_wp,  11.8_wp,&
                         11.8_wp,  11.8_wp,  11.8_wp /)           ,&
                     ! ( ECMWF: 3.0, 3.0, 3.0, 3.4, 3.6, 3.8, 3.2, 3.2,
                     !          2.4, 2.2, 2.0, 2.0, 2.0, 2.5, 3.0
                     !          used for sea stations only )
                     !
                     ! (root of) SYNOP height error variance (land)
                     !           ------------
    oezsynp (nerlev) = (/ 7.0_wp,   8.0_wp,   8.6_wp,  12.1_wp,&
                         14.9_wp,  18.8_wp,  25.4_wp,  27.7_wp,&
                         32.4_wp,  39.4_wp,  50.3_wp,  59.3_wp,&
                         69.8_wp,  96.0_wp, 114.2_wp /)           ,&
                     !
                     ! (root of) GPS height error variance (land)
                     !           ----------
                     !                       (values not yet specific for gps) !
    oezgps  (nerlev) = (/ 7.0_wp,   8.0_wp,   8.6_wp,  12.1_wp,&
                         14.9_wp,  18.8_wp,  25.4_wp,  27.7_wp,&
                         32.4_wp,  39.4_wp,  50.3_wp,  59.3_wp,&
                         69.8_wp,  96.0_wp, 114.2_wp /)           ,&
                     !
                     ! (root of) SATOB wind error variance
                     !           ----------
    oesatob (nerlev) = (/ 7.0_wp,   7.0_wp,   7.0_wp,   7.0_wp,&
                          7.0_wp,   7.0_wp,   7.0_wp,   7.0_wp,&
                          7.0_wp,   7.0_wp,   7.0_wp,   7.0_wp,&
                          7.0_wp,   7.0_wp,   7.0_wp /)           ,&
                     !
    oevscat =  7.0_wp  ,& ! (root of) SCATT wind   error variance
    oevdrib =  5.4_wp  ,& ! (root of) DRIBU wind   error variance
    oezdrib = 14.0_wp  ,& ! (root of) DRIBU height error variance
    oezship = 14.0_wp  ,& ! (root of) SHIP (sea SYNOP) height error variance
    rherr1  = 10.0_wp  ,& ! (root of) fixed    / normal conditions
    rherr2  = 15.0_wp  ,& ! relative humidity <  if temperature below 233K
    rherr3  = 20.0_wp  ,& ! error variances    \ if rel. humidity below 20%
    oezpd   = 0.003_wp    ! GPS Zenith Path Delay default error variance
 

!-------------------------------------------------------------------------------
! Section 5 : Different kinds of limits
!-------------------------------------------------------------------------------

!         5.0    Horizontal assignment limits (moved to 'obs_assign_gridpt'
!                ----------------------------
!   rtempmx =  1.414_wp ,& ! horizontal search radius to assign a TEMP /
!                              !   PILOT to a grid pt. (if (rtempmx < 0), search
!                              !   is limited to the 4 neighbouring grid pts.)
!   rsypmx  =  1.414_wp ,& ! as 'rtempmx', but for surface-level reports

!         5.1    Redundancy check limits
!                -----------------------
!               (The values for 'rhzlim', 'rhzlair' are replaced by a value [1m]
!                smaller than the minimum grid length used, since the horizontal
!                redundancy check is replaced by a check for common grid points
!                assigned to the reports.)

  REAL (KIND=wp)           , PARAMETER  :: &
    rtmlim  = .150_wp ,& !  time limit, all reports except TEMP, AIREP [hrs]
    rtmlrs  = .751_wp ,& !  time limit for radiosondes (TEMP, PILOT)   [hrs]
    rtmlsy  = .350_wp ,& !  time limit for SYNOP (> 20 min)            [hrs]
    rtmlair = .250_wp ,& !  time limit for reports of obs type 'AIREP' [hrs]
                             !   (time of lowest level of multi-level ODR)
    rhzlim  =  0.2_wp ,& !  horiz.  dist. limit for obs sta. (\ AIREP)  [km]
    rhzlshp =  9.9_wp ,& !  horiz.  dist. limit for ship/buoy/SYNOP     [km]
    rhzlair = .001_wp ,& !  horizont. distance limit for AIREP reports  [km]
    rvtlim  = 10.0_wp ,& !  vertic. dist. limit for obs sta. (\ AIREP)   [m]
    rvtlsy  = 99.0_wp ,& !  vertic. dist. limit (sta.ht.) for SYNOP      [m]
    rvtlshp =199.0_wp ,& !  vertic. dist. limit (sta.ht.) for ship/buoy  [m]
    rvtlair =  4.9_wp ,& !  vertical  distance limit for AIREP reports [hpa]
    rdplim  = 24.9_wp ,& !  vertic. dist. limit within multi-level rep.[hpa]
    rprlim  =  1.1_wp    !  vertic. dist. limit for replacing 'missing
                             !    data' within multi-level reports,        [hpa]
                             !    and for removing colocated levels

!         5.2    Temperature / humidity / pressure / height / fog limits
!                -------------------------------------------------------

  REAL (KIND=wp)           , PARAMETER  :: &
    rttlim =    233._wp ,& ! temperature limit below which rel. humidity
                               !   observation error is increased
    rrhlim =     0.2_wp ,& ! relative humidity limit below which its
                               !   observation error is increased
    rhtsat  =   0.96_wp ,& ! relative humidity threshold beyond which
                               !   observed humidity it set to saturation
    rtshlm =    1.02_wp ,& ! gross error upper limit for relative humidity
    rerrpp = 106000._wp ,& ! msl pressure above which observed pressure is
                               !   assumed to be erroneous
    pminsigt= 10000._wp ,& ! significant TEMP/PILOT level neglected unless
    pminsigv= 20000._wp ,& !   p-obs >= pminsigt and temperat. obs exists or
                               !   p-obs >= pminsigv and wind obs exists
    pminsig = 10000._wp ,& ! AOF read only: neglect significant-level TEMP /
                               !   PILOT data at pressure <= 'pminsig' (100 hPa)
    pqmin   = 29999._wp ,& ! pressure [pa] of level above which moisture
                               !   data are not used
    rpplim  =   900._wp ,& ! pressure level obove which obs. are not used
    vfoglim =   500._wp ,& ! visibility threshold [m] below which the
                               !   existence of low cloud (fog) is assumed
                               !   in the presence of precipitation
    fflim_vad  =  3._wp    ! lower limit for accepting VAD wind speed
 
!         5.3    Height / pressure limits for reporting practice check
!                -----------------------------------------------------

  REAL (KIND=wp)           , PARAMETER  :: &
    rh300  =    300._wp ,& !  300m station height limit applied as:
                               !   - min. height if station reports 900mb level
    rh800  =    800._wp ,& !  800m station height limit applied as:
                               !   - max. height if station reports pmsl
                               !   - max. height if station reports 1000mb level
                               !   - min. allowed if station reports 850mb level
    rh1700 =   1700._wp ,& ! 1700m station height limit applied as:
                               !   - max. height if station reports 900mb level
    rh2300 =   2300._wp ,& ! 2300m station height limit applied as:
                               !   - max. height if station reports 850mb level
                               !   - min. height if station reports 700mb level
    rh3700 =   3700._wp    ! 3700m station height limit applied as:
                               !   - max. height if station reports 700mb level

!         5.4    Limits for the directional wind shear check
!                -------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    nnqcdd  =  7      ! number of levels in the quality control threshold tables

  INTEGER (KIND=iintegers) , PARAMETER  :: &
                      ! thresholds for directional shear, as funct. of speed
    nqcdd  (nnqcdd  ) =          (/  30 , 40 , 50 , 60 , 70 , 80 , 90 /)      ,&
                      ! limit values for sum of wind speeds (for 'nqcdd')
    nqcddff(nnqcdd,2) = RESHAPE( (/  72 , 61 , 57 , 53 , 49 , 46 , 41          &
                                  , 110 , 84 , 77 , 70 , 63 , 52 , 50 /)       &
                               , (/ nnqcdd,2 /) )

  REAL (KIND=wp)           , PARAMETER  :: &
                      ! upper limit of wind speed ('x1'), for which the speed
    qcfddff(2) = (/ & ! shear check is never passed whenever the directional
                      ! shear check is not passed (i.e. for such small wind
      4.5625_wp,& ! speeds, only the speed shear check needs to be done)
      7.8250_wp/) ! (for any directional shear check, the sum of wind speeds
                      !  '(x1 + x2)' must exceed 'nqcddff(?,nnqcdd)'.
                      !  --> conditions for speed shear check
                      !      |x1 - x2| = a + b *(x1 + x2)  , a=20.6  ,  b=0.275
                      !       x1 + x2  = c                 , c=nqcddff(?,nnqcdd)
                      !      ==>  x1 = 0.5* (c*(1-b) - a)  , x1, x2: wind speed)


!-------------------------------------------------------------------------------
! Section 6 : Variables used for the production of aircraft multi-level reports
!-------------------------------------------------------------------------------

!         6.1     Limits
!                 ------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    minslml =  4             ! min. no. of levels in multi-level aircraft report

  REAL (KIND=wp)           , PARAMETER  :: &
    thairt  =   0.25_wp ,& ! maximum temporal distance [h] between the
                               !   lowest report and any single-level report to
                               !   be added to a multi-level AIRCRAFT report
    thairv  =  5500._wp    ! maximum vertical distance [Pa] between two
                               !   successive levels within a multi-level
                               !   AIRCRAFT report

!         6.2     Array sizes
!                 -----------
  
  INTEGER (KIND=iintegers) , PARAMETER , PRIVATE  :: &
    ilstid  =  9      ! character length of the station identity
                      ! Note: Must be equal to 'ilstid' in modules
                      !       'data_obs_record' and 'data_nudge_gather' !!!

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    maxaid  =  150    ! max. number of single-level aircraft reports with the
                      ! same station id.

  INTEGER (KIND=iintegers)  :: &
    nairls            ! length of list 'iairls' (see below)

!         6.3     Temporary ODR and lists
!                 -----------------------
  
  INTEGER (KIND=iintegers) , ALLOCATABLE :: &
    mzmlbd (:,:,:) ,& ! body of temporary aircraft multi-level ODR
    mzmlhd   (:,:) ,& ! header of temporary aircraft multi-level ODR
    mzsgbd   (:,:) ,& ! body of temporary aircraft single-level ODR
    mzsghd   (:,:) ,& ! header of temporary aircraft single-level ODR
    mzslon     (:) ,& ! list for single-level reports, as part of the temporary
                      ! single-level ODR, containing the 'long' list indices
                      ! (see below) of the single-level reports
    mzmlon   (:,:) ,& ! list for multi -level reports, as part of the temporary
                      ! multi -level ODR, containing the 'long' list indices
                      ! (see below) of all single-level reports, which the
                      ! multi-level reports are made of
    iarlls   (:,:) ,& ! 'long' list, containing all aircraft single-level
                      ! reports, that are processed at the current timestep
    iairls (:,:,:) ,& ! other lists with aircraft single-level report meeting
                      ! various requirements
    ismls      (:)    ! list for multi-level reports, that contains the report
                      ! index in the temporary multi-level report array of the
                      ! reports with the current station id

  REAL    (KIND=wp)        , ALLOCATABLE :: &
    zmlbdy (:,:,:) ,& ! body of temporary aircraft multi-level ODR
    zmlhed   (:,:) ,& ! header of temporary aircraft multi-level ODR
    zsgbdy   (:,:) ,& ! body of temporary aircraft single-level ODR
    zsghed   (:,:) ,& ! header of temporary aircraft single-level ODR
    rarlls     (:) ,& ! model layer thickness at reports in 'long' list
    rairls     (:)    ! model layer thickness at reports in other list 1

  CHARACTER (LEN= ilstid ) , ALLOCATABLE :: &
    yarlls     (:) ,& ! station id's in 'long' list
    yairls   (:,:)    ! station id's in other lists
 

!-------------------------------------------------------------------------------
! Section 7 :  For reporting rejection of data: Output buffer, size and formats
!-------------------------------------------------------------------------------

  INTEGER (KIND = iintegers), ALLOCATABLE :: &
    outbuf   (:)    ! buffer containing output for a single node

  INTEGER (KIND=iintegers)  ::  &
    nacout       ,& ! actual number of records stored in the output buffer
    nmxoln       ,& ! maximum length of output buffer for events
    istrej          ! length of strings (station id) in output buffer

  INTEGER (KIND=iintegers) , PARAMETER  :: &
                    ! output format numbers (see 'obs_cdf_print_reject'):
    nfmt1  =  1  ,& ! no pressure
    nfmt2  =  2  ,& ! excess of precipitation
    nfmt3  =  3  ,& ! no accepted data
    nfmt4  =  4  ,& ! excess of levels
    nfmt5  =  5  ,& ! several surface levels
    nfmt6  =  6  ,& ! excess of pressure tendency
    nfmt7  =  7  ,& ! excess of lapse rate
    nfmt8  =  8  ,& ! excess of wind speed shear
    nfmt9  =  9  ,& ! excess of directional shear
    nfmt10 = 10  ,& ! redundancy of surface-level report
    nfmt11 = 11  ,& ! redundancy of multi-level report
    nfmt12 = 12  ,& ! redundancy of aircraft report
    nfmt13 = 13  ,& ! redundancy of wind
    nfmt14 = 14  ,& ! redundancy of temperature
    nfmt15 = 15  ,& ! redundancy of humidity
    nfmt16 = 16  ,& ! redundancy of pressure / height
    nfmt17 = 17  ,& ! thinning of aircraft reports
    nfmt18 = 18  ,& ! exaggerated flight colocation
    nfmt19 = 19  ,& ! flight track error
    nfmt20 = 20  ,& ! message only: fog and precipitation
    nfmt21 = 21  ,& ! message only: fog and invisible sky
    nfmt22 = 22  ,& ! message only: fog and no cloud base
    nfmt23 = 23  ,& ! message only: cloud and no cloud base or fog
    nfmt24 = 24  ,& ! report (partly) blacklisted
    nfmt25 = 25  ,& ! report not on whitelist
    nfmt26 = 26     ! suspicious aircraft identity
 

!-------------------------------------------------------------------------------
! Section 8 : Temporary global model fields
!-------------------------------------------------------------------------------

  REAL    (KIND = wp)       , ALLOCATABLE :: &
    hsurf_tot (:,:) ,& ! total array of model surface height
    fland_tot (:,:)    ! total array of fraction of land in each grid element
 
 
!-------------------------------------------------------------------------------

  END MODULE data_obs_cdfin
