!+ Source module for the observation processing in the data assimilation mode
!-------------------------------------------------------------------------------

MODULE src_obs_proc_cdf

!-------------------------------------------------------------------------------
! Description:
!   This module organizes the reading and observation processing of observation
!   reports from NetCDF observation input files.
!   Special tasks are:
!    - determining from which time interval observations of which type should
!      be read
!    - preparing and calling the interface to a collection of modules (called
!      'library') which do not explicitly use variables from the model
!      environment (e.g. model fields, namelist parameters, etc.) and which
!      actually do the reading and observation processing of the reports from
!      NetCDF input files; the preparation includes the setting of all the
!      required variables which depend explicitly on namelist parameters
!    - calling routines from the library which
!      - do reading the reports from NetCDF files, selecting, assigning them
!        to model grid points and distributing them to the nodes (precessors)
!        according to the sub-domain which they lie on, etc.
!      - derive piecewise vertical profiles from single-level aircraft reports
!      - do multi-level gross error checking (wind shear, stability)
!      - print some diagnostic (ASCII) output
!
! Method:
!   This module contains the following module procedures:
!    - obs_org_cdf_proc        : organize reading obs from NetCDF input files
!    - obs_cdf_prep_interface  : prepare interface to NetCDF obs input file read
!    - obs_cdf_prep_cma_rules  : prepare rules for the use of obs and code types
!    - obs_cdf_mark_old_rep    : flag 'old' report to be deleted in the ODR
!    - obs_cdf_print_caution   : print alerting messages
!    - obs_cdf_print_caut_neff : print alerting messages
!    - obs_cdf_print_eventwarn : print alerting messages
!    - interface_1dvar         : 1dvar interface to COSMO model environment
!    - obs_1dvar_mark_old_rep  : delete 'old' 1dvar reports in the ODR
!    - obs_read_gps            : read + pre-process GPS reps. from COST ASCII f.
!    - obs_gps_init            : check existence + open COST ASCII GPS obs file
!   Driver routine 'obs_org_cdf_proc' is called in 'organize_nudging' of module
!   src_obs_use_org.f90. All other procedures are called in 'obs_org_cdf_proc'.
!
!   It uses from:
!    - src_obs_cdfin_org:   - obs_cdf_read_org
!                           - obs_cdf_interface
!                           - obs_cdf_proc_init
!                           - obs_cdfin_open_files
!                           - obs_fof_check_files
!                           - obs_copy_err_status
!                           - obs_mult_qualicheck
!    - src_obs_cdfin_print: - obs_cdf_print_statist
!                           - obs_cdf_print_events
!                           - obs_cdf_print_reject
!                           - obs_cdf_print_odr
!                           - obs_print_number_rep
!    - src_obs_cdfin_blk:   - obs_read_blacklist
!    - src_obs_proc_air:    - obs_air_org_mult
!  \ - src_obs_1dvar_org:   - organize_1dvar
!    - src_obs_cdfin_util:  - f_p2dp
!                           - obs_assign_gridpt
!                           - get_global_surface
!    - environment:         - exchg_boundaries, comm_barrier, model_abort
!    - time_utilities:      - get_timings
!    - parallel_utilities:  - global_values
!                           - scatterv_values
!                           - distribute_values
!    - utilities:           - convert_month
!                           - phi2phirot, rla2rlarot
!                           - to_upper
!                           - diff_minutes
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V4_5         2008/09/10 Christoph Schraff
!  Initial release
! V4_7         2008/12/12 Christoph Schraff
!  Bug corrections to use the whitelist, on control output for white- and black-
!  listed reports, and at distributing values from multi-level AMDAR reports.
!  Check for consistency between reported station pressure and reduced pressure
!  introduced (if inconsistent then use reduced instead of station pressure).
! V4_9         2009/07/16 Ulrich Schaettler, Kathleen Helmert
!  Use full path for opening blacklist file
!  Replaced some quotes with double quotes (problems with IBM preprocessor)
!  Non-consideration of small VAD-winds (Kathleen Helmert)
! V4_10        2009/09/11 Davide Cesari
!  Add characters after a backslash in comments; g95 interprets that as
!  continuation line
! V4_12        2010/05/11 Ulrich Schaettler, Oli Fuhrer
!  Renamed t0 to t0_melt because of conflicting names
!  Some bug fixes
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_14        2010/06/14 Christoph Schraff
!  Bug fix to avoid model crash at processing of mobile TEMPs during COPS.
!  Bug fix for writing correctly the aircraft flight phase to the ODR.
! V4_15        2010/11/19 Klaus Stephan
!  Call of new sub routine "obs_cdf_raso_rh_bias" for correction for
!  Vaisala RS92 radiosonde humidity.
!  The type of bias correction can be selected by new namelist param. mqcorr92.
! V4_21        2011/12/06 Oliver Fuhrer
!  Correct print statements in SR obs_cdf_read_common (around lines 3353,3356)
!  Comments: Variable ncc in line 5566 is uninitialized (obs_cdf_store_multilev)
!            Variable kobtyp2 in line 8806 is uninitialized (obs_cdf_redundancy)
! V4_22        2012/01/31 Christoph Schraff
!  - New modules src_obs_cdfin_org.f90, src_obs_cdfin_blk.f90,
!    src_obs_cdfin_comhead.f90, src_obs_cdfin_mult.f90, src_obs_cdfin_sing.f90,
!    src_obs_cdfin_print.f90, src_obs_cdfin_util.f90, src_obs_print_vof.f90
!    extracted from this module. In the above modules + in src_obs_proc_air.f90,
!    the only data modules now used are data_parameters.f90, data_obs_cdfin.f90,
!    data_obs_lib_cosmo.f90, and data_obs_record.f90, and the additional modules
!    used are environment.f90, parallel_utilities.f90, utilities.f90 and netcdf. 
!    To allow for this, the interface to these modules has been completely
!    revised by extending subroutine argument lists and by introducing routines
!    obs_cdf_prep_interface, obs_cdf_prep_cma_rules, obs_cdf_print_caution, 
!    obs_cdf_print_eventwarn, and obs_cdf_print_caut_neff.
!  - Minor bug fix: Obs are read by (tconbox-dt)/2 beyond nstop because the
!    last increments may be valid for a time as late as nstop+(tconbox-dt)/2.
!  - Karolin Eichler: use of GPS data introduced, discriminated by processing
!    centre.
!  - New routines for interfacing a 1DVAR pre-processing for satellite radiance
!    data. The interface allows to limit the use of data modules to those also
!    used for the reading and pre-processing of conventional obs from NetCDF
!    files, except for the additional use of 'data_1dvar.f90'.
!    However, the 1DVAR pre-processing is NOT included in this official COSMO
!    version - an experimental version based on V4_18 is available from
!    christoph.schraff_at_dwd.de .
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer, Ulrich Schaettler
!  Replaced qx-variables by using them from the tracer module
!  Introduced nexch_tag for MPI boundary exchange tag to replace ntstep
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  Introduced MESSy interface
! V4_28        2013/07/12 Christoph Schraff
!  Optional capability added for reading obs from feedobs files.
!  Bug fix in routine 'obs_cdf_prep_cma_rules': set GPS obs type in 'cma'.  
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
! V5_1         2014-11-28 Christoph Schraff, Oliver Fuhrer
!  Extensions for reading / processing Mode-S aircraft observations.
!  Replaced ireals by wp (working precision) (OF)
! V5_4         2016-03-10 Christoph Schraff
!  - Removal of AOF interface as alternative file format for observation input.
!  - Reading and pre-processing of GPS-IWV data from COST-716 (V2.0 or V1.0)
!    ASCII files moved here from former module 'src_obs_processing'.
!    For these GPS data, a proper check for the station height (incl. its
!    difference to the model orography) is introduced, and for writing to VOF,
!    the passive bit is set in station characteristics if height flag is set.
!    Diagnostic output on processed obs added to standard output.
!  - Dimension of 'neventr' and 'neventd' reduced from 3 to 2.
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
!   irealgrib    ! KIND-type parameter for real variables in the grib library

!-------------------------------------------------------------------------------

USE data_obs_lib_cosmo, ONLY :  &

    c0         ,& ! standard real constant 0.0
    c1         ,& ! standard real constant 1.0
    c2         ,& ! standard real constant 2.0
    c05        ,& ! standard real constant 0.5
    c3600      ,& ! standard real constant 3600.0
    epsy       ,& ! = 1.E-8_wp : commonly used very small value > 0
    i1         ,& ! standard integer constant 1

! 2b. Pointers for arrays used for reading/processing from NetCDF files
! ---------------------------------------------------------------------

    i_subpos   ,& ! positions of the subdomains in the total domain
    i_gpscen   ,& ! array of processing centres of GPS reports used actively
    r_p        ,& ! pressure (at main levels)
    r_hhl      ,& ! height   (at half levels)
    r_t_ll     ,& ! temperature at lowest model (main) level
    r_ps       ,& ! surface pressure
    r_frland   ,& ! land fraction

! 2c. Pointers for arrays used in 1DVAR
! -------------------------------------

    r_t        ,& ! temperature                                   (  K  )
    r_qv       ,& ! specific water vapor content                  (kg/kg)
    r_t_g      ,& ! weighted surface temperature                  (  K  )
    r_t_2m     ,& ! temperature in 2m                             (  K  )
    r_qv_2m    ,& ! specific water vapor content in 2m            (kg/kg)
    r_u_10m    ,& ! zonal wind in 10m                             ( m/s )
    r_v_10m    ,& ! meridional wind in 10m                        ( m/s )

! 4. I/O device numbers for obs processing / nudging and file names
! -----------------------------------------------------------------

    ! file names
    yucautn    ,& ! caution messages if too many obs for ODR size
    yurejct    ,& ! direct reporting of rejected obs. reports
    yuobsdr    ,& ! observations stored in the observation data record
    yuprint    ,& ! all the remaining information
    ! device numbers
    nucautn    ,& ! caution messages if too many obs for current ODR size
    nustat     ,& ! statistics of processed reports
    nurej      ,& ! direct reporting of rejected obs. reports
    nuodr      ,& ! observations stored in the observation data record
    nupr       ,& ! all the remaining information 
    lopen_odr  ,& ! .true. if yuobsdr is open
    lopen_rej     ! .true. if yurejct is open

USE data_obs_lib_cosmo, ONLY :  &

! 5. CMA observation type and code type numbers
! ---------------------------------------------

    nsynop     ,& ! SYNOP report
    nairep     ,& ! AIREP report (all aircraft reports)
    nsatob     ,& ! SATOB report
    ndribu     ,& ! DRIBU report
    ntemp      ,& ! TEMP  report
    npilot     ,& ! PILOT report
    nsatem     ,& ! SATEM report
    nsattv     ,& ! SATEM report
    ngps       ,& ! GPS   report
    nscatt     ,& ! SCATT report (from NetCDF only)
    nsrscd     ,& !   synop surface report
    natscd     ,& !   automatic synop surface report
    nshscd     ,& !   ship synop report
    nabscd     ,& !   ship synop abbreviated report
    nshred     ,& !   shred report
    natshs     ,& !   automatic ship synop report
    nmetar     ,& !   Metar
    naircd     ,& !   aircraft report
    ncodar     ,& !   codar report
    ncolba     ,& !   colba report
    namdar     ,& !   amdar report
    nacar      ,& !   acar  report
    nmodes     ,& !   mode-s report
    nstbcd     ,& !   satob report
    nhrvis     ,& !   high-res. VIS wind report
    namv       ,& !   AMV   report
    nsst       ,& !   sst report
    ndrbcd     ,& !   dribu report
    nbathy     ,& !   bathy report
    ntesac        !   tesac report

USE data_obs_lib_cosmo, ONLY :  &
    nldtcd     ,& !   temp land   report
    nshtcd     ,& !   temp ship   report
    nmotcd     ,& !   temp mobile report
    ntdrop     ,& !   temp drop   report
    nrocob     ,& !   rocob      report
    nrocsh     ,& !   rocob ship report
    nldpcd     ,& !   pilot land   report
    nshpcd     ,& !   pilot ship   report
    nmopcd     ,& !   pilot mobile report
    nwp_eu     ,& !   European wind profiler report
    nra_eu     ,& !   European SODAR/RASS report
    nravad     ,& !   Radar VAD wind report
    npr_us     ,& !   US Wind Profiler/RASS report
    nstmcd     ,& !   satem report
    nstovs     ,& !   high resolution ATOVS satellite data
    nsmsg1     ,& !   MSG_1 (METEOSAT-9) satellite retrieval
    nsmsg2     ,& !   MSG_2 (METEOSAT-9) satellite retrieval
    nnoa15     ,& !   NOAA15 satellite retrieval  (nnoa15 == 200 + idnoaa15)
    nnoa16     ,& !   NOAA16 satellite retrieval  (nnoa16 == 200 + idnoaa16)
    nnoa17     ,& !   NOAA17 satellite retrieval  (nnoa17 == 200 + idnoaa17)
    nnoa18     ,& !   NOAA18 satellite retrieval
    nascat     ,& !   ASCAT scatterometer report
    nqscat        !   QuickScat scatterometer report
!   ngpgfz     ,& !   GPS report processed by GFZ
!   ngpasi     ,& !   GPS report processed by ASI
!   ngpbkg     ,& !   GPS report processed by BKG
!   ngpgop     ,& !   GPS report processed by GOP
!   ngpieec    ,& !   GPS report processed by IEEC
!   ngpsgn     ,& !   GPS report processed by SGN
!   ngplpt     ,& !   GPS report processed by LPT
!   ngpmet     ,& !   GPS report processed by MET
!   ngprob     ,& !   GPS report processed by ROB
!   ngpige     ,& !   GPS report processed by IGE
!   ngpknmi    ,& !   GPS report processed by KNMI
!   ngpnga        !   GPS report processed by NGA

! 5. CMA obs type and code types
! ------------------------------

USE data_obs_lib_cosmo, ONLY :  &
    n_cma      ,& ! number of CMA obs and code types
    t_cmatyp   ,& ! data type for information on CMA observation and code types
    cma        ,& ! array of meta data on CMA observation and code types
    i_cma      ,& ! function to determine the index of 'cma'
                  ! referring to a given CMA observation and code type
    n_gpc      ,& ! number of GPS processing centres
    t_gpscen   ,& ! data type for information on GPS processing centres
    gpc        ,& ! array of meta data on GPS processing centres
    n_gpc_offs ,& ! index offset of GPS code type elements in array 'cma'
    cs         ,& ! southern latitude limit
    cn         ,& ! northern latitude limit
    cw         ,& ! western longitude limit
    ce            ! eastern longitude limit

! end of data_obs_lib_cosmo

!-------------------------------------------------------------------------------

USE data_obs_cdfin, ONLY :  &

!         1 : NetCDF Observation Input File formats
!             (this section is used ONLY IF obs data are read from NetCDF files)
!         1.1   internal attributes of the different NetCDF input files
!               -------------------------------------------------------

    mxtotin        ,& ! total maximum number of NetCDF input files
    icdfinlen      ,& ! maximum length of NetCDF observation input file name
    ncdf_temp      ,& ! indicator for processing of NetCDF TEMP       input
    ncdf_tempship  ,& ! indicator for processing of NetCDF TEMPSHIP   input
    ncdf_tempdrop  ,& ! indicator for proc. NetCDF TEMP Dropsonde     input
    ncdf_pilot     ,& ! indicator for proc. NetCDF PILOT (z-levels)   input
    ncdf_pilot_p   ,& ! indicator for proc. NetCDF PILOT (p-levels)   input
    ncdf_amdar_ml  ,& ! indicator for proc. NetCDF AMDAR multi-level  input
    ncdf_amdar_vp  ,& ! indicator for proc. NetCDF AMDAR vert.profile input
    ncdf_amdar     ,& ! indicator for proc. NetCDF AMDAR single-level input
    ncdf_acars     ,& ! indicator for proc. NetCDF ACARS single-level input
    ncdf_modes     ,& ! indicator for proc. NetCDF MODE-S KNMI format input
    ncdf_modes_acr ,& ! indicator for proc. NetCDF MODE-S ACARS fmt.  input
    ncdf_wprof     ,& ! indicator for proc. NetCDF wind profiler      input
    ncdf_rass      ,& ! indicator for proc. NetCDF RASS profiler      input
    ncdf_radar_vad ,& ! indicator for proc. NetCDF radar wind prof.   input
    ncdf_synop     ,& ! indicator for proc. NetCDF SYNOP              input
    ncdf_synop_mob ,& ! indicator for proc. NetCDF SYNOP mobile       input
    ncdf_ship      ,& ! indicator for proc. NetCDF SHIP               input
    ncdf_buoy      ,& ! indicator for proc. NetCDF BUOY               input
    ncdf_metar     ,& ! indicator for proc. NetCDF METAR sfc aviation input
    ncdf_gps_zenith,& ! indicator for proc. NetCDF GPS (ZPD / IWV)    input
    ncdf_ascat     ,& ! indicator for proc. NetCDF ASCAT scatterom.   input
    ncdf_qscat     ,& ! indicator for proc. NetCDF QuickScat scatter. input
    ncdf_satob     ,& ! indicator for proc. NetCDF SATOB wind         input
    ncdf_acars_uk  ,& ! indicator for proc. NetCDF ACARS UK + Canada  input
    ncdf_acars_us  ,& ! indicator for proc. NetCDF ACARS US w. humid. input
    ncdf_fof       ,& ! indicator for proc. 'fof' feedobs (feedback) file input
    n_cdfin        ,& ! number of existing NetCDF observation input files
    n_fofin        ,& ! number of existing NetCDF feedobs (feedback) input files
    icdfin         ,& ! obs file type of NetCDF observation input files 
    dimids            ! dimension IDs in NetCDF files

USE data_obs_cdfin, ONLY :  &

!         3.1     Format of event counters
!                 ------------------------
    mxreve     ,& ! length of report event counter array
    nesodr     ,& ! report number exceeding size of ODR (==> adjust Namelist)
    nenoml     ,& ! multi-levl report not made due to ODR array size
    mxdeve     ,& ! length of data event counter array

!         3.2    Event counter arrays
!                --------------------
    neventr    ,& ! counter array of report events
    neventd    ,& ! counter array of data events

!         5.1    Redundancy check limits
!                -----------------------
    rtmlim     ,& ! time limit for all reports except AIREP      [hrs]
    rtmlrs     ,& ! time limit for radiosondes (TEMP, PILOT)     [hrs]
    rtmlair    ,& ! time limit for reports of obs type 'AIREP'   [hrs]
                  !  (time of lowest level of multi-level ODR)

!         7.   For reporting rejection of data: Output buffer, size and formats
!              ----------------------------------------------------------------
    outbuf        ! buffer containing output for a single node

! end of data_obs_cdfin

!-------------------------------------------------------------------------------

USE data_obs_record, ONLY :   &

!       1.1    ODR header format
!              -----------------
!       1.1.1  Header formats of ODR reports: 'omlhed' and 'osghed'
!              ----------------------------------------------------
    nhtime     ,& ! (exact) time of observation in forecast hours
    nh1wte     ,& ! 1dvar retriev.: influence radius of temporal w.: future
    nh1tip     ,& ! 1dvar retriev.: temporal interval for re-doing minimizat.

!       1.1.2  Header formats of ODR reports: 'momlhd' and 'mosghd'
!              ----------------------------------------------------
    nhio       ,& ! (local) x-coord. of grid pt. assigned to obs
    nhjo       ,& ! (local) y-coord. of grid pt. assigned to obs
    nhobtp     ,& ! observation type
    nhcode     ,& ! code type
    nhpass     ,& ! flag for report being set to 'passive'       (see 1.1.4)
    nhqcfw     ,& ! threshold quality control flags for pressure, status of
                  ! verification output

!       1.3    ODR body format
!              ---------------
    nbtp       ,& ! pressure [Pa]   (multi-level reports: 'omlbdy'
    nbsp       ,& ! pressure        (single-level report: 'osgbdy')    [Pa]
    nbgp       ,& ! pressure [Pa]   (GPS reports: 'ogpbdy')
    nbvp       ,& ! pressure [Pa]   (sat retrievals: 'otvbdy' (mandatory))

!       1.5    Further quantities related to ODR
!              ---------------------------------
    ntotml     ,& ! tot. number of stored multi-level reports
    ntotsg     ,& ! tot. number of stored single-level reports
    ntotgp     ,& ! tot. number of stored GPS reports
    ntottv     ,& ! tot. number of stored satellite retrievals

!       2.     Observation data records (ODR)
!       -------------------------------------
    omlbdy     ,& ! body   of multi-level ODR
    omlhed     ,& ! header of multi-level ODR
    osgbdy     ,& ! body   of single-level ODR
    osghed     ,& ! header of single-level ODR
    ogpbdy     ,& ! body   of GPS ODR
    ogphed     ,& ! header of GPS ODR
    otvbdy     ,& ! body   of satellite retrieval ODR
    otvhed     ,& ! header of satellite retrieval ODR
    momlhd     ,& ! header of multi-level ODR
    mosghd     ,& ! header of single-level ODR
    mogphd     ,& ! header of GPS ODR
    motvhd     ,& ! header of satellite retrieval ODR
    yomlhd     ,& ! header of multi-level ODR
    yosghd     ,& ! header of single-level ODR
    yogphd     ,& ! header of GPS ODR
    yotvhd        ! header of satellite retrieval ODR

! end of data_obs_record

!-------------------------------------------------------------------------------

USE data_nudge_all, ONLY :   &

! 1. parameters and related variables
! -----------------------------------

    icdfdirlen ,& ! max. length of name of directory where
                  !   NetCDF observation input files reside
    mxgpc      ,& ! max. number of GPS processing centres used
    lwonl      ,& ! .true. for the node (sub-domain) at which the file with the
                  !        unit number 'nupr' is open
                  !       (i.e. the sub-domain where grid pt. (ionl,jonl) lies)
    acthr      ,& ! actual forecast hour


! 2. namelist variables controlling the data assimilation
! -------------------------------------------------------

! 2.1  temporal variables
!      ------------------

    lverif     ,& ! on - off switch for verification
    lverpas    ,& ! .t. : on - off switch for verif. also of passive reports
    l1dvar     ,& ! .f. : on - off switch for 1DVar
    nudgend    ,& ! 0   : end of nudging period in timesteps
    nverend    ,& ! 0   : end of verification period in timesteps
    hversta    ,& ! 0   : start of verification period in 'model integ. hours'
    hverend    ,& ! 0   : end of verification period in 'model integr. hours'
    tconbox    ,& ! 6*dt: timestep [s] for computing analysis increments
                  !       (i.e. time box of constant analysis increments)
    tipolmx    ,& ! max. time span  for linear interpolat. for upper-air data
    tipmxsu    ,& ! max. time span  for linear interpolat.of surf.-level data
    wtukrsa    ,& ! temporal radius of infl. towards the past for TEMP/PILOT
    wtukrse    ,& ! temporal radius of infl. towards the future for TEMP/PILOT
    wtukara    ,& ! temporal radius of infl. towards the past for aircraft obs
    wtukare    ,& ! temporal radius of infl. towards the future for aircr.data
    wtuksua    ,& ! temporal radius of infl. towards the past for surface data
    wtuksue       ! temporal radius of influence towards the future
    
USE data_nudge_all, ONLY :   &

! 5.    Spatial weights : Spreading of observational information
! --------------------------------------------------------------

    lsvcorl    ,& ! adjustment of vertical correlation scales
    msprpar    ,& ! switch specifying the surfaces along which observation
                  ! increments of upper-air data are (primarily) spread
    vcorls     ,& ! square of the vertical correlation scale
    rhinfl     ,& ! constant part of the 'correlation scale of the autore-
                  ! gressive horizontal correlation function' (='COSAC') [km]
    rhvfac     ,& ! multipl. factor to the vertically varying part of 'COSAC'

! 7.    Quality control and quality weights
! -----------------------------------------

    dtqc       ,& ! 720.: timestep (in [s]) for the threshold quality control

! 8.    Observation processing
! ----------------------------  

    irun_osse  ,& !  model run to derive obs values from file yfofin='fof'
    losse_fg   ,& !  1st guess check flag from 'fof' converted to 'dataset flag'
    thairh     ,& !  maximum horizontal distance for combining single level
                  !  aircraft reports to a multi-level report
    ycdfdir    ,& !  directory where NetCDF obs input + blacklist files reside

    doromx     ,& !  cut-off and gaussian radius of height differences 
                  !  between model orography and station height
    altopsu    ,& !  SYNOP obs. above height 'altopsu' are not assimilated
    fperturb   ,& !  factor to obs error variances to define size of random
                  !  perturbations added to the obs (only from yfofin='fof')
    mqcorr92   ,& !  switch for bias correction for Vaisala RS92 sonde humidity
    maxmlo     ,& !      max. number of multi-level reports within total domain
    maxsgo     ,& !      max. number of (surface-level + upper-air single-level
                  !        reports within the total domain
    maxgpo     ,& !      max. number of GPS reports within total domain
    maxtvo     ,& !      max. number of sat retrievals within total domain
    maxmlv     ,& !  100   : max. number of obs levels in multi-level (m-l) ODR
    mxfrep     ,& !   -1   : max. number of reports in NetCDF feedobs file
    mxfobs     ,& !   -1   : max. number of observations in feedobs file
    nolbc      ,& !  number of grid rows at lateral boundaries
                  !  where obs are neglected
    iseed         !    0   : external seed for random number generator

USE data_nudge_all, ONLY :   &

    lsynop     ,& ! .t.    : .t. if SYNOP data is used
    laircf     ,& ! .t.    : .t. if AIREP data is used (aircraft)
    lsatob     ,& ! .false.: .t. if SATOB data is used
    ldribu     ,& ! .t.    : .t. if BUOY  data is used (drifting buoy)
    ltemp      ,& ! .t.    : .t. if TEMP  data is used
    lpilot     ,& ! .t.    : .t. if PILOT data is used
    lsatem     ,& ! .false.: .t. if SATEM data is used
    lgps       ,& ! .false.: .t. if GPS   data is used
    lscatt     ,& ! .t.    : .t. if SCATT data is used (scatterometer)
    lprodr     ,& ! .t. for diagnostic print of obs data records ODR
    lcd011     ,& ! .t.    : .t. if synop code  11 data is used (land synop)
    lcd014     ,& ! .t.    : .t. if synop code  14 data is used (automatic)
    lcd021     ,& ! .t.    : .t. if synop code  21 data is used (ship)
    lcd022     ,& ! .t.    : .t. if synop code  22 data is used (ship abbrev.)
    lcd023     ,& ! .t.    : .t. if synop code  23 data is used (shred)
    lcd024     ,& ! .t.    : .t. if synop code  24 data is used (autom. ship)
    lcd140     ,& ! .t.    : .t. if synop code 140 data is used (metar)
    lcd041     ,& ! .t.    : .t. if airep code  41 data is used (codar)
    lcd141     ,& ! .t.    : .t. if airep code 141 data is used (airep)
    lcd241     ,& ! .t.    : .t. if airep code 241 data is used (colba)
    lcd144     ,& ! .t.    : .t. if airep code 144 data is used (amdar)
    lcd146     ,& ! .t.    : .t. if airep code 146 data is used (mode-s)
    lcd244     ,& ! .t.    : .t. if airep code 244 data is used (acars)
    lcd088     ,& ! .t.    : .t. if satob code  88 data is used (satob)
    lcd188     ,& ! .false.: .t. if satob code 188 data is used (sst)
    lcd063     ,& ! .t.    : .t. if dribu code  63 data is used (bathy)
    lcd064     ,& ! .t.    : .t. if dribu code  64 data is used (tesac)
    lcd165        ! .t.    : .t. if dribu code 165 data is used (drift. buoy)

USE data_nudge_all, ONLY :   &

    lcd035     ,& ! .t.    : .t. if temp  code  35 data is used (land temp)
    lcd036     ,& ! .t.    : .t. if temp  code  36 data is used (temp ship)
    lcd037     ,& ! .t.    : .t. if temp  code  37 data is used (mobile)
    lcd135     ,& ! .t.    : .t. if temp  code 135 data is used (dropsonde)
    lcd039     ,& ! .t.    : .t. if temp  code  39 data is used (rocob)
    lcd040     ,& ! .t.    : .t. if temp  code  40 data is used (rocob ship)
    lcd032     ,& ! .t.    : .t. if pilot code  32 data is used (land pilot)
    lcd033     ,& ! .t.    : .t. if pilot code  33 data is used (pilot ship)
    lcd038     ,& ! .t.    : .t. if pilot code  38 data is used (mobile)
    lcd132     ,& ! .t.    : .t. if pilot code 132 data is used (win-prof eu)
    lcd133     ,& ! .t.    : .t. if pilot code 133 data is used (sod/rass eu)
    lcd136     ,& ! .t.    : .t. if pilot code 136 data is used (pro/rass us)
    lcd137     ,& ! .t.    : .t. if pilot code 137 data is used (Radar VAD)
    lcd086     ,& ! .t.    : .t. if satem code  86 data is used (satem)
    lcd186     ,& ! .t.    : .t. if atovs code 186 data is used (hi-res ATOVS)
    lcd122     ,& ! .t.    : .t. if scatt code 122 data is used (QuickScat)
    lcd123     ,& ! .t.    : .t. if scatt code 123 data is used (ASCAT)
!   lcd096     ,& ! .t.    : .t. if gps data from COST ASCII file is used
    igpscen    ,& ! array of GPS processing centres used actively
    mcdmsg1    ,& !  processing / use of MSG1   code  71 data
    mcdmsg2    ,& !  processing / use of MSG2   code  72 data
    mcdno15    ,& !  processing / use of NOAA15 code 206 data
    mcdno16    ,& !  processing / use of NOAA16 code 207 data
    mcdno17    ,& !  processing / use of NOAA17 code 208 data
    mcdno18       !  processing / use of NOAA18 code 209 data

USE data_nudge_all, ONLY :   &

    obnlat     ,& !  northern boundary of observation area
    obslat     ,& !  southern boundary of observation area
    obwlon     ,& !  western boundary of observation area
    obelon     ,& !  eastern boundary of observation area
    exnlat     ,& !  northern boundary for exclusion area
    exslat     ,& !  southern boundary for exclusion area
    exwlon     ,& !  western boundary for exclusion area
    exelon     ,& !  eastern boundary for exclusion area

    ionl       ,& !  / grid point coordinates
    jonl       ,& !  \ for standard output
    ionl2      ,& !  / 2nd grid pt coordinates
    jonl2      ,& !  \ for other standard output

!      3.2  Device numbers
!           --------------
!   nugps      ,& ! GPS observation file unit

!      5.1  Switch for observation processing
!           ---------------------------------
!   lexigps    ,& ! .TRUE if GPS input file exists
    tmaxbox    ,& ! maximum interval [s] for computing analysis increments
    cqcwbox       ! (maximum) half interval [s] within which observations are
                  ! quality controlled and written on the feedobs file

! end of data_nudge_all

!-------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 1. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------

    ie           ,& ! number of grid points in zonal direction
    je           ,& ! number of grid points in meridional direction
    ke           ,& ! number of grid points in vertical direction
    ie_tot       ,& ! number of grid points in zonal direction
    je_tot       ,& ! number of grid points in meridional direction
    jstartpar,    & ! start index for computations in the parallel program
    jendpar,      & ! end index for computations in the parallel program

! 4. constants for the horizontal rotated grid and related variables
! ------------------------------------------------------------------

    pollon       ,& ! longitude of the rotated north pole (in degrees, E>0)
    pollat       ,& ! latitude of the rotated north pole (in degrees, N>0)
    polgam       ,& ! angle between the north poles of the systems
    dlon         ,& ! grid point distance in zonal direction (in degrees)
    dlat         ,& ! grid point distance in meridional direction (in degrees)
    startlat_tot ,& ! transformed latitude of the lower left grid point
                    ! of the total domain (in degrees, N>0)
    startlon_tot ,& ! transformed longitude of the lower left grid point
                    ! of the total domain (in degrees, E>0)
    degrad       ,& ! factor for transforming degree to rad

! 5. variables for the time discretization and related variables
! --------------------------------------------------------------

    dt           ,& ! long time-step
    dtdeh        ,& ! dt / 3600 seconds

! 8. Organizational variables to handle the COSMO humidity tracers
! ----------------------------------------------------------------
    idt_qv

! end of data_modelconfig

!-------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 2. physical constants and related variables
! -------------------------------------------

    t0_melt      ,& ! melting temperature of ice   
    r_d          ,& ! gas constant for dry air
    rdv          ,& ! r_d / r_v
    rdocp        ,& ! r_d / cp_d
    g            ,& ! acceleration due to gravity

! 3. constants for parametrizations
! ---------------------------------

    b1           ,& ! variables for computing the saturation vapour pressure
    b2w          ,& ! over water (w) and ice (i)
    b2i          ,& !               -- " --
    b3           ,& !               -- " --
    b4w          ,& !               -- " --
    b4i             !               -- " --

! end of data_constants

!-------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------
    dp0          ,& ! reference pressure thickness of layers          (  Pa )
    p0           ,& ! reference pressure at main levels               (  Pa )
    hhl          ,& ! geometical height of half model levels          (  m  )

! 2. external parameter fields                                        (unit)
! ----------------------------
!   hsurf        ,& ! height of surface topography                    ( m   )
    fr_land      ,& ! fraction of land in a grid element              --

! 3. prognostic variables
! -----------------------
    pp           ,& ! deviation from the reference pressure           (  Pa )
    t            ,& ! temperature                                     (  K  )

! 5. fields for surface values and soil model variables               (unit )
! -----------------------------------------------------
    ps           ,& ! surface pressure                                (  Pa )
    t_g          ,& ! weighted surface temperature                    (  K  )

! 7. fields for model output and diagnostics                          (unit)
! ------------------------------------------
    t_2m         ,& ! temperature in 2m                               (  K  )
    qv_2m        ,& ! specific water vapor content in 2m              (kg/kg)
    u_10m        ,& ! zonal wind in 10m                               ( m/s )
    v_10m           ! meridional wind in 10m                          ( m/s )

! end of data_fields

!-------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    nstop        ,& ! last time step of the forecast
    ntstep       ,& ! actual time step
                    ! indices for permutation of three time levels
    nnew         ,& ! corresponds to ntstep + 1

! 3. controlling the physics
! --------------------------
    itype_gscp   ,& ! type of grid-scale precipitaiton physics

! 7. additional control variables
! -------------------------------
    ltime,        & ! detailled timings of the program are given
    lperi_x,      & ! if lartif_data=.TRUE.: periodic boundary conditions in x-dir.
                    ! or with Davies conditions (.FALSE.)
    lperi_y,      & ! if lartif_data=.TRUE.: periodic boundary conditions in y-dir.
                    ! or with Davies conditions (.FALSE.)
    l2dim           ! 2 dimensional runs

! end of data_runcontrol

!-------------------------------------------------------------------------------

USE data_parallel,      ONLY :  &
    num_compute  ,& ! number of compute PEs
    nboundlines  ,& ! number of overlapping boundary lines of the subdomains
    my_cart_id   ,& ! rank of this subdomain in the cartesian communicator
    icomm_cart   ,& ! communicator for the virtual cartesian topology
    isubpos      ,& ! positions of the subdomains in the total domain
                    ! (i- and j-indices of the lower left and upper right grid
                    !  point in the order:  i_ll, j_ll, i_ur, j_ur  ; only the
                    !  domain interior is considered, not the boundary lines)
    imp_reals    ,& ! determines the correct REAL type used for MPI
    imp_integers ,& ! determines the correct INTEGER type used for MPI
    imp_character,& ! determines the correct CHARACTER type used for MPI
                    ! --> for 'exchg_boundaries':
    ncomm_type,   & ! type of communication
    my_cart_neigh,& ! neighbors of this subdomain in the cartesian grid
    nexch_tag,    & ! tag to be used for MPI boundary exchange
                    !  (in calls to exchg_boundaries)
    sendbuf,      & ! sending buffer for boundary exchange:
                    ! 1-4 are used for sending, 5-8 are used for receiving
    isendbuflen     ! length of one column of sendbuf

! end of data_parallel

!-------------------------------------------------------------------------------

USE data_io,            ONLY :  &

    ydate_ini       ! start of the forecast
                    ! yyyymmddhh (year, month, day, hour)

! end of data_io

!-------------------------------------------------------------------------------

 USE environment,              ONLY :  &
    comm_barrier    ,& ! explicit synchronization point
    exchg_boundaries,& ! performs the boundary exchange between
                       ! neighboring processors
    model_abort        ! aborts the program in case of errors

!-------------------------------------------------------------------------------

 USE parallel_utilities,       ONLY :  &
    global_values      ! computes global values by operating on local arrays

!-------------------------------------------------------------------------------

 USE time_utilities,           ONLY :  get_timings, i_barrier_waiting_nud,     &
                                       i_communications_nud, i_obs_processing

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_util,       ONLY :  &
    f_p2dp             ! get approx. model layer thickness at given pressure

!-------------------------------------------------------------------------------

 USE src_obs_proc_air,         ONLY :  & 
    obs_air_org_mult   ! production of multi-level aircraft reports from
                       !   single-level reports

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_print,      ONLY :  & 
    obs_cdf_print_statist,& ! prints summary of statistics on processed reports
    obs_cdf_print_events ,& ! prints summary of report and data events
    obs_cdf_print_reject ,& ! prints messages on report / data rejection
    obs_cdf_print_odr    ,& ! prints a part of the ODR
    obs_print_number_rep    ! prints statistics of processed reports per node

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_org,        ONLY :  & 
    obs_cdf_read_org     ,& ! organizes reading, process. + distribution of obs
    obs_cdf_interface    ,& ! transfers input values from the model environment
                            ! to the observation processing library
    obs_cdf_proc_init    ,& ! allocates long term storage arrays + opens files
    obs_cdfin_open_files ,& ! opens NetCDF observation input files
    obs_fof_check_files  ,& ! checks existence of feedobs (feedback) input files
    obs_cdf_raso_rh_bias ,& ! bias correction of Vaisala RS92 r-sonde humidity
    obs_cdf_mult_qualicheck ! model-independent multi-level gross error checking

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_blk,        ONLY :  &
    obs_read_blacklist   ! read blacklist and whitelist from file and store them

!-------------------------------------------------------------------------------

!USE src_obs_processing,       ONLY :  &
!   obs_read_gps      ,& ! read and pre-process GPS reports from COST ASCII file
!   obs_gps_init         ! check existence and open COST ASCII GPS obs file

!-------------------------------------------------------------------------------

USE src_tracer,                ONLY: trcr_get, trcr_errorstr

!-------------------------------------------------------------------------------

! XXX src_obs_1dvar_org,         ONLY:  &
!   organize_1dvar       ! produces satellite retrievals from radiances by 1DVAR

!-------------------------------------------------------------------------------
! Modules used for additional routines 'obs_read_gps' and 'obs_gps_init':
!-------------------------------------------------------------------------------

USE data_obs_cdfin, ONLY :  &

!         3.1.1   Report event counter array format
!                 ---------------------------------
    netime     ,& ! time out of range (too old) (or time missing)
    nenoal     ,& ! no station altitude
    neloca     ,& ! station location out of range
    nezdif     ,& ! distance 'model orography - station altitude' too large
    neobct     ,& ! observation or code type excluded on area with sta. location
    nesodr     ,& ! report number exceeding size of ODR (==> adjust Namelist)
    neredn     ,& ! redundancy between 2 multi-level, or 2 single-level reports

!         5.1    Redundancy check limits
!                -----------------------
    rhzlim     ,& !  horiz.  dist. limit for obs sta. (\ airep)  [km]
    rvtlim     ,& !  vertic. dist. limit for obs sta. (\ airep)   [m]

!         7.   For reporting rejection of data: Output buffer, size and formats
!              ----------------------------------------------------------------
    nacout     ,& !  actual number of records stored in the output buffer
    nmxoln     ,& !  maximum length of output buffer
    istrej     ,& !  length of strings (station id) in output buffer
    nfmt10     ,& ! redundancy of surface-level report

!         8.   Temporary global model fields
!              -----------------------------
    hsurf_tot  ,& ! total array of model surface height
    fland_tot     ! total array of fraction of land in each grid element

! end of data_obs_cdfin

!-------------------------------------------------------------------------------

USE data_obs_lib_cosmo, ONLY :  &

!         1.   General parameters
!              ------------------
    rmdi       ,& ! =-1.E31_wp : commonly used missing data indicator
    rmdich     ,& ! =-1.E30_wp : commonly used check value for miss data

!         2.   Variables and parameters obtained by calling 'obs_cdf_interface'
!              ----------------------------------------------------------------
    maxgpl     ,& ! size (report dimension) of the (ground-based) GPS  ODR

!         4.   I/O device numbers for obs processing and file names
!              ----------------------------------------------------
    yustats    ,& ! statistics of processed reports
    nustat     ,& ! statistics of processed reports

!         5.   CMA observation type and code type numbers
!              ------------------------------------------
    ngpgfz        !   GPS report processed by GFZ

! end of data_obs_lib_cosmo

!-------------------------------------------------------------------------------

USE data_nudge_all, ONLY :   &

!         2.3  namelist variables, for observation processing
!              ----------------------------------------------
    lcd096     ,& !  .t. if gps   code  96 data is used

!         5.   I/O device numbers and file names for nudging
!              ---------------------------------------------
    nugps      ,& ! GPS observation file unit
    yugps      ,& ! input gps observation file

!         6.   Switch for observation processing and varia
!              -------------------------------------------
    lexigps    ,& ! .TRUE if GPS input file exists
    liwvssc       ! .t. : spatial consistency check of IWV performed

! end of data_nudge_all

!-------------------------------------------------------------------------------

USE data_obs_record, ONLY :   &

!       1.1    ODR header format
!       ------------------------
!       1.1.1  Header formats of ODR reports:'omlhed','osghed','ogphed','otvhed'
!              -----------------------------------------------------------------
    mxghed     ,& ! header length of GPS reports
    nhilon     ,& ! longitude of observing station
    nhjlat     ,& ! latitude  of observing station
    nhalt      ,& ! station altitude [m]
    nhtime     ,& ! time of observat. in forecast hours
    nhsurf     ,& ! height of model grid pt. to which obs. is assigned
    nhzio      ,& ! longitude of obs. station (or lowest datum) in grid pt. unit
    nhzjo      ,& ! latitude  of obs. station in grid pt. units
!   nhsolz     ,& ! solar zenith angle [deg]

!       1.1.2  Header formats of ODR reports:'momlhd','mosghd','mopghd','motvhd'
!              -----------------------------------------------------------------
    mxghdf     ,& ! header length of GPS reports
    nhio       ,& ! (local) x-coord. of grid pt. assigned to obs
    nhjo       ,& ! (local) y-coord. of grid pt. assigned to obs
    nhitot     ,& ! global x-coord. of grid pt. assigned to obs
    nhjtot     ,& ! global y-coord. of grid pt. assigned to obs
    nhobtp     ,& ! observation type
    nhcode     ,& ! code type
    nhschr     ,& ! station characteristics                      (see 1.1.4)
    nhpass     ,& ! flag for report being set to 'passive'       (see 1.1.4)
    nhqcfw     ,& ! threshold quality control flags for pressure, status of
                  ! verification output
    nhflag     ,& ! report flags (obs type, surf., alt., sta ID) (see 1.2.1)
    nhcorr     ,& ! update sequence number (station correction indicator)
    nhcat      ,& ! data     category (from BUFR Section 1) 
    nhcats     ,& ! data sub-category (from BUFR Section 1)
    nhkz       ,& ! DWD internal classification number (observation type)
    nhcent     ,& ! originating centre  +  (1000* sub-centre)
    nhstid     ,& ! station identity number 
    nhdate     ,& ! absolute exact observation date [yyyymmdd]
    nhhrmn     ,& ! absolute exact observation time [hhmm]
    nhsyhr     ,& ! absolute nominal (synoptic) observation time [yymmddhh]
    mxghdf     ,& ! header length of GPS reports

!       1.1.3  Header formats of ODR reports:'yomlhd','yosghd','yopghd','yotvhd'
!              -----------------------------------------------------------------
    ilstid     ,& ! character length of the station identity
    ilstidp       ! char. length used for printing the station ID
                  ! Note: (ilstid >= ilstidg >= ilstidp), cf. data_nudge_gather

USE data_obs_record, ONLY :   &

!       1.2    Bit patterns for packed information in ODR (and VOF) header
!       ------------------------------------------------------------------
    nvpsbp     ,& ! bit pos. for report set passive since at least  nvschr
                  !              1 of the following 5 flags or
                  !              the flight track flag applies
    nvobbp     ,& ! bit pos. for flag: 'station location out of       "
                  !                     user-specified area'
    nvalbp     ,& ! bit pos. for flag: 'distance model orography -    "
                  !                     station altitude too large'
    nvexbp     ,& ! bit pos. for flag: 'observation or code type      "
                  !                     excluded at station location'
    nvrdbp     ,& ! bit pos. for flag: 'redundant report'             "
    nvsebp     ,& ! bit pos. for report located at sea grid pt.       "
    nvinbp     ,& ! bit pos. for instrument specification indicator   "

!       1.3    ODR body format
!              ---------------
!       1.3.5  Body format of ODR of GPS reports: 'ogpbdy'
!              -------------------------------------------
    mxgbdy     ,& ! body length of GPS reports
    nbgtze     ,& ! error in total zenith delay [mm]
    nbgzpd     ,& ! zenith path delay (total zenith delay)             [mm]
    nbgzwd     ,& ! zenith wet delay [mm]
    nbgiwv     ,& ! integrated water vapour [mm]
!   nbgp       ,& ! pressure [Pa]
    nbgt       ,& ! temperature [K]
    nbgrh      ,& ! relative humidity [/]
    nbgbia     ,& ! bias correction to integrated water vapour [mm]
    nbgiwa     ,& ! adjusted (bias corrected) integrated water vapour [mm]

!       1.3.6  Body format of ODR of GPS report flags: 'mogpbd'
!              -------------------------------------------------
    mxgbdf     ,& ! body length of GPS reports
    nbgflg     ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbgerr     ,& ! status flag word        (bit pattern, see below: 'nb?err')
    nbgqcf     ,& ! threshold quality control flags      (see below: 'nb?qcf')
    nbglid     ,& ! level identity          (bit pattern, see below: 'nb?lid')

!       1.4.2  Other bit patt. for packed info in ODR (VOF) body, general words
!              ----------------------------------------------------------------
    nvriwv     ,& ! bit pos. for status 'active' for IWV            nb?err
    nvlidp     ,& ! level id. bit pattern                           nb?lid

!       1.5    Further quantities related to ODR
!              ---------------------------------
    imdi       ,& ! missing data indicator for ODR integers (2^31-1)
    fdoro      ,& ! scaling factor to vertical distances betw. model

!       2.     Observation data records (ODR)
!       -------------------------------------
    mogpbd        ! body of GPS ODR

! end of data_obs_record

!-------------------------------------------------------------------------------

  USE mo_fdbk_tables, ONLY :   &

    FL_OBSTYPE ,& ! passive report type (at obs. location)
    FL_AREA    ,& ! location not in valid area (i.e. outside user-defined area)
    FL_HEIGHT     ! location not in valid height range

! end of mo_fdbk_tables

!-------------------------------------------------------------------------------

 USE parallel_utilities,       ONLY :  &
    scatterv_values, & ! scatters batches of variable length to other nodes
    distribute_values  ! distributes a set of values from one node to all others

!-------------------------------------------------------------------------------

 USE utilities,                ONLY :  &
    convert_month,   & ! number of month converted from a 3-character string
    phi2phirot,      & ! converts phi from the real to the rotated system
    rla2rlarot,      & ! converts lambda from the real to the rotated system
    to_upper,        & ! converts alphabetic characters from lower to upper case
    diff_minutes       ! compute difference in minutes between 2 dates/times

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_util,       ONLY :  &
    obs_assign_gridpt ,& ! horiz. assignment of an obs report to a grid pt.
    get_global_surface   ! get orography / land fraction on global domain

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
!+ Module procedure in "src_obs_proc_cdf" for reading and distributing reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_org_cdf_proc ( aiwthr , nlastproc )

!-------------------------------------------------------------------------------
! Description:
!   This procedure of module "src_obs_proc_cdf" organizes the reading,
!   pre-selection and distribution of observation reports.
!
! Method:
!   This routine is called at each timestep when new analysis increments are
!   computed. In section 2, it is determined whether and for which time interval
!   new observations have to be read from which NetCDF observation input file.
!   Then, the reading is done, the observations are assigned to model grid 
!   points, distributed to the appropriate nodes (sub-domains), and stored in
!   the internal storage arrays ODR.
!   The pre-processing of observations includes redundancy checking, gross
!   error checking, thinning and aggregation of single-level aircraft reports
!   to piecewise vertical profiles, flight track checking, assignment of
!   observation errors, etc.
!   Writing of a range of monitoring / diagnostic (ASCII) output is performed
!   or prepared.
!
! Written by        : DWD, Christoph Schraff
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  REAL (KIND=wp)           , INTENT (IN)  ::  &
    aiwthr              ! model time [hrs] for which the analysis increments are
                        ! computed (for which the temporal weights are valid)

  INTEGER (KIND=iintegers) , INTENT (IN)  ::  &
    nlastproc           ! timestep at which this routine is called the last time


! Local parameters: None
! ----------------

  REAL    (KIND=wp   )     , PARAMETER  ::  &
    thdlnpt = 28.6_wp    ! mean vertical potential temperature gradient in
                             ! the troposphere in  log(p) units

! Local scalars:
! -------------

  LOGICAL                  , SAVE  ::  &
    lobpfrs = .TRUE.  ,& ! variable for 'first time of obs. processing'
    lneedob  (mxtotin)   ! new observation need to be read at current timestep

  INTEGER (KIND=iintegers) , SAVE  ::  &
    min_end  (mxtotin),& ! end of reading time interval [min]
    min_last (mxtotin),& ! last obs time [min] to be used ever in this model run
    ncdf_tbox(mxtotin),& ! length of time box [min] to be read
    nexceold             ! number of multi-level reports rejected at this and
                         ! the previous reading time due to ODR size

  REAL (KIND=wp)           , SAVE  ::  &
    cdf_wta (mxtotin),& ! max time range of influence towards the past for obs
    tairobox            ! time box for reading aircraft observations [h]

  REAL (KIND=wp)           , TARGET  ::  &
    zp     (ie,je,ke)   ! model pressure field

  INTEGER (KIND=iintegers) ::  &
    n_totin          ,& ! total number of input files ('cdfin' + feedobs 'fof')
    nodrtot  (3)     ,& ! total number of multi-level / single-level / GPS rep.
    nodrold  (3)     ,& ! number of multi-level / single-level / GPS reports
                        ! before having read new data at the current timestep
                        !   (can be modified in 'obs_cdf_del_old_rep')
    imaxl    (4)     ,& ! size (report dimension) of the (4 types of) ODR
    nexceed  (5)     ,& ! number of reports in excess of ODR array size
    nexcedml         ,& ! number of rejected multi-level reports due to ODR size
    nexceair (2)     ,& ! number of multi-level reports derived from single-
                        !           level reports but in excess of array size
    min_sta (mxtotin),& ! start of reading time interval [min]
    min_box          ,& ! minimum length [min] of time box to be read
    min_need         ,& ! minimum obs time needed to be read now
    icdf             ,& ! obs file type of NetCDF observation input file
    ilcf, irep, ivar ,& ! loop indices
    io   , jo        ,& ! local horizontal coordinates of observation
    nstat, istat        ! error status variables

  REAL (KIND=wp)           ::  &
    cdf_wte (mxtotin),& ! max time range of influence towards the future for obs
    cdf_wtuke        ,& ! time range of influence to the future for a single obs
    cdf_last         ,& ! last obs time [hr] to be used ever in this model run
    ctshift          ,& ! time [hrs] by which an obs may be shifted (in the
                        !   redundancy check / making multi-level aircraft rep)
    rtmlmt           ,& ! colocation threshold for time
    zpob             ,& ! pressure at observation
    vscale       (4) ,& ! scale (in [ln(p)] units) of the static Gaussian
                        ! vertical correlation function (see obs_air_org_mult)
                        ! --> usual and reasonable values: 0.577,0.577, 0.2, 0.2
    tscale           ,& ! scale [hrs] of the temporal Gaussian weight
                        ! --> usual and reasonable value : 0.25
    fsize        (4)    ! ratio of namelist parameter (max??o)
                        !       to ODR size (max??l)

  LOGICAL                  ::  &
    ldo_airmult      ,& ! do call 'obs_air_org_mult' because new air-obs read
    ldo_pr_out          ! do call 'obs_cdf_print_odr' because new obs have been
 
  CHARACTER (LEN=20)       ::  &
    yroutine            ! name of this subroutine

  INTEGER (KIND=iintegers) , SAVE  ::  &
    imaxmll          ,& ! size (report dimension) of the  multi-level (m-l)  ODR
    imaxsgl          ,& ! size (report dimension) of the single-level (s-l)  ODR
    imaxgpl          ,& ! size (report dimension) of the (ground-based) GPS  ODR
    imaxtvl          ,& ! size (report dimension) of the satellite retrieval ODR
    madj_hum            ! mode switch on adjusting observed humidity to model
 
! Local arrays:
! ------------

  REAL (KIND=wp)           , ALLOCATABLE  ::       &
    odp    (:)          ! approx. model layer thickness at the obs level

  INTEGER (KIND=iintegers) , ALLOCATABLE  ::       &
    ievnt (:)           ! index for data event counters
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_org_cdf_proc
!-------------------------------------------------------------------------------

  yroutine = 'obs_org_cdf_proc'

!-------------------------------------------------------------------------------
! Section 0: Initalisation of the (MPI) environment and domain decomposion
!            if not already done
!            (see 'organize_setup' in module src_setup.f90 of the COSMO model;
!            'domain_decomposition', 'check_decomposition', and 'grid_constants'
!            also reside src_setup.f90,
!            while 'init_environment' and 'init_procgrid' are in environment.f90
!            and 'init_par_utilities' in parallel_utilities.f90)
!-------------------------------------------------------------------------------

! CALL init_environment (nproc, my_world_id, icomm_world, igroup_world,       &
!                        imp_integers, imp_reals, imp_grib, imp_byte,         &
!                        imp_character, imp_logical, iexch_req, irealgrib,    &
!                        yzerrmsg, izerror )

! ALLOCATE ( isubpos(0:num_compute-1,4)   , STAT=istat )
! isubpos(:,:) = 0

! ! Initialize the cartesian processor grid
! CALL init_procgrid (                                                        &
!       nproc, nprocx, nprocy, nprocio, lperi, lreproduce, lreorder,          &
!       igroup_world, my_world_id, icomm_compute, icomm_cart,                 &
!       igroup_cart, my_cart_id, my_cart_pos, my_cart_neigh,                  &
!       my_peri_neigh, icomm_row, lcompute_pe, yzerrmsg, izerror)

! IF ( lcompute_pe )  CALL domain_decomposition

! ! Compute constants related to the grid
! CALL grid_constants

! ! Initialize the utility module parallel_utilities
! CALL init_par_utilities                                                     &
!  (ie, je, ke, ie_tot, je_tot, ke_tot, ie_max, je_max,                       &
!   istartpar, iendpar, jstartpar, jendpar, nproc, nprocx, nprocy, nprocio,   &
!   isubpos, nboundlines, icomm_cart, my_cart_id, imp_reals, imp_integers)

! ! In debug-mode, control information about the decomposition is printed
! IF ( (idbg_level > 5) .AND. (num_compute > 1) .AND. (lcompute_pe))           &
!   CALL check_decomposition (yzerrmsg, izerror)

!-------------------------------------------------------------------------------
! Section 1: Interface for the observation pre-processing:
!            make the 'model environment' available to the obs process library
!             (this basically means: copying values from the 'model environment'
!              into variables from data modules that are part of the observation
!              pre-processing + 1DVAR library
!              ('library' denotes here simply a set of modules))
!            Remarks:
!            - The subroutines called here are part of this module, except for
!              'obs_cdf_interface' itself which is part of the obs proc library.
!            - The values of all the other variables which are set outside the
!              obs pre-proc library but which must be available within a certain
!              1DVAR-related subroutine are transfered directly from the 'model'
!              environment through the argument list of the respective routine.
!            - All the other variables which must be commonly available within
!              the obs pre-proc library but can be set within the library are
!              stored in one of the data modules that are part of the library.
!-------------------------------------------------------------------------------

! prepare the rules for the use of CMA observation and code types
! by filling 'cma' of type 't_cmatyp' as a function of namelist input
! -------------------------------------------------------------------

  IF (lobpfrs)  CALL obs_cdf_prep_cma_rules
!               ===========================

! prepare those input variables for 'obs_cdf_interface'
! which are not yet known here:
! - some constant variables depending on namelist input
! - model 3-D pressure field
! -----------------------------------------------------

  CALL obs_cdf_prep_interface ( ie, je, ke , zp , madj_hum                     &
                              , imaxtvl, imaxmll, imaxsgl, imaxgpl )
! ===========================

! make the 'model' environment available to the obs proc library
! --------------------------------------------------------------
  ! interface for scalars: make a copy by calling 'obs_cdf_interface'

  CALL obs_cdf_interface ( ie, je, ke, ie_tot, je_tot                          &
                         , num_compute, nboundlines, my_cart_id, icomm_cart    &
                         , imp_reals, imp_integers, imp_character              &
                         , pollon, pollat, polgam, dlon, dlat                  &
                         , startlat_tot, startlon_tot, degrad                  &
                         , imaxmll, imaxsgl, imaxgpl, imaxtvl, maxmlv          &
                         , nolbc, madj_hum, doromx, altopsu, acthr             &
                         , irun_osse, losse_fg, fperturb, iseed                &
                         , g, t0_melt, r_d, rdv, rdocp                         &
                         , b1, b2w, b2i, b3, b4w, b4i                          &
                         , lverpas, lwonl, ydate_ini, mxgpc )
! ======================

  ! interface for arrays: set 'library 1' pointers to target arrays
  ! (target arrays must have 'target' attribute !)
  ! (this avoids the need to make copies of these arrays, by (de-)allocation)
                                 ! target arrays (must) have dimensions:
  r_p       =>  zp               ! (ie,je,ke)
  r_hhl     =>  hhl              ! (ie,je,ke+1)
  r_t_ll    =>  t(:,:,ke,nnew)   ! (ie,je)
  r_ps      =>  ps(:,:,nnew)     ! (ie,je)
  r_frland  =>  fr_land          ! (ie,je)
  i_subpos  =>  isubpos          ! (0:num_compute-1,4)
  i_gpscen  =>  igpscen          ! (mxgpc)

!-------------------------------------------------------------------------------
! Section 2: Initialise observation processing, including opening of files,
!            reading of blacklist, and allocation of long term storage arrays
!            and of short term storage rejection messaging arrays ('outbuf').
!            (The subroutines called here are part of the obs precess. library.)
!-------------------------------------------------------------------------------

! array allocation

  CALL obs_cdf_proc_init ( lobpfrs , n_cma )
! ======================

  IF (lobpfrs)  THEN

    CALL obs_read_blacklist ( icdfdirlen , ycdfdir )
!   =======================

    n_cdfin = 0
! check how many NetCDF observation input files exist for which obs type;
! get file annexes and open files (and get their file unit numbers)

    CALL obs_cdfin_open_files ( icdfdirlen , ycdfdir )
!   =========================

    n_fofin = 0
! check how many NetCDF feedobs (feedback) input files exist; get file annexes

    CALL obs_fof_check_files ( icdfdirlen , ycdfdir )
!   ========================

! write a CAUTION to YUCAUTN only if not even 1 NetCDF input file
! (or feedobs input file) exists
    IF ((n_cdfin == 0) .AND. (n_fofin == 0) .AND. (my_cart_id == 0)) THEN
      OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'          &
                                 , POSITION='APPEND', IOSTAT=nstat)
      WRITE( nucautn,'("CAUTION !!!!! NO NetCDF OBSERVATION INPUT FILES IN:"   &
                     &,A)' ) ycdfdir(1:LEN_TRIM(ycdfdir))
      PRINT          '("CAUTION !!!!! NO NetCDF OBSERVATION INPUT FILES IN:"   &
                     &,A)' , ycdfdir(1:LEN_TRIM(ycdfdir))
      CLOSE (nucautn)
    ENDIF

! initialise reading of GPS IWV obs from a COST ASCII format file

    CALL obs_gps_init
!   =================

! Set initial number of processed reports to zero
    ntotml = 0
    ntotsg = 0
    ntotgp = 0
    nexceold = 0
  ENDIF

!-------------------------------------------------------------------------------
! Section 3: Determine whether and for which time interval observations have to
!            be read now from which type of NetCDF observation input file
!-------------------------------------------------------------------------------

  IF (lobpfrs) THEN
! minimum length [min] of time box to be read
    min_box = INT( tconbox/60._wp -epsy ) + 1

    DO ilcf = 1 , n_cdfin + n_fofin
      IF (ilcf <= n_cdfin)  icdf  =  icdfin(ilcf)
      IF (ilcf >  n_cdfin)  icdf  =  ncdf_fof
      rtmlmt  =  rtmlim
! nudging time windows which determine how much to the past and the future
! observations are required for exactly the current timestep only
      IF (     (icdf == ncdf_temp ) .OR. (icdf == ncdf_tempship)               &
                                    .OR. (icdf == ncdf_tempdrop)               &
          .OR. (icdf == ncdf_pilot) .OR. (icdf == ncdf_pilot_p )) THEN
        cdf_wta (ilcf) = MAX( wtukrsa , tipolmx )
        cdf_wte (ilcf) = MAX( wtukrse , tipolmx )
        cdf_wtuke      =      wtukrse
        rtmlmt         =      rtmlrs
      ELSEIF (     (icdf == ncdf_amdar)    .OR. (icdf == ncdf_acars)           &
              .OR. (icdf == ncdf_acars_uk) .OR. (icdf == ncdf_acars_us)        &
              .OR. (icdf == ncdf_modes)    .OR. (icdf == ncdf_modes_acr)       &
              .OR. (icdf == ncdf_amdar_ml) .OR. (icdf == ncdf_amdar_vp)) THEN
        cdf_wta (ilcf) =      wtukara
        cdf_wte (ilcf) =      wtukare
        cdf_wtuke      =      wtukare
        rtmlmt         =      rtmlair
      ELSEIF (     (icdf == ncdf_wprof) .OR. (icdf == ncdf_rass)               &
              .OR. (icdf == ncdf_radar_vad)) THEN
!   this must be the same as in Section 1 of routine 'local_sort_reports'
!   in module 'src_sing_local.f90' where the temporal weights are hard-coded
!   for these observation code types !
        cdf_wta (ilcf) = MAX( MIN( 0.5_wp, wtukrsa )                      &
                            , MIN( 0.5_wp, tipolmx ) )
        cdf_wte (ilcf) = MAX( MIN( 0.2_wp, wtukrse )                      &
                            , MIN( 0.5_wp, tipolmx ) )
        cdf_wtuke      =      MIN( 0.2_wp, wtukrse )
      ELSEIF (icdf == ncdf_gps_zenith) THEN
        cdf_wta (ilcf) = MAX( wtuksua , tipmxsu )
        cdf_wte (ilcf) = MAX( wtuksue , tipmxsu )
        cdf_wtuke      =      wtuksue
      ELSEIF (icdf == ncdf_fof) THEN
        ! feedobs file case
        cdf_wta (ilcf) = MAX( wtukrsa , wtukara , wtuksua , tipolmx , tipmxsu )
        cdf_wte (ilcf) = MAX( wtukrse , wtukare , wtuksue , tipolmx , tipmxsu )
        cdf_wtuke      = MAX( wtukrse , wtukare , wtuksue )
      ELSE
        ! surface obs
        cdf_wta (ilcf) = MAX( wtuksua , tipmxsu )
        cdf_wte (ilcf) = MAX( wtuksue , tipmxsu )
        cdf_wtuke      =      wtuksue
      ENDIF

! the obs time box has to encompass at least
! - the (max.) half interval 'cqcwbox' at which the obs are quality controlled
!   and written on the feedobs file !
!   (because 'cqcwbox' has to be used as a condition in 'local_sort_reports')
! - the colocation threshold for time in the redundancy check (since the obs
!   can be shifted in time by this threshold)
      ctshift        = MAX( cqcwbox /c3600 , rtmlmt )
      cdf_wta (ilcf) = MAX( cdf_wta(ilcf) , ctshift )
      cdf_wte (ilcf) = MAX( cdf_wte(ilcf) , ctshift )

! last observation time to be used ever in this model run
!   ( 0.5*tmaxbox is added because the increments that are computed (nearly)
!     at nstop are in fact valid (nearly) for time nstop+(tconbox-dt)/2 )
      cdf_last        = MIN( nstop*dtdeh + cdf_wta(ilcf)  -epsy                &
                                         + c05* tmaxbox /c3600                 &
                           , MAX( nudgend*dtdeh - cdf_wtuke                    &
                                , hverend ) + ctshift +epsy )
!                               , nverend*dtdeh )         +epsy )
      min_last (ilcf) = INT( 60._wp* cdf_last )

! length of time box [min] to be read each time when reading is required
! (this time box must be at least as large as 'tconbox' which is the timestep
!  for computing new analysis increments and for calling this subroutine)
!     IF ((icdf == ncdf_ascat) .OR. (icdf == ncdf_qscat)) THEN
!       ! special real-time treatment for scatterometer is not needed as long as
!       ! pre-processing (choosing 1 out of 4 wind values) is done outside COSMO
!       ncdf_tbox (ilcf) = 1
!     ELSE
        ncdf_tbox (ilcf) = 60
!   time box [h] for aircraft obs
        tairobox         = REAL( ncdf_tbox(ilcf) , wp ) / 60._wp
!     ENDIF
      ncdf_tbox (ilcf) = MAX( ncdf_tbox(ilcf) , min_box )

! time interval of observations which are used for computing the analysis
! increments at the current timestep
      min_sta (ilcf) = INT( MIN( aiwthr - (cdf_wte(ilcf) -epsy)                &
                               , hversta - ctshift )             *60._wp )
      min_end (ilcf) = INT(     (aiwthr + (cdf_wta(ilcf) -epsy)) *60._wp )
! end of reading time interval which is sufficient to compute the analysis
! increments within the next 'ncdf_tbox(ilcf)' minutes of model time
      min_end (ilcf) = MIN( min_end(ilcf) + ncdf_tbox(ilcf) , min_last(ilcf) )
      lneedob (ilcf) = .TRUE.
    ENDDO

  ELSEIF (.NOT. lobpfrs) THEN
    DO ilcf = 1 , n_cdfin + n_fofin
! decide if new observations need to be read
      min_need = INT( (aiwthr + (cdf_wta(ilcf) -epsy)) *60._wp )
      lneedob (ilcf) =       (min_need      >  min_end (ilcf))                 &
                       .AND. (min_need      <= min_last(ilcf)+c05*tconbox)     &
                       .AND. (min_end(ilcf) <  min_last(ilcf))
      IF (lneedob(ilcf)) THEN
! determine time interval of observations to be read now
        min_sta (ilcf) =      min_end(ilcf) + 1
        min_end (ilcf) = MIN( min_end(ilcf) + ncdf_tbox(ilcf) , min_last(ilcf) )
      ENDIF
    ENDDO
  ENDIF

! mark those reports which are too old to be used any more
! (in order to be deleted subseqently in 'obs_cdf_read_org')

  IF (.NOT. lobpfrs)  CALL obs_cdf_mark_old_rep ( aiwthr )
!                     =========================

!-------------------------------------------------------------------------------
! Section 4.0: Read GPS IWV obs from a COST ASCII format file
!-------------------------------------------------------------------------------

! get numbers of reports which have been read prior to the current timestep
  nodrold (1) = ntotml
  nodrold (2) = ntotsg
  nodrold (3) = ntotgp

  nexceed     =  0
  ldo_airmult = .FALSE.
  ldo_pr_out  = .FALSE.

  IF ((lobpfrs) .AND. (lexigps)) THEN

    CALL obs_read_gps
!   =================
  
    IF (my_cart_id == 0)  CLOSE (nugps)
    ldo_pr_out  = .TRUE.
  ENDIF

!-------------------------------------------------------------------------------
! Section 4: Read from NetCDF files, select, assign to grid points, and
!            distribute observational reports to nodes (processors) according
!            to the sub-domains which contain the assigned grid points;
!            then pre-process and store the reports in the ODR arrays;
!            finally check for redundancy of reports; also clean up the ODR
!            arrays from reports which are not used any more
!-------------------------------------------------------------------------------

  n_totin  =  n_cdfin + n_fofin

  IF (n_totin >= 1) THEN

    CALL obs_cdf_read_org ( lneedob(1:n_totin), min_sta(1:n_totin)             &
                          , min_end(1:n_totin), ycdfdir, icdfdirlen            &
                          , nodrold, nexceed(1:3), ldo_airmult, ldo_pr_out )
!   =====================

! bias correction dep on solar zenith angle for Vaisala RS92 radiosonde humidity

!   CALL obs_cdf_raso_rh_bias ( mqcorr92, nodrold(1), nodrcdf(1) )
    CALL obs_cdf_raso_rh_bias ( mqcorr92, nodrold(1), ntotml )
!   =========================

  ENDIF

!-------------------------------------------------------------------------------
! Section 5: Produce vertical profiles (multi-level reports) from single-level
!            aircraft reports, and save total number of stored ODR reports
!-------------------------------------------------------------------------------

  nexceair (1) = 0
  nexceair (2) = 0

  IF (ldo_airmult) THEN
!   PRINT *,'callair0 ', lsvcorl, nodrold(1), nodrold(2), thairh, ntotml, ntotsg

! this call is not required since 'obs_cdf_interface' has already been called
!   CALL obs_air_interface (num_compute, my_cart_id, icomm_cart, imp_reals     &
!                          ,imp_integers,degrad, imaxmll, imaxsgl, acthr, lwonl)

    IF (lsvcorl) THEN

      ALLOCATE (odp (MAX(ntotsg,1)) , STAT=istat )
      DO irep = 1 , ntotsg
        io   = mosghd (irep,nhio)
        jo   = mosghd (irep,nhjo)
        zpob = osgbdy (irep,nbsp)
!       odp (irep)  =  p2dp ( zpob , io , jo )
!                      ====
        odp (irep)  =  f_p2dp ( zpob , ke , zp(io,jo,:) , dp0(io,jo,:) )
!                      ======
      ENDDO
      tscale   = MAX( MIN( wtukara , wtukare )                                 &
                    , (wtukara + wtukare)/ 10.0_wp ) / c2

      ! get vertical correlation scale as by the namelist,
      ! convert it to ln(p) units
      DO ivar = 1 , 4
        vscale (ivar) = SQRT( vcorls(ivar) )
        IF (msprpar == 2) vscale (ivar) = vscale(ivar) / thdlnpt
      ENDDO
!     PRINT *,'callair2 ', MAXVAL( odp ), MAXVAL( -odp )
!     PRINT *,'callair3 ', vscale
!     PRINT *,'callair4 ', tscale

      CALL obs_air_org_mult ( lsvcorl, nodrold(2), nodrold(1), thairh          &
                            , lobpfrs, isubpos , nexceair                      &
                            , odp , vscale , tscale , rhinfl , rhvfac)
!     =====================

      DEALLOCATE (odp , STAT=istat )

    ELSE

      CALL obs_air_org_mult ( lsvcorl, nodrold(2), nodrold(1), thairh          &
                            , lobpfrs, isubpos , nexceair )
!     =====================

    ENDIF
!   PRINT *,'callair7 ', nodrold(1), nodrold(2), ntotml, ntotsg

    nexceed (1)  =  nexceed(1) + nexceair(1)

  ENDIF

!-------------------------------------------------------------------------------
! Section 6: Model-independent quality control check for multi-level data
!-------------------------------------------------------------------------------

  IF (ntotml > nodrold(1)) THEN
    ALLOCATE ( ievnt (ntotml+1) , STAT=istat )
    DO irep = nodrold(1)+1 , ntotml
      ievnt (irep) = i_cma ( momlhd(irep,nhobtp), momlhd(irep,nhcode) )
!                    =====
    ENDDO

    CALL obs_cdf_mult_qualicheck ( nodrold(1), ntotml, maxmlv, ievnt, rdocp )
!   ============================

    DEALLOCATE ( ievnt , STAT=istat )
  ENDIF

! required only before library 2 on obs operators is ready:

! IF (ldo_pr_out)  CALL obs_copy_err_status ( 1, nodrold, ntotml, ntotsg )
!                  ========================

!-------------------------------------------------------------------------------
! Section 7: Print out ODR and diagnostic output (+ de-allocate 'outbuf')
!-------------------------------------------------------------------------------

  imaxl (1)    =  imaxmll
  imaxl (2)    =  imaxsgl
  imaxl (3)    =  imaxgpl
  imaxl (4)    =  imaxtvl
  nodrtot (1)  =  ntotml
  nodrtot (2)  =  ntotsg
  nodrtot (3)  =  ntotgp

  IF (ldo_pr_out) THEN

    CALL obs_cdf_print_reject ( lwonl , num_compute , my_cart_id , icomm_cart  &
                              , imp_integers )
!   =========================

    CALL obs_cdf_print_odr ( nodrtot , nodrold , num_compute , my_cart_id      &
                           , icomm_cart , imp_integers , 0 )
!   ======================

    CALL obs_print_number_rep ( nodrtot , nodrold , imaxl , lwonl, num_compute &
                              , my_cart_id , icomm_cart , imp_integers )
!   =========================

  ELSE
    DEALLOCATE ( outbuf , STAT=istat )
  ENDIF

!-------------------------------------------------------------------------------
! Section 8: Process satellite radiances and produce retrievals by 1DVAR
!-------------------------------------------------------------------------------

  IF ((my_cart_id == 0) .AND. (l1dvar)) THEN
    !   just for testing !!!!!!!!!!!!!!
    IF (lopen_odr)  CLOSE (nuodr)
    IF (lopen_rej)  CLOSE (nurej)
    lopen_odr = .FALSE.
    lopen_rej = .FALSE.
  ENDIF

  IF (l1dvar)  CALL interface_1dvar ( aiwthr , nlastproc , nexceed(4) )
!              ====================

  ! caution: pointers 'r_p' and 'r_ps' have been reset inside 'interface_1dvar'
  !          and can therefore not be used hereafter without being reset again

  IF (my_cart_id == 0) THEN
    IF (lopen_odr)  CLOSE (nuodr)
    IF (lopen_rej)  CLOSE (nurej)
    lopen_odr = .FALSE.
    lopen_rej = .FALSE.
  ENDIF

!-------------------------------------------------------------------------------
! Section 9: Write alerting ('CAUTION') messages if the ODR array size is
!            insufficient to accommodate all observation reports;
!            the messages include a guess by how much certain namelist
!            parameters should be increased to render the ODR size large enough.
!            (Due to the namelist dependency (e.g. 'maxmlo'), the called
!             routines are part of this module rather than the cdfin-library.)
!-------------------------------------------------------------------------------

! write alerting messages to a specific file (yucautn), if required
! ----------------------------------------------------

  nexcedml     =             nexceed(1)
  nexceed (5)  =  nexceold + nexceed(1)
  nexceold     =  nexcedml

  fsize (1)    =  REAL( maxmlo ,wp ) / MAX( REAL( imaxmll ,wp ) ,c1 )
  fsize (2)    =  REAL( maxsgo ,wp ) / MAX( REAL( imaxsgl ,wp ) ,c1 )
  fsize (3)    =  REAL( maxgpo ,wp ) / MAX( REAL( imaxgpl ,wp ) ,c1 )
  fsize (4)    =  REAL( maxtvo ,wp ) / MAX( REAL( imaxtvl ,wp ) ,c1 )

  CALL obs_cdf_print_caution ( nexceed , nexceair , imaxl , fsize )
! ==========================

!-------------------------------------------------------------------------------
! Section 10: Printing diagnostics: Statistics on processed reports and events
!-------------------------------------------------------------------------------

! IF (lastproc) THEN
  IF (ntstep == nlastproc) THEN

! Print statistics on the processed observational reports (and open 'nustat')

    CALL obs_cdf_print_statist ( 1 , num_compute                               &
                                   , my_cart_id , icomm_cart , imp_integers )
!   ==========================

! further summarise warning messages (to unit number 'nustat')

    CALL obs_cdf_print_eventwarn ( mxreve , neventr , fsize )
!   ============================

! Print summary of report and data events (and finally close 'nustat')

    CALL obs_cdf_print_events (  0 , 0 , mxreve , neventr , num_compute        &
                                       , my_cart_id , icomm_cart , imp_integers)
!   =========================
    CALL obs_cdf_print_events (  0 , 1 , mxdeve , neventd , num_compute        &
                                       , my_cart_id , icomm_cart , imp_integers)
!   =========================
    CALL obs_cdf_print_events ( -1 , 2 , mxdeve , neventd , num_compute        &
                                       , my_cart_id , icomm_cart , imp_integers)
!   =========================

    DEALLOCATE ( neventr   , STAT = istat )
    DEALLOCATE ( neventd   , STAT = istat )
  ENDIF

! flush diagnostic file with unit 'nupr'
! --------------------------------------
! (this is the only file used here which is open
!  when entering and leaving this routine)

! IF (lwonl) THEN
!   CLOSE (nupr)
!   OPEN  (nupr,FILE=yuprint,FORM='FORMATTED',STATUS='OLD'                     &
!                           ,POSITION='APPEND',IOSTAT=nstat)
!   IF (nstat /= 0) yerrmsg = 'OPENING OF FILE yuprint FAILED'
!   IF (nstat /= 0) CALL model_abort (my_cart_id, 1008, yerrmsg, yroutine)
! ENDIF

  lobpfrs = .FALSE.

! required only before library 2 on obs operators is ready:

! IF (ldo_pr_out)  CALL obs_copy_err_status ( 1, nodrold, ntotml, ntotsg )
!                  ========================

!-------------------------------------------------------------------------------
! End Subroutine obs_org_cdf_proc
!-------------------------------------------------------------------------------
END SUBROUTINE obs_org_cdf_proc


!===============================================================================
!+ Module procedure in "src_obs_proc_cdf" for reading and distributing reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_prep_interface ( ie, je, ke , zp, madj_hum                  &
                                  , imaxtvl, imaxmll, imaxsgl, imaxgpl )

!-------------------------------------------------------------------------------
! Description:
!   This procedure of module "src_obs_proc_cdf" prepares those input variables
!   for the interface 'obs_cdf_interface' to the observation processing library
!   which are not yet known:
!    - some constant variables depending on namelist input
!    - some time-dependent model fields
!
! Method:
!   Use namelist parameters and model fields from the 'model' environment.
!
! Written by        : DWD, Christoph Schraff  (original version: 16.01.09)
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
    ie , je , ke        ! dimension of model fields 'zp'

  REAL (KIND=wp)           , INTENT (OUT)    ::  &
    zp    (ie,je,ke)    ! model pressure field

  INTEGER (KIND=iintegers) , INTENT (OUT)    ::  &
    imaxtvl          ,& ! size (report dimension) of the satellite retrieval ODR
    madj_hum            ! = 1 : adjust observed humidity (by ratio of saturation
                        !       vapour pressure over water to the one over ice,
                        !       applied if cloud ice is not a state variable

  INTEGER (KIND=iintegers) , INTENT (OUT) , OPTIONAL   ::  &
    imaxmll          ,& ! size (report dimension) of the  multi-level (m-l)  ODR
    imaxsgl          ,& ! size (report dimension) of the single-level (s-l)  ODR
    imaxgpl             ! size (report dimension) of the (ground-based) GPS  ODR

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  LOGICAL                  , SAVE  ::  &
    lprepfrs = .TRUE.   ! variable for 'first time of obs. processing'

  REAL (KIND=wp)           ::  &
    zfsize           ,& ! factor used to set array sizes
    zfproc              ! additional factor dependent on number of processors
 
! CHARACTER (LEN=20)       ::  &
!   yroutine            ! name of this subroutine
! CHARACTER (LEN=30)       ::  &
!   yerr                ! error message
 
! Local arrays:
! ------------
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_prep_interface
!-------------------------------------------------------------------------------

! yroutine = 'obs_cdf_prep_interface'

!-------------------------------------------------------------------------------
! Section 1: Constant values, dependent on namelist parameters
!-------------------------------------------------------------------------------

  IF (lprepfrs) THEN

! determine the report dimensions of the Observation data record (ODR)
! --------------------------------------------------------------------
! (i.e. the maximum numbers (imaxmll, imaxsgl, imaxpgl, imaxtvl) of reports
!  in a subdomain as a function of the number of nodes, the temporal
!  weight functions, and the total number of reports as specified in
!  the namelist input (by: maxmlo, maxsgo, maxpgo, maxtvo))

    zfproc  = MIN( c1 , c1/10.0_wp + 9.0_wp/(10.0_wp *SQRT(SQRT( c1*num_compute ))) )
!   zfproc  = MIN( c1 , c1/ 5.0_wp  + 4.0_wp/( 5.0_wp *    SQRT( c1*num_compute ) ) )

!   to be adjusted:
    zfsize  = c2
    imaxtvl =  MAX( i1 , NINT( maxtvo * zfproc * zfsize ) )

    IF (      (PRESENT( imaxsgl )) .AND. (PRESENT( imaxmll ))                  &
        .AND. (PRESENT( imaxgpl ))) THEN
      zfsize  = MIN( (  MAX( wtukrsa, wtukara, wtuksua, tipolmx, tipmxsu )     &
                      + MAX( wtuksue, tipmxsu ) + c1)                          &
                    /   MAX( wtuksua+wtuksue, c2*tipmxsu, c1 )  ,  4.0_wp )
      imaxsgl = MAX( i1 , NINT( maxsgo * zfproc * zfsize ) )

      zfsize  = MAX( wtukrsa+wtukrse, wtukara+wtukare, c2*tipolmx, c1 )
!     zfsize  = (zfsize + naofbox) / zfsize
      zfsize  = (zfsize + c1)      / zfsize
      imaxmll = MAX( i1 , NINT( maxmlo * zfproc * zfsize ) )

      zfsize  = MAX( c1 , REAL(MIN( nstop, nudgend ), wp) *dtdeh )
      imaxgpl = MAX( i1 , NINT( maxgpo * zfproc * zfsize ) )
      IF (lwonl) THEN
        WRITE( nupr,'(" MAXIMUM NUMBER OF REPORTS IN THE TOTAL DOMAIN:",       &
                     &"   MAXMLO=",I5,"  MAXSGO=",I5, /,                       &
                     &" MAXIMUM NUMBER OF REPORTS IN SUBDOMAINS:      ",       &
                     &"   MAXMLL=",I5,"  MAXSGL=",I5)' )                       &
               maxmlo, maxsgo, imaxmll, imaxsgl
        WRITE( nupr,'(" MAXIMUM NUMBER OF REPORTS IN THE TOTAL DOMAIN:",       &
                     &"   MAXGPO=",I5,"  MAXTVO=",I5, /,                       &
                     &" MAXIMUM NUMBER OF REPORTS IN SUBDOMAINS:      ",       &
                     &"   MAXGPL=",I5,"  MAXTVL=",I5)' )                       &
               maxgpo, maxtvo, imaxgpl, imaxtvl
      ENDIF
    ELSEIF (lwonl) THEN
        WRITE( nupr,'(" MAXIMUM NUMBER OF RETRIEVALS IN THE TOTAL DOMAIN:",    &
                     &"   MAXTVO=",I5,/,                                       &
                     &" MAXIMUM NUMBER OF RETRIEVALS IN SUBDOMAINS:      ",    &
                     &"   MAXTVL=",I5)' )   maxtvo, imaxtvl
    ENDIF

! decide whether humidity observation values need to be adjusted
! to become compatible with the model
! (this is the case if cloud ice is not a model state variable)
! --------------------------------------------------------------------

    madj_hum  =  0
    IF (itype_gscp <= 2)  madj_hum  =  1

  ENDIF

!-------------------------------------------------------------------------------
! Section 2: Variable quantities
!-------------------------------------------------------------------------------

  zp (:,:,:)  =  p0(:,:,:) + pp(:,:,:,nnew)

  lprepfrs = .FALSE.

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_prep_interface
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_prep_interface


!===============================================================================
!+ Module procedure in "src_obs_proc_cdf" for reading and distributing reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_prep_cma_rules

!-------------------------------------------------------------------------------
! Description:
!   This procedure of module "src_obs_proc_cdf" prepares rules for the use
!   observation and code types.
!
! Method:
!   Array 'cma' of type 't_cmatyp' is updated by the rules that depend on
!   latitute and longitude are written, using namelist parameters.
!
! Written by        : DWD, Christoph Schraff  (original version: 02.01.09)
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments: None
! --------------------

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    icma, ic, icn       ! loop indices
 
  LOGICAL                  ::  &
    lgpc  (0:n_gpc)     ! .TRUE. if GPS processing centre shall be used actively
!   lcd097           ,&
!   lcd098           ,&
!   lcd099           ,&
!   lcd100           ,&
!   lcd101           ,&
!   lcd102           ,&
!   lcd103           ,&
!   lcd104           ,&
!   lcd105           ,&
!   lcd106           ,&
!   lcd107           ,&
!   lcd108

  CHARACTER (LEN=12)       ::  &
    ybuf                ! name of this subroutine
! CHARACTER (LEN=20)       ::  &
!   yroutine            ! name of this subroutine
! CHARACTER (LEN=30)       ::  &
!   yerr                ! error message
     
! Local arrays:
! ------------
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_prep_cma_rules
!-------------------------------------------------------------------------------

! yroutine = 'obs_cdf_prep_cma_rules'

  ! complete table 'cma' (for GPS code types, i.e. sub-centres)
  DO ic = 1 , n_gpc
    cma (n_gpc_offs+ic) %obtyp  =  cma (n_gpc_offs) %obtyp
    cma (n_gpc_offs+ic) %cdtyp  =  gpc(ic) %cdtyp
    ybuf  =  cma (n_gpc_offs) %name (1:12)
    cma (n_gpc_offs+ic) %name   =  ybuf(1:LEN_TRIM(ybuf)) //  ' by '           &
                                                          // gpc(ic) %name
  ENDDO

  ! determine which centres in 'gpc' shall be used actively
  lgpc (0)       =  lgps
  lgpc (1:n_gpc) = .FALSE.
  DO icn = 1 , mxgpc
    IF (igpscen(icn) >= 0) THEN
      DO ic = 1 , n_gpc
        IF (gpc(ic) %gpscen == igpscen(icn))  lgpc (ic) = .TRUE.
      ENDDO
    ENDIF
  ENDDO


  ! preparation of use of GPS processing centres: if centre is given in array
  ! 'igpscen', a corresponding logical lcdxxx is set true for active use;
  ! if centre is not in 'igpscen' then all data of this centre are set passive
! lcd108=.false.
! lcd097=.false.
! lcd098=.false.
! lcd099=.false.
! lcd100=.false.
! lcd101=.false.
! lcd102=.false. 
! lcd103=.false.
! lcd104=.false.
! lcd105=.false.
! lcd106=.false.
! lcd107=.false.
! DO ic = 1 , mxgpc
!   IF (igpscen(ic) == 23)                           lcd108=.TRUE.
!   IF (igpscen(ic) == 21)                           lcd097=.TRUE.
!   IF (igpscen(ic) == 30)                           lcd098=.TRUE.
!   IF (igpscen(ic) == 24)                           lcd099=.TRUE.
!   IF (igpscen(ic) == 25)                           lcd100=.TRUE.
!   IF (igpscen(ic) == 29)                           lcd101=.TRUE.
!   IF (igpscen(ic) == 26)                           lcd102=.TRUE.
!   IF (igpscen(ic) ==  0)                           lcd103=.TRUE.
!   IF((igpscen(ic) == 32) .OR. (igpscen(ic) == 37)) lcd104=.TRUE.
!   IF (igpscen(ic) == 35)                           lcd105=.TRUE.
!   IF (igpscen(ic) == 33)                           lcd106=.TRUE.
!   IF (igpscen(ic) == 34)                           lcd107=.TRUE.
! ENDDO

!-------------------------------------------------------------------------------

  DO icma = 1 , n_cma
    cma(icma)%obslat = MAX( cma(icma)%obslat , obslat )
    cma(icma)%obnlat = MIN( cma(icma)%obnlat , obnlat )
    cma(icma)%obwlon = MAX( cma(icma)%obwlon , obwlon )
    cma(icma)%obelon = MIN( cma(icma)%obelon , obelon )
    IF (cma(icma)%obtyp == nsynop) THEN
      IF (     (.NOT. lsynop)                                                  &
          .OR. ((cma(icma)%cdtyp == nsrscd) .AND. (.NOT. lcd011))              &
          .OR. ((cma(icma)%cdtyp == natscd) .AND. (.NOT. lcd014))              &
          .OR. ((cma(icma)%cdtyp == nshscd) .AND. (.NOT. lcd021))              &
          .OR. ((cma(icma)%cdtyp == nabscd) .AND. (.NOT. lcd022))              &
          .OR. ((cma(icma)%cdtyp == nshred) .AND. (.NOT. lcd023))              &
          .OR. ((cma(icma)%cdtyp == natshs) .AND. (.NOT. lcd024))              &
          .OR.  (cma(icma)%cdtyp == nmetar)                      ) THEN
        cma(icma)%exslat = MIN( cma(icma)%exslat , exslat )
        cma(icma)%exnlat = MAX( cma(icma)%exnlat , exnlat )
        cma(icma)%exwlon = MIN( cma(icma)%exwlon , exwlon )
        cma(icma)%exelon = MAX( cma(icma)%exelon , exelon )
      ENDIF
    ELSEIF (cma(icma)%obtyp == nairep) THEN
      IF (     (.NOT. laircf)                                                  &
          .OR. ((cma(icma)%cdtyp == ncodar) .AND. (.NOT. lcd041))              &
          .OR. ((cma(icma)%cdtyp == naircd) .AND. (.NOT. lcd141))              &
          .OR. ((cma(icma)%cdtyp == namdar) .AND. (.NOT. lcd144))              &
          .OR. ((cma(icma)%cdtyp == nmodes) .AND. (.NOT. lcd146))              &
          .OR. ((cma(icma)%cdtyp == nacar ) .AND. (.NOT. lcd244))              &
          .OR. ((cma(icma)%cdtyp == ncolba) .AND. (.NOT. lcd241))) THEN
        cma(icma)%exslat = MIN( cma(icma)%exslat , exslat )
        cma(icma)%exnlat = MAX( cma(icma)%exnlat , exnlat )
        cma(icma)%exwlon = MIN( cma(icma)%exwlon , exwlon )
        cma(icma)%exelon = MAX( cma(icma)%exelon , exelon )
      ENDIF
    ELSEIF (cma(icma)%obtyp == nsatob) THEN
      IF (     (.NOT. lsatob)                                                  &
          .OR. ((cma(icma)%cdtyp == nstbcd) .AND. (.NOT. lcd088))              &
          .OR.  (cma(icma)%cdtyp == nhrvis)                                    &
          .OR.  (cma(icma)%cdtyp == namv  )                                    &
          .OR. ((cma(icma)%cdtyp == nsst  ) .AND. (.NOT. lcd188))) THEN
        cma(icma)%exslat = MIN( cma(icma)%exslat , exslat )
        cma(icma)%exnlat = MAX( cma(icma)%exnlat , exnlat )
        cma(icma)%exwlon = MIN( cma(icma)%exwlon , exwlon )
        cma(icma)%exelon = MAX( cma(icma)%exelon , exelon )
      ENDIF
    ELSEIF (cma(icma)%obtyp == ndribu) THEN
      IF (     (.NOT. ldribu)                                                  &
          .OR. ((cma(icma)%cdtyp == nbathy) .AND. (.NOT. lcd063))              &
          .OR. ((cma(icma)%cdtyp == ntesac) .AND. (.NOT. lcd064))              &
          .OR. ((cma(icma)%cdtyp == ndrbcd) .AND. (.NOT. lcd165))) THEN
        cma(icma)%exslat = MIN( cma(icma)%exslat , exslat )
        cma(icma)%exnlat = MAX( cma(icma)%exnlat , exnlat )
        cma(icma)%exwlon = MIN( cma(icma)%exwlon , exwlon )
        cma(icma)%exelon = MAX( cma(icma)%exelon , exelon )
      ENDIF
    ELSEIF (cma(icma)%obtyp == ntemp ) THEN
      IF (     (.NOT. ltemp )                                                  &
          .OR. ((cma(icma)%cdtyp == nldtcd) .AND. (.NOT. lcd035))              &
          .OR. ((cma(icma)%cdtyp == nshtcd) .AND. (.NOT. lcd036))              &
          .OR. ((cma(icma)%cdtyp == nmotcd) .AND. (.NOT. lcd037))              &
          .OR. ((cma(icma)%cdtyp == ntdrop) .AND. (.NOT. lcd135))              &
          .OR. ((cma(icma)%cdtyp == nrocob) .AND. (.NOT. lcd039))              &
          .OR. ((cma(icma)%cdtyp == nrocsh) .AND. (.NOT. lcd040))) THEN
        cma(icma)%exslat = MIN( cma(icma)%exslat , exslat )
        cma(icma)%exnlat = MAX( cma(icma)%exnlat , exnlat )
        cma(icma)%exwlon = MIN( cma(icma)%exwlon , exwlon )
        cma(icma)%exelon = MAX( cma(icma)%exelon , exelon )
      ENDIF
    ELSEIF (cma(icma)%obtyp == npilot) THEN
      IF (     (.NOT. lpilot)                                                  &
          .OR. ((cma(icma)%cdtyp == nldpcd) .AND. (.NOT. lcd032))              &
          .OR. ((cma(icma)%cdtyp == nshpcd) .AND. (.NOT. lcd033))              &
          .OR. ((cma(icma)%cdtyp == nmopcd) .AND. (.NOT. lcd038))              &
          .OR. ((cma(icma)%cdtyp == nwp_eu) .AND. (.NOT. lcd132))              &
          .OR. ((cma(icma)%cdtyp == nra_eu) .AND. (.NOT. lcd133))              &
          .OR. ((cma(icma)%cdtyp == npr_us) .AND. (.NOT. lcd136))              &
          .OR. ((cma(icma)%cdtyp == nravad) .AND. (.NOT. lcd137))) THEN
        cma(icma)%exslat = MIN( cma(icma)%exslat , exslat )
        cma(icma)%exnlat = MAX( cma(icma)%exnlat , exnlat )
        cma(icma)%exwlon = MIN( cma(icma)%exwlon , exwlon )
        cma(icma)%exelon = MAX( cma(icma)%exelon , exelon )
      ENDIF
    ELSEIF (cma(icma)%obtyp == nsattv) THEN
      IF (     (.NOT. lsatem)                                                  &
          .OR. ((cma(icma)%cdtyp == nsmsg1) .AND. (mcdmsg1 <= 1))              &
          .OR. ((cma(icma)%cdtyp == nsmsg2) .AND. (mcdmsg2 <= 1))              &
          .OR. ((cma(icma)%cdtyp == nnoa15) .AND. (mcdno15 <= 1))              &
          .OR. ((cma(icma)%cdtyp == nnoa16) .AND. (mcdno16 <= 1))              &
          .OR. ((cma(icma)%cdtyp == nnoa17) .AND. (mcdno17 <= 1))              &
          .OR. ((cma(icma)%cdtyp == nnoa18) .AND. (mcdno18 <= 1))) THEN
        cma(icma)%exslat = MIN( cma(icma)%exslat , exslat )
        cma(icma)%exnlat = MAX( cma(icma)%exnlat , exnlat )
        cma(icma)%exwlon = MIN( cma(icma)%exwlon , exwlon )
        cma(icma)%exelon = MAX( cma(icma)%exelon , exelon )
      ENDIF
      IF (         ((cma(icma)%cdtyp == nsmsg1) .AND. (mcdmsg1 >= 1))          &
              .OR. ((cma(icma)%cdtyp == nsmsg2) .AND. (mcdmsg2 >= 1))) THEN
        cma(icma) %attr = 1
      ELSEIF (     ((cma(icma)%cdtyp == nnoa15) .AND. (mcdno15 >= 1))          &
              .OR. ((cma(icma)%cdtyp == nnoa16) .AND. (mcdno16 >= 1))          &
              .OR. ((cma(icma)%cdtyp == nnoa17) .AND. (mcdno17 >= 1))          &
              .OR. ((cma(icma)%cdtyp == nnoa18) .AND. (mcdno18 >= 1))) THEN
        cma(icma) %attr = 2
      ENDIF
    ELSEIF (cma(icma)%obtyp == ngps  ) THEN
      IF ((.NOT. lgps) .OR. (.NOT. lgpc(icma-n_gpc_offs))) THEN
!     IF (     (.NOT. lgps  )                                                  &
!         .OR. ((cma(icma)%cdtyp == ngpgfz ) .AND. (.NOT. lcd108 ))            &
!         .OR. ((cma(icma)%cdtyp == ngpasi ) .AND. (.NOT. lcd097 ))            &
!         .OR. ((cma(icma)%cdtyp == ngpbkg ) .AND. (.NOT. lcd098 ))            &
!         .OR. ((cma(icma)%cdtyp == ngpgop ) .AND. (.NOT. lcd099 ))            &
!         .OR. ((cma(icma)%cdtyp == ngpieec) .AND. (.NOT. lcd100 ))            &
!         .OR. ((cma(icma)%cdtyp == ngpsgn ) .AND. (.NOT. lcd101 ))            &
!         .OR. ((cma(icma)%cdtyp == ngplpt ) .AND. (.NOT. lcd102 ))            &
!         .OR. ((cma(icma)%cdtyp == ngpmet ) .AND. (.NOT. lcd103 ))            &
!         .OR. ((cma(icma)%cdtyp == ngprob ) .AND. (.NOT. lcd104 ))            &
!         .OR. ((cma(icma)%cdtyp == ngpige ) .AND. (.NOT. lcd105 ))            &
!         .OR. ((cma(icma)%cdtyp == ngpknmi) .AND. (.NOT. lcd106 ))            &
!         .OR. ((cma(icma)%cdtyp == ngpnga ) .AND. (.NOT. lcd107 ))) THEN
        cma(icma)%exslat = MIN( cma(icma)%exslat , exslat )
        cma(icma)%exnlat = MAX( cma(icma)%exnlat , exnlat )
        cma(icma)%exwlon = MIN( cma(icma)%exwlon , exwlon )
        cma(icma)%exelon = MAX( cma(icma)%exelon , exelon )
      ENDIF
    ELSEIF (cma(icma)%obtyp == nscatt) THEN
      IF (     (.NOT. lscatt)                                                  &
          .OR. ((cma(icma)%cdtyp == nascat) .AND. (.NOT. lcd123))              &
          .OR. ((cma(icma)%cdtyp == nqscat) .AND. (.NOT. lcd122))) THEN
        cma(icma)%exslat = MIN( cma(icma)%exslat , exslat )
        cma(icma)%exnlat = MAX( cma(icma)%exnlat , exnlat )
        cma(icma)%exwlon = MIN( cma(icma)%exwlon , exwlon )
        cma(icma)%exelon = MAX( cma(icma)%exelon , exelon )
      ENDIF
    ELSE
      cma(icma)%exslat = cs
      cma(icma)%exnlat = cn
      cma(icma)%exwlon = cw
      cma(icma)%exelon = ce
    ENDIF
  ENDDO

!     some code types have to be added:
!       METAR, TEMP MOBIL, PILOT MOBIL, etc.
!     following code types are not used any more:
!       kcdtyp == nabscd .AND. .NOT.lcd022    .OR.                             &
!       kcdtyp == nshred .AND. .NOT.lcd023    .OR.                             &
!       kcdtyp == ncodar .AND. .NOT.lcd141    .OR.                             &
!       kcdtyp == ncolba .AND. .NOT.lcd241    .OR.                             &
!       kcdtyp == nacar  .AND. .NOT.lcd244    .OR.                             &
!       kcdtyp == nsst   .AND. .NOT.lcd188    .OR.                             &
!       kcdtyp == nbathy .AND. .NOT.lcd063    .OR.                             &
!       kcdtyp == nrocob .AND. .NOT.lcd039    .OR.                             &
!       kcdtyp == nrocsh .AND. .NOT.lcd040    .OR.                             &

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_prep_cma_rules
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_prep_cma_rules


!===============================================================================
!+ Module procedure in "src_obs_proc_cdf" for deleting 'old' reports in the ODR
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_mark_old_rep ( hract )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_proc_cdf" prepares the deletion of
!   past reports which are not used any more (for nudging or verification) by
!   marking these report through setting 'mo??hd(,.nhpass) = -1'.
!
! Method:
!   Depends on temporal weight functions defined by nudging namelist parameters.
!
! Written by        : DWD, Christoph Schraff
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  REAL    (KIND=wp   )     , INTENT (IN)     ::  &
    hract            ! hour for which analysis increments are re-calculated
 
! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    nmlob , nsgob ,& ! loop indices
    ngpob         ,& ! loop indices
    istat            ! error status variable

  REAL (KIND=wp)           ::  &
    wtuke         ,& ! temporal radius of influence
    sract         ,& ! time [sec] for which analysis increm. are re-calculated
    c3600r           ! 1 / 3600.0_wp 

  LOGICAL                  ::  &
    ldel             ! delete report

! Local arrays:
! ------------

  LOGICAL                 , ALLOCATABLE :: &
    lmsg      (:)    ! caution message needs to be printed for current report
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_mark_old_rep
!-------------------------------------------------------------------------------

  ALLOCATE ( lmsg  (MAX( ntotsg, ntotml, ntotgp, 1 )) , STAT=istat )
  c3600r = c1 / c3600
  sract  = hract * c3600

!-------------------------------------------------------------------------------
!  Section 1: Multi-level reports
!-------------------------------------------------------------------------------

  DO nmlob = 1 , ntotml
    IF (momlhd(nmlob,nhobtp) == nairep) THEN
      wtuke = wtukare
    ELSEIF (      (momlhd(nmlob,nhobtp) == npilot)                             &
            .AND. (     (momlhd(nmlob,nhcode) == nwp_eu)                       &
                   .OR. (momlhd(nmlob,nhcode) == nra_eu)                       &
                   .OR. (momlhd(nmlob,nhcode) == npr_us)                       &
                   .OR. (momlhd(nmlob,nhcode) == nravad))) THEN
!       this must be the same as in Section 1 of routine 'local_sort_reports'
!       in module 'src_sing_local.f90' where the temporal weights are hard-coded
!       for these observation code types !
      wtuke = MAX( MIN( 0.2_wp, wtukrse ) , 0.5_wp )
    ELSE
      wtuke = MAX( wtukrse , tipolmx )
    ENDIF
    ldel  = (omlhed(nmlob,nhtime) + wtuke  <=  hract)
    lmsg (nmlob) = .FALSE.
    IF ((lverif) .AND. (momlhd(nmlob,nhqcfw) > -99)                            &
                 .AND. (omlhed(nmlob,nhtime) >= hversta-epsy)                  &
                 .AND. (omlhed(nmlob,nhtime) <= hverend+epsy)) THEN
!     (routine is now called if (lgetai))
!     IF (omlhed(nmlob,nhtime) + (dtqc+tconbox)*c3600r  <=  ntstep*dtdeh) THEN
!     IF (ABS( omlhed(nmlob,nhtime)*c3600 +dtqc - ntstep*dt ) <= dt*c05) THEN
      IF (ABS( omlhed(nmlob,nhtime)*c3600 +dtqc -sract ) <= tmaxbox*c05) THEN
        lmsg (nmlob) = .TRUE.
      ELSE
        ldel  = .FALSE.
      ENDIF
    ENDIF
    IF (ldel) momlhd(nmlob,nhpass) = -1
  ENDDO

  DO nmlob = 1 , ntotml
    IF (lmsg(nmlob))                                                           &
      PRINT '(2X,"NOTE: m-l report ",A ," ,",F5.2," [hrs],",F7.0               &
            &," [Pa] maybe not written to YUVERIF")' ,                         &
            yomlhd(nmlob), omlhed(nmlob,nhtime), omlbdy(nmlob,1,nbtp)
  ENDDO

!-------------------------------------------------------------------------------
!  Section 2: Single-level reports
!-------------------------------------------------------------------------------

  DO  nsgob = 1 , ntotsg
    wtuke = MAX( wtuksue , tipmxsu )
    IF (     (mosghd(nsgob,nhobtp) == nairep)                                  &
        .OR. (mosghd(nsgob,nhobtp) == nsatob)) wtuke = wtukare
    ldel = (osghed(nsgob,nhtime) + wtuke  <=  hract)
    lmsg (nsgob) = .FALSE.
    IF ((lverif) .AND. (mosghd(nsgob,nhqcfw) > -99)                            &
                 .AND. (osghed(nsgob,nhtime) >= hversta-epsy)                  &
                 .AND. (osghed(nsgob,nhtime) <= hverend+epsy)) THEN
      IF (ABS( osghed(nsgob,nhtime)*c3600 +dtqc -sract ) <= tmaxbox*c05) THEN
        lmsg (nsgob) = .TRUE.
      ELSE
        ldel = .FALSE.
      ENDIF
    ENDIF
    IF (ldel) mosghd(nsgob,nhpass) = -1
  ENDDO

  DO  nsgob = 1 , ntotsg
    IF (lmsg(nsgob))                                                           &
      PRINT '(2X,"NOTE: s-l report ",A ," ,",F5.2," [hrs],",F7.0               &
            &," [Pa] maybe not written to YUVERIF")' ,                         &
            yosghd(nsgob), osghed(nsgob,nhtime), osgbdy(nsgob,nbsp)
  ENDDO

!-------------------------------------------------------------------------------
!  Section 3: GPS reports
!-------------------------------------------------------------------------------

  DO  ngpob = 1 , ntotgp
    wtuke = MAX( wtuksue , tipmxsu )
    ldel = (ogphed(ngpob,nhtime) + wtuke  <=  hract)
    lmsg (ngpob) = .FALSE.
    IF ((lverif) .AND. (mogphd(ngpob,nhqcfw) > -99)                            &
                 .AND. (ogphed(ngpob,nhtime) >= hversta-epsy)                  &
                 .AND. (ogphed(ngpob,nhtime) <= hverend+epsy)) THEN
      IF (ABS( ogphed(ngpob,nhtime)*c3600 +dtqc -sract ) <= tmaxbox*c05) THEN
        lmsg (ngpob) = .TRUE.
      ELSE
        ldel = .FALSE.
      ENDIF
    ENDIF
    IF (ldel) mogphd(ngpob,nhpass) = -1
  ENDDO

  DO  ngpob = 1 , ntotgp
    IF (lmsg(ngpob))                                                           &
      PRINT '(2X,"NOTE: GPS report ",A ," ,",F5.2," [hrs],",F7.0               &
            &," [Pa] maybe not written to YUVERIF")' ,                         &
            yogphd(ngpob), ogphed(ngpob,nhtime), ogpbdy(ngpob,nbgp)
  ENDDO

  DEALLOCATE ( lmsg  , STAT=istat )

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_mark_old_rep
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_mark_old_rep


!===============================================================================
!+ Module procedure in "src_obs_proc_cdf" for printing alerting messages
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_print_caution ( nexceed , nexceair , imaxl , fsize )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "src_obs_proc_cdf" writes alerting
!   ('CAUTION') messages to a specific file (yucautn) if the ODR array size
!   is insufficient to accommodate all observations;
!   the messages include an estimation by how much certain namelist parameters
!   should be increased to render the ODR size large enough.
!
! Method:
!   If the program is run in parallel mode, the maximum over all nodes of the
!   local excess is determined and used for the estimation.
! 
! Written by        : DWD, Christoph Schraff  (original version: 16.01.09)
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

  INTEGER (KIND=iintegers), INTENT (INOUT)  ::  &
    nexceed  (5)     ,& ! number of reports in excess of ODR array size
    nexceair (2)        ! number of multi-level reports derived from single-
                        !           level reports but in excess of array size

  INTEGER (KIND=iintegers), INTENT (IN)     ::  &
    imaxl    (4)        ! size (report dimension) of the (4 types of) ODR

  REAL    (KIND=wp)       , INTENT (IN)         ::       &
    fsize    (4)        ! ratio of namelist parameter (max??o)
                        !       to ODR size (max??l)

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    ierr  , nstat       ! error status
    
  CHARACTER (LEN=21)       ::  &
    yroutine            ! name of this subroutine
  CHARACTER (LEN=30)       ::  &
    yerr                ! error message

! Local arrays: 
! ------------ 
! 
!------------ End of header ----------------------------------------------------
    
  yroutine = 'obs_cdf_print_caution'

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_print_caution
!-------------------------------------------------------------------------------

! collect counters
! ----------------

  IF (num_compute > 1) THEN

    CALL global_values ( nexceed, 5, 'MAX',imp_integers,icomm_cart, 0,yerr,ierr)
!   ==================

    IF (ierr == 0)                                                             &

      CALL global_values (nexceair, 2,'MAX',imp_integers,icomm_cart,0,yerr,ierr)
!     ==================

    IF (ierr /= 0)  CALL model_abort (my_cart_id, 11015, yerr, yroutine)
!                   ----------------
  ENDIF

! write messages
! --------------

  IF ((my_cart_id == 0) .AND. (     (MAXVAL(nexceed(1:4)) > 0)                 &
                               .OR. (MAXVAL(nexceair(1:2)) > 0))) THEN
    OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'            &
                               , POSITION='APPEND', IOSTAT=nstat)
    IF (nstat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
    IF (nstat /= 0) CALL model_abort (my_cart_id, 1409, yerr, yroutine)
    IF (nexceed(2) > 0) THEN
      WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," LOCAL SINGLE-LEVEL OBS." &
                     &," BEYOND maxsgl ",I5)' ) ntstep, nexceed(2), imaxl(2)
      nexceed(2)  =  NINT( nexceed(2) *fsize(2) ) + 1
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxsgo BY AT LEAST" &
                     &,I6)' ) nexceed(2)
    ENDIF
    ! write alerting messages if array size is insufficient to accommodate all
    ! multi-level reports that can be derived from single-lev (aircraft) reports 
    IF (nexceair(2) > 0) THEN
      WRITE( nucautn,'("CAUTION !!!!! t=",F6.3,":",I5," NEW MULTI-LEV. AIRCR"  &
                     &,"AFT REPS BEYOND ARRAY SIZE ")' ) acthr, nexceair(2)
      nexceair(2)  =  NINT( nexceair(2) *fsize(1) ) + 1
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxmlo BY AT LEAST" &
                     &,I6)' ) nexceair(2)
    ENDIF
    IF (nexceed(1) > 0) THEN
      IF (nexceair(1) > 0)                                                     &
        WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," LOC MULTI-LEV. AIRCR"  &
                       &,"AFT REPORTS BEYOND ARRAY SIZE ")') ntstep, nexceair(1)
      WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," LOCAL MULTI-LEVEL OBS."  &
                     &," BEYOND maxmll ",I5)' ) ntstep, nexceed(1), imaxl(1)
      nexceed(1)  =  NINT( nexceed(1) *fsize(1) ) + 1
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxmlo BY AT LEAST" &
                     &,I6)' ) nexceed(1)
      IF (nexceed(5) > nexceed(1)) THEN
        nexceed(5)  =  NINT( nexceed(5) *fsize(1) ) + 1
        WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxmlo BY ABOUT " &
                       &,I6)' ) nexceed(5)
        WRITE( nucautn,'("     OR THESE AIRCRAFT REPORTS ARE ASSIMILATED AS S" &
                       &,"INGLE-LEVEL REPORTS")' )
      ENDIF
    ENDIF
    IF (nexceed(3) > 0) THEN
      WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," LOCAL GPS (IWV) OBS."    &
                     &," BEYOND maxgpl ",I5)' ) ntstep, nexceed(3), imaxl(3)
      nexceed(3)  =  NINT( nexceed(3) *fsize(3) ) + 1
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxgpo BY AT LEAST" &
                     &,I6)' ) nexceed(3)
    ENDIF
    IF (nexceed(4) > 0) THEN
      WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," SATELLITE RETRIEVALS"   &
                     &," BEYOND maxtvl ",I5)' ) ntstep, nexceed(4), imaxl(4)
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxtvo BY AT LEAST" &
                     &,I6)' ) nexceed(4)
    ENDIF
    CLOSE (nucautn)
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure obs_cdf_print_caution
!-------------------------------------------------------------------------------

END SUBROUTINE obs_cdf_print_caution


!===============================================================================
!+ Module procedure in "src_obs_proc_cdf" for printing alerting messages
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_print_caut_neff ( nexce_rep, nexce_bdy, max_rep, max_body )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "src_obs_proc_cdf" writes alerting
!   ('CAUTION') messages to a specific file (yucautn) if the NetCDF feedobs file
!   size is insufficient to accommodate all reports or observations;
!   the messages include a guess by how much certain namelist parameters should
!   be increased to render the ODR size large enough.
!
! Method:
! 
! Written by        : DWD, Christoph Schraff  (original version: 03.01.11)
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

  INTEGER (KIND=iintegers), INTENT (INOUT)  ::  &
    nexce_rep        ,& ! number of reports exceeding NetCDF feedobs file size
    nexce_bdy        ,& ! number of obs     exceeding NetCDF feedobs file size
    max_rep          ,& ! max. number of reports in NetCDF feedobs file (FOF)
    max_body            ! max. number of obs in NetCDF feedobs file (FOF)

! Local parameters: None
! ----------------

! Local varibles:
! --------------

  INTEGER (KIND=iintegers) ::  &
    nexceed  (2)     ,& ! number of reports / obs in excess of feedobs file size
    nall             ,& ! total number of reports / obs to be written to file
    ierr  , nstat       ! error status
    
  CHARACTER (LEN=23)       ::  &
    yroutine            ! name of this subroutine
  CHARACTER (LEN=30)       ::  &
    yerr                ! error message

! Local arrays: 
! ------------ 
! 
!------------ End of header ----------------------------------------------------
    
  yroutine = 'obs_cdf_print_caut_neff'

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_print_caut_neff
!-------------------------------------------------------------------------------

! collect counters
! ----------------

  nexceed (1) = nexce_rep
  nexceed (2) = nexce_bdy
  IF (num_compute > 1) THEN

    CALL global_values ( nexceed, 2, 'MAX',imp_integers,icomm_cart, 0,yerr,ierr)
!   ==================
    IF (ierr /= 0)  CALL model_abort (my_cart_id, 11016, yerr, yroutine)
!                   ----------------
  ENDIF

! write messages
! --------------

  IF ((my_cart_id == 0) .AND. (MAXVAL( nexceed ) > 0)) THEN
    OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'            &
                               , POSITION='APPEND', IOSTAT=nstat)
    IF (nstat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
    IF (nstat /= 0) CALL model_abort (my_cart_id, 1410, yerr, yroutine)
    IF (nexceed(1) > 0) THEN
      nall  = max_rep  + nexceed(1)
      WRITE( nucautn,'("CAUTION !!!!! total number of reports ",I6             &
                     &," > FOF size max_rep =",I6)' )  nall, max_rep
      IF (mxfrep <= 0) THEN
        nall  = nexceed(1) / MAX( 1 , NINT( hverend - hversta ) )
        WRITE( nucautn,'("   ==>  INCREASE SUM OF NAMELIST VARIABLES maxmlo "  &
                       &,"+ maxsgo + 2*maxgpo + 2*maxtvo")' )
        WRITE( nucautn,'(8X,"BY AT LEAST",I7," FOR NetCDF FEEDOBS FILE fof_*")'&
             ) nall
      ELSE
        WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE mxfrep BY "       &
                       &,"AT LEAST",I7," FOR FEEDOBS FILE fof_*")' ) nexceed(1)
      ENDIF
    ENDIF
    IF (nexceed(2) > 0) THEN
      nall  = max_body + nexceed(2)
      WRITE( nucautn,'("CAUTION !!!!! total number of obs ",I10                &
                     &," > FOF size max_body =",I10)' )  nall, max_body
      IF (mxfobs <= 0) THEN
        nall  = nexceed(2) / MAX( 1 , NINT( hverend - hversta ) )
        WRITE( nucautn,'("   ==>  INCREASE SUM OF NAMELIST VARIABLES "         &
                       &,"1.6*maxmlo*maxmlv + 8*maxsgo")' )
        WRITE( nucautn,'("        + 6*maxgpo + 4*ke*maxtvo BY AT LEAST",I7     &
                       &," FOR NetCDF FEEDOBS FILE")' )  nall
      ELSE
        WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE mxfobs BY "       &
                       &,"AT LEAST",I7," FOR FEEDOBS FILE fof_*")' ) nexceed(2)
      ENDIF
    ENDIF
    CLOSE (nucautn)
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure obs_cdf_print_caut_neff
!-------------------------------------------------------------------------------

END SUBROUTINE obs_cdf_print_caut_neff


!===============================================================================
!+ Module procedure in "src_obs_proc_cdf" for printing warning messages
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_print_eventwarn ( mxeve , nevent , fsize )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "src_obs_proc_cdf" writes summarising
!   warning messages to the file with unit number 'nustat', if the size of the
!   Observation Data Record (ODR) as specified by namelist parameters has been
!   too small to accommodate all observation reports.
!
! Method:
!   If the program is run in parallel mode, the diagnostic arrays summed up
!   over all nodes.
! 
! Written by        : DWD, Christoph Schraff  (original version: 05.01.09)
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

  INTEGER (KIND=iintegers), INTENT (IN)         ::       &
    mxeve               ,& ! length of first dimension of record 'nevent'
    nevent (mxeve,n_cma)   ! record containing the event counters

  REAL    (KIND=wp)       , INTENT (IN)         ::       &
    fsize  (4)             ! ratio of namelist parameter (max??o)
                           !       to ODR size (max??l)

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    nwarn  (4)       ,& ! local max of events due to ODR size limit
    ima              ,& ! loop index over array 'cma'
    izerror             ! error status
    
  CHARACTER (LEN=80)       ::  yzerrmsg

! Local arrays: 
! ------------ 
! 
!------------ End of header ----------------------------------------------------
    
  izerror  = 0
  yzerrmsg = '   '
    
!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_print_eventwarn
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1: Writing of messages which warn of insufficient array size
!-------------------------------------------------------------------------------

! sum up counters for warning messages
! ------------------------------------

  nwarn (:) = 0
  DO ima = 1 , n_cma
    IF (     (cma(ima)%obtyp == nsynop) .OR. (cma(ima)%obtyp == ndribu)        &
        .OR. (cma(ima)%obtyp == nsatob) .OR. (cma(ima)%obtyp == nscatt)) THEN
      nwarn (2)  =  nwarn(2)  +  nevent(nesodr,ima)
    ELSEIF   (cma(ima)%obtyp == nairep)                                  THEN
      nwarn (2)  =  nwarn(2)  +  nevent(nesodr,ima)
      nwarn (1)  =  nwarn(1)  +  nevent(nenoml,ima)
    ELSEIF ( (cma(ima)%obtyp == ntemp ) .OR. (cma(ima)%obtyp == npilot)) THEN
      nwarn (1)  =  nwarn(1)  +  nevent(nesodr,ima)
    ELSEIF   (cma(ima)%obtyp == ngps  )                                  THEN
      nwarn (3)  =  nwarn(3)  +  nevent(nesodr,ima)
    ELSEIF   (cma(ima)%obtyp == nsatem)                                  THEN
      nwarn (4)  =  nwarn(4)  +  nevent(nesodr,ima)
    ENDIF
  ENDDO
  IF (num_compute > 1) THEN

    CALL global_values ( nwarn, 4, 'MAX', imp_integers                         &
                       , icomm_cart, 0, yzerrmsg, izerror )
!   ===================

  ENDIF

! print warning messages
! ----------------------
  
  IF (my_cart_id == 0) THEN
    nwarn (1) = NINT( REAL( nwarn(1), wp ) * fsize(1) + 0.49_wp )
    nwarn (2) = NINT( REAL( nwarn(2), wp ) * fsize(2) + 0.49_wp )
    nwarn (3) = NINT( REAL( nwarn(3), wp ) * fsize(3) + 0.49_wp )
    nwarn (4) = NINT( REAL( nwarn(4), wp ) * fsize(4) + 0.49_wp )
!   IF (MAXVAL( nwarn ) > 0) THEN 
    IF (MAX( nwarn(1),nwarn(2),nwarn(3),nwarn(4) ) > 0) THEN
      WRITE( nustat,'(''1'')' ) 
      WRITE( nustat,'(''0'')' )
      WRITE( nustat,'(''+     !!! CAUTION !!!!! CAUTION !!!!! ''               &
                    &,''CAUTION !!!!! CAUTION !!!!! CAUTION !!!!!'')' )
      WRITE( nustat,'(''+     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~''               &
                    &,''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'')' )
    ENDIF
    IF (nwarn(1) > 0)                                                          &
        WRITE( nustat,'(''0     WARNING: array size for multi-level ''         &
                    &,''observations is too small'', /                         &
                    &,''      =======  to accommodate all observations'', /    &
                    &,''       --> Increase "MAXMLO" (namelist) by'',I5        &
                    &,'': usually ok for local obs. array'', /                 &
                    &,''           =================                   ''      &
                    &,'' (possibly still insufficient for'', /                 &
                    &,''                                               ''      &
                    &,'' the global obs increment array!)'')' )        nwarn(1)
    IF (nwarn(2) > 0)                                                          &
      WRITE( nustat,'(''0     WARNING: array size for single-level ''          &
                    &,''observations is too small'', /                         &
                    &,''      =======  to accommodate all observations'', /    &
                    &,''       --> Increase "MAXSGO" (namelist) by'',I5        &
                    &,'': usually ok for local obs. array'', /                 &
                    &,''           =================                   ''      &
                    &,'' (possibly still insufficient for'', /                 &
                    &,''                                               ''      &
                    &,'' the global obs increment array!)'')' )        nwarn(2)
    IF (nwarn(3) > 0)                                                          &
        WRITE( nustat,'(''0     WARNING: array size for GPS ''                 &
                    &,''observations is too small'', /                         &
                    &,''      =======  to accommodate all observations'', /    &
                    &,''       --> Increase "MAXGPO" (namelist) by'',I5        &
                    &,'': usually ok for local obs. array'', /                 &
                    &,''           =================                   ''      &
                    &,'' (possibly still insufficient for'', /                 &
                    &,''                                               ''      &
                    &,'' the global obs increment array!)'')' )        nwarn(3)
    IF (nwarn(4) > 0)                                                          &
        WRITE( nustat,'(''0     WARNING: array size for SAT ''                 &
                    &,''retrievals   is too small'', /                         &
                    &,''      =======  to accommodate all observations'', /    &
                    &,''       --> Increase "MAXTVO" (namelist) by'',I5        &
                    &,'': usually ok for local obs. array'', /                 &
                    &,''           =================                   ''      &
                    &,'' (possibly still insufficient for'', /                 &
                    &,''                                               ''      &
                    &,'' the global obs increment array!)'')' )        nwarn(4)
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure obs_cdf_print_eventwarn
!-------------------------------------------------------------------------------

END SUBROUTINE obs_cdf_print_eventwarn


!===============================================================================
!+ Module procedure in "src_obs_proc_cdf" for deleting 'old' reports in the ODR
!-------------------------------------------------------------------------------

SUBROUTINE obs_1dvar_mark_old_rep ( hract )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_proc_cdf" prepares the deletion of
!   past retrievals which are not used any more (for nudging or verification) by
!   marking these report through setting 'motvhd(,.nhpass) = -1'.
!
! Method:
!   Depends on temporal weight functions defined in 'data_1dvar.f90'.
!
! Written by        : DWD, Christoph Schraff
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  REAL    (KIND=wp   )     , INTENT (IN)     ::  &
    hract            ! hour for which analysis increments are re-calculated
 
! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    nrtvo         ,& ! loop index
    istat            ! error status variable

  REAL (KIND=wp)           ::  &
    wtuke         ,& ! temporal radius of influence
    sract            ! time [sec] for which analysis increm. are re-calculated

  LOGICAL                  ::  &
    ldel             ! delete report

! Local arrays:
! ------------

  LOGICAL                 , ALLOCATABLE :: &
    lmsg      (:)    ! caution message needs to be printed for current report
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_1dvar_mark_old_rep
!-------------------------------------------------------------------------------

  ALLOCATE ( lmsg  (MAX( ntottv, 1 )) , STAT=istat )
  sract  = hract * c3600

  DO  nrtvo = 1 , ntottv
    wtuke = MAX( otvhed(nrtvo,nh1wte) , otvhed(nrtvo,nh1tip) )
    ldel    = (otvhed(nrtvo,nhtime) + wtuke  <=  hract)
    lmsg (nrtvo) = .FALSE.
    IF ((lverif) .AND. (motvhd(nrtvo,nhqcfw) > -99)                            &
                 .AND. (otvhed(nrtvo,nhtime) >= hversta-epsy)                  &
                 .AND. (otvhed(nrtvo,nhtime) <= hverend+epsy)) THEN
!     IF (otvhed(nrtvo,nhtime) + (dtqc+tconbox)/c3600  <= acthr) THEN
      IF (ABS( otvhed(nrtvo,nhtime)*c3600 +dtqc -sract) <= tmaxbox*c05) THEN
        lmsg (nrtvo) = .TRUE.
      ELSE
        ldel = .FALSE.
      ENDIF
    ENDIF
    IF (ldel) motvhd(nrtvo,nhpass) = -1
  ENDDO

  DO  nrtvo = 1 , ntottv
    IF (lmsg(nrtvo))                                                           &
      PRINT '(2X,"NOTE: sat retrieval ",A ," ,",F5.2," [hrs],",F7.0            &
            &," [Pa] maybe not written to YUVERIF")' ,                         &
            yotvhd(nrtvo), otvhed(nrtvo,nhtime), otvbdy(nrtvo,1,nbvp)
  ENDDO

  DEALLOCATE ( lmsg  , STAT=istat )

!-------------------------------------------------------------------------------
! End Subroutine obs_1dvar_mark_old_rep
!-------------------------------------------------------------------------------
END SUBROUTINE obs_1dvar_mark_old_rep


!===============================================================================
!+ Module procedure in "src_obs_proc_cdf" for interfacing the 1dvar / retrievals
!-------------------------------------------------------------------------------

SUBROUTINE interface_1dvar ( tactread , nlastproc , nexceed )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_proc_cdf" just contains the
!   interface to subroutine 'organize_1dvar'. This allows for having only one
!   single call of subroutine 'organize_1dvar' in the total COSMO code and
!   thus having only one single point with linking option 'ifdef' for
!   ex-/including the 1DVAR modules in a COSMO binary.
!
! Written by        :  Christoph Schraff, DWD  (original version: 29.01.10)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!===============================================================================
!
! Declarations:
!
!===============================================================================

IMPLICIT NONE

!===============================================================================

! Subroutine arguments:
! --------------------

  REAL    (KIND=wp)        , INTENT (IN)     ::  &
    tactread         ! reference time [h] to define which sat obs need to be
                     ! read now (tactread = aiwthr,
                     !           where aiwthr is the time for which the next
                     !           analysis increments are valid)

  INTEGER (KIND=iintegers) , INTENT (IN)     ::  &
    nlastproc        ! timestep at which this routine is called the last time

  INTEGER (KIND=iintegers) , INTENT (OUT)    ::  &
    nexceed          ! number of reports in excess of ODR array size

! Local parameters: None
! ----------------

! Local variables: None
! ----------------

  INTEGER (KIND=iintegers) , SAVE  ::  &
    imaxtvl       ,& ! size (report dimension) of the satellite retrieval ODR
    madj_hum         ! mode switch on adjusting observed humidity to model
 
  LOGICAL                  , SAVE  ::  &
    lobpfrs = .TRUE. ! variable for 'first time of obs. processing'

  INTEGER (KIND=iintegers) ::  &
    imaxmll       ,& ! size (report dimension) of the  multi-level (m-l)  ODR
    imaxsgl       ,& ! size (report dimension) of the single-level (s-l)  ODR
    imaxgpl          ! size (report dimension) of the (ground-based) GPS  ODR

  REAL (KIND=wp)           , TARGET  ::  &
    zp  (ie,je,ke)   ! model pressure field

  REAL (KIND=wp)           ::  &
    hstop         ,& ! end of the forecast in hours
    hmaxbox       ,& ! maximum interval for current routine being called [hrs]
                     !   = max. interval for computing analysis increments
    hlastproc        ! hour at which this routine is called the last time
 
  CHARACTER (LEN=icdfdirlen)       ::  &
    ysatdir          ! directory for satellite datasets and parameters
  CHARACTER (LEN=20)       ::  &
    yroutine         ! name of this subroutine
  CHARACTER (LEN=40)       ::  &
    yerrmsg          ! error message

  INTEGER (KIND=iintegers)  :: kzdims(24)
  INTEGER (KIND=iintegers)  :: izerror, istat
  CHARACTER (LEN=255)       :: yzerrmsg

! Tracer pointers:
! -----------------
  REAL (KIND=wp)     , POINTER  :: &
    qv(:,:,:) => NULL()            ! QV at nnew
!
!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine interface_1dvar
!-------------------------------------------------------------------------------

  yroutine = 'interface_1dvar'

! Retrieve the required microphysics tracers (at nnew)
  CALL trcr_get( izerror, idt_qv, ptr_tlev = nnew, ptr = qv )
  IF (izerror /= 0_iintegers) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort( my_cart_id, izerror, yzerrmsg, 'interface_1dvar' )
  ENDIF

!-------------------------------------------------------------------------------
! Section 1: Preparations for the 1DVAR pre-processing
!-------------------------------------------------------------------------------

! first exchange 1 row of boundary data (of local domains)
! (required for bi-linear interpolation from neighbouring 4 grid pts,
!  because the values within the boundary zone of the sub-domains
!  - may not be updated  ( --> no reproducibility )  or
!  - may even be missing ( --> causes minimisation to fail )
! -------------------------------------------------------------------

  IF (num_compute > 1) THEN
    IF (ltime) THEN
      CALL get_timings (i_obs_processing     , ntstep, dt, izerror)
      CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
      CALL get_timings (i_barrier_waiting_nud, ntstep, dt, izerror)
    ENDIF
    kzdims( :  ) = 0_iintegers
    kzdims(1:24) = (/ke,ke,ke,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)

    CALL exchg_boundaries ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart,    &
                           num_compute, ie, je, kzdims, jstartpar, jendpar, 1, &
                           nboundlines, my_cart_neigh, lperi_x, lperi_y, l2dim,&
                           20000+nexch_tag,.FALSE.,ncomm_type,izerror,yzerrmsg,&
                           pp(:,:,:,nnew), t(:,:,:,nnew), qv(:,:,:),           &
                           t_g(:,:,nnew) , ps(:,:,nnew) ,                      &
                           t_2m(:,:), qv_2m(:,:), u_10m(:,:), v_10m(:,:) )
!   =====================

    IF (izerror /= 0) CALL model_abort (my_cart_id, 12345, yzerrmsg            &
                                       ,'CREATE_FIRST_GUESS')
    IF (ltime) CALL get_timings (i_communications_nud , ntstep, dt, izerror)
  ENDIF

! mark those retrievals which are too old to be used any more
!   (so that they can be deleted subseqently in 'obs_del_store_retriev'
!    called within 'organize_1dvar', see module 'src_obs_1dvar_org.f90')
! ----------------------------------------------------------------------

  IF (.NOT. lobpfrs)  CALL obs_1dvar_mark_old_rep ( tactread )
!                     ===========================

!-------------------------------------------------------------------------------
! Section 2: Interface for the 1DVAR pre-processing:
!            make the 'model environment' available to the 1DVAR process modules
!             (this basically means: copying values from the 'model environment'
!              into variables from data modules that are part of the observation
!              pre-processing + 1DVAR library
!              ('library' denotes here simply a set of modules))
!            Remarks:
!            - The values of all the other variables which are set outside the
!              obs pre-proc library but which must be available within a certain
!              1DVAR-related subroutine are transfered directly from the 'model'
!              environment through the argument list of the respective routine.
!            - All the other variables which must be commonly available within
!              the obs pre-proc library but can be set within the library are
!              stored in one of the data modules that are part of the library.
!-------------------------------------------------------------------------------

  ! 'r_p, r_ps' must be re-determined after 'exchg_boundaries'

  CALL obs_cdf_prep_interface ( ie, je, ke, zp, madj_hum, imaxtvl )
! ===========================

  r_p       =>  zp               ! (ie,je,ke)
  r_ps      =>  ps(:,:,nnew)     ! (ie,je)

    ! 2.: environment used only for 1DVAR
    ! 2a: interface for scalars: none

    ! 2b: interface for arrays: set 'library 1' pointers (residing in
    !     data_obs_lib_cosmo.f90) to target arrays (from 'model environment')
    !     (target arrays must have 'target' attribute !)
    !     (this avoids the need to make copies of these arrays,
    !      i.e. no memory (de-)allocation required)
                                   ! target arrays (must) have dimensions:
  r_t       =>  t (:,:,:,nnew)     ! (ie,je,ke)
  r_qv      =>  qv(:,:,:)          ! (ie,je,ke)
  r_t_g     =>  t_g (:,:,nnew)     ! (ie,je)
  r_t_2m    =>  t_2m               ! (ie,je)
  r_qv_2m   =>  qv_2m              ! (ie,je)
  r_u_10m   =>  u_10m              ! (ie,je)
  r_v_10m   =>  v_10m              ! (ie,je)

!-------------------------------------------------------------------------------
! Section 3: Perform the reading of the satellite radiances and
!            the 1DVAR pre-processing by calling the master routine for 1DVAR
!-------------------------------------------------------------------------------

  hmaxbox   = tmaxbox / c3600
  ! careful with that axe, Eugene !
  !   hstop must be computed here from nstop, since
  !   hstop from data_modelconfig.f90 is correct only if (my_cart_id == 0) !!
  hstop     = nstop     * dtdeh
  hlastproc = nlastproc * dtdeh
  nexceed   = 0
  ysatdir   = ycdfdir

! CALL organize_1dvar ( tactread, hstop, hmaxbox, hlastproc                    &
!                     , icdfdirlen, ysatdir , l1dvar, nexceed )
! ===================

!-------------------------------------------------------------------------------
! Section 4: Clean up the 1DVAR pre-processing and its interface
!-------------------------------------------------------------------------------

  !   flush file YUPRINT
  IF (lwonl) THEN
    CLOSE (nupr)
    OPEN  (nupr,FILE=yuprint,FORM='FORMATTED',STATUS='OLD'                     &
                            ,POSITION='APPEND',IOSTAT=istat)
    IF (istat /= 0) yerrmsg = 'OPENING of file YUPRINT failed'
    IF (istat /= 0) CALL model_abort (my_cart_id, 1009, yerrmsg, yroutine)
  ENDIF

  lobpfrs = .FALSE.

!-------------------------------------------------------------------------------
! End of module procedure interface_1dvar
!-------------------------------------------------------------------------------

END SUBROUTINE interface_1dvar


!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_proc_cdf" for reading GPS obs from an ASCII file
!-------------------------------------------------------------------------------

SUBROUTINE  obs_read_gps

!-------------------------------------------------------------------------------
! Description:
!   This module procedure performs GPS-IWV report processing.
!
! Method:
!   GPS data extraction from ascii files in COST-716 format.
!   Reference: "Format specification for COST-716 processed data."
!   COST-716 Document. Version 1.0-3, 28 September 2000.
!   More data files can be joined in the input file.
!   Values less than 0.0 are considered missing values.
!   The grid point assignement is perfomed as for synop observations.
!   A redundancy check is also included.
!
!   Modifications:
!   - V3.3  2003/04/22 Maria Tomassini + C. Schraff:   Initial version.
!   - V3.18 2006/03/03 C. Schraff:        GPS processing adapted for IWV
!                                         spatial consistency check.
!   - V4_22 2012/01/31 Karolin Eichler:   Reprocessed GPS data can be read from
!                                         ASCII files in COST-716 V2.0 (or V1.0)
!                                         format.
!
! Written by        : DWD, Maria Tomassini + Christoph Schraff (2001)
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!===============================================================================
!
! Declarations:
!
!===============================================================================

IMPLICIT NONE

!===============================================================================

! Subroutine arguments: None
! --------------------

! Local parameters:
! -----------------

  INTEGER (KIND=iintegers) ::  &
    nlengps = 29        ! report header length in send buffer

  REAL    (KIND = wp)        , PARAMETER  :: &
    c1000   =  1000.0_wp     ,& !
    c1000r  =   0.001_wp        !

! Local scalars:
! -------------

  LOGICAL                  ::  &
    lprcs               ,& ! obs to be processed (if lverif or liwvscc)
    lsfcob              ,& ! false, if no surface data are to be used
                           ! from a temp or pilot
    lseaobs             ,& ! true, if observation taken on sea
    lcoloc                 ! .TRUE. for spatial or temporal colocation

  INTEGER (KIND=iintegers) ::  &
    ianyy     ,imoyy    ,& ! year     / of the date and time
    ianmm     ,imomm    ,& ! month   /  of the current
    iandd     ,imodd    ,& ! day    <   observation
    ianhh     ,imohh    ,& ! hour    \  resp. of the initial
    ianmin    ,imomin   ,& ! minute   \ model date and time
               imosec   ,& ! second    \ model date and time
    imindif             ,& ! time difference [min] between these two times
    ierrf                  ! error status for time difference routine

  INTEGER (KIND=iintegers) ::  &
    nmxbln,              & ! max. length of buffer containing all reports
    ngpob,               & ! GPS report loop index
    icma,                & ! indices for statistics ('cma', events)
    istat,     nstat,    & ! error status variables
    izmplcode, ierr,     & ! for MPI error code
    izlocstat,           & ! for space allocation error
    npe,                 & ! rank of the receiving PE
    nobbuf,              & ! number of reports in report buffers
    ibuflen,  iallen,    & ! total length of report buffers
    inodcnt,             & ! length of local report buffer
    inodlen(num_compute),& ! total length of report buffer at each node
    inodpos(num_compute),& ! displacement of local report buffers in 'iallbuf'
    ipos,                & ! start address -1 of an observation in report buffer
    irplen,              & ! report length
    ilen,                & ! length of report
    iobdat,              & ! observation date as yymmdd or yyyymodd
    ntime,               & ! actual time of observation
    nelevt,              & ! station altitude
    nstachr,             & ! station characteristics
    nkflag,              & ! processing flag
    mpassiv,             & ! passive flag assigned to report
    mrepflg,             & ! report flag word as part of station charcteristics
    nbufbit,             & ! buffer of bit pattern
    ncentre,             & ! data centre       (WMO Common Code Table C11 ???)
    nexceed,             & ! number of reports in excess of array size
    iob,     job,        & ! grid point asigned to observation (sub-domain)
    iob_tot, job_tot,    & ! grid points assigned to obs (total model area)
    iminmin, iminmax    ,& ! min./ max. time difference [min] obs - model
    istcor,              & ! == 1, if active report is a station correction
    ngplo ,              & ! ODR index of GPS report
    nodract, nodrrej,    & ! index of active / rejected ODR report
    nobtyp , ncdtyp ,    & ! observation/code type
    nobtypr, ncdtypr,    & ! observation/code type of rejected report
    irnode,              & ! process an observation belongs to
    inode,   inobs,      & ! loop indices
    icl,     ii            ! loop indices

  REAL (KIND=wp)           ::  &
    zobtim,              & ! observation time in [hrs], rel. to model init time
    timlim,              & ! limit for time difference in forecast hour units
    rscale,              & ! scaling factor
    rdegkm2,             & ! factor for conversion of degree to km (**2)
    cosolat,             & ! cos (latitude)
    rhzob,               & ! horizontal distance in km
    zsurf,               & ! height of model grid pt. to which obs. is assigned
    fisd,                & ! (scaled) height difference between model orography
                           !    and station altitude
    zio_tot, zjo_tot,    & ! observation position in g.p. units (total area)
    roblon,  roblat,     & ! latitude and longitude in the geographical system
    rlon,    rlat,       & ! latitude and longitude in the rotated system
    tzderr,              & !      / variables
    tzd,                 & !     /  for
    zwd,                 & !    /   temporary storage
    iwv,                 & !   <    of GPS measurements,
    surfpres,            & !    \   as
    surftemp,            & !     \  in
    surfrh                 !      \ GPSbufr
 
 
  CHARACTER (LEN=ilstidp)  ::  &
    ystid   , ystidr       ! station identity
  CHARACTER (LEN=12)       :: &
    yobdat                 ! date and time of current observation
  CHARACTER (LEN=20)       :: &
    yroutine               ! name of this subroutine
  CHARACTER (LEN=30)       :: &
    yerrmsg                ! error message
  CHARACTER (LEN=40)       ::  &
    yerr                   ! error message

! GPS variables:

! Parameter
  CHARACTER (LEN=4),        PARAMETER :: &
      COSTfmtver = "V1.0" ! COST format version supported here
  CHARACTER (LEN=4),        PARAMETER :: &
      COSTfmtve2 = "V2.0" ! COST format version supported here
  CHARACTER (LEN=10),       PARAMETER :: &
      SOVFmark   = "COST-716 V"    ! Start of vfile marker
  INTEGER (KIND=iintegers), PARAMETER :: &
      EOVFstat    = -999            ! End of vfile status code
  CHARACTER (LEN=100),      PARAMETER :: &
      EOVFmark    = REPEAT("-",100) ! End   of vfile marker
  INTEGER (KIND=iintegers), PARAMETER :: &
      MaxSamples  = 288             ! Max. no. of data samples per vfile
  INTEGER (KIND=iintegers), PARAMETER :: & 
      MaxSlants   = 24              ! Max. no. of slant obs per data sample

! COST File Header info

  TYPE COSTheader
    CHARACTER (LEN=4)  :: FmtVer        ! File format ident code (Vm.n)
    CHARACTER (LEN=4)  :: StnID         ! Station ident code
    CHARACTER (LEN=20) :: RcvrType      ! Receiver type
    CHARACTER (LEN=20) :: AntType       ! Antenna/Radome types
    REAL (KIND=wp)     :: StnLat        ! Station latitude  wrt WGS-84 (deg)
    REAL (KIND=wp)     :: StnLon        ! Station longitude wrt WGS=84 (deg)
    REAL (KIND=wp)     :: StnHt         ! Station height    wrt WGS-84 (m)
    REAL (KIND=wp)     :: StnHtAMSL     ! Station height    wrt EGM96  (m)
    CHARACTER (LEN=20) :: DateTimeFirst ! Date/time of 1st sample(dd-MMM-yyyy)
    CHARACTER (LEN=20) :: DateTimeProc  ! Date/time of process.  (dd-MMM-yyyy)
    CHARACTER (LEN=20) :: ProcCentre    ! Processing centre name
    CHARACTER (LEN=20) :: ProcMethod    ! Processing software type
    CHARACTER (LEN=20) :: OrbitType     ! GNSS orbit source and type
    CHARACTER (LEN=20) :: MetSource     ! Source of surface meteorological data
    INTEGER (KIND=iintegers) :: TimeInc    ! Sample time incr./resolut. (min.)
    INTEGER (KIND=iintegers) :: UpdateInt  ! Batch update interval (minutes)
    INTEGER (KIND=iintegers) :: BatchLen   ! Tot. length of batch data (minutes)
    INTEGER (KIND=iintegers) :: PCDH       ! Header Product Confidence Data
    INTEGER (KIND=iintegers) :: NumSamples ! Number of data samples in file
  END TYPE
  TYPE (COSTheader) Header 

! Slant delay & associated parameters

  TYPE COSTslant
    CHARACTER (LEN=4)  :: SatID         ! GNSS satellite idents
    REAL (KIND=wp)     :: TSD           ! Total Slant delays (mm)
    REAL (KIND=wp)     :: TSDErr        ! Estimates errors in TSD (mm)
    REAL (KIND=wp)     :: Azim          ! Slant azimuth angles wrt North(deg)
    REAL (KIND=wp)     :: Elev          ! Slant elevations wrt horizon(deg)
  END TYPE

! File Data parameters

  TYPE COSTsample
    INTEGER (KIND=iintegers), DIMENSION(3) :: & 
                          TimeStamp     ! Data sample time stamp(hr,min,sec)
    INTEGER (KIND=iintegers)               :: &
                          PCDD          ! Data Product Confidence data
    REAL (KIND=wp)     :: TZD           ! Total Zenith Delay (mm)
    REAL (KIND=wp)     :: TZDErr        ! Estimated error in TZD (mm)
    REAL (KIND=wp)     :: ZWD           ! Zenith Wet Delay (mm)
    REAL (KIND=wp)     :: IWV           ! Integrated Water Vapour (mm)
    REAL (KIND=wp)     :: SurfPres      ! Surface pressure (hPa)
    REAL (KIND=wp)     :: SurfTemp      ! Surface temperature (K)
    REAL (KIND=wp)     :: SurfRH        ! Surface relative humidity (%)
    REAL (KIND=wp)     :: GradientN     ! N/S gradient (mm)
    REAL (KIND=wp)     :: GradientE     ! E/W gradient (mm)
    REAL (KIND=wp)     :: GrdErrN       ! Estimated error in N/S gradient(mm)
    REAL (KIND=wp)     :: GrdErrE       ! Estimated error in E/W gradient(mm)
    REAL (KIND=wp)     :: LogTEC        ! (Log10 of) total electron content
    INTEGER (KIND=iintegers)              :: &
                          NumSlants     ! Number of slants this data sample
    TYPE (COSTslant)   Slant(MaxSlants) ! Define up to MaxSlants slants
  END TYPE
  TYPE (COSTsample)       Sample 

! Input format definitions
  CHARACTER (LEN=5),  PARAMETER :: hfmt1="(A20)"
  CHARACTER (LEN=65), PARAMETER :: hfmt2=           &
            "(A4/A20,5X,A20/2F12.5,2F12.3/2(A20,5X)/4(A20,5X)/3I5/Z8.8/I4)"
  CHARACTER (LEN=3),  PARAMETER :: dfmt1="(A)"
  CHARACTER (LEN=27), PARAMETER :: dfmt2="(3I3,Z9.8,7F7.1,4F7.2,F7.3)"
  CHARACTER (LEN=4),  PARAMETER :: dfmt3="(I4)"
  CHARACTER (LEN=10), PARAMETER :: dfmt4="(A4,4F7.1)"

  CHARACTER (LEN=20)            :: hline      ! header line
  CHARACTER (LEN=110)           :: dline      ! data line
  INTEGER  (KIND=iintegers)     :: nsamp      ! sample loop counter
  INTEGER  (KIND=iintegers)     :: nslan      ! slant  loop counter

  INTEGER  (KIND=iintegers)     :: nconv      ! used for date conversion

  LOGICAL, SAVE                 :: leof  = .false.   ! end-of-file flag

! Local arrays:
! ------------
  INTEGER (KIND=iintegers), ALLOCATABLE :: &
    nodbuf (:)          ,& ! buffer containing all reports for a single node
    iallbuf(:)          ,& ! buffer containing observations from all PEs
    iposbuf(:)          ,& ! position of the reports in the report buffer
    inodbuf(:)             ! node the reports belong to
!
!------------- End of header ---------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine obs_read_gps
!-------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! Section 1.0: Initialization of organizational variables
  !-----------------------------------------------------------------------------

  yroutine = 'obs_read_gps'

  lprcs    = (lverif) .OR. (liwvssc)
  ngpob = ntotgp
! Find maximum length of the buffer containing all reports for all PE's
  nmxbln  = ( maxgpo * nlengps )
  nexceed = 0

  ALLOCATE( nodbuf (nmxbln)             , STAT=izlocstat )
  nodbuf (:)   = 0

  ibuflen      = 0
  npe          = 0
  IF (num_compute > 1) THEN
    ALLOCATE ( iallbuf (nmxbln)           , STAT=izlocstat )
    ALLOCATE ( inodbuf (nmxbln/nlengps+1) , STAT=izlocstat )
    ALLOCATE ( iposbuf (nmxbln/nlengps+1) , STAT=izlocstat )
! (nlengps is the minimum length of a report in the buffer 'iallbuf') 
    iallbuf (:)  = 0
    inodbuf (:)  = 0
    iposbuf (:)  = 0
    nobbuf       = 0
  ENDIF 

! Initial date and time of the forecast run
  READ (ydate_ini,'(I4,5I2)') imoyy, imomm, imodd, imohh, imomin, imosec

  iminmin =  999
  iminmax = -999

  !-------------------------------------------------------------------------
  ! Section 1.1  Compose full fields of model surface height / land fraction
  !-------------------------------------------------------------------------

  ALLOCATE ( hsurf_tot (ie_tot,je_tot)      , STAT=izlocstat )
  ALLOCATE ( fland_tot (ie_tot,je_tot)      , STAT=izlocstat )

  CALL get_global_surface ( num_compute     , my_cart_id                       &
                          , ie     , je     , hhl(:,:,ke+1) , fr_land          &
                          , ie_tot , je_tot , hsurf_tot     , fland_tot )
! =======================

  !-----------------------------------------------------------------------------
  ! Section 2.    Read GPS data
  !-----------------------------------------------------------------------------

  IF (my_cart_id == 0) THEN
! =========================

!CS: open files if reading other obs from NetCDF files
    OPEN (nurej  , FILE=yurejct, FORM='FORMATTED', STATUS='UNKNOWN'            &
                               , POSITION='APPEND', IOSTAT=nstat)
    IF (nstat /= 0) THEN
      yerrmsg = 'OPENING OF FILE yurejct FAILED'
      CALL model_abort (my_cart_id, 1010, yerrmsg, yroutine)
    ENDIF
    IF (nstat == 0)                                                            &
      OPEN (nustat , FILE=yustats, FORM='FORMATTED', STATUS='UNKNOWN'          &
                                 , POSITION='APPEND', IOSTAT=nstat)
    IF (nstat /= 0) THEN
      yerrmsg = 'OPENING OF FILE yustats FAILED'
      CALL model_abort (my_cart_id, 1011, yerrmsg, yroutine)
    ENDIF

    read_gps : DO WHILE (.NOT.leof)
!   -------------------------------

! Read (next) header

    hline   = ' ' 
    istat = 0
    DO WHILE ( hline(1:10) /= SOVFmark .AND. istat == 0 ) 
      READ ( UNIT=nugps, FMT=hfmt1, IOSTAT=istat ) hline
    ENDDO

! If OK, read rest of header info
    IF ( istat == 0 ) THEN
      Header%FmtVer = hline(10:13)
      IF ((Header%FmtVer /= COSTfmtver).AND.(Header%FmtVer /= COSTfmtve2)) THEN
         PRINT *, "COST format version ", Header%FmtVer, " not supported"
         PRINT *, "**** NO GPS DATA READ FROM INPUT FILE ****"
         EXIT read_gps
      ENDIF
 
      READ ( UNIT=nugps, FMT=hfmt2, IOSTAT=istat ) &
                        Header%StnID,         &
                        Header%RcvrType,      &
                        Header%AntType,       &
                        Header%StnLat,        &
                        Header%StnLon,        &
                        Header%StnHt,         &
                        Header%StnHtAMSL,     &
                        Header%DateTimeFirst, &
                        Header%DateTimeProc,  &
                        Header%ProcCentre,    &
                        Header%ProcMethod,    &
                        Header%OrbitType,     &
                        Header%MetSource,     &
                        Header%TimeInc,       &
                        Header%UpdateInt,     &
                        Header%BatchLen,      &
                        Header%PCDH,          &
                        Header%NumSamples
! Ensure strings are upper case

      CALL to_upper ( Header%StnID )
      CALL to_upper ( Header%RcvrType )
      CALL to_upper ( Header%AntType )
      CALL to_upper ( Header%DateTimeFirst )
      CALL to_upper ( Header%DateTimeProc )
      CALL to_upper ( Header%ProcCentre )
      CALL to_upper ( Header%ProcMethod )
      CALL to_upper ( Header%OrbitType )
      CALL to_upper ( Header%MetSource )

! If header ok, read samples
      nsamp = 0
      IF ( Header%NumSamples < 0 ) Header%NumSamples = MaxSamples

      DO WHILE ( nsamp < Header%NumSamples .AND. istat==0 )
        READ ( UNIT=nugps, FMT=dfmt1, IOSTAT=istat ) dline

        IF ( istat == 0 ) THEN
          IF ( dline(1:100) /= EOVFmark ) THEN
            READ (dline, FMT=dfmt2, IOSTAT=istat ) &
                       Sample%TimeStamp, &
                       Sample%PCDD,      &
                       Sample%TZD,       &
                       Sample%TZDErr,    &
                       Sample%ZWD,       &
                       Sample%IWV,       &
                       Sample%SurfPres,  &
                       Sample%SurfTemp,  &
                       Sample%SurfRH,    &
                       Sample%GradientN, &
                       Sample%GradientE, &
                       Sample%GrdErrN,   &
                       Sample%GrdErrE,   &
                       Sample%LogTEC
            nsamp = nsamp + 1        ! (Another) sample read successfully
  
! Next read no. of slant samples and input that number of slants
            READ ( UNIT=nugps, FMT=dfmt3, IOSTAT=istat ) &
                        Sample%NumSlants

            DO nslan = 1, MIN(Sample%NumSlants,MaxSlants)
              READ ( UNIT=nugps, FMT=dfmt4, IOSTAT=istat ) &
                         Sample%Slant(nslan)%SatID,  &
                         Sample%Slant(nslan)%TSD,    &
                         Sample%Slant(nslan)%TSDErr, &
                         Sample%Slant(nslan)%Azim,   &
                         Sample%Slant(nslan)%Elev
            END DO ! on slants

! %%%% no correct allignement below
  !-----------------------------------------------------------------------------
  ! Section 2.1   Pre-selection of GPS reports
  !-----------------------------------------------------------------------------

 
  ! Skip empty report, i.e. if negative or zero (?) TZD
  !-----------------------------------------------------------------------
    IF (Sample%TZD <= c0)  CYCLE read_gps

    nkflag  = 0

  ! Station id
  !-----------------------------------------------------------------------
    WRITE( ystid(1:4),'(A4)' ) Header%StnID
    DO icl = 5, ilstidp
      ystid (icl:icl) = ' '
    ENDDO

  ! Observation and code type
  !-----------------------------------------------------------------------
    nobtyp  = ngps
    ncdtyp  = ngpgfz

!   get diagnostic array position
!   -----------------------------
    icma  = i_cma ( nobtyp , ncdtyp )
    cma(icma)%cnt_pr = cma(icma)%cnt_pr + 1
!   pre-set passive flag to zero
    mpassiv = 0
    mrepflg = 0

  ! blacklists (to be done)
  !-----------------------------------------------------------------------

  ! Observation time
  !-----------------------------------------------------------------------

! Extract date from header
! Year
    READ(Header%DateTimeFirst(8:11),'(I4)') nconv
!CS: time cheat
!   nconv = 2009
!CS: time cheat
    iobdat = nconv * 10000
! Month (trasform month name in month number)
    CALL convert_month ( Header%DateTimeFirst(4:6), nconv , 1)
!CS: time cheat
!   nconv = 2
!CS: time cheat
    iobdat = iobdat + nconv*100
! Day
    READ(Header%DateTimeFirst(1:2),'(I2)') nconv
    iobdat = iobdat + nconv
! Extract time from sample
    ntime  = Sample%TimeStamp(1)*100 + Sample%TimeStamp(2)

    WRITE( yobdat,'(I8, I4)' )  iobdat, ntime
    READ ( yobdat,'(I4,4I2)' )  ianyy, ianmm, iandd, ianhh, ianmin

!   Time difference in minutes between observation time and initial model time
    CALL diff_minutes ( imoyy, imomm, imodd, imohh, imomin,                    &
                        ianyy, ianmm, iandd, ianhh, ianmin, 0 , imindif, ierrf )
!   =================

!   observation time in forecast hour units  
    zobtim   = imindif/60.0_wp
    iminmin  = MIN( iminmin , imindif )
    iminmax  = MAX( iminmax , imindif )

!   reject report if too old
!   ------------------------
    timlim  =  zobtim + MAX( wtuksue , tipmxsu )
!   PRINT *,'GPS yobdat ', yobdat, zobtim, timlim
    IF (ntstep*dt/3600._wp > timlim ) THEN
      neventr (netime,icma) = neventr(netime,icma) + 1
      cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
      WRITE(nurej,'(" GPS STA ",A ," : OBSERVATION TOO OLD: ",F6.1,            &
                &" [FORECAST HRS]")' ) ystid, zobtim
      CYCLE read_gps
    ENDIF

  ! Observation position
  !-----------------------------------------------------------------------

!   observation height (above Mean Sea Level)
!   -----------------------------------------
    nelevt     =   NINT(Header%StnHtAMSL)
    IF ((nelevt < -400) .OR. (nelevt > 9000)) THEN
      neventr (nenoal,icma) = neventr(nenoal,icma) +1
      cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
      WRITE( nurej,'(" STATION ",A ," : GPS STATION HEIGHT MISSING"            &
                   &,F6.1," HRS")' ) ystid, zobtim
      CYCLE read_gps
    ENDIF

!   observation location
!   --------------------
    roblat     =   Header%StnLat
    IF ( Header%StnLon > 180._wp ) THEN
      roblon =  Header%StnLon - 360._wp
    ELSE
      roblon =  Header%StnLon
    ENDIF

    rlon = rla2rlarot (roblat, roblon, pollat, pollon, polgam)
    rlat = phi2phirot (roblat, roblon, pollat, pollon)
!   observation position in grid point units (total area)
    zio_tot    = 1._wp + (rlon - startlon_tot) /dlon
    zjo_tot    = 1._wp + (rlat - startlat_tot) /dlat

!   Treat GPS as synop observation to
!   assign observation to grid point
!   --------------------------------
    nobtyp  = nsynop
    lseaobs =    .FALSE.

    CALL obs_assign_gridpt ( zio_tot, zjo_tot, nelevt, .FALSE., .FALSE.        &
                           , ie_tot, je_tot, hsurf_tot, fland_tot              &
                           , dlat, dlon, startlat_tot, degrad                  &
                           , doromx, fdoro, imdi, nboundlines , nolbc          &
                           , lseaobs , iob_tot, job_tot, zsurf, lsfcob )
!   ======================
!   reassign correct observation type
    nobtyp  = ngps
    IF (lseaobs)  mrepflg = IBSET( mrepflg , nvsebp )

!   reject report if observation location outside model domain
!   ----------------------------------------------------------
    IF   (iob_tot == 0) THEN
      WRITE(nurej,'(" GPS STA ",A ," : OBS. LOCATION OUT OF DOMAIN ",2F7.1     &
                  &," , ",F6.1," HRS")' ) ystid, roblon, roblat, zobtim
      IF (mpassiv == 0) neventr (neloca,icma) = neventr (neloca, icma) + 1
      IF (mpassiv >  0) cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
      cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
      CYCLE read_gps
    ENDIF

!   PRINT '("ZG ",A,I10,2x,L,2I5,2x,L,F5.0,2F7.2)', ystid, nelevt, lseaobs     &
!                                           , iob_tot,job_tot, lsfcob, zsurf   &
!                                           , roblat, roblon

!   observation location outside user-defined area
!   ----------------------------------------------
    IF (     (roblat > obnlat) .OR. (roblat < obslat)                          &
        .OR. ((obwlon >= obelon) .AND. (roblon < obwlon)                       &
                                 .AND. (roblon > obelon))                      &
        .OR. ((obwlon <  obelon) .AND. (     (roblon < obwlon)                 &
                                        .OR. (roblon > obelon)))) THEN
      WRITE( nurej,'(" GPS STA ",A ," : OBS. LOCATION OUT OF USER-SPECIFIED"   &
                   &," AREA ",2F7.1," , ",F6.1," HRS")' )                      &
             ystid, roblon, roblat, zobtim
      IF (mpassiv == 0) neventr (neloca,icma) = neventr (neloca,icma) + 1
      IF (mpassiv >  0) cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
      IF (.NOT. lprcs ) cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
      IF (.NOT. lprcs ) CYCLE read_gps
 
      cma(icma)%cnt_ps = cma(icma)%cnt_ps + 1
      mpassiv = 2
      mrepflg = IBSET( mrepflg , nvobbp )
      nkflag  = IBSET( nkflag  , FL_AREA )
    ENDIF

!   distance "model orography - station altitude" too large
!     (but it's ok if GPS obs is over sea)
!   -------------------------------------------------------
!   as for SYNOP
    IF (iob_tot >  0) THEN
      fisd   = REAL( nelevt, wp ) - zsurf
      fisd = ((fdoro(2)-c1)/c2 + SIGN( (fdoro(2)+c1)/c2 , fisd) ) * fisd
      IF ((nelevt > NINT(altopsu(2))) .OR. (fisd > doromx(2)))                 &
        iob_tot = - iob_tot
    ENDIF
!   IF ((iob_tot <  0) .OR. (lseaobs)) THEN
    IF (iob_tot <  0) THEN
      WRITE( nurej,'(" GPS STA ",A ," : HEIGHT ",I5," DIFF. TO MODEL ",F5.0    &
                   &," TOO LARGE, ",F7.1," HRS")') ystid, nelevt, zsurf, zobtim
      IF (mpassiv == 0) neventr (nezdif,icma) = neventr (nezdif,icma) + 1
      IF (mpassiv >  0) cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
      IF (.NOT. lprcs ) cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
      IF (.NOT. lprcs ) CYCLE read_gps

      cma(icma)%cnt_ps = cma(icma)%cnt_ps + 1
      mpassiv = 2
      mrepflg = IBSET( mrepflg , nvalbp )
      nkflag  = IBSET( nkflag  , FL_HEIGHT )
      iob_tot = ABS( iob_tot )
      job_tot = ABS( job_tot )
    ENDIF

  ! Observation type, code type
  !-----------------------------------------------------------------------

    IF (     (nobtyp == ngps   .AND. .NOT.lgps  )                              &
        .OR. (ncdtyp == ngpgfz .AND. .NOT.lcd096)) THEN
      IF (      (roblat <= exnlat+epsy) .AND. (roblat >= exslat-epsy)          &
          .AND. (roblon >= exwlon-epsy) .AND. (roblon <= exelon+epsy)) THEN
        IF (mpassiv == 0) neventr (neobct,icma) = neventr(neobct,icma) + 1
        IF (mpassiv >  0) cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
        IF (.NOT. lprcs ) cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
        IF (.NOT. lprcs ) CYCLE read_gps
        cma(icma)%cnt_ps = cma(icma)%cnt_ps + 1
        mpassiv = 2
        mrepflg = IBSET( mrepflg , nvexbp )
        nkflag  = IBSET( nkflag  , FL_OBSTYPE )
      ENDIF
    ENDIF

!csc
    IF (TRIM(Header%ProcCentre) == 'GFZ') THEN
      ncentre = 23023
    ELSE
      ncentre = imdi
    ENDIF
    nstachr = 0

  !-----------------------------------------------------------------------
  ! Section  2.2      Find the appropriate node for this report
  !-----------------------------------------------------------------------

    irnode = imdi
    Nodes: DO  inode = 0, num_compute-1
      IF ( isubpos(inode,1) <= iob_tot .AND.                                   &
           isubpos(inode,3) >= iob_tot .AND.                                   &
           isubpos(inode,2) <= job_tot .AND.                                   &
           isubpos(inode,4) >= job_tot ) THEN
        irnode = inode
        EXIT Nodes
      ENDIF
    ENDDO Nodes
    IF (irnode == imdi) THEN
 
      IF (mpassiv == 0) neventr (neloca,icma) = neventr(neloca,icma) + 1
      IF (mpassiv >  0) cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
      cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
      WRITE( nurej,'(" GPS STA ",a ," : NO APPROPRIATE NODE FOUND ",           &
                   &2F7.1," , ",F6.1," HRS")' ) ystid, roblon, roblat, zobtim
      CYCLE read_gps
    ENDIF

!   count accepted reports
    IF (mpassiv == 0) cma(icma)%cnt_ac = cma(icma)%cnt_ac + 1
  !----------------------------------------------------------------------
  ! Section  2.3  Store this report in buffer 'nodbuf' ( --> 'iallbuf')
  !               together with the already extracted header information
  !----------------------------------------------------------------------

! passive reports are rejected, if lverpas = .FALSE.
    IF ((.NOT. lverpas) .AND. (mpassiv > 0)) THEN
      cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
      cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1  
      CYCLE read_gps
    ENDIF

!  GPS observations to be saved 
    tzderr    = Sample%TZDErr
    tzd       = Sample%TZD
    zwd       = Sample%ZWD
    iwv       = Sample%IWV
    surfpres  = Sample%SurfPres
    surftemp  = Sample%SurfTemp
    surfrh    = Sample%SurfRH
!CS
!   PRINT '(" OBS_READ_GPS: ",A,2F6.1,4F8.1,3I5,F5.0)',  &
!   ystid,zobtim,iwv,tzd,surfpres,surftemp,surfrh,iob_tot,job_tot,nelevt,zsurf

!   actual buffer length of array 'nodbuf'
    ipos   = ibuflen
!   report length of current observation: 
!          header information + values in GPS report buffer gpsbuf
    irplen = nlengps 
    IF (ipos+irplen > nmxbln)       THEN
      nexceed  =  nexceed + 1
      WRITE( nustat ,'(A20,": CAUTION: maxgpo TOO SMALL ==> REPORT REJEC"      &
                     &,"TED", / ,26X,"SINCE ACTUAL / MAXIMUM BUFFER LENGTH:"   &
                     &,I7," / ",I7)' )                  yroutine, ipos, nmxbln
      PRINT          '(A20,": CAUTION: maxgpo TOO SMALL ==> REPORT REJEC"      &
                     &,"TED", / ,26X,"SINCE ACTUAL / MAXIMUM BUFFER LENGTH:"   &
                     &,I7," / ",I7)' ,                  yroutine, ipos, nmxbln
! insufficient ODR size is the only report event, which is updated even for
! reports set to passive previously
!     IF (mpassiv == 0)                                                        &
      neventr (nesodr,icma) = neventr(nesodr,icma) + 1
      IF (mpassiv == 0) cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
      IF (mpassiv >  0) cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
      cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
      CYCLE read_gps
    ENDIF
    IF (mpassiv == 2) mrepflg = IBSET( mrepflg , nvpsbp )
    nodbuf (ipos+ 1) = irplen
    nodbuf (ipos+ 2) = nobtyp
    nodbuf (ipos+ 3) = ncdtyp
    nodbuf (ipos+ 4) = nstachr
    nodbuf (ipos+ 5) = nkflag
    nodbuf (ipos+ 6) = ncentre
    nodbuf (ipos+ 7) = iobdat
    nodbuf (ipos+ 8) = ntime
    nodbuf (ipos+ 9) = imindif
    nodbuf (ipos+10) = nelevt
    nodbuf (ipos+11) = NINT( zsurf   *c1000 )
    nodbuf (ipos+12) = NINT( zio_tot *c1000 )
    nodbuf (ipos+13) = NINT( zjo_tot *c1000 )
    nodbuf (ipos+14) = NINT( roblon  *c1000 )
    nodbuf (ipos+15) = NINT( roblat  *c1000 )
    nodbuf (ipos+16) = iob_tot
    nodbuf (ipos+17) = job_tot
    nodbuf (ipos+18) = mrepflg
    nodbuf (ipos+19) = NINT( tzderr  *c1000 )
    nodbuf (ipos+20) = NINT( tzd     *c1000 )
    nodbuf (ipos+21) = NINT( zwd     *c1000 )
    nodbuf (ipos+22) = NINT( iwv     *c1000 )
    nodbuf (ipos+23) = NINT( surfpres*c1000 )
    nodbuf (ipos+24) = NINT( surftemp*c1000 )
    nodbuf (ipos+25) = NINT( surfrh  *c1000 )
    nodbuf (ipos+26) = ICHAR(ystid(1:1))
    nodbuf (ipos+27) = ICHAR(ystid(2:2))
    nodbuf (ipos+28) = ICHAR(ystid(3:3))
    nodbuf (ipos+29) = ICHAR(ystid(4:4))

    IF (num_compute > 1) THEN
!     Update number of reports
      nobbuf           = nobbuf + 1
!     Node the report belongs to
      inodbuf (nobbuf) = irnode
!     Position of the report in 'nodbuf'
      iposbuf (nobbuf) = ipos
    ENDIF
!   Update buffer length
    ibuflen          = ibuflen + irplen

! %%%% no correct alligment above

          ELSE
            istat = EOVFstat

          END IF  ! on dline

        END IF    ! istat on Sample
      END DO    ! on Samples

    ELSEIF ( istat < 0 ) THEN
      leof = .true.         ! end-of-file

    ELSEIF ( istat > 0 ) THEN
!   I/O error
      WRITE( nurej,'(" I/O ERROR DETECTED ON GPS FILE")')
      yerrmsg = 'I/O ERROR DETECTED ON GPS FILE'
      CALL model_abort (my_cart_id, 1001, yerrmsg, yroutine)

    ENDIF ! istat of header

    ENDDO read_gps
!   --------------

    !   diagnostic printout on standard output
    inobs = ibuflen / nlengps 
    PRINT '(" processing",I6," reports from",I5," - ",I4," [min] from file ",A)' &
          ,inobs, iminmin, iminmax, yugps(1:LEN_TRIM(yugps))

  !----------------------------------------------------------------------
  ! Section  2.4  Sort reports according to the nodes they belong to
  !----------------------------------------------------------------------

    IF (num_compute > 1) THEN
      iallen      = 0
      inodlen (:) = 0
      Sort: DO inode = 0 , num_compute-1
        DO inobs = 1 , nobbuf
          IF (inodbuf(inobs) == inode) THEN
            ipos = iposbuf (inobs)
            ilen = nodbuf (ipos+1)
            DO ii = 1 , ilen
              iallbuf (iallen+ii) = nodbuf (ipos+ii)
            ENDDO
            inodlen (inode+1) = inodlen(inode+1) + ilen
            iallen            = iallen           + ilen
          ENDIF
        ENDDO
      ENDDO Sort
    ENDIF  ! num_compute > 1

    CLOSE (nurej)
    CLOSE (nustat)

    IF (nexceed > 0) THEN
      OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'          &
                                 , POSITION='APPEND', IOSTAT=istat)
      IF (istat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
      IF (istat /= 0) CALL model_abort (my_cart_id, 1408, yerr, yroutine)
      WRITE( nucautn,'("CAUTION !!!!! t=",I5,": TOO SMALL BUFFER LENGTH FOR"   &
                     &," READING GPS:",I9)' ) ntstep, nmxbln
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxgpo BY AT LEAST" &
                     &,I6)' ) nexceed
      CLOSE (nucautn)
    ENDIF
  ENDIF ! my_cart_id == 0 
 
  !----------------------------------------------------------------------
  ! Section  2.5  Scatter the reports to the nodes they belong to
  !----------------------------------------------------------------------

  IF (num_compute > 1) THEN

! Broadcast the number of elements for each node

    CALL distribute_values ( inodlen, num_compute, 0, imp_integers, icomm_cart &
                           , izmplcode )
!   ======================

! Determine the displacement of the first of the elements for each node
! within array 'iallbuf', the number of elements for the current node

    inodpos (1) = 0
    DO inode = 2 , num_compute
      inodpos (inode) = inodpos(inode-1) + inodlen(inode-1)
    ENDDO
    inodcnt = inodlen(my_cart_id+1)

! Determine 'nodbuf' by scattering the elements to the nodes they belong to

    nodbuf (:) = 0
    IF (inodpos(num_compute)+inodlen(num_compute) > 0) THEN

      CALL scatterv_values ( iallbuf, nmxbln, 0, inodlen                       &
                           , nodbuf , nmxbln, my_cart_id, num_compute          &
                           , imp_integers, icomm_cart, yerrmsg, izmplcode )
!     ====================

      IF (izmplcode /= 0)                                                      &
        CALL model_abort (my_cart_id, 11177, yerrmsg, yroutine, izmplcode)
    ENDIF

  ENDIF  ! num_compute > 1


  IF (lwonl)                                                                   &
    WRITE(nupr,'(" OBS_READ_GPS:  length of first report received",            &
               &" at node number ",i3," :  ",i4," words")') my_cart_id,nodbuf(1)

  !--------------------------------------------------------------------------
  ! Section 3.  Deallocate local arrays
  !--------------------------------------------------------------------------

  IF (num_compute > 1) THEN
    DEALLOCATE ( iallbuf   , STAT = istat )
    DEALLOCATE ( iposbuf   , STAT = istat )
    DEALLOCATE ( inodbuf   , STAT = istat )
  ENDIF
  DEALLOCATE ( hsurf_tot   , STAT = istat )
  DEALLOCATE ( fland_tot   , STAT = istat )

  !--------------------------------------------------------------------------
  ! Section 4.  Fill GPS observation data arrays
  !--------------------------------------------------------------------------

  nexceed = 0

! start index -1 of first observation in report buffer

  ipos = 0

  get_gps_report :   DO
! ~~~~~~~~~~~~~~~~~~~~~

! get length of next report
    IF  (ipos +1 > nmxbln)                                   EXIT get_gps_report
    irplen   = nodbuf(ipos+1)
    IF  (irplen == 0)                                        EXIT get_gps_report

! read GPS buffer
    irplen   = nodbuf(ipos + 1)
    nobtyp   = nodbuf(ipos + 2)
    ncdtyp   = nodbuf(ipos + 3)
    nstachr  = nodbuf(ipos + 4)
    nkflag   = nodbuf(ipos + 5)
    ncentre  = nodbuf(ipos + 6)
    iobdat   = nodbuf(ipos + 7)
    ntime    = nodbuf(ipos + 8)
    zobtim   = REAL ( nodbuf(ipos + 9), wp ) /60._wp
    nelevt   = nodbuf(ipos +10)
    zsurf    = REAL ( nodbuf(ipos +11), wp ) * c1000r
    zio_tot  = REAL ( nodbuf(ipos +12), wp ) * c1000r
    zjo_tot  = REAL ( nodbuf(ipos +13), wp ) * c1000r
    roblon   = REAL ( nodbuf(ipos +14), wp ) * c1000r
    roblat   = REAL ( nodbuf(ipos +15), wp ) * c1000r
    iob_tot  = nodbuf(ipos +16)
    job_tot  = nodbuf(ipos +17)
    mrepflg  = nodbuf(ipos +18)
    mpassiv  = IBITS( mrepflg , nvpsbp , 1 ) * 2
    tzderr   = REAL ( nodbuf(ipos +19), wp ) * c1000r
    tzd      = REAL ( nodbuf(ipos +20), wp ) * c1000r
    zwd      = REAL ( nodbuf(ipos +21), wp ) * c1000r
    iwv      = REAL ( nodbuf(ipos +22), wp ) * c1000r
    surfpres = REAL ( nodbuf(ipos +23), wp ) * c1000r
    surftemp = REAL ( nodbuf(ipos +24), wp ) * c1000r
    surfrh   = REAL ( nodbuf(ipos +25), wp ) * c1000r
    ystid(1:1)=ACHAR( nodbuf(ipos +26) )
    ystid(2:2)=ACHAR( nodbuf(ipos +27) )
    ystid(3:3)=ACHAR( nodbuf(ipos +28) )
    ystid(4:4)=ACHAR( nodbuf(ipos +29) )
    DO icl = 5, ilstidp
      ystid (icl:icl) = ' '
    ENDDO

! Assigned grid point in the processor domain
    iob     = iob_tot - isubpos (my_cart_id,1) + 1 + nboundlines
    job     = job_tot - isubpos (my_cart_id,2) + 1 + nboundlines

! check whether there is space in the ODR
! ---------------------------------------
    IF (ngpob >= maxgpl) THEN
      nexceed  =  nexceed + 1
      IF (lwonl)                                                               &
        WRITE( nupr,'(" CAUTION !!!!! STATION ",A ,": ",I5,                    &
                     &"th GPS REPORT EXCEEDS ODR SIZE MAXGPL ",I5)')           &
               ystid, ngpob, maxgpl
      PRINT         '(" CAUTION !!!!! STATION ",A ,": ",I5,                    &
                     &"th GPS REPORT EXCEEDS ODR SIZE MAXGPL")' ,              &
               ystid, ngpob
      icma  = i_cma ( nobtyp , ncdtyp )
! insufficient ODR size is the only report event, which is updated even for
! reports set to passive previously
!     IF (mpassiv == 0)                                                        &
         neventr (nesodr,icma) = neventr(nesodr,icma) + 1
      IF (mpassiv == 0) cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
      IF (mpassiv >  0) cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
      cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
      ipos = ipos + irplen
      CYCLE get_gps_report
    ENDIF

! create ODR station characteristic word, incl. flags and instr. specification
! ----------------------------------------------------------------------------
    nstachr = mrepflg
    IF (iwv >  epsy)  nbufbit = 36
    IF (iwv <= epsy)  nbufbit = 38
!   IF ((iwv >  epsy) .AND. (.NOT. lgpsbias))  nbufbit = 36
!   IF ((iwv >  epsy) .AND. (      lgpsbias))  nbufbit = 37
!   IF ((iwv <= epsy) .AND. (.NOT. lgpsbias))  nbufbit = 38
!   IF ((iwv <= epsy) .AND. (      lgpsbias))  nbufbit = 39
!   nstachr = insert( nstachr , nbufbit , nvinbp )
    nstachr = IOR( nstachr , ISHFT( nbufbit , nvinbp ) )

! fill GPS observation data arrays
!------------------------------------------------------------------------------

    ngpob = ngpob + 1
    ogphed (ngpob,nhilon) = roblon
    ogphed (ngpob,nhjlat) = roblat
    ogphed (ngpob,nhalt ) = REAL ( nelevt, wp )
    ogphed (ngpob,nhtime) = zobtim
    ogphed (ngpob,nhsurf) = zsurf
    ogphed (ngpob,nhzio ) = zio_tot
    ogphed (ngpob,nhzjo ) = zjo_tot
 
    mogphd (ngpob,nhio  ) = iob
    mogphd (ngpob,nhjo  ) = job
    mogphd (ngpob,nhitot) = iob_tot
    mogphd (ngpob,nhjtot) = job_tot
    mogphd (ngpob,nhobtp) = nobtyp
    mogphd (ngpob,nhcode) = ncdtyp
    mogphd (ngpob,nhschr) = nstachr
    mogphd (ngpob,nhpass) = mpassiv
    mogphd (ngpob,nhqcfw) = -1

    mogphd (ngpob,nhflag) = nkflag
!   mogphd (ngpob,nhcorr) = imdi
    mogphd (ngpob,nhcorr) = 0
    mogphd (ngpob,nhcat ) = 0
    mogphd (ngpob,nhcats) = 99
    mogphd (ngpob,nhkz  ) = imdi
    mogphd (ngpob,nhcent) = ncentre
    mogphd (ngpob,nhstid) = imdi
    mogphd (ngpob,nhdate) = iobdat
    mogphd (ngpob,nhhrmn) = ntime
    mogphd (ngpob,nhsyhr) = MOD( iobdat,1000000 ) *100 + ntime/100
 
    yogphd (ngpob) = ystid
    DO icl = ilstidp, ilstid
      yogphd(ngpob) (icl:icl) = ' '
    ENDDO
 
! Write values if  >= 0 on the ODR (tzd has been already checked)
    ogpbdy (ngpob,nbgzpd) = tzd
    IF ( tzderr   >= c0 )  ogpbdy (ngpob,nbgtze) = tzderr
    IF ( zwd      >= c0 )  ogpbdy (ngpob,nbgzwd) = zwd
    IF ( iwv      >= c0 )  ogpbdy (ngpob,nbgiwv) = iwv
    IF ( surfpres >= c0 )  ogpbdy (ngpob,nbgp)   = surfpres *100._wp
    IF ( surftemp >= c0 )  ogpbdy (ngpob,nbgt)   = surftemp
! Temperature before March 2001 was in C
!   IF ( surftemp >= -99._wp ) ogpbdy (ngpob,nbgt)   = surftemp
    IF ( surfrh   > epsy)  ogpbdy (ngpob,nbgrh)  = surfrh /100._wp
 
    mogpbd (ngpob,nbgflg) = 0
    mogpbd (ngpob,nbgqcf) = 0
    mogpbd (ngpob,nbglid) = IBSET( 0, nvlidp(7) )
!CS
    mogpbd (ngpob,nbgerr) = 0
    IF ( tzderr   >= c0 )  mogpbd (ngpob,nbgerr) = IBSET( 0, nvriwv )
    ogpbdy (ngpob,nbgbia) = c0
!CS

! Update ipos for next report
    ipos = ipos + irplen

  !--------------------------------------------------------------------------
  ! Section 6.  Redundancy check
  !--------------------------------------------------------------------------

    rscale   =  100._wp*(1._wp+epsy)
    rdegkm2  =  111._wp **2
    cosolat  =  COS( (startlat_tot+(job_tot-1)*dlat) * degrad )

!  Preset ODR index for check against all prior reports
    ngplo  = 0

!  Get next report as 'second' report
!  ----------------------------------
    next_gps_report: DO
!   ~~~~~~~~~~~~~~~~~~~
      ngplo  = ngplo + 1
!  Leave 'loop' if no redundancy with present report
      IF (ngplo == ngpob)                                   EXIT next_gps_report
!  Check if current / second report is (partly) redundant
!  ------------------------------------------------------

!  Check time difference
      lcoloc = (ABS( ogphed(ngplo,nhtime) - ogphed(ngpob,nhtime) )  <=  rtmlim)
      IF (.NOT. lcoloc)                                    CYCLE next_gps_report
!  Check horizontal (and vertical) (quasi-) colocation and station id
      rhzob  =   (cosolat *(mogphd(ngplo,nhitot) -iob_tot) * dlon)**2          &
               +          ((mogphd(ngplo,nhjtot) -job_tot) * dlat)**2
      lcoloc =      (rhzob *rdegkm2 <= rhzlim *rhzlim)                         &
              .AND. (ABS( ogphed(ngplo,nhalt)-ogphed(ngpob,nhalt) ) <=  rvtlim)
      IF (.NOT. lcoloc)                                    CYCLE next_gps_report
      IF (.NOT. (yogphd(ngplo) == yogphd(ngpob)))          CYCLE next_gps_report
!  Decide which report is redundant
!  --------------------------------

      istcor   =  0
      IF ((mogphd(ngpob,nhpass) < 2) .AND. (mogphd(ngplo,nhpass) == 2)) THEN
        nodract  =  ngpob
      ELSEIF ((mogphd(ngpob,nhpass) == 2) .AND. (mogphd(ngplo,nhpass) < 2)) THEN
        nodract  =  ngplo
      ELSEIF (mogphd(ngpob,nhcorr) > mogphd(ngplo,nhcorr)) THEN
        nodract  =  ngpob
        istcor   =  1
      ELSEIF (mogphd(ngpob,nhcorr) < mogphd(ngplo,nhcorr)) THEN
        nodract  =  ngplo
        istcor   =  1
      ELSEIF ((tzderr > epsy) .AND. (     (ogpbdy(ngplo,nbgtze) > tzderr)      &
                                     .OR. (ogpbdy(ngplo,nbgtze) < rmdich))) THEN
        nodract  =  ngpob
      ELSE
        nodract  =  ngplo
      ENDIF
      nodrrej  =  ngpob + ngplo - nodract

      nobtypr  =  mogphd(nodrrej,nhobtp)
      ncdtypr  =  mogphd(nodrrej,nhcode)
      mpassiv  =  mogphd(nodrrej,nhpass)
      IF (lprodr)   ystidr   =  yogphd(nodrrej)

!  Reject redundant data
!  ---------------------

      IF ((lprcs ) .AND. (lverpas)) THEN
        mogphd (nodrrej,nhpass)  = 2
        mogphd (nodrrej,nhschr)  = IBSET( mogphd(nodrrej,nhschr) , nvrdbp )
        mogphd (nodrrej,nhschr)  = IBSET( mogphd(nodrrej,nhschr) , nvpsbp )
      ELSE
        IF (nodract == ngpob) THEN
          ogpbdy (ngplo ,1:mxgbdy) = ogpbdy (nodract,1:mxgbdy)
          mogpbd (ngplo ,1:mxgbdf) = mogpbd (nodract,1:mxgbdf)
          ogphed (ngplo ,1:mxghed) = ogphed (nodract,1:mxghed)
          mogphd (ngplo ,1:mxghdf) = mogphd (nodract,1:mxghdf)
          yogphd (ngplo )          = yogphd (nodract)
        ENDIF
        ogpbdy (ngpob ,1:mxgbdy) = rmdi
        mogpbd (ngpob ,1:mxgbdf) = imdi
        ogphed (ngpob ,1:mxghed) = rmdi
        mogphd (ngpob ,1:mxghdf) = imdi
        yogphd (ngpob )          = '        '
      ENDIF

!  Update statistics and report events and print message
!  -----------------------------------------------------

      icma  = i_cma ( nobtypr , ncdtypr )
      IF (mpassiv < 2) THEN
        cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
        neventr (neredn,icma) = neventr(neredn,icma) + 1
      ELSE
        cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
      ENDIF
      IF ((lprcs ) .AND. (lverpas)) THEN
        cma(icma)%cnt_ps = cma(icma)%cnt_ps + 1
      ELSE
        cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
      ENDIF

      IF (lprodr) THEN
        ilen = 15 + 2*istrej
        IF (nacout+ilen <= nmxoln) THEN
          outbuf(nacout+ 1) = ilen
          outbuf(nacout+ 2) = nfmt10
          outbuf(nacout+ 3) = ncdtypr
          DO icl = 1, 5
            outbuf(nacout+ 3+icl) = imdi
          ENDDO
          IF (ogpbdy(nodract,nbgiwv) > rmdich)                                 &
            outbuf(nacout+ 4) = INT(ogpbdy(nodract,nbgiwv)*rscale)
          IF (ogpbdy(nodract,nbgzwd) > rmdich)                                 &
            outbuf(nacout+ 5) = INT(ogpbdy(nodract,nbgzwd)*rscale)
          IF (ogpbdy(nodract,nbgzpd) > rmdich)                                 &
            outbuf(nacout+ 6) = INT(ogpbdy(nodract,nbgzpd)*rscale)
          IF (ogpbdy(nodract,nbgtze) > rmdich)                                 &
            outbuf(nacout+ 7) = INT(ogpbdy(nodract,nbgtze)*rscale)
          IF (ogpbdy(nodract,nbgp  ) > rmdich)                                 &
            outbuf(nacout+ 8) = INT(ogpbdy(nodract,nbgp  )*rscale)
          outbuf(nacout+ 9) = mogphd(nodract,nhcode)
          outbuf(nacout+10) = INT(ogphed(nodract,nhilon)*rscale)
          outbuf(nacout+11) = INT(ogphed(nodract,nhjlat)*rscale)
          outbuf(nacout+12) = 1 + MAX( 0_iintegers, MIN( i1, nodract - nodrrej))
          outbuf(nacout+13) = istcor
          outbuf(nacout+14) = 0
          outbuf(nacout+15) = INT(ogphed(nodract,nhtime)*rscale)
 
          DO icl = 1 , istrej
            outbuf(nacout+15       +icl) = ICHAR( ystidr          (icl:icl) )
            outbuf(nacout+15+istrej+icl) = ICHAR( yogphd(nodract) (icl:icl) )
          ENDDO
          nacout  = nacout + ilen
        ENDIF
      ENDIF
! Adjust total number of GPS reports

      IF (.NOT. ((lprcs ) .AND. (lverpas)))  ngpob = ngpob - 1
      IF (.NOT. ((lprcs ) .AND. (lverpas)))                 EXIT next_gps_report

    ENDDO next_gps_report
!   ~~~~~~~~~~~~~~~~~~~~~

  ENDDO get_gps_report
! ~~~~~~~~~~~~~~~~~~~~
ntotgp = ngpob
  DEALLOCATE ( nodbuf    , STAT = izlocstat )

! number of GPS obs not used due to insufficient ODR size
! -------------------------------------------------------

  IF (num_compute > 1) THEN
    CALL global_values ( nexceed, 1,'MAX', imp_integers,icomm_cart, 0,yerr,ierr)
!   ==================
    IF (ierr /= 0)  CALL model_abort (my_cart_id, 11014, yerr, yroutine)
!                   ----------------
  ENDIF

  IF ((my_cart_id == 0) .AND. (nexceed > 0)) THEN
    OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'            &
                               , POSITION='APPEND', IOSTAT=istat)
    IF (istat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
    IF (istat /= 0) CALL model_abort (my_cart_id, 1408, yerr, yroutine)
    WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," LOCAL UP-AIR SING-L."      &
                   &," OBS. INCR. BEYOND maxgpl ",I5)' ) ntstep, nexceed, maxgpl
    WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxgpo BY AT LEAST"   &
                   &,I6)' ) nexceed
    CLOSE (nucautn)
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure obs_read_gps
!-------------------------------------------------------------------------------
END SUBROUTINE obs_read_gps


!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_proc_cdf" to check for and open ASCII GPS file
!-------------------------------------------------------------------------------

SUBROUTINE obs_gps_init

!-------------------------------------------------------------------------------
!
! Description :
!   This module procedure initialises the reading from a dedicated ASCII GPS
!   data input file:
!   1. Check existence of GPS input file ('lexigps')
!   2: Open GPS input file, if it exists.
!
! Written by        : DWD, Christoph Schraff (2001)
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!-------------------------------------------------------------------------------
  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Declarations:
!
!-------------------------------------------------------------------------------

! Subroutine arguments: None
! --------------------

! Local variables:
! ---------------

  INTEGER (KIND=iintegers) ::  &
    nstat, ierr, jerr      ! error status variables

  CHARACTER (LEN=70)       ::  &
    yerrmsl                ! error message
  CHARACTER (LEN=20)       ::  &
    yroutine               ! name of this subroutine
  CHARACTER (LEN=40)       ::  &
    yerr                   ! error message

!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine obs_gps_init
!-------------------------------------------------------------------------------

  yroutine = 'obs_gps_init'

  ierr = 0
  jerr = 0
  IF (my_cart_id == 0) THEN
    INQUIRE (FILE=yugps, EXIST=lexigps)
    IF (lexigps) jerr = 1
  ENDIF
  IF (num_compute > 1)                                                         &
    CALL global_values( jerr, 1,'MAX',imp_integers,icomm_cart, -1, yerr,ierr )
!   ==================

  lexigps = (jerr >= 1)

  IF (lexigps) THEN
! GPS input file exists
    IF (my_cart_id == 0) THEN
      OPEN (nugps,FILE=yugps,FORM='FORMATTED',STATUS='OLD',IOSTAT=nstat)
      IF ((nstat /= 0) .OR. (ierr /= 0)) THEN
        WRITE( yerrmsl,'("OPENING OF FILE yugps FAILED: 0 /= ",I5,",",I5)' )   &
               nstat, ierr
        CALL model_abort (my_cart_id, 1014, yerrmsl, yroutine)
      ENDIF
      PRINT *,"OPEN gps COST-716 ASCII file: ", yugps
      REWIND nugps
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure obs_gps_init
!-------------------------------------------------------------------------------
END SUBROUTINE obs_gps_init

!-------------------------------------------------------------------------------

END MODULE src_obs_proc_cdf
