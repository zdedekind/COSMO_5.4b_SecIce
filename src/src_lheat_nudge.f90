!+ Module for latent heat nudging 
!-------------------------------------------------------------------------------

MODULE src_lheat_nudge

!-------------------------------------------------------------------------------
!
! Description:
!   The module "lheat_nudge" performs the latent heat nudging (lhn).
!   The lhn adds temperature increments to the prognostic variable t
!   so that the total temperature increase due to latent heat release
!   in the current timestep corresponds to the amount of analyzed
!   (or observed) precipitation.
!   The temperature increments added due to lhn are derived from the 
!   model heating rate profiles (large scale condensation and convective 
!   heating) scaled by the ratio of analyzed to modelled precipitation
!   (total precipitation: rain and snow from large scale and 
!   convective processes). The analyzed precipitation is based on radar
!   data merged with the model (total) precipitation fields.
!    
!   The module contains as an organizational unit the subroutine
!   "organize_lhn" which is called from the module organize_assimilation.
!   Further module procedures (subroutines) called by organize_lhn:
!   -> lhn_obs_prep : reading+preparing the precip radar data
!      |
!      |--> lhn_obs_open  : open the radar data file and read general header
!      |--> lhn_obs_read  : read a record (header + data) from radar data file
!      |--> distribute_field : distribute field to all PE's
!
!   -> lhn_skill_scores   : verification of precipitation model against radar
!   -> lhn_pr_ana : 'analysis' of precipitation (weighting of model, obs precip)
!
!   -> lhn_t_inc  : derivation of temperature increments by scaling of
!      |           model latent heating profiles 
!      |--> lhn_search    : search for appropriate model heating profile
!      |--> lhn_filt      : vertical filtering of local tinc_lhn profile
!      |--> lhn_limit     : limiting of the tinc_lhn
!      |--> lhn_relax     : horizontal filtering of tinc_lhn
!           |--> hor_filt
!                |--> init_horizontal_filtering
!                |--> horizontal_filtering
!                |--> exchange_boundaries
!   -> lhn_q_inc  : adjust humidity (i.e. qv) to new temperature (t+tinc_lhn)
!   -> or lhn_satad  : do saturation adjustment instead of humidity adjustment
!
!   Note: The names of input/output variables/arrays defined only once 
!         in the module declaration section but used and "filled" by the 
!         different subroutines are documented in the description parts
!         of each procedure for clarity.
!
! -> possibilitiew for speed-up : i
!    include distance weighting (section3 of lhn_obs_prep) as subroutine calls
!       in time weighting (section 2) to avoid IF(obs_dist > radmax) which
!       is redundant to IF(obs(.,.,ilast) >= rmiss .and. obs(.,.,inext)>=rmiss)
!
!
! Current Code Owner: DWD, Klaus Stephan
!    phone: +49  69  8062 2689
!    fax:   +49  69  8062 3721
!    email: klaus.stephan@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! <VERSION>  <DATE>     Christina Ko"pken
!  Initial code
! 3.13       2004/12/03 S. Klink/ K. Stephan
!  Initial code
! 3.14       2005/01/25 Klaus Stephan
!  Update of the scheme
! 3.18       2006/03/03 Klaus Stephan
!  Another update of the scheme
! 3.21       2006/12/04 Klaus Stephan
!  Another update of the scheme
! V3_24        2007/04/26 Klaus Stephan
!  Bug correction for scale_fac: limitation of fac_lhn_up only after LOG
!  Bug correction for summing up temperatures: nnow_nnold must be used
!  Changes in lhn_blacklist_dx_open_and_read: no check for date
!  ktop_lhn, kbot_lhn were replaced by 1 and ke
!  Introduced another section: vertical restriction of increments
!  A change in the diagnostic output
!  Eliminated calls to subroutines hydci_pp_gr and hydci_gr_lhn
!  Eliminated Namelist variables lhn_diagprec, rlhn_scale_dp
! V3_28        2007/07/18 Klaus Stephan
!  lhn_prof_search: new criteria for choosing a point as valuable
!  reference precipitaion: integegration of precipitation flux startes
!                          at layer k which is the first to be greater
!                          or equal than the maximum of the flux within
!                          the column
!  vertical restriction: definition of the uppest layer by indication
!                        the temperature of the layer
!  diagnosis: introduction of rain histogram
!  IO:        bugfix of generation of INPUT-filename (D. Leuenberger)
!             reproducibility check, possibility to ingnore date by
!             namelist noobs_date
! V4_2         2007/12/05 Klaus Stephan
!  - Establishment of possibility to use radar beam height as an additional
!    information for the scheme, switched on/off by lhn_height. An additional grib file
!    is required which contains the radar beam height. It allows a bright band
!    detection and a more sophisticated validation against radar.
!  - Establishment of bright band detection algorithm switch on/off by lhn_bright.
!    Therefore the radar beam height is compared ot the height of the model's
!    freezing level. If the radar measures within the frezzing level and shows locally
!    significantly higher values than in the surroundig the point is marked as bright
!    band and will be unaccounted for LHN.
!  - Enhancement of online verification to the model precipitation at radar beam height
!    (pr_mod_at_dx, pr_mod_at_dx_sum, dbz_at_dx)
!  - minor changes in timing routine
! V4_4         2008/07/16 Ulrich Schaettler, Daniel Leuenberger
!  Eliminated timing variables and other unused variables
!  Enclosed barrier calls with ltime_barrier
!  Changed NL parameter lyear_360 to itype_calendar, to have several options
!  Adapted interface of get_timings
!  Bug Fix for the initialization of blacklist
! V4_5         2008/09/10 Ulrich Schaettler
!  Allocation of blacklist_tot moved to outer subroutine, because values are
!  needed for the whole model run
! V4_8         2009/02/16 Klaus Stephan
!  Elimination of some unused variables
! V4_10        2009/09/11 Davide Cesari
!  Add characters after a backslash in comments; g95 interprets that as
!  continuation line
! V4_12        2010/05/11 Klaus Stephan, Daniel Leuenberger
!  - replaced field dist by field spqual
!  - replaced field wobs_dist by field wobs_space
!  - replaced variable delta_t_obs with lhn_dt_obs
!  - introduced the local variable zlhigh_freq_data (true if radar data is
!    available in high temporal frequency, hard coded to 10min)
!  - new misdat handling: only one missval is now used for the obs and spqual
!    data
!  - more statistics about input data is written to YULHN (taking into account
!    missval handling)
!  - new statistics on spatial weight is written to YULHN
!  - for ease of parsing, warnings are now consistently marked with a
!    leading "WARNING: " in YULHN
!  - bugfix in lhn_sumrad: change 12._ireals by real(ndata), which is the number
!     of datasets read by the model. In case of nudgecast the number will be
!     lower than 12
!  - change simple "nradar" to be a namelist variable
!  - get a CAUTION in case of less radar height information
!  Renamed t0 to t0_melt because of conflicting names (Ulrich Schaettler)
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_15        2010/11/19 Klaus Stephan
!  - Bug correction of bright band detection:
!     - Allocation of brightband now only once at the beginning
!     - change the hight interval, where bright band is considered
!  - minor changes for ktop_temp and ttm_cv
! V4_17        2011/02/24 Ulrich Blahak
!  Adapted interface of exchg_boundaries; corrected kzdims(1:20) -> kzdims(1:24);
!  eliminated my_peri_neigh
! V4_19        2011/08/01 Ulrich Schaettler
!  Introduced conditional compilation for GRIBDWD
! V4_20        2011/08/31 Klaus Stephan
!  Bug correction: To consider both, blacklist flag and bright band flag for
!      calculation of pr_obs and the weight, it should be an .AND. instead of .OR.
! V4_21        2011/12/06 KLaus Stephan
!  Another modification to above bug fix
!  Delete unused variables numnext, numlast
! V4_23        2012/05/10 Ulrich Schaettler
!  Removed switch lprogprec
!  Adapted call to SR distribute_fields (added sender PE) (Uli Schaettler)
! V4_24        2012/06/22 Hendrik Reich, Ulrich Schaettler
!  Adapted length of strings for date variables (HR)
!  Introduced variable izdebug for debug output (to be passed to SR radar_lm_ray)(US)
! V4_25        2012/09/28 Klaus Stephan, Daniel Leuenberger,
!                         Anne Roches, Oliver Fuhrer, Ulrich Blahak
!  Bug fix for variable min0_lhn (should be known at each processor)
!  Changes to enable subhourly model starting points (LETKF approach):
!   - use fully information of new get_utc routine (hh_rad and min0_lhn)
!   - additional call of SR lhn_sumrad if first observation is within current hour (lhm1)
!   - modifications of obs_time, next_obs_time
!  Correction of array size of diagnostic variables with respect to local processors.
!  Replaced qx-variables by using them from the tracer module (AR; OF)
!  UB:
!  Implemented internal switch "l2mom_satads" to be able to switch on/off the 
!    extra saturation adjustments outside the microphysics parts.
!  Added radar_sb_ray for 2-moment microphysics.
!  Corrected interface to radar_lm_ray depending on itype_gscp, added
!   missing qc and multiplication with rho for the input q's.
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  Introduction of MESSy interface
! V4_28        2013/07/12 Klaus Stephan, Ulrich Schaettler
!  Revision of data reading procedure, to avoid the former rewind of the radar data files.
!  Now, the data of one file are read only once, shortly after opening the data file.
!  The time assignment is done in a next step, as well as the identification as precipitation or
!  quality data. The change also affects the SR lhn_sumrad, where no reading has to be
!  performed anymore. In addition this SR has been changed a bit, which is slightly changing the results. 
!    - There is a small inconsistency between former and new SR considering the presicion of 
!      KIND=ireals and KIND=irealgrib. In the former version a observation point was added
!      if "ds_rad(ind) > 0.0_irealgrib".
!      In the new version the observation is already converted to ireals (SR lhn_obs_read), therefore the 
!      condition is changed to "datafield_all(i,j,n) >= 0.0_ireals"
!    - The former normalisation term "ndata" is now space dependend, which better accounted for data outage of
!      single radar stations (ie. "datafield_all(i,j,n) < 0.0_ireals").
!  All input procedures have been extended to work also with grib_api (GRIB1)  (US)
!  Use parameters for vertical grid from module vgrid_refatm_utils (US)
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler, Klaus Stephan
!  Unification of MESSy interfaces and COSMO Tracer structure
!  For the COSMO-Model only use vcoord from vgrid_refatm_utils
!  Set missing value for grib_api (KS)
!  Error and EOF check in lhn_obs_read: if EOF (ilen=0), then do not check the error
!   (otherwise it could go in an endless CYCLE loop).
! V4_30        2013/11/08 Ulrich Schaettler
!  Modifications to read input data also in GRIB2
! V5_1         2014-11-28 Klaus Stephan, Ulrich Schaettler, Oliver Fuhrer
!  - Further modifications to read input data also in GRIB2, i.e. introduce correct short-names
!  - make grib file reading more consistent with other parts of COSMO code and 
!    fixing some bug fixes for correct reading of grib bitmap section
!  - improvement of data file handling of blacklist and height file
!  - set bright band detection to .false. in case of wrong height infomation
!    stored in the height file (ltlhnbright)
!  - bug fixes in determining observation times in case of dt > 60s
!  - bug fix in determining correct day of observation in case of Februar
!  Bugfix of bright band detection in according to verification results. 
!   Now, less points will be detected as bright band, which seems to be more realistic.
!  Correct some minor bug in data file handling:
!   - give a more suitable error message when input file is missing
!   - restrict number of records to be read to expected number of records
!  Adapted interface to function compute_grib_intbuffer_length (US)
!  Changed argument vcoord%vert_coord in SR calhzero to the reference profile hhl_prof (US)
!  Replaced ireals by wp (working precision) (OF)
! V5_3a        2015-11-24 Klaus Stephan
!  Added update of tt_lheat and re-initialization to 0 for timelevel nnew 
!  (was in src_relaxation before, which is now called before the assimilation)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!===============================================================================

! Modules used:

#ifdef GRIBAPI
USE grib_api
#endif

USE data_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    iintegers, & ! KIND-type parameter for standard integer variables
    irealgrib, & ! KIND-type parameter for real variables in the grib library
    iwlength,  & ! length of an integer word in byte
    intgribf,  & ! KIND-type parameter for fortran files in the grib library
    intgribc     ! KIND-type parameter for C files in the grib library

!-------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 1. vertical coordinate parameters and related variables
! -------------------------------------------------------

    hhl_prof,     & ! a special hhl-profile

! 2. horizontal and vertical sizes of the fields and related variables
!---------------------------------------------------------------------
    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ie_tot,       & ! number of grid points in zonal direction
    je_tot,       & ! number of grid points in meridional direction
    ie_max,       & ! Max. of ie on all processors
    je_max,       & ! Max. of je on all processors
    ke,           & ! number of grid points in vertical direction
    ieje,         & ! ie*je
    startlat_tot, & ! transformed latitude of the lower left grid point
    startlon_tot, & ! transformed longitude of the lower left grid point
    dlon,         & ! grid point distance in zonal direction (in degrees)
    dlat,         & ! grid point distance in meridional direction (in degrees)

! 3. start- and end-indices for the computations in the horizontal layers
!------------------------------------------------------------------------
    istart,       & ! start index for the forecast of w, t, qv, qc and pp
    iend,         & ! end index for the forecast of w, t, qv, qc and pp
    jstart,       & ! start index for the forecast of w, t, qv, qc and pp
    jend,         & ! end index for the forecast of w, t, qv, qc and pp
    istartpar,    & ! start index for the forecast of w, t, qv, qc and pp
    iendpar,      & ! end index for the forecast of w, t, qv, qc and pp
    jstartpar,    & ! start index for the forecast of w, t, qv, qc and pp
    jendpar,      & ! end index for the forecast of w, t, qv, qc and pp

! 5. variables for the time discretization and related variables
!----------------------------------------------------------------

    dt,           & ! long time-step
    dt2,          & ! 2 * dt
    klv950,       & ! k index of the LM-mainlevel, on 950 HPa
    klv850,       & ! k index of the LM-mainlevel, on 850 HPa
    klv700,       & ! k index of the LM-mainlevel, on 700 HPa

! 8. Organizational variables to handle the COSMO humidity tracers
! ----------------------------------------------------------------
    idt_qv,  idt_qc,  idt_qi,  idt_qr,  idt_qs,  idt_qg,  idt_qh,  &
             idt_qnc, idt_qni, idt_qnr, idt_qns, idt_qng, idt_qnh

! end of data_modelconfig
!-------------------------------------------------------------------------------

USE data_constants,  ONLY : &

    rdv,          & ! r_d / r_v
    o_m_rdv,      & ! 1 - r_d/r_v
    rvd_m_o,      & !
    t0_melt,      & !
    lh_v,         & !
    cp_d,         & !
    cpdr,         & !
    b1,           & !
    b2w,          & !
    b3,           & !
    b4w,          & !
    b234w,        & !
    repsilon        ! precision of 1.0 in current floating point format

USE data_constants,  ONLY : & ! for radar_lm_ray
    pi,           & ! circle constant
    rho_w,        & ! density of liquid water
    rho_ice,      & ! density of ice          (kg/m^3)
    K_w,          & ! dielectric constant for water
    K_ice           ! dielectric constant for ice

USE data_io, ONLY        : &

! Variables for handling the Gribfile I/O:
    ydate_ini,     & ! start of the forecast (yyyymmddhh (year, month, day, hour, min, sec))
    ydirini,       & ! catalog-name of the file
    yform_read,    & ! format of the (read) files
    idwdednr,      & ! grib edition number for DWD library
    npds,          & ! Dimension for product definition section (pds)
    ngds,          & ! Dimension for grid description section (gds)
    nbms,          & ! Dimension for bit map section (bms)
    nbds,          & ! Dimension for binary data section
    ndsup,         & ! Dimension for dsup
    undefgrib,     & ! value for "undefined" in the grib routines
    ndims            ! Dimension for idims (contains all dimensions)

!-------------------------------------------------------------------------------

USE data_runcontrol, ONLY :  &

! 1. start and end of the forecast
! --------------------------------

    ntstep,       & ! actual time step
    nold,         & ! corresponds to ntstep - 1
    nnow,         & ! corresponds to ntstep
    nnew,         & ! corresponds to ntstep + 1

! 3. controlling the physics
! --------------------------

    lconv,        & ! forecast with convection
    itype_gscp,   & ! type of grid-scale precipitation physics

! 5. additional control variables
! -------------------------------

    l2tls,        & ! time integration by two timelevel RK-scheme (.TRUE.)
                    ! or by default three-time level KW-scheme (.FALSE.)
    lcond,        & ! forecast with condensation/evaporation
    ltime,        & ! detailled timings of the program are given

! 7. additional control variables
! -------------------------------
    l2mom_satads, & ! in case of 2-moment scheme, do all the satads
                    ! (like for the 1-moment schemes), not just the
                    ! satad after the microphysics at the end of the timestep.

! 9. Variables for Ascii file handling, time measuring, ...
! ---------------------------------------------------------
    itype_calendar, & ! for specifying the calendar used
    lperi_x,        & ! if lartif_data=.TRUE.: periodic boundary conditions in x-dir.
                      ! or with Davies conditions (.FALSE.)
    lperi_y,        & ! if lartif_data=.TRUE.: periodic boundary conditions in y-dir.
                      ! or with Davies conditions (.FALSE.)
    l2dim,          & ! 2 dimensional runs

! 12. controlling verbosity of debug output
! -----------------------------------------
    idbg_level,    & ! to control the verbosity of debug output
    ldebug_ass,    & ! if .TRUE., debug output for assimilation
    lprintdeb_all    ! .TRUE.:  all tasks print debug output
                     ! .FALSE.: only task 0 prints debug output

! end of data_runcontrol

USE data_fields, ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------

    p0,             & ! reference pressure at full levels             ( Pa  )
    dp0,            & ! reference pressure thickness of layers        ( Pa  )

! 2. external parameter fields                                        (unit)
! ----------------------------

    hsurf,          & ! height of surface topography                  ( m   )
    hhl,            &

! 3. prognostic variables                                             (unit)
! -----------------------

    t,              & ! temperature                                   (  k  )
    pp,             & ! deviation from the reference pressure         ( pa  )
    rho,            & ! total density of air                          (kg/m3)

! 6. fields that are computed in the parametrization and dynamics     (unit )
!-----------------------------------------------------------------

!   fields of the convection
    tt_conv        ,& ! temperature tendency due to convection        ( K/s  )
    prr_con        ,& ! precipitation rate of rain, convective        (kg/m2*s)  
    prs_con        ,& ! precipitation rate of snow, convective        (kg/m2*s)

!   fields of the precipitation (large scale)
    prr_gsp        ,& ! precipitation rate of rain, grid-scale        (kg/m2*s)
    prs_gsp        ,& ! precipitation rate of snow, grid-scale        (kg/m2*s)
    prg_gsp        ,& ! precipitation rate of graupel, grid-scale     (kg/m2*s)
    prh_gsp        ,& ! precipitation rate of hail, grid-scale        (kg/m2*s)

! for control output at diagnostic grid points only :
    u              ,&
    v              ,&
    w              ,& ! vertical wind speed (defined on half levels)  ( m/s )
    ps             ,& ! surface pressure                              ( pa  )
    clch           ,& ! cloud cover with high clouds                    --
    clcm           ,& ! cloud cover with medium clouds                  --
    clcl           ,& ! cloud cover with low clouds                     --
    clct              ! total cloud cover                               --

! end of data_fields
!-------------------------------------------------------------------------------


USE data_lheat_nudge, ONLY  :   &

! 1. Namelist variables controlling the latent heat nudging
! ---------------------------------------------------------

    llhn         ,& ! main switch for LHN
    llhnverif    ,& ! switch for verification against radar
    lhn_search   ,& ! search for appropriate nearby model heating profile
    lhn_filt     ,& ! vertical filtering of lhn t-increments
    lhn_relax    ,& ! horizontal filtering of lhn t-increments
    lhn_limit    ,& ! impose an absolute limit on lhn t-increments (abs_lhn_lim)
    lhn_hum_adj  ,& ! apply a humidity adjustment along with t-increment
    lhn_spqual   ,& ! switch for the use of a spatial quality function
    lhn_black    ,& ! use blacklist for radar data
    lhn_incloud  ,& ! 
    lhn_diag     ,& ! produce more detailed diagnostic output during lhn
    lhn_qrs      ,& ! calculate integrated precipitation flux
    lhn_logscale ,& ! apply logarithmic scaling factors
    lhn_wweight  ,& ! apply weighting of increments with respect to mean horizontal wind
    
    nlhn_start   ,& ! start of latent heat nudging period in timesteps
    nlhn_end     ,& ! end of latent heat nudging period in timesteps
    nlhnverif_start   ,& ! start of verification of latent heat nudging approach in timesteps
    nlhnverif_end     ,& ! end of verification of latent heat nudging approach in timesteps
    rlhn_search  ,& ! radius (gridpoints) for profiles search (if lhn_search)
    ktop_lhn     ,& ! index for uppest model layer for which lhn is performed
    kbot_lhn     ,& ! index for lowest model layer for which lhn is performed
    ktop_temp    ,& ! temperature of uppest model layer for which lhn is performed

    noobs_date   ,& ! dates without data
    n_noobs      ,& ! number of dates without data

    lhn_dt_obs   ,& ! observational increment (different for DX/PI)
    nradar       ,& ! max. number of radar stations within input data

    nlhn_relax   ,& ! number of interations of horizontal filtering
    abs_lhn_lim  ,& ! absolute limit for lhn t-increments (imposed if lhn_limit)
    fac_lhn_search  ,& ! factor for search nearby profiles
    fac_lhn_up   ,& ! limiting factor for upscaling of model heating profile
    fac_lhn_down ,& ! limiting factor for downscaling of model heating profile
    rad_wobs_lhn ,& ! max. distance to radar for full observation weigh
    rqrsgmax     ,& ! ratio of maximum of qrsgflux, needed for reference precipitation
    radar_in     ,& ! directory for radar input data
    blacklist_file  ,& ! blacklist_file_name

    lhn_coef     ,& ! factor for reduction of lhn t-increments
    thres_lhn    ,& ! threshold of rain rates to be consinderd within lhn approach

    nulhn        ,& ! unit of lhn output file
    yulhn           ! name of lhn output file

USE data_lheat_nudge, ONLY  :   &
    lhn_height   ,& ! use height infos for radar data
    lhn_bright   ,& ! apply bright band detection
    height_file     ! dxheight_file_name

USE data_lheat_nudge, ONLY  :   &

! 2. fields and related variables                                       (unit)
! -------------------------------

    tt_lheat    ,& ! profile of t-increments due to latent heating   ( K/s )
                   ! (stored for current and previous timestep)
    tinc_lhn    ,& ! temperature increments due to lhn               (K/s)
    qrsflux     ,& ! total precipitation flux
    obs         ,& ! observations on model grid at six observation time levels
    spqual      ,& ! spatial quality function on model grid at two observation time levels
!   fields for testing purposes only:
    tminc_lhn   ,& ! cumulated temperature increments due to lhn     (K)
    ttm_lheat   ,& ! cumulated latent heating (grid scale + conv)    (K)
    pr_obs_sum  ,& ! cumulated precipitation (hourly)
    pr_mod_sum  ,& ! cumulated precipitation (hourly)
    pr_ref_sum  ,& ! cumulated precipitation (hourly)
    ttm_cv      ,& ! array for test output of diverse 2D fields
    blacklist   ,& ! blacklist for DX radar data
    dxheight    ,& ! DX radar heights
    brightband  ,& ! array for brightband flacs

    iblock_rad  ,& ! array for gribed data
    ibmap_rad   ,& ! array for bit map
    dsup_rad    ,& ! array for special data
    ds_rad         ! array for unpacked data


! end of data_lheat_nudge
!-------------------------------------------------------------------------------

USE data_parallel, ONLY :   &

! 1. Information about the processor grid
! ---------------------------------------
    ltime_barrier,   & ! if .TRUE.: use additional barriers for determining the
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    imp_integers,    & ! determines the correct INTEGER type used in the model
                       ! for MPI
    imp_logical,     & ! determines the correct LOGICAL type used in the model
                       ! for MPI
    num_compute     ,& !
    my_cart_id      ,& ! rank of this subdomain in the cartesian communicator
    my_cart_neigh,   & ! neighbors of this subdomain in the cartesian grid
    sendbuf,         & ! sending buffer for boundary exchange:
                       ! 1-4 are used for sending, 5-8 are used for receiving
    isendbuflen     ,& ! length of one column of sendbuf


! 2. Further information for MPI
! ------------------------------
    icomm_cart,      & ! communicator for the virtual cartesian topology
    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains
    isubpos            ! positions of the subdomains in the total domain. Given
                       ! are the i- and the j-indices of the lower left and the
                       ! upper right grid point in the order
                       !                  i_ll, j_ll, i_ur, j_ur.
                       ! Only the interior of the domains are considered, not
                       ! the boundary lines.

! end of data_parallel

!-------------------------------------------------------------------------------

USE io_utilities,             ONLY :  &
    open_file, close_file, compute_grib_intbuffer_length

!-------------------------------------------------------------------------------

USE utilities, ONLY :   &
    get_utc_date

!-------------------------------------------------------------------------------

USE vgrid_refatm_utils, ONLY :   &
    vcoord

!-------------------------------------------------------------------------------

USE parallel_utilities, ONLY :   &
    global_values      ,&
    distribute_values  ,& ! subr. to distribute values to all nodes
    distribute_field   ,& ! subr. distributing a total model field to all nodes
    gather_field       ,& ! subr. gathering a total field, sending it to all PEs
    i_global,j_global  ,& ! functions determining global coord. of point i/j
    ij_local              ! function determining at which node to find point i/j

USE parallel_utilities, ONLY :   &
    exchange_profiles     ! subr. for the exchange of some profiles between nodes

!-------------------------------------------------------------------------------

USE meteo_utilities,     ONLY :  satad

!-------------------------------------------------------------------------------

USE pp_utilities,     ONLY :  calhzero, radar_lm_ray
#ifdef TWOMOM_SB
USE pp_utilities,     ONLY :  radar_sb_ray
#endif

!-------------------------------------------------------------------------------

USE environment, ONLY :     &
    exchg_boundaries,       & ! performs the boundary exchange between
                              ! neighboring processors
    model_abort,            & ! aborts the program in case of errors
    comm_barrier              ! explicit synchronization point

!-------------------------------------------------------------------------------

USE time_utilities, ONLY :  get_timings, i_lhn_computations, i_lhn_obs_prep, &
       i_lhn_relax, i_lhn_t_inc, i_lhn_q_inc, i_barrier_waiting_lhn,         &
       i_communications_lhn, i_lhn_search

!-------------------------------------------------------------------------------

USE data_lhn_diag , ONLY :  &
    tt_lheat_o,tinc_lhn_o,ttm_cv_o

!-------------------------------------------------------------------------------

USE src_tracer,     ONLY : trcr_get, trcr_errorstr
USE data_tracer,    ONLY : T_ERR_NOTFOUND

!-------------------------------------------------------------------------------

USE data_obs_lib_cosmo, ONLY: epsy

USE src_lheating,               ONLY :  &
      get_gs_lheating           ! storage of grid scale latent heating for lhn

!===============================================================================

IMPLICIT NONE

!===============================================================================

! Local scalars:
!---------------

  INTEGER (KIND=iintegers) ::  &
    izlocstat         ,& ! error status on allocation of fields
    istat             ,& ! summed errors due to allocation of fields
    izerror           ,& ! error status variable
    ntreat            ,& ! number of grid points to be treated by lhn
    ntreat_tot        ,& ! number of grid points to be treated by lhn
    nnow_nold         ,& ! is set to nold or nnow, depending on the chosen
                         ! integration scheme: leapfrog or 2-timelevel
      
    idummy            

! variable for GRIB-input
! first for obs_read_dx
  INTEGER (KIND=intgribc)  ::  &
    maxlenc, ilenc, ierrc ! variables for C routines

  INTEGER (KIND=intgribf)  ::  &
    ipds  (npds)  ,& ! product definition section for radardata GRIB file
    igds  (ngds)  ,& ! grid description section for radardata GRIB file
    ibms  (nbms)  ,& ! bit map section for radardata GRIB file
    ibds  (nbds)  ,& ! binary data section for radardat GRIB file
    lds           ,&
    lbm           ,&
    lfd

  INTEGER (KIND=intgribf)            ::       &
    idims (ndims)    ! array for all dimensions

  REAL  (KIND=wp)                  ::           &
    zdt                  ,& ! valid time step for integration
    radmax=200.0_wp      ,& ! maximal radius of data coverage around radar ( km )
    sec_per_hr=3600.0_wp ,& ! seconds per hour
    sec_per_hr_inv          ! inverse of seconds per hour

  REAL (KIND=wp)                   ::           &
    missval = -0.1_wp ! missing value of radar data (precip rates)


  CHARACTER (LEN=30)    ::  yroutine    ! name of the subroutine
  CHARACTER (LEN=255)   ::  yerrmsg     ! text message for model abort

! Local arrays  
!--------------
  INTEGER (KIND=iintegers), ALLOCATABLE ::  &
    iexchange(:)         ! array for exchange between PE's (used by global_values)

  REAL  (KIND=wp),     ALLOCATABLE ::           &
    pr_mod(:,:)     ,& ! total model precipitation rate              (kg/m2*s)
    pr_obs(:,:)     ,& ! observed (radar) precipitation rate         (kg/m2*s)
    pr_mod_nofilt(:,:),& !
    pr_obs_nofilt(:,:),& !
    pr_ana(:,:)     ,& ! analyzed precipitation rate                 (kg/m2*s)
    pr_at_dx(:,:)   ,& ! precipitation rate at observation altitude
    pr_at_dx_sum(:,:),& ! precipitation sum at observation altitude
    dbz_at_dx(:,:)   ,& ! reflectivity at observation altitude
    z_radar(:,:,:)   ,& ! reflectivity
    wobs_space(:,:) ,& ! weights (spatial) for the precip obs          ( 1 )
    wobs_time(:,:)     ! weights (temporal) for the precip obs         ( 1 )

  REAL  (KIND=wp),     ALLOCATABLE ::           &
    tt_lheat_int(:,:) ,& ! vert. integrated tt_lheat
    tt_clim_int(:,:)  ,& ! integrated climatolgical heating profile
    qrsflux_int(:,:)  ,& ! vert. integrated qrsflux
    scale_diag(:,:)   ,& ! global distribution of scale_fac
    treat_diag(:,:)   ,& ! diagnose of treatment
    windcor_diag(:,:) ,& ! weight with respect to the mean wind
    ktop_diag(:,:)       ! diagnosis of ktop

  REAL (KIND=wp),     ALLOCATABLE ::       &
    blacklist_tot(:,:)        ! complete field of blacklist information


  INTEGER (KIND=iintegers), ALLOCATABLE ::  &
    i_treat(:)      ,& ! i indeces of grid points to be treated by lhn
    j_treat(:)         ! j indeces of grid points to be treated by lhn

  LOGICAL, ALLOCATABLE :: &
    scale_fac_index(:,:)

  LOGICAL  :: &
    ltlhn           ,& ! true if latent heat nudging is active
    ltlhnverif      ,& ! true if radar verification is active
    zlhigh_freq_data   ! true if radar data is available in high temporal
                       ! frequency

  LOGICAL,SAVE  :: &
    ltlhnbright        ! true if brightband can be detected within time step

  INTEGER (KIND=iintegers) ::  &
    hfwidth, hfw_m_nb, ie_hf, je_hf

  INTEGER (KIND=iintegers) ::  &
    step_of_hour       ! timestep of current forecast hour

  REAL (KIND=wp),     POINTER :: &
    qv     (:,:,:)=> NULL(),&  ! QV at tlev=nnew
    qc     (:,:,:)=> NULL(),&  ! QC at tlev=nnew
    qi     (:,:,:)=> NULL(),&  ! QI at tlev=nnew
    qi_nn  (:,:,:)=> NULL(),&  ! QI at tlev=nnow_nold
    qg_nn  (:,:,:)=> NULL(),&  ! QG at tlev=nnow_nold
    qr_nn  (:,:,:)=> NULL(),&  ! QR at tlev=nnow_nold
    qs_nn  (:,:,:)=> NULL(),&  ! QS at tlev=nnow_nold
    qc_nn  (:,:,:)=> NULL()    ! QC at tlev=nnow_nold
#ifdef TWOMOM_SB
  REAL (KIND=wp),     POINTER :: &
    qh_nn  (:,:,:)=> NULL() ,& ! QH at tlev=nnow_nold
    qni_nn  (:,:,:)=> NULL(),& ! QI at tlev=nnow_nold
    qng_nn  (:,:,:)=> NULL(),& ! QG at tlev=nnow_nold
    qnr_nn  (:,:,:)=> NULL(),& ! QR at tlev=nnow_nold
    qns_nn  (:,:,:)=> NULL(),& ! QS at tlev=nnow_nold
    qnc_nn  (:,:,:)=> NULL(),& ! QC at tlev=nnow_nold
    qnh_nn  (:,:,:)=> NULL()   ! QH at tlev=nnow_nold
#endif

!-------------------------------------------------------------------------------
! End of declarations.    Public and private subroutines :
!-------------------------------------------------------------------------------

!===============================================================================

CONTAINS

!===============================================================================
!+ Module procedure "organize_lhn" in "lheat_nudge" to organize the
!  different program steps/subroutines for the latent heat nudging
!-------------------------------------------------------------------------------

SUBROUTINE organize_lhn

!-------------------------------------------------------------------------------
!
! Description:
! This subroutine is the main routine of LHN and is called by lmorg. It contains
! the organization of the LHN:
!
! 1. Read and prepare radar precipitation data
! 2. Get total model latent heating profiles: add convective latent heating 
!    contributions to large scale latent heating (in terms of temperature 
!    tendency as K/s)
! 3. Determine total model precipitation rate, analyze precipitation, i.e. 
!    combine model and observation values
! 4. Determine the latent heat nudging temperature increment by
!    scaling the model profiles; do profile search if requested
! 5. Apply the latent heat nudging temperature increment
! 6. Apply a corresponding humidity increment to maintain or produce
!    (nearly) saturation ("humidity enhancement") or reduce qv in case
!    of negative temperature increments (leave rel. hum. unchanged)
! 7. Output of values and profiles at diagnostic grid points
!    each timestep
!
!-------------------------------------------------------------------------------

! Local parameters, scalars, arrays:
!-------------------------------------------------------------------------------
! Local scalars:
!---------------

  LOGICAL, SAVE                    ::           &
    lfirst=.TRUE.        ! execute some parts only at first call of lhn routine
 
  INTEGER (KIND=iintegers)         ::           &
    i,j,k,             & ! loop indices
    izdebug              ! for debug output

  INTEGER (KIND=iintegers)         ::           &
    kqrs, &              ! upper layer with qrs_flux > 0.0
    k_dx

  REAL (KIND=wp)                   ::           &
    zdt_1                ! inverse of the timestep for physics ( = 1/dt )

  REAL (KIND=wp)                   ::           &
    tnow                 ! current model time (min; relative to start of frct)

  INTEGER (KIND=iintegers)         ::           &
    nstat             ,& ! IOSTAT of nulhn
    ntdeallow            ! last time in LHN

  INTEGER (KIND=iintegers), SAVE   ::           &
    steps_per_hour       ! timestep per forecast hour

! Local arrays:
!--------------

  REAL  (KIND=wp)           :: &
    prtot_gsp(ie,je)   ,& ! local array to compute a mask used in WHERE function
    prtot_con(ie,je)      ! local array to compute a mask used in WHERE function

  REAL (KIND=wp)     ::       &
    zprmod       (ie,je)  ,&
    zprmod_qrs   (ie,je)  ,&
    zprmod_ref   (ie,je)  ,&
    zprmodatdx   (ie,je)  ,&
    zprrad       (ie,je)  ,&
    zprmod_ref_f (ie,je)  ,&
    zprrad_f     (ie,je)

  REAL (KIND=wp)             ::       &
    vcoordsum,qrsgmax,qrsgthres,      &
    ztlhn   (ie,je,ke)   ! temporary variable for LHN

  INTEGER (KIND=iintegers) , SAVE :: &
    k_dx_min,k_dx_max

!- End of header
!===============================================================================
!-------------------------------------------------------------------------------
! Begin of subroutine
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Section 0 : Preliminaries : 
!             Check if lhn should be executed at the current timestep
!             Allocate space for fields; determine dt
!-------------------------------------------------------------------------------
  yroutine='organize_lhn'

  ! Initialize debug output variable
  IF (ldebug_ass) THEN
    IF (lprintdeb_all) THEN
      izdebug = idbg_level
    ELSE
      IF (my_cart_id == 0) THEN
        izdebug = idbg_level
      ELSE
        izdebug = 0
      ENDIF
    ENDIF
  ELSE
    izdebug = 0
  ENDIF

  ltlhn      = (llhn) .AND. (ntstep >= nlhn_start) .AND. (ntstep <= nlhn_end)
  ltlhnverif = (llhnverif) .AND. (ntstep >= nlhnverif_start) .AND. (ntstep <= nlhnverif_end)

 IF ((ltlhn) .OR. (ltlhnverif)) THEN
! allocate space for observational fields (long term storage)
  istat = 0
  IF (lfirst) THEN
     ALLOCATE (blacklist_tot(ie_tot,je_tot) , STAT=izlocstat );                                istat = istat + izlocstat

     ALLOCATE (obs         (ie,je,6)     , STAT=izlocstat);                                    istat = istat + izlocstat
     ALLOCATE (spqual      (ie,je,6)     , STAT=izlocstat);                                    istat = istat + izlocstat
     ALLOCATE (blacklist   (ie,je)       , STAT=izlocstat);  blacklist   (:,:)   = 0.0_wp; istat = istat + izlocstat
     ALLOCATE (dxheight    (ie,je,nradar), STAT=izlocstat);  dxheight    (:,:,:) = 0.0_wp; istat = istat + izlocstat
     ALLOCATE (tminc_lhn   (ie,je,ke)    , STAT=izlocstat);  tminc_lhn   (:,:,:) = 0.0_wp; istat = istat + izlocstat
     ALLOCATE (ttm_lheat   (ie,je,ke)    , STAT=izlocstat);  ttm_lheat   (:,:,:) = 0.0_wp; istat = istat + izlocstat
     ALLOCATE (pr_obs_sum  (ie,je)       , STAT=izlocstat);  pr_obs_sum  (:,:)   = 0.0_wp; istat = istat + izlocstat
     ALLOCATE (pr_mod_sum  (ie,je)       , STAT=izlocstat);  pr_mod_sum  (:,:)   = 0.0_wp; istat = istat + izlocstat
     ALLOCATE (pr_ref_sum  (ie,je)       , STAT=izlocstat);  pr_ref_sum  (:,:)   = 0.0_wp; istat = istat + izlocstat
     ALLOCATE (pr_at_dx_sum(ie,je)       , STAT=izlocstat);  pr_at_dx_sum(:,:)   = 0.0_wp; istat = istat + izlocstat
     ALLOCATE (brightband  (ie,je)       , STAT=izlocstat);  brightband  (:,:)   = 0.0_wp; istat = istat + izlocstat

! opening of files used in the nudging (and not only in the obs. processing)
! --------------------------------------------------------------------------

     IF (my_cart_id == 0) THEN
       OPEN (nulhn,FILE=yulhn,FORM='FORMATTED',IOSTAT=nstat)
       IF (nstat /= 0) THEN
          yerrmsg = 'OPENING OF FILE yulhn FAILED'
          CALL model_abort (my_cart_id, 6401, yerrmsg, yroutine)
       ENDIF
       REWIND nulhn
     ENDIF
 
     steps_per_hour = INT(sec_per_hr/dt)
     step_of_hour = 0_iintegers
  ENDIF   ! lfirst

! allocate space for the fields needed by several of the subroutines called
! (space deallocated again after termination of subroutine organize_lhn)

  ALLOCATE (tinc_lhn     (ie,je,ke), STAT=izlocstat);  tinc_lhn     (:,:,:) = 0.0_wp; istat = istat + izlocstat
  ALLOCATE (pr_mod       (ie,je)   , STAT=izlocstat);  pr_mod       (:,:)   = 0.0_wp; istat = istat + izlocstat
  ALLOCATE (pr_obs       (ie,je)   , STAT=izlocstat);  pr_obs       (:,:)   = 0.0_wp; istat = istat + izlocstat
  ALLOCATE (pr_mod_nofilt(ie,je)   , STAT=izlocstat);  pr_mod_nofilt(:,:)   = 0.0_wp; istat = istat + izlocstat
  ALLOCATE (pr_obs_nofilt(ie,je)   , STAT=izlocstat);  pr_obs_nofilt(:,:)   = 0.0_wp; istat = istat + izlocstat
  ALLOCATE (pr_ana       (ie,je)   , STAT=izlocstat);  pr_ana       (:,:)   = 0.0_wp; istat = istat + izlocstat
  ALLOCATE (pr_at_dx     (ie,je)   , STAT=izlocstat);  pr_at_dx     (:,:)   = 0.0_wp; istat = istat + izlocstat
  ALLOCATE (dbz_at_dx    (ie,je)   , STAT=izlocstat);  dbz_at_dx    (:,:)   = 0.0_wp; istat = istat + izlocstat
  ALLOCATE (z_radar      (ie,je,ke), STAT=izlocstat);  z_radar      (:,:,:) = 0.0_wp; istat = istat + izlocstat
  ALLOCATE (ttm_cv       (ie,je,ke), STAT=izlocstat);  ttm_cv       (:,:,:) = 0.0_wp; istat = istat + izlocstat
  ALLOCATE (wobs_space   (ie,je)   , STAT=izlocstat);  wobs_space   (:,:)   = 0.0_wp; istat = istat + izlocstat
  ALLOCATE (wobs_time    (ie,je)   , STAT=izlocstat);  wobs_time    (:,:)   = 0.0_wp; istat = istat + izlocstat
  ALLOCATE (i_treat      (ieje)    , STAT=izlocstat);  i_treat      (:)     = 0_iintegers;istat = istat + izlocstat
  ALLOCATE (j_treat      (ieje)    , STAT=izlocstat);  j_treat      (:)     = 0_iintegers;istat = istat + izlocstat

  ALLOCATE (tt_lheat_int (ie,je)   , STAT=izlocstat);  tt_lheat_int (:,:)   = 0.0_wp; istat = istat + izlocstat
  ALLOCATE (tt_clim_int  (ie,je)   , STAT=izlocstat);  tt_clim_int  (:,:)   = 0.0_wp; istat = istat + izlocstat
  ALLOCATE (windcor_diag (ie,je)   , STAT=izlocstat);  windcor_diag (:,:)   = 1.0_wp; istat = istat + izlocstat
  ALLOCATE (ktop_diag    (ie,je)   , STAT=izlocstat);  ktop_diag    (:,:)   = 1.0_wp; istat = istat + izlocstat
  ALLOCATE (scale_diag   (ie,je)   , STAT=izlocstat);  scale_diag   (:,:)   = 0.0_wp; istat = istat + izlocstat
  ALLOCATE (treat_diag   (ie,je)   , STAT=izlocstat);  treat_diag   (:,:)   =-1.0_wp; istat = istat + izlocstat

  IF (lhn_qrs) THEN
    ALLOCATE (qrsflux_int(ie,je)   , STAT=izlocstat);  qrsflux_int  (:,:)   = 0.0_wp; istat = istat + izlocstat
  ENDIF

  IF (num_compute > 1) THEN
    ALLOCATE ( iexchange(rlhn_search+2+11) , STAT=izlocstat )
    istat = istat + izlocstat
  ENDIF

  IF (istat /= 0) THEN
     yerrmsg =' ERROR  *** allocation of space for lhn - fields failed'
     CALL model_abort (my_cart_id,6402,yerrmsg,yroutine)
  ENDIF

! select timelevel and timestep for calculations
  IF (l2tls) THEN
    zdt = dt
    nnow_nold=nnew  ! nnow -> nnew !ks
  ELSE
    zdt = dt2
    nnow_nold=nnew
  ENDIF
  zdt_1 = 1.0_wp/zdt
  sec_per_hr_inv = 1.0_wp/sec_per_hr
  tnow  = ntstep*dt/60.0_wp

!ks>>
   tt_lheat_o (:,:,:) = 0.0_wp
   tinc_lhn_o (:,:,:) = 0.0_wp
   ttm_cv_o (:,:,:)   = 0.0_wp
!ks<<

!  IF(my_cart_id == 0) PRINT *,'LHN : relevant time step/time now : ',zdt,tnow

  ! retrieve the required microphysics tracers
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnew, ptr = qv)
  IF (izerror /= 0) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nnew, ptr = qc)
  IF (izerror /= 0) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF

  CALL trcr_get(izerror, idt_qi, ptr_tlev = nnew, ptr = qi)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qi, ptr_tlev = nnow_nold, ptr = qi_nn)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qr, ptr_tlev = nnow_nold, ptr = qr_nn)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nnow_nold, ptr = qc_nn)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  IF (itype_gscp >= 1) THEN
    CALL trcr_get(izerror, idt_qs, ptr_tlev = nnow_nold, ptr = qs_nn)
    IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
      yerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
    ENDIF
  ENDIF
  IF (itype_gscp >= 4) THEN
    CALL trcr_get(izerror, idt_qg, ptr_tlev = nnow_nold, ptr = qg_nn)
    IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
      yerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
    ENDIF
  END IF
#ifdef TWOMOM_SB
  IF (itype_gscp >= 100) THEN
    CALL trcr_get(izerror, idt_qh, ptr_tlev = nnow_nold, ptr = qh_nn)
    IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
      yerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
    ENDIF
    CALL trcr_get(izerror, idt_qnc, ptr_tlev = nnow_nold, ptr = qnc_nn)
    IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
      yerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
    ENDIF
    CALL trcr_get(izerror, idt_qnr, ptr_tlev = nnow_nold, ptr = qnr_nn)
    IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
      yerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
    ENDIF
    CALL trcr_get(izerror, idt_qni, ptr_tlev = nnow_nold, ptr = qni_nn)
    IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
      yerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
    ENDIF
    CALL trcr_get(izerror, idt_qns, ptr_tlev = nnow_nold, ptr = qns_nn)
    IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
      yerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
    ENDIF
    CALL trcr_get(izerror, idt_qng, ptr_tlev = nnow_nold, ptr = qng_nn)
    IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
      yerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
    ENDIF
    CALL trcr_get(izerror, idt_qnh, ptr_tlev = nnow_nold, ptr = qnh_nn)
    IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
      yerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
    ENDIF
  END IF
#endif

! ************************************************************************
! WARNING: The tracers QI, QG, QR, QS do not necessarily exist (depending
!   on the settings for the microphysics) and the code below will fail
!   very nastily in this case. This should be checked and handled grace-
!   fully in the near future. Specifically the call to radar_lm_ray()
!   should be modified or internally check for NULL pointers. QI is also
!   used at another place.
! ************************************************************************

  IF (ltime) THEN
     CALL get_timings (i_lhn_computations, ntstep, dt, izerror)
  ENDIF

!-------------------------------------------------------------------------------
! Section 1 : Read and prepare radar precipitation data 
!             - Determine if new data have to be read
!             - Read new data, project the observation data onto the model grid
!             - Distribute gridded observations to the PE's (if parallel)
!             - Determine spatial and temporal weights for observations 
!               (includes interpolation in time between consecutive obs)
!-------------------------------------------------------------------------------

     CALL lhn_obs_prep

     IF (ltime) THEN
        CALL get_timings (i_lhn_obs_prep, ntstep, dt, izerror)
     ENDIF

!-------------------------------------------------------------------------------
! Section 2 : Determine total model precipitation rate
!             Analyze precipitation, i.e. combine model and observation values
!-------------------------------------------------------------------------------

! total model precipitation: scale rain + scale snow + conv. rain + conv. snow
   prtot_gsp = prr_gsp + prs_gsp

   IF (itype_gscp >= 4 ) &
      prtot_gsp = prtot_gsp + prg_gsp

#ifdef TWOMOM_SB
   IF (itype_gscp >= 2000 ) &
      prtot_gsp = prtot_gsp + prh_gsp
#endif

   IF (lconv) THEN
      prtot_con = prr_con + prs_con
   ELSE
      prtot_con = 0.0_wp
   ENDIF
   pr_mod = prtot_gsp + prtot_con

   zprmod = pr_mod

 ENDIF
! ------------------------------------------------------------------------------
! Section 3: calculate skill scores
! ------------------------------------------------------------------------------
! IF (ltlhnverif) THEN

!     CALL gather_field (pr_mod, ie, je, zprmod_tot, ie_tot, je_tot, 0, izerror)
!     CALL gather_field (pr_obs, ie, je, zprrad_tot, ie_tot, je_tot, 0, izerror)
!     CALL gather_field (wobs_space, ie, je, wobs_space_tot, ie_tot, je_tot, 0, izerror)

!     CALL lhn_skill_scores

! ENDIF

! ------------------------------------------------------------------------------
! Section 4: get reference precipition for comparison of radar and model
! ------------------------------------------------------------------------------

 IF ((ltlhn) .OR. (ltlhnverif)) THEN 

!ks>> a new idea: take the vertikal integral of the precipitation flux as reference.
!   It is computed in src_gscp.hydci_pp or src_gscp.hydci_pp_gr

   IF (lhn_qrs) THEN
     DO j=jstart,jend
      DO i=istart,iend
       qrsflux_int(i,j) = 0.0_wp
       qrsgmax=MAXVAL(qrsflux(i,j,1:ke))
       qrsgthres=max(thres_lhn,rqrsgmax*qrsgmax)
       vcoordsum=0.0_wp
       kqrs=ke+1_iintegers
       DO k=1,ke
          IF ( qrsflux(i,j,k) >= qrsgthres ) then
             kqrs=k
             EXIT
          ENDIF
       ENDDO
       DO k=kqrs,ke
           qrsflux_int(i,j) = qrsflux_int(i,j) + qrsflux(i,j,k)  &
                            * (hhl(i,j,k)-hhl(i,j,k+1))
           vcoordsum=vcoordsum+(hhl(i,j,k)-hhl(i,j,k+1))
       ENDDO
       IF (vcoordsum /= 0.0_wp) qrsflux_int(i,j) = qrsflux_int(i,j) / vcoordsum
       ttm_cv(i,j,ke-22) = vcoordsum            ! ive: 26
      ENDDO
     ENDDO
     pr_mod = qrsflux_int
     zprmod_qrs = pr_mod
   ENDIF

   IF (ltime) THEN
      CALL get_timings (i_lhn_computations, ntstep, dt, izerror)
   ENDIF

!ks>> horizontal filtering of precipitation:
   pr_obs_nofilt(:,:) = pr_obs(:,:)
   pr_mod_nofilt(:,:) = pr_mod(:,:)
   IF (lhn_relax) THEN
      CALL hor_filt(pr_obs,nlhn_relax,1,1,.TRUE.)
      CALL hor_filt(pr_mod,nlhn_relax,1,1,.TRUE.)
      IF (ltime) THEN
          CALL get_timings (i_lhn_relax, ntstep, dt, izerror)
      ENDIF
   ENDIF
  
   IF (lhn_height) THEN
      IF (lfirst) THEN
         k_dx_min=ke+2
         k_dx_max=0
         DO j=1,je
            DO i=1,ie
               DO k_dx=1,nradar
                  IF ( dxheight(i,j,k_dx) > 0._wp ) THEN
                     DO k=1,ke
                        IF (dxheight(i,j,k_dx) < hhl(i,j,k) .AND. dxheight(i,j,k_dx) >= hhl (i,j,k+1)) THEN
                            k_dx_min=min(k_dx_min,k)
                            k_dx_max=max(k_dx_max,k)
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      IF (k_dx_max >= k_dx_min) THEN
         pr_at_dx(istart:iend,jstart:jend) =                        &
                    MAXVAL(qrsflux(istart:iend,jstart:jend,k_dx_min:k_dx_max),3)
      ENDIF
      IF (step_of_hour == steps_per_hour - 1 .OR. lfirst) THEN
        IF (itype_gscp == 3) THEN
          CALL radar_lm_ray (ie,je,ke, pi, rho_w, rho_ice, K_w, K_ice, t0_melt,    &
              klv850, my_cart_id, itype_gscp, izdebug, t(:,:,:,nnow_nold),         &
              qc_nn*rho, qr_nn*rho, qi_nn*rho, qs_nn*rho, z_radar = z_radar)
        ELSEIF (itype_gscp == 4) THEN
          CALL radar_lm_ray (ie,je,ke, pi, rho_w, rho_ice, K_w, K_ice,  t0_melt,   &
              klv850, my_cart_id, itype_gscp, izdebug, t(:,:,:,nnow_nold),         &
              qc_nn*rho, qr_nn*rho, qi_nn*rho, qs_nn*rho, q_grau=qg_nn*rho,        &
              z_radar = z_radar)
#ifdef TWOMOM_SB
        ELSEIF (itype_gscp >= 100) THEN
          CALL radar_sb_ray (ie,je,ke, pi, klv850, my_cart_id, t(:,:,:,nnow_nold), &
               qc_nn*rho, qr_nn*rho, qi_nn*rho, qs_nn*rho, qg_nn*rho, qh_nn*rho,   &
               qnc_nn*rho, qnr_nn*rho, qni_nn*rho, qns_nn*rho, qng_nn*rho,         &
               qnh_nn*rho, z_radar = z_radar )
#endif
        ELSE
          yerrmsg(:) = ' '
          yerrmsg =' ERROR  *** specified itype_gscp not implemented for radar_lm_ray'
          CALL model_abort (my_cart_id,6427,yerrmsg,yroutine)
        END IF
        IF (k_dx_max >= k_dx_min) &
          dbz_at_dx(istart:iend,jstart:jend)=MAXVAL(z_radar(istart:iend,jstart:jend,k_dx_min:k_dx_max),3)
      ENDIF
   ENDIF

 ENDIF

 IF (ltlhnverif) THEN
   zprmod_ref  (:,:) = pr_mod_nofilt(:,:)
   zprrad      (:,:) = pr_obs_nofilt(:,:)
   zprmod_ref_f(:,:) = pr_mod       (:,:)
   zprrad_f    (:,:) = pr_obs       (:,:)
   zprmodatdx  (istart:iend,jstart:jend) = pr_at_dx (istart:iend,jstart:jend)

 !USUS:  die Felder pr_mod*, pr_obs* werden nur fuer ltlhnverif angelegt:    IF (ltlhnverif) THEN
   !ks>> write diagnostic output of pr_rad, pr_mod and pr_mod_ref
   CALL lhn_verification ('SW',zprmod,zprmod_ref,zprrad,zprmodatdx,zprmod_ref_f,zprrad_f)
 ELSE
   zprmod_ref  (:,:) = 0.0_wp
   zprrad      (:,:) = 0.0_wp
   zprmod_ref_f(:,:) = 0.0_wp
   zprrad_f    (:,:) = 0.0_wp
   zprmodatdx  (:,:) = 0.0_wp
 ENDIF


! combine model and observed precipitation to analyzed precipitation field
! determine points (and number ntreat) for which LHN-modifs are necessary
! ------------------------------------------------------------------------------
! Section 5: get analized precipitation and all points which are treated within LHN
! ------------------------------------------------------------------------------

 IF (ltlhn) THEN

     CALL lhn_pr_ana(ntreat)
  
!-------------------------------------------------------------------------------
! Section 6: Get total model latent heating profiles
!            (in terms of temperature tendency as K/s)
!-------------------------------------------------------------------------------
! assure well defined large scale heating rates if current timestep is the first
! model time step (i.e. start of model integration), only values for nnew 
   IF (lfirst) THEN
      IF (.NOT.l2tls) tt_lheat(:,:,:,nold) = tt_lheat(:,:,:,nnew)
      tt_lheat(:,:,:,nnow) = tt_lheat(:,:,:,nnew)
      IF (my_cart_id == 0)  &
       PRINT *,'organize_lhn first called : set tt_lheat(.,.,.,old/now)=(.,new)'
   ENDIF

   DO  k=1,ke
!       WHERE ( prtot_gsp <= 0.0_wp )
!         tt_lheat(:,:,k,nnow_nold) = 0.0_wp
!       ELSEWHERE
         tt_lheat(:,:,k,nnow_nold) = tt_lheat(:,:,k,nnow_nold) * zdt_1
!       ENDWHERE
   ENDDO

   IF (ltime) THEN
      CALL get_timings (i_lhn_computations, ntstep, dt, izerror)
   ENDIF

! horizontal filtering of modelled latent heating rate

   IF (lhn_relax) THEN
      DO k=1,ke
         ttm_cv(:,:,ke-20)=ttm_cv(:,:,ke-20)+tt_lheat(:,:,k,nnow_nold)
      ENDDO
   
      CALL hor_filt(tt_lheat(:,:,:,nnow_nold),nlhn_relax,1,ke,.FALSE.)

      DO k=1,ke
         ttm_cv(:,:,ke-21)=ttm_cv(:,:,ke-21)+tt_lheat(:,:,k,nnow_nold)
      ENDDO
      IF (ltime) THEN
         CALL get_timings (i_lhn_relax, ntstep, dt, izerror)
      ENDIF
   ENDIF

!-------------------------------------------------------------------------------
! Section 7 : Determine the latent heat nudging temperature increment by
!             scaling the model profiles; do profile search if requested
!             (This step needs communication between nodes to exchange the
!             heating profiles that are found outside the PE's subdomain)
!-------------------------------------------------------------------------------

   ALLOCATE (scale_fac_index(ie,je),STAT=izlocstat)
   IF (izlocstat /= 0) THEN
     yerrmsg =' ERROR  *** allocation of space for lhn - fields failed'
     CALL model_abort (my_cart_id,6403,yerrmsg,yroutine)
   ENDIF
   scale_fac_index = .FALSE.

   CALL lhn_t_inc(ntreat)

!-------------------------------------------------------------------------------
! Section 8 : Apply the latent heat nudging temperature increment
!-------------------------------------------------------------------------------

   DO     k=1,ke
     DO   j=jstart,jend
       DO i=istart,iend
          IF (ltlhnbright .AND. NINT(brightband(i,j)) > 0_iintegers) tinc_lhn(i,j,k) = 0.0_wp
          tinc_lhn(i,j,k) = lhn_coef * tinc_lhn(i,j,k)
          t(i,j,k,nnow_nold) = t(i,j,k,nnow_nold) + tinc_lhn(i,j,k) * zdt
       ENDDO
     ENDDO
   ENDDO

   IF (ltime) THEN
      CALL get_timings (i_lhn_t_inc, ntstep, dt, izerror)
   ENDIF

!-------------------------------------------------------------------------------
! Section 9 : Apply a corresponding humidity increment to maintain or produce
!             (nearly) saturation ("humidity enhancement") or reduce qv in case
!             of negative temperature increments (leave rel. hum. unchanged)
!             or
!             Apply saturation adjustment to balance the thermodynamic fields
!             after the LHN
!-------------------------------------------------------------------------------

IF (itype_gscp < 100 .OR. l2mom_satads) THEN

  IF (lhn_hum_adj) THEN
   CALL lhn_q_inc(ntreat)
  ELSE
   CALL lhn_satad(ntreat)
  ENDIF

ENDIF

  IF (ltime) THEN
     CALL get_timings (i_lhn_q_inc, ntstep, dt, izerror)
  ENDIF

  DEALLOCATE (scale_fac_index, STAT=izlocstat)

! rescale of tt_lheat to unit of K, this should be done, because tt_lheat has a time dimension
! and could be used one time step later in the same way as now

  tt_lheat(:,:,:,nnow_nold) = tt_lheat(:,:,:,nnow_nold) * zdt

 ENDIF

 IF ((ltlhn) .OR. (ltlhnverif)) THEN
!-------------------------------------------------------------------------------
! Section 10 : Diagnostic procedure...
! a) get accumulated total model heating and heating increments
! b) integrate observed precipitation rates over one hour
! c) control output for some variables via ttm_cv
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Section 10a : get accumulated total model heating and heating increments
!-------------------------------------------------------------------------------

! for control output: summation of temperature increments due to lhn
   tminc_lhn(:,:,:) = tminc_lhn(:,:,:) + tinc_lhn(:,:,:) * zdt

! for control output: accumulated model heating : -> ttm_lheat
   ttm_lheat(:,:,:) = ttm_lheat(:,:,:) + tt_lheat(:,:,:,nnow_nold)
  

!-------------------------------------------------------------------------------
! Section 10b : integrate observed precipitation rates over one hour
!-------------------------------------------------------------------------------
! for verification integrate observed precipitation rates over one hour:

   WHERE (zprrad > 0.0_wp)     pr_obs_sum(:,:) = pr_obs_sum(:,:)       &
                                                   + zprrad(:,:) * dt
   WHERE (zprmod > 0.0_wp)     pr_mod_sum(:,:) = pr_mod_sum(:,:)       &
                                                   + zprmod(:,:) * dt
   WHERE (zprmod_ref > 0.0_wp) pr_ref_sum(:,:) = pr_ref_sum(:,:)       &
                                                   + zprmod_ref(:,:) * dt
   WHERE (pr_at_dx > 0.0_wp)   pr_at_dx_sum(:,:) = pr_at_dx_sum(:,:)   &
                                                   + pr_at_dx(:,:) * dt

!-------------------------------------------------------------------------------
! Section 10c : control output for some variables via ttm_cv
!-------------------------------------------------------------------------------
! control output of pr_mod, pr_obs, pr_ana in ttm_cv (upper 3 levels)
! ouput as mm/h , since pr_mod is in kg/(s*m**2) and scaling is *3600.

   WHERE (pr_obs >= 0.0_wp)
      ttm_cv(:,:,ke) = pr_obs
   ELSEWHERE                                 ! ive: 1
      ttm_cv(:,:,ke) = -1.0_wp
   ENDWHERE
   ttm_cv(:,:,ke-1) = pr_mod                 ! ive: 2
   ttm_cv(:,:,ke-2) = pr_ana                 ! ive: 3
   ttm_cv(:,:,ke-3) = wobs_space             ! ive: 4
   ttm_cv(:,:,ke-4) = wobs_time              ! ive: 5
   WHERE (pr_obs > thres_lhn .AND. pr_mod > thres_lhn)
      ttm_cv(:,:,ke-5) = 1.0_wp          ! ive: 6
   ENDWHERE
   WHERE (pr_obs > thres_lhn .AND. pr_mod <= thres_lhn)
      ttm_cv(:,:,ke-5) = 2.0_wp          ! ive: 6
   ENDWHERE
   WHERE (pr_obs <= thres_lhn .AND. pr_mod > thres_lhn)
      ttm_cv(:,:,ke-5) = 3.0_wp          ! ive: 6
   ENDWHERE
   WHERE (pr_obs > pr_mod)
      ttm_cv(:,:,ke-6) = 1.0_wp          ! ive: 7
   ELSEWHERE
      ttm_cv(:,:,ke-6) = 0.0_wp          ! ive: 7
   ENDWHERE
   WHERE (pr_obs < pr_mod)
      ttm_cv(:,:,ke-6) = -1.0_wp         ! ive: 7
   ELSEWHERE
      ttm_cv(:,:,ke-6) = 0.0_wp          ! ive: 7
   ENDWHERE
   WHERE (MINVAL(tinc_lhn,3) == 0._wp .AND. MAXVAL(tinc_lhn,3) == 0._wp )
      ttm_cv(:,:,ke-7) = 0.0_wp          ! ive: 8
   ELSEWHERE
      ttm_cv(:,:,ke-7) = 1.0_wp          ! ive: 8
   ENDWHERE
   ttm_cv(:,:,ke-8) = zprmod                 ! ive: 9
   ttm_cv(:,:,ke-9) = zprmod_ref             ! ive: 10
   ttm_cv(:,:,ke-10) = zprrad                ! ive: 11
   ttm_cv(:,:,ke-11) = pr_obs_sum            ! ive: 12
   ttm_cv(:,:,ke-12) = pr_mod_sum            ! ive: 13
   ttm_cv(:,:,ke-13) = pr_ref_sum            ! ive: 14
   ttm_cv(:,:,ke-14) = zprmod_ref_f          ! ive: 15
   ttm_cv(:,:,ke-15) = zprrad_f              ! ive: 16
   ttm_cv(:,:,ke-16) = tt_lheat_int          ! ive: 17
   IF (lhn_qrs) ttm_cv(:,:,ke-17) = qrsflux_int ! ive: 18
   ttm_cv(:,:,ke-18) = scale_diag(:,:)       ! ive: 19
   ttm_cv(:,:,ke-19) = treat_diag            ! ive: 20
!   ttm_cv(:,:,ke-20) = tt_lheat_nofilt(:,:,:,:) ! ive: 21 used for unfiltered tt_lheat (see above)
!   ttm_cv(:,:,ke-21) = tt_lheat_filt(:,:,:,:)   ! ive: 22 used for filtered tt_lheat (see above)
!   ttm_cv(:,:,ke-22) = vcoordsum            ! ive: 23 already defined (see above)
   IF (lhn_qrs) ttm_cv(:,:,ke-25) = zprmod_qrs  ! ive: 26
   ttm_cv(:,:,ke-26) = windcor_diag          ! ive: 27
   ttm_cv(:,:,ke-27) = tt_clim_int           ! ive: 28
   ttm_cv(:,:,ke-28) = ktop_diag             ! ive: 29
   ttm_cv(:,:,ke-29) = brightband            ! ive: 30
   ttm_cv(:,:,ke-30) = pr_at_dx              ! ive: 31
   ttm_cv(:,:,ke-31) = pr_at_dx_sum          ! ive: 32
   ttm_cv(:,:,ke-32) = dbz_at_dx             ! ive: 33
!     ttm_cv(:,:,ke-33) = height_interval(:,:) ---> defined in SR detect_bright_band
!     ttm_cv(:,:,ke-34) = minheight(:,:)       ---> defined in SR detect_bright_band
!     ttm_cv(:,:,ke-35) = maxheight(:,:)       ---> defined in SR detect_bright_band

   tt_lheat_o (:,:,:) = tt_lheat (:,:,:,nnow_nold)
   tinc_lhn_o (:,:,:) = tinc_lhn (:,:,:)
   ttm_cv_o (:,:,:)   = ttm_cv(:,:,:)

!   steps_per_hour = INT(sec_per_hr/dt)
   step_of_hour   = step_of_hour + 1_iintegers
   IF (step_of_hour == steps_per_hour) THEN
      IF (ltlhnverif) THEN
         CALL lhn_verification ('HR',pr_mod_sum,pr_ref_sum,pr_obs_sum,pr_at_dx_sum)
      ENDIF
      step_of_hour     = 0_iintegers
      ttm_lheat(:,:,:) = 0.0_wp
      pr_obs_sum(:,:)  = 0.0_wp
      pr_mod_sum(:,:)  = 0.0_wp
      pr_ref_sum(:,:)  = 0.0_wp
      pr_at_dx_sum(:,:)  = 0.0_wp
   ENDIF

!-------------------------------------------------------------------------------
! Section 11 : Deallocate fields needed only during the lhn - step
!-------------------------------------------------------------------------------

   DEALLOCATE ( tinc_lhn, STAT=izlocstat )
   DEALLOCATE ( pr_mod, STAT=izlocstat )
   DEALLOCATE ( pr_obs, STAT=izlocstat )
   DEALLOCATE ( pr_mod_nofilt, STAT=izlocstat )
   DEALLOCATE ( pr_obs_nofilt, STAT=izlocstat )
   DEALLOCATE ( pr_ana, STAT=izlocstat )
   DEALLOCATE ( pr_at_dx, STAT=izlocstat )
   DEALLOCATE ( dbz_at_dx, STAT=izlocstat )
   DEALLOCATE ( z_radar, STAT=izlocstat )
   DEALLOCATE ( ttm_cv, STAT=izlocstat )
   DEALLOCATE ( wobs_space, STAT=izlocstat )
   DEALLOCATE ( wobs_time, STAT=izlocstat )
   DEALLOCATE ( i_treat, STAT=izlocstat )
   DEALLOCATE ( j_treat, STAT=izlocstat )
   IF (lhn_qrs)  DEALLOCATE ( qrsflux_int, STAT=izlocstat )
   DEALLOCATE ( tt_lheat_int, STAT=izlocstat )
   DEALLOCATE ( tt_clim_int, STAT=izlocstat )
   DEALLOCATE ( scale_diag , STAT=izlocstat )
   DEALLOCATE ( treat_diag , STAT=izlocstat )
   DEALLOCATE ( windcor_diag , STAT=izlocstat )
   DEALLOCATE ( ktop_diag , STAT=izlocstat )

   IF(num_compute > 1) DEALLOCATE ( iexchange, STAT=izlocstat )

   ! Added update of tt_lheat and re-initialization to 0 for timelevel nnew 
   ! (was in src_relaxation before, which is now called before the assimilation)
   DO k = 1, ke
      IF ( .NOT.l2tls ) tt_lheat(:,:,k,nold) = tt_lheat(:,:,k,nnow)
      tt_lheat(:,:,k,nnow) = tt_lheat(:,:,k,nnew)
   ENDDO
   ztlhn = 0.0_wp
   CALL get_gs_lheating ('new',1,ke,ztlhn) ! set tt_lheat(:,:,:,nnew) = 0.0

   ntdeallow = MAX(nlhn_end,nlhnverif_end)
   IF(ntstep == ntdeallow) THEN
    DEALLOCATE ( obs , STAT=izlocstat )
    DEALLOCATE ( spqual, STAT=izlocstat )
    DEALLOCATE ( blacklist, STAT=izlocstat )
    DEALLOCATE ( dxheight, STAT=izlocstat )
    DEALLOCATE ( tminc_lhn, STAT=izlocstat )
    DEALLOCATE ( ttm_lheat, STAT=izlocstat )
    DEALLOCATE ( pr_obs_sum, STAT=izlocstat )
    DEALLOCATE ( pr_mod_sum, STAT=izlocstat )
    DEALLOCATE ( pr_ref_sum, STAT=izlocstat )
    DEALLOCATE ( pr_at_dx_sum, STAT=izlocstat )
    DEALLOCATE ( brightband, STAT=izlocstat )
   
    IF (my_cart_id == 0) CLOSE (nulhn)
   ENDIF

! set switch for first call to false

   lfirst = .FALSE.

  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure organize_lhn
!-------------------------------------------------------------------------------

END SUBROUTINE organize_lhn

!===============================================================================
!+ Module procedure in "lheat_nudge" preparing radar precip data for lhn
!-------------------------------------------------------------------------------

SUBROUTINE lhn_obs_prep

!-------------------------------------------------------------------------------
!
! Description:
!   This procedure of the module "lheat_nudge" is called from organize_lhn.
!   It reads in observations (rain rates in mm/h) from the radar data input
!   file and prepares them for further use:
!   The observations are interpolation in time
!   between consecutive observation times is done and spatial and temporal
!   weighting factors are assigned to the resulting observations at each
!   grid point.
!
!   Namelist parameters used: 
!            rad_wobs_lhn, lhn_black, blacklist_file, radar_in
!   Input arrays : none. (Use of general information on model grid)
!   Output arrays: pr_obs,wobs_space,wobs_time
!
! Method:
!   1a: Determine if new observations have to be read at current timestep
!   1b: Read the observational data from binary file 
!   1c: Distribute the field with new radar rain rates + spatial quality information
!       to different PE's
!   1d: Update time informaton for observation fields and send to PE's
!    2: Interpolate observation data in time and assign temporal weight:
!       (depending on the local variables dtobsint and dtinflmax)
!    3: Determine the spatial normalized weight (depending on namelist
!       parameters rad_wobs_lhn and local variables radmax, wobs_prmin)
!
! Note : Only sections 1c and 1d contain communication between nodes.
!
! Input files:
!    Radar data input file radar blacklist file read in subroutine lhn_obs_read. See lhn_obs_read
!    for further information.
!
!-------------------------------------------------------------------------------

! Subroutine / Function arguments: None
!-------------------------------------------------------------------------------

! Local parameters, scalars, arrays :
!-------------------------------------------------------------------------------
! Local scalars:
!---------------
! 1. Variables for organizing the code
  LOGICAL, SAVE                        ::       &
    lfirst=.TRUE.  ,&   ! indicator for steps to be executed at first call only
    lhn_obs_read_status ! indicator for steps to be executed only for data input

  LOGICAL                              ::       &
    lhm1           ,&   ! first observation 1 hour older than model start
    lobsread,lqualread

! 2. Variables for input of observational data

  INTEGER (KIND=iintegers), SAVE             ::       &
    nulhnrad        ,&  ! unit number for radardata GRIB file
    nulhnblackdx    ,&  ! unit number for blacklist file (DX)
    nulhnh              ! unit number for height file (DX)

  INTEGER (KIND=iintegers)                   ::       &
    nrec_max      ,& ! number of maximal expected data records to be read
    nobs          ,& ! number of currently read data records
    ierr          ,& ! error flag
    istat,izlocstat,&!
    obs_time_s, n

  INTEGER (KIND=iintegers), SAVE             ::       &
    next_obs_time ,& !
    nrec          ,& ! number of currently read records
    field_num     ,& ! number of data record in memory field
    obs_time(6)      !

! variables for time difference subroutine 
  INTEGER (KIND=iintegers),SAVE  ::       &
    cc_rad    ,& ! century   /
    yy_rad    ,& ! year     / current
    mm_rad    ,& ! month   /  date and time
    dd_rad    ,& ! day    <   of the
    hh_rad    ,& ! hour    \  radar
    min_rad   ,& ! minute   \ data
    cc_black  ,&
    yy_black  ,&
    mm_black  ,&
    dd_black  ,&
    hh_black  ,&
    min_black ,&
    min0_lhn  ,& ! minute of first LHN step after model start
    nactday   ,& ! day of the year
    day_of_month(12) ! number of days of month

  REAL (KIND=wp),     SAVE             ::       &
    pr_time_limit ,& !
    hh0_lhn       ,& ! hour    \  of first LHN step 
    delta_t_file     ! nominal time interval of consecutive radar data files

  CHARACTER (LEN=14), SAVE    ::  actdate   ! current date
  CHARACTER (LEN=12), SAVE    ::  obsdate   ! current date
  CHARACTER (LEN=28)    ::  dum1

! 3. Variables used for calculating the temporal weights of the observations

  INTEGER (KIND=iintegers), SAVE             ::       &
     iread          ,& ! number of data fields in time cache
     center_time    ,& ! time of actual radar observation
     next_time_1    ,& ! time of next observation
     next_time_2    ,& ! time of next observation
     next_time_3    ,& ! time of next observation
     prev_time_1    ,& ! time of previos observation
     prev_time_2    ,& ! time of previos observation
     weight_index_0 ,& !
     weight_index_p1,& !
     weight_index_p2,& !
     weight_index_p3,& !
     weight_index_m1,& !
     weight_index_m2   !


  INTEGER (KIND=iintegers)                   ::       &
    i,j                ! loop counter

  INTEGER (KIND=iintegers)                   ::       &
    inoobs         ,& ! loop counter for noobs_dates
    delta_t0          ! is used as integer below!

  REAL (KIND=wp)                       ::       &
    tnow           ,& ! current model time (min; relative to start of forecast)
    twlast         ,& ! weight for last obs
    twnext            ! weight for next obs

!    tnext             ! model time of next model step 
                      ! (min, relative to start of forecast)

  INTEGER (KIND=iintegers)             ::       &
! concerning the time interpolation of DX data
    num1delta_t_obs ,& ! number of points with obs-dist 1 delta_t
    num2delta_t_obs ,& ! number of points with obs-dist 2 delta_t
    num3delta_t_obs ,& ! number of points with obs-dist 3 delta_t
    num4delta_t_obs ,& ! number of points with obs-dist 4 delta_t
    numnone         ,& ! number of points without obs
    numblack        ,& ! number of points which are blacklisted
! concerning the spatial interpolation of radar data
    numfull         ,& ! number of points with obs with full spatial weight
    numred          ,& ! number of points with reduced weight
    numzero         ,& ! number of points with obs zero spatial weight
    znmissval          ! number of missing values

  INTEGER (KIND=iintegers), SAVE             ::       &
    ntbright
! Local arrays:
!--------------
! 1. temporary arrays for observational input (all data held on one processor)

  REAL  (KIND=wp)                    ::       &
    obs_tot      (ie_tot,je_tot), & ! complete field of gridded observations 
    spqual_tot   (ie_tot,je_tot), & ! complete field of gridded spatial quality function
    dxheight_tot (ie_tot,je_tot,nradar) , & ! complete field of height information
    sumrad_tot   (ie_tot,je_tot), sumrad(ie,je)

  REAL  (KIND=wp),    SAVE,ALLOCATABLE      ::       &
    obs_tot_all  (:,:,:)            ! complete field of gridded observations

  INTEGER (KIND=iintegers),SAVE,ALLOCATABLE ::       &
    obs_var_read  (:),&
    obs_tab_read  (:)

  CHARACTER (LEN=12),SAVE,ALLOCATABLE :: obs_date_read (:)

! Local variables, arrays for communication of reading PE to all PE's :
!----------------------------------------------------------------------
  INTEGER (KIND=iintegers), PARAMETER ::   &
    ibuflen = 5          ! length of buffer 

  INTEGER (KIND=iintegers)            ::   &
    nzbytes              ! bytes per word for grib packing

!- End of header
!===============================================================================
!-------------------------------------------------------------------------------
! Begin of subroutine
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Section 0: Initialization of some variables
!            At first call: Get start date and time of this model run
!                           Open data file and read general header information
!            At each call : Allocate space for temporal fields
!
!-------------------------------------------------------------------------------
  yroutine='lhn_obs_prep'

  lhn_obs_read_status=.FALSE.  
  day_of_month = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  inoobs=1
  nrec_max = INT(60._wp/lhn_dt_obs)
  IF (lhn_spqual ) nrec_max = 2 * nrec_max
  lqualread=.true.


  ! define if radar data is high frequency data, i.e. if lhn_dt_obs is below
  ! or equal to a specific time interval. We define it to be 10 minutes.
  IF (lhn_dt_obs <= 10.0_wp) THEN
     zlhigh_freq_data = .TRUE.
  ELSE
     zlhigh_freq_data = .FALSE.
  ENDIF

  tnow  = REAL(ntstep*dt,wp)/60.0_wp

  IF (lfirst) THEN
     
     ! set defaults to be assigned to gridded arrays
     pr_time_limit= -0.01_wp
     delta_t_file = 60.0_wp
     ntbright=0
     sumrad=0.0_wp
     ALLOCATE (obs_var_read(nrec_max),obs_tab_read(nrec_max),obs_date_read(nrec_max))
     ALLOCATE (obs_tot_all(ie_tot,je_tot,nrec_max))
     obs_var_read = 0_iintegers
     obs_tab_read = 0_iintegers
     obs_date_read = "190001010100"
     obs_tot_all  = missval

     IF (zlhigh_freq_data) THEN
        iread = 6_iintegers
        delta_t0   = 10_iintegers
     ELSE
       iread = 2_iintegers
       delta_t0 = 0_iintegers
     ENDIF

     CALL get_utc_date (ntstep,ydate_ini,dt,itype_calendar,actdate,dum1,nactday,hh0_lhn)
     READ(actdate,'(6I2)') cc_rad,yy_rad,mm_rad,dd_rad,hh_rad,min0_lhn
     IF (my_cart_id == 0) THEN
! initialize variables  
       lds = ie_tot*je_tot
       lbm = lds / BIT_SIZE (1_intgribf) + 1_iintegers 
       nzbytes = 8        ! for grib packing: bytes per word (to be on the safe side)
       lfd = compute_grib_intbuffer_length(ie_tot, je_tot, nzbytes, iwlength)

       istat=0
       ALLOCATE ( iblock_rad(lfd), ibmap_rad(lbm), STAT=izlocstat )
       istat = istat + izlocstat
       ALLOCATE ( ds_rad(lds), dsup_rad(ndsup), STAT=izlocstat )
       istat = istat + izlocstat

       IF (istat /= 0) THEN
          yerrmsg =' ERROR  *** allocation of space for lhn - input fields failed'
          CALL model_abort (my_cart_id,6404,yerrmsg,yroutine)
       ENDIF

! read in blacklist
       cc_black=20_iintegers
       yy_black=04_iintegers
       mm_black=07_iintegers
       dd_black=19_iintegers
       hh_black=23_iintegers
       min_black=45_iintegers
       blacklist_tot(:,:)=0.0_wp
       IF (lhn_black) &
          CALL lhn_blacklist_dx_open_and_read (nulhnblackdx,blacklist_file,cc_black,&
                                               yy_black,mm_black,dd_black,hh_black,min_black,blacklist_tot,ierr)

       dxheight_tot(:,:,:)=0.0_wp
       IF (lhn_height) &
           CALL lhn_dx_height_open_and_read (nulhnh,height_file,dxheight_tot,ierr)

       ltlhnbright=(lhn_bright .AND. MAXVAL(dxheight_tot) > 0.0_wp)

       tnow  = REAL(ntstep*dt,wp)/60.0_wp

       ! determine time for the first radar image
       min0_lhn=NINT(REAL(min0_lhn,wp)/(lhn_dt_obs))*NINT(lhn_dt_obs)
       min_rad=min0_lhn-delta_t0
       obs_time(1)=NINT(tnow/lhn_dt_obs)*NINT(lhn_dt_obs)-delta_t0

       lhm1=.false.
       IF (min_rad < 0) THEN
          min_rad = min_rad + 60
          min_rad=NINT(REAL(min_rad, wp)/(lhn_dt_obs))*NINT(lhn_dt_obs)
          hh_rad  = hh_rad  - 1
          lhm1=.true.

          IF (hh_rad < 0) THEN
             hh_rad  = hh_rad + 24
             dd_rad  = dd_rad - 1

             IF (dd_rad == 0) THEN
                mm_rad = mm_rad - 1

                IF (mm_rad == 0) THEN
                   mm_rad = 12
                   dd_rad = day_of_month(mm_rad)
                   yy_rad = yy_rad - 1
                   IF (yy_rad < 0) THEN
                      yy_rad = yy_rad + 100
                      cc_rad = cc_rad - 1
                   ENDIF
                ELSE IF (mm_rad == 2) THEN
                   IF ((yy_rad == 0 .AND. MOD(cc_rad,4)==0) .OR. (MOD(yy_rad,4) == 0)) THEN
                        dd_rad = day_of_month(mm_rad) + 1
                   ELSE
                        dd_rad = day_of_month(mm_rad)
                   ENDIF
                ELSE
                   dd_rad = day_of_month(mm_rad)
                ENDIF
            ENDIF
         ENDIF
       ENDIF

       WRITE(obsdate,'(6i2.2)') cc_rad,yy_rad,mm_rad,dd_rad,hh_rad,min_rad

! Open first data file
       CALL lhn_obs_open (nulhnrad,obsdate,ierr)
       IF (ierr == 0 ) THEN
         CALL lhn_obs_read(nulhnrad,nrec_max,obs_tot_all,obs_date_read,obs_var_read,obs_tab_read,nrec,ierr)
         IF (ltlhnbright .AND. .NOT. lhm1) THEN
           nobs=nrec
           IF (lhn_spqual) nobs=NINT(REAL(nrec,wp)*0.5_wp)
           CALL lhn_sumrad(nobs,obs_tot_all,sumrad_tot,blacklist_tot,ierr)
           ntbright=1_iintegers
         ENDIF
         CALL close_file( nulhnrad, yform_read, 0, 0, 1, .FALSE., idbg_level, yerrmsg, ierr)
         IF (ierr /= 0) THEN
           WRITE(nulhn,*) 'Something wrong by closing radar file with ('//yform_read//')'
         ENDIF
       ENDIF
     ENDIF   ! my_cart_id == 0

     IF (num_compute > 1) &
         CALL distribute_values (ltlhnbright,1,0,imp_logical,icomm_cart,ierr)

     DO i=1,iread
        field_num=i
        IF (my_cart_id == 0) THEN
           IF (obsdate == noobs_date(inoobs)) THEN
               WRITE(nulhn, *)'no data are used for ',noobs_date(inoobs)
               obs_tot(:,:)=missval
               spqual_tot(:,:)=missval
               inoobs=inoobs+1
               WRITE(nulhn,*)'next noobs:',inoobs,noobs_date(inoobs)
           ELSE
             obs_tot(:,:)=missval
             spqual_tot(:,:)=missval
             lobsread=.false.
             IF (lhn_spqual) lqualread=.false.
             DO n=1,nrec
               IF ( obsdate == obs_date_read(n)) THEN
                 IF ( obs_var_read(n) == 61 .AND. obs_tab_read(n) == 2) THEN
                    obs_tot(:,:) = obs_tot_all(:,:,n)
                    WRITE(nulhn,*)'next noobs:',inoobs,noobs_date(inoobs)
                    lobsread=.true.
                 ELSE IF ( obs_var_read(n) == 162 .AND. obs_tab_read(n) == 250) THEN
                    spqual_tot(:,:) = obs_tot_all(:,:,n)
                    lqualread=.true.
                 ELSE
                    WRITE(nulhn, *) &
                    'WARNING: observation file does not include precipiation data (ee=61) nor quality data (ee=162) for ',obsdate
                    EXIT
!                    obs_tot(:,:)    = missval
!                    spqual_tot(:,:) = missval
                 ENDIF
!               ELSE
!                 WRITE(nulhn, *) &
!                 'WARNING: observation file does not include precipiation data (ee=61) nor quality data (ee=162) for ',obsdate
!                 obs_tot(:,:)    = missval
!                 spqual_tot(:,:) = missval
                 IF (lobsread .AND. lqualread) EXIT
               ENDIF
             ENDDO
           ENDIF

           IF (lhn_diag) THEN
             ! write availability and statistics of radar data
             znmissval = count(obs_tot/=missval)
             IF (znmissval == 0) THEN
                WRITE(nulhn, *) 'data available: 0 ', obsdate
                WRITE(nulhn, '("Stat. of obs    at ",a12," (min,max,mean): ",3a12)') obsdate,&
                     'undef','undef','undef'
             ELSE
                WRITE(nulhn, *) 'data available: 1 ', obsdate
                WRITE(nulhn, '("Stat. of obs    at ",a12," (min,max,mean): ",3f12.5)') obsdate,&
                     MINVAL(obs_tot,MASK=obs_tot/=missval),&
                     MAXVAL(obs_tot,MASK=obs_tot/=missval),&
                     SUM(obs_tot,MASK=obs_tot/=missval)/REAL(znmissval,wp)
                WRITE(nulhn,*) "Number of missval in obs: ",ie_tot*je_tot-znmissval
             END IF

             IF (lhn_spqual) THEN
                ! write statistics of spatial quality data
                znmissval = count(spqual_tot/=missval)
                IF (znmissval == 0) THEN
                   WRITE(nulhn, '("Stat. of spqual at ",a12," (min,max,mean): ",3a12)') obsdate,&
                        'undef','undef','undef'
                ELSE
                   WRITE(nulhn, '("Stat. of spqual at ",a12," (min,max,mean): ",3f12.5)') obsdate,&
                        MINVAL(spqual_tot,MASK=spqual_tot/=missval),&
                        MAXVAL(spqual_tot,MASK=spqual_tot/=missval),&
                        SUM(spqual_tot,MASK=spqual_tot/=missval)/REAL(znmissval,wp)
                   WRITE(nulhn,*) "Number of missval in spqual: ",ie_tot*je_tot-znmissval
                END IF
             END IF
           ENDIF

           IF (lhn_spqual) THEN

              ! replace missval of spqual with zero
              WHERE(spqual_tot==missval) spqual_tot = 0.0_wp

              ! Check spqual for values larger than one and smaller than zero
              znmissval = count(spqual_tot>1.0_wp)
              IF (znmissval > 0) THEN
                 WRITE(nulhn,*) "WARNING: values larger than one in spqual, correcting them"
                 where (spqual_tot > 1.0_wp) spqual_tot = 1.0_wp
              END IF
              znmissval = count(spqual_tot<0.0_wp)
              IF (znmissval > 0) THEN
                 WRITE(nulhn,*) "WARNING: values smaller than zero in spqual, correcting them"
                 where (spqual_tot < 0.0_wp) spqual_tot = 0.0_wp
              END IF

           END IF



           min_rad=min_rad+NINT(lhn_dt_obs)
           IF (min_rad >= 60) THEN

              min_rad=min_rad - 60
              min_rad=NINT(REAL(min_rad,wp)/(lhn_dt_obs))*NINT(lhn_dt_obs)
              hh_rad = hh_rad + 1
              IF (hh_rad >= 24) THEN
                 hh_rad = 0
                 dd_rad = dd_rad + 1
                 IF (dd_rad > day_of_month(mm_rad) ) THEN
                    dd_rad = 1
                    mm_rad = mm_rad + 1
                    IF (mm_rad > 12) THEN
                       mm_rad = 1
                       yy_rad = yy_rad + 1
                       IF (yy_rad > 99) THEN
                          yy_rad = yy_rad - 100
                          cc_rad = cc_rad + 1
                       ENDIF
                   ENDIF
                 ENDIF
              ENDIF
              WRITE(obsdate,'(6i2.2)') cc_rad,yy_rad,mm_rad,dd_rad,hh_rad,min_rad
              CALL lhn_obs_open (nulhnrad,obsdate,ierr)
              IF (ierr == 0 ) THEN
                CALL lhn_obs_read(nulhnrad,nrec_max,obs_tot_all,obs_date_read,obs_var_read,obs_tab_read,nrec,ierr)
                IF (ltlhnbright) THEN
                   nobs=nrec
                   IF (lhn_spqual) nobs=NINT(REAL(nrec,wp)*0.5_wp)
                   CALL lhn_sumrad(nobs,obs_tot_all,sumrad_tot,blacklist_tot,ierr)
                   ntbright=1_iintegers
                ENDIF
                CALL close_file( nulhnrad, yform_read, 0, 0, 1, .FALSE., idbg_level, yerrmsg, ierr)
                IF (ierr /= 0) THEN
                  WRITE(nulhn,*) 'Something wrong by closing radar file with ('//yform_read//')'
                ENDIF
              ENDIF
           ENDIF !min_rad>=60
           WRITE(obsdate,'(6i2.2)') cc_rad,yy_rad,mm_rad,dd_rad,hh_rad,min_rad

           IF (field_num > 1) &
               obs_time(field_num)=obs_time(field_num-1)+NINT(lhn_dt_obs)
        ENDIF ! my_cart_id == 0

        IF (ltime) THEN
          CALL get_timings (i_lhn_obs_prep, ntstep, dt, izerror)
          IF (ltime_barrier) THEN
            CALL comm_barrier (icomm_cart, izerror, yerrmsg)
            CALL get_timings (i_barrier_waiting_lhn, ntstep, dt, izerror)
          ENDIF
        ENDIF

        IF (num_compute > 1) THEN
           CALL distribute_field (obs_tot(1:ie_tot,1:je_tot), &
                ie_tot,je_tot,obs(1:ie,1:je,field_num),ie,je,0,ierr)
           IF (lhn_spqual) &
              CALL distribute_field (spqual_tot(1:ie_tot,1:je_tot), &
                   ie_tot,je_tot,spqual(1:ie,1:je,field_num),ie,je,0,ierr)
        ELSE 
           ! we are running on one PE
           obs(:,:,field_num) = obs_tot(:,:)
           IF (lhn_spqual) &
              spqual(:,:,field_num) = spqual_tot(:,:)
        ENDIF

        IF (num_compute > 1) THEN
           ! obs_time is integer!
           CALL distribute_values (obs_time(field_num),1,0,imp_integers,icomm_cart,ierr)
        ENDIF


        IF (ltime) THEN
            CALL get_timings (i_communications_lhn, ntstep, dt, izerror)
        ENDIF

     ENDDO ! i=iread
   

     IF (ltime) THEN
       CALL get_timings (i_lhn_obs_prep, ntstep, dt, izerror)
       IF (ltime_barrier) THEN
         CALL comm_barrier (icomm_cart, izerror, yerrmsg)
         CALL get_timings (i_barrier_waiting_lhn, ntstep, dt, izerror)
       ENDIF
     ENDIF

! distribute blacklist
     IF (lhn_black) THEN
        IF (num_compute > 1) THEN
            CALL distribute_field (blacklist_tot(1:ie_tot,1:je_tot), &
                 ie_tot,je_tot,blacklist(1:ie,1:je),ie,je,0,ierr)
        ELSE
           ! we are running on one PE
            blacklist(:,:) = blacklist_tot(:,:)
        ENDIF

     ENDIF

! distribute heights
     IF (lhn_height) THEN
        IF (num_compute > 1) THEN
           DO n=1,nradar
            CALL distribute_field (dxheight_tot(1:ie_tot,1:je_tot,n), &
                    ie_tot,je_tot,dxheight(1:ie,1:je,n),ie,je,0,ierr)
           ENDDO
        ELSE
           ! we are running on one PE
           dxheight(:,:,:) = dxheight_tot(:,:,:)
        ENDIF
     ENDIF

     CALL distribute_values (ntbright,1,0,imp_integers,icomm_cart,ierr)

     IF (ltlhnbright .AND. ntbright == 1_iintegers) THEN
        IF (num_compute > 1) THEN
            CALL distribute_field (sumrad_tot(1:ie_tot,1:je_tot), &
                 ie_tot,je_tot,sumrad(1:ie,1:je),ie,je,0,ierr)
        ELSE
            sumrad(:,:)=sumrad_tot(:,:) 
        ENDIF
     ENDIF

     IF (ltime) THEN
         CALL get_timings (i_communications_lhn, ntstep, dt, izerror)
     ENDIF

     IF (ltlhnbright .AND. ntbright == 1_iintegers) THEN
        CALL detect_bright_band (sumrad)
        ntbright = 0_iintegers
     ENDIF

     lhn_obs_read_status=.TRUE.
     next_obs_time=NINT(tnow/lhn_dt_obs)*NINT(lhn_dt_obs)+NINT(lhn_dt_obs)

! set switch for first call to false
     lfirst=.FALSE.
  ENDIF ! lfirst

!-------------------------------------------------------------------------------
! Section 1: Input of observations, distribute field to PE's
!-------------------------------------------------------------------------------
! Section 1a: Determine if new observations have to be read at current timestep
!             Allocate space for observations
!            (remark: The unit used for all time variables is minutes,
!                     time is relative to start of forecast)
!-------------------------------------------------------------------------------

! determine current model time as minutes since start of forecast
  tnow  = ntstep*dt/60.0_wp

  IF (.NOT.lfirst .AND. (tnow >= next_obs_time)) THEN
     lhn_obs_read_status=.TRUE.
     next_obs_time=NINT(tnow/lhn_dt_obs)*NINT(lhn_dt_obs)+NINT(lhn_dt_obs)
!    next_obs_time=next_obs_time+NINT(lhn_dt_obs) ??? 

     IF (field_num == iread) THEN
        field_num = 1_iintegers
     ELSE
        field_num = field_num + 1_iintegers
     ENDIF
     
     IF (my_cart_id == 0) THEN
 
        IF (zlhigh_freq_data) THEN
           obs_time_s=NINT(tnow/lhn_dt_obs)*NINT(lhn_dt_obs)
           obs_time_s=obs_time_s+INT(3*lhn_dt_obs) 
        ELSE
           obs_time_s=NINT(tnow/lhn_dt_obs)*NINT(lhn_dt_obs)
           obs_time_s=obs_time_s+INT(lhn_dt_obs) 
        ENDIF
        
        IF (obsdate == noobs_date(inoobs)) THEN
          WRITE(nulhn, *)'no data are used for ',noobs_date(inoobs)
          obs_tot(:,:)=missval
          inoobs=inoobs+1
          WRITE(nulhn,*)'next noobs:',inoobs,obsdate,'',noobs_date(inoobs-1),'',noobs_date(inoobs)
        ELSE
          obs_tot(:,:)=missval
          spqual_tot(:,:)=missval
          lobsread=.false.
          IF (lhn_spqual) lqualread=.false.
          DO n=1,nrec
            IF ( obsdate == obs_date_read(n)) THEN
              IF ( obs_var_read(n) == 61 .AND. obs_tab_read(n) == 2) THEN
                obs_tot(:,:) = obs_tot_all(:,:,n)
                WRITE(nulhn,*)'next noobs:',inoobs,noobs_date(inoobs)
                lobsread=.true.
              ELSE IF ( obs_var_read(n) == 162 .AND. obs_tab_read(n) == 250) THEN
                spqual_tot(:,:) = obs_tot_all(:,:,n)
                lqualread=.true.
              ELSE
                WRITE(nulhn, *) &
                 'WARNING: observation file does not include precipiation data (ee=61) nor quality data (ee=162) for ',obsdate
                EXIT
!                obs_tot(:,:)    = missval
!                spqual_tot(:,:) = missval
              ENDIF
!            ELSE
!              WRITE(nulhn, *) &
!              'WARNING: observation file does not include precipiation data (ee=61) nor quality data (ee=162) for ',obsdate
!              obs_tot(:,:)    = missval
!              spqual_tot(:,:) = missval
              IF (lobsread .AND. lqualread) EXIT
            ENDIF
          ENDDO
        ENDIF

        IF (lhn_diag) THEN
           ! write availability and statistics of radar data
           znmissval = count(obs_tot/=missval)
           IF (znmissval == 0) THEN
              WRITE(nulhn, *) 'data available: 0 ', obsdate
              WRITE(nulhn, '("Stat. of obs    at ",a12," (min,max,mean): ",3a12)') obsdate,&
                   'undef','undef','undef'
           ELSE
              WRITE(nulhn, *) 'data available: 1 ', obsdate
              WRITE(nulhn, '("Stat. of obs    at ",a12," (min,max,mean): ",3f12.5)') obsdate,&
                   MINVAL(obs_tot,MASK=obs_tot/=missval),&
                   MAXVAL(obs_tot,MASK=obs_tot/=missval),&
                   SUM(obs_tot,MASK=obs_tot/=missval)/REAL(znmissval,wp)
           END IF

           ! write statistics of spatial quality data
           IF (lhn_spqual) THEN
                znmissval = count(spqual_tot/=missval)
                IF (znmissval == 0) THEN
                   WRITE(nulhn, '("Stat. of spqual at ",a12," (min,max,mean): ",3a12)') obsdate,&
                        'undef','undef','undef'
                ELSE
                   WRITE(nulhn, '("Stat. of spqual at ",a12," (min,max,mean): ",3f12.5)') obsdate,&
                        MINVAL(spqual_tot,MASK=spqual_tot/=missval),&
                        MAXVAL(spqual_tot,MASK=spqual_tot/=missval),&
                        SUM(spqual_tot,MASK=spqual_tot/=missval)/REAL(znmissval,wp)
                   WRITE(nulhn,*) "Number of missval in spqual: ",ie_tot*je_tot-znmissval
                END IF
           END IF
        ENDIF

        IF (lhn_spqual) THEN

           ! replace missval of spqual with zero
           WHERE(spqual_tot==missval) spqual_tot = 0.0_wp

           ! Check spqual for values larger than one and smaller than zero
           znmissval = count(spqual_tot>1.0_wp)
           IF (znmissval > 0) THEN
             WRITE(nulhn,*) "WARNING: values larger than one in spqual, correcting them"
             WHERE (spqual_tot > 1.0_wp) spqual_tot = 1.0_wp
           END IF
           znmissval = COUNT(spqual_tot<0.0_wp)
           IF (znmissval > 0) THEN
             WRITE(nulhn,*) "WARNING: values smaller than zero in spqual, correcting them"
             WHERE (spqual_tot < 0.0_wp) spqual_tot = 0.0_wp
           END IF

        END IF

        min_rad=min_rad+NINT(lhn_dt_obs)
        IF (min_rad >= 60) THEN

           min_rad=min_rad - 60
           min_rad=NINT(REAL(min_rad,wp)/(lhn_dt_obs))*NINT(lhn_dt_obs)
           hh_rad = hh_rad + 1
           IF (hh_rad >= 24) THEN
               hh_rad = 0
               dd_rad = dd_rad + 1
               IF (dd_rad > day_of_month(mm_rad) ) THEN
                   dd_rad = 1
                   mm_rad = mm_rad + 1
                   IF (mm_rad > 12) THEN
                       mm_rad = 1
                       yy_rad = yy_rad + 1
                       IF (yy_rad > 99) THEN
                           yy_rad = yy_rad - 100
                           cc_rad = cc_rad + 1
                        ENDIF
                   ENDIF
               ENDIF
           ENDIF
           WRITE(obsdate,'(6i2.2)') cc_rad,yy_rad,mm_rad,dd_rad,hh_rad,min_rad
           CALL lhn_obs_open (nulhnrad,obsdate,ierr)
           IF (ierr == 0 ) THEN
             CALL lhn_obs_read(nulhnrad,nrec_max,obs_tot_all,obs_date_read,obs_var_read,obs_tab_read,nrec,ierr)
             IF (ltlhnbright) THEN
                nobs=nrec
                IF (lhn_spqual) nobs=NINT(REAL(nrec,wp)*0.5_wp)
                CALL lhn_sumrad(nobs,obs_tot_all,sumrad_tot,blacklist_tot,ierr)
                ntbright=1_iintegers
             ENDIF
             CALL close_file( nulhnrad, yform_read, 0, 0, 1, .FALSE., idbg_level, yerrmsg, ierr)
             IF (ierr /= 0) THEN
               WRITE(nulhn,*) 'Something wrong by closing radar file with ('//yform_read//')'
             ENDIF
           ENDIF
        ENDIF
        WRITE(obsdate,'(6i2.2)') cc_rad,yy_rad,mm_rad,dd_rad,hh_rad,min_rad


        obs_time(field_num)=obs_time_s

     ENDIF ! my_cart_id = 0
!-------------------------------------------------------------------------------
! Section 1c: Distribute the field with new obs to different PE's
!-------------------------------------------------------------------------------

     IF (ltime) THEN
       CALL get_timings (i_lhn_obs_prep, ntstep, dt, izerror)
       IF (ltime_barrier) THEN
         CALL comm_barrier (icomm_cart, izerror, yerrmsg)
         CALL get_timings (i_barrier_waiting_lhn, ntstep, dt, izerror)
       ENDIF
     ENDIF

     IF (num_compute > 1) THEN
        CALL distribute_field (obs_tot(1:ie_tot,1:je_tot), &
             ie_tot,je_tot,obs(1:ie,1:je,field_num),ie,je,0,ierr)
        IF (lhn_spqual) &
           CALL distribute_field (spqual_tot(1:ie_tot,1:je_tot), &
                ie_tot,je_tot,spqual(1:ie,1:je,field_num),ie,je,0,ierr)
        CALL distribute_values (ntbright,1,0,imp_integers,icomm_cart,ierr)
        IF (ltlhnbright .AND. ntbright == 1_iintegers) &
           CALL distribute_field (sumrad_tot(1:ie_tot,1:je_tot), &
                ie_tot,je_tot,sumrad(1:ie,1:je),ie,je,0,ierr)
     ELSE 
        ! we are running on one PE
        obs(:,:,field_num) = obs_tot(:,:)
        IF (lhn_spqual) &
           spqual(:,:,field_num) = spqual_tot(:,:)
        IF (ltlhnbright .AND. ntbright == 1_iintegers) &
           sumrad(:,:)=sumrad_tot(:,:) 
     ENDIF
        
!-------------------------------------------------------------------------------
! Section 1d: Update time informaton for observation fields if new record has
!             been read (distribute time of new data record (if parallel mode)
!-------------------------------------------------------------------------------

! distribute the time information to all PE's
     IF (num_compute > 1) THEN
        CALL distribute_values (obs_time(field_num),1,0,imp_integers,icomm_cart,ierr)
     ENDIF
     
     IF (ltime) THEN
         CALL get_timings (i_communications_lhn, ntstep, dt, izerror)
     ENDIF

     IF (ltlhnbright .AND. ntbright == 1_iintegers) THEN
         CALL detect_bright_band (sumrad)
         ntbright = 0_iintegers
     ENDIF
  
  ENDIF   ! tnow >= next_obs_time

! some control output on observations, observation times
  IF (my_cart_id == 0) THEN
     WRITE(nulhn, *)' current time of model run : tnow [min after start] = ',tnow
     WRITE(nulhn, *)' current date of model     : actdate       = ',actdate
     WRITE(nulhn, *)' preset value for time interval '
     WRITE(nulhn, *)'           of observations : lhn_dt_obs   = ',lhn_dt_obs
     WRITE(nulhn, *)' time of next observations : next_obs_time = ',next_obs_time
  ENDIF

!-------------------------------------------------------------------------------
! Section 2: Interpolate observation data in time and assign temporal weight;
!            spatial weights are interpolated likewise (wobs_space)
!            - Interpolate data linearly if consecutive obs are close enough 
!              in time (dtobs < dtobsint) and both obs are present
!              full time weighting factor : wobs_time = 1.
!            - If one observation is missing either at last or next observation
!              time, the available obs is taken but receives a reduced time 
!              weighting factor depending on the distance of the obs to the
!              current timestep and the specified limit of influence (dtinfl)
!            - If both data are present but dtobs > dtobsint both obs are
!              weighted according to the previous case, the total time weight
!              is a weighted mean of both reduced weights
! !! Remark: some computation time may be saved if obs(i,j,inext)-obs(i,j,ilast)
!            are held in and array obsdif(i,j)
!-------------------------------------------------------------------------------

  IF (lhn_obs_read_status) THEN

     DO i=1,iread
        IF (((tnow-obs_time(i)) >= 0_iintegers)&
             .AND.((tnow-obs_time(i)) < lhn_dt_obs)) THEN
           weight_index_0=i
           center_time=obs_time(i)
           next_time_1=center_time+1_iintegers*INT(lhn_dt_obs)
           IF (iread > 2) THEN
            next_time_2=center_time+2_iintegers*INT(lhn_dt_obs)
            next_time_3=center_time+3_iintegers*INT(lhn_dt_obs)
            prev_time_1=center_time-1_iintegers*INT(lhn_dt_obs)
            prev_time_2=center_time-2_iintegers*INT(lhn_dt_obs)
           ENDIF
           EXIT
        ENDIF
     ENDDO
     
     DO i=1,iread
        IF (obs_time(i) == next_time_1) weight_index_p1=i
        IF (iread > 2) THEN
         IF (obs_time(i) == next_time_2) weight_index_p2=i
         IF (obs_time(i) == next_time_3) weight_index_p3=i
         IF (obs_time(i) == prev_time_1) weight_index_m1=i
         IF (obs_time(i) == prev_time_2) weight_index_m2=i
        ENDIF
     ENDDO
  ENDIF

! reset counters
  num1delta_t_obs = 0_iintegers
  num2delta_t_obs = 0_iintegers
  num3delta_t_obs = 0_iintegers
  num4delta_t_obs = 0_iintegers
  numblack        = SUM(NINT(blacklist(istart:iend,jstart:jend)))
  numfull         = 0_iintegers
  numred          = 0_iintegers
  numzero         = 0_iintegers
  numnone         = 0_iintegers

! If the data is in high frequency take into account observations that
! are within the interval [-2,3]*lhn_dt_obs fore the time interpolation
! of obs and wobs_space
IF (zlhigh_freq_data) THEN
  DO   j=jstart,jend
     DO i=istart,iend
!        IF (NINT(blacklist(i,j)) /= 1_iintegers .AND. NINT(brightband(i,j)) /= 1_iintegers) THEN
           IF ((obs(i,j,weight_index_0) >= pr_time_limit) .AND. &
               (obs(i,j,weight_index_p1) >= pr_time_limit)) THEN
              ! observation is valid at t=0 and t=+lhn_dt_obs
              IF (obs(i,j,weight_index_0)  < 0.0_wp) &
                   obs(i,j,weight_index_0) = 0.0_wp
              IF (obs(i,j,weight_index_p1) < 0.0_wp) &
                   obs(i,j,weight_index_p1)= 0.0_wp
              pr_obs(i,j)    = obs(i,j,weight_index_0)                           &
                            + (obs(i,j,weight_index_p1)-obs(i,j,weight_index_0)) &
                            * (tnow-obs_time(weight_index_0))/ lhn_dt_obs
              pr_obs(i,j) = pr_obs(i,j)*sec_per_hr_inv
              wobs_time(i,j) = 1.0_wp
              num1delta_t_obs = num1delta_t_obs + 1
              IF (.NOT.lhn_spqual) CYCLE
              wobs_space(i,j) = spqual(i,j,weight_index_0)                                   &
                             + (spqual(i,j,weight_index_p1)-spqual(i,j,weight_index_0))      &
                             * (tnow-obs_time(weight_index_0))/ lhn_dt_obs
              CYCLE
           ELSEIF((obs(i,j,weight_index_m1) >= pr_time_limit) .AND. &
                  (obs(i,j,weight_index_p1) >= pr_time_limit)) THEN
              ! observation is valid at t=-lhn_dt_obs and t=+lhn_dt_obs
              IF (obs(i,j,weight_index_m1) < 0.0_wp) &
                   obs(i,j,weight_index_m1)= 0.0_wp
              IF (obs(i,j,weight_index_p1) < 0.0_wp) &
                   obs(i,j,weight_index_p1)= 0.0_wp
              pr_obs(i,j)     = obs(i,j,weight_index_m1)                               &
                             + (obs(i,j,weight_index_p1)-obs(i,j,weight_index_m1))  &
                             * (tnow-obs_time(weight_index_m1))/ &
                               (2.0_wp*lhn_dt_obs)
              pr_obs(i,j) = pr_obs(i,j)*sec_per_hr_inv
              wobs_time(i,j) = 0.75_wp   
              num2delta_t_obs = num2delta_t_obs + 1
              IF (.NOT.lhn_spqual) CYCLE
              wobs_space(i,j) = spqual(i,j,weight_index_m1)                               &
                             + (spqual(i,j,weight_index_p1)-spqual(i,j,weight_index_m1))  &
                             * (tnow-obs_time(weight_index_m1))/ &
                               (2.0_wp*lhn_dt_obs)
              CYCLE
           ELSEIF((obs(i,j,weight_index_0) >= pr_time_limit) .AND. &
                  (obs(i,j,weight_index_p2)>= pr_time_limit)) THEN
              ! observation is valid at t=0 and t=+2lhn_dt_obs
              IF (obs(i,j,weight_index_0) < 0.0_wp) &
                   obs(i,j,weight_index_0)= 0.0_wp
              IF (obs(i,j,weight_index_p2) < 0.0_wp) &
                   obs(i,j,weight_index_p2)= 0.0_wp
              pr_obs(i,j) = obs(i,j,weight_index_0)                                   &
                            + (obs(i,j,weight_index_p2)-obs(i,j,weight_index_0))      &
                            * (tnow-obs_time(weight_index_0))/ &
                              (2.0_wp*lhn_dt_obs)
              pr_obs(i,j) = pr_obs(i,j)*sec_per_hr_inv
              wobs_time(i,j) = 0.75_wp
              num2delta_t_obs = num2delta_t_obs + 1
              IF (.NOT.lhn_spqual) CYCLE
              wobs_space(i,j) = spqual(i,j,weight_index_0)                                   &
                             + (spqual(i,j,weight_index_p2)-spqual(i,j,weight_index_0))      &
                             * (tnow-obs_time(weight_index_0))/ &
                               (2.0_wp*lhn_dt_obs)
              CYCLE
           ELSEIF((obs(i,j,weight_index_m1) >= pr_time_limit) .AND. &
                  (obs(i,j,weight_index_p2) >= pr_time_limit)) THEN
              ! observation is valid at t=-lhn_dt_obs and t=+2lhn_dt_obs
              IF (obs(i,j,weight_index_m1) < 0.0_wp) &
                   obs(i,j,weight_index_m1)= 0.0_wp
              IF (obs(i,j,weight_index_p2) < 0.0_wp) &
                   obs(i,j,weight_index_p2)= 0.0_wp
              pr_obs(i,j)    = obs(i,j,weight_index_m1)                                  &
                            + (obs(i,j,weight_index_p2)-obs(i,j,weight_index_m1))     &
                            * (tnow-obs_time(weight_index_m1))/ &
                              (3.0_wp*lhn_dt_obs)
              pr_obs(i,j) = pr_obs(i,j)*sec_per_hr_inv
              wobs_time(i,j) = 0.5_wp
              num3delta_t_obs = num3delta_t_obs + 1
              IF (.NOT.lhn_spqual) CYCLE
              wobs_space(i,j)= spqual(i,j,weight_index_m1)                                  &
                            + (spqual(i,j,weight_index_p2)-spqual(i,j,weight_index_m1))     &
                            * (tnow-obs_time(weight_index_m1))/ &
                              (3.0_wp*lhn_dt_obs)
              CYCLE
           ELSEIF((obs(i,j,weight_index_m2) >= pr_time_limit) .AND. &
                  (obs(i,j,weight_index_p1) >= pr_time_limit)) THEN
              ! observation is valid at t=-2lhn_dt_obs and t=+lhn_dt_obs
              IF (obs(i,j,weight_index_m2) < 0.0_wp) &
                   obs(i,j,weight_index_m2)= 0.0_wp
              IF (obs(i,j,weight_index_p1) < 0.0_wp) &
                   obs(i,j,weight_index_p1)= 0.0_wp
              pr_obs(i,j) = obs(i,j,weight_index_m2)                                  &
                            + (obs(i,j,weight_index_p1)-obs(i,j,weight_index_m2))     &
                            * (tnow-obs_time(weight_index_m2))/ &
                              (3.0_wp*lhn_dt_obs)
              pr_obs(i,j) = pr_obs(i,j)*sec_per_hr_inv
              wobs_time(i,j) = 0.5_wp
              num3delta_t_obs = num3delta_t_obs + 1
              IF (.NOT.lhn_spqual) CYCLE
              wobs_space(i,j) = spqual(i,j,weight_index_m2)                                  &
                             + (spqual(i,j,weight_index_p1)-spqual(i,j,weight_index_m2))     &
                             * (tnow-obs_time(weight_index_m2))/ &
                               (3.0_wp*lhn_dt_obs)
              CYCLE
           ELSEIF((obs(i,j,weight_index_0) >= pr_time_limit) .AND. &
                  (obs(i,j,weight_index_p3)>= pr_time_limit)) THEN
              ! observation is valid at t=0 and t=+3lhn_dt_obs
              IF (obs(i,j,weight_index_0) < 0.0_wp) &
                   obs(i,j,weight_index_0)= 0.0_wp
              IF (obs(i,j,weight_index_p3) < 0.0_wp) &
                   obs(i,j,weight_index_p3)= 0.0_wp
              pr_obs(i,j) = obs(i,j,weight_index_0)                                  &
                            + (obs(i,j,weight_index_p3)-obs(i,j,weight_index_0))     &
                            * (tnow-obs_time(weight_index_0))/ &
                              (3.0_wp*lhn_dt_obs)
              pr_obs(i,j) = pr_obs(i,j)*sec_per_hr_inv
              wobs_time(i,j) = 0.5_wp
              num3delta_t_obs = num3delta_t_obs + 1
              IF (.NOT.lhn_spqual) CYCLE
              wobs_space(i,j) = spqual(i,j,weight_index_0)                                  &
                             + (spqual(i,j,weight_index_p3)-spqual(i,j,weight_index_0))     &
                             * (tnow-obs_time(weight_index_0))/ &
                               (3.0_wp*lhn_dt_obs)
              CYCLE
           ELSEIF((obs(i,j,weight_index_m2) >= pr_time_limit) .AND. &
                  (obs(i,j,weight_index_p2) >= pr_time_limit)) THEN
              ! observation is valid at t=-2lhn_dt_obs and t=+2lhn_dt_obs
              IF (obs(i,j,weight_index_m2) < 0.0_wp) &
                   obs(i,j,weight_index_m2)= 0.0_wp
              IF (obs(i,j,weight_index_p2) < 0.0_wp) &
                   obs(i,j,weight_index_p2)= 0.0_wp
              pr_obs(i,j) = obs(i,j,weight_index_m2)                                  &
                            + (obs(i,j,weight_index_p2)-obs(i,j,weight_index_m2))     &
                            * (tnow-obs_time(weight_index_m2))/ &
                              (4.0_wp*lhn_dt_obs)
              pr_obs(i,j) = pr_obs(i,j)*sec_per_hr_inv
              wobs_time(i,j) = 0.25_wp
              num4delta_t_obs = num4delta_t_obs + 1
              IF (.NOT.lhn_spqual) CYCLE
              wobs_space(i,j) = spqual(i,j,weight_index_m2)                                  &
                             + (spqual(i,j,weight_index_p2)-spqual(i,j,weight_index_m2))     &
                             * (tnow-obs_time(weight_index_m2))/ &
                               (4.0_wp*lhn_dt_obs)
              CYCLE
           ELSEIF((obs(i,j,weight_index_m1) >= pr_time_limit) .AND. &
                  (obs(i,j,weight_index_p3) >= pr_time_limit)) THEN
              ! observation is valid at t=-lhn_dt_obs and t=+3lhn_dt_obs
              IF (obs(i,j,weight_index_m1) < 0.0_wp) &
                   obs(i,j,weight_index_m1)= 0.0_wp
              IF (obs(i,j,weight_index_p3) < 0.0_wp) &
                   obs(i,j,weight_index_p3)= 0.0_wp
              pr_obs(i,j) = obs(i,j,weight_index_m1)                                  &
                            + (obs(i,j,weight_index_p3)-obs(i,j,weight_index_m1))     &
                            * (tnow-obs_time(weight_index_m1))/ &
                              (4.0_wp*lhn_dt_obs)
              pr_obs(i,j) = pr_obs(i,j)*sec_per_hr_inv
              wobs_time(i,j) = 0.25_wp
              num4delta_t_obs = num4delta_t_obs + 1
              IF (.NOT.lhn_spqual) CYCLE
              wobs_space(i,j) = spqual(i,j,weight_index_m1)                                  &
                             + (spqual(i,j,weight_index_p3)-spqual(i,j,weight_index_m1))     &
                             * (tnow-obs_time(weight_index_m1))/ &
                              (4.0_wp*lhn_dt_obs)
              CYCLE
           ELSE
              ! observation is not valid
              pr_obs(i,j) = pr_time_limit
              wobs_space(i,j) = 0.0_wp
              wobs_time(i,j) = 0.0_wp
              numnone = numnone + 1
           ENDIF
!        ELSE
!           ! observation is blacklisted or bright band
!           pr_obs(i,j) = pr_time_limit
!           wobs_space(i,j) = 0.0_wp
!           wobs_time(i,j) = 0.0_wp
!           IF (NINT(blacklist(i,j)) == 1_iintegers ) &
!             numblack = numblack + 1_iintegers
!        ENDIF
     ENDDO
  ENDDO

ELSE

  ! use only this and the next observation
  ! determine some constant variables/factors needed in the following loops
  ! time weighting factors when one or both observations have reduced influence
  twlast = MAX (0._wp, (obs_time(weight_index_p1) - tnow)/ &
               (obs_time(weight_index_p1)-obs_time(weight_index_0)) ) ! weight for last obs
  twnext = MAX (0._wp, (tnow - obs_time(weight_index_0) )/ &
               (obs_time(weight_index_p1)-obs_time(weight_index_0)) ) ! weight for next obs

! linear interpolation between last and next obs (if close enough in time)
    DO   j=jstart,jend
      DO i=istart,iend
        IF ((obs(i,j,weight_index_0) >= pr_time_limit) .AND. &
            (obs(i,j,weight_index_p1) >= pr_time_limit)) THEN
           IF (obs(i,j,weight_index_0)  < 0.0_wp) &
                obs(i,j,weight_index_0) = 0.0_wp
           IF (obs(i,j,weight_index_p1) < 0.0_wp) &
                obs(i,j,weight_index_p1)= 0.0_wp
           pr_obs(i,j) = obs(i,j,weight_index_0)                                   &
                         + (obs(i,j,weight_index_p1)-obs(i,j,weight_index_0))      &
                         * (tnow-obs_time(weight_index_0))/ lhn_dt_obs
           pr_obs(i,j) = pr_obs(i,j)*sec_per_hr_inv
           wobs_time(i,j) = 1.0_wp
           num1delta_t_obs = num1delta_t_obs + 1
           IF (.NOT.lhn_spqual) CYCLE
           wobs_space(i,j) = spqual(i,j,weight_index_0)                               &
                          + (spqual(i,j,weight_index_p1)-spqual(i,j,weight_index_0))             &
                          * (tnow-obs_time(weight_index_0))/ lhn_dt_obs
           CYCLE
         ELSEIF(obs(i,j,weight_index_0) >= pr_time_limit) THEN
           IF (obs(i,j,weight_index_0)  < 0.0_wp) &
               obs(i,j,weight_index_0) = 0.0_wp
           pr_obs(i,j)    = obs(i,j,weight_index_0) 
           wobs_time(i,j) = twlast   
           IF (.NOT.lhn_spqual) CYCLE
           wobs_space(i,j) = spqual(i,j,weight_index_0)
           CYCLE
         ELSEIF(obs(i,j,weight_index_p1) >= pr_time_limit) THEN
            pr_obs(i,j) = obs(i,j,weight_index_p1)
            wobs_time(i,j) = twnext
            IF (.NOT.lhn_spqual) CYCLE
            wobs_space(i,j) = spqual(i,j,weight_index_p1)
            CYCLE
         ELSE
            pr_obs(i,j) = pr_time_limit
            wobs_space(i,j) = 0.0_wp
            wobs_time(i,j) = 0.0_wp
            numnone = numnone + 1
         ENDIF
      ENDDO
    ENDDO
ENDIF

! if no spatial quality function is used, set wobs_space constant to one
IF (.NOT.lhn_spqual) wobs_space(istart:iend,jstart:jend) = 1.0_wp

! determine statistics about spatial weights

numfull = count(  wobs_space(istart:iend,jstart:jend) == 1.0_wp)
numred  = count( (wobs_space(istart:iend,jstart:jend) <  1.0_wp) &
           .and. (wobs_space(istart:iend,jstart:jend) >  0.0_wp) )
numzero = count(  wobs_space(istart:iend,jstart:jend) == 0.0_wp)

! some control output on time and spatial weighting
!-------------------------------------------------------------------------------

  IF (lhn_diag) THEN
   IF (num_compute > 1) THEN
     ! get summed diagnostics at PE0 from all PE's using collect_values 
     ! (collect_values needs real vector as input)
     iexchange = 0
     iexchange( 1)= num1delta_t_obs
     iexchange( 2)= num2delta_t_obs
     iexchange( 3)= num3delta_t_obs
     iexchange( 4)= num4delta_t_obs
     iexchange( 5)= numnone
     iexchange( 6)= numblack
     iexchange( 7)= numfull
     iexchange( 8)= numred
     iexchange( 9)= numzero

     IF (ltime) THEN
       CALL get_timings (i_lhn_obs_prep, ntstep, dt, izerror)
       IF (ltime_barrier) THEN
         CALL comm_barrier (icomm_cart, izerror, yerrmsg)
         CALL get_timings (i_barrier_waiting_lhn, ntstep, dt, izerror)
       ENDIF
     ENDIF

     CALL global_values (iexchange,9, 'SUM',imp_integers,icomm_cart, 0,yerrmsg, izerror)

     IF (ltime) THEN
        CALL get_timings (i_communications_lhn, ntstep, dt, izerror)
     ENDIF

   ENDIF

   IF (my_cart_id == 0) THEN
     IF (num_compute > 1) THEN 
        num1delta_t_obs    = iexchange( 1)
        num2delta_t_obs    = iexchange( 2)
        num3delta_t_obs    = iexchange( 3)
        num4delta_t_obs    = iexchange( 4)
        numnone            = iexchange( 5)        
        numblack           = iexchange( 6)
        numfull            = iexchange( 7)
        numred             = iexchange( 8)
        numzero            = iexchange( 9)
     ENDIF

   ! remark : total number of points is (ie_tot-2*nboundlines)*(je_tot-2*nboundlines)
     WRITE(nulhn, *) 
     WRITE(nulhn, *)' Diagnostics of RADAR obs time interpolation, lhn_dt_obs = ',lhn_dt_obs
     WRITE(nulhn, *)' number of treated obs points in time weighting : ',  &
          num1delta_t_obs+num2delta_t_obs+num3delta_t_obs+num4delta_t_obs+numnone
     WRITE(nulhn, *)' n of points with obs-dist 1 delta_t  : num1delta_t_obs = ',num1delta_t_obs
     WRITE(nulhn, *)' n of points with obs-dist 2 delta_t  : num2delta_t_obs = ',num2delta_t_obs
     WRITE(nulhn, *)' n of points with obs-dist 3 delta_t  : num3delta_t_obs = ',num3delta_t_obs
     WRITE(nulhn, *)' n of points with obs-dist 4 delta_t  : num4delta_t_obs = ',num4delta_t_obs
     WRITE(nulhn, *)' n of points without obs              : numnone = ',numnone     
     WRITE(nulhn, *)' n of points which are blacklisted    : numblack = ',numblack
     WRITE(nulhn, *)
     WRITE(nulhn, *)' Diagnostics of RADAR obs space weighting'
     WRITE(nulhn, *)' n of points with full    spatial weight : numfull = ',numfull
     WRITE(nulhn, *)' n of points with reduced spatial weight : numred  = ',numred
     WRITE(nulhn, *)' n of points with zero    spatial weight : numzero = ',numzero
     WRITE(nulhn, *)
   ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! Section 5: Deallocate space of temporal arrays
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! End of subroutine 
!-------------------------------------------------------------------------------

END SUBROUTINE lhn_obs_prep

!===============================================================================
!+ Module procedure in "lheat_nudge" reading radar precip data from input file
!-------------------------------------------------------------------------------
 
SUBROUTINE lhn_obs_open (nulhnrad,yobs_date,ierr)

!-------------------------------------------------------------------------------
!
! Description:
!   This procedure in the module "lheat_nudge" is called by "lhn_obs_prep" and
!   opens the radar data input file yradar_inpt and reads the general header
!   information.
!
! Method:
!   Call of GRIB I/O routines of DWDLIB or grib_api.
!
! Input files:
!   GRIB file.
!   
!-------------------------------------------------------------------------------

! Subroutine scalars (intent: out):
!----------------------------------

  INTEGER (KIND=iintegers), INTENT(OUT)            ::       &
    nulhnrad         ,& ! unit number for radardata GRIB file
    ierr             ! status/error code

  CHARACTER (LEN=12), INTENT(IN)                   ::       &
    yobs_date

! Local scalars, arrays:
!-------------------------------------------------------------------------------

  CHARACTER (LEN=100)    ::  ypath
  CHARACTER (LEN= 14)    ::  yradar_inpt

  LOGICAL                ::  lzexist1=.FALSE., lzexist2=.FALSE.

!- End of header
!===============================================================================
!-------------------------------------------------------------------------------
! Begin of subroutine
!-------------------------------------------------------------------------------
  yroutine = 'lhn_obs_open'

! Open input file and read the header information
!-------------------------------------------------------------------------------
! open input file
  nulhnrad=0_iintegers
  ierr=0_iintegers

! Check, which file exists: either '.grib1'  (for backwards compatibility) or only '.grib'
! (at DWD the following convention is used:
!  - for yform_read='grb1': the preprocessing job prepares a filename with '.grib1'
!  - for yform_read='apix': the preprocessing job prepares a filename with '.grib'  )

  yradar_inpt = yobs_date(3:10)//'.grib1'
  ypath       = TRIM(radar_in)//TRIM(yradar_inpt)

  INQUIRE (FILE=ypath, EXIST=lzexist1)

  IF (.NOT. lzexist1) THEN
    ! check wether a .grib2 file exists
    yradar_inpt = yobs_date(3:10)//'.grib'
    ypath       = TRIM(radar_in)//yradar_inpt(1:13)

    INQUIRE (FILE=TRIM(ypath), EXIST=lzexist2)
    IF (.NOT. lzexist2) THEN
      PRINT *,  ' *** ERROR LHN: no radar input file available:  ', ypath
      PRINT *,  ' Please note file name convention:'
      PRINT *,  '     using DWDLIB: *.grib1'
      PRINT *,  '     using GRIB_API: *.grib'
      WRITE(nulhn, *) ' *** ERROR LHN: no radar input file available:  ', ypath
      WRITE(nulhn, *) ' Please note file name convention for LHN input files:'
      WRITE(nulhn, *) '     using DWDLIB (only grib1): *.grib1'
      WRITE(nulhn, *) '     using GRIB_API (any case): *.grib'
!      yerrmsg = 'Error in LHN: no radar input file available'
!      CALL model_abort(my_cart_id,6010,yerrmsg,yroutine)
      ierr=1
      RETURN
    ENDIF
  ENDIF


! if this part is reached, the file ypath exists and can be opened (whether .grib1 or .grib)

  PRINT *, 'LHN:  I open input file: ',TRIM(ypath)
  WRITE(nulhn, *)'I open input file: ',TRIM(ypath)

  IF (yform_read == 'grb1' .AND. .NOT. lzexist1) THEN
    ! the radar input file could be in GRIB2, then the reading with DWDLIB will fail
    PRINT *,  ' *** WARNING LHN: GRIB format not known: could cause problems with DWDLIB'
    WRITE(nulhn, *) ' *** WARNING LHN: GRIB format not known: could cause problems with DWDLIB'
  END IF

  CALL open_file( nulhnrad, ypath, 'r  ', yform_read, 0, 0, 1, .FALSE., idbg_level, yerrmsg, ierr )
  IF (ierr /= 0) THEN
     PRINT *, 'Error in opening the radar data file: ', yradar_inpt
     WRITE(nulhn, *)'Error in opening the radar data file: ', yradar_inpt
! fuo: we don't want cosmo to stop if no LHN files are being found
!     yerrmsg='Error in opening the radar data file.'
!     CALL model_abort(my_cart_id,6010,yerrmsg,yroutine)
  ENDIF

!-------------------------------------------------------------------------------
! End of subroutine 
!-------------------------------------------------------------------------------

END SUBROUTINE lhn_obs_open

!===============================================================================
!+ Module procedure in "lheat_nudge" reading radar precip data from input file
!-------------------------------------------------------------------------------
 
!SUBROUTINE lhn_obs_read (nulhnrad,itype,imocc,imoyy,imomm,imodd,imohh,imomin,datafield,ierr)
SUBROUTINE lhn_obs_read (nulhnrad,nrecs,datafield,cdate_rad,ivar,itab,nread,ierr_cnt)

!-------------------------------------------------------------------------------
!
! Description:
!   This procedure in the module "lheat_nudge" is called by "lhn_obs_prep" and
!   reads one header and one data record data from the radar data input file.
!
! Method:
!   Simple read from sequential unformatted file.
!
! Input files:
!   Sequential unformatted file YRADAR.
!
!-------------------------------------------------------------------------------

! Subroutine arguments: 
!-------------------------------------------------------------------------------
! Subroutine scalars (intent: in):
!-----------------------------------

  INTEGER (KIND=iintegers), INTENT(IN)              ::       &
    nrecs !

  INTEGER (KIND=iintegers), INTENT(INOUT)          ::       &
    nulhnrad            ! unit number for radardata GRIB file

! Subroutine scalars (intent: out):
!----------------------------------

!  INTEGER (KIND=iintegers), INTENT(IN)            ::       &
!    itype         ,& ! type of variable. 0 for obs, 1 for spatial quality
!    imocc         ,& ! date of requested observations : century
!    imoyy         ,& ! date of requested observations : year
!    imomm         ,& ! date of requested observations : month
!    imodd         ,& ! date of requested observations : day
!    imohh         ,& ! date of requested observations : hour
!    imomin           ! date of requested observations : minute

  CHARACTER (LEN=12), INTENT(OUT)           ::        &
    cdate_rad(nrecs)     ! array of date of record

  INTEGER (KIND=iintegers), INTENT(OUT)           ::       &
    ivar(nrecs)       ,& ! array of grib element numbers
    itab(nrecs)       ,& ! array of grib table numbers
    nread             ,& ! index of read record
    ierr_cnt             ! status/error code

! Subroutine arrays (intent: out):
!---------------------------------

  REAL (KIND=wp),     INTENT(OUT)                   ::       &
    datafield(ie_tot,je_tot,nrecs)      ! field for decoded data

! Local parameters, scalars, arrays :
!-------------------------------------------------------------------------------
! Local scalars :
!----------------
  INTEGER (KIND=iintegers)             ::       &
    i,j             ,& ! loop counter
    iyy,imm,idd     ,& ! date of observational data
    ihh,imin        ,& ! time of observational data
    icc,iccyy       ,& ! variables for the computation of the year of observations
    ind             ,& ! index
    ilen            ,& !
    ierr
!   zgrib_nr        ,& ! grib_nr of expected data record
!   ztable_nr          ! table_nr of expected data record

  INTEGER (KIND=intgribf)              ::       &
    igribid, ieditionnr, ie_grib, je_grib, igriblen

  REAL (KIND=wp)                       ::       &
    zstartlon_grib, zendlon_grib, zstartlat_grib, zendlat_grib

  CHARACTER (LEN= 30)                  ::       &
    ydate, ytime, yzshortname

!  LOGICAL                              ::       &
!    lmatch_date   = .false.,& ! true if grib record has right date
!    lmatch_var    = .false.,& ! true if grib record has right variable
!    lmatch_domain = .false.,& ! true if grib record has right domain dimensions
!    lrewind

!  LOGICAL, SAVE                        ::       &
!    lfirst=.TRUE.       ! indicator for steps to be executed at first call only

!- End of header
!===============================================================================

!-------------------------------------------------------------------------------
! Begin of subroutine
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Section 0: Initialization of some variables
!-------------------------------------------------------------------------------
  yroutine = 'lhn_obs_read'
  ierr     = 0_iintegers

  datafield = missval
  ierr_cnt  = 0_iintegers
  nread     = 0_iintegers
  ivar      = 0_iintegers
  itab      = 0_iintegers
  cdate_rad = "190001010100"

!  lrewind  = .false.

  ! define grib specification of radar data
!  IF (itype == 0) THEN
!     ! grib specification of observations
!     zgrib_nr  = 61
!     ztable_nr = 2
!  ELSEIF (itype == 1) THEN
!     ! grib specification of quality function
!     zgrib_nr  = 162
!     ztable_nr = 250
!  ELSE
!!     WRITE(nulhn,'("ERROR: unexpected itype: ",i3,", expected 0 or 1")')& &,itype
!     WRITE(nulhn,*)'ERROR: unexpected itype: ',itype,', expected 0 or 1'
!     ierr = 1
!     RETURN
!  ENDIF

  IF (nulhnrad == 0) RETURN !there is no data file open yet
    ! what about grib_api, could that be 0 there?

! Read data for one data record
!-------------------------------------------------------------------------------

! read time information
! read one record of the GRIB file

  maxlenc=INT(lfd*iwlength, intgribc)
  yzshortname(:) = ' '


! Initializations for the grib library
! (input product definition section)

     idims=0_iintegers
     idims(1) = npds
     idims(2) = ngds
     idims(3) = nbms
     idims(4) = nbds
     idims(5) = -1
     idims(6) = ndsup
     idims(7) = lds
     idims(8) = lfd

  readloop: DO !until no more record
     
     SELECT CASE (yform_read)
     CASE ('grb1')
#ifdef GRIBDWD
       CALL cuegin(nulhnrad, maxlenc, iblock_rad, ilenc, ierrc)

       ilen=INT(ilenc, iintegers)
       ierr=INT(ierrc, iintegers)
       IF (ilen /= 0 .AND. ierr /= 0) THEN
          WRITE(nulhn, *) 'WARNING: error in reading from the GRIB radar data file: ntstep: ', ntstep
          IF (ierr == -2) THEN
             WRITE(nulhn, *)'Field length of the GRIB field not correct.'
          ENDIF
          CYCLE readloop
       ENDIF

#endif

     CASE ('apix')
#ifdef GRIBAPI
       CALL grib_new_from_file (nulhnrad, igribid, ierrc)

       IF ( ierrc /= GRIB_SUCCESS ) THEN
         IF (ierrc == GRIB_END_OF_FILE) THEN
           !EOF is reached 
           ilen = 0
         ELSE
           ! cycle readloop analogous as for DWDLIB
           CYCLE readloop
           !US ierr = 2
         ENDIF
       ELSE
         ! get edition number
         CALL grib_get (igribid, 'editionNumber', ieditionnr, ierrc)
         IF ( ierrc /= GRIB_SUCCESS ) THEN
           WRITE (nulhn, *) 'Error in reading edition number with grib_get'
           ierr = 3
         ENDIF

         ! get total length of message
         CALL grib_get (igribid, 'totalLength', ilen, ierrc)
         IF ( ierrc /= GRIB_SUCCESS ) THEN
           WRITE (nulhn, *) 'Error in reading totalLength with grib_get'
           ierr = 3
         ENDIF
       ENDIF
#endif
     END SELECT

     IF (ilen == 0_iintegers) THEN
        WRITE(nulhn, *)'end of obs file reached: ntstep: ', ntstep
#ifdef GRIBAPI
        IF (yform_read == 'apix') THEN
          ! The handle can be released now
          CALL grib_release(igribid)
        ENDIF
#endif
        EXIT readloop
     ENDIF
     nread=nread+1

     IF (nread > nrecs) THEN
        PRINT *, 'WARNING (LHN): number of records in radar grib_file exceed expected number',nrecs,'!'
        PRINT *, '               continue model run without further reading!'

        WRITE(nulhn, *) 'WARNING (LHN): number of records in radar grib_file exceed expected number',nrecs,'!'
        WRITE(nulhn, *) '               continue model run without further reading!'
        RETURN
     ENDIF

     SELECT CASE (yform_read)
     CASE ('grb1')
#ifdef GRIBDWD
       CALL grbin1(idwdednr, REAL(missval,irealgrib), ndims, idims, iblock_rad, ibmap_rad, &
                   ipds, igds, ibms, ibds, dsup_rad, ds_rad, ierr)
       IF (ierr /= 0) THEN
          WRITE(nulhn, *)'WARNING: error in decoding the record and ipds of the GRIB radar data file'
          ierr_cnt = ierr_cnt + ierr
       ENDIF

       ! Get metadata for checking the date and time
       icc  = ipds(22)-1
       iyy  = ipds(11)
       iccyy= iyy + icc*100
       imm  = ipds(12)
       idd  = ipds(13)
       ihh  = ipds(14)
       imin = ipds(15)
       WRITE(cdate_rad(nread),'(6i2.2)') icc,iyy,imm,idd,ihh,imin

       ! Get metadata for checking the variable
       ivar(nread)=INT(ipds(7),iintegers)
       itab(nread)=INT(ipds(2),iintegers)

       ! Get metadata for checking the grid
       ie_grib        =      igds( 5)
       je_grib        =      igds( 6)
       zstartlat_grib = REAL(igds( 7), wp)*0.001_wp
       zstartlon_grib = REAL(igds( 8), wp)*0.001_wp
       zendlat_grib   = REAL(igds(10), wp)*0.001_wp
       zendlon_grib   = REAL(igds(11), wp)*0.001_wp
#endif

     CASE ('apix')
#ifdef GRIBAPI
       ! Get metadata for checking the date and time
       CALL grib_get (igribid, 'dataDate', ydate, ierr)
       CALL grib_get (igribid, 'dataTime', ytime, ierr)
       cdate_rad(nread) = ydate(1:8)//ytime(1:4)

       IF     (ieditionnr == 1) THEN
         ! Get metadata for checking the variable
         CALL grib_get (igribid, 'indicatorOfParameter', ivar(nread), ierr)
         CALL grib_get (igribid, 'table2Version',        itab(nread), ierr)
       ELSEIF (ieditionnr == 2) THEN

         CALL grib_get (igribid, 'shortName',            yzshortname, ierr)
         IF     (TRIM(yzshortname) == 'RAD_PRECIP') THEN
           ivar(nread) = 61
           itab(nread) =  2
         ELSEIF (TRIM(yzshortname) == 'RAD_QUAL') THEN
           ivar(nread) = 162
           itab(nread) = 250
         ELSE
           PRINT *, 'ERROR: *** this shortName is not recognized in Latent Heat Nudging ***'
           CALL model_abort (my_cart_id,6411,'wrong shortname for LHN', 'lhn_obs_read')
         ENDIF
       ENDIF

       ! Get metadata for checking the grid
       CALL grib_get (igribid, 'Ni',               ie_grib, ierr)
       CALL grib_get (igribid, 'Nj',               je_grib, ierr)

       CALL grib_get (igribid, 'latitudeOfFirstGridPointInDegrees',  zstartlat_grib,  ierr)
       CALL grib_get (igribid, 'longitudeOfFirstGridPointInDegrees', zstartlon_grib,  ierr)
       CALL grib_get (igribid, 'latitudeOfLastGridPointInDegrees',   zendlat_grib,    ierr)
       CALL grib_get (igribid, 'longitudeOfLastGridPointInDegrees',  zendlon_grib,    ierr)

       ! Careful: for grib2, longitudes are between    0...360 degrees,
       !          for grib1, longitudes are between -180...180
       ! the COSMO convention still is -180...180
       IF (ieditionnr == 2) THEN
         IF (zstartlon_grib >= 180.0_wp) THEN
           zstartlon_grib = zstartlon_grib - 360.0_wp
         ENDIF
         IF (zendlon_grib >  180.0_wp) THEN
           zendlon_grib = zendlon_grib - 360.0_wp
         ENDIF
       ENDIF
#endif
     END SELECT

     ! perform some checks of correctness of radar grib input
     !-------------------------------------------------------

     ! check for geographical domain
     IF (                                                                                  &
          ! lower left corner
          (ABS(zstartlat_grib - startlat_tot)                     <= 1.0E-3_wp) .AND.  &
          (ABS(zstartlon_grib - startlon_tot)                     <= 1.0E-3_wp) .AND.  &
          ! upper right corner
          (ABS(zendlat_grib   - (startlat_tot + (je_tot-1)*dlat)) <= 1.0E-3_wp) .AND.  &
          (ABS(zendlon_grib   - (startlon_tot + (ie_tot-1)*dlon)) <= 1.0E-3_wp) .AND.  &
          ! dimensions
          (ie_tot == INT(ie_grib, iintegers))                                       .AND.  &
          (je_tot == INT(je_grib, iintegers))                                              &
          ) THEN
!       lmatch_domain = .true.
     ELSE
!       lmatch_domain = .false.
        WRITE(nulhn, *) "WARNING: Unexpected domain of radar observations"
        WRITE(nulhn, '("WARNING:   start_lat =",f9.5,", expected ",f9.5)') zstartlat_grib, startlat_tot
        WRITE(nulhn, '("WARNING:   start_lon =",f9.5,", expected ",f9.5)') zstartlon_grib, startlon_tot
        WRITE(nulhn, '("WARNING:     end_lat =",f9.5,", expected ",f9.5)') zendlat_grib  , startlat_tot + (je_tot-1)*dlat
        WRITE(nulhn, '("WARNING:     end_lon =",f9.5,", expected ",f9.5)') zendlon_grib  , startlon_tot + (ie_tot-1)*dlon
        WRITE(nulhn, '("WARNING:      ie_tot =",I5  ,", expected ",I5)'  ) ie_grib       , ie_tot
        WRITE(nulhn, '("WARNING:      je_tot =",I5  ,", expected ",I5)'  ) je_grib       , ie_tot
     END IF

     ! if all matches, then read out data
!     IF (lmatch_date .AND. lmatch_var .AND. lmatch_domain) THEN
#ifdef GRIBAPI
       IF (yform_read == 'apix') THEN
         ! Get size of data and data itself
         ! Set missing_value before
         CALL grib_set (igribid, 'missingValue', missval)
         CALL grib_get_size (igribid, 'values', igriblen, ierr)

         IF (igriblen > lds) THEN
           WRITE (nulhn,*) ' *** ERROR: size of message is too big for allocated field: ', igriblen, lds
           ierr = 6
         ENDIF

         CALL grib_get (igribid, 'values', ds_rad, ierr)
       ENDIF
#endif
        DO j=1,je_tot
           DO i=1,ie_tot
              ind=(j-1)*ie_tot+i
              IF (ds_rad(ind) == REAL(missval,irealgrib)) THEN
                 datafield(i,j,nread) = missval
              ELSE
                 datafield(i,j,nread) = REAL(ds_rad(ind),wp)
              ENDIF
           ENDDO
        ENDDO
!        RETURN
!     END IF

#ifdef GRIBAPI
     IF (yform_read == 'apix') THEN
       ! The handle can be released now
       CALL grib_release(igribid)
     ENDIF    
#endif

  ENDDO readloop

  WRITE(nulhn, *)nread,' obs-fields found: ntstep: ', ntstep

!-------------------------------------------------------------------------------
! End of subroutine 
!-------------------------------------------------------------------------------

END SUBROUTINE lhn_obs_read
 
SUBROUTINE lhn_sumrad (nrec,datafield_all,datafield,blacklist,ierr)

!-------------------------------------------------------------------------------
!
! Description:
!   This procedure in the module "lheat_nudge" is called by "lhn_obs_prep" and
!   reads one header and all data record data from the radar data input file.
!   It calculates the sum of all records and a local ratio between the local sum
!   and the areal sum.
!   This procedure is required to run the bright_band_detection
!
! Method:
!   Simple read from sequential unformatted file.
!
! Input files:
!   Sequential unformatted file YRADAR.
!
!-------------------------------------------------------------------------------

! Subroutine arguments: 
!-------------------------------------------------------------------------------
! Subroutine scalars (intent: in):
!-----------------------------------

  INTEGER (KIND=iintegers), INTENT(INOUT)          ::       &
    nrec               ! number for read data records

! Subroutine scalars (intent: out):
!----------------------------------

  INTEGER (KIND=iintegers), INTENT(OUT)           ::       &
    ierr             ! status/error code

! Subroutine arrays (intent: out):          
!---------------------------------

  REAL (KIND=wp),     INTENT(IN)                   ::       &
    datafield_all(ie_tot,je_tot,nrec),blacklist(ie_tot,je_tot)        ! field for decoded data

  REAL (KIND=wp),     INTENT(OUT)                   ::       &
    datafield(ie_tot,je_tot)        ! field for decoded data

! Local parameters, scalars, arrays :
!-------------------------------------------------------------------------------
! Local scalars :
!----------------
  INTEGER (KIND=iintegers)             ::       &
    i,j,n,                & ! loop counter
    nsum,                 & ! counters
    ndata(ie_tot,je_tot)

  REAL (KIND=wp)                       ::       &
    datasum

!- End of header
!===============================================================================
!-------------------------------------------------------------------------------
! Begin of subroutine
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Section 0: Initialization of some variables
!-------------------------------------------------------------------------------
  yroutine = 'lhn_sumrad'
  ierr     = 0_iintegers

  datafield = 0.0_wp
  datasum   = 0.0_wp
  nsum      = 0_iintegers
  ndata     = 0_iintegers

  IF (nrec == 0) RETURN

  DO j=1,je_tot
    DO i=1,ie_tot
      IF (NINT(blacklist(i,j), iintegers) /= 1_iintegers ) THEN
        DO n=1,nrec
            IF (datafield_all(i,j,n) >= 0.0_wp) THEN
              datafield(i,j) = datafield(i,j) + datafield_all(i,j,n)
              ndata(i,j) = ndata(i,j) + 1_iintegers
            ENDIF
        ENDDO
        IF (ndata(i,j) < 1) CYCLE
        datafield(i,j)=datafield(i,j)/REAL(ndata(i,j),wp)
        IF (datafield(i,j) > 2.5_wp) THEN
           datasum=datasum+datafield(i,j)
           nsum=nsum+1
        ENDIF
      ENDIF
    ENDDO
  ENDDO

  IF (nsum < 1_iintegers ) THEN
     datafield(:,:)=0.0_wp
  ELSE  ! nsum > 0
    datasum=datasum / REAL(nsum,wp)
    IF (datasum > 0.0_wp) datafield(:,:)=datafield(:,:)/datasum
  ENDIF

!-------------------------------------------------------------------------------
! End of subroutine 
!-------------------------------------------------------------------------------

END SUBROUTINE lhn_sumrad

!===============================================================================
!+ Module procedure in "lheat_nudge" reading blacklisted radar pixel
!-------------------------------------------------------------------------------

SUBROUTINE lhn_blacklist_dx_open_and_read (nulhnblackdx,blacklist_file,imocc,&
                                           imoyy,imomm,imodd,imohh,imomin,blacklist_tot,ierr)
!-------------------------------------------------------------------------------
!
! Description:
!   This procedure in the module "lheat_nudge" is called by "lhn_obs_prep" and
!   opens and reads the blacklist file for the DX data.
!
! Method:
!   Call of GRIB I/O routines of DWDLIB.
!
! Input files:
!   GRIB file.
!  
!-------------------------------------------------------------------------------

! Subroutine scalars (intent: out):
!----------------------------------

  INTEGER (KIND=iintegers), INTENT(OUT)            ::       &
    nulhnblackdx         ,& ! unit number for radardata GRIB file
    ierr             ! status/error code

  CHARACTER (LEN=*), INTENT(IN)                   ::       &
    blacklist_file


! Subroutine scalars (intent: out):
!----------------------------------

  INTEGER (KIND=iintegers), INTENT(IN)            ::       &
    imocc         ,& ! date of requested observations : century
    imoyy         ,& ! date of requested observations : year
    imomm         ,& ! date of requested observations : month
    imodd         ,& ! date of requested observations : day
    imohh         ,& ! date of requested observations : hour
    imomin           ! date of requested observations : minute

! Subroutine arrays (intent: out):
!---------------------------------

  REAL (KIND=wp),     INTENT(OUT)                   ::       &
    blacklist_tot(ie_tot,je_tot)        ! field for decoded data

! Local parameters, scalars, arrays :
!-------------------------------------------------------------------------------
! Local scalars :
!----------------
  INTEGER (KIND=iintegers)             ::       &
    i,j             ,& ! loop counter
    ind             ,& ! index
    igribid, iee    ,& ! 
    ieditionnr      ,& ! 
    igriblen        ,& ! 
    ilen               !

  CHARACTER (LEN=100)    ::  path
  CHARACTER (LEN= 30)    ::  yzshortname

  LOGICAL :: lb_exist

!- End of header
!===============================================================================

!-------------------------------------------------------------------------------
! Begin of subroutine
!-------------------------------------------------------------------------------
  yroutine = 'lhn_blacklist_dx_open_and_read'

! Open input file and read the header information
!-------------------------------------------------------------------------------
! open input file
  nulhnblackdx=0_iintegers
  yzshortname(:) = ' '

! open GRIB file

  path=radar_in(1:LEN_TRIM(radar_in))//blacklist_file(1:LEN_TRIM(blacklist_file))

  INQUIRE(file=path,EXIST=lb_exist)
  IF (.NOT.lb_exist) THEN
    WRITE(*,*) "LHN: blacklist file (DX) does not exist: ",path
    WRITE(nulhn,*) "blacklist file (DX) does not exist: ",path
    RETURN
  ENDIF

  PRINT *,'LHN: open blacklist file (DX): ',path

  CALL open_file( nulhnblackdx, path, 'r  ', yform_read, 0, 0, 1, .FALSE., idbg_level, yerrmsg, ierr )
  IF (nulhnblackdx == 0) RETURN ! there is no data file open yet
  IF (ierr /= 0) THEN
     WRITE(nulhn, *) 'Error in opening the blacklist file (DX): ',blacklist_file
     RETURN
  ENDIF

! Read data for one data record
!-------------------------------------------------------------------------------

! read time information
! read one record of the GRIB file

  maxlenc=INT(lfd*iwlength, intgribc)

! Initializations for the grib library
! (input product definition section)

  idims=0_iintegers
  idims(1) = npds
  idims(2) = ngds
  idims(3) = nbms
  idims(4) = nbds
  idims(5) = -1
  idims(6) = ndsup
  idims(7) = lds
  idims(8) = lfd

  readloop: DO !until no more record

     SELECT CASE (yform_read)
     CASE ('grb1')
#ifdef GRIBDWD
       CALL cuegin(nulhnblackdx, maxlenc, iblock_rad, ilenc, ierrc)

       ierr=INT(ierrc, iintegers)
       IF (ierr /= 0) THEN
          WRITE(nulhn, *) 'Error in reading the first record of the GRIB blacklist file (DX).'
          IF (ierr == -2) THEN
             WRITE(nulhn, *)'Field length of the GRIB field not correct.'
          ENDIF
       ENDIF

       ilen=INT(ilenc, iintegers)
#endif
     CASE ('apix')
#ifdef GRIBAPI
       CALL grib_new_from_file (nulhnblackdx, igribid, ierrc)

       IF ( ierrc /= GRIB_SUCCESS ) THEN
         IF (ierrc == GRIB_END_OF_FILE) THEN
           !EOF is reached 
           ilen = 0
         ELSE
           ierr = 2
         ENDIF
       ELSE
         ! get edition number
         CALL grib_get (igribid, 'editionNumber', ieditionnr, ierrc)
         IF ( ierrc /= GRIB_SUCCESS ) THEN
           WRITE (nulhn, *) 'Error in reading edition number with grib_get'
           ierr = 3
         ENDIF

         ! get total length of message
         CALL grib_get (igribid, 'totalLength', ilen, ierrc)
         IF ( ierrc /= GRIB_SUCCESS ) THEN
           WRITE (nulhn, *) 'Error in reading totalLength with grib_get'
           ierr = 3
         ENDIF
       ENDIF
#endif
     END SELECT

     IF (ilen == 0_iintegers) THEN
        WRITE(nulhn, *)'end of blacklist file reached'
#ifdef GRIBAPI
        IF (yform_read == 'apix') THEN
          ! The handle can be released now
          CALL grib_release(igribid)
        ENDIF    
#endif
        EXIT readloop
     ENDIF

     SELECT CASE (yform_read)
     CASE ('grb1')
#ifdef GRIBDWD
       CALL grbin1(idwdednr, undefgrib, ndims, idims, iblock_rad, ibmap_rad, &
                   ipds, igds, ibms, ibds, dsup_rad, ds_rad, ierr)
       IF (ierr /= 0) THEN
          WRITE(nulhn, *)'Error in decoding the first record and ipds of the GRIB blacklist file (DX)'
       ENDIF

       ! get element number
       iee = ipds(7)
#endif
     CASE ('apix')
#ifdef GRIBAPI
       IF     (ieditionnr == 1) THEN
         ! Get element number
         CALL grib_get (igribid, 'indicatorOfParameter', iee, ierr)
       ELSEIF (ieditionnr == 2) THEN
!US      this is only a very crude first guess: this should be refined soon
         CALL grib_get (igribid, 'shortName',            yzshortname, ierr)
         IF     (TRIM(yzshortname) == 'RAD_BL') THEN
           iee = 63
         ENDIF
       ENDIF

       ! Get size of data and data itself
       ! Set missing_value before
       CALL grib_set (igribid, 'missingValue', missval)
       CALL grib_get_size (igribid, 'values', igriblen, ierr)
       IF (igriblen > lds) THEN
         WRITE (nulhn,*) ' *** ERROR: size of message is too big for allocated field: ', igriblen, lds
         ierr = 6
       ENDIF

         CALL grib_get (igribid, 'values', ds_rad, ierr)
#endif
     END SELECT

     IF (iee == 63_iintegers ) THEN
       DO j=1,je_tot
         DO i=1,ie_tot
           ind=(j-1)*ie_tot+i
           blacklist_tot(i,j)=REAL(ds_rad(ind), wp)     ! ds_rad is irealgrib
         ENDDO
       ENDDO
#ifdef GRIBAPI
       IF (yform_read == 'apix') THEN
         ! The handle can be released now
         CALL grib_release(igribid)
       ENDIF    
#endif

       EXIT readloop
     ENDIF

#ifdef GRIBAPI
    IF (yform_read == 'apix') THEN
      ! The handle can be released now
      CALL grib_release(igribid)
    ENDIF    
#endif

  ENDDO readloop

  CALL close_file( nulhnblackdx, yform_read, 0, 0, 1, .FALSE., idbg_level, yerrmsg, ierr)
  IF (ierr /= 0) THEN
    WRITE(nulhn,*) 'Something wrong by closing blacklist file ('//yform_read//')'
  ENDIF

!-------------------------------------------------------------------------------
! End of subroutine
!-------------------------------------------------------------------------------

END SUBROUTINE lhn_blacklist_dx_open_and_read

SUBROUTINE lhn_dx_height_open_and_read (nulhnh,height_file, height_tot, ierr)
!-------------------------------------------------------------------------------
!
! Description:
!   This procedure in the module "lheat_nudge" is called by "lhn_obs_prep" and
!   opens and reads the heights for the DX data.
!
! Method:
!   Call of GRIB I/O routines of DWDLIB.
!
! Input files:
!   GRIB file.
!  
!-------------------------------------------------------------------------------

! Subroutine scalars (intent: out):
!----------------------------------

  INTEGER (KIND=iintegers), INTENT(OUT)            ::       &
    nulhnh         ,& ! unit number for radardata GRIB file
    ierr             ! status/error code

  CHARACTER (LEN=*), INTENT(IN)                   ::       &
    height_file


! Subroutine scalars (intent: out):
!----------------------------------

! Subroutine arrays (intent: out):
!---------------------------------

  REAL (KIND=wp),     INTENT(OUT)                   ::       &
    height_tot(ie_tot,je_tot,nradar)        ! field for decoded data

! Local parameters, scalars, arrays :
!-------------------------------------------------------------------------------
! Local scalars :
!----------------
  INTEGER (KIND=iintegers)             ::       &
    i,j             ,& ! loop counter
    ind,n           ,& ! index
    igribid         ,& ! 
    ieditionnr      ,& ! 
    igriblen        ,& ! 
    ilen               !

  CHARACTER (LEN=100)    ::  path
  LOGICAL :: lh_exist

!- End of header
!===============================================================================

!-------------------------------------------------------------------------------
! Begin of subroutine
!-------------------------------------------------------------------------------
  yroutine = 'lhn_dx_height_open_and_read'

! Open input file and read the header information
!-------------------------------------------------------------------------------
! open input file
  nulhnh=0_iintegers

! open GRIB file

  path=radar_in(1:LEN_TRIM(radar_in))//height_file(1:LEN_TRIM(height_file))

  INQUIRE(file=path,EXIST=lh_exist)
  IF (.NOT.lh_exist) THEN
    WRITE(*,*) "LHN: radar beam height file (DX) does not exist: ",path
    WRITE(nulhn,*) "radar beam height file (DX) does not exist: ",path
    RETURN
  ENDIF

  WRITE(*,*)'LHN: open height file (DX): ',path
  WRITE(nulhn,*)'LHN: open height file (DX): ',path

  CALL open_file( nulhnh, path, 'r  ', yform_read, 0, 0, 1, .FALSE., idbg_level, yerrmsg, ierr )

  IF (nulhnh == 0) RETURN !there is no data file open yet
  IF (ierr /= 0) THEN
     WRITE(nulhn, *)'Error in opening the height file (DX): ',height_file
     RETURN
  ENDIF

! Read data for one data record
!-------------------------------------------------------------------------------

! read time information
! read one record of the GRIB file

  maxlenc=INT(lfd*iwlength, intgribc) 
! Initializations for the grib library
! (input product definition section)

     idims=0_iintegers
     idims(1) = npds
     idims(2) = ngds
     idims(3) = nbms
     idims(4) = nbds
     idims(5) = -1
     idims(6) = ndsup
     idims(7) = lds
     idims(8) = lfd

  height_tot=-999.9_wp
  n=0

  readloop: DO !until no more record

     SELECT CASE (yform_read)
     CASE ('grb1')
#ifdef GRIBDWD
       CALL cuegin(nulhnh, maxlenc, iblock_rad, ilenc, ierrc)

       ierr=INT(ierrc, iintegers)
       IF (ierr /= 0) THEN
          WRITE(nulhn, *) 'Error in reading the first record of the GRIB height file (DX).'
          IF (ierr == -2) THEN
             WRITE(nulhn, *)'Field length of the GRIB field not correct.'
          ENDIF
       ENDIF

       ilen=INT(ilenc, iintegers)
#endif
     CASE ('apix')
#ifdef GRIBAPI
       CALL grib_new_from_file (nulhnh, igribid, ierrc)

       IF ( ierrc /= GRIB_SUCCESS ) THEN
         IF (ierrc == GRIB_END_OF_FILE) THEN
           !EOF is reached 
           ilen = 0
         ELSE
           ierr = 2
         ENDIF
       ELSE
         ! get edition number
         CALL grib_get (igribid, 'editionNumber', ieditionnr, ierrc)
         IF ( ierrc /= GRIB_SUCCESS ) THEN
           WRITE (nulhn, *) 'Error in reading edition number with grib_get'
           ierr = 3
         ENDIF

         ! get total length of message
         CALL grib_get (igribid, 'totalLength', ilen, ierrc)
         IF ( ierrc /= GRIB_SUCCESS ) THEN
           WRITE (nulhn, *) 'Error in reading totalLength with grib_get'
           ierr = 3
         ENDIF
       ENDIF
#endif
     END SELECT

     IF (ilen == 0_iintegers) THEN
        WRITE(nulhn, *)'end of height file reached'
#ifdef GRIBAPI
        IF (yform_read == 'apix') THEN
          ! The handle can be released now
          CALL grib_release(igribid)
        ENDIF    
#endif
        EXIT readloop
     ENDIF

     SELECT CASE (yform_read)
     CASE ('grb1')
#ifdef GRIBDWD
       ds_rad(:)=undefgrib
       CALL grbin1(idwdednr, undefgrib, ndims, idims, iblock_rad, ibmap_rad, &
                   ipds, igds, ibms, ibds, dsup_rad, ds_rad, ierr)
       IF (ierr /= 0) THEN
          WRITE(nulhn, *)'Error in decoding the first record and ipds of the GRIB height file (DX)'
       ENDIF
#endif
     CASE ('apix')
#ifdef GRIBAPI
       ! Get size of data and data itself
       CALL grib_get_size (igribid, 'values', igriblen, ierr)
       IF (igriblen > lds) THEN
         WRITE (nulhn,*) ' *** ERROR: size of message is too big for allocated field: ', igriblen, lds
         ierr = 6
       ENDIF

       CALL grib_get (igribid, 'values', ds_rad, ierr)
#endif
     END SELECT


     n=n+1
     IF ( n > nradar) THEN
        yerrmsg = 'ERROR  *** number of radar height records exceeds nradar (namelist NUDGING)!'
        CALL model_abort (my_cart_id,6411,yerrmsg,yroutine)
     ENDIF

     DO j=1,je_tot
       DO i=1,ie_tot
           ind=(j-1)*ie_tot+i
           height_tot(i,j,n)=REAL(ds_rad(ind), wp)     ! ds_rad is irealgrib
       ENDDO
     ENDDO

#ifdef GRIBAPI
     IF (yform_read == 'apix') THEN
       ! The handle can be released now
       CALL grib_release(igribid)
     ENDIF    
#endif

  ENDDO readloop

  IF ( n /= nradar ) THEN
     WRITE(nulhn,*) 'CAUTION: number of radar height records does not match number of radars: ',n,nradar
     WRITE(nulhn,*) '         bright band detection can not performed precisely!'
     WRITE(*,*) 'LHN: CAUTION: number of radar height records does not match number of radars: ',n,nradar
     WRITE(*,*) '         bright band detection can not performed precisely!'
  ENDIF

  CALL close_file( nulhnh, yform_read, 0, 0, 1, .FALSE., idbg_level, yerrmsg, ierr)
  IF (ierr /= 0) THEN
    WRITE(nulhn,*) 'Something wrong by closing height file with ('//yform_read//')'
  ENDIF

!-------------------------------------------------------------------------------
! End of subroutine
!-------------------------------------------------------------------------------

END SUBROUTINE lhn_dx_height_open_and_read


SUBROUTINE detect_bright_band(sumrad)
!-------------------------------------------------------------------------------
!
! Description:
!   This procedure in the module "lheat_nudge" is called by "lhn_obs_prep" and
!   detects all grid points which are possibly influenced by bright band effects
!  
!-------------------------------------------------------------------------------

  REAL (KIND=wp),     INTENT(IN)                   ::       &
   sumrad(ie,je)

  INTEGER (KIND=iintegers) :: &
    nbright

  INTEGER (KIND=iintegers) :: &
    i,j,nh,anzheight(ie,je)

  REAL  (KIND=wp)                    ::       &
    hzero_mod(ie,je),height_interval(ie,je),minheight(ie,je),maxheight(ie,je)

!- End of header
!===============================================================================
!-------------------------------------------------------------------------------
! Begin of subroutine
!-------------------------------------------------------------------------------
  yroutine = 'detect_bright_band'

     !-------------------------------------------------------------------------------
     ! Evaluate radar height with resprect to temperature
     !-------------------------------------------------------------------------------

     anzheight         = 0_iintegers
     nbright           = 0_iintegers

     CALL calhzero( hzero_mod(:,:), t(:,:,:,nnow_nold), hhl, hhl_prof, &
                    ie, je, ke, t0_melt)
     minheight(:,:) = 9999.9_wp
     maxheight(:,:) = -9999.9_wp
     height_interval(:,:) = 0.0_wp
     brightband(:,:) = 0.0_wp
     DO i=istart,iend
        DO j=istart,jend
                   
           IF(NINT(blacklist(i,j)) /= 1_iintegers) THEN
              DO nh=1,nradar
                 minheight(i,j)=min(dxheight(i,j,nh),minheight(i,j))
                 maxheight(i,j)=max(dxheight(i,j,nh),maxheight(i,j))
              ENDDO
              height_interval(i,j)=(maxheight(i,j)-minheight(i,j))/1000.0_wp
              DO nh=1,nradar
                 IF  ( dxheight(i,j,nh) > 0.0_wp                   &
               .AND.  (hzero_mod(i,j)-dxheight(i,j,nh)) >= -100_wp &
               .AND.  (hzero_mod(i,j)-dxheight(i,j,nh)) <= 1000_wp &
!               .AND.   sumrad(i,j) > 8.5_wp) THEN
               .AND.   sumrad(i,j) > 1.0_wp) THEN
                       brightband(i,j)=1.0_wp
                 ENDIF
              ENDDO
           ENDIF
           IF (brightband(i,j) > 0.0_wp) nbright=nbright+1
        ENDDO
     ENDDO

     ttm_cv(:,:,ke-33) = height_interval(:,:)
     ttm_cv(:,:,ke-34) = minheight(:,:)
     ttm_cv(:,:,ke-35) = maxheight(:,:)

     CALL global_values (nbright,1, 'SUM',imp_integers,icomm_cart, 0,yerrmsg, izerror)
     IF (my_cart_id == 0) &
     WRITE(nulhn, *)' n of points which are bright band    : numbright = ',nbright

!-------------------------------------------------------------------------------
! End of subroutine
!-------------------------------------------------------------------------------

END SUBROUTINE detect_bright_band

!===============================================================================
!+ Module procedure in "lheat_nudge" merging model and observed precipitation
!-------------------------------------------------------------------------------
 
SUBROUTINE lhn_pr_ana(ntreat)

!-------------------------------------------------------------------------------
!
! Description:
!   This procedure of the module "lheat_nudge" analyzes the precipitation
!   rate (which is later used for determining the total latent heat increment
!   to be applied).
!
!   Input arrays : pr_mod,pr_obs
!                  wobs_space,wobs_time
!   Output arrays: pr_ana
!                  i_treat,j_treat 
!
! Method:
!   A simple weighting of the observations and the model first guess value
!   is done at each grid point taking spatial and temporal weighting factors
!   for the observation into account. The weighting factors have been 
!   determined earlier (in procedure lhn_obs_proc).
!
!-------------------------------------------------------------------------------

! Subroutine arguments :
!-------------------------------------------------------------------------------
! Scalar arguments, intent(out) :
!-------------------------------
  INTEGER   (KIND=iintegers), INTENT(OUT)     ::       &
    ntreat             ! number of grid points to be treated by lhn
 
! Local parameters, scalars, arrays :
!-------------------------------------------------------------------------------
! Local scalars:
!---------------

  INTEGER   (KIND=iintegers)       ::       &
    i,j           ! loop indeces

  REAL (KIND=wp)                   ::           &
    w           ! final (normalized) weight for an observation

!- End of header
!===============================================================================
!-------------------------------------------------------------------------------
! Begin of subroutine
!-------------------------------------------------------------------------------

! analyze precipitation by weighting observation and model values

  ntreat = 0_iintegers
  DO   j=jstart,jend
    DO i=istart,iend
       w = wobs_time(i,j) * wobs_space(i,j)
       IF ( w <= 0.0_wp                                                  &
         .OR. ((pr_obs(i,j) < thres_lhn) .AND. (pr_mod(i,j) < thres_lhn) )   &
         .OR. NINT(blacklist(i,j)) == 1_iintegers                            &
         .OR. NINT(brightband(i,j)) == 1_iintegers                           &
          ) THEN
          pr_ana(i,j) = pr_mod(i,j)
       ELSE
          IF ( w < 0.0_wp .OR. w > 1.0_wp ) &
             PRINT *,'WARNING : lhn_pr_ana w unvalid : ',w,i,j,my_cart_id
          pr_ana(i,j) = w * pr_obs(i,j) + (1.0_wp-w) * pr_mod(i,j)
          ntreat = ntreat +1
          i_treat(ntreat) = i
          j_treat(ntreat) = j
          treat_diag(i,j)=0.0_wp
       ENDIF
    ENDDO
  ENDDO

!-------------------------------------------------------------------------------
! End of subroutine 
!-------------------------------------------------------------------------------

END SUBROUTINE lhn_pr_ana

!===============================================================================
!+ Module procedure in "lheat_nudge" determining T - increments due to LHN
!-------------------------------------------------------------------------------
 
SUBROUTINE lhn_t_inc (ntreat)

!-------------------------------------------------------------------------------
!
! Description:
!
!   Namelist parameters used: lhn_search,rlhn_search
!                             fac_lhn_up,fac_lhn_down
!                             lhn_relax,nlhn_relax
!                             lhn_limit,abs_lhn_lim
!                             lhn_filt,lhn_diag
!                             lhn_incloud
!                             
!   Input arrays : tt_lheat, 
!                  pr_ana,pr_mod
!                  i_treat,j_treat
!   Output arrays: tinc_lhn
!
! Method:
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Subroutine arguments :
!-------------------------------------------------------------------------------
! Scalar arguments, intent(in) :
!-------------------------------
  INTEGER   (KIND=iintegers), INTENT(IN)     ::       &
    ntreat             ! number of grid points to be treated by lhn

!-------------------------------------------------------------------------------
! Local parameters, scalars, arrays :
!-------------------------------------------------------------------------------
! Local scalars:
!---------------
  LOGICAL, SAVE                     ::           &
    lfirst=.TRUE.      ! flag for first call of this subroutine

  LOGICAL                           ::           &
    lfound(ntreat)             ! flag for a succesful search for a nearby profile

  INTEGER   (KIND=iintegers)                 ::       &
    nprof_max       ,& ! maximal number of profiles needed from adjacent nodes
    nlev            ,& ! number of levels of profiles to be treated
    n_fetch         ,& ! number of profiles needed from adjacent nodes
    i_near          ,& ! i-index for nearby profile
    j_near          ,& ! j-index for nearby profile
    i_near_tot      ,& ! i-index (in total model grid) for nearby profile
    j_near_tot      ,& ! j-index (in total model grid) for nearby profile
    i,j,k           ,& ! indeces for position in local subdomain
    i_tot,j_tot     ,& ! indeces for position in total model grid
    n,ipoint        ,& ! loop indeces
    n_local         ,& ! diagnostic : number of points with local profiles
    n_search        ,& ! diagnostic : number of points with search calls
    n_success       ,& ! diagnostic : number of points with successful search
    n_clim          ,& ! diagnostic : number of points with climatic profile
    n_up_lim        ,& ! diagnostic : number of points with limited upsclaing
    n_down_lim      ,& ! diagnostic : number of points with limited downscaling
    n_ex_lim_p      ,& ! diagnostic : number of pts with inc. above limit
    n_ex_lim_n      ,& ! diagnostic : number of pts with inc. below (neg) limit
    ierror             ! error flag
    
! for check only, delete later
  INTEGER   (KIND=iintegers)                 ::       &
    icheck,jcheck,nodecheck,bottom

  REAL (KIND=wp)                             ::       &
    eps=0.2_wp/3600.0_wp ,& ! limit used for profile filterering (flags/
                              ! eliminates values below 0.2 deg heating/h)
    pr_quot         ,& ! ratio of analyzed (observed) to model precipitation
    pr_clim         ,& ! precip corresponding to climatolgical heating profile
    abs_lim_neg     ,& ! negative absolut limit for increments (= -abs_lhn_lim)
    slope           ,&
    prmax           ,&
    prmax_th

  REAL (KIND=wp)                             ::       &
    pr_quot_clim    ,&      ! ratio of analyzed (observed) to clim. precipitation rate
    pr_clim_ratio(ie,je) ,& ! precip corresponding to climatolgical heating profile
    pr_clim_tot(ie_tot,je_tot)      ! precip corresponding to climatolgical heating profile

! Local arrays:
!--------------

! choice of profiles
  INTEGER   (KIND=iintegers)                              ::       &
    i_prof(ntreat)   ,& ! i index of profile to scale (coord. at node concerned)
    j_prof(ntreat)   ,& ! j index of profile to scale (coord. at node concerned)
    node_prof(ntreat)   !  node number at which profile to scale is found

  INTEGER   (KIND=iintegers), ALLOCATABLE                 ::       &
    i_fetch(:)    ,& ! i-indices of profiles to get from  other nodes
    j_fetch(:)    ,& ! j-indices of profiles to get from  other nodes
    node_fetch(:) ,& ! node from which a profile has to be transferred
    num_fetch(:)     ! treated point belonging to eached fetched profile

  LOGICAL                             ::       &
    lfirstsearch  ,& ! .true. for the first profile search at a timestep
    lelim         ,& ! =1 for elimination of isolated peaks in profile filter
    lsmooth          ! =1 for smoothing in profile filter

  INTEGER   (KIND=iintegers)                              ::       &
    nradius(rlhn_search+2)    ,& ! info on successful search at different radii
    isearch(4*rlhn_search*(rlhn_search+1),2)    ! template holding i,j-positions
                                ! of searched points relative to treated point

  INTEGER   (KIND=iintegers)                              ::       &
    nelimosc, nelimiso, nsmooth ,& ! diagnostic variables for filter_prof
    nelimosc_proc, nelimiso_proc, nsmooth_proc, &
    n_windcor, n_windcor0

  REAL (KIND=wp),     ALLOCATABLE                         ::       &
    tt_lh_prof(:,:)  ! array for heating profiles received from other nodes

  REAL (KIND=wp)                             ::       &
    pr_mod_tot(ie_tot,je_tot)   ,& ! complete field of total model precipitation
    hsurf_tot(ie_tot,je_tot)    ,& ! complete field of surface height
    scale_fac(ntreat)           ,& ! scaling factor determined for each profile
    prof_filt(ke)               ,& ! array for filtered vertical heating profile
    tt_lhn(ke)                  ,& ! scaled total t-change due to latant heat release
    tt_clim_35_orig(35)         ,& ! climatological heating profile for 35 levels
    halflevel_35_orig(36)       ,& !
    level_new(ke)               ,&
    level_35_orig(35)

  REAL (KIND=wp)        ::       &
    tt_clim(ke)                    ! climatological heating profile

  REAL (KIND=wp)        ::       &
    umean, vmean, zvb, zvb_llim, zvb_ulim, w950, w850, w700, wind_corr, &
    dpsum


  DATA tt_clim_35_orig / 0.00000000_wp,0.00000000_wp,0.00000000_wp,0.00000000_wp,0.00000000_wp, & !K/s
                         0.00000000_wp,0.00000000_wp,0.00000000_wp,0.00001111_wp,0.00001944_wp, &
                         0.00006389_wp,0.00013889_wp,0.00033333_wp,0.00047500_wp,0.00062222_wp, &
                         0.00079167_wp,0.00105000_wp,0.00125278_wp,0.00144444_wp,0.00150000_wp, &
                         0.00147500_wp,0.00133889_wp,0.00103056_wp,0.00081111_wp,0.00064167_wp, &
                         0.00054444_wp,0.00045556_wp,0.00035556_wp,0.00029444_wp,0.00020833_wp, &
                         0.00011944_wp,0.00004167_wp,0.00001389_wp,0.00000000_wp,0.00000000_wp/

! corresponds to the following level_heights (half levels!!) (predefined in src_artifdata.f90)
  DATA halflevel_35_orig / 23580.4414_wp, 21773.3691_wp, 18455.5488_wp, 16559.7266_wp, 14970.4609_wp, &
                           13623.1367_wp, 12477.7646_wp, 11478.0508_wp, 10561.6865_wp,  9712.8486_wp, &
                            8919.9834_wp,  8154.0073_wp,  7431.7842_wp,  6764.6855_wp,  6144.3555_wp, &
                            5564.2437_wp,  5019.1157_wp,  4504.7178_wp,  4017.5457_wp,  3554.6763_wp, &
                            3125.2937_wp,  2725.8110_wp,  2353.2324_wp,  2015.3969_wp,  1719.1958_wp, &
                            1461.1770_wp,  1238.5348_wp,  1039.6235_wp,   853.8964_wp,   680.6733_wp, &
                             519.3526_wp,   378.1536_wp,   256.2476_wp,   152.9479_wp,    67.6839_wp, &
                               0.0000_wp /

  DATA pr_clim/0.0008333_wp/        ! corresponds to precip of 3 mm/h

!- End of header
!===============================================================================
!-------------------------------------------------------------------------------
! Begin of subroutine
!-------------------------------------------------------------------------------
  yroutine='lhn_t_inc'
  ierror  = 0_iintegers
!-------------------------------------------------------------------------------
! Section 0 : Preliminaries : 
!             - determine (approx.) max number of profiles from neighbour nodes
!             Initialize fields
!-------------------------------------------------------------------------------

! calculate the height of the levels (middle of layers), location, 
! where lh-values are defined
  
  DO i=1,ke
     level_new(i)=(vcoord%vert_coord(i)+vcoord%vert_coord(i+1))*0.5_wp
  ENDDO
  DO i=1,35
     level_35_orig(i)=(halflevel_35_orig(i)+halflevel_35_orig(i+1))*0.5_wp
  ENDDO

! calculate new profile
  tt_clim(1)=0.0_wp
  DO i=2,ke
! between which levels of the 35 layer model lies the level of the current model
     DO j=1,34
        bottom=1
        IF ( (level_new(i) < level_35_orig(j)) .AND. &
             (level_new(i) > level_35_orig(j+1)) ) THEN
           bottom=j
           EXIT
        ENDIF
     ENDDO


! calculate the slope between the 35 layer model levels
     slope=(tt_clim_35_orig(bottom+1)-tt_clim_35_orig(bottom))/ &
           (level_35_orig(bottom+1)-level_35_orig(bottom))

! calculate the tt_clim(i) for the current model layers

     tt_clim(i)=tt_clim_35_orig(bottom)+slope*(level_new(i)-level_35_orig(bottom))
  ENDDO

!ks>> tt_clim_int as criteria for search_profile
  DO j=jstart,jend
  DO i=istart,iend
     dpsum=0.0_wp
     tt_clim_int(i,j) = 0.0_wp
     tt_lheat_int(i,j) = 0.0_wp
     pr_clim_ratio(i,j) = 0.0_wp
  DO k=1,ke
     tt_clim_int(i,j)  = tt_clim_int(i,j) + tt_clim(k)  &
                       * dp0(i,j,k)
     tt_lheat_int(i,j) = tt_lheat_int(i,j) + tt_lheat(i,j,k,nnow_nold)  &
                       * dp0(i,j,k)
     dpsum=dpsum+dp0(i,j,k)
  ENDDO
     tt_clim_int(i,j) = tt_clim_int(i,j) / dpsum
     tt_lheat_int(i,j) = tt_lheat_int(i,j) / dpsum
     IF (tt_lheat_int(i,j) > 0.0_wp) pr_clim_ratio(i,j)= pr_clim * tt_lheat_int(i,j) / tt_clim_int(i,j)
  ENDDO
  ENDDO
!ks<<

! set temperature increments to be determined to zero
  tinc_lhn = 0.0_wp

! initializations of flag, arrays for profile search
  lfirstsearch = .TRUE.
  nprof_max = rlhn_search * 2*(ie+je)
  isearch = 0
  nradius = 0

! flags for profile filter
  lelim=.TRUE.
  lsmooth=.TRUE.

! initializations of counters
  n_local    = 0
  n_search   = 0
  n_success  = 0
  n_clim     = 0
  n_up_lim   = 0
  n_down_lim = 0
  n_windcor  = 0
  n_windcor0 = 0

! other initializations (assume that profiles are at local i,j,node first)
  scale_fac = 1._wp
  i_prof = -1_iintegers
  j_prof = -1_iintegers
  node_prof = -1_iintegers

  nlev = ke

!-------------------------------------------------------------------------------
! Section 1 : Assemble the complete field of model precipitation and topo at 
!             each PE if profile search is requested (lhn_search=.true.) and if
!             the search radius extends to the domain of adjacent processors.
!-------------------------------------------------------------------------------
  IF (num_compute > 1) THEN
     IF (lhn_search) THEN

        IF (ltime) THEN
          CALL get_timings (i_lhn_t_inc, ntstep, dt, izerror)
          IF (ltime_barrier) THEN
            CALL comm_barrier (icomm_cart, izerror, yerrmsg)
            CALL get_timings (i_barrier_waiting_lhn, ntstep, dt, izerror)
          ENDIF
        ENDIF

        CALL gather_field (pr_mod, ie, je, pr_mod_tot, ie_tot, je_tot, -1, ierror)
        CALL gather_field (hsurf,  ie, je, hsurf_tot,  ie_tot, je_tot, -1, ierror)
        CALL gather_field (pr_clim_ratio,  ie, je, pr_clim_tot,  ie_tot, je_tot, -1, ierror)

        IF (ltime) THEN
           CALL get_timings (i_communications_lhn, ntstep, dt, izerror)
        ENDIF

     ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! Section 2 : Select the latent heating profile to scale and set scaling factor:
!             (only points with pr_ana != pr_mod are treated)
!             - determine if the local or a nearby heating profile is taken
!             - search for a nearby heating profile (if lhn_search)
!             - retain the indices of the selected profile and assign them to a
!               special list if it belongs to a subdomain held at another node
!             - set the scaling factor (depending on the ratio of analyzed to
!               model precipitation and respecting the scaling limits)
!-------------------------------------------------------------------------------

  DO ipoint=1,ntreat
    i=i_treat(ipoint)
    j=j_treat(ipoint)
    pr_quot = (pr_ana(i,j)+repsilon) / (pr_mod(i,j)+repsilon)
    pr_quot_clim = pr_ana(i,j) / pr_clim

   IF(lhn_logscale) THEN
     pr_quot = 1._wp + LOG(pr_quot)
   ENDIF

! local model precip is within range [ fac_lhn_down*pr_ana , fac_lhn_up*pr_ana ]
! -> scaling of a local profile (upscaling or downscaling)
    IF ( pr_quot >= fac_lhn_down .AND. pr_quot <= fac_lhn_up ) THEN
       scale_fac(ipoint) = pr_quot
       i_prof(ipoint) = i
       j_prof(ipoint) = j
       node_prof(ipoint) = my_cart_id
       n_local = n_local + 1
       treat_diag(i,j)=1.0_wp

! local model precip is too large -> limited downscaling of local profile
    ELSEIF ( pr_quot < fac_lhn_down ) THEN
       scale_fac(ipoint) = fac_lhn_down
       i_prof(ipoint) = i
       j_prof(ipoint) = j
       node_prof(ipoint) = my_cart_id
       n_down_lim = n_down_lim + 1
       treat_diag(i,j)=2.0_wp

! local model precip is too small and search for nearby profile allowed -> search
! for a nearby profile to scale
    ELSEIF ( pr_quot > fac_lhn_search .AND. lhn_search ) THEN

       IF (ltime) THEN
          CALL get_timings (i_lhn_t_inc, ntstep, dt, izerror)
       ENDIF

!ks: Check if there isn't already enough tt_lheat within the profile

!!sk20060818       IF (tt_lheat_int(i,j) < pr_quot_clim*tt_clim_int) THEN 
          IF (num_compute > 1) THEN
             i_tot = i_global(i)
             j_tot = j_global(j)
          ELSE
             i_tot = i
             j_tot = j
          ENDIF
          CALL lhn_prof_search (i_tot,j_tot,pr_ana(i,j),hsurf(i,j)              &
                               ,pr_mod_tot,hsurf_tot,pr_clim_tot                &
                               ,ie_tot,je_tot,rlhn_search,lfirstsearch,isearch  &
                               ,lfound(ipoint),i_near_tot,j_near_tot,nradius)

          n_search = n_search + 1
          ! if a profile is found: localize the node and the coordinates at that node
          ! where the nearby profile can be found
          IF (lfound(ipoint)) THEN
             CALL ij_local (i_near_tot,j_near_tot,icheck,jcheck,nodecheck,izerror)
             i_prof(ipoint) = icheck
             j_prof(ipoint) = jcheck
             node_prof(ipoint) = nodecheck

             prmax    = pr_mod(i,j)*fac_lhn_up
             prmax_th = thres_lhn*fac_lhn_up
             prmax_th = MIN(prmax_th,pr_ana(i,j))
             prmax    = MAX(prmax,prmax_th)

             scale_fac(ipoint) = 1._wp + (prmax-pr_mod(i,j))/pr_mod_tot(i_near_tot,j_near_tot)
             IF (lhn_logscale) scale_fac(ipoint) = 1._wp + LOG(scale_fac(ipoint))
             scale_fac(ipoint) = MIN(scale_fac(ipoint),fac_lhn_up)
            
             treat_diag(i,j)=3.0_wp

             IF (ltime) THEN
                CALL get_timings (i_lhn_search, ntstep, dt, izerror)
             ENDIF

          ELSE IF (pr_clim >= pr_mod(i,j)) THEN
             ! no appropriate profile is found -> scaling of climatological profile

             prmax    = pr_mod(i,j)*fac_lhn_up
             prmax_th = thres_lhn*fac_lhn_up*3 !ks-fuddle_faktor (-;
             prmax_th = MIN(prmax_th,pr_ana(i,j))
             prmax    = MAX(prmax,prmax_th)

             scale_fac(ipoint) = 1._wp + (prmax-pr_mod(i,j))/pr_clim
             IF (lhn_logscale) scale_fac(ipoint) = 1._wp + LOG(scale_fac(ipoint))
             scale_fac(ipoint) = MIN(scale_fac(ipoint),fac_lhn_up)

             i_prof(ipoint) = 0
             j_prof(ipoint) = 0
             node_prof(ipoint) = my_cart_id
             n_clim = n_clim + 1
             treat_diag(i,j)=4.0_wp

             IF (ltime) THEN
                CALL get_timings (i_lhn_search, ntstep, dt, izerror)
             ENDIF

          ELSE
             scale_fac(ipoint) = fac_lhn_up
             i_prof(ipoint) = i
             j_prof(ipoint) = j
             node_prof(ipoint) = my_cart_id
             n_up_lim = n_up_lim + 1
             treat_diag(i,j)=5.0_wp
          ENDIF
!ks: If there is already a high amount of tt_lheat within the profile do not any further scaling
!!sk20060818       ELSE   
!!sk20060818          scale_fac(ipoint) = fac_lhn_up
!!sk20060818          i_prof(ipoint) = i
!!sk20060818          j_prof(ipoint) = j
!!sk20060818          node_prof(ipoint) = my_cart_id
!!sk20060818          n_up_lim = n_up_lim + 1
!!sk20060818          treat_diag(i,j)=6.0_wp
!!sk20060818       ENDIF
    ELSEIF ( pr_quot > fac_lhn_search ) THEN

       IF (ltime) THEN
          CALL get_timings (i_lhn_t_inc, ntstep, dt, izerror)
       ENDIF

       prmax    = pr_mod(i,j)*fac_lhn_up
       prmax_th = thres_lhn*fac_lhn_up*3 !ks-fuddle_faktor (-;
       prmax_th = MIN(prmax_th,pr_ana(i,j))
       prmax    = MAX(prmax,prmax_th)

       scale_fac(ipoint) = 1._wp + (prmax-pr_mod(i,j))/pr_clim
       IF (lhn_logscale) scale_fac(ipoint) = 1._wp + LOG(scale_fac(ipoint))
       scale_fac(ipoint) = MIN(scale_fac(ipoint),fac_lhn_up)

       i_prof(ipoint) = 0
       j_prof(ipoint) = 0
       node_prof(ipoint) = my_cart_id
       n_clim = n_clim + 1
       treat_diag(i,j)=4.0_wp

       IF (ltime) THEN
          CALL get_timings (i_lhn_search, ntstep, dt, izerror)
       ENDIF

! local model precip is too small and no search for nearby profile allowed -> limited 
! upscaling of local profile
    ELSE
       scale_fac(ipoint) = fac_lhn_up
       i_prof(ipoint) = i
       j_prof(ipoint) = j
       node_prof(ipoint) = my_cart_id
       n_up_lim = n_up_lim + 1
       treat_diag(i,j)=7.0_wp
    ENDIF
    
!ks>> 1.0 changed to fac_lhn_up
    IF (pr_quot > fac_lhn_up) THEN
       scale_fac_index(i,j)=.TRUE.
    ELSE
       scale_fac_index(i,j)=.FALSE.
    ENDIF

  ENDDO

!-------------------------------------------------------------------------------
! Section 3 : Exchange of heating profiles between different nodes (in case
!             that heating profiles from nearby gridpoints in adjacent
!             subdomains are needed).
!             Allocate arrays needed for profile exchange between processors
!-------------------------------------------------------------------------------

  IF (lhn_search .AND.  num_compute > 1) THEN
! allocation of space for arrays needed for 'order lists' of profiles
     istat = 0
     ALLOCATE ( num_fetch(nprof_max) , STAT=izlocstat )
     istat = istat + izlocstat
     num_fetch = -1
     ALLOCATE ( node_fetch(nprof_max) , STAT=izlocstat )
     istat = istat + izlocstat
     node_fetch = -1
     ALLOCATE ( i_fetch(nprof_max) , STAT=izlocstat )
     istat = istat + izlocstat
     i_fetch = -1
     ALLOCATE ( j_fetch(nprof_max) , STAT=izlocstat )
     istat = istat + izlocstat
     j_fetch = -1
     IF (istat /= 0) THEN
        yerrmsg =' ERROR  *** Allocation of arrays for profile request failed!'
        CALL model_abort (my_cart_id,6405,yerrmsg,yroutine)
     ENDIF

! determine the profiles needed from other nodes
     n_fetch = 0
     DO ipoint = 1,ntreat
        IF (node_prof(ipoint) /= my_cart_id) THEN
        ! if profile is located at another node, THEN
        ! fill fetch-arrays with info of profiles to fetch
            n_fetch = n_fetch + 1
            i_fetch(n_fetch) = i_prof(ipoint)
            j_fetch(n_fetch) = j_prof(ipoint)
            node_fetch(n_fetch) = node_prof(ipoint)
            num_fetch(n_fetch) = ipoint
        ENDIF
     ENDDO

! exchange the profiles
     ALLOCATE ( tt_lh_prof(ke,n_fetch) , STAT=izlocstat ) ; tt_lh_prof = 0.0_wp
     IF (izlocstat /= 0) THEN
        yerrmsg =' ERROR  *** Allocation of arrays for profile request failed!'
        CALL model_abort (my_cart_id,6406,yerrmsg,yroutine)
     ENDIF
     CALL exchange_profiles (nlev,n_fetch                              &
                            ,i_fetch,j_fetch,node_fetch                &
                            ,tt_lheat(1:ie,1:je,1:ke,nnow_nold)        &
                            ,1,ke                                      &
                            ,tt_lh_prof(1:ke,1:n_fetch))   
  ENDIF

!-------------------------------------------------------------------------------
! Section 4 : Scale the heating profiles with the predetermined factors
!             Insert the climatological profile if no suitable profile was found
!             -> Get the temperature increment due _only_ to lhn - correction
!             tt_lheat (= dT due to saturation adjustment has already been added to T
!             (src_leapfrog) so only the part due to LHN has to be added now:
!             tinc_lhn = tt_lheat * (scale_fac - 1)
!             Filter the profiles to exclude spurious peaks/noise
!-------------------------------------------------------------------------------

! filter and scale the local profiles and nearby profiles from this node
! or climatological profiles
  DO ipoint=1,ntreat
    IF (node_prof(ipoint) == my_cart_id) THEN
       i = i_treat(ipoint)
       j = j_treat(ipoint)
       i_near = i_prof(ipoint)
       j_near = j_prof(ipoint)
       tt_lhn = 0.0_wp

       scale_diag(i,j)=scale_fac(ipoint)

!      treat a local/nearby profile
!      i_near and j_near are defined to be zero, if a climatological profile is chosen
       IF (i_near /= 0) THEN
          tt_lhn(1:ke) = tt_lheat(i_near,j_near,1:ke,nnow_nold) &
                                    * scale_fac(ipoint)
!      insert a climatological profile
       ELSE
          tt_lhn(1:ke) = tt_clim(1:ke)*scale_fac(ipoint)
       ENDIF

!      determine the temperature correction due to lhn

       IF (i_near /= 0) THEN
          tinc_lhn(i,j,1:ke) = tt_lhn(1:ke) - tt_lheat(i_near,j_near,1:ke,nnow_nold)

!ks 13.04.05: new definition for lhn_incloud
         IF (lhn_incloud) THEN
            WHERE (tt_lheat(i_near,j_near,1:ke,nnow_nold) < 0.0_wp) &
                  tinc_lhn(i,j,1:ke) = 0.0_wp
         ENDIF

       ELSE
          tinc_lhn(i,j,1:ke) = tt_lhn(1:ke) - tt_clim(1:ke)
       ENDIF

    ENDIF
  ENDDO

  IF (lhn_search .AND. num_compute > 1) THEN
     DO n = 1,n_fetch

       ipoint = num_fetch(n)
       i      = i_treat(ipoint)
       j      = j_treat(ipoint)
       i_near = i_prof(ipoint)
       j_near = j_prof(ipoint)

       scale_diag(i,j)=scale_fac(ipoint)

!ks 12.12.05: change because of repeatability reasons
       tt_lhn(1:ke)       = tt_lh_prof(1:ke,n) * scale_fac(ipoint)
       tinc_lhn(i,j,1:ke) = tt_lhn(1:ke) - tt_lh_prof(1:ke,n)

!ks 13.04.05: new definition for lhn_incloud 
     IF (lhn_incloud) THEN
        WHERE (tinc_lhn(i,j,1:ke) < 0.0_wp) &
               tinc_lhn(i,j,1:ke) = 0.0_wp
         ENDIF
     ENDDO
  ENDIF


!-------------------------------------------------------------------------------
! Section 5 : Impose absolute limit on increments if requested (if lhn_limit)
!-------------------------------------------------------------------------------

! absolut limit to increments
  n_ex_lim_p = 0_iintegers
  n_ex_lim_n = 0_iintegers
  IF (lhn_limit) THEN
    abs_lim_neg = -1.0_wp * abs_lhn_lim
    DO ipoint=1,ntreat
      i = i_treat(ipoint)
      j = j_treat(ipoint)
      DO k=1,ke
        IF (tinc_lhn(i,j,k) > abs_lhn_lim) THEN
          tinc_lhn(i,j,k) = abs_lhn_lim
          n_ex_lim_p = n_ex_lim_p + 1
        ELSEIF (tinc_lhn(i,j,k) < abs_lim_neg) THEN
          tinc_lhn(i,j,k) = abs_lim_neg
          n_ex_lim_n = n_ex_lim_n + 1
        ENDIF
      ENDDO
    ENDDO
  ENDIF

!-------------------------------------------------------------------------------
! Section 5 : Weighting of the temperature increment with respect to
!             the mean horizontal wind within the column
!             The mean horizontal wind is defined as linear combination of
!             0.5*wind(950 hPa)+0.25*(wind(850 hPa)+wind(700 hPa))
!
!-------------------------------------------------------------------------------

  IF (lhn_wweight) THEN
     DO ipoint=1,ntreat
        i = i_treat(ipoint)
        j = j_treat(ipoint)
        w950 = 0.5_wp
        w850 = 0.25_wp
        w700 = 0.25_wp
        umean = w950 * u(i,j,klv950,nnow_nold)       &
              + w850 * u(i,j,klv850,nnow_nold)       &
              + w700 * u(i,j,klv700,nnow_nold)
        vmean = w950 * v(i,j,klv950,nnow_nold)       &
              + w850 * v(i,j,klv850,nnow_nold)       &
              + w700 * v(i,j,klv700,nnow_nold)
        zvb   = SQRT(umean * umean + vmean * vmean)
        zvb_llim = 20.0_wp
        zvb_ulim = 30.0_wp

        IF (zvb <= zvb_llim) THEN
           wind_corr=1.0_wp
        ELSE IF (zvb <= zvb_ulim) THEN
           wind_corr=1.0_wp - (1.0_wp/(zvb_ulim-zvb_llim))   &
                                  * (zvb - zvb_llim)
           n_windcor=n_windcor+1
        ELSE
           wind_corr=0.0_wp
           n_windcor0=n_windcor0+1
        ENDIF
        tinc_lhn(i,j,:) = tinc_lhn(i,j,:) * wind_corr
        windcor_diag(i,j)=wind_corr
     ENDDO
  ENDIF

!-------------------------------------------------------------------------------
! Section 6 : vertical restriction of increments
!-------------------------------------------------------------------------------

  DO ipoint=1,ntreat
     i = i_treat(ipoint)
     j = j_treat(ipoint)
     IF ( ktop_temp > 100_wp ) THEN
        DO k=1,ke
           IF (T(i,j,k,nnow_nold) >= ktop_temp) THEN
              ktop_lhn=k
              exit
           ENDIF
        ENDDO
     ENDIF
     ktop_diag(i,j)=REAL(ktop_lhn,wp)
     DO k=1,ke
        IF (k < ktop_lhn .OR. k > kbot_lhn ) &
            tinc_lhn(i,j,k)=0.0_wp
     ENDDO
  ENDDO

!-------------------------------------------------------------------------------
! Section 7 : Vertical filtering of increments (if lhn_filt)
!-------------------------------------------------------------------------------

  IF (lhn_filt) THEN
    nelimosc_proc = 0_iintegers
    nelimiso_proc = 0_iintegers
    nsmooth_proc  = 0_iintegers
    DO ipoint=1,ntreat
      i = i_treat(ipoint)
      j = j_treat(ipoint)
         CALL filter_prof (tinc_lhn(i,j,1:ke),prof_filt,1,ke,      &
                        eps,lelim,lsmooth,nelimosc,nelimiso,nsmooth)
      tinc_lhn(i,j,1:ke) = prof_filt(1:ke)
      nelimosc_proc = nelimosc_proc + nelimosc
      nelimiso_proc = nelimiso_proc + nelimiso
      nsmooth_proc  = nsmooth_proc  + nsmooth
    ENDDO

! Diagostic output from filter_prof
!-------------------------------------------------------------------------------
     IF (lhn_diag) THEN
        IF (num_compute > 1) THEN
           iexchange    = 0
           iexchange( 1)= nelimosc_proc
           iexchange( 2)= nelimiso_proc
           iexchange( 3)= nsmooth_proc

           IF (ltime) THEN
             CALL get_timings (i_lhn_t_inc, ntstep, dt, izerror)
             IF (ltime_barrier) THEN
               CALL comm_barrier (icomm_cart, izerror, yerrmsg)
               CALL get_timings (i_barrier_waiting_lhn, ntstep, dt, izerror)
             ENDIF
           ENDIF

           CALL global_values (iexchange(1:3),3, 'SUM',imp_integers,icomm_cart, &
                               0,yerrmsg, izerror)

           IF (ltime) THEN
              CALL get_timings (i_communications_lhn, ntstep, dt, izerror)
           ENDIF

        ENDIF
        IF (my_cart_id == 0) THEN
           IF (num_compute > 1) THEN
              nelimosc_proc = iexchange(1)
              nelimiso_proc = iexchange(2)
              nsmooth_proc  = iexchange(3)
           ENDIF
           WRITE(nulhn, *)' Vert. Filtering : n points eliminate oscillations: ', &
                          nelimosc_proc
           WRITE(nulhn, *)' Vert. Filtering : n points eliminate isolate peaks: ',&
                          nelimiso_proc
           WRITE(nulhn, *)' Vert. Filtering : n points smoothed : ',nsmooth_proc
        ENDIF
     ENDIF
   ENDIF

!-------------------------------------------------------------------------------
! Section 8 : Diagnostic output on lhn - increments and profile search
!             Summing information at PE0 for printout
!-------------------------------------------------------------------------------

  IF (lhn_diag) THEN

     IF (num_compute > 1) THEN
     ! get summed diagnostics at PE0 from all PE's using collect_values 
     ! (collect_values needs real vector as input)
        iexchange = 0
        iexchange( 1) = ntreat
        iexchange( 2) = n_local
        iexchange( 3) = n_up_lim
        iexchange( 4) = n_down_lim
        iexchange( 5) = n_ex_lim_p
        iexchange( 6) = n_ex_lim_n
        iexchange( 7) = n_search
        iexchange( 8) = n_clim
        iexchange( 9) = n_fetch
        iexchange(10) = n_windcor
        iexchange(11) = n_windcor0

        DO i=1,rlhn_search+2
           iexchange(11+i)=nradius(i)
        ENDDO

        IF (ltime) THEN
          CALL get_timings (i_lhn_t_inc, ntstep, dt, izerror)
          IF (ltime_barrier) THEN
            CALL comm_barrier (icomm_cart, izerror, yerrmsg)
            CALL get_timings (i_barrier_waiting_lhn, ntstep, dt, izerror)
          ENDIF
        ENDIF

        CALL global_values (iexchange,rlhn_search+2+11, 'SUM',imp_integers,  &
                            icomm_cart, 0,yerrmsg, izerror)

        IF (ltime) THEN
           CALL get_timings (i_communications_lhn, ntstep, dt, izerror)
        ENDIF

     ENDIF
       
! printout of diagnostics
   IF (my_cart_id == 0) THEN
      IF (num_compute > 1) THEN
         ntreat_tot = iexchange( 1)
         n_local    = iexchange( 2)
         n_up_lim   = iexchange( 3)
         n_down_lim = iexchange( 4)
         n_ex_lim_p = iexchange( 5)
         n_ex_lim_n = iexchange( 6)
         n_search   = iexchange( 7)
         n_clim     = iexchange( 8)
         n_fetch    = iexchange( 9)
         n_windcor  = iexchange(10)
         n_windcor0 = iexchange(11)

         DO i=1,rlhn_search+2
            nradius(i) = iexchange(11+i)
         ENDDO
      ELSE
         ntreat_tot = ntreat
      ENDIF
      WRITE(nulhn, *)
      WRITE(nulhn, *)' Diagnostics of LHN - nudging scheme, subroutine lhn_t_inc'
      IF (lfirst) THEN
         WRITE(nulhn, *)' parameters set for LHN :'
         WRITE(nulhn, *)' Profile search :           lhn_search = ',lhn_search
         WRITE(nulhn, *)' Vertical Filtering of increments  : lhn_filt   = ',lhn_filt
         WRITE(nulhn, *)' Horizontal Filtering of increments : lhn_relax  = ',lhn_relax,nlhn_relax
         WRITE(nulhn, *)' Absolute limit of incs.  :  lhn_limit = ',lhn_limit,abs_lhn_lim &
              ,' (K/second)'
         WRITE(nulhn, *)' Humidity enhancement :     lhn_hum_adj = ',lhn_hum_adj
         WRITE(nulhn, *)' Diagnostic output :        lhn_diag    = ',lhn_diag
         WRITE(nulhn, *)' Number of points treated (domain minus 2*nboundlines): ', &
              (ie_tot-2*nboundlines)*(je_tot-2*nboundlines)
      ENDIF

!      IF (ntreat_tot > 0) THEN
         WRITE(nulhn, *)'Diagnostics of LHN, lhn_t_inc, timestep : ',ntstep
         WRITE(nulhn, *)' n of points with increments         : ',ntreat_tot,ntstep
         WRITE(nulhn, *)' n of points with local profiles     : ',n_local,ntstep
         WRITE(nulhn, *)' n of points with clim. prof         : ',n_clim,ntstep
         WRITE(nulhn, *)' n of points with limited upscaling  : ',n_up_lim,ntstep
         WRITE(nulhn, *)' n of points with limited downscaling: ',n_down_lim,ntstep
         WRITE(nulhn, *)' n of points with profile search     : ',n_search,ntstep
         DO i=1,rlhn_search
            WRITE(nulhn, *)'    successful search, radius = ',i,': ',nradius(i)
            n_success = n_success + nradius(i)
         ENDDO
         WRITE(nulhn, *)' n of points with successful search  : ',n_success,ntstep
         WRITE(nulhn, *)'             with failed     search  : ',n_search-n_success
         WRITE(nulhn, *)'             with failure due to dz  : ',nradius(rlhn_search+2)
         WRITE(nulhn, *)' total number of points searched     : ',nradius(rlhn_search+1)
         WRITE(nulhn, *)' number of profiles from other nodes : ',n_fetch
         WRITE(nulhn, *)' points with imposed positive limit  : ',n_ex_lim_p
         WRITE(nulhn, *)' points with imposed negative limit  : ',n_ex_lim_n
         WRITE(nulhn, *)' points with wind weighting < 1 and > 0 : ',n_windcor
         WRITE(nulhn, *)' points with wind weighting equal 0     : ',n_windcor0
         WRITE(nulhn, *)
!      ENDIF
   ENDIF

!    DO k=ktop_lhn,kbot_lhn
!       WRITE(nulhn, *)' max t_inc at level ',k,MAXVAL(tinc_lhn(1:ie,1:je,k))
!    ENDDO
!    DO k=ktop_lhn,kbot_lhn
!       WRITE(nulhn, *)' min t_inc at level ',k,MINVAL(tinc_lhn(1:ie,1:je,k))
!    ENDDO
!   ENDIF

  ENDIF

!-------------------------------------------------------------------------------
! Section 9 : Deallocate space
!-------------------------------------------------------------------------------

! space has only been allocated if profile search is performed
  IF (lhn_search) THEN
    DEALLOCATE ( num_fetch  , STAT=izlocstat )
    DEALLOCATE ( node_fetch , STAT=izlocstat )
    DEALLOCATE ( i_fetch , STAT=izlocstat )
    DEALLOCATE ( j_fetch , STAT=izlocstat )
    IF (num_compute > 1) DEALLOCATE ( tt_lh_prof , STAT=izlocstat )
  ENDIF

! set switch for first call to false
  lfirst = .FALSE.

!-------------------------------------------------------------------------------
! End of subroutine 
!-------------------------------------------------------------------------------

END SUBROUTINE lhn_t_inc

!===============================================================================
!+ Module procedure in "lheat_nudge" determining T - increments due to LHN !NEW!
!-------------------------------------------------------------------------------
 
!===============================================================================
!+ Module procedure in "lheat_nudge" adjusting humidity to LHN T - increments
!-------------------------------------------------------------------------------
 
SUBROUTINE lhn_q_inc(ntreat)

!-------------------------------------------------------------------------------
!
! Description:
!   Subroutine which computes q-increments for the previously added temperature
!   increments resulting from LHN.
!
! Method:
!   In areas of negative T-increments, qv is adjusted so that the relative
!   humidity from before is retained unaltered.
!   In areas of positive  T-increments, qv is raised, so that saturation is
!   reached. No adjustment to qc is made, since LHN assumes that the
!   condensed water vapour has precipitated.
!
!   relhum == const ---> relhum(told) = relhum(tnew)
!
!              rdv * (relhum(told)*esat(tnew))
!       qv  = ---------------------------------------
!              p - o_m_rdv * (relhum(told)*esat(tnew))
!
!
!   Used fields : tinc_lhn, t, qv, p0, pp
!   Modified fields : qv(.,.,.,nnew)
!
!-------------------------------------------------------------------------------
 
! Scalar arguments, intent(in) :
!-------------------------------
  INTEGER   (KIND=iintegers), INTENT(IN)     ::       &
    ntreat             ! number of grid points to be treated by lhn

! Local parameters, scalars, arrays :
!------------------------------------

! Local parameters:
  REAL    (KIND=wp   ),     SAVE ::  &
   delt_minn= -3.E-6_wp ,& ! minimal T-change before applying T-adjustment
   delt_minp= 3.E-6_wp  ,& ! minimal T-change before applying T-adjustment
   f_raise  = 1._wp     ,& ! relative humidity in positive adjustment areas
   tau_nudge = 1._wp/1800._wp     ,& ! time weight for nudging of the humidity
                               ! increment (increment spread over 30. min=1800 sec.)
   fac_q_max= 2._wp        ! maximal factor allowed in change of qv

! Local scalars:
  REAL    (KIND=wp   )     ::  &
   f_esat  ,& ! Name of satement function (saturation water vapour pressure)
   f_qv    ,& ! Name of satement function (specific humidity from p,e)
   f_e     ,& ! Name of satement function (water vapour pressure from q,p)
   zt      ,& ! Dummy argument for statement functions
   ze      ,& ! ...
   zqv     ,& ! ...
   zp      ,& ! ...
   esat    ,& ! saturation water vapour pressure
   relhum     ! relative humidity

  INTEGER (KIND=iintegers)  ::  &
   i,j,k,       & ! loop indices
   nred,ninc   ,& ! number of points where humidity qv is reduced/increased
   ninc2          ! number of points where humidity qv is reduced/increased

!- End of header
!===============================================================================
!-------------------------------------------------------------------------------
! Begin of subroutine
!-------------------------------------------------------------------------------
! STATEMENT FUNCTIONS (identical to those used in meteo_utilities (satad) )
  f_esat(zt)     = b1*EXP( b2w*(zt-b3)/(zt-b4w) )      ! Magnus formula
  f_qv(ze,zp)    = rdv*ze/( zp - o_m_rdv*ze )          ! spec. hum. qv from e, p
  f_e (zqv,zp)   = MAX( epsy , zqv ) * zp / (rdv + zqv*o_m_rdv)  ! e from  qv, p

!-------------------------------------------------------------------------------
! 1. Increase / decrease specific humidity
!-------------------------------------------------------------------------------

  nred = 0
  ninc = 0
  ninc2 = 0

  zt = tau_nudge * zdt

 DO i=istart,iend
  DO j=jstart,jend

    DO   k=1,ke

     IF ( tinc_lhn(i,j,k) < delt_minn .OR. tinc_lhn(i,j,k) > delt_minp) THEN
       zp = p0(i,j,k)+pp(i,j,k,nnew)

!       !saturation pressure before temperature increment
       esat = f_esat ( t(i,j,k,nnew) - tinc_lhn(i,j,k) * zdt)

!       ! relhum before temperature increment
       relhum = f_e ( qv(i,j,k) , zp ) / esat
       relhum = MIN ( f_raise,relhum)

!       ! specific humidity after temperature increment so that relhum is unchanged
       qv(i,j,k) = f_qv ( relhum * f_esat(t(i,j,k,nnew)) , zp )

       IF ( tinc_lhn(i,j,k) > delt_minp ) THEN
        
         ninc = ninc + 1
!ks: if criteria changed
         IF ( (scale_fac_index(i,j)) .AND. (qc(i,j,k)+qi(i,j,k) <= epsy) ) THEN
!! add an additional increment to qv at gridpoints where the precipitation rate should
!! be increased and f has not reached 100% so far!
!! Attention: do not add these increments at points with a positive temperature increment
!! generally, because positive temperature increments can also occur at gridpoints where
!! the precipitation rate should be decreased (pos. temp. inc. below clouds, higher 
!! evaporation)

          ninc2= ninc2 + 1
          zqv = f_qv ( f_raise * f_esat (t(i,j,k,nnew)), zp )
          zqv = MIN ( zqv , fac_q_max * qv(i,j,k))
          zqv = zqv - qv(i,j,k)
           !! All grid points with qc+qi>0 will not be adjusted!
          qv(i,j,k) = qv(i,j,k) + zqv*zt

         ENDIF
! Diagnostics Output

       ELSE IF ( tinc_lhn(i,j,k) < delt_minn ) THEN
         nred = nred + 1
       ENDIF
     ENDIF

    ENDDO
   ENDDO
  ENDDO

!-------------------------------------------------------------------------------
! 2. Diagostic output
!-------------------------------------------------------------------------------
 IF (lhn_diag) THEN
  IF (num_compute > 1) THEN
     iexchange    = 0
     iexchange( 1)= ninc
     iexchange( 2)= nred
     iexchange( 3)= ninc2

     IF (ltime) THEN
       CALL get_timings (i_lhn_q_inc, ntstep, dt, izerror)
       IF (ltime_barrier) THEN
         CALL comm_barrier (icomm_cart, izerror, yerrmsg)
         CALL get_timings (i_barrier_waiting_lhn, ntstep, dt, izerror)
       ENDIF
     ENDIF

     CALL global_values (iexchange,3, 'SUM',imp_integers,icomm_cart, 0,yerrmsg, izerror)

     IF (ltime) THEN
        CALL get_timings (i_communications_lhn, ntstep, dt, izerror)
     ENDIF

  ENDIF
  IF (my_cart_id == 0) THEN
     IF (num_compute > 1) THEN
     ninc = iexchange(1)
     nred = iexchange(2)
     ninc2 = iexchange(3)
     ENDIF
     WRITE(nulhn, *)' Humidity adjustment : n points increased : ',ninc
     WRITE(nulhn, *)' Humidity adjustment : n points nudged to saturation: ',ninc2
     WRITE(nulhn, *)' Humidity adjustment : n points decreased : ',nred
  ENDIF
 ENDIF
!-------------------------------------------------------------------------------
! End of subroutine 
!-------------------------------------------------------------------------------

END SUBROUTINE lhn_q_inc

!===============================================================================
!+ Saturation Adjustment after changing the temperature without adjusting the moisture
!-------------------------------------------------------------------------------

SUBROUTINE lhn_satad(ntreat)

!-------------------------------------------------------------------------------
! Description:
! Changing the temperature yields a changing of the moistures fields too. If no 
! humidity adjustment is used within the LHN temperature and moisture are not 
! in thermodynamic balance any more. For consistency reasons a SATuration ADjustment 
! has to be done.
!
!-------------------------------------------------------------------------------

! Subroutine arguments
!-------------------------------------------------------------------------------
  INTEGER   (KIND=iintegers), INTENT(IN)     ::       &
    ntreat             ! number of grid points to be treated by lhn

! Local parameters, scalars, arrays :
!-------------------------------------------------------------------------------
! Local scalars:
!---------------

REAL    (KIND=wp   )     ::  &
    zte     (ie,je)  ,& ! Auxiliary fields for subr. *satad*
    zqve    (ie,je)  ,& !      |
    zqce    (ie,je)  ,& !      |
    zphfe   (ie,je)  ,& !      |
    ztstart (ie,je)  ,& !    \ | /
    zr1     (ie,je)  ,& !     \|/
    zr2     (ie,je)  ,& !      V
    zr3     (ie,je)  ,& !
    zr4     (ie,je)  ,& !
    zr5     (ie,je)  ,& !
    zr6     (ie,je)  ,& !
    zr7     (ie,je)  ,& !
    zr8     (ie,je)     !

INTEGER (KIND=iintegers) :: &
    i,j,k,kitera

!- End of header
!===============================================================================
!-------------------------------------------------------------------------------
! Begin of subroutine
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
 kitera = 1_iintegers
DO  k=1,ke
 DO i=istart,iend
  DO j=jstart,jend

   zphfe  (i,j) =  p0(i,j,k) + pp(i,j,k,nnew)
   zqve   (i,j) =  qv(i,j,k)
   zqce   (i,j) =  qc(i,j,k)
   zte    (i,j) =  t (i,j,k,nnew)
   ztstart(i,j) =  t (i,j,k,nnow)
  ENDDO
 ENDDO


  CALL satad ( kitera, zte, zqve, zqce, ztstart, zphfe,                        &
               zr1, zr2, zr3, zr4, zr5, zr6, zr7, zr8,                         &
               b1, b2w, b3, b4w, b234w, rdv, o_m_rdv, rvd_m_o, lh_v,cpdr,cp_d, &
               ie, je, istart, iend, jstart, jend  )

 DO i=istart,iend
  DO j=jstart,jend
   t (i,j,k,nnew) = zte (i,j)
   qv(i,j,k     ) = zqve(i,j)
   qc(i,j,k     ) = zqce(i,j)
  ENDDO
 ENDDO
ENDDO
!-------------------------------------------------------------------------------
! End of subroutine 
!-------------------------------------------------------------------------------

END SUBROUTINE lhn_satad

!===============================================================================
!+ Search of a nearby grid point with an appropriate precipitation rate
!-------------------------------------------------------------------------------
 
SUBROUTINE lhn_prof_search (i_tot,j_tot,prana,topo               &
                           ,prmod,topo_tot,prclim_tot,idim,jdim  &
                           ,nrad,lfirstsearch,isearch            &
                           ,lfound,i_found,j_found,nradius)

!-------------------------------------------------------------------------------
!
! Description:
!   This subroutine of the module "src_lheat_nudge" performs the search for
!   a grid point with a suitable precipitation rate within a specified maximal
!   search radius around the treated point. The subroutine is called by every PE.
!   
!
! Method:
!   The surrounding grid points are tested at increasing radius ranges. 
!   A suitable grid point has to have a sufficiently high precip rate and 
!   a topographic height within +/- 200m from that of the treated point. 
!   The nearest grid point satisfying this criterium is chosen. If several points
!   at the same range satisfy this limit, the point with the best fitting 
!   precipitation rate compared to the analyzed precip rate (prana) is selected.
!
!-------------------------------------------------------------------------------

! Subroutine arguments
!-------------------------------------------------------------------------------
! Scalar arguments with intent(in):
!---------------------------------
  INTEGER   (KIND=iintegers), INTENT(IN)       ::           &
    idim,jdim       ,& ! dimensions of the total model domain
    i_tot,j_tot     ,& ! coordinates of the point to treat (in total grid)
    nrad               ! radius for profile search (in grid points)
   
  REAL  (KIND=wp),     INTENT(IN)              ::           &
    prana             ,& ! precipitation rate at grid point to treat
    topo                 ! height of topography at grid point to treat
    

  LOGICAL, INTENT(INOUT)                  ::       &
    lfirstsearch     ! indicator for steps to be executed at first call only

! Array arguments with intent(in):
!---------------------------------
  REAL  (KIND=wp),     INTENT(IN)              ::           &
    prmod(idim,jdim)     ,& ! precipitation rate for the total domain
    topo_tot(idim,jdim)  ,& ! Model topography for the total domain
    prclim_tot(idim,jdim)   ! pr_clim scaled with actual latent heat release for the total domain

! Array arguments with intent(inout):
!-------------------------------------

  INTEGER   (KIND=iintegers), INTENT(INOUT)      ::           &
    isearch (4*nrad*(nrad+1),2) ,&! template holding i,j-positions
                                  ! of searched points relative to treated point
    nradius(nrad+2)        ! information of successful search at different radii
 
! Scalar arguments with intent(out):
!-------------------------------------
  LOGICAL, INTENT(OUT)                         ::           &
    lfound                   ! true if a suitable profile was found

  INTEGER   (KIND=iintegers), INTENT(OUT)      ::           &
    i_found,j_found    ! coordinates of the (best, nearest) profile found


!-------------------------------------------------------------------------------
! Local parameters, scalars, arrays :
!-------------------------------------------------------------------------------
! Local scalars:
!---------------
  LOGICAL          :: ltopo     ! flag for search failure due to height-diff

  INTEGER    (KIND=iintegers), SAVE           ::           &
    nsearch       ,& ! number of grid points tested for suitable precip rate
    ntopo            ! number of grid points with failure due to delta-topo

  REAL  (KIND=wp)                      ::           &
    ratio         ,& ! ratio of precipitation rate at treated and tested point
    best          ,& ! best precipitation ratio found so far
    pr_eps=0.5_wp    ! ratio of tt_lheat_int at a searching point to the tt_lheat_int
                     ! at the local point

  INTEGER    (KIND=iintegers)                 ::           &
    i_min,i_max   ,& ! min,max i in total domain (avoids outer nboundlines)
    j_min,j_max   ,& ! min,max j in total domain (avoids outer nboundlines)
    jn,jr         ,& ! loop counters
    i,j              ! coordinates of a grid point


!- End of header
!===============================================================================
!-------------------------------------------------------------------------------
! Begin of subroutine
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Section 1 : Set up isearch template at first call a timestep, 
!             Initialize variables
!-------------------------------------------------------------------------------

  lfound = .FALSE.
  ltopo  = .FALSE.
  i_found = 0
  j_found = 0
  i_min = 1 + nboundlines
  i_max = ie_tot - nboundlines
  j_min = 1 + nboundlines
  j_max = je_tot - nboundlines

  IF (lfirstsearch) THEN

     DO jr=1,nrad
       DO jn=1,2*jr
          isearch(4*(jr-1)*jr+jn,1) = -jr-1+jn
          isearch(4*(jr-1)*jr+jn,2) = -jr
          isearch(4*(jr-1)*jr+jn+2*jr,1) = jr
          isearch(4*(jr-1)*jr+jn+2*jr,2) = -jr-1+jn
          isearch(4*(jr-1)*jr+jn+4*jr,1) = jr+1-jn
          isearch(4*(jr-1)*jr+jn+4*jr,2) = jr
          isearch(4*(jr-1)*jr+jn+6*jr,1) = -jr
          isearch(4*(jr-1)*jr+jn+6*jr,2) = jr+1-jn
        ENDDO     ! jn
      ENDDO       ! jr

! initialize nradius
      nsearch = 0
      ntopo   = 0
      nradius(:) = 0

      lfirstsearch = .FALSE.
  ENDIF

!-------------------------------------------------------------------------------
! Section 2 : Search for a suitable point, i.e. where the modelled precipitation
!             is close to the analyzed/observed one
!-------------------------------------------------------------------------------

  ratio = 0.0_wp
  best = 0.0_wp

! Loop over ranges
   DO jr=1,nrad

! Loop over points at this range if no suitable profile was found so far
      IF (.NOT.lfound) THEN
         DO jn=1,8*jr
            i = i_tot + isearch( 4*jr*(jr-1) +jn , 1)
            j = j_tot + isearch( 4*jr*(jr-1) +jn , 2)
! check that the point is within the model domain minus outer nboundlines lines
            IF ((i < i_min) .OR. (i > i_max) .OR. (j < j_min) .OR. (j > j_max)) CYCLE

! count number of points searched
            nsearch = nsearch + 1
 
! test if difference in grid point height is not too large
            IF ( ABS(topo - topo_tot(i,j)) > 200.0_wp ) THEN
               ltopo = .TRUE.
               CYCLE
            ENDIF
            
! test if precipitation and tt_lheat_int is high enough
            IF ( prmod(i,j) > (1.0_wp - pr_eps) * prclim_tot(i,j) .AND. &
                 prmod(i,j) < (1.0_wp + pr_eps) * prclim_tot(i,j) ) THEN
              IF ( prana >= (fac_lhn_down*prmod(i,j)) .AND. &
                   prana <= (fac_lhn_up * prmod(i,j)) ) THEN
                 ratio = prana / (prmod(i,j)+repsilon)
                 IF ( ratio > 1.0_wp ) ratio = 1.0_wp/ratio

  
! keep record of the best match at this range
                 IF (ratio > best) THEN
                     best = ratio
                     i_found = i
                     j_found = j
                     IF (.NOT.lfound) nradius(jr)=nradius(jr) + 1
                     lfound = .TRUE.
                 ENDIF ! ratio test
              ENDIF ! ratio test
            ENDIF    ! precip rate test
         ENDDO       ! jn points loop
      ENDIF          ! found one test
   ENDDO             ! jr range loop

   IF (ltopo .AND. .NOT.lfound) ntopo = ntopo + 1
   nradius(nrad+1) = nsearch
   nradius(nrad+2) = ntopo

!-------------------------------------------------------------------------------
! End of subroutine 
!-------------------------------------------------------------------------------

END SUBROUTINE lhn_prof_search

!===============================================================================
!+ Filtering of vertical (heating) profiles, elimination of isolated peaks
!-------------------------------------------------------------------------------

SUBROUTINE filter_prof (prof,prof_filt,kup,klow,eps,lelim,lsmooth, &
                        nelimosc,nelimiso,nsmooth)

!-------------------------------------------------------------------------------
! Description:
!  This subroutines filters a vertical profile (e.g. heating profile to be used
!  in latent heat nudging to elimiate computational noise).
!
! Method: 
!  lelim : eliminate isolated peaks of small vertical extent
!    a) value on one level below and above is below specified eps
!    b) value two levels above is below eps and less than one level below the
!       profile value was idntified as very small or isolated peak
!
!  lsmooth : apply a simple one-dimensional shapiro-filter with S=1/2
!            to levels where the value is above eps
!
!-------------------------------------------------------------------------------
 
! Subroutine arguments, intent=in and intent=inout
  LOGICAL, INTENT(IN)       ::    &
    lelim             ,& ! flag 0/1 for elimination of isolated peaks
    lsmooth              ! flag 0/1 for smoothing

  INTEGER (KIND=iintegers), INTENT(IN)       ::    &
    kup               ,& ! array dimensions of profile prof
    klow                 ! array dimensions of profile prof

  REAL (KIND=wp),        INTENT(IN)       ::      &
    prof(kup:klow)    ,& ! profile array
    eps                  ! limits above which values are modified

  REAL (KIND=wp),        INTENT(OUT)       ::      &
    prof_filt(kup:klow)      ! array for filtered profile

! Local subroutine variables and arrays
  LOGICAL                                   ::    &
    lflag(kup:klow)          ! array to flag levels to be set to zero
  INTEGER   (KIND=iintegers)                ::    &
    k,nheat, &               ! loop index, acceptable level counter
    nelimosc,nelimiso,nsmooth ! diag. output
  REAL (KIND=wp)               ::      &
    proffilt(kup:klow)       ! profile array

!- End of header
!-------------------------------------------------------------------------------
! Begin Subroutine filter_prof
!-------------------------------------------------------------------------------

   proffilt(:)=prof(:)
   prof_filt(:)=prof(:)

   nelimosc = 0
   nelimiso = 0
   nsmooth  = 0

! eliminate isolated peaks
   IF (lelim) THEN

     lflag(:) = .FALSE.

!    eliminate oscillations around zero between adjacent levels
     IF ( ABS(prof(klow)) < eps ) lflag(klow) = .TRUE.
     IF ( ABS(prof(kup)) < eps ) lflag(kup) = .TRUE.
     DO k=klow-1,kup+1,-1
        IF ( (prof(k-1)*prof(k)) <= 0.0_wp .AND.   &
             (prof(k)*prof(k+1)) <= 0.0_wp         ) THEN
           lflag(k) = .TRUE.
           nelimosc = nelimosc + 1
        ENDIF
     ENDDO
     WHERE ( lflag ) proffilt = 0.0_wp

!    eliminate isolated peaks of small vertical extent
!    a) value on one level below and above is below specified eps
!    b) value two levels above is below eps and less than one level below the
!       profile value was idntified as very small or isolated peak
!    nheat : counter of lower levels with accepted heating rate (profile value)
!          - is reset to zero when an isolated peak is diagnosed or the heating
!            at the level considered is below specified eps
     nheat = 0
     lflag(:) = .FALSE.
     DO k=klow-1,kup+2,-1
        nheat=nheat+1
        IF ( ABS(proffilt(k-1))  <= eps .AND. ABS(proffilt(k+1)) <= eps ) THEN
              lflag(k) = .TRUE.
              nheat = 0
           nelimiso = nelimiso + 1
        ENDIF
        IF ( nheat < 2 .AND. ABS(proffilt(k-2)) <= eps ) THEN
              lflag(k-1) = .TRUE.
              lflag(k) = .TRUE.
              nheat = 0
           nelimiso = nelimiso + 1
        ENDIF
     ENDDO
     WHERE ( lflag ) proffilt = 0.0_wp

     prof_filt(:)=proffilt(:)

   ENDIF

! smooth profile
   IF (lsmooth) THEN

      DO k=klow-1,kup+1,-1
         IF ( ABS(proffilt(k)) >= eps ) THEN
            prof_filt(k) = 0.5_wp  * proffilt(k) + 0.25_wp * (proffilt(k+1)+proffilt(k-1))
            nsmooth = nsmooth + 1
         ENDIF
      ENDDO
      IF ( ABS(proffilt(klow)) >= eps) THEN
         prof_filt(klow)=0.66_wp * proffilt(klow) + 0.33_wp * proffilt(klow-1)
         nsmooth = nsmooth + 1
      ENDIF
      IF ( ABS(proffilt(kup)) >= eps) THEN
         prof_filt(kup)=0.66_wp * proffilt(kup) + 0.33_wp * proffilt(kup+1)
         nsmooth = nsmooth + 1
      ENDIF
   ENDIF

!-------------------------------------------------------------------------------

END SUBROUTINE filter_prof

!===============================================================================

SUBROUTINE hor_filt(field,nfilt,kup,klow,pos_def)

!-------------------------------------------------------------------------------
! Description:
!
!-------------------------------------------------------------------------------

! Subroutine / Function arguments
!-------------------------------------------------------------------------------

REAL (KIND=wp)     ,INTENT(INOUT) :: &
  field(ie,je,ke)

INTEGER (KIND=iintegers) ,INTENT(IN) :: &
  kup,klow,nfilt

LOGICAL, INTENT(IN)  :: &
  pos_def

! Local scalars:
! -------------
INTEGER (KIND=iintegers) ::  &
  izerror,kzdims(24),  &
  i,j,k,n, &
  hfjstartpar,hfjendpar

REAL (KIND=wp)     , ALLOCATABLE :: &
  field_tmp(:,:,:)

CHARACTER (LEN=200)        ::  &
  yzerrmsg    ! error message for error handling


  izerror   = 0
  yroutine  = 'hor_filt'
  yzerrmsg  = '     '

  ! width of the stencil for the horizontal filter
  hfwidth  = 4
  hfw_m_nb = hfwidth - nboundlines
  ie_hf = ie + 2*hfw_m_nb
  je_hf = je + 2*hfw_m_nb
  ALLOCATE( field_tmp(ie_hf,je_hf,kup:klow), STAT = izlocstat )
  IF (izlocstat /= 0) THEN
     yerrmsg =' ERROR  *** allocation of space for lhn - fields failed'
     CALL model_abort (my_cart_id,6407,yerrmsg,yroutine)
  ENDIF

  IF (nfilt > 0_iintegers) THEN
     CALL init_horizontal_filtering_lh( field(:,:,kup:klow), field_tmp(:,:,kup:klow),  &
                                        kup, klow, hfjstartpar, hfjendpar )
  ENDIF

  DO n = 1, nfilt

    CALL horizontal_filtering_lh( field_tmp(:,:,kup:klow), kup, klow )

    ! exchange boundaries after filtering
    IF (num_compute > 1) THEN
    kzdims(1:24) =(/ klow-kup+1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,        &
                     0,0,0,0 /)
    CALL exchg_boundaries                                            &
      ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie_hf,       &
        je_hf, kzdims, hfjstartpar, hfjendpar, hfwidth, hfwidth,     &
        my_cart_neigh, lperi_x, lperi_y, l2dim, 2000, .FALSE., 1, izerror, yzerrmsg,          &
        field_tmp(:,:,kup:klow) )
    ENDIF

  END DO

  IF (nfilt > 0_iintegers) THEN
     DO k = kup, klow
       DO j = 1, je
         DO i = 1, ie
           field(i,j,k) = field_tmp(i+hfw_m_nb,j+hfw_m_nb,k)
           IF (pos_def) THEN
              IF (field(i,j,k) < 0.0_wp) field(i,j,k) = 0.0_wp
           ENDIF
         END DO
       END DO
     END DO
  ENDIF

  DEALLOCATE( field_tmp, STAT = izlocstat )

END SUBROUTINE hor_filt

!===============================================================================

SUBROUTINE init_horizontal_filtering_lh ( field, field_hf, kstart, kend, &
                                          hfjstartpar, hfjendpar )
!-------------------------------------------------------------------------------
!
! Description:
!
!-------------------------------------------------------------------------------
! Subroutine / Function arguments
!-------------------------------------------------------------------------------

INTEGER (KIND=iintegers), INTENT(IN) :: &
     kstart, kend
INTEGER (KIND=iintegers), INTENT(OUT) :: &
     hfjstartpar, hfjendpar
REAL    (KIND=wp   ),     INTENT (IN) ::  &

     field(ie,je,kstart:kend)
REAL    (KIND=wp   ),     INTENT (OUT) ::  &
     field_hf(ie_hf,je_hf,kstart:kend)

! Local scalars:
! -------------
INTEGER (KIND=iintegers) ::  &
     izerror, i, j, k, kzdims(24) !,  &
!      istart_hf, iend_hf, jstart_hf, jend_hf

CHARACTER (LEN=200)        ::  &
  yzerrmsg    ! error message for error handling


  izerror   = 0
  yzerrmsg  = '   '


  field_hf(:,:,:) = 0.0_wp

  DO k = kstart, kend
    DO j = 1, je
      DO i = 1, ie
        field_hf(i+hfw_m_nb,j+hfw_m_nb,k) = field(i,j,k)
      END DO
    END DO
  END DO
  ! Determine start- and end-indices for communication and
  ! set values in boundline frame
  ! west
  IF (my_cart_neigh(1) == -1) THEN
     DO k = kstart, kend
        DO j = 1,je
           DO i= 1, hfw_m_nb
              field_hf(i,j+hfw_m_nb,k) = field(1,j,k)
           ENDDO
        ENDDO
     ENDDO
    ! southwest corner
     IF (my_cart_neigh(4) == -1) THEN
        DO k = kstart, kend
           DO j = 1, hfw_m_nb
              DO i= 1, hfw_m_nb
                 field_hf(i,j,k) = field(1,1,k)
              ENDDO
           ENDDO
        ENDDO
     ENDIF
    ! northwest corner
     IF (my_cart_neigh(2) == -1) THEN
        DO k = kstart, kend
           DO j = je+1, je+hfw_m_nb
              DO i = 1, hfw_m_nb
                 field_hf(i,j+hfw_m_nb,k) = field(1,je,k)
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDIF
  ! east
  IF (my_cart_neigh(3) == -1) THEN
     DO k = kstart, kend
        DO j = 1, je
           DO i = ie+1, ie+hfw_m_nb
              field_hf(i+hfw_m_nb,j+hfw_m_nb,k) = field(ie,j,k)
           ENDDO
        ENDDO
     ENDDO
    ! southeast corner
     IF (my_cart_neigh(4) == -1) THEN
        DO k = kstart, kend
           DO j = 1, hfw_m_nb
              DO i = ie+1, ie+hfw_m_nb 
                 field_hf(i+hfw_m_nb,j,k) = field(ie,1,k)
              ENDDO
           ENDDO
        ENDDO
     ENDIF
    ! northeast corner
     IF (my_cart_neigh(2) == -1) THEN
        DO k = kstart, kend
           DO j = je+1, je+hfw_m_nb
              DO i = ie+1, ie+hfw_m_nb
                 field_hf(i+hfw_m_nb,j+hfw_m_nb,k) = field(ie,je,k)
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDIF
  ! south
  IF (my_cart_neigh(4) == -1) THEN
     hfjstartpar = 1
     DO k = kstart, kend
        DO j = 1, hfw_m_nb
           DO i = 1, ie
              field_hf(i+hfw_m_nb,j,k) = field(i,1,k)
           ENDDO
        ENDDO
     ENDDO
  ELSE
     hfjstartpar = 1 + hfw_m_nb + nboundlines
  ENDIF
  ! north
  IF (my_cart_neigh(2) == -1) THEN
     hfjendpar   = je_hf
     DO k = kstart, kend
        DO j = je+1, je+hfw_m_nb
           DO i = 1, ie
              field_hf(i+hfw_m_nb,j+hfw_m_nb,k) = field(i,je,k)
           ENDDO
        ENDDO
     ENDDO
  ELSE
     hfjendpar   = je_hf - nboundlines + hfw_m_nb
  ENDIF

  IF (num_compute > 1) THEN
    kzdims(1:24) = (/ kend-kstart+1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  &
                      0,0,0,0 /)
    CALL exchg_boundaries                                                   &
      ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie_hf,    &
        je_hf, kzdims, hfjstartpar, hfjendpar, hfwidth, hfwidth,  &
        my_cart_neigh, lperi_x, lperi_y, l2dim, 2000, .FALSE., 1, izerror, yzerrmsg, field_hf )
  ENDIF

END SUBROUTINE init_horizontal_filtering_lh

!===============================================================================

SUBROUTINE horizontal_filtering_lh (field_hf, kstart, kend)

! -------------------------------------------------------------------------
! Description:
!
!-------------------------------------------------------------------------------

! Subroutine / Function arguments
!-------------------------------------------------------------------------------

INTEGER (KIND=iintegers), INTENT(IN) :: &
  kstart, kend
REAL    (KIND=wp   ),     INTENT (INOUT) ::  &
  field_hf(ie_hf,je_hf,kstart:kend)

! Local scalars:
! -------------
INTEGER (KIND=iintegers) ::  &
  ilow, iup, istart, iend,   & !
  jlow, jup, jstart, jend,   & !
  i, j, k                      !  Loop indices

! Local (automatic) arrays:
! -------------------------
REAL    (KIND=wp   )     ::  &
  field_tmp(ie_hf,je_hf,kstart:kend), &
  zfw9p(-hfwidth:hfwidth),   & ! filter weights for 9-point filter
  zfw3p(-1:1)                  ! filter weights for 3-point filter

  hfw_m_nb = hfwidth - nboundlines
  istart = 1 + nboundlines
  iend   = ie_hf - 2*hfw_m_nb - nboundlines
  jstart = 1 + nboundlines
  jend   = je_hf - 2*hfw_m_nb - nboundlines

  ! init working array
  field_tmp(:,:,:) = field_hf(:,:,:)

  ! coefficiens as in dfilt4
!!$  zfw9p = (/  -0.00390625,              &
!!$               0.03125,                 &
!!$              -0.109375,                &
!!$               0.21875,                 &
!!$               0.7265625,               &
!!$               0.21875,                 &
!!$              -0.109375,                &
!!$               0.03125,                 &
!!$              -0.00390625 /)

  zfw9p = (/   0.0_wp,                 &
               0.0_wp,                 &
               0.0_wp,                 &
               0.25_wp,                &
               0.5_wp,                 &
               0.25_wp,                &
               0.0_wp,                 &
               0.0_wp,                 &
               0.0_wp /)

  zfw3p = (/ 0.25_wp, 0.5_wp, 0.25_wp /)

! west
  IF (my_cart_neigh(1) == -1) THEN
     ilow = 1 + 2*hfwidth
  ELSE
     ilow = istart + hfw_m_nb
  END IF
  ! east
  IF (my_cart_neigh(3) == -1) THEN
     iup = iend - nboundlines
  ELSE
     iup = iend + hfw_m_nb
  END IF
  ! south
  IF (my_cart_neigh(4) == -1) THEN
     jlow = 1 + 2*hfwidth
  ELSE
     jlow = jstart + hfw_m_nb
  END IF
  ! north
  IF (my_cart_neigh(2) == -1) THEN
     jup = jend - nboundlines
  ELSE
     jup = jend + hfw_m_nb
  END IF

  !
  ! apply 9-point-filter in x-direction
  !
  DO k = kstart, kend
     DO j = 1, je_hf
        DO i = ilow, iup
           field_tmp(i,j,k) = 0.0_wp
!US        DO l = -hfwidth, hfwidth
!US           field_tmp(i,j,k) = field_tmp(i,j,k) + zfw9p(l)*field_hf(i+l,j,k)
!US        ENDDO

           field_tmp(i,j,k) = field_tmp(i,j,k) + zfw9p(-1)*field_hf(i-1,j,k)    &
                                               + zfw9p( 0)*field_hf(i  ,j,k)    &
                                               + zfw9p( 1)*field_hf(i+1,j,k)
        ENDDO
     ENDDO
  ENDDO
  ! apply 3-point-filter in x-direction at west boundary
  IF (my_cart_neigh(1) == -1) THEN
     DO k = kstart, kend
        DO j = 1, je_hf
           DO i = hfw_m_nb+1, ilow-1
              field_tmp(i,j,k) = 0.0_wp
!US           DO l = -1, 1
!US              field_tmp(i,j,k) = field_tmp(i,j,k) + zfw3p(l)*field_hf(i+l,j,k)
!US           ENDDO

              field_tmp(i,j,k) = field_tmp(i,j,k) + zfw3p(-1)*field_hf(i-1,j,k)    &
                                                  + zfw3p( 0)*field_hf(i  ,j,k)    &
                                                  + zfw3p( 1)*field_hf(i+1,j,k)
           ENDDO
        ENDDO
     ENDDO
  ENDIF

  ! apply 3-point-filter in x-direction at east boundary
  IF (my_cart_neigh(3) == -1) THEN
     DO k = kstart, kend
        DO j = 1, je_hf
           DO i = iup+1, ie+hfw_m_nb
              field_tmp(i,j,k) = 0.0_wp
!US           DO l = -1, 1
!US              field_tmp(i,j,k) = field_tmp(i,j,k) + zfw3p(l)*field_hf(i+l,j,k)
!US           ENDDO

              field_tmp(i,j,k) = field_tmp(i,j,k) + zfw3p(-1)*field_hf(i-1,j,k)    &
                                                  + zfw3p( 0)*field_hf(i  ,j,k)    &
                                                  + zfw3p( 1)*field_hf(i+1,j,k)
           ENDDO
        ENDDO
     ENDDO
  ENDIF


  !
  ! apply 9-point-filter in y-direction
  !
  DO k = kstart, kend
     DO j = jlow, jup
        DO i = 1, ie_hf
           field_hf(i,j,k) = 0.0_wp
!US        DO l = -hfwidth, hfwidth
!US           field_hf(i,j,k) = field_hf(i,j,k) + zfw9p(l)*field_tmp(i,j+l,k)
!US        ENDDO

           field_hf(i,j,k) = field_hf(i,j,k) + zfw9p(-1)*field_tmp(i,j-1,k)    &
                                             + zfw9p( 0)*field_tmp(i,j  ,k)    &
                                             + zfw9p( 1)*field_tmp(i,j+1,k)
        ENDDO
     ENDDO
  ENDDO


  ! apply 3-point-filter in y-direction at south boundary
  IF (my_cart_neigh(4) == -1) THEN
     DO k = kstart, kend
        DO j = hfw_m_nb+1, jlow-1
           DO i = 1, ie_hf
              field_hf(i,j,k) = 0.0_wp
!US           DO l = -1, 1
!US              field_hf(i,j,k) = field_hf(i,j,k) + zfw3p(l)*field_tmp(i,j+l,k)
!US           ENDDO

           field_hf(i,j,k) = field_hf(i,j,k) + zfw3p(-1)*field_tmp(i,j-1,k)    &
                                             + zfw3p( 0)*field_tmp(i,j  ,k)    &
                                             + zfw3p( 1)*field_tmp(i,j+1,k)
           ENDDO
        ENDDO
     ENDDO
  ENDIF

  ! apply 3-point-filter in y-direction at north boundary
  IF (my_cart_neigh(2) == -1) THEN
     DO k = kstart, kend
        DO j = jup+1, je+hfw_m_nb
           DO i = 1, ie_hf
              field_hf(i,j,k) = 0.0_wp
!US           DO l = -1, 1
!US              field_hf(i,j,k) = field_hf(i,j,k) + zfw3p(l)*field_tmp(i,j+l,k)
!US           ENDDO

           field_hf(i,j,k) = field_hf(i,j,k) + zfw3p(-1)*field_tmp(i,j-1,k)    &
                                             + zfw3p( 0)*field_tmp(i,j  ,k)    &
                                             + zfw3p( 1)*field_tmp(i,j+1,k)
           ENDDO
        ENDDO
     ENDDO
  ENDIF

!------------------------------------------------------------------------------
! End of subroutine horizontal_filtering_lh
!------------------------------------------------------------------------------

END SUBROUTINE horizontal_filtering_lh

!===============================================================================

SUBROUTINE lhn_verification (ytime,zprmod,zprmod_ref,zprrad,zprmodatdx,zprmod_ref_f,zprrad_f)

!-------------------------------------------------------------------------------
!
! Description:
! This subroutine calculates different parameter for verification of model precipitaion
! against radar observations. The values are stored in the LHN output file YULHN.
!-------------------------------------------------------------------------------

! Subroutine / Function arguments
!-------------------------------------------------------------------------------

 CHARACTER (LEN=*), INTENT(IN)  ::       &
   ytime

 REAL (KIND=wp),     INTENT(IN)  ::       &
   zprmod(ie,je),                    &
   zprmod_ref(ie,je),                &
   zprmodatdx(ie,je),                &
   zprrad(ie,je)

 REAL (KIND=wp),     INTENT(IN), OPTIONAL  ::       &
   zprmod_ref_f(ie,je),              &
   zprrad_f(ie,je)

! Local scalars:
! -------------

 REAL (KIND=wp)             ::       &
   zflar,                            & !sum of diagnostic/prognostic precipitation
   timefac,                          &
   zprmod_s,                         &
   zprmod_ref_s,                     &
   zprrad_s,                         &
   zprmodatdx_s,                     &
   zprmod_ref_f_s,                   &
   zprrad_f_s,                       &
   zprmod_ano,                       &
   zprmod_ref_ano,                   &
   zprrad_ano,                       &
   zprmod_var,                       &
   zprmod_ref_var,                   &
   zprrad_var,                       &
   zprcovar,                         &
   zprcovar_ref,                     &
   zprdiff,                          &
   zprmse,                           &
   zprmae,                           &
   zprcorrel,                        &
   zprcorrel_ref,                    &
   zprcme

 REAL (KIND=wp)             ::       &
   realbuf      (7)         ! for communication


 INTEGER (KIND=iintegers) :: &
   i,j,n,nrealbuf,i_ver(ie*je),j_ver(ie*je)

  INTEGER (KIND=iintegers)         ::           &
   zpranz, zprcount

   IF ( ytime == 'SW') THEN
      nrealbuf=6
      timefac=3600.0_wp
   ELSE
      nrealbuf=4
      timefac=1.0_wp
   ENDIF


   zprmod_s       = 0.0_wp
   zprmod_ref_s   = 0.0_wp
   zprmod_ref_f_s = 0.0_wp
   zprrad_s       = 0.0_wp
   zprmodatdx_s   = 0.0_wp
   zprrad_f_s     = 0.0_wp
   zprmod_ano     = 0.0_wp
   zprmod_ref_ano = 0.0_wp
   zprrad_ano     = 0.0_wp
   zprmod_var     = 0.0_wp
   zprmod_ref_var = 0.0_wp
   zprrad_var     = 0.0_wp
   zprcovar       = 0.0_wp
   zprcovar_ref   = 0.0_wp
   zprdiff        = 0.0_wp
   zprmse         = 0.0_wp
   zprmae         = 0.0_wp
   zprcorrel      = 0.0_wp
   zprcorrel_ref  = 0.0_wp
   zprcme         = 0.0_wp
   zprcount       = 0_iintegers
   zpranz         = 0_iintegers

   DO j = jstart,jend
     DO i = istart,iend
       IF (wobs_space(i,j) > 0.75_wp            .AND.  &
           NINT(blacklist(i,j)) /= 1_iintegers  .AND.  &
           NINT(brightband(i,j)) /= 1_iintegers .AND.  &
           zprrad(i,j) >= 0.0_wp) THEN
           zprmod_s        = zprmod_s        + zprmod(i,j)
           zprmod_ref_s    = zprmod_ref_s    + zprmod_ref(i,j)
           zprrad_s        = zprrad_s        + zprrad(i,j)
           zprmodatdx_s    = zprmodatdx_s    + zprmodatdx(i,j)
           IF (PRESENT (zprmod_ref_f)) zprmod_ref_f_s  = zprmod_ref_f_s  + zprmod_ref_f(i,j)
           IF (PRESENT (zprrad_f)) zprrad_f_s      = zprrad_f_s      + zprrad_f(i,j)
           zprcount        = zprcount        + 1_iintegers
           i_ver(zprcount) = i
           j_ver(zprcount) = j
       ENDIF
     ENDDO
   ENDDO

   IF (num_compute > 1) THEN
      realbuf    = 0.0_wp
      realbuf(1) = zprmod_s
      realbuf(2) = zprmod_ref_s
      realbuf(3) = zprrad_s
      realbuf(4) = zprmodatdx_s
      IF (PRESENT (zprmod_ref_f)) realbuf(5) = zprmod_ref_f_s
      IF (PRESENT (zprrad_f))     realbuf(6) = zprrad_f_s

      IF (ltime) THEN
        CALL get_timings (i_lhn_computations, ntstep, dt, izerror)
        IF (ltime_barrier) THEN
          CALL comm_barrier (icomm_cart, izerror, yerrmsg)
          CALL get_timings (i_barrier_waiting_lhn, ntstep, dt, izerror)
        ENDIF
      ENDIF

      CALL global_values (realbuf, nrealbuf, 'SUM', imp_reals, icomm_cart, -1,    &
          yerrmsg, izerror)

      IF (ltime) THEN
         CALL get_timings (i_communications_lhn, ntstep, dt, izerror)
      ENDIF

      zprmod_s       = realbuf(1)
      zprmod_ref_s   = realbuf(2)
      zprrad_s       = realbuf(3)
      zprmodatdx_s   = realbuf(4)
      IF (PRESENT (zprmod_ref_f)) zprmod_ref_f_s = realbuf(5)
      IF (PRESENT (zprrad_f))     zprrad_f_s     = realbuf(6)

      iexchange      = 0
      iexchange( 1)  = zprcount

      IF (ltime) THEN
        CALL get_timings (i_lhn_computations, ntstep, dt, izerror)
        IF (ltime_barrier) THEN
          CALL comm_barrier (icomm_cart, izerror, yerrmsg)
          CALL get_timings (i_barrier_waiting_lhn, ntstep, dt, izerror)
        ENDIF
      ENDIF

      CALL global_values (iexchange, 1, 'SUM', imp_integers, icomm_cart, -1,    &
          yerrmsg, izerror)

      IF (ltime) THEN
         CALL get_timings (i_communications_lhn, ntstep, dt, izerror)
      ENDIF

      zpranz = iexchange( 1)

   ENDIF

  IF (zpranz > 0) THEN
      zflar          = 1.0_wp / REAL (zpranz,wp)
      zprmod_s       = zprmod_s       * zflar * timefac
      zprmod_ref_s   = zprmod_ref_s   * zflar * timefac
      zprrad_s       = zprrad_s       * zflar * timefac
      zprmodatdx_s   = zprmodatdx_s   * zflar * timefac
      IF (PRESENT (zprmod_ref_f)) zprmod_ref_f_s = zprmod_ref_f_s * zflar * timefac
      IF (PRESENT (zprrad_f))     zprrad_f_s     = zprrad_f_s     * zflar * timefac

   DO n = 1,zprcount
       i=i_ver(n)
       j=j_ver(n)
       zprdiff         = zprmod(i,j)*timefac - zprrad(i,j)*timefac
       zprmse          = zprmse          + (zprdiff*zprdiff)
       zprmae          = zprmae          + ABS(zprdiff)
       zprmod_ano      = zprmod(i,j)*timefac - zprmod_s
       zprmod_var      = zprmod_var      + (zprmod_ano*zprmod_ano)
       zprmod_ref_ano  = zprmod_ref(i,j)*timefac - zprmod_ref_s
       zprmod_ref_var  = zprmod_ref_var  + (zprmod_ref_ano*zprmod_ref_ano)
       zprrad_ano      = zprrad(i,j)*timefac - zprrad_s
       zprrad_var      = zprrad_var      + (zprrad_ano*zprrad_ano)
       zprcovar        = zprcovar        + (zprmod_ano*zprrad_ano)
       zprcovar_ref    = zprcovar_ref    + (zprmod_ref_ano*zprrad_ano)
   ENDDO

   IF (num_compute > 1) THEN
      realbuf    = 0.0_wp
      realbuf(1) = zprmod_var
      realbuf(2) = zprmod_ref_var
      realbuf(3) = zprrad_var
      realbuf(4) = zprcovar
      realbuf(5) = zprcovar_ref
      realbuf(6) = zprmse
      realbuf(7) = zprmae

      IF (ltime) THEN
        CALL get_timings (i_lhn_computations, ntstep, dt, izerror)
        IF (ltime_barrier) THEN
          CALL comm_barrier (icomm_cart, izerror, yerrmsg)
          CALL get_timings (i_barrier_waiting_lhn, ntstep, dt, izerror)
        ENDIF
      ENDIF

      CALL global_values (realbuf, 7, 'SUM', imp_reals, icomm_cart, 0,    &
          yerrmsg, izerror)

      IF (ltime) THEN
         CALL get_timings (i_communications_lhn, ntstep, dt, izerror)
      ENDIF

      zprmod_var     = realbuf(1)
      zprmod_ref_var = realbuf(2)
      zprrad_var     = realbuf(3)
      zprcovar       = realbuf(4)
      zprcovar_ref   = realbuf(5)
      zprmse         = realbuf(6)
      zprmae         = realbuf(7)

   ENDIF
  ENDIF

  IF (my_cart_id == 0) THEN
      zprmod_var     = zprmod_var     * zflar
      zprmod_ref_var = zprmod_ref_var * zflar
      zprrad_var     = zprrad_var     * zflar
      zprcovar       = zprcovar       * zflar
      zprcovar_ref   = zprcovar_ref   * zflar
      zprmse         = zprmse         * zflar
      zprmae         = zprmae         * zflar
      IF ((zprmod_var + zprrad_var) /= 0.0_wp) &
         zprcorrel      = zprcovar/EXP(0.5_wp*LOG(zprmod_var + zprrad_var))
      IF ((zprmod_ref_var + zprrad_var) /= 0.0_wp) &
         zprcorrel_ref  = zprcovar_ref/EXP(0.5_wp*LOG(zprmod_ref_var + zprrad_var))
      IF ((zprrad_var) /= 0.0_wp) &
         zprcme         = 1.0_wp - (zprmse/zprrad_var)
!   ENDIF

!   IF (my_cart_id == 0) THEN
      WRITE(nulhn, *)'Verification:'
      WRITE(nulhn, '(a,a3,i6,3f8.3,2f12.3,f8.3)')'Modell (mod,ref,filt,var,var_ref,modatdx)',ytime,ntstep,zprmod_s,zprmod_ref_s &
           ,zprmod_ref_f_s,zprmod_var,zprmod_ref_var,zprmodatdx_s
      WRITE(nulhn, '(a,a3,i6,2f8.3,f12.3)')'Radar (obs,filt,var)',ytime,ntstep,zprrad_s,zprrad_f_s,zprrad_var
      WRITE(nulhn, '(a,a3,i6,9f12.3)')'Statist (bias,mse,rmse,mae,covar,correl,cme,covar_ref,correl_ref)' &
      ,ytime,ntstep,zprmod_s-zprrad_s,zprmse,SQRT(zprmse),zprmae,zprcovar,zprcorrel,zprcme,zprcovar_ref,zprcorrel_ref

   ENDIF

   CALL lhn_skill_scores (ytime,zprmod,zprrad)
!   CALL lhn_skill_scores (ytime,zprmodatdx,zprrad)



END SUBROUTINE lhn_verification

!===============================================================================

SUBROUTINE lhn_skill_scores (ytime,zprmod,zprrad)

!-------------------------------------------------------------------------------
!
! Description:
! This subroutine calculates different skill scores for verification of model precipitaion
! against radar observations. The scores are stored in the LHN output file YULHN.
!-------------------------------------------------------------------------------

! Subroutine / Function arguments
!-------------------------------------------------------------------------------
 CHARACTER (LEN=*), INTENT(IN)            ::       &
   ytime

 REAL (KIND=wp),     INTENT(IN)  ::       &
   zprmod(ie,je),                    &
   zprrad(ie,je)


! Local scalars:
! -------------
  LOGICAL, SAVE                     ::           &
   lfirst=.TRUE.      ! flag for first call of this subroutine

  INTEGER (KIND=iintegers) ::  &
   i,j,ass,bss,css,dss,zss,ierror       ! table of contengency

  INTEGER (KIND=iintegers) ::  &
   ii,jj,  &
   nthre ,&
   itab(7,7),i1,i2,ith, &
   histmod(7),histobs(7),anzobs,anzmod

  REAL (KIND=wp)     ::  &
   rass,rbss,rcss,rdss     ! table of contengency as real

  REAL (KIND=wp)                   ::           &
   timefac,hr, far, fr, fbi, ts, rets, rhss, ets, pai_plus, pai_minus, pod, tss, hss ! skill scores

  REAL (KIND=wp)                   ::           &
   thr(6)

   nthre=6

   IF ( ytime == 'SW') THEN
      timefac=3600.0_wp
   ELSE
      timefac=1.0_wp
   ENDIF


   hr=0._wp
   fr=0._wp
   far=0._wp
   fbi=0._wp
   ts=0._wp
   ets=0._wp
   hss=0._wp
   pai_plus=0._wp
   pai_minus=0._wp
   tss=0._wp
   pod=0._wp
   rhss=0.0_wp
   rets=0.0_wp
   thr = (/ 0.1_wp, 0.2_wp, 0.5_wp, 1.0_wp, 2.0_wp, 5.0_wp /)
   thr = thr / timefac
   histobs(:) = 0_iintegers
   histmod(:) = 0_iintegers
   anzobs=0_iintegers
   anzmod=0_iintegers

! contingence table:
!             |   Observed    |
!       -----------------------
!       |     |  yes  |   no  |
! -----------------------------
! Mod   | yes |  ass  |  bss  |
!       -----------------------
! elled | no  |  css  |  dss  |
! -----------------------------

   DO jj=1,nthre+1
    DO ii=1,nthre+1
     itab(ii,jj)=0
    ENDDO
   ENDDO
   DO j=jstart,jend
      DO i=istart,iend
       IF (wobs_space(i,j) > 0.75_wp            .AND. &
           NINT(blacklist(i,j)) /= 1_iintegers  .AND. &
           NINT(brightband(i,j)) /= 1_iintegers) THEN
           i1=1
           i2=1
           IF (zprrad(i,j) > 0.0_wp) THEN
              histobs(1)=histobs(1)+1_iintegers
           ENDIF
           IF (zprmod(i,j) > 0.0_wp) THEN
              histmod(1)=histmod(1)+1_iintegers
           ENDIF
           DO ith=1,nthre
            IF (zprrad(i,j).GE.thr(ith)) THEN
             i1=i1+1
             histobs(ith+1)=histobs(ith+1)+1_iintegers
            ENDIF
            IF (zprmod(i,j).GE.thr(ith)) THEN
             i2=i2+1
             histmod(ith+1)=histmod(ith+1)+1_iintegers
            ENDIF
           ENDDO
           itab(i1,i2)=itab(i1,i2)+1
           anzobs=anzobs+1_iintegers
           anzmod=anzmod+1_iintegers
       ENDIF
     ENDDO
   ENDDO

   DO ith=1,nthre

     ass=0_iintegers
     bss=0_iintegers
     css=0_iintegers
     dss=0_iintegers
     zss=0_iintegers
!  Observation no / Forecast no
     DO j=1,ith
      DO i=1,ith
       dss=dss+itab(i,j)
      ENDDO
     ENDDO
!  Observation no / Forecast yes
     DO j=ith+1,nthre+1
      DO i=1,ith
       bss=bss+itab(i,j)
      ENDDO
     ENDDO
!  Observation yes / Forecast no
     DO j=1,ith
      DO i=ith+1,nthre+1
       css=css+itab(i,j)
      ENDDO
     ENDDO
!  Observation yes / Forecast yes
     DO  j=ith+1,nthre+1
      DO  i=ith+1,nthre+1
       ass=ass+itab(i,j)
      ENDDO
     ENDDO


  zss=ass+bss+css+dss

! calculate skill scores

     IF (num_compute > 1) THEN
      iexchange = 0
      iexchange( 1)= ass
      iexchange( 2)= bss
      iexchange( 3)= css
      iexchange( 4)= dss
      iexchange( 5)= zss
      iexchange( 6)= histmod(ith)
      iexchange( 7)= histobs(ith)

      IF (ltime) THEN
        CALL get_timings (i_lhn_computations, ntstep, dt, izerror)
        IF (ltime_barrier) THEN
          CALL comm_barrier (icomm_cart, izerror, yerrmsg)
          CALL get_timings (i_barrier_waiting_lhn, ntstep, dt, izerror)
        ENDIF
      ENDIF

      CALL global_values (iexchange,7, 'SUM',imp_integers,icomm_cart, 0,yerrmsg, ierror)

      IF (ltime) THEN
         CALL get_timings (i_communications_lhn, ntstep, dt, izerror)
      ENDIF
     ENDIF

  IF(my_cart_id == 0) THEN
     IF (num_compute > 1) THEN
        ass    = iexchange( 1)
        bss    = iexchange( 2)
        css    = iexchange( 3)
        dss    = iexchange( 4)
        zss    = iexchange( 5)
        histmod(ith) = iexchange( 6)
        histobs(ith) = iexchange( 7)
     ENDIF
   rass=REAL(ass,wp)
   rbss=REAL(bss,wp)
   rcss=REAL(css,wp)
   rdss=REAL(dss,wp)

   ! new nomenclature according to U. Damrath in "Die neue Modellkette des DWD II"
   ! hitrate or percent correct score
   IF ((rass+rbss+rcss+rdss) > 0._wp) &
    hr  = 100.0_wp * (rass+rdss)/(rass+rbss+rcss+rdss)
   ! Probability of detection
   IF ((rass+rcss) > 0._wp) &
    pod  = 100.0_wp * (rass) / (rass+rcss)
   ! false alarm rate
   IF ((rbss+rdss) > 0._wp) &
    fr = 100.0_wp * (rbss) / (rdss+rbss)
   ! false alarm ratio
   IF ((rass+rbss) > 0._wp) &
    far = 100.0_wp * (rbss) / (rass+rbss)
   ! frequency bias
   IF ((rass+rcss) > 0._wp) &
    fbi = (rass+rbss) / (rass+rcss)
   ! threat score
   IF ((rbss+rcss+rass) > 0._wp) &
    ts  = 100.0_wp * (rass) / (rbss+rcss+rass)
   ! equitable threat score
   ! randomly correct forecasted wet points
   IF ((rass+rbss+rcss+rdss) > 0._wp) &
    rets = ((rass+rbss) * (rass+rcss)) / ((rass+rbss+rcss+rdss))
   IF ((rass+rbss+rcss-rets) > 0._wp) &
    ets = 100.0_wp * (rass-rets) / (rass+rbss+rcss-rets)
   ! True skill statistics or Hanssen-Kuipers discriminant or Kuipers score
   IF ((rass+rcss) > 0._wp .AND. (rbss+rdss) > 0._wp) &
    tss = ( (rass)/(rass+rcss) + (rdss)/(rbss+rdss) - 1._wp) * 100_wp
   ! precipitation area index +
   IF ((rbss+rdss) > 0._wp) &
    pai_plus  = (rcss+rdss) / (rbss+rdss)
   ! precipitation area index -
   IF ((rass+rcss) > 0._wp) &
    pai_minus = (rass+rbss) / (rass+rcss)
   ! Heidke skill score
   IF ((rass+rbss+rcss+rdss) > 0._wp) &
    rhss = (((rass+rbss)*(rass+rcss)+(rcss+rdss)*(rbss+rdss))/(rass+rbss+rcss+rdss))
   IF ((rass+rbss+rcss+rdss-rhss) > 0._wp) &
    hss = 100.0_wp * (rass+rdss-rhss)/(rass+rbss+rcss+rdss-rhss)

   WRITE(nulhn,'(a25,a3,6i7,e12.3)')'skill scores (a,b,c,d):',ytime,ntstep,ass,bss,css,dss,zss,thr(ith)
   IF (lfirst) &
      WRITE(nulhn,*) '#skill scores: ',ytime,' ntstep, hit rate, probability of detection, false alarm ratio, &
    &  false alarm rate, frequency bias, threat score, equitable threat score, heidtke skill score, &
    &  true skill statistics'
   WRITE(nulhn,'(a14,a3,i6,9f8.2)')'skill scores:',ytime,ntstep,hr,pod,far,fr,fbi,ts,ets,hss,tss

  ENDIF

  lfirst=.FALSE.
 ENDDO ! loop over thresholds
 IF (num_compute > 1) THEN
    iexchange = 0
    iexchange( 1)= histmod(7)
    iexchange( 2)= histobs(7)
    iexchange( 3)= anzobs
    iexchange( 4)= anzmod
    CALL global_values (iexchange, 4, 'SUM',imp_integers,icomm_cart, 0,yerrmsg, ierror)
 ENDIF

 IF(my_cart_id == 0) THEN
    IF (num_compute > 1) THEN
       histmod(7) = iexchange( 1)
       histobs(7) = iexchange( 2)
       anzobs = iexchange( 3)
       anzmod = iexchange( 4)
    ENDIF
    DO ith=1,6
       histobs(ith)=histobs(ith)-histobs(ith+1)
       histmod(ith)=histmod(ith)-histmod(ith+1)
    ENDDO
    WRITE(nulhn,'(a17,a3,i6,8i12)')'Histogramm model:',ytime,ntstep,anzmod,(histmod(i),i=1,7)
    WRITE(nulhn,'(a17,a3,i6,8i12)')'Histogramm radar:',ytime,ntstep,anzobs,(histobs(i),i=1,7)
 ENDIF

END SUBROUTINE lhn_skill_scores

!===============================================================================
! End of module
!-------------------------------------------------------------------------------

END MODULE src_lheat_nudge
