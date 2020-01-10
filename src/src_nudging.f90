!+ Source module for organizing the nudging and computing analysis increments
!-------------------------------------------------------------------------------

MODULE src_nudging

!-------------------------------------------------------------------------------
!
! Description:
!   The module "src_nudging" contains the routines which perform the final steps
!   of the data assimilation by nudging towards direct observations.
!   Namely, the analysis increments are computed (using the previous results of
!   the spreading of the observation increments to the target grid points), and
!   explicit balancing steps for the analysis increment fields are optinally 
!   performed.
!   In detail:
!    - The direct nudging analysis increments for horizontal wind, temperature,
!      and specific humidity are computed using the net (spreaded) observation
!      increments and net nudging weights for wind, temperature, and relative
!      humidity at the target grid points as input.
!    - A temperature correction balances the (previously computed) pressure
!      analysis increments at the lowest model level (from nudging surface
!      pressure data) according to pressure - temperature error correlations
!      valid for the mesoscale.
!    - The mass field and the wind field are explicitly coupled (balanced) by
!       - a geostrophic wind correction which balances mass field increments
!         from nudging surface pressure data,
!       - a geostrophic (surface) pressure correction Which balances wind
!         increments from nudging 10-m wind data (e.g. from scatterometers)
!    - Upper-air pressure increments are derived by hydrostatic balancing of
!      the final pressure analysis increments at the lowest model levels
!      and the final temperature and humidity analysis increments at all
!      model levels.
!    - Finally, the prognostic fields are updated by adding the analysis
!      increments (i.e. by execution of the nudging equations).
!
! Method:
!   This module contains the following procedures:
!    - nudge_horiz_wind   : updating the horizontal wind fields by nudging
!    - nudge_humid_mass   : updating the mass fields by nudging
!    - ps_temperatur_corr : temper. correction balancing surface pressure incr.
!    - geostroph_ps_corr  : geostrophic pressure corr. from surfacee wind incr.
!    - poisson_coarse     : core to solve Poisson eq. for 'geostroph_ps_corr'
!   Driving routine is "organize_nudging" of module "src_obs_use_org.f90".
!
!   This module also contains elemental functions, formerly statement functions:
!    - fpvsw     : Magnus formula for water: saturation vapour pressure from T
!    - fpv2q     : specific humidity from water vapour pressure and air pressure
!    - fq2pv     : vapour pressure from specific humidity and air pressure
!    - rmod      : MOD function for positive reals
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
! 1.13       1998/10/22 Christoph Schraff
!  Initial release
! 1.19       1998/12/11 Christoph Schraff
!  Introduction of the verification mode (writing observations to the VOF).
!  Improved diagnostics e.g. on the geostrophic wind correction.
!  ANSI violations removed.
! 1.20       1999/01/07 Guenther Doms
!  Renaming of some global variables
! 1.29       1999/05/11 Ulrich Schaettler
!  Adapted interfaces to utility-modules and prepared use of MPE_IO
! 1.31       1999/07/01 Christoph Schraff
!  quantities related to MPI communicator 'icomm_world' replaced; input to
!  'global_values' adjusted.
! 1.34       1999/12/10 Ulrich Schaettler
!  Changed calls to timing routines
! 1.38       2000/04/06 Christoph Schraff
!  Bug correction when computing wind analysis increments at each timestep.
! 1.39       2000/05/03 Ulrich Schaettler
!  Changed some variable names and adapted call to timing routine get_timings.
! 2.4        2001/01/29 Christoph Schraff
!  Addition of precedure calls for the spatial consistency check for surface
!  pressure data. Preservation of specific instead of relative humidity accord-
!  ing to namelist var. 'khumbal', when nudging temperature. If 'khumbal' > 99,
!  specific humidity is also preserved in the 'temperature correction'.
!  Exchange of boundaries prior to using model values for the nudging.
! 2.5        2001/06/01 Christoph Schraff
!  Re-organisation of 'organize_nudging'. Savety test at array allocations.
! 2.17       2002/05/08 Christoph Schraff
!  Minor bug corrections at temperature correction, effect on increments < 0.2%.
! 2.18       2002/07/16 Ulrich Schaettler
!  Added a test for variable psanai in line 1940f before a division
!  (because of problems on the NEC)
! 2.19       2002/10/24 Christoph Schraff
!  Revision to geostrophic wind correction: height-dependent geostrophic factor,
!  increased upper limit to increments, bug correction (with very small effect).
! 3.5        2003/09/02 Ulrich Schaettler
!  Adapted interface for routine exchg_boundaries
! 3.6        2003/12/11 Christoph Schraff
!  Bug correction for print-out at geostrophic wind correction.
! 3.7        2004/02/18 Ulrich Schaettler
!  Renamed cphi (crlat), acphir (acrlat)
! 3.12       2004/09/15 Christoph Schraff
!  Extension to (prepare to) include assimilation of satellite retrievals.
! 3.13       2004/12/03 Klaus Stephan
!  Modifications to run with latent heat nudging
! 3.18       2006/03/03 Christoph Schraff
!  New option for alternative weighting for multiple observations.
!  Option for separate weighting for different observation types.
!  Introduction of spatial consistency check of IWV from TEMP humidity and GPS.
!  Preparation for use of real-data 1DVar satellite retrievals (MSG, NOAA15-18).
!  Flushing of output files. To allow for reproducibility, replace pp(:,nnew) by 
!  pp(:,nnow) for the scaling of the geostrophic wind increments.
!    and by Klaus Stephan
!  LHN namelist parameter moved to data_lheat_nudge to avoid to many dependencies
!  Introduced temperature increment due to latent heat (Jochen Foerstner)
! 3.21       2006/12/04 Christoph Schraff, Klaus Stephan
!  Time integration of analysis increments introduced for diagnostic purposes.
!  Call get_gs_lheating also for llhnverif
! V3_23        2007/03/30 Ulrich Schaettler
!  Introduced ldump_ascii for flushing the ASCII files
! V4_4         2008/07/16 Ulrich Schaettler
!  Eliminated timing variables which are unused
!  Adapted interface of get_timings
! V4_5         2008/09/10 Christoph Schraff
!  Extension to read observations from NetCDF files instead of an AOF file.
!  Statistics of wind direction analysis increments corrected
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_17        2011/02/24 Ulrich Blahak
!  Adapted interface of exchg_boundaries; corrected kzdims(1:20) -> kzdims(1:24);
!  eliminated my_peri_neigh
! V4_22        2012/01/31 Christoph Schraff
!  - Master routine 'organize_nudging' moved to new module src_obs_use_org.f90.
!  - Routine introduced to derive geostrophically balanced surface pressure
!    increments from surface wind increments by solving the Poisson equation,
!    to improve the assimilation of e.g. scatterometer data
!    (Poisson solver by Heinz-Werner Bitzer).
!  - Latitude dependent reduction of geostrophic wind correction introduced,
!    in order to get reasonably small geostrophic increments near the equator.
!  - Refined option for computation of net increments and weights separately for
!    different sets of observing systems, and determination of final analysis
!    increments by (re-)using the weighting function for multiple increments.
!  - Frequency of opening and closing files (for flushing) reduced.
!  - Some variables moved from module 'data_nudge_all' to 'data_obs_lib_cosmo'.
!  - (Potentially fatal) bug corrected (diagnostic output on geost. wind corr).
! V4_23        2012/05/10 Ulrich Schaettler
!  Adapted call to SR distribute_fields (added sender PE)
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer, Ulrich Blahak, Ulrich Schaettler
!  Replaced qx-variables by using them from the tracer module
!  UB:
!  Implemented internal switch "l2mom_do_extra_satads" to be able
!   to switch on/off the extra saturation adjustments outside the microphysics parts.
!  Introduced nexch_tag for MPI boundary exchange tag to replace ntstep (HJP)
! V4_27        2013/03/19 Christoph Schraff, Astrid Kerkweg, Ulrich Schaettler
!  Header lines in control output 'yustats' adapted to sub-hourly intervals
!  and 14-digit absolute time declaration.
!  Introduced MESSy interface (AK)
! V4_28        2013/07/12 Christoph Schraff, Ulrich Schaettler
!  Bug fix to re-establish (after 4_22) bit reproducibility irrespective of
!  domain decomposition (by extending loops for computing 'zabgeou', 'zabgeov').
!  Statement functions replaced by elemental functions. Improved comments.
!  Use parameters for vertical grid from module vgrid_refatm_utils (US)
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
!  For the COSMO-Model only use vcoord and refatm from vgrid_refatm_utils
! V4_30        2013/11/08 Christoph Schraff
!  Bug fix, to avoid using de-allocated array 'ztcorr' if ntpscor=1.
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
! V5_4         2016-03-10 Christoph Schraff
!  Removal of if statements for AOF interface.
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
    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction

! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the 
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from 
!    the other ones because of the use of the staggered Arakawa-B-grid.
!    
!   zonal direction
    istart,       & ! start index for the forecast of w, t, qv, qc and pp
    iend,         & ! end index for the forecast of w, t, qv, qc and pp
    istartu,      & ! start index for the forecast of u
    iendu,        & ! end index for the forecast of u
    istartv,      & ! start index for the forecast of v
    iendv,        & ! end index for the forecast of v

!   meridional direction
    jstart,       & ! start index for the forecast of w, t, qv, qc and pp
    jend,         & ! end index for the forecast of w, t, qv, qc and pp
    jstartu,      & ! start index for the forecast of u
    jendu,        & ! start index for the forecast of u
    jstartv,      & ! start index for the forecast of v
    jendv,        & ! end index for the forecast of v
    jstartpar,    & ! start index for computations in the parallel program
    jendpar,      & ! end index for computations in the parallel program

! 4. constants for the horizontal rotated grid and related variables
! ------------------------------------------------------------------

    dlon,         & ! grid point distance in zonal direction (in degrees)
    dlat,         & ! grid point distance in meridional direction (in degrees)
    startlat_tot, & ! rotated latitude  \  of the lower left grid point of the
!   startlon_tot, & ! rotated longitude /  total domain (in degrees, N,E>0)
    eddlon,       & ! 1 / dlon
    eddlat,       & ! 1 / dlat
    edadlat,      & ! 1 / (radius of the earth * dlat)
    degrad,       & ! factor for transforming degree to rad

! 4. variables for the time discretization and related variables
! --------------------------------------------------------------

    dt,           & ! long time-step
!   ed2dt,        & ! 1 / (2 * dt)
    dt2,          & ! 2 * dt
    dtdeh,        & ! dt / 3600 seconds

! 8. Organizational variables to handle the COSMO humidity tracers
! ----------------------------------------------------------------
    idt_qv,  idt_qc

! end of data_modelconfig

!-------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 1. mathematical constants
! -------------------------

    pi,           & ! circle constant

! 2. physical constants and related variables
! -------------------------------------------

    g,            & ! acceleration due to gravity
    r_d,          & ! gas constant for dry air
    rdv,          & ! r_d / r_v
    o_m_rdv,      & ! 1 - r_d/r_v
    rvd_m_o,      & ! r_v/r_d - 1
    cp_d,         & ! specific heat of dry air at constant pressure
    cpdr,         & ! 1 / cp_d
    lh_v,         & ! latent heat of vapourization
    r_earth,      & ! mean radius of the earth

! 3. constants for parametrizations
! ---------------------------------

    b1,           & ! variables for computing the saturation vapour pressure
    b2w,          & ! over water (w) and ice (i)
    b3,           & !               -- " --
    b4w,          & !               -- " --
    b234w           ! b2w * (b3 - b4w)

! end of data_constants

!-------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------

    rho0        , & ! reference density at the full model levels      (kg/m3)
    dp0         , & ! reference pressure thickness of layers          ( Pa)
    p0          , & ! reference pressure at main levels               ( Pa)
    hhl         , & ! geometical height of half model levels          ( m )

! 2. external parameter fields                                        (unit)
! ----------------------------

    rlat        , & ! geographical latitude                           ( rad )
    fc          , & ! coriolis-parameter                              ( 1/s )
    crlat       , & ! cosine of transformed latitude
    acrlat      , & ! 1 / ( crlat * radius of the earth )             ( 1/m )

! 3. prognostic variables                                             (unit)
! -----------------------

    u           , & ! zonal wind speed                                ( m/s )
    v           , & ! meridional wind speed                           ( m/s )
!   w           , & ! vertical wind speed (defined on half levels)    ( m/s )
    t           , & ! temperature                                     (  k  )
    pp          , & ! deviation from the reference pressure           ( pa  )

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------

    rho         , & ! total density of air                            (kg/m3)
    tinc_lh     , & ! temperature increment due to latent heat        (  K  )
    qrs         , & ! precipitation water (water loading)             (kg/kg)

! 10. analysis increment fields
! -----------------------------

    ff_anai     , & ! wind velocity                                   ( m/s )
    dd_anai     , & ! wind direction                                  ( rad )
    t_anai      , & ! temperature                                     (  k  )
    p_anai      , & ! deviation from the reference pressure           ( Pa  )
    qv_anai     , & ! specific water vapor content                    (kg/kg)
    qc_anai         ! specific cloud water content (via saturation adjustm)

! end of data_fields

!-------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------

!   nstart,       & ! first time step of the forecast
    nstop,        & ! last time step of the forecast
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
!   nold,         & ! corresponds to ntstep - 1
    nnow,         & ! corresponds to ntstep
    nnew,         & ! corresponds to ntstep + 1

! 7. additional control variables
! -------------------------------

    itype_gscp,   & ! type of grid-scale precipitation physics
    l2tls,        & ! forecast with 2-TL integration scheme
    lcond,        & ! forecast with condensation/evaporation
    ldiabf_satad, & ! include diabatic forcing due to saturation adjustment
    l2mom_satads, & ! in case of 2-moment scheme, do all the satads
                    ! (like for the 1-moment schemes), not just the
                    ! satad after the microphysics at the end of the timestep.
    ltime,        & ! detailled timings of the program are given
    lreproduce,   & ! the results are reproducible in parallel mode
    ldump_ascii,  & ! for flushing (close and re-open) the ASCII files
    lperi_x,      & ! if lartif_data=.TRUE.: periodic boundary conditions in x-dir.
                    ! or with Davies conditions (.FALSE.)
    lperi_y,      & ! if lartif_data=.TRUE.: periodic boundary conditions in y-dir.
                    ! or with Davies conditions (.FALSE.)
    l2dim           ! 2 dimensional runs

! end of data_runcontrol 

!-------------------------------------------------------------------------------


USE data_parallel,      ONLY :  &
    num_compute,   & ! number of compute PEs
    nboundlines,   & ! number of boundary lines of the domain for which
                     ! no forecast is computed = overlapping boundary
                     ! lines of the subdomains

!   ltime_barrier, & ! if .TRUE.: use additional barriers for determining the
                     ! load-imbalance
    ncomm_type,    & ! type of communication
    my_cart_id,    & ! rank of this subdomain in the cartesian communicator
    my_cart_neigh, & ! neighbors of this subdomain in the cartesian grid
    isubpos,       & ! positions of the subdomains in the total domain. Given
                     ! are the i- and the j-indices of the lower left and the
                     ! upper right grid point in the order
                     !                  i_ll, j_ll, i_ur, j_ur.
                     ! Only the interior of the domains are considered, not
                     ! the boundary lines.
    icomm_cart,    & ! communicator for the virtual cartesian topology
    imp_reals,     & ! determines the correct REAL type used in the model
                     ! for MPI
    nexch_tag,     & ! tag to be used for MPI boundary exchange
                     !  (in calls to exchg_boundaries)
    sendbuf,       & ! sending buffer for boundary exchange:
                     ! 1-4 are used for sending, 5-8 are used for receiving
    isendbuflen,   & ! length of one column of sendbuf
    imp_integers     ! determines the correct INTEGER type used for MPI


! end of data_parallel 

!-------------------------------------------------------------------------------

USE data_io,            ONLY :  &

    ydate_ini              ! start of the forecast
                           ! yyyymmddhhmmss (year, month, day, hour, min., sec.)

! end of data_io

!-------------------------------------------------------------------------------

USE data_nudge_all , ONLY :   &

! 1. parameters and related variables
! -----------------------------------

    lwonl        ,& ! .TRUE if grid pt. (ionl ,jonl ) lies in the local domain
    lwonl2       ,& ! .TRUE if grid pt. (ionl2,jonl2) lies in the local domain
    lcroot       ,& ! .TRUE if (my_cart_id == 0) (for print on std. output)
    lfirst       ,& ! .TRUE if 'organize_nudging' is called the first time
    onl_cart_id  ,& ! 'my_cart_id' of node with area containing (ionl,jonl)

! 2. namelist variables controlling the data assimilation
! -------------------------------------------------------

    nwtyp        ,& ! 1   : if > 1 then compute net obs. increments for 'nwtyp'
                    !       different sets of observing systems separately
    kwtyp        ,& ! 1   : function for weights W for multiple observations
    nudgend      ,& ! 0   : end of nudging period in timesteps
    nverend      ,& ! 0   : end of verification period in timesteps
    tconbox      ,& ! 6*dt: timestep [s] for computing analysis increments
                    !       (i.e. time box of constant analysis increments)
    luvgcor      ,& ! .t. : .t. ==> geostrophic wind correction applied
    mpsgcor      ,& ! 1   : mode to apply geostrophic pressure correction
    ntpscor      ,& ! 1   : switch for temperature correction when the 'surface'
                    !       pressure (i.e. p. at lowest model level) is nudged
    khumbal      ,& ! 100 : range around convectively precipitating grid pts, at
                    !       which specified (not relative) humidity is preserved
    ptpstop      ,& ! 400.: pressure [hPa] at upper boundary of the temperature
                    !       correction for 'surface' pressure nudging
    qgeo         ,& ! .3  : factor to the geostrophic wind increments at 1000hPa
    qgeotop      ,& ! .5  : factor to the geostrophic wind increments at ptpstop
    qgeops       ,& ! 1.  : factor to the geostrophic pressure increments
    gnudg        ,& ! 6,12,6,6*10^-4: nudging coefficients for TEMP / PILOT data
    ionl         ,& ! 167    : / grid point coordinates
    jonl         ,& ! 103    : \ for standard output on nudging
    ionl2        ,& ! 167    : / 2nd grid pt coordinates
    jonl2        ,& ! 103    : \ for other standard output on nudging

! 6. Miscellany
! -------------

    zconbas      ,& ! height of base of convectively instable region
    zcontop         ! height of top  of convectively instable region

! end of data_nudge_all

!-------------------------------------------------------------------------------

USE data_obs_lib_cosmo , ONLY :   & 

! 1. General parameters
! ---------------------

    c0         ,& ! standard real constant 0.0
    c1         ,& ! standard real constant 1.0
    c2         ,& ! standard real constant 2.0
    c05        ,& ! standard real constant 0.5
    epsy       ,& ! = 1.E-8_wp : commonly used very small value > 0

! 4. I/O device numbers for obs processing / nudging and file names
! -----------------------------------------------------------------

    yustats    ,& ! file name for statistics of processed reports
    yuprint    ,& ! file name for all the remaining information
    nustat     ,& ! unit number of file for statistics of processed reports
    nupr          ! unit number of file for all the remaining information

! end of data_obs_lib_cosmo

!-------------------------------------------------------------------------------

USE data_nudge_spread , ONLY :   &

! 1. Analysis increment fields (to be kept in the long-term storage)
! ------------------------------------------------------------------

    uanai        ,& ! analysis increments of zonal wind
    vanai        ,& ! analysis increments of meridional wind
    tanai        ,& ! analysis increments of temperature
    qanai        ,& ! analysis increments of specific humidity
    psanai       ,& ! analysis increments of pressure at the lowest model level
                    !   from nudging of surface pressure observations
    psaigeo      ,& ! analysis increments of pressure at the lowest model level
                    !   from geostrophic balancing of wind analysis increments
    taips        ,& ! temperature analysis increments due to pressure nudging

! 2. Analysis increment fields
! ----------------------------

    zpai         ,& ! analysis incr. of pressure, condens./evapor. not included
    zroi         ,& ! analysis incr. of density , condens./evapor. not included
    ztwips       ,& ! temperature increments implied (physic'ly or statistic.)
                    ! by pressure analysis increments at the lowest model level
    ztpgeo       ,& ! temperature increments implied (physic'ly or statistic.)
                    ! by geostrophic (surface) pressure analysis increments

! 3. Output of spreading procedures
! ---------------------------------

    omy          ,& ! sum of  spatial * temporal * quality 'spreading weights'
    om2          ,& ! sum of  squares of 'spreading weights'
    zwi          ,& ! sum of weighted (observation) increments (weights are
                    !        squares of 'spreading weights')

! 4. Further fields used to compute analysis increments
! -----------------------------------------------------

    faclbc       ,& ! reduction of nudging weights (near lateral boundaries)

! 10. Other variables
! -------------------

    wvgeo        ,& ! vertical weights to the geostrophic wind correction
    nhisto          ! data density histogram for surface pressure obs 

! end of data_nudge_spread

!-------------------------------------------------------------------------------

USE data_lheat_nudge,    ONLY :  &
    llhn,         & ! main switch for latent heat nudging
    llhnverif       ! verification for latent heat nudging

!-------------------------------------------------------------------------------

USE environment,              ONLY :  &
      exchg_boundaries, & ! performs the boundary exchange between 
                          ! neighboring processors
      comm_barrier,     & ! explicit synchronization point
      model_abort         ! aborts the program in case of errors

!-------------------------------------------------------------------------------

USE parallel_utilities,       ONLY :  &
      global_values    ,& ! computes global values by operating on local arrays
      gather_field     ,& ! gathers the parts of a total field from all
                          !   subdomains
      distribute_field    ! distributes the parts of a total field to all
                          !   subdomains

!-------------------------------------------------------------------------------

USE utilities,                ONLY :  &
      uv2df               ! converts wind components into direction and speed
!     uv2df_vec           ! converts wind components into direction and speed

!-------------------------------------------------------------------------------

USE pp_utilities,             ONLY :  &
      calpmsl             ! computes the mean sea-level pressusre

!-------------------------------------------------------------------------------

USE meteo_utilities,          ONLY :  &
!     calps,            & ! surface pressure
      satad               ! saturation adjustment

!-------------------------------------------------------------------------------

USE time_utilities,           ONLY :  &
      get_timings,   i_barrier_waiting_nud,  &
      i_communications_nud,   i_nudging_corr

!-------------------------------------------------------------------------------

USE vgrid_refatm_utils, ONLY :   &
    vcoord

!-------------------------------------------------------------------------------

USE src_lheating,           ONLY:    &
  get_gs_lheating          ! storage of grid scale latent heating for lhn

!-------------------------------------------------------------------------------

USE src_tracer,             ONLY: trcr_get, trcr_errorstr

!-------------------------------------------------------------------------------
! Local declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE

! 1. Variables
! ------------

  INTEGER (KIND = iintegers)              , PRIVATE :: &
    ktoptps           ! uppermost full level influenced by the temperature
                      ! correction

  REAL    (KIND = wp)                     , PRIVATE :: &
    zpaionl        ,& ! analysis increment of (upper-air) pressure (hydrostati-
                      ! cally derived) at diagnostic grid pt (for diagn. output)
    ztaionl        ,& ! analysis increm. of temperature (incl. condens./evapor.)
    zqaionl        ,& ! analysis increm. of humidity (in the form of the density
                      ! increment due to moisture, condens./evapor. included)
    zrealdiff         ! elapsed time since latest measuring

! For error handling
! ------------------
  INTEGER (KIND=iintegers) ::  izerror

  CHARACTER (LEN=80)       ::  yzerrmsg

!===============================================================================

!-------------------------------------------------------------------------------
! Public and Private Subroutines
!-------------------------------------------------------------------------------

CONTAINS


!-------------------------------------------------------------------------------
!+ Module procedure for updating the horizontal wind fields by nudging
!-------------------------------------------------------------------------------

SUBROUTINE nudge_horiz_wind ( lgetai , lconai , k )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure updates the horizontal wind fields by execution of
!   nudging equations.
!   Optionally, the update includes geostrophic wind increments, which are
!   deduced here from the analysis increments of pressure and density resulting
!   from nudging the pressure at the lowest model layer and applying the
!   temperature correction.
!   The resulting total analysis increments of horizontal wind are stored if
!   identical analysis increments are used during several timesteps.
!
! Method:
!   1.  For the horizontal wind components, 'net (spreaded) observation incre-
!       ments' (squared weighted mean of observation increments) and 'net
!       nudging weights' (weighted mean of nudging weights) at the (center(!) of
!       the) target grid points are derived from the sum of weights, the sum of
!       squared weights, and the sum of squared weighted observation increments
!       as available from the spreading procedures. This derivation of weighting
!       of multiple data is based on Benjamin & Seaman (MWR, 1985) and accounts
!       for data density at the model grid points.
!   2a. If the program is run in parallel mode, one row of these net observation
!       increments and net nudging weights are exchanged at the boundaries of
!       the sub-domains between the nodes. If geostrophic wind increments are
!       also to be computed, this exchange includes the analysis increments of
!       pressure and density (as from nudging the pressure at the lowest model
!       layer and applying the temperature correction).
!   2b. The net increments and weights are shifted from the center of the
!       Arakawa-C grid points to the locations where the model wind components
!       are defined, and the weights are reduced near the lateral boundaries of
!       the total domain.
!   3a. Geostrophic wind increments are computed from the analysis increments
!       of pressure and density. These increments are determined prior to the
!       latent heat / sensible heat transfer (cf. proc. *nudge_humid_mass*) to
!       prevent unrealistic small-scale wind increments in saturated regions,
!       yet include a change in specific humidity to keep the relative humidity
!       constant.
!   3b. The geostrophic wind increments are reduced depending on parameters.
!       The reduction is sub-divided in a global reduction (to account for the
!       limited validity of the geostrophic approximation for the scales of
!       interest), a reduction near the lateral boudaries (where pressure and
!       density analysis increments are reduced, resulting in unrealistic hori-
!       zontal gradients), and a reduction depending on the vertical model level
!       (to account for strong ageostophy near the ground, and a large impact
!       of biased temperature data to geostrophic wind in the upper atmosphere).
!       Furthermore, a general upper limit is imposed, since the geostrophic
!       wind correction (integrated in time) may get too large, where nudging
!       increases geopotential gradients in a curved, significantly ageostrophic
!       flow (e.g. at mesoscale cyclones).
!       (Note this information again needs to be exchanged between nodes here.)
!   4.  Analysis increments for the horizontal wind components are computed from
!       the shifted 'net (spreaded) observation increments', the 'net nudging
!       weights', and the geostrophic wind increments, and are added to the
!       model fields.
!       (--> "Nudging of horizontal wind")
!   The current procedure is called by proc. 'organize_nudging' within a verti-
!   cal loop over model levels. However, if the program is run in parallel mode,
!   or (and) if geostrophic wind increments are to be computed, only step '1.'
!   is applied, except when the vertical loop has arrived at its last level 
!   (i.e. when all the spreading and mass field nudging has already been com-
!   pleted at the call of the current procedure). Then, steps 2 to 4 are applied
!   consecutively for all vertical model levels in a loop on its own.
!   The reasons are that firstly in the parallel mode, this allows each node to 
!   be independent during all the (expensive) spreading and the mass field
!   nudging, and communication between nodes (apart from the exchange of the
!   local observational information prior to the spreading) is only needed at
!   the very end of the nudging scheme for steps 2 to 4. And secondly for the
!   computation of geostrophic wind increments, pressure analysis increments
!   need not only be known at the current model level, but also at the adjacent
!   model levels.
!
! Written by        :  Christoph Schraff, DWD  (original version: 10.06.97)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
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

  LOGICAL                 , INTENT (IN)         ::       &
    lgetai           ,& ! .TRUE if new analysis increment are to be computed
                        !          and used at the current timestep
    lconai              ! TRUE if analysis increments are const during 'tconbox'

  INTEGER (KIND=iintegers), INTENT (IN)         ::       &
    k                   ! index of current vertical model level

! Local parameters:
! ----------------

  REAL    (KIND=wp   )     , PARAMETER ::  &
    pmaxgeo= 20.0_wp ,& ! upper limit to abs. of geostrophic wind increment
                            ! vector [m/s] per hour. The value of this limit is
                            ! imposed to each increment vector at each timestep
    plbcno = 10.0_wp ,& ! Number of grid rows at the lateral boundaries of
                            ! the total domain with zero weight for the geo-
                            ! strophic wind increments
    plbcipr= 0.125_wp,& ! Number of grid rows within which the horizontal
                            ! weight to these increments increase linearly from
                            ! zero at the lateral boundaries to one in the inner
                            ! part of the total domain  (reciproke of idem)
    ddfflim=  3.0_wp    ! wind speed limit below which analysis increments
                            ! of wind direction are not added to time sum

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    izerror          ,& ! error index
    i, j, kzdims(24) ,& ! loop indices in horizontal direction
    kuwi             ,& ! vertical index for net increments and weights
    kwi              ,& ! vertical index for analysis increments
    kv               ,& ! loop index for inner vertical loop
    kva    , kve     ,& ! range of loop index for inner vertical loop
    jstadz , jenddz  ,& ! horiz. index range used for the mass weighted Coriolis
    istadz , ienddz  ,& !   parameter, or for vertical gradients
    jstaum1, istavm1 ,& ! jstartu-1 , istartv-1
    kp1    , km1     ,& ! indices of model levels used for vertical gradients
    ityp             ,& ! index of set of obs system in weighted increment array
    ixuv             ,& ! index for grid pt. location for control output
    ix     , ixe     ,& ! indices used for control output
    kk               ,& ! vertical loop index
    ipr(2) , jpr(2)  ,& ! loop indices in horizontal direction for printing
    ilpr             ,& ! loop index
    navai            ,& ! number of grid points in (sub-)domain
    nulast           ,& ! timestep at which statistics is printed
    kdpr   , kepr    ,& ! indices defining printed vertical levels
    nzerr  , nstat      ! indicator of status for (de-)allocation of fields

  REAL    (KIND=wp   )     ::  &
    fmultw           ,& ! 0. or 1., depending on multiple observation weighting
    zdt              ,& ! timestep relevant for the nudging
    zgnudg           ,& ! basic nudging coefficient for horizontal wind
    zdt2m            ,& ! arithmetic weight (0.5) times twice the timestep
    zwvgeo           ,& ! vertically dependent weight to geostrophic increments
    zrtv             ,& ! R * virtual temperature (time level 'nnow')
    z2rhomr          ,& ! half the inverse of the horizonally averaged density
    zprefac          ,& ! factor to the horizontal pressure gradient
    zfppdz           ,& ! factor to the vertical pressure gradient (term) as
                        ! part of the horizontal pressure gradient
    zabsgeo          ,& ! square of geostrophic wind increment vector
    zalatu  ,zalatv  ,& ! ABS( true geogr. latitude of Arakawa-C u-/v-grid pt. )
    ziec_tot,zjec_tot,& ! location of center of domain in global coord.
    ziec    ,zjec    ,& ! location of center of domain in local  coord.
    zlbciy  ,zlbcjy  ,& ! distance from center with geost. horiz. weight = 0.5
    zwlbc            ,& ! horizontal weight for geostrophic wind increments
    zmaxgeo2         ,& ! (upper limit to abs. of geost. wind incr. vector) **2
    zfbuoyt          ,& ! factor in the buoyancy term in the w-equation
    zfbuoyb          ,& ! as 'zfbuoyt' but for for a model level further below
    zdp01dt ,zdp0dt1 ,& ! other factors in the buoyancy term in the w-equation
    zqgeo            ,& ! pressure dependent grade of geostrophy
    zfqgeo           ,& ! interpol. factor for pressure dep. grade of geostrophy
    zpgref           ,& ! reference pressure, at which 'qgeo' is valid (1000hPa)
    zdpgrt           ,& ! pressure difference related to 'qgeo', 'qgeotop'
    zpkemin          ,& ! for printing: minimum pressure in local domain
    zp1000           ,& ! for printing: for p(ke) == zp1000 then ps = 1000 hPa
    zdp1000          ,& ! for printing: minimum departure from 'zp1000'
    zmass            ,& ! for printing: mass of base state
    zdanai              ! analysis incr. of wind direction (-180 < zdanai<= 180)

  REAL    (KIND=wp   )     ::  &
    zumn             ,& ! zonal  comp. of wind after  adding analysis increments
    zumo             ,& ! zonal  comp. of wind before adding analysis increments
    zvmn             ,& ! merid. comp. of wind after  adding analysis increments
    zvmo             ,& ! merid. comp. of wind before adding analysis increments
    zddn             ,& ! wind direction after  adding analysis increments
    zddo             ,& ! wind direction before adding analysis increments
    zffn             ,& ! wind speed     after  adding analysis increments
    zffo                ! wind speed     before adding analysis increments

  LOGICAL                  ::  &
    lgetgeo          ,& ! geostrophic wind correct. is computed at current time
    lprgeo           ,& ! printing for control
    lwrite              ! printout for control

  CHARACTER (LEN=11)       ::  yuvgeo
  CHARACTER (LEN=5 )       ::  ychh     ! hour of integration time
  CHARACTER (LEN=255)      ::  yzerrmsg
  CHARACTER (LEN=25)       ::  yzroutine

! Local (automatic) arrays:
! -------------------------

  REAL    (KIND=wp   )     , ALLOCATABLE , SAVE ::       &
    zuwi    (:,:,:)  ,& ! net (spreaded) observation increment of zonal wind
    zvwi    (:,:,:)  ,& ! net (spreaded) observation increment of merid. wind
    zumy    (:,:,:)     ! net nudging weight for horizontal wind

  REAL    (KIND=wp   )     , ALLOCATABLE , SAVE ::  &
    ztpsai      (:)  ,& ! temperature increment from temperature correction
    ztnuai      (:)  ,& ! temperature incr. from temp. corr. & simple nudging
    ztecai      (:)  ,& ! temperature incr. including condensation & evaporation
    zqnuai      (:)  ,& ! specific humidity increment from simple nudging
    zqecai      (:)  ,& ! specific humidity incr. incl. condensat. & evaporation
    zpecai      (:)     ! pressure analysis increment after complete nudging
 
  REAL    (KIND=wp   )     ::  &
    zwsum   (ie,je,2),& ! sum of weights assigned to the different obs. types
    zuwipos (ie,je)  ,& ! net obs. increment of zonal wind shifted to u-grid pt.
    zvwipos (ie,je)  ,& ! net obs. increment of meri. wind shifted to v-grid pt.
    zumydt2 (ie,je)  ,& ! 2*timestep times net nudging weight shifted to u-g.p.
    zvmydt2 (ie,je)  ,& ! 2*timestep times net nudging weight shifted to v-g.p.
    zuwigeo (ie,je)  ,& ! zonal  comp. of geostrophic wind increments
    zvwigeo (ie,je)  ,& ! merid. comp. of geostrophic wind increments
    zuwigef (ie,je)  ,& ! zonal  comp. of full geostrophic wind increments
    zvwigef (ie,je)  ,& ! merid. comp. of full geostrophic wind increments
    zabgeou (ie,je)  ,& ! square of geostrophic wind increment vector at u-g.pt.
    zabgeov (ie,je)  ,& ! square of geostrophic wind increment vector at v-g.pt.
    zfcmw   (ie,je)  ,& ! mass weighted Coriolis parameter
    zrho    (ie,je)  ,& ! density (time level 'nnow')
    zppdz   (ie,je)  ,& ! vertical gradient of the pressure perturbation (*2)
    zpaidz  (ie,je)  ,& ! vertical gradient of the pressure analysis incr. (*2)
    zdp0mz  (ie,je)  ,& ! vertical gradient of the base state pressure (*2)
    zdpgxrho(ie,je)  ,& ! perturbation of west-east pressure gradient due to
                        ! pressure perturbations
    zpgxdrho(ie,je)  ,& ! perturbation of west-east pressure gradient due to
                        ! density  perturbations
    zdpgyrho(ie,je)  ,& ! perturbation of north-south pressure gradient due to
                        ! pressure perturbations
    zpgydrho(ie,je)  ,& ! perturbation of north-south pressure gradient due to
                        ! density  perturbations
    zpp     (ie,je)  ,& ! pressure at mass grid points
    zppu    (ie,je)  ,& ! pressure at zonal wind grid points
    zppv    (ie,je)  ,& ! pressure at meridional wind grid points
    zedacdl (   je)  ,& ! geometric factor to horizontal pressure gradient
    zqcrlat (   je)  ,& ! averaging factor to merid. geostrophic wind
    zuvmx   (    4)     ! max. geostrophic increments at u-, v- grid pts.

  INTEGER (KIND=iintegers) ::  &
    ijuvmx  ( 2, 4)     ! location with max. geostrophic incr. at u-, v- grid p.

  REAL    (KIND=wp   )     ::  &
    avai_ff (   ke)  ,& ! horiz. average of wind speed     analysis increments
    avai_dd (   ke)  ,& ! horiz. average of wind direction analysis increments
    avai_t  (   ke)  ,& ! horiz. average of temperature    analysis increments
    avai_p  (   ke)  ,& ! horiz. average of pressure       analysis increments
    avai_fi (   ke)  ,& ! horiz. average of geopotential   analysis increments
    avai_qv (   ke)  ,& ! horiz. average of specific humidity   analysis increm.
    avai_qc (   ke)  ,& ! horiz. average of cloud water content analysis increm.
    avai_tqv         ,& ! horiz. average of integrated vapour analysis increm.
    avai_tqc         ,& ! horiz. average of total cloud water analysis increm.
    avai_dse         ,& ! domain average of dry static energy analysis increm.
    zugeai  (   ke)  ,& ! zonal  comp. of geostrophic wind incr. (for printout)
    zvgeai  (   ke)  ,& ! merid. comp. of geostrophic wind incr. (for printout)
    zpnoai  (   ke)  ,& ! pressure analysis incr, model temperature unchanged
    zppsai  (   ke)  ,& ! pressure analysis incr, temperature correction applied
    zpnuai  (   ke)  ,& ! pressure analysis incr, temper. corr & nudging applied
    zpnew   (   ke)     ! pressure

  CHARACTER (LEN=42)       ::  ypr(2)       ! output for control

  REAL    (KIND=wp   )     , ALLOCATABLE , SAVE ::       &
    zxmax       (:)     ! max. increments of horizontal wind components

  INTEGER (KIND=iintegers) , ALLOCATABLE , SAVE ::       &
    ijxmx     (:,:)     ! coord. with max. increments of horizontal wind compon.

! Tracer pointers:
! -----------------
  REAL (KIND=wp)     , POINTER  :: &
    qv(:,:,:) => NULL()       ,&   ! QV at nnow
    qc(:,:,:) => NULL()            ! QC at nnow

!
!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine nudge_horiz_wind
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1: Preliminaries (incl. allocation of memory)
!-------------------------------------------------------------------------------

  izerror   = 0_iintegers
  yzerrmsg  = ''
  yzroutine = 'nudge_horiz_wind'

  IF (l2tls) THEN
    zdt = dt
  ELSE
    zdt = dt2
  ENDIF
  kzdims(:) = 0_iintegers
 
  lgetgeo  =  (luvgcor) .AND. (lgetai)

  IF ( (num_compute > 1) .OR. (lgetgeo) ) THEN
    kuwi = k
    IF (k == ke) THEN
      IF (ALLOCATED( zuwi )) THEN
        PRINT '("CAUTION in src_nudging: zuwi is already allocated "           &
              &,"at time / PE",I6,I5)', ntstep, my_cart_id
        DEALLOCATE ( zuwi , zvwi , zumy , STAT=nzerr )
      ENDIF
      ALLOCATE ( zuwi(ie,je,ke) , stat=nzerr )
      ALLOCATE ( zvwi(ie,je,ke) , stat=nzerr )
      ALLOCATE ( zumy(ie,je,ke) , stat=nzerr )
    ENDIF
  ELSE
    kuwi = 1
    IF (ALLOCATED( zuwi )) THEN
      PRINT '("CAUTION in src_nudging: zuwi is already allocated "             &
            &,"at time / PE",I6,I5)', ntstep, my_cart_id
      DEALLOCATE ( zuwi , zvwi , zumy , STAT=nzerr )
    ENDIF
    ALLOCATE ( zuwi(ie,je,1) , stat=nzerr )
    ALLOCATE ( zvwi(ie,je,1) , stat=nzerr )
    ALLOCATE ( zumy(ie,je,1) , stat=nzerr )
  ENDIF

  ! retrieve the required microphysics tracers (at nnow)
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnow, ptr = qv)
  IF (izerror /= 0_iintegers) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nnow, ptr = qc)
  IF (izerror /= 0_iintegers) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

!-------------------------------------------------------------------------------
!  Section 2a: Computation of a
!              1):    squared weighted mean   =    net (spreaded) obs. increment
!                     of obs. increments           at the target grid pt.
!              2):    weighted mean           =    net nudging weight
!                     of nudging weights           at the target grid pt.
!              Note that these quantities are still defined at the center of
!              the Arakawa-C grid points.
!-------------------------------------------------------------------------------
 
  IF (lgetai) THEN
    lwrite  =  (lwonl) .AND. (k == ke-5) .AND. (ntstep <= 5)
    IF (lwrite) WRITE( nupr,'(''omu: weight/ sqr/ weighted incr'',I3,4F9.4)' ) &
                       k, omy(ionl,jonl,1,1), om2(ionl,jonl,1,1)               &
                        , zwi(ionl,jonl,1,1), zwi(ionl,jonl,2,1) 
! Inclusion of nudging coefficient (the MAX(.,epsy) is required for nudging of
! surface-level data, if the nudging coeff. for TEMP / PILOT is set to zero)
    zgnudg = MAX( gnudg(1) , epsy )

    DO   j = jstart, jend
      DO i = istart, iend
        zwsum (i,j,1)     =  c0
        zwsum (i,j,2)     =  c0
        zuwi  (i,j,kuwi)  =  c0
        zvwi  (i,j,kuwi)  =  c0
      ENDDO
    ENDDO
    DO ityp = 1 , nwtyp
      fmultw = REAL( kwtyp(ityp) - 1, wp )
      DO   j = jstart, jend
        DO i = istart, iend
          IF ((om2(i,j,1,ityp) > epsy) .AND. (omy(i,j,1,ityp) > epsy)) THEN
! averaged increment for observation type 'ityp'
            zwi (i,j,1,ityp)  =   zwi(i,j,1,ityp)                              &
                               / (om2(i,j,1,ityp) + fmultw *omy(i,j,1,ityp))
            zwi (i,j,2,ityp)  =   zwi(i,j,2,ityp)                              &
                               / (om2(i,j,1,ityp) + fmultw *omy(i,j,1,ityp))
! individual weight 'w(ityp)' assigned to the observation type 'ityp'
            omy (i,j,4,ityp)  =  (om2(i,j,1,ityp) + fmultw *omy(i,j,1,ityp))   &
                               / (omy(i,j,1,ityp) + fmultw)
! sum of the individual weights 'w(ityp)' assigned to the different obs. types
            zwsum (i,j,1)     =  omy(i,j,4,ityp)  +  zwsum(i,j,1)
          ELSE
            zwi (i,j,1,ityp)  =  c0
            zwi (i,j,2,ityp)  =  c0
            omy (i,j,4,ityp)  =  c0
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    fmultw = REAL( kwtyp(1) - 1, wp )
    IF (nwtyp >= 2)  fmultw = REAL( kwtyp(nwtyp+1) - 1, wp )
    DO ityp = 1 , nwtyp
      DO   j = jstart, jend
        DO i = istart, iend
          IF (zwsum(i,j,1) > epsy) THEN
! total multiple weight 'W(ityp)' assigned to observation type 'ityp'
            omy (i,j,1,ityp)  =  omy(i,j,4,ityp) * (omy(i,j,4,ityp) + fmultw)  &
                                                 / (zwsum(i,j,1)    + fmultw)
! weighted averaged increment for observation type 'ityp'
            zwi (i,j,1,ityp)  =  zwi(i,j,1,ityp) * omy(i,j,1,ityp)
            zwi (i,j,2,ityp)  =  zwi(i,j,2,ityp) * omy(i,j,1,ityp)
! sum of total multiple weights 'W(ityp)' assigned to the different obs. types
            zwsum (i,j,2)     =  omy(i,j,1,ityp)  +  zwsum(i,j,2)
          ELSE
            omy (i,j,1,ityp)  =  c0
            zwi (i,j,1,ityp)  =  c0
            zwi (i,j,2,ityp)  =  c0
          ENDIF
! sum of the weighted averaged increments for the different observation types
          zuwi (i,j,kuwi)  =  zwi(i,j,1,ityp)  +  zuwi(i,j,kuwi)
          zvwi (i,j,kuwi)  =  zwi(i,j,2,ityp)  +  zvwi(i,j,kuwi)
        ENDDO
      ENDDO
    ENDDO
    DO   j = jstart, jend
      DO i = istart, iend
        IF (zwsum(i,j,2) > epsy) THEN
! double-averaged total observation increment
          zuwi (i,j,kuwi)  =  zuwi(i,j,kuwi) / zwsum(i,j,2)
          zvwi (i,j,kuwi)  =  zvwi(i,j,kuwi) / zwsum(i,j,2)
! final weight, including nudging coefficient
          zumy (i,j,kuwi)  =  zwsum(i,j,2)  * zgnudg
        ELSE
          zuwi (i,j,kuwi)  =  c0
          zvwi (i,j,kuwi)  =  c0
          zumy (i,j,kuwi)  =  c0
        ENDIF
      ENDDO
    ENDDO

!   DO   j = jstart, jend
!     DO i = istart, iend
!       IF (om2(i,j,1) > epsy) THEN
!         zuwi (i,j,kuwi) = zwi(i,j,1) / om2(i,j,1)
!         zvwi (i,j,kuwi) = zwi(i,j,2) / om2(i,j,1)
!         zumy (i,j,kuwi) = om2(i,j,1) / omy(i,j,1) * zgnudg
!       ELSE
!         zuwi (i,j,kuwi) = c0
!         zvwi (i,j,kuwi) = c0
!         zumy (i,j,kuwi) = c0
!       ENDIF
!     ENDDO
!   ENDDO
    IF (lwrite) WRITE( nupr,'(''omu: weight     / weighted incr'',I3,F12.7,6X  &
                            &,2F9.4)' )   k, zumy(ionl,jonl,kuwi)              &
                     , zuwi(ionl,jonl,kuwi), zvwi(ionl,jonl,kuwi) 
  ENDIF

!-------------------------------------------------------------------------------
! Section 3: Exchange of 1 row of boundary data (of local domains)
!-------------------------------------------------------------------------------

IF ( (num_compute > 1) .AND. (lgetai) .AND. (k == 1) ) THEN
  IF (ltime) THEN
    CALL get_timings (i_nudging_corr, ntstep, dt, izerror)
    CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
    CALL get_timings (i_barrier_waiting_nud, ntstep, dt, izerror)
  ENDIF

  IF (lgetgeo) THEN

    kzdims(1:24)=(/ke,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)

    CALL exchg_boundaries ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart,       &
                           num_compute, ie, je, kzdims, jstartpar, jendpar, 1,    &
                           nboundlines, my_cart_neigh, lperi_x, lperi_y, l2dim,   &
                           20000+nexch_tag,.FALSE., ncomm_type, izerror, yzerrmsg,&
                           zuwi(:,:,:), zvwi(:,:,:), zumy(:,:,:),                 &
                           zpai(:,:,:), zroi(:,:,:) )
!   =====================

  ELSE

    kzdims(1:24)=(/ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)

    CALL exchg_boundaries ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart,       &
                           num_compute, ie, je, kzdims, jstartpar, jendpar, 1,    &
                           nboundlines, my_cart_neigh, lperi_x, lperi_y, l2dim,   &
                           20000+nexch_tag,.FALSE., ncomm_type, izerror, yzerrmsg,&
                           zuwi(:,:,:), zvwi(:,:,:), zumy(:,:,:) )
!   =====================

  ENDIF

  IF (ltime) CALL get_timings (i_communications_nud, ntstep, dt, izerror)
ENDIF

!-------------------------------------------------------------------------------
! Section 4: Setting up of vertical loop over model levels, and
!            computation (shift) of net increments and nudging weights at the
!            u-grid pts. and v-grid pts. of the Arakawa-C grid respectively
!-------------------------------------------------------------------------------

  IF ( (num_compute > 1) .OR. (lgetgeo) ) THEN
! Wind nudging is executed only once (if (k == 1)), but for all levels
! ==> no vertical loop except for (k == 1) when loop reaches from ke to 1
    kva  = 0
    kve  = 1
    IF (k == 1) kva = ke
  ELSE
! Wind nudging is executed each time for the one level k
! ==> vertical loop is only over one model level k
    kva  = k
    kve  = k
    kuwi = 1
  ENDIF

loop_over_model_levels:  DO kv = kva , kve , -1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IF (lgetai) THEN
    IF ( (num_compute > 1) .OR. (lgetgeo) ) kuwi = kv
    zdt2m  = c05 * zdt

    DO   j = jstartu, jendu
      DO i = istartu, iendu
        zuwipos (i,j) =  (zuwi(i,j,kuwi) + zuwi(i+1,j,kuwi)) * c05
! Inclusion of twice the timestep, & of weight reduction near lateral boundary
        zumydt2 (i,j) =  (zumy(i,j,kuwi) + zumy(i+1,j,kuwi)) * zdt2m           &
                       * faclbc(i,j,2)
      ENDDO
    ENDDO
    DO   j = jstartv, jendv
      DO i = istartv, iendv
        zvwipos (i,j) =  (zvwi(i,j,kuwi) + zvwi(i,j+1,kuwi)) * c05
        zvmydt2 (i,j) =  (zumy(i,j,kuwi) + zumy(i,j+1,kuwi)) * zdt2m           &
                       * faclbc(i,j,3)
      ENDDO
    ENDDO
  ENDIF


!-------------------------------------------------------------------------------
! Section 5: Computation of (full) geostrophic wind increments
!-------------------------------------------------------------------------------

  zwvgeo = c0

! vertical weight (depending on model layer only) to geostrophic wind correction

  IF (lgetgeo) THEN
    IF (kv >= ktoptps) zwvgeo = wvgeo(kv)
    zpgref  = 100000.0_wp
    zdpgrt  = c1 / (zpgref - ptpstop)
  ENDIF

  IF (zwvgeo > epsy) THEN

! --> need values within   (istart-1,iend+1) , (jstart-1,jend+1) , (kv-1,kv+1)
!     of zroi, zpai, pp(:,:,:,nnow), p0, zrho
! (Note that 'zrho' is defined in *slow_tendencies* within 
!                   (istart-1,iend+1) , (jstart-1,jend+1) , (1,ke)
!          , 'pp(,,,nnow), p0 are defined on the total domain)

! preset geostrophic increments notably at the boundaries of the total domain

    DO    j = 1 , je
      DO  i = 1 , ie
        zuwigeo (i,j) = c0
        zvwigeo (i,j) = c0
        zuwigef (i,j) = c0
        zvwigef (i,j) = c0
      ENDDO
    ENDDO

    DO    j = jstartv, jendv
      zedacdl (j) = eddlon * acrlat(j,1)
      zqcrlat (j) = c05 *(crlat(j,1) + crlat(j+1,1)) / crlat(j,2)
    ENDDO
    zedacdl (jendv+1) = eddlon * acrlat(jendv+1,1)

! Coriolis parameter per mass, at the corners of Arakawa-C grid pts.

    jstaum1 = jstartu-1
    istavm1 = istartv-1
    jstadz  = MIN( jstaum1 , jstartv )
    jenddz  = MAX( jendu   , jendv   )
    istadz  = MIN( istartu , istavm1 )
    ienddz  = MAX( iendu   , iendv   )
    DO    j = jstadz, jenddz
      DO  i = istadz, ienddz
        zfcmw (i,j) = fc(i,j) * 4.0_wp / ( dp0(i  ,j,kv) + dp0(i  ,j+1,kv) &
                                             + dp0(i+1,j,kv) + dp0(i+1,j+1,kv))
      ENDDO
    ENDDO

! density and various vertical gradients (*2)
! (used for the horizontal pressure gradient perturbation)

    kp1     = MIN( ke, INT( kv+1 ,iintegers) )
    km1     = MAX( 1 ,      kv-1 )
    jstadz  = MIN( jstaum1   , jstartv   )
    jenddz  = MAX( jendu  +1 , jendv  +1 )
    istadz  = MIN( istartu   , istavm1   )
    ienddz  = MAX( iendu  +1 , iendv  +1 )
    DO    j = jstadz, jenddz
      DO  i = istadz, ienddz
        zrtv         = r_d * t(i,j,kv,nnow) *( c1 + rvd_m_o*qv(i,j,kv)   &
                                           -qc(i,j,kv) -qrs(i,j,kv))
        zrho   (i,j) = (p0(i,j,kv) + pp(i,j,kv,nnow)) / zrtv
        zppdz  (i,j) = pp(i,j,kp1,nnow) - pp(i,j,km1,nnow)
        zpaidz (i,j) = zpai(i,j,kp1)    - zpai(i,j,km1)
        zdp0mz (i,j) = p0(i,j,kp1)      - p0(i,j,km1)
        zpp    (i,j) = p0(i,j,kv)  + pp(i,j,kv,nnow)
      ENDDO
    ENDDO

! perturbation of the east-west pressure gradient at (Arakawa-C) u-grid pts.

! --> need values within   (istartv-1,iendv+1) , (jstartv  ,jendv+1)
    DO    j = jstartv  , jendv+1
      DO  i = istartv-1, iendv
        z2rhomr =   c1 / (zrho(i+1,j) + zrho(i,j))
        zprefac =   zedacdl(j) * c2 * z2rhomr
        zfppdz  =  (p0(i+1,j,kv) - p0(i,j,kv)) / (zdp0mz(i+1,j) + zdp0mz(i,j))
        zpgxdrho (i,j) =  zprefac *( (pp(i+1,j,kv,nnow) - pp(i,j,kv,nnow))     &
                                    -(zppdz  (i+1,j) + zppdz  (i,j)) *zfppdz)  &
                        * z2rhomr *  (zroi(i+1,j,kv) + zroi(i,j,kv))
        zdpgxrho (i,j) =  zprefac *( (zpai(i+1,j,kv) - zpai(i,j,kv))           &
                                    -(zpaidz (i+1,j) + zpaidz (i,j)) *zfppdz)
      ENDDO
    ENDDO

! perturbation of the north-south pressure gradient at (Arakawa-C) v-grid pts.

! --> need values within   (istartu  ,iendu+1) , (jstartu-1,jendu+1)
    DO    j = jstartu-1, jendu
      DO  i = istartu  , iendu+1
        z2rhomr =   c1 / (zrho(i,j+1) + zrho(i,j))
        zprefac =   edadlat    * c2 * z2rhomr
        zfppdz  =  (p0(i,j+1,kv) - p0(i,j,kv)) / (zdp0mz(i,j+1) + zdp0mz(i,j))
        zpgydrho (i,j) =  zprefac *( (pp(i,j+1,kv,nnow) - pp(i,j,kv,nnow))     &
                                    -(zppdz  (i,j+1) + zppdz  (i,j)) *zfppdz)  &
                        * z2rhomr *  (zroi(i,j+1,kv) + zroi(i,j,kv))
        zdpgyrho (i,j) =  zprefac *( (zpai(i,j+1,kv) - zpai(i,j,kv))           &
                                    -(zpaidz (i,j+1) + zpaidz (i,j)) *zfppdz)
      ENDDO
    ENDDO

! full geostrophic wind increments

    DO   j = jstartu, jendu
      DO i = istartu, iendu
        zuwigef (i,j) =   (+ zpgydrho(i,j-1) + zpgydrho(i+1,j-1)               &
                           + zpgydrho(i,j  ) + zpgydrho(i+1,j  )               &
                           - zdpgyrho(i,j-1) - zdpgyrho(i+1,j-1)               &
                           - zdpgyrho(i,j  ) - zdpgyrho(i+1,j  ))              &
                        / (zfcmw(i,j-1)  + zfcmw(i,j))                         &
                        / (dp0(i+1,j,kv) + dp0(i,j,kv))
        zppu    (i,j) = c05 *(zpp(i+1,j) + zpp(i,j))
      ENDDO
    ENDDO
    DO   j = jstartv, jendv
      DO i = istartv, iendv
        zvwigef (i,j) =   (- zpgxdrho(i-1,j) - zpgxdrho(i-1,j+1)               &
                           - zpgxdrho(i  ,j) - zpgxdrho(i  ,j+1)               &
                           + zdpgxrho(i-1,j) + zdpgxrho(i-1,j+1)               &
                           + zdpgxrho(i  ,j) + zdpgxrho(i  ,j+1))              &
                        / (zfcmw(i-1,j)  + zfcmw(i,j))                         &
                        / (dp0(i,j+1,kv) + dp0(i,j,kv)) * zqcrlat(j)
! preparation for section 6
        zppv    (i,j) = c05 *(zpp(i,j+1) + zpp(i,j))
      ENDDO
    ENDDO


!-------------------------------------------------------------------------------
! Section 6: Reduction of geostrophic wind increments depending on parameters
!-------------------------------------------------------------------------------

!   zqwvgeo  = qgeo * zwvgeo
    zmaxgeo2 = (pmaxgeo * dtdeh) **2
    IF (.NOT. l2tls)  zmaxgeo2 = zmaxgeo2 * c2 * c2

! reduce geostrophic wind increments pressure dependently and close to ground

    DO   j = jstartu, jendu
      DO i = istartu, iendu
        zfqgeo        = (zpgref - zppu(i,j)) * zdpgrt
        zqgeo         = zfqgeo *qgeotop  +  (c1 - zfqgeo) *qgeo
        zuwigeo (i,j) = zuwigef(i,j) * zqgeo * zwvgeo
      ENDDO
    ENDDO
    DO   j = jstartv, jendv
      DO i = istartv, iendv
        zfqgeo        = (zpgref - zppv(i,j)) * zdpgrt
        zqgeo         = zfqgeo *qgeotop  +  (c1 - zfqgeo) *qgeo
        zvwigeo (i,j) = zvwigef(i,j) * zqgeo * zwvgeo
      ENDDO
    ENDDO

! print for control output

!   relative contribution of geopotential gradient vs. the density gradient term
!   ziec     = c0 ;   zjec     = c0 ;   ziec_tot = c0 ;   zjec_tot = c0
!   DO   j = jstart, jend
!     DO i = istart, iend
!       ziec      =  MAX( ziec , ABS( zdpgxrho(i,j) ) , ABS( zdpgyrho(i,j) ) )
!       zjec      =  MAX( zjec , ABS( zpgxdrho(i,j) ) , ABS( zpgydrho(i,j) ) )
!       ziec_tot  =  ziec_tot + ABS( zdpgxrho(i,j) )
!       zjec_tot  =  zjec_tot + ABS( zpgxdrho(i,j) )
!     ENDDO
!   ENDDO
!   ziec_tot =10000._wp * ziec_tot /(jend - jstart + 1) /(iend - istart + 1)
!   zjec_tot =10000._wp * zjec_tot /(jend - jstart + 1) /(iend - istart + 1)
!   ziec     =10000._wp * ziec
!   zjec     =10000._wp * zjec
!   IF ((MOD(kv,5) == 0) .AND. (MOD( ntstep,40 ) == 0))                        &
!     PRINT '("PE",I4,": GEOS",I4,I3,4F10.6)'   , my_cart_id                   &
!           , ntstep, kv, ziec, zjec, ziec_tot, zjec_tot

    IF (kv == kva)  lprgeo = (ntstep < 9) .AND. (.NOT. lfirst) .AND. (lwonl)   &
                                          .AND. (ntpscor > 0)
    IF ((lprgeo) .AND. (kv == ke)) THEN
      zpkemin = 109900.0_wp
      zdp1000 = 109900.0_wp
      zp1000  = 100000.0_wp * c05 *(vcoord%sigm_coord(ke+1) +              &
                                        vcoord%sigm_coord(ke))
      ipr (1) = 0
      jpr (1) = 0
      ipr (2) = 0
      jpr (2) = 0
      DO   j = jstart , jend
        DO i = istart , iend
          IF (ABS( psanai(i,j) ) > epsy) THEN
            IF (zpp(i,j) < zpkemin) THEN
              ipr(2) = i
              jpr(2) = j
              zpkemin = zpp(i,j)
            ENDIF
            IF (ABS( zpp(i,j)-zp1000 ) < zdp1000) THEN
              ipr(1) = i
              jpr(1) = j
              zdp1000 = ABS( zpp(i,j) - zp1000 )
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      IF (zpkemin > 109800.0_wp) lprgeo = .FALSE.
    ENDIF
    IF ((lprgeo) .AND. (kv == ke)) THEN
      ix = NINT( ptpstop / 100.0_wp )
      WRITE( nupr,'(" ")' )
      WRITE( nupr,'("Geostrophic wind correction for (surface) pr"             &
                  &,"essure nudging: top:",I4,"hPa, pts: (",I3,",",I3,")")' )  &
             ix, ipr(1), jpr(1)
      WRITE( nupr,'("============================================"             &
                  &,"================================= (",I3,",",I3,")")' )    &
                 ipr(2), jpr(2)
      WRITE( nupr,'("height  p weight u-geo  v-p v-w v-geo levl |"             &
                  &,"height  p weight u-geo  v-p v-w v-geo levl")' )
      WRITE( nupr,'("   m   hPa        cm/s  hPa      cm/s      |"             &
                  &,"   m   hPa        cm/s  hPa      cm/s     ")' )
    ENDIF
    IF (lprgeo) THEN
      DO ilpr = 1 , 2
        i = ipr(ilpr)
        j = jpr(ilpr)
        i = MAX( istartu, istartv, MIN( iendu, iendv, i ) )
        j = MAX( jstartu, jstartv, MIN( jendu, jendv, j ) )
        ixuv   = NINT( c05 * c05 * (  hhl(i  ,j,kv) + hhl(i  ,j,kv+1)          &
                                    + hhl(i+1,j,kv) + hhl(i+1,j,kv+1)) )
        zfqgeo = (zpgref - zppu(i,j)) * zdpgrt
        zqgeo  = zwvgeo * (zfqgeo *qgeotop  +  (c1 - zfqgeo) *qgeo)
        ix     = NINT( zppu(i,j) * .01_wp )
        ziec   = zuwigeo(i,j) * 100.0_wp
        zfqgeo = (zpgref - zppv(i,j)) * zdpgrt
        zfqgeo = zwvgeo * (zfqgeo *qgeotop  +  (c1 - zfqgeo) *qgeo)
        ixe    = NINT( zppv(i,j) * .01_wp )
        zjec   = zvwigeo(i,j) * 100.0_wp
        WRITE( ypr(ilpr),'(I5,2(I5,F5.2,F7.3),I3)' )                           &
               ixuv, ix, zqgeo, ziec, ixe, zfqgeo, zjec, kv
      ENDDO
      IF (lwonl)  WRITE( nupr,'(A42,'' |'',A42)' ) ypr(1), ypr(2)
    ENDIF

! exchange boundaries of geostrophic wind increments for their general limitat.

    IF (num_compute > 1) THEN
      IF (ltime) THEN
        CALL get_timings (i_nudging_corr, ntstep, dt, izerror)
        CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
        CALL get_timings (i_barrier_waiting_nud, ntstep, dt, izerror)
      ENDIF

      IF (lfirst) THEN

        kzdims(1:24)=(/1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)

        CALL exchg_boundaries ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart,   &
                        num_compute, ie, je, kzdims, jstartpar, jendpar, 1,       &
                        nboundlines, my_cart_neigh, lperi_x, lperi_y, l2dim,      &
                        20000+nexch_tag, .FALSE., ncomm_type , izerror, yzerrmsg, &
                        zuwigeo(:,:), zvwigeo(:,:), zuwigef(:,:), zvwigef(:,:) )
!       =====================

      ELSE

        kzdims(1:24)=(/1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)

        CALL exchg_boundaries ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart,   &
                        num_compute, ie, je, kzdims, jstartpar, jendpar, 1,       &
                        nboundlines, my_cart_neigh, lperi_x, lperi_y, l2dim,      &
                        20000+nexch_tag, .FALSE., ncomm_type , izerror, yzerrmsg, &
                        zuwigeo(:,:), zvwigeo(:,:) )
!       =====================

      ENDIF
      IF (ltime) CALL get_timings (i_communications_nud,ntstep, dt, izerror)
    ENDIF

! reduce geostrophic wind increments depending on geographic latitute
! (for ABS( latitude ) < 45 deg. ;   the SIN-function approximates
!  the latitude dependency obtained by the NMC method for GME / 3DVAR)

    !  note: for 'zabgeov', 'zuwigeo' is required from istartu-1, to jendu+1
    DO   j = jstartu  , jendu+1
      DO i = istartu-1, iendu
        zalatu = c05 * ABS( rlat(i,j) + rlat(i+1,j) )
        IF (zalatu < pi*0.25_wp)                                       &
          zuwigeo (i,j) = zuwigeo(i,j) * SIN( c2* zalatu )
      ENDDO
    ENDDO
    !  note: for 'zabgeou', 'zvwigeo' is required from jstartv-1, to iendv+1
    DO   j = jstartv-1, jendv
      DO i = istartv  , iendv+1
        zalatv = c05 * ABS( rlat(i,j) + rlat(i,j+1) )
        IF (zalatv < pi*0.25_wp)                                       &
          zvwigeo (i,j) = zvwigeo(i,j) * SIN( c2* zalatv )
      ENDDO
    ENDDO

! prepare reduction of geostrophic wind increments near lateral boundaries

    ziec_tot = c05 * ie_tot
    zjec_tot = c05 * je_tot
    IF (num_compute > 1) THEN
      ziec = ziec_tot - isubpos(my_cart_id,1) + nboundlines + 1
      zjec = zjec_tot - isubpos(my_cart_id,2) + nboundlines + 1
    ELSE
      ziec = ziec_tot
      zjec = zjec_tot
    ENDIF
    zlbciy = ziec_tot - plbcno - c05/ plbcipr
    zlbcjy = zjec_tot - plbcno - c05/ plbcipr

! get maximum geostrophic wind increments and
! the global coordinates of their location for printout for control

!   IF ((ntstep >= 11) .AND. (ntstep <= 11+INT( tconbox/dt-epsy ))) THEN
    IF (lfirst) THEN
      DO ix = 1 , 4
        zuvmx    (ix) = c0
        ijuvmx (1,ix) =  0
        ijuvmx (2,ix) =  0
      ENDDO
      DO   j = jstartu, jendu
        DO i = istartu, iendu
          zabsgeo =  zuwigef(i,j) **2                                          &
                   + 0.125_wp * (  zvwigef(i,j-1) + zvwigef(i+1,j-1)       &
                                     + zvwigef(i,j  ) + zvwigef(i+1,j  )) **2 
          IF (zabsgeo > zuvmx(1)) THEN
            zuvmx    (1) = zabsgeo
            ijuvmx (1,1) = i
            ijuvmx (2,1) = j
          ENDIF
          zwlbc =  MAX( c0 , c1 - MAX( c0, ABS(ziec-i-c05)-zlbciy ) *plbcipr ) &
                 * MAX( c0 , c1 - MAX( c0, ABS(zjec-j    )-zlbcjy ) *plbcipr )
          IF ((zabsgeo > zuvmx(3)) .AND. (zwlbc > 0.9_wp)) THEN
            zuvmx    (3) = zabsgeo
            ijuvmx (1,3) = i
            ijuvmx (2,3) = j
          ENDIF
        ENDDO
      ENDDO
      DO   j = jstartv, jendv
        DO i = istartv, iendv
          zabsgeo =  0.125_wp * (  zuwigef(i-1,j) + zuwigef(i-1,j+1)       &
                                     + zuwigef(i  ,j) + zuwigef(i  ,j+1)) **2  &
                   + zvwigef(i,j) **2
          IF (zabsgeo > zuvmx(2)) THEN
            zuvmx    (2) = zabsgeo
            ijuvmx (1,2) = i
            ijuvmx (2,2) = j
          ENDIF
          zwlbc =  MAX( c0 , c1 - MAX( c0, ABS(ziec-i    )-zlbciy ) * plbcipr) &
                 * MAX( c0 , c1 - MAX( c0, ABS(zjec-j-c05)-zlbcjy ) * plbcipr)
          IF ((zabsgeo > zuvmx(4)) .AND. (zwlbc > 0.9_wp)) THEN
            zuvmx    (4) = zabsgeo
            ijuvmx (1,4) = i
            ijuvmx (2,4) = j
          ENDIF
        ENDDO
      ENDDO
      DO ix = 1 , 4
        zuvmx (ix) = SQRT( zuvmx(ix) )
      ENDDO
      IF (num_compute > 1) THEN
        IF (ltime) THEN
          CALL get_timings (i_nudging_corr, ntstep, dt, izerror)
          CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
          CALL get_timings (i_barrier_waiting_nud, ntstep, dt, izerror)
        ENDIF

        CALL global_values ( zuvmx, 4, 'MAX', imp_reals, icomm_cart            &
                           , onl_cart_id, ijuvmx, yzerrmsg, izerror)
!       ==================

        IF (ltime) CALL get_timings (i_communications_nud,ntstep,dt,izerror)
      ENDIF
    ENDIF

! reduce geostrophic wind increments by imposing a general upper limit,
!                                    and near the lateral boundaries

    DO   j = jstartu, jendu
      DO i = istartu, iendu
        zabgeou (i,j) =  zuwigeo(i,j) **2                                      &
                       + 0.0625_wp *(  zvwigeo(i,j-1) + zvwigeo(i+1,j-1)   &
                                         + zvwigeo(i,j  ) + zvwigeo(i+1,j  ))**2 
      ENDDO
    ENDDO
    DO   j = jstartv, jendv
      DO i = istartv, iendv
        zabgeov (i,j) =  zvwigeo(i,j) **2                                      &
                       + 0.0625_wp *(  zuwigeo(i-1,j) + zuwigeo(i-1,j+1)   &
                                         + zuwigeo(i  ,j) + zuwigeo(i  ,j+1))**2
      ENDDO
    ENDDO
    zabsgeo = c0
    DO   j = jstartu, jendu
      DO i = istartu, iendu
        IF (zabgeou(i,j) > zmaxgeo2)                                           &
          zuwigeo (i,j) = zuwigeo(i,j) * SQRT( zmaxgeo2 / zabgeou(i,j) )
! reduce geostrophic wind increments near lateral boundaries
        zwlbc =  MAX( c0 , c1 - MAX( c0, ABS(ziec-i-c05)-zlbciy ) * plbcipr )  &
               * MAX( c0 , c1 - MAX( c0, ABS(zjec-j    )-zlbcjy ) * plbcipr )
        zuwigeo (i,j) = zuwigeo(i,j) * zwlbc
        zabsgeo  =  MAX( zabsgeo , zabgeou(i,j) *zwlbc *zwlbc )
      ENDDO
    ENDDO
    DO   j = jstartv, jendv
      DO i = istartv, iendv
        IF (zabgeov(i,j) > zmaxgeo2)                                           &
          zvwigeo (i,j) = zvwigeo(i,j) * SQRT( zmaxgeo2 / zabgeov(i,j) )
! reduce geostrophic wind increments near lateral boundaries
        zwlbc =  MAX( c0 , c1 - MAX( c0, ABS(ziec-i    )-zlbciy ) * plbcipr )  &
               * MAX( c0 , c1 - MAX( c0, ABS(zjec-j-c05)-zlbcjy ) * plbcipr )
        zvwigeo (i,j) = zvwigeo(i,j) * zwlbc
        zabsgeo  =  MAX( zabsgeo , zabgeov(i,j) *zwlbc *zwlbc )
!       IF ((zabgeov(i,j) > zmaxgeo2) .AND. (zwlbc > 0.95_wp)              &
!                                     .AND. (MOD(ntstep,45) <= 3))             &
!         PRINT '("Large GEO-v:",I4,2F7.3,F6.3,3I4)'                           &
!                , my_cart_id, zabgeov(i,j), zmaxgeo2, zwlbc, i, j, kv
      ENDDO
    ENDDO
!   IF (zabsgeo > zmaxgeo2*0.5_wp)                                         &
!     PRINT '("Large GEO:",2I4,2F7.3,I4)',my_cart_id,ntstep,zabsgeo, zmaxgeo2,kv

!-------------------------------------------------------------------------------
! Section 7: Printout for control of geostrophic wind increments
!-------------------------------------------------------------------------------

! printout for control

!   IF ((ntstep >= 11) .AND. (ntstep <= 11+INT( tconbox/dt-epsy ))) THEN
    IF (lfirst) THEN
      ixuv  = 1
      IF (zuvmx(2) > zuvmx(1)) ixuv = 2
      IF (lwonl) THEN
        IF (ixuv == 1) yuvgeo = ', u-gr.pt: '
        IF (ixuv == 2) yuvgeo = ', v-gr.pt: '
        WRITE( nupr,'(4I4,'': max. unreduced geost. incr. on total domain''    &
                    &,A11,F7.2)' )                                             &
               ntstep, kv, ijuvmx(1,ixuv), ijuvmx(2,ixuv), yuvgeo, zuvmx(ixuv)
        ixuv  = 3
        IF (zuvmx(4) > zuvmx(3)) ixuv = 4
        IF (ixuv == 3) yuvgeo = ', u-gr.pt: '
        IF (ixuv == 4) yuvgeo = ', v-gr.pt: '
        WRITE( nupr,'(4I4,'': max. unreduced geost. incr. on inner domain''    &
                    &,A11,F7.2)' )                                             &
               ntstep, kv, ijuvmx(1,ixuv), ijuvmx(2,ixuv), yuvgeo, zuvmx(ixuv)
        WRITE( nupr,'(4I4,F5.2,'': reduced geost. incr. at i/jonl:'',2F8.3)' ) &
               ntstep, kv, i, j, zwvgeo, zuwigeo(ionl,jonl), zvwigeo(ionl,jonl)
      ENDIF
    ENDIF
! prepare further printout (in section 8)
    IF ((lwonl) .AND. (ntstep < 121)) THEN
      zugeai(kv) = zuwigeo(ionl,jonl)
      zvgeai(kv) = zvwigeo(ionl,jonl)
    ENDIF
  ELSEIF ((lwonl) .AND. (ntstep < 121)) THEN
    zugeai(kv) = c0
    zvgeai(kv) = c0

  ENDIF

 
!-------------------------------------------------------------------------------
! Section 8: Analysis increments of horizontal wind,   updating the
!            prognostic horizontal wind field by adding these increments
!-------------------------------------------------------------------------------

  kwi = kv
  IF (.NOT. lconai) kwi = 1

  IF ((lgetai) .AND. (zwvgeo > epsy)) THEN
    DO   j = jstartu, jendu
      DO i = istartu, iendu
        uanai (i,j,kwi)  =  zumydt2(i,j) / (c1 + zumydt2(i,j)) * zuwipos(i,j)  &
                          +      c1      / (c1 + zumydt2(i,j)) * zuwigeo(i,j)
      ENDDO
    ENDDO
    DO   j = jstartv, jendv
      DO i = istartv, iendv
        vanai (i,j,kwi)  =  zvmydt2(i,j) / (c1 + zvmydt2(i,j)) * zvwipos(i,j)  &
                          +      c1      / (c1 + zvmydt2(i,j)) * zvwigeo(i,j)
      ENDDO
    ENDDO
  ELSEIF (lgetai) THEN
    DO   j = jstartu, jendu
      DO i = istartu, iendu
        uanai (i,j,kwi)  =  zumydt2(i,j) / (c1 + zumydt2(i,j)) * zuwipos(i,j)
      ENDDO
    ENDDO
    DO   j = jstartv, jendv
      DO i = istartv, iendv
        vanai (i,j,kwi)  =  zvmydt2(i,j) / (c1 + zvmydt2(i,j)) * zvwipos(i,j)
      ENDDO
    ENDDO
  ENDIF

  DO   j = jstartu, jendu
    DO i = istartu, iendu
      u (i,j,kv,nnew) =  u(i,j,kv,nnew) + uanai(i,j,kwi)
    ENDDO
  ENDDO
  DO   j = jstartv, jendv
    DO i = istartv, iendv
      v (i,j,kv,nnew) =  v(i,j,kv,nnew) + vanai(i,j,kwi)
    ENDDO
  ENDDO

  IF ((lwonl2) .AND. (kv == ke-5) .AND. (ntstep*dtdeh <= c2))                  &
    WRITE( nupr,'(''nudge_horiz_wind: uv-ana.incr. '',3I4, 2F9.4)' )           &
           ionl2, jonl2, kv, uanai(ionl2,jonl2,kwi), vanai(ionl2,jonl2,kwi)

! control output of max. (analysis or net spreaded observation) increments

  IF (lgetai) THEN
    IF (kv == ke) THEN
      IF (ALLOCATED( zxmax )) THEN
        PRINT '("CAUTION in src_nudging: zxmax is already allocated "          &
              &,"at time / PE",I6,I5)', ntstep, my_cart_id
        DEALLOCATE ( zxmax , ijxmx , STAT=nzerr )
      ENDIF
      ALLOCATE ( zxmax (  12*ke) , stat=nzerr )
      ALLOCATE ( ijxmx (2,12*ke) , stat=nzerr )
    ENDIF
    ixe = kv + 11*ke
    DO  ix = kv , ixe , ke
      zxmax   (ix) = c0
      ijxmx (1,ix) = 0
      ijxmx (2,ix) = 0
    ENDDO
    DO   j = jstartu, jendu
!for the NEC
!CDIR NODEP
      DO i = istartu, iendu
        ixe = kv
        IF (ABS(uanai(i,j,kwi)) > ABS(zxmax(ixe))) THEN
          zxmax   (ixe) = uanai(i,j,kwi)
          ijxmx (1,ixe) = i
          ijxmx (2,ixe) = j
        ENDIF
        ixe = kv + 2*ke
        IF (ABS(zuwipos(i,j)) > ABS(zxmax(ixe))) THEN
          zxmax   (ixe) = zuwipos(i,j)
          ijxmx (1,ixe) = i
          ijxmx (2,ixe) = j
        ENDIF
        IF (zwvgeo > epsy) THEN
          ixe = kv + 4*ke
          IF (ABS(zuwigeo(i,j)) > ABS(zxmax(ixe))) THEN
            zxmax   (ixe) = zuwigeo(i,j)
            ijxmx (1,ixe) = i
            ijxmx (2,ixe) = j
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    DO   j = jstartv, jendv 
!for the NEC
!CDIR NODEP
      DO i = istartv, iendv
        ixe = kv + ke
        IF (ABS(vanai(i,j,kwi)) > ABS(zxmax(ixe))) THEN
          zxmax   (ixe) = vanai(i,j,kwi)
          ijxmx (1,ixe) = i
          ijxmx (2,ixe) = j
        ENDIF
        ixe = kv + 3*ke
        IF (ABS(zvwipos(i,j)) > ABS(zxmax(ixe))) THEN
          zxmax   (ixe) = zvwipos(i,j)
          ijxmx (1,ixe) = i
          ijxmx (2,ixe) = j
        ENDIF
        IF (zwvgeo > epsy) THEN
          ixe = kv + 5*ke
          IF (ABS(zvwigeo(i,j)) > ABS(zxmax(ixe))) THEN
            zxmax   (ixe) = zvwigeo(i,j)
            ijxmx (1,ixe) = i
            ijxmx (2,ixe) = j
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    IF ( (num_compute > 1) .AND. (kv == 1)) THEN
      DO  kk = 1 , ke
        DO  ix = kk , kk+5*ke , ke
          ixe = ix + 6*ke
          zxmax   (ixe) = - zxmax  (ix)
          ijxmx (1,ixe) =   ijxmx(1,ix)
          ijxmx (2,ixe) =   ijxmx(2,ix)
        ENDDO
      ENDDO
      ixe = 12 * ke
      IF (ltime) THEN
        CALL get_timings (i_nudging_corr, ntstep, dt, izerror)
        CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
        CALL get_timings (i_barrier_waiting_nud, ntstep, dt, izerror)
      ENDIF

      CALL global_values ( zxmax, ixe, 'MAX', imp_reals, icomm_cart            &
                         , onl_cart_id, ijxmx, yzerrmsg, izerror)
!     ==================

      IF (ltime) CALL get_timings (i_communications_nud,ntstep, dt, izerror)
      IF (lwonl) THEN
        DO  kk = 1 , ke
          DO  ix = kk , kk+5*ke , ke
            ixe = ix + 6*ke
            IF (ABS(zxmax(ixe)) > ABS(zxmax(ix))) THEN
              zxmax   (ix) = - zxmax  (ixe)
              ijxmx (1,ix) =   ijxmx(1,ixe)
              ijxmx (2,ix) =   ijxmx(2,ixe)
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDIF

    IF ((lwonl) .AND. (kv == 1)) THEN
      WRITE( nupr,'(''max. (ABS of) increments of horizontal wind field ''     &
                  &,''components in [m/s]'')')
      WRITE( nupr,'(''(ana incr: analysis increments'')' )
      WRITE( nupr,'('' net incr: net observation increments'')' )
      WRITE( nupr,'('' geo corr: geostrophic wind correction)'')' )
      WRITE( nupr,'(''le- u-ana  global  v-ana  global  u-net  global''        &
                  &,''  v-net  global  u-geo  global  v-geo  global'')' )
      WRITE( nupr,'(''vel  incr  coord.   incr  coord.   incr  coord.''        &
                  &,''   incr  coord.   corr  coord.   corr  coord.'')' )
      DO kk = 1 , ke
        ixe = kk + 5*ke
        WRITE( nupr,'(I2, 2(F7.2,2I4), 2(F7.1,2I4), 2(F7.2,2I4))' )   kk       &
                   , (zxmax(ix), ijxmx(1,ix), ijxmx(2,ix) , ix=kk,ixe,ke)
      ENDDO
    ENDIF
    IF (kv == 1) THEN
      DEALLOCATE ( zxmax , stat=nzerr )
      DEALLOCATE ( ijxmx , stat=nzerr )
    ENDIF
  ENDIF


ENDDO loop_over_model_levels
!~~~~~~~~~~~~~~~~~~~~~~~~~~~

!-------------------------------------------------------------------------------
! Section 9: Time integration of analysis increments
!            (done at k=1 as it requires exchg_boundaries)
!-------------------------------------------------------------------------------

  IF ((k == 1) .AND. (num_compute > 1)) THEN
    IF (ltime) THEN
      CALL get_timings (i_nudging_corr, ntstep, dt, izerror)
      CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
      CALL get_timings (i_barrier_waiting_nud, ntstep, dt, izerror)
    ENDIF

! Exchange of 1 row of boundary data (of local domains)

    kzdims(1:24)=(/ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)

    CALL exchg_boundaries ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart,       &
                           num_compute, ie, je, kzdims, jstartpar, jendpar, 1,    &
                           nboundlines, my_cart_neigh, lperi_x, lperi_y, l2dim,   &
                           20000+nexch_tag,.FALSE., ncomm_type, izerror, yzerrmsg,&
                           u(:,:,:,nnew), v(:,:,:,nnew),                          &
                           uanai(:,:,:) , vanai(:,:,:) )
!   =====================

    IF (ltime) CALL get_timings (i_communications_nud, ntstep, dt, izerror)
  ENDIF

  IF (k == 1) THEN
    DO  kv = 1 , ke
      DO   j = jstart, jend
        DO i = istart, iend
          zumn   =         c05* (u(i,j,kv,nnew) + u(i-1,j,kv,nnew))
          zumo   =  zumn - c05* (uanai(i,j,kv ) + uanai(i-1,j,kv ))
          zvmn   =         c05* (v(i,j,kv,nnew) + v(i,j-1,kv,nnew))
          zvmo   =  zvmn - c05* (vanai(i,j,kv ) + vanai(i,j-1,kv ))

! for differences of wind speed and direction, transformation from
! rotated to geographical coordinates is not necessary
!     CALL uv2df_vec (zumn, zvmn, zddn, zffn, ie, je)
!     CALL uv2df_vec (zumo, zvmo, zddo, zffo, ie, je)
          ! call uv2df inside loop (instead of call uv2df_vec outside loop)
          ! allows inlining (faster on NEC)
          CALL uv2df (zumn, zvmn, zddn, zffn)
          CALL uv2df (zumo, zvmo, zddo, zffo)
!         ==========
          ff_anai (i,j,kv) = ff_anai(i,j,kv) + zffn - zffo
          IF (MIN( zffo,zffn ) > ddfflim) THEN
!           dd_anai (i,j,kv) = dd_anai(i,j,kv) + zddn - zddo
            zdanai  =  zddn - zddo
            IF (zdanai >  180._wp)  zdanai  =  zdanai - 360._wp
            IF (zdanai < -180._wp)  zdanai  =  zdanai + 360._wp
            dd_anai (i,j,kv) = dd_anai(i,j,kv) + zdanai
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDIF


!-------------------------------------------------------------------------------
! Section 10: General printout for control
!-------------------------------------------------------------------------------

! vertical profiles of horizontally averaged analysis increments
! --------------------------------------------------------------

! nulast must be the same as at call of obs_org_print_stat (due to nustat)
  nulast  =  MIN( INT( MAX( nverend , nudgend ) ,iintegers) , nstop )

  IF ((k == 1) .AND. (     (MOD( NINT((ntstep+c05)*dt),3600 ) < NINT(dt))      &
                      .OR. (ntstep == nulast))) THEN
    avai_tqv  =  c0
    avai_tqc  =  c0
    avai_dse  =  c0
    zmass     =  c0
    navai     = (jend - jstart + 1) * (iend - istart + 1)
    DO kv = ke , 1 , -1
      avai_ff (kv) = c0
      avai_dd (kv) = c0
      avai_t  (kv) = c0
      avai_p  (kv) = c0
      avai_qv (kv) = c0
      avai_qc (kv) = c0
      avai_fi (kv) = c0
      DO   j = jstart, jend
        DO i = istart, iend
          avai_ff (kv) = avai_ff(kv) + ff_anai(i,j,kv)
          avai_dd (kv) = avai_dd(kv) + dd_anai(i,j,kv)
          avai_t  (kv) = avai_t (kv) +  t_anai(i,j,kv)
          avai_p  (kv) = avai_p (kv) +  p_anai(i,j,kv)
          avai_qv (kv) = avai_qv(kv) + qv_anai(i,j,kv)
          avai_qc (kv) = avai_qc(kv) + qc_anai(i,j,kv)
          avai_fi (kv) = avai_fi(kv) +  p_anai(i,j,kv) / rho(i,j,kv)
          avai_tqv     = avai_tqv    + qv_anai(i,j,kv) * rho(i,j,kv)           &
                                        * (hhl(i,j,kv) - hhl(i,j,kv+1)) 
          avai_tqc     = avai_tqc    + qc_anai(i,j,kv) * rho(i,j,kv)           &
                                        * (hhl(i,j,kv) - hhl(i,j,kv+1)) 
          avai_dse     = avai_dse    +  t_anai(i,j,kv) * dp0(i,j,kv) * cp_d 
          IF (kv == 1) zmass = zmass +    ( p0(i,j,ke) + dp0(i,j,ke) *c05 )    &
                                     -    ( p0(i,j,1 ) - dp0(i,j,1 ) *c05 )
        ENDDO
      ENDDO
    ENDDO
    avai_dse  =  avai_dse / zmass

    ixe = 7*ke + 4
    IF (ALLOCATED( zxmax )) THEN
      PRINT '("CAUTION in src_nudging: zxmax is already allocated "          &
            &,"at time / PE",I6,I5)', ntstep, my_cart_id
      DEALLOCATE ( zxmax , STAT=nzerr )
    ENDIF
    ALLOCATE ( zxmax (ixe) , stat=nzerr )
    DO kv = ke , 1 , -1
      zxmax(kv     ) = avai_ff(kv)
      zxmax(kv+1*ke) = avai_dd(kv)
      zxmax(kv+2*ke) = avai_t (kv)
      zxmax(kv+3*ke) = avai_p (kv)
      zxmax(kv+4*ke) = avai_qv(kv)
      zxmax(kv+5*ke) = avai_qc(kv)
      zxmax(kv+6*ke) = avai_fi(kv)
    ENDDO
    zxmax(1+7*ke) = avai_tqv
    zxmax(2+7*ke) = avai_tqc
    zxmax(3+7*ke) = avai_dse
    zxmax(4+7*ke) = REAL( navai , wp )
    IF (num_compute > 1) THEN
      IF (ltime) THEN
        CALL get_timings (i_nudging_corr, ntstep, dt, izerror)
        CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
        CALL get_timings (i_barrier_waiting_nud, ntstep, dt, izerror)
      ENDIF

      CALL global_values ( zxmax, ixe, 'SUM', imp_reals, icomm_cart            &
                         , 0, yzerrmsg, izerror)
!     ==================

      IF (ltime) CALL get_timings (i_communications_nud,ntstep, dt, izerror)
    ENDIF

    IF (lcroot) THEN
      navai  =  NINT( zxmax(4+7*ke) )
      DO kv = ke , 1 , -1
        avai_ff(kv) = zxmax(kv     ) / navai
        avai_dd(kv) = zxmax(kv+1*ke) / navai * 57.296_wp
        avai_t (kv) = zxmax(kv+2*ke) / navai
        avai_p (kv) = zxmax(kv+3*ke) / navai *   0.01_wp
        avai_qv(kv) = zxmax(kv+4*ke) / navai * 1000.0_wp
        avai_qc(kv) = zxmax(kv+5*ke) / navai * 1000000.0_wp
        avai_fi(kv) = zxmax(kv+6*ke) / navai
      ENDDO
      avai_tqv = zxmax(1+7*ke) / navai * 1000.0_wp
      avai_tqc = zxmax(2+7*ke) / navai * 1000000.0_wp
      avai_dse = zxmax(3+7*ke) / navai

      OPEN (nustat,FILE=yustats,FORM='FORMATTED',STATUS='UNKNOWN'              &
                               ,POSITION='APPEND',IOSTAT=nstat)
      IF (nstat /= 0) PRINT'("OPENING OF FILE yustats FAILS, nudge_horiz_wind")'
      WRITE( ychh,'(F5.2)' ) ntstep *dtdeh
      WRITE( nustat, '(" ")' )
      WRITE( nustat,'("1",5X,"Diagnostics on time-integrated analysis "        &
                    &,"increments")' )
      WRITE( nustat,'(    6X,A14," +",A5," h:")' )  ydate_ini, ychh
      WRITE( nustat,'(    6X,"horizontal mean of ",A5,"-hourly sums of analy"  &
                    &,"sis increments")' )  ychh
      WRITE( nustat,'(24X,"TQV [g/m2]    TQC [mg/m2]       Dry Static Energy " &
                    &,"[J/kg]")' )
      WRITE( nustat,'(22X,F9.3,6X,F9.3,4X,7X,F9.3," (volumn mean)")' )         &
                      avai_tqv, avai_tqc, avai_dse
      WRITE( nustat,'(6X,"level p0(k)    |v|    wind dir     T        p   "    &
                    &,"     FI      qv       qc")' )
      WRITE( nustat,'(6X,"      hPa      m/s       deg       K       hPa  "    &
                    &,"   m2/s2    g/kg    mg/kg")' )
      DO kk = 1 , ke
        WRITE( nustat,'(5X,I3,F9.2,F9.3,F9.2,F9.3,F9.2,F9.1,2F9.3)' )          &
         kk, c05*(vcoord%sigm_coord(kk+1)+vcoord%sigm_coord(kk))*1000._wp  &
             , avai_ff(kk), avai_dd(kk), avai_t(kk), avai_p(kk), avai_fi(kk)   &
             , avai_qv(kk), avai_qc(kk)  
      ENDDO
      CLOSE( nustat )
    ENDIF
    DEALLOCATE ( zxmax , stat=nzerr )
  ENDIF

! other output
! ------------

  IF ((lgetai) .AND. (lwonl) .AND. (ntstep < 121)) THEN
    IF (k == ke) THEN
      IF (ALLOCATED( ztpsai )) THEN
        PRINT '("CAUTION in src_nudging: ztpsai is already allocated "         &
              &,"at time / PE",I6,I5)', ntstep, my_cart_id
        DEALLOCATE ( ztpsai, ztnuai, ztecai, zqnuai, zqecai, zpecai, STAT=nzerr)
      ENDIF
      ALLOCATE ( ztpsai(ke) , stat=nzerr )
      ALLOCATE ( ztnuai(ke) , stat=nzerr )
      ALLOCATE ( ztecai(ke) , stat=nzerr )
      ALLOCATE ( zqnuai(ke) , stat=nzerr )
      ALLOCATE ( zqecai(ke) , stat=nzerr )
      ALLOCATE ( zpecai(ke) , stat=nzerr )
    ENDIF
    ztpsai(k) =         c1 /(c1 + zdt *omy(ionl,jonl,2,1)) * ztwips(ionl,jonl,k)
    ztnuai(k) = ztpsai(k) +       zdt *omy(ionl,jonl,2,1)                      &
                           /(c1 + zdt *omy(ionl,jonl,2,1)) * zwi(ionl,jonl,3,1)
    zqnuai(k) = rvd_m_o   *       zdt *omy(ionl,jonl,3,1)                      &
                           /(c1 + zdt *omy(ionl,jonl,3,1)) * zwi(ionl,jonl,4,1)
    ztecai(k) = ztaionl
    zqecai(k) = zqaionl
    zpecai(k) = zpaionl
  ENDIF

  IF ((lgetai) .AND. (lwonl) .AND. (k == 1) .AND. (ntstep < 121)) THEN
    IF (     (ntstep < 11)                                                     &
        .OR. (MOD(ntstep,20) < INT( rmod(tconbox,dt,epsy)-epsy ))) THEN
      zpnoai(ke) = zpecai(ke)
      zppsai(ke) = zpecai(ke)
      zpnuai(ke) = zpecai(ke)
      zpnew (ke) = pp(ionl,jonl,ke,nnew) + p0(ionl,jonl,ke)
      DO kk = ke-1 , 1 , -1
        zfbuoyt   = r_d * rho0(ionl,jonl,kk  ) * t(ionl,jonl,kk  ,nnow)
        zfbuoyb   = r_d * rho0(ionl,jonl,kk+1) * t(ionl,jonl,kk+1,nnow)
        zfbuoyt   = c1 / (c2 + dp0(ionl,jonl,kk+1) /zfbuoyt)
        zfbuoyb   =      (c2 - dp0(ionl,jonl,kk  ) /zfbuoyb)
        zdp01dt   = dp0(ionl,jonl,kk+1) / t(ionl,jonl,kk  ,nnow)
        zdp0dt1   = dp0(ionl,jonl,kk  ) / t(ionl,jonl,kk+1,nnow)
        zpnoai(kk) = zfbuoyt *(  zpnoai(kk+1) * zfbuoyb )
        zppsai(kk) = zfbuoyt *(  zppsai(kk+1) * zfbuoyb                        &
                               + ztpsai(kk+1) * zdp0dt1 + ztpsai(kk) * zdp01dt )
        zpnuai(kk) = zfbuoyt *(  zpnuai(kk+1) * zfbuoyb                        &
                               + ztnuai(kk+1) * zdp0dt1 + ztnuai(kk) * zdp01dt &
                               + zqnuai(kk  ) * dp0(ionl,jonl,kk+1)            &
                               + zqnuai(kk+1) * dp0(ionl,jonl,kk) )
        zpnew (kk) = pp(ionl,jonl,kk,nnew) + p0(ionl,jonl,kk)
        zqnuai(kk+1) = zqnuai(kk+1) * 100000.0_wp
        zqecai(kk+1) = zqecai(kk+1) * 100000.0_wp
        IF (.NOT. lgetgeo) THEN
          zugeai (kk+1) = c0
          zvgeai (kk+1) = c0
        ENDIF
      ENDDO
      zqnuai(1) = zqnuai(1) * 100000.0_wp
      zqecai(1) = zqecai(1) * 100000.0_wp
      kdpr = (ke-3) / 4
      kepr = 1 + 4*kdpr
      WRITE( nupr,'(''nudging, pressure                   '',5F7.0,F8.0)' )    &
           ( zpnew (kk), kk=1,kepr,kdpr ) , zpnew (ke)
      WRITE( nupr,'(''nudging, p-incr, no T-change        '',5F7.3,F8.3)' )    &
           ( zpnoai(kk), kk=1,kepr,kdpr ) , zpnoai(ke)
      WRITE( nupr,'(''nudging, p-incr, T-corr, no T-nudge '',5F7.3,F8.3)' )    &
           ( zppsai(kk), kk=1,kepr,kdpr ) , zppsai(ke)
      WRITE( nupr,'(''nudging, p-incr, T-nudge, no evapo. '',5F7.3,F8.3)' )    &
           ( zpnuai(kk), kk=1,kepr,kdpr ) , zpnuai(ke)
      WRITE( nupr,'(''nudging, p-incr, T-nudge complete   '',5F7.3,F8.3)' )    &
           ( zpecai(kk), kk=1,kepr,kdpr ) , zpecai(ke)
      WRITE( nupr,'(''nudging, T-incr, T-corr, no T-nudge '',5F7.3,F8.3)' )    &
           ( ztpsai(kk), kk=1,kepr,kdpr ) , ztpsai(ke)
      WRITE( nupr,'(''nudging, T-incr, T-nudge, no evapo. '',5F7.3,F8.3)' )    &
           ( ztnuai(kk), kk=1,kepr,kdpr ) , ztnuai(ke)
      WRITE( nupr,'(''nudging, T-incr, T-nudge complete   '',5F7.3,F8.3)' )    &
           ( ztecai(kk), kk=1,kepr,kdpr ) , ztecai(ke)
      WRITE( nupr,'(''nudging, q-incr, RH-nudge, no evapo.'',5F7.3,F8.3)' )    &
           ( zqnuai(kk), kk=1,kepr,kdpr ) , zqnuai(ke)
      WRITE( nupr,'(''nudging, q-incr, T+RH-nudge complete'',5F7.3,F8.3)' )    &
           ( zqecai(kk), kk=1,kepr,kdpr ) , zqecai(ke)
      IF (lgetgeo) THEN
        kwi = 1 + 3*kdpr
        IF (.NOT. lconai) kwi = 1
        WRITE( nupr,'(''nudging, u-incr, geostrophic        '',5F7.3,F8.3)' )  &
             ( zugeai(kk), kk=1,kepr,kdpr ) , zugeai(ke)
        WRITE( nupr,'(''nudging, v-incr, geostrophic        '',5F7.3,F8.3)' )  &
             ( zvgeai(kk), kk=1,kepr,kdpr ) , zvgeai(ke)
!            ( zvgeai(kk), kk=1,kepr,kdpr ) , vanai(ionl,jonl,kwi)
      ENDIF
    ENDIF
  ENDIF


!-------------------------------------------------------------------------------
! Section 11: Clean-up
!-------------------------------------------------------------------------------

  IF ((.NOT. ((num_compute > 1) .OR. (lgetgeo))) .OR. (k == 1)) THEN
    DEALLOCATE ( zuwi , stat=nzerr )
    DEALLOCATE ( zvwi , stat=nzerr )
    DEALLOCATE ( zumy , stat=nzerr )
  ENDIF

  IF ((lgetai) .AND. (lwonl) .AND. (ntstep < 121) .AND. (k == 1)) THEN
    DEALLOCATE ( ztpsai , stat=nzerr )
    DEALLOCATE ( ztnuai , stat=nzerr )
    DEALLOCATE ( ztecai , stat=nzerr )
    DEALLOCATE ( zqnuai , stat=nzerr )
    DEALLOCATE ( zqecai , stat=nzerr )
    DEALLOCATE ( zpecai , stat=nzerr )
  ENDIF

  IF (ldump_ascii) THEN
    ! flush YUPRINT file
    ! Note: 'FLUSH' is standard in Fortran2003, but not in Fortran95.
    !       'Call flush' is available in the NEC compiler library,
    !       but not on every platform;
    !       it would be more efficient than close and open.
    IF ((lwonl) .AND. ((k == 1) .OR. ((k == ke) .AND. (lgetai))                &
                                .OR. ((k >= ke-5) .AND. (lfirst)))) THEN
      CLOSE (nupr)
      OPEN  (nupr,FILE=yuprint,FORM='FORMATTED',STATUS='OLD'                   &
                              ,POSITION='APPEND',IOSTAT=nstat)
      IF (nstat /= 0) PRINT '("OPENING OF FILE yuprint FAILED, nudge_horiz_wind")'
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure nudge_horiz_wind
!-------------------------------------------------------------------------------

END SUBROUTINE nudge_horiz_wind


!-------------------------------------------------------------------------------
!+ Module procedure for updating the mass fields by nudging
!-------------------------------------------------------------------------------

SUBROUTINE nudge_humid_mass ( lgetai , lconai , k )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure updates the temperature, humidity and pressure fields
!   by execution of nudging equations, and by a hydrostatic adjustment for the
!   upper-level pressure.
!   The following analysis increments are stored:
!   - temperature and humidity excluding effects of condensation and evaporation
!     if identical analysis increments are used during several timesteps
!   - temperature, humidity and pressure, including effects of condensation and
!     evaporation prior to the updating of the pressure
!   - pressure and density (from nudging the pressure at the lowest model and
!     applying the temperature correction, excluding effects of condensation and
!     evaporation), if geostrophic wind increments are to be computed
!
! Method:
!   1.  The pressure at the lowest model level is updated by adding previously 
!       computed pressure analysis increments
!       (--> "Nudging of 'surface' pressure").
!   2.  For temperature, 'net (spreaded) observation increments' (squared
!       weighted mean of observation increments) and 'net nudging weights'
!       (weighted mean of nudging weights) at the target grid points are
!       derived from the sum of weights, the sum of squared weights, and the
!       sum of squared weighted observation increments as available from the
!       spreading procedures. This derivation of weighting of multiple data
!       is based on Benjamin & Seaman (MWR, 1985) and accounts for data density
!       at the model grid points.
!       Analysis increments for temperature are then computed from these
!       'net increments' and 'net weights', and are added to the model fields.
!       (--> "Nudging of temperature")
!     \ To prepare '3.', alternative temperature analysis increments are
!     \ derived analogously, but having set the observation increments to zero
!     \ (in 'src_mult_spread', 'src_sing_spread') for those temperature obser-
!     \ vations, for which specific instead of relative humidity shoud be
!     \ preserved.
!   3.  The specific humidity is adjusted so that the relative humidity is as
!       before nudging the temperature (apart from the effect of those tempera-
!       ture observations, for which specific humidity shoud be preserved).
!   4.  'Net nudging weights' and 'net (spreaded) observation increments' of
!       relative humidity are computed, and are (approximately) converted into
!       'net observation increments' of specific humidity  -  errors due to the
!       applied approximation are usually far less than 0.1% relative humidity.
!   5.  Analysis increments for specific water vapour content are computed from
!       the 'net increments' and weights, and are added to the model fields.
!       (--> "Nudging of humidity")
!   6.  The saturation adjustment technique is applied for condensation and 
!       evaporation. 
!   7.  Analysis increments of pressure at upper levels are derived hydrostatic-
!       ally from a) pressure analysis increments at the lowest model level and
!                 b) analysis increments of temperature and humidity, 
!       and added to the model fields. Optionally, this is followed by another
!       application of the saturation adjustment technique.
!
! Written by        :  Christoph Schraff, DWD  (original version: 10.06.97)
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

  LOGICAL                 , INTENT (IN)         ::       &
    lgetai           ,& ! .TRUE if new analysis increment are to be computed 
                        !          and used at the current timestep
    lconai              ! TRUE if analysis increments are const during 'tconbox'

  INTEGER (KIND=iintegers), INTENT (IN)         ::       &
    k                   ! index of current vertical model level

! Local parameters:
! ----------------
  LOGICAL , PARAMETER      ::  &
    lpsatad = .FALSE.   !

! REAL    (KIND=wp   )     , PARAMETER ::  &
!   fswtp = 0.0_wp  ! saturation factor, which deminishes the weight given
                        ! to the temperature correction for surface pressure
                        ! nudging, depending on the influence of temperature
                        ! observations (--> NAMELIST !!?)

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    izerror          ,& ! error index
    i      , j       ,& ! loop indices in horizontal direction
    kk               ,& ! vertical index
    kwi              ,& ! vertical index for analysis increments
    kp1              ,& ! k + 1
    ityp             ,& ! index of set of obs system in weighted increment array
    ix     , ixe     ,& ! indices used for control output
    kitera           ,& ! number of iterations in the saturation adjustment
                        ! after nudging updating at time level nnew
    nzerr  , nstat      ! status of memory (de-)allocation

  REAL    (KIND=wp   )     ::  &
    fmultw           ,& ! 0. or 1., depending on multiple observation weighting
    zdt              ,& ! timestep relevant for the nudging
    zgnudgt, zgnudgq ,& ! basic nudging coefficient for temp., rel. humidity
    zrtv             ,& ! R * virtual temperature
    zqmy                ! total nudging weight for spreaded humidity obs. incr.

  LOGICAL                  ::  &
    lgetgeo          ,& ! geostrophic wind correct. is computed at current time
    lwrite              ! printout for control

  CHARACTER(LEN=255)       :: yzerrmsg
  CHARACTER(LEN=25)        :: yzroutine

! Local (automatic) arrays:
! -------------------------

  REAL    (KIND=wp   )     , ALLOCATABLE , SAVE ::       &
    paninc    (:,:)  ,& ! analysis increment of (upper-air) pressure
                        ! (hydrostatically derived)
    taninc    (:,:)  ,& ! analysis incr. of temperature (incl. condens./evapor.)
    qaninc    (:,:)  ,& ! analysis incr. of humidity (in the form of the density
                        ! increment due to moisture, condens./evapor. included)
!   taieva    (:,:)  ,& ! analysis incr. of temperature (excl. condens./evapor.)
    qaieva    (:,:)     ! analysis incr. of humidity (in the form of the density
                        ! increment due to moisture, condens./evapor. excluded)

  REAL    (KIND=wp   )     ::  &
    zwsum   (ie,je,3),& ! sum of weights assigned to the different obs. types
    taninb  (ie,je)  ,& ! analysis increments of temperature (incl. condens.
                        ! /evapor.) at level k+1
    qaninb  (ie,je)  ,& ! analysis increments of humidity at level k+1 (in the
                        ! form of the correction of density due to moisture, 
                        ! condens./evapor. included)
    zpkl    (ie,je)  ,& ! pressure
    zrhmod  (ie,je)  ,& ! relative humidity (model value)
    zewsat  (ie,je)  ,& ! saturation vapour pressure
    zrhni   (ie,je)  ,& ! net spreaded observat. increment of relative humidity
    qanait  (ie,je)  ,& ! specific humidity increment due to nudging temperature
    qaigeo  (ie,je)  ,& ! specific humidity increment due to temperature correc.
    zfbuoyt (ie,je)  ,& ! factor in the buoyancy term in the w-equation
    zfbuoyb (ie,je)  ,& ! as 'zfbuoyt' but for a model level further below

    zqce    (ie,je)  ,& ! Auxiliary fields for subr. *satad*
    zphfe   (ie,je)  ,& !      |
    zr1     (ie,je)  ,& !      | 
    zr2     (ie,je)  ,& !    \ | /
    zr3     (ie,je)  ,& !     \|/
    zr4     (ie,je)  ,& !      V
    zr5     (ie,je)  ,& !
    zr6     (ie,je)  ,& !
    zr7     (ie,je)  ,& !
    zr8     (ie,je)     !

  REAL    (KIND=wp   )     , ALLOCATABLE , SAVE ::       &
    zxmax       (:)     ! max. analysis increments of mass fields

  INTEGER (KIND=iintegers) , ALLOCATABLE , SAVE ::       &
    ijxmx     (:,:)     ! coord. of max. analysis increments of mass fields

! Tracer pointers:
! -----------------
  REAL (KIND=wp)     , POINTER  :: &
    qv(:,:,:) => NULL()        ,& ! QV at nnew
    qc(:,:,:) => NULL()           ! QC at nnew
!
!------------ End of header ----------------------------------------------------


 
!-------------------------------------------------------------------------------
! Begin Subroutine nudge_humid_mass
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 0: Preliminaries (e.g. prepar. storage of analysis increments)
!-------------------------------------------------------------------------------

izerror   = 0_iintegers
yzerrmsg  = ''
yzroutine = 'nudge_humid_mass'

! Retrieve the required microphysics tracers (at nnew)
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnew, ptr = qv)
  IF (izerror /= 0_iintegers) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nnew, ptr = qc)
  IF (izerror /= 0_iintegers) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  lgetgeo  =  (luvgcor) .AND. (lgetai)

  IF (l2tls) THEN
    zdt = dt
  ELSE
    zdt = dt2
  ENDIF
 
  IF (k == ke) THEN
    IF (ALLOCATED( paninc )) THEN
      PRINT '("CAUTION in src_nudging: paninc is already allocated "           &
            &,"at time / PE",I6,I5)', ntstep, my_cart_id
      DEALLOCATE ( paninc , taninc , qaninc , STAT=nzerr )
    ENDIF
    ALLOCATE ( paninc(ie,je) , stat=nzerr )
    ALLOCATE ( taninc(ie,je) , stat=nzerr )
    ALLOCATE ( qaninc(ie,je) , stat=nzerr )
    paninc = c0
    taninc = c0
    qaninc = c0
  ENDIF
  IF ((k == ke) .AND. (lgetgeo)) THEN
    IF (ALLOCATED( qaieva )) THEN
      PRINT '("CAUTION in src_nudging: qaieva is already allocated "           &
            &,"at time / PE",I6,I5)', ntstep, my_cart_id
      DEALLOCATE ( qaieva , STAT=nzerr )
    ENDIF
!   ALLOCATE ( taieva(ie,je) , stat=nzerr )
    ALLOCATE ( qaieva(ie,je) , stat=nzerr )
!   taieva = c0
    qaieva = c0
  ENDIF

  kwi = k
  IF (.NOT. lconai) kwi = 1


  IF (k < ke) THEN
    DO   j = jstart, jend
      DO i = istart, iend
        taninb (i,j) = taninc (i,j)
        qaninb (i,j) = qaninc (i,j)
      ENDDO
    ENDDO
  ENDIF
  DO   j = jstart, jend
    DO i = istart, iend
      taninc (i,j) = t (i,j,k,nnew)
      qaninc (i,j) = rvd_m_o*qv(i,j,k) - qc(i,j,k) - qrs(i,j,k)
    ENDDO
  ENDDO


! prepare storage of density incr. for computing the geostrophic wind correction

  IF (lgetgeo) THEN
    DO   j = jstart, jend
      DO i = istart, iend
        zroi(i,j,k) =  (p0(i,j,k) + pp(i,j,k,nnew))                            &
                     / (r_d * t(i,j,k,nnew) * (c1 + qaninc(i,j)))
      ENDDO
    ENDDO
  ENDIF

!-------------------------------------------------------------------------------
!  Section 1:  Addition of analysis increments to pressure at lowest model level
!-------------------------------------------------------------------------------
 
  IF (k == ke) THEN
    DO   j = jstart, jend
      DO i = istart, iend
        paninc (i,j)     =  psanai(i,j)    + psaigeo(i,j)
        pp (i,j,k,nnew)  =  pp(i,j,k,nnew) + paninc(i,j)
        p_anai (i,j,k)   =  p_anai(i,j,k)  + paninc(i,j) - psaigeo(i,j)
      ENDDO
    ENDDO
  ENDIF
  DO   j = jstart, jend
    DO i = istart, iend
      zpkl (i,j) = p0(i,j,k) + pp(i,j,k,nnew)
    ENDDO
  ENDDO
 
 
!-------------------------------------------------------------------------------
!  Section 2 : Updating of temperature
!-------------------------------------------------------------------------------
!  Section 2a: Computation of a
!              1):    squared weighted mean   =    net (spreaded) obs. increment
!                     of obs. increments           at the target grid pt.
!              2):    weighted mean           =    net nudging weight
!                     of nudging weights           at the target grid pt.
!-------------------------------------------------------------------------------
 
  IF (lgetai) THEN

    lwrite  =  (lwonl) .AND. (k == ke-5) .AND. (ntstep <= 5)
    IF (lwrite) WRITE( nupr,'(''omt: weight/ sqr/ weighted incr'',I3,3F9.4)' ) &
                       k, omy(ionl,jonl,2,1), om2(ionl,jonl,2,1)               &
                        , zwi(ionl,jonl,3,1) 
! the MAX(.,epsy) is required for nudging of surface-level data,
! if the nudging coeff. for TEMP / PILOT is set to zero
    zgnudgt = MAX( gnudg(3) , epsy )

    zwsum  =  c0
    DO ityp = 1 , nwtyp
      fmultw = REAL( kwtyp(ityp) - 1, wp )
      DO   j = jstart, jend
        DO i = istart, iend
          IF ((om2(i,j,2,ityp) > epsy) .AND. (omy(i,j,2,ityp) > epsy)) THEN
! averaged increment for observation type 'ityp'
            zwi (i,j,3,ityp)  =   zwi(i,j,3,ityp)                              &
                               / (om2(i,j,2,ityp) + fmultw *omy(i,j,2,ityp))
! individual weight 'w(ityp)' assigned to the observation type 'ityp'
            omy (i,j,4,ityp)  =  (om2(i,j,2,ityp) + fmultw *omy(i,j,2,ityp))   &
                               / (omy(i,j,2,ityp) + fmultw)
! sum of the individual weights 'w(ityp)' assigned to the different obs. types
            zwsum (i,j,1)     =  omy(i,j,4,ityp)  +  zwsum(i,j,1)
          ELSE
            zwi (i,j,3,ityp)  =  c0
            omy (i,j,4,ityp)  =  c0
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    fmultw = REAL( kwtyp(1) - 1, wp )
    IF (nwtyp >= 2)  fmultw = REAL( kwtyp(nwtyp+1) - 1, wp )
    DO ityp = 1 , nwtyp
      DO   j = jstart, jend
        DO i = istart, iend
          IF (zwsum(i,j,1) > epsy) THEN
! total multiple weight 'W(ityp)' assigned to observation type 'ityp'
            omy (i,j,2,ityp)  =  omy(i,j,4,ityp) * (omy(i,j,4,ityp) + fmultw)  &
                                                 / (zwsum(i,j,1)    + fmultw)
! weighted averaged increment for observation type 'ityp'
            zwi (i,j,3,ityp)  =  zwi(i,j,3,ityp) * omy(i,j,2,ityp)
! sum of total multiple weights 'W(ityp)' assigned to the different obs. types
            zwsum (i,j,2)     =  omy(i,j,2,ityp)  +  zwsum(i,j,2)
          ELSE
            omy (i,j,2,ityp)  =  c0
            zwi (i,j,3,ityp)  =  c0
          ENDIF
! sum of the weighted averaged increments for the different observation types
          zwsum (i,j,3)       =  zwi(i,j,3,ityp)  +  zwsum(i,j,3)
        ENDDO
      ENDDO
    ENDDO
    DO   j = jstart, jend
      DO i = istart, iend
        IF (zwsum(i,j,2) > epsy) THEN
! double-averaged total observation increment
          zwi (i,j,3,1)  =  zwsum(i,j,3)  / zwsum(i,j,2)
! final weight, including nudging coeff., effect of lateral boundary, timestep !
          omy (i,j,2,1)  =  zwsum(i,j,2)  * zgnudgt * faclbc(i,j,1)  *  zdt
        ELSE
          zwi (i,j,3,1)  =  c0
          omy (i,j,2,1)  =  c0
        ENDIF
      ENDDO
    ENDDO

!   DO   j = jstart, jend
!     DO i = istart, iend
!       IF (om2(i,j,2,ityp) > epsy) THEN
!         zwi (i,j,3)  =  zwi(i,j,3) / om2(i,j,2)
! Inclusion of nudging coefficient, & of weight reduction near lateral boundary
!         omy (i,j,2)  =  om2(i,j,2) / omy(i,j,2) * zgnudgt * faclbc(i,j,1)
!       ELSE
!         zwi (i,j,3)  =  c0
!       ENDIF
!     ENDDO
!   ENDDO

! weighted mean of increments used to balance (specific / relative) humidity

!!  DO   j = jstart, jend
!!    DO i = istart, iend
!!      IF (om2(i,j,2,ityp) > epsy) THEN
!!        zwi (i,j,5)  =  zwi(i,j,5) / om2(i,j,2)
!!      ELSE
!!        zwi (i,j,5)  =  c0
!!      ENDIF
!!    ENDDO
!!  ENDDO

    IF (lwrite) WRITE( nupr,'(''omt: weight     / weighted incr'',I3,F12.7,6X  &
                            &,F9.4)' ) k, omy(ionl,jonl,2,1), zwi(ionl,jonl,3,1)
    IF ((lwonl) .AND. (k == ke-5) .AND. (ntstep <= 40))                        &
      WRITE( nupr,'(''ntstep='',I4,'', k='',I2,'', T / qv / qc :''             &
                  &,F7.2,2F9.6,2I4,'', ztwi / omyt:'',F7.2,F8.5)' )            &
             ntstep, k, t(ionl,jonl,k,nnew), qv(ionl,jonl,k)                   &
           , qc(ionl,jonl,k), ionl,jonl, zwi(ionl,jonl,3,1)                    &
                                            , omy(ionl,jonl,2,1)

  ENDIF

 
!-------------------------------------------------------------------------------
!  Section 2b: Addition of nudging terms to prognostic temperature
!              (i.e. computation of temperature analysis increments, updating
!               the prognostic temperature field by adding these increments)
!-------------------------------------------------------------------------------
 
  IF (lgetai) THEN
    DO   j = jstart, jend
      DO i = istart, iend
! analysis increment of temperature,
! condensation / evaporation not yet taken into account
        tanai (i,j,kwi)  =  omy(i,j,2,1) / (c1 + omy(i,j,2,1)) * zwi(i,j,3,1)  &
                          +     c1       / (c1 + omy(i,j,2,1))                 &
!!                                       / (c1 + fswtp*om2(i,j,2,1))           &
                                         * (ztwips(i,j,k) + ztpgeo(i,j,k))
      ENDDO
    ENDDO
    DO   j = jstart, jend
      DO i = istart, iend
! temperature analysis increment used to balance humidity
! (conditional preservation of relative humidity)
!       tairh (i,j,kwi)  =  omy(i,j,2,1) / (c1 + omy(i,j,2,1)) * zwi(i,j,5,1)
! temperature analysis increment due to pressure nudging
! (unconditional preservation of relative humidity)
          taips (i,j,kwi) =     c1       / (c1 + omy(i,j,2,1)) *(  ztwips(i,j,k) &
                                                                 + ztpgeo(i,j,k) )
      ENDDO
    ENDDO
  ENDIF

! new (nudged) value of 'T', 'qv' (before condensation / evaporation)
  DO   j = jstart, jend
    DO i = istart, iend
      t (i,j,k,nnew)  =  t (i,j,k,nnew) + tanai(i,j,kwi)
      t_anai (i,j,k)  =  t_anai(i,j,k)  + tanai(i,j,kwi)
    ENDDO
  ENDDO

 
!-------------------------------------------------------------------------------
!  Section 3:  Updating of the specific humidity so that the relative or,
!              depending on convective activity, specific humidity is as before
!              nudging the temperature (if 'khumbal' < 99);
!              preparation of the specific humidity increment due to the
!              temperature correction (if 'khumbal' < 100).
!-------------------------------------------------------------------------------

  IF (lgetgeo)         qaigeo = c0
  IF (khumbal >= 100)  qanait = c0

  IF (khumbal < 100) THEN
    DO   j = jstart, jend
      DO i = istart, iend
! specific humidity prior to the nudging of temperature
        qanait (i,j)  = qv(i,j,k)
! relative humidity prior to (and after) the nudging of (part of) temperature
        zrhmod (i,j)  =   fq2pv( qv(i,j,k), zpkl(i,j), rdv )                   &
                        / fpvsw( taninc(i,j), b1, b2w, b3, b4w )
      ENDDO
    ENDDO
  ENDIF
  IF ((lgetgeo) .AND. (khumbal < 100)) THEN
    DO   j = jstart, jend
      DO i = istart, iend
! specific humidity increment due to the temperature correction
        IF (ABS( psanai(i,j) ) > epsy) THEN
          qaigeo (i,j) = fpv2q( fpvsw( taninc(i,j) + ztwips(i,j,k)             &
                                     , b1, b2w, b3, b4w ) * zrhmod(i,j)        &
                              , zpkl(i,j) , rdv )  -  qv(i,j,k)
          qaigeo (i,j) = qaigeo(i,j) * rvd_m_o
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  IF (khumbal <= -1) THEN
    DO   j = jstart, jend
      DO i = istart, iend
! saturation vapour pressure after the nudging of (part of the) temperature
        zewsat  (i,j)  =  fpvsw( taninc(i,j) + tanai(i,j,kwi), b1, b2w, b3,b4w )
! specific humidity after the nudging of temperature
        qv (i,j,k)     =  fpv2q( zewsat(i,j) * zrhmod(i,j) , zpkl(i,j) , rdv )
! specific humidity increment due to the nudging of temperature
        qanait  (i,j)  =  qv(i,j,k) - qanait(i,j)
      ENDDO
    ENDDO
  ELSEIF (khumbal < 99) THEN
    DO   j = jstart, jend
      DO i = istart, iend
! saturation vapour pressure after the nudging of (part of the) temperature
        IF (      (hhl(i,j,k  ) > zconbas(i,j))                                &
            .AND. (hhl(i,j,k+1) < zcontop(i,j))) THEN
          zewsat  (i,j)  =  fpvsw( taninc(i,j) + taips(i,j,kwi), b1,b2w,b3,b4w )
        ELSE
          zewsat  (i,j)  =  fpvsw( taninc(i,j) + tanai(i,j,kwi), b1,b2w,b3,b4w )
        ENDIF
! specific humidity after the nudging of temperature
        qv (i,j,k)       =  fpv2q( zewsat(i,j) * zrhmod(i,j) , zpkl(i,j) , rdv )
! specific humidity increment due to the nudging of temperature
        qanait    (i,j)  =  qv(i,j,k) - qanait(i,j)
      ENDDO
    ENDDO
  ENDIF


!-------------------------------------------------------------------------------
!  Section 4:  Computation of a net (spreaded) obs. increments of relative
!              humidity and of net nudging weights at the target grid pt.
!              (cf. section 2a).
!              Approximate conversion of the net obs. increments of relative
!              humidity into net obs. increments of specific humidity.
!-------------------------------------------------------------------------------
 
  IF (lgetai) THEN

!   lwrite  =  (lwonl) .AND. (k == ke-5) .AND. (ntstep <= 5)
    IF (lwrite) WRITE( nupr,'(''omq: weight/ sqr/ weighted incr'',I3,3F9.4)' ) &
                       k, omy(ionl,jonl,3,1), om2(ionl,jonl,3,1)               &
                        , zwi(ionl,jonl,4,1) 
! the MAX(.,epsy) is required for nudging of surface-level data,
! if the nudging coeff. for TEMP / PILOT is set to zero
    zgnudgq = MAX( gnudg(4) , epsy )

    zwsum  =  c0
    DO ityp = 1 , nwtyp
      fmultw = REAL( kwtyp(ityp) - 1, wp )
      DO   j = jstart, jend
        DO i = istart, iend
          IF ((om2(i,j,3,ityp) > epsy) .AND. (omy(i,j,3,ityp) > epsy)) THEN
! averaged relative humidity increment for observation type 'ityp'
            zwi (i,j,4,ityp)  =   zwi(i,j,4,ityp)                              &
                               / (om2(i,j,3,ityp) + fmultw *omy(i,j,3,ityp))
! individual weight 'w(ityp)' assigned to the observation type 'ityp'
            omy (i,j,4,ityp)  =  (om2(i,j,3,ityp) + fmultw *omy(i,j,3,ityp))   &
                               / (omy(i,j,3,ityp) + fmultw)
! sum of the individual weights 'w(ityp)' assigned to the different obs. types
            zwsum (i,j,1)     =  omy(i,j,4,ityp)  +  zwsum(i,j,1)
          ELSE
            zwi (i,j,4,ityp)  =  c0
            omy (i,j,4,ityp)  =  c0
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    DO ityp = 1 , nwtyp
      fmultw = REAL( kwtyp(1) - 1, wp )
      IF (nwtyp >= 2)  fmultw = REAL( kwtyp(nwtyp+1) - 1, wp )
      DO   j = jstart, jend
        DO i = istart, iend
          IF (zwsum(i,j,1) > epsy) THEN
! total multiple weight 'W(ityp)' assigned to observation type 'ityp'
            omy (i,j,3,ityp)  =  omy(i,j,4,ityp) * (omy(i,j,4,ityp) + fmultw)  &
                                                 / (zwsum(i,j,1)    + fmultw)
! weighted averaged relative humidity increment for observation type 'ityp'
            zwi (i,j,4,ityp)  =  zwi(i,j,4,ityp) * omy(i,j,3,ityp)
! sum of total multiple weights 'W(ityp)' assigned to the different obs. types
            zwsum (i,j,2)     =  omy(i,j,3,ityp)  +  zwsum(i,j,2)
          ELSE
            omy (i,j,3,ityp)  =  c0
            zwi (i,j,4,ityp)  =  c0
          ENDIF
! sum of the weighted averaged increments for the different observation types
          zwsum (i,j,3)       =  zwi(i,j,4,ityp)  +  zwsum(i,j,3)
        ENDDO
      ENDDO
    ENDDO
    DO   j = jstart, jend
      DO i = istart, iend
        IF (zwsum(i,j,2) > epsy) THEN
! double-averaged total observation increment  =  net increment of rel. humidity
          zrhni (i,j)    =  zwsum(i,j,3)  / zwsum(i,j,2)
! net nudging weight, incl. nudging coeff, effect of lateral boundary, timestep!
          omy (i,j,3,1)  =  zwsum(i,j,2)  * zgnudgq * faclbc(i,j,1)  *  zdt
        ELSE
          zrhni (i,j)    =  c0
          omy (i,j,3,1)  =  c0
        ENDIF
      ENDDO
    ENDDO

    DO   j = jstart, jend
      DO i = istart, iend
        IF (omy(i,j,3,1) > epsy) THEN
!       IF (om2(i,j,3)   > epsy) THEN
! Net obs. increment of relative humidity
!         zrhni (i,j)    =  zwi(i,j,4) / om2(i,j,3)
! Net obs. increment of vapour pressure (d(RH) * new saturation vapour pressure)
          zwi (i,j,4,1)  =  zrhni(i,j) * fpvsw( t(i,j,k,nnew), b1, b2w, b3,b4w )
! Bogus observation of vapour pressure
          zwi (i,j,4,1)  =  zwi(i,j,4,1) + fq2pv( qv(i,j,k), zpkl(i,j), rdv )
! Bogus observation of specific humidity (approximation)
          zwi (i,j,4,1)  =                 fpv2q( zwi(i,j,4,1), zpkl(i,j), rdv )
! Net obs. increment of specific humidity (approximation)
          zwi (i,j,4,1)  =  zwi(i,j,4,1)  -  qv(i,j,k)
        ELSE
          zwi (i,j,4,1)  =  c0
!         zrhni (i,j)    =  c0
        ENDIF
      ENDDO
    ENDDO

!   DO   j = jstart, jend
!     DO i = istart, iend
!       IF (om2(i,j,3) > epsy)                                                 &
! Net nudging weight, with inclusion of nudging coefficient,
!                                and of weight reduction near lateral boundaries
!         omy (i,j,3) = om2(i,j,3) / omy(i,j,3) * zgnudgq * faclbc(i,j,1)
!     ENDDO
!   ENDDO

    IF (lwrite) WRITE( nupr,'(''omq: weight     / weighted incr'',I3,F12.7,6X  &
                            &,F9.4)' ) k, omy(ionl,jonl,3,1), zwi(ionl,jonl,4,1)

  ENDIF

 
!-------------------------------------------------------------------------------
!  Section 5:  Analysis increments of specific humidity,
!              updating the prognostic humidity field by adding these increments
!              - Further preparations
!-------------------------------------------------------------------------------
 
  IF (lgetai) THEN
    DO   j = jstart, jend
      DO i = istart, iend
! analysis increment of specific water content 'delta-q' = 'delta-qv'
        zqmy            = omy(i,j,3,1) / (c1 + omy(i,j,3,1))
        qanai (i,j,kwi) = zqmy * zwi(i,j,4,1)
        zrhni (i,j)     = zqmy * zrhni(i,j)
      ENDDO
    ENDDO
  ENDIF

! new (nudged) value of 'qv' (before condensation / evaporation)
  DO   j = jstart, jend
    DO i = istart, iend
      qv (i,j,k)       =  qv(i,j,k)       +  qanai(i,j,kwi)
      qv_anai (i,j,k)  =  qv_anai(i,j,k)  +  qanai(i,j,kwi) + qanait(i,j)
    ENDDO
  ENDDO

! preparation of computation of upper-air analysis increments of pressure

  IF (k < ke) THEN
    DO   j = jstart, jend
      DO i = istart, iend
        zfbuoyt (i,j)  =  r_d * rho0(i,j,k  ) * t(i,j,k  ,nnow)
        zfbuoyb (i,j)  =  r_d * rho0(i,j,k+1) * t(i,j,k+1,nnow)
      ENDDO
    ENDDO
  ENDIF

! do / prepare the storage of the pressure / density analysis increments excl.
! condensation + evaporation for the comput. of the geostrophic wind correction

  IF (lgetgeo) THEN
    IF (k < ke) THEN
      kp1 = k + 1
      DO   j = jstart, jend
        DO i = istart, iend
!         zpai (i,j,k) =          c1       / (c2 + dp0(i,j,k+1) /zfbuoyt(i,j)) &
!                       *(   zpai(i,j,kp1) * (c2 - dp0(i,j,k  ) /zfbuoyb(i,j)) &
!                         +  tanai (i,j,kwi) * dp0(i,j,k+1) / t(i,j,k  ,nnow)  &
!                         +  taieva(i,j)     * dp0(i,j,k  ) / t(i,j,k+1,nnow)  &
!                         + (qanai (i,j,kwi) + qanait(i,j)) * rvd_m_o          &
!                                            * dp0(i,j,k+1)                    &
!                         +  qaieva(i,j)     * dp0(i,j,k  )  )
          zpai (i,j,k) =          c1       / (c2 + dp0(i,j,kp1) /zfbuoyt(i,j)) &
                        *(   zpai(i,j,kp1) * (c2 - dp0(i,j,k  ) /zfbuoyb(i,j)) &
                          +  ztwips(i,j,k  ) * dp0(i,j,kp1) / t(i,j,k  ,nnow)  &
                          +  ztwips(i,j,kp1) * dp0(i,j,k  ) / t(i,j,kp1,nnow)  &
                          +  qaigeo(i,j)     * dp0(i,j,kp1)                    &
                          +  qaieva(i,j)     * dp0(i,j,k  )  )
        ENDDO
      ENDDO
    ENDIF
    DO   j = jstart, jend
      DO i = istart, iend
        qaieva (i,j) =  qaigeo(i,j)
!       qaieva (i,j) = (qanai(i,j,kwi) + qanait(i,j)) * rvd_m_o
!       taieva (i,j) =  tanai(i,j,kwi)
!       zrtv         =  r_d * t(i,j,k,nnew) * (c1 + qaninc(i,j) + qaieva(i,j))
        zrtv         =  r_d * (taninc(i,j) + ztwips(i,j,k))                    &
                                          * (c1 + qaninc(i,j) + qaieva(i,j))
        zroi (i,j,k) =  (zpkl(i,j) + zpai(i,j,k)) /zrtv  -  zroi(i,j,k)
      ENDDO
    ENDDO
    IF (k == ke) THEN
      DO   j = jstart, jend
        DO i = istart, iend
          zpai (i,j,k)  =  psanai(i,j)
        ENDDO
      ENDDO
    ENDIF
  ENDIF


!-------------------------------------------------------------------------------
!  Section 6:  Saturation adjustment for the updated variables
!-------------------------------------------------------------------------------

IF (lcond) THEN

 IF (itype_gscp < 100 .OR. l2mom_satads) THEN

  kitera       = 1

  IF ( ldiabf_satad ) THEN
    ! initialize temperature increment due to latent heat
    tinc_lh(:,:,k) = tinc_lh(:,:,k) - t(:,:,k,nnew)
  ENDIF

  ! store T before saturation adjustment for calculation of latent heating rate
  IF (llhn .OR. llhnverif) CALL get_gs_lheating ('add',k,k)
!                          ====================

  DO   j = jstart, jend
    DO i = istart, iend
      zphfe  (i,j) =  p0(i,j,k) + pp(i,j,k,nnew)
      zqce   (i,j) =  qc(i,j,k)
    ENDDO
  ENDDO

  CALL satad ( kitera, t(:,:,k,nnew), qv(:,:,k), zqce, t(:,:,k,nnow)           &
             , zphfe, zr1, zr2, zr3, zr4, zr5, zr6, zr7, zr8                   &
             , b1, b2w, b3, b4w, b234w, rdv, o_m_rdv, rvd_m_o, lh_v,cpdr,cp_d  &
             , ie, je, istart, iend, jstart, jend  )
! ==========

  DO j = jstart, jend
    DO i = istart, iend
      qc_anai (i,j,k)  =  qc_anai(i,j,k) + zqce(i,j) - qc(i,j,k)
      qc (i,j,k)       =  zqce(i,j)
    ENDDO
  ENDDO

  IF ( ldiabf_satad ) THEN
    ! compute temperature increment due to latent heat
    tinc_lh(:,:,k) = tinc_lh(:,:,k) + t(:,:,k,nnew)
  ENDIF

! calculate gridscale latent heating rate from the saturation adjustment
  IF (llhn .OR. llhnverif) CALL get_gs_lheating ('inc',k,k)
!                          ====================

 ENDIF

ENDIF
 
!-------------------------------------------------------------------------------
!  Section 6b: Storage of analysis increments of temperature and humidity
!              (in the form of the correction of density due to moisture)
!              in which condensation / evaporation is taken into account.
!              Printout for control.
!-------------------------------------------------------------------------------

  DO   j = jstart, jend
    DO i = istart, iend
      taninc (i,j) = t(i,j,k,nnew) - taninc(i,j)
      qaninc (i,j) = rvd_m_o*qv(i,j,k) - qc(i,j,k) - qrs(i,j,k) - qaninc(i,j)
    ENDDO
  ENDDO

  IF ((lwonl) .AND. (k == ke-5) .AND. (ntstep <= 40))                          &
    WRITE( nupr,'(''ntstep='',I4,'', k='',I2,'', T / qv / qc :''               &
                &,F7.2,2F9.6,2I4,'', lgetai: '',L1)' )                         &
           ntstep, k, t(ionl,jonl,k,nnew), qv(ionl,jonl,k)                     &
         , qc(ionl,jonl,k), ionl, jonl, lgetai
  IF ((lwonl2) .AND. (k == ke-5) .AND. (ntstep*dtdeh <= c2))                   &
    WRITE( nupr,'(''Tqnudge: T,q-ana.incr. '',3I4, F9.4, 2F10.6)' )            &
           ionl2, jonl2, k, tanai(ionl2,jonl2,k), qanait(ionl2,jonl2)          &
                          , qanai(ionl2,jonl2,k)

!-------------------------------------------------------------------------------
!  Section 7:  - Analysis increments of pressure, hydrostatically derived from
!                1. pressure analysis increments at adjacent model level below
!                2. analysis increments of temperature and humidity at current
!                   level and at adjacent level below.
!              - Addition of analysis increments to prognostic model pressure
!              - Saturation adjustment not applied
!-------------------------------------------------------------------------------

  IF (k < ke) THEN

    DO   j = jstart, jend
      DO i = istart, iend

        paninc (i,j) =          c1      / (c2 + dp0(i,j,k+1) /zfbuoyt(i,j))    &
                      * (   paninc(i,j) * (c2 - dp0(i,j,k  ) /zfbuoyb(i,j))    &
                         +  taninc(i,j) * dp0(i,j,k+1) / t(i,j,k  ,nnow)       &
                         +  taninb(i,j) * dp0(i,j,k  ) / t(i,j,k+1,nnow)       &
                         +  qaninc(i,j) * dp0(i,j,k+1)                         &
                         +  qaninb(i,j) * dp0(i,j,k  )  )

        pp (i,j,k,nnew)  =  pp(i,j,k,nnew) + paninc(i,j)
        p_anai (i,j,k)   =  p_anai(i,j,k)  + paninc(i,j)

      ENDDO
    ENDDO

    IF ((lcond) .AND. (lpsatad)) THEN
      kitera       = 1
      IF ( ldiabf_satad ) THEN
        ! initialize temperature increment due to latent heat
        tinc_lh(:,:,k) = tinc_lh(:,:,k) - t(:,:,k,nnew)
      ENDIF
      ! store T before saturation adjustment for calculating latent heating rate
      IF (llhn .OR. llhnverif) CALL get_gs_lheating ('add',k,k)
!                              ====================
      DO   j = jstart, jend
        DO i = istart, iend
          zphfe  (i,j) =  p0(i,j,k) + pp(i,j,k,nnew)
          zqce   (i,j) =  qc(i,j,k)
        ENDDO
      ENDDO
      CALL satad ( kitera, t(:,:,k,nnew), qv(:,:,k), zqce, t(:,:,k,nnow)       &
                 , zphfe, zr1, zr2, zr3, zr4, zr5, zr6, zr7, zr8 , b1, b2w     &
                 , b3, b4w, b234w, rdv, o_m_rdv, rvd_m_o, lh_v, cpdr, cp_d     &
                 , ie, je, istart, iend, jstart, jend  )
!     ==========
      DO j = jstart, jend
        DO i = istart, iend
          qc_anai (i,j,k)  =  qc_anai(i,j,k) + zqce(i,j) - qc(i,j,k)
          qc (i,j,k)       =  zqce(i,j)
        ENDDO
      ENDDO
      IF ( ldiabf_satad ) THEN
        ! compute temperature increment due to latent heat
        tinc_lh(:,:,k) = tinc_lh(:,:,k) + t(:,:,k,nnew)
      ENDIF
      ! calculate gridscale latent heating rate from the saturation adjustment
      IF (llhn .OR. llhnverif) CALL get_gs_lheating ('inc',k,k)
!                              ====================
    ENDIF

  ENDIF

! store pressure & density increments for comput. the geostrophic wind correct.

! IF (lgetgeo) THEN
!   DO   j = jstart, jend
!     DO i = istart, iend
!       zroi (i,j,k)  = (p0(i,j,k) + pp(i,j,k,nnew)) / zrtv(i,j)  -  zroi(i,j,k)
!!      zpai (i,j,k)  =  paninc(i,j)
!     ENDDO
!   ENDDO
! ENDIF

! store analysis increments for diagnostic output

  IF (lwonl) THEN
    zpaionl = paninc (ionl,jonl)
    ztaionl = taninc (ionl,jonl)
    zqaionl = qaninc (ionl,jonl)
  ENDIF

! control output of max. analysis increments

  IF (lgetai) THEN
    IF (k == ke) THEN
      IF (ALLOCATED( zxmax )) THEN
        PRINT '("CAUTION in src_nudging: zxmax is already allocated "          &
              &,"at time / PE",I6,I5)', ntstep, my_cart_id
        DEALLOCATE ( zxmax , ijxmx , STAT=nzerr )
      ENDIF
      ALLOCATE ( zxmax (  12*ke) , stat=nzerr )
      ALLOCATE ( ijxmx (2,12*ke) , stat=nzerr )
    ENDIF
    ixe = k + 11*ke
    DO  ix = k , ixe , ke
      zxmax   (ix) = c0
      ijxmx (1,ix) = 0
      ijxmx (2,ix) = 0
    ENDDO
    DO   j = jstart, jend
!for the NEC
!CDIR NODEP
      DO i = istart, iend
        ixe = k
        IF (ABS(ztwips(i,j,k)) > ABS(zxmax(ixe))) THEN
          zxmax   (ixe) = ztwips(i,j,k)
          ijxmx (1,ixe) = i
          ijxmx (2,ixe) = j
        ENDIF
        ixe = k + ke
        IF (ABS(tanai (i,j,kwi)) > ABS(zxmax(ixe))) THEN
          zxmax   (ixe) = tanai (i,j,kwi)
          ijxmx (1,ixe) = i
          ijxmx (2,ixe) = j
        ENDIF
        ixe = k + 2*ke
        IF (ABS(qanait(i,j)) > ABS(zxmax(ixe))) THEN
          zxmax   (ixe) = qanait(i,j)
          ijxmx (1,ixe) = i
          ijxmx (2,ixe) = j
        ENDIF
        ixe = k + 3*ke
        IF (ABS(qanai (i,j,kwi)) > ABS(zxmax(ixe))) THEN
          zxmax   (ixe) = qanai (i,j,kwi)
          ijxmx (1,ixe) = i
          ijxmx (2,ixe) = j
        ENDIF
        ixe = k + 4*ke
        IF (ABS(zrhni (i,j)) > ABS(zxmax(ixe))) THEN
          zxmax   (ixe) = zrhni (i,j)
          ijxmx (1,ixe) = i
          ijxmx (2,ixe) = j
        ENDIF
        ixe = k + 5*ke
        IF (ABS(paninc(i,j)) > ABS(zxmax(ixe))) THEN
          zxmax   (ixe) = paninc(i,j)
          ijxmx (1,ixe) = i
          ijxmx (2,ixe) = j
        ENDIF
      ENDDO
    ENDDO

    IF ((num_compute > 1 ) .AND. (k == 1)) THEN
      DO  kk = 1 , ke
        DO  ix = kk , kk+5*ke , ke
          ixe = ix + 6*ke
          zxmax   (ixe) = - zxmax  (ix)
          ijxmx (1,ixe) =   ijxmx(1,ix)
          ijxmx (2,ixe) =   ijxmx(2,ix)
        ENDDO
      ENDDO
      ixe = 12 * ke
      IF (ltime) THEN
        CALL get_timings (i_nudging_corr, ntstep, dt, izerror)
        CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
        CALL get_timings (i_barrier_waiting_nud, ntstep, dt, izerror)
      ENDIF

      CALL global_values ( zxmax, ixe, 'MAX', imp_reals, icomm_cart            &
                         , onl_cart_id, ijxmx, yzerrmsg, izerror)
!     ==================

      CALL global_values ( nhisto, 15, 'SUM', imp_integers, icomm_cart         &
                         , onl_cart_id, yzerrmsg, izerror)
!     ==================

      IF (ltime) CALL get_timings (i_communications_nud, ntstep, dt,izerror)
      DO  kk = 1 , ke
        DO  ix = kk , kk+5*ke , ke
          ixe = ix + 6*ke
          IF (ABS(zxmax(ixe)) > ABS(zxmax(ix))) THEN
            zxmax   (ix) = - zxmax  (ixe)
            ijxmx (1,ix) =   ijxmx(1,ixe)
            ijxmx (2,ix) =   ijxmx(2,ixe)
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    zxmax (k+2*ke) = zxmax(k+2*ke) * 1000.0_wp
    zxmax (k+3*ke) = zxmax(k+3*ke) * 1000.0_wp
    zxmax (k+4*ke) = zxmax(k+4*ke) *  100.0_wp
    ixe = k + 5*ke
    IF ((lwonl) .AND. (     (ABS(zxmax(k     )) >  10.0_wp)                &
                       .OR. (ABS(zxmax(k+  ke)) >  10.0_wp)                &
                       .OR. (ABS(zxmax(k+2*ke)) >  10.0_wp)                &
                       .OR. (ABS(zxmax(k+3*ke)) >  10.0_wp)                &
                       .OR. (ABS(zxmax(k+4*ke)) > 100.0_wp)                &
                       .OR. (ABS(zxmax(k+5*ke)) > 500.0_wp)))              &
      WRITE( nupr,'(''Tq-AI'',I3,2(F7.0,2I4),2(F9.2,2I4),F6.0,2I4,F7.0,2I4)' ) &
             k   , (zxmax(ix), ijxmx(1,ix), ijxmx(2,ix) , ix=k,ixe,ke)

    IF ((lwonl) .AND. (k == 1)) THEN
      WRITE( nupr,'(''Analysis increments at timestep'',I5)' )  ntstep
      WRITE( nupr,'(''max. (ABS of) analysis increments of mass field'')')
      WRITE( nupr,'(''(prior to condensation / evaporation except for''        &
                  &,'' pressure increments)'')' )
      WRITE( nupr,'(''le- T-corr global T-incr  global qv-incr due  gl.''      &
                  &,''  qv-incr  global RH-incr global  p-incr  global'')' )
      WRITE( nupr,'(''vel   [K]  coord.    [K]  coord. to T-incr coord.''      &
                  &,''   [g/kg]  coord.    [%]  coord.    [Pa]  coord.'')' )
      DO kk = 1 , ke
        ixe = kk + 5*ke
        WRITE( nupr,'(I2,F7.3,2I4,F7.2,2I4, 2(F9.4,2I4), F7.2,2I4, F8.3,2I4)') &
               kk  , (zxmax(ix), ijxmx(1,ix), ijxmx(2,ix) , ix=kk,ixe,ke)
      ENDDO
      WRITE( nupr,'("data density histogram for surface pressure obs:")' )
      IF (ntstep <= NINT(tconbox/dt)+1)                                        &
        WRITE( nupr,'("class 0-1  1-2  2-3  3-4  4-5  5-6  6-7  7-8  8-9"      &
                    &," 9-10 10-15  <20  <30  <40  >40")' )
      WRITE( nupr,'("Hit:",10I5,I6,4I5)' ) nhisto
    ENDIF
    IF (k == 1) THEN
      DEALLOCATE ( zxmax , stat=nzerr )
      DEALLOCATE ( ijxmx , stat=nzerr )
    ENDIF
  ENDIF

! clean-up

  IF (k == 1) THEN
    DEALLOCATE ( paninc , stat=nzerr )
    DEALLOCATE ( taninc , stat=nzerr )
    DEALLOCATE ( qaninc , stat=nzerr )
  ENDIF
  IF ((k == 1) .AND. (lgetgeo)) THEN
!   DEALLOCATE ( taieva , stat=nzerr )
    DEALLOCATE ( qaieva , stat=nzerr )
  ENDIF

! flush YUPRINT file
! IF ((ldump_ascii) .AND. (lwonl) .AND. (lgetai)
  IF ((ldump_ascii) .AND. (lwonl) .AND. (lfirst)                               &
                    .AND. ((k == 1) .OR. (k == ke))) THEN
    CLOSE (nupr)
    OPEN  (nupr,FILE=yuprint,FORM='FORMATTED',STATUS='OLD'                     &
                            ,POSITION='APPEND',IOSTAT=nstat)
    IF (nstat /= 0) PRINT '("OPENING OF FILE yuprint FAILED, nudge_humid_mass")'
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure nudge_humid_mass
!-------------------------------------------------------------------------------

END SUBROUTINE nudge_humid_mass


!-------------------------------------------------------------------------------
!+ Module procedure for temperature correction balancing surface pressure incr.
!-------------------------------------------------------------------------------

SUBROUTINE ps_temperatur_corr ( nactio )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure computes gridded temperature increments
!   ('temperature correction') to account for density changes which are
!   (physically or statistically) consistent with the pressure analysis
!   increments at the lowest model level resulting from the assimilation of
!   surface pressure data.
!
! Method:
!   For the 1st type of temperature correction, upper-level pressure increments
!   are specified in terms of the 'surface-level' pressure analysis increments
!   by use of a pressure-pressure ('p-p') correlation function 'cpp.1' .
!   For the 2nd type of temperature correction, the 'p-p' correlation is derived
!   approximately from a geopotential-geopotential correlation which is similar
!   to the one that is used by Bell et al. (1996) of UKMO for the AC-scheme.
!     - cpp.1(k) = ck(k)^2 * EXP( (1.-ck(k)^3) /8.)
!     - cfi.2(k) = .5 *ck(k) *(1.+ck(k))     (as approx. to the funct. in AC-s.)
!       cpp.2(k) = cfi.2(k) * p(k) / p(ke)
!     where      ck(k) = (p(k) - ptpstop) / (p(ke) - ptpstop)   for p > ptpstop
!                ck(k) = 0.                                     for p < ptpstop
!   From the upper-air pressure increments, upper-air virtual temperature
!   increments are computed locally and hydrostatically assuming constant
!   density within model layers.
!   In order to split these increments into temperature and humidity increments
!   so that the relative humidity remains unchanged, a quasi-iterative procedure
!   is devised:
!    1.: temperature increments are computed from the virtual temperature
!        increments assuming constant specific humidity
!    2.: specific humidity increments are computed from the temperature incre-
!        ments assuming constant relative humdity
!    3.: a residual pressure increment at the top level of the temperature
!        correction is computed by bottom-up integration of the hydrostatic
!        equation using the pressure increments at the lowest model level and
!        the upper-air temperature and specific humidity increments
!    4.: the upper-air residual pressure increment is converted into a residual
!        pressure increment at the lowest model level (by applying an appropri-
!        ate factor, or, equivalently, by top-down integration of the hydro-
!        static equation using pressure increments only)
!    5.: from this residual pressure increment, temperature increments can be
!        computed (by use of the 'p-p' correlation) (1st iteration), and added
!        to the original temperature increments  -
!        this equivalent to multiplying the original temperature increments by a
!        factor : (analysis + residual pressure incr.) / (analysis press. incr.)
!
!   Note: - cpp.2 is an approximation to the function used in the AC-scheme.
!         - Using cpp.1(k), the resulting surface pressure - temperature cor-
!           relation is in reasonable agreement with the statistical correlation
!           found by Anderson (personal communication, 1996) for the mesoscale
!           part (wavenumbers n=45, n=55) of the ECMWF model. 
!
! Written by        :  Christoph Schraff, DWD  (original version: 24.06.97)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
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
    nactio   ! = 0  : for determination of 'ktop' only (top full level of the
             !        temperature correction), including 3D- hydrostatic
             !        pressure and virtual temperature
             ! = 1  : for computing the temperature correction using 'ktop';
             !        without any communication between PE's

! Local parameters: None
! ----------------
! Local scalars:
! -------------

  LOGICAL                  , SAVE   ::  &
    lfst_prtps = .TRUE. ! print output of temperature correction once

  INTEGER (KIND=iintegers) ::  &
    izerror          ,& ! error index
    i      , j       ,& ! loop indices in horizontal direction
    kk               ,& ! loop index in the vertical direction
    khp1   , kp1     ,& ! index of lower half resp. main level
    ihpa   , ihpb    ,& ! cycling half level indices
    ipr(2) , jpr(2)  ,& ! loop indices in horizontal direction for printing
    ilpr             ,& ! loop index
    nzerr  , nstat      ! status of memory (de-)allocation

  REAL    (KIND=wp   )     ::  &
    rdg              ,& ! r_d / g
    zck              ,& ! vertical coordinate used for geopotential correlat.
    zptop            ,& ! pressure level of top of temperature correction
    zrhmod           ,& ! relative humidity (model value)
    zpkemin          ,& ! for printing: minimum pressure in local domain
    zp1000           ,& ! for printing: for p(ke) == zp1000 then ps = 1000 hPa
    zdp1000          ,& ! for printing: minimum departure from 'zp1000'
    zprh             ,& ! for printing: height at (ipr,jpr,kk)
    zprp             ,& ! for printing: pressure at (ipr,jpr,kk)
    zprhp            ,& ! for printing: hydrostatic pressure at (ipr,jpr,kk)
    zprcz            ,& ! for printing: height correlation at (ipr,jpr,kk)
    zprdz            ,& ! for printing: height correction at (ipr,jpr,kk)
    zprdt            ,& ! for printing: normalized temperature correction
    zpaipr  (ke,2)      ! for printing: residual relative pressure increment

  LOGICAL                  ::  &
    lgotop           ,& !
    lprtps              ! printing for control

  CHARACTER (LEN=39)       ::  ypr(2)       ! output for control
  CHARACTER (LEN=255)      ::  yzerrmsg
  CHARACTER (LEN=25)       ::  yzroutine

  INTEGER (KIND=iintegers) , SAVE ::  &
    ktop                ! uppermost full level influenced by the
                        ! temperature correction

! Local (automatic) arrays:
! -------------------------

  REAL    (KIND = wp)       , ALLOCATABLE :: &
    ps_fld  (:,:)    ,& ! input surface pressure field
    ztcorr  (:,:,:)  ,& ! temperature corrction balancing surface pressure field
    ztvrhl  (:,:)    ,& ! virtual temperature at lower half level
    zpdenom (:,:)    ,& ! factor used for the geopotental correlation
    zpker   (:,:)    ,& ! reciproke of hydrostatic pressure at lowest full level
    zfbuoyt (:,:)    ,& ! factor in the buoyancy term in the w-equation
    zfbuoyb (:,:)    ,& ! as 'zfbuoyt' but for a model level further below
    zqitps  (:,:)    ,& ! humidity increment due to the temperature correction
                        ! (in the form of the density increment due to moisture
                        ! perturbations, condens./evapor. excluded)
    zqitpb  (:,:)    ,& ! as 'zqitps', but at the model layer below
    zpainc  (:,:)    ,& ! pressure change at the top of the atmosphere from the
                        ! surface pressure change without temperature correction
!   zfisai  (:,:,:)  ,& ! factor used for the conversion of the geopotential
                        ! correlation into height increments
    zcfi    (:,:)    ,& ! vertical correlation of geopotential
                        ! (as array for more convenient printing only)
    zdhl    (:,:,:)  ,& ! height increment (ntpscor >= 3) or pressure increment
                        ! at half levels
    zcfim   (:,:)    ,& ! as 'zcfi', but on full model levels
    zdhlm   (:,:)       ! as 'zdhl', but on full model levels

  REAL    (KIND=wp   )     , ALLOCATABLE , SAVE ::       &
    ztvrml  (:,:,:)  ,& ! virtual temperature at full levels
    zhyphl  (:,:,:)     ! hydrostatic pressure at half levels

! Tracer pointers:
! -----------------
  REAL (KIND=wp)     , POINTER  :: &
    qv(:,:,:) => NULL()         ,& ! QV at nnew
    qc(:,:,:) => NULL()            ! QC at nnew

!
!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine ps_temperatur_corr
!-------------------------------------------------------------------------------

izerror   = 0_iintegers
yzerrmsg  = ''
yzroutine = 'ps_temperatur_corr'


! Retrieve the required microphysics tracers (at nnew)
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnew, ptr = qv)
  IF (izerror /= 0_iintegers) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nnew, ptr = qc)
  IF (izerror /= 0_iintegers) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

!-------------------------------------------------------------------------------
!  Section 0: Determination of the top full level of the temperature correction
!-------------------------------------------------------------------------------

IF (nactio == 0) THEN

! ptpttop = MIN( ptpstop , p0(:,:,ke) - 10000.0_wp )

  rdg = r_d / g

  ktop  = ke + 1

  IF (ntpscor > 0) THEN

!   IF (ntpscor == 2) ntpscor = 4 ! use dfi-dfis -correlat., see *input_nudging*
!   IF (ntpscor == 3) ntpscor = 1 ! use dp-dps -correlation

! virtual temperature at full levels and hydrostatic pressure at half levels

    IF (ALLOCATED( ztvrml )) THEN
      PRINT '("CAUTION in src_nudging: ztvrml is already allocated "           &
            &,"at time / PE",I6,I5)', ntstep, my_cart_id
      DEALLOCATE ( ztvrml , zhyphl , STAT=nzerr )
    ENDIF
!    PRINT*,'DEBUG1',my_cart_id,istart,iend,jstart,jend,ke
    ALLOCATE ( ztvrml (istart:iend,jstart:jend,ke  ) , STAT = nzerr )
    ALLOCATE ( zhyphl (istart:iend,jstart:jend,ke+1) , STAT = nzerr )
!    PRINT*,'DEBUG2',my_cart_id,istart,iend,jstart,jend,ke,nzerr

    DO   j = jstart, jend
      DO i = istart, iend
        ztvrml (i,j,1) =  t(i,j,1,nnew) *( c1 + rvd_m_o*qv(i,j,1)              &
                                              - qc(i,j,1)-qrs(i,j,1))
        zhyphl (i,j,1) =  p0(i,j,1) + pp(i,j,1,nnew)
        zhyphl (i,j,2) = zhyphl(i,j,1) * (c1 + c05 *(hhl(i,j,1) - hhl(i,j,2))  &
                                                   /(rdg *ztvrml(i,j,1)))
        zhyphl (i,j,1) = c2 *zhyphl(i,j,1) - zhyphl(i,j,2)
      ENDDO
    ENDDO
    DO kk = 2 , ke
      khp1 = kk + 1
      DO   j = jstart, jend
        DO i = istart, iend
          ztvrml (i,j,kk)   = t(i,j,kk,nnew) *(c1 + rvd_m_o*qv(i,j,kk)         &
                                                  - qc(i,j,kk)-qrs(i,j,kk))
          zhyphl (i,j,khp1) = zhyphl(i,j,kk) +  (hhl(i,j,kk) - hhl(i,j,khp1))  &
                                               *(p0(i,j,kk) + pp(i,j,kk,nnew)) &
                                               /(rdg * ztvrml(i,j,kk))
        ENDDO
      ENDDO
    ENDDO

! top full level of temperature correction

    zptop = ptpstop - epsy
    ktop  = 1
    i = (istart + iend) / 2
    j = (jstart + jend) / 2
    vertical: DO kk = 2 , ke
      IF ( MAX( zhyphl(istart,jstart,kk) , zhyphl(iend  ,jend  ,kk)            &
              , zhyphl(istart,jend  ,kk) , zhyphl(iend  ,jstart,kk)            &
              , zhyphl(i     ,j     ,kk) )  >  zptop)              EXIT vertical
      ktop = kk
    ENDDO vertical
    lgotop = .FALSE.
    DO WHILE ((.NOT. lgotop) .AND. (ktop > 1))
      lgotop = .TRUE.
      kk = MIN( ktop , ke )
      DO   j = jstart, jend
        DO i = istart, iend
          IF (zhyphl(i,j,kk)  >  zptop)   lgotop = .FALSE.
        ENDDO
      ENDDO
      IF (.NOT. lgotop) ktop = ktop - 1
    ENDDO

! make this top level identical in all subdomains if necessary

    IF ((num_compute > 1) .AND. ((lreproduce) .OR. (luvgcor))) THEN
      IF (ltime) THEN
        CALL get_timings (i_nudging_corr, ntstep, dt, nzerr)
        CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
        CALL get_timings (i_barrier_waiting_nud, ntstep, dt, nzerr)
      ENDIF

      CALL global_values ( ktop, 1, 'MIN', imp_integers, icomm_cart, -1        &
                         , yzerrmsg, izerror )
!     ==================

      IF (ltime) CALL get_timings (i_communications_nud, ntstep, dt, nzerr)
    ENDIF

  ENDIF

!-------------------------------------------------------------------------------
!  Section 1b: Allocation of memory
!-------------------------------------------------------------------------------

ELSEIF (nactio >= 1) THEN

  IF (ntpscor > 0) THEN
    ALLOCATE ( ps_fld (istart:iend,jstart:jend     ) , STAT = nzerr )
    ALLOCATE ( ztcorr (istart:iend,jstart:jend,ke  ) , STAT = nzerr )
 
    IF (nactio == 1) THEN
      DO   j = jstart , jend
        DO i = istart , iend
          ps_fld (i,j) = psanai(i,j)
        ENDDO
      ENDDO
    ELSE  !  IF (nactio == 2)
      DO   j = jstart , jend
        DO i = istart , iend
          ps_fld (i,j) = psaigeo(i,j)
        ENDDO
      ENDDO
    ENDIF
  ENDIF
 
!-------------------------------------------------------------------------------
!  Section 1a: Preparation and header of printout for control
!-------------------------------------------------------------------------------

  lprtps = (lfst_prtps) .AND. ((lcroot) .OR. (lwonl)) .AND. (ntpscor > 0)
  IF (lprtps) THEN
    zpkemin = c2 * 109900.0_wp
    zdp1000 = c2 * 109900.0_wp
    zp1000  = c2 * 100000.0_wp * c05 *(vcoord%sigm_coord(ke+1) +    &
                                           vcoord%sigm_coord(ke))
    ipr (1) = 0
    jpr (1) = 0
    ipr (2) = 0
    jpr (2) = 0
    DO   j = jstart , jend
      DO i = istart , iend
        IF (ABS( ps_fld(i,j) ) > epsy) THEN
          IF (zhyphl(i,j,ke)+zhyphl(i,j,ke+1) < zpkemin) THEN
            ipr(2) = i
            jpr(2) = j
            zpkemin = zhyphl(i,j,ke) + zhyphl(i,j,ke+1)
          ENDIF
          IF (ABS( zhyphl(i,j,ke)+zhyphl(i,j,ke+1)-zp1000 ) < zdp1000) THEN
            ipr(1) = i
            jpr(1) = j
            zdp1000 = ABS( zhyphl(i,j,ke) + zhyphl(i,j,ke+1) - zp1000 )
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    IF (zpkemin > c2*109800.0_wp) lprtps = .FALSE.
    IF ((lprtps) .AND. (lcroot)) THEN
      kk = NINT( ptpstop / 100.0_wp )
      PRINT *
      PRINT '(''Temperature correction for (surface) pres''                    &
            &,''sure nudging: top at'',I4,''hPa in layer'',I3)' , kk, ktop
      PRINT '(''=========================================''                    &
            &,''======================================='')'
      PRINT '(''                               (*): relat''                    &
            &,''ive to the surface pressure increment  '')'
      PRINT '(''height T-corr pressure  z-z    dp level |''                    &
            &,''height T-corr pressure  z-z    dp level'')'
      PRINT '(''       K/hPa full hydro correl (*)      |''                    &
            &,''       K/hPa full hydro correl (*)     '')'
    ENDIF
    IF ((lprtps) .AND. (lwonl)) THEN
      kk = NINT( ptpstop / 100.0_wp )
      WRITE( nupr,'('' '')' )
      WRITE( nupr,'(''Temperature correction for (surface) pres''              &
                  &,''sure nudging: top at'',I4,''hPa in layer'',I3)' ) kk, ktop
      WRITE( nupr,'(''=========================================''              &
                  &,''======================================='')' )
      WRITE( nupr,'(''(surface) pressure analysis increments:''                &
                  &,F8.2,'' at'',2I4,'' ,'',F8.2,'' at'',2I4)' )               &
             ps_fld(ipr(1),jpr(1)), ipr(1), jpr(1)                             &
           , ps_fld(ipr(2),jpr(2)), ipr(2), jpr(2)
      WRITE( nupr,'(''                               (*): relat''              &
                  &,''ive to the surface pressure increment  '')' )
      WRITE( nupr,'(''height T-corr pressure  z-z    dp level |''              &
                  &,''height T-corr pressure  z-z    dp level'')' )
      WRITE( nupr,'(''       K/hPa full hydro correl (*)      |''              &
                  &,''       K/hPa full hydro correl (*)     '')' )
    ENDIF
    IF ((.NOT. lprtps) .AND. (lcroot)) THEN
      PRINT *,'No pressure increments for temperature correction'
    ELSEIF ((.NOT. lprtps) .AND. (lwonl)) THEN
      WRITE( nupr,'(''No pressure increments for temperature correction'')' )
    ENDIF
    lfst_prtps = .FALSE.
  ENDIF

!-------------------------------------------------------------------------------
!  Section 1b: Allocation of memory
!-------------------------------------------------------------------------------

  IF (ntpscor > 0) THEN
    ALLOCATE ( ztvrhl (istart:iend,jstart:jend     ) , STAT = nzerr )
    ALLOCATE ( zpdenom(istart:iend,jstart:jend     ) , STAT = nzerr )
    ALLOCATE ( zpker  (istart:iend,jstart:jend     ) , STAT = nzerr )
    ALLOCATE ( zfbuoyt(istart:iend,jstart:jend     ) , STAT = nzerr )
    ALLOCATE ( zfbuoyb(istart:iend,jstart:jend     ) , STAT = nzerr )
    ALLOCATE ( zqitps (istart:iend,jstart:jend     ) , STAT = nzerr )
    ALLOCATE ( zqitpb (istart:iend,jstart:jend     ) , STAT = nzerr )
    ALLOCATE ( zpainc (istart:iend,jstart:jend     ) , STAT = nzerr )
!   ALLOCATE ( zfisai (istart:iend,jstart:jend,2   ) , STAT = nzerr )
    ALLOCATE ( zcfi   (istart:iend,jstart:jend     ) , STAT = nzerr )
    ALLOCATE ( zdhl   (istart:iend,jstart:jend,2   ) , STAT = nzerr )
    ALLOCATE ( zcfim  (istart:iend,jstart:jend     ) , STAT = nzerr )
    ALLOCATE ( zdhlm  (istart:iend,jstart:jend     ) , STAT = nzerr )
    rdg   = r_d / g
    zptop = ptpstop - epsy
  ENDIF

!-------------------------------------------------------------------------------
!  Section 2: Computation of the temperature increments
!-------------------------------------------------------------------------------

  IF (ktop <= ke) THEN
    DO kk = 1 , ktop-1
      DO   j = jstart, jend
        DO i = istart, iend
!         ztwips (i,j,kk) = c0
          ztcorr (i,j,kk) = c0
        ENDDO
      ENDDO
    ENDDO

! setting scaling quantities, and computing the geopotential change at the
! lowest model level as approximation to the geopotential change at 'ptpstop'

    khp1 = ke + 1
    DO   j = jstart, jend
      DO i = istart, iend
        zpker   (i,j)  =  c2  / (zhyphl(i,j,ke) + zhyphl(i,j,khp1))
!       zpdenom (i,j)  =  zhyphl(i,j,khp1) - ptpstop
        zpdenom (i,j)  =  c05 *(zhyphl(i,j,ke) + zhyphl(i,j,khp1)) - ptpstop
!       zfisai  (i,j,1)=  rdg * ztvrml(i,j,ke) * zpker(i,j) * psanai(i,j)
      ENDDO
    ENDDO
  ENDIF

! initializing vertical indices

  ihpa = 1
  ihpb = 2

! compute: - hydrostatic pressure at the half level topping layer 'ktop'
!          - height or pressure change at the half level topping layer 'ktop'
!          - virtual temperature at full level 'ktop'
! ---------------------------------------------------------------------------

  IF (ktop == 1) THEN
    DO   j = jstart, jend
      DO i = istart, iend
        zdhl   (i,j,ihpb) = c0
        IF (zhyphl(i,j,1) > zptop) THEN
! height or pressure change at the top half level is non-zero;
! virtual temperature at top half level: constant lapse rate within model layer
! (i.e. the gradient at level 1/2 equals that at level 3/2, c.f. (4.13) for pp)
          ztvrhl (i,j) = ztvrml(i,j,1) + (ztvrml(i,j,1) - ztvrml(i,j,2))       &
                                        *dp0(i,j,1) /(dp0(i,j,1) + dp0(i,j,2))
! geopotential correlation at current lower half level
          zck          = MAX( c0 , (zhyphl(i,j,1) - ptpstop) / zpdenom(i,j) )
          IF (MOD(ntpscor,2) == 1) THEN
            zcfi (i,j) = zck**2 * EXP( (c1 - zck**3) * 0.125_wp )
!           zcfi (i,j) = zck**3 * EXP( (c1 - zck**4) * 0.25_wp )
          ELSEIF (MOD(ntpscor,2) == 0) THEN
            zcfi (i,j) = zck * c05 * (c1 + zck)
          ENDIF
! geopotential or pressure change at current lower half level
!         IF (ntpscor >= 3) THEN
!           zdhl (i,j,ihpb) =  zcfi(i,j) * zfisai(i,j,1)
!         ELSE
            zdhl (i,j,ihpb) =  zcfi(i,j) * ps_fld(i,j)
!           IF (MOD(ntpscor,2) == 0)                                           &
            IF (ntpscor >= 3)                                                  &
              zdhl (i,j,ihpb) =  zdhl(i,j,ihpb) * zpker(i,j) * zhyphl(i,j,1)   &
                                                * ztvrml(i,j,ke) / ztvrml(i,j,1)
!         ENDIF
        ENDIF
      ENDDO
    ENDDO
  ELSEIF (ktop <= ke) THEN
! height or pressure change at half level topping layer 'ktop' is always zero
    DO kk = 1 , ktop-1
      ihpa = 3 - ihpa
      ihpb = 3 - ihpb
      DO   j = jstart, jend
        DO i = istart, iend
          zdhl (i,j,ihpb) = c0
        ENDDO
      ENDDO
    ENDDO
  ENDIF

! compute the temperature correction for layers k <= ktop
! -------------------------------------------------------

  DO kk = ktop , ke
    khp1 =      kk + 1
    kp1  = MIN( khp1 , ke )
    ihpa = 3 - ihpa
    ihpb = 3 - ihpb
    DO   j = jstart, jend
      DO i = istart, iend
! virtual temperature at current lower half level
        ztvrhl (i,j) =  ( dp0(i,j,kp1) *ztvrml(i,j,kk )                        &
                         +dp0(i,j,kk ) *ztvrml(i,j,kp1))                       &
                      / ( dp0(i,j,kk) +dp0(i,j,kp1))
! geopotential correlation at current lower half level
        zck          = MAX( c0 , (zhyphl(i,j,khp1) - ptpstop) / zpdenom(i,j) )
        IF (MOD(ntpscor,2) == 1) THEN
          zcfi (i,j) = zck**2 * EXP( (c1 - zck**3) * 0.125_wp )
!         zcfi (i,j) = zck**3 * EXP( (c1 - zck**4) * 0.25_wp )
        ELSEIF (MOD(ntpscor,2) == 0) THEN
          zcfi (i,j) = zck * c05 * (c1 + zck)
        ENDIF
        zck          = MAX( c0 , (  c05 *(zhyphl(i,j,kk) + zhyphl(i,j,khp1))   &
                                  - ptpstop) / zpdenom(i,j) )
        IF (MOD(ntpscor,2) == 1) THEN
          zcfim (i,j) = zck**2 * EXP( (c1 - zck**3) * 0.125_wp )
!         zcfim (i,j) = zck**3 * EXP( (c1 - zck**4) * 0.25_wp )
        ELSEIF (MOD(ntpscor,2) == 0) THEN
          zcfim (i,j) = zck * c05 * (c1 + zck)
        ENDIF
      ENDDO
    ENDDO
    DO   j = jstart, jend
      DO i = istart, iend
! geopotential or pressure change at current lower half level and
! temperature correction at current full level
!       IF (ntpscor >= 3) THEN
!         zdhl (i,j,ihpb) =  zcfi(i,j) * zfisai(i,j,1)
!         ztwips (i,j,kk) =  (zdhl  (i,j,ihpa) - zdhl  (i,j,ihpb))             &
!                          * (zhyphl(i,j,khp1) + zhyphl(i,j,kk)) * c05         &
!                          /((zhyphl(i,j,khp1) - zhyphl(i,j,kk)) * rdg)
!       ELSE
          zdhl (i,j,ihpb) =  zcfi(i,j) * ps_fld(i,j)
          zdhlm  (i,j)    =  zcfim(i,j) * ps_fld(i,j)
!         IF (MOD(ntpscor,2) == 0) THEN
          IF (ntpscor >= 3) THEN
            zdhl (i,j,ihpb) =  zdhl(i,j,ihpb) * zpker(i,j) * zhyphl(i,j,khp1)  &
                             * ztvrml(i,j,ke) / ztvrhl(i,j)
            zdhlm  (i,j)    =  zdhlm(i,j)     * zpker(i,j)                     &
                             * c05 * (zhyphl(i,j,kk) + zhyphl(i,j,khp1))       &
                             * ztvrml(i,j,ke) / ztvrml(i,j,kk)
          ENDIF
!         ztwips (i,j,kk) =  ( (hhl(i,j,kk) - hhl(i,j,khp1)) /rdg *zdhlm(i,j)  &
          ztcorr (i,j,kk) =  ( (hhl(i,j,kk) - hhl(i,j,khp1)) /rdg *zdhlm(i,j)  &
                              -(zdhl  (i,j,ihpb) - zdhl  (i,j,ihpa))           &
                               *ztvrml(i,j,kk))                                &
                            /(  zhyphl(i,j,khp1) - zhyphl(i,j,kk)              &
                              + zdhl  (i,j,ihpb) - zdhl  (i,j,ihpa))
!       ENDIF
!       ztwips (i,j,kk) = ztwips(i,j,kk) * t(i,j,kk,nnew) / ztvrml(i,j,kk)
        ztcorr (i,j,kk) = ztcorr(i,j,kk) * t(i,j,kk,nnew) / ztvrml(i,j,kk)
      ENDDO
    ENDDO
  ENDDO


!-------------------------------------------------------------------------------
!  Section 2b: 2nd iteration, if increments should vanish above a certain height
!-------------------------------------------------------------------------------

  IF ((ktop > 1) .AND. (ktop <= ke)) THEN
    DO   j = jstart, jend
      DO i = istart, iend
        zfbuoyt (i,j)  =  r_d * rho0(i,j,ke) * t(i,j,ke,nnow)
!       zpainc  (i,j)  =  psanai(i,j)
        zpainc  (i,j)  =  ps_fld(i,j)
      ENDDO
    ENDDO
  ENDIF
  IF ((ktop > 1) .AND. (ktop <= ke) .AND. (khumbal < 100)) THEN
    DO   j = jstart, jend
      DO i = istart, iend
        zrhmod       =  fq2pv( qv(i,j,ke) , p0(i,j,ke) +pp(i,j,ke,nnew) , rdv) &
                      / fpvsw( t (i,j,ke,nnew) , b1, b2w, b3, b4w )
! specific humidity increment due to the temperature correction
!       zqitps (i,j) =  fpv2q( fpvsw( t(i,j,ke,nnew) + ztwips(i,j,ke)          &
!                                   , b1, b2w, b3, b4w ) * zrhmod              &
        zqitps (i,j) =  fpv2q( fpvsw( t(i,j,ke,nnew) + ztcorr(i,j,ke)          &
                                    , b1, b2w, b3, b4w ) * zrhmod              &
                             , p0(i,j,ke) +pp(i,j,ke,nnew) , rdv )  - qv(i,j,ke)
        zqitps (i,j) =  zqitps(i,j) * rvd_m_o
      ENDDO
    ENDDO
    DO kk = ke-1 , ktop , -1
      kp1  = kk + 1
      khp1 = kk + 1
      DO   j = jstart, jend
        DO i = istart, iend
          zfbuoyb(i,j) =  zfbuoyt(i,j)
          zfbuoyt(i,j) =  r_d * rho0(i,j,kk) * t(i,j,kk,nnow)
          zqitpb (i,j) =  zqitps (i,j)
          zrhmod       = fq2pv( qv(i,j,kk) , p0(i,j,kk)+pp(i,j,kk,nnew) , rdv) &
                        /fpvsw( t (i,j,kk,nnew) , b1, b2w, b3, b4w )
          zqitps (i,j) = fpv2q( fpvsw( t(i,j,kk,nnew) + ztcorr(i,j,kk)         &
                                     , b1, b2w, b3, b4w ) * zrhmod             &
                              , p0(i,j,kk) +pp(i,j,kk,nnew) , rdv ) - qv(i,j,kk)
          zqitps (i,j) =  zqitps(i,j) * rvd_m_o
!         IF (zhyphl(i,j,khp1) > zptop)                                        &
            zpainc (i,j) =         c1     / (c2 + dp0(i,j,kp1) /zfbuoyt(i,j))  &
                          *(  zpainc(i,j) * (c2 - dp0(i,j,kk ) /zfbuoyb(i,j))  &
                            + ztcorr(i,j,kk ) * dp0(i,j,kp1) / t(i,j,kk ,nnow) &
                            + ztcorr(i,j,kp1) * dp0(i,j,kk ) / t(i,j,kp1,nnow) &
                            + zqitps(i,j)     * dp0(i,j,kp1)                   &
                            + zqitpb(i,j)     * dp0(i,j,kk )  )
! (Note: - (zhyphl(i,j,kk) > zptop) only if (kk == ktop == 1) .
!        - the denominator is equivalent to integrating 'zpainc' hydrostatically
!          from level 'ktop' ('kk') down to level 'ke' (without T,q-increments))
!         IF (       (zhyphl(i,j,khp1) >  zptop)                               &
!             .AND. ((zhyphl(i,j,kk  ) <= zptop) .OR. (kk == ktop))            &
!             .AND.  (ABS(ps_fld(i,j)) >  epsy))                               &
!           zpainc (i,j)   = c1 +  zpainc(i,j) / ps_fld(i,j)                   &
!                                /( c05 *(zhyphl(i,j,kk) + zhyphl(i,j,khp1))   &
!                                  *zpker (i,j))
!           zfisai (i,j,2) =  zpainc(i,j) * rdg * ztvrml(i,j,kk)               &
!                                         / (p0(i,j,kk) + pp(i,j,kk,nnew))
        ENDDO
      ENDDO
    ENDDO
  ELSEIF ((ktop > 1) .AND. (ktop <= ke)) THEN
    DO kk = ke-1 , ktop , -1
      kp1  = kk + 1
      khp1 = kk + 1
      DO   j = jstart, jend
        DO i = istart, iend
          zfbuoyb(i,j) =  zfbuoyt(i,j)
          zfbuoyt(i,j) =  r_d * rho0(i,j,kk) * t(i,j,kk,nnow)
!         IF (zhyphl(i,j,khp1) > zptop)                                        &
            zpainc (i,j) =         c1     / (c2 + dp0(i,j,kp1) /zfbuoyt(i,j))  &
                          *(  zpainc(i,j) * (c2 - dp0(i,j,kk ) /zfbuoyb(i,j))  &
                            + ztcorr(i,j,kk ) * dp0(i,j,kp1) / t(i,j,kk ,nnow) &
                            + ztcorr(i,j,kp1) * dp0(i,j,kk ) / t(i,j,kp1,nnow) )
! (Note: - (zhyphl(i,j,kk) > zptop) only if (kk == ktop == 1) .
!        - the denominator is equivalent to integrating 'zpainc' hydrostatically
!          from level 'ktop' ('kk') down to level 'ke' (without T,q-increments))
!         IF (       (zhyphl(i,j,khp1) >  zptop)                               &
!             .AND. ((zhyphl(i,j,kk  ) <= zptop) .OR. (kk == ktop))            &
!             .AND.  (ABS(ps_fld(i,j)) >  epsy))                               &
!           zpainc (i,j)   = c1 +  zpainc(i,j) / ps_fld(i,j)                   &
!                                /( c05 *(zhyphl(i,j,kk) + zhyphl(i,j,khp1))   &
!                                  *zpker (i,j))
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  IF ((ktop > 1) .AND. (ktop <= ke)) THEN
    DO   j = jstart, jend
      DO i = istart, iend
!   the correction to the denominator is equivalent to integrating 'zpainc'
!   hydrostatically from level 'ktop' (kk) down to level ke (without T,q-incr.))
        IF (ABS(ps_fld(i,j)) > epsy) THEN
          zpainc (i,j) =  ps_fld(i,j)                                          &
                        /(ps_fld(i,j) -  zpainc(i,j)                           &
                                       /(c05 *( zhyphl(i,j,ktop)               &
                                               +zhyphl(i,j,ktop+1))*zpker(i,j)))
        ENDIF
      ENDDO
    ENDDO
    DO kk = ktop , ke
      DO   j = jstart, jend
        DO i = istart, iend
!         IF (ABS(zfisai(i,j,1)) > epsy)                                       &
!           ztcorr (i,j,kk) =  ztcorr(i,j,kk)                                  &
!                            * (c1 + zfisai(i,j,2) / zfisai(i,j,1))
          IF (ABS(ps_fld(i,j)) > epsy)                                         &
            ztcorr (i,j,kk) =  ztcorr(i,j,kk)  *  zpainc(i,j)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

! storing of temperature increments on the appropriate arrays

  IF (ntpscor > 0) THEN
    IF (nactio == 1) THEN
      DO kk = 1 , ke
        DO   j = jstart , jend
          DO i = istart , iend
            ztwips (i,j,kk) = ztcorr(i,j,kk)
          ENDDO
        ENDDO
      ENDDO
    ELSE  !  IF (nactio == 2)
      DO kk = 1 , ke
        DO   j = jstart , jend
          DO i = istart , iend
            ztpgeo (i,j,kk) = ztcorr(i,j,kk)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
!   zmaxai = MAX(   MAXVAL(  ztcorr(:,:,ke) ) ,   MAXVAL(  ztcorr(:,:,ke-1) ) )
!   zminai = MIN( - MAXVAL( -ztcorr(:,:,ke) ) , - MAXVAL( -ztcorr(:,:,ke-1) ) )
!   IF ((lwonl) .OR. (MAX( zmaxai,ABS(zminai) ) > c1))                         &
!     PRINT *, 'max(ztcorr):', nactio, zmaxai, zminai
  ENDIF

!-------------------------------------------------------------------------------
!  Section 2c: Printout for control
!-------------------------------------------------------------------------------

  IF (lprtps) THEN
    DO ilpr = 1 , 2
      IF (ipr(ilpr) > 0) THEN
        i = ipr(ilpr)
        j = jpr(ilpr)
        zpaipr (ke,ilpr) =  ps_fld(i,j)
        zfbuoyt(i,j) =  r_d * rho0(i,j,ke) * t(i,j,ke,nnow)
      ENDIF
      IF ((ipr(ilpr) > 0) .AND. (khumbal < 100)) THEN
        zrhmod       =  fq2pv( qv(i,j,ke) , p0(i,j,ke) +pp(i,j,ke,nnew) , rdv) &
                      / fpvsw( t (i,j,ke,nnew) , b1, b2w, b3, b4w )
! specific humidity increment due to the temperature correction
        zqitps (i,j) =  fpv2q( fpvsw( t(i,j,ke,nnew) + ztcorr(i,j,ke)          &
                                    , b1, b2w, b3, b4w ) * zrhmod              &
                             , p0(i,j,ke) +pp(i,j,ke,nnew) , rdv )  - qv(i,j,ke)
        zqitps (i,j) =  zqitps(i,j) * rvd_m_o
        DO kk = ke-1 , ktop , -1
          kp1  = kk + 1
          zfbuoyb(i,j) =  zfbuoyt(i,j)
          zfbuoyt(i,j) =  r_d * rho0(i,j,kk) * t(i,j,kk,nnow)
          zqitpb (i,j) =  zqitps (i,j)
          zrhmod       = fq2pv( qv(i,j,kk) , p0(i,j,kk)+pp(i,j,kk,nnew) , rdv) &
                        /fpvsw( t (i,j,kk,nnew) , b1, b2w, b3, b4w )
          zqitps (i,j) = fpv2q( fpvsw( t(i,j,kk,nnew) + ztcorr(i,j,kk)         &
                                     , b1, b2w, b3, b4w ) * zrhmod             &
                              , p0(i,j,kk) +pp(i,j,kk,nnew) , rdv ) - qv(i,j,kk)
          zqitps (i,j) =  zqitps(i,j) * rvd_m_o
!         IF (zhyphl(i,j,khp1) > zptop)                                        &
            zpaipr (kk,ilpr) =     c1     / (c2 + dp0(i,j,kp1) /zfbuoyt(i,j))  &
                     *(  zpaipr(kp1,ilpr) * (c2 - dp0(i,j,kk ) /zfbuoyb(i,j))  &
                            + ztcorr(i,j,kk ) * dp0(i,j,kp1) / t(i,j,kk ,nnow) &
                            + ztcorr(i,j,kp1) * dp0(i,j,kk ) / t(i,j,kp1,nnow) &
                            + zqitps(i,j)     * dp0(i,j,kp1)                   &
                            + zqitpb(i,j)     * dp0(i,j,kk )  )
        ENDDO
      ELSEIF (ipr(ilpr) > 0) THEN
        DO kk = ke-1 , ktop , -1
          kp1  = kk + 1
          zfbuoyb(i,j) =  zfbuoyt(i,j)
          zfbuoyt(i,j) =  r_d * rho0(i,j,kk) * t(i,j,kk,nnow)
!         IF (zhyphl(i,j,khp1) > zptop)                                        &
            zpaipr (kk,ilpr) =     c1     / (c2 + dp0(i,j,kp1) /zfbuoyt(i,j))  &
                     *(  zpaipr(kp1,ilpr) * (c2 - dp0(i,j,kk ) /zfbuoyb(i,j))  &
                            + ztcorr(i,j,kk ) * dp0(i,j,kp1) / t(i,j,kk ,nnow) &
                            + ztcorr(i,j,kp1) * dp0(i,j,kk ) / t(i,j,kp1,nnow) )
        ENDDO
      ENDIF
    ENDDO
    DO kk = ktop , ke
      khp1 = kk + 1
      DO ilpr = 1 , 2
        i = ipr(ilpr)
        j = jpr(ilpr)
        zprh  = c05 * (hhl(i,j,kk) + hhl(i,j,khp1))
        zprhp = (zhyphl(i,j,kk) + zhyphl(i,j,khp1)) / 200.0_wp
        zprp  = (p0(i,j,kk)     + pp(i,j,kk,nnew) ) / 100.0_wp
        zck   = MAX( c0 , (zprhp*100.0_wp - ptpstop) / zpdenom(i,j) )
        IF (MOD(ntpscor,2) == 1) THEN
          zprcz = zck**2 * EXP( (c1 - zck**3) * 0.125_wp )
!!        zprcz = zck**3 * EXP( (c1 - zck**4) * 0.25_wp )
          zprcz = zprcz /( zpker(i,j) *c05 *(zhyphl(i,j,kk) +zhyphl(i,j,khp1)) &
                          *ztvrml(i,j,ke) / ztvrml(i,j,kk))
        ELSEIF (MOD(ntpscor,2) == 0) THEN
          zprcz = zck * c05 * (c1 + zck)
        ENDIF
!!      zprdz = zprcz *(zfisai(i,j,1)+zfisai(i,j,2)) / ps_fld(i,j) *100.0_wp
!       zprdz = zprcz
!       IF (MOD(ntpscor,2) == 0)                                               &
!         zprdz = zprcz * zpker(i,j) *c05 *(zhyphl(i,j,kk) + zhyphl(i,j,khp1)) &
!                       * ztvrml(i,j,ke) / ztvrml(i,j,kk)
        zprdz = zpaipr(kk,ilpr) / ps_fld(i,j)
        zprdt = ztwips(i,j,kk) / ps_fld(i,j) *100.0_wp
        WRITE( ypr(ilpr),'(F6.0,F7.3,2F5.0,F6.3,F7.3,I3)' )                    &
               zprh, zprdt, zprp, zprhp, zprcz, zprdz, kk
      ENDDO
      IF (lcroot) PRINT       '(A39,'' |'',A39)' , ypr(1), ypr(2)
      IF (lwonl)  WRITE( nupr,'(A39,'' |'',A39)' ) ypr(1), ypr(2)
      IF ((lwonl) .AND. (kk == ke)) WRITE( nupr,'('' '')' )
    ENDDO
  ENDIF

! clean-up

  IF (ntpscor > 0) THEN
    IF (((nactio == 1) .AND. (mpsgcor == 0)) .OR. (nactio == 2)) THEN
      DEALLOCATE ( ztvrml , STAT = nzerr )
      DEALLOCATE ( zhyphl , STAT = nzerr )
    ENDIF

    DEALLOCATE ( ps_fld , STAT = nzerr )
    DEALLOCATE ( ztcorr , STAT = nzerr )

    DEALLOCATE ( ztvrhl , STAT = nzerr )
    DEALLOCATE ( zpdenom, STAT = nzerr )
    DEALLOCATE ( zpker  , STAT = nzerr )
    DEALLOCATE ( zfbuoyt, STAT = nzerr )
    DEALLOCATE ( zfbuoyb, STAT = nzerr )
    DEALLOCATE ( zpainc , STAT = nzerr )
    DEALLOCATE ( zqitps , STAT = nzerr )
    DEALLOCATE ( zqitpb , STAT = nzerr )
!   DEALLOCATE ( zfisai , STAT = nzerr )
    DEALLOCATE ( zcfi   , STAT = nzerr )
    DEALLOCATE ( zdhl   , STAT = nzerr )
    DEALLOCATE ( zcfim  , STAT = nzerr )
    DEALLOCATE ( zdhlm  , STAT = nzerr )
  ENDIF

  ktoptps = ktop
  IF (ktop == ke+1) ktoptps = 1

  ! flush YUPRINT file
  IF ((ldump_ascii) .AND. (lwonl) .AND. (lfirst)) THEN
    CLOSE (nupr)
    OPEN  (nupr,FILE=yuprint,FORM='FORMATTED',STATUS='OLD'                     &
                            ,POSITION='APPEND',IOSTAT=nstat)
    IF (nstat /= 0) PRINT '("OPENING OF FILE yuprint FAIL, ps_temperatur_corr")'
  ENDIF

ENDIF


!-------------------------------------------------------------------------------
! End of module procedure ps_temperatur_corr
!-------------------------------------------------------------------------------

END SUBROUTINE ps_temperatur_corr


!-------------------------------------------------------------------------------
!+ Module procedure for geostrophic pressure correction from surfacee wind incr.
!-------------------------------------------------------------------------------

SUBROUTINE geostroph_ps_corr

!-------------------------------------------------------------------------------
!
! Description:
!   This procedure produces (near-)surface pressure increments which are in
!   geostrophic balance with (the rotational part of) the (near-)surface
!   analysis increment field derived from (selected types of) surface-level
!   wind observations (e.g. scatterometer wind data).
!
! Method:
!   1.  Computation of wind analysis increments at the center of the Arakawa-C
!       grid points at the lowest model level from the quantities (sums of
!       weights and weighted observation increments) prepared in the spreading
!       subroutine 'surf_upair_spread'. These increments are derived from
!       selected observation types which may have different specified weighting.
!       The lateral spreading for each of the observation increments is
!       isotropic, i.e. orographic effects are excluded (to obtain smooth wind
!       increments and decrease ageostrophic parts). In this sense, the
!       increments are regarded as being representative for the mean sea level.
!   2.  One boundary line is exchanged of the wind increment field, and
!       the rotation of the field ( rot(v) ) is computed, which is one part
!       of the input for the Poisson solver.
!   3.  To pre-condition the Poisson solver, a smooth first guess geopotential
!       increment field is derived. The applied algorithms are not readily
!       parallelised and therefore solved at one node for the whole domain on
!       a coarse grid:
!     3a. The input wind increments fields are collected by node 0.
!     3b. The collected fields are (bi-linearly) interpolated to a coarse grid.
!     3c. To obtain the balanced geopotential increment field, the Poisson
!         equation is solved iteratively such that the direction of the
!         integration is alternated from iteration to iteration.
!     3d. The coarse geopotential increment field is (bi-linearly) interpolated
!         to the (fine) model grid.
!     3e. The global geopotential field is distributed to the local sub-domains.
!   4.  The parallelised Poisson solver interatively derives fine-scale
!       geopotential increment fields using ( rot(v) ) of the original wind
!       field and the smooth geopotential first guess as input.
!   5.  The geopotential increments representative for the mean sea level are
!       converted to pressure increments valid at the lowest model level
!       and scaled according to parameter 'qgeops'.
!     5a. Conversion to pressure increments at mean sea level, by extrapolating
!         model temperature and pressure to mean sea level.
!     5b. Extrapolation of the pressure increments from mean sea level to the
!         lowest model level by applying the correction factor used in the
!         lateral spreading of pressure observation increments according to
!         the temperature correction.
!     5c. Reduction of pressure increments according to parameter 'qgeops'.
!
! Written by        :  Christoph Schraff + Heinz-Werner Bitzer, DWD
!                                                   (original version: 20.04.07)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
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

! Local parameters:
! ----------------

! REAL    (KIND=wp   )     , PARAMETER ::  &
!   plbcno = 10.0_wp ,& ! Number of grid rows at the lateral boundaries of
!                           ! the total domain with zero weight for the geo-
!                           ! strophic wind increments
!   plbcipr= 0.125_wp   ! Number of grid rows within which the horizontal
!                           ! weight to these increments increase linearly from
!                           ! zero at the lateral boundaries to one in the inner
!                           ! part of the total domain  (reciproke of idem)

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    i, j, kzdims(24) ,& ! loop indices in horizontal direction
    jm1              ,& ! index MAX( j-1, 1 )
    ityp             ,& ! index of set of obs system in weighted increment array
    irm              ,& ! error status variable
                        ! --> variables related to Poisson solvers:
    ntpois, ntpois_c ,& ! number of iterations used in the Poisson solvers
    itpois           ,& ! loop over iterations used in the Poisson solvers
    ige, jge         ,& ! number of grid points in the coarse grid (each direc.)
    ig , jg , jgp1   ,& ! loop indices on the coarse intermediate grid
    in1, jn1, jnd    ,& ! indices on the coarse intermediate grid
    in2, jn2         ,& ! indices on the coarse intermediate grid
    np1, np2         ,& ! indices for direction in Poisson solver (coarse grid)
    it0, it1, it2    ,& ! cycling indices in the Poisson solver (fine grid)
    ierror  , istat(9)  ! error status variables

  REAL    (KIND=wp   )     ::  &
    fmultw           ,& ! 0. or 1., depending on multiple observation weighting
    zdt              ,& ! timestep relevant for the nudging
    zgnudg           ,& ! basic nudging coefficient for horizontal wind
    zlapse           ,& ! constant lapse rate to extrapolate to mean sea level
    r_d_r            ,& ! 1 / r_d
    zaimx            ,& ! maximum absolute near-surface wind analysis increment
    zmaxai , zminai  ,& !
    zmaxai_1,zminai_1,& !
    zdfi_max,zdfi_min,& !
    zqmass           ,& ! mass affected by 'temp. corr.' devided by 'zmassps'
    zqmass2, zqmass3 ,& ! zqmass  **2    ;   zqmass **3
                        ! --> variables related to Poisson solvers:
    zfpois           ,& ! tuning parameter in the Poisson solver
    dlatg   , dlong  ,& ! mesh width (in degrees) of coarse grid 
    odadlatg         ,& ! 1 / (meridional mesh width (in km) of coarse grid)
    dinkx   , dinky  ,& ! ratio of mesh width between original and target grid
    e1 , e2 , e3 , e4,& ! geometrical factors used in the Poisson solver
    x3      , y3     ,& ! position of target point in reference grid
    xn      , yn     ,& ! position of target point within grid cell
                        ! diagnostics:
    umax             ,& ! max. geostrophic wind from geopotential increm. field
    divmax  , dvmax  ,& ! max. differences of geostr. wind between iterations
    dv               ,& ! differences of geostr. wind between iterations
    q11              ,& ! relative max. diff. of geostr. wind between iterations
    sum     , xnum   ,& ! other variables for diagnostics
    epsv                ! threshold on 'q11' to stop iteration (not implemented)
    
  LOGICAL                  ::  &
    lwrite              ! printout for control

  CHARACTER (LEN=20)       ::  &
    yroutine            ! name of this subroutine
  CHARACTER (LEN=40)       ::  &
    yerrmsg             ! error message
  CHARACTER (LEN=50)       :: yerror

! Local (automatic) arrays:
! -------------------------

  REAL    (KIND=wp   )     ::  &
    zwsum   (ie,je,2),& ! sum of weights assigned to the different obs. types
    zaiu    (ie,je)  ,& ! analysis increment of zonal  wind on local sub-domain
    zaiv    (ie,je)  ,& ! analysis increment of merid. wind on local sub-domain
    zumy    (ie,je)  ,& ! 2*timestep times net nudging weight shifted to u-g.p.
    zrotv   (ie,je)  ,& ! rotation of the wind increment field ( rot(dv) )
    zfcc    (ie,je)  ,& ! Coriolis parameter at mass points of C-grid
    zdfi    (ie,je,2),& ! geopotential increment field
    ztmsl   (ie,je)  ,& ! near-surface temperature reduced to mean sea level
    zpmsl   (ie,je)  ,& ! near-surface pressure    reduced to mean sea level
    zsprcor (ie,je)  ,& ! correction factor depending on orographic height
    zpke_v  (ie,je)  ,& ! model pressure at lowest model level
    cfitop_v(ie,je)  ,& ! weighted conveyor of geopot. correl. to p- correlation
    zdxr       (je)  ,& ! reciprocal x-grid spacing (1/( a cos(phi) dlam )
    zdyr       (je)     ! 1/( a cos(phi) dlat )

!! REAL    (KIND=wp   )     ::  &
!!  avai_ff (   ke)  ,& ! horiz. average of wind speed     analysis increments
!!  avai_dd (   ke)  ,& ! horiz. average of wind direction analysis increments
!!  zugeai  (   ke)  ,& ! zonal  comp. of geostrophic wind incr. (for printout)
!!  zvgeai  (   ke)  ,& ! merid. comp. of geostrophic wind incr. (for printout)

  REAL    (KIND=wp   )     , ALLOCATABLE ::       &
    zaiu_tot  (:,:)  ,& ! analysis increment of zonal  wind on global domain
    zaiv_tot  (:,:)  ,& ! analysis increment of merid. wind on global domain
    zdfi_tot  (:,:)     ! global first guess geopotential increment field
                        !   (divided by Coriolis parameter)

  REAL    (KIND=wp   )     , ALLOCATABLE ::       &
                         ! following variables are defined on the coarse grid:
    zdxg(:) , zdyg(:) ,& ! 1 / (mesh width (in km) of coarse grid * COS( lat ))
    ug(:,:) , vg(:,:) ,& ! wind vector interpolated to the coarse grid
    zrotv_c     (:,:) ,& ! rotation of the wind increment field ( rot(dv) )
    fi_0        (:,:) ,& ! geopotential increment fields derived from 'zrotv_c'
                         !   (divided by Coriolis parameter)
    ug2(:,:), vg2(:,:),& ! wind vector derived from geopotential increment field
    crlatg(:)         ,& ! COS( latitude )
                         ! following variables are defined on the model grid:
    upn(:,:), vpn(:,:),& ! wind vector from geopot. increm. field, new iteration
    upo(:,:), vpo(:,:)   ! wind vector from geopot. increm. field, old iteration
!
!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine geostroph_ps_corr
!-------------------------------------------------------------------------------

  kzdims(:) = 0_iintegers

  yroutine = 'geostroph_ps_corr'
  irm = 0

  DO   j = 1, je
    DO i = 1, ie
      zaiu  (i,j)    =  c0
      zaiv  (i,j)    =  c0
      zrotv (i,j)    =  c0
    ENDDO
  ENDDO

!-------------------------------------------------------------------------------
!  Section 1 : Computation of wind analysis increments at the center of the
!              Arakawa-C grid points at the lowest model level.
!              These increments are derived from selected observation types
!              with different specified weighting.
!-------------------------------------------------------------------------------
!  Section 1a: Computation of a
!              1):    squared weighted mean   =    net (spreaded) obs. increment
!                     of obs. increments           at the target grid pt.
!              2):    weighted mean           =    net nudging weight
!                     of nudging weights           at the target grid pt.
!        Note: The calculations in this section are as in section 2a of
!              subr. nudge_horiz_wind.
!-------------------------------------------------------------------------------

  IF (l2tls) THEN
    zdt = dt
  ELSE
    zdt = dt2
  ENDIF

  lwrite  =  (lwonl) .AND. (ntstep <= 5)
  IF (lwrite) WRITE( nupr,'("wup: weight/ sqr/ weighted incr",4F9.4)' )        &
                     omy(ionl,jonl,5,1), om2(ionl,jonl,4,1)                    &
                   , zwi(ionl,jonl,5,1), zwi(ionl,jonl,6,1)
! Inclusion of nudging coefficient (the MAX(.,epsy) is required for nudging of
! surface-level data, if the nudging coeff. for TEMP / PILOT is set to zero)
  zgnudg = MAX( gnudg(1) , epsy )

  DO   j = jstart, jend
    DO i = istart, iend
      zwsum (i,j,1)  =  c0
      zwsum (i,j,2)  =  c0
    ENDDO
  ENDDO
  DO ityp = 1 , nwtyp
    fmultw = REAL( kwtyp(ityp) - 1, wp )
    DO   j = jstart, jend
      DO i = istart, iend
        IF ((om2(i,j,4,ityp) > epsy) .AND. (omy(i,j,5,ityp) > epsy)) THEN
! averaged increment for observation type 'ityp'
          zwi (i,j,5,ityp)  =   zwi(i,j,5,ityp)                                &
                             / (om2(i,j,4,ityp) + fmultw *omy(i,j,5,ityp))
          zwi (i,j,6,ityp)  =   zwi(i,j,6,ityp)                                &
                             / (om2(i,j,4,ityp) + fmultw *omy(i,j,5,ityp))
! individual weight 'w(ityp)' assigned to the observation type 'ityp'
          omy (i,j,4,ityp)  =  (om2(i,j,4,ityp) + fmultw *omy(i,j,5,ityp))     &
                             / (omy(i,j,5,ityp) + fmultw)
! sum of the individual weights 'w(ityp)' assigned to the different obs. types
          zwsum (i,j,1)     =  omy(i,j,4,ityp)  +  zwsum(i,j,1)
        ELSE
          zwi (i,j,5,ityp)  =  c0
          zwi (i,j,6,ityp)  =  c0
          omy (i,j,4,ityp)  =  c0
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  fmultw = REAL( kwtyp(1) - 1, wp )
  IF (nwtyp >= 2)  fmultw = REAL( kwtyp(nwtyp+1) - 1, wp )
  DO ityp = 1 , nwtyp
    DO   j = jstart, jend
      DO i = istart, iend
        IF (zwsum(i,j,1) > epsy) THEN
! total multiple weight 'W(ityp)' assigned to observation type 'ityp'
          omy (i,j,5,ityp)  =  omy(i,j,4,ityp) * (omy(i,j,4,ityp) + fmultw)    &
                                               / (zwsum(i,j,1)    + fmultw)
! weighted averaged increment for observation type 'ityp'
          zwi (i,j,5,ityp)  =  zwi(i,j,5,ityp) * omy(i,j,5,ityp)
          zwi (i,j,6,ityp)  =  zwi(i,j,6,ityp) * omy(i,j,5,ityp)
! sum of total multiple weights 'W(ityp)' assigned to the different obs. types
          zwsum (i,j,2)     =  omy(i,j,5,ityp)  +  zwsum(i,j,2)
        ELSE
          omy (i,j,5,ityp)  =  c0
          zwi (i,j,5,ityp)  =  c0
          zwi (i,j,6,ityp)  =  c0
        ENDIF
! sum of the weighted averaged increments for the different observation types
        zaiu (i,j)  =  zwi(i,j,5,ityp)  +  zaiu(i,j)
        zaiv (i,j)  =  zwi(i,j,6,ityp)  +  zaiv(i,j)
      ENDDO
    ENDDO
  ENDDO
  DO   j = jstart, jend
    DO i = istart, iend
      IF (zwsum(i,j,2) > epsy) THEN
! double-averaged total observation increment
        zaiu (i,j)  =  zaiu(i,j) / zwsum(i,j,2)
        zaiv (i,j)  =  zaiv(i,j) / zwsum(i,j,2)
! final weight, including nudging coefficient
        zumy (i,j)  =  zwsum(i,j,2)  * zgnudg
      ELSE
        zaiu (i,j)  =  c0
        zaiv (i,j)  =  c0
        zumy (i,j)  =  c0
      ENDIF
    ENDDO
  ENDDO
  IF (lwrite) WRITE( nupr,'("wup: weight    / weighted incr",F12.7,6X,2F9.4)') &
                     zumy(ionl,jonl), zaiu(ionl,jonl), zaiv(ionl,jonl)

!-------------------------------------------------------------------------------
!  Section 1b: Computation of the final wind analysis increments
!-------------------------------------------------------------------------------

  DO   j = jstart, jend
    DO i = istart, iend
      zumy (i,j)  =  zumy(i,j) * zdt * faclbc(i,j,1)
      zaiu (i,j)  =  zumy(i,j) / (c1 + zumy(i,j)) * zaiu(i,j)
      zaiv (i,j)  =  zumy(i,j) / (c1 + zumy(i,j)) * zaiv(i,j)
    ENDDO
  ENDDO

  IF ((lwonl2) .AND. (ntstep*dtdeh <= c2))                                     &
    WRITE( nupr,'(''geostroph_ps_corr: uv-ana.incr. '',2I4, 2F9.4)' )          &
           ionl2, jonl2, zaiu(ionl2,jonl2), zaiv(ionl2,jonl2)

!-------------------------------------------------------------------------------
!  Section 1c: Check whether wind analysis increments are non-zero
!              and the balancing has to be computed or not
!-------------------------------------------------------------------------------

  zaimx  =  MAXVAL( ABS( zaiu ) )  +  MAXVAL( ABS( zaiv ) )

  ! exchange of 2 boundary lines of the wind increments (!)
  IF (num_compute > 1) THEN
    IF (ltime) THEN
      CALL get_timings (i_nudging_corr, ntstep, dt, izerror)
      CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
      CALL get_timings (i_barrier_waiting_nud, ntstep, dt, izerror)
    ENDIF

    CALL global_values ( zaimx, 1,'MAX', imp_reals,icomm_cart, -1,yerror,ierror)
!   ==================

    IF (ltime) CALL get_timings (i_communications_nud, ntstep, dt, izerror)
  ENDIF

  IF (zaimx < epsy) THEN
    psaigeo (:,:) = c0
                                                                          RETURN
  ENDIF

!-------------------------------------------------------------------------------
! Section 2: Computation of the rotation of the wind analysis increment field
!            ( rot(v) )    (multiplied by the Coriolis parameter)
!-------------------------------------------------------------------------------

  ! Coriolis parameter at the mass grid points
  DO j = 1, je
    jm1 = MAX( j-1, 1 )
    zfcc (1,j) = c05 *(fc(1,j) + fc(1,jm1))
    DO i = 2, ie
      zfcc (i,j) = 0.25_wp *(fc(i,j) + fc(i-1,j) + fc(i,jm1) + fc(i-1,jm1))
    ENDDO
  ENDDO

  ! exchange of 2 boundary lines of the wind increments (!)
  IF (num_compute > 1) THEN
    IF (ltime) THEN
      CALL get_timings (i_nudging_corr, ntstep, dt, izerror)
      CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
      CALL get_timings (i_barrier_waiting_nud, ntstep, dt, izerror)
    ENDIF

    kzdims(1:24)=(/1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)

    CALL exchg_boundaries ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart,       &
                           num_compute, ie, je, kzdims, jstartpar, jendpar, 2,    &
                           nboundlines, my_cart_neigh, lperi_x, lperi_y, l2dim,   &
                           20000+nexch_tag,.FALSE., ncomm_type, izerror, yzerrmsg,&
                           zaiu(:,:), zaiv(:,:) )
!   =====================

    IF (ltime) CALL get_timings (i_communications_nud, ntstep, dt, izerror)
  ENDIF

  ! metrical factors for horizontal discretization
  DO  j = 1 , je
    zdxr  (j) =     acrlat(j,1)*eddlon
    zdyr  (j) =     acrlat(j,1)*eddlat
  ENDDO

  ! rotation (vorticity) of the wind increment field
  DO   j = MAX( jstart-1 , 2 ) , MIN( jend+1 , je-1 )
    DO i = MAX( istart-1 , 2 ) , MIN( iend+1 , ie-1 )
      zrotv(i,j) =   zdxr(j) *c05 *(zaiv(i+1,j) - zaiv(i-1,j))                 &
                   - zdyr(j) *c05 *(  crlat(j+1,1) *zaiu(i,j+1)                &
                                    - crlat(j-1,1) *zaiu(i,j-1))
    ENDDO
  ENDDO

!-------------------------------------------------------------------------------
! Section 3:  Derivation of a smooth first guess geopotential increment field
!             (by solving the Poisson equation on an intermediate coarse grid)
!-------------------------------------------------------------------------------
! Section 3a: Collect the input wind increments fields (f*v) at node 0
!-------------------------------------------------------------------------------

  IF (irm == 0)  ALLOCATE ( zaiu_tot (ie_tot,je_tot)      , STAT=irm )
  IF (irm == 0)  ALLOCATE ( zaiv_tot (ie_tot,je_tot)      , STAT=irm )
  IF (irm /= 0) WRITE( yerrmsg,'(''ERROR at alloc zaiu/v_tot'',I5)' ) irm
  IF (irm /= 0) CALL model_abort (my_cart_id, 11931, yerrmsg, yroutine, irm )

  IF (num_compute > 1) THEN

    CALL gather_field ( zaiu , ie, je, zaiu_tot , ie_tot, je_tot, 0, irm )
!   =================

    IF (irm /= 0) WRITE( yerrmsg,'(''ERROR at gather_field'',I6)' ) irm
    IF (irm /= 0) CALL model_abort (my_cart_id, 11941, yerrmsg, yroutine, irm )

    CALL gather_field ( zaiv , ie, je, zaiv_tot , ie_tot, je_tot, 0, irm )
!   =================

    IF (irm /= 0) WRITE( yerrmsg,'(''ERROR at gather_field'',I6)' ) irm
    IF (irm /= 0) CALL model_abort (my_cart_id, 11942, yerrmsg, yroutine, irm )
  ELSE
    zaiu_tot = zaiu
    zaiv_tot = zaiv
  ENDIF

  ! initialise 'zdfi_tot'
  istat    = 0
  ALLOCATE ( zdfi_tot (ie_tot,je_tot) , stat=istat(1) )
  zdfi_tot = c0

!-------------------------------------------------------------------------------
! Section 3b: Bi-linear interpolation of the collected fields to a coarse grid
!-------------------------------------------------------------------------------

  IF (my_cart_id == 0) THEN

    ! set parameters for field size of the coarse grid
!   ige    = 70
!   jge    = 69
    ige    = NINT( REAL( ie_tot-1, wp) * dlon * 1.5_wp ) + 1
    jge    = NINT( REAL( je_tot-1, wp) * dlat * 1.5_wp ) + 1
    ige    = MAX( 5 , (ie_tot-1)/20 , MIN( (ie_tot-1)/5 , ige-1 ) ) + 1
    jge    = MAX( 5 , (je_tot-1)/20 , MIN( (je_tot-1)/5 , jge-1 ) ) + 1
    IF (ntstep == 0)  PRINT *,'coarse grid size for Poisson solver = ', ige, jge

    ! set other parameters
    ntpois_c = 1640
    zfpois   = 13._wp / 12._wp

    ALLOCATE( zrotv_c (ige,jge)      , fi_0 (ige,jge)           , stat=istat(2))
    ALLOCATE( zdxg(jge) , zdyg(jge)  , crlatg (0:jge+1)         , stat=istat(3))
    ALLOCATE( ug(ige,jge),vg(ige,jge), ug2(ige,jge),vg2(ige,jge), stat=istat(4))
!   PRINT*,'istat_2,3,4 = ', istat(2), istat(3), istat(4)

    ! initialisation of fields
    ug = c0 ;   vg   = c0 ;   zrotv_c = c0 ;   fi_0 = c0

    ! initialisation of variables
    dinkx    = REAL( ie_tot-1, wp) /REAL( ige-1, wp)
    dinky    = REAL( je_tot-1, wp) /REAL( jge-1, wp)
    dlatg    = dinky *dlat
    dlong    = dinkx *dlon
    odadlatg = c1/ (r_earth *dlat *degrad *dinky)
 
!   WRITE(*,'(a,5(e12.5,2x))') 'dinkx dinky dlatg dlong = '                    &
!                              ,dinkx,dinky,dlatg,dlong
 
    ! auxilliary fields
    DO jg = 0, jge+1
      y3         = startlat_tot + REAL( jg-1, wp) *dlatg
      crlatg(jg) = COS( y3 *degrad )
    ENDDO
    DO jg = 1, jge
      zdxg(jg)   = c1/ (r_earth *crlatg(jg) *dlong *degrad)
!CS: is 'crlatg' really needed for 'zdyg' ?
      zdyg(jg)   = c1/ (r_earth *crlatg(jg) *dlatg *degrad)
!     e1         = c1/ zdxg(jg)
!     e2         = c1/ zdyg(jg)
!     WRITE(*,'(a,3(e12.5,2x))') 'cos izdx izdy = ',crlatg(jg),e1,e2
    ENDDO
 
! do the bi-linear interpolation
! ------------------------------
 
    dvmax = c0
 
    DO jg = 1, jge
      y3 = dinky *REAL( jg-1, wp) + c1
      jn1 = MIN( INT( y3 +epsy, iintegers ) , je_tot-1 )
      jn2 = jn1+1
 
      DO ig = 1, ige
        x3 = dinkx *REAL( ig-1, wp) + c1
        in1 = MIN( INT( x3 +epsy, iintegers ) , ie_tot-1 )
        in2 = in1+1
 
        yn = y3 - REAL( jn1, wp)
        xn = x3 - REAL( in1, wp)
 
        ug (ig,jg) =  (c1-yn)* ( (c1-xn)* zaiu_tot(in1,jn1)                    &
                                +    xn * zaiu_tot(in2,jn1))                   &
                     +    yn * ( (c1-xn)* zaiu_tot(in1,jn2)                    &
                                +    xn * zaiu_tot(in2,jn2))
        vg (ig,jg) =  (c1-yn)* ( (c1-xn)* zaiv_tot(in1,jn1)                    &
                                +    xn * zaiv_tot(in2,jn1))                   &
                     +    yn * ( (c1-xn)* zaiv_tot(in1,jn2)                    &
                                +    xn * zaiv_tot(in2,jn2))

        dv        = SQRT( ug(ig,jg) *ug(ig,jg) + vg(ig,jg) *vg(ig,jg) )
        dvmax     = MAX( dvmax,dv )

      ENDDO  !ig
    ENDDO  !jg

!-------------------------------------------------------------------------------
! Section 3c: Compute balanced geopotential(*) increments on the coarse grid
!             (by solving the Poisson equation with an iterative algorithm)
!             (*): in fact: geopotential divided by Coriolis parameter
!-------------------------------------------------------------------------------

! compute rotation (vorticity) of the wind increment field on the coarse grid
! ---------------------------------------------------------------------------

    DO j = 2, jge-1
      DO i = 2, ige-1
        zrotv_c(i,j) =   zdxg(j) *c05 *(vg(i+1,j) - vg(i-1,j))                 &
                       - zdyg(j) *c05 *(  crlatg(j+1) *ug(i,j+1)               &
                                        - crlatg(j-1) *ug(i,j-1))
      ENDDO
    ENDDO

!   PRINT*,'zfpois = ',zfpois
    ! initialisation of variables used for diagnostics
    ! could also be used for a dynamic criterion to stop the iterations)
    ug  = c0
    vg  = c0
    ug2 = c0
    vg2 = c0

! solve Poisson equation iteratively, 
! looping in an alternative direction from one to the next iteration
! ------------------------------------------------------------------

    DO itpois = 1 , ntpois_c
      np1 = MOD( itpois, 4 )
 
    ! (np1 == 1): DO j = 2,jge-1    ; DO i = 2,ige-1    ! inner loop west-east
    ! (np1 == 2): DO j = jge-1,2,-1 ; DO i = ige-1,2,-1 ! inner loop east-west
    ! (np1 == 3): DO i = 2,ige-1    ; DO j = 2,jge-1    ! inner loop south-north
    ! (np1 == 0): DO i = ige-1,2,-1 ; DO j = jge-1,2,-1 ! inner loop north-south
      np2 = MOD( np1,2 )
      jn1 = np2* 2       + (1-np2)* (jge-1)
      jn2 = np2* (jge-1) + (1-np2)* 2
      in1 = np2* 2       + (1-np2)* (ige-1)
      in2 = np2* (ige-1) + (1-np2)* 2
      jnd = np2          - (1-np2)
 
      IF ((np1 == 1) .OR. (np1 == 2)) THEN
        DO j = jn1, jn2, jnd
          DO i = in1, in2, jnd
            e1 = zdxg(j)*zdxg(j)
            e2 = c05*(crlatg(j+1)+crlatg(j)) *zdyg(j) *odadlatg
            e3 = c05*(crlatg(j-1)+crlatg(j)) *zdyg(j) *odadlatg
            e4 = c1/(c2*e1+e2+e3)
 
            CALL poisson_coarse
!           ===================
          ENDDO
        ENDDO
      ELSEIF ((np1 == 3) .OR. (np1 == 0)) THEN
        DO i = in1, in2, jnd
          DO j = jn1, jn2, jnd
            e1 = zdxg(j)*zdxg(j)
            e2 = c05*(crlatg(j+1)+crlatg(j)) *zdyg(j) *odadlatg
            e3 = c05*(crlatg(j-1)+crlatg(j)) *zdyg(j) *odadlatg
            e4 = c1/(c2*e1+e2+e3)
 
            CALL poisson_coarse
!           ===================
          ENDDO
        ENDDO
      ENDIF
 
! for diagnostics
! ---------------
 
      DO j = 2, jge-1
      DO i = 2, ige-1
        ug2(i,j) = -c05* (fi_0(i  ,j+1) - fi_0(i  ,j-1)) *odadlatg
        vg2(i,j) =  c05* (fi_0(i+1,j)   - fi_0(i-1,j)  ) *zdxg(j)
      ENDDO; ENDDO
 
      divmax = c0
      umax   = c0
 
      DO j = 2, jge-1
      DO i = 2, ige-1
        dv     = SQRT( (ug2(i,j)-ug(i,j))**2 + (vg2(i,j)-vg(i,j))**2 )
        divmax = MAX( divmax,dv )
        dv     = SQRT( ug2(i,j)**2 + vg2(i,j)**2 )
        umax   = MAX( umax,dv )
      ENDDO; ENDDO
 
      ug = ug2
      vg = vg2
 
    ENDDO ! itpois
 
!-------------------------------------------------------------------------------
! Section 3d: Interpolate geopotential increments to the fine model grid
!-------------------------------------------------------------------------------

    dinkx    = REAL( ige-1, wp) / REAL( ie_tot-1, wp)
    dinky    = REAL( jge-1, wp) / REAL( je_tot-1, wp)

    DO j = 1, je_tot
      y3 = dinky *REAL(j-1,wp) + c1
      jg = MIN( INT( y3 +epsy, iintegers ) , jge-1 )
      jgp1 = jg + 1
 
      DO i = 1, ie_tot
        x3 = dinkx *REAL(i-1,wp) + c1
        ig = MIN( INT( x3 +epsy, iintegers ) , ige-1 )
        yn = y3 - REAL( jg, wp)
        xn = x3 - REAL( ig, wp)
 
        zdfi_tot (i,j) =  (c1-yn)* ( (c1-xn)* fi_0(ig  ,jg  )                  &
                                    +    xn * fi_0(ig+1,jg  ))                 &
                         +    yn * ( (c1-xn)* fi_0(ig  ,jgp1)                  &
                                    +    xn * fi_0(ig+1,jgp1))
      ENDDO  ! i
    ENDDO  ! j
 
!   IF ((lwonl2) .AND. (MOD( ntstep,60 ) < 6)) THEN
!     DO i = 1,12
!       WRITE( nupr,*) 'zdfi_tot = ',zdfi_tot(i,329), zaiu_tot(i,329)          &
!                                  , zaiv_tot(i,329)
!     ENDDO
!     DO i = 1,4
!       WRITE( nupr,*) 'fi_0  = ',fi_0(i,33), fi_0(i,34), fi_0(i,35)
!       WRITE( nupr,*) 'zrotv_c = ',zrotv_c(i,33), zrotv_c(i,34), zrotv_c(i,35)
!     ENDDO
!   ENDIF
 
    DEALLOCATE( ug, vg, ug2, vg2, zrotv_c, fi_0,  stat=istat(5))
    DEALLOCATE( zdxg, zdyg, crlatg, stat=istat(6) )
 
  ENDIF  !  my_cart_id == 0
 
!-------------------------------------------------------------------------------
! Section 3e: Distribute the global geopotential field to the local sub-domains
!-------------------------------------------------------------------------------

  DEALLOCATE ( zaiu_tot , STAT=irm )
  DEALLOCATE ( zaiv_tot , STAT=irm )

  CALL distribute_field ( zdfi_tot(:,:)  , ie_tot, je_tot                      &
                        , zdfi    (:,:,1), ie,     je,     0,   istat(7) )
! =====================

  DEALLOCATE ( zdfi_tot , STAT=irm )

!-------------------------------------------------------------------------------
! Section 4: Computation of the geopotential increment field by a Poisson solver
!-------------------------------------------------------------------------------

  ! set paramters
  ntpois = 50
  epsv   = 1.667e-5_wp

  ! initialisations
  ALLOCATE ( upo(ie,je), vpo(ie,je), upn(ie,je), vpn(ie,je), stat=istat(8) )

  upo = c0 ;   vpo = c0 ;   dvmax       = c0
  upn = c0 ;   vpn = c0 ;   zdfi(:,:,2) = zdfi(:,:,1)
  it0 = 0  ;   it1 = 1  ;   it2 = 2
  xnum   = REAL( (ie-2)*(je-2), wp)
 
  ! get first guess geostrophic wind field (for diagnostics)
  DO j = 2, je-1
  DO i = 2, ie-1
    upo (i,j) = -c05* (zdfi(i  ,j+1,it1) - zdfi(i  ,j-1,it1)) *edadlat
    vpo (i,j) =  c05* (zdfi(i+1,j  ,it1) - zdfi(i-1,j  ,it1)) *zdxr(j)
  ENDDO; ENDDO
 
! compute balanced geopotential(*) increments on the model grid
! (by solving the Poisson equation with an iterative algorithm)
! (*): in fact: geopotential divided by Coriolis parameter
! -------------------------------------------------------------
 
  DO itpois = 1, ntpois
 
    DO j = 2, je-1
      e1 = zdxr(j)*zdxr(j)
      e2 = crlat(j  ,2) *zdyr(j) *edadlat
      e3 = crlat(j-1,2) *zdyr(j) *edadlat
      e4 = c1 / (c2*e1 + e2 + e3)
      DO i = 2, ie-1
        zdfi(i,j,it2) =   e1*e4*(zdfi(i+1,j,it1) + zdfi(i-1,j,it1))            &
                        + e2*e4* zdfi(i,j+1,it1)                               &
                        + e3*e4* zdfi(i,j-1,it1)                               &
                        -    e4* zrotv(i,j)
      ENDDO
    ENDDO
 
    kzdims(1:24)=(/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)

    CALL exchg_boundaries ( 1, sendbuf, isendbuflen, imp_reals, icomm_cart,       &
                           num_compute, ie, je, kzdims, jstartpar, jendpar, 1,    &
                           nboundlines, my_cart_neigh, lperi_x, lperi_y, l2dim,   &
                           20000+nexch_tag,.FALSE., ncomm_type, izerror, yzerrmsg,&
                           zdfi )
!   =====================
 
! diagnostics
! -----------
 
    DO j = 2, je-1
    DO i = 2, ie-1
      upn (i,j) = -c05*(zdfi(i,j+1,it2) - zdfi(i,j-1,it2)) *edadlat
      vpn (i,j) =  c05*(zdfi(i+1,j,it2) - zdfi(i-1,j,it2)) *zdxr(j)
    ENDDO; ENDDO
 
    sum   = 0.0_wp
    dvmax = 0.0_wp
    umax  = 0.0_wp
 
    DO j = 2,je-1
    DO i = 2,ie-1
      dv    = SQRT( (upn(i,j)-upo(i,j))**2 + (vpn(i,j)-vpo(i,j))**2 )
      dvmax = MAX( dvmax,dv )
      sum   = sum + dv
      dv    = SQRT( upn(i,j)**2 + vpn(i,j)**2 )
      umax  = MAX( umax,dv )
    ENDDO; ENDDO
 
    sum = sum/xnum
 
    CALL global_values (dvmax,1,'MAX',imp_reals,icomm_cart,-1,yerror,ierror)
    CALL global_values (sum  ,1,'MAX',imp_reals,icomm_cart,-1,yerror,ierror)
    CALL global_values (umax ,1,'MAX',imp_reals,icomm_cart,-1,yerror,ierror)
!   ==================

    q11 = dvmax *1000._wp /(umax +1.e-12_wp)
    IF (lwonl2) THEN
      WRITE(nupr,'(a,i4,2x,2x,2(f10.6,2x))') 'itera q11 umax = ',itpois,q11,umax
    ENDIF
 
    upo = upn
    vpo = vpn
 
    it0 = it2
    it2 = it1
    it1 = it0
 
!   CALL global_values (q11  ,1,'MAX',imp_reals,icomm_cart,-1,yerror,ierror)
!   IF ((q11 < epsv) .AND. (itpois > 20)) EXIT
  ENDDO  !itpois
 
  DO j = 1,je
  DO i = 1,ie
    zdfi(i,j,it1) = zdfi(i,j,it1) *zfcc(i,j)
    zdfi(i,j,it2) = zdfi(i,j,it1)
  ENDDO; ENDDO
 
  IF (lwonl2) THEN
    zdfi_max = MAXVAL( zdfi(:,:,it2) )
    zdfi_min = MINVAL( zdfi(:,:,it2) )
    WRITE( nupr,*) 'zdfi_max = ',my_cart_id, zdfi_max, zdfi_min
  ENDIF
 
  DEALLOCATE ( upn, vpn, upo, vpo, stat=istat(9) )
 
!-------------------------------------------------------------------------------
! Section 5:  Conversion of geopotential increments to pressure increments.
!             Note: Since the non-isotropic part of the lateral weights, which
!                   reflect mainly the orography, has not been included in the
!                   wind analysis increments, the resulting geostrophic
!                   geopotential increments should be regarded as being
!                   representative for the mean sea level. They have to be
!                   converted to pressure increments and then extrapolated
!                   to the model orography level.
!-------------------------------------------------------------------------------
! Section 5a: Conversion of geopotential increments to pressure increments
!             at 'mean sea level' (more precisely: at half of the thickness
!             of the lowest model level).
!             This requires temperature and pressure at 'mean sea level'.
!-------------------------------------------------------------------------------

  psaigeo (:,:) = c0

  CALL calpmsl( zpmsl(:,:), pp(:,:,ke,nnew)+p0(:,:,ke), t(:,:,ke,nnew)         &
              , rho0(:,:,ke), dp0(:,:,ke), hhl(:,:,ke+1), ie, je, g, r_d )
! ============

  zlapse = 0.0065_wp
  r_d_r  = c1 / r_d
  DO   j = jstart, jend
    DO i = istart, iend
! extrapolate temperature to mean sea level using 'zlapse' (as in calpmsl)
      ztmsl   (i,j)  =  t(i,j,ke,nnew) + zlapse *hhl(i,j,ke+1)
! convert geopotential increments into pressure increments
      psaigeo (i,j)  =  zdfi(i,j,it2) *zpmsl(i,j) /ztmsl(i,j) *r_d_r
    ENDDO
  ENDDO

!-------------------------------------------------------------------------------
! Section 5b: Extrapolation of the pressure increments from mean sea level
!             to the lowest model level by applying the correction factor
!             used in the lateral spreading of pressure observation increments
!             according to the temperature correction.
!-------------------------------------------------------------------------------

  zmaxai_1 = -1.e9_wp
  zminai_1 =  1.e9_wp
  DO   j = jstart, jend
    DO i = istart, iend
      zpke_v  (i,j) = p0(i,j,ke) + pp(i,j,ke,nnew)
      zsprcor (i,j) = zpke_v(i,j) / zpmsl(i,j)
      cfitop_v(i,j) = c1
      IF (ntpscor >= 3) THEN
        cfitop_v(i,j) = zpke_v(i,j) / zpmsl(i,j) * ztmsl(i,j) / t(i,j,ke,nnew)
      ENDIF
      IF (ntpscor >= 1) THEN
        zqmass  = (zpke_v(i,j) - ptpstop) / (zpmsl(i,j) - ptpstop)
        IF (MOD(ntpscor,2) == 1) THEN
          zqmass2 = zqmass *zqmass
          zqmass3 = zqmass *zqmass2
          zsprcor (i,j) =  cfitop_v(i,j) * zqmass2                         &
                         * EXP( (c1-zqmass3)*0.125_wp / MAX(c1,zqmass3))
        ELSE
          zsprcor (i,j) =  cfitop_v(i,j) * zqmass * MAX(c1 , zqmass)       &
                         * (c05*(c1+zqmass)) **NINT(SIGN( c1, c1-zqmass ))
        ENDIF
      ENDIF
      psaigeo (i,j)  =  psaigeo(i,j) * zsprcor(i,j)
      IF (lwonl2) THEN
        zmaxai_1 = MAX( psaigeo(i,j), zmaxai_1 )
        zminai_1 = MIN( psaigeo(i,j), zminai_1 )
      ENDIF

!-------------------------------------------------------------------------------
! Section 5c: Reduction of pressure increments according to parameter 'qgeops'.
!-------------------------------------------------------------------------------

      psaigeo (i,j)  =  psaigeo(i,j) * qgeops
    ENDDO
  ENDDO

  zmaxai =   MAXVAL(   psaigeo(:,:) )
  zminai = - MAXVAL( - psaigeo(:,:) )
  IF (lwonl2) THEN
    WRITE( nupr,*) 'max(psaigeo_1): ',my_cart_id, zmaxai_1, zminai_1
    WRITE( nupr,*) 'max(psaigeo_2): ',my_cart_id, zmaxai  , zminai
  ELSEIF (MAX( zmaxai,ABS(zminai) ) > c2) THEN
    PRINT       *, 'max(psaigeo_2): ',my_cart_id, zmaxai  , zminai
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure geostroph_ps_corr
!-------------------------------------------------------------------------------
!
CONTAINS
!
!-------------------------------------------------------------------------------
SUBROUTINE poisson_coarse
!
      fi_0(i,j)   =   (c1-zfpois)*fi_0(i,j)                                    &
                    +     zfpois *(  e1*e4*(fi_0(i+1,j)+fi_0(i-1,j))           &
                                   + e2*e4* fi_0(i,j+1)                        &
                                   + e3*e4* fi_0(i,j-1)                        &
                                   -    e4* zrotv_c(i,j))
!
END SUBROUTINE poisson_coarse
!-------------------------------------------------------------------------------

END SUBROUTINE geostroph_ps_corr

!===============================================================================

ELEMENTAL REAL (KIND=wp) FUNCTION fpvsw  ( zt, b1, b2w, b3, b4w )
  !----------------------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  zt, b1, b2w, b3, b4w
  !----------------------------------------------------------------------------
  ! Magnus formula for water:  input  'zt'   : temperature
  !                            output 'fpvsw': saturation water vapour pressure
  ! Magnus formula for ice  :  if constants 'b2i', 'b4i' for ice are used for
  !                            'b2w', 'b4w' in the call
  !----------------------------------------------------------------------------
  !
  fpvsw  =  b1 * EXP( b2w*(zt-b3)/(zt-b4w) )
  !
END FUNCTION fpvsw

!-------------------------------------------------------------------------------

ELEMENTAL REAL (KIND=wp) FUNCTION fpv2q  ( zpv, zp, rdv )
  !--------------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  zpv, zp, rdv
  !---------------------------------------------------------------------------
  ! specific humidity from water vapour pressure 'zpv' and air pressure 'zp'
  !   (rdv = r_d / r_v
  !        = gas constant for dry air / gas constant for water vapour )
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
  ! from specific humidity 'zqv' and air pressure 'zp'
  !   (rdv = r_d / r_v
  !        = gas constant for dry air / gas constant for water vapour )
  !---------------------------------------------------------------------------
  !
  fq2pv  =  MAX( epsy , zqv ) * zp / (rdv + zqv*(c1-rdv))
  !
END FUNCTION fq2pv

!-------------------------------------------------------------------------------

ELEMENTAL REAL (KIND=wp) FUNCTION rmod  ( ztti, ziv, epsy )
  !----------------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  ztti, ziv, epsy
  !---------------------------------------------------------------------------
  ! MOD function for positive reals
  !---------------------------------------------------------------------------
  !
  rmod  =  ztti  -  ziv * REAL(INT( ztti/ziv + epsy ), wp) + epsy
  !
END FUNCTION rmod

!-------------------------------------------------------------------------------

END MODULE src_nudging
