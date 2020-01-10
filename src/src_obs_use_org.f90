!+ Source module for organizing the tasks which make use of observations
!-------------------------------------------------------------------------------

MODULE src_obs_use_org

!-------------------------------------------------------------------------------
!
! Description:
!   The module "src_obs_use_org" organizes the tasks which make use of
!   observations (except for the gridded observations from Grib files,
!                 as used in latent heat nudging) :
!    - data assimilation by nudging towards direct observations;
!    - application of observation operators for KENDA (LETKF);
!    - writing NetCDF feedobs files for KENDA (LETKF) and for verification;
!    - writing YUVERIF files for verification;
!    - preparation of surface analyses (by doing the observation processing).
!   All these tasks require reading and processing of observations which is
!   also organised here.
!
!   This module contains the following procedures:
!    - organize_nudging :  called by 'organize_assimilation'
!
!   This module also contains an elemental function, formerly statement funct.:
!    - rmod             :  MOD function for positive reals
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
!  Initial release, by extraction from src_nudging.f90 , plus modifications:
!  - Modified interfaces of the calls for the observation (pre-)processing
!    and for writing the VOF (yuverif) file.
!  - New routine call for writing a NetCDF feedobs file (for KENDA / verif.).
!  - 'nexcess' introduced for local array bound checks.
!  - Frequency of opening and closing files (for flushing) reduced.
!  - Some variables moved from module 'data_nudge_all' to 'data_obs_lib_cosmo'.
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer, Ulrich Schaettler
!  Replaced qx-variables by using them from the tracer module
!  Introduced nexch_tag for MPI boundary exchange tag to replace ntstep
! V4_27        2013/03/19 Christoph Schraff, Astrid Kerkweg, Ulrich Schaettler
!  'lrefend' for reference time of model run added to argument list of
!  'obs_print_vof'.
!  Introduced MESSy interface (AK)
! V4_28        2013/07/12 Christoph Schraff
!  Statement function 'rmod' replaced by elemental function.
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
! V5_1         2014-11-28 Christoph Schraff, Oliver Fuhrer
!  (mruntyp == 2) added as condition for 'lrefend'.
!  Output 'verifn' done also if (.not.lnudge) but (lverif).
!  Replaced ireals by wp (working precision) (OF)
! V5_4         2016-03-10 Christoph Schraff
!  Removal of (all dependencies on) the interface to read obs from AOF files.
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
!   istart,       & ! start index for the forecast of w, t, qv, qc and pp
!   iend,         & ! end index for the forecast of w, t, qv, qc and pp

!   meridional direction
!   jstart,       & ! start index for the forecast of w, t, qv, qc and pp
!   jend,         & ! end index for the forecast of w, t, qv, qc and pp
    jstartpar,    & ! start index for computations in the parallel program
    jendpar,      & ! end index for computations in the parallel program

! 4. constants for the horizontal rotated grid and related variables
! ------------------------------------------------------------------

    pollon,       & ! longitude of the rotated north pole (in degrees, E>0)
    pollat,       & ! latitude of the rotated north pole (in degrees, N>0)
    dlon,         & ! grid point distance in zonal direction (in degrees)
    dlat,         & ! grid point distance in meridional direction (in degrees)
    startlat_tot, & ! rotated latitude  \  of the lower left grid point of the
    startlon_tot, & ! rotated longitude /  total domain (in degrees, N,E>0)

! 4. variables for the time discretization and related variables
! --------------------------------------------------------------

    dt,           & ! long time-step
!   ed2dt,        & ! 1 / (2 * dt)
    dt2,          & ! 2 * dt
    dtdeh,        & ! dt / 3600 seconds
    idt_qv, idt_qc

! end of data_modelconfig

!-------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 3. prognostic variables                                             (unit)
! -----------------------

    u           , & ! zonal wind speed                                ( m/s )
    v           , & ! meridional wind speed                           ( m/s )
    t           , & ! temperature                                     (  k  )
    pp              ! deviation from the reference pressure           ( pa  )

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
!   nnow,         & ! corresponds to ntstep
    nnew,         & ! corresponds to ntstep + 1

! 7. additional control variables
! -------------------------------

    nvers,        & ! version number of experiment for documentation
    l2tls,        & ! forecast with 2-TL integration scheme
    ltime,        & ! detailled timings of the program are given
!   lreproduce,   & ! the results are reproducible in parallel mode
!   ldump_ascii,  & ! for flushing (close and re-open) the ASCII files
    lperi_x,      & ! if lartif_data=.TRUE.: periodic boundary conditions in x-dir.
                    ! or with Davies conditions (.FALSE.) 
    lperi_y,      & ! if lartif_data=.TRUE.: periodic boundary conditions in y-dir.
                    ! or with Davies conditions (.FALSE.) 
    l2dim,        & ! 2 dimensional runs

! 11. controlling ensemble mode (EPS) 
! -----------------------------------

    leps,         & ! switch ensemble mode on/off
    iepsmem         ! ID of the member in the ensemble (ID >= 0)

! end of data_runcontrol 

!-------------------------------------------------------------------------------

USE data_parallel,      ONLY :  &
    num_compute,   & ! number of compute PEs
    nboundlines,   & ! number of boundary lines of the domain for which
                     ! no forecast is computed = overlapping boundary
                     ! lines of the subdomains

!   ldatatypes,    & ! if .TRUE.: use MPI-Datatypes for some communications
!   ltime_barrier, & ! if .TRUE.: use additional barriers for determining the
                     ! load-imbalance
    ncomm_type,    & ! type of communication
    my_cart_id,    & ! rank of this subdomain in the cartesian communicator
    my_cart_neigh, & ! neighbors of this subdomain in the cartesian grid
!   isubpos,       & ! positions of the subdomains in the total domain. Given
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

    ydate_ini           ,& ! start of the forecast
                           ! yyyymmddhh (year, month, day, hour)
    ! NetCDF global attributes
    yncglob_institution ,& ! originating center name
    yncglob_source      ,& ! program name and version
!   yncglob_title       ,& ! title string for the output
!   yncglob_references  ,& ! URL, report etc.
    ! for directory name of grib output
    root                   ! pointer to the root of gribout namelists

! end of data_io

!-------------------------------------------------------------------------------

USE data_nudge_all , ONLY :   &

! 1. parameters and related variables
! -----------------------------------

    lwonl        ,& ! .TRUE if grid pt. (ionl ,jonl ) lies in the local domain
    lcroot       ,& ! .TRUE if (my_cart_id == 0) (for print on std. output)
    lfirst       ,& ! .TRUE if 'organize_nudging' is called the first time
    acthr        ,& ! actual forecast hour
    aiwthr       ,& ! model time [hrs] for which temporal weights are valid
    ltnudge      ,& ! nudging to be done at current timestep
    ltverif         ! verification to be done at current timestep

USE data_nudge_all , ONLY :   &

! 2. namelist variables controlling the data assimilation
! -------------------------------------------------------

    lnudge       ,& ! .f. : on - off switch for nudging
    lverif       ,& ! .f. : on - off switch for verification
    mruntyp      ,& ! -1  : type of current model run used for increments in VOF
    mveripr      ,& ! 3   : type of verification/observation file(s) written
    nudgsta      ,& ! 0   : start of nudging period in timesteps
    nudgend      ,& ! 0   : end of nudging period in timesteps
    nversta      ,& ! 0   : start of verification period in timesteps
    nverend      ,& ! 0   : end of verification period in timesteps
    hversta      ,& ! start of verification period in 'model integration hours'
    hverend      ,& ! end of verification period in 'model integration hours'
    tconbox      ,& ! 6*dt: timestep [s] for computing analysis increments
                    !       (i.e. time box of constant analysis increments)
    mpsgcor      ,& ! 1   : mode to apply geostrophic pressure correction
    khumbal      ,& ! 100 : range around convectively precipitating grid pts, at
                    !       which specified (not relative) humidity is preserved
    dtqc         ,& ! timestep (in [s]) for the threshold quality control
    qcc          ,& ! constant parts of the quality control thresholds (='QCT')
    qcvf         ,& ! multiplication factor to vertically varying part of QCT
    qccsu        ,& ! constant parts of the quality control thresholds (='QCT')
    maxmlo       ,& !  max. number of multi-level reports in the total domain
    maxsgo       ,& !  max. number of (surface-level and upper-air) single-level
                    !                 reports within the total domain
    maxgpo       ,& !  max. number of GPS reports within total domain
    maxtvo       ,& !  max. number of sat retrievals within total domain
    maxmlv       ,& !  100   : max. number of obs levels in multi-level reports
    mxfrep       ,& !   -1   : max. number of reports in NetCDF feedobs file
    mxfobs       ,& !   -1   : max. number of observations in feedobs file
    ionl         ,& ! 167    : / grid point coordinates
    jonl            ! 103    : \ for standard output on nudging

USE data_nudge_all , ONLY :   &

! 6. Miscellany
! -------------

!   ltvprcs      ,& ! .TRUE if 1dvar       processing at current timestep
    liwvssc      ,& ! .t. : spatial consistency check of IWV performed
    qctfpr          ! for VOF output only: time factor for QC thresholds

! end of data_nudge_all

!-------------------------------------------------------------------------------

USE data_obs_lib_cosmo , ONLY :   & 

! 1. General parameters
! ---------------------

    c0         ,& ! standard real constant 0.0
    c1         ,& ! standard real constant 1.0
    c2         ,& ! standard real constant 2.0
    c05        ,& ! standard real constant 0.5
    c3600      ,& ! standard real constant 3600.0
    epsy       ,& ! = 1.E-8_wp : commonly used very small value > 0

! 4. I/O device numbers for obs processing / nudging and file names
! -----------------------------------------------------------------

    nupr          ! unit number of file for all the remaining information

! end of data_obs_lib_cosmo

!-------------------------------------------------------------------------------

USE environment,              ONLY :  &
      exchg_boundaries, & ! performs the boundary exchange between 
                          ! neighboring processors
      comm_barrier,     & ! explicit synchronization point
      model_abort         ! aborts the program in case of errors

!-------------------------------------------------------------------------------

USE time_utilities,           ONLY :  get_timings, i_barrier_waiting_nud,  &
     i_communications_nud, i_nudging_corr, i_local_info, i_obs_processing, &
     i_spread_ps_pre, i_spread_multi_level, i_spread_upper_air,            &
     i_spread_surface

!-------------------------------------------------------------------------------

! USE src_obs_cdfin_util,       ONLY :  &
!   obs_gather_buffer    ! gathers elements from buffer arrays from all nodes

!-------------------------------------------------------------------------------

USE src_obs_proc_cdf,       ONLY:    &
  obs_org_cdf_proc, obs_cdf_print_caut_neff

USE src_obs_print_vof,      ONLY:    &
  obs_print_vof

USE src_obs_cdfout_feedobs, ONLY:    &
  obs_org_cdfout_feedobs

USE src_gather_info,        ONLY:    &
  gather_local_info, gather_spread_aux, gather_nudge_aux, gather_varia

USE src_sing_local,         ONLY:    &
  local_info_aux, local_sort_reports,                                          &
  ps_local_info, surf_local_info, upair_local_info, ps_spatial_check

USE src_mult_local,         ONLY:    &
  mult_org_localinfo, mult_iwv_horicheck

USE src_sing_spread,        ONLY:    &
  ps_spreading, surf_org_spread, upair_org_spread

USE src_mult_spread,        ONLY:    &
  mult_org_spread

USE src_nudging,            ONLY:    &
  nudge_horiz_wind, nudge_humid_mass, ps_temperatur_corr, geostroph_ps_corr

USE src_tracer,             ONLY: trcr_get, trcr_errorstr

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
!+ Module procedure in "src_obs_use_org" for organizing the use of observations
!-------------------------------------------------------------------------------

SUBROUTINE organize_nudging

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "src_obs_use_org" is the driving routine
!   for the tasks which make use of observations
!   (except for the gridded observations from Grib files, as used in LHN) :
!    - data assimilation by nudging towards direct observations;
!    - application of observation operators for KENDA (LETKF);
!    - writing NetCDF feedobs files for KENDA (LETKF) and for verification;
!    - writing YUVERIF files for verification;
!    - preparation of surface analyses (by doing the observation processing).
!   All these tasks require reading and processing of observations which is
!   also organised here.
!
! Method:
!   Basically, this is the driving routine for the data assimilation by
!   observation nudging. All sections in this routine are required for the
!   nudging (plus writing of files for verification).
!   For other tasks, only some of the sections are required and run.
!   1. Reading observations from AOF (if selected) and observation processing,
!      only when (usually once per hour) there may exist data, which have not
!      yet been read but would start influencing the current timestep.
!   2. Checking, if new analysis increments are to be computed at the current
!      timestep. If so, determination of the time for which the temporal
!      weights are valid, else omission of steps 3 - 7 .
!   3. Reading observations from NetCDF Observation Input files (if selected)
!      and observation processing.
!      (In the called routine, obs processing typically only once per hour.)
!   4. Computation of the observation increments by applying the observation
!      forward operators, and of all the other required local information
!      on the observations and from their location, including (the results of)
!      the threshold quality control ('first guess check').
!   5. Writing obserations, simulated observations, flags etc. onto
!      NetCDF feedobs and / or ASCII VOF (YUVERIF) files.
!   6. Nudging only: Distribution of the local information from each node
!                    to all nodes (if the program is run in parallel mode).
!   7. Nudging only: Spreading (extrapolation and weighting) of the observation
!                    increments to the target model grid points.
!   8. Nudging only: Computation of the analysis increments, and execution of
!                    the nudging equations by updating the prognostic variables.
!   Steps 7 and 8 are done model level by model level bottom up.
!   Steps (sections) 3, 4, 7 do not contain any communication between nodes.
!
! Current Code Owner:  Christoph Schraff
! Written by        :  Christoph Schraff, DWD  (original version: 13.06.97)
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

! Subroutine arguments: None
! --------------------

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND = iintegers) , SAVE :: &
    nexce_rep  =  0  ,& ! number of reports exceeding NetCDF feedback file size
    nexce_bdy  =  0     ! number of obs     exceeding NetCDF feedback file size

  REAL    (KIND=wp   )     ::  &
    ailast           ,& ! time of the last timestep within the current AITB
                        !   (AITB = analysis increment time box)
    ainext              ! earliest possible starting time of the next AITB

  LOGICAL                  ::  &
    lgetai           ,& ! .TRUE if new analysis increment are to be computed
                        !          and used at the current timestep
    lconai           ,& ! TRUE if analysis increments are const during 'tconbox'
    lastprc          ,& ! .TRUE if at the current timestep, the observation
                        !          processing (obs_org_cdf_proc) is called
                        !          for the last time in this COSMO run
    lrefend             ! .TRUE if verification reference time in NetCDF feedobs
                        !          file is set to end of model integration time

  CHARACTER (LEN=80)       ::  &
    yaitext             ! diagnostics on constant analysis increments

  INTEGER (KIND=iintegers) ::  &
    nulast           ,& ! timestep at which statistics is printed
    nexproc          ,& ! timestep at which the observation processing is called
                        !   for the next time
    nailast          ,& ! last timestep at which analysis increments are re-
                        !   computed (and obs processing (for NetCDF) is called)
    nla              ,& ! smallest possible value for 'nailast' acc. to tconbox 
    nli              ,& ! loop index
    iensmem          ,& ! ensemble member ( -1: deterministic)
    nexcess (5) = 0  ,& ! number of reports (obs incr.) in excess of array size
    max_rep          ,& ! max. number of reports in NetCDF feedobs file (FOF)
    max_body         ,& ! max. number of obs in NetCDF feedobs file (FOF)
    kk               ,& ! index of current vertical model level
    kzdims(24)       ,& ! vertical dimensions of the fields for exchange
!   i,j, itot, jtot  ,& ! for printout
    izerror             ! error status

  CHARACTER (LEN=255)      ::  yzerrmsg
  CHARACTER (LEN=25)       ::  yzroutine

! INTEGER (KIND=iintegers) ::  ncount, iicount, ircount, ijtot
! INTEGER (KIND=iintegers), ALLOCATABLE ::  ibufloc(:)
! INTEGER (KIND=iintegers), POINTER     ::  ibuftot(:)
! REAL    (KIND=wp)       , ALLOCATABLE ::  rbufloc(:)
! REAL    (KIND=wp)       , POINTER     ::  rbuftot(:)

! Local (automatic) arrays: None
! -------------------------
!
! Tracer pointers:
! -----------------
  REAL (KIND=wp)     , POINTER  :: &
    qv(:,:,:)=> NULL()         ,&   ! QV at nnew
    qc(:,:,:)=> NULL()              ! QC at nnew
!
!------------ End of header ----------------------------------------------------

 
!-------------------------------------------------------------------------------
! Begin Subroutine organize_nudging
!-------------------------------------------------------------------------------

izerror   = 0_iintegers
yzerrmsg  = ''
yzroutine = 'organize_nudging'

kzdims(:) = 0_iintegers

!-------------------------------------------------------------------------------
! Section 1: Initializations
!-------------------------------------------------------------------------------

IF ((ntstep == 0) .AND. (.NOT. l2tls)) THEN
  dt    = c2 * dt
  dt2   = c2 * dt
  dtdeh = dt / c3600
ENDIF

IF (lfirst)  CALL local_info_aux ( 1 )
!            ===================

ltnudge  =  (lnudge) .AND. (ntstep >= nudgsta) .AND. (ntstep <= nudgend)
ltverif  =  (lverif) .AND. (ntstep >= nversta) .AND. (ntstep <= nverend)

IF ((ltnudge) .OR. (ltverif)) THEN

  IF (num_compute > 1) THEN
    IF (ltime) THEN
      CALL get_timings (i_local_info, ntstep, dt, izerror)
      CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
      CALL get_timings (i_barrier_waiting_nud, ntstep, dt, izerror)
    ENDIF

! Exchange of 1 row of boundary data (of local domains)

    kzdims(1:24)=(/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)

    CALL exchg_boundaries ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart,       &
                           num_compute, ie, je, kzdims, jstartpar, jendpar, 1,    &
                           nboundlines, my_cart_neigh, lperi_x, lperi_y, l2dim,   &
                           20000+nexch_tag,.FALSE., ncomm_type, izerror, yzerrmsg,&
                           u(:,:,:,nnew), v(:,:,:,nnew) )
!   =====================

    IF (ltime) CALL get_timings (i_communications_nud, ntstep, dt, izerror)
  ENDIF

  acthr   = ntstep *dtdeh
  lconai  = (tconbox > dt+epsy)

  CALL gather_varia ( 1 , lconai )
! =================

! CALL local_info_aux ( 1 )
! ===================

  ! retrieve the required microphysics tracers (at nnew)
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnew, ptr = qv)
  IF (izerror /= 0_iintegers) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nnew, ptr = qc )
  IF (izerror /= 0_iintegers) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

!-------------------------------------------------------------------------------
! Section 2: Check if new analysis increments are required at the current time-
!            step, and if so, determination of the time for which the temporal
!            weights are valid
!-------------------------------------------------------------------------------

  IF (lconai) THEN
    lgetai  =   (lfirst) .OR. (rmod( ntstep*dt,tconbox ) < dt)
    IF (lgetai) THEN
! starting time [hrs] of the next (regular) 'analysis increment time box' (AITB)
      ainext = REAL(INT( c3600 *acthr /tconbox + epsy + c1 ),wp) * tconbox / c3600
! time [hrs] of the last timestep within the present 'analysis incr. time box'
      ailast = REAL(INT( c3600 *ainext /dt     - epsy      ),wp) * dt      / c3600
! mean time [hrs] of all timesteps in the present 'analysis increment time box'
      aiwthr = c05 * (acthr + ailast)
      IF (((lcroot) .OR. (lwonl)) .AND. (lfirst)) THEN
!       hconbox = tconbox / c3600
        WRITE( yaitext,'('' !! Analysis Increments ("AI") held constant ''     &
                       &,''during time boxes of ca'',F6.3,'' hours'')' )       &
               tconbox / c3600
        IF (lcroot) PRINT       *
        IF (lcroot) PRINT       '(A)' , yaitext
        IF (lwonl ) WRITE( nupr,'(A)' ) yaitext
      ENDIF
      IF (((lcroot) .OR. (lwonl)) .AND. (acthr <= c2)) THEN
        WRITE( yaitext,'(4X,''hour: AI-box: ['',F5.3,'','',F5.3,''], mean:''   &
                       &,F6.3,'', next: AI:'',F5.2)') &
               acthr, ailast, aiwthr, ainext
        IF (lcroot) PRINT       '(A)' , yaitext
        IF (lwonl ) WRITE( nupr,'(A)' ) yaitext
      ENDIF
! next timestep when analysis increments are computed again
      nexproc = INT( c3600 *ainext /dt - epsy + c1 )
    ENDIF
  ELSE
    lgetai = .TRUE.
    aiwthr = acthr
    nexproc = ntstep + 1
  ENDIF

!-------------------------------------------------------------------------------
! Section 3: Observation Processing (if reading from NetCDF files)
!-------------------------------------------------------------------------------

  IF (lgetai) THEN
    nulast  =  MIN( INT( MAX( nverend , nudgend ) ,iintegers) , nstop )
    lastprc =  (nexproc > nulast)
  ENDIF

  IF (lgetai) THEN

    IF (ltime) CALL get_timings (i_local_info, ntstep, dt, izerror)

    ! last timestep 'nailast' when lgetai =.true.
    ! (i.e. when analysis increments are re-computed and obs processing called)
      ! nailast must be at nulast or < (tconbox/dt) timesteps before nulast
    nla     = nulast - INT( tconbox/dt -epsy )
    nailast = ntstep
    DO nli = nla , nulast
      IF (rmod( nli*dt,tconbox ) < dt)  nailast = nli
    ENDDO

    CALL obs_org_cdf_proc ( aiwthr , nailast )
!   =====================

    IF (ltime) CALL get_timings (i_obs_processing, ntstep, dt, izerror)
  ENDIF

!-------------------------------------------------------------------------------
! Section 4: Computation of all the required local information
!            on the observations and parameters from their location
!-------------------------------------------------------------------------------
!            Part I: Main part including computation of observation increments,
!                              without communication between PE's
!-------------------------------------------------------------------------------

  IF (lgetai) THEN

! preparations and sorting of local reports
! -----------------------------------------

    CALL local_info_aux ( 2 )
!   ===================

    CALL local_sort_reports
!   =======================

!-------------------------------------------------------------------------------
! surface pressure observations   (Note: ps_local_info has to be called prior to
! -----------------------------          the other *local_info because of usage
!                                        of mo??hd(.,nhqcfw) )

    CALL ps_local_info ( nexcess(4) )
!   ==================

! upper-air single-level observations
! -----------------------------------

    CALL upair_local_info ( nexcess(2) )
!   =====================

! surface-level observations (without surface pressure observations)
! --------------------------

    CALL surf_local_info ( nexcess(3) )
!   ====================

! multi-level observations
! ------------------------

    CALL mult_org_localinfo ( nexcess(1) , nexcess(5) )
!   =======================

! printing control output
! -----------------------

    CALL local_info_aux ( 3 )
!   ===================

!-------------------------------------------------------------------------------
!             Part II: Secondary part with communication between PE's
!-------------------------------------------------------------------------------

! spatial consistency checks
! --------------------------

    CALL gather_spread_aux ( 0 )
!   ======================

    IF (ltime) CALL get_timings (i_local_info, ntstep, dt, izerror)

    IF (num_compute > 1)  CALL gather_local_info ( 0 , nexcess )
!                         ======================

    CALL gather_varia ( 2 )
!   =================

    CALL ps_spatial_check
!   =====================

    IF (liwvssc)  CALL mult_iwv_horicheck
!                 =======================

!-------------------------------------------------------------------------------
! Section 5: Wrting observations (incl. simulated observations, flags etc.)
!            onto NetCDF feedobs and / or ASCII VOF (YUVERIF) files.
!            Clean-up of section 4.  Requires communication.
!-------------------------------------------------------------------------------

! printing of observations on ASCII file for verification
! -------------------------------------------------------

    IF ((lverif) .AND. (mveripr/2  >= 1)) THEN
      !   ygrbt = end   of model run : for NUDGE runs only
      !   ygrbt = start of model run : for FORECAST (LETKF) or NUDGECAST
      lrefend  =  (lnudge) .AND. (nverend >= nstop) .AND. (mruntyp == 2)

      CALL obs_print_vof ( lastprc, lrefend, ydate_ini, hversta, hverend       &
                         , mruntyp, mveripr                                    &
                         , num_compute, my_cart_id, icomm_cart, imp_integers   &
                         , dtqc, qcc, qcvf, qccsu, qctfpr                      &
                         , ie_tot, je_tot, ke, pollat, pollon, dlat, dlon      &
                         , startlat_tot, startlon_tot )
!     ==================

    ENDIF

! printing of observations on NetCDF file for Ensemble Kalman Filter or verif.
! ----------------------------------------------------------------------------

    IF ((lverif) .AND. (MOD( mveripr, 2 ) == 1)) THEN
      !   veri_ref_time = end   of model run : for NUDGE runs only
      !   veri_ref_time = start of model run : for FORECAST (LETKF) or NUDGECAST
      lrefend  =  (lnudge) .AND. (nverend >= nstop) .AND. (mruntyp == 2)
      iensmem  =  iepsmem
      IF (.NOT. leps) iensmem = -1
      max_rep  =  mxfrep
      max_body =  mxfobs
      IF (mxfrep <= 0)  max_rep  = MAX( 1 , NINT( hverend - hversta ) )        &
                                   * (maxmlo + maxsgo + 2*maxgpo + 2*maxtvo)
      IF (mxfobs <= 0)  max_body = MAX( 1 , NINT( hverend - hversta ) )        &
                                   * (    maxmlo *5 *maxmlv/3   +   maxsgo *8  &
                                      + 2*maxtvo *2 *ke         + 2*maxgpo *3 )

      CALL obs_org_cdfout_feedobs ( lastprc, lrefend, ydate_ini                &
                                  , hversta, hverend                           &
                                  , mruntyp, max_rep, max_body                 &
                                  , num_compute, my_cart_id, icomm_cart        &
                                  , imp_integers, imp_reals                    &
                                  , ie_tot, je_tot, ke, pollat, pollon         &
                                  , dlat, dlon, startlat_tot, startlon_tot     &
                                  , iensmem, nvers                             &
                                  , yncglob_institution, yncglob_source        &
                                  , root% ydir                                 &
                                  , nexce_rep, nexce_bdy )
!     ===========================

      IF (lastprc)  CALL obs_cdf_print_caut_neff ( nexce_rep, nexce_bdy        &
                                                 , max_rep, max_body )
!                   ============================

    ENDIF

! printing results of threshold quality control and single-level data
! -------------------------------------------------------------------

    CALL local_info_aux ( 4 )
!   ===================

    IF (num_compute > 1)  CALL gather_local_info ( 1 )
!                         ======================

    CALL gather_varia ( 3 )
!   =================

  ENDIF

  CALL local_info_aux ( 5 )
! ===================

  IF (ltime) CALL get_timings (i_local_info, ntstep, dt, izerror)

ENDIF   !   (ltnudge) .OR. (ltverif)

!-------------------------------------------------------------------------------
! Section 6: Prepare the spreading and the nudging by distributing
!            the required local information from each PE to all PE's
!            (this includes communication between PE's)
!-------------------------------------------------------------------------------

IF (ltnudge) THEN

! prepare humidity balancing for temperature nudging
! --------------------------------------------------

  IF (khumbal < 99)  CALL gather_nudge_aux
!                    =====================

  IF (ltime) CALL get_timings (i_nudging_corr, ntstep, dt, izerror)

  IF (lgetai) THEN

! (distribute) top level of temperature correction
! ------------------------------------------------

    IF (ltnudge) CALL ps_temperatur_corr ( 0 )
!                =======================

    IF (ltime) CALL get_timings (i_nudging_corr, ntstep, dt, izerror)

! distribute the local information related to the active observations
! -------------------------------------------------------------------

    IF ((num_compute > 1) .AND. (ltnudge))  CALL gather_local_info ( 2 )
!                                           ======================

!-------------------------------------------------------------------------------
! Section 7: Spreading (extrapolation) of the observation increments
!            (this section is without communication between PE's)
!-------------------------------------------------------------------------------

! preparations
! ------------

    CALL gather_varia ( 4 )
!   =================

    CALL gather_spread_aux ( 1 )
!   ======================

! spreading 'surface' pressure data (incl. computation of analysis increments)
! ---------------------------------

    CALL ps_spreading
!   =================

    IF (ltime) CALL get_timings (i_spread_ps_pre, ntstep, dt, izerror)

! temperature correction 'implied' by 'surface' pressure analysis increments
! --------------------------------------------------------------------------

    CALL ps_temperatur_corr ( 1 )
!   =======================

    IF (ltime) CALL get_timings (i_nudging_corr, ntstep, dt, izerror)

  ENDIF   !   (lgetai)

! bottom-up vertical loop over the model levels
! ---------------------------------------------
! ______________________________________________
  loop_over_vertical_levels:  DO kk = ke , 1 , -1
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IF (lgetai) THEN

      CALL gather_spread_aux ( 2 , kk )
!     ======================

      IF (ltime) CALL get_timings (i_spread_ps_pre,ntstep, dt, izerror)

! spreading multi-level data
! --------------------------

      CALL mult_org_spread ( kk )
!     ====================

      IF (ltime) CALL get_timings (i_spread_multi_level,ntstep, dt, izerror)

! spreading upper-air single-level data
! -------------------------------------

      CALL upair_org_spread ( kk )
!     =====================

      IF (ltime) CALL get_timings (i_spread_upper_air,ntstep, dt, izerror)

! spreading surface-level data
! ----------------------------

      CALL surf_org_spread ( kk )
!     ====================

      IF (ltime) CALL get_timings (i_spread_surface,ntstep, dt, izerror)

    ENDIF


!-------------------------------------------------------------------------------
! Section 8: Computation of the analysis increments,
!            and update of the prognostic variables
!            (this includes communication between PE's only if (kk == 1))
!-------------------------------------------------------------------------------

! compute geostrophic pressure increments from (10-m) wind analysis increments
! ----------------------------------------------------------------------------

    IF ((mpsgcor >= 1) .AND. (lgetai) .AND. (kk == ke)) THEN

      CALL geostroph_ps_corr
!     ======================

! compute temperature correction for geostrophic pressure increments
! ------------------------------------------------------------------

      CALL ps_temperatur_corr ( 2 )
!     =======================

    ENDIF
!-------------------------------------------------------------------------------

! nudging of mass field, including humidity
! -----------------------------------------

    CALL nudge_humid_mass ( lgetai , lconai , kk )
!   =====================

! nudging of horizontal wind, includ. computation of geostrophic wind increments
! ------------------------------------------------------------------------------

    CALL nudge_horiz_wind ( lgetai , lconai , kk )
!   =====================

    IF (ltime) CALL get_timings (i_nudging_corr, ntstep, dt, izerror)

! _______________________________
  ENDDO loop_over_vertical_levels
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ENDIF   !   (ltnudge)


IF ((lwonl) .AND. (ntstep >= nstop)) THEN
  ! retrieve the required microphysics tracers (at nnew)
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnew, ptr = qv)
  IF (izerror /= 0_iintegers) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  DO kk = ke , 1 , -1
    IF ((kk == ke/2) .OR. (kk == ke) .OR. (kk == 1)) THEN
      PRINT*,'verifn: Tq ',ntstep,kk, t(ionl,jonl,kk,nnew),qv(ionl,jonl,kk)
      PRINT*,'verifn: pu ',ntstep,kk,pp(ionl,jonl,kk,nnew), u(ionl,jonl,kk,nnew)
    ENDIF
  ENDDO
ENDIF


! IF (ntstep >= nstop   ) THEN

!   ijtot  =  (ke/5+1) *(ie_tot/10+1) *(je_tot/10+1)
!   ALLOCATE( ibufloc (ijtot*3) )
!   ALLOCATE( rbufloc (ijtot*6) )
!   ibufloc = 0
!   rbufloc = 0.0_wp
!   ncount   = 0

!!  DO  kk = ke    , 1    , -1
!   DO  kk = ke    , 1    , -5
!     DO   j = jstart, jend
!       DO i = istart, iend
!!        IF (num_compute > 1) THEN
!           itot  =  i + isubpos(my_cart_id,1) - nboundlines - 1
!           jtot  =  j + isubpos(my_cart_id,2) - nboundlines - 1
!!        ELSE
!!          itot  =  i
!!          jtot  =  j
!!        ENDIF
!         IF ((MOD( itot,10 ) == 0) .AND. (MOD( jtot,10 ) == 0)) THEN
!!          PRINT *,'XYT', itot, jtot, kk,  t(i,j,kk,nnew), qv(i,j,kk)
!!          PRINT *,'XYP', itot, jtot, kk, pp(i,j,kk,nnew), qc(i,j,kk)
!!          PRINT *,'XYV', itot, jtot, kk,  u(i,j,kk,nnew),  v(i,j,kk,nnew)
!           ncount = ncount + 1
!           ibufloc (3*ncount-2) = itot
!           ibufloc (3*ncount-1) = jtot
!           ibufloc (3*ncount  ) = kk
!           rbufloc (6*ncount-5) =  t(i,j,kk,nnew)
!           rbufloc (6*ncount-4) = qv(i,j,kk)
!           rbufloc (6*ncount-3) = pp(i,j,kk,nnew)
!           rbufloc (6*ncount-2) = qc(i,j,kk)
!           rbufloc (6*ncount-1) =  u(i,j,kk,nnew)
!           rbufloc (6*ncount  ) =  v(i,j,kk,nnew)
!         ENDIF
!       ENDDO
!     ENDDO
!   ENDDO

!   CALL obs_gather_buffer ( ncount*3, ijtot*3, ibufloc, iicount, ibuftot, 0   &
!                          , num_compute, my_cart_id, icomm_cart, imp_integers )
! ! ======================

!   CALL obs_gather_buffer ( ncount*6, ijtot*6, rbufloc, ircount, rbuftot, 0   &
!                          , num_compute, my_cart_id, icomm_cart, imp_reals )
! ! ======================
!   DEALLOCATE( ibufloc )
!   DEALLOCATE( rbufloc )

!   IF (my_cart_id == 0) THEN
!     ncount = iicount / 3
!     PRINT *,'YYcount ', iicount, ircount, ncount
!     IF (ncount /= ircount/6) PRINT *,'CAUTION ', ncount, ircount, iicount
!     DO j = 1, ncount
!       itot = ibuftot(3*j-2)
!       jtot = ibuftot(3*j-1)
!       kk   = ibuftot(3*j  )
!       PRINT *,'XYT', itot, jtot, kk, rbuftot(6*j-5), rbuftot(6*j-4)
!       PRINT *,'XYP', itot, jtot, kk, rbuftot(6*j-3), rbuftot(6*j-2)
!       PRINT *,'XYV', itot, jtot, kk, rbuftot(6*j-1), rbuftot(6*j  )
!     ENDDO
!   ENDIF
!   DEALLOCATE ( ibuftot )
!   DEALLOCATE ( rbuftot )
! ENDIF

!-------------------------------------------------------------------------------
! Section 9: Cleanup  (without communication between PE's)
!-------------------------------------------------------------------------------

IF ((ltnudge) .OR. (ltverif)) THEN

  IF (lgetai)  CALL gather_spread_aux ( 3 )
!              ======================

  CALL gather_varia ( 5 )
! =================

  IF (ltime) CALL get_timings (i_spread_ps_pre, ntstep, dt, izerror)

ENDIF   !   (ltnudge) .OR. (ltverif)

IF ((ltnudge) .OR. (ltverif))  lfirst  = .FALSE.

IF ((ntstep == 0) .AND. (.NOT. l2tls)) THEN
  IF (my_cart_id == 0) PRINT '(5X,"STEP 0: dt,dt2,dtdeh: within the nudging "  &
                              &,2F5.0, F7.4)' , dt, dt2, dtdeh
  dt    = c05 * dt
  dt2   = c2 * dt
  dtdeh = dt / c3600
  IF (my_cart_id == 0) PRINT '(5X,"STEP 0: dt,dt2,dtdeh: outside the nudging"  &
                              &,2F5.0, F7.4)' , dt, dt2, dtdeh
ENDIF

!-------------------------------------------------------------------------------
! End of module procedure organize_nudging
!-------------------------------------------------------------------------------

END SUBROUTINE organize_nudging

!===============================================================================

ELEMENTAL REAL (KIND=wp) FUNCTION rmod  ( ztti, ziv )
  !----------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  ztti, ziv
  !  uses parameter:                           epsy
  !---------------------------------------------------------------------------
  ! MOD function for positive REALS
  !---------------------------------------------------------------------------
  !
  rmod  =  ztti  -  ziv * REAL(INT( ztti/ziv + epsy ),wp) + epsy
  !
END FUNCTION rmod

!-------------------------------------------------------------------------------

END MODULE src_obs_use_org
