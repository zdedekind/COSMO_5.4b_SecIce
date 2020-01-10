!+ Subroutine (called from lmorg) to organize the model initialization
!------------------------------------------------------------------------------

SUBROUTINE dfi_initialization (lbd_frame, undef, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
! A dynamic initialization of the inital data is performed using the
! digital filtering method.
!
! Method:
!  After Peter Lynch, 1996.
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  Ulrich.Schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.4        1998/05/22 Guenther Doms
!  Initial release
! 1.7        1998/07/16 Guenther Doms
!  Change in the calling list of routine 'calps'.
! 1.8        1998/08/03 Ulrich Schaettler
!  Dependency on module data_io added because of parameterlist for org_radiation.
! 1.9        1998/09/16 Guenther Doms
!  Include USE of varibale 'l2tls' and routine 'vertical_diffusion'.
! 1.20       1999/01/07 Guenther Doms
!  Renaming of some global variables.
! 1.29       1999/05/11 Ulrich Schaettler
!  Changed USE-list for data_parallel and utility-modules and excluded 
!  organize_physics
! 1.30       1999/06/24 Matthias Raschendofer
!  Additional Initialization of tketens and ntke
! 1.34       1999/12/10 Ulrich Schaettler
!  Changed timings and CALLs to utility- and organizing routines
! 1.39       2000/05/03 Ulrich Schaettler
!  Included data from former module data_filter. 
!  Eliminated direct reference to module src_relaxation.
!  Read in Namelist Input and include subroutine input_inictl.
!  Adapted Calls to timing routines.
! 2.8        2001/07/06 Ulrich Schaettler
!  Eliminated non-necessary variables and dependencies from the USE-lists
! 2.11       2001/09/28 Ulrich Schaettler
!  Introduced filtering of cloud ice qi
! 2.13       2002/01/18 Ulrich Schaettler
!  Bug correction in data exchange for the cloud ice scheme
! 2.18       2002/07/16 Reinhold Schrodin
!  Eliminated variable rhde, use cf_snow instead
!  Use lmulti_layer to define boundary values for old soil model correctly
! 3.5        2003/09/02 Ulrich Schaettler
!  Corrections for ndfi=0/1; Corrections to boundary exchange for variables
!  from convection; Adaptation of interface for exchg_boundaries.
! 3.6        2003/12/11 Ulrich Schaettler
!  Limit the boundary values for the humidity variables to MIN=0.0;
!  Introduced use of frames for the boundary variables
! 3.7        2004/02/18 Ulrich Schaettler
!  Replaced type_gscp by lprog_qi; Adaptations for 2 timelevel scheme
! 3.15       2005/03/03 Ulrich Schaettler
!  Replaced FLOAT by REAL
! 3.18       2006/03/03 Ulrich Schaettler
!  Reset the sunshine-hours after initialization
! 3.21       2006/12/04 Burkhardt Rockel
!  Renamed variable sunshhrs to dursun
! V3_23        2007/03/30 Lucio Torrisi
!  Bug fix (for re-constructing boundary fields after backward step)
!  Adaptation to Runge-Kutta scheme
! V4_1         2007/12/04 Lucio Torrisi
!  Initialize t_so, w_so in case of multi layer soil model
! V4_4         2008/07/16 Ulrich Schaettler, Lucio Torrisi
!  Eliminated timing variables which are unused
!  Adapted interface of get_timings
!  Also initialize rho_snow in case of multi layer soil model
! V4_9         2009/07/16 Ulrich Schaettler
!  Eliminate reduction of qv in case lana_qi, llb_qi = .TRUE.
! V4_11        2009/11/30 Ekaterina Machulskaya, Lucio Torrisi
!  Modifications to run the DFI initialization with multi-layer snow model (EM)
!  Initialization of alhfl_bs, alhfl_pl (LT)
! V4_12        2010/05/11 Michael Baldauf
!  Corrected dimension of variable qvsflx
!  Adapted interface parameter t_s for SR tgcom to account for seaice (JPS)
!  Commented unused statement function fpvsi (to avoid compiler warnings; US)
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_17        2011/02/24 Ulrich Blahak
!  Adapted interface of exchg_boundaries; corrected kzdims(1:20) -> kzdims(1:24);
!  eliminated my_peri_neigh
! V4_18        2011/05/26 Ulrich Schaettler
!  In case of ndfi=1:
!   Re-initialize start events for meanvalues, grid points and synsat pictures
!   Switched on soil processes and re-initialize soil fields after dfi with
!    mean values between start and end of dfi
!  Changed the code owner
! V4_20        2011/08/31 Ulrich Schaettler
!  Limited values for ndfi to 1/2.
!  In case of l2tls and ndfi==2: Abort. The Runge-Kutta scheme cannot be 
!  integrated backwards
! V4_23        2012/05/10 Ulrich Schaettler, Oli Fuhrer, CLM
!  Removed switch lprogprec: switching off transport of qr,qs,qg in adiabatic
!   (Leapfrog) backward integration will be possible again with tracer module
!  Removed src_2timelevel and related exchange routine
!  Removed SR exchange_l2dim, because this is now done within exchg_boundaries
!  Removed field qvt_diff from the COSMO-Model (Oli Fuhrer)
!  Put computation of next SynSat output into IF (luse_rttov) (Oli Fuhrer)
!  Added new field snow_melt and nullify it during Cleanup (CLM)
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer, Hans-Juergen Panitz
!  Replaced all global qx-variables by local ones and get their values from 
!  the corresponding tracers
!  Introduced nexch_tag for MPI boundary exchange tag to replace ntstep
! V4_26        2012/12/06 Ulrich Schaettler
!  Adapted variable names of multi-layer snow model to corresponding
!   short names for I/O
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  MESSy interface introduced
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Eliminated MESSy interface, because MESSy cannot be used with DFI
! V5_1         2014-11-28 Ulrich Blahak, Oliver Fuhrer
!  Changed the format of some YUSPECIF entries for the CLM namelist tool.
!  Implemented F2003 IOMSG-mechanism for better namelist error messages. (UB)
!  Replaced ireals by wp (working precision) (OF)
! V5_2         2015-05-21 Ulrich Schaettler
!  Reorganized numbering of MPI datatypes ID; Communication for humidity tracers 
!  cannot be done with MPI datatypes any more, because the corresponding pointers 
!  are local for this SR and not in long term storage
! V5_4         2016-03-10 Lucio Torrisi
!  Better treatment of clouds and precipitation during diabatic DFI step
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Modules used:

!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 2. physical constants and related variables
! -------------------------------------------
    r_d,          & ! gas constant for dry air
    r_v,          & ! gas constant for water vapor
    rdv,          & ! r_d / r_v
    o_m_rdv,      & ! 1 - r_d/r_v
    rvd_m_o,      & ! r_v/r_d - 1
    g,            & ! acceleration due to gravity

! 3. constants for parametrizations
! ---------------------------------
    b1,           & ! variables for computing the saturation vapour pressure
    b2w,          & ! over water (w) and ice (i)
    b2i,          & !               -- " --
    b3,           & !               -- " --
    b4w,          & !               -- " --
    b4i,          & !               -- " --
    aks2,         & ! variable for horizontal diffusion of second order
    aks4            ! variable for horizontal diffusion of fourth order

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------
    rho0       ,    & ! reference density at the full model levels    (kg/m3)
    dp0        ,    & ! reference pressure thickness of layers        ( Pa)
    p0         ,    & ! reference pressure at main levels             ( Pa)

! 2. external parameter fields                                        (unit)
! ----------------------------
    llandmask  ,    & ! landpoint mask

! 3. prognostic variables                                             (unit)
! -----------------------
    u          ,    & ! zonal wind speed                              ( m/s )
    v          ,    & ! meridional wind speed                         ( m/s )
    w          ,    & ! vertical wind speed (defined on half levels)  ( m/s )
    t          ,    & ! temperature                                   (  k  )
    pp                ! deviation from the reference pressure         ( pa  )

USE data_fields     , ONLY :   &

! 4. tendency fields for the prognostic variables                     (unit )
! -----------------------------------------------
    utens        ,  & ! u-tendency without sound-wave terms           ( m/s2)
    vtens        ,  & ! v-tendency without sound-wave terms           ( m/s2)
    wtens        ,  & ! w-tendency without sound-wave terms           ( m/s2
                      ! (defined on half levels )
    ttens        ,  & ! t-tendency without sound-wave terms           ( m/s2)
    pptens       ,  & ! pp-tendency without sound-wave terms          ( m/s2)
    tketens      ,  & ! tke-tendency (defined on half levels)         ( m/s2)

! 5. fields for surface values and soil model variables               (unit )
! -----------------------------------------------------
    ps        ,     & ! surface pressure                              ( pa  )
    t_snow    ,     & ! temperature of the snow-surface               (  k  )
    t_snow_mult,    & ! temperature of the snow-surface               (  k  )
    dzh_snow_mult,  & !
    w_snow_mult,    & ! total water content of snow                   (m H2O)
    wliq_snow ,     & ! liquid water content of snow                  (m H2O)
    t_s       ,     & ! temperature of the ground surface             (  k  )
    t_g       ,     & ! weighted surface temperature                  (  k  )
    qv_s      ,     & ! specific water vapor content on the surface   (kg/kg)
    t_m       ,     & ! temperature between upper and medium layer    (  K  )
    t_so      ,     & ! multi-layer soil temperature                  (  K  )
    w_so      ,     & ! multi-layer soil moisture                     (m H2O)
    w_snow    ,     & ! water content of snow                         (m H2O)
    rho_snow  ,     & ! snow density                                  ()
    freshsnow ,     & ! indicator for age of snow in top of snow layer(  -  )
    w_i       ,     & ! water content of interception water           (m H2O)
    w_g1      ,     & ! water content of the upper soil layer         (m H2O)
    w_g2      ,     & ! water content of the medium soil layer        (m H2O)
    w_g3      ,     & ! water content of the lower soil layer         (m H2O)
    w_cl      ,     & ! water content of the climatological layer     (m H2O)

!   fields for prognostic variables of the lake model FLake or ocean
!   variables
    t_ice      ,    & ! temperature of ice/water surface              (  K  )
    h_ice      ,    & ! lake/sea ice thickness                        (  m  )

! prognostic variables of the lake model                              (unit )
! -----------------------------------------------------
    t_mnw_lk   ,    & ! mean temperature of the water column          (  K  )
    t_wml_lk   ,    & ! mixed-layer temperature                       (  K  )
    t_bot_lk   ,    & ! temperature at the water-bottom sediment
                      ! interface                                     (  K  )
    t_b1_lk    ,    & ! temperature at the bottom of the upper layer
                      ! of the sediments                              (  K  )
    c_t_lk     ,    & ! shape factor with respect to the
                      ! temperature profile in lake thermocline       (  -  )
    h_ml_lk    ,    & ! thickness of the mixed-layer                  (  m  )
    h_b1_lk           ! thickness of the upper layer
                      ! of bottom sediments                           (  m  )

USE data_fields     , ONLY :   &

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------
    rho        ,    & ! density of moist air
    qrs        ,    & ! precipitation water (water loading)           (kg/kg)
    dqvdt      ,    & ! threedimensional moisture convergence         ( 1/s )
    tke        ,    & ! SQRT(2 * turbulent kinetik energy)            ( m/s )
    qvsflx     ,    & ! surface flux of water vapour                  (kg/m2s)
    aumfl_s    ,    & ! average u-momentum flux (surface)             ( N/m2)
    avmfl_s    ,    & ! average v-momentum flux (surface)             ( N/m2)
    ashfl_s    ,    & ! average sensible heat flux (surface)          ( W/m2)
    alhfl_s    ,    & ! average latent heat flux (surface)            ( W/m2)

! 7. fields for model output and diagnostics                          (unit )
! ------------------------------------------
    rain_gsp   ,    & ! amount of rain from grid-scale precip. (sum)  (kg/m2)
    snow_gsp   ,    & ! amount of snow from grid-scale precip. (sum)  (kg/m2)
    rain_con   ,    & ! amount of rain from convective precip. (sum)  (kg/m2)
    snow_con   ,    & ! amount of snow from convective precip. (sum)  (kg/m2)
    snow_melt  ,    & ! amount of snow melt (sum)                     (kg/m2)
    asob_s     ,    & ! average solar radiation budget (surface)      ( W/m2)
    athb_s     ,    & ! average thermal radiation budget (surface)    ( W/m2)
    apab_s     ,    & ! average photosynthetic active radiation (sfc) ( W/m2)
    asob_t     ,    & ! average solar radiation budget (model top)    ( W/m2)
    athb_t     ,    & ! average thermal radiation budget (model top)  ( W/m2)
    dursun     ,    & ! sunshine duration                             (  s  )
    alhfl_bs   ,    & ! average latent heat flux from bare soil evap. ( W/m2)
    alhfl_pl   ,    & ! average latent heat flux from plants          ( W/m2)

! 8. fields for the boundary values                                   (unit )
! ---------------------------------
    u_bd          , & ! boundary field for u                          ( m/s )
    v_bd          , & ! boundary field for v                          ( m/s )
    t_bd          , & ! boundary field for t                          (  k  )
    w_bd          , & ! boundary field for w                          ( m/s )
    pp_bd             ! boundary field for pp                         (  pa )

!------------------------------------------------------------------------------

USE data_io,          ONLY :   &
    lana_rho_snow   ! if .TRUE., take rho_snow-values from analysis file
                    ! else, it is set in the model

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------
    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction
    ke1,          & ! ke+1
    ke_soil,      & ! number of layers in multi-layer soil model

! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
    istartpar,    & ! start index for computations in the parallel program
    iendpar,      & ! end index for computations in the parallel program
    jstartpar,    & ! start index for computations in the parallel program
    jendpar,      & ! end index for computations in the parallel program

! 4. variables for the time discretization and related variables
! --------------------------------------------------------------
    dt,           & ! long time-step
    ed2dt,        & ! 1 / (2 * dt)
    dt2,          & ! 2 * dt
    dtdeh,        & ! dt / 3600 seconds
    nehddt,       & ! 3600 s / dt
    istart, iend, & !
    jstart, jend, & !
    epsass,       & ! eps for the Asselin-filter

! 8. Organizational variables to handle the COSMO humidity tracers
! ----------------------------------------------------------------
    idt_qv, idt_qc, idt_qi, idt_qr, idt_qs, idt_qg

!------------------------------------------------------------------------------

USE data_parallel,      ONLY :  &
    num_compute,     & ! number of compute PEs
    my_cart_neigh,   & ! neighbors of this subdomain in the cartesian grid
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    icomm_cart,      & ! communicator for the virtual cartesian topology
    iexch_req,       & ! stores the sends requests for the neighbor-exchange
                       ! that can be used by MPI_WAIT to identify the send
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    imp_integers,    & ! determines the correct INTEGER type used in the
                       ! model for MPI
    imp_logical,     & ! determines the correct LOGICAL type used in the
                       ! model for MPI
    nexch_tag,       & ! tag to be used for MPI boundary exchange
                       !  (in calls to exchg_boundaries)
    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains
    ldatatypes,      & ! if .TRUE.: use MPI-Datatypes for some communications
    ltime_barrier,   & ! if .TRUE.: use additional barriers for determining the
                       ! load-imbalance
    ncomm_type,      & ! type of communication
    sendbuf,         & ! sending buffer for boundary exchange:
                       ! 1-4 are used for sending, 5-8 are used for receiving
    isendbuflen,     & ! length of one column of sendbuf
    nproc, realbuf, intbuf, logbuf

!------------------------------------------------------------------------------

USE data_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    nstart,       & ! first time step of the forecast
    nstop,        & ! last time step of the forecast
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
    nold,         & ! corresponds to ntstep - 1
    nnow,         & ! corresponds to ntstep
    nnew,         & ! corresponds to ntstep + 1
    ntke,         & ! timestep for tke

! 2. boundary definition and update
! ---------------------------------
    nlastbound,   & ! time step of the last boundary update
    nincbound,    & ! time step increment of boundary update
    nbd1,         & ! indices for permutation of the
    nbd2            ! two boundary time levels

USE data_runcontrol , ONLY :   &

! 3. controlling the physics
! --------------------------
    itype_gscp,   & ! type of grid-scale precipitaiton physics
    nlgw,         & ! number of prognostic soil water levels
    lphys,        & ! forecast with physical parametrizations
    lgsp,         & ! forecast with grid scale precipitation
    lsoil,        & ! forecast with soil model
    ltur,         & ! forecast with vertical diffusion
    itype_turb,   & ! type of turbulent diffusion parametrization
    lmulti_layer, & ! run multi-layer soil model
    lmulti_snow,  & ! run multi-layer snow model
    llake,        & ! forecast with lake model
    lseaice,      & ! forecast with sea ice model
    lprog_tke,    & ! prognostic treatment of TKE (for itype_turb=5/7)
    nincconv,     & ! time step increment for running the convection scheme

! 5. additional control variables
! -------------------------------
    l2tls,        & ! forecast with 2-TL integration scheme
    lcond,        & ! forecast with condensation/evaporation
    ltime,        & ! detailled timings of the program are given
    lconv,        & ! forecast with convection
    l2dim,        & !
    lw_freeslip,  & !
    lperi_x,      & ! if lgen=.TRUE.: periodic boundary conditions (.TRUE.) in x-dir.
                    ! or with Davies conditions (.FALSE.)
    lperi_y,      & ! if lgen=.TRUE.: periodic boundary conditions (.TRUE.) in y-dir.
                    ! or with Davies conditions (.FALSE.)
    l2dim,        & ! 2 dimensional runs
    luse_rttov,   & !
    nbl_exchg,    & !
    lconf_avg,    & ! average convective forcings in case of massflux closure
    idbg_level,   & !
    lprintdeb_all,& !
    nuspecif        ! unit number for file YUSPECIF

!------------------------------------------------------------------------------

USE data_soil       , ONLY :   &
    lsoilinit_dfi,& ! initialize soil after dfi forward launching
    cf_snow

!------------------------------------------------------------------------------

#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
USE data_satellites , ONLY :  &
    nsat_steps, nsat_next
#endif

!------------------------------------------------------------------------------

USE src_meanvalues,     ONLY:  n0meanval, nextmeanval
USE src_gridpoints,     ONLY:  n0gp, nnextgp

!------------------------------------------------------------------------------
! Imported utilities
!------------------------------------------------------------------------------

USE utilities       ,   ONLY :  dolph
USE meteo_utilities ,   ONLY :  calps, tgcom, calrho
USE environment     ,   ONLY :  exchg_boundaries, comm_barrier, model_abort
USE time_utilities  ,   ONLY :  get_timings, i_add_computations,              &
                        i_dyn_computations, i_phy_computations, i_relaxation, &
                        i_barrier_waiting_dyn, i_communications_dyn
USE parallel_utilities, ONLY :  distribute_values

USE data_tracer     ,   ONLY :  T_ERR_NOTFOUND, T_LBC_ID, T_LBC_FILE,         &
                                T_INI_ID, T_INI_FILE
USE src_tracer      ,   ONLY :  trcr_get, trcr_errorstr, trcr_meta_get

!==============================================================================
 
IMPLICIT NONE
 
!==============================================================================

! Variables for the digital filtering
  REAL    (KIND=wp)          ::  &
    tspan, & ! Time-span (in seconds) for the adiabatic and diabatic
             ! stages of the initialization
    taus,  & ! Cuttoff period (in seconds) for the filter
             ! (Taus is the stop-band edge for the Dolph filter)
    dtbak, & ! Time-step for the backcast filtering stage (sec)
    dtfwd, & ! Time-step for the forecast filtering stage (sec)
    thsbak,& ! Cuttoff frequency for the hindcast stage (Hz)
    thsfwd   ! Cuttoff frequency for the forecast stage (Hz)

  INTEGER (KIND=iintegers)   ::  &
    ndfi,  & ! indicator for kind of filtering
             ! = 0 : No filtering (corresponds to LDFI=.F.)
             ! = 1 : Launching (forward stage only)
             ! = 2 : Full two stage filtering (default)
    nfilt, & ! indicator for method of filtering
             ! = 1 : Dolph-Chebyshev filter (default)
             !       (no other filter is implemented yet)
    nbak,  & ! Number of time-steps for the adiabatic backcast
    nfwd     ! Number of time-steps for the  diabatic forecast

  INTEGER (KIND=iintegers)   ::  &
    itype_dfi_qx, &  ! treatment of qx variables in dfi (active only for ndfi=1)
                     ! = 0 : all qx variables are filtered in forward stage
                     ! = 1 : qv filtered, other variables initialised with
                     !       istantaneous values valid at nhalf step of forward stage
                     ! = 2 : qv filtered but correct to mantain saturated points,
                     !       other variables initialised with
                     !       istantaneous values valid at nhalf step of forward stage
    itype_dfi_soil   ! treatment of soil variables in dfi (active only for ndfi=1)
                     ! = 0 : all soil variables are derived from the average of
                     !       values at the initial and final step of forward stage
                     ! = 1 : qv filtered, other variables initialised with
                     !       istantaneous values valid at nhalf step of forward stage

! Global (supporting) arrays
  REAL    (KIND=wp),     ALLOCATABLE      ::  &
    hdfi(:)          ! Filter coefficients of DFI-scheme, which are calculated
                     ! in SUBROUTINE dolph

!==============================================================================

! Parameter list:
! ---------------

LOGICAL,                  INTENT(IN)    ::         &
  lbd_frame       ! if .TRUE., boundary data are on a frame

REAL (KIND=wp),           INTENT(IN)     ::        &
  undef           ! value for undefined grid points

INTEGER (KIND=iintegers), INTENT(OUT)    ::        &
  ierror       ! error status

CHARACTER (LEN=  *),      INTENT(OUT)    ::        &
  yerrmsg      ! error message

! Local scalars:
! -------------

INTEGER (KIND=iintegers)   ::        &
  nsp, i,j,k, nhalf, nmax, nx, & ! support variable
  nzntstep, nzincbound,        & ! to store the original settings
  iztype_gscp, nzstop,         & ! to store the original settings
  nzlastbound, istat, izerrstat, nuin, kzdims(24), nx0, izdebug

INTEGER (KIND=iintegers)   ::        &
  izlbc_qi,                    & ! lateral BC type for QI
  izlbc_qr,                    & ! lateral BC type for QR
  izlbc_qs,                    & ! lateral BC type for QS
  izlbc_qg,                    & ! lateral BC type for QG
  izic_qi                        ! IC type for QI

REAL (KIND=wp)             ::        &
  z1, z2, zst, zsge, zsp,      & ! factors for time interpolation
  fpvsw,        fqvs,          & ! for the statement functions
!        fpvsi,                & ! for the statement functions
  zqvt,                        & ! for computing qv_s
  zaks2, zaks4, zdt, zepsass     ! to store the original settings

LOGICAL                    ::        &
  lzphys, lzcond, lzsoil,      & ! to store the original settings
  lzlake, lzseaice,            & ! to store the original settings
  lzconv                         ! to determine the communication

CHARACTER (LEN= 9)         ::   yinput

! Local Arrays: 
! -------------

REAL (KIND=wp)             ::        &
  hfu (ie,je,ke), & ! Summation array for weighted filtered values of U
  hfv (ie,je,ke), & ! Summation array for weighted filtered values of V
  hfw (ie,je,ke1),& ! Summation array for weighted filtered values of W
  hft (ie,je,ke), & ! Summation array for weighted filtered values of T
  hfrh(ie,je,ke), & ! Summation array for weighted filtered values of RH
  hfqv(ie,je,ke), & ! Summation array for weighted filtered values of QV
  hfqc(ie,je,ke), & ! Summation array for weighted filtered values of QC
  hfqi(ie,je,ke), & ! Summation array for weighted filtered values of QI
  hfqr(ie,je,ke), & ! Summation array for weighted filtered values of QR
  hfqs(ie,je,ke), & ! Summation array for weighted filtered values of QS
  hfqg(ie,je,ke), & ! Summation array for weighted filtered values of QG
  hfpp(ie,je,ke), & ! Summation array for weighted filtered values of PP
  qvb (ie,je,ke), & ! for saving the original boundary values
  qcb (ie,je,ke), & !
  qrb (ie,je,ke), & !
  qsb (ie,je,ke), & !
  qgb (ie,je,ke), & !
  qib (ie,je,ke)    !

REAL (KIND=wp),     ALLOCATABLE         ::        &
  zrh      (:,:,:),   & !
  zqc      (:,:,:),   & !
  zqv      (:,:,:),   & !
  zqi      (:,:,:),   & !
  zqr      (:,:,:),   & !
  zqs      (:,:,:),   & !
  zqg      (:,:,:)      !

! fields to save the initial values for soil variables
! for later averaging, if only forward launching is done
REAL (KIND=wp),     ALLOCATABLE         ::        &
  zts      (:,:),     & !
  zt_m     (:,:),     & !
  zw_g1    (:,:),     & !
  zw_g2    (:,:),     & !
  zw_g3    (:,:),     & !
  zt_so    (:,:,:),   & !
  zw_so    (:,:,:),   & !
  zfreshsnw(:,:),     & !
  zrhosnw  (:,:),     & !
  zqv_s    (:,:),     & !
  zw_snow  (:,:),     & !
  zt_snow  (:,:),     & !
  ztmnw_lk (:,:),     & !
  ztwml_lk (:,:),     & !
  ztbot_lk (:,:),     & !
  zct_lk   (:,:),     & !
  zhml_lk  (:,:),     & !
  ztice    (:,:),     & !
  zhice    (:,:)        !

REAL (KIND=wp)             ::        &
  zt_s(ie,je)       ! = t_s   on land and sea
                    ! = t_ice on sea ice (if present)

REAL (KIND=wp),     ALLOCATABLE         ::        &
  cheby(:), time(:), time2(:), weight(:), weight2(:)

! Local pointers:
!----------------
REAL (KIND=wp),     POINTER :: &
  qv_now (:,:,:)  => NULL(),      &  ! QV at nnow
  qv_new (:,:,:)  => NULL(),      &  ! QV at nnew
  qv_old (:,:,:)  => NULL(),      &  ! QV at nold
  qv_nx  (:,:,:)  => NULL(),      &  ! QV at nx0
  qv_bd  (:,:,:,:)=> NULL(),      &  ! QV boundaries
  qv_tens(:,:,:)  => NULL(),      &  ! QV tendency
  qc_now (:,:,:)  => NULL(),      &  ! QC at nnow
  qc_new (:,:,:)  => NULL(),      &  ! QC at nnew
  qc_old (:,:,:)  => NULL(),      &  ! QC at nold
  qc_nx  (:,:,:)  => NULL(),      &  ! QC at nx0
  qc_bd  (:,:,:,:)=> NULL(),      &  ! QC boundaries
  qc_tens(:,:,:)  => NULL(),      &  ! QC tendency
  qi_now (:,:,:)  => NULL(),      &  ! QI at nnow
  qi_new (:,:,:)  => NULL(),      &  ! QI at nnew
  qi_old (:,:,:)  => NULL(),      &  ! QI at nold
  qi_nx  (:,:,:)  => NULL(),      &  ! QI at nx0
  qi_bd  (:,:,:,:)=> NULL(),      &  ! QI boundaries
  qi_tens(:,:,:)  => NULL(),      &  ! QI tendency
  qr_now (:,:,:)  => NULL(),      &  ! QR at nnow
  qr_new (:,:,:)  => NULL(),      &  ! QR at nnew
  qr_old (:,:,:)  => NULL(),      &  ! QR at nold
  qr_nx  (:,:,:)  => NULL(),      &  ! QR at nx0
  qr_bd  (:,:,:,:)=> NULL(),      &  ! QR boundaries
  qs_now (:,:,:)  => NULL(),      &  ! QS at nnow
  qs_new (:,:,:)  => NULL(),      &  ! QS at nnew
  qs_old (:,:,:)  => NULL(),      &  ! QS at nold
  qs_nx  (:,:,:)  => NULL(),      &  ! QS at nx0
  qs_bd  (:,:,:,:)=> NULL(),      &  ! QS boundaries
  qg_now (:,:,:)  => NULL(),      &  ! QG at nnow
  qg_new (:,:,:)  => NULL(),      &  ! QG at nnew
  qg_old (:,:,:)  => NULL(),      &  ! QG at nold
  qg_nx  (:,:,:)  => NULL(),      &  ! QG at nx
  qg_bd  (:,:,:,:)=> NULL()          ! QG boundaries

!US  qr_tens(:,:,:)  => NULL(),      &  ! QR tendency
!US  qs_tens(:,:,:)  => NULL(),      &  ! QS tendency
!US  qg_tens(:,:,:)  => NULL()          ! QG tendency

! For error handling
! ------------------
  INTEGER (KIND=iintegers) :: izerror

  CHARACTER (LEN=255)      :: yzerrmsg
  CHARACTER (LEN=25)       :: yzroutine

!- End of header
!==============================================================================
 
!------------------------------------------------------------------------------
! Statement Functions
!------------------------------------------------------------------------------
 
  fpvsw(zst)       = b1 * EXP( b2w * (zst-b3)/(zst-b4w) )
! fpvsi(zst)       = b1 * EXP( b2i * (zst-b3)/(zst-b4i) )
  fqvs (zsge, zsp) = rdv * zsge / ( zsp - o_m_rdv * zsge )

!------------------------------------------------------------------------------
! Begin Subroutine dfi_initialization
!------------------------------------------------------------------------------

  izerror   = 0_iintegers
  yzerrmsg  = '   '
  ierror    = 0_iintegers
  yerrmsg   = '   '
  yzroutine = 'dfi_initialization'
  izerrstat = 0_iintegers
  kzdims(:) = 0_iintegers

  ! Initialize, whether debug output shall be done
  IF (lprintdeb_all) THEN
    izdebug = idbg_level
  ELSE
    IF (my_cart_id == 0) THEN
      izdebug = idbg_level
    ELSE
      izdebug = 0
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 1.1: Read in NAMELIST input
!------------------------------------------------------------------------------

  ! Open NAMELIST-INPUT file
  IF (my_cart_id == 0) THEN
    PRINT *,'    INPUT OF THE NAMELISTS FOR INITIALIZATION'
    yinput   = 'INPUT_INI'
    nuin     =  1

    OPEN(nuin   , FILE=yinput  , FORM=  'FORMATTED', STATUS='UNKNOWN',  &
         IOSTAT=izerrstat)
    IF(izerrstat /= 0) THEN
      yerrmsg  = ' ERROR    *** Error while opening file INPUT_INI *** '
      ierror   = 2
      RETURN
    ENDIF
  ENDIF

  ! read the NAMELIST group inictl:
  CALL input_inictl (nuspecif, nuin, izerrstat)

  IF (my_cart_id == 0) THEN
    ! Close file for input of the NAMELISTS
    CLOSE (nuin    , STATUS='KEEP')
  ENDIF

  IF (izerrstat < 0) THEN
    yerrmsg = 'ERROR *** while reading NAMELIST Group /INICTL/ ***'
    ierror  = 3
    RETURN
  ELSEIF (izerrstat > 0) THEN
    yerrmsg = 'ERROR *** Wrong values occured in NAMELIST INPUT_DYN ***'
    ierror  = 4
    RETURN
  ENDIF

!------------------------------------------------------------------------------
! Section 1.2: Initializations
!------------------------------------------------------------------------------
  
  ! store original settings
  nzntstep    = ntstep
  nzstop      = nstop
  nzincbound  = nincbound
  nzlastbound = nlastbound
  iztype_gscp = itype_gscp
  zdt         = dt
  zepsass     = epsass
  zaks2       = aks2
  zaks4       = aks4
  lzphys      = lphys
  lzcond      = lcond
  lzsoil      = lsoil
  lzlake      = llake
  lzseaice    = lseaice
 
  ! initialize new time level for the surface fields: these are held 
  ! constant during the backward and the forward integration
  IF (.NOT. l2tls) THEN
    nold    = 1
    nnow    = 2
    nnew    = 3
    nx0     = nold
  ELSE
    nnow    = 2
    nnew    = 1
    nx0     = nnow
  ENDIF

  CALL update_trcr_pointers()

  ! Retrieve the metadata for the tracers
  !------------------------------------------
  IF (ASSOCIATED(qi_now)) THEN
    CALL trcr_meta_get(izerror, idt_qi, T_LBC_ID, izlbc_qi)
    IF (izerror /= 0_iintegers) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_meta_get(izerror, idt_qi, T_INI_ID, izic_qi)
    IF (izerror /= 0_iintegers) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
  ENDIF
  IF (ASSOCIATED(qr_now)) THEN
    CALL trcr_meta_get(izerror, idt_qr, T_LBC_ID, izlbc_qr)
    IF (izerror /= 0_iintegers) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
  ENDIF
  IF (ASSOCIATED(qs_now)) THEN
    CALL trcr_meta_get(izerror, idt_qs, T_LBC_ID, izlbc_qs)
    IF (izerror /= 0_iintegers) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
  ENDIF
  IF (ASSOCIATED(qg_now)) THEN
    CALL trcr_meta_get(izerror, idt_qg, T_LBC_ID, izlbc_qg)
    IF (izerror /= 0_iintegers) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
  ENDIF

  ! Initializations for the surface variables
  DO j = 1,je
    DO i = 1,ie
      IF ( llandmask(i,j) .EQV. .TRUE. ) THEN
        t_snow(i,j,nnew) = t_snow(i,j,nx0)
        t_s   (i,j,nnew) = t_s   (i,j,nx0)
        qv_s  (i,j,nnew) = qv_s  (i,j,nx0)
        w_snow(i,j,nnew) = w_snow(i,j,nx0)
      ELSE
        t_snow(i,j,nnew) = t_snow(i,j,nx0)
        t_s   (i,j,nnew) = t_s   (i,j,nx0)
        zqvt             = fqvs ( fpvsw ( t_s(i,j,nx0) ) , ps (i,j,nx0) )
        qv_s  (i,j,:)    = zqvt
      ENDIF
    ENDDO
  ENDDO

  ! Initializations for the soil variables
  IF (.NOT. lmulti_layer) THEN
    DO j = 1,je
      DO i = 1,ie
        IF ( llandmask(i,j) .EQV. .TRUE. ) THEN
          t_m   (i,j,nnew) = t_m   (i,j,nx0)
          w_g1  (i,j,nnew) = w_g1  (i,j,nx0)
          w_g2  (i,j,nnew) = w_g2  (i,j,nx0)
          IF ( nlgw == 3 ) THEN
            w_g3(i,j,nnew) = w_g3(i,j,nx0)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ELSE
    DO j = 1,je
      DO i = 1,ie
        IF ( llandmask(i,j) .EQV. .TRUE. ) THEN
          t_so(i,j,:,nnew)   = t_so(i,j,:,nx0)
          w_so(i,j,:,nnew)   = w_so(i,j,:,nx0)
          rho_snow(i,j,nnew) = rho_snow(i,j,nx0)
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  ! compute interface variable zt_s depending on actual settings
  DO j = 1,je
    DO i = 1,ie
      zt_s(i,j) = t_s(i,j,nnew)
      IF (lseaice) THEN
        IF (.NOT. llandmask(i,j) .AND. h_ice(i,j,nx0) > 0.0_wp) THEN
          zt_s(i,j) = t_ice(i,j,nx0)
        ENDIF
      ENDIF
    ENDDO
  ENDDO

  IF (lmulti_snow) THEN
    ! Initializations for the multi-layer snow-model
    DO j = 1,je
      DO i = 1,ie
        IF ( llandmask(i,j) .EQV. .TRUE. ) THEN
          t_snow_mult   (i,j,:,nnew) = t_snow_mult   (i,j,:,nx0)
          w_snow_mult   (i,j,:,nnew) = w_snow_mult   (i,j,:,nx0)
          wliq_snow     (i,j,:,nnew) = wliq_snow     (i,j,:,nx0)
          dzh_snow_mult (i,j,:,nnew) = dzh_snow_mult (i,j,:,nx0)
        ELSE
          t_snow_mult   (i,j,:,nnew) = t_snow_mult   (i,j,:,nx0)
        ENDIF
      ENDDO
    ENDDO

    CALL tgcom ( t_g(:,:,nnew), t_snow_mult(:,:,1,nnew), &
                 zt_s(:,:)    , w_snow(:,:,nnew), &
                 llandmask(:,:) , ie, je, cf_snow,&
                 istartpar, iendpar, jstartpar, jendpar )
  ELSE
    CALL tgcom ( t_g(:,:,nnew), t_snow(:,:,nnew), &
                 zt_s(:,:)    , w_snow(:,:,nnew), &
                 llandmask(:,:) , ie, je, cf_snow,&
                 istartpar, iendpar, jstartpar, jendpar )
  ENDIF
  
  ! Allocate array hdfi
  nmax = MAX (nbak, nfwd)
  ALLOCATE (hdfi   (0:nmax), STAT=istat)
  ALLOCATE (cheby  (0:nmax), STAT=istat)
  ALLOCATE (time   (0:nmax), STAT=istat)
  ALLOCATE (time2  (0:nmax), STAT=istat)
  ALLOCATE (weight (0:nmax), STAT=istat)
  ALLOCATE (weight2(0:nmax), STAT=istat)

!------------------------------------------------------------------------------
! Section 2: Adiabatic backward integration
!------------------------------------------------------------------------------

  ! Adiabatic backward integration is only possible with Leapfrog Scheme.
  ! For Runge-Kutta (l2tls=.TRUE.), only ndfi=1 is accepted
  IF (ndfi == 2) THEN
    ! settings for backward integration
    IF ( .NOT.l2tls ) THEN
      nold = 1
      nnow = 2
      nnew = 3
      dt   = -dtbak * 0.5_wp   ! only in the first step
!   is not possible for ndfi=2
!   ELSE
!     nnow = 2
!     nnew = 1
!     dt   = -dtbak
!     nstop = nbak - 1
    END IF

    CALL update_trcr_pointers()

    nincbound  = NINT ( nzincbound * zdt / dtbak )
    aks2       = -zaks2
    aks4       = -zaks4
    lphys      = .FALSE.
    lcond      = .TRUE.
    itype_gscp = 2
  
    ! Calculate the filter weights
    nhalf   = nbak / 2
    CALL dolph (dtbak, taus, nhalf, hdfi, cheby, time, time2, weight, weight2)
  
    ! Copy the weighted initial fields into the DFI work arrays in COMDFI
    hfu  (:,:,:) = hdfi(0) * u (:,:,:,nnow)
    hfv  (:,:,:) = hdfi(0) * v (:,:,:,nnow)
    hfw  (:,:,:) = hdfi(0) * w (:,:,:,nnow)
    hft  (:,:,:) = hdfi(0) * t (:,:,:,nnow)
    hfqc (:,:,:) = hdfi(0) * qc_now(:,:,:)
    hfqv (:,:,:) = hdfi(0) * qv_now(:,:,:)
    hfpp (:,:,:) = hdfi(0) * pp(:,:,:,nnow)
  
    ! In the backward integration the boundary relaxation is done between
    ! 0 and -nincbound. The boundary values for variable psi and for
    ! -nincbound are calculated by:
    !   psi (nbd2) = psi(nbd1) - (psi(nbd2)-psi(nbd1)) = 2*psi(nbd1) - psi(nbd2)
    ! For the humidity variables the MAX of this value and 0.0 is taken.
    ! At the end of the forward integration the original values are calculated
    ! by the inverse:
    !  psi (nbd2) = 0.5 * (psi(nbd1) + psi(nbd2))
    ! For the humidity variables the original values have to be saved, because
    ! with the MAX-function this computation is not reversible
    WHERE (t_bd (:,:,:,nbd2) /= undef)
      qvb  (:,:,:)      = qv_bd(:,:,:,nbd2)
      qcb  (:,:,:)      = qc_bd(:,:,:,nbd2)
      u_bd (:,:,:,nbd2) = 2 * u_bd (:,:,:,nbd1) - u_bd (:,:,:,nbd2)
      v_bd (:,:,:,nbd2) = 2 * v_bd (:,:,:,nbd1) - v_bd (:,:,:,nbd2)
      t_bd (:,:,:,nbd2) = 2 * t_bd (:,:,:,nbd1) - t_bd (:,:,:,nbd2)
      pp_bd(:,:,:,nbd2) = 2 * pp_bd(:,:,:,nbd1) - pp_bd(:,:,:,nbd2)
      qv_bd(:,:,:,nbd2) = MAX (2*qv_bd(:,:,:,nbd1)-qv_bd(:,:,:,nbd2),0.0_wp)
      qc_bd(:,:,:,nbd2) = MAX (2*qc_bd(:,:,:,nbd1)-qc_bd(:,:,:,nbd2),0.0_wp)
    ENDWHERE

    ! Time loop for adiabatic backward integration
    DO ntstep = 0, nbak-1

      ! Set nexch_tag dependend on the time step
      nexch_tag = MOD (ntstep, 24*3600/INT(dt))

      ! variables concerned with the time step
      dt2    = 2.0_wp * dt
      ed2dt  = 1.0_wp / dt2
      dtdeh  = dt / 3600.0_wp
      nehddt = NINT ( 3600.0_wp / dt )

      ! preset values for time level nnew
      ! factors for linear time interpolation
      z2 = REAL (ntstep+1-nlastbound, wp) / REAL (nincbound, wp)
      z1 = 1.0_wp - z2

      ! fields of the atmosphere
      IF (lbd_frame) THEN
        WHERE (t_bd (:,:,:,nbd2) /= undef)
          u (:,:,:,nnew) = z1 * u_bd (:,:,:,nbd1) + z2 * u_bd (:,:,:,nbd2)
          v (:,:,:,nnew) = z1 * v_bd (:,:,:,nbd1) + z2 * v_bd (:,:,:,nbd2)
          t (:,:,:,nnew) = z1 * t_bd (:,:,:,nbd1) + z2 * t_bd (:,:,:,nbd2)
          pp(:,:,:,nnew) = z1 * pp_bd(:,:,:,nbd1) + z2 * pp_bd(:,:,:,nbd2)
          qv_new(:,:,:)  = z1 * qv_bd(:,:,:,nbd1) + z2 * qv_bd(:,:,:,nbd2)
          qc_new(:,:,:)  = z1 * qc_bd(:,:,:,nbd1) + z2 * qc_bd(:,:,:,nbd2)
        ELSEWHERE
          u (:,:,:,nnew) = u (:,:,:,nnow)
          v (:,:,:,nnew) = v (:,:,:,nnow)
          t (:,:,:,nnew) = t (:,:,:,nnow)
          pp(:,:,:,nnew) = pp(:,:,:,nnow)
          qv_new(:,:,:)  = qv_now(:,:,:)
          qc_new(:,:,:)  = qc_now(:,:,:)
        END WHERE
      ELSE
        u (:,:,:,nnew) = z1 * u_bd (:,:,:,nbd1) + z2 * u_bd (:,:,:,nbd2)
        v (:,:,:,nnew) = z1 * v_bd (:,:,:,nbd1) + z2 * v_bd (:,:,:,nbd2)
        t (:,:,:,nnew) = z1 * t_bd (:,:,:,nbd1) + z2 * t_bd (:,:,:,nbd2)
        pp(:,:,:,nnew) = z1 * pp_bd(:,:,:,nbd1) + z2 * pp_bd(:,:,:,nbd2)
        qv_new(:,:,:)  = z1 * qv_bd(:,:,:,nbd1) + z2 * qv_bd(:,:,:,nbd2)
        qc_new(:,:,:)  = z1 * qc_bd(:,:,:,nbd1) + z2 * qc_bd(:,:,:,nbd2)
      ENDIF


      ! initialize tendency fields with 0.0
      utens   (:,:,:) = 0.0_wp
      vtens   (:,:,:) = 0.0_wp
      wtens   (:,:,:) = 0.0_wp
      ttens   (:,:,:) = 0.0_wp
      qv_tens (:,:,:) = 0.0_wp
      qc_tens (:,:,:) = 0.0_wp
      pptens  (:,:,:) = 0.0_wp
      IF ( lphys .AND. ltur .AND. itype_turb==3 ) THEN
         tketens (:,:,:) = 0.0_wp
      END IF

      ! Calculate surface pressure for new time level and density of moist air
      CALL calps ( ps(:,:   ,nnew), pp(:,:,ke,nnew), t(:,:,ke,nnew),     &
                   qv_new(:,:,ke), qc_new(:,:,ke), qrs(:,:,ke)   ,       &
                   rho0(:,:,ke), p0(:,:,ke), dp0(:,:,ke),                &
                   ie, je, rvd_m_o, r_d,                                 &
                   istartpar, iendpar, jstartpar, jendpar )
      CALL calrho ( t(:,:,:,nnow), pp(:,:,:,nnow), qv_now(:,:,:),        &
                   qc_now(:,:,:), qrs, p0, rho, ie, je, ke, r_d, rvd_m_o)

      IF (ltime) CALL get_timings (i_add_computations, 0, dt, izerror)

      ! perform the dynamics and the relaxation
      CALL organize_dynamics ('compute', izerror, yzerrmsg, -zdt, .TRUE.)
      IF (ltime) CALL get_timings (i_dyn_computations, 0, dt, izerror)

      CALL organize_dynamics ('relaxation', izerror, yzerrmsg, -zdt, .TRUE.)
      IF (ltime) CALL get_timings (i_relaxation, 0, dt, izerror)

      WHERE (qc_new (:,:,:) < 1.0E-15_wp)
        qc_new (:,:,:) = 0.0_wp
      ENDWHERE
      WHERE (qv_new (:,:,:) < 1.0E-15_wp)
        qv_new (:,:,:) = 0.0_wp
      ENDWHERE
      WHERE (qrs(:,:,:) < 1.0E-15_wp)
        qrs(:,:,:) = 0.0_wp
      ENDWHERE

      IF (ltime_barrier) THEN
        CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
        IF (ltime) CALL get_timings (i_barrier_waiting_dyn, 0, dt, izerror)
      ENDIF

      kzdims(1:24) =                                                       &
         (/ke,ke,ke,ke,ke1,ke1,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                &
         (nnew+30, sendbuf, isendbuflen, imp_reals, icomm_cart,            &
          num_compute, ie, je, kzdims, jstartpar, jendpar, 2, nboundlines, &
          my_cart_neigh, lperi_x, lperi_y, l2dim,                          &
          100+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,        &
          u (:,:,:,nnow), u (:,:,:,nnew), v (:,:,:,nnow), v (:,:,:,nnew),  &
          w (:,:,:,nnow), w (:,:,:,nnew), t (:,:,:,nnow), t (:,:,:,nnew),  &
          pp(:,:,:,nnow), pp(:,:,:,nnew))

      kzdims(1:24) =                                                       &
         (/ke,ke,ke,ke,ke1,ke1,ke,ke,ke,ke,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                &
         ( 0     , sendbuf, isendbuflen, imp_reals, icomm_cart,            &
          num_compute, ie, je, kzdims, jstartpar, jendpar, 2, nboundlines, &
          my_cart_neigh, lperi_x, lperi_y, l2dim,                          &
          100+nexch_tag, .FALSE.   , ncomm_type, izerror, yzerrmsg,        &
          qv_now(:,:,:) , qv_new(:,:,:) , qc_now(:,:,:) , qc_new(:,:,:))

      IF (ltime) CALL get_timings (i_communications_dyn, 0, dt, izerror)

      ! calculate surface pressure
      CALL calps ( ps(:,:,nnow),    pp(:,:,ke,nnow), t  (:,:,ke,nnow),   &
                   qv_now(:,:,ke) , qc_now(:,:,ke) , qrs(:,:,ke),        &
                   rho0(:,:,ke), p0(:,:,ke), dp0(:,:,ke),                &
                   ie, je, rvd_m_o, r_d, 1, ie, 1, je)
      CALL calps ( ps(:,:,nnew),    pp(:,:,ke,nnew), t  (:,:,ke,nnew),   &
                   qv_new(:,:,ke) , qc_new(:,:,ke) , qrs(:,:,ke),        &
                   rho0(:,:,ke), p0(:,:,ke), dp0(:,:,ke),                &
                   ie, je, rvd_m_o, r_d, 1, ie, 1, je)

      ! double dt after first step
      IF (ntstep == 0 .AND. .NOT.l2tls) THEN
        dt = 2.0_wp * dt
      ENDIF
  
      ! accumulate weighted fields in work arrays
      hfu  (:,:,:) = hfu (:,:,:) + hdfi(ntstep+1) * u (:,:,:,nnew)
      hfv  (:,:,:) = hfv (:,:,:) + hdfi(ntstep+1) * v (:,:,:,nnew)
      hfw  (:,:,:) = hfw (:,:,:) + hdfi(ntstep+1) * w (:,:,:,nnew)
      hft  (:,:,:) = hft (:,:,:) + hdfi(ntstep+1) * t (:,:,:,nnew)
      hfqc (:,:,:) = hfqc(:,:,:) + hdfi(ntstep+1) * qc_new(:,:,:)
      hfqv (:,:,:) = hfqv(:,:,:) + hdfi(ntstep+1) * qv_new(:,:,:)
      hfpp (:,:,:) = hfpp(:,:,:) + hdfi(ntstep+1) * pp(:,:,:,nnew)
  
      ! cyclic change of time levels
      IF ( .NOT.l2tls ) THEN
        nsp  = nold
        nold = nnow
        nnow = nnew
        nnew = nsp
!     not possible for ndfi=2
!     ELSE
!       nnow = 3-nnow
!       nnew = 3-nnew
      END IF

      CALL update_trcr_pointers()

    ENDDO  ! adiabatic backward integration
  
    ! copy filtered fields to prognostic fields as initial values 
    ! (time index 1 and 2)
    DO nx = 1, 2

      CALL trcr_get(izerror, idt_qv, ptr_tlev = nx, ptr = qv_nx)
      IF (izerror /= 0_iintegers) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF
      CALL trcr_get(izerror, idt_qc, ptr_tlev = nx, ptr = qc_nx)
      IF (izerror /= 0_iintegers) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF

      u  (:,:,:,nx) = hfu (:,:,:)
      v  (:,:,:,nx) = hfv (:,:,:)
      w  (:,:,:,nx) = hfw (:,:,:)
      t  (:,:,:,nx) = hft (:,:,:)
      pp (:,:,:,nx) = hfpp(:,:,:)
      qc_nx(:,:,:)  = hfqc(:,:,:)
      qv_nx(:,:,:)  = hfqv(:,:,:)

      ! calculate surface pressure
      CALL calps ( ps(:,:,nx), pp(:,:,ke,nx), t(:,:,ke,nx), qv_nx(:,:,ke),    &
                   qc_nx(:,:,ke), qrs(:,:,ke), rho0(:,:,ke), p0(:,:,ke),      &
                   dp0(:,:,ke), ie, je, rvd_m_o, r_d, 1, ie, 1, je)
    ENDDO
  
    ! reset settings
    aks2 = zaks2
    aks4 = zaks4

    ! Reset boundary values for nbd2 by the inverse:
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Bugfix in 3.23: this was wrong before
    !! WRONG   psi (nbd2) = 0.5 * (psi(nbd1) + psi(nbd2)) !! WRONG
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  psi (nbd2) (+1h) = 2. * psi(nbd1) (0h)  - psi(nbd2) (-1h)
    WHERE (t_bd (:,:,:,nbd2) /= undef)
      u_bd (:,:,:,nbd2) = 2.0_wp * u_bd (:,:,:,nbd1) - u_bd (:,:,:,nbd2)
      v_bd (:,:,:,nbd2) = 2.0_wp * v_bd (:,:,:,nbd1) - v_bd (:,:,:,nbd2)
      t_bd (:,:,:,nbd2) = 2.0_wp * t_bd (:,:,:,nbd1) - t_bd (:,:,:,nbd2)
      pp_bd(:,:,:,nbd2) = 2.0_wp * pp_bd(:,:,:,nbd1) - pp_bd(:,:,:,nbd2)
      qv_bd(:,:,:,nbd2) = qvb(:,:,:)
      qc_bd(:,:,:,nbd2) = qcb(:,:,:)
    ENDWHERE

  ENDIF ! ndfi == 2

!------------------------------------------------------------------------------
! Section 3: Diabatic forward integration
!------------------------------------------------------------------------------

  ! settings for forward integration
  IF ( .NOT.l2tls ) THEN
    nold = 1
    nnow = 2
    nnew = 3
    dt   = dtfwd * 0.5_wp   ! only in the first step
    nx0  = nold
  ELSE
    nnow = 2
    nnew = 1
    dt   = dtfwd
    nx0  = nnow
    nstop = nfwd - 1
  END IF

  CALL update_trcr_pointers()

  ! Special timelevel nx0

  CALL trcr_get(izerror, idt_qv, ptr_tlev = nx0, ptr = qv_nx)
  IF (izerror /= 0_iintegers) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nx0, ptr = qc_nx)
  IF (izerror /= 0_iintegers) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qi, ptr_tlev = nx0, ptr = qi_nx)
  IF (izerror /= 0_iintegers .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  IF (itype_dfi_qx > 0) THEN
    IF (itype_dfi_qx == 2) ALLOCATE(zrh(ie,je,ke),STAT=istat)
    ALLOCATE(zqc(ie,je,ke),STAT=istat)
    ALLOCATE(zqi(ie,je,ke),STAT=istat)
    ALLOCATE(zqr(ie,je,ke),STAT=istat)
    ALLOCATE(zqs(ie,je,ke),STAT=istat)
    ALLOCATE(zqg(ie,je,ke),STAT=istat)
  ENDIF

  ! Resetting of control switches
  ! New in V4_18: Soil-Model, FLake and Seaice scheme are now used for the
  !               forward integration, if ndfi=1 is chosen
  lphys      = lzphys
  lcond      = lzcond
  lsoil      = lzsoil
  llake      = lzlake
  lseaice    = lzseaice
  itype_gscp = iztype_gscp

  IF (ndfi == 2) THEN

    ! When using backward and forward integration, no soil processes are
    ! taken into account
    lsoil      = .FALSE.
    llake      = .FALSE.
    lseaice    = .FALSE.

    ! In the forward integration the boundary relaxation is done between
    ! -nincound and +nincbound. The boundary values for variable psi and for
    ! -nincbound are calculated by:
    !  psi (nbd1) = psi(nbd1) - (psi(nbd2)-psi(nbd1)) = 2*psi(nbd1) - psi(nbd2)
    ! For the humidity variables the MAX from these values and 0.0 are taken.
    ! At the end of the forward integration the original values are calculated
    ! by the inverse:
    !  psi (nbd1) = 0.5 * (psi(nbd1) + psi(nbd2))
    ! For the humidity variables the original values have to be saved, because
    ! with the MAX-function this computation is not reversible
    nincbound  =   NINT ( 2.0_wp * nzincbound * zdt / dtfwd)
    nlastbound = - NINT ( (nzincbound*zdt - nbak*dtbak*0.5_wp) / dtfwd )
    WHERE (t_bd (:,:,:,nbd2) /= undef)
      qvb  (:,:,:)      = qv_bd(:,:,:,nbd1)
      qcb  (:,:,:)      = qc_bd(:,:,:,nbd1)

      u_bd (:,:,:,nbd1) = 2.0_wp * u_bd (:,:,:,nbd1) - u_bd (:,:,:,nbd2)
      v_bd (:,:,:,nbd1) = 2.0_wp * v_bd (:,:,:,nbd1) - v_bd (:,:,:,nbd2)
      t_bd (:,:,:,nbd1) = 2.0_wp * t_bd (:,:,:,nbd1) - t_bd (:,:,:,nbd2)
      pp_bd(:,:,:,nbd1) = 2.0_wp * pp_bd(:,:,:,nbd1) - pp_bd(:,:,:,nbd2)
      qv_bd(:,:,:,nbd1) = MAX (2.0_wp * qv_bd (:,:,:,nbd1) - qv_bd (:,:,:,nbd2),   &
                               0.0_wp)
      qc_bd(:,:,:,nbd1) = MAX (2.0_wp * qc_bd (:,:,:,nbd1) - qc_bd (:,:,:,nbd2),   &
                               0.0_wp)
    ENDWHERE

    IF (ASSOCIATED(qi_now)) THEN
      IF (izlbc_qi == T_LBC_FILE) THEN
        WHERE (qi_bd (:,:,:,nbd2) /= undef)
          qib  (:,:,:)      = qi_bd(:,:,:,nbd1)
          qi_bd(:,:,:,nbd1) = MAX (2.0_wp * qi_bd(:,:,:,nbd1) - qi_bd(:,:,:,nbd2), &
                                   0.0_wp)
        ENDWHERE
      ENDIF
    ENDIF
    IF (ASSOCIATED(qr_now)) THEN
      IF (izlbc_qr == T_LBC_FILE) THEN
        WHERE (qr_bd (:,:,:,nbd2) /= undef)
          qrb  (:,:,:)      = qr_bd(:,:,:,nbd1)
          qr_bd(:,:,:,nbd1) = MAX (2.0_wp * qr_bd(:,:,:,nbd1) - qr_bd(:,:,:,nbd2), &
                                   0.0_wp)
        END WHERE
      ENDIF
    ENDIF
    IF (ASSOCIATED(qs_now)) THEN
      IF (izlbc_qs == T_LBC_FILE) THEN
        WHERE (qs_bd (:,:,:,nbd2) /= undef)
          qsb  (:,:,:)      = qs_bd(:,:,:,nbd1)
          qs_bd(:,:,:,nbd1) = MAX (2.0_wp * qs_bd(:,:,:,nbd1) - qs_bd(:,:,:,nbd2), &
                                   0.0_wp)
        END WHERE
      ENDIF
    ENDIF
    IF (ASSOCIATED(qg_now)) THEN
      IF (izlbc_qg == T_LBC_FILE) THEN
        WHERE (qg_bd (:,:,:,nbd2) /= undef)
          qgb  (:,:,:)      = qg_bd(:,:,:,nbd1)
          qg_bd(:,:,:,nbd1) = MAX (2.0_wp * qg_bd(:,:,:,nbd1) - qg_bd(:,:,:,nbd2), &
                                   0.0_wp)
        END WHERE
      ENDIF
    ENDIF

  ELSEIF (ndfi == 1) THEN

    ! only a forward integration is performed between 0 and nincbound. 
    ! The boundary values need not be changed for that
    nincbound  = NINT ( nzincbound * zdt / dtfwd )
    nlastbound = 0_iintegers

    ! Allocate space and save initial variables for soil and surface variables
    IF (itype_dfi_soil == 0) THEN
      ALLOCATE (zts          (ie,je), STAT=izerror)
      IF (izerror == 0) THEN
        zts(:,:)  = t_s(:,:,nx0)
      ENDIF

      ALLOCATE (zqv_s        (ie,je), STAT=izerror)
      IF (izerror == 0) THEN
        zqv_s(:,:) = qv_s(:,:,nx0)
      ENDIF

      ALLOCATE (zw_snow      (ie,je), STAT=izerror)
      IF (izerror == 0) THEN
        zw_snow(:,:) = w_snow(:,:,nx0)
      ENDIF

      ALLOCATE (zt_snow      (ie,je), STAT=izerror)
      IF (izerror == 0) THEN
        zt_snow(:,:) = t_snow(:,:,nx0)
      ENDIF

      IF (lsoil) THEN
        IF (lmulti_layer) THEN
          ALLOCATE (zt_so    (ie,je,0:ke_soil+1), STAT=izerror)
          IF (izerror == 0) THEN
            zt_so(:,:,:)  = t_so(:,:,:,nx0)
          ENDIF

          ALLOCATE (zw_so    (ie,je,ke_soil+1), STAT=izerror)
          IF (izerror == 0) THEN
            zw_so(:,:,:)  = w_so(:,:,:,nx0)
          ENDIF

          ALLOCATE (zfreshsnw(ie,je), STAT=izerror)
          IF (izerror == 0) THEN
            zfreshsnw(:,:)  = freshsnow(:,:)
          ENDIF

          IF (lana_rho_snow) THEN
            ALLOCATE (zrhosnw(ie,je), STAT=izerror)
            IF (izerror == 0) THEN
              zrhosnw(:,:)  = rho_snow(:,:,nx0)
            ENDIF
          ENDIF
        ELSE
          ALLOCATE (zt_m     (ie,je), STAT=izerror)
          IF (izerror == 0) THEN
            zt_m(:,:)  = t_m(:,:,nx0)
          ENDIF

          ALLOCATE (zw_g1    (ie,je), STAT=izerror)
          IF (izerror == 0) THEN
            zw_g1(:,:)  = w_g1(:,:,nx0)
          ENDIF

          ALLOCATE (zw_g2    (ie,je), STAT=izerror)
          IF (izerror == 0) THEN
            zw_g2(:,:)  = w_g2(:,:,nx0)
          ENDIF

          IF ( nlgw == 3 ) THEN
            ALLOCATE (zw_g3  (ie,je), STAT=izerror)
            IF (izerror == 0) THEN
              zw_g3(:,:)  = w_g3(:,:,nx0)
            ENDIF
          ENDIF
        ENDIF
      ENDIF

      IF (llake) THEN
        ALLOCATE (ztmnw_lk   (ie,je), STAT=izerror)
        IF (izerror == 0) THEN
          ztmnw_lk(:,:)  = t_mnw_lk(:,:,nx0)
        ENDIF

        ALLOCATE (ztwml_lk   (ie,je), STAT=izerror)
        IF (izerror == 0) THEN
          ztwml_lk(:,:)  = t_wml_lk(:,:,nx0)
        ENDIF

        ALLOCATE (ztbot_lk   (ie,je), STAT=izerror)
        IF (izerror == 0) THEN
          ztbot_lk(:,:)  = t_bot_lk(:,:,nx0)
        ENDIF

        ALLOCATE (zct_lk     (ie,je), STAT=izerror)
        IF (izerror == 0) THEN
          zct_lk(:,:)  = c_t_lk(:,:,nx0)
        ENDIF

        ALLOCATE (zhml_lk    (ie,je), STAT=izerror)
        IF (izerror == 0) THEN
          zhml_lk(:,:)  = h_ml_lk(:,:,nx0)
        ENDIF
      ENDIF

      IF (lseaice .OR. llake) THEN
        ALLOCATE (ztice      (ie,je), STAT=izerror)
        IF (izerror == 0) THEN
          ztice(:,:)  = t_ice(:,:,nx0)
        ENDIF

        ALLOCATE (zhice      (ie,je), STAT=izerror)
        IF (izerror == 0) THEN
          zhice(:,:)  = h_ice(:,:,nx0)
        ENDIF
      ENDIF

      IF (izerror /= 0) THEN
        ierror = 10
        yerrmsg = 'Problem allocating space for ndfi=1'
        RETURN
      ENDIF
    ENDIF
  ENDIF


  ! Calculate the filter weights
  nhalf   = nfwd / 2
  CALL dolph (dtfwd, taus, nhalf, hdfi, cheby, time, time2, weight, weight2)

  ! If the ice-scheme is turned on, but no initial values are given for qi,
  ! initial values have to be set according to lmorg:
  IF (ASSOCIATED(qi_now)) THEN
    IF (izic_qi /= T_INI_FILE) THEN
      DO k = 1, ke
        qi_now(:,:,k) = 0.0_wp
        IF (.NOT. l2tls) THEN
          qi_old(:,:,k) = 0.0_wp
        ENDIF

        IF ( MINVAL(t(:,:,k,nnow)) < 248.15_wp ) THEN
          DO j = 1, je
            DO i = 1, ie
              IF ( t(i,j,k,nnow) < 248.15_wp ) THEN
                qi_now(i,j,k) = qc_now(i,j,k)
                qc_now(i,j,k) = 0.0_wp
! this reduction is only useful for GME data without cloud ice
! put all these computations into the INT2LM later
!               qv(i,j,k,nnow) = qv(i,j,k,nnow) &
!                  *fqvs( fpvsi(t(i,j,k,nnow)), p0(i,j,k)+pp(i,j,k,nnow) ) &
!                  /fqvs( fpvsw(t(i,j,k,nnow)), p0(i,j,k)+pp(i,j,k,nnow) )

                IF (.NOT. l2tls) THEN
                  qi_old(i,j,k) = qi_now(i,j,k)
                  qc_old(i,j,k) = 0.0_wp
                  qv_old(i,j,k) = qv_now(i,j,k)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDIF

  ENDIF

  ! Copy the weighted initial fields into the DFI work arrays in COMDFI
  hfu  (:,:,:) = hdfi(0) * u (:,:,:,nnow)
  hfv  (:,:,:) = hdfi(0) * v (:,:,:,nnow)
  hfw  (:,:,:) = hdfi(0) * w (:,:,:,nnow)
  hft  (:,:,:) = hdfi(0) * t (:,:,:,nnow)
  hfpp (:,:,:) = hdfi(0) * pp(:,:,:,nnow)

  IF (itype_dfi_qx == 0) THEN
    ! filter all variables 
    hfqv (:,:,:) = hdfi(0) * qv_now(:,:,:)
    hfqc (:,:,:) = hdfi(0) * qc_now(:,:,:)
    IF (ASSOCIATED(qi_now)) THEN
      hfqi (:,:,:) = hdfi(0) * qi_now(:,:,:)
    ENDIF
    IF (ASSOCIATED(qr_now)) THEN
      hfqr (:,:,:) =  hdfi(0) * qr_now(:,:,:)
    ENDIF
    IF (ASSOCIATED(qs_now)) THEN
      hfqs (:,:,:) =  hdfi(0) * qs_now(:,:,:)
    ENDIF
    IF (ASSOCIATED(qg_now)) THEN
      hfqg (:,:,:) = hdfi(0) * qg_now(:,:,:)
    ENDIF
  ELSE
    ! filter only qv
    hfqv (:,:,:) = hdfi(0) * qv_now(:,:,:)
  ENDIF

  ! Time loop for diabatic forward integration
  DO ntstep = 0, nfwd-1

    IF ( izdebug > 1 ) THEN
      PRINT *, 'DFI STEP ',ntstep
    ENDIF

    ! Set nexch_tag dependend on the time step
    nexch_tag = MOD (ntstep, 24*3600/INT(dt))

    ! variables concerned with the time step
    dt2    = 2.0_wp * dt
    ed2dt  = 1.0_wp / dt2
    dtdeh  = dt / 3600.0_wp
    nehddt = NINT ( 3600.0_wp / dt )

    ! preset values for time level nnew
    ! factors for linear time interpolation
    z2 = REAL (ntstep+1-nlastbound, wp) / REAL (nincbound, wp)
    z2 = MIN ( 1.0_wp , z2 )
    z1 = 1.0_wp - z2

    ! fields of the atmosphere
    IF (lbd_frame) THEN
      WHERE (t_bd (:,:,:,nbd2) /= undef)
        u (:,:,:,nnew) = z1 * u_bd (:,:,:,nbd1) + z2 * u_bd (:,:,:,nbd2)
        v (:,:,:,nnew) = z1 * v_bd (:,:,:,nbd1) + z2 * v_bd (:,:,:,nbd2)
        t (:,:,:,nnew) = z1 * t_bd (:,:,:,nbd1) + z2 * t_bd (:,:,:,nbd2)
        pp(:,:,:,nnew) = z1 * pp_bd(:,:,:,nbd1) + z2 * pp_bd(:,:,:,nbd2)
        qv_new(:,:,:)  = z1 * qv_bd(:,:,:,nbd1) + z2 * qv_bd(:,:,:,nbd2)
        qc_new(:,:,:)  = z1 * qc_bd(:,:,:,nbd1) + z2 * qc_bd(:,:,:,nbd2)
      ELSEWHERE
        u (:,:,:,nnew) = u (:,:,:,nnow)
        v (:,:,:,nnew) = v (:,:,:,nnow)
        t (:,:,:,nnew) = t (:,:,:,nnow)
        pp(:,:,:,nnew) = pp(:,:,:,nnow)
        qv_new(:,:,:)  = qv_now(:,:,:)
        qc_new(:,:,:)  = qc_now(:,:,:)
      END WHERE
      IF (.NOT. lw_freeslip) THEN
        WHERE (t_bd (:,:,:,nbd2) /= undef)
          w(:,:,:,nnew) = z1 *  w_bd(:,:,:,nbd1) + z2 *  w_bd(:,:,:,nbd2)
        ELSEWHERE
          w(:,:,:,nnew) =  w(:,:,:,nnow)
        ENDWHERE
      ENDIF
    ELSE
      u (:,:,:,nnew) = z1 * u_bd (:,:,:,nbd1) + z2 * u_bd (:,:,:,nbd2)
      v (:,:,:,nnew) = z1 * v_bd (:,:,:,nbd1) + z2 * v_bd (:,:,:,nbd2)
      t (:,:,:,nnew) = z1 * t_bd (:,:,:,nbd1) + z2 * t_bd (:,:,:,nbd2)
      pp(:,:,:,nnew) = z1 * pp_bd(:,:,:,nbd1) + z2 * pp_bd(:,:,:,nbd2)
      qv_new(:,:,:)  = z1 * qv_bd(:,:,:,nbd1) + z2 * qv_bd(:,:,:,nbd2)
      qc_new(:,:,:)  = z1 * qc_bd(:,:,:,nbd1) + z2 * qc_bd(:,:,:,nbd2)
      IF (.NOT. lw_freeslip) THEN
        w (:,:,:,nnew) = z1 * w_bd (:,:,:,nbd1) + z2 * w_bd (:,:,:,nbd2)
      ENDIF
    ENDIF

    IF (ASSOCIATED(qi_now)) THEN
      IF (izlbc_qi == T_INI_FILE) THEN
        IF (lbd_frame) THEN
          WHERE (t_bd (:,:,:,nbd2) /= undef)
            qi_new(:,:,:) = z1 * qi_bd(:,:,:,nbd1) + z2 * qi_bd(:,:,:,nbd2)
          ELSEWHERE
            qi_new(:,:,:) = qi_now(:,:,:)
          ENDWHERE
        ELSE
          qi_new(:,:,:) = z1 * qi_bd(:,:,:,nbd1) + z2 * qi_bd(:,:,:,nbd2)
        ENDIF
      ELSE
        ! treatment of qi analogous to initialize_loop in lmorg
        DO k = 1, ke
          qi_new(:,:,k) = 0.0_wp
          IF ( MINVAL(t(:,:,k,nnew)) < 248.15_wp ) THEN
            DO j = 1, je
              DO i = 1, ie
                IF ( t(i,j,k,nnew) < 248.15_wp ) THEN
                  qi_new(i,j,k) = qc_new(i,j,k)
                  qc_new(i,j,k) = 0.0_wp
! this reduction is only useful for GME data without cloud ice
! put all these computations into the INT2LM later
!                 qv_new(i,j,k) = qv_new(i,j,k) &
!                        *fqvs( fpvsi(t(i,j,k,nnew)), p0(i,j,k)+pp(i,j,k,nnew) ) &
!                        /fqvs( fpvsw(t(i,j,k,nnew)), p0(i,j,k)+pp(i,j,k,nnew) )
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF
    ENDIF

  ! the initial values are reformulated similarily
    IF (ASSOCIATED(qi_now)) THEN
      IF ( (ntstep == 0) .AND. (izic_qi /= T_INI_FILE) ) THEN
      DO k = 1, ke
        IF ( MINVAL(t(:,:,k,nnow)) < 248.15_wp ) THEN
          DO j = 1, je
            DO i = 1, ie
              IF ( t(i,j,k,nnow) < 248.15_wp ) THEN
                qi_now(i,j,k) = qc_now(i,j,k)
                qi_nx (i,j,k) = qi_now(i,j,k)
                qc_now(i,j,k) = 0.0_wp
                qc_nx (i,j,k) = 0.0_wp
! this reduction is only useful for GME data without cloud ice
! put all these computations into the INT2LM later
!               qv_now(i,j,k) = qv_now(i,j,k) &
!                      *fqvs( fpvsi(t(i,j,k,nnow)), p0(i,j,k)+pp(i,j,k,nnow) ) &
!                      /fqvs( fpvsw(t(i,j,k,nnow)), p0(i,j,k)+pp(i,j,k,nnow) )
!               qv_nx(i,j,k)  = qv_now(i,j,k)
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      ENDIF
    ENDIF

  ! prognostic precipitation fields
    IF (ASSOCIATED(qr_new) .AND. ASSOCIATED(qs_new)) THEN
      IF (izlbc_qr == T_LBC_FILE .AND. izlbc_qs == T_LBC_FILE) THEN
        IF (lbd_frame) THEN
          WHERE (qr_bd (:,:,:,nbd2) /= undef)
            qr_new(:,:,:) = z1 * qr_bd(:,:,:,nbd1) + z2 * qr_bd(:,:,:,nbd2)
            qs_new(:,:,:) = z1 * qs_bd(:,:,:,nbd1) + z2 * qs_bd(:,:,:,nbd2)
          ELSEWHERE
            qr_new(:,:,:) = qr_now(:,:,:)
            qs_new(:,:,:) = qs_now(:,:,:)
          END WHERE
        ELSE
          qr_new(:,:,:) = z1 * qr_bd(:,:,:,nbd1) + z2 * qr_bd(:,:,:,nbd2)
          qs_new(:,:,:) = z1 * qs_bd(:,:,:,nbd1) + z2 * qs_bd(:,:,:,nbd2)
        ENDIF
      ENDIF
    ENDIF
    IF (ASSOCIATED(qg_new)) THEN
      IF (izlbc_qg == T_LBC_FILE) THEN
        IF (lbd_frame) THEN
          WHERE (qg_bd (:,:,:,nbd2) /= undef)
            qg_new(:,:,:) = z1 * qg_bd(:,:,:,nbd1) + z2 * qg_bd(:,:,:,nbd2)
          ELSEWHERE
            qg_new(:,:,:) = qg_now(:,:,:)
          END WHERE
        ELSE
          qg_new(:,:,:) = z1 * qg_bd(:,:,:,nbd1) + z2 * qg_bd(:,:,:,nbd2)
        ENDIF
      ENDIF
    ENDIF

    ! initialize tendency fields with 0.0
    utens  (:,:,:) = 0.0_wp
    vtens  (:,:,:) = 0.0_wp
    wtens  (:,:,:) = 0.0_wp
    ttens  (:,:,:) = 0.0_wp
    qv_tens (:,:,:)= 0.0_wp
    qc_tens (:,:,:)= 0.0_wp
    IF(ASSOCIATED(qi_tens)) qi_tens (:,:,:)= 0.0_wp
    pptens (:,:,:) = 0.0_wp
    IF ( lphys .AND. ltur .AND. itype_turb==3 ) THEN
       tketens (:,:,:) = 0.0_wp
    END IF

    ! Calculate surface pressure for new time level and density of moist air
    CALL calps ( ps(:,:   ,nnew), pp(:,:,ke,nnew), t(:,:,ke,nnew),     &
                 qv_new(:,:,ke ), qc_new(:,:,ke) , qrs(:,:,ke)   ,     &
                 rho0(:,:,ke), p0(:,:,ke), dp0(:,:,ke),                &
                 ie, je, rvd_m_o, r_d,                                  &
                 istartpar, iendpar, jstartpar, jendpar )
    CALL calrho ( t(:,:,:,nnow), pp(:,:,:,nnow), qv_now(:,:,:),        &
                 qc_now(:,:,:), qrs, p0, rho, ie, je, ke, r_d, rvd_m_o)

    IF (ltime) CALL get_timings (i_add_computations, 0, dt, izerror)

    ! perform the physics, dynamics, relaxation and boundary exchange
    IF (lphys) CALL organize_physics ('compute', izerror, yzerrmsg)

    IF (ltime) CALL get_timings (i_phy_computations, 0, dt, izerror)

    CALL organize_dynamics ('compute', izerror, yzerrmsg, zdt, .TRUE.)
    IF (ltime) CALL get_timings (i_dyn_computations, 0, dt, izerror)

    IF (ASSOCIATED(qr_new)) THEN
        DO  k = 1, ke
          ! western boundary
          IF (my_cart_neigh(1) == -1) THEN
            DO  j = jstart, jend
              DO i = 1, nboundlines
                qr_new(i,j,k) = qr_new(istart,j,k)
              ENDDO
            ENDDO
          ENDIF
          ! eastern boundary
          IF (my_cart_neigh(3) == -1) THEN
            DO  j = jstart, jend
              DO i = ie-nboundlines+1, ie
                qr_new(i,j,k) = qr_new(iend  ,j,k)
              ENDDO
            ENDDO
          ENDIF
          ! southern boundary
          IF (my_cart_neigh(4) == -1) THEN
            DO  j = 1, nboundlines
              qr_new(:,j,k) = qr_new(:,jstart,k)
            ENDDO
          ENDIF
          ! northern boundary
          IF (my_cart_neigh(2) == -1) THEN
            DO  j = je-nboundlines+1, je
              qr_new(:,j,k) = qr_new(:,jend  ,k)
            ENDDO
          ENDIF
        ENDDO
    ENDIF
      ! qs
    IF (ASSOCIATED(qs_new)) THEN
        DO  k = 1, ke
          ! western boundary
          IF (my_cart_neigh(1) == -1) THEN
            DO  j = jstart, jend
              DO i = 1, nboundlines
                qs_new(i,j,k) = qs_new(istart,j,k)
              ENDDO
            ENDDO
          ENDIF
          ! eastern boundary
          IF (my_cart_neigh(3) == -1) THEN
            DO  j = jstart, jend
              DO i = ie-nboundlines+1, ie
                qs_new(i,j,k) = qs_new(iend  ,j,k)
              ENDDO
            ENDDO
          ENDIF
          ! southern boundary
          IF (my_cart_neigh(4) == -1) THEN
            DO  j = 1, nboundlines
              qs_new(:,j,k) = qs_new(:,jstart,k)
            ENDDO
          ENDIF
          ! northern boundary
          IF (my_cart_neigh(2) == -1) THEN
            DO  j = je-nboundlines+1, je
              qs_new(:,j,k) = qs_new(:,jend  ,k)
            ENDDO
          ENDIF
        ENDDO
    ENDIF

    IF (ASSOCIATED(qg_new)) THEN
        DO  k = 1, ke
          ! western boundary
          IF (my_cart_neigh(1) == -1) THEN
            DO  j = jstart, jend
              DO i = 1, nboundlines
                qg_new(i,j,k) = qg_new(istart,j,k)
              ENDDO
            ENDDO
          ENDIF
          ! eastern boundary
          IF (my_cart_neigh(3) == -1) THEN
            DO  j = jstart, jend
              DO i = ie-nboundlines+1, ie
                qg_new(i,j,k) = qg_new(iend  ,j,k)
              ENDDO
            ENDDO
          ENDIF
          ! southern boundary
          IF (my_cart_neigh(4) == -1) THEN
            DO  j = 1, nboundlines
              qg_new(:,j,k) = qg_new(:,jstart,k)
            ENDDO
          ENDIF
          ! northern boundary
          IF (my_cart_neigh(2) == -1) THEN
            DO  j = je-nboundlines+1, je
              qg_new(:,j,k) = qg_new(:,jend  ,k)
            ENDDO
          ENDIF
        ENDDO
    ENDIF

    CALL organize_dynamics ('relaxation', izerror, yzerrmsg, zdt, .TRUE.)
    IF (ltime) CALL get_timings (i_relaxation, 0, dt, izerror)

    ! Final update of temperature and humidity variables due to
    ! cloud microphysics in case of the cloud ice scheme
    IF (lphys) CALL organize_physics ('finish_compute', izerror, yzerrmsg)
    IF (ltime) CALL get_timings (i_phy_computations, 0, dt, izerror)

    IF (ASSOCIATED(qr_new)) THEN
      WHERE (qr_new (:,:,:) < 1.0E-15_wp)
        qr_new (:,:,:) = 0.0_wp
      ENDWHERE
    ENDIF
    IF (ASSOCIATED(qs_new)) THEN
      WHERE (qs_new (:,:,:) < 1.0E-15_wp)
        qs_new (:,:,:) = 0.0_wp
      ENDWHERE
    ENDIF
    IF (ASSOCIATED(qg_new)) THEN
      WHERE (qg_new (:,:,:) < 1.0E-15_wp)
        qg_new (:,:,:) = 0.0_wp
      ENDWHERE
    ENDIF
    IF (ASSOCIATED(qi_new)) THEN
      WHERE (qi_new (:,:,:) < 1.0E-15_wp)
        qi_new (:,:,:) = 0.0_wp
      ENDWHERE
    ENDIF
    WHERE (qc_new (:,:,:) < 1.0E-15_wp)
      qc_new (:,:,:) = 0.0_wp
    ENDWHERE
    WHERE (qv_new (:,:,:) < 1.0E-15_wp)
      qv_new (:,:,:) = 0.0_wp
    ENDWHERE
    WHERE (qrs(:,:,:) < 1.0E-15_wp)
      qrs(:,:,:) = 0.0_wp
    ENDWHERE

    ! Check, whether additional communication for the convection is
    ! necessary
    lzconv = lconv .AND. lconf_avg .AND.                         &
              ((ntstep < 1) .OR. (MOD(ntstep+2,nincconv)==0))

    IF (ltime_barrier) THEN
      CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
      IF (ltime) CALL get_timings (i_barrier_waiting_dyn, 0, dt, izerror)
    ENDIF

    IF     ( l2tls ) THEN
      CALL exchange_runge_kutta
    ELSE ! Leapfrog:
      CALL exchange_leapfrog
    ENDIF
    IF (ltime) CALL get_timings (i_communications_dyn, 0, dt, izerror)

    ! Using ndfi == 1 the lfff00000000/z/p are not produced, but for many
    ! operational applications these files are needed.

    IF (ntstep == 0 .AND. ndfi == 1) THEN

      CALL near_surface (nnow)
      IF (ltime) CALL get_timings (i_add_computations, ntstep, dt, izerror)

#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
      IF (luse_rttov) THEN
        CALL organize_satellites ('compute', izerror, yzerrmsg)
        IF (izerror /= 0_iintegers) THEN
          CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,             &
                                       'organize_satellites: compute')
        ENDIF
      ENDIF
#endif

      CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
      IF (ltime) CALL get_timings (i_add_computations, ntstep, dt, izerror)

      CALL organize_data ('result', ntstep, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                     'result: input-init')
      ENDIF
    ENDIF

    ! calculate surface pressure
    CALL calps ( ps(:,:,nnow),    pp(:,:,ke,nnow), t  (:,:,ke,nnow),   &
                 qv_now(:,:,ke) , qc_now(:,:,ke) , qrs(:,:,ke),        &
                 rho0(:,:,ke), p0(:,:,ke), dp0(:,:,ke),                &
                 ie, je, rvd_m_o, r_d, 1, ie, 1, je)
    CALL calps ( ps(:,:,nnew),    pp(:,:,ke,nnew), t  (:,:,ke,nnew),   &
                 qv_new(:,:,ke) , qc_new(:,:,ke) , qrs(:,:,ke),        &
                 rho0(:,:,ke), p0(:,:,ke), dp0(:,:,ke),                &
                 ie, je, rvd_m_o, r_d, 1, ie, 1, je)

    ! double dt after first step
    IF (ntstep == 0 .AND. .NOT.l2tls) THEN
      dt = 2.0_wp * dt
    ENDIF

    ! accumulate weighted fields in work arrays
    hfu  (:,:,:) = hfu (:,:,:) + hdfi(ntstep+1) * u (:,:,:,nnew)
    hfv  (:,:,:) = hfv (:,:,:) + hdfi(ntstep+1) * v (:,:,:,nnew)
    hfw  (:,:,:) = hfw (:,:,:) + hdfi(ntstep+1) * w (:,:,:,nnew)
    hft  (:,:,:) = hft (:,:,:) + hdfi(ntstep+1) * t (:,:,:,nnew)
    hfpp (:,:,:) = hfpp(:,:,:) + hdfi(ntstep+1) * pp(:,:,:,nnew)
    hfqv (:,:,:) = hfqv(:,:,:) + hdfi(ntstep+1) * qv_new(:,:,:)

    IF (itype_dfi_qx == 0) THEN

      ! filter all other qx-variables
      hfqc (:,:,:) = hfqc(:,:,:) + hdfi(ntstep+1) * qc_new(:,:,:)
      IF (ASSOCIATED(qi_new)) THEN
        hfqi (:,:,:) = hfqi(:,:,:) + hdfi(ntstep+1) * qi_new(:,:,:)
      ENDIF
      IF (ASSOCIATED(qr_new)) THEN
        hfqr (:,:,:) = hfqr(:,:,:) + hdfi(ntstep+1) * qr_new(:,:,:)
      ENDIF
      IF (ASSOCIATED(qs_new)) THEN
        hfqs (:,:,:) = hfqs(:,:,:) + hdfi(ntstep+1) * qs_new(:,:,:)
      ENDIF
      IF (ASSOCIATED(qg_new)) THEN
        hfqg (:,:,:) = hfqg(:,:,:) + hdfi(ntstep+1) * qg_new(:,:,:)
      ENDIF

    ELSE

      ! save values for all other qx-variables from nhalf time step
      IF ( ntstep == nhalf ) THEN
        IF(izdebug > 1) PRINT*,'DFI: Memorize qx var. for nhalf step ', ntstep

        IF (itype_dfi_qx == 2) THEN
          DO k = 1, ke
            DO j = 1, je
              DO i = 1, ie
                zrh(i,j,k) = MIN(1._wp,MAX(0._wp,(qv_new(i,j,k))/ &
!                zrh(i,j,k) = MIN(1._wp,MAX(0._wp,(qv_new(i,j,k)+qc_new(i,j,k))/ &
                 fqvs ( fpvsw ( t(i,j,k,nnew) ) , p0(i,j,k) + pp(i,j,k,nnew)  )))
!                       (p0(i,j,k) + pp(i,j,k,nnew) - O_m_rdv* &
!                       fpvsw(t (i,j,k,nnew)))/(Rdv*fpvsw(t (i,j,k,nnew))))) 
                ! Memorize rh for saturated grid points
                IF(.NOT.((zrh(i,j,k) < 1._wp .AND. (qc_new(i,j,k) > 1.E-6_wp &
                  .OR.qi_new(i,j,k) > 1.E-6_wp)) .OR. zrh(i,j,k) == 1._wp)) THEN
!                  IF( .NOT.zrh(i,j,k) == 1._wp) THEN
                  zrh(i,j,k)  = undef
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        zqc(:,:,:) = qc_new(:,:,:)
        IF (ASSOCIATED(qi_new)) zqi(:,:,:) = qi_new(:,:,:)
        IF (ASSOCIATED(qr_new)) zqr(:,:,:) = qr_new(:,:,:)
        IF (ASSOCIATED(qs_new)) zqs(:,:,:) = qs_new(:,:,:)
        IF (ASSOCIATED(qg_new)) zqg(:,:,:) = qg_new(:,:,:)
      ENDIF
    ENDIF

    IF (itype_dfi_soil == 1) THEN
      IF ( ntstep == nhalf ) THEN
        IF(izdebug > 1) PRINT*,'DFI: Memorize soil var. for nhalf step ', ntstep
        ALLOCATE (zts          (ie,je), STAT=izerror)
        IF (izerror == 0) THEN
          zts(:,:)  = t_s(:,:,nnew)
        ENDIF

        ALLOCATE (zqv_s        (ie,je), STAT=izerror)
        IF (izerror == 0) THEN
          zqv_s(:,:) = qv_s(:,:,nnew)
        ENDIF

        ALLOCATE (zw_snow      (ie,je), STAT=izerror)
        IF (izerror == 0) THEN
          zw_snow(:,:) = w_snow(:,:,nnew)
        ENDIF

        ALLOCATE (zt_snow      (ie,je), STAT=izerror)
        IF (izerror == 0) THEN
          zt_snow(:,:) = t_snow(:,:,nnew)
        ENDIF

        IF (lsoil) THEN
          IF (lmulti_layer) THEN
            ALLOCATE (zt_so    (ie,je,0:ke_soil+1), STAT=izerror)
            IF (izerror == 0) THEN
              zt_so(:,:,:)  = t_so(:,:,:,nnew)
            ENDIF

            ALLOCATE (zw_so    (ie,je,ke_soil+1), STAT=izerror)
            IF (izerror == 0) THEN
              zw_so(:,:,:)  = w_so(:,:,:,nnew)
            ENDIF

            ALLOCATE (zfreshsnw(ie,je), STAT=izerror)
            IF (izerror == 0) THEN
              zfreshsnw(:,:)  = freshsnow(:,:)
            ENDIF

            IF (lana_rho_snow) THEN
              ALLOCATE (zrhosnw(ie,je), STAT=izerror)
              IF (izerror == 0) THEN
                zrhosnw(:,:)  = rho_snow(:,:,nnew)
              ENDIF
            ENDIF
          ELSE
            ALLOCATE (zt_m     (ie,je), STAT=izerror)
            IF (izerror == 0) THEN
              zt_m(:,:)  = t_m(:,:,nnew)
            ENDIF

            ALLOCATE (zw_g1    (ie,je), STAT=izerror)
            IF (izerror == 0) THEN
              zw_g1(:,:)  = w_g1(:,:,nnew)
            ENDIF

            ALLOCATE (zw_g2    (ie,je), STAT=izerror)
            IF (izerror == 0) THEN
              zw_g2(:,:)  = w_g2(:,:,nnew)
            ENDIF

            IF ( nlgw == 3 ) THEN
              ALLOCATE (zw_g3  (ie,je), STAT=izerror)
              IF (izerror == 0) THEN
                zw_g3(:,:)  = w_g3(:,:,nnew)
              ENDIF
            ENDIF
          ENDIF
        ENDIF

        IF (llake) THEN
          ALLOCATE (ztmnw_lk   (ie,je), STAT=izerror)
          IF (izerror == 0) THEN
            ztmnw_lk(:,:)  = t_mnw_lk(:,:,nnew)
          ENDIF

          ALLOCATE (ztwml_lk   (ie,je), STAT=izerror)
          IF (izerror == 0) THEN
            ztwml_lk(:,:)  = t_wml_lk(:,:,nnew)
          ENDIF

          ALLOCATE (ztbot_lk   (ie,je), STAT=izerror)
          IF (izerror == 0) THEN
            ztbot_lk(:,:)  = t_bot_lk(:,:,nnew)
          ENDIF

          ALLOCATE (zct_lk     (ie,je), STAT=izerror)
          IF (izerror == 0) THEN
            zct_lk(:,:)  = c_t_lk(:,:,nnew)
          ENDIF

          ALLOCATE (zhml_lk    (ie,je), STAT=izerror)
          IF (izerror == 0) THEN
            zhml_lk(:,:)  = h_ml_lk(:,:,nnew)
          ENDIF
        ENDIF

        IF (lseaice .OR. llake) THEN
          ALLOCATE (ztice      (ie,je), STAT=izerror)
          IF (izerror == 0) THEN
            ztice(:,:)  = t_ice(:,:,nnew)
          ENDIF

          ALLOCATE (zhice      (ie,je), STAT=izerror)
          IF (izerror == 0) THEN
            zhice(:,:)  = h_ice(:,:,nnew)
          ENDIF
        ENDIF

        IF (izerror /= 0) THEN
          ierror = 10
          yerrmsg = 'Problem allocating space for ndfi=1'
          RETURN
        ENDIF
      ENDIF

    ENDIF

    ! cyclic change of time levels
    IF ( .NOT.l2tls ) THEN
      nsp  = nold
      nold = nnow
      nnow = nnew
      nnew = nsp
      nx0  = nnow
    ELSE
      nnow = 3-nnow
      nnew = 3-nnew
      nx0  = nnow
    END IF

    CALL update_trcr_pointers()

    ! nx0:  time level to compute mean values for soil variables later on

  ENDDO  ! diabatic forward integration

  ! copy filtered fields to prognostic fields as initial values 
  ! (time index 1 and 2)
  DO nx = 1, 2
    CALL trcr_get(izerror, idt_qv, ptr_tlev = nx, ptr = qv_nx)
    IF (izerror /= 0_iintegers) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get(izerror, idt_qc, ptr_tlev = nx, ptr = qc_nx)
    IF (izerror /= 0_iintegers) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get(izerror, idt_qi, ptr_tlev = nx, ptr = qi_nx)
    IF (izerror /= 0_iintegers .AND. izerror /= T_ERR_NOTFOUND) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get(izerror, idt_qr, ptr_tlev = nx, ptr = qr_nx)
    IF (izerror /= 0_iintegers .AND. izerror /= T_ERR_NOTFOUND) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get(izerror, idt_qs, ptr_tlev = nx, ptr = qs_nx)
    IF (izerror /= 0_iintegers .AND. izerror /= T_ERR_NOTFOUND) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get(izerror, idt_qg, ptr_tlev = nx, ptr = qg_nx)
    IF (izerror /= 0_iintegers .AND. izerror /= T_ERR_NOTFOUND) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF

    u  (:,:,:,nx) = hfu (:,:,:)
    v  (:,:,:,nx) = hfv (:,:,:)
    w  (:,:,:,nx) = hfw (:,:,:)
    t  (:,:,:,nx) = hft (:,:,:)
    pp (:,:,:,nx) = hfpp(:,:,:)
    qv_nx (:,:,:) = hfqv(:,:,:)

    IF (itype_dfi_qx == 0) THEN
      qc_nx (:,:,:) = hfqc(:,:,:)
      IF (ASSOCIATED(qi_nx)) qi_nx(:,:,:) = hfqi(:,:,:)
      IF (ASSOCIATED(qr_nx)) qr_nx(:,:,:) = hfqr(:,:,:)
      IF (ASSOCIATED(qs_nx)) qs_nx(:,:,:) = hfqs(:,:,:)
      IF (ASSOCIATED(qg_nx)) qg_nx(:,:,:) = hfqg(:,:,:)
    ELSE
      IF (itype_dfi_qx == 2) THEN
        DO k = 1, ke
          DO j = 1, je
            DO i = 1, ie
              IF(zrh(i,j,k) /= undef) THEN
                qv_nx (i,j,k) = MAX(qv_nx(i,j,k),zrh(i,j,k)* &
                 fqvs ( fpvsw ( t(i,j,k,nx) ) , p0(i,j,k) + pp(i,j,k,nx)  ))
!                qv_nx (i,j,k) = MAX(qv_nx(i,j,k),zrh(i,j,k)*Rdv*fpvsw(t(i,j,k,nx)) &
!                   /(p0(i,j,k) + pp(i,j,k,nx)  - O_m_rdv*fpvsw(t(i,j,k,nx))))
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        IF(nx == 2) DEALLOCATE(zrh)
      ENDIF
      qc_nx (:,:,:) = zqc(:,:,:)
      IF(nx == 2) DEALLOCATE(zqc)
      IF (ASSOCIATED(qi_nx)) THEN
        qi_nx(:,:,:) = zqi(:,:,:)
        IF(nx == 2) DEALLOCATE(zqi)
      ENDIF
      IF (ASSOCIATED(qr_nx)) THEN
        qr_nx(:,:,:) = zqr(:,:,:)
        IF(nx == 2) DEALLOCATE(zqr)
      ENDIF
      IF (ASSOCIATED(qs_nx)) THEN
        qs_nx(:,:,:) = zqs(:,:,:)
        IF(nx == 2) DEALLOCATE(zqs)
      ENDIF
      IF (ASSOCIATED(qg_nx)) THEN
        qg_nx(:,:,:) = zqg(:,:,:)
        IF(nx == 2) DEALLOCATE(zqg)
      ENDIF
    ENDIF

    ! calculate surface pressure
    CALL calps ( ps(:,:,nx), pp(:,:,ke,nx), t(:,:,ke,nx), qv_nx(:,:,ke),     &
                 qc_nx(:,:,ke), qrs(:,:,ke), rho0(:,:,ke), p0(:,:,ke),    &
                 dp0(:,:,ke), ie, je, rvd_m_o, r_d, 1, ie, 1, je)
  ENDDO

  IF (ndfi == 2) THEN

    ! Reset boundary values for nbd1 by the inverse:
    !  psi (nbd1) = 0.5 * (psi(nbd1) + psi(nbd2))
    WHERE (t_bd (:,:,:,nbd2) /= undef)
      u_bd (:,:,:,nbd1) = 0.5_wp * (u_bd (:,:,:,nbd1) + u_bd (:,:,:,nbd2))
      v_bd (:,:,:,nbd1) = 0.5_wp * (v_bd (:,:,:,nbd1) + v_bd (:,:,:,nbd2))
      t_bd (:,:,:,nbd1) = 0.5_wp * (t_bd (:,:,:,nbd1) + t_bd (:,:,:,nbd2))
      pp_bd(:,:,:,nbd1) = 0.5_wp * (pp_bd(:,:,:,nbd1) + pp_bd(:,:,:,nbd2))
      qv_bd(:,:,:,nbd1) = qvb(:,:,:)
      qc_bd(:,:,:,nbd1) = qcb(:,:,:)
    ENDWHERE

    IF (ASSOCIATED(qi_now)) THEN
      IF (izlbc_qi == T_LBC_FILE) THEN
        WHERE (qi_bd (:,:,:,nbd2) /= undef)
          qi_bd(:,:,:,nbd1) = qib(:,:,:)
        ENDWHERE
      ENDIF
    ENDIF
    IF (ASSOCIATED(qr_now) .AND. ASSOCIATED(qs_now)) THEN
      IF (izlbc_qr == T_LBC_FILE .AND. izlbc_qs == T_LBC_FILE) THEN
        WHERE (qr_bd (:,:,:,nbd2) /= undef)
          qr_bd(:,:,:,nbd1) = qrb  (:,:,:)
          qs_bd(:,:,:,nbd1) = qsb  (:,:,:)
        END WHERE
      ENDIF
    ENDIF
    IF (ASSOCIATED(qg_now)) THEN
      IF (izlbc_qg == T_LBC_FILE) THEN
        WHERE (qg_bd (:,:,:,nbd2) /= undef)
          qg_bd(:,:,:,nbd1) = qgb  (:,:,:)
        END WHERE
      ENDIF
    ENDIF

  ENDIF

  ! reset settings
  IF (ndfi == 1) THEN

    ntstep = nzntstep + nhalf
    nstart = nzntstep + nhalf

    ! also reset n0meanval and nsat_next
    IF (n0meanval   < nstart) n0meanval   = nstart
    IF (nextmeanval < nstart) nextmeanval = n0meanval

    IF (n0gp        < nstart) n0gp        = nstart
    IF (nnextgp     < nstart) nnextgp     = nstart

#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
    IF (luse_rttov) THEN
      nsat_next = 1
      DO WHILE ( (nsat_steps(nsat_next) <  nstart) .AND. (nsat_steps(nsat_next) >= 0) )
        nsat_next = nsat_next + 1
      ENDDO
    ENDIF
#endif

    IF (itype_dfi_soil == 0) THEN
      ! Average soil variables to nhalf
      t_s(:,:,1) = 0.5_wp * (zts(:,:) + t_s(:,:,nx0))
      t_s(:,:,2) = t_s(:,:,1)
      DEALLOCATE (zts)

      qv_s(:,:,1) = 0.5_wp * (zqv_s(:,:) + qv_s(:,:,nx0))
      qv_s(:,:,2) = qv_s(:,:,1)
      DEALLOCATE (zqv_s)

      w_snow(:,:,1) = 0.5_wp * (zw_snow(:,:) + w_snow(:,:,nx0))
      w_snow(:,:,2) = w_snow(:,:,1)
      DEALLOCATE (zw_snow)

      t_snow(:,:,1) = 0.5_wp * (zt_snow(:,:) + t_snow(:,:,nx0))
      t_snow(:,:,2) = t_snow(:,:,1)
      DEALLOCATE (zt_snow)

      IF (lsoil) THEN
        IF (lmulti_layer) THEN
          t_so(:,:,:,1) = 0.5_wp * (zt_so(:,:,:) + t_so(:,:,:,nx0))
          t_so(:,:,:,2) = t_so(:,:,:,1)
          DEALLOCATE (zt_so)

          w_so(:,:,:,1) = 0.5_wp * (zw_so(:,:,:) + w_so(:,:,:,nx0))
          w_so(:,:,:,2) = w_so(:,:,:,1)
          DEALLOCATE (zw_so)

          freshsnow(:,:) = 0.5_wp * (zfreshsnw(:,:) + freshsnow(:,:))
          DEALLOCATE (zfreshsnw)

          IF (lana_rho_snow) THEN
            rho_snow(:,:,1) = 0.5_wp * (zrhosnw(:,:) + rho_snow(:,:,nx0))
            rho_snow(:,:,2) = rho_snow(:,:,1)
            DEALLOCATE (zrhosnw)
          ENDIF
        ELSE
          t_m(:,:,1) = 0.5_wp * (zt_m(:,:) + t_m(:,:,nx0))
          t_m(:,:,2) = t_m(:,:,1)
          DEALLOCATE (zt_m)

          w_g1(:,:,1) = 0.5_wp * (zw_g1(:,:) + w_g1(:,:,nx0))
          w_g1(:,:,2) = w_g1(:,:,1)
          DEALLOCATE (zw_g1)

          w_g2(:,:,1) = 0.5_wp * (zw_g2(:,:) + w_g2(:,:,nx0))
          w_g2(:,:,2) = w_g2(:,:,1)
          DEALLOCATE (zw_g2)

          IF ( nlgw == 3 ) THEN
            w_g3(:,:,1) = 0.5_wp * (zw_g3(:,:) + w_g3(:,:,nx0))
            w_g3(:,:,2) = w_g3(:,:,1)
            DEALLOCATE (zw_g3)
          ENDIF
        ENDIF
      ENDIF

      IF (llake) THEN
        t_mnw_lk(:,:,1) = 0.5_wp * (ztmnw_lk(:,:) + t_mnw_lk(:,:,nx0))
        t_mnw_lk(:,:,2) = t_mnw_lk(:,:,1)
        DEALLOCATE (ztmnw_lk)

        t_wml_lk(:,:,1) = 0.5_wp * (ztwml_lk(:,:) + t_wml_lk(:,:,nx0))
        t_wml_lk(:,:,2) = t_wml_lk(:,:,1)
        DEALLOCATE (ztwml_lk)

        t_bot_lk(:,:,1) = 0.5_wp * (ztbot_lk(:,:) + t_bot_lk(:,:,nx0))
        t_bot_lk(:,:,2) = t_bot_lk(:,:,1)
        DEALLOCATE (ztbot_lk)

        c_t_lk(:,:,1) = 0.5_wp * (zct_lk(:,:) + c_t_lk(:,:,nx0))
        c_t_lk(:,:,2) = c_t_lk(:,:,1)
        DEALLOCATE (zct_lk)

        h_ml_lk(:,:,1) = 0.5_wp * (zhml_lk(:,:) + h_ml_lk(:,:,nx0))
        h_ml_lk(:,:,2) = h_ml_lk(:,:,1)
        DEALLOCATE (zhml_lk)
      ENDIF

      IF (lseaice .OR. llake) THEN
        t_ice(:,:,1) = 0.5_wp * (ztice(:,:) + t_ice(:,:,nx0))
        t_ice(:,:,2) = t_ice(:,:,1)
        DEALLOCATE (ztice)

        h_ice(:,:,1) = 0.5_wp * (zhice(:,:) + h_ice(:,:,nx0))
        h_ice(:,:,2) = h_ice(:,:,1)
        DEALLOCATE (zhice)
      ENDIF

    ELSE
      t_s(:,:,1) = zts(:,:)
      t_s(:,:,2) = t_s(:,:,1)
      DEALLOCATE (zts)

      qv_s(:,:,1) = zqv_s(:,:)
      qv_s(:,:,2) = qv_s(:,:,1)
      DEALLOCATE (zqv_s)

      w_snow(:,:,1) = zw_snow(:,:)
      w_snow(:,:,2) = w_snow(:,:,1)
      DEALLOCATE (zw_snow)

      t_snow(:,:,1) = zt_snow(:,:)
      t_snow(:,:,2) = t_snow(:,:,1)
      DEALLOCATE (zt_snow)

      IF (lsoil) THEN
        IF (lmulti_layer) THEN
          t_so(:,:,:,1) = zt_so(:,:,:)
          t_so(:,:,:,2) = t_so(:,:,:,1)
          DEALLOCATE (zt_so)

          w_so(:,:,:,1) = zw_so(:,:,:)
          w_so(:,:,:,2) = w_so(:,:,:,1)
          DEALLOCATE (zw_so)

          freshsnow(:,:) = zfreshsnw(:,:)
          DEALLOCATE (zfreshsnw)

          IF (lana_rho_snow) THEN
            rho_snow(:,:,1) = zrhosnw(:,:)
            rho_snow(:,:,2) = rho_snow(:,:,1)
            DEALLOCATE (zrhosnw)
          ENDIF
        ELSE
          t_m(:,:,1) = zt_m(:,:)
          t_m(:,:,2) = t_m(:,:,1)
          DEALLOCATE (zt_m)

          w_g1(:,:,1) = zw_g1(:,:)
          w_g1(:,:,2) = w_g1(:,:,1)
          DEALLOCATE (zw_g1)

          w_g2(:,:,1) = zw_g2(:,:)
          w_g2(:,:,2) = w_g2(:,:,1)
          DEALLOCATE (zw_g2)

          IF ( nlgw == 3 ) THEN
            w_g3(:,:,1) = zw_g3(:,:)
            w_g3(:,:,2) = w_g3(:,:,1)
            DEALLOCATE (zw_g3)
          ENDIF
        ENDIF
      ENDIF

      IF (llake) THEN
        t_mnw_lk(:,:,1) = ztmnw_lk(:,:)
        t_mnw_lk(:,:,2) = t_mnw_lk(:,:,1)
        DEALLOCATE (ztmnw_lk)

        t_wml_lk(:,:,1) = ztwml_lk(:,:)
        t_wml_lk(:,:,2) = t_wml_lk(:,:,1)
        DEALLOCATE (ztwml_lk)

        t_bot_lk(:,:,1) = ztbot_lk(:,:)
        t_bot_lk(:,:,2) = t_bot_lk(:,:,1)
        DEALLOCATE (ztbot_lk)

        c_t_lk(:,:,1) = zct_lk(:,:)
        c_t_lk(:,:,2) = c_t_lk(:,:,1)
        DEALLOCATE (zct_lk)

        h_ml_lk(:,:,1) = zhml_lk(:,:)
        h_ml_lk(:,:,2) = h_ml_lk(:,:,1)
        DEALLOCATE (zhml_lk)
      ENDIF

      IF (lseaice .OR. llake) THEN
        t_ice(:,:,1) = ztice(:,:)
        t_ice(:,:,2) = t_ice(:,:,1)
        DEALLOCATE (ztice)

        h_ice(:,:,1) = zhice(:,:)
        h_ice(:,:,2) = h_ice(:,:,1)
        DEALLOCATE (zhice)
      ENDIF

    ENDIF

  ELSEIF (ndfi == 2) THEN

    ntstep = nzntstep

  ENDIF

  nincbound  = nzincbound
  nlastbound = nzlastbound
  dt         = zdt
  epsass     = zepsass
  lphys      = lzphys
  lcond      = lzcond
  lsoil      = lzsoil
  llake      = lzlake
  lseaice    = lseaice
  itype_gscp = iztype_gscp

  IF ( .NOT.l2tls ) THEN
    nold = 3
    nnow = 1
    nnew = 2
  ELSE
    nnow = 2
    nnew = 1
    nstop = nzstop
  ENDIF
  ntke       = 0

  CALL update_trcr_pointers()

!------------------------------------------------------------------------------
! Section 4: Cleanup
!------------------------------------------------------------------------------

  ! Reset diagnostic fields
  rain_gsp = 0.0_wp
  snow_gsp = 0.0_wp
  rain_con = 0.0_wp
  snow_con = 0.0_wp
  snow_melt= 0.0_wp
  asob_t   = 0.0_wp
  athb_t   = 0.0_wp
  asob_s   = 0.0_wp
  athb_s   = 0.0_wp
  apab_s   = 0.0_wp
  aumfl_s  = 0.0_wp
  avmfl_s  = 0.0_wp
  ashfl_s  = 0.0_wp
  alhfl_s  = 0.0_wp
  alhfl_bs = 0.0_wp
  alhfl_pl = 0.0_wp
  dursun   = 0.0_wp

  ! Deallocate hdfi
  DEALLOCATE (hdfi, cheby, time, time2, weight, weight2)

  ! Last time measurement for the digital filtering
  IF (ltime) CALL get_timings (i_add_computations, 0, dt, izerror)

!------------------------------------------------------------------------------
! Internal procedure
!------------------------------------------------------------------------------

CONTAINS

!==============================================================================
!+ Internal procedure in "dfi_initialization" for the retrieval of tracers
!------------------------------------------------------------------------------

SUBROUTINE update_trcr_pointers

  ! Retrieve the required microphysics tracers
  !------------------------------------------

  ! QV
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnew, ptr = qv_new, ptr_bd = qv_bd, ptr_tens = qv_tens)
  IF (izerror /= 0_iintegers) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnow, ptr = qv_now)
  IF (izerror /= 0_iintegers) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  ! QC
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nnew, ptr = qc_new, ptr_bd = qc_bd, ptr_tens = qc_tens)
  IF (izerror /= 0_iintegers) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nnow, ptr = qc_now)
  IF (izerror /= 0_iintegers) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  ! QI
  CALL trcr_get(izerror, idt_qi, ptr_tlev = nnew, ptr = qi_new, ptr_bd = qi_bd, ptr_tens = qi_tens)
  IF (izerror /= 0_iintegers .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qi, ptr_tlev = nnow, ptr = qi_now)
  IF (izerror /= 0_iintegers .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  ! QR
  CALL trcr_get(izerror, idt_qr, ptr_tlev = nnew, ptr = qr_new, ptr_bd = qr_bd)
  IF (izerror /= 0_iintegers .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qr, ptr_tlev = nnow, ptr = qr_now)
  IF (izerror /= 0_iintegers .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  ! QS
  CALL trcr_get(izerror, idt_qs, ptr_tlev = nnew, ptr = qs_new, ptr_bd = qs_bd)
  IF (izerror /= 0_iintegers .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qs, ptr_tlev = nnow, ptr = qs_now)
  IF (izerror /= 0_iintegers .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  ! QG
  CALL trcr_get(izerror, idt_qg, ptr_tlev = nnew, ptr = qg_new, ptr_bd = qg_bd)
  IF (izerror /= 0_iintegers .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qg, ptr_tlev = nnow, ptr = qg_now)
  IF (izerror /= 0_iintegers .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  ! additional timelevel for LF
  IF (.NOT. l2tls) THEN
    CALL trcr_get(izerror, idt_qv, ptr_tlev = nold, ptr = qv_old)
    IF (izerror /= 0_iintegers) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get(izerror, idt_qc, ptr_tlev = nold, ptr = qc_old)
    IF (izerror /= 0_iintegers) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get(izerror, idt_qi, ptr_tlev = nold, ptr = qi_old)
    IF (izerror /= 0_iintegers .AND. izerror /= T_ERR_NOTFOUND ) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get(izerror, idt_qr, ptr_tlev = nold, ptr = qr_old)
    IF (izerror /= 0_iintegers .AND. izerror /= T_ERR_NOTFOUND ) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get(izerror, idt_qs, ptr_tlev = nold, ptr = qs_old)
    IF (izerror /= 0_iintegers .AND. izerror /= T_ERR_NOTFOUND ) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get(izerror, idt_qg, ptr_tlev = nold, ptr = qg_old)
    IF (izerror /= 0_iintegers .AND. izerror /= T_ERR_NOTFOUND ) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE update_trcr_pointers

!==============================================================================
!==============================================================================
!+ Internal procedure in "dfi_initialization" for the input of NAMELIST inictl
!------------------------------------------------------------------------------

SUBROUTINE input_inictl (nuspecif, nuin, ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine organizes the input of the NAMELIST-group inictl. 
!   The group inictl contains variables for controlling initialization
!   (kind of filtering, backward and foreward time step)
!
! Method:
!   All variables are initialized with default values and then read in from
!   the file INPUT. The input values are checked for errors and for
!   consistency. If wrong input values are detected the program prints
!   an error message. The program is not stopped in this routine but an
!   error code is returned to the calling routine that aborts the program after
!   reading in all other namelists.
!   In parallel mode, the variables are distributed to all nodes with the
!   environment-routine distribute_values.   
!   Both, default and input values are written to the file YUSPECIF
!   (specification of the run).
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
  INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
    nuspecif,     & ! Unit number for protocolling the task
    nuin            ! Unit number for Namelist INPUT file

  INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
    ierrstat        ! error status variable

! Variables for default values
  INTEGER (KIND=iintegers)   ::       &
    ndfi_d, & ! indicator for kind of filtering
              ! = 0 : No filtering (corresponds to LDFI=.F.)
              ! = 1 : Launching (forward stage only)
              ! = 2 : Full two stage filtering (default)
    nfilt_d   !  indicator for method of filtering
              ! = 1 : Dolph-Chebyshev filter (default)

  INTEGER (KIND=iintegers)   ::       &
    itype_dfi_qx_d,   & ! treatment of qx variables in dfi (active only for ndfi=1)
    itype_dfi_soil_d    ! treatment of soil variables in dfi (active only for ndfi=1)

  REAL (KIND=wp)             ::       &
    tspan_d,      & ! Time-span (in seconds) for the adiabatic and diabatic
                    ! stages of the initialization
    dtbak_d,      & ! Time-step for the backcast filtering stage (sec)
    dtfwd_d,      & ! Time-step for the forecast filtering stage (sec)
    taus_d          ! Cuttoff period (in seconds) for the filter
                    ! (Taus is the stop-band edge for the Dolph filter)

! Other Variables
  INTEGER (KIND=iintegers)   ::       &
    ierr, nhalff, nhalfb, iz_err

  REAL (KIND=wp)             ::       &
    zpi    ! 

  CHARACTER(LEN=250)         :: iomsg_str

! Define the namelist group
  NAMELIST /inictl/ ndfi, nfilt, tspan, dtbak, dtfwd, taus,                  &
                    itype_dfi_qx, itype_dfi_soil

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE input_inictl
!------------------------------------------------------------------------------

iz_err = 0_iintegers

IF (my_cart_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 1: Initialize the default variables
!------------------------------------------------------------------------------

  ndfi_d           = 0_iintegers
  nfilt_d          = 1_iintegers
  itype_dfi_qx_d   = 0_iintegers
  itype_dfi_soil_d = 0_iintegers

  tspan_d          = 3600.0_wp
  dtbak_d          = 90.0_wp
  dtfwd_d          = 90.0_wp
  taus_d           = 3600.0_wp

!------------------------------------------------------------------------------
!- Section 2: Initialize variables with defaults
!------------------------------------------------------------------------------

  ndfi             = ndfi_d
  nfilt            = nfilt_d
  itype_dfi_qx     = itype_dfi_qx_d
  itype_dfi_soil   = itype_dfi_soil_d

  tspan            = tspan_d
  dtbak            = dtbak_d
  dtfwd            = dtfwd_d
  taus             = taus_d

!------------------------------------------------------------------------------
!- Section 3: Input of the namelist values
!------------------------------------------------------------------------------

  iomsg_str(:) = ' '
  READ (nuin, inictl, IOSTAT=iz_err, IOMSG=iomsg_str)

  IF (iz_err /= 0) WRITE (*,'(A,A)') 'Namelist-ERROR INICTL: ', TRIM(iomsg_str)
ENDIF

IF (nproc > 1) THEN
  ! distribute error status to all processors
  CALL distribute_values  (iz_err, 1, 0, imp_integers,  icomm_cart, ierr)
ENDIF

IF (iz_err /= 0) THEN
  ierrstat = -1
  RETURN
ENDIF

IF (my_cart_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 4: Check values for errors and consistency
!------------------------------------------------------------------------------

  ! Check, which kind of filtering has been chosen
  IF ((ndfi < 1) .OR. (ndfi > 2)) THEN
    IF (my_cart_id == 0) THEN
      PRINT *, ' ERROR  *** wrong initialization type ndfi = ', ndfi, ' ***'
      PRINT *, '        *** only ndfi = 1 / 2 is possible       ***'
    ENDIF
    RETURN
  ENDIF

  ! Check whether backward integration (ndfi=2) is used for the Runge-Kutta
  ! scheme: This cannot be done
  IF ( (ndfi==2) .AND. l2tls) THEN
    PRINT *,' ERROR    *** Runge-Kutta cannot be integrated backwards ***'
    PRINT *,'          *** only ndfi=1 possible with l2tls=.TRUE.     ***'
    ierrstat = 1002
  ENDIF

  IF (nfilt /= 1) THEN
    PRINT *,' WARNING  ***  only nfilt = 1 is implemented   *** '
    PRINT *,' WARNING  *** program continues with nfilt = 1 *** '
    nfilt = 1
  ENDIF

  ! Compute and check additional values
  nbak   = NINT ( tspan / dtbak )
  nfwd   = NINT ( tspan / dtfwd )
  nhalfb = nbak / 2
  nhalff = nfwd / 2
  IF (2*nhalfb /= nbak) THEN
    PRINT *,' ERROR    ***  nbak must be even in initialization *** '
    ierrstat = 1002
  ENDIF
  IF (2*nhalff /= nfwd) THEN
    PRINT *,' ERROR    ***  nfwd must be even in initialization *** '
    ierrstat = 1002
  ENDIF

  ! Check whether nincbound is larger than the time span for the filtering
  ! in the initialization:
  IF ( nincbound*dt < tspan ) THEN
    PRINT *,' ERROR    *** tspan must be less than ',nincbound*dt, ' *** '
    ierrstat = 1002
  ENDIF

  zpi    = 4.0_wp * ATAN (1.0_wp)
  thsbak = 2.0_wp * zpi * dtbak / taus
  thsfwd = 2.0_wp * zpi * dtfwd / taus

ENDIF

!------------------------------------------------------------------------------
!- Section 5: Distribute variables to all nodes
!------------------------------------------------------------------------------

IF (nproc > 1) THEN

  IF (my_cart_id == 0) THEN
    intbuf  ( 1) = ndfi
    intbuf  ( 2) = nfilt
    intbuf  ( 3) = nbak
    intbuf  ( 4) = nfwd
    intbuf  ( 5) = itype_dfi_qx
    intbuf  ( 6) = itype_dfi_soil
    realbuf ( 1) = tspan
    realbuf ( 2) = dtbak 
    realbuf ( 3) = dtfwd
    realbuf ( 4) = taus
    realbuf ( 5) = thsbak
    realbuf ( 6) = thsfwd
  ENDIF

  CALL distribute_values (intbuf, 6, 0, imp_integers, icomm_cart, ierr)
  CALL distribute_values (realbuf,6, 0, imp_reals,    icomm_cart, ierr)

  IF (my_cart_id /= 0) THEN
    ndfi           = intbuf  ( 1)
    nfilt          = intbuf  ( 2)
    nbak           = intbuf  ( 3)
    nfwd           = intbuf  ( 4)
    itype_dfi_qx   = intbuf  ( 5)
    itype_dfi_soil = intbuf  ( 6)
    tspan          = realbuf ( 1)
    dtbak          = realbuf ( 2)
    dtfwd          = realbuf ( 3)
    taus           = realbuf ( 4)
    thsbak         = realbuf ( 5)
    thsfwd         = realbuf ( 6)
  ENDIF

ENDIF

!------------------------------------------------------------------------------
!- Section 6: Output of the namelist variables and their default values
!------------------------------------------------------------------------------

IF (my_cart_id == 0) THEN

  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(A23)') '0     NAMELIST:  inictl'
  WRITE (nuspecif, '(A23)') '      -----------------'
  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(T7,A,T33,A,T51,A,T71,A)') 'Variable', 'Actual Value',   &
                                               'Default Value', 'Format'

  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                 'ndfi       ',ndfi,       ndfi_d     ,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                 'nfilt      ',nfilt,     nfilt_d     ,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                          'itype_dfi_qx ',itype_dfi_qx,itype_dfi_qx_d ,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                    'itype_dfi_soil ',itype_dfi_soil,itype_dfi_soil_d ,' I '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                                          'tspan', tspan, tspan_d     ,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                                          'dtbak', dtbak, dtbak_d     ,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                                          'dtfwd', dtfwd, dtfwd_d     ,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                                          'taus' , taus,  taus_d      ,' R '
  WRITE (nuspecif, '(A2)')  '  '

ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE input_inictl

!==============================================================================

SUBROUTINE exchange_leapfrog

  IF     (lzconv) THEN
    kzdims(1:24) =                                                         &
       (/ke,ke,ke,ke,ke1,ke1,ke,ke,ke,ke,ke,ke,1,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                  &
       (39+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, &
        ie, je, kzdims, jstartpar, jendpar, 2, nboundlines,                &
        my_cart_neigh, lperi_x, lperi_y, l2dim,                            &
        200+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,          &
        u (:,:,:,nnow), u (:,:,:,nnew), v (:,:,:,nnow), v (:,:,:,nnew),    &
        w (:,:,:,nnow), w (:,:,:,nnew), t (:,:,:,nnow), t (:,:,:,nnew),    &
        pp(:,:,:,nnow), pp(:,:,:,nnew), qrs(:,:,:)    , dqvdt(:,:,:)  ,    &
        qvsflx(:,:) )
  ELSE
    kzdims(1:24) =                                                         &
       (/ke,ke,ke,ke,ke1,ke1,ke,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                  &
       (36+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, &
        ie, je, kzdims, jstartpar, jendpar, 2, nboundlines,                &
        my_cart_neigh, lperi_x, lperi_y, l2dim,                            &
        200+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,          &
        u (:,:,:,nnow), u (:,:,:,nnew), v (:,:,:,nnow), v (:,:,:,nnew),    &
        w (:,:,:,nnow), w (:,:,:,nnew), t (:,:,:,nnow), t (:,:,:,nnew),    &
        pp(:,:,:,nnow), pp(:,:,:,nnew), qrs(:,:,:)    )
  ENDIF

  ! Exchange of humidity tracers according to itype_gscp
  IF     (itype_gscp >= 4) THEN
    kzdims(1:24) =                                                         &
       (/ke,ke,ke,ke,ke,ke,ke,ke,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                  &
       ( 0     , sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, &
        ie, je, kzdims, jstartpar, jendpar, 2, nboundlines,                &
        my_cart_neigh, lperi_x, lperi_y, l2dim,                            &
        250+nexch_tag, .FALSE.   , ncomm_type, izerror, yzerrmsg,          &
        qv_now(:,:,:), qv_new(:,:,:), qc_now(:,:,:), qc_new(:,:,:),        &
        qr_now(:,:,:), qr_new(:,:,:), qs_now(:,:,:), qs_new(:,:,:),        &
        qi_now(:,:,:), qi_new(:,:,:), qg_now(:,:,:), qg_new(:,:,:))
  ELSEIF (itype_gscp==3) THEN
    kzdims(1:24) =                                                         &
       (/ke,ke,ke,ke,ke,ke,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                  &
       ( 0     , sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, &
        ie, je, kzdims, jstartpar, jendpar, 2, nboundlines,                &
        my_cart_neigh, lperi_x, lperi_y, l2dim,                            &
        250+nexch_tag, .FALSE.   , ncomm_type, izerror, yzerrmsg,          &
        qv_now(:,:,:), qv_new(:,:,:), qc_now(:,:,:), qc_new(:,:,:),        &
        qr_now(:,:,:), qr_new(:,:,:), qs_now(:,:,:), qs_new(:,:,:),        &
        qi_now(:,:,:), qi_new(:,:,:))
  ELSEIF (itype_gscp==2) THEN
    kzdims(1:24) =                                                         &
       (/ke,ke,ke,ke,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                  &
       ( 0     , sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, &
        ie, je, kzdims, jstartpar, jendpar, 2, nboundlines,                &
        my_cart_neigh, lperi_x, lperi_y, l2dim,                            &
        250+nexch_tag, .FALSE.   , ncomm_type, izerror, yzerrmsg,          &
        qv_now(:,:,:), qv_new(:,:,:), qc_now(:,:,:), qc_new(:,:,:),        &
        qr_now(:,:,:), qr_new(:,:,:), qs_now(:,:,:), qs_new(:,:,:))
  ELSEIF (itype_gscp==1) THEN
    kzdims(1:24) =                                                         &
       (/ke,ke,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                  &
       ( 0     , sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, &
        ie, je, kzdims, jstartpar, jendpar, 2, nboundlines,                &
        my_cart_neigh, lperi_x, lperi_y, l2dim,                            &
        250+nexch_tag, .FALSE.   , ncomm_type, izerror, yzerrmsg,          &
        qv_now(:,:,:), qv_new(:,:,:), qc_now(:,:,:), qc_new(:,:,:),        &
        qr_now(:,:,:), qr_new(:,:,:))
  ENDIF

END SUBROUTINE exchange_leapfrog

!==============================================================================
!==============================================================================

SUBROUTINE exchange_runge_kutta

INTEGER(KIND=iintegers) :: zntke


  kzdims(1:24)=(/ke,ke,ke1,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
  CALL exchg_boundaries                                                &
    (50+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,&
     ie, je, kzdims, jstartpar, jendpar, nbl_exchg, nboundlines,       &
     my_cart_neigh, lperi_x, lperi_y, l2dim,                           &
     200+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,         &
     u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),   &
     pp(:,:,:,nnew), qrs(:,:,:) )

  IF ( lzconv ) THEN
    IF ( lprog_tke ) THEN
      IF (itype_turb /= 3 .OR. ntke == 0) THEN
        zntke = nnew
      ELSE
        zntke = ntke
      ENDIF
      kzdims(1:24)=(/ke,1,ke1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                  &
        (52+zntke, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, &
         ie, je, kzdims, jstartpar, jendpar, nbl_exchg, nboundlines,         &
         my_cart_neigh, lperi_x, lperi_y, l2dim,                             &
         400+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,           &
         dqvdt(:,:,:), qvsflx(:,:), tke(:,:,:,zntke) )
    ELSE
      kzdims(1:24)=(/ke,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                  &
        (8, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,        &
         ie, je, kzdims, jstartpar, jendpar, nbl_exchg, nboundlines,         &
         my_cart_neigh, lperi_x, lperi_y, l2dim,                             &
         400+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,           &
         dqvdt(:,:,:), qvsflx(:,:) )
    END IF
  ELSE
    IF ( lprog_tke ) THEN
      kzdims(1:24)=(/ke1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      IF (itype_turb /= 3 .OR. ntke == 0) THEN
        zntke = nnew
      ELSE
        zntke = ntke
      ENDIF
      CALL exchg_boundaries                                                  &
        (54+zntke, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, &
         ie, je, kzdims, jstartpar, jendpar, nbl_exchg, nboundlines,         &
         my_cart_neigh, lperi_x, lperi_y, l2dim,                             &
         400+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,           &
         tke(:,:,:,zntke) )
    END IF
  END IF


  ! Exchange of tracers
  IF     (itype_gscp >= 4) THEN
    kzdims(1:24)=(/ke,ke,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                &
      ( 0     , sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,&
       ie, je, kzdims, jstartpar, jendpar, nbl_exchg, nboundlines,       &
       my_cart_neigh, lperi_x, lperi_y, l2dim,                           &
       200+nexch_tag, .FALSE.   , ncomm_type, izerror, yzerrmsg,         &
       qv_new(:,:,:), qc_new(:,:,:), qr_new(:,:,:), qs_new(:,:,:),       &
       qi_new(:,:,:), qg_new(:,:,:))
  ELSEIF (itype_gscp == 3) THEN
    kzdims(1:24)=(/ke,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                &
      ( 0     , sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,&
       ie, je, kzdims, jstartpar, jendpar, nbl_exchg, nboundlines,       &
       my_cart_neigh, lperi_x, lperi_y, l2dim,                           &
       200+nexch_tag, .FALSE.   , ncomm_type, izerror, yzerrmsg,         &
       qv_new(:,:,:), qc_new(:,:,:), qr_new(:,:,:), qs_new(:,:,:),       &
       qi_new(:,:,:))
  ELSEIF (itype_gscp == 2) THEN
    kzdims(1:24)=(/ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                &
      ( 0     , sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,&
       ie, je, kzdims, jstartpar, jendpar, nbl_exchg, nboundlines,       &
       my_cart_neigh, lperi_x, lperi_y, l2dim,                           &
       300+nexch_tag, .FALSE.   , ncomm_type, izerror, yzerrmsg,         &
       qv_new(:,:,:), qc_new(:,:,:), qr_new(:,:,:), qs_new(:,:,:))
  ELSE ! kessler_pp:
    kzdims(1:24)=(/ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                &
      ( 0     , sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,&
       ie, je, kzdims, jstartpar, jendpar, nbl_exchg, nboundlines,       &
       my_cart_neigh, lperi_x, lperi_y, l2dim,                           &
       300+nexch_tag, .FALSE.   , ncomm_type, izerror, yzerrmsg,         &
       qv_new(:,:,:), qc_new(:,:,:), qr_new(:,:,:))
  ENDIF

END SUBROUTINE exchange_runge_kutta

!==============================================================================

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE dfi_initialization
