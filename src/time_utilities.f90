!+ Utility module for timing routines
!------------------------------------------------------------------------------

MODULE  time_utilities

!------------------------------------------------------------------------------
!
! Description:
!   This module provides service utilities for timing different parts of the
!   model. 
!
!   Routines (module procedures) currently contained:
!
!     - init_timings:
!       initialize the variables for time-measuring
!
!     - get_timings:
!       Gets the time used for a special part of the model
!
!     - collect_timings:
!       At the end of the program all timings are collected to one node and 
!       printed.
!
! Method:
!   get_timings has to be called for special parts of the program. The special
!   parts are characterized by strings, which are passed to get_timings. 
!   35 parts are preset (see table below). If a string is passed that is not
!   recognized by get_timings, an additional entry is created and the times for
!   that entry are also stored (alltogether up to 70 entries are possible).
!   Before get_timings is called the variables necessary for time measuring 
!   have to be initialized by calling init_timings.
!   At the end of the program collect_timings collects all measurements to
!   node 0 and prints them to disk.
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.34       1999/12/10 Ulrich Schaettler
!  Initial release
! 1.39       2000/05/03 Ulrich Schaettler
!  Eliminate model_step as private variable and introduce dt as subroutine
!  argument in get_timings.
! 3.7        2004/02/18 Ulrich Schaettler
!  Changed treatment of unit number for ascii output file
! 3.18       2006/03/03 Klaus Stephan
!  Add LHN components for time measuring
! 3.21       2006/12/04 Ulrich Schaettler
!  Adapted table for LMK timings
! V4_3         2008/02/25 Ulrich Schaettler
!  Corrected start of time measurement in case of restart runs;
!  Eliminated T3E Flop counting
! V4_4         2008/07/16 Ulrich Schaettler
!  Introduced itype_timing, to determine, how the timings shall be handled
! V4_9         2009/07/16 Ulrich Schaettler
!  Implement timings for COSMO_ART
!  Corrected computations of sums, which started an index too early
!  Corrected computation of ifirsthour and storing initialization and cleanup
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_18        2011/05/26 Ulrich Schaettler
!  Implemented additional timings for COSMO-ART (Christoph Knote)
! V4_23        2012/05/10 Ulrich Schaettler
!  Modified SR collect_timings to save memory: 
!    for itype_timing=2/4 the sum over all hours is now computed in the PEs; 
!    for itype_timing=1/3, only the hours are gathered to PE 0 to save memory;
!          and a loop over all hours has been implemented
! V4_25        2012/09/28 Ulrich Blahak, Carlos Osuna
!  Introduced new logical variable l_2mom in interface to init_timings
!  Implemented additional timings for 2-mom-microphysics
!  Add i_asynio_wait to time that compute PE has to wait for asyn IO PE
!   to finish writing all the pending data before it can shutdown (CO)
! V4_28        2013/07/12 Ulrich Schaettler
!  Implemented additional output timers for meta data handling:
!    i_meta_data_r, i_meta_data_w
! V5_1         2014-11-28 Ulrich Blahak, Oliver Fuhrer
!  Implemented additional timings for radar forward operator
!  Replaced ireals by wp (working precision) (OF)
!  Added timings for copying data to / from block structure (US)
! V5_2         2015-05-21 Ulrich Blahak
!  Added timer i_radar_barrier for MPI barrier waiting in the
!   radar forward operator
! V5_3         2015-10-09 Oliver Fuhrer
!  Added timers for cpp-dycore
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Declarations:

! Modules used:
USE kind_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    i4, i8       ! KIND-type parameter for standard/double integer variables

USE data_parallel,    ONLY :   my_cart_id, num_compute, icomm_cart,        &
                               imp_reals, nc_asyn_io, nprocx, nprocy

!==============================================================================

IMPLICIT NONE

!==============================================================================

! include statements
INCLUDE "mpif.h"

!==============================================================================

! Private Declarations
  REAL (KIND=wp),     ALLOCATABLE, PRIVATE  ::      &
    timings (:,:)    ! for storing the times for different parts of the program

  INTEGER, PUBLIC          ::      &
    i_dyn_computations    , i_horizontal_diffusion, i_horizontal_advection,  &
    i_add_tend_moisture   , i_complete_tendencies , i_slow_tendencies     ,  &
    i_relaxation          , i_spectr_nudging      , i_semi_implicit       ,  &
    i_fast_waves          , i_dyn_09              ,                          &
#ifdef CPP_DYCORE
    i_cppdycore_copyin    , i_cppdycore_step      , i_cppdycore_copyout   ,  &
    i_cppdycore_swap      ,                                                  &
#endif
    i_phy_computations    , i_precipitation       , i_radiation           ,  &
    i_turbulence          , i_convection          , i_soil_model          ,  &
    i_sso,     i_phy_17   , i_phy_18              , i_phy_19              ,  &
    i_nud_computations    , i_obs_processing      , i_local_info          ,  &
    i_spread_multi_level  , i_spread_upper_air    , i_spread_surface      ,  &
    i_spread_ps_pre       , i_nudging_corr        , i_nud_28, i_nud_29    ,  &
    i_lhn_computations    , i_lhn_obs_prep        , i_lhn_t_inc           ,  &
    i_lhn_q_inc           , i_lhn_relax           , i_lhn_search          ,  &
    i_lhn_36,  i_lhn_37   , i_lhn_38              , i_lhn_39              ,  &
    i_add_computations    , i_input               , i_read_data           ,  &
    i_computations_I      , i_inp_44              , i_output              ,  &
    i_computations_O      , i_write_data          , i_out_48              ,  &
    i_meta_data_r         , i_meta_data_w         ,                          &
    i_initializations     , i_cleanup             , i_asynio_wait         ,  &
    i_semi_impl_comm      , i_semi_impl_barrier   , i_dyn_comm_01         ,  &
    i_fast_waves_comm     , i_fast_waves_barrier  ,                          &
    i_global_communi_dyn  , i_barrier_globcom_dyn ,                          &
    i_communications_dyn  , i_barrier_waiting_dyn , i_dyn_comm_02         ,  &
    i_communications_phy  , i_barrier_waiting_phy , i_phy_comm_01         ,  &
    i_communications_nud  , i_barrier_waiting_nud , i_nud_comm_01         ,  &
    i_communications_lhn  , i_barrier_waiting_lhn , i_lhn_comm_01         ,  &
    i_distribute_data     , i_gather_data         , i_copyblocks          ,  &
    ! timing for radar simulator
    i_radarsim            , i_radar_ini           , i_radar_compgrid      ,  &
    i_radar_comm          , i_radar_comppolar     , i_radar_ongeom        ,  &
    i_radar_out           , i_radar_barrier       ,                          &
    !
    i_art_staub_1         , i_art_staub_2         , i_art_staub_3         ,  &
    i_art_staub_4         , i_art_staub_5         , i_art_vorher          ,  &
    ! CK 20101117 additional timings for ART
    i_art_chemie          , i_art_aerosol         ,                          &
    i_art_emissions_I     , i_art_emissions       ,                          &
    i_art_bounds_I        , i_art_bounds          ,                          &
    i_art_initializations , i_art                 ,                          &
    ! CK end
    i_2mom_start          , i_2mom_init           , i_2mom_todens         ,  &
    i_2mom_clouds         , i_2mom_tospecif       , i_2mom_sedi           ,  &
    i_2mom_cleanup        , i_2mom                , i_2mom_init_0         ,  &
    i_2mom_cgps


  LOGICAL                                   ::      &
    ltu_phys,      & ! whether physical parameterizations are computed
    ltu_dfi,       & ! whether the initialization was done
    ltu_storedfi,  & ! for storing the times of dfi
    ltu_useobs,    & ! whether nudging was done
    ltu_2tls,      & ! time integration by two timelevel RK-scheme
    ltu_semi_imp,  & ! semi-implicit or split-explicit scheme
    ltu_art,       & ! run COSMO_ART
    ltu_2mom,      & ! run 2-moment microphysics
    ltu_radarsim     ! run radar simulator / forward operator

  CHARACTER (LEN=22), ALLOCATABLE, PRIVATE  ::      &
    ytable  (:)

  CHARACTER (LEN= 8) :: yutiming = 'YUTIMING'

  INTEGER ::         &
    itu_timing,      & ! to determine type of timing
    nutiming,        & ! unit number for output file
    idim_table,      & ! dimension of ytable
    idim_table_base, & ! dimension of ytable
    ifirsthour,      & ! first forecast hour
    ilasthour          ! last  forecast hour
    
  ! this variable really has to be the standard integer
  INTEGER, PRIVATE                          ::      &
    icountsold      ! number of counts in the last call

!==============================================================================

CONTAINS

!==============================================================================
!------------------------------------------------------------------------------

SUBROUTINE init_timings (nstart, nstop, dt, itype_timing, ldfi, lphys,     &
                         luseobs, l2tls, lsemi_imp, l_cosmo_art, l_2mom,   &
                         l_radar, istatus)

!------------------------------------------------------------------------------
!
! Description:
!   Initializes the variables necessary for measuring the different parts
!   of the model.
!
!------------------------------------------------------------------------------
!
! Parameter List:
! ---------------

INTEGER                 , INTENT (IN)             ::  &
  nstart,       & ! first timestep of the model
  nstop,        & ! last timestep of the model
  itype_timing    ! determines how to handle the timing

REAL    (KIND=wp),        INTENT (IN)             ::  &
  dt              ! time step

LOGICAL,                  INTENT(IN)              ::  &
  lphys,      & ! whether physical parameterizations are computed
  ldfi,       & ! whether the initialization was done
  luseobs,    & ! whether nudging was done
  l2tls,      & ! time integration by two timelevel RK-scheme
  lsemi_imp,  & ! semi-implicit or split-explicit scheme
  l_cosmo_art,& ! run COSMO_ART
  l_2mom,     & ! run 2-moment microphysics
  l_radar       ! run radar forward operator

INTEGER                 , INTENT (OUT), OPTIONAL  ::  &
  istatus         ! optional argument for error value

! Local Variables:
! ----------------

INTEGER           :: ir, im    ! arguments to SYSTEM_CLOCK
INTEGER           :: istat, itabidx
INTEGER (KIND=i8) :: i1, i2
LOGICAL           :: lpres     ! if optional argument is present

!------------ End of header ---------------------------------------------------

! Begin subroutine init_timings

  ! Set a default value for idim_table_base
  idim_table_base = 80

  idim_table    = idim_table_base
#ifdef CPP_DYCORE
  idim_table    = idim_table + 10
#endif
  IF (l_cosmo_art) idim_table = idim_table + 20
  IF (l_2mom)      idim_table = idim_table + 10
  itu_timing    = itype_timing
  ltu_phys      = lphys
  ltu_dfi       = ldfi
  ltu_storedfi  = ldfi
  ltu_useobs    = luseobs
  ltu_2tls      = l2tls
  ltu_semi_imp  = lsemi_imp
  ltu_art       = l_cosmo_art
  ltu_2mom      = l_2mom
  ltu_radarsim  = l_radar

  lpres = PRESENT (istatus)
  IF (lpres) THEN
    istatus = 0
  ENDIF

  ! Check, whether system clock is present and initialize icountsold
  CALL SYSTEM_CLOCK ( COUNT=icountsold, COUNT_RATE=ir, COUNT_MAX=im )
  IF ( ir == 0 ) THEN
    ! no system clock present: set error value
    IF ( lpres ) THEN
      istatus = 1
    ENDIF
  ENDIF

  ! Allocate character variable for storing the table of different
  ! parts of the model and fill the first entries
  ALLOCATE (ytable(idim_table), STAT=istat)

  ! Timings for Computation
  ytable(:)  = '                      '
  ytable( 1) = 'Dyn. Computations     ';  i_dyn_computations     =  1
  ytable( 2) = '  Horizontal Diffusion';  i_horizontal_diffusion =  2
  ytable( 3) = '  Horizontal Advection';  i_horizontal_advection =  3
  IF (ltu_2tls) THEN
  ytable( 4) = '  Add. Tend+Moist     ';  i_add_tend_moisture    =  4
  ELSE
  ytable( 4) = '  Complete Tendencies ';  i_complete_tendencies  =  4
  ENDIF
  ytable( 5) = '  Slow Tendencies     ';  i_slow_tendencies      =  5
  IF (ltu_semi_imp) THEN
  ytable( 6) = '  Semi Implicit Scheme';  i_semi_implicit        =  6
  ELSE
  ytable( 6) = '  Fast Waves          ';  i_fast_waves           =  6
  ENDIF
  ytable( 7) = '  Relaxation          ';  i_relaxation           =  7
  ytable( 8) = '  Spectr. Nudging     ';  i_spectr_nudging       =  8

  ytable(10) = 'Phy. Computations     ';  i_phy_computations     = 10
  ytable(11) = '  Precipitation       ';  i_precipitation        = 11
  ytable(12) = '  Radiation           ';  i_radiation            = 12
  ytable(13) = '  Turbulence          ';  i_turbulence           = 13
  ytable(14) = '  Convection          ';  i_convection           = 14
  ytable(15) = '  Soil Model          ';  i_soil_model           = 15
  ytable(16) = '  Sub-grid Scale Oro. ';  i_sso                  = 16
  ytable(19) = 'Copy data for blocks  ';  i_copyblocks           = 19

  ytable(20) = 'Nudg. Computations    ';  i_nud_computations     = 20
  ytable(21) = '  Obs. Processing     ';  i_obs_processing       = 21
  ytable(22) = '  Local Info          ';  i_local_info           = 22
  ytable(23) = '  Spread. Multi-Level ';  i_spread_multi_level   = 23
  ytable(24) = '  Spread. Upper-Air   ';  i_spread_upper_air     = 24
  ytable(25) = '  Spread. Surface     ';  i_spread_surface       = 25
  ytable(26) = '  Spread. ps + pre    ';  i_spread_ps_pre        = 26
  ytable(27) = '  Nudging & Corr.     ';  i_nudging_corr         = 27

  ytable(30) = 'LHN Computations      ';  i_lhn_computations     = 30
  ytable(31) = '  lhn_obs_prep        ';  i_lhn_obs_prep         = 31
  ytable(32) = '  lhn_t_inc           ';  i_lhn_t_inc            = 32
  ytable(33) = '  lhn_q_inc           ';  i_lhn_q_inc            = 33
  ytable(34) = '  lhn_relax           ';  i_lhn_relax            = 34
  ytable(35) = '  lhn_search          ';  i_lhn_search           = 35

  ytable(40) = 'Add. Computations     ';  i_add_computations     = 40

  ytable(41) = 'Input                 ';  i_input                = 41
  ytable(42) = '  read data           ';  i_read_data            = 42
  ytable(43) = '  meta data           ';  i_meta_data_r          = 43
  ytable(44) = '  computations I      ';  i_computations_I       = 44

  ytable(45) = 'Output                ';  i_output               = 45
  ytable(46) = '  computations O      ';  i_computations_O       = 46
  ytable(47) = '  write data          ';  i_write_data           = 47
  ytable(48) = '  meta  data          ';  i_meta_data_w          = 48

  ytable(49) = 'Initializations       ';  i_initializations      = 49
  ytable(50) = 'Cleanup               ';  i_cleanup              = 50

  IF (ltu_radarsim) THEN
  ytable(51) = 'Radar Simulator       ';  i_radarsim             = 51
  ytable(52) = '  Init./const. geom.  ';  i_radar_ini            = 52
  ytable(53) = '  Grid point values   ';  i_radar_compgrid       = 53
  ytable(54) = '  MPI Communications  ';  i_radar_comm           = 54
  ytable(55) = '  Online beam propag. ';  i_radar_ongeom         = 55
  ytable(56) = '  Comp. on polar grid ';  i_radar_comppolar      = 56
  ytable(57) = '  Output              ';  i_radar_out            = 57
  ytable(58) = '  Barrier Waiting     ';  i_radar_barrier        = 58
  END IF

  ! Timings for Communication
  IF (ltu_semi_imp) THEN
  ytable(61) = 'Semi Impl. Comm.      ';  i_semi_impl_comm       = 61
  ytable(62) = 'Semi Impl. Barrier    ';  i_semi_impl_barrier    = 62
  ELSE
  ytable(61) = 'Fast Waves Comm.      ';  i_fast_waves_comm      = 61
  ytable(62) = 'Fast Waves Barrier    ';  i_fast_waves_barrier   = 62
  ENDIF
  ytable(63) = 'Communications  Dyn   ';  i_communications_dyn   = 63
  ytable(64) = 'Barrier Waiting Dyn   ';  i_barrier_waiting_dyn  = 64
  ytable(65) = 'Global Comm     Dyn   ';  i_global_communi_dyn   = 65
  ytable(66) = 'Barrier GlobCom Dyn   ';  i_barrier_globcom_dyn  = 66

  ytable(67) = 'Communications  Phy   ';  i_communications_phy   = 67
  ytable(68) = 'Barrier Waiting Phy   ';  i_barrier_waiting_phy  = 68

  ytable(70) = 'Communications  Nud   ';  i_communications_nud   = 70
  ytable(71) = 'Barrier Waiting Nud   ';  i_barrier_waiting_nud  = 71

  ytable(73) = 'Communications  LHN   ';  i_communications_lhn   = 73
  ytable(74) = 'Barrier Waiting LHN   ';  i_barrier_waiting_lhn  = 74

  ytable(76) = '  distribute data     ';  i_distribute_data      = 76
  ytable(77) = '  gather data         ';  i_gather_data          = 77
  ytable(78) = 'asynIO wait           ';  i_asynio_wait          = 78

  itabidx = idim_table_base

#ifdef CPP_DYCORE
  ytable(itabidx+ 1) = 'wrapper copy in       ';  i_cppdycore_copyin    = itabidx+ 1
  ytable(itabidx+ 2) = 'wrapper dycore step   ';  i_cppdycore_step      = itabidx+ 2
  ytable(itabidx+ 3) = 'wrapper copy out      ';  i_cppdycore_copyout   = itabidx+ 3
  ytable(itabidx+ 4) = 'wrapper copy swap     ';  i_cppdycore_swap      = itabidx+ 4
  itabidx = itabidx + 10
#endif

  IF (ltu_art) THEN
    ytable(itabidx+ 1) = 'ART                   ';  i_art                 =  itabidx+ 1
    ytable(itabidx+ 2) = '  initializations     ';  i_art_initializations =  itabidx+ 2
    ytable(itabidx+ 3) = '  dust 1              ';  i_art_staub_1         =  itabidx+ 3
    ytable(itabidx+ 4) = '  dust 2              ';  i_art_staub_2         =  itabidx+ 4
    ytable(itabidx+ 5) = '  dust 3              ';  i_art_staub_3         =  itabidx+ 5
    ytable(itabidx+ 6) = '  dust 4              ';  i_art_staub_4         =  itabidx+ 6
    ytable(itabidx+ 7) = '  dust 5              ';  i_art_staub_5         =  itabidx+ 7
    ytable(itabidx+ 8) = '  before              ';  i_art_vorher          =  itabidx+ 8
    ytable(itabidx+ 9) = '  chemistry           ';  i_art_chemie          =  itabidx+ 9
    ytable(itabidx+10) = '  aerosol             ';  i_art_aerosol         =  itabidx+10
    ytable(itabidx+11) = '  emissions I         ';  i_art_emissions_I     =  itabidx+11
    ytable(itabidx+12) = '  emissions           ';  i_art_emissions       =  itabidx+12
    ytable(itabidx+13) = '  boundaries I        ';  i_art_bounds_I        =  itabidx+13
    ytable(itabidx+14) = '  boundaries          ';  i_art_bounds          =  itabidx+14
    itabidx = itabidx + 20
  ENDIF

  IF (ltu_2mom) THEN
    ytable(itabidx+ 0) = 'dummy to start clock  ';  i_2mom_start         = itabidx+ 0
    ytable(itabidx+ 1) = '2-moment-microph.     ';  i_2mom               = itabidx+ 1
    ytable(itabidx+ 2) = '  read lookup table   ';  i_2mom_init_0        = itabidx+ 2
    ytable(itabidx+ 3) = '  initializations     ';  i_2mom_init          = itabidx+ 3
    ytable(itabidx+ 4) = '  grid point search   ';  i_2mom_cgps          = itabidx+ 4
    ytable(itabidx+ 5) = '  conv. to dens.      ';  i_2mom_todens        = itabidx+ 5
    ytable(itabidx+ 6) = '  microphys.          ';  i_2mom_clouds        = itabidx+ 6
    ytable(itabidx+ 7) = '  conv. to specif.    ';  i_2mom_tospecif      = itabidx+ 7
    ytable(itabidx+ 8) = '  sedimentation       ';  i_2mom_sedi          = itabidx+ 8
    ytable(itabidx+ 9) = '  cleanup             ';  i_2mom_cleanup       = itabidx+ 9
    itabidx = itabidx + 10
  END IF

  ! Allocate variables for time measuring
  i1 = NINT(nstart * dt, i8)
  i2 = NINT(nstop  * dt, i8)
  IF (i2 == 0) THEN
    IF (ldfi) THEN
      ifirsthour = 0
    ELSE
      ifirsthour = 1
    ENDIF
    ilasthour  = 1
  ELSE
    IF (ldfi) THEN
      ifirsthour = 0
    ELSE
      IF (i1 == 0) THEN
        ifirsthour = 1
      ELSE
        IF(MOD (i1, 3600_i8) == 0_i8) THEN
          ifirsthour = INT (i1 / 3600_i8)
        ELSE
          ifirsthour = INT (i1 / 3600_i8) + 1
        ENDIF
      ENDIF
    ENDIF
    IF(MOD (i2,3600_i8) == 0_i8) THEN
      ilasthour = INT (i2 / 3600_i8)
    ELSE
      ilasthour = INT (i2 / 3600_i8) + 1
    ENDIF
  ENDIF
  ALLOCATE ( timings (idim_table, ifirsthour:ilasthour+1), STAT=istat )
  timings (:,:) = 0.0_wp

END SUBROUTINE init_timings

!==============================================================================
!==============================================================================
!------------------------------------------------------------------------------

SUBROUTINE get_timings (ipartname, nstep, dt, istatus)

!------------------------------------------------------------------------------
!
! Description:
!   Gets the elapsed wall-clock time in seconds since the last call. The value
!   is stored in the array "timings" with an entry corresponding to the entry
!   of ypartname in the table of entries "ytable". If ypartname does not 
!   correspond to an entry of that table, an additional entry is created 
!   (up to 70 entries on the whole).
!   The optional parameter istatus will return the following values:
!      0: on successful completion
!      1: no system clock present
!      2: no more entry available
!
! Method:
!   The intrinsic function SYSTEM_CLOCK is used, that returns the number of
!   clock counts since some system dependent event in the past (e.g. midnight
!   for a 24-hour system clock). The difference of clock counts since the last
!   call is determined and converted into seconds. The variable "icountsold" 
!   (see below) has to be SAVEd for the next call.
!
!------------------------------------------------------------------------------
!
! Parameter List:
! ---------------

INTEGER                 , INTENT (IN)             ::  &
  ipartname,    & ! index of the entry
  nstep           ! actual time step

REAL    (KIND=wp),        INTENT (IN)             ::  &
  dt              ! length of time step

INTEGER                 , INTENT (OUT), OPTIONAL  ::  &
  istatus         ! optional argument for error value

! Local Variables:
! ----------------

REAL (KIND=wp)     :: realtimedif       ! timing in seconds
INTEGER            :: ihours

INTEGER (KIND=i8)  :: i1

INTEGER            :: icountsnew,     & ! number of counts in this call
                      ir, im            ! other arguments to SYSTEM_CLOCK

LOGICAL                  ::           &
  lpres             ! if optional argument is present

!------------ End of header ---------------------------------------------------

! Begin subroutine get_time

  lpres = PRESENT (istatus)
  IF (lpres) THEN
    istatus = 0
  ENDIF

  IF ( (ipartname == i_initializations) .AND. (nstep == 0) .AND. ltu_dfi) THEN
    ! digital filtering has ended
    ltu_storedfi = .FALSE.
  ENDIF

  ! Determination of ihours for time-measurement
  IF (ltu_storedfi) THEN
    ihours = 0
  ELSE
    i1     = NINT(nstep * dt, i8)
    IF (i1 == 0) THEN
      ihours = 1
    ELSE
      IF (MOD (i1,3600_i8) == 0_i8) THEN
        ihours = INT(i1 / 3600)
      ELSE
        ihours = INT(i1 / 3600) + 1
      ENDIF
    ENDIF
  ENDIF

  ! Get timings
  CALL SYSTEM_CLOCK ( COUNT=icountsnew, COUNT_RATE=ir, COUNT_MAX=im )

  IF ( ir == 0 ) THEN
    ! system clock is not present
    IF (lpres) THEN
      istatus = 1
    ENDIF
  ELSE
    ! convert the clock counts to seconds
    IF ( icountsnew >= icountsold ) THEN
      realtimedif = ( REAL (icountsnew - icountsold, wp) )      &
                    / REAL (ir,wp)
    ELSE
      realtimedif = REAL (im- (icountsold-icountsnew ), wp)     &
                    / REAL (ir, wp)
    ENDIF
    icountsold = icountsnew

    ! Store value in the appropriate entry of timings:
    IF ( (1 <= ipartname) .AND. (ipartname <= idim_table) ) THEN
      timings(ipartname, ihours) = timings(ipartname, ihours) + realtimedif
    ELSE
      print *, ' *** WARNING:  Unknown table entry for time measurement *** '
      istatus = 2
    ENDIF
  ENDIF

END SUBROUTINE get_timings

!==============================================================================
!==============================================================================
!+ Subroutine that collects the timings from all nodes and prints them
!------------------------------------------------------------------------------

SUBROUTINE collect_timings

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine collects the timings from all nodes (if running in 
!   parallel mode) and prints them to an ASCII file YUTIMING. The measured
!   real times for the computations (dynamics, physics, diagnostics) and
!   the communications (small time step, large time step, physics, other
!   communications) are given for every forecast hour.
!
! Method:
!
!------------------------------------------------------------------------------

! Local variables

  REAL (KIND=wp),     ALLOCATABLE       ::                    &
    ztotal_times(:,:)        ! for gathering the times from all nodes

  REAL (KIND=wp)                        ::                    &
    zavgtime(4, idim_table),                                  &
    ztotal(ifirsthour:ilasthour+1), zsumtimes

  INTEGER                    ::                               &
    istat, izmplcode, niostat, ntab,                          &
    npr, nfpr, nlpr, nzsendcount, nzrecvcount,                &
    nzroot, ihours, n, ientrysum, iz, izstartprint, izendprint

  CHARACTER (LEN=25) yroutine
  CHARACTER (LEN=75) yerrmsg

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE collect_timings
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------

  ! Initializations
  yroutine    = 'collect_timings'
  izmplcode   = 0
  nzroot      = 0
  nzsendcount = idim_table         !* (ilasthour-ifirsthour+2)
  nzrecvcount = nzsendcount
  ientrysum   = ilasthour+1
  zavgtime(:,:) = 0.0_wp

  ! Allocate the buffers
  ALLOCATE ( ztotal_times (idim_table, 0:num_compute-1),   STAT = istat )
  ztotal_times(:,:) = 0.0_wp

  ! Compute sum over all hours in all tasks
  DO ihours = ifirsthour, ilasthour
    DO ntab = 1, idim_table
      timings(ntab, ientrysum) = timings(ntab, ientrysum) + timings(ntab, ihours)
    ENDDO
  ENDDO

  IF (my_cart_id == 0) THEN
    ! Calculate sums for computations
    zsumtimes = 0.0_wp
    DO ihours = ifirsthour, ilasthour + 1
      ! Calculate the sum of all timings for this hour "ihours"
      ztotal(ihours) = 0.0_wp
      DO ntab = 1, idim_table
        ztotal(ihours) = ztotal(ihours) + timings(ntab, ihours)
      ENDDO
      IF (ihours == ifirsthour) THEN
        ! subtract time for initialization
        ztotal(ifirsthour) = ztotal(ifirsthour) - timings(49,ifirsthour)
      ENDIF
      IF (ihours == ilasthour) THEN
        ! subtract time for cleanup
        ztotal(ilasthour) = ztotal(ilasthour) - timings(50,ilasthour)
      ENDIF
      zsumtimes = zsumtimes + ztotal(ihours)
    ENDDO
  ENDIF
      
  DO ihours = ifirsthour, ilasthour + 1
  
    ! Dyn. Computations:  sum( 1- 9)
    DO iz = 2, 9
      timings( i_dyn_computations,ihours) = timings( i_dyn_computations,ihours) + &
                                            timings(                 iz,ihours)
    ENDDO

    IF (ltu_2mom) THEN
      ! 2-moment microphysics: sum(102-109 to 101)
      DO iz = i_2mom_init_0, i_2mom_cleanup
            timings( i_2mom, ihours) =                             &
                           timings( i_2mom         ,ihours) +      &
                           timings(              iz,ihours)
      ENDDO
      ! Set the "original" i_precipitation (=11) correctly:
      ! This has to be done before the summation of the physical computations
      timings( i_precipitation, ihours ) =                         &
                           timings( i_precipitation ,ihours) +     &
                           timings( i_2mom          ,ihours)
    ENDIF

    ! Phy. Computations:  sum(10-18)
    ! 19 is the time for the copying to / from block structure
    DO iz = 11, 18
      timings( i_phy_computations,ihours) = timings( i_phy_computations,ihours) + &
                                            timings(                 iz,ihours)
    ENDDO

    ! Nud. Computations:  sum(20-29)
    DO iz = 21, 29
      timings( i_nud_computations,ihours) = timings( i_nud_computations,ihours) + &
                                            timings(                 iz,ihours)
    ENDDO

    ! LHN Computations:   sum(30-39)
    DO iz = 31, 39
      timings( i_lhn_computations,ihours) = timings( i_lhn_computations,ihours) + &
                                            timings(                 iz,ihours)
    ENDDO

    ! Input:              sum(41-44)
    DO iz = 42, 44
      timings( i_input           ,ihours) = timings( i_input           ,ihours) + &
                                            timings(                 iz,ihours)
    ENDDO
    ! add timings for global communications
    timings(i_input,ihours)             = timings(i_input,ihours) + &
                                          timings(i_distribute_data,ihours)

    ! Output:             sum(45-48)
    DO iz = 46, 48
      timings( i_output          ,ihours) = timings( i_output          ,ihours) + &
                                            timings(                 iz,ihours)
    ENDDO
    ! add timings for global communications
    timings(i_output,ihours)            = timings(i_output,ihours) + &
                                          timings(i_gather_data,ihours)

    IF (ltu_radarsim) THEN
      ! Radar forward operator: sum(52-58)
      DO iz = i_radar_ini, i_radar_barrier
        timings( i_radarsim      ,ihours) = timings( i_radarsim        ,ihours) + &
                                            timings(                 iz,ihours)
      ENDDO
    END IF

    IF (ltu_art) THEN
      ! ART:             sum(81-95)
      DO iz = i_art_initializations, i_art_bounds
        timings( i_art           ,ihours) = timings( i_art             ,ihours) + &
                                            timings(                 iz,ihours)
      ENDDO
    ENDIF

  ENDDO

!------------------------------------------------------------------------------
!- Section 2: Print headers to the file YUTIMING
!------------------------------------------------------------------------------

  IF (my_cart_id == 0) THEN

    OPEN(nutiming, FILE=yutiming, FORM=  'FORMATTED', STATUS='UNKNOWN',  &
         IOSTAT=niostat)
    IF(niostat /= 0) THEN
      yerrmsg = ' ERROR    *** opening file YUTIMING failed *** '
    ENDIF

    ! Just prepare the output
    SELECT CASE (itu_timing)

    CASE (1)       ! detailed timings per task and hour

      ! Print a headline in file YUDEBUG
      WRITE (nutiming, '(A)') '  '
      WRITE (nutiming, '(A)') '    Detailed timings for the forecast per task:'
      WRITE (nutiming, '(A)') '    ==========================================='
      WRITE (nutiming, '(A)') '  '

      ! Print the information from all processes
      WRITE (nutiming, '(A,F9.2)')                                            &
        'Time for the setup of the model:    ',                               &
                                timings(i_initializations,ifirsthour)
      WRITE (nutiming, '(A,F9.2)')                                            &
        'Time for the cleanup of the model:  ',                               &
                                timings(i_cleanup,ilasthour)
      WRITE (nutiming, '(A)')  '  '

!     izstartprint = ifirsthour
!     izendprint   = ilasthour

    CASE (2)       ! detailed timings per task, sum over all hours

      ! Print a headline in file YUDEBUG
      WRITE (nutiming, '(A)') '  '
      WRITE (nutiming, '(A)') '    Sums of  timings for the forecast per task:'
      WRITE (nutiming, '(A)') '    ==========================================='
      WRITE (nutiming, '(A)') '  '

      ! Print the information from all processes
      WRITE (nutiming, '(A,F9.2)')                                            &
        'Time for the setup of the model:    ',                               &
                                timings(i_initializations,ifirsthour)
      WRITE (nutiming, '(A,F9.2)')                                            &
        'Time for the cleanup of the model:  ',                               &
                                timings(i_cleanup,ilasthour)
      WRITE (nutiming, '(A)')  '  '

!     izstartprint = ientrysum
!     izendprint   = ientrysum

    CASE (3)       ! mean values for all tasks and hours

      ! Print a headline in file YUDEBUG
      WRITE (nutiming, '(A)') '  '
      WRITE (nutiming, '(A)') '    Mean values of the timings for all tasks:'  
      WRITE (nutiming, '(A)') '    ========================================='
      WRITE (nutiming, '(A)') '  '
      WRITE (nutiming, '(A,I5,A,I5,A)')                                       &
              'Timings for ',nprocx,' x ',nprocy,':  '
      WRITE (nutiming, '(A)')  '  '

      ! Print the information from all processes
      WRITE (nutiming, '(A,F9.2)')                                            &
        'Time for the setup of the model:    ',                               &
                                timings(i_initializations,ifirsthour)
      WRITE (nutiming, '(A,F9.2)')                                            &
        'Time for the cleanup of the model:  ',                               &
                                timings(i_cleanup,ilasthour)
      WRITE (nutiming, '(A,F9.2)')                                            &
        'Time for the total model run:       ', ztotal(ientrysum)
      WRITE (nutiming, '(A)')  '  '

      WRITE (nutiming, '(A,A)') '                    ',           &
                                '         Min.         Avg.         Max.        Total'

!     izstartprint = ifirsthour
!     izendprint   = ilasthour

    CASE (4)       ! mean values for all tasks, sum over all hours 

      ! Print a headline in file YUDEBUG
      WRITE (nutiming, '(A)') '  '
      WRITE (nutiming, '(A)') '    Mean values of the timings for all tasks:'  
      WRITE (nutiming, '(A)') '    ========================================='
      WRITE (nutiming, '(A)') '  '
      WRITE (nutiming, '(A,I5,A,I5,A)')                                       &
              'Timings for ',nprocx,' x ',nprocy,':  '
      WRITE (nutiming, '(A)')  '  '

      ! Print the information from all processes
      WRITE (nutiming, '(A,F9.2)')                                            &
        'Time for the setup of the model:    ',                               &
                                timings(i_initializations,ifirsthour)
      WRITE (nutiming, '(A,F9.2)')                                            &
        'Time for the cleanup of the model:  ',                               &
                                timings(i_cleanup,ilasthour)
      WRITE (nutiming, '(A,F9.2)')                                            &
        'Time for the total model run:       ', ztotal(ientrysum)
      WRITE (nutiming, '(A)')  '  '

      WRITE (nutiming, '(A,A)') '                    ',           &
                                '         Min.         Avg.         Max.        Total'

!     izstartprint = ientrysum
!     izendprint   = ientrysum

    END SELECT     ! type of timing

  ENDIF

  ! All processors have to determine izstartprint, izendprint
  SELECT CASE (itu_timing)
  CASE (1,3)       ! print outs per hour
      izstartprint = ifirsthour
      izendprint   = ilasthour
  CASE (2,4)       ! sum over all hours
      izstartprint = ientrysum
      izendprint   = ientrysum
  END SELECT     ! type of timing

!------------------------------------------------------------------------------
!- Section 3:  Gather data from all PEs on hourly basis and do the printing
!------------------------------------------------------------------------------

  ! And now do the printing
  SELECT CASE (itu_timing)

  CASE (1,2)

    DO ihours = izstartprint, izendprint

      IF (num_compute > 1) THEN
        CALL MPI_GATHER    (timings(:,ihours), nzsendcount, imp_reals,         &
                            ztotal_times,      nzrecvcount, imp_reals,         &
                            nzroot, icomm_cart, izmplcode)
        IF (izmplcode /= 0) THEN
          print *, 'somethings wrong here'
          RETURN
        ENDIF
      ELSE
        ztotal_times(:,0) = timings(:,ihours)
      ENDIF

      IF (my_cart_id == 0) THEN
        ! in this case it must really be checked whether ihours = 0
        ! only then I want to have this string
        IF (ihours == 0) THEN
          WRITE (nutiming, '(A,F10.2)')                                      &
            'Time for the initialization in node 0:',ztotal(ihours)
        ELSEIF (ihours <= ilasthour) THEN
          WRITE (nutiming, '(A,I6,A,F10.2)')                                 &
            'Time for hour ',ihours,' in node 0:',ztotal(ihours)
        ENDIF
        WRITE (nutiming, '(A)')  '  '
  
        DO npr = 0 , num_compute-1 , 8
          nfpr = npr
          nlpr = MIN (nfpr+7 , num_compute-1)
 
          ! Timings for the dynamics
          DO ntab = 1, 8
            WRITE (nutiming, '(A,8F10.2)') ytable(                 ntab)(1:22), &
                                    (ztotal_times(                 ntab,n),n=nfpr,nlpr)
          ENDDO
          WRITE   (nutiming, '(A,8F10.2)') ytable(    i_fast_waves_comm)(1:22), &
                                    (ztotal_times(    i_fast_waves_comm,n),n=nfpr,nlpr)
          WRITE   (nutiming, '(A,8F10.2)') ytable( i_fast_waves_barrier)(1:22), &
                                    (ztotal_times( i_fast_waves_barrier,n),n=nfpr,nlpr)
          WRITE   (nutiming, '(A,8F10.2)') ytable( i_communications_dyn)(1:22), &
                                    (ztotal_times( i_communications_dyn,n),n=nfpr,nlpr)
          WRITE   (nutiming, '(A,8F10.2)') ytable(i_barrier_waiting_dyn)(1:22), &
                                    (ztotal_times(i_barrier_waiting_dyn,n),n=nfpr,nlpr)
          WRITE   (nutiming, '(A,8F10.2)') ytable( i_global_communi_dyn)(1:22), &
                                    (ztotal_times( i_global_communi_dyn,n),n=nfpr,nlpr)
          WRITE   (nutiming, '(A,8F10.2)') ytable(i_barrier_globcom_dyn)(1:22), &
                                    (ztotal_times(i_barrier_globcom_dyn,n),n=nfpr,nlpr)
          WRITE   (nutiming, '(A)')  '  '

#ifdef CPP_DYCORE
          ! Timings for the wrapper
          WRITE (nutiming, '(A,8F10.2)') ytable(  i_cppdycore_copyin)(1:22), &
                                    (ztotal_times(i_cppdycore_copyin,n),n=nfpr,nlpr)
          WRITE (nutiming, '(A,8F10.2)') ytable(    i_cppdycore_step)(1:22), &
                                (ztotal_times(      i_cppdycore_step,n),n=nfpr,nlpr)
          WRITE (nutiming, '(A,8F10.2)') ytable( i_cppdycore_copyout)(1:22), &
                                   (ztotal_times(i_cppdycore_copyout,n),n=nfpr,nlpr)
          WRITE (nutiming, '(A,8F10.2)') ytable(    i_cppdycore_swap)(1:22), &
                                      (ztotal_times(i_cppdycore_swap,n),n=nfpr,nlpr)
#endif

          ! Timings for the physics
          IF (ltu_phys) THEN
            DO ntab = 10, 16
            WRITE (nutiming, '(A,8F10.2)') ytable(                 ntab)(1:22), &
                                    (ztotal_times(                 ntab,n),n=nfpr,nlpr)
            ENDDO
            WRITE (nutiming, '(A,8F10.2)') ytable( i_copyblocks        )(1:22), &
                                    (ztotal_times( i_copyblocks        ,n),n=nfpr,nlpr)
            WRITE (nutiming, '(A,8F10.2)') ytable( i_communications_phy)(1:22), &
                                    (ztotal_times( i_communications_phy,n),n=nfpr,nlpr)
            WRITE (nutiming, '(A,8F10.2)') ytable(i_barrier_waiting_phy)(1:22), &
                                    (ztotal_times(i_barrier_waiting_phy,n),n=nfpr,nlpr)
            WRITE (nutiming, '(A)')  '  '
          ENDIF

          ! Timings for the nudging
          IF (ltu_useobs) THEN
            DO ntab = 20, 27
            WRITE (nutiming, '(A,8F10.2)') ytable(                 ntab)(1:22), &
                                    (ztotal_times(                 ntab,n),n=nfpr,nlpr)
            ENDDO
            WRITE (nutiming, '(A,8F10.2)') ytable( i_communications_nud)(1:22), &
                                    (ztotal_times( i_communications_nud,n),n=nfpr,nlpr)
            WRITE (nutiming, '(A,8F10.2)') ytable(i_barrier_waiting_nud)(1:22), &
                                    (ztotal_times(i_barrier_waiting_nud,n),n=nfpr,nlpr)
            WRITE (nutiming, '(A)')  '  '
          ENDIF

          ! Timings for the LHN
          IF (ltu_useobs) THEN
            DO ntab = 30, 35
            WRITE (nutiming, '(A,8F10.2)') ytable(                 ntab)(1:22), &
                                    (ztotal_times(                 ntab,n),n=nfpr,nlpr)
            ENDDO
            WRITE (nutiming, '(A,8F10.2)') ytable( i_communications_lhn)(1:22), &
                                    (ztotal_times( i_communications_lhn,n),n=nfpr,nlpr)
            WRITE (nutiming, '(A,8F10.2)') ytable(i_barrier_waiting_lhn)(1:22), &
                                    (ztotal_times(i_barrier_waiting_lhn,n),n=nfpr,nlpr)
            WRITE (nutiming, '(A)')  '  '
          ENDIF

          ! Timings for the additional computations
          WRITE   (nutiming, '(A,8F10.2)') ytable(   i_add_computations)(1:22), &
                                    (ztotal_times(   i_add_computations,n),n=nfpr,nlpr)

          ! Timings for the input
          WRITE   (nutiming, '(A,8F10.2)') ytable(              i_input)(1:22), &
                                    (ztotal_times(              i_input,n),n=nfpr,nlpr)
          WRITE   (nutiming, '(A,8F10.2)') ytable(          i_read_data)(1:22), &
                                    (ztotal_times(          i_read_data,n),n=nfpr,nlpr)
          WRITE   (nutiming, '(A,8F10.2)') ytable(        i_meta_data_r)(1:22), &
                                    (ztotal_times(        i_meta_data_r,n),n=nfpr,nlpr)
          WRITE   (nutiming, '(A,8F10.2)') ytable(     i_computations_I)(1:22), &
                                    (ztotal_times(     i_computations_I,n),n=nfpr,nlpr)
          WRITE   (nutiming, '(A,8F10.2)') ytable(    i_distribute_data)(1:22), &
                                    (ztotal_times(    i_distribute_data,n),n=nfpr,nlpr)
          WRITE   (nutiming, '(A)')  '  '

          ! Timings for the output
          WRITE   (nutiming, '(A,8F10.2)') ytable(             i_output)(1:22), &
                                    (ztotal_times(        i_output,n),n=nfpr,nlpr)
          WRITE   (nutiming, '(A,8F10.2)') ytable(     i_computations_O)(1:22), &
                                    (ztotal_times(i_computations_O,n),n=nfpr,nlpr)
          WRITE   (nutiming, '(A,8F10.2)') ytable(        i_meta_data_w)(1:22), &
                                    (ztotal_times(  i_meta_data_w,n),n=nfpr,nlpr)
          WRITE   (nutiming, '(A,8F10.2)') ytable(         i_write_data)(1:22), &
                                    (ztotal_times(    i_write_data,n),n=nfpr,nlpr)
          WRITE   (nutiming, '(A,8F10.2)') ytable(        i_gather_data)(1:22), &
                                    (ztotal_times(   i_gather_data,n),n=nfpr,nlpr)
          IF( nc_asyn_io > 0 ) THEN
            WRITE   (nutiming, '(A,8F10.2)') ytable(    i_asynio_wait)(1:22),   &
                                    (ztotal_times(   i_asynio_wait,n),n=nfpr,nlpr)
          ENDIF

          IF (ltu_radarsim) THEN
            ! Timings for the radar simulator
            WRITE   (nutiming, '(A,8F10.2)') ytable(       i_radarsim)(1:22),&
                           (ztotal_times(       i_radarsim,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(      i_radar_ini)(1:22),&
                           (ztotal_times(      i_radar_ini,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable( i_radar_compgrid)(1:22),&
                           (ztotal_times( i_radar_compgrid,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable( i_radar_comm)(1:22),&
                           (ztotal_times( i_radar_comm,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable( i_radar_ongeom)(1:22),&
                           (ztotal_times( i_radar_ongeom,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(i_radar_comppolar)(1:22),&
                           (ztotal_times(i_radar_comppolar,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(      i_radar_out)(1:22),&
                           (ztotal_times(      i_radar_out,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(      i_radar_barrier)(1:22),&
                           (ztotal_times(      i_radar_barrier,n),n=nfpr,nlpr)
          END IF

          ! Timings for ART
          IF (ltu_art) THEN
            WRITE   (nutiming, '(A)')  '  '
            WRITE   (nutiming, '(A,8F10.2)') ytable(                i_art)(1:22), &
                                            (ztotal_times(i_art,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(i_art_initializations)(1:22), &
                                            (ztotal_times(i_art_initializations,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(        i_art_staub_1)(1:22), &
                                            (ztotal_times(i_art_staub_1,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(        i_art_staub_2)(1:22), &
                                            (ztotal_times(i_art_staub_2,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(        i_art_staub_3)(1:22), &
                                            (ztotal_times(i_art_staub_3,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(        i_art_staub_4)(1:22), &
                                            (ztotal_times(i_art_staub_4,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(        i_art_staub_5)(1:22), &
                                            (ztotal_times(i_art_staub_5,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(         i_art_vorher)(1:22), &
                                            (ztotal_times(i_art_vorher,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(         i_art_chemie)(1:22), &
                                            (ztotal_times(i_art_chemie,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(        i_art_aerosol)(1:22), &
                                            (ztotal_times(i_art_aerosol,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(    i_art_emissions_I)(1:22), &
                                            (ztotal_times(i_art_emissions_I,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(      i_art_emissions)(1:22), &
                                            (ztotal_times(i_art_emissions,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(       i_art_bounds_I)(1:22), &
                                            (ztotal_times(i_art_bounds_I,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(         i_art_bounds)(1:22), &
                                            (ztotal_times(i_art_bounds,n),n=nfpr,nlpr)
          ENDIF

          ! Timings for 2-moment microphysics:
          IF (ltu_2mom) THEN
            WRITE   (nutiming, '(A)')  '  '
            WRITE   (nutiming, '(A,8F10.2)') ytable(               i_2mom)(1:22),&
                                             (ztotal_times(i_2mom,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(        i_2mom_init_0)(1:22),&
                                             (ztotal_times(i_2mom_init_0,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(          i_2mom_init)(1:22),&
                                             (ztotal_times(i_2mom_init,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(          i_2mom_cgps)(1:22),&
                                             (ztotal_times(i_2mom_cgps,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(        i_2mom_todens)(1:22),&
                                             (ztotal_times(i_2mom_todens,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(        i_2mom_clouds)(1:22),&
                                             (ztotal_times(i_2mom_clouds,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(      i_2mom_tospecif)(1:22),&
                                             (ztotal_times(i_2mom_tospecif,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(          i_2mom_sedi)(1:22),&
                                             (ztotal_times(i_2mom_sedi,n),n=nfpr,nlpr)
            WRITE   (nutiming, '(A,8F10.2)') ytable(       i_2mom_cleanup)(1:22),&
                                             (ztotal_times(i_2mom_cleanup,n),n=nfpr,nlpr)
          ENDIF

          WRITE   (nutiming, '(A)')  '  '
          WRITE   (nutiming, '(A)')  '  '
        ENDDO
      ENDIF
    ENDDO

  CASE (3,4)

    DO ihours = izstartprint, izendprint

      IF (num_compute > 1) THEN
        CALL MPI_GATHER    (timings(:,ihours), nzsendcount, imp_reals,         &
                            ztotal_times,      nzrecvcount, imp_reals,         &
                            nzroot, icomm_cart, izmplcode)
        IF (izmplcode /= 0) THEN
          print *, 'somethings wrong here'
          RETURN
        ENDIF
      ELSE
        ztotal_times(:,0) = timings(:,ihours)
      ENDIF

      IF (my_cart_id == 0) THEN
        IF (ihours == 0) THEN
          WRITE (nutiming, '(A,F10.2)')                                      &
            'Time for the initialization in node 0:',ztotal(ihours)
        ELSEIF (ihours <= ilasthour) THEN
          WRITE (nutiming, '(A,I6,A,F10.2)')                                 &
            'Time for hour ',ihours,' in node 0:',ztotal(ihours)
        ENDIF
        WRITE (nutiming, '(A)')  '  '
 
        DO ntab = 1, idim_table
          IF (LEN_TRIM(ytable(ntab)) > 0) THEN
            zavgtime(1,ntab) = MINVAL (ztotal_times(ntab,:))
            zavgtime(3,ntab) = MAXVAL (ztotal_times(ntab,:))
            zavgtime(4,ntab) = 0.0_wp
            DO n=0,num_compute-1
              zavgtime(4,ntab) = zavgtime(4,ntab) + ztotal_times(ntab,n)
            ENDDO
            zavgtime(2,ntab) = zavgtime(4,ntab) / num_compute
          ENDIF
        ENDDO

        ! Timings for the dynamics
        DO ntab = 1, 8
          WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                         zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
        ENDDO
        ntab = i_fast_waves_comm
        WRITE   (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                         zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
        ntab = i_fast_waves_barrier
        WRITE   (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                         zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
        ntab = i_communications_dyn
        WRITE   (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                         zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
        ntab = i_barrier_waiting_dyn
        WRITE   (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                         zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
        ntab = i_global_communi_dyn
        WRITE   (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                         zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
        ntab = i_barrier_globcom_dyn
        WRITE   (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                         zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
        WRITE   (nutiming, '(A)')  '  '

#ifdef CPP_DYCORE
        ! Timings for the wrapper
        ntab = i_cppdycore_copyin
        WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                       zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
        ntab = i_cppdycore_step
        WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                       zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
        ntab = i_cppdycore_copyout
        WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                       zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
        ntab = i_cppdycore_swap
        WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                       zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
#endif

        ! Timings for the physics
        IF (ltu_phys) THEN
          DO ntab = 10, 16
          WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                         zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
          ENDDO
          ntab = i_copyblocks
          WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                         zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
          ntab = i_communications_phy
          WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                         zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
          ntab = i_barrier_waiting_phy
          WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                         zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
          WRITE (nutiming, '(A)')  '  '
        ENDIF

        ! Timings for the nudging
        IF (ltu_useobs) THEN
          DO ntab = 20, 27
          WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                         zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
          ENDDO
          ntab = i_communications_nud
          WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                         zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
          ntab = i_barrier_waiting_nud
          WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                         zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
          WRITE (nutiming, '(A)')  '  '
        ENDIF

        ! Timings for the LHN
        IF (ltu_useobs) THEN
          DO ntab = 30, 35
          WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                         zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
          ENDDO
          ntab = i_communications_lhn
          WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                         zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
          ntab = i_barrier_waiting_lhn
          WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                         zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
          WRITE (nutiming, '(A)')  '  '
        ENDIF

        ! Timings for the additional computations
        ntab = i_add_computations
        WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                       zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)

        ! Timings for the input, output
        ntab = i_input
        WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                       zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
        ntab = i_read_data
        WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                       zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
        ntab = i_meta_data_r
        WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                       zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
        ntab = i_computations_I
        WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                       zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
        ntab = i_distribute_data
        WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                       zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)

        ! Timings for the output
        ntab = i_output
        WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                       zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
        ntab = i_computations_O
        WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                       zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
        ntab = i_meta_data_w
        WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                       zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
        ntab = i_write_data
        WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                       zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
        ntab = i_gather_data
        WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                       zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
        IF( nc_asyn_io > 0 ) THEN
          ntab = i_asynio_wait
          WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab),  &
                       zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
        ENDIF

        ! Timings for the radar simulator
        IF (ltu_radarsim) THEN
          WRITE (nutiming, '(A)')  '  '
          DO ntab = i_radarsim, i_radar_barrier
          WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                        zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
          END DO
        END IF

        IF (ltu_art) THEN
          WRITE (nutiming, '(A)')  '  '
          DO ntab = i_art, i_art_bounds
          WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                         zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
          ENDDO
        ENDIF

        IF (ltu_2mom) THEN
          WRITE (nutiming, '(A)')  '  '
          DO ntab = i_2mom, i_2mom_cleanup
            WRITE (nutiming, '(A,8F13.2)') ytable(ntab)(1:22), zavgtime(1,ntab), &
                           zavgtime(2,ntab), zavgtime(3,ntab), zavgtime(4,ntab)
          ENDDO
        ENDIF

        WRITE (nutiming, '(A)')  '  '
        WRITE (nutiming, '(A)')  '  '
      ENDIF
    ENDDO

  END SELECT     ! type of timing

!------------------------------------------------------------------------------
!- Section 4: Cleanup
!------------------------------------------------------------------------------

  IF (my_cart_id == 0) THEN
    CLOSE (nutiming, STATUS='KEEP')
  ENDIF

  DEALLOCATE ( ztotal_times, STAT=istat)
  DEALLOCATE ( timings,      STAT=istat)
  DEALLOCATE ( ytable,       STAT=istat)

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------
  
END SUBROUTINE collect_timings

!==============================================================================
!==============================================================================
!------------------------------------------------------------------------------
! End of the module
!------------------------------------------------------------------------------

END MODULE time_utilities
