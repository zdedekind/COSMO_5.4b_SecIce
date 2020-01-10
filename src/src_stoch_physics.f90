!+ Source module for stochastics physics utility routines
!==============================================================================

MODULE src_stoch_physics

!==============================================================================
!
! Description:
!   This module provides service utilities for stochastic perturbation of
!   physics tendencies (SPPT) calculations mainly based on Buizza et al (1999)
!   and Palmer et al. (2009).
!   Options for "smoothed" perturbations in time and space are also provided.
!
!     A random number rn is drawn from a uniform distribution or a Gaussian
!     distribution (lgauss_rn) for every point of an externally specified
!     horizontal coarse grid (box) having (ie_rn, je_rn) points and
!     (dlon_rn, dlat_rn) grid spacings.
!     The same random number is used for every point of the vertical grid
!     column (except if vertical tapering is applied). 
!
!     All the model grid points contained in the same coarse grid point have
!     the same random number (except if horizontal interpolation is applied). 
!
!     A new random number is drawn every hinc_rn hours.
!
!     The random number stream is defined by a date/time, the ensemble member
!     number, and an external seed number (nseed_rn) set by namelist.
!
!     In case of uniform distribution (0<rn<1) the magnitude of stochastics 
!     physics perturbation is defined as: 
!          ppt/pt=1+a*(1-2*rn)
!          where:
!           pt  : tendency
!           ppt : perturbed tendency
!           a   : max perturbation amplitude (range_rn)
!           rn  : random number (0<rn<1)         
!
!     In case of gaussian distribution (lgauss_rn=.TRUE.) the magnitude of
!       stochastics physics perturbation is defined as: 
!          ppt/pt=1+MIN(a, std*rn) if rn > 0
!          ppt/pt=1-MIN(a,-std*rn) if rn < 0
!       where:
!          std : perturbation standard deviation (stdv_rn)
!  
!     Thus, depending on a (range_rn), the relative size of the tendency 
!     perturbation can be:   (1-a) < ppt/pt < (1+a)
! 
!     The same random number is applied to perturb T, u, v, qv 
!     physical tendencies.
!
!     itype_vtaper_rn: type of random number tapering near surface/stratosphere
!                  itype_vtaper_rn=0 (no vert. tapering)
!                  itype_vtaper_rn=1 (prescribed: sfc and stratosph.)
!                  itype_vtaper_rn=2 (prescribed: only stratosph.)
!                  itype_vtaper_rn=3 (specified from namelist - vtaper_rn(1:ke))
!
!     If lhorint_rn=T the random numbers are horizontally interpolated on
!       the model grid
! 
!     If ltimeint_rn=T the random numbers are linearly interpolated in time 
!       every time step between the previous and the next drawn value (every
!       hinc_rn hours) 
!
!     If npattern_rn > 1 it is possible to specify more random number patterns  
!
!     imode_rn    0:  for each random number pattern, only one stream of
!                     random numbers is created and used for all time levels,
!                     and initial model date/time is used for the initial
!                     seed of the random number stream;
!                     the same offset of the coarse grid relative to the lower
!                     left corner of the COSMO grid is then used for all random
!                     number time steps (for a given pattern)
!                 1:  for each random number pattern, a new stream of random
!                     numbers is created for each random number time step, and
!                     for the seed of the random number stream, the date/time
!                     is used for which the respective random number field is
!                     valid;
!                     this option should be used for a data assimilation cycle
!                     and subsequent forecasts to ensure temporal correlations;
!                     the offset of the coarse grid relative to the lower left
!                     corner of the COSMO grid is then different for each
!                     random number time step
!
!     itype_qxpert_rn define which hum variables tend. are perturbed
!                 0:  Only qv tendencies are perturbed
!                 1:  qv,qc,qi tendencies are perturbed
!                 2:  qv,qc,qi,qr,qs,qg tendencies are perturbed
!
!     itype_qxlim_rn  type of reduction/removal of the perturbation 
!                     in case of negative (qv, qc, qi) or
!                     supersaturated (qv) values
!
! Routines (module procedures) currently contained:
!   - init_stoch_phys    : initializes rn generator, produces first rn field
!   - compute_stoch_phys : computes a 3-D random field for SPPT
!   - gen_rand_field     : produces a 3-D random field (for 1 time level)
!   - prep_init_rand_numb: prepares production of random number fields
!   - init_rn_generator  : initialise random number generator
!   - gen_rand_numb      : creates a set of random numbers
!   - int_rand_numb      : interpolates rn from coarse grid to COSMO grid
!   - set_seed_rand_numb : sets the seed of a rn stream
!   - apply_tqx_tend_adj : computes T, qv perturbation checking humidity limits
!
! Current Code Owner: CNMCA, Lucio Torrisi
!  phone:  +39 06 91293890
!  fax:    +39 06 91293889
!  email:  torrisi@meteoam.it
!
! History:
! Version      Date       Name
! ----------   ---------- ----
! V5_1         2014-11-28 Lucio Torrisi, Christoph Schraff
!  Initial release.
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:
!------------------------------------------------------------------------------

USE kind_parameters , ONLY :   &
  wp,        & ! KIND-type parameters for real variables in working precision
  dp           ! KIND-type parameters for double precision variables

!------------------------------------------------------------------------------

USE vgrid_refatm_utils, ONLY :   &
    refatm,         & ! reference pressure at sea level
    vcoord            ! sigma-coordinate refering to PMSL

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &
    dt,           & ! long time-step
    dt2,          & ! 2*long time-step

! horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------
    dlon,         & ! grid point distance in zonal direction (in degrees)
    dlat,         & ! grid point distance in meridional direction (in degrees)
!   startlon_tot, & ! transformed longitude of the lower left grid point
                    ! of the total domain (in degrees, E>0)
!   startlat_tot, & ! transformed latitude of the lower left grid point
                    ! of the total domain (in degrees, N>0)
    ie_tot,       & ! number of grid points in zonal direction
    je_tot,       & ! number of grid points in meridional direction
    ke_tot,       & ! number of grid points in vertical direction
    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction
!   zonal direction
    istart,       & ! start index for the forecast of w, t, qv, qc and pp
    iend,         & ! end index for the forecast of w, t, qv, qc and pp
!   istartu,      & ! start index for the forecast of u
!   iendu,        & ! end index for the forecast of u
!   istartv,      & ! start index for the forecast of v
!   iendv,        & ! end index for the forecast of v

!   meridional direction
    jstart,       & ! start index for the forecast of w, t, qv, qc and pp
    jend            ! end index for the forecast of w, t, qv, qc and pp
!   jstartu,      & ! start index for the forecast of u
!   jendu,        & ! end index for the forecast of u
!   jstartv,      & ! start index for the forecast of v
!   jendv           ! end index for the forecast of v

USE data_runcontrol,    ONLY :  &
    ntstep,        & ! actual time step
    nstart,        & ! first time step of the forecast
    l2tls,         & ! forecast with 2-TL integration scheme
    itype_calendar,& ! for specifying the calendar used
    iepsmem,       & ! ID of the member in the ensemble (ID >= 0) 
    iepstot,       & ! total number of EPS members

    lhorint_rn,    & ! random numbers (defined on a rn horiz. coarse grid)
                     !   horizontally interpolated on model grid (otherwise
                     !   model grid points contained in the same rn coarse grid
                     !   point have the same random number value)
    ltimeint_rn,   & ! random numbers (defined on a rn horiz. coarse grid) are
                     !   interpolated in time
    lgauss_rn,     & ! use a gaussian distribution of random numbers
                     !   (otherwise a uniform distribution is used)
    imode_rn,      & ! = 0: use only 1 stream of rn for all rn time steps
                     ! = 1: use a new stream of rn for every rn time step,
                     !      this enables temporal corrations in DA cycles
    npattern_rn,   & ! number of rn patterns with different scale/correlation
    hinc_rn,       & ! time increment (in [hrs]) for drawing a new field of
                     !   random numbers
    hstep_rn,      & ! time [   hrs   ] for which latest rn field is valid
    nseed_rn,      & ! external part of seed for random number generation
    nseed_rn2,     & ! external seed to generate rn at ntstep>0 (if imode_rn>0)
    dlat_rn,       & ! random number coarse grid point distance in
                     !   meridional direction (in degrees)
    dlon_rn,       & ! random number coarse grid point distance in
                     !   zonal direction (in degrees)
    stdv_rn,       & ! standard deviation of the gaussian distribution of 
                     !   random numbers
    range_rn,      & ! max magnitude of random numbers 
    rvtaper_rn,    & ! externally specified function for vertical tapering of
                     !  the random number
    itype_vtaper_rn,&! type of tapering near surface and in stratosphere
    n1_rn,         & ! / indices for permutation of the
    n2_rn,         & ! \ two random number time levels
    ie_rn,         & ! number of horiz. coarse grid points in zonal direction
    je_rn            ! number of horiz. coarse grid points in meridional direct.
                     !   where random numbers are defined

USE data_io,            ONLY :  &
    ydate_ini        ! start of the forecast 
                     ! yyyymmddhh (year, month, day, hour)

USE data_parallel,      ONLY :  &
    icomm_compute, & ! communicator for the group of compute PEs
    imp_reals,     & ! correct type for MPI
    imp_integers,  & ! correct type for MPI
    num_compute,   & ! number of compute PEs
    my_cart_id       ! rank of this subdomain in the cartesian communicator
                     ! model for MPI
   
USE data_constants  , ONLY :   &
    !   physical constants and related variables
    rdv,          & ! r_d / r_v
    !   constants for parametrizations
    b1,           & ! variables for computing the saturation vapour pressure
    b2w,          & ! over water (w) and ice (i)
    b3,           & !               -- " --
    b4w             !               -- " --
    
USE parallel_utilities,      ONLY :  &
    i_global        ,& ! function for global i-index from local i-index
    j_global        ,& ! function for global j-index from local j-index
    distribute_values  ! distributes values from 1 node to all the others

USE utilities,               ONLY :  &
    get_utc_date       ! actual date of the forecast in different forms

USE mo_random,               ONLY :  &
    random_state_t  ,& ! derived type: random generator state
    construct       ,& ! constructs new random generator (state) from seed
    random_number   ,& ! random number generator
    random_gauss       ! generator for Gaussian random number distribution

!==============================================================================

IMPLICIT NONE

!==============================================================================

!------------------------------------------------------------------------------
! Declarations
!------------------------------------------------------------------------------

REAL    (KIND=wp), ALLOCATABLE, PUBLIC  :: &
     pertstoph  (:,:,:)        ! stochastic multiplier of physics tendencies

REAL    (KIND=wp), ALLOCATABLE, PRIVATE :: &
     taper      (:)         ,& ! vertical tapering function  
     rlatwei_rn (:,:,:)     ,& ! weights of for horizontal interpolation of rn
     rlonwei_rn (:,:,:)     ,& !          
     randnumb   (:,:,:,:)   ,& ! random numbers on coarse grid
     randnumbint(:,:,:,:,:)    ! random numbers on model grid

REAL    (KIND=wp), PRIVATE :: &
     range_tot ,& ! max. allowed value of random numbers
     zdt          ! time step

INTEGER, PRIVATE :: &
     inisec    ,& ! initial time of model run in terms of daytime [sec]
     k50       ,& ! indexes for vertical tapering
     k100      ,& !
     ksfc      ,& !
     k870         ! 

INTEGER, ALLOCATABLE, PRIVATE :: &
     jlatbox_rn(:,:,:), & ! / Indexes of nearest rn coarse grid points for each
     ilonbox_rn(:,:,:)    ! \ each model grid point used in rn horiz. interpol.

REAL    (KIND=wp),        PARAMETER,   PRIVATE  :: &
     epsy_hr  =  1.0E-4_wp        ! safety tolerance for time checks, must be
                                  !   < 1sec in units [hrs], assuming that the
                                  !   model timestep 'dt' is a multiple of 1 sec

SAVE
 
TYPE (random_state_t),    ALLOCATABLE, PRIVATE :: yd_random_stream(:)

!------------------------------------------------------------------------------

CONTAINS


!==============================================================================
!+ Module procedure in "src_stoch_phys" to initialise creation of rn fields
!------------------------------------------------------------------------------

SUBROUTINE init_stoch_phys

!------------------------------------------------------------------------------
!
! Description:
!   Initialise random number generator, generate first set of random numbers,
!   distribute it to the other processors, and create a first field of random
!   numbers on the COSMO grid (one for each random number pattern).
!   In case of restart, it reproduces the same rn (if (imode_rn == 0) then
!   'gen_rand_numb' is called the same number of times)
!
!------------------------------------------------------------------------------

IMPLICIT NONE

REAL (KIND=wp),      PARAMETER   :: c3600 = 3600.0_wp

INTEGER    ::  &
    imode_seed,& ! = 1 : seed for first random field (for ydate <= ydate_ini)
                 ! = 2 : seed for subsequent random field
    mtstep    ,& ! loop index for timestep prior to initial timestep
    j  , n    ,& ! loop indices
    kss       ,& ! second
    irefsec   ,& ! time of day (in [sec]) of first random field
    nzday     ,& ! Julian day
    ierr         ! error control

REAL    (KIND=wp)        ::  &
    htim_ref  ,& ! (<=0): time (in forecast hours) of first random field of run
    zhr            ! hour

REAL (KIND=wp),      ALLOCATABLE :: zhelp (:)

CHARACTER (LEN=14)  :: ydate_rn ! actual date in the form   yyyymmddhh
CHARACTER (LEN=28)  :: yakdat2  ! actual date, form:   wd   dd.mm.yy  hh UTC
!------------------------------------------------------------------------------

  n1_rn = 1
  n2_rn = 2

  !   prepare initialisation of random number (rn) generator, and
  !   get time of day (in [sec]) of initial time of model run 'inisec'

  CALL prep_init_rand_numb
! ========================

! create a first set of random number fields (for each rn pattern)
! ------------------------------------------

  DO n=1,npattern_rn
    IF (imode_rn == 0) THEN
      !   construct random number generator, compute shift of rn coarse grid
      CALL init_rn_generator ( n , ydate_ini , 1 )
    ! ======================
      !   time (in timesteps) of first required random field
      hstep_rn(n)  =  0.0_wp
      !   for restart runs: generate all rn before nstart to get reproducibility
      IF (nstart > 0) THEN
        DO mtstep = 0, nstart-1
          IF (mtstep*dt/c3600-epsy_hr > hstep_rn(n) + hinc_rn(n)) THEN
            hstep_rn(n) = hstep_rn(n) + hinc_rn(n)
            IF (my_cart_id == 0)  CALL gen_rand_numb ( n1_rn(n), n )
                                ! ==================
          ENDIF
        ENDDO
      ENDIF
    ELSEIF (imode_rn >  0) THEN
      !   time of day (in [sec]) of first random field (for current rn pattern)
      irefsec   =  (inisec / NINT( hinc_rn(n)*3600.0_wp ))                  &
                  *          NINT( hinc_rn(n)*3600.0_wp )
      !   time (in forecast hours) of first random field of run: htim_ref <= 0 !
      htim_ref  =  REAL( irefsec - inisec, wp ) / c3600
      !   time (in forecast hours) of first required random field
      !   !! it is assumed that (MOD( 24 hrs , hinc_rn(n) ) == 0) for all 'n' !!
      IF (nstart == 0) THEN
        hstep_rn(n) = htim_ref
      ELSE
        hstep_rn(n) = htim_ref + hinc_rn(n)* INT( (ntstep*dt/c3600 - htim_ref) &
                                                  /hinc_rn(n) +epsy_hr )
      ENDIF
      !   date/time of first required random field
      kss = NINT( hstep_rn(n) *3600.0_wp)
      CALL get_utc_date ( kss, ydate_ini, 1.0_wp, itype_calendar, ydate_rn  &
                        , yakdat2, nzday, zhr )
    ! =================
      !   check consistency of namelist input for seeds of rn generator
      IF ((nseed_rn2(n) /= nseed_rn(n)) .AND. (ABS( htim_ref ) > epsy_hr)) THEN
        PRINT *,"CAUTION ! seed of random number generator at t=0 is different "
        PRINT *,"          from seed at the end of the previous run in DA cycle"
      ENDIF
      !   construct random number generator, compute shift of rn coarse grid
      IF (ABS( hstep_rn(n) - htim_ref ) <= epsy_hr)  imode_seed = 1
      IF (ABS( hstep_rn(n) - htim_ref ) >  epsy_hr)  imode_seed = 2
      CALL init_rn_generator ( n , ydate_rn , imode_seed )
    ! ======================
    ENDIF

    !   get first random number field on coarse grid
    IF (my_cart_id == 0)  CALL gen_rand_numb ( n1_rn(n), n )
                        ! ==================

    !   distribute rn field, and get first random number field on model grid
    CALL gen_rand_field ( n1_rn(n), n )
  ! ===================

    !   make sure that 'randnumbint(.,n2_rn,.) is defined (e.g. for ntstep=0)
    !   even if it is not computed within the first call of 'compute_stoch_phys'
    randnumbint(:,:,:,n2_rn(n),n) = randnumbint(:,:,:,n1_rn(n),n)
  ENDDO

END SUBROUTINE init_stoch_phys


!==============================================================================
!+ Module procedure in "src_stoch_phys" for computing a random field for SPPT
!------------------------------------------------------------------------------

SUBROUTINE compute_stoch_phys

!------------------------------------------------------------------------------
!
! Description:
!   Compute stochastic physics perturbations. 
!   Model grid points contained in the same box (i.e. coarse grid point) have
!   the same random number.
!   If requested:
!     - horizontal interpolation of random numbers from the coarse grid to the
!       model grid every hinc_rn hours
!     - time interpolation of random numbers every time step 
!
!------------------------------------------------------------------------------

IMPLICIT NONE

REAL (KIND=wp) , PARAMETER   :: c3600 = 3600.0_wp

INTEGER            ::  i, j, k, n, mtstep ! loop indices
INTEGER            ::  kss, nzday         ! second / Julian day
REAL    (KIND=wp)  ::  zwgt1, zwgt2       ! factors for time interpolation
REAL    (KIND=wp)  ::  ppert              ! SPPT factor from 1 rn pattern
REAL    (KIND=wp)  ::  zhr                ! hour
CHARACTER (LEN=14) :: ydate_rn ! actual date in the form   yyyymmddhh
CHARACTER (LEN=28) :: yakdat2  ! actual date, form:   wd   dd.mm.yy  hh UTC
!------------------------------------------------------------------------------

  pertstoph(:,:,:) = 0.0_wp

! compute 3-D random fields for new time levels, if required
! ----------------------------------------------------------

  DO n=1,npattern_rn

    !   compute a new random field every 'hinc_rn' hours
    IF (ntstep*dt/c3600-epsy_hr > hstep_rn(n)) THEN
      hstep_rn(n) = hstep_rn(n) + hinc_rn(n)
      IF (imode_rn >  0) THEN
        !   create a new seed    
        kss = NINT( hstep_rn(n) *3600.0_wp )
        CALL get_utc_date ( kss, ydate_ini, 1.0_wp, itype_calendar          &
                          , ydate_rn, yakdat2, nzday, zhr )
      ! =================
        CALL init_rn_generator ( n , ydate_rn , 2 )
      ! ======================
      ENDIF
      n1_rn(n) = 3 - n1_rn(n)
      n2_rn(n) = 3 - n2_rn(n)
      IF (my_cart_id == 0)  CALL gen_rand_numb ( n2_rn(n), n )
                          ! ==================
      CALL gen_rand_field ( n2_rn(n), n )
    ! ===================
    ENDIF

! compute a 3-D random field valid for current timestep for SPPT
! --------------------------------------------------------------

    !   calculation of coefficients for time interpolation if requested
    IF (ltimeint_rn) THEN
      zwgt2 = (ntstep*dt/c3600 - (hstep_rn(n)-hinc_rn(n)))/(hinc_rn(n))
      zwgt1 = 1.0_wp - zwgt2
    ELSE
      zwgt1 = 0.0_wp
      zwgt2 = 1.0_wp
    ENDIF
!   IF ((my_cart_id == 0) .AND. ((ntstep <=  1) .OR. (MOD( ntstep,36 ) == 0))) &
!     PRINT*,'SPPTt ', ntstep, my_cart_id, n, range_tot, zwgt1, zwgt2, n2_rn(n)

    DO k = 1, ke
      DO j = jstart, jend+1
        DO i = istart, iend+1
          ppert = (  zwgt1*randnumbint(i,j,k,n1_rn(n),n)                       &
                   + zwgt2*randnumbint(i,j,k,n2_rn(n),n))

          !  get sum over perturbation patterns and limit it to max/min values
          pertstoph(i,j,k) = pertstoph(i,j,k) + ppert
          IF (n == npattern_rn .AND. n > 1) THEN
            IF(pertstoph(i,j,k) > range_tot) THEN
              pertstoph(i,j,k) =  range_tot
            ELSEIF(pertstoph(i,j,k) < -range_tot) THEN
              pertstoph(i,j,k) = - range_tot
            ENDIF
          ENDIF
        END DO
      END DO
    END DO
  ENDDO

  !   perturbations are centered around 1
  pertstoph(:,:,:) = pertstoph(:,:,:) + 1.0_wp

! IF ((MOD(my_cart_id,24) == 0) .AND.((ntstep <= 1).OR.(MOD(ntstep,36) == 0))) &
!   PRINT*,'SPPTp ', ntstep, my_cart_id, MAXVAL(pertstoph(:,:,:))              &
!                                      , MINVAL(pertstoph(:,:,:))              &
!                                      , SUM   (pertstoph(:,:,:))/(ie*je*ke)

END SUBROUTINE compute_stoch_phys


!==============================================================================
!+ Module procedure in "src_stoch_phys" for generating a 3-D random field
!------------------------------------------------------------------------------

SUBROUTINE gen_rand_field ( n_rn, np )

!------------------------------------------------------------------------------
!
! Description: Generate one 3-D random field on model grid
!
!------------------------------------------------------------------------------

INTEGER, INTENT (IN)       :: n_rn  ,& ! time index for rn field
                              np       ! index of random pattern

INTEGER                    :: j, ierr
REAL(KIND=wp), ALLOCATABLE :: zhelp (:)
!------------------------------------------------------------------------------

  !   distribution of coarse grid random field to all processors
  IF (num_compute > 1) THEN
    ALLOCATE ( zhelp(ie_rn(np)) )
    DO j=1,je_rn(np)
      zhelp(:) = randnumb(1:ie_rn(np),j,n_rn,np)
      CALL distribute_values ( zhelp, ie_rn(np),0,imp_reals,icomm_compute,ierr )
    ! ~~~~~~~~~~~~~~~~~~~~~~
      IF (my_cart_id /= 0)  randnumb(1:ie_rn(np),j,n_rn,np) = zhelp(:)
    ENDDO
    DEALLOCATE ( zhelp )
  ENDIF

  !   horizontal interpolation to model grid and vertical tapering
  CALL int_rand_numb ( istart, iend+1, jstart, jend+1, 1, ke                   &
                     , n_rn, np, randnumbint(:,:,:,n_rn,np) )
! ==================

END SUBROUTINE gen_rand_field


!==============================================================================
!+ Module procedure in "src_stoch_phys" for preparing random field production
!------------------------------------------------------------------------------

SUBROUTINE prep_init_rand_numb

!------------------------------------------------------------------------------
!
! Description: Allocate required fields, compute vertical tapering function
!              and rn coarse grid and its parameters
!
!------------------------------------------------------------------------------

INTEGER       ::  khh, kmi  ,& ! hour, minute
                  k  , n       ! loop indices
REAL(KIND=wp) ::  p50, p100, p870, psfc, zpno, zpnu 
!------------------------------------------------------------------------------

  IF ( l2tls ) THEN
    zdt   = dt
  ELSE
    zdt   = dt2
  ENDIF

! memory allocation
! -----------------

  ALLOCATE (yd_random_stream(npattern_rn))
  ALLOCATE (rlatwei_rn (1+1:ie_tot-1, 1+1:je_tot-1, npattern_rn),              &
            rlonwei_rn (1+1:ie_tot-1, 1+1:je_tot-1, npattern_rn))
  ALLOCATE (jlatbox_rn (1+1:ie_tot-1, 1+1:je_tot-1, npattern_rn),              &
            ilonbox_rn (1+1:ie_tot-1, 1+1:je_tot-1, npattern_rn))
  ALLOCATE (taper(1:ke_tot))
  taper=1.0_wp

! ALLOCATE (randnumbint (1:ie, 1:je, ke, 2, npattern_rn))
  ALLOCATE (randnumbint (istart:iend+1, jstart:jend+1, ke, 2, npattern_rn))
  randnumbint= 0.0_wp

  ALLOCATE (pertstoph (1:ie, 1:je, ke))
  pertstoph  = 1.0_wp

! vertical tapering function definition
! -------------------------------------

  IF(itype_vtaper_rn > 0.AND.itype_vtaper_rn <= 2)THEN
    zpno=0.0_wp
    DO k = 1,ke_tot
      zpno = refatm%p0sl*vcoord%sigm_coord(k  )
      zpnu = refatm%p0sl*vcoord%sigm_coord(k+1)
      IF ( (zpno <= 870.0E2_wp) .AND. (870.0E2_wp < zpnu) ) THEN
        k870 = k
        p870 = zpno
      ENDIF
      IF ( (zpno <= 100.0E2_wp) .AND. (100.0E2_wp < zpnu) ) THEN
        k100 = k
        p100 = zpno
      ENDIF
      IF ( (zpno <= 50.0E2_wp) .AND. (50.0E2_wp < zpnu) )  THEN
        k50 = k
        p50 = zpno
      ENDIF
    ENDDO
    ksfc = ke_tot+1
    psfc = refatm%p0sl*vcoord%sigm_coord(ksfc)
    DO k= 1,ke_tot
      IF(k <= k50)taper(k)=0.0_wp
      IF(k > k50 .AND. k < k100)taper(k)= &
         (refatm%p0sl*0.5_wp*(vcoord%sigm_coord(k)+ &
          vcoord%sigm_coord(k+1))-p50)/(p100-p50)
      IF(k >= k100.AND. k <= k870)taper(k)=1.0_wp
      IF (itype_vtaper_rn == 2)THEN
        IF(k > k870)taper(k)=1.0_wp
      ELSE
        IF(k > k870 .AND. k < ksfc)taper(k)= &
          (refatm%p0sl*0.5_wp*(vcoord%sigm_coord(k)+ &
           vcoord%sigm_coord(k+1))-psfc)/(p870-psfc)
!        IF(k >= ksfc) taper(k)=0.0_wp
      ENDIF
    ENDDO
  ELSEIF(itype_vtaper_rn == 3)THEN
    k50=0
    ksfc=ke_tot+1
    DO k= 1,ke_tot
      taper(k)=rvtaper_rn(k)
    ENDDO
  ELSE
    k50=0
    ksfc=ke_tot+1
  ENDIF

  IF (my_cart_id==0 .AND. itype_vtaper_rn > 0)                                &
    PRINT*,' SPPT vertical tapering at model levels:', taper

! various further preparations
! ----------------------------

  range_tot=0.0_wp
  DO n=1,npattern_rn
!
    ! number of grid pts. of coarse grid of random numbers:
    !   Coarse grid has to cover COSMO grid with an extra half grid spacing
    !   due to the wind staggering and without two grid spacing due to the
    !   external boundary frame that have values defined by lateral BC
!   dim_rn = 2 + INT(((idim_tot-1)*ddim + 0.5*ddim - 2.0 * ddim)/ddim_rn(n))
    ie_rn(n) = 2+INT(((REAL(ie_tot,wp)-2.5_wp)*dlon)/dlon_rn(n))+1
    je_rn(n) = 2+INT(((REAL(je_tot,wp)-2.5_wp)*dlat)/dlat_rn(n))+1

    !   max. allowed value of random numbers
    range_tot = range_tot + (range_rn(n)**2)
  ENDDO
  range_tot = SQRT(range_tot)
  IF (my_cart_id == 0)  PRINT '(" Total range of random numbers for SPPT is " &
                              &,F8.4)' , range_tot

  ALLOCATE (randnumb(MAXVAL(ie_rn),MAXVAL(je_rn),2,npattern_rn))

  !   time of day (in [sec]) of initial time of model run
  IF (imode_rn >  0) THEN
    READ ( ydate_ini(9:10), '(I2.2)' ) khh 
    READ ( ydate_ini(11:12),'(I2.2)' ) kmi 
    inisec = khh*3600 + kmi*60
  ENDIF

END SUBROUTINE prep_init_rand_numb


!==============================================================================
!+ Module procedure in "src_stoch_phys" for initialising random numb. generator
!------------------------------------------------------------------------------

SUBROUTINE init_rn_generator ( n , ydate , imode_seed )

!----------------------------------------------------------------------------
!
! Description: Set the seed of a random number stream,
!              and compute parameters of rn coarse grid
!
!----------------------------------------------------------------------------

IMPLICIT NONE

CHARACTER  (LEN=14), INTENT(IN) :: ydate ! reference date of pattern 'n'
INTEGER, INTENT(IN)             :: n     ! index of random pattern
INTEGER, INTENT(IN)             :: imode_seed
   ! mode of seed: = 1 : seed for first random field (for ydate <= ydate_ini)
   !               = 2 : seed for subsequent random field

INTEGER            :: i, j, ier, kconseed, ijshift(2)
INTEGER            :: iseed       ! seed for random number generator
!REAL(KIND=wp)     :: zshift(2)
REAL (KIND=dp)     :: zshift(2)
!----------------------------------------------------------------------------

! initialise random number generator
! ----------------------------------

  IF (my_cart_id == 0) THEN
    IF (imode_seed == 1)  iseed = nseed_rn (n)
    IF (imode_seed == 2)  iseed = nseed_rn2(n)
    kconseed = iepsmem + (n-1)*MAX( iepstot,2000 )

    CALL set_seed_rand_numb ( ydate, kconseed , iseed )
  ! =======================
    CALL construct ( yd_random_stream(n), seed=iseed )
  ! ==============
    PRINT*,'SPPT random pattern',n,': seed:',iseed,' for date ',ydate, kconseed
  ENDIF

! (definition of correspondence between rn coarse and model grid:)
! compute random shift of SW corner of the coarse grid
! ----------------------------------------------------

  IF (my_cart_id == 0) THEN
    !   draw random numbers in [0,1] (with uniform distribution)
    CALL random_number (zshift(1:2), yd_random_stream(n))
  ! ==================
    WHERE (zshift == 1.0_wp)
      zshift = 0.0_wp
    END WHERE
    ijshift(1) = INT( zshift(1)*dlat_rn(n) /dlat )
    ijshift(2) = INT( zshift(2)*dlon_rn(n) /dlon )
  ENDIF
  IF (num_compute > 1)                                                         &
    CALL distribute_values ( ijshift, 2, 0, imp_integers, icomm_compute, ier )
  ! ~~~~~~~~~~~~~~~~~~~~~~

  DO j=2,je_tot-1
    jlatbox_rn(:,j,n) = INT( REAL(j+ijshift(1)-2,wp)*dlat/dlat_rn(n) ) + 1
    rlatwei_rn(:,j,n) = MOD( REAL(j+ijshift(1)-2,wp)*dlat, dlat_rn(n) )    &
                       /dlat_rn(n)
  ENDDO
  DO i=2,ie_tot-1
    ilonbox_rn(i,:,n) = INT( REAL(i+ijshift(2)-2,wp)*dlon/dlon_rn(n) ) + 1
    rlonwei_rn(i,:,n) = MOD( REAL(i+ijshift(2)-2,wp)*dlon, dlon_rn(n) )    &
                       /dlon_rn(n)
  ENDDO
 
  IF (my_cart_id == 0) THEN
    PRINT '(" SPPT random pattern",I2," coarse grid: size",I3," x",I3          &
          &,"; SW corner shifted by",I3," /",I3," g.p.")'                      &
           , n, MAXVAL(ilonbox_rn(:,:,n)), MAXVAL(jlatbox_rn(:,:,n)), ijshift
  ENDIF

END SUBROUTINE init_rn_generator


!==============================================================================
!+ Module procedure in "src_stoch_phys" for drawing random numbers 
!------------------------------------------------------------------------------

SUBROUTINE gen_rand_numb ( n_rn, np )

!----------------------------------------------------------------------------
!
! Description: Generate a set of random numbers 
!
!----------------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT (IN) :: n_rn  ,& ! time index for rn field
                        np       ! index of random pattern

INTEGER                               :: j, i
!REAL   (KIND=wp), ALLOCATABLE        :: zglptsu(:,:)
REAL    (KIND=dp), ALLOCATABLE        :: zglptsu(:,:)
CHARACTER (LEN=15)                    :: yformat
!----------------------------------------------------------------------------
 
  ALLOCATE(zglptsu(ie_rn(np),je_rn(np)))
  zglptsu=0.0_wp

  !   a random number is defined for each coarse grid point (box), so that
  !   model grid points inside the same box have the same random number
  DO j=1,je_rn(np)
    IF (lgauss_rn) THEN

      !   N(0,1) random numbers (zglptsu) .
      CALL random_gauss (zglptsu(:,j), yd_random_stream(np))
    ! =================

      !   transformation to N(0,stdv_rn) random numbers
      zglptsu(:,j) = zglptsu(:,j) *stdv_rn(np)

      !   limitation of maximum/minimum values (+/- range_rn) for each pattern
      WHERE (zglptsu(:,j) > range_rn(np))
        zglptsu(:,j) = range_rn(np)      
      ELSEWHERE (zglptsu(:,j) < -range_rn(np))
        zglptsu(:,j) = -range_rn(np)      
      END WHERE
    ELSE
      !   random numbers with uniform distribution in [0,1]
      CALL random_number (zglptsu(:,j), yd_random_stream(np))
    ! ==================
      zglptsu(:,j) = range_rn(np)-2.0_wp*range_rn(np)*zglptsu(:,j)
    ENDIF
  ENDDO

  randnumb(1:ie_rn(np),1:je_rn(np),n_rn,np) = zglptsu(:,:)

  PRINT '(" SPPT random pattern",I2,": min, max, mean:",3F8.4)'                &
         , np, MINVAL(randnumb(1:ie_rn(np),1:je_rn(np),n_rn,np))               &
             , MAXVAL(randnumb(1:ie_rn(np),1:je_rn(np),n_rn,np))               &
             , SUM(zglptsu(:,:))/(ie_rn(np)*je_rn(np))

  DEALLOCATE(zglptsu)

! WRITE( yformat,'("(A4,I4,F6.2,2I4,",I2,"F5.2)")' ) ie_rn(np)
! DO j =1 ,je_rn(np)
!   PRINT (yformat), 'SPPTrn ', np, hstep_rn(np), n_rn, j                      &
!                             , (randnumb(i,j,n_rn,np),i=1,ie_rn(np))
! ENDDO

END SUBROUTINE gen_rand_numb

!==============================================================================
!+ Module procedure in "src_stoch_phys" for interpolating random number fields
!------------------------------------------------------------------------------

SUBROUTINE int_rand_numb ( ki1sc, ki1ec, ki2sc, ki2ec, ki3s, ki3e              &
                         , kn_rn, np , prandnumbint )

!------------------------------------------------------------------------------
!
! Description: Horizontal interpolation and vertical tapering of random numbers 
!
!
!------------------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN) :: &
              ki1sc  , &! first  dimension of calculation, start index
              ki1ec  , &! first  dimension of calculation, end   index
              ki2sc  , &! second dimension of calculation, start index
              ki2ec  , &! second dimension of calculation, end   index
              ki3s   , &! third dimension of calculation,  start index
              ki3e      ! third dimension of calculation,  end   index
INTEGER,        INTENT(IN)   :: kn_rn,np       
REAL (KIND=wp), INTENT(OUT)  :: &
              prandnumbint(ki1sc:ki1ec,ki2sc:ki2ec,ki3s:ki3e)
!REAL   (KIND=wp),        INTENT(OUT)  :: prandnumbint(1:ie,1:je,1:ke)

REAL    (KIND=wp)        :: zrand_help(ki1sc:ki1ec,ki2sc:ki2ec)
INTEGER                  :: i,j,k
INTEGER                  :: ii,jj
!------------------------------------------------------------------------------

  zrand_help(ki1sc:ki1ec,ki2sc:ki2ec)=0.0_wp

! Horizontal interpolation of the random numbers (smoothed values)
  IF (lhorint_rn) THEN
    DO jj=ki2sc,ki2ec
      DO ii=ki1sc,ki1ec
        i=i_global(ii)
        j=j_global(jj)
        zrand_help(ii,jj) =                                                    &
              randnumb(ilonbox_rn(i,j,np)  ,jlatbox_rn(i,j,np)  ,kn_rn,np)     &
              *(1.0_wp-rlonwei_rn(i,j,np))*(1.0_wp-rlatwei_rn(i,j,np))         &
            + randnumb(ilonbox_rn(i,j,np)+1,jlatbox_rn(i,j,np)  ,kn_rn,np)     &
              *   rlonwei_rn(i,j,np)*(1.0_wp-rlatwei_rn(i,j,np))               &
            + randnumb(ilonbox_rn(i,j,np)  ,jlatbox_rn(i,j,np)+1,kn_rn,np)     &
              *(1.0_wp-rlonwei_rn(i,j,np))*  rlatwei_rn(i,j,np)                &
            + randnumb(ilonbox_rn(i,j,np)+1,jlatbox_rn(i,j,np)+1,kn_rn,np)     &
              *   rlonwei_rn(i,j,np)*   rlatwei_rn(i,j,np)
      ENDDO
    ENDDO
  ELSE
! No horizontal interpolation of the random numbers (constant values in boxes)
    DO jj=ki2sc,ki2ec
      DO ii=ki1sc,ki1ec
        i=i_global(ii)
        j=j_global(jj)
        zrand_help(ii,jj)= &
          randnumb(ilonbox_rn(i,j,np),jlatbox_rn(i,j,np),kn_rn,np)
      ENDDO
    ENDDO
  ENDIF

! Vertical tapering of the random numbers as prescribed
  IF(itype_vtaper_rn > 0)THEN
    DO j=ki2sc,ki2ec
      DO i=ki1sc,ki1ec
        DO k=ki3s,ki3e
          prandnumbint(i,j,k) = zrand_help(i,j)*taper(k)
        ENDDO
      ENDDO
    ENDDO
  ELSE
! No vertical tapering of the random numbers 
    DO j=ki2sc,ki2ec
      DO i=ki1sc,ki1ec
        DO k=ki3s,ki3e
          prandnumbint(i,j,k)= zrand_help(i,j)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

END SUBROUTINE int_rand_numb

!==============================================================================
!+ Module procedure in "src_stoch_phys" for set the seed of random numbers
!==============================================================================

SUBROUTINE set_seed_rand_numb ( ydate, kconseed , iseed )

!------------------------------------------------------------------------------
!
! Description:
!      From IFS library of ECMWF (setran.F)
!      *set_seed_rand_numb* - Sets the seed for a random number stream
!      as a function of initial date and member number  
!      
! Method:
!      A seed 'iseed' is set as a function of 'ydate' (which can be e.g. the
!      initial time of the model run), of integer 'kconseed' (which can depend
!      e.g. on the ensemble member number), and of the input value of 'iseed'
!      itself (which may depend e.g. on a namelist value).
!      For dates chosen at random, the seeds are approximately uniformly
!      distributed between 1 and HUGE(0). A highly nonlinear function is used
!      to reduce the possibility of correlations between random sequences
!      generated for different initial dates.
!
!------------------------------------------------------------------------------
 
IMPLICIT NONE
 
CHARACTER(LEN=14), INTENT(IN)    :: ydate
INTEGER,           INTENT(IN)    :: kconseed
INTEGER,           INTENT(INOUT) ::  &
  iseed       ! seed constructed 'kconseed', 'ydate', and input value of 'iseed'

INTEGER, PARAMETER     ::  &
  ibtshf = 11 ! number of bits to shift 'kconseed'
INTEGER   :: idigits, iradix, is, jdigit, iscale, &
             inseed, keseed, nindat, nsssss
INTEGER   :: ndd, nmm, nccaa, naa, namd, ncth, nzzaa,   &
             nzzmm, ncent, nyearc 
INTEGER   :: kgrdat, ksec, kaaaa, kmm, kdd, kss
INTEGER   :: kcent, kyearc, kmonth, kday, khh, kmi
REAL (KIND=dp)     :: zirr1, zs, zt, ztim
REAL (KIND=dp)     :: rjudat
!INTEGER,PARAMETER :: NINDAT=19880826 , NSSSSS=0
 
ndd(kgrdat)   = MOD(kgrdat,100)
nmm(kgrdat)   = MOD((kgrdat-ndd(kgrdat))/100,100)
nccaa(kgrdat) = kgrdat/10000
nzzaa(kaaaa,kmm) = kaaaa-( (1-SIGN(1,kmm-3))/2 )
nzzmm(kmm) = kmm+6*(1-SIGN(1,kmm-3))
rjudat(kaaaa,kmm,kdd) = 1720994.5_dp + REAL(2-nzzaa(kaaaa,kmm)/100  &
                        + (nzzaa(kaaaa,kmm)/100)/4 &
                        + INT(365.25_dp * REAL(nzzaa(kaaaa,kmm),dp))&
                        + INT(30.601_dp * REAL(nzzmm(kmm)+1,dp))    &
                        + kdd,dp)

! End of header
!------------------------------------------------------------------------------
  inseed = iseed
  keseed = kconseed

  !   produce a number from 'iseed' and 'kconseed'
  keseed  = IOR( ISHFT( IBITS( iseed , 0 , 30-ibtshf ) , ibtshf ) , keseed )

  READ ( ydate(1:4) , '(I4.4)' ) kaaaa
  READ ( ydate(5:6) , '(I2.2)' ) kmm
  READ ( ydate(7:8) , '(I2.2)' ) kdd
  READ ( ydate(9:10), '(I2.2)' ) khh 
  READ ( ydate(11:12),'(I2.2)' ) kmi 
  nindat = kaaaa*10000+kmm*100+kdd
  nsssss = khh*3600 + kmi*60
  iradix  = RADIX(ztim)
  idigits = DIGITS(ztim)
  zirr1 = 0.5_dp*(SQRT(5.0_dp)-1.0_dp)

!--- generate a unique number from the date and the input keseed
 
  ztim =   rjudat(nccaa(nindat),nmm(nindat),ndd(nindat))              &
         - 1720994.5_dp +  REAL(nsssss,dp)/86400.0_dp             &
         - 2581470.3_dp * keseed  
 
!--- multiply by an irrational number to randomize the bits and scale
!--- to between 0 and 1.
 
  ztim = FRACTION(zirr1*ABS(ztim))
 
!--- reverse the bits
 
  zs = 0.0_dp
  zt = ztim
  DO jdigit = 1, idigits
    zt = zt*iradix
    is = int(zt)
    zt = zt-is
    zs = (zs+is)/iradix
  ENDDO

!--- Scale to an odd number between 0 and HUGE-100000000
!--- (Allow some headroom in order to use set_seed_rand_numb to set an initial 
!---  seed and then generate new seeds by incrementing.)
 
  iscale = (HUGE(iseed)-100000000)/2
  iseed = 1 + 2*INT( iscale*zs )
 
END SUBROUTINE set_seed_rand_numb


!==============================================================================
!+ Module procedure in "src_stoch_phys" for checking neg./sup.sat. humidity 
!------------------------------------------------------------------------------

SUBROUTINE apply_tqx_tend_adj ( iitype_qxpert_rn, iitype_qxlim_rn, zpf         &
                              , zt, zqv, zqc, zqi, zqr, zqs, zpertu , zttens   &
                              , zqvtens, zqctens, zqitens, zqrtens, zqstens    &
                              , lnopertu , zqg , zqgtens )
!------------------------------------------------------------------------------
!
!    itype_qxpert_rn define which hum variables tend. are perturbed
!                 0:    Only qv tendencies are perturbed
!                 1:    qv,qc,qi tendencies are perturbed
!                 2:    qv,qc,qi,qr,qs,qg tendencies are perturbed
!
!
!    itype_qxlim_rn  type of reduction/removal of the perturbation 
!                    in case of negative (qv, qc, qi) or 
!                    supersaturated (qv) values
!                 0:    No limitation of perturbed tendencies
!                 1:    If new qv values are negative or super-sat. -> T and qv
!                       tendencies are not perturbed
!                       If new qx (qc,qi,qr,qs,qg values) are negative -> qx
!                       tendencies are not perturbed
!                 2:    (currently discarded:)
!                       Compute the perturbed tendencies of T and qv, such that
!                       the new values do not exceed the limits for qv, and the
!                       limitation procedure does not introduce any bias by
!                       definition.  -  The lower limit for qv is set to zero,
!                       the upper limit to the specific water vapour content
!                       at saturation over water.
!------------------------------------------------------------------------------

IMPLICIT NONE

INTEGER, INTENT(IN) :: &
  iitype_qxlim_rn ,& ! defines which hydrometeor tend. are perturbed (see above)
  iitype_qxpert_rn   ! defines limiter to tendencies (see above)

REAL   (KIND=wp),    INTENT(IN) :: &
  zpertu,  & ! stochastic multiplier of physics tendencies
  zt,      & ! model field: temperature
  zqv,     & ! model field: specific water vapour content
  zqc,     & ! model field: specific cloud water  content
  zqi,     & ! model field: specific cloud  ice   content
  zqr,     & ! model field: specific rain  water  content
  zqs,     & ! model field: specific snow  water  content
  zpf        ! model field: pressure

REAL (KIND=wp),      INTENT(INOUT) :: &
  zttens,  & ! tendency of temperature
  zqvtens, & ! tendency of specific water vapour content
  zqctens, & ! tendency of specific cloud water  content
  zqitens, & ! tendency of specific cloud  ice   content
  zqrtens, & ! tendency of specific rain  water  content
  zqstens    ! tendency of specific snow  water  content

LOGICAL,             INTENT(OUT)   :: &
  lnopertu   ! .true.: no tendency perturbation applied (at current grid pt.)

REAL (KIND=wp),      INTENT(IN)   , OPTIONAL ::  &
  zqg        ! model field: specific graupel water content

REAL (KIND=wp),      INTENT(INOUT), OPTIONAL ::  &
  zqgtens    ! tendency of specific graupel water content

REAL (KIND=wp)     :: zqvsp, ztten, zqvten, zqcten, zqiten, zqvn, ztn, zqcn,   &
                      zqrten, zqsten, zqgten, zqrn, zqsn, zqgn, zqin, zpv

!REAL (KIND=wp),      PARAMETER :: c0 = 0.0_wp, c1 = 1.0_wp
!REAL (KIND=wp) ::  &
! zlim_low ,& ! lower limit to new value of qv
! zlim_upp ,& ! upper limit to new value of qv
! zapertu  ,& ! test perturbation factor applied to tendency
! zfpert      ! perturbation factor applied to tendency, replaces 'zpertu'
!------------------------------------------------------------------------------

  lnopertu = .false.

  IF(iitype_qxlim_rn == 0)THEN
    zttens = zttens * zpertu
    zqvtens = zqvtens * zpertu
    IF(iitype_qxpert_rn >= 1)THEN
      zqctens = zqctens * zpertu
      zqitens = zqitens * zpertu
    ENDIF
    IF(iitype_qxpert_rn == 2)THEN
      zqrtens = zqrtens * zpertu
      zqstens = zqstens * zpertu
      IF (PRESENT(zqg)) THEN
        zqgtens = zqgtens * zpertu
      ENDIF
    ENDIF

! Original IFS humidity check
  ELSEIF(iitype_qxlim_rn == 1)THEN
    ztten = zttens * zpertu
    zqvten = zqvtens * zpertu
    zqvn=zqv+zdt*zqvten
    ztn=zt+zdt*ztten
    zpv   = fpvsw( ztn, b1, b2w, b3, b4w )
    zqvsp = fpv2q( zpv, zpf, rdv )
!   zqvsp = fpv2q( fpvsw( ztn, b1, b2w, b3, b4w ), zpf, rdv )
!   IF (zqvn < 0. .OR. (zqvn > zqvsp .AND.ztn > 248.15)) THEN
    IF (zqvn < 0. .OR. zqvn > zqvsp ) THEN
        lnopertu= .true.
    ELSE
        zttens = ztten
        zqvtens = zqvten
    ENDIF

    IF(iitype_qxpert_rn >= 1)THEN
      zqcten = zqctens * zpertu
      zqcn=zqc+zdt*zqcten
      IF (zqcn >= 0.0_wp) THEN
        zqctens = zqcten
      ENDIF

      zqiten = zqitens * zpertu
      zqin=zqi+zdt*zqiten
      IF (zqin >= 0.) THEN
        zqitens = zqiten
      ENDIF
    ENDIF

    IF(iitype_qxpert_rn == 2)THEN
      zqrten = zqrtens * zpertu
      zqrn=zqr+zdt*zqrten
      IF (zqrn >= 0.0_wp) THEN
        zqrtens = zqrten
      ENDIF

      zqsten = zqstens * zpertu
      zqsn=zqs+zdt*zqsten
      IF (zqsn >= 0.0_wp) THEN
        zqstens = zqsten
      ENDIF

      IF (PRESENT(zqg)) THEN
        zqgten = zqgtens * zpertu
        zqgn=zqg+zdt*zqgten
        IF (zqgn >= 0.0_wp) THEN
          zqgtens = zqgten
        ENDIF
      ENDIF
    ENDIF

! (by Christoph Schraff: This code should avoid a dry bias of the model,
!                        but still appears to have some problems)
! ELSEIF(iitype_qxlim_rn == 2)THEN

!   zfpert   = zpertu

    !  limit the tendency perturbation such that if the modulus of the limited
    !  tendency perturbation is added to the original tendency, then the limits
    !  are not exceeded  (it is assumed that 'zpertu > 0')
      !   test perturbation factor = 1.0 + modulus of perturbation random number
!     zapertu = 1.0_wp + ABS( zpertu - 1.0_wp )
      !   lower limit on qv : 0.0
!     zlim_low = c0
      !   upper limit on qv : qv_sat over water (not only imposed for T > -25 C)
      !   (qv_sat is related to the minimum of new T with or w/o perturbation;
      !   unless the safe approximation of taking the minimum was used here,
      !   an iterative procedure would be needed)
!     zlim_upp = fpv2q( fpvsw( zt + zdt *MIN( zttens, zttens *zapertu )        &
!                           , b1, b2w, b3, b4w ), zpf, rdv )
      !   new qv-value without perturbation
!     zqvn = zqv + zdt *zqvtens
      !   if the unperturbed new qv-value exceeds the limits then do not perturb
!     IF ((zqvn <= zlim_low) .OR. (zqvn >= zlim_upp)) THEN
!       zfpert   = c1
!       lnopertu = .true.
      !   check if adding the test tendency lets 'qv_n' exceed the upper limit
!     ELSEIF ((zqv + zdt*zqvtens* zapertu > zlim_upp) .AND. (zqvtens > c0)) THEN
        !   zlim_upp - zqvn  : modulus of limited (adjusted) tendency perturbat.
!       zfpert  =  c1  +  SIGN( (zlim_upp - zqvn) / zqvtens , zpertu - c1 )
      !   check if adding the test tendency lets 'qv_n' go below the lower limit
!     ELSEIF ((zqv + zdt*zqvtens* zapertu < zlim_low) .AND. (zqvtens < c0)) THEN
!       zfpert  =  c1  +  SIGN( (zqvn - zlim_low) / zqvtens , zpertu - c1 )
!     ENDIF
!   zttens  = zttens  * zfpert
!   zqvtens = zqvtens * zfpert

!   IF (iitype_qxpert_rn >= 1) THEN
!     IF (     (zqc + zdt*zqctens         >= zlim_low)                         &
!         .AND.(zqc + zdt*zqctens* zapertu < zlim_low) .AND.(zqctens < c0)) THEN
!       zqctens = zqctens + SIGN( zqc + zdt *zqctens - zlim_low , zpertu - c1 )
!     ELSEIF (zqc + zdt*zqctens* zapertu >= zlim_low) THEN
!       zqctens = zqctens * zpertu
!     ENDIF
!     IF (     (zqi + zdt*zqitens         >= zlim_low)                         &
!         .AND.(zqi + zdt*zqitens* zapertu < zlim_low) .AND.(zqitens < c0)) THEN
!       zqitens = zqitens + SIGN( zqi + zdt *zqitens - zlim_low , zpertu - c1 )
!     ELSEIF (zqi + zdt*zqitens* zapertu >= zlim_low) THEN
!       zqitens = zqitens * zpertu
!     ENDIF
!   ENDIF
    
!   IF (iitype_qxpert_rn == 2) THEN
!     IF (     (zqr + zdt*zqrtens         >= zlim_low)                         &
!         .AND.(zqr + zdt*zqrtens* zapertu < zlim_low) .AND.(zqrtens < c0)) THEN
!       zqrtens = zqrtens + SIGN( zqr + zdt *zqrtens - zlim_low , zpertu - c1 )
!     ELSEIF (zqr + zdt*zqrtens* zapertu >= zlim_low) THEN
!       zqrtens = zqrtens * zpertu
!     ENDIF
!     IF (     (zqs + zdt*zqstens         >= zlim_low)                         &
!         .AND.(zqs + zdt*zqstens* zapertu < zlim_low) .AND.(zqstens < c0)) THEN
!       zqstens = zqstens + SIGN( zqs + zdt *zqstens - zlim_low , zpertu - c1 )
!     ELSEIF (zqs + zdt*zqstens* zapertu >= zlim_low) THEN
!       zqstens = zqstens * zpertu
!     ENDIF
!     IF (PRESENT(zqg)) THEN
!       IF (   (zqg + zdt*zqgtens         >= zlim_low)                         &
!         .AND.(zqg + zdt*zqgtens* zapertu < zlim_low) .AND.(zqgtens < c0)) THEN
!         zqgtens = zqgtens + SIGN( zqg + zdt *zqgtens - zlim_low, zpertu - c1 )
!       ELSEIF (zqg + zdt*zqgtens* zapertu >= zlim_low) THEN
!         zqgtens = zqgtens * zpertu
!       ENDIF
!     ENDIF
!   ENDIF
  ENDIF

END SUBROUTINE apply_tqx_tend_adj

!-------------------------------------------------------------------------------

ELEMENTAL REAL FUNCTION fpvsw  ( zt, b1, b2w, b3, b4w )
  !----------------------------------------------------
  REAL    (KIND=wp), INTENT (IN)  ::  zt, b1, b2w, b3, b4w
  !---------------------------------------------------------------------------
  ! Magnus formula for water:  input  'zt'   : temperature
  !                            output 'fpvsw': saturation water vapour pressure
  !        b1, b2w, b3, b4w :  constants of Magnus formula for water
  !---------------------------------------------------------------------------
  !
  fpvsw  =  b1 * EXP( b2w *(zt-b3) /(zt-b4w) )
  !
END FUNCTION fpvsw

!-------------------------------------------------------------------------------

ELEMENTAL REAL FUNCTION fpv2q  ( zpv, zp, rdv )
  !--------------------------------------------
  REAL    (KIND=wp), INTENT (IN)  ::  zpv, zp, rdv
  !---------------------------------------------------------------------------
  ! specific humidity from water vapour pressure 'zpv' and air pressure 'zp'
  !   (rdv = r_d / r_v
  !        = gas constant for dry air / gas constant for water vapour )
  !---------------------------------------------------------------------------
  !
  fpv2q  =  rdv * zpv / (zp - (1.0_wp-rdv)*zpv)
  !
END FUNCTION fpv2q

!------------------------------------------------------------------------------

END MODULE src_stoch_physics
