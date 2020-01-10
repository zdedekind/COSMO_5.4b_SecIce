!+ Module containing a routine for computing the horizontal diffusion
!------------------------------------------------------------------------------

MODULE hori_diffusion

!------------------------------------------------------------------------------
!
! Description:
!   This module provides a routine to calculate the horizontal diffusion 
!   for the atmospheric variables.
!
!   Routine(s) currently contained:
!
!       comp_hori_diff
!
! Current Code Owner: DWD, Michael Baldauf
!  phone:  +49  69  8062 2733
!  fax:    +49  69  8062 3721
!  email:  Michael.Baldauf@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 3.2        2003/02/07 Ulrich Schaettler
!  Initial release
! 3.5        2003/09/02 Ulrich Schaettler
!  Adaptation of interface for exchg_boundaries
! 3.7        2004/02/18 Jochen Foerstner / Ulrich Schaettler
!  Adaptations for horizontal diffusion of moisture quantities with orographic
!  limiter (for new Runge-Kutta scheme)
!  Renamed cphi (crlat)
! 3.16       2005/07/22 Ulrich Schaettler
!  Changes for the use in the Runge-Kutta scheme (converted SR to MODULE)
! 3.18       2006/03/03 Jochen Foerstner
!  Avoid unnecessary computations in case of hd_corr_q=0
! 3.21       2006/12/04 Jochen Foerstner
!  Introduction of 3D mask array hd_mask_dcoeff
! V4_4         2008/07/16 Ulrich Schaettler
!  Adapted interface of get_timings
! V4_8         2009/02/16 Oliver Fuhrer
!  Initialize fields for limiters, if new NL switch linit_fields is set
! V4_9         2009/07/16 Ulrich Schaettler, Heike Vogel, Christian Bollmann
!  Introduce first version of COSMO_ART
!  Adaptations for 3D versions of routines in numeric_utilities
!  Implemented different mask fields for t, q, and u-fields
! V4_12        2010/05/11 Oli Fuhrer
!  Introduced option to treat horizontal diffusion differently in the interior
!  and the boundary zone of the domain; pressure is now also treated differently 
!  than temperature; Option 3 from itype_hdiff has been eliminated
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_17        2011/02/24 Ulrich Blahak
!  Adapted interface of exchg_boundaries; corrected kzdims(1:20) -> kzdims(1:24);
!  eliminated my_peri_neigh;
!  changed ldatatypes to .FALSE. for all communication for
!  COSMO-ART and POLLEN
! V4_18        2011/05/26 Ulrich Schaettler
!  Changed the code owner
! V4_21        2011/12/06 Michael Baldauf
!  new routines for Smagorinsky diffusion of u and v, which need also
!  additional boundary exchanges
! V4_23        2012/05/10 Ulrich Schaettler
!  Removed src_2timelevel and related stuff
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer, Hans-Juergen Panitz
!  Replaced qx-variables by using them from the tracer module
!  Introduced nexch_tag for MPI boundary exchange tag to replace ntstep (HJP)
! V4_26        2012/12/06 Anne Roches
!  Replacement of hd_corr_q_XXX by hd_corr_trcr_XXX in order to be consistent
!  also with the naming of other switches (e.g. ltrcr_trilin, lef_adv_trcr_notpd)
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  MESSy interface introduced: horizontal diffusion of tracers
! V4_28        2013/07/12 KIT, Ulrich Schaettler
!  Changes to adapt COSMO-ART to new tracer module: all dependencies to
!  COSMOART and POLLEN deleted, because this is now handled by the tracer module
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
! V5_1         2014-11-28 Michael Baldauf, Oliver Fuhrer
!  New subr. 'targeted_diffusion_cold_pools': performs 2nd order diffusion on
!  temperature where deviation from T0 is too big
!  Replaced ireals by wp (working precision) (OF)
!  Fixed allocation of [i,j,k]_liste: (ie*je*(dk+1)) instead of (ie*je*dk) (UB)
! V5_2         2015-05-21 Ulrich Schaettler
!  Reorganized numbering of MPI datatypes ID
! V5_3         2015-10-09 Andrea Arteaga, Oliver Fuhrer
!  Bug fix in calculation of ztha (for targeted diffusion) (A. Arteaga)
!  Bug fix in indexing of Smagorinsky-diffusion (A. Arteaga)
!  Excessive 'warm pools' now are removed by the targeted diffusion for 
!    cold pools (OF)
!  Implemented serialization comments (OF)
! V5_4         2016-03-10 Oliver Fuhrer, Pascal Spoerri, Xavier Lapillonne
!  Allow the targeted cold pool diffusion threshold to be controlled by the
!  thresh_cold_pool variable.
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

USE data_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &


! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------
    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction
    ke1,          & ! ke+1

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
    eddlon,       & ! 1 / dlon
    eddlat,       & ! 1 / dlat

! 5. variables for the time discretization and related variables
! --------------------------------------------------------------
    dt              ! long time-step

! end of data_modelconfig

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------
    rho0       ,    & ! reference density at the full model levels    (kg/m3)
    t0         ,    & ! reference temperature                         ( K   )
    p0         ,    & ! reference pressure at main levels             ( Pa)
    hhl        ,    & ! geometical height of half model levels        ( m )

! 2. external parameter fields                                        (unit)
! ----------------------------
    hd_mask_dcoeff_p, & ! 3D-domain mask for horizontal diffusion * dcoeff--
    hd_mask_dcoeff_t, & ! 3D-domain mask for horizontal diffusion * dcoeff--
    hd_mask_dcoeff_q, & ! 3D-domain mask for horizontal diffusion * dcoeff--
    hd_mask_dcoeff_u, & ! 3D-domain mask for horizontal diffusion * dcoeff--
    ofa_hdx    ,    & !
    ofa_hdy    ,    & !
    crlat      ,    & ! cosine of transformed latitude
    acrlat     ,    & ! 1 / ( crlat * radius of the earth )           ( 1/m )

! 3. prognostic variables                                             (unit)
! -----------------------
    u          ,    & ! zonal wind speed                              ( m/s )
    v          ,    & ! meridional wind speed                         ( m/s )
    w          ,    & ! vertical wind speed (defined on half levels)  ( m/s )
    t          ,    & ! temperature                                   (  k  )
    pp                ! deviation from the reference pressure         ( pa  )

! end of data_fields

!------------------------------------------------------------------------------

USE data_parallel,      ONLY :  &
    my_cart_id,&
    num_compute,     & ! number of compute PEs
    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains
    ldatatypes,      & ! if .TRUE.: use MPI-Datatypes for some communications
    ltime_barrier,   & ! if .TRUE.: use additional barriers for determining the
                       ! load-imbalance
    ncomm_type,      & ! type of communication
    my_cart_neigh,   & ! neighbors of this subdomain in the cartesian grid
    icomm_cart,      & ! communicator for the virtual cartesian topology
    iexch_req,       & ! stores the sends requests for the neighbor-exchange
                       ! that can be used by MPI_WAIT to identify the send
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    nexch_tag,       & ! tag to be used for MPI boundary exchange
                       !  (in calls to exchg_boundaries)
    sendbuf,         & ! sending buffer for boundary exchange:
                       ! 1-4 are used for sending, 5-8 are used for receiving
    isendbuflen        ! length of one column of sendbuf

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
    nnew,         & ! corresponds to ntstep + 1
    nnow,         & ! corresponds to ntstep (required for MESSy)

! 7. additional control variables
! -------------------------------
    lperi_x,      & ! if lartif_data=.TRUE.: periodic boundary conditions (.TRUE.) in x-dir.
                    ! or with Davies conditions (.FALSE.)
    lperi_y,      & ! if lartif_data=.TRUE.: periodic boundary conditions (.TRUE.) in y-dir.
                    ! or with Davies conditions (.FALSE.)
    l2dim,        & ! 2 dimensional runs
    ltime,        & ! detailled timings of the program are given
    l2tls,        & ! switch for 2-timelevel scheme
    l_cosmo_art,  & ! if .TRUE., run the COSMO_ART
    l_pollen,     & ! of pollen
                    ! correction factor for horizontal diffusion fluxes of
    hd_corr_u_bd   , & ! u,v,w      in boundary zone
    hd_corr_t_bd   , & ! t          in boundary zone
    hd_corr_trcr_bd, & ! tracers    in boundary zone
    hd_corr_p_bd   , & ! p          in boundary zone
    hd_corr_u_in   , & ! u,v,w      in domain
    hd_corr_t_in   , & ! t          in domain
    hd_corr_trcr_in, & ! tracers    in domain
    hd_corr_p_in   , & ! p          in domain
    hd_dhmax,     & ! maximum gridpoint height difference for applying
                    ! horizontal diffusion fluxes between them
    l_diff_Smag,  & ! use Smagorinsky-diffusion for u and v
    l_diff_cold_pools,& ! use targeted diffusion for cold pools
    thresh_cold_pool, & ! threshold for targeted diffusion for cold pools

! 12. controlling verbosity of debug output and miscellaneous
! -----------------------------------------------------------
    linit_fields    ! to initialize also local variables with a default value


!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 2. physical constants and related variables
! -------------------------------------------
    r_earth,      & ! mean radius of the earth (m)
    r_d             ! gas constant for dry air

!------------------------------------------------------------------------------

USE data_tracer      , ONLY :  T_DIFF_ID, T_DIFF_ON, T_MISSING

!------------------------------------------------------------------------------

USE src_tracer       , ONLY : trcr_get_ntrcr, trcr_get, trcr_meta_get,     &
                              trcr_errorstr

!------------------------------------------------------------------------------

USE environment      , ONLY :  exchg_boundaries, comm_barrier, model_abort
USE time_utilities   , ONLY :  get_timings, i_dyn_computations,            &
           i_barrier_waiting_dyn, i_communications_dyn, i_horizontal_diffusion
USE numeric_utilities, ONLY :  lap_4a, lap_4am, lap_4aml, lap_2

!==============================================================================

IMPLICIT NONE

!==============================================================================

CONTAINS

!==============================================================================

SUBROUTINE comp_hori_diff (itype)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine computes the horizontal diffusion for the atmospheric
!   variables. Depending on the type, different schemes are used for different
!   variables.
!
!   For itype=1/2: the diffusion for the prognostic variables u,v,w,pp,t,qv
!                  and qc at time level n+1 (nnew).
!
!------------------------------------------------------------------------------


! Subroutine arguments:
! --------------------
  INTEGER (KIND=iintegers), INTENT (IN) ::  &
    itype      ! type of horizontal diffusion used

! Local scalars:
! -----------------------------

  INTEGER (KIND=iintegers) ::  &
    i,  j,  k,           & ! Loop indices
    iztrcr,              & ! Index for tracer looping
    nztrcr_diff,         & ! Number of tracers undergoing horiz. diffusion
    kzdims(24)             ! Vertical dimensions for exchg_datatypes

! Local (automatic) arrays:
! -----------------------------
  REAL    (KIND=wp   )     ::  &
    ztha   (ie,je,ke),   & ! temperature deviation from reference
    zlap   (ie,je,ke),   & ! work array
    zcrlato(   je   ),   & !
    zcrlatu(   je   ),   & !
    zcrlavo(   je   ),   & !
    zcrlavu(   je   )      !

  REAL (KIND=wp),     POINTER :: &
    ztrcr  (:,:,:) => NULL()   ! tracer variable at tlev=nnew

  REAL    (KIND=wp   ),     ALLOCATABLE ::  &
     k_diff_Smag_u(:,:,:),  &  ! dim.-less coefficients for Smagorinsky diffusion
     k_diff_Smag_v(:,:,:)

  REAL    (KIND=wp   ),     ALLOCATABLE ::  &
    zlap_3d       (:,:,:)  , & ! Laplacian for 1 selected variable (temperature)
    zlap_3d_trcr  (:,:,:,:)    ! Laplacian for tracer variables

  INTEGER  (KIND=iintegers), ALLOCATABLE :: &
    iztrcr_diff(:)             ! Index list of tracers which undergo artificial
                               ! hyperdiffusion

  REAL (KIND=wp   ),     ALLOCATABLE :: &
    Tp(:,:,:)      ! temperature deviation from reference state T0
  REAL (KIND=wp   )     ::   &
    T_threshold    ! threshold to trigger targeted diffusion against cold pools
  INTEGER :: dk    ! number of levels above ground for the targeted diffusion

! For error handling
! ------------------
  INTEGER (KIND=iintegers) ::  izerror

  CHARACTER (LEN=255)      ::  yzerrmsg
  CHARACTER (LEN=25)       ::  yzroutine

! End of header
!==============================================================================

  izerror  = 0_iintegers
  yzerrmsg = '    '
  yzroutine= 'comp_hori_diff'
  kzdims(:)= 0_iintegers

  ALLOCATE (zlap_3d(ie,je,ke), STAT=izerror)

  ALLOCATE ( iztrcr_diff(trcr_get_ntrcr() ), STAT = izerror)
  iztrcr_diff = T_MISSING

  ! Meta-data handling for tracers
  CALL trcr_meta_get(izerror, T_DIFF_ID, iztrcr_diff)
  IF (izerror /= 0_iintegers) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  ! Identification of the tracers undergoing hori. diffusion
  nztrcr_diff = 0_iintegers
  DO  iztrcr = 1, trcr_get_ntrcr()
    IF ( iztrcr_diff(iztrcr) == T_DIFF_ON ) THEN
      nztrcr_diff = nztrcr_diff + 1_iintegers
      iztrcr_diff(nztrcr_diff) = iztrcr
    ENDIF
  ENDDO

  ! Memory allocation for variables used by the tracers
  ALLOCATE ( zlap_3d_trcr(ie, je, ke, nztrcr_diff), STAT=izerror)
  zlap_3d_trcr = 0.0_wp

!------------------------------------------------------------------------------
! Section 1: Boundary Exchange of prognostic variables
!------------------------------------------------------------------------------


  IF (ltime) CALL get_timings (i_dyn_computations, ntstep, dt, izerror)
  IF (ltime_barrier) THEN
    CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
    IF (ltime) CALL get_timings (i_barrier_waiting_dyn, ntstep, dt, izerror)
  ENDIF

  IF ( .NOT.l2tls ) THEN
    IF ( hd_corr_trcr_in /= 0.0_wp .OR. hd_corr_trcr_bd /= 0.0_wp ) THEN
      kzdims(1:24) =                                                        &
              (/ke,ke,ke1,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                 &
       (22+nnew,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, &
        ie, je, kzdims, jstartpar, jendpar, 2, nboundlines,                 &
        my_cart_neigh, lperi_x, lperi_y, l2dim,                             &
        1000+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,          &
        u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),     &
        pp(:,:,:,nnew) )

      ! Exchange tracers undergoing hori. diffusion
      DO iztrcr = 1, nztrcr_diff

        ! get pointer to tracer (at nnew)
        CALL trcr_get(izerror, iztrcr_diff(iztrcr), ptr_tlev=nnew, ptr=ztrcr)
        IF (izerror /= 0_iintegers) THEN
          yzerrmsg = trcr_errorstr(izerror)
          CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
        ENDIF 

        ! halo-update
        kzdims(1:24) =                                                        &
                (/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                 &
         ( 4     ,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, &
          ie, je, kzdims, jstartpar, jendpar, 2, nboundlines,                 &
          my_cart_neigh, lperi_x, lperi_y, l2dim,                             &
          2000+nexch_tag+iztrcr_diff(iztrcr), ldatatypes, ncomm_type,         &
          izerror, yzerrmsg, ztrcr(:,:,:) )

      ENDDO

    ELSE   ! hd_corr_q_in/bd = 0: only mass, no humidities or tracers
      kzdims(1:24) =                                                         &
              (/ke,ke,ke1,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                  &
       (22+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
        ie, je, kzdims, jstartpar, jendpar, 2, nboundlines,                  &
        my_cart_neigh, lperi_x, lperi_y, l2dim,                              &
        10000+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,          &
        u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),      &
        pp(:,:,:,nnew) )
    END IF
  ELSE   ! l2tls
    ! mass variables are exchanged in src_runge_kutta
    IF ( hd_corr_trcr_in /= 0.0_wp .OR. hd_corr_trcr_bd /= 0.0_wp ) THEN

      ! Exchange tracers undergoing hori. diffusion
      DO iztrcr = 1, nztrcr_diff

        ! get pointer to tracer (at nnew)
        CALL trcr_get(izerror, iztrcr_diff(iztrcr), ptr_tlev=nnew, ptr=ztrcr)
        IF ( izerror /= 0_iintegers ) THEN
          yzerrmsg = trcr_errorstr(izerror)
          CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
        ENDIF

        ! halo-update
        kzdims(1:24) =                                                       &
                (/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                &
         ( 4     ,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,&
          ie, je, kzdims, jstartpar, jendpar, 2, nboundlines,                &
          my_cart_neigh, lperi_x, lperi_y, l2dim,                            &
          12000+nexch_tag+iztrcr_diff(iztrcr), ldatatypes, ncomm_type,       &
          izerror, yzerrmsg, ztrcr(:,:,:) )

      ENDDO

    ENDIF
  ENDIF

  IF (ltime) CALL get_timings (i_communications_dyn, ntstep, dt, izerror)

!------------------------------------------------------------------------------
! Section 2: Preparations
!------------------------------------------------------------------------------

  DO j = 2, je-1
    zcrlato(j) = crlat(j  ,2)/crlat(j,1)
    zcrlatu(j) = crlat(j-1,2)/crlat(j,1)
    zcrlavo(j) = crlat(j+1,1)/crlat(j,2)
    zcrlavu(j) = crlat(j  ,1)/crlat(j,2)
  ENDDO
  zcrlato(1)   = crlat(1  ,2)/crlat(1,1)

!------------------------------------------------------------------------------
! Section 3: Now do the horizontal diffusion
!------------------------------------------------------------------------------

  !$ser savepoint HorizontalDiffusionUnittest.DoStep-in LargeTimeStep=ntstep
  !$ser data t=t(:,:,:,nnew)-t0 pp=pp(:,:,:,nnew) w=w(:,:,:,nnew)  &
  !$ser&     u=u(:,:,:,nnew) v=v(:,:,:,nnew)
  !$ser tracer %all@nnew

  !$ser savepoint HorizontalDiffusionUnittest.ColdPool-in LargeTimeStep=ntstep
  !$ser data t=t(:,:,:,nnew)-t0

  IF ( l_diff_cold_pools ) THEN

    ! extra treatment for 'cold pools'
    ! which can occur in narrow valleys

    ALLOCATE( Tp(ie,je,ke), STAT=izerror )

    dk = 3    ! number of levels above the surface which are treated

    DO k=ke-dk, ke
      DO j=1, je
        DO i=1, ie
          Tp(i,j,k) = T(i,j,k,nnew) - T0(i,j,k)
        END DO
      END DO
    END DO

    CALL targeted_diffusion_cold_pools( Tp(:,:,:), thresh_cold_pool, dk )

    DO k=ke-dk, ke
      DO j=1, je
        DO i=1, ie
          T(i,j,k,nnew) = Tp(i,j,k) + T0(i,j,k)
        END DO
      END DO
    END DO

    DEALLOCATE( Tp )

  END IF

  !MB: Bug fix by Andrea Arteaga (04.11.2014): shift of calc of ztha ...
  ztha (:,:,:) = t(:,:,:,nnew) - p0(:,:,:)/(r_d*rho0(:,:,:))

  !$ser savepoint HorizontalDiffusionUnittest.ColdPool-out LargeTimeStep=ntstep
  !$ser data t=t(:,:,:,nnew)-t0

  SELECT CASE(itype)

!------------------------------------------------------------------------------
! Section 3a: First type of horizontal diffusion
!------------------------------------------------------------------------------

  CASE(1)

    DO  k = 1, ke

      ! Apply Laplace-operator to all fields which are defined at the
      ! scalar grid point and update the variables accordingly

      CALL lap_4a(t (:,:,k,nnew), ztha(:,:,k),    zlap, zcrlato, zcrlatu,   &
                  hd_mask_dcoeff_t(:,:,k), ie, je,                          &
                  istart, iend, jstart, jend)

      IF ( hd_corr_trcr_in /= 0.0_wp .OR. hd_corr_trcr_bd /= 0.0_wp ) THEN

        ! Laplace-operator O(4) on the tracers undergoing hori. diffusion
        DO iztrcr = 1, nztrcr_diff

          ! get pointer to tracer (at nnew)
          CALL trcr_get(izerror, iztrcr_diff(iztrcr), ptr_tlev=nnew, ptr=ztrcr)
          IF ( izerror /= 0_iintegers ) THEN
            yzerrmsg = trcr_errorstr(izerror)
            CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
          ENDIF

          ! compute diffusion operator
          CALL lap_4a( ztrcr(:,:,k), ztrcr(:,:,k), zlap, zcrlato, zcrlatu,    &
                       hd_mask_dcoeff_q(:,:,k), ie, je, istart, iend,         &
                       jstart, jend )
        ENDDO

      ENDIF

      CALL lap_4a(pp(:,:,k,nnew), pp(:,:,k,nnew), zlap, zcrlato, zcrlatu,   &
                  hd_mask_dcoeff_p(:,:,k), ie, je,                          &
                  istart, iend, jstart, jend)

      ! Apply Laplace operator to the horizontal (vertical) wind components which are
      ! defined at the u- and v- (w-) gridpoint locations of the c-grid
      ! and update the variables accordingly. Also, the diffusion coefficient
      ! is reduced by a constant factor hd_corr_u.

      IF (k > 1) THEN
        CALL lap_4a(w(:,:,k,nnew), w(:,:,k,nnew), zlap, zcrlato, zcrlatu,   &
                    hd_mask_dcoeff_u(:,:,k), ie, je,                        &
                    istart, iend, jstart, jend)
      ENDIF

      CALL lap_4a(u(:,:,k,nnew), u(:,:,k,nnew), zlap, zcrlato, zcrlatu,     &
                  hd_mask_dcoeff_u(:,:,k), ie, je,                          &
                  istartu, iendu, jstartu, jendu)

      CALL lap_4a(v(:,:,k,nnew), v(:,:,k,nnew), zlap, zcrlavo, zcrlavu,     &
                  hd_mask_dcoeff_u(:,:,k), ie, je,                          &
                  istartv, iendv, jstartv, jendv)

    ENDDO

!------------------------------------------------------------------------------
! Section 3b: Horizontal diffusion with orographic limiter
!------------------------------------------------------------------------------

  CASE(2)

    ! For t, qv and qc, a monotonic operator with orographic limiting
    ! is applied. Also, the diffusion coefficient is reduced by a
    ! constant factor hd_corr_t/q.

    ! The monotonic limiter requires the Laplacian on the total domain.
    ! Thus, a data exchange between processors is required.
    IF ( hd_corr_trcr_in /= 0.0_wp .OR. hd_corr_trcr_bd /= 0.0_wp ) THEN

      ! Laplace-operator O(2) on the tracers undergoing hori. diffusion
      DO iztrcr = 1, nztrcr_diff

        ! get pointer to tracer (at nnew)
        CALL trcr_get(izerror, iztrcr_diff(iztrcr), ptr_tlev=nnew, ptr=ztrcr)
        IF ( izerror /= 0_iintegers ) THEN
          yzerrmsg = trcr_errorstr(izerror)
          CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
        ENDIF

        ! compute Laplacian
        CALL lap_2( ztrcr(:,:,:), zlap_3d_trcr(:,:,:,iztrcr),&
                    zcrlato, zcrlatu, ie, je, ke )

      ENDDO

    ENDIF

    CALL lap_2 (ztha(:,:,:), zlap_3d(:,:,:), zcrlato, zcrlatu, ie, je, ke )

    IF (ltime) CALL get_timings (i_horizontal_diffusion, ntstep, dt, izerror)
    IF (ltime_barrier) THEN
      CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
      IF (ltime) CALL get_timings (i_barrier_waiting_dyn, ntstep,dt, izerror)
    ENDIF

    IF ( hd_corr_trcr_in /= 0.0_wp .OR. hd_corr_trcr_bd /= 0.0_wp ) THEN
      kzdims(1:24)=(/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                  &
         ( 4, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,      &
          ie, je, kzdims, jstartpar, jendpar, 2, nboundlines,                &
          my_cart_neigh, lperi_x, lperi_y, l2dim,                            &
          15000+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,        &
          zlap_3d(:,:,:) )

        ! Exchange tracers undergoing hori. diffusion
        DO iztrcr = 1, nztrcr_diff

          ! halo-update
          kzdims(1:24)=(/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
          CALL exchg_boundaries                                              &
             ( 4      , sendbuf, isendbuflen, imp_reals, icomm_cart,         &
               num_compute, ie, je, kzdims, jstartpar, jendpar, 2,           &
               nboundlines, my_cart_neigh, lperi_x, lperi_y, l2dim,          &
               16000+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,   &
               zlap_3d_trcr(:,:,:,iztrcr) )

        ENDDO

    ELSE
      kzdims(1:24)=(/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                  &
         ( 4, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,      &
          ie, je, kzdims, jstartpar, jendpar, 2, nboundlines,                &
          my_cart_neigh, lperi_x, lperi_y, l2dim,                            &
          22000+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,        &
          zlap_3d(:,:,:) )
    END IF
    IF (ltime) CALL get_timings (i_communications_dyn, ntstep, dt, izerror)

    CALL lap_4aml(t (:,:,:,nnew), ztha(:,:,:), zlap_3d(:,:,:),             &
                  ofa_hdx, ofa_hdy, zcrlato, zcrlatu,                      &
                  hd_mask_dcoeff_t(:,:,:),                                 &
                  ie, je, ke, istart, iend, jstart, jend )

    IF ( hd_corr_trcr_in /= 0.0_wp .OR. hd_corr_trcr_bd /= 0.0_wp ) THEN

      !Laplace-operator O(4) with monotonic flux limiter on tracers
      !undergoing hori. diffusion
      DO iztrcr = 1, nztrcr_diff

        ! get pointer to tracer (at nnew)
        CALL trcr_get(izerror, iztrcr_diff(iztrcr), ptr_tlev=nnew, ptr=ztrcr)
        IF ( izerror /= 0_iintegers ) THEN
          yzerrmsg = trcr_errorstr(izerror)
          CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
        ENDIF

        ! compute diffusion operator
        CALL lap_4aml( ztrcr(:,:,:), ztrcr(:,:,:),                            &
                       zlap_3d_trcr(:,:,:,iztrcr),                            &
                       ofa_hdx, ofa_hdy, zcrlato, zcrlatu,                    &
                       hd_mask_dcoeff_q(:,:,:), ie, je, ke, istart, iend,     &
                       jstart, jend )

      ENDDO

    ENDIF

    CALL lap_4am(pp(:,:,:,nnew), pp(:,:,:,nnew), zlap, zcrlato, zcrlatu,   &
                 hd_mask_dcoeff_p(:,:,:),                                  &
                 ie, je, ke, istart, iend, jstart, jend, 1, ke)

    ! Apply Laplace operator to the horizontal (vertical) wind components which are
    ! defined at the u- and v- (w-) gridpoint locations of the c-grid
    ! and update the variables accordingly. Also, the diffusion coefficient
    ! is reduced by a constant factor hd_corr_u.

    CALL lap_4am(w(:,:,:,nnew),  w(:,:,:,nnew), zlap, zcrlato, zcrlatu,   &
                 hd_mask_dcoeff_u(:,:,:),                                 &
                 ie, je, ke, istart, iend, jstart, jend, 2, ke)

    CALL lap_4am(u(:,:,:,nnew),  u(:,:,:,nnew), zlap, zcrlato, zcrlatu,    &
                 hd_mask_dcoeff_u(:,:,:),                                  &
                 ie, je, ke, istartu, iendu, jstartu, jendu, 1, ke)

    CALL lap_4am(v(:,:,:,nnew),  v(:,:,:,nnew), zlap, zcrlavo, zcrlavu,    &
                 hd_mask_dcoeff_u(:,:,:),                                  &
                 ie, je, ke, istartv, iendv, jstartv, jendv, 1, ke)

  END SELECT

  kzdims(1:24)=(/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
  CALL exchg_boundaries                                                    &
     ( 25+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
       ie, je, kzdims, jstartpar, jendpar, 2, nboundlines, my_cart_neigh,  &
       lperi_x, lperi_y, l2dim,                                            &
       18000+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,         &
       u(:,:,:,nnew), v(:,:,:,nnew) )


  IF ( l_diff_Smag ) THEN

    ! Smagorinsky-diffusion

    ALLOCATE ( k_diff_Smag_u (ie,je,ke), STAT=izerror)
    ALLOCATE ( k_diff_Smag_v (ie,je,ke), STAT=izerror)

    CALL smagorinsky_coeff( k_diff_Smag_u(:,:,:),    k_diff_Smag_v(:,:,:),     &
      &                     hd_mask_dcoeff_u(:,:,:), hd_mask_dcoeff_u(:,:,:),  &
      &                     u(:,:,:,nnew),           v(:,:,:,nnew),      dt )

    CALL lap_2( u(:,:,:,nnew), zlap(:,:,:), zcrlato, zcrlatu, ie, je, ke )
    u(:,:,:,nnew) = u(:,:,:,nnew) + k_diff_Smag_u(:,:,:) * zlap(:,:,:)

    !MB: Bug fix von Andrea Arteaga (13.11.2014)
    CALL lap_2( v(:,:,:,nnew), zlap(:,:,:), zcrlavo, zcrlavu, ie, je, ke )
    v(:,:,:,nnew) = v(:,:,:,nnew) + k_diff_Smag_v(:,:,:) * zlap(:,:,:)

    kzdims(1:24)=(/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                  &
     ( 25+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
       ie, je, kzdims, jstartpar, jendpar, 2, nboundlines, my_cart_neigh,  &
       lperi_x, lperi_y, l2dim,                                            &
       19000+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,         &
       u(:,:,:,nnew), v(:,:,:,nnew) )


    DEALLOCATE ( k_diff_Smag_u )
    DEALLOCATE ( k_diff_Smag_v )

  END IF

  DEALLOCATE (zlap_3d)

  IF ( ALLOCATED(zlap_3d_trcr) ) THEN
    DEALLOCATE ( zlap_3d_trcr )
  ENDIF

  IF (ltime) CALL get_timings (i_horizontal_diffusion, ntstep, dt, izerror)

!------------------------------------------------------------------------------
!  End of the module procedure comp_hori_diff
!------------------------------------------------------------------------------

END SUBROUTINE comp_hori_diff

!==============================================================================

SUBROUTINE smagorinsky_coeff( k_diff_Smag_u, k_diff_Smag_v,  &
        hd_mask_dcoeff_u,  hd_mask_dcoeff_v,                 &
        u, v, dt)

!------------------------------------------------------------------------------
! Description:
!   calculate the diffusion coefficient for the
!   horizontal nonlinear Smagorinsky diffusion
!
! Input:
!   u(:,:,:), v(:,:,:): for the calculation of tension and shear strain
!
!   hd_mask_dcoeff_u, hd_mask_dcoeff_v: are used to approximately
!        avoid double counting with the 4th order diffusion
!
! Output:
!   k_diff_Smag_u, k_diff_Smag_v: dimensionless diffusion
!      coefficients at the u- and v-position
!
!------------------------------------------------------------------------------

  REAL (KIND=wp),     INTENT(IN)  :: hd_mask_dcoeff_u(ie, je, ke)
  REAL (KIND=wp),     INTENT(IN)  :: hd_mask_dcoeff_v(ie, je, ke)

  REAL (KIND=wp),     INTENT(IN)  :: u(ie, je, ke)
  REAL (KIND=wp),     INTENT(IN)  :: v(ie, je, ke)
  REAL (KIND=wp),     INTENT(IN)  :: dt

  REAL (KIND=wp),     INTENT(OUT) :: k_diff_Smag_u(ie, je, ke)
  REAL (KIND=wp),     INTENT(OUT) :: k_diff_Smag_v(ie, je, ke)

  REAL (KIND=wp),     ALLOCATABLE :: T_sqr_s (:,:) ! tension^2 (at scalar position)
  REAL (KIND=wp),     ALLOCATABLE :: S_sqr_uv(:,:) ! shear^2   (at uv position)

  REAL (KIND=wp)     :: frac_1_dx, frac_1_dy
  REAL (KIND=wp)     :: T_s, S_uv
  REAL (KIND=wp)     :: T_sqr_u, T_sqr_v   ! tension^2 at u- / v- position
  REAL (KIND=wp)     :: S_sqr_u, S_sqr_v   ! shear^2   at u- / v- position
  REAL (KIND=wp)     :: c_Smag        ! some similarity with the Smagorinsky-constant
  REAL (KIND=wp)     :: tau_Smag      ! time-scale
  REAL (KIND=wp)     :: weight_K_4th  ! avoid double counting with the 4th order diffusion

  REAL (KIND=wp)     :: mean_u, mean_level_u, maxi_u
  REAL (KIND=wp)     :: mean_v, mean_level_v, maxi_v

  INTEGER :: i,j,k

  ALLOCATE( T_sqr_s (ie,je) )
  ALLOCATE( S_sqr_uv(ie,je) )

  frac_1_dy = eddlat / r_earth

  ! dimensionless parameter (approximately the Smagorinsky constant):
  !c_Smag = 0.1_wp
  c_Smag = 0.03_wp
  !c_Smag = 0.01_wp

  ! a weight to avoid double counting with the 4th order diffusion:
  weight_K_4th = 0.5_wp

  tau_Smag = c_Smag * dt    ! time-scale

  DO k=1, ke

    DO j=2, je-1
      frac_1_dx = acrlat(j,1) * eddlon

      DO i=2, ie-1

        ! tension T at scalar position:
        T_s   = ( u(i,j,k) - u(i-1,j,k) ) * frac_1_dx  &
          &   - ( v(i,j,k) - v(i,j-1,k) ) * frac_1_dy

        ! shear S at uv-position (i+1/2, j+1/2):
        S_uv  = ( u(i,j+1,k) - u(i,j,k) ) * frac_1_dy  &
          &   + ( v(i+1,j,k) - v(i,j,k) ) * frac_1_dx

        ! remark: for very large model areas, metric correction terms should be
        ! included (Smagorinsky, 1993)

        ! calculate squares BEFORE averaging to u- / v- grid positions:
        ! (to prevent cancellation of neighbouring positive and negative values):
        T_sqr_s (i,j) = T_s  * T_s
        S_sqr_uv(i,j) = S_uv * S_uv

      END DO
    END DO

    T_sqr_s ( 1,2:je-1) = T_sqr_s (   2,2:je-1)
    T_sqr_s (ie,2:je-1) = T_sqr_s (ie-1,2:je-1)

    S_sqr_uv( 1,2:je-1) = S_sqr_uv(   2,2:je-1)
    S_sqr_uv(ie,2:je-1) = S_sqr_uv(ie-1,2:je-1)

    T_sqr_s (2:ie-1, 1) = T_sqr_s (2:ie-1,   2)
    T_sqr_s (2:ie-1,je) = T_sqr_s (2:ie-1,je-1)

    S_sqr_uv(2:ie-1, 1) = S_sqr_uv(2:ie-1,   2)
    S_sqr_uv(2:ie-1,je) = S_sqr_uv(2:ie-1,je-1)

    ! (remark: corner points (1,1,:), (1,je,:), (ie,1,:), (ie,je,:)
    ! are not set, but they are never used in the following)

    DO j=2, je-1
      DO i=2, ie-1

        T_sqr_u = 0.5_wp * ( T_sqr_s(i+1,j) + T_sqr_s(i,j ) )
        T_sqr_v = 0.5_wp * ( T_sqr_s(i,j+1) + T_sqr_s(i,j ) )

        S_sqr_u = 0.5_wp * ( S_sqr_uv(i,j) + S_sqr_uv(i,j-1) )
        S_sqr_v = 0.5_wp * ( S_sqr_uv(i,j) + S_sqr_uv(i-1,j) )

        k_diff_Smag_u(i,j,k) = tau_Smag * SQRT( T_sqr_u + S_sqr_u )
        k_diff_Smag_v(i,j,k) = tau_Smag * SQRT( T_sqr_v + S_sqr_v )


        ! approximately avoid double counting with the 4th order
        ! horizontal diffusion:

        k_diff_Smag_u(i,j,k) = k_diff_Smag_u(i,j,k)   &
          &                  - weight_K_4th * hd_mask_dcoeff_u(i,j,k)

        k_diff_Smag_v(i,j,k) = k_diff_Smag_v(i,j,k)   &
          &                  - weight_K_4th * hd_mask_dcoeff_v(i,j,k)

        ! clip negative values:
        k_diff_Smag_u(i,j,k) = MAX( 0.0_wp, k_diff_Smag_u(i,j,k) )
        k_diff_Smag_v(i,j,k) = MAX( 0.0_wp, k_diff_Smag_v(i,j,k) )


        ! avoid numerical instability of the diffusion:
        k_diff_Smag_u(i,j,k) = MIN( 0.5_wp, k_diff_Smag_u(i,j,k) )
        k_diff_Smag_v(i,j,k) = MIN( 0.5_wp, k_diff_Smag_v(i,j,k) )

      END DO
    END DO

  END DO

  ! set boundary values:
  k_diff_Smag_u(1 ,: ,:) = k_diff_Smag_u(2   ,:   ,:)
  k_diff_Smag_u(ie,: ,:) = k_diff_Smag_u(ie-1,:   ,:)

  k_diff_Smag_v(1 ,: ,:) = k_diff_Smag_v(2   ,:   ,:)
  k_diff_Smag_v(ie,: ,:) = k_diff_Smag_v(ie-1,:   ,:)

  k_diff_Smag_u(:, 1 ,:) = k_diff_Smag_u(:   ,2   ,:)
  k_diff_Smag_u(:, je,:) = k_diff_Smag_u(:   ,je-1,:)

  k_diff_Smag_v(:, 1 ,:) = k_diff_Smag_v(:   ,2   ,:)
  k_diff_Smag_v(:, je,:) = k_diff_Smag_v(:   ,je-1,:)

  ! remark: a proper (i.e. rotation invariant) treatment of
  ! corner points would additionally need:
  ! k_diff_Smag_u(1 ,1 ,:) = k_diff_Smag_u(2   ,2   ,:)
  ! k_diff_Smag_v(1 ,1 ,:) = k_diff_Smag_v(2   ,2   ,:)
  ! k_diff_Smag_u(1 ,je,:) = k_diff_Smag_u(2   ,je-1,:)
  ! k_diff_Smag_v(1 ,je,:) = k_diff_Smag_v(2   ,je-1,:)
  ! k_diff_Smag_u(ie,1 ,:) = k_diff_Smag_u(ie-1,2   ,:)
  ! k_diff_Smag_v(ie,1 ,:) = k_diff_Smag_v(ie-1,2   ,:)
  ! k_diff_Smag_u(ie,je,:) = k_diff_Smag_u(ie-1,je-1,:)
  ! k_diff_Smag_v(ie,je,:) = k_diff_Smag_v(ie-1,je_1,:)
  ! but this seems not to be very important

  !tracer(:,:,:,1) = k_diff_Smag(:,:,:)


  ! some statistics output:
  IF ( .FALSE. ) THEN
    mean_u = 0.0_wp
    mean_v = 0.0_wp

    DO k=1, ke
      mean_level_u = 0.0_wp
      mean_level_v = 0.0_wp

      DO j=1,je
        DO i=1,ie
          mean_level_u = mean_level_u + k_diff_Smag_u(i,j,k)
          mean_level_v = mean_level_v + k_diff_Smag_v(i,j,k)
        END DO
      END DO

      mean_u = mean_u + mean_level_u
      mean_v = mean_v + mean_level_v
    END DO

    mean_u = mean_u / ( ie * je * ke )
    mean_v = mean_v / ( ie * je * ke )

    maxi_u = MAXVAL ( k_diff_Smag_u(:,:,:) )
    maxi_v = MAXVAL ( k_diff_Smag_v(:,:,:) )

    WRITE(*,'(A,2(A,F12.6,A,F12.6))') "k_diff_Smag: ",   &
      &     "mean_u=", mean_u, ", max_u=", maxi_u,       &
      &   ", mean_v=", mean_v, ", max_v=", maxi_v
  END IF

  DEALLOCATE( T_sqr_s  )
  DEALLOCATE( S_sqr_uv )

END SUBROUTINE smagorinsky_coeff

!==============================================================================

SUBROUTINE targeted_diffusion_cold_pools( Tp, T_threshold, dk )

  ! avoid unrealistic undershootings in a field Tp 
  ! (by more than T_threshold compared to the average value of the surrounding) 
  ! near the bottom (i.e. in levels ke, ke-1, ..., ke-dk)
  ! by targeted diffusion.

  ! It is assumed that this occurs only in a few grid points.
  ! Therefore the targeted diffusion is not implemtend as a field
  ! operation but only on selected grid points.

  REAL (KIND=wp),     INTENT(INOUT) :: Tp(ie, je, ke)
  REAL (KIND=wp),     INTENT(IN)    :: T_threshold
  INTEGER, INTENT(IN)               :: dk  

  INTEGER, ALLOCATABLE :: i_liste(:)
  INTEGER, ALLOCATABLE :: j_liste(:)
  INTEGER, ALLOCATABLE :: k_liste(:)

  INTEGER :: nmbr_points
  INTEGER :: i,j,k,l

  REAL (KIND=wp) :: alpha_x, alpha_y, c1
  REAL (KIND=wp) :: Tp_environ
  REAL (KIND=wp), ALLOCATABLE :: T_new(:)

  LOGICAL :: is_add_exchange_necessary

  INTEGER (KIND=iintegers) ::  &
    kzdims(24)          ! Vertical dimensions for exchg_datatypes

  INTEGER :: ibegin, istop
  INTEGER :: jbegin, jstop

  INTEGER (KIND=iintegers) ::  izerror
  CHARACTER (LEN=80)       ::  yzerrmsg


!MB WRITE(*,*) "Subr. [targeted_diffusion_cold_pools] ..."

  ALLOCATE( i_liste(ie*je*(dk+1)) )
  ALLOCATE( j_liste(ie*je*(dk+1)) )
  ALLOCATE( k_liste(ie*je*(dk+1)) )
  ALLOCATE( T_new  (ie*je*(dk+1)) )

  ! diffusion 'constants':

  alpha_x = 0.125_wp  ! dim.less diffusion coeff. x-direction
  alpha_y = 0.125_wp  ! dim.less diffusion coeff. y-direction

  c1 = 1.0_wp - 2.0_wp * ( alpha_x + alpha_y)

  ! find critical grid points

  nmbr_points = 0

  ! western boundary:
  IF ( my_cart_neigh(1) == -1 ) THEN
    ibegin = istart
  ELSE 
    ibegin = MAX(1,istart-2)
  END IF

  ! eastern boundary
  IF ( my_cart_neigh(3) == -1 ) THEN
    istop  = iend
  ELSE
    istop  = MIN(ie,iend+2)
  END IF

  ! northern boundary
  IF ( my_cart_neigh(2) == -1 ) THEN
    jstop = jend
  ELSE
    jstop = MIN(je,jend+2)
  END IF

  ! southern boundary
  IF ( my_cart_neigh(4) == -1 ) THEN
    jbegin = jstart
  ELSE
    jbegin = MAX(1,jstart-2)
  END IF


  DO k= ke-dk, ke
    DO j=jbegin, jstop
      DO i=ibegin, istop

        ! remark: the index limits shall guarantee reproducibility
        ! (= independence of domain decomposition)
        ! for the subsequent 4th order artificial diffusion
        ! without an additional exchange of Tp.

        !Tp_environ = MIN( Tp(i+1,j  ,k), &
        !  &               Tp(i-1,j  ,k), &
        !  &               Tp(i  ,j+1,k), &
        !  &               Tp(i  ,j-1,k) )

        Tp_environ = 0.25_wp * ( Tp(i+1,j  ,k)   &
          &                    + Tp(i-1,j  ,k)   &
          &                    + Tp(i  ,j+1,k)   &
          &                    + Tp(i  ,j-1,k) )

        !MB: Vorschlag von Oliver Fuhrer, Andrea Arteaga (04.11.2014):
        IF ( ABS(Tp(i,j,k) - Tp_environ) > T_threshold ) THEN
          nmbr_points = nmbr_points + 1
          i_liste(nmbr_points) = i
          j_liste(nmbr_points) = j
          k_liste(nmbr_points) = k
        END IF

      END DO
    END DO
  END DO

  ! targeted diffusion

  if ( nboundlines >= 3 ) THEN
    is_add_exchange_necessary = .FALSE.
    ! Two additional rows and columns are calculated above.
    ! Therefore, no exchange is necessary for the following
    ! 4th order horizontal diffusion.
  ELSE
    is_add_exchange_necessary = .TRUE.
    ! Is only necessary, if a detected grid point lies at the border
    ! (see commented if-statement below).
    ! However, to guarantee reproducibility for any domain 
    ! decomposition, an exchange is forced here.
  END IF
    
  DO l=1, nmbr_points
    i = i_liste(l)
    j = j_liste(l)
    k = k_liste(l)

    T_new(l) =    c1 * Tp(i,j,k)                       &
      &  + alpha_x * ( Tp(i-1,j  ,k) + Tp(i+1,j  ,k) ) &
      &  + alpha_y * ( Tp(i  ,j-1,k) + Tp(i  ,j+1,k) )

    !IF (    ( i == istart ) .OR. ( i == iend )         &
    !  & .OR.( j == jstart ) .OR. ( j == jend ) ) THEN
    !  ! grid point lies at the border of the computational domain
    !  ! additional exchange is necessary
    !  is_add_exchange_necessary = .TRUE.
    ! ! IMPORTANT: is_add_exchange_necessary must be globally exchanged!
    !END IF

  END DO

  DO l=1, nmbr_points
    i = i_liste(l)
    j = j_liste(l)
    k = k_liste(l)

    Tp(i,j,k) = T_new(l)

  END DO

  IF ( is_add_exchange_necessary ) THEN

    kzdims(1:24)=(/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                 &
      ( 4, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,      &
      ie, je, kzdims, jstartpar, jendpar, 2, nboundlines, my_cart_neigh,  &
      lperi_x, lperi_y, l2dim,                                            &
      18000+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,         &
      Tp(:,:,:) )

  END IF

  !MB>>
  ! Statistik
!MB  WRITE(*,'(A,I6)') "  Anzahl Gitterpunkte: ", nmbr_points
!MB  DO l=1, MIN( nmbr_points, 10)
!MB    i = i_liste(l)
!MB    j = j_liste(l)
!MB    k = k_liste(l)
!MB
!MB    WRITE(*,'(4I5,F13.6,A,4F13.6)')  my_cart_id, i, j, k,    &
!MB       &  Tp(i,j,k), " -- ", Tp(i-1,j,k), Tp(i+1,j,k), Tp(i,j-1,k), Tp(i,j+1,k)
!MB  END DO

  !PRINT*, Tp(144,45,ke), Tp(143,45,ke), Tp(145,45,ke), Tp(144,44,ke), Tp(144,46,ke)
  !MB<<

  DEALLOCATE( i_liste, j_liste, k_liste, T_new )
END SUBROUTINE targeted_diffusion_cold_pools

!==============================================================================

END MODULE hori_diffusion
