!+ Source module for computing slow tendencies of Runge-Kutta dynamics
!------------------------------------------------------------------------------

MODULE src_slow_tendencies_rk

!------------------------------------------------------------------------------
!
! Description:
!   This module contains subroutines which compute the slow tendencies for
!   the Runge-Kutta version.
!
!   These routines have been in module src_runge_kutta.f90 before, which 
!   has been splitted.
!
!   Routines currenctly included are:
!    - complete_tendencies_init
!    - complete_tendencies_uvwtpp
!    - complete_tendencies_uvwtpp_eva
!    - compute_moisture_divergence
!    - complete_tendencies_tke
!    - complete_tendencies_trcr
!    - explicit_horizontal_diffusion
!    - implicit_vert_diffusion_uvwt
!
! Current Code Owner: DWD, Michael Baldauf
!  phone:  +49  69  8062 2733
!  fax:    +49  69  8062 3721
!  email:  michael.baldauf@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 3.21       2006/12/04 Ulrich Schaettler
!  Initial release
! V3_23        2007/03/30 Jochen Foerstner
!  Corrections for l3dturb_metr
! V4_5         2008/09/10 Guenther Zaengl
!  Possibility to switch off the surface momentum fluxes in idealized cases
!  by setting lfreeslip_sfc=.TRUE.
! V4_9         2009/07/16 Ulrich Schaettler, Christian Bollmann (NEC)
!  Implemented NEC ON_ADB Directives; changed some j- and k-loops
! V4_10        2009/09/11 Ulrich Schaettler, Christian Bollmann (NEC)
!  Changed syntax of on_adb option for NEC because of new compiler
! V4_11        2009/11/30 Guenther Zaengl
!  Implemented use of cloud ice tendencies from convection scheme
! V4_12        2010/05/11 Michael Baldauf
!  Replace local variables zMlambda_s, zMphi_s by new metric coefficients
!  dzeta_dlam, dzeta_dphi, resp.
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_15        2010/11/19 Michael Baldauf
!  new SUBROUTINE complete_tend_uvwtpp_CN3Crow
! V4_17        2011/02/24 Ulrich Blahak
!  Eliminated lfreeslip_surf (free-slip BC and/or no-surface-heat/moisture-flux
!    conditions can now be imposed by switches lnosurffluxes_m/h in namelist ARTIFCTL);
!  fixed use of lkge2_prog_tke in horizontal turbulent diffusion of TKE;
!  fixed various bugs in loop index boundaries in horizontal turbulent diffusion;
!  correction for ztch: introduced lower velocity limit vel_min (COMMENTED OUT FOR NOW)
! V4_18        2011/05/26 Ulrich Schaettler
!  Implementation of new SR complete_tendencies_cosmo_art, complete_tendencies_pollen
!  for COSMO-ART and Pollen (Christoph Knote)
!  Added comments for ztvb, ztch regarding the use of t_g instead of t_s and
!   the clipping of zvbke by vel_min. These things have to be tested further
!   and made consistent in the rest of the model in the future. (by Uli Blahak)
! V4_20        2011/08/31 Ulrich Schaettler
!  Activated lower velocity limit vel_min
! V4_23        2012/05/10 Ulrich Schaettler, Oliver Fuhrer
!  Use fields sqrtg_r_*, dzeta_* from new module grid_metrics_utilities
!  Introduced new SR compute_moisture_divergence and splitted the computations
!   for dqvdt from SR complete_tendencies_qvqcqi_tke (Oliver Fuhrer)
!  Removed field qvt_diff (Oliver Fuhrer)
!  Bug fix for vertical diffusion for cloud ice in SR complete_tendencies_qvqcqi_tke;
!    it must only be called, if ltur=.TRUE.  (Oliver Fuhrer)
! V4_24        2012/06/22 Oliver Fuhrer
!  Order of summations for humidity tracer qvtens in advection: was different 
!  than for all the other humidity tracers and has been adapted now
!    This changes the results because of numerical reasons
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer
!  Replaced qx-variables by using them from the tracer module
! V4_26        2012/12/06 Anne Roches
!  Changes and technical adaptations to the tracer handling
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  MESSy interface introduced
! V4_28        2013/07/12 KIT, Ulrich Schaettler
!  Changes to adapt COSMO-ART to new tracer module: all dependencies to
!  COSMOART and POLLEN deleted, because this is now handled by the tracer module
! V4_29        2013/10/04 Ulrich Blahak, Astrid Kerkweg
!  Bugfixes in explicit_horizontal_diffusion (3D-Turbulence):
!   - some loop index bounds were too small by 1 gridpoint
!   - corrected wrong metrical term in ztaud23
!   (This will change results only if l3dturb=.TRUE. and l3dturb_metr=.TRUE.)
!  Unification of MESSy interfaces and COSMO Tracer structure (AK)
! V5_1         2014-11-28 Michael Baldauf, Ulrich Blahak, Matthias Raschendorfer
!                         Oliver Fuhrer, Anne Roches, Lucio Torrisi
!  Introduction of a new variable p0ref_recip (for 1/p0ref) (MB)
!  Adaptations for lprog_tke=.true. in order to apply
!   explicit_horizontal_diffusion also for the
!   case itype_turb=3 (standarad turbulence according Raschendorfer).
!  Bug fixes for some cut-and-paste errors in explicit_horizontal_diffusion (l3dturb):
!   - in the advection parts for TKE (wrong indices, wrong grid factors)
!   - in the 2.nd flux components (wrong indices for dzeta_dphi)
!   - eliminated doubled tracer-loop in the uppermost level (k=1) in the 
!     calculation of ztrcr_tens
!  Replaced ireals by wp (working precision) (OF)
!  Removed hacks for the tracers (AR)
!  Stochastic perturbation of physics tendencies (SPPT) added for tracers. (LT)
!  Use t_g instead of t_s for flux calculations (UB)
! V5_3         2015-10-09 Michael Baldauf, Ulrich Blahak, Oliver Fuhrer
!  Complete redesign of the 3D diffusion subroutine 'explicit_horizontal_diffusion'
!    (and partially of complete_tendencies_trcr)
!   and a suitable subroutine 'complete_tendencies_uvwt_3D':
!   now all terms, which can be treated with a tridiagonal solver are handled so.
!   This strongly improves stability in steep terrain. (MB)
!  For computation of lhfl_s use the correct lh_v or lh_s depending on 
!   interfacial temperature t_g (UB)
!  Implemented serialization comments (OF)
! V5_4         2016-03-10 Pascal Spoerri
!  Update of serialization statements
! V5_4a        2016-05-10 Matthias Raschendorfer
!  Use new module turb_data instead of data_turbulence
!  Use imode_turb now from turb_data
!  The lower boundary condition is now controlled by new switch lsflcnd
!      (and not imode_turb==0 any more: because meaning of imode_turb has changed)
!  The new switch (itype_vdif < 0) replaces (itype_turb/=3 .OR. imode_turb < 2)
!  The auxiliary fields za1t, za2t are no more used: 
!    they were only used for lvertdiff=.TRUE. and in that case were set to a1t, a2t
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

USE data_parameters , ONLY :   &
  wp,        & ! KIND-type parameters for real variables
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

  ! 4. constants for the horizontal rotated grid and related variables
  ! ------------------------------------------------------------------

  eddlon,       & ! 1 / dlon
  eddlat,       & ! 1 / dlat

  ! 5. variables for the time discretization and related variables
  ! --------------------------------------------------------------

  dt,           & ! long time-step

! 8. Organizational variables to handle the COSMO humidity tracers
! ----------------------------------------------------------------
  idt_qv,  idt_qc

! end of data_modelconfig

!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

  ! 2. physical constants and related variables
  ! -------------------------------------------

  t0_melt,      & ! absolute zero for temperature
  r_d,          & ! gas constant for dry air
  rvd_m_o,      & ! r_v/r_d - 1
  cp_d,         & ! specific heat of dry air at constant pressure
  rdocp,        & ! r_d / cp_d
  lh_v,         & ! latent heat of vapourization
  lh_s,         & ! latent heat of sublimation
  g,            & ! acceleration due to gravity (for COSMO_ART)
  r_earth,      & ! mean radius of the earth

  ! 3. constants for parametrizations
  ! ---------------------------------

  p0ref_recip     ! reciprocal of ref. pressure for Exner-fct. (Pa^-1)

! end of data_constants

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

  ! 1. constant fields for the reference atmosphere                 (unit)
  ! -----------------------------------------------
  p0         ,    & ! reference pressure at main levels             ( Pa  )
  hhl        ,    & ! geometical height of half model levels        (  m  )

  ! 2. external parameter fields                                    (unit)
  ! ----------------------------
  crlat      ,    & ! cosine of transformed latitude
  acrlat     ,    & ! 1 / ( crlat * radius of the earth )           ( 1/m )

  ! 3. prognostic variables                                         (unit)
  ! -----------------------
  u          ,    & ! zonal wind speed                              ( m/s )
  v          ,    & ! meridional wind speed                         ( m/s )
  w          ,    & ! vertical wind speed (defined on half levels)  ( m/s )
  t          ,    & ! temperature                                   (  k  )
  tke        ,    & ! turbulent kinetic energy (on half levels)     (m2/s2)
  pp         ,    & ! deviation from the reference pressure         ( pa  )

  ! 4. tendency fields for the prognostic variables                 (unit )
  ! -----------------------------------------------
  !    timely deviation  by diabatic and adiabatic processes 
  !    without sound-wave terms
  utens,          & ! u-tendency without sound-wave terms           ( m/s2)
  vtens,          & ! v-tendency without sound-wave terms           ( m/s2)
  wtens,          & ! w-tendency without sound-wave terms           ( m/s2)
  ttens,          & ! t-tendency without sound-wave terms           ( K/s )
  tketens,        & ! tke-tendency                                  (m2/s3)
  pptens            ! pp-tendency without sound-wave terms          (Pa/s )

USE data_fields     , ONLY :   &

  ! 5. fields for surface values and soil model variables           (unit )
  ! -----------------------------------------------------
  ps        ,     & ! surface pressure                              ( pa  )
! t_s       ,     & ! temperature of the ground surface             (  k  )
  t_g       ,     & ! weighted surface temperature                  (  k  )
  qv_s      ,     & ! specific water vapor content on the surface   (kg/kg)

  ! 6. fields that are computed in the parametrization and dynamics (unit )
  ! ---------------------------------------------------------------
  rho,          & ! density of moist air
  a1t, a2t  ,   & ! implicit weight of vertical diffusion

  !   turbulent coefficients in the atmosphere
  !   (defined on half levels)
  !   vertical turbulent diffusion coefficients
  tkvm     ,      & ! ... for momentum                               (m2/s)
  tkvh     ,      & ! ... for heat and moisture                      (m2/s)
  ! horizontal turbulent diffusion coefficients
  tkhm     ,      & ! ... for momentum                               (m2/s)
  tkhh     ,      & ! ... for heat and moisture                      (m2/s)

  !   turbulent coefficients at the surface
  tcm      ,    & ! turbulent diffusion coefficients for momentum   --
  tch      ,    & ! turbulent diffusion coefficients for heat       --
                  ! and moisture
  !   fields that are computed in the dynamics
  dqvdt    ,    & ! threedimensional moisture convergence         (  1/s)
  umfl_s   ,    & ! average u-momentum flux (surface)             ( N/m2)
  vmfl_s   ,    & ! average v-momentum flux (surface)             ( N/m2)
  shfl_s   ,    & ! average sensible heat flux (surface)          ( W/m2)
  qvsflx   ,    & ! surface flux of water vapour                  (kg/m2s)
  lhfl_s   ,    & ! average latent heat flux (surface)            ( W/m2)
  wcon     ,    & ! contravariant vertical velocity
  uadvt    ,    & ! advective tendency of u
  vadvt    ,    & ! advective tendency of v
  wadvt    ,    & ! advective tendency of w
  ppadvt   ,    & ! advective tendency of pp
  tadvt           ! advective tendency of t

! end of data_fields

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &
  nnow,           & ! corresponds to ntstep
  nnew,           & ! corresponds to ntstep + 1
  ntstep,         & ! current timestep (for serialization)
  ntke,           & ! current time index for tke in case of itype_turb=3
  ltur,           & ! forecast with turbulent diffusion
  itype_turb,     & ! type of turbulent diffusion parametrization
  itype_vdif,     & ! type of vertical diffusion calculation
  l3dturb,        & ! 3D-turbulence (additional horizontal diffusion)
  l3dturb_metr,   & ! switch on/off additional metric terms for 3D-turbulence
  lprog_tke,      & ! prognostic treatment of TKE (for itype_turb=5/7)
  lsppt,          & ! stochastic perturbed physics tendencies
  itype_qxpert_rn,& ! define which hum variables tend. are perturbed
  l_cosmo_art,    & ! if .TRUE., run the COSMO_ART
  l_pollen,       & ! of pollen
  idbg_level,     &
  lprintdeb_all,  &
  ldebug_dyn

!------------------------------------------------------------------------------

USE turb_data, ONLY :   &
  lsflcnd,        & ! lower surface flux condition
  imode_turb,     & ! mode of turbulent diffusion parametrization
  vel_min           ! minimal velocity scale [m/s]

!------------------------------------------------------------------------------

USE grid_metrics_utilities,   ONLY :  &
  sqrtg_r_s  ,    & ! reciprocal square root of G at skalar points  ( 1/m )
  sqrtg_r_u  ,    & ! reciprocal square root of G at u points       ( 1/m )
  sqrtg_r_v  ,    & ! reciprocal square root of G at v points       ( 1/m )
  sqrtg_r_w  ,    & ! reciprocal square root of G at w points       ( 1/m )
  dzeta_dlam ,    & ! metric coefficient                            ( 1   )
  dzeta_dphi ,    & ! metric coefficient                            ( 1   )
  dzeta_dlam_at_u  , & !
  dzeta_dlam_at_w  , & !
  dzeta_dlam_at_uw , & !
  dzeta_dphi_at_v  , & !
  dzeta_dphi_at_w  , & !
  dzeta_dphi_at_vw

!------------------------------------------------------------------------------

USE data_tracer         ,     ONLY :  &
  T_TURB_ID   , T_TURB_1D     , T_TURB_3D       , T_NAME_ID      ,            &
  T_BBC_ID    , T_BBC_ZEROVAL , T_BBC_ZEROFLUX  , T_BBC_SURF_VAL ,            &
  T_CLP_ID    , T_CLP_ON      , T_ERR_NOTFOUND

!------------------------------------------------------------------------------

USE data_parallel       ,     ONLY :  &
  my_cart_id

!------------------------------------------------------------------------------

USE environment         ,     ONLY :  &
  model_abort

!------------------------------------------------------------------------------

USE numeric_utilities_rk,     ONLY :  &
  clipping                    !

!------------------------------------------------------------------------------

USE src_tracer          ,     ONLY :  &
  trcr_get, trcr_meta_get, trcr_errorstr, trcr_get_ntrcr, trcr_calc, trcr_get_index

!------------------------------------------------------------------------------

USE src_stoch_physics,        ONLY :  &
  pertstoph

!==============================================================================

IMPLICIT NONE

!==============================================================================

!------------------------------------------------------------------------------
! Declarations
!------------------------------------------------------------------------------

! coefficients used in complete_tendencies
! ----------------------------------------
REAL    (KIND=wp)   ::  &
  za1t_surf, za2t_surf    !

REAL (KIND=wp),     ALLOCATABLE ::  &
  zsqrtgrho_r_s(:,:,:), & ! reciprocal square root of G * rho at skalar points
  zsqrtgrho_r_u(:,:,:), & ! reciprocal square root of G * rho at u points
  zsqrtgrho_r_v(:,:,:), & ! reciprocal square root of G * rho at v points
  zsqrtgrho_r_w(:,:,:), & ! reciprocal square root of G * rho at w points
  zpia(:,:,:),          & !
  zpianf(:,:,:),        & !
  ztheta(:,:,:),        & !
  ztheta_l(:,:,:),      & !
  ztmkvm(:,:,:),        & !
  ztmkvh(:,:),          & !
  ztch(:,:),            & !
  ztcm(:,:),            & !
  zkh(:,:,:),           & !
  ztmkvw(:,:,:),        & !
  ztmkvtke(:,:,:)         !

REAL (KIND=wp), ALLOCATABLE, DIMENSION(:,:,:), PRIVATE ::  &
  Kvm_rho_at_s,    &
  Kvm_rho_at_uw,   &
  Kvm_rho_at_vw,   &
  Kvh_rho_at_s,    &
  Kvh_rho_at_w,    &
  Kvh_rho_pi_at_w


LOGICAL                  ::  &
  lvertdiff,            & ! .TRUE. if vertical diffusion is calculated here
  lvertdiff_w,          & ! .TRUE. for vertical diffusion of w
  lmoist_turb,          & ! .TRUE. if moist turb. param. (Potsdam) is used
  lmassf_diffusion        !

! implicit weights:
REAL    (KIND=wp), ALLOCATABLE, PRIVATE ::  &
  beta_uu1(:),   &
  beta_uu2(:),   &
  beta_uu3(:),   &
  beta_vv1(:),   &
  beta_vv2(:),   &
  beta_vv3(:),   &
  beta_ww1(:),   &
  beta_ww2(:),   &
  beta_1(:),     &
  beta_2(:)

LOGICAL :: &
  l_diffusion_tracers  ! =.TRUE., if at least one tracer is diffused

INTEGER, PRIVATE :: izdebug

!==============================================================================

CONTAINS

!==============================================================================
!+ init procedure for completing the time step precalculate some coeffs.
!------------------------------------------------------------------------------

SUBROUTINE complete_tendencies_init

!------------------------------------------------------------------------------
!
! Description:
!   This procedure provides some coefficients needed in the computation of 
!   the final slow tendencies for the dynamic variables.
!
!------------------------------------------------------------------------------

! Declarations:

! Local scalars:
! -------------
INTEGER (KIND=iintegers) ::  &
  i,  j,  k              !  Loop indices in lon., lat. and vert. direction


REAL    (KIND=wp)  ::  &
  ztmkvtke_tmp(ie,je), & !
  ztvb, zvbke,         & !
  zdr, zpp,            & !
  zfpi, zppa,          & !
  zlhvocp,             & !
  ztkvz, ztkvl

INTEGER (KIND=iintegers) :: izerror
CHARACTER (LEN=255)      :: yzerrmsg
CHARACTER (LEN=50)       :: yzroutine


! Tracer pointers:
!----------------
REAL (KIND=wp),     POINTER :: &
  qc (:,:,:) => NULL()

! Statement function zfpi for calculation of the exner function
! where the dummy argument zppa is pressure

zfpi(zppa) = (zppa*p0ref_recip)**rdocp

! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine complete_tendencies_init
!------------------------------------------------------------------------------

  !$ser savepoint VerticalDiffusionUnittest.PrepareStep-in LargeTimeStep=ntstep
  !$ser data u_nnow=u(:,:,:,nnow) v_nnow=v(:,:,:,nnow)     &
  !$ser&     p_s_nnew=ps(:,:,nnew) p_s_nnow=ps(:,:,nnow) t_g=t_g(:,:,nnow) rho=rho &
  !$ser&     tkvm=tkvm tkvh=tkvh tch=tch tcm=tcm
  !$ser tracer QV#surf@nnow

!------------------------------------------------------------------------------
! Section 1: Some preparations
!------------------------------------------------------------------------------

  izerror   = 0_iintegers
  yzerrmsg  = ''
  yzroutine = 'complete_tendencies_init'

  IF (ldebug_dyn) THEN
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

  IF (izdebug >= 10) THEN
    WRITE(*,*) "Subr. [", TRIM(yzroutine), "]..."
  END IF

! Precalculate some help variables to avoid divisions in subsequent code
! zedsqrtgrho = 1/(sqrt(G)*rho)
!------------------------------------------------------------------------------

  DO  k = 1, ke
    DO j = jstart, jend
      DO i = istart, iend
        zsqrtgrho_r_s(i,j,k) = sqrtg_r_s(i,j,k) / rho(i,j,k)
      ENDDO
    ENDDO
  ENDDO

  DO  k = 1, ke
    DO j = jstartu, jendu
      DO i = istartu, iendu
        zsqrtgrho_r_u(i,j,k) = 2.0_wp*sqrtg_r_u(i,j,k)              &
          / ( rho(i,j,k)+rho(i+1,j,k) )
      ENDDO
    ENDDO
  ENDDO

  DO  k = 1, ke
    DO j = jstartv, jendv
      DO i = istartv, iendv
        zsqrtgrho_r_v(i,j,k) = 2.0_wp*sqrtg_r_v(i,j,k)              &
          / ( rho(i,j,k)+rho(i,j+1,k) )
      ENDDO
    ENDDO
  ENDDO

  IF ( lvertdiff_w ) THEN
    DO  k = 2, ke
      DO j = jstart, jend
        DO i = istart, iend
          zsqrtgrho_r_w(i,j,k) = 2.0_wp*sqrtg_r_w(i,j,k)            &
            / ( rho(i,j,k)+rho(i,j,k-1) )
        ENDDO
      ENDDO
    ENDDO
  END IF

  ! Setting of parameters for horizontal averaging of vertical diffusion
  ! coefficients: ztkvz, ztkvl
  ztkvz = 0.9_wp
  ztkvl = (1.0_wp-ztkvz)*0.25_wp

  ! In order to switch off the calculation of vertical diffusion,
  ! we use a local copy for the implicit weights a1t and a2t with values
  ! set to zero. This has to be done if the tendencies have already been
  ! calculated in turbdiff. An exception is the cloud ice, the tendencies
  ! have not been calculated in turbdiff.
  ! This will be removed when a final form of time-integration has been found.

!MR
  ! The calculation of vertical diffusion is now turned off by setting 
  ! lvertdiff to .FALSE. The local copies za1t, za2t for the implicit weights 
  ! a1t, a2t are no more needed therefore.

!MR  (itype_turb /= 3 .OR. imode_turb < 2) is now replaced by itype_vdif < 0
  IF ( ltur .AND. (itype_vdif < 0) ) THEN
    ! vertical diffusion needs to be calculated
    lvertdiff = .TRUE.
  ELSE
    ! no turbulence active or vertical diffusion already calculated
    lvertdiff = .FALSE.
  END IF

  ! Selection of lower boundary conditions

!MR (imode_turb == 0)  is now replaced by (.NOT. lscflnd)
  IF (.NOT.lsflcnd) THEN
    ! condition of the lower boundary is the surface concentration
    za1t_surf = a1t(ke1)  ! implicit weight for ke-value
    za2t_surf = a2t(ke1)  ! explicit weight for ke-value
  ELSE
    ! condition of the lower boundary is the explicit surface mass flux
    za1t_surf = 0.0_wp ! no implicit weight
    za2t_surf = 1.0_wp ! full explicit weight
  ENDIF

  ! Calculation of Exner-pressure and potential temperature
  ! -------------------------------------------------------
  DO  k = 1 , ke
    DO  j = 1, je
      DO  i = 1, ie
        zpp           = p0(i,j,k) + pp(i,j,k,nnow)
        zpia(i,j,k)   = zfpi( zpp )
        ztheta(i,j,k) = t(i,j,k,nnow) / zpia(i,j,k)
      ENDDO
    ENDDO
  ENDDO

  ! Some preparations for moist turbulence formulation of sensible heat flux
  ! ------------------------------------------------------------------------
  IF ( lmoist_turb .AND. l3dturb ) THEN
    
    zlhvocp = lh_v / cp_d

    ! Retrieve the required microphysics tracers
    CALL trcr_get(izerror, idt_qc, ptr_tlev = nnow, ptr = qc)
    IF (izerror /= 0) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF

    ! liquid water potential temperature
    DO  k = 1 , ke
      DO  j = 1, je
        DO  i = 1, ie
          ztheta_l(i,j,k) = ztheta(i,j,k)                               &
                          - zlhvocp * qc(i,j,k) / zpia(i,j,k)
        ENDDO
      ENDDO
    ENDDO
    
  END IF
  
  ! Calculation of modified transfer coefficients 
  ! ---------------------------------------------
  DO   j  =  jstart, jend+1
    DO i  =  istart, iend+1
      zvbke      = 0.5_wp*SQRT ( (u(i,j,ke,nnow) + u(i-1,j,ke,nnow))**2    &
                               + (v(i,j,ke,nnow) + v(i,j-1,ke,nnow))**2 )
! UB>> muss hier nicht t_g anstelle von t_s hin???
!      ztvb       = t_s (i,j,nnow)*(1.0_wp + rvd_m_o*qv_s(i,j,nnow))
      ztvb       = t_g (i,j,nnow)*(1.0_wp + rvd_m_o*qv_s(i,j,nnow))
! UB<<
      ztcm(i,j)  = tcm(i,j)*zvbke*ps(i,j,nnow)/(r_d*ztvb)
! UB>> limit the horizontal velocity zvbke to a max. of vel_min, to prevent 0 heat/moisture fluxes at 0 windspeed:
!  (NOT ACTIVATED FOR NOW) --> should also be done consistently in the rest
!  of the code, i.e., src_soil_multlay.f90
!!!      ztch(i,j)  = tch(i,j)*zvbke*ps(i,j,nnow)/(r_d*ztvb)
      ztch(i,j)  = tch(i,j)*MAX(zvbke,vel_min)*ps(i,j,nnow)/(r_d*ztvb)
! UB<<
    ENDDO
  ENDDO


  ! Calculation of the modified turbulent diffusion coefficients
  ! ------------------------------------------------------------

  DO  k = 2 , ke
    DO  j = 1, je
      DO  i = 1, ie
        zdr    = 0.5_wp*( rho(i,j,k) + rho(i,j,k-1) ) * sqrtg_r_w(i,j,k)
        ztmkvh(i,j  ) = tkvh(i,j,k) * zdr
        ztmkvm(i,j,k) = tkvm(i,j,k) * zdr
        zpianf(i,j,k) = 0.5_wp * ( zpia(i,j,k) + zpia(i,j,k-1) )
      ENDDO
    ENDDO
    DO    j = jstart, jend
      DO  i = istart, iend
        zkh(i,j,k) = ztkvz*ztmkvh(i,j)                       &
                   + ztkvl*( ztmkvh(i-1,j) + ztmkvh(i+1,j)   &
                           + ztmkvh(i,j-1) + ztmkvh(i,j+1) )
        ztmkvw  (i,j,k) = ztkvz*ztmkvm(i,j,k)                                 &
                        + ztkvl*( ztmkvm(i+1,j,k) + ztmkvm(i,j+1,k)           &
                                + ztmkvm(i-1,j,k) + ztmkvm(i,j-1,k) )
      ENDDO
    ENDDO
  ENDDO

  IF ( lprog_tke .AND. (itype_turb /= 3) ) THEN
    DO  k = 2 , ke
      DO  j = jstart-1, jend+1
        DO  i = istart-1, iend+1
          ztmkvtke_tmp(i,j) = tkvm(i,j,k) * sqrtg_r_w(i,j,k)
        ENDDO
      ENDDO
      DO    j = jstart, jend
        DO  i = istart, iend
          ztmkvtke(i,j,k) = ztkvz*ztmkvtke_tmp(i,j)                             &
                          + ztkvl*( ztmkvtke_tmp(i+1,j) + ztmkvtke_tmp(i,j+1)   &
                                  + ztmkvtke_tmp(i-1,j) + ztmkvtke_tmp(i,j-1) )
        ENDDO
      ENDDO
    ENDDO
  END IF
  
  DO    j = 1, je
    DO  i  = 1, ie
      zpianf(i,j,ke1) = zfpi( ps(i,j,nnow) )
    ENDDO
  ENDDO

!$ser savepoint VerticalDiffusionUnittest.PrepareStep-out LargeTimeStep=ntstep
!$ser data ztch=ztch ztcm=ztcm zkh=zkh ztmkvm=ztmkvm ztmkvw=ztmkvw zsqrtgrho_r_s=zsqrtgrho_r_s

END SUBROUTINE complete_tendencies_init

!==============================================================================
!==============================================================================
!+ procedure for completing the time step for u, v, w, t and pp 
!+ in combination with calculation of the diffusion at the beginning
!------------------------------------------------------------------------------

SUBROUTINE complete_tendencies_uvwtpp( nadv, irk, dtadv, opt_ark, opt_brk )

!------------------------------------------------------------------------------
!
! Description:
!   This procedure calculates the vertical tendencies for the 
!   prognostic variables u,v,w,pp and T.
!
! Method:
!   Using the previously calculated tendencies from horizontal advection,
!   which are stored on the tendency arrays, and tendencies due to 
!   physics and adiabatic processes of the dynamic variables, 
!   the vertical advection is solved by a vertically implicit 
!   scheme (modified Crank-Nicolson).
!
!------------------------------------------------------------------------------

! Declarations:

! Subroutine arguments:
! ---------------------

INTEGER (KIND=iintegers), INTENT(IN) ::  &
  nadv, irk

REAL (KIND=wp),     INTENT(IN) ::  &
  dtadv

REAL (KIND=wp),     INTENT(IN), OPTIONAL ::  &
  opt_ark, opt_brk


! Local scalars:
! -------------
INTEGER (KIND=iintegers) ::  &
  i,  j,  k              !  Loop indices in lon., lat. and vert. direction


REAL    (KIND=wp)  ::  &
  zbetav, zbetp, zbetm,& !
  zgav  , zgcv ,       & !
  zag, zas,            & !
  zcg, zcs,            & !
  zbg, zdg,            & !
  zd1g, zd2g,          & !
  znew, zz,            & !
  zzdtr,               & !
  zark, zbrk             !

! Local (automatic) arrays:
! ------------------------
REAL    (KIND=wp)    ::  &
  zc      (ie,je,ke),    & ! Upper main diagonal of ...
  zd1     (ie,je,ke),    & ! Right hand side of ...
  zd2     (ie,je,ke),    & ! Right hand side of ...
  ze1     (ie,je,ke),    & ! Soluton vector of ...
  ze2     (ie,je,ke)       ! Soluton vector of ...

! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine complete_tendencies_uvwtpp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Some preparations
!------------------------------------------------------------------------------

  ! test if optional parameters are present -
  ! used in combination with TVD-variant of Runge-Kutta scheme
  IF ( PRESENT(opt_ark) .AND. PRESENT(opt_brk) ) THEN
    zark = opt_ark
    zbrk = opt_brk
  ELSE
    zark = 1.0_wp
    zbrk = 0.0_wp
  END IF

  ! Setting of parameters for implicit calculation of vertical advection       
  ! (zbetav: weight for the n+1 time level           
  !          zbetav=0: centered, +1: full implicit,-1: explicit)
  zbetav = 0.0_wp
  zbetp  = 0.5_wp*( 1.0_wp + zbetav )
  zbetm  = 0.5_wp*( 1.0_wp - zbetav )

  ! Setting of reciprocal time step
  zzdtr = 1.0_wp/dtadv

! Precalculate some help variables to avoid divisions in subsequent code
!-----------------------------------------------------------------------

!NEC_CB Included below
!  DO  k = 1, ke
!    DO j = jstart, jend
!      DO i = istart, iend
!        zgavx (i,j,k) = - 0.5_wp*wcon(i,j,k  )
!        zgcvx (i,j,k) =   0.5_wp*wcon(i,j,k+1)
!      ENDDO
!    ENDDO
!  ENDDO

!------------------------------------------------------------------------------
! Section 2: Setup of tridiagonal matrix systems resulting from the implicit
!            numerical formulation of advection.
!            -  pressure perturbation and temperature
!------------------------------------------------------------------------------

  !$ser savepoint VerticalAdvectionUnittest.DoUVW-in LargeTimeStep=ntstep RKStageNumber=irk
  !$ser data dtadv=dtadv wcon=wcon
  !$ser data u_pos=zark*u(:,:,:,nnow)+zbrk*u(:,:,:,nadv) u=u(:,:,:,nadv) u_tens=utens uadvt=uadvt
  !$ser data v_pos=zark*v(:,:,:,nnow)+zbrk*v(:,:,:,nadv) v=v(:,:,:,nadv) v_tens=vtens vadvt=vadvt
  !$ser data w_pos=zark*w(:,:,:,nnow)+zbrk*w(:,:,:,nadv) w=w(:,:,:,nadv) w_tens=wtens wadvt=wadvt

  !$ser savepoint VerticalAdvectionUnittest.DoPPTP-in LargeTimeStep=ntstep RKStageNumber=irk
  !$ser data dtadv=dtadv wcon=wcon
  !$ser data pp_pos=zark*pp(:,:,:,nnow)+zbrk*pp(:,:,:,nadv) pp=pp(:,:,:,nadv) pp_tens=pptens ppadvt=ppadvt
  !$ser data t_pos=zark*t(:,:,:,nnow)+zbrk*t(:,:,:,nadv) t=t(:,:,:,nadv) t_tens=ttens tadvt=tadvt

  DO j = jstart, jend

    ! Top layer       
!!CDIR ON_ADB(pp)
!!CDIR ON_ADB(t)
!!CDIR ON_ADB(zc)
!CDIR ON_ADB(zd1)
!CDIR ON_ADB(zd2)
    DO i = istart, iend
      zgcv =   0.5_wp*wcon(i,j,2)
      zcg  = zgcv*zbetp
      zcs  = zgcv*zbetm
      zbg  = zzdtr - zcg
      zd1g = zzdtr * ( zark*pp(i,j,1,nnow) + zbrk*pp(i,j,1,nadv) )   &
           + pptens(i,j,1) + ppadvt(i,j,1)                           &
           - zcs * ( pp(i,j,2,nadv) - pp(i,j,1,nadv) )
      zd2g = zzdtr * ( zark*t (i,j,1,nnow) + zbrk*t (i,j,1,nadv) )   &
           + ttens (i,j,1) + tadvt (i,j,1)                           &
           - zcs * ( t (i,j,2,nadv) - t (i,j,1,nadv) )
      zd1(i,j,1) = zd1g / zbg
      zd2(i,j,1) = zd2g / zbg
      zc(i,j,1)  = zcg  / zbg
    ENDDO

    ! The layers from k=2 to k=ke-1
    DO k = 2, ke-1
!CDIR ON_ADB(pp)
!CDIR ON_ADB(t)
!CDIR ON_ADB(zc)
!CDIR ON_ADB(zd1)
!CDIR ON_ADB(zd2)
      DO i = istart, iend
        zgav = - 0.5_wp*wcon(i,j,k  )
        zgcv =   0.5_wp*wcon(i,j,k+1)
        zag  = zgav*zbetp
        zas  = zgav*zbetm
        zcg  = zgcv*zbetp
        zcs  = zgcv*zbetm
        zbg  = zzdtr - zag - zcg
        zd1g = zzdtr * ( zark*pp(i,j,k,nnow) + zbrk*pp(i,j,k,nadv) )   &
             + pptens(i,j,k) + ppadvt(i,j,k)                           &
             - zas * ( pp(i,j,k-1,nadv) - pp(i,j,k,nadv) )             &
             - zcs * ( pp(i,j,k+1,nadv) - pp(i,j,k,nadv) )
        zd2g = zzdtr * ( zark*t (i,j,k,nnow) + zbrk*t (i,j,k,nadv) )   &
             + ttens (i,j,k) + tadvt (i,j,k)                           &
             - zas * ( t (i,j,k-1,nadv) - t (i,j,k,nadv) )             &
             - zcs * ( t (i,j,k+1,nadv) - t (i,j,k,nadv) )
        zz   = 1.0_wp / ( zbg - zag*zc(i,j,k-1) )
        zc (i,j,k) = zcg * zz
        zd1(i,j,k) = ( zd1g -zag*zd1(i,j,k-1) ) * zz
        zd2(i,j,k) = ( zd2g -zag*zd2(i,j,k-1) ) * zz
      ENDDO
    ENDDO

    ! The bottom layer
!CDIR ON_ADB(pp)
!CDIR ON_ADB(t)
!CDIR ON_ADB(zc)
!CDIR ON_ADB(zd1)
!CDIR ON_ADB(zd2)
!CDIR ON_ADB(ze1)
!CDIR ON_ADB(ze2)
    DO i = istart, iend
      zgav = - 0.5_wp*wcon(i,j,ke)
      zag  = zgav*zbetp 
      zas  = zgav*zbetm 
      zbg  = zzdtr - zag
      zd1g = zzdtr * ( zark*pp(i,j,ke,nnow) + zbrk*pp(i,j,ke,nadv) )   &
           + pptens(i,j,ke) + ppadvt(i,j,ke)                           &
           - zas * ( pp(i,j,ke-1,nadv) - pp(i,j,ke,nadv) )
      zd2g = zzdtr * ( zark*t (i,j,ke,nnow) + zbrk*t (i,j,ke,nadv) )   &
           + ttens (i,j,ke) + tadvt (i,j,ke)                           &
           - zas * ( t (i,j,ke-1,nadv) - t (i,j,ke,nadv) )
      zz   = 1.0_wp / ( zbg - zag*zc(i,j,ke-1) )
      znew = ( zd1g - zag*zd1(i,j,ke-1) ) * zz
      ppadvt(i,j,ke) =  &
        ( znew - zark*pp(i,j,ke,nnow) - zbrk*pp(i,j,ke,nadv) ) * zzdtr
      ze1(i,j,ke) = znew
      znew = ( zd2g - zag*zd2(i,j,ke-1) ) * zz
      tadvt (i,j,ke) =  &
        ( znew - zark*t (i,j,ke,nnow) - zbrk*t (i,j,ke,nadv) ) * zzdtr
      ze2(i,j,ke) = znew
    ENDDO

    ! Backsubstitution and storage of the complete slow tendencies

    DO k = ke-1, 1, -1
!CDIR ON_ADB(ze1)
!CDIR ON_ADB(ze2)
      DO i = istart, iend
        ze1(i,j,k)     = zd1(i,j,k) - zc(i,j,k)*ze1(i,j,k+1)
        ppadvt(i,j,k) =  &
          ( ze1(i,j,k) - zark*pp(i,j,k,nnow) - zbrk*pp(i,j,k,nadv) ) * zzdtr
        ze2(i,j,k)     = zd2(i,j,k) - zc(i,j,k)*ze2(i,j,k+1)
        tadvt (i,j,k) =  &
          ( ze2(i,j,k) - zark*t (i,j,k,nnow) - zbrk*t (i,j,k,nadv) ) * zzdtr
      ENDDO
    ENDDO

  ENDDO ! j-loop

!------------------------------------------------------------------------------
! Section 3: Setup of tridiagonal matrix systems resulting from the implicit
!            numerical formulation of advection.
!            -  vertical velocity
!------------------------------------------------------------------------------

!NEC_CB Included below
!  DO  k = 2, ke
!    ! Precalculate some help variables to avoid divisions in subsequent code
!    DO j = jstart, jend
!      DO i = istart, iend
!        zgavx(i,j,k) =  &
!          - 0.25_wp*( wcon(i,j,k) + wcon(i,j,k-1) )
!        zgcvx(i,j,k) =  &
!          + 0.25_wp*( wcon(i,j,k) + wcon(i,j,k+1) )
!      ENDDO
!    ENDDO
!  ENDDO

  ! Top layer       
  DO j = jstart, jend
!CDIR ON_ADB(w)
!CDIR ON_ADB(zc)
!CDIR ON_ADB(zd1)
    DO i = istart, iend
      zgav       = - 0.25_wp*( wcon(i,j,2) + wcon(i,j,1) )
      zgcv       = + 0.25_wp*( wcon(i,j,2) + wcon(i,j,3) )
      zcg        = zgcv*zbetp
      zbg        = zzdtr - zcg - zgav*zbetp
      zdg        = zzdtr * ( zark*w(i,j,2,nnow) + zbrk*w(i,j,2,nadv) )   &
        + wtens(i,j,2) + wadvt(i,j,2)                                    &
        - zgcv*zbetm * ( w(i,j,3,nadv) - w(i,j,2,nadv) )         &
        + zgav*zbetm*w(i,j,2,nadv) - zgav*w(i,j,1,nadv)
      zd1(i,j,2) = zdg / zbg
      zc (i,j,2) = zcg / zbg
    ENDDO

    ! The layers from k=3 to k=ke-1
    DO k = 3, ke-1
!CDIR ON_ADB(w)
!CDIR ON_ADB(zc)
!CDIR ON_ADB(zd1)
      DO i = istart, iend
        zgav       = - 0.25_wp*( wcon(i,j,k) + wcon(i,j,k-1) )
        zgcv       = + 0.25_wp*( wcon(i,j,k) + wcon(i,j,k+1) )
        zag        = zgav*zbetp
        zas        = zgav*zbetm
        zcg        = zgcv*zbetp
        zcs        = zgcv*zbetm
        zbg        = zzdtr - zag - zcg
        zdg        = zzdtr * ( zark*w(i,j,k,nnow) + zbrk*w(i,j,k,nadv) )  &
          + wtens(i,j,k) + wadvt(i,j,k)                                   &
          - zas * ( w(i,j,k-1,nadv) - w(i,j,k,nadv) )                     &
          - zcs * ( w(i,j,k+1,nadv) - w(i,j,k,nadv) )
        zz         = 1.0_wp/( zbg - zag*zc(i,j,k-1) )
        zc (i,j,k) = zcg * zz
        zd1(i,j,k) = ( zdg - zag*zd1(i,j,k-1) ) * zz
      ENDDO
    ENDDO

    ! The bottom layer
!CDIR ON_ADB(ze1)
    DO i = istart, iend
      zgav       = - 0.25_wp*( wcon(i,j,ke) + wcon(i,j,ke-1) )
      zgcv       = + 0.25_wp*( wcon(i,j,ke) + wcon(i,j,ke+1) )
      zag        = zgav*zbetp 
      zbg        = zzdtr - zag - zgcv*zbetp
      zdg        = zzdtr * ( zark*w(i,j,ke,nnow) + zbrk*w(i,j,ke,nadv) )  &
        + wtens(i,j,ke) + wadvt(i,j,ke)                                   &
        - zgav*zbetm*( w(i,j,ke-1,nadv) - w(i,j,ke,nadv) )       &
        + zgcv*zbetm*w(i,j,ke,nadv)                              &
        - zgcv*w(i,j,ke1,nadv)
      znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
      wadvt(i,j,ke) =  &
        ( znew - zark*w(i,j,ke,nnow) - zbrk*w(i,j,ke,nadv) ) * zzdtr
      ze1(i,j,ke) = znew
    ENDDO

    ! Backsubstitution and storage of the complete slow tendencies

    DO k = ke-1, 2, -1
!CDIR ON_ADB(ze1)
      DO i = istart, iend
        ze1(i,j,k)     = zd1(i,j,k) - zc(i,j,k)*ze1(i,j,k+1)
        wadvt (i,j,k) =  &
          ( ze1(i,j,k) - zark*w(i,j,k,nnow) - zbrk*w(i,j,k,nadv) ) * zzdtr
      ENDDO
    ENDDO
  ENDDO ! j-loop

!------------------------------------------------------------------------------
! Section 4: Setup of tridiagonal matrix systems resulting from the implicit
!            numerical formulation of advection.
!            -  horizontal wind velocity u
!------------------------------------------------------------------------------

  ! Top layer:   k = 1
  DO j = jstartu, jendu
!CDIR ON_ADB(wcon)
!CDIR ON_ADB(u)
!CDIR ON_ADB(zc)
!CDIR ON_ADB(zd1)
    DO i = istartu, iendu
      zgcv       = 0.25_wp*( wcon(i,j,2)+wcon(i+1,j,2) )
      zcg        = zgcv*zbetp
      zbg        = zzdtr - zcg
      zdg        = zzdtr * ( zark*u(i,j,1,nnow) + zbrk*u(i,j,1,nadv) )  &
        + utens(i,j,1) + uadvt(i,j,1)                                   &
        - zgcv*zbetm * ( u(i,j,2,nadv) - u(i,j,1,nadv) )
      zd1(i,j,1) = zdg / zbg
      zc (i,j,1) = zcg / zbg
    ENDDO

    ! The layers from k=2 to k=ke-1
    DO k = 2, ke-1
!CDIR ON_ADB(wcon)
!CDIR ON_ADB(u)
!CDIR ON_ADB(zc)
!CDIR ON_ADB(zd1)
      DO i = istartu, iendu
        zgav       =  &
          -0.25_wp*( wcon(i,j,k  )+wcon(i+1,j,k  ) )
        zgcv       =  &
           0.25_wp*( wcon(i,j,k+1)+wcon(i+1,j,k+1) )
        zag        = zgav*zbetp
        zas        = zgav*zbetm
        zcg        = zgcv*zbetp
        zcs        = zgcv*zbetm
        zbg        = zzdtr - zag - zcg
        zdg        = zzdtr * ( zark*u(i,j,k,nnow) + zbrk*u(i,j,k,nadv) )  &
          + utens(i,j,k) + uadvt(i,j,k)                                   &
          - zas * ( u(i,j,k-1,nadv) - u(i,j,k,nadv) )                     &
          - zcs * ( u(i,j,k+1,nadv) - u(i,j,k,nadv) )
        zz         = 1.0_wp / ( zbg - zag*zc(i,j,k-1) )
        zc (i,j,k) = zcg * zz
        zd1(i,j,k) = ( zdg -zag*zd1(i,j,k-1) ) * zz
      ENDDO
    ENDDO

    ! The bottom layer:  k = ke
!CDIR ON_ADB(ze1)
    DO i = istartu, iendu
      zgav       = - 0.25_wp*( wcon(i,j,ke)+wcon(i+1,j,ke) )
      zag        = zgav*zbetp
      zas        = zgav*zbetm
      zbg        = zzdtr - zag
      zdg        = zzdtr * ( zark*u(i,j,ke,nnow) + zbrk*u(i,j,ke,nadv) )  &
        + utens(i,j,ke) + uadvt(i,j,ke)                                   &
        - zas * ( u(i,j,ke-1,nadv) - u(i,j,ke,nadv) )
      znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
      uadvt (i,j,ke) =  &
        ( znew - zark*u(i,j,ke,nnow) - zbrk*u(i,j,ke,nadv) ) * zzdtr
      ze1(i,j,ke) = znew
    ENDDO

    ! Backsubstitution and storage of the complete slow tendencies

    DO k = ke-1, 1, -1
!CDIR ON_ADB(ze1)
      DO i = istartu, iendu
        ze1(i,j,k) =  zd1(i,j,k) - zc(i,j,k)*ze1(i,j,k+1)
        uadvt(i,j,k) =  &
          ( ze1(i,j,k) - zark*u(i,j,k,nnow) - zbrk*u(i,j,k,nadv) ) * zzdtr
      ENDDO
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section 5: Setup of tridiagonal matrix systems resulting from the implicit
!            numerical formulation of advection.
!            -  horizontal wind velocity v
!------------------------------------------------------------------------------

  ! Top layer  k=1
  DO j = jstartv, jendv
!CDIR ON_ADB(wcon)
!CDIR ON_ADB(v)
!CDIR ON_ADB(zc)
!CDIR ON_ADB(zd1)
    DO i = istartv, iendv
      zgcv       = 0.25_wp*( wcon(i,j,2)+wcon(i,j+1,2) )
      zcg        = zgcv*zbetp
      zbg        = zzdtr - zcg
      zdg        = zzdtr * ( zark*v(i,j,1,nnow) + zbrk*v(i,j,1,nadv) )  &
        + vtens(i,j,1) + vadvt(i,j,1)                                   &
        - zgcv*zbetm * ( v(i,j,2,nadv) - v(i,j,1,nadv) )
      zd1(i,j,1) = zdg / zbg
      zc (i,j,1) = zcg / zbg
    ENDDO

    ! The layers from k=2 to k=ke-1
    DO  k = 2, ke-1
!CDIR ON_ADB(wcon)
!CDIR ON_ADB(v)
!CDIR ON_ADB(zc)
!CDIR ON_ADB(zd1)
      DO i = istartv, iendv
        zgav       = -0.25_wp*( wcon(i,j,k  )+wcon(i,j+1,k  ) )
        zgcv       =  0.25_wp*( wcon(i,j,k+1)+wcon(i,j+1,k+1) )
        zag        = zgav*zbetp
        zas        = zgav*zbetm
        zcg        = zgcv*zbetp
        zcs        = zgcv*zbetm 
        zbg        = zzdtr - zag - zcg
        zdg        = zzdtr * ( zark*v(i,j,k,nnow) + zbrk*v(i,j,k,nadv) )  &
          + vtens(i,j,k) + vadvt(i,j,k)                                   &
          - zas * ( v(i,j,k-1,nadv) - v(i,j,k,nadv) )                     &
          - zcs * ( v(i,j,k+1,nadv) - v(i,j,k,nadv) )
        zz         = 1.0_wp / ( zbg - zag*zc(i,j,k-1) )
        zc (i,j,k) = zcg * zz
        zd1(i,j,k) = ( zdg -zag*zd1(i,j,k-1) ) * zz
      ENDDO
    ENDDO

    ! The bottom layer k=ke
!CDIR ON_ADB(ze1)
    DO i = istartv, iendv
      zgav       = - 0.25_wp*( wcon(i,j,ke)+wcon(i,j+1,ke) )
      zag        = zgav*zbetp
      zas        = zgav*zbetm
      zbg        = zzdtr - zag
      zdg        = zzdtr * ( zark*v(i,j,ke,nnow) + zbrk*v(i,j,ke,nadv) )  &
        + vtens(i,j,ke) + vadvt(i,j,ke)                                   &
        - zas * ( v(i,j,ke-1,nadv) - v(i,j,ke,nadv) )
      znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
      vadvt (i,j,ke) =  &
        ( znew - zark*v(i,j,ke,nnow) - zbrk*v(i,j,ke,nadv) ) * zzdtr
      ze1(i,j,ke) = znew
    ENDDO

    ! Backsubstitution and storage of the complete slow tendencies

    DO k = ke-1, 1, -1
!CDIR ON_ADB(ze1)
      DO i = istartv, iendv
        ze1(i,j,k) = zd1(i,j,k) - zc(i,j,k)*ze1(i,j,k+1)
        vadvt(i,j,k) =  &
          ( ze1(i,j,k) - zark*v(i,j,k,nnow) - zbrk*v(i,j,k,nadv) ) * zzdtr
      ENDDO
    ENDDO
  ENDDO

  !$ser savepoint VerticalAdvectionUnittest.DoUVW-out LargeTimeStep=ntstep RKStageNumber=irk
  !$ser data uadvt=uadvt vadvt=vadvt wadvt=wadvt

  !$ser savepoint VerticalAdvectionUnittest.DoPPTP-out LargeTimeStep=ntstep RKStageNumber=irk
  !$ser data ppadvt=ppadvt tadvt=tadvt

!------------------------------------------------------------------------------
! End of subroutine complete_tendencies_uvwtpp
!------------------------------------------------------------------------------

END SUBROUTINE complete_tendencies_uvwtpp

!==============================================================================
!==============================================================================
!+ procedure for completing the time step for u, v, w, t and pp
!+ in combination with vertical explicit advection and
!+ calculation of the vertical diffusion at the beginning
!------------------------------------------------------------------------------

!option! -pvctl _on_adb
SUBROUTINE complete_tendencies_uvwtpp_eva

!------------------------------------------------------------------------------
!
! Description:
!   This procedure provides the final slow tendencies for the dynamic 
!   variables.
!
!------------------------------------------------------------------------------

! Declarations:

! Local scalars:
! -------------
INTEGER (KIND=iintegers) ::  &
  i,  j,  k              !  Loop indices in lon., lat. and vert. direction

! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine complete_tendencies_uvwtpp_eva
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section : Pressure perturbation
!------------------------------------------------------------------------------

  DO k = 1, ke
    DO j = jstart, jend
      DO i = istart, iend
        ppadvt(i,j,k) = ppadvt(i,j,k) + pptens(i,j,k)
      ENDDO
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section : Vertical velocity
!------------------------------------------------------------------------------

  DO k = 2, ke
    DO j = jstart, jend
      DO i = istart, iend
        wadvt(i,j,k) = wadvt(i,j,k) + wtens(i,j,k)
      ENDDO
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section : Temperature 
!------------------------------------------------------------------------------

  DO k = 1, ke
    DO j = jstart, jend
      DO i = istart, iend
        tadvt(i,j,k) = tadvt(i,j,k) + ttens(i,j,k)
      ENDDO
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section : Horizontal wind velocity u 
!------------------------------------------------------------------------------

  DO k = 1, ke
    DO j = jstartu, jendu
      DO i = istartu, iendu
        uadvt(i,j,k) = uadvt(i,j,k) + utens(i,j,k)
      ENDDO
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section : Horizontal wind velocity v
!------------------------------------------------------------------------------

  DO k = 1, ke
    DO j = jstartv, jendv
      DO i = istartv, iendv
        vadvt(i,j,k) = vadvt(i,j,k) + vtens(i,j,k)
      ENDDO
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! End of subroutine complete_tendencies_uvwtpp_eva
!------------------------------------------------------------------------------

END SUBROUTINE complete_tendencies_uvwtpp_eva

!==============================================================================
!==============================================================================
!+ procedure for computing the moisture divergence
!------------------------------------------------------------------------------

SUBROUTINE compute_moisture_divergence

!------------------------------------------------------------------------------
!
! Description:
!   This procedure is only used in src_runge_kutta and calculates the
!   moisture divergence (used for the convection scheme
!
! Method:
!   Using the previously calculated updates (on the timelevel nnew) for the
!   moisture variables from advection, physics and adiabatic processes, the 
!   vertical diffusion is solved by a vertically implicit scheme (modified 
!   Crank-Nicolson).
!
!------------------------------------------------------------------------------

! Declarations:

! Local scalars:
! -------------
INTEGER (KIND=iintegers) ::  &
  i,  j,  k              !  Loop indices in lon., lat. and vert. direction

INTEGER (KIND=iintegers) :: izerror

CHARACTER(LEN=255)       :: yzerrmsg
CHARACTER(LEN=50)        :: yzroutine

REAL    (KIND=wp)    ::  &
  zgat, zgct,            & !
  zag, zas,              & !
  zcg, zcs,              & !
  zbg, zdg,              & !
  zz, zzdtr                !

! Local (automatic) arrays:
! ------------------------
REAL    (KIND=wp)    ::  &
  zc      (ie,je,ke),    & ! Upper main diagonal of ...
  zd      (ie,je,ke),    & ! Right hand side of ...
  ze      (ie,je,ke)       ! Soluton vector of ...

! Tracer pointers
!----------------
REAL (KIND=wp),     POINTER :: &
  qv_now  (:,:,:) => NULL(),        & ! QV at nnow
  qv_new  (:,:,:) => NULL()           ! QV at nnew

! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine compute_moisture_divergence
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Some preparations
!------------------------------------------------------------------------------

  izerror   = 0_iintegers
  yzerrmsg  = ''
  yzroutine = 'compute_moisture_div'

  ! Setting of reciprocal time step
  zzdtr = 1.0_wp / dt

  CALL trcr_calc(1, idt_qv)

  ! Retrieve the required microphysics tracers (at specified timelevel)
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnow, ptr = qv_now)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnew, ptr = qv_new)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  IF ( lvertdiff ) THEN
    !US means:  (ltur .AND. (itype_vdif < 0) )

    !--------------------------------------------------------------------------
    ! Section 2a: Setup of tridiagonal matrix systems resulting from the
    !             implicit numerical formulation of diffusion of qv
    !--------------------------------------------------------------------------

    ! First, the matrix elements a(k), b(k) and c(k) of the coefficient
    ! matrix are set (these are the same for qv, qc and qi).
    ! The right hand side is stored on d(k).

    ! Top layer
    DO j = jstart, jend
!CDIR ON_ADB(zc)
!CDIR ON_ADB(zd)
      DO i = istart, iend
        zgct       = - zkh(i,j,2)*zsqrtgrho_r_s(i,j,1)
        zcg        = zgct*a1t(2)
        zcs        = zgct*a2t(2)
        zbg        = zzdtr - zcg
        zdg        = zzdtr * qv_new(i,j,1)                             &
                   - zcs * ( qv_now(i,j,2) - qv_now(i,j,1) )
        zc (i,j,1) = zcg / zbg
        zd (i,j,1) = zdg / zbg
      ENDDO
    ENDDO

    ! The layers from k=2 to k=ke-1
    DO k = 2, ke-1
      DO j = jstart, jend
!CDIR ON_ADB(zc)
!CDIR ON_ADB(zd)
        DO i = istart, iend
          zgat       = - zkh(i,j,k  )*zsqrtgrho_r_s(i,j,k)
          zgct       = - zkh(i,j,k+1)*zsqrtgrho_r_s(i,j,k)
          zag        = zgat*a1t(k)
          zas        = zgat*a2t(k)
          zcg        = zgct*a1t(k+1)
          zcs        = zgct*a2t(k+1)
          zbg        = zzdtr - zag - zcg
          zdg        = zzdtr * qv_new(i,j,k)                           &
                     - zas * ( qv_now(i,j,k-1) - qv_now(i,j,k) )       &
                     - zcs * ( qv_now(i,j,k+1) - qv_now(i,j,k) )
          zz         = 1.0_wp / ( zbg - zag*zc(i,j,k-1) )
          zc (i,j,k) = zcg * zz
          zd (i,j,k) = ( zdg -zag*zd(i,j,k-1) ) * zz
        ENDDO
      ENDDO
    ENDDO

    ! the bottom layer
    DO j = jstart, jend
!CDIR ON_ADB(zc)
!CDIR ON_ADB(zd)
!CDIR ON_ADB(ze)
      DO i = istart, iend
        zgat       = - zkh (i,j,ke)*zsqrtgrho_r_s(i,j,ke)
        zgct       = - ztch(i,j   )*zsqrtgrho_r_s(i,j,ke)
        zag        = zgat*a1t(ke)
        zas        = zgat*a2t(ke)
        zcg        = zgct*za1t_surf
        zcs        = zgct*za2t_surf
        zbg        = zzdtr - zag - zcg
        zdg        = zzdtr * qv_new(i,j,ke  )                          &
                   - zas * ( qv_now(i,j,ke-1) - qv_now(i,j,ke) )       &
                   + zcs *   qv_now(i,j,ke  ) - zgct * qv_s(i,j,nnow)
        zz         = 1.0_wp / ( zbg  - zag*zc(i,j,ke-1) )
        ze(i,j,ke) = ( zdg - zag*zd(i,j,ke-1) ) * zz
        dqvdt(i,j,ke) =  &
          ( MAX(0.0_wp, ze(i,j,ke)) - qv_now(i,j,ke) )*zzdtr
      ENDDO
    ENDDO

    ! Backsubstitution and storage of the complete slow tendencies (convergence):
    DO k = ke-1, 1, -1
      DO j = jstart, jend
!CDIR ON_ADB(ze)
        DO i = istart, iend
          ze(i,j,k) = zd(i,j,k) - zc(i,j,k)*ze(i,j,k+1)
          dqvdt(i,j,k) =  &
            ( MAX(0.0_wp,ze(i,j,k)) - qv_now(i,j,k) ) * zzdtr
        ENDDO
      ENDDO
    ENDDO

!MR: new case
  ELSEIF (ltur .AND. (itype_vdif > 0) ) THEN
    ! vertical qv-diffusion is already contained in 'dqvdt' and needs to be
    ! incremented by the part of qv-flux-convergence, which has already been 
    ! added to 'qv_new': 

    DO k = 1, ke
      DO j = jstart, jend
        DO i = istart, iend
          dqvdt(i,j,k) =  dqvdt(i,j,k) + &
            ( MAX(0.0_wp, qv_new(i,j,k) ) - qv_now(i,j,k) )*zzdtr
        ENDDO
      ENDDO
    ENDDO

  ELSE  !US means:  (.NOT. ltur  .OR. (itype_vdif == 0))
    ! qv-flux-convergence needs to be with the part, which already has been 
    ! added to 'qv_new' and may later be completed with vertical qv-diffusion 
    ! by the call of the turbulence-scheme:

    DO k = 1, ke
      DO j = jstart, jend
        DO i = istart, iend
          dqvdt(i,j,k) =  &
            ( MAX(0.0_wp, qv_new(i,j,k) ) - qv_now(i,j,k) )*zzdtr
        ENDDO
      ENDDO
    ENDDO
    
  END IF

  CALL trcr_calc(-1, idt_qv)

  !----------------------------------------------------------------------------
  ! End of subroutine compute_moisture_divergence 
  !----------------------------------------------------------------------------

END SUBROUTINE compute_moisture_divergence

!==============================================================================
!==============================================================================
!+ procedure for completing the time step for tke
!------------------------------------------------------------------------------

SUBROUTINE complete_tendencies_tke

!------------------------------------------------------------------------------
!
! Description:
!   This procedure is only used in src_runge_kutta and calculates the
!   final updates for the prognostic tke at time level n+1 (nnew).
!   The vertical diffusion as the last slow tendency is computed here for tke.
!
! Method:
!   Using the previously calculated updates (on the timelevel nnew) for the
!   tke variable from advection, physics and adiabatic processes, the 
!   vertical diffusion is solved by a vertically implicit scheme (modified 
!   Crank-Nicolson). The tke variable is updated directly.
!
!------------------------------------------------------------------------------

! Declarations:

! Local scalars:
! -------------
INTEGER (KIND=iintegers) ::  &
  i,  j,  k              !  Loop indices in lon., lat. and vert. direction

INTEGER (KIND=iintegers) :: &
  km1, kp1
  
REAL    (KIND=wp)    ::  &
  zag, zas, z2ag,        & !
  zcg, zcs, z2cg ,       & !
  zbg, zdg,              & !
  zz, zzdtr                !

! Local (automatic) arrays:
! ------------------------
REAL    (KIND=wp)    ::  &
  zgatz   (ie,je,ke),    & !
  zgctz   (ie,je,ke),    & !
  zc      (ie,je,ke),    & ! Upper main diagonal of ...
  zd1     (ie,je,ke)       ! Right hand side of ...

REAL    (KIND=wp)    ::  &
  rho_inv_at_w,   &
  r_earth_recip

CHARACTER(LEN=50) :: yzroutine

! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine complete_tendencies_tke  
!------------------------------------------------------------------------------

  yzroutine = "complete_tendencies_tke"

  IF ( izdebug >= 10 ) THEN
    WRITE(*,*) "subr. [", TRIM(yzroutine), "] ..."
  END IF

!------------------------------------------------------------------------------
! Section 1: Some preparations
!------------------------------------------------------------------------------

  ! Setting of reciprocal time step
  zzdtr = 1.0_wp / dt

  r_earth_recip = 1.0_wp / r_earth

  !----------------------------------------------------------------------------
  ! Section 2: Do vertical diffusion for turbulent kinetic energy (tke),
  !            if required.
  !----------------------------------------------------------------------------

  IF ( lprog_tke .AND. (itype_turb /= 3) ) THEN
    
    DO  k = 2, ke
      
      km1 = MAX( 2, k-1 )
      kp1 = MIN( ke, k+1 )
      
      ! Precalculate some help variables to avoid divisions in subsequent code
      DO j = jstart, jend
        DO i = istart, iend
          zgatz(i,j,k) = - sqrtg_r_w(i,j,k)   &
                         * ( ztmkvtke(i,j,k)+ztmkvtke(i,j,km1) )
          zgctz(i,j,k) = - sqrtg_r_w(i,j,k)   &
                         * ( ztmkvtke(i,j,kp1)+ztmkvtke(i,j,k) )
        ENDDO
      ENDDO
      
    ENDDO

    ! Top layer       
    DO j = jstart, jend
      DO i = istart, iend
        zag        = 0.5_wp*(a1t(2)+a1t(1)) * zgatz(i,j,2)
        zas        = 0.5_wp*(a2t(2)+a2t(1)) * zgatz(i,j,2)
        zcg        = 0.5_wp*(a1t(3)+a1t(2)) * zgctz(i,j,2)
        zcs        = 0.5_wp*(a2t(3)+a2t(2)) * zgctz(i,j,2)
        zbg        = zzdtr - zag - zcg
        zdg        = zzdtr * tke(i,j,2,nnew) + tketens(i,j,2)                &
                     + zas * tke(i,j,2,nnow)                                 &
                     - zcs * ( tke(i,j,3,nnow) - tke(i,j,2,nnow) )
        zd1(i,j,2) = zdg/zbg
        zc (i,j,2) = zcg/zbg
      ENDDO
    ENDDO

    IF ( l3dturb_metr ) THEN

      ! The layers from k=3 to k=ke-1
      DO k = 3, ke-1
        DO j = jstart, jend
          DO i = istart, iend

            rho_inv_at_w = 2.0_wp / ( rho(i,j,k) + rho(i,j,k-1) )

            ! the explicit counterparts of z2ag, z2cg are in 'explicit_horizontal_diffusion':
            z2ag = rho_inv_at_w * Kvh_rho_at_s(i,j,k) * beta_1(k)        &
              &  * (   acrlat(j,1)   * dzeta_dlam_at_w(i,j,k)            &
              &      * acrlat(j,1)   * dzeta_dlam(i,j,k)                 &
              &    +   r_earth_recip * dzeta_dphi_at_w(i,j,k)            &
              &      * r_earth_recip * dzeta_dphi(i,j,k) )

            z2cg = rho_inv_at_w * Kvh_rho_at_s(i,j,k+1) * beta_1(k+1)    &
              &  * (   acrlat(j,1)   * dzeta_dlam_at_w(i,j,k)            &
              &      * acrlat(j,1)   * dzeta_dlam(i,j,k+1)               &
              &    +   r_earth_recip * dzeta_dphi_at_w(i,j,k)            &
              &      * r_earth_recip * dzeta_dphi(i,j,k+1) )

            zag        = 0.5_wp*(a1t(k)+a1t(k-1)) * zgatz(i,j,k) - z2ag
            zas        = 0.5_wp*(a2t(k)+a2t(k-1)) * zgatz(i,j,k)
            zcg        = 0.5_wp*(a1t(k+1)+a1t(k)) * zgctz(i,j,k) - z2cg
            zcs        = 0.5_wp*(a2t(k+1)+a2t(k)) * zgctz(i,j,k)
            zbg        = zzdtr - zag - zcg
            zdg        = zzdtr * tke(i,j,k,nnew) + tketens(i,j,k)              &
                       - zas * ( tke(i,j,k-1,nnow) - tke(i,j,k,nnow) )         &
                       - zcs * ( tke(i,j,k+1,nnow) - tke(i,j,k,nnow) )
            zz         = 1.0_wp/( zbg - zag*zc(i,j,k-1) )
            zc (i,j,k) = zcg*zz
            zd1(i,j,k) = ( zdg - zag*zd1(i,j,k-1) )*zz
          ENDDO
        ENDDO
      ENDDO

    ELSE

      ! The layers from k=3 to k=ke-1
      DO k = 3, ke-1
        DO j = jstart, jend
          DO i = istart, iend
            zag        = 0.5_wp*(a1t(k)+a1t(k-1)) * zgatz(i,j,k)
            zas        = 0.5_wp*(a2t(k)+a2t(k-1)) * zgatz(i,j,k)
            zcg        = 0.5_wp*(a1t(k+1)+a1t(k)) * zgctz(i,j,k)
            zcs        = 0.5_wp*(a2t(k+1)+a2t(k)) * zgctz(i,j,k)
            zbg        = zzdtr - zag - zcg
            zdg        = zzdtr * tke(i,j,k,nnew) + tketens(i,j,k)              &
                       - zas * ( tke(i,j,k-1,nnow) - tke(i,j,k,nnow) )         &
                       - zcs * ( tke(i,j,k+1,nnow) - tke(i,j,k,nnow) )
            zz         = 1.0_wp/( zbg - zag*zc(i,j,k-1) )
            zc (i,j,k) = zcg*zz
            zd1(i,j,k) = ( zdg - zag*zd1(i,j,k-1) )*zz
          ENDDO
        ENDDO
      ENDDO

    END IF

    ! The bottom layer: We assume TKE=0 at the lowest level!
    DO j = jstart, jend
      DO i = istart, iend
        zag        = 0.5_wp*(a1t(ke)+a1t(ke-1)) * zgatz(i,j,ke)
        zas        = 0.5_wp*(a2t(ke)+a2t(ke-1)) * zgatz(i,j,ke)
        zcg        = 0.5_wp*(a1t(ke1)+a1t(ke))  * zgctz(i,j,ke)
        zcs        = 0.5_wp*(a2t(ke1)+a2t(ke))  * zgctz(i,j,ke)
        zbg        = zzdtr - zag - zcg
        zdg        = zzdtr * tke(i,j,ke,nnew) + tketens(i,j,ke)               &
                     - zas * ( tke(i,j,ke-1,nnow) - tke(i,j,ke,nnow) )        &
                     + zcs * tke(i,j,ke,nnow)
        zz         = 1.0_wp / ( zbg  - zag*zc(i,j,ke-1) )
        tke(i,j,ke,nnew) = ( zdg -zag*zd1(i,j,ke-1) ) * zz
      ENDDO
    ENDDO

    ! Backsubstitution
    DO k = ke-1, 2, -1
      DO j = jstart, jend
        DO i = istart, iend
          tke(i,j,k,nnew) = zd1(i,j,k) - zc(i,j,k)*tke(i,j,k+1,nnew)
        ENDDO
      ENDDO
    ENDDO

    CALL clipping( tke(:,:,1:ke,nnew), ie, je, ke )

  END IF

  !----------------------------------------------------------------------------
  ! End of subroutine complete_tendencies_tke 
  !----------------------------------------------------------------------------

END SUBROUTINE complete_tendencies_tke

!==============================================================================
!==============================================================================
!+ procedure for completing the time step for the tracer species
!------------------------------------------------------------------------------

SUBROUTINE complete_tendencies_trcr

!------------------------------------------------------------------------------
!
! Description:
!   This procedure is only used in src_runge_kutta and calculates the
!   final updates for the prognostic tracer variables 
!   at time level n+1 (nnew). 
!   The vertical diffusion as the last slow tendency 
!   is computed here for all tracer variables the vertical diffusion acts on.
!
! Method:
!   Using the previously calculated updates (on the timelevel nnew) for the
!   tracer variables from advection, physics and adiabatic processes, the 
!   vertical diffusion is solved by a vertically implicit scheme (modified 
!   Crank-Nicolson). The tracer variables
!   are updated directly. Also, the surface fluxes of tracers are stored.
!
!------------------------------------------------------------------------------

! Declarations:

! Local scalars:
! -------------
INTEGER (KIND=iintegers) ::  &
  i,  j,  k, iztrcr            ! Loop indices in lon., lat. and vert. direction 
                               ! and for tracers

INTEGER (KIND=iintegers) ::  & 
  izerror,izqv, izqc, izqi, izqr, izqs, izqg

CHARACTER (LEN=255)      ::  &
  yzerrmsg

CHARACTER (LEN=25)       ::  yzroutine


REAL    (KIND=wp)        ::  &
  zgat, zgct,                & !
  zag, zas,                  & !
  z2ag, z2cg,                & !
  zcg, zcs,                  & !
  zbg, zdg,                  & !
  zz, zzdtr                    !

REAL (KIND=wp)        ::  &
  rho_inv,                &
  r_earth_recip

INTEGER ::  &
  k_s2w_up, k_s2w_low

REAL    (KIND=wp)        :: help1,help2


! Local (automatic) arrays:
! ------------------------
REAL    (KIND=wp)        ::  &
  ztrcor  (ie,je   )      ,  & !
  tmp     (ie,je,ke)      ,  & ! Temporary field in case of perturbed physics
  zc      (ie,je,ke)      ,  & ! Upper main diagonal of ...
  zd      (ie,je,ke)           ! Right hand side of ...

INTEGER (KIND=iintegers) ::  &
  izturb(trcr_get_ntrcr()),  & ! Turbulent mixing for tracers
  izbbc (trcr_get_ntrcr()),  & ! Bottom BC for tracers
  izclip(trcr_get_ntrcr())     ! Clipping for tracers

LOGICAL                  ::  &
  lzst_pert(trcr_get_ntrcr())  ! Stochastic perturbation for this tracer

! Tracer pointers:
!-----------------
REAL (KIND=wp),     POINTER :: &
  ztrcr_now  (:,:,:) => NULL(),        & ! tracer variable at nnow
  ztrcr_new  (:,:,:) => NULL(),        & ! tracer variable at nnew
  ztrcr_tens (:,:,:) => NULL(),        & ! tracer tendency 
  ztrcr_surf (:,:,:) => NULL(),        & ! tracer surface data
  qv_now     (:,:,:) => NULL(),        & ! QV at nnow
  qv_new     (:,:,:) => NULL()           ! QV at nnew

! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine complete_tendencies_trcr
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 0: Check that lvertdiff = TRUE
!------------------------------------------------------------------------------


  izerror   = 0_iintegers
  yzerrmsg  = ''
  yzroutine = 'complete_tendencies_trcr'

  IF ( izdebug >= 10 ) THEN
    WRITE(*,*) "subr. [", TRIM(yzroutine), "] ..."
  END IF

  ! lvertdiff = FALSE is not implemented for the tracers
!MR:  only necessary until the case "itype_vdif=-2" is removed
  IF ( (itype_turb == 3) .AND. (imode_turb >= 2) .AND.     &
       (itype_vdif == -2)                        .AND.     &
       (trcr_get_ntrcr() > 0_iintegers) )                   THEN
    yzerrmsg =  'This vertical diffusion type is not yet ' // &
                'implemented for tracers'
    izerror  =  5000
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

!------------------------------------------------------------------------------
! Section 1: Some preparations
!------------------------------------------------------------------------------

  ! Setting of reciprocal time step
  zzdtr = 1.0_wp / dt

  r_earth_recip = 1.0_wp / r_earth

  ! Retrieve the required metadata
  CALL trcr_meta_get(izerror, T_TURB_ID, izturb)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_meta_get(izerror, T_BBC_ID, izbbc)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_meta_get(izerror, T_CLP_ID, izclip)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

 ! Some more preparations for stochastic perturbed physics

  lzst_pert(:) = .FALSE.   ! no stochastic perturbation for all tracers

  IF (lsppt) THEN          ! do the correct setting for the tracers
                           ! requiring stochastic perturbation
    ! Remark: the following settings for the perturbation
    !         are active only if the turbulent mixing is switched on
    !         for these tracers (check this in the calls to trcr_new in
    !         organize_physics.f90)

    ! Perturb qv in any case
    CALL trcr_get_index(izerror, 'QV', izqv)
    IF (izerror /= 0) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
    lzst_pert(izqv) = .TRUE. ! stochastic perturbation for QV

    ! Perturb also qc and qi in case itype_qxpert_rn is 1 or 2
    IF (itype_qxpert_rn == 1 .OR. itype_qxpert_rn == 2) THEN
      ! Perturb qc
      CALL trcr_get_index(izerror, 'QC', izqc)
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF
      lzst_pert(izqc) = .TRUE. ! stochastic perturbation for QC

      ! Perturb qi
      CALL trcr_get_index(izerror, 'QI', izqi)
      IF (izerror /= 0) THEN
        IF (izerror == T_ERR_NOTFOUND) THEN
          yzerrmsg = "qi cannot be perturbed because is not present in the model."
        ELSE
          yzerrmsg = trcr_errorstr(izerror)
        ENDIF
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF
      lzst_pert(izqi) = .TRUE. ! stochastic perturbation for QI
    ENDIF

    ! Perturb also qr, qs and qg in case itype_qxpert_rn is 2
    IF (itype_qxpert_rn == 2) THEN

      ! Perturb qr
      CALL trcr_get_index(izerror, 'QR', izqr)
      IF (izerror /= 0) THEN
        IF (izerror == T_ERR_NOTFOUND) THEN
          yzerrmsg = "qr cannot be perturbed because is not present in the model."
        ELSE
          yzerrmsg = trcr_errorstr(izerror)
        ENDIF
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF
      lzst_pert(izqr) = .TRUE. ! stochastic perturbation for QR

      ! Perturb qs
      CALL trcr_get_index(izerror, 'QS', izqs)
      IF (izerror /= 0) THEN
        IF (izerror == T_ERR_NOTFOUND) THEN
          yzerrmsg = "qs cannot be perturbed because is not present in the model."
        ELSE
          yzerrmsg = trcr_errorstr(izerror)
        ENDIF
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF
      lzst_pert(izqs) = .TRUE. ! stochastic perturbation for QS

      ! Perturb qg
      CALL trcr_get_index(izerror, 'QG', izqg)
      IF (izerror /= 0) THEN
        IF (izerror == T_ERR_NOTFOUND) THEN
          yzerrmsg = "qg cannot be perturbed because is not present in the model."
          lzst_pert(izqg) = .FALSE. ! no stochastic perturbation for QG, because not present
        ELSE
          ! some other error occured: abort
          yzerrmsg = trcr_errorstr(izerror)
          CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
        ENDIF
      ELSE
        lzst_pert(izqg) = .TRUE. ! stochastic perturbation for QG
      ENDIF
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Setup of tridiagonal matrix systems resulting from the
!             implicit numerical formulation of diffusion of the tracers
!            (if required).
!------------------------------------------------------------------------------

  ! loop over tracers

  DO iztrcr = 1, trcr_get_ntrcr()

    ! get pointer to tracer (at nnew)
    CALL trcr_get(izerror, iztrcr, ptr_tlev=nnew, ptr=ztrcr_new,         &
                  ptr_tens=ztrcr_tens)
    IF (izerror /= 0) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF

    IF ( ( izturb(iztrcr) == T_TURB_1D .OR.                              &
           izturb(iztrcr) == T_TURB_3D ) .AND. lvertdiff ) THEN

      ! get pointer to tracer (at nnow)
      CALL trcr_get(izerror, iztrcr, ptr_tlev=nnow, ptr=ztrcr_now)
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF
    
      ! get pointer to associated surface field (if present)
      IF ( izbbc(iztrcr) == T_BBC_SURF_VAL  ) THEN
        CALL trcr_meta_get(izerror, iztrcr, 'SURF_FIELD',ztrcr_surf)
        IF (izerror /= 0) THEN
          yzerrmsg = trcr_errorstr(izerror)
          CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
        ENDIF
      ENDIF

      ! Computation of vertical diffusion
      ! ---------------------------------
      ! First, the matrix elements a(k), b(k) and c(k) of the coefficient
      ! matrix are set (these are the same for all tracers).
      ! The right hand side is stored on d(k).

      DO j = jstart, jend

        IF (lzst_pert(iztrcr)) THEN   ! for SPPT tendencies perturbation
          DO k = 1, ke
            DO i = istart, iend
              tmp(i,j,k) = ztrcr_new(i,j,k)
            ENDDO
          ENDDO
        ENDIF

        ! Top layer
        DO i = istart, iend
          zgct      = - zkh(i,j,2) * zsqrtgrho_r_s(i,j,1)
          zcg       = zgct*a1t(2)
          zcs       = zgct*a2t(2)
          zbg       = zzdtr - zcg
          zdg       = zzdtr * ztrcr_new(i,j,1) + ztrcr_tens(i,j,1)            &
                    - zcs * ( ztrcr_now(i,j,2) - ztrcr_now(i,j,1) )
          zc(i,j,1) = zcg / zbg
          zd(i,j,1) = zdg / zbg
        ENDDO
   
     
        ! The layers from k=2 to k=ke-1
        IF ( l3dturb_metr ) THEN
          DO k = 2, ke-1
            DO i = istart, iend

              rho_inv = 1.0_wp / rho(i,j,k)

              ! the explicit counterparts of z2ag, z2cg are in ttens (via 'explicit_horizontal_diffusion'):
              z2ag = rho_inv * Kvh_rho_at_w(i,j,k)   * beta_1(k) *                                       &
                 &     ( acrlat(j,1)   * dzeta_dlam(i,j,k) * acrlat(j,1)   * dzeta_dlam_at_w(i,j,k)      &
                 &     + r_earth_recip * dzeta_dphi(i,j,k) * r_earth_recip * dzeta_dphi_at_w(i,j,k) )

              z2cg = rho_inv * Kvh_rho_at_w(i,j,k+1) * beta_1(k+1) *                                     &
                 &     ( acrlat(j,1)   * dzeta_dlam(i,j,k) * acrlat(j,1)   * dzeta_dlam_at_w(i,j,k+1)    &
                 &     + r_earth_recip * dzeta_dphi(i,j,k) * r_earth_recip * dzeta_dphi_at_w(i,j,k+1) )

              zgat      = - zkh(i,j,k  ) * zsqrtgrho_r_s(i,j,k)
              zgct      = - zkh(i,j,k+1) * zsqrtgrho_r_s(i,j,k)

              zag       = zgat*a1t(k)   - z2ag
              zas       = zgat*a2t(k)
              zcg       = zgct*a1t(k+1) - z2cg
              zcs       = zgct*a2t(k+1)
              zbg       = zzdtr - zag - zcg
              zdg       = zzdtr * ztrcr_new(i,j,k) + ztrcr_tens(i,j,k)          &
                 - zas * ( ztrcr_now(i,j,k-1) - ztrcr_now(i,j,k) )       &
                 - zcs * ( ztrcr_now(i,j,k+1) - ztrcr_now(i,j,k) )
              zz        = 1.0_wp / ( zbg - zag*zc(i,j,k-1) )
              zc(i,j,k) = zcg * zz
              zd(i,j,k) = ( zdg -zag*zd(i,j,k-1) ) * zz
            ENDDO
          ENDDO
        ELSE
          DO k = 2, ke-1
            DO i = istart, iend

              zgat      = - zkh(i,j,k  ) * zsqrtgrho_r_s(i,j,k)
              zgct      = - zkh(i,j,k+1) * zsqrtgrho_r_s(i,j,k)

              zag       = zgat*a1t(k)
              zas       = zgat*a2t(k)
              zcg       = zgct*a1t(k+1)
              zcs       = zgct*a2t(k+1)
              zbg       = zzdtr - zag - zcg
              zdg       = zzdtr * ztrcr_new(i,j,k) + ztrcr_tens(i,j,k)          &
                        - zas * ( ztrcr_now(i,j,k-1) - ztrcr_now(i,j,k) )       &
                        - zcs * ( ztrcr_now(i,j,k+1) - ztrcr_now(i,j,k) )
              zz        = 1.0_wp / ( zbg - zag*zc(i,j,k-1) )
              zc(i,j,k) = zcg * zz
              zd(i,j,k) = ( zdg -zag*zd(i,j,k-1) ) * zz
            ENDDO
          ENDDO
        END IF
   
        ! the bottom layer (different version for different lower BC, only
        !   the line zdg=... actually changes, but it's necessary to pull
        !   out the IF statement for vectorization)
        IF     ( izbbc(iztrcr) ==  T_BBC_ZEROFLUX ) THEN
          DO i = istart, iend
            zgat      = - zkh (i,j,ke)*zsqrtgrho_r_s(i,j,ke)
            zgct      = - ztch(i,j   )*zsqrtgrho_r_s(i,j,ke)
            zag       = zgat*a1t(ke)
            zas       = zgat*a2t(ke)
            zcg       = zgct*za1t_surf
            zcs       = zgct*za2t_surf
            zbg       = zzdtr - zag - zcg
            zdg     = zzdtr * ztrcr_new(i,j,ke  ) + ztrcr_tens(i,j,ke)        &
                    - zas * ( ztrcr_now(i,j,ke-1) - ztrcr_now(i,j,ke) )          
            zz        = 1.0_wp / ( zbg  - zag*zc(i,j,ke-1) )
            ztrcr_new(i,j,ke) = ( zdg - zag*zd(i,j,ke-1) ) * zz
          ENDDO
        ELSEIF ( izbbc(iztrcr) == T_BBC_ZEROVAL   ) THEN
          DO i = istart, iend
            zgat      = - zkh (i,j,ke)*zsqrtgrho_r_s(i,j,ke)
            zgct      = - ztch(i,j   )*zsqrtgrho_r_s(i,j,ke)
            zag       = zgat*a1t(ke)
            zas       = zgat*a2t(ke)
            zcg       = zgct*za1t_surf
            zcs       = zgct*za2t_surf
            zbg       = zzdtr - zag - zcg
            zdg     = zzdtr * ztrcr_new(i,j,ke  ) + ztrcr_tens(i,j,ke)        &
                    - zas * ( ztrcr_now(i,j,ke-1) - ztrcr_now(i,j,ke) )       &
                    + zcs *   ztrcr_now(i,j,ke  )
            zz        = 1.0_wp / ( zbg  - zag*zc(i,j,ke-1) )
            ztrcr_new(i,j,ke) = ( zdg - zag*zd(i,j,ke-1) ) * zz
          ENDDO
        ELSEIF ( izbbc(iztrcr) == T_BBC_SURF_VAL  ) THEN
          DO i = istart, iend
            zgat      = - zkh (i,j,ke)*zsqrtgrho_r_s(i,j,ke)
            zgct      = - ztch(i,j   )*zsqrtgrho_r_s(i,j,ke)
            zag       = zgat*a1t(ke)
            zas       = zgat*a2t(ke)
            zcg       = zgct*za1t_surf
            zcs       = zgct*za2t_surf
            zbg       = zzdtr - zag - zcg
            zdg     = zzdtr * ztrcr_new(i,j,ke  ) + ztrcr_tens(i,j,ke)        &
                    - zas * ( ztrcr_now(i,j,ke-1) - ztrcr_now(i,j,ke) )       &
                    + zcs *   ztrcr_now(i,j,ke  ) - zgct * ztrcr_surf(i,j,nnow)
            zz        = 1.0_wp / ( zbg  - zag*zc(i,j,ke-1) )
            ztrcr_new(i,j,ke) = ( zdg - zag*zd(i,j,ke-1) ) * zz
          ENDDO
        ELSE
          yzerrmsg =  'This BBC option is not valid for the tracers'      
          izerror  =  5000 
          CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
        ENDIF
    
        ! Backsubstitution and storage of the complete slow tendencies
        DO k = ke-1, 1, -1
          DO i = istart, iend
            ztrcr_new(i,j,k) = zd(i,j,k) - zc(i,j,k) * ztrcr_new(i,j,k+1)
          ENDDO
        ENDDO
    
        IF (lzst_pert(iztrcr)) THEN   ! SPPT tendencies perturbation
          DO k = ke, 1, -1
!CDIR ON_ADB(ztrcr_new)
            DO i = istart, iend
              ztrcr_new(i,j,k) = tmp(i,j,k) + (  ztrcr_new(i,j,k)             &
                                               - tmp(i,j,k) ) *pertstoph(i,j,k)
            ENDDO
          ENDDO
        ENDIF
      ENDDO

    ELSE

      ! no vertical diffusion, simply add the tendency
      DO k = 1, ke
        DO j = jstart, jend
          DO i = istart, iend
            ztrcr_new(i,j,k) = ztrcr_new(i,j,k) + dt * ztrcr_tens(i,j,k)
          ENDDO
        ENDDO
      ENDDO

    ENDIF

!MR: clipping also if ".NOT.lvertdiff"
!US: this can now be executed, if ltur=.FALSE.: this behaviour is different than before
    ! Clipping if required
    IF ( ( izturb(iztrcr) == T_TURB_1D .OR.                              &
           izturb(iztrcr) == T_TURB_3D ) .AND. &
           izclip(iztrcr) == T_CLP_ON  ) THEN
       CALL clipping( ztrcr_new(:,:,:), ie, je, ke )
    ENDIF

  ENDDO ! loop over tracers

  !----------------------------------------------------------------------------
  ! Section 3: Calculation of the surface moisture flux 'qvsflx'
  !            The latent heat flux is integrated in time
  !----------------------------------------------------------------------------

! FUO TODO: this could be moved out of complete_tendencies_trcr since
!           it is currently not a diagnostic which does not generally apply to all
!           tracers.

!MR:<
  IF (lvertdiff) THEN
!MR>

  ! retrieve the required tracers (at specified timelevel)
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnew, ptr = qv_new)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnow, ptr = qv_now)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  DO j =jstart, jend
    DO i = istart , iend
      qvsflx(i,j)  =  - ztch(i,j) * &
                     ( za2t_surf*(qv_s(i,j,nnow) - qv_now(i,j,ke)) +  &
                       za1t_surf*(qv_s(i,j,nnew) - qv_new(i,j,ke)) )
      ! For the latent heat flux, we do not take into account the
      !  implicit time weighting in case of imode_turb=0 for simplicity,
      !  so no "za2t_surf" and "za1t_surf" and no "nnew"-branch:
      IF (t_g(i,j,nnow) >= t0_melt) THEN
        lhfl_s(i,j) = lh_v * qvsflx(i,j)
      ELSE
        lhfl_s(i,j) = lh_s * qvsflx(i,j)
      END IF
!MR:
!The inclusion of a freezing rate is not consistent, as the real latent heat flux
!depends on the freezed surface fraction!

    ENDDO
  ENDDO

  ENDIF  ! lvertdiff

  !----------------------------------------------------------------------------
  ! End of subroutine complete_tendencies_trcr 
  !----------------------------------------------------------------------------

END SUBROUTINE complete_tendencies_trcr

!==============================================================================
!==============================================================================

SUBROUTINE calc_K_rho_at_y

!------------------------------------------------------------------------------
!
! Description:
!   allocate and calculate those averaged diffusion coefficient fields
!   that are used both in explicit_horizontal_diffusion 
!   and in complete_tendencies_uvwt_3D
!   (this is only needed for 3D-diffusion)
!
! output:
!   3D fields
!     Kvm_rho_at_s, Kvm_rho_at_uw, Kvm_rho_at_vw,
!     Kvh_rho_pi_at_w
!   switch
!     l_diffusion_tracers
!
!------------------------------------------------------------------------------

  USE src_tracer          ,     ONLY :  &
    trcr_meta_get, trcr_get_ntrcr

  INTEGER :: i, j, k

  INTEGER ::  &
     k_s2w_up, k_s2w_low

  ! Auxiliary fields:
  REAL    (KIND=wp), ALLOCATABLE    ::  &  
     Kvm_rho_at_w(:,:,:),  &
     rho_pi(:,:)

  REAL    (KIND=wp) ::  z_pi
  
  INTEGER (KIND=iintegers) ::  &
    izturb(trcr_get_ntrcr())    ! Turbulent mixing for tracers

  INTEGER :: &
    iztrcr, izerror


  IF ( izdebug >= 10 ) THEN
    WRITE(*,*) "subr. [calc_K_rho_at_y] ..."
  END IF

  ALLOCATE( Kvm_rho_at_s    (ie, je, ke  ) )
  ALLOCATE( Kvm_rho_at_uw   (ie, je, ke+1) )
  ALLOCATE( Kvm_rho_at_vw   (ie, je, ke+1) )
  ALLOCATE( Kvh_rho_pi_at_w (ie, je, ke+1) )

  ALLOCATE( Kvm_rho_at_w    (ie, je, ke+1) )   ! here only auxiliary field


  ! --- determine Kvm_rho_at_s, ..._uw, ..._vw ---

  DO k=1, ke
    ! for averaging from the scalar position to the w-position
    k_s2w_up  = MAX( k-1, 1  )
    k_s2w_low = MIN( k  , ke )

    ! Kvm * rho at the two grid positions 'scalar' and 'w' need the least averaging:
    DO j = 1, je
      DO  i = 1, ie
        Kvm_rho_at_s (i,j,k) = rho(i,j,k) * 0.5_wp * ( tkvm(i,j,k) + tkvm(i,j,k+1) ) 
      END DO
    END DO

    DO j=1, je
      DO i=1, ie
        Kvm_rho_at_w (i,j,k) = 0.5_wp * ( rho(i,j,k_s2w_up) + rho(i,j,k_s2w_low) ) * tkvm(i,j,k)
      END DO
    END DO

  END DO
  Kvm_rho_at_w (:,:,ke+1) = rho(:,:,ke) * tkvm(:,:,ke+1)

  DO k=1, ke1
    DO j = 1, je
      DO  i = 1, ie-1
        Kvm_rho_at_uw(i,j,k) = 0.5_wp * ( Kvm_rho_at_w (i,j,k) + Kvm_rho_at_w (i+1,j  ,k) )
      END DO
    END DO
    Kvm_rho_at_uw(ie,:,k) = 0.0_wp

    DO j = 1, je-1
      DO  i = 1, ie
        Kvm_rho_at_vw(i,j,k) = 0.5_wp * ( Kvm_rho_at_w (i,j,k) + Kvm_rho_at_w (i  ,j+1,k) )
      END DO
    END DO
    Kvm_rho_at_vw(:,je,k) = 0.0_wp

  END DO

  DEALLOCATE( Kvm_rho_at_w  )


  ! --- determine Kvh_rho_pi_at_w ---

  ALLOCATE( rho_pi(ie, ke) )

  DO j = 1, je

    DO k=1, ke
      DO  i = 1, ie
        z_pi = ( ( p0(i,j,k) + pp(i,j,k,nnow) ) * p0ref_recip )**rdocp
        rho_pi(i,k) = rho(i,j,k) * z_pi
      END DO
    END DO

    DO k=1, ke+1
      ! for averaging from the scalar position to the w-position
      k_s2w_up  = MAX( k-1, 1  )
      k_s2w_low = MIN( k  , ke )

      DO  i = 1, ie
        Kvh_rho_pi_at_w (i,j,k) = 0.5_wp * ( rho_pi(i,k_s2w_up) + rho_pi(i,k_s2w_low) ) * tkvh(i,j,k)
      END DO
    END DO

  END DO

  DEALLOCATE( rho_pi )


  ! --- if diffusion is applied to tracers, then determine Kvh_rho_at_w ---

  l_diffusion_tracers = .FALSE.

  CALL trcr_meta_get(izerror, T_TURB_ID, izturb)

  DO iztrcr = 1, trcr_get_ntrcr()

    IF ( ( izturb(iztrcr) == T_TURB_1D .OR.                              &
           izturb(iztrcr) == T_TURB_3D )    .AND. lvertdiff ) THEN
      l_diffusion_tracers = .TRUE.
    END IF
  END DO


  IF ( l_diffusion_tracers ) THEN

    ALLOCATE( Kvh_rho_at_w (ie, je, ke+1) )

    DO k=1, ke+1
      ! for averaging from the w-position to the scalar position:
      !k_w2s_up  = MAX( 1,  k   )
      !k_w2s_low = MIN( ke, k+1 )

      ! for averaging from the scalar position to the w-position
      k_s2w_up  = MAX( k-1, 1  )
      k_s2w_low = MIN( k  , ke )

      DO j = 1, je
        DO  i = 1, ie
          Kvh_rho_at_w(i,j,k) = 0.5_wp * ( rho(i,j,k_s2w_up) + rho(i,j,k_s2w_low) ) * tkvh(i,j,k)
        ENDDO
      ENDDO
    ENDDO

  END IF


  IF ( l_diffusion_tracers .OR. lprog_tke ) THEN

    ALLOCATE( Kvh_rho_at_s (ie, je, ke) )

    DO k=1, ke
      DO j = 1, je
        DO  i = 1, ie
          Kvh_rho_at_s (i,j,k) = rho(i,j,k) * 0.5_wp * ( tkvh(i,j,k) + tkvh(i,j,k+1) ) 
        END DO
      END DO
    END DO

  END IF


END SUBROUTINE calc_K_rho_at_y

!==============================================================================
!==============================================================================

SUBROUTINE deallocate_K_rho_at_y

  IF ( izdebug >= 10 ) THEN
    WRITE(*,*) "subr. [deallocate_K_rho_at_y] ..."
  END IF

  DEALLOCATE( Kvm_rho_at_s  )
  DEALLOCATE( Kvm_rho_at_uw )
  DEALLOCATE( Kvm_rho_at_vw )

  DEALLOCATE( Kvh_rho_pi_at_w )

  IF ( l_diffusion_tracers ) THEN
    DEALLOCATE( Kvh_rho_at_w )
  END IF

  IF ( l_diffusion_tracers .OR. lprog_tke ) THEN
    DEALLOCATE( Kvh_rho_at_s )
  END IF

END SUBROUTINE deallocate_K_rho_at_y

!==============================================================================
!==============================================================================

SUBROUTINE set_beta_xxx_for_3D_diffusion

!------------------------------------------------------------------------------
!
! Description:
! allocate and set all implicit weighting coefficients beta_xxx
!
!------------------------------------------------------------------------------

  INTEGER :: k

  ALLOCATE( beta_uu1(ke1) )
  ALLOCATE( beta_uu2(ke1) )
  ALLOCATE( beta_uu3(ke1) )
  ALLOCATE( beta_vv1(ke1) )
  ALLOCATE( beta_vv2(ke1) )
  ALLOCATE( beta_vv3(ke1) )
  ALLOCATE( beta_ww1(ke1) )
  ALLOCATE( beta_ww2(ke1) )
  ALLOCATE( beta_1(ke1) )
  ALLOCATE( beta_2(ke1) )

  DO k=1, ke+1
    beta_uu1(k) = a1t(k)
    beta_uu2(k) = a1t(k)
    beta_uu3(k) = a1t(k)
    beta_vv1(k) = a1t(k)
    beta_vv2(k) = a1t(k)
    beta_vv3(k) = a1t(k)
    beta_ww1(k) = a1t(k)
    beta_ww2(k) = a1t(k)
    beta_1(k)   = a1t(k)
    beta_2(k)   = a1t(k)
  END DO

END SUBROUTINE set_beta_xxx_for_3D_diffusion

!==============================================================================
!==============================================================================

SUBROUTINE deallocate_beta_xxx_for_3D_diffusion

  DEALLOCATE( beta_uu1 )
  DEALLOCATE( beta_uu2 )
  DEALLOCATE( beta_uu3 )
  DEALLOCATE( beta_vv1 )
  DEALLOCATE( beta_vv2 )
  DEALLOCATE( beta_vv3 )
  DEALLOCATE( beta_ww1 )
  DEALLOCATE( beta_ww2 )
  DEALLOCATE( beta_1 )
  DEALLOCATE( beta_2 )

END SUBROUTINE deallocate_beta_xxx_for_3D_diffusion

!==============================================================================
!==============================================================================


!option! -pvctl _on_adb
SUBROUTINE explicit_horizontal_diffusion

  !------------------------------------------------------------------------------
  !
  ! Description:
  ! calculates tendencies due to the calculation of 
  ! horizontal flux-divergences and optionally metric corrections
  !
  ! output:
  ! utens, vtens, wtens, ttens (and possibly tendencies for tracer variables or tketens)
  !
  ! input:
  ! utens, vtens, wtens, ttens (and possibly tendencies for tracer variables or tketens)
  ! Kvm, Kvh
  !
  ! references:
  ! M. Baldauf (2005): "The Coordinate Transformations of the 3-Dimensional
  !        Turbulent Diffusion in LMK", COSMO-Newsletter No. 5, 132-140
  ! M. Baldauf (2006): "Implementation of the 3D-Turbulence Metric Terms in LMK"
  !           COSMO-Newsletter No. 6, 44-50
  !
  !------------------------------------------------------------------------------

  USE data_runcontrol, ONLY: l3dturb_metr

  IMPLICIT NONE


  ! Local scalars:
  !---------------
  !REAL    (KIND=wp)     ::  &
  ! zhfdx, zhfdy, zhfkh

  INTEGER (KIND=iintegers) :: &
    i, j, k,                &
    k_s2w_up, k_s2w_low,    &
    iztrcr, nztrcr_3dturb,  &
    icount,                 &
    izerror

  CHARACTER (LEN=255)      :: &
    yzerrmsg

  CHARACTER (LEN=50)       :: yzroutine

  !LOGICAL                  :: &
  !  lkge2_prog_tke

  REAL    (KIND=wp), ALLOCATABLE, DIMENSION(:,:,:)     ::  &
    Kvh_rho_at_u,   &
    Kvh_rho_at_v,   &
    Kvh_rho_at_uw,  &
    Kvh_rho_at_vw


  REAL    (KIND=wp), ALLOCATABLE, DIMENSION(:,:,:)     ::  &
    Kvh_rho_pi_at_s,  &
    Kvh_rho_pi_at_u,  &
    Kvh_rho_pi_at_v,  &
    work

  REAL    (KIND=wp), ALLOCATABLE, DIMENSION(:,:,:)     ::  &
    ztke, ztketens

  REAL    (KIND=wp)   ::  &
    Kvm_rho_at_uv, &
    z_pi

  REAL    (KIND=wp), ALLOCATABLE, DIMENSION(:,:,:)     ::  &
    rho_inv, rho_inv_at_u, rho_inv_at_v, rho_inv_at_w 

  REAL    (KIND=wp), ALLOCATABLE, DIMENSION(:,:)     ::    &
    work_ik, work_jk

  REAL    (KIND=wp)   ::                          &
    ds_dlam_at_u, ds_dlam_at_w,                   &
    ds_dphi_at_v, ds_dphi_at_w,                   &
    ds_dzeta, ds_dzeta_at_u,  ds_dzeta_at_v, ds_dzeta_at_w 

  REAL    (KIND=wp)   ::                          &
    dtke_dlam_at_s,  dtke_dlam_at_uw,             &
    dtke_dphi_at_s,  dtke_dphi_at_vw,             &
    dtke_dzeta_at_s, dtke_dzeta_at_uw, dtke_dzeta_at_vw

  REAL    (KIND=wp)   ::                          &
    du_dlam_at_s,  du_dlam_at_uw, du_dzeta,       &
    du_dphi_at_uv, du_dphi_at_uw, du_dphi_at_vw,  &
    du_dzeta_at_s, du_dzeta_at_uv, du_dzeta_at_uw, du_dzeta_at_vw

  REAL    (KIND=wp)   ::                          &
    dv_dlam_at_uv, dv_dlam_at_vw, dv_dlam_at_uw,  &
    dv_dphi_at_s,  dv_dphi_at_vw,                 &
    dv_dzeta, dv_dzeta_at_s, dv_dzeta_at_uw,      &
    dv_dzeta_at_vw, dv_dzeta_at_uv

  REAL    (KIND=wp)   ::            &
    dw_dlam_at_s, dw_dlam_at_uw,    &
    dw_dphi_at_s, dw_dphi_at_vw,    &
    dw_dzeta_at_uw, dw_dzeta_at_vw, dw_dzeta_at_s

  REAL    (KIND=wp)   ::                                  &
    dzeta_dlam_at_v, dzeta_dlam_at_uv, dzeta_dlam_at_vw,  &
    dzeta_dphi_at_u, dzeta_dphi_at_uv, dzeta_dphi_at_uw,  &
    dzeta_dz_at_uw  , dzeta_dz_at_vw

  ! Local arrays:
  !--------------
  REAL    (KIND=wp)     ::         &
    ztau12_at_uv (ie,je,ke),       & !
    ztau13_at_uw (ie,je,ke+1),     & !
    ztau23_at_vw (ie,je,ke+1)        !

  INTEGER (KIND=iintegers) ::           &
    iztrcr_3dturb(trcr_get_ntrcr())   , &
    izturb       (trcr_get_ntrcr())

  REAL    (KIND=wp)     ::  &
    r_earth_recip

  CHARACTER (LEN=50)       :: y_tracer_name

  ! Tracer pointers:
  ! ----------------

  REAL (KIND=wp),     POINTER ::        &
    ztrcr_now  (:,:,:) => NULL(),       & ! tracer variable at nnow
    ztrcr_tens (:,:,:) => NULL()          ! tracer tendency

  !USUS should be ok for single and double precision
  REAL(KIND=wp),     PARAMETER :: ztkemin = 1.0E-5_wp  ! lower security limit of sqrt(2*TKE)

  ! End of header
  !==============================================================================

  ! -------- (0) preparations --------------------------------------------

  izerror   = 0_iintegers
  yzerrmsg  = ''
  yzroutine = 'explicit_horizontal_diffusion'

  IF ( izdebug >= 10 ) THEN
    WRITE(*,*) "subr. [", TRIM(yzroutine), "] ..."
  END IF

  nztrcr_3dturb = 0_iintegers
  icount        = 1_iintegers

  ! retrieve the required metadata
  CALL trcr_meta_get(izerror, T_TURB_ID, izturb)
  IF (izerror /= 0_iintegers) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  r_earth_recip = 1.0_wp / r_earth

  ALLOCATE( work_ik( ie, ke+1) )
  ALLOCATE( work_jk( je, ke+1) )


  ! === utens =========================================

  ALLOCATE( rho_inv_at_u (ie, je, ke) )

  DO k = 1, ke
    DO j = jstartu, jendu
      DO  i = istartu, iendu
        rho_inv_at_u(i,j,k) = 2.0_wp / ( rho(i,j,k) + rho(i+1,j,k) )
      END DO
    END DO
  END DO

  ! --- utens: contributions by tau_11 ---

  DO j = jstartu, jendu

    DO k = 1, ke
      DO  i = istartu, iendu+1

        du_dlam_at_s = eddlon * ( u(i,j,k,nnow) - u(i-1,j,k,nnow) )

        ! tau11 at s-position:
        work_ik(i,k) = -2.0_wp * Kvm_rho_at_s(i,j,k) * acrlat(j,1) * du_dlam_at_s
      END DO
    END DO

    IF ( l3dturb_metr ) THEN

      k = 1
      DO  i = istartu, iendu+1

        du_dzeta_at_s = 0.5_wp * ( ( u(i-1,j,k+1,nnow) - u(i-1,j,k,nnow) )    &
          &                      + ( u(i  ,j,k+1,nnow) - u(i  ,j,k,nnow) ) ) 

        ! tau11 at s-position:
        work_ik(i,k) = work_ik(i,k)                                         &
          &      - 2.0_wp * Kvm_rho_at_s(i,j,k)                             &
          &            * acrlat(j,1) * dzeta_dlam(i,j,k) * du_dzeta_at_s
      END DO

      DO k = 2, ke-1
        DO  i = istartu, iendu+1

          du_dzeta_at_s = 0.25_wp * ( ( u(i-1,j,k+1,nnow) - u(i-1,j,k-1,nnow) )    &
            &                       + ( u(i  ,j,k+1,nnow) - u(i  ,j,k-1,nnow) ) ) 

          ! tau11 at s-position:
          work_ik(i,k) = work_ik(i,k)                                         &
            &      - 2.0_wp * Kvm_rho_at_s(i,j,k)                             &
            &            * acrlat(j,1) * dzeta_dlam(i,j,k) * du_dzeta_at_s
        END DO
      END DO

      k = ke
      DO  i = istartu, iendu+1

        du_dzeta_at_s = 0.5_wp * ( ( u(i-1,j,k,nnow) - u(i-1,j,k-1,nnow) )    &
          &                      + ( u(i  ,j,k,nnow) - u(i  ,j,k-1,nnow) ) ) 

        ! tau11 at s-position:
        work_ik(i,k) = work_ik(i,k)                                       &
          &      - 2.0_wp * Kvm_rho_at_s(i,j,k)                           &
          &            * acrlat(j,1) * dzeta_dlam(i,j,k) * du_dzeta_at_s
      END DO

    END IF

    !MB  DO k = 1, ke
    DO k = 2, ke-1
      DO i = istartu, iendu

        utens(i,j,k)= utens(i,j,k) - rho_inv_at_u(i,j,k)     &
          &    * acrlat(j,1) * ( work_ik(i+1,k) - work_ik(i,k) ) * eddlon

      END DO
    END DO

  END DO   ! DO j


  IF ( l3dturb_metr ) THEN

    DO j = jstartu, jendu

      ! tau11 at uw-position:
      k=1
      DO  i = istartu, iendu

        du_dlam_at_uw = 0.5_wp * ( u(i+1,j,k,nnow) - u(i-1,j,k,nnow) ) * eddlon

        du_dzeta = u(i,j,k+1,nnow) - u(i,j,k,nnow)       ! ??????

        ! tau11 at uw-position:
        work_ik(i,k) = -2.0_wp * Kvm_rho_at_uw(i,j,k) * acrlat(j,1)          &
          &  * ( du_dlam_at_uw                                               &
          &    + (1.0_wp-beta_uu1(k)) * dzeta_dlam_at_uw(i,j,k) * du_dzeta )
      END DO

      DO k=2, ke
        DO  i = istartu, iendu

          du_dlam_at_uw = 0.25_wp * eddlon                           &
            &       * ( ( u(i+1,j,k-1,nnow) - u(i-1,j,k-1,nnow) )    &
            &         + ( u(i+1,j,k  ,nnow) - u(i-1,j,k  ,nnow) ) )

          du_dzeta = u(i,j,k,nnow) - u(i,j,k-1,nnow)

          ! tau11 at uw-position:
          work_ik(i,k) = -2.0_wp * Kvm_rho_at_uw(i,j,k) * acrlat(j,1)          &
            &  * ( du_dlam_at_uw                                               &
            &    + (1.0_wp-beta_uu1(k)) * dzeta_dlam_at_uw(i,j,k) * du_dzeta )
        END DO
      END DO

      !MB  DO k=1, ke-1
      DO k=2, ke-1
        DO i = istartu, iendu

          utens(i,j,k)= utens(i,j,k) - rho_inv_at_u(i,j,k)                      &
            &  * acrlat(j,1) * dzeta_dlam_at_u(i,j,k) * ( work_ik(i,k+1) - work_ik(i,k) )

        END DO
      END DO

    END DO    ! DO j

  END IF

  ! --- utens: contributions by tau_12 ---

  ! tau12 at uv-position (can be used also for v-tens!):
  DO k=1, ke
    DO j = MIN(jstartu-1, jstartv), MAX(jendu, jendv)
      DO i = MIN(istartu, istartv-1), MAX(iendu, iendv)

        dv_dlam_at_uv = ( v(i+1,j  ,k,nnow) - v(i,j,k,nnow) ) * eddlon
        du_dphi_at_uv = ( u(i  ,j+1,k,nnow) - u(i,j,k,nnow) ) * eddlat

        Kvm_rho_at_uv = 0.25_wp * ( Kvm_rho_at_s (i,j  ,k) + Kvm_rho_at_s (i+1,j  ,k)  &
          &                       + Kvm_rho_at_s (i,j+1,k) + Kvm_rho_at_s (i+1,j+1,k) )

        ztau12_at_uv(i,j,k) = - Kvm_rho_at_uv           &
          &     * (    acrlat(j,2) * dv_dlam_at_uv      &
          &        + r_earth_recip * du_dphi_at_uv )
      END DO
    END DO
  END DO

  IF ( l3dturb_metr ) THEN

    k=1
    DO j = MIN(jstartu-1, jstartv), MAX(jendu, jendv)
      DO i = MIN(istartu, istartv-1), MAX(iendu, iendv)

        dv_dzeta_at_uv = 0.5_wp * ( ( v(i  ,j,k+1,nnow) - v(i  ,j,k,nnow) )    &
          &                       + ( v(i+1,j,k+1,nnow) - v(i+1,j,k,nnow) ) )

        du_dzeta_at_uv = 0.5_wp * ( ( u(i,j  ,k+1,nnow) - u(i,j  ,k,nnow) )    &
          &                       + ( u(i,j+1,k+1,nnow) - u(i,j+1,k,nnow) ) )

        dzeta_dlam_at_uv = 0.25_wp * ( dzeta_dlam(i,j  ,k) + dzeta_dlam(i+1,j  ,k)   &
          &                          + dzeta_dlam(i,j+1,k) + dzeta_dlam(i+1,j+1,k) )

        dzeta_dphi_at_uv = 0.25_wp * ( dzeta_dphi(i,j  ,k) + dzeta_dphi(i+1,j  ,k)   &
          &                          + dzeta_dphi(i,j+1,k) + dzeta_dphi(i+1,j+1,k) )

        Kvm_rho_at_uv = 0.25_wp * ( Kvm_rho_at_s (i,j  ,k) + Kvm_rho_at_s (i+1,j  ,k)   &
          &                      + Kvm_rho_at_s (i,j+1,k) + Kvm_rho_at_s (i+1,j+1,k) )

        ztau12_at_uv(i,j,k) = ztau12_at_uv(i,j,k) - Kvm_rho_at_uv        &
          &    * ( dzeta_dlam_at_uv * acrlat(j,2)   * dv_dzeta_at_uv     &
          &      + dzeta_dphi_at_uv * r_earth_recip * du_dzeta_at_uv )  

      END DO
    END DO


    DO k=2, ke-1
      DO j = MIN(jstartu-1, jstartv), MAX(jendu, jendv)
        DO i = MIN(istartu, istartv-1), MAX(iendu, iendv)

          dv_dzeta_at_uv = 0.25_wp * ( ( v(i  ,j,k+1,nnow) - v(i  ,j,k-1,nnow) )    &
            &                        + ( v(i+1,j,k+1,nnow) - v(i+1,j,k-1,nnow) ) )

          du_dzeta_at_uv = 0.25_wp * ( ( u(i,j  ,k+1,nnow) - u(i,j  ,k-1,nnow) )    &
            &                        + ( u(i,j+1,k+1,nnow) - u(i,j+1,k-1,nnow) ) )

          dzeta_dlam_at_uv = 0.25_wp * ( dzeta_dlam(i,j  ,k) + dzeta_dlam(i+1,j  ,k)   &
            &                          + dzeta_dlam(i,j+1,k) + dzeta_dlam(i+1,j+1,k) )

          dzeta_dphi_at_uv = 0.25_wp * ( dzeta_dphi(i,j  ,k) + dzeta_dphi(i+1,j  ,k)   &
            &                          + dzeta_dphi(i,j+1,k) + dzeta_dphi(i+1,j+1,k) )

          Kvm_rho_at_uv = 0.25_wp * ( Kvm_rho_at_s (i,j  ,k) + Kvm_rho_at_s (i+1,j  ,k)   &
            &                      + Kvm_rho_at_s (i,j+1,k) + Kvm_rho_at_s (i+1,j+1,k) )

          ztau12_at_uv(i,j,k) = ztau12_at_uv(i,j,k) - Kvm_rho_at_uv        &
            &    * ( dzeta_dlam_at_uv * acrlat(j,2)   * dv_dzeta_at_uv     &
            &      + dzeta_dphi_at_uv * r_earth_recip * du_dzeta_at_uv )  

        END DO
      END DO
    END DO

    k=ke
    DO j = MIN(jstartu-1, jstartv), MAX(jendu, jendv)
      DO i = MIN(istartu, istartv-1), MAX(iendu, iendv)

        dv_dzeta_at_uv = 0.5_wp * ( ( v(i  ,j,k,nnow) - v(i  ,j,k-1,nnow) )    &
          &                       + ( v(i+1,j,k,nnow) - v(i+1,j,k-1,nnow) ) )

        du_dzeta_at_uv = 0.5_wp * ( ( u(i,j  ,k,nnow) - u(i,j  ,k-1,nnow) )    &
          &                       + ( u(i,j+1,k,nnow) - u(i,j+1,k-1,nnow) ) )

        dzeta_dlam_at_uv = 0.25_wp * ( dzeta_dlam(i,j  ,k) + dzeta_dlam(i+1,j  ,k)   &
          &                          + dzeta_dlam(i,j+1,k) + dzeta_dlam(i+1,j+1,k) )

        dzeta_dphi_at_uv = 0.25_wp * ( dzeta_dphi(i,j  ,k) + dzeta_dphi(i+1,j  ,k)   &
          &                          + dzeta_dphi(i,j+1,k) + dzeta_dphi(i+1,j+1,k) )

        Kvm_rho_at_uv = 0.25_wp * ( Kvm_rho_at_s (i,j  ,k) + Kvm_rho_at_s (i+1,j  ,k)   &
          &                      + Kvm_rho_at_s (i,j+1,k) + Kvm_rho_at_s (i+1,j+1,k) )

        ztau12_at_uv(i,j,k) = ztau12_at_uv(i,j,k) - Kvm_rho_at_uv        &
          &    * ( dzeta_dlam_at_uv * acrlat(j,2)   * dv_dzeta_at_uv     &
          &      + dzeta_dphi_at_uv * r_earth_recip * du_dzeta_at_uv )  

      END DO
    END DO

  END IF


  !MB   DO k = 1, ke
  DO k = 2, ke-1
    DO j = jstartu, jendu
      DO i = istartu, iendu

        utens(i,j,k)= utens(i,j,k) - rho_inv_at_u(i,j,k) *                   &
          &    r_earth_recip * ( ztau12_at_uv(i,j,k) - ztau12_at_uv(i,j-1,k) ) * eddlat

      END DO
    END DO
  END DO


  ! tau12 at uw-position:

  IF ( l3dturb_metr ) THEN

    DO j = jstartu, jendu

      DO k=2, ke
        DO i = istartu, iendu

          dv_dlam_at_uw = 0.25_wp * eddlon *                       &
            &    ( ( v(i+1,j  ,k,  nnow) - v(i,j  ,k  ,nnow) )     &
            &    + ( v(i+1,j-1,k,  nnow) - v(i,j-1,k  ,nnow) )     &
            &    + ( v(i+1,j  ,k-1,nnow) - v(i,j  ,k-1,nnow) )     &
            &    + ( v(i+1,j-1,k-1,nnow) - v(i,j-1,k-1,nnow) ) )

          du_dphi_at_uw = 0.25_wp * eddlat *                       &
            &    ( ( u(i,j+1,k  ,nnow) - u(i,j-1,k  ,nnow) )       &
            &    + ( u(i,j+1,k-1,nnow) - u(i,j-1,k-1,nnow) ) )

          dv_dzeta_at_uw = 0.25_wp * ( ( v(i  ,j  ,k,nnow) - v(i  ,j  ,k-1,nnow) )    &
            &                        + ( v(i+1,j  ,k,nnow) - v(i+1,j  ,k-1,nnow) )    &
            &                        + ( v(i  ,j-1,k,nnow) - v(i  ,j-1,k-1,nnow) )    &
            &                        + ( v(i+1,j-1,k,nnow) - v(i+1,j-1,k-1,nnow) ) )

          du_dzeta_at_uw = u(i,j,k,nnow) - u(i,j,k-1,nnow)

          dzeta_dphi_at_uw = 0.5_wp * ( dzeta_dphi_at_w(i,j,k) + dzeta_dphi_at_w(i+1,j,k) )

          ! tau12 at uw-position:
          work_ik(i,k) = -Kvm_rho_at_uw(i,j,k) *                                                 &
            &     (  acrlat(j,1)   * ( dv_dlam_at_uw                                             &
            &                        + dzeta_dlam_at_uw(i,j,k) * dv_dzeta_at_uw )                &
            &      + r_earth_recip * ( du_dphi_at_uw                                             &
            &                        + (1.0_wp-beta_uu3(k)) * dzeta_dphi_at_uw * du_dzeta_at_uw ) )

        END DO
      END DO

      !MB  DO k=2, ke
      DO k=2, ke-1
        DO i = istartu, iendu

          dzeta_dphi_at_u = 0.5_wp * ( dzeta_dphi(i,j,k) + dzeta_dphi(i+1,j,k) )

          utens(i,j,k)= utens(i,j,k) - rho_inv_at_u(i,j,k) *                         &
            &    r_earth_recip * dzeta_dphi_at_u * ( work_ik(i,k+1) - work_ik(i,k) )

        END DO
      END DO

    END DO

  END IF

  ! --- utens: contributions by tau_13 ---

  ! tau13 at uw-position (can be used also for w-tens!):
  DO k=1, ke+1
    DO j = MIN(jstartu, jstart), MAX(jendu, jend)
      DO i = MIN(istartu, istart-1), MAX(iendu, iend)

        dw_dlam_at_uw   = eddlon * ( w(i+1,j,k,nnow) - w(i,j,k,nnow) )

        ztau13_at_uw(i,j,k) = - Kvm_rho_at_uw(i,j,k) * acrlat(j,1) * dw_dlam_at_uw

      END DO
    END DO
  END DO

  IF ( l3dturb_metr ) THEN

    k=1
    DO j = MIN(jstartu, jstart), MAX(jendu, jend)
      DO i = MIN(istartu, istart-1), MAX(iendu, iend)

        dw_dzeta_at_uw  = 0.5_wp *                                 &
          &        ( ( w(i  ,j,k+1,nnow) - w(i  ,j,k,nnow) )       &
          &        + ( w(i+1,j,k+1,nnow) - w(i+1,j,k,nnow) ) )

        ztau13_at_uw(i,j,k) = ztau13_at_uw(i,j,k) - Kvm_rho_at_uw(i,j,k)     &
          &      * acrlat(j,1) * dzeta_dlam_at_uw(i,j,k) * dw_dzeta_at_uw
      END DO
    END DO

    DO k=2, ke-1
      DO j = MIN(jstartu, jstart), MAX(jendu, jend)
        DO i = MIN(istartu, istart-1), MAX(iendu, iend)

          dw_dzeta_at_uw  = 0.25_wp *                                 &
            &        ( ( w(i  ,j,k+1,nnow) - w(i  ,j,k-1,nnow) )      &
            &        + ( w(i+1,j,k+1,nnow) - w(i+1,j,k-1,nnow) ) )

          ztau13_at_uw(i,j,k) = ztau13_at_uw(i,j,k) - Kvm_rho_at_uw(i,j,k)    &
            &     * acrlat(j,1) * dzeta_dlam_at_uw(i,j,k) * dw_dzeta_at_uw
        END DO
      END DO
    END DO

    k=ke
    DO j = MIN(jstartu, jstart), MAX(jendu, jend)
      DO i = MIN(istartu, istart-1), MAX(iendu, iend)

        dw_dzeta_at_uw = 0.5_wp *                                  &
          &        ( ( w(i  ,j,k,nnow) - w(i  ,j,k-1,nnow) )       &
          &        + ( w(i+1,j,k,nnow) - w(i+1,j,k-1,nnow) ) )

        ztau13_at_uw(i,j,k) = ztau13_at_uw(i,j,k) - Kvm_rho_at_uw(i,j,k)    &
          &      * acrlat(j,1) * dzeta_dlam_at_uw(i,j,k) * dw_dzeta_at_uw
      END DO
    END DO

  END IF

  !MB   DO k=2, ke
  DO k=2, ke-1
    DO j = jstartu, jendu
      DO i = istartu, iendu

        utens(i,j,k)= utens(i,j,k) - rho_inv_at_u(i,j,k) *                   &
          &    (-sqrtg_r_u(i,j,k)) * ( ztau13_at_uw(i,j,k+1) - ztau13_at_uw(i,j,k) )

      END DO
    END DO
  END DO

  DEALLOCATE( rho_inv_at_u )


  ! The following contribution to tau13 is treated in 'implicit_vert_diffusion_uvwt_3D'.
  ! Nevertheless, it is needed later on for wtens
  k=1
  DO j = MIN(jstartu, jstart), MAX(jendu, jend)
    DO i = MIN(istartu, istart-1), MAX(iendu, iend)

      du_dzeta_at_uw = u(i,j,k+1,nnow) - u(i,j,k,nnow)       ! ??????????????

      dzeta_dz_at_uw = - 0.5_wp * ( sqrtg_r_w(i,j,k) + sqrtg_r_w(i+1,j,k) )

      ztau13_at_uw(i,j,k) = ztau13_at_uw(i,j,k)                      &
        &    - Kvm_rho_at_uw(i,j,k) * dzeta_dz_at_uw * du_dzeta_at_uw

    END DO
  END DO

  DO k=2, ke
    DO j = MIN(jstartu, jstart), MAX(jendu, jend)
      DO i = MIN(istartu, istart-1), MAX(iendu, iend)

        du_dzeta_at_uw = u(i,j,k,nnow) - u(i,j,k-1,nnow)

        dzeta_dz_at_uw = -0.5_wp * ( sqrtg_r_w(i,j,k) + sqrtg_r_w(i+1,j,k) )

        ztau13_at_uw(i,j,k) = ztau13_at_uw(i,j,k)     &
          &   - Kvm_rho_at_uw(i,j,k) *  dzeta_dz_at_uw * du_dzeta_at_uw 
      END DO
    END DO
  END DO

  k=ke+1
  DO j = MIN(jstartu, jstart), MAX(jendu, jend)
    DO i = MIN(istartu, istart-1), MAX(iendu, iend)

      du_dzeta_at_uw = u(i,j,k-1,nnow) - u(i,j,k-2,nnow)    ! ??????

      dzeta_dz_at_uw = - 0.5_wp * ( sqrtg_r_w(i,j,k) + sqrtg_r_w(i+1,j,k) )

      ztau13_at_uw(i,j,k) = ztau13_at_uw(i,j,k)     &
        &   - Kvm_rho_at_uw(i,j,k) *  dzeta_dz_at_uw * du_dzeta_at_uw 
    END DO
  END DO


  ! === vtens =========================================

  ALLOCATE( rho_inv_at_v (ie, je, ke) )

  DO k = 1, ke
    DO j = jstartv, jendv
      DO  i = istartv, iendv
        rho_inv_at_v(i,j,k) = 2.0_wp / ( rho(i,j,k) + rho(i,j+1,k) )
      END DO
    END DO
  END DO

  ! --- vtens: contributions by tau_22 ---

  DO  i = istartv, iendv

    DO k = 1, ke
      DO j = jstartv, jendv+1

        dv_dphi_at_s = eddlat * ( v(i,j,k,nnow) - v(i,j-1,k,nnow) )

        ! tau22 at s-position:
        work_jk(j,k) = -2.0_wp * Kvm_rho_at_s(i,j,k) * r_earth_recip * dv_dphi_at_s
      END DO
    END DO

    IF ( l3dturb_metr ) THEN

      k = 1
      DO j = jstartv, jendv+1

        dv_dzeta_at_s = 0.5_wp * ( ( v(i,j-1,k+1,nnow) - v(i,j-1,k,nnow) )    &
          &                      + ( v(i,j  ,k+1,nnow) - v(i,j  ,k,nnow) ) ) 

        ! tau22 at s-position:
        work_jk(j,k) = work_jk(j,k)                                           &
          &      - 2.0_wp * Kvm_rho_at_s(i,j,k)                               &
          &            * r_earth_recip * dzeta_dphi(i,j,k) * dv_dzeta_at_s
      END DO

      DO k = 2, ke-1
        DO j = jstartv, jendv+1

          dv_dzeta_at_s = 0.25_wp * ( ( v(i,j-1,k+1,nnow) - v(i,j-1,k-1,nnow) )    &
            &                       + ( v(i,j  ,k+1,nnow) - v(i,j  ,k-1,nnow) ) ) 

          ! tau22 at s-position:
          work_jk(j,k) = work_jk(j,k)                                          &
            &      - 2.0_wp * Kvm_rho_at_s(i,j,k)                              &
            &            * r_earth_recip * dzeta_dphi(i,j,k) * dv_dzeta_at_s
        END DO
      END DO

      k = ke
      DO j = jstartv, jendv+1

        dv_dzeta_at_s = 0.5_wp * ( ( v(i,j-1,k,nnow) - v(i,j-1,k-1,nnow) )    &
          &                      + ( v(i,j  ,k,nnow) - v(i,j  ,k-1,nnow) ) ) 

        ! tau22 at s-position:
        work_jk(j,k) = work_jk(j,k)                                          &
          &      - 2.0_wp * Kvm_rho_at_s(i,j,k)                              &
          &            * r_earth_recip * dzeta_dphi(i,j,k) * dv_dzeta_at_s
      END DO

    END IF

    !MB  DO k = 1, ke
    DO k = 2, ke-1
      DO j = jstartv, jendv

        vtens(i,j,k)= vtens(i,j,k) - rho_inv_at_v(i,j,k)       &
          &    * r_earth_recip * ( work_jk(j+1,k) - work_jk(j,k) ) * eddlat

      END DO
    END DO

  END DO   ! DO i


  IF ( l3dturb_metr ) THEN

    DO  i = istartv, iendv

      ! tau22 at vw-position:
      k=1
      DO j = jstartv, jendv

        dv_dphi_at_vw = 0.5_wp * ( v(i,j+1,k,nnow) - v(i,j-1,k,nnow) ) * eddlat

        dv_dzeta = v(i,j,k+1,nnow) - v(i,j,k,nnow)       ! ??????

        ! tau22 at vw-position:
        work_jk(j,k) = -2.0_wp * Kvm_rho_at_vw(i,j,k) * r_earth_recip           &
          &    * ( dv_dphi_at_vw                                                &
          &      + (1.0_wp-beta_vv1(k)) * dzeta_dphi_at_vw(i,j,k) * dv_dzeta )
      END DO

      DO k=2, ke
        DO j = jstartv, jendv

          dv_dphi_at_vw = 0.25_wp * eddlat                           &
            &       * ( ( v(i,j+1,k-1,nnow) - v(i,j-1,k-1,nnow) )    &
            &         + ( v(i,j+1,k  ,nnow) - v(i,j-1,k  ,nnow) ) )

          dv_dzeta = v(i,j,k,nnow) - v(i,j,k-1,nnow)

          ! tau22 at vw-position:
          work_jk(j,k) = -2.0_wp * Kvm_rho_at_vw(i,j,k) * r_earth_recip          &
            &    * ( dv_dphi_at_vw                                               &
            &      + (1.0_wp-beta_vv1(k)) * dzeta_dphi_at_vw(i,j,k) * dv_dzeta )
        END DO
      END DO

      !MB  DO k=1, ke-1
      DO k=2, ke-1
        DO j = jstartv, jendv

          vtens(i,j,k)= vtens(i,j,k) - rho_inv_at_v(i,j,k)                      &
            &  * r_earth_recip * dzeta_dphi_at_v(i,j,k) * ( work_jk(j,k+1) - work_jk(j,k) )

        END DO
      END DO

    END DO    ! DO i

  END IF

  ! --- vtens: contributions by tau_12 ---

  ! (use tau12 at uv-position from above)

  !MB   DO k = 1, ke
  DO k = 2, ke-1
    DO j = jstartv, jendv
      DO i = istartv, iendv

        vtens(i,j,k)= vtens(i,j,k) - rho_inv_at_v(i,j,k) *                   &
          &    acrlat(j,2) * ( ztau12_at_uv(i,j,k) - ztau12_at_uv(i-1,j,k) ) * eddlon

      END DO
    END DO
  END DO

  ! tau12 at vw-position:

  IF ( l3dturb_metr ) THEN

    DO i = istartv, iendv

      DO k=2, ke
        DO j = jstartv, jendv

          du_dphi_at_vw = 0.25_wp * eddlat *                       &
            &    ( ( u(i  ,j+1,k  ,nnow) - u(i  ,j,k  ,nnow) )     &
            &    + ( u(i-1,j+1,k  ,nnow) - u(i-1,j,k  ,nnow) )     &
            &    + ( u(i  ,j+1,k-1,nnow) - u(i  ,j,k-1,nnow) )     &
            &    + ( u(i-1,j+1,k-1,nnow) - u(i-1,j,k-1,nnow) ) )

          dv_dlam_at_vw = 0.25_wp * eddlon *                       &
            &    ( ( v(i+1,j,k  ,nnow) - v(i-1,j,k  ,nnow) )       &
            &    + ( v(i+1,j,k-1,nnow) - v(i-1,j,k-1,nnow) ) )

          du_dzeta_at_vw = 0.25_wp * ( ( u(i  ,j  ,k,nnow) - u(i  ,j  ,k-1,nnow) )    &
            &                        + ( u(i  ,j+1,k,nnow) - u(i  ,j+1,k-1,nnow) )    &
            &                        + ( u(i-1,j  ,k,nnow) - u(i-1,j  ,k-1,nnow) )    &
            &                        + ( u(i-1,j+1,k,nnow) - u(i-1,j+1,k-1,nnow) ) )

          dv_dzeta_at_vw = v(i,j,k,nnow) - v(i,j,k-1,nnow)

          dzeta_dlam_at_vw = 0.5_wp * ( dzeta_dlam_at_w(i,j,k) + dzeta_dlam_at_w(i,j+1,k) )

          ! tau12 at vw-position:
          work_jk(j,k) = - Kvm_rho_at_vw(i,j,k)                                                 &
            &   * ( acrlat(j,2) * ( dv_dlam_at_vw                                               &
            &                     + (1.0_wp-beta_vv3(k)) * dzeta_dlam_at_vw  * dv_dzeta_at_vw ) &
            &     + r_earth_recip * ( du_dphi_at_vw                                             &
            &                     + dzeta_dphi_at_vw(i,j,k) * du_dzeta_at_vw ) )


        END DO
      END DO

      !MB   DO k=2, ke
      DO k=2, ke-1
        DO j = jstartv, jendv

          dzeta_dlam_at_v = 0.5_wp * ( dzeta_dlam(i,j,k) + dzeta_dlam(i,j+1,k) )

          vtens(i,j,k)= vtens(i,j,k) - rho_inv_at_v(i,j,k) *                         &
            &    acrlat(j,2) * dzeta_dlam_at_v * ( work_jk(j,k+1) - work_jk(j,k) )

        END DO
      END DO

    END DO

  END IF

  ! --- vtens: contributions by tau_32 ---

  ! tau32 at vw-position (can be used also for w-tens!):
  DO k=1, ke+1
    DO j = MIN(jstartv, jstart-1), MAX(jendv, jend)
      DO i = MIN(istartv, istart), MAX(iendv, iend)

        dw_dphi_at_vw   = eddlat * ( w(i,j+1,k,nnow) - w(i,j,k,nnow) )

        ztau23_at_vw(i,j,k) = - Kvm_rho_at_vw(i,j,k) * r_earth_recip * dw_dphi_at_vw

      END DO
    END DO
  END DO

  IF ( l3dturb_metr ) THEN

    k=1
    DO j = MIN(jstartv, jstart-1), MAX(jendv, jend)
      DO i = MIN(istartv, istart), MAX(iendv, iend)

        dw_dzeta_at_vw  = 0.5_wp *                                 &
          &        ( ( w(i,j  ,k+1,nnow) - w(i,j  ,k,nnow) )       &
          &        + ( w(i,j+1,k+1,nnow) - w(i,j+1,k,nnow) ) )

        ztau23_at_vw(i,j,k) = ztau23_at_vw(i,j,k) - Kvm_rho_at_vw(i,j,k) *      &
          &          r_earth_recip * dzeta_dphi_at_vw(i,j,k) * dw_dzeta_at_vw
      END DO
    END DO

    DO k=2, ke-1
      DO j = MIN(jstartv, jstart-1), MAX(jendv, jend)
        DO i = MIN(istartv, istart), MAX(iendv, iend)

          dw_dzeta_at_vw  = 0.25_wp *                                 &
            &        ( ( w(i,j  ,k+1,nnow) - w(i,j  ,k-1,nnow) )      &
            &        + ( w(i,j+1,k+1,nnow) - w(i,j+1,k-1,nnow) ) )

          ztau23_at_vw(i,j,k) = ztau23_at_vw(i,j,k) - Kvm_rho_at_vw(i,j,k)       &
            &     * r_earth_recip * dzeta_dphi_at_vw(i,j,k) * dw_dzeta_at_vw
        END DO
      END DO
    END DO

    k=ke
    DO j = MIN(jstartv, jstart-1), MAX(jendv, jend)
      DO i = MIN(istartv, istart), MAX(iendv, iend)

        dw_dzeta_at_vw = 0.5_wp *                                  &
          &        ( ( w(i,j  ,k,nnow) - w(i,j  ,k-1,nnow) )       &
          &        + ( w(i,j+1,k,nnow) - w(i,j+1,k-1,nnow) ) )

        ztau23_at_vw(i,j,k) = ztau23_at_vw(i,j,k) - Kvm_rho_at_vw(i,j,k)      &
          &     * r_earth_recip * dzeta_dphi_at_vw(i,j,k) * dw_dzeta_at_vw
      END DO
    END DO

  END IF

  !MB   DO k=2, ke
  DO k=2, ke-1
    DO j = jstartv, jendv
      DO i = istartv, iendv

        vtens(i,j,k)= vtens(i,j,k) - rho_inv_at_v(i,j,k) *                   &
          &    (-sqrtg_r_v(i,j,k)) * ( ztau23_at_vw(i,j,k+1) - ztau23_at_vw(i,j,k) )

      END DO
    END DO
  END DO

  DEALLOCATE( rho_inv_at_v )


  ! The following contribution to tau23 is treated in 'implicit_vert_diffusion_uvwt_3D'.
  ! Nevertheless, it is needed later on for wtens
  k=1
  DO j = MIN(jstartv, jstart-1), MAX(jendv, jend)
    DO i = MIN(istartv, istart), MAX(iendv, iend)

      dv_dzeta_at_vw  = v(i,j,k+1,nnow) - v(i,j,k,nnow)       ! ??????????????

      dzeta_dz_at_vw = - 0.5_wp * ( sqrtg_r_w(i,j,k) + sqrtg_r_w(i,j+1,k) )

      ztau23_at_vw(i,j,k) = ztau23_at_vw(i,j,k)       &
        &     - Kvm_rho_at_vw(i,j,k) * dzeta_dz_at_vw * dv_dzeta_at_vw

    END DO
  END DO

  DO k=2, ke
    DO j = MIN(jstartv, jstart-1), MAX(jendv, jend)
      DO i = MIN(istartv, istart), MAX(iendv, iend)

        dv_dzeta_at_vw  = v(i,j,k,nnow) - v(i,j,k-1,nnow)

        dzeta_dz_at_vw = - 0.5_wp * ( sqrtg_r_w(i,j,k) + sqrtg_r_w(i,j+1,k) )

        ztau23_at_vw(i,j,k) = ztau23_at_vw(i,j,k)       &
          &    - Kvm_rho_at_vw(i,j,k) * dzeta_dz_at_vw * dv_dzeta_at_vw
      END DO
    END DO
  END DO

  k=ke+1
  DO j = MIN(jstartv, jstart-1), MAX(jendv, jend)
    DO i = MIN(istartv, istart), MAX(iendv, iend)

      dv_dzeta_at_vw  = v(i,j,k-1,nnow) - v(i,j,k-2,nnow)       ! ????????

      dzeta_dz_at_vw = - 0.5_wp * ( sqrtg_r_w(i,j,k) + sqrtg_r_w(i,j+1,k) )

      ztau23_at_vw(i,j,k) = ztau23_at_vw(i,j,k)      &
        &    - Kvm_rho_at_vw(i,j,k) * dzeta_dz_at_vw * dv_dzeta_at_vw 

    END DO
  END DO


  ! === wtens =========================================

  ALLOCATE( rho_inv_at_w (ie, je, ke+1) )

  DO k = 2, ke
    DO j = jstart, jend
      DO  i = istart, iend
        rho_inv_at_w(i,j,k) = 2.0_wp / ( rho(i,j,k) + rho(i,j,k-1) )
      END DO
    END DO
  END DO
  rho_inv_at_w(:,:,1   ) = 1.0_wp / rho(:,:,1)
  rho_inv_at_w(:,:,ke+1) = 1.0_wp / rho(:,:,ke)

  ! --- wtens: contributions by tau_13 ---

  !MB    DO k = 2, ke
  DO k = 3, ke-1
    DO j = jstart, jend
      DO  i = istart, iend

        wtens(i,j,k)= wtens(i,j,k) - rho_inv_at_w(i,j,k)                      &
          &    * acrlat(j,1) * ( ztau13_at_uw(i,j,k) - ztau13_at_uw(i-1,j,k) ) * eddlon
      END DO
    END DO
  END DO

  IF ( l3dturb_metr ) THEN

    DO j = jstart, jend

      k=1
      DO i = istart, iend

        du_dzeta_at_s  = 0.5_wp * ( ( u(i  ,j,k+1,nnow) - u(i  ,j,k,nnow) )   &
          &                       + ( u(i-1,j,k+1,nnow) - u(i-1,j,k,nnow) ) )

        dw_dlam_at_s   = 0.25_wp * eddlon *                              &
          &               ( ( w(i+1,j,k  ,nnow) - w(i-1,j,k  ,nnow) )    &
          &               + ( w(i+1,j,k+1,nnow) - w(i-1,j,k+1,nnow) ) )

        dw_dzeta_at_s  = w(i,j,k+1,nnow) - w(i,j,k,nnow)

        ! tau13 at scalar position
        work_ik(i,k) = - Kvm_rho_at_s(i,j,k) *              &
          &       ( - sqrtg_r_s(i,j,k) * du_dzeta_at_s      &
          &         + acrlat(j,1) * ( dw_dlam_at_s          &
          &                         + (1.0_wp-beta_ww1(k)) * dzeta_dlam(i,j,k) * dw_dzeta_at_s ) )
      END DO

      DO k=2, ke-1
        DO i = istart, iend

          du_dzeta_at_s  = 0.25_wp * ( ( u(i  ,j,k+1,nnow) - u(i  ,j,k-1,nnow) )    &
            &                       + ( u(i-1,j,k+1,nnow) - u(i-1,j,k-1,nnow) ) )

          dw_dlam_at_s   = 0.25_wp * eddlon *                              &
            &               ( ( w(i+1,j,k  ,nnow) - w(i-1,j,k  ,nnow) )    &
            &               + ( w(i+1,j,k+1,nnow) - w(i-1,j,k+1,nnow) ) )

          dw_dzeta_at_s  = w(i,j,k+1,nnow) - w(i,j,k,nnow)

          ! tau13 at scalar position
          work_ik(i,k) = - Kvm_rho_at_s(i,j,k) *              &
            &       ( - sqrtg_r_s(i,j,k) * du_dzeta_at_s      &
            &         + acrlat(j,1) * ( dw_dlam_at_s          &
            &                         + (1.0_wp-beta_ww1(k)) * dzeta_dlam(i,j,k) * dw_dzeta_at_s ) )
        END DO
      END DO

      k=ke
      DO i = istart, iend

        du_dzeta_at_s  = 0.5_wp * ( ( u(i  ,j,k,nnow) - u(i  ,j,k-1,nnow) )    &
          &                       + ( u(i-1,j,k,nnow) - u(i-1,j,k-1,nnow) ) )

        dw_dlam_at_s   = 0.25_wp * eddlon *                              &
          &               ( ( w(i+1,j,k  ,nnow) - w(i-1,j,k  ,nnow) )    &
          &               + ( w(i+1,j,k+1,nnow) - w(i-1,j,k+1,nnow) ) )

        dw_dzeta_at_s  = w(i,j,k+1,nnow) - w(i,j,k,nnow)

        ! tau13 at scalar position
        work_ik(i,k) = - Kvm_rho_at_s(i,j,k) *              &
          &       ( - sqrtg_r_s(i,j,k) * du_dzeta_at_s      &
          &         + acrlat(j,1) * ( dw_dlam_at_s          &
          &                         + (1.0_wp-beta_ww1(k)) * dzeta_dlam(i,j,k) * dw_dzeta_at_s ) )
      END DO

      !MB   DO k = 2, ke
      DO k = 3, ke-1
        DO  i = istart, iend

          wtens(i,j,k)= wtens(i,j,k) - rho_inv_at_w(i,j,k)           &
            &   * acrlat(j,1) * dzeta_dlam_at_w(i,j,k) * ( work_ik(i,k) - work_ik(i,k-1) )

        END DO
      END DO

    END DO     ! DO j

  END IF

  ! --- wtens: contributions by tau_23 ---

  !MB   DO k = 2, ke
  DO k = 3, ke-1
    DO j = jstart, jend
      DO  i = istart, iend

        wtens(i,j,k)= wtens(i,j,k) - rho_inv_at_w(i,j,k)                        &
          &    * r_earth_recip * ( ztau23_at_vw(i,j,k) - ztau23_at_vw(i,j-1,k) ) * eddlat
      END DO
    END DO
  END DO

  IF ( l3dturb_metr ) THEN

    DO j = jstart, jend

      k=1
      DO i = istart, iend

        dv_dzeta_at_s  = 0.5_wp * ( ( v(i,j  ,k+1,nnow) - v(i,j  ,k,nnow) )     &
          &                       + ( v(i,j-1,k+1,nnow) - v(i,j-1,k,nnow) ) )

        dw_dphi_at_s   = 0.25_wp * eddlat *                              &
          &               ( ( w(i,j+1,k  ,nnow) - w(i,j-1,k  ,nnow) )    &
          &               + ( w(i,j+1,k+1,nnow) - w(i,j-1,k+1,nnow) ) )

        dw_dzeta_at_s  = w(i,j,k+1,nnow) - w(i,j,k,nnow)

        ! tau23 at scalar position
        work_ik(i,k) = - Kvm_rho_at_s(i,j,k) *                &
          &       ( - sqrtg_r_s(i,j,k) * dv_dzeta_at_s        &
          &         + r_earth_recip * ( dw_dphi_at_s          &
          &                           + (1.0_wp-beta_ww1(k)) * dzeta_dphi(i,j,k) * dw_dzeta_at_s ) )
      END DO

      DO k=2, ke-1
        DO i = istart, iend

          dv_dzeta_at_s  = 0.25_wp * ( ( v(i,j  ,k+1,nnow) - v(i,j  ,k-1,nnow) )   &
            &                        + ( v(i,j-1,k+1,nnow) - v(i,j-1,k-1,nnow) ) )

          dw_dphi_at_s   = 0.25_wp * eddlat *                              &
            &               ( ( w(i,j+1,k  ,nnow) - w(i,j-1,k  ,nnow) )    &
            &               + ( w(i,j+1,k+1,nnow) - w(i,j-1,k+1,nnow) ) )

          dw_dzeta_at_s  = w(i,j,k+1,nnow) - w(i,j,k,nnow)

          ! tau23 at scalar position
          work_ik(i,k) = - Kvm_rho_at_s(i,j,k) *                &
            &       ( - sqrtg_r_s(i,j,k) * dv_dzeta_at_s        &
            &         + r_earth_recip * ( dw_dphi_at_s          &
            &                           + (1.0_wp-beta_ww1(k)) * dzeta_dphi(i,j,k) * dw_dzeta_at_s ) )
        END DO
      END DO

      k=ke
      DO i = istart, iend

        dv_dzeta_at_s  = 0.5_wp * ( ( v(i,j  ,k,nnow) - v(i,j  ,k-1,nnow) )   &
          &                       + ( v(i,j-1,k,nnow) - v(i,j-1,k-1,nnow) ) )

        dw_dphi_at_s   = 0.25_wp * eddlat *                              &
          &               ( ( w(i,j+1,k  ,nnow) - w(i,j-1,k  ,nnow) )    &
          &               + ( w(i,j+1,k+1,nnow) - w(i,j-1,k+1,nnow) ) )

        dw_dzeta_at_s  = w(i,j,k+1,nnow) - w(i,j,k,nnow)

        ! tau23 at scalar position
        work_ik(i,k) = - Kvm_rho_at_s(i,j,k) *                &
          &       ( - sqrtg_r_s(i,j,k) * dv_dzeta_at_s        &
          &         + r_earth_recip * ( dw_dphi_at_s          &
          &                           + dzeta_dphi(i,j,k) * dw_dzeta_at_s ) )
      END DO

      !MB  DO k = 2, ke
      DO k = 3, ke-1
        DO  i = istart, iend

          wtens(i,j,k)= wtens(i,j,k) - rho_inv_at_w(i,j,k)           &
            &   * r_earth_recip * dzeta_dphi_at_w(i,j,k) * ( work_ik(i,k) - work_ik(i,k-1) )

        END DO
      END DO

    END DO   ! DO j

  END IF


  ! --- wtens: contributions by tau_33 ---

  ! the following tendency is already treated in 'implicit_vert_diffusion_uvwt_3D'

  !DO j = jstart, jend

  !  DO k=1, ke-1
  !    DO i = istart, iend

  !      dw_dzeta_at_s  = w(i,j,k+1,nnow) - w(i,j,k,nnow)

  !      ! tau33 at 2-position:
  !      work_ik(i,k) = - 2.0_wp * Kvm_rho_at_s(i,j,k) *               &
  !         &       (1.0_wp-beta_ww2(k)) * (-sqrtg_r_s(i,j,k)) * dw_dzeta_at_s
  !    END DO
  !  END DO

  !  k=ke
  !  DO i = istart, iend

  !    dw_dzeta_at_s  = w(i,j,k,nnow) - w(i,j,k-1,nnow)

  !    ! tau33 at 2-position:
  !    work_ik(i,k) = - 2.0_wp * Kvm_rho_at_s(i,j,k) *               &
  !       &       (1.0_wp-beta_ww2(k)) * (-sqrtg_r_s(i,j,k)) * dw_dzeta_at_s
  !  END DO

  !  DO k = 2, ke
  !    DO  i = istart, iend

  !      wtens(i,j,k)= wtens(i,j,k) - rho_inv_at_w(i,j,k) *          &
  !         &   (-sqrtg_r_w(i,j,k)) * ( work_ik(i,k) - work_ik(i,k-1) )

  !    END DO
  !  END DO

  !END DO


  ! ========== scalar diffusion: T ========================


  ALLOCATE( rho_inv (ie, je, ke) )

  DO k = 1, ke
    DO j = jstartu-1, jendv+1
      DO  i = istartu-1, iendu+1
        rho_inv(i,j,k) = 1.0_wp / rho(i,j,k)
      END DO
    END DO
  END DO


  ! === ttens =========================================

  ALLOCATE( Kvh_rho_pi_at_s (ie, je, ke) )
  ALLOCATE( Kvh_rho_pi_at_u (ie, je, ke) )
  ALLOCATE( Kvh_rho_pi_at_v (ie, je, ke) )

  ! determine diffusion coefficients at different positions:

  ALLOCATE( work(ie, je, ke) )

  DO k=1, ke
    DO j = jstart-2, jend+2
      DO  i = istart-2, iend+2
        z_pi = ( ( p0(i,j,k) + pp(i,j,k,nnow) ) * p0ref_recip )**rdocp
        work(i,j,k) = rho(i,j,k) * z_pi
      END DO
    END DO
  END DO

  DO k=1, ke
    ! for averaging from the w-position to the scalar position:
    !k_w2s_up  = MAX( 1,  k   )
    !k_w2s_low = MIN( ke, k+1 )

    ! for averaging from the scalar position to the w-position
    k_s2w_up  = MAX( k-1, 1 )
    k_s2w_low = MIN( k  , ke )

    DO j = jstart-1, jend+1
      DO  i = istart-1, iend+1

        ! Kvh * rho * pi at the two grid positions 'scalar' and 'w' need the least averaging:
        Kvh_rho_pi_at_s (i,j,k) = work(i,j,k) * 0.5_wp * ( tkvh(i,j,k) + tkvh(i,j,k+1) ) 

      END DO
    END DO
  END DO

  DO k=1, ke
    DO j = jstart, jend
      DO  i = istart-1, iend
        Kvh_rho_pi_at_u(i,j,k) = 0.5_wp * ( Kvh_rho_pi_at_s(i,j,k) + Kvh_rho_pi_at_s(i+1,j,k) ) 
      END DO
    END DO

    DO j = jstart-1, jend
      DO  i = istart, iend
        Kvh_rho_pi_at_v(i,j,k) = 0.5_wp * ( Kvh_rho_pi_at_s(i,j,k) + Kvh_rho_pi_at_s(i,j+1,k) ) 
      END DO
    END DO
  END DO

  DEALLOCATE( work )

  ! --- ttens: contributions by H1 ---

  DO j = jstart, jend

    DO k = 1, ke
      DO  i = istart-1, iend

        ds_dlam_at_u = eddlon * ( ztheta(i+1,j,k) - ztheta(i,j,k) )

        ! H1 at u-position:
        work_ik(i,k) = - Kvh_rho_pi_at_u(i,j,k) * acrlat(j,1) * ds_dlam_at_u
      END DO
    END DO

    IF ( l3dturb_metr ) THEN

      k = 1
      DO  i = istart-1, iend

        ds_dzeta_at_u = 0.5_wp * ( ( ztheta(i+1,j,k+1) - ztheta(i+1,j,k) )    &
          &                      + ( ztheta(i  ,j,k+1) - ztheta(i  ,j,k) ) ) 

        ! H1 at u-position:
        work_ik(i,k) = work_ik(i,k)   - Kvh_rho_pi_at_u(i,j,k)                  &
          &            * acrlat(j,1) * dzeta_dlam_at_u(i,j,k) * ds_dzeta_at_u
      END DO

      DO k = 2, ke-1
        DO  i = istart-1, iend

          ds_dzeta_at_u = 0.25_wp * ( ( ztheta(i+1,j,k+1) - ztheta(i+1,j,k-1) )    &
            &                       + ( ztheta(i  ,j,k+1) - ztheta(i  ,j,k-1) ) ) 

          ! H1 at u-position:
          work_ik(i,k) = work_ik(i,k)   - Kvh_rho_pi_at_u(i,j,k)                 &
            &            * acrlat(j,1) * dzeta_dlam_at_u(i,j,k) * ds_dzeta_at_u

        END DO
      END DO

      k = ke
      DO  i = istart-1, iend

        ds_dzeta_at_u = 0.5_wp * ( ( ztheta(i+1,j,k) - ztheta(i+1,j,k-1) )    &
          &                      + ( ztheta(i  ,j,k) - ztheta(i  ,j,k-1) ) ) 

        ! H1 at u-position:
        work_ik(i,k) = work_ik(i,k) - Kvh_rho_pi_at_u(i,j,k)                   &
          &            * acrlat(j,1) * dzeta_dlam_at_u(i,j,k) * ds_dzeta_at_u
      END DO

    END IF

    !MB  DO k = 1, ke
    DO k = 2, ke-1
      DO i = istart, iend

        ttens(i,j,k)= ttens(i,j,k) - rho_inv(i,j,k)     &
          &    * acrlat(j,1) * ( work_ik(i,k) - work_ik(i-1,k) ) * eddlon

      END DO
    END DO

  END DO


  IF ( l3dturb_metr ) THEN

    DO j = jstart, jend

      ! H1 at w-position:
      k=1
      DO  i = istart, iend

        ds_dlam_at_w = 0.5_wp * ( ztheta(i+1,j,k) - ztheta(i-1,j,k) ) * eddlon

        ds_dzeta = ztheta(i,j,k+1) - ztheta(i,j,k)     ! ????????

        ! H1 at w-position:
        work_ik(i,k) = - Kvh_rho_pi_at_w(i,j,k) * acrlat(j,1)              &
          &    * ( ds_dlam_at_w                                            &
          &      + (1.0_wp-beta_1(k)) * dzeta_dlam_at_w(i,j,k) * ds_dzeta )
      END DO

      DO k=2, ke
        DO  i = istart, iend

          ds_dlam_at_w = 0.25_wp * ( ( ztheta(i+1,j,k-1) - ztheta(i-1,j,k-1) )    &
            &                      + ( ztheta(i+1,j,k  ) - ztheta(i-1,j,k  ) ) ) * eddlon

          ds_dzeta = ztheta(i,j,k) - ztheta(i,j,k-1)

          ! H1 at w-position:
          work_ik(i,k) = - Kvh_rho_pi_at_w(i,j,k) * acrlat(j,1)              &
            &    * ( ds_dlam_at_w                                            &
            &      + (1.0_wp-beta_1(k)) * dzeta_dlam_at_w(i,j,k) * ds_dzeta )
        END DO
      END DO

      !MB   DO k=1, ke-1
      DO k=2, ke-1
        DO i = istart, iend

          ttens(i,j,k)= ttens(i,j,k) - rho_inv(i,j,k)                      &
            &  * acrlat(j,1) * dzeta_dlam(i,j,k) * ( work_ik(i,k+1) - work_ik(i,k) )

        END DO
      END DO

    END DO

  END IF


  ! --- ttens: contributions by H2 ---

  DO  i = istart, iend

    DO k = 1, ke
      DO j = jstart-1, jend

        ds_dphi_at_v = eddlat * ( ztheta(i,j+1,k) - ztheta(i,j,k) )

        ! H2 at v-position:
        work_jk(j,k) = - Kvh_rho_pi_at_v(i,j,k) * r_earth_recip * ds_dphi_at_v
      END DO
    END DO

    IF ( l3dturb_metr ) THEN

      k = 1
      DO j = jstart-1, jend

        ds_dzeta_at_v = 0.5_wp * ( ( ztheta(i,j+1,k+1) - ztheta(i,j+1,k) )    &
          &                      + ( ztheta(i,j  ,k+1) - ztheta(i,j  ,k) ) ) 

        ! H2 at v-position:
        work_jk(j,k) = work_jk(j,k)  - Kvh_rho_pi_at_v(i,j,k)                      &
          &            * r_earth_recip * dzeta_dphi_at_v(i,j,k) * ds_dzeta_at_v
      END DO

      DO k = 2, ke-1
        DO j = jstart-1, jend

          ds_dzeta_at_v = 0.25_wp * ( ( ztheta(i,j+1,k+1) - ztheta(i,j+1,k-1) )    &
            &                       + ( ztheta(i,j  ,k+1) - ztheta(i,j  ,k-1) ) ) 

          ! H2 at v-position:
          work_jk(j,k) = work_jk(j,k)  - Kvh_rho_pi_at_v(i,j,k)                     &
            &            * r_earth_recip * dzeta_dphi_at_v(i,j,k) * ds_dzeta_at_v
        END DO
      END DO

      k = ke
      DO j = jstart-1, jend

        ds_dzeta_at_v = 0.5_wp * ( ( ztheta(i,j+1,k) - ztheta(i,j+1,k-1) )    &
          &                      + ( ztheta(i,j  ,k) - ztheta(i,j  ,k-1) ) ) 

        ! H2 at v-position:
        work_jk(j,k) = work_jk(j,k)  - Kvh_rho_pi_at_v(i,j,k)                     &
          &            * r_earth_recip * dzeta_dphi_at_v(i,j,k) * ds_dzeta_at_v
      END DO

    END IF

    !MB   DO k = 1, ke
    DO k = 2, ke-1
      DO j = jstart, jend

        ttens(i,j,k)= ttens(i,j,k) - rho_inv(i,j,k)     &
          &    * r_earth_recip * ( work_jk(j,k) - work_jk(j-1,k) ) * eddlat

      END DO
    END DO

  END DO


  IF ( l3dturb_metr ) THEN

    DO  i = istart, iend

      ! H2 at w-position:
      k=1
      DO j = jstart, jend

        ds_dphi_at_w = 0.5_wp * ( ztheta(i,j+1,k) - ztheta(i-1,j,k) ) * eddlat

        ds_dzeta = ztheta(i,j,k+1) - ztheta(i,j,k)     ! ????????

        ! H2 at w-position:
        work_jk(j,k) = - Kvh_rho_pi_at_w(i,j,k) * r_earth_recip            &
          &    * (                                         ds_dphi_at_w    &
          &      + (1.0_wp-beta_1(k)) * dzeta_dphi_at_w(i,j,k) * ds_dzeta )
      END DO

      DO k=2, ke
        DO j = jstart, jend

          ds_dphi_at_w = 0.25_wp * ( ( ztheta(i,j+1,k-1) - ztheta(i,j-1,k-1) )    &
            &                      + ( ztheta(i,j+1,k  ) - ztheta(i,j-1,k  ) ) ) * eddlat

          ds_dzeta = ztheta(i,j,k) - ztheta(i,j,k-1)

          ! H2 at w-position:
          work_jk(j,k) = - Kvh_rho_pi_at_w(i,j,k) * r_earth_recip             &
            &    * ( ds_dphi_at_w                                             &
            &      + (1.0_wp-beta_1(k)) * dzeta_dphi_at_w(i,j,k) * ds_dzeta )
        END DO
      END DO

      !MB   DO k=1, ke-1
      DO k=2, ke-1
        DO j = jstart, jend

          ttens(i,j,k)= ttens(i,j,k) - rho_inv(i,j,k)                      &
            &  * r_earth_recip * dzeta_dphi(i,j,k) * ( work_jk(j,k+1) - work_jk(j,k) )

        END DO
      END DO

    END DO

  END IF

  ! --- ttens: contributions by H3 ---

  ! the following tendency is already treated in 'implicit_vert_diffusion_uvwt_3D'

  !DO j = jstart, jend

  !  DO k=1, ke-1
  !    DO i = istart, iend

  !      ds_dzeta_at_w  = ztheta(i,j,k) - ztheta(i,j,k-1)

  !      ! H3 at w-position:
  !      work_ik(i,k) = - Kvh_rho_pi_at_w(i,j,k) *               &
  !         &       (1.0_wp-beta_2(k)) * (-sqrtg_r_w(i,j,k)) * ds_dzeta_at_w
  !    END DO
  !  END DO

  !  k=ke
  !  DO i = istart, iend

  !    ds_dzeta_at_w  = ztheta(i,j,k) - ztheta(i,j,k-1)

  !    ! H3 at w-position:
  !    work_ik(i,k) = - Kvh_rho_pi_at_w(i,j,k) *               &
  !       &       (1.0_wp-beta_2(k)) * (-sqrtg_r_w(i,j,k)) * ds_dzeta_at_w
  !  END DO

  !  DO k = 1, ke-1
  !    DO  i = istart, iend
  !
  !      ttens(i,j,k)= ttens(i,j,k) - rho_inv(i,j,k) *          &
  !         &   (-sqrtg_r_s(i,j,k)) * ( work_ik(i,k+1) - work_ik(i,k) )

  !    END DO
  !  END DO

  !END DO

  DEALLOCATE( Kvh_rho_pi_at_s )
  DEALLOCATE( Kvh_rho_pi_at_u )
  DEALLOCATE( Kvh_rho_pi_at_v )


  ! ========== scalar diffusion: tracers =============================


  ! loop over tracers
  DO  iztrcr = 1, trcr_get_ntrcr()

    ! check for each tracer if 3D turbulent mixing is required
    IF ( izturb(iztrcr) == T_TURB_3D ) THEN
      ! if so, increase the number of tracers undergoing 3D turb. mix.
      nztrcr_3dturb = nztrcr_3dturb + 1_iintegers
      iztrcr_3dturb(icount) = iztrcr
      icount = icount + 1
    END IF
  END DO


  IF ( nztrcr_3dturb > 0 ) THEN

    ALLOCATE( Kvh_rho_at_u (ie, je, ke) )
    ALLOCATE( Kvh_rho_at_v (ie, je, ke) )

    DO k=1, ke
      DO j = jstart, jend
        DO  i = istart-1, iend
          Kvh_rho_at_u(i,j,k) = 0.5_wp * ( Kvh_rho_at_s(i,j,k) + Kvh_rho_at_s(i+1,j,k) ) 
        END DO
      END DO

      DO j = jstart-1, jend
        DO  i = istart, iend
          Kvh_rho_at_v(i,j,k) = 0.5_wp * ( Kvh_rho_at_s(i,j,k) + Kvh_rho_at_s(i,j+1,k) ) 
        END DO
      END DO
    END DO

  END IF

  DO iztrcr = 1, nztrcr_3dturb

    CALL trcr_get( izerror, iztrcr_3dturb(iztrcr), ptr=ztrcr_now, ptr_tlev=nnow, ptr_tens = ztrcr_tens )

    ! --- ztrcr_tens: contributions by H1 ---

    IF ( my_cart_id==0 ) THEN
      CALL trcr_meta_get( izerror, iztrcr_3dturb(iztrcr), T_NAME_ID, y_tracer_name )
      WRITE(*,*) "control output (only done for step 0):"
      WRITE(*,*) "  explicit_horizontal_diffusion for tracer ", y_tracer_name
    END IF

    DO j = jstart, jend

      DO k = 1, ke
        DO  i = istart-1, iend

          ds_dlam_at_u = eddlon * ( ztrcr_now(i+1,j,k) - ztrcr_now(i,j,k) )

          ! H1 at u-position:
          work_ik(i,k) = - Kvh_rho_at_u(i,j,k) * acrlat(j,1) * ds_dlam_at_u
        END DO
      END DO

      IF ( l3dturb_metr ) THEN

        k = 1
        DO  i = istart-1, iend

          ds_dzeta_at_u = 0.5_wp * ( ( ztrcr_now(i+1,j,k+1) - ztrcr_now(i+1,j,k) )    &
            &                      + ( ztrcr_now(i  ,j,k+1) - ztrcr_now(i  ,j,k) ) ) 

          ! H1 at u-position:
          work_ik(i,k) = work_ik(i,k) - Kvh_rho_at_u(i,j,k)                         &
            &            * acrlat(j,1) * dzeta_dlam_at_u(i,j,k) * ds_dzeta_at_u
        END DO

        DO k = 2, ke-1
          DO  i = istart-1, iend

            ds_dzeta_at_u = 0.25_wp * ( ( ztrcr_now(i+1,j,k+1) - ztrcr_now(i+1,j,k-1) )    &
              &                       + ( ztrcr_now(i  ,j,k+1) - ztrcr_now(i  ,j,k-1) ) ) 

            ! H1 at u-position:
            work_ik(i,k) = work_ik(i,k) - Kvh_rho_at_u(i,j,k)                     &
              &            * acrlat(j,1) * dzeta_dlam_at_u(i,j,k) * ds_dzeta_at_u

          END DO
        END DO

        k = ke
        DO  i = istart-1, iend

          ds_dzeta_at_u = 0.5_wp * ( ( ztrcr_now(i+1,j,k) - ztrcr_now(i+1,j,k-1) )    &
            &                      + ( ztrcr_now(i  ,j,k) - ztrcr_now(i  ,j,k-1) ) ) 

          ! H1 at u-position:
          work_ik(i,k) = work_ik(i,k) - Kvh_rho_at_u(i,j,k)                       &
            &            * acrlat(j,1) * dzeta_dlam_at_u(i,j,k) * ds_dzeta_at_u
        END DO

      END IF

      !MB  DO k = 1, ke
      DO k = 2, ke-1
        DO i = istart, iend

          ztrcr_tens(i,j,k)= ztrcr_tens(i,j,k) - rho_inv(i,j,k)     &
            &    * acrlat(j,1) * ( work_ik(i,k) - work_ik(i-1,k) ) * eddlon

        END DO
      END DO

    END DO


    IF ( l3dturb_metr ) THEN

      DO j = jstart, jend

        ! H1 at w-position:
        k=1
        DO  i = istart, iend

          ds_dlam_at_w = 0.5_wp * ( ztrcr_now(i+1,j,k) - ztrcr_now(i-1,j,k) ) * eddlon

          ds_dzeta = ztrcr_now(i,j,k+1) - ztrcr_now(i,j,k)     ! ????????

          ! H1 at w-position:
          work_ik(i,k) = - Kvh_rho_at_w(i,j,k) * acrlat(j,1)                 &
            &    * ( ds_dlam_at_w                                            &
            &      + (1.0_wp-beta_1(k)) * dzeta_dlam_at_w(i,j,k) * ds_dzeta )
        END DO

        DO k=2, ke
          DO  i = istart, iend

            ds_dlam_at_w = 0.25_wp * ( ( ztrcr_now(i+1,j,k-1) - ztrcr_now(i-1,j,k-1) )    &
              &                     + ( ztrcr_now(i+1,j,k  ) - ztrcr_now(i-1,j,k  ) ) ) * eddlon

            ds_dzeta = ztrcr_now(i,j,k) - ztrcr_now(i,j,k-1)

            ! H1 at w-position:
            work_ik(i,k) = - Kvh_rho_at_w(i,j,k) * acrlat(j,1)                 &
              &    * ( ds_dlam_at_w                                            &
              &      + (1.0_wp-beta_1(k)) * dzeta_dlam_at_w(i,j,k) * ds_dzeta )
          END DO
        END DO

        !MB   DO k=1, ke-1
        DO k=2, ke-1
          DO i = istart, iend

            ztrcr_tens(i,j,k)= ztrcr_tens(i,j,k) - rho_inv(i,j,k)                      &
              &  * acrlat(j,1) * dzeta_dlam(i,j,k) * ( work_ik(i,k+1) - work_ik(i,k) )

          END DO
        END DO

      END DO

    END IF


    ! --- ztrcr_tens: contributions by H2 ---

    DO  i = istart, iend

      DO k = 1, ke
        DO j = jstart-1, jend

          ds_dphi_at_v = eddlat * ( ztrcr_now(i,j+1,k) - ztrcr_now(i,j,k) )

          ! H2 at v-position:
          work_jk(j,k) = - Kvh_rho_at_v(i,j,k) * r_earth_recip * ds_dphi_at_v
        END DO
      END DO

      IF ( l3dturb_metr ) THEN

        k = 1
        DO j = jstart-1, jend

          ds_dzeta_at_v = 0.5_wp * ( ( ztrcr_now(i,j+1,k+1) - ztrcr_now(i,j+1,k) )   &
            &                      + ( ztrcr_now(i,j  ,k+1) - ztrcr_now(i,j  ,k) ) ) 

          ! H2 at v-position:
          work_jk(j,k) = work_jk(j,k)  - Kvh_rho_at_v(i,j,k)                         &
            &            * r_earth_recip * dzeta_dphi_at_v(i,j,k) * ds_dzeta_at_v
        END DO

        DO k = 2, ke-1
          DO j = jstart-1, jend

            ds_dzeta_at_v = 0.25_wp * ( ( ztrcr_now(i,j+1,k+1) - ztrcr_now(i,j+1,k-1) )   &
              &                       + ( ztrcr_now(i,j  ,k+1) - ztrcr_now(i,j  ,k-1) ) ) 

            ! H2 at v-position:
            work_jk(j,k) = work_jk(j,k)  - Kvh_rho_at_v(i,j,k)                       &
              &            * r_earth_recip * dzeta_dphi_at_v(i,j,k) * ds_dzeta_at_v
          END DO
        END DO

        k = ke
        DO j = jstart-1, jend

          ds_dzeta_at_v = 0.5_wp * ( ( ztrcr_now(i,j+1,k) - ztrcr_now(i,j+1,k-1) )    &
            &                      + ( ztrcr_now(i,j  ,k) - ztrcr_now(i,j  ,k-1) ) ) 

          ! H2 at v-position:
          work_jk(j,k) = work_jk(j,k)  - Kvh_rho_at_v(i,j,k)                  &
            &            * r_earth_recip * dzeta_dphi_at_v(i,j,k) * ds_dzeta_at_v
        END DO

      END IF

      !MB   DO k = 1, ke
      DO k = 2, ke-1
        DO j = jstart, jend

          ztrcr_tens(i,j,k)= ztrcr_tens(i,j,k) - rho_inv(i,j,k)     &
            &    * r_earth_recip * ( work_jk(j,k) - work_jk(j-1,k) ) * eddlat

        END DO
      END DO

    END DO


    IF ( l3dturb_metr ) THEN

      DO  i = istart, iend

        ! H2 at w-position:
        k=1
        DO j = jstart, jend

          ds_dphi_at_w = 0.5_wp * ( ztrcr_now(i,j+1,k) - ztrcr_now(i-1,j,k) ) * eddlat

          ds_dzeta = ztrcr_now(i,j,k+1) - ztrcr_now(i,j,k)     ! ????????

          ! H2 at w-position:
          work_jk(j,k) = - Kvh_rho_at_w(i,j,k) * r_earth_recip               &
            &    * (                                         ds_dphi_at_w    &
            &      + (1.0_wp-beta_1(k)) * dzeta_dphi_at_w(i,j,k) * ds_dzeta )
        END DO

        DO k=2, ke
          DO j = jstart, jend

            ds_dphi_at_w = 0.25_wp * ( ( ztrcr_now(i,j+1,k-1) - ztrcr_now(i,j-1,k-1) )    &
              &                      + ( ztrcr_now(i,j+1,k  ) - ztrcr_now(i,j-1,k  ) ) ) * eddlat

            ds_dzeta = ztrcr_now(i,j,k) - ztrcr_now(i,j,k-1)

            ! H2 at w-position:
            work_jk(j,k) = - Kvh_rho_at_w(i,j,k) * r_earth_recip                &
              &    * ( ds_dphi_at_w                                             &
              &      + (1.0_wp-beta_1(k)) * dzeta_dphi_at_w(i,j,k) * ds_dzeta )
          END DO
        END DO

        !MB   DO k=1, ke-1
        DO k=2, ke-1
          DO j = jstart, jend

            ztrcr_tens(i,j,k)= ztrcr_tens(i,j,k) - rho_inv(i,j,k)                      &
              &  * r_earth_recip * dzeta_dphi(i,j,k) * ( work_jk(j,k+1) - work_jk(j,k) )

          END DO
        END DO

      END DO

    END IF

    ! --- ztrcr_tens: contributions by H3 ---

    ! the following tendency is already treated in 'implicit_vert_diffusion_uvwt_3D'

    !DO j = jstart, jend

    !  DO k=1, ke-1
    !    DO i = istart, iend

    !      ds_dzeta_at_w  = ztrcr_now(i,j,k) - ztrcr_now(i,j,k-1)

    !      ! H3 at w-position:
    !      work_ik(i,k) = - Kvh_rho_at_w(i,j,k) *               &
    !         &       (1.0_wp-beta_2(k)) * (-sqrtg_r_w(i,j,k)) * ds_dzeta_at_w
    !    END DO
    !  END DO

    !  k=ke
    !  DO i = istart, iend

    !    ds_dzeta_at_w  = ztrcr_now(i,j,k) - ztrcr_now(i,j,k-1)

    !    ! H3 at w-position:
    !    work_ik(i,k) = - Kvh_rho_at_w(i,j,k) *               &
    !       &       (1.0_wp-beta_2(k)) * (-sqrtg_r_w(i,j,k)) * ds_dzeta_at_w
    !  END DO

    !  DO k = 1, ke-1
    !    DO  i = istart, iend
    !
    !      ztrcr_tens(i,j,k)= ztrcr_tens(i,j,k) - rho_inv(i,j,k) *          &
    !         &   (-sqrtg_r_s(i,j,k)) * ( work_ik(i,k+1) - work_ik(i,k) )

    !    END DO
    !  END DO

    !END DO

  END DO

  IF ( nztrcr_3dturb > 0 ) THEN
    DEALLOCATE( Kvh_rho_at_u )
    DEALLOCATE( Kvh_rho_at_v )
  END IF

  ! === tketens =========================================

  IF ( lprog_tke ) THEN

    ALLOCATE( ztke    (1:ie, 1:je, 1:ke1) )
    ALLOCATE( ztketens(1:ie, 1:je, 1:ke1) )

    ztketens = 0.0_wp
    IF (itype_turb == 3) THEN !tke contains turb. velocity scale SQRT(2*TKE)
      ztke(:,:,:) = 0.5_wp * MAX( tke(:,:,:,ntke), ztkemin )**2
    ELSE
      ztke(:,:,:) = tke(:,:,:,nnow)
    END IF


    ALLOCATE( Kvh_rho_at_uw (ie, je, ke+1) )
    ALLOCATE( Kvh_rho_at_vw (ie, je, ke+1) )

    DO k=1, ke
      DO j = jstart, jend
        DO  i = istart, iend
          Kvh_rho_at_uw(i,j,k) = 0.5_wp * ( Kvh_rho_at_w(i,j,k) + Kvh_rho_at_w(i+1,j,k) ) 
          Kvh_rho_at_vw(i,j,k) = 0.5_wp * ( Kvh_rho_at_w(i,j,k) + Kvh_rho_at_w(i,j+1,k) ) 
        END DO
      END DO
    END DO


    ! --- tketens: contributions by H1 ---

    DO j = jstart, jend

      ! H1 at uw-position:
      DO k=1, ke
        DO i = istart, iend

          dtke_dlam_at_uw   = eddlon * ( ztke(i+1,j,k) - ztke(i,j,k) )

          work_ik(i,k) = - Kvh_rho_at_uw(i,j,k) * acrlat(j,1) * dtke_dlam_at_uw

        END DO
      END DO


      IF ( l3dturb_metr ) THEN

        k=1
        DO i = istart, iend

          dtke_dzeta_at_uw  = 0.5_wp *                         &
            &        ( ( ztke(i  ,j,k+1) - ztke(i  ,j,k) )     &
            &        + ( ztke(i+1,j,k+1) - ztke(i+1,j,k) ) )

          work_ik(i,k) = work_ik(i,k) - Kvh_rho_at_uw(i,j,k)     &
            &      * acrlat(j,1) * dzeta_dlam_at_uw(i,j,k) * dtke_dzeta_at_uw
        END DO

        DO k=2, ke-1
          DO i = istart, iend

            dtke_dzeta_at_uw  = 0.25_wp *                         &
              &        ( ( ztke(i  ,j,k+1) - ztke(i  ,j,k-1) )    &
              &        + ( ztke(i+1,j,k+1) - ztke(i+1,j,k-1) ) )

            work_ik(i,k) = work_ik(i,k) - Kvh_rho_at_uw(i,j,k)    &
              &     * acrlat(j,1) * dzeta_dlam_at_uw(i,j,k) * dtke_dzeta_at_uw
          END DO
        END DO

        k=ke
        DO i = istart, iend

          dtke_dzeta_at_uw = 0.5_wp *                          &
            &        ( ( ztke(i  ,j,k) - ztke(i  ,j,k-1) )     &
            &        + ( ztke(i+1,j,k) - ztke(i+1,j,k-1) ) )

          work_ik(i,k) = work_ik(i,k) - Kvh_rho_at_uw(i,j,k)    &
            &      * acrlat(j,1) * dzeta_dlam_at_uw(i,j,k) * dtke_dzeta_at_uw
        END DO

      END IF

      !MB    DO k = 2, ke
      DO k = 3, ke-1
        DO  i = istart, iend

          ztketens(i,j,k)= ztketens(i,j,k) - rho_inv_at_w(i,j,k)            &
            &   * acrlat(j,1) * ( work_ik(i,k) - work_ik(i-1,k) ) * eddlon
        END DO
      END DO

    END DO


    IF ( l3dturb_metr ) THEN

      DO j = jstart, jend

        k=1
        DO i = istart, iend

          dtke_dlam_at_s   = 0.25_wp * eddlon *                        &
            &               ( ( ztke(i+1,j,k  ) - ztke(i-1,j,k  ) )    &
            &               + ( ztke(i+1,j,k+1) - ztke(i-1,j,k+1) ) )

          dtke_dzeta_at_s  = ztke(i,j,k+1) - ztke(i,j,k)

          ! H1 at scalar position
          work_ik(i,k) = - Kvh_rho_at_s(i,j,k) *                &
            &       ( acrlat(j,1) * ( dtke_dlam_at_s            &
            &                       + (1.0_wp-beta_1(k)) * dzeta_dlam(i,j,k) * dtke_dzeta_at_s ) )
        END DO

        DO k=2, ke-1
          DO i = istart, iend

            dtke_dlam_at_s   = 0.25_wp * eddlon *                        &
              &               ( ( ztke(i+1,j,k  ) - ztke(i-1,j,k  ) )    &
              &               + ( ztke(i+1,j,k+1) - ztke(i-1,j,k+1) ) )

            dtke_dzeta_at_s  = ztke(i,j,k+1) - ztke(i,j,k)

            ! H1 at scalar position
            work_ik(i,k) = - Kvh_rho_at_s(i,j,k) *                &
              &       ( acrlat(j,1) * ( dtke_dlam_at_s            &
              &                       + (1.0_wp-beta_1(k)) * dzeta_dlam(i,j,k) * dtke_dzeta_at_s ) )
          END DO
        END DO

        k=ke
        DO i = istart, iend

          dtke_dlam_at_s   = 0.25_wp * eddlon *                        &
            &               ( ( ztke(i+1,j,k  ) - ztke(i-1,j,k  ) )    &
            &               + ( ztke(i+1,j,k+1) - ztke(i-1,j,k+1) ) )

          dtke_dzeta_at_s  = ztke(i,j,k+1) - ztke(i,j,k)

          ! tau13 at scalar position
          work_ik(i,k) = - Kvh_rho_at_s(i,j,k) *                &
            &       ( acrlat(j,1) * ( dtke_dlam_at_s            &
            &                       + (1.0_wp-beta_1(k)) * dzeta_dlam(i,j,k) * dtke_dzeta_at_s ) )
        END DO

        !MB   DO k = 2, ke
        DO k = 3, ke-1
          DO  i = istart, iend

            ztketens(i,j,k)= ztketens(i,j,k) - rho_inv_at_w(i,j,k)           &
              &   * acrlat(j,1) * dzeta_dlam_at_w(i,j,k) * ( work_ik(i,k) - work_ik(i,k-1) )

          END DO
        END DO

      END DO

    END IF

    ! --- tketens: contributions by H2 ---

    DO i = istart, iend

      ! H2 at vw-position:
      DO k=1, ke
        DO j = jstart, jend

          dtke_dphi_at_vw   = eddlat * ( ztke(i,j+1,k) - ztke(i,j,k) )

          work_jk(j,k) = - Kvh_rho_at_vw(i,j,k) * r_earth_recip * dtke_dphi_at_vw

        END DO
      END DO

      IF ( l3dturb_metr ) THEN

        k=1
        DO j = jstart, jend

          dtke_dzeta_at_vw  = 0.5_wp *                         &
            &        ( ( ztke(i,j  ,k+1) - ztke(i,j  ,k) )     &
            &        + ( ztke(i,j+1,k+1) - ztke(i,j+1,k) ) )

          work_jk(j,k) = work_jk(j,k) - Kvh_rho_at_vw(i,j,k) *      &
            &          r_earth_recip * dzeta_dphi_at_vw(i,j,k) * dtke_dzeta_at_vw
        END DO

        DO k=2, ke-1
          DO j = jstart, jend

            dtke_dzeta_at_vw  = 0.25_wp *                        &
              &        ( ( ztke(i,j  ,k+1) - ztke(i,j  ,k-1) )   &
              &        + ( ztke(i,j+1,k+1) - ztke(i,j+1,k-1) ) )

            work_jk(j,k) = work_jk(j,k) - Kvh_rho_at_vw(i,j,k)       &
              &     * r_earth_recip * dzeta_dphi_at_vw(i,j,k) * dtke_dzeta_at_vw
          END DO
        END DO

        k=ke
        DO j = jstart, jend

          dtke_dzeta_at_vw = 0.5_wp *                         &
            &        ( ( ztke(i,j  ,k) - ztke(i,j  ,k-1) )    &
            &        + ( ztke(i,j+1,k) - ztke(i,j+1,k-1) ) )

          work_jk(j,k) = work_jk(j,k) - Kvh_rho_at_vw(i,j,k)      &
            &     * r_earth_recip * dzeta_dphi_at_vw(i,j,k) * dtke_dzeta_at_vw
        END DO

      END IF

      !MB   DO k = 2, ke
      DO k = 3, ke-1
        DO j = jstart, jend

          ztketens(i,j,k)= ztketens(i,j,k) - rho_inv_at_w(i,j,k)                   &
            &    * r_earth_recip * ( work_jk(j,k) - work_jk(j-1,k) ) * eddlat
        END DO
      END DO

    END DO

    IF ( l3dturb_metr ) THEN

      DO i = istart, iend

        k=1
        DO j = jstart, jend

          dtke_dphi_at_s   = 0.25_wp * eddlat *                       &
            &               ( ( ztke(i,j+1,k  ) - ztke(i,j-1,k  ) )   &
            &               + ( ztke(i,j+1,k+1) - ztke(i,j-1,k+1) ) )

          dtke_dzeta_at_s  = ztke(i,j,k+1) - ztke(i,j,k)

          ! H2 at scalar position
          work_jk(j,k) = - Kvh_rho_at_s(i,j,k) *                &
            &       ( r_earth_recip * ( dtke_dphi_at_s          &
            &                         + (1.0_wp-beta_1(k)) * dzeta_dphi(i,j,k) * dtke_dzeta_at_s ) )
        END DO

        DO k=2, ke-1
          DO j = jstart, jend

            dtke_dphi_at_s   = 0.25_wp * eddlat *                       &
              &               ( ( ztke(i,j+1,k  ) - ztke(i,j-1,k  ) )   &
              &               + ( ztke(i,j+1,k+1) - ztke(i,j-1,k+1) ) )

            dtke_dzeta_at_s  = ztke(i,j,k+1) - ztke(i,j,k)

            ! H2 at scalar position
            work_jk(j,k) = - Kvh_rho_at_s(i,j,k) *                  &
              &       ( r_earth_recip * ( dtke_dphi_at_s            &
              &                         + (1.0_wp-beta_1(k)) * dzeta_dphi(i,j,k) * dtke_dzeta_at_s ) )
          END DO
        END DO

        k=ke
        DO j = jstart, jend

          dtke_dphi_at_s   = 0.25_wp * eddlat *                        &
            &               ( ( ztke(i,j+1,k  ) - ztke(i,j-1,k  ) )    &
            &               + ( ztke(i,j+1,k+1) - ztke(i,j-1,k+1) ) )

          dtke_dzeta_at_s  = ztke(i,j,k+1) - ztke(i,j,k)

          ! H2 at scalar position
          work_jk(j,k) = - Kvh_rho_at_s(i,j,k) *                &
            &       ( r_earth_recip * ( dtke_dphi_at_s          &
            &                         + dzeta_dphi(i,j,k) * dtke_dzeta_at_s ) )
        END DO

        !MB  DO k = 2, ke
        DO k = 3, ke-1
          DO  j = jstart, jend

            ztketens(i,j,k)= ztketens(i,j,k) - rho_inv_at_w(i,j,k)           &
              &   * r_earth_recip * dzeta_dphi_at_w(i,j,k) * ( work_jk(j,k) - work_jk(j,k-1) )

          END DO
        END DO

      END DO

    END IF


    ! --- tketens: contributions by H3 ---

    ! the following tendency is already treated in 'implicit_vert_diffusion_uvwt_3D'

    !DO j = jstart, jend

    !  DO k=1, ke-1
    !    DO i = istart, iend

    !      dtke_dzeta_at_s  = ztke(i,j,k+1) - ztke(i,j,k)

    !      ! H3 at 2-position:
    !      work_ik(i,k) = - Kvh_rho_at_s(i,j,k) *               &
    !         &       (1.0_wp-beta_2(k)) * (-sqrtg_r_s(i,j,k)) * dtke_dzeta_at_s
    !    END DO
    !  END DO

    !  k=ke
    !  DO i = istart, iend

    !    dtke_dzeta_at_s  = ztke(i,j,k) - ztke(i,j,k-1)

    !    ! H3 at 2-position:
    !    work_ik(i,k) = - Kvh_rho_at_s(i,j,k) *               &
    !       &       (1.0_wp-beta_2(k)) * (-sqrtg_r_s(i,j,k)) * dtke_dzeta_at_s
    !  END DO

    !  DO k = 2, ke
    !    DO  i = istart, iend

    !      ztketens(i,j,k)= ztketens(i,j,k) - rho_inv_at_w(i,j,k) *          &
    !         &   (-sqrtg_r_w(i,j,k)) * ( work_ik(i,k) - work_ik(i,k-1) )

    !    END DO
    !  END DO

    !END DO


    IF (itype_turb == 3) THEN
      DO k=2, ke
        DO j=jstart, jend
          DO i=istart, iend
            ! Here, tke is sqrt(2*tke) and ztketens the "real" TKE tendency.
            ! The sqrt(2*tke) tendency is obtained by dividing
            ! the "real" TKE tendency by SQRT(2*tke):
            ztketens(i,j,k) = ztketens(i,j,k) / tke(i,j,k,ntke)
            ! limit the tendency so that SQRT(2*tke) stays >= 0.0:
            ztketens(i,j,k) = MAX( -tke(i,j,k,ntke) / dt, ztketens(i,j,k) )
            tketens(i,j,k)  = tketens(i,j,k) + ztketens(i,j,k)
          END DO
        END DO
      END DO
    ELSE
      DO k=2, ke
        DO j=jstart, jend
          DO i=istart, iend
            ! Here, tke is the "real" tke and ztketens the "real" TKE tendency.
            ! limit the tendency so that tke stays >= 0.0:
            ztketens(i,j,k) = MAX( -tke(i,j,k,nnow) / dt, ztketens(i,j,k))
            tketens(i,j,k)  = tketens(i,j,k) + ztketens(i,j,k)
          END DO
        END DO
      END DO
    END IF

    DEALLOCATE( ztke )
    DEALLOCATE( ztketens )

    DEALLOCATE( Kvh_rho_at_uw )
    DEALLOCATE( Kvh_rho_at_vw )

  END IF   ! IF ( lprog_tke ) 


  DEALLOCATE( rho_inv_at_w )

  DEALLOCATE( rho_inv )

  DEALLOCATE( work_ik )
  DEALLOCATE( work_jk )

END SUBROUTINE explicit_horizontal_diffusion


!==============================================================================
!==============================================================================

SUBROUTINE implicit_vert_diffusion_uvwt

  !----------------------------------------------------------------------------
  !
  ! Description:
  !   The vertical diffusion as a slow tendency 
  !   is computed here for all variables the vertical diffusion acts on. 
  !
  ! Method:
  !   The vertical diffusion is solved by a vertically implicit 
  !   scheme (modified Crank-Nicolson).
  !
  !----------------------------------------------------------------------------

  ! Declarations:

  ! Local scalars:
  ! -------------
  INTEGER (KIND=iintegers) ::  &
    i,  j,  k              !  Loop indices in lon., lat. and vert. direction

  INTEGER (KIND=iintegers) :: &
    km1, kp1

  INTEGER (KIND=iintegers) :: izerror
  CHARACTER (LEN=255)      :: yzerrmsg
  CHARACTER (LEN=50)       :: yzroutine

  REAL    (KIND=wp)        ::  &
    zgat, zgct,          & !
    zag, zas,            & !
    zcg, zcs,            & !
    zbg, zdg,            & !
    znew, zz,            & !
    zzdtr, ztmcmq,       & !
    zlhvocp                !

  ! Local (automatic) arrays:
  ! ------------------------
  REAL    (KIND=wp)        ::  &
    zkm     (ie,je   ),    & !
    zgatz   (ie,je,ke),    & ! 
    zgctz   (ie,je,ke),    & ! 
    zc      (ie,je,ke),    & ! Upper main diagonal of ...
    zd1     (ie,je,ke),    & ! Right hand side of ...
    ze      (ie,je,ke)       ! Soluton vector of ...

  ! Tracer pointers:
  !-----------------
  REAL (KIND=wp),     POINTER :: &
    qc  (:,:,:)  => NULL()             ! QC at tlev=nnow

  ! End of header
  !============================================================================

  !----------------------------------------------------------------------------
  ! Begin Subroutine implicit_vert_diffusion_uvwt
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Section 1: Some preparations
  !----------------------------------------------------------------------------

  izerror   = 0_iintegers
  yzerrmsg  = ''
  yzroutine = 'implicit_vert_diffusion_uvwt'

  IF ( izdebug >= 10 ) THEN
    WRITE(*,*) "subr. [", TRIM(yzroutine), "] ..."
  END IF 

  ! Setting of reciprocal time step
  zzdtr = 1.0_wp / dt

  !$ser savepoint VerticalDiffusionUnittest.DoUVWT-in LargeTimeStep=ntstep
  !$ser data u_nnow=u(:,:,:,nnow) u_tens=utens &
  !$ser&     v_nnow=v(:,:,:,nnow) v_tens=vtens &
  !$ser&     w_nnow=w(:,:,:,nnow) w_tens=wtens &
  !$ser&     t_nnow=t(:,:,:,nnow) t_tens=ttens &
  !$ser&     pp_nnow=pp(:,:,:,nnow) rho=rho    &
  !$ser&     p_s_nnew=ps(:,:,nnew) p_s_nnow=ps(:,:,nnow) &
  !$ser&     t_g_nnew=t_g(:,:,nnew) t_g_nnow=t_g(:,:,nnow) &
  !$ser&     tkvm=tkvm tkvh=tkvh tch=tch tcm=tcm
  !$ser tracer QV#surf@nnow

  IF ( lvertdiff ) THEN
    
    !--------------------------------------------------------------------------
    ! Section 2: Temperature 
    !--------------------------------------------------------------------------

    IF ( lmoist_turb ) THEN
      
      ! Retrieve the required microphysics tracers
      CALL trcr_get(izerror, idt_qc, ptr_tlev = nnow, ptr = qc)
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF

      zlhvocp = lh_v / cp_d
      
      ! Top layer
      DO j = jstart, jend
        DO i = istart, iend
          zgct       = - zkh(i,j,2)*zsqrtgrho_r_s(i,j,1)
          zcg        = zgct*a1t(2)*zpianf(i,j,2)/zpia(i,j,2)
          zbg        = zzdtr - zgct*a1t(2)*zpianf(i,j,2)/zpia(i,j,1)
          zdg        = zzdtr * t(i,j,1,nnow) + ttens(i,j,1)
          zdg        = zdg - zgct*zpianf(i,j,2) * (                          &
                          a2t(2)*  ( ztheta(i,j,2) - ztheta(i,j,1) )         &
                       - zlhvocp * ( qc(i,j,2)/zpia(i,j,2)                   &
                                   - qc(i,j,1)/zpia(i,j,1) ) )
          zd1(i,j,1) = zdg/zbg
          zc(i,j,1)  = zcg/zbg
        ENDDO
      ENDDO

      ! The layers from k=2 to k=ke-1
      DO k = 2, ke-1
        DO j = jstart, jend
          DO i = istart, iend
            zgat       = - zkh(i,j,k  )*zsqrtgrho_r_s(i,j,k)
            zgct       = - zkh(i,j,k+1)*zsqrtgrho_r_s(i,j,k)
            zag        = zgat*a1t(k  )*zpianf(i,j,k  )/zpia(i,j,k-1)
            zcg        = zgct*a1t(k+1)*zpianf(i,j,k+1)/zpia(i,j,k+1)
            zbg        = zzdtr - zgat*a1t(k  )*zpianf(i,j,k  )/zpia(i,j,k)  &
                               - zgct*a1t(k+1)*zpianf(i,j,k+1)/zpia(i,j,k)
            zdg        = zzdtr * t(i,j,k,nnow) + ttens(i,j,k)
            zdg        = zdg - zgat * zpianf(i,j,k) *  (                     &
                   a2t(k)*  ( ztheta(i,j,k-1) - ztheta(i,j,k) )              &
                - zlhvocp * ( qc(i,j,k-1)/zpia(i,j,k-1)                      &
                            - qc(i,j,k  )/zpia(i,j,k  ) ) )
            zdg        = zdg - zgct * zpianf(i,j,k+1)* (                     &
                   a2t(k+1)*( ztheta(i,j,k+1) - ztheta(i,j,k) )              &
                - zlhvocp * ( qc(i,j,k+1)/zpia(i,j,k+1)                      &
                            - qc(i,j,k  )/zpia(i,j,k  ) ) )
            zz         = 1.0_wp/( zbg - zag*zc(i,j,k-1) )
            zc (i,j,k) = zcg*zz
            zd1(i,j,k) = ( zdg -zag*zd1(i,j,k-1) )*zz
          ENDDO
        ENDDO
      ENDDO

      ! The bottom layer
      DO j = jstart, jend
        DO i = istart, iend
          zgat       = - zkh (i,j,ke)*zsqrtgrho_r_s(i,j,ke)
          zgct       = - ztch(i,j   )*zsqrtgrho_r_s(i,j,ke)
          zag        = zgat*a1t(ke)*zpianf(i,j,ke)/zpia(i,j,ke-1)
          zbg        = zzdtr - zgat*a1t(ke )*zpianf(i,j,ke )/zpia(i,j,ke)   &
                             - zgct*za1t_surf*zpianf(i,j,ke1)/zpia(i,j,ke)
          zdg        = zzdtr * t(i,j,ke,nnow) + ttens(i,j,ke)
          zdg        = zdg - zgat*zpianf(i,j,ke ) * (                        &
                 a2t(ke)* ( ztheta(i,j,ke-1) - ztheta(i,j,ke) )              &
              - zlhvocp * ( qc(i,j,ke-1)/zpia(i,j,ke-1)                      &
                          - qc(i,j,ke  )/zpia(i,j,ke  ) ) )
          zdg        = zdg + zgct*zpianf(i,j,ke1)* (                         &
                za2t_surf * ztheta(i,j,ke)                                   &
              - zlhvocp *   qc(i,j,ke)/zpia(i,j,ke  ) )                      &
                           - zgct*t_g(i,j,nnow)
          znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
          ttens(i,j,ke) = ( znew - t(i,j,ke,nnow) ) * zzdtr
          ze   (i,j,ke) = znew
        ENDDO
      ENDDO

    ELSE

      ! Top layer
      DO j = jstart, jend
        DO i = istart, iend
          zgct       = - zkh(i,j,2)*zsqrtgrho_r_s(i,j,1)
          zcg        = zgct*a1t(2)*zpianf(i,j,2)/zpia(i,j,2)
          zbg        = zzdtr - zgct*a1t(2)*zpianf(i,j,2)/zpia(i,j,1)
          zdg        = zzdtr * t(i,j,1,nnow) + ttens(i,j,1)            &
                      - zgct*a2t(2)*zpianf(i,j,2) *                    &
                       ( ztheta(i,j,2) - ztheta(i,j,1) )
          zd1(i,j,1) = zdg/zbg
          zc(i,j,1)  = zcg/zbg
        ENDDO
      ENDDO

      ! The layers from k=2 to k=ke-1
      DO k = 2, ke-1
        DO j = jstart, jend
          DO i = istart, iend
            zgat       = - zkh(i,j,k  )*zsqrtgrho_r_s(i,j,k)
            zgct       = - zkh(i,j,k+1)*zsqrtgrho_r_s(i,j,k)
            zag        = zgat*a1t(k  )*zpianf(i,j,k  )/zpia(i,j,k-1)
            zcg        = zgct*a1t(k+1)*zpianf(i,j,k+1)/zpia(i,j,k+1)
            zbg        = zzdtr - zgat *a1t(k  )*zpianf(i,j,k  )/zpia(i,j,k)  &
                       - zgct*a1t(k+1)*zpianf(i,j,k+1)/zpia(i,j,k)
            zdg        = zzdtr * t(i,j,k,nnow) + ttens(i,j,k)                &
              - zgat* a2t(k)*zpianf(i,j,k) *                                 &
                ( ztheta(i,j,k-1) - ztheta(i,j,k) )                          &
              - zgct* a2t(k+1)*zpianf(i,j,k+1)*                              &
                ( ztheta(i,j,k+1) - ztheta(i,j,k) )
            zz         = 1.0_wp/( zbg - zag*zc(i,j,k-1) )
            zc (i,j,k) = zcg*zz
            zd1(i,j,k) = ( zdg -zag*zd1(i,j,k-1) )*zz
          ENDDO
        ENDDO
      ENDDO

      ! The bottom layer
      DO j = jstart, jend
        DO i = istart, iend
          zgat       = - zkh (i,j,ke)*zsqrtgrho_r_s(i,j,ke)
          zgct       = - ztch(i,j   )*zsqrtgrho_r_s(i,j,ke)
          zag        = zgat*a1t(ke)*zpianf(i,j,ke)/zpia(i,j,ke-1)
          zbg        = zzdtr - zgat* a1t(ke )*zpianf(i,j,ke )/zpia(i,j,ke)   &
            - zgct*za1t_surf*zpianf(i,j,ke1)/zpia(i,j,ke)
          zdg        = zzdtr * t(i,j,ke,nnow) + ttens(i,j,ke)                &
            - zgat* a2t(ke)*zpianf(i,j,ke ) *                                &
              ( ztheta(i,j,ke-1) - ztheta(i,j,ke) )                          &
            + zgct*za2t_surf *zpianf(i,j,ke1)*                               &
              ztheta(i,j,ke) - zgct*t_g(i,j,nnow)
          znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
          ttens(i,j,ke) = ( znew - t(i,j,ke,nnow) ) * zzdtr
          ze   (i,j,ke) = znew
        ENDDO
      ENDDO

    END IF
  
    ! Backsubstitution and storage of the complete slow tendencies

    DO k = ke-1, 1, -1
      DO j = jstart, jend
        DO i = istart, iend
          ze   (i,j,k) = zd1(i,j,k) - zc(i,j,k)*ze(i,j,k+1)
          ttens(i,j,k) = ( ze(i,j,k) - t(i,j,k,nnow) ) * zzdtr
        ENDDO
      ENDDO
    ENDDO

    ! Calculation of the sensible heat flux at the surface
    ! This flux is integrated in time

    DO j = jstart, jend
      DO  i  = istart , iend
        shfl_s (i,j) = - ztch(i,j)*cp_d*                  &
          ( za2t_surf*(t_g(i,j,nnow) - zpianf(i,j,ke1)*   &
                       t  (i,j,ke,nnow)/zpia(i,j,ke) ) +  &
            za1t_surf*(t_g(i,j,nnew) - zpianf(i,j,ke1)*   &
                       ze(i,j,ke)/zpia(i,j,ke) ) )
      ENDDO
    ENDDO

    !--------------------------------------------------------------------------
    ! Section : Horizontal wind velocity u 
    !--------------------------------------------------------------------------

    ! Top layer:   k = 1
    DO j = jstartu, jendu
      DO i = istartu, iendu
        zkm(i,j)   = 0.5_wp*( ztmkvm(i,j,2) + ztmkvm(i+1,j,2) )
        zgct       = - zkm(i,j)*zsqrtgrho_r_u(i,j,1)
        zcg        = zgct*a1t(2)
        zbg        = zzdtr - zcg
        zdg        = zzdtr * u(i,j,1,nnow) + utens(i,j,1)    &
          - a2t(2)*zgct * ( u(i,j,2,nnow) - u(i,j,1,nnow) )
        zd1(i,j,1) = zdg/zbg
        zc (i,j,1) = zcg/zbg
      ENDDO
    ENDDO

    ! The layers from k=2 to k=ke-1
    DO k = 2, ke-1
      DO j = jstartu, jendu
        DO i = istartu, iendu
          zgat       = - zkm(i,j)*zsqrtgrho_r_u(i,j,k)
          zkm(i,j)   = 0.5_wp*( ztmkvm(i,j,k+1) + ztmkvm(i+1,j,k+1) )
          zgct       = - zkm(i,j)*zsqrtgrho_r_u(i,j,k)
          zag        = zgat*a1t(k)
          zas        = zgat*a2t(k)
          zcg        = zgct*a1t(k+1)
          zcs        = zgct*a2t(k+1)
          zbg        = zzdtr - zag - zcg
          zdg        = zzdtr * u(i,j,k,nnow) + utens(i,j,k)  &
            - zas * ( u(i,j,k-1,nnow) - u(i,j,k,nnow) )      &
            - zcs * ( u(i,j,k+1,nnow) - u(i,j,k,nnow) )
          zz         = 1.0_wp/( zbg - zag*zc(i,j,k-1) )
          zc (i,j,k)  = zcg*zz
          zd1(i,j,k)  = ( zdg -zag*zd1(i,j,k-1) )*zz
        ENDDO
      ENDDO
    ENDDO

    ! The bottom layer:  k = ke
    ! Including the calculation of the u-momentum flux at the surface
    DO j = jstartu, jendu
      DO i = istartu, iendu
        zgat       = - zkm(i,j)*zsqrtgrho_r_u(i,j,ke)
        ztmcmq     = 0.5_wp*( ztcm(i,j) + ztcm(i+1,j) )
        zgct       = - ztmcmq*zsqrtgrho_r_u(i,j,ke)
        zag        = zgat*a1t(ke)
        zas        = zgat*a2t(ke)
        zcg        = a1t(ke1)*zgct
        zcs        = a2t(ke1)*zgct
        zbg        = zzdtr - zag - zcg
        zdg        = zzdtr * u(i,j,ke,nnow) + utens(i,j,ke)  &
          - zas * ( u(i,j,ke-1,nnow) - u(i,j,ke,nnow) )      &
          + zcs * u(i,j,ke,nnow)
        znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
        utens(i,j,ke) = ( znew - u(i,j,ke,nnow) ) * zzdtr
        umfl_s(i,j  ) = ztmcmq*( a2t(ke1)*u(i,j,ke,nnow) + a1t(ke1)*znew )
        ze   (i,j,ke) = znew
      ENDDO
    ENDDO

    ! Backsubstitution and storage of the complete slow tendencies

    DO k = ke-1, 1, -1
      DO j = jstartu, jendu
        DO i = istartu, iendu
          ze   (i,j,k) = zd1(i,j,k) - zc(i,j,k)*ze(i,j,k+1)
          utens(i,j,k) = ( ze(i,j,k) - u(i,j,k,nnow) ) * zzdtr
        ENDDO
      ENDDO
    ENDDO

    !--------------------------------------------------------------------------
    ! Section : Horizontal wind velocity v
    !--------------------------------------------------------------------------

    ! Top layer  k=1
    DO j = jstartv, jendv
      DO i = istartv, iendv
        zkm(i,j)   = 0.5_wp*( ztmkvm(i,j,2) + ztmkvm(i,j+1,2) )
        zgct       = - zkm(i,j)*zsqrtgrho_r_v(i,j,1)
        zcg        = zgct*a1t(2)
        zbg        = zzdtr - zcg
        zdg        = zzdtr * v(i,j,1,nnow) + vtens(i,j,1)     &
          - a2t(2)*zgct * ( v(i,j,2,nnow) - v(i,j,1,nnow) )
        zd1(i,j,1) = zdg/zbg
        zc (i,j,1) = zcg/zbg
      ENDDO
    ENDDO

    ! The layers from k=2 to k=ke-1
    DO  k = 2, ke-1
      DO j = jstartv, jendv
        DO i = istartv, iendv
          zgat       = - zkm(i,j)*zsqrtgrho_r_v(i,j,k)
          zkm(i,j)   = 0.5_wp*( ztmkvm(i,j,k+1) + ztmkvm(i,j+1,k+1) )
          zgct       = - zkm(i,j)*zsqrtgrho_r_v(i,j,k)
          zag        = zgat*a1t(k)
          zas        = zgat*a2t(k)
          zcg        = zgct*a1t(k+1)
          zcs        = zgct*a2t(k+1)
          zbg        = zzdtr - zag - zcg
          zdg        = zzdtr * v(i,j,k,nnow) + vtens(i,j,k)  &
            - zas * ( v(i,j,k-1,nnow) - v(i,j,k,nnow) )      &
            - zcs * ( v(i,j,k+1,nnow) - v(i,j,k,nnow) )
          zz         = 1.0_wp/( zbg - zag*zc(i,j,k-1) )
          zc (i,j,k) = zcg*zz
          zd1(i,j,k) = ( zdg -zag*zd1(i,j,k-1) )*zz
        ENDDO
      ENDDO
    ENDDO

    ! The bottom layer k=ke
    ! Including the calculation of the v-momentum flux at the surface
    DO j = jstartv, jendv
      DO i = istartv, iendv
        zgat       = - zkm(i,j)*zsqrtgrho_r_v(i,j,ke)
        ztmcmq     = 0.5_wp*( ztcm(i,j) + ztcm(i,j+1) )
        zgct       = - ztmcmq*zsqrtgrho_r_v(i,j,ke)
        zag        = zgat*a1t(ke)
        zas        = zgat*a2t(ke)
        zcg        = a1t(ke1)*zgct
        zcs        = a2t(ke1)*zgct
        zbg        = zzdtr - zag - zcg
        zdg        = zzdtr * v(i,j,ke,nnow) + vtens(i,j,ke)  &
          - zas * ( v(i,j,ke-1,nnow) - v(i,j,ke,nnow) )      &
          + zcs * v(i,j,ke,nnow)
        znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
        vtens(i,j,ke) = ( znew - v(i,j,ke,nnow) ) * zzdtr
        vmfl_s(i,j  ) = ztmcmq*( a2t(ke1)*v(i,j,ke,nnow) + a1t(ke1)*znew )
        ze   (i,j,ke) = znew
      ENDDO
    ENDDO

    ! Backsubstitution and storage of the complete slow tendencies

    DO k = ke-1, 1, -1
      DO j = jstartv, jendv
        DO i = istartv, iendv
          ze   (i,j,k) = zd1(i,j,k) - zc(i,j,k)*ze(i,j,k+1)
          vtens(i,j,k) = ( ze(i,j,k) - v(i,j,k,nnow) ) * zzdtr
        ENDDO
      ENDDO
    ENDDO

  END IF
  
  !----------------------------------------------------------------------------
  ! Section : Vertical wind velocity w
  !----------------------------------------------------------------------------

  IF ( lvertdiff_w ) THEN
    
    DO j = jstart, jend

      DO  k = 2, ke
        km1 = MAX( 2, k-1 )
        kp1 = MIN( ke, k+1 )
      
        ! Precalculate some help variables to avoid divisions in subsequent code
        DO i = istart, iend
          zgatz(i,j,k) = - zsqrtgrho_r_w(i,j,k)   &
                         * ( ztmkvw(i,j,k)+ztmkvw(i,j,km1) )
          zgctz(i,j,k) = - zsqrtgrho_r_w(i,j,k)   &
                         * ( ztmkvw(i,j,kp1)+ztmkvw(i,j,k) )
        ENDDO
      ENDDO
      
      ! Top layer       
      DO i = istart, iend
        zag        = 0.5_wp*(a1t(2)+a1t(1)) * zgatz(i,j,2)
        zas        = 0.5_wp*(a2t(2)+a2t(1)) * zgatz(i,j,2)
        zcg        = 0.5_wp*(a1t(3)+a1t(2)) * zgctz(i,j,2)
        zcs        = 0.5_wp*(a2t(3)+a2t(2)) * zgctz(i,j,2)
        zbg        = zzdtr - zag - zcg
        zdg        = zzdtr * w(i,j,2,nnow) + wtens(i,j,2)                    &
                     + zas * w(i,j,2,nnow)                                   &
                     - zcs * ( w(i,j,3,nnow) - w(i,j,2,nnow) )
        zd1(i,j,2)  = zdg/zbg
        zc (i,j,2)  = zcg/zbg
      ENDDO

      ! The layers from k=3 to k=ke-1
      DO k = 3, ke-1
        DO i = istart, iend
          zag        = 0.5_wp*(a1t(k)+a1t(k-1)) * zgatz(i,j,k)
          zas        = 0.5_wp*(a2t(k)+a2t(k-1)) * zgatz(i,j,k)
          zcg        = 0.5_wp*(a1t(k+1)+a1t(k)) * zgctz(i,j,k)
          zcs        = 0.5_wp*(a2t(k+1)+a2t(k)) * zgctz(i,j,k)
          zbg        = zzdtr - zag - zcg
          zdg        = zzdtr * w(i,j,k,nnow) + wtens(i,j,k)                  &
                     - zas * ( w(i,j,k-1,nnow) - w(i,j,k,nnow) )             &
                     - zcs * ( w(i,j,k+1,nnow) - w(i,j,k,nnow) )
          zz         = 1.0_wp/( zbg - zag*zc(i,j,k-1) )
          zc (i,j,k)  = zcg*zz
          zd1(i,j,k)  = ( zdg - zag*zd1(i,j,k-1) )*zz
        ENDDO
      ENDDO

      ! The bottom layer
      DO i = istart, iend
        zag        = 0.5_wp*(a1t(ke)+a1t(ke-1)) * zgatz(i,j,ke)
        zas        = 0.5_wp*(a2t(ke)+a2t(ke-1)) * zgatz(i,j,ke)
        zcg        = 0.5_wp*(a1t(ke1)+a1t(ke))  * zgctz(i,j,ke)
        zcs        = 0.5_wp*(a2t(ke1)+a2t(ke))  * zgctz(i,j,ke)
        zbg        = zzdtr - zag - zcg
        zdg        = zzdtr * w(i,j,ke,nnow) + wtens(i,j,ke)                 &
                     - zas * ( w(i,j,ke-1,nnow) - w(i,j,ke,nnow) )          &
                     - zcs * ( w(i,j,ke1,nnow) - w(i,j,ke,nnow))            &
                     - zcg * w(i,j,ke1,nnow) 
        znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
        wtens(i,j,ke) = ( znew - w(i,j,ke,nnow) ) * zzdtr
        ze   (i,j,ke) = znew
      ENDDO

      ! Backsubstitution and storage of the complete slow tendencies

      DO k = ke-1, 2, -1
        DO i = istart, iend
          ze    (i,j,k) = zd1(i,j,k) - zc(i,j,k)*ze(i,j,k+1)
          wtens (i,j,k) = ( ze(i,j,k) - w(i,j,k,nnow) ) * zzdtr
        ENDDO
      ENDDO
    ENDDO

  END IF

!$ser savepoint VerticalDiffusionUnittest.DoUVWT-out LargeTimeStep=ntstep
!$ser data u_tens=utens v_tens=vtens w_tens=wtens t_tens=ttens &
!$ser&     shfl_s=shfl_s umfl_s=umfl_s vmfl_s=vmfl_s

  !----------------------------------------------------------------------------
  ! End of subroutine implicit_vert_diffusion_uvwt
  !----------------------------------------------------------------------------

END SUBROUTINE implicit_vert_diffusion_uvwt


!==============================================================================
!==============================================================================


SUBROUTINE implicit_vert_diffusion_uvwt_3D

  !----------------------------------------------------------------------------
  !
  ! Description:
  !   The vertical diffusion as a slow tendency 
  !   is computed here for all variables the vertical diffusion acts on. 
  !
  ! output:
  !   utens, vtens, wtens, ttens (and possible tracer tendencies)
  !
  ! input:
  !   utens, vtens, wtens, ttens (and possible tracer tendencies)
  !   Kvm, Kvh, 
  !
  ! Method:
  !   The vertical diffusion is solved by a vertically implicit 
  !   scheme (modified Crank-Nicolson).
  !
  !----------------------------------------------------------------------------

  ! Declarations:

  ! Local scalars:
  ! -------------
  INTEGER (KIND=iintegers) ::  &
    i,  j,  k,      &    !  Loop indices in lon., lat. and vert. direction
    k_s2w_up, k_s2w_low

  INTEGER (KIND=iintegers) :: &
    km1, kp1

  INTEGER (KIND=iintegers) :: izerror
  CHARACTER (LEN=255)      :: yzerrmsg
  CHARACTER (LEN=50)       :: yzroutine

  REAL    (KIND=wp)        ::  &
    zgat, zgct,          & !
    zag, zas,            & !
    zcg, zcs,            & !
    zbg, zdg,            & !
    z2ag, z2cg,          & !
    znew, zz,            & !
    zzdtr, ztmcmq,       & !
    zlhvocp                !

  ! Local (automatic) arrays:
  ! ------------------------
  REAL    (KIND=wp)        ::  &
    zkm     (ie,je   ),    & !
    zgatz   (ie,je,ke),    & ! 
    zgctz   (ie,je,ke),    & ! 
    zc      (ie,je,ke),    & ! Upper main diagonal of ...
    zd1     (ie,je,ke),    & ! Right hand side of ...
    ze      (ie,je,ke)       ! Soluton vector of ...

  ! Tracer pointers:
  !-----------------
  REAL (KIND=wp),     POINTER :: &
    qc  (:,:,:)  => NULL()             ! QC at tlev=nnow

  REAL (KIND=wp)        ::  &
    rho_inv, rho_inv_at_u,    rho_inv_at_v,    rho_inv_at_w

  REAL (KIND=wp)        ::  &
    r_earth_recip,    &
    z_pi

  REAL (KIND=wp)        ::  &
    dzeta_dlam_at_v, dzeta_dphi_at_u

  REAL (KIND=wp), ALLOCATABLE, DIMENSION(:,:,:) ::    &
    dzeta_dlam_at_vw,  &
    dzeta_dphi_at_uw


  ! End of header
  !============================================================================

  !----------------------------------------------------------------------------
  ! Begin Subroutine implicit_vert_diffusion_uvwt_3D
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Section 1: Some preparations
  !----------------------------------------------------------------------------

  izerror   = 0_iintegers
  yzerrmsg  = ''
  yzroutine = 'implicit_vert_diffusion_uvwt_3D'

  IF ( izdebug >= 10 ) THEN
    WRITE(*,*) "subr. [", TRIM(yzroutine), "] ..."
  END IF

  ! Setting of reciprocal time step
  zzdtr = 1.0_wp / dt

  r_earth_recip = 1.0_wp / r_earth


  ALLOCATE( dzeta_dlam_at_vw (ie,je,ke+1) )
  ALLOCATE( dzeta_dphi_at_uw (ie,je,ke+1) )

  DO k = 1, ke+1
    DO j = jstartv, jendv
      DO i = istartv, iendv
        dzeta_dlam_at_vw(i,j,k) = 0.5_wp * ( dzeta_dlam_at_w(i,j,k) + dzeta_dlam_at_w(i,j+1,k) )
      ENDDO
    ENDDO
  ENDDO

  DO k = 1, ke+1
    DO j = jstartu, jendu
      DO i = istartu, iendu
        dzeta_dphi_at_uw(i,j,k) = 0.5_wp * ( dzeta_dphi_at_w(i,j,k) + dzeta_dphi_at_w(i+1,j,k) )
      ENDDO
    ENDDO
  ENDDO



  IF ( lvertdiff ) THEN

    !--------------------------------------------------------------------------
    ! Section 2: Temperature 
    !--------------------------------------------------------------------------

    IF ( lmoist_turb ) THEN

      ! Retrieve the required microphysics tracers
      CALL trcr_get(izerror, idt_qc, ptr_tlev = nnow, ptr = qc)
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF

      zlhvocp = lh_v / cp_d

      ! Top layer
      DO j = jstart, jend
        DO i = istart, iend
          zgct       = - zkh(i,j,2)*zsqrtgrho_r_s(i,j,1)
          zcg        = zgct*a1t(2)*zpianf(i,j,2)/zpia(i,j,2)
          zbg        = zzdtr - zgct*a1t(2)*zpianf(i,j,2)/zpia(i,j,1)
          zdg        = zzdtr * t(i,j,1,nnow) + ttens(i,j,1)
          zdg        = zdg - zgct*zpianf(i,j,2) * (                          &
            &             a2t(2)*  ( ztheta(i,j,2) - ztheta(i,j,1) )         &
            &          - zlhvocp * ( qc(i,j,2)/zpia(i,j,2)                   &
            &                      - qc(i,j,1)/zpia(i,j,1) ) )
          zd1(i,j,1) = zdg/zbg
          zc(i,j,1)  = zcg/zbg
        ENDDO
      ENDDO

      ! The layers from k=2 to k=ke-1
      DO k = 2, ke-1
        DO j = jstart, jend
          DO i = istart, iend

            rho_inv = 1.0_wp / rho(i,j,k)

            ! the explicit counterparts of z2ag, z2cg are in ttens (via 'explicit_horizontal_diffusion'):
            z2ag = rho_inv * Kvh_rho_pi_at_w(i,j,k)   * beta_1(k) *                                   &
              &     ( acrlat(j,1)   * dzeta_dlam(i,j,k) * acrlat(j,1)   * dzeta_dlam_at_w(i,j,k)      &
              &     + r_earth_recip * dzeta_dphi(i,j,k) * r_earth_recip * dzeta_dphi_at_w(i,j,k) )

            z2cg = rho_inv * Kvh_rho_pi_at_w(i,j,k+1) * beta_1(k+1) *                                 &
              &     ( acrlat(j,1)   * dzeta_dlam(i,j,k) * acrlat(j,1)   * dzeta_dlam_at_w(i,j,k+1)    &
              &     + r_earth_recip * dzeta_dphi(i,j,k) * r_earth_recip * dzeta_dphi_at_w(i,j,k+1) )

            zgat       = - zkh(i,j,k  ) * zsqrtgrho_r_s(i,j,k)
            zgct       = - zkh(i,j,k+1) * zsqrtgrho_r_s(i,j,k)

            zag        = zgat*a1t(k  )*zpianf(i,j,k  )/zpia(i,j,k-1)        &
              &          - z2ag
            zcg        = zgct*a1t(k+1)*zpianf(i,j,k+1)/zpia(i,j,k+1)        &
              &          - z2cg
            zbg        = zzdtr - zgat*a1t(k  )*zpianf(i,j,k  )/zpia(i,j,k)  &
              &                - zgct*a1t(k+1)*zpianf(i,j,k+1)/zpia(i,j,k)  &
              &                + z2ag + z2cg

            zdg        = zzdtr * t(i,j,k,nnow) + ttens(i,j,k)
            zdg        = zdg - zgat * zpianf(i,j,k) *  (                       &
              &                   a2t(k)*  ( ztheta(i,j,k-1) - ztheta(i,j,k) ) &
              &                - zlhvocp * ( qc(i,j,k-1)/zpia(i,j,k-1)         &
              &                             - qc(i,j,k  )/zpia(i,j,k  ) ) )
            zdg        = zdg - zgct * zpianf(i,j,k+1)* (                       &
              &                   a2t(k+1)*( ztheta(i,j,k+1) - ztheta(i,j,k) ) &
              &                - zlhvocp * ( qc(i,j,k+1)/zpia(i,j,k+1)         &
              &                            - qc(i,j,k  )/zpia(i,j,k  ) ) )

            zz         = 1.0_wp/( zbg - zag*zc(i,j,k-1) )
            zc (i,j,k) = zcg*zz
            zd1(i,j,k) = ( zdg -zag*zd1(i,j,k-1) )*zz
          ENDDO
        ENDDO
      ENDDO

      ! The bottom layer
      DO j = jstart, jend
        DO i = istart, iend
          zgat       = - zkh (i,j,ke)*zsqrtgrho_r_s(i,j,ke)
          zgct       = - ztch(i,j   )*zsqrtgrho_r_s(i,j,ke)
          zag        = zgat*a1t(ke)*zpianf(i,j,ke)/zpia(i,j,ke-1)
          zbg        = zzdtr - zgat*a1t(ke )*zpianf(i,j,ke )/zpia(i,j,ke)   &
            &          - zgct*za1t_surf*zpianf(i,j,ke1)/zpia(i,j,ke)
          zdg        = zzdtr * t(i,j,ke,nnow) + ttens(i,j,ke)
          zdg        = zdg - zgat*zpianf(i,j,ke ) * (                        &
            &    a2t(ke)* ( ztheta(i,j,ke-1) - ztheta(i,j,ke) )              &
            & - zlhvocp * ( qc(i,j,ke-1)/zpia(i,j,ke-1)                      &
            & - qc(i,j,ke  )/zpia(i,j,ke  ) ) )
          zdg        = zdg + zgct*zpianf(i,j,ke1)* (                         &
            &   za2t_surf * ztheta(i,j,ke)                                   &
            & - zlhvocp *   qc(i,j,ke)/zpia(i,j,ke  ) )                      &
            & - zgct*t_g(i,j,nnow)
          znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
          ttens(i,j,ke) = ( znew - t(i,j,ke,nnow) ) * zzdtr
          ze   (i,j,ke) = znew
        ENDDO
      ENDDO

    ELSE

      ! Top layer
      DO j = jstart, jend
        DO i = istart, iend
          zgct       = - zkh(i,j,2)*zsqrtgrho_r_s(i,j,1)
          zcg        = zgct*a1t(2)*zpianf(i,j,2)/zpia(i,j,2)
          zbg        = zzdtr - zgct*a1t(2)*zpianf(i,j,2)/zpia(i,j,1)
          zdg        = zzdtr * t(i,j,1,nnow) + ttens(i,j,1)            &
            &        - zgct* a2t(2)*zpianf(i,j,2) *                    &
            &            ( ztheta(i,j,2) - ztheta(i,j,1) )
          zd1(i,j,1) = zdg/zbg
          zc(i,j,1)  = zcg/zbg
        ENDDO
      ENDDO

      ! The layers from k=2 to k=ke-1
      DO k = 2, ke-1
        DO j = jstart, jend
          DO i = istart, iend

            rho_inv = 1.0_wp / rho(i,j,k)

            ! the explicit counterparts of z2ag, z2cg are in ttens (via 'explicit_horizontal_diffusion'):
            z2ag = rho_inv * Kvh_rho_pi_at_w(i,j,k)   * beta_1(k) *                                   &
              &     ( acrlat(j,1)   * dzeta_dlam(i,j,k) * acrlat(j,1)   * dzeta_dlam_at_w(i,j,k)      &
              &     + r_earth_recip * dzeta_dphi(i,j,k) * r_earth_recip * dzeta_dphi_at_w(i,j,k) )

            z2cg = rho_inv * Kvh_rho_pi_at_w(i,j,k+1) * beta_1(k+1) *                                 &
              &     ( acrlat(j,1)   * dzeta_dlam(i,j,k) * acrlat(j,1)   * dzeta_dlam_at_w(i,j,k+1)    &
              &     + r_earth_recip * dzeta_dphi(i,j,k) * r_earth_recip * dzeta_dphi_at_w(i,j,k+1) )

            zgat = - zkh(i,j,k  ) * zsqrtgrho_r_s(i,j,k)
            zgct = - zkh(i,j,k+1) * zsqrtgrho_r_s(i,j,k)

            zag  = zgat*a1t(k  )*zpianf(i,j,k  )/zpia(i,j,k-1)        &
              &   - z2ag
            zcg  = zgct*a1t(k+1)*zpianf(i,j,k+1)/zpia(i,j,k+1)        &
              &   - z2cg
            zbg  = zzdtr - zgat*a1t(k  )*zpianf(i,j,k  )/zpia(i,j,k)  &
              &          - zgct*a1t(k+1)*zpianf(i,j,k+1)/zpia(i,j,k)  &
              &          + z2ag + z2cg

            zdg  = zzdtr * t(i,j,k,nnow) + ttens(i,j,k)                &
              &    - zgat * a2t(k)   * zpianf(i,j,k)   *               &
              &         ( ztheta(i,j,k-1) - ztheta(i,j,k) )            &
              &    - zgct * a2t(k+1) * zpianf(i,j,k+1) *               &
              &         ( ztheta(i,j,k+1) - ztheta(i,j,k) )

            zz         = 1.0_wp/( zbg - zag*zc(i,j,k-1) )
            zc (i,j,k) = zcg*zz
            zd1(i,j,k) = ( zdg -zag*zd1(i,j,k-1) )*zz
          ENDDO
        ENDDO
      ENDDO

      ! The bottom layer
      DO j = jstart, jend
        DO i = istart, iend
          zgat       = - zkh (i,j,ke)*zsqrtgrho_r_s(i,j,ke)
          zgct       = - ztch(i,j   )*zsqrtgrho_r_s(i,j,ke)

          zag        = zgat*a1t(ke)*zpianf(i,j,ke)/zpia(i,j,ke-1)
          zbg        = zzdtr - zgat* a1t(ke )*zpianf(i,j,ke )/zpia(i,j,ke)   &
            &                - zgct*za1t_surf*zpianf(i,j,ke1)/zpia(i,j,ke)
          zdg        = zzdtr * t(i,j,ke,nnow) + ttens(i,j,ke)                &
            &         - zgat* a2t(ke)  *zpianf(i,j,ke ) *                    &
            &                   ( ztheta(i,j,ke-1) - ztheta(i,j,ke) )        &
            &         + zgct*za2t_surf *zpianf(i,j,ke1)*                     &
            &                     ztheta(i,j,ke)                             &
            &         - zgct*t_g(i,j,nnow)
          znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
          ttens(i,j,ke) = ( znew - t(i,j,ke,nnow) ) * zzdtr
          ze   (i,j,ke) = znew
        ENDDO
      ENDDO

    END IF

    ! Backsubstitution and storage of the complete slow tendencies

    DO k = ke-1, 1, -1
      DO j = jstart, jend
        DO i = istart, iend
          ze   (i,j,k) = zd1(i,j,k) - zc(i,j,k)*ze(i,j,k+1)
          ttens(i,j,k) = ( ze(i,j,k) - t(i,j,k,nnow) ) * zzdtr
        ENDDO
      ENDDO
    ENDDO

    ! Calculation of the sensible heat flux at the surface
    ! This flux is integrated in time

    DO j = jstart, jend
      DO  i  = istart , iend
        shfl_s (i,j) = - ztch(i,j)*cp_d*                  &
          ( za2t_surf*(t_g(i,j,nnow) - zpianf(i,j,ke1)*   &
          t  (i,j,ke,nnow)/zpia(i,j,ke) ) +  &
          za1t_surf*(t_g(i,j,nnew) - zpianf(i,j,ke1)*   &
          ze(i,j,ke)/zpia(i,j,ke) ) )
      ENDDO
    ENDDO

    !--------------------------------------------------------------------------
    ! Section : Horizontal wind velocity u 
    !--------------------------------------------------------------------------

    ! Top layer:   k = 1
    DO j = jstartu, jendu
      DO i = istartu, iendu
        zkm(i,j)   = 0.5_wp*( ztmkvm(i,j,2) + ztmkvm(i+1,j,2) )
        zgct       = - zkm(i,j)*zsqrtgrho_r_u(i,j,1)
        zcg        = zgct*a1t(2)
        zbg        = zzdtr - zcg
        zdg        = zzdtr * u(i,j,1,nnow) + utens(i,j,1)    &
          - a2t(2)*zgct * ( u(i,j,2,nnow) - u(i,j,1,nnow) )
        zd1(i,j,1) = zdg/zbg
        zc (i,j,1) = zcg/zbg
      ENDDO
    ENDDO

    ! The layers from k=2 to k=ke-1
    DO k = 2, ke-1
      DO j = jstartu, jendu
        DO i = istartu, iendu

          rho_inv_at_u = 2.0_wp / ( rho(i,j,k) + rho(i+1,j,k) )

          dzeta_dphi_at_u = 0.5_wp * ( dzeta_dphi(i,j,k) + dzeta_dphi(i+1,j,k) )

          ! the explicit counterparts of z2ag, z2cg are in 'explicit_horizontal_diffusion':
          z2ag = rho_inv_at_u * Kvm_rho_at_uw(i,j,k)                                 &
            &  * ( acrlat(j,1)   * dzeta_dlam_at_u(i,j,k)                            &
            &      * beta_uu1(k) * 2.0_wp * acrlat(j,1) * dzeta_dlam_at_uw(i,j,k)    &
            &    + r_earth_recip * dzeta_dphi_at_u                                   &
            &      * beta_uu3(k) *        r_earth_recip * dzeta_dphi_at_uw(i,j,k) )

          z2cg = rho_inv_at_u * Kvm_rho_at_uw(i,j,k+1)                                 &
            &  * ( acrlat(j,1)   * dzeta_dlam_at_u(i,j,k)                              &
            &      * beta_uu1(k+1) * 2.0_wp * acrlat(j,1) * dzeta_dlam_at_uw(i,j,k+1)  &
            &    + r_earth_recip * dzeta_dphi_at_u                                     &
            &      * beta_uu3(k+1) *        r_earth_recip * dzeta_dphi_at_uw(i,j,k+1) )

          zgat       = - zkm(i,j)*zsqrtgrho_r_u(i,j,k)
          zkm(i,j)   = 0.5_wp*( ztmkvm(i,j,k+1) + ztmkvm(i+1,j,k+1) )
          zgct       = - zkm(i,j)*zsqrtgrho_r_u(i,j,k)

          zag        = zgat*a1t(k)   - z2ag
          zas        = zgat*a2t(k)
          zcg        = zgct*a1t(k+1) - z2cg
          zcs        = zgct*a2t(k+1)
          zbg        = zzdtr - zag - zcg
          zdg        = zzdtr * u(i,j,k,nnow) + utens(i,j,k)  &
            - zas * ( u(i,j,k-1,nnow) - u(i,j,k,nnow) )      &
            - zcs * ( u(i,j,k+1,nnow) - u(i,j,k,nnow) )
          zz         = 1.0_wp/( zbg - zag*zc(i,j,k-1) )
          zc (i,j,k)  = zcg*zz
          zd1(i,j,k)  = ( zdg -zag*zd1(i,j,k-1) )*zz
        ENDDO
      ENDDO
    ENDDO

    ! The bottom layer:  k = ke
    ! Including the calculation of the u-momentum flux at the surface
    DO j = jstartu, jendu
      DO i = istartu, iendu
        zgat       = - zkm(i,j)*zsqrtgrho_r_u(i,j,ke)
        ztmcmq     = 0.5_wp*( ztcm(i,j) + ztcm(i+1,j) )
        zgct       = - ztmcmq*zsqrtgrho_r_u(i,j,ke)
        zag        = zgat*a1t(ke)
        zas        = zgat*a2t(ke)
        zcg        = a1t(ke1)*zgct
        zcs        = a2t(ke1)*zgct
        zbg        = zzdtr - zag - zcg
        zdg        = zzdtr * u(i,j,ke,nnow) + utens(i,j,ke)  &
          - zas * ( u(i,j,ke-1,nnow) - u(i,j,ke,nnow) )      &
          + zcs * u(i,j,ke,nnow)
        znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
        utens(i,j,ke) = ( znew - u(i,j,ke,nnow) ) * zzdtr
        umfl_s(i,j  ) = ztmcmq*( a2t(ke1)*u(i,j,ke,nnow) + a1t(ke1)*znew )
        ze   (i,j,ke) = znew
      ENDDO
    ENDDO

    ! Backsubstitution and storage of the complete slow tendencies

    DO k = ke-1, 1, -1
      DO j = jstartu, jendu
        DO i = istartu, iendu
          ze   (i,j,k) = zd1(i,j,k) - zc(i,j,k)*ze(i,j,k+1)
          utens(i,j,k) = ( ze(i,j,k) - u(i,j,k,nnow) ) * zzdtr
        ENDDO
      ENDDO
    ENDDO

    !--------------------------------------------------------------------------
    ! Section : Horizontal wind velocity v
    !--------------------------------------------------------------------------

    ! Top layer  k=1
    DO j = jstartv, jendv
      DO i = istartv, iendv
        zkm(i,j)   = 0.5_wp*( ztmkvm(i,j,2) + ztmkvm(i,j+1,2) )
        zgct       = - zkm(i,j)*zsqrtgrho_r_v(i,j,1)
        zcg        = zgct*a1t(2)
        zbg        = zzdtr - zcg
        zdg        = zzdtr * v(i,j,1,nnow) + vtens(i,j,1)     &
          - a2t(2)*zgct * ( v(i,j,2,nnow) - v(i,j,1,nnow) )
        zd1(i,j,1) = zdg/zbg
        zc (i,j,1) = zcg/zbg
      ENDDO
    ENDDO

    ! The layers from k=2 to k=ke-1
    DO  k = 2, ke-1
      DO j = jstartv, jendv
        DO i = istartv, iendv

          rho_inv_at_v = 2.0_wp / ( rho(i,j,k) + rho(i,j+1,k) )

          dzeta_dlam_at_v = 0.5_wp * ( dzeta_dlam(i,j,k) + dzeta_dlam(i,j+1,k) )

          ! the explicit counterparts of z2ag, z2cg are in 'explicit_horizontal_diffusion':
          z2ag = rho_inv_at_v * Kvm_rho_at_vw(i,j,k)                                     &
            &  * ( acrlat(j,2)   * dzeta_dlam_at_v                                       &
            &      * beta_vv3(k) *          acrlat(j,2)   * dzeta_dlam_at_vw(i,j,k)      &
            &    + r_earth_recip * dzeta_dphi_at_v(i,j,k)                                &
            &      * beta_vv1(k) * 2.0_wp * r_earth_recip * dzeta_dphi_at_vw(i,j,k) )

          z2cg = rho_inv_at_v * Kvm_rho_at_vw(i,j,k+1)                                   &
            &  * ( acrlat(j,2)   * dzeta_dlam_at_v                                       &
            &      * beta_vv3(k+1) *          acrlat(j,2)   * dzeta_dlam_at_vw(i,j,k+1)  &
            &    + r_earth_recip * dzeta_dphi_at_v(i,j,k)                                &
            &      * beta_vv1(k+1) * 2.0_wp * r_earth_recip * dzeta_dphi_at_vw(i,j,k+1) )

          zgat       = - zkm(i,j)*zsqrtgrho_r_v(i,j,k)
          zkm(i,j)   = 0.5_wp*( ztmkvm(i,j,k+1) + ztmkvm(i,j+1,k+1) )
          zgct       = - zkm(i,j)*zsqrtgrho_r_v(i,j,k)

          zag        = zgat*a1t(k)   - z2ag
          zas        = zgat*a2t(k)
          zcg        = zgct*a1t(k+1) - z2cg
          zcs        = zgct*a2t(k+1)
          zbg        = zzdtr - zag - zcg
          zdg        = zzdtr * v(i,j,k,nnow) + vtens(i,j,k)  &
            - zas * ( v(i,j,k-1,nnow) - v(i,j,k,nnow) )      &
            - zcs * ( v(i,j,k+1,nnow) - v(i,j,k,nnow) )
          zz         = 1.0_wp/( zbg - zag*zc(i,j,k-1) )
          zc (i,j,k) = zcg*zz
          zd1(i,j,k) = ( zdg -zag*zd1(i,j,k-1) )*zz
        ENDDO
      ENDDO
    ENDDO

    ! The bottom layer k=ke
    ! Including the calculation of the v-momentum flux at the surface
    DO j = jstartv, jendv
      DO i = istartv, iendv
        zgat       = - zkm(i,j)*zsqrtgrho_r_v(i,j,ke)
        ztmcmq     = 0.5_wp*( ztcm(i,j) + ztcm(i,j+1) )
        zgct       = - ztmcmq*zsqrtgrho_r_v(i,j,ke)
        zag        = zgat*a1t(ke)
        zas        = zgat*a2t(ke)
        zcg        = a1t(ke1)*zgct
        zcs        = a2t(ke1)*zgct
        zbg        = zzdtr - zag - zcg
        zdg        = zzdtr * v(i,j,ke,nnow) + vtens(i,j,ke)  &
          - zas * ( v(i,j,ke-1,nnow) - v(i,j,ke,nnow) )      &
          + zcs * v(i,j,ke,nnow)
        znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
        vtens(i,j,ke) = ( znew - v(i,j,ke,nnow) ) * zzdtr
        vmfl_s(i,j  ) = ztmcmq*( a2t(ke1)*v(i,j,ke,nnow) + a1t(ke1)*znew )
        ze   (i,j,ke) = znew
      ENDDO
    ENDDO

    ! Backsubstitution and storage of the complete slow tendencies

    DO k = ke-1, 1, -1
      DO j = jstartv, jendv
        DO i = istartv, iendv
          ze   (i,j,k) = zd1(i,j,k) - zc(i,j,k)*ze(i,j,k+1)
          vtens(i,j,k) = ( ze(i,j,k) - v(i,j,k,nnow) ) * zzdtr
        ENDDO
      ENDDO
    ENDDO

  END IF

  !----------------------------------------------------------------------------
  ! Section : Vertical wind velocity w
  !----------------------------------------------------------------------------

  IF ( lvertdiff_w ) THEN

    DO j = jstart, jend

      DO  k = 2, ke
        km1 = MAX( 2, k-1 )
        kp1 = MIN( ke, k+1 )

        ! Precalculate some help variables to avoid divisions in subsequent code
        DO i = istart, iend
          zgatz(i,j,k) = - zsqrtgrho_r_w(i,j,k)   &
            * ( ztmkvw(i,j,k)+ztmkvw(i,j,km1) )
          zgctz(i,j,k) = - zsqrtgrho_r_w(i,j,k)   &
            * ( ztmkvw(i,j,kp1)+ztmkvw(i,j,k) )
        ENDDO
      ENDDO

      ! Top layer       
      DO i = istart, iend
        zag        = 0.5_wp*(a1t(2)+a1t(1)) * zgatz(i,j,2)
        zas        = 0.5_wp*(a2t(2)+a2t(1)) * zgatz(i,j,2)
        zcg        = 0.5_wp*(a1t(3)+a1t(2)) * zgctz(i,j,2)
        zcs        = 0.5_wp*(a2t(3)+a2t(2)) * zgctz(i,j,2)
        zbg        = zzdtr - zag - zcg
        zdg        = zzdtr * w(i,j,2,nnow) + wtens(i,j,2)         &
          + zas * w(i,j,2,nnow)                                   &
          - zcs * ( w(i,j,3,nnow) - w(i,j,2,nnow) )
        zd1(i,j,2)  = zdg/zbg
        zc (i,j,2)  = zcg/zbg
      ENDDO

      ! The layers from k=3 to k=ke-1
      DO k = 3, ke-1
        DO i = istart, iend

          rho_inv_at_w = 2.0_wp / ( rho(i,j,k) + rho(i,j,k-1) )

          ! the explicit counterparts of z2ag, z2cg are in 'explicit_horizontal_diffusion':
          z2ag = rho_inv_at_w * Kvm_rho_at_s(i,j,k) * beta_ww1(k)      &
            &  * (   acrlat(j,1)   * dzeta_dlam_at_w(i,j,k)            &
            &      * acrlat(j,1)   * dzeta_dlam(i,j,k)                 &
            &    +   r_earth_recip * dzeta_dphi_at_w(i,j,k)            &
            &      * r_earth_recip * dzeta_dphi(i,j,k) )

          z2cg = rho_inv_at_w * Kvm_rho_at_s(i,j,k+1) * beta_ww1(k+1)  &
            &  * (   acrlat(j,1)   * dzeta_dlam_at_w(i,j,k)            &
            &      * acrlat(j,1)   * dzeta_dlam(i,j,k+1)               &
            &    +   r_earth_recip * dzeta_dphi_at_w(i,j,k)            &
            &      * r_earth_recip * dzeta_dphi(i,j,k+1) )

          zag        = 0.5_wp*(a1t(k)+a1t(k-1)) * zgatz(i,j,k) - z2ag
          zas        = 0.5_wp*(a2t(k)+a2t(k-1)) * zgatz(i,j,k)
          zcg        = 0.5_wp*(a1t(k+1)+a1t(k)) * zgctz(i,j,k) - z2cg
          zcs        = 0.5_wp*(a2t(k+1)+a2t(k)) * zgctz(i,j,k)
          zbg        = zzdtr - zag - zcg
          zdg        = zzdtr * w(i,j,k,nnow) + wtens(i,j,k)         &
            - zas * ( w(i,j,k-1,nnow) - w(i,j,k,nnow) )             &
            - zcs * ( w(i,j,k+1,nnow) - w(i,j,k,nnow) )
          zz         = 1.0_wp/( zbg - zag*zc(i,j,k-1) )
          zc (i,j,k)  = zcg*zz
          zd1(i,j,k)  = ( zdg - zag*zd1(i,j,k-1) )*zz
        ENDDO
      ENDDO

      ! The bottom layer
      DO i = istart, iend
        zag        = 0.5_wp*(a1t(ke)+a1t(ke-1)) * zgatz(i,j,ke)
        zas        = 0.5_wp*(a2t(ke)+a2t(ke-1)) * zgatz(i,j,ke)
        zcg        = 0.5_wp*(a1t(ke1)+a1t(ke))  * zgctz(i,j,ke)
        zcs        = 0.5_wp*(a2t(ke1)+a2t(ke))  * zgctz(i,j,ke)
        zbg        = zzdtr - zag - zcg
        zdg        = zzdtr * w(i,j,ke,nnow) + wtens(i,j,ke)      &
          - zas * ( w(i,j,ke-1,nnow) - w(i,j,ke,nnow) )          &
          - zcs * ( w(i,j,ke1,nnow) - w(i,j,ke,nnow))            &
          - zcg * w(i,j,ke1,nnow) 
        znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
        wtens(i,j,ke) = ( znew - w(i,j,ke,nnow) ) * zzdtr
        ze   (i,j,ke) = znew
      ENDDO

      ! Backsubstitution and storage of the complete slow tendencies

      DO k = ke-1, 2, -1
        DO i = istart, iend
          ze    (i,j,k) = zd1(i,j,k) - zc(i,j,k)*ze(i,j,k+1)
          wtens (i,j,k) = ( ze(i,j,k) - w(i,j,k,nnow) ) * zzdtr
        ENDDO
      ENDDO
    ENDDO

  END IF


  DEALLOCATE( dzeta_dlam_at_vw )
  DEALLOCATE( dzeta_dphi_at_uw )

  !----------------------------------------------------------------------------
  ! End of subroutine implicit_vert_diffusion_uvwt
  !----------------------------------------------------------------------------

END SUBROUTINE implicit_vert_diffusion_uvwt_3D


!==============================================================================
!==============================================================================


SUBROUTINE complete_tend_uvwtpp_CN3Crow( dt_adv, nstar )

  !----------------------------------------------------------------------------
  !
  ! Description:
  !   This procedure calculates the tendencies of the vertical advection 
  !   (uadvt, vadvt, wadvt, ppadvt and tadvt)
  !   for the prognostic variables u,v,w,pp and T.
  !
  ! Method:
  !   the vertical advection is solved by a vertically implicit 
  !   scheme (modified Crank-Nicolson 3. order) starting at the states
  !   u(:,:,:,nstar), v(:,:,:,now), w(..,nstar), pp(..), T(..)
  !
  !----------------------------------------------------------------------------

  USE numeric_utilities, ONLY: solve_5banddiag, solve_5banddiag_vec

  USE data_parallel,      ONLY :  &
    my_cart_id     ! rank of this subdomain in the cartesian communicator


  IMPLICIT NONE

  ! Declarations:

  ! Subroutine arguments:
  ! ---------------------

  INTEGER (KIND=iintegers), INTENT(in) ::  &
    nstar

  REAL (KIND=wp),     INTENT(in) ::  &
    dt_adv


  ! Local scalars:
  ! -------------
  INTEGER (KIND=iintegers) ::  &
    i,  j,  k              !  Loop indices in lon., lat. and vert. direction

  REAL    (KIND=wp)        ::  &
    z_beta_v, z_beta_p, z_beta_m,     &  ! Crank-Nicholson-weights
    zC2, &           ! square of local Courant number
    zh1, zh2, &      ! help variable for storage of intermediate results
    z_dt_half, &     ! = dt_adv / 2
    z_dt_quart, &    ! = dt_adv / 4
    z_dt_recip, &    ! = 1 / dt
    z_beta_p_half, & ! = z_beta_p / 2
    z_beta_m_half    ! = z_beta_m / 2


  ! Local (automatic) arrays:
  ! ------------------------

  ! Arrays for the LGS with a  5-band-diagonal-matrix:
  REAL    (KIND=wp),    ALLOCATABLE ::  &
    z_Cour    (:,:,:),   &
    z_lgs     (:,:,:,:), &  ! the 5 banddiagonals of the lin. eq. system
    z_lgs_rhs (:,:,:),   &
    zh        (:,:,:)

  REAL    (KIND=wp),    ALLOCATABLE ::  &
    z_lgs_store (:,:,:,:) ! copy of the 5 banddiagonals of the lin. eq. system

  ! nicht sehr elegant, dass man hier extra Felder deklarieren muss:
  REAL    (KIND=wp),    ALLOCATABLE ::   &
    z_Cour_w    (:,:,:),   &
    z_lgs_w     (:,:,:,:), & ! the 5 banddiagonals of the lin. eq. system
    z_lgs_rhs_w (:,:,:),   &
    zh_w        (:,:,:)

  INTEGER (KIND=iintegers) ::  istat

  LOGICAL :: flag_vector_version  ! call Numerical Recipes routines, 
                                  ! which are optimized for vector computers

  ! End of header
  !============================================================================

  !----------------------------------------------------------------------------
  ! Begin Subroutine complete_tend_uvwtpp_CN3Crow
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Section 1: Some preparations
  !----------------------------------------------------------------------------

  flag_vector_version = .TRUE.

  ALLOCATE( z_Cour     (1:ie, 1:je, 1:ke), STAT=istat )
  ALLOCATE( z_lgs      (1:ie, 1:je, 1:ke, 1:5), STAT=istat )
  ALLOCATE( z_lgs_rhs  (1:ie, 1:je, 1:ke), STAT=istat )
  ALLOCATE( zh         (1:ie, 1:je, 1:ke), STAT=istat )

  ALLOCATE( z_Cour_w   (1:ie, 1:je, 1:ke1), STAT=istat )
  ALLOCATE( z_lgs_w    (1:ie, 1:je, 1:ke1, 1:5), STAT=istat )
  ALLOCATE( z_lgs_rhs_w(1:ie, 1:je, 1:ke1), STAT=istat )
  ALLOCATE( zh_w       (1:ie, 1:je, 1:ke1), STAT=istat )

  IF ( izdebug >= 10 ) THEN
    WRITE(*,*) "subr. complete_tend_uvwtpp_CN3Crow ..."
    WRITE(*,*) "  flag_vector_version = ", flag_vector_version
  END IF


  !MB: for debugging:
  !z_Cour     = 0.0
  !z_lgs      = 0.0
  !z_lgs_rhs  = 0.0
  !zh         = 0.0

  !z_Cour_w     = 0.0
  !z_lgs_w      = 0.0
  !z_lgs_rhs_w  = 0.0
  !zh_w         = 0.0

  ! precalculated factors for optimization 
  z_dt_half  = 0.5_wp  * dt_adv 
  z_dt_quart = 0.25_wp * dt_adv 
  z_dt_recip = 1.0_wp / dt_adv

  ! Setting of parameters for implicit calculation of vertical advection       
  ! (z_beta_v: weight for the n+1 time level           
  !          z_beta_v=0: centered, =+1: full implicit,=-1: explicit)
  z_beta_v = 0.0_wp
  z_beta_p = 0.5_wp*( 1.0_wp + z_beta_v )
  z_beta_m = 0.5_wp*( 1.0_wp - z_beta_v )

  z_beta_p_half = 0.5_wp * z_beta_p
  z_beta_m_half = 0.5_wp * z_beta_m

  !----------------------------------------------------------------------------
  ! Section 2: Setup of tridiagonal matrix systems resulting from the implicit
  !            numerical formulation of advection.
  !            - variables at scalar positions (T and pp)
  !----------------------------------------------------------------------------

  ! Band-diagonal-matrix

  ! Top layer k=1 (one-sided diff.)  
  DO j = jstart, jend
    DO i = istart, iend

      ! remark: wcon contains 'contravar. vertical velocity * sqrt(G)':
      z_Cour(i,j,1) = sqrtg_r_s(i,j,1) * ( wcon(i,j,1) + wcon(i,j,2) ) * z_dt_half 

      IF ( z_Cour(i,j,1) < 0.0_wp ) THEN
      z_lgs(i,j,1,1) = 0.0_wp
      z_lgs(i,j,1,2) = 0.0_wp
      z_lgs(i,j,1,3) = 1.0_wp - z_Cour(i,j,1) * z_beta_p
      z_lgs(i,j,1,4) =            + z_Cour(i,j,1) * z_beta_p
      z_lgs(i,j,1,5) = 0.0_wp
      ELSE
        ! use explicit upwind:
        z_lgs(i,j,1,1) = 0.0_wp
        z_lgs(i,j,1,2) = 0.0_wp
        z_lgs(i,j,1,3) = 1.0_wp
        z_lgs(i,j,1,4) = 0.0_wp
        z_lgs(i,j,1,5) = 0.0_wp
      END IF

    ENDDO
  ENDDO

  ! 2. Layer k=2 (centered diff. 2. order)
  DO j = jstart, jend
    DO i = istart, iend

      z_Cour(i,j,2) = sqrtg_r_s(i,j,2) * ( wcon(i,j,2) + wcon(i,j,3) ) * z_dt_half 

      z_lgs(i,j,2,1) = 0.0_wp
      z_lgs(i,j,2,2) = - z_Cour(i,j,2) * z_beta_p_half
      z_lgs(i,j,2,3) = 1.0_wp
      z_lgs(i,j,2,4) = - z_lgs(i,j,2,2)
      z_lgs(i,j,2,5) = 0.0_wp

    ENDDO
  ENDDO

  ! The layers from k=3 to k=ke-2
  DO k = 3, ke-2
    DO j = jstart, jend
      DO i = istart, iend

        z_Cour(i,j,k) = sqrtg_r_s(i,j,k) * ( wcon(i,j,k) + wcon(i,j,k+1) ) * z_dt_half 
        zC2 = z_Cour(i,j,k) * z_Cour(i,j,k)
        zh2 = z_Cour(i,j,k) * z_beta_p

        IF ( z_Cour(i,j,k) > 0 ) THEN
          z_lgs(i,j,k,1) = zh2 * ( 1.0_wp/6.0_wp + zC2/12.0_wp )
          z_lgs(i,j,k,2) = zh2 * ( - 1.0_wp      - zC2/4.0_wp  )
          z_lgs(i,j,k,3) = zh2 * ( 1.0_wp/2.0_wp + zC2/4.0_wp  ) + 1.0_wp
          z_lgs(i,j,k,4) = zh2 * ( 1.0_wp/3.0_wp - zC2/12.0_wp )
          z_lgs(i,j,k,5) = 0.0_wp
        ELSE 
          z_lgs(i,j,k,1) = 0.0_wp
          z_lgs(i,j,k,2) = - zh2 * ( 1.0_wp/3.0_wp - zC2/12.0_wp )
          z_lgs(i,j,k,3) = - zh2 * ( 1.0_wp/2.0_wp + zC2/4.0_wp  ) + 1.0_wp
          z_lgs(i,j,k,4) = - zh2 * ( - 1.0_wp      - zC2/4.0_wp  )
          z_lgs(i,j,k,5) = - zh2 * ( 1.0_wp/6.0_wp + zC2/12.0_wp )
        END IF

      ENDDO
    ENDDO
  ENDDO

  ! Layer k=ke-1 (centered diff. 2. order)
  DO j = jstart, jend
    DO i = istart, iend

      z_Cour(i,j,ke-1) = sqrtg_r_s(i,j,ke-1) * ( wcon(i,j,ke-1) + wcon(i,j,ke) ) * z_dt_half 

      z_lgs(i,j,ke-1,1) = 0.0_wp
      z_lgs(i,j,ke-1,2) = - z_Cour(i,j,ke-1) * z_beta_p_half
      z_lgs(i,j,ke-1,3) = 1.0_wp
      z_lgs(i,j,ke-1,4) = - z_lgs(i,j,ke-1,2)
      z_lgs(i,j,ke-1,5) = 0.0_wp

    ENDDO
  ENDDO

  ! The bottom layer k=ke
  DO j = jstart, jend
    DO i = istart, iend

      z_Cour(i,j,ke) = sqrtg_r_s(i,j,ke) * ( wcon(i,j,ke) + wcon(i,j,ke+1) ) * z_dt_half  

      IF ( z_Cour(i,j,k) > 0 ) THEN
        z_lgs(i,j,ke,1) = 0.0_wp
        z_lgs(i,j,ke,2) =            - z_Cour(i,j,ke) * z_beta_p
        z_lgs(i,j,ke,3) = 1.0_wp + z_Cour(i,j,ke) * z_beta_p
        z_lgs(i,j,ke,4) = 0.0_wp
        z_lgs(i,j,ke,5) = 0.0_wp
      ELSE
        ! use explicit upwind:
        z_lgs(i,j,ke,1) = 0.0_wp
        z_lgs(i,j,ke,2) = 0.0_wp
        z_lgs(i,j,ke,3) = 1.0_wp
        z_lgs(i,j,ke,4) = 0.0_wp
        z_lgs(i,j,ke,5) = 0.0_wp
      END IF
    ENDDO
  ENDDO

  ! --- right hand side for T ---

  ! Top layer (k=1)
  DO j = jstart, jend
    DO i = istart, iend
      IF ( z_Cour(i,j,1) < 0.0_wp ) THEN
      z_lgs_rhs(i,j,1) = T(i,j,1,nstar)                                   &
        &  - z_Cour(i,j,1) * z_beta_m * ( T(i,j,2,nstar)-T(i,j,1,nstar) ) 
      ELSE
        ! use explicit upwind:
        z_lgs_rhs(i,j,1) = T(i,j,1,nstar)                                   &
          &  - z_Cour(i,j,1) * 1.0_wp * ( T(i,j,2,nstar)-T(i,j,1,nstar) ) 
      END IF
    ENDDO
  ENDDO

  ! 2. layer k=2
  DO j = jstart, jend
    DO i = istart, iend
      z_lgs_rhs(i,j,2) = T(i,j,2,nstar)                                          &
        &  - z_Cour(i,j,2) * ( T(i,j,3,nstar) - T(i,j,1,nstar) ) * z_beta_m_half 

    ENDDO
  ENDDO

  ! The layers from k=3 to k=ke-2
  DO k = 3, ke-2
    DO j = jstart, jend
      DO i = istart, iend

        zC2 = z_Cour(i,j,k) * z_Cour(i,j,k) 

        ! store advection-operator in help-variable zh1:
        IF ( z_Cour(i,j,k) > 0 ) THEN
          zh1 =   ( 1.0_wp/6.0_wp + zC2/12.0_wp ) * T(i,j,k-2,nstar)  &
              & + ( - 1.0_wp      - zC2/4.0_wp  ) * T(i,j,k-1,nstar)  &
              & + ( 1.0_wp/2.0_wp + zC2/4.0_wp  ) * T(i,j,k  ,nstar)  &
              & + ( 1.0_wp/3.0_wp - zC2/12.0_wp ) * T(i,j,k+1,nstar)
        ELSE 
          zh1 = - ( 1.0_wp/6.0_wp + zC2/12.0_wp ) * T(i,j,k+2,nstar)  &
              & - ( - 1.0_wp      - zC2/4.0_wp  ) * T(i,j,k+1,nstar)  &
              & - ( 1.0_wp/2.0_wp + zC2/4.0_wp  ) * T(i,j,k  ,nstar)  &
              & - ( 1.0_wp/3.0_wp - zC2/12.0_wp ) * T(i,j,k-1,nstar)
        END IF
        z_lgs_rhs(i,j,k) = T(i,j,k,nstar)                  &
          &   - z_Cour(i,j,k) * zh1 * z_beta_m  

      ENDDO
    ENDDO
  ENDDO

  ! layer k=ke-1
  DO j = jstart, jend
    DO i = istart, iend
      z_lgs_rhs(i,j,ke-1) = T(i,j,ke-1,nstar)                                          &
        &   - z_Cour(i,j,ke-1) * ( T(i,j,ke,nstar)-T(i,j,ke-2,nstar) ) * z_beta_m_half
    ENDDO
  ENDDO

  ! Bottom layer k=ke
  DO j = jstart, jend
    DO i = istart, iend
      IF ( z_Cour(i,j,k) > 0 ) THEN
        z_lgs_rhs(i,j,ke) = T(i,j,ke,nstar)              &
          &  - z_Cour(i,j,ke) * z_beta_m * ( T(i,j,ke,nstar)-T(i,j,ke-1,nstar) )
      ELSE
        ! use explicit upwind:
        z_lgs_rhs(i,j,ke) = T(i,j,ke,nstar)              &
          &  - z_Cour(i,j,ke) * 1.0_wp * ( T(i,j,ke,nstar)-T(i,j,ke-1,nstar) )
      END IF
    ENDDO
  ENDDO

  IF ( flag_vector_version ) THEN 
    ALLOCATE( z_lgs_store( 1:ie, 1:je, 1:ke, 1:5 ), STAT=istat )

    z_lgs_store(:,:,:,:) = z_lgs(:,:,:,:)

    CALL solve_5banddiag_vec( z_lgs_store, z_lgs_rhs, zh, &
         ie, je, istart, iend, jstart, jend, ke ) 
    ! (remark: z_lgs_store now has changed)

    DEALLOCATE( z_lgs_store )
  ELSE

    CALL solve_5banddiag( z_lgs, z_lgs_rhs, zh, &
         ie, je, istart, iend, jstart, jend, ke ) 
  END IF

  DO k = 1, ke
    DO j = jstart, jend
      DO i = istart, iend
        tadvt(i,j,k) = ( zh(i,j,k) - T(i,j,k, nstar) ) * z_dt_recip  
      ENDDO
    ENDDO
  ENDDO

  ! --- right hand side for pp ---

  ! Top layer (k=1)      
  DO j = jstart, jend
    DO i = istart, iend
      IF ( z_Cour(i,j,1) < 0.0_wp ) THEN
        z_lgs_rhs(i,j,1) = pp(i,j,1,nstar)                                   &
          &  - z_Cour(i,j,1) * z_beta_m * ( pp(i,j,2,nstar)-pp(i,j,1,nstar) )
      ELSE
        ! use explicit upwind:
        z_lgs_rhs(i,j,1) = pp(i,j,1,nstar)                                   &
          &  - z_Cour(i,j,1) * 1.0_wp * ( pp(i,j,2,nstar)-pp(i,j,1,nstar) )
      END IF
    ENDDO
  ENDDO

  ! 2. layer k=2
  DO j = jstart, jend
    DO i = istart, iend
      z_lgs_rhs(i,j,2) = pp(i,j,2,nstar)                                           &
        &  - z_Cour(i,j,2) * ( pp(i,j,3,nstar) - pp(i,j,1,nstar) ) * z_beta_m_half

    ENDDO
  ENDDO

  ! The layers from k=3 to k=ke-2
  DO k = 3, ke-2
    DO j = jstart, jend
      DO i = istart, iend

        zC2 = z_Cour(i,j,k) * z_Cour(i,j,k) 

        ! store advection-operator in help-variable zh1:
        IF ( z_Cour(i,j,k) > 0 ) THEN
          zh1 =   ( 1.0_wp/6.0_wp + zC2/12.0_wp ) * pp(i,j,k-2,nstar)  &
              & + ( - 1.0_wp      - zC2/4.0_wp  ) * pp(i,j,k-1,nstar)  &
              & + ( 1.0_wp/2.0_wp + zC2/4.0_wp  ) * pp(i,j,k  ,nstar)  &
              & + ( 1.0_wp/3.0_wp - zC2/12.0_wp ) * pp(i,j,k+1,nstar)
        ELSE 
          zh1 = - ( 1.0_wp/6.0_wp + zC2/12.0_wp ) * pp(i,j,k+2,nstar)  &
              & - ( - 1.0_wp      - zC2/4.0_wp  ) * pp(i,j,k+1,nstar)  &
              & - ( 1.0_wp/2.0_wp + zC2/4.0_wp  ) * pp(i,j,k  ,nstar)  &
              & - ( 1.0_wp/3.0_wp - zC2/12.0_wp ) * pp(i,j,k-1,nstar)
        END IF
        z_lgs_rhs(i,j,k) = pp(i,j,k,nstar)                  &
          &   - z_Cour(i,j,k) * zh1 * z_beta_m 

      ENDDO
    ENDDO
  ENDDO

  ! layer k=ke-1
  DO j = jstart, jend
    DO i = istart, iend
      z_lgs_rhs(i,j,ke-1) = pp(i,j,ke-1,nstar)          &
        &   - z_Cour(i,j,ke-1) * ( pp(i,j,ke,nstar)-pp(i,j,ke-2,nstar) ) * z_beta_m_half
    ENDDO
  ENDDO

  ! Bottom layer k=ke
  DO j = jstart, jend
    DO i = istart, iend
      IF ( z_Cour(i,j,k) > 0 ) THEN
        z_lgs_rhs(i,j,ke) = pp(i,j,ke,nstar)              &
          &  - z_Cour(i,j,ke) * z_beta_m * ( pp(i,j,ke,nstar)-pp(i,j,ke-1,nstar) )
      ELSE
        ! use explicit upwind:
        z_lgs_rhs(i,j,ke) = pp(i,j,ke,nstar)              &
          &  - z_Cour(i,j,ke) * 1.0_wp * ( pp(i,j,ke,nstar)-pp(i,j,ke-1,nstar) )
      END IF
    ENDDO
  ENDDO

  IF ( flag_vector_version ) THEN 
    CALL solve_5banddiag_vec( z_lgs, z_lgs_rhs, zh, &
        ie, je, istart, iend, jstart, jend, ke ) 
    ! attention: z_lgs has now changed

  ELSE
    CALL solve_5banddiag( z_lgs, z_lgs_rhs, zh, &
        ie, je, istart, iend, jstart, jend, ke ) 
  END IF

  DO k = 1, ke
    DO j = jstart, jend
      DO i = istart, iend
        ppadvt(i,j,k) = ( zh(i,j,k) - pp(i,j,k, nstar) ) * z_dt_recip  
      ENDDO
    ENDDO
  ENDDO

  !----------------------------------------------------------------------------
  ! Section 3: Setup of tridiagonal matrix systems resulting from the implicit
  !            numerical formulation of advection.
  !            -  variables at w-positions (w)
  !----------------------------------------------------------------------------


  ! Band-diagonal-matrix

  ! Top layer k=1 (one-sided diff.)  
  DO j = jstart, jend
    DO i = istart, iend

      z_Cour_w(i,j,1) = sqrtg_r_w(i,j,1) * wcon(i,j,1) * dt_adv 

      IF ( z_Cour(i,j,1) < 0.0_wp ) THEN
        z_lgs_w(i,j,1,1) = 0.0_wp
        z_lgs_w(i,j,1,2) = 0.0_wp
        z_lgs_w(i,j,1,3) = 1.0_wp - z_Cour_w(i,j,1) * z_beta_p
        z_lgs_w(i,j,1,4) =            + z_Cour_w(i,j,1) * z_beta_p
        z_lgs_w(i,j,1,5) = 0.0_wp
      ELSE
        ! use explicit upwind
        z_lgs_w(i,j,1,1) = 0.0_wp
        z_lgs_w(i,j,1,2) = 0.0_wp
        z_lgs_w(i,j,1,3) = 1.0_wp
        z_lgs_w(i,j,1,4) = 0.0_wp
        z_lgs_w(i,j,1,5) = 0.0_wp
      END IF

    ENDDO
  ENDDO

  ! 2. Layer k=2 (centered diff. 2. order)
  DO j = jstart, jend
    DO i = istart, iend

      z_Cour_w(i,j,2) = sqrtg_r_w(i,j,2) * wcon(i,j,2) * dt_adv 

      z_lgs_w(i,j,2,1) = 0.0_wp
      z_lgs_w(i,j,2,2) = - z_Cour_w(i,j,2) * z_beta_p_half
      z_lgs_w(i,j,2,3) = 1.0_wp
      z_lgs_w(i,j,2,4) = - z_lgs_w(i,j,2,2)
      z_lgs_w(i,j,2,5) = 0.0_wp

    ENDDO
  ENDDO

  ! The layers from k=3 to k=ke-1
  DO k = 3, ke-1
    DO j = jstart, jend
      DO i = istart, iend

        z_Cour_w(i,j,k) = sqrtg_r_w(i,j,k) * wcon(i,j,k) * dt_adv 
        zC2 = z_Cour_w(i,j,k) * z_Cour_w(i,j,k)
        zh2 = z_Cour_w(i,j,k) * z_beta_p

        IF ( z_Cour_w(i,j,k) > 0 ) THEN
          z_lgs_w(i,j,k,1) = zh2 * ( 1.0_wp/6.0_wp + zC2/12.0_wp )
          z_lgs_w(i,j,k,2) = zh2 * ( - 1.0_wp      - zC2/4.0_wp  )
          z_lgs_w(i,j,k,3) = zh2 * ( 1.0_wp/2.0_wp + zC2/4.0_wp  ) + 1.0_wp
          z_lgs_w(i,j,k,4) = zh2 * ( 1.0_wp/3.0_wp - zC2/12.0_wp )
          z_lgs_w(i,j,k,5) = 0.0_wp
        ELSE 
          z_lgs_w(i,j,k,1) = 0.0_wp
          z_lgs_w(i,j,k,2) = - zh2 * ( 1.0_wp/3.0_wp - zC2/12.0_wp )
          z_lgs_w(i,j,k,3) = - zh2 * ( 1.0_wp/2.0_wp + zC2/4.0_wp  ) + 1.0_wp
          z_lgs_w(i,j,k,4) = - zh2 * ( - 1.0_wp      - zC2/4.0_wp  )
          z_lgs_w(i,j,k,5) = - zh2 * ( 1.0_wp/6.0_wp + zC2/12.0_wp ) 
        END IF

      ENDDO
    ENDDO
  ENDDO

  ! Layer k=ke (centered diff. 2. order)
  DO j = jstart, jend
    DO i = istart, iend

      z_Cour_w(i,j,ke) = sqrtg_r_w(i,j,ke) * wcon(i,j,ke) * dt_adv 

      z_lgs_w(i,j,ke,1) = 0.0_wp
      z_lgs_w(i,j,ke,2) = - z_Cour_w(i,j,ke) * z_beta_p_half
      z_lgs_w(i,j,ke,3) = 1.0_wp
      z_lgs_w(i,j,ke,4) = - z_lgs_w(i,j,ke,2)
      z_lgs_w(i,j,ke,5) = 0.0_wp

    ENDDO
  ENDDO

  ! The bottom layer k=ke+1
  DO j = jstart, jend
    DO i = istart, iend

      z_Cour_w(i,j,ke1) = sqrtg_r_w(i,j,ke1) * wcon(i,j,ke1) * dt_adv  

      IF ( z_Cour(i,j,k) > 0 ) THEN
        z_lgs_w(i,j,ke1,1) = 0.0_wp
        z_lgs_w(i,j,ke1,2) =            - z_Cour_w(i,j,ke1) * z_beta_p
        z_lgs_w(i,j,ke1,3) = 1.0_wp + z_Cour_w(i,j,ke1) * z_beta_p
        z_lgs_w(i,j,ke1,4) = 0.0_wp
        z_lgs_w(i,j,ke1,5) = 0.0_wp
      ELSE
        z_lgs_w(i,j,ke1,1) = 0.0_wp
        z_lgs_w(i,j,ke1,2) = 0.0_wp
        z_lgs_w(i,j,ke1,3) = 1.0_wp
        z_lgs_w(i,j,ke1,4) = 0.0_wp
        z_lgs_w(i,j,ke1,5) = 0.0_wp
      END IF
    ENDDO
  ENDDO

  ! --- right hand side for w ---

  ! Top layer (k=1)
  DO j = jstart, jend
    DO i = istart, iend
      IF ( z_Cour(i,j,1) < 0.0_wp ) THEN
        z_lgs_rhs_w(i,j,1) = w(i,j,1,nstar)                                   &
          &  - z_Cour_w(i,j,1) * z_beta_m * ( w(i,j,2,nstar)-w(i,j,1,nstar) )
      ELSE
        ! use explicit upwind:
        z_lgs_rhs_w(i,j,1) = w(i,j,1,nstar)                                   &
          &  - z_Cour_w(i,j,1) * 1.0_wp * ( w(i,j,2,nstar)-w(i,j,1,nstar) )
      END IF
    ENDDO
  ENDDO

  ! 2. layer k=2
  DO j = jstart, jend
    DO i = istart, iend
      z_lgs_rhs_w(i,j,2) = w(i,j,2,nstar)                                          &
        &  - z_Cour_w(i,j,2) * ( w(i,j,3,nstar) - w(i,j,1,nstar) ) * z_beta_m_half

    ENDDO
  ENDDO

  ! The layers from k=3 to k=ke-1
  DO k = 3, ke-1
    DO j = jstart, jend
      DO i = istart, iend

        zC2 = z_Cour_w(i,j,k) * z_Cour_w(i,j,k) 

        ! store advection-operator in help-variable zh1:
        IF ( z_Cour_w(i,j,k) > 0 ) THEN
          zh1 =   ( 1.0_wp/6.0_wp + zC2/12.0_wp ) * w(i,j,k-2,nstar)  &
              & + ( - 1.0_wp      - zC2/4.0_wp  ) * w(i,j,k-1,nstar)  &
              & + ( 1.0_wp/2.0_wp + zC2/4.0_wp  ) * w(i,j,k  ,nstar)  &
              & + ( 1.0_wp/3.0_wp - zC2/12.0_wp ) * w(i,j,k+1,nstar)
        ELSE 
          zh1 = - ( 1.0_wp/6.0_wp + zC2/12.0_wp ) * w(i,j,k+2,nstar)  &
              & - ( - 1.0_wp      - zC2/4.0_wp  ) * w(i,j,k+1,nstar)  &
              & - ( 1.0_wp/2.0_wp + zC2/4.0_wp  ) * w(i,j,k  ,nstar)  &
              & - ( 1.0_wp/3.0_wp - zC2/12.0_wp ) * w(i,j,k-1,nstar)
        END IF
        z_lgs_rhs_w(i,j,k) = w(i,j,k,nstar)                  &
          &   - z_Cour_w(i,j,k) * zh1 * z_beta_m

      ENDDO
    ENDDO
  ENDDO

  ! layer k=ke
  DO j = jstart, jend
    DO i = istart, iend
      z_lgs_rhs_w(i,j,ke) = w(i,j,ke,nstar)                                              &
        &   - z_Cour_w(i,j,ke) * ( w(i,j,ke+1,nstar)-w(i,j,ke-1,nstar) ) * z_beta_m_half
    ENDDO
  ENDDO

  ! Bottom layer k=ke+1
  DO j = jstart, jend
    DO i = istart, iend
      IF ( z_Cour(i,j,k) > 0 ) THEN
        z_lgs_rhs_w(i,j,ke1) = w(i,j,ke1,nstar)              &
          &  - z_Cour_w(i,j,ke1) * z_beta_m * ( w(i,j,ke1,nstar)-w(i,j,ke,nstar) )
      ELSE
        ! use explicit upwind:
        z_lgs_rhs_w(i,j,ke1) = w(i,j,ke1,nstar)              &
          &  - z_Cour_w(i,j,ke1) * 1.0_wp * ( w(i,j,ke1,nstar)-w(i,j,ke,nstar) )
      END IF
    ENDDO
  ENDDO

  IF ( flag_vector_version ) THEN 
    CALL solve_5banddiag_vec( z_lgs_w, z_lgs_rhs_w, zh_w, &
        ie, je, istart, iend, jstart, jend, ke1 ) 
    ! attention: z_lgs_w has now changed
  ELSE 
    CALL solve_5banddiag( z_lgs_w, z_lgs_rhs_w, zh_w, &
        ie, je, istart, iend, jstart, jend, ke1 ) 
  END IF

  DO k = 1, ke1
    DO j = jstart, jend
      DO i = istart, iend
        wadvt(i,j,k) = ( zh_w(i,j,k) - w(i,j,k, nstar) ) * z_dt_recip  
      ENDDO
    ENDDO
  ENDDO

  !----------------------------------------------------------------------------
  ! Section 4: Setup of tridiagonal matrix systems resulting from the implicit
  !            numerical formulation of advection.
  !            - variables at u-positions (u)
  !----------------------------------------------------------------------------

  ! Band-diagonal-matrix

  ! Top layer k=1 (one-sided diff.)  
  DO j = jstartu, jendu
    DO i = istartu, iendu

      z_Cour(i,j,1) = sqrtg_r_u(i,j,1) *                         &
        &             ( wcon(i,  j,1) + wcon(i,  j,2)            &
        &             + wcon(i+1,j,1) + wcon(i+1,j,2) ) * z_dt_quart 

      IF ( z_Cour(i,j,1) < 0.0_wp ) THEN
        z_lgs(i,j,1,1) = 0.0_wp
        z_lgs(i,j,1,2) = 0.0_wp
        z_lgs(i,j,1,3) = 1.0_wp - z_Cour(i,j,1) * z_beta_p
        z_lgs(i,j,1,4) =            + z_Cour(i,j,1) * z_beta_p
        z_lgs(i,j,1,5) = 0.0_wp
      ELSE
        ! use explicit upwind:
        z_lgs(i,j,1,1) = 0.0_wp
        z_lgs(i,j,1,2) = 0.0_wp
        z_lgs(i,j,1,3) = 1.0_wp
        z_lgs(i,j,1,4) = 0.0_wp
        z_lgs(i,j,1,5) = 0.0_wp
      END IF
    ENDDO
  ENDDO

  ! 2. Layer k=2 (centered diff. 2. order)
  DO j = jstartu, jendu
    DO i = istartu, iendu

      z_Cour(i,j,2) =  sqrtg_r_u(i,j,2) *                        &
        &             ( wcon(i,  j,2) + wcon(i,  j,3)            &
        &             + wcon(i+1,j,2) + wcon(i+1,j,3) ) * z_dt_quart 

      z_lgs(i,j,2,1) = 0.0_wp
      z_lgs(i,j,2,2) = - z_Cour(i,j,2) * z_beta_p_half
      z_lgs(i,j,2,3) = 1.0_wp
      z_lgs(i,j,2,4) = - z_lgs(i,j,2,2)
      z_lgs(i,j,2,5) = 0.0_wp

    ENDDO
  ENDDO

  ! The layers from k=3 to k=ke-2
  DO k = 3, ke-2
    DO j = jstartu, jendu
      DO i = istartu, iendu

        z_Cour(i,j,k) =  sqrtg_r_u(i,j,k) *                         &
          &             ( wcon(i,  j,k) + wcon(i,  j,k+1)           &
          &             + wcon(i+1,j,k) + wcon(i+1,j,k+1) ) * z_dt_quart 

        zC2 = z_Cour(i,j,k) * z_Cour(i,j,k)
        zh2 = z_Cour(i,j,k) * z_beta_p

        IF ( z_Cour(i,j,k) > 0 ) THEN
          z_lgs(i,j,k,1) = zh2 * ( 1.0_wp/6.0_wp + zC2/12.0_wp )
          z_lgs(i,j,k,2) = zh2 * ( - 1.0_wp      - zC2/4.0_wp  )
          z_lgs(i,j,k,3) = zh2 * ( 1.0_wp/2.0_wp + zC2/4.0_wp  ) + 1.0_wp
          z_lgs(i,j,k,4) = zh2 * ( 1.0_wp/3.0_wp - zC2/12.0_wp )
          z_lgs(i,j,k,5) = 0.0_wp
        ELSE 
          z_lgs(i,j,k,1) = 0.0_wp
          z_lgs(i,j,k,2) = - zh2 * ( 1.0_wp/3.0_wp - zC2/12.0_wp )
          z_lgs(i,j,k,3) = - zh2 * ( 1.0_wp/2.0_wp + zC2/4.0_wp  ) + 1.0_wp
          z_lgs(i,j,k,4) = - zh2 * ( - 1.0_wp      - zC2/4.0_wp  )
          z_lgs(i,j,k,5) = - zh2 * ( 1.0_wp/6.0_wp + zC2/12.0_wp )
        END IF

      ENDDO
    ENDDO
  ENDDO

  ! Layer k=ke-1 (centered diff. 2. order)
  DO j = jstartu, jendu
    DO i = istartu, iendu

      z_Cour(i,j,ke-1) =  sqrtg_r_u(i,j,ke-1) *                        &
        &                ( wcon(i,  j,ke-1)+ wcon(i,  j,ke)            &
        &                + wcon(i+1,j,ke-1)+ wcon(i+1,j,ke) ) * z_dt_quart  

      z_lgs(i,j,ke-1,1) = 0.0_wp
      z_lgs(i,j,ke-1,2) = - z_Cour(i,j,ke-1) * z_beta_p_half
      z_lgs(i,j,ke-1,3) = 1.0_wp
      z_lgs(i,j,ke-1,4) = - z_lgs(i,j,ke-1,2)
      z_lgs(i,j,ke-1,5) = 0.0_wp

    ENDDO
  ENDDO

  ! The bottom layer k=ke
  DO j = jstartu, jendu
    DO i = istartu, iendu

      z_Cour(i,j,ke) =  sqrtg_r_u(i,j,ke) *                          &
        &              ( wcon(i,  j,ke)+ wcon(i,  j,ke+1)            &
        &              + wcon(i+1,j,ke)+ wcon(i+1,j,ke+1) ) * z_dt_quart  

      IF ( z_Cour(i,j,k) > 0 ) THEN
        z_lgs(i,j,ke,1) = 0.0_wp
        z_lgs(i,j,ke,2) =            - z_Cour(i,j,ke) * z_beta_p
        z_lgs(i,j,ke,3) = 1.0_wp + z_Cour(i,j,ke) * z_beta_p
        z_lgs(i,j,ke,4) = 0.0_wp
        z_lgs(i,j,ke,5) = 0.0_wp
      ELSE
        ! use explicit upwind:
        z_lgs(i,j,ke,1) = 0.0_wp
        z_lgs(i,j,ke,2) = 0.0_wp
        z_lgs(i,j,ke,3) = 1.0_wp
        z_lgs(i,j,ke,4) = 0.0_wp
        z_lgs(i,j,ke,5) = 0.0_wp
      END IF
    ENDDO
  ENDDO

  ! --- right hand side for u ---

  ! Top layer (k=1)
  DO j = jstartu, jendu
    DO i = istartu, iendu
      IF ( z_Cour(i,j,1) < 0.0_wp ) THEN
        z_lgs_rhs(i,j,1) = u(i,j,1,nstar)                                   &
          &  - z_Cour(i,j,1) * z_beta_m * ( u(i,j,2,nstar)-u(i,j,1,nstar) ) 
      ELSE
        ! use explicit upwind:
        z_lgs_rhs(i,j,1) = u(i,j,1,nstar)                                   &
          &  - z_Cour(i,j,1) * 1.0_wp * ( u(i,j,2,nstar)-u(i,j,1,nstar) ) 
      END IF
    ENDDO
  ENDDO

  ! 2. layer k=2
  DO j = jstartu, jendu
    DO i = istartu, iendu
      z_lgs_rhs(i,j,2) = u(i,j,2,nstar)                                          &
        &  - z_Cour(i,j,2) * ( u(i,j,3,nstar) - u(i,j,1,nstar) ) * z_beta_m_half 

    ENDDO
  ENDDO

  ! The layers from k=3 to k=ke-2
  DO k = 3, ke-2
    DO j = jstartu, jendu
      DO i = istartu, iendu

        zC2 = z_Cour(i,j,k) * z_Cour(i,j,k) 

        ! store advection-operator in help-variable zh1:
        IF ( z_Cour(i,j,k) > 0 ) THEN
          zh1 =   ( 1.0_wp/6.0_wp + zC2/12.0_wp ) * u(i,j,k-2,nstar)  &
              & + ( - 1.0_wp      - zC2/4.0_wp  ) * u(i,j,k-1,nstar)  &
              & + ( 1.0_wp/2.0_wp + zC2/4.0_wp  ) * u(i,j,k  ,nstar)  &
              & + ( 1.0_wp/3.0_wp - zC2/12.0_wp ) * u(i,j,k+1,nstar)
        ELSE 
          zh1 = - ( 1.0_wp/6.0_wp + zC2/12.0_wp ) * u(i,j,k+2,nstar)  &
              & - ( - 1.0_wp      - zC2/4.0_wp  ) * u(i,j,k+1,nstar)  &
              & - ( 1.0_wp/2.0_wp + zC2/4.0_wp  ) * u(i,j,k  ,nstar)  &
              & - ( 1.0_wp/3.0_wp - zC2/12.0_wp ) * u(i,j,k-1,nstar)
        END IF
        z_lgs_rhs(i,j,k) = u(i,j,k,nstar)                  &
          &   - z_Cour(i,j,k) * zh1 * z_beta_m 

      ENDDO
    ENDDO
  ENDDO

  ! layer k=ke-1
  DO j = jstartu, jendu
    DO i = istartu, iendu
      z_lgs_rhs(i,j,ke-1) = u(i,j,ke-1,nstar)                                          &
        &   - z_Cour(i,j,ke-1) * ( u(i,j,ke,nstar)-u(i,j,ke-2,nstar) ) * z_beta_m_half
    ENDDO
  ENDDO

  ! Bottom layer k=ke
  DO j = jstartu, jendu
    DO i = istartu, iendu
      IF ( z_Cour(i,j,k) > 0 ) THEN
        z_lgs_rhs(i,j,ke) = u(i,j,ke,nstar)              &
          &  - z_Cour(i,j,ke) * z_beta_m * ( u(i,j,ke,nstar)-u(i,j,ke-1,nstar) )
      ELSE
        ! use explicit upwind:
        z_lgs_rhs(i,j,ke) = u(i,j,ke,nstar)              &
          &  - z_Cour(i,j,ke) * 1.0_wp * ( u(i,j,ke,nstar)-u(i,j,ke-1,nstar) )
      END IF
    ENDDO
  ENDDO

  IF ( flag_vector_version ) THEN 
    CALL solve_5banddiag_vec( z_lgs, z_lgs_rhs, zh,  &
       ie, je, istartu, iendu, jstartu, jendu, ke ) 
    ! attention: z_lgs has now changed
  ELSE
    CALL solve_5banddiag( z_lgs, z_lgs_rhs, zh,  &
       ie, je, istartu, iendu, jstartu, jendu, ke ) 
  END IF

  DO k = 1, ke
    DO j = jstartu, jendu
      DO i = istartu, iendu
        uadvt(i,j,k) = ( zh(i,j,k) - u(i,j,k, nstar) ) * z_dt_recip  
      ENDDO
    ENDDO
  ENDDO

  !----------------------------------------------------------------------------
  ! Section 5: Setup of tridiagonal matrix systems resulting from the implicit
  !            numerical formulation of advection.
  !            - variables at v-positions (v)
  !----------------------------------------------------------------------------

  ! Band-diagonal-matrix

  ! Top layer k=1 (one-sided diff.)  
  DO j = jstartv, jendv
    DO i = istartv, iendv

      z_Cour(i,j,1) =  sqrtg_r_v(i,j,1) *                         &
        &             ( wcon(i,j,  1) + wcon(i,j,  2)             &
        &             + wcon(i,j+1,1) + wcon(i,j+1,2) ) * z_dt_quart

      IF ( z_Cour(i,j,1) < 0.0_wp ) THEN
        z_lgs(i,j,1,1) = 0.0_wp
        z_lgs(i,j,1,2) = 0.0_wp
        z_lgs(i,j,1,3) = 1.0_wp - z_Cour(i,j,1) * z_beta_p
        z_lgs(i,j,1,4) =            + z_Cour(i,j,1) * z_beta_p
        z_lgs(i,j,1,5) = 0.0_wp
      ELSE
        ! use explicit upwind:
        z_lgs(i,j,1,1) = 0.0_wp
        z_lgs(i,j,1,2) = 0.0_wp
        z_lgs(i,j,1,3) = 1.0_wp
        z_lgs(i,j,1,4) = 0.0_wp
        z_lgs(i,j,1,5) = 0.0_wp
      END IF

    ENDDO
  ENDDO

  ! 2. Layer k=2 (centered diff. 2. order)
  DO j = jstartv, jendv
    DO i = istartv, iendv

      z_Cour(i,j,2) =  sqrtg_r_v(i,j,2) *                         &
        &             ( wcon(i,j,  2) + wcon(i,j,  3)             &
        &             + wcon(i,j+1,2) + wcon(i,j+1,3) ) * z_dt_quart

      z_lgs(i,j,2,1) = 0.0_wp
      z_lgs(i,j,2,2) = - z_Cour(i,j,2) * z_beta_p_half
      z_lgs(i,j,2,3) = 1.0_wp
      z_lgs(i,j,2,4) = - z_lgs(i,j,2,2)
      z_lgs(i,j,2,5) = 0.0_wp

    ENDDO
  ENDDO

  ! The layers from k=3 to k=ke-2
  DO k = 3, ke-2
    DO j = jstartv, jendv
      DO i = istartv, iendv

        z_Cour(i,j,k) =  sqrtg_r_v(i,j,k) *                         &
          &             ( wcon(i,j,  k) + wcon(i,j,  k+1)           &
          &             + wcon(i,j+1,k) + wcon(i,j+1,k+1) ) * z_dt_quart

        zC2 = z_Cour(i,j,k) * z_Cour(i,j,k)
        zh2 = z_Cour(i,j,k) * z_beta_p

        IF ( z_Cour(i,j,k) > 0 ) THEN
          z_lgs(i,j,k,1) = zh2 * ( 1.0_wp/6.0_wp + zC2/12.0_wp )
          z_lgs(i,j,k,2) = zh2 * ( - 1.0_wp      - zC2/4.0_wp  )
          z_lgs(i,j,k,3) = zh2 * ( 1.0_wp/2.0_wp + zC2/4.0_wp  ) + 1.0_wp
          z_lgs(i,j,k,4) = zh2 * ( 1.0_wp/3.0_wp - zC2/12.0_wp )
          z_lgs(i,j,k,5) = 0.0_wp
        ELSE 
          z_lgs(i,j,k,1) = 0.0_wp
          z_lgs(i,j,k,2) = - zh2 * ( 1.0_wp/3.0_wp - zC2/12.0_wp )
          z_lgs(i,j,k,3) = - zh2 * ( 1.0_wp/2.0_wp + zC2/4.0_wp  ) + 1.0_wp
          z_lgs(i,j,k,4) = - zh2 * ( - 1.0_wp      - zC2/4.0_wp  )
          z_lgs(i,j,k,5) = - zh2 * ( 1.0_wp/6.0_wp + zC2/12.0_wp ) 
        END IF

      ENDDO
    ENDDO
  ENDDO

  ! Layer k=ke-1 (centered diff. 2. order)
  DO j = jstartv, jendv
    DO i = istartv, iendv
      z_Cour(i,j,ke-1) =  sqrtg_r_v(i,j,ke-1) *                      &
        &                ( wcon(i,j,  ke-1)+ wcon(i,j,  ke)          &
        &                + wcon(i,j+1,ke-1)+ wcon(i,j+1,ke) ) * z_dt_quart
 
      z_lgs(i,j,ke-1,1) = 0.0_wp
      z_lgs(i,j,ke-1,2) = - z_Cour(i,j,ke-1) * z_beta_p_half
      z_lgs(i,j,ke-1,3) = 1.0_wp
      z_lgs(i,j,ke-1,4) = - z_lgs(i,j,ke-1,2)
      z_lgs(i,j,ke-1,5) = 0.0_wp

    ENDDO
  ENDDO

  ! The bottom layer k=ke
  DO j = jstartv, jendv
    DO i = istartv, iendv

      z_Cour(i,j,ke) =  sqrtg_r_v(i,j,ke) *                        &
        &              ( wcon(i,j,  ke)+ wcon(i,j,  ke+1)          &
        &              + wcon(i,j+1,ke)+ wcon(i,j+1,ke+1) ) * z_dt_quart
 
      IF ( z_Cour(i,j,k) > 0 ) THEN
        z_lgs(i,j,ke,1) = 0.0_wp
        z_lgs(i,j,ke,2) =            - z_Cour(i,j,ke) * z_beta_p
        z_lgs(i,j,ke,3) = 1.0_wp + z_Cour(i,j,ke) * z_beta_p
        z_lgs(i,j,ke,4) = 0.0_wp
        z_lgs(i,j,ke,5) = 0.0_wp
      ELSE
        ! use explicit upwind:
        z_lgs(i,j,ke,1) = 0.0_wp
        z_lgs(i,j,ke,2) = 0.0_wp
        z_lgs(i,j,ke,3) = 1.0_wp
        z_lgs(i,j,ke,4) = 0.0_wp
        z_lgs(i,j,ke,5) = 0.0_wp
      END IF
    ENDDO
  ENDDO

  ! --- right hand side for v ---

  ! Top layer (k=1)      
  DO j = jstartv, jendv
    DO i = istartv, iendv
      IF ( z_Cour(i,j,1) < 0.0_wp ) THEN
        z_lgs_rhs(i,j,1) = v(i,j,1,nstar)                                   &
          &  - z_Cour(i,j,1) * z_beta_m * ( v(i,j,2,nstar)-v(i,j,1,nstar) ) 
      ELSE
        ! use explicit upwind:
        z_lgs_rhs(i,j,1) = v(i,j,1,nstar)                                   &
          &  - z_Cour(i,j,1) * 1.0_wp * ( v(i,j,2,nstar)-v(i,j,1,nstar) ) 
      END IF
    ENDDO
  ENDDO

  ! 2. layer k=2
  DO j = jstartv, jendv
    DO i = istartv, iendv
      z_lgs_rhs(i,j,2) = v(i,j,2,nstar)                                          &
        &  - z_Cour(i,j,2) * ( v(i,j,3,nstar) - v(i,j,1,nstar) ) * z_beta_m_half 

    ENDDO
  ENDDO

  ! The layers from k=3 to k=ke-2
  DO k = 3, ke-2
    DO j = jstartv, jendv
      DO i = istartv, iendv

        zC2 = z_Cour(i,j,k) * z_Cour(i,j,k) 

        ! store advection-operator in help-variable zh1:
        IF ( z_Cour(i,j,k) > 0 ) THEN
          zh1 =   ( 1.0_wp/6.0_wp + zC2/12.0_wp ) * v(i,j,k-2,nstar)  &
              & + ( - 1.0_wp      - zC2/4.0_wp  ) * v(i,j,k-1,nstar)  &
              & + ( 1.0_wp/2.0_wp + zC2/4.0_wp  ) * v(i,j,k  ,nstar)  &
              & + ( 1.0_wp/3.0_wp - zC2/12.0_wp ) * v(i,j,k+1,nstar)
        ELSE 
          zh1 = - ( 1.0_wp/6.0_wp + zC2/12.0_wp ) * v(i,j,k+2,nstar)  &
              & - ( - 1.0_wp      - zC2/4.0_wp  ) * v(i,j,k+1,nstar)  &
              & - ( 1.0_wp/2.0_wp + zC2/4.0_wp  ) * v(i,j,k  ,nstar)  &
              & - ( 1.0_wp/3.0_wp - zC2/12.0_wp ) * v(i,j,k-1,nstar)
        END IF
        z_lgs_rhs(i,j,k) = v(i,j,k,nstar)                  &
          &   - z_Cour(i,j,k) * zh1 * z_beta_m 

      ENDDO
    ENDDO
  ENDDO

  ! layer k=ke-1
  DO j = jstartv, jendv
    DO i = istartv, iendv
      z_lgs_rhs(i,j,ke-1) = v(i,j,ke-1,nstar)                                          &
        &   - z_Cour(i,j,ke-1) * ( v(i,j,ke,nstar)-v(i,j,ke-2,nstar) ) * z_beta_m_half 
    ENDDO
  ENDDO

  ! Bottom layer k=ke
  DO j = jstartv, jendv
    DO i = istartv, iendv
      IF ( z_Cour(i,j,k) > 0 ) THEN
        z_lgs_rhs(i,j,ke) = v(i,j,ke,nstar)              &
          &  - z_Cour(i,j,ke) * z_beta_m * ( v(i,j,ke,nstar)-v(i,j,ke-1,nstar) ) 
      ELSE
        ! use explicit upwind:
        z_lgs_rhs(i,j,ke) = v(i,j,ke,nstar)              &
          &  - z_Cour(i,j,ke) * 1.0_wp * ( v(i,j,ke,nstar)-v(i,j,ke-1,nstar) ) 
      END IF
    ENDDO
  ENDDO

  IF ( flag_vector_version ) THEN 
    CALL solve_5banddiag_vec( z_lgs, z_lgs_rhs, zh,  &
        ie, je, istartv, iendv, jstartv, jendv, ke ) 
    ! attention: z_lgs has now changed
  ELSE
    CALL solve_5banddiag( z_lgs, z_lgs_rhs, zh,  &
        ie, je, istartv, iendv, jstartv, jendv, ke ) 
  END IF

  DO k = 1, ke
    DO j = jstartv, jendv
      DO i = istartv, iendv
        vadvt(i,j,k) = ( zh(i,j,k) - v(i,j,k, nstar) ) * z_dt_recip  
      ENDDO
    ENDDO
  ENDDO

  DEALLOCATE( z_Cour, z_lgs, z_lgs_rhs, zh )

  DEALLOCATE( z_Cour_w, z_lgs_w, z_lgs_rhs_w, zh_w )

  !----------------------------------------------------------------------------
  ! End of subroutine complete_tend_uvwtpp_CN3Crow
  !----------------------------------------------------------------------------

END SUBROUTINE complete_tend_uvwtpp_CN3Crow


END MODULE src_slow_tendencies_rk
