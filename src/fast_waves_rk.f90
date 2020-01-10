!+ Module containing procedures to compute fast wave terms for RK-scheme
!------------------------------------------------------------------------------

MODULE fast_waves_rk

!------------------------------------------------------------------------------
!
! Description:
!   This module contains procedure(s) (for the use in src_runge_kutta) to
!   calculate the new values of the prognostic variables
!   u, v, w, pp and T at time level n+1 (nnew).
!
! Current Code Owner: DWD, Michael Baldauf    
!  phone:  +49  69  8062 2733
!  fax:    +49  69  8062 3721
!  email:  michael.baldauf@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 3.7        2004/02/18 Jochen Foerstner
!  Initial release
! 3.8        2004/03/23 Jochen Foerstner
!  Some optimizations and corrections for l2dim and lperi; editorial changes
! 3.13       2004/12/03 Jochen Foerstner
!  Use of new variables for reciprocal sqrt(G); Test of radiative lateral
!  boundary condition
! 3.14       2005/01/25 Jochen Foerstner
!  Upwind-1st order and centered-2nd order operator in advection form
! 3.16       2005/07/22 Jochen Foerstner
!  Converted former external subroutine to a module with a module procedure
! 3.18       2006/03/03 Jochen Foerstner
!  Updates
! 3.21       2006/12/04 Jochen Foerstner
!  More Updates
! V3_23        2007/03/30 Lucio Torrisi
!  Deallocation of memory at the end of time stepping (for DFI)
! V4_4         2008/07/16 Ulrich Schaettler
!  Adapted interface of get_timings
!  Consideration of the sign of dt for DFI (Lucio Torrisi)
! V4_5         2008/09/10 Guenther Zaengl
! Implementation of new sponge layer near upper model boundary
! Finish implementation of radiative lateral boundary condition
! More efficient w lower boundary condition
! V4_8_1_1     2009/04/21 Ulrich Schaettler
!  Activate old w lower boundary condition again for Bott scheme
! V4_9         2009/07/16 Guenther Zaengl, Christian Bollmann, Ulrich Schaettler
!  More accurate discretization of metric terms
!  Option for using potential temperature as advected variable
!  Avoid vectorization over small loops in treatment at boundary (Bollmann)
!  Implemented ON_ADB directives (Christian Bollmann)
! V4_11        2009/11/30 Guenther Zaengl
!  Implemented switch for choosing bottom boundary for w (itype_bbc_w)
!  Lower boundary condition 2nd order for metric terms
!  Partly vectorized radiative boundary condition
! V4_12        2010/05/11 Guenther Zaengl
!  Use 2nd order bottom boundary condition only for itype_bbc_w=0/2/4
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_15        2010/11/19 Michael Baldauf
!  Summarize calculation of weighting coefficients for vertical interpolation in a
!  new subroutine 'weighting_factors_full2half'  (G. Zaengl)
! V4_17        2011/02/24 Ulrich Blahak
!  Adapted interface of exchg_boundaries;
!  corrected kzdims(1:20) -> kzdims(1:24);
!  eliminated my_peri_neigh and all if-branches for single processor 2D-runs 
!    and/or periodic BCs (now part of exchg_boundaries() in all cases);
!  eliminated lfreeslip_surf (free-slip BC and/or no-surface-heat/moisture-flux
!    conditions may now be imposed by new switches lnosurffluxes_m/h in namelist IDEAL);
!  introduced lnosurffluxes_m instead of lfreeslip_surf when limiting the 
!    extrapolated values of U, V to the surface; 
!  added new dependency on src_artifdata because of lnosurffluxes_m.
! V4_23        2012/05/10 Michael Baldauf
!  Shifted SR weighting_factors_full2half and fields wgtfac* to new module
!    grid_metrics_utilities
!  Cleaned up the code and removed unnecessary (internal) options
! V4_24        2012/06/22 Michael Baldauf
!  Boundary exchange of suten, svten now takes place in org_runge_kutta
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer, Hans-Juergen Panitz
!  Replaced qx-variables by using them from the tracer module
!  Introduced nexch_tag for MPI boundary exchange tag to replace ntstep (HJP)
! V4_27        2013/03/19 Michael Baldauf, Astrid Kerkweg, Ulrich Schaettler
!  Removed unused variables iup, jup;
!  Introduced error message, if wrong value for itype_bbc_w is used
!  MESSy interface introduced (AK)
! V4_28        2013/07/12 Ulrich Schaettler
!  Pass kflat as parameter to avoid dependency on vgrid_refatm_utils
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
! V5_2         2015-05-21 Ulrich Schaettler
!  Reorganized numbering of MPI datatypes ID
! V5_3         2015-10-09 Oliver Fuhrer
!  Implemented serialization comments
! V5_4b        2016-07-12 Oliver Fuhrer
!  Integrate new lateral boundary condition module for FW-solvers
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

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
     ke_rd,        & ! lowest level with Rayleigh-damping

! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from
!    the other ones because of the use of the staggered Arakawa-B-grid.
!
!    zonal direction
     istart,       & ! start index for the forecast of w, t, qv, qc and pp
     iend,         & ! end index for the forecast of w, t, qv, qc and pp
     istartu,      & ! start index for the forecast of u
     iendu,        & ! end index for the forecast of u
     istartv,      & ! start index for the forecast of v
     iendv,        & ! end index for the forecast of v
     istartpar,    & ! start index for computations in the parallel program
     iendpar,      & ! end index for computations in the parallel program

!    meridional direction
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
     edadlat,      & ! 1 / (radius of the earth * dlat)

! 5. variables for the time discretization and related variables
! --------------------------------------------------------------

     dt,           & ! long time-step
     betasw,       & ! beta-variable for treatment of soundwaves    in w
     betagw,       & ! beta-variable for treatment of gravity-waves in w
     beta2sw,      & ! beta-variable for treatment of soundwaves    in p*, T*
     beta2gw,      & ! beta-variable for treatment of gravity-waves in p*, T*

! 8. Organizational variables to handle the COSMO humidity tracers
! ----------------------------------------------------------------
    idt_qv, idt_qc

! end of data_modelconfig

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                      (unit)
! -----------------------------------------------
     rho0       ,    & ! reference density at the full model levels    (kg/m3)
     p0         ,    & ! reference pressure at main levels             ( Pa)
     dt0dz      ,    & ! temperature gradient of reference atmosphere  ( K/m )
     t0         ,    & ! reference temperature                         ( K   )
     hhl        ,    & ! geometical height of half model levels        ( m )

! 2. external parameter fields                                         (unit)
! ----------------------------
     crlat      ,    & ! cosine of transformed latitude
     acrlat     ,    & ! 1 / ( crlat * radius of the earth )           ( 1/m )

! 3. prognostic variables                                              (unit)
! -----------------------
     u          ,    & ! zonal wind speed                              ( m/s )
     v          ,    & ! meridional wind speed                         ( m/s )
     w          ,    & ! vertical wind speed (defined on half levels)  ( m/s )
     t          ,    & ! temperature                                   (  k  )
     pp         ,    & ! deviation from the reference pressure         ( pa  )

! 6. fields that are computed in the parametrization and dynamics      (unit )
! ---------------------------------------------------------------
     rdcoef     ,    & ! Rayleigh damping coefficients
     qrs        ,    & ! precipitation water (water loading)           (kg/kg)
     rho               ! density of moist air

! end of data_fields

!------------------------------------------------------------------------------

USE data_parallel,      ONLY :  &
     num_compute,     & ! number of compute PEs
     nboundlines,     & ! number of boundary lines of the domain for which
                        ! no forecast is computed = overlapping boundary
                        ! lines of the subdomains
     ldatatypes,      & ! if .TRUE.: use MPI-Datatypes for some communications
     ltime_barrier,   & ! if .TRUE.: use additional barriers for determining
                        ! the load-imbalance
     ncomm_type,      & ! type of communication
     my_cart_id,      & ! rank of this subdomain in the cartesian communicator
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
     nstart,       & ! first time step of the forecast
     nstop,        & ! last time step of the forecast
     ntstep,       & ! actual time step
                     ! indices for permutation of three time levels
     nnow,         & ! corresponds to ntstep
     nnew,         & ! corresponds to ntstep + 1

! 7. additional control variables
! -------------------------------
     llm,          & ! if .TRUE., running with a lowered upper boundary
     lperi_x,      & ! if lartif_data=.TRUE.: periodic boundary conditions in x-dir.
                     !               .FALSE.:   or with Davies conditions
     lperi_y,      & ! if lartif_data=.TRUE.: periodic boundary conditions in y-dir.
                     !               .FALSE.:   or with Davies conditions
     lradlbc,      & ! if lartif_data=.TRUE.: radiative lateral boundary conditions
                     !               .FALSE.:   or with Davies conditions
     l2dim,        & ! if lartif_data=.TRUE.: 2dimensional model version or
                     !               .FALSE.:   full 3dimensional version
     ltime,        & ! detailled timings of the program are given
     irk_order,    & ! order of the Runge-Kutta scheme in dyn. core
     iadv_order,   & ! order of the advection scheme in dyn. core
     itheta_adv,   & ! =0: use T' (perturbation temperature) for advection
                     ! =1: use theta' (perturbation potential temperature)
                     ! =2: use theta (full potential temperature)
     itype_bbc_w,  & ! Bottom boundary condition for vertical wind
                     ! =0/1: RK-like method following iadv_order
                     ! =2/3: differencing following iadv_order without RK stepping
                     ! =4/5: Fourth-order centered differences
                     ! 0/2/4: include quadratic extrapolation of horizontal wind to sfc
                     ! 1/3/5: no extrapolation of horizontal wind to sfc
     lw_freeslip,  & ! if .TRUE.: with free slip lateral boundary condition and
                     ! if .FALSE. specified lateral boundary values for w
     ldyn_bbc,     & ! if .TRUE., dynamical bottom boundary condition
     lspubc ,      & ! if .TRUE., use Rayleigh damping in the upper levels
     ltur ,        & ! if .TRUE., turbulence parameterization is used
     itype_spubc     ! type of Rayleigh damping in the upper levels

!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 2. physical constants and related variables
! -------------------------------------------
     r_d,          & ! gas constant for dry air
     rvd_m_o,      & ! r_v/r_d - 1
     cpdr,         & ! 1 / cp_d
     gamma,        & ! 1 / (1 - rdocp)   ( = cp_d/cv_d)
     g               ! acceleration due to gravity

!------------------------------------------------------------------------------

USE src_tracer,       ONLY: trcr_get, trcr_errorstr

!------------------------------------------------------------------------------

USE environment     , ONLY :  exchg_boundaries, comm_barrier, model_abort
USE time_utilities  , ONLY :  get_timings, i_barrier_waiting_dyn, i_fast_waves,&
                              i_fast_waves_comm, i_fast_waves_barrier

USE grid_metrics_utilities, ONLY: &
     sqrtg_r_s  ,    & ! reciprocal square root of G at skalar points  ( 1/m )
     sqrtg_r_u  ,    & ! reciprocal square root of G at u points       ( 1/m )
     sqrtg_r_v  ,    & ! reciprocal square root of G at v points       ( 1/m )
     sqrtg_r_w  ,    & ! reciprocal square root of G at w points       ( 1/m )
     wgtfac,         & ! weighting factor for vertical interpolation
     wgtfac_u,       & ! weighting factor for vertical interpolation
     wgtfac_v,       & ! weighting factor for vertical interpolation
     wgtfacq,        & ! weighting factor for vertical interpolation
     wgtfacq_u,      & ! weighting factor for vertical interpolation
     wgtfacq_v         ! weighting factor for vertical interpolation

USE meteo_utilities , ONLY :  calrho_tp_pp

USE numeric_utilities_rk, ONLY:  &
     udsdx,        &  ! udsdx_* interface
     udsdx_up5_xy     ! optimized version for default UP-5th hori. adv.

USE src_artifdata,    ONLY : lnosurffluxes_m

USE src_lbc, ONLY : lbc_value, lbc_tendency, lbc_zerograd, lbc_copy

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Local scalars:
! -----------------------------
LOGICAL, PRIVATE ::                   &
  llin_tdyn_rk           ! flag to switch between linearization with nnow
                         ! and intermediate TD-values of Runge-Kutta scheme

INTEGER (KIND=iintegers), PRIVATE ::  &
  im, ip,              & ! index-boundaries for advection stencil
  imipmx,              & ! maximum of the stencil indices
  ilowu, jlowv           !  start- and end-indices for subdomains with
                         !  or without neighbours

REAL    (KIND=wp   ),     PRIVATE ::  &
  zfy,                 & !
  zbsp, zbsm,          & !
  zbgp, zbgm,          & !
  zb2sp, zb2sm,        & !
  zb2gp, zb2gm,        & !
  zb2smob2sp,          & !
  zb2gmob2gp

! Local (allocatable) arrays:
! -----------------------------
REAL(KIND=wp),     DIMENSION(:),     ALLOCATABLE, PRIVATE ::  &
  zfx, zfydn, zfyd, zfyds                              !

REAL(KIND=wp),     DIMENSION(:,:,:), ALLOCATABLE, PRIVATE ::  &
  za_opt, zb_opt, zc_opt, zd, zcw1, zcw2,            & !
  zcp1_h, zct1_h, zcp1, zcp2, zct1, zct2,            & !
  zalphab, zalphat, zbetab, zbetat,                  & !
  zbalphab, zbalphat, zwtens,                        & !
  zrhoqx_i, zrhoqy_i, zdzdx, zdzdy                     !

! Variables for dyn. bottom BC
!-----------------------------
REAL (KIND=wp),     DIMENSION(:,:), ALLOCATABLE, PRIVATE ::  &
  xdzdx, xdzdy, xdzdx_v, xdzdy_u,  &
  xrhsx, xrhsy, xrhsz,             &
  xlhsx, xlhsy, zcwp


! For radiative lateral boundary condition
! ----------------------------------------
LOGICAL, PRIVATE ::                   &
  lradlbc_u, lradlbc_v, lradlbc_w,    & ! flags for radiative lateral BC
  lradlbc_pp, lradlbc_t

INTEGER (KIND=iintegers), PRIVATE ::  &
  islb, ielb, jslb, jelb,             & ! indices of boundary lines to set
  islbu, ielbu, jslbv, jelbv            ! dito for u / v

! Weighting factors
REAL (KIND=wp),     ALLOCATABLE, PRIVATE :: wrdfac(:)

! Factors for theta advection
REAL (KIND=wp)     :: cvovcp, cvovcp2

!==============================================================================

CONTAINS

!==============================================================================

SUBROUTINE fast_waves_runge_kutta (nsmsteps, dts, xkd, suten, svten, swten,   &
                         stpten, sppten, t_bd_tend, pp_bd_tend, u_bd_tend,    &
                         v_bd_tend, w_bd_tend, irk, a_rk, b_rk, kzflat )

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure calculates the new values of the prognostic variables
!   u, v, w, pp and T at time level n+1 (nnew).
!
! Method:
!   The Runge-Kutta/split explicit time integration scheme is used: Within
!   each Runge-Kutta step the slow tendencies of the prognostic variables
!   due to advection and diffusion (which have been calculated before and
!   are stored in the corresponding tendency fieleds) are held constant
!   wheras the fast terms, which describe sound and gravity wave propagation,
!   are stepped forward by using a small timestep dts.
!   Horizontal wave propagation is treated explicitly, for vertical wave
!   propagation an implicit scheme is used.
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------
INTEGER (KIND=iintegers), INTENT (IN) ::  &
  nsmsteps,           & ! number of small time steps for integration
  irk,                & ! counter for Runge-Kutta steps
  kzflat                ! level-index where half-levels bcome flat

REAL    (KIND=wp   ),     INTENT (IN) ::  &
  dts,                & ! small time step
  xkd                   ! constant factor for divergence-type damping

REAL    (KIND=wp   ),     INTENT (INOUT) ::  &
  suten (ie,je,ke),   & ! slow tendency of u-velocity (*dts)
  svten (ie,je,ke),   & ! slow tendency of v-velocity (*dts)
  swten (ie,je,ke1),  & ! slow tendency of w-velocity (*dts)
  stpten(ie,je,ke),   & ! slow tendency of temperature (*dts)
  sppten(ie,je,ke)      ! slow tendency of pressure perturbation (*dts)

REAL    (KIND=wp   ), INTENT (IN) ::  &
  t_bd_tend(ie,je,ke), &
  pp_bd_tend(ie,je,ke),&
  u_bd_tend(ie,je,ke), &
  v_bd_tend(ie,je,ke), &
  w_bd_tend(ie,je,ke1)

REAL    (KIND=wp   ),     INTENT(IN) ::  &
  a_rk, b_rk

! Local scalars:
! -----------------------------
INTEGER (KIND=iintegers) ::  &
  kzdims(24),          & ! Vertical dimensions for exchg_boundaries
  n_rk,                & !
  i,  j,  k,           & !  Loop indices
  izts,                & !  small timestep loop index
  kp1, km1,            & !  k+1, k-1
  izstata,             & ! error status at allocation
  izstatd                ! error status at deallocation

REAL    (KIND=wp   )     ::  &
  zrofac, zptot,       & !
  zbqt, zbqb,          & !
  zdpdx, zdpdy,        & !
  zdzpz, zfact,        & !
  zdpdzu, zdpdzv,      & !
  zpgradx, zpgrady,    & !
  zaw, zcw, zwdiv,     & !
  zwdz_s, zwavg_s,     & !
  zrhoqx, zrhoqy,      & !
  zdenom, zdts_q         !

! Local (allocatable) arrays:
! -----------------------------

REAL(KIND=wp),     ALLOCATABLE ::  zhml(:,:,:)

! Local (automatic) arrays:
! -----------------------------
REAL    (KIND=wp   )     ::  &
  zpi    (ie,je,ke)  ,      & !
  za     (ie,je,2:ke),      & !
  zb     (ie,je,2:ke),      & !
  zc     (ie,je,2:ke),      & !
  ztdiv  (ie,je,ke)  ,      & !
  zppten (ie,ke)     ,      & !
  ztpten (ie,ke)     ,      & !
  zttotr (ie,je,ke)  ,      & !
  zt0otp0(ie,je,ke)  ,      & !
  zdpdz  (ie,je)     ,      & !
  zphl   (ie,je,ke1) ,      & !
  zuhl   (ie,je,ke1) ,      & !
  zvhl   (ie,je,ke1)          !

 ! Variables for dyn. bottom BC
 !-----------------------------
REAL (KIND=wp)      :: &
  zxrhsx, zxrhsy, zxrhsz, zxdzdx, zxdzdy

! For error handling
! ------------------
INTEGER (KIND=iintegers) ::  izerror
CHARACTER (LEN=255)      ::  yzerrmsg
CHARACTER (LEN=25)       ::  yzroutine


! For radiative lateral boundary condition
! ----------------------------------------
REAL (KIND=wp),     ALLOCATABLE :: &
  zu_lbx(:,:,:),  zu_lby(:,:,:),   & ! arrays for radiative lateral BC
  zv_lbx(:,:,:),  zv_lby(:,:,:),   & ! to store values at the boundaries
  zw_lbx(:,:,:),  zw_lby(:,:,:),   & !
  zpp_lbx(:,:,:), zpp_lby(:,:,:),  & !
  zt_lbx(:,:,:),  zt_lby(:,:,:)      !

! For microphysics tracers
!-------------------------
REAL (KIND=wp),     POINTER :: &
  qv_new  (:,:,:)=> NULL(),   & ! QV at tlev=nnew
  qv_now  (:,:,:)=> NULL(),   & ! QV at tlev=nnow
  qc_new  (:,:,:)=> NULL(),   & ! QC at tlev=nnew
  qc_now  (:,:,:)=> NULL()      ! QC at tlev=nnow

!
! End of header
!==============================================================================

yzroutine = 'fast_waves_runge_kutta'

!------------------------------------------------------------------------------
! Begin Subroutine fast_waves_runge_kutta
!------------------------------------------------------------------------------

  IF ( ntstep == nstart .AND. irk == 1 ) THEN

    IF (.NOT. ALLOCATED(wrdfac))  ALLOCATE(wrdfac(ke))
    wrdfac(:) = 1.0_wp

    ! Damping factor for new sponge condition (Klemp et al., 2008, MWR)
    IF (lspubc .AND. (itype_spubc == 3) ) THEN
      DO k = 1, ke_rd
        wrdfac(k) = 1.0_wp/(1.0_wp + 10.0_wp*rdcoef(k)*dts)
      ENDDO
    ENDIF


    llin_tdyn_rk = .FALSE.  ! linearization with intermediate RK-values

    ! Factors needed for theta advection
    cvovcp = 1.0_wp - r_d*cpdr
    cvovcp2 = 0.5_wp*cvovcp*(cvovcp-1.0_wp)

    IF (.NOT.( lperi_x .OR. lperi_y) ) THEN

      ! some preparations for radiative lateral boundary condition
      ! ----------------------------------------------------------
      IF (lradlbc) THEN
        lradlbc_u  = .TRUE.
        lradlbc_v  = .TRUE.
        lradlbc_w  = .TRUE.
        lradlbc_pp = .TRUE.
        lradlbc_t  = .TRUE.
      ELSE
        lradlbc_u  = .FALSE.
        lradlbc_v  = .FALSE.
        lradlbc_w  = .FALSE.
        lradlbc_pp = .FALSE.
        lradlbc_t  = .FALSE.
      ENDIF

      islb  = istart
      ielb  = iend
      jslb  = jstart
      jelb  = jend
      islbu = istartu-1
      ielbu = iendu+1
      jslbv = jstartv-1
      jelbv = jendv+1

    END IF

    ALLOCATE(               &
      zfx   (je),           &
      zfydn (je),           &
      zfyd  (je),           &
      zfyds (je),           &
      zcp2  (ie,je,ke),     &
      zct2  (ie,je,ke),     &
      zhml  (ie,je,ke),     &
      zdzdx (ie,je,ke),     &
      zdzdy (ie,je,ke),     &
      STAT=izstata )

    ! some preparations to simplify the choice of
    ! different advection shemes 
    ! setting of the indices of the advected stencil
    SELECT CASE(iadv_order)
    CASE(1,2)
      im = -1_iintegers
      ip = 1_iintegers
    CASE(3,4)
      im = -2_iintegers
      ip = 2_iintegers
    CASE(5,6)
      im = -3_iintegers
      ip = 3_iintegers
    END SELECT
    imipmx = MAX( 2_iintegers, ABS(im), ABS(ip) )

!------------------------------------------------------------------------------
!  Section 1: Setup of parameters for time integration
!------------------------------------------------------------------------------

    !$ser savepoint FastWavesRKUnittest.Init-in LargeTimeStep=0 RKStageNumber=0 SmallTimeStep=0
    !$ser data dts=dts

    ! betasw  is the Ikawa beta parameter for soundwaves    in w
    ! betagw            "                 for gravity-waves in w
    ! beta2sw           "                 for soundwaves    in p*, T*
    ! beta2gw           "                 for gravity-waves in p*, T*
    ! (=0: centered, =1 backward, =-1 forward)
    zbsp  = ( 1.0_wp + betasw  ) * 0.5_wp
    zbsm  = ( 1.0_wp - betasw  ) * 0.5_wp
    zbgp  = ( 1.0_wp + betagw  ) * 0.5_wp
    zbgm  = ( 1.0_wp - betagw  ) * 0.5_wp
    zb2sp = ( 1.0_wp + beta2sw ) * 0.5_wp
    zb2sm = ( 1.0_wp - beta2sw ) * 0.5_wp
    zb2gp = ( 1.0_wp + beta2gw ) * 0.5_wp
    zb2gm = ( 1.0_wp - beta2gw ) * 0.5_wp

    IF (my_cart_neigh(1) == -1) THEN
      ilowu = istartu
    ELSE
      ilowu = istartu-1
    ENDIF

    IF (my_cart_neigh(4) == -1) THEN
      jlowv = jstartv
    ELSE
      jlowv = jstartv-1
    ENDIF

    ! Factors for the setup of the tridiagonal gaussian matrix and some
    ! local arrays for horizontal and vertical discretization

    zfy     = edadlat

    DO j = jstart - 1, jend + 1
      zfx(j)   = eddlon*acrlat(j,1)
      zfydn(j) = eddlat*acrlat(j,1)*crlat(j  ,2)
      zfyds(j) = eddlat*acrlat(j,1)*crlat(j-1,2)
      zfyd(j)  = eddlat*acrlat(j,1)*crlat(j  ,1)
    ENDDO

    ! Setup of some constant arrays
    IF ( itheta_adv /= 2 ) THEN
      ! gravity-waves are treated (partly) implicit in p*- and T*-equation
      DO  k = 1, ke
        DO  j = jstart, jend
          DO  i = istart, iend
            zcp2(i,j,k) =   zb2gp * rho0(i,j,k)  * 0.5_wp*g
            zct2(i,j,k) = - zb2gp * dt0dz(i,j,k) * 0.5_wp
          ENDDO
        ENDDO
      ENDDO
    ELSE IF ( itheta_adv == 2 ) THEN
      DO  k = 1, ke
        DO  j = jstart, jend
          DO  i = istart, iend
            zcp2(i,j,k) =   zb2gp * rho0(i,j,k)  * 0.5_wp*g
            zct2(i,j,k) = 0.0_wp
          ENDDO
        ENDDO
      ENDDO
    END IF

    ! Compute height of main levels
    DO  k = 1, ke
      DO  j = jstart-1, jend+1
!CDIR ON_ADB(hhl)
        DO  i = istart-1, iend+1
          zhml(i,j,k) = 0.5_wp*( hhl(i,j,k) + hhl(i,j,k+1) )
        ENDDO
      ENDDO
    ENDDO

    DO  k = 1, ke
      zfact = 0.5_wp
!?    IF ( k==1 .OR. k==ke ) zfact = 0.5_wp
      DO  j = jstart-1, jend
!CDIR ON_ADB(zhml)
        DO  i = istart-1, iend
          zdzdx(i,j,k) = zfact * sqrtg_r_u(i,j,k) *           &
                        ( zhml(i+1,j,k) - zhml(i,j,k) )
          zdzdy(i,j,k) = zfact * sqrtg_r_v(i,j,k) *           &
                        ( zhml(i,j+1,k) - zhml(i,j,k) )
        ENDDO
      ENDDO
    ENDDO

    IF ( ldyn_bbc ) THEN
      
      ALLOCATE(           &
        xrhsx(ie,je),     &
        xrhsy(ie,je),     &
        xrhsz(ie,je),     &
        xlhsx(ie,je),     &
        xlhsy(ie,je),     &
        xdzdx(ie,je),     &
        xdzdy(ie,je),     &
        xdzdx_v(ie,je),   &
        xdzdy_u(ie,je),   &
        STAT=izstata )
#ifdef _SERIALZE
      xrhsx = 0.0_wp
      xrhsy = 0.0_wp
      xrhsz = 0.0_wp
      xlhsx = 0.0_wp
      xlhsy = 0.0_wp
      xdzdx = 0.0_wp
      xdzdy = 0.0_wp
      xdzdx_v = 0.0_wp
      xdzdy_u = 0.0_wp
#endif
    
      DO  j = jlowv, jendv+1
!CDIR ON_ADB(zhml)
        DO  i = istartv-1, iendv
          xdzdx(i,j) = ( zhml(i+1,j,ke) - zhml(i,j,ke) ) * zfx(j)
        ENDDO
      ENDDO
      DO  j = jstartu-1, jendu
!CDIR ON_ADB(zhml)
!CDIR ON_ADB(xdzdy)
        DO  i = ilowu, iendu+1
          xdzdy(i,j) = ( zhml(i,j+1,ke) - zhml(i,j,ke) ) * zfy
        ENDDO
      ENDDO
      
      DO j = jstartu, jendu
!CDIR ON_ADB(xdzdy)
        DO i = ilowu, iendu
          xdzdy_u(i,j) = 0.25_wp * ( xdzdy(i+1,j)   + xdzdy(i,j)      &
                                       + xdzdy(i+1,j-1) + xdzdy(i,j-1) )
          xlhsx(i,j) = 1.0_wp / ( xdzdx(i,j)**2 + xdzdy_u(i,j)**2     &
                                                    + 1.0_wp )
        ENDDO
      ENDDO

      DO j = jlowv, jendv
!CDIR ON_ADB(xdzdx)
        DO i = istartv, iendv
          xdzdx_v(i,j) = 0.25_wp * ( xdzdx(i,j+1)   + xdzdx(i,j)      &
                                       + xdzdx(i-1,j+1) + xdzdx(i-1,j) )
          xlhsy(i,j) = 1.0_wp / ( xdzdx_v(i,j)**2 + xdzdy(i,j)**2     &
                                                      + 1.0_wp )
        ENDDO
      ENDDO
        
    END IF

    !$ser savepoint FastWavesRKUnittest.Init-out LargeTimeStep=0 RKStageNumber=0 SmallTimeStep=0
    !$ser data cp2=zcp2 ct2=zct2 wgtfacq=wgtfacq wgtfac=wgtfac
    !$ser data zdzdx=zdzdx zdzdy=zdzdy
    !$ser data xdzdx=xdzdx xdzdy=xdzdy xlhsx=xlhsx xlhsy=xlhsy if (ldyn_bbc)

    DEALLOCATE( zhml, STAT=izstatd )

  END IF

  IF (.NOT.( lperi_x .OR. lperi_y)) THEN

    ! some preparations for radiative lateral boundary condition
    ! ----------------------------------------------------------
    IF ( lradlbc_u )  THEN
      ALLOCATE( zu_lbx(6,je,ke), zu_lby(ie,6,ke), STAT=izstata )
      zu_lbx(:,:,:) = 0.0_wp
      zu_lby(:,:,:) = 0.0_wp
    END IF
    IF ( lradlbc_v )  THEN
      ALLOCATE( zv_lbx(6,je,ke), zv_lby(ie,6,ke), STAT=izstata )
      zv_lbx(:,:,:) = 0.0_wp
      zv_lby(:,:,:) = 0.0_wp
    END IF
    IF ( lradlbc_w )  THEN
      ALLOCATE( zw_lbx(6,je,ke1),zw_lby(ie,6,ke1),STAT=izstata )
      zw_lbx(:,:,:) = 0.0_wp
      zw_lby(:,:,:) = 0.0_wp
    END IF
    IF ( lradlbc_pp ) THEN
      ALLOCATE( zpp_lbx(6,je,ke),zpp_lby(ie,6,ke),STAT=izstata )
      zpp_lbx(:,:,:) = 0.0_wp
      zpp_lby(:,:,:) = 0.0_wp
    END IF
    IF ( lradlbc_t )  THEN
      ALLOCATE( zt_lbx(6,je,ke), zt_lby(ie,6,ke), STAT=izstata )
      zt_lbx(:,:,:) = 0.0_wp
      zt_lby(:,:,:) = 0.0_wp
    END IF

  END IF


  izerror  = 0
  yzerrmsg = '   '
  kzdims(:)= 0_iintegers

  ! LBC for divergence is zero, although these values should
  ! not be used in the computations for the fast-waves solver
  CALL lbc_value( zpi, 0.0_wp,                                      &
                  ie, je, ke, istart, iend, jstart, jend, 1, ke)

!------------------------------------------------------------------------------
!  Section 2: Some preparations for time integration
!             and for optimization
!------------------------------------------------------------------------------

  IF ( irk == 1 .OR. .NOT.llin_tdyn_rk ) THEN
    n_rk = nnow
  ELSE
    n_rk = nnew
  END IF

  ALLOCATE( zd(ie,je,ke1), STAT=izstata )

  ! Get needed microphysics tracers out of the tracer structure
  CALL trcr_get (izerror, idt_qv, ptr_tlev = nnew, ptr = qv_new)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort (my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get (izerror, idt_qv, ptr_tlev = nnow, ptr = qv_now)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort (my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get (izerror, idt_qc, ptr_tlev = nnew, ptr = qc_new)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort (my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get (izerror, idt_qc, ptr_tlev = nnow, ptr = qc_now)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort (my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  IF ( irk == 1 .OR. llin_tdyn_rk ) THEN

    ALLOCATE(                 &
      zrhoqx_i (ie,je,ke),    &
      zrhoqy_i (ie,je,ke),    &
      zcw1     (ie,je,ke),    &
      zcw2     (ie,je,ke),    &
      zcp1_h   (ie,je,ke),    &
      zct1_h   (ie,je,ke),    &
      zcp1     (ie,je,ke),    &
      zct1     (ie,je,ke),    &
      zalphab  (ie,je,ke),    &
      zalphat  (ie,je,ke),    &
      zbalphab (ie,je,ke),    &
      zbalphat (ie,je,ke),    &
      zbetab   (ie,je,ke),    &
      zbetat   (ie,je,ke),    &
      zwtens   (ie,je,ke1),   &
      za_opt   (ie,je,2:ke),  &
      zb_opt   (ie,je,2:ke),  &
      zc_opt   (ie,je,2:ke),  &
      STAT=izstata )

    IF ( ldyn_bbc ) THEN
       ALLOCATE( zcwp(ie,je), STAT=izstata )
       !$ser zero zcwp
    ENDIF

    !$ser savepoint FastWavesRKUnittest.PrepareStep-in LargeTimeStep=ntstep RKStageNumber=irk SmallTimeStep=0
    !$ser data t_nnow=t(:,:,:,nnow) pp_nnow=pp(:,:,:,nnow) qrs=qrs(:,:,:) rho=rho
    !$ser tracer QV@nnow QC@nnow

    ! Avoid division by zrhoqx and zrhoqy in the time integration
    DO  k = 1, ke
      DO  j = jstart-1, jend+1
!CDIR ON_ADB(rho)
        DO  i = istart-1, iend
          zrhoqx          = 0.5_wp*( rho(i,j,k) + rho(i+1,j,k) )
          zrhoqx_i(i,j,k) = zfx(j) / zrhoqx
        ENDDO
      ENDDO
      DO  j = jstart-1, jend
!CDIR ON_ADB(rho)
        DO  i = istart-1, iend+1
          zrhoqy          = 0.5_wp*( rho(i,j,k) + rho(i,j+1,k) )
          zrhoqy_i(i,j,k) = zfy / zrhoqy
        ENDDO
      ENDDO
    ENDDO

    ! Factors for the setup of the tridiagonal gaussian matrix and some
    ! local arrays for horizontal and vertical discretization

    DO  k = 2, ke
      DO j = jstart, jend
!CDIR ON_ADB(rho)
!CDIR ON_ADB(rho0)
!CDIR ON_ADB(qc_now)
!CDIR ON_ADB(qrs)
!CDIR ON_ADB(qv_now)
        DO  i = istart, iend
          zrofac      = 1._wp / ( rho(i,j,k) + rho(i,j,k-1) )
          zcw1(i,j,k) = 2.0_wp * sqrtg_r_w(i,j,k) * zrofac
          zrofac      = ( rho0(i,j,k) + rho0(i,j,k-1) ) * zrofac
          zcw2(i,j,k) = 0.5_wp * g * zrofac
          ! Vertical wind tendency due to bouyancy and water loading
          zbqb  =  rvd_m_o*qv_now(i,j,k  ) - qc_now(i,j,k  ) - qrs(i,j,k  )
          zbqt  =  rvd_m_o*qv_now(i,j,k-1) - qc_now(i,j,k-1) - qrs(i,j,k-1)
          zwtens(i,j,k) = zcw2(i,j,k)*( zbqt + zbqb )
        ENDDO
      ENDDO
    ENDDO

    ! zwtens(:,:,ke1) is used to store the w-tendency (buoyancy eff.)
    !                 at the full level ke needed to the dynamical
    !                 bottom boundary condition
    IF ( ldyn_bbc ) THEN
      DO j = jstart-1, jend+1
!CDIR ON_ADB(rho)
!CDIR ON_ADB(rho0)
!CDIR ON_ADB(qc_now)
!CDIR ON_ADB(qrs)
!CDIR ON_ADB(qv_now)
        DO i = istart-1, iend+1
          zrofac = rho0(i,j,ke) / rho (i,j,ke) * g
          zwtens(i,j,ke1) = zrofac * ( rvd_m_o*qv_now(i,j,ke)           &
                                             - qc_now(i,j,ke)           &
                                             - qrs(i,j,ke) )
        ENDDO
      ENDDO
    ENDIF

    IF (itheta_adv /= 2) THEN
      DO  k = 1, ke
        DO  j = jstart-1, jend+1
          DO  i = istart-1, iend+1
            zttotr(i,j,k)  = 1.0_wp / ( t0(i,j,k) + t(i,j,k,n_rk) )
          ENDDO
        ENDDO
      ENDDO
    ELSE
      DO  k = 1, ke
        DO  j = jstart-1, jend+1
          DO  i = istart-1, iend+1
            zttotr(i,j,k)  = 1.0_wp / t(i,j,k,n_rk)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    DO  k = 1, ke
      DO  j = jstart, jend
        DO  i = istart, iend
          zptot  = p0(i,j,k) + pp(i,j,k,n_rk)
          zcp1_h(i,j,k)  = gamma * zptot
          zct1_h(i,j,k)  = zcp1_h(i,j,k) * cpdr/rho(i,j,k)
        ENDDO
      ENDDO
    ENDDO

    IF (itheta_adv == 0) THEN
      DO  k = 1, ke
        DO  j = jstart, jend
          DO  i = istart, iend
            zt0otp0(i,j,k) = t0(i,j,k) * zttotr(i,j,k) / p0(i,j,k)
          ENDDO
        ENDDO
      ENDDO
    ELSE
      DO  k = 1, ke
        DO  j = jstart, jend
          DO  i = istart, iend
            zt0otp0(i,j,k) = (cvovcp + cvovcp2 * pp(i,j,k,n_rk)/p0(i,j,k) ) &
                              * t0(i,j,k) * zttotr(i,j,k) / p0(i,j,k)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF ( ldyn_bbc .AND. (itheta_adv == 0)) THEN
      DO j = jstart-1, jend+1
        DO i = istart-1, iend+1
          zcwp(i,j) = zttotr(i,j,ke) / ( r_d*rho0(i,j,ke) )
        ENDDO
      ENDDO
    ELSE IF ( ldyn_bbc .AND. (itheta_adv /= 0)) THEN
      DO j = jstart-1, jend+1
        DO i = istart-1, iend+1
          zcwp(i,j) = zttotr(i,j,ke)*t0(i,j,ke)/p0(i,j,ke)
        ENDDO
      ENDDO
    ENDIF

    ! soundwaves are treated (partly) implicit in p*- and T*-equation
    DO  k = 1, ke
      DO  j = jstart, jend
        DO  i = istart, iend
          zcp1(i,j,k) = zb2sp * zcp1_h(i,j,k) * sqrtg_r_s(i,j,k)
          zct1(i,j,k) = zb2sp * zct1_h(i,j,k) * sqrtg_r_s(i,j,k)
        ENDDO
      ENDDO
    ENDDO


    ! Setup of the lower (za), the main (zb) and the upper (zc) diagonals
    ! of the Gauss matrix
    ! (only in the first Runge-Kutta step, if llin_tdyn_rk=.FALSE.)

    ! soundwaves and gravity-waves are treated (partly) implicit
    ! in p*- and T*-equation

    IF (itheta_adv == 0) THEN
      ! gravity-waves are treated (partly) implicit in w-equation
      DO  k = 2, ke
        DO  j = jstart, jend
!CDIR ON_ADB(zt0otp0)
!CDIR ON_ADB(zcp1)
!CDIR ON_ADB(zcp2)
!CDIR ON_ADB(zct1)
!CDIR ON_ADB(zct2)
!CDIR ON_ADB(zttotr)
          DO  i = istart, iend

            zaw             = zcw2(i,j,k) * zt0otp0(i,j,k-1)
            zalphat(i,j,k)  =      zcw1(i,j,k) +      zaw
            zbalphat(i,j,k) = zbsp*zcw1(i,j,k) + zbgp*zaw

            zcw             = zcw2(i,j,k) * zt0otp0(i,j,k  )
            zalphab(i,j,k)  =      zcw1(i,j,k) -      zcw
            zbalphab(i,j,k) = zbsp*zcw1(i,j,k) - zbgp*zcw

            za_opt(i,j,k) =   &
              - zbalphat(i,j,k) * ( zcp1(i,j,k-1) - zcp2(i,j,k-1) )
            zb_opt(i,j,k) =   &
                zbalphab(i,j,k) * ( zcp1(i,j,k  ) - zcp2(i,j,k  ) )   &
              + zbalphat(i,j,k) * ( zcp1(i,j,k-1) + zcp2(i,j,k-1) )
            zc_opt(i,j,k) =   &
              - zbalphab(i,j,k) * ( zcp1(i,j,k  ) + zcp2(i,j,k  ) )

            zbetab(i,j,k) = zcw2(i,j,k) * zttotr(i,j,k  )
            zbetat(i,j,k) = zcw2(i,j,k) * zttotr(i,j,k-1)

            zaw       = zbetat(i,j,k)*zbgp
            zcw       = zbetab(i,j,k)*zbgp
            za_opt(i,j,k) = za_opt(i,j,k)                         &
              + zaw * ( zct1(i,j,k-1) - zct2(i,j,k-1) )
            zb_opt(i,j,k) = zb_opt(i,j,k)                         &
              + zcw * ( zct1(i,j,k  ) - zct2(i,j,k  ) )           &
                - zaw * ( zct1(i,j,k-1) + zct2(i,j,k-1) )
            zc_opt(i,j,k) = zc_opt(i,j,k)                         &
              - zcw * ( zct1(i,j,k  ) + zct2(i,j,k  ) )

          ENDDO
        ENDDO
      ENDDO
    ELSE IF (itheta_adv == 1) THEN
      ! version for perturbation potential temperature
      DO  k = 2, ke
        DO  j = jstart, jend
!CDIR ON_ADB(zt0otp0)
!CDIR ON_ADB(zcp1)
!CDIR ON_ADB(zcp2)
!CDIR ON_ADB(zct1)
!CDIR ON_ADB(zct2)
!CDIR ON_ADB(zttotr)
          DO  i = istart, iend

            zaw             = zcw2(i,j,k) * zt0otp0(i,j,k-1)
            zalphat(i,j,k)  =      zcw1(i,j,k) +      zaw
            zbalphat(i,j,k) = zbsp*zcw1(i,j,k) + zbgp*zaw

            zcw             = zcw2(i,j,k) * zt0otp0(i,j,k  )
            zalphab(i,j,k)  =      zcw1(i,j,k) -      zcw
            zbalphab(i,j,k) = zbsp*zcw1(i,j,k) - zbgp*zcw

            za_opt(i,j,k) =   &
              - zbalphat(i,j,k) * ( zcp1(i,j,k-1) - zcp2(i,j,k-1) )
            zb_opt(i,j,k) =   &
                zbalphab(i,j,k) * ( zcp1(i,j,k  ) - zcp2(i,j,k  ) )   &
              + zbalphat(i,j,k) * ( zcp1(i,j,k-1) + zcp2(i,j,k-1) )
            zc_opt(i,j,k) =   &
              - zbalphab(i,j,k) * ( zcp1(i,j,k  ) + zcp2(i,j,k  ) )

            zbetab(i,j,k) = zcw2(i,j,k) * zttotr(i,j,k  )
            zbetat(i,j,k) = zcw2(i,j,k) * zttotr(i,j,k-1)

            zaw       = zbetat(i,j,k)*zbgp
            zcw       = zbetab(i,j,k)*zbgp
            za_opt(i,j,k) = za_opt(i,j,k) - zaw * zct2(i,j,k-1)
            zb_opt(i,j,k) = zb_opt(i,j,k)                         &
              - zcw * zct2(i,j,k) - zaw * zct2(i,j,k-1)
            zc_opt(i,j,k) = zc_opt(i,j,k) - zcw * zct2(i,j,k  )

          ENDDO
        ENDDO
      ENDDO
    ELSE IF (itheta_adv == 2) THEN
      ! version for full potential temperature
      DO  k = 2, ke
        DO  j = jstart, jend
!CDIR ON_ADB(zt0otp0)
!CDIR ON_ADB(zcp1)
!CDIR ON_ADB(zcp2)
!CDIR ON_ADB(zct1)
!CDIR ON_ADB(zct2)
!CDIR ON_ADB(zttotr)
          DO  i = istart, iend

            zaw             = zcw2(i,j,k) * zt0otp0(i,j,k-1)
            zalphat(i,j,k)  =      zcw1(i,j,k) +      zaw
            zbalphat(i,j,k) = zbsp*zcw1(i,j,k) + zbgp*zaw

            zcw             = zcw2(i,j,k) * zt0otp0(i,j,k  )
            zalphab(i,j,k)  =      zcw1(i,j,k) -      zcw
            zbalphab(i,j,k) = zbsp*zcw1(i,j,k) - zbgp*zcw

            za_opt(i,j,k) =   &
              - zbalphat(i,j,k) * ( zcp1(i,j,k-1) - zcp2(i,j,k-1) )
            zb_opt(i,j,k) =   &
                zbalphab(i,j,k) * ( zcp1(i,j,k  ) - zcp2(i,j,k  ) )   &
              + zbalphat(i,j,k) * ( zcp1(i,j,k-1) + zcp2(i,j,k-1) )
            zc_opt(i,j,k) =   &
              - zbalphab(i,j,k) * ( zcp1(i,j,k  ) + zcp2(i,j,k  ) )

            zbetab(i,j,k) = zcw2(i,j,k) * zttotr(i,j,k  )
            zbetat(i,j,k) = zcw2(i,j,k) * zttotr(i,j,k-1)

          ENDDO
        ENDDO
      ENDDO
    END IF

    !$ser savepoint FastWavesRKUnittest.PrepareStep-out LargeTimeStep=ntstep RKStageNumber=irk SmallTimeStep=0
    !$ser data cw1=zcw1 cw2=zcw2 bwlwten=zwtens ttotr=zttotr t0otp0=zt0otp0 cp1_h=zcp1_h ct1_h=zct1_h  &
    !$ser&     cp1=zcp1 ct1=zct1 cwp=zcwp

  END IF      !!! IF ( irk == 1 .OR. llin_tdyn_rk ) THEN ...

  zdts_q = dts*dts

  ! Setup of the lower (za), main (zb) and the upper (zc) diagonals
  ! of the Gauss matrix
  DO  k = 2, ke
    DO  j = jstart, jend
      DO  i = istart, iend
        za(i,j,k) = zdts_q * za_opt(i,j,k)
        zc(i,j,k) = zdts_q * zc_opt(i,j,k)
        zb(i,j,k) = 1.0_wp + zdts_q * zb_opt(i,j,k)
      ENDDO
    ENDDO
  ENDDO

  !
  ! LU-Decomposition of the matrix
  !
  DO  j = jstart, jend
!CDIR ON_ADB(za)
!CDIR ON_ADB(zb)
    DO  i = istart, iend
      zdenom     =  zb(i,j,ke)
      zb(i,j,ke) = 1.0_wp / zdenom
      za(i,j,ke) = -za(i,j,ke) * zb(i,j,ke)
    ENDDO
  ENDDO
  DO  k = ke-1, 2, -1
    DO  j = jstart, jend
!CDIR ON_ADB(za)
      DO  i = istart, iend
        zdenom    =  zc(i,j,k) * za(i,j,k+1) + zb(i,j,k)
        zb(i,j,k) = 1.0_wp / zdenom
        za(i,j,k) = -za(i,j,k) * zb(i,j,k)
      ENDDO
    ENDDO
  ENDDO

  !$ser savepoint FastWavesRKUnittest.PrepareStage-in LargeTimeStep=ntstep RKStageNumber=irk SmallTimeStep=0
  !$ser data dts=dts swten=swten bwlwten=zwtens

  !-----------------------------------------------------------------
  ! add vertical wind tendency due to bouyancy and water loading
  !-----------------------------------------------------------------
  DO  k = 2, ke
    DO j = jstart, jend
      DO  i = istart, iend
        swten(i,j,k) = swten(i,j,k) + zwtens(i,j,k) * dts
      ENDDO
    ENDDO
  ENDDO

  IF ( ldyn_bbc ) THEN
    DO j = jstart-1, jend+1
      DO i = istart-1, iend+1
        swten(i,j,ke1) = swten(i,j,ke1) + zwtens(i,j,ke1) * dts
      ENDDO
    ENDDO
  ENDIF


  IF ( irk == irk_order .OR. llin_tdyn_rk ) THEN
    DEALLOCATE(                &
      zcw1, zcw2, zwtens,      &
      za_opt, zb_opt, zc_opt,  &
      STAT=izstatd )
  END IF


  !-----------------------------------------------------------------
  ! Calculate the start values of the variables for the integration
  !-----------------------------------------------------------------

  ! NOTE: the initialization below overwrites the BC for all of the fields
  !       at timelevel nnew and thus undoes the BCs which have been setup
  !       in initialize_loop() -> we need to setup correct BCs during the
  !       fast-waves solver again (see below)!

  DO k = 1, ke
    DO j = 1, je
      DO i = 1, ie
        u (i,j,k,nnew) = a_rk * u(i,j,k,nnow)  + b_rk * u(i,j,k,nnew)
        v (i,j,k,nnew) = a_rk * v(i,j,k,nnow)  + b_rk * v(i,j,k,nnew)
        w (i,j,k,nnew) = a_rk * w(i,j,k,nnow)  + b_rk * w(i,j,k,nnew)
        pp(i,j,k,nnew) = a_rk * pp(i,j,k,nnow) + b_rk * pp(i,j,k,nnew)
        t (i,j,k,nnew) = a_rk * t(i,j,k,nnow)  + b_rk * t(i,j,k,nnew)
      END DO
    END DO
  END DO
  w(:,:,ke1,nnew) = a_rk * w(:,:,ke1,nnow) + b_rk * w(:,:,ke1,nnew)

  ! Initialize the working pressure perturbation including
  ! one boundary line, which is used for the pressure gradient
  ! (unlike the u and v fields, the BC for zpi carrying pp is
  ! not updated on every small timestep)

  ! NOTE: the boundary line is initialized irrespective of whether
  !       we are actually sitting at the global border with this PE
  !       in order to avoid an additional exchange. This assumes that
  !       pp(nnow) and pp(nnew) are updated in at least one line of the
  !       halo-region when entering fast_waves_rk()

  ! NOTE: is this really correct, should we keep the border at
  !       a_rk*now + b_rk*new or keep it constant over all RK stages
  !       at nnow?

  DO  k = 1, ke
    zpi(istart-1:iend+1,jstart-1:jend+1,k) =                     &
      pp(istart-1:iend+1,jstart-1:jend+1,k,nnew)
  ENDDO
  CALL lbc_copy( zpi, pp(:,:,:,nnew),                            &
                 ie, je, ke, istart, iend, jstart, jend, 1, ke, nlines=1)

  !$ser savepoint FastWavesRKUnittest.PrepareStage-out LargeTimeStep=ntstep RKStageNumber=irk SmallTimeStep=0
  !$ser data swten=swten wbbctens_stage=swten(:,:,ke1)

!------------------------------------------------------------------------------
!     ***************************************************************
!     *   begin of time integration using small time steps          *
!     ***************************************************************
!------------------------------------------------------------------------------

  loop_over_small_timesteps:  DO  izts = 1, nsmsteps

  !$ser savepoint TimeIntegratorUnittest.DoSmallStep-in LargeTimeStep=ntstep RKStageNumber=irk SmallTimeStep=izts
  !$ser data u=u(:,:,:,nnew) v=v(:,:,:,nnew) w=w(:,:,:,nnew) t=t(:,:,:,nnew) pp=pp(:,:,:,nnew) zpi=zpi

    ! Set boundary values for t, pp, u, v and w from the boundary tendencies

    IF ( .NOT. lradlbc_t ) THEN
      CALL lbc_tendency( t(:,:,:,nnew), t_bd_tend, dts,                     &
                         ie, je, ke, istart, iend, jstart, jend, 1, ke )
    END IF
    IF ( .NOT. lradlbc_pp ) THEN
      CALL lbc_tendency( pp(:,:,:,nnew), pp_bd_tend, dts,                   &
                         ie, je, ke, istart, iend, jstart, jend, 1, ke )
    END IF
    IF ( .NOT. lradlbc_u ) THEN
      CALL lbc_tendency( u(:,:,:,nnew), u_bd_tend, dts,                     &
                         ie, je, ke, istartu, iendu, jstartu, jendu, 1, ke )
    END IF
    IF ( .NOT. lradlbc_v ) THEN
      CALL lbc_tendency( v(:,:,:,nnew), v_bd_tend, dts,                     &
                         ie, je, ke, istartv, iendv, jstartv, jendv, 1, ke )
    END IF
    IF ( .NOT. lradlbc_w ) THEN
      IF ( lw_freeslip ) THEN
        CALL lbc_zerograd( w(:,:,:,nnew),                                   &
                           ie, je, ke1, istart, iend, jstart, jend, 1, ke1)
      ELSE
        CALL lbc_tendency( w(:,:,:,nnew), w_bd_tend, dts,                   &
                           ie, je, ke1, istart, iend, jstart, jend, 1, ke1)
      ENDIF
    END IF

    IF ( .NOT.( lperi_x .OR. lperi_y ) ) THEN

      ! radiative lateral boundary condition (init)
      ! -----------------------------------------------
      IF ( lradlbc_u ) &
        CALL init_radiative_lbc(u(:,:,:,nnew),zu_lbx(:,:,:),zu_lby(:,:,:),  &
                                islbu,ielbu,jslb,jelb,ke)
      IF ( lradlbc_v ) &
        CALL init_radiative_lbc(v(:,:,:,nnew),zv_lbx(:,:,:),zv_lby(:,:,:),  &
                                islb,ielb,jslbv,jelbv,ke)
      IF ( lradlbc_w ) &
        CALL init_radiative_lbc(w(:,:,:,nnew),zw_lbx(:,:,:),zw_lby(:,:,:),  &
                                islb,ielb,jslb,jelb,ke1)
      IF ( lradlbc_pp ) &
        CALL init_radiative_lbc(pp(:,:,:,nnew),zpp_lbx(:,:,:),zpp_lby(:,:,:), &
                                islb,ielb,jslb,jelb,ke)
      IF ( lradlbc_t ) &
        CALL init_radiative_lbc(t(:,:,:,nnew),zt_lbx(:,:,:),zt_lby(:,:,:),  &
                                islb,ielb,jslb,jelb,ke)

    ENDIF

    
!------------------------------------------------------------------------------
!  Section 3: Integration of the wind components u and v
!------------------------------------------------------------------------------

    ! Periodic boundary conditions on pressure (zpi) within the time loop !

    IF (lperi_x .OR. lperi_y .OR. l2dim) THEN
      IF (ltime) THEN
        CALL get_timings (i_fast_waves, ntstep, dt, izerror)
        IF (ltime_barrier) THEN
          CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
          CALL get_timings (i_fast_waves_barrier, ntstep, dt, izerror)
        ENDIF
      ENDIF
      
      kzdims(1:24)=(/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                 &
           ( 0,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
            ie, je, kzdims, jstartpar, jendpar, 2, nboundlines,             &
            my_cart_neigh, lperi_x, lperi_y, l2dim,                         &
            5000+nexch_tag, .FALSE., ncomm_type, izerror, yzerrmsg,         &
            zpi(:,:,:) )
      IF (ltime) CALL get_timings (i_fast_waves_comm, ntstep, dt, izerror)

      ! NOTE: BCs on zpi have already been set above

    END IF

    !$ser savepoint FastWavesRKUnittest.UV-in LargeTimeStep=ntstep RKStageNumber=irk SmallTimeStep=izts
    !$ser data dts=dts rho=rho zpi=zpi u=u(:,:,:,nnew) v=v(:,:,:,nnew)
    !$ser data suten=suten svten=svten wbbctens_stage=swten(:,:,ke1) wgtfac=wgtfac
    !$ser data cwp=zcwp xdzdx=xdzdx xdzdy=xdzdy xlhsx=xlhsx xlhsy=xlhsy

    ! Precompute p at half levels for efficiency improvement
    ! This is done by linear interpolation in height
    DO  k = MAX(2,kzflat), ke
!CDIR ON_ADB(zpi)
      zphl(istart-1:iend+1,jstart-1:jend+1,k)     =                    &
        zpi(istart-1:iend+1,jstart-1:jend+1,k)    *                    &
        wgtfac(istart-1:iend+1,jstart-1:jend+1,k) +                    &
        zpi(istart-1:iend+1,jstart-1:jend+1,k-1)  *                    &
        (1.0_wp-wgtfac(istart-1:iend+1,jstart-1:jend+1,k))
    ENDDO
    ! Quadratic extrapolation for bottom level
    IF (.NOT. ldyn_bbc) THEN
!CDIR ON_ADB(zpi)
      zphl(istart-1:iend+1,jstart-1:jend+1,ke+1)    =                  &
        zpi(istart-1:iend+1,jstart-1:jend+1,ke)     *                  &
        wgtfacq(istart-1:iend+1,jstart-1:jend+1,1)  +                  &
        zpi(istart-1:iend+1,jstart-1:jend+1,ke-1)   *                  &
        wgtfacq(istart-1:iend+1,jstart-1:jend+1,2)  +                  &
        zpi(istart-1:iend+1,jstart-1:jend+1,ke-2)   *                  &
        wgtfacq(istart-1:iend+1,jstart-1:jend+1,3)
    ENDIF
    ! Linear extrapolation for top level (typically not used)
    IF (kzflat == 1) THEN
      zphl(istart-1:iend+1,jstart-1:jend+1,1)     =                    &
        zpi(istart-1:iend+1,jstart-1:jend+1,1)    *                    &
        (1.0_wp+wgtfac(istart-1:iend+1,jstart-1:jend+1,1)) -       &
        zpi(istart-1:iend+1,jstart-1:jend+1,2)    *                    &
        wgtfac(istart-1:iend+1,jstart-1:jend+1,1)
    ENDIF

    loop_over_levels_3:  DO  k = 1, ke

      kp1 = MIN ( ke, k+1 )
      km1 = MAX ( 1 , k-1 )

      IF ( k < kzflat ) THEN
        
        DO  j = jstartu, jendu
!CDIR ON_ADB(zpi)
          DO  i = ilowu, iendu   
            zpgradx =  (zpi(i+1,j,k) - zpi(i,j,k)) * zrhoqx_i(i,j,k) * dts
            u(i,j,k,nnew) = u(i,j,k,nnew) - zpgradx + suten(i,j,k)
          ENDDO
        ENDDO
        DO  j = jlowv, jendv
!CDIR ON_ADB(zpi)
          DO  i = istartv, iendv   
            zpgrady =  (zpi(i,j+1,k) - zpi(i,j,k)) * zrhoqy_i(i,j,k) * dts
            v(i,j,k,nnew) = v(i,j,k,nnew) - zpgrady + svten(i,j,k)
          ENDDO
        ENDDO
        
      ELSEIF ( (.NOT. ldyn_bbc .AND.              k >= kzflat) .OR.    &
               (      ldyn_bbc .AND. k < ke .AND. k >= kzflat) ) THEN

!CDIR ON_ADB(zphl)
!CDIR ON_ADB(zdpdz)
        zdpdz(istart-1:iend+1,jstart-1:jend+1) =           &
          zphl(istart-1:iend+1,jstart-1:jend+1,k+1) -      &
          zphl(istart-1:iend+1,jstart-1:jend+1,k)

        DO  j = jstartu, jendu
!CDIR ON_ADB(zdpdz)
!CDIR ON_ADB(zpi)
          DO  i = ilowu, iendu   
            zdpdzu  =  zdpdz(i+1,j) + zdpdz(i,j)
            zdzpz   =  zdpdzu * zdzdx(i,j,k)
            zdpdx   =  zpi(i+1,j,k) - zpi(i,j,k)
            zpgradx =  zdpdx + zdzpz
            u(i,j,k,nnew) = u(i,j,k,nnew) - zpgradx*zrhoqx_i(i,j,k) * dts  &
                                          + suten(i,j,k)
          ENDDO
        ENDDO
        DO  j = jlowv, jendv
!CDIR ON_ADB(zdpdz)
!CDIR ON_ADB(zpi)
          DO  i = istartv, iendv   
            zdpdzv  =  zdpdz(i,j+1) + zdpdz(i,j)
            zdzpz   =  zdpdzv * zdzdy(i,j,k)
            zdpdy   =  zpi(i,j+1,k) - zpi(i,j,k)
            zpgrady =  zdpdy + zdzpz
            v(i,j,k,nnew) = v(i,j,k,nnew) - zpgrady*zrhoqy_i(i,j,k) * dts  &
                                          + svten(i,j,k)
          ENDDO
        ENDDO

      ELSE
        
        ! Dynamical Bottom boundary condition (derived from Almut's formulation)
        !
        !  swten(:,:,ke1) is used to store the w-tendency (hor. adv. and buoyancy eff.)
        !                 at the full level ke (computed in src_runge_kutta.f90)
        !

        DO j = jlowv, jendv+1
!CDIR ON_ADB(zpi)
!CDIR ON_ADB(xrhsx)
          DO i = istartv-1, iendv
            xrhsx(i,j) = - zrhoqx_i(i,j,ke)*( zpi(i+1,j,ke)-zpi(i,j,ke) )  &
                         + suten(i,j,ke) / dts
          ENDDO
        ENDDO

        DO j = jstartu-1, jendu
!CDIR ON_ADB(zpi)
!CDIR ON_ADB(xrhsy)
          DO i = ilowu, iendu+1
            xrhsy(i,j) = - zrhoqy_i(i,j,ke)*( zpi(i,j+1,ke)-zpi(i,j,ke) )  &
                         + svten(i,j,ke) / dts
          ENDDO
        ENDDO

        IF (itheta_adv == 0) THEN
          DO j = jstart-1, jend+1
!CDIR ON_ADB(zpi)
!CDIR ON_ADB(xrhsz)
            DO i = istart-1, iend+1
              xrhsz(i,j) = rho0(i,j,ke)/rho(i,j,ke)*g * ( 1.0_wp         &
                         - ( p0(i,j,ke)+zpi(i,j,ke) )*zcwp(i,j) )            &
                         + swten(i,j,ke1) / dts
            ENDDO
          ENDDO
        ELSE
          DO j = jstart-1, jend+1
!CDIR ON_ADB(zpi)
!CDIR ON_ADB(xrhsz)
            DO i = istart-1, iend+1
              xrhsz(i,j) = rho0(i,j,ke)/rho(i,j,ke)*g * ( 1.0_wp         &
                         - ( p0(i,j,ke)+cvovcp*zpi(i,j,ke) )  &
                         * zcwp(i,j) )+ swten(i,j,ke1) / dts
            ENDDO
          ENDDO
        ENDIF

        DO j = jstartu, jendu
!CDIR ON_ADB(xrhsx)
!CDIR ON_ADB(xrhsy)
!CDIR ON_ADB(xrhsz)
          DO i = ilowu, iendu

            zxdzdx = xdzdx(i,j)
            zxdzdy = xdzdy_u(i,j)
            zxrhsy = 0.25_wp * ( xrhsy(i+1,j)   + xrhsy(i,j)      &
                                   + xrhsy(i+1,j-1) + xrhsy(i,j-1) )
            zxrhsz = 0.5_wp * ( xrhsz(i+1,j) + xrhsz(i,j) )

            u(i,j,ke,nnew) = u(i,j,ke,nnew) + dts * ( xlhsx(i,j) * (  &
                                  - zxdzdx**2     * xrhsx(i,j)        &
                                  - zxdzdx*zxdzdy * zxrhsy            &
                                  + zxdzdx        * zxrhsz            &
                                                ) + xrhsx(i,j) )

          ENDDO
        ENDDO

        DO j = jlowv, jendv
!CDIR ON_ADB(xrhsx)
!CDIR ON_ADB(xrhsy)
!CDIR ON_ADB(xrhsz)
          DO i = istartv, iendv

            zxdzdx = xdzdx_v(i,j)
            zxdzdy = xdzdy(i,j)
            zxrhsx = 0.25_wp * ( xrhsx(i,j)   + xrhsx(i-1,j)      &
                                   + xrhsx(i,j+1) + xrhsx(i-1,j+1) )
            zxrhsz = 0.5_wp * ( xrhsz(i,j+1) + xrhsz(i,j) )

            v(i,j,ke,nnew) = v(i,j,ke,nnew) + dts * ( xlhsy(i,j) * (  &
                                  - zxdzdx*zxdzdy * zxrhsx            &
                                  - zxdzdy**2     * xrhsy(i,j)        &
                                  + zxdzdy        * zxrhsz            &
                                                ) + xrhsy(i,j) )

          ENDDO
        ENDDO

      ENDIF

    ENDDO loop_over_levels_3

    !$ser savepoint FastWavesRKUnittest.UV-out LargeTimeStep=ntstep RKStageNumber=irk SmallTimeStep=izts
    !$ser data u=u(:,:,:,nnew) v=v(:,:,:,nnew)
    
    ! Exchange u and v for periodic boundary conditions
    IF ( lperi_x .OR. lperi_y .OR. l2dim) THEN
      IF (ltime) THEN
        CALL get_timings (i_fast_waves, ntstep, dt, izerror)
        IF (ltime_barrier) THEN
          CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
          CALL get_timings (i_fast_waves_barrier, ntstep, dt, izerror)
        ENDIF
      ENDIF
      
      kzdims(1:24)=(/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                   &
           (0      , sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,&
            ie, je, kzdims, jstartpar, jendpar, 2, nboundlines,               &
            my_cart_neigh, lperi_x, lperi_y, l2dim,                           &
            6000+nexch_tag, .FALSE.   , ncomm_type, izerror, yzerrmsg,        &
            u(:,:,:,nnew), v(:,:,:,nnew))
      IF (ltime) CALL get_timings (i_fast_waves_comm, ntstep, dt, izerror)
      
      ! NOTE: only haloexchg and no BC necessary for periodic boundary conditions since the
      !       BCs on u and v have been set inside the small timestep loop

    ENDIF


!------------------------------------------------------------------------------
!  Section 4: Integration of vertical velocity, temperature and
!             pertubation pressure using a vertical implicit scheme
!------------------------------------------------------------------------------

    !$ser savepoint FastWavesRKUnittest.DivergenceWPPTP-in LargeTimeStep=ntstep RKStageNumber=irk SmallTimeStep=izts
    !$ser data u=u(:,:,:,nnew) v=v(:,:,:,nnew) w=w(:,:,:,nnew) pp=pp(:,:,:,nnew) t=t(:,:,:,nnew)
    !$ser data swten=swten sppten=sppten stpten=stpten
    !$ser data wgtfacq=wgtfacq wrdfac=wrdfac dts=dts

    zb2smob2sp = zb2sm / zb2sp

    ! Precompute u and v at half levels for efficiency improvement
    ! This is done by linear interpolation in height
    DO  k = MAX(2,kzflat), ke
!CDIR ON_ADB(u)
      zuhl(istart-1:iend,jstart-1:jend,k)       =                    &
        u(istart-1:iend,jstart-1:jend,k,nnew)   *                    &
        wgtfac_u(istart-1:iend,jstart-1:jend,k) +                    &
        u(istart-1:iend,jstart-1:jend,k-1,nnew) *                    &
        (1.0_wp-wgtfac_u(istart-1:iend,jstart-1:jend,k))
!CDIR ON_ADB(v)
      zvhl(istart-1:iend,jstart-1:jend,k)       =                    &
        v(istart-1:iend,jstart-1:jend,k,nnew)   *                    &
        wgtfac_v(istart-1:iend,jstart-1:jend,k) +                    &
        v(istart-1:iend,jstart-1:jend,k-1,nnew) *                    &
        (1.0_wp-wgtfac_v(istart-1:iend,jstart-1:jend,k))
    ENDDO
    SELECT CASE (itype_bbc_w)
    CASE(0,2,4)
      ! Quadratic extrapolation for bottom level if extrapolated wind
      ! is also used for the bottom boundary condition of w
      DO j = jstart-1, jend
!CDIR ON_ADB(u)
!CDIR ON_ADB(v)
        DO i = istart-1, iend
          zuhl(i,j,ke1) = u(i,j,ke  ,nnew)*wgtfacq_u(i,j,1)            &
                        + u(i,j,ke-1,nnew)*wgtfacq_u(i,j,2)            &
                        + u(i,j,ke-2,nnew)*wgtfacq_u(i,j,3)
          zvhl(i,j,ke1) = v(i,j,ke  ,nnew)*wgtfacq_v(i,j,1)            &
                        + v(i,j,ke-1,nnew)*wgtfacq_v(i,j,2)            &
                        + v(i,j,ke-2,nnew)*wgtfacq_v(i,j,3)

          ! Limit extrapolated value to the absolute value at the lowermost main 
          ! level; this is needed to avoid spurious overshoots for underresolved
          ! katabatic flows
!!$ UB>>          IF (ltur .AND. .NOT. lfreeslip_sfc) THEN
          IF (ltur .AND. .NOT. lnosurffluxes_m) THEN
            IF (u(i,j,ke,nnew) > 0.0_wp) THEN
              zuhl(i,j,ke1) = MIN(zuhl(i,j,ke1),u(i,j,ke,nnew))
            ELSE
              zuhl(i,j,ke1) = MAX(zuhl(i,j,ke1),u(i,j,ke,nnew))
            ENDIF
            IF (v(i,j,ke,nnew) > 0.0_wp) THEN
              zvhl(i,j,ke1) = MIN(zvhl(i,j,ke1),v(i,j,ke,nnew))
            ELSE
              zvhl(i,j,ke1) = MAX(zvhl(i,j,ke1),v(i,j,ke,nnew))
            ENDIF
          ENDIF
        ! Ensure that extrapolation does not overshoot to opposite sign
        ! In that case, surface wind is set to zero
        IF (zuhl(i,j,ke1)*u(i,j,ke,nnew) <= 0._wp) &
          zuhl(i,j,ke1) = 0._wp
        IF (zvhl(i,j,ke1)*v(i,j,ke,nnew) <= 0._wp) &
          zvhl(i,j,ke1) = 0._wp
        ENDDO
      ENDDO
    CASE(1,3,5)
      ! Linear extrapolation for bottom level if u/v(ke) is done here
      ! but IS NOT USED for the bottom boundary condition of w
      DO j = jstart-1, jend
!CDIR ON_ADB(u)
!CDIR ON_ADB(v)
        DO i = istart-1, iend
          zuhl(i,j,ke1) = u(i,j,ke,nnew) + wgtfac_u(i,j,ke1) &
                        * (u(i,j,ke,nnew) - u(i,j,ke-1,nnew))
          zvhl(i,j,ke1) = v(i,j,ke,nnew) + wgtfac_v(i,j,ke1) &
                        * (v(i,j,ke,nnew) - v(i,j,ke-1,nnew))
        ENDDO
      ENDDO
    CASE DEFAULT
      yzerrmsg="this case for itype_bbc_w is not available if itype_fast_waves=1!"
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'bottom_BC_for_w')
    END SELECT

    ! Linear extrapolation for top level (typically not used because kzflat>1)
    IF (kzflat == 1) THEN
      zuhl(istart-1:iend,jstart-1:jend,1)        =                    &
        u(istart-1:iend,jstart-1:jend,1,nnew)    *                    &
        (1.0_wp+wgtfac_u(istart-1:iend,jstart-1:jend,1)) -        &
        u(istart-1:iend,jstart-1:jend,2,nnew)    *                    &
        wgtfac_u(istart-1:iend,jstart-1:jend,1)
      zvhl(istart-1:iend,jstart-1:jend,1)        =                    &
        v(istart-1:iend,jstart-1:jend,1,nnew)    *                    &
        (1.0_wp+wgtfac_v(istart-1:iend,jstart-1:jend,1)) -        &
        v(istart-1:iend,jstart-1:jend,2,nnew)    *                    &
        wgtfac_v(istart-1:iend,jstart-1:jend,1)
    ENDIF

    DO k = 1, ke
      
      kp1 = MIN( ke, k+1 )
      km1 = MAX( 1 , k-1 )
!       zfact = 1.0_wp
!       IF ( k==1 .OR. k==ke ) zfact = 0.5_wp
      IF ( k < kzflat ) THEN
        DO j = jstart, jend
!CDIR ON_ADB(u)
!CDIR ON_ADB(v)
!CDIR ON_ADB(zpi)
          DO i = istart, iend
            zpi(i,j,k) = zfx(j)*( u(i,j,k,nnew) - u(i-1,j,k,nnew) )        &
                       + zfydn(j)*v(i,j,k,nnew) - zfyds(j)*v(i,j-1,k,nnew)
          ENDDO
        ENDDO
      ELSE
        DO j = jstart, jend
!CDIR ON_ADB(u)
!CDIR ON_ADB(v)
!CDIR ON_ADB(zdzdx)
!CDIR ON_ADB(zdzdy)
!CDIR ON_ADB(zpi)
!CDIR ON_ADB(zuhl)
!CDIR ON_ADB(zvhl)
          DO i = istart, iend
            zpi(i,j,k) =                                                     &
              zfx(j)*( u(i,j,k,nnew) - u(i-1,j,k,nnew)                       &
              + ( zuhl(i-1,j,k+1)-zuhl(i-1,j,k) )*zdzdx(i-1,j,k)             &
              + ( zuhl(i  ,j,k+1)-zuhl(i  ,j,k) )*zdzdx(i  ,j,k) )           &
              + zfydn(j)*(  v(i,j  ,k,nnew)                                  &
              + ( zvhl(i,j  ,k+1)-zvhl(i,j  ,k) )*zdzdy(i,j  ,k) )           &
              + zfyds(j)*(- v(i,j-1,k,nnew)                                  &
              + ( zvhl(i,j-1,k+1)-zvhl(i,j-1,k) )*zdzdy(i,j-1,k) )
          ENDDO
        ENDDO
      ENDIF
      
      IF (itheta_adv == 0) THEN
        DO j = jstart, jend
!CDIR ON_ADB(zpi)
          DO  i = istart, iend

            ! Calculation of the explicit part of the 3-d divergence
            zwdiv = ( w(i,j,k+1,nnew) - w(i,j,k,nnew) )*zb2smob2sp

            ! 3-d divergence for temperature
            ztdiv(i,j,k) = ( zwdiv * zct1(i,j,k) - zpi(i,j,k) * zct1_h(i,j,k) ) &
                       * dts

            ! Store the 3-d divergence on zpi
            zpi(i,j,k)   = ( zwdiv * zcp1(i,j,k) - zpi(i,j,k) * zcp1_h(i,j,k) ) &
                       * dts
          ENDDO
        ENDDO
      ELSE
        DO j = jstart, jend
!CDIR ON_ADB(zpi)
          DO  i = istart, iend

            ! Calculation of the explicit part of the 3-d divergence
            zwdiv = ( w(i,j,k+1,nnew) - w(i,j,k,nnew) )*zb2smob2sp

            ztdiv(i,j,k) = 0._wp   ! disappears with theta-advection

            ! Store the 3-d divergence on zpi
            zpi(i,j,k)   = ( zwdiv * zcp1(i,j,k) - zpi(i,j,k) * zcp1_h(i,j,k) ) &
                       * dts
          ENDDO
        ENDDO
      ENDIF
    ENDDO

    zb2gmob2gp = zb2gm / zb2gp

    loop_from_south_to_north_2:  DO  j = jstart, jend
    !------------------------------------------------

      IF (itheta_adv == 0) THEN
        DO   k = 1, ke
!CDIR ON_ADB(zppten)
!CDIR ON_ADB(ztpten)
          DO  i = istart, iend

            zwavg_s = ( w(i,j,k+1,nnew)+w(i,j,k,nnew) )*zb2gmob2gp

            ! Explicit part of total pressure tendency
            zppten(i,k) = sppten(i,j,k) + zpi(i,j,k)     &
                        + dts * zcp2(i,j,k) * zwavg_s

            ! Explicit part of total temperature tendency
            ztpten(i,k) = stpten(i,j,k) + ztdiv(i,j,k)   &
                        + dts * zct2(i,j,k) * zwavg_s

          ENDDO
        ENDDO
      ELSE IF (itheta_adv == 1) THEN
        DO k = 1, ke
!CDIR ON_ADB(zppten)
!CDIR ON_ADB(ztpten)
          DO  i = istart, iend

            zwavg_s = ( w(i,j,k+1,nnew)+w(i,j,k,nnew) )*zb2gmob2gp

            ! Explicit part of total pressure tendency
            zppten(i,k) = sppten(i,j,k) + zpi(i,j,k)     &
                        + dts * zcp2(i,j,k) * zwavg_s

            ! Explicit part of total temperature tendency
            ztpten(i,k) = stpten(i,j,k) + dts * zct2(i,j,k) * zwavg_s

          ENDDO
        ENDDO
      ELSE
        DO k = 1, ke
!CDIR ON_ADB(zppten)
!CDIR ON_ADB(ztpten)
          DO  i = istart, iend

            zwavg_s = ( w(i,j,k+1,nnew)+w(i,j,k,nnew) )*zb2gmob2gp

            ! Explicit part of total pressure tendency
            zppten(i,k) = sppten(i,j,k) + zpi(i,j,k)     &
                        + dts * zcp2(i,j,k) * zwavg_s

            ! Explicit part of total temperature tendency
            ztpten(i,k) = stpten(i,j,k)

          ENDDO
        ENDDO
      ENDIF

      ! Calculate the righthand side of the matrix eqation
!!!          IF ( k > 1 ) THEN
        IF (itheta_adv /= 2) THEN
          DO   k = 2, ke
!CDIR ON_ADB(zppten)
!CDIR ON_ADB(ztpten)
!CDIR ON_ADB(pp)
!CDIR ON_ADB(t)
            DO  i = istart, iend

              zd(i,j,k) = w(i,j,k,nnew) + swten(i,j,k)                 &
                + dts * (   zbalphab(i,j,k) * zppten(i,k  )            &
                          + zalphab(i,j,k)  * pp(i,j,k  ,nnew)         &
                          - zbalphat(i,j,k) * zppten(i,k-1)            &
                          - zalphat(i,j,k)  * pp(i,j,k-1,nnew) )

              zd(i,j,k) = zd(i,j,k)                                    &
                + dts * (   zbetab(i,j,k) * (   zbgp*ztpten(i,k  )     &
                                              + t(i,j,k  ,nnew) )      &
                +           zbetat(i,j,k) * (   zbgp*ztpten(i,k-1)     &
                                              + t(i,j,k-1,nnew) ) )

            ENDDO
          ENDDO
        ELSE
          DO   k = 2, ke
!CDIR ON_ADB(zppten)
!CDIR ON_ADB(ztpten)
!CDIR ON_ADB(pp)
!CDIR ON_ADB(t)
            DO  i = istart, iend

              zd(i,j,k) = w(i,j,k,nnew) + swten(i,j,k)                 &
                + dts * (   zbalphab(i,j,k) * zppten(i,k  )            &
                          + zalphab(i,j,k)  * pp(i,j,k  ,nnew)         &
                          - zbalphat(i,j,k) * zppten(i,k-1)            &
                          - zalphat(i,j,k)  * pp(i,j,k-1,nnew) )

              zd(i,j,k) = zd(i,j,k)                                    &
                + dts * (   zbetab(i,j,k) * (   zbgp*ztpten(i,k  )     &
                          + (t(i,j,k  ,nnew) - t0(i,j,k)) )            &
                +           zbetat(i,j,k) * (   zbgp*ztpten(i,k-1)     &
                          + (t(i,j,k-1,nnew) - t0(i,j,k-1)) ) )

            ENDDO
          ENDDO
        ENDIF
!!!          ENDIF


      DO   k = 1, ke
!CDIR ON_ADB(zppten)
!CDIR ON_ADB(ztpten)
!CDIR ON_ADB(pp)
!CDIR ON_ADB(t)
        DO  i = istart, iend

          ! Explicit part of pp(n+1)
          pp(i,j,k,nnew) = pp(i,j,k,nnew) + zppten(i,k)

          ! Explicit part of tp(n+1)
          t (i,j,k,nnew) = t (i,j,k,nnew) + ztpten(i,k)

        ENDDO
      END DO

      ! Set the upper (w=0) and bottom (w=u*dz/dx) boundary conditions on w
      IF (llm) THEN
        DO  i = istart, iend
          w (i,j,  1,nnew) = 0.0_wp
          w (i,j,ke1,nnew) = 0.0_wp
          zd(i,j,ke1)      = w(i,j,ke1,nnew)
        ENDDO
      ENDIF

    ENDDO loop_from_south_to_north_2

    IF ( .NOT.llm ) THEN

      SELECT CASE (itype_bbc_w) 
      CASE (0)
        IF ( iadv_order /= 5 ) THEN
          CALL w_bbc_rk( hhl(:,:,ke1), zuhl(:,:,ke1), zvhl(:,:,ke1), &
                         w(:,:,ke1,nnew), zfx(:), zfyd(:) )
        ELSE
          CALL w_bbc_rk_up5( hhl(:,:,ke1), zuhl(:,:,ke1), zvhl(:,:,ke1), &
                             w(:,:,ke1,nnew), zfx(:), zfyd(:) )
        END IF
      CASE (1)
        IF ( iadv_order /= 5 ) THEN
          CALL w_bbc_rk( hhl(:,:,ke1), u(:,:,ke,nnew), v(:,:,ke,nnew), &
                         w(:,:,ke1,nnew), zfx(:), zfyd(:) )
        ELSE
          CALL w_bbc_rk_up5( hhl(:,:,ke1), u(:,:,ke,nnew), v(:,:,ke,nnew), &
                             w(:,:,ke1,nnew), zfx(:), zfyd(:) )
        END IF
      CASE(2)
        CALL w_bbc_var(zuhl(:,:,ke1), zvhl(:,:,ke1), w(:,:,:,nnew), zfx, zfyd)
      CASE(3)
        CALL w_bbc_var(u(:,:,ke,nnew), v(:,:,ke,nnew), w(:,:,:,nnew), zfx, zfyd)
      CASE(4)
        CALL w_bbc_cen4(zuhl(:,:,ke1), zvhl(:,:,ke1), w(:,:,:,nnew), zfx, zfyd)
      CASE(5)
        CALL w_bbc_cen4(u(:,:,ke,nnew), v(:,:,ke,nnew), w(:,:,:,nnew), zfx, zfyd)
      END SELECT

      DO  j = jstart, jend
        DO  i = istart, iend
          zd(i,j,ke1)     = w(i,j,ke1,nnew)
          w(i,j,1,nnew)   = 0.0_wp
        ENDDO
      ENDDO

    ENDIF

    DO  j = jstart, jend

      ! Up-down algorithm to solve the tridiagonal matrix equation for w

      DO  k = ke, 2, -1
!CDIR ON_ADB(zd)
        DO  i = istart, iend
          zd(i,j,k) = ( zd(i,j,k) - zd(i,j,k+1)*zc(i,j,k) )*zb(i,j,k)
        ENDDO
      ENDDO

      DO  k = 1, ke-1
!CDIR ON_ADB(w)
        DO  i = istart, iend
          w(i,j,k+1,nnew) = za(i,j,k+1)*w(i,j,k,nnew) + zd(i,j,k+1)
        ENDDO
      ENDDO

    ENDDO

    ! Rayleigh damping on w according to Klemp et al. (2008; MWR)
    IF (lspubc .AND. (itype_spubc == 3) ) THEN
      DO k = 1, ke_rd
        DO  j = jstart, jend
          DO  i = istart, iend
            ! wrdfac is defined on top for efficiency reasons
            w(i,j,k,nnew) = w(i,j,k,nnew)*wrdfac(k)  
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    ! Calculation of pp(n+1) and t(n+1)
    ! Completition of the divergence damping term zpi

    ! soundwaves and gravity-waves are treated (partly) implicit
    ! in p*- and T*-equation
    IF (itheta_adv == 0) THEN
      DO k = 1, ke
        DO  j = jstart, jend
!CDIR ON_ADB(w)
          DO  i = istart, iend
            zwdz_s          = w(i,j,k+1,nnew) - w(i,j,k,nnew)
            zwavg_s         = w(i,j,k+1,nnew) + w(i,j,k,nnew)
            ! pp(n+1)
            zwdiv           = dts * zcp1(i,j,k) * zwdz_s
            pp (i,j,k,nnew) = pp(i,j,k,nnew) + zwdiv         &
                              + dts * zcp2(i,j,k) * zwavg_s
            zpi(i,j,k)      = zpi(i,j,k) + zwdiv
            ! tp(n+1)
            zwdiv           = dts * zct1(i,j,k) * zwdz_s
            t  (i,j,k,nnew) = t(i,j,k,nnew)  + zwdiv         &
                              + dts * zct2(i,j,k) * zwavg_s
          ENDDO
        ENDDO
      ENDDO
    ELSE IF (itheta_adv == 1) THEN
      DO k = 1, ke
        DO  j = jstart, jend
!CDIR ON_ADB(w)
          DO  i = istart, iend
            zwdz_s          = w(i,j,k+1,nnew) - w(i,j,k,nnew)
            zwavg_s         = w(i,j,k+1,nnew) + w(i,j,k,nnew)
            ! pp(n+1)
            zwdiv           = dts * zcp1(i,j,k) * zwdz_s
            pp (i,j,k,nnew) = pp(i,j,k,nnew) + zwdiv         &
                            + dts * zcp2(i,j,k) * zwavg_s
            zpi(i,j,k)      = zpi(i,j,k) + zwdiv
            ! tp(n+1)
            t  (i,j,k,nnew) = t(i,j,k,nnew) + dts * zct2(i,j,k) * zwavg_s
          ENDDO
        ENDDO
      ENDDO
    ELSE IF (itheta_adv == 2) THEN
      DO k = 1, ke
        DO  j = jstart, jend
!CDIR ON_ADB(w)
          DO  i = istart, iend
            zwdz_s          = w(i,j,k+1,nnew) - w(i,j,k,nnew)
            zwavg_s         = w(i,j,k+1,nnew) + w(i,j,k,nnew)
            ! pp(n+1)
            zwdiv           = dts * zcp1(i,j,k) * zwdz_s
            pp (i,j,k,nnew) = pp(i,j,k,nnew) + zwdiv         &
                            + dts * zcp2(i,j,k) * zwavg_s
            zpi(i,j,k)      = zpi(i,j,k) + zwdiv
          ENDDO
        ENDDO
      ENDDO
    ENDIF


    IF (izts < nsmsteps) THEN
      ! Addition of 3-d divergence to perturbation pressure
      ! for divergence damping
      DO  k = 1, ke
        zpi(istart:iend,jstart:jend,k) =                                 &
                   pp(istart:iend,jstart:jend,k,nnew)                    &
            + xkd*zpi(istart:iend,jstart:jend,k)
      ENDDO

      !$ser savepoint FastWavesRKUnittest.DivergenceWPPTP-out LargeTimeStep=ntstep RKStageNumber=irk SmallTimeStep=izts
      !$ser verbatim ! recalculate bottom boundary condition for serialization
      !$ser verbatim CALL w_bbc_var(zuhl(:,:,ke1), zvhl(:,:,ke1), zphl(:,:,:), zfx, zfyd)
      !$ser data w=w(:,:,:,nnew) pp=pp(:,:,:,nnew) t=t(:,:,:,nnew)
      !$ser data zpi=zpi wbbc=zphl(:,:,ke1)

      ! Exchange one boundary line of zpi
      IF (ltime) THEN
        CALL get_timings (i_fast_waves, ntstep, dt, izerror)
        IF (ltime_barrier) THEN
          CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
          CALL get_timings (i_fast_waves_barrier, ntstep, dt, izerror)
        ENDIF
      ENDIF
      kzdims(1:24)=(/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                 &
           ( 0,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
            ie, je, kzdims, jstartpar, jendpar, 1, nboundlines,             &
            my_cart_neigh, lperi_x, lperi_y, l2dim,                         &
            10000+nexch_tag, .FALSE.   , ncomm_type, izerror, yzerrmsg,     &
            zpi(:,:,:) )

      ! NOTE: zpi still carries one line of BC values which have been initialized
      !       before the small timestep loop, so nothing has to be done here

      IF (ltime) CALL get_timings (i_fast_waves_comm, ntstep, dt, izerror)

    ENDIF

    !$ser verbatim ! serialize after last of small timesteps
    !$ser verbatim IF (izts == nsmsteps) THEN
    !$ser savepoint FastWavesRKUnittest.DoSmallStep-out LargeTimeStep=ntstep RKStageNumber=irk SmallTimeStep=izts
    !$ser data u=u(:,:,:,nnew) v=v(:,:,:,nnew) w=w(:,:,:,nnew) pp=pp(:,:,:,nnew) t=t(:,:,:,nnew)
    !$ser verbatim ENDIF

    IF ( .NOT.( lperi_x .OR. lperi_y ) ) THEN
      
      ! radiative lateral boundary condition
      ! -----------------------------------------------
      IF ( lradlbc_u ) &
        CALL radiative_lbc (zu_lbx(:,:,:), zu_lby(:,:,:), u(:,:,:,nnew),     &
             u(:,:,:,nnow), v(:,:,:,nnow), ke, 'u', dts)
      IF ( lradlbc_v ) &
        CALL radiative_lbc (zv_lbx(:,:,:), zv_lby(:,:,:), v(:,:,:,nnew),     &
             u(:,:,:,nnow), v(:,:,:,nnow), ke, 'v', dts)
      IF ( lradlbc_w ) &
        CALL radiative_lbc (zw_lbx(:,:,:), zw_lby(:,:,:), w(:,:,:,nnew),     &
             u(:,:,:,nnow), v(:,:,:,nnow), ke1, 'w', dts)
      IF ( lradlbc_pp ) &
        CALL radiative_lbc (zpp_lbx(:,:,:), zpp_lby(:,:,:), pp(:,:,:,nnew),  &
             u(:,:,:,nnow), v(:,:,:,nnow), ke, 's', 0.0_wp)
      IF ( lradlbc_t ) &
        CALL radiative_lbc (zt_lbx(:,:,:), zt_lby(:,:,:), t(:,:,:,nnew),     &
             u(:,:,:,nnow), v(:,:,:,nnow), ke, 's', dts)

    END IF
    
    !$ser savepoint TimeIntegratorUnittest.DoSmallStep-out LargeTimeStep=ntstep RKStageNumber=irk SmallTimeStep=izts
    !$ser data u=u(:,:,:,nnew) v=v(:,:,:,nnew) w=w(:,:,:,nnew) pp=pp(:,:,:,nnew) t=t(:,:,:,nnew)

!------------------------------------------------------------------------------
!     ***************************************************************
!     *   end of time integration using small time steps            *
!     ***************************************************************
!------------------------------------------------------------------------------

  ENDDO loop_over_small_timesteps

  
  IF ( .NOT.( lperi_x .OR. lperi_y) ) THEN
    
    ! radiative lateral boundary condition
    IF ( lradlbc_u )  DEALLOCATE( zu_lbx,  zu_lby,  STAT=izstatd )
    IF ( lradlbc_v )  DEALLOCATE( zv_lbx,  zv_lby,  STAT=izstatd )
    IF ( lradlbc_w )  DEALLOCATE( zw_lbx,  zw_lby,  STAT=izstatd )
    IF ( lradlbc_pp ) DEALLOCATE( zpp_lbx, zpp_lby, STAT=izstatd )
    IF ( lradlbc_t )  DEALLOCATE( zt_lbx,  zt_lby,  STAT=izstatd )
    
  END IF


  IF ( irk == irk_order .OR. llin_tdyn_rk ) THEN
    DEALLOCATE(                &
      zrhoqx_i, zrhoqy_i,      &
      zcp1_h, zct1_h,          &
      zcp1, zct1,              &
      zalphab, zalphat,        &
      zbalphab, zbalphat,      &
      zbetab, zbetat,          &
      STAT=izstatd )
    IF ( ldyn_bbc ) DEALLOCATE( zcwp, STAT=izstatd )
  END IF

  DEALLOCATE( zd, STAT=izstatd )

!------------------------------------------------------------------------------
!  Section 5a: Set free-slip lateral boundary conditions on w
!              (also for periodic boundary conditions; these are overwritten later)
!              (applies also for single-processor and/or 2D case)
!------------------------------------------------------------------------------

  ! Apply a freeslip (zero-gradient) boundary condition to w
  IF ( lw_freeslip ) THEN
    CALL lbc_zerograd( w(:,:,:,nnew), ie, je, ke1, istart, iend, jstart, jend, 1, ke1)
  ENDIF

!------------------------------------------------------------------------------
!  Section 5b: Set periodic lateral boundary conditions on w and all
!              all other variables if required
!------------------------------------------------------------------------------

  ! Boundary exchange:
  IF ( lperi_x .OR. lperi_y .OR. l2dim) THEN
    IF (ltime) THEN
      CALL get_timings (i_fast_waves, ntstep, dt, izerror)
      IF (ltime_barrier) THEN
        CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
        CALL get_timings (i_fast_waves_barrier, ntstep, dt, izerror)
      ENDIF
    ENDIF
    
    ! NOTE: Is the exchange of qv_new and qc_new here really needed?
    kzdims(1:24)=(/ke1,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                     &
         ( 0      ,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,&
          ie, je, kzdims, jstartpar, jendpar, 2, nboundlines,                 &
          my_cart_neigh, lperi_x, lperi_y, l2dim,                             &
          15000+nexch_tag, .FALSE.   , ncomm_type, izerror, yzerrmsg,         &
          w (:,:,:,nnew), t (:,:,:,nnew), qv_new(:,:,:), qc_new(:,:,:),       &
          pp(:,:,:,nnew) )
    IF (ltime) CALL get_timings (i_fast_waves_comm, ntstep, dt, izerror)
    
    ! NOTE: u and v not required here since they have already been updated
    !       inside the small timestep loop in caes of periodic boundary conditions
    ! NOTE: Why is there a exchg_boundaries() call here??? The exchg_boundaries() is
    !       also done after returning to src_runge_kutta.f90 for u,v,w,t,pp!!! Why
    !       is there an exchangen for qv and qc here, if they are not modified???

  END IF

  IF ( llin_tdyn_rk ) THEN
    ! Compute density of moist air for start values of tp and pp
    ! ----------------------------------------------------------------
    CALL calrho_tp_pp( t(:,:,:,nnew), pp(:,:,:,nnew),    &
      qv_now(:,:,:), qc_now(:,:,:), qrs,                 &
      t0, p0, rho, ie, je, ke, r_d, rvd_m_o )
  END IF

  IF (ltime) CALL get_timings (i_fast_waves, ntstep, dt, izerror)

  IF   ( irk == irk_order .AND. ntstep == nstop) THEN
    DEALLOCATE( zfx, zfydn, zfyd, zfyds, zcp2, zct2, zdzdx, zdzdy,   &
                wrdfac, STAT=izstatd )
  END IF

!------------------------------------------------------------------------------
!  End of the module procedure fast_waves_runge_kutta
!------------------------------------------------------------------------------

END SUBROUTINE fast_waves_runge_kutta


!==============================================================================
!==============================================================================

SUBROUTINE w_bbc_rk( q_in, u_sfc, v_sfc, q_out, fadux, fadvy )

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "fast_waves_rk" computes the bottom
!   boundary for the vertical velocity in the procedure fast_waves_runge_kutta.
!
! Method:
!   The same time integration and advection scheme as in src_runge_kutta 
!   is used to compute w=u*dz/dx.
!
!------------------------------------------------------------------------------

! Declarations:

! Subroutine arguments:
! ---------------------
REAL (KIND=wp),     INTENT(IN) :: &
  q_in(ie,je),     & ! input field which is advected ( normally hhl )
  u_sfc(ie,je),    & ! zonal wind component at surface level
  v_sfc(ie,je),    & ! meridional wind component at surface level
  fadux(je),       & !
  fadvy(je)          !

REAL (KIND=wp),     INTENT(OUT) :: &
  q_out(ie,je)       ! output field: bbc for w ( = u*dz/dx )

! Local variables:
! ----------------

INTEGER (KIND=iintegers) ::  &
  kzdims(24),      & ! Vertical dimensions for exchg_boundaries
  i, j, irk_bbc      !  Loop indices

REAL (KIND=wp)     :: &
  zqten(ie,je),    & ! field to store the tendency
  q_tmp(ie,je,1),  & ! temporary working array
  zui(ie,je),      & ! advection velocity in x-direction
  zvi(ie,je)         ! advection velocity in y-direction

REAL (KIND=wp)     :: &
  zdtrk, &              ! time step in rk-loop
  zsign

INTEGER (KIND=iintegers) ::  izerror
CHARACTER (LEN=80)       ::  yzerrmsg

!------------------------------------------------------------------------------

  izerror  = 0
  yzerrmsg = '   '
  kzdims(:)= 0_iintegers
  zsign = SIGN(1._wp,dt)

  ! compute advection velocities
  DO  j = jstart, jend
!CDIR ON_ADB(u_sfc)
!CDIR ON_ADB(v_sfc)
    DO  i = istart, iend
      zui(i,j) = 0.5_wp*( u_sfc(i,j) + u_sfc(i-1,j) )
      zvi(i,j) = 0.5_wp*( v_sfc(i,j) + v_sfc(i,j-1) )
    ENDDO
  ENDDO

  ! save input field ( old values )
  q_tmp(:,:,1) = q_in(:,:)

  runge_kutta_loop: DO irk_bbc = 1, irk_order

    IF ( irk_bbc == irk_order ) THEN

      ! last RK-step: save tendency in q_out
      DO  j = jstart, jend
        DO  i = istart, iend
          
          ! u*dz/dx + v*dz/dy
          q_out(i,j) = udsdx( 0, iadv_order, q_tmp(1,1,1),   &
                              ie, je, 1, i, j, 1, im, ip,    &
                              zui(i,j), fadux(j),            &
                              zvi(i,j), fadvy(j), zsign)

        END DO
      END DO

    ELSE

      DO  j = jstart, jend
!CDIR ON_ADB(zqten)
        DO  i = istart, iend

          ! u*dz/dx + v*dz/dy
          zqten(i,j) = udsdx( 0, iadv_order, q_tmp(1,1,1),   &
                              ie, je, 1, i, j, 1, im, ip,    &
                              zui(i,j), fadux(j),            &
                              zvi(i,j), fadvy(j), zsign)

        END DO
      END DO

      IF (irk_bbc == 1 .AND. irk_order == 3) THEN
! JF:           zdtrk = dt/4.0_wp
        zdtrk = dt/3.0_wp
      ELSE
        zdtrk = dt/2.0_wp
      END IF

      ! intermediate RK-steps: update q_tmp
      DO  j = jstart, jend
!CDIR ON_ADB(zqten)
        DO  i = istart, iend
          q_tmp(i,j,1) = q_in(i,j) + zdtrk * zqten(i,j)
        END DO
      END DO

      ! Boundary exchange in parallel mode
      IF (ltime) THEN
        CALL get_timings (i_fast_waves, ntstep, dt, izerror)
        IF (ltime_barrier) THEN
          CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
          CALL get_timings (i_fast_waves_barrier, ntstep, dt, izerror)
        ENDIF
      ENDIF
      
      kzdims(1:24)=(/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                  &
           ( 0,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
            ie, je, kzdims, jstartpar, jendpar, imipmx, nboundlines,         &
            my_cart_neigh, lperi_x, lperi_y, l2dim,                          &
            20000+nexch_tag, .FALSE.   , ncomm_type, izerror, yzerrmsg,      &
            q_tmp(:,:,1) )

      ! NOTE: which type of BC should be applied here for q_tmp?
 
      IF (ltime) CALL get_timings (i_fast_waves_comm, ntstep, dt, izerror)

    ENDIF

  ENDDO runge_kutta_loop

END SUBROUTINE w_bbc_rk

!==============================================================================

SUBROUTINE w_bbc_rk_up5( q_in, u_sfc, v_sfc, q_out, fadux, fadvy )

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "fast_waves_rk" computes the bottom
!   boundary for the vertical velocity in the procedure fast_waves_runge_kutta.
!   This is the version with explicit call of FUNCTION udsdx_up5_xy for
!   5th-UP advection.  
!
! Method:
!   The same time integration and advection scheme as in src_runge_kutta 
!   is used to compute w=u*dz/dx.
!
!------------------------------------------------------------------------------

! Declarations:

! Subroutine arguments:
! ---------------------
REAL (KIND=wp),     INTENT(IN) :: &
  q_in(ie,je),     & ! input field which is advected ( normally hhl )
  u_sfc(ie,je),    & ! zonal wind component at surface level
  v_sfc(ie,je),    & ! meridional wind component at surface level
  fadux(je),       & !
  fadvy(je) 

REAL (KIND=wp),     INTENT(OUT) :: &
  q_out(ie,je)       ! output field: bbc for w ( = u*dz/dx )

! Local variables:
! ----------------

INTEGER (KIND=iintegers) ::  &
  kzdims(24),      & ! Vertical dimensions for exchg_boundaries
  i,  j,  irk_bbc    !  Loop indices

REAL (KIND=wp)     :: &
  zqten(ie,je),    & ! field to store the tendency
  q_tmp(ie,je,1),  & ! temporary working array
  zui(ie,je),      & ! advection velocity in x-direction
  zvi(ie,je)         ! advection velocity in y-direction

REAL (KIND=wp)     :: &
  zdtrk, &              ! time step in rk-loop
  zsign

INTEGER (KIND=iintegers) ::  izerror
CHARACTER (LEN=80)       ::  yzerrmsg

!------------------------------------------------------------------------------

  izerror  = 0
  yzerrmsg = '   '
  kzdims(:)= 0_iintegers
  zsign = SIGN(1._wp,dt)

  ! compute advection velocities
  DO  j = jstart, jend
!CDIR NONEIGHBORS
!CDIR ON_ADB(u_sfc)
!CDIR ON_ADB(v_sfc)
    DO  i = istart, iend
      zui(i,j) = 0.5_wp*( u_sfc(i,j) + u_sfc(i-1,j) )
      zvi(i,j) = 0.5_wp*( v_sfc(i,j) + v_sfc(i,j-1) )
    ENDDO
  ENDDO

  ! save input field ( old values )
  q_tmp(:,:,1) = q_in(:,:)

  runge_kutta_loop_up5: DO irk_bbc = 1, irk_order

    IF ( irk_bbc == irk_order ) THEN

      ! last RK-step: save tendency in q_out
      DO  j = jstart, jend
!CDIR NONEIGHBORS
        DO  i = istart, iend
          
          ! u*dz/dx + v*dz/dy
          q_out(i,j) = udsdx_up5_xy( q_tmp(1,1,1),           &
                                     ie, je, 1, i, j, 1,     &
                                     zui(i,j), fadux(j),     &
                                     zvi(i,j), fadvy(j), zsign)

        END DO
      END DO

    ELSE

      DO  j = jstart, jend
!CDIR NONEIGHBORS
!CDIR ON_ADB(zqten)
        DO  i = istart, iend

          ! u*dz/dx + v*dz/dy
          zqten(i,j) = udsdx_up5_xy( q_tmp(1,1,1),           &
                                     ie, je, 1, i, j, 1,     &
                                     zui(i,j), fadux(j),     &
                                     zvi(i,j), fadvy(j), zsign)

        END DO
      END DO

      IF (irk_bbc == 1 .AND. irk_order == 3) THEN
! JF:           zdtrk = dt/4.0_wp
        zdtrk = dt/3.0_wp
      ELSE
        zdtrk = dt/2.0_wp
      END IF

      ! intermediate RK-steps: update q_tmp
      DO  j = jstart, jend
!CDIR ON_ADB(zqten)
        DO  i = istart, iend
          q_tmp(i,j,1) = q_in(i,j) + zdtrk * zqten(i,j)
        END DO
      END DO

      ! Boundary exchange in parallel mode
      IF (ltime) THEN
        CALL get_timings (i_fast_waves, ntstep, dt, izerror)
        IF (ltime_barrier) THEN
          CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
          CALL get_timings (i_fast_waves_barrier, ntstep, dt, izerror)
        ENDIF
      ENDIF
      
      kzdims(1:24)=(/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                  &
           ( 0,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
            ie, je, kzdims, jstartpar, jendpar, imipmx, nboundlines,         &
            my_cart_neigh, lperi_x, lperi_y, l2dim,                          &
            20000+nexch_tag, .FALSE.   , ncomm_type, izerror, yzerrmsg,      &
            q_tmp(:,:,1) )

      ! NOTE: which type of BC should be applied here for q_tmp?

      IF (ltime) CALL get_timings (i_fast_waves_comm, ntstep, dt, izerror)

    ENDIF

  ENDDO runge_kutta_loop_up5

END SUBROUTINE w_bbc_rk_up5

!==============================================================================

SUBROUTINE w_bbc_cen4(u_sfc, v_sfc, w, fadux, fadvy )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine computes the bottom boundary condition for vertical velocity
!   using fourth-order centred differences
!
! Declarations:
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! ---------------------
REAL (KIND=wp),     INTENT(IN) :: &
  u_sfc(ie,je),    & ! zonal wind component at surface level
  v_sfc(ie,je),    & ! meridional wind component at surface level
  fadux(je),       & !
  fadvy(je)          !

REAL (KIND=wp),     INTENT(OUT) :: &
  w(ie,je,ke1)       ! output field: bbc for w ( = u*dz/dx )

! Local variables:
! ----------------

INTEGER (KIND=iintegers) :: i,  j    !  Loop indices

REAL (KIND=wp)     :: &
  zui(ie,je),      & ! advection velocity in x-direction
  zvi(ie,je)       ! advection velocity in y-direction
REAL (KIND=wp)     :: fac1,fac2

fac1 = 4.0_wp/3.0_wp
fac2 = 1.0_wp/6.0_wp

  ! average horizontal surface wind components to w points
  DO  j = jstart, jend
!CDIR ON_ADB(u_sfc)
!CDIR ON_ADB(v_sfc)
    DO  i = istart, iend
      zui(i,j) = 0.5_wp*( u_sfc(i,j) + u_sfc(i-1,j) )
      zvi(i,j) = 0.5_wp*( v_sfc(i,j) + v_sfc(i,j-1) )
    ENDDO
  ENDDO

  ! diagnose w(sfc) using fourth-order centered differences
  DO  j = jstart, jend
!CDIR ON_ADB(hhl)
    DO  i = istart, iend
      w(i,j,ke1) = 0.5_wp*(fadux(j)*zui(i,j)*(fac1*(hhl(i+1,j,ke1)-hhl(i-1,j,ke1))- &
                   fac2*(hhl(i+2,j,ke1)-hhl(i-2,j,ke1)))+fadvy(j)*zvi(i,j)*             &
                   (fac1*(hhl(i,j+1,ke1)-hhl(i,j-1,ke1))-fac2*(hhl(i,j+2,ke1)-hhl(i,j-2,ke1))))
    ENDDO
  ENDDO


END SUBROUTINE w_bbc_cen4

!==============================================================================

SUBROUTINE w_bbc_var(u_sfc, v_sfc, w, fadux, fadvy )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine computes the bottom boundary condition for vertical velocity
!   using the same differencing scheme as for advection
!
! Declarations:
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! ---------------------
REAL (KIND=wp),     INTENT(IN) :: &
  u_sfc(ie,je),    & ! zonal wind component at surface level
  v_sfc(ie,je),    & ! meridional wind component at surface level
  fadux(je),       & !
  fadvy(je)          !

REAL (KIND=wp),     INTENT(OUT) :: &
  w(ie,je,ke1)       ! output field: bbc for w ( = u*dz/dx )

! Local variables:
! ----------------

INTEGER (KIND=iintegers) :: i,  j    !  Loop indices

REAL (KIND=wp)     :: &
  zui(ie,je),      & ! advection velocity in x-direction
  zvi(ie,je),      & ! advection velocity in y-direction
  zsign              !

  zsign = SIGN(1._wp,dt)

  ! average horizontal surface wind components to w points
  DO  j = jstart, jend
!CDIR ON_ADB(u_sfc)
!CDIR ON_ADB(v_sfc)
    DO  i = istart, iend
      zui(i,j) = 0.5_wp*( u_sfc(i,j) + u_sfc(i-1,j) )
      zvi(i,j) = 0.5_wp*( v_sfc(i,j) + v_sfc(i,j-1) )
    ENDDO
  ENDDO

  ! diagnose w(sfc) using the same differencing scheme as for advection
  DO  j = jstart, jend
!CDIR ON_ADB(hhl)
    DO  i = istart, iend

      w(i,j,ke1) = udsdx( 0, iadv_order, hhl(1:ie,1:je,ke1),   &
                          ie, je, 1, i, j, 1, im, ip,          &
                          zui(i,j), fadux(j),                  &
                          zvi(i,j), fadvy(j), zsign)

    ENDDO
  ENDDO

END SUBROUTINE w_bbc_var

!==============================================================================

SUBROUTINE init_radiative_lbc( q_in, q_lbx, q_lby,             &
                               lb_w, lb_e, lb_s, lb_n, kdim )

!------------------------------------------------------------------------------
!
! Description: save actual values at lateral boundaries for radiative BC
!
!------------------------------------------------------------------------------

! Declarations:

! Subroutine arguments:
! ---------------------
INTEGER (KIND=iintegers), INTENT(IN) :: &
  lb_w,              & ! index of western boundary
  lb_e,              & ! index of eastern boundary
  lb_s,              & ! index of southern boundary
  lb_n,              & ! index of northern boundary
  kdim                 ! vertical dimension of fields

REAL (KIND=wp),     INTENT(IN) :: &
  q_in(ie,je,kdim)     ! input field

REAL (KIND=wp),     INTENT(OUT) :: &
  q_lbx(6,je,kdim),  & ! west and east boundary values
  q_lby(ie,6,kdim)     ! south and north boundary values

!------------------------------------------------------------------------------

  ! save boundary values in output fields

  ! west boundary
  IF ( my_cart_neigh(1) == -1 ) THEN
    q_lbx(1,1:je,1:kdim) = q_in(lb_w,1:je,1:kdim)
    q_lbx(2,1:je,1:kdim) = q_in(lb_w+1,1:je,1:kdim)
  END IF

  
  ! east boundary
  IF ( my_cart_neigh(3) == -1 ) THEN
    q_lbx(5,1:je,1:kdim) = q_in(lb_e-1,1:je,1:kdim)
    q_lbx(6,1:je,1:kdim) = q_in(lb_e,1:je,1:kdim)
  END IF
  
  IF (.NOT. l2dim) THEN
    ! south boundary
    IF ( my_cart_neigh(4) == -1 ) THEN
      q_lby(1:ie,1,1:kdim) = q_in(1:ie,lb_s,1:kdim)
      q_lby(1:ie,2,1:kdim) = q_in(1:ie,lb_s+1,1:kdim)
    END IF
  
    ! north boundary
    IF ( my_cart_neigh(2) == -1 ) THEN
      q_lby(1:ie,5,1:kdim) = q_in(1:ie,lb_n-1,1:kdim)
      q_lby(1:ie,6,1:kdim) = q_in(1:ie,lb_n,1:kdim)
    END IF
  ENDIF

END SUBROUTINE init_radiative_lbc


SUBROUTINE radiative_lbc( q_lbx, q_lby, q_out, u_old, v_old, kdim, var_typ, dts )

!------------------------------------------------------------------------------
!
! Description: radiative lateral boundary conditions for dynamic variables
!
! Method:      radiate / transport perturbations out of domain using a
!              fixed sound wave speed as advective velocity (Durran, 1993)  
!
!------------------------------------------------------------------------------

! Declarations:

! Subroutine arguments:
! ---------------------
INTEGER (KIND=iintegers), INTENT(IN) :: &
  kdim                 ! vertical dimension of fields

CHARACTER(LEN=1), INTENT(IN) :: &
  var_typ              ! 's' for values at scalar positions
                       ! 'u', 'v', 'w' for velocity components

REAL (KIND=wp),     INTENT(IN) :: &
  dts,               & ! small time step
  q_lbx(6,je,kdim),  & ! west and east input boundary values
  q_lby(ie,6,kdim),  & ! south and north input boundary values
  u_old(ie,je,kdim), v_old(ie,je,kdim)  ! Advecting wind speed

REAL (KIND=wp),     INTENT(OUT) :: &
  q_out(ie,je,kdim)    ! output field (updated at lateral boundaries)


! Local variables:
! ----------------
INTEGER (KIND=iintegers) :: &
  i,j,k,k1            !  Loop indices

REAL (KIND=wp)     :: &
  zminspd, zfac1, zfac2, zfac3, zfac4

REAL (KIND=wp)     :: &
  zcs,               & ! assumed speed of propagating gravity wave modes (25 m/s)
  zfadsxw(je,kdim),zfadsxe(je,kdim),     & ! Advecting speeds multiplied by dts and 1/dx
  zfadsys(ie,kdim),zfadsyn(ie,kdim),     & !
  zfaduxw(je,kdim),zfaduxe(je,kdim),     & !
  zfadvxw(je,kdim),zfadvxe(je,kdim),     & !
  zfaduys(ie,kdim),zfaduyn(ie,kdim),     & !
  zfadvys(ie,kdim),zfadvyn(ie,kdim),     & !
  zqouti (ie,kdim),zqoutj (je,kdim)
! sound wave speed
!   zcs  = SQRT( gamma*r_d*303.0_wp )

!  Use phase speed typical for gravity waves
zcs = 25._wp
zminspd = 2._wp

! Vertical staggering for w is not implemented but likely unnecessary
! GZ: changed loop order because the k loop is usually longer than the j loop
!CDIR NOVECTOR
DO  j = jstart , jend
  zfac1 = dts*acrlat(j,1)*eddlon
  DO k = 1, ke
! x-direction - western boundary
    zfadsxw(j,k) = MAX(zcs-u_old(istartu,j,k),zminspd)*zfac1
    zfaduxw(j,k) = zfadsxw(j,k)
! x-direction - eastern boundary
    zfadsxe(j,k) = MAX(zcs+u_old(iendu,j,k),zminspd)*zfac1
    zfaduxe(j,k) = zfadsxe(j,k)
  ENDDO
ENDDO
!CDIR NOVECTOR
DO  j = jstartv , jendv
  zfac2 = dts*acrlat(j,2)*eddlon
  DO k = 1, ke
! x-direction - western boundary
    zfadvxw(j,k) = MAX(zcs-u_old(istartu,j,k),zminspd)*zfac2
! x-direction - eastern boundary
    zfadvxe(j,k) = MAX(zcs+u_old(iendu,j,k),zminspd)*zfac2
  ENDDO
ENDDO

IF (kdim == ke1) THEN
  k = ke1
  k1 = ke
  DO  j = jstart , jend
! x-direction - western boundary
    zfadsxw(j,k) = MAX(zcs-u_old(istartu,j,k1),zminspd)*dts*acrlat(j,1)*eddlon
    zfaduxw(j,k) = zfadsxw(j,k)
! x-direction - eastern boundary
    zfadsxe(j,k) = MAX(zcs+u_old(iendu,j,k1),zminspd)*dts*acrlat(j,1)*eddlon
    zfaduxe(j,k) = zfadsxe(j,k)
  ENDDO
  DO  j = jstartv , jendv
! x-direction - western boundary
    zfadvxw(j,k) = MAX(zcs-u_old(istartu,j,k1),zminspd)*dts*acrlat(j,2)*eddlon
! x-direction - eastern boundary
    zfadvxe(j,k) = MAX(zcs+u_old(iendu,j,k1),zminspd)*dts*acrlat(j,2)*eddlon
  ENDDO
ENDIF

IF (.NOT.l2dim) THEN
  zfac1 = dts*acrlat(jslb,1)*eddlat*crlat(jslb,1)
  zfac2 = dts*acrlat(jslbv,2)*eddlat*crlat(jslbv,2)
  zfac3 = dts*acrlat(jelb,1)*eddlat*crlat(jelb,1)
  zfac4 = dts*acrlat(jelbv,2)*eddlat*crlat(jelbv,2)
  DO k = 1, kdim
    k1 = MIN(k,ke)
    DO  i = istart , iend
! y-direction - southern boundary
      zfadsys(i,k) = MAX(zcs-v_old(i,jstartv,k1),zminspd)*zfac1
      zfadvys(i,k) = MAX(zcs-v_old(i,jstartv,k1),zminspd)*zfac2
! y-direction - northern boundary
      zfadsyn(i,k) = MAX(zcs+v_old(i,jendv,k1),zminspd)*zfac3
      zfadvyn(i,k) = MAX(zcs+v_old(i,jendv,k1),zminspd)*zfac4
    ENDDO
    DO  i = istartu , iendu
! y-direction - southern boundary
      zfaduys(i,k) = MAX(zcs-v_old(i,jstartv,k1),zminspd)*zfac1
! y-direction - northern boundary
      zfaduyn(i,k) = MAX(zcs+v_old(i,jendv,k1),zminspd)*zfac3
    ENDDO
  ENDDO
ENDIF

  ! west boundary
  IF ( my_cart_neigh(1) == -1 ) THEN
    SELECT CASE(var_typ)
    CASE ('s', 'w')
      zqoutj(jstart:jend,1:kdim) = q_lbx(1,jstart:jend,1:kdim) + zfadsxw(jstart:jend,1:kdim) * &
                                  (q_lbx(2,jstart:jend,1:kdim) - q_lbx(1,jstart:jend,1:kdim))
!CDIR NOLOOPCHG
      DO i=1,islb
        q_out(i,jstart:jend,1:kdim) = zqoutj(jstart:jend,1:kdim)
      END DO
    CASE ('u')
      zqoutj(jstart:jend,1:kdim) = q_lbx(1,jstart:jend,1:kdim) + zfaduxw(jstart:jend,1:kdim) * &
                                  (q_lbx(2,jstart:jend,1:kdim) - q_lbx(1,jstart:jend,1:kdim))
!CDIR NOLOOPCHG
      DO i=1,islbu
        q_out(i,jstart:jend,1:kdim) = zqoutj(jstart:jend,1:kdim)
      END DO
    CASE ('v')
      zqoutj(jstartv:jendv,1:kdim) = q_lbx(1,jstartv:jendv,1:kdim) + zfadvxw(jstartv:jendv,1:kdim) * &
                                    (q_lbx(2,jstartv:jendv,1:kdim) - q_lbx(1,jstartv:jendv,1:kdim))
!CDIR NOLOOPCHG
      DO i=1,islb
        q_out(i,jstartv:jendv,1:kdim) = zqoutj(jstartv:jendv,1:kdim)
      END DO
    END SELECT
  END IF

  ! east boundary
  IF ( my_cart_neigh(3) == -1 ) THEN
    SELECT CASE(var_typ)
    CASE ('s', 'w')
      zqoutj(jstart:jend,1:kdim) = q_lbx(6,jstart:jend,1:kdim) - zfadsxe(jstart:jend,1:kdim) * &
                                  (q_lbx(6,jstart:jend,1:kdim) - q_lbx(5,jstart:jend,1:kdim))
!CDIR NOLOOPCHG
      DO i=ielb,ie
        q_out(i,jstart:jend,1:kdim) = zqoutj(jstart:jend,1:kdim)
      END DO
    CASE ('u')
      zqoutj(jstart:jend,1:kdim) = q_lbx(6,jstart:jend,1:kdim) - zfaduxe(jstart:jend,1:kdim) * &
                                  (q_lbx(6,jstart:jend,1:kdim) - q_lbx(5,jstart:jend,1:kdim))
!CDIR NOLOOPCHG
      DO i=ielbu,ie
        q_out(i,jstart:jend,1:kdim) = zqoutj(jstart:jend,1:kdim)
      END DO
    CASE ('v')
      zqoutj(jstartv:jendv,1:kdim) = q_lbx(6,jstartv:jendv,1:kdim) - zfadvxe(jstartv:jendv,1:kdim) * &
                                    (q_lbx(6,jstartv:jendv,1:kdim) - q_lbx(5,jstartv:jendv,1:kdim))
!CDIR NOLOOPCHG
      DO i=ielb,ie
        q_out(i,jstartv:jendv,1:kdim) = zqoutj(jstartv:jendv,1:kdim)
      END DO
    END SELECT
  END IF

IF (.NOT.l2dim) THEN
  ! south boundary
  IF ( my_cart_neigh(4) == -1 ) THEN
    SELECT CASE(var_typ)
    CASE ('s', 'w')
      zqouti(istart:iend,1:kdim) = q_lby(istart:iend,1,1:kdim) + zfadsys(istart:iend,1:kdim) * &
                                  (q_lby(istart:iend,2,1:kdim) - q_lby(istart:iend,1,1:kdim))
!CDIR NOLOOPCHG
      DO j=1,jslb
        q_out(istart:iend,j,1:kdim) = zqouti(istart:iend,1:kdim)
      END DO
    CASE ('u')
      zqouti(istartu:iendu,1:kdim) = q_lby(istartu:iendu,1,1:kdim) + zfaduys(istartu:iendu,1:kdim) * &
                                    (q_lby(istartu:iendu,2,1:kdim) - q_lby(istartu:iendu,1,1:kdim))
!CDIR NOLOOPCHG
      DO j=1,jslb
        q_out(istartu:iendu,j,1:kdim) = zqouti(istartu:iendu,1:kdim)
      END DO
    CASE ('v')
      zqouti(istart:iend,1:kdim) = q_lby(istart:iend,1,1:kdim) + zfadvys(istart:iend,1:kdim) * &
                                  (q_lby(istart:iend,2,1:kdim) - q_lby(istart:iend,1,1:kdim))
!CDIR NOLOOPCHG
      DO j=1,jslbv
        q_out(istart:iend,j,1:kdim) = zqouti(istart:iend,1:kdim)
      END DO
    END SELECT
  END IF

  ! north boundary
  IF ( my_cart_neigh(2) == -1 ) THEN
    SELECT CASE(var_typ)
    CASE ('s', 'w')
      zqouti(istart:iend,1:kdim) = q_lby(istart:iend,6,1:kdim) - zfadsyn(istart:iend,1:kdim) * &
                                  (q_lby(istart:iend,6,1:kdim) - q_lby(istart:iend,5,1:kdim))
!CDIR NOLOOPCHG
      DO j=jelb,je
        q_out(istart:iend,j,1:kdim) = zqouti(istart:iend,1:kdim)
      END DO
    CASE ('u')
      zqouti(istartu:iendu,1:kdim) = q_lby(istartu:iendu,6,1:kdim) - zfaduyn(istartu:iendu,1:kdim) * &
                                    (q_lby(istartu:iendu,6,1:kdim) - q_lby(istartu:iendu,5,1:kdim))
!CDIR NOLOOPCHG
      DO j=jelb,je
        q_out(istartu:iendu,j,1:kdim) = zqouti(istartu:iendu,1:kdim)
      END DO
    CASE ('v')
      zqouti(istart:iend,1:kdim) = q_lby(istart:iend,6,1:kdim) - zfadvyn(istart:iend,1:kdim) * &
                                  (q_lby(istart:iend,6,1:kdim) - q_lby(istart:iend,5,1:kdim))
!CDIR NOLOOPCHG
      DO j=jelbv,je
        q_out(istart:iend,j,1:kdim) = zqouti(istart:iend,1:kdim)
      END DO
    END SELECT
  END IF
ENDIF

  ! north-west corner
  IF ( my_cart_neigh(2) == -1 .AND. my_cart_neigh(1) == -1 ) THEN
    SELECT CASE(var_typ)
    CASE ('s', 'w')
      DO k = 1, kdim
        q_out(islb,jelb,k) = q_lbx(2,jelb-1,k)
      END DO
    CASE ('u')
      DO k = 1, kdim
        q_out(islbu,jelb,k) = q_lbx(2,jelb-1,k)
      END DO
    CASE ('v')
      DO k = 1, kdim
        q_out(islb,jelbv,k) = q_lbx(2,jelbv-1,k)
      END DO
    END SELECT
  END IF

  ! north-east corner
  IF ( my_cart_neigh(2) == -1 .AND. my_cart_neigh(3) == -1 ) THEN
    SELECT CASE(var_typ)
    CASE ('s', 'w')
      DO k = 1, kdim
        q_out(ielb,jelb,k) = q_lbx(5,jelb-1,k)
      END DO
    CASE ('u')
      DO k = 1, kdim
        q_out(ielbu,jelb,k) = q_lbx(5,jelb-1,k)
      END DO
    CASE ('v')
      DO k = 1, kdim
        q_out(ielb,jelbv,k) = q_lbx(5,jelbv-1,k)
      END DO
    END SELECT
  END IF

  ! south-west corner
  IF ( my_cart_neigh(4) == -1 .AND. my_cart_neigh(1) == -1 ) THEN
    SELECT CASE(var_typ)
    CASE ('s', 'w')
      DO k = 1, kdim
        q_out(islb,jslb,k) = q_lbx(2,jslb+1,k)
      END DO
    CASE ('u')
      DO k = 1, kdim
        q_out(islbu,jslb,k) = q_lbx(2,jslb+1,k)
      END DO
    CASE ('v')
      DO k = 1, kdim
        q_out(islb,jslbv,k) = q_lbx(2,jslbv+1,k)
      END DO
    END SELECT
  END IF

  ! south-east corner
  IF ( my_cart_neigh(4) == -1 .AND. my_cart_neigh(3) == -1 ) THEN
    SELECT CASE(var_typ)
    CASE ('s', 'w')
      DO k = 1, kdim
        q_out(ielb,jslb,k) = q_lbx(5,jslb+1,k)
      END DO
    CASE ('u')
      DO k = 1, kdim
        q_out(ielbu,jslb,k) = q_lbx(5,jslb+1,k)
      END DO
    CASE ('v')
      DO k = 1, kdim
        q_out(ielb,jslbv,k) = q_lbx(5,jslbv+1,k)
      END DO
    END SELECT
  END IF

END SUBROUTINE radiative_lbc

!==============================================================================

!------------------------------------------------------------------------------
!  End of the module fast_waves_rk
!------------------------------------------------------------------------------

END MODULE fast_waves_rk
