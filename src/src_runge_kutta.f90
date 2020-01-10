!+ Source module for 2-timelevel Runge-Kutta version of the dynamics
!------------------------------------------------------------------------------

MODULE src_runge_kutta

!------------------------------------------------------------------------------
!
! Description:
!   The module "src_runge_kutta" performs one time step of the integration
!   of the spatially discretized thermodynamical equations using one of  
!   several selectable 2-timelevel Runge-Kutta schemes.
!   Normal RK-variants of 2nd- or 3rd-order (in principle it is also possible
!   to use a simple 1st-order Euler forward scheme) as well as a slightly
!   more sophisticated 3rd-order TVD-Runge-Kutta (total variation diminishing)
!   scheme are implemented.
!   The type of scheme is set using the namelist-Parameter irunge_kutta
!   (1: normal RK-scheme, 2: TVD-RK-scheme) and
!   the order of the scheme (makes only sense for case 1) is set with
!   the namelist-Parameter irk_order (value: 1, 2 or 3).
!   In addition, the order of the operator for the horizontal advection of
!   the variables u, v, w, pp an T is choosen via the namelist-Parameter
!   iadv_order (value: 3, 4, 5 or 6).
!   Best choice COSMO_DE: irunge_kutta=1, irk_order=3, iadv_order=5.
!               COSMO_EU: irunge_kutta=1, irk_order=3, iadv_order=3.
!   In contrary to Almuts scheme the RK-loop contains the effects
!   of the complete slow tendencies - especially the whole 3D advection is
!   computed in each RK-step - and the integration of the small time steps
!   done in the routine "fast_waves_runge_kutta".
!   Driving routine is the model procedure "org_runge_kutta", which
!   calls the routines required. 
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
!  Technical changes and some optimizations
! 3.13       2004/12/03 Jochen Foerstner
!  Adaptations for calls to boundary exchange; more modifications
!  Modifications to run with latent heat nudging (Klaus Stephan et.al.)
! 3.14       2005/01/25 Jochen Foerstner
!  Enhanced formulation of moist turbulence for 3D turbulence scheme
! 3.16       2005/07/22 Jochen Foerstner
!  New formulation of advection for the moisture quantities
!  (Courant number independent Eulerian and Semi-Lagrangian variants)
!  Optimizations
! 3.18       2006/03/03 Jochen Foerstner / Klaus Stephan
!  Conservation form of Bott Scheme is now standard;
!  Changed interfaces and some additional optimizations
!  LHN namelist parameter moved to data_lheat_nudge to avoid to many dependencies
! 3.19       2006/04/25 Michael Baldauf
!  Corrections for computing small time step
! 3.21       2006/12/04 Jochen Foerstner / Ulrich Schaettler
!  More updates
!  Added additional terms for deep atmosphere (R. Petrik)
!  Put subroutines for slow tendencies (to src_slow_tendencies_rk) and 
!  advection (to src_advection_rk) to extra modules
! V3_23        2007/03/30 Jochen Foerstner, Michael Baldauf, Lucio Torrisi
!  Bug correction for computing CFL criterion in SR check_cfl_horiz_advection
!  Apply Rayleigh damping only, where necessary
!  Eliminated SR exchange_boundaries_nnow
!  Adaptations to work with the DFI (by Lucio Torrisi)
! V4_1         2007/12/04 Ulrich Schaettler
!  Compute MAXVAL for zwmax only for istartpar:iendpar,jstartpar:jendpar
!  to get reproducible results
! V4_4         2008/07/16 Ulrich Schaettler, Michael Baldauf
!  Adapted interface of get_timings
!  Adaptations to changes in src_advection_rk
!  Replaced lkainfri by itype_conv
!  Introduced tendencies of sub-grid scale orography scheme (J-P Schulz)
! V4_5         2008/09/10 Guenther Zaengl, Dmitrii Mironov, Michael Baldauf
!  Adaptations for new reference atmosphere (Guenther Zaengl)
!  Deliver tendencies of qi and qc of the new convection scheme for the
!    cloud microphysics scheme. (D. Mironov, M. Baldauf)
! V4_8         2009/02/16 Michael Baldauf
!  Stable integration of coriolis terms (is still hardcoded to false; 
!  for using it, set l_Coriolis_every_RK_substep=.TRUE.)
!  Bug fix in the computation of dt0dz for irefatm=2
!    (hhl has to be used on main levels)
! V4_9         2009/07/16 Guenther Zaengl
!  Option for using potential temperature as prognostic variable in RK part
!  Check calls to get_timings with ltime (U. Schaettler)
!  Inserted compiler directives
!  Inserted call to SR collapse to change loop indices
! V4_11        2009/11/30 Guenther Zaengl
!  Implemented optional use of irefatm=3 with constant BV frequency
! V4_12        2010/05/11 Oli Fuhrer
!  Replaced DOUBLE PRECISION by KIND-definition in SR coriolis
!  Adapted call to SR comp_hori_diffusion and eliminated option itype_hdiff=3
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_15        2010/11/19 Michael Baldauf
!  Possibility to call a 3rd order implicit vertical advection.
!  New subroutine 'calc_wcon_sqrtg'
! V4_17        2011/02/24 Ulrich Blahak
!  Adapted interface of exchg_boundaries; corrected kzdims(1:20) -> kzdims(1:24);
!  eliminated my_peri_neigh; added possibility to output max. V_h and min./max. W
!  to stdout for diagnostic purposes, if (ldebug_dyn .AND. idbg_level > 3).
! V4_18        2011/05/26 Ulrich Schaettler, Michael Baldauf
!  Introduced conditional compilation for Nudging
!  Bug fix: Subr. 'coriolis' is now called with the correct timelevels
!     in the case l_Coriolis_every_RK_substep=.TRUE. (Michael Baldauf)
!  Implementation of COSMO-ART interfaces (Christoph Knote)
! V4_20        2011/08/31 Michael Baldauf
!  Modified number of small time steps for the Runge-Kutta substeps
!  Call SR coriolis only, if lcori is true
!  Initialize some variables in any case, not only for itheta_avd >= 1 (Oli Fuhrer)
!  Print debug output only in case of ldebug_dyn
! V4_23        2012/05/10 Michael Baldauf, Ulrich Schaettler, Oliver Fuhrer
!  Shifted SR calc_wcon_sqrtg from src_runge_kutta to here
!  Allow new variants of the Bott-Advection schemes:
!    BOTT2_STRANG_B, BOTT4_STRANG_B: Strang Splitting only at the bottom (5 levels)
!    BOTT2_XYZYX, BOTT4_XYZYX: modified sequence 'xyzyx' compared to the current
!  Removed src_2timelevel and related stuff
!  Call to new SR compute_moisture_divergence from src_slow_tendencies_rk
!   for dqvdt (Oliver Fuhrer)
!  Eliminated field qvt_diff (Oliver Fuhrer)
!  Removed computations of total physical tendencies and moved them to
!    organize_physics (Oli Fuhrer)
! V4_24        2012/06/22 Michael Baldauf
!  Optional call of the new subroutine fast_waves_strong_conserv
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer, Ulrich Blahak
!                         Hans-Juergen Panitz
!  Replaced qx-variables by using them from the tracer module
!  UB:
!  Implemented internal switch "l2mom_satads" to be able
!   to switch on the extra saturation adjustments outside the microphysics parts.
!  Introduced nexch_tag for MPI boundary exchange tag to replace ntstep (HJP)
! V4_26        2012/12/06 Anne Roches
!  Replacement of hd_corr_q_XXX by hd_corr_trcr_XXX in order to be
!    consistent also with the naming of other switches (e.g. ltrcr_trilin,
!  Renaming of T_CLP_POSDEF to T_CLP_ON since only on and off are available
!    for the moment.
! V4_27        2013/03/19 Michael Baldauf, Astrid Kerkweg
!  Introduced call to new subroutine init_fast_waves_sc_3
!     (from module fast_waves_sc.f90)
!  MESSy interface introduced: complete tendencies for tracers (AK)
! V4_28        2013/07/12 Michael Baldauf, KIT, Ulrich Schaettler
!  Proper calling of SR 'l_calc_lhs_at_1st_RKstep' in the case of 3D divergence damping.
!     (changes only, if 3D div.-damping is used) (Michael Baldauf)
!  Changes to adapt COSMO-ART to new tracer module: some dependencies to
!  COSMOART and POLLEN deleted, because this is now handled by the tracer module
!  Use subroutines and variables for vertical grid and reference atmospheres
!    from module vgrid_refatm_utils
!  BUG FIXES:
!   - use time level nnow (instead of nnew) for the moisture variables
!     in the water loading contribution of the buoyancy (only for itype_fast_waves=2)
!   - exchange of wadvt in the case ldyn_bbc=.TRUE. necessary
! V4_29        2013/10/04 Davide Cesari, Ulrich Schaettler, Astrid Kerkweg
!  Implemented new SR finalize_runge_kutta to do a proper clean up of the RK method
!  Unification of MESSy interfaces and COSMO Tracer structure
!  For the COSMO-Model only use vcoord and refatm from vgrid_refatm_utils
! V5_1         2014-11-28 Michael Baldauf, Oliver Fuhrer, Lucio Torrisi
!  Introduction of a new variable p0ref_recip (for 1/p0ref)
!  Replaced ireals by wp (working precision) (OF)
!  Stochastic perturbation of physics tendencies added to T-,u-,v- tendencies
! V5_2         2015-05-21 Ulrich Schaettler
!  Reorganized numbering of MPI datatypes ID
! V5_3         2015-10-09 Michael Baldauf, Andrea Arteaga
!                         Ulrich Blahak, Oliver Fuhrer
!  New Namelist switch l_euler_dynamics (MB)
!  more efficient calculation of qcond now at the beginning of the RK-loop (AA)
!  Added storage of LH by satad on ttens_diab for output of pure diabatic tendencies (UB)
!  Moved initialization of dt0dz (done for ntstep=nstart) to init_dynamics (OF)
!  Implemented serialization comments (OF)
! V5_4         2016-03-10 Oliver Fuhrer, Pascal Spoerri, Xavier Lapillonne
!  Added simple clipping for SL3 advection scheme.
!  Hardcoded number of saturation adjustment iterations to 3 for C++ dycore
!  Update serialization statements
! V5_4a        2016-05-10 Pascal Spoerri, Matthias Raschendorfer
!  Integrated a namelist switch (l_satad_dyn_iter) to allow for a dynamic number 
!  of saturation adjustments. (PS)
!  Eliminated allocation/deallocation of za1t, za2t (are no more used 
!   in src_slow_tendencies_rk) (MR)
! V5_4b        2016-07-12 Jochen Foerstner, Pascal Spoerri, Oliver Fuhrer
!  Save the temperature tendencies for the Tiedtke-Bechtold convection scheme 
!  in ttens_conv  (JF)
!  Integrated the new lateral boundary condition module for FW-solvers (OF, PS)
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
  istartpar,    & ! start index for computations in the parallel program
  iendpar,      & ! end index for computations in the parallel program

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

  ! 5. variables for the time discretization and related variables
  ! --------------------------------------------------------------

  dt,           & ! long time-step
  betasw,       & ! beta-variable for treatment of soundwaves    in w
  betagw,       & ! beta-variable for treatment of gravity-waves in w
  beta2sw,      & ! beta-variable for treatment of soundwaves    in p*, T*
  beta2gw,      & ! beta-variable for treatment of gravity-waves in p*, T*

  ! 8. Organizational variables to handle the COSMO humidity tracers
  ! ----------------------------------------------------------------
  idt_qv, idt_qc, idt_qr, idt_qs, idt_qi, idt_qg, idt_qh

! end of data_modelconfig

!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

  pi,           & ! circle constant

  ! 2. physical constants and related variables
  ! -------------------------------------------

  r_d,          & ! gas constant for dry air
  rdv,          & ! r_d / r_v
  o_m_rdv,      & ! 1 - r_d/r_v
  rvd_m_o,      & ! r_v/r_d - 1
  cp_d,         & ! specific heat of dry air at constant pressure
  cpdr,         & ! 1 / cp_d
  gamma,        & ! = cp_d / cv_d
  lh_v,         & ! latent heat of vapourization
  g,            & ! acceleration due to gravity
  r_earth,      & ! mean radius of the earth


  ! 3. constants for parametrizations
  ! ---------------------------------

  p0ref_recip,  & ! reciprocal of ref. pressure for Exner-fct. (Pa^-1)
  b1,           & ! variables for computing the saturation vapour pressure
  b2w,          & ! over water (w) and ice (i)
  b3,           & !               -- " --
  b4w,          & !               -- " --
  b234w,        & ! b2w * (b3 - b4w)
  aks4            ! variable for horizontal diffusion of fourth order

! end of data_constants

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

  ! 1. constant fields for the reference atmosphere                 (unit)
  ! -----------------------------------------------
  p0         ,    & ! reference pressure at main levels             ( Pa  )
  t0         ,    & ! reference temperature                         ( K   )
  hhl        ,    & ! geometical height of half model levels        ( m )

  ! 2. external parameter fields                                    (unit)
  ! ----------------------------
  fc         ,    & ! coriolis-parameter                            ( 1/s )
  fccos      ,    & ! coriolis-parameter mit cosinus                ( 1/s )
  crlat      ,    & ! cosine of transformed latitude

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
  pptens            ! pp-tendency without sound-wave terms          ( Pa/s )

USE data_fields     , ONLY :   &

  ! 6. fields that are computed in the parametrization and dynamics (unit )
  ! ---------------------------------------------------------------
  rho,          & ! density of moist air         (serialize)
  rdcoef,       & ! Rayleigh damping coefficient (serialize)
  dqvdt,        & ! threedimensional moisture convergence(serialize)( 1/s )
  qvsflx,       & ! surface flux of water vapour                    (kg/m2s)
  lhfl_s,       & ! latent heat flux (surface) (serialize)          ( W/m2)
  qrs,          & ! precipitation water (water loading) (serialize) (kg/kg)
  tinc_lh,      & ! temperature increment due to latent heat        (  K  )
  ttens_diab,   & ! temperature tendency due to pure diabatic processes for output (  K  )
  t_g,          & ! weighted surface temperature (serialize)        (  k  )

  wcon     ,    & ! contravariant vertical velocity
  uadvt    ,    & ! advective tendency of u
  vadvt    ,    & ! advective tendency of v
  wadvt    ,    & ! advective tendency of w
  ppadvt   ,    & ! advective tendency of pp
  tadvt    ,    & ! advective tendency of t
  ttens_conv      ! to save the temperature tendencies for Tiedtke-Bechtold scheme

! end of data_fields

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

  ! 1. start and end of the forecast
  ! --------------------------------
  nstart,        & ! first time step of the forecast
  ntstep,        & ! actual time step
  ! indices for permutation of three time levels
  nnow,          & ! corresponds to ntstep
  nnew,          & ! corresponds to ntstep + 1

  ! 3. controlling the physics
  ! --------------------------
  itype_conv,    & ! type of convection parametrization
  itype_gscp,    & ! type of microphys. parametrization
  itype_turb,    & ! type of turbulent diffusion parametrization
  l3dturb,       & ! 3D-turbulence (additional horizontal diffusion)
  l3dturb_metr,  & ! switch on/off additional metric terms for 3D-turbulence
  lprog_tke,     & ! prognostic treatment of TKE (for itype_turb=5/7)
  l_cosmo_art,   & ! if .TRUE., run the COSMO_ART
  l_pollen,      & ! if .TRUE., run the Pollen
  lsppt,         & ! switch stoch. physics tend. perturbation
  l_dzeta_d_needed, & !

  ! 4. controlling the dynamics
  ! ---------------------------
  l_euler_dynamics, & ! on/off Euler-solver; 
                      ! however, Coriolis force and tracer advection may be performed
  lcori,            & ! if =.TRUE.:  with Coriolis force
  lcori_deep,       & ! if =.TRUE.: take cos(phi) coriolis terms into account
  itype_fast_waves, & ! Type of fast waves solver for Runge-Kutta dynamics (1=old, 2=new)
  ldyn_bbc,         & ! dynamic bottom boundary condition
  l_satad_dyn_iter, & ! Dynamic number of iterations for the saturation adjustment
                      ! in the runge kutta dynamical core. 

  ! 6. controlling the upper boundary condition
  ! -------------------------------------------
  lrubc,            & ! with radiative upper boundary condition

  ! 7. additional control variables
  ! -------------------------------
  lw_freeslip         ! if .TRUE.: with free slip lateral boundary condition and
                      ! if .FALSE. specified lateral boundary values for w


USE data_runcontrol , ONLY :   &

  ! 7. additional control variables
  ! -------------------------------
  lcond,         & ! forecast with condensation/evaporation
  loutput_diab,  & ! internal switch for output of temperature increment due to diabatic processes
  ldiabf_lh,     & ! include diabatic forcing due to latent heat in RK-scheme
  ldiabf_satad,  & ! include diabatic forcing due to saturation adjustment
  l2mom_satads,  & ! in case of 2-moment scheme, do all the satads
                   ! (like for the 1-moment schemes), not just the
                   ! satad after the microphysics at the end of the timestep.
  lperi_x,       & ! if lartif_data=.TRUE.: periodic boundary conditions in x-dir.
                   ! or with Davies conditions (.FALSE.)
  lperi_y,       & ! if lartif_data=.TRUE.: periodic boundary conditions in y-dir.
                   ! or with Davies conditions (.FALSE.)
  l2dim,         & ! 2 dimensional runs
  ltime,         & ! detailed timings of the program are given
  y_scalar_advect,&! type of scalar advection scheme
                   ! "SL3_SC", "SL3_MF", "SL3_SFD", "Bott2", "Bott4"
                   ! "Bott2_Strang", "Bott4_Strang", "vanLeer", "PPM"
  y_vert_adv_dyn,& ! switch to choose between explicit and several implicit schemes
                   ! for the vertical advection of dynamic variables
  lreproduce,    & ! the results are reproducible in parallel mode
  lhordiff,      & ! running with horizontal diffusion
  nbl_exchg,     & ! number of boundlines to exchange
  irunge_kutta,  & ! =1: use new RK scheme,
                   ! =2: use new TVD-RK scheme
  irk_order,     & ! order of the Runge-Kutta time-integration scheme
  iadv_order,    & ! order of the horizontal advection scheme in dyn. core
  ieva_order,    & ! order of the explicit vertical adv. scheme in dyn. core
  itheta_adv,    & ! =0: use T' (perturbation temperature) for advection
                   ! =1: use theta' (perturbation potential temperature)
                   ! =2: use theta (full potential temperature)
  xkd,           & ! coefficient for divergence damping
  itype_hdiff,   & ! type of horizontal diffusion (=1: 4th order linear),
                   ! =2: 4th order linear monotonic with orographic limit)
  hd_corr_trcr_bd,  & ! correction factor for horizontal diffusion flux of tracers in boundary zone
  hd_corr_trcr_in,  & ! correction factor for horizontal diffusion flux of tracers in domain
  lartif_data,   & 

! 12. controlling verbosity of debug output
! -----------------------------------------
  idbg_level,    & ! to control the verbosity of debug output
  ldebug_dyn,    & ! if .TRUE., debug output for dynamics
  lprintdeb_all    ! .TRUE.:  all tasks print debug output
                   ! .FALSE.: only task 0 prints debug output
! end of data_runcontrol 

!------------------------------------------------------------------------------

USE data_parallel,      ONLY :  &
  num_compute,   & ! number of compute PEs
  nboundlines,   & ! number of boundary lines of the domain for which
                   ! no forecast is computed = overlapping boundary
                   ! lines of the subdomains
  ldatatypes,    & ! if .TRUE.: use MPI-Datatypes for some communications
  ltime_barrier, & ! if .TRUE.: use additional barriers for determining the
                   ! load-imbalance
  ncomm_type,    & ! type of communication
  my_cart_id,    & ! rank of this subdomain in the cartesian communicator
  my_cart_neigh, & ! neighbors of this subdomain in the cartesian grid
  icomm_cart,    & ! communicator for the virtual cartesian topology
                   ! that can be used by MPI_WAIT to identify the send
  imp_reals,     & ! determines the correct REAL type used in the model
                   ! for MPI
  imp_integers,  &
  nexch_tag,     & ! tag to be used for MPI boundary exchange
                   !  (in calls to exchg_boundaries)
  sendbuf,       & ! sending buffer for boundary exchange:
                   ! 1-4 are used for sending, 5-8 are used for receiving
  isendbuflen      ! length of one column of sendbuf

!------------------------------------------------------------------------------

! coefficients used in complete_tendencies
! ----------------------------------------

USE src_slow_tendencies_rk, ONLY:  &
  zsqrtgrho_r_s,        & ! reciprocal square root of G * rho at skalar points
  zsqrtgrho_r_u,        & ! reciprocal square root of G * rho at u points
  zsqrtgrho_r_v,        & ! reciprocal square root of G * rho at v points
  zsqrtgrho_r_w,        & ! reciprocal square root of G * rho at w points
  zpia,                 & !
  zpianf,               & !
  ztheta,               & !
  ztheta_l,             & !
  ztmkvm,               & !
  ztmkvh,               & !
  ztch,                 & !
  ztcm,                 & !
  zkh,                  & !
  ztmkvw,               & !
  ztmkvtke,             & !
  lvertdiff_w,          & ! .TRUE. for vertical diffusion of w
  lmoist_turb,          & ! .TRUE. if moist turb. param. (Postdam) is used
  lmassf_diffusion        !

USE src_slow_tendencies_rk, ONLY:  &
  compute_moisture_divergence,     &
  complete_tendencies_init,        &
  complete_tendencies_uvwtpp,      &
  complete_tendencies_uvwtpp_eva,  &
  complete_tendencies_tke,         &
  complete_tendencies_trcr,        &
  complete_tend_uvwtpp_CN3Crow,    &
  explicit_horizontal_diffusion,   &
  implicit_vert_diffusion_uvwt,    &
  implicit_vert_diffusion_uvwt_3D, &
  calc_K_rho_at_y,                 &
  deallocate_K_rho_at_y,           &
  set_beta_xxx_for_3D_diffusion,   &
  deallocate_beta_xxx_for_3D_diffusion

! coefficients used in advection
! ------------------------------

USE src_advection_rk, ONLY:  &
  lef_adv_trcr_notpd,   & ! .TRUE. if Euler forward adv. scheme for tracers
                          ! is NOT positive definite
  ltrcr_trilin,         & ! .TRUE. if trilin. interpol. is used for tracers
  ltrcr_conserv_form,   & ! .TRUE. if tracer transport in conservation form
  lcalrho_advprog,      & ! .TRUE. if rho is advected prognostically
  implicit_sedim,       & !
  advection,            & !
  advection_pd,         & !
  calc_wcon_sqrtg         !

!------------------------------------------------------------------------------

USE environment,              ONLY :  &
  exchg_boundaries,       & ! performs the boundary exchange between 
                            ! neighboring processors
  comm_barrier,           & ! 
  model_abort, collapse

!------------------------------------------------------------------------------

USE time_utilities,           ONLY :  get_timings, i_dyn_computations,      &
     i_slow_tendencies, i_horizontal_diffusion, i_horizontal_advection,     &
     i_fast_waves, i_fast_waves_barrier, i_fast_waves_comm,                 &
     i_add_tend_moisture, i_communications_dyn, i_barrier_waiting_dyn,      &
     i_global_communi_dyn, i_barrier_globcom_dyn


USE meteo_utilities,          ONLY :  &
  satad                     !

USE parallel_utilities,       ONLY :  &
  i_global, j_global,     & !
  global_values             !

USE vgrid_refatm_utils,       ONLY :  &
  refatm, vcoord

USE hori_diffusion,           ONLY :  &
  comp_hori_diff

USE numeric_utilities_rk,     ONLY :  &
  clipping,               & !
  init_bott_coeffs          !

USE grid_metrics_utilities,   ONLY :  &
  calc_dzeta_dxx_at_y,  &
  deallocate_dzeta_dxx_at_y

USE fast_waves_rk,            ONLY :  &
  fast_waves_runge_kutta    ! fast waves solver for Runge-Kutta core

USE fast_waves_sc,            ONLY :  &
  alloc_fast_waves_sc,           &
  init_fast_waves_sc_2,          &
  init_fast_waves_sc_3,          &
  dealloc_fast_waves_sc,         &
  finalize_fast_waves_sc,        &
  calc_total_T_p_rho,            &
  init_div_damping_coeff,        &
  alpha_div_h, alpha_div_v_to_h, &
  lhs_of_tridiag_system_for_w,   &
  l_3D_div_damping,              &
  l_calc_lhs_at_1st_rkstep,      &
  fast_waves_strong_conserv        ! fast waves solver for Runge-Kutta core
                                   ! in strong conservation form

!------------------------------------------------------------------------------

USE src_stoch_physics,        ONLY : pertstoph

!------------------------------------------------------------------------------
  
USE src_lbc,                  ONLY : lbc_compute_tendency

!------------------------------------------------------------------------------

USE data_tracer ,             ONLY :  &
  T_DIFF_ID     , T_DIFF_ON   , T_CLP_ID   , T_CLP_ON, T_ERR_NOTFOUND

!------------------------------------------------------------------------------

USE src_tracer ,              ONLY :  &
  trcr_get, trcr_meta_get, trcr_errorstr, trcr_get_ntrcr, trcr_get_block

!------------------------------------------------------------------------------

#ifdef NUDGING
USE data_lheat_nudge,    ONLY :  &
    llhn,         & ! main switch for latent heat nudging (lhn)
    llhnverif,    & ! main switch for latent heat nudging (lhn)
    tt_lheat                ! profile of t-increments due to latent heating   ( K/s )
                            ! (stored for current and previous timestep)

!------------------------------------------------------------------------------

USE src_lheating,             ONLY :  &
    get_gs_lheating         ! storage of grid scale latent heating for lhn
#endif

!------------------------------------------------------------------------------

! CK 20101204 ART RK interfaces
#ifdef POLLEN
USE data_pollen,         ONLY :     &
    isp_pollen   , & ! number of pollen species
    trcr_idx_pollen, &
    vsed_p           ! sedimentation velocity aerosols                  (m/s)
#endif

#ifdef COSMOART
USE data_cosmo_art,      ONLY :     &
    laero        , & ! of aerosols
    isp_aerotrans, & ! number of transported aerosol variables
    isp_aero     , & ! number of transported aerosol variables
    trcr_idx_aero, &
    vd           , & ! deposition velocity gas phase                    (m/s)
    vdepa        , & ! deposition velocity aerosols                     (m/s)
    vseda            ! sedimentation velocity aerosols                  (m/s)
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================

!------------------------------------------------------------------------------
! Declarations
!------------------------------------------------------------------------------

REAL (KIND = wp)                  :: &
  alpha_rk,       & ! 1. parameter for RK-scheme
  beta_rk,        & ! 2. parameter for RK-scheme
  gamma_rk,       & ! 3. parameter for RK-scheme
  dts,            & ! small time step for time splitting
  dtsmall,        & ! small time step in last step
  dtsmall_1,      & ! small time step first intermediate step
  dtsmall_2,      & ! small time step first intermediate step
  cfl_eva,        & ! CFL-value (for explicit vertical advection)
  epsray            ! Rayleigh damping coefficient for vhmx_vol > vhmx_cfl

INTEGER (KIND=iintegers)          :: &
  ismtstep,       & ! number of small time steps in last step
  ismtstep_1,     & ! number of small time steps in first intermediate step
  ismtstep_2        ! number of small time steps in second intermediate step

INTEGER (KIND=iintegers)          :: &
  j2dim,          & ! middle j index for 2D-runs
  ny_2dim,        & ! number of gridpoints in y-direction for 2D-runs
  im, ip,         & ! index-boundaries for advection stencil
  imipmx            ! maximum of the stencil for explicit vertical advection

REAL (KIND=wp), ALLOCATABLE :: &
   t_bd_tend (:,:,:), &
  pp_bd_tend (:,:,:), &
   u_bd_tend (:,:,:), &
   v_bd_tend (:,:,:), &
   w_bd_tend (:,:,:)

!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure in "src_runge_kutta" for organizing the time stepping  
!------------------------------------------------------------------------------

SUBROUTINE org_runge_kutta

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "src_runge_kutta" is the driving routine 
!   of the module, i.e. it acts as interface to the organizing routine
!   organize_dynamics
!
! Method:
!
!------------------------------------------------------------------------------

! Local scalars:
! -------------
INTEGER (KIND=iintegers) ::  &
  i,j,k,irk,           & !  loop indices
#ifdef COSMOART
  isp,                 & !  loop indices
#endif
  kzdims(24),          & !  vertical dimensions for exchg_boundaries
  iztrcr,              & !
  kitpro,              & !  number of iterations in the saturation adjustment
  izstata,             & !  error status at allocation
  izstatd,             & !  error status at deallocation
  izerror,             & !
  izs,                 & ! sum of error status
  ismts,               & ! actual number of small time steps
  nadv,                & ! actual time level to use in advection
  izdebug

REAL    (KIND=wp   )     ::  &
  zdtr

CHARACTER (LEN=255)      ::  &
  yzerrmsg

CHARACTER (LEN=25)       :: yzroutine

LOGICAL :: lapply_Rayleigh_damping
LOGICAL :: l_Coriolis_every_RK_substep

! For potential temperature advection
REAL    (KIND=wp   )     ::  &
  rovcp, govcp, exner(ie,je,ke)

! Local allocatable arrays:
! -------------------------
REAL (KIND = wp),     ALLOCATABLE :: &  
  u_half(:,:,:)  ,& ! 1. velocity component after 
                    ! half of the small time steps
  v_half(:,:,:)  ,& ! 2. velocity component ...
  w_half(:,:,:)  ,& ! 3. velocity component ...
  epsray(:,:,:)     ! Rayleigh-damping-coefficient

#ifdef COSMOART
REAL (KIND = wp),     ALLOCATABLE :: &   !  COSMO_ART
  zcaero(:,:,:,:)
REAL (KIND=wp),     POINTER ::                        &
  caero_new   (:,:,:,:)  => NULL()                     ! caero at nnew
#endif

#ifdef POLLEN
REAL (KIND = wp),     ALLOCATABLE :: &   !  COSMO_ART
  zcpollen(:,:,:,:)
REAL (KIND=wp),     POINTER ::                        &
  cpollen_new   (:,:,:,:)  => NULL()                     ! caero at nnew
#endif

REAL (KIND = wp),     ALLOCATABLE :: &
  q_cond(:,:,:)

! Local (automatic) arrays:
! -------------------------
REAL    (KIND=wp   )     ::  &
  zwmax(ke)            ! maximum vertical velocity in a horizontal layer

INTEGER (KIND=iintegers) ::  &
  izdiff(trcr_get_ntrcr()) ,& !
  izclip(trcr_get_ntrcr())    !

! Tracer pointers:
!-----------------

REAL (KIND = wp),     POINTER     :: &
  ztrcr(:,:,:)=> NULL()  ,& ! tracer variable at nnew
  qv   (:,:,:)=> NULL()  ,& ! QV at nnew
  qc   (:,:,:)=> NULL()  ,& ! QC at nnew
  qi   (:,:,:)=> NULL()  ,& ! QI at nnew
  qr   (:,:,:)=> NULL()  ,& ! QR at nnew
  qs   (:,:,:)=> NULL()  ,& ! QS at nnew
  qg   (:,:,:)=> NULL()  ,& ! QG at nnew
  qh   (:,:,:)=> NULL()     ! QH at nnew

REAL (KIND = wp),     POINTER     :: &
  qv_nnow   (:,:,:)=> NULL()  ,& ! QV at nnow
  qc_nnow   (:,:,:)=> NULL()  ,& ! QC at nnow
  qi_nnow   (:,:,:)=> NULL()  ,& ! QI at nnow
  qr_nnow   (:,:,:)=> NULL()  ,& ! QR at nnow
  qs_nnow   (:,:,:)=> NULL()  ,& ! QS at nnow
  qg_nnow   (:,:,:)=> NULL()  ,& ! QG at nnow
  qh_nnow   (:,:,:)=> NULL()     ! QH at nnow

! For debug output of MAX(U) and MIN/MAX(W) every timestep to stdout:
REAL    (KIND=wp   )     ::  &
     vhmx_loc, wmx_loc, wmn_loc

! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine org_runge_kutta
!------------------------------------------------------------------------------

  yzroutine = 'org_runge_kutta'

! Initialize, whether debug output shall be done
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

! l_Coriolis_every_RK_substep = .TRUE.      ! stable variant
  l_Coriolis_every_RK_substep = .FALSE.

  IF ( .NOT. l_euler_dynamics ) THEN
    l_Coriolis_every_RK_substep = .FALSE.
    ! (remark: Coriolis force is possible, but cannot be treated in the RK-scheme)
  END IF

  SELECT CASE( TRIM(y_scalar_advect) )
  CASE( "PPM", "PPM_STRANG" )
    lef_adv_trcr_notpd = .TRUE.
  CASE DEFAULT
    lef_adv_trcr_notpd = .FALSE.
  END SELECT

  ! Tri-linear or tri-cubic interpolation for tracers
  ltrcr_trilin = .FALSE.

  ! Transport of tracers in conservation form on or off
  ltrcr_conserv_form = .TRUE.

  ! Compute advection of rho prognostically simultaneously to qx transport
  lcalrho_advprog = .TRUE.

  ! Vertical diffusion of vertical velocity on or off
  lvertdiff_w = .TRUE.

  ! Vertical mass redistribution scheme for cloud water on or off
  lmassf_diffusion = .TRUE. ! after calculation of the diffusion

  ! Include diabatic forcing due to latent heat of saturation adjustment
  ldiabf_satad = .FALSE.


  ! Prog. advection of rho is only an option if conservation form for tracers
  ! is chosen.
  IF ( .NOT.ltrcr_conserv_form ) lcalrho_advprog = .FALSE.

  ! Set ldiabf_satad according to NAMELIST parameter ldiabf_lh
  IF ( .NOT.ldiabf_lh ) ldiabf_satad = .FALSE.

  ! Set constants for potential temperature advection
  rovcp = r_d/cp_d
  govcp = g/cp_d

  ! Retrieve the required metadata
  CALL trcr_meta_get(izerror, T_DIFF_ID, izdiff)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_meta_get(izerror, T_CLP_ID, izclip)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

!------------------------------------------------------------------------------
! Section 1: Some preparations at the initial time step(s):
!------------------------------------------------------------------------------

  IF ( ntstep == nstart ) THEN

    IF (izdebug > 5) THEN
      PRINT *, '         Initializations for nstart'
    ENDIF

    SELECT CASE( TRIM(y_scalar_advect) )
    CASE( "BOTT2", "BOTT2_STRANG", "BOTT2_XYZYX", "BOTT2_STRANG_B",   &
          "BOTT4", "BOTT4_STRANG", "BOTT4_XYZYX", "BOTT4_STRANG_B" )
      CALL init_bott_coeffs
    END SELECT

    CALL calc_small_timestep( dt, ismtstep, ismtstep_1, ismtstep_2,         &
                              dtsmall, dtsmall_1, dtsmall_2 )

    IF ( itype_fast_waves == 2 ) CALL init_fast_waves_sc_3( dtsmall )

    ! some preparations to simplify the choice of
    ! different advection shemes 
    ! setting of the indices of the advected stencil
    SELECT CASE(iadv_order)
    CASE(1)
      im = -1
      ip = 1
    CASE(2)
      im = -1
      ip = 1
    CASE(3)
      im = -2
      ip = 2
    CASE(4)
      im = -2
      ip = 2
    CASE(5)
      im = -3
      ip = 3
    CASE(6)
      im = -3
      ip = 3
    END SELECT
    IF ( y_vert_adv_dyn == "expl" ) THEN
      cfl_eva = 1.0_wp
      SELECT CASE(ieva_order)
      CASE(1)
        imipmx = 1
      CASE(2)
        imipmx = 1
      CASE(3)
        imipmx = 2
        IF ( irk_order == 3 ) cfl_eva = 1.4_wp
      CASE(4)
        imipmx = 2
        cfl_eva = 0.9_wp
      CASE(5)
        imipmx = 3
        cfl_eva = 0.95_wp
      CASE(6)
        imipmx = 3
        cfl_eva = 0.75_wp
      END SELECT
    END IF

    j2dim = nboundlines + 1
    ny_2dim = nboundlines + j2dim

    IF (my_cart_id == 0) THEN
      PRINT *,' LEVELINDEX WHERE LEVELS BECOME FLAT, KFLAT = ', vcoord%kflat
      PRINT *,' BETA_SW  = ', betasw,  ' BETA_GW  = ', betagw
      PRINT *,' BETA2_SW = ', beta2sw, ' BETA2_GW = ', beta2gw
      PRINT *,' XKD = ', xkd
      IF (lrubc) THEN
        PRINT *,' CAUTION!!: RADIATIVE UPPER BOUNDARY CONDITION NOT &
          &IMPLEMENTED'
        PRINT *,' CAUTION!!: MODEL RUN USES RIGID UPPER LID'
      ENDIF
    ENDIF

    ! Note: the calculation of dt0dz which was here before has been
    !       moved to the init_dynamcis in organize_dynamics.f90

  ENDIF

  IF (ldebug_dyn .AND. idbg_level > 3) THEN
    ! Determine the min/.max. W in the model domain and the max |vh| and write to stdout
    ! for diagnostic purposes:
    ! all processes have to go to this loop because of global communication
    vhmx_loc = SQRT(MAXVAL( u(:,:,:,nnow)**2 + v(:,:,:,nnow)**2 ))
    wmx_loc  = MAXVAL(w(:,:,:,nnow))
    wmn_loc  = MINVAL(w(:,:,:,nnow))

    IF (ltime) THEN
      CALL get_timings (i_dyn_computations, ntstep, dt, izerror)
      IF (ltime_barrier) THEN
        CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
        CALL get_timings (i_barrier_globcom_dyn, ntstep, dt, izerror)
      ENDIF
    ENDIF

    IF (num_compute > 1) THEN
      CALL global_values (vhmx_loc, 1, 'MAX', imp_reals, icomm_cart, -1,    &
           yzerrmsg, izerror)
      CALL global_values (wmx_loc, 1, 'MAX', imp_reals, icomm_cart, -1,    &
           yzerrmsg, izerror)
      CALL global_values (wmn_loc, 1, 'MIN', imp_reals, icomm_cart, -1,    &
           yzerrmsg, izerror)
    ENDIF

    IF (ltime) CALL get_timings (i_global_communi_dyn, ntstep, dt, izerror)

    IF ( my_cart_id == 0 ) THEN
      WRITE (*,'(a,es12.5)') 'Maximum V_h : ', vhmx_loc
      WRITE (*,'(a,es12.5,a,es12.5)') 'Maximum W : ', wmx_loc, '     Minimum W : ', wmn_loc
    END IF
  END IF

  ALLOCATE( epsray(1:ie, 1:je, 1:ke), STAT=izstata )
  IF ( izstata /= 0 ) THEN
    yzerrmsg="allocation of epsray"
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, yzroutine)
  END IF

  ALLOCATE( q_cond(1:ie, 1:je, 1:ke), STAT=izstata )
  IF ( izstata /= 0 ) THEN
    yzerrmsg="allocation of q_cond"
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, yzroutine)
  END IF

  IF ( l_euler_dynamics ) THEN
    CALL check_cfl_horiz_advection ( lapply_Rayleigh_damping, epsray )
  ELSE
    lapply_Rayleigh_damping = .FALSE.
  END IF

  ! get required tracers (all qx at nnew)
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnew, ptr = qv)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nnew, ptr = qc)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qi, ptr_tlev = nnew, ptr = qi)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qr, ptr_tlev = nnew, ptr = qr)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qs, ptr_tlev = nnew, ptr = qs)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qg, ptr_tlev = nnew, ptr = qg)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  ! try to get any possible moisture variable (at time level nnow):
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnow, ptr = qv_nnow)
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nnow, ptr = qc_nnow)
  CALL trcr_get(izerror, idt_qi, ptr_tlev = nnow, ptr = qi_nnow)
  CALL trcr_get(izerror, idt_qr, ptr_tlev = nnow, ptr = qr_nnow)
  CALL trcr_get(izerror, idt_qs, ptr_tlev = nnow, ptr = qs_nnow)
  CALL trcr_get(izerror, idt_qg, ptr_tlev = nnow, ptr = qg_nnow)
  CALL trcr_get(izerror, idt_qh, ptr_tlev = nnow, ptr = qh_nnow)

#ifdef TWOMOM_SB
  CALL trcr_get(izerror, idt_qh, ptr_tlev = nnew, ptr = qh)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
#endif

!------------------------------------------------------------------------------
! Section 2: Allocate and store boundary tendencies for u and v in case of
!            non-periodic boundary conditions.
!------------------------------------------------------------------------------

  IF (izdebug > 5) THEN
    PRINT *, '         Store boundary tendencies for u and v'
  ENDIF

  ! These tendencies are required to compute the wind divergence along the
  ! boundary lines in the fast-wave solver.
  ! The boundary tendencies are premultiplied by the small time step
  ! The nnew-values are the boundary values provided in the beginning of the
  ! time step.

  ALLOCATE ( t_bd_tend(ie,je,ke), pp_bd_tend(ie,je,ke), &
             u_bd_tend(ie,je,ke), v_bd_tend(ie,je,ke),  &
             w_bd_tend(ie,je,ke1), STAT=izstata )
  IF ( izstata /= 0 ) THEN
    yzerrmsg="allocation of t_bd_tend, ..."
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, yzroutine)
  END IF

  ! NOTE: computation of lateral boundary condition tendency using src_lbc
  !       (this is used in the fast-waves solver to update boundary values on
  !       each small timestep)
  CALL lbc_compute_tendency( t_bd_tend, t(:,:,:,nnew), t(:,:,:,nnow), dt,        &
                             ie, je, ke, istart, iend, jstart, jend, 1, ke )
  CALL lbc_compute_tendency( pp_bd_tend, pp(:,:,:,nnew), pp(:,:,:,nnow), dt,     &
                             ie, je, ke, istart, iend, jstart, jend, 1, ke )
  CALL lbc_compute_tendency( u_bd_tend, u(:,:,:,nnew), u(:,:,:,nnow), dt,        &
                             ie, je, ke, istartu, iendu, jstartu, jendu, 1, ke )
  CALL lbc_compute_tendency( v_bd_tend, v(:,:,:,nnew), v(:,:,:,nnow), dt,        &
                             ie, je, ke, istartv, iendv, jstartv, jendv, 1, ke )
  IF ( .NOT. lw_freeslip ) THEN
    CALL lbc_compute_tendency( w_bd_tend, w(:,:,:,nnew), w(:,:,:,nnow), dt,      &
                               ie, je, ke1, istart, iend, jstart, jend, 1, ke1 )
  ENDIF

  ! Allocation of the separate tendency arrays for advection
  izs = 0
  ALLOCATE ( uadvt (ie,je,ke),  STAT=izstata ); izs = izs + izstata
  ALLOCATE ( vadvt (ie,je,ke),  STAT=izstata ); izs = izs + izstata
  ALLOCATE ( wadvt (ie,je,ke1), STAT=izstata ); izs = izs + izstata
  ALLOCATE ( ppadvt(ie,je,ke),  STAT=izstata ); izs = izs + izstata
  ALLOCATE ( tadvt (ie,je,ke),  STAT=izstata ); izs = izs + izstata

  IF ( izs /= 0 ) THEN
    yzerrmsg="allocation of xadvt"
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, yzroutine)
  END IF

  ! NOTE: tendencies are set to zero in the boundary zone
  uadvt (:,:,:) = 0.0_wp
  vadvt (:,:,:) = 0.0_wp
  wadvt (:,:,:) = 0.0_wp
  ppadvt(:,:,:) = 0.0_wp
  tadvt (:,:,:) = 0.0_wp

  ALLOCATE ( wcon (ie,je,ke1), STAT=izstata )
  ! NOTE: wcon is set to zero in the boundary zone
  wcon(:,:,:) = 0.0_wp

  IF (ltime) CALL get_timings (i_dyn_computations, ntstep, dt, izerror)

!------------------------------------------------------------------------------
! Section 3: Do forecast using the two time level scheme
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 3a: Complete the slow tendencies from physical process by some
!             adiabatic forcings (Coriolis force, water loading etc.)
!             and add them to the advective tendencies.
!------------------------------------------------------------------------------

  IF (izdebug > 5) THEN
    PRINT *, '         Complete slow tendencies'
  ENDIF

  ! Add Rayleigh friction to the velocity field.
  ! ---------------------------------------------------------------------------

  ! Set the time level for Rayleigh friction and add damping terms

  IF ( lapply_Rayleigh_damping ) THEN
    DO  k = 1 , ke
      DO   j = jstartu, jendu
        DO i = istartu, iendu
          utens(i,j,k) = utens(i,j,k) - epsray(i,j,k) * u(i,j,k,nnow)
        ENDDO
      ENDDO
      DO   j = jstartv, jendv
        DO i = istartv, iendv
          vtens(i,j,k) = vtens(i,j,k) - epsray(i,j,k) * v(i,j,k,nnow)
        ENDDO
      ENDDO
    ENDDO
  END IF

  DEALLOCATE( epsray, STAT=izstatd )
  IF ( izstatd /= 0 ) THEN
    yzerrmsg="deallocation of epsray"
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, yzroutine)
  END IF

  ! add estimated temperature tendency due to latent heat
  ! release in previous time step
  IF ( ldiabf_lh ) THEN
    zdtr = 1.0_wp / dt
    DO k = 1, ke
      DO j = jstart, jend 
        DO i = istart, iend
          ttens(i,j,k) = ttens(i,j,k) + zdtr * tinc_lh(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! precalculate some coefficients for complete_tendencies_...
  ! ----------------------------------------------------------------

  ALLOCATE ( zsqrtgrho_r_s(ie,je,ke), zsqrtgrho_r_u(ie,je,ke),           &
    zsqrtgrho_r_v(ie,je,ke), STAT=izstata )
  IF ( izstata /= 0 ) THEN
    yzerrmsg="allocation of zsqrtgrho_r_s, ..."
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, yzroutine)
  END IF

  ALLOCATE ( zpia(ie,je,ke), zpianf(ie,je,ke1),                          &
    ztmkvm(ie,je,ke), ztmkvh(ie,je), ztch(ie,je), ztcm(ie,je),           &
    ztheta(ie,je,ke), zkh(ie,je,ke), STAT=izstata )
  IF ( izstata /= 0 ) THEN
    yzerrmsg="allocation of zpia, zpianf"
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, yzroutine)
  END IF

  IF ( lvertdiff_w ) THEN
    ALLOCATE ( ztmkvw(ie,je,ke),                        &
      zsqrtgrho_r_w(ie,je,ke), STAT=izstata )
    IF ( izstata /= 0 ) THEN
      yzerrmsg="allocation of ztmkvw, ..."
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, yzroutine)
    END IF
  END IF

  IF ( lprog_tke ) THEN
    ALLOCATE ( ztmkvtke(ie,je,ke), STAT=izstata )
    IF ( izstata /= 0 ) THEN
      yzerrmsg="allocation of ztmkvtke"
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, yzroutine)
    END IF
  END IF

  IF ( itype_turb == 6 .OR. itype_turb == 8 ) THEN
    lmoist_turb = .TRUE.
    IF ( l3dturb ) THEN
      ALLOCATE ( ztheta_l(ie,je,ke), STAT=izstata )
      IF ( izstata /= 0 ) THEN
        yzerrmsg="allocation of ztheta_l"
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, yzroutine)
      END IF
    END IF
  ELSE
    lmoist_turb = .FALSE.
  END IF

  CALL complete_tendencies_init

  IF (ltime) CALL get_timings (i_slow_tendencies, ntstep, dt, izerror)

  ! initialize temperature tendencies field for IFS convection scheme
  IF ( itype_conv == 2 ) THEN
    ttens_conv(:,:,:) = ttens_conv(:,:,:) - ttens(:,:,:)
  END IF

  ! NOTE: how many values in boundary zone required for
  !       explicit_horizontal_diffusion???

  ! Horizontal diffusion
  ! ----------------------------------------------------------------

  IF ( l3dturb ) THEN
    IF (  l_dzeta_d_needed ) THEN
      CALL calc_dzeta_dxx_at_y( izerror )
    END IF
    CALL calc_K_rho_at_y
    CALL set_beta_xxx_for_3D_diffusion

    CALL explicit_horizontal_diffusion
  END IF

  ! Vertical diffusion
  ! ----------------------------------------------------------------

  ! NOTE: no boundary values required as all computations are in columns for
  !       implicit_vert_diffusion_uvwt

  IF ( l3dturb .AND. l3dturb_metr ) THEN
    CALL implicit_vert_diffusion_uvwt_3D
  ELSE
    CALL implicit_vert_diffusion_uvwt
  END IF

  ! save temperature tendencies due to diffusion for IFS convection scheme
  ! in case diffusion is calculated here
  IF ( itype_conv == 2 ) THEN
    ttens_conv(:,:,:) = ttens_conv(:,:,:) + ttens(:,:,:)
  END IF

  IF (lsppt) THEN

    !$ser savepoint VerticalDiffusionUnittest.DoSPPT-in LargeTimeStep=ntstep
    !$ser data t_tens=ttens u_tens=utens v_tens=vtens &
    !$ser&     pertstoph=pertstoph tinc_lh=tinc_lh

    DO k = 1, ke
      IF ( ldiabf_lh ) THEN
        DO j = jstart, jend
          DO i = istart, iend
            ttens  (i,j,k) = (ttens(i,j,k) - zdtr * tinc_lh(i,j,k)) * &
                           pertstoph(i,j,k) + zdtr * tinc_lh(i,j,k)
          ENDDO
        ENDDO
      ELSE
        DO j = jstart, jend
          DO i = istart, iend
            ttens  (i,j,k) = ttens(i,j,k) * pertstoph(i,j,k)
          ENDDO
        ENDDO
      ENDIF
! spec. hum. perturbed in src_slow_tendencies_rk
      DO    j = jstartu, jendu
        DO  i = istartu, iendu
          utens (i,j,k) = utens (i,j,k) * 0.5_wp * (pertstoph(i,j,k)+  &
                          pertstoph(i+1,j,k))
        ENDDO
      ENDDO
      DO    j = jstartv, jendv
        DO  i = istartv, iendv
          vtens (i,j,k) = vtens (i,j,k) * 0.5_wp * (pertstoph(i,j,k)+  &
                          pertstoph(i,j+1,k))
        ENDDO
      ENDDO
    ENDDO

    !$ser savepoint VerticalDiffusionUnittest.DoSPPT-out LargeTimeStep=ntstep
    !$ser data t_tens=ttens u_tens=utens v_tens=vtens
  END IF

  ! substract temperature tendencies field for IFS convection scheme
  ! for calculation of temperature tendencies due to dynamics at the end
  IF ( itype_conv == 2 ) THEN
    ttens_conv(:,:,:) = ttens_conv(:,:,:) - ttens(:,:,:)
    IF ( ldiabf_lh ) ttens_conv(:,:,:) = ttens_conv(:,:,:) + zdtr * tinc_lh(:,:,:)
  END IF

  IF (ltime) CALL get_timings (i_horizontal_diffusion, ntstep, dt, izerror)

  !$ser savepoint ConvertTemperatureUnittest.DoTP-in LargeTimeStep=ntstep
  !$ser data t_nnow=t(:,:,:,nnow) t_nnew_bd=t(:,:,:,nnew)

  ! NOTE: Beware, transformation to perturbation form of T is done also in boundary zone!
  !       --> All LBCs on t() should now be be in transformed T!!!

  IF (itheta_adv == 0) THEN
    ! calculate temperature perturbation tp=t-t0
    ! (temporarily override array t)
!CDIR COLLAPSE
    DO k = 1, ke
      DO j = 1, je
        DO i = 1, ie
          t(i,j,k,nnow) = t(i,j,k,nnow) - t0(i,j,k)
!$ser verbatim t(i,j,k,nnew) = t(i,j,k,nnew) - t0(i,j,k) ! not needed
                                                         ! but useful for serialization
        END DO
      END DO
    END DO
  ELSE IF (itheta_adv == 1) THEN
    ! Convert temperature into potential temperature and compute perturbation
!CDIR COLLAPSE
    DO k = 1, ke
      DO j = 1, je
        DO i = 1, ie
          exner(i,j,k) = EXP(-rovcp*LOG((p0(i,j,k)+pp(i,j,k,nnow))*p0ref_recip))
          t(i,j,k,nnow) = t(i,j,k,nnow)*exner(i,j,k) - t0(i,j,k)
          ttens(i,j,k)  = ttens(i,j,k)*exner(i,j,k)
!           t(i,j,k,nnew) = t(i,j,k,nnew)*EXP(-rovcp*LOG((p0(i,j,k)+ &
!                           pp(i,j,k,nnew))*p0ref_recip)) - t0(i,j,k) ! not needed
        END DO
      END DO
    END DO
  ELSE IF (itheta_adv == 2) THEN
    ! Convert temperature into potential temperature
!CDIR COLLAPSE
    DO k = 1, ke
      DO j = 1, je
        DO i = 1, ie
          exner(i,j,k) = EXP(-rovcp*LOG((p0(i,j,k)+pp(i,j,k,nnow))*p0ref_recip))
          t(i,j,k,nnow) = t(i,j,k,nnow)*exner(i,j,k)
          ttens(i,j,k)  = ttens(i,j,k)*exner(i,j,k)
!           t(i,j,k,nnew) = t(i,j,k,nnew)*EXP(-rovcp*LOG((p0(i,j,k)+ &
!                           pp(i,j,k,nnew))*p0ref_recip)) ! not needed
        END DO
      END DO
    END DO
  ENDIF
  IF ( ANY( itheta_adv == (/1,2/) ) ) THEN
    t_bd_tend(:,:,:) = t_bd_tend(:,:,:) * exner(:,:,:)
  ENDIF

  !$ser savepoint ConvertTemperatureUnittest.DoTP-out LargeTimeStep=ntstep
  !$ser data t_nnow=t(:,:,:,nnow) t_nnew_bd=t(:,:,:,nnew)

  IF ( lcori .AND. (.NOT. l_Coriolis_every_RK_substep) ) THEN
    CALL coriolis( utens, vtens, wtens, nnow )
  END IF

  IF ( itype_fast_waves == 2 ) THEN
    ! q_cond should contain ALL condensated moisture fractions.
    ! Use of humidity tracers was with nadv before; but at this point,
    ! nadv is always nnew

    IF (               ASSOCIATED(qr_nnow)     &
      &   .AND. (.NOT. ASSOCIATED(qi_nnow))    &
      &   .AND. (.NOT. ASSOCIATED(qs_nnow))    &
      &   .AND. (.NOT. ASSOCIATED(qg_nnow)) ) THEN
      ! e.g. if the Kessler warm rain scheme is used
      q_cond(:,:,:) = qc_nnow(:,:,:) + qr_nnow(:,:,:)

    ELSE IF (          ASSOCIATED(qi_nnow)     &
      &   .AND.        ASSOCIATED(qr_nnow)     &
      &   .AND.        ASSOCIATED(qs_nnow)     &
      &   .AND. (.NOT. ASSOCIATED(qg_nnow)) ) THEN
      ! e.g. if the 5-class microphysics scheme is used
      q_cond(:,:,:) = qc_nnow(:,:,:) + qi_nnow(:,:,:) + qr_nnow(:,:,:) + qs_nnow(:,:,:)

    ELSE IF (          ASSOCIATED(qi_nnow)     &
      &   .AND.        ASSOCIATED(qr_nnow)     &
      &   .AND.        ASSOCIATED(qs_nnow)     &
      &   .AND.        ASSOCIATED(qg_nnow) ) THEN
      ! e.g. if the 6-class graupel microphysics scheme is used
      q_cond(:,:,:) = qc_nnow(:,:,:) + qi_nnow(:,:,:) + qr_nnow(:,:,:) + qs_nnow(:,:,:) + qg_nnow(:,:,:)

    ELSE
      q_cond(:,:,:) = qc_nnow(:,:,:)
      IF ( ASSOCIATED(qi_nnow) )  q_cond(:,:,:) = q_cond(:,:,:) + qi_nnow(:,:,:)
      IF ( ASSOCIATED(qr_nnow) )  q_cond(:,:,:) = q_cond(:,:,:) + qr_nnow(:,:,:)
      IF ( ASSOCIATED(qs_nnow) )  q_cond(:,:,:) = q_cond(:,:,:) + qs_nnow(:,:,:)
      IF ( ASSOCIATED(qg_nnow) )  q_cond(:,:,:) = q_cond(:,:,:) + qg_nnow(:,:,:)
      IF ( (my_cart_id == 0) .AND. (ntstep==0) ) THEN
        WRITE(*,*) "HINT in subroutine org_runge_kutta:"
        WRITE(*,*) "  this is a generic version to calculate q_cond."
        WRITE(*,*) "  You can improve its efficiency analogous to the programming lines above."
      END IF
    END IF

#ifdef TWOMOM_SB
    IF ( ASSOCIATED( qh_nnow ) )  q_cond(:,:,:) = q_cond(:,:,:) + qh_nnow(:,:,:)
#endif

  END IF

  !----------------------------------------------------------------------------
  ! Runge-Kutta time integration ...
  !----------------------------------------------------------------------------

  ! set time level for advection in first RK step
  nadv = nnow

  alpha_rk = 1.0_wp
  beta_rk  = 0.0_wp

  !$ser savepoint TimeIntegratorUnittest.DoStep-in LargeTimeStep=ntstep
  !$ser data  u_nnow=u(:,:,:,nnow)   u_tens=utens
  !$ser data  v_nnow=v(:,:,:,nnow)   v_tens=vtens
  !$ser data  w_nnow=w(:,:,:,nnow)   w_tens=wtens
  !$ser data  t_nnow=t(:,:,:,nnow)   t_tens=ttens
  !$ser data pp_nnow=pp(:,:,:,nnow) pp_tens=pptens
  !$ser data qrs=qrs rho=rho tinc_lh=tinc_lh rdcoef=rdcoef
  !$ser tracer QV@nnew QC@nnew
  !$ser tracer QV@nnow QC@nnow QI@nnow QR@nnow QS@nnow QG@nnow

  IF ( l_euler_dynamics ) THEN

    runge_kutta_loop: DO irk = 1, irk_order

      IF (izdebug > 5) THEN
        PRINT *, '         Runge-Kutta loop irk = ', irk
      ENDIF

      !--------------------------------------------------------------------------
      ! Section 3b: Compute slow tendencies by first evaluating the advective
      !             and then providing slow tendencies including physical ones.
      !--------------------------------------------------------------------------

      IF (irk_order == 3) THEN
        SELECT CASE(irunge_kutta)
        CASE(1)
          SELECT CASE(irk)
          CASE(1)
            gamma_rk = 1.0_wp/3.0_wp
            ismts = ismtstep_1   ! set number of small time steps
            dts   = dtsmall_1
          CASE(2)
            gamma_rk = 0.5_wp
            ismts = ismtstep_2   ! set number of small time steps
            dts   = dtsmall_2
          CASE(3)
            gamma_rk = 1.0_wp
            ismts = ismtstep     ! set number of small time steps
            dts   = dtsmall
          END SELECT
        CASE(2)
          SELECT CASE(irk)
          CASE(1)
            alpha_rk = 1.0_wp
            beta_rk  = 0.0_wp
            gamma_rk = 1.0_wp
            ismts = ismtstep_1   ! set number of small time steps
            dts   = dtsmall_1
          CASE(2)
            alpha_rk = 0.75_wp
            beta_rk  = 0.25_wp
            gamma_rk = 0.25_wp
            ismts = ismtstep_2   ! set number of small time steps
            dts   = dtsmall_2
          CASE(3)
            alpha_rk = 1.0_wp/3.0_wp
            beta_rk  = 2.0_wp/3.0_wp
            gamma_rk = 2.0_wp/3.0_wp
            ismts = ismtstep     ! set number of small time steps
            dts   = dtsmall
          END SELECT
        END SELECT
      ELSE
        IF (irk == irk_order) THEN
          gamma_rk = 1.0_wp
          ismts = ismtstep     ! set number of small time steps
          dts   = dtsmall
        ELSE
          gamma_rk = 0.5_wp
          ismts = ismtstep_1   ! set number of small time steps
          dts   = dtsmall_1
        END IF
      END IF

      !$ser savepoint TimeIntegratorUnittest.DoRKStage-in LargeTimeStep=ntstep RKStageNumber=irk
      !$ser data u=u(:,:,:,nadv) v=v(:,:,:,nadv) w=w(:,:,:,nadv) pp=pp(:,:,:,nadv) t=t(:,:,:,nadv) &
      !$ser&     u_tens=utens(:,:,:) v_tens=vtens(:,:,:) w_tens=wtens(:,:,:) t_tens=ttens(:,:,:)   &
      !$ser&     pp_tens=pptens(:,:,:) dtadv=gamma_rk*dt dts=dts

      !--------------------------------------------------------------------------
      ! Section 3b1: Compute the time tendencies due to horizontal advection (and Coriolis terms)
      !              and store them on separate tendency fields (uadvt,...)
      !--------------------------------------------------------------------------

      ! Compute advection for dynamic variables
      ! NOTE: advection requires nbl_exchg number of valid boundary gridpoints
      !       for each of the advected fields (pp, t, u, v, w) at timelevel
      !       nadv
      CALL advection( nadv, irk, im, ip, imipmx, cfl_eva, j2dim, ny_2dim)


      IF ( lcori .AND. l_Coriolis_every_RK_substep ) THEN
        ! Compute Coriolisterms for u, v (and perhaps w)
        ! NOTE: coriolis requires 1 line of valid boundary gridpoints
        !       for u and v at timelevel nadv
        CALL coriolis( uadvt, vadvt, wadvt, nadv )
      END IF

      IF (ltime) CALL get_timings (i_horizontal_advection, ntstep, dt, izerror)

      !--------------------------------------------------------------------------
      ! Section 3b2: Compute the slow tendencies
      !              and store them on separate tendency fields (uadvt,...)
      !--------------------------------------------------------------------------

      ! NOTE: no boundary points required as all computations are within the columns

      IF ( y_vert_adv_dyn == "impl2" ) THEN

        IF ( irk_order == 3 .AND. irunge_kutta == 2 ) THEN
          CALL complete_tendencies_uvwtpp( nadv, irk, gamma_rk*dt, alpha_rk, beta_rk )
        ELSE
          CALL complete_tendencies_uvwtpp( nadv, irk, gamma_rk*dt )
        END IF

      ELSE IF ( y_vert_adv_dyn == "expl" ) THEN

        CALL complete_tendencies_uvwtpp_eva

      ELSE IF ( y_vert_adv_dyn == "impl3" ) THEN  

        ! NOTE: maybe reuse & rename complete_tendencies_uvwtpp_eva as this does nothing else than
        !       what the below code does!!!

        ! tendencies are not added in this implicit scheme, therefore do it here explicitely:
        ! NOTE: this computation could be done (istart:iend,jstart:jend) instead of full domain
        !       as utens, vtens, wtens, pptens, ttens should be zero in boundary region
        uadvt(:,:,:)  = uadvt(:,:,:)  + utens(:,:,:)
        vadvt(:,:,:)  = vadvt(:,:,:)  + vtens(:,:,:)
        wadvt(:,:,:)  = wadvt(:,:,:)  + wtens(:,:,:)
        ppadvt(:,:,:) = ppadvt(:,:,:) + pptens(:,:,:)
        tadvt(:,:,:)  = tadvt(:,:,:)  + ttens(:,:,:)

      ELSE
        yzerrmsg="false value for y_vert_adv_dyn"
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, yzroutine)
      END IF

      ! set time level for advection in remaining RK step(s)
      nadv = nnew

      !--------------------------------------------------------------------------
      ! Runge-Kutta step(s) ...
      ! Section 3c: Compute the intermediate values of the dynamic fields
      !             (using the time level nnew for storage) by integrating
      !             the sound and gravity wave terms, including the slow
      !             tendencies with small time steps for the time interval
      !             (dt/4,) (dt/2,) dt with the split-explicit scheme.
      !             Exchange all boundary fields afterwards.
      !--------------------------------------------------------------------------

      ! Exchange one boundary line of uadvt and vadvt:

      IF (ltime) THEN
        CALL get_timings (i_slow_tendencies, ntstep, dt, izerror)
        IF (ltime_barrier) THEN
          CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
          CALL get_timings (i_barrier_waiting_dyn, ntstep, dt, izerror)
        ENDIF
      ENDIF

      IF ( ldyn_bbc ) THEN

        kzdims(1:24) =                  &
          (/  ke,  ke,   2,   0,   0,   &
               0,   0,   0,   0,   0,   &
               0,   0,   0,   0,   0,   &
               0,   0,   0,   0,   0,   &
               0,   0,   0,   0  /)

        CALL exchg_boundaries                                                    &
          ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je, &
            kzdims, jstartpar, jendpar, 2, nboundlines, my_cart_neigh,           &
            lperi_x, lperi_y, l2dim,                                             &
            0, .FALSE., ncomm_type, izerror, yzerrmsg,                           &
            uadvt(:,:,:),          &
            vadvt(:,:,:),          &
            wadvt(:,:,ke:ke1)  )

        ! NOTE: here some BC should probably be applied

      ELSE

        kzdims(1:24) =                     &
          (/  ke,  ke,   0,   0,   0,      &
               0,   0,   0,   0,   0,      &
               0,   0,   0,   0,   0,      & 
               0,   0,   0,   0,   0,      &
               0,   0,   0,   0  /)

        CALL exchg_boundaries                                                    &
          ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je, &
            kzdims, jstartpar, jendpar, 2, nboundlines, my_cart_neigh,           &
            lperi_x, lperi_y, l2dim,                                             &
            0, .FALSE., ncomm_type, izerror, yzerrmsg,                           &
            uadvt(:,:,:),          &
            vadvt(:,:,:)  )

        ! NOTE: here some BC should probably be applied

      END IF


      IF ( izerror /= 0_iintegers ) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, yzroutine)
      END IF

      IF (ltime) CALL get_timings (i_fast_waves_comm, ntstep, dt, izerror)

      IF ( itype_fast_waves == 2 ) THEN

        IF ( irunge_kutta == 1 ) THEN
          u (:,:,:,nnew) = u (:,:,:,nnow)
          v (:,:,:,nnew) = v (:,:,:,nnow)
          w (:,:,:,nnew) = w (:,:,:,nnow)
          T (:,:,:,nnew) = T (:,:,:,nnow)
          pp(:,:,:,nnew) = pp(:,:,:,nnow)
        ELSE
          ! valid for both RK-variants (irunge_kutta=1/2):
          DO k = 1, ke
            DO j = 1, je
              DO i = 1, ie
                u (i,j,k,nnew) = alpha_rk * u (i,j,k,nnow) + beta_rk * u (i,j,k,nnew)
                v (i,j,k,nnew) = alpha_rk * v (i,j,k,nnow) + beta_rk * v (i,j,k,nnew)
                w (i,j,k,nnew) = alpha_rk * w (i,j,k,nnow) + beta_rk * w (i,j,k,nnew)
                pp(i,j,k,nnew) = alpha_rk * pp(i,j,k,nnow) + beta_rk * pp(i,j,k,nnew)
                t (i,j,k,nnew) = alpha_rk * t (i,j,k,nnow) + beta_rk * t (i,j,k,nnew)
              END DO
            END DO
          END DO
          w(:,:,ke1,nnew) = alpha_rk * w(:,:,ke1,nnow) + beta_rk * w(:,:,ke1,nnew)
        END IF

        IF ( l_calc_lhs_at_1st_RKstep .AND. ( irk==1 ) ) THEN
          !$ser savepoint FastWavesSCUnittest.LHS-in LargeTimeStep=ntstep
          !$ser data pp=pp(:,:,:,nnew) t=T(:,:,:,nnew) dts=dts
          !$ser tracer QV@nnow QC@nnow QI@nnow QR@nnow QS@nnow QG@nnow

          ! it is sufficient to calculate the coefficients for
          ! the implicit equation for w only once per timestep.
          ! Unfortunately, some fields must be allocateted already now:
          CALL alloc_fast_waves_sc
          CALL init_fast_waves_sc_2
          CALL calc_total_T_p_rho( T(:,:,:,nnew), pp(:,:,:,nnew), qv_nnow, q_cond )
          CALL init_div_damping_coeff( xkd, dts, alpha_div_h, alpha_div_v_to_h )
          CALL lhs_of_tridiag_system_for_w( dts, l_3D_div_damping )
        END IF

        ! fast_waves_sc updates the values of u,v,w,T,pp at time level 'nnew':
        CALL fast_waves_strong_conserv( ismts, dts, xkd,                              &
          u(:,:,:,nnew), v(:,:,:,nnew), w(:,:,:,nnew), T(:,:,:,nnew), pp(:,:,:,nnew), &
          qv_nnow, q_cond, uadvt, vadvt, wadvt, tadvt, ppadvt,                        &
          t_bd_tend, pp_bd_tend, u_bd_tend, v_bd_tend, w_bd_tend, irk )

      ELSE IF ( itype_fast_waves == 1 ) THEN

        ! premultiply tendencies with small time step
        ! NOTE: this is inefficient, should do this on-the-fly when the fields are
        !       actually being used for the first time (to avoid read/writing 5 full 3D fields)
        uadvt  = uadvt * dts
        vadvt  = vadvt * dts
        wadvt  = wadvt * dts
        tadvt  = tadvt * dts
        ppadvt = ppadvt * dts

        ! NOTE: see LBC pre- and post-requirements in fast_waves_rk.f90
        CALL fast_waves_runge_kutta( ismts, dts, xkd,                          &
          uadvt, vadvt, wadvt, tadvt, ppadvt,                                  &
          t_bd_tend, pp_bd_tend, u_bd_tend, v_bd_tend,                         &
          w_bd_tend, irk, alpha_rk, beta_rk, vcoord%kflat )

      ELSE
        yzerrmsg = "false value in itype_fast_waves"
        CALL model_abort (my_cart_id, 1234, yzerrmsg, yzroutine)
      END IF

      ! NOTE: This should only update the calculation domain of the field,
      !           excluding the boundary
      IF ( irk == irk_order .AND. ldiabf_lh ) THEN
        ! subtract estimated temperature increment due to latent heat 
        ! release in previous time step
        IF (itheta_adv == 0) THEN
          t(istart:iend,jstart:jend,:,nnew) = t(istart:iend,jstart:jend,:,nnew) - tinc_lh(istart:iend,jstart:jend,:)
          !t(:,:,:,nnew) = t(:,:,:,nnew) - tinc_lh(:,:,:)
        ELSE
          t(istart:iend,jstart:jend,:,nnew) = t(istart:iend,jstart:jend,:,nnew) - &
          tinc_lh(istart:iend,jstart:jend,:) * exner(istart:iend,jstart:jend,:)
          !t(:,:,:,nnew) = t(:,:,:,nnew) - tinc_lh(:,:,:) * exner(:,:,:)
        ENDIF
      ENDIF
    
      IF (ltime) THEN
        CALL get_timings (i_fast_waves, ntstep, dt, izerror)
        IF (ltime_barrier) THEN
          CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
          CALL get_timings (i_fast_waves_barrier, ntstep, dt, izerror)
        ENDIF
      ENDIF

      kzdims(1:24)=(/ke,ke,ke1,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                &
        (56+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,&
        ie, je, kzdims, jstartpar, jendpar, nbl_exchg, nboundlines,        &
        my_cart_neigh, lperi_x, lperi_y, l2dim,                            &
        20000+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,        &
        u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),    &
        pp(:,:,:,nnew))

      !$ser savepoint TimeIntegratorUnittest.DoRKStage-out LargeTimeStep=ntstep RKStageNumber=irk
      !$ser data u=u(:,:,:,nnew) v=v(:,:,:,nnew) w=w(:,:,:,nnew) t=t(:,:,:,nnew) pp=pp(:,:,:,nnew)

      IF (ltime) CALL get_timings (i_fast_waves_comm, ntstep, dt, izerror)

    END DO runge_kutta_loop

    !$ser savepoint TimeIntegratorUnittest.DoStep-out LargeTimeStep=ntstep
    !$ser data u_nnow=u(:,:,:,nnow)    u=u(:,:,:,nnew)   u_tens=utens(:,:,:)  
    !$ser data v_nnow=v(:,:,:,nnow)    v=v(:,:,:,nnew)   v_tens=vtens(:,:,:)  
    !$ser data w_nnow=w(:,:,:,nnow)    w=w(:,:,:,nnew)   w_tens=wtens(:,:,:)  
    !$ser data t_nnow=t(:,:,:,nnow)    t=t(:,:,:,nnew)   t_tens=ttens(:,:,:)  
    !$ser data pp_nnow=pp(:,:,:,nnow) pp=pp(:,:,:,nnew) pp_tens=pptens(:,:,:)  

    IF ( ( itype_fast_waves == 2 ) .AND. l_calc_lhs_at_1st_RKstep ) THEN
      CALL dealloc_fast_waves_sc
    END IF

  ELSE  ! i.e. l_euler_dynamics == .FALSE.

    IF ( my_cart_id == 0 ) THEN
      WRITE(*,*) "WARNING: l_euler_dynamics == .FALSE.!"
    END IF

    ! only Euler forward with physics tendencies:
    ! NOTE: this computation could be done (istart:iend,jstart:jend) instead of full domain
    !           as utens, vtens, wtens, pptens, ttens should be zero in boundary region
    u (:,:,:,nnew) = u (:,:,:,nnow) + dt * utens(:,:,:)
    v (:,:,:,nnew) = v (:,:,:,nnow) + dt * vtens(:,:,:)
    w (:,:,:,nnew) = w (:,:,:,nnow) + dt * wtens(:,:,:)
    T (:,:,:,nnew) = T (:,:,:,nnow) + dt * ttens(:,:,:)
    pp(:,:,:,nnew) = pp(:,:,:,nnow) + dt * pptens(:,:,:)

  END IF

  IF ( y_vert_adv_dyn == "impl3" ) THEN
    ! complete operatorsplitting between hor-Adv./RK  and vertical advection

    CALL calc_wcon_sqrtg( u(:,:,:,nnew), v(:,:,:,nnew), w(:,:,:,nnew), wcon(:,:,:) )
    ! attention: now  wcon  contains  'sqrtg * contravar. vertical velocity' !
 
    ! we need 1 boundary line of wcon since the computation of the
    ! courant number for vertical u,v-advection required wcon at i,i+1
    ! and j,j+1 positions.
    IF (ltime) THEN
      CALL get_timings (i_slow_tendencies, ntstep, dt, izerror)
      IF (ltime_barrier) THEN
        CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
        CALL get_timings (i_barrier_waiting_dyn, ntstep, dt, izerror)
      ENDIF
    ENDIF

    kzdims(1:24)=(/ke1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                &
      ( 7     , sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,&
      ie, je, kzdims, jstartpar, jendpar, 1, nboundlines,                &
      my_cart_neigh, lperi_x, lperi_y, l2dim,                            &
      20000+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,        &
      wcon)

    ! NOTE: here some LBC should probably be applied

    IF (ltime) CALL get_timings (i_communications_dyn, ntstep, dt, izerror)

    CALL complete_tend_uvwtpp_CN3Crow( dt, nnew )

    ! NOTE: this computation could be done (istart:iend,jstart:jend) instead of full domain
    !       as uadvt, vadvt, wadvt, ppadvt, tadvt should be zero in boundary region
    u (:,:,:,nnew) = u (:,:,:,nnew) + dt * uadvt (:,:,:)
    v (:,:,:,nnew) = v (:,:,:,nnew) + dt * vadvt (:,:,:)
    w (:,:,:,nnew) = w (:,:,:,nnew) + dt * wadvt (:,:,:)
    pp(:,:,:,nnew) = pp(:,:,:,nnew) + dt * ppadvt(:,:,:)
    t (:,:,:,nnew) = t (:,:,:,nnew) + dt * tadvt (:,:,:)

    ! the updated u,v,w,pp,t need to be updated in the halo
    ! zone in order to be consistent with the other vertical
    ! advection schemes which come out of  the RK-stage-loop
    ! with a final halo-update
    IF (ltime) THEN
      CALL get_timings (i_slow_tendencies, ntstep, dt, izerror)
      IF (ltime_barrier) THEN
        CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
        CALL get_timings (i_barrier_waiting_dyn, ntstep, dt, izerror)
      ENDIF
    ENDIF

    kzdims(1:24)=(/ke,ke,ke1,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                &
      (56+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,&
      ie, je, kzdims, jstartpar, jendpar, nbl_exchg, nboundlines,        &
      my_cart_neigh, lperi_x, lperi_y, l2dim,                            &
      20000+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,        &
      u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),    &
      pp(:,:,:,nnew))

    ! NOTE: here some BC should probably be applied

    IF (ltime) CALL get_timings (i_communications_dyn, ntstep, dt, izerror)

  END IF

  !$ser savepoint ConvertTemperatureUnittest.DoT-in LargeTimeStep=ntstep
  !$ser data t=t(:,:,:,nnew)

  ! NOTE: Beware, back transformation from perturbation form of T is done also in boundary zone!
  !       --> All LBCs on t() should now be in absolute temperatures again

  IF (itheta_adv == 0) THEN
    ! calculate absolute temperature t=tp+t0
!CDIR COLLAPSE
    DO k = 1, ke
      DO j = 1, je
        DO i = 1, ie
          t(i,j,k,nnow) = t(i,j,k,nnow) + t0(i,j,k)
          t(i,j,k,nnew) = t(i,j,k,nnew) + t0(i,j,k)
        END DO
      END DO
    END DO
  ELSE IF (itheta_adv == 1) THEN
    ! convert theta perturbation back to temperature
!CDIR COLLAPSE
    DO k = 1, ke
      DO j = 1, je
        DO i = 1, ie
          t(i,j,k,nnow) = (t(i,j,k,nnow) + t0(i,j,k)) / exner(i,j,k)
          t(i,j,k,nnew) = (t(i,j,k,nnew) + t0(i,j,k)) * &
                          EXP(rovcp*LOG((p0(i,j,k)+pp(i,j,k,nnew))*p0ref_recip))
        END DO
      END DO
    END DO
  ELSE IF (itheta_adv == 2) THEN
    ! convert theta back to temperature
!CDIR COLLAPSE
    DO k = 1, ke
      DO j = 1, je
        DO i = 1, ie
          t(i,j,k,nnow) = t(i,j,k,nnow) / exner(i,j,k)
          t(i,j,k,nnew) = t(i,j,k,nnew) * &
                          EXP(rovcp*LOG((p0(i,j,k)+pp(i,j,k,nnew))*p0ref_recip))
        END DO
      END DO
    END DO
  ENDIF

  !$ser savepoint ConvertTemperatureUnittest.DoT-out LargeTimeStep=ntstep
  !$ser data t=t(:,:,:,nnew)

  DEALLOCATE( q_cond, STAT=izstatd )
  IF ( izstatd /= 0 ) THEN
    yzerrmsg="deallocation of q_cond"
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, yzroutine)
  END IF

  ! Deallocate the private arrays for the advective tendencies
  izs = 0
  DEALLOCATE ( uadvt , STAT=izstatd ); izs = izs + izstatd
  DEALLOCATE ( vadvt , STAT=izstatd ); izs = izs + izstatd
  DEALLOCATE ( wadvt , STAT=izstatd ); izs = izs + izstatd
  DEALLOCATE ( ppadvt, STAT=izstatd ); izs = izs + izstatd
  DEALLOCATE ( tadvt , STAT=izstatd ); izs = izs + izstatd
  IF ( izstatd /= 0 ) THEN
    yzerrmsg="deallocation of uadvt, ..."
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, yzroutine)
  END IF

  ! Deallocate the boundary tendency fields
  DEALLOCATE ( t_bd_tend, pp_bd_tend, u_bd_tend, v_bd_tend,  &
               w_bd_tend, STAT=izstata )
  IF ( izstata /= 0 ) THEN
    yzerrmsg="deallocation of t_bd_tend, ..."
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, yzroutine)
  END IF

  !----------------------------------------------------------------------------
  ! Section 4a: Advection of scalars: tracers and tke
  !----------------------------------------------------------------------------

  IF (izdebug > 5) THEN
    PRINT *, '         Advection of scalars'
  ENDIF

  IF ( y_scalar_advect /= "none" ) THEN

    ! Allocate working arrays
    ALLOCATE ( u_half(ie,je,ke), v_half(ie,je,ke), w_half(ie,je,ke1), &
               STAT=izstata )
    IF ( izstata /= 0 ) THEN
      yzerrmsg="allocation of u_half, ..."
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, yzroutine)
    END IF

    ! Calculate estimates of the velocities for time n+1/2 
    ! -------------------------------------------------------------------------
    ! Integrate the velocities
!CDIR COLLAPSE
    u_half(:,:,:)  = 0.5_wp*( u(:,:,:,nnow)+u(:,:,:,nnew) )
!CDIR COLLAPSE
    v_half(:,:,:)  = 0.5_wp*( v(:,:,:,nnow)+v(:,:,:,nnew) )
!CDIR COLLAPSE
    w_half(:,:,:)  = 0.5_wp*( w(:,:,:,nnow)+w(:,:,:,nnew) )

! Note: the dycore computes half level in the advection part
!       the results just computed are not serializedL the dycore
!       computes the from scratch (i.e. from u_in and u_nnow etc.)

  !$ser savepoint AdvectionPDBottUnittest.DoTracers-in LargeTimeStep=ntstep
  !$ser data u=u(:,:,:,nnew) u_nnow=u(:,:,:,nnow)                                  &
  !$ser&     v=v(:,:,:,nnew) v_nnow=v(:,:,:,nnow)                                  &
  !$ser&     w=w(:,:,:,nnew) w_nnow=w(:,:,:,nnow)                                  &
  !$ser&     rho=rho(:,:,:)
  !$ser tracer %all@nnow

    ! NOTE: Requires variable number of boundary lines depending on advection
    !       scheme used. See also inside advection_pd() for specifics of BC
    !       treatment.

    !
    ! Perform positive-definite advection of moisture variables + TKE
    !
    CALL advection_pd( u_half(:,:,:), v_half(:,:,:), w_half(:,:,:), nnow, dt, &
                       im, ip, j2dim, ny_2dim)

!$ser savepoint AdvectionPDBottUnittest.DoTracers-out LargeTimeStep=ntstep
!$ser data rho=rho
!$ser tracer %all@nnew

    ! Deallocate working arrays
    DEALLOCATE ( u_half, v_half, w_half, STAT=izstatd )
    IF ( izstatd /= 0 ) THEN
      yzerrmsg="deallocation of u_half, ..."
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, yzroutine)
    END IF

  ELSE ! if (y_scalar_advect == "none" )
    IF ( my_cart_id == 0 ) THEN
      PRINT*, "WARNING: tracer-advection is switched off!"
    END IF
  END IF  ! y_scalar_advect /= "none"

  IF (ltime) CALL get_timings (i_add_tend_moisture, ntstep, dt, izerror)

  !$ser savepoint SedimentationUnittest.DoTracers-in LargeTimeStep=ntstep
  !$ser tracer %all@nnew 
  !PS TODO Add %all#sedimvel above

#ifdef COSMOART
  IF (l_cosmo_art) THEN
    IF (laero) THEN

      CALL trcr_get_block(izerror, idx_start=trcr_idx_aero(1), idx_end=trcr_idx_aero(isp_aero), &
                 ptr_tlev = nnew, ptr = caero_new)
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF

      ALLOCATE ( zcaero(ie,je,ke,isp_aerotrans),STAT=izstata )
      DO isp=1,isp_aerotrans

        DO k = 1, ke
          DO j = 1, je
            DO i = 1, ie
              zcaero(i,j,k,isp)=caero_new(i,j,k,isp)
              vseda(i,j,k,isp) = (-1.0_wp)*vseda(i,j,k,isp)
            ENDDO
          ENDDO
        ENDDO
        CALL implicit_sedim( caero_new(:,:,:,isp), zcaero(:,:,:,isp), &
                             vseda(:,:,:,isp), vseda(:,:,:,isp))
      ENDDO
      DEALLOCATE ( zcaero )
    ENDIF
  ENDIF
#endif

#ifdef POLLEN
  ! Sedimentation of pollen
  IF (l_pollen) THEN
    ALLOCATE ( zcpollen(ie,je,ke,isp_pollen),STAT=izstata )

    CALL trcr_get_block(izerror, idx_start=trcr_idx_pollen(1), idx_end=trcr_idx_pollen(isp_pollen), &
               ptr_tlev = nnew, ptr = cpollen_new)
    IF (izerror /= 0) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF

    zcpollen(:,:,:,:)=cpollen_new(:,:,:,:)
    DO isp=1,isp_pollen
      DO k = 1,ke
        DO j = 1,je
          DO i = 1,ie
            vsed_p(i,j,k,isp) = -1.0_wp * vsed_p(i,j,k,isp)
          ENDDO
        ENDDO
      ENDDO
      CALL implicit_sedim( cpollen_new(:,:,:,isp), zcpollen(:,:,:,isp), &
                           vsed_p(:,:,:,isp), vsed_p(:,:,:,isp) )
    ENDDO
    DEALLOCATE ( zcpollen )
  ENDIF
#endif

  !$ser savepoint SedimentationUnittest.DoTracers-out LargeTimeStep=ntstep
  !$ser tracer %all@nnew

  !----------------------------------------------------------------------------
  ! Section 4b: Horizontal diffusion and tendencies of tracers and tke
  !----------------------------------------------------------------------------

  IF (izdebug > 5) THEN
    PRINT *, '         Horizontal diffusion and tendencies of tracers and tke'
  ENDIF

#ifdef COSMOART
  ! CK 20101117 there have been calls to comp_hori_diff here,
  ! but in the current version the diffusion and tendencies are
  ! calculated directly in the comp_hori_diff call (20 lines below)
  ! so I assume this is unnecessary...
#endif

  IF (ltime) CALL get_timings (i_horizontal_diffusion, ntstep, dt, izerror)

  !$ser savepoint VerticalDiffusionUnittest.DoTracers-in LargeTimeStep=ntstep
  !$ser data rho=rho(:,:,:) u_nnow=u(:,:,:,nnow) v_nnow=v(:,:,:,nnow)
  !$ser data ztch=ztch zsqrtgrho_r_s=zsqrtgrho_r_s zkh=zkh pertstoph=pertstoph
  !$ser data t_g=t_g(:,:,nnow)
  !$ser tracer %all@nnow %all@nnew %all#tens %all#surf@nnow

  ! Complete timestep, if required by vertical diffusion
  ! NOTE: no boundary points required as all computations are within the columns
  CALL compute_moisture_divergence
  CALL complete_tendencies_tke
  CALL complete_tendencies_trcr

  !$ser savepoint VerticalDiffusionUnittest.DoTracers-out LargeTimeStep=ntstep
  !$ser tracer %all@nnew
  !$ser data dqvdt=dqvdt qvsflx=qvsflx lhfl_s=lhfl_s

  IF ( l3dturb ) THEN
    IF (  l_dzeta_d_needed ) THEN
      CALL deallocate_dzeta_dxx_at_y
    END IF
    CALL deallocate_K_rho_at_y
    CALL deallocate_beta_xxx_for_3D_diffusion
  END IF

  DEALLOCATE ( wcon, STAT=izstatd )
  DEALLOCATE ( zsqrtgrho_r_s, zsqrtgrho_r_u, zsqrtgrho_r_v, STAT=izstatd )
  DEALLOCATE ( zpia, zpianf, ztmkvm, ztmkvh, ztch, ztcm, zkh,              &
               ztheta, STAT=izstatd )
  IF ( lvertdiff_w )  DEALLOCATE ( ztmkvw, zsqrtgrho_r_w, STAT=izstatd )
  IF ( lprog_tke )    DEALLOCATE ( ztmkvtke, STAT=izstatd )
  IF ( lmoist_turb .AND. l3dturb )  DEALLOCATE ( ztheta_l, STAT=izstatd )

  IF (ltime) CALL get_timings (i_add_tend_moisture, ntstep, dt, izerror)

  !----------------------------------------------------------------------------
  ! Section 5: Horizontal diffusion at end of time step
  !----------------------------------------------------------------------------

  IF ( lhordiff ) THEN
    
    IF (izdebug > 5) THEN
      PRINT *, '         Horizontal diffusion at end of time step'
    ENDIF

    ! NOTE: Requires variable number of boundary lines depending on diffusion
    !       scheme used. See also inside comp_hori_diff() for specifics of LBC
    !       treatment.

    CALL comp_hori_diff( itype_hdiff )

    IF (ltime) CALL get_timings (i_horizontal_diffusion, ntstep, dt, izerror)

    IF ( hd_corr_trcr_bd /= 0.0_wp .OR. hd_corr_trcr_in /= 0.0_wp ) THEN

      IF (ltime) THEN
        IF (ltime_barrier) THEN
          CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
          CALL get_timings (i_barrier_waiting_dyn, ntstep, dt, izerror)
        ENDIF
      ENDIF

      kzdims(1:24)=(/ke,ke,ke1,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                 &
        (22+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, &
         ie, je, kzdims, jstartpar, jendpar, 2, nboundlines, my_cart_neigh, &
        lperi_x, lperi_y, l2dim,                                            &
        20000+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,         &
        u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),     &
        pp(:,:,:,nnew))

      ! NOTE: What about BC for u, v, w, t, qv, qc, pp at timelevel nnew?
      !       What about BC for cpollen at timelevel nnew?

        ! Loop over tracers
        DO iztrcr = 1,trcr_get_ntrcr()

          ! check for each tracer that horizontal diffusion was needed 
          ! (and thus that exchange is required to communicate the 
          !  value of the tracer at tlev=nnew updated in comp_hori_diff)
          IF (izdiff(iztrcr) == T_DIFF_ON) THEN

            ! get pointer to tracer
            CALL trcr_get(izerror, iztrcr, ptr_tlev=nnew, ptr= ztrcr)
            IF (izerror /= 0) THEN
              yzerrmsg = trcr_errorstr(izerror)
              CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
            ENDIF

            ! halo-update
            kzdims(1:24) =                                                    &
              (/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
            CALL exchg_boundaries                                             &
              ( 4     , sendbuf, isendbuflen, imp_reals, icomm_cart,          &
              num_compute, ie, je, kzdims, jstartpar, jendpar, 2, nboundlines,&
              my_cart_neigh, lperi_x, lperi_y, l2dim, 20000+nexch_tag,        &
              ldatatypes, ncomm_type, izerror, yzerrmsg, ztrcr(:,:,:))

            ! NOTE: here some BC should be applied. Which type? Should respect user's setting...

            ! clipping of eventual negative values created by horiz. diff.
            IF ( izclip(iztrcr) == T_CLP_ON ) THEN
              CALL clipping( ztrcr(:,:,:), ie, je, ke )
            ENDIF

          ENDIF
        ENDDO

    ELSE
      kzdims(1:24)=(/ke,ke,ke1,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                 &
        (22+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, &
         ie, je, kzdims, jstartpar, jendpar, 2, nboundlines, my_cart_neigh, &
        lperi_x, lperi_y, l2dim,                                            &
        20000+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,         &
        u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),     &
        pp(:,:,:,nnew))
    ENDIF
    
    !$ser savepoint HorizontalDiffusionUnittest.DoStep-out LargeTimeStep=ntstep
    !$ser data t=t(:,:,:,nnew)-t0(:,:,:) pp=pp(:,:,:,nnew) u=u(:,:,:,nnew)   &
    !$ser&     v=v(:,:,:,nnew) w=w(:,:,:,nnew)
    !$ser tracer %all@nnew

    ! NOTE: what about BC for u, v, w, t, pp at timelevel nnew?

    IF (ltime) CALL get_timings (i_communications_dyn, ntstep, dt, izerror)

  ENDIF

  !----------------------------------------------------------------------------
  ! Section 6: First saturation adjustment for t, qv, qc at time level nnew
  !----------------------------------------------------------------------------

  IF (izdebug > 5) THEN
    PRINT *, '         First saturation adjustment'
  ENDIF

  IF ( ldiabf_lh ) THEN
    ! initialize temperature increment due to latent heat
    IF ( ldiabf_satad ) THEN
      tinc_lh(:,:,:) = - t(:,:,:,nnew)
    ELSE
      tinc_lh(:,:,:) = 0.0_wp
    END IF
  END IF
  
  IF (lcond) THEN

    IF (itype_gscp < 100 .OR. l2mom_satads) THEN

#ifdef NECSX
  CALL collapse(.TRUE., ie, istartpar, iendpar, jstartpar, jendpar)
#endif

      IF ( loutput_diab ) THEN
        ! prepare storing of the latent heat release by the following saturation adjustment:
        ttens_diab(:,:,:) = ttens_diab(:,:,:) - t(:,:,:,nnew) / dt
      ENDIF

      ! The array wcon is allocated again here but just as intermediate storage
      ! for routine satad.
      ALLOCATE (wcon(ie,je,9), STAT=izstata )
      IF ( izstata /= 0 ) THEN
        yzerrmsg="allocation of wcon"
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, yzroutine)
      END IF

      IF (l_satad_dyn_iter) THEN
        ! Determine the maxima of w(:,:,k,:) for reproducible results in the
        ! parallel program
        DO k = 1,ke
          zwmax(k) = MAXVAL ( ABS(w(istartpar:iendpar,jstartpar:jendpar,k,nnew)) )
        ENDDO

        IF ( (lreproduce) .AND. (num_compute > 1) ) THEN
          IF (ltime) THEN
            CALL get_timings (i_dyn_computations, ntstep, dt, izerror)
            IF (ltime_barrier) THEN
              CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
              CALL get_timings (i_barrier_globcom_dyn, ntstep, dt, izerror)
            ENDIF
          ENDIF

          CALL global_values (zwmax, ke, 'MAX', imp_reals, icomm_cart, -1,    &
            yzerrmsg, izerror)

          IF (ltime) CALL get_timings (i_global_communi_dyn, ntstep, dt, izerror)
        ENDIF
      ENDIF

#ifdef NUDGING
      IF (llhn .OR. llhnverif) THEN
        ! For calculation of latent heating rate 
        ! T has to be stored before saturation adjustment
        !$ser savepoint LatentHeatingUnittest.DoSetTt_lheat-in LargeTimeStep=ntstep
        !$ser data tt_lheat_nnew=tt_lheat(:,:,:,nnew) tt_lheat_nnow=tt_lheat(:,:,:,nnow)
        IF (itype_gscp>=3) THEN
          ! in case of itype_gscp==3 (hydci or hydci_pp) get latent heat rate 
          ! of the last time step and store it on the next time step. Both 
          ! Routines are runnig at the end of the time step and therefore after 
          ! the latent heat nudging.
          tt_lheat(:,:,:,nnew)=tt_lheat(:,:,:,nnow)
        ENDIF
        !$ser savepoint LatentHeatingUnittest.DoSetTt_lheat-out LargeTimeStep=ntstep
        !$ser data tt_lheat_nnew=tt_lheat(:,:,:,nnew)

        !$ser savepoint LatentHeatingUnittest.DoGsLheatingAdd-in LargeTimeStep=ntstep
        !$ser data tt_lheat_nnew=tt_lheat(:,:,:,nnew) t=t(:,:,:,nnew)
        CALL get_gs_lheating ('add', 1, ke)
        !$ser savepoint LatentHeatingUnittest.DoGsLheatingAdd-out LargeTimeStep=ntstep
        !$ser data tt_lheat_nnew=tt_lheat(:,:,:,nnew)
      ENDIF
#endif

!US because all humidity variables are now needed in the Runge-Kutta loop,
!   these tracers have to fetched before
!   ! get required tracers (qv and qc at nnew)
!   CALL trcr_get(izerror, 'QV', ptr_tlev = nnew, ptr = qv)
!   IF (izerror /= 0) THEN
!     yzerrmsg = trcr_errorstr(izerror)
!     CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
!   ENDIF
!   CALL trcr_get(izerror, 'QC', ptr_tlev = nnew, ptr = qc)
!   IF (izerror /= 0) THEN
!     yzerrmsg = trcr_errorstr(izerror)
!     CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
!   ENDIF

      DO  k = 1, ke
        DO j = jstartpar, jendpar
          DO i = istartpar, iendpar
            wcon(i,j,1) =  p0(i,j,k) + pp(i,j,k,nnew)
          ENDDO
        ENDDO

        IF (l_satad_dyn_iter) THEN
          ! compute number of iterations according to maximum updraft velocity
          kitpro = 1
          IF (zwmax(k) >=  2.0_wp) kitpro = 2
          IF (zwmax(k) >= 10.0_wp) kitpro = 3
        ELSE
           ! fix number of iteration to 3. Note: In case this number is changed, 
           ! the C++ dycore should be adapted accordingly
          kitpro = 3
        ENDIF

        ! NOTE: saturation adjustment only done in compute domain

        CALL satad ( kitpro, t(:,:,k,nnew), qv(:,:,k),            &
          qc(:,:,k), t(:,:,k,nnow), wcon(:,:,1),                  &
          wcon(:,:,2), wcon(:,:,3), wcon(:,:,4), wcon(:,:,5),     &
          wcon(:,:,6), wcon(:,:,7), wcon(:,:,8), wcon(:,:,9),     &
          b1, b2w, b3, b4w, b234w, rdv, o_m_rdv,                  &
          rvd_m_o, lh_v, cpdr, cp_d,                              &
          ie, je, istartpar, iendpar, jstartpar, jendpar  )
      ENDDO

      DEALLOCATE (wcon, STAT=izstatd)

#ifdef NUDGING
      ! calculate gridscale latent heating rate from the saturation adjustment
      IF (llhn .OR. llhnverif) THEN
        !$ser savepoint LatentHeatingUnittest.DoGsLheatingInc-in LargeTimeStep=ntstep
        !$ser data tt_lheat_nnew=tt_lheat(:,:,:,nnew) t=t(:,:,:,nnew)
        !XL_CHANGE: 'inc' only in compute domain
        CALL get_gs_lheating ('inc', 1, ke)
        !$ser savepoint LatentHeatingUnittest.DoGsLheatingInc-out LargeTimeStep=ntstep
        !$ser data tt_lheat_nnew=tt_lheat(:,:,:,nnew)
      ENDIF
#endif

#ifdef NECSX
      CALL collapse(.FALSE., ie, istartpar, iendpar, jstartpar, jendpar)
#endif

      IF ( loutput_diab ) THEN
        ! complete the storing of the latent heat release by the above saturation adjustment:
        ttens_diab(:,:,:) = ttens_diab(:,:,:) + t(:,:,:,nnew) / dt
      ENDIF

    ENDIF

  ENDIF

  IF ( ldiabf_satad ) THEN
    ! compute temperature increment due to latent heat
    tinc_lh(:,:,:) = tinc_lh(:,:,:) + t(:,:,:,nnew)
  END IF
  
  ! save temperature tendencies due to dynamics for IFS convection scheme
  IF ( itype_conv == 2 ) THEN
    ttens_conv(:,:,:) = ttens_conv(:,:,:) + ( t(:,:,:,nnew) - t(:,:,:,nnow) ) / dt
  END IF

  IF (ltime) CALL get_timings (i_dyn_computations, ntstep, dt, izerror)

!------------------------------------------------------------------------------
! End of module procedure org_runge_kutta
!------------------------------------------------------------------------------

END SUBROUTINE org_runge_kutta

!==============================================================================
!+ Module procedure in "src_runge_kutta" for computing small time step
!------------------------------------------------------------------------------

SUBROUTINE calc_small_timestep( dt, ismtstep, ismtstep_1, ismtstep_2, &
                                dtsmall, dtsmall_1, dtsmall_2 )

!----------------------------------------------------------------------------
!
! Description:
!   Calculation of the small time step sizes dtsmall, _1, _2 and the number of
!   small steps (ismtstep, _1, _2) for each Runge-Kutta-substep
!
! Method:
!
!----------------------------------------------------------------------------

  ! Subroutine arguments:
  ! ---------------------

  REAL    (KIND=wp),        INTENT(IN)  :: dt
  INTEGER (KIND=iintegers), INTENT(OUT) :: ismtstep, ismtstep_1, ismtstep_2
  REAL    (KIND=wp),        INTENT(OUT) :: dtsmall,  dtsmall_1,  dtsmall_2

  ! Local scalars:
  ! -------------

  REAL    (KIND=wp)     ::  &
    z_cs, z_dtsmax, z_dx, z_dy, z_crlat_min, z_dt

  INTEGER (KIND=iintegers) :: &
    z_ismtstep_max, z_ismtstep_sum, izerror

  CHARACTER (LEN=80) ::  &
    yzerrmsg

  INTEGER (KIND=iintegers) :: icalc_version  ! 1 = as before, 0 = Standard RK

  !----------------------------------------------------------------------------
  ! Begin Subroutine calc_small_timestep
  !----------------------------------------------------------------------------

  icalc_version = 0
  !WRITE(*,*) "Subr [calc_small_timestep]: icalc_version =", icalc_version

  z_crlat_min = MINVAL( crlat(:,1) )

  IF (num_compute > 1) THEN
    IF (ltime) THEN
      CALL get_timings (i_dyn_computations, ntstep, dt, izerror)
      IF (ltime_barrier) THEN
        CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
        CALL get_timings (i_barrier_globcom_dyn, ntstep, dt, izerror)
      ENDIF
    ENDIF

    CALL global_values (z_crlat_min, 1, 'MIN', imp_reals, icomm_cart, -1, &
        yzerrmsg, izerror)

    IF (ltime) CALL get_timings (i_global_communi_dyn, ntstep, dt, izerror)
  ENDIF

  z_dt = ABS(dt)   ! in DFI backward stepping, the time step is negative
  z_dx = r_earth * 3.1415927_wp / 180.0_wp * dlon * z_crlat_min
  z_dy = r_earth * 3.1415927_wp / 180.0_wp * dlat

  z_cs        = SQRT( gamma * r_d * 303.0_wp )
  z_dtsmax    = 0.9_wp / z_cs / SQRT( 1.0_wp/z_dx**2 + 1.0_wp/z_dy**2 )
  z_dtsmax    = MIN( z_dtsmax, 32.0_wp )      !MB: reason?

  SELECT CASE(irk_order)
  CASE(1)
    ! = Euler-Forward-scheme
    ismtstep = INT( z_dt/z_dtsmax, iintegers) + 1
    IF ( icalc_version == 1 ) THEN
      IF ( ismtstep == 1 ) ismtstep = 2
      dtsmall  = z_dt/ismtstep
    ELSE
      dtsmall  = z_dt/ismtstep
    ENDIF
    z_ismtstep_max = ismtstep
    z_ismtstep_sum = ismtstep
  CASE(2)
    ismtstep_1 = INT( ( z_dt/2.0_wp) / z_dtsmax, iintegers) + 1
    IF ( ismtstep_1 == 1 ) ismtstep_1 = 2
    ismtstep   = INT( z_dt/z_dtsmax, iintegers) + 1
    IF ( ismtstep == 1 ) ismtstep = 2
    dtsmall_1  = (z_dt/2.0_wp) / ismtstep_1
    dtsmall    = z_dt / ismtstep
    z_ismtstep_max = ismtstep
    z_ismtstep_sum = ismtstep_1 + ismtstep
  CASE(3)
    SELECT CASE(irunge_kutta)
    CASE(1)
      IF ( icalc_version == 1 ) THEN
        ismtstep_1 = INT( (z_dt/3.0_wp) / z_dtsmax, iintegers) + 1
        IF ( ismtstep_1 == 1 ) ismtstep_1 = 2
        ismtstep_2 = INT( (z_dt/2.0_wp) / z_dtsmax, iintegers) + 1
        IF ( ismtstep_2 == 1 ) ismtstep_2 = 2
        ismtstep   = INT( z_dt/z_dtsmax, iintegers) + 1
        IF ( ismtstep == 1 ) ismtstep = 2
        dtsmall_1  = (z_dt/3.0_wp) / ismtstep_1
        dtsmall_2  = (z_dt/2.0_wp) / ismtstep_2
        dtsmall    = z_dt / ismtstep
        z_ismtstep_max = ismtstep
      ELSE
        ! generate a ratio dt/dtsmall so, that numbers 2/3/6, or 4/6/12, ...
        ! of small steps in the 3 RK substeps are used
        ! (remark: 2/3/6 is more efficient than higher numbers,
        ! therefore it could even be advisable to reduce the BIG timestep dt
        ! a bit, to achieve these lower numbers)

        ismtstep = INT( z_dt / z_dtsmax, iintegers)
        ! round up to multiples of 6:
        ismtstep = 6 * INT( ismtstep / 6 + 1, iintegers )

        dtsmall = z_dt / ismtstep

        ismtstep_1 = ismtstep / 3
        dtsmall_1 = dtsmall

        ismtstep_2 = ismtstep / 2
        dtsmall_2 = dtsmall

        z_ismtstep_max = ismtstep
      ENDIF
    CASE(2)
      ismtstep_1 = INT( z_dt / z_dtsmax, iintegers) + 1
      IF ( ismtstep_1 == 1 ) ismtstep_1 = 2
      ismtstep_2 = INT( ( z_dt/4.0_wp) / z_dtsmax, iintegers) + 1
      IF ( ismtstep_2 == 1 ) ismtstep_2 = 2
      ismtstep   = INT( (2.0_wp * z_dt/3.0_wp) / z_dtsmax, iintegers) + 1
      IF ( ismtstep == 1 ) ismtstep = 2
      dtsmall_1  = z_dt / ismtstep_1
      dtsmall_2  = (z_dt/4.0_wp) / ismtstep_2
      dtsmall    = (2.0_wp*z_dt/3.0_wp) / ismtstep
      z_ismtstep_max = ismtstep_1
    END SELECT
    z_ismtstep_sum = ismtstep_1 + ismtstep_2 + ismtstep
  END SELECT

  ! Adapt the sign for the small steps (for DFI)
  IF (dt < 0.0_wp) THEN
    dtsmall   = - dtsmall
    dtsmall_1 = - dtsmall_1
    dtsmall_2 = - dtsmall_2
  ENDIF

  IF (my_cart_id == 0) THEN
    PRINT *,"CS=",z_cs," DT=",dt," DTSMAX=" ,z_dtsmax
    IF      ( irk_order == 1 ) THEN
      PRINT *,"SHORT TIME STEP ", dtsmall, ", number of steps: ", ismtstep
    ELSE IF ( irk_order == 2 ) THEN
      PRINT *,"SHORT TIME STEP of 1. RK step ", dtsmall_1, ", number of steps: ", ismtstep_1
      PRINT *,"SHORT TIME STEP of 2. RK step ", dtsmall,   ", number of steps: ", ismtstep
    ELSE IF ( irk_order == 3 ) THEN
      PRINT *,"SHORT TIME STEP of 1. RK step ", dtsmall_1, ", number of steps: ", ismtstep_1
      PRINT *,"SHORT TIME STEP of 2. RK step ", dtsmall_2, ", number of steps: ", ismtstep_2
      PRINT *,"SHORT TIME STEP of 3. RK step ", dtsmall,   ", number of steps: ", ismtstep
    END IF
    PRINT *,"ismtstep_sum = ",z_ismtstep_sum," ismtstep_max = ", z_ismtstep_max
  ENDIF

END SUBROUTINE calc_small_timestep

!==============================================================================
!==============================================================================

SUBROUTINE check_cfl_horiz_advection ( lapply_Rayleigh_damping, epsray )

!------------------------------------------------------------------------------
!
! Description:
! 1.) Check, if the CFL-criterion for 2-dim. horizontal advection at
!     timelevel nnow is violated.
! 2.) If this is the case then set the Rayleigh damping coefficients
!
!------------------------------------------------------------------------------

!============================================================================

! Declarations:

! Subroutine arguments:
! ---------------------

REAL(KIND=wp),     INTENT(OUT) :: epsray(1:ie, 1:je, 1:ke)
      ! Rayleigh damping coefficient, if CFL-crit. is violated
LOGICAL, INTENT(OUT) ::  lapply_Rayleigh_damping

! Local scalars:
! -------------

INTEGER           :: i, j, k
REAL(KIND=wp)     :: cfl_limit_1D     ! maximum CFL in the 1-dim. case
REAL(KIND=wp)     :: cfl, cfl_max
REAL(KIND=wp)     :: dt_dx, dt_dy
REAL(KIND=wp)     :: epsray0
INTEGER           :: izerror, izdebug
CHARACTER(LEN=80) :: yzerrmsg
LOGICAL           :: lz_crit_point_already_printed

!----------------------------------------------------------------------------
! Begin Subroutine check_cfl_advection
!----------------------------------------------------------------------------

! Initialize, whether debug output shall be done
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

  dt_dx = ABS(dt) / ( Pi * r_earth * dlon / 180.0_wp ) 
          ! correction with 'crlat' see below
  dt_dy = ABS(dt) / ( Pi * r_earth * dlat / 180.0_wp )

  ! 1D-CFL-criteria for different advection schemes
  ! (e.g. Wicker, Skamarock (2002)
  SELECT CASE ( irk_order )
  CASE ( 3 )

    SELECT CASE ( iadv_order )
    CASE(1)
      cfl_limit_1D = 1.24_wp
    CASE(2)
      cfl_limit_1D = 1.72_wp
    CASE(3)
      cfl_limit_1D = 1.61_wp
    CASE(4)
      cfl_limit_1D = 1.25_wp
    CASE(5)
      cfl_limit_1D = 1.42_wp
    CASE(6)
      cfl_limit_1D = 1.08_wp
    CASE DEFAULT
      cfl_limit_1D = 1.0_wp
    END SELECT

  CASE default
    cfl_limit_1D = 1.0_wp
  END SELECT

  ! calculate CFL-criterion
  cfl_max = 0.0_wp
  DO k = 1,ke
    DO j = jstartv,jendv
      DO i = istartu,iendu
        cfl =   ABS( u(i,j,k,nnow) ) * dt_dx / crlat(j,1)  &
              + ABS( v(i,j,k,nnow) ) * dt_dy
        cfl_max = MAX( cfl_max, cfl )
      ENDDO
    ENDDO
  ENDDO

  IF (num_compute > 1) THEN
    IF (ltime) THEN
      CALL get_timings (i_dyn_computations, ntstep, dt, izerror)
      IF (ltime_barrier) THEN
        CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
        CALL get_timings (i_barrier_globcom_dyn, ntstep, dt, izerror)
      ENDIF
    ENDIF

    CALL global_values ( cfl_max, 1, 'MAX', imp_reals, icomm_cart, -1,  &
      yzerrmsg, izerror)

    IF (ltime) CALL get_timings (i_global_communi_dyn, ntstep, dt, izerror)
  ENDIF

  IF (izdebug > 2) THEN
    PRINT *, 'CFL=', cfl_max, ', Max (|u|+|v|) ~ ', cfl_max / dt_dx
    ! this is only a rough estimation for the velocity,
    ! assumption: dt_dx = dt_dy = const. is used
  END IF

  IF ( cfl_max > 0.95_wp * cfl_limit_1D ) THEN

    IF ( my_cart_id == 0 ) THEN
      PRINT*,' !!!!*** WARNING ***!!! CFL-criterion for horizontal &
               &advection is violated'
      PRINT*,' !!!! Max ( |cfl_x| + |cfl_y| ) = ',cfl_max,  &
                            " ( > 0.95 * ", cfl_limit_1D, " )"
    ENDIF

    ! calculate the Rayleigh-damping-coefficient
    ! and search again the point of the CFL-violation for debugging

    lapply_Rayleigh_damping = .TRUE.

    epsray0 = 1.0_wp / (1000.0_wp * ABS(dt))
    epsray0 = epsray0 * (cfl_max - 0.95_wp * cfl_limit_1D) /  &
                        (0.05_wp * cfl_limit_1D)
    ! stability constraint
    epsray0 = SIGN(1.0_wp,dt) * MIN( epsray0, 1.0_wp / (2.0_wp * ABS(dt) ) )

    lz_crit_point_already_printed = .FALSE.

    epsray(:,:,:)=0.0_wp

    DO k = 1,ke
      DO j = jstartv,jendv
        DO i = istartu,iendu

          cfl =   ABS( u(i,j,k,nnow) ) * dt_dx / crlat(j,1)  &
                + ABS( v(i,j,k,nnow) ) * dt_dy

          ! reduce only high velocities:
          IF ( cfl > 0.90_wp * cfl_limit_1D ) THEN
            epsray(i,j,k) = epsray0
          ELSE
            epsray(i,j,k) = 0.0_wp
          END IF

          IF ( .NOT.(lz_crit_point_already_printed) .AND.   &
               (cfl >= 0.9999_wp * cfl_max )   ) THEN
            PRINT *, ' global indices: i=', i_global(i),                  &
                                    ', j=', j_global(j), ', k=', k
            PRINT *, "|cfl_x| + |cfl_y| = ", cfl, ", u=", u(i,j,k,nnow),  &
                                                  ", v=", v(i,j,k,nnow)
            PRINT *, "epsray0 = ", epsray0
            lz_crit_point_already_printed = .TRUE.
               ! printing of only one grid point with
               ! extremely high velocities is sufficient
          END IF

        ENDDO
      ENDDO
    ENDDO

  ELSE

    lapply_Rayleigh_damping = .FALSE.
    !epsray(:,:,:) = 0.0_wp

  ENDIF

END SUBROUTINE check_cfl_horiz_advection

!==============================================================================
!==============================================================================

SUBROUTINE coriolis( utens, vtens, wtens, nn ) 

!-----------------------------------------------------
!
! Description:
! Add Coriolis force to u- and v-tendencies
!
! Input:
!   u(:,:,:,nn), v(:,:,:,nn), w(:,:,:,nn)  at timelevel nn
!   fc, fccos
!
! Input- and Output:
!   utens, vtens, wtens
!
!------------------------------------------------------
  REAL(KIND=wp),     INTENT(INOUT) :: utens(ie,je,ke), vtens(ie,je,ke)
  REAL(KIND=wp),     INTENT(INOUT) :: wtens(ie,je,ke1)
  INTEGER (KIND=iintegers), INTENT(IN) :: nn

  INTEGER          :: i,j,k
  REAL(KIND=wp)    :: z_fv_north, z_fv_south, zfq
  REAL(KIND=wp)    :: z_fcw_west, z_fcw_east, zfcw 
  REAL(KIND=wp)    :: z_fu_west,  z_fu_east
  REAL(KIND=wp)    :: zfcu

!$ser savepoint CoriolisUnittest.Apply-in LargeTimeStep=ntstep
!$ser data u_nnow=u(:,:,:,nnow) v_nnow=v(:,:,:,nnow) u_tens=utens v_tens=vtens

  DO k = 1 , ke
    DO j = jstartu, jendu
!CDIR ON_ADB(fc)
!CDIR ON_ADB(v)
      DO i = istartu, iendu
        z_fv_north = fc(i,j)   * ( v(i,j  ,k,nn) + v(i+1,j  ,k,nn) )
        z_fv_south = fc(i,j-1) * ( v(i,j-1,k,nn) + v(i+1,j-1,k,nn) )
        zfq = 0.25_wp * ( z_fv_north + z_fv_south )
        utens(i,j,k) = utens(i,j,k) + zfq
      ENDDO
    ENDDO

    DO j = jstartv, jendv
!CDIR ON_ADB(fc)
!CDIR ON_ADB(u)
      DO i = istartv, iendv
        z_fu_east = fc(i,j)   * ( u(i  ,j,k,nn) + u(i  ,j+1,k,nn) )
        z_fu_west = fc(i-1,j) * ( u(i-1,j,k,nn) + u(i-1,j+1,k,nn) )
        zfq = 0.25_wp * ( z_fu_east + z_fu_west )
        vtens(i,j,k) = vtens(i,j,k) - zfq 
      ENDDO
    ENDDO
  END DO

  IF (lcori_deep) THEN
    ! additional 2*Omega*cos(phi) * w -term in the momentum equation for u
    DO k = 1 , ke
      DO j = jstartu, jendu
        DO i = istartu, iendu
          z_fcw_west = fccos(i  ,j) * ( w(i  ,j,k,nn) + w(i  ,j,k+1,nn) )
          z_fcw_east = fccos(i+1,j) * ( w(i+1,j,k,nn) + w(i+1,j,k+1,nn) )
          zfcw = 0.25_wp * ( z_fcw_west + z_fcw_east )
          utens(i,j,k) = utens(i,j,k) - zfcw
        ENDDO
      ENDDO
    END DO

    ! additional 2*Omega*cos(phi) * u -term in the momentum equation for w
    DO k = 2 , ke
      ! no calculations for wtens(:,:,1/ke) due to boundary conditions
      DO j = jstart, jend
        DO i = istart, iend
          zfcu = 0.25_wp * fccos(i,j)                   &
               * ( u(i-1,j,k-1,nn) + u(i,j,k-1,nn) +        &
                   u(i-1,j,k  ,nn) + u(i,j,k  ,nn) )
          wtens(i,j,k) = wtens(i,j,k) + zfcu
        ENDDO
      ENDDO
    ENDDO
  ENDIF

!$ser savepoint CoriolisUnittest.Apply-out LargeTimeStep=ntstep
!$ser data u_tens=utens v_tens=vtens

END SUBROUTINE coriolis

!==============================================================================
!==============================================================================

SUBROUTINE finalize_runge_kutta

!-----------------------------------------------------
!
! Description:
!   Do a proper clean up of the Runge-Kutta method
!
!------------------------------------------------------

  IF ( itype_fast_waves == 2 ) CALL finalize_fast_waves_sc


END SUBROUTINE finalize_runge_kutta

!==============================================================================

END MODULE src_runge_kutta
