!+ Source module for the inland lake water model
!------------------------------------------------------------------------------

MODULE src_flake

!------------------------------------------------------------------------------
!
! Description:
!
!  The main program unit of the lake model FLake,
!  containing most of the FLake procedures.
!  Most FLake variables and local parameters are declared.
!
!  FLake (Fresh-water Lake) is a lake model capable of predicting the surface
!  temperature in lakes of various depth on the time scales from a few hours
!  to a year.
!  The model is based on a two-layer parametric representation of the evolving
!  temperature profile, where the structure of the stratified layer between the
!  upper mixed layer and the basin bottom, the lake thermocline, is described
!  using the concept of self-similarity of the temperature-depth curve.
!  The concept was put forward by Kitaigorodskii and Miropolsky (1970) to
!  describe the vertical temperature structure of the oceanic seasonal
!  thermocline. It has been successfully used in geophysical applications.
!  The concept of self-similarity of the evolving temperature profile is also
!  used to describe the vertical structure of the thermally active upper layer
!  of bottom sediments and of the ice and snow cover.
!
!  The lake model incorporates the heat budget equations for the four layers
!  in question, viz., snow, ice, water and bottom sediments, developed with
!  due regard for the vertically distributed character of solar radiation
!  heating. The entrainment equation that incorporates the spin-up term
!  is used to compute the depth of a convectively-mixed layer.
!  A relaxation-type equation is used to compute the wind-mixed layer depth in
!  stable and neutral stratification, where a multi-limit formulation for the
!  equilibrium mixed-layer depth accounts for the effects of the earth's rotation,
!  of the surface buoyancy flux and of the static stability in the thermocline.
!  The equations for the mixed-layer depth are developed with due regard for
!  the volumetric character of the radiation heating.
!  Simple thermodynamic arguments are invoked to develop
!  the evolution equations for the ice and snow thickness.
!  The heat flux through the water-bottom sediment interface is computed,
!  using a parameterization proposed by Golosov et al. (1998).
!  The heat flux trough the air-water interface
!  (or through the air-ice or air-snow interface)
!  is provided by the driving atmospheric model.
!
!  Empirical constants and parameters of the lake model are estimated, using
!  independent empirical and numerical data. They should not be re-evaluated
!  when the model is applied to a particular lake. The only lake-specific
!  parameters are the lake depth, the optical characteristics of lake water,
!  the temperature at the lower boundary of the thermally active layer of bottom
!  sediments, and the depth of that layer.
!
!  A detailed description of the lake model FLake is given in
!
!  Mironov, D. V., 2008: 
!  Parameterization of lakes in numerical weather prediction. Description of a lake model. 
!  COSMO Technical Report, No. 11, Deutscher Wetterdienst, Offenbach am Main, Germany, 41 pp.
!
!  Mironov, D., E. Heise, E. Kourzeneva, B. Ritter, N. Schneider, and A. Terzhevik, 2010:
!  Implementation of the lake parameterisation scheme FLake
!  into the numerical weather prediction model COSMO.
!  Boreal Env. Res., 15, 218-230.
!
!  In the present configuration, the bottom sediment module of FLake is switched off
!  and the heat flux at the water-bottom sediment interface is set to zero. 
!  Although the snowfall rate is provided by COSMO, snow over lake ice is not considered explicitly.
!  The effect of snow is accounted for implicitly (parametrically) 
!  through changes in the surface albedo with respect to solar radiation.
!
!
!  Interface to the FLake Model:
!  This module contains two subroutines, "flake_init" and "flake_interface",
!  called by the COSMO routine "organize_physics".
!  In "flake_init", the FLake variables are initialized.
!  The routine "flake_interface" is the communication routine between 
!  the COSMO model and "flake_driver" (see "flake_driver" for further comments).
!
!  Lines embraced with "!_tmp" contain temporary parts of the code.
!  Lines embraced/marked with "!_dev" may be replaced
!  as improved parameterizations are developed and tested.
!  Lines embraced/marked with "!_dm" are DM's comments
!  that may be helpful to a user.
!  Lines embraced/marked with "!_dbg" are used
!  for debugging purposes only.
!
!
! Current Code Owner: DWD, Dmitrii Mironov
!  Phone:  +49-69-8062 2705
!  Fax:    +49-69-8062 3721
!  E-mail: dmitrii.mironov@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 3.18       2006/03/03 Dmitrii Mironov
!  Initial release
! 3.21       2006/12/04 Dmitrii Mironov, Ulrich Schaettler
!  General update of the scheme
! V4_11        2009/11/30 Ekaterina Machulskaya, Jan-Peter Schulz
!  Adaptations for multi-layer snow model
!  Eliminate option for cold start (is put to INT2LM) (US)
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_15        2010/11/19 Ulrich Schaettler
!  The field t_s_lake is removed again after adaptations in the SST Analysis
! V4_16        2010/12/07 Dmitrii Mironov
!  Adaptation to the two time level scheme, update of I/O (SUBROUTINE flake_init)
! V4_21        2011/12/06 Dmitrii Mironov
!  Changes for better vectorization
! V4_23        2012/05/10 Burkhardt Rockel (CLM)
!  Restrict some writing information to "debug output"
!  Set  t_snow_mult = t_snow only for grid points depth_lk > 0.0,
!   otherwise program may crash later on in organize_dynamics
! V4_24        2012/06/22 Dmitrii Mironov
!  Added consistency checks at the beginning, because values might be inconsistent
!    due to grib packing
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer
!  Replaced qx-variables by using them from the tracer module
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  Introduced MESSy interface
! V4_28        2013/07/12 Dmitrii Mironov
!  Modified some comments
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
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

USE data_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! physical constants and related variables
! -------------------------------------------
    lh_v,         & ! latent heat of vapourization
    lh_f,         & ! latent heat of fusion
    cp_d,         & ! specific heat of dry air at constant pressure
    rdocp           ! r_d / cp_d

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------
    ie         ,  & ! number of grid points in zonal direction
    je         ,  & ! number of grid points in meridional direction
    ke         ,  & ! number of grid points in vertical direction
    ke1        ,  & ! ke + 1
    ke_soil    ,  & ! number of layers in multi-layer soil model
    ke_snow    ,  & ! number of layers in multi-layer snow model

! start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from
!    the other ones because of the use of the staggered Arakawa-C-grid.
!
    istartpar  ,  & ! start index for computations in the parallel program
    iendpar    ,  & ! end index for computations in the parallel program
    jstartpar  ,  & ! start index for computations in the parallel program
    jendpar    ,  & ! end index for computations in the parallel program

! variables for the time discretization and related variables
! --------------------------------------------------------------
    dt         ,  & ! long time-step
    dt2        ,  & ! dt*2.

! 8. Organizational variables to handle the COSMO humidity tracers
! ----------------------------------------------------------------
    idt_qv

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! external parameter fields                                           (unit)
! ----------------------------
    p0         ,    & ! base state pressure                           ( pa  )
    fc         ,    & ! coriolis-parameter                            ( 1/s )
    llandmask  ,    & ! landpoint mask                                (  -  )

! external parameter fields of the lake model                         (unit)
! ----------------------------
    fr_lake    ,    & ! lake fraction in a grid element [0,1]         (  -  )
    depth_lk   ,    & ! lake depth                                    (  m  )
    fetch_lk   ,    & ! wind fetch over lake                          (  m  )
    dp_bs_lk   ,    & ! depth of the thermally active layer
                      ! of bottom sediments                           (  m  )
    t_bs_lk    ,    & ! climatological temperature at the bottom of
                      ! the thermally active layer of sediments       (  K  )
    gamso_lk   ,    & ! attenuation coefficient for
                      ! solar radiation in lake water                 ( 1/m )
! prognostic variables                                                (unit)
! -----------------------
    u          ,    & ! zonal wind speed                              ( m/s )
    v          ,    & ! meridional wind speed                         ( m/s )
    t          ,    & ! temperature                                   (  k  )
    pp         ,    & ! deviation from the reference pressure         ( pa  )

! four-dimensional fields (x,y,z,t)
! ---------------------------------
    t_so       ,    & ! multi-layer soil temperature                  (  K  )

! fields for surface values and soil model variables                  (unit )
! -----------------------------------------------------
    ps         ,    & ! surface pressure                              ( pa  )
    t_s        ,    & ! temperature of the ground surface             (  K  )
    t_snow     ,    & ! temperature of the snow-surface               (  k  )
    t_snow_mult,    & ! temperature of the snow-surface               (  k  )
    t_ice      ,    & ! temperature at the snow-ice or
                      ! air-ice interface                             (  K  )
    h_snow     ,    & ! height of snow                                (  m  )
    h_ice      ,    & ! ice thickness                                 (  m  )
    t_g        ,    & ! weighted surface temperature                  (  K  )
    qv_s       ,    & ! specific water vapor content at the surface   (kg/kg)

! fields computed by the turbulence scheme                            (unit )
! -----------------------------------------------------
    tch        ,    & ! turbulent transfer coefficient for heat       ( -- )
    tcm        ,    & ! turbulent transfer coefficient for momentum   ( -- )

! fields computd by the radiation scheme                              (unit )
! -----------------------------------------------------
    sobs       ,    & ! solar radiation at the ground                 ( W/m2)
    thbs       ,    & ! thermal radiation at the ground               ( W/m2)

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

!------------------------------------------------------------------------------

USE data_flake,         ONLY : &

  !  flake_configure:

  lflk_botsed_use  , & ! .TRUE. indicates that bottom-sediment scheme is used
  rflk_depth_bs_ref, & ! Reference value of depth of thermally active layer 
                       ! of bottom sediments

  ! derived types
  opticpar_medium  , & ! derived type 

  ! parameters
  c_cbl_1          , & ! Constant in the CBL entrainment equation
  c_cbl_2          , & ! Constant in the CBL entrainment equation
  c_sbl_ZM_n       , & ! Constant in the ZM1996 equation for the equilibrium 
                       ! SBL depth
  c_sbl_ZM_s       , & ! Constant in the ZM1996 equation for the equilibrium 
                       ! SBL depth
  c_sbl_ZM_i       , & ! Constant in the ZM1996 equation for the equilibrium 
                       ! SBL depth
  c_relax_h        , & ! Constant in the relaxation equation for the SBL depth
  c_relax_C        , & ! Constant in the relaxation equation for the shape 
                       !   factor with respect to the temperature profile in 
                       !   the thermocline

  C_T_min          , & ! Minimum value of the shape factor C_T (thermocline)
  C_T_max          , & ! Maximum value of the shape factor C_T (thermocline)
  Phi_T_pr0_1      , & ! Constant in the expression for the T shape-function derivative
  Phi_T_pr0_2      , & ! Constant in the expression for the T shape-function derivative
  C_TT_1           , & ! Constant in the expression for C_TT (thermocline)
  C_TT_2           , & ! Constant in the expression for C_TT (thermocline)
  C_B1             , & ! Shape factor (upper layer of bottom sediments)
  C_B2             , & ! Shape factor (lower layer of bottom sediments)
  Phi_B1_pr0       , & ! B1 shape-function derivative
  C_S_lin          , & ! Shape factor (linear temperature profile in the snow layer)
  Phi_S_pr0_lin    , & ! S shape-function derivative (linear profile)
  C_I_lin          , & ! Shape factor (linear temperature profile in the ice layer)
  Phi_I_pr0_lin    , & ! I shape-function derivative (linear profile)
  Phi_I_pr1_lin    , & ! I shape-function derivative (linear profile)
  Phi_I_ast_MR     , & ! Constant in the MR2004 expression for I shape factor
  C_I_MR           , & ! Constant in the MR2004 expression for I shape factor
  H_Ice_max        , & ! Maximum ice tickness in the Mironov/Ritter ice model [m]

  h_Snow_min_flk   , & ! Minimum snow thickness [m]
  h_Ice_min_flk    , & ! Minimum ice thickness [m]
  h_ML_min_flk     , & ! Minimum mixed-layer depth [m]
  h_ML_max_flk     , & ! Maximum mixed-layer depth [m]
  H_B1_min_flk     , & ! Minimum thickness of the upper layer of bottom sediments [m]
  u_star_min_flk   , & ! Minimum value of the surface friction velocity [m s^{-1}]

  c_small_flk      , & ! A small number
  c_maxearg_flk        ! Maximum value of the EXP function argument [-]

USE data_flake,         ONLY : &

  tpl_grav         , & ! Acceleration due to gravity [m s^{-2}]
  tpl_T_r          , & ! Temperature of maximum density of fresh water [K]
  tpl_T_f          , & ! Fresh water freezing point [K]
  tpl_a_T          , & ! Constant in the fresh-water equation of state [K^{-2}]
  tpl_rho_w_r      , & ! Maximum density of fresh water [kg m^{-3}]
  tpl_rho_I        , & ! Density of ice [kg m^{-3}]
  tpl_rho_S_min    , & ! Minimum snow density [kg m^{-3}]
  tpl_rho_S_max    , & ! Maximum snow density [kg m^{-3}]
  tpl_Gamma_rho_S  , & ! Empirical parameter [kg m^{-4}] in the expression for 
                       !    the snow density
  tpl_L_f          , & ! Latent heat of fusion     [J kg^{-1}]
  tpl_c_w          , & ! Specific heat of water    [J kg^{-1} K^{-1}]
  tpl_c_I          , & ! Specific heat of ice      [J kg^{-1} K^{-1}]
  tpl_c_S          , & ! Specific heat of snow     [J kg^{-1} K^{-1}]
  tpl_kappa_w      , & ! Molecular heat conductivity of water 
                       !   [J m^{-1} s^{-1} K^{-1}]
  tpl_kappa_I      , & ! Molecular heat conductivity of ice
                       !   [J m^{-1} s^{-1} K^{-1}]
  tpl_kappa_S_min  , & ! Minimum molecular heat conductivity of snow
                       !   [J m^{-1} s^{-1} K^{-1}]
  tpl_kappa_S_max  , & ! Maximum molecular heat conductivity of snow
                       !   [J m^{-1} s^{-1} K^{-1}]
  tpl_Gamma_kappa_S, & ! Empirical parameter in expression for the 
                       !   snow heat conductivity
                       !   [J m^{-2} s^{-1} K^{-1}]

  ! Optical characteristics for water, ice and snow.
  opticpar_water_ref      , & ! Water (reference)
! opticpar_water_trans    , & ! Transparent Water (two-band)
! opticpar_whiteice_ref   , & ! White ice
! opticpar_blueice_ref    , & ! Blue ice
! opticpar_drysnow_ref    , & ! Dry snow
! opticpar_meltingsnow_ref, & ! Melting snow
  opticpar_ice_opaque     , & ! Opaque ice
  opticpar_snow_opaque        ! Opaque snow

!------------------------------------------------------------------------------

USE data_parallel   , ONLY :   &
  my_cart_id       ! rank of this subdomain in the cartesian communicator

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! controlling the physics
! --------------------------
    lmulti_layer   , & ! run multi-layer soil model
    lmulti_snow    , & ! run multi-layer snow model
    lseaice        , & ! forecast with sea ice model

! start and end of the forecast
! --------------------------------
    ntstep       ,  & ! actual time step
    nold         ,  & ! corresponds to ntstep - 1
    nnow         ,  & ! corresponds to ntstep
    nnew         ,  & ! corresponds to ntstep + 1

! controlling the dynamics
! ---------------------------
    l2tls        ,  & ! time integration by two time level RK-scheme (.TRUE.)
                      ! or by three time level KW-scheme (.FALSE.)

! 12. controlling verbosity of debug output
! -----------------------------------------
    idbg_level,   & ! to control the verbosity of debug output
    ldebug_soi,   & ! if .TRUE., debug output for soil and surface
    lprintdeb_all   ! .TRUE.:  all tasks print debug output
                    ! .FALSE.: only task 0 prints debug output

!------------------------------------------------------------------------------

USE src_flake_sfcflx, ONLY:    &
  SfcFlx_rhoair                   ! Function, returns the air density

!------------------------------------------------------------------------------
  
USE environment,      ONLY: model_abort
  
!------------------------------------------------------------------------------
  
USE src_tracer,       ONLY: trcr_get, trcr_errorstr

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations
!
!  The variables declared below are accessible to all program units of the
!  MODULE flake.
!  These are basically the quantities computed by FLake.
!  All variables declared below have a suffix "flk".

!  FLake variables of type REAL

!  Temperatures at the previous time step ("p")
!  and the updated temperatures ("n")
REAL (KIND = wp)     ::           &
  T_mnw_p_flk, T_mnw_n_flk      , & ! Mean temperature of the water column [K]
  T_snow_p_flk, T_snow_n_flk    , & ! Temperature at the air-snow interface [K]
  T_ice_p_flk, T_ice_n_flk      , & ! Temperature at the snow-ice or air-ice
                                    ! interface [K]
  T_wML_p_flk, T_wML_n_flk      , & ! Mixed-layer temperature [K]
  T_bot_p_flk, T_bot_n_flk      , & ! Temperature at the water-bottom sediment
                                    ! interface [K]
  T_B1_p_flk, T_B1_n_flk            ! Temperature at the bottom of the upper
                                    ! layer of the sediments [K]

!  Thickness of various layers at the previous time step ("p")
!  and the updated values ("n")
REAL (KIND = wp)     ::           &
  h_snow_p_flk, h_snow_n_flk    , & ! Snow thickness [m]
  h_ice_p_flk, h_ice_n_flk      , & ! Ice thickness [m]
  h_ML_p_flk, h_ML_n_flk        , & ! Thickness of the mixed-layer [m]
  H_B1_p_flk, H_B1_n_flk            ! Thickness of the upper layer of bottom
                                    !  sediments [m]

!  The shape factor(s) at the previous time step ("p")
!  and the updated value(s) ("n")
REAL (KIND = wp)     ::           &
  C_T_p_flk, C_T_n_flk          , & ! Shape factor (thermocline)
  C_TT_flk                      , & ! Dimensionless parameter (thermocline)
  C_Q_flk                       , & ! Shape factor with respect to the heat
                                    ! flux (thermocline)
  C_I_flk                       , & ! Shape factor (ice)
  C_S_flk                           ! Shape factor (snow)

!  Derivatives of the shape functions
REAL (KIND = wp)     ::           &
  Phi_T_pr0_flk                 , & ! d\Phi_T(0)/d\zeta   (thermocline)
  Phi_I_pr0_flk                 , & ! d\Phi_I(0)/d\zeta_I (ice)
  Phi_I_pr1_flk                 , & ! d\Phi_I(1)/d\zeta_I (ice)
  Phi_S_pr0_flk                     ! d\Phi_S(0)/d\zeta_S (snow)

!  Heat and radiation fluxes
REAL (KIND = wp)     ::           &
  Q_snow_flk                    , & ! Heat flux through the air-snow interface
                                    ! [W m^{-2}]
  Q_ice_flk                     , & ! Heat flux through the snow-ice or air-ice
                                    ! interface [W m^{-2}]
  Q_w_flk                       , & ! Heat flux through the ice-water or
                                    ! air-water interface [W m^{-2}]
  Q_bot_flk                     , & ! Heat flux through the water-bottom
                                    ! sediment interface [W m^{-2}]
  I_atm_flk                     , & ! Radiation flux at the lower boundary of
                                    ! the atmosphere [W m^{-2}],
                                    ! i.e. the incident radiation flux with no
                                    ! regard for the surface albedo.
  I_snow_flk                    , & ! Radiation flux through the air-snow
                                    ! interface [W m^{-2}]
  I_ice_flk                     , & ! Radiation flux through the snow-ice or
                                    ! air-ice interface [W m^{-2}]
  I_w_flk                       , & ! Radiation flux through the ice-water or
                                    ! air-water interface [W m^{-2}]
  I_h_flk                       , & ! Radiation flux through the mixed-layer-
                                    ! thermocline interface [W m^{-2}]
  I_bot_flk                     , & ! Radiation flux through the water-bottom
                                    ! sediment interface [W m^{-2}]
  I_intm_0_h_flk                , & ! Mean radiation flux over the mixed layer
                                    ! [W m^{-2}]
  I_intm_h_D_flk                , & ! Mean radiation flux over the thermocline
                                    ! [W m^{-2}]
  Q_star_flk                        ! A generalized heat flux scale [W m^{-2}]

!  Velocity scales
REAL (KIND = wp)     ::           &
  u_star_w_flk                  , & ! Friction velocity in the surface layer of
                                    ! lake water [m s^{-1}]
  w_star_sfc_flk                    ! Convective velocity scale, using a
                                    ! generalized heat flux scale [m s^{-1}]

!  The rate of snow accumulation
REAL (KIND = wp)     ::           &
  dMsnowdt_flk                      ! The rate of snow accumulation
                                    ! [kg m^{-2} s^{-1}]

!==============================================================================

CONTAINS

!==============================================================================
!  Initialize the lake model
!------------------------------------------------------------------------------

SUBROUTINE flake_init

!------------------------------------------------------------------------------
!
! Description:
!
!  This routine initializes prognostic variables of the lake model FLake.
!  At present, the bottom-sediment module of FLake is switched off so that 
!  the external parameters of the bottom-sediment module,
!  namely, the depth of the thermally active layer of the sediment
!  and the climatological temperature at the bottom of that layer,
!  are not required.
!  These parameters are set to their reference values and 
!  are not part of the input (not read from the two-dimensional GRIB files).
!  The same applies to the attenuation coefficient with respect to solar 
!  radiation and to the wind fetch. Default values of these variables are used.
!  Three FLake prognostic variables, namely, 
!  the snow thickness (always zero as snow over lake ice is not considered 
!  explicitly), thickness of the upper layer of bottom sediment and 
!  the temperature at the outer edge of that layer 
!  (these are equal to their reference values as the bottom-sediment module is 
!  switched off), are not included into the COSMO IO list.
!  These variable are handled internally.
!
!==============================================================================
!
! Declarations

!  Local variables of type INTEGER
INTEGER (KIND = iintegers) :: &
  i, j, ks, nt              , & ! Loop indices
  istarts                   , & ! Start index for x-direction
  iends                     , & ! End   index for x-direction
  jstarts                   , & ! Start index for y-direction
  jends                     , & ! End   index for y-directin
  izdebug                   , & ! for debug output
  nztlev                    , & ! number of time-levels for prognostic variables
  nx                            ! Time-level for initialization (nnow or nnew)

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

! Set method of debug output
IF (ldebug_soi) THEN
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

! Set start and end indices for DO loops in the horizontal directions

istarts = istartpar
iends   = iendpar
jstarts = jstartpar
jends   = jendpar

! Set time level for initialization
! Notice that for two time level scheme (Runge-Kutta) nx=nnew should be used
IF(l2tls) THEN
  nx = nnew
  nztlev = 2
ELSE
  nx = nnow
  nztlev = 3
ENDIF

! Lake fraction and lake depth are input values (read from GRIB files).
! Set other FLake external parameters to reference values.
DO j = jstarts, jends
!CDIR NOMOVE
DO i = istarts, iends
  fetch_lk(i,j) = 1.E4_wp           ! Use a constant long fetch
  dp_bs_lk(i,j) = rflk_depth_bs_ref     ! Reference value
  t_bs_lk (i,j) = tpl_T_r               ! Reference value
  gamso_lk(i,j) = opticpar_water_ref%extincoef_optic(1)
                                        ! Reference value
END DO
END DO

! Since a loss of accuracy may occur during IO
! (e.g. due to GRIB encoding/decoding),
! FLake prognostic variables must be checked for consistency.
! That is, the relations between
! "t_wml_lk", "t_mnw_lk", "t_bot_lk", "h_ml_lk", "c_t_lk", "h_ice" and "t_ice",
! suggested by the assumed shape (self-similarity) of the temperature-depth curve,
! must hold or should be enforced.
DO j = jstarts, jends
!CDIR NOMOVE
DO i = istarts, iends
  IF(depth_lk(i,j) > 0.0_wp) THEN
    ! Limit the mixed-layer depth and temperature-profile shape factor
    h_ml_lk(i,j,nx) = MIN(depth_lk(i,j), MAX(h_ml_lk(i,j,nx), 0._wp))
    c_t_lk(i,j,nx)  = MIN(C_T_max, MAX(c_t_lk(i,j,nx), C_T_min))
    IF(h_ice(i,j,nx) >= h_Ice_min_flk) THEN
    ! Ice-covered lake-type grid boxes
      h_ice(i,j,nx)    = MIN(h_ice(i,j,nx), H_Ice_max)
      t_ice(i,j,nx)    = MIN(t_ice(i,j,nx), tpl_T_f)
      t_wml_lk(i,j,nx) = tpl_T_f
      t_mnw_lk(i,j,nx) = MAX(tpl_T_f, MIN(t_mnw_lk(i,j,nx), tpl_T_r))
      t_bot_lk(i,j,nx) = MAX(tpl_T_f, MIN(t_bot_lk(i,j,nx), tpl_T_r))
      ! Mixing down to the lake bottom or "inverse" stratification,
      ! reset mixed-layer depth, shape factor, mean temperature, and bottom temperature.
      IF( (h_ml_lk(i,j,nx) >= depth_lk(i,j)-h_ML_min_flk) .OR.  &
          (t_bot_lk(i,j,nx) <= t_wml_lk(i,j,nx)) ) THEN
        h_ml_lk(i,j,nx)  = 0._wp
        c_t_lk(i,j,nx)   = C_T_min
        t_mnw_lk(i,j,nx) = t_wml_lk(i,j,nx)
        t_bot_lk(i,j,nx) = t_wml_lk(i,j,nx)
      END IF
    ELSE
    ! Ice-free lake-type grid boxes
      h_ice(i,j,nx)    = 0.0_wp
      t_ice(i,j,nx)    = tpl_T_f
      t_wml_lk(i,j,nx) = MAX(tpl_T_f, t_wml_lk(i,j,nx))
      t_mnw_lk(i,j,nx) = MAX(tpl_T_f, t_mnw_lk(i,j,nx))
      t_bot_lk(i,j,nx) = MAX(tpl_T_f, t_bot_lk(i,j,nx))
      ! Avoid temperature of maximum density crossover
      IF( ((t_wml_lk(i,j,nx) <= tpl_T_r) .AND. (t_bot_lk(i,j,nx) > tpl_T_r)) .OR.  &
           ((t_wml_lk(i,j,nx) >= tpl_T_r) .AND. (t_bot_lk(i,j,nx) < tpl_T_r)) )    &
        t_bot_lk(i,j,nx) = tpl_T_r
      ! Mixing down to the lake bottom or "inverse" stratification,
      ! reset mixed-layer depth, shape factor, mised-layer temperature, and bottom temperature.
      IF( (h_ml_lk(i,j,nx) >= depth_lk(i,j)-h_ML_min_flk) .OR.                               &
           ((t_wml_lk(i,j,nx) <= tpl_T_r) .AND. (t_bot_lk(i,j,nx) < t_wml_lk(i,j,nx))) .OR.  &
           ((t_wml_lk(i,j,nx) >= tpl_T_r) .AND. (t_bot_lk(i,j,nx) > t_wml_lk(i,j,nx))) ) THEN
        h_ml_lk(i,j,nx)  = depth_lk(i,j)
        c_t_lk(i,j,nx)   = C_T_min
        t_wml_lk(i,j,nx) = t_mnw_lk(i,j,nx)
        t_bot_lk(i,j,nx) = t_mnw_lk(i,j,nx)
      END IF
    END IF
    ! Snow over lake ice is not considered explicitly,
    ! set "h_snow" to zero and "t_snow" to the ice surface temperature.
    h_snow(i,j,nx)   = 0.0_wp
    t_snow(i,j,nx)   = t_ice(i,j,nx)
  END IF
END DO
END DO


! The outermost loop is over the time levels
DO nt = 1, nztlev
DO j = jstarts, jends
!CDIR NOMOVE
DO i = istarts, iends
  Lake_points: IF(depth_lk(i,j) > 0.0_wp) THEN
    ! At lake points, set "t_s", "t_g" and "t_so" to FLake values
    t_s(i,j,nt) = t_s(i,j,nx)
    t_g(i,j,nt) = t_s(i,j,nx)
    IF (lmulti_layer) THEN
      ! Save lake surface temperature in "t_so(:,:,0,:)"
      ! if the multi-layer soil model is used
      t_so(i,j,0,nt) = t_s(i,j,nx)
    ENDIF
    ! Set FLake variables at all time levels
    ! ("nnew, "nnow", and for three time level scheme "nold")
    ! to their values at "nx" (nx=nnew for two time level scheme, see above).
    t_wml_lk(i,j,nt) = t_wml_lk(i,j,nx)
    t_mnw_lk(i,j,nt) = t_mnw_lk(i,j,nx)
    t_bot_lk(i,j,nt) = t_bot_lk(i,j,nx)
    c_t_lk  (i,j,nt) = c_t_lk  (i,j,nx)
    h_ml_lk (i,j,nt) = h_ml_lk (i,j,nx)
    h_ice   (i,j,nt) = h_ice   (i,j,nx)
    t_ice   (i,j,nt) = t_ice   (i,j,nx)
    ! As snow over lake ice is not considered explicitly,
    ! "h_snow" is zero
    h_snow  (i,j,nt) = h_snow  (i,j,nx)
    ! and "t_snow" is equal to the ice surface temperature
    t_snow  (i,j,nt) = t_snow  (i,j,nx)
  ELSE Lake_points
    ! Set FLake variables at non-lake (sea and land) points at their reference
    ! values at all time levels. Then, all FLake variables at land and sea points
    ! remain equal to their initial values during the entire COSMO run
    ! independent of permutation of time indices.
    t_wml_lk(i,j,nt) = tpl_T_r
    t_mnw_lk(i,j,nt) = tpl_T_r
    t_bot_lk(i,j,nt) = tpl_T_r
    c_t_lk  (i,j,nt) = C_T_min
    h_ml_lk (i,j,nt) = 0.0_wp
    ! If sea ice scheme is not used,
    ! "t_ice" and "h_ice" are only used if "llake" is set true.
    ! Initial values of "t_ice" and "h_ice" are set through the surface analysis.
    IF(.NOT.lseaice) THEN
      h_ice (i,j,nt) = h_ice (i,j,nx)
      t_ice (i,j,nt) = t_ice (i,j,nx)
    ENDIF
  ENDIF Lake_points
! Variables are not part of the input (not read from GRIB files),
! set them to their reference values at all points.
  t_b1_lk(i,j,nt) = tpl_T_r          ! Reference value, bottom-sediment module
                                     !     is switched off
  h_b1_lk(i,j,nt) = rflk_depth_bs_ref! Reference value, bottom-sediment module
                                     !     is switched off
END DO
END DO
END DO

! If multi-layer snow model is used, set "t_snow_mult" at lake points
! to the ice surface temperature at all vertical levels.
IF (lmulti_snow) THEN
  DO nt = 1, nztlev   ! DO loop over time levels
  DO ks = 0, ke_snow  ! DO loop over snow layers
    WHERE (depth_lk(:,:) > 0.0_wp) t_snow_mult(:,:,ks,nt) = t_snow(:,:,nx)
  END DO
  END DO
ENDIF

! If multi-layer soil model is used, set "t_so(:,:,:,:)" at lake points
! to the fresh-water temperature of maximum density at all levels except 0.
! Notice that the climatological soil temperature,
! i.e. t_so(:,:,ke_soil+1,:), is left intact.
IF (lmulti_layer) THEN
  DO nt = 1, nztlev   ! DO loop over time levels
  DO ks = 1, ke_soil  ! DO loop over soil layers
    WHERE (depth_lk(:,:) > 0.0_wp) t_so(:,:,ks,nt) = tpl_T_r
  END DO
  END DO
ENDIF

!_dbg
  IF (izdebug > 5) THEN
    WRITE(*,*) 'FLake: WARM START initialization completed. '
  ENDIF
!_dbg

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END SUBROUTINE flake_init

!------------------------------------------------------------------------------
!  End of lake model initialization
!==============================================================================

!==============================================================================
!  FLake interface. Advances the lake model variables one time step ahead.
!------------------------------------------------------------------------------

SUBROUTINE flake_interface 

!------------------------------------------------------------------------------
!
! Description:
!
!  The FLake interface is a communication routine between "flake_driver" and 
!  COSMO (or any other model that uses FLake).
!  It assigns the FLake variables at the previous time step 
!  to their input values given by COSMO,
!  calls a number of routines to compute the heat and radiation fluxes,
!  calls "flake_driver",
!  and returns the updated FLake variables to COSMO.
!  The "flake_interface" does not contain any FLake physics. 
!  It only serves as a convenient means to organize calls of "flake_driver"
!  and of external routines that compute heat and radiation fluxes.
!  These calls are executed in a DO loop over all horizontal grid-points.
!  For sea and land grid-points, no action is taken so that
!  FLake variables remain equal to their reference values. 
!
!==============================================================================

! Declarations

!  Local variables of type INTEGER

INTEGER (KIND = iintegers) :: &
  i, j                      , & ! Loop indices
  ksnow                     , & ! Loop index (snow layers)
  im1, jm1                  , & ! i-1, j-1
  nx                        , & ! Time-level for integration (nold or nnow)
  izerror                   , & ! Error number
  izdebug                   , & ! for debug output
  istarts                   , & ! Start index for x-direction
  iends                     , & ! End   index for x-direction
  jstarts                   , & ! Start index for y-direction
  jends                         ! End   index for y-directin

!  Local variables of type REAL

REAL (KIND = wp)     :: &
  del_time            , & ! The time step used by FLake [s]
  depth_w             , & ! The lake depth [m]
  fetch               , & ! Typical wind fetch [m]
  depth_bs            , & ! Depth of the thermally active layer of the bottom 
                          ! sediments [m]
  T_bs                , & ! Temperature at the outer edge of 
                          ! the thermally active layer of sediments [K]
  par_Coriolis        , & ! The Coriolis parameter [s^{-1}]
  U_a_in              , & ! Wind speed at the first COSMO main level above the 
                          ! surface [m s^{-1}]
  T_a_in              , & ! Air temperature at the first COSMO main level above 
                          ! the surface [K]
  q_a_in              , & ! Air specific humidity at the first COSMO main level 
                          ! above the surface [-]
  q_s_in              , & ! Specific humidity at the surface [-]
  P_a_in              , & ! Surface air pressure [N m^{-2} = kg m^{-1} s^{-2}]
  rho_a               , & ! Air density [kg m^{-3}]
  T_sfc_p             , & ! Surface temperature at the previous time step [K]  
  T_sfc_n             , & ! Updated surface temperature [K]  
  Q_momentum          , & ! Momentum flux [N m^{-2}]
  Q_sensible          , & ! Sensible heat flux [W m^{-2}]
  Q_latent                ! Latent heat flux [W m^{-2}]

REAL (KIND = wp)     :: &
  albedo_water        , & ! Water surface albedo with respect to the short-wave
                          ! radiation
  albedo_ice          , & ! Ice surface albedo with respect to the short-wave
                          ! radiation
  albedo_snow             ! Snow surface albedo with respect to the short-wave
                          ! radiation

CHARACTER (LEN=80)   :: yzerrmsg

!  Local variables of derived type OPTICPAR_MEDIUM

TYPE (opticpar_medium) :: & 
  opticpar_water        , & ! Optical characteristics of water
  opticpar_ice          , & ! Optical characteristics of ice
  opticpar_snow             ! Optical characteristics of snow 

!  Local parameters of type REAL 

REAL (KIND = wp),     PARAMETER :: &
  Vel_min = 1.0E-02_wp       ! Minimum wind speed [m s^{-1}] 
                             ! Used by MR to compute transfer coefficient for 
                             ! momentum 
                             
!  Tracer pointers

REAL (KIND=wp),     POINTER :: &
  qv  (:,:,:) => NULL()      ! QV at tlev=nx

CHARACTER(LEN=25) :: yzroutine = 'flake_interface'

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

! Set method of debug output
IF (ldebug_soi) THEN
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

! Set start and end indices for DO loops in the horizontal directions
istarts = istartpar
iends   = iendpar
jstarts = jstartpar
jends   = jendpar

! Set time step and the previous time level for the two time level scheme
IF ((ntstep==0).AND.(.NOT.l2tls)) THEN
  ! For three time level scheme, use the original dt and not dt/2
  del_time = 2._wp*dt
ELSE
  del_time = dt
END IF
nx = nnow

!------------------------------------------------------------------------------
! Retrieve pointer to required tracers
!------------------------------------------------------------------------------

CALL trcr_get(izerror, idt_qv, ptr_tlev = nx, ptr = qv)
IF (izerror /= 0) THEN
  yzerrmsg = trcr_errorstr(izerror)
  CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
ENDIF

!------------------------------------------------------------------------------
!  Set optical characteristics of the lake water, lake ice and snow
!------------------------------------------------------------------------------

! Use default values
opticpar_water = opticpar_water_ref    ! Reference value, 95% of radiation is 
                                       !    absorbed within the upper 1 m
opticpar_ice   = opticpar_ice_opaque   ! Opaque ice
opticpar_snow  = opticpar_snow_opaque  ! Opaque snow

!------------------------------------------------------------------------------
!  Set/compute albedo of the lake water
!------------------------------------------------------------------------------

! Use default value
!_dm The net solar radiation flux from COSMO already accounts for the surface albedo.
!_nu albedo_water = albedo_water_ref
albedo_water = 0._wp
!_dm In case "albedo_water_ref" is used, 
!_dm changes should be made in COSMO where appropriate.
!_dm Otherwise, the radiation fluxes used by FLake become inconsistent 
!_dm with the COSMO radiation calculations.

!------------------------------------------------------------------------------
!  Compute/set albedos of the lake ice and snow
!------------------------------------------------------------------------------

! Use interpolation formula similar to that in the operational GME (c/o DM&BR)
!_dm The net solar radiation flux from COSMO already accounts for the ice surface albedo.
!_dm Important! Check if the sea/lake ice albedo is handled consistently
!_dm in the COSMO radiation routines. 
!_dm Introduce MR's approximation for the ice albedo as needed.
!_nu albedo_ice   = EXP(-c_albice_MR*(tpl_T_f-T_sfc_p)/tpl_T_f)
!_nu albedo_ice   = albedo_whiteice_ref*(1._wp-albedo_ice) + 
!_nu                                            albedo_blueice_ref*albedo_ice
    albedo_ice   = 0._wp
    albedo_snow  = albedo_ice     ! Snow is not considered
!_dm In case "albedo_ice" and "albedo_snow" defined in FLake are used, 
!_dm changes should be made in COSMO where appropriate.
!_dm Otherwise, the radiation fluxes used by FLake become inconsistent 
!_dm with the COSMO radiation calculations.

!------------------------------------------------------------------------------
!  DO loops over COSMO horizontal grid points
!------------------------------------------------------------------------------

DO j = jstarts, jends
!CDIR NOMOVE
DO i = istarts, iends

  IF(depth_lk(i,j) > 0.0_wp) THEN   
    ! Lake depth is positive, advance FLake variables. 
    ! Otherwise, take no action.

    depth_w      = depth_lk(i,j)
    fetch        = fetch_lk(i,j)      
    depth_bs     = dp_bs_lk(i,j)
    T_bs         = t_bs_lk (i,j)
    par_Coriolis = fc      (i,j)

!------------------------------------------------------------------------------
!  Set initial values
!------------------------------------------------------------------------------

! "t_snow" is always defined independent of the use of
! multi-layer snow model
    T_snow_p_flk = t_snow  (i,j,nx)
    T_ice_p_flk  = t_ice   (i,j,nx)  
    T_mnw_p_flk  = t_mnw_lk(i,j,nx) 
    T_wML_p_flk  = t_wml_lk(i,j,nx)
    T_bot_p_flk  = t_bot_lk(i,j,nx)
    T_B1_p_flk   = t_b1_lk (i,j,nx)
    C_T_p_flk    = c_t_lk  (i,j,nx)
    h_snow_p_flk = h_snow  (i,j,nx)
    h_ice_p_flk  = h_ice   (i,j,nx)
    h_ML_p_flk   = h_ml_lk (i,j,nx)
    H_B1_p_flk   = h_b1_lk (i,j,nx)
    T_sfc_p      = t_s     (i,j,nx)    ! Use COSMO variable "t_s"

!------------------------------------------------------------------------------
!  Set the rate of snow accumulation
!------------------------------------------------------------------------------

    dMsnowdt_flk = 0._wp      ! Snow is not considered

!------------------------------------------------------------------------------
!  Compute short-wave radiation fluxes (positive downward)
!------------------------------------------------------------------------------

    I_atm_flk = sobs(i,j)        ! Net solar radiation flux from COSMO
    CALL flake_radflux ( depth_w, albedo_water, albedo_ice, albedo_snow, &
                         opticpar_water, opticpar_ice, opticpar_snow )

!------------------------------------------------------------------------------
!  Compute long-wave radiation fluxes (positive downward)
!------------------------------------------------------------------------------

    Q_w_flk = thbs(i,j)          ! Net long-wave radiation flux from COSMO
                                 ! "thbs" already includes the upward 
                                 ! radiation from the surface

!------------------------------------------------------------------------------
!  Compute the surface friction velocity and fluxes of sensible and latent heat 
!------------------------------------------------------------------------------

!_dm DM's flux calculation routines are not used.
!_dm The fluxes of momentum and of sensible and latent heat are computed
!_dm using transfer coefficients from the COSMO turbulence routines (c/o MR). 

    im1  = MAX(1,i-1)             ! Compute wind speed at the first COSMO main 
                                  ! level above the surface
    jm1  = MAX(1,j-1)             ! as in the soil model 
    U_a_in       =  0.5_wp * SQRT                                &
                    ( (u(i,j,ke,nx) + u(im1,j,ke,nx))**2_iintegers   &
                     +(v(i,j,ke,nx) + v(i,jm1,ke,nx))**2_iintegers )

    q_a_in       = qv(i,j,ke)     ! Specific humidity at the first COSMO main 
                                  ! level above the surface
    q_s_in       = qv_s(i,j,nx)   ! Specific humidity at the surface 
    T_a_in       = t(i,j,ke,nx)*(ps(i,j,nx)/(p0(i,j,ke)+pp(i,j,ke,nx)))**rdocp 
                                  ! T_a_in is the potential temperature at the 
                                  ! first COSMO main level
                                  ! The surface temperature is the potential 
                                  ! temperature by definition
    P_a_in       = ps(i,j,nx)     ! Surface pressure
    rho_a = SfcFlx_rhoair(T_sfc_p, q_s_in, P_a_in) ! Air density at the surface 

    Q_momentum   = MAX(Vel_min, U_a_in)            ! Wind speed scale used by MR
    Q_sensible   = -rho_a*cp_d*tch(i,j)*Q_momentum*(T_a_in-T_sfc_p)   
    Q_latent     = lh_v
    IF(h_ice_p_flk.GE.h_Ice_min_flk) Q_latent = Q_latent + lh_f   
                                  ! Add latent heat of fusion over ice
    Q_latent     = -rho_a*Q_latent*tch(i,j)*Q_momentum*(q_a_in-q_s_in)
    Q_momentum   = -rho_a*tcm(i,j)*Q_momentum*U_a_in
    u_star_w_flk = SQRT(-Q_momentum/tpl_rho_w_r)

!------------------------------------------------------------------------------
!  Compute heat fluxes Q_snow_flk, Q_ice_flk, Q_w_flk
!------------------------------------------------------------------------------

    Q_w_flk = Q_w_flk - Q_sensible - Q_latent  
                     ! Add sensible and latent heat fluxes (notice the signs)

    IF(h_ice_p_flk.GE.h_Ice_min_flk) THEN        ! Ice exists
      IF(h_snow_p_flk.GE.h_Snow_min_flk) THEN    ! There is snow above the ice
        Q_snow_flk = Q_w_flk
        Q_ice_flk  = 0._wp
        Q_w_flk    = 0._wp
      ELSE                                       ! No snow above the ice
        Q_snow_flk = 0._wp
        Q_ice_flk  = Q_w_flk
        Q_w_flk    = 0._wp
      END IF
    ELSE                                         ! No ice-snow cover
        Q_snow_flk = 0._wp
        Q_ice_flk  = 0._wp
    END IF

!------------------------------------------------------------------------------
!  Advance FLake variables
!------------------------------------------------------------------------------

    CALL flake_driver ( depth_w, depth_bs, T_bs, par_Coriolis,         &
                        opticpar_water%extincoef_optic(1),             &
                        del_time, T_sfc_p, T_sfc_n, izdebug )

!------------------------------------------------------------------------------
!  Set output values
!------------------------------------------------------------------------------

! FLake variables
! "t_snow" is always defined independent of the use of
! multi-layer snow model
    t_snow  (i,j,nnew) = T_snow_n_flk
    t_ice   (i,j,nnew) = T_ice_n_flk      
    t_mnw_lk(i,j,nnew) = T_mnw_n_flk     
    t_wml_lk(i,j,nnew) = T_wML_n_flk    
    t_bot_lk(i,j,nnew) = T_bot_n_flk   
    t_b1_lk (i,j,nnew) = T_B1_n_flk    
    c_t_lk  (i,j,nnew) = C_T_n_flk     
    h_snow  (i,j,nnew) = h_snow_n_flk   
    h_ice   (i,j,nnew) = h_ice_n_flk    
    h_ml_lk (i,j,nnew) = h_ML_n_flk     
    h_b1_lk (i,j,nnew) = H_B1_n_flk    
    t_s     (i,j,nnew) = T_sfc_n       

! Weighted surface temperature of the COSMO grid point
! is equal to the lake surface temperature
    t_g     (i,j,nnew) = t_s (i,j,nnew)

! Save lake surface temperature in t_so(:,:,0,:)
! in case the multi-layer soil model is used.
    IF (lmulti_layer) THEN
      t_so(i,j,0,nnew) = t_s (i,j,nnew)
    ENDIF

  ENDIF    

END DO          
END DO             ! End of DO loops over COSMO horizontal grid points

IF (lmulti_snow) THEN
  ! If multi-layer snow model is used, save updated snow temperature
  ! also in "t_snow_mult" at all vertical levels
  DO ksnow = 0, ke_snow  ! DO loop over snow layers
    WHERE (depth_lk(:,:) > 0.0_wp) t_snow_mult(:,:,ksnow,nnew) = t_snow(:,:,nnew)
  END DO
ENDIF

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END SUBROUTINE flake_interface

!------------------------------------------------------------------------------
!  End of FLake interface
!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================
!------------------------------------------------------------------------------

SUBROUTINE flake_radflux ( depth_w, albedo_water, albedo_ice, albedo_snow, & 
                           opticpar_water, opticpar_ice, opticpar_snow )       

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the radiation fluxes 
!  at the snow-ice, ice-water, air-water, 
!  mixed layer-thermocline and water column-bottom sediment interfaces,
!  the mean radiation flux over the mixed layer,
!  and the mean radiation flux over the thermocline.
!
!------------------------------------------------------------------------------

! Declarations

!  Input (procedure arguments)

REAL (KIND = wp),     INTENT(IN) ::   &
  depth_w                           , & ! The lake depth [m]
  albedo_water                      , & ! Albedo of the water surface 
  albedo_ice                        , & ! Albedo of the ice surface
  albedo_snow                           ! Albedo of the snow surface

TYPE (opticpar_medium), INTENT(IN) :: & 
  opticpar_water                    , & ! Optical characteristics of water
  opticpar_ice                      , & ! Optical characteristics of ice
  opticpar_snow                         ! Optical characteristics of snow 


!  Local variables of type INTEGER
!_nu INTEGER (KIND = iintegers) :: & ! Help variable(s)
!_nu   i                             ! DO loop index

!==============================================================================

!     An n-band aproximation of the exponential decay law for the solar 
!     radiation flux, where n can be any value from 1 to 10, 
!     is not used due to bad vectorization. 
!     In the new version of "flake_radflux",
!     the number of wave-length bands is restricted to two and 
!     the summation over bands is perfomed explicitly (no DO loops).

!  Original version (!_nu = unused) but should remain in the file for a while

!_nu !==============================================================================
!_nu !  Start calculations
!_nu !------------------------------------------------------------------------------
!_nu
!_nu   IF(h_ice_p_flk >= h_Ice_min_flk) THEN            
!_nu     ! Ice exists
!_nu     IF(h_snow_p_flk >= h_Snow_min_flk) THEN        
!_nu       ! There is snow above the ice
!_nu       I_snow_flk = I_atm_flk*(1._wp-albedo_snow) 
!_nu       I_bot_flk = 0._wp
!_nu !CDIR EXPAND=opticpar_water%nband_optic
!_nu       DO i=1, opticpar_snow%nband_optic
!_nu         I_bot_flk = I_bot_flk +                    & 
!_nu         opticpar_snow%frac_optic(i)*EXP(-opticpar_snow%extincoef_optic(i)  &
!_nu                                    *h_snow_p_flk) 
!_nu       END DO 
!_nu       I_ice_flk  = I_snow_flk*I_bot_flk
!_nu     ELSE                                           
!_nu       ! No snow above the ice 
!_nu       I_snow_flk = I_atm_flk  
!_nu       I_ice_flk  = I_atm_flk*(1._wp-albedo_ice)
!_nu     END IF 
!_nu     I_bot_flk = 0._wp
!_nu     DO i=1, opticpar_ice%nband_optic
!_nu       I_bot_flk = I_bot_flk +                      & 
!_nu       opticpar_ice%frac_optic(i)*EXP(-opticpar_ice%extincoef_optic(i)      &
!_nu                                 *h_ice_p_flk) 
!_nu     END DO 
!_nu     I_w_flk      = I_ice_flk*I_bot_flk
!_nu   ELSE                                             
!_nu     ! No ice-snow cover
!_nu     I_snow_flk   = I_atm_flk  
!_nu     I_ice_flk    = I_atm_flk
!_nu     I_w_flk      = I_atm_flk*(1._wp-albedo_water)
!_nu   END IF 
!_nu 
!_nu   IF(h_ML_p_flk >= h_ML_min_flk) THEN
!_nu     ! Radiation flux at the bottom of the mixed layer
!_nu     I_bot_flk = 0._wp
!_nu     DO i=1, opticpar_water%nband_optic
!_nu       I_bot_flk = I_bot_flk +            & 
!_nu       opticpar_water%frac_optic(i)*EXP(-opticpar_water%extincoef_optic(i)  &
!_nu                                   *h_ML_p_flk) 
!_nu     END DO 
!_nu     I_h_flk = I_w_flk*I_bot_flk
!_nu   ELSE
!_nu     ! Mixed-layer depth is less then a minimum value
!_nu     I_h_flk = I_w_flk
!_nu   END IF
!_nu 
!_nu   ! Radiation flux at the lake bottom
!_nu   I_bot_flk = 0._wp
!_nu   DO i=1, opticpar_water%nband_optic
!_nu     I_bot_flk = I_bot_flk +              & 
!_nu     opticpar_water%frac_optic(i)*EXP(-opticpar_water%extincoef_optic(i)    &
!_nu                                 *depth_w) 
!_nu   END DO 
!_nu   I_bot_flk = I_w_flk*I_bot_flk
!_nu 
!_nu   IF(h_ML_p_flk.GE.h_ML_min_flk) THEN
!_nu     ! Integral-mean radiation flux over the mixed layer
!_nu     I_intm_0_h_flk = 0._wp
!_nu     DO i=1, opticpar_water%nband_optic
!_nu       I_intm_0_h_flk = I_intm_0_h_flk +                                    &
!_nu       opticpar_water%frac_optic(i)/opticpar_water%extincoef_optic(i)*      &
!_nu       (1._wp - EXP(-opticpar_water%extincoef_optic(i)*h_ML_p_flk))
!_nu     END DO 
!_nu     I_intm_0_h_flk = I_w_flk*I_intm_0_h_flk/h_ML_p_flk
!_nu   ELSE
!_nu     I_intm_0_h_flk = I_h_flk
!_nu   END IF
!_nu 
!_nu   IF(h_ML_p_flk.LE.depth_w-h_ML_min_flk) THEN
!_nu     ! Integral-mean radiation flux over the thermocline
!_nu     I_intm_h_D_flk = 0._wp 
!_nu     DO i=1, opticpar_water%nband_optic
!_nu       I_intm_h_D_flk = I_intm_h_D_flk +                                    &
!_nu       opticpar_water%frac_optic(i)/opticpar_water%extincoef_optic(i)*      &
!_nu       ( EXP(-opticpar_water%extincoef_optic(i)*h_ML_p_flk)                 &
!_nu       - EXP(-opticpar_water%extincoef_optic(i)*depth_w) )
!_nu     END DO 
!_nu     I_intm_h_D_flk = I_w_flk*I_intm_h_D_flk/(depth_w-h_ML_p_flk)
!_nu   ELSE
!_nu     I_intm_h_D_flk = I_h_flk
!_nu   END IF
!_nu 
!_nu !------------------------------------------------------------------------------
!_nu !  End calculations
!_nu !------------------------------------------------------------------------------

!  New version

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

  IF(h_ice_p_flk >= h_Ice_min_flk) THEN            
    ! Ice exists
    IF(h_snow_p_flk >= h_Snow_min_flk) THEN        
      ! There is snow above the ice
      I_snow_flk = I_atm_flk*(1._wp-albedo_snow) 
      I_bot_flk =                                                             &
      opticpar_snow%frac_optic(1)                                             &
      *EXP(-MIN(opticpar_snow%extincoef_optic(1)*h_snow_p_flk, c_maxearg_flk)) +  &
      opticpar_snow%frac_optic(2)                                             &
      *EXP(-MIN(opticpar_snow%extincoef_optic(2)*h_snow_p_flk, c_maxearg_flk))
      I_ice_flk  = I_snow_flk*I_bot_flk
    ELSE                                           
      ! No snow above the ice 
      I_snow_flk = I_atm_flk  
      I_ice_flk  = I_atm_flk*(1._wp-albedo_ice)
    END IF 
    I_bot_flk =                                                               &
    opticpar_ice%frac_optic(1)                                                &
    *EXP(-MIN(opticpar_ice%extincoef_optic(1)*h_ice_p_flk, c_maxearg_flk)) +  &
    opticpar_ice%frac_optic(2)                                                &
    *EXP(-MIN(opticpar_ice%extincoef_optic(2)*h_ice_p_flk, c_maxearg_flk))
    I_w_flk      = I_ice_flk*I_bot_flk
  ELSE                                             
    ! No ice-snow cover
    I_snow_flk   = I_atm_flk  
    I_ice_flk    = I_atm_flk
    I_w_flk      = I_atm_flk*(1._wp-albedo_water)
  END IF 

  IF(h_ML_p_flk >= h_ML_min_flk) THEN
    ! Radiation flux at the bottom of the mixed layer
    I_bot_flk =                                                                &
    opticpar_water%frac_optic(1)                                               &
    *EXP(-MIN(opticpar_water%extincoef_optic(1)*h_ML_p_flk, c_maxearg_flk)) +  &
    opticpar_water%frac_optic(2)                                               &
    *EXP(-MIN(opticpar_water%extincoef_optic(2)*h_ML_p_flk, c_maxearg_flk))
    I_h_flk = I_w_flk*I_bot_flk
  ELSE
    ! Mixed-layer depth is less then a minimum value
    I_h_flk = I_w_flk
  END IF

  ! Radiation flux at the lake bottom
  I_bot_flk =                                                             &
  opticpar_water%frac_optic(1)                                            &
  *EXP(-MIN(opticpar_water%extincoef_optic(1)*depth_w, c_maxearg_flk)) +  &
  opticpar_water%frac_optic(2)                                            &
  *EXP(-MIN(opticpar_water%extincoef_optic(2)*depth_w, c_maxearg_flk))
  I_bot_flk = I_w_flk*I_bot_flk

  IF(h_ML_p_flk >= h_ML_min_flk) THEN
    ! Integral-mean radiation flux over the mixed layer
    I_intm_0_h_flk =                                                                        &
    opticpar_water%frac_optic(1)/opticpar_water%extincoef_optic(1)*                         &
    (1._wp - EXP(-MIN(opticpar_water%extincoef_optic(1)*h_ML_p_flk, c_maxearg_flk))) +  &
    opticpar_water%frac_optic(2)/opticpar_water%extincoef_optic(2)*                         &
    (1._wp - EXP(-MIN(opticpar_water%extincoef_optic(2)*h_ML_p_flk, c_maxearg_flk)))
    I_intm_0_h_flk = I_w_flk*I_intm_0_h_flk/h_ML_p_flk
  ELSE
    I_intm_0_h_flk = I_h_flk
  END IF

  IF(h_ML_p_flk <= depth_w-h_ML_min_flk) THEN
    ! Integral-mean radiation flux over the thermocline
    I_intm_h_D_flk =                                                           &
    opticpar_water%frac_optic(1)/opticpar_water%extincoef_optic(1)*            &
    ( EXP(-MAX(opticpar_water%extincoef_optic(1)*h_ML_p_flk, c_maxearg_flk))   &
    - EXP(-MAX(opticpar_water%extincoef_optic(1)*depth_w, c_maxearg_flk)) ) +  &
    opticpar_water%frac_optic(2)/opticpar_water%extincoef_optic(2)*            &
    ( EXP(-MAX(opticpar_water%extincoef_optic(2)*h_ML_p_flk, c_maxearg_flk))   &
    - EXP(-MAX(opticpar_water%extincoef_optic(2)*depth_w, c_maxearg_flk)) )
    I_intm_h_D_flk = I_w_flk*I_intm_h_D_flk/(depth_w-h_ML_p_flk)
  ELSE
    I_intm_h_D_flk = I_h_flk
  END IF

!------------------------------------------------------------------------------
!  End calculations
!------------------------------------------------------------------------------

END SUBROUTINE flake_radflux

!==============================================================================
!==============================================================================
!------------------------------------------------------------------------------

SUBROUTINE flake_driver ( depth_w, depth_bs, T_bs, par_Coriolis,       &
                          extincoef_water_typ,                         &
                          del_time, T_sfc_p, T_sfc_n, idbg )

!------------------------------------------------------------------------------
!
! Description:
!
!  The main driving routine of the lake model FLake 
!  where computations are performed.
!  Advances the surface temperature
!  and other FLake variables one time step.
!  At the moment, the Euler explicit scheme is used.
!
!------------------------------------------------------------------------------

!  Input (procedure arguments)

REAL (KIND = wp),     INTENT(IN) ::   &
  depth_w                           , & ! The lake depth [m]
  depth_bs                          , & ! Depth of the thermally active layer 
                                        ! of bottom sediments [m]
  T_bs                              , & ! Temperature at the outer edge of 
                                        ! the thermally active layer of bottom 
                                        ! sediments [K]
  par_Coriolis                      , & ! The Coriolis parameter [s^{-1}]
  extincoef_water_typ               , & ! "Typical" extinction coefficient of 
                                        ! the lake water [m^{-1}], used to 
                                        ! compute the equilibrium CBL depth
  del_time                          , & ! The model time step [s]
  T_sfc_p                               ! Surface temperature at the previous 
                                        ! time step [K]  
                                        ! (equal to either T_ice, T_snow 
                                        !  or T_wML)

!  Output (procedure arguments)

REAL (KIND = wp),     INTENT(OUT) ::  &
  T_sfc_n                               ! Updated surface temperature [K] 
                                        ! (equal to the updated value of 
                                        !  either T_ice, T_snow or T_wML)

INTEGER (KIND=iintegers), INTENT(IN) :: idbg   ! for debug output

!  Local variables of type LOGICAL
LOGICAL ::          &
  l_ice_create    , & ! .TRUE. = ice does not exist but should be created
  l_snow_exists   , & ! .TRUE. = there is snow above the ice
  l_ice_meltabove     ! .TRUE. = snow/ice melting from above takes place

!  Local variables of type INTEGER
!_nu INTEGER (KIND = iintegers) :: &
!_nu   i                             ! Loop index

!  Local variables of type REAL
REAL (KIND = wp)     ::    &
  d_T_mnw_dt             , & ! Time derivative of T_mnw [K s^{-1}] 
  d_T_ice_dt             , & ! Time derivative of T_ice [K s^{-1}] 
  d_T_bot_dt             , & ! Time derivative of T_bot [K s^{-1}] 
  d_T_B1_dt              , & ! Time derivative of T_B1 [K s^{-1}] 
  d_h_snow_dt            , & ! Time derivative of h_snow [m s^{-1}]
  d_h_ice_dt             , & ! Time derivative of h_ice [m s^{-1}]
  d_h_ML_dt              , & ! Time derivative of h_ML [m s^{-1}]
  d_H_B1_dt              , & ! Time derivative of H_B1 [m s^{-1}]
  d_C_T_dt                   ! Time derivative of C_T [s^{-1}]

!  Local variables of type REAL
REAL (KIND = wp)     :: &
  N_T_mean        , & ! The mean buoyancy frequency in the thermocline [s^{-1}] 
  ZM_h_scale      , & ! The ZM96 equilibrium SBL depth scale [m] 
  conv_equil_h_scale  ! The equilibrium CBL depth scale [m]

!  Local variables of type REAL
REAL (KIND = wp)     :: &
  h_ice_threshold, & ! If h_ice<h_ice_threshold, use quasi-equilibrium ice model 
  flk_str_1      , & ! Help storage variable
  flk_str_2      , & ! Help storage variable
  R_H_icesnow    , & ! Dimensionless ratio, used to store intermediate results
  R_rho_c_icesnow, & ! Dimensionless ratio, used to store intermediate results
  R_TI_icesnow   , & ! Dimensionless ratio, used to store intermediate results
  R_Tstar_icesnow    ! Dimensionless ratio, used to store intermediate results

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

!_dm 
! Security. Set time-rate-of-change of prognostic variables to zero.
! Set prognostic variables to their values at the previous time step.
! (This is to avoid spurious changes of prognostic variables 
! when FLake is used within a 3D model, e.g. to avoid spurious generation of ice 
! at the neighbouring lake points as noticed by Burkhardt Rockel.)
!_dm 

d_T_mnw_dt   = 0._wp 
d_T_ice_dt   = 0._wp 
d_T_bot_dt   = 0._wp 
d_T_B1_dt    = 0._wp 
d_h_snow_dt  = 0._wp 
d_h_ice_dt   = 0._wp 
d_h_ML_dt    = 0._wp 
d_H_B1_dt    = 0._wp 
d_C_T_dt     = 0._wp 
T_snow_n_flk = T_snow_p_flk   
T_ice_n_flk  = T_ice_p_flk    
T_wML_n_flk  = T_wML_p_flk   
T_mnw_n_flk  = T_mnw_p_flk     
T_bot_n_flk  = T_bot_p_flk  
T_B1_n_flk   = T_B1_p_flk      
h_snow_n_flk = h_snow_p_flk 
h_ice_n_flk  = h_ice_p_flk   
h_ML_n_flk   = h_ML_p_flk    
H_B1_n_flk   = H_B1_p_flk   
C_T_n_flk    = C_T_p_flk    

!------------------------------------------------------------------------------
!  Compute fluxes, using variables from the previous time step.
!------------------------------------------------------------------------------

!_dm
! At this point, the heat and radiation fluxes, namely,
! Q_snow_flk, Q_ice_flk, Q_w_flk, 
! I_atm_flk, I_snow_flk, I_ice_flk, I_w_flk, I_h_flk, I_bot_flk,     
! the mean radiation flux over the mixed layer, I_intm_0_h_flk, 
! and the mean radiation flux over the thermocline, I_intm_h_D_flk, 
! should be known.
! They are computed within "flake_interface" (or within the driving model)
! and are available to "flake_driver"
! through the above variables declared in the MODULE "flake".
! In case a lake is ice-covered, Q_w_flk is re-computed below.
!_dm

! Heat flux through the ice-water interface
IF(h_ice_p_flk >= h_Ice_min_flk) THEN    
  ! Ice exists 
  IF(h_ML_p_flk <= h_ML_min_flk) THEN    
    ! Mixed-layer depth is zero, compute flux 

    ! Flux with linear T(z) 
    Q_w_flk = -tpl_kappa_w*(T_bot_p_flk-T_wML_p_flk)/depth_w  

    ! d\Phi(0)/d\zeta (thermocline)
    Phi_T_pr0_flk = Phi_T_pr0_1*C_T_p_flk-Phi_T_pr0_2         

    ! Account for an increased d\Phi(0)/d\zeta 
    Q_w_flk = Q_w_flk*MAX(Phi_T_pr0_flk, 1._wp)           

  ELSE                    

    ! Mixed-layer depth is greater than zero, set flux to zero
    Q_w_flk = 0._wp                  

  END IF   
END IF   

! A generalized heat flux scale 
Q_star_flk = Q_w_flk + I_w_flk + I_h_flk - 2._wp*I_intm_0_h_flk

! Heat flux through the water-bottom sediment interface
IF(lflk_botsed_use) THEN
  Q_bot_flk = -tpl_kappa_w                                                   &
           *(T_B1_p_flk-T_bot_p_flk)/MAX(H_B1_p_flk, H_B1_min_flk)*Phi_B1_pr0
ELSE  
  Q_bot_flk = 0._wp   ! The bottom-sediment scheme is not used
END IF


!------------------------------------------------------------------------------
!  Check if ice exists or should be created.
!  If so, compute the thickness and the temperature of ice and snow.
!------------------------------------------------------------------------------

!_dm
! Notice that a quasi-equilibrium ice-snow model is used 
! to avoid numerical instability when the ice is thin.
! This is always the case when new ice is created.
!_dm

!_dev
! The dependence of snow density and of snow heat conductivity 
! on the snow thickness is accounted for parametrically.
! That is, the time derivatives of \rho_S and \kappa_S are neglected.
! The exception is the equation for the snow thickness 
! in case of snow accumulation and no melting, 
! where d\rho_S/dt is incorporated.
! Furthermore, some (presumably small) correction terms incorporating 
! the snow density and the snow heat conductivity are dropped out.
! Those terms may be included as better formulations 
! for \rho_S and \kappa_S are available.
!_dev

! Default values
l_ice_create    = .FALSE.  
l_ice_meltabove = .FALSE.  

Ice_exist: IF(h_ice_p_flk < h_Ice_min_flk) THEN   
  ! Ice does not exist 

  l_ice_create = (T_wML_p_flk <= (tpl_T_f+c_small_flk)) .AND.    &
                 (Q_w_flk     <   0._wp)

  IF(l_ice_create) THEN                            
    ! Ice does not exist but should be created
    d_h_ice_dt = -Q_w_flk/tpl_rho_I/tpl_L_f                                  

    ! Advance h_ice 
    h_ice_n_flk = h_ice_p_flk + d_h_ice_dt*del_time                          

    ! Ice temperature
    T_ice_n_flk = tpl_T_f + h_ice_n_flk*Q_w_flk/tpl_kappa_I/Phi_I_pr0_lin    
    d_h_snow_dt = dMsnowdt_flk/tpl_rho_S_min 

    ! Advance h_snow
    h_snow_n_flk = h_snow_p_flk + d_h_snow_dt*del_time                       

    ! d\Phi_I(1)/d\zeta_I (ice)
    Phi_I_pr1_flk = Phi_I_pr1_lin                                    & 
                  + Phi_I_ast_MR*MIN(1._wp, h_ice_n_flk/H_Ice_max)       

    R_H_icesnow = Phi_I_pr1_flk/Phi_S_pr0_lin                        &
                * tpl_kappa_I/flake_snowheatconduct(h_snow_n_flk)    &
                * h_snow_n_flk/MAX(h_ice_n_flk, h_Ice_min_flk)

    ! Snow temperature
    T_snow_n_flk = T_ice_n_flk + R_H_icesnow*(T_ice_n_flk-tpl_T_f)           
  END IF

ELSE Ice_exist                                     
  ! Ice exists

  ! Check if there is snow above the ice
  l_snow_exists = h_snow_p_flk >= h_Snow_min_flk   

  Melting: IF(T_snow_p_flk >= (tpl_T_f-c_small_flk)) THEN  
    ! T_sfc = T_f, check for melting from above
    ! T_snow = T_ice if snow is absent 

    IF (l_snow_exists) THEN   
      ! There is snow above the ice
      flk_str_1 = Q_snow_flk + I_snow_flk - I_ice_flk    ! Atmospheric forcing
      IF(flk_str_1 >= 0._wp) THEN  ! Melting of snow and ice from above
        l_ice_meltabove = .TRUE.
        d_h_snow_dt = (-flk_str_1/tpl_L_f+dMsnowdt_flk)/             &
                                         flake_snowdensity(h_snow_p_flk)
        d_h_ice_dt  = -(I_ice_flk - I_w_flk - Q_w_flk)/tpl_L_f/tpl_rho_I 
      END IF 
    ELSE                     
      ! No snow above the ice
      ! Atmospheric forcing + heating from the water
      flk_str_1 = Q_ice_flk + I_ice_flk - I_w_flk - Q_w_flk  

      IF(flk_str_1.GE.0._wp) THEN  
        ! Melting of ice from above, snow accumulation may occur
        l_ice_meltabove = .TRUE.
        d_h_ice_dt  = -flk_str_1/tpl_L_f/tpl_rho_I 
        d_h_snow_dt = dMsnowdt_flk/tpl_rho_S_min
      END IF 
    END IF 
    IF(l_ice_meltabove) THEN  ! Melting from above takes place
      h_ice_n_flk  = h_ice_p_flk  + d_h_ice_dt *del_time  ! Advance h_ice
      h_snow_n_flk = h_snow_p_flk + d_h_snow_dt*del_time  ! Advance h_snow
      T_ice_n_flk  = tpl_T_f            ! Set T_ice to the freezing point
      T_snow_n_flk = tpl_T_f            ! Set T_snow to the freezing point
    END IF

  END IF Melting

  No_Melting: IF(.NOT.l_ice_meltabove) THEN     ! No melting from above

    d_h_snow_dt = flake_snowdensity(h_snow_p_flk)  
    IF(d_h_snow_dt.LT.tpl_rho_S_max) THEN    ! Account for d\rho_S/dt
     flk_str_1 = h_snow_p_flk*tpl_Gamma_rho_S/tpl_rho_w_r
     flk_str_1 = flk_str_1/(1._wp-flk_str_1)
    ELSE                                     ! Snow density is equal to its 
                                             ! maximum value, d\rho_S/dt=0
     flk_str_1 = 0._wp
    END IF

    ! Snow accumulation
    d_h_snow_dt  = dMsnowdt_flk/d_h_snow_dt/(1._wp+flk_str_1)       

    ! Advance h_snow
    h_snow_n_flk = h_snow_p_flk + d_h_snow_dt*del_time                         
    
    ! h_ice relative to its maximum value
    Phi_I_pr0_flk = h_ice_p_flk/H_Ice_max                              

    ! Shape factor (ice)
    C_I_flk = C_I_lin - C_I_MR*(1._wp+Phi_I_ast_MR)*Phi_I_pr0_flk  

    ! d\Phi_I(1)/d\zeta_I (ice)
    Phi_I_pr1_flk = Phi_I_pr1_lin + Phi_I_ast_MR*Phi_I_pr0_flk         

    ! d\Phi_I(0)/d\zeta_I (ice)
    Phi_I_pr0_flk = Phi_I_pr0_lin - Phi_I_pr0_flk                      

    h_ice_threshold = MAX(1._wp, 2._wp*C_I_flk*tpl_c_I*              &
                                               (tpl_T_f-T_ice_p_flk)/tpl_L_f)
    h_ice_threshold = Phi_I_pr0_flk/C_I_flk*tpl_kappa_I/tpl_rho_I/tpl_c_I*   &
                                                              h_ice_threshold
    ! Threshold value of h_ice
    h_ice_threshold = SQRT(h_ice_threshold*del_time)                   

    ! h_ice(threshold) < 0.9*H_Ice_max
    h_ice_threshold = MIN(0.9_wp*H_Ice_max,                              &
                          MAX(h_ice_threshold, h_Ice_min_flk))

    IF(h_ice_p_flk < h_ice_threshold) THEN  
      ! Use a quasi-equilibrium ice model

      IF(l_snow_exists) THEN   ! Use fluxes at the air-snow interface
        flk_str_1 = Q_snow_flk + I_snow_flk - I_w_flk
      ELSE                     ! Use fluxes at the air-ice interface
        flk_str_1 = Q_ice_flk + I_ice_flk - I_w_flk
      END IF
      d_h_ice_dt = -(flk_str_1-Q_w_flk)/tpl_L_f/tpl_rho_I
      ! Advance h_ice
      h_ice_n_flk = h_ice_p_flk + d_h_ice_dt *del_time                         
      ! Ice temperature
      T_ice_n_flk = tpl_T_f + h_ice_n_flk*flk_str_1/tpl_kappa_I/Phi_I_pr0_flk

    ELSE                                     
      ! Use a complete ice model

      d_h_ice_dt = tpl_kappa_I*(tpl_T_f-T_ice_p_flk)/h_ice_p_flk*Phi_I_pr0_flk
      d_h_ice_dt = (Q_w_flk+d_h_ice_dt)/tpl_L_f/tpl_rho_I
      ! Advance h_ice
      h_ice_n_flk = h_ice_p_flk  + d_h_ice_dt*del_time

      ! Dimensionless parameters
      R_TI_icesnow = tpl_c_I*(tpl_T_f-T_ice_p_flk)/tpl_L_f
      R_Tstar_icesnow = 1._wp - C_I_flk

      IF(l_snow_exists) THEN  
        ! There is snow above the ice
        R_H_icesnow = Phi_I_pr1_flk/Phi_S_pr0_lin                         &
                        * tpl_kappa_I/flake_snowheatconduct(h_snow_p_flk) &
                        * h_snow_p_flk/h_ice_p_flk
        R_rho_c_icesnow = flake_snowdensity(h_snow_p_flk)*                &
                          tpl_c_S/tpl_rho_I/tpl_c_I 
!_dev 
!_dm 
! These terms should be included as an improved understanding of the snow 
! scheme is gained, of the effect of snow density in particular. 
!_dm 
!_nu        R_Tstar_icesnow = R_Tstar_icesnow                               &
!_nu                        + (1._wp+C_S_lin*h_snow_p_flk/h_ice_p_flk)  &
!_nu                        *  R_H_icesnow*R_rho_c_icesnow
!_dev

        ! Dimensionless parameter
        R_Tstar_icesnow = R_Tstar_icesnow*R_TI_icesnow   

!_dev
!_nu        R_Tstar_icesnow = R_Tstar_icesnow                               &
!_nu                        + (1._wp-R_rho_c_icesnow)*tpl_c_I *         &
!_nu                          T_ice_p_flk/tpl_L_f
!_dev
        ! Atmospheric fluxes
        flk_str_2 = Q_snow_flk+I_snow_flk-I_w_flk                  
        flk_str_1  = C_I_flk*h_ice_p_flk + (1._wp+C_S_lin*R_H_icesnow)  &
                            *R_rho_c_icesnow*h_snow_p_flk

        ! Effect of snow accumulation
        d_T_ice_dt = -(1._wp-2._wp*C_S_lin)*R_H_icesnow             &
                        *(tpl_T_f-T_ice_p_flk)                              &
                        * tpl_c_S*dMsnowdt_flk                          
      ELSE
        ! No snow above the ice

        ! Dimensionless parameter
        R_Tstar_icesnow = R_Tstar_icesnow*R_TI_icesnow  

        ! Atmospheric fluxes
        flk_str_2 = Q_ice_flk+I_ice_flk-I_w_flk                    
        flk_str_1  = C_I_flk*h_ice_p_flk
        d_T_ice_dt = 0._wp

      END IF 
      ! Add flux due to heat conduction
      d_T_ice_dt = d_T_ice_dt + tpl_kappa_I*(tpl_T_f-T_ice_p_flk)/h_ice_p_flk&
                    * Phi_I_pr0_flk * (1._wp-R_Tstar_icesnow)                     
      ! Add flux from water to ice
      d_T_ice_dt = d_T_ice_dt - R_Tstar_icesnow*Q_w_flk            

      ! Add atmospheric fluxes
      d_T_ice_dt = d_T_ice_dt + flk_str_2                          

      ! Total forcing
      d_T_ice_dt = d_T_ice_dt/tpl_rho_I/tpl_c_I                    

      ! dT_ice/dt 
      d_T_ice_dt = d_T_ice_dt/flk_str_1                            

      ! Advance T_ice
      T_ice_n_flk = T_ice_p_flk + d_T_ice_dt*del_time                          
    END IF

    ! h_ice relative to its maximum value
    Phi_I_pr1_flk = MIN(1._wp, h_ice_n_flk/H_Ice_max)          

    ! d\Phi_I(1)/d\zeta_I (ice)
    Phi_I_pr1_flk = Phi_I_pr1_lin + Phi_I_ast_MR*Phi_I_pr1_flk     
    R_H_icesnow = Phi_I_pr1_flk/Phi_S_pr0_lin                             &
                  * tpl_kappa_I/flake_snowheatconduct(h_snow_n_flk)       &
                  * h_snow_n_flk/MAX(h_ice_n_flk, h_Ice_min_flk)

    ! Snow temperature
    T_snow_n_flk = T_ice_n_flk + R_H_icesnow*(T_ice_n_flk-tpl_T_f)             

  END IF No_Melting

END IF Ice_exist   

! Security, limit h_ice by its maximum value
h_ice_n_flk = MIN(h_ice_n_flk, H_Ice_max)      

! Security, limit the ice and snow temperatures by the freezing point 
T_snow_n_flk = MIN(T_snow_n_flk, tpl_T_f)  
T_ice_n_flk =  MIN(T_ice_n_flk,  tpl_T_f)    

! Security, avoid too low values 
  T_snow_n_flk = MAX(T_snow_n_flk, 73.15_wp)  
  T_ice_n_flk =  MAX(T_ice_n_flk,  73.15_wp)    

! Remove too thin ice and/or snow
IF(h_ice_n_flk < h_Ice_min_flk)  THEN        ! Check ice
  h_ice_n_flk = 0._wp       ! Ice is too thin, remove it, and
  T_ice_n_flk = tpl_T_f         ! set T_ice to the freezing point.
  h_snow_n_flk = 0._wp      ! Remove snow when there is no ice, and
  T_snow_n_flk = tpl_T_f        ! set T_snow to the freezing point.
  l_ice_create = .FALSE.        ! "Exotic" case, ice has been created but 
                                ! proved to be too thin
ELSE IF(h_snow_n_flk < h_Snow_min_flk) THEN  ! Ice exists, check snow
  h_snow_n_flk = 0._wp      ! Snow is too thin, remove it, 
  T_snow_n_flk = T_ice_n_flk    ! and set the snow temperature equal to the 
                                ! ice temperature.
END IF

!------------------------------------------------------------------------------
!  Compute the mean temperature of the water column.
!------------------------------------------------------------------------------

IF (l_ice_create) THEN
  ! Ice has just been created, set Q_w to zero
  Q_w_flk = 0._wp     
ENDIF

d_T_mnw_dt = (Q_w_flk - Q_bot_flk + I_w_flk - I_bot_flk) /         &
                                              tpl_rho_w_r/tpl_c_w/depth_w
! Advance T_mnw
T_mnw_n_flk = T_mnw_p_flk + d_T_mnw_dt*del_time   

! Limit T_mnw by the freezing point 
T_mnw_n_flk = MAX(T_mnw_n_flk, tpl_T_f)           

!------------------------------------------------------------------------------
!  Compute the mixed-layer depth, the mixed-layer temperature, 
!  the bottom temperature and the shape factor
!  with respect to the temperature profile in the thermocline. 
!  Different formulations are used, depending on the regime of mixing. 
!------------------------------------------------------------------------------

HTC_Water: IF (h_ice_n_flk >= h_Ice_min_flk) THEN    ! Ice exists

  ! Limit the mean temperature under the ice by T_r 
  T_mnw_n_flk = MIN(T_mnw_n_flk, tpl_T_r) 

  ! The mixed-layer temperature is equal to the freezing point 
  T_wML_n_flk = tpl_T_f                   

  IF(l_ice_create) THEN                  ! Ice has just been created 
    IF(h_ML_p_flk.GE.depth_w-h_ML_min_flk) THEN ! h_ML=D when ice is created 
      h_ML_n_flk = 0._wp             ! Set h_ML to zero 
      C_T_n_flk = C_T_min                ! Set C_T to its minimum value 
    ELSE                                 ! h_ML<D when ice is created 
      h_ML_n_flk = h_ML_p_flk            ! h_ML remains unchanged 
      C_T_n_flk = C_T_p_flk              ! C_T (thermocline) remains unchanged 
    END IF 
    T_bot_n_flk = T_wML_n_flk - (T_wML_n_flk-T_mnw_n_flk)/C_T_n_flk/        &
                                             (1._wp-h_ML_n_flk/depth_w)
                                         ! Update the bottom temperature 

  ELSE IF(T_bot_p_flk.LT.tpl_T_r) THEN   ! Ice exists and T_bot < T_r, 
                                         ! molecular heat transfer 
    h_ML_n_flk = h_ML_p_flk              ! h_ML remains unchanged 
    C_T_n_flk = C_T_p_flk                ! C_T (thermocline) remains unchanged 
    T_bot_n_flk = T_wML_n_flk - (T_wML_n_flk-T_mnw_n_flk)/C_T_n_flk/        &
                                             (1._wp-h_ML_n_flk/depth_w)
                                         ! Update the bottom temperature 

  ELSE                                   ! Ice exists and T_bot = T_r, 
                                         !convection due to bottom heating 
    T_bot_n_flk = tpl_T_r                ! T_bot is equal to the temperature 
                                         ! of maximum density 
    IF(h_ML_p_flk.GE.c_small_flk) THEN   ! h_ML > 0 
      C_T_n_flk = C_T_p_flk              ! C_T (thermocline) remains unchanged 
      h_ML_n_flk = depth_w*(1._wp-(T_wML_n_flk-T_mnw_n_flk)/            &
                                       (T_wML_n_flk-T_bot_n_flk)/C_T_n_flk)
      h_ML_n_flk = MAX(h_ML_n_flk, 0._wp)   ! Update the mixed-layer depth  
    ELSE                                 ! h_ML = 0 
      h_ML_n_flk = h_ML_p_flk            ! h_ML remains unchanged 
      C_T_n_flk = (T_wML_n_flk-T_mnw_n_flk)/(T_wML_n_flk-T_bot_n_flk) 
      C_T_n_flk = MIN(C_T_max, MAX(C_T_n_flk, C_T_min)) 
                                      ! Update the shape factor (thermocline)
    END IF 
  END IF 

  ! Security, limit the bottom temperature by T_r 
  T_bot_n_flk = MIN(T_bot_n_flk, tpl_T_r)    

ELSE HTC_Water                                      ! Open water

! Generalised buoyancy flux scale and convective velocity scale
  flk_str_1 = flake_buoypar(T_wML_p_flk)*Q_star_flk/tpl_rho_w_r/tpl_c_w                    
  IF(flk_str_1.LT.0._wp) THEN       
    ! Convection     
    w_star_sfc_flk = (-flk_str_1*h_ML_p_flk)**(1._wp/3._wp)
  ELSE 
    ! Neutral or stable stratification
    w_star_sfc_flk = 0._wp
  END IF 

!_dm
! The equilibrium depth of the CBL due to surface cooling with the volumetric
! heating is not computed as a solution to the transcendental equation.
! Instead, an algebraic formula is used
! that interpolates between the two asymptotic limits.
!_dm

  conv_equil_h_scale = -Q_w_flk/MAX(I_w_flk, c_small_flk)
  IF(conv_equil_h_scale.GT.0._wp .AND. conv_equil_h_scale.LT.1._wp  &
    .AND. T_wML_p_flk.GT.tpl_T_r) THEN   
    ! The equilibrium CBL depth scale is only used above T_r
    conv_equil_h_scale = SQRT(6._wp*conv_equil_h_scale)                 &
              + 2._wp*conv_equil_h_scale/(1._wp-conv_equil_h_scale)
    conv_equil_h_scale = MIN(depth_w, conv_equil_h_scale/extincoef_water_typ)
  ELSE
    ! Set the equilibrium CBL depth to zero
    conv_equil_h_scale = 0._wp
  END IF

! Mean buoyancy frequency in the thermocline
  N_T_mean = flake_buoypar(0.5_wp*(T_wML_p_flk+T_bot_p_flk)) *      &
                                      (T_wML_p_flk-T_bot_p_flk)
  IF(h_ML_p_flk.LE.depth_w-h_ML_min_flk) THEN
    N_T_mean = SQRT(N_T_mean/(depth_w-h_ML_p_flk))  ! Compute N                   
  ELSE 
    N_T_mean = 0._wp                            ! h_ML=D, set N to zero
  END IF 

! The rate of change of C_T
  d_C_T_dt = MAX(w_star_sfc_flk, u_star_w_flk, u_star_min_flk)**2_iintegers

  ! Relaxation time scale for C_T
  d_C_T_dt = N_T_mean*(depth_w-h_ML_p_flk)**2_iintegers       &
           / c_relax_C/d_C_T_dt                               

  ! Rate-of-change of C_T 
  d_C_T_dt = (C_T_max-C_T_min)/MAX(d_C_T_dt, c_small_flk)     

! Compute the shape factor and the mixed-layer depth, 
! using different formulations for convection and wind mixing

  ! C_TT, using C_T at the previous time step
  C_TT_flk = C_TT_1*C_T_p_flk-C_TT_2      

  ! C_Q using C_T at the previous time step
  C_Q_flk = 2._wp*C_TT_flk/C_T_p_flk  

  Mixing_regime: IF(flk_str_1.LT.0._wp) THEN  ! Convective mixing 

    ! Update C_T, assuming dh_ML/dt>0
    C_T_n_flk = C_T_p_flk + d_C_T_dt*del_time           

    ! Limit C_T 
    C_T_n_flk = MIN(C_T_max, MAX(C_T_n_flk, C_T_min))   

    ! Re-compute dC_T/dt
    d_C_T_dt = (C_T_n_flk-C_T_p_flk)/del_time           

    IF(h_ML_p_flk.LE.depth_w-h_ML_min_flk) THEN       ! Compute dh_ML/dt
      IF(h_ML_p_flk.LE.h_ML_min_flk) THEN    
        ! Use a reduced entrainment equation (spin-up)
        d_h_ML_dt = c_cbl_1/c_cbl_2*MAX(w_star_sfc_flk, c_small_flk)

!_dbg
!       IF (idbg > 10) THEN
!         PRINT *, ' FLake: reduced entrainment eq. D_time*d_h_ML_dt  = ', &
!                                                 d_h_ML_dt*del_time
!         PRINT *, '         w_*       = ', w_star_sfc_flk
!         PRINT *, '         \beta*Q_* = ', flk_str_1
!       ENDIF
!_dbg

      ELSE
        ! Use a complete entrainment equation 
        R_H_icesnow     = depth_w/h_ML_p_flk
        R_rho_c_icesnow = R_H_icesnow-1._wp
        R_TI_icesnow    = C_T_p_flk/C_TT_flk
        R_Tstar_icesnow = (R_TI_icesnow/2._wp-1._wp)*R_rho_c_icesnow &
                          + 1._wp
        d_h_ML_dt = -Q_star_flk*(R_Tstar_icesnow*(1._wp+c_cbl_1)-1._wp) &
                    - Q_bot_flk
        ! Q_* and Q_b flux terms
        d_h_ML_dt = d_h_ML_dt/tpl_rho_w_r/tpl_c_w

        flk_str_2 = (depth_w-h_ML_p_flk)*(T_wML_p_flk-T_bot_p_flk)*C_TT_2/   &
                                                              C_TT_flk*d_C_T_dt 

        ! Add dC_T/dt term
        d_h_ML_dt = d_h_ML_dt + flk_str_2

        flk_str_2 = I_bot_flk + (R_TI_icesnow-1._wp)*I_h_flk -           &
                                 R_TI_icesnow*I_intm_h_D_flk
        flk_str_2 = flk_str_2 + (R_TI_icesnow-2._wp)*R_rho_c_icesnow *   &
                                (I_h_flk-I_intm_0_h_flk)
        flk_str_2 = flk_str_2/tpl_rho_w_r/tpl_c_w

        ! Add radiation terms
        d_h_ML_dt = d_h_ML_dt + flk_str_2
        flk_str_2 = -c_cbl_2*R_Tstar_icesnow*Q_star_flk/tpl_rho_w_r/tpl_c_w/ &
                                   MAX(w_star_sfc_flk, c_small_flk)
        flk_str_2 = flk_str_2 + C_T_p_flk*(T_wML_p_flk-T_bot_p_flk)

        ! dh_ML/dt = r.h.s.
        d_h_ML_dt = d_h_ML_dt/flk_str_2
      END IF 

!_dm
! Notice that dh_ML/dt may appear to be negative  
! (e.g. due to buoyancy loss to bottom sediments and/or
! the effect of volumetric radiation heating),
! although a negative generalized buoyancy flux scale indicates 
! that the equilibrium CBL depth has not yet been reached
! and convective deepening of the mixed layer should take place. Physically, 
! this situation reflects an approximate character of the lake model.
! Using the self-similar temperature profile in the thermocline, 
! there is always communication between the mixed layer, the thermocline 
! and the lake bottom. As a result, the rate of change of the CBL depth
! is always dependent on the bottom heat flux and the radiation heating of the 
! thermocline. In reality, convective mixed-layer deepening may be completely 
! decoupled from the processes underneath. In order to account for this fact,
! the rate of CBL deepening is set to a small value
! if dh_ML/dt proves to be negative.
! This is "double insurance" however, 
! as a negative dh_ML/dt is encountered very rarely.
!_dm

!_dbg
!     IF (idbg > 10) THEN
!       IF(d_h_ML_dt.LT.0._wp) THEN 
!         PRINT *, 'FLake: negative d_h_ML_dt during convection, = ', d_h_ML_dt
!         PRINT *, '                d_h_ML_dt*del_time = ',        &
!                                          MAX(d_h_ML_dt, c_small_flk)*del_time
!         PRINT *, '         u_*       = ', u_star_w_flk   
!         PRINT *, '         w_*       = ', w_star_sfc_flk
!         PRINT *, '         h_CBL_eqi = ', conv_equil_h_scale
!         PRINT *, '         ZM scale  = ', ZM_h_scale
!         PRINT *, '        h_ML_p_flk = ', h_ML_p_flk
!       END IF
!       PRINT *, 'FLake: Convection, = ', d_h_ML_dt
!       PRINT *, '         Q_*       = ', Q_star_flk
!       PRINT *, '         \beta*Q_* = ', flk_str_1
!     ENDIF
!_dbg

      d_h_ML_dt = MAX(d_h_ML_dt, c_small_flk)    

      ! Update h_ML 
      h_ML_n_flk = h_ML_p_flk + d_h_ML_dt*del_time

      ! Security, limit h_ML
      h_ML_n_flk = MAX(h_ML_min_flk, MIN(h_ML_n_flk, depth_w))
    ELSE
      ! Mixing down to the lake bottom
      h_ML_n_flk = depth_w
    END IF

  ELSE Mixing_regime                              ! Wind mixing

    ! The surface friction velocity
    d_h_ML_dt  = MAX(u_star_w_flk, u_star_min_flk)

    ZM_h_scale = (ABS(par_Coriolis)/c_sbl_ZM_n + N_T_mean/c_sbl_ZM_i)*     &
                                                       d_h_ML_dt**2_iintegers
    ZM_h_scale = ZM_h_scale + flk_str_1/c_sbl_ZM_s
    ZM_h_scale = MAX(ZM_h_scale, c_small_flk)
    ZM_h_scale = d_h_ML_dt**3_iintegers/ZM_h_scale 

    ! The ZM96 SBL depth scale 
    ZM_h_scale = MAX(h_ML_min_flk, MIN(ZM_h_scale, h_ML_max_flk))

    ! Equilibrium mixed-layer depth 
    ZM_h_scale = MAX(ZM_h_scale, conv_equil_h_scale)

!_dm 
! In order to avoid numerical discretization problems,
! an analytical solution to the evolution equation 
! for the wind-mixed layer depth is used.
! That is, an exponential relaxation formula is applied
! over the time interval equal to the model time step.
!_dm 

    d_h_ML_dt = c_relax_h*d_h_ML_dt/ZM_h_scale*del_time
    ! Update h_ML 
    h_ML_n_flk = ZM_h_scale - (ZM_h_scale-h_ML_p_flk)*EXP(-d_h_ML_dt)

    ! Limit h_ML 
    h_ML_n_flk = MAX(h_ML_min_flk, MIN(h_ML_n_flk, depth_w))

    ! Re-compute dh_ML/dt
    d_h_ML_dt = (h_ML_n_flk-h_ML_p_flk)/del_time

    IF(h_ML_n_flk.LE.h_ML_p_flk)           &
      ! Mixed-layer retreat or stationary state, dC_T/dt<0
      d_C_T_dt = -d_C_T_dt

    C_T_n_flk = C_T_p_flk + d_C_T_dt*del_time           ! Update C_T
    C_T_n_flk = MIN(C_T_max, MAX(C_T_n_flk, C_T_min))   ! Limit C_T 
    d_C_T_dt = (C_T_n_flk-C_T_p_flk)/del_time           ! Re-compute dC_T/dt


!_dbg
!   IF (idbg > 10) THEN
!     PRINT *, 'FLake: wind mixing: d_h_ML_dt*del_time = ', d_h_ML_dt*del_time
!     PRINT *, '              h_CBL_eqi = ', conv_equil_h_scale
!     PRINT *, '              ZM scale  = ', ZM_h_scale
!     PRINT *, '              w_*       = ', w_star_sfc_flk
!     PRINT *, '              u_*       = ', u_star_w_flk
!     PRINT *, '             h_ML_p_flk = ', h_ML_p_flk
!   ENDIF
!_dbg

  END IF Mixing_regime

! Compute the time-rate-of-change of the the bottom temperature, 
! depending on the sign of dh_ML/dt 
! Update the bottom temperature and the mixed-layer temperature

  IF(h_ML_n_flk <= depth_w-h_ML_min_flk) THEN       
    ! Mixing did not reach the bottom 

    IF(h_ML_n_flk > h_ML_p_flk) THEN   
      ! Mixed-layer deepening 
      R_H_icesnow     = h_ML_p_flk/depth_w
      R_rho_c_icesnow = 1._wp-R_H_icesnow 
      R_TI_icesnow    = 0.5_wp*C_T_p_flk*R_rho_c_icesnow+C_TT_flk*     &
                                           (2._wp*R_H_icesnow-1._wp)
      R_Tstar_icesnow = (0.5_wp+C_TT_flk-C_Q_flk)/R_TI_icesnow
      R_TI_icesnow    = (1._wp-C_T_p_flk*R_rho_c_icesnow)/R_TI_icesnow
     
      d_T_bot_dt = (Q_w_flk-Q_bot_flk+I_w_flk-I_bot_flk)/tpl_rho_w_r/tpl_c_w
      d_T_bot_dt = d_T_bot_dt - C_T_p_flk*(T_wML_p_flk-T_bot_p_flk)*d_h_ML_dt

      ! Q+I fluxes and dh_ML/dt term
      d_T_bot_dt = d_T_bot_dt*R_Tstar_icesnow/depth_w

      flk_str_2 = I_intm_h_D_flk - (1._wp-C_Q_flk)*I_h_flk - C_Q_flk * &
                                                                 I_bot_flk
      flk_str_2 = flk_str_2*R_TI_icesnow/(depth_w-h_ML_p_flk)/tpl_rho_w_r/ &
                                                                 tpl_c_w
      ! Add radiation-flux term
      d_T_bot_dt = d_T_bot_dt + flk_str_2

      flk_str_2 = (1._wp-C_TT_2*R_TI_icesnow)/C_T_p_flk
      flk_str_2 = flk_str_2*(T_wML_p_flk-T_bot_p_flk)*d_C_T_dt

      ! Add dC_T/dt term
      d_T_bot_dt = d_T_bot_dt + flk_str_2
      
    ELSE
      ! Mixed-layer retreat or stationary state
      ! dT_bot/dt=0
      d_T_bot_dt = 0._wp                                            
    END IF

    ! Update T_bot  
    T_bot_n_flk = T_bot_p_flk + d_T_bot_dt*del_time   

    ! Security, limit T_bot by the freezing point
    T_bot_n_flk = MAX(T_bot_n_flk, tpl_T_f)

    flk_str_2 = (T_bot_n_flk-tpl_T_r)*flake_buoypar(T_mnw_n_flk)

    ! Security, avoid T_r crossover 
    IF(flk_str_2.LT.0._wp) T_bot_n_flk = tpl_T_r  

    T_wML_n_flk = C_T_n_flk*(1._wp-h_ML_n_flk/depth_w)
    T_wML_n_flk = (T_mnw_n_flk-T_bot_n_flk*T_wML_n_flk)/(1._wp-T_wML_n_flk)

    ! Security, limit T_wML by the freezing point
    T_wML_n_flk = MAX(T_wML_n_flk, tpl_T_f)

  ELSE
    ! Mixing down to the lake bottom 

    h_ML_n_flk = depth_w
    T_wML_n_flk = T_mnw_n_flk
    T_bot_n_flk = T_mnw_n_flk
    C_T_n_flk = C_T_min

  END IF

END IF HTC_Water

!------------------------------------------------------------------------------
!  Compute the depth of the upper layer of bottom sediments
!  and the temperature at that depth.
!------------------------------------------------------------------------------

Use_sediment: IF(lflk_botsed_use) THEN   ! The bottom-sediment scheme is used
  
  IF (H_B1_p_flk >= depth_bs-H_B1_min_flk) THEN  
    ! No T(z) maximum (no thermal wave) 
    H_B1_p_flk = 0._wp               ! Set H_B1_p to zero
    T_B1_p_flk = T_bot_p_flk             ! Set T_B1_p to the bottom temperature
  END IF 

  flk_str_1 = 2._wp*Phi_B1_pr0/(1._wp-C_B1)*tpl_kappa_w/tpl_rho_w_r/ &
                                                          tpl_c_w*del_time
  ! Threshold value of H_B1
  h_ice_threshold = SQRT(flk_str_1)

  ! Limit H_B1
  h_ice_threshold = MIN(0.9_wp*depth_bs, h_ice_threshold)    

  flk_str_2 = C_B2/(1._wp-C_B2)*(T_bs-T_B1_p_flk)/(depth_bs-H_B1_p_flk)

  IF (H_B1_p_flk < h_ice_threshold) THEN
    ! Use a truncated equation for H_B1(t)
    H_B1_n_flk = SQRT(H_B1_p_flk**2_iintegers+flk_str_1)  ! Advance H_B1
    d_H_B1_dt = (H_B1_n_flk-H_B1_p_flk)/del_time          ! Re-compute dH_B1/dt 
  ELSE
    ! Use a full equation for H_B1(t)
    flk_str_1 = (Q_bot_flk+I_bot_flk)/H_B1_p_flk/tpl_rho_w_r/tpl_c_w
    flk_str_1 = flk_str_1 - (1._wp-C_B1)*(T_bot_n_flk-T_bot_p_flk)/del_time
    d_H_B1_dt = (1._wp-C_B1)*(T_bot_p_flk-T_B1_p_flk)/H_B1_p_flk +    &
                                                              C_B1*flk_str_2
    d_H_B1_dt = flk_str_1/d_H_B1_dt
    H_B1_n_flk = H_B1_p_flk + d_H_B1_dt*del_time          ! Advance H_B1
  END IF 
  d_T_B1_dt = flk_str_2*d_H_B1_dt
  T_B1_n_flk = T_B1_p_flk + d_T_B1_dt*del_time            ! Advance T_B1


!_dbg
! IF (idbg > 10) THEN
!   PRINT *, 'BS module: '
!   PRINT *, '  Q_bot   = ', Q_bot_flk
!   PRINT *, '  d_H_B1_dt = ', d_H_B1_dt
!   PRINT *, '  d_T_B1_dt = ', d_T_B1_dt
!   PRINT *, '  H_B1    = ', H_B1_n_flk
!   PRINT *, '    T_bot = ', T_bot_n_flk
!   PRINT *, '  T_B1    = ', T_B1_n_flk
!   PRINT *, '    T_bs  = ',  T_bs
! ENDIF
!_dbg

!_nu  
! Use a very simplistic procedure, where only the upper layer profile is used, 
! H_B1 is always set to depth_bs, and T_B1 is always set to T_bs. Then, the 
! time derivatives are zero, and the sign of the bottom heat flux depends on 
! whether T_bot is smaller or greater than T_bs.
! This is, of course, an oversimplified scheme.
!_nu  d_H_B1_dt = 0._wp
!_nu  d_T_B1_dt = 0._wp
!_nu  H_B1_n_flk = H_B1_p_flk + d_H_B1_dt*del_time   ! Advance H_B1
!_nu  T_B1_n_flk = T_B1_p_flk + d_T_B1_dt*del_time   ! Advance T_B1
!_nu  

  l_snow_exists = (H_B1_n_flk >= depth_bs-H_B1_min_flk)       &
                                             ! H_B1 reached depth_bs, or
             .OR. (H_B1_n_flk <  H_B1_min_flk)                &
                                             ! H_B1 decreased to zero, or
             .OR. ((T_bot_n_flk-T_B1_n_flk)*(T_bs-T_B1_n_flk) <= 0._wp)
                                             ! there is no T(z) maximum
  IF(l_snow_exists) THEN      
    H_B1_n_flk = depth_bs ! Set H_B1 to the depth of the thermally active layer
    T_B1_n_flk = T_bs     ! Set T_B1 to the climatological temperature 
  END IF

ELSE Use_sediment
  ! The bottom-sediment scheme is not used

  ! H_B1 is set to a reference value 
  H_B1_n_flk = rflk_depth_bs_ref   

  ! T_B1 is set to the temperature of maximum density
  T_B1_n_flk = tpl_T_r

END IF Use_sediment

!------------------------------------------------------------------------------
!  Impose additional constraints.
!------------------------------------------------------------------------------

! In case of unstable stratification, force mixing down to the bottom
flk_str_2 = (T_wML_n_flk-T_bot_n_flk)*flake_buoypar(T_mnw_n_flk)
IF(flk_str_2.LT.0._wp) THEN 

!_dbg
!IF (idbg > 10) THEN
!  PRINT *, 'FLake: inverse (unstable) stratification !!! '
!  PRINT *, '       Mixing down to the bottom is forced.'
!  PRINT *, '  T_wML_p, T_wML_n ', T_wML_p_flk-tpl_T_f, T_wML_n_flk-tpl_T_f
!  PRINT *, '  T_mnw_p, T_mnw_n ', T_mnw_p_flk-tpl_T_f, T_mnw_n_flk-tpl_T_f
!  PRINT *, '  T_bot_p, T_bot_n ', T_bot_p_flk-tpl_T_f, T_bot_n_flk-tpl_T_f
!  PRINT *, '  h_ML_p,  h_ML_n  ', h_ML_p_flk,          h_ML_n_flk
!  PRINT *, '  C_T_p,   C_T_n   ', C_T_p_flk,           C_T_n_flk
!ENDIF
!_dbg

  h_ML_n_flk = depth_w
  T_wML_n_flk = T_mnw_n_flk
  T_bot_n_flk = T_mnw_n_flk
  C_T_n_flk = C_T_min

END IF


!------------------------------------------------------------------------------
!  Update the surface temperature.
!------------------------------------------------------------------------------

IF     (h_snow_n_flk >= h_Snow_min_flk) THEN   
  ! Snow exists, use the snow temperature
  T_sfc_n = T_snow_n_flk
ELSEIF (h_ice_n_flk >= h_Ice_min_flk) THEN
  ! Ice exists but there is no snow, use the ice temperature
  T_sfc_n = T_ice_n_flk
ELSE 
  ! No ice-snow cover, use the mixed-layer temperature
  T_sfc_n = T_wML_n_flk
END IF

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END SUBROUTINE flake_driver

!==============================================================================
!==============================================================================
!------------------------------------------------------------------------------

REAL (KIND = wp)     FUNCTION flake_buoypar (T_water)

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the buoyancy parameter,
!  using a quadratic equation of state for the fresh-water.
!  
!------------------------------------------------------------------------------

!  Input (function argument) 
REAL (KIND = wp),     INTENT(IN) :: &
  T_water                             ! Water temperature [K]

!------------------------------------------------------------------------------
!  Start calculations
!------------------------------------------------------------------------------

! Buoyancy parameter [m s^{-2} K^{-1}]

  flake_buoypar = tpl_grav*tpl_a_T*(T_water-tpl_T_r)

!------------------------------------------------------------------------------
!  End calculations
!------------------------------------------------------------------------------

END FUNCTION flake_buoypar

!==============================================================================
!==============================================================================
!------------------------------------------------------------------------------

REAL (KIND = wp)     FUNCTION flake_snowdensity (hz_snow)

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the snow density,
!  using an empirical approximation from Heise et al. (2003).
!  
!------------------------------------------------------------------------------

!  Input (function argument) 
REAL (KIND = wp),     INTENT(IN) :: &
  hz_snow                              ! Snow thickness [m]

!------------------------------------------------------------------------------
!  Start calculations
!------------------------------------------------------------------------------

! Snow density [kg m^{-3}]

! Security. Ensure that the expression in () does not become negative at a 
! very large hz_snow.
  flake_snowdensity = MAX( c_small_flk,                                 &
                          (1._wp - hz_snow*tpl_Gamma_rho_S/tpl_rho_w_r) )
  flake_snowdensity = MIN( tpl_rho_S_max, tpl_rho_S_min/flake_snowdensity )

!------------------------------------------------------------------------------
!  End calculations
!------------------------------------------------------------------------------

END FUNCTION flake_snowdensity 

!==============================================================================
!==============================================================================
!------------------------------------------------------------------------------

REAL (KIND = wp)     FUNCTION flake_snowheatconduct (hz_snow)

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the snow heat conductivity,
!  using an empirical approximation from Heise et al. (2003).
!  
!------------------------------------------------------------------------------
 
!  Input (function argument) 
REAL (KIND = wp),     INTENT(IN) :: &
  hz_snow                              ! Snow thickness [m]

!------------------------------------------------------------------------------
!  Start calculations
!------------------------------------------------------------------------------

! Snow heat conductivity [J m^{-1} s^{-1} K^{-1} = kg m s^{-3} K^{-1}]

  ! Compute snow density
  flake_snowheatconduct = flake_snowdensity( hz_snow )   

  flake_snowheatconduct = MIN( tpl_kappa_S_max, tpl_kappa_S_min            &
              + hz_snow*tpl_Gamma_kappa_S*flake_snowheatconduct/tpl_rho_w_r )

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END FUNCTION flake_snowheatconduct

!==============================================================================

END MODULE src_flake

