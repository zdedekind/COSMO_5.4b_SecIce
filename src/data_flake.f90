!+ Data module for variables concerned with the FLake Model
!------------------------------------------------------------------------------

MODULE data_flake

!------------------------------------------------------------------------------
!
! Description:
!  This module contains common definitions, variables and constants for the
!  FLake Model. These are:
!
!  - flake_albedo_ref:
!   "reference" values of albedo for the lake water, lake ice and snow. 
!    As in "flake_paramoptic_ref", two ice categories, viz. white ice and blue
!    ice, and two snow categories, viz. dry snow and melting snow, are used.  
!
!  - flake_configure:
!    Switches and reference values of parameters 
!    that configure the lake model FLake are set.
!
!  - flake_derivedtypes:
!    Derived type(s) is(are) defined.
!
!  - flake_parameters:
!    Values of empirical constants of the lake model FLake 
!    and of several thermodynamic parameters are set.
!
!  - flake_paramoptic_ref:
!    This module contains "reference" values of the optical characteristics
!    of the lake water, lake ice and snow.
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
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_24        2012/06/22 Dmitrii Mironov
!  Introduced a maximum value for EXP arguments
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Modules used:

USE data_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations
 
!  flake_albedo_ref:
!  -----------------

!  Albedo for water, ice and snow.
REAL (KIND = wp),     PARAMETER ::    &
  albedo_water_ref       = 0.07_wp  , & ! Water
  albedo_whiteice_ref    = 0.60_wp  , & ! White ice
  albedo_blueice_ref     = 0.10_wp  , & ! Blue ice
  albedo_drysnow_ref     = 0.60_wp  , & ! Dry snow 
  albedo_meltingsnow_ref = 0.10_wp      ! Melting snow 

!  Empirical parameters.
REAL (KIND = wp),     PARAMETER ::    &
  c_albice_MR = 95.6_wp                 ! Constant in the interpolation formula for
                                        ! the ice albedo (Mironov and Ritter 2004)

!------------------------------------------------------------------------------

!  flake_configure:
!  ----------------

! Switches and reference values of parameters
LOGICAL, PARAMETER :: &
  lflk_botsed_use   = .FALSE.      

  ! .TRUE. indicates that the bottom-sediment scheme is used to compute the
  ! depth penetrated by the thermal wave, the temperature at this depth and
  ! the bottom heat flux.  Otherwise, the heat flux at the water-bottom
  ! sediment interface is set to zero, the depth penetrated by the thermal
  ! wave is set to a reference value defined below, and the temperature at
  ! this depth is set to the temperature of maximum density of the fresh water.

REAL (KIND = wp),     PARAMETER :: &
  rflk_depth_bs_ref = 10.0_wp  

  ! Reference value of the depth of the thermally active layer of bottom
  ! sediments [m].  This value is used to (formally) define the depth
  ! penetrated by the thermal wave in case the bottom-sediment scheme
  ! is not used.

!------------------------------------------------------------------------------

!  flake_derived_types:
!  --------------------

!  Derived type(s) is(are) defined.

!  Maximum value of the wave-length bands 
!  in the exponential decay law for the radiation flux.
!  A storage for a ten-band approximation is allocated,
!  although a smaller number of bands is actually used.

INTEGER (KIND = iintegers), PARAMETER :: & 
  nband_optic_max = 10_iintegers

!  Define TYPE "opticpar_medium"
TYPE opticpar_medium

  INTEGER (KIND = iintegers)          :: & 
    nband_optic                            ! Number of wave-length bands

  REAL (KIND = wp)                    :: & 
    frac_optic       (nband_optic_max),  & ! Fractions of total radiation flux 
    extincoef_optic  (nband_optic_max)     ! Extinction coefficients 

END TYPE opticpar_medium

!------------------------------------------------------------------------------

! flake_parameters:
! -----------------

!  Values of empirical constants of the lake model FLake 
!  and of several thermodynamic parameters are set.

!  Dimensionless constants in the equations for the mixed-layer depth and for
!  the shape factor with respect to the temperature profile in the thermocline

REAL (KIND = wp),     PARAMETER ::      &
  c_cbl_1    = 0.17_wp , & ! Constant in the CBL entrainment equation
  c_cbl_2    = 1.0_wp  , & ! Constant in the CBL entrainment equation
  c_sbl_ZM_n = 0.5_wp  , & ! Constant in the ZM1996 equation for the 
                           !       equilibrium SBL depth
  c_sbl_ZM_s = 10.0_wp , & ! Constant in the ZM1996 equation for the 
                           !       equilibrium SBL depth
  c_sbl_ZM_i = 20.0_wp , & ! Constant in the ZM1996 equation for the 
                           !       equilibrium SBL depth
  c_relax_h  = 0.010_wp, & ! Constant in the relaxation equation for the 
                           !       SBL depth
  c_relax_C  = 0.0030_wp   ! Constant in the relaxation equation for the 
                           !       shape factor with respect to the 
                           !       temperature profile in the thermocline

!  Parameters of the shape functions 
!  Indices refer to T - thermocline, S - snow, I - ice,
!  B1 - upper layer of the bottom sediments, 
!  B2 - lower layer of the bottom sediments.
!  "pr0" and "pr1" denote zeta derivatives of the corresponding shape function 
!  at "zeta=0" ad "zeta=1", respectively.

REAL (KIND = wp), PARAMETER ::   &
  C_T_min       = 0.5_wp       , & ! Minimum value of the shape factor C_T (thermocline)
  C_T_max       = 0.8_wp       , & ! Maximum value of the shape factor C_T (thermocline)
  Phi_T_pr0_1   = 40._wp/3._wp , & ! Constant in the expression for the T shape-function derivative
  Phi_T_pr0_2   = 20._wp/3._wp , & ! Constant in the expression for the T shape-function derivative
  C_TT_1        = 11._wp/18._wp, & ! Constant in the expression for C_TT (thermocline)
  C_TT_2        = 7._wp/45._wp , & ! Constant in the expression for C_TT (thermocline)
  C_B1          = 2._wp/3._wp  , & ! Shape factor (upper layer of bottom sediments)
  C_B2          = 3._wp/5._wp  , & ! Shape factor (lower layer of bottom sediments)
  Phi_B1_pr0    = 2._wp        , & ! B1 shape-function derivative
  C_S_lin       = 0.5_wp       , & ! Shape factor (linear temperature profile in the snow layer)
  Phi_S_pr0_lin = 1._wp        , & ! S shape-function derivative (linear profile)
  C_I_lin       = 0.5_wp       , & ! Shape factor (linear temperature profile in the ice layer)
  Phi_I_pr0_lin = 1._wp        , & ! I shape-function derivative (linear profile)
  Phi_I_pr1_lin = 1._wp        , & ! I shape-function derivative (linear profile)
  Phi_I_ast_MR  = 2._wp        , & ! Constant in the MR2004 expression for I shape factor
  C_I_MR        = 1._wp/12._wp , & ! Constant in the MR2004 expression for I shape factor
  H_Ice_max     = 3._wp            ! Maximum ice tickness in the Mironov and Ritter ice model [m]

!  Security constants
REAL (KIND = wp), PARAMETER ::   &
  h_Snow_min_flk = 1.0E-5_wp   , & ! Minimum snow thickness [m]
  h_Ice_min_flk  = 1.0E-9_wp   , & ! Minimum ice thickness [m]
  h_ML_min_flk   = 1.0E-2_wp   , & ! Minimum mixed-layer depth [m]
  h_ML_max_flk   = 1.0E+3_wp   , & ! Maximum mixed-layer depth [m]
  H_B1_min_flk   = 1.0E-3_wp   , & ! Minimum thickness of the upper layer of bottom sediments [m]
  u_star_min_flk = 1.0E-6_wp       ! Minimum value of the surface friction velocity [m s^{-1}]

!  Security constant(s)
REAL (KIND = wp), PARAMETER ::   &
  c_small_flk    = 1.0E-10_wp  , & ! A small number
  c_maxearg_flk  = 1.0E+02_wp      ! Maximum value of the EXP function argument [-]

!  Thermodynamic parameters
REAL (KIND = wp), PARAMETER ::       &
  tpl_grav          = 9.81_wp      , & ! Acceleration due to gravity [m s^{-2}]
  tpl_T_r           = 277.13_wp    , & ! Temperature of maximum density of fresh water [K]
  tpl_T_f           = 273.15_wp    , & ! Fresh water freezing point [K]
  tpl_a_T           = 1.6509E-05_wp, & ! Constant in the fresh-water equation of state [K^{-2}]
  tpl_rho_w_r       = 1.0E+03_wp   , & ! Maximum density of fresh water [kg m^{-3}]
  tpl_rho_I         = 9.1E+02_wp   , & ! Density of ice [kg m^{-3}]
  tpl_rho_S_min     = 1.0E+02_wp   , & ! Minimum snow density [kg m^{-3}]
  tpl_rho_S_max     = 4.0E+02_wp   , & ! Maximum snow density [kg m^{-3}]
  tpl_Gamma_rho_S   = 2.0E+02_wp   , & ! Empirical parameter [kg m^{-4}] in
                                       !   the expression for the snow density 
  tpl_L_f           = 3.3E+05_wp   , & ! Latent heat of fusion      [J kg^{-1}]
  tpl_c_w           = 4.2E+03_wp   , & ! Specific heat of water     [J kg^{-1} K^{-1}]
  tpl_c_I           = 2.1E+03_wp   , & ! Specific heat of ice       [J kg^{-1} K^{-1}]
  tpl_c_S           = 2.1E+03_wp   , & ! Specific heat of snow      [J kg^{-1} K^{-1}]
  tpl_kappa_w       = 5.46E-01_wp  , & ! Molecular heat conductivity of water 
                                       !   [J m^{-1} s^{-1} K^{-1}]
  tpl_kappa_I       = 2.29_wp      , & ! Molecular heat conductivity of ice   
                                       !   [J m^{-1} s^{-1} K^{-1}]
  tpl_kappa_S_min   = 0.2_wp       , & ! Minimum molecular heat conductivity  
                                       !   of snow [J m^{-1} s^{-1} K^{-1}]
  tpl_kappa_S_max   = 1.5_wp       , & ! Maximum molecular heat conductivity  
                                       !   of snow [J m^{-1} s^{-1} K^{-1}]
  tpl_Gamma_kappa_S = 1.3_wp           ! Empirical parameter in expression for the 
                                       !   snow heat conductivity [J m^{-2} s^{-1} K^{-1}] 

!------------------------------------------------------------------------------

! flake_paramoptic_ref:
! ---------------------

!  This module contains "reference" values of the optical characteristics
!  of the lake water, lake ice and snow. These reference values may be used 
!  if no information about the optical characteristics of the lake in question 
!  is available. An exponential decay law for the solar radiation flux is assumed.
!  In the simplest one-band approximation,
!  the extinction coefficient for water is set to a large value,
!  leading to the absorption of 95% of the incoming radiation 
!  within the uppermost 1 m of the lake water. 
!  The extinction coefficients for ice and snow are taken from 
!  Launiainen and Cheng (1998). The estimates for the ice correspond 
!  to the uppermost 0.1 m of the ice layer and to the clear sky conditions 
!  (see Table 2 in op. cit.).
!  Very large values of the extinction coefficients for ice and snow ("opaque")
!  can be used to prevent penetration of the solar radiation 
!  through the snow-ice cover.
 
INTEGER (KIND = iintegers), PRIVATE :: & ! Help variable(s)
  i                                      ! DO loop index

!  Optical characteristics for water, ice and snow.
!  The simplest one-band approximation is used as a reference.
TYPE (opticpar_medium), PARAMETER ::                           & 
  ! Water (reference)
  opticpar_water_ref = opticpar_medium(1,                      &
    (/1._wp, (0._wp,i=2,nband_optic_max)/),                    &
    (/3._wp, (1.E+10_wp,i=2,nband_optic_max)/))              , &

  ! Transparent Water (two-band)
  opticpar_water_trans = opticpar_medium(2,                    & 
    (/0.10_wp, 0.90_wp, (0._wp,i=3,nband_optic_max)/),         &
    (/2.0_wp, 0.20_wp, (1.E+10_wp,i=3,nband_optic_max)/))    , &
  ! Transparent Water (one-band)
!_nu  opticpar_water_trans = opticpar_medium(1,                   & 
!_nu    (/1._wp, (0._wp,i=2,nband_optic_max)/),                   &
!_nu    (/0.30_wp, (1.E+10_wp,i=2,nband_optic_max)/))           , &
  ! White ice
  opticpar_whiteice_ref = opticpar_medium(1,                   & 
    (/1._wp, (0._wp,i=2,nband_optic_max)/),                    &   
    (/17.1_wp, (1.E+10_wp,i=2,nband_optic_max)/))            , &
  ! Blue ice
  opticpar_blueice_ref = opticpar_medium(1,                    & 
    (/1._wp, (0._wp,i=2,nband_optic_max)/),                    &
    (/8.4_wp, (1.E+10_wp,i=2,nband_optic_max)/))             , &
  ! Dry snow 
  opticpar_drysnow_ref = opticpar_medium(1,                    & 
    (/1._wp, (0._wp,i=2,nband_optic_max)/),                    &
    (/25.0_wp, (1.E+10_wp,i=2,nband_optic_max)/))            , &
  ! Melting snow 
  opticpar_meltingsnow_ref = opticpar_medium(1,                & 
    (/1._wp, (0._wp,i=2,nband_optic_max)/),                    &
    (/15.0_wp, (1.E+10_wp,i=2,nband_optic_max)/))            , &
  ! Opaque ice
  opticpar_ice_opaque = opticpar_medium(1,                     & 
    (/1._wp, (0._wp,i=2,nband_optic_max)/),                    &
    (/1.0E+07_wp, (1.E+10_wp,i=2,nband_optic_max)/))         , &
  ! Opaque snow
  opticpar_snow_opaque = opticpar_medium(1,                    & 
    (/1._wp, (0._wp,i=2,nband_optic_max)/),                    &
    (/1.0E+07_wp, (1.E+10_wp,i=2,nband_optic_max)/)) 

!==============================================================================

END MODULE data_flake

