!+ Source module  "gscp"
!------------------------------------------------------------------------------

MODULE src_gscp    

!------------------------------------------------------------------------------
!
! Description:
!
!   The module "gscp" calculates the rates of change of temperature,
!   cloud condensate and water vapor due to cloud microphysical processes
!   related to the formation of grid scale clouds and precipitation.
!   The microphysical subroutines are either called from "organize_gscp"
!   (organizational subroutine in this module, called from organize_physics)
!   or from "organize_physics" itself.
!
!   The following subroutines may be called:
!
!   With prognostic precipitation turned on (from Version 4.23 only this is possible)
!
!      called form organize_physics at the end of a time step:
!        (1) a warm rain Kessler parameterization (itype_gscp=1) named "kessler_pp"
!        (2) a scheme including snow              (itype_gscp=2) named "hydor_pp"
!        (3) a scheme with prognostic qr,qs,qi    (itype_gscp=3) named "hydci_pp"
!        (4) a scheme with prognostic qr,qs,qi,qg (itype_gscp=4) named "hydci_pp_gr"
!
!      For diagnostic initialization of qr and qs (ldiniprec=.TRUE)
!       called from organize_gscp (at the beginning of a time step):
!        (3) a cloud ice scheme (itype_gscp=3) named "hydci" (with linionly)
!
! Current Code Owner: DWD, Axel Seifert
!  phone:  +49  69  8062 2729
!  fax:    +49  69  8062 3721
!  email:  Axel.Seifert@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Guenther Doms     
!  Initial release
! 1.3        1998/04/15 Guenther Doms
!  Some minor modifications in Hydor and Kessler for optimization
! 1.4        1998/05/22 Guenther Doms
!  Inclusion of control parameter l2tls to select time levels
!  according to the time integration scheme used.
! 1.20       1999/01/07 Guenther Doms
!  Renaming of some global variables
! 1.24       1999/03/01 Guenther Doms
!  Inclusion of the new optional cloud ice scheme ('hydci').
! 1.29       1999/05/11 Ulrich Schaettler
!  Corrections to intrinsic calls for use of non-Cray machines
! 2.8        2001/07/06 Ulrich Schaettler
!  Introduced optimization for vectorization
! 2.14       2002/02/15 Guenther Doms
!  Reduction of the autoconversion coefficient of cloud water (zccau)
!  from 7.E-4 to 4.E-4 in the cloud ice scheme, to be conform with GME.
! 2.17       2002/05/08 Ulrich Schaettler
!  Corrections to get reproducible results in the cloud-ice scheme
! 2.18       2002/07/16 Guenther Doms
!  Evaporation of melting snow is also considered in the cloud ice scheme.
!  Bug correction in hydor and changed some coefficients in hydor and hydci
!  (by Almut Gassmann)
!  New Routines (hydorprog, hydciprog) for prognostic precipitation
!  (by Almut Gassmann)
! 2.19       2002/10/24 Ulrich Schaettler
!  Eliminated evaporation of melting snow again because of problems in
!  the tropical regions.
!  Changes to the prognostic precipitation routines by Almut Gassmann.
! 3.2        2003/02/07 Ulrich Schaettler
!  Set cloud ice threshold (qi0) for autoconversion to 0 (in hydci).
! 3.6        2003/12/11 Ulrich Schaettler
!  Replaced most power-functions by EXP (... * LOG(...)) for optimization
! 3.7        2004/02/18 Guenther Doms
!  New Subroutine hydci_pplf for prognostic treatment of precipitation
!  (only possible, if running with cloud ice)
!  Eliminated a CYCLE-Statement in Section 2 of hydorprog to get reproducible
!  results  (by Jochen Foerstner)
! 3.8        2004/03/23 Jochen Foerstner / Thorsten Reinhardt
!  Technical changes for using hydci_pp also for the Runge-Kutta scheme
!  Bug-fix and additional optimizations in hydci_pp
! 3.12       2004/09/15 Guenther Doms / Christoph Schraff
!  Diagnostic initialisation of prognostic rain and snow.
! 3.13       2004/12/03 Ulrich Schaettler
!  Include some variables and routines for the Latent Heat Nudging,
!  (which calls e.g. hydci also in case of lprogprec=.TRUE.)
!  Get diagnostic precipitation for LHN in case of lprogprec == .true.
!                                              (Klaus Stephan)
!  Introduction of graupel scheme; New routines hydor_pp and kessler_pp
!  for prognostic treatment of precipitation als in these cases
!                                              (Thorsten Reinhardt)
! 3.14       2005/01/25 Thorsten Reinhardt
!  Adaptations for latent heat nudging also in hydor_pp, kessler_pp, hydci_pp_gr
! 3.16       2005/07/22 Thorsten Reinhardt
!  Some editorial changes
! 3.18       2006/03/03 Klaus Stephan / Mauro Ballabio / Thorsten Reinhardt
!  LHN namelist parameter moved to data_lheat_nudge to avoid to many dependencies
!  Include some new variables and routines for the Latent Heat Nudging
!  Initialization of some variables in the Graupel scheme for vectorization
!  Introduction of temperature increment due to latent heat (in RK scheme)
! 3.21       2006/12/04 Ulrich Schaettler
!  Use new NL variables qc0, qi0 from data_constants
!  Changes in hydci_pp_gr and hydci_gr_lhn (by T. Reinhardt, K. Stephan)
! 3.22       2007/01/24 Axel Seifert
!  New version of hydci_pp with major changes in the parameterization of
!  snow: Variable intercept parameter based on Field et al. (2005),
!  temperature-dependent sticking efficiency, changes in geometry and
!  fall speed of snow.
!  Use of Seifert and Beheng (2001) warm rain scheme with a constant
!  cloud droplet number concentration.
! V3_23        2007/03/30 Axel Seifert, Jochen Foerstner, Lucio Torrisi
!  Added initialization of several scalar variables
!  Added correction of tinc_lh for ldiabf_lh in hydci_pp
!  Some corrections for Latent Heat Nudging (Klaus Stephan)
! V3_24        2007/04/26 Axel Seifert
!  Adapted hydci_pp_gr to the new version hydci_pp
!  Eliminated subroutines hydci_lhn and hydci_gr_lhn, which are not used 
!  any more
! V3_27        2007/06/29 Ulrich Schaettler
!  Editorial changes
! V4_5         2008/09/10 Ulrich Schaettler, Jens-Olaf Beismann
!  Vector optimizations in SR hydci_pp and hydci_pp_gr
!  Moved declaration of mu_rain and cloud_num (before: mur, zcnum) to data_gscp
!  Some variables are determined depending on value of mu_rain (Christoph Gebhardt)
! V4_9         2009/07/16 Ulrich Schaettler, Christian Bollmann
!  Added VOB (vector overtake barrier) in NEC Compiler Directives
!  Additional vector optimizations
! V4_10        2009/09/11 Axel Seifert
!  Bug fix in SR hydci_pp_gr: Variable scac was set twice
!  Added compiler directive to use option _on_adb for NEC
! V4_12        2010/05/11 Ulrich Schaettler
!  Renamed t0 to t0_melt because of conflicting names
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_14        2010/06/14 Axel Seifert
!  Density correction of terminal velocity for snow flakes
! V4_18        2011/05/26 Ulrich Schaettler
!  Introduced conditional compilation for NUDGING
! V4_20        2011/08/31 Axel Seifert
!  Correction of size distribution for rain droplets;
!  generalization of implementation for arbitrary mu_rain;
!  bugfix in density correction of fall speeds
! V4_21        2011/12/06 Axel Seifert, Uli Blahak
!  Introduction of new factor rain_n0_factor
!  Bug fix in hydci_pp, hydci_pp_gr: a factor of rho has been forgotten in a 
!    computation
! V4_23        2012/05/10 Ulrich Schaettler
!  Removed src_2timelevel and related stuff (hydorprog, hydciprog)
!  Removed switch lprogprec and corresponding subroutines
!    that were called only if lprogprec=.FALSE. (kessler, hydor)
!    Remaining subroutines: kessler_pp, hydor_pp, hydci_pp, hydci_pp_gr
!                and        hydci (for linionly=.TRUE.)
!    Also removed organize_gscp, because there is nothing left in that SR
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer
!  Replaced qx-variables by using them from the tracer module
! V4_27        2013/03/19 Felix Rieper
!  Modified several physical coefficients which were wrongly initialized:
!   cloud ice scheme hydci_pp:
!     ccsdxp: now use variable zv1s instead of zbms for initialization
!     zbsmel, ccsdep, zbev: factor sqrt(0.5) must be eliminated
!   graupel scheme hydci_pp_gr:
!     ccsdxp: now use variable zv1s instead of zbms for initialization
!  Adapted values of variables zcsdep, zsamel to lately modified diffusion 
!  coefficient zdv in the cloud ice scheme hydci_pp
!  MESSy interface introduced
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
! V5_1         2014-11-28 Oliver Fuhrer, Lucio Torrisi
!  Replaced ireals by wp (working precision) (OF)
!  Call to 'apply_tqx_tend_adj' added for SPPT (stochastic perturbation of
!    physics tendencies) in routines hydci_pp and hydci_pp_gr
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
    wp        ,& ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &


! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------
    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ieje,         & ! ie * je
    ke,           & ! number of grid points in vertical direction
    ke1,          & ! ke+1

! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the 
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from 
!    the other ones because of the use of the staggered Arakawa-C-grid.

    istartpar,    & ! start index for computations in the parallel program
    iendpar,      & ! end index for computations in the parallel program
    jstartpar,    & ! start index for computations in the parallel program
    jendpar,      & ! end index for computations in the parallel program
    jstart, jend, istart, iend, &

! 4. variables for the time discretization and related variables
! --------------------------------------------------------------
    dt,           & ! timestep
    dt2,          & ! 2 * dt

! 8. Organizational variables to handle the COSMO humidity tracers
! ----------------------------------------------------------------
    idt_qv, idt_qc, idt_qs, idt_qr, idt_qi, idt_qg

! end of data_modelconfig

!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &
    pi,           & ! 

! 2. physical constants and related variables
! -------------------------------------------
    t0_melt,      & ! melting temperature of ice    
    r_d,          & ! gas constant for dry air
    r_v,          & ! gas constant for water vapour
    rdv,          & ! r_d / r_v
    o_m_rdv,      & ! 1 - r_d/r_v
    rvd_m_o,      & ! r_v/r_d - 1
    cp_d,         & ! specific heat of dry air
    cpdr,         & ! 1 / cp_d
    lh_v,         & ! latent heat of vapourization
    lh_f,         & ! latent heat of fusion
    lh_s,         & ! latent heat of sublimation
    g,            & ! acceleration due to gravity

! 3. constants for parametrizations
! ---------------------------------
    b1,           & ! variables for computing the saturation vapour pressure
    b2w,          & ! over water (w) and ice (i)
    b2i,          & !               -- " --
    b3,           & !               -- " --
    b4w,          & !               -- " --
    b234w,        & !               -- " --
    b4i,          & !               -- " --

! 4. tuning constants for radiation, cloud physics, turbulence
! ------------------------------------------------------------

    qi0,          & ! cloud ice threshold for autoconversion
    qc0,          & ! cloud water threshold for autoconversion

! 5. Precision-dependent security parameters (epsilons)
! ------------------------------------------------------

    repsilon        ! precision of 1.0 in current floating point format
! end of data_constants

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------
    rho0       ,    & ! reference density at the full model levels    (kg/m3)
    dp0        ,    & ! pressure thickness of layers                  (Pa)
    p0         ,    & ! reference pressure at full levels             (Pa)
    hhl        ,    & ! height of model half levels                   (m )

! 3. prognostic variables                                             (unit)
! -----------------------
    t          ,    & ! temperature                                   (  k  )
    pp         ,    & ! deviation from the reference pressure         ( pa  )

! 4. tendency fields for the prognostic variables                     (unit )
! -----------------------------------------------
!    timely deviation  by diabatic and adiabatic processes 
!    without sound-wave terms
    ttens        ,  & ! t-tendency without sound-wave terms           ( K/s )


! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------
    tinc_lh    ,    & ! temperature increment due to latent heat      (  K  )
    rho        ,    & ! density of moist air
!   fields of the precipitation
    qrs        ,    & ! precipitation water (water loading)           (kg/kg)
    prr_gsp    ,    & ! precipitation rate of rain, grid-scale        (kg/m2*s)
    prs_gsp    ,    & ! precipitation rate of snow, grid-scale        (kg/m2*s)
    prg_gsp           ! precipitation rate of graupel, grid-scale     (kg/m2*s)

! end of data_fields

!------------------------------------------------------------------------------

USE data_gscp       ! all variables are used here

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    ntstep,         & ! actual time step
    nstart,         & ! first time step of the forecast
    nold,           & ! corresponds to ntstep - 1
    nnow,           & ! corresponds to ntstep 
    nnew,           & ! corresponds to ntstep + 1

! 3. controlling the physics      
! --------------------------
    itype_gscp,     & ! type of grid-scale precipitation physics   
    ldiniprec,      & ! diagnostic initialisation of prognostic precip (qr, qs)
    lsppt,          & ! switch stoch. physics tend. perturbation
    itype_qxpert_rn,& ! define which hum variables tend. are perturbed
    itype_qxlim_rn, & ! type of reduction/removal of the perturbation 
                      ! in case of negative (qv, qc, qi) or 
                      ! supersaturated (qv) values

! 3. controlling the dynamics     
! ---------------------------
    l2tls,          & ! forecast with 2-TL integration scheme

! 7. additional control variables 
! --------------------------------
    ldiabf_lh,      & ! include diabatic forcing due to latent heat in RK-scheme

! 12. controlling verbosity of debug output
! -----------------------------------------
    idbg_level,     & ! to control the verbosity of debug output
    ldebug_gsp,     & ! if .TRUE., debug output for grid scale precipitation
    lprintdeb_all     ! .TRUE.:  all tasks print debug output
                      ! .FALSE.: only task 0 prints debug output

! end of data_runcontrol 
!------------------------------------------------------------------------------

USE data_parallel,            ONLY : &
    my_cart_id      ! rank of this subdomain in the cartesian communicator

!------------------------------------------------------------------------------

USE environment,              ONLY : collapse, model_abort
USE pp_utilities,             ONLY : gamma_fct
USE meteo_utilities,          ONLY : satad        ! saturation adjustment
USE src_stoch_physics,        ONLY : pertstoph, apply_tqx_tend_adj

!------------------------------------------------------------------------------

USE src_tracer,               ONLY : trcr_get, trcr_errorstr

#ifdef MESSY
USE messy_main_data_bi,       ONLY: precmelt_ls, rainform_bave, snowform_bave &
                                  , precflxno_ls, snowflxno_ls, precflx_ls    &
                                  , snowflx_ls,   preccover_ls, sediice_ls
#endif

!------------------------------------------------------------------------------

#ifdef NUDGING
USE data_lheat_nudge,           ONLY :  &
    llhn,         & ! main switch for latent heat nudging
    llhnverif,    & ! main switch for latent heat nudging
    lhn_qrs,      & ! use integrated precipitaion flux as reference
    tt_lheat,     & ! latent heat release of the model
    qrsflux         ! total precipitation flux

!------------------------------------------------------------------------------

USE src_lheating,             ONLY :  &
    get_gs_lheating            ! storage of grid scale latent heating for lhn
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================

REAL (KIND=wp), PARAMETER :: &
  eps_div = repsilon ! small, precision-dependent value to be used in divisions
                     ! to avoid division by zero, e.g. a/MAX(b,eps_div)

!==============================================================================
! Module procedures
!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure "organize_gscp" in "gscp" to select the parameterization
!  scheme for grid-scale clouds and precipitation    
!------------------------------------------------------------------------------

!SUBROUTINE organize_gscp       

!------------------------------------------------------------------------------
!
! Description:
!   This module selects the parameterization scheme to calculate the effects
!   of grid-scale clouds and precipitation. Depending on the namelist input
!   parameter "itype_gscp", the following routines may be called
!
!   With prognostic precipiation (only this is possible now)
!     For diagnostic initialization of qr and qs (ldiniprec=.TRUE)
!     (3) a cloud ice scheme (itype_gscp=3)
!         named "hydci" (with argument linionly)
!
!   The other microphysical subroutines kessler_pp, hydor_pp,
!   hydci_pp and hydci_pp_gr are called after
!   the Dynamics from the driving subroutine organize_physics.
!
!------------------------------------------------------------------------------

! Subroutine arguments: None
! --------------------
! Local parameters,sclars and fields: None
! ----------------------------------

! Select for the parameterization chosen by the user:

!------------------------------------------------------------------------------
! End of module procedure organize_gscp
!------------------------------------------------------------------------------

!END SUBROUTINE organize_gscp

!==============================================================================
!==============================================================================
!+ Module procedure "hydci" in "gscp" for computing effects of grid scale 
!  precipitation including cloud water, rain and snow
!------------------------------------------------------------------------------

SUBROUTINE hydci (linionly, yerrmsg, ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure calculates the rates of change of temperature,
!   cloud water, cloud ice  and water vapor due to cloud microphysical 
!   processe related to the formation of grid scale precipitation. The 
!   tendencies are added to the global tendency fields.
!   The precipitation fluxes of rain and snow are calculated diagnostically
!   by integrating the corresponding budget equations from the top of
!   the atmosphere down to the surface. The resulting precipitation fluxes
!   at the surface are stored on the corresponding global fields.
!
! Method:
!   The tendencies involving cloud water and cloud ice are calculated with  
!   a quasi-implicit scheme.
!   For the budget equations of rain and snow stationarity and horizontal
!   homogenity are assumed. With respect to the integration of this equations
!   the precipitation fluxes of rain and snow are defined at half levels.
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------
  LOGICAL                 , INTENT(IN)           ::  &
    linionly             ! only diagnostic initialisation of qr, qs

  CHARACTER(LEN=*)        , INTENT(OUT)          ::  &
    yerrmsg

  INTEGER(KIND=iintegers) , INTENT(OUT)          ::  &
    ierrstat

! Local parameters:
! ----------------
  REAL    (KIND=wp   ),     PARAMETER ::  &
    ! Parameters for autoconversion of cloud water and cloud ice 
    zccau = 4.0E-4_wp, & ! autoconversion coefficient (cloud water to rain)
    zciau = 1.0E-3_wp, & ! autoconversion coefficient (cloud ice   to snow)

    ! Other process coefficients.
    ! These coefficients have been precalculated on the basis of the
    ! following basic parameters:
    ! N0R = 8.0 E6 : Parameter in the size distrubution function for rain
    ! N0S = 8.0 E5 : Parameter in the size distrubution function for snow
    ! ECR = 0.8    : Collection efficiency for rain collecting cloud water
    ! ECS = 0.9    : Collection efficiency for snow collecting cloud water
    ! EIS = 0.5    : Collection efficiency for snow collecting cloud ice
    ! EIR = 0.8    : Collection efficiency for rain collecting cloud ice
    ! V0R = 130.0  : Factor in the terminal velocity for raindrops
    !                VTR(D) = V0R*D**(1/2)
    ! V0S = 4.9    : Factor in the terminal velocity for snow paricles 
    !                VTS(DS) = V0S*DS**(1/4)
    ! AMS = 0.038  : Formfactor in the mass-size relation of snowparticles
    !                m(DS) = AMS*DS**2
    ! AMI = 0.038  : Formfactor in the mass-size relation of cloud ice crystals
    !                m(DI) = AMI*DS**3
    ! ETA  = 1.75 E-5 : Viscosity of air    
    ! DIFF = 2.22 E-5 : Molecular diffusion coefficient for water vapour
    ! LHEAT= 2.40 E-2 : Heat conductivity of air
    ! RHON = 1.0      : Density of air 
    ! RHOW = 1.0 E3   : Density of water 
    ! AR  = PI*1.0E3*N0R
    ! BR  = v0r*Gamma(4.5)/6
    ! AS  = 2*N0S*AMS
    ! BS  = N0S*V0S*AMS*Gamma(3.25)
    zcrim = 18.6_wp   , & ! 0.25*PI*ECS/AMS
    zcagg = 10.3_wp   , & ! 0.25*PI*EIS/AMS
    zcac  = 0.24_wp   , & ! 3.0*ECR*((ZAR*ZBR)**(2/9)/(7.0*E3)
    zcicri= 0.24_wp   , & ! 3.0*EIR*((ZAR*ZBR)**(2/9)/(7.0*E3)
    zcrcri= 3.20E-5_wp, & ! 0.25*PI*EIR*5.5*4.5*(ar*br)**(-4/9)
    zcev  = 1.0E-3_wp , & ! 2*PI*DIFF*HW*N0R*(AR*BR)**(-4/9)
    zbev  = 5.9_wp    , & ! 0.26*sqrt(rho*v0r/eta)*gam(2.75)*(AR*BR)**(-1/6)
    zcsdep= 1.8E-2_wp , & ! 4*DIFF*HW*N0S*(BS**(-8/13))
    zbsdep= 12.3_wp   , & ! 0.26*sqrt(rho*v0s/eta)*gam(21/8)*(BS)**(-5/26)
    zcsmel= 8.43E-5_wp, & ! 4*N0S*(BS**(-8/13))*xlheat/rho/lh_f
    zbsmel= 12.05_wp  , & ! 0.26*sqrt(rho*v0s/eta)*gam(21/8)*(BS)**(-5/26)
    zasmel= 2.31E3_wp , & ! diff*lh_v/lh_f*rho*lh_f/xlheat
    zcidep= 1.3E-5_wp , & ! 4*diff*hw/(ami**(1/3))
    zcrfrz= 3.75E-2_wp, & ! coefficient for raindrop freezing

    ! Additional parameters
    zfqrp = 0.105_wp  , & ! ar**(1/9)*br**(-8/9), to calculate qr from pr
    zfqsp = 0.430_wp  , & ! as*bs**(-12/13.), to calculate qs from ps
    zthn  = 236.15_wp , & ! temperature for hom. freezing of cloud water
    ztrfrz= 271.15_wp , & ! threshold temperature for heterogeneous
                              ! freezing of raindrops
    zmi0  = 1.0E-12_wp, & ! initial crystal mass for cloud ice nucleation
    zmimax= 1.0E-9_wp , & ! maximum mass of cloud ice crystals   
    zmsmin= 3.0E-9_wp , & ! initial mass of snow crystals        

    ! Constant exponents in the transfer rate equations
    x12o13 = 12.0_wp/13.0_wp  ,  x1o6   =  1.0_wp/ 6.0_wp, &
    x4o9   =  4.0_wp/ 9.0_wp  ,  x7o9   =  7.0_wp/ 9.0_wp, &
    x8o13  =  8.0_wp/13.0_wp  ,  x5o26  =  5.0_wp/26.0_wp, &
    x8o9   =  8.0_wp/ 9.0_wp  ,  x13o9  = 13.0_wp/ 9.0_wp

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    k, i, j, ij,       & ! loop indees
    ncc   ,            & ! number of cloud/precipitation points
    nx                   ! time-level for computation

  REAL    (KIND=wp   )     ::  &
    fpvsw, fpvsi, fqvs,& ! name of statement functions
    fxna ,             & ! statement function for ice crystal number
    ztx  , zpv  , zpx ,& ! dummy arguments for statement functions
    znimax,            & ! maximum number of cloud ice crystals
    zpvsw0,            & ! sat.vap. pressure at melting temperature
    zqvsw0,            & ! sat.specific humidity at melting temperature
    zdtc  , zdti,      & ! timestep for integration (water / ice )
    zdtrc ,            & ! reciprocal of timestep for integration of water
    zdtri ,            & ! reciprocal of timestep for integration of ice
    zscau , zscac  , zscrim , zscshe, zsnuc , & ! local values of the
    zsiau , zsagg  , zsidep , zsicri, zsrcri, & ! transfer rates
    zsdau , zssdep , zssmelt, zsimelt, zsev,  & ! defined below
    zsrfrz, zcorr,                            & !
    zscsum, zscmax,    & ! terms for limiting  total cloud water depletion
    zsrsum, zsrmax,    & ! terms for limiting  total rain water depletion
    znin,              & ! number of cloud ice crystals at nucleation
    znid,              & ! number of cloud ice crystals for deposition
    zmi ,              & ! mass of a cloud ice crystal
    zsvidep, zsvisub,  & ! deposition, sublimation of cloud ice
    zsimax , zsisum ,  & ! terms for limiting total
    zsvmax,            & ! cloud ice depletion
    zrhodz,            & ! rho*dz
    zqvsi, zqvsw,      & ! sat. specitic humidity at ice and water saturation
    zztau, zxfac, zprx, zpsx, zx1, zx2, &! help variables
    zrlogprr, zrlogprs, zeln7o9prr, zeln13o9prr, zeln5o26prs, zeln8o13prs

  LOGICAL :: &
    lmassb_tot,        & ! switch for doing a mass-budget calculation
    lcwexist,          & ! switch for existence of cloud water
    lciexist,          & ! switch for existence of cloud ice  
    lprexist,          & ! switch for existence of rain       
    lpsexist,          & ! switch for existence of snow       
    ldini                ! only diagnostic initialisation of qr, qs

! Local (automatic) arrays:
! -------------------------
  REAL    (KIND=wp   )     ::  &
    scau  (ie,je) ,& ! transfer rate due to autoconversion of cloud water
    scac  (ie,je) ,& ! transfer rate due to accretion of cloud water
    snuc  (ie,je) ,& ! transfer rate due nucleation of cloud ice    
    scfrz (ie,je) ,& ! transfer rate due homogeneous freezing of cloud water
    simelt(ie,je) ,& ! transfer rate due melting of cloud ice
    sidep (ie,je) ,& ! transfer rate due depositional growth of cloud ice
    ssdep (ie,je) ,& ! transfer rate due depositional growth of snow      
    sdau  (ie,je) ,& ! transfer rate due depositional cloud ice autoconversion
    srim  (ie,je) ,& ! transfer rate due riming of snow
    sshed (ie,je) ,& ! transfer rate due shedding
    sicri (ie,je) ,& ! transfer rate due cloud ice collection by rain (sink qi)
    srcri (ie,je) ,& ! transfer rate due cloud ice collection by rain (sink qr)
    sagg  (ie,je) ,& ! transfer rate due aggregation of snow and cloud ice
    siau  (ie,je) ,& ! transfer rate due autoconversion of cloud ice      
    ssmelt(ie,je) ,& ! transfer rate due melting of snow
    sev   (ie,je) ,& ! transfer rate due evaporation of rain
    srfrz (ie,je) ,& ! transfer rate due to rainwater freezing
    zpr   (ie,je) ,& ! precipitation flux of rain entering from above
    zps   (ie,je) ,& ! precipitation flux of snow entering from above
    zqvt  (ie,je) ,& ! layer tendency of specific humidity
    zqct  (ie,je) ,& ! layer tendency of cloud water       
    zqit  (ie,je) ,& ! layer tendency of cloud ice         
    zqrt  (ie,je) ,& ! layer tendency of rain              
    zqst  (ie,je) ,& ! layer tendency of snow              

    zqv   (ieje ) ,& ! specific humidity    
    zqc   (ieje ) ,& ! cloud water content
    zqi   (ieje ) ,& ! cloud ice content
    zt    (ieje ) ,& ! temperature
    zph   (ieje ) ,& ! pressure
    zrho  (ieje ) ,& ! air density
    zprr  (ieje ) ,& ! precipitation flux of rain entering from above
    zprs  (ieje )    ! precipitation flux of snow entering from above

  REAL    (KIND=wp   )     ::  &
    zconvqit(ke)  ,& !
    zconvqct(ke)  ,& !
    zconvqvt(ke)  ,& !
    zconvqrt(ke)  ,& !
    zconvqst(ke)  ,& !
    zconvtot(ke)  ,& ! total of conversion rates in a layer k
    zconptot(ke)  ,& ! total of conversion rates prod. precip. in a layer k
    zprecipp(ke)  ,& ! total of precipitation divergence in a layer k
    scau_sum  (ke),& ! = SUM(scau)
    scac_sum  (ke),& ! = SUM(scac)
    snuc_sum  (ke),& ! = SUM(snuc)
    scfrz_sum (ke),& ! = SUM(scfrz)
    simelt_sum(ke),& ! = SUM(simelt)
    sidep_sum (ke),& ! = SUM(sidep)
    ssdep_sum (ke),& ! = SUM(ssdep)
    sdau_sum  (ke),& ! = SUM(sdau)
    srim_sum  (ke),& ! = SUM(srim)
    sshed_sum (ke),& ! = SUM(sshed)
    sicri_sum (ke),& ! = SUM(sicri)
    srcri_sum (ke),& ! = SUM(srcri)
    sagg_sum  (ke),& ! = SUM(sagg )
    siau_sum  (ke),& ! = SUM(siau)
    ssmelt_sum(ke),& ! = SUM(ssmelt)
    sev_sum   (ke)   ! = SUM(sev)


  INTEGER (KIND=iintegers) ::  &
    ic (ieje),         & !  i-index of cloud/precipitation points
    jc (ieje)            !  j-index of cloud/precipitation points
 

  INTEGER (KIND=iintegers) :: izerror
  CHARACTER (LEN=255)      :: yzerrmsg

! Tracer pointers
!----------------
  REAL (KIND=wp),     POINTER :: &
    qv_nx   (:,:,:)=> NULL(),   & ! QV at nx
    qv      (:,:,:)=> NULL(),   & ! QV at nnew
    qv_tens (:,:,:)=> NULL(),   & ! QV tendency
    qc_nx   (:,:,:)=> NULL(),   & ! QC at nx
    qc      (:,:,:)=> NULL(),   & ! QC at nnew
    qi_nx   (:,:,:)=> NULL(),   & ! QI at nx
    qr_now  (:,:,:)=> NULL(),   & ! QR at nnow
    qr_old  (:,:,:)=> NULL(),   & ! QR at nold
    qs_now  (:,:,:)=> NULL(),   & ! QS at nnow
    qs_old  (:,:,:)=> NULL()      ! QS at nold

  CHARACTER (LEN=25)       :: yzroutine = 'hydci'

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine hydci
!------------------------------------------------------------------------------

! Statement functions
! -------------------

! saturation vapour pressure over water (fpvsw), over ice (fpvsi)
! and specific humidity at vapour saturation (fqvs)
  fpvsw(ztx)     = b1*EXP( b2w*(ztx-b3)/(ztx-b4w) )
  fpvsi(ztx)     = b1*EXP( b2i*(ztx-b3)/(ztx-b4i) )
  fqvs (zpv,zpx) = rdv*zpv/( zpx - o_m_rdv*zpv )

! Number of activate ice crystals;  ztx is temperature
  fxna(ztx)   = 1.0E2_wp * EXP(0.2_wp * (t0_melt - ztx))

!------------------------------------------------------------------------------
!  Section 1: Initial setting of local and global variables
!------------------------------------------------------------------------------

  yerrmsg = '                                  '
  ierrstat = 0_iintegers

  ! This statement in principle is obsolet now, because hydci can only be 
  ! called with linionly=.TRUE. now that only prognostic precipitation is done
  IF (linionly) THEN
    ldini = linionly
  ELSE
    ldini = .FALSE.
    ierrstat = 1
    yerrmsg  = 'hydci can only be called with linionly=.TRUE.'
    RETURN
  END IF

! Mass budget calculation on or off
  lmassb_tot = .FALSE.

! Some constant coefficients
  znimax = fxna(zthn) ! Maximum number of cloud ice crystals
  zpvsw0 = fpvsw(t0_melt)  ! sat. vap. pressure for t = t0_melt

! Delete mixing ratio of rain and snow and precipitation fluxes from previous
! timestep
  qrs   (:,:,:) = 0.0_wp
  prr_gsp (:,:) = 0.0_wp
  prs_gsp (:,:) = 0.0_wp

! select timelevel and timestep for calculations
  IF ( l2tls ) THEN
    nx    = nnew
    zdtc  = dt 
    zdti  = dt 
    zdtrc = 1.0_wp / dt
    zdtri = 1.0_wp / dt
  ELSE
    nx    = nnew
    IF (ldini) nx = nnow
    zdtc  = dt2
    zdti  = dt
    zdtrc = 1.0_wp / dt2
    zdtri = 1.0_wp / dt 
  ENDIF

  ! retrieve the required microphysics tracers (at corresponding timelevel)
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnew, ptr = qv, ptr_tens = qv_tens)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nx, ptr = qv_nx)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nx, ptr = qc_nx)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nnew, ptr = qc)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qi, ptr_tlev = nx, ptr = qi_nx)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  IF ( ldini ) THEN
    CALL trcr_get(izerror, idt_qr, ptr_tlev = nnow, ptr = qr_now)
    IF (izerror /= 0) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
    IF (.NOT. l2tls) THEN
      CALL trcr_get(izerror, idt_qr, ptr_tlev = nold, ptr = qr_old)
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF
    ENDIF
    CALL trcr_get(izerror, idt_qs, ptr_tlev = nnow, ptr = qs_now)
    IF (izerror /= 0) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
    IF (.NOT. l2tls) THEN
      CALL trcr_get(izerror, idt_qs, ptr_tlev = nold, ptr = qs_old)
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF
    ENDIF
  ENDIF

! *********************************************************************
! Loop from the top of the model domain to the surface to calculate the
! transfer rates and the precipitation fluxes.     
! *********************************************************************

DO  k = 1, ke
 
  ! Check for existence of cloud water, cloud ice, rain and snow
  lcwexist = MAXVAL(qc_nx(istartpar:iendpar,jstartpar:jendpar,k)) > 0.0_wp
  lciexist = MAXVAL(qi_nx(istartpar:iendpar,jstartpar:jendpar,k)) > 0.0_wp 
  IF (k > 1 ) THEN
    lprexist = MAXVAL(prr_gsp(istartpar:iendpar,jstartpar:jendpar)) > 0.0_wp 
    lpsexist = MAXVAL(prs_gsp(istartpar:iendpar,jstartpar:jendpar)) > 0.0_wp 
  ELSE
    lprexist = .FALSE.
    lpsexist = .FALSE.
  ENDIF
 
  ! Initialize the conversion rates with zeros in every layer
  ! ---------------------------------------------------------
  scau  (:,:) = 0.0_wp
  scac  (:,:) = 0.0_wp
  snuc  (:,:) = 0.0_wp
  scfrz (:,:) = 0.0_wp
  simelt(:,:) = 0.0_wp
  sidep (:,:) = 0.0_wp
  ssdep (:,:) = 0.0_wp
  sdau  (:,:) = 0.0_wp
  srim  (:,:) = 0.0_wp
  sshed (:,:) = 0.0_wp
  sicri (:,:) = 0.0_wp
  srcri (:,:) = 0.0_wp
  sagg  (:,:) = 0.0_wp
  siau  (:,:) = 0.0_wp
  ssmelt(:,:) = 0.0_wp
  sev   (:,:) = 0.0_wp
  srfrz (:,:) = 0.0_wp
 
  ! Deposition nucleation for low temperatures below a threshold
  IF ( MINVAL(t (istartpar:iendpar,jstartpar:jendpar,k,nx)) < 248.15_wp  &
       .AND.                                                                 &
       MAXVAL(qv_nx(istartpar:iendpar,jstartpar:jendpar,k)) > 8.E-6_wp ) THEN
    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar
        zx1 = p0(i,j,k) + pp(i,j,k,nx)
        zqvsi     = fqvs( fpvsi(t(i,j,k,nx)), zx1 )
        IF( t(i,j,k,nx) < 248.15_wp .AND. qv_nx(i,j,k) > zqvsi .AND. &
           qv_nx(i,j,k) > 8.E-6_wp  .AND. qi_nx(i,j,k) <= 0.0_wp) THEN
        znin  = MIN( fxna(t(i,j,k,nx)), znimax )
        zsnuc = zmi0 / rho(i,j,k) * znin * zdtri
        snuc(i,j) = zsnuc
      ENDIF
      ENDDO
    ENDDO
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 2: Search for cloudy grid points with cloud water and
  !            calculation of the conversion rates involving qc
  !----------------------------------------------------------------------------
  IF (lcwexist) THEN
  ncc = 0
  DO j = jstartpar, jendpar
    DO i = istartpar, iendpar
      IF ( qc_nx(i,j,k) > 0.0_wp ) THEN
        ncc = ncc + 1
        ic(ncc) = i
        jc(ncc) = j
        zqc (ncc) = qc_nx(i,j,k)
        zqi (ncc) = qi_nx(i,j,k)
        zt  (ncc) = t (i,j,k,nx)
        zrho(ncc) = rho(i,j,k)
        zprr(ncc) = prr_gsp(i,j)
        zprs(ncc) = prs_gsp(i,j)
      ENDIF
    ENDDO
  ENDDO

  IF ( ncc > 0 ) THEN ! conversion rates involving cloud water
    DO ij = 1, ncc
      zscau  = zccau * MAX( zqc(ij) - qc0, 0.0_wp )
      IF (zprr(ij) == 0.0_wp) THEN
        zscac  = 0.0_wp
      ELSE
        zscac  = zcac  * EXP(x7o9 * LOG(zprr(ij))) * zqc(ij)
      ENDIF
      zscrim = zcrim * zqc(ij) * zprs(ij)
      zscshe = 0.0_wp
      IF( zt(ij) >= t0_melt ) THEN
        zscshe = zscrim
        zscrim = 0.0_wp
      ENDIF
      ! Check for maximum depletion of cloud water and adjust the
      ! transfer rates accordingly
      zscmax = zqc(ij)*zdtrc 
      zscsum = zscau + zscac + zscrim + zscshe 
      zcorr  = zscmax / MAX( zscmax, zscsum )
      IF( zt(ij) <= zthn ) THEN
        scfrz(ic(ij),jc(ij)) = zscmax
      ELSE
        scau (ic(ij),jc(ij)) = zcorr*zscau
        scac (ic(ij),jc(ij)) = zcorr*zscac
        srim (ic(ij),jc(ij)) = zcorr*zscrim
        sshed(ic(ij),jc(ij)) = zcorr*zscshe 
      ENDIF
      ! Calculation of heterogeneous nucleation of cloud ice.
      ! This is done in this section, because we require water saturation
      ! for this process (i.e. the existence of cloud water) to exist.
      ! Hetrogeneous nucleation is assumed to occur only when no
      ! cloud ice is present and the temperature is below a nucleation
      ! threshold.
      IF( zt(ij) <= 267.15_wp .AND. zqi(ij) <= 0.0_wp ) THEN
        znin  = MIN( fxna(zt(ij)), znimax )
        zsnuc = zmi0 / zrho(ij) * znin * zdtri
        snuc(ic(ij),jc(ij)) = zsnuc
      ENDIF
      ! Calculation of in-cloud rainwater freezing
      IF ( zt(ij) < ztrfrz ) THEN
        zsrfrz = zcrfrz*SQRT( ( (ztrfrz-zt(ij))*zprr(ij) )**3 )
        srfrz(ic(ij),jc(ij)) = zsrfrz
      ENDIF
     
    ENDDO
  ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 3: Search for cold grid points with cloud ice and/or snow and
  !            calculation of the conversion rates involvin qi and ps
  !----------------------------------------------------------------------------
  IF (lciexist.OR.lpsexist) THEN
  ncc = 0
  DO j = jstartpar, jendpar
    DO i = istartpar, iendpar
      IF ( qi_nx(i,j,k) + prs_gsp(i,j) > 0.0_wp .AND. T(i,j,k,nx) <= t0_melt ) THEN
        ncc = ncc + 1
        ic(ncc) = i
        jc(ncc) = j
        zt  (ncc) = t (i,j,k,nx)
        zph (ncc) = p0(i,j,k) + pp(i,j,k,nx)
        zqv (ncc) = qv_nx(i,j,k)
        zqi (ncc) = qi_nx(i,j,k)
        zrho(ncc) = rho(i,j,k)
        zprr(ncc) = prr_gsp(i,j)
        zprs(ncc) = prs_gsp(i,j)
      ENDIF
    ENDDO
  ENDDO

  IF ( ncc > 0 ) THEN  ! Calculation of the conversion rates
    DO ij = 1, ncc

      IF (zprr(ij) /= 0.0_wp) THEN
        zrlogprr    = LOG (zprr(ij))
        zeln7o9prr  = EXP (x7o9  * zrlogprr)
        zeln13o9prr = EXP (x13o9 * zrlogprr)
      ELSE
        zeln7o9prr  = 0.0_wp
        zeln13o9prr = 0.0_wp
      ENDIF
      IF (zprs(ij) /= 0.0_wp) THEN
        zrlogprs    = LOG (zprs(ij))
        zeln5o26prs = EXP (x5o26 * zrlogprs)
        zeln8o13prs = EXP (x8o13 * zrlogprs)
      ELSE
        zeln5o26prs = 0.0_wp
        zeln8o13prs = 0.0_wp
      ENDIF

      zqvsi     = fqvs( fpvsi(zt(ij)), zph(ij) )
      znin      = MIN( fxna(zt(ij)), znimax )
      zmi       = MIN( zrho(ij)*zqi(ij)/znin, zmimax )
      zmi       = MAX( zmi0, zmi )
      ! maximum of cloud ice deposition/sublimation rate (zsvmax)
      ! is limited to achieve ice saturation zqvsip at the next time step
      ! (counted negative, i.e ! qv(n+1) = q(n)-zsvmax*dt ).
      !zpvsidt   = lh_s/(r_v*zt(ij)**2)
      !zdqvidt   = zpvsidt*( 1.0_wp + rvd_m_o*zqvsi ) * zqvsi 
      !zqvsip    = zqvsi + zdqvidt*lh_s*cpdr * ( zqv(ij) - zqvsi ) &
      !            / (1. + lh_s*cpdr*zdqvidt)
      !zsvmax    = (zqv(ij) - zqvsip) * zdtrc 
      zsvmax    = (zqv(ij) - zqvsi) * zdtrc 
      zsiau     = zciau * MAX( zqi(ij)- qi0, 0.0_wp )
      zsagg     = zcagg * zqi(ij) * zprs(ij)
      znid      = zrho(ij)*zqi(ij)/zmi
      zsidep    = zcidep * znid * EXP(0.33_wp*LOG(zmi)) * ( zqv(ij) - zqvsi )
      zsvidep   = 0.0_wp
      zsvisub   = 0.0_wp
      zsimax    = zqi(ij)*zdtri 
      IF( zsidep > 0.0_wp ) THEN
        zsvidep = MIN(   zsidep, zsvmax )
      ELSEIF (zsidep < 0.0_wp ) THEN
        zsvisub = - MAX(-zsimax, zsvmax )
      END IF
      zztau     = 1.5_wp*( EXP(0.66_wp*LOG(zmsmin/zmi)) - 1.0_wp)
      zsdau     = zsvidep/zztau
      zsicri    = zcicri * zqi(ij) * zeln7o9prr
      zsrcri    = zcrcri * (zqi(ij)/zmi) * zeln13o9prr
      zxfac     = 1.0_wp + zbsdep* zeln5o26prs
      zssdep    = zcsdep*zxfac*( zqv(ij) - zqvsi )* zeln8o13prs
      ! Check for maximal depletion of cloud ice
      ! No check is done for depositional autoconversion because
      ! this is a always a fraction of the gain rate due to
      ! deposition (i.e the sum of this rates is always positive)
      zsisum = zsiau + zsagg + zsicri + zsvisub
      zcorr  = 0.0_wp
      IF( zsimax > 0.0_wp ) zcorr  = zsimax / MAX( zsimax, zsisum )
      sidep(ic(ij),jc(ij))  = zsvidep - zcorr*zsvisub
      sdau (ic(ij),jc(ij))  = zsdau
      siau (ic(ij),jc(ij))  = zcorr*zsiau
      sagg (ic(ij),jc(ij))  = zcorr*zsagg
      ssdep(ic(ij),jc(ij))  = zssdep
      srcri(ic(ij),jc(ij))  = zsrcri
      sicri(ic(ij),jc(ij))  = zcorr*zsicri
    ENDDO
  ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 4: Search for warm grid points with cloud ice and/or snow and
  !            calculation of the melting rates of qi and ps
  !----------------------------------------------------------------------------
  IF (lciexist.OR.lpsexist) THEN
  ncc = 0
  DO j = jstartpar, jendpar
    DO i = istartpar, iendpar
      IF ( qi_nx(i,j,k) + prs_gsp(i,j) > 0.0_wp .AND. t(i,j,k,nx) > t0_melt ) THEN
        ncc = ncc + 1
        ic(ncc) = i
        jc(ncc) = j
        zt  (ncc) = t (i,j,k,nx)
        zph (ncc) = p0(i,j,k) + pp(i,j,k,nx)
        zqv (ncc) = qv_nx(i,j,k)
        zqi (ncc) = qi_nx(i,j,k)
        zprs(ncc) = prs_gsp(i,j)
      ENDIF
    ENDDO
  ENDDO
  IF ( ncc > 0 ) THEN ! Calculation of the conversion rates
    DO ij = 1, ncc
      zsimelt  = zqi(ij)*zdtri 
      zqvsw0    = fqvs( zpvsw0, zph(ij))
      zx1      = (zt(ij) - t0_melt) + zasmel*(zqv(ij) - zqvsw0)

      IF (zprs(ij) /= 0.0_wp) THEN
        zrlogprs  = LOG (zprs(ij))
        zx2       = 1.0_wp + zbsmel*EXP(x5o26*zrlogprs)
        zssmelt  = zcsmel * zx1 * zx2 * EXP(x8o13*zrlogprs)
      ELSE
        zssmelt  = 0.0_wp
      ENDIF
      ssmelt(ic(ij),jc(ij)) = MAX( zssmelt, 0.0_wp )
      simelt(ic(ij),jc(ij)) = zsimelt
    ENDDO
  ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 5: Search for grid points with rain in subsaturated areas
  !            and calculation of the evaporation rate of rain
  !----------------------------------------------------------------------------
  IF (lprexist) THEN
  ncc = 0
  DO j = jstartpar, jendpar
    DO i = istartpar, iendpar
      IF ( prr_gsp(i,j) > 0.0_wp .AND. qc_nx(i,j,k) <= 0.0_wp ) THEN
        ncc = ncc + 1
        ic(ncc) = i
        jc(ncc) = j
        zprr(ncc) = prr_gsp(i,j)
        zqv (ncc) = qv_nx(i,j,k)
        zph (ncc) = p0(i,j,k) + pp(i,j,k,nx)
        zt  (ncc) = t (i,j,k,nx)
      ENDIF
    ENDDO
  ENDDO
  IF ( ncc > 0 ) THEN ! Calculate evaporation of rain
    DO ij = 1, ncc
      zqvsw  = fqvs( fpvsw(zt(ij)), zph(ij) )

      IF (zprr(ij) /= 0.0_wp) THEN
        zrlogprr  = LOG (zprr(ij))
        zx1     = 1.0_wp + zbev* EXP(x1o6*zrlogprr)
        zsev = zcev*zx1*(zqvsw - zqv(ij))* EXP(x4o9*zrlogprr)
      ELSE
        zsev = 0.0_wp
      ENDIF

      sev(ic(ij),jc(ij)) = MAX( zsev, 0.0_wp )
      ! Calculation of below-cloud rainwater freezing
      IF ( zt(ij) < ztrfrz ) THEN
        zsrfrz = zcrfrz*SQRT( ( (ztrfrz-zt(ij))*zprr(ij) )**3 )
        srfrz(ic(ij),jc(ij)) = zsrfrz
      ENDIF
    ENDDO
  ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 6: Calculate the total tendencies of the prognostic variables.
  !            Perform the vertical integration of the precipitation fluxes.
  !            Update the prognostic variables in the interior domain.
  !----------------------------------------------------------------------------
  DO j = jstartpar, jendpar
    DO i = istartpar, iendpar
      zrhodz = rho(i,j,k)*(hhl(i,j,k)-hhl(i,j,k+1))
      zsrmax = prr_gsp(i,j)/zrhodz
      zsrsum = sev(i,j) + srfrz(i,j) + srcri(i,j)
      zcorr  = 1.0_wp
      IF(zsrsum > 0) THEN
        zcorr  = zsrmax / MAX( zsrmax, zsrsum )
      ENDIF
      sev  (i,j) = zcorr*sev  (i,j)
      srfrz(i,j) = zcorr*srfrz(i,j)
      srcri(i,j) = zcorr*srcri(i,j)

      ssmelt(i,j) = MIN(ssmelt(i,j), prs_gsp(i,j)/zrhodz)
      IF (ssdep(i,j) < 0.0_wp ) THEN
        ssdep(i,j) = MAX(ssdep(i,j), - prs_gsp(i,j)/zrhodz)
      ENDIF
      zqvt(i,j) =   sev  (i,j) - sidep(i,j) - ssdep (i,j) - snuc(i,j)
      zqct(i,j) = - scau (i,j) - scfrz(i,j) + simelt(i,j) - scac(i,j)   &
                  - sshed(i,j) - srim (i,j)
      zqit(i,j) =   snuc (i,j) + scfrz(i,j) - simelt(i,j) - sicri (i,j) &
                  + sidep(i,j) - sdau (i,j) - sagg  (i,j) - siau  (i,j)
      zqrt(i,j) =   scau (i,j) + sshed(i,j) + scac  (i,j) + ssmelt(i,j) &
                  - sev  (i,j) - srcri(i,j) - srfrz (i,j)
      zqst(i,j) =   siau (i,j) + sdau (i,j) + sagg  (i,j) - ssmelt(i,j) &
                  + sicri(i,j) + srcri(i,j) + srim  (i,j) + ssdep (i,j) &
                  + srfrz(i,j)
      zpr (i,j) = prr_gsp(i,j)
      zps (i,j) = prs_gsp(i,j)
      prr_gsp(i,j) = prr_gsp(i,j) + zqrt(i,j)*zrhodz
      prs_gsp(i,j) = prs_gsp(i,j) + zqst(i,j)*zrhodz
    ENDDO
  ENDDO

  DO j = jstartpar, jendpar
    DO i = istartpar, iendpar
      IF ( prr_gsp(i,j) < 1.0E-15_wp ) prr_gsp(i,j) = 0.0_wp
      IF ( prs_gsp(i,j) < 1.0E-15_wp ) prs_gsp(i,j) = 0.0_wp
    ENDDO
  ENDDO

  IF (ldini) THEN
    ! calculate the diagnostic rain amount without feedback to the
    ! other variables

    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar
        zprx = 0.5_wp*( zpr(i,j) + prr_gsp(i,j) )
        zpsx = 0.5_wp*( zps(i,j) + prs_gsp(i,j) )
        IF ( zprx > 0.0_wp ) THEN
          qrs(i,j,k) = 1.0_wp/rho(i,j,k)*zfqrp*EXP(x8o9*LOG(zprx))
          IF (ldini) THEN
            qr_now(i,j,k) = qrs(i,j,k)
            IF (.NOT.l2tls) qr_old(i,j,k) = qrs(i,j,k)
          ENDIF
        ENDIF
        IF ( zpsx > 0.0_wp ) THEN
          qrs(i,j,k) = qrs(i,j,k)+1.0_wp/rho(i,j,k)*zfqsp*EXP(x12o13*LOG(zpsx))
          IF (ldini) THEN
            qs_now(i,j,k) = 1.0_wp/rho(i,j,k)*zfqsp*EXP(x12o13*LOG(zpsx))
            IF (.NOT.l2tls) qs_old(i,j,k) = qs_now(i,j,k)
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    DO j = jstart, jend
      DO i = istart, iend
        qrs(i,j,k)   = qrs(i,j,k) + qi_nx(i,j,k)
      ENDDO
    ENDDO

    ! Do a mass-budget calculation for the total domain if required
    IF (lmassb_tot) THEN
      zconvqit(k) = 0.0_wp
      zconvqct(k) = 0.0_wp
      zconvqvt(k) = 0.0_wp
      zconvqrt(k) = 0.0_wp
      zconvqst(k) = 0.0_wp
      zprecipp(k) = 0.0_wp
      DO j = jstartpar, jendpar
        DO i = istartpar, iendpar
          zrhodz      = rho(i,j,k)*(hhl(i,j,k)-hhl(i,j,k+1))
          zconvqit(k) = zconvqit(k) + zqit(i,j)*zrhodz
          zconvqct(k) = zconvqct(k) + zqct(i,j)*zrhodz
          zconvqvt(k) = zconvqvt(k) + zqvt(i,j)*zrhodz
          zconvqrt(k) = zconvqrt(k) + zqrt(i,j)*zrhodz
          zconvqst(k) = zconvqst(k) + zqst(i,j)*zrhodz
          zprecipp(k) = zprecipp(k) + ( prr_gsp(i,j) + prs_gsp(i,j) &
                                      - zpr(i,j) - zps(i,j) )
        ENDDO
      ENDDO
      zconptot(k) = zconvqit(k)+zconvqct(k)+zconvqvt(k)
      zconvtot(k) = zconptot(k)+zconvqrt(k)+zconvqst(k)
      ! sums of the conversion rates
      scau_sum  (k) = SUM(scau)
      scac_sum  (k) = SUM(scac)
      snuc_sum  (k) = SUM(snuc)
      scfrz_sum (k) = SUM(scfrz)
      simelt_sum(k) = SUM(simelt)
      sidep_sum (k) = SUM(sidep)
      ssdep_sum (k) = SUM(ssdep)
      sdau_sum  (k) = SUM(sdau)
      srim_sum  (k) = SUM(srim)
      sshed_sum (k) = SUM(sshed)
      sicri_sum (k) = SUM(sicri)
      srcri_sum (k) = SUM(srcri)
      sagg_sum  (k) = SUM(sagg )
      siau_sum  (k) = SUM(siau)
      ssmelt_sum(k) = SUM(ssmelt)
      sev_sum   (k) = SUM(sev)
    ENDIF

  ENDIF  ! (ldini)
ENDDO

!------------------------------------------------------------------------------
! End of module procedure hydci   
!------------------------------------------------------------------------------

END SUBROUTINE hydci  

!==============================================================================
!==============================================================================
!+ Module procedure "kessler_pp" in "gscp" for computing effects of grid scale 
!  precipitation including cloud water and rain in 
!  context with the Leapfrog and the Runge-Kutta time-integration
!------------------------------------------------------------------------------

SUBROUTINE kessler_pp

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure calculates the rates of change of temperature, 
!   cloud water, water vapor, and rain due to cloud microphysical processes 
!   related to the formation of grid scale precipitation. The variables will 
!   be updated in this subroutine.
!   The precipitation fluxes at the 
!   surface are stored on the corresponding global fields.
!   This subroutine relies on conversion rates used in the subroutine 
!   hydorprog.
!
! Method:
!   The tendencies involving cloud water are calculated using a full implicit 
!   scheme whereas evaporation is computed explicitly.
!   Rain is a prognostic variable and the sedimentation term is 
!   computed implicitly.
!
!------------------------------------------------------------------------------
!
! Declarations:

! Subroutine arguments: None
! --------------------
!
! Local parameters:
! ----------------
  REAL    (KIND=wp   ),     PARAMETER ::  &
    ! basic constants of the parameterization scheme
    zaau   = 1.0_wp/1000.0_wp,         & ! coef. for autoconversion
    zaac   = 1.72_wp,                      & ! coef. for accretion (neu)
    zbev   = 9.1_wp,                       & ! coef. for drop ventilation

    ! constants for the process rates
    z1d8   = 1.0_wp/8.0_wp,               & !
    z7d8   = 7.0_wp/8.0_wp,               & !
    z3d16  = 3.0_wp/16.0_wp,              & !

    ! constants for sedimentation
    zvz0r = 12.63_wp,                         & !

    ! to avoid illegal operations in power expressions
    znull = 1.E-20_wp

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    k     ,            & ! loop index in vertical direction              
    i     ,            & ! loop index in x-direction              
    j     ,            & ! loop index in y-direction
    nu                   ! time level


  REAL    (KIND=wp   )     ::  &
    ztx   ,            & ! 
    zpx   ,            & !
    zgex  ,            & !  
    fspw  ,            & ! function for equilibrium vapour pressure over water
    fsqv  ,            & ! function for specific humidity at 
    fsa3  ,            & !
    zspw  ,            & ! equilibrium vapour pressure over water
    zsa3  ,            & !
    zc3   ,            & !
    zc1c  ,            & !
    zc1   ,            & !
    zx    ,            & !
    zsrmax

  REAL    (KIND=wp   )     ::  & 
    zphf  ,            & !  pressure in a k-layer
    zsqvw ,            & !  specific humidity at water saturation
    zqvts ,            & !  qv-tendency in a layer
    zqcts ,            & !  qc-tendency in a layer
    zqrts ,            & !  qr-tendency in a layer
    ztts  ,            & !  t -tendency in a layer
    zswra ,            & !  autoconversion rate
    zswrk ,            & !  accretion rate
    zsrd                 !  evaporation rate
 

  REAL    (KIND=wp   )     ::  &
    zvzr(ie,je),       & !
    zpres (ie,je) ,&     ! pressure
    zprvr(ie,je),      & !
    zpkr(ie,je),       & !
    zpkm1r(ie,je),     & !
    zdummy(ie,je,8),   & !
    zzar,              & !
    zdh,               & !
    zdtdh,             & !
    zqrk,              & !
    lnzqrk,            & !
    zdt,               & ! 
    zdtr,              & ! 
    zimr,              & !
    qcg,               & !
    tg,                & !
    qvg,               & !
    qrg,               & !
    rhog,              & !
    rhogr

  INTEGER (KIND=iintegers)    :: izerror
  CHARACTER (LEN=255)         :: yzerrmsg
    
! Tracer pointers
!----------------
  REAL (KIND=wp),     POINTER :: &
    qv  (:,:,:)=> NULL(),    & ! QV at nu
    qc  (:,:,:)=> NULL(),    & ! QC at nu
    qr  (:,:,:)=> NULL()       ! QR at nu
    
    
  CHARACTER(LEN=25) :: yzroutine = 'kessler_pp'

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine kessler_pp
!------------------------------------------------------------------------------

! Statement functions
! -------------------
 
fspw(ztx)     = b1*EXP( b2w*(ztx-b3)/(ztx-b4w) ) ! sat.vap. pressure over water
fsqv(zgex,zpx) = rdv*zgex/( zpx - o_m_rdv*zgex ) ! specific humidity at sat.
  
! coefficients for conversion rates
fsa3(ztx) = 3.86E-3_wp - 9.41E-5_wp*(ztx-t0_melt) 

!------------------------------------------------------------------------------
!  Section 1: Initial setting of local variables
!
!------------------------------------------------------------------------------


zdummy(:,:,:)=0.0_wp

! select timelevel and timestep for calculations
nu    = nnew
IF ( l2tls ) THEN
  zdt   = dt
ELSE
  zdt   = dt2
ENDIF
zdtr  = 1.0_wp / zdt


zpkr(:,:) = 0.0_wp
zprvr(:,:) = 0.0_wp
zvzr(:,:) = 0.0_wp

! retrieve the required microphysics tracers (at corresponding timelevel)
CALL trcr_get(izerror, idt_qv, ptr_tlev = nu, ptr = qv)
IF (izerror /= 0) THEN
  yzerrmsg = trcr_errorstr(izerror)
  CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
ENDIF
CALL trcr_get(izerror, idt_qc, ptr_tlev = nu, ptr = qc)
IF (izerror /= 0) THEN
  yzerrmsg = trcr_errorstr(izerror)
  CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
ENDIF
CALL trcr_get(izerror, idt_qr, ptr_tlev = nu, ptr = qr)
IF (izerror /= 0) THEN
  yzerrmsg = trcr_errorstr(izerror)
  CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
ENDIF

#ifdef NUDGING
  ! add part of latent heating calculated in subroutine kessler_pp to model latent
  ! heating field: subtract temperature from model latent heating field
    IF (llhn) &
       CALL get_gs_lheating ('add',1,ke)
#endif

loop_over_levels: DO k = 1, ke

  IF ( ldiabf_lh ) THEN
    ! initialize temperature increment due to latent heat
    tinc_lh(:,:,k) = tinc_lh(:,:,k) - t(:,:,k,nu)
  ENDIF

  !----------------------------------------------------------------------------
  !  Section 2: Test for clouds and precipitation in present layer.
  !             If no cloud or precipitation points have been found
  !             go to the next layer.
  !----------------------------------------------------------------------------


  loopj: DO j = jstartpar, jendpar
    loopi: DO i = istartpar, iendpar
      IF(qr (i,j,k) < 1.0e-15_wp) qr (i,j,k) = 0.0_wp
      IF(qc (i,j,k) < 1.0e-15_wp) qc (i,j,k) = 0.0_wp
    ENDDO loopi
  ENDDO loopj


  DO j = jstartpar, jendpar
    DO i = istartpar, iendpar

      qcg  = qc(i,j,k)
      qvg  = qv(i,j,k)
      qrg  = qr(i,j,k)
      tg   = t(i,j,k,nu)
      rhog = rho(i,j,k)
      rhogr=  1.0_wp / rhog
      zphf = p0(i,j,k) + pp(i,j,k,nu)

      zqrk = qrg * rhog

      zdh = hhl(i,j,k)-hhl(i,j,k+1) 
      zdtdh = 0.5_wp * zdt / zdh 

      zpkm1r(i,j) = zpkr(i,j)  
      zzar        = zqrk/zdtdh + zprvr(i,j) + zpkm1r(i,j) 
!      zpkr(i,j)   = zqrk * zvz0r * EXP(z1d8*LOG(MAX(zqrk,znull)))
      IF (zqrk.GT.znull) THEN
        zpkr(i,j)   = zqrk * zvz0r * EXP(z1d8*LOG(zqrk))
      ELSE
        zpkr(i,j)   = 0.0_wp
      ENDIF  
      zpkr(i,j)   = MIN( zpkr(i,j) , zzar )
      zzar        = zdtdh * (zzar-zpkr(i,j))

      IF( zvzr(i,j) == 0.0_wp ) zvzr(i,j) = zvz0r * EXP(z1d8*LOG(MAX(0.5_wp*zqrk,znull)))
 
      zimr        = 1.0_wp / (1.0_wp + zvzr(i,j) * zdtdh)
      zqrk        = zzar*zimr

!      lnzqrk      = LOG(MAX(zqrk,znull))  ! Optimization
      IF (zqrk.GT.znull) THEN
        lnzqrk      = LOG(zqrk)
      ELSE
        lnzqrk      = 0.0_wp
      ENDIF
        
      ztts        = 0.0_wp
      zqvts       = 0.0_wp
      zqcts       = 0.0_wp

      !------------------------------------------------------------------------
      !  Section 5: Calculation of cloud microphysics for cloud case
      !             ( qc > 0)
      !------------------------------------------------------------------------

      IF (qcg > 0.0_wp) THEN

        ! Calculate conversion rates

        ! Coefficients

        zc1c  = zaac *  EXP(lnzqrk*z7d8)
        zc1   = zaau + zc1c
        zx    = qcg / (1.0_wp + zc1*zdt)
   
        ! Conversion rates
        zswra  = zaau * zx
        zswrk  = zc1c * zx
 
        ! Store tendencies
        zqcts  = - zswra - zswrk
        zqrts  =   zswra + zswrk
 
        ! Update values

        qrg = MAX(0.0_wp,(zzar*rhogr + zqrts*zdt)*zimr)

      !------------------------------------------------------------------------
      !  Section 7: Calculation of cloud microphysics for
      !             precipitation case without cloud ( qc = 0 )
      !------------------------------------------------------------------------

      ELSEIF ( (zqrk) > 0.0_wp .AND. qcg <= 0.0_wp )    THEN

        zspw  = fspw(tg)
        zsqvw = fsqv(zspw, zphf)

        ! Coefficients
        zsrmax = zzar*rhogr * zdtr
        zsa3   = fsa3(tg)
        zc3  = zsa3 * SQRT(zqrk) * (1.0_wp + zbev * EXP(lnzqrk*z3d16))
 
        ! Conversion rates
        zsrd   = -zc3 * (qvg - zsqvw)
        zsrd   = MIN (zsrmax, zsrd)

 
        ! Store tendencies
        zqvts =   zsrd
        ztts  = - lh_v*zsrd*cpdr
        zqrts = - zsrd

        ! Update values

        qrg = MAX(0.0_wp,(zzar*rhogr + zqrts*zdt)*zimr)
 
      ENDIF         

      !------------------------------------------------------------------------
      ! Section 8: Complete time step
      !------------------------------------------------------------------------

      IF ( k /= ke ) THEN
        ! Store precipitation fluxes and sedimentation velocities for the next level
        zprvr(i,j) = qrg*rhog*zvzr(i,j)
        zvzr(i,j)  = zvz0r * EXP(z1d8 * LOG(MAX((qrg+qr(i,j,k+1))*0.5_wp*rhog,znull)))
      ELSE
        ! Precipitation flux at the ground
        prr_gsp(i,j) = 0.5_wp * (qrg*rhog*zvzr(i,j) + zpkr(i,j))
      ENDIF

      ! Update of prognostic variables or tendencies
      qr (i,j,k   ) = qrg
      qrs(i,j,k   ) = qrg
      t  (i,j,k,nu) = t (i,j,k,nu) + ztts*zdt 
      qv (i,j,k   ) = MAX ( 0.0_wp, qv(i,j,k) + zqvts*zdt )
      qc (i,j,k   ) = MAX ( 0.0_wp, qc(i,j,k) + zqcts*zdt )

    ENDDO
  ENDDO

  ! Do a final saturation adjustment for new values of t, qv and qc
  zpres(:,:) = p0(:,:,k) + pp(:,:,k,nu)

  CALL satad ( 1, t(:,:,k,nu), qv(:,:,k),                 &
               qc(:,:,k), t(:,:,k,nu), zpres,             &
               zdummy(:,:,1),zdummy(:,:,2),zdummy(:,:,3), &
               zdummy(:,:,4),zdummy(:,:,5),zdummy(:,:,6), &
               zdummy(:,:,7),zdummy(:,:,8),               &
               b1, b2w, b3, b4w, b234w, rdv, o_m_rdv,     &
               rvd_m_o, lh_v, cpdr, cp_d,                 &
               ie, je, istartpar, iendpar, jstartpar, jendpar)

  IF ( ldiabf_lh ) THEN
    ! compute temperature increment due to latent heat
    tinc_lh(:,:,k) = tinc_lh(:,:,k) + t(:,:,k,nu)
  ENDIF
    
ENDDO loop_over_levels

#ifdef NUDGING
! add part of latent heating calculated in subroutine kessler_pp to model latent
! heating field: add temperature to model latent heating field
IF (llhn) &
 CALL get_gs_lheating ('inc',1,ke)
#endif

!------------------------------------------------------------------------------
! End of module procedure kessler_pp
!------------------------------------------------------------------------------

END SUBROUTINE kessler_pp


!==============================================================================
!+ Module procedure "hydor_pp" in "gscp" for computing effects of grid scale 
!  precipitation including cloud water, rain and snow in 
!  context with the Leapfrog and the Runge-Kutta time-integration
!------------------------------------------------------------------------------

SUBROUTINE hydor_pp

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure calculates the rates of change of temperature, cloud 
!   water, water vapor, rain and snow due to cloud microphysical processes 
!   related to the formation of grid scale precipitation. The variables will 
!   be updated in this subroutine.
!   The precipitation fluxes at the 
!   surface are stored on the corresponding global fields.
!   This subroutine relies on conversion rates used in the subroutine hydorprog.
!
! Method:
!   The tendencies involving cloud water are calculated using a full implicit 
!   scheme whereas evaporation and deposition are computed explicitly.
!   Rain and snow are prognostic variables end the sedimentation term is 
!   computed implicitly.
!
!------------------------------------------------------------------------------
!
! Declarations:

! Subroutine arguments: None
! --------------------
!
! Local parameters:
! ----------------
  REAL    (KIND=wp   ),     PARAMETER ::  &
    ! basic constants of the parameterization scheme
    zamc   = 0.08_wp,                         & !
    zamv   = 0.02_wp,                         & !
    zt1    = 253.15_wp,                       & !
    zt2    = 235.15_wp,                       & !
    zaau   = 1.0_wp/10000.0_wp,           & ! coef. for autoconversion
    zaac   = 1.72_wp,                         & ! coef. for accretion (neu)
    zanuc  = 1.0_wp/1000.0_wp,            & ! coef. for nucleation
    zarim  = 1.97_wp,                         & ! coef. for riming (neu) 
    zamelt = 7.2E-6_wp,                       & ! coef. for melting (neu)
    zbev   = 9.1_wp,                          & ! coef. for drop ventilation
    zbdep  = 13.0_wp,                         & ! coef. for ice ventilation
    zbmelt = 13.0_wp,                         & ! coef. for melting ice

    ! constants for the process rates
    z03d4  = 3.0_wp/40.0_wp,              & !
    z09d4  = 9.0_wp/40.0_wp,              & !
    z43d4  = 43.0_wp/40.0_wp,             & !
    z7d4   = 7.0_wp/4.0_wp,               & !
    z1d8   = 1.0_wp/8.0_wp,               & !
    z5d8   = 5.0_wp/8.0_wp,               & !
    z7d8   = 7.0_wp/8.0_wp,               & !
    z13d8  = 13.0_wp/8.0_wp,              & !
    z3d16  = 3.0_wp/16.0_wp,              & !

    ! constants for rain freezing processes
    zfnpib = 1.0E4_wp,                        & !
    zfnpig = 2.3E5_wp,                        & !
    zfsaif = 9.95E-5_wp,                      & !
    zfsakf = 1.55E-3_wp*5.0E-3_wp,        & !

    ! constants for sedimentation
    zvz0r = 12.63_wp,                         & !
    zvz0s = 2.87_wp,                          & !

    ! to avoid illegal operations in power expressions
    znull = 1.E-20_wp

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    k     ,            & ! loop index in vertical direction              
    i     ,            & ! loop index in x-direction              
    j     ,            & ! loop index in y-direction              
    nu                   ! time level


  REAL    (KIND=wp   )     ::  &
    ztx   ,            & ! 
    zpx   ,            & !
    zgex  ,            & !  
    fspw  ,            & ! function for equilibrium vapour pressure over water
    fspi  ,            & ! function for equilibrium vapour pressure over ice
    fsqv  ,            & ! function for specific humidity at 
    fam   ,            & !
    feps  ,            & !
    fsa3  ,            & !
    fsa6  ,            & !
    zspi  ,            & ! equilibrium vapour pressure over ice
    zspw  ,            & ! equilibrium vapour pressure over water
    zam   ,            & !
    zeps  ,            & !
    zw2am ,            & !
    zw4am ,            & !
    zamcr2,            & !
    zamelt_am,         & !
    zarim_am,          & !
    zbmelt_am,         & !
    zbdep_am,          & !
    zsa3  ,            & !
    zsa6  ,            & !
    zc1a  ,            & !
    zc1b  ,            & !
    zc1c  ,            & !
    zc1d  ,            & !
    zc1   ,            & !
    zc3   ,            & !
    zc6   ,            & !
    zx    ,            & !
    zsrmax,            & !
    zsemax,            & !
    zzcool,            & !
    zznpi 

  REAL    (KIND=wp   )     ::  & 
    zphf  ,            & !  pressure in a k-layer
    zamr  ,            & !  reciprocal of form factor am
    zsqvw ,            & !  specific humidity at water saturation
    zsqvi ,            & !  specific humidity at ice saturation
    zqvts ,            & !  qv-tendency in a layer
    zqcts ,            & !  qc-tendency in a layer
    zqrts ,            & !  qr-tendency in a layer
    zqsts ,            & !  qs-tendency in a layer
    ztts  ,            & !  t -tendency in a layer
    zswra ,            & !  autoconversion rate
    zswrk ,            & !  accretion rate
    zswen ,            & !  nucleation rate
    zsweg ,            & !  riming rate
    zswegr,            & !  shedding rate
    zsde  ,            & !  deposition rate
    zser  ,            & !  melting rate
    zsrd  ,            & !  evaporation rate
    zsreif,            & !  immersion freezing rate
    zsrekf               !  contact freezing rate

  REAL    (KIND=wp   )     ::  &
    zvzr(ie,je),       & !
    zvzs(ie,je),       & !
    zprvr(ie,je),      & !
    zprvs(ie,je),      & !
    zpkr(ie,je),       & !
    zpks(ie,je),       & !
    zpkm1r(ie,je),     & !
    zpkm1s(ie,je),     & !
    zzar,              & !
    zzas,              & !
    zdh,               & !
    zdtdh,             & !
    zqrk,              & !
    lnzqrk,            & !
    zqsk,              & !
    lnzqsk,            & !
    zdt,               & ! 
    zdtr,              & ! 
    zimr,              & !
    zims,              & !
    qcg,               & !
    tg,                & !
    qvg,               & !
    qrg,               & !
    qsg,               & !
    rhog,              & !
    rhogr

!  LOGICAL                  ::  &
!    lcpexist             !

! Local (automatic) arrays:
! -------------------------
  REAL    (KIND=wp   )     ::  &
    zpres(ie,je)     ,& ! 
    zdummy(ie,je,8)      !

! Tracer pointers
!----------------
  REAL (KIND=wp),     POINTER :: &
    qv  (:,:,:) => NULL(),      & ! QV at nu
    qc  (:,:,:) => NULL(),      & ! QC at nu
    qr  (:,:,:) => NULL(),      & ! QR at nu
    qs  (:,:,:) => NULL()         ! QS at nu

  INTEGER (KIND=iintegers) :: izerror
  CHARACTER (LEN=255)      :: yzerrmsg

  CHARACTER (LEN=25)       :: yzroutine = 'hydor_pp'
!
!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine hydor_pp
!------------------------------------------------------------------------------

! Statement functions
! -------------------
 
fspw(ztx)     = b1*EXP( b2w*(ztx-b3)/(ztx-b4w) ) ! sat.vap. pressure over water
fspi(ztx)     = b1*EXP( b2i*(ztx-b3)/(ztx-b4i) ) ! sat.vap. pressure over ice
fsqv(zgex,zpx) = rdv*zgex/( zpx - o_m_rdv*zgex ) ! specific humidity at sat.
 
! interpreted ice/water ratio 
feps(ztx) = 0.5_wp * (1.0_wp+SIN(pi*(t0_melt-19.0_wp-ztx)/38.0_wp)) 
 
! form factor of ice  
fam(ztx) = zamc - zamv*(1.0_wp + COS(0.1_wp*pi*(ztx-t0_melt+10.0_wp))) 

! coefficients for conversion rates
fsa3(ztx) = 3.86E-3_wp - 9.41E-5_wp*(ztx-t0_melt) 
fsa6(ztx) = 1.09E-3_wp - 3.34E-5_wp*(ztx-t0_melt)
 

!------------------------------------------------------------------------------
!  Section 1: Initial setting of local variables
!
!------------------------------------------------------------------------------

! Precalculated coefficients for melting and shedding and deposition
zamcr2    = 1.0_wp/SQRT(zamc)
zamelt_am = zamelt*zamcr2    
zbmelt_am = zbmelt*SQRT(zamcr2)
zbdep_am  = zbdep *SQRT(zamcr2)
zarim_am  = zarim/zamc

zdummy(:,:,:)=0.0_wp

! Delete mixing ratio of rain and snow and precipitation fluxes from previous
! timestep
prr_gsp (:,:) = 0.0_wp
prs_gsp (:,:) = 0.0_wp

! select timelevel and timestep for calculations

nu    = nnew
IF ( l2tls ) THEN
  zdt   = dt
ELSE
  zdt   = dt2
ENDIF
zdtr  = 1.0_wp / zdt


zpkr(:,:) = 0.0_wp
zpks(:,:) = 0.0_wp
zprvr(:,:) = 0.0_wp
zprvs(:,:) = 0.0_wp
zvzr(:,:) = 0.0_wp
zvzs(:,:) = 0.0_wp

! retrieve the required microphysics tracers (at corresponding timelevel)
CALL trcr_get(izerror, idt_qv, ptr_tlev = nu, ptr = qv)
IF (izerror /= 0) THEN
  yzerrmsg = trcr_errorstr(izerror)
  CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
ENDIF
CALL trcr_get(izerror, idt_qc, ptr_tlev = nu, ptr = qc)
IF (izerror /= 0) THEN
  yzerrmsg = trcr_errorstr(izerror)
  CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
ENDIF
CALL trcr_get(izerror, idt_qr, ptr_tlev = nu, ptr = qr)
IF (izerror /= 0) THEN
  yzerrmsg = trcr_errorstr(izerror)
  CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
ENDIF
CALL trcr_get(izerror, idt_qs, ptr_tlev = nu, ptr = qs)
IF (izerror /= 0) THEN
  yzerrmsg = trcr_errorstr(izerror)
  CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
ENDIF

#ifdef NUDGING
  ! add part of latent heating calculated in subroutine hydor_pp to model latent
  ! heating field: subtract temperature from model latent heating field
    IF (llhn) &
       CALL get_gs_lheating ('add',1,ke)
#endif

loop_over_levels: DO k = 1, ke

  IF ( ldiabf_lh ) THEN
    ! initialize temperature increment due to latent heat
    tinc_lh(:,:,k) = tinc_lh(:,:,k) - t(:,:,k,nu)
  ENDIF

  !----------------------------------------------------------------------------
  !  Section 2: Test for clouds and precipitation in present layer.
  !             If no cloud or precipitation points have been found
  !             go to the next layer.
  !----------------------------------------------------------------------------

!  lcpexist = .FALSE.

!  loopj: DO j = jstartpar, jendpar
!    loopi: DO i = istartpar, iendpar
!      IF(qr (i,j,k,nu) < 1.0e-15_wp) qr (i,j,k,nu) = 0.0_wp
!      IF(qs (i,j,k,nu) < 1.0e-15_wp) qs (i,j,k,nu) = 0.0_wp
!     IF(qc (i,j,k,nu) < 1.0e-15_wp) qc (i,j,k,nu) = 0.0_wp
!     IF( (qc(i,j,k,nu) > 0.0_wp)   .OR. &
!          (qr(i,j,k,nu) > 0.0_wp)   .OR. &
!         (qs(i,j,k,nu) > 0.0_wp)   .OR. &
!          (zprvr(i,j)   > 0.0_wp)   .OR. &
!          (zprvs(i,j)   > 0.0_wp)   .OR. &
!          (zpkr(i,j)    > 0.0_wp)   .OR. &
!          (zpks(i,j)    > 0.0_wp) ) THEN
!          lcpexist = .TRUE.
!      ENDIF
!    ENDDO loopi
!    IF ( lcpexist ) EXIT loopj
!  ENDDO loopj

  !!! The following CYCLE has been eliminated, because the results are
  !!! not reproducible otherwise !!! (Should be investigated further on)

  !!! IF (.NOT. lcpexist) CYCLE loop_over_levels

  DO j = jstartpar, jendpar
    DO i = istartpar, iendpar

      qcg  = qc(i,j,k)
      qvg  = qv(i,j,k)
      qrg  = qr(i,j,k)
      qsg  = qs(i,j,k)
      tg   = t(i,j,k,nu)
      rhog = rho(i,j,k)
      rhogr=  1.0_wp / rhog
      zphf = p0(i,j,k) + pp(i,j,k,nu)

      zqrk = qrg * rhog
      zqsk = qsg * rhog

      zdh = hhl(i,j,k)-hhl(i,j,k+1) 
      zdtdh = 0.5_wp * zdt / zdh 

      zpkm1r(i,j) = zpkr(i,j) 
      zpkm1s(i,j) = zpks(i,j) 
      zzar        = zqrk/zdtdh + zprvr(i,j) + zpkm1r(i,j) 
      zzas        = zqsk/zdtdh + zprvs(i,j) + zpkm1s(i,j)
      zpkr(i,j)   = zqrk * zvz0r * EXP(z1d8*LOG(MAX(zqrk,znull)))
      zpks(i,j)   = zqsk * zvz0s * EXP(z03d4*LOG(MAX(zqsk,znull)))
      zpkr(i,j)   = MIN( zpkr(i,j) , zzar ) 
      zpks(i,j)   = MIN( zpks(i,j) , zzas ) 
      zzar        = zdtdh * (zzar-zpkr(i,j))
      zzas        = zdtdh * (zzas-zpks(i,j))

      IF( zvzr(i,j) == 0.0_wp ) zvzr(i,j) = zvz0r * EXP(z1d8*LOG(MAX(0.5_wp*zqrk,znull)))
      IF( zvzs(i,j) == 0.0_wp ) zvzs(i,j) = zvz0s * EXP(z03d4*LOG(MAX(0.5_wp*zqsk,znull)))
      zimr        = 1.0_wp / (1.0_wp + zvzr(i,j) * zdtdh)
      zims        = 1.0_wp / (1.0_wp + zvzs(i,j) * zdtdh)
      zqrk        = zzar*zimr
      zqsk        = zzas*zims

      lnzqrk      = LOG(MAX(zqrk,znull))  ! Optimization
      lnzqsk      = LOG(MAX(zqsk,znull))  ! Optimization

      ztts        = 0.0_wp
      zqvts       = 0.0_wp
      zqcts       = 0.0_wp
        
      !------------------------------------------------------------------------
      !  Section 4: Calculation of cloud microphysics for cold cloud case
      !             ( qc > 0, t < t0_melt )
      !------------------------------------------------------------------------

      IF (qcg > 0.0_wp .AND. tg <  t0_melt ) THEN

        zspi      = fspi(tg)
        zsqvi     = fsqv(zspi, zphf)
        zam  = zamc
        IF (tg > zt1)  zam  = fam (tg)
        zamr = 1.0_wp/zam
  
        ! Coefficients
        zeps  = 1.0_wp
        IF (tg > zt2)  zeps = feps(tg)
        zw2am = SQRT(zamr)
        zw4am = SQRT(zw2am)
        zsa6  = fsa6(tg)
 
        zc1a  = zanuc * zeps 
        zc1b  = zarim * zamr *  EXP(lnzqsk*z43d4)
        zc1c  = zaau * (1.0_wp - zeps)
        zc1d  = zaac * (1.0_wp - zeps) *  EXP(lnzqrk*z7d8)
        zc1   = zc1a + zc1b + zc1c + zc1d
        zc6   = zsa6*zw2am* EXP(lnzqsk*z5d8) * (1.0_wp+zbdep*zw4am*EXP(lnzqsk*z09d4))
        zx    = qcg / (1.0_wp + zc1*zdt)
 
        ! Conversion rates
        zswra   = zc1c * zx
        zswrk   = zc1d * zx
        zswen   = zc1a * zx
        zsweg   = zc1b * zx
        zsde    = zc6 * MAX( 0.0_wp, qvg-zsqvi )
        zsde    = MIN(zsde, qcg*zdtr)
        zsde    = zsde/(1.0_wp+ zc1*zdt)
        zsrmax  = zzar*rhogr*zdtr
        zzcool  = t0_melt - tg
        zsreif  = zfsaif*(EXP(0.66_wp*zzcool)-1.0_wp) * EXP(lnzqrk*z7d4)
        zsreif  = MIN(zsrmax, zsreif)
        zsrmax  = zsrmax - zsreif
        zznpi   = zfnpib + MAX(zfnpig*(zphf-5.E4_wp) &
                               /5.E4_wp, 0.0_wp)
        zzcool  = MAX(270.16_wp - tg, 0.0_wp)
        zsrekf  = zfsakf*zznpi*EXP(1.3_wp*LOG(zzcool+znull)) * EXP(lnzqrk*z13d8)
        zsrekf  = MIN(zsrmax, zsrekf)

        ! Store Tendencies 
        zqcts =     - zswra - zswrk - zswen - zsweg - zsde
        ztts  = lh_f * (  zswen + zsweg + zsde + zsreif + zsrekf )*cpdr
        zqrts =  zswra + zswrk - zsreif - zsrekf 
        zqsts =  zswen + zsweg + zsde + zsreif + zsrekf

        ! Update values 

        qrg = MAX(0.0_wp,(zzar*rhogr + zqrts*zdt)*zimr)
        qsg = MAX(0.0_wp,(zzas*rhogr + zqsts*zdt)*zims)

      !------------------------------------------------------------------------
      !  Section 5: Calculation of cloud microphysics for warm cloud case
      !             ( qc > 0, t >= t0_melt )
      !------------------------------------------------------------------------

      ELSEIF (qcg > 0.0_wp .AND. tg >= t0_melt) THEN

        ! Calculate conversion rates and precipitation fluxes               

        ! Coefficients
        zc1a  = zarim_am *  EXP(lnzqsk*z43d4)
        zc1b  = zaau
        zc1c  = zaac *  EXP(lnzqrk*z7d8)
        zc1   = zc1a + zc1b + zc1c
        zx    = qcg / (1.0_wp + zc1*zdt)
   
        ! Conversion rates
        zswra  = zc1b * zx
        zswrk  = zc1c * zx
        zswegr = zc1a * zx
        zser   = zamelt_am*EXP(lnzqsk*z5d8) * (1.0_wp+zbmelt_am*EXP(lnzqsk*z09d4)) * (tg -t0_melt)
        zsemax = zzas*rhogr *zdtr
        zser   = MIN(zser, zsemax)
 
        ! Store tendencies and integrate precipitation fluxes
        zqcts  = - zswra - zswrk - zswegr
        ztts   = - lh_f * zser * cpdr
        zqrts  =   zswra + zswrk + zswegr + zser
        zqsts  = - zser
  
        ! Update values

        qrg = MAX(0.0_wp,(zzar*rhogr + zqrts*zdt)*zimr)
        qsg = MAX(0.0_wp,(zzas*rhogr + zqsts*zdt)*zims)
 
      !------------------------------------------------------------------------
      !  Section 6: Calculation of cloud microphysics for
      !             cold precipitation case
      !             without cloud ( qc = 0, t <  t0_melt )
      !------------------------------------------------------------------------

      ELSEIF ( (zqrk+zqsk) > 0.0_wp .AND. qcg <= 0.0_wp &
                                        .AND. tg  <  t0_melt )    THEN

        zspw  = fspw(tg)
        zsqvw = fsqv(zspw, zphf)
        zspi  = fspi(tg)
        zsqvi = fsqv(zspi, zphf)
        zam       = zamc
        IF (tg > zt1) zam = fam (tg)
        zamr  = 1.0_wp/zam

        ! Coefficients
        zsrmax = zzar*rhogr * zdtr
        zsemax = zzas*rhogr * zdtr
        zw2am  = SQRT(zamr)
        zw4am  = SQRT(zw2am)
        zsa6   = fsa6(tg)
        zsa3   = fsa3(tg)
        zc6    = zsa6 * zw2am *  EXP(lnzqsk*z5d8)*(1.0_wp + zbdep * zw4am * EXP(lnzqsk*z09d4))
        zc3    = zsa3 * SQRT(zqrk) * (1.0_wp + zbev * EXP(lnzqrk*z3d16))
 
        ! Conversion rates
        zsde    =  zc6 * (qvg - zsqvi)
        zsde    = MAX (-zsemax, zsde)
        zsrd    = -zc3 * (qvg - zsqvw)
        zsrd    = MIN ( zsrmax, zsrd)
        zsrmax  = zsrmax - zsrd
        zzcool  = t0_melt - tg
        zsreif  = zfsaif*(EXP(0.66_wp*zzcool)-1.0_wp) * EXP(lnzqrk*z7d4)
        zsreif  = MIN(zsrmax, zsreif)
        zsrmax  = zsrmax - zsreif
        zznpi   = zfnpib + MAX(zfnpig*(zphf-5.E4_wp)/5.E4_wp, 0.0_wp)
        zzcool  = MAX(270.16_wp - tg, 0.0_wp)
        zsrekf  = zfsakf*zznpi*EXP(1.3_wp*LOG(MAX(zzcool,znull))) * EXP(lnzqrk*z13d8)
        zsrekf  = MIN(zsrmax, zsrekf)

        ! Store tendencies and integrate precipitation fluxes
        zqvts = - zsde + zsrd
        ztts  =  ( lh_f*(zsde + zsreif + zsrekf) + lh_v*(zsde - zsrd) ) * cpdr
        zqrts = - (zsrd + zsreif + zsrekf)
        zqsts =    zsde + zsreif + zsrekf
 
        ! Update values

        qrg = MAX(0.0_wp,(zzar*rhogr + zqrts*zdt)*zimr)
        qsg = MAX(0.0_wp,(zzas*rhogr + zqsts*zdt)*zims)


      !------------------------------------------------------------------------
      !  Section 7: Calculation of cloud microphysics for
      !             warm precipitation case
      !             without cloud ( qc = 0, t >= t0_melt )
      !------------------------------------------------------------------------

      ELSEIF ( (zqrk+zqsk) > 0.0_wp .AND. qcg <= 0.0_wp &
                                        .AND. tg  >= t0_melt )    THEN

        zspw  = fspw(tg)
        zsqvw = fsqv(zspw, zphf)

        ! Coefficients
        zsrmax = zzar*rhogr * zdtr
        zsemax = zzas*rhogr * zdtr
        zsa6   = fsa6(tg)
        zsa3   = fsa3(tg)
        zc6  = zsa6 * zamcr2 * EXP(lnzqsk*z5d8) * (1.0_wp + zbdep_am *  EXP(lnzqsk*z09d4))
        zc3  = zsa3 * SQRT(zqrk) * (1.0_wp + zbev * EXP(lnzqrk*z3d16))
 
        ! Conversion rates
        zsde   =  zc6 * (qvg - zsqvw)
        zsde   = MAX (-zsemax, zsde)
        zsemax = zsemax + zsde
        zsrd   = -zc3 * (qvg - zsqvw)
        zsrd   = MIN (zsrmax, zsrd)

        zser   = 0.0_wp
        IF (zqsk > 0.0_wp)  THEN
          zser = zamelt_am * EXP(lnzqsk*z5d8) * (1.0_wp + zbmelt_am * EXP(lnzqsk*z09d4))*(tg-t0_melt)
          zser = MIN (zser, zsemax)
        ENDIF
 
        ! Store tendencies and integrate precipitation fluxes
        zqvts =     -  zsde + zsrd
        ztts  = ( lh_f*(zsde-zser) + lh_v*(zsde-zsrd) )*cpdr
        zqrts = zser - zsrd
        zqsts = zsde - zser

        ! Update values

        qrg = MAX(0.0_wp,(zzar*rhogr + zqrts*zdt)*zimr)
        qsg = MAX(0.0_wp,(zzas*rhogr + zqsts*zdt)*zims)
 
      ENDIF         

      !------------------------------------------------------------------------
      ! Section 8: Complete time step
      !------------------------------------------------------------------------

      IF ( k /= ke ) THEN
        ! Store precipitation fluxes and sedimentation velocities for the next level
        zprvr(i,j) = qrg*rhog*zvzr(i,j)
        zprvs(i,j) = qsg*rhog*zvzs(i,j)
        zvzr(i,j)  = zvz0r * EXP(z1d8 * LOG(MAX((qrg+qr(i,j,k+1))*0.5_wp*rhog,znull)))
        zvzs(i,j)  = zvz0s * EXP(z03d4 * LOG(MAX((qsg+qs(i,j,k+1))*0.5_wp*rhog,znull)))
      ELSE
        ! Precipitation fluxes at the ground
        prr_gsp(i,j) = 0.5_wp * (qrg*rhog*zvzr(i,j) + zpkr(i,j))
        prs_gsp(i,j) = 0.5_wp * (qsg*rhog*zvzs(i,j) + zpks(i,j))
      ENDIF

      ! Update of prognostic variables
      qr (i,j,k   ) = qrg
      qs (i,j,k   ) = qsg
      qrs(i,j,k   ) = qrg+qsg
      t  (i,j,k,nu) = t (i,j,k,nu) + ztts*zdt 
      qv (i,j,k   ) = MAX ( 0.0_wp, qv(i,j,k) + zqvts*zdt )
      qc (i,j,k   ) = MAX ( 0.0_wp, qc(i,j,k) + zqcts*zdt )
 
    ENDDO
  ENDDO

  ! Do a final saturation adjustment for new values of t, qv and qc
  zpres(:,:) = p0(:,:,k) + pp(:,:,k,nu)

  CALL satad ( 1, t(:,:,k,nu), qv(:,:,k),                 &
               qc(:,:,k), t(:,:,k,nu), zpres,             &
               zdummy(:,:,1),zdummy(:,:,2),zdummy(:,:,3), &
               zdummy(:,:,4),zdummy(:,:,5),zdummy(:,:,6), &
               zdummy(:,:,7),zdummy(:,:,8),               &
               b1, b2w, b3, b4w, b234w, rdv, o_m_rdv,     &
               rvd_m_o, lh_v, cpdr, cp_d,                 &
               ie, je, istartpar, iendpar, jstartpar, jendpar)

  IF ( ldiabf_lh ) THEN
    ! compute temperature increment due to latent heat
    tinc_lh(:,:,k) = tinc_lh(:,:,k) + t(:,:,k,nu)
  ENDIF
    
ENDDO loop_over_levels

#ifdef NUDGING
! add part of latent heating calculated in subroutine hydor_pp to model latent
! heating field: add temperature to model latent heating field
IF (llhn) &
 CALL get_gs_lheating ('inc',1,ke)
#endif

!------------------------------------------------------------------------------
! End of module procedure hydor_pp
!------------------------------------------------------------------------------

END SUBROUTINE hydor_pp

!==============================================================================
!+ Module procedure "hydci_pp" in "gscp" for computing effects of grid scale
!  precipitation including cloud water, cloud ice, rain and snow in
!  context with the Leapfrog and the Runge-Kutta time-integration
!------------------------------------------------------------------------------

SUBROUTINE hydci_pp

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure calculates the rates of change of temperature, cloud
!   water, cloud ice, water vapor, rain and snow due to cloud microphysical
!   processe related to the formation of grid scale precipitation. The
!   variables are updated in this subroutine. Rain and snow are prognostic
!   variables. The precipitation fluxes at the surface are stored on the
!   corresponding global fields.
!   The subroutine relies on conversion terms used in hydci.
!
! Method:
!   The sedimentation of rain and snow is computed implicitly.
!
!------------------------------------------------------------------------------
! Declarations:
!
!------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

! Subroutine arguments: None
! --------------------
!
! Local parameters:
! ----------------

  INTEGER (KIND=iintegers), PARAMETER ::  &
    iautocon       = 1,&
    isnow_n0temp   = 2

  REAL    (KIND=wp   ),     PARAMETER ::  &
    zqmin = 1.0E-15_wp,& ! threshold for computations

    ! Parameters for autoconversion of cloud water and cloud ice 
    zccau  = 4.0E-4_wp, & ! autoconversion coefficient (cloud water to rain)
    zciau  = 1.0E-3_wp, & ! autoconversion coefficient (cloud ice   to snow)

    zkcau  = 9.44e+09_wp, & ! kernel coeff for SB2001 autoconversion
    zkcac  = 5.25e+00_wp, & ! kernel coeff for SB2001 accretion
    zcnue  = 2.00e+00_wp, & ! gamma exponent for cloud distribution
    zxstar = 2.60e-10_wp, & ! separating mass between cloud and rain
    zkphi1 = 6.00e+02_wp, & ! constant in phi-function for autoconversion
    zkphi2 = 0.68e+00_wp, & ! exponent in phi-function for autoconversion
    zkphi3 = 5.00e-05_wp, & ! exponent in phi-function for accretion

    ! Other process coefficients.
    ! These coefficients have been precalculated on the basis of the
    ! following basic parameters (original hydci):
    ! N0R = 8.0 E6 : Parameter in the size distrubution function for rain
    ! ECR = 0.8    : Collection efficiency for rain collecting cloud water
    ! ECS = 0.9    : Collection efficiency for snow collecting cloud water
    ! EIR = 0.8    : Collection efficiency for rain collecting cloud ice
    ! V0R = 130.0  : Factor in the terminal velocity for raindrops
    !                VTR(D) = V0R*D**(1/2)
    ! AMS = 0.038  : Formfactor in the mass-size relation of snowparticles
    !                m(DS) = AMS*DS**2
    ! AMI = 130.0  : Formfactor in the mass-size relation of cloud ice crystals
    !                m(DI) = AMI*DI**3
    ! ETA  = 1.75 E-5 : Viscosity of air    
    ! DIFF = 2.22 E-5 : Molecular diffusion coefficient for water vapour
    ! LHEAT= 2.40 E-2 : Heat conductivity of air
    ! RHO  = 1.0      : Density of air 
    ! RHOW = 1.0 E3   : Density of water 
    ! lh_v            : latent heat of vapourization
    ! lh_f            : latent heat of fusion
    ! AR  = PI*RHOW*N0R
    ! AS  = 2*N0S*AMS
    ! HW              : Howell factor ( =(1/(1+H_w)) or =(1/(1+H_i))) 

    zhw   = 2.270603_wp,     & ! Howell factor
    zecs  = 0.9_wp,   & ! Collection efficiency for snow collecting cloud water

    zadi  = 0.217_wp, & ! Formfactor in the size-mass relation of ice particles
    zbdi  = 0.302_wp, & ! Exponent in the size-mass relation of ice particles
    zams  = 0.069_wp, & ! Formfactor in the mass-size relation of snow particles
    zbms  = 2.000_wp, & ! Exponent in the mass-size relation of snow particles

    zv1s  = 0.50_wp,  & ! Exponent in the terminal velocity for snow

    zami  = 130.0_wp, & ! Formfactor in the mass-size relation of cloud ice
    zn0s0 = 8.0E5_wp, & ! 
    zn0s1 = 13.5_wp * 5.65E5_wp, & ! parameter in N0S(T)
    zn0s2 = -0.107_wp , & ! parameter in N0S(T), Field et al
    zcac  = 1.72_wp   , & ! (15/32)*(PI**0.5)*(ECR/RHOW)*V0R*AR**(1/8)
    zcicri= 1.72_wp   , & ! (15/32)*(PI**0.5)*(EIR/RHOW)*V0R*AR**(1/8)
    zcrcri= 1.24E-3_wp, & ! (PI/24)*EIR*V0R*Gamma(6.5)*AR**(-5/8)
!    zcev  = 3.1E-3_wp , & ! 2*PI*DIFF*HW*N0R*AR**(-1/2)
!    zbev  = 9.0_wp    , & ! 0.26*sqrt(0.5*RHO*v0r/eta)*Gamma(2.75)*AR**(-3/16)
    zcsmel= 1.48E-4_wp, & ! 4*LHEAT*N0S*AS**(-2/3)/(RHO*lh_f)

!FR old    
!   zbsmel= 14.37_wp  , & ! 0.26_wp*sqrt(0.5_wp*RHO*v0s/eta)*Gamma(21.0_wp/8.0_wp)*AS**(-5.0_wp/24.0_wp)
    zbsmel= 20.32_wp  , & ! 0.26_wp*sqrt(       RHO*v0s/eta)*Gamma(21.0_wp/8.0_wp)*AS**(-5.0_wp/24.0_wp)

!FR old    
!   zasmel= 2.31E3_wp , & ! DIFF*lh_v*RHO/LHEAT
    zasmel= 2.43E3_wp , & ! DIFF*lh_v*RHO/LHEAT

    zcrfrz= 1.68_wp   , & ! coefficient for raindrop freezing

    zrho0 = 1.225e+0_wp, & ! reference air density
    zrhow = 1.000e+3_wp, & ! density of liquid water

    zdv    = 2.22e-5_wp, & ! molecular diffusion coefficient for water vapour
    zlheat = 2.40E-2_wp, & ! thermal conductivity of dry air
    zeta   = 1.75e-5_wp    ! kinematic viscosity of air 

  REAL    (KIND=wp   ),     PARAMETER ::  &
    ! Additional parameters
    zthet = 248.15_wp , & ! temperature for het. nuc. of cloud ice
    zthn  = 236.15_wp , & ! temperature for hom. freezing of cloud water
    ztrfrz= 271.15_wp , & ! threshold temperature for heterogeneous
                              ! freezing of raindrops
    zmi0   = 1.0E-12_wp, & ! initial crystal mass for cloud ice nucleation
    zmimax = 1.0E-9_wp , & ! maximum mass of cloud ice crystals   
    zmsmin = 3.0E-9_wp , & ! initial mass of snow crystals        

    ! Constant exponents in the transfer rate equations
    x1o12  =  1.0_wp/12.0_wp  ,  x3o16  =  3.0_wp/16.0_wp, &
    x7o8   =  7.0_wp/ 8.0_wp  ,  x2o3   =  2.0_wp/ 3.0_wp, &
    x5o24  =  5.0_wp/24.0_wp  ,  x1o8   =  1.0_wp/ 8.0_wp, &
    x13o8  = 13.0_wp/ 8.0_wp  ,  x13o12 = 13.0_wp/12.0_wp, &
    x27o16 = 27.0_wp/16.0_wp  ,  x1o3   =  1.0_wp/ 3.0_wp, &
    x1o2   =  1.0_wp/ 2.0_wp

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    k, i, j,           & ! loop indees
    nx, nn, izdebug      ! time-level for computation

  REAL    (KIND=wp   )     ::  &
    fpvsw, fpvsi, fqvs,& ! name of statement functions
    fxna ,             & ! statement function for ice crystal number
    ztx  , zpv  , zpx ,& ! dummy arguments for statement functions
    znimax,            & ! maximum number of cloud ice crystals
    zpvsw0,            & ! sat.vap. pressure at melting temperature
    zqvsw0,            & ! sat.specific humidity at melting temperature
    zdt,               & ! timestep for integration (water / ice )
    zdtr ,             & ! reciprocal of timestep for integration
    zscau , zscac  , zscrim , zscshe, zsnuc , & ! local values of the
    zsiau , zsagg  , zsidep , zsicri, zsrcri, & ! transfer rates
    zsdau , zssdep , zssmelt,                 & ! defined below
    zsrfrz,                                   & !
    zscsum, zscmax, zcorr,  & ! terms for limiting  total cloud water depletion
    zsrsum, zsrmax,    & ! terms for limiting  total rain water depletion
    zssmax,            & ! term for limiting snow depletion
    znin,              & ! number of cloud ice crystals at nucleation
    znid,              & ! number of cloud ice crystals for deposition
    zmi ,              & ! mass of a cloud ice crystal
    zsvidep, zsvisub,  & ! deposition, sublimation of cloud ice
    zsimax , zsisum , zsvmax,   & ! terms for limiting total
    zqvsi, zqvsw,      & ! sat. specitic humidity at ice and water saturation
    zztau, zxfac, zx1, zx2, ztt, &   ! some help variables
    ztau, zphi, zhi, zdvtp, ztc

  REAL    (KIND=wp   )     ::  &
    zqct   ,& ! layer tendency of cloud water
    zqvt   ,& ! layer tendency of water vapour
    zqit   ,& ! layer tendency of cloud ice
    zqrt   ,& ! layer tendency of rain
    zqst      ! layer tendency of snow

  REAL (KIND=wp)             ::       &
    mma(10), mmb(10)
    
  REAL    (KIND=wp   )     ::  &
    zdh, zlnqrk,zlnqsk,     & !
    zlnlogmi, zar, zn0r,    & !
    qcg,tg,qvg,qrg, qsg,qig,rhog,ppg, alf,bet,m2s,m3s,hlp

  LOGICAL :: &
    llqr,llqs,llqc,llqi  !   switch for existence of qr, qs, qc, qi

! Local (automatic) arrays:
! -------------------------
  REAL    (KIND=wp   )     ::  &
    zpres       (ie,je),     & ! pressure
    zvzr        (ie,je),     & !
    zvzs        (ie,je),     & !
    zpkr        (ie,je),     & !
    zpks        (ie,je),     & !
    zprvr       (ie,je),     & !
    zprvs       (ie,je),     & !
    zdummy      (ie,je,8),   & !
    zcsdep      (ie,je),     & !
    zcidep      (ie,je),     & !
    zvz0s       (ie,je),     & !
    zcrim       (ie,je),     & !
    zcagg       (ie,je),     & !
    zbsdep      (ie,je),     & !
    zcslam      (ie,je),     & !
    zn0s        (ie,je),     & !
    zimr        (ie,je),     & !
    zims        (ie,je),     & !
    zzar        (ie,je),     & !
    zzas        (ie,je),     & !
    zqrk        (ie,je),     & !
    zqsk        (ie,je),     & !
    zdtdh       (ie,je),     & !
    z1orhog     (ie,je),     & ! 1/rhog
    zrho1o2     (ie,je),     & ! (rho0/rhog)**1/2
    zeln7o8qrk  (ie,je),     & !
    zeln27o16qrk(ie,je),     & !
    zeln13o8qrk (ie,je),     & !
!!    zeln3o16qrk (ie,je),     & !
    zeln13o12qsk(ie,je),     & !
    zeln5o24qsk (ie,je),     & !
    zeln2o3qsk  (ie,je)        !

  REAL    (KIND=wp   )     ::  &
    scau   (ie,je), & ! transfer rate due to autoconversion of cloud water
    scac   (ie,je), & ! transfer rate due to accretion of cloud water
    snuc   (ie,je), & ! transfer rate due nucleation of cloud ice
    scfrz  (ie,je), & ! transfer rate due homogeneous freezing of cloud water
    simelt (ie,je), & ! transfer rate due melting of cloud ice
    sidep  (ie,je), & ! transfer rate due depositional growth of cloud ice
    ssdep  (ie,je), & ! transfer rate due depositional growth of snow
    sdau   (ie,je), & ! transfer rate due depositional cloud ice autoconversion
    srim   (ie,je), & ! transfer rate due riming of snow
    sshed  (ie,je), & ! transfer rate due shedding
    sicri  (ie,je), & ! transfer rate due cloud ice collection by rain (sink qi)
    srcri  (ie,je), & ! transfer rate due cloud ice collection by rain (sink qr)
    sagg   (ie,je), & ! transfer rate due aggregation of snow and cloud ice
    siau   (ie,je), & ! transfer rate due autoconversion of cloud ice
    ssmelt (ie,je), & ! transfer rate due melting of snow
    sev    (ie,je), & ! transfer rate due evaporation of rain
    srfrz  (ie,je)    ! transfer rate due to rainwater freezing

  ! Integer arrays for a better vectorization
  INTEGER (KIND=iintegers) ::  &
    idx1(ie*je),       jdx1(ie*je),   & !
    idx2(ie*je),       jdx2(ie*je),   & !
    idx3(ie*je),       jdx3(ie*je),   & !
    idx4(ie*je),       jdx4(ie*je),   & !
    idx5(ie*je),       jdx5(ie*je),   & !
    idx6(ie*je),       jdx6(ie*je)      !

  ! Dimensions and loop counter for storing the indices
  INTEGER (KIND=iintegers) ::  &
    ic1, ic2, ic3, ic4, ic5, ic6, i1d

! Tracer pointers
!----------------
  REAL (KIND=wp),     POINTER :: &
    qv  (:,:,:) => NULL(),      & ! QV at nx
    qc  (:,:,:) => NULL(),      & ! QC at nx
    qi  (:,:,:) => NULL(),      & ! QI at nx
    qr  (:,:,:) => NULL(),      & ! QR at nx
    qs  (:,:,:) => NULL()         ! QS at nx

  INTEGER(KIND=iintegers)     :: izerror
  CHARACTER (LEN=255)         :: yzerrmsg

  CHARACTER(LEN=25)           :: yzroutine = 'hydci_pp'

  LOGICAL :: ldum

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine hydci_pp
!------------------------------------------------------------------------------

! Statement functions
! -------------------

! saturation vapour pressure over water (fpvsw), over ice (fpvsi)
! and specific humidity at vapour saturation (fqvs)
  fpvsw(ztx)     = b1*EXP( b2w*(ztx-b3)/(ztx-b4w) )
  fpvsi(ztx)     = b1*EXP( b2i*(ztx-b3)/(ztx-b4i) )
  fqvs (zpv,zpx) = rdv*zpv/( zpx - o_m_rdv*zpv )

! Number of activate ice crystals;  ztx is temperature
  fxna(ztx)   = 1.0E2_wp * EXP(0.2_wp * (t0_melt - ztx))

! Coeffs for moment relation based on 2nd moment (Field 2005)
  mma = (/   5.065339_wp, -0.062659_wp, -3.032362_wp, 0.029469_wp, -0.000285_wp, &
             0.312550_wp,  0.000204_wp,  0.003199_wp, 0.000000_wp, -0.015952_wp /)
  mmb = (/   0.476221_wp, -0.015896_wp,  0.165977_wp, 0.007468_wp, -0.000141_wp, &
             0.060366_wp,  0.000079_wp,  0.000594_wp, 0.000000_wp, -0.003577_wp /)

IF (ldebug_gsp) THEN
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

!------------------------------------------------------------------------------
!  Section 1: Initial setting of local and global variables
!------------------------------------------------------------------------------

!!  mu_rain = 0.5  ! is a namelist parameter

  IF (ntstep == nstart) THEN
    zconst = zkcau / (20.0_wp*zxstar*cloud_num*cloud_num) &
              * (zcnue+2.0_wp)*(zcnue+4.0_wp)/(zcnue+1.0_wp)**2
    ccsrim = 0.25_wp*pi*zecs*v0snow*gamma_fct(zv1s+3.0_wp)
    ccsagg = 0.25_wp*pi*v0snow*gamma_fct(zv1s+3.0_wp)

!FR old    
!   ccsdep = 0.26_wp*gamma_fct((zv1s+5.0_wp)/2.0_wp)*SQRT(0.5_wp/zeta)
    ccsdep = 0.26_wp*gamma_fct((zv1s+5.0_wp)/2.0_wp)*SQRT(1.0_wp/zeta)

    ccsvxp = -(zv1s/(zbms+1.0_wp)+1.0_wp)
    ccsvel = zams*v0snow*gamma_fct(zbms+zv1s+1.0_wp)&
      & *(zams*gamma_fct(zbms+1.0_wp))**ccsvxp 
    ccsvxp = ccsvxp + 1.0_wp
    ccslam = zams*gamma_fct(zbms+1.0_wp)
    ccslxp = 1.0_wp / (zbms+1.0_wp)
    ccswxp = zv1s*ccslxp
    ccsaxp = -(zv1s+3.0_wp)

!FR old    
!   ccsdxp = -(zbms+1.0_wp)/2.0_wp
    ccsdxp = -(zv1s+1.0_wp)/2.0_wp

    ccshi1 = lh_s*lh_s/(zlheat*r_v)

!FR old    
!   ccdvtp = 2.11E-5_wp * t0_melt**(-1.94_wp) * 101325.0_wp
    ccdvtp = 2.22E-5_wp * t0_melt**(-1.94_wp) * 101325.0_wp

    ccidep = 4.0_wp * zami**(-x1o3)
    ! empirical relation adapted from Ulbrich (1983)
    zn0r   = 8.0E6_wp * EXP(3.2_wp*mu_rain) * (0.01_wp)**(-mu_rain)  
    ! to tune the zn0r variable
    zn0r   = zn0r * rain_n0_factor

    zar    = pi*zrhow/6.0_wp * zn0r * gamma_fct(mu_rain+4.0_wp) ! pre-factor in lambda
    zcevxp = (mu_rain+2.0_wp)/(mu_rain+4.0_wp)
    zcev   = 2.0_wp*pi*zdv/zhw*zn0r*zar**(-zcevxp) * gamma_fct(mu_rain+2.0_wp)
    zbevxp = (2.0_wp*mu_rain+5.5_wp)/(2.0_wp*mu_rain+8.0_wp)-zcevxp

!FR old 
!   zbev   = 0.26_wp * SQRT(0.5_wp*zrho0*130.0_wp/zeta)*zar**(-zbevxp)
    zbev   = 0.26_wp * SQRT(       zrho0*130.0_wp/zeta)*zar**(-zbevxp) &
           * gamma_fct((2.0_wp*mu_rain+5.5_wp)/2.0_wp) / gamma_fct(mu_rain+2.0_wp)

    zvzxp  = 0.5_wp/(mu_rain+4.0_wp)
    zvz0r  = 130.0_wp*gamma_fct(mu_rain+4.5_wp)/gamma_fct(mu_rain+4.0_wp)*zar**(-zvzxp)
    IF (izdebug > 10) THEN
      WRITE (*,*) 'SRC_GSCP: Initialized hydci_pp'
      WRITE (*,'(A,E10.3)') '      ccslam = ',ccslam
      WRITE (*,'(A,E10.3)') '      ccsvel = ',ccsvel
      WRITE (*,'(A,E10.3)') '      ccsrim = ',ccsrim
      WRITE (*,'(A,E10.3)') '      ccsagg = ',ccsagg
      WRITE (*,'(A,E10.3)') '      ccsdep = ',ccsdep
      WRITE (*,'(A,E10.3)') '      ccslxp = ',ccslxp
      WRITE (*,'(A,E10.3)') '      ccidep = ',ccidep
      WRITE (*,'(A,E10.3)') '      mu_r   = ',mu_rain
      WRITE (*,'(A,E10.3)') '      zn0r   = ',zn0r
      WRITE (*,'(A,E10.3)') '      zbev   = ',zbev
      WRITE (*,'(A,E10.3)') '      zbevxp = ',zbevxp
      WRITE (*,'(A,E10.3)') '      zcev   = ',zcev
      WRITE (*,'(A,E10.3)') '      zcevxp = ',zcevxp
      WRITE (*,'(A,E10.3)') '      zvzxp  = ',zvzxp
      WRITE (*,'(A,E10.3)') '      zvz0r  = ',zvz0r
    ENDIF
  ENDIF

! Some constant coefficients
  znimax = fxna(zthn) ! Maximum number of cloud ice crystals
  zpvsw0 = fpvsw(t0_melt)  ! sat. vap. pressure for t = t0_melt

! Delete precipitation fluxes from previous timestep
!CDIR BEGIN COLLAPSE
  prr_gsp (:,:) = 0.0_wp
  prs_gsp (:,:) = 0.0_wp
  zpkr    (:,:) = 0.0_wp
  zpks    (:,:) = 0.0_wp
  zprvr   (:,:) = 0.0_wp
  zprvs   (:,:) = 0.0_wp
  zvzr    (:,:) = 0.0_wp
  zvzs    (:,:) = 0.0_wp
! zdummy(:,:,:) = 0.0_wp
!CDIR END
  
! select timelevel and timestep for calculations
  nx    = nnew
  IF ( l2tls ) THEN
    zdt   = dt
  ELSE
    zdt   = dt2
  ENDIF
  zdtr  = 1.0_wp / zdt

  ! retrieve the required microphysics tracers (at corresponding timelevel)
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nx, ptr = qv)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nx, ptr = qc)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qi, ptr_tlev = nx, ptr = qi)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qr, ptr_tlev = nx, ptr = qr)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qs, ptr_tlev = nx, ptr = qs)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

#ifdef NUDGING
  ! add part of latent heating calculated in subroutine hydci_pp to model latent
  ! heating field: subtract temperature from model latent heating field
  IF (llhn .OR. llhnverif) THEN
    IF (lhn_qrs) THEN
!CDIR COLLAPSE
      qrsflux(:,:,:) = 0.0_wp
    ENDIF
    CALL get_gs_lheating ('add',1,ke)
  ENDIF
#endif

! *********************************************************************
! Loop from the top of the model domain to the surface to calculate the
! transfer rates  and sedimentation terms    
! *********************************************************************

  loop_over_levels: DO  k = 1, ke

    IF ( ldiabf_lh ) THEN
      ! initialize temperature increment due to latent heat
      tinc_lh(:,:,k) = tinc_lh(:,:,k) - t(:,:,k,nx)
    ENDIF

  !----------------------------------------------------------------------------
  ! Section 2: Check for existence of rain and snow
  !            Initialize microphysics and sedimentation scheme
  !----------------------------------------------------------------------------

!CDIR BEGIN COLLAPSE
    zcrim (:,:) = 0.0_wp
    zcagg (:,:) = 0.0_wp
    zbsdep(:,:) = 0.0_wp
    zvz0s (:,:) = 0.0_wp
    zn0s  (:,:) = zn0s0
!CDIR END

    ic1 = 0
    ic2 = 0

    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar

        qrg = qr(i,j,k)
        qsg = qs(i,j,k)
        qvg = qv(i,j,k)
        qcg = qc(i,j,k)
        qig = qi(i,j,k)
        tg  = t(i,j,k,nx)
        ppg = p0(i,j,k) + pp(i,j,k,nx)
        rhog = rho(i,j,k)

        !..for density correction of fall speeds
        z1orhog(i,j) = 1.0_wp/rhog
        zrho1o2(i,j) = EXP(LOG(zrho0*z1orhog(i,j))*x1o2)

        zqrk(i,j) = qrg * rhog
        zqsk(i,j) = qsg * rhog

        llqr = zqrk(i,j) > zqmin
        llqs = zqsk(i,j) > zqmin

        zdh = hhl(i,j,k)-hhl(i,j,k+1)
        zdtdh(i,j) = 0.5_wp * zdt / zdh

        zzar(i,j)   = zqrk(i,j)/zdtdh(i,j) + zprvr(i,j) + zpkr(i,j)
        zzas(i,j)   = zqsk(i,j)/zdtdh(i,j) + zprvs(i,j) + zpks(i,j)

        IF (llqs) THEN
          ic1 = ic1 + 1
          idx1(ic1) = i
          jdx1(ic1) = j
        ENDIF
        IF (llqr) THEN
          ic2 = ic2 + 1
          idx2(ic2) = i
          jdx2(ic2) = j
        ENDIF
      ENDDO
    ENDDO

!CDIR NODEP,VOVERTAKE,VOB
    DO i1d = 1, ic1
      i = idx1(i1d)
      j = jdx1(i1d)

      qsg = qs(i,j,k)
      tg  = t(i,j,k,nx)

      IF (isnow_n0temp == 1) THEN
        ! Calculate n0s using the temperature-dependent
        ! formula of Field et al. (2005)
        ztc = tg - t0_melt
        ztc = MAX(MIN(ztc,0.0_wp),-40.0_wp)          
        zn0s(i,j) = zn0s1*EXP(zn0s2*ztc)
        zn0s(i,j) = MIN(zn0s(i,j),1e9_wp)
        zn0s(i,j) = MAX(zn0s(i,j),1e6_wp)
      ELSEIF (isnow_n0temp == 2) THEN
        ! Calculate n0s using the temperature-dependent moment 
        ! relations of Field et al. (2005)
        ztc = tg - t0_melt
        ztc = MAX(MIN(ztc,0.0_wp),-40.0_wp)
        
        nn  = 3
        hlp = mma(1)      +mma(2)*ztc      +mma(3)*nn       +mma(4)*ztc*nn+mma(5)*ztc**2 &
            + mma(6)*nn**2+mma(7)*ztc**2*nn+mma(8)*ztc*nn**2+mma(9)*ztc**3+mma(10)*nn**3
        alf = 10.0_wp**hlp
        bet = mmb(1)      +mmb(2)*ztc      +mmb(3)*nn       +mmb(4)*ztc*nn+mmb(5)*ztc**2 &
            + mmb(6)*nn**2+mmb(7)*ztc**2*nn+mmb(8)*ztc*nn**2+mmb(9)*ztc**3+mmb(10)*nn**3

        ! Here is the exponent bms=2.0 hardwired! not ideal! (Uli Blahak)
        m2s = qsg * rho(i,j,k) / zams   ! UB rho added as bugfix
        m3s = alf*EXP(bet*LOG(m2s))

        hlp  = zn0s1*EXP(zn0s2*ztc)
        ! UB: the 13.5 is actually 3^3 / gamma(3) ...
        zn0s(i,j) = 13.50_wp * m2s*(m2s/m3s)**3
        zn0s(i,j) = MAX(zn0s(i,j),0.5_wp*hlp)
        zn0s(i,j) = MIN(zn0s(i,j),1e2_wp*hlp)
        zn0s(i,j) = MIN(zn0s(i,j),1e9_wp)
        zn0s(i,j) = MAX(zn0s(i,j),1e6_wp)
      ELSE
        ! Old constant n0s
        zn0s(i,j) = 8.0e5_wp
      ENDIF
      zcrim (i,j) = ccsrim*zn0s(i,j)
      zcagg (i,j) = ccsagg*zn0s(i,j)
      zbsdep(i,j) = ccsdep*SQRT(v0snow)
      zvz0s (i,j) = ccsvel*EXP(ccsvxp * LOG(zn0s(i,j)))

      IF (zvzs(i,j) == 0.0_wp) THEN
        zvzs(i,j) = zvz0s(i,j) * EXP (ccswxp * LOG (0.5_wp*zqsk(i,j))) * zrho1o2(i,j)
      ENDIF
    ENDDO

!CDIR NODEP,VOVERTAKE,VOB
    DO i1d = 1, ic2
      i = idx2(i1d)
      j = jdx2(i1d)

      IF (zvzr(i,j) == 0.0_wp) THEN
        zvzr(i,j) = zvz0r * EXP (zvzxp * LOG (0.5_wp*zqrk(i,j))) * zrho1o2(i,j)
      ENDIF
    ENDDO

  !----------------------------------------------------------------------------
  ! Section 3: 
  !----------------------------------------------------------------------------

!CDIR BEGIN COLLAPSE
    zeln7o8qrk   (:,:) = 0.0_wp
    zeln27o16qrk (:,:) = 0.0_wp
    zeln13o8qrk  (:,:) = 0.0_wp
!!    zeln3o16qrk  (:,:) = 0.0_wp
    zeln13o12qsk (:,:) = 0.0_wp
    zeln5o24qsk  (:,:) = 0.0_wp
    zeln2o3qsk   (:,:) = 0.0_wp

!FR old  
!   zcsdep       (:,:) = 3.2E-2_wp
    zcsdep       (:,:) = 3.367E-2_wp

    zcidep       (:,:) = 1.3E-5_wp
    zcslam       (:,:) = 1e10_wp

    scau         (:,:) = 0.0_wp
    scac         (:,:) = 0.0_wp
    snuc         (:,:) = 0.0_wp
    scfrz        (:,:) = 0.0_wp
    simelt       (:,:) = 0.0_wp
    sidep        (:,:) = 0.0_wp
    ssdep        (:,:) = 0.0_wp
    sdau         (:,:) = 0.0_wp
    srim         (:,:) = 0.0_wp
    sshed        (:,:) = 0.0_wp
    sicri        (:,:) = 0.0_wp
    srcri        (:,:) = 0.0_wp
    sagg         (:,:) = 0.0_wp
    siau         (:,:) = 0.0_wp
    ssmelt       (:,:) = 0.0_wp
    sev          (:,:) = 0.0_wp
    srfrz        (:,:) = 0.0_wp
!CDIR END

    ic1 = 0
    ic2 = 0
    ic3 = 0
    ic4 = 0
    ic5 = 0
    ic6 = 0
    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar

        qrg  = qr(i,j,k)
        qsg  = qs(i,j,k)
        qvg  = qv(i,j,k)
        qcg  = qc(i,j,k)
        qig  = qi(i,j,k)
        tg   =  t(i,j,k,nx)
        ppg  = p0(i,j,k) + pp(i,j,k,nx)
        rhog = rho(i,j,k)

        llqr = zqrk(i,j) > zqmin
        llqs = zqsk(i,j) > zqmin

        IF (llqr) THEN
!US reported by Thorsten Reinhardt: Multiplication with zrho1o2 was missing
          zpkr(i,j) = zqrk(i,j) * zvz0r * EXP (zvzxp * LOG (zqrk(i,j))) * zrho1o2(i,j)
        ELSE
          zpkr(i,j) = 0.0_wp
        ENDIF

        IF (llqs) THEN
!US reported by Thorsten Reinhardt: Multiplication with zrho1o2 was missing
          zpks(i,j) = zqsk (i,j) * zvz0s(i,j) * EXP (ccswxp * LOG (zqsk(i,j))) * zrho1o2(i,j)
        ELSE
          zpks(i,j) = 0.0_wp
        ENDIF

        zpkr(i,j)   = MIN( zpkr(i,j) , zzar(i,j) )
        zpks(i,j)   = MIN( zpks(i,j) , zzas(i,j) )

        zzar(i,j)   = zdtdh(i,j) * (zzar(i,j)-zpkr(i,j))
        zzas(i,j)   = zdtdh(i,j) * (zzas(i,j)-zpks(i,j))

        zimr(i,j)   = 1.0_wp / (1.0_wp + zvzr(i,j) * zdtdh(i,j))
        zims(i,j)   = 1.0_wp / (1.0_wp + zvzs(i,j) * zdtdh(i,j))

        zqrk(i,j)   = zzar(i,j)*zimr(i,j)
        zqsk(i,j)   = zzas(i,j)*zims(i,j)

        llqr = zqrk(i,j) > zqmin
        llqs = zqsk(i,j) > zqmin
        llqc =       qcg > zqmin
        llqi =       qig > zqmin

        IF (llqr) THEN
          ic1 = ic1 + 1
          idx1(ic1) = i
          jdx1(ic1) = j
        ENDIF
        IF (llqs) THEN
          ic2 = ic2 + 1
          idx2(ic2) = i
          jdx2(ic2) = j
        ENDIF
        IF (llqi .OR. llqs) THEN
          ic3 = ic3 + 1
          idx3(ic3) = i
          jdx3(ic3) = j
        ENDIF
        IF ( tg < zthet .AND. qvg >  8.E-6_wp &
                        .AND. qig <= 0.0_wp ) THEN
          ic4 = ic4 + 1
          idx4(ic4) = i
          jdx4(ic4) = j
        ENDIF
        IF (llqc) THEN
          ic5 = ic5 + 1
          idx5(ic5) = i
          jdx5(ic5) = j
        ENDIF
        IF (llqr .AND. qcg <= 0.0_wp) THEN
          ic6 = ic6 + 1
          idx6(ic6) = i
          jdx6(ic6) = j
        ENDIF
      ENDDO
    ENDDO

! ic1
!CDIR NODEP,VOVERTAKE,VOB
    DO i1d =1, ic1
      i = idx1(i1d)
      j = jdx1(i1d)

      qcg  = qc(i,j,k)
      qig  = qi(i,j,k)
      tg   =  t(i,j,k,nx)
      llqi =  qig > zqmin

      zlnqrk       = LOG (zqrk(i,j))
      IF ( qig+qcg > zqmin ) THEN
        zeln7o8qrk(i,j)   = EXP (x7o8   * zlnqrk)
      ENDIF
      IF ( tg < ztrfrz ) THEN
        zeln27o16qrk(i,j) = EXP (x27o16 * zlnqrk)
      ENDIF
      IF (llqi) THEN
        zeln13o8qrk(i,j)  = EXP (x13o8  * zlnqrk)
      ENDIF
!!      IF (qcg <= 0.0_wp ) THEN
!!        zeln3o16qrk(i,j)  = EXP (x3o16  * zlnqrk)
!!      ENDIF
    ENDDO

! ic2
!CDIR NODEP,VOVERTAKE,VOB
    DO i1d =1, ic2
      i = idx2(i1d)
      j = jdx2(i1d)

      qcg = qc(i,j,k)
      qig = qi(i,j,k)

      zlnqsk       = LOG (zqsk(i,j))
      IF (qig+qcg > zqmin) THEN
        zeln13o12qsk(i,j) = EXP (x13o12 * zlnqsk)
      ENDIF
      zeln5o24qsk(i,j)  = EXP (x5o24  * zlnqsk)
      zeln2o3qsk(i,j)   = EXP (x2o3   * zlnqsk)
    ENDDO

! ic3
!CDIR NODEP,VOVERTAKE,VOB
    DO i1d =1, ic3
      i = idx3(i1d)
      j = jdx3(i1d)

      tg   =   t(i,j,k,nx)
      ppg  =  p0(i,j,k) + pp(i,j,k,nx)
      rhog = rho(i,j,k)
      llqs = zqsk(i,j) > zqmin

      zqvsi  = fqvs( fpvsi(tg), ppg )
      zdvtp  = ccdvtp * EXP(1.94_wp * LOG(tg)) / ppg          
      zhi    = ccshi1*zdvtp*rhog*zqvsi/(tg*tg)
      hlp    = zdvtp / (1.0_wp + zhi)
      zcidep(i,j) = ccidep * hlp
      IF (llqs) THEN
        zcslam(i,j) = EXP(ccslxp * LOG(ccslam * zn0s(i,j) / zqsk(i,j) ))
        zcslam(i,j) = MIN(zcslam(i,j),1e15_wp)
        zcsdep(i,j) = 4.0_wp * zn0s(i,j) * hlp
      ENDIF
    ENDDO

    !--------------------------------------------------------------------------
    ! Section 4: Initialize the conversion rates with zeros in every layer
    !            Deposition nucleation for low temperatures
    !--------------------------------------------------------------------------

    ! Deposition nucleation for low temperatures below a threshold

! ic4
!CDIR NODEP,VOVERTAKE,VOB
    DO i1d = 1, ic4
      i = idx4(i1d)
      j = jdx4(i1d)

      qvg  =  qv(i,j,k)
      tg   =   t(i,j,k,nx)
      ppg  =  p0(i,j,k) + pp(i,j,k,nx)
      rhog = rho(i,j,k)

      zqvsi   = fqvs( fpvsi(tg), ppg )
      IF( qvg > zqvsi ) THEN
        znin  = MIN( fxna(tg), znimax )
        zsnuc = zmi0 * z1orhog(i,j) * znin * zdtr
        snuc(i,j) = zsnuc
      ENDIF
    ENDDO

    !--------------------------------------------------------------------------
    ! Section 5: Search for cloudy grid points with cloud water and
    !            calculation of the conversion rates involving qc
    !--------------------------------------------------------------------------

! ic5
!CDIR NODEP,VOVERTAKE,VOB
    DO i1d =1, ic5
      i = idx5(i1d)
      j = jdx5(i1d)

      qrg  =   qr(i,j,k)
      qvg  =   qv(i,j,k)
      qcg  =   qc(i,j,k)
      qig  =   qi(i,j,k)
      tg   =    t(i,j,k,nx)
      ppg  =   p0(i,j,k) + pp(i,j,k,nx)
      rhog =  rho(i,j,k)
      llqs = zqsk(i,j) > zqmin

      IF (iautocon == 0) THEN
        ! Kessler (1969) autoconversion rate
        zscau  = zccau * MAX( qcg - qc0, 0.0_wp )
        zscac  = zcac  * qcg * zeln7o8qrk(i,j)
      ELSEIF (iautocon == 1) THEN
        ! Seifert and Beheng (2001) autoconversion rate
        ! with constant cloud droplet number concentration cloud_num
        IF (qcg > 1e-6) THEN
#ifndef MESSY
          ztau   = MIN(1.0_wp-qcg/(qcg+qrg),0.9_wp)
#else
          ztau   = MAX(MIN(1.0_wp-qcg/(qcg+qrg),0.9_wp) , 0._wp)
#endif
          zphi   = zkphi1 * ztau**zkphi2 * (1.0_wp - ztau**zkphi2)**3
          zscau  = zconst * qcg*qcg*qcg*qcg &
            &    * (1.0_wp + zphi/(1.0_wp - ztau)**2)
          zphi   = (ztau/(ztau+zkphi3))**4
          zscac  = zkcac * qcg * qrg * zphi !* zrho1o2(i,j)
        ELSE
          zscau  = 0.0_wp
          zscac  = 0.0_wp
        ENDIF
      ENDIF
      IF (llqs) THEN
        zscrim = zcrim(i,j) * EXP(ccsaxp * LOG(zcslam(i,j))) * qcg !* zrho1o2(i,j)
      ELSE
        zscrim = 0.0_wp
      ENDIF

      zscshe = 0.0_wp
      IF( tg >= t0_melt ) THEN
        zscshe = zscrim
        zscrim = 0.0_wp
      ENDIF
      ! Check for maximum depletion of cloud water and adjust the
      ! transfer rates accordingly
      zscmax = qcg*zdtr 
      zscsum = zscau + zscac + zscrim + zscshe 
      zcorr  = zscmax / MAX( zscmax, zscsum )
      IF( tg <= zthn ) THEN
        scfrz(i,j) = zscmax
      ELSE
        scau (i,j) = zcorr*zscau
        scac (i,j) = zcorr*zscac
        srim (i,j) = zcorr*zscrim
        sshed(i,j) = zcorr*zscshe 
      ENDIF

      ! Calculation of heterogeneous nucleation of cloud ice.
      ! This is done in this section, because we require water saturation
      ! for this process (i.e. the existence of cloud water) to exist.
      ! Hetrogeneous nucleation is assumed to occur only when no
      ! cloud ice is present and the temperature is below a nucleation
      ! threshold.
      IF( tg <= 267.15_wp .AND. qig <= 0.0_wp ) THEN
        znin  = MIN( fxna(tg), znimax )
        zsnuc = zmi0 * z1orhog(i,j) * znin * zdtr
        snuc(i,j) = zsnuc
      ENDIF
      ! Calculation of in-cloud rainwater freezing
      IF ( tg < ztrfrz ) THEN
        zsrfrz = zcrfrz*SQRT( (ztrfrz-tg)**3 )* zeln27o16qrk(i,j)
        srfrz(i,j) = zsrfrz
      ENDIF
    ENDDO
   
      !------------------------------------------------------------------------
      ! Section 6: Search for cold grid points with cloud ice and/or snow and
      !            calculation of the conversion rates involving qi and ps
      !------------------------------------------------------------------------

! also ic3
!CDIR NODEP,VOVERTAKE,VOB
    DO i1d =1, ic3
      i = idx3(i1d)
      j = jdx3(i1d)

      qvg  =  qv(i,j,k)
      qig  =  qi(i,j,k)
      qsg  =  qs(i,j,k)
      tg   =   t(i,j,k,nx)
      ppg  =  p0(i,j,k) + pp(i,j,k,nx)
      rhog = rho(i,j,k)
      llqi =  qig > zqmin

      IF (tg<=t0_melt) THEN
        zqvsi   = fqvs( fpvsi(tg), ppg )
        znin    = MIN( fxna(tg), znimax )
        zmi     = MIN( rhog*qig/znin, zmimax )
        zmi     = MAX( zmi0, zmi )
        zsvmax  = (qvg - zqvsi) * zdtr
        zsagg   = zcagg(i,j) * EXP(ccsaxp*LOG(zcslam(i,j))) * qig
        zsagg   = MAX( zsagg, 0.0_wp ) & !* zrho1o2(i,j) &
          * MAX(0.2_wp,MIN(EXP(0.09_wp*(tg-t0_melt)),1.0_wp))
        znid      = rhog * qig/zmi
        IF (llqi) THEN
          zlnlogmi= LOG (zmi)
          zsidep    = zcidep(i,j) * znid * EXP(0.33_wp * zlnlogmi)   &
                        * ( qvg - zqvsi )
        ELSE
          zsidep = 0.0_wp
        ENDIF
        zsvidep   = 0.0_wp
        zsvisub   = 0.0_wp
        zsimax    = qig*zdtr 
        IF( zsidep > 0.0_wp ) THEN
          zsvidep = MIN( zsidep, zsvmax )
        ELSEIF (zsidep < 0.0_wp ) THEN
          zsvisub = - MAX(-zsimax, zsvmax )
        ENDIF
        zsiau = zciau * MAX( qig - qi0, 0.0_wp ) &
          * MAX(0.2_wp,MIN(EXP(0.09_wp*(tg-t0_melt)),1.0_wp))
        IF (llqi) THEN
          zlnlogmi = LOG(zmsmin/zmi)
          zztau    = 1.5_wp*( EXP(0.66_wp*zlnlogmi) - 1.0_wp)
          zsdau    = zsvidep/MAX(zztau,eps_div)
        ELSE
          zsdau    =  0.0_wp
        ENDIF
        zsicri    = zcicri * qig * zeln7o8qrk(i,j)
        zsrcri    = zcrcri * (qig/zmi) * zeln13o8qrk(i,j)
        zxfac     = 1.0_wp + zbsdep(i,j) * EXP(ccsdxp*LOG(zcslam(i,j)))
        zssdep    = zcsdep(i,j) * zxfac * ( qvg - zqvsi ) / (zcslam(i,j)+eps_div)**2

        ! Check for maximal depletion of vapor by sdep
        IF (zssdep > 0.0_wp) zssdep = MIN(zssdep, zsvmax-zsvidep)
        ! Check for maximal depletion of snow by sdep
        IF (zssdep < 0.0_wp) zssdep = MAX(zssdep, -qsg*zdtr)

        zsisum = zsiau + zsdau + zsagg + zsicri + zsvisub
        zcorr  = 0.0_wp
        IF( zsimax > 0.0_wp ) zcorr  = zsimax / MAX( zsimax, zsisum )
        sidep(i,j)  = zsvidep - zcorr*zsvisub
        sdau (i,j)  = zcorr*zsdau
        siau (i,j)  = zcorr*zsiau
        sagg (i,j)  = zcorr*zsagg
        ssdep(i,j)  = zssdep
        srcri(i,j)  = zsrcri
        sicri(i,j)  = zcorr*zsicri

      !------------------------------------------------------------------------
      ! Section 7: Search for warm grid points with cloud ice and/or snow and
      !            calculation of the melting rates of qi and ps
      !------------------------------------------------------------------------

      ELSE ! tg > 0
        simelt(i,j) = qig*zdtr
        zqvsw0      = fqvs( zpvsw0, ppg)
        zx1         = (tg - t0_melt) + zasmel*(qvg - zqvsw0)
        zx2         = 1.0_wp + zbsmel * zeln5o24qsk(i,j)
        zssmelt     = zcsmel * zx1 * zx2 * zeln2o3qsk(i,j)
        ssmelt(i,j) = MAX( zssmelt, 0.0_wp )
      ENDIF ! tg
    ENDDO

    !--------------------------------------------------------------------------
    ! Section 8: Search for grid points with rain in subsaturated areas
    !            and calculation of the evaporation rate of rain
    !--------------------------------------------------------------------------

! ic6
!CDIR NODEP,VOVERTAKE,VOB
    DO i1d =1, ic6
      i = idx6(i1d)
      j = jdx6(i1d)

      qvg = qv(i,j,k)
      tg  =  t(i,j,k,nx)
      ppg = p0(i,j,k) + pp(i,j,k,nx)

      zlnqrk      = LOG (zqrk(i,j))
      zqvsw       = fqvs( fpvsw(tg), ppg )
      zx1         = 1.0_wp + zbev * EXP (zbevxp  * zlnqrk)
      sev(i,j)    = zcev*zx1*(zqvsw - qvg) * EXP (zcevxp  * zlnqrk)
!      zqvsw    = fqvs( fpvsw(tg), ppg )
!      zx1      = 1.0_wp + zbev* zeln3o16qrk(i,j)
!      zsev     = zcev*zx1*(zqvsw - qvg)*SQRT(zqrk(i,j))
!      sev(i,j) = MAX( zsev, 0.0_wp )

      ! Calculation of below-cloud rainwater freezing
      IF ( tg < ztrfrz ) THEN
        zsrfrz = zcrfrz*SQRT( (ztrfrz-tg)**3 ) * zeln27o16qrk(i,j)
        srfrz(i,j)  = zsrfrz
      ENDIF
    ENDDO

    !--------------------------------------------------------------------------
    ! Section 9: Calculate the total tendencies of the prognostic variables.
    !            Update the prognostic variables in the interior domain.
    !--------------------------------------------------------------------------

    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar

        qrg = qr(i,j,k)
        qsg = qs(i,j,k)
        qig = qi(i,j,k)
        rhog = rho(i,j,k)
        ppg = p0(i,j,k) + pp(i,j,k,nx)

        zsrmax = zzar(i,j)*z1orhog(i,j)*zdtr
        zssmax = zzas(i,j)*z1orhog(i,j)*zdtr
        zsrsum = sev(i,j) + srfrz(i,j) + srcri(i,j)
        zcorr  = 1.0_wp
        IF(zsrsum > 0) THEN
          zcorr  = zsrmax / MAX( zsrmax, zsrsum )
        ENDIF
        sev  (i,j) = zcorr*sev(i,j)
        srfrz(i,j) = zcorr*srfrz(i,j)
        srcri(i,j) = zcorr*srcri(i,j)

        ssmelt(i,j) = MIN(ssmelt(i,j), zssmax)
        IF (ssdep(i,j) < 0.0_wp ) THEN
          ssdep(i,j) = MAX(ssdep(i,j), - zssmax)
        ENDIF
        zqvt = sev(i,j)   - sidep(i,j) - ssdep(i,j)  - snuc(i,j)
        zqct = simelt(i,j)- scau(i,j)  - scfrz(i,j)  - scac(i,j)   - sshed(i,j) - srim(i,j) 
        zqit = snuc(i,j)  + scfrz(i,j) - simelt(i,j) - sicri(i,j)  + sidep(i,j) - sdau(i,j)  - sagg(i,j) - siau(i,j)
        zqrt = scau(i,j)  + sshed(i,j) + scac(i,j)   + ssmelt(i,j) - sev(i,j)   - srcri(i,j) - srfrz(i,j) 
        zqst = siau(i,j)  + sdau(i,j)  + sagg(i,j)   - ssmelt(i,j) + sicri(i,j) + srcri(i,j) + srim(i,j)       &
                                                                                + ssdep(i,j) + srfrz(i,j)
        ztt = cpdr*( lh_v*(zqct+zqrt) + lh_s*(zqit+zqst) )

        IF (lsppt)THEN
          IF(ntstep>0.AND.i>=istart.and.i<=iend.and.j>=jstart.and.j<=jend) THEN
            CALL apply_tqx_tend_adj(itype_qxpert_rn,itype_qxlim_rn,ppg, &
                t(i,j,k,nx),qv(i,j,k), qc(i,j,k),qi(i,j,k),qr(i,j,k),qs(i,j,k),&
                 pertstoph(i,j,k),ztt,zqvt,zqct,zqit,zqrt,zqst,ldum)
          ENDIF
        ENDIF

        ! Update variables and add qi to qrs for water loading 
        qig = MAX ( 0.0_wp, qig + zqit*zdt)
        qrg = MAX ( 0.0_wp, (zzar(i,j)*z1orhog(i,j) + zqrt*zdt)*zimr(i,j))
        qsg = MAX ( 0.0_wp, (zzas(i,j)*z1orhog(i,j) + zqst*zdt)*zims(i,j))

#ifdef MESSY
        ! convert kg/(kg*s) into kg/(m2*s)
        precmelt_ls(i,j,k)   = (ssmelt(i,j) + simelt(i,j)) &
                               * rhog * (hhl(i,j,k)-hhl(i,j,k+1))
        ! rain formation in kg/kg
        rainform_bave(i,j,k) = zqrt * zdt
        ! snow formation in kg/kg
        snowform_bave(i,j,k) = zqst * zdt
        ! rain flux  in kg/(m2*s)
        precflx_ls(i,j,k)  = zqrt  * rhog * (hhl(i,j,k)-hhl(i,j,k+1))
        ! snow flux in  kg/(m2*s)
        snowflx_ls(i,j,k)  = zqst  * rhog * (hhl(i,j,k)-hhl(i,j,k+1))
        ! rain flux without in-cloud production kg/(m2*s)
        ! in-coming rain + melting snow - evaporation - freezing of rain
        precflxno_ls(i,j,k)  = zprvr(i,j)                               &
                               + (ssmelt(i,j) - sev(i,j) - srfrz(i,j))  &
                               * rhog * (hhl(i,j,k)-hhl(i,j,k+1))
        ! snow flux without in-cloud production kg/(m2*s)
        ! incoming snow + freezing of rain - sublimation - melting snow
        snowflxno_ls(i,j,k)  = zprvs(i,j)                               &
                               + (srfrz(i,j) - ssdep(i,j)- ssmelt(i,j)) &
                               * rhog * (hhl(i,j,k)-hhl(i,j,k+1))
        ! make a crude assumtion about the precipitating cloud cover
        ! if rain or snow exists, the cover is 1
        IF (qrg > 1.e-15_wp .OR. qsg > 1.e-15_wp) THEN
           preccover_ls(i,j,k) = 1._wp
        ELSE
           preccover_ls(i,j,k) = 0._wp
        ENDIF
        sediice_ls(i,j,k) = 0._wp
#endif

        !----------------------------------------------------------------------
        ! Section 10: Complete time step
        !----------------------------------------------------------------------

        IF ( k /= ke) THEN
          ! Store precipitation fluxes and sedimentation velocities 
          ! for the next level
          zprvr(i,j) = qrg*rhog*zvzr(i,j)
          zprvs(i,j) = qsg*rhog*zvzs(i,j)
          IF (zprvr(i,j) <= zqmin) zprvr(i,j)=0.0_wp
          IF (zprvs(i,j) <= zqmin) zprvs(i,j)=0.0_wp

#ifdef NUDGING
          ! for the latent heat nudging
          IF ((llhn .OR. llhnverif) .AND. lhn_qrs) THEN
            qrsflux(i,j,k) = zprvr(i,j)+zprvs(i,j)
            qrsflux(i,j,k) = 0.5_wp*(qrsflux(i,j,k)+zpkr(i,j)+zpks(i,j))
          ENDIF
#endif

          IF (qrg+qr(i,j,k+1) <= zqmin) THEN
            zvzr(i,j)= 0.0_wp
          ELSE
            zvzr(i,j)= zvz0r * EXP(zvzxp*LOG((qrg+qr(i,j,k+1))*0.5_wp*rhog)) * zrho1o2(i,j)
          ENDIF
          IF (qsg+qs(i,j,k+1) <= zqmin) THEN
            zvzs(i,j)= 0.0_wp
          ELSE
            zvzs(i,j)= zvz0s(i,j) * EXP(zv1s/(zbms+1.0_wp)*LOG((qsg+qs(i,j,k+1))*0.5_wp*rhog)) * zrho1o2(i,j)
          ENDIF
        ELSE
          ! Precipitation fluxes at the ground
          prr_gsp(i,j) = 0.5_wp * (qrg*rhog*zvzr(i,j) + zpkr(i,j))
          prs_gsp(i,j) = 0.5_wp * (qsg*rhog*zvzs(i,j) + zpks(i,j))

#ifdef NUDGING
          ! for the latent heat nudging
          IF ((llhn .OR. llhnverif) .AND. lhn_qrs)        &
            qrsflux(i,j,k) = prr_gsp(i,j)+prs_gsp(i,j)
#endif

        ENDIF

        ! Update of prognostic variables or tendencies
        qr (i,j,k   ) = qrg
        qs (i,j,k   ) = MAX ( 0.0_wp, qsg )
        qi (i,j,k   ) = qig
        qrs(i,j,k   ) = qrg+qsg+qig
        t  (i,j,k,nx) = t (i,j,k,nx) + ztt*zdt 
        qv (i,j,k   ) = MAX ( 0.0_wp, qv(i,j,k) + zqvt*zdt )
        qc (i,j,k   ) = MAX ( 0.0_wp, qc(i,j,k) + zqct*zdt )
      ENDDO
    ENDDO

  IF (izdebug > 10) THEN
    ! Check for negative values
    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar
        IF (qr(i,j,k) < 0) THEN
          WRITE (*,'(a)') ' WARNING: hydci_pp, negative value in qr'
        ENDIF
        IF (qc(i,j,k) < 0) THEN
          WRITE (*,'(a)') ' WARNING: hydci_pp, negative value in qc'
        ENDIF
        IF (qi(i,j,k) < 0) THEN
          WRITE (*,'(a)') ' WARNING: hydci_pp, negative value in qi'
        ENDIF
        IF (qs(i,j,k) < 0) THEN
          WRITE (*,'(a)') ' WARNING: hydci_pp, negative value in qs'
        ENDIF
        IF (qv(i,j,k) < 0) THEN
          WRITE (*,'(a)') ' WARNING: hydci_pp, negative value in qv'
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  ! Do a final saturation adjustment for new values of t, qv and qc
!CDIR COLLAPSE
  zpres(:,:) = p0(:,:,k) + pp(:,:,k,nx)

  CALL satad ( 1, t(:,:,k,nx), qv(:,:,k),                 &
               qc(:,:,k), t(:,:,k,nx), zpres,             &
               zdummy(:,:,1),zdummy(:,:,2),zdummy(:,:,3), &
               zdummy(:,:,4),zdummy(:,:,5),zdummy(:,:,6), &
               zdummy(:,:,7),zdummy(:,:,8),               &
               b1, b2w, b3, b4w, b234w, rdv, o_m_rdv,     &
               rvd_m_o, lh_v, cpdr, cp_d,                 &
               ie, je, istartpar, iendpar, jstartpar, jendpar )

  IF ( ldiabf_lh ) THEN
    ! compute temperature increment due to latent heat
!CDIR COLLAPSE
    tinc_lh(:,:,k) = tinc_lh(:,:,k) + t(:,:,k,nx)
  ENDIF

ENDDO loop_over_levels

#ifdef NUDGING
! add part of latent heating calculated in subroutine hydci to model latent
! heating field: add temperature to model latent heating field
IF (llhn .OR. llhnverif) &
 CALL get_gs_lheating ('inc',1,ke)
#endif

!------------------------------------------------------------------------------
! End of module procedure hydci_pp
!------------------------------------------------------------------------------

END SUBROUTINE hydci_pp

!==============================================================================

!option! -pvctl ifopt
!option! -pvctl _on_adb
SUBROUTINE hydci_pp_gr

!------------------------------------------------------------------------------
! Description:
!   This module procedure calculates the rates of change of temperature, cloud
!   water, cloud ice, water vapor, rain, snow, and graupel due to cloud
!   microphysical processes related to the formation of grid scale
!   precipitation.
!   The variables are updated in this subroutine. Rain, snow, and graupel are
!   prognostic variables. The precipitation fluxes at the surface are stored
!   on the corresponding global fields.
!   The subroutine relies mostly on conversion terms used in hydci_pp.
!   In contrast to hydci_pp, qc0 = 0.0002 (instead of 0.0) is used!
!   This is version G29TK21.
!
! Method:
!   The sedimentation of rain, snow, and graupel is computed implicitly.
!
! Vectorization:
!   Most computations in this routine are grouped in IF-clauses. But the IFs
!   inside DO-loops often hinder or even disables vectorization. 
!   For the big IF-chunks, the condition is now checked at the beginning of 
!   the subroutine and the corresponding indices are stored in extra index
!   arrays. The following loops then are running only over these indices,
!   avoiding at least some IF-clauses inside the loops.
!
!------------------------------------------------------------------------------
! Declarations:
!
!------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------
!
! Subroutine arguments: None
! --------------------
!
! Local parameters:
! ----------------

  INTEGER (KIND=iintegers), PARAMETER ::  &
    iautocon       = 1, &
    isnow_n0temp   = 2

  REAL    (KIND=wp   ),     PARAMETER ::  &
    zqmin = 1.0E-15_wp,& ! threshold for computations

    ! Parameters for autoconversion of cloud water and cloud ice
    zccau = 4.0E-4_wp, & ! autoconversion coefficient (cloud water to rain)
    zciau = 1.0E-3_wp, & ! autoconversion coefficient (cloud ice   to snow)

    zkcau  = 9.44e+09_wp, & ! kernel coeff for SB2001 autoconversion
    zkcac  = 5.25e+00_wp, & ! kernel coeff for SB2001 accretion
    zcnue  = 2.00e+00_wp, & ! gamma exponent for cloud distribution
    zxstar = 2.60e-10_wp, & ! separating mass between cloud and rain
    zkphi1 = 6.00e+02_wp, & ! constant in phi-function for autoconversion
    zkphi2 = 0.68e+00_wp, & ! exponent in phi-function for autoconversion
    zkphi3 = 5.00e-05_wp, & ! exponent in phi-function for accretion

    ! Other process coefficients.
    ! These coefficients have been precalculated on the basis of the
    ! following basic parameters:
    !
    ! N0R = 8.0 E6 : Parameter in the size distribution function for rain
    ! N0S = 8.0 E5 : Parameter in the size distribution function for snow
    ! N0G = 4.0 E6 : Parameter in the size distribution function for graupel
    ! ECR = 0.8    : Collection efficiency for rain collecting cloud water
    ! ECS = 0.9    : Collection efficiency for snow/gr. collecting cloud water
    ! EIS = 0.5    : Collection efficiency for snow/gr. collecting cloud ice
    ! EIR = 0.8    : Collection efficiency for rain collecting cloud ice
    ! V0R = 130.0  : Factor in the terminal velocity for raindrops
    !                VTR(D) = V0R*D**(1/2)
    ! V0S = 4.9    : Factor in the terminal velocity for snow particles
    !                VTS(DS) = V0S*DS**(1/4)
    ! V0G = 442.0  : Factor in the terminal velocity for graupel particles
    !                VTG(DG) = V0G*DG**(0.89)
    ! AMS = 0.038  : Formfactor in the mass-size relation of snow particles
    !                m(DS) = AMS*DS**2
    ! AMG = 169.6  : Formfactor in the mass-size relation of graupel particles
    !                m(DG) = AMG*DG**(3.1)
    ! AMI = 0.038  : Formfactor in the mass-size relation of cloud ice crystals
    !                m(DI) = AMI*DI**3
    ! ETA  = : Viscosity of air (= 1.72e-5*393.0*(ztx/273.0)**1.5/(ztx+120.0))
    ! DIFF = 2.22 E-5 : Molecular diffusion coefficient for water vapour
    ! For calculation of fits for deposition/melting, a temperature/pressure
    ! dependent diffusion coefficient formula is used:
    ! DIFF(ztx,zp)=pnull/zp*(2.22e-5+1.46e-7*(ztx-t0_melt))  (t>t0_melt)!L.-B.
    ! DIFF(ztx,zp)=pnull/zp*(2.22e-5+1.25e-7*(ztx-t0_melt))  (t<t0_melt)!L.-B.
    ! (DIFF formula taken from Landolt-Boernstein), pnull=100000.0
    ! LHEAT= 2.40 E-2 : Heat conductivity of air
    ! RHOW = 1.0 E3   : Density of water
    ! lh_v            : latent heat of vapourization
    ! lh_f            : latent heat of fusion
    ! AR  = PI*RHOW*N0R
    ! AS  = 2*N0S*AMS
    ! HW              : Howell factor ( =(1/(1+H_w)) or =(1/(1+H_i)))

    zhw   = 2.270603_wp,     & ! Howell factor
    zecs  = 0.9_wp,   & ! Collection efficiency for snow collecting cloud water

    zadi  = 0.217_wp, & ! Formfactor in the size-mass relation of ice particles
    zbdi  = 0.302_wp, & ! Exponent in the size-mass relation of ice particles
    zams  = 0.038_wp, & ! Formfactor in the mass-size relation of snow particles
    zbms  = 2.000_wp, & ! Exponent in the mass-size relation of snow particles

    zv1s  = 0.50_wp,  & ! Exponent in the terminal velocity for snow

    zami  = 130.0_wp, & ! Formfactor in the mass-size relation of cloud ice
    zn0s0 = 8.0E5_wp, & !
    zn0s1 = 13.5_wp * 5.65E5_wp, & ! parameter in N0S(T)
    zn0s2 = -0.107_wp, & ! parameter in N0S(T), Field et al

    zcsg=0.5_wp,        & !coefficient for snow-graupel conversion by riming
    zcrim_g=4.43_wp,    & !
    zrimexp_g=0.94878_wp, &

    zrho0 = 1.225e+0_wp, & ! reference air density
    zrhow = 1.000e+3_wp, & ! density of liquid water

    zcagg_g = 2.46_wp , & !
    zcac  = 1.72_wp   , & ! (15/32)*(PI**0.5)*(ECR/RHOW)*V0R*AR**(1/8)
    zcicri= 1.72_wp   , & ! (15/32)*(PI**0.5)*(EIR/RHOW)*V0R*AR**(1/8)
    zcrcri= 1.24E-3_wp    ! (PI/24)*EIR*V0R*Gamma(6.5)*AR**(-5/8)

! These are now in data_gscp
!  REAL    (KIND=wp   )      ::  &
!    ! these are no parameters any more but depend on mu_rain
!    zcev                  , & ! 2*PI*DIFF*HW*N0R*AR**(-1/2)
!    zbev                  , & ! 0.26*sqrt(0.5*RHO*v0r/eta)*Gamma(2.75)*AR**(-3/16)
!    zcevxp                , & !
!    zbevxp                , & !
!    zvzxp                 , & !
!    zvz0r                     ! coefficient of sedimentation velocity for rain

  REAL    (KIND=wp   ),     PARAMETER ::  &
    zasmel= 2.95E3_wp , & ! DIFF*lh_v*RHO/LHEAT
    zcrfrz= 1.68_wp   , & ! coefficient for raindrop freezing
    zexpsedg=0.217_wp,  & ! exponent for graupel sedimentation

    zdv    = 2.22e-5_wp, & ! molecular diffusion coefficient for water vapour
    zlheat = 2.40E-2_wp, & ! thermal conductivity of dry air
    zeta   = 1.75e-5_wp, & ! kinematic viscosity of air

    ! Additional parameters
    zthn  = 236.15_wp , & ! temp. for hom. freezing of cloud/rain water
    ztrfrz= 271.15_wp , & ! threshold temperature for heterogeneous
                              ! freezing of raindrops
    zmi0  = 1.0E-12_wp, & ! initial crystal mass for cloud ice nucleation
    zmimax= 1.0E-9_wp , & ! maximum mass of cloud ice crystals
    zmsmin= 3.0E-9_wp , & ! initial mass of snow crystals
    zvz0g = 12.24_wp  , & ! coefficient of sedimentation velocity for graupel
    ztcrit=3339.5_wp  , & ! factor in calculation of critical temperature

    ! Constant exponents in the transfer rate equations
    x1o12  =  1.0_wp/12.0_wp  ,  x3o16  =  3.0_wp/16.0_wp, &
    x7o8   =  7.0_wp/ 8.0_wp  ,                                    &
    x1o8   =  1.0_wp/ 8.0_wp  ,  x1o3   =  1.0_wp/ 3.0_wp, &
    x13o8  = 13.0_wp/ 8.0_wp  ,  x13o12 = 13.0_wp/12.0_wp, &
    x27o16 = 27.0_wp/16.0_wp  ,  x3o4   =  0.75_wp,            &
    x1o2   =  1.0_wp/ 2.0_wp

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    izdebug,           & ! for debug output
    k, i, j,           & ! loop indices
    nx,                & ! time-level for computation
    nn

  REAL    (KIND=wp   )     ::  &
    fpvsw, fpvsi, fqvs,& ! name of statement functions
    fxna ,             & ! statement function for ice crystal number
    ztx  , zpv  , zpx ,& ! dummy arguments for statement functions
    znimax,            & ! maximum number of cloud ice crystals
    zpvsw0,            & ! sat.vap. pressure at melting temperature
    zqvsw0,            & ! sat.specific humidity at melting temperature
    zqvsw0diff,        & ! qv-zqvsw0
    zdt,               & ! timestep for cloud physics
    zdtr ,             & ! reciprocal of timestep for cloud physics
    zscsum, zscmax, zcorr,  & ! terms for limiting  total cloud water depletion
    zsrsum,            & ! terms for limiting total rain water depletion
    znin,              & ! number of cloud ice crystals at nucleation
    znid,              & ! number of cloud ice crystals for deposition
    zmi ,              & ! mass of a cloud ice crystal
    zsvidep, zsvisub,  & ! deposition, sublimation of cloud ice
    zsimax , zsisum , zsvmax,& ! terms for limiting total cloud ice depletion
    zqvsi, zqvsw,      & ! sat. specitic humidity at ice and water saturation
    zqvsidiff,         & ! qv-zqvsi
    ztfrzdiff,         & ! ztrfrz-t
    zztau, zx1, ztt,   & ! help variables
    ztau, zphi, zhi, zdvtp, ztc, zxfac, &
    zeff      ! a lot more help variables

  REAL    (KIND=wp   )     ::  &
    zqct   ,& ! layer tendency of cloud water
    zqvt   ,& ! layer tendency of water vapour
    zqit   ,& ! layer tendency of cloud ice
    zqrt   ,& ! layer tendency of rain
    zqgt   ,& ! layer tendency of graupel
    zqst      ! layer tendency of snow

  REAL (KIND=wp)             ::       &
    mma(10), mmb(10)

  REAL    (KIND=wp   )     ::  &
    zdh, zlnqrk, &
    zln1o2, zlnqsk, zlnqgk, zlnlogmi, zar, zn0r, &
    qcg, tg, qvg, qrg, qsg, qgg, qig, rhog, ppg, alf,bet,m2s,m3s,hlp

  LOGICAL ::  &
    llqs,llqc,llqi,llqg     !  switches for existence of qr, qs, qc, qi, qg

  LOGICAL ::  &
    llqr(ie,je)

! Local (automatic) arrays:
! -------------------------
  REAL    (KIND=wp   )     ::  &
    zpres (ie,je),     & ! pressure
    zvzr  (ie,je),     & !
    zvzs  (ie,je),     & !
    zvzg  (ie,je),     & !
    zqrk  (ie,je),     & !
!   zlnqrk(ie,je),     & !
    zqsk  (ie,je),     & !
    zqgk  (ie,je),     & !
    zpkr  (ie,je),     & !
    zpks  (ie,je),     & !
    zpkg  (ie,je),     & !
    zprvr (ie,je),     & !
    zprvs (ie,je),     & !
    zprvg (ie,je),     & !
    zdtdh (ie,je),     & !
    zn0s  (ie,je),     & !
    zdummy(ie,je,8)

  REAL    (KIND=wp   )     ::  &
    scau  (ie,je),     & ! transfer rate due to autoconversion of cloud water
    scac  (ie,je),     & ! transfer rate due to accretion of cloud water
    snuc  (ie,je),     & ! transfer rate due to nucleation of cloud ice
    scfrz (ie,je),     & ! transfer rate due to homogeneous freezing of cloud water
    simelt(ie,je),     & ! transfer rate due to melting of cloud ice
    sidep (ie,je),     & ! transfer rate due to depositional growth of cloud ice
    ssdep (ie,je),     & ! transfer rate due to depositional growth of snow
    sgdep (ie,je),     & ! transfer rate due to depositional growth of graupel
    sdau  (ie,je),     & ! transfer rate due to depositional cloud ice autoconversion
    srim  (ie,je),     & ! transfer rate due to riming of snow
    srim2 (ie,je),     & ! transfer rate due to riming of graupel
    sconsg(ie,je),     & ! transfer rate due to conversion from snow to graupel by riming
    sshed (ie,je),     & ! transfer rate due to shedding
    sicri (ie,je),     & ! transfer rate due to cloud ice collection by rain (sink qi)
    srcri (ie,je),     & ! transfer rate due to cloud ice collection by rain (sink qr)
    sagg  (ie,je),     & ! transfer rate due to aggregation of snow and cloud ice
    sagg2 (ie,je),     & ! transfer rate due to aggregation of graupel and cloud ice
    siau  (ie,je),     & ! transfer rate due to autoconversion of cloud ice
    ssmelt(ie,je),     & ! transfer rate due to melting of snow
    sgmelt(ie,je),     & ! transfer rate due to melting of graupel
    sev   (ie,je),     & ! transfer rate due to evaporation of rain
    sconr (ie,je),     & ! transfer rate due to condensation on melting snow/graupel
    srfrz (ie,je)        ! transfer rate due to rainwater freezing

  REAL    (KIND=wp   )     ::  &
    zeln7o8qrk  (ie,je), & !
    zeln27o16qrk(ie,je), & !
    zeln13o8qrk (ie,je), & !
    zeln3o4qsk  (ie,je), & !
    zeln8qsk    (ie,je), & !
    zeln6qgk    (ie,je), & !
    zelnrimexp_g(ie,je), & !
    zcrim       (ie,je), & !
    zcagg       (ie,je), & !
    zbsdep      (ie,je), & !
    z1orhog     (ie,je), & ! 1/rhog
    zrho1o2     (ie,je), & ! (rho0/rhog)**1/2
    zvz0s       (ie,je), & !
    zzar        (ie,je), & !
    zzas        (ie,je), & !
    zzag        (ie,je), & !
    zimr        (ie,je), & !
    zims        (ie,je), & !
    zimg        (ie,je), & !
    zcsdep      (ie,je), & !
    zcidep      (ie,je), & !
    zcslam      (ie,je), & !
    zsrmax      (ie,je), & !
    zssmax      (ie,je), & !
    zsgmax      (ie,je)

!NEC_CB
  REAL    (KIND=wp   )     ::  &
    helpr(ie*je), helps(ie*je), helpg(ie*je)

  ! Integer arrays for a better vectorization
  INTEGER (KIND=iintegers) ::  &
    idx1(ie*je),       jdx1(ie*je),   & !
    idx2(ie*je),       jdx2(ie*je),   & !
    idx3(ie*je),       jdx3(ie*je),   & !
    idx4(ie*je),       jdx4(ie*je),   & !
    idx5(ie*je),       jdx5(ie*je),   & !
    idx6(ie*je),       jdx6(ie*je),   & !
    idx7(ie*je),       jdx7(ie*je),   & !
    idx8(ie*je),       jdx8(ie*je)

  ! Dimensions and loop counter for storing the indices
  INTEGER (KIND=iintegers) ::  &
    ic1, ic2, ic3, ic4, ic5, ic6, ic7, ic8, i1d

! Tracer pointers
!----------------
  REAL (KIND=wp),     POINTER :: &
    qv  (:,:,:) => NULL(),         & ! QV at nx
    qc  (:,:,:) => NULL(),         & ! QC at nx
    qi  (:,:,:) => NULL(),         & ! QI at nx
    qg  (:,:,:) => NULL(),         & ! QG at nx
    qr  (:,:,:) => NULL(),         & ! QR at nx
    qs  (:,:,:) => NULL()            ! QS at nx

  INTEGER (KIND=iintegers)    :: izerror
  CHARACTER (LEN=255)         :: yzerrmsg
    
  CHARACTER (LEN=25)          :: yzroutine = 'hydci_pp_gr'

  LOGICAL :: ldum

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine hydci_pp_gr
!------------------------------------------------------------------------------

! Statement functions
! -------------------

! saturation vapour pressure over water (fpvsw), over ice (fpvsi)
! and specific humidity at vapour saturation (fqvs)
  fpvsw(ztx)     = b1*EXP( b2w*(ztx-b3)/(ztx-b4w) )
  fpvsi(ztx)     = b1*EXP( b2i*(ztx-b3)/(ztx-b4i) )
  fqvs (zpv,zpx) = rdv*zpv/( zpx - o_m_rdv*zpv )

! Number of activate ice crystals;  ztx is temperature
  fxna(ztx)   = 1.0E2_wp * EXP(0.2_wp * (t0_melt - ztx))

! Coeffs for moment relation based on 2nd moment (Field 2005)
  mma = (/   5.065339_wp, -0.062659_wp, -3.032362_wp, 0.029469_wp, -0.000285_wp, &
             0.312550_wp,  0.000204_wp,  0.003199_wp, 0.000000_wp, -0.015952_wp /)
  mmb = (/   0.476221_wp, -0.015896_wp,  0.165977_wp, 0.007468_wp, -0.000141_wp, &
             0.060366_wp,  0.000079_wp,  0.000594_wp, 0.000000_wp, -0.003577_wp /)

!------------------------------------------------------------------------------
!  Section 1: Initial setting of local and global variables
!------------------------------------------------------------------------------

IF (ldebug_gsp) THEN
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

#ifdef NECSX
  CALL collapse(.TRUE., ie, istartpar, iendpar, jstartpar, jendpar)
#endif

  IF (ntstep == nstart) THEN
    zconst = zkcau / (20.0_wp*zxstar*cloud_num*cloud_num) &
              * (zcnue+2.0_wp)*(zcnue+4.0_wp)/(zcnue+1.0_wp)**2
    ccsrim = 0.25_wp*pi*zecs*v0snow*gamma_fct(zv1s+3.0_wp)
    ccsagg = 0.25_wp*pi*v0snow*gamma_fct(zv1s+3.0_wp)
    ccsdep = 0.26_wp*gamma_fct((zv1s+5.0_wp)/2.0_wp)*SQRT(0.5_wp/zeta)
    ccsvxp = -(zv1s/(zbms+1.0_wp)+1.0_wp)
    ccsvel = zams*v0snow*gamma_fct(zbms+zv1s+1.0_wp)&
      & *(zams*gamma_fct(zbms+1.0_wp))**ccsvxp
    ccsvxp = ccsvxp + 1.0_wp
    ccslam = zams*gamma_fct(zbms+1.0_wp)
    ccslxp = 1.0_wp / (zbms+1.0_wp)
    ccswxp = zv1s*ccslxp
    ccsaxp = -(zv1s+3.0_wp)

!US old    ccsdxp = -(zbms+1.0_wp)/2.0_wp
    ccsdxp = -(zv1s+1.0_wp)/2.0_wp

    ccshi1 = lh_s*lh_s/(zlheat*r_v)
    ccdvtp = 2.11E-5_wp * t0_melt**(-1.94_wp) * 101325.0_wp
    ccidep = 4.0_wp * zami**(-x1o3)
    ! fall speed of rain is approximated as v = 130 D**(1/2), see above
    ! empirical relation adapted from Ulbrich (1983)
    zn0r   = 8.0E6_wp * EXP(3.2_wp*mu_rain) * (0.01_wp)**(-mu_rain)
    ! to tune the zn0r variable
    zn0r   = zn0r * rain_n0_factor

    zar    = pi*zrhow/6.0_wp * zn0r * gamma_fct(mu_rain+4.0_wp) ! pre-factor in lambda
    zcevxp = (mu_rain+2.0_wp)/(mu_rain+4.0_wp)
    zcev   = 2.0_wp*pi*zdv/zhw*zn0r*zar**(-zcevxp) * gamma_fct(mu_rain+2.0_wp)
    zbevxp = (2.0_wp*mu_rain+5.5_wp)/(2.0_wp*mu_rain+8.0_wp)-zcevxp
    zbev   = 0.26_wp * SQRT(0.5_wp*zrho0*130.0_wp/zeta)*zar**(-zbevxp) &
      &    * gamma_fct((2.0_wp*mu_rain+5.5_wp)/2.0_wp) / gamma_fct(mu_rain+2.0_wp)
    zvzxp  = 0.5_wp/(mu_rain+4.0_wp)
    zvz0r  = 130.0_wp*gamma_fct(mu_rain+4.5_wp)/gamma_fct(mu_rain+4.0_wp)*zar**(-zvzxp)

    IF (izdebug > 10) THEN
      WRITE (*,*) 'SRC_GSCP: Initialized hydci_pp_gr'
      WRITE (*,'(A,E10.3)') '      ccslam = ',ccslam
      WRITE (*,'(A,E10.3)') '      ccsvel = ',ccsvel
      WRITE (*,'(A,E10.3)') '      ccsrim = ',ccsrim
      WRITE (*,'(A,E10.3)') '      ccsagg = ',ccsagg
      WRITE (*,'(A,E10.3)') '      ccsdep = ',ccsdep
      WRITE (*,'(A,E10.3)') '      ccslxp = ',ccslxp
      WRITE (*,'(A,E10.3)') '      ccidep = ',ccidep
      WRITE (*,'(A,E10.3)') '      mu_r   = ',mu_rain
      WRITE (*,'(A,E10.3)') '      zn0r   = ',zn0r
      WRITE (*,'(A,E10.3)') '      zbev   = ',zbev
      WRITE (*,'(A,E10.3)') '      zbevxp = ',zbevxp
      WRITE (*,'(A,E10.3)') '      zcev   = ',zcev
      WRITE (*,'(A,E10.3)') '      zcevxp = ',zcevxp
      WRITE (*,'(A,E10.3)') '      zvzxp  = ',zvzxp
      WRITE (*,'(A,E10.3)') '      zvz0r  = ',zvz0r
    ENDIF
  ENDIF

! Some constant coefficients
  znimax = fxna(zthn) ! Maximum number of cloud ice crystals
  zpvsw0 = fpvsw(t0_melt)  ! sat. vap. pressure for t = t0_melt

! select timelevel and timestep for calculations
nx    = nnew
IF ( l2tls ) THEN
  zdt   = dt
ELSE
  zdt   = dt2
ENDIF
zdtr  = 1.0_wp / zdt


!CDIR BEGIN COLLAPSE
zpkr (:,:) = 0.0_wp
zpks (:,:) = 0.0_wp
zpkg (:,:) = 0.0_wp
zprvr(:,:) = 0.0_wp
zprvs(:,:) = 0.0_wp
zprvg(:,:) = 0.0_wp
zvzr (:,:) = 0.0_wp
zvzs (:,:) = 0.0_wp
zvzg (:,:) = 0.0_wp
!CDIR END

! retrieve the required microphysics tracers (at corresponding timelevel)
CALL trcr_get(izerror, idt_qv, ptr_tlev = nx, ptr = qv)
IF (izerror /= 0) THEN
  yzerrmsg = trcr_errorstr(izerror)
  CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
ENDIF
CALL trcr_get(izerror, idt_qc, ptr_tlev = nx, ptr = qc)
IF (izerror /= 0) THEN
  yzerrmsg = trcr_errorstr(izerror)
  CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
ENDIF
CALL trcr_get(izerror, idt_qi, ptr_tlev = nx, ptr = qi)
IF (izerror /= 0) THEN
  yzerrmsg = trcr_errorstr(izerror)
  CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
ENDIF
CALL trcr_get(izerror, idt_qg, ptr_tlev = nx, ptr = qg)
IF (izerror /= 0) THEN
  yzerrmsg = trcr_errorstr(izerror)
  CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
ENDIF
CALL trcr_get(izerror, idt_qr, ptr_tlev = nx, ptr = qr)
IF (izerror /= 0) THEN
  yzerrmsg = trcr_errorstr(izerror)
  CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
ENDIF
CALL trcr_get(izerror, idt_qs, ptr_tlev = nx, ptr = qs)
IF (izerror /= 0) THEN
  yzerrmsg = trcr_errorstr(izerror)
  CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
ENDIF

#ifdef NUDGING
! add part of latent heating calculated in subroutine hydci_pp_gr to model latent
! heating field: subtract temperature from model latent heating field
IF ((llhn .OR. llhnverif)) THEN
  IF (lhn_qrs) THEN
!CDIR COLLAPSE
    qrsflux(:,:,:) = 0.0_wp
  ENDIF
  CALL get_gs_lheating ('add',1,ke)
ENDIF
#endif

! *********************************************************************
! Loop from the top of the model domain to the surface to calculate the
! transfer rates and sedimentation terms    
! *********************************************************************

loop_over_levels: DO  k = 1, ke

  IF ( ldiabf_lh ) THEN
    ! initialize temperature increment due to latent heat
    tinc_lh(:,:,k) = tinc_lh(:,:,k) - t(:,:,k,nx)
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 2: Implicit sedimentation
  !----------------------------------------------------------------------------

!CDIR BEGIN COLLAPSE
  zcrim (:,:) = 0.0_wp
  zcagg (:,:) = 0.0_wp
  zbsdep(:,:) = 0.0_wp
  zvz0s (:,:) = 0.0_wp
  zn0s  (:,:) = zn0s0
!CDIR END

  !----------------------------------------------------------------------------
  ! 2.1: Preparations for computations and to check the different conditions
  !----------------------------------------------------------------------------

  ! Nullify counters
  ic1 = 0
  ic2 = 0
  ic3 = 0

!CDIR ON_ADB(idx1)
!CDIR ON_ADB(jdx1)
!CDIR ON_ADB(idx2)
!CDIR ON_ADB(jdx2)
!CDIR ON_ADB(idx3)
!CDIR ON_ADB(jdx3)
  DO j = jstartpar, jendpar
    DO i = istartpar, iendpar

      qrg = qr(i,j,k)
      qsg = qs(i,j,k)
      qgg = qg(i,j,k) 
      qvg = qv(i,j,k)
      qcg = qc(i,j,k)
      qig = qi(i,j,k)
      tg  = t(i,j,k,nx)
      ppg = p0(i,j,k) + pp(i,j,k,nx)
      rhog = rho(i,j,k)

      !..for density correction of fall speeds
      z1orhog(i,j) = 1.0_wp/rhog
      zrho1o2(i,j) = EXP(LOG(zrho0*z1orhog(i,j))*x1o2)

      zpres(i,j)= ppg
      zqrk(i,j) = qrg * rhog
      zqsk(i,j) = qsg * rhog
      zqgk(i,j) = qgg * rhog

      llqr(i,j)  = (zqrk(i,j) > zqmin)
      llqs       = (zqsk(i,j) > zqmin)
      llqg       = (zqgk(i,j) > zqmin)

      zdh        = hhl(i,j,k)-hhl(i,j,k+1)
      zdtdh(i,j) = 0.5_wp * zdt / zdh

!NEC_CB Never used, only in the three lines below
!NEC_CB      zpkm1r(i,j) = zpkr(i,j)
!NEC_CB      zpkm1s(i,j) = zpks(i,j)
!NEC_CB      zpkm1g(i,j) = zpkg(i,j)
      zzar  (i,j) = zqrk(i,j)/zdtdh(i,j) + zprvr(i,j) + zpkr(i,j)
      zzas  (i,j) = zqsk(i,j)/zdtdh(i,j) + zprvs(i,j) + zpks(i,j)
      zzag  (i,j) = zqgk(i,j)/zdtdh(i,j) + zprvg(i,j) + zpkg(i,j)

      zpkr(i,j)   = MIN( zpkr(i,j) , zzar(i,j) )
      ! Check the different conditions
      IF (llqs) THEN
        ic1 = ic1 + 1
        idx1(ic1) = i
        jdx1(ic1) = j
      ENDIF
      IF (llqr(i,j)) THEN
        ic2 = ic2 + 1
        idx2(ic2) = i
        jdx2(ic2) = j
      ENDIF
      IF (llqg) THEN
        ic3 = ic3 + 1
        idx3(ic3) = i
        jdx3(ic3) = j
      ENDIF

    ENDDO
  ENDDO

  !----------------------------------------------------------------------------
  ! 2.2: Calculate snow intercept and derived parameters only for chosen points
  !----------------------------------------------------------------------------

!CDIR BEGIN COLLAPSE
  zpkr  (:,:) = 0.0_wp
  zpks  (:,:) = 0.0_wp
  zpkg  (:,:) = 0.0_wp
!CDIR END

  zln1o2=EXP (ccswxp * LOG (0.5_wp))
!CDIR NODEP,VOVERTAKE,VOB
!CDIR ON_ADB(idx1)
!CDIR ON_ADB(jdx1)
  ic1_loop1: DO i1d = 1,ic1
    i = idx1(i1d)
    j = jdx1(i1d)

    qsg       = qs(i,j,k)
    tg        = t(i,j,k,nx)

    IF (isnow_n0temp == 1) THEN
      ! Calculate n0s using the temperature-dependent
      ! formula of Field et al. (2005)
      ztc = tg - t0_melt
      ztc = MAX(MIN(ztc,0.0_wp),-40.0_wp)
      zn0s(i,j) = zn0s1*EXP(zn0s2*ztc)
      zn0s(i,j) = MIN(zn0s(i,j),1e9_wp)
      zn0s(i,j) = MAX(zn0s(i,j),1e6_wp)
    ELSEIF (isnow_n0temp == 2) THEN
      ! Calculate n0s using the temperature-dependent moment
      ! relations of Field et al. (2005)
      ztc = tg - t0_melt
      ztc = MAX(MIN(ztc,0.0_wp),-40.0_wp)

      nn  = 3
      hlp = mma(1)      +mma(2)*ztc      +mma(3)*nn       +mma(4)*ztc*nn+mma(5)*ztc**2 &
          + mma(6)*nn**2+mma(7)*ztc**2*nn+mma(8)*ztc*nn**2+mma(9)*ztc**3+mma(10)*nn**3
      alf = 10.0_wp**hlp
      bet = mmb(1)      +mmb(2)*ztc      +mmb(3)*nn       +mmb(4)*ztc*nn+mmb(5)*ztc**2 &
          + mmb(6)*nn**2+mmb(7)*ztc**2*nn+mmb(8)*ztc*nn**2+mmb(9)*ztc**3+mmb(10)*nn**3

      ! Here is the exponent bms=2.0 hardwired! not ideal! (Uli Blahak)
      m2s = qsg * rho(i,j,k) / zams  ! UB rho added as bugfix
      m3s = alf*EXP(bet*LOG(m2s))

      hlp  = zn0s1*EXP(zn0s2*ztc)
      ! UB: the 13.5 is actually 3^3 / gamma(3) ...
      zn0s(i,j) = 13.50_wp * m2s*(m2s/m3s)**3
      zn0s(i,j) = MAX(zn0s(i,j),0.5_wp*hlp) 
      zn0s(i,j) = MIN(zn0s(i,j),1e2_wp*hlp)
      zn0s(i,j) = MIN(zn0s(i,j),1e9_wp)
      zn0s(i,j) = MAX(zn0s(i,j),1e6_wp)
    ELSE
      ! Old constant n0s
      zn0s(i,j) = 8.0e5_wp
    ENDIF
    zcrim (i,j) = ccsrim*zn0s(i,j)
    zcagg (i,j) = ccsagg*zn0s(i,j)
    zbsdep(i,j) = ccsdep*SQRT(v0snow)
    zvz0s (i,j) = ccsvel*EXP(ccsvxp * LOG(zn0s(i,j)))
    zlnqsk = zvz0s(i,j) * EXP (ccswxp * LOG (zqsk(i,j))) * zrho1o2(i,j)
    zpks  (i,j) = zqsk(i,j) * zlnqsk
    IF (zvzs(i,j) == 0.0_wp) THEN
      zvzs(i,j) = zlnqsk * zln1o2
    ENDIF
  ENDDO ic1_loop1

  ! sedimentation fluxes
  zln1o2=EXP (zvzxp * LOG (0.5_wp))
!CDIR NODEP,VOVERTAKE,VOB
!CDIR ON_ADB(idx2)
!CDIR ON_ADB(jdx2)
  ! llqr:
  ic2_loop1: DO i1d = 1,ic2
    i = idx2(i1d)
    j = jdx2(i1d)
    zlnqrk = zvz0r * EXP (zvzxp * LOG (zqrk(i,j))) * zrho1o2(i,j)
    zpkr(i,j) = zqrk(i,j) * zlnqrk
    IF (zvzr(i,j) == 0.0_wp) THEN
      zvzr(i,j) = zlnqrk * zln1o2
    ENDIF
  ENDDO ic2_loop1

  ! the other computation from here has been added above:  llqs:

  zln1o2=EXP (zexpsedg * LOG (0.5_wp))
!CDIR NODEP,VOVERTAKE,VOB
!CDIR ON_ADB(idx3)
!CDIR ON_ADB(jdx3)
  ! llqg:
  ic3_loop1: DO i1d = 1,ic3
    i = idx3(i1d)
    j = jdx3(i1d)
    zlnqgk = zvz0g * EXP (zexpsedg * LOG (zqgk(i,j))) * zrho1o2(i,j)
    zpkg(i,j) = zqgk(i,j) * zlnqgk
    IF (zvzg(i,j) == 0.0_wp) THEN
      zvzg(i,j) = zlnqgk * zln1o2
    ENDIF
  ENDDO ic3_loop1

  !----------------------------------------------------------------------------
  ! 2.3: Second part of preparations
  !----------------------------------------------------------------------------

!CDIR BEGIN COLLAPSE
  zeln7o8qrk   (:,:) = 0.0_wp
  zeln27o16qrk (:,:) = 0.0_wp
  zeln13o8qrk  (:,:) = 0.0_wp
  zsrmax       (:,:) = 0.0_wp
  zeln3o4qsk   (:,:) = 0.0_wp
  zeln8qsk     (:,:) = 0.0_wp
  zssmax       (:,:) = 0.0_wp
  zelnrimexp_g (:,:) = 0.0_wp
  zeln6qgk     (:,:) = 0.0_wp
  zsgmax       (:,:) = 0.0_wp
  zcsdep       (:,:) = 3.2E-2_wp
  zcidep       (:,:) = 1.3E-5_wp
  zcslam       (:,:) = 1e10_wp
!CDIR END

  ! Nullify counters again
  ic1 = 0
  ic2 = 0
  ic3 = 0
  ic4 = 0
  ic5 = 0
  ic6 = 0
  ic7 = 0
  ic8 = 0

  DO j = jstartpar, jendpar
    DO i = istartpar, iendpar

      qcg       = qc(i,j,k)
      qvg       = qv(i,j,k)
      qig       = qi(i,j,k)
      qrg       = qr(i,j,k)
      qsg       = qs(i,j,k)
      qgg       = qg(i,j,k)
      tg        =  t(i,j,k,nx)
      ppg       = zpres(i,j) !NEC_CB p0(i,j,k) + pp(i,j,k,nx)
      rhog      = rho(i,j,k)
!     zqsk(i,j) = qsg * rhog
!     zqgk(i,j) = qgg * rhog
!NEC_CB      llqs      = (zqsk(i,j) > zqmin)
!NEC_CB      llqr(i,j) = (zqrk(i,j) > zqmin)
!NEC_CB      llqg      = (zqgk(i,j) > zqmin)

      zpkr(i,j)   = MIN( zpkr(i,j) , zzar(i,j) )
      zpks(i,j)   = MIN( zpks(i,j) , zzas(i,j) )
      zpkg(i,j)   = MIN( zpkg(i,j) , zzag(i,j) )
      zzar(i,j)   = zdtdh(i,j) * (zzar(i,j)-zpkr(i,j))
      zzas(i,j)   = zdtdh(i,j) * (zzas(i,j)-zpks(i,j))
      zzag(i,j)   = zdtdh(i,j) * (zzag(i,j)-zpkg(i,j))

!NEC_CB Moved to section 2.2 to reuse log(zq?k(i,j))
!NEC_CB      IF (zvzr(i,j) == 0.0_wp .AND. llqr(i,j)) THEN
!NEC_CB        zvzr(i,j) = zvz0r * EXP (zvzxp * LOG (0.5_wp*zqrk(i,j)))
!NEC_CB      ENDIF
!NEC_CB      IF (zvzs(i,j) == 0.0_wp .AND. llqs) THEN
!NEC_CB        zvzs(i,j) = zvz0s(i,j) * EXP (ccswxp * LOG (0.5_wp*zqsk(i,j)))
!NEC_CB      ENDIF
!NEC_CB      IF (zvzg(i,j) == 0.0_wp .AND. llqg) THEN
!NEC_CB        zvzg(i,j) = zvz0g * EXP (zexpsedg * LOG (0.5_wp*zqgk(i,j)))
!NEC_CB      ENDIF
      zimr(i,j)   = 1.0_wp / (1.0_wp + zvzr(i,j) * zdtdh(i,j))
      zims(i,j)   = 1.0_wp / (1.0_wp + zvzs(i,j) * zdtdh(i,j))
      zimg(i,j)   = 1.0_wp / (1.0_wp + zvzg(i,j) * zdtdh(i,j))

      zqrk(i,j)   = zzar(i,j) * zimr(i,j)
      zqsk(i,j)   = zzas(i,j) * zims(i,j)
      zqgk(i,j)   = zzag(i,j) * zimg(i,j)

      llqr(i,j) = (zqrk(i,j) > zqmin)
      llqs = (zqsk(i,j) > zqmin)
      llqg = (zqgk(i,j) > zqmin)
      llqc = (      qcg > zqmin)
      llqi = (      qig > zqmin)
      
      ! Initialize index arrays for subsequent IF-clauses
      IF (llqr(i,j)) THEN
       ic1 = ic1 + 1
       idx1(ic1) = i
       jdx1(ic1) = j
      ENDIF
      IF (llqs) THEN
       ic2 = ic2 + 1
       idx2(ic2) = i
       jdx2(ic2) = j
      ENDIF
      IF (llqg) THEN
       ic3 = ic3 + 1
       idx3(ic3) = i
       jdx3(ic3) = j
      ENDIF
      IF (llqi .OR. llqs) THEN
       ic4 = ic4 + 1
       idx4(ic4) = i
       jdx4(ic4) = j
      ENDIF
      IF ( (tg < 248.15_wp) .AND. (qvg > 8.E-6_wp) .AND. (.NOT.llqi) ) THEN
       ic5 = ic5 + 1
       idx5(ic5) = i
       jdx5(ic5) = j
      ENDIF
      IF ( llqc ) THEN
       ic6 = ic6 + 1
       idx6(ic6) = i
       jdx6(ic6) = j
      ENDIF
      IF (llqi .OR. llqs .OR. llqg) THEN
       ic7 = ic7 + 1
       idx7(ic7) = i
       jdx7(ic7) = j
      ENDIF
      IF (llqr(i,j) .AND. (qcg <= 0.0_wp))  THEN
       ic8 = ic8 + 1
       idx8(ic8) = i
       jdx8(ic8) = j
      ENDIF
    ENDDO
  ENDDO

  !----------------------------------------------------------------------------
  ! 2.4: IF (llqr): ic1
  !----------------------------------------------------------------------------

!CDIR NODEP,VOVERTAKE,VOB
  ic1_loop3: DO i1d = 1,ic1
    i = idx1(i1d)
    j = jdx1(i1d)

    rhog = rho(i,j,k)
    qcg  =  qc(i,j,k)
    qig  =  qi(i,j,k)
    tg   =   t(i,j,k,nx)
    llqi = ( qig > zqmin)

    zlnqrk         = LOG (zqrk(i,j))
    zsrmax(i,j)    = zzar(i,j)/rhog*zdtr
    IF (qig+qcg > zqmin ) THEN
      zeln7o8qrk  (i,j) = EXP (x7o8   * zlnqrk)
    ENDIF
    IF ( tg < ztrfrz ) THEN
      zeln27o16qrk(i,j) = EXP (x27o16 * zlnqrk)
    ENDIF
    IF (llqi) THEN
      zeln13o8qrk (i,j) = EXP (x13o8  * zlnqrk)
    ENDIF
    ! IF (qcg .LE. 0.0_wp ) THEN    ! zeln3o16qrk is not used below !!
    !   zeln3o16qrk (i,j) = EXP (x3o16  * zlnqrk)
    ! ENDIF
  ENDDO ic1_loop3

  !----------------------------------------------------------------------------
  ! 2.5: IF (llqs): ic2
  !----------------------------------------------------------------------------

!CDIR NODEP,VOVERTAKE,VOB
  ic2_loop2: DO i1d=1,ic2
    i = idx2(i1d)
    j = jdx2(i1d)

    rhog = rho(i,j,k)
    qcg  = qc(i,j,k)
    qig  = qi(i,j,k)
    zlnqsk        = LOG (zqsk(i,j))
    zssmax(i,j)   = zzas(i,j) / rhog*zdtr
    IF (qig+qcg > zqmin) THEN
      zeln3o4qsk(i,j) = EXP (x3o4 *zlnqsk)
    ENDIF
    zeln8qsk(i,j) = EXP (0.8_wp *zlnqsk)
  ENDDO ic2_loop2

  !----------------------------------------------------------------------------
  ! 2.6: IF (llqg): ic3
  !----------------------------------------------------------------------------

!CDIR NODEP,VOVERTAKE,VOB
  ic3_loop2: DO i1d=1,ic3
    i = idx3(i1d) 
    j = jdx3(i1d)

    rhog = rho(i,j,k)
    qcg  = qc(i,j,k)
    qig  = qi(i,j,k)

    zlnqgk        = LOG (zqgk(i,j))
    zsgmax(i,j)   = zzag(i,j) / rhog*zdtr
    IF (qig+qcg > zqmin) THEN
      zelnrimexp_g(i,j) = EXP (zrimexp_g * zlnqgk)
    ENDIF
    zeln6qgk(i,j) = EXP (0.6_wp *zlnqgk)
  ENDDO ic3_loop2

  !----------------------------------------------------------------------------
  ! 2.7: slope of snow PSD and coefficients for depositional growth (lqi,lqs)
  !----------------------------------------------------------------------------

!CDIR NODEP,VOVERTAKE,VOB
  ic4_loop1: DO i1d=1,ic4
    i = idx4(i1d)
    j = jdx4(i1d)

    rhog = rho(i,j,k)
    qcg  = qc(i,j,k)
    qig  = qi(i,j,k)
    ppg  = zpres(i,j) !NEC_CB p0(i,j,k) + pp(i,j,k,nx)
    tg   = t(i,j,k,nx)
    llqs = (zqsk(i,j) > zqmin)

    zqvsi       = fqvs( fpvsi(tg), ppg )
    zdvtp       = ccdvtp * EXP(1.94_wp * LOG(tg)) / ppg
    zhi         = ccshi1*zdvtp*rhog*zqvsi/(tg*tg)
    hlp         = zdvtp / (1.0_wp + zhi)
    zcidep(i,j) = ccidep * hlp

    IF (llqs) THEN
      zcslam(i,j) = EXP(ccslxp * LOG(ccslam * zn0s(i,j) / zqsk(i,j) ))
      zcslam(i,j) = MIN(zcslam(i,j),1e15_wp)
      zcsdep(i,j) = 4.0_wp * zn0s(i,j) * hlp
    ENDIF

  ENDDO ic4_loop1

  !----------------------------------------------------------------------------
  ! Section 3: Initialize the conversion rates with zeros in every layer
  !            Deposition nucleation for low temperatures
  !----------------------------------------------------------------------------

!CDIR BEGIN COLLAPSE
  scau   (:,:) = 0.0_wp
  scac   (:,:) = 0.0_wp
  snuc   (:,:) = 0.0_wp
  scfrz  (:,:) = 0.0_wp
  simelt (:,:) = 0.0_wp
  sidep  (:,:) = 0.0_wp
  ssdep  (:,:) = 0.0_wp
  sgdep  (:,:) = 0.0_wp
  sdau   (:,:) = 0.0_wp
  srim   (:,:) = 0.0_wp
  srim2  (:,:) = 0.0_wp
  sconsg (:,:) = 0.0_wp
  sshed  (:,:) = 0.0_wp
  sicri  (:,:) = 0.0_wp
  srcri  (:,:) = 0.0_wp
  sagg   (:,:) = 0.0_wp
  sagg2  (:,:) = 0.0_wp
  siau   (:,:) = 0.0_wp
  ssmelt (:,:) = 0.0_wp
  sgmelt (:,:) = 0.0_wp
  sev    (:,:) = 0.0_wp
  sconr  (:,:) = 0.0_wp
  srfrz  (:,:) = 0.0_wp
!CDIR END

  ! Deposition nucleation for low temperatures below a threshold
!CDIR NODEP,VOVERTAKE,VOB
  ic5_loop1: DO i1d = 1,ic5
    i = idx5(i1d)
    j = jdx5(i1d)

    rhog = rho(i,j,k)
    qvg  = qv(i,j,k)
    ppg  = zpres(i,j) !NEC_CB p0(i,j,k) + pp(i,j,k,nx)
    tg   = t(i,j,k,nx)

    zqvsi     = fqvs( fpvsi(tg), ppg )
    IF( qvg > zqvsi ) THEN
      znin      = MIN( fxna(tg), znimax )
      snuc(i,j) = zmi0 / rhog * znin * zdtr
    ENDIF
  ENDDO ic5_loop1

  !----------------------------------------------------------------------------
  ! Section 4: Search for cloudy grid points with cloud water and
  !  (ic6)     calculation of the conversion rates involving qc, qs and qg
  !----------------------------------------------------------------------------

!CDIR NODEP,VOVERTAKE,VOB
  ic6_loop1: DO i1d = 1,ic6
    i = idx6(i1d)
    j = jdx6(i1d)

    rhog = rho(i,j,k)
    qcg  =  qc(i,j,k)
    qrg  =  qr(i,j,k)
    qig  =  qi(i,j,k)
    ppg  =  zpres(i,j) !NEC_CB p0(i,j,k) + pp(i,j,k,nx)
    tg   =   t(i,j,k,nx)
    llqs = (zqsk(i,j) > zqmin)
    llqi = (      qig > zqmin)

    zscmax = qcg*zdtr
    IF( tg > zthn ) THEN
      IF (iautocon == 0) THEN
        ! Kessler (1969) autoconversion rate
        scau(i,j) = zccau * MAX( qcg - qc0, 0.0_wp )
        scac(i,j) = zcac  * qcg * zeln7o8qrk(i,j)
      ELSEIF (iautocon == 1) THEN
        ! Seifert and Beheng (2001) autoconversion rate
        ! with constant cloud droplet number concentration cloud_num
        IF (qcg > 1e-6) THEN
          ztau  = MIN(1.0_wp-qcg/(qcg+qrg),0.9_wp)
          zphi  = zkphi1 * ztau**zkphi2 * (1.0_wp - ztau**zkphi2)**3
          scau(i,j) = zconst * qcg*qcg*qcg*qcg &
                      * (1.0_wp + zphi/(1.0_wp - ztau)**2)
          zphi      = (ztau/(ztau+zkphi3))**4
          scac(i,j) = zkcac * qcg * qrg * zphi
        ENDIF
      ENDIF
      IF (llqr(i,j)) THEN
        ! Calculation of in-cloud rainwater freezing
        IF ( tg < ztrfrz ) THEN
          ztfrzdiff=ztrfrz-tg
          srfrz(i,j) = zcrfrz*ztfrzdiff*SQRT(ztfrzdiff)* zeln27o16qrk(i,j)
        ENDIF
      ENDIF
      IF (llqs) THEN
        srim(i,j) = zcrim(i,j) * qcg *  EXP(ccsaxp * LOG(zcslam(i,j)))
      ENDIF
      srim2(i,j)  = zcrim_g * qcg * zelnrimexp_g(i,j)
      IF( tg >= t0_melt ) THEN
        sshed(i,j) = srim(i,j)+srim2(i,j)
        srim (i,j) = 0.0_wp
        srim2(i,j) = 0.0_wp
      ELSE
        IF (qcg.GE.qc0) THEN
          sconsg(i,j) = zcsg * qcg * zeln3o4qsk(i,j)
        ENDIF
      ENDIF
      ! Check for maximum depletion of cloud water and adjust the
      ! transfer rates accordingly
      zscsum = scau(i,j) + scac(i,j) + srim(i,j) + srim2(i,j) + sshed(i,j)
      zcorr  = zscmax / MAX( zscmax, zscsum )
      scau   (i,j) = zcorr*scau(i,j)
      scac   (i,j) = zcorr*scac(i,j)
      srim   (i,j) = zcorr*srim(i,j)
      srim2  (i,j) = zcorr*srim2(i,j)
      sshed  (i,j) = zcorr*sshed(i,j)
      sconsg (i,j) = MIN (sconsg(i,j), srim(i,j)+zssmax(i,j))
    ELSE !tg >= tg: ! hom. freezing of cloud and rain water
      scfrz(i,j) = zscmax
      srfrz(i,j) = zsrmax(i,j)
    ENDIF
    ! Calculation of heterogeneous nucleation of cloud ice.
    ! This is done in this section, because we require water saturation
    ! for this process (i.e. the existence of cloud water) to exist.
    ! Heterogeneous nucleation is assumed to occur only when no
    ! cloud ice is present and the temperature is below a nucleation
    ! threshold.
    IF( tg <= 267.15_wp .AND. .NOT.llqi ) THEN
      znin      = MIN( fxna(tg), znimax )
      snuc(i,j) = zmi0 / rhog * znin * zdtr
    ENDIF
  ENDDO ic6_loop1

!CDIR NODEP,VOVERTAKE,VOB
  ic7_loop1: DO i1d = 1,ic7
    i = idx7(i1d)
    j = jdx7(i1d)

    rhog = rho(i,j,k)
    qcg  =  qc(i,j,k)
    qrg  =  qr(i,j,k)
    qig  =  qi(i,j,k)
    qvg  =  qv(i,j,k)
    ppg  =  zpres(i,j) !NEC_CB p0(i,j,k) + pp(i,j,k,nx)
    tg   =   t(i,j,k,nx)
    llqs = (zqsk(i,j) > zqmin)
    llqi = (      qig > zqmin)

    IF (tg<=t0_melt) THEN           ! cold case

      !------------------------------------------------------------------------
      ! Section 5: Search for cold grid points with cloud ice and/or
      !            snow/graupel and calculation of the conversion rates
      !            involving qi, qs and qg
      !------------------------------------------------------------------------

      zqvsi     = fqvs( fpvsi(tg), ppg )
      zqvsidiff = qvg-zqvsi
      zsvmax    = zqvsidiff * zdtr
      IF (llqi) THEN
        zeff       = MAX(0.2_wp,MIN(EXP(0.09_wp*(tg-t0_melt)),1.0_wp))
        sagg(i,j)  = zeff * qig * zcagg(i,j) * EXP(ccsaxp*LOG(zcslam(i,j)))
        sagg2(i,j) = zeff * qig * zcagg_g * zelnrimexp_g(i,j)
        siau(i,j)  = zeff * zciau * MAX( qig - qi0, 0.0_wp )
        znin       = MIN( fxna(tg), znimax )
        zmi        = MIN( rhog*qig/znin, zmimax )
        zmi        = MAX( zmi0, zmi )
        znid       = rhog * qig/zmi
        zlnlogmi   = LOG (zmi)
        sidep(i,j) = zcidep(i,j) * znid * EXP(0.33_wp * zlnlogmi) * zqvsidiff
        zsvidep    = 0.0_wp
        zsvisub    = 0.0_wp
        zsimax     = qig*zdtr
        IF( sidep(i,j) > 0.0_wp ) THEN
          zsvidep = MIN( sidep(i,j), zsvmax )
        ELSEIF ( sidep(i,j) < 0.0_wp ) THEN
          zsvisub  =   MAX ( sidep(i,j),  zsvmax)
          zsvisub  = - MAX (    zsvisub, -zsimax)
        ENDIF
        zlnlogmi   = LOG  (zmsmin/zmi)
        zztau      = 1.5_wp*( EXP(0.66_wp*zlnlogmi) - 1.0_wp)
        sdau(i,j)  = zsvidep/zztau
        sicri(i,j) = zcicri * qig * zeln7o8qrk(i,j)
        srcri(i,j) = zcrcri * (qig/zmi) * zeln13o8qrk(i,j)
      ELSE
        zsimax    =  0.0_wp
        zsvidep   =  0.0_wp
        zsvisub   =  0.0_wp
      ENDIF

      zxfac      = 1.0_wp + zbsdep(i,j) * EXP(ccsdxp*LOG(zcslam(i,j)))
      ssdep(i,j) = zcsdep(i,j) * zxfac * zqvsidiff / (zcslam(i,j)+eps_div)**2
      !ssdep     = (2.91955_wp-0.0109928_wp*tg                   &
      !             + 15871.3_wp/ppg + 1.74744E-6_wp*ppg) *      &
      !             zqvsidiff * zeln8qsk
      sgdep(i,j) = (0.398561_wp-0.00152398_wp*tg                 &
                    + 2554.99_wp/ppg+ 2.6531E-7_wp*ppg) *        &
                     zqvsidiff * zeln6qgk(i,j)
      ! Check for maximal depletion of cloud ice
      ! No check is done for depositional autoconversion because
      ! this is a always a fraction of the gain rate due to
      ! deposition (i.e the sum of this rates is always positive)
      zsisum = siau(i,j) + sagg(i,j) + sagg2(i,j) + sicri(i,j) + zsvisub
      zcorr  = 0.0_wp
      IF( zsimax > 0.0_wp ) zcorr  = zsimax / MAX( zsimax, zsisum )
      sidep(i,j)  = zsvidep - zcorr*zsvisub
      siau (i,j)  = zcorr*siau(i,j)
      sagg (i,j)  = zcorr*sagg(i,j)
      sagg2(i,j)  = zcorr*sagg2(i,j)
      sicri(i,j)  = zcorr*sicri(i,j)
      IF ( zqvsidiff < 0.0_wp ) THEN
        ssdep(i,j) = MAX(ssdep(i,j), - zssmax(i,j))
        sgdep(i,j) = MAX(sgdep(i,j), - zsgmax(i,j))
      ENDIF
    ELSE          ! tg - warm case:

      !------------------------------------------------------------------------
      ! Section 6: Search for warm grid points with cloud ice and/or 
      !            snow/graupel and calculation of the melting rates
      !            of qi, qs and qg
      !------------------------------------------------------------------------

      ! cloud ice melts instantaneously
      simelt(i,j) = qig*zdtr

      zqvsw0     = fqvs( zpvsw0, ppg)
      zqvsw0diff = qvg-zqvsw0

      IF ( tg > (t0_melt-ztcrit*zqvsw0diff) ) THEN
        !calculate melting rate
        zx1         = (tg - t0_melt) + zasmel*zqvsw0diff
        ssmelt(i,j) = (79.6863_wp/ppg+0.612654E-3_wp)* zx1 * zeln8qsk(i,j)
        ssmelt(i,j) = MIN (ssmelt(i,j),zssmax(i,j))
        sgmelt(i,j) = (12.31698_wp/ppg+7.39441e-05_wp)* zx1 * zeln6qgk(i,j)
        sgmelt(i,j) = MIN (sgmelt(i,j), zsgmax(i,j))
        !deposition + melting, ice particle temperature: t0_melt
        !calculation without howell-factor!
        ssdep(i,j)  = (31282.3_wp/ppg+0.241897_wp)       &
                      * zqvsw0diff * zeln8qsk(i,j)
        sgdep(i,j)  = (0.153907_wp-ppg*7.86703e-07_wp)  &
                      * zqvsw0diff * zeln6qgk(i,j)
        IF (zqvsw0diff < 0.0_wp) THEN
          !melting + evaporation of snow/graupel
          ssdep(i,j) = MAX (-zssmax(i,j),ssdep(i,j))
          sgdep(i,j) = MAX (-zsgmax(i,j),sgdep(i,j))
          !melt water evaporates
          ssmelt(i,j) = ssmelt(i,j)+ssdep(i,j)
          sgmelt(i,j) = sgmelt(i,j)+sgdep(i,j)
          ssmelt(i,j) = MAX( ssmelt(i,j), 0.0_wp )
          sgmelt(i,j) = MAX( sgmelt(i,j), 0.0_wp )
        ELSE
          !deposition on snow/graupel is interpreted as increase
          !in rain water ( qv --> qr, sconr)
          !therefore,  sconr=(zssdep+zsgdep)
          sconr(i,j)=ssdep(i,j)+sgdep(i,j)
          ssdep(i,j)=0.0_wp
          sgdep(i,j)=0.0_wp
        ENDIF
      ELSE
        !if t<t_crit
        !no melting, only evaporation of snow/graupel
        zqvsw      = fqvs( fpvsw(tg), ppg )
        zqvsidiff  = qvg-zqvsw
        ssdep(i,j) = (0.28003_wp-ppg*0.146293E-6_wp) &
                       * zqvsidiff * zeln8qsk(i,j)
        sgdep(i,j) = (0.0418521_wp-ppg*4.7524E-8_wp) &
                       * zqvsidiff *zeln6qgk(i,j)
        ssdep(i,j) = MAX(-zssmax(i,j) ,ssdep(i,j) )
        sgdep(i,j) = MAX(-zsgmax(i,j) ,sgdep(i,j) )
     ENDIF !t_crit
    ENDIF !tg

  ENDDO ic7_loop1

  !----------------------------------------------------------------------------
  ! Section 7: Search for grid points with rain in subsaturated areas
  !  (ic8)     and calculation of the evaporation rate of rain
  !----------------------------------------------------------------------------

!CDIR NODEP,VOVERTAKE,VOB
  ic8_loop1: DO i1d = 1,ic8
    i = idx8(i1d)
    j = jdx8(i1d)

    qvg = qv(i,j,k)
    tg  =  t(i,j,k,nx)
    ppg = zpres(i,j) !NEC_CB p0(i,j,k) + pp(i,j,k,nx)

    zlnqrk      = LOG (zqrk(i,j))
    zqvsw       = fqvs( fpvsw(tg), ppg )
    zx1         = 1.0_wp + zbev * EXP (zbevxp  * zlnqrk)
    sev(i,j)    = zcev*zx1*(zqvsw - qvg) * EXP (zcevxp  * zlnqrk)

    IF( tg > zthn ) THEN
      ! Calculation of below-cloud rainwater freezing
      IF ( tg < ztrfrz ) THEN
        ztfrzdiff = ztrfrz-tg
        srfrz(i,j)= zcrfrz*ztfrzdiff*SQRT(ztfrzdiff) * zeln27o16qrk(i,j)
      ENDIF
    ELSE ! Hom. freezing of rain water
      srfrz(i,j) = zsrmax(i,j)
    ENDIF
  ENDDO ic8_loop1

  !----------------------------------------------------------------------------
  ! Section 9: Calculate the total tendencies of the prognostic variables.
  !            Update the prognostic variables in the interior domain.
  !----------------------------------------------------------------------------

  ic1=0
  ic2=0
  ic3=0
!CDIR ON_ADB(idx1)
!CDIR ON_ADB(jdx1)
!CDIR ON_ADB(idx2)
!CDIR ON_ADB(jdx2)
!CDIR ON_ADB(idx3)
!CDIR ON_ADB(jdx3)
  DO j = jstartpar, jendpar
    DO i = istartpar, iendpar

      rhog = rho(i,j,k)
      qig  =  qi(i,j,k)
      ppg       = zpres(i,j) !NEC_CB p0(i,j,k) + pp(i,j,k,nx)

      zsrsum = sev(i,j) + srfrz(i,j) + srcri(i,j)
      zcorr  = 1.0_wp
      IF(zsrsum > 0) THEN
        zcorr  = zsrmax(i,j) / MAX( zsrmax(i,j), zsrsum )
      ENDIF
      sev(i,j)   = zcorr*sev(i,j)
      srfrz(i,j) = zcorr*srfrz(i,j)
      srcri(i,j) = zcorr*srcri(i,j)

      zqvt =   sev   (i,j) - sidep (i,j) - ssdep (i,j) - sgdep (i,j)         &
             - snuc  (i,j) - sconr (i,j)
      zqct =   simelt(i,j) - scau  (i,j) - scfrz (i,j) - scac  (i,j)         &
             - sshed (i,j) - srim  (i,j) - srim2 (i,j)
      zqit =   snuc  (i,j) + scfrz (i,j) - simelt(i,j) - sicri (i,j)         &
             + sidep (i,j) - sdau  (i,j) - sagg  (i,j) - sagg2 (i,j) - siau(i,j)
      zqrt =   scau  (i,j) + sshed (i,j) + scac  (i,j) + ssmelt(i,j)         &
             + sgmelt(i,j) - sev   (i,j) - srcri (i,j) - srfrz (i,j) + sconr(i,j)
      zqst =   siau  (i,j) + sdau  (i,j) - ssmelt(i,j) + srim  (i,j)         &
             + ssdep (i,j) + sagg  (i,j) - sconsg(i,j)
      zqgt =   sagg2 (i,j) - sgmelt(i,j) + sicri (i,j) + srcri (i,j)         &
             + sgdep (i,j) + srfrz (i,j) + srim2 (i,j) + sconsg(i,j)

      ztt = cpdr * ( lh_v*(zqct+zqrt) + lh_s*(zqit+zqst+zqgt) )

      IF (lsppt)THEN
        IF(ntstep>0.AND.i>=istart.and.i<=iend.and.j>=jstart.and.j<=jend) THEN
          CALL apply_tqx_tend_adj(itype_qxpert_rn,itype_qxlim_rn,ppg, &
                t(i,j,k,nx),qv(i,j,k), qc(i,j,k),qi(i,j,k),qr(i,j,k),qs(i,j,k),&
                pertstoph(i,j,k),ztt,zqvt,zqct,zqit,zqrt,zqst,ldum,qg(i,j,k),zqgt)
        ENDIF
      ENDIF

      ! Update variables
      qig = MAX ( 0.0_wp, qig + zqit*zdt)
      qrg = MAX ( 0.0_wp, (zzar(i,j)/rhog + zqrt*zdt)*zimr(i,j))
      qsg = MAX ( 0.0_wp, (zzas(i,j)/rhog + zqst*zdt)*zims(i,j))
      qgg = MAX ( 0.0_wp, (zzag(i,j)/rhog + zqgt*zdt)*zimg(i,j))

#ifdef MESSY
        ! convert kg/(kg*s) into kg/(m2*s)
        precmelt_ls(i,j,k)   = (ssmelt(i,j) + simelt(i,j) + sgmelt(i,j)) &
                               * rhog * (hhl(i,j,k)-hhl(i,j,k+1))
        ! rain formation in kg/kg
        rainform_bave(i,j,k) = zqrt * zdt
        ! snow formation in kg/kg
        snowform_bave(i,j,k) = ( zqst + zqgt ) * zdt
        ! rain flux  in kg/(m2*s)
        precflx_ls(i,j,k)  = zqrt  * rhog * (hhl(i,j,k)-hhl(i,j,k+1))
        ! snow flux in  kg/(m2*s)
        snowflx_ls(i,j,k)  = (zqst + zqgt) * rhog * (hhl(i,j,k)-hhl(i,j,k+1))
        ! rain flux without incloud production kg/(m2*s)
        ! in-coming rain + melting snow - evaporation - freezing of rain
        precflxno_ls(i,j,k)    = zprvr(i,j) +  &
                    (ssmelt(i,j) + sgmelt(i,j) - sev(i,j) - srfrz(i,j))  &
                     * rhog * (hhl(i,j,k)-hhl(i,j,k+1))
        ! snow flux without incloud production kg/(m2*s)
        ! incoming snow + freezing of rain - sublimation - melting snow
        snowflxno_ls(i,j,k)    = zprvs(i,j) +  &
                    (srfrz(i,j) - ssdep(i,j)  - sgdep(i,j)   &
                                - ssmelt(i,j) - sgmelt(i,j)) &
                    * rhog * (hhl(i,j,k)-hhl(i,j,k+1))
        ! make a crude assumtion about the precipitating cloud cover
        ! if rain or snow exists, the cover is 1
        IF (qrg > 1.e-15_wp .OR. qsg > 1.e-15_wp) THEN
           preccover_ls(i,j,k) = 1._wp
        ELSE 
           preccover_ls(i,j,k) = 0._wp
        ENDIF
        sediice_ls(i,j,k) = 0._wp
#endif

      !------------------------------------------------------------------------
      ! Section 10: Complete time step
      !------------------------------------------------------------------------

      IF ( k /= ke ) THEN
        ! Store precipitation fluxes and sedimentation velocities ! for the next level
        zprvr(i,j) = qrg*rhog*zvzr(i,j)
        zprvs(i,j) = qsg*rhog*zvzs(i,j)
        zprvg(i,j) = qgg*rhog*zvzg(i,j)
        IF (zprvr(i,j) .LE. zqmin) zprvr(i,j)=0.0_wp
        IF (zprvs(i,j) .LE. zqmin) zprvs(i,j)=0.0_wp
        IF (zprvg(i,j) .LE. zqmin) zprvg(i,j)=0.0_wp

#ifdef NUDGING
        ! for the latent heat nudging
        IF ((llhn .OR. llhnverif) .AND. lhn_qrs) THEN
          qrsflux(i,j,k) = zprvr(i,j)+zprvs(i,j)+zprvg(i,j)
          qrsflux(i,j,k) = 0.5_wp*(qrsflux(i,j,k)+zpkr(i,j)+zpks(i,j)+zpkg(i,j))
        ENDIF
#endif

        IF (qrg+qr(i,j,k+1) .GT. zqmin) THEN
          ic1=ic1+1
          idx1(ic1)=i
          jdx1(ic1)=j
          helpr(ic1)=(qrg+qr(i,j,k+1))*0.5_wp*rhog
        ENDIF
        IF (qsg+qs(i,j,k+1) .GT. zqmin) THEN
          ic2=ic2+1
          idx2(ic2)=i
          jdx2(ic2)=j
          helps(ic2)=(qsg+qs(i,j,k+1))*0.5_wp*rhog
        ENDIF
        IF (qgg+qg(i,j,k+1) .GT. zqmin) THEN
          ic3=ic3+1
          idx3(ic3)=i
          jdx3(ic3)=j
          helpg(ic3)=(qgg+qg(i,j,k+1))*0.5_wp*rhog
        ENDIF
      ELSE
        ! Precipitation fluxes at the ground
        prr_gsp(i,j) = 0.5_wp * (qrg*rhog*zvzr(i,j) + zpkr(i,j))
        prs_gsp(i,j) = 0.5_wp * (qsg*rhog*zvzs(i,j) + zpks(i,j))
        prg_gsp(i,j) = 0.5_wp * (qgg*rhog*zvzg(i,j) + zpkg(i,j))

#ifdef NUDGING
        ! for the latent heat nudging
        IF ((llhn.OR.llhnverif) .AND. lhn_qrs)                     &
          qrsflux(i,j,k) = prr_gsp(i,j)+prs_gsp(i,j)+prg_gsp(i,j)
#endif
      ENDIF

      ! Update of prognostic variables or tendencies
      ! and add qi to qrs for water loading
      qr (i,j,k   ) = qrg
      qs (i,j,k   ) = qsg
      qg (i,j,k   ) = qgg
      qi (i,j,k   ) = qig
      qrs(i,j,k   ) = qrg+qsg+qig+qgg
      t  (i,j,k,nx) = t (i,j,k,nx) + ztt*zdt
      qv (i,j,k   ) = MAX ( 0.0_wp, qv(i,j,k) + zqvt*zdt )
      qc (i,j,k   ) = MAX ( 0.0_wp, qc(i,j,k) + zqct*zdt )

    ENDDO
  ENDDO

  IF ( k /= ke ) THEN
!CDIR BEGIN COLLAPSE
    zvzr(:,:)= 0.0_wp
    zvzs(:,:)= 0.0_wp
    zvzg(:,:)= 0.0_wp
!CDIR END
  END IF

!CDIR NODEP,VOVERTAKE,VOB
!CDIR ON_ADB(idx1)
!CDIR ON_ADB(jdx1)
  DO i1d = 1,ic1
    i = idx1(i1d)
    j = jdx1(i1d)
    zvzr(i,j)= zvz0r * EXP(zvzxp *LOG(helpr(i1d))) * zrho1o2(i,j)
  END DO
!CDIR NODEP,VOVERTAKE,VOB
!CDIR ON_ADB(idx2)
!CDIR ON_ADB(jdx2)
  DO i1d = 1,ic2
    i = idx2(i1d)
    j = jdx2(i1d)
    zvzs(i,j)= zvz0s(i,j) * EXP(ccswxp*LOG(helps(i1d))) * zrho1o2(i,j)
  END DO
!CDIR NODEP,VOVERTAKE,VOB
!CDIR ON_ADB(idx3)
!CDIR ON_ADB(jdx3)
  DO i1d = 1,ic3
    i = idx3(i1d)
    j = jdx3(i1d)
    zvzg(i,j)= zvz0g * EXP(zexpsedg*LOG(helpg(i1d))) * zrho1o2(i,j)
  END DO

  ! Do a final saturation adjustment for new values of t, qv and qc
!NEC_CB !CDIR COLLAPSE
!NEC_CB   zpres(:,:) = p0(:,:,k) + pp(:,:,k,nx)
!NEC_CB better do this in the very beginning instead of multiple times

  CALL satad ( 1, t(:,:,k,nx), qv(:,:,k),                 &
               qc(:,:,k), t(:,:,k,nx), zpres,             &
               zdummy(:,:,1),zdummy(:,:,2),zdummy(:,:,3), &
               zdummy(:,:,4),zdummy(:,:,5),zdummy(:,:,6), &
               zdummy(:,:,7),zdummy(:,:,8),               &
               b1, b2w, b3, b4w, b234w, rdv, o_m_rdv,     &
               rvd_m_o, lh_v, cpdr, cp_d,                 &
               ie, je, istartpar, iendpar, jstartpar, jendpar )

  IF ( ldiabf_lh ) THEN
    ! compute temperature increment due to latent heat
!CDIR COLLAPSE
    tinc_lh(:,:,k) = tinc_lh(:,:,k) + t(:,:,k,nx)
  ENDIF

ENDDO loop_over_levels

#ifdef NUDGING
! add part of latent heating calculated in subroutine hydci_pp_gr to model latent
! heating field: add temperature to model latent heating field
IF (llhn .OR. llhnverif) &
 CALL get_gs_lheating ('inc',1,ke)
#endif

#ifdef NECSX
  CALL collapse(.FALSE., ie, istartpar, iendpar, jstartpar, jendpar)
#endif
!------------------------------------------------------------------------------
! End of module procedure hydci_pp_gr
!------------------------------------------------------------------------------

END SUBROUTINE hydci_pp_gr

!==============================================================================

END MODULE src_gscp    
