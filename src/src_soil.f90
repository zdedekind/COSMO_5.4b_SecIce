!+ Source module  "src_soil"
!------------------------------------------------------------------------------

MODULE src_soil     

!------------------------------------------------------------------------------
!
! Description:
!   The module "src_soil" performs calculations related to the parameterization
!   of soil processes. It contains the two soubroutines terra1.inc and 
!   terra2.inc which comprise the soil model of LM. All parametric scalar and 
!   array data for these soil model routines are defined in the data module 
!   data_soil.f.
!
!   All global variables of the model that are used by the soil model routines
!   terra1.inc and terra2.inc are imported by USE statements below.
!
!   The parameterization package has been provided by B. Ritter in a
!   Plug-compatible Fortran77-Version, which is based on the EM/DM soil model 
!   by E. Heise. Some technical modifications have been done for the 
!   F90 and the parallel Version:
!   Internal communication by common-blocks is replaced by module parameters,
!   scalars and arrays defined in module data_soil. The plug compatible
!   I/O lists of the subroutines have been replaced by the Module interface
!   defined by Use lists below.
!
! Current Code Owner: DWD, Juergen Helmert
!  phone:  +49  69  8062 2704
!  fax:    +49  69  8062 3721
!  email:  Juergen.Helmert@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Guenther Doms     
!  Initial release
! 1.4        1998/05/22 Guenther Doms
!  Inclusion of the control parameter l2tls to select the
!  time levels according the the time integration scheme used.
! 1.20       1999/01/07 Guenther Doms
!  Renaming of some global variables
! 1.30       1999/06/24 Erdmann Heise
!  Implementation of simplified BATS-version for the computation of
!  evapotranspiration
! 1.33       1999/10/14 Matthias Raschendorfer
!  Removal of a LOGICAL namelist-parameter (lbats).
!  Introduction of 2 INTEGER-namelist-parameters controlling the evaporation:
!  These are: itype_(trvg, evsl), rat_lam.
! 2.2        2000/08/18 Matthias Raschendorfer
!  Declaration of the fields (sai, eai and tai).
!  USE of the INTEGER-namelist-parameter itype_tran.
! 2.8        2001/07/06 Ulrich Schaettler
!  Put include-files to this file and introduced optimization for 
!  vectorization in subroutine terra1
! 2.17       2002/05/08 Ulrich Schaettler
!  Replaced definition of snow and water covered fraction analogous to 
!  organize_radiation
! 2.18       2002/07/16 Reinhold Schrodin
!  Eliminated variable rhde, use cf_snow instead
! 3.18       2006/03/03 Ulrich Schaettler
!  Adaptations to do some initializations also for restart runs (nstart > 0)
! 3.21       2006/12/04 Ulrich Schaettler
!  crsmin, rat_lam put to data_soil; eliminated variables that are not used
!  Additional use of graupel, if present (Thorsten Reinhardt)
! V3_23        2007/03/30 Matthias Raschendorfer
!  Introduction of 'tfv' to consider different laminar resistance for heat
!  and water vapour; thus 'rat_lam' is not used here any longer.
! V4_12        2010/05/11 Ulrich Schaettler
!  Renamed t0 to t0_melt because of conflicting names
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_18        2011/05/26 Ulrich Schaettler
!  Changed the code owner
! V4_20        2011/08/31 Matthias Raschendorfer
!  Warning if SCLM mode is used with surface flux forcing (with ifdef)
! V4_23        2012/05/10 H.-J. Panitz, IMK/TRO, U. Boehm
!  Added new field snow_melt
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer, Ulrich Blahak
!  Replaced qx-variables by using them from the tracer module
!  Added hail rate (in case of two-moment microphysics) to the
!   precipitation quantities at the ground where it seems necessary.
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  MESSy interface introduced
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
!  Use new global fields: fr_snow, fr_wi (were local zf_snow, zf_wi before) and ustar_fv
! V5_1         2014-11-28 Matthias Raschendorfer, Oliver Fuhrer
!  Formal Introduction of the SC-framework
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

USE data_modelconfig, ONLY :   &

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------
    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction

! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the 
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from 
!    the other ones because of the use of the staggered Arakawa-C-grid.
!    
    istartpar,    & ! start index for computations in the parallel program
    iendpar,      & ! end index for computations in the parallel program
    jstartpar,    & ! start index for computations in the parallel program
    jendpar,      & ! end index for computations in the parallel program

! 4. variables for the time discretization and related variables
! --------------------------------------------------------------
    dt,           & ! long time-step
    dt2,          & ! dt*2.            

! 8. Organizational variables to handle the COSMO humidity tracers
! ----------------------------------------------------------------
    idt_qv

! end of data_modelconfig

!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 1. mathematical constants
! -------------------------
    pi,           & ! circle constant

! 2. physical constants and related variables
! -------------------------------------------

    t0_melt,      & ! absolute zero for temperature
    r_d,          & ! gas constant for dry air    
    rdv,          & ! r_d / r_v
    o_m_rdv,      & ! 1 - r_d/r_v
    rvd_m_o,      & ! r_v/r_d -1  
    cp_d,         & ! specific heat of dry air at constant pressure
    rdocp,        & ! r_d / cp_d
    lh_v,         & ! latent heat of vapourization
    lh_f,         & ! latent heat of fusion         
    lh_s,         & ! latent heat of sublimation    
    g,            & ! acceleration due to gravity
    sigma,        & ! Boltzmann-constant

! 3. constants for parametrizations
! ---------------------------------
    b1,           & ! variables for computing the saturation vapour pressure
    b2w,          & ! over water (w) and ice (i)
    b2i,          & !               -- " --
    b3,           & !               -- " --
    b4w,          & !               -- " --
    b4i,          & !               -- " --
    rho_w           ! density of liquid water  

! end of data_constants

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------
    p0         ,    & ! base state pressure                           (Pa) 

! 2. external parameter fields                                        (unit)
! ----------------------------
    soiltyp    ,    & ! type of the soil (keys 0-9)                     --
    plcov      ,    & ! fraction of plant cover                         --
    rootdp     ,    & ! depth of the roots                            ( m  )
    sai        ,    & ! surface area index                              --
    tai        ,    & ! transpiration area index                        --
    eai        ,    & ! earth area (evaporative surface area) index     --
    llandmask  ,    & ! landpoint mask

! 3. prognostic variables                                             (unit)
! -----------------------
    u          ,    & ! zonal wind speed                              ( m/s )
    v          ,    & ! meridional wind speed                         ( m/s )
    t          ,    & ! temperature                                   (  k  )
    pp         ,    & ! deviation from the reference pressure         ( pa  )

! 5. fields for surface values and soil model variables               (unit )
! -----------------------------------------------------
    ps        ,     & ! surface pressure                              ( pa  )
    t_snow    ,     & ! temperature of the snow-surface               (  K  )
    t_s       ,     & ! temperature of the ground surface             (  K  )
    t_g       ,     & ! weighted surface temperature                  (  K  )
    qv_s      ,     & ! specific humidity at the surface              (  K  )
    t_m       ,     & ! temperature between upper and medium layer    (  K  )
    t_cl      ,     & ! temperature beween medium and lowest layer    (  K  )
    w_snow    ,     & ! water content of snow                         (m H2O)
    fr_snow   ,     & ! surface fraction covered by snow           (  -  )
    fr_wi     ,     & ! surface fraction covered by interception water (-)
    ustar_fv  ,     & ! friction velocity (ustar)                  ( m/s )
    w_i       ,     & ! water content of interception water           (m H2O)
    w_g1      ,     & ! water content of the upper soil layer         (m H2O)
    w_g2      ,     & ! water content of the medium soil layer        (m H2O)
    w_g3      ,     & ! water content of the lower soil layer         (m H2O)
    w_cl              ! water content of the climatological layer     (m H2O)

USE data_fields     , ONLY :   &

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------
!   fields of convective and grid-scale precipitation
    prr_con     ,   & ! precipitation rate of rain, convective        (kg/m2*s)
    prs_con     ,   & ! precipitation rate of snow, convective        (kg/m2*s)
    prr_gsp     ,   & ! precipitation rate of rain, grid-scale        (kg/m2*s)
    prs_gsp     ,   & ! precipitation rate of snow, grid-scale        (kg/m2*s)
    prg_gsp     ,   & ! precipitation rate of graupel, grid-scale     (kg/m2*s)
    prh_gsp     ,   & ! precipitation rate of hail, grid-scale        (kg/m2*s)

!   fields of the turbulence scheme
    tch         ,   & ! turbulent transfer coefficient for heat       ( -- )
    tcm         ,   & ! turbulent transfer coefficient for momentum   ( -- )
    tfv         ,   & ! laminar reduction factor for evaporation      ( -- )

!   fields of the radiation
    sobs        ,   & ! solar radiation at the ground                 ( W/m2)
    thbs        ,   & ! thermal radiation at the ground               ( W/m2)
    pabs        ,   & ! photosynthetic active radiation               ( W/m2)

! 7. fields for model output and diagnostics                          (unit )
! ---------------------------------------------------------------
    runoff_s    ,   & ! surface water runoff; summ over forecast      (kg/m2)
    runoff_g    ,   & ! soil water runoff; summ over forecast         (kg/m2)
    rstom       ,   & ! stomata resistance                           ( s/m )
    snow_melt         ! snow melt amount; summ over forecast          (kg/m2)

! end of data_fields

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    nstart,       & ! first time step of the forecast
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
    nold  ,       & ! corresponds to ntstep - 1
    nnow  ,       & ! corresponds to ntstep 
    nnew  ,       & ! corresponds to ntstep + 1

! 3. controlling the physics
! --------------------------
    itype_gscp,   & ! type of grid-scale precipitation physics
    nlgw ,        & ! number of prognostic soil water levels
    itype_trvg,   & ! type of vegetation transpiration parameterization
    itype_evsl,   & ! type of parameterization of bare soil evaporation
    itype_tran,   & ! type of surface to atmospher transfer

! 5. additional control variables
! --------------------------
     lscm,        &  ! SCM-run
     l2tls           ! forecast with 2-TL integration scheme 

! end of data_runcontrol 

!------------------------------------------------------------------------------

USE data_soil       ! All variables from data module "data_soil are used by
                    ! this module. These variables start with letter "c"
                    ! (except for epsilons, which start with "eps_").
! end of data_soil       

!------------------------------------------------------------------------------

! External routines uses by this module (in subroutine terra2)

USE meteo_utilities, ONLY : tgcom

!------------------------------------------------------------------------------

USE data_parallel,   ONLY : my_cart_id

!------------------------------------------------------------------------------

USE environment,     ONLY : model_abort

!------------------------------------------------------------------------------

USE src_tracer,      ONLY: trcr_get, trcr_errorstr

!------------------------------------------------------------------------------

#ifdef SCLM
USE data_1d_global, ONLY : &

    lertflu, i_cal, i_upd, i_mod, im, jm, &

    SNC, INC, & !snow_cover and intc_cover
    SHF, LHF    !surface heat fluxes
#endif

!------------------------------------------------------------------------------
! Declarations
!------------------------------------------------------------------------------

IMPLICIT NONE

! The following arrays are declared for communication between the subroutines
! terra1 and terra2 which comprise the soil model. They are allocated and
! precalculated in terra1. In module subroutine terra2, these arrays are
! modified and finally deallocated.

REAL (KIND = wp),     ALLOCATABLE, PRIVATE  :: &
  xdzs     (:,:)   , & ! 
  xalam    (:,:)   , & ! heat conductivity
  xdqvts   (:,:)   , & ! 
  xdqvtsnow(:,:)   , & ! 
  xdwidt   (:,:)   , & ! interception store tendency 
  xdwsndt  (:,:)   , & ! snow water contenttendency
  xesoil   (:,:)   , & ! evaporation from bare soil
  xrr      (:,:)   , & ! total rain rate including formation of dew
  xrs      (:,:)   , & ! total snow rate including formation of rime
  xrhoch   (:,:)   , & ! transfer coefficient*rho*g
  xrocg    (:,:)   , & ! heat capacity of first soil layer 
  xrocm    (:,:)   , & ! heat capacity of sec. soil layer
  xrocs    (:,:)   , & ! heat capacity of snow
  xth_low  (:,:)   , & ! 
  xtrang   (:,:,:)     ! transpiration from various ground layers


!     Definition of temperature related variables in the soil model
!
!
!     --------------T_snow--------
!                                | 
!           fr_snow              |            1 - fr_snow
!                                |      
!     ______________T_s__________|________________T_s_____________
!     ////////////////////////////////////////////////////////////
!
!
!
!     _________________________________T_m________________________
!
!
!
!
!
!
!
!
!
!     _________________________________T_cl_______________________
!
!
!
!     The surface temperature T_g is the snow fraction (fr_snow)
!     weighted sum of T_snow and T_s!
!==============================================================================

CONTAINS

!==============================================================================
!+ Computation of the first part of the soil parameterization scheme
!------------------------------------------------------------------------------

SUBROUTINE terra1             

!------------------------------------------------------------------------------
!
! Description:
!   The module procedure terra1 performs the first part of the soil processes
!   parameterization scheme. In this part, evaporation from the surface is 
!   calculated and a number of array variables is provided for later use in 
!   the module procedure terra2. The communication between these modules is 
!   by module arrays which are allocated in terra1 and are deallocated at the 
!   end of module procedure terra2. Following the execution of terra1, a 
!   moisture convergence based convection scheme may be called. The resulting 
!   convective precipitation rate is then used in terra2 for the final 
!   updating of all prognostic variables of the soil model.
!
! Method:
!   The soil model is based on the "Extendend force restore" (EFR) method for
!   the prediction of soil temperature in a two-layer scheme following Jacobsen
!   and Heise (1982, Contr.Atm.Phys.,55). The prediction of soil water is 
!   calculated by a) budget equations with a two-layer or (optional) three-layer
!   soil moisture scheme, or b) by a simplified BATS-version with two or three
!   soil layers.
!
!   For a detailed description of the physical and mathematical methods,
!   see the EM/DM Documentation.
!
!------------------------------------------------------------------------------
! Subroutine arguments: None
! --------------------
!
! Local parameters:
! ----------------

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    izerror,        & ! error number
    nx     ,        & ! Time-level for integration
    kwl    ,        & ! loop index for soil moisture layers           
    i      ,        & ! loop index in x-direction              
    j      ,        & ! loop index in y-direction              
    im1, jm1,       & ! i-1, j-1
    jb     ,        & ! loop index for soil-type               
    mstyp  ,        & ! soil type index
    istat  ,        & ! error status for allocation      
    istarts,        & ! start index for x-direction      
    iends  ,        & ! end   index for x-direction     
    jstarts,        & ! start index for y-direction     
    jends             ! end   index for y-directin      

  REAL    (KIND=wp   )     ::  &
    zdt    ,        & ! integration time step
    zwmmin ,        & !
    zwgmin ,zwgmax, & !
    zwmmax ,        & !
    zx     ,        & !
    zomb   ,        & !
    zzdlam ,        & !
    zwmean ,        & !
    zrhoc  ,        & !
    zzalam ,        & !
    zdtdrhw,        & ! 
    zrhwddt,        & !
    z1d2dt ,        & !  1./2*timestep
    zuv    ,        & ! wind velocity     in lowest atmospheric layer
    ztvs   ,        & !
    zplow  ,        & ! pressure          of lowest atmospheric layer
    zqvlow ,        & ! specific humidity of lowest atmospheric layer
    zqs    ,        & ! saturation humidty at T_s and ps
    zdqs   ,        & ! derivative of zqs with respect to T_s
    zqsnow ,        & ! saturation humidty at T_snow and ps
    zdqsnow,        & ! derivative of zqs with respect to T_snow
    zwqg   , zwqm,  & !
    z4wdpv ,        & !
    zrss   , zrww,  & !
    zrhosn , zbeta, & !
    zice   ,        & ! indicator of soil type ice
    zevap  ,        & !
    ztlpmwp           !

  REAL    (KIND=wp   )     ::  &
    ztgt0  , zroot, & !
    zr_root,        & !
    ze_sum ,        & ! sum of all contributions to evapotranspiration
    zsf_heav,       & ! Statement function: Heaviside function
    zsf_psat_iw,    & ! Saturation water vapour pressure over ice or water 
                      ! depending on temperature "zstx"
    zsf_qsat,       & ! Specific humidity at saturation pressure 
                      ! (depending on the saturation water vapour pressure 
                      !  "zspsatx" and the air pressure "zspx")
    zsf_dqvdt_iw,   & ! Statement function: First derivative of specific  
                      ! saturation humidity 
                      ! with respect to temperature (depending on temperature 
                      ! "zstx" and saturation specific humidity pressure 
                      ! "zsqsatx")
    zstx   ,        & ! dummy argument for Stmt. function
    zspx   ,        & ! dummy argument for Stmt. function
    zspsatx,        & ! dummy argument for Stmt. function
    zsqsatx,        & ! dummy argument for Stmt. function
    z2iw   ,        & ! dummy argument for Stmt. function
    z4iw   ,        & ! dummy argument for Stmt. function
    z234iw            ! dummy argument for Stmt. function

  REAL    (KIND=wp   )     ::  &
!   New Local scalars added for BATS-scheme
    zbf1   ,        & ! auxiliary variable
    zbf2   ,        & ! auxiliary variable
    zdmax  ,        & ! auxiliary variable
    zd     ,        & ! auxiliary variable
    zck    ,        & ! auxiliary variable
    zs1    ,        & ! auxiliary variable
    zfqmax ,        & ! maximum sustainable flux in uppermost soil layer
    zevapor,        & ! evaporation rate of bare soil in BATS-scheme
    zrla   ,        & ! atmospheric resistance
    zpar   ,        & ! PAR (interim value)
    zfill  ,        & ! = 1, if layer is totally filled by roots, = 0 otherwise
    zpart  ,        & ! = 1, if layer is partly filled by roots, = 0 otherwise
    zropart,        & ! fraction of layer filled by roots
    zf_rad ,        & ! radiation function for stomatal resistance
    zf_wat ,        & ! soil water function for stomatal resistance
    zf_tem ,        & ! temperature function for stomatal resistance
    zf_sat ,        & ! saturation deficit function for stomatal resistance
    zepsat ,        & ! saturation vapour pressure at near surface temperature
    zepke  ,        & ! near surface water vapour pressure
    zrstom ,        & ! stomata resistance
    zedrstom,       & ! 1./zrstom
    ztrabpf,        & ! area average transpiration rate
    zr1, zr2          ! for intermediate storage

! Local (automatic) arrays:
! -------------------------
  REAL    (KIND=wp   )     ::  &
    zadp     (ie,je)  , & ! air dryness point
    zbwt     (ie,je)  , & ! root depth (with artificial minimum value)
    zrock    (ie,je)  , & ! ice/rock-indicator
    zporv    (ie,je)  , & ! 
    zts      (ie,je)  , & ! 
    zts_pm   (ie,je)  , & ! 
    ztsnow   (ie,je)  , & ! 
    ztsnow_pm(ie,je)  , & ! 
    zwin     (ie,je)  , & ! water cont. of interception store
    zwsnow   (ie,je)  , & ! 
    zfcap    (ie,je)  , & ! 
    zpwp     (ie,je)  , & ! 
    ztmch    (ie,je)  , & ! 
    zwg_fr   (ie,je,4), & ! fractional water content of layers in ground
    zdlam    (ie,je)  , & ! 
    zep_s    (ie,je)  , & ! potential evaporation for T_s    
    zep_snow (ie,je)  , & ! potential evaporation for T_snow 
    zdzwsu (4)        , & ! sum over hydrological layer thicknesses
    zdzwg  (4)        , & ! 
    znlgw1f(3)        , & !
!   New Local arrays added for BATS-scheme
    zk0di    (ie,je)  , & ! surface type dependent parameter
    zbedi    (ie,je)  , & ! surface type dependent parameter
    zsnull   (ie,je)  , & ! mean relativ (zporv) water content of the soil
    zcatm    (ie,je)  , & ! atmospheric transfer velocity
    zrveg    (ie,je)  , & ! additional resistance of the vegetation 
    zsdef_qv (ie,je)  , & ! saturation deficit
    ztraleav (ie,je)  , & ! transpiration rate of dry leaves
    zwroot   (ie,je)  , & ! mean water content over root depth
    zrootfc  (ie,je,3)    ! distributes transpiration to soil layers
  INTEGER  (KIND=iintegers ) ::  &
     m_styp (ie,je)       ! soil type

  CHARACTER(LEN=255)         ::  &
    yzerrmsg

  CHARACTER(LEN=80)          ::  &
    yzroutine

  REAL (KIND=wp),     POINTER :: &
    qv  (:,:,:)  => NULL()          ! QV at nx

!- End of header
!=======================================================================

!------------------------------------------------------------------------------
! Begin Subroutine terra1              
!------------------------------------------------------------------------------

! Declaration of STATEMENT-FUNCTIONS

  zsf_heav     (zstx                    ) = 0.5_wp+SIGN( 0.5_wp, zstx )
  zsf_psat_iw  (zstx,z2iw   ,z4iw       )                                     &
                   = b1*EXP(z2iw*(zstx - b3)/(zstx - z4iw))
  zsf_qsat     (zspsatx, zspx           )                                     &
                   = rdv*zspsatx/(zspx-o_m_rdv*zspsatx)
  zsf_dqvdt_iw (zstx,zsqsatx,z4iw,z234iw)                                     &
                   = z234iw*(1.0_wp+rvd_m_o*zsqsatx)*zsqsatx/(zstx-z4iw)**2

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  yzroutine = 'terra1'

  ! Allocation of arrays for the communication with the second part of the 
  ! soil model (terra2); these arrays have been defined in the module 
  ! declaration section

  ALLOCATE ( xdzs     (ie,je)  , STAT = istat ) ! 
  ALLOCATE ( xalam    (ie,je)  , STAT = istat ) ! heat conductivity
  ALLOCATE ( xdqvts   (ie,je)  , STAT = istat ) ! 
  ALLOCATE ( xdqvtsnow(ie,je)  , STAT = istat ) ! 
  ALLOCATE ( xdwidt   (ie,je)  , STAT = istat ) ! interception store tendency 
  ALLOCATE ( xdwsndt  (ie,je)  , STAT = istat ) ! snow water contenttendency
  ALLOCATE ( xesoil   (ie,je)  , STAT = istat ) ! evaporation from bare soil
  ALLOCATE ( xrr      (ie,je)  , STAT = istat ) ! total rain rate including 
                                                ! formation of dew
  ALLOCATE ( xrs      (ie,je)  , STAT = istat ) ! total snow rate including 
                                                ! formation of rime
  ALLOCATE ( xrhoch   (ie,je)  , STAT = istat ) ! transfer coefficient*rho*g
  ALLOCATE ( xrocg    (ie,je)  , STAT = istat ) ! heat capacity of first 
                                                ! soil layer 
  ALLOCATE ( xrocm    (ie,je)  , STAT = istat ) ! heat capacity of sec. 
                                                ! soil layer
  ALLOCATE ( xrocs    (ie,je)  , STAT = istat ) ! heat capacity of snow
  ALLOCATE ( xth_low  (ie,je)  , STAT = istat ) ! 
  ALLOCATE ( xtrang   (ie,je,3), STAT = istat ) ! transpiration from various 
                                                ! ground layers

! select timelevel and timestep for calculations
  IF ( l2tls ) THEN
    nx  = nnow
    zdt = dt
  ELSE
    nx  = nold
    zdt = dt2
  ENDIF

! Horizontal domain for computation
  istarts = istartpar
  iends   = iendpar
  jstarts = jstartpar
  jends   = jendpar

  ! retrieve the required microphysics tracers (at timelevel nx)
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nx, ptr = qv)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

! Computation of derived constants

  z1d2dt   = 1.0_wp/zdt         ! 1./2*timestep

! layer thickness for hydrological part
  IF (nlgw == 2) THEN    ! (2-layer model)
     zdzwg(1)  = cdzw12 
     zdzwg(2)  = cdzw22
  ELSE                      ! (3-layer model)
     zdzwg(1)  = cdzw13
     zdzwg(2)  = cdzw23
     zdzwg(3)  = cdzw33
  END IF
  zdzwg(nlgw+1)= zdzwg(nlgw)

  zdtdrhw = zdt/rho_w   ! 2*timestep/density of liquid water
  zrhwddt = rho_w/zdt   ! inverse of ...

  ! time constant for infiltration of water from interception store must not 
  ! be less than dt2
  ctau_i        = MAX(ctau_i,zdt)

! Specification of depth of soil layer centers
  zdzwsu(1) = 0.5_wp*zdzwg(1)
  DO kwl    = 1, nlgw
     znlgw1f(kwl)  = 0.0_wp
     IF (kwl < nlgw) THEN
        zdzwsu(kwl+1)= zdzwsu(kwl) + 0.5_wp*(zdzwg(kwl) + zdzwg(kwl+1))
     ENDIF
  ENDDO
  zdzwsu(nlgw+1)  = zdzwsu(nlgw) + 0.5_wp*zdzwg(nlgw)
  znlgw1f(1)         = 1.0_wp

! Prepare basic surface properties (for land-points only)
 
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF(llandmask(i,j)) THEN        ! for land-points only
        mstyp       = NINT(soiltyp(i,j))        ! soil type
        m_styp(i,j) = mstyp                     ! array for soil type
        xdzs (i,j)  = cdsmin                    ! minimum snow depth
        zporv(i,j)  = cporv(mstyp)              ! pore volume
        zpwp (i,j)  = cpwp (mstyp)              ! plant wilting point
        zadp (i,j)  = cadp (mstyp)              ! air dryness point
        zfcap(i,j)  = cfcap(mstyp)              ! field capacity
        zrock(i,j)  = crock(mstyp)              ! rock or ice indicator
        xrocg(i,j)  = crhoc(mstyp)              ! heat capacity
        zdlam(i,j)  = cala1(mstyp)-cala0(mstyp) ! heat conductivity parameter
        xalam(i,j)  = cala0(mstyp)              ! heat conductivity parameter
        ! New arrays for BATS-scheme
        zk0di(i,j)  = ck0di(mstyp)              !
        zbedi(i,j)  = cbedi(mstyp)              !
        zbwt(i,j)   = MAX(0.1_wp,rootdp(i,j))! Artificial minimum value for 
                                                ! root depth
      ENDIF
    ENDDO
  ENDDO

! For ntstep=0 : a) determine layer thickness of all soil types for EFR-method
!                b) determine constant clgk0 for BATS-scheme
! ntstep=nstart: c) confine water content to range between air dryness point 
!                   and pore volume

  IF (ntstep == nstart) THEN

!--------------------------------------------------------------------------------------
#ifdef SCLM
    IF (lscm .AND. (LHF%mod(0)%ist.EQ.i_mod .OR. SHF%mod(0)%ist.EQ.i_mod)) THEN
       WRITE (*,*) 'ATTENTION: SURFACE FLUX FORCING IS NOT YET IMPLEMENTED IN THE '// &
                   '2-LAYER FORCE RESTORE SOIL MODEL!'
    END IF
#endif
!SCLM ---------------------------------------------------------------------------------

    zx         = SQRT(ctau1/ctau2)
    zomb       = 2.0_wp*pi/(ctau1*86400.0_wp)
    DO jb      = 1, 10
      zzdlam   = cala1(jb) - cala0(jb)
      zwmean   = 0.5_wp*(cfcap(jb) + cpwp(jb))
      z4wdpv   = 4.0_wp*zwmean/cporv(jb)
      zrhoc    = crhoc (jb) + rho_w*zwmean*chc_w
      zzalam   = cala0(jb) + zzdlam*(0.25_wp + 0.30_wp*zzdlam   &
                 /(1.0_wp + 0.75_wp*zzdlam)) &
                 *MIN(z4wdpv,1.0_wp+(z4wdpv-1.0_wp)  &
                 *(1.0_wp + 0.35_wp*zzdlam) &
                 /(1.0_wp + 1.95_wp*zzdlam) )
      cdz1(jb) = SQRT (2.0_wp*zzalam/(zrhoc*zomb))/(1.0_wp+zx)
      clgk0(jb) = LOG10(MAX(eps_soil,ck0di(jb)/ckrdi))
    ENDDO

  ENDIF

  ! it has to be checked, if it is the first step after a forward
  ! launching digital filtering initialization: then lsoilinit_dfi
  ! is set to TRUE in dfi_initialization
  IF ((ntstep == 0) .OR. lsoilinit_dfi) THEN
    lsoilinit_dfi = .FALSE.

    DO   j     = jstarts, jends
      DO i     = istarts, iends
        IF (llandmask(i,j)) THEN              ! for land-points only
          zwgmin           = zadp (i,j)*zdzwg(1)    ! lower-limit for 1. layer
          zwgmax           = zporv(i,j)*zdzwg(1)    ! upper-limit for 1. layer
          zwmmin           = zadp (i,j)*zdzwg(2)    ! lower-limit for 2. layer
          zwmmax           = zporv(i,j)*zdzwg(2)    ! upper-limit for 2. layer
          w_g1(i,j,nx  )   = MAX( zwgmin, MIN(zwgmax,w_g1(i,j,nx  )) )
          w_g1(i,j,nnow)   = MAX( zwgmin, MIN(zwgmax,w_g1(i,j,nnow)) )
          w_g2(i,j,nx  )   = MAX( zwmmin, MIN(zwmmax,w_g2(i,j,nx  )) )
          w_g2(i,j,nnow)   = MAX( zwmmin, MIN(zwmmax,w_g2(i,j,nnow)) )
          IF (nlgw == 2) THEN
            w_cl(i,j)      = MAX(zwmmin,MIN(zwmmax,w_cl(i,j)    ))
          ELSE
            zwmmin         = zadp (i,j)*zdzwg(3)    ! lower-limit for 3. layer
            zwmmax         = zporv(i,j)*zdzwg(3)    ! upper-limit for 3. layer
            w_g3(i,j,nx  ) = MAX( zwmmin, MIN(zwmmax,w_g3(i,j,nx  )) )
            w_g3(i,j,nnow) = MAX( zwmmin, MIN(zwmmax,w_g3(i,j,nnow)) )
            w_cl(i,j)      = MAX( zwmmin, MIN(zwmmax,w_cl(i,j)     ) )
          END IF
        END IF
      END DO
    END DO

  END IF             ! End of 'initial time step'-calculations

! conversion of tch to tmch 

  DO   j = jstarts, jends
    jm1  = MAX( 1, j-1 )
    DO i = istarts, iends
       im1        = MAX( 1, i-1)
       zuv        = 0.5_wp*SQRT ( (u(i,j,ke,nx) + u(im1,j,ke,nx))**2 &
                                 +(v(i,j,ke,nx) + v(i,jm1,ke,nx))**2 )
       ztvs       = t_s (i,j,nx)*(1.0_wp + rvd_m_o*qv_s(i,j,nx))
       ztmch(i,j) = tch(i,j)*zuv*g*ps(i,j,nx)/(r_d*ztvs)
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section 2: temperatures, water contents (in mH2O), surface pressure, 
!------------------------------------------------------------------------------

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN     ! for land-points only

         IF (w_snow(i,j,nx) > 0.0_wp) THEN                           
            ! existence of snow 
            ! --> no water in interception store and t_snow < t0_melt
            w_i   (i,j,nx) = 0.0_wp
            t_snow(i,j,nx) = MIN (t0_melt - eps_soil, t_snow(i,j,nx) )
         ELSE IF (t_snow(i,j,nx) >= t0_melt) THEN                      
            ! no snow and t_snow >= t0_melt --> ts > t0_melt and t_snow = ts
            t_s   (i,j,nx) = MAX (t0_melt + eps_soil, t_s(i,j,nx) )
            t_snow(i,j,nx) = t_s(i,j,nx)                        
         ELSE
            ! no snow and  t_snow < t0_melt 
            ! --> Ts < Tmelt and T_snow = ts and no water in interception store
            t_s   (i,j,nx) = MIN (t0_melt - eps_soil, t_s(i,j,nx) )
            t_snow(i,j,nx) = t_s(i,j,nx)
            w_i   (i,j,nx) = 0.0_wp
         END IF

         IF (nlgw == 2) THEN   ! two soil layers for hydrology
            zwg_fr(i,j,1) = w_g1(i,j,nx)/zdzwg(1)
            zwg_fr(i,j,2) = w_g2(i,j,nx)/zdzwg(2)
            zwg_fr(i,j,3) = w_cl (i,j     )/zdzwg(2)
         ELSE                     ! three soil layers for hydrology
            zwg_fr(i,j,1) = w_g1(i,j,nx)/zdzwg(1)
            zwg_fr(i,j,2) = w_g2(i,j,nx)/zdzwg(2)
            zwg_fr(i,j,3) = w_g3(i,j,nx)/zdzwg(3)
            zwg_fr(i,j,4) = w_cl (i,j     )/zdzwg(3)
         END IF           

      END IF             
    ENDDO
  ENDDO
      
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN   ! for land-points only
        zplow          = p0(i,j,ke) + pp(i,j,ke,nx)
        ztsnow   (i,j) = t_snow(i,j,nx)
        zts      (i,j) = t_s   (i,j,nx)
        zts_pm   (i,j) = zsf_heav(zts   (i,j) - t0_melt)
        ztsnow_pm(i,j) = zsf_heav(ztsnow(i,j) - t0_melt)
        zwsnow   (i,j) = w_snow(i,j,nx)
        zwin     (i,j) = w_i   (i,j,nx)

        ! moisture and potential temperature of lowest layer
        zqvlow         = qv(i,j,ke)
        xth_low (i,j)  =  t(i,j,ke,nx) *( (ps(i,j,nx)/zplow)**rdocp )

        ! density*transfer coefficient*wind velocity 
        xrhoch(i,j)    = ztmch(i,j)*(1.0_wp/g) + eps_soil

        ! saturation specific humidity for ts and t_snow and first derivative
        z2iw        = zts_pm(i,j)*b2w + (1.0_wp - zts_pm(i,j))*b2i
        z4iw        = zts_pm(i,j)*b4w + (1.0_wp - zts_pm(i,j))*b4i
        z234iw      = z2iw*(b3 - z4iw)
        zqs         = zsf_qsat( zsf_psat_iw(zts(i,j), z2iw,z4iw), ps(i,j,nx) )
        xdqvts(i,j) = zsf_dqvdt_iw (zts(i,j), zqs, z4iw, z234iw)
        zdqs        = zqvlow - zqs   
        IF (ABS(zdqs).Lt.eps_soil) zdqs = 0.0_wp
        z2iw        =         ztsnow_pm(i,j) *b2w + (1.0_wp - ztsnow_pm(i,j))*b2i
        z4iw        =         ztsnow_pm(i,j) *b4w + (1.0_wp - ztsnow_pm(i,j))*b4i
        z234iw      = z2iw*(b3 - z4iw)
        zqsnow      = zsf_qsat(zsf_psat_iw(ztsnow(i,j),z2iw,z4iw), ps(i,j,nx))
        xdqvtsnow(i,j)= zsf_dqvdt_iw(ztsnow(i,j), zqsnow, z4iw,z234iw)
        zdqsnow     = zqvlow - zqsnow
        IF (ABS(zdqsnow).Lt.eps_soil) zdqsnow = 0.0_wp

        ! potential evaporation at T_snow and Ts
        zep_snow(i,j) = (1.0_wp-ztsnow_pm(i,j))*tfv(i,j)*xrhoch(i,j)*zdqsnow
        zep_s   (i,j) =                         tfv(i,j)*xrhoch(i,j)*zdqs
      END IF               
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section 3: heat conductivity, frozen fraction, snow and 
!            water covered fraction snow height, volume heat content
!------------------------------------------------------------------------------

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only
        zwqg         = 0.5_wp*(zfcap(i,j) + zpwp(i,j))  
        zwqm         = zwqg
        z4wdpv       = 4.0_wp*zwqg/zporv(i,j)

        ! heat conductivity
        xalam(i,j)   = xalam(i,j) + zdlam(i,j)                                &
                      * (0.25_wp + 0.30_wp*zdlam(i,j) / (1.0_wp+0.75_wp*zdlam(i,j)))  &
                      * MIN (z4wdpv, 1.0_wp + (z4wdpv-1.0_wp)    &
                      *(1.0_wp+0.35_wp*zdlam(i,j)) &
                      /(1.0_wp+1.95_wp*zdlam(i,j)))

        ! volume heat content (layer thickness not included yet)
        xrocm(i,j)   = xrocg(i,j) + rho_w*zwqm*chc_w    ! attention: do not
        xrocg(i,j)   = xrocg(i,j) + rho_w*zwqg*chc_w    ! swap statements  

        ! snow and water covered fraction
        ! modified at 15.02.2002 (adapted to GME formulation, identical with
        ! zsnow = in organize_radiation.f90)
        zrss = MAX(0.01_wp, MIN(1.0_wp,zwsnow(i,j)/cf_snow) )
        zrww =       ztsnow_pm(i,j)  * MAX(0.01_wp,                       &
                    1.0_wp - EXP ( MAX (-5.0_wp, -zwin  (i,j)/cf_w) ) )
        IF (zrss > 0.99_wp) zrss = 1.0_wp
        IF (zrww > 0.99_wp) zrww = 1.0_wp
        fr_snow(i,j) = zrss*zsf_heav(zwsnow(i,j) - 0.5_wp*eps_soil)
        fr_wi  (i,j) = zrww*zsf_heav(zwin  (i,j) - 0.5_wp*eps_soil)

        ! density of snow as function of snow water content and snow height
        !  confined to range cdsmin to 1.5m for temperature calculation
        zrhosn     = MAX( crhosmin,                                           &
                          MIN( crhosmax, crhosmin + crhos_dw*zwsnow(i,j) ) )
        xdzs (i,j) = MIN(1.5_wp, MAX(cdsmin,                              &
                zwsnow(i,j) * rho_w / zrhosn / MAX(0.01_wp,fr_snow(i,j))))
        xrocs(i,j) = chc_i*xdzs(i,j)*zrhosn
      END IF         
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section 4: Hydrology, 1.Section
!------------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Section 4.1: Evaporation
  !----------------------------------------------------------------------------
  
  ! Evaporation from interception store and snow
  ! Formation of dew and rime
  ! Evaporation and transpiration are negative, dew and rime
  ! positive quantities, since positive sign indicates a flux
  ! directed towards the earth's surface!
  
    DO   j = jstarts, jends
      DO i = istarts, iends
  
        xdwidt (i,j)     = 0.0_wp ! Initialisation of all
        xdwsndt(i,j)     = 0.0_wp ! evaporation quantities
        xesoil (i,j)     = 0.0_wp ! to be equal to zero
        xrr    (i,j)     = 0.0_wp ! 
        xrs    (i,j)     = 0.0_wp ! 
        xtrang (i,j,1)   = 0.0_wp ! Initialisation of 
        xtrang (i,j,2)   = 0.0_wp ! transpiration
        xtrang (i,j,3)   = 0.0_wp ! quantities    
  
        IF (llandmask(i,j)) THEN  ! land points only
  
          ! Evaporation from interception store if it contains water (wi>0) and
          ! if zep_s<0 indicates potential evaporation for temperature Ts
          ! amount of water evaporated is limited to total content of store
          xdwidt(i,j) = zsf_heav(-zep_s(i,j))      &
                        * MAX(-zrhwddt*zwin(i,j), fr_wi(i,j)*zep_s(i,j))
  
          ! Evaporation of snow, if snow exists (wsnow>0) and if zep_snow<0
          ! indicates potential evaporation for temperature Tsnow
          xdwsndt(i,j) = zsf_heav(-zep_snow(i,j))  &
                         * MAX(-zrhwddt*zwsnow(i,j), fr_snow(i,j)*zep_snow(i,j))
  
          ! Formation of dew or rime, if zep_s > 0 . distinction between
          ! dew or rime is only controlled by sign of surface temperature
          ! and not effected by presence of snow !
          xrr(i,j) = zsf_heav(zep_s   (i,j))*zep_s   (i,j)*    zts_pm(i,j)
          xrs(i,j) = zsf_heav(zep_snow(i,j))*zep_snow(i,j)*(1.0_wp-zts_pm(i,j))
        END IF   
      ENDDO
    ENDDO
  
  
  !----------------------------------------------------------------------------
  ! Section 4.2a: Bare soil evaporation, bucket version
  !----------------------------------------------------------------------------
  
  IF (itype_evsl.EQ.1) THEN   ! Bucket version
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN       ! land-points only
          IF (zep_s(i,j) < 0.0_wp) THEN  ! upwards directed potential evaporation 
            ! reduction factor for evaporation based on water content of 
            ! first soil layer, air dryness point and field capacity
            zbeta = MAX( 0.0_wp, MIN( 1.0_wp,                        &
                           (zwg_fr(i,j,1)-zadp(i,j)) / (zfcap(i,j)-zadp(i,j))))
            zbeta = zbeta**2
  
            ! if first soil layer is frozen, allow evaporation at potential
            ! rate; if soil type is rock, evaporation is not allowed at all
            zice        = zsf_heav(1.5_wp - m_styp(i,j)) ! =1 only for ice
            zevap       = zrock(i,j) + zice          ! =1 for all, but rock
            zbeta       = zbeta + (1.0_wp - zbeta)*zice
            xesoil(i,j) = zevap*zbeta                  & ! reduction
                          * zep_s(i,j)                 & ! evaporation
                          * (1.0_wp - fr_wi  (i,j))    & ! not water covered
                          * (1.0_wp - fr_snow(i,j))    & ! not snow covered
                          * eai(i,j)/sai(i,j)        ! relative source surface of the bare soil
          END IF ! upwards directed potential evaporation
        END IF   ! land points
      END DO
    END DO
  END IF         ! Bucket version
  
  !----------------------------------------------------------------------------
  ! Section 4.2b: Bare soil evaporation, BATS version
  !----------------------------------------------------------------------------

  IF (itype_evsl.EQ.2) THEN
    ! Calculation of bare soil evaporation after Dickinson (1984)
    ! Determination of mean water content relative to volume of voids
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN   ! land points only
          zsnull(i,j)   = 0.0_wp
        ENDIF
      ENDDO
    ENDDO
    DO kwl = 1,nlgw
      DO   j = jstarts, jends
        DO i = istarts, iends
          IF (llandmask(i,j)) THEN   ! land points only
            zsnull(i,j) = zsnull(i,j) + zwg_fr(i,j,kwl)*zdzwg(kwl)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN   ! land points only
          zsnull(i,j)   = zsnull(i,j)/(zporv(i,j)*zdzwsu(nlgw + 1))
        END IF
      END DO
    END DO

    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN       ! land points only
          IF (zep_s(i,j) < 0.0_wp) THEN   ! upwards directed potential evaporation
            ! Treatment of Ice (m_styp=1) and Rocks (m_styp=2)
            zice   = zsf_heav(1.5_wp - m_styp(i,j)) ! = 1 only for ice
            zevap  = zrock(i,j) + zice                  ! = 1 for all soil types
                                                        !   but rock (=0)
            zbeta  = 0.0_wp
            IF (m_styp(i,j).ge.3) THEN ! Computations not for ice and rocks
              ! auxiliary quantities
              zbf1   = 5.5_wp - 0.8_wp* zbedi(i,j)*                &
                      (1.0_wp + 0.1_wp*(zbedi(i,j) - 4.0_wp)*  &
                       clgk0(m_styp(i,j)) )
              zbf2   = (zbedi(i,j) - 3.7_wp + 5.0_wp/zbedi(i,j))/  &
                      (5.0_wp + zbedi(i,j))
              zdmax  = zbedi(i,j)*cfinull*zk0di(i,j)/crhowm
              zs1    = zwg_fr(i,j,1)/zporv(i,j)
!             zd     = 1.02_wp*zdmax*zs1**zbedi(i,j)*zsnull(i,j)**2 *  &
              zd     = 1.02_wp*zdmax*zs1**(zbedi(i,j) + 2) *           &
                       (zsnull(i,j)/zs1)**zbf1
              zck    = (1.0_wp + 1550.0_wp*cdmin/zdmax)*zbf2
              ! maximum sustainable moisture flux in the uppermost surface layer
              ! in kg/(s*m**2)
!             zfqmax = - rho_w*zck*zd*zsnull(i,j)/SQRT(zdzwsu(nlgw+1)*zdzwg(1))
              zfqmax = - rho_w*zck*zd*zs1/SQRT(zdzwsu(nlgw+1)*zdzwg(1))
              ! consideration of frozen ground (flux reduction to 10%)
              zfqmax = zfqmax*(1.0_wp - 0.9_wp*    &
                            (1.0_wp - zts_pm(i,j)))
              zevapor= MAX(zep_s(i,j),zfqmax)
              zbeta  = zevapor/MIN(zep_s(i,j),-eps_div)
            END IF ! Computations not for ice and rocks
            zbeta  = zbeta + (1.0_wp - zbeta)*zice
            ! zbeta=1 (ice), zbeta=0 (rocks), zbeta unchanged for all other
            ! soil types
            ! consideration of plant or snow/water cover
            xesoil(i,j) = zevap*zbeta                   & ! reduction
                          * zep_s(i,j)                  & ! evaporation
                          * (1.0_wp - fr_wi  (i,j)) & ! not water covered
                          * (1.0_wp - fr_snow(i,j)) & ! not snow covered
                          * eai(i,j)/sai(i,j)
                                   ! relative source surface of the bare soil
          END IF  ! upwards directed potential evaporation
        END IF    ! land points
      END DO
    END DO
  END IF ! BATS version

  !----------------------------------------------------------------------------
  ! Section 4.3a: transpiration by plants, bucket version
  !----------------------------------------------------------------------------

  IF (itype_trvg.EQ.1) THEN   ! Bucket version
  ! for lower layers, even for water saturated soil and and maximum
  ! root extension, transpiration is limited to the potential evaporation rate
  ! the consideration of a root depth (zroot) allows the maximum
  ! transpiration, if the active soil layer is saturated and completely
  ! penetrated by roots
  
    DO kwl = 1,nlgw       ! loop over soil layers 
      DO   j = jstarts, jends
        DO i = istarts, iends
          IF (llandmask(i,j)) THEN     ! land-points only
            IF (zep_s(i,j) < 0.0_wp  &        ! potential evaporation and
               .AND. m_styp(i,j) >= 3) THEN   ! neither ice nor rock
  
              ! turgor-loss-point minus plant wilting point
              ztlpmwp = (zfcap(i,j) - zpwp(i,j)) &
                       *(0.81_wp + 0.121_wp*ATAN(-86400.0_wp*zep_s(i,j) - 4.75_wp))
              ! determine whether neither surface nor lower boundary of
              ! second layer are below freezing point
              ztgt0   = zts_pm(i,j)*zsf_heav(t_m(i,j,nx) - t0_melt)
   
              ! determine reduction factor beta for this layer
              zbeta   = ztgt0 * MAX(0.0_wp, MIN(1.0_wp,               &
                            (zwg_fr(i,j,kwl) - zpwp(i,j)) / ztlpmwp))
              zbeta   = zbeta**2
  
              ! consider the effect of root depth
              zroot = MIN ( zdzwg(kwl), MAX(0.0_wp,                       &
                         zbwt(i,j) - (zdzwsu(kwl) - 0.5_wp*zdzwg(kwl))))
              zr_root = zroot/zdzwsu(nlgw+1)
              xtrang(i,j,kwl) = zrock(i,j)*zr_root*zbeta       & ! reduction
                                   * zep_s(i,j)                & ! transpiration
                                   * (1.0_wp - fr_wi(i,j))     & ! non-water
                                   * (1.0_wp - fr_snow(i,j))   & ! non-snow
                                   * tai(i,j)/sai(i,j)       ! transp. surface
                                     ! relative source surface of the plants
            END IF  ! upwards directed potential evaporation .AND. m_styp > 2
          END IF    ! land-points only
        END DO
      END DO
    END DO          ! loop over soil layers
  END IF            ! Bucket version
  
  !----------------------------------------------------------------------------
  ! Section 4.3b: transpiration by plants, BATS version
  !----------------------------------------------------------------------------

  IF (itype_trvg.EQ.2) THEN   ! BATS version
    ! This version is based on Dickinson's (1984) BATS scheme, simplified by
    ! neglecting the water and energy transports between the soil and the plant
    ! canopy. This leads to a Monteith combination formula for the computation
    ! of plant transpiration.
    !** In the present stage, values of the lowest model level are used      ***
    !** to replace values in 2 m height and in 10 m height, respectively.    ***

    ! Determination of the transfer functions CA, CF, and CV

    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j) .AND. m_styp(i,j).ge.3) THEN      ! land points only,
                                                             ! not ice or rocks
          IF (zep_s(i,j) < 0.0_wp) THEN  ! upwards directed potential evaporation
            ! Soil water function
            zwroot(i,j)= 0.0_wp
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    DO kwl = 1,nlgw
      zr1 = zdzwsu(kwl)+0.5_wp*zdzwg(kwl)
      zr2 = zdzwsu(kwl)-0.5_wp*zdzwg(kwl)
      DO   j = jstarts, jends
        DO i = istarts, iends
          IF (llandmask(i,j) .AND. m_styp(i,j).ge.3) THEN      ! land points only,
                                                               ! not ice or rocks
            IF (zep_s(i,j) < 0.0_wp) THEN  ! upwards directed potential evaporation

              zfill    =                        zsf_heav(zbwt(i,j) - zr1)
              zpart    = (1.0_wp - zfill) * zsf_heav(zbwt(i,j) - zr2)
              zropart  = zfill*zdzwg(kwl) +     zpart * (zbwt(i,j) - zr2)
              zwroot(i,j) = zwroot(i,j) + zwg_fr(i,j,kwl)*zropart/zbwt(i,j)
!             zrootfc(i,j,kwl) =  zwg_fr(i,j,kwl)*zropart/zbwt(i,j)
              zrootfc(i,j,kwl) =  zwg_fr(i,j,kwl)*zropart/zdzwsu(nlgw+1)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j) .AND. m_styp(i,j).ge.3) THEN      ! land points only,
                                                             ! not ice or rocks
#ifndef MESSY
          IF (zep_s(i,j) < 0.0_wp) THEN  ! upwards directed potential evaporation
#endif
            zuv        = SQRT(u(i,j,ke,nx)**2 + v(i,j,ke,nx)**2)
            zcatm(i,j) = tch(i,j)*zuv           ! Function CA

            ustar_fv(i,j) = zuv*SQRT(tcm(i,j))

            IF (itype_tran.eq.1) THEN
!US               zustar     = zuv*SQRT(tcm(i,j))
               zrla       = 1.0_wp/MAX(cdash*SQRT(ustar_fv(i,j)),eps_div)
            ELSE
               zrla       = 0.0_wp
            ENDIF

            ! to compute CV, first the stomatal resistance has to be determined
            ! this requires the determination of the F-functions:
            ! Radiation function
            zpar       = pabs(i,j)  !  PAR
            zf_rad     = MAX(0.0_wp,MIN(1.0_wp,zpar/cparcrit))
            ztlpmwp    = (zfcap(i,j) - zpwp(i,j)) &
                         *(0.81_wp + 0.121_wp*ATAN(-86400.0_wp*zep_s(i,j) - 4.75_wp))

            ! zwroot now is the integral mean of fractional water content over
            ! the root depth, and zrootfc is the contribution of the single
            ! layers to this integral mean
            zf_wat     = MAX(0.0_wp,MIN(1.0_wp,(zwroot(i,j) - zpwp(i,j))/ztlpmwp))
            ! Temperature function
            zf_tem     = MAX(0.0_wp,MIN(1.0_wp,4.0_wp*(t(i,j,ke,nx) - t0_melt)*      &
                         (ctend - t(i,j,ke,nx))/(ctend - t0_melt)**2))
            ! Saturation deficit function (function not used, but computations
            ! necessary for determination of  slope of the saturation curve)
            z2iw       = zts_pm(i,j)*b2w + (1.0_wp - zts_pm(i,j))*b2i
            z4iw       = zts_pm(i,j)*b4w + (1.0_wp - zts_pm(i,j))*b4i
            zepsat     = zsf_psat_iw(t(i,j,ke,nx),z2iw,z4iw)
            zepke      = qv(i,j,ke)*ps(i,j,nx)/                   &
                                      (rdv + o_m_rdv*qv(i,j,ke))
            zf_sat     = MAX(0.0_wp,MIN(1.0_wp,1.0_wp - (zepsat - zepke)/csatdef))
            zsdef_qv(i,j) = rdv*zepsat/(ps(i,j,nx) - o_m_rdv*zepsat)   &
                         - qv(i,j,ke) ! saturation deficit qv for later use
            zf_sat     = 1.0_wp
            zedrstom   = 1.0_wp/crsmax + (1.0_wp/crsmin - 1.0_wp/crsmax)*         &
                          zf_rad*zf_wat*zf_tem*zf_sat
            zrstom     = 1.0_wp/zedrstom                 ! stomatal resistance

            rstom(i,j) = zrstom
            zrveg(i,j) = zrla+zrstom

#ifndef MESSY
          END IF  ! negative potential evaporation only
#endif
        END IF    ! land points .AND. m_styp > 2
      END DO
    END DO

    ! Determination of the transpiration rate of dry leaves
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j) .AND. m_styp(i,j).ge.3) THEN      !land points only,
                                                             ! not ice or rocks
          IF (zep_s(i,j) < 0.0_wp) THEN  ! upwards potential evaporation

            ztraleav(i,j) = zep_s(i,j)*tai(i,j) / (sai(i,j) + zrveg(i,j)*zcatm(i,j))

          END IF  ! upwards directed potential evaporation only
        END IF    ! land points .AND. m_styp > 2
      END DO
    END DO
    ! Consideration of water and snow coverage, distribution to the different
    ! soil layers
    DO     kwl       = 1,nlgw
      DO   j         = jstarts, jends
        DO i         = istarts, iends
          IF (llandmask(i,j) .AND. m_styp(i,j).ge.3) THEN ! land points only,
                                                          ! not ice or rocks
            IF (zep_s(i,j) < 0.0_wp) THEN          ! upwards potential evaporation
              ztrabpf  = ztraleav(i,j)*             & 
                         (1.0_wp - fr_wi(i,j))*        & ! not water covered
                         (1.0_wp - fr_snow(i,j))         ! not snow covered
              xtrang(i,j,kwl) = ztrabpf*zrootfc(i,j,kwl)/MAX(eps_div,zwroot(i,j))
            END IF  ! upwards directed potential evaporation only
          END IF    ! land points .AND. m_styp > 2
        END DO
      END DO
    END DO          ! loop over soil layers

  END IF ! BATS version

  !----------------------------------------------------------------------------
  ! Section 4.4: total evapotranspiration and 
  !              associated ficticious soil humidity qv_s
  !----------------------------------------------------------------------------
  
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN   ! land-points only
          ze_sum = xdwsndt(i,j  )  & ! evaporation of snow
                 + xdwidt (i,j  )  & ! evaporation from interception store
                 + xesoil (i,j  )  & ! evaporation from bare soil
                 + xtrang (i,j,1)  & ! transpiration from first soil layer
                 + xtrang (i,j,2)  & ! transpiration from second soil layer
                 + xtrang (i,j,3)  & ! transpiration from third soil layer
                 + xrr    (i,j  )  & ! formation of dew
                 + xrs    (i,j  )    ! formation of rime
          qv_s(i,j,nx  ) = qv (i,j,ke) - ze_sum /(xrhoch(i,j) + eps_div)
          qv_s(i,j,nnew) = qv_s(i,j,nx)
        END IF     ! land points
      END DO
    END DO

!------------------------------------------------------------------------------
! End of module procedure terra1
!------------------------------------------------------------------------------

END SUBROUTINE terra1 

!==============================================================================
!+ Computation of the second part of the soil parameterization scheme
!------------------------------------------------------------------------------

SUBROUTINE terra2

!------------------------------------------------------------------------------
!
! Description:
!   The module procedure terra2 performs the second part of the soil processes
!   parametrization scheme.
!   In this part, the parameterization of soil processes is completed and all
!   prognostic variables of the soil model are finally updated. The calulation
!   uses the precipitation rates from a previous call of the convection 
!   scheme and data from module arrays provided by terra1. These module arrays
!   are deallocated at the end of terra2.
!
! Method:
!   The soil model is based on the "Extendend force restore" (EFR) method for
!   the prediction of soil temperature in a two-layer scheme following Jacobsen
!   and Heise (1982, Contr.Atm.Phys.,55). The prediction of soil water is 
!   calculated by usual budget equations with a two-layer or three-layer soil 
!   moisture scheme (optional).
!
!   For a detailed description of the physical and mathematical methods,
!   see the EM/DM Documentation.
!
!------------------------------------------------------------------------------
! Subroutine arguments: None
! --------------------
!
! Local parameters:
! ----------------
  REAL    (KIND=wp)       , PARAMETER ::  &
    zimpli = 1.0_wp   ! weight for implicit calculation 
                      ! (1 means full implicit, 0 means full exlplicit)

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    nx     ,        & ! time-level index
    kwl    ,        & ! loop index for soil moisture layers           
    i      ,        & ! loop index in x-direction              
    j      ,        & ! loop index in y-direction              
    mstyp  ,        & ! soil type index
    msr_off,        & ! number of layers contributing to surface run off
    istat  ,        & ! error statur for deallocation    
    istarts,        & ! start index for x-direction      
    iends  ,        & ! end   index for x-direction     
    jstarts,        & ! start index for y-direction     
    jends             ! end   index for y-directin      

  REAL    (KIND=wp   )     ::  &
    ! various local variables for intermediate storage
    zeisb     , zgu       , zgm      , zabmbm  , zomb    , zomkon   ,  &
    zdmddb    , zomhqx    , zbbmam   , zdtdrhw , zrhwddt , zroffdt  ,  &
    z1d2dt    , zx        , zalfbs   , zalfms  , zbetss  ,             &
    zbetms    , zwgmin    , zwgmax   , zgstr   , zwinstr ,             &
    zinfmx    , zwimax    , zalf     , zkwil   , zvers   , zro_inf  ,  &
    zdwsndtt  , zwsnstr   , zdwseps  , zdwidtt , zdwieps , zro_wi   ,  &
    zro_sfak  , zro_gfak  , zfmb_fak , ztm_pm  , zdwil   , zdwg     ,  & 
    zredfu    , zro,  zwgn, zro2     , zkorr   , zrnet_s ,             &
    zshfl_s   , zlhfl_s   , zsprs    , zalas   , zgsb    , zfor_s   ,  &
    zrnet_snow, zdtmdtt   , zts_im   , zfak    , zwsneu  , ztsnow_im,  &
    ztsnew    , ztsnownew , zenm_snow, zenm_g  , zenm    , zensm    ,  &
    zfr_melt  , zdtsmx    , zdwsnm   , zdwgme  , zwinew  ,             &
    zshfl_snow, zlhfl_snow, zfor_snow, zdt     ,                       &
    zsf_heav,       & ! Heaviside function
    zstx              ! dummy argument for Heaviside function

! Local (automatic) arrays:
! -------------------------
  REAL    (KIND=wp   )     ::  &
    zadp     (ie,je)  , & ! air dryness point
    zrock    (ie,je)  , & ! ice/rock-indicator
    zporv    (ie,je)  , & ! 
    zts      (ie,je)  , & ! 
    zts_pm   (ie,je)  , & ! 
    zthsnw   (ie,je)  , & ! 
    zthsoi   (ie,je)  , & ! 
    ztsnow   (ie,je)  , & ! 
    ztsnow_pm(ie,je)  , & ! 
    zwin     (ie,je)  , & ! water cont. of interception store
    zwsnow   (ie,je)  , & ! 
    zfcap    (ie,je)  , & ! 
    zdw      (ie,je)  , & ! 
    zdzg     (ie,je)  , & ! 
    zdzm     (ie,je)  , & ! 
    zik2     (ie,je)  , & ! 
    zkw      (ie,je)  , & ! 
    ztm      (ie,je)  , & ! 
    zdw1     (ie,je)  , & ! 
    zkw1     (ie,je)  , & ! 
    zwg_fr   (ie,je,4), & ! fractional water content of layers in ground
    zdtsdt   (ie,je)  , & ! 
    zdtmdt   (ie,je)  , & ! 
    zdtsnowdt(ie,je)  , & ! 
    zroff    (ie,je)  , & ! 
    zroff_s  (ie,je)  , & ! 
    zwlarc   (ie,je)  , & ! 
    zinfil   (ie,je)  , & ! infiltration rate
    ztsn     (ie,je)  , & ! 
    ztmn     (ie,je)  , & ! 
    ztsnown  (ie,je)  , & ! 
    zwsnn    (ie,je)  , & ! 
    zverbo   (ie,je)  , & ! 
    zversn   (ie,je)  , & ! 
    zdwgdt   (ie,je,3), & ! 
    zflmg    (ie,je,4), & ! 
    zdzwsu (4)        , & ! sum over hydrological layer thicknesses
    zdzwg  (4)        , & ! 
    znlgw1f(3)            !

  INTEGER  (KIND=iintegers ) ::  &
     m_styp (ie,je)       ! soil type
!
!- End of header
!=======================================================================

!------------------------------------------------------------------------------
! Begin Subroutine terra2              
!------------------------------------------------------------------------------

! Declaration of STATEMENT-FUNCTIONS
  zsf_heav(zstx)  = 0.5_wp + SIGN( 0.5_wp, zstx)

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

! select timelevel and timestep for calculations
  IF ( l2tls ) THEN
    nx  = nnow
    zdt = dt
  ELSE
    nx  = nold
    zdt = dt2
  ENDIF

  ! Horizontal domain for computation
  istarts = istartpar
  iends   = iendpar
  jstarts = jstartpar
  jends   = jendpar

  ! Computation of derived constants
  z1d2dt = 1.0_wp/zdt      ! 1./2*timestep

  ! Number of layers contributing to surface run-off  
  msr_off  = 0
  IF (nlgw ==3) msr_off = 0
 
  zx      = SQRT(ctau1/ctau2)
  zomb    = 2.0_wp*pi/(ctau1*86400.0_wp)
  zalfbs  = 1.0_wp + zx + zx*zx
  zalfms  = zx
  zbetss  =    SQRT(1.0_wp + zx*zx)*EXP(  zx/(1.0_wp + zx))
  zbetms  = zx*SQRT(1.0_wp + zx*zx)*EXP( -zx/(1.0_wp + zx))
  zdmddb  = zalfbs/zbetms - 1.0_wp
  zomhqx  = SQRT(0.5_wp*zomb)/(1.0_wp + zx)
  zomkon  = zomhqx*zx*zx/zbetms
  zabmbm  = zomhqx*(zalfbs - zbetms)
  zbbmam  = zomhqx*(zbetss - zalfms)

  ! layer thickness for hydrological part
  IF (nlgw == 2) THEN  ! (2-layer model)
     zdzwg(1)  = cdzw12 
     zdzwg(2)  = cdzw22
  ELSE                    ! (3-layer model)
     zdzwg(1)  = cdzw13
     zdzwg(2)  = cdzw23
     zdzwg(3)  = cdzw33
  END IF
  zdzwg(nlgw+1)= zdzwg(nlgw)

  zdtdrhw = zdt/rho_w   ! 2*timestep/density of liquid water
  zrhwddt = rho_w/zdt   ! inverse of ...
  zroffdt = 0.5_wp*zdt     ! time step for run-off computation

  ! time constant for infiltration of water from interception store
  ! must not be less than 2*time step
  ctau_i  = MAX( ctau_i, zdt )

  ! Specification of depth of soil layer centers
  zdzwsu(1)       = 0.5_wp*zdzwg(1)
  DO kwl     = 1,nlgw
    znlgw1f(kwl)  = 0.0_wp
    IF (kwl < nlgw) THEN
      zdzwsu(kwl+1) = zdzwsu(kwl) + 0.5_wp*(zdzwg(kwl)+zdzwg(kwl+1))
    ENDIF
  END DO
  zdzwsu(nlgw+1)  = zdzwsu(nlgw) + 0.5_wp*zdzwg(nlgw)
  znlgw1f(1)         = 1.0_wp

!------------------------------------------------------------------------------
! Section 2: Prepare basic surface properties and create some local
!            arrays of surface related quantities (for land-points only)
!------------------------------------------------------------------------------
 
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN                 ! land-points

        zts      (i,j) = t_s(i,j,nx)
        zts_pm   (i,j) = zsf_heav(zts   (i,j) - t0_melt)
        ztm      (i,j) = t_m(i,j,nx)
        ztsnow   (i,j) = t_snow(i,j,nx)
        ztsnow_pm(i,j) = zsf_heav(ztsnow(i,j) - t0_melt)
        zwsnow   (i,j) = w_snow(i,j,nx)
        zwin     (i,j) = w_i   (i,j,nx)

        IF (nlgw == 2) THEN   ! two soil layers for hydrology
           zwg_fr(i,j,1)  = w_g1(i,j,nx)/zdzwg(1)
           zwg_fr(i,j,2)  = w_g2(i,j,nx)/zdzwg(2)
           zwg_fr(i,j,3)  = w_cl(i,j   )/zdzwg(2)
        ELSE                  ! three soil layers for hydrology
           zwg_fr(i,j,1)  = w_g1(i,j,nx)/zdzwg(1)
           zwg_fr(i,j,2)  = w_g2(i,j,nx)/zdzwg(2)
           zwg_fr(i,j,3)  = w_g3(i,j,nx)/zdzwg(3)
           zwg_fr(i,j,4)  = w_cl (i,j   )/zdzwg(3)
        END IF                ! number of soil layers for hydrology

        mstyp        = NINT(soiltyp(i,j))     ! soil type
        m_styp(i,j)  = mstyp
        zdzg  (i,j)  = cdz1  (mstyp)      ! top layer thickness
        zdzm  (i,j)  = zdzg  (i,j)*zdmddb ! second layer thickness
        zporv (i,j)  = cporv (mstyp)      ! pore volume
        zadp  (i,j)  = cadp  (mstyp)      ! air dryness point
        zfcap (i,j)  = cfcap (mstyp)      ! field capacity
        zrock (i,j)  = crock (mstyp)      ! rock or ice indicator
        zdw   (i,j)  = cdw0  (mstyp)      ! hydrological diff.parameter
        zdw1  (i,j)  = cdw1  (mstyp)      ! hydrological diff.parameter
        zkw   (i,j)  = ckw0  (mstyp)      ! hydrological cond.parameter
        zkw1  (i,j)  = ckw1  (mstyp)      ! hydrological cond.parameter
        zik2  (i,j)  = cik2  (mstyp)      ! minimum infiltration rate

      END IF
    END DO
  END DO

!------------------------------------------------------------------------------
! Section 3: Estimate thermal surface fluxes
!------------------------------------------------------------------------------
 
  ! Estimate thermal surface fluxes over snow covered and snow free
  ! part of surface based on area mean values calculated in radiation
  ! code (positive = downward)

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN
        zgstr = - sigma*(1.0_wp - Ctalb) * ( (1.0_wp - fr_snow(i,j))*zts(i,j) &
                  + fr_snow(i,j)*ztsnow(i,j) )**4 - thbs(i,j)
        zthsnw(i,j) = - sigma*(1.0_wp - Ctalb)*ztsnow(i,j)**4 - zgstr
        zthsoi(i,j) = - sigma*(1.0_wp - Ctalb)*zts(i,j)**4 - zgstr
        ! the estimation of the solar component would require the availability
        ! of the diffuse and direct components of the solar flux
      END IF
    END DO
  END DO

 
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF(llandmask(i,j))THEN     ! land-points only

        ! store forcing terms due to evapotranspiration, formation of dew
        ! and rime for later use
        zverbo(i,j) = xdwidt(i,j) + xesoil(i,j) + xtrang(i,j,1) + xtrang(i,j,2)&
                      + xtrang(i,j,3) + (1.0_wp-fr_snow(i,j))*(xrr(i,j) + xrs(i,j))
        zversn(i,j) = xdwsndt(i,j) + xrs(i,j)                                  &
                                   + zsf_heav (zwsnow(i,j) - eps_soil) * xrr(i,j)

        ! add grid scale and convective precipitation (and graupel, if present)
        ! to dew and rime
        xrr(i,j) = xrr(i,j) + prr_con(i,j) + prr_gsp(i,j)
        IF ( itype_gscp >= 2000 ) THEN
          xrs(i,j) = xrs(i,j) + prs_con(i,j) + prs_gsp(i,j) + prg_gsp(i,j) + prh_gsp(i,j)
        ELSEIF ( itype_gscp >= 4 ) THEN
          xrs(i,j) = xrs(i,j) + prs_con(i,j) + prs_gsp(i,j) + prg_gsp(i,j)
        ELSE
          xrs(i,j) = xrs(i,j) + prs_con(i,j) + prs_gsp(i,j)
        ENDIF

        ! infiltration and surface run-off
 
        ! subtract evaporation from interception store to avoid negative
        ! values due to sum of evaporation+infiltration
        zwinstr = zwin(i,j) + xdwidt(i,j)*zdtdrhw
        zwinstr = MAX(0.0_wp,zwinstr)

        ! maximum infiltration rate (rock/ich/water-exclusion
        zinfmx = zrock(i,j)*zts_pm(i,j)*csvoro &
                 *( cik1*MAX(0.5_wp,plcov(i,j))*&
             MAX(0.0_wp, zporv(i,j)-zwg_fr(i,j,1))/zporv(i,j) + zik2(i,j) )      
        zwimax = cwimax*(1.0_wp + plcov(i,j)*5.0_wp)
        zalf   = SQRT(MAX(0.0_wp,1.0_wp - zwinstr/zwimax))

        ! infiltration of water from interception store (if Ts above freezing)
        zvers           = zts_pm(i,j)*zwinstr*rho_w/ctau_i
      
        ! possible contribution of rain to infiltration
        IF (xrr(i,j)-eps_soil > 0.0_wp) THEN
          zalf = MAX( zalf, (zrhwddt*MAX(0.0_wp, zwimax-zwinstr) + zvers)/xrr(i,j) )
          zalf = MAX( 0.01_wp, MIN(1.0_wp, zalf) )
          zalf = zts_pm(i,j)*zalf + (1.0_wp - zts_pm(i,j))

          ! if rain falls onto snow, all rain is considered for infiltration
          ! ?  wieso, der regen kann doch auch fallen, wenn schnee 
          ! ?  nur partiell liegt 
          IF (zwsnow(i,j) > 0.0_wp) zalf = 0.0_wp
        END IF

        ! add rain contribution to infiltration rate
        zvers = zvers + (1.0_wp - zalf)*xrr(i,j)

        ! final infiltration rate limited by maximum value
        zinfil(i,j) = MIN(zinfmx,zvers)

        ! surface run-off (residual of potential minus actual infiltration)
        zro_inf       = MAX( 0.0_wp, zvers - zinfil(i,j) )
        runoff_s(i,j) = runoff_s(i,j) + zro_inf*zroffdt
        zflmg (i,j,1) = - zinfil(i,j)

        ! change of snow water and interception water store
        ! (negligible residuals are added to the run-off)

        ! snow store
        zdwsndtt = xrs(i,j) + xdwsndt(i,j)
        zwsnstr  = zwsnow(i,j) + zdwsndtt*zdtdrhw
        zwsnstr  = MAX(0.0_wp, zwsnstr)    ! avoid negative values
        zdwseps  = 0.0_wp
        IF (zwsnstr > 0.0_wp .AND. zwsnstr < eps_soil) THEN
          zdwseps    = zwsnstr*zrhwddt
          runoff_s(i,j) = runoff_s(i,j) + zdwseps*zroffdt
          zdwsndtt   = - zsf_heav(zwsnow(i,j) - 0.1_wp*eps_soil) &
                         *zwsnow(i,j)*zrhwddt
          zwsnstr    = 0.0_wp   ! ? what for
        END IF
        xdwsndt(i,j) = zdwsndtt
        zdwidtt      = xrr(i,j) + xdwidt(i,j)-zinfil(i,j)-zro_inf
        zwinstr      = zwin(i,j) + zdwidtt*zdtdrhw
        zwinstr      = MAX( 0.0_wp, zwinstr )
 
        ! interception store
        zdwieps      = 0.0_wp
        IF (zwinstr > 0.0_wp .AND. zwinstr < eps_soil) THEN
          zdwieps    = zwinstr*zrhwddt
          runoff_s(i,j)= runoff_s(i,j) + zdwieps*zroffdt
          zdwidtt    = - zsf_heav(zwin(i,j) - 0.1_wp*eps_soil) &
                         *zwin(i,j)*zrhwddt
          zwinstr    = 0.0_wp
        END IF
        ! add excess over zwimax to runoff  
        zro_wi       = zts_pm(i,j)*zrhwddt*MAX( 0.0_wp, zwinstr-zwimax )
        zdwidtt      = zdwidtt - zro_wi
        xdwidt(i,j)  = zdwidtt
        runoff_s(i,j)= runoff_s(i,j) + zro_wi*zroffdt
        zroff_s(i,j) = 0.0_wp
        zroff  (i,j) = 0.0_wp
      END IF            ! land-points only
    END DO
  END DO

!------------------------------------------------------------------------------
! Section 4: Loop over soil layers
!------------------------------------------------------------------------------
 
  DO  kwl = 1,nlgw

    ! utility variables used to avoid if-constructs in following loops
    zro_sfak = zsf_heav(0.5_wp + msr_off - kwl)  
    zro_gfak = 1.0_wp - zro_sfak
    zfmb_fak = zsf_heav(0.5_wp + kwl - nlgw)

    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN      ! land-points only
          zdwgdt(i,j,kwl)   = 0.0_wp
          zflmg (i,j,kwl+1) = 0.0_wp
          ztm_pm            = zsf_heav(ztm(i,j) - t0_melt)

          ! sedimentation and capillary transport in soil

          IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
            ! minimum and maximum fractional water content of actual layer and
            ! the next layer below it
            zwgmin = MIN( zwg_fr(i,j,kwl), zwg_fr(i,j,kwl+1) )
            zwgmax = MAX( zwg_fr(i,j,kwl), zwg_fr(i,j,kwl+1) )
            zdwil  = zdw(i,j) * EXP( zdw1(i,j) * (zporv(i,j)                 &
                                       - cakw * zwgmin - (1.0_wp-cakw) * zwgmax) &
                                       / (zporv(i,j) - zadp(i,j)) )
            zkwil  = zkw(i,j) * EXP( zkw1(i,j) * (zporv(i,j)                 &
                                       - cakw * zwgmin - (1.0_wp-cakw) * zwgmax) &
                                       / (zporv(i,j) - zadp(i,j)) )
            zflmg(i,j,kwl+1) = zrock(i,j)*ztm_pm*rho_w   &
                               *( zdwil*(zwg_fr(i,j,kwl+1) - zwg_fr(i,j,kwl)) &
                                  /(0.5_wp*(zdzwg(kwl)+zdzwg(kwl+1))) - zkwil ) 
            ! limit flux at lower boundary of first layer to infiltration rate
            IF (kwl == 1 .AND. zflmg(i,j,2) < 0.0_wp .AND. zinfil(i,j) > 0.0_wp) THEN
              zflmg(i,j,2) = MAX( zflmg(i,j,2), -zinfil(i,j) )
            ENDIF
            ! 10% of pore volume is allowed as maximum downward flux 
            ! per timestep
            zflmg(i,j,kwl+1) = MAX ( zflmg(i,j,kwl+1),                         &
                                     -0.1_wp * zrhwddt * zdzwg(kwl) * zporv(i,j) )
            ! Artificial reduction of upward flux at the lower boundary of first
            ! layer to 50 % (to account for the normally different structure 
            ! of this layer).
            IF (kwl.eq.1 .AND. zflmg(i,j,2) .gt.0.0_wp)                           &
                zflmg(i,j,2) = 0.5_wp*zflmg(i,j,2)

            ! set flux at lower boundary of soil model to zero
            zflmg(i,j,nlgw+1) = 0.0_wp

            ! first run_off calculation without consideration of 
            ! evapotranspiration    
            zdwg   = zflmg(i,j,kwl+1) - zflmg(i,j,kwl)
            zredfu =  MAX( 0.0_wp, MIN( 1.0_wp,            &
              (zwg_fr(i,j,kwl)-zfcap(i,j))/MAX(zporv(i,j)-zfcap(i,j),eps_div)) )
            zredfu = zsf_heav(zdwg)*zredfu
            zro    = zdwg*zredfu
            zdwg   = zdwg*(1.0_wp - zredfu)

            ! add evaporation (znlgw1f: first layer only) 
            ! and transpiration (for each layer) 
            zdwg   = zdwg + znlgw1f(kwl) * xesoil(i,j) * zrock(i,j)           &
                          + xtrang (i,j,kwl)
            zwgn   = zwg_fr(i,j,kwl) + zdtdrhw*zdwg/zdzwg(kwl)
            zro2   = zrhwddt*zdzwg(kwl)*MAX(0.0_wp, zwgn - zporv(i,j))*zrock(i,j)
            zkorr  = zrhwddt*zdzwg(kwl)*MAX(0.0_wp, zadp(i,j) - zwgn )*zrock(i,j)
            zdwgdt(i,j,kwl)= zdwg + zkorr - zro2
            zro    = zro      + zro2
            zflmg(i,j,kwl+1)= zflmg(i,j,kwl+1) + zkorr
            zroff_s (i,j) = zroff_s (i,j) + zro*zro_sfak
            runoff_s(i,j) = runoff_s(i,j) + zro*zro_sfak*zroffdt
            zroff   (i,j) = zroff   (i,j) + zro*zro_gfak
            runoff_g(i,j) = runoff_g(i,j) + zro*zro_gfak*zroffdt
            runoff_g(i,j) = runoff_g(i,j) - zflmg(i,j,nlgw+1)*zroffdt*zfmb_fak
          END IF          ! ice/rock-exclusion

        END IF   ! land-points only
      END DO  
    END DO  
  END DO         ! end loop over soil layers           

!------------------------------------------------------------------------------
! Section 5: Thermal part
!------------------------------------------------------------------------------

  ! freezing of soil water, heat fluxes and temperature tendency 
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN  ! land-points only
        zwlarc(i,j) = SQRT(xalam(i,j)*xrocg(i,j))
        ! multiplication of xrocg and xrocm by the layer thickness, after
        ! zwlarc has been multiplied by the proper value vor the volume 
        ! heat content
        xrocg(i,j)  = xrocg(i,j)*zdzg(i,j)
        xrocm(i,j)  = xrocm(i,j)*zdzm(i,j)
      END IF
    END DO  
  END DO  

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only

        ! heat fluxes 
        zgu   = zwlarc(i,j)* zomkon*(t_cl(i,j)-ztm(i,j))
        zgm   = zwlarc(i,j)*(zabmbm*(t_cl(i,j)-zts(i,j))                        &
                                            - zbbmam*(t_cl(i,j)-ztm(i,j)))

        ! forcing for snow-free soil
        ! (evaporation, transpiration, formation of dew and rime are already
        !  weighted by correspondind surface fraction)
        ! net radiation, sensible and latent heat flux
        zrnet_s = (1.0_wp - fr_snow(i,j))*(sobs(i,j)+zthsoi(i,j)) 
        zshfl_s = (1.0_wp - fr_snow(i,j))*cp_d*xrhoch(i,j)*(xth_low(i,j) - zts(i,j))
        zlhfl_s = (zts_pm(i,j)*lh_v + (1.0_wp-zts_pm(i,j))*lh_s)*zverbo(i,j)

        zsprs      = 0.0_wp
        ! thawing of snow falling on soil with Ts > T0
        IF (ztsnow_pm(i,j)*xrs(i,j) > 0.0_wp) THEN
          ! snow fall on soil with T>T0, snow water content increases 
          ! interception store water content
          zsprs        = - lh_f*xrs(i,j)
          xdwidt (i,j) = xdwidt (i,j) + xrs(i,j)
          xdwsndt(i,j) = xdwsndt(i,j) - xrs(i,j)

          ! avoid overflow of interception store, add possible excess to
          ! water content of first layer from where it may run-off in next 
          ! time-step
          zwimax       = cwimax*(1.0_wp + plcov(i,j)*5.0_wp)
          zwinstr      = zwin(i,j) + xdwidt(i,j)*zdtdrhw
          IF (zwinstr > zwimax) THEN  ! overflow of interception store
            zvers      = (zwinstr - zwimax)*zrhwddt
            zinfil(i,j)= zinfil(i,j) + zvers
            xdwidt(i,j)= xdwidt(i,j) - zvers
            zdwgdt(i,j,1)= zdwgdt(i,j,1) + zvers
          ENDIF                       ! overflow of interception store

        ! freezing of rain falling on soil with Ts < T0  (black-ice !!!)
        ELSEIF (zwsnow(i,j) == 0.0_wp .AND. (1.0_wp-ztsnow_pm(i,j))*xrr(i,j) > 0.0_wp)THEN
          zsprs   = lh_f*xrr(i,j)
          xdwidt (i,j) = xdwidt (i,j) - xrr(i,j)
          xdwsndt(i,j) = xdwsndt(i,j) + xrr(i,j)
        END IF

        ! heat conductivity of snow as funtion of water content
        zalas  = MAX(calasmin,MIN(calasmax, calasmin + calas_dw*zwsnow(i,j)))
        zgsb   = zalas*(ztsnow(i,j) - zts(i,j))/xdzs(i,j)

        ! total forcing for uppermost soil layer
        zfor_s = zrnet_s + zshfl_s + zlhfl_s + zsprs                          &
                         + fr_snow(i,j) * (1.0_wp-ztsnow_pm(i,j)) * zgsb

        ! forcing contributions for snow formation of dew and rime are 
        ! contained in ze_ges, heat fluxes must not be multiplied by 
        ! snow covered fraction
        zrnet_snow = sobs(i,j) + zthsnw(i,j)
        zshfl_snow = xrhoch(i,j)*cp_d*(xth_low(i,j) - ztsnow(i,j))
        zlhfl_snow = lh_s*zversn(i,j)
        zfor_snow  = zrnet_snow + zshfl_snow + zlhfl_snow

        ! temperature forecast
        ! (in the following two sections a correct formulation for the
        !  impact of the change in water content is still missing) 

        zdtmdtt   = zdt*2.0_wp*(zgu-zgm)/xrocm(i,j)
        ztmn(i,j) = ztm(i,j) + zdtmdtt                        ! Tm(t+2*dt)
        ztsn(i,j) = zts(i,j) + zdt*2.0_wp*(zgm+zfor_s)/xrocg(i,j)-zdtmdtt 
                                                              ! Ts(t+2*dt)

        ! implicit treatment of Ts-forecast 
        ! (sensible and latent heat flux, zgm and zgsb (for snow temperature) 
        ! are considered)
        zts_im = - (1.0_wp - fr_snow(i,j)) * xrhoch(i,j)                          &
                   * (cp_d + xdqvts(i,j) * (zts_pm(i,j) * lh_f                    &
                          + (1.0_wp - zts_pm(i,j)) * lh_s))                       &
                   - zwlarc (i,j) * zabmbm - fr_snow(i,j) * zalas/xdzs(i,j)
        zfak   = MAX (eps_soil, 1.0_wp - zdt * zimpli * zts_im / xrocg(i,j))
        ztsn   (i,j) = zts (i,j) + (ztsn(i,j)-zts(i,j))/zfak
        ztsnown(i,j) = ztsn(i,j)
        zwsneu       = zwsnow(i,j) + xdwsndt(i,j)*zdtdrhw

        ! forecast of snow temperature Tsnow
        IF (ztsnow(i,j) < t0_melt .AND. zwsneu > eps_soil) THEN
          ztsnown(i,j) = ztsnow(i,j) + zdt*2.0_wp*(zfor_snow - zgsb)/xrocs(i,j)  &
                                    - ( ztsn(i,j) - zts(i,j) )
          ! implicit formulation  
          ztsnow_im    = - xrhoch(i,j) * (cp_d + xdqvtsnow(i,j) * lh_s)       &
                                       - zalas/xdzs(i,j)
          zfak  = MAX(eps_soil,1.0_wp - zdt*zimpli*ztsnow_im/xrocs(i,j))
          ztsnown(i,j) = ztsnow(i,j) + (ztsnown(i,j)-ztsnow(i,j))/zfak   
        ELSE
          ztsnown(i,j) = ztsn(i,j)
        END IF

        zdtsnowdt(i,j) = (ztsnown(i,j) - ztsnow(i,j))*z1d2dt
        zdtsdt   (i,j) = (ztsn   (i,j) - zts   (i,j))*z1d2dt
        zdtmdt   (i,j) = (ztmn   (i,j) - ztm   (i,j))*z1d2dt
        zwsnn    (i,j) =  zwsnow (i,j) + zdtdrhw*xdwsndt(i,j)

      END IF          ! land-points only
    END DO  
  END DO  

!------------------------------------------------------------------------------
! Section 6: Melting of snow
!------------------------------------------------------------------------------

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only

        zeisb         = zsf_heav(m_styp(i,j) - 1.5_wp)
        zwsneu        = zwsnow(i,j) + zdtdrhw*xdwsndt(i,j)
        zroff_s(i,j)  = 0.0_wp
        zroff  (i,j)  = 0.0_wp

        IF (ztsnow(i,j) < t0_melt .AND. zwsneu > 0.0_wp) THEN
          ztsnew    = ztsn(i,j)         ! initial
          ztsnownew = ztsnown(i,j)      ! settings

          ! a)soil surface temperature below freezing, but snow temperature
          !   above freezing leads to melting at the top; melted water
          !   penetrates the snow layer, where it freezes again, i.e. the
          !   heat content is only redistributed in the snow layer
          IF (ztsnown(i,j) > t0_melt) THEN   ! new snow temperature > t0_melt
            ztsnownew = t0_melt - eps_soil   ! increase of soil temperature
            ztsnew    = ztsn(i,j) + zeisb*(ztsnown(i,j) - ztsnownew)/ &
                                    (1.0_wp+xrocg(i,j)/xrocs(i,j))
            zdtsnowdt(i,j) = zdtsnowdt(i,j) + (ztsnownew-ztsnown(i,j))*z1d2dt
            zdtsdt   (i,j) = zdtsdt   (i,j) + (ztsnew   -   ztsn(i,j))*z1d2dt
          END IF

          ! b)melting at the bottom of the snow layer
          IF (ztsnew > t0_melt) THEN
            ! energy available from change of snow mean layer temperature to
            ! t0_melt - epsilon
            zenm_snow = xrocs(i,j)*( 0.5_wp*(ztsnownew+ztsnew) - (t0_melt-eps_soil) )
            ! energy available in uppermost soil layer for melting
            zenm_g    = 0.5_wp*xrocg(i,j)*(ztsnew - (t0_melt + eps_soil))
            ! total available energy for melting of snow     
            zenm      = zenm_g + zenm_snow

            ! redistribution of heat only, if zenm < 0
            IF (zenm <= 0.0_wp) THEN
              zdtsnowdt(i,j) = zdtsnowdt(i,j) + z1d2dt*(ztsnew - (t0_melt - eps_soil)) &
                                              * (1.0_wp + xrocg(i,j)/xrocs(i,j))
              zdtsdt   (i,j) = zdtsdt   (i,j) - z1d2dt*(ztsnew - (t0_melt - eps_soil)) &
                                              * zeisb
            ELSE    ! melting
              zensm    = lh_f*zwsnn(i,j)*rho_w  ! latent energy content of snow
              zfr_melt = MIN(zensm,zenm)/zensm  ! fraction of snow to be melted
              IF (zfr_melt > 0.9999_wp) zfr_melt = 1.0_wp

              ! reduction of Ts due to melting; Tsnow follows from new Ts
              ! and the mean snow layer temperature t0_melt - epsilon
              zdtsmx = (ztsnew - (t0_melt + eps_soil))*z1d2dt
              zdtsdt(i,j) = zdtsdt(i,j) - zeisb * MIN (zdtsmx,               &
                                    2.0_wp* (zensm-zenm_snow) * z1d2dt/xrocg(i,j))
              ztsnew      = zts(i,j) + zdt*zdtsdt(i,j)
              ztsnownew   = 2.0_wp*(t0_melt - eps_soil) - ztsnew

              ! melted snow is allowed to penetrate the soil (up to field 
              ! capacity), if the soil type is neither ice nor rock; else it 
              ! contributes to surface run-off;
              ! fractional water content of the first soil layer determines
              ! a reduction factor which controls additional run-off
              zdwsnm       = zfr_melt*zwsneu*z1d2dt*rho_w  ! available water
              snow_melt(i,j) = snow_melt(i,j) + zdwsnm/z1d2dt

              zdwgme       = zdwsnm*zrock(i,j)            ! contribution to wg
              zroff_s(i,j) = (1.0_wp - zrock(i,j))*zdwsnm     ! surface run_off
              zredfu       = MAX( 0.0_wp,  MIN( 1.0_wp, &
               (zwg_fr(i,j,1)-zfcap(i,j))/ MAX(zporv(i,j)-zfcap(i,j), eps_div)))
              zroff(i,j)   = zdwgme*zredfu             ! excess in soil
              zinfil(i,j)  = zinfil(i,j  ) + zdwgme*(1.0_wp - zredfu)
              zdwgdt(i,j,1)= zdwgdt(i,j,1) + zdwgme*(1.0_wp - zredfu)
              xdwsndt(i,j) = xdwsndt (i,j) - zdwsnm
              IF (zfr_melt > 0.9999_wp) xdwsndt(i,j)= -zwsnow(i,j)*zrhwddt

              ! if all snow melts, Tsnow must be set equal to Ts
              ztsnew = zts(i,j) + zdt*zdtsdt(i,j)
              zfak   = zsf_heav(zenm - zensm)    ! =1, if all snow melts
              ztsnownew      = (1.0_wp - zfak ) * ztsnownew                       &
                                     + zfak * MAX( t0_melt+eps_soil, ztsnew )
              zdtsnowdt(i,j) = (ztsnownew - ztsnow(i,j))*z1d2dt
              runoff_s (i,j) = runoff_s(i,j)+(zroff_s(i,j)+zroff(i,j))*zroffdt

            END IF     ! redistribution of heat or true net melting
          END IF       ! melting from bottom of snow layer
        END IF

        ! freezing of water in interception store is not accounted for,
        ! instead the corresponding water content is added to the run-off
        zwinew          = zwin(i,j) + zdtdrhw*xdwidt(i,j)
        IF (ztsnow(i,j) > t0_melt .AND. ztsnown(i,j) < t0_melt .AND. zwinew > 0.0_wp) THEN
          zroff_s (i,j) = zwinew*zrhwddt
          runoff_s(i,j) = runoff_s(i,j)+zroff_s(i,j)*zroffdt
          xdwidt (i,j) = xdwidt(i,j) - zroff_s(i,j)
        END IF

      ENDIF          ! land-points only
    ENDDO  
  ENDDO  

!------------------------------------------------------------------------------
! Section 7: Final updating of prognostic values 
!------------------------------------------------------------------------------

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN  ! for landpoints only

        t_snow(i,j,nnew) = t_snow(i,j,nx) + zdt*zdtsnowdt(i,j)
        t_s   (i,j,nnew) = t_s   (i,j,nx) + zdt*zdtsdt   (i,j)
        t_m   (i,j,nnew) = t_m   (i,j,nx) + zdt*zdtmdt   (i,j)
        w_snow (i,j,nnew) = w_snow (i,j,nx) + zdt*xdwsndt(i,j) / rho_w  
        w_i (i,j,nnew) = w_i (i,j,nx) + zdt*xdwidt (i,j)  /rho_w  
        w_g1(i,j,nnew) = w_g1(i,j,nx) + zdt*zdwgdt (i,j,1)/rho_w  
        w_g2(i,j,nnew) = w_g2(i,j,nx) + zdt*zdwgdt (i,j,2)/rho_w  

        IF (nlgw == 3) THEN
          w_g3(i,j,nnew) = w_g3(i,j,nx) + zdt*zdwgdt(i,j,3)/rho_w  
        END IF
        IF (w_snow(i,j,nnew) <= eps_soil) w_snow(i,j,nnew) = 0.0_wp
        IF (w_i(i,j,nnew) <= eps_soil) w_i(i,j,nnew) = 0.0_wp
        IF (w_snow(i,j,nnew) >  0.0_wp .AND. t_snow(i,j,nnew) > t0_melt) THEN
          w_i   (i,j,nnew) = w_snow(i,j,nnew)
          w_snow(i,j,nnew) = 0.0_wp
        END IF

      END IF          ! land-points only
    END DO  
  END DO  

! computation of the temperature at the boundary soil/snow-atmosphere

  CALL tgcom ( t_g(:,:,nnew), t_snow(:,:,nnew), t_s(:,:,nnew),      &
               w_snow(:,:,nnew), llandmask(:,:), ie, je, cf_snow,   &
               istarts, iends, jstarts, jends )

!------------------------------------------------------------------------------
! Section 7: Dealloate module arrays
!------------------------------------------------------------------------------

  DEALLOCATE ( xdzs      , STAT = istat )
  DEALLOCATE ( xalam     , STAT = istat )
  DEALLOCATE ( xdqvts    , STAT = istat )
  DEALLOCATE ( xdqvtsnow , STAT = istat ) 
  DEALLOCATE ( xdwidt    , STAT = istat )
  DEALLOCATE ( xdwsndt   , STAT = istat )
  DEALLOCATE ( xesoil    , STAT = istat )
  DEALLOCATE ( xrr       , STAT = istat )
  DEALLOCATE ( xrs       , STAT = istat )
  DEALLOCATE ( xrhoch    , STAT = istat )
  DEALLOCATE ( xrocg     , STAT = istat )
  DEALLOCATE ( xrocm     , STAT = istat )
  DEALLOCATE ( xrocs     , STAT = istat )
  DEALLOCATE ( xth_low   , STAT = istat )
  DEALLOCATE ( xtrang    , STAT = istat )

!------------------------------------------------------------------------------
! End of module procedure terra2
!------------------------------------------------------------------------------

END SUBROUTINE terra2 

!==============================================================================

!------------------------------------------------------------------------------
! End of module src_soil
!------------------------------------------------------------------------------

END MODULE src_soil     
