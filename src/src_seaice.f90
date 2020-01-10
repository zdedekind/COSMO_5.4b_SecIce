!+ Source module for the sea ice model
!-------------------------------------------------------------------------------

MODULE src_seaice

!-------------------------------------------------------------------------------
!
! Description:
!   The module "src_seaice" performs calculations related to the
!   parameterization of sea ice. It contains the subroutine seaice.f90 with the
!   present sea ice model which is based on work by D. Mironov and B. Ritter
!   (DWD).
!
!   All global variables of the model that are used by the sea ice routine are
!   imported by USE statements below.
!
!   The parameterization package has been extracted from the DWD global model
!   GME (V2_20 of the DWD gmtri library). Some modifications have been made in
!   the code: Internal communication by common-blocks is replaced by module
!   parameters, scalars and arrays defined in this module.
!
! Current Code Owner: DWD, Jan-Peter Schulz
!  phone:  +49  69  8062 2751
!  fax:    +49  69  8062 3721
!  email:  jan-peter.schulz@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V4_10        2009/09/11 Jan-Peter Schulz
!  Initial release
! V4_12        2010/05/11 Ulrich Schaettler, Jan-Peter Schulz
!  Renamed t0 to t0_melt because of conflicting names (Uli)
!  Introduced call to SR tgcom (JPS)
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer
!  Replaced qx-variables by using them from the tracer module
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  Introduced MESSy interface
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
! V5_4a        2016-05-10 Matthias Raschendorfer
!  Use new module turb_data instead of data_turbulence
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
    wp         ,    & ! KIND-type parameter for real variables
    iintegers         ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------
    ! number of grid points for this domain
    ie         ,    & ! number of grid points in zonal direction
    je         ,    & ! number of grid points in meridional direction
    ke         ,    & ! number of grid points in vertical direction

! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from
!    the other ones because of the use of the staggered Arakawa-C-grid.
!
    istartpar  ,    & ! start index for computations in the parallel program
    iendpar    ,    & ! end index for computations in the parallel program
    jstartpar  ,    & ! start index for computations in the parallel program
    jendpar    ,    & ! end index for computations in the parallel program

! 5. variables for the time discretization and related variables
! --------------------------------------------------------------
    dt         ,    & ! long time-step

! 8. Organizational variables to handle the COSMO humidity tracers
! ----------------------------------------------------------------
    idt_qv

! end of data_modelconfig

!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 2. physical constants and related variables
! -------------------------------------------
    t0_melt    ,    & ! melting temperature of ice
    r_d        ,    & ! gas constant for dry air
    rdv        ,    & ! r_d / r_v
    o_m_rdv    ,    & ! 1 - r_d/r_v
    rvd_m_o    ,    & ! r_v/r_d - 1
    cp_d       ,    & ! specific heat of dry air at constant pressure
    rdocp      ,    & ! r_d / cp_d
    lh_f       ,    & ! latent heat of fusion
    lh_s       ,    & ! latent heat of sublimation
    rho_ice    ,    & ! density of ice          (kg/m^3)

! 3. constants for parametrizations
! ---------------------------------
    b1         ,    & ! variables for computing the saturation vapour pressure
    b2i        ,    & ! over water (w) and ice (i) according to Teten's formula
    b3         ,    & !               -- " --
    b4i               !               -- " --

! end of data_constants

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------
    p0         ,    & ! reference pressure at full levels             ( Pa  )

! 2. external parameter fields                                        (unit)
! ----------------------------
    depth_lk   ,    & ! lake depth                                    (  m  )
    llandmask  ,    & ! landpoint mask
    lseamask   ,    & ! ocean point mask, i.e. water but not lake

! 3. prognostic variables                                             (unit)
! -----------------------
    u          ,    & ! zonal wind speed                              ( m/s )
    v          ,    & ! meridional wind speed                         ( m/s )
    t          ,    & ! temperature                                   (  k  )
    pp         ,    & ! deviation from the reference pressure         ( pa  )

! 5. fields for surface values and soil/canopy model variables        (unit )
! ------------------------------------------------------------
    ps         ,    & ! surface pressure                              ( pa  )
    t_snow     ,    & ! temperature of the snow-surface               (  K  )
    t_snow_mult,    & ! temperature of the snow-surface               (  K  )
    t_s        ,    & ! temperature of the ground surface             (  K  )
    t_g        ,    & ! weighted surface temperature                  (  K  )
    w_snow     ,    & ! water content of snow                         (m H2O)

!   fields for prognostic variables of the lake model FLake or ocean
!   variables
    t_ice      ,    & ! temperature of ice/water surface              (  K  )
    h_ice      ,    & ! lake/sea ice thickness                        (  m  )

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------
!   turbulent coefficients at the surface
    tch        ,    & ! transfer coefficient for heat and moisture    ( -- )
    tfv        ,    & ! laminar reduction factor for evaporation      ( -- )

!   fields from the radiation scheme
    sobs       ,    & ! solar radiation at the ground                 ( W/m2)
    thbs              ! thermal radiation at the ground               ( W/m2)

! end of data_fields

!------------------------------------------------------------------------------

USE turb_data       , ONLY :   &
    zt_ice            ! freezing temperature of sea ice          (-1.7 deg C)

!------------------------------------------------------------------------------

USE data_soil       , ONLY :   &
    cf_snow           ! parameter for the calculation of the
                      ! fractional snow coverage

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    nstart     ,    & ! first time step of the forecast
    ntstep     ,    & ! actual time step
                      ! indices for permutation of three time levels
    nnow       ,    & ! corresponds to ntstep
    nnew       ,    & ! corresponds to ntstep + 1

! 3. controlling the physics
! --------------------------
    lmulti_snow,    & ! run multi-layer snow model
    llake      ,    & ! forecast with lake model

! 4. controlling the dynamics
! ---------------------------
    l2tls             ! time integration by two timelevel RK-scheme (.TRUE.)
                      ! or by default three-time level KW-scheme (.FALSE.)

! end of data_runcontrol

!------------------------------------------------------------------------------

USE data_parallel   , ONLY : my_cart_id

!------------------------------------------------------------------------------
! Imported utilities
!------------------------------------------------------------------------------

USE meteo_utilities , ONLY : tgcom

!------------------------------------------------------------------------------

USE environment     , ONLY : model_abort

!------------------------------------------------------------------------------

USE src_tracer      , ONLY : trcr_get, trcr_errorstr

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Declarations

!==============================================================================
! Module procedure in "src_seaice"
!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure in "src_seaice"
!------------------------------------------------------------------------------

SUBROUTINE seaice

!------------------------------------------------------------------------------
!
! Description: Parameterisation of the evolution of sea ice thickness and
!              temperature.
!
! Method:      The scheme uses the atmospheric forcing (fluxes of latent and
!              sensible heat, solar and thermal radiation) and the internal
!              heat flux through the ice sheet to simulate the evolution of
!              ice thickness and temperature at the top of the ice layer.
!              Depending on the sign of the atmospheric flux and the
!              temperature difference between the top and the bottom of the
!              ice layer either melting or freezing are simulated. The
!              temperature at the top of the ice layer is adjusted accordingly.
!              The temperature at the bottom of the ice layer is constant and
!              set to the freezing temperature assumed for (sea) water.
!
!------------------------------------------------------------------------------
!
! Declarations:
!
! Local parameters:
! ----------------

REAL    (KIND=wp   ),     PARAMETER ::  &
  zepsi    = 1.0E-6_wp,   & ! security constant
  zH_imin  = 0.05_wp  ,   & ! minimum ice thickness        [m]
  zc_i     = 2100.0_wp,   & ! ice heat capacity            [J/(kg K)]
  zk_i     = 1.2E-06_wp     ! ice temperature conductivity [m**2/s]   

! Local scalars:
! -------------

INTEGER (KIND=iintegers) ::  &

! Indices
  izerror                , & ! error index
  nx                     , & ! time level for integration
  i                      , & ! loop index in x-direction
  j                      , & ! loop index in y-direction
  im1, jm1                   ! i-1, j-1

REAL    (KIND=wp   )     ::  &
  zdt                    , & ! integration time step

! Connection to the atmosphere
  zuv                    , & ! wind velocity in lowest atmospheric layer
  zpatm                  , & ! pressure in lowest atmospheric layer
  ztatm                  , & ! potential temperature in lowest atmospheric layer
  zqsat                  , & ! saturation specific humidity at the surface
  ztvs                   , & ! virtual temperature at the surface
  zrho                       ! air density at the surface

REAL    (KIND=wp   )     ::  &
  zrhoc_i                , & ! derived constants
  zcidl,zrhoiLf,zkicidLf , & ! derived constants
  zthfl_a,zhfl_i         , & ! atmospheric and ice internal heat flux [W/m**2]
  zshfl,zlhfl            , & ! sensible and latent heat flux          [W/m**2]
  zq_i                   , & ! ice internal temperature flux          [K m/s]
  zT_i,zT_im1            , & ! ice surface temperature                [K]
  zH_i,zH_im1            , & ! ice thickness                          [m]
  zT_f                   , & ! freezing temperature for sea ice       [K]
  zT_f0                  , & ! melting temperature at top of thick ice layer
  zT_ftop                , & ! variable melting temperature at top
  zCshape_i              , & ! temperature profile shape factor
  zPhi_i_pr0             , & ! temperature shape function derivative at the ice-water interface
  zPhi_i_pr1             , & ! temperature shape function derivative at the air-ice interface
  zCshape_i_lin          , & ! shape factor for linear temperature profile
  zPhi_i_pr0_lin         , & ! shape function derivative at the ice-water interface for linear profile
  zPhi_i_pr1_lin         , & ! shape function derivative at the air-ice interface for linear profile
  zCshape_tunpar1        , & ! tuning parameter (used to compute the shape factor)
  zCshape_tunpar2        , & ! tuning parameter (used to compute the shape factor)
  zH_i_max                   ! maximum ice thickness [m] (tuning parameter)

REAL    (KIND=wp   )     ::  &
  zpsat_i, zqvap         , & ! statement functions and
  ztemp, zpvap, zpres        ! their formal arguments

CHARACTER (LEN=255)      :: yzerrmsg
CHARACTER (LEN=25)       :: yzroutine

! Local arrays:
! -------------

REAL    (KIND=wp   )     ::  &
  zt_s(ie,je)                ! = t_s   on land and sea
                             ! = t_ice on sea ice (if present)

! Tracer pointers
!----------------
REAL (KIND=wp),     POINTER :: &
  qv  (:,:,:)  => NULL()

!- End of header
!------------------------------------------------------------------------------

! Definition of statement functions:

! Saturation vapour pressure over ice
  zpsat_i(ztemp) = b1 * EXP( b2i*(ztemp-b3)/(ztemp-b4i) )

! Specific humidity
  zqvap(zpvap,zpres) = rdv * zpvap / ( zpres - o_m_rdv*zpvap )

!------------------------------------------------------------------------------

  yzroutine = 'seaice'

! Select time level and time step for calculations.
! A 2-time level scheme is used for the sea ice model. Therefore, it is not
! necessary to distinguish between leapfrog and Runge-Kutta scheme. But the
! first time step must not be done with dt/2, as in the leapfrog scheme.

  nx = nnow
  IF ( (ntstep == 0) .AND. (.NOT. l2tls) ) THEN
    ! Use the original dt and not dt/2
    zdt = 2 * dt
  ELSE
    zdt = dt
  END IF

  ! retrieve the required microphysics tracers (at timelevel nx)
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nx, ptr = qv)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

! Construct mask of ocean points, i.e. water points which are not lake points.

  IF (ntstep == nstart) THEN
    DO j = jstartpar,jendpar
      DO i = istartpar,iendpar

        lseamask(i,j) = .NOT. llandmask(i,j)

        IF (llake) THEN
          IF (depth_lk(i,j) > 0.0_wp) lseamask(i,j) = .FALSE.
        END IF

      END DO
    END DO
  END IF

! Physical 'constants'

  zT_f         = t0_melt + zt_ice              ! 273.15 - 1.7 K
  zT_f0        = t0_melt                       ! 273.15 K

! The ice shape function parameters (constants)
!  The shape factor and the shape function derivatives at the ice-water and at the air-ice interface
!  depend on the ice thickness. A linear temperature profile is recovered in case of thin ice.

  zCshape_i_lin    = 0.5_wp  ! zCshape_i for linear profile
  zPhi_i_pr0_lin   = 1.0_wp  ! zPhi_i_pr0 for linear profile
  zPhi_i_pr1_lin   = 1.0_wp  ! zPhi_i_pr1 for linear profile
  zH_i_max         = 3.0_wp  ! Maximum ice thickness [m]
                                 ! (recommended value)
  zCshape_tunpar1  = 2.0_wp  ! Shape factor tuning parameter
                                 ! (recommended value, allowable values range from -1 to 5)
  zCshape_tunpar2  = 1.0_wp/12.0_wp ! Shape factor tuning parameter
                                            ! (fixed  value, should not be altered)

! Derived 'constants'

  zrhoc_i      = rho_ice*zc_i          ! 1.932E+06 J/(m**3 K)
  zcidl        = zc_i/lh_f             ! 5.988E-03 1/K
  zrhoiLf      = rho_ice*lh_f
  zkicidLf     = zk_i*zc_i/lh_f        ! 7.545E-09 m**2/(K s)

! Compute the ice shape-function parameters (functions of the ice thickness).
! At the moment, these are set to constant values
! corresponding to a linear temperature profile within the ice.

  zPhi_i_pr0 = zPhi_i_pr0_lin
  zPhi_i_pr1 = zPhi_i_pr1_lin
  zCshape_i  = zCshape_i_lin

  DO j = jstartpar,jendpar
    DO i = istartpar,iendpar

SEA:  IF (lseamask(i,j)) THEN
        zH_im1 = h_ice(i,j,nx) ! Previous time step ice thickness
        zT_im1 = t_ice(i,j,nx) ! Previous time step ice surface temperature

! Compute the ice shape function parameters (functions of the ice thickness).
!_  Not used at the moment. When the first version of the ice model
!_  which uses a linear temperature profile within the ice goes operational,
!_  these changes should be incorporated into the COSMO model and tested.
!
!_nu  zPhi_i_pr0 = zH_im1/zH_i_max
!_nu  zPhi_i_pr1 = zPhi_i_pr1_lin + zCshape_tunpar1*zPhi_i_pr0
!_nu  zCshape_i = zCshape_i_lin - zCshape_tunpar2*(1.0_wp+zCshape_tunpar1)*zPhi_i_pr0
!_nu  zPhi_i_pr0 = zPhi_i_pr0_lin - zPhi_i_pr0

ICE:    IF (zH_im1 > 0.0_wp) THEN    ! Ice exists

!   Compute ice thickness dependent ice top temperature (Not Used)
!   zT_ftop = zT_f + (zT_f0-zT_f)*(1.-exp(-zH_im1/0.5_wp))

          zT_ftop = zT_f0

!   Compute surface turbulent heat fluxes

          jm1   = MAX(1,j-1)
          im1   = MAX(1,i-1)
          zuv   = 0.5_wp * SQRT( (u(i,j,ke,nx) + u(im1,j,ke,nx))**2 &
                                +(v(i,j,ke,nx) + v(i,jm1,ke,nx))**2 )

          zpatm = p0(i,j,ke) + pp(i,j,ke,nx)
          ztatm = t(i,j,ke,nx) * ((ps(i,j,nx)/zpatm)**rdocp)
          zqsat = zqvap(zpsat_i(zT_im1),ps(i,j,nx))

          ztvs  = zT_im1*(1.0_wp + rvd_m_o*zqsat)
          zrho  = ps(i,j,nx)/(r_d*ztvs)

          zshfl = zrho * tch(i,j) * zuv            * cp_d * (ztatm         - zT_im1)
          zlhfl = zrho * tch(i,j) * zuv * tfv(i,j) * lh_s * (qv(i,j,ke)    - zqsat )

!   Total atmospheric forcing / heat flux  [W/m**2]

          zthfl_a = sobs(i,j) + thbs(i,j) + zshfl + zlhfl

!   Ice internal temperature and heat fluxes, using linear temperature profile

          zq_i       = - zk_i*(zT_im1 - zT_f)/zH_im1
          zhfl_i     = zq_i*zrhoc_i

          IF (zT_im1>=zT_ftop .AND. zthfl_a>=0.0_wp) THEN ! Melting from above and below
            zH_i     = zH_im1 - zdt*(zthfl_a+(zPhi_i_pr1-zPhi_i_pr0)*zhfl_i)/zrhoiLf
            zT_i     = zT_ftop
          ELSE                                           ! Freezing or melting from below

            zH_i     = zH_im1 - zdt*zkicidLf*zPhi_i_pr0*(zT_im1-zT_f)/zH_im1
            zT_i     = zT_im1 + (zdt/zCshape_i/zH_im1)*(zthfl_a/zrhoc_i+zq_i*zPhi_i_pr0)
            zH_i     = min(zH_i, zH_i_max)     ! Security
            zT_i     = min(zT_i, zT_ftop)      ! Security
          END IF

!   Surface temperature will only be updated if originally existing sea ice melts completely

          IF (zH_i < zH_imin) THEN   ! Ice thickness less than threshold
            zH_i = 0.0_wp        ! Remove left overs and set
            zT_i = zT_f              ! temperature to freezing value
            t_s(i,j,nnew) = zT_i     ! Set the surface temperature to freezing value
          END IF

        ELSE  ICE ! No ice exists and we are not going to create any
          t_s(i,j,nnew) = t_s(i,j,nx)
          zT_i          = t_s(i,j,nnew)
          zH_i          = 0.0_wp
        END IF  ICE

        t_ice(i,j,nnew) = zT_i
        h_ice(i,j,nnew) = zH_i

      END IF SEA

    END DO
  END DO

!   Computation of the temperature at the boundary sea/ice - atmosphere

  DO j = jstartpar,jendpar
    DO i = istartpar,iendpar
      zt_s(i,j) = t_s(i,j,nnew)
      IF (.NOT. llandmask(i,j) .AND. h_ice(i,j,nnew) > 0.0_wp) THEN
        zt_s(i,j) = t_ice(i,j,nnew)
      ENDIF
    ENDDO
  ENDDO

  IF (lmulti_snow) THEN
    CALL tgcom ( t_g(:,:,nnew), t_snow_mult(:,:,1,nnew), &
                 zt_s(:,:)    , w_snow(:,:,nnew), &
                 llandmask(:,:) , ie, je, cf_snow,&
                 istartpar, iendpar, jstartpar, jendpar )
  ELSE
    CALL tgcom ( t_g(:,:,nnew), t_snow(:,:,nnew), &
                 zt_s(:,:)    , w_snow(:,:,nnew), &
                 llandmask(:,:) , ie, je, cf_snow,&
                 istartpar, iendpar, jstartpar, jendpar )
  ENDIF

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE seaice

!------------------------------------------------------------------------------
! End of module src_seaice
!------------------------------------------------------------------------------

END MODULE src_seaice
