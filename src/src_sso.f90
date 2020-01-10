!+ Source module for the sub-grid scale orography
!------------------------------------------------------------------------------

MODULE src_sso

!------------------------------------------------------------------------------
!
! Description:
!   The module "src_sso" performs calculations related to the parameterization
!   of sub-grid scale orographic (SSO) effects. The present SSO scheme is based
!   on Lott and Miller (1997).
!
!   All global variables of the model that are used by the SSO routines are
!   imported by USE statements below. The interface of the SSO routines and the
!   model is provided by the organizational routine "organize_sso".
!
!   The parameterization package has been extracted from the DWD global model
!   GME (V2_15 of the DWD gmtri library) and is originally based on an earlier
!   version by M. Miller and F. Lott (ECMWF). Some modifications have been made
!   in the code: Internal communication by common-blocks is replaced by module
!   parameters, scalars and arrays defined in this module.
!
! Current Code Owner: DWD, Jan-Peter Schulz
!  phone:  +49  69  8062 2742
!  fax:    +49  69  8062 3721
!  email:  jan-peter.schulz@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V4_4         2008/07/16 Jan-Peter Schulz
!  Initial release
! V4_8         2009/02/16 Guenther Zaengl
!  Use p0hl (reference pressure at half levels) for full consistency with
!  new reference atmosphere implementation
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
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
    ke1        ,    & ! ke + 1

! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the 
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from 
!    the other ones because of the use of the staggered Arakawa-C-grid.
!    
!   zonal direction
    istart     ,    & ! start index for the forecast of w, t, qd, qw and pp
    iend       ,    & ! end index for the forecast of w, t, qd, qw and pp
    istartu    ,    & ! start index for the forecast of u
    iendu      ,    & ! end index for the forecast of u
    istartv    ,    & ! start index for the forecast of v
    iendv      ,    & ! end index for the forecast of v
    istartpar  ,    & ! start index for computations in the parallel program
    iendpar    ,    & ! end index for computations in the parallel program

!   meridional direction
    jstart     ,    & ! start index for the forecast of w, t, qd, qw and pp
    jend       ,    & ! end index for the forecast of w, t, qd, qw and pp
    jstartu    ,    & ! start index for the forecast of u
    jendu      ,    & ! end index for the forecast of u
    jstartv    ,    & ! start index for the forecast of v
    jendv      ,    & ! end index for the forecast of v
    jstartpar  ,    & ! start index for computations in the parallel program
    jendpar    ,    & ! end index for computations in the parallel program

! 5. variables for the time discretization and related variables
! --------------------------------------------------------------
    dt         ,    & ! long time-step
    dt2               ! dt*2.

! end of data_modelconfig

!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 2. physical constants and related variables
! -------------------------------------------
    r_d        ,    & ! gas constant for dry air
    cp_d       ,    & ! specific heat of dry air at constant pressure
    g                 ! acceleration due to gravity

! end of data_constants

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------
    p0         ,    & ! reference pressure at full levels             ( Pa  )
    p0hl       ,    & ! reference pressure at half levels             ( Pa  )
    dp0        ,    & ! reference pressure thickness of layer         ( Pa  )
    hhl        ,    & ! geometrical height of half levels             ( m   )

! 2. external parameter fields                                        (unit)
! ----------------------------
    hsurf      ,    & ! height of surface topography                  ( m   )
    sso_stdh   ,    & ! standard deviation of sub-grid scale orography( m   )
    sso_gamma  ,    & ! anisotropy of sub-grid scale orography          --
    sso_theta  ,    & ! angle betw. principal axis of orography and E ( rad )
    sso_sigma  ,    & ! mean slope of sub-grid scale orography          --

! 3. prognostic variables                                             (unit)
! -----------------------
    u          ,    & ! zonal wind speed                              ( m/s )
    v          ,    & ! meridional wind speed                         ( m/s )
    t          ,    & ! temperature                                   (  k  )
    pp         ,    & ! deviation from the reference pressure         ( pa  )

! 5. fields for surface values and soil model variables               (unit )
! -----------------------------------------------------
    ps         ,    & ! surface pressure                              ( pa  )

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------
!   fields from the sub-grid scale orography scheme
    ut_sso     ,    & ! u-tendency due to SSO                         ( m/s2)
    vt_sso     ,    & ! v-tendency due to SSO                         ( m/s2)
    tt_sso     ,    & ! temperature tendency due to SSO               ( K/s )
    ustr_sso   ,    & ! u-stress (surface momentum flux) due to SSO   ( N/m2)
    vstr_sso   ,    & ! v-stress (surface momentum flux) due to SSO   ( N/m2)
    vdis_sso          ! vert. int. dissipation of kin. en. due to SSO ( W/m2)

! end of data_fields

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    ntstep     ,    & ! actual time step
                      ! indices for permutation of three time levels
    nold       ,    & ! corresponds to ntstep - 1
    nnow       ,    & ! corresponds to ntstep

! 5. additional control variables
! -------------------------------
    l2tls             ! forecast with 2-TL integration scheme

! end of data_runcontrol 

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Declarations

! The following parameters are tunable constants for the sub-grid scale
! orography scheme.

REAL (KIND = wp)     ::     &

! Tunable parameters
! ------------------
! Gsigcr  = 0.80_wp   , &   ! top layer for low level drag
! Gkdrag  = 0.30          , &   ! gw drag constant (Original ECMWF value)
  Gkdrag  = 0.075_wp  , &   ! gw drag constant
  Gkwake  = 0.50_wp   , &   ! low level wake drag constant
! Gkdrag  = Gkdrag_read   , &   ! Gkdrag_read read in or set in gme_tuning_constants
! Gkwake  = Gkwake_read   , &   ! Gkwake_read read in or set in gme_tuning_constants
  Grcrit  = 0.25_wp   , &   ! critical Richardson number
  Gfrcrit = 0.50_wp   , &   ! critical Froude number

! Security constants
! ------------------
  Gvsec   = 0.10_wp   , &   ! to secure the projection calculation
  Gssec   = 1.E-12_wp , &   ! to secure stability
  Gtsec   = 1.E-07_wp       ! to secure the stress calculation

!==============================================================================
! Module procedures in "src_sso"
!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure in "src_sso"
!------------------------------------------------------------------------------

SUBROUTINE organize_sso

!------------------------------------------------------------------------------
!
! Description:
!
!   The module procedure organize_sso is the interface of the model to the
!   parameterization package for sub-grid scale orographic effects.
!
! Externals:
!
!   sso: Lott and Miller SSO scheme
!
!------------------------------------------------------------------------------

! Local scalars and automatic arrays (also to be used in lower level routines):
! -----------------------------------------------------------------------------

  ! Input for the SSO routine "sso"
  ! -------------------------------
  REAL    (KIND=wp   )      ::  &
    zt   (ie,je,ke),    & ! temperature at full levels
    zu   (ie,je,ke),    & ! zonal wind component
    zv   (ie,je,ke),    & ! meridional wind component
    zfif (ie,je,ke),    & ! geopotential at full levels
    zpf  (ie,je,ke),    & ! pressure at full levels
    zph  (ie,je,ke1),   & ! pressure at half levels
    zdt                   ! actual timestep for 2 or 3 TL-scheme     
  LOGICAL   ::  &
    ldebug                ! debug indicator

  INTEGER (KIND=iintegers) ::  &
    i, j, k, km1, nx

! For error handling
! ------------------
  INTEGER (KIND=iintegers) ::  &
    izerror

  CHARACTER (LEN=80)       ::  &
    yzerrmsg
!
! End of header
!==============================================================================

  izerror  = 0
  yzerrmsg = '   '
  ldebug   = .FALSE.

  ! Select timelevel and timestep of the computation
  IF ( l2tls ) THEN
    nx  = nnow
    zdt = dt
  ELSE
    nx  = nold
    zdt = dt2
  ENDIF

  ! In order to save some CPU-time, the SSO scheme is called at fixed time
  ! increments nincsso (e.g. every 10 timsteps) and stores the SSO tendencies
  ! and other fields on global arrays which are held fixed in the intermediate
  ! steps. The time increment nincsso can be set on NAMELIST input.

  ! Reset the SSO tendencies, the surface momentum fluxes due to SSO and the
  ! dissipation due to SSO each time when the scheme is called.
! ut_sso  (:,:,:) = 0.0_wp
! vt_sso  (:,:,:) = 0.0_wp
! tt_sso  (:,:,:) = 0.0_wp
! ustr_sso(:,:  ) = 0.0_wp
! vstr_sso(:,:  ) = 0.0_wp
! vdis_sso(:,:  ) = 0.0_wp

  ! Prepare the input arrays for the SSO scheme.
  DO k = 1, ke 
    km1 = MAX ( 1, k-1 )
    DO j = jstart, jend
      DO i = istart, iend
        zt   (i,j,k) = t(i,j,k,nx)
        zu   (i,j,k) = 0.5_wp*( u(i,j,k,nx) + u(i-1,j,k,nx) )
        zv   (i,j,k) = 0.5_wp*( v(i,j,k,nx) + v(i,j-1,k,nx) )
        zfif (i,j,k) = 0.5_wp*g*( hhl(i,j,k) + hhl(i,j,k+1) )
        zpf  (i,j,k) = p0(i,j,k) + pp(i,j,k,nx)
        zph  (i,j,k) = p0hl(i,j,k)  &
                     + 0.5_wp*(pp(i,j,k,nx) + pp(i,j,km1,nx))
      ENDDO
    ENDDO
  ENDDO
  DO j = jstart, jend
    DO i = istart, iend
      zph  (i,j,ke1) = ps(i,j,nx)
    ENDDO
  ENDDO

  ! Call to sso.

  CALL sso (                                                 &
       zpf   , zph   , zfif  , zt      , zu , zv , g*hsurf , &
       sso_stdh, sso_gamma, sso_theta, sso_sigma,            &
       zdt   , ntstep,                                       &
       ldebug,                                               &
       ut_sso, vt_sso, tt_sso, ustr_sso, vstr_sso, vdis_sso  )

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE organize_sso

!==============================================================================
!+ Module procedure in "src_sso"
!------------------------------------------------------------------------------

SUBROUTINE sso (                                                       &
           ppf    , pph    , pfif   , pt       , pu , pv  , pfis     , &
           psso_stdh, psso_gamma, psso_theta, psso_sigma,              &
           pdt    , knstep ,                                           &
           ldebug ,                                                    &
           pdu_sso, pdv_sso, pdt_sso, pustr_sso, pvstr_sso, pvdis_sso  )

!------------------------------------------------------------------------------
!
! Purpose:   The module procedure sso performs the parameterisation of
!            sub-grid scale orographic effects according to Lott and
!            Miller (1997). It computes the physical tendencies of the
!            prognostic variables u, v and T due to vertical transports
!            by sub-grid scale orographically excited gravity waves.
!
! Method:    The scheme consists of three parts, i.e.
!             - the calculation of lowlevel drag due to wake effects of
!               the sub-grid scale orography
!             - the gravity wave stress and
!             - the stress profile in the vertical
!            The stress is computed based on the low level wind, the 
!            static stability and subgrid orographic parameters (i.e.
!            standard deviation of height, angle of the main orographic
!            axis relative to east, anisotropy and slope of the orography).
!            At each level a wave Richardson number is computed and by
!            requiring that its value is never less than a critical one,
!            a value of stress can be determined for each model level.
!            The critical Froude number of the flow determines the depth
!            of the low level drag.
!
! Reference: Paper on SSO scheme by Lott and Miller (1997).
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------

!     Input variables
!     ---------------

!     Grid scale variables
!     --------------------
      REAL(KIND=wp)     :: pt   (ie,je,ke)
!                            temperature at full levels         (K)
      REAL(KIND=wp)     :: pu   (ie,je,ke)
!                            zonal wind component               (m/s)
      REAL(KIND=wp)     :: pv   (ie,je,ke)
!                            meridional wind component          (m/s)
      REAL(KIND=wp)     :: pfif (ie,je,ke)
!                            geopotential at full levels        (m**2/s**2)
      REAL(KIND=wp)     :: pfis (ie,je)
!                            geopotential at surface            (m**2/s**2)
      REAL(KIND=wp)     :: pph  (ie,je,ke1)
!                            pressure at half levels            (Pa)
      REAL(KIND=wp)     :: ppf  (ie,je,ke)
!                            pressure at full levels            (Pa)

!     GWD-Parameters for each model grid point
!     (standard deviation, anisotropy, angle and slope)
!     -------------------------------------------------
      REAL(KIND=wp)     :: psso_stdh  (ie,je)
      REAL(KIND=wp)     :: psso_gamma (ie,je)
      REAL(KIND=wp)     :: psso_theta (ie,je)
      REAL(KIND=wp)     :: psso_sigma (ie,je)

      REAL(KIND=wp)     :: pdt           ! time step

      INTEGER knstep ! time step of integration

      LOGICAL ldebug ! debug control switch

!     Output variables
!     ----------------

!     Tendencies of T, u and v
!
      REAL(KIND=wp)     :: pdt_sso(ie,je,ke)
      REAL(KIND=wp)     :: pdv_sso(ie,je,ke)
      REAL(KIND=wp)     :: pdu_sso(ie,je,ke)

!     Surface (u,v) momentum fluxes and vertically integrated dissipation
!
      REAL(KIND=wp)     :: pustr_sso(ie,je)
      REAL(KIND=wp)     :: pvstr_sso(ie,je)
      REAL(KIND=wp)     :: pvdis_sso(ie,je)

!     Arrays and variables local to *sso* or used for communication
!     with higher level subroutines
!     -----------------------------

!     Indicators of significant levels for the operation of the scheme
      INTEGER  mcrit  (ie,je) 
      INTEGER  mkcrith(ie,je) 
      INTEGER  mkenvh (ie,je) 
      INTEGER  mknu   (ie,je) 
      INTEGER  mknu2  (ie,je) 

      LOGICAL  lo_sso (ie,je)

      REAL(KIND=wp)     :: zfi    (ie,je,ke) 
      ! geopotential minus surface geopotential (m**2/s**2)
      REAL(KIND=wp)     :: ztau   (ie,je,ke1) 
      ! gravity wave stress               (Pa)
      REAL(KIND=wp)     :: zstrdu (ie,je,ke1)
      ! flux of u-momentum (GWD and Blocking)
      REAL(KIND=wp)     :: zstrdv (ie,je,ke1)
      ! flux of v-momentum (GWD and Blocking) 
      REAL(KIND=wp)     :: zstab  (ie,je,ke1)
      ! squared Brunt-vaisala frequency   (1/s**2)
      REAL(KIND=wp)     :: zvph   (ie,je,ke1)
      ! wind profile projected onto plane of low level wind
      REAL(KIND=wp)     :: zrho   (ie,je,ke1)
      ! density at half levels            (kg/m**3)
      REAL(KIND=wp)     :: zri    (ie,je,ke1)
      ! mean flow Richardson number       (-) 
      REAL(KIND=wp)     :: zpsi   (ie,je,ke1)
      ! angle between flow and main axis of orography
      !
      REAL(KIND=wp)     :: zzdep  (ie,je,ke)
      !
      REAL(KIND=wp)     :: zdudt  (ie,je)
      ! sso-tendency for zonal wind    (m/s**2)
      REAL(KIND=wp)     :: zdvdt  (ie,je)
      ! sso-tendency for merid.wind    (m/s**2)
      REAL(KIND=wp)     :: zdtdt  (ie,je)
      ! sso-tendency for temperature   (K/s)
      REAL(KIND=wp)     :: zulow  (ie,je)
      ! u-component of low level wind  (m/s)
      REAL(KIND=wp)     :: zvlow  (ie,je)
      ! v-component of low level wind  (m/s)
      ! directional parameters (see *sso_setup*)
      REAL(KIND=wp)     :: zvidis (ie,je)
      REAL(KIND=wp)     :: zd1    (ie,je) 
      REAL(KIND=wp)     :: zd2    (ie,je) 
      REAL(KIND=wp)     :: zdmod  (ie,je)
!
!     Utility variables
!     -----------------
      REAL(KIND=wp)     :: zgdph,zcons1
      REAL(KIND=wp)     :: zdedt,zdis,zdelp,ztemp,zb,zc,zcs,zss,zconb,zabsv,zzd1
      REAL(KIND=wp)     :: zratio,zbet,zust,zvst,zdt2

      INTEGER j1,j2,j3      ! loop indices

!     Timestep is already set for 2TL or 3TL scheme, respectively,
!     in calling routine organize_sso.
!     zdt2  = 2._wp* pdt
      zdt2  = pdt

!     zcons1=1._wp/(G*pdt*2._wp)
      zcons1=1._wp/(G*zdt2)

!     Initialize tendencies and compute geopotential above ground
!     ===========================================================
      DO j3=1,ke
        DO j2=jstart,jend
          DO j1=istart,iend
            pdu_sso(j1,j2,j3) = 0.0_wp
            pdv_sso(j1,j2,j3) = 0.0_wp
            pdt_sso(j1,j2,j3) = 0.0_wp
            zfi    (j1,j2,j3) = pfif(j1,j2,j3)-pfis(j1,j2)
          END DO
        END DO
      END DO

!     Control operation of scheme by selection of points with standard
!     deviation of sub-grid scale orography > 10 m only
!     =================================================
      DO j2=jstart,jend
        DO j1=istart,iend
          IF (psso_stdh(j1,j2).GT.10._wp) THEN
            lo_sso(j1,j2)=.TRUE.
          ELSE
            lo_sso(j1,j2)=.FALSE.
          ENDIF
        END DO
      END DO

! ========================================================
!     Computation of basic state variables in *sso_setup*
! ========================================================

      CALL sso_setup (                                    &
         pph ,ppf ,pu ,pv ,pt ,zfi,                       &
         psso_stdh,psso_theta,psso_gamma,                 &
         knstep,                                          &
         lo_sso, ldebug,                                  &
         zrho  , zri   , zstab,ztau ,zvph ,zpsi , zzdep,  &
         zulow , zvlow ,zd1  ,zd2  ,zdmod,                &
                 mkcrith, mcrit, mkenvh,mknu,mknu2)

! ========================================================
!     Surface gravity wave stress amplitude
! ========================================================
 
      CALL gw_stress (                                &
         zrho,zstab,zvph,psso_stdh,psso_sigma,zdmod,  &
         knstep,                                      &
         lo_sso,ldebug,                               &
         ztau) 
 
! ========================================================
!     Gravity wave stress profile
! ========================================================

      CALL gw_profil(                                  &
         pph, zrho,zstab,zvph,zri,ztau,                &
         zdmod,psso_sigma,psso_gamma,psso_stdh,        &
                 mkcrith, mcrit , mkenvh, mknu,mknu2,  &
         lo_sso,ldebug )
 
! ========================================================
!     Computation of SSO effects' tendencies
! ========================================================

!     Initialisation of tendencies for ALL grid points
!     ------------------------------------------------
      DO j2=jstart,jend
        DO j1=istart,iend
        zvidis (j1,j2)=0.0_wp
        zdudt  (j1,j2)=0.0_wp
        zdvdt  (j1,j2)=0.0_wp
        zdtdt  (j1,j2)=0.0_wp
        END DO
      END DO
 
!     Compute and add low level drag tendencies to the GWD ones
!     ---------------------------------------------------------
      DO j3=1,ke
       DO j2=jstart,jend
        DO j1=istart,iend

        IF (lo_sso(j1,j2)) THEN

!       Gravity wave drag (cf. documentation EQ.4.13)
!       ---------------------------------------------
        zdelp = pph(j1,j2,j3+1)-pph(j1,j2,j3)
        ztemp = -G*(ztau(j1,j2,j3+1)-ztau(j1,j2,j3))                    &
                  /(zvph(j1,j2,ke1)*zdelp)
        zdudt(j1,j2)=(zulow(j1,j2)*zd1(j1,j2)-zvlow(j1,j2)*zd2(j1,j2))  &
                                    *ztemp/zdmod(j1,j2)
        zdvdt(j1,j2)=(zvlow(j1,j2)*zd1(j1,j2)+zulow(j1,j2)*zd2(j1,j2))  &
                                    *ztemp/zdmod(j1,j2)
        IF (j3 < 4) THEN
         zdudt(j1,j2)= SIGN(MIN(ABS(zdudt(j1,j2)),20._wp/3600._wp),zdudt(j1,j2))
         zdvdt(j1,j2)= SIGN(MIN(ABS(zdvdt(j1,j2)),20._wp/3600._wp),zdvdt(j1,j2))
        ENDIF

!       Low level drag ('blocking') (cf. documentation EQ.4.14 ff.)
!       -----------------------------------------------------------
        IF (j3.GE.mkenvh(j1,j2)) THEN
         zb  = 1.0_wp-0.18_wp*psso_gamma(j1,j2)-0.04_wp*psso_gamma(j1,j2)**2
         zc  = 0.48_wp*psso_gamma(j1,j2)+0.3_wp*psso_gamma(j1,j2)**2
         zcs = COS(zpsi(j1,j2,j3))**2
         zss = 1.0_wp-zcs
         zzd1  =zb*zcs+zc*zss
         zconb =zdt2*Gkwake*psso_sigma(j1,j2)/(2._wp*psso_stdh(j1,j2))
         zabsv =0.5_wp*SQRT(pu(j1,j2,j3)**2+pv(j1,j2,j3)**2)
         zratio=(zcs+psso_gamma(j1,j2)*zss)/(psso_gamma(j1,j2)*zcs+zss)
         zbet  =MAX(0._wp,2._wp-1._wp/zratio)*zconb*zzdep(j1,j2,j3)*zzd1*zabsv
!        Partially implicit tendency calculation
!        ---------------------------------------
         zdudt(j1,j2)=-pu(j1,j2,j3)/zdt2*(zbet/(1._wp+zbet)) 
         zdvdt(j1,j2)=-pv(j1,j2,j3)/zdt2*(zbet/(1._wp+zbet)) 
        END IF

        pdu_sso(j1,j2,j3)=zdudt(j1,j2)
        pdv_sso(j1,j2,j3)=zdvdt(j1,j2)
        zust             = pu(j1,j2,j3)+zdt2*zdudt(j1,j2)
        zvst             = pv(j1,j2,j3)+zdt2*zdvdt(j1,j2)
        zdis=0.5_wp*(pu(j1,j2,j3)**2+pv(j1,j2,j3)**2-zust**2-zvst**2)
        zdedt            = zdis/zdt2
        zdtdt(j1,j2)     = zdedt       /Cp_d
        pdt_sso(j1,j2,j3)= zdtdt(j1,j2)
        zvidis(j1,j2)=zvidis(j1,j2)+zdis*zdelp    ! de-activated

        ENDIF

        END DO
       END DO
      END DO     ! loop over vertical layers

! ======================================================================
!     Flux computations of original code of *GWDRAG* are not required in
!     COSMO model, but may be reactivated by uncommenting the following
!     statements and declaration of the appropriate arrays
! ======================================================================

!     Stress components and dissipation
!     ---------------------------------

      DO j2=jstart,jend
        DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
            pvdis_sso(j1,j2)=zcons1*zvidis(j1,j2)
          ELSE
            pvdis_sso(j1,j2)=0._wp
          ENDIF   
        END DO
      END DO

!     Initialize flux at top
!     ----------------------
      DO j2=jstart,jend
        DO j1=istart,iend
          zstrdu(j1,j2,1)=0._wp
          zstrdv(j1,j2,1)=0._wp
        END DO
      END DO

!     Increment flux based on tendency in each layer
!     ----------------------------------------------
      DO j3=1,ke
        DO j2=jstart,jend
          DO j1=istart,iend
            zgdph=-G  /(pph(j1,j2,j3+1)-pph(j1,j2,j3))
            zstrdu(j1,j2,j3+1)=pdu_sso(j1,j2,j3)/zgdph + zstrdu(j1,j2,j3)
            zstrdv(j1,j2,j3+1)=pdv_sso(j1,j2,j3)/zgdph + zstrdv(j1,j2,j3)
          END DO
        END DO
      END DO

!     Store flux at surface
!     ---------------------
      DO j2=jstart,jend
        DO j1=istart,iend
        pustr_sso(j1,j2)=zstrdu(j1,j2,ke1)
        pvstr_sso(j1,j2)=zstrdv(j1,j2,ke1)
        END DO
      END DO

!     Control printout
!     ----------------
      IF (ldebug) THEN
        DO j2=jstart,jend
          DO j1=istart,iend
          IF (j1.EQ.55 .and. j2.EQ.64) THEN
            PRINT *, ' '
            PRINT *, ' '
            PRINT *, ' Diagnosis SSO scheme j1=55 j2=64'
            PRINT *, ' '
            PRINT *, ' fis      : ', pfis       (j1,j2)
            PRINT *, ' sso_stdh : ', psso_stdh  (j1,j2)
            PRINT *, ' sso_gamma: ', psso_gamma (j1,j2)
            PRINT *, ' sso_theta: ', psso_theta (j1,j2)
            PRINT *, ' sso_sigma: ', psso_sigma (j1,j2)
            PRINT *, ' '
            PRINT *, '  j3    ph (Pa)     pf (Pa)     u(m/s)    ',  &
           'v (m/s)     T (K)      fif (m^2/s^2)     du_sso     ',  &
           'dv_sss     dt_sso'
          ENDIF
          ENDDO
        ENDDO
  
        DO j2=jstart,jend
          DO j1=istart,iend
          IF (j1.EQ.55 .and. j2.EQ.64) THEN
            DO j3=1,ke
            WRITE (*,'(i3, 9E13.6)') j3, pph(j1,j2,j3+1), ppf(j1,j2,j3),  &
            pu(j1,j2,j3), pv(j1,j2,j3), pt(j1,j2,j3), pfif(j1,j2,j3),     &
            pdu_sso(j1,j2,j3), pdv_sso(j1,j2,j3), pdt_sso(j1,j2,j3)
            ENDDO
          ENDIF
          ENDDO
        ENDDO
      ENDIF

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE sso

!==============================================================================
!+ Module procedure in "src_sso"
!------------------------------------------------------------------------------

SUBROUTINE sso_setup (                                      &
           pph ,ppf ,pu ,pv ,pt ,pfi,                       &
           psso_stdh,psso_theta,psso_gamma,                 &
           knstep,                                          &
           lo_sso, ldebug,                                  &
           prho  , pri   , pstab,ptau ,pvph ,ppsi , pzdep,  &
           pulow , pvlow ,pd1  ,pd2  ,pdmod,                &
                   kkcrith, kcrit, kkenvh,kknu,kknu2)

!------------------------------------------------------------------------------
!
! Purpose: The module procedure sso_setup sets up parameters for the
!          parameterization of sub-grid scale orographic effects:
!
!          - definition of various reference model levels
!          - Brunt-Vaisala frequency on half levels
!          - mean wind components in the layer between one and two
!            standard deviations of sso above ground
!          - geometry factors, relating the orientation of sso and wind
!          - Phillips parameters
!          - vertical wind profile in plane of gravity wave stress
!          - basic flow Richardson number
!          - calculation of depth of 'blocked' layer
!          - calculation of layer in which low level wave braking occurs
!          - calculation of assumed vertical profile of sso
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------

!     Input arrays and variables
!     ==========================
!
      REAL(KIND=wp)     :: pph (ie,je,ke1)
      REAL(KIND=wp)     :: ppf (ie,je,ke)
!
      REAL(KIND=wp)     :: pu  (ie,je,ke)
      REAL(KIND=wp)     :: pv  (ie,je,ke)
      REAL(KIND=wp)     :: pt  (ie,je,ke)
      REAL(KIND=wp)     :: pfi (ie,je,ke)

!     subgrid scale orography parameters
      REAL(KIND=wp)     :: psso_stdh (ie,je)
      REAL(KIND=wp)     :: psso_theta(ie,je) 
      REAL(KIND=wp)     :: psso_gamma(ie,je) 

      INTEGER knstep ! time step of integration

      LOGICAL lo_sso(ie,je)
      LOGICAL ldebug ! debug control switch

!     Output arrays
!     =============

      REAL(KIND=wp)     :: prho (ie,je,ke1)
!     density on half levels          (kg/m**3)
      REAL(KIND=wp)     :: pri  (ie,je,ke1)
!     mean flow Richardson number     (-)
      REAL(KIND=wp)     :: pstab(ie,je,ke1)
!     squared Brunt-Vaisala frequency (1/s**2)
      REAL(KIND=wp)     :: ptau (ie,je,ke1)
!     gravity wave stress profile     (Pa)
      REAL(KIND=wp)     :: pvph (ie,je,ke1) 
!     projected flow on half levels   (m/s)
      REAL(KIND=wp)     :: ppsi (ie,je,ke1)
!     angle between orography and blocked flow (1:ke)
!                            or low level flow (ke1)
      REAL(KIND=wp)     :: pzdep(ie,je,ke)
!     height dependency factor for 'blocking' tendency
      REAL(KIND=wp)     :: pulow(ie,je) 
!     low level zonal wind            (m/s)
      REAL(KIND=wp)     :: pvlow(ie,je)
!     low level meridional wind       (m/s)

!     directional parameters
      REAL(KIND=wp)     :: pd1  (ie,je)
      REAL(KIND=wp)     :: pd2  (ie,je)
      REAL(KIND=wp)     :: pdmod(ie,je)

      INTEGER kkcrith(ie,je)
!     maximum level for wave breaking

      INTEGER kcrit  (ie,je) 
!     critical level

      INTEGER kkenvh (ie,je) 
!     index of top of 'envelope' / blocking level

      INTEGER kknu   (ie,je)
!     model level at 4*stdh

      INTEGER kknu2  (ie,je)
!     model level at 3*stdh

!     local arrays and variables
!     ==========================

      REAL(KIND=wp)     :: zvpf   (ie,je,ke)
      ! projected flow on full levels (m/s)
      REAL(KIND=wp)     :: zdp    (ie,je,ke)
      ! pressure difference between layers
      REAL(KIND=wp)     :: zsqst  (ie,je,ke)
      REAL(KIND=wp)     :: znorm  (ie,je)
      REAL(KIND=wp)     :: znup   (ie,je)
      REAL(KIND=wp)     :: znum   (ie,je)
      REAL(KIND=wp)     :: znu    (ie,je) 

      REAL(KIND=wp)     :: zcons1, zcons2 ! utility constants
      REAL(KIND=wp)     :: zu             ! security for low level zonal wind
      REAL(KIND=wp)     :: zb             ! Phillips parameter B
      REAL(KIND=wp)     :: zc             ! Phillips parameter C
      REAL(KIND=wp)     :: zdelp          ! pressure thickness of layers
      REAL(KIND=wp)     :: zvt1,zvt2      ! utility variables for flow projection
      REAL(KIND=wp)     :: zst            ! utility variable for stability calculation
      REAL(KIND=wp)     :: zdwind         ! utility variable for wind shear calculation
      REAL(KIND=wp)     :: zwind          ! utility variable for proj. wind calculation
      REAL(KIND=wp)     :: zggeenv,zggeo,zgvar ! geopotential utility variables
      REAL(KIND=wp)     :: zhcrit 

      INTEGER mknub(ie,je)
      INTEGER mknul(ie,je)

      INTEGER mi3h              ! vertical loop limit
      INTEGER j1,j2,j3          ! loop variables   

      LOGICAL lo1  (ie,je,ke1)
      LOGICAL llo               ! utility switch

!     The following parameter is a tunable constant for the sub-grid scale
!     orography scheme
!     ================

      INTEGER (KIND=iintegers) :: Nktopg                    
                                ! number of topmost layer used to define low level
                                ! flow in case of high vertical model resolution

!-------------------------------------------------------------------------------------

! Begin subroutine

      Nktopg = ke               ! number of topmost layer used to defined low level

!     computational constants
!     =======================

!     mi3h =(ki3e-ki3s+1)/3
      mi3h =ke/3

      zcons1=1._wp/R_d
      zcons2=G**2/Cp_d
 
!C*         2.1     DEFINE LOW LEVEL WIND, PROJECT WINDS IN PLANE OF
!C*                 LOW LEVEL WIND, DETERMINE SECTOR IN WHICH TO TAKE
!C*                 THE VARIANCE AND SET INDICATOR FOR CRITICAL LEVELS.
!C
!     security on anisotropy factor and presetting of critical levels

      DO j2=jstart,jend
        DO j1=istart,iend
        psso_gamma(j1,j2) = MAX(psso_gamma(j1,j2),Gtsec)
        kknu      (j1,j2) = ke
        kknu2     (j1,j2) = ke
        mknub     (j1,j2) = ke
        mknul     (j1,j2) = ke
        lo1(j1,j2,ke1)    =.FALSE.    ! Initialize control variable
        END DO
      END DO

!
!!!!  define top of low level drag calculation layer (*kkcrit*)
!     and other critical levels
!     ============================================================
!
      DO j3=ke,mi3h,-1     ! vertical loop   
        DO j2=jstart,jend
          DO j1=istart,iend
          zhcrit          =4._wp*psso_stdh(j1,j2)
          lo1(j1,j2,j3)=((pfi(j1,j2,j3)/G).GT.zhcrit          )
            IF(lo1(j1,j2,j3).NEQV.lo1(j1,j2,j3+1)) THEN
            kknu(j1,j2)=j3  ! first layer with height > 4*stdh
            ENDIF
          END DO
        END DO
      END DO                ! end of vertical loop

      DO j3=ke,mi3h,-1    ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          zhcrit          =3._wp*psso_stdh(j1,j2)
          lo1(j1,j2,j3)=((pfi(j1,j2,j3)/G).GT.zhcrit          )
            IF(lo1(j1,j2,j3).NEQV.lo1(j1,j2,j3+1)) THEN
            kknu2(j1,j2)=j3 ! first layer with height > 3*stdh
            ENDIF
          END DO
        END DO
      END DO                ! end of vertical loop

      DO j3=ke,mi3h,-1    ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          zhcrit          =2._wp*psso_stdh(j1,j2)
          lo1(j1,j2,j3)=((pfi(j1,j2,j3)/G).GT.zhcrit          )
            IF(lo1(j1,j2,j3).NEQV.lo1(j1,j2,j3+1)) THEN
            mknub(j1,j2)=j3  ! first layer with height > 2*stdh
            ENDIF
          END DO
        END DO
      END DO                ! end of vertical loop

      DO j3=ke,mi3h,-1    ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          zhcrit          =psso_stdh(j1,j2)
          lo1(j1,j2,j3)=((pfi(j1,j2,j3)/G).GT.zhcrit          )
            IF(lo1(j1,j2,j3).NEQV.lo1(j1,j2,j3+1)) THEN
            mknul(j1,j2)=j3 ! first layer with height > 1*stdh
            ENDIF
          END DO
        END DO
      END DO                ! end of vertical loop
 
!     Confine critical level indices to be less or equal to Nktopg
!     ============================================================
      DO j2=jstart,jend
        DO j1=istart,iend
        kknu(j1,j2) =MIN(kknu(j1,j2),Nktopg)
        mknub(j1,j2)=MIN(mknub(j1,j2),Nktopg)
        IF(mknub(j1,j2).EQ.Nktopg) mknul(j1,j2)=ke
        IF(mknub(j1,j2).EQ.mknul(j1,j2)) mknub(j1,j2) = mknub(j1,j2) - 1
        END DO
      END DO

!     Initialize various arrays
!     =========================
      DO j2=jstart,jend
        DO j1=istart,iend
        prho (j1,j2,ke1) = 0.0_wp
        pstab(j1,j2,1  ) = 0.0_wp
        pstab(j1,j2,ke1) = 0.0_wp
        pri  (j1,j2,1  ) = 0.0_wp
        pri  (j1,j2,ke1) = 9999.0_wp
        ppsi (j1,j2,ke1) = 0.0_wp
        pvph (j1,j2,1)   = 0.0_wp
        pulow(j1,j2)     = 0.0_wp
        pvlow(j1,j2)     = 0.0_wp
        kkcrith(j1,j2)   = ke
        kkenvh(j1,j2)    = ke ! default for top of envelope layer
        kcrit(j1,j2)     = 1  ! default for critical level
        znu  (j1,j2)     = 0.0_wp
        znum (j1,j2)     = 0.0_wp
        lo1  (j1,j2,ke1) = .FALSE.
        END DO
      END DO
 
!     pressure thickness, density and Brunt-Vaisala frequency (squared)
!     =================================================================
      DO j3=ke,2,-1        ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
          zdp (j1,j2,j3) = ppf(j1,j2,j3)-ppf(j1,j2,j3-1)
!         density on half levels
          prho (j1,j2,j3) = 2._wp*pph(j1,j2,j3)*zcons1                 &
     &                     /(pt(j1,j2,j3)+pt(j1,j2,j3-1))
!         squared Brunt-Vaisala frequency on half levels
          pstab(j1,j2,j3)= 2._wp*zcons2/(pt(j1,j2,j3)+pt(j1,j2,j3-1))  &
     &                    *( 1._wp-Cp_d*prho(j1,j2,j3)                 &
     &                         *(pt(j1,j2,j3)-pt(j1,j2,j3-1))              &
     &                             / zdp(j1,j2,j3) )
!         security on Brunt-Vaisala frequency
          pstab(j1,j2,j3)=MAX(pstab(j1,j2,j3),gssec)
          zsqst(j1,j2,j3)=SQRT(pstab(j1,j2,j3))
          ENDIF
          END DO
        END DO
      END DO                ! end of vertical loop
 
!     Definition of blocked flow
!     ==========================
      DO j3=ke,mi3h,-1          ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
          IF(j3.GE.mknub(j1,j2).AND.j3.LE.mknul(j1,j2)) THEN
!         only layers with height between one and two stdh contribute 
          pulow(j1,j2) = pulow(j1,j2)                                     &
     &                  +pu   (j1,j2,j3)*(pph(j1,j2,j3+1)-pph(j1,j2,j3))
          pvlow(j1,j2) = pvlow(j1,j2)                                     &
     &                  +pv   (j1,j2,j3)*(pph(j1,j2,j3+1)-pph(j1,j2,j3))
          END IF
          END IF
          END DO
        END DO
      END DO                ! end of vertical loop
 
!     Division by pressure thickness of contributing layers and
!     determination of wind speed of blocked flow
      DO j2=jstart,jend
        DO j1=istart,iend
        IF(lo_sso(j1,j2)) THEN
        pulow(j1,j2) = pulow(j1,j2) /                                    &
     &              (pph(j1,j2,mknul(j1,j2)+1)-pph(j1,j2,mknub(j1,j2)))
        pvlow(j1,j2) = pvlow(j1,j2) /                                    &
     &              (pph(j1,j2,mknul(j1,j2)+1)-pph(j1,j2,mknub(j1,j2)))
        znorm(j1,j2)=MAX(SQRT(pulow(j1,j2)**2+pvlow(j1,j2)**2),Gvsec)
        pvph(j1,j2,ke1)=znorm(j1,j2)  ! projected flow at lower boundary
        END IF
        END DO
      END DO
 
!     Axes of subgrid scale orography and plane of profiles 
!     =====================================================
      DO j2=jstart,jend      
        DO j1=istart,iend   
        IF (lo_sso(j1,j2)) THEN
        llo=(pulow(j1,j2).LT.Gvsec).AND.(pulow(j1,j2).GE.-Gvsec)
          IF(llo) THEN
          ZU=pulow(j1,j2)+2._wp*Gvsec
          ELSE
          ZU=pulow(j1,j2)
          ENDIF
!       angle between principal axis of orog. and low-level wind direction
!       ------------------------------------------------------------------
        ppsi (j1,j2,ke1) = psso_theta(j1,j2)-ATAN(pvlow(j1,j2)/ZU)
!       Phillips parameters B and C
!       ---------------------------
        zb           = 1._wp-0.18_wp*psso_gamma(j1,j2)                      &
     &                   -0.04_wp*psso_gamma(j1,j2)**2
        zc           = 0.48_wp*psso_gamma(j1,j2)+0.3_wp*psso_gamma(j1,j2)**2
!       projection parameters D1 and D2 (see documentation)
        pd1  (j1,j2) = zb-(zb-zc)*(SIN(ppsi(j1,j2,ke1))**2)
        pd2  (j1,j2) = (zb-zc)*SIN(ppsi(j1,j2,ke1))                   &
     &                        *COS(ppsi(j1,j2,ke1))
        pdmod(j1,j2) = SQRT(pd1(j1,j2)**2+pd2(j1,j2)**2)
        END IF
        END DO
      END DO
 
!     projection of flow into plane of low level stress  (eq.4.7)
!     ===========================================================
      DO j3=1,ke                  ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
          zvt1 = pulow(j1,j2)*pu(j1,j2,j3)+pvlow(j1,j2)*pv(j1,j2,j3)
          zvt2 =-pvlow(j1,j2)*pu(j1,j2,j3)+pulow(j1,j2)*pv(j1,j2,j3)
          zvpf(j1,j2,j3)=(zvt1*pd1(j1,j2)+zvt2*pd2(j1,j2))             &
     &                      /(znorm(j1,j2)*pdmod(j1,j2))
          ENDIF

!      initialize stress array *ptau*, depth of blocked layer *pzdep*
!      and angle array *ppsi* (not for j3=ke1) and reset control
!      variable *llo1*
!      -----------------------------------------------------------------
          ptau(j1,j2,j3)  =0.0_wp
          pzdep(j1,j2,j3) =0.0_wp
          ppsi(j1,j2,j3)  =0.0_wp
          lo1(j1,j2,j3)   =.FALSE.
          END DO
        END DO
      END DO                ! end of vertical loop
!C
!     linear interpolation of projected flow to half levels
!     ========================================================
      DO j3=2,ke     ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
          pvph(j1,j2,j3)=                                        &
     &     ((pph(j1,j2,j3)-ppf(j1,j2,j3-1))*zvpf(j1,j2,j3) +     &
     &      (ppf(j1,j2,j3)-pph(j1,j2,j3  ))*zvpf(j1,j2,j3-1))    &
     &                         /zdp(j1,j2,j3)
            IF(pvph(j1,j2,j3).LT.Gvsec) THEN  ! critical layer
            pvph(j1,j2,j3)=Gvsec
            if (j3.lt.mknub(j1,j2)) kcrit(j1,j2  )=j3
            ENDIF
          ENDIF
          END DO
        END DO
      END DO                ! end of vertical loop
 
!     Brunt-Vaisala frequency and density for lowest level
!     ====================================================
      DO j3=mi3h,ke   ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
            IF(j3.GE.(mknub(j1,j2)+1).AND.j3.LE.mknul(j1,j2)) THEN
            zst=zcons2/pt(j1,j2,j3)*(1._wp-Cp_d*prho(j1,j2,j3)*     &
     &                    (pt(j1,j2,j3)-pt(j1,j2,j3-1))/zdp(j1,j2,j3))
            pstab(j1,j2,ke1)=pstab(j1,j2,ke1)+zst*zdp(j1,j2,j3)
            pstab(j1,j2,ke1)=MAX(pstab(j1,j2,ke1),Gssec)
            prho (j1,j2,ke1)= prho(j1,j2,ke1)                     &
     &                    +pph(j1,j2,j3)*2._wp*zdp(j1,j2,j3)  &
     &                    *zcons1/(pt(j1,j2,j3)+pt(j1,j2,j3-1))
            ENDIF
          ENDIF
          END DO
        END DO
      END DO                ! end of vertical loop
 
!     normalization
!     -------------
      DO j2=jstart,jend
        DO j1=istart,iend
        IF(lo_sso(j1,j2)) THEN
        pstab(j1,j2,ke1)=pstab(j1,j2,ke1)                    &
     &      /(ppf(j1,j2,mknul(j1,j2))-ppf(j1,j2,mknub(j1,j2)))
        pstab(j1,j2,ke1)=MAX(pstab(j1,j2,ke1),gssec)
        prho (j1,j2,ke1)=prho(j1,j2,ke1)                     &
     &      /(ppf(j1,j2,mknul(j1,j2))-ppf(j1,j2,mknub(j1,j2)))
        END IF
        END DO
      END DO
 
!     mean flow Richardson number on half levels
!     ==========================================
      DO j3=2,ke     ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
          zdwind=MAX(ABS(zvpf(j1,j2,j3)-zvpf(j1,j2,j3-1)),Gvsec)
          pri(j1,j2,j3)=pstab(j1,j2,j3)*(zdp(j1,j2,j3)                  &
     &            /(G*prho(j1,j2,j3)*zdwind))**2
          pri(j1,j2,j3)=MAX(pri(j1,j2,j3),Grcrit)
          ENDIF
          END DO
        END DO
      END DO                ! end of vertical loop

!     define top of 'envelope' layer (cf. eq.4.8)
!     ===========================================
      DO j3=2,ke-1     ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
            IF (j3.GE.kknu2(j1,j2)) THEN
            znum (j1,j2)=znu(j1,j2)
            zwind= (pulow(j1,j2)*pu(j1,j2,j3)                           &
     &             +pvlow(j1,j2)*pv(j1,j2,j3))/                         &
     &         MAX(SQRT(pulow(j1,j2)**2+pvlow(j1,j2)**2),Gvsec)
            zwind=MAX(ABS(zwind),Gvsec)
            zdelp=pph(j1,j2,j3+1)-pph(j1,j2,j3)
!           vertically integrated left side of eq.4.8
            znu(j1,j2) = znu(j1,j2) + (zdelp/G)*                             &
     &              ( (zsqst(j1,j2,j3+1)/prho(j1,j2,j3+1)                    &
     &                +zsqst(j1,j2,j3  )/prho(j1,j2,j3  ) )/2._wp)/zwind
            IF((znum(j1,j2).LE.gfrcrit).AND.(znu(j1,j2).GT.gfrcrit)          &
     &                          .AND.(kkenvh(j1,j2).EQ.ke))                &
     &      kkenvh(j1,j2)=j3     
            ENDIF
          ENDIF
          END DO
        END DO
      END DO                ! end of vertical loop

!     dynamical mixing height for the breaking of gravity waves
!     =========================================================
      DO j2=jstart,jend
        DO j1=istart,iend
        znup(j1,j2)=0.0_wp
        znum(j1,j2)=0.0_wp
        END DO
      END DO

      DO j3=ke-1,2,-1  ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
            IF (j3.LT.kkenvh(j1,j2)) THEN    ! only above envelope height
            znum(j1,j2)=znup(j1,j2)
            zwind=(pulow(j1,j2)*pu(j1,j2,j3)+pvlow(j1,j2)*pv(j1,j2,j3))/    &
     &            MAX(SQRT(pulow(j1,j2)**2+pvlow(j1,j2)**2),Gvsec)
            zwind=MAX(ABS(zwind),Gvsec)
            zdelp=pph(j1,j2,j3+1)-pph(j1,j2,j3)
            znup(j1,j2) = znup(j1,j2) + (zdelp/G)*                          &
     &            ( (zsqst(j1,j2,j3+1)/prho(j1,j2,j3+1)                     &
     &              +zsqst(j1,j2,j3  )/prho(j1,j2,j3  ) )/2._wp)/zwind
            IF((znum(j1,j2).LE.1.5_wp).AND.(znup(j1,j2).GT.1.5_wp)                &
     &                          .AND.(kkcrith(j1,j2).EQ.ke))              &
     &      kkcrith(j1,j2)=j3
            ENDIF
          ENDIF
          END DO
        END DO
      END DO                ! end of vertical loop

!     allow low level wave breaking only above height of 4*stdh 
      DO j2=jstart,jend
        DO j1=istart,iend
        kkcrith(j1,j2)=MIN(kkcrith(j1,j2),kknu(j1,j2))
        END DO
      END DO
 
!     directional information for flow blocking
!     =========================================
      DO j3=mi3h,ke     ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
          IF(j3.GE.kkenvh(j1,j2)) THEN  ! only within envelope layer
          llo=(pu(j1,j2,j3).LT.Gvsec).AND.(pu(j1,j2,j3).GE.-Gvsec)
            IF(llo) THEN
            ZU=pu(j1,j2,j3)+2._wp*Gvsec
            ELSE
            ZU=pu(j1,j2,j3)
            ENDIF
          ppsi(j1,j2,j3)=psso_theta(j1,j2)-ATAN(pv(j1,j2,j3)/ZU)
          ENDIF
          ENDIF
          END DO
        END DO
      END DO                ! end of vertical loop
 
!     assumed vertical profile of sso for blocking calculations
!     =========================================================
      DO j3=mi3h,ke     ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
          IF(j3.GE.kkenvh(j1,j2)) THEN   ! only within envelope layer
          zggeenv= MAX(1._wp,                                             &
     &     (pfi(j1,j2,kkenvh(j1,j2))+pfi(j1,j2,kkenvh(j1,j2)-1))/2._wp)
          zggeo  = MAX(pfi(j1,j2,j3),1._wp)
          zgvar  = MAX(psso_stdh(j1,j2)*G,1._wp)
          pzdep(j1,j2,j3)=SQRT((zggeenv-zggeo)/(zggeo+zgvar))
          END IF
          END IF
          END DO
        END DO
      END DO                ! end of vertical loop

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE sso_setup

!==============================================================================
!+ Module procedure in "src_sso"
!------------------------------------------------------------------------------

SUBROUTINE gw_stress (                                  &
           prho,pstab,pvph,psso_stdh,psso_sigma,pdmod,  &
           knstep,                                      &
           lo_sso, ldebug,                              &
           ptau )

!------------------------------------------------------------------------------
!
! Purpose: The module procedure gw_stress computes the gravity stress
!          amplitude following eq.4.11 of the documentation.
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------

!     Input
!     =====
!
      REAL(KIND=wp)     :: prho (ie,je,ke1)
      ! density on half levels    (kg/m**3)
      REAL(KIND=wp)     :: pstab(ie,je,ke1)
      ! squared Brunt-Vaisala frequency  (1/s**2)
      REAL(KIND=wp)     :: pvph (ie,je,ke1)
      ! wind on half levels projected on plane of low level wind (m/s)
      REAL(KIND=wp)     :: psso_stdh (ie,je)
      ! standard deviation of sso-height (m)
      REAL(KIND=wp)     :: psso_sigma(ie,je)
      ! mean sso-slope                   (-)
      REAL(KIND=wp)     :: pdmod(ie,je)
      ! projection parameter = SQRT(D1**2+D2**2)    cf. eq.4.7

      INTEGER knstep ! time step of integration

      LOGICAL ldebug ! debug control switch
      LOGICAL lo_sso(ie,je)
      ! 

!     Output
!     ======
      REAL(KIND=wp)     :: ptau(ie,je,ke1)
      ! gravity wave stress amplitude   (Pa)

!     local variables
!     ===============
                                  ! utility variables, which may be used to modify
      REAL(KIND=wp)     :: zblock ! the magnitude of the subgrid scale standard
      REAL(KIND=wp)     :: zeff   ! deviation which enters the stress amplitude
                                  ! calculation (zblock=0.0 in operational version)

      INTEGER j1,j2  ! loop variables

!     gravity wave stress amplitude (eq.4.11)
!     =======================================

      DO j2=jstart,jend
        DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
            zblock=0.0_wp 
            zeff=MAX(0._wp,2._wp*psso_stdh(j1,j2)-zblock)
            ptau(j1,j2,ke1)=Gkdrag*prho(j1,j2,ke1)*psso_sigma(j1,j2)  &
     &                     *zeff**2/4._wp                         &
     &                     /psso_stdh(j1,j2)*pvph(j1,j2,ke1)          &
     &                     *pdmod(j1,j2)*SQRT(pstab(j1,j2,ke1))
          ELSE
            ptau(j1,j2,ke1)=0.0_wp
          ENDIF
        END DO
      END DO

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE gw_stress

!==============================================================================
!+ Module procedure in "src_sso"
!------------------------------------------------------------------------------

SUBROUTINE gw_profil(                                    &
           pph, prho,pstab,pvph,pri,ptau,                &
           pdmod,psso_sigma,psso_gamma,psso_stdh,        &
                   kkcrith, kcrit , kkenvh, kknu,kknu2,  &
           lo_sso,ldebug )

!------------------------------------------------------------------------------
!
! Purpose: The module procedure gw_profil computes the vertical profile of
!          gravity wave stress.
!
! Method:  The stress profile for gravity waves is computed as follows:
!
!          - it is constant (no gravity wave drag) at all levels
!            between the ground and the top of the blocked layer kkenvh
!          - it decreases linearly with height from the top of the
!            blocked layer to 3 * stdh (kknu), to simulate lee waves or
!            non-linear gravity wave breaking
!          - at levels above (kknu) it is constant, except when the
!            wave encounters a critical level (kcrit) or when it breaks
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------

!     Input
!     =====
!
      REAL(KIND=wp)     :: pph  (ie,je,ke1)
      ! half level pressure           (Pa)
      REAL(KIND=wp)     :: prho (ie,je,ke1)
      ! density on half levels    (kg/m**3)
      REAL(KIND=wp)     :: pstab(ie,je,ke1)
      ! squared Brunt-Vaisala frequency  (1/s**2)
      REAL(KIND=wp)     :: pri  (ie,je,ke1)
      ! mean flow Richardson number      ( - )    
      REAL(KIND=wp)     :: pvph (ie,je,ke1)
      ! wind on half levels projected on plane of low level wind (m/s)
      REAL(KIND=wp)     :: ptau (ie,je,ke1)
      !gravity wave stress profile
      REAL(KIND=wp)     :: pdmod(ie,je)
      ! projection parameter = SQRT(D1**2+D2**2)    cf. eq.4.7
      REAL(KIND=wp)     :: psso_stdh (ie,je)
      ! standard deviation of sso-height (m)
      REAL(KIND=wp)     :: psso_sigma(ie,je)
      ! mean sso-slope                   (-)
      REAL(KIND=wp)     :: psso_gamma(ie,je)
      ! anisotropy factor of sso         (-)

!     various significant levels
      INTEGER kkcrith(ie,je)
      INTEGER kcrit  (ie,je)
      INTEGER kkenvh (ie,je)
      INTEGER kknu   (ie,je)
      INTEGER kknu2  (ie,je)

      LOGICAL lo_sso (ie,je)
      LOGICAL ldebug
     
!     local arrays and variables
!     ==========================

      REAL(KIND=wp)     :: zdz2   (ie,je,ke)
      REAL(KIND=wp)     :: ztau   (ie,je,ke1)
      REAL(KIND=wp)     :: znorm  (ie,je)
      REAL(KIND=wp)     :: zoro   (ie,je)

      REAL(KIND=wp)     :: zb,zdelp,zdelpt,zalpha,zalfa,zdel,zriw    ! utitility variables
      REAL(KIND=wp)     :: zsqri,zdz2n                               ! utitility variables

      INTEGER j1,j2,j3                 ! loop indices

      DO j2=jstart,jend
        DO j1=istart,iend
        IF(lo_sso(j1,j2)) THEN
        zoro(j1,j2) = psso_sigma(j1,j2)*pdmod(j1,j2)  &
     &               /(4._wp*MAX(psso_stdh(j1,j2),1.0_wp))
!HF  &               /4._wp/MAX(psso_stdh(j1,j2),1.0_wp)
        ztau(j1,j2,kknu(j1,j2)+1) = ptau(j1,j2,kknu(j1,j2)+1)
        ztau(j1,j2,ke1          ) = ptau(j1,j2,ke1          )
        ENDIF   
        END DO
      END DO

      DO j3=ke,2,-1     ! vertical loop

!     constant stress up to top of blocking layer
!     ===========================================
      DO j2=jstart,jend
        DO j1=istart,iend
        IF(lo_sso(j1,j2)) THEN
          IF(j3.GE.kknu2(j1,j2)) THEN
          ptau(j1,j2,j3)=ztau(j1,j2,ke1)
          ENDIF
        ENDIF   
        END DO
      END DO

!     wave displacement at next level
!     ===============================
      DO j2=jstart,jend
        DO j1=istart,iend
        IF(lo_sso(j1,j2)) THEN
          IF(j3.LT.kknu2(j1,j2)) THEN
          znorm(j1,j2)=Gkdrag*prho(j1,j2,j3)*SQRT(pstab(j1,j2,j3))  &
     &              *pvph(j1,j2,j3)*zoro(j1,j2)
          zdz2(j1,j2,j3)=ptau(j1,j2,j3+1)/MAX(znorm(j1,j2),Gssec)
          ENDIF
        ENDIF   
        END DO
      END DO

!     wave Richardson number, new wave displacement and stress
!     breaking evaluation and critical level
!     ========================================================
      DO j2=jstart,jend
        DO j1=istart,iend
        IF(lo_sso(j1,j2)) THEN
          IF(j3.LT.kknu2(j1,j2)) THEN    ! only above blocking layer
            IF((ptau(j1,j2,j3+1).LT.Gtsec).OR.(j3.LE.kcrit(j1,j2))) THEN
            ptau(j1,j2,j3)=0.0_wp    ! above critical levels
            ELSE
            zsqri=SQRT(pri(j1,j2,j3))
            zalfa=SQRT(pstab(j1,j2,j3)*zdz2(j1,j2,j3))/pvph(j1,j2,j3)
!HF         zriw=pri(j1,j2,j3)*(1._wp-zalfa)/(1+zalfa*zsqri)**2
            zriw=pri(j1,j2,j3)*(1._wp-zalfa)/(1._wp+zalfa*zsqri)**2
              IF(zriw.LT.Grcrit) THEN      ! breaking occurs
              zdel=4._wp/zsqri/Grcrit+1._wp/Grcrit**2+4._wp/Grcrit
              zb=1._wp/Grcrit+2._wp/zsqri
              zalpha=0.5_wp*(-zb+SQRT(zdel))
              zdz2n=(pvph(j1,j2,j3)*zalpha)**2/pstab(j1,j2,j3)
              ptau(j1,j2,j3)=znorm(j1,j2)*zdz2n
              ELSE
              ptau(j1,j2,j3)=znorm(j1,j2)*zdz2(j1,j2,j3)
              ENDIF
            ptau(j1,j2,j3)=MIN(ptau(j1,j2,j3),ptau(j1,j2,j3+1))
            ENDIF
          ENDIF
        ENDIF    
        END DO
      END DO

      END DO       ! end of vertical loop

!     reorganisation of stress profile, if breaking occurs at low levels
!     ==================================================================
      DO j2=jstart,jend
        DO j1=istart,iend
        IF(lo_sso(j1,j2)) THEN
        ztau(j1,j2,kkenvh(j1,j2)) =ptau(j1,j2,kkenvh(j1,j2))
        ztau(j1,j2,kkcrith(j1,j2))=ptau(j1,j2,kkcrith(j1,j2))
        ENDIF   
        END DO
      END DO

!     linear decrease between kkenvh and kkcrith
      DO j3=1,ke      ! vertical loop
        DO j2=jstart,jend
          DO j1=istart,iend
          IF(lo_sso(j1,j2)) THEN
          IF(j3.GT.kkcrith(j1,j2).AND.j3.LT.kkenvh(j1,j2))THEN
          zdelp=pph(j1,j2,j3)-pph(j1,j2,kkenvh(j1,j2))
          zdelpt=pph(j1,j2,kkcrith(j1,j2))-pph(j1,j2,kkenvh(j1,j2))
          ptau(j1,j2,j3)=ztau(j1,j2,kkenvh(j1,j2)) +                    &
     &        (ztau(j1,j2,kkcrith(j1,j2))-ztau(j1,j2,kkenvh(j1,j2)) )*  &
     &            zdelp/zdelpt
          ENDIF
          ENDIF            
          END DO
        END DO
      END DO       ! end of vertical loop

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE gw_profil

!------------------------------------------------------------------------------
! End of module src_sso
!------------------------------------------------------------------------------

END MODULE src_sso
