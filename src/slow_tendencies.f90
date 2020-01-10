!+ External procedure for computing the slow tendencies
!------------------------------------------------------------------------------

SUBROUTINE slow_tendencies ( wcon )

!------------------------------------------------------------------------------
!
! Description:
!   This external procedure is only used in src_leapfrog and calculates the
!   final slow tendencies for the prognostic variables u,v,w,pp and T and
!   updates the moisture variables at time level n+1 (nnew).
!
! Method:
!   Using the previously calculated tendencies from advection, physics and
!   adiabatic processes, which are stored on the tendency arrays, the
!   vertical advection and vertical diffusion is solved by a vertically
!   implicit scheme (modified Crank-Nicolson).
!   The resulting final tendencies are restored on the tendency arrays.
!   The water vapour, cloud water and cloud ice variables are updated directly.
!   Additionally, the moisture convertgence is stored on dqvdt.
!   Also, the surface fluxes of momentum, heat and moisture are stored.
!
! Current Code Owner: DWD, Michael Baldauf
!  phone:  +49  69  8062 2733
!  fax:    +49  69  8062 3721
!  email:  Michael.Baldauf@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 2.10       2001/07/24 Guenther Doms
!  Initial release
! 2.13       2002/01/18 Christoph Schraff
!  Bug correction: increase dimension of 'za1t','za2t' from ke to ke1.
! 2.14       2002/02/15 Guenther Doms
!  Bug correction in the tridiagonal matrix for w: loop now from
!  k = 3, ke-1 instead of accidentially from k = 2, ke-1.
! 2.18       2002/07/16 Ulrich Schaettler
!  adaptations to use slow_tendencies only in the 3 time level leapfrog scheme
!  the new 2 time level scheme does not need it (work by Almut Gassmann).
! 3.7        2004/02/18 Michael Baldauf
!  Selection of correct time level for cloud ice qi, depending on lprogprec
! 3.14       2005/01/25 Ulrich Schaettler
!  Correction in the computation of latent and sensible heat flux (to be 
!  consistent with lower boundary condition)
! 3.18       2006/03/03 Dmitrii Mironov
!  Add variables and fields used by the lake model FLake;
!  changed computation of the surface latent heat flux with regard to lake ice.
! 3.21       2006/12/04 MeteoSwiss / Ulrich Schaettler
!  Modifications to use Bechtold convection scheme
! V4_4         2008/07/16 Ulrich Schaettler
!  qitens is now also used by the Tiedtke scheme: therefore some comments 
!  have been removed.
! V4_9         2009/07/16 Ulrich Schaettler, Heike Vogel
!  Introduced treatment of variables for COSMO-ART
! V4_12        2010/05/11 Michael Baldauf
!  Corrected dimension of qvsflx
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_17        2011/02/24 Ulrich Blahak
!  Correction for ztch: introduced lower velocity limit vel_min (COMMENTED OUT FOR NOW)
! V4_18        2011/05/26 Ulrich Schaettler
!  Unified treatment of deposition, added surface concentration (Christoph Knote)
! V4_20        2011/08/31 Ulrich Schaettler
!  Activated lower velocity limit vel_min
! V4_23        2012/05/10 Oliver Fuhrer, Ulrich Schaettler
!  Splitted the computations for moisture divergence dqvdt into a different code block
!  Removed thereby (silly) computations " + a ... - a", which now gives numerically
!    different results
!  Removed switch lprogprec
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer
!  Replaced qx-variables by using them from the tracer module
! V4_26        2012/12/06 Anne Roches
!  Changes and technical adaptations to the tracer handling
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  MESSy interface introduced
! V4_28        2013/07/12 KIT, Ulrich Schaettler
!  Changes to adapt COSMO-ART to new tracer module: all dependencies to
!  COSMOART and POLLEN deleted, because this is now handled by the tracer module
!  Eliminated reference to kflat (not used) (US)
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
! V5_1         2014-11-28 Oliver Fuhrer, Anne Roches, Ulrich Blahak
!  Replaced ireals by wp (working precision) (OF)
!  Removed hacks for the tracer module  (AR)
!  Use t_g instead of t_s for flux calculations (UB)
! V5_4a        2016-05-10 Matthias Raschendorfer
!  Use new module turb_data instead of data_turbulence
!  Use imode_turb, lcpfluc now from turb_data
!  The lower boundary condition is now controlled by new switch lsflcnd
!      (and not imode_turb==0 any more: because meaning of imode_turb has changed)
!  The new switch (itype_vdif < 0) replaces (itype_turb/=3 .OR. imode_turb < 2)
!  Corrections for the treatment of the lower boundary condition and for the 
!      case that vertical turbulent diffusion is calculated by the new common 
!      routine in module 'turb_turbdiff'.
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


! 4. constants for the horizontal rotated grid and related variables
! ------------------------------------------------------------------

    degrad,       & ! factor for transforming degree to rad
    eddlon,       & ! 1 / dlon
    eddlat,       & ! 1 / dlat
    edadlat,      & ! 1 / (radius of the earth * dlat)

! 4. variables for the time discretization and related variables
! --------------------------------------------------------------

    dt,           & ! long time-step
    ed2dt,        & ! 1 / (2 * dt)
    dt2,          & ! 2 * dt
    epsass,       & ! eps for the Asselin-filter
    betasw,       & ! beta-variable for treatment of soundwaves
    vhmx_vol,     & ! maximum absolute horizontal wind in total model domain
    vhmx_cfl,     & ! maximum absolute horizontal wind velocity from CFL

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

    r_d,          & ! gas constant for dry air
    r_v,          & ! gas constant for water vapor
    rdv,          & ! r_d / r_v
    o_m_rdv,      & ! 1 - r_d/r_v
    rvd_m_o,      & ! r_v/r_d - 1
    cp_d,         & ! specific heat of dry air at constant pressure
    rdocp,        & ! r_d / cp_d
    rcpv,         & ! wcp_d/wcp_v - 1
    rcpl,         & ! wcp_d/wcp_l - 1
    lh_v,         & ! latent heat of vapourization
    lh_f,         & ! latent heat of fusion
    lh_s,         & ! latent heat of sublimation
    g,            & ! acceleration due to gravity
    gq,           & ! g*g
    gr              ! 1 / g

! end of data_constants

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------

    rho0       ,    & ! reference density at the full model levels    (kg/m3)
    dp0        ,    & ! reference pressure thickness of layers        ( Pa)
    p0         ,    & ! reference pressure at main levels             ( Pa)

! 3. prognostic variables                                             (unit)
! -----------------------

    u          ,    & ! zonal wind speed                              ( m/s )
    v          ,    & ! meridional wind speed                         ( m/s )
    w          ,    & ! vertical wind speed (defined on half levels)  ( m/s )
    t          ,    & ! temperature                                   (  k  )
    pp         ,    & ! deviation from the reference pressure         ( pa  )

! 4. tendency fields for the prognostic variables                     (unit )
! -----------------------------------------------
!    timely deviation  by diabatic and adiabatic processes
!    without sound-wave terms

    utens        ,  & ! u-tendency without sound-wave terms           ( m/s2)
    vtens        ,  & ! v-tendency without sound-wave terms           ( m/s2)
    wtens        ,  & ! w-tendency without sound-wave terms           ( m/s2

                      ! (defined on half levels )
    ttens        ,  & ! t-tendency without sound-wave terms           ( K/s )
    pptens       ,  & ! pp-tendency without sound-wave terms          (Pa/s )

! 5. fields for surface values and soil model variables               (unit )
! -----------------------------------------------------

    ps        ,     & ! surface pressure                              ( pa  )
!   t_s       ,     & ! temperature of the ground surface             (  k  )
    t_g       ,     & ! weighted surface temperature                  (  k  )
    qv_s      ,     & ! specific water vapor content on the surface   (kg/kg)
    fr_lake   ,     & ! lake fraction in a grid element [0,1]         (  -  )
    depth_lk  ,     & ! lake depth                                    (  m  )
    h_ice             ! ice thickness                                 (  m  )

USE data_fields     , ONLY :   &

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------

    rho       ,     & ! density of moist air

    a1t, a2t  ,     & ! implicit weight of vertical diffusion

!   turbulent coefficients in the atmosphere
!   (defined on half levels)
    tkvm       ,    & ! turbulent diffusion coefficients for momentum (m/s2 )
    tkvh       ,    & ! turbulent diffusion coefficients for heat     (m/s2 )
                      ! and moisture

!   turbulent coefficients at the surface
    tcm      ,      & ! turbulent diffusion coefficients for momentum   --
    tch      ,      & ! turbulent diffusion coefficients for heat       --
                      ! and moisture

!   fields that are computed in the dynamics
    dqvdt      ,    & ! threedimensional moisture convergence         (  1/s)
    qvsflx     ,    & ! surface flux of water vapour                  (kg/m2s)
    umfl_s     ,    & ! average u-momentum flux (surface)             ( N/m2)
    vmfl_s     ,    & ! average v-momentum flux (surface)             ( N/m2)
    shfl_s     ,    & ! average sensible heat flux (surface)          ( W/m2)
    lhfl_s            ! average latent heat flux (surface)            ( W/m2)

! end of data_fields

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    nstart,       & ! first time step of the forecast
    nstop,        & ! last time step of the forecast
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
    nold,         & ! corresponds to ntstep - 1
    nnow,         & ! corresponds to ntstep
    nnew,         & ! corresponds to ntstep + 1

! 3. controlling the physics
! --------------------------
    lphys,        & ! forecast with physical parametrizations
    llake,        & ! forecst with lake model FLake
    lgsp,         & ! forecast with grid-scale precipitation
    ltur,         & ! forecast with turbulent diffusion
    itype_turb,   & ! type of turbulent diffusion parametrization
    itype_vdif,   & ! type of vertical diffusion calculation
    l_cosmo_art,  & ! if .TRUE., run the COSMO_ART
    l_pollen        ! of pollen

! end of data_runcontrol

!------------------------------------------------------------------------------

USE data_flake,       ONLY:  &
    h_Ice_min_flk   ! Minimum ice thickness ( m )

!------------------------------------------------------------------------------

USE turb_data,        ONLY :   &
    lsflcnd,      & ! lower surface flux condition
    imode_turb,   & ! mode of turbulent diffusion parametrization
    lcpfluc,      & ! consideration of fluctuations of the heat capacity of air
    vel_min         ! minimal velocity scale [m/s]

!------------------------------------------------------------------------------

USE environment,      ONLY :   &
  model_abort

!------------------------------------------------------------------------------

USE data_parallel,    ONLY :   &
  my_cart_id

!------------------------------------------------------------------------------

USE numeric_utilities_rk,     ONLY :  &
  clipping

!------------------------------------------------------------------------------

USE data_tracer ,     ONLY :                                                  &
  T_ADV_ID      ,  T_ADV_ON       ,                                           &
  T_TURB_ID     ,  T_TURB_1D      ,                                           &
  T_CLP_ID      ,  T_CLP_ON       ,                                           &
  T_BBC_ID      ,  T_BBC_ZEROFLUX ,  T_BBC_ZEROVAL  , T_BBC_SURF_VAL

!------------------------------------------------------------------------------

USE src_tracer ,      ONLY :   &
  trcr_get_ntrcr, trcr_get, trcr_meta_get, trcr_errorstr

!------------------------------------------------------------------------------


!=======================================================================
 
IMPLICIT NONE
 
!=======================================================================
 
! Subroutine arguments:
! ---------------------
  REAL    (KIND=wp   ),     INTENT (IN)    ::  &
    wcon(ie,je,ke1)       ! contravariant vertical velocity

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    i,  j,  k, isp,      & !  Loop indices in lon., lat. and vert. direction
    iztrcr,              & !
    izt1, izt2,          & !
    izerror                !

  REAL    (KIND=wp   )     ::  &
    zadvfac, zturfac,    & !
    zbetav, zbetp, zbetm,& !
    ztkvz , ztkvl,       & !
    ztvw  , zpp  ,       & !
    ztvb  , zvbke,       & !
    zrho0s, zqx ,        & !
    zswb  , zswt ,       & !
    zgav  , zgat ,       & !
    zgcv  , zgct ,       & !
    zag   , zas  ,       & !
    zcg   , zcs  ,       & !
    zbg   , zdg  ,       & !
    zd1g  , zd2g ,       & !
    zdtrcrg,             & !
    zfac_surf,           & !
    za1t_surf, za2t_surf,& !
    znew  , zz   ,       & !
    ztmcmq, zdr  ,       & !
    zfpi  , zppa ,       & !
    zrwcv , zrwcl, zws,  & !
    zeddpb, zdp0r, zdtr    !

  LOGICAL                  ::  &
    lmassf,              & !
    lvertdiff

  CHARACTER (LEN=80)       :: yzerrmsg
  CHARACTER (LEN=25)       :: yzroutine

! Local (automatic) arrays:
! ------------------------
  REAL    (KIND=wp   )     ::  &
    za1t    (      ke1),   & !
    za2t    (      ke1),   & !
    zpia    (ie,je,ke),    & !
    zpianf  (ie,je,ke1),   & !
!   zcpa    (ie,je,ke),    & !
!   zcpanf  (ie,je,ke1),   & !
    ztmkvm  (ie,je,ke),    & !
    zeddpq  (ie,je,ke),    & !
    ztmkvh  (ie,je   ),    & !
    ztcm    (ie,je   ),    & !
    ztch    (ie,je   ),    & !
    zkh     (ie,je,ke),    & !
    zkm     (ie,je   ),    & !
    ztrcor  (ie,je   ),    &
    zgavx   (ie,je,ke),    & !
    zgcvx   (ie,je,ke),    & !
    zc      (ie,je,ke),    & ! Upper main diagonal of ...
    zd1     (ie,je,ke),    & ! Right hand side of ...
    zd2     (ie,je,ke),    & ! Right hand side of ...
    zdtrcr  (ie,je,ke),    & ! Right hand side of ...
    ze      (ie,je,ke),    & ! Soluton vector of ...
    zdphr   (ie,je)

  INTEGER(KIND=iintegers)  ::  &
    izadv   (trcr_get_ntrcr())    , & ! advection for all tracers
    izturb  (trcr_get_ntrcr())    , & ! turb mix for all tracers
    izbbc   (trcr_get_ntrcr())    , & ! bottom BC for all tracers
    izclip  (trcr_get_ntrcr())    , & ! clipping for all tracers
    izsp_adv_lf (trcr_get_ntrcr())    ! SP_ADV_LF for all tracers

  LOGICAL                  ::  &
    lzmss_flx(trcr_get_ntrcr())       ! MASSFLX_CLP for all tracers

! Tracer pointers:
!-----------------
  REAL (KIND=wp),     POINTER :: &
    ztrcr      (:,:,:) => NULL(),& ! tracer variable at tlev=nnew
    ztrcr_t1   (:,:,:) => NULL(),& ! tracer variable at tlev=izt1 (new/old)
    ztrcr_t2   (:,:,:) => NULL(),& ! tracer variable at tlev=izt2 (now/old)
    ztrcr_tens (:,:,:) => NULL(),& ! tracer tendency variable
    ztrcr_surf (:,:,:) => NULL(),& ! tracer surface variable
    qv_new     (:,:,:) => NULL(),& ! QV at nnew
    qv_old     (:,:,:) => NULL()   ! QV at nold

! Statement function zfpi for calculation of the exner function
! where the dummy argument zppa is pressure

  zfpi(zppa) = (1.E-5_wp*zppa)**rdocp
 
! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine slow_tendencies     
!------------------------------------------------------------------------------

  izerror  = 0_iintegers
  yzerrmsg = ''
  yzroutine= 'slow_tendencies'

!------------------------------------------------------------------------------
! Section 1a: Some preparations
!------------------------------------------------------------------------------

  ! Retrieve the required microphysics tracers
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnew, ptr = qv_new)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  CALL trcr_get(izerror, idt_qv, ptr_tlev = nold, ptr = qv_old)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  ! Retrieve the required metadata
  CALL trcr_meta_get(izerror, T_ADV_ID, izadv)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
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
  CALL trcr_meta_get(izerror, 'SP_ADV_LF', izsp_adv_lf)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_meta_get(izerror, 'MASSFLX_CLP', lzmss_flx)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  ! Setting of parameters for implicit calculation of vertical advection       
  ! (zbetav: weight for the n+1 time level           
  !          zbetav=0: centered, +1: full implicit,-1: explicit)
  zbetav = 0.0_wp
  zbetp  = 0.5_wp*( 1.0_wp + zbetav )
  zbetm  = 0.5_wp*( 1.0_wp - zbetav )

  IF (lcpfluc .EQV. .TRUE.) THEN
     zrwcv=rcpv
     zrwcl=rcpl
  ELSE
     zrwcv=0.0_wp
     zrwcl=0.0_wp
  END IF

  ! Setting of parameters for horizontal averaging of vertical diffusion
  ! coefficients: ztkvz, ztkvl
  ztkvz = 0.9_wp
  ztkvl = (1.0_wp-ztkvz)*0.25_wp

  ! Vertical mass redistribution scheme for cloud water on or off
  lmassf = .TRUE.

  ! In order to switch off the calculation of vertical diffusion,
  ! we use a local copy for the implicit weights a1t and a2t with values
  ! set to zero. This has to be done if the tendencies have alreay been
  ! calculated in turbdiff. An exception is the cloud ice, the tendencies
  ! have not beec calculated in turbdiff.
  ! This will be removed when a final form of time-integration has been found.

!MR  (itype_turb /= 3 .OR. imode_turb < 2) is now replaced by itype_vdif < 0
  IF (ltur .AND. (itype_vdif < 0)) THEN
    DO k = 1, ke1
      za1t(k) = a1t(k)
      za2t(k) = a2t(k)
    ENDDO
    lvertdiff =.TRUE.
    zfac_surf = 1.0_wp
  ELSE
    DO k = 1, ke1
      za1t(k) = 0.0_wp
      za2t(k) = 0.0_wp
    ENDDO
    lvertdiff =.FALSE.
    zfac_surf = 0.0_wp
  END IF

  ! lvertdiff = .FALSE. is not implemented for the tracers
!MR:  only necessary until the case "itype_vdif=-2" is removed
  IF ((itype_turb ==  3) .AND. (imode_turb >= 2)      .AND.            &
      (itype_vdif == -2)                              .AND.            &
      (trcr_get_ntrcr().GT.0_iintegers)                      ) THEN
    yzerrmsg =  'This vertical diffusion type is not yet ' // &
                'implemented for tracers'
    izerror  =  5000
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  ! Selection of lower boundary conditions

  IF (.NOT. lsflcnd) THEN
    ! condition of the lower boundary is the surface concentration
    za1t_surf = za1t(ke1)  ! implicit weight for ke-value
    za2t_surf = za2t(ke1)  ! explicit weight for ke-value
  ELSE
    ! condition of the lower boundary is the explicit surface mass flux
    za1t_surf = 0.0_wp    ! no implicit weight
    za2t_surf = zfac_surf ! full explicit weight
  ENDIF

  ! Calculation of Exner-pressure and cp-factor         
  ! -------------------------------------------
  DO  k = 1 , ke
    DO  j = jstart-1, jend+1
      DO  i = istart-1, iend+1
        zpp         = p0(i,j,k) + pp(i,j,k,nnow)
        zpia(i,j,k) = zfpi( zpp )
!       zcpa(i,j,k) = 1_wp-zrwcv*qv(i,j,k,nnow)+zrwcl*qc(i,j,k,nnow)
      ENDDO     
    ENDDO
  ENDDO 

  ! Calculation of modified transfer coefficients 
  ! ---------------------------------------------
  DO   j  =  jstart, jend+1
    DO i  =  istart, iend+1
     zvbke      = 0.5_wp*SQRT ( (u(i,j,ke,nold) + u(i-1,j,ke,nold))**2     &
                            +(v(i,j,ke,nold) + v(i,j-1,ke,nold))**2 )
!    ztvb       = t_s (i,j,nold)*(1.0_wp + rvd_m_o*qv_s(i,j,nold))
     ztvb       = t_g (i,j,nold)*(1.0_wp + rvd_m_o*qv_s(i,j,nold))
     ztcm(i,j)  = tcm(i,j)*zvbke*g*ps(i,j,nold)/(r_d*ztvb)
     ztch(i,j)  = tch(i,j)*MAX(zvbke,vel_min)*g*ps(i,j,nold)/(r_d*ztvb)
    ENDDO       
  ENDDO


  ! Calculation of the modified turbulent diffusion coefficients,
  ! of the Exner function on half levels (zpianf) and
  ! of the cp-factor on half levels (zpianf)
  ! ------------------------------------------------------------

  DO  k = 2 , ke
    DO    j = jstart-1, jend+1
      DO  i = istart-1, iend+1
        zdr    = 1.0_wp / ( dp0(i,j,k-1)+dp0(i,j,k) )
        zswb   = dp0(i,j,k-1)*zdr
        zswt   = dp0(i,j,k  )*zdr
        zrho0s =   zswb*rho0(i,j,k  )*rho(i,j,k  )          &
                 + zswt*rho0(i,j,k-1)*rho(i,j,k-1)
        zdr    = gq*zrho0s*2.0_wp * zdr
        ztmkvh(i,j  ) = tkvh(i,j,k)*zdr
        ztmkvm(i,j,k) = tkvm(i,j,k)*zdr
        zpianf(i,j,k) = zswb*zpia(i,j,k) + zswt*zpia(i,j,k-1)
!       zcpanf(i,j,k) = zswb*zcpa(i,j,k) + zswt*zcpa(i,j,k-1)
      ENDDO     
    ENDDO       
    DO    j = jstart, jend
      DO  i = istart, iend
          zkh(i,j,k) = ztkvz*ztmkvh(i,j)                     &
                     + ztkvl*( ztmkvh(i-1,j) + ztmkvh(i+1,j) &
                             + ztmkvh(i,j-1) + ztmkvh(i,j+1) ) 
      ENDDO     
    ENDDO       
  ENDDO         
  DO    j = jstart-1, jend+1
    DO  i  = istart-1, iend+1
      zpianf(i,j,ke1) = zfpi( ps(i,j,nnow) )
!     zcpanf(i,j,ke1) = 1_wp-zrwcv*qv_s(i,j,nold)+zrwcl*qc(i,j,ke,nold)
    ENDDO        
  ENDDO        

  ! Precalculate some help variables to avoid divisions in subsequent code
  !-----------------------------------------------------------------------
  DO  k = 1, ke
    DO j = jstart, jend
      DO i = istart, iend
        zdr  = 1.0_wp / dp0(i,j,k)
        zeddpq(i,j,k) = rho0(i,j,k)*zdr / rho(i,j,k)
        zgavx (i,j,k) = - 0.5_wp*wcon(i,j,k  )*zdr
        zgcvx (i,j,k) =   0.5_wp*wcon(i,j,k+1)*zdr
      ENDDO
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section 1b: Computation of the moisture divergence
!------------------------------------------------------------------------------

  ! First, the matrix elements a(k), b(k) and c(k) of the coefficient matrix
  ! are set (these are the same for qv and qc ).
  ! The right hand side is stored on d (k).

  ! Top layer
  DO j = jstart, jend
    DO i = istart, iend
      zgct       = - zkh(i,j,2)*zeddpq(i,j,1)
      zcg        = zgcvx(i,j,1)*zbetp + zgct*za1t(2)
      zcs        = zgcvx(i,j,1)*zbetm + zgct*za2t(2)
      zbg        = ed2dt - zcg
      zd1(i,j,1)  = ed2dt * qv_old(i,j,1)                                     &
                    - zcs *(qv_old(i,j,2) - qv_old(i,j,1))                    &
                    + dqvdt(i,j,1)
      zd1(i,j,1)  = zd1(i,j,1)/zbg
      zc (i,j,1)  = zcg/zbg
    ENDDO
  ENDDO

  ! The layers from k=2 to k=ke-1
  DO k = 2, ke-1
    DO j = jstart, jend
      DO i = istart, iend
        zgat       = - zkh(i,j,k  )*zeddpq(i,j,k)
        zgct       = - zkh(i,j,k+1)*zeddpq(i,j,k)
        zag        = zgavx(i,j,k)*zbetp + zgat*za1t(k)
        zas        = zgavx(i,j,k)*zbetm + zgat*za2t(k)
        zcg        = zgcvx(i,j,k)*zbetp + zgct*za1t(k+1)
        zcs        = zgcvx(i,j,k)*zbetm + zgct*za2t(k+1)
        zbg        = ed2dt - zag - zcg
        zd1g       = ed2dt *  qv_old(i,j,k  )                                 &
                     - zas * (qv_old(i,j,k-1) - qv_old(i,j,k))                &
                     - zcs * (qv_old(i,j,k+1) - qv_old(i,j,k))                &
                     + dqvdt(i,j,k)
        zz         = 1.0_wp/( zbg - zag*zc(i,j,k-1) )
        zc (i,j,k) = zcg*zz
        zd1(i,j,k) = ( zd1g -zag*zd1(i,j,k-1) )*zz
      ENDDO
    ENDDO
  ENDDO

  ! the bottom layer
  DO j = jstart, jend
    DO i = istart, iend
      zgat       = - zkh (i,j,ke)*zeddpq(i,j,ke)
      zgct       = - ztch(i,j   )*zeddpq(i,j,ke)
      zag        = zgavx(i,j,ke)*zbetp + zgat*za1t(ke)
      zas        = zgavx(i,j,ke)*zbetm + zgat*za2t(ke)
      zcg        = zgct*za1t_surf
      zcs        = zgct*za2t_surf
      zbg        = ed2dt - zag - zcg
      zd1g       = ed2dt *  qv_old(i,j,ke  )                                  &
                   - zas * (qv_old(i,j,ke-1) - qv_old(i,j,ke))                &
!MR: Bug???        + zcs *  qv_old(i,j,ke  ) - zgct * qv_s(i,j,nold)          &
                   + zcs * (qv_old(i,j,ke  ) - qv_s(i,j,nold)) - zcg*qv_s(i,j,nold) &
                   + dqvdt(i,j,ke)
      zz         = 1.0_wp / ( zbg  - zag*zc (i,j,ke-1) )
      znew            = ( zd1g - zag*zd1(i,j,ke-1) ) * zz
      dqvdt(i,j,ke)   = ( MAX(0.0_wp, znew) - qv_old(i,j,ke) )*ed2dt
      ze   (i,j,ke)   = znew
    ENDDO
  ENDDO

  ! Backsubstitution and storage of the complete slow tendencies

  DO k = ke-1, 1, -1
    DO j =jstart, jend
      DO i = istart, iend
        ze(i,j,k     ) = zd1(i,j,k) - zc(i,j,k)*ze(i,j,k+1)
      ENDDO
    ENDDO
  ENDDO

!MR: extra treatment for dqvdt
  IF (.NOT.ltur .OR. (itype_vdif <= 0)) THEN
    ! 'dqvdt' needs to be loeaded by qv-convergence due to vertical advection
    !    (if (.NOT.ltur .OR. itype_vdif == 0)) 
    ! or due to vertical advection and diffusion
    !    (if (ltur .AND. itype_vdif < 0)):

    DO k = ke, 1, -1
      DO j =jstart, jend
        DO i = istart, iend
          dqvdt(i,j,k  ) = (MAX(0.0_wp,ze(i,j,k)) - qv_old(i,j,k))*ed2dt
        ENDDO
      ENDDO
    ENDDO

  ELSE   !  if (ltur .AND. itype_vdif > 0)
    ! vertical qv-diffusion is already contained in 'dqvdt' and needs to be
    ! incremented by the part due to vertical advection:

    DO k = ke, 1, -1
      DO j =jstart, jend
        DO i = istart, iend
          dqvdt(i,j,k  ) = dqvdt(i,j,k) + &
                           (MAX(0.0_wp,ze(i,j,k)) - qv_old(i,j,k))*ed2dt
        ENDDO
      ENDDO
    ENDDO
  ENDIF

!------------------------------------------------------------------------------
! Section 2a: Setup of tridiagonal matrix systems resulting from the implicit
!             numerical formulation of advection and diffusion of the tracers.
!------------------------------------------------------------------------------

  ! loop over tracers
  DO iztrcr = 1, trcr_get_ntrcr()

    ! switch on/off vertical diffusion and vertical advection
    IF ( izturb(iztrcr) == T_TURB_1D .AND. lvertdiff ) THEN
      zturfac = 1.0_wp
    ELSE
      zturfac = 0.0_wp
    ENDIF
    IF ( izadv(iztrcr) == T_ADV_ON ) THEN
      zadvfac = 1.0_wp
    ELSE
      zadvfac = 0.0_wp
    ENDIF
    
    ! get pointer to tracer (at nnew) as well as tendency
    CALL trcr_get(izerror, iztrcr, ptr_tlev=nnew, ptr=ztrcr,            &
                  ptr_tens=ztrcr_tens)
    IF (izerror /= 0) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine) 
    ENDIF
    
    ! izt1 is timelevel used for time-discretization (i.e. nnew = izt1 * dt2 + RHS)
    ! izt2 is the timelevel used for the turb. diffusions term in the RHS
    ! NOTE: izsp_adv_lf is typically only =2 for QI
    IF ( izsp_adv_lf(iztrcr) == 2_iintegers ) THEN
      izt1 = nnew
      izt2 = nold
    ELSE
      izt1 = nold
      izt2 = nold
    ENDIF
      
    CALL trcr_get(izerror, iztrcr, ptr_tlev=izt1, ptr=ztrcr_t1)
    IF (izerror /= 0) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
    
    ! in LF core vertical turbulent diffusion and implicit vertical advection is not
    ! computed for qr, qs, qg species which are advected using the fully 3D semi-Lagrangian
    ! scheme (izsp_adv_lf = 3)
    IF ( izsp_adv_lf(iztrcr) /= 3_iintegers ) THEN

      CALL trcr_get(izerror, iztrcr, ptr_tlev=izt2, ptr=ztrcr_t2)
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF
      
      IF (izbbc(iztrcr) == T_BBC_SURF_VAL) THEN
        CALL trcr_meta_get(izerror, iztrcr, 'SURF_FIELD', ztrcr_surf)
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
      
      ! Top layer
      DO j = jstart, jend
        DO i = istart, iend
          zgct          = - zkh(i,j,2)*zeddpq(i,j,1)
          IF ( izsp_adv_lf(iztrcr) == 2_iintegers ) THEN
            zcg           = zturfac*zgct* a1t(2)
            zcs           = zturfac*zgct* a2t(2) 
          ELSE 
            zcg           = zadvfac*zgcvx(i,j,1)*zbetp + zturfac*zgct*a1t(2) 
            zcs           = zadvfac*zgcvx(i,j,1)*zbetm + zturfac*zgct*a2t(2)
          ENDIF
          zbg             = ed2dt - zcg
          zdtrcr(i,j,1)   = ed2dt * ztrcr_t1(i,j,1) + ztrcr_tens(i,j,1)        &
                          - zcs *(ztrcr_t2(i,j,2) - ztrcr_t2(i,j,1))
          zdtrcr(i,j,1)   = zdtrcr(i,j,1)/zbg
          zc (i,j,1)      = zcg/zbg
        ENDDO
      ENDDO
      
      ! The layers from k=2 to k=ke-1
      DO k = 2, ke-1
        DO j = jstart, jend
          DO i = istart, iend
            zgat       = - zkh(i,j,k  )*zeddpq(i,j,k)
            zgct       = - zkh(i,j,k+1)*zeddpq(i,j,k)
            IF ( izsp_adv_lf(iztrcr) == 2_iintegers ) THEN
              zag        = zturfac*zgat*a1t(k)
              zas        = zturfac*zgat*a2t(k)
              zcg        = zturfac*zgct*a1t(k+1)
              zcs        = zturfac*zgct*a2t(k+1)
            ELSE
              zag        = zadvfac*zgavx(i,j,k)*zbetp + zturfac*zgat*a1t(k)
              zas        = zadvfac*zgavx(i,j,k)*zbetm + zturfac*zgat*a2t(k)
              zcg        = zadvfac*zgcvx(i,j,k)*zbetp + zturfac*zgct*a1t(k+1)
              zcs        = zadvfac*zgcvx(i,j,k)*zbetm + zturfac*zgct*a2t(k+1)
            ENDIF
            zbg        = ed2dt - zag - zcg
            zdtrcrg    = ed2dt * ztrcr_t1(i,j,k) + ztrcr_tens(i,j,k)          &
                       - zas * (ztrcr_t2(i,j,k-1) - ztrcr_t2(i,j,k))          &
                       - zcs * (ztrcr_t2(i,j,k+1) - ztrcr_t2(i,j,k))
            zz         = 1.0_wp/( zbg - zag*zc(i,j,k-1) )
            zc (i,j,k) = zcg*zz
            zdtrcr(i,j,k) = ( zdtrcrg -zag*zdtrcr(i,j,k-1) )*zz
          ENDDO
        ENDDO
      ENDDO
      
      ! the bottom layer
      DO j = jstart, jend
        DO i = istart, iend
          zgat       = - zkh (i,j,ke)*zeddpq(i,j,ke)
          zgct       = - ztch(i,j   )*zeddpq(i,j,ke)
          IF ( izsp_adv_lf(iztrcr) == 2_iintegers ) THEN
            zag        = zturfac*zgat*a1t(ke)
            zas        = zturfac*zgat*a2t(ke)
          ELSE
            zag        = zadvfac*zgavx(i,j,ke)*zbetp + zturfac*zgat*a1t(ke)
            zas        = zadvfac*zgavx(i,j,ke)*zbetm + zturfac*zgat*a2t(ke)
          ENDIF
          zcg        = zturfac*zgct*za1t_surf
          zcs        = zturfac*zgct*za2t_surf
          zbg        = ed2dt - zag - zcg
          IF (izbbc(iztrcr) == T_BBC_ZEROFLUX) THEN
!MR: last line added
            zdtrcrg      = ed2dt * ztrcr_t1(i,j,ke) + ztrcr_tens(i,j,ke)      &
                         - zas * (ztrcr_t2(i,j,ke-1) - ztrcr_t2(i,j,ke))      &
                         - zcg *                       ztrcr_t2(i,j,ke)
          ELSEIF (izbbc(iztrcr) == T_BBC_ZEROVAL) THEN
            zdtrcrg      = ed2dt * ztrcr_t1(i,j,ke) + ztrcr_tens(i,j,ke)      &
                         - zas * (ztrcr_t2(i,j,ke-1) - ztrcr_t2(i,j,ke))      &
                         + zcs *  ztrcr_t2(i,j,ke) 
          ELSEIF (izbbc(iztrcr) == T_BBC_SURF_VAL) THEN
!MR: last 2 lines changed
            zdtrcrg      = ed2dt * ztrcr_t1(i,j,ke) + ztrcr_tens(i,j,ke)      &
                         - zas * (ztrcr_t2(i,j,ke-1) - ztrcr_t2(i,j,ke))      &
                         + zcs * (ztrcr_t2(i,j,ke) - ztrcr_surf(i,j,nold))    &
                         - zcg *                     ztrcr_surf(i,j,nold)
          ELSE
            WRITE(yzerrmsg,'(i3)') iztrcr
            yzerrmsg =  'This BBC option is not valid for the tracers'//      &
                        ' (TRCR='//TRIM(yzerrmsg)//')'      
            izerror  =  5000 
            CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
          ENDIF
          zz         = 1.0_wp / ( zbg  - zag*zc (i,j,ke-1) )
          ztrcr(i,j,ke) = ( zdtrcrg - zag*zdtrcr(i,j,ke-1) ) * zz
        ENDDO
      ENDDO
      
      ! Backsubstitution and storage of the complete slow tendencies
      
      DO k = ke-1, 1, -1
        DO j =jstart, jend
          DO i = istart, iend
            ztrcr(i,j,k) = zdtrcr(i,j,k) - zc(i,j,k)*ztrcr(i,j,k+1)
          ENDDO
        ENDDO
      ENDDO
      

      ! Clipping if required
      IF ( izclip(iztrcr) == T_CLP_ON ) THEN

        ! massflux correction
        IF ( lzmss_flx(iztrcr) ) THEN      
          DO j = jstart , jend
            DO i = istart , iend
              ztrcor(i,j) = 0.0_wp
            ENDDO
          ENDDO
          DO k = 1 , ke
            DO j = jstart , jend
              DO i = istart , iend
                ztrcor(i,j) = ztrcr(i,j,k) + ztrcor(i,j)*zeddpq(i,j,k)
                IF( ztrcor(i,j) < 0.0_wp ) THEN
                  ztrcr(i,j,k) = 0.0_wp
                  ztrcor(i,j)      = ztrcor(i,j)/zeddpq(i,j,k)
                ELSE
                  ztrcr(i,j,k) = ztrcor(i,j)
                  ztrcor(i,j)    = 0.0_wp
                ENDIF
              ENDDO
            ENDDO
          ENDDO

        ! call normal clipping routine
        ELSE
          CALL clipping(ztrcr(:,:,:), ie, je, ke)  
        ENDIF

      ENDIF

! fuo: this does not seem to work, although something similar
!      is done in the RK core (not important as it is only active
!      for izsp_adv_lf = 3)
!
!    ELSE
!
!      ! no vertical diffusion, simply add the tendency
!      DO k = 1, ke
!        DO j = jstart, jend
!          DO i = istart, iend
!            ztrcr(i,j,k) = ztrcr_t1(i,j,k) + dt2 * ztrcr_tens(i,j,k)
!          ENDDO
!        ENDDO
!      ENDDO

    ENDIF  

  ENDDO ! loop over tracers

!-------------------------------------------------------------------------------
! Section 3: Calculation of the surface moisture flux 'qvsflx'
!            The latent heat flux is integrated in time
!-------------------------------------------------------------------------------

!MR: Only to be done for lvertdiff

  IF (lvertdiff) THEN
    DO j =jstart, jend
      DO i = istart , iend
        qvsflx(i,j)  =  - ztch(i,j)*gr                                    &
                       * (za2t_surf*(qv_s(i,j,nold) - qv_old(i,j,ke))  +  &
                          za1t_surf*(qv_s(i,j,nnew) - qv_new(i,j,ke)))

        lhfl_s(i,j) = lh_v*qvsflx(i,j)

        ! Account for lake ice (FLake is used AND lake point AND lake ice exists)
        IF (llake) THEN
          IF( (depth_lk(i,j) > 0.0_wp) .AND.                          &
              (h_ice(i,j,nnow) >= h_Ice_min_flk) ) THEN
            ! Add latent heat of fusion
            lhfl_s(i,j) = lhfl_s(i,j) + lh_f * qvsflx(i,j)

            ! No change otherwise, but is that true for sea ice?
          END IF
!MR:
! The inclusion of a freezing rate is not consistent, as the real latent heat flux
! depends on the freezed surface fraction!
        ENDIF

      ENDDO
    ENDDO
  ENDIF

!------------------------------------------------------------------------------
! Section 4: Setup of tridiagonal matrix systems resulting from the implicit
!            numerical formulation of advection and diffusion.
!            -  pressure perturbation (no vertical diffusion)
!------------------------------------------------------------------------------ 
 
  ! Top layer       
  DO j = jstart, jend
    DO i = istart, iend
      zcg        = zgcvx(i,j,1)*zbetp
      zcs        = zgcvx(i,j,1)*zbetm
      zbg        = ed2dt - zcg
      zd1g       = ed2dt*pp(i,j,1,nold) + pptens(i,j,1)   &
                  - zcs * ( pp(i,j,2,nold) - pp(i,j,1,nold) )
      zd1(i,j,1) = zd1g/zbg
      zc(i,j,1)  = zcg /zbg
    ENDDO
  ENDDO

  ! The layers from k=2 to k=ke-1
  DO k = 2, ke-1
    DO j = jstart, jend
      DO i = istart, iend
        zag        = zgavx(i,j,k)*zbetp
        zas        = zgavx(i,j,k)*zbetm
        zcg        = zgcvx(i,j,k)*zbetp
        zcs        = zgcvx(i,j,k)*zbetm
        zbg        = ed2dt - zag - zcg
        zdg        = ed2dt * pp(i,j,k,nold) + pptens(i,j,k)  &
                   - zas*(pp(i,j,k-1,nold)-pp(i,j,k,nold))  &
                   - zcs*(pp(i,j,k+1,nold)-pp(i,j,k,nold))
        zz         = 1.0_wp/( zbg - zag*zc(i,j,k-1) )
        zc (i,j,k) = zcg*zz
        zd1(i,j,k) = ( zdg -zag*zd1(i,j,k-1) )*zz
      ENDDO
    ENDDO
  ENDDO

  ! The bottom layer
  DO j = jstart, jend
    DO i = istart, iend
      zag        = zgavx(i,j,ke)*zbetp 
      zas        = zgavx(i,j,ke)*zbetm 
      zbg        = ed2dt - zag
      zdg        = ed2dt * pp(i,j,ke,nold) + pptens(i,j,ke)         &
                   - zas *( pp(i,j,ke-1,nold) - pp(i,j,ke,nold) ) 
      znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
      pptens(i,j,ke) = ( znew - pp(i,j,ke,nold) ) * ed2dt
      ze    (i,j,ke) = znew
    ENDDO
  ENDDO

  ! Backsubstitution and storage of the complete slow tendencies

  DO k = ke-1, 1, -1
    DO j = jstart, jend
      DO i = istart, iend
        ze    (i,j,k) =   zd1(i,j,k) - zc(i,j,k)*ze(i,j,k+1)
        pptens(i,j,k) = ( ze(i,j,k) - pp(i,j,k,nold) ) * ed2dt
      ENDDO
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section 5: Setup of tridiagonal matrix systems resulting from the implicit
!            numerical formulation of advection and diffusion.
!            -  temperature
!
!------------------------------------------------------------------------------ 

  ! Top layer       
  DO j = jstart, jend
    DO i = istart, iend
      zgct       = - zkh(i,j,2)*zeddpq(i,j,1)
!     zgct       = - zkh(i,j,2)*zeddpq(i,j,1)*zcpanf(i,j,2)/zcpa(i,j,1)

      zcg        =   zgcvx(i,j,1)*zbetp &
                    + zgct*za1t(2)*zpianf(i,j,2)/zpia(i,j,2)
      zbg        = ed2dt - zgcvx(i,j,1)*zbetp &
                   - zgct*za1t(2)*zpianf(i,j,2)/zpia(i,j,1)
      zdg        = ed2dt*t(i,j,1,nold) + ttens(i,j,1) - zgct*za2t(2)* &
                   ( t(i,j,2,nold)*zpianf(i,j,2)/zpia(i,j,2)         &
                   - t(i,j,1,nold)*zpianf(i,j,2)/zpia(i,j,1) )       &
                   - zgcvx(i,j,1)*zbetm*( t(i,j,2,nold) - t(i,j,1,nold) )
      zd1(i,j,1) = zdg/zbg
      zc(i,j,1)  = zcg /zbg
    ENDDO
  ENDDO    

  ! The layers from k=2 to k=ke-1
  DO k = 2, ke-1
    DO j = jstart, jend
      DO i = istart, iend
        zgat       = - zkh(i,j,k  )*zeddpq(i,j,k)
        zgct       = - zkh(i,j,k+1)*zeddpq(i,j,k)
!       zgat       = - zkh(i,j,k  )*zeddpq(i,j,k)*zcpanf(i,j,k)/zcpa(i,j,k)
!       zgct       = - zkh(i,j,k+1)*zeddpq(i,j,k)*zcpanf(i,j,k+1)/zcpa(i,j,k)

        zag        = zgavx(i,j,k)*zbetp + &
                     zgat*za1t(k  )*zpianf(i,j,k  )/zpia(i,j,k-1)
        zcg        = zgcvx(i,j,k)*zbetp + &
                     zgct*za1t(k+1)*zpianf(i,j,k+1)/zpia(i,j,k+1)
        zbg        = ed2dt - zgavx(i,j,k)*zbetp - zgcvx(i,j,k)*zbetp   &
                           - zgat*za1t(k  )*zpianf(i,j,k  )/zpia(i,j,k) &
                           - zgct*za1t(k+1)*zpianf(i,j,k+1)/zpia(i,j,k)
        zdg        = ed2dt * t(i,j,k,nold) + ttens(i,j,k) &
                     - zgat*za2t(k  ) *zpianf(i,j,k  )*      &
                (t(i,j,k-1,nold)/zpia(i,j,k-1) - t(i,j,k,nold)/zpia(i,j,k)) & 
                     - zgct*za2t(k+1) *zpianf(i,j,k+1)*      &
                (t(i,j,k+1,nold)/zpia(i,j,k+1) - t(i,j,k,nold)/zpia(i,j,k))
        zdg        = zdg                                                 &
                     - zgavx(i,j,k)*zbetm*(t(i,j,k-1,nold)-t(i,j,k,nold)) &
                     - zgcvx(i,j,k)*zbetm*(t(i,j,k+1,nold)-t(i,j,k,nold))
        zz         = 1.0_wp/( zbg - zag*zc(i,j,k-1) )
        zc (i,j,k) = zcg*zz
        zd1(i,j,k) = ( zdg -zag*zd1(i,j,k-1) )*zz
      ENDDO
    ENDDO
  ENDDO

  ! The bottom layer
  DO j = jstart, jend
    DO i = istart, iend
      zgat       = - zkh (i,j,ke)*zeddpq(i,j,ke)
      zgct       = - ztch(i,j   )*zeddpq(i,j,ke)
!     zgat       = - zkh (i,j,ke)*zeddpq(i,j,ke)*zcpanf(i,j,ke )/zcpa(i,j,ke)
!     zgct       = - ztch(i,j   )*zeddpq(i,j,ke)*zcpanf(i,j,ke1)/zcpa(i,j,ke)

      zag        = zgavx(i,j,ke)*zbetp &
                   + zgat*za1t(ke)*zpianf(i,j,ke)/zpia(i,j,ke-1)
      zbg        = ed2dt - zgavx(i,j,ke)*zbetp                      &
                     - zgat*za1t(ke )*zpianf(i,j,ke )/zpia(i,j,ke)  &
                     - zgct*za1t_surf*zpianf(i,j,ke1)/zpia(i,j,ke)
!MR: changed computation of zdg: bug fix
      zdg        = ed2dt * t(i,j,ke,nold) + ttens(i,j,ke) &
                   - zgat*za2t(ke ) * zpianf(i,j,ke )     &
                      * (t(i,j,ke-1,nold)/zpia(i,j,ke-1) - t(i,j,ke,nold)/zpia(i,j,ke)) &
                   + zgct*za2t_surf &
                         *(zpianf(i,j,ke1)*t(i,j,ke,nold)/zpia(i,j,ke) - t_g(i,j,nold)) &
                   - zgct*za1t_surf*                                     t_g(i,j,nold)
!MR
      zdg        = zdg        - zgavx(i,j,ke)*zbetm*      &
                      ( t(i,j,ke-1,nold) - t(i,j,ke,nold) )
      znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
      ttens (i,j,ke) = ( znew - t(i,j,ke,nold) ) * ed2dt
      ze    (i,j,ke) = znew
    ENDDO
  ENDDO

  ! Backsubstitution and storage of the complete slow tendencies

  DO k = ke-1, 1, -1
    DO j = jstart, jend
      DO i = istart, iend
        ze   (i,j,k) =   zd1(i,j,k) - zc(i,j,k)*ze(i,j,k+1)
        ttens(i,j,k) = ( ze(i,j,k) - t(i,j,k,nold) ) * ed2dt
      ENDDO
    ENDDO
  ENDDO

  ! Calculation of the sensible heat flux at the surface
  ! This flux is integrated in time
!MR: Only to be done for lvertdiff

  IF (lvertdiff) THEN
    DO j = jstart, jend
      DO  i  = istart , iend
        shfl_s (i,j) = - ztch(i,j)*gr*cp_d*                          &
                      (za2t_surf*(t_g(i,j,nold) - zpianf(i,j,ke1)*   &
                                  t  (i,j,ke,nold)/zpia(i,j,ke) ) +  &
                       za1t_surf*(t_g(i,j,nnew) - zpianf(i,j,ke1)*   &
                                  ze(i,j,ke)/zpia(i,j,ke) ) )
      ENDDO
    ENDDO
  ENDIF

!------------------------------------------------------------------------------
! Section 6: Setup of tridiagonal matrix systems resulting from the implicit
!            numerical formulation of advection and diffusion.
!            -  vertical velocity (diffusion omitted! )
!------------------------------------------------------------------------------ 
  DO  k = 2, ke
    ! Precalculate some help variables to avoid divisions in subsequent code
    DO j = jstart, jend
      DO i = istart, iend
        zdr  = 1.0_wp / ( dp0(i,j,k) + dp0(i,j,k-1) )
        zgav       = - ( wcon(i,j,k) + wcon(i,j,k-1) ) * zdr
        zgcv       = + ( wcon(i,j,k) + wcon(i,j,k+1) ) * zdr
        zgavx(i,j,k) = zgav*dp0(i,j,k  )*zdr
        zgcvx(i,j,k) = zgcv*dp0(i,j,k-1)*zdr
      ENDDO
    ENDDO
  ENDDO

  ! Top layer       
  DO j = jstart, jend
    DO i = istart, iend
      zcg        = zgcvx(i,j,2)*zbetp
      zbg        = ed2dt - zcg - zgavx(i,j,2)*zbetp
      zdg        = ed2dt*w(i,j,2,nold) + wtens(i,j,2)   &
                   - zgcvx(i,j,2)*zbetm*( w(i,j,3,nold) - w(i,j,2,nold) )  &
                   + zgavx(i,j,2)*zbetm*w(i,j,2,nold) &
                   - zgavx(i,j,2)*w(i,j,1,nold)
      zd1(i,j,2)  = zdg/zbg
      zc (i,j,2)  = zcg/zbg
    ENDDO
  ENDDO

  ! The layers from k=3 to k=ke-1
  DO k = 3, ke-1
    DO j = jstart, jend
      DO i = istart, iend
        zag        = zgavx(i,j,k)*zbetp
        zas        = zgavx(i,j,k)*zbetm
        zcg        = zgcvx(i,j,k)*zbetp
        zcs        = zgcvx(i,j,k)*zbetm
        zbg        = ed2dt - zag - zcg
        zdg        = ed2dt * w(i,j,k,nold) + wtens(i,j,k)            &
                   - zas * ( w(i,j,k-1,nold) - w(i,j,k,nold) )  &
                   - zcs * ( w(i,j,k+1,nold) - w(i,j,k,nold) )
        zz         = 1.0_wp/( zbg - zag*zc(i,j,k-1) )
        zc (i,j,k)  = zcg*zz
        zd1(i,j,k)  = ( zdg - zag*zd1(i,j,k-1) )*zz
      ENDDO
    ENDDO
  ENDDO

  ! The bottom layer
  DO j = jstart, jend
    DO i = istart, iend
      zag        = zgavx(i,j,ke)*zbetp 
      zbg        = ed2dt - zag - zgcvx(i,j,ke)*zbetp
      zdg        = ed2dt * w(i,j,ke,nold) + wtens(i,j,ke)         &
                   - zgavx(i,j,ke)*zbetm*( w(i,j,ke-1,nold) - w(i,j,ke,nold) ) &
                   + zgcvx(i,j,ke)*zbetm*w(i,j,ke,nold) &
                   - zgcvx(i,j,ke)*w(i,j,ke1,nold)
      znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
      wtens(i,j,ke) = ( znew - w(i,j,ke,nold) ) * ed2dt
      ze   (i,j,ke) = znew
    ENDDO
  ENDDO

  ! Backsubstitution and storage of the complete slow tendencies

  DO k = ke-1, 2, -1
    DO j = jstart, jend
      DO i = istart, iend
        ze    (i,j,k) =   zd1(i,j,k) - zc(i,j,k)*ze(i,j,k+1)
        wtens (i,j,k) = ( ze(i,j,k) - w(i,j,k,nold) ) * ed2dt
      ENDDO
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section 7: Setup of tridiagonal matrix systems resulting from the implicit
!            numerical formulation of advection and diffusion.
!            -  horizontal wind velocity u
!------------------------------------------------------------------------------ 

  ! Precalculate some help variables to avoid divisions in subsequent code
  ! zeddpq = rho0/rho is used for the u- and the v-equations
  DO  k = 1, ke
    DO j = jstart-1, jend+1
      DO i = istart-1, iend+1
        zeddpq(i,j,k) = rho0(i,j,k) / rho(i,j,k)
      ENDDO
    ENDDO
  ENDDO

  ! Top layer:   k = 1
  DO j = jstartu, jendu
    DO i = istartu, iendu
      zdp0r      = 1.0_wp / (dp0(i,j,1) + dp0(i+1,j,1))
      zeddpb     = (zeddpq(i,j,1) + zeddpq(i+1,j,1)) * zdp0r
      zkm(i,j)   = 0.5_wp*( ztmkvm(i,j,2) + ztmkvm(i+1,j,2) )
      zgct       = - zkm(i,j)*zeddpb
      zgcv       = 0.5_wp*(wcon(i,j,2)+wcon(i+1,j,2))*zdp0r
      zcg        = zgcv*zbetp + zgct*za1t(2)
      zbg        = ed2dt - zcg
      zdg        = ed2dt*u(i,j,1,nold) + utens(i,j,1)  &
                   - ( za2t(2)*zgct + zbetm*zgcv )*        &
                     ( u(i,j,2,nold) - u(i,j,1,nold) )
      zd1(i,j,1) = zdg/zbg
      zc (i,j,1) = zcg/zbg
    ENDDO
  ENDDO

  ! The layers from k=2 to k=ke-1
  DO k = 2, ke-1
    DO j = jstartu, jendu
      DO i = istartu, iendu
        zdp0r      = 1.0_wp / (dp0(i,j,k) + dp0(i+1,j,k))
        zeddpb     = (zeddpq(i,j,k) + zeddpq(i+1,j,k)) * zdp0r
        zgat       = - zkm(i,j)*zeddpb
        zkm(i,j)   = 0.5_wp*( ztmkvm(i,j,k+1) + ztmkvm(i+1,j,k+1) )
        zgct       = - zkm(i,j)*zeddpb
        zgav       = -0.5_wp*(wcon(i,j,k  )+wcon(i+1,j,k  ))*zdp0r
        zgcv       =  0.5_wp*(wcon(i,j,k+1)+wcon(i+1,j,k+1))*zdp0r
        zag        = zgav*zbetp + zgat*za1t(k)
        zas        = zgav*zbetm + zgat*za2t(k)
        zcg        = zgcv*zbetp + zgct*za1t(k+1)
        zcs        = zgcv*zbetm + zgct*za2t(k+1)
        zbg        = ed2dt - zag - zcg
        zdg        = ed2dt*u(i,j,k,nold) + utens(i,j,k)  &
                     - zas * (u(i,j,k-1,nold) - u(i,j,k,nold))  &
                     - zcs * (u(i,j,k+1,nold) - u(i,j,k,nold))
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
      zdp0r      = 1.0_wp / (dp0(i,j,ke) + dp0(i+1,j,ke))
      zeddpb     = (zeddpq(i,j,ke) + zeddpq(i+1,j,ke)) * zdp0r
      zgat       = - zkm(i,j) *zeddpb
      ztmcmq     = 0.5_wp*( ztcm(i,j) + ztcm(i+1,j) )
      zgct       = - ztmcmq*zeddpb
      zgav       = - 0.5_wp*(wcon(i,j,ke)+wcon(i+1,j,ke))*zdp0r
      zag        = zgav*zbetp + zgat*za1t(ke)
      zas        = zgav*zbetm + zgat*za2t(ke)
      zcg        = za1t(ke1)*zgct
      zcs        = za2t(ke1)*zgct
      zbg        = ed2dt - zag - zcg
      zdg        = ed2dt*u(i,j,ke,nold) + utens(i,j,ke)  &
                   - zas * ( u(i,j,ke-1,nold) - u(i,j,ke,nold) ) &
                   + zcs * u(i,j,ke,nold)
      znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
      utens (i,j,ke) = ( znew - u(i,j,ke,nold) ) * ed2dt
!MR: compute umfl_s only for lvertdiff below
!MR   umfl_s(i,j   ) = ztmcmq*gr*( a2t(ke1)*u(i,j,ke,nold) + a1t(ke1)*znew )
      ze    (i,j,ke) = znew
    ENDDO
  ENDDO

!MR:
  IF (lvertdiff) THEN
    DO j = jstartu, jendu
      DO i = istartu, iendu
       umfl_s(i,j) = zkm(i,j)*( a2t(ke1)*u(i,j,ke,nold) + a1t(ke1)*ze(i,j,ke) )
      ENDDO
    ENDDO
  ENDIF

  ! Backsubstitution and storage of the complete slow tendencies

  DO k = ke-1, 1, -1
    DO j = jstartu, jendu
      DO i = istartu, iendu
        ze   (i,j,k) =   zd1(i,j,k) - zc(i,j,k)*ze(i,j,k+1)
        utens(i,j,k) = ( ze (i,j,k) - u(i,j,k,nold) ) * ed2dt
      ENDDO
    ENDDO
  ENDDO


!------------------------------------------------------------------------------
! Section 8: Setup of tridiagonal matrix systems resulting from the implicit
!            numerical formulation of advection and diffusion.
!            -  horizontal wind velocity v
!------------------------------------------------------------------------------ 

  ! Top layer  k=1
  DO j = jstartv, jendv
    DO i = istartv, iendv
      zdp0r      = 1.0_wp / (dp0(i,j,1)+dp0(i,j+1,1))
      zeddpb     = (zeddpq(i,j,1) + zeddpq(i,j+1,1)) * zdp0r
      zkm(i,j)   = 0.5_wp*( ztmkvm(i,j,2) + ztmkvm(i,j+1,2) )
      zgct       = - zkm(i,j)*zeddpb
      zgcv       = 0.5_wp*(wcon(i,j,2)+wcon(i,j+1,2))*zdp0r
      zcg        = zgcv*zbetp + zgct*za1t(2)
      zbg        = ed2dt - zcg
      zdg        = ed2dt*v(i,j,1,nold) + vtens(i,j,1)  &
                   - ( za2t(2)*zgct + zbetm*zgcv )*        &
                     ( v(i,j,2,nold) - v(i,j,1,nold) )
      zd1(i,j,1) = zdg/zbg
      zc (i,j,1) = zcg/zbg
    ENDDO
  ENDDO

  ! The layers from k=2 to k=ke-1
  DO  k = 2, ke-1
    DO j = jstartv, jendv
      DO i = istartv, iendv
        zdp0r      = 1.0_wp / (dp0(i,j,k)+dp0(i,j+1,k))
        zeddpb     = (zeddpq(i,j,k) + zeddpq(i,j+1,k)) * zdp0r
        zgat       = - zkm(i,j)*zeddpb
        zkm(i,j)   = 0.5_wp*( ztmkvm(i,j,k+1) + ztmkvm(i,j+1,k+1) )
        zgct       = - zkm(i,j)*zeddpb
        zgav       = -0.5_wp*(wcon(i,j,k  )+wcon(i,j+1,k  ))*zdp0r
        zgcv       =  0.5_wp*(wcon(i,j,k+1)+wcon(i,j+1,k+1))*zdp0r
        zag        = zgav*zbetp + zgat*za1t(k)
        zas        = zgav*zbetm + zgat*za2t(k)
        zcg        = zgcv*zbetp + zgct*za1t(k+1)
        zcs        = zgcv*zbetm + zgct*za2t(k+1)
        zbg        = ed2dt - zag - zcg
        zdg        = ed2dt*v(i,j,k,nold) + vtens(i,j,k)  &
                     - zas * (v(i,j,k-1,nold) - v(i,j,k,nold))  &
                     - zcs * (v(i,j,k+1,nold) - v(i,j,k,nold))
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
      zdp0r      = 1.0_wp / (dp0(i,j,ke)+dp0(i,j+1,ke))
      zeddpb     = (zeddpq(i,j,ke) + zeddpq(i,j+1,ke)) * zdp0r
      zgat       = - zkm(i,j) *zeddpb
      ztmcmq     = 0.5_wp*( ztcm(i,j) + ztcm(i,j+1) )
      zgct       = - ztmcmq*zeddpb
      zgav       = - 0.5_wp*(wcon(i,j,ke)+wcon(i,j+1,ke))*zdp0r
      zag        = zgav*zbetp + zgat*za1t(ke)
      zas        = zgav*zbetm + zgat*za2t(ke)
      zcg        = za1t(ke1)*zgct
      zcs        = za2t(ke1)*zgct
      zbg        = ed2dt - zag - zcg
      zdg        = ed2dt*v(i,j,ke,nold) + vtens(i,j,ke)  &
                   - zas * ( v(i,j,ke-1,nold) - v(i,j,ke,nold) ) &
                   + zcs * v(i,j,ke,nold)
      znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
      vtens (i,j,ke) = ( znew - v(i,j,ke,nold) ) * ed2dt
!MR: compute vmfl_s only for lvertdiff below
!MR   vmfl_s(i,j   ) = ztmcmq*gr*( a2t(ke1)*v(i,j,ke,nold) + a1t(ke1)*znew )
      ze    (i,j,ke) = znew
    ENDDO
  ENDDO

!MR
  IF (lvertdiff) THEN
    DO j = jstartv, jendv
      DO i = istartv, iendv
        vmfl_s(i,j) = zkm(i,j)*( a2t(ke1)*v(i,j,ke,nold) + a1t(ke1)*ze(i,j,ke) )
      ENDDO
    ENDDO
  ENDIF

  ! Backsubstitution and storage of the complete slow tendencies

  DO k = ke-1, 1, -1
    DO j = jstartv, jendv
      DO i = istartv, iendv
        ze   (i,j,k) =   zd1(i,j,k) - zc(i,j,k)*ze(i,j,k+1)
        vtens(i,j,k) = ( ze(i,j,k) - v(i,j,k,nold) ) * ed2dt
      ENDDO
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! End of subroutine slow_tendencies
!------------------------------------------------------------------------------

END SUBROUTINE slow_tendencies
