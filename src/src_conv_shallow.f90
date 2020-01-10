!+ Source module  "src_conv_shallow"
!------------------------------------------------------------------------------

MODULE src_conv_shallow

!------------------------------------------------------------------------------
!
! Description:
!   The module "src_conv_shallow" performs calculations related to the 
!   parameterization of subgrid-scale moist convection for use in 
!   high-resolution runs on the meso-gamma scale (LMK). The present 
!   approach is based on a simple mass-flux scheme with moisture convergence
!   closure.
!
!   All global variables of the model that are used by the convection routines
!   are imported by USE statements below. The interface of the convection
!   routines and the model is provided by the organizational routine
!   "organize_conv_shallow".
!
! Current Code Owner: DWD, Dmitrii Mironov
!  phone:  +49  69  8062 2705
!  fax:    +49  69  8062 3721
!  email:  Dmitrii.Mironov@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 3.16       2005/07/22 Jochen Foerstner
!  Initial release
! 3.18       2006/03/03 Jochen Foerstner
!  Corrections for writing instantaneous values (Jochen Foerstner)
! V4_5         2008/09/10 Ulrich Schaettler
!  Moved declaration of entr_sc (before: entrscv) to new module data_convection
! V4_8         2009/02/16 Guenther Zaengl
!  Use p0hl (reference pressure at half levels) for full consistency with
!  new reference atmosphere implementation
! V4_12        2010/05/11 Ulrich Schaettler
!  Removed t0(_melt), qvsflx
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_18        2011/05/26 Martin Koehler
!  Replaced hard coded value by Namelist variable: thick_sc
! V4_20        2011/08/31 Matthias Raschendorfer
!  Introducing calculation of convective buoyant TKE production 'tket_conv'
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer
!  Replaced qx-variables by using them from the tracer module
! V4_26        2012/12/06 Matthias Raschendorfer
!  Allowing only non negative buoyant production terms of TKE
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  Introduced MESSy interface
! V4_28        2013/07/12 Ulrich Schaettler
!  Eliminated reference to sigmr (not used)
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
! V5_1         2014-11-28 Matthias Raschendorfer, Oliver Fuhrer, Michael Baldauf
!  Calculation of 'tket_conv' corrected.
!  Replaced ireals by wp (working precision) (OF)
!  Removed option lexpl_lbc=.FALSE.  (MB)
! V5_3         2015-10-09 Ulrich Blahak
!  Added computation of ttdiab_conv, the pure diabatic tendency due to convection
! V5_4b        2016-07-12 Ulrich Schaettler
!  Removed dp0 from the USE lists (not used here)
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

USE data_convection, ONLY :   &
    entr_sc,      & ! mean entrainment rate for shallow convection
    thick_sc        ! limit for convective clouds to be "shallow" (in Pa)

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------
    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction
    ke1,          & ! ke + 1
    ieke,         & ! ie*ke

! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the 
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from 
!    the other ones because of the use of the staggered Arakawa-C-grid.
!    
    istart   ,    & !
    iend     ,    & !
    jstart   ,    & !
    jend     ,    & !
    istartpar,    & ! start index for computations in the parallel program
    iendpar,      & ! end index for computations in the parallel program
    jstartpar,    & ! start index for computations in the parallel program
    jendpar,      & ! end index for computations in the parallel program


! 4. constants for the horizontal rotated grid and related variables
! ------------------------------------------------------------------
    edadlat,      & ! 1 / (radius of the earth * dlat)

! 5. variables for the time discretization and related variables
! --------------------------------------------------------------
    dt,           & ! long time-step
    dt2,          & ! dt*2.            

! 8. Organizational variables to handle the COSMO humidity tracers
! ----------------------------------------------------------------
    idt_qv

! end of data_modelconfig

!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 2. physical constants and related variables
! -------------------------------------------

    r_d,          & ! gas constant for dry air    
    r_v,          & ! gas constant of water vapour
    rdv,          & ! r_d / r_v
    rvd_m_o,      & ! r_v/r_d - 1  
    cp_d,         & ! specific heat of dry air at constant pressure
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
    b4i             !               -- " --

! end of data_constants

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------
    p0         ,    & ! base state pressure                           (Pa) 
    p0hl       ,    & ! base state pressure on half levels            (Pa)
    hhl        ,    & ! geometrical height of half levels             ( m ) 

! 2. external parameter fields                                        (unit)
! ----------------------------
    rmy        ,    & ! Davis-parameter for boundary relaxation         --


! 3. prognostic variables                                             (unit)
! -----------------------
    u          ,    & ! zonal wind speed                              ( m/s )
    v          ,    & ! meridional wind speed                         ( m/s )
    w          ,    & ! vertical velocity                             ( m/s )
    t          ,    & ! temperature                                   (  k  )
    pp         ,    & ! deviation from the reference pressure         ( pa  )

! 5. fields for surface values and soil model variables               (unit )
! -----------------------------------------------------
    ps        ,     & ! surface pressure                              ( pa  )
    t_g       ,     & ! weighted surface temperature                  (  K  )

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------
!   fields for convective subgrid-scale precipitation
    clc_con     ,   & ! cloud cover due to convection                   --
    clw_con     ,   & ! cloud liquid water due to convection            --
    bas_con     ,   & ! level index of convective cloud base            -- 
    top_con     ,   & ! level index of convective cloud base            --
    tt_conv     ,   & ! temperature tendency due to convection        ( K/s  )
    ttdiab_conv ,   & ! pure diabatic temperature tendency due to convection ( K/s  )
    qvt_conv    ,   & ! humidity    tendency due to convection        ( 1/s  )
    tket_conv   ,   & ! TKE-tendency due to convective buoyancy       ( m2/s3)
    dqvdt       ,   & ! threedimendional moisture convergence         (1/s)
    mflx_con    ,   & ! cloud base massflux                           (kg/m2*s)
    cape_con    ,   & ! convective available energy                   (   J/kg)
    qcvg_con          ! moisture convergence for Kuo-type closure     (    1/s)

! end of data_fields

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
    nold,         & ! corresponds to ntstep - 1
    nnow,         & ! corresponds to ntstep

! 3. controlling the physics
! --------------------------
    lconf_avg,    & ! average convective forcings in case of massflux closure
    lcape,        & ! convection with CAPE closure
    lconv_inst,   & ! output of instantaneous values of top_con/bas_con
                    ! instead of min/max for an output interval

! 5. additional control variables
! -------------------------------
    l2tls,        & ! forecast with 2-TL integration scheme
    loutput_diab    ! internal switch for output of temperature tendencies due to pure diabatic processes

! end of data_runcontrol 

!------------------------------------------------------------------------------

USE data_parallel , ONLY:    &
  my_cart_id         ! rank of this subdomain in the cartesian communicator
                     ! that can be used by MPI_WAIT to identify the send

!------------------------------------------------------------------------------

USE environment   , ONLY: model_abort

!------------------------------------------------------------------------------

USE src_tracer,     ONLY: trcr_get, trcr_errorstr

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Declarations

! The following parameters are tunable constants for the cumulus convection
! scheme. 

LOGICAL, PRIVATE ::       &
  lmfmid        = .TRUE.     ! switch for mid-level convection
  !lmfmid        = .FALSE.     ! switch for mid-level convection

REAL (KIND = wp),     PARAMETER, PRIVATE  :: &
  entrmid       = 0.00010_wp, & ! mean entrainment rate for mid-level convection

  cmfcmax       = 1.0_wp    , & ! maximum mass flux
  cmfcmin       = 1.E-10_wp , & ! minimum mass flux  (security)
  cmfctop       = 0.33_wp       ! relative mass flux above level of non-buoyancy

REAL (KIND = wp),     PRIVATE  :: &
    ctmelt               , & ! tripel point
    cc2                  , & ! definition of utitility constants
    c5hlccp              , & ! for saturation humidity
    c5hlscp              , & !
    chlcdcp              , & !
    chlsdcp                  !

! The logical variable LFIRSTI is used to initialize additional parameters     
! on the first CALL of Subroutine organize_conv_tiedtke
LOGICAL, PRIVATE ::       &
  lfirsti = .TRUE.           ! switch for initialization

!==============================================================================
! Module procedures in src_conv_shallow
!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure for organization
!------------------------------------------------------------------------------

SUBROUTINE organize_conv_shallow

!------------------------------------------------------------------------------
!
! Description:
!
!   The module procedure organize_conv_tiedtke is the interface of the model
!   to the parameterization package for moist convection.
!   At present, only one parameterization scheme, the Tiedtke mass-flux scheme,
!   is available.
!
! Externals
!   cu_shallow  : Shallow-convection mass-flux scheme             
!
!------------------------------------------------------------------------------

! Local scalars and automatic arrays (also to be used in lower level routines):
! -----------------------------------------------------------------------------

  ! Input for the convection routine "cu_shallow"
  ! -----------------------------------------
  REAL    (KIND=wp   )      ::  &
     zt      (ie,ke),    & ! temperature at full levels
     zqv     (ie,ke),    & ! specific humidiy at full levels
     zu      (ie,ke),    & ! zonal wind component
     zv      (ie,ke),    & ! meridional wind component
     zw      (ie,ke),    & ! vertical velocity in z-system
     zfif    (ie,ke  ),  & ! geopotential at full levels
     zfih    (ie,ke1),   & ! geopotential at half levels
     zpf     (ie,ke  ),  & ! pressure at full levels 
     zph     (ie,ke1),   & ! pressure at half levels 
     zt_g    (ie    ),   & ! surface temperature
     zdqvdt  (ie,ke),    & ! moisture tencency
     zcumflx (ie),       & ! cu_base massflux
     zcucape (ie),       & ! convective available energy
     zcuconv (ie),       & ! moisture convergence to be used in cu_base massflux
     zdt                   ! actual timestep for 2 or 3 TL-scheme     

  ! Output from the convection routine "cu_shallow"
  ! -----------------------------------------
  REAL    (KIND=wp   )     ::  &
     zmyfac   (ie),     & ! Davies-Relaxation factor
     zdt_con  (ie,ke),  & ! convective tendency of temperature
     zdtdiab_con(ie,ke),& ! conv T-tend. due to LH exchanges
     zdqv_con (ie,ke),  & ! convective tendency of specific humidity
     zdtke_con(ie,ke),  & ! convective buoyant TKE production on half levels
     zco_clw  (ie,ke)     ! convective cloud water

  INTEGER (KIND=iintegers) ::  &
     mbas_con (ie)      , & ! cloud base level index
     mtop_con (ie)          ! cloud top  level index

  LOGICAL                  ::  &
     locum (ie)             ! indicator for convection at gridpoints
  
  ! Other local arrays and scalars  
  REAL    (KIND=wp   )     ::  &
     zbas, ztop,            &
     zdttop, zdtbot, zdepth_max, zdepth_now

 REAL    (KIND=wp   )     ::  &  ! Weights for horizontal averaging
     zcent = 0.2500_wp, & ! centre weight in a nine point stencil
     zside = 0.1250_wp, & ! weight for side points 
     zedge = 0.0625_wp    ! weight for edge points

  INTEGER (KIND=iintegers) ::  &
     i, j, k, km1, nx,      & !
     mtop, mbas

! For error handling
! ------------------
  INTEGER (KIND=iintegers) ::  &
    izerror

  CHARACTER (LEN=255)      ::  &
    yzerrmsg

  CHARACTER (LEN=25)       ::  &
    yzroutine

! Tracer pointers
!----------------
  REAL (KIND=wp),     POINTER :: &
    qv  (:,:,:) => NULL()            ! QV at nx

!
! End of header
!==============================================================================

  izerror  = 0
  yzerrmsg = '   '
  yzroutine= 'organize_conv_shallow'

  ! Select timelevel and timestep of the computation
  IF ( l2tls ) THEN
    nx  = nnow
    zdt = dt
  ELSE
    nx  = nold
    zdt = dt2
  ENDIF


  ! Some initializations
  IF (lfirsti) THEN
    ctmelt  = 273.16_wp                ! tripel point
    cc2     = b1 * rdv                 ! definition of utitility constants
    c5hlccp = b2w*(b3-b4w)*lh_v/cp_d   ! for saturation humidity
    c5hlscp = b2i*(b3-b4i)*lh_s/cp_d   !
    chlcdcp = lh_v/cp_d                !
    chlsdcp = lh_s/cp_d                !
    lfirsti = .FALSE.
  ENDIF

  ! Retrieve the required microphysics tracers
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nx, ptr = qv)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  ! In order to save come CPU-time, the convection scheme in a former version
  ! of the model has been called only for those points, which have been 
  ! identified as convecively unstable in the previous time step. The present
  ! version calls the convection routine only at fixed time increments nincconv
  ! (e.g. every 10 timsteps) and stores the concective tendencies and the 
  ! precipitation rates on global arrays which are held fixed in the 
  ! intermediate steps.  The time increment nincconv can be set on NAMELIST  
  ! input.

  ! Reset the  convective cloud cover each time when the convection scheme 
  ! is called

  clc_con (:,:,:) = 0.0_wp

 
  ! In the present version of the model, the vertical layer index of the
  ! convective cloud base and cloud top are stored as hourly maximum
  ! values and after each hour of integration these fields are reset
  ! to zero in routine "near_surface".
  ! the arrays top_con and bas_con are reset in routine "near_surface"
  ! or here for instantaneous values.
  IF ( lconv_inst ) THEN
    bas_con (:,:) = 0.0_wp
    top_con (:,:) = 0.0_wp
  ENDIF

  !--------------------------------------------------------------
  ! Loop from south to north. The convection scheme is called for
  ! eache j-slice.
  !--------------------------------------------------------------

  DO  j = jstart, jend

    ! Preset the output arrays of the convection scheme
    zdt_con  (:,:) = 0.0_wp
    zdtdiab_con(:,:)= 0.0_wp
    zdqv_con (:,:) = 0.0_wp
    zdtke_con(:,:) = 0.0_wp
    zco_clw  (:,:) = 0.0_wp
    mbas_con (:)   = 0
    mtop_con (:)   = 0
    locum    (:)   = .FALSE.

    ! Preset the additional output fields for convection
    zcumflx(:) = 0.0_wp
    zcucape(:) = 0.0_wp
    zcuconv(:) = 0.0_wp

    ! Prepare the input arrays for the convection scheme
    DO i = istart, iend
      zmyfac(i) = 1.0_wp - rmy(i,j,1)
    ENDDO

    DO  k = 1, ke 
      km1 = MAX ( 1, k-1 )
      DO i = istart, iend
        zt    (i,k)  = t (i,j,k,nx)
        zqv   (i,k)  = qv(i,j,k   )
        zu    (i,k)  = 0.5_wp*( u(i,j,k,nx) + u(i-1,j,k,nx) )
        zv    (i,k)  = 0.5_wp*( v(i,j,k,nx) + v(i,j-1,k,nx) )
        zw    (i,k)  = 0.5_wp*( w(i,j,k,nx) + w(i,j,k+1,nx) )
        zfif  (i,k)  = 0.5_wp*g*( hhl(i,j,k) + hhl(i,j,k+1) )
        zfih  (i,k)  = g* hhl(i,j,k) 
        zdqvdt(i,k)  = dqvdt(i,j,k)*zmyfac(i)
        zpf   (i,k)  = p0(i,j,k) + pp(i,j,k,nx)
        zph   (i,k)  = p0hl(i,j,k) + 0.5_wp*(pp(i,j,k,nx) + pp(i,j,km1,nx))
      ENDDO
    ENDDO
    DO i = istart, iend
      zfih  (i,ke1)  = g* hhl(i,j,ke1) 
      zph   (i,ke1)  = ps(i,j,nx)
      zt_g  (i)      = t_g  (i,j,nx)
    ENDDO     

    IF (lconf_avg) THEN ! Replace convective forcings by area average
      DO  k = 1, ke
        km1 = MAX ( 1, k-1 )
        DO i = istart, iend
          zw    (i,k)  = 0.5_wp * &
                         (  zcent*( w(i,j,k,nx) + w(i,j,k+1,nx)  )           &
                          + zside*( w(i-1,j  ,k  ,nx) + w(i+1,j  ,k  ,nx)    &
                                  + w(i,j-1  ,k  ,nx) + w(i  ,j+1,k  ,nx)  ) &
                          + zedge*( w(i-1,j-1,k  ,nx) + w(i+1,j-1,k  ,nx)    &
                                  + w(i-1,j+1,k  ,nx) + w(i+1,j+1,k  ,nx)  ) &
                          + zside*( w(i-1,j  ,k+1,nx) + w(i+1,j  ,k+1,nx)    &
                                  + w(i  ,j-1,k+1,nx) + w(i  ,j+1,k+1,nx)  ) &
                          + zedge*( w(i-1,j-1,k+1,nx) + w(i+1,j-1,k+1,nx)    &
                                  + w(i-1,j+1,k+1,nx) + w(i+1,j+1,k+1,nx)  ) )
          zdqvdt(i,k)  = zmyfac(i) * ( zcent*dqvdt(i,j,k)                  &
                          + zside*( dqvdt(i-1,j,k  ) + dqvdt(i+1,j,k  )    &
                                  + dqvdt(i,j-1,k  ) + dqvdt(i,j+1,k  )  ) &
                          + zedge*( dqvdt(i-1,j-1,k) + dqvdt(i+1,j-1,k)    &
                                  + dqvdt(i-1,j+1,k) + dqvdt(i+1,j+1,k)  ) )
        ENDDO
      ENDDO
    ENDIF
 
    CALL cu_shallow (                                          &
         zt     , zqv     , zu      , zv    , zw    ,          &
         zfif   , zfih    , zpf     , zph   ,                  &
         zt_g   , zdqvdt  ,                                    &
         ie     , ke      , istart  , iend  ,                  &
         zcumflx, zcuconv , zco_clw,                           &
         zdt_con, zdtdiab_con, zdqv_con, zdtke_con, mbas_con,  &
         mtop_con, locum)

    ! Store the output from the convection scheme on the corresponding
    ! global arrays
 
    DO  k = 1, ke
      DO  i = istart, iend
        IF( locum(i) ) THEN
          tt_conv  (i,j,k) = zdt_con (i,k) 
          qvt_conv (i,j,k) = zdqv_con(i,k)
          tket_conv(i,j,k) = zdtke_con(i,k)
          clw_con  (i,j,k) = zco_clw (i,k)
        ELSE
          tt_conv  (i,j,k) = 0.0_wp
          qvt_conv (i,j,k) = 0.0_wp
          tket_conv(i,j,k) = 0.0_wp
          clw_con  (i,j,k) = 0.0_wp
        ENDIF
      ENDDO     
    ENDDO     
    IF (loutput_diab) THEN
      DO  k = 1, ke
        DO  i = istart, iend
          IF (locum(i)) THEN
            ttdiab_conv(i,j,k) = zdtdiab_con (i,k)
          ELSE
            ttdiab_conv(i,j,k) = 0.0_wp
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    DO  i = istart, iend
      IF( locum(i) .AND. mtop_con(i) > 0 ) THEN
        mtop_con (i) = MAX   ( mtop_con(i)-2, 2 )
        zdepth_max   = bas_con(i,j) - top_con(i,j)
        zdepth_now   = REAL ( mbas_con(i) - mtop_con(i), wp )
        IF ( zdepth_max < zdepth_now ) THEN
          top_con(i,j) = REAL ( mtop_con(i), wp )
          bas_con(i,j) = REAL ( mbas_con(i), wp )
        ENDIF
      ENDIF
    ENDDO     

    ! Calculate the convective cloud cover by a simple empirical
    ! relation (following B.Ritter, FE14). Anvils are assumed for
    ! a temperature increase at top level
    DO  i = istart, iend
      IF( locum(i) .AND. mtop_con(i) > 0 ) THEN
        mtop = mtop_con(i)
        mbas = mbas_con(i)
        zbas = 0.5_wp*( hhl(i,j,mbas) + hhl(i,j,mbas+1) )
        ztop = 0.5_wp*( hhl(i,j,mtop) + hhl(i,j,mtop+1) )
        DO  k = mtop, mbas-1
          clc_con(i,j,k) = 0.35_wp*(ztop-zbas)/5000.0_wp 
          IF ( k == mtop ) THEN
            zdtbot = t(i,j,k+1,nx) - t(i,j,k  ,nx)
            zdttop = t(i,j,k  ,nx) - t(i,j,k-1,nx)
            IF ( zdtbot > 0.0_wp .AND. zdttop <= 0.0_wp ) THEN
              clc_con(i,j,k) = 2.0_wp*clc_con(i,j,k)
            ENDIF
          ENDIF
          clc_con(i,j,k) = MIN ( 1.0_wp, MAX(0.05_wp, clc_con(i,j,k)) )
        ENDDO        
      ENDIF
    ENDDO        
 
    DO i = istart, iend
      IF (locum(i)) THEN
        mflx_con(i,j) = zcumflx(i)
        cape_con(i,j) = zcucape(i)
        qcvg_con(i,j) = zcuconv(i)
      ELSE
        mflx_con(i,j) = 0.0_wp
        cape_con(i,j) = 0.0_wp
        qcvg_con(i,j) = 0.0_wp
      END IF
    END DO

  !--------------------------------------------------------------
  ! End of the loop from south to north. 
  !--------------------------------------------------------------
  ENDDO       

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------
END SUBROUTINE organize_conv_shallow

!==============================================================================

!+ Module procedure in "Convection" 
!------------------------------------------------------------------------------

SUBROUTINE cu_shallow (                                          &
           pt     , pqv     , pu      , pv    , pw    ,          &
           pfif   , pfih    , ppf     , pph   ,                  &
           pt_g   , pdqvdt  ,                                    &
           idim   , kdim    , isc     , iec   ,                  &
           pcumflx, pcuconv , pco_clw,                           &
           pdt_con, pdtdiab_con, pdqv_con, pdtke_con,            &
           mbas_con, mtop_con,  locum)

!------------------------------------------------------------------------------
!
! Description:
!
!   The module procedure cu_shallow organizes the massflux (shallow) cumulus 
!   parameterization scheme.
!   This routine computes the physical tendencies of the prognostic variables 
!   t and qv due to convective processes. These are calculated from
!     - convective fluxes due to updrafts
!     - condensation of cloud water in the updrafts 
!     - evaporation of cloud water in the environment
!
!   Input for the scheme are the grid scale values of T, qv, qc, w, p, fi 
!   and the moisture tendency due to turbulent mixing and 3-d advection 
! 
!   Output of the scheme are 
!     - tendencies of T and qv
!     - cloud base and cloud top model level indices and hights
!
!
! Method:
!
!   The parameterisation is based on a mass flux representation of
!   convective processes and proceeds along the following steps:
!
!     - definition of constants and parameters
!     - specification of half level values and initialization of
!       updraft values 
!     - determination of cloud base (if existing) 
!       and specification of cloud base mass flux from PBL moisture
!       budget (lcape = .false.), or from convective available
!       energy (lcape = .true. )
!     - cloud ascent calculations (no downdrafts are considered)
!     - final adjustments to convective fluxes 
!     - calculation of the tendencies of T and qv 
!
! Externals
!
!   cu_asc  : cloud ascent for entraining plume
!   cu_cond : condensation/evaporation calculations
!
! Switches 
!
!   lmfmid = .T.  MIDLEVEL CONVECTION IS SWITCHED ON
!
! Model parameters 
!
!   Have been defined within the module
!
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

! Input data
! ----------
  INTEGER (KIND=iintegers), INTENT (IN) ::  &
     idim ,       & ! array dimension in zonal direction
     kdim ,       & ! array dimension in vertical direction 
     isc  ,       & ! start index for first  array computation
     iec            ! end   index for first  array computation

  REAL    (KIND=wp   ),     INTENT (IN) ::  &
     pt      (idim,kdim),     & ! temperature at full levels
     pqv     (idim,kdim),     & ! specific humidiy at full levels
     pu      (idim,kdim),     & ! zonal wind component
     pv      (idim,kdim),     & ! meridional wind component
     pw      (idim,kdim),     & ! vertical velocity in z-system (full levles)
     pfif    (idim,kdim  ),   & ! geopotential at full levels
     pfih    (idim,kdim+1),   & ! geopotential at half levels
     ppf     (idim,kdim  ),   & ! pressure at full levels 
     pph     (idim,kdim+1),   & ! pressure at half levels 
     pt_g    (idim       ),   & ! surface temperature
     pdqvdt  (idim,kdim)        ! moisture tencency

! Output data
! -----------
  REAL    (KIND=wp   ),     INTENT (OUT) ::  &
     pdt_con  (idim, kdim), & ! convective tendency of temperature
     pdtdiab_con (idim, kdim),& ! conv T-tend. due to LH exchanges
     pdqv_con (idim, kdim), & ! convective tendency of specific humidity
     pdtke_con(idim, kdim), & ! convective boyant TKE prpoduction on half levels
     pco_clw  (idim, kdim), & ! convective cloud liquid water
     pcumflx  (idim),       & ! cu_base massflux
     pcuconv  (idim)          ! moisture convergence used for cu_base massflux
     ! These output-arrays have to be set to zero on call of this routine !

  INTEGER (KIND=iintegers), INTENT (OUT) ::  &
     mbas_con (idim)      , & ! cloud base level index
     mtop_con (idim)          ! cloud top  level index

  LOGICAL                 , INTENT (OUT) ::  &
     locum(idim)          ! indicatior for convection at gridpoints

! Local scalars and automatic arrays (also to be used in lower level routines):
! ----------------------------------
  INTEGER (KIND=iintegers) ::  &
    i, k,                  & ! loop indices over spatial dimensions
    mkb,                   & !
    mlab   (idim,kdim),    & !
    mtype  (idim)            !

  REAL    (KIND=wp   )     ::  &
    ztu    (idim,kdim) ,      & !
    zqu    (idim,kdim) ,      & !
    zlu    (idim,kdim) ,      & !
    zlude  (idim,kdim) ,      & !
    zmfu   (idim,kdim) ,      & !
    zmfus  (idim,kdim) ,      & !
    zmfuq  (idim,kdim) ,      & !
    zmful  (idim,kdim) ,      & !


    zqsen  (idim,kdim) ,      & ! saturation specific humidiy at full levels
    ztenh  (idim,kdim) ,      & !
    zqenh  (idim,kdim) ,      & !
    zqsenh (idim,kdim) ,      & !
    zrhgdz (idim,kdim)          ! rho*g*dz

  REAL    (KIND=wp   )     ::  &
    z1dp      (idim)   ,      & !
    zrho      (idim)   ,      & !
    zqold     (idim)   ,      & !
    zentr     (idim)   ,      & !
    zmfub     (idim)   ,      & !
    zdqpbl    (idim)   ,      & !
    zdmfen    (idim)   ,      & !
    zdmfde    (idim)            !

  REAL    (KIND=wp   )     ::  &
    zbuo, zc3, zc4, zzs, zqumqe, zdqmin, zpbmpt, zzp, zl,         & !
    zdmfmax, zdprho, zzentr, zqeen, zseen, zscde, zqude,          & !
    zmfusk, zmfuqk, zmfulk, zzdmf

  LOGICAL                  ::  &
    loflag (idim), llo1          !              
 
!------------ End of header ---------------------------------------------------
 
!------------------------------------------------------------------------------
! Begin Subroutine cu_shallow             
!------------------------------------------------------------------------------
 
  ! the highest level of these variables is not set below!
  pdt_con  (:,1)   = 0.0_wp
  pdqv_con (:,1)   = 0.0_wp
  pdtke_con(:,1)   = 0.0_wp
  pdtdiab_con(:,1) = 0.0_wp
 
!******************************************************************************
!*                                                                            *
!* Section 1: Initialisation of some variables and computation of additional  *
!*            half level properties ( former subroutine cu_ini )              *
!*                                                                            *
!******************************************************************************
 
  ! Initialization of arrays needed in the parameterisation of
  ! convection; in particular interpolation of large-scale
  ! (environmental) fields to model half levels and determination
  ! of the level of maximum (upward) vertical velocity
 
  ! 1. Specify saturation specific humidity, large scale parameter at half 
  !    levels, adjust temperature if statically unstable, find level of maximum
  !    vertical velocity.
  !    Half level geopotential set to central value between full layers !
  !    First guess for half level temperature from adiabatic interpolation
  !    of adjacent full layer values (maximum of two values is taken)
  !    First guess for half level saturation humidity is saturation
  !    humidity at upper full level

  DO k = 1, kdim
    DO i = isc, iec
      IF( pt(i,k) - ctmelt > 0.0_wp) THEN     
          zc3    = b2w
          zc4    = b4w
        ELSE                            
          zc3    = b2i
          zc4    = b4i
        END IF
        zqsen(i,k) = cc2*EXP( zc3*(pt(i,k)-B3)/(pt(i,k)-zc4) )/ppf(i,k)
        zqsen(i,k) = zqsen(i,k) / ( 1.0_wp - rvd_m_o*zqsen(i,k) )
    ENDDO
  ENDDO

  DO k = 2, kdim
    DO i = isc, iec
      ztenh(i,k) = ( MAX (cp_d*pt(i,k-1)+pfif(i,k-1),cp_d*pt(i,k)+pfif(i,k) ) &
                    - pfih(i,k) )/cp_d  
      zqsenh(i,k)= zqsen(i,k-1)
      z1dp  (i)  = 1.0_wp/pph(i,k)
      loflag(i)  = .TRUE.
    END DO
 
    ! adjust temperature and humidity at half level to value for moist
    ! adiabatic interpolation from neighbour levels (i.e. consider 
    ! condensation/evaporation process)

    CALL cu_cond ( ztenh(:,k), zqsenh(:,k), z1dp,              &
                   loflag , .TRUE.  , .TRUE.,                  &
                   idim   , isc     , iec            )
 
    ! interpolation of specific humidity to half levels, following
    ! a moist adiabate from the upper full level and avoiding supersaturation
    DO i = isc, iec
      zqenh(i,k) = MIN( pqv(i,k-1), zqsen(i,k-1) ) + zqsenh(i,k) - zqsen(i,k-1) 
      zqenh(i,k) = MAX( zqenh(i,k), 0.0_wp )
    END DO
   
  END DO    ! vertical loop
 
  ! lowest layer
  DO i = isc, iec
!   avoid 'convective drizzle' by skipping the next two statements
!   ztenh (i,kdim) = (cp_d*pt(i,kdim) + pfif(i,kdim) - pfih(i,kdim))/cp_d  
!   zqenh (i,kdim) = pqv(i,kdim)
    ztenh (i,1   ) = pt  (i,1)
    zqenh (i,1   ) = pqv (i,1)
  END DO

  ! avoid unstable structure in half level temperature profile 
  DO k = kdim-1, 2, -1
    DO i = isc, iec
       zzs = MAX( cp_d*ztenh(i,k)+pfih(i,k), cp_d*ztenh(i,k+1)+pfih(i,k+1) )
       ztenh(i,k) = ( zzs - pfih(i,k) )/cp_d  ! ztenh may become incompatible
    END DO                                    ! with zqenh
  END DO
 
  ! Initialize updraft and downdraft values 
  DO k = 1, kdim
    DO i = isc, iec
      ztu    (i,k) = ztenh (i,k)
      zqu    (i,k) = zqenh (i,k)
      zlu    (i,k) = 0.0_wp
      zmfu   (i,k) = 0.0_wp
      zmfus  (i,k) = 0.0_wp
      zmfuq  (i,k) = 0.0_wp
      zlude  (i,k) = 0.0_wp
      mlab   (i,k) = 0
      zrhgdz (i,k) = ppf(i,k)*( pfih(i,k) - pfih(i,k+1) )/(r_d*pt(i,k))
    END DO
  END DO

!******************************************************************************
!*                                                                            *
!* Section 2: Cloud base calculations                                         *
!*                                                                            *
!******************************************************************************

  ! a) Determination of cloud-base values
  !     Method
  !       humidity, pressure and geopotential at half levels this
  !       routine computes corresponding cloud base values and the
  !       following flag indices:
  !             mlab=1 as indicator of sub-cloud levels
  !             mlab=3 as indicator of unstable in-cloud levels
  !     - to define a cloud base, surface air is lifted dry-adiaba-
  !       tically up to cloud base (non-entraining plume, i.e. for
  !       constant massflux)
  !     - temperature and humidity are adjusted inside of the cloud
  !       to take into account condensation effects

  ! Initialize lifting level variables
  DO i = isc, iec
    mlab (i,kdim ) = 1        ! lowest layer is below cloud base
    mbas_con(i)    = kdim-1   ! cloud base index
    mtop_con(i)    = kdim-1   ! cloud base index
    locum(i)       =.false.
  END DO

  DO k = kdim-1, 2, -1     ! Vertical loop over layers

    DO i = isc, iec
      z1dp (i) = 1.0_wp/pph(i,k)
      IF ( mlab(i,k+1).EQ.1 ) THEN
        loflag(i) = .TRUE.   ! ascent continues
      ELSE
        loflag(i) = .FALSE.  ! ascent complete
      ENDIF
    END DO

    DO i = isc, iec
      IF( loflag(i) ) THEN        ! cloud base not found yet
        zqu(i,k) = zqu(i,k+1) ! retain parcel's humdity
                      ! parcel temperature after ascent to next layer
        ztu(i,k) = ( cp_d*ztu(i,k+1) + pfih(i,k+1) - pfih(i,k) )/cp_d  

        ! difference between parcel (virtual) temperature and environmental
        ! (virtual) temperature determines buoancy of parcel (0.5K added)
        zbuo =   ztu  (i,k) * ( 1.0_wp + rvd_m_o*zqu  (i,k) )  &
               - ztenh(i,k) * ( 1.0_wp + rvd_m_o*zqenh(i,k) ) + 0.5_wp
        IF(zbuo > 0.0_wp) mlab(i,k) = 1 ! sub-cloud indicator for positve buoyancy
        zqold(i) = zqu(i,k) ! store parcel humidity in local variable
      END IF               
    END DO

    !Check for condensation and adjust parcel temperature and humidity
    !if condensation/sublimation occurs

    CALL cu_cond ( ztu(:,k), zqu(:,k), z1dp  ,                 &
                   loflag , .TRUE.  , .FALSE.,                 &
                   idim   , isc     , iec             )

!     If ascent calculations are still active and parcel humidity
!     changed due to condensation:
!      
    DO i = isc, iec
      IF( loflag(i) .AND. zqu(i,k).NE.zqold(i) ) THEN
        mlab(i,k) = 2 ! Index for lifting condensation level
        zlu(i,k) = zlu(i,k) + zqold(i) - zqu(i,k)
        zbuo =   ztu  (i,k)*( 1.0_wp + rvd_m_o*zqu  (i,k) )  &
               - ztenh(i,k)*( 1.0_wp + rvd_m_o*zqenh(i,k) ) + 0.5_wp
        IF(zbuo > 0.0_wp) THEN ! Test for buoyancy
          mbas_con(i) = k        ! define cloud base index
          locum(i)    = .TRUE.   ! indicate existence of unstable cloud base
          mlab(i,k) = 3          ! index for unstabale LCL
        END IF
      END IF
    END DO
 
  END DO   ! Vertical loop

   
  ! b) total moisture convergence and decision on type of convection
 
  DO i = isc, iec
    zdqpbl(i) = 0.0_wp
  ENDDO

  DO k = 2, kdim
    DO i = isc, iec
      IF ( k >= mbas_con(i) )  THEN
        zdqpbl(i) = zdqpbl(i) + pdqvdt(i,k)*zrhgdz(i,k)
      ENDIF
    ENDDO
  ENDDO

  DO i = isc, iec
    mtype(i) = 2     ! only shallow convection allowed here
  ENDDO
  

!******************************************************************************
!*                                                                            *
!* Section 3: Moisture supply in boundary layer and preliminary cloud base    *
!*            mass flux (excluding downdraft effects)                         *
!*                                                                            *
!******************************************************************************
 

   ! Moisture convergence closure
 
    DO i = isc, iec
      mkb = mbas_con(i)
      zqumqe = zqu(i,mkb) + zlu(i,mkb) - zqenh(i,mkb)
      zdqmin = MAX( 0.01_wp*zqenh(i,mkb), 1.E-10_wp )
!     avoid 'convective drizzle' by using a minimum moisture convergence
!     llo1   = zdqpbl(i) > 0.0_wp   & ! positive moisture convergence
      llo1   = zdqpbl(i) > (1.E-3_wp*g*zqumqe) & ! minimum positive moist. conv.
              .AND. zqumqe > zdqmin & ! parcel humidity exceeds env.hum.
              .AND. locum(i)          ! convective grid point
      IF ( llo1 ) THEN
        zmfub(i) = zdqpbl(i) / (g*MAX(zqumqe, zdqmin))
      ELSE
        zmfub(i) = 0.0_wp
        locum(i) = .FALSE.            ! GP switched off
      ENDIF
      zmfub(i) = MIN( zmfub(i), cmfcmax)   !  massflux upper limit
      pcumflx(i)= zmfub(i)
      IF (zdqpbl(i).gt.0.0001_wp) THEN
        pcuconv(i) = zdqpbl(i)
      ELSE
        pcuconv(i) = 0.0_wp
      END IF
      zentr(i) = entr_sc
    ENDDO


!******************************************************************************
!*                                                                            *
!* Section 5: Cloud ascent for entraining plume                               *
!*                                                                            *
!******************************************************************************

! a) cloud ascent calcualation for an entraining plume to obtain the
!    vertical in-cloud profiles of T, q and l
!    Surface air is lifted dry-adiabatically to cloud base from there, a moist 
!    ascent is calculated for an entraining/detraining plume; 
!    cases where shallow convection does not occur, a check for possible 
!    existence of mid-level (shallow) convection is made

! Security parameter
  zdmfmax = cmfcmax /( kdim/4 )
 
! Set defaults

  DO i = isc, iec
    IF( .NOT.locum(i) ) mtype(i) = 0   ! type of convection
  ENDDO

  DO k = 1, kdim
    DO i = isc, iec
  !   zlu   (i,k) = 0.0_wp
      zmfu  (i,k) = 0.0_wp
      zmfus (i,k) = 0.0_wp
      zmfuq (i,k) = 0.0_wp
  !   zmful (i,k) = 0.0_wp
      zlude (i,k) = 0.0_wp
    ENDDO
  ENDDO
 
  DO i = isc, iec
    IF (.NOT.locum(i) ) THEN
      DO k = 1, kdim 
        mlab(i,k) = 0
      ENDDO
      mbas_con(i) = kdim -1
      mtop_con(i) = kdim -1
      zmfub(i)     = 0.0_wp
    ENDIF
  ENDDO
 
! Initialize values at lifting level
 
  DO i = isc, iec
    IF ( locum(i) ) THEN
      mkb = mbas_con(i)
      zmfu (i,mkb) = zmfub(i)
      zmfus(i,mkb) = zmfub(i)*( cp_d*ztu(i,mkb) + pfih(i,mkb) )
      zmfuq(i,mkb) = zmfub(i)*zqu(i,mkb)
      zmful(i,mkb) = zmfub(i)*zlu(i,mkb)
      mtop_con(i)  = mbas_con(i)
    ENDIF
  ENDDO

 
! Perform ascent:  a) dry adiabatic lifting and 
!                  b) allowing for condensation
!                  c) check buoancy and set flags
!                     (mlab=1 :sub-cloud layer, mlab=2 :cloud layer)
  
  DO k = kdim-1, 2, -1
 
    DO i = isc, iec
      IF( mlab(i,k+1) == 0 ) mlab(i,k) = 0
      IF( mlab(i,k+1) == 3 ) THEN
        loflag(i)= .TRUE.   ! active (unstable) grid point  in layer below
      ELSE
        loflag(i)= .FALSE.  ! inactive grid point                          
      ENDIF
    ENDDO

    ! calculation of entrainement/detrainement rates
    DO i = isc, iec
      zdmfen(i) = 0.0_wp
      zdmfde(i) = 0.0_wp
      zrho  (i) = pph(i,k+1)/(r_d*ztenh(i,k+1))
      z1dp  (i) = 1.0_wp/pph(i,k)
    ENDDO
 
    DO i = isc, iec
      IF( loflag(i) ) THEN
        zdprho = ( pfih(i,k) - pfih(i,k+1) ) / g
        zzentr = zentr(i)*zmfu(i,k+1)*zdprho
        IF( k < mbas_con(i) ) THEN
          zdmfde(i) = zzentr
          zdmfen(i) = zzentr
        ENDIF
      ENDIF
    ENDDO
 
    ! adiabatic ascent for entraining/detraining plume

    DO i = isc, iec
      IF(loflag(i)) THEN
        zdmfen(i) = MIN( zdmfen(i), zdmfmax )
        zdmfde(i) = MIN( zdmfde(i), 0.75_wp*zmfu(i,k+1) )
        zmfu(i,k)= zmfu(i,k+1) + zdmfen(i) - zdmfde(i)
        zqeen = zqenh(i,k+1)*zdmfen(i)
        zseen = ( cp_d*ztenh (i,k+1) + pfih(i,k+1) )*zdmfen(i)
        zscde = ( cp_d*ztu(i,k+1) + pfih(i,k+1) )*zdmfde(i)
        zqude = zqu(i,k+1)*zdmfde(i)
        zlude(i,k) = zlu(i,k+1)*zdmfde(i)
        zmfusk = zmfus(i,k+1) + zseen - zscde
        zmfuqk = zmfuq(i,k+1) + zqeen - zqude
        zmfulk = zmful(i,k+1) - zlude(i,k)
        zlu (i,k) =  zmfulk*( 1.0_wp/MAX(cmfcmin,zmfu(i,k)) )
        zqu (i,k) =  zmfuqk*( 1.0_wp/MAX(cmfcmin,zmfu(i,k)) )
        ztu (i,k) = (zmfusk*( 1.0_wp/MAX(cmfcmin,zmfu(i,k)) ) -      &
                                  pfih(i,k))/cp_d  
        ztu(i,k)  = MAX( 100.0_wp, ztu(i,k) )
        ztu(i,k)  = MIN( 400.0_wp, ztu(i,k) )
        zqold(i)  = zqu(i,k)
      ENDIF
    ENDDO
 
    ! corrections for moist ascent by adjusting T,q and l
    ! calulation of condensation and corresponding adjustment of T and q

    CALL cu_cond ( ztu(:,k), zqu(:,k), z1dp  ,              &
                   loflag , .TRUE.  , .FALSE.,                 &
                   idim   , isc     , iec            )

    DO i = isc, iec
      IF( loflag(i) .AND.  zqu(i,k).NE.zqold(i) ) THEN
        zlu (i,k) = zlu(i,k) + zqold(i) - zqu(i,k)
        zbuo      =     ztu(i,k)*(1.0_wp+rvd_m_o*zqu(i,k))  &
                    - ztenh(i,k)*(1.0_wp+rvd_m_o*zqenh(i,k))
        IF(mlab(i,k+1)==1) zbuo = zbuo + 0.5_wp
        IF( zbuo > 0.0_wp .AND. zmfu(i,k) >= 0.1_wp*zmfub(i) ) THEN
          mtop_con(i) = k
          mlab(i,k)   = 3
        ELSE
          mlab(i,k)  = 0
          zmfu(i,k) = 0.0_wp
        ENDIF
      ENDIF
    ENDDO

    DO i = isc, iec
      IF( loflag(i) ) THEN
        zmful(i,k) = zlu(i,k)*zmfu(i,k)
        zmfus(i,k) = ( cp_d*ztu(i,k) + pfih(i,k) )*zmfu(i,k)
        zmfuq(i,k) = zqu(i,k)*zmfu(i,k)
      ENDIF
    ENDDO

  ENDDO         ! vertical loop
 
!     convective fluxes above non-buoancy level
!     
!        (NOTE: CLOUD VARIABLES LIKE T,Q and L ARE NOT
!               AFFECTED BY DETRAINMENT and ARE ALREADY KNOWN
!               FROM PREVIOUS CALCULATIONS ABOVE)
 
  DO i = isc, iec 
    IF( mtop_con(i) == kdim-1 ) locum(i) = .FALSE.  
    mbas_con(i) = MAX( mbas_con(i), mtop_con(i) )
  ENDDO

  DO i = isc, iec
    IF( locum(i) ) THEN
      k          = mtop_con(i) - 1
      zzdmf      = cmfctop
      zdmfde(i)  = ( 1.0_wp - zzdmf )*zmfu(i,k+1)
      zlude(i,k) = zdmfde(i)*zlu(i,k+1)
      zmfu(i,k) = zmfu(i,k+1) - zdmfde(i)
      zmfus(i,k)  = ( cp_d*ztu(i,k) + pfih(i,k) )*zmfu(i,k)
      zmfuq(i,k)  = zqu(i,k)*zmfu(i,k)
      zmful(i,k)  = zlu(i,k)*zmfu(i,k)
      zlude(i,k-1) = zmful(i,k)
    ENDIF
  ENDDO


  ! Additional checks 
  ! Check cloud depth and delete deep convection
 
  DO i=isc,iec
    IF (locum(i)) THEN
      ! cloud thickness (in Pa)
      zpbmpt = pph(i,mbas_con(i)) - pph(i,mtop_con(i))
      IF (zpbmpt > thick_sc) THEN
        locum(i) = .FALSE.
        mbas_con(i) = 0
        mtop_con(i) = 0
      ENDIF
    ENDIF
  ENDDO
 
  ! Store in-cloud water content
  DO k = 1, kdim
    DO i = isc, iec
      IF (locum(i)) THEN
        pco_clw(i,k) = MAX ( 0.0_wp, MIN(100.0_wp, zlu(i,k)) )
      ENDIF
    ENDDO
  ENDDO

!******************************************************************************
!*                                                                            *
!* Section 7: Determine final convective fluxes (former subroutine cu_flux)   *
!*            and conside evaporation of rain in sub-cloud layers             *
!*                                                                            *
!******************************************************************************

  DO i = isc, iec
     IF( .NOT.locum(i))  mtype(i) = 0
  ENDDO

  DO k = 2, kdim
    DO i = isc, iec

      IF( locum(i) .AND. k >= mtop_con(i)-1 ) THEN
         zmfus(i,k) = zmfus(i,k) - zmfu (i,k)*(cp_d*ztenh(i,k)+pfih(i,k))
         zmfuq(i,k) = zmfuq(i,k) - zmfu (i,k)*zqenh(i,k)
      ELSE              ! above cloud and non-convective points
         zmfu  (i,k  ) = 0.0_wp
         zmfus (i,k  ) = 0.0_wp
         zmfuq (i,k  ) = 0.0_wp
         zmful (i,k  ) = 0.0_wp
         zlude (i,k-1) = 0.0_wp
      ENDIF
    ENDDO
  ENDDO 

  DO k = 2, kdim
    DO i = isc, iec
      IF( locum(i) .AND. k > mbas_con(i) ) THEN
        mkb = mbas_con(i)
        zzp = ( pph(i,kdim+1) - pph(i,k) )/( pph(i,kdim+1) - pph(i,mkb) ) 
        IF ( mtype(i).EQ.3 ) zzp=zzp**2
        zmfus(i,k) = zmfus(i,mkb) * zzp
        zmfuq(i,k) = zmfuq(i,mkb) * zzp
        zmful(i,k) = zmful(i,mkb) * zzp
      ENDIF
    ENDDO
  ENDDO      
 
 
!******************************************************************************
!*                                                                            *
!* Section 8: Compute the final tendencies for grid scale variables T and qv  *
!*            (former subroutine cu_dtdq)                                     *
!*                                                                            *
!******************************************************************************

  DO k = 2, kdim - 1  ! above lowest model layer
    DO i = isc, iec
      IF(locum(i)) THEN ! for convective grid points only
        llo1 = (pt(i,k)-ctmelt) > 0.0_wp  !
        IF ( llo1 ) THEN
          zl = lh_v
        ELSE
          zl = lh_s
        ENDIF
        pdt_con(i,k) =    g/zrhgdz(i,k)/cp_d  &
                      * (  zmfus (i,k+1) - zmfus (i,k) &
                         - zl*( zmful(i,k+1)-zmful(i,k) ) )
                               
        pdtdiab_con(i,k) =    g/zrhgdz(i,k)/cp_d  &
                      * ( - zl*( zmful(i,k+1)-zmful(i,k) ) )

        pdqv_con(i,k)=   g/zrhgdz(i,k)                   &
                       *(  zmfuq(i,k+1) - zmfuq(i,k)   &
                         + zmful(i,k+1) - zmful(i,k)   )
      END IF
    END DO
  ENDDO

  k = kdim      ! special case: lowest model layer

  DO i=isc,iec
    IF(locum(i)) THEN      ! for convective grid points only
      llo1 = (pt(i,k)-ctmelt) > 0.0_wp
      IF ( llo1 ) THEN
         zl = lh_v
      ELSE
         zl = lh_s
      ENDIF
      pdt_con(i,k) = - g /zrhgdz(i,k)/cp_d * ( zmfus(i,k) -zl*zmful(i,k) )
      pdtdiab_con(i,k) = - g /zrhgdz(i,k)/cp_d * ( -zl*zmful(i,k) )
      pdqv_con(i,k)= - g /zrhgdz(i,k) * ( zmfuq(i,k)  + zmful(i,k) )
    END IF
  END DO

  DO k = 2, kdim ! for all model half levels (except the lower model boundary)
    DO i = isc, iec
      IF (locum(i)) THEN ! for convective grid points only
        ! RA: Computation of the convective buoyant TKE production:
        pdtke_con(i,k) = MAX( 0.0_wp, g*r_d/pph(i,k) *    &
                         ( (1.0_wp+rvd_m_o*zqenh(i,k))*zmfus(i,k)/cp_d + &
                           rvd_m_o*ztenh(i,k)*zmfuq(i,k) ))

        ! RA
      END IF
    END DO
  END DO

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------
END SUBROUTINE cu_shallow

!==============================================================================
!+ Module procedure in "Convection"
!------------------------------------------------------------------------------

SUBROUTINE cu_cond (                                               &
           pt     , pqv     , p1op   , plflag , pcflag , peflag,   &
           idim   , isc     , iec                                  )

!-----------------------------------------------------------------------------!
! Description:
!
!   The module procedure cu_cond does a saturation adjustment for
!   temperature and specific humidity.
!
!   Method:    Thermodynamic adjustment by instantaneous condensation
!              at constant pressure using a double iteration method.
!              Release of latent heat of condendation and of deposition
!              is considered depending on temperature.
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------
! Subroutine arguments:
! --------------------
! Input data
! ----------
  INTEGER (KIND=iintegers), INTENT (IN) ::  &
     idim ,        & ! array dimension in zonal direction
     isc  ,        & ! start index for first  array computation
     iec             ! end   index for first  array computation

  REAL     (KIND=wp   ),     INTENT (IN) ::  &
     p1op  (idim)    ! reciprocal of pressure, 1.0/p

  LOGICAL                 ,  INTENT (IN) ::  &
     plflag (idim),& ! switch for points where adjustment shall be made
     pcflag,       & ! condensation only (.TRUE)
     peflag          ! evaporation only  (.TRUE)

! Input/Output data
! -----------
  REAL     (KIND=wp   ),     INTENT (INOUT) ::  &
     pt    (idim), & ! temperature on input, adjusted on output
     pqv   (idim)    ! specific humidity on input, adjusted on ouput

! Local scalars and automatic arrays:
! ----------------------------------
  INTEGER (KIND=iintegers) ::  &
    i                ! loop indix

  REAL    (KIND=wp   )     ::  &
    zcond(idim)      ! condensation amount

  REAL    (KIND=wp   )     ::  &
    zc3, zc4, zc5, zhldcp, zqsat, zcor, zcond1, zfacc, zface  ! local storage

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine cu_cond
!------------------------------------------------------------------------------
    zfacc = 0.0_wp
    zface = 0.0_wp
    IF (pcflag) zfacc = 1.0_wp
    IF (peflag) zface = 1.0_wp

    DO i = isc, iec
      zcond(i) = 0.0_wp     ! Initialize condensation variable
      IF(plflag(i)) THEN        ! only, if ascent still continues
        IF( pt(i) - ctmelt > 0.0_wp) THEN     ! condensation
          zc3    = b2w
          zc4    = b4w
          zc5    = c5hlccp
          zhldcp = chlcdcp
        ELSE                                       ! deposition
          zc3    = b2i
          zc4    = b4i
          zc5    = c5hlscp
          zhldcp = chlsdcp
        END IF
        zqsat = cc2*EXP( zc3*(pt(i)-b3)/(pt(i)-zc4) )*p1op(i)
        zqsat = MIN( 0.5_wp, zqsat )
        zcor  = 1.0_wp / ( 1.0_wp - rvd_m_o*zqsat )
        zqsat = zqsat*zcor
        zcond1   = (pqv(i)-zqsat)/(1.0_wp+zc5*zqsat*zcor/(pt(i)-zc4)**2)
        zcond(i) = zfacc*MAX( zcond1, 0.0_wp )  + &
                   zface*MIN( zcond1, 0.0_wp )
        pt(i)    = pt(i) + zhldcp*zcond(i)
        pqv(i)   = pqv(i) - zcond(i)
      END IF
    END DO
    !Second iteration
    DO i = isc, iec
      IF( plflag(i) .AND. zcond(i).NE.0.0_wp) THEN  !saturation adjustment
        IF( pt(i) - ctmelt > 0.0_wp ) THEN
          zc3    = b2w
          zc4    = b4w
          zc5    = c5hlccp
          zhldcp = chlcdcp
        ELSE
          zc3    = b2i
          zc4    = b4i
          zc5    = c5hlscp
          zhldcp = chlsdcp
        END IF
        zqsat = cc2*EXP( zc3*(pt(i)-b3)/(pt(i)-zc4) )*p1op(i)
        zqsat = MIN( 0.5_wp, zqsat )
        zcor  = 1.0_wp / ( 1.0_wp - rvd_m_o*zqsat )
        zqsat = zqsat * zcor
        zcond1=(pqv(i)-zqsat)/(1.0_wp+zc5*zqsat*zcor/(pt(i)-zc4)**2)
        pt(i)   = pt(i) + zhldcp*zcond1
        pqv(i)  = pqv(i) - zcond1
      END IF
    END DO

!-----------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------
END SUBROUTINE cu_cond

!==============================================================================

!------------------------------------------------------------------------------
! End of module src_conv_shallow
!------------------------------------------------------------------------------

END MODULE src_conv_shallow
