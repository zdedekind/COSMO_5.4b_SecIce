!+ Source module for the interface between organize_physics and turbulence
!----------------------------------------------------------------------------------

MODULE turbulence_interface

!-------------------------------------------------------------------------------
!
! Description:
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  Ulrich.Schaettler@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V4_20        2011/08/31 Ulrich Schaettler
!  Initial release
! V4_23        2012/05/10 Oliver Fuhrer
!  Eliminated qvt_diff from interface to turbulence_diff
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer
!  Replaced qx-variables by using them from the tracer module
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  MESSy interface introduced
! V4_28        2013/07/12 Oliver Fuhrer
!  Removed memory leak: if h_ice_p is allocated, it has to be deallocated again
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler, Ulrich Blahak
!  Unification of MESSy interfaces and COSMO Tracer structure
!  For artifdata: specify tch, t_g and qv_s in a way that predefined surface
!    sensible and latent heat fluxes will result (UB)
! V5_1         2014-11-28 Ulrich Blahak, Matthias Raschendorfer, Oliver Fuhrer
!  Avoiding the CALL of horizontal_diffcoeffs in case of itype_turb=3.
!   These are calculated now in turbdiff() depending on itype_sher. They are 0 for
!   itype_sher<2, isotropic for itype_sher=2 and have an additional
!   horizontal mode for itype_sher=3.
!  Introduced tket_adv for advection of tke in case of itype_turb=3.
!  Replaced ireals by wp (working precision) (OF)
! V5_3         2015-10-09 Ulrich Schaettler
!  Added subroutine turbulence_organize, which organizes the calling of the
!  blocked version of the turbulence
! V5_4a        2016-05-10 Matthias Raschendorfer
!  Use several variables from turb_data, which have been in data_runcontrol before
!  Removed interface for new scheme, which is in a separate module now
! V5_4b        2016-07-12 Ulrich Blahak
!  Added switch lsensiflux_fix
!
! Code Description:
! Language: Fortran 90 
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Modules used:

USE data_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!-------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------

    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction
    ke1,          & ! ke + 1
    kcm,          & ! index of the lowest model layer higher than the canopy

! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the 
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from 
!    the other ones because of the use of the staggered Arakawa-B-grid.
!    
!   zonal direction
    istart,       & ! start index for the forecast of w, t, qd, qw and pp
    iend,         & ! end index for the forecast of w, t, qd, qw and pp
    istartu,      & ! start index for the forecast of u
    iendu,        & ! end index for the forecast of u
    istartv,      & ! start index for the forecast of v
    iendv,        & ! end index for the forecast of v
    istartpar,    & ! start index for computations in the parallel program
    iendpar,      & ! end index for computations in the parallel program

!   meridional direction
    jstart,       & ! start index for the forecast of w, t, qd, qw and pp
    jend,         & ! end index for the forecast of w, t, qd, qw and pp
    jstartu,      & ! start index for the forecast of u
    jendu,        & ! start index for the forecast of u
    jstartv,      & ! start index for the forecast of v
    jendv,        & ! end index for the forecast of v
    jstartpar,    & ! start index for computations in the parallel program
    jendpar,      & ! end index for computations in the parallel program

! 4. constants for the horizontal rotated grid and related variables
! ------------------------------------------------------------------

!   dlam,         & ! grid point distance in zonal direction (in degrees)
!   dphi,         & ! grid point distance in meridional direction (in degrees)
    eddlon,       & ! 1 / dlon
    eddlat,       & ! 1 / dlat
    edadlat,      & ! 1 / (radius of the earth * dlat)

! 5. variables for the time discretization and related variables
! --------------------------------------------------------------

    dt,           & ! long time-step
    dt2,          & ! 2 * dt
    dtdeh,        & ! dt / 3600 seconds
    lalloc_h_ice, & !

! 8. Organizational variables to handle the COSMO humidity tracers
! ----------------------------------------------------------------
    idt_qv, idt_qc

! end of data_modelconfig

!-------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 2. physical constants and related variables
! -------------------------------------------

    r_earth         ! mean radius of the earth

! end of data_constants

!-------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------

    hhl ,           & ! height of model half levels                   ( m )
    p0  ,           & ! reference pressure at main levels             (pa )
    dp0 ,           & ! pressure thickness of layer                   (pa )

! 2. external parameter fields                                        (unit)
! ----------------------------
    gz0,            & ! roughness length * g of the vertically not
                      ! resolved canopy                               (m2/s2)
    fr_land,        & ! land portion of a grid point area             ( 1 )
    depth_lk ,      & ! lake depth                                    ( m )
    sai,            & ! surface area index                            ( 1 )
    lai,            & ! leaf area index                               ( 1 )
    tai,            & ! transpiration area index                      ( 1 )
    eai,            & ! (evaporative) earth area index                ( 1 )
    plcov,          & ! fraction of plant cover                       ( 1 )
    h_can,          & ! hight of the vertically resolved canopy       ( m )
    d_pat,          & ! horizontal pattern length scale               ( m )
    c_big,          & ! effective drag coefficient of canopy elements 
                      ! larger than or equal to the tubulent length 
                      ! scale                                         (1/m)
    c_sml,          & ! effective drag coefficient of canopy elements
                      ! smaller than the tubulent length scale        (1/m) 
    r_air,          & ! log of air containing fraction of a gridbox inside 
                      ! the canopy                                    ( 1 ) 

! 3. prognostic variables                                             (unit)
! -----------------------

    u          ,    & ! zonal wind speed                              ( m/s )
    v          ,    & ! meridional wind speed                         ( m/s )
    w          ,    & ! vertical wind speed (defined on half levels)  ( m/s )
    t          ,    & ! temperature                                   (  k  )
    pp         ,    & ! deviation from the reference pressure         ( pa  )
    tke               ! SQRT(2*TKE); TKE='turbul. kin. energy'        ( m/s )
                      ! (defined on half levels)

USE data_fields     , ONLY :   &

! 4. tendency fields for the prognostic variables                     (unit )
! -----------------------------------------------
!    timely deviation of diffused variables

    ut_sso       ,  & ! u-tendency due to the SSO-Scheme (mass-pos.)  ( m/s2)
    vt_sso       ,  & ! v-tendency due to the SSO-Scheme (mass-pos.)  ( m/s2)
    utens        ,  & ! u-tendency without sound-wave terms           ( m/s2)
    vtens        ,  & ! v-tendency without sound-wave terms           ( m/s2)
    ttens        ,  & ! t-tendency without sound-wave terms           ( k/s )

! defined on half-levels:

    tketens      ,  & ! non-advective tendency of SQRT(2*TKE)         (m/s2)
    tket_adv     ,  & ! pure advective tendency of SQRT(2*TKE)        (m/s2)
    edr          ,  & ! eddy dissipation rate of TKE (EDR)            (m2/s3)
    tket_sso     ,  & ! TKE-tendency due to SSO wake production       (m2/s3)
    tket_hshr    ,  & ! TKE-tendency due to (sep.) horiz. shear       (m2/s3)
    tket_conv    ,  & ! TKE-tendency due to convective buoyancy       (m2/s3)

! 5. fields for surface values and soil/canopy model variables        (unit )
! ------------------------------------------------------------

    ps        ,     & ! surface pressure                              ( pa  )
    h_ice     ,     & ! ice thickness                                 (  m  )
    t_g       ,     & ! weighted surface temperature                  (  k  )
    qv_s      ,     & ! specific water vapor content on the surface   (kg/kg)

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------

    rho        ,    & ! total density of air                          (kg/m3)
    tkvm       ,    & ! turbulent diffusion coefficients for momentum (m/s2 )
    tkvh       ,    & ! turbulent diffusion coefficients for heat     (m/s2 )
                      ! and moisture
                      ! horizontal turbulent diffusion coefficients
    tkhm       ,    & ! ... for momentum            (l3dturb)         (m2/s)
    tkhh       ,    & ! ... for heat and moisture   (l3dturb)         (m2/s)
    rcld       ,    & ! standard deviation of the saturation deficit        
                      ! (as input and output)                             
                      ! fractional cloud cover (in turbdiff)            --
    tcm        ,    & ! turbulent transfer coefficients for momentum    --
    tch        ,    & ! turbulent transfer coefficients for heat        --

    tfm        ,    & ! factor of laminar transfer of momentum          --
    tfh        ,    & ! factor of laminar transfer of scalars           --
    tfv         ,   & ! laminar reduction factor for evaporation        --

! 7. fields for model output and diagnostics                          (unit )
! ------------------------------------------

    t_2m       ,    & ! temperature in 2m                             (  K  )
    qv_2m      ,    & ! specific water vapor content in 2m            (kg/kg)
    td_2m      ,    & ! dew-point in 2m                               (  K  )
    rh_2m      ,    & ! relative humidity in 2m                       (  %  )
    u_10m      ,    & ! zonal wind in 10m                             ( m/s )
    v_10m      ,    & ! meridional wind in 10m                        ( m/s )

! 8. external parameter fields
! ------------------------------------------

    acrlat            ! 1 / ( crlat * radius of the earth )           ( 1/m )

! end of data_fields

!-------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. Time indices                  
! ---------------                  
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
    nold,         & ! corresponds to ntstep - 1
    nnow,         & ! corresponds to ntstep
    nnew,         & ! corresponds to ntstep + 1
                    ! indices for permutation of the tke time levels
    ntke,         & ! time step index of current TKE values
                    ! (corresponds to nold for "leap-frog" or nnow for 2 time levels)
    ninctura,     & ! time step increment for running the vertical diffusion
    l3dturb,      & ! 3D-turbulence: CALL explicit_horizontal_diffusion (RK)
    lartif_data,  & ! forecast with self-defined artificial data
    l2tls,        & ! time integration by two timelevel RK-scheme (.TRUE.)
                    ! or by default three-time level KW-scheme (.FALSE.)
    lperi_x, lperi_y, l2dim, nbl_exchg, &

! 3. controlling the physics
! --------------------------
    itype_tran,   & ! type of surface-atmosphere transfer
    itype_turb,   & ! type of turbulent diffusion parametrization

    lprintdeb_all, idbg_level, ldebug_tur

USE turb_data, ONLY : &

    itype_sher,   & ! type of shear production for TKE
    ltkeshs,      & ! consider separ. horiz. shear production of TKE
    loutshs         ! consider separ. horiz. shear production of TKE for output

! end of data_runcontrol 

!-------------------------------------------------------------------------------

#ifndef SCLM
USE src_artifdata,      ONLY:  tkvhfix, tkhhfix, tkvmfix, tkhmfix, &
                               lnosurffluxes_m, lnosurffluxes_h, &
                               set_idealized_surffluxes, lsensiflux_fix
#endif

!-------------------------------------------------------------------------------

USE data_parallel,      ONLY :  &
    my_cart_id         ! rank of this subdomain in the cartesian communicator

!-------------------------------------------------------------------------------

USE turbulence_utilities, ONLY :  &
    l_hori

!-------------------------------------------------------------------------------

USE turbulence_tran, ONLY :  &
    turbtran                ! transfer scheme of Raschendorfer turbulence

USE turbulence_diff, ONLY :  &
    turbdiff                ! turbulent diffusion scheme of Raschendorfer turbulence

USE src_turbulence, ONLY: &
    parturs, partura, horizontal_diffcoeffs, prankolmo_rk

!-------------------------------------------------------------------------------

USE src_tracer,       ONLY : trcr_get, trcr_errorstr, trcr_meta_get, &
                             trcr_get_ntrcr !number of all tracers

!-------------------------------------------------------------------------------

USE environment,      ONLY :   model_abort

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

CONTAINS

!===============================================================================

!+ organize_turbulence calls the different schemes depending on options set
!  it is working on the ijk (non-blocked) version of the physics

SUBROUTINE organize_turbulence (ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!  This routine organizes the calls to the different turbulence schemes and
!  acts as interface between the COSMO-Model (routine organize_physics) and
!  the different turbulence packages. It uses data from the corresponding
!  COSMO-Modules and passes these on the parameterizations.
!
! Method:
!
!------------------------------------------------------------------------------

! Formal Parameters:

INTEGER (KIND=iintegers), INTENT(OUT)           ::                      &
  ierror       ! error status

CHARACTER (LEN= *),       INTENT(OUT)           ::                      &
  yerrmsg      ! error message

!------------------------------------------------------------------------------

! Local Variables:
INTEGER (KIND=iintegers)    ::   &
  iini,         & ! type of initialization 
                  !      0: no initialization
                  !      1: separate before the time loop (only ICON)
                  !      2: within the first time step
  itim,         & ! any time step index
  nprv,         & ! previous time step index for tke
  ntur,         & ! current new time step index of tke
  ntim,         & ! number of TKE time levels
  nvor,         & ! another index to have something to copy
  nvar,         & ! indicating time level for physics
  tkestep(3),   & ! max. 3 timelevels
  izerror,      & !
  izdebug         ! for debug printout

LOGICAL                     ::   &
  lstfnct         ! calculation of a new stability function in sub. turbdiff

CHARACTER (LEN=255):: yzerrmsg
CHARACTER (LEN=25) :: yzroutine

REAL (KIND=wp)     :: &
  zdt             ! time step to be used (depends on Leapfrog / Runge-Kutta)

REAL (KIND=wp),     POINTER :: &
  h_ice_p(:,:)
 
REAL (KIND=wp)              :: &
  pres(ie,je,ke)  ! full pressure                                 ( pa  )

! Tracer pointers:
REAL (KIND=wp),     POINTER :: &
  qv  (:,:,:)    => NULL(),     &    ! QV at nvar
  qv_tens (:,:,:)=> NULL(),     &    ! QV tendency
  qc  (:,:,:)    => NULL(),     &    ! QC at nvar
  qc_tens (:,:,:)=> NULL()           ! QC tendency

!-End of Header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  ierror    = 0_iintegers
  yerrmsg   = '     '
  yzroutine = 'organize_turbulence'

  IF (ldebug_tur) THEN
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

  ! check, whether atmospheric turbulence / stability function has to be 
  ! computed
  IF ( (ntstep < 10) .OR. ( MOD(ntstep+1,ninctura) == 0 ) ) THEN
    ! calculate turbulence or at least stability functions
    lstfnct = .TRUE.
  ELSE
    ! don't calculate stability functions or even turbulence
    lstfnct = .FALSE.
  ENDIF

  IF (ntke == 0) THEN 
    ! only in the very fist time step an initialization is performed
    iini=2
  ELSE 
    ! in any other time step no initialization
    iini=0
  END IF

  ! set timelevels, ntim, nprv, ntur, which are passed on to routines
  ! number of TKE time levels

  ! get the correct timelevel for dynamic variables
  IF (l2tls) THEN
    nvar = nnow
    zdt  = dt
  ELSE
    nvar = nold
    zdt  = dt2
  ENDIF

  ! retrieve the required microphysics tracers (at specified timelevel)
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nvar, ptr = qv)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qv, ptr_tens = qv_tens )
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nvar, ptr = qc)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tens = qc_tens)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  ntim=UBOUND(tke, 4) 

  ! for testing purposes only
  ! ntim=1

  IF (ntim == 1) THEN
    nprv = 1
    ntur = 1
  ELSEIF (itype_turb == 3) THEN 
    ! independent 2-timlevel stepping for 'turbdiff'
    DO itim = 1, ntim
      tkestep(itim) = MOD(ntke+itim, ntim)+1
    END DO

    ! time step index of present tke  (corresponds to nold/nnow)
    nprv = tkestep(ntim-1) 
    ! time step index of updated tke  (corresponds to nnew)
    ntur = tkestep(ntim)   
  ELSE 
    ! any other possible turbulence scheme with prognostic TKE 
    ! using at least 'turbtran'
    ! time step index of present tke is that of the other present variables
    nprv = nvar   ! nnow / nnold, depending on Leapfrog / Runge-Kutta
    ! time step index of updated tke is that of the other updated variables
    ntur = nnew
  END IF

  nvor = nprv

  ! compute full pressure
  pres(:,:,:) = p0(:,:,:) + pp(:,:,:,nvar)

  ! compute additional variables for the schemes
  l_hori = r_earth / SQRT(eddlon*eddlat)

  IF (lalloc_h_ice) THEN
    h_ice_p => h_ice(:,:,nnow)
  ELSE
    ALLOCATE (h_ice_p(ie,je))
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Surface transfer scheme
!------------------------------------------------------------------------------

  IF (izdebug > 5) THEN
    PRINT *, '  new TURBULENCE with itype_tran, lstfnct = ', itype_tran, lstfnct
  ENDIF

  IF     (itype_tran == 1) THEN
    ! former Louis scheme
    CALL parturs

  ELSEIF (itype_tran == 2) THEN
    ! Default scheme for surface-to-atmosphere transfer by Raschendorfer based 
    ! on an iterative solution of the turbulent budget equations 
    ! for the surface layer:

!US old scheme    CALL turbtran(nvar,dt)

    CALL turbtran(dt, lstfnct, iini, nprv, ntur, ntim, nvor,                   &
                  ie, je, ke, ke1, kcm,     1_iintegers,                       &
                  istartpar, iendpar, jstartpar, jendpar,                      &
                  hhl, fr_land, depth_lk, sai,                                 &
                  u(:,:,:,nvar), v(:,:,:,nvar), t(:,:,:,nvar), qv(:,:,:),      &
                  qc(:,:,:)    , pres(:,:,:)  , ps(:,:,nvar) , qv_s(:,:,nvar), &
                  t_g(:,:,nvar), h_ice_p,                                      &
                  gz0, tcm, tch, tfm, tfh, tfv,                                &
                  tke, tkvm, tkvh, rcld,                                       &
                  t_2m, qv_2m, td_2m, rh_2m, u_10m, v_10m,                     &
                  edr,                                                         &
                  yerrormsg=yerrmsg, ierrstat=ierror)
  ! optional parameters shfl_s, lhfl_s not used yet

  ENDIF

  ! special option for idealized test cases
#ifndef SCLM
  IF (lartif_data) THEN

    ! Specify tch, t_g and qv_s in a way that predefined surface
    ! sensible and latent heat fluxes result, based on
    ! some NAMELIST parameters:
    IF (lsensiflux_fix) THEN
      CALL set_idealized_surffluxes()
    ENDIF

    ! Set surface diffusion coefficients and therefore surface fluxes 
    ! of momentum, heat and moisture to 0:
    IF (lnosurffluxes_m) THEN
      tcm(:,:) = 0.0_wp
    ENDIF
    IF (lnosurffluxes_h) THEN
      tch(:,:) = 0.0_wp
    ENDIF
  ENDIF
#endif

!------------------------------------------------------------------------------
! Section 3: Atmospheric diffusion scheme
!------------------------------------------------------------------------------

  IF (izdebug > 5) THEN
    PRINT *, '  new TURBULENCE with itype_turb, lstfnct = ', itype_turb, lstfnct
  ENDIF

  SELECT CASE( itype_turb )

  CASE( 1 )

    ! former Mueller scheme with diagnostic TKE
    IF (lstfnct) THEN
      CALL partura
      IF ( l3dturb ) CALL horizontal_diffcoeffs
    ENDIF

  CASE( 3 )

    ! Default scheme for atmospheric diffusion by Raschendorfer based on a 
    ! generalized prognostic TKE equation and turbulent budget equations for 
    ! quasi conservative variables:

    IF (izdebug > 5) THEN
      PRINT *, '      TURBULENCE with iini, nprv, ntur:  ', iini, nprv, ntur
    ENDIF

!US    CALL turbdiff(nvar,zdt,dt,lstfnct)

    CALL turbdiff(zdt, dt, lstfnct, iini, ntur, ntim, nvor,                   &
                  ie, je, ke, ke1, kcm,   1,                                  &
                  istartpar, iendpar, jstartpar, jendpar,                     &
                  istart   , iend   , jstart   , jend   ,                     &
                  istartu  , iendu  , jstartu  , jendu  ,                     &
                  istartv  , iendv  , jstartv  , jendv  ,                     &
                  hhl=hhl, fr_land=fr_land, depth_lk=depth_lk,                &
                  d_pat=d_pat, c_big=c_big, c_sml=c_sml, r_air=r_air, dp0=dp0,&
                  u=u(:,:,:,nvar), v=v(:,:,:,nvar), w=w(:,:,:,nvar),          &
                  t=t(:,:,:,nvar), qv=qv(:,:,:)   , qc=qc(:,:,:)   ,          &
                  prs=pres(:,:,:), ps=ps(:,:,nvar), qv_s=qv_s(:,:,nvar),      &
                  t_g=t_g(:,:,nvar), h_ice=h_ice_p,                           &
                  gz0=gz0, tcm=tcm, tch=tch, tfm=tfm, tfh=tfh, tfv=tfv,       &
                  tke=tke, rcld=rcld,                                         &
                  tkvm=tkvm, tkvh=tkvh, tkhm=tkhm, tkhh=tkhh,                 &
                  u_tens=utens, v_tens=vtens, t_tens=ttens, qv_tens=qv_tens,  &
                  qc_tens=qc_tens, tketens=tketens, tketens_adv=tket_adv,     &
                  ut_sso=ut_sso, vt_sso=vt_sso, tket_conv=tket_conv,          &
                  tket_sso=tket_sso, tket_hshr=tket_hshr, edr=edr,            &
                  yerrormsg=yerrmsg, ierrstat=ierror)

  CASE( 5:8 )

    CALL prankolmo_rk( lstfnct )

#ifndef SCLM
  CASE( 100 )

    ! special option for idealized test cases
    IF (lartif_data) THEN
      ! Set constant diffusion coefficients from namelist parameters:
      WRITE (*,'(a,3(es8.1,","),es8.1,a)')   &
             '!!!!!!! ATTENTION: tkvh, tkhh, tkvm, tkhm = const = ', &
             tkvhfix,tkhhfix,tkvmfix,tkhmfix,' !!!!!!!!!!!!!!'
      tkvh = tkvhfix
      tkvm = tkvmfix
      IF (l3dturb) THEN
        tkhh = tkhhfix
        tkhm = tkhmfix
      END IF
    ELSE
      yerrmsg  = ' ERROR    *** itype_turb = 100 (const. diff. coeff.) only '//&
           'if lartif_data=.TRUE. *** '
      ierror   = 2798
      RETURN
    END IF
#endif

  END SELECT

  ! Update ntke:
  ntke = ntur   ! update index of current TKE time level after 'turbdiff' has passed

  ! Deallocate h_ice_p again, if it has been allocated
  IF (.NOT. lalloc_h_ice) THEN
    DEALLOCATE (h_ice_p)
  ENDIF

  ! Special option for idealized test cases
  ! Set surface diffusion coefficients and therefore surface fluxes 
  ! of momentum, heat and moisture to 0:
#ifndef SCLM
  IF (lartif_data) THEN
    IF (lnosurffluxes_m) THEN
      tcm(:,:) = 0.0_wp
    END IF
    IF (lnosurffluxes_h) THEN
      tch(:,:) = 0.0_wp
    END IF
  END IF
#endif

END SUBROUTINE organize_turbulence

!===============================================================================

END MODULE turbulence_interface
