!+ Some utility routines for the radiation
!------------------------------------------------------------------------------

MODULE radiation_utilities

!------------------------------------------------------------------------------
!
! Description: Routines included
!
!  - compute_sunshine_conditions (
!  - surface_albedo              (
!  - cloudiness_and_humidity     (
!  - co2_scenarios_and_absorbers (
!  - calc_rad_corrections        (
!
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  Ulrich.Schaettler@dwd.de
!
! History:
! Version      Date       Name
! ----------   ---------- ----
! V5_3         2015-10-09 Ulrich Schaettler
!  Initial release
! V5_3a        2015-11-24 Ulrich Schaettler
!  Adaptations in ACC code to run new radiation on the GPUs
! V5_4         2016-03-10 Xavier Lapillonne
!  Bugfix: added loop indices for variable zcosphi in SR compute_sunshine_conditions
!  Adaptations for running radiation on GPUs
! V5_4a        2016-05-10 Ulrich Schaettler, Matthias Raschendorfer
!  Use new module turb_data instead of data_turbulence
!  Use itype_wcld now from turb_data (instead of runcontrol)
! @VERSION@    @DATE@     Ulrich Schaettler
!  Use isoiltyp instead of soiltyp now
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

!USE parallel_utilities, ONLY :  global_values
USE data_parallel,      ONLY :  &
  nprocx,        & !
  icomm_cart,    & ! communicator for the virtual cartesian topology
                   ! that can be used by MPI_WAIT to identify the send
  imp_reals        ! determines the correct REAL type used in the model
                   ! for MPI

USE kind_parameters, ONLY :   &
    wp,             & ! KIND-type parameter for real variables
    sp,             & ! KIND-type parameter for real variables (single precision)
    dp                ! KIND-type parameter for real variables (double precision)

USE data_constants  , ONLY :   &
    pi,             & ! circle constant
    solc,           & ! solar constant
    t0_melt,        & ! melting temperature of ice
    rdv,            & ! r_d / r_v
    o_m_rdv,        & ! 1 - r_d/r_v
    rvd_m_o,        & ! r_v/r_d - 1
    cpdr,           & ! 1 / cp_d
    rdocp,          & ! r_d / cp_d
    lh_v,           & ! latent heat of vapourization
    lhocp,          & ! lh_v/cp_d
    b1,             & ! variables for computing the saturation vapour pressure
    b2w,            & ! over water (w) and ice (i)
    b2i,            & !               -- " --
    b3,             & !               -- " --
    b4w,            & !               -- " --
    b4i,            & !               -- " --
    b234w,          & ! b2w * (b3 - b4w)
    uc1,            & ! variable for computing the rate of cloud cover in 
    uc2,            & ! the unsaturated case
    ucl,            & !               -- " --
    rprecision        !

USE data_fields,     ONLY :   &
    ! co2_scenarios_and_absorbers
    p0hl       ,    & ! reference pressure at half levels             ( Pa  )
    dp0        ,    & ! pressure thickness
    vio3       ,    & ! vertical integrated ozone contents            (pa O3)
    hmo3       ,    & ! ozone maximum                                 ( pa  )
    aer_su     ,    & ! monthly aerosol climatology sulfate drops     (0 - 1)
    aer_du     ,    & ! mon. aerosol climatology total dust           (0 - 1)
    aer_or     ,    & ! mon. aerosol climatology organic (water sol.) (0 - 1)
    aer_bc     ,    & ! mon. aerosol climatology black carbon         (0 - 1)
    aer_ss     ,    & ! mon. aerosol climatology sea salt             (0 - 1)
    aerlan     ,    & ! aerosol-distribution for rural areas            --
    aerurb     ,    & ! aerosol-distribution for urban areas            --
    aerdes     ,    & ! aerosol-distribution for desert areas           --
    aersea     ,    & ! aerosol-distribution for sea                    --

    ! cloudiness_and_humidity: input
    llandmask  ,    & ! landpoint mask
    t_g        ,    & ! weighted surface temperature                  (  k  )
    t          ,    & ! temperature                                   (  k  )
    p0         ,    & ! reference pressure at full levels             ( pa  )
    pp         ,    & ! deviation from the reference pressure         ( pa  )
    ps         ,    & ! surface pressure                              ( pa  )
    rcld       ,    & ! standard deviation of the saturation deficit    --
    clc_con    ,    & ! cloud cover due to convection                   --     

    ! cloudiness_and_humidity: output
    qc_rad     ,    & ! subgrid-scale specific cloud water content    (kg/kg)
    qi_rad     ,    & ! subgrid-scale specific ice water content      (kg/kg)
    clc_sgs    ,    & ! subgrid-scale stratiform cloud cover            --
    clch       ,    & ! cloud cover with high clouds                    --   
    clcm       ,    & ! cloud cover with medium clouds                  --   
    clcl       ,    & ! cloud cover with low clouds                     --   
    clct       ,    & ! total cloud cover                               --   

    ! surface_albedo  (also uses llandmask, t_g, defined above)
    depth_lk   ,    & ! lake depth                                    (  m  )
    for_e      ,    & ! ground fraction covered by evergreen forest     --
    for_d      ,    & ! ground fraction covered by deciduous forest     --
    plcov      ,    & ! fraction of plant cover                         --
    isoiltyp   ,    & ! type of the soil (keys 0-9)                     --
    alb_dry    ,    & ! surface albedo for dry soil                   (     )
    alb_sat    ,    & ! surface albedo for saturated soil             (     )
    alb_dif    ,    & ! solar surface albedo                          (     )
    emis_rad   ,    & ! thermal surface emissivity                    (0 - 1)
    h_ice      ,    & ! lake/sea ice thickness                        (  m  )
    t_ice      ,    & ! temperature of ice/water surface              (  K  )
    t_s        ,    & ! surface temperature                           (  K  )
    t_snow     ,    & ! temperature of the snow-surface               (  k  )
    w_snow     ,    & ! water content of snow                         (m H2O)
    t_snow_mult,    & ! temperature of the snow-surface               (  k  )
    freshsnow  ,    & ! weighting function indicating 'freshness' of snow in
                      ! upper few centimeters of snow cover           ( --- )
    w_so              ! multi-layer soil moisture                     (m H2O)

USE data_flake, ONLY : &
    ! flake_parameters
    h_Ice_min_flk             , & ! Minimum ice thickness [m]
    tpl_T_f                   , & ! Fresh water freezing point [K]
    ! flake_albedo_ref
    albedo_whiteice_ref       , & ! White ice
    albedo_blueice_ref        , & ! Blue ice
    c_albice_MR                   ! Constant in the interpolation formula for

USE data_modelconfig, ONLY :   &
    ie,             & ! number of grid points in zonal direction
    je,             & ! number of grid points in meridional direction
    istartpar,      & ! start index for computations in the parallel program
    iendpar,        & ! end index for computations in the parallel program
    jstartpar,      & ! start index for computations in the parallel program
    jendpar,        & ! end index for computations in the parallel program
    dt,             & ! long time-step
    degrad,         & ! factor for transforming degree to rad
    pollon,         & ! longitude of the rotated north pole (in degrees, E>0)
    pollat,         & ! latitude of the rotated north pole (in degrees, N>0)
    czmls,          & ! depth of the soil layers in meters
    klv800,         & ! k index of the LM-mainlevel, on 800 HPa
    klv500,         & ! k index of the LM-mainlevel, on 500 HPa
    klv400            ! k index of the LM-mainlevel, on 400 HPa

USE data_runcontrol, ONLY :   &
    lreproduce,   & ! the results are reproducible in parallel mode
    nnow,           & ! corresponds to ntstep 
    itype_aerosol,  & ! type of aerosol map internal/external
    icldm_rad,      & ! mode of cloud representation in radiation  parametr.
    lsoil,          & ! forecast with soil model
    lforest,        & ! if .true., run with forest (evergreen and deciduous)
    lemiss,         & ! external emissivity map
    lmulti_snow,    & ! run multi-layer snow model
    llake,          & ! forecst with lake model FLake
    lseaice,        & ! forecast with sea ice model
    itype_albedo,   & ! type of solar surface albedo
    l_cosmo_art,    & ! if .TRUE., run the COSMO_ART
#ifdef TWOMOM_SB
    iradpar_cloud,  & ! type of parameterization of radiative transfer parameters (ext., sing. alb., asym.)  
                      ! for grid- and subgrid scale clouds (cloud water, ice water)
                      !   1 = old method
                      !   2 = based on eff. radius (code from Elias Zubler, ETHZ)
#endif
    lperi_x,        & ! if lartif_data=.TRUE.: periodic boundary conditions in x-dir.
                      !               .FALSE.:   or with Davies conditions
    lperi_y,        & ! if lartif_data=.TRUE.: periodic boundary conditions in y-dir.
                      !               .FALSE.:   or with Davies conditions
    l2dim,          & ! 2 dimensional runs
    lscm,           & ! SCM-run
    itype_calendar,&! for specifying the calendar used
    ico2_rad,       & ! type of CO2 concentration in radiation parameterization
    lco2_stab,      & ! use CO2 stabilisation
    iy_co2_stab       ! default year of CO2 stabilisation


USE data_soil       , ONLY :   &
    csalb         , & !
    ctalb         , & !
    csalb_p       , & !
    csalb_snow    , & !
    csalb_snow_min, & ! min. solar albedo of snow for forest free surfaces
    csalb_snow_max, & ! max. solar albedo of snow for forest free surfaces
    csalb_snow_fd , & ! solar albedo of snow for surfaces with deciduous forest
    csalb_snow_fe , & ! solar albedo of snow for surfaces with evergreen forest
    cf_snow       , & !
    cporv         , & !
    cadp

USE turb_data     , ONLY :   &
    itype_wcld,     & ! type of water cloud diagnosis
    clc_diag,       & ! cloud cover at saturation in statistical cloud diagnostic
    q_crit            ! critical value for normalized over-saturation

USE radiation_data,  ONLY :   &
    rad_csalbw    , & !
    pvdaes        , & ! normalized vertical distribution (sea)
    pvdael        , & ! normalized vertical distribution (land)
    pvdaeu        , & ! normalized vertical distribution (urban)
    pvdaed        , & ! normalized vertical distrubution (desert)
    ptrbga        , & ! b. optical depths divided by pressure (tropospheric)
    pvobga        , & ! b. optical depths divided by pressure (volcanic)
    pstbga        , & ! b. optical depths divided by pressure (stratospheric)
    ptrpt         , & ! temperature exponent for the stratosperic definition
    paeops        , & ! total optical depths for vertical varying aerosols (sea)
    paeopl        , & ! total optical depths for vertical varying aerosols (land)
    paeopu        , & ! total optical depths for vertical varying aerosols (urban)
    paeopd        , & ! total optical depths for vertical varying aerosols (desert)
    mind_ilon_rad , & ! correspondance between blocked (coarse) radiation
    mind_jlat_rad     ! grid and the usual COSMO grid (i,j,k)

USE meteo_utilities,    ONLY :  cloud_diag
USE utilities,          ONLY :  get_utc_date

#ifdef COSMOART
USE data_cosmo_art,     ONLY :     &
    lrad_dust         ! mineral dust aerosols

#ifdef TWOMOM_SB
USE src_cloud_opt_reff, ONLY :   &
   reff_for_rad,                 &
#endif
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================


CONTAINS

!==============================================================================
! Module procedure to compute sunshine conditions
!------------------------------------------------------------------------------

SUBROUTINE compute_sunshine_conditions(                                     &
              ydate_ini, idbg, ista_comp, iend_comp, jsta_comp, jend_comp,  &
              ibc, tsc, tec, rlat, rlon, itaja_zsct_previous,               &
              zsct_save,zdekcos_save, zdeksin_save, zdtzgl_save,            &
              sun_azi,sun_el, zsct, zeit0, zsmu0, nz_cosmu0pos,             &
              zmaxmu02d, zeitrad, zsinphi, zcosphi, lsolar)

!------------------------------------------------------------------------------
!
! Description:
!
!   Compute the sunshine conditions over a certain timestep interval.
!   To compute the conditions for the computation of the radiative
!   equilibrium in fesft, the current and next full radiation timesteps
!   shall be passed as timestep loop index arguments. 
!   The mean and maximum zsmu0 (cosine of the solar zenith angle) are
!   computed. 
!   A solar factor is set, to indicate a sunshine condition per grid point:
!     0: sun below zenith: no solar radiation
!     1: sun above zenith: do solar radiation
!   If the sun is shining at one grid point in a subdomain, the solar radiation
!   is computed for all grid points, but the results are masked out afterwards
!   using the solar factor (which gives reproducible results at the end).
!
!   To compute the conditions at the current timestep for the computation
!   of the radiative fluxes, which is done at every timestep, the current
!   timestep is passed as both start and end index of the timestep loop.
!   In this case, also the sun azimuth and elevation angles are computed. !RUS_BUG?
!
!------------------------------------------------------------------------------

! Input data
! ----------

  CHARACTER (LEN=14), INTENT(IN) :: &
    ydate_ini   ! start of the forecast yyyymmddhh (year, month, day, hour, min, sec)

  INTEGER, INTENT(IN) :: &
    idbg,     &
    tsc, tec, & ! timestep loop indices
    ista_comp,& ! loop start / end indices
    iend_comp,& !
    jsta_comp,& !
    jend_comp,& !
    ibc         ! current block to work on

  REAL (KIND=wp), INTENT(IN) :: &
    rlat(:,:), &
    rlon(:,:)

! Input-output data
! -----------------

  INTEGER, INTENT(INOUT) :: &
    itaja_zsct_previous

  REAL (KIND=wp), INTENT(INOUT) :: &
    zsct_save,        &
    zdtzgl_save,      &
    zdeksin_save,     &
    zdekcos_save       

! Output data
! -----------

  REAL (KIND=wp), INTENT(OUT)   :: &
    sun_azi  (ie,je), &
    sun_el   (ie,je)

  REAL (KIND=wp), INTENT(OUT) :: &
    zsct,             &
    zeit0,            &
    zsmu0    (:,:)      ! Cosine of zenith angle 

  LOGICAL, INTENT(OUT) :: &
    lsolar              ! control switch for solar calculations

! Local arrays, passed through argument list (declared as INTENT(OUT))
! --------------------------------------------------------------------

  INTEGER,        INTENT(OUT) :: &
    nz_cosmu0pos(:,:)

  REAL (KIND=wp), INTENT(OUT) :: &
    zmaxmu02d(:,:)             , &   ! is not really OUT, but only a work array
    zeitrad  (:,:)             , &
    zsinphi  (:,:)             , &
    zcosphi  (:,:)

! Local parameters: 
! ----------------

  REAL    (KIND=wp)       , PARAMETER ::  &
    z_1d7200 = 1._wp/7200._wp ,&
    zcent  = 0.2500_wp, & ! centre weight in a nine point stencil !T.R.
    zside  = 0.1250_wp, & ! weight for side points !T.R.
    zedge  = 0.0625_wp, & ! weight for edge points !T.R.
    zepclc = MAX(1.0E-8_wp,rprecision), &
                          ! avoids cloud cover =1.0 and = 0.0 
    zeph2o = 1.0E-9_wp, & ! minimum value for specific humidity
    zepemu = 1.0E-9_wp, & ! avoids cosine of zenith angle = 0.0
    zclwcm = 1.0E-9_wp, & ! avoids cloud water content = 0.0
    rtod = 57.2957795_wp  ! conversion from radians to degrees

! Local scalars:
! -------------

  CHARACTER (LEN=255) &
    yzerrmsg      ! for error message

  CHARACTER (LEN=14) :: &
    yrad1         ! output from routine get_utc_dat

  CHARACTER (LEN=28) :: &
    yrad2         ! output from routine get_utc_dat

  INTEGER            :: &
    nz_zsct,    &
    nzrad,      &
    itaja,      &
    jj,         &
    i, j,       & ! loop indices over spatial dimensions
    ip, jp,     & ! loop indices over spatial dimensions for blocked data structure
    jmu0,       &
    izerror

  REAL (KIND=wp) :: &
    zsct_h,     &
    zcosthi,    &
    zstunde,    &
    zmaxmu0,    &
    phi_s,      &
    x1, x2,     &
    zsmu0_loc

  REAL (KIND=wp) :: &
    zdek,       &
    ztho,       &
    ztwo,       &
    zsocof

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Start GPU data region
!------------------------------------------------------------------------------

  !$acc data &
  !---- Argument arrays - intent(in)
  !$acc present ( rlat, rlon                              ) &
  !---- Argument arrays - intent(out)
  !$acc present ( sun_azi, sun_el, zsmu0, zmaxmu02d       ) &
  !$acc present ( zeitrad, zsinphi, zcosphi, nz_cosmu0pos ) &
  !---- radiation_data
  !$acc present ( mind_ilon_rad, mind_jlat_rad            )

!------------------------------------------------------------------------------
! Begin subroutine compute_sunshine_condition
!------------------------------------------------------------------------------

  nz_zsct           = 0
  zsct_h            = 0.0_wp
  zmaxmu0           = 0.0_wp

  !$acc parallel
  !$acc loop gang vector collapse(2)
  DO  jp = jsta_comp, jend_comp
    DO  ip = ista_comp, iend_comp
      nz_cosmu0pos(ip,jp) = 0
      zsmu0(ip,jp)        = 0.0_wp
      zmaxmu02d(ip,jp)    = 0.0_wp
    ENDDO
  ENDDO
  !$acc end parallel

  timesteploop: DO jmu0=tsc,tec

    CALL get_utc_date ( jmu0, ydate_ini, dt, itype_calendar, yrad1, yrad2,  &
                        itaja, zstunde )
    READ (yrad1(1:4),'(I4)') jj

    IF (itaja /= itaja_zsct_previous) THEN ! new calendar day
      itaja_zsct_previous = itaja
  
      ztwo    = 0.681_wp + 0.2422_wp*(jj-1949)-(jj-1949)/4
      ztho    = 2.0_wp*pi*( REAL(itaja, wp) -1.0_wp + ztwo )/365.2422_wp
      zdtzgl_save  = 0.000075_wp + 0.001868_wp*COS(   ztho) - 0.032077_wp*SIN(   ztho) &
        - 0.014615_wp*COS(2.0_wp*ztho) - 0.040849_wp*SIN(2.0_wp*ztho)
      zdek    = 0.006918_wp - 0.399912_wp*COS(   ztho) + 0.070257_wp*SIN(   ztho) &
        - 0.006758_wp*COS(2.0_wp*ztho) + 0.000907_wp*SIN(2.0_wp*ztho) &
        - 0.002697_wp*COS(3.0_wp*ztho) + 0.001480_wp*SIN(3.0_wp*ztho)
  
      zdeksin_save = SIN (zdek)
      zdekcos_save = COS (zdek)
  
      zsocof  = 1.000110_wp + 0.034221_wp*COS(   ztho) + 0.001280_wp*SIN(   ztho) &
        + 0.000719_wp*COS(2.0_wp*ztho) + 0.000077_wp*SIN(2.0_wp*ztho)
  
      ! Update the solar constant.
      zsct_save = zsocof*solc

      ! Update the sum of the solar constant over all timesteps.
      zsct_h = zsct_h + zsct_save
      nz_zsct = nz_zsct + 1
    ENDIF

    zstunde = zstunde + dt*z_1d7200 !add half a timestep

    zeit0   = pi*(zstunde-12._wp)/12._wp + zdtzgl_save

    ! Compute zeitrad, zsinphi, and zcosphi for idealized simulations.
    IF ((.NOT. lscm) .AND. lperi_x) THEN
      ! In case of lperi_x=.true., use the solar time of the model reference point
      ! (as implied by pollon,pollat) to avoid boundary disturbances.
      ! The problem is that the "true" sun is not periodic, but has to be
      ! "artificially forced" to be periodic for periodic BCs.

      !$acc parallel
      !$acc loop gang vector collapse(2)
      DO jp=jsta_comp, jend_comp
        DO ip=ista_comp, iend_comp
          IF (pollat >= 0.0_wp) THEN
            zeitrad(ip,jp) = zeit0 + degrad*(pollon-SIGN(1.0_wp,pollon)*180.0_wp)
          ELSE
            zeitrad(ip,jp) = zeit0 + degrad*pollon
          ENDIF
        ENDDO
      ENDDO
      !$acc end parallel

      ! the ELSE case depends on j and is computed within the j-loop
    ENDIF

    IF ((.NOT. lscm) .AND. (lperi_y .OR. l2dim)) THEN
      ! Similar thing for lperi_y=.true. or l2dim=.true.:
      ! Use the geogr. latitude of the model reference point:

      !$acc parallel
      !$acc loop gang vector collapse(2)
      DO jp=jsta_comp, jend_comp
        DO ip=ista_comp, iend_comp
          zsinphi(ip,jp) = SIN (degrad*(90.0_wp-ABS(pollat)))
          zcosphi(ip,jp) = COS (degrad*(90.0_wp-ABS(pollat)))
        ENDDO
      ENDDO
      !$acc end parallel
    ! the ELSE case depends on j and is computed within the j-loop
    ENDIF

    !$acc parallel
    !$acc loop gang vector collapse(2)
    DO  jp = jsta_comp, jend_comp
      DO  ip = ista_comp, iend_comp

        i = mind_ilon_rad(ip,jp,ibc)
        j = mind_jlat_rad(ip,jp,ibc)

        ! Compute zeitrad, zsinphi, and zcosphi for full, non-idealized simulations.
        IF (lscm .OR. (.NOT.lperi_x)) THEN
          zeitrad(ip,jp) = zeit0 + rlon(i,j)
        ENDIF
        IF (lscm .OR. (.NOT.(lperi_y .OR. l2dim))) THEN
          zsinphi(ip,jp) = SIN (rlat(i,j))
!         zcosphi(ip,jp) = COS (rlat(i,j))
        ENDIF

        zcosphi(ip,jp) = SQRT(1.0_wp - zsinphi(ip,jp)**2)
        zcosthi        = zdeksin_save * zsinphi(ip,jp) + zdekcos_save * zcosphi(ip,jp) * COS(zeitrad(ip,jp))

        ! Compute the cosine of the solar zenith angle.
        IF ( zcosthi > zepemu ) THEN
          zsmu0(ip,jp) = zsmu0(ip,jp) + zcosthi ! Update the sum over all timesteps.
          nz_cosmu0pos(ip,jp) = nz_cosmu0pos(ip,jp) + 1
          zmaxmu02d(ip,jp) = MAX (zcosthi, zmaxmu02d(ip,jp))
        ENDIF

      ENDDO
    ENDDO
    !$acc end parallel

  ENDDO timesteploop !jmu0

  ! Compute the mean zsmu0 over all timesteps from the sum determined in the jmu0 loop.
  !$acc parallel
  !$acc loop gang vector collapse(2) reduction(max:zmaxmu0)
  DO  jp = jsta_comp, jend_comp
    DO  ip = ista_comp, iend_comp
      IF ( nz_cosmu0pos(ip,jp) > 0 ) THEN
        zsmu0(ip,jp) = zsmu0(ip,jp) / REAL(nz_cosmu0pos(ip,jp), wp)
      ELSE
        zsmu0(ip,jp) = zepemu
      ENDIF
      zmaxmu0 = MAX (zmaxmu0, zmaxmu02d(ip,jp))
    ENDDO
  ENDDO
  !$acc end parallel

  IF (tsc == tec) THEN
    ! Sun azimuth and sun elevation (for computation of relative sunshine duration) 
    ! are only computed for the update of the zenith angle, which is done in 
    ! every time step. The criterion is tsc==tec (i.e. no radiation interval, 
    ! but only a single time step)

    !$acc parallel
    !$acc loop gang vector collapse(2)
    DO  jp = jsta_comp, jend_comp
      DO  ip = ista_comp, iend_comp
        i = mind_ilon_rad(ip,jp,ibc)
        j = mind_jlat_rad(ip,jp,ibc)

        sun_el(i,j) = ASIN(zsmu0(ip,jp))
        x1 = zdekcos_save * SIN(zeitrad(ip,jp)) / COS(sun_el(i,j))
        x2 = ( SIN(rlat(i,j)) * zdekcos_save * COS(zeitrad(ip,jp)) - &
               COS(rlat(i,j)) * zdeksin_save ) / COS(sun_el(i,j))
        IF (x2 < -1.0_wp) x2 = -1.0_wp
        IF (x2 >  1.0_wp) x2 = 1.0_wp
        phi_s = ACOS(x2)
        IF (x1 < 0) phi_s = - phi_s
        sun_azi(i,j) = rtod*(phi_s + pi)
        sun_el(i,j) = rtod*sun_el(i,j)
      ENDDO
    ENDDO
    !$acc end parallel
  ENDIF


  ! Compute the mean solar constant from the sum determined in the jmu0 loop.
  IF ( nz_zsct > 0 ) THEN
    zsct = zsct_h/REAL(nz_zsct,wp)
  ELSE
    zsct = zsct_save
  ENDIF

#ifdef COSMOART
  !T.R. Should this be shifted after the calculation of the _actual_ zenith angle
  !     (and be done _every_ timestep)?

!US: what to do about the indices here?
! IF(l_cosmo_art) THEN
!   !$acc parallel present_or_copy(mmy)
!   !$acc loop gang
!   DO  jp = jsta_comp, jend_comp
!     !$acc loop vector
!     DO  ip = ista_comp, iend_comp
!       mmy(ip,jp) = zsmu0(ip,jp)
!     ENDDO
!   ENDDO
!   !$acc end parallel
! ENDIF
#endif

  ! Determine whether it's necessary to compute the radiation.
  ! zmaxmu0 now is the maximum over the whole subdomain, so the solar calculations
  ! can be switched off for the whole subdomain, if it is "Night on Earth"
  IF (zmaxmu0 > zepemu) THEN
    lsolar = .TRUE.
  ELSE
    lsolar = .FALSE.
  ENDIF

!------------------------------------------------------------------------------
! End of the subroutine 
!------------------------------------------------------------------------------

  !$acc end data

END SUBROUTINE compute_sunshine_conditions

!==============================================================================
!==============================================================================
!+ Module procedure to calculate solar and thermal surface albedo 
!------------------------------------------------------------------------------

SUBROUTINE surface_albedo              (                                     &
                        ista_comp, iend_comp, jsta_comp, jend_comp, ibc,     &
                        ntl, ralso, ralth)                     

!------------------------------------------------------------------------------
!
! Description:
!  Calculation of surface albedo taking soil type,              
!  vegetation and snow/ice conditions into account
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER, INTENT(IN)         ::  &
    ista_comp,            & ! start and end incides of computation
    iend_comp,            & !
    jsta_comp,            & !
    jend_comp,            & !
    ibc,                  & ! block to compute in blocked data structure
    ntl                     ! time level used for computations

  REAL (KIND=wp), INTENT(OUT) ::  &
    ralso(:,:),           & !
    ralth(:,:)              !


! Local arrays and scalars:
! -------------------------
  INTEGER                   ::                &
    i, j, ist, ip, jp

  REAL (KIND=wp)            ::                &
    t_test, zwetfrac, zsalb_snow, zsnow_alb, zsnow, zvege


!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Start GPU data region
!------------------------------------------------------------------------------

  !$acc data &
  !---- Argument arrays - intent(out)
  !$acc present ( ralso, ralth                                      ) &
  !---- data_fields
  !$acc present ( llandmask, isoiltyp, plcov, emis_rad, depth_lk    ) &
  !$acc present ( alb_dry, alb_sat, alb_dif, for_e, for_d           ) &
  !$acc present ( t_g, t_s, t_snow, t_snow_mult, w_so, w_snow       ) &
  !$acc present ( freshsnow, h_ice, t_ice                           ) &

  !---- data_modelconfig
  !$acc present ( czmls                                             ) &
  !---- data_soil
  !$acc present ( csalb, cporv, cadp                                ) &
  !---- radiation_data
  !$acc present ( mind_ilon_rad, mind_jlat_rad, rad_csalbw          )

!------------------------------------------------------------------------------
! Begin subroutine surface_albedo
!------------------------------------------------------------------------------

  DO jp = jsta_comp, jend_comp       ! loop over nradcoarse points in j-direction
                                     !   is loop from 1 to 1 without coarse grid
    !$acc parallel
    !$acc loop gang vector
    DO  ip = ista_comp, iend_comp    ! loop over grid points in the block in i-direction

      ! get i/j indices for 2D COSMO data structure
      i = mind_ilon_rad(ip,jp,ibc)
      j = mind_jlat_rad(ip,jp,ibc)

      IF (lemiss) THEN
        ralth(ip,jp) = 1.0_wp-emis_rad(i,j)  ! geographical dependent thermal albedo
      ELSE
        ralth(ip,jp) = ctalb
      ENDIF

      ist = 10

      ! In the following IF statement, t_snow has been used up to now.
      ! In NetCDF files, t_snow is undefined (-1E20) where no snow exists.
      ! This leads to ice-points over the whole sea. t_g could be used instead,
      ! but this changes the results and has to be tested more intensively.
      ! As an intermediate solution, we use t_snow, where it is defined,
      ! otherwise t_g (in grib-files, t_snow is defined as t_s, where no snow
      ! exists.

      IF ( lmulti_snow ) THEN
        IF (t_snow_mult(i,j,1,ntl) < 0.0_wp) THEN
          t_test = t_g   (i,j,ntl)
        ELSE
          t_test = t_snow_mult(i,j,1,ntl)
        ENDIF
      ELSE
        IF (t_snow(i,j,ntl) < 0.0_wp) THEN
          t_test = t_g   (i,j,ntl)
        ELSE
          t_test = t_snow(i,j,ntl)
        ENDIF
      ENDIF

      IF ( llandmask(i,j) .OR. t_test >= t0_melt - 1.7_wp ) THEN
        ist = isoiltyp(i,j) ! water (ist=9) and sea ice (ist=10) included
      ENDIF
      ralso(ip,jp) = csalb(ist)

      IF (lsoil .AND. llandmask(i,j)) THEN

        IF     (itype_albedo == 1) THEN
          ralso(ip,jp) = csalb(ist) - rad_csalbw(ist)*w_so(i,j,1,ntl)

        ELSEIF (itype_albedo == 2) THEN
          zwetfrac = (0.5_wp*w_so(i,j,1,ntl)/czmls(1)-cadp(ist))/ &
                     (cporv(ist)-cadp(ist))
          zwetfrac = MIN(1.0_wp,MAX(0.0_wp,zwetfrac))
          ralso(ip,jp) = (1.0_wp-zwetfrac)*alb_dry(i,j) + zwetfrac*alb_sat(i,j)

        ELSEIF (itype_albedo == 3) THEN
          ! some points use soiltyp albedo, predefined in EXTPAR
          ralso(ip,jp) = alb_dif(i,j)  ! MODIS background albedo
        ENDIF

      ENDIF  ! lsoil, llandmask

    ENDDO    ! i-loop
    !$acc end parallel

    IF (lseaice) THEN
      !$acc parallel copyin(nnow)
      !$acc loop gang vector
      DO ip = ista_comp, iend_comp
        ! get i/j indices for 2D COSMO data structure
        i = mind_ilon_rad(ip,jp,ibc)
        j = mind_jlat_rad(ip,jp,ibc)

        ! In case the sea ice model is used AND water point AND ice is present,
        ! compute ice albedo for water points with an empirical formula taken from GME.
        ! The ice albedo is the lower the warmer, and therefore wetter, the ice is.
        ! Use ice temperature at time level nnow (2-time level scheme in sea ice model).

        IF ((.NOT. llandmask(i,j)) .AND. (h_ice(i,j,nnow) > 0.0_wp))               &
          ralso(ip,jp) = (1.0_wp-0.3846_wp*EXP(-0.35_wp*(t0_melt-t_ice(i,j,nnow)))) &
          * csalb(10)
      ENDDO
      !$acc end parallel
    ENDIF

    IF (llake) THEN
      !$acc parallel copyin(nnow)
      !$acc loop gang vector
      DO  ip = ista_comp, iend_comp
        ! get i/j indices for 2D COSMO data structure
        i = mind_ilon_rad(ip,jp,ibc)
        j = mind_jlat_rad(ip,jp,ibc)

        IF((depth_lk(i,j)      >  0.0_wp) .AND.    &
          (h_ice   (i,j,nnow) >= h_Ice_min_flk) ) THEN
          !  In case the lake model FLake is used AND lake point AND ice is present,
          !  compute ice albedo for lake points with an empirical formulation 
          !  proposed by Mironov and Ritter (2004) for use in GME 
          !  [ice_albedo=function(ice_surface_temperature)].
          !  Use surface temperature at time level "nnow".

          ralso(ip,jp) = EXP(-c_albice_MR*(tpl_T_f-t_s(i,j,nnow))/tpl_T_f)
          ralso(ip,jp) = albedo_whiteice_ref * (1._wp-ralso(ip,jp)) +      &
            albedo_blueice_ref  * ralso(ip,jp)
        ENDIF
      ENDDO
      !$acc end parallel
    ENDIF

    ! Snow cover and vegetation
    ! -------------------------

    IF (lsoil) THEN
      !$acc parallel
      !$acc loop gang vector
      DO  ip = ista_comp, iend_comp    ! istartpar, iendpar
        ! get i/j indices for 2D COSMO data structure
        i = mind_ilon_rad(ip,jp,ibc)
        j = mind_jlat_rad(ip,jp,ibc)

        zvege= 0.0_wp
        zsnow= 0.0_wp
        IF (llandmask(i,j)) THEN 
          ! consider effects of aging on solar snow albedo
          zsalb_snow = csalb_snow_min + &
            freshsnow(i,j)*(csalb_snow_max-csalb_snow_min)
          IF (lforest) THEN
            zsnow_alb = zsalb_snow*(1._wp-for_e(i,j)-for_d(i,j))       &
              + csalb_snow_fe * for_e(i,j)                       &
              + csalb_snow_fd * for_d(i,j)
          ELSE
            zsnow_alb = zsalb_snow
          ENDIF

          ! account for snow cover and plant cover and compute final solar
          ! snow albedo
          zvege = plcov(i,j)
          IF (w_snow(i,j,ntl) > 0.0_wp)                              &
            zsnow = MIN(1.0_wp, w_snow(i,j,ntl)/cf_snow)

          IF     ( (itype_albedo == 1) .OR. (itype_albedo == 2) ) THEN
            ralso(ip,jp) = zsnow * zsnow_alb +                             &
              (1.0_wp - zsnow) * (zvege * csalb_p + (1.0_wp - zvege) * ralso(ip,jp))

          ELSEIF   (itype_albedo == 3) THEN
            ralso(ip,jp) = zsnow * zsnow_alb +                             &
              (1.0_wp - zsnow) *                                           ralso(ip,jp)

          ELSEIF   (itype_albedo == 4) THEN
!US for_e, for_d are used, but lforest is not checked???
            ralso(ip,jp)= zsnow*zsnow_alb +                                &
                     (1._wp-zsnow)*                                  &
                     ( zvege*(for_e(i,j)*0.10_wp     +               &
                              for_d(i,j)*0.15_wp     +               &
                     (1._wp-for_e(i,j)-for_d(i,j))*0.20_wp)+     &
                     (1._wp-zvege)*ralso(ip,jp))
          ENDIF

        ENDIF !   llandmask
      ENDDO
      !$acc end parallel
    ENDIF
  ENDDO

  !$acc end data

END SUBROUTINE surface_albedo

!==============================================================================
!==============================================================================
!+ Module procedure to calculate cloudiness and humidity
!------------------------------------------------------------------------------

SUBROUTINE cloudiness_and_humidity     (qv, qc, qi,                            &
                        ista_comp, iend_comp, jsta_comp, jend_comp, ibc, kend, &
                        ntl, lprog_qi, rsw, rwv, rclwc, rciwc, rclc,           &
                        zse, zclcmax, zclcmin, zclcm1)

!------------------------------------------------------------------------------
!
! Description:
!  Depending on the namelist variable icldm_rad, the cloud cover is computed
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER, INTENT(IN)         ::  &
    ista_comp,            & ! start and end incides of computation
    iend_comp,            & !
    jsta_comp,            & !
    jend_comp,            & !
    ibc,                  & ! block to compute
    kend,                 & ! vertical dimension of fields
    ntl                     ! time level used for computations

  LOGICAL, INTENT(IN)         ::  &
    lprog_qi

  REAL (KIND=wp), INTENT(IN)  ::  &
    qv(:,:,:),            & ! water vapour
    qc(:,:,:),            & ! cloud water
    qi(:,:,:)               ! cloud ice

  REAL (KIND=wp), INTENT(OUT) ::  &
    rsw  (:,:,:),         & !
    rwv  (:,:,:),         & !
    rclwc(:,:,:),         & !
    rciwc(:,:,:),         & !
    rclc (:,:,:)            !


! Local arrays, passed through argument list (declared as INTENT(OUT))
! --------------------------------------------------------------------

  REAL (KIND=wp), INTENT(OUT) ::  &
    zse    (:,:,:),         & !
    zclcmax(:,:,:),         & !
    zclcmin(:,:,:),         & !
    zclcm1 (:,  :)

! Local scalars:
! --------------
  INTEGER                   ::                &
    i, j, k  , ip, jp


  REAL (KIND=wp)            ::                &
    fgew, fgee, fgqv,          & ! name of statement functions
    ztt,  zzpv, zzpa,          & ! dummy arguments of stat. func.
    zclwfs, zclwcs,            & !
    zclics, zclws ,            & !
    zclick, zclwck,  zclwk,    & !
    zclwfk, zck, zcs, zdthdz,  & !
    zqdw, zsex, zf_ice, zpim,  & !
    zpio, zpiu, zsigma, zthvo, & !
    zthvu, zuc,                & !
    zt_ice1, zt_ice2, zph

  REAL (KIND=wp), PARAMETER ::  &
    zepclc = MAX(1.0E-8_wp,rprecision), & ! avoids cloud cover =1.0 and = 0.0 
    zeph2o = 1.0E-9_wp,                 & ! minimum value for specific humidity
    zclwcm = 1.0E-9_wp                    ! avoids cloud water content = 0.0

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine cloudiness and humidity
!------------------------------------------------------------------------------

! statement function to calculate saturation vapour pressure over water
  fgew(ztt)       = b1 * EXP( b2w*(ztt - b3)/(ztt - b4w) ) ! ztt: temperature

! statement function to calculate saturation vapour pressure over ice
  fgee(ztt)       = b1 * EXP( b2i*(ztt - b3)/(ztt - b4i) ) ! ztt: temperature

! statement function to calculate specific humitdity  
  fgqv(zzpv,zzpa) = rdv*zzpv/(zzpa - o_m_rdv*zzpv)   ! zzpv: vapour pressure

!------------------------------------------------------------------------------
! Start GPU data region
!------------------------------------------------------------------------------

  !$acc data &
  !---- Argument arrays - intent(in)
  !$acc present ( qv, qc, qi                                        ) &
  !---- Argument arrays - intent(out)
  !$acc present ( rsw, rwv, rclwc, rciwc, rclc                      ) &
  !$acc present ( zse, zclcmax, zclcmin, zclcm1                     ) &
  !---- data_fields
  !$acc present ( llandmask, p0, pp, ps, t, t_g, clc_sgs, clc_con  ) &
  !$acc present ( qc_rad, qi_rad, clch, clct, clcm, clcl           ) &
  !---- radiation_data
  !$acc present ( mind_ilon_rad, mind_jlat_rad                      )

!------------------------------------------------------------------------------
! Section 1: Calculate water vapour saturation mixing ratios over water and ice 
!------------------------------------------------------------------------------

  ! maximum (in-)cloud water content: 0.5% of specific humidity at saturation
  zclwfs=0.005_wp !0.5% of specific humidity at saturation in stratiform clouds
  zclwfk=0.010_wp !1.0% of specific humidity at saturation in convective clouds

  zt_ice1= t0_melt -  5.0_wp
  zt_ice2= t0_melt - 25.0_wp

  !$acc parallel
  !$acc loop gang vector collapse(3)
  DO k = 1, kend
    DO jp = jsta_comp, jend_comp        ! from 1 to 1 without coarse grid
      DO  ip = ista_comp, iend_comp

        ! get i/j indices for 2D COSMO data structure
        i = mind_ilon_rad(ip,jp,ibc)
        j = mind_jlat_rad(ip,jp,ibc)

        ! specific humidity (rwv) specific total water content (zqdw), 
        ! specific humidity at saturation
        ! over water (rsw ) and ice (zse)
        zph       = p0(i,j,k) + pp(i,j,k,ntl)
        zse (ip,k,jp) = fgqv ( fgee(t(i,j,k,ntl)), zph)
        rsw (ip,k,jp) = fgqv ( fgew(t(i,j,k,ntl)), zph)
      ENDDO
    ENDDO
  ENDDO
  !$acc end parallel

!------------------------------------------------------------------------------
! Section 2: Calculate stratiform cloud cover (non-convective)
!------------------------------------------------------------------------------

  IF     ( icldm_rad == 0 ) THEN

    !--------------------------------------------------------------------------
    ! Section 2.0: No interpretation of clouds at all for radiative calculations
    !--------------------------------------------------------------------------

    !$acc parallel
    !$acc loop gang vector collapse(3)
    DO  k = 1, kend
      DO jp = jsta_comp, jend_comp
        DO  ip = ista_comp, iend_comp

          rclwc(ip,k,jp) = 0.0_wp
          rciwc(ip,k,jp) = 0.0_wp
          rclc (ip,k,jp) = 0.0_wp
        ENDDO
      ENDDO
    ENDDO
    !$acc end parallel

  ELSEIF ( icldm_rad == 1 ) THEN

    !--------------------------------------------------------------------------
    ! Section 2.1: Only grid-sale water clouds are passed to the radiation routine
    !--------------------------------------------------------------------------

    !$acc parallel
    !$acc loop gang vector collapse(3)
    DO  k = 1, kend
      DO jp = jsta_comp, jend_comp
!CDIR NODEP
        DO  ip = ista_comp, iend_comp

          ! get i/j indices for 2D COSMO data structure
          i = mind_ilon_rad(ip,jp,ibc)
          j = mind_jlat_rad(ip,jp,ibc)

          rclwc(ip,k,jp) = qc(i,j,k)
          IF (lprog_qi) THEN
            IF ( qc(i,j,k)+qi(i,j,k) > 0.0_wp ) THEN
              rclc(ip,k,jp) = 1.0_wp
            ELSE
              rclc(ip,k,jp) = 0.0_wp
            END IF
            rciwc(ip,k,jp) = qi(i,j,k)
          ELSE
            IF ( qc(i,j,k) > 0.0_wp ) THEN
              rclc(ip,k,jp) = 1.0_wp
            ELSE
              rclc(ip,k,jp) = 0.0_wp
            END IF
            rciwc(ip,k,jp) = 0.0_wp
          ENDIF
          clc_sgs(i,j,k) = rclc(ip,k,jp)
        ENDDO
      ENDDO
    ENDDO
    !$acc end parallel

  ELSEIF (icldm_rad == 2) THEN

    !--------------------------------------------------------------------------
    ! Section 2.2: Cloud cover and water content from statistical diagnosis
    !--------------------------------------------------------------------------

!   CALL cloud_diag(rclc,rclwc,                                   &
!     ista_comp, iend_comp, jsta_comp, jend_comp, 1, kend,        &
!     1, ie, 1, je,1, kend,                                       &
!     ie, je, kend, kend+1,                                       &
!     rdv,o_m_rdv,rvd_m_o,lhocp,t0_melt,                          &
!     b1,b2w,b3,b4w,b234w,b2i,b4i,                                &
!     uc1,uc2,ucl, clc_diag, q_crit,                              &
!     t(:,:,:,ntl),qv(:,:,:),qc(:,:,:),                           &
!     pp(:,:,:,ntl)+p0(:,:,:),rcld,ps(:,:,ntl),                   &
!     itype_wcld)

!   DO  k = 1, kend
!     DO j = jsta_comp, jend_comp        ! jstartpar, jendpar
!CDIR NODEP
!       DO  i = ista_comp, iend_comp     ! istartpar, iendpar
!         ! convective (in-)cloud water content 
!         ! as a function of specific humidity at saturation
!         IF ( t(i,j,k,ntl) >= t0_melt )  THEN
!           zclwck = rsw(i,k,j)*zclwfk  !                           
!         ELSE
!           zclwck = zse(i,k,j)*zclwfk
!         ENDIF

!         ! cloud cover of the non convective part of the grid box and cloud ice
!         zcs = rclc(i,k,j)
!         rciwc(i,k,j) = 0.0_wp

!         IF (lprog_qi) THEN
!           ! if there is a grid scale cloud with cloud ice,
!           ! even there might has been diagnosed subgrid scale water clouds,
!           ! their water is thought to be distributed over the 
!           ! whole grid volume:
!           IF ( qi(i,j,k) > 0.0_wp ) THEN
!             zcs = 1.0_wp
!           ENDIF
!           rciwc(i,k,j) = qi(i,j,k)
!         ENDIF
!         clc_sgs(i,j,k) = zcs

!         ! convective cloud cover
!         zck = clc_con(i,j,k)

!         ! grid scale cloud cover and water contend
!         rclc (i,k,j) = zcs + zck*(1.0_wp-zcs)
!         rclwc(i,k,j) = rclwc(i,k,j)*(1.0_wp-zck) + zclwck*zck
!       ENDDO
!     ENDDO
!   ENDDO

  ELSEIF ( icldm_rad == 4 .OR. icldm_rad == 3 ) THEN

    !--------------------------------------------------------------------------
    ! Section 2.3/4: Standard diagnosis
    !--------------------------------------------------------------------------

    !$acc parallel
    !$acc loop gang vector collapse(3)
    DO  k = 1, kend
      DO jp = jsta_comp, jend_comp
!CDIR NODEP
        DO  ip = ista_comp, iend_comp

          ! get i/j indices for 2D COSMO data structure
          i = mind_ilon_rad(ip,jp,ibc)
          j = mind_jlat_rad(ip,jp,ibc)

          ! Critical relative humidity as function of thermal stability
          zph      = p0(i,j,k) + pp(i,j,k,ntl)
          zsigma   = zph / ps(i,j,ntl)
          zdthdz   = 0.0_wp

          zsex      = rsw(ip,k,jp)
          zqdw      = qv(i,j,k) + qc(i,j,k)
          IF (lprog_qi) THEN
            zf_ice      = 1.0_wp - MIN( 1.0_wp, MAX( 0.0_wp, &
              (t(i,j,k,ntl)-zt_ice2)/(zt_ice1-zt_ice2) ) )
            zqdw        = zqdw      + qi(i,j,k)
            zsex        = rsw(ip,k,jp) * (1.0_wp - zf_ice) + zse(ip,k,jp)*zf_ice
          ENDIF

          IF(k == kend) THEN
            zpio    = ( 1.0E-5_wp *( p0(i,j,k)+pp(i,j,k,ntl) ) )**rdocp
            zpiu    = ( 1.0E-5_wp * ps(i,j,ntl)  )**rdocp
            zpim    = 0.5_wp*(zpio+zpiu)
            zthvo   = t  (i,j,k  ,ntl)/zpio
            zthvu   = t_g(i,j,    ntl)/zpiu
            zdthdz  = zthvo - zthvu
          ELSE IF(zsigma > 0.95_wp) THEN
            zpio    = ( 1.0E-5_wp *( p0(i,j,k  )+pp(i,j,k  ,ntl) ) )**rdocp
            zpiu    = ( 1.0E-5_wp *( p0(i,j,k+1)+pp(i,j,k+1,ntl) ) )**rdocp
            zpim    = 0.5_wp*(zpio+zpiu)
            zthvo   = t(i,j,k  ,ntl)/zpio
            zthvu   = t(i,j,k+1,ntl)/zpiu
            zdthdz  = zthvo - zthvu + (lh_v*cpdr/zpim)*(qv(i,j,k)-qv(i,j,k+1))
          ENDIF

          ! grid scale cloud cover as function of relative humidity
          zuc     = 0.95_wp - uc1*zsigma*(1.0_wp-zsigma)*(1.0_wp+uc2*(zsigma-0.5_wp))
          zcs     = MAX ( 0.0_wp, MIN ( 1.0_wp, (zqdw/zsex-zuc)/(ucl-zuc) ) )**2

          ! Corrections and limitations
          IF ( (zsigma > 0.95_wp) .AND. (zdthdz < 0.0_wp) ) THEN
            zcs = 0.0_wp  ! no cloud cover in unstable stratification
          ENDIF
          IF ( qc(i,j,k) > 0.0_wp ) THEN  ! grid scale clouds
            IF ( llandmask(i,j) .AND. k < kend ) zcs = 1.0_wp
          ENDIF
          IF (lprog_qi) THEN
            IF (qi(i,j,k) > 1.0E-7_wp) THEN
              zcs = 1.0_wp ! grid scale clouds with cloud ice
            ENDIF
          ENDIF

          ! store grid-scale cloud cover on global array
          clc_sgs(i,j,k) = zcs

          ! Maximum in-cloud water content:  1.0% of specific humidity at saturation
          !                                  except for convective clouds (fixed)
          ! Standard diagnosis

          IF (lprog_qi) THEN
            zclws  = 0.005_wp*zsex
            zclwcs = zclws*(1.0_wp-zf_ice)
            zclics = zclws*zf_ice

            ! Check for grid-scale water or ice-clouds
            ! Now change in zclwcs only if qc(i,j,k,ntl) > 0.0
            IF ( qc(i,j,k) > 0.0_wp ) THEN      ! grid scale cloud water
              zclwcs = MAX( zclwcs, 0.5_wp*qc(i,j,k) )
            ENDIF
            ! Now change in zclics only if qi(i,j,k,ntl) > 1.0E-7
            IF ( qi(i,j,k) > 1.0E-7_wp ) THEN  ! grid scale cloud ice
              zclics = MAX( zclics, 0.5_wp*qi(i,j,k) )
            ENDIF

            ! Convective cloud water / ice content
            zclwk  = MAX( 2.0_wp*zclws, 0.0002_wp )
            zclwck = zclwk*(1.0_wp-zf_ice)
            zclick = zclwk*zf_ice

            ! Reduce the cloud cover of ice clouds in the upper troposphere
            ! for the diagnosis of clch and clct
            ! Changed by Axel: 
            IF ( (zclwcs <= 1.0E-10_wp) .AND. (zclics  > 0.0_wp) ) THEN
              clc_sgs(i,j,k) = clc_sgs(i,j,k)*MIN( 1._wp, MAX(0.0_wp, &
                                 ( LOG(zclics)       - LOG(1.E-7_wp) )/ &
                                 ( LOG(8.E-6_wp) - LOG(1.E-7_wp) )) )
            ENDIF
            !IF ((k <= klv500) .AND. (zclwcs <= 1.0E-10_wp) .AND. &
            !  (zclics  > 0.0_wp) ) THEN
            !  clc_sgs(i,j,k) = clc_sgs(i,j,k)*MIN( 1._wp, MAX(0.2_wp, &
            !    ( LOG(zclics)       - LOG(1.E-7_wp) )/ &
            !    ( LOG(5.E-5_wp) - LOG(1.E-7_wp) )) )
            !ENDIF

            ! set area-average cloud water/ice content
            rclwc(ip,k,jp) = zclwck*clc_con(i,j,k) + &
              zclwcs*clc_sgs(i,j,k)*(1.0_wp-clc_con(i,j,k))
            rciwc(ip,k,jp) = zclick*clc_con(i,j,k) + &
              zclics*clc_sgs(i,j,k)*(1.0_wp-clc_con(i,j,k))

          ELSE

            zclwcs = 0.005_wp*rsw(ip,k,jp)
            zclwck = MAX( zclwcs, 0.0002_wp )
            IF ( qc(i,j,k) > 0.0_wp ) THEN  ! grid scale clouds
              zclwcs = MAX( zclwcs, 0.5_wp*qc(i,j,k) )
            ENDIF

            ! set area-average cloud water/ice content
            rclwc(ip,k,jp) = zclwck*clc_con(i,j,k) + &
              zclwcs*clc_sgs(i,j,k)*(1.0_wp-clc_con(i,j,k))

            ! set area average cloud ice content (in-cloud)
            rciwc(ip,k,jp) = 0.0_wp
          ENDIF

          ! calculate combined cloud cover
          rclc (ip,k,jp) = clc_sgs(i,j,k) + &
              clc_con(i,j,k)*( 1.0_wp - clc_sgs(i,j,k) )

        ENDDO
      ENDDO
    ENDDO
    !$acc end parallel

  ENDIF  ! icldm_rad

  ! Restrictions for radiative calculations
  ! ------------------------------------

  !$acc parallel
  !$acc loop gang vector collapse(3)
  DO  k = 1, kend
    DO jp = jsta_comp, jend_comp
!CDIR NODEP
      DO  ip = ista_comp, iend_comp

        ! get i/j indices for 2D COSMO data structure
        i = mind_ilon_rad(ip,jp,ibc)
        j = mind_jlat_rad(ip,jp,ibc)

        rwv  (ip,k,jp) = MIN( MAX(zeph2o,qv(i,j,k)), rsw(ip,k,jp) )
        rclc (ip,k,jp) = MAX( zepclc, MIN(1.0_wp-zepclc,rclc(ip,k,jp)) )
        rclwc(ip,k,jp) = MAX( zclwcm, rclwc(ip,k,jp) )
        rciwc(ip,k,jp) = MAX( zclwcm, rciwc(ip,k,jp) )

        ! set qc_rad, qi_rad
        qc_rad(i,j,k) = rclwc(ip,k,jp)
        qi_rad(i,j,k) = rciwc(ip,k,jp)
      ENDDO
    ENDDO
  ENDDO
  !$acc end parallel

#if defined TWOMOM_SB && defined COSMO_ART
  ! Calculate cloud optical properties based on cloud effective radii
  ! ------------------------------------------------------------------
  IF(iradpar_cloud > 1) THEN
    CALL reff_for_rad(rclc,rciwc,iradpar_cloud)
  ENDIF
#endif

  !----------------------------------------------------------------------------
  ! Section 3:  Calculate and store total cloud cover for 3 integral layers
  !             (high, medium, low)
  !----------------------------------------------------------------------------

  DO jp = jsta_comp, jend_comp      ! loop from 1 to 1 without coarse radiation

    !$acc parallel
    !$acc loop gang vector
    DO  ip = ista_comp, iend_comp
      zclcm1(ip,jp) = 1.0_wp - rclc(ip,1,jp)
    ENDDO
    !$acc end parallel

    !$acc parallel
    !$acc loop gang vector collapse(2)
    DO  k  = 2, kend
      DO ip = ista_comp, iend_comp
        zclcmax(ip,k,jp) = 1.0_wp - MAX(rclc(ip,k,jp), rclc(ip,k-1,jp))
        zclcmin(ip,k,jp) = 1.0_wp / (1.0_wp - rclc(ip,k-1,jp))
      ENDDO
    ENDDO
    !$acc end parallel

    !$acc parallel
    !$acc loop seq
    DO  k  = 2, klv400
      !$acc loop gang vector
      DO ip = ista_comp, iend_comp
        zclcm1(ip,jp) = zclcm1(ip,jp) * zclcmax(ip,k,jp) * zclcmin(ip,k,jp)
      ENDDO
    ENDDO
    !$acc end parallel

    !$acc parallel
    !$acc loop gang vector
    DO  ip = ista_comp, iend_comp
      ! get i/j indices for 2D COSMO data structure
      i = mind_ilon_rad(ip,jp,ibc)
      j = mind_jlat_rad(ip,jp,ibc)
      clch  (i,j) = 1.0_wp - zclcm1(ip,jp) - zepclc
    ENDDO
    !$acc end parallel

    !$acc parallel
    !$acc loop seq
    DO  k  = klv400+1, kend
      !$acc loop gang vector
      DO  ip = ista_comp, iend_comp
        zclcm1(ip,jp) = zclcm1(ip,jp)* zclcmax(ip,k,jp) * zclcmin(ip,k,jp)
      ENDDO
    ENDDO
    !$acc end parallel

    !$acc parallel
    !$acc loop gang vector
    DO  ip = ista_comp, iend_comp
      ! get i/j indices for 2D COSMO data structure
      i = mind_ilon_rad(ip,jp,ibc)
      j = mind_jlat_rad(ip,jp,ibc)
      clct  (i,j)   = 1.0_wp - zclcm1(ip,jp) - zepclc
      zclcm1(ip,jp) = 1.0_wp - rclc(ip,klv400+1,jp)
    ENDDO
    !$acc end parallel

    !$acc parallel
    !$acc loop seq
    DO  k  = klv400+2,klv800
      !$acc loop gang vector
      DO  ip = ista_comp, iend_comp
        zclcm1(ip,jp) = zclcm1(ip,jp) * zclcmax(ip,k,jp) * zclcmin(ip,k,jp)
      ENDDO
    ENDDO
    !$acc end parallel

    !$acc parallel
    !$acc loop gang vector
    DO  ip = ista_comp, iend_comp
      ! get i/j indices for 2D COSMO data structure
      i = mind_ilon_rad(ip,jp,ibc)
      j = mind_jlat_rad(ip,jp,ibc)
      clcm  (i,j)   = 1.0_wp - zclcm1(ip,jp) - zepclc
      zclcm1(ip,jp) = 1.0_wp - rclc(ip,klv800+1,jp)
    ENDDO
    !$acc end parallel

    !$acc parallel
    !$acc loop seq
    DO  k  = klv800+2,kend
      !$acc loop gang vector
      DO  ip = ista_comp, iend_comp
        ! get i/j indices for 2D COSMO data structure
        i = mind_ilon_rad(ip,jp,ibc)
        j = mind_jlat_rad(ip,jp,ibc)
        zclcm1(ip,jp) = zclcm1(ip,jp) * zclcmax(ip,k,jp) * zclcmin(ip,k,jp)
      ENDDO
    ENDDO
    !$acc end parallel

    !$acc parallel
    !$acc loop gang vector
    DO  ip = ista_comp, iend_comp
      ! get i/j indices for 2D COSMO data structure
      i = mind_ilon_rad(ip,jp,ibc)
      j = mind_jlat_rad(ip,jp,ibc)
      clcl  (i,j) = 1.0_wp - zclcm1(ip,jp) - zepclc
    ENDDO
    !$acc end parallel

  ENDDO

  !$acc end data

END SUBROUTINE cloudiness_and_humidity

!==============================================================================
!==============================================================================
!+ Module procedure to handle CO2 scenarios and aerosols for radiation
!------------------------------------------------------------------------------

SUBROUTINE co2_scenarios_and_absorbers (iday_of_year, iyear,                   &
                        ista_comp, iend_comp, jsta_comp, jend_comp, ibc, kend, &
                        idebug, rti,                                           &
                        rduco2f, rduo3f, raeq1, raeq2, raeq3, raeq4, raeq5,    &
                        zqcfo, zo3h, zqofo, zaeqdo, zaequo, zaeqlo, zaeqso,    &
                        zaetr_top)
                
!------------------------------------------------------------------------------
!
! Description:
!  This subroutine chooses a CO2 scenario and computes amounts of absorbers
!  (aerosols, CO2, O3)
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER, INTENT(IN)         ::  &
    iday_of_year,         & ! number of the day of the year
    iyear,                & ! actual year
    ista_comp,            & ! start and end incides of computation
    iend_comp,            & !
    jsta_comp,            & !
    jend_comp,            & !
    ibc,                  & ! block to compute
    kend,                 & ! vertical dimension of fields
                            ! (varies for regular or coarse radiation grid)
    idebug                  ! for debug output


  REAL (KIND=wp), INTENT(IN)  ::  &
    rti    (:,:,:)          ! temperature at layer boundaries

  REAL (KIND=wp), INTENT(OUT) ::  &
    rduco2f(:,:,:),       & !
    rduo3f (:,:,:),       & !
    raeq1  (:,:,:),       & !
    raeq2  (:,:,:),       & !
    raeq3  (:,:,:),       & !
    raeq4  (:,:,:),       & !
    raeq5  (:,:,:)          !


! Local arrays, passed through argument list (declared as INTENT(OUT))
! --------------------------------------------------------------------

  REAL (KIND=wp), INTENT(OUT) ::  &
    zqcfo    (:,:),               &
    zo3h     (:,:),               &
    zqofo    (:,:),               &
    zaeqdo   (:,:),               &
    zaequo   (:,:),               &
    zaeqlo   (:,:),               &
    zaeqso   (:,:),               &
    zaetr_top(:,:)

! Local scalars:
! --------------
  INTEGER         ::        &
    i, j, k  , ip, jp

  REAL (KIND=wp)  ::        &
    zqcfn, zqofn, zaeqdn, zaequn, zaeqln, zaeqsn, zaetr_bot, zaetr,  &
    zyear, zyearmax, zyearmin, zqco2, zdpo, zdpn

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Start GPU data region
!------------------------------------------------------------------------------

  !$acc data &
  !---- Argument arrays - intent(in)
  !$acc present ( rti                                               ) &
  !---- Argument arrays - intent(out)
  !$acc present ( rduco2f, rduo3f, raeq1, raeq2, raeq3, raeq4, raeq5) &
  !$acc present ( zqcfo, zo3h, zqofo, zaeqdo, zaequo, zaeqlo, zaeqso) &
  !$acc present ( zaetr_top                                         ) &
  !---- data_fields
  !$acc present ( dp0, p0hl, vio3, hmo3                            ) &
  !$acc present ( aersea, aerlan, aerurb, aerdes                   ) &
  !$acc present ( aer_su, aer_du, aer_or, aer_bc, aer_ss           ) &
  !---- radiation_data
  !$acc present ( mind_ilon_rad, mind_jlat_rad                     ) &
  !$acc present ( pvdaes, pvdael, pvdaeu, pvdaed                   )

!------------------------------------------------------------------------------
! Section 1: Choose a CO2 scenario
!------------------------------------------------------------------------------

  IF (idebug > 10) THEN
    PRINT *, '           choose CO2 scenario and compute amount of absorbers'
  ENDIF

  IF (idebug >= 2 .AND. lco2_stab) THEN
    PRINT *, '           CO2 stabilisation active!'
  ENDIF

  ! Set zqco2, dependent on the chosen CO2-Scenario
  zyear = REAL(iyear,wp) + REAL(iday_of_year,wp)/365.0_wp

  ! Define upper limits of fitted scenarios polynoms to avoid eloping
  zyearmin = REAL (1950, wp)
  IF (ico2_rad < 7) THEN ! SRES fitting
    zyearmax = REAL (2100, wp)
  ELSE                   ! RCP fitting valid until 2150
    zyearmax = REAL (2150, wp)
  ENDIF

  ! zyear can be modified, because it is only used for calculation of CO2 trends
  IF (zyear < zyearmin) THEN
    zyear = zyearmin
  ELSEIF (zyear > zyearmax) THEN
    zyear = zyearmax
  ENDIF

  ! CO2 stabilisation
  IF (lco2_stab) THEN
    IF (zyear > REAL(iy_co2_stab,wp)) THEN
      zyear = REAL(iy_co2_stab,wp)
    ENDIF
  ENDIF

  SELECT CASE (ico2_rad)
  CASE (0)
    ! specific CO2 content of the atmosphere (=330 PPM) (default for DWD)
    zqco2  = 0.5014E-3_wp

  CASE (1)
    ! time dependent CO2 content (fits of IPCC scenario values, taken from ECHAM5)
    ! A1B scenario (for 1950 <= zyear <= 2100)
    !   only CO2
    zqco2 = (- 2.2915249519766070E+07_wp                &
      + 45714.032150104744_wp      * zyear       &
      - 34.178190922262594_wp      * zyear*zyear &
      + 0.01134997524110402_wp     * zyear**3    &
      - 1.4124678138498344E-06_wp  * zyear**4) * 1.519E-06_wp

  CASE (2)
    ! A1B scenario (for 1950 <= zyear <= 2100)
    !   eff. CO2 (i.e. CO2 & CH4 & N2O)
    zqco2 = ( -2.131843263017098E+07_wp                &
      + 42697.69425574343_wp      * zyear       &
      - 32.04969808544885_wp      * zyear*zyear &
      + 0.010685253016710392_wp   * zyear**3    &
      - 1.3349801856070718E-06_wp * zyear**4) * 1.519E-06_wp

  CASE (3)
    ! B1 scenario (for 1950 <= zyear <= 2100)
    !   only CO2
    zqco2 = (- 1.0401357268181011E+07_wp               &
      + 21152.707545487563_wp     * zyear       &
      - 16.116691528852456_wp     * zyear*zyear &
      + 0.005452554505141226_wp   * zyear**3    &
      - 6.910849734430986E-07_wp  * zyear**4) * 1.519E-06_wp

  CASE (4)
    ! B1 scenario (for 1950 <= zyear <= 2100)
    !   eff. CO2 (i.e. CO2 & CH4 & N2O)
    zqco2 = (- 7.716609874947305E+06_wp                &
      + 15881.335647116388_wp     * zyear       &
      - 12.239258629216023_wp     * zyear*zyear &
      + 0.0041862325463834565_wp  * zyear**3    &
      - 5.361489502050553E-07_wp  * zyear**4) * 1.519E-06_wp

  CASE (5)
    ! A2 scenario (for 1950 <= zyear <= 2100)
    !   only CO2
    zqco2 = (  3.682237592956851E06_wp                 &
      - 7547.069807360021_wp      * zyear       &
      + 5.8133367065151145_wp     * zyear*zyear &
      - 0.001994454601121309_wp   * zyear**3    &
      + 2.571600007798381E-07_wp  * zyear**4 ) * 1.519E-06_wp

  CASE (6)
    ! A2 scenario (for 1950 <= zyear <= 2100)
    !   eff. CO2 (i.e. CO2 & CH4 & N2O)
    zqco2 = ( - 340960.0590212098_wp                    &
      + 403.20639583857496_wp     * zyear       &
      - 0.074859345260926_wp      * zyear*zyear &
      - 0.00005743139714985962_wp * zyear**3    &
      + 1.837122734626407E-08_wp  * zyear**4) * 1.519E-06_wp

  CASE (7)
    ! RCP2.6 scenario (for 1950 <= zyear <= 2150)
    !   eff. CO2 (all GreeHouseGases (GHG) considered)
    zqco2 = ( + 5.8284208232E+08_wp        &
              - 1.4124858918E+06_wp * zyear       &
              + 1.3686382349E+03_wp * zyear*zyear &
              - 6.6279390807E-01_wp * zyear**3    &
              + 1.6041979084E-04_wp * zyear**4    &
              - 1.5524630613E-08_wp * zyear**5) * 1.519E-06_wp

  CASE (8)
    ! RCP4.5 scenario (for 1950 <= zyear <= 2150)
    !   eff. CO2 (all GreeHouseGases (GHG) considered)
    zqco2 = ( + 1.9955662739E+07_wp               &
              - 3.8768658589E+04_wp * zyear       &
              + 2.8220059919E+01_wp * zyear*zyear &
              - 9.1219963074E-03_wp * zyear**3    &
              + 1.1048642039E-06_wp * zyear**4) * 1.519E-06_wp

  CASE (9)
    ! RCP6 scenario (for 1950 <= zyear <= 2150)
    !   eff. CO2 (all GreeHouseGases (GHG) considered)
    zqco2 = ( - 2.1182177462E+07_wp               &
              + 4.1828994948E+04_wp * zyear       &
              - 3.0962236444E+01_wp * zyear*zyear &
              + 1.0181182525E-02_wp * zyear**3    &
              - 1.2547431825E-06_wp * zyear**4) * 1.519E-06_wp

  CASE (10)
    ! RCP8.5 scenario (for 1950 <= zyear <= 2150)
    !   eff. CO2 (all GreeHouseGases (GHG) considered)
    zqco2 = ( - 4.0501145412E+07_wp               &
              + 7.9386473439E+04_wp * zyear       &
              - 5.8292720579E+01_wp * zyear*zyear &
              + 1.9002921793E-02_wp * zyear**3    &
              - 2.3202412328E-06_wp * zyear**4) * 1.519E-06_wp

  END SELECT

!------------------------------------------------------------------------------
! Section 2:  Calculate amounts of absorbers (CO2, O3, Aerosol)
!------------------------------------------------------------------------------

#ifdef COSMOART
  ! Change climatology for dust
  IF (l_cosmo_art .AND. lrad_dust) THEN
    DO jp = jsta_comp, jend_comp           ! loop over (coarse) points in j-direction
      !$acc parallel
      !$acc loop gang vector
      DO  ip = ista_comp, iend_comp        ! loop in i-direction
        ! get i/j indices for 2D COSMO data structure
        i = mind_ilon_rad(ip,jp,ibc)
        j = mind_jlat_rad(ip,jp,ibc)

        IF (itype_aerosol == 1) THEN
          aerdes(i,j) = 0.0_wp
        ELSEIF (itype_aerosol == 2) THEN
          aer_du(i,j) = 0.0_wp
        ENDIF
      ENDDO
      !$acc end parallel
    ENDDO
  ENDIF
#endif

  IF (itype_aerosol == 1) THEN
    DO jp = jsta_comp, jend_comp           ! loop over (coarse) points in j-direction
      !$acc parallel
      !$acc loop gang vector
      DO  ip = ista_comp, iend_comp        ! loop in i-direction
        ! get i/j indices for 2D COSMO data structure
        i = mind_ilon_rad(ip,jp,ibc)
        j = mind_jlat_rad(ip,jp,ibc)

        zdpo           = p0hl(i,j,1)
        zo3h     (ip,jp) = SQRT(hmo3(i,j))**3
        zqcfo    (ip,jp) = zqco2*zdpo
        zaeqso   (ip,jp) = paeops*aersea(i,j)*pvdaes(1)
        zaeqlo   (ip,jp) = paeopl*aerlan(i,j)*pvdael(1)
        zaequo   (ip,jp) = paeopu*aerurb(i,j)*pvdaeu(1)
        zaeqdo   (ip,jp) = paeopd*aerdes(i,j)*pvdaed(1)
        zaetr_top(ip,jp) = 1.0_wp
        zqofo    (ip,jp) = vio3(i,j)*SQRT(zdpo**3)/(SQRT(zdpo**3) + zo3h(ip,jp))
      ENDDO
      !$acc end parallel
    ENDDO
  ELSEIF (itype_aerosol == 2) THEN
    ! new Tegen aerosol climatology: no multiplication with tau(max) as climatology not normalised!
    DO jp = jsta_comp, jend_comp           ! loop over (coarse) points in j-direction
      !$acc parallel
      !$acc loop gang vector
      DO  ip = ista_comp, iend_comp        ! loop in i-direction
        ! get i/j indices for 2D COSMO data structure
        i = mind_ilon_rad(ip,jp,ibc)
        j = mind_jlat_rad(ip,jp,ibc)

        zdpo           = p0hl(i,j,1)
        zo3h     (ip,jp) = SQRT(hmo3(i,j))**3
        zqcfo    (ip,jp) = zqco2*zdpo
        zaeqso   (ip,jp) =  aer_ss(i,j)              *pvdaes(1)
        zaeqlo   (ip,jp) =( aer_or(i,j)+aer_su(i,j) )*pvdael(1)
        zaequo   (ip,jp) =  aer_bc(i,j)              *pvdaeu(1)
        zaeqdo   (ip,jp) =  aer_du(i,j)              *pvdaed(1)
        zaetr_top(ip,jp) =  1.0_wp
        zqofo    (ip,jp) = vio3(i,j)*SQRT(zdpo**3)/(SQRT(zdpo**3) + zo3h(ip,jp))
      ENDDO
      !$acc end parallel
    ENDDO
  ENDIF

  IF (itype_aerosol == 1) THEN
    DO jp = jsta_comp, jend_comp             ! loop over (coarse) points in j-direction
!CDIR UNROLL=10
      !$acc parallel
      !$acc loop seq
      DO k = 1, kend
        !$acc loop gang vector
        DO  ip = ista_comp, iend_comp        ! loop in i-direction
          ! get i/j indices for 2D COSMO data structure
          i = mind_ilon_rad(ip,jp,ibc)
          j = mind_jlat_rad(ip,jp,ibc)

          zdpn             = p0hl(i,j,k+1)
          zaeqsn           = paeops*aersea(i,j)*pvdaes(k+1)
          zaeqln           = paeopl*aerlan(i,j)*pvdael(k+1)
          zaequn           = paeopu*aerurb(i,j)*pvdaeu(k+1)
          zaeqdn           = paeopd*aerdes(i,j)*pvdaed(k+1)
          zaetr_bot        = zaetr_top(ip,jp) * ( MIN (1.0_wp, rti(ip,k,jp)/rti(ip,k+1,jp)) )**ptrpt
          zqcfn            = zqco2 * zdpn
          zqofn            = vio3(i,j)*SQRT(zdpn**3)/(SQRT(zdpn**3) + zo3h(ip,jp))
          rduco2f(ip,k,jp) = zqcfn-zqcfo(ip,jp)
          rduo3f (ip,k,jp) = zqofn-zqofo(ip,jp)
          zaetr            = SQRT(zaetr_bot*zaetr_top(ip,jp))
          raeq1(ip,k,jp)   = (1.0_wp-zaetr) * (ptrbga*dp0(i,j,k)+zaeqln-zaeqlo(ip,jp)+zaeqdn-zaeqdo(ip,jp))
          raeq2(ip,k,jp)   = (1.0_wp-zaetr) * ( zaeqsn-zaeqso(ip,jp) )
          raeq3(ip,k,jp)   = (1.0_wp-zaetr) * ( zaequn-zaequo(ip,jp) )
          raeq4(ip,k,jp)   =         zaetr  *   pvobga*dp0(i,j,k)
          raeq5(ip,k,jp)   =         zaetr  *   pstbga*dp0(i,j,k)

          zqcfo(ip,jp)     = zqcfn
          zqofo(ip,jp)     = zqofn
          zaetr_top(ip,jp) = zaetr_bot
          zaeqso(ip,jp)    = zaeqsn
          zaeqlo(ip,jp)    = zaeqln
          zaequo(ip,jp)    = zaequn
          zaeqdo(ip,jp)    = zaeqdn
        ENDDO
      ENDDO
      !$acc end parallel
    ENDDO
  ELSEIF (itype_aerosol == 2) THEN
    ! new Tegen aerosol climatology: no multiplication with tau(max) as climatology not normalised!
    DO jp = jsta_comp, jend_comp           ! loop over (coarse) points in j-direction
!CDIR UNROLL=10
      !$acc parallel
      !$acc loop seq
      DO  k = 1, kend
       !$acc loop gang vector
        DO  ip = ista_comp, iend_comp        ! loop in i-direction
          ! get i/j indices for 2D COSMO data structure
          i = mind_ilon_rad(ip,jp,ibc)
          j = mind_jlat_rad(ip,jp,ibc)

          zdpn             = p0hl(i,j,k+1)
          zaeqsn           =   aer_ss(i,j)              * pvdaes(k+1)
          zaeqln           =  (aer_or(i,j)+aer_su(i,j)) * pvdael(k+1)
          zaequn           =   aer_bc(i,j)              * pvdaeu(k+1)
          zaeqdn           =   aer_du(i,j)              * pvdaed(k+1)
          zaetr_bot        = zaetr_top(ip,jp) * ( MIN (1.0_wp, rti(ip,k,jp)/rti(ip,k+1,jp)) )**ptrpt
          zqcfn            = zqco2 * zdpn
          zqofn            = vio3(i,j)*SQRT(zdpn**3)/(SQRT(zdpn**3) + zo3h(ip,jp))
          rduco2f(ip,k,jp) = zqcfn-zqcfo(ip,jp)
          rduo3f (ip,k,jp) = zqofn-zqofo(ip,jp)
          zaetr            = SQRT(zaetr_bot*zaetr_top(ip,jp))

          raeq1(ip,k,jp)   = (1.0_wp-zaetr)*( ptrbga*dp0(i,j,k) + zaeqln - zaeqlo(ip,jp) )
          raeq2(ip,k,jp)   = (1.0_wp-zaetr)*(zaeqsn-zaeqso(ip,jp))
          raeq3(ip,k,jp)   = (1.0_wp-zaetr)*(zaeqdn-zaeqdo(ip,jp))
          raeq4(ip,k,jp)   = (1.0_wp-zaetr)*(zaequn-zaequo(ip,jp))
          raeq5(ip,k,jp)   =    zaetr * pstbga*dp0(i,j,k)

          zqcfo(ip,jp)     = zqcfn
          zqofo(ip,jp)     = zqofn
          zaetr_top(ip,jp) = zaetr_bot
          zaeqso(ip,jp)    = zaeqsn
          zaeqlo(ip,jp)    = zaeqln
          zaequo(ip,jp)    = zaequn
          zaeqdo(ip,jp)    = zaeqdn
        ENDDO
      ENDDO
      !$acc end parallel
    ENDDO
  ENDIF

  !$acc end data

END SUBROUTINE co2_scenarios_and_absorbers

!==============================================================================
!==============================================================================

SUBROUTINE calc_rad_corrections (thetain, phi, horizonte, smu0, rlati,       &
                                 rloni, deksini, dekcosi, zeit0i, fcor,      &
                                 phi_sun, theta_sun, theta,                  &
                                 idim, jdim, nhordim, ibc,                   &
                                 ista_comp, iend_comp, jsta_comp, jend_comp, &
                                 idbg, ierror, yerror)

!------------------------------------------------------------------------------
!
! Description:
!   Compute parameters swdir_cor needed for the grid scale topographic 
!   correction of direct solar radiation at the surface.
!   If correction option is chosen, this subroutine is called before fesft.
!   Following Mueller and Sherrer (2005), MWR
!
!------------------------------------------------------------------------------

! Variables of the parameter list
! -------------------------------

INTEGER                 ,  INTENT (IN)  ::   &
  idim, jdim, nhordim,           & ! field dimensions
  ista_comp, iend_comp,          & ! loop start / end indices
  jsta_comp, jend_comp,          & !
  ibc,                           & ! block to compute
  idbg                             ! for debug output

REAL    (KIND=wp)       ,  INTENT (IN)     ::   &
  thetain  (:,:),                & ! slope angle
  phi      (:,:),                & ! slope aspect
! horizonte(idim,jdim,nhordim),  & ! horizont
  horizonte(:,:,:),              & ! horizont
  rlati    (:,:),                & ! latitude (geogr) [rad]
  rloni    (:,:),                & ! longitude(geogr) [rad]
  smu0     (:,:)                   ! sun zenith angle

REAL    (KIND=wp)       ,  INTENT (IN)     ::   &
  deksini,                       & ! sin of sun declination angle
  dekcosi,                       & ! cos of sun declination angle
  zeit0i                           !T.R.

! Output of the routine
REAL    (KIND=wp)       ,  INTENT (OUT) ::   &
  fcor     (:,:)

INTEGER,                   INTENT (OUT) ::   &
  ierror

CHARACTER(LEN=*),          INTENT (OUT) ::   &
  yerror

! Local arrays, passed through argument list (declared as INTENT(OUT))
! --------------------------------------------------------------------

REAL    (KIND=wp),         INTENT (OUT) ::   &
  phi_sun  (:,:),                & ! sun azimuth angle [rad]
  theta_sun(:,:),                & ! sun elevation angle [rad]
  theta    (:,:)

! Local parameters and variables
! ------------------------------

REAL (KIND=wp),     PARAMETER :: &
  zepemu = 1.0E-07_wp,              &
  rtod = 57.2957795_wp                ! radiantas to degrees

REAL    (KIND=wp)       ::       &
  zeitrad,                       & ! T.R.
  phi_s,                         & !
  x1,x2,ha_sun                     !


INTEGER                  ::        &
  ii,shadow,i,k,j, ip, jp

LOGICAL                  ::        &
  lshade, lslope_aspect            !switches

!------------------------------------------------------------------------------

ierror        = 0
yerror(:)     = ' '

lshade        = .TRUE.
lslope_aspect = .TRUE.

!$acc data &
!---- Argument arrays - intent(in)
!$acc present ( thetain, phi, horizonte, smu0, rlati, rloni  ) &
!---- Argument arrays - intent(out)
!$acc present ( fcor                                         ) &
!$acc present ( phi_sun,theta_sun,theta                      )

!$acc parallel 
!$acc loop gang vector collapse(2)
DO jp = jsta_comp, jend_comp
  DO ip = ista_comp, iend_comp

    ! get i/j indices for 2D COSMO data structure
    i = mind_ilon_rad(ip,jp,ibc)
    j = mind_jlat_rad(ip,jp,ibc)

    ! sun elevation angle
    theta_sun(ip,jp) = ASIN(smu0(ip,jp))

    ! sun azimuth angle
    zeitrad = zeit0i + rloni(i,j)  !T.R.
    x1 = dekcosi * SIN(zeitrad) / COS(theta_sun(ip,jp))
    x2 = ( SIN(rlati(i,j)) * dekcosi * COS(zeitrad) - &
                        COS(rlati(i,j)) * deksini ) / COS(theta_sun(ip,jp))

    IF (x2 < -1.0_wp) x2 = -1.0_wp
    IF (x2 >  1.0_wp) x2 =  1.0_wp
    phi_s = COS(x2)
    IF (x1 < 0) phi_s = - phi_s
    phi_sun(ip,jp) = phi_s + pi

    ! sun elevation angle corrected by refraction (empiric formula, A.H.Siemer(1988))
    theta_sun(ip,jp) = theta_sun(ip,jp) + (1.569000_wp - theta_sun(ip,jp)) / &
                                    (185.5_wp + 3620.0_wp * theta_sun(ip,jp))

    ! night or day
    IF (theta_sun(ip,jp) < 0.0_wp) THEN
       theta_sun(ip,jp) = 0.0_wp
    ENDIF

    ! compute shadow
    ! the horizon has a spatial resolution of 360/nhordim degrees. a distance weighted
    ! linear interpolation is done between the two neighboring points.
    ii = INT(rtod*phi_sun(ip,jp)/(360.0_wp/nhordim))
    IF (ii >= nhordim) THEN
      ii = nhordim - 1
    ENDIF

    k  = MOD(ii+1,24)

    IF (horizonte(i,j,k+1) > 90.0_wp .OR. horizonte(i,j,k+1) < 0.0_wp) THEN
#ifndef _OPENACC
      PRINT *,'!!ERROR!!, horizon_angle > 90deg or < 0deg',horizonte(i,j,k+1)
      !US there shall be no stop in a parallel program!!!  STOP
      yerror = '*** ERROR:  Error with horizon_angle in SR calc_rad_corrections'
#endif
      !$acc atomic write      !XL:This may not work with PGI
      ierror = 2777
      !$acc end atomic
    ENDIF

    ha_sun = (horizonte(i,j,k+1)*(rtod*phi_sun(ip,jp)-15*ii)+ &
              horizonte(i,j,ii+1)*(15*(ii+1)-rtod*phi_sun(ip,jp)))/15.0_wp

    ! compute shadowmask
    shadow = 1
    IF (rtod*theta_sun (ip,jp) < ha_sun .AND. lshade) THEN
      shadow = 0
    ENDIF

    ! compute fcor
    ! slope angle and aspect switched off
    IF (.NOT. lslope_aspect) THEN
      theta(ip,jp) = 0.0_wp
    ELSE
      theta(ip,jp) = thetain(i,j)
    ENDIF

    IF (theta_sun(ip,jp) > 0.01_wp) THEN
      ! Mueller and Scherrer (2005) formula (MWR)
      ! fcor(i,j) = shadow * (1 + ( tan(theta(ip,jp)) / tan(theta_sun(ip,jp)) )*&
      !              cos( phi_sun(ip,jp) - phi(i,j)) )

      ! New formula (lower correction, theoretically correct derived)
      fcor(ip,jp) = REAL(shadow, wp) * ( COS(theta(ip,jp)) + (SIN(theta(ip,jp))/TAN(theta_sun(ip,jp)) )*&
                    COS( phi_sun(ip,jp) - phi(i,j)) )
    ELSE
      fcor(ip,jp) = 1.0_wp
    ENDIF

    ! Consistency check to avoid negative corrections:
    ! active in situations with low sun elevation (slope angles > sun elevation)
    ! when the slope aspect is greater than the daily maxima or smaller than the
    ! daily minima of the sun azimuth angle (during the sunshine time, a kind 
    ! of self shading effect).
    IF (fcor(ip,jp) < 0.0_wp) THEN
      fcor(ip,jp) = 0.0_wp
    ENDIF

#ifndef _OPENACC
    IF ( (idbg > 20) .AND. (i == 4) ) THEN
      PRINT *, '   calc_rad_corrections:  debug point:  ', i, j
      PRINT *, '   deksini, dekcosi, zeitrad = ', deksini, dekcosi, zeitrad
      PRINT *, '   rlat, rlon, theta, phi, smu0 = ', rlati(i,j), rloni(i,j), &
                                                      theta(ip,jp), phi(i,j), smu0(i,j)
      PRINT *, '   ii, k, horizon = ', ii, k, horizonte(i,j,ii+1), horizonte(i,j,k+1)
      PRINT *, '   ha_sun, theta_sun, phi_sun = ', ha_sun, rtod*theta_sun(ip,jp), &
                                                           rtod*phi_sun(ip,jp)

      IF ( (shadow == 0) .AND. (theta_sun(ip,jp) > 0.01_wp) ) THEN
        PRINT *, '    calc_rad_corrections:  ', 'DAY-SHADOW'
      ENDIF
      IF ( (shadow == 1) .AND. (theta_sun(ip,jp) > 0.01_wp) ) THEN
        PRINT *, '    calc_rad_corrections:  ', 'DAY-SUN'
      ENDIF
      IF (theta_sun(ip,jp) <= 0.01_wp ) then
        PRINT *, '    calc_rad_corrections:  ', 'NIGHT'
      ENDIF
      PRINT *, '    calc_rad_corrections: fcor  ', fcor(ip,jp)
    ENDIF
#endif

  ENDDO
ENDDO
!$acc end parallel

#ifdef _OPENACC
IF (ierror>0)  yerror = '*** ERROR:  Error in calc_rad_corrections'
#endif

!$acc end data


END SUBROUTINE calc_rad_corrections

!==============================================================================

END MODULE radiation_utilities
