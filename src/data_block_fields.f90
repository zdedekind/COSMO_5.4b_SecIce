!+ Data module for all fields in block data structure
!------------------------------------------------------------------------------

MODULE data_block_fields

!------------------------------------------------------------------------------
!
! Description:
!
! This module declares all fields which are used in the physics and have
! the block data structure field(nproma,ke)
!
! All fields are declared as allocatable arrays
!
! Convention: 
!  - name of fields are the same as in data_fields with a "_b", for block, at the end
!  - if timelevel nnew is used, the suffix is "_new_b"
!
! Current Code Owner: MeteoSwiss, Xavier Lapillonne
!  phone: +41 58 460 9237
!  fax  : +41 58 460 9278
!  email:  xavier.lapillonne@meteoswiss.ch
!
! History:
! Version      Date       Name
! ----------   ---------- ----
! V5_1         2014-11-28 Xavier Lapillonne
!  Initial release
! V5_3         2015-10-09 Ulrich Schaettler, Xavier Lapillonne
!  Adaptations to MCH (POMPA) version
! V5_4         2016-03-10 Oliver Fuhrer
!  Updated code owner information
! V5_4a        2016-05-10 Ulrich Schaettler
!  Included new fields for TKE scheme
! V5_4b        2016-07-12 Ulrich Blahak, Xavier Lapillonne, Ulrich Schaettler
!  Included new field h0noise_b for lartif_data (UB)
!  Added new fields for convection schemes (XL)
!  Added new fields for COSMO-ICON versions of surface schemes (US)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Declarations:
!
! Modules used:
USE kind_parameters, ONLY :  wp    

!==============================================================================

IMPLICIT NONE

!==============================================================================

! 1. constant fields for the reference atmosphere                        (unit)
! -----------------------------------------------

  REAL (KIND=wp), ALLOCATABLE, TARGET     ::   &
    p0_b   (:,:),      & ! reference pressure at full levels             ( Pa  )
    p0hl_b (:,:),      & ! reference pressure at half levels             ( Pa  )
    rho0_b (:,:),      & ! reference density at the full model levels    (kg/m3)
    dp0_b  (:,:),      & ! reference pressure thickness of layers        ( Pa  )
    hhl_b  (:,:),      & ! geometrical height of half levels             ( m   )
    dz_b   (:,:)         ! vertical difference of geometrical height of half levels  ( m   )

! 2. external parameter fields                                           (unit)
! -----------------------------                           

  REAL (KIND=wp), ALLOCATABLE, TARGET     ::   &
    hsurf_b      (:),  & ! height of surface topography                  (  m  )
    fr_land_b    (:),  & ! fraction of land in a grid element              --
    plcov_b      (:),  & ! fraction of plant cover                         --
    rootdp_b     (:),  & ! depth of the roots                            (  m  )
    sso_stdh_b   (:),  & ! standard deviation of sub-grid scale orography(  m  )
    sai_b        (:),  & ! surface area index                              --
    tai_b        (:),  & ! transpiration area index                        --
    eai_b        (:),  & ! (evaporative) earth area index                  --
    l_pat_b      (:),  & ! effective length scale of circulation patterns  --
    rsmin2d_b    (:),  & ! minimum stomata resistance                    ( s/m )
    fr_lake_b    (:),  & ! lake fraction in a grid element               ([0,1])
    depth_lk_b   (:),  & ! lake depth                                    (  m  )
    fc_b         (:)     ! coriolis-parameter                            ( 1/s )

  LOGICAL, ALLOCATABLE                    ::   &
    llandmask_b  (:),  & ! land-sea-mask
    llakemask_b  (:)     ! mask for fresh-water lakes

  INTEGER, ALLOCATABLE, TARGET            ::   &
    isoiltyp_b   (:)     ! type of the soil (keys 0-9)                     --

! 3. prognostic variables                                               (unit)
! -----------------------
  REAL (KIND=wp), ALLOCATABLE, TARGET     ::   &
    ptot_b     (:,:),  & ! total pressure p0+pp
    pp_b       (:,:),  & ! deviation from the reference pressure         ( pa  )
    t_b        (:,:),  & ! temperature                                   (  k  )
    u_m_b      (:,:),  & ! zonal wind speed at mass point                ( m/s )
    v_m_b      (:,:),  & ! meridional wind speed at mass point           ( m/s )
    w_b        (:,:),  & ! vertical wind speed (defined on half levels)  ( m/s )
    rho_b      (:,:),  & ! total density of air                          (kg/m3)
    qv_b       (:,:),  & ! specific water vapor content                  (kg/kg) 
    qc_b       (:,:),  & ! specific cloud liquid water content           (kg/kg)
    qi_b       (:,:),  & ! specific cloud ice content                    (kg/kg)
    qr_b       (:,:),  & ! specific rain content                         (kg/kg)
    qs_b       (:,:),  & ! specific snow content                         (kg/kg)
    qg_b       (:,:),  & ! specific graupel content                      (kg/kg)

! fields of the turbulent scheme defined on half-levels:
    tke_b    (:,:,:),  & ! SQRT(2*TKE); TKE='turbul. kin. energy'        ( m/s )          

! time level new
    t_new_b    (:,:),  & ! temperature                                   (  k  )
    pp_new_b   (:,:),  & ! deviation from the reference pressure         ( pa  )
    qv_new_b   (:,:),  & ! specific water vapor content                  (kg/kg) 
    qc_new_b   (:,:),  & ! specific cloud liquid water content           (kg/kg)
    qi_new_b   (:,:),  & ! specific cloud ice content                    (kg/kg)
    qr_new_b   (:,:),  & ! specific rain content                         (kg/kg)
    qs_new_b   (:,:),  & ! specific snow content                         (kg/kg)
    qg_new_b   (:,:)     ! specific graupel content                      (kg/kg)

! 4. tendency fields for the prognostic variables                     (unit )
! -----------------------------------------------
  REAL (KIND=wp), ALLOCATABLE, TARGET     ::   &
    utens_b    (:,:),  & ! u-tendency without sound-wave terms           ( m/s2)
    vtens_b    (:,:),  & ! v-tendency without sound-wave terms           ( m/s2)
    wtens_b    (:,:),  & ! w-tendency without sound-wave terms           ( m/s2)
                         ! (defined on half levels)
    ttens_b    (:,:),  & ! t-tendency without sound-wave terms           ( K/s )
    pptens_b   (:,:),  & ! pp-tendency without sound-wave terms          (Pa/s )
    tketens_b  (:,:),  & ! non-advective tke-tendency (on half-levels)   ( m/s )
    qvtens_b   (:,:),  & ! qv-tendency                                   (     )
    qctens_b   (:,:),  & ! qc-tendency                                   (     )
    qitens_b   (:,:)     ! qi-tendency                                   (     )


! 5. fields for surface values                                           (unit )
! ----------------------------
  REAL (KIND=wp), ALLOCATABLE, TARGET     ::   &
    ps_b       (:),    & ! surface pressure                      ( pa  )
    t_s_b      (:),    & ! surface temperature                           (  K  )
    t_s_new_b  (:),    & ! surface temperature at time level nnew        (  K  )
    t_g_b      (:),    & ! weighted surface temperature                  (  K  )
    t_g_new_b  (:),    & ! weighted surface temperature at nnew          (  K  )
    qv_s_b     (:),    & ! specific water vapor content at the surface   (kg/kg)
    qv_s_new_b (:)       ! specific water vapor content at the surface   (kg/kg)

! 6. fields that are computed in the parametrization and dynamics        (unit )
! ---------------------------------------------------------------

! from microphysics
! -----------------
  REAL (KIND=wp), ALLOCATABLE, TARGET     ::   &
    tinc_lh_b  (:,:),  & ! temperature increment due to latent heat      (  K  )
    qrs_b      (:,:),  & ! precipitation water content (water loading)   (kg/kg)
    prr_gsp_b  (:)  ,  & ! precipitation rate of rain, grid-scale        (kg/m2*s)
    prs_gsp_b  (:)  ,  & ! precipitation rate of snow, grid-scale        (kg/m2*s)
    prg_gsp_b  (:)  ,  & ! precipitation rate of graupel, grid-scale     (kg/m2*s)

    !latent heat nudging
    tt_lheat_b (:,:),  & ! profile of t-increments due to latent heating ( K/s )
                         ! (stored for current and previous timestep)
    tt_lheat_new_b(:,:) ,& ! profile of t-increments due to latent heating  ( K/s )
                           ! (stored for current and previous timestep)
    qrsflux_b(:,:)         ! total precipitation flux
     
! from radiation scheme (needed in surface TERRA)                        (unit )
! -----------------------------------------------
  REAL (KIND=wp), ALLOCATABLE, TARGET     ::   &
    sobs_b       (:),  & ! solar radiation at the ground                 ( W/m2)
    thbs_b       (:),  & ! thermal radiation at the ground               ( W/m2)
    pabs_b       (:)     ! photosynthetic active radiation at the ground ( W/m2)


! from the sub-grid scale orography
! ---------------------------------
  REAL (KIND=wp), ALLOCATABLE, TARGET     ::   &
    ut_sso_b   (:,:),  & ! u-tendency due to SSO                         (m/s2)
    vt_sso_b   (:,:),  & ! v-tendency due to SSO                         (m/s2)
    tket_sso_b (:,:)     ! TKE-tendency due to SSO wake production       (m2/s3)


! from turbulence scheme                                                 (unit )
! ----------------------
  REAL (KIND=wp), ALLOCATABLE, TARGET     ::   &
    l_hori_b     (:),  & ! horizontal grid spacing (location dependent horizontal scale)
    gz0_b        (:),  & ! surface roughness * g                         (m2/s2)

!   turbulent transfer coefficients at the surface
    tcm_b        (:),  & ! for momentum                                  ( -- )
    tch_b        (:),  & ! for heat and moisture                         ( -- )
!   turbulent transfer factors for laminar- and roughness-layer transfer
    tfm_b        (:),  & ! of momentum                                     --
    tfh_b        (:),  & ! of scalars                                      --
    tfv_b        (:),  & ! of water vapor compared to heat                 --
!   turbulent transfer velocities at the surface
    tvm_b        (:),  & ! for momentum                                  ( m/s)
    tvh_b        (:),  & ! for heat and moisture                         ( m/s)
    tkr_b        (:),  & ! reference surf. diff. coeff. (l*Ustar)        (m2/s)
    tkred_sfc_b  (:),  & ! reduction factor for minimum diffusion coefficients near the surface
!   turbulence statistics in the atmosphere
    rcld_b     (:,:),  & ! standard deviation of the saturation deficit    --
!   coefficients for turbulent diffusion in the atmosphere
!   (defined on half levels)
                         ! vertical   turbulent diffusion coefficients
    tkvm_b     (:,:),  & ! ... for momentum                              (m2/s)
    tkvh_b     (:,:),  & ! ... for heat and moisture                     (m2/s)
                         ! horizontal turbulent diffusion coefficients
    tkhm_b     (:,:),  & ! ... for momentum                              (m2/s)
    tkhh_b     (:,:),  & ! ... for heat and moisture                     (m2/s)

    ! variables needed for turb_prepare and other turbulence issues
    hdef2_b    (:,:),  & ! horizontal deformation square
    hdiv_b     (:,:),  & ! horizontal divergence
    dwdx_b     (:,:),  & ! horizontal gradients of vertical wind
    dwdy_b     (:,:),  & ! horizontal gradients of vertical wind

    edr_b      (:,:),  & ! eddy dissipation rate of TKE (EDR)            (m2/s3)
    tket_hshr_b(:,:),  & ! TKE-tendency due to (sep.) horiz. shear       (m2/s3)
    tket_adv_b (:,:),  & ! pure advective tke-tendency (on half-levels)  (m/s )
    dqvdt_b    (:,:),  & ! threedimensional moisture convergence         (1/s )
    ut_turb_b  (:,:),  & ! u-tendency due to turbulent vert. diffusion   (m/s2)
    vt_turb_b  (:,:)     ! v-tendency due to turbulent vert. diffusion   (m/s2)


! from the surface schemes                                                 (unit )
! ------------------------
  REAL (KIND=wp), ALLOCATABLE, TARGET     ::   &
    t_snow_b           (:), & ! snow temperature                           (  K  )
    t_snow_new_b       (:), & !
    w_snow_b           (:), & ! water content of snow                      (m H2O)
    w_snow_new_b       (:), & !
    h_snow_b           (:), & ! snow height                                (  m  )
    h_snow_new_b       (:), & !
    rho_snow_b         (:), & ! prognostic snow density                    (kg/m3)
    rho_snow_new_b     (:), & !
    freshsnow_b        (:), & ! weighting function indicating 'freshness' of snow
    fr_snow_b          (:), & ! surface fraction covered by snow           (  -  )
    w_i_b              (:), & ! water content of interception water        (m H2O)
    w_i_new_b          (:), & !
    w_p_b              (:), & ! water content of interception water pond   (m H2O)
    w_p_new_b          (:), & !
    w_s_b              (:), & ! water content of interception snow water   (m H2O)
    w_s_new_b          (:), & !
    t_so_b           (:,:), & ! multi-layer soil temperature               (  K  )
    t_so_new_b       (:,:), & !
    w_so_b           (:,:), & ! multi-layer soil moisture                  (m H2O)
    w_so_new_b       (:,:), & !
    w_so_ice_b       (:,:), & ! multi-layer soil ice                       (m H2O)
    w_so_ice_new_b   (:,:), & !

    t_snow_mult_b    (:,:), & ! temperature of the snow-surface            (  K  )
    t_snow_mult_new_b(:,:), & !
    w_snow_mult_b    (:,:), & ! total (liquid+solid) water content of snow (m H2O)
    w_snow_mult_new_b(:,:), & !
    wliq_snow_b      (:,:), & ! liquid water content in snow               (m H2O)
    wliq_snow_new_b  (:,:), & !
    rho_snow_mult_b  (:,:), & ! prognostic snow density                    (kg/m3)
    rho_snow_mult_new_b(:,:),&!
    dzh_snow_mult_b  (:,:), & ! layer thickness between half levels in snow(  m  )
    dzh_snow_mult_new_b(:,:),&!

    runoff_s_b         (:), & ! surface water runoff; sum over forecast
    runoff_g_b         (:), & ! soil water runoff; sum over forecast
    lhfl_bs_b          (:), & ! latent heat flux from bare soil evaporation
    lhfl_pl_b        (:,:), & ! average latent heat flux from plants
    rstom_b            (:)    ! stomata resistance                   

! fields for prognostic variables of the lake model FLake or ocean variables
  REAL (KIND=wp), ALLOCATABLE, TARGET     ::   &
    qmomflux_b       (:), & ! momentum flux at the surface                 ( N/m2)
    fetch_lk_b       (:), & ! wind fetch over lake                         (  m  )
    dp_bs_lk_b       (:), & ! depth of the thermally active layer
                            ! of bottom sediments                          (  m  )
    t_bs_lk_b        (:), & ! climatological temperature at the bottom of
                            ! the thermally active layer of sediments      (  K  )
    gamso_lk_b       (:), & ! attenuation coefficient for
                            ! solar radiation in lake water                ( 1/m )
    t_mnw_lk_b       (:), & ! mean temperature of the water column         (  K  )
    t_mnw_lk_new_b   (:), & !
    t_wml_lk_b       (:), & ! mixed-layer temperature                      (  K  )
    t_wml_lk_new_b   (:), & !
    t_bot_lk_b       (:), & ! temperature at the water-bottom sediment
                            ! interface                                    (  K  )
    t_bot_lk_new_b   (:), & !
    t_b1_lk_b        (:), & ! temperature at the bottom of the upper layer
                            ! of the sediments                             (  K  )
    t_b1_lk_new_b    (:), & !
    c_t_lk_b         (:), & ! shape factor with respect to the
                            ! temperature profile in lake thermocline      (  -  )
    c_t_lk_new_b     (:), & !
    h_ml_lk_b        (:), & ! thickness of the mixed-layer                 (  m  )
    h_ml_lk_new_b    (:), & !
    h_b1_lk_b        (:), & ! thickness of the upper layer
                            ! of bottom sediments                          (  m  )
    h_b1_lk_new_b    (:), & !

    t_ice_b          (:), & ! temperature of ice/water surface             (  K  )
    t_ice_new_b      (:), & !
    h_ice_b          (:), & ! lake/sea ice thickness                       (  m  )
    h_ice_new_b      (:)    !

! from the convection schemes (Tiedtke, Tiedtke-Bechthold, shallow)         (unit ) 
! ---------------------------
  REAL (KIND=wp), ALLOCATABLE, TARGET     ::   &
    clc_con_b      (:,:), & !
    clw_con_b      (:,:), & !
    prr_con_b        (:), & ! precipitation rate of rain, convective        (kg/m2*s)
    prs_con_b        (:), & ! precipitation rate of snow, convective        (kg/m2*s)
    prne_con_b       (:), & ! precipitation rate no evaporation, convective (kg/m2*s)
    bas_con_b        (:), & !
    top_con_b        (:), & !
    tt_conv_b      (:,:), & !
    ttdiab_conv_b  (:,:), & ! pure diabatic temperature tendency due to convection ( K/s  )
    ttens_conv_b   (:,:), & ! temperature tendency for the Tiedtke-Bechtold convection ( K/s  )
    qvt_conv_b     (:,:), & !
    qct_conv_b     (:,:), & !
    qit_conv_b     (:,:), & !
     ut_conv_b     (:,:), & !
     vt_conv_b     (:,:), & !
   tket_conv_b     (:,:), & ! TKE-tendency due to convective buoyancy       (m2/s3)
    mflx_con_b       (:), & ! convective mass flux
    cape_con_b       (:), & !
    qcvg_con_b       (:), & !
     tke_con_b       (:), & ! turbulent kinetic energy convective !MR: not defined anywhere
   vgust_con_b       (:)    !

  ! Additional fields for the convection interface
  REAL (KIND=wp), ALLOCATABLE, TARGET     ::   &
    fif_b          (:,:), & !
    fih_b          (:,:), & !
    pf_b           (:,:), & !
    ph_b           (:,:), & !
    w_conv_b       (:,:), & ! (averaged) w for convection
    dqvdt_conv_b   (:,:), & ! (averaged) dqvdt for convection
    qhfl_b           (:), & ! (averaged) qhfl (qvsflx) for convection
    ! for Tiedtke-Bechtold
    cloud_num_b    (:,:), & !
    trop_mask_b    (:,:), & !
    mtn_mask_b     (:,:)    !

  INTEGER, ALLOCATABLE, TARGET ::  &
    conv_k850_b    (:,:), & !
    conv_k950_b    (:,:)    !


!   fields that are computed in the dynamics and / or physics
  REAL (KIND=wp), ALLOCATABLE, TARGET     ::   &
    shfl_s_b         (:), & ! sensible heat flux (surface)                  ( W/m2)
    lhfl_s_b         (:), & ! latent heat flux (surface)                    ( W/m2)
    qvsflx_b         (:)    ! surface flux of water vapour                  (kg/m2s)

! 7. fields for model output and diagnostics                          (unit )
! ------------------------------------------
  REAL (KIND=wp), ALLOCATABLE, TARGET     ::   &
    t_2m_b           (:), & ! temperature in 2m                             (  K  )
    qv_2m_b          (:), & ! specific water vapor content in 2m            (kg/kg)
    td_2m_b          (:), & ! dew-point in 2m                               (  K  )
    rh_2m_b          (:), & ! relative humidity in 2m                       (  %  )
    u_10m_b          (:), & ! zonal wind in 10m                             ( m/s )
    v_10m_b          (:)    ! meridional wind in 10m                        ( m/s )


! Miscellaneous fields used:
!---------------------------
  REAL (KIND=wp), ALLOCATABLE, TARGET     ::   &
     pertstoph_b(:,:), & ! stochastic multiplier of physics tendencies
     h0noise_b  (:)      ! random noise for idealized surface fluxes (lartifdata=.true.)  

! Additional pointer for the microphysics (needed, if microphysics is called 
! at the beginning of the time loop
!---------------------------------
  REAL (KIND=wp), POINTER     :: &
    t_nx_b     (:,:),  & ! temperature                                  (  k  )
    pp_nx_b    (:,:),  & ! deviation from the reference pressure        ( pa  )
    qv_nx_b    (:,:),  & ! specific water vapor content                 (kg/kg) 
    qc_nx_b    (:,:),  & ! specific cloud liquid water content          (kg/kg)
    qi_nx_b    (:,:),  & ! specific cloud ice content                   (kg/kg)
    qr_nx_b    (:,:),  & ! specific rain content                        (kg/kg)
    qs_nx_b    (:,:),  & ! specific snow content                        (kg/kg)
    qg_nx_b    (:,:),  & ! specific graupel content                     (kg/kg)
    tt_lheat_nx_b(:,:)

END MODULE data_block_fields
