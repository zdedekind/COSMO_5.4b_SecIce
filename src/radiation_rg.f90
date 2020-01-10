!+ Source module with all subroutines for Ritter-Geleyn radiation scheme
!------------------------------------------------------------------------------

MODULE radiation_rg

!------------------------------------------------------------------------------
!
! Description:
!
!   The program package for the Ritter-Geleyn parameterization of radiative 
!   transfer consists of following routines:
!         fesft, coe_so, coe_th, inv_so, inv_th, opt_so and opt_th.
!   In addition, two subroutines are provided to allocate / deallocate 
!   necessary working arrays:  radiation_rg_wkarr_alloc, radiation_rg_wkarr_dealloc
!
!   All parametric data that are required by these routines are defined in the
!   data module radiation_data.
!
!   Additionally, the routine radiation_init has to be called once before the 
!   first call of the driving routine fesft for the radiation package in the 
!   driving routine radiation_organize. Aerdis is a small help routine to 
!   receive some parameters for the vertical distribution of background 
!   aerosol in the driving routine.
!
! Current Code Owner: DWD, Bodo Ritter
!  phone:  +49  69  8062 2703
!  fax:    +49  69  8062 3721
!  email:  Bodo.Ritter@dwd.de
!
! History:
! Version      Date       Name
! ----------   ---------- ----
! V5_3         2015-10-09 Stefan Ruedisuehli, Ulrich Schaettler
!  Initial release (based on src_radiation from V5_2)
! V5_4         2016-03-10 Xavier Lapillonne, Heike Vogel (KIT)
!  Some fixes for running radiation on GPUs (XL)
!  Adapted COSMO-ART code to blocked data structure (HV)
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

USE kind_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    sp,        & ! KIND-type parameter for real variables (single precision)
    dp           ! KIND-type parameter for real variables (double precision)

!------------------------------------------------------------------------------

USE radiation_data, ONLY : &
  coali,                   & !
  cobti,                   & !
  jpsol   , & ! Number of solar spectral intervals
  jpther  , & ! Number of thermal spectral intervals
  jpspec  , & ! =jpsol+jpther (Total number of spectral intervals)

  solant  , & ! Fraction of solar energy at TOA contained in individual
              ! solar spectral intervals

  planck  , & ! coefficients for the description of the fraction of the total
              ! black body radiation contained in an thermal spectral interval
              ! as a function (2.order polynomial) of temperature
  zketypa , & ! (r) following ROBERTS ET AL. 1976,(a) implicitly derived 
              ! from the AFGL spectral data
  ztetypa , & ! (r) following ROBERTS ET AL. 1976,(a) implicitly derived 
              ! from the AFGL spectral data
  zteref  , & ! reference temperature

  ! absorption properties of atmospheric gases
  ncgas   , & ! number of coefficients for each interval and gas (maximum=7)
  nfast   , & ! control variable for choice between ESFT/FESFT method in 
              ! each interval
  coai    , & ! weigthing coefficients
  cobi    , & ! absorption coefficients

! array data variables for opt_th and opt_so with intent (in):
  zaea, zaes, zaeg, zaef,         & !
  zlwe, zlwemn, zlwemx,           & !
  zlww, zlwg,                     & !
  ziwe, ziwemn, ziwemx,           & !
  ziww, ziwg,                     & !
  zrsc                              !

!------------------------------------------------------------------------------

#ifdef COSMOART
USE data_runcontrol,    ONLY : l_cosmo_art

USE data_cosmo_art,        ONLY: &
    lgas         ,               & ! with gas phase chemistry
    laero        ,               & ! with aerosol
    lrad_dust    ,               & ! mineral dust aerosols
    lrad_seas    ,               & ! sea salt aerosols
    lrad_aero                      ! anthropogenic aerosols
USE data_block_fields_art, ONLY: &
    Eup_b        ,               &
    Edown_b      ,               &
    Edir_b

#ifdef TWOMOM_SB
USE src_cloud_opt_reff, ONLY :   &
   reff_for_rad,                 &
   calc_cloud_opt,               &
   calc_optthick_cloud,          &
   zlwoda_so_prefac,             &
   zlwods_so_prefac,             &
   zlwb0_so,                     &
   zlwb_so_prefac,               &
   ziwoda_so_prefac,             &
   ziwods_so_prefac,             &
   ziwb0_so,                     &
   ziwb_so_prefac,               &
   zlwoda_th_prefac,             &
   zlwods_th_prefac,             &
   zlwb0_th,                     &
   zlwb_th,                      &
   ziwoda_th_prefac,             &
   ziwods_th_prefac,             &
   ziwb0_th
#endif
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================
! Working arrays - fesft_dp
!==============================================================================

  REAL    (KIND=dp), ALLOCATABLE        ::  &
    zketyp (:) ,      & !
    ztetyp (:)          !

! Local (automatic) arrays:
! ------------------------
! Arrays local to *fesft* or required for communication with
! subroutines called from *fesft*

  REAL    (KIND=dp), ALLOCATABLE        ::  &
    ! 'Grey' and gaseous fluxes for individual spectral intervals
    ! "_c" means: corrected if lradtopo
    zflux    (:,:)  , & ! (W/m**2)
    zflux_c  (:,:)  , & ! (W/m**2)
    zfluxi   (:,:)  , & ! 1./(W/m**2)
    zfluxu   (:,:)  , & ! (W/m**2)
    zfluxu_c (:,:)  , & ! (W/m**2)
    zfluxui  (:,:)  , & ! 1./(W/m**2)
    zfluxd   (:,:)  , & ! (W/m**2)
    zfluxd_c (:,:)  , & ! (W/m**2)
    zfluxdi  (:,:)  , & ! 1./(W/m**2)
    zfgas    (:,:)  , & ! (W/m**2)
    zfgasu   (:,:)  , & ! (W/m**2)
    zfgasd   (:,:)      ! (W/m**2)

  REAL    (KIND=dp), ALLOCATABLE        ::  &
    pbbr     (:,:)  , & ! (W/m**2) Black body radiation at layer boundaries
    pflpt    (:)    , & ! Solar flux at TOA
    palp     (:)    , & ! Solar surface albedo for parallel radiation
    pqsmu0   (:)    , & ! Inverse of cosine of zenith angle
    palogt   (:,:)  , & ! ln T
    palogp   (:,:)  , & ! ln p
    papra    (:)    , & ! (Pa) pressure at one level
    pduh2oc  (:,:)  , & ! layer water vapour content (Pa), cloudy
    pduh2of  (:,:)  , & ! layer water vapour content (Pa), cloud-free)
    pdulwc   (:,:)  , & ! (Pa H2O-liquid) layer incloud liquid water content
    pduiwc   (:,:)  , & ! (Pa H2O-ice) layer incloud ice content
    prholwc  (:,:)  , & ! (kg/m**3)
    prhoiwc  (:,:)  , & ! (kg/m**3)
    zduetpc  (:,:)  , & ! water vapour e-type contribution (cloudy)
    zduetpf  (:,:)  , & ! water vapour e-type contribution (cloud-free)
    zapre    (:)    , & ! (Pa) pressure at one level

    ! layer mean temperature, water vapour mixing ratio, utility arrays
    ztm      (:)    , & !
    zzwv     (:)    , & !
    zcpo     (:)    , & !
    zcpn     (:)    , & !
    zcmo     (:)    , & !
    zcmn     (:)        !

  REAL    (KIND=dp), ALLOCATABLE        ::  &
    ! Output data from opt_th/opt_so
    podac   (:,:)   , & ! absorption optical depth
    podaf   (:,:)   , & ! in cloudy and free part
    podsc   (:,:)   , & ! scattering optical depth
    podsf   (:,:)   , & ! in cloudy and free part
    pbsfc   (:,:)   , & ! backscattering fraction 
    pbsff   (:,:)   , & ! in cloudy and free part
    pusfc   (:,:)   , & ! upscattering   fraction 
    pusff   (:,:)   , & ! in cloudy and free part

    ! cloud geometry factors
    pca1    (:,:)   , & !
    pcb1    (:,:)   , & !
    pcc1    (:,:)   , & !
    pcd1    (:,:)   , & !
    pca2    (:,:)   , & !
    pcb2    (:,:)   , & !
    pcc2    (:,:)   , & !
    pcd2    (:,:)   , & !

    !fluxes calculated in inv_th/inv_so
    pflfd   (:,:) , & ! (W/m**2)
    pflfu   (:,:) , & ! (W/m**2)
    pflfp   (:,:) , & ! (W/m**2)
    pflcd   (:,:) , & ! (W/m**2)
    pflcu   (:,:) , & ! (W/m**2)
    pflcp   (:,:)     ! (W/m**2)

!==============================================================================
! Working arrays - fesft_sp
!==============================================================================

  REAL    (KIND=dp), ALLOCATABLE ::  &
    ! Temperature at layer boundaries, pressure thickness of layers,
    ! cloud cover in each layer, water vapour and saturation water
    ! vapour mixing ratio, liquid water and ice mixing ratio
    ! layer CO2 content and O3 content
    pti_dp    (:,:)  , & ! (K)
    pdp_dp    (:,:)  , & ! (Pa)
    pclc_dp   (:,:)  , & ! (1)
    pwv_dp    (:,:)  , & ! (kg/kg)
    psw_dp    (:,:)  , & ! (kg/kg)
    pqlwc_dp  (:,:)  , & ! (kg/kg)
    pqiwc_dp  (:,:)  , & ! (kg/kg)
    pduco2_dp (:,:)  , & ! (Pa CO2)
    pduo3_dp  (:,:)  , & ! (Pa O3)

    ! Aerorsole optical depth at 0.55 micrometer for five types
    paeq1_dp  (:,:)  , & ! (1)
    paeq2_dp  (:,:)  , & ! (1)
    paeq3_dp  (:,:)  , & ! (1)
    paeq4_dp  (:,:)  , & ! (1)
    paeq5_dp  (:,:)  , & ! (1)

    ! Cosine of zenith angle, thermal and solar surface albedo
    psmu0_dp  (:)    , & ! (1)
    palth_dp  (:)    , & ! (1)
    palso_dp  (:)    , & ! (1)
  
    ! External data and radiation corrections
    pskyview_dp(:)   , & !
    pfcor_dp   (:)

  REAL    (KIND=dp), ALLOCATABLE ::  &
    ! Surface pressure
    papre_dp  (:)

! Output data
! -----------
  REAL    (KIND=dp), ALLOCATABLE ::  &
    ! Thermal and solar radiative fluxes at each layer boundary
    ! dito for cloud-free conditions (TOA and surface only)
    pflt_dp      (:,:) , & ! (W/m**2)
    pfls_dp      (:,:) , & ! (W/m**2)

    ! corrected thermal and solar surface flux
    ! (if not lradtopo, just pflt(ke1) and pfls(ke1)
    pflt_s_dp    (:)   , & ! (W/m**2)
    pfls_s_dp    (:)   , & ! (W/m**2)

    ! and for the Climate LM Version: solar direct downward radiative flux 
    pflsdir_dp   (:,:) , & ! (W/m**2)

    ! components of solar and thermal fluxes at surface
    ! (influenced by lradtopo for topographic corrections)
    pfltd_dp     (:)   , & ! (W/m**2)
    pfltu_dp     (:)   , & ! (W/m**2)
    pflsd_dp     (:)   , & ! (W/m**2)
    pflsu_dp     (:)   , & ! (W/m**2)
    pflsp_dp     (:)   , & ! (W/m**2)

    ! surface flux of photosynthetic active radiation and components
    pflpar_dp    (:)   , & ! (W/m**2)
    pflsu_par_dp (:)   , & ! (W/m**2)
    pflsd_par_dp (:)   , & ! (W/m**2)
    pflsp_par_dp (:)       ! (W/m**2)

!==============================================================================
! Working arrays - inv_th/so
!==============================================================================

  REAL    (KIND=dp), ALLOCATABLE ::  &

    ! layer properties calculated in *coe_th*
    pa1c (:)  , & !
    pa1f (:)  , & ! 
    pa2c (:)  , & ! 
    pa2f (:)  , & ! 
    pa3c (:)  , & ! 
    pa3f (:)  , & ! 
    pa4c (:)  , & ! 
    pa4f (:)  , & ! 
    pa5c (:)  , & ! 
    pa5f (:)  , & ! 

    ! Utility arrays
    ztu1 (:,:), & !
    ztu2 (:,:), & !
    ztu3 (:,:), & !
    ztu4 (:,:), & !
    ztu5 (:,:), & !
    ztu6 (:,:), & !
    ztu7 (:,:), & !
    ztu8 (:,:), & !
    ztu9 (:,:)    !

!==============================================================================
! Working arrays - opt_th/so
!==============================================================================

!    optical properties (absorption, scattering, back- and upscatter fraction)
!    for liquid water, ice and total aerosole
!    for Rayleigh scattering, only the optical depth is stored as array, since
!    back- and upscatter fractions are constant (=0.5)

  REAL    (KIND=dp), ALLOCATABLE ::  &
     zlwoda (:), & !
     zlwods (:), & ! 
     zlwb0  (:), & ! 
     zlwb   (:), & ! 
     ziwoda (:), & ! 
     ziwods (:), & ! 
     ziwb0  (:), & ! 
     ziwb   (:), & ! 
     zaeoda (:), & ! 
     zaeods (:), & ! 
     zaeb0  (:), & ! 
     zaeb   (:), & ! 
     zraods (:)    ! 

#ifdef COSMOART
  REAL    (KIND=dp), ALLOCATABLE ::  &
    asym_dust     (:,:,:), & !
    tau_abs_dust  (:,:,:), & !
    tau_scat_dust (:,:,:), & !
    asym_seas     (:,:,:), & !
    tau_abs_seas  (:,:,:), & !
    tau_scat_seas (:,:,:), & ! 
    asym_aero     (:,:,:), & ! 
    tau_abs_aero  (:,:,:), & !
    tau_scat_aero (:,:,:)    !
#endif

!==============================================================================

! indices for a debug point (define it here instead of doing it in all routines)
! Check the size of your domain, before activating debug output
  INTEGER, PARAMETER ::  &
     j1b    = 50,          & ! debug point index nproma dimension
     jb     = 4              ! debug point index block  dimension

! Interface Blocks
INTERFACE fesft
  MODULE PROCEDURE        &
    fesft_sp,             &
    fesft_dp
END INTERFACE

!==============================================================================

CONTAINS

!==============================================================================
!==============================================================================
!+ Module procedure in "Radiation" 
!------------------------------------------------------------------------------

SUBROUTINE fesft_sp(                                                     &
       pti   , pdp  , pclc  , pwv  , psw  , pqlwc, pqiwc, pduco2, pduo3, &
       paeq1 , paeq2, paeq3 , paeq4, paeq5,                              &
       papre , psmu0, palso , palth, pskyview, pfcor,                    &
       psig  , psct ,                                                    &
       ki1sd , ki1ed,                ki3sd, ki3ed,                       &
       ki1sc , ki1ec,                ki3sc, ki3ec,                       &
       lsolar, lcrf , lradtopo, idebug, jindex,                          &
       pflt  , pfls , pflt_s, pfls_s, pflsdir,                           &
       pfltd , pfltu, pflsd , pflsu, pflsp,                              &
       pflpar, pflsu_par, pflsd_par, pflsp_par,                          &
       ierror, yerrmsg)


!------------------------------------------------------------------------------
!
! Description:
!
!   This is the specific subroutine call for fesft with single precision 
!   variables. It transforms all input variables to double precision and then 
!   calls fesft_dp
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

! Input data
! ----------
  INTEGER, INTENT (IN) ::  &
    ! indices for array dimensioning
    ki1sd,       & ! start index for first  array dimension
    ki1ed,       & ! end   index for first  array dimension
    ki3sd,       & ! start index for third  array dimension
    ki3ed,       & ! end   index for third  array dimension

    ! and the same for the computations
    ki1sc,       & ! start index for first  array computation
    ki1ec,       & ! end   index for first  array computation
    ki3sc,       & ! start index for third  array computation
    ki3ec,       & ! end   index for third  array computation

    ! for activating debug output
    jindex,      & ! actual j-index computed
    idebug         ! debug control switch       

  LOGICAL                 , INTENT (IN) ::  &
    lcrf,        & ! control switch for cloud-free calcul.
    lsolar,      & ! control switch for solar calculations
    lradtopo       ! control switch for topographic corrections

  REAL    (KIND=sp)       , INTENT (IN) ::  &
    psig,        & ! Stefan-Boltzman constant 
    psct           ! solar constant (at time of year)

  REAL    (KIND=sp)       , INTENT (IN) ::  &
    ! Temperature at layer boundaries, pressure thickness of layers,
    ! cloud cover in each layer, water vapour and saturation water
    ! vapour mixing ratio, liquid water and ice mixing ratio
    ! layer CO2 content and O3 content
    pti    (ki1sd:ki1ed,            ki3sd:ki3ed+1)    , & ! (K)
    pdp    (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (Pa)
    pclc   (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (1)
    pwv    (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (kg/kg)
    psw    (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (kg/kg)
    pqlwc  (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (kg/kg)
    pqiwc  (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (kg/kg)
    pduco2 (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (Pa CO2)
    pduo3  (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (Pa O3)

    ! Aerorsole optical depth at 0.55 micrometer for five types
    paeq1  (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (1)
    paeq2  (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (1)
    paeq3  (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (1)
    paeq4  (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (1)
    paeq5  (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (1)

    ! Cosine of zenith angle, thermal and solar surface albedo
    psmu0  (ki1sd:ki1ed)                              , & ! (1)
    palth  (ki1sd:ki1ed)                              , & ! (1)
    palso  (ki1sd:ki1ed)                              , & ! (1)
  
    ! External data and radiation corrections
    pskyview(ki1sd:ki1ed)                             , & !
    pfcor   (ki1sd:ki1ed)            

  REAL    (KIND=sp)       , INTENT (IN)  ::  &
    ! Surface pressure
    papre  (ki1sd:ki1ed)                                  ! (Pa)

! Output data
! -----------
  REAL    (KIND=sp)       , INTENT (OUT) ::  &
    ! Thermal and solar radiative fluxes at each layer boundary
    ! dito for cloud-free conditions (TOA and surface only)
    pflt      (ki1sd:ki1ed,            ki3sd:ki3ed+1) , & ! (W/m**2)
    pfls      (ki1sd:ki1ed,            ki3sd:ki3ed+1) , & ! (W/m**2)

    ! corrected thermal and solar surface flux
    ! (if not lradtopo, just pflt(ke1) and pfls(ke1)
    pflt_s    (ki1sd:ki1ed)                           , & ! (W/m**2)
    pfls_s    (ki1sd:ki1ed)                           , & ! (W/m**2)

    ! and for the Climate LM Version: solar direct downward radiative flux 
    pflsdir   (ki1sd:ki1ed,            ki3sd:ki3ed+1) , & ! (W/m**2)

    ! components of solar and thermal fluxes at surface
    ! (influenced by lradtopo for topographic corrections)
    pfltd     (ki1sd:ki1ed)                           , & ! (W/m**2)
    pfltu     (ki1sd:ki1ed)                           , & ! (W/m**2)
    pflsd     (ki1sd:ki1ed)                           , & ! (W/m**2)
    pflsu     (ki1sd:ki1ed)                           , & ! (W/m**2)
    pflsp     (ki1sd:ki1ed)                           , & ! (W/m**2)

    ! surface flux of photosynthetic active radiation and components
    pflpar    (ki1sd:ki1ed)                           , & ! (W/m**2)
    pflsu_par (ki1sd:ki1ed)                           , & ! (W/m**2)
    pflsd_par (ki1sd:ki1ed)                           , & ! (W/m**2)
    pflsp_par (ki1sd:ki1ed)                               ! (W/m**2)

  INTEGER,          INTENT(OUT)  ::     ierror       ! error code
  CHARACTER(LEN=*), INTENT(OUT)  ::     yerrmsg      ! error message

!------------------------------------------------------------------------------

! Local variables in double precision
  REAL    (KIND=dp)                     ::  &
    psig_dp,     & ! Stefan-Boltzman constant 
    psct_dp        ! solar constant (at time of year)

  INTEGER :: ip,k

!------------------------------------------------------------------------------
! Start GPU data region
!------------------------------------------------------------------------------

  !$acc data &
  !---- Argument arrays - intent(in)
  !$acc present ( pti,pdp,pclc,pwv,psw,pqlwc,pqiwc,pduco2,pduo3  ) &
  !$acc present ( paeq1,paeq2,paeq3,paeq4,paeq5,psmu0,palso      ) &
  !$acc present ( palth,pskyview,pfcor                           ) &
  !---- Argument arrays - intent(in)  (formerly inout)           
  !$acc present ( papre                                          ) &
  !---- Argument arrays - intent(out)                            
  !$acc present ( pflt,pfls,pflt_s,pfls_s,pflsdir,pfltd,pfltu    ) &
  !$acc present ( pflsd,pflsu,pflsp,pflpar,pflsu_par,pflsd_par   ) &
  !$acc present ( pflsp_par                                      ) &
  !---- Local automatic arrays
  !$acc present ( pti_dp,pdp_dp,pclc_dp,pwv_dp,psw_dp,pqlwc_dp   ) &
  !$acc present ( pqiwc_dp,pduco2_dp,pduo3_dp,paeq1_dp,paeq2_dp  ) &
  !$acc present ( paeq3_dp,paeq4_dp,paeq5_dp,psmu0_dp,palth_dp   ) &
  !$acc present ( palso_dp,pskyview_dp,pfcor_dp,papre_dp,pflt_dp ) &
  !$acc present ( pfls_dp,pflt_s_dp,pfls_s_dp,pflsdir_dp         ) &
  !$acc present ( pfltd_dp,pfltu_dp,pflsd_dp,pflsu_dp,pflsp_dp   ) &
  !$acc present ( pflpar_dp,pflsu_par_dp,pflsd_par_dp            ) &
  !$acc present ( pflsp_par_dp                                   )

!------------------------------------------------------------------------------
! Begin of the subroutine
!------------------------------------------------------------------------------

  ierror     = 0
  yerrmsg(:) = ' '

! Variables with intent in
  psig_dp  = REAL(psig, dp)
  psct_dp  = REAL(psct, dp)

  !$acc parallel
  !$acc loop gang vector collapse(2)
  DO k=ki3sd,ki3ed+1
    DO ip=ki1sd,ki1ed
      pti_dp (ip,k) = 0.0_dp
      pflt   (ip,k) = 0.0_sp
      pfls   (ip,k) = 0.0_sp
      pflsdir(ip,k) = 0.0_sp
    ENDDO
  ENDDO
  !$acc end parallel

  !$acc parallel
  !$acc loop gang vector collapse(2)
  DO k=ki3sd,ki3ed
    DO ip=ki1sd,ki1ed
      pdp_dp    (ip,k) = 0.0_dp
      pclc_dp   (ip,k) = 0.0_dp
      pwv_dp    (ip,k) = 0.0_dp
      psw_dp    (ip,k) = 0.0_dp
      pqlwc_dp  (ip,k) = 0.0_dp
      pqiwc_dp  (ip,k) = 0.0_dp
      pduco2_dp (ip,k) = 0.0_dp
      pduo3_dp  (ip,k) = 0.0_dp
      paeq1_dp  (ip,k) = 0.0_dp
      paeq2_dp  (ip,k) = 0.0_dp
      paeq3_dp  (ip,k) = 0.0_dp
      paeq4_dp  (ip,k) = 0.0_dp
      paeq5_dp  (ip,k) = 0.0_dp
    ENDDO
  ENDDO
  !$acc end parallel

  !$acc parallel
  !$acc loop gang vector
  DO ip=ki1sd,ki1ed
    psmu0_dp   (ip) = 0.0_dp
    palth_dp   (ip) = 0.0_dp
    palso_dp   (ip) = 0.0_dp
    pskyview_dp(ip) = 0.0_dp
    pfcor_dp   (ip) = 0.0_dp
    pflt_s     (ip) = 0.0_sp
    pfls_s     (ip) = 0.0_sp
    pfltd      (ip) = 0.0_sp
    pfltu      (ip) = 0.0_sp
    pflsd      (ip) = 0.0_sp
    pflsu      (ip) = 0.0_sp
    pflsp      (ip) = 0.0_sp
    pflpar     (ip) = 0.0_sp
    pflsu_par  (ip) = 0.0_sp
    pflsd_par  (ip) = 0.0_sp
    pflsp_par  (ip) = 0.0_sp
  ENDDO
  !$acc end parallel

  !$acc parallel
  !$acc loop gang vector collapse(2)
  DO k=ki3sc,ki3ec+1
    DO ip=ki1sc,ki1ec
      pti_dp (ip,k) = REAL(pti (ip,k),dp)
    ENDDO
  ENDDO
  !$acc end parallel

  !$acc parallel
  !$acc loop gang vector collapse(2)
  DO k=ki3sc,ki3ec
    DO ip=ki1sc,ki1ec
      pdp_dp    (ip,k) = REAL(pdp    (ip,k),dp)
      pclc_dp   (ip,k) = REAL(pclc   (ip,k),dp)
      pwv_dp    (ip,k) = REAL(pwv    (ip,k),dp)
      psw_dp    (ip,k) = REAL(psw    (ip,k),dp)
      pqlwc_dp  (ip,k) = REAL(pqlwc  (ip,k),dp)
      pqiwc_dp  (ip,k) = REAL(pqiwc  (ip,k),dp)
      pduco2_dp (ip,k) = REAL(pduco2 (ip,k),dp)
      pduo3_dp  (ip,k) = REAL(pduo3  (ip,k),dp)

      paeq1_dp  (ip,k) = REAL(paeq1  (ip,k),dp)
      paeq2_dp  (ip,k) = REAL(paeq2  (ip,k),dp)
      paeq3_dp  (ip,k) = REAL(paeq3  (ip,k),dp)
      paeq4_dp  (ip,k) = REAL(paeq4  (ip,k),dp)
      paeq5_dp  (ip,k) = REAL(paeq5  (ip,k),dp)
    ENDDO
  ENDDO
  !$acc end parallel

  !$acc parallel
  !$acc loop gang vector
  DO ip=ki1sc,ki1ec
    psmu0_dp   (ip) = REAL(psmu0   (ip),dp)
    palth_dp   (ip) = REAL(palth   (ip),dp)
    palso_dp   (ip) = REAL(palso   (ip),dp)

    pskyview_dp(ip) = REAL(pskyview(ip),dp)
    pfcor_dp   (ip) = REAL(pfcor   (ip),dp)

  ! Variable with intent in (was inout formerly)
    papre_dp   (ip)  = REAL(papre  (ip),dp)
  ENDDO
  !$acc end parallel

! Now call double precision fesft
  CALL fesft_dp                                                               &
      (pti_dp,      pdp_dp,      pclc_dp,     pwv_dp,      psw_dp,            &
       pqlwc_dp,    pqiwc_dp,    pduco2_dp,   pduo3_dp,                       &
       paeq1_dp,    paeq2_dp,    paeq3_dp,    paeq4_dp,    paeq5_dp,          &
       papre_dp,    psmu0_dp,    palso_dp,    palth_dp,    pskyview_dp,       &
       pfcor_dp,    psig_dp,     psct_dp,                                     &
       ki1sd,       ki1ed,                                 ki3sd,      ki3ed, &
       ki1sc,       ki1ec,                                 ki3sc,      ki3ec, &
       lsolar,      lcrf,        lradtopo,    idebug,      jindex,            &
       pflt_dp,     pfls_dp,     pflt_s_dp,   pfls_s_dp,   pflsdir_dp,        &
       pfltd_dp,    pfltu_dp,    pflsd_dp,    pflsu_dp,    pflsp_dp,          &
       pflpar_dp,   pflsu_par_dp, pflsd_par_dp, pflsp_par_dp,                 &
       ierror,      yerrmsg)

! Transform output variables back to single precision
  !$acc parallel
  !$acc loop gang vector collapse(2)
  DO k=ki3sc,ki3ec+1
    DO ip=ki1sc,ki1ec
      pflt   (ip,k) = REAL(pflt_dp   (ip,k),sp)
      pfls   (ip,k) = REAL(pfls_dp   (ip,k),sp)
      pflsdir(ip,k) = REAL(pflsdir_dp(ip,k),sp)
    ENDDO
  ENDDO
  !$acc end parallel

  !$acc parallel
  !$acc loop gang vector
  DO ip=ki1sc,ki1ec
    pflt_s   (ip) = REAL(pflt_s_dp   (ip),sp)
    pfls_s   (ip) = REAL(pfls_s_dp   (ip),sp)
    pfltd    (ip) = REAL(pfltd_dp    (ip),sp)
    pfltu    (ip) = REAL(pfltu_dp    (ip),sp)
    pflsd    (ip) = REAL(pflsd_dp    (ip),sp)
    pflsu    (ip) = REAL(pflsu_dp    (ip),sp)
    pflsp    (ip) = REAL(pflsp_dp    (ip),sp)
    pflpar   (ip) = REAL(pflpar_dp   (ip),sp)
    pflsu_par(ip) = REAL(pflsu_par_dp(ip),sp)
    pflsd_par(ip) = REAL(pflsd_par_dp(ip),sp)
    pflsp_par(ip) = REAL(pflsp_par_dp(ip),sp)
  ! Variable with intent inout
  ! papre    (ip) = REAL(papre_dp    (ip),sp)
  ENDDO
  !$acc end parallel

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

  !$acc end data

END SUBROUTINE fesft_sp

!==============================================================================
!==============================================================================
!+ Module procedure in "Radiation" 
!------------------------------------------------------------------------------

SUBROUTINE fesft_dp(                                                     &
       pti   , pdp  , pclc  , pwv  , psw  , pqlwc, pqiwc, pduco2, pduo3, &
       paeq1 , paeq2, paeq3 , paeq4, paeq5,                              &
       papre , psmu0, palso , palth, pskyview, pfcor,                    &
       psig  , psct ,                                                    &
       ki1sd , ki1ed,                ki3sd, ki3ed,                       &
       ki1sc , ki1ec,                ki3sc, ki3ec,                       &
       lsolar, lcrf , lradtopo, idebug, jindex,                          &
       pflt  , pfls , pflt_s, pfls_s, pflsdir,                           &
       pfltd , pfltu, pflsd , pflsu, pflsp,                              &
       pflpar, pflsu_par, pflsd_par, pflsp_par,                          &
       ierror, yerrmsg)


!------------------------------------------------------------------------------
!
! Description:
!
!   The module procedure fesft organizes the radiative transfer calculations.
!
! Method:
!
! This routine organizes the radiative transfer calculations by
! calling a set of dedicated routines for the calculation of
! basic optical properties (opt_th/opt_so), the derivation of
! layer coefficients (coe_th/coe_so) for an implicit delta-two-
! stream scheme (cf.Ritter and Geleyn, 1992) and the inversion
! (inv_th/inv_so) of the corresponding system matrix. These 
! operations are performed seperately for thermal and solar parts
! of the spectrum and are embedded in loops over the various
! spectral intervals. Within each interval, a data-controlled
! decision is taken, whether the so-called ESFT or FESFT approach
! is used for the handling of gaseous absorption (cf. Ritter and
! Geleyn, 1992).
! Controlled by the logical input variable LCRF, the calculation
! of radiative fluxes in cloud-free conditions can be done in
! addition to the results for the given atmospheric cloud structure.
! (not implemented yet)
! Before the actual flux calculation starts, some preliminary steps
! provide utility arrays which are applicable to all spectral inter-
! vals (e.g. cloud geometry factors, integrated layer water content, etc.)
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

! Input data
! ----------
  INTEGER, INTENT (IN) ::  &
    ! indices for array dimensioning
    ki1sd,       & ! start index for first  array dimension
    ki1ed,       & ! end   index for first  array dimension
    ki3sd,       & ! start index for third  array dimension
    ki3ed,       & ! end   index for third  array dimension

    ! and the same for the computations
    ki1sc,       & ! start index for first  array computation
    ki1ec,       & ! end   index for first  array computation
    ki3sc,       & ! start index for third  array computation
    ki3ec,       & ! end   index for third  array computation

    ! for activating debug output
    jindex,      & ! actual j-index computed
    idebug         ! debug control switch       

  LOGICAL                 , INTENT (IN) ::  &
    lcrf,        & ! control switch for cloud-free calcul.
    lsolar,      & ! control switch for solar calculations
    lradtopo       ! control switch for topographic corrections

  REAL    (KIND=dp)       , INTENT (IN) ::  &
    psig,        & ! Stefan-Boltzman constant 
    psct           ! solar constant (at time of year)

  REAL    (KIND=dp)       , INTENT (IN) ::  &
    ! Temperature at layer boundaries, pressure thickness of layers,
    ! cloud cover in each layer, water vapour and saturation water
    ! vapour mixing ratio, liquid water and ice mixing ratio
    ! layer CO2 content and O3 content
    pti    (ki1sd:ki1ed,            ki3sd:ki3ed+1)    , & ! (K)
    pdp    (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (Pa)
    pclc   (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (1)
    pwv    (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (kg/kg)
    psw    (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (kg/kg)
    pqlwc  (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (kg/kg)
    pqiwc  (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (kg/kg)
    pduco2 (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (Pa CO2)
    pduo3  (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (Pa O3)

    ! Aerorsole optical depth at 0.55 micrometer for five types
    paeq1  (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (1)
    paeq2  (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (1)
    paeq3  (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (1)
    paeq4  (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (1)
    paeq5  (ki1sd:ki1ed,            ki3sd:ki3ed)      , & ! (1)

    ! Cosine of zenith angle, thermal and solar surface albedo
    psmu0  (ki1sd:ki1ed)                              , & ! (1)
    palth  (ki1sd:ki1ed)                              , & ! (1)
    palso  (ki1sd:ki1ed)                              , & ! (1)
  
    ! External data and radiation corrections
    pskyview(ki1sd:ki1ed)                             , & !
    pfcor   (ki1sd:ki1ed)                

  REAL    (KIND=dp)       , INTENT (IN)  ::  &
    ! Surface pressure
    papre  (ki1sd:ki1ed)                                  ! (Pa)

! Output data
! -----------
  REAL    (KIND=dp)       , INTENT (OUT) ::  &
    ! Thermal and solar radiative fluxes at each layer boundary
    ! dito for cloud-free conditions (TOA and surface only)
    pflt      (ki1sd:ki1ed,            ki3sd:ki3ed+1) , & ! (W/m**2)
    pfls      (ki1sd:ki1ed,            ki3sd:ki3ed+1) , & ! (W/m**2)

    ! corrected thermal and solar surface flux
    ! (if not lradtopo, just pflt(ke1) and pfls(ke1)
    pflt_s    (ki1sd:ki1ed)                           , & ! (W/m**2)
    pfls_s    (ki1sd:ki1ed)                           , & ! (W/m**2)

    ! and for the Climate LM Version: solar direct downward radiative flux 
    pflsdir   (ki1sd:ki1ed,            ki3sd:ki3ed+1) , & ! (W/m**2)

    ! components of solar and thermal fluxes at surface
    ! (influenced by lradtopo for topographic corrections)
    pfltd     (ki1sd:ki1ed)                           , & ! (W/m**2)
    pfltu     (ki1sd:ki1ed)                           , & ! (W/m**2)
    pflsd     (ki1sd:ki1ed)                           , & ! (W/m**2)
    pflsu     (ki1sd:ki1ed)                           , & ! (W/m**2)
    pflsp     (ki1sd:ki1ed)                           , & ! (W/m**2)

    ! surface flux of photosynthetic active radiation and components
    pflpar    (ki1sd:ki1ed)                           , & ! (W/m**2)
    pflsu_par (ki1sd:ki1ed)                           , & ! (W/m**2)
    pflsd_par (ki1sd:ki1ed)                           , & ! (W/m**2)
    pflsp_par (ki1sd:ki1ed)                               ! (W/m**2)

  INTEGER,          INTENT(OUT)  ::     ierror       ! error code
  CHARACTER(LEN=*), INTENT(OUT)  ::     yerrmsg      ! error message

! Local parameters: 
! ----------------
  REAL    (KIND=dp)       , PARAMETER ::  &
    zepflx = 1.0E-8_dp,         & ! Minimum 'grey' flux to avoid 1./0.
    zrd    = 287.05_dp,         & ! Ra (gas constant of dry air)

    zrvdm1 = 461.51_dp/287.05_dp-1.0_dp,  & ! Rv/Ra - 1
    zrvd   = 461.51_dp/287.05_dp,         & ! Rv/Ra
    zepai  = 0.0_dp               ! Could be used to save computing time for 
                                  ! 'unimportant' gaseous absorption coefficients

! Local scalars:
! -------------
  INTEGER  ::  &
    icrf, igase,           & !
    igasm1, igasz,         & !
    j1,j3,                 & ! loop indices over spatial dimensions
    jc,jh2o,jco2,jo3,      & ! loop indices over gaseous coefficients
    jspec,jspect,          & ! loop indices over spectrum
    jg,jjg ,               & ! loop indices over gases      
    icc,ih2o,ico2,io3        ! loop limit for gaseous absorption

  LOGICAL                  ::  &
    ldebug_th      ,       & ! debug control switch for thermal     
    ldebug_so      ,       & ! debug control switch for solar       
    ldebug_opt_th  ,       & ! debug control switch for opt_th      
    ldebug_opt_so  ,       & ! debug control switch for opt_so      
    ldebug_inv_th  ,       & ! debug control switch for inv_th      
    ldebug_inv_so            ! debug control switch for inv_so      

  REAL    (KIND=dp)        ::  &
    zet ,zaiprod,          & !
    zcoai,zcobi,           & !
    zemissivity, zalbedo     !

! Local arrays:
! -------------
  INTEGER  ::  &
    icgas (3)                !

  CHARACTER (LEN=  5)  yzroutine

!- End of header
!==============================================================================
 
!------------------------------------------------------------------------------
! Start GPU data region
!------------------------------------------------------------------------------

  !$acc data                                                      &
  !---- Argument arrays - intent(in)
  !$acc present ( pti,pdp,pclc,pwv,psw,pqlwc,pqiwc,pduco2,pduo3 ) &
  !$acc present ( paeq1,paeq2,paeq3,paeq4,paeq5,psmu0,palth     ) &
  !$acc present ( palso,pskyview,pfcor                          ) &
  !---- Argument arrays - intent(inout)
  !$acc present ( papre                                         ) &
  !---- Argument arrays - intent(in) (formerly inout)
  !$acc present ( pflt,pfls,pflt_s,pfls_s,pflsdir,pfltd,pfltu   ) &
  !$acc present ( pflsd,pflsu,pflsp,pflpar,pflsu_par,pflsd_par  ) &
  !$acc present ( pflsp_par                                     ) &
  !---- Local automatic arrays
  !$acc present ( zketyp,ztetyp,zflux,zflux_c,zfluxi,zfluxu     ) &
  !$acc present ( zfluxu_c,zfluxui,zfluxd,zfluxd_c,zfluxdi      ) &
  !$acc present ( zfgas,zfgasu,zfgasd,pbbr,pflpt,palp,pqsmu0    ) &
  !$acc present ( palogt,palogp,papra,pduh2oc,pduh2of,pdulwc    ) &
  !$acc present ( pduiwc,prholwc,prhoiwc,zduetpc,zduetpf,zapre  ) &
  !$acc present ( ztm,zzwv,zcpo,zcpn,zcmo,zcmn,podac,podaf,podsc) &
  !$acc present ( podsf,pbsfc,pbsff,pusfc,pusff,pca1,pcb1,pcc1  ) &
  !$acc present ( pcd1,pca2,pcb2,pcc2,pcd2,pflfd,pflfu,pflfp    ) &
  !$acc present ( pflcd,pflcu,pflcp                             )

!------------------------------------------------------------------------------
! Begin Subroutine fesft               
!------------------------------------------------------------------------------

  yzroutine  = 'fesft'
  yerrmsg(:) = ' '
  ierror     = 0

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  icrf  = 0
 
  ! Debug switches for lower level subroutines
  ldebug_th     = .FALSE.
  ldebug_so     = .FALSE.
  ldebug_opt_th = .FALSE.
  ldebug_opt_so = .FALSE.
  ldebug_inv_th = .FALSE.
  ldebug_inv_so = .FALSE.

  IF (idebug > 15) THEN
    PRINT *,' **** FESFT  *********************** debug point : ',j1b,jb
  ENDIF 

  ! Preset output arrays
  !$acc parallel
  !$acc loop gang vector collapse(2)
  DO j3 = ki3sc, ki3ec+1
    DO j1 = ki1sc, ki1ec
      pflt   (j1,   j3) = 0.0_dp
      pfls   (j1,   j3) = 0.0_dp
      ! and for Climate-LM Version
      pflsdir(j1,   j3) = 0.0_dp
    ENDDO
  ENDDO
  !$acc end parallel

#ifdef COSMOART
  IF(l_cosmo_art .AND. lgas) THEN
    !$acc parallel
    !$acc loop vector vector collapse(2)
    DO j3 = ki3sc, ki3ec+1
      DO j1 = ki1sc, ki1ec  
        Edir_b (j1,j3) =  0.0_dp
        Edown_b(j1,j3) =  0.0_dp
        Eup_b  (j1,j3) =  0.0_dp
      ENDDO
    ENDDO
    !$acc end parallel
  ENDIF
#endif

  !$acc parallel
  !$acc loop gang vector
  DO j1 = ki1sc, ki1ec
    pflpar   (j1) = 0.0_dp
    pflsu_par(j1) = 0.0_dp
    pflsd_par(j1) = 0.0_dp
    pflsp_par(j1) = 0.0_dp
    pfls_s   (j1) = 0.0_dp
    pflt_s   (j1) = 0.0_dp
    pfltu    (j1) = 0.0_dp
    pfltd    (j1) = 0.0_dp
    pflsu    (j1) = 0.0_dp
    pflsd    (j1) = 0.0_dp
    pflsp    (j1) = 0.0_dp
  ENDDO
  !$acc end parallel
 
  ! Choice of e-type-absorption and temperature correction coefficients
  DO jspec = 1, jpther     
     zketyp(jspec) = zketypa(jspec)
     ztetyp(jspec) = ztetypa(jspec)
  ENDDO
  ! Initialized on CPU to avoid copying zketypa,ztetypa to GPU,
  ! but zketyp,ztetyp subsequently used on both CPU and GPU.
  !$acc update device (zketyp,ztetyp)


  ! cloud geometry factors
 
  ! first part for top layer
  !$acc parallel
  !$acc loop gang vector
  DO j1 = ki1sc, ki1ec
    pca2(j1,ki3sc)   = 1.0_dp - pclc(j1   ,ki3sc)
    pcd2(j1,ki3sc)   = 1.0_dp
    zcpn(j1)         = MAX(pclc(j1   ,ki3sc),pclc(j1   ,ki3sc+1))
    zcmn(j1)         = MIN(pclc(j1   ,ki3sc),pclc(j1   ,ki3sc+1))
    pca2(j1,ki3sc+1) = (1.0_dp-zcpn(j1   ))/pca2(j1   ,ki3sc)
    pcd2(j1,ki3sc+1) =      zcmn(j1   ) /pclc(j1   ,ki3sc)
  ENDDO
  !$acc end parallel
 
  ! first part for inner layers
 
  !$acc parallel
  !$acc loop seq
  DO j3 = ki3sc+1, ki3ec-1
    !$acc loop gang vector
    DO j1 = ki1sc, ki1ec
      zcpo(j1)      = zcpn(j1   )
      zcmo(j1)      = zcmn(j1   )
      zcpn(j1)      = MAX (pclc(j1   ,j3),pclc(j1   ,j3+1))
      zcmn(j1)      = MIN (pclc(j1   ,j3),pclc(j1   ,j3+1))
      pca2(j1,j3+1) = (1.0_dp-zcpn(j1   ))/(1.0_dp-pclc(j1   ,j3))
      pca1(j1,j3-1) = (1.0_dp-zcpo(j1   ))/(1.0_dp-pclc(j1   ,j3))
      pcd2(j1,j3+1) = zcmn(j1   )/pclc(j1   ,j3)
      pcd1(j1,j3-1) = zcmo(j1   )/pclc(j1   ,j3)
    ENDDO
  ENDDO
  !$acc end parallel
 
  ! first part for lowest layer
  !$acc parallel
  !$acc loop gang vector
  DO j1 = ki1sc, ki1ec
    pca1(j1,ki3ec-1) = (1.0_dp-zcpn(j1   ))/(1.0_dp-pclc(j1   ,ki3ec))
    pcd1(j1,ki3ec-1) = zcmn(j1   )/pclc(j1   ,ki3ec)
    pca1(j1,ki3ec  ) = 1.0_dp
    pcd1(j1,ki3ec  ) = 1.0_dp
  ENDDO
  !$acc end parallel
 
  ! second part of geometry factors
  !$acc parallel
  !$acc loop gang vector collapse(2)
  DO j3 = ki3sc, ki3ec
    DO j1 = ki1sc, ki1ec
      pcb1(j1,j3) = 1.0_dp-pca1(j1   ,j3)
      pcc1(j1,j3) = 1.0_dp-pcd1(j1   ,j3)
      pcb2(j1,j3) = 1.0_dp-pca2(j1   ,j3)
      pcc2(j1,j3) = 1.0_dp-pcd2(j1   ,j3)
    ENDDO
  ENDDO
  !$acc end parallel
 
  ! Optically relevant layer constituents
  ! (Note: CO2 and O3 amounts are provided by calling routine)
  !$acc parallel
  !$acc loop gang vector
  DO j1 = ki1sc, ki1ec
    papra(j1) = papre(j1   )   ! surface pressure
    zapre(j1) = papre(j1   )   ! surface pressure
  ENDDO
  !$acc end parallel
 
  ! water vapour, liquid water and ice content, logarithm of layer
  ! mean temperature and pressure, absorber amount for e-type absorption
 
  !$acc parallel
  !$acc loop seq
  DO j3 = ki3ec, ki3sc,-1    ! Bottom to top
    !$acc loop gang vector
    DO j1 = ki1sc, ki1ec
      ztm   (j1)    = 0.5_dp*(pti(j1   ,j3)+pti(j1   ,j3+1))
      papra (j1)    = papra(j1   ) - 0.5_dp*pdp(j1   ,j3)
      palogt(j1,j3) = LOG  (ztm  (j1   ))
      palogp(j1,j3) = LOG  (papra(j1   ))
  
      ! cloud-free:  water vapour and e-type absorber amount
      zzwv   (j1) = MAX( (pwv(j1   ,j3)-pclc(j1   ,j3)*psw(j1   ,j3)) &
                            /(1.0_dp-pclc(j1   ,j3)) , 0.0_dp)
      pduh2of(j1,j3) = pdp(j1   ,j3)*zzwv(j1   )
      zduetpf(j1,j3) = pduh2of(j1   ,j3)*pduh2of(j1   ,j3) &
                         *papra(j1   )*zrvd/pdp(j1   ,j3)
  
      ! cloudy:  water vapour, e-type absorber amount, liquid water and ice
      pdulwc (j1,j3) = pdp(j1   ,j3) * (pqlwc(j1   ,j3)/pclc(j1   ,j3))
      pdulwc (j1,j3) = MAX( pdulwc(j1   ,j3), 0.0_dp )
      pduiwc (j1,j3) = pdp(j1   ,j3) * (pqiwc(j1   ,j3)/pclc(j1   ,j3))
      pduiwc (j1,j3) = MAX( pduiwc(j1   ,j3), 0.0_dp )
      pduh2oc(j1,j3) = pdp(j1   ,j3) * psw(j1   ,j3)
      zduetpc(j1,j3) = pduh2oc(j1   ,j3) * pduh2oc(j1   ,j3) *           &
                                papra(j1   ) * zrvd / pdp(j1   ,j3)
      prholwc(j1,j3) = (pqlwc(j1   ,j3) / pclc(j1   ,j3)) * papra(j1   ) &
                            / (zrd*ztm(j1   ) * (1.0_dp+zrvdm1*psw(j1   ,j3)))
      prhoiwc(j1,j3) = (pqiwc(j1   ,j3) / pclc(j1   ,j3)) * papra(j1   ) &
                            / (zrd*ztm(j1   ) * (1.0_dp+zrvdm1*psw(j1   ,j3)))
      ! Secure minium for ice density for use in empirical function with ALOG
      prhoiwc(j1,j3) = MAX (prhoiwc(j1   ,j3),1.0E-06_dp) 
       
      papra  (j1)    = papra(j1) - 0.5_dp * pdp(j1,j3)

    ENDDO
  ENDDO     ! End of vertical loop
  !$acc end parallel

  ! Identify *zapre* with top of model pressure (for Rayleigh scattering)
  !$acc parallel
  !$acc loop gang vector
  DO j1 = ki1sc, ki1ec
    zapre(j1) = papra(j1)
  ENDDO
  !$acc end parallel
 
1 CONTINUE  ! Address for backward jump to perform cloud-free calculations

!------------------------------------------------------------------------------
! Section 2: Thermal radiative flux calculations 
!------------------------------------------------------------------------------

  ! Loop over thermal spectral intervals
  !================================================================
  thermal_spectral_loop:  DO jspec=jpsol+1,jpspec
  !================================================================
      
  !----------------------------------------------------------------------------
  ! Section 2.1: Initializations
  !----------------------------------------------------------------------------

    jspect = jspec - jpsol
 
    ! Black body radiation at layer boundaries in spectral interval
    !$acc parallel
    !$acc loop gang
    DO j3 = ki3sc, ki3ec+1
      !$acc loop vector
      DO j1 = ki1sc, ki1ec
        pbbr(j1,j3)= ( planck(1,jspect) + pti(j1,j3)                      &
                   * ( planck(2,jspect) + pti(j1,j3)*planck(3,jspect) ) ) &
                   * psig * (pti(j1,j3)**2)**2
      ENDDO
    ENDDO
    !$acc end parallel
 
    ! Optical properties of non-gaseous constituents        
    IF (idebug > 10 ) THEN
       print *,' FESFT    Call to opt_th for jspec: ',jspec
    ENDIF 

    CALL opt_th( prholwc, pdulwc, prhoiwc, pduiwc,              &
                 paeq1  , paeq2 , paeq3  , paeq4 , paeq5 ,      &
                 ki1sd  , ki1ed , ki3sd  , ki3ed , jspec ,      &
                 ki1sc  , ki1ec , ki3sc  , ki3ec ,              &
                 ldebug_opt_th  , jindex ,                      &
                 podac  , podaf , podsc  , podsf , pbsfc , pbsff)

    ! Addition of e-type contribution
    IF (zketyp(jspect) /= 0.0_dp) THEN
      zet   = 1.0_dp/EXP(ztetyp(jspect)/zteref)
      !$acc parallel
      !$acc loop seq
      DO j3 = ki3sc,ki3ec
        !$acc loop gang vector
        DO j1 = ki1sc, ki1ec
          ztm(j1)      = 0.5_dp*(pti(j1,j3)+pti(j1,j3+1))
          podaf(j1,j3) = podaf(j1,j3) + zduetpf(j1,j3) &
                       * zet * EXP(ztetyp(jspect)/ztm(j1)) * zketyp(jspect)
          podac(j1,j3) = podac(j1,j3) + zduetpc(j1,j3) &
                       * zet * EXP(ztetyp(jspect)/ztm(j1)) * zketyp(jspect)
        ENDDO  
      ENDDO  
      !$acc end parallel
    ENDIF

  !----------------------------------------------------------------------------
  ! Section 2.2: Selection of ESFT or FESFT method for interval considered
  !----------------------------------------------------------------------------

    !--------------------------------------------------------------
    IF (nfast(jspec) == 0) THEN           ! ESFT method
    !--------------------------------------------------------------

      !$acc parallel
      !$acc loop gang vector collapse(2)
      DO j3 = ki3sc, ki3ec+1
        DO j1 = ki1sc, ki1ec
          zflux(j1,j3) = 0.0_dp  ! Preset flux in spectral interval
        ENDDO
      ENDDO
      !$acc end parallel
 
      ! Loop over various absorption coefficients of each of the three gases
 
      ih2o = ncgas(jspec,1)
      ico2 = ncgas(jspec,2)
      io3  = ncgas(jspec,3)
 
      DO jh2o=    1,ih2o    ! Loop over H2O coefficients
        DO jco2=  1,ico2    ! Loop over CO2 coefficients
          DO jo3= 1,io3     ! Loop over O3  coefficients
 
            zaiprod = coai(jh2o,jspec,1)*coai(jco2,jspec,2)*coai(jo3,jspec,3)
 
            IF (icrf.eq.0) THEN      ! partially cloudy atmosphere
              CALL       inv_th (                                              &
                pclc   ,pca1   ,pca2  ,pcb1  ,pcb2  ,pcc1  ,pcc2 ,pcd1 ,pcd2 , &
                pduh2oc,pduh2of,pduco2,pduo3 ,palogp,palogt,                   &
                podsc  ,podsf  ,podac ,podaf ,pbsfc ,pbsff ,                   &
                pbbr   ,palth,                                                 &
                jspec  ,jh2o   ,jco2  ,jo3   ,                                 &
                ki1sd  ,ki1ed  ,ki3sd ,ki3ed ,ki1sc ,ki1ec ,ki3sc ,ki3ec ,     &
                ldebug_inv_th  ,jindex,                                        &
                pflcu  ,pflfu  ,pflcd ,pflfd)
            ELSE                     ! 'cloud-free' atmosphere    
              print *,' CRF not yet implemented'
            END IF
 
            ! Incrementation of flux in spectral interval

            !$acc parallel
            !$acc loop gang
            DO j3 = ki3sc, ki3ec+1
              !$acc loop vector
              DO j1 = ki1sc, ki1ec
                zflux(j1,j3) = zflux(j1,j3)                            &
                   + zaiprod * ( pflfu(j1,j3) + pflcu(j1,j3)           &
                               - pflfd(j1,j3) - pflcd(j1,j3) )
              ENDDO
            ENDDO
            !$acc end parallel

          ENDDO       ! Loop over O3 absorption coefficients
        ENDDO         ! Loop over CO2 absorption coefficients
      ENDDO           ! Loop over H2O absorption coefficients

    !--------------------------------------------------------------
    ELSE                                  ! FESFT method
    !--------------------------------------------------------------

      igase = 0
      DO jg=1,3
        icgas (jg) = 0
        IF (ncgas(jspec,jg).GT.1) THEN
          igase     = igase + 1
        ENDIF
      ENDDO
      igasm1    = igase -1
 
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF (igase.le.1) THEN  !(no 'grey' fluxes required)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        !$acc parallel
        !$acc loop gang
        DO j3 = ki3sc, ki3ec+1
          !$acc loop vector
          DO j1 = ki1sc, ki1ec
            zfluxi(j1,j3) = 1.0_dp
          ENDDO
        ENDDO
        !$acc end parallel

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE   ! more than 1 gas --> 'grey' fluxes required
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        IF (icrf.EQ.0) THEN    ! partially cloudy atmosphere
          CALL       inv_th (                                               &
             pclc   ,pca1   ,pca2  ,pcb1  ,pcb2  ,pcc1  ,pcc2 ,pcd1 ,pcd2 , &
             pduh2oc,pduh2of,pduco2,pduo3 ,palogp,palogt,                   &
             podsc  ,podsf  ,podac ,podaf ,pbsfc ,pbsff ,                   &
             pbbr   ,palth,                                                 &
             jspec  ,0      ,0     ,0     ,                                 &
             ki1sd  ,ki1ed  ,ki3sd ,ki3ed ,ki1sc ,ki1ec ,ki3sc ,ki3ec ,     &
             ldebug_inv_th  ,jindex,                                        &
             pflcu  ,pflfu  ,pflcd ,pflfd)
        ELSE                     ! 'cloud-free' atmosphere    
          print *,' CRF not yet implemented'
        ENDIF

        ! Storage of 'grey' fluxes and their inverse (**igasm1)
        !$acc parallel
        !$acc loop gang
        DO j3 = ki3sc, ki3ec+1
          !$acc loop vector
          DO j1 = ki1sc, ki1ec
            zflux (j1,j3) = pflfu(j1,j3) + pflcu(j1,j3)          &
                          - pflfd(j1,j3) - pflcd(j1,j3)
            zfluxi(j1,j3) = 1.0_dp / MIN( -zepflx, zflux(j1,j3) )**igasm1
          ENDDO
        ENDDO
        !$acc end parallel

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ENDIF    ! No.of relevant gases
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      igasz = 0
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO jg = 3, 1, -1   !     Loop over gases
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
        IF (ncgas(jspec,jg).GT.1) THEN   ! include gas only, if necessary
          igasz = igasz + 1
 
          DO jjg = 1,3
            icgas(jjg) = 0      !   Set absorption coefficient index for all
          ENDDO                 !   gases to zero
 
          !$acc parallel
          !$acc loop gang
          DO j3 = ki3sc, ki3ec+1
            !$acc loop vector
            DO j1 = ki1sc, ki1ec
              zfgas(j1,j3) = 0.0_dp    ! Preset 'gaseous' flux
            ENDDO
          ENDDO
          !$acc end parallel
 
          icc = ncgas(jspec,jg) ! No.of relevant coefficients    
 
          ! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          DO jc = icc,1,-1        ! Loop over absorption coefficients
          ! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

            zcoai = coai(jc,jspec,jg)
            zcobi = cobi(jc,jspec,jg)
 
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            IF ( ((zcoai.GE.zepai).AND.(zcobi.GT.0.0_dp)) .OR. (igase.EQ.1) ) THEN   
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

              ! Solve linear system, if necessary
              icgas(jg) = jc
              IF (icrf.EQ.0) THEN  ! partially cloudy atmosphere

                CALL       inv_th (                                              & 
                  pclc   ,pca1   ,pca2  ,pcb1  ,pcb2  ,pcc1  ,pcc2 ,pcd1 ,pcd2 , &
                  pduh2oc,pduh2of,pduco2,pduo3 ,palogp,palogt,                   &
                  podsc  ,podsf  ,podac ,podaf ,pbsfc ,pbsff ,                   &
                  pbbr   ,palth,                                                 &
                  jspec  ,icgas(1),icgas(2),icgas(3),                            &
                  ki1sd  ,ki1ed  ,ki3sd ,ki3ed ,ki1sc ,ki1ec ,ki3sc ,ki3ec ,     &
                  ldebug_inv_th  ,jindex,                                        &
                  pflcu  ,pflfu  ,pflcd ,pflfd)
              ELSE                     ! 'cloud-free' atmosphere    
                print *,' CRF not yet implemented'
              ENDIF
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ELSE               ! use 'grey' fluxes directly
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

              ! Solve linear system, if necessary
              !$acc parallel
              !$acc loop gang vector collapse(2)
              DO j3 = ki3sc, ki3ec+1
                DO j1 = ki1sc, ki1ec
                  pflfu(j1,j3) = zflux(j1,j3)
                  pflfd(j1,j3) = 0.0_dp
                  pflcu(j1,j3) = 0.0_dp
                  pflcd(j1,j3) = 0.0_dp
                ENDDO
              ENDDO
              !$acc end parallel
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ENDIF            ! Necessity to calculate fluxes
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
            !$acc parallel
            !$acc loop gang
            DO j3 = ki3sc, ki3ec+1
              !$acc loop vector
              DO j1 = ki1sc, ki1ec
                zfgas(j1,j3) = zfgas(j1,j3) + zcoai * ( pflfu(j1,j3) &
                             + pflcu(j1,j3) - pflfd(j1,j3) - pflcd(j1,j3) )
              ENDDO
            ENDDO
            !$acc end parallel

#ifndef _OPENACC
            !?????????????????????????????????????????????????????????????????
            IF (ldebug_th .AND. jb==jindex) THEN
               print *,' FESFT in debug mode for thermal fluxes'
               print *,' only one interval/coefficient considered '
               print *,'zfgas(j1b,jb,ki3sc): ',zfgas(j1b    ,ki3sc)
               EXIT thermal_spectral_loop
            ENDIF 
            !?????????????????????????????????????????????????????????????????
#endif
          ! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          ENDDO      ! Loop over absorption coefficients
          ! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
 
          ! Combination of inverse of 'grey' fluxes and actual gaseous flux
          !$acc parallel
          !$acc loop gang
          DO j3 = ki3sc, ki3ec+1
            !$acc loop vector
            DO j1 = ki1sc, ki1ec
              zfluxi(j1,j3) = zfluxi(j1,j3)*zfgas(j1,j3)
            ENDDO
          ENDDO
          !$acc end parallel
 
          IF (igasz.eq.igasm1) THEN    ! Avoid unphysical pseudo-transmission
            !$acc parallel
            !$acc loop gang
            DO j3 = ki3sc, ki3ec+1
              !$acc loop vector
              DO j1 = ki1sc, ki1ec
                zfluxi(j1,j3) = MIN( 1.0_dp, MAX( 0.0_dp, zfluxi(j1,j3) ) )
              ENDDO
            ENDDO
            !$acc end parallel
          ENDIF
    
        ENDIF                   ! Test, whether gas needs to be included
    
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END DO          ! End of loop over gases
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
      ! Store FESFT result in zflux
      !$acc parallel
      !$acc loop gang
      DO j3 = ki3sc, ki3ec+1
        !$acc loop vector
        DO j1 = ki1sc, ki1ec
          zflux(j1,j3) = zfluxi(j1,j3)
        ENDDO
      ENDDO
      !$acc end parallel

    !--------------------------------------------------------------
    END IF                             ! ESFT/FESFT-Selection 
    !--------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Section 2.3: Addition of flux for spectral interval to total thermal flux
  !----------------------------------------------------------------------------
 
    IF (icrf.EQ.0) THEN     ! Add flux at all levels
      !$acc parallel
      !$acc loop gang
      DO j3 = ki3sc, ki3ec+1
        DO j1 = ki1sc, ki1ec
          pflt(j1,j3) = pflt(j1,j3) + zflux(j1,j3)
        ENDDO
      ENDDO
      !$acc end parallel
    END IF
 
  ! End of spectral loop
  ! ===============================================================
  ENDDO  thermal_spectral_loop
  ! ===============================================================

  !----------------------------------------------------------------------------
  ! Section 2.4: Storage of components for thermal radiative surface flux
  !----------------------------------------------------------------------------

  !$acc parallel
  !$acc loop gang vector
  DO j1 = ki1sc, ki1ec
    ! Recompute surface thermal flux components based on lower boundary 
    ! condition (cf. Ritter and Geleyn (1992))
    zemissivity =  1.0_dp-palth(j1) ! surface emissivity
    pfltd (j1) = (pflt(j1,ki3ec+1) + zemissivity * psig * pti(j1,ki3ec+1)**4) &
               / zemissivity
    pfltu (j1) = pfltd(j1) - pflt(j1,ki3ec+1)
    pflt_s(j1) = pflt(j1,ki3ec+1)
  ENDDO
  !$acc end parallel

#ifndef _OPENACC
!NEC_CB Moved Debug-Prints of out the loop
  IF (idebug > 15) THEN
    IF ((jb==jindex).AND.(j1b>=ki1sc).AND.(j1b<=ki1ec)) THEN
      WRITE (*,'(A32,2F16.6)') 'FESFT: zemissivity, palth',              &
           zemissivity, palth(j1b)  
      WRITE (*,'(A60, F16.6)')                                           &
          'FESFT: thermal fluxes before and after correction, skyview',  &
           pskyview(j1b)
      WRITE (*,'(A32,2I4,3F16.6)') 'FESFT th before: net, down, up:',    &
           j1b, jb, pflt(j1b,ki3ec+1), pfltd(j1b), pfltu(j1b)
    ENDIF
  ENDIF
#endif

  IF (lradtopo) THEN
    !$acc parallel
    !$acc loop gang vector
    DO j1 = ki1sc, ki1ec
      ! corrected thermal balance
      ! correction as in Mueller and Scherrer (2005)
      pfltd (j1) = pfltd(j1) *           pskyview(j1)  +   &
                   pfltu(j1) * (1.0_dp - pskyview(j1))
      pflt_s(j1) = pfltd(j1)-pfltu(j1)
    ENDDO
    !$acc end parallel

#ifndef _OPENACC
!NEC_CB Moved Debug-Prints of out the loop
    IF (idebug > 15) THEN
      IF ((jb==jindex).AND.(j1b>=ki1sc).AND.(j1b<=ki1ec)) THEN
        WRITE (*,'(A32,2I4,3F16.6)') 'FESFT th after: net,down,up:',     &
                j1b, jb, pflt_s(j1b  ), pfltd(j1b  ), pfltu(j1b  )
      END IF
    END IF
#endif
  END IF

!------------------------------------------------------------------------------
! Section 3: Solar flux calculations, if required
!------------------------------------------------------------------------------

IF (lsolar) THEN

  !----------------------------------------------------------------------------
  ! Section 3.1: Initializations
  !----------------------------------------------------------------------------

  ! Inverse of cosine of zenith angle and surface albedo for
  ! parallel radiation

  !$acc parallel
  !$acc loop gang vector
  DO j1 = ki1sc, ki1ec
    pqsmu0(j1) = 1.0_dp / psmu0(j1)
    palp  (j1) = (1.0_dp +                                        &
             0.5_dp * (psmu0(j1) * (1.0_dp/palso(j1) - 1.0_dp)))  &
          / (1.0_dp + (psmu0(j1) * (1.0_dp/palso(j1) - 1.0_dp)))**2
  ENDDO
  !$acc end parallel

  ! Loop over solar spectral intervals
  ! ===============================================================
  solar_spectral_loop:  DO jspec = 1, jpsol
  ! ===============================================================

    ! Preset flux in spectral interval
    !$acc parallel
    !$acc loop gang vector collapse(2)
    DO j3 = ki3sc, ki3ec+1
      DO j1 = ki1sc, ki1ec
        zflux (j1,j3) = 0.0_dp
        zfluxd(j1,j3) = 0.0_dp
        zfluxu(j1,j3) = 0.0_dp
      ENDDO
    ENDDO
    !$acc end parallel
 
    ! Upper boundary condition and reference pressure for Rayleigh sc.
    !$acc parallel
    !$acc loop gang vector
    DO j1 = ki1sc, ki1ec
      pflpt(j1) = psct * solant(jspec) * psmu0(j1)
      papra(j1) = zapre(j1) * pqsmu0(j1)
    ENDDO
    !$acc end parallel
 
    ! Optical properties of non-gaseous constituents        
 
    IF (idebug > 10) THEN
       print *,' FESFT    Call to opt_so for jspec: ',jspec
    ENDIF   

    CALL opt_so ( prholwc,pdulwc,prhoiwc,pduiwc,               &
                  paeq1  ,paeq2 ,paeq3  ,paeq4 , paeq5 ,       &
                  pdp    ,papra ,psmu0  ,pqsmu0,               &
                  ki1sd  ,ki1ed ,ki3sd  ,ki3ed , jspec ,       &
                  ki1sc  ,ki1ec ,ki3sc  ,ki3ec ,               &
                  ldebug_opt_so ,jindex ,                      &
                  podac  ,podaf ,podsc  ,podsf , pbsfc ,pbsff ,&
                  pusfc  ,pusff  )

  !----------------------------------------------------------------------------
  ! Section 3.2: Selection of ESFT or FESFT method for interval considered
  !----------------------------------------------------------------------------

    !--------------------------------------------------------------
    IF (nfast(jspec).eq.0) THEN           ! ESFT method
    !--------------------------------------------------------------

      ih2o = ncgas(jspec,1)
      ico2 = ncgas(jspec,2)
      io3  = ncgas(jspec,3)
 
      DO jh2o    = 1, ih2o    ! Loop over H2O coefficients
        DO jco2  = 1, ico2    ! Loop over CO2 coefficients
          DO jo3 = 1, io3     ! Loop over O3  coefficients
 
            zaiprod=coai(jh2o,jspec,1)*coai(jco2,jspec,2)*coai(jo3,jspec,3)
 
            IF (icrf.eq.0) THEN      ! partially cloudy atmosphere
              CALL inv_so (                                                     &
                 pclc   ,pca1   ,pca2  ,pcb1  ,pcb2  ,pcc1  ,pcc2 ,pcd1 ,pcd2 , &
                 pflpt  ,psmu0  ,pqsmu0,palp  ,palso ,                          &
                 pduh2oc,pduh2of,pduco2,pduo3 ,palogp,palogt,                   &
                 podsc  ,podsf  ,podac ,podaf ,pbsfc ,pbsff ,pusfc,pusff,       &
                 jspec  ,jh2o   ,jco2  ,jo3   ,                                 &
                 ki1sd  ,ki1ed  ,ki3sd ,ki3ed ,ki1sc ,ki1ec ,ki3sc,ki3ec,       &
                 ldebug_inv_so  ,jindex,                                        &
                 pflcu  ,pflfu  ,pflcd ,pflfd ,pflcp ,pflfp)
 
            ELSE                     ! cloud-free calculation
              print *,' CRF-Code not implemented yet'
            ENDIF
 
            ! Incrementation of flux in spectral interval
            !$acc parallel
            !$acc loop gang
            DO j3 = ki3sc, ki3ec+1
              !$acc loop vector
              DO j1 = ki1sc, ki1ec
                zflux (j1,j3) = zflux (j1,j3)                           &
                              + zaiprod * (pflfp(j1,j3) + pflcp(j1,j3))
                zfluxd(j1,j3) = zfluxd(j1,j3)                           &
                              + zaiprod * (pflfd(j1,j3) + pflcd(j1,j3))
                zfluxu(j1,j3) = zfluxu(j1,j3)                           &
                              + zaiprod * (pflfu(j1,j3) + pflcu(j1,j3))
              ENDDO
            ENDDO
            !$acc end parallel
 
          ENDDO       ! Loop over O3 absorption coefficients
        ENDDO         ! Loop over CO2 absorption coefficients
      ENDDO           ! Loop over H2O absorption coefficients

    !--------------------------------------------------------------
    ELSE                                  ! FESFT method
    !--------------------------------------------------------------

      igase = 0
      DO jg = 1, 3
         icgas(jg) = 0
         IF (ncgas(jspec,jg).GT.1) THEN
            igase = igase + 1
         ENDIF
      ENDDO

      igasm1    = igase -1
 
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF (igase.le.1) THEN  !(no 'grey' fluxes required)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        !$acc parallel
        !$acc loop gang vector collapse(2)
        DO j3 = ki3sc, ki3ec+1
          DO j1 = ki1sc, ki1ec
            zflux  (j1,j3) = 1.0_dp
            zfluxd (j1,j3) = 1.0_dp
            zfluxu (j1,j3) = 1.0_dp
            zfluxi (j1,j3) = 1.0_dp
            zfluxdi(j1,j3) = 1.0_dp
            zfluxui(j1,j3) = 1.0_dp
          ENDDO
        ENDDO
        !$acc end parallel

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE   ! more than 1 gas --> 'grey' fluxes required
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        IF (icrf.eq.0) THEN      ! partially cloudy atmosphere
          CALL inv_so (                                                    &
             pclc   ,pca1   ,pca2  ,pcb1  ,pcb2  ,pcc1  ,pcc2 ,pcd1 ,pcd2 ,&
             pflpt  ,psmu0  ,pqsmu0,palp  ,palso ,                         &
             pduh2oc,pduh2of,pduco2,pduo3 ,palogp,palogt,                  & 
             podsc  ,podsf  ,podac ,podaf ,pbsfc ,pbsff ,pusfc,pusff,      &
             jspec  ,0      ,0     ,0     ,                                &
             ki1sd  ,ki1ed  ,ki3sd ,ki3ed ,ki1sc ,ki1ec ,ki3sc,ki3ec,      &
             ldebug_inv_so  ,jindex,                                       &
             pflcu  ,pflfu  ,pflcd ,pflfd ,pflcp ,pflfp)
        ELSE                     ! cloud-free calculation
          yerrmsg = ' CRF-Code not implemented yet !!!'
          ierror  = 2616
#ifndef _OPENACC
!Return statement not allowed within data region
          RETURN
#endif
        ENDIF

        ! Storage of 'grey' fluxes and their inverse (**igasm1)
        !$acc parallel
        !$acc loop gang
        DO j3 = ki3sc, ki3ec+1
          !$acc loop vector
          DO j1 = ki1sc, ki1ec
            zfluxi (j1,j3) = 1.0_dp                                  &
                           / MAX (pflfp(j1,j3)+pflcp(j1,j3), zepflx) ** igasm1
            zfluxdi(j1,j3) = 1.0_dp                                  &
                           / MAX (pflfd(j1,j3)+pflcd(j1,j3), zepflx) ** igasm1
            zfluxui(j1,j3) = 1.0_dp                                  &
                           / MAX (pflfu(j1,j3)+pflcu(j1,j3), zepflx) ** igasm1
            zflux  (j1,j3) = pflfp(j1,j3) + pflcp(j1,j3)
            zfluxd (j1,j3) = pflfd(j1,j3) + pflcd(j1,j3)
            zfluxu (j1,j3) = pflfu(j1,j3) + pflcu(j1,j3)
          ENDDO
        ENDDO
        !$acc end parallel

#ifndef _OPENACC
        IF (ldebug_so .AND. jb==jindex) THEN
           print *,' FESFT in debug mode for solar fluxes for js = ', jindex
           print *,' Grey fluxes   '
           DO j3=ki3sc,ki3ec+1
              print *,'par/down/up    : ',zflux (j1b    ,j3) &
                                         ,zfluxd(j1b    ,j3) &
                                         ,zfluxu(j1b    ,j3),j3
           ENDDO
        ENDIF 
#endif
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ENDIF    ! No.of relevant gases
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      igasz = 0
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO jg = 3, 1, -1   !     Loop over gases
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
      IF (ncgas(jspec,jg).GT.1) THEN   ! include gas only, if necessary
   
        igasz = igasz + 1
 
        DO jjg = 1,3
           icgas(jjg) = 0
        ENDDO 
 
        ! Initialize 'gaseous' fluxes
        !$acc parallel
        !$acc loop gang vector collapse(2)
        DO j3 = ki3sc, ki3ec+1
          DO j1 = ki1sc, ki1ec
            zfgas (j1,j3) = 0.0_dp
            zfgasd(j1,j3) = 0.0_dp
            zfgasu(j1,j3) = 0.0_dp
          ENDDO
        ENDDO
        !$acc end parallel
 
        ! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
        icc = ncgas(jspec,jg)
        DO jc = icc,1,-1        ! Loop over absorption coefficients
        ! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          zcoai = coai(jc,jspec,jg)
          zcobi = cobi(jc,jspec,jg)
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          IF ( ((zcoai.GE.zepai).AND.(zcobi.GT.0.0_dp)) .OR. (igase.EQ.1) ) THEN 
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! Solve linear system, if necessary
            icgas(jg) = jc
            IF (icrf.eq.0) THEN      ! partially cloudy atmosphere
              CALL inv_so (                                                     &
                 pclc   ,pca1   ,pca2  ,pcb1  ,pcb2  ,pcc1  ,pcc2 ,pcd1 ,pcd2 , &
                 pflpt  ,psmu0  ,pqsmu0,palp  ,palso ,                          &
                 pduh2oc,pduh2of,pduco2,pduo3 ,palogp,palogt,                   &
                 podsc  ,podsf  ,podac ,podaf ,pbsfc ,pbsff ,pusfc,pusff,       &
                 jspec  ,icgas(1),icgas(2),icgas(3),                            &
                 ki1sd  ,ki1ed  ,ki3sd ,ki3ed ,ki1sc ,ki1ec ,ki3sc,ki3ec,       &
                 ldebug_inv_so  ,jindex,                                        &
                 pflcu  ,pflfu  ,pflcd ,pflfd ,pflcp ,pflfp)
            ELSE                     ! cloud-free clculations
              yerrmsg = ' CRF-Code not implemented yet !!!'
              ierror  = 2616
#ifndef _OPENACC
              !Return statement not allowed within data region
              RETURN
#endif
            ENDIF

          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ELSE               ! use 'grey' fluxes directly
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            !$acc parallel
            !$acc loop gang
            DO j3 = ki3sc, ki3ec+1
              !$acc loop vector
              DO j1 = ki1sc, ki1ec
                pflfp(j1,j3) = zflux (j1,j3)
                pflcp(j1,j3) = 0.0_dp
                pflfd(j1,j3) = zfluxd(j1,j3)
                pflcd(j1,j3) = 0.0_dp
                pflfu(j1,j3) = zfluxu(j1,j3)
                pflcu(j1,j3) = 0.0_dp
              ENDDO
            ENDDO
            !$acc end parallel

          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ENDIF            ! Necessity to calculate fluxes
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
          !$acc parallel
          !$acc loop gang
          DO j3 = ki3sc, ki3ec+1
            !$acc loop vector
            DO j1 = ki1sc, ki1ec
              zfgas (j1,j3) = zfgas (j1,j3)                          &
                            + zcoai * (pflfp(j1,j3) + pflcp(j1,j3))
              zfgasd(j1,j3) = zfgasd(j1,j3)                          &
                            + zcoai * (pflfd(j1,j3) + pflcd(j1,j3))
              zfgasu(j1,j3) = zfgasu(j1,j3)                          &
                            + zcoai * (pflfu(j1,j3) + pflcu(j1,j3))
            ENDDO
          ENDDO
          !$acc end parallel

#ifndef _OPENACC
          !?????????????????????????????????????????????????????????????????
          IF (ldebug_so .AND. jb==jindex) THEN
             print *,' FESFT in debug mode for solar fluxes for js = ', jindex
             print *,' only one interval/coefficient considered '
             print *,' zcoai = ',zcoai
             DO j3=ki3sc,ki3ec+1
                print *,'zfgas(j1b,jb,j3): ',zfgas(j1b    ,j3)
             ENDDO 
             EXIT solar_spectral_loop
          ENDIF 
          !?????????????????????????????????????????????????????????????????
#endif
        ! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
        ENDDO      ! Loop over absorption coefficients
        ! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
 
        ! Combination of inverse of 'grey' fluxes and actual gaseous flux
        !$acc parallel
        !$acc loop gang vector collapse(2)
        DO j3 = ki3sc, ki3ec+1
          DO j1 = ki1sc, ki1ec
            zfluxi (j1,j3) = zfluxi (j1,j3) * zfgas (j1,j3)
            zfluxdi(j1,j3) = zfluxdi(j1,j3) * zfgasd(j1,j3)
            zfluxui(j1,j3) = zfluxui(j1,j3) * zfgasu(j1,j3)
          ENDDO
        ENDDO
        !$acc end parallel
 
        IF (igasz.eq.igasm1) THEN    !     Avoid unphysical pseudo-transmission
          !$acc parallel
          !$acc loop gang vector collapse(2)
          DO j3 = ki3sc, ki3ec+1
            DO j1 = ki1sc, ki1ec
              zfluxi (j1,j3) = MIN( 1.0_dp, MAX( 0.0_dp, zfluxi (j1,j3)) )
              zfluxdi(j1,j3) = MIN( 1.0_dp, MAX( 0.0_dp, zfluxdi(j1,j3)) )
              zfluxui(j1,j3) = MIN( 1.0_dp, MAX( 0.0_dp, zfluxui(j1,j3)) )
            ENDDO
          ENDDO
          !$acc end parallel
        END IF

      END IF                   ! Test, whether gas needs to be included
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END DO          ! End of loop over gases
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
      ! Store FESFT result in zflux
      !$acc parallel
      !$acc loop gang vector collapse(2)
      DO j3 = ki3sc, ki3ec+1
        DO j1 = ki1sc, ki1ec
          zflux (j1,j3) = zfluxi (j1,j3)
          zfluxd(j1,j3) = zfluxdi(j1,j3)
          zfluxu(j1,j3) = zfluxui(j1,j3)
        ENDDO
      ENDDO
      !$acc end parallel
 
    !--------------------------------------------------------------
    END IF                             ! ESFT/FESFT-Selection 
    !--------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Section 3.3: Compute corrected fluxes, if lradtopo
  !----------------------------------------------------------------------------

    IF (lradtopo) THEN
#ifndef _OPENACC
      IF (idebug > 15 .AND. jb==jindex) THEN
        WRITE (*,'(A60,2F16.6)')                                             &
                  'FESFT: solar fluxes before and after, skyview, fcor',     &
                          pskyview(j1b    ), pfcor(j1b    )
        WRITE (*,'(A32,2I4,3F16.6)')  'FESFT: zfluxd, zflux, zfluxu',        &
                j1b, jb, zfluxd(j1b    ,ki3ec+1), zflux(j1b    ,ki3ec+1),   &
                          zfluxu(j1b    ,ki3ec+1)
        zalbedo = zfluxu(j1b    ,ki3ec+1) /                                  &
                          (zflux(j1b    ,ki3ec+1)+zfluxd(j1b    ,ki3ec+1))
        WRITE (*,'(A32,2F16.6)') 'albedo, palso', zalbedo, palso(j1b    )
      ENDIF
#endif

      !$acc parallel
      !$acc loop gang vector
      DO j1 = ki1sc, ki1ec

        zalbedo = zfluxu(j1,ki3ec+1) /                                    &
                      (zflux(j1,ki3ec+1)+zfluxd(j1,ki3ec+1))

        ! direct down corrected
        zflux_c (j1,ki3ec+1) = pfcor(j1) * zfluxi(j1,ki3ec+1)

        ! diffuse down corrected
        zfluxd_c(j1,ki3ec+1) =                                            &
                  zfluxdi(j1,ki3ec+1) *         pskyview(j1)              &
                + zfluxui(j1,ki3ec+1) * (1.0_dp-pskyview(j1))

        ! diffuse up adapted to new other components
        zfluxu_c(j1,ki3ec+1) =                                            &
                 (zflux_c(j1,ki3ec+1) + zfluxd_c(j1,ki3ec+1)) * zalbedo

      ENDDO
      !$acc end parallel

#ifndef _OPENACC
      IF (idebug > 15 .AND. jb==jindex) THEN
        WRITE (*,'(A42,2I4,3F16.6)')                                         &
                      'FESFT corrected: zfluxd, zflux, zfluxu', j1b, jb,    &
                      zfluxd_c(j1b    ,ki3ec+1), zflux_c(j1b    ,ki3ec+1),   &
                      zfluxu_c(j1b    ,ki3ec+1)
      ENDIF
#endif
    ENDIF

  !----------------------------------------------------------------------------
  ! Section 3.4: Addition of flux for spectral interval to total solar flux
  !----------------------------------------------------------------------------
 
    IF (icrf == 0) THEN     ! Add flux at all levels

      !$acc parallel
      !$acc loop gang
      DO j3 = ki3sc, ki3ec+1
        !$acc loop vector
        DO j1 = ki1sc, ki1ec
          pfls(j1,j3) = pfls (j1,j3)                                 &
                      + zflux(j1,j3) + zfluxd(j1,j3) - zfluxu(j1,j3)
          ! for the Climate-LM Version
          IF ( (.NOT. lradtopo) .OR. (pfcor(j1) /= 0.0_dp) ) THEN
            pflsdir(j1,j3) = pflsdir (j1,j3) + zflux(j1,j3)
          ! the else part just lets pflsdir untouched
          ENDIF
        ENDDO
      ENDDO
      !$acc end parallel

      ! Store individual components of solar flux at surface
      IF (lradtopo) THEN
        !$acc parallel
        !$acc loop gang vector
        DO j1 = ki1sc, ki1ec
          pflsu (j1) = pflsu (j1) + zfluxu_c(j1,ki3ec+1)
          pflsd (j1) = pflsd (j1) + zfluxd_c(j1,ki3ec+1)
          pflsp (j1) = pflsp (j1) + zflux_c (j1,ki3ec+1)
          pfls_s(j1) = pfls_s(j1) + zflux_c (j1,ki3ec+1)            &
                     + zfluxd_c(j1,ki3ec+1) - zfluxu_c(j1,ki3ec+1)
        ENDDO
        !$acc end parallel
      ELSE
        !$acc parallel
        !$acc loop gang vector
        DO j1 = ki1sc, ki1ec
          pflsu (j1) = pflsu (j1) + zfluxu (j1,ki3ec+1)
          pflsd (j1) = pflsd (j1) + zfluxd (j1,ki3ec+1)
          pflsp (j1) = pflsp (j1) + zflux  (j1,ki3ec+1)
          pfls_s(j1) = pfls_s(j1) + zflux  (j1,ki3ec+1)             &
                     + zfluxd (j1,ki3ec+1) - zfluxu (j1,ki3ec+1)
        ENDDO
        !$acc end parallel
      ENDIF

#ifndef _OPENACC
      IF (idebug > 15 .AND. jb==jindex) THEN
        WRITE (*,'(A40,3F16.6)') 'FESFT: diff_up, diff_down, dir_down',      &
            pflsu(j1b    ), pflsd(j1b    ), pflsp(j1b    )
      ENDIF
#endif

      IF (jspec == 3) THEN   ! Photosynthetic active radiation
        IF (lradtopo) THEN   ! T.R.
          !$acc parallel
          !$acc loop gang vector
          DO j1 = ki1sc, ki1ec
            pflsu_par(j1) = pflsu_par (j1) + zfluxu_c(j1,ki3ec+1)
            pflsd_par(j1) = pflsd_par (j1) + zfluxd_c(j1,ki3ec+1)
            pflsp_par(j1) = pflsp_par (j1) + zflux_c (j1,ki3ec+1)
            pflpar   (j1) = pflpar    (j1) + zflux_c (j1,ki3ec+1)       &
                          + zfluxd_c(j1,ki3ec+1) - zfluxu_c(j1,ki3ec+1)
          ENDDO
          !$acc end parallel
        ELSE
          !$acc parallel
          !$acc loop gang vector
          DO j1 = ki1sc, ki1ec
            pflsu_par(j1) = pflsu_par (j1) + zfluxu(j1,ki3ec+1)
            pflsd_par(j1) = pflsd_par (j1) + zfluxd(j1,ki3ec+1)
            pflsp_par(j1) = pflsp_par (j1) + zflux (j1,ki3ec+1)
            pflpar   (j1) = pflpar    (j1) + zflux (j1,ki3ec+1)       &
                          + zfluxd (j1,ki3ec+1) - zfluxu (j1,ki3ec+1)
          ENDDO
          !$acc end parallel
        ENDIF
      ENDIF !(jspec == 3)

    ENDIF
 
#ifdef COSMOART
  IF(l_cosmo_art .AND. lgas) THEN
    IF (jspec == 3) THEN
      IF (lradtopo) THEN
        !$acc parallel
        !$acc loop gang
        DO j3 = ki3sc, ki3ec+1
          !$acc loop vector
          DO j1 = ki1sc, ki1ec
            Edir_b (j1,j3) =  zflux_c (j1,j3)
            Edown_b(j1,j3) =  zfluxd_c(j1,j3)
            Eup_b  (j1,j3) =  zfluxu_c(j1,j3)
          ENDDO
        ENDDO
        !$acc end parallel
      ELSE
        !$acc parallel
        !$acc loop gang
        DO j3 = ki3sc, ki3ec+1
          !$acc loop vector
          DO j1 = ki1sc, ki1ec
            Edir_b (j1,j3) =  zflux (j1,j3)
            Edown_b(j1,j3) =  zfluxd(j1,j3)
            Eup_b  (j1,j3) =  zfluxu(j1,j3)
          ENDDO
        ENDDO
        !$acc end parallel
      ENDIF
    ENDIF
  ENDIF
#endif

  ! End of solar spectral loop
  ! ===============================================================
  ENDDO solar_spectral_loop 
  ! ===============================================================
 
ENDIF     ! Test, whether solar calculation or not
 
!------------------------------------------------------------------------------
! Section 4: Repeat calculations for cloud-free fluxes
!------------------------------------------------------------------------------

  ! Repeat calculations for cloud-free fluxes if switch for CRF
  ! is set to .true. and cloud-free fluxes have not yet been
  ! computed

  IF (.NOT. lcrf) THEN
!T.R.: pfltf & pflsf removed
!    DO j2 = ki2sc, ki2ec
!      DO j1 = ki1sc, ki1ec
!        pflsf(j1   ,1)  = pfls(j1   ,ki3sc  )
!        pflsf(j1   ,2)  = pfls(j1   ,ki3ec+1)
!        pfltf(j1   ,1)  = pflt(j1   ,ki3sc  )
!        pfltf(j1   ,2)  = pflt(j1   ,ki3ec+1)
!      ENDDO
!    ENDDO
  ELSE IF (icrf.eq.0) THEN  ! Branch to cloud-free calculations only once
    icrf = 1     
    GO TO 1
  ENDIF
 
!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

  !$acc end data

END SUBROUTINE fesft_dp

!==============================================================================
!==============================================================================
!+ Module procedure in "radiation_rg"
!------------------------------------------------------------------------------

SUBROUTINE coe_th (                                                    &
       pduh2oc,pduh2of,pduco2 ,pduo3  ,palogp ,palogt ,                &
       podsc  ,podsf  ,podac  ,podaf  ,pbsfc  ,pbsff  ,                &
       ki3    ,kspec  ,kh2o   ,kco2   ,ko3    ,                        &
       ki1sd  ,ki1ed  ,ki3sd  ,ki3ed  ,ki1sc  ,ki1ec  ,                &
       ldebug ,jindex ,                                                &
       pa1c   ,pa1f   ,pa2c   ,pa2f   ,pa3c   ,pa3f)

!------------------------------------------------------------------------------
!
! Description:
!
!   The module procedure coe_th calculates the optical effects of atmospheric 
!   layers on thermal radiation based on basic optical properties of non-gaseous 
!   constituents and gaseous absorption coefficients selected through the 
!   corresponding control variables in the argument list.
!   This routine computes layer effects (transmissivity, reflectivity
!   and emmisivity) in the thermal part of the radiative spectrum
!   both for the cloud-free and the cloudy part of a model layer.
!   The calculation is based on the implicit delt-two-stream equations
!   (cf. Ritter and Geleyn, 1992) and uses basic optical properties
!   (i.e. absorption and scattering optical depth and backscattered
!   fraction for non-gaseous atmospheric constituents as well as 
!   gaseous absorption properties) as input. 
!
! Method:
!
! - addition of individual gaseous absorption effects to the optical
!   properties of the non-gaseous constituents
! - determination of layer effects (cf. Zdunkowski et al., 1982, 1986
!   and Ritter and Geleyn, 1992)
!     
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

! Input data
! ----------
  INTEGER, INTENT (IN) ::  &
     ki1sd,       & ! start index for first  array dimension
     ki1ed,       & ! end   index for first  array dimension
     ki3sd,       & ! start index for third  array dimension
     ki3ed,       & ! end   index for third  array dimension

   ! and the same for the computations
     ki1sc,       & ! start index for first  array computation
     ki1ec,       & ! end   index for first  array computation
     ki3  ,       & ! vertical layer considered    
     kspec,       & ! spectral interval considered
     kh2o ,       & ! table index for h2o absorption properties
     kco2 ,       & ! table index for co2 absorption properties
     ko3  ,       & ! table index for o3  absorption properties
     jindex         ! index for j-loop

  LOGICAL                 , INTENT (IN) ::  &
     ldebug         ! debug control switch       

  REAL    (KIND=dp)       , INTENT (IN) ::  &

     ! optical relevant gas quantities (Pa)
     pduh2oc(ki1sd:ki1ed,ki3sd:ki3ed), & ! h2o inside cloud
     pduh2of(ki1sd:ki1ed,ki3sd:ki3ed), & ! h2o out of cloud
     pduco2 (ki1sd:ki1ed,ki3sd:ki3ed), & ! co2 content 
     pduo3  (ki1sd:ki1ed,ki3sd:ki3ed), & ! o3  content 

     ! Logarithm of layer mean temperature and pressure
     palogt (ki1sd:ki1ed,ki3sd:ki3ed), & ! ln T
     palogp (ki1sd:ki1ed,ki3sd:ki3ed), & ! ln p

     ! Optical properties of non-gaseous constituents (..c=cloudy; ..f=free)  
     podsc  (ki1sd:ki1ed,ki3sd:ki3ed), & ! 
     podsf  (ki1sd:ki1ed,ki3sd:ki3ed), & ! 
     podac  (ki1sd:ki1ed,ki3sd:ki3ed), & ! 
     podaf  (ki1sd:ki1ed,ki3sd:ki3ed), & ! 
     pbsfc  (ki1sd:ki1ed,ki3sd:ki3ed), & ! 
     pbsff  (ki1sd:ki1ed,ki3sd:ki3ed)    ! 

! Output data
! -----------
  REAL    (KIND=dp)       , INTENT (OUT) ::  &
     pa1c  (ki1sd:ki1ed)             , & ! transmissivity in cloud   
     pa1f  (ki1sd:ki1ed)             , & ! transmissivity cloud-free  
     pa2c  (ki1sd:ki1ed)             , & ! reflectivity in cloud    
     pa2f  (ki1sd:ki1ed)             , & ! reflectivity cloud-free      
     pa3c  (ki1sd:ki1ed)             , & ! emissivity in cloud    
     pa3f  (ki1sd:ki1ed)                 ! emissivity cloud-free       

! Local parameters: 
! ----------------
  REAL    (KIND=dp)       , PARAMETER ::  &
     zargli  = 80.0_dp     , &  ! argument limit for EXP 
     ztsec   = 1.0E-35_dp  , &  ! (=exp(-zargli) avoids ALOG(0.0)
     zodmax  = 1.0E+6_dp   , &  ! maximum allowed optical depth
     zudiff  = 2.0_dp      , &  ! Diffusivity factors for gases and other constituents
     zangfa  = 1.648721271_dp   ! exp(0.5)

! Local scalars:
! -------------
  INTEGER  ::  &
    j1,j3                    ! loop indices over spatial dimensions

  REAL    (KIND=dp)        ::  &
    zeps, ztau, zrho, zodgf, zodgc, zod1, zod2
 
! End of header
!==============================================================================
 
!------------------------------------------------------------------------------
! Begin Subroutine coe_th              
!------------------------------------------------------------------------------

  j3     = ki3

  IF (ldebug .AND. jb==jindex) THEN
    DO j1 = ki1sc, ki1ec
      IF (j1 == ki1sc) THEN
        print *,'**** coe_th ******************************'
        print *,'**** debug point j1 = ',j1, ' in level j3 = ', j3
        print *,'**** coe_th kspec=',kspec
        print *,'**** coe_th j3   =',j3   
        print *,'**** coe_th kh2o =',kh2o 
        print *,'**** coe_th kco2 =',kco2 
        print *,'**** coe_th ko3  =',ko3  
        print *,'**** pduh2of(j1,j3)=',pduh2of(j1,j3)
        print *,'**** pduh2oc(j1,j3)=',pduh2oc(j1,j3)
        print *,'**** pduco2 (j1,j3)=',pduco2 (j1,j3)
        print *,'**** pduo3  (j1,j3)=',pduo3  (j1,j3)
        print *,'**** palogp (j1,j3)=',palogp (j1,j3)
        print *,'**** palogt (j1,j3)=',palogt (j1,j3)
        if (kh2o > 0) print *,'**** cobi (kh2o,kspec,1)    =',cobi (kh2o,kspec,1)
        if (kco2 > 0) print *,'**** cobi (kco2,kspec,2)    =',cobi (kco2,kspec,2)
        if (ko3  > 0) print *,'**** cobi (ko3 ,kspec,3)    =',cobi (ko3 ,kspec,3)
        if (kh2o > 0) print *,'**** coali(kh2o,kspec,1)    =',coali (kh2o,kspec,1)
        if (kco2 > 0) print *,'**** coali(kco2,kspec,2)    =',coali (kco2,kspec,2)
        if (ko3  > 0) print *,'**** coali(ko3 ,kspec,3)    =',coali (ko3 ,kspec,3)
        if (kh2o > 0) print *,'**** cobti(kh2o,kspec,1)    =',cobti(kh2o,kspec,1)
        if (kco2 > 0) print *,'**** cobti(kco2,kspec,2)    =',cobti(kco2,kspec,2)
        if (ko3  > 0) print *,'**** cobti(ko3 ,kspec,3)    =',cobti(ko3 ,kspec,3)
      ENDIF
    ENDDO
  ENDIF 

  ! Optical depth of gases

  DO j1 = ki1sc, ki1ec
    zodgf = 0.0_dp     ! Initialisation
 
    IF (kco2.ne.0) then     ! Include CO2 contribution
      zodgf = zodgf + pduco2(j1,j3) * (cobi(kco2,kspec,2)        &
                  * EXP ( coali(kco2,kspec,2) * palogp(j1,j3)    &
                         -cobti(kco2,kspec,2) * palogt(j1,j3)))
    ENDIF                  ! CO2

    IF (ko3 /= 0) THEN  ! Include O3 contribution
      zodgf = zodgf + pduo3 (j1,j3) * (cobi(ko3 ,kspec,3)*       &
                  EXP ( coali(ko3 ,kspec,3) * palogp(j1,j3)      &
                       -cobti(ko3 ,kspec,3) * palogt(j1,j3)))
    ENDIF

    ! Cloudy = cloud free for CO2 and O3 :
    zodgc = zodgf

    IF (kh2o /= 0) THEN  ! Include H2O contribution
      zodgf = zodgf + pduh2of(j1,j3)* (cobi(kh2o,kspec,1)*       &
                  EXP ( coali(kh2o,kspec,1) * palogp(j1,j3)      &
                       -cobti(kh2o,kspec,1) * palogt(j1,j3)))
      zodgc = zodgc + pduh2oc(j1,j3)* (cobi(kh2o,kspec,1)*       &
                  EXP ( coali(kh2o,kspec,1) * palogp(j1,j3)      &
                       -cobti(kh2o,kspec,1) * palogt(j1,j3)))
    ENDIF

    zodgf = MIN (zodgf, zodmax)
    zodgc = MIN (zodgc, zodmax)

    ! Pseudo-optical depth in cloud-free part of layer

    zod2 = zudiff * pbsff(j1,j3) * podsf(j1,j3)
    zod1 = zod2 + zudiff * podaf(j1,j3)
    zod1 = zod1 + zangfa * zodgf

    ! Layer coefficients in cloud-free part of layer

    zeps=SQRT(zod1*zod1-zod2*zod2)
! for a better vectorization
!   IF (zeps.LT.zargli) THEN
!     ztau = EXP  (-zeps)
      ztau = EXP  (-MIN(zeps,zargli))
!   ELSE
!     ztau = ztsec
!   END IF
    zrho = zod2/(zod1+zeps)
    pa1f(j1)=ztau*(1.0_dp-(zrho**2))*(1.0_dp/(1.0_dp-(zrho**2)*(ztau**2)))
    pa2f(j1)=zrho*(1.0_dp-(ztau**2))*(1.0_dp/(1.0_dp-(zrho**2)*(ztau**2)))
    pa3f(j1)=(1.0_dp-pa1f(j1)+pa2f(j1))/(zod1+zod2)

    ! Pseudo-optical depth in cloudy part of layer
    zod2 = zudiff * pbsfc(j1,j3) * podsc(j1,j3)
    zod1 = zod2 + zudiff * podac(j1,j3)
    zod1 = zod1 + zangfa * zodgc

    ! Layer coefficients in cloudy part of layer

    zeps=SQRT(zod1*zod1-zod2*zod2)
! for a better vectorization
!   IF (zeps.LT.zargli) THEN
!     ztau = EXP  (-zeps)
      ztau = EXP  (-MIN(zeps,zargli))
!   ELSE
!     ztau = ztsec
!   END IF
    zrho = zod2/(zod1+zeps)
    pa1c(j1)=ztau*(1.0_dp-(zrho**2))*(1.0_dp/(1.0_dp-(zrho**2)*(ztau**2)))
    pa2c(j1)=zrho*(1.0_dp-(ztau**2))*(1.0_dp/(1.0_dp-(zrho**2)*(ztau**2)))
    pa3c(j1)=(1.0_dp-pa1c(j1)+pa2c(j1))/(zod1+zod2)

  ENDDO
 
!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE coe_th

!==============================================================================
!==============================================================================
!+ Module procedure in "radiation_rg"
!------------------------------------------------------------------------------

SUBROUTINE coe_so (                                                     &
       pduh2oc,pduh2of,pduco2 ,pduo3  ,palogp ,palogt ,                 &
       podsc  ,podsf  ,podac  ,podaf  ,pbsfc  ,pbsff  ,pusfc ,pusff ,   &
       psmu0  ,pqsmu0 ,                                                 &
       ki3    ,kspec  ,kh2o   ,kco2   ,ko3    ,                         &
       ki1sd  ,ki1ed  ,ki3sd  ,ki3ed  ,ki1sc  ,ki1ec  ,                 &
       ldebug ,jindex ,                                                 &
       pa1c   ,pa1f   ,pa2c   ,pa2f   ,pa3c   ,pa3f ,                   &
       pa4c   ,pa4f   ,pa5c   ,pa5f )

!------------------------------------------------------------------------------
!
! Description:
!
!   The module procedure coe_so calculates the optical effects of atmospheric
!   layers on solar radiation based on basic optical properties of non-gaseous
!   constituents and gaseous absorption coefficients selected through the
!   corresponding control variables.
!   This routine computes layer effects (transmissivity, reflectivity)
!   for diffuse and direct solar radiation both for the cloudy and the
!   cloud-free part of a model layer.
!   The calculation is based on the implicit delt-two-stream equations
!   (cf. Ritter and Geleyn, 1992) and uses basic optical properties
!   (i.e. absorption and scattering optical depth, backscattered and
!   upscattered fraction for non-gaseous atmospheric constituents and
!   gaseous absorption properties) as input.
!
! Method:
!
! - addition of individual gaseous absorption effects to the optical
!   properties of the non-gaseous constituents
!   (optical depth multiplied by alpha1 to alpha4)
!
! - determination of layer effects (cf. Zdunkowski et al., 1982, 1986
!   and Ritter and Geleyn, 1992)
! - the resonance case for those effects related to the direct solar
!   radiation is avoided by a small displacement of the local inverse
!   of the cosine of the zenith angle (if necessary)
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

! Input data
! ----------
  INTEGER, INTENT (IN) ::  &
     ki1sd,       & ! start index for first  array dimension
     ki1ed,       & ! end   index for first  array dimension
     ki3sd,       & ! start index for third  array dimension
     ki3ed,       & ! end   index for third  array dimension

   ! and the same for the computations
     ki1sc,       & ! start index for first  array computation
     ki1ec,       & ! end   index for first  array computation
     ki3  ,       & ! vertical layer considered
     kspec,       & ! spectral interval considered
     kh2o ,       & ! table index for h2o absorption properties
     kco2 ,       & ! table index for co2 absorption properties
     ko3  ,       & ! table index for o3  absorption properties
     jindex         ! index for j-loop

  LOGICAL                 , INTENT (IN) ::  &
     ldebug         ! debug control switch

  REAL    (KIND=dp)       , INTENT (IN) ::  &

     ! opticall relevant gas quantities (Pa)
     pduh2oc(ki1sd:ki1ed,ki3sd:ki3ed), & ! h2o inside cloud
     pduh2of(ki1sd:ki1ed,ki3sd:ki3ed), & ! h2o out of cloud
     pduco2 (ki1sd:ki1ed,ki3sd:ki3ed), & ! co2 content
     pduo3  (ki1sd:ki1ed,ki3sd:ki3ed), & ! o3  content

     ! Logarithm of layer mean temperature and pressure
     palogt (ki1sd:ki1ed,ki3sd:ki3ed), & ! ln T
     palogp (ki1sd:ki1ed,ki3sd:ki3ed), & ! ln p

     ! Optical properties of non-gaseous constituents (..c=cloudy; ..f=free)
     podsc  (ki1sd:ki1ed,ki3sd:ki3ed), & !
     podsf  (ki1sd:ki1ed,ki3sd:ki3ed), & !
     podac  (ki1sd:ki1ed,ki3sd:ki3ed), & !
     podaf  (ki1sd:ki1ed,ki3sd:ki3ed), & !
     pbsfc  (ki1sd:ki1ed,ki3sd:ki3ed), & !
     pbsff  (ki1sd:ki1ed,ki3sd:ki3ed), & !
     pusfc  (ki1sd:ki1ed,ki3sd:ki3ed), & !
     pusff  (ki1sd:ki1ed,ki3sd:ki3ed), & !

     psmu0  (ki1sd:ki1ed)            , & ! cosine of zenith angle
     pqsmu0 (ki1sd:ki1ed)                ! inverse of cosine ...

! Output data
! -----------
  REAL    (KIND=dp)       , INTENT (OUT) ::  &
     pa1c  (ki1sd:ki1ed)             , & ! direct radiation transmis-
     pa1f  (ki1sd:ki1ed)             , & ! sivity cloudy/cloud-free
     pa2c  (ki1sd:ki1ed)             , & ! direct radition downward
     pa2f  (ki1sd:ki1ed)             , & ! scattering cloudy/cloud-free
     pa3c  (ki1sd:ki1ed)             , & ! direct radiation upward
     pa3f  (ki1sd:ki1ed)             , & ! scattering cloudy/cloud-free
     pa4c  (ki1sd:ki1ed)             , & ! diffuse flux transmissivity
     pa4f  (ki1sd:ki1ed)             , & ! cloudy/cloud-free
     pa5c  (ki1sd:ki1ed)             , & ! diffuse flux reflectivity
     pa5f  (ki1sd:ki1ed)                 ! cloudy/cloud-free

! Local parameters:
! ----------------
  REAL    (KIND=dp)       , PARAMETER ::  &
     zargli  = 80.0_dp     , &  ! argument limit for EXP
     ztsec   = 1.0E-35_dp  , &  ! (=exp(-zargli) avoids ALOG(0.0)
     zodmax  = 1.0E+6_dp   , &  ! maximum allowed optical depth
     zepres  = 1.0E-7_dp   , &  ! for resonance case avoidance
                                ! 32bit-accuracy (1.E-14 for 64bit machine)
     zudiff  = 2.0_dp      , &  ! Diffusivity factors for gases and other constituents
     zangfa  = 1.648721271_dp   ! exp(0.5)

! Local scalars:
! -------------
  INTEGER  ::  &
    j1, j3                    ! loop indices over spatial dimensions

  REAL    (KIND=dp)        ::  &
    zeps,                      & !
    ze,zm,zg1,zg2,ze1mwf,zmu0if  !

  REAL    (KIND=dp)        ::  &
     zodgf, zodgc, zod1, zod2, zod3, zod4, zod5
 
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine coe_so
!------------------------------------------------------------------------------

  j3     = ki3

  IF (ldebug .AND. jb==jindex) THEN
     print *,'**** coe_so ******************************'
     print *,'**** debug point index : ',j1b,jb, ' in level ', j3
     print *,'**** coe_so kspec=',kspec
     print *,'**** coe_so j3   =',j3
     print *,'**** coe_so kh2o =',kh2o
     print *,'**** coe_so kco2 =',kco2
     print *,'**** coe_so ko3  =',ko3
     print *,'**** pduh2of(j1b,jb,j3)=',pduh2of(j1b,j3)
     print *,'**** pduh2oc(j1b,jb,j3)=',pduh2oc(j1b,j3)
     print *,'**** pduco2 (j1b,jb,j3)=',pduco2 (j1b,j3)
     print *,'**** pduo3  (j1b,jb,j3)=',pduo3  (j1b,j3)
     print *,'**** psmu0  (j1b,jb)   =',psmu0  (j1b)
  ENDIF

  ! Optical depth of gases

  DO j1 = ki1sc, ki1ec
    zodgf = 0.0_dp          ! Initialisation

    IF (kco2 /= 0) THEN     ! Include CO2 contribution
      zodgf = zodgf + pduco2(j1,j3)* (cobi(kco2,kspec,2)*       &
                  EXP ( coali(kco2,kspec,2) * palogp(j1,j3)     &
                       -cobti(kco2,kspec,2) * palogt(j1,j3)))
    ENDIF                  ! CO2

    IF (ko3 /= 0) THEN     ! Include O3 contribution
      zodgf = zodgf + pduo3 (j1,j3)* (cobi(ko3 ,kspec,3)*       &
                  EXP ( coali(ko3 ,kspec,3) * palogp(j1,j3)     &
                       -cobti(ko3 ,kspec,3) * palogt(j1,j3)))
    ENDIF

    ! Cloudy = cloud free for CO2 and O3 :
    zodgc = zodgf

    IF (kh2o /= 0) THEN    ! Include H2O contribution
      zodgf = zodgf + pduh2of(j1,j3)* (cobi(kh2o,kspec,1)*       &
                  EXP ( coali(kh2o,kspec,1) * palogp(j1,j3)      &
                       -cobti(kh2o,kspec,1) * palogt(j1,j3)))
      zodgc = zodgc + pduh2oc(j1,j3)* (cobi(kh2o,kspec,1)*       &
                  EXP ( coali(kh2o,kspec,1) * palogp(j1,j3)      &
                       -cobti(kh2o,kspec,1) * palogt(j1,j3)))
    ENDIF

    zodgf = MIN (zodgf, zodmax)
    zodgc = MIN (zodgc, zodmax)

    ! Pseudo-optical depth in cloud-free part of layer

    zod2 = zudiff * pbsff(j1,j3) * podsf(j1,j3)
    zod1 = zod2 + zudiff * podaf(j1,j3)
    zod3 = pusff(j1,j3) * podsf(j1,j3)
    zod4 = podsf(j1,j3) - zod3
    zod5 = podsf(j1,j3) + podaf(j1,j3)
    zod1 = zod1 + zangfa * zodgf
    zod5 = zod5 + zodgf

    ! Layer coefficients in cloud-free part of layer

    zeps=SQRT(zod1*zod1-zod2*zod2)
! for a better vectorization
!   IF (zeps.LT.zargli) THEN
!     ze = EXP  (-zeps)
      ze = EXP  (-MIN(zeps,zargli))
!   ELSE
!     ze = ztsec
!   END IF
    zm = zod2/(zod1+zeps)
    pa4f(j1)=ze*(1.0_dp-(zm**2))*(1.0_dp/(1.0_dp-(zm**2)*(ze**2)))
    pa5f(j1)=zm*(1.0_dp-(ze**2))*(1.0_dp/(1.0_dp-(zm**2)*(ze**2)))

    ze1mwf = zeps / zod5
    zmu0if = ze1mwf + SIGN ( MAX(ABS(pqsmu0(j1)-ze1mwf),zepres) &
                            ,(pqsmu0(j1)-ze1mwf) )
    zod3 = zod3 * zmu0if
    zod4 = zod4 * zmu0if
    zod5 = zod5 * zmu0if
! for a better vectorization
!   IF (zod5.LT.zargli) THEN
!     pa1f(j1) = EXP  (-zod5)
      pa1f(j1) = EXP  (-MIN(zod5,zargli))
!   ELSE
!     pa1f(j1) = ztsec
!   END IF
    zg1 = ( zod3*(zod5-zod1) -zod2*zod4) /(zod5*zod5 - zeps*zeps)
    zg2 =-( zod4*(zod5+zod1) +zod2*zod3) /(zod5*zod5 - zeps*zeps)
    pa2f(j1) = zg2*(pa1f(j1)-pa4f(j1)) -zg1*pa5f(j1)*pa1f(j1)
    pa3f(j1) = zg1*(1.0_dp-pa4f(j1)*pa1f(j1)) -zg2*pa5f(j1)

    ! Pseudo-optical depth in cloudy part of layer

    zod2 = zudiff * pbsfc(j1,j3) * podsc(j1,j3)
    zod1 = zod2 + zudiff * podac(j1,j3)
    zod3 = pusfc(j1,j3) * podsc(j1,j3)
    zod4 = podsc(j1,j3) - zod3
    zod5 = podsc(j1,j3) + podac(j1,j3)
    zod1 = zod1 + zangfa * zodgc
    zod5 = zod5 + zodgc

    ! Layer coefficients in cloudy part of layer

    zeps=SQRT(zod1*zod1-zod2*zod2)
! for a better vectorization
!   IF (zeps.LT.zargli) THEN
!     ze = EXP  (-zeps)
      ze = EXP  (-MIN(zeps,zargli))
!   ELSE
!     ze = ztsec
!   END IF
    zm = zod2/(zod1+zeps)
    pa4c(j1)=ze*(1.0_dp-(zm**2))*(1.0_dp/(1.0_dp-(zm**2)*(ze**2)))
    pa5c(j1)=zm*(1.0_dp-(ze**2))*(1.0_dp/(1.0_dp-(zm**2)*(ze**2)))

    ze1mwf = zeps / zod5
    zmu0if = ze1mwf + SIGN ( MAX(ABS(pqsmu0(j1)-ze1mwf),zepres) &
                            ,(pqsmu0(j1)-ze1mwf) )
    zod3 = zod3 * zmu0if
    zod4 = zod4 * zmu0if
    zod5 = zod5 * zmu0if
! for a better vectorization
!   IF (zod5.LT.ZARGLI) THEN
!      pa1c(j1) = EXP  (-zod5)
       pa1c(j1) = EXP  (-MIN(zod5,zargli))
!   ELSE
!      pa1c(j1) = ztsec
!   END IF
    zg1 = ( zod3*(zod5-zod1) -zod2*zod4) /(zod5*zod5 - zeps*zeps)
    zg2 =-( zod4*(zod5+zod1) +zod2*zod3) /(zod5*zod5 - zeps*zeps)
    pa2c(j1) = zg2*(pa1c(j1)-pa4c(j1)) -zg1*pa5c(j1)*pa1c(j1)
    pa3c(j1) = zg1*(1.0_dp-pa4c(j1)*pa1c(j1)) -zg2*pa5c(j1)
  ENDDO

  IF (ldebug .AND. jb==jindex) THEN
      print *,'**** cloudy     pa1c (j1b,jb)    =',pa1c (j1b)
      print *,'**** cloudy     pa2c (j1b,jb)    =',pa2c (j1b)
      print *,'**** cloudy     pa3c (j1b,jb)    =',pa3c (j1b)
      print *,'**** cloudy     pa4c (j1b,jb)    =',pa4c (j1b)
      print *,'**** cloudy     pa5c (j1b,jb)    =',pa5c (j1b)
  ENDIF

!-------------------------------------------------------------------------------
! End of the subroutine
!-------------------------------------------------------------------------------

END SUBROUTINE coe_so

!==============================================================================
!==============================================================================
!+ Module procedure in "radiation_rg"
!------------------------------------------------------------------------------

SUBROUTINE coe_th_gpu (                                 &
       pduh2oc,pduh2of,pduco2 ,pduo3  ,palogp ,palogt , &
       podsc  ,podsf  ,podac  ,podaf  ,pbsfc  ,pbsff  , &
       kspec  ,kh2o   ,kco2   ,ko3    ,                 &
       pa1c   ,pa1f   ,pa2c   ,pa2f   ,pa3c   ,pa3f)

!------------------------------------------------------------------------------
!
! Description:
!
!   The module procedure coe_th_gpu calculates the optical effects of atmospheric 
!   layers on thermal radiation based on basic optical properties of non-gaseous 
!   constituents and gaseous absorption coefficients selected through the 
!   corresponding control variables in the argument list.
!   This routine computes layer effects (transmissivity, reflectivity
!   and emmisivity) in the thermal part of the radiative spectrum
!   both for the cloud-free and the cloudy part of a model layer.
!   The calculation is based on the implicit delt-two-stream equations
!   (cf. Ritter and Geleyn, 1992) and uses basic optical properties
!   (i.e. absorption and scattering optical depth and backscattered
!   fraction for non-gaseous atmospheric constituents as well as 
!   gaseous absorption properties) as input. 
!
! Method:
!
! - addition of individual gaseous absorption effects to the optical
!   properties of the non-gaseous constituents
! - determination of layer effects (cf. Zdunkowski et al., 1982, 1986
!   and Ritter and Geleyn, 1992)
!     
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

! Input data
! ----------
  INTEGER, INTENT(IN) :: &
     kspec,       & ! spectral interval considered
     kh2o ,       & ! table index for h2o absorption properties
     kco2 ,       & ! table index for co2 absorption properties
     ko3            ! table index for o3  absorption properties
     
  REAL    (KIND=dp   ), INTENT (IN) ::  &

     ! opticall relevant gas quantities (Pa)
     pduh2oc, & ! h2o inside cloud
     pduh2of, & ! h2o out of cloud
     pduco2 , & ! co2 content 
     pduo3  , & ! o3  content 

     ! Logarithm of layer mean temperature and pressure
     palogt , & ! ln T
     palogp , & ! ln p

     ! Optical properties of non-gaseous constituents (..c=cloudy; ..f=free)  
     podsc  , & ! 
     podsf  , & ! 
     podac  , & ! 
     podaf  , & ! 
     pbsfc  , & ! 
     pbsff      ! 

! Output data
! -----------
  REAL    (KIND=dp   ), INTENT (OUT) ::  &
     pa1c  , & ! transmissivity in cloud   
     pa1f  , & ! transmissivity cloud-free  
     pa2c  , & ! reflectivity in cloud    
     pa2f  , & ! reflectivity cloud-free      
     pa3c  , & ! emissivity in cloud    
     pa3f      ! emissivity cloud-free       

! Local parameters: 
! ----------------
  REAL    (KIND=dp   ), PARAMETER ::  &
     zargli  = 80.0_dp     , &  ! argument limit for EXP 
     ztsec   = 1.0E-35_dp  , &  ! (=exp(-zargli) avoids ALOG(0.0)
     zodmax  = 1.0E+6_dp   , &  ! maximum allowed optical depth
     zudiff  = 2.0_dp      , &  ! Diffusivity factors for gases and other constituents
     zangfa  = 1.648721271_dp   ! exp(0.5)

  INTEGER, PARAMETER ::  &
     j1b    = 1,           & ! debug point index first dimension
     jb     = 1              ! debug point index second dimension

  REAL    (KIND=dp   ) ::  &
    zeps, ztau, zrho, zodgf, zodgc, zod1, zod2

! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine coe_th_gpu              
!------------------------------------------------------------------------------

  ! Optical depth of gases
      zodgf = 0.0_dp     ! Initialisation
 
      IF (kco2.ne.0) then     ! Include CO2 contribution
        zodgf = zodgf + pduco2 * (REAL(cobi(kco2,kspec,2),dp) &
                    * EXP ( REAL(coali(kco2,kspec,2),dp) * palogp    &
                           -REAL(cobti(kco2,kspec,2),dp) * palogt))
      ENDIF                  ! CO2
      !US IF (ldebug) print *,'**** zodgf(CO2)        =',zodgf
 
      IF (ko3 /= 0) THEN  ! Include O3 contribution
        zodgf = zodgf + pduo3  * (REAL(cobi(ko3 ,kspec,3),dp)* &
                    EXP ( REAL(coali(ko3 ,kspec,3),dp) * palogp      &
                         -REAL(cobti(ko3 ,kspec,3),dp) * palogt))
      ENDIF
      !US IF (ldebug) print *,'**** zodgf(CO2+O3)     =',zodgf
 
      ! Cloudy = cloud free for CO2 and O3 :
      zodgc = zodgf
 
      IF (kh2o /= 0) THEN  ! Include H2O contribution
        zodgf = zodgf + pduh2of* (REAL(cobi(kh2o,kspec,1),dp)* &
                    EXP ( REAL(coali(kh2o,kspec,1),dp) * palogp      &
                         -REAL(cobti(kh2o,kspec,1),dp) * palogt))
        zodgc = zodgc + pduh2oc* (REAL(cobi(kh2o,kspec,1),dp)* &
                    EXP ( REAL(coali(kh2o,kspec,1),dp) * palogp      &
                         -REAL(cobti(kh2o,kspec,1),dp) * palogt))
      ENDIF
!------------------------------------------------------------------------------
      !US IF (ldebug) print *,'**** zodgf(CO2+O3+H2O) =',zodgf
      !US IF (ldebug) print *,'**** zodgc(CO2+O3+H2O) =',zodgc
 
      zodgf = MIN (zodgf, zodmax)
      zodgc = MIN (zodgc, zodmax)

      !US IF (ldebug) print *,'**** nach securit auf optical depth '
      !US IF (ldebug) print *,'**** zodgf(CO2+O3+H2O) =',zodgf
      !US IF (ldebug) print *,'**** zodgc(CO2+O3+H2O) =',zodgc
 
      ! Pseudo-optical depth in cloud-free part of layer

      zod2 = zudiff * pbsff * podsf
      zod1 = zod2 + zudiff * podaf
      zod1 = zod1 + zangfa * zodgf

      !US IF (ldebug) THEN
      !US   print *,'**** cloud-free zod1 (j1b,jb)=',zod1
      !US   print *,'**** cloud-free zod2 (j1b,jb)=',zod2
      !US ENDIF 
 
      ! Layer coefficients in cloud-free part of layer
 
      zeps=SQRT(zod1*zod1-zod2*zod2)
      IF (zeps.LT.zargli) THEN
        ztau = EXP  (-zeps)
      ELSE
        ztau = ztsec
      END IF
      zrho = zod2/(zod1+zeps)
      pa1f=ztau*(1._dp-(zrho**2))*(1._dp/(1._dp-(zrho**2)*(ztau**2)))
      pa2f=zrho*(1._dp-(ztau**2))*(1._dp/(1._dp-(zrho**2)*(ztau**2)))
      pa3f=(1._dp-pa1f+pa2f)/(zod1+zod2)

      !US IF (ldebug) THEN
      !US   print *,'**** cloud-free pa1f (j1b,jb)=',pa1f (j1b,jb)
      !US   print *,'**** cloud-free pa2f (j1b,jb)=',pa2f (j1b,jb)
      !US   print *,'**** cloud-free pa3f (j1b,jb)=',pa3f (j1b,jb)
      !US ENDIF
 
      ! Pseudo-optical depth in cloudy part of layer
      zod2 = zudiff * pbsfc * podsc
      zod1 = zod2 + zudiff * podac
      zod1 = zod1 + zangfa * zodgc
 
      ! Layer coefficients in cloudy part of layer
 
      zeps=SQRT(zod1*zod1-zod2*zod2)
      IF (zeps.LT.zargli) THEN
        ztau = EXP  (-zeps)
      ELSE
        ztau = ztsec
      END IF
      zrho = zod2/(zod1+zeps)
      pa1c=ztau*(1._dp-(zrho**2))*(1._dp/(1._dp-(zrho**2)*(ztau**2)))
      pa2c=zrho*(1._dp-(ztau**2))*(1._dp/(1._dp-(zrho**2)*(ztau**2)))
      pa3c=(1._dp-pa1c+pa2c)/(zod1+zod2)

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE coe_th_gpu

!==============================================================================
!==============================================================================
!+ Module procedure in "radiation_rg"
!------------------------------------------------------------------------------

SUBROUTINE coe_so_gpu (                                                 &
       pduh2oc,pduh2of,pduco2 ,pduo3  ,palogp ,palogt ,                 &
       podsc  ,podsf  ,podac  ,podaf  ,pbsfc  ,pbsff  ,pusfc ,pusff ,   &
       psmu0  ,pqsmu0 ,                                                 &
       kspec  ,kh2o   ,kco2   ,ko3    ,                                 &
       pa1c   ,pa1f   ,pa2c   ,pa2f   ,pa3c   ,pa3f ,                   &
       pa4c   ,pa4f   ,pa5c   ,pa5f )

!------------------------------------------------------------------------------
!
! Description:
!
!   The module procedure coe_so calculates the optical effects of atmospheric
!   layers on solar radiation based on basic optical properties of non-gaseous
!   constituents and gaseous absorption coefficients selected through the
!   corresponding control variables.
!   This routine computes layer effects (transmissivity, reflectivity)
!   for diffuse and direct solar radiation both for the cloudy and the
!   cloud-free part of a model layer.
!   The calculation is based on the implicit delt-two-stream equations
!   (cf. Ritter and Geleyn, 1992) and uses basic optical properties
!   (i.e. absorption and scattering optical depth, backscattered and
!   upscattered fraction for non-gaseous atmospheric constituents and
!   gaseous absorption properties) as input.
!
! Method:
!
! - addition of individual gaseous absorption effects to the optical
!   properties of the non-gaseous constituents
!   (optical depth multiplied by alpha1 to alpha4)
!
! - determination of layer effects (cf. Zdunkowski et al., 1982, 1986
!   and Ritter and Geleyn, 1992)
! - the resonance case for those effects related to the direct solar
!   radiation is avoided by a small displacement of the local inverse
!   of the cosine of the zenith angle (if necessary)
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

! Input data
! ----------
  INTEGER, INTENT (IN) ::  &
     kspec,       & ! spectral interval considered
     kh2o ,       & ! table index for h2o absorption properties
     kco2 ,       & ! table index for co2 absorption properties
     ko3            ! table index for o3  absorption properties

  REAL    (KIND=dp   ), INTENT (IN) ::  &

     ! opticall relevant gas quantities (Pa)
     pduh2oc, & ! h2o inside cloud
     pduh2of, & ! h2o out of cloud
     pduco2 , & ! co2 content
     pduo3  , & ! o3  content

     ! Logarithm of layer mean temperature and pressure
     palogt , & ! ln T
     palogp , & ! ln p

     ! Optical properties of non-gaseous constituents (..c=cloudy; ..f=free)
     podsc  , & !
     podsf  , & !
     podac  , & !
     podaf  , & !
     pbsfc  , & !
     pbsff  , & !
     pusfc  , & !
     pusff  , & !

     psmu0              , & ! cosine of zenith angle
     pqsmu0                 ! inverse of cosine ...

! Output data
! -----------
  REAL    (KIND=dp   ), INTENT (OUT) ::  &
     pa1c  , & ! direct radiation transmis-
     pa1f  , & ! sivity cloudy/cloud-free
     pa2c  , & ! direct radition downward
     pa2f  , & ! scattering cloudy/cloud-free
     pa3c  , & ! direct radiation upward
     pa3f  , & ! scattering cloudy/cloud-free
     pa4c  , & ! diffuse flux transmissivity
     pa4f  , & ! cloudy/cloud-free
     pa5c  , & ! diffuse flux reflectivity
     pa5f      ! cloudy/cloud-free

! Local parameters:
! ----------------
  REAL    (KIND=dp   ), PARAMETER ::  &
     zargli  = 80.0_dp     , &  ! argument limit for EXP
     ztsec   = 1.0E-35_dp  , &  ! (=exp(-zargli) avoids ALOG(0.0)
     zodmax  = 1.0E+6_dp   , &  ! maximum allowed optical depth
     zepres  = 1.0E-7_dp   , &  ! for resonance case avoidance
                             ! 32bit-accuracy (1.E-14 for 64bit machine)
     zudiff  = 2.0_dp      , &  ! Diffusivity factors for gases and other constituents
     zangfa  = 1.648721271_dp   ! exp(0.5)

  INTEGER, PARAMETER ::  &
     j1b    = 1,           & ! debug point index first dimension
     jb     = 1              ! debug point index second dimension

! Local scalars:
! -------------
  INTEGER  ::  &
    j1,j2,j3                 ! loop indices over spatial dimensions

  REAL    (KIND=dp   ) ::  &
    zeps,                      & !
    ze,zm,zg1,zg2,ze1mwf,zmu0if  !

  REAL    (KIND=dp   ) ::  &
     zodgf, zodgc, zod1, zod2, zod3, zod4, zod5

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine coe_so_gpu
!------------------------------------------------------------------------------

  ! Optical depth of gases

      zodgf = 0.0_dp        ! Initialisation

      IF (kco2 /= 0) THEN     ! Include CO2 contribution
        zodgf = zodgf + pduco2* (REAL(cobi(kco2,kspec,2),dp)*       &
                    EXP ( REAL(coali(kco2,kspec,2),dp) * palogp     &
                         -REAL(cobti(kco2,kspec,2),dp) * palogt))
      ENDIF                  ! CO2
      !US IF (ldebug) print *,'**** zodgf(CO2)        =',zodgf

      IF (ko3 /= 0) THEN     ! Include O3 contribution
        zodgf = zodgf + pduo3 * (REAL(cobi(ko3 ,kspec,3),dp)*       &
                    EXP ( REAL(coali(ko3 ,kspec,3),dp) * palogp     &
                         -REAL(cobti(ko3 ,kspec,3),dp) * palogt))
      ENDIF
      !US IF (ldebug) print *,'**** zodgf(CO2+O3)     =',zodgf

      ! Cloudy = cloud free for CO2 and O3 :
      zodgc = zodgf

      IF (kh2o /= 0) THEN    ! Include H2O contribution
        zodgf = zodgf + pduh2of* (REAL(cobi(kh2o,kspec,1),dp)*       &
                    EXP ( REAL(coali(kh2o,kspec,1),dp) * palogp      &
                         -REAL(cobti(kh2o,kspec,1),dp) * palogt))
        zodgc = zodgc + pduh2oc* (REAL(cobi(kh2o,kspec,1),dp)*       &
                    EXP ( REAL(coali(kh2o,kspec,1),dp) * palogp      &
                         -REAL(cobti(kh2o,kspec,1),dp) * palogt))
      ENDIF
      !US IF (ldebug) print *,'**** zodgf(CO2+O3+H2O) =',zodgf
      !US IF (ldebug) print *,'**** zodgc(CO2+O3+H2O) =',zodgc

      zodgf = MIN (zodgf, zodmax)
      zodgc = MIN (zodgc, zodmax)
      !US IF (ldebug) print *,'**** nach securit auf optical depth '
      !US IF (ldebug) print *,'**** zodgf(CO2+O3+H2O) =',zodgf
      !US IF (ldebug) print *,'**** zodgc(CO2+O3+H2O) =',zodgc

      ! Pseudo-optical depth in cloud-free part of layer

      zod2 = zudiff * pbsff * podsf
      zod1 = zod2 + zudiff * podaf
      zod3 = pusff * podsf
      zod4 = podsf - zod3
      zod5 = podsf + podaf
      zod1 = zod1 + zangfa * zodgf
      zod5 = zod5 + zodgf
      !US IF (ldebug) THEN
      !US   print *,'**** cloud-free zod1 (j1b,jb)    =',zod1
      !US   print *,'**** cloud-free zod2 (j1b,jb)    =',zod2
      !US   print *,'**** cloud-free zod3 (j1b,jb)    =',zod3
      !US   print *,'**** cloud-free zod4 (j1b,jb)    =',zod4
      !US   print *,'**** cloud-free zod5 (j1b,jb)    =',zod5
      !US ENDIF

      ! Layer coefficients in cloud-free part of layer

      zeps=SQRT(zod1*zod1-zod2*zod2)
      IF (zeps.LT.zargli) THEN
        ze = EXP  (-zeps)
      ELSE
        ze = ztsec
      END IF
      zm = zod2/(zod1+zeps)
      pa4f=ze*(1._dp-(zm**2))*(1._dp/(1._dp-(zm**2)*(ze**2)))
      pa5f=zm*(1._dp-(ze**2))*(1._dp/(1._dp-(zm**2)*(ze**2)))

      ze1mwf = zeps / zod5
      zmu0if = ze1mwf + SIGN ( MAX(ABS(pqsmu0-ze1mwf),zepres) &
                              ,(pqsmu0-ze1mwf) )
      zod3 = zod3 * zmu0if
      zod4 = zod4 * zmu0if
      zod5 = zod5 * zmu0if
      IF (zod5.LT.zargli) THEN
        pa1f = EXP  (-zod5)
      ELSE
        pa1f = ztsec
      END IF
      zg1 = ( zod3*(zod5-zod1) -zod2*zod4) /(zod5*zod5 - zeps*zeps)
      zg2 =-( zod4*(zod5+zod1) +zod2*zod3) /(zod5*zod5 - zeps*zeps)
      pa2f = zg2*(pa1f-pa4f) -zg1*pa5f*pa1f
      pa3f = zg1*(1._dp-pa4f*pa1f) -zg2*pa5f

      !US IF (ldebug) THEN
      !US   print *,'**** cloud-free pa1f (j1b,jb)    =',pa1f (j1b,jb)
      !US   print *,'**** cloud-free pa2f (j1b,jb)    =',pa2f (j1b,jb)
      !US   print *,'**** cloud-free pa3f (j1b,jb)    =',pa3f (j1b,jb)
      !US   print *,'**** cloud-free pa4f (j1b,jb)    =',pa4f (j1b,jb)
      !US   print *,'**** cloud-free pa5f (j1b,jb)    =',pa5f (j1b,jb)
      !US ENDIF

      ! Pseudo-optical depth in cloudy part of layer

      zod2 = zudiff * pbsfc * podsc
      zod1 = zod2 + zudiff * podac
      zod3 = pusfc * podsc
      zod4 = podsc - zod3
      zod5 = podsc + podac
      zod1 = zod1 + zangfa * zodgc
      zod5 = zod5 + zodgc

      !US IF (ldebug) THEN
      !US   print *,'**** cloudy     zod1 (j1b,jb)    =',zod1
      !US   print *,'**** cloudy     zod2 (j1b,jb)    =',zod2
      !US   print *,'**** cloudy     zod3 (j1b,jb)    =',zod3
      !US   print *,'**** cloudy     zod4 (j1b,jb)    =',zod4
      !US   print *,'**** cloudy     zod5 (j1b,jb)    =',zod5
      !US ENDIF

      ! Layer coefficients in cloudy part of layer

      zeps=SQRT(zod1*zod1-zod2*zod2)
      IF (zeps.LT.zargli) THEN
        ze = EXP  (-zeps)
      ELSE
        ze = ztsec
      END IF
      zm = zod2/(zod1+zeps)
      pa4c=ze*(1._dp-(zm**2))*(1._dp/(1._dp-(zm**2)*(ze**2)))
      pa5c=zm*(1._dp-(ze**2))*(1._dp/(1._dp-(zm**2)*(ze**2)))

      ze1mwf = zeps / zod5
      zmu0if = ze1mwf + SIGN ( MAX(ABS(pqsmu0-ze1mwf),zepres) &
                              ,(pqsmu0-ze1mwf) )
      zod3 = zod3 * zmu0if
      zod4 = zod4 * zmu0if
      zod5 = zod5 * zmu0if
      IF (zod5.LT.ZARGLI) THEN
         pa1c = EXP  (-zod5)
      ELSE
         pa1c = ztsec
      END IF
      zg1 = ( zod3*(zod5-zod1) -zod2*zod4) /(zod5*zod5 - zeps*zeps)
      zg2 =-( zod4*(zod5+zod1) +zod2*zod3) /(zod5*zod5 - zeps*zeps)
      pa2c = zg2*(pa1c-pa4c) -zg1*pa5c*pa1c
      pa3c = zg1*(1._dp-pa4c*pa1c) -zg2*pa5c

!-------------------------------------------------------------------------------
! End of the subroutine
!-------------------------------------------------------------------------------

END SUBROUTINE coe_so_gpu

!==============================================================================
!==============================================================================
!+ Module procedure in "radiation_rg"
!------------------------------------------------------------------------------

SUBROUTINE inv_th (                                                   &
       pclc   ,pca1   ,pca2  ,pcb1  ,pcb2  ,pcc1  ,pcc2 ,pcd1 ,pcd2 , &
       pduh2oc,pduh2of,pduco2,pduo3 ,palogp,palogt,                   &
       podsc  ,podsf  ,podac ,podaf ,pbsfc ,pbsff ,                   &
       pbbr   ,palth,                                                 &
       kspec  ,kh2o   ,kco2  ,ko3   ,                                 &
       ki1sd  ,ki1ed  ,ki3sd ,ki3ed ,ki1sc ,ki1ec ,ki3sc ,ki3ec ,     &
       ldebug ,jindex ,                                               &
       pflcu  ,pflfu  ,pflcd ,pflfd)

!------------------------------------------------------------------------------
!
! Description:
!
!   The module procedure inv_th solves a linear equation system for thermal 
!   fluxes using a Gaussian elimination-backsubstitution algorithm dedicated 
!   to the specific structure of the system matrix.
!
! Method:
!
! - setting of the RHS of the system using the layer boundary black
!   body radiation and allowing for partial cloud cover in each layer
! - solution of the equation system including the lower boundary
!   condition
! - matrix coefficients are calculated in the course of the elimination
!   step for one layer at a time through a call to routine *coe_th*
! - the final result, i.e. the so-called black body flux differences
!   (cf.Ritter and Geleyn, 1992) are stored seperately for cloudy and
!   cloud-free part of each layer boundary
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

! Input data
! ----------
  INTEGER, INTENT (IN) ::  &
     ki1sd,       & ! start index for first  array dimension
     ki1ed,       & ! end   index for first  array dimension
     ki3sd,       & ! start index for third  array dimension
     ki3ed,       & ! end   index for third  array dimension

   ! and the same for the computations
     ki1sc,       & ! start index for first  array computation
     ki1ec,       & ! end   index for first  array computation
     ki3sc,       & ! start index for third  array computation
     ki3ec,       & ! end   index for third  array computation
     kspec,       & ! spectral interval considered
     kh2o ,       & ! table index for h2o absorption properties
     kco2 ,       & ! table index for co2 absorption properties
     ko3  ,       & ! table index for o3  absorption properties
     jindex         ! index for j-loop

  LOGICAL                 , INTENT (IN) ::  &
     ldebug         ! debug control switch       

  REAL    (KIND=dp)       , INTENT (IN) ::  &
     pclc  (ki1sd:ki1ed,            ki3sd:ki3ed),   & ! cloud cover
     pca1  (ki1sd:ki1ed,            ki3sd:ki3ed),   & ! cloud geometry factor  
     pca2  (ki1sd:ki1ed,            ki3sd:ki3ed),   & ! cloud geometry factor  
     pcb1  (ki1sd:ki1ed,            ki3sd:ki3ed),   & ! cloud geometry factor  
     pcb2  (ki1sd:ki1ed,            ki3sd:ki3ed),   & ! cloud geometry factor  
     pcc1  (ki1sd:ki1ed,            ki3sd:ki3ed),   & ! cloud geometry factor  
     pcc2  (ki1sd:ki1ed,            ki3sd:ki3ed),   & ! cloud geometry factor  
     pcd1  (ki1sd:ki1ed,            ki3sd:ki3ed),   & ! cloud geometry factor  
     pcd2  (ki1sd:ki1ed,            ki3sd:ki3ed),   & ! cloud geometry factor  
     pbbr  (ki1sd:ki1ed,            ki3sd:ki3ed+1), & ! black body radiation   
     palth (ki1sd:ki1ed)                              ! surface albedo

  ! Input data to be passed to *coe_th*
  REAL    (KIND=dp)       , INTENT (IN) ::  &

     ! layer gas contents (cloudy and cloud-free, if distinction necessary)
     pduh2oc(ki1sd:ki1ed,            ki3sd:ki3ed), & ! h2o-vapour cloudy      
     pduh2of(ki1sd:ki1ed,            ki3sd:ki3ed), & ! h2o-vapour cloud-free  
     pduco2 (ki1sd:ki1ed,            ki3sd:ki3ed), & ! co2
     pduo3  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! o3 
     ! optical properties of 'grey' constituents (cloudy and cloud-free)
     podsc  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! scattering optical depth
     podsf  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! scattering optical depth
     podac  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! absorption optical depth
     podaf  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! absorption optical depth
     pbsfc  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! backscatter fraction
     pbsff  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! backscatter fraction

     palogp (ki1sd:ki1ed,            ki3sd:ki3ed), & ! ln(p)
     palogt (ki1sd:ki1ed,            ki3sd:ki3ed)    ! ln(T)
 
! Output data
! -----------
  REAL    (KIND=dp)       , INTENT (OUT) ::  &
     pflcu (ki1sd:ki1ed,            ki3sd:ki3ed+1), & ! flux up   cloudy
     pflfu (ki1sd:ki1ed,            ki3sd:ki3ed+1), & ! flux up   cloud-free 
     pflcd (ki1sd:ki1ed,            ki3sd:ki3ed+1), & ! flux down cloudy     
     pflfd (ki1sd:ki1ed,            ki3sd:ki3ed+1)    ! flux down cloud-free 

! Local scalars:
! -------------
  INTEGER  ::  &
    j1,j3                    ! loop indices over spatial dimensions

  LOGICAL                  ::  &
    ldebug_coe_th            ! debug switch for *coe_th*

  REAL    (KIND=dp)        ::  &
    ztd1 ,ztd2 ,ztd3 ,ztd4 ,ztd5 ,ztd6 , ztd7,  & !
    ztds1,ztds2,ztds3,ztus1                       !
 
!- End of header
!==============================================================================
 
!------------------------------------------------------------------------------
! Start GPU data region
!------------------------------------------------------------------------------

  !$acc data                                                            &
  !---- Argument arrays - intent(in)
  !$acc present ( pclc,pca1,pca2,pcb1,pcb2,pcc1,pcc2,pcd1,pcd2,pbbr   ) &
  !$acc present ( palth,pduh2oc,pduh2of,pduco2,pduo3,podsc,podsf      ) &
  !$acc present ( podac,podaf,pbsfc,pbsff                             ) &
  !---- Argument arrays - intent(out)
  !$acc present ( pflcu,pflfu,pflcd,pflfd                             ) &
  !---- Local automatic arrays
  !$acc present ( pa1c,pa1f,pa2c,pa2f,pa3c,pa3f                       ) &
  !$acc present ( ztu1,ztu2,ztu3,ztu4,ztu5,ztu6,ztu7,ztu8,ztu9        ) &
  !---- Module arrays
  !$acc present ( cobi,coali,cobti                                    )

!------------------------------------------------------------------------------
! Begin Subroutine inv_th              
!------------------------------------------------------------------------------

  ldebug_coe_th = .FALSE. 

! Upper boundary condition
 
  !$acc parallel
  !$acc loop gang vector collapse(2)
  DO j3 = ki3sc, ki3ec+1
    DO j1 = ki1sc, ki1ec
      pflfd(j1,j3) = pbbr(j1,j3)
      pflcd(j1,j3) = 0.0_dp
    ENDDO
  ENDDO
  !$acc end parallel

#ifndef _OPENACC
  IF (ldebug .AND. jb==jindex) THEN
     print *,' *** INV_TH **************************'
     print *,' *** debug point : ',j1b,jb
     do j3 = ki3sc, ki3ec+1
       print *,'pflfd(j1b,',j3,') : ',pflfd(j1b,j3)
       print *,'pflcd(j1b,',j3,') : ',pflcd(j1b,j3)
     enddo
  ENDIF
#endif
 
! Determine effects of first layer in *coe_th*
#ifdef _OPENACC
  !$acc parallel
  !$acc loop gang vector
  DO j1 = ki1sc, ki1ec
    CALL coe_th_gpu(pduh2oc (j1,ki3sc), pduh2of (j1,ki3sc), &
                    pduco2  (j1,ki3sc), pduo3   (j1,ki3sc), &
                    palogp  (j1,ki3sc), palogt  (j1,ki3sc), &
                    podsc   (j1,ki3sc), podsf   (j1,ki3sc), &
                    podac   (j1,ki3sc), podaf   (j1,ki3sc), &
                    pbsfc   (j1,ki3sc), pbsff   (j1,ki3sc), &
                    kspec   , kh2o    , kco2    , ko3     , &
                    pa1c(j1), pa1f(j1), pa2c(j1),           &
                    pa2f(j1), pa3c(j1), pa3f(j1)            )
  ENDDO
  !$acc end parallel
#else
  CALL coe_th ( pduh2oc,pduh2of,pduco2 ,pduo3  ,palogp ,palogt , &
                podsc  ,podsf  ,podac  ,podaf  ,pbsfc  ,pbsff  , &
                ki3sc  ,kspec  ,kh2o   ,kco2   ,ko3    ,         &
                ki1sd  ,ki1ed  ,ki3sd  ,ki3ed  ,ki1sc  ,ki1ec  , &
                ldebug_coe_th  , jindex ,                        &
                pa1c   ,pa1f   ,pa2c   ,pa2f   ,pa3c   ,pa3f)
#endif
 
! Set RHS
  !$acc parallel
  !$acc loop gang vector
  DO j1 = ki1sc, ki1ec
    pflfu(j1,ki3sc)   = (1.0_dp-pclc(j1,ki3sc))*pa3f(j1)* &
                         (pbbr(j1,ki3sc)-pbbr(j1,ki3sc+1))
    pflcu(j1,ki3sc)   =   pclc(j1,ki3sc)*pa3c(j1)* &
                         (pbbr(j1,ki3sc)-pbbr(j1,ki3sc+1))
    pflfd(j1,ki3sc+1) = -pflfu(j1,ki3sc)
    pflcd(j1,ki3sc+1) = -pflcu(j1,ki3sc)
  ENDDO
  !$acc end parallel
 
! Elimination for first layer
  !$acc parallel
  !$acc loop gang vector
  DO j1 = ki1sc, ki1ec
    pflfu(j1,ki3sc)   = pflfu(j1,ki3sc  ) + pa2f(j1) *            &
                       (pca2 (j1,ki3sc)   * pflfd(j1,ki3sc))
    pflfd(j1,ki3sc+1) = pflfd(j1,ki3sc+1) + pa1f(j1) *            &
                       (pca2 (j1,ki3sc)   * pflfd(j1,ki3sc))
    pflcu(j1,ki3sc)   = pflcu(j1,ki3sc  ) + pa2c(j1) *            &
                       (pcb2 (j1,ki3sc)   * pflfd(j1,ki3sc))
    pflcd(j1,ki3sc+1) = pflcd(j1,ki3sc+1) + pa1c(j1) *            &
                       (pcb2 (j1,ki3sc)   * pflfd(j1,ki3sc))
  ENDDO
  !$acc end parallel

#ifndef _OPENACC
  IF (ldebug .AND. jb==jindex) THEN
     print *,' after elimination'
     print *,'pflfd(j1b,jb,ki3sc+1) : ',pflfd(j1b    ,ki3sc+1)
  ENDIF
#endif
 
! Store some utitlity variables for first layer

  !$acc parallel
  !$acc loop gang vector
  DO j1 = ki1sc, ki1ec
    ztu1(j1,ki3sc) = 0.0_dp
    ztu2(j1,ki3sc) = pca1(j1,ki3sc)*pa1f(j1)
    ztu3(j1,ki3sc) = pcc1(j1,ki3sc)*pa1f(j1)
    ztu4(j1,ki3sc) = pcb1(j1,ki3sc)*pa1c(j1)
    ztu5(j1,ki3sc) = pcd1(j1,ki3sc)*pa1c(j1)
    ztu6(j1,ki3sc) = pca1(j1,ki3sc)*pa2f(j1)
    ztu7(j1,ki3sc) = pcc1(j1,ki3sc)*pa2f(j1)
    ztu8(j1,ki3sc) = pcb1(j1,ki3sc)*pa2c(j1)
    ztu9(j1,ki3sc) = pcd1(j1,ki3sc)*pa2c(j1)
  ENDDO
  !$acc end parallel
 
! Vertical loop
 
!RUS: Turn loop structure with multiple ip loops inside a
!RUS: single k loop into perfectly nested k-ip loop on GPU.
#ifdef _OPENACC
  !$acc parallel
  !$acc loop gang vector
  DO j1 = ki1sc, ki1ec
    !$acc loop seq
    DO j3 = ki3sc+1, ki3ec

    ! Determine effect of the layer in *coe_th*
      CALL coe_th_gpu(pduh2oc    (j1,j3), pduh2of    (j1,j3), &
                      pduco2     (j1,j3), pduo3      (j1,j3), &
                      palogp     (j1,j3), palogt     (j1,j3), &
                      podsc      (j1,j3), podsf      (j1,j3), &
                      podac      (j1,j3), podaf      (j1,j3), &
                      pbsfc      (j1,j3), pbsff      (j1,j3), &
                      kspec   , kh2o    , kco2    , ko3     , &
                      pa1c(j1), pa1f(j1), pa2c(j1),           &
                      pa2f(j1), pa3c(j1), pa3f(j1)            )

      pflfu(j1,j3  ) = (1.0_dp - pclc(j1,j3)) * pa3f(j1)      &
                     * (pbbr(j1,j3) - pbbr(j1,j3+1))
      pflcu(j1,j3  ) =  pclc(j1,j3) * pa3c(j1)                &
                     * (pbbr(j1,j3) - pbbr(j1,j3+1))
      pflfd(j1,j3+1) = -pflfu(j1,j3)
      pflcd(j1,j3+1) = -pflcu(j1,j3)

  ! IF (ldebug .AND. jb==jindex) THEN
  !    print *,' in vertical loop j3=',j3
  !    print *,'pflfd(j1b,jb,j3+1) : ',pflfd(j1b    ,j3+1)
  ! ENDIF 

    ! Elimination and storage of utility variables
      ztd1 = 1.0_dp/(1.0_dp-pa2f(j1)*( pca2(j1,j3)*ztu6(j1,j3-1)    &
                                 +pcc2(j1,j3)*ztu8(j1,j3-1)))
      pflfu(j1,j3) = ztd1*( pflfu(j1,j3) +                          &
                     pa2f(j1)*( pca2(j1,j3)*pflfd(j1,j3)            &
                               +pcc2(j1,j3)*pflcd(j1,j3)))
      ztu1 (j1,j3) = ztd1*pa2f(j1)*( pca2(j1,j3)*ztu7(j1,j3-1)      &
                                    +pcc2(j1,j3)*ztu9(j1,j3-1))
      ztu2 (j1,j3) = ztd1*pa1f(j1)*pca1(j1,j3)
      ztu3 (j1,j3) = ztd1*pa1f(j1)*pcc1(j1,j3)
      ztd2 = pa2c(j1   )*(  pcb2(j1,j3)*ztu6(j1,j3-1)               &
                          + pcd2(j1,j3)*ztu8(j1,j3-1))
      ztd3 = 1.0_dp/(1.0_dp-pa2c(j1)*( pcb2(j1,j3)*ztu7(j1,j3-1)    &
                                      +pcd2(j1,j3)*ztu9(j1,j3-1))   &
                                      -ztd2*ztu1(j1,j3))
      pflcu(j1,j3) = ztd3*( pflcu(j1,j3) +                          &
                     pa2c(j1)*( pcb2(j1,j3)*pflfd(j1,j3)            &
                               +pcd2(j1,j3)*pflcd(j1,j3))           &
                       + ztd2*pflfu(j1,j3))
      ztu4 (j1,j3) = ztd3*( pa1c(j1)*pcb1(j1,j3)+ztd2*ztu2(j1,j3))
      ztu5 (j1,j3) = ztd3*( pa1c(j1)*pcd1(j1,j3)+ztd2*ztu3(j1,j3))
      ztd4 = pa1f(j1)*( pca2(j1,j3)*ztu6(j1,j3-1)                   &
                       +pcc2(j1,j3)*ztu8(j1,j3-1))
      ztd5 = pa1f(j1)*( pca2(j1,j3)*ztu7(j1,j3-1)                   &
                       +pcc2(j1,j3)*ztu9(j1,j3-1))
      pflfd(j1,j3+1) = pflfd(j1,j3+1)                               &
                       +pa1f(j1)*( pca2(j1,j3)*pflfd(j1,j3)         &
                                  +pcc2(j1,j3)*pflcd(j1,j3))        &
                       +ztd4*pflfu(j1,j3) + ztd5*pflcu(j1,j3)
      ztu6 (j1,j3) = pa2f(j1)*pca1(j1,j3)                           &
                    +ztd4*ztu2(j1,j3)+ztd5*ztu4(j1,j3)
      ztu7 (j1,j3) = pa2f(j1)*pcc1(j1,j3)                           &
                    +ztd4*ztu3(j1,j3)+ztd5*ztu5(j1,j3)
      ztd6 = pa1c(j1)*( pcb2(j1,j3)*ztu6(j1,j3-1)                   &
                       +pcd2(j1,j3)*ztu8(j1,j3-1))
      ztd7 = pa1c(j1)*( pcb2(j1,j3)*ztu7(j1,j3-1)                   &
                       +pcd2(j1,j3)*ztu9(j1,j3-1))
      pflcd(j1,j3+1) = pflcd(j1,j3+1)                               &
                      +pa1c(j1)*( pcb2(j1,j3)*pflfd(j1,j3)          &
                                 +pcd2(j1,j3)*pflcd(j1,j3))         &
                      +ztd6*pflfu(j1,j3) + ztd7*pflcu(j1,j3)
      ztu8(j1,j3) = pa2c(j1)*pcb1(j1,j3)                            &
                   +ztd6*ztu2(j1,j3)+ztd7*ztu4(j1,j3)
      ztu9(j1,j3) = pa2c(j1)*pcd1(j1,j3)                            &
                   +ztd6*ztu3(j1,j3)+ztd7*ztu5(j1,j3)
    ENDDO
  ! IF (ldebug .AND. jb==jindex) THEN
  !    print *,' after elimination in vertical loop j3=',j3
  !    print *,'pflfd(j1b,jb,j3+1) : ',pflfd(j1b    ,j3+1)
  ! ENDIF
  ENDDO     ! End of vertical loop over layers
  !$acc end parallel
#else
  DO j3 = ki3sc+1, ki3ec
 
    ! Determine effect of the layer in *coe_th*
    CALL coe_th ( pduh2oc,pduh2of,pduco2 ,pduo3  ,palogp ,palogt , &
                  podsc  ,podsf  ,podac  ,podaf  ,pbsfc  ,pbsff  , &
                  j3     ,kspec  ,kh2o   ,kco2   ,ko3    ,         &
                  ki1sd  ,ki1ed  ,ki3sd  ,ki3ed  ,ki1sc  ,ki1ec  , &
                  ldebug_coe_th  ,jindex ,                         &
                  pa1c   ,pa1f   ,pa2c   ,pa2f   ,pa3c   ,pa3f)
 
    ! Set RHS
    DO j1 = ki1sc, ki1ec
      pflfu(j1,j3  ) = (1.0_dp - pclc(j1,j3)) * pa3f(j1)      &
                     * (pbbr(j1,j3) - pbbr(j1,j3+1))
      pflcu(j1,j3  ) =  pclc(j1,j3) * pa3c(j1)                &
                     * (pbbr(j1,j3) - pbbr(j1,j3+1))
      pflfd(j1,j3+1) = -pflfu(j1,j3)
      pflcd(j1,j3+1) = -pflcu(j1,j3)
    ENDDO
 
    IF (ldebug .AND. jb==jindex) THEN
       print *,' in vertical loop j3=',j3
       print *,'pflfd(j1b,jb,j3+1) : ',pflfd(j1b    ,j3+1)
    ENDIF

    ! Elimination and storage of utility variables
 
    DO j1 = ki1sc, ki1ec
      ztd1 = 1.0_dp/(1.0_dp-pa2f(j1)*( pca2(j1,j3)*ztu6(j1,j3-1)    &
                                      +pcc2(j1,j3)*ztu8(j1,j3-1)))
      pflfu(j1,j3) = ztd1*( pflfu(j1,j3) +                          &
                     pa2f(j1)*( pca2(j1,j3)*pflfd(j1,j3)            &
                               +pcc2(j1,j3)*pflcd(j1,j3)))
      ztu1 (j1,j3) = ztd1*pa2f(j1)*( pca2(j1,j3)*ztu7(j1,j3-1)      &
                                    +pcc2(j1,j3)*ztu9(j1,j3-1))
      ztu2 (j1,j3) = ztd1*pa1f(j1)*pca1(j1,j3)
      ztu3 (j1,j3) = ztd1*pa1f(j1)*pcc1(j1,j3)
      ztd2 = pa2c(j1)*( pcb2(j1,j3)*ztu6(j1,j3-1)                   &
                      + pcd2(j1,j3)*ztu8(j1,j3-1))
      ztd3 = 1.0_dp/(1.0_dp-pa2c(j1)*( pcb2(j1,j3)*ztu7(j1,j3-1)    &
                                      +pcd2(j1,j3)*ztu9(j1,j3-1))   &
                                      -ztd2*ztu1(j1,j3))
      pflcu(j1,j3) = ztd3*( pflcu(j1,j3) +                          &
                    pa2c(j1)*( pcb2(j1,j3)*pflfd(j1,j3)             &
                              +pcd2(j1,j3)*pflcd(j1,j3))            &
                              + ztd2*pflfu(j1,j3))
      ztu4 (j1,j3) = ztd3*( pa1c(j1)*pcb1(j1,j3)+ztd2*ztu2(j1,j3))
      ztu5 (j1,j3) = ztd3*( pa1c(j1)*pcd1(j1,j3)+ztd2*ztu3(j1,j3))
      ztd4 = pa1f(j1)*( pca2(j1,j3)*ztu6(j1,j3-1)                   &
                       +pcc2(j1,j3)*ztu8(j1,j3-1))
      ztd5 = pa1f(j1)*( pca2(j1,j3)*ztu7(j1,j3-1)                   &
                       +pcc2(j1,j3)*ztu9(j1,j3-1))
      pflfd(j1,j3+1) = pflfd(j1,j3+1)                               &
                      +pa1f(j1)*( pca2(j1,j3)*pflfd(j1,j3)          &
                                 +pcc2(j1,j3)*pflcd(j1,j3))         &
                      + ztd4*pflfu(j1,j3) + ztd5*pflcu(j1,j3)
      ztu6 (j1,j3) = pa2f(j1)*pca1(j1,j3)                           &
                    +ztd4*ztu2(j1,j3)+ztd5*ztu4(j1,j3)
      ztu7 (j1,j3) = pa2f(j1)*pcc1(j1,j3)                           &
                    +ztd4*ztu3(j1,j3)+ztd5*ztu5(j1,j3)
      ztd6 = pa1c(j1)*( pcb2(j1,j3)*ztu6(j1,j3-1)                   &
                       +pcd2(j1,j3)*ztu8(j1,j3-1))
      ztd7 = pa1c(j1)*( pcb2(j1,j3)*ztu7(j1,j3-1)                   &
                       +pcd2(j1,j3)*ztu9(j1,j3-1))
      pflcd(j1,j3+1) = pflcd(j1,j3+1)                               &
                      +pa1c(j1)*( pcb2(j1,j3)*pflfd(j1,j3)          &
                                 +pcd2(j1,j3)*pflcd(j1,j3))         &
                      +ztd6*pflfu(j1,j3) + ztd7*pflcu(j1,j3)
      ztu8(j1,j3) = pa2c(j1)*pcb1(j1,j3)                            &
                   +ztd6*ztu2(j1,j3)+ztd7*ztu4(j1,j3)
      ztu9(j1,j3) = pa2c(j1)*pcd1(j1,j3)                            &
                   +ztd6*ztu3(j1,j3)+ztd7*ztu5(j1,j3)
    ENDDO

    IF (ldebug .AND. jb==jindex) THEN
       print *,' after elimination in vertical loop j3=',j3
       print *,'pflfd(j1b,jb,j3+1) : ',pflfd(j1b    ,j3+1)
    ENDIF
 
  ENDDO     ! End of vertical loop over layers
#endif
 
  ! Elimination and backsubstitution at the surface

  !$acc parallel
  !$acc loop gang vector
  DO j1 = ki1sc, ki1ec
    ztds1    =1.0_dp/(1.0_dp-palth(j1)*ztu6(j1,ki3ec))
    pflfu(j1,ki3ec+1)= ztds1 *palth(j1)*pflfd(j1,ki3ec+1)
    ztus1    =ztds1 *palth(j1)*ztu7(j1,ki3ec)
    ztds2    =palth(j1)*ztu8(j1,ki3ec)
    ztds3    =1.0_dp/(1.0_dp-palth(j1)*ztu9(j1,ki3ec)-ztds2*ztus1)
    pflcu(j1,ki3ec+1)=ztds3*( palth(j1)*pflcd(j1,ki3ec+1) &
                                 +ztds2*pflfu(j1,ki3ec+1))
    pflfu(j1,ki3ec+1)=pflfu(j1,ki3ec+1)+ztus1*pflcu(j1,ki3ec+1)
  ENDDO
  !$acc end parallel

!     Layer-by-layer backsubstitution
 
  !$acc parallel
  !$acc loop seq
  DO j3 =ki3ec,ki3sc,-1
    !$acc loop gang vector
    DO j1 = ki1sc, ki1ec
      pflcd(j1,j3+1) = pflcd(j1,j3+1) + ztu8 (j1,j3)               &
                 * pflfu(j1,j3+1) + ztu9 (j1,j3) * pflcu(j1,j3+1)
      pflfd(j1,j3+1) = pflfd(j1,j3+1) + ztu6 (j1,j3)               &
                 * pflfu(j1,j3+1) + ztu7 (j1,j3) * pflcu(j1,j3+1)
      pflcu(j1,j3  ) = pflcu(j1,j3  ) + ztu4 (j1,j3)               &
                 * pflfu(j1,j3+1) + ztu5 (j1,j3) * pflcu(j1,j3+1)
      pflfu(j1,j3  ) = pflfu(j1,j3  ) + ztu2 (j1,j3)               &
                 * pflfu(j1,j3+1) + ztu3 (j1,j3) * pflcu(j1,j3+1)  &
                                  + ztu1 (j1,j3) * pflcu(j1,j3)
    ENDDO

#ifndef _OPENACC
    IF (ldebug .AND. jb==jindex) THEN
      print *,' after backsubst.  in vertical loop j3=',j3
      print *,'pflfd(j1b,jb,j3+1) : ',pflfd(j1b    ,j3+1)
    ENDIF
#endif
  ENDDO
  !$acc end parallel
 
!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

  !$acc end data

END SUBROUTINE inv_th

!==============================================================================
!==============================================================================
!+ Module procedure in "Radiation" 
!------------------------------------------------------------------------------

SUBROUTINE inv_so (                                                    &
       pclc   ,pca1   ,pca2  ,pcb1  ,pcb2  ,pcc1  ,pcc2 ,pcd1 ,pcd2 ,  &
       pflpt  ,psmu0  ,pqsmu0,palp  ,palso ,                           &
       pduh2oc,pduh2of,pduco2,pduo3 ,palogp,palogt,                    &
       podsc  ,podsf  ,podac ,podaf ,pbsfc ,pbsff ,pusfc,pusff,        &
       kspec  ,kh2o   ,kco2  ,ko3   ,                                  &
       ki1sd  ,ki1ed  ,ki3sd ,ki3ed ,ki1sc ,ki1ec ,ki3sc ,ki3ec ,      &
       ldebug ,jindex ,                                                &
       pflcu  ,pflfu  ,pflcd ,pflfd ,pflcp ,pflfp)

!------------------------------------------------------------------------------
!
! Description:
!
!   The module procedure inv_so solves the linear system of equations for 
!   solar fluxes.
!   The routine solves a linear equation system for solar fluxes using
!   a Gaussian elimination-backsubstitution algorithm dedicated to the
!   specific structure of the system matrix.
!
! Method:
!
! - setting of the RHS of the system using the parallel solar radiation
!   at the top of the atmosphere and allowing for partial cloud cover
! - solution of the equation system including the lower boundary
!   condition
! - matrix coefficients are calculated in the course of the elimination
!   step for one layer at a time through a call to routine *coe_so*
! - the final result, i.e. upward and downward diffuse and parallel 
!   solar fluxes are stored seperately for cloudy and cloud-free parts
!   of each layer boundary
!     
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

! Input data
! ----------
  INTEGER, INTENT (IN) ::  &
     ki1sd,       & ! start index for first  array dimension
     ki1ed,       & ! end   index for first  array dimension
     ki3sd,       & ! start index for third  array dimension
     ki3ed,       & ! end   index for third  array dimension

   ! and the same for the computations
     ki1sc,       & ! start index for first  array computation
     ki1ec,       & ! end   index for first  array computation
     ki3sc,       & ! start index for third  array computation
     ki3ec,       & ! end   index for third  array computation
     kspec,       & ! spectral interval considered
     kh2o ,       & ! table index for h2o absorption properties
     kco2 ,       & ! table index for co2 absorption properties
     ko3  ,       & ! table index for o3  absorption properties
     jindex         ! index for j-loop

  LOGICAL                 , INTENT (IN) ::  &
     ldebug         ! debug control switch       

  REAL    (KIND=dp)       , INTENT (IN) ::  &
     pclc  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! cloud cover
     pca1  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! cloud geometry factor  
     pca2  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! cloud geometry factor  
     pcb1  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! cloud geometry factor  
     pcb2  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! cloud geometry factor  
     pcc1  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! cloud geometry factor  
     pcc2  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! cloud geometry factor  
     pcd1  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! cloud geometry factor  
     pcd2  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! cloud geometry factor  

     pflpt (ki1sd:ki1ed)            , &  ! parallel solar flux at TOA
     palp  (ki1sd:ki1ed)            , &  ! surface albedo for parallel
     palso (ki1sd:ki1ed)                 ! and for diffuse radiation  

  ! Input data to be passed to *coe_so*
  REAL    (KIND=dp)       , INTENT (IN) ::  &

     ! layer gas contents (cloudy and cloud-free, if distinction necessary)
     pduh2oc(ki1sd:ki1ed,            ki3sd:ki3ed), & ! h2o-vapour cloudy      
     pduh2of(ki1sd:ki1ed,            ki3sd:ki3ed), & ! h2o-vapour cloud-free  
     pduco2 (ki1sd:ki1ed,            ki3sd:ki3ed), & ! co2
     pduo3  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! o3 
     ! optical properties of 'grey' constituents (cloudy and cloud-free)
     podsc  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! scattering optical depth
     podsf  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! scattering optical depth
     podac  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! absorption optical depth
     podaf  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! absorption optical depth
     pbsfc  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! backscatter fraction
     pbsff  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! backscatter fraction
     pusfc  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! upscatter   fraction
     pusff  (ki1sd:ki1ed,            ki3sd:ki3ed), & ! upscatter   fraction

     palogp (ki1sd:ki1ed,            ki3sd:ki3ed), & ! ln(p)
     palogt (ki1sd:ki1ed,            ki3sd:ki3ed), & ! ln(T)

     psmu0 (ki1sd:ki1ed)            , & ! cosine of zenith angle
     pqsmu0(ki1sd:ki1ed)                ! 1./cosine of zenith angle

! Output data
! -----------
  REAL    (KIND=dp)       , INTENT (OUT) ::  &
     pflcu (ki1sd:ki1ed,            ki3sd:ki3ed+1), & ! flux up   cloudy
     pflfu (ki1sd:ki1ed,            ki3sd:ki3ed+1), & ! flux up   cloud-free 
     pflcd (ki1sd:ki1ed,            ki3sd:ki3ed+1), & ! flux down cloudy     
     pflfd (ki1sd:ki1ed,            ki3sd:ki3ed+1), & ! flux down cloud-free 
     pflcp (ki1sd:ki1ed,            ki3sd:ki3ed+1), & ! flux par. cloudy     
     pflfp (ki1sd:ki1ed,            ki3sd:ki3ed+1)    ! flux par. cloud-free 

! Local scalars:
! -------------
  INTEGER  ::  &
    j1,j3                    ! loop indices over spatial dimensions

  LOGICAL                  ::  &
    ldebug_coe_so            ! debug switch for *coe_so*

  REAL    (KIND=dp)        ::  &
    ztd1 ,ztd2 ,ztd3 ,ztd4 ,ztd5 ,ztd6 , ztd7,  & !
    ztds1,ztds2,ztds3,ztus1                       !
 
!- End of header
!==============================================================================
 
!------------------------------------------------------------------------------
! Start GPU data region
!------------------------------------------------------------------------------

  !$acc data                                                            &
  !---- Argument arrays - intent(in)
  !$acc present ( pclc,pca1,pca2,pcb1,pcb2,pcc1,pcc2,pcd1,pcd2,pflpt  ) &
  !$acc present ( palp,palso,pduh2oc,pduh2of,pduco2,pduo3,podsc,podsf ) &
  !$acc present ( podac,podaf,pbsfc,pbsff,pusfc,pusff,palogp,palogt   ) &
  !$acc present ( psmu0,pqsmu0                                        ) &
  !---- Argument arrays - intent(out)
  !$acc present ( pflcu,pflfu,pflcd,pflfd,pflcp,pflfp                 ) &
  !---- Local automatic arrays
  !$acc present ( pa1c,pa1f,pa2c,pa2f,pa3c,pa4f,pa4c,pa5f,pa5c,pa3f   ) &
  !$acc present ( ztu1,ztu2,ztu3,ztu4,ztu5,ztu6,ztu7,ztu8,ztu9        ) &
  !---- Module arrays
  !$acc present ( cobi,coali,cobti                                    )

!------------------------------------------------------------------------------
! Begin Subroutine inv_so              
!------------------------------------------------------------------------------

  ldebug_coe_so = .FALSE. 

  ! Upper boundary condition
 
  !$acc parallel
  !$acc loop gang vector
  DO j1 = ki1sc, ki1ec
    pflfp(j1,ki3sc) = pflpt(j1)
    pflcp(j1,ki3sc) = 0.0_dp
    pflfd(j1,ki3sc) = 0.0_dp
    pflcd(j1,ki3sc) = 0.0_dp
  ENDDO
  !$acc end parallel

#ifndef _OPENACC
  IF (ldebug .AND. jb==jindex) THEN
     print *,' *** INV_SO **************************'
     print *,' *** Debug point: ',j1b,jb
     print *,'pflfp(j1b,jb,ki3sc) : ',pflfp(j1b    ,ki3sc)
     print *,'pflcp(j1b,jb,ki3sc) : ',pflcp(j1b    ,ki3sc)
     print *,'pflfp(j1b,jb,ki3sc) : ',pflfp(j1b    ,ki3sc)
     print *,'pflcd(j1b,jb,ki3sc) : ',pflcd(j1b    ,ki3sc)
  ENDIF
#endif

  ! Determine effects of first layer in *coe_so*
#ifdef _OPENACC
  !$acc parallel
  !$acc loop gang vector
  DO j1 = ki1sc, ki1ec
    CALL coe_so_gpu(pduh2oc (j1,ki3sc), pduh2of (j1,ki3sc), &
                    pduco2  (j1,ki3sc), pduo3   (j1,ki3sc), &
                    palogp  (j1,ki3sc), palogt  (j1,ki3sc), &
                    podsc   (j1,ki3sc), podsf   (j1,ki3sc), &
                    podac   (j1,ki3sc), podaf   (j1,ki3sc), &
                    pbsfc   (j1,ki3sc), pbsff   (j1,ki3sc), &
                    pusfc   (j1,ki3sc), pusff   (j1,ki3sc), &
                    psmu0   (j1)      , pqsmu0  (j1)      , &
                    kspec   , kh2o    , kco2    , ko3     , &
                    pa1c(j1), pa1f(j1), pa2c(j1),           &
                    pa2f(j1), pa3c(j1), pa3f(j1),           &
                    pa4c(j1), pa4f(j1), pa5c(j1), pa5f(j1)  )
  ENDDO
  !$acc end parallel
#else
  CALL  coe_so (                                                     &
      pduh2oc,pduh2of,pduco2 ,pduo3  ,palogp ,palogt ,               &
      podsc  ,podsf  ,podac  ,podaf  ,pbsfc  ,pbsff  ,pusfc ,pusff , &
      psmu0  ,pqsmu0 ,                                               &
      ki3sc  ,kspec  ,kh2o   ,kco2   ,ko3    ,                       &
      ki1sd  ,ki1ed  ,ki3sd  ,ki3ed  ,ki1sc  ,ki1ec  ,               &
      ldebug_coe_so  ,jindex ,                                       &
      pa1c   ,pa1f   ,pa2c   ,pa2f   ,pa3c   ,pa3f ,                 &
      pa4c   ,pa4f   ,pa5c   ,pa5f )
#endif
 
  ! Top layer elimination
 
  !$acc parallel
  !$acc loop gang vector
  DO j1 = ki1sc, ki1ec
    pflfu(j1,ki3sc  ) = pa3f(j1) * pca2(j1,ki3sc) *pflfp(j1,ki3sc)
    pflfp(j1,ki3sc+1) = pa1f(j1) * pca2(j1,ki3sc) *pflfp(j1,ki3sc)
    pflfd(j1,ki3sc+1) = pa2f(j1) * pca2(j1,ki3sc) *pflfp(j1,ki3sc)
    pflcu(j1,ki3sc  ) = pa3c(j1) * pcb2(j1,ki3sc) *pflfp(j1,ki3sc)
    pflcp(j1,ki3sc+1) = pa1c(j1) * pcb2(j1,ki3sc) *pflfp(j1,ki3sc)
    pflcd(j1,ki3sc+1) = pa2c(j1) * pcb2(j1,ki3sc) *pflfp(j1,ki3sc)
  ENDDO
  !$acc end parallel

#ifndef _OPENACC
  IF (ldebug .AND. jb==jindex) THEN
     print *,' *** INV_SO **************************'
     print *,'pflfu(j1b,jb,ki3sc)  : ',pflfu(j1b    ,ki3sc)
     print *,'pflcu(j1b,jb,ki3sc)  : ',pflcu(j1b    ,ki3sc)
     print *,'pflfd(j1b,jb,ki3sc+1): ',pflfd(j1b    ,ki3sc+1)
     print *,'pflcd(j1b,jb,ki3sc+1): ',pflcd(j1b    ,ki3sc+1)
     print *,'pa1f (j1b,jb)        : ',pa1f (j1b    )        
     print *,'pa1c (j1b,jb)        : ',pa1c (j1b    )        
     print *,'pa2f (j1b,jb)        : ',pa2f (j1b    )        
     print *,'pa2c (j1b,jb)        : ',pa2c (j1b    )        
     print *,'pa3f (j1b,jb)        : ',pa3f (j1b    )        
     print *,'pa3c (j1b,jb)        : ',pa3c (j1b    )        
  ENDIF 
#endif
 
  ! Storage of utility arrays for the top layer

  !$acc parallel
  !$acc loop gang vector
  DO j1 = ki1sc, ki1ec
    ztu1(j1,1) = 0.0_dp
    ztu2(j1,1) = pca1(j1,1) * pa4f(j1)
    ztu3(j1,1) = pcc1(j1,1) * pa4f(j1)
    ztu4(j1,1) = pcb1(j1,1) * pa4c(j1)
    ztu5(j1,1) = pcd1(j1,1) * pa4c(j1)
    ztu6(j1,1) = pca1(j1,1) * pa5f(j1)
    ztu7(j1,1) = pcc1(j1,1) * pa5f(j1)
    ztu8(j1,1) = pcb1(j1,1) * pa5c(j1)
    ztu9(j1,1) = pcd1(j1,1) * pa5c(j1)
  ENDDO
  !$acc end parallel
 
  ! Successive layer-by-layer elimination
 
!RUS: Turn loop structure with multiple ip loops inside a
!RUS: single k loop into perfectly nested k-ip loop on GPU.
#ifdef _OPENACC
  !$acc parallel
  !$acc loop gang vector
  DO j1 = ki1sc, ki1ec
    !$acc loop seq
    DO j3 = ki3sc+1, ki3ec         ! Loop over vertical

    ! Determine effects of layer in *coe_so*
      CALL coe_so_gpu(pduh2oc (j1,j3)   , pduh2of (j1,j3)   , &
                      pduco2  (j1,j3)   , pduo3   (j1,j3)   , &
                      palogp  (j1,j3)   , palogt  (j1,j3)   , &
                      podsc   (j1,j3)   , podsf   (j1,j3)   , &
                      podac   (j1,j3)   , podaf   (j1,j3)   , &
                      pbsfc   (j1,j3)   , pbsff   (j1,j3)   , &
                      pusfc   (j1,j3)   , pusff   (j1,j3)   , &
                      psmu0   (j1)      , pqsmu0  (j1)      , &
                      kspec   , kh2o    , kco2    , ko3     , &
                      pa1c(j1), pa1f(j1), pa2c(j1),           &
                      pa2f(j1), pa3c(j1), pa3f(j1),           &
                      pa4c(j1), pa4f(j1), pa5c(j1), pa5f(j1)  )

    ! Elimination

      ztd1 = 1.0_dp/(1.0_dp-pa5f(j1)*(pca2(j1,j3)*ztu6(j1,j3-1)          &
                                     +pcc2(j1,j3)*ztu8(j1,j3-1)))
      pflfu(j1,j3) = ztd1*( pa3f(j1)*( pca2(j1,j3)*pflfp(j1,j3)          &
                                      +pcc2(j1,j3)*pflcp(j1,j3) )        &
                           +pa5f(j1)*( pca2(j1,j3)*pflfd(j1,j3)          &
                                      +pcc2(j1,j3)*pflcd(j1,j3) )   )
      ztu1 (j1,j3) = ztd1*pa5f(j1)* ( pca2(j1,j3)*ztu7(j1,j3-1)          &
                                     +pcc2(j1,j3)*ztu9(j1,j3-1))
      ztu2(j1,j3)  =  ztd1*pa4f(j1)*pca1(j1,j3)
      ztu3(j1,j3)  =  ztd1*pa4f(j1)*pcc1(j1,j3)
      ztd2         = pa5c(j1)*( pcb2(j1,j3)*ztu6(j1,j3-1)                &
                               +pcd2(j1,j3)*ztu8(j1,j3-1) )
      ztd3         = 1.0_dp/( 1.0_dp-pa5c(j1)*(pcb2(j1,j3)*ztu7(j1,j3-1) &
                                              +pcd2(j1,j3)*ztu9(j1,j3-1))&
                             -ztd2*ztu1(j1,j3) )
      pflcu(j1,j3) = ztd3 *(                                             &
                              pa3c(j1)*( pcb2(j1,j3)*pflfp(j1,j3)        &
                                        +pcd2(j1,j3)*pflcp(j1,j3) )      &
                             +pa5c(j1)*( pcb2(j1,j3)*pflfd(j1,j3)        &
                                        +pcd2(j1,j3)*pflcd(j1,j3) )      &
                             +ztd2*pflfu(j1,j3)   )
      ztu4(j1,j3)    = ztd3 *( pa4c(j1)*pcb1(j1,j3)+ztd2*ztu2(j1,j3) )
      ztu5(j1,j3)    = ztd3 *( pa4c(j1)*pcd1(j1,j3)+ztd2*ztu3(j1,j3) )
      pflfp(j1,j3+1) = pa1f(j1)*(pca2(j1,j3)*pflfp(j1,j3)                &
                                +pcc2(j1,j3)*pflcp(j1,j3))
      pflcp(j1,j3+1) = pa1c(j1)*(pcb2(j1,j3)*pflfp(j1,j3)                &
                                +pcd2(j1,j3)*pflcp(j1,j3))
      ztd4 = pa4f(j1)*( pca2(j1,j3)*ztu6(j1,j3-1)                        &
                       +pcc2(j1,j3)*ztu8(j1,j3-1) )
      ztd5 = pa4f(j1)*( pca2(j1,j3)*ztu7(j1,j3-1)                        &
                       +pcc2(j1,j3)*ztu9(j1,j3-1) )
      pflfd(j1,j3+1) = pa2f(j1)*(pca2(j1,j3)*pflfp(j1,j3)                &
                                +pcc2(j1,j3)*pflcp(j1,j3))               &
                      +pa4f(j1)*(pca2(j1,j3)*pflfd(j1,j3)                &
                                +pcc2(j1,j3)*pflcd(j1,j3))               &
                      +ztd4*pflfu(j1,j3) + ztd5*pflcu(j1,j3)
      ztu6(j1,j3) = pa5f(j1)*pca1(j1,j3)                                 &
                      +ztd4*ztu2(j1,j3) + ztd5*ztu4(j1,j3)
      ztu7(j1,j3) = pa5f(j1)*pcc1(j1,j3)                                 &
                      +ztd4*ztu3(j1,j3) + ztd5*ztu5(j1,j3)
      ztd6 = pa4c(j1)*( pcb2(j1,j3)*ztu6(j1,j3-1)                        &
                       +pcd2(j1,j3)*ztu8(j1,j3-1) )
      ztd7 = pa4c(j1)*( pcb2(j1,j3)*ztu7(j1,j3-1)                        &
                        +pcd2(j1,j3)*ztu9(j1,j3-1) )
      pflcd(j1,j3+1) =  pa2c(j1)*(pcb2(j1,j3)*pflfp(j1,j3)               &
                                 +pcd2(j1,j3)*pflcp(j1,j3))              &
                         + pa4c(j1)*(pcb2(j1,j3)*pflfd(j1,j3)            &
                                 +pcd2(j1,j3)*pflcd(j1,j3))              &
                         +ztd6*pflfu(j1,j3) + ztd7*pflcu(j1,j3)
      ztu8(j1,j3) = pa5c(j1)*pcb1(j1,j3)                                 &
                   +ztd6*ztu2(j1,j3) + ztd7*ztu4(j1,j3)
      ztu9(j1,j3) = pa5c(j1)*pcd1(j1,j3)                                 &
                   +ztd6*ztu3(j1,j3) + ztd7*ztu5(j1,j3)
    ENDDO
  END DO       ! Vertical loop
  !$acc end parallel
#else
  DO j3 = ki3sc+1, ki3ec         ! Loop over vertical

    ! Determine effects of layer in *coe_so*
    CALL  coe_so (                                                       &
          pduh2oc,pduh2of,pduco2 ,pduo3  ,palogp ,palogt ,               &
          podsc  ,podsf  ,podac  ,podaf  ,pbsfc  ,pbsff  ,pusfc ,pusff , &
          psmu0  ,pqsmu0 ,                                               &
          j3     ,kspec  ,kh2o   ,kco2   ,ko3    ,                       &
          ki1sd  ,ki1ed  ,ki3sd  ,ki3ed  ,ki1sc  ,ki1ec  ,               &
          ldebug_coe_so  ,jindex ,                                       &
          pa1c   ,pa1f   ,pa2c   ,pa2f   ,pa3c   ,pa3f ,                 &
          pa4c   ,pa4f   ,pa5c   ,pa5f )
 
    ! Elimination
 
    DO j1 = ki1sc, ki1ec
      ztd1 = 1.0_dp/(1.0_dp-pa5f(j1)*(pca2(j1,j3)*ztu6(j1,j3-1)          &
                                     +pcc2(j1,j3)*ztu8(j1,j3-1)))
      pflfu(j1,j3) = ztd1*( pa3f(j1)*( pca2(j1,j3)*pflfp(j1,j3)          &
                                      +pcc2(j1,j3)*pflcp(j1,j3) )        &
                           +pa5f(j1)*( pca2(j1,j3)*pflfd(j1,j3)          &
                                      +pcc2(j1,j3)*pflcd(j1,j3) )   )
      ztu1 (j1,j3) = ztd1*pa5f(j1)* ( pca2(j1,j3)*ztu7(j1,j3-1)          &
                                     +pcc2(j1,j3)*ztu9(j1,j3-1))
      ztu2(j1,j3)  =  ztd1*pa4f(j1)*pca1(j1,j3)
      ztu3(j1,j3)  =  ztd1*pa4f(j1)*pcc1(j1,j3)
      ztd2         = pa5c(j1)*( pcb2(j1,j3)*ztu6(j1,j3-1)                &
                               +pcd2(j1,j3)*ztu8(j1,j3-1) )
      ztd3         = 1.0_dp/( 1.0_dp-pa5c(j1)*(pcb2(j1,j3)*ztu7(j1,j3-1) &
                                              +pcd2(j1,j3)*ztu9(j1,j3-1))&
                             - ztd2*ztu1(j1,j3) )
      pflcu(j1,j3) = ztd3 *(   pa3c(j1)*( pcb2(j1,j3)*pflfp(j1,j3)       &
                                         +pcd2(j1,j3)*pflcp(j1,j3) )     &
                              +pa5c(j1)*( pcb2(j1,j3)*pflfd(j1,j3)       &
                                         +pcd2(j1,j3)*pflcd(j1,j3) )     &
                              +ztd2*pflfu(j1,j3)   )
      ztu4(j1,j3)  = ztd3 *( pa4c(j1)*pcb1(j1,j3)+ztd2*ztu2(j1,j3) )
      ztu5(j1,j3)  = ztd3 *( pa4c(j1)*pcd1(j1,j3)+ztd2*ztu3(j1,j3) )
      pflfp(j1,j3+1) = pa1f(j1)*(pca2(j1,j3)*pflfp(j1,j3)                &
                                +pcc2(j1,j3)*pflcp(j1,j3))
      pflcp(j1,j3+1) = pa1c(j1)*(pcb2(j1,j3)*pflfp(j1,j3)                &
                                +pcd2(j1,j3)*pflcp(j1,j3))
      ztd4 = pa4f(j1)*( pca2(j1,j3)*ztu6(j1,j3-1)                        &
                       +pcc2(j1,j3)*ztu8(j1,j3-1) )
      ztd5 = pa4f(j1)*( pca2(j1,j3)*ztu7(j1,j3-1)                        &
                       +pcc2(j1,j3)*ztu9(j1,j3-1) )
      pflfd(j1,j3+1) = pa2f(j1)*(pca2(j1,j3)*pflfp(j1,j3)                &
                                +pcc2(j1,j3)*pflcp(j1,j3))               &
                      +pa4f(j1)*(pca2(j1,j3)*pflfd(j1,j3)                &
                                +pcc2(j1,j3)*pflcd(j1,j3))               &
                         +ztd4*pflfu(j1,j3) + ztd5*pflcu(j1,j3)
      ztu6(j1,j3) = pa5f(j1)*pca1(j1,j3)                                 &
                      + ztd4*ztu2(j1,j3) + ztd5*ztu4(j1,j3)
      ztu7(j1,j3) = pa5f(j1)*pcc1(j1,j3)                                 &
                      + ztd4*ztu3(j1,j3) + ztd5*ztu5(j1,j3)
      ztd6 = pa4c(j1)*( pcb2(j1,j3)*ztu6(j1,j3-1)                        &
                       +pcd2(j1,j3)*ztu8(j1,j3-1) )
      ztd7 = pa4c(j1)*( pcb2(j1,j3)*ztu7(j1,j3-1)                        &
                       +pcd2(j1,j3)*ztu9(j1,j3-1) )
      pflcd(j1,j3+1) =  pa2c(j1)*(pcb2(j1,j3)*pflfp(j1,j3)               &
                                 +pcd2(j1,j3)*pflcp(j1,j3))              &
                      + pa4c(j1)*(pcb2(j1,j3)*pflfd(j1,j3)               &
                                 +pcd2(j1,j3)*pflcd(j1,j3))              &
                         +ztd6*pflfu(j1,j3) + ztd7*pflcu(j1,j3)
      ztu8(j1,j3) = pa5c(j1)*pcb1(j1,j3)                                 &
                      + ztd6*ztu2(j1,j3) + ztd7*ztu4(j1,j3)
      ztu9(j1,j3) = pa5c(j1)*pcd1(j1,j3)                                 &
                      + ztd6*ztu3(j1,j3) + ztd7*ztu5(j1,j3)
    ENDDO

    IF (ldebug .AND. jb==jindex) THEN
      print *,' inv_so j3=',j3, '; jb=',jb
      print *,'pflfu(j1b,jb,j3)  : ',pflfu(j1b,    j3)
      print *,'pflcu(j1b,jb,j3)  : ',pflcu(j1b,    j3)
      print *,'pflfd(j1b,jb,j3+1): ',pflfd(j1b,    j3+1)
      print *,'pflcd(j1b,jb,j3+1): ',pflcd(j1b,    j3+1)
      print *,'pflfp(j1b,jb,j3+1): ',pflfp(j1b,    j3+1)
      print *,'pflcp(j1b,jb,j3+1): ',pflcp(j1b,    j3+1)
      print *,'ztu1 (j1b,jb,j3)  : ',ztu1 (j1b,    j3)
      print *,'ztu2 (j1b,jb,j3)  : ',ztu2 (j1b,    j3)
      print *,'ztu3 (j1b,jb,j3)  : ',ztu3 (j1b,    j3)
      print *,'ztu4 (j1b,jb,j3)  : ',ztu4 (j1b,    j3)
      print *,'ztu5 (j1b,jb,j3)  : ',ztu5 (j1b,    j3)
      print *,'ztu6 (j1b,jb,j3)  : ',ztu6 (j1b,    j3)
      print *,'ztu7 (j1b,jb,j3)  : ',ztu7 (j1b,    j3)
      print *,'ztu8 (j1b,jb,j3)  : ',ztu8 (j1b,    j3)
      print *,'ztu9 (j1b,jb,j3)  : ',ztu9 (j1b,    j3)
      print *,' .....'
    ENDIF

  END DO       ! Vertical loop
#endif
 
  ! Elimination and back-substitution at surface
 
  !$acc parallel
  !$acc loop gang vector
  DO j1 = ki1sc, ki1ec
    ztds1  = 1.0_dp/(1.0_dp-palso(j1)*ztu6(j1,ki3ec))
    pflfu(j1,ki3ec+1)= ztds1 *(palp (j1)*pflfp(j1,ki3ec+1) &
                              +palso(j1)*pflfd(j1,ki3ec+1))
    ztus1  =  ztds1*palso(j1)*ztu7(j1,ki3ec)
    ztds2  =        palso(j1)*ztu8(j1,ki3ec)
    ztds3  = 1.0_dp/(1.0_dp-palso(j1)*ztu9(j1,ki3ec)-ztds2*ztus1)
    pflcu(j1,ki3ec+1) = ztds3 *(palp (j1)*pflcp(j1,ki3ec+1) &
                               +palso(j1)*pflcd(j1,ki3ec+1) &
                               +ztds2    *pflfu(j1,ki3ec+1))
    pflfu(j1,ki3ec+1) = pflfu(j1,ki3ec+1) +ztus1*pflcu(j1,ki3ec+1)
  ENDDO
  !$acc end parallel

#ifndef _OPENACC
  IF (ldebug .AND. jb==jindex) THEN
     print *,' inv_so surface ------------------------------'
     print *,'pflfu(j1b,jb,ki3ec+1) : ',pflfu(j1b,    ki3ec+1)
     print *,'pflcu(j1b,jb,ki3ec+1) : ',pflcu(j1b,    ki3ec+1)
     print *,'palp (j1b,jb): ',palp (j1b    )
     print *,'palso(j1b,jb): ',palso(j1b    )
     print *,'ztds1               ',ztds1
     print *,'ztds2               ',ztds2
     print *,'ztds3               ',ztds3
     print *,'ztus1               ',ztus1
     print *,' .....'
  ENDIF
#endif

  ! Layer-by-layer backsubstitution

  IF (ldebug) print *,' inv_so BACKSUBSTITUTION'

!------------------------------------------------------------------------------
  !$acc parallel
  !$acc loop seq
  DO j3 = ki3ec, ki3sc, -1
    !$acc loop gang vector
    DO j1 = ki1sc, ki1ec
      pflcd(j1,j3+1) =  pflcd(j1,j3+1)                             &
                              + ztu8(j1,j3)*pflfu(j1,j3+1)         &
                              + ztu9(j1,j3)*pflcu(j1,j3+1)
      pflfd(j1,j3+1) =  pflfd(j1,j3+1)                             &
                              + ztu6(j1,j3)*pflfu(j1,j3+1)         &
                              + ztu7(j1,j3)*pflcu(j1,j3+1)
      pflcu(j1,j3  ) =  pflcu(j1,j3  )                             &
                              + ztu4(j1,j3)*pflfu(j1,j3+1)         &
                              + ztu5(j1,j3)*pflcu(j1,j3+1)
      pflfu(j1,j3  ) =  pflfu(j1,j3  )                             &
                              + ztu2(j1,j3)*pflfu(j1,j3+1)         &
                              + ztu3(j1,j3)*pflcu(j1,j3+1)         &
                              + ztu1(j1,j3)*pflcu(j1,j3)
    ENDDO

#ifndef _OPENACC
    IF (ldebug .AND. jb==jindex) THEN
      print *,' inv_so j3=',j3
      print *,'pflfu(j1b,jb,j3)  : ',pflfu(j1b    ,j3)
      print *,'pflcu(j1b,jb,j3)  : ',pflcu(j1b    ,j3)
      print *,'pflfd(j1b,jb,j3+1): ',pflfd(j1b    ,j3+1)
      print *,'pflcd(j1b,jb,j3+1): ',pflcd(j1b    ,j3+1)
    ENDIF  
#endif
  ENDDO  
  !$acc end parallel
 
!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

  !$acc end data

END SUBROUTINE inv_so

!==============================================================================
!==============================================================================
!+ Module procedure in "Radiation" 
!------------------------------------------------------------------------------

SUBROUTINE opt_th(prholwc  ,pdulwc  ,prhoiwc  ,pduiwc  ,                  &
                  paeq1    ,paeq2   ,paeq3    ,paeq4   , paeq5   ,        &
                  ki1sd    ,ki1ed   ,ki3sd    ,ki3ed   ,                  &
                  kspec    ,ki1sc   ,ki1ec    ,ki3sc   ,ki3ec    ,        &
                  ldebug   ,jindex  ,                                     &
                  podac    ,podaf   ,podsc    ,podsf   , pbsfc   ,pbsff   )

!------------------------------------------------------------------------------
!
! Description:
!
!   The module procedure opt_th calculates the optical properties of the 
!   non-gaseous constituents for one spectral interval in the thermal part 
!   of the spectrum.
!   Purpose is the calculation of (absorption and scattering) optical
!   depth and backward scattered fraction of diffuse radiation (excluding 
!   the contribution by gaseous constituents).
!
! Method:
!
! - determination of optical properies (i.e. extinction coefficient,
!   single scattering albedo and asymetry factor of the phase function)
!   using approximate relations for cloud water and cloud ice and
!   combination of five type of aerosols
!
! - calculation of optical depth (scattering and absorption) and back-
!   scattered fraction suitable for use in an implicit delta-two-stream
!   scheme
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

! Input data
! ----------
  INTEGER, INTENT (IN) ::  &
     ki1sd,       & ! start index for first  array dimension
     ki1ed,       & ! end   index for first  array dimension
     ki3sd,       & ! start index for third  array dimension
     ki3ed,       & ! end   index for third  array dimension
     kspec,       & ! selected spectral interval

   ! and the same for the computations
     ki1sc,       & ! start index for first  array computation
     ki1ec,       & ! end   index for first  array computation
     ki3sc,       & ! start index for third  array computation
     ki3ec,       & ! end   index for third  array computation
     jindex         ! index for j-loop

  LOGICAL                 , INTENT (IN) ::  &
     ldebug         ! debug control switch       

  REAL    (KIND=dp)       , INTENT (IN) ::  &
           !  Liquid and ice water density and content within for the cloudy
           !  part of each layer
     prholwc(ki1sd:ki1ed,            ki3sd:ki3ed), &
     prhoiwc(ki1sd:ki1ed,            ki3sd:ki3ed), &
     pdulwc (ki1sd:ki1ed,            ki3sd:ki3ed), &
     pduiwc (ki1sd:ki1ed,            ki3sd:ki3ed), &

           !  Aerosole contents (optical depths at 0.55 micrometer) for 5 types
     paeq1  (ki1sd:ki1ed,            ki3sd:ki3ed), &
     paeq2  (ki1sd:ki1ed,            ki3sd:ki3ed), &
     paeq3  (ki1sd:ki1ed,            ki3sd:ki3ed), &
     paeq4  (ki1sd:ki1ed,            ki3sd:ki3ed), &
     paeq5  (ki1sd:ki1ed,            ki3sd:ki3ed)   

! Output data
! -----------
  REAL    (KIND=dp)       , INTENT (OUT) ::  &
     podac (ki1sd:ki1ed,            ki3sd:ki3ed), & ! absorption optical depth
     podaf (ki1sd:ki1ed,            ki3sd:ki3ed), & ! in cloudy and free part
     podsc (ki1sd:ki1ed,            ki3sd:ki3ed), & ! scattering optical depth
     podsf (ki1sd:ki1ed,            ki3sd:ki3ed), & ! in cloudy and free part
     pbsfc (ki1sd:ki1ed,            ki3sd:ki3ed), & ! backscattering fraction 
     pbsff (ki1sd:ki1ed,            ki3sd:ki3ed)    ! in cloudy and free part

! Local parameters: 
! ----------------
  REAL    (KIND=dp)       , PARAMETER ::  &
     z1dg   = 1.0_dp/9.80665_dp, & ! 1./g
     z1d8   = 0.125_dp         , & ! 1./8 
     zepopd = 1.0E-6_dp        , & ! Security constant for optical depth
     zepssa = 1.0E-6_dp            ! Security constant for single scattering albedo

! Local scalars:
! -------------
  INTEGER  ::  &
    j1,j3,                 & ! loop indices over spatial dimensions
    ja                       ! local loop index

  REAL    (KIND=dp)        ::  &
    ! individual optical properties of liquid water and ice 
    z_lwe, z_iwe,          & ! extinction coefficient
    z_lww, z_iww,          & ! single scattering coefficient
    z_lwg, z_iwg,          & ! asymetry factor
    z_lwf, z_iwf,          & ! forward scattered fraction 
    zzg
 
!- End of header
!==============================================================================
 
!------------------------------------------------------------------------------
! Start GPU data region
!------------------------------------------------------------------------------

  !$acc data                                                              &
  !---- Argument arrays - intent(in)
  !$acc present ( prholwc,pdulwc,prhoiwc,pduiwc,paeq1,paeq2,paeq3,paeq4 ) &
  !$acc present ( paeq5                                                 ) &
  !---- Argument arrays - intent(out)
  !$acc present ( podac,podaf,podsc,podsf,pbsfc,pbsff                   ) &
  !---- Data module arrays
  !$acc present ( zaeg,zaef,zlwg,zlww,zlwe,zlwemn,zlwemx,ziwg,ziww,ziwe ) &
  !$acc present ( ziwemn,ziwemx,zaea,zaes,zaeg,zrsc                     ) &
  !---- Local automatic arrays
  !$acc present ( zlwoda,zlwods,zlwb0,ziwoda,ziwods,ziwb0,zaeoda,zaeods ) &
  !$acc present ( zaeb0                                                 )

!------------------------------------------------------------------------------
! Begin Subroutine opt_th              
!------------------------------------------------------------------------------


  IF (ldebug .AND. jb==jindex) THEN
     print *,' **** opt-th   start ********************'
     print *,' **** debug point : ',j1b,jb
  ENDIF     
 
  !$acc parallel
  !$acc loop gang vector
  DO ja=1,5
    zaef(kspec,ja)  = zaeg(kspec,ja)**2 ! forward sc.fraction f.aerosols
  ENDDO
  !$acc end parallel

#ifndef _OPENACC
  IF (ldebug .AND. jb==jindex) THEN
     DO ja=1,5
       print *,'ja, zaef(kspec,ja): ',ja,zaef(kspec,ja)
     ENDDO
  ENDIF 
#endif

  IF (ldebug) print *,' In opt-th   vor vertical loop'

! Vertical loop
! ------------- 

  DO j3 = ki3sc, ki3ec
 
    IF (ldebug) print *,' In opt-th   j3 = ',j3  

    ! Optical properties of liquid water and ice as function of the specific
    ! liquid water and ice content
 
#if defined TWOMOM_SB && defined COSMO_ART
!TS
    IF(iradpar_cloud == 1) THEN
!TS
#endif

    !$acc parallel
    !$acc loop gang vector
    DO j1 = ki1sc, ki1ec

      ! Liquid water

      z_lwg      = zlwg(1,kspec) + zlwg(2,kspec)*prholwc(j1,j3)
      z_lwg      = MAX (0.0_dp,MIN(1.0_dp,z_lwg))
      z_lwf      = z_lwg*z_lwg
      z_lww      = zlww(1,kspec) + zlww(2,kspec)*prholwc(j1,j3)
      z_lww      = MAX(zepssa,MIN(1.0_dp,z_lww)) 
      z_lwe      = z1dg * (zlwe(1,kspec) + zlwe(2,kspec)/ &
                   (zlwe(3,kspec)*prholwc(j1,j3)+zlwe(4,kspec)))
      z_lwe      = MAX(zlwemn(kspec),MIN(zlwemx(kspec),z_lwe))

      zlwoda(j1) = z_lwe *pdulwc(j1,j3) * (1.0_dp-z_lww) 
      zlwods(j1) = z_lwe *pdulwc(j1,j3) *         z_lww  * (1.0_dp-z_lwf)
      zlwb0 (j1) = z1d8*(4.0_dp+z_lwg)/(1.0_dp+z_lwg)
 
      ! Ice
 
      z_iwg      = ziwg(1,kspec) + ziwg(2,kspec)* LOG(prhoiwc(j1,j3))
      z_iwg      = MAX(0.0_dp,MIN(1.0_dp,z_iwg))
      z_iwf      = z_iwg*z_iwg
      z_iww      = ziww(1,kspec) + ziww(2,kspec)* LOG(prhoiwc(j1,j3))
      z_iww      = MAX(1.E-12_dp, MIN (1.0_dp, z_iww) )
      z_iwe      = z1dg * (ziwe(1,kspec) + ziwe(2,kspec)/  &
                   (ziwe(3,kspec)*prhoiwc(j1,j3)+ziwe(4,kspec)))
      z_iwe      = MAX(ziwemn(kspec),MIN(ziwemx(kspec),z_iwe  ))
 
      ziwoda(j1) = z_iwe * pduiwc(j1,j3)*(1.0_dp-z_iww)
      ziwods(j1) = z_iwe * pduiwc(j1,j3)*        z_iww  *(1.0_dp-z_iwf)
      ziwb0 (j1) = z1d8*(4.0_dp+z_iwg)/(1.0_dp+z_iwg)
    END DO
    !$acc end parallel

#ifndef _OPENACC
    IF (ldebug .AND. jb==jindex) THEN
      print *,' prholwc (j1b,jb) :',prholwc (j1b    ,j3)
      print *,' pdulwc  (j1b,jb) :',pdulwc  (j1b    ,j3)
      print *,' zlwoda  (j1b,jb) :',zlwoda  (j1b    )
      print *,' zlwods  (j1b,jb) :',zlwods  (j1b    )
      print *,' zlwb0   (j1b,jb) :',zlwb0   (j1b    )
      print *,' z_lwg                :',z_lwg                
      print *,' z_lwf                :',z_lwf                
      print *,' z_lww                :',z_lww                
      print *,' z_lwe                :',z_lwe                
      print *,' prhoiwc (j1b,jb) :',prhoiwc (j1b    ,j3)
      print *,' pduiwc  (j1b,jb) :',pduiwc  (j1b    ,j3)
      print *,' ziwoda  (j1b,jb) :',ziwoda  (j1b    )
      print *,' ziwods  (j1b,jb) :',ziwods  (j1b    )
      print *,' ziwb0   (j1b,jb) :',ziwb0   (j1b    )
    ENDIF   
#endif
 
#if defined TWOMOM_SB && defined COSMO_ART
!TS cloud optical properties
    ELSE

!US now it gets critical: because we removed the second index here!
      !$acc parallel
      !$acc loop gang vector
      DO j1 = ki1sc, ki1ec

        ! Use cloud optical properties based on effective radii
        zlwoda(j1) = zlwoda_th_prefac(j1,j3,kspec) * pdulwc(j1,j3)
        zlwods(j1) = zlwods_th_prefac(j1,j3,kspec) * pdulwc(j1,j3)
        zlwb0 (j1) = zlwb0_th(j1,j3,kspec)

        ziwoda(j1) = ziwoda_th_prefac(j1,j3,kspec) * pduiwc(j1,j3)
        ziwods(j1) = ziwods_th_prefac(j1,j3,kspec) * pduiwc(j1,j3)
        ziwb0 (j1) = ziwb0_th(j1,j3,kspec)

      END DO
      !$acc end parallel
    ENDIF
#endif

    ! Aerosoles

#ifdef COSMOART
    IF(.NOT. l_cosmo_art) THEN
#endif
 
      !$acc parallel
      !$acc loop gang vector
      DO j1 = ki1sc, ki1ec
        zaeoda(j1) =  paeq1(j1,j3)*zaea(kspec,1) &
                    + paeq2(j1,j3)*zaea(kspec,2) &
                    + paeq3(j1,j3)*zaea(kspec,3) &
                    + paeq4(j1,j3)*zaea(kspec,4) &
                    + paeq5(j1,j3)*zaea(kspec,5)
        zaeods(j1) =                                                    &
               ( paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1)) ) &
              +( paeq2(j1,j3)*zaes(kspec,2)*(1.0_dp-zaef(kspec,2)) ) &
              +( paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3)) ) &
              +( paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4)) ) &
              +( paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5)) )
        zzg=                                                                   &
            ((paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1)))*zaeg(kspec,1)  &
            +(paeq2(j1,j3)*zaes(kspec,2)*(1.0_dp-zaef(kspec,2)))*zaeg(kspec,2)  &
            +(paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3)))*zaeg(kspec,3)  &
            +(paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4)))*zaeg(kspec,4)  &
            +(paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5)))*zaeg(kspec,5)) &
            / MAX( zaeods(j1),zepopd)
        zaeb0 (j1) = z1d8*(4.0_dp+zzg)/(1.0_dp+zzg)
      ENDDO
      !$acc end parallel

#ifdef COSMOART
!US also here: 2nd index removed!!
    ELSE   ! l_cosmo_art
      IF ((.NOT.lrad_dust) .AND. (.NOT. lrad_seas) .AND. (.NOT. lrad_aero)) THEN
        !$acc parallel
        !$acc loop gang vector
        DO j1 = ki1sc, ki1ec
          zaeoda(j1) =  paeq1(j1,j3)*zaea(kspec,1) &
                      + paeq2(j1,j3)*zaea(kspec,2) &
                      + paeq3(j1,j3)*zaea(kspec,3) &
                      + paeq4(j1,j3)*zaea(kspec,4) &
                      + paeq5(j1,j3)*zaea(kspec,5)
          zaeods(j1) =                                                    &
                    ( paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1)) ) &
                   +( paeq2(j1,j3)*zaes(kspec,2)*(1.0_dp-zaef(kspec,2)) ) &
                   +( paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3)) ) &
                   +( paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4)) ) &
                   +( paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5)) )
          zzg=                                                                    &
              ((paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1)))*zaeg(kspec,1)  &
              +(paeq2(j1,j3)*zaes(kspec,2)*(1.0_dp-zaef(kspec,2)))*zaeg(kspec,2)  &
              +(paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3)))*zaeg(kspec,3)  &
              +(paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4)))*zaeg(kspec,4)  &
              +(paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5)))*zaeg(kspec,5)) &
              / MAX( zaeods(j1),zepopd)
          zaeb0 (j1) = z1d8*(4.0_dp+zzg)/(1.0_dp+zzg)
        ENDDO
        !$acc end parallel

      ELSEIF ((lrad_dust) .AND. (.NOT. lrad_seas) .AND. (.NOT. lrad_aero)) THEN

        !$acc parallel
        !$acc loop gang vector
        DO j1 = ki1sc, ki1ec
          zaeoda(j1) =   paeq1(j1,j3)*zaea(kspec,1) + tau_abs_dust(j1,j3,kspec)   &
                       + paeq2(j1,j3)*zaea(kspec,2)                               &
                       + paeq3(j1,j3)*zaea(kspec,3)                               &
                       + paeq4(j1,j3)*zaea(kspec,4)                               &
                       + paeq5(j1,j3)*zaea(kspec,5)

          zaeods(j1) = paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1))+             &
                  tau_scat_dust(j1,j3,kspec)*(1.0_dp-(asym_dust(j1,j3,kspec)**2))     &
                     + paeq2(j1,j3)*zaes(kspec,2)*(1.0_dp-zaef(kspec,2))              &
                     + paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3))              &
                     + paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))              &
                     + paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))

          zzg=( paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1))*zaeg(kspec,1) +       &
                       tau_scat_dust(j1,j3,kspec)*                                      &
                       (1.0_dp-(asym_dust(j1,j3,kspec)**2))*asym_dust(j1,j3,kspec)      &
                     + paeq2(j1,j3)*zaes(kspec,2)*(1.0_dp-zaef(kspec,2))*zaeg(kspec,2)  &
                     + paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3))*zaeg(kspec,3)  &
                     + paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))*zaeg(kspec,4)  &
                     + paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))*zaeg(kspec,5)) &
                          / MAX( zaeods(j1),zepopd)
          zaeb0 (j1) = z1d8*(4.0_dp+zzg)/(1.0_dp+zzg)
        ENDDO
        !$acc end parallel

      ELSEIF ((lrad_dust) .AND. (.NOT. lrad_seas) .AND. (lrad_aero)) THEN

        !$acc parallel
        !$acc loop gang vector
        DO j1 = ki1sc, ki1ec
          zaeoda(j1) = paeq1(j1,j3)*zaea(kspec,1) + tau_abs_dust(j1,j3,kspec) &
                + paeq2(j1,j3)*zaea(kspec,2)                                  &
                + paeq3(j1,j3)*zaea(kspec,3) + tau_abs_aero(j1,j3,kspec)      &
                + paeq4(j1,j3)*zaea(kspec,4)                                  &
                + paeq5(j1,j3)*zaea(kspec,5)

          zaeods(j1) = paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1)) +           &
                  tau_scat_dust(j1,j3,kspec)*(1.0_dp-(asym_dust(j1,j3,kspec)**2))    &
                 +paeq2(j1,j3)*zaes(kspec,2)*(1.0_dp-zaef(kspec,2))                  &
                 +paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3)) +                &
                  tau_scat_aero(j1,j3,kspec)*(1.0_dp-(asym_aero(j1,j3,kspec)**2))    &
                 +paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))                  &
                 +paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))

          zzg= ( paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1))*zaeg(kspec,1) +                          & 
               tau_scat_dust(j1,j3,kspec)*(1.0_dp-(asym_dust(j1,j3,kspec)**2))*asym_dust(j1,j3,kspec)       &
             + paeq2(j1,j3)*zaes(kspec,2)*(1.0_dp-zaef(kspec,2))*zaeg(kspec,2)                              &
             + paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3))*zaeg(kspec,3)  +                           &
               tau_scat_aero(j1,j3,kspec)*(1.0_dp-(asym_aero(j1,j3,kspec)**2.0_dp))*asym_aero(j1,j3,kspec)  &
             + paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))*zaeg(kspec,4)                              &
             + paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))*zaeg(kspec,5) )                            &
               / MAX( zaeods(j1),zepopd)
          zaeb0 (j1) = z1d8*(4.0_dp+zzg)/(1.0_dp+zzg)
        ENDDO
        !$acc end parallel

      ELSEIF ((.NOT. lrad_dust) .AND. (.NOT. lrad_seas) .AND. (lrad_aero)) THEN

        !$acc parallel
        !$acc loop gang vector
        DO j1 = ki1sc, ki1ec
          zaeoda(j1) = paeq1(j1,j3)*zaea(kspec,1)                         &
                + paeq2(j1,j3)*zaea(kspec,2)                              &
                + paeq3(j1,j3)*zaea(kspec,3) + tau_abs_aero(j1,j3,kspec)  &
                + paeq4(j1,j3)*zaea(kspec,4)                              &
                + paeq5(j1,j3)*zaea(kspec,5)

          zaeods(j1) = paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1))                &
                 +paeq2(j1,j3)*zaes(kspec,2)*(1.0_dp-zaef(kspec,2))                     &
                 +paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3))+                    &
                  tau_scat_aero(j1,j3,kspec)*(1.0_dp-(asym_aero(j1,j3,kspec)**2.0_dp))  &
                 +paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))                     &
                 +paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))

          zzg=( paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1))*zaeg(kspec,1)                              &
               +paeq2(j1,j3)*zaes(kspec,2)*(1.0_dp-zaef(kspec,2))*zaeg(kspec,2)                              &
               +paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3))*zaeg(kspec,3)  +                           &
                tau_scat_aero(j1,j3,kspec)*(1.0_dp-(asym_aero(j1,j3,kspec)**2.0_dp))*asym_aero(j1,j3,kspec)  &
               +paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))*zaeg(kspec,4)                              &
               +paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))*zaeg(kspec,5) )                            &
                   / MAX( zaeods(j1),zepopd)
          zaeb0 (j1) = z1d8*(4.0_dp+zzg)/(1.0_dp+zzg)
        ENDDO
        !$acc end parallel

      ELSEIF ((lrad_dust) .AND. (lrad_seas) .AND. (.NOT.lrad_aero)) THEN

        !$acc parallel
        !$acc loop gang vector
        DO j1 = ki1sc, ki1ec
          zaeoda(j1) =   paeq1(j1,j3)*zaea(kspec,1) + tau_abs_dust(j1,j3,kspec)   &
                  + tau_abs_seas(j1,j3,kspec)                                     &
                  + paeq3(j1,j3)*zaea(kspec,3)                                    &
                  + paeq4(j1,j3)*zaea(kspec,4)                                    &
                  + paeq5(j1,j3)*zaea(kspec,5)

          zaeods(j1) = paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1))+             &
                  tau_scat_dust(j1,j3,kspec)*(1.0_dp-(asym_dust(j1,j3,kspec)**2))     &
                + tau_scat_seas(j1,j3,kspec)*(1.0_dp-(asym_seas(j1,j3,kspec)**2))     &
                + paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3))                   &
                + paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))                   &
                + paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))

          zzg=( paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1))*zaeg(kspec,1) +                       &
                tau_scat_dust(j1,j3,kspec)*(1.0_dp-(asym_dust(j1,j3,kspec)**2))*asym_dust(j1,j3,kspec)  &
               +tau_scat_seas(j1,j3,kspec)*(1.0_dp-(asym_seas(j1,j3,kspec)**2))*asym_seas(j1,j3,kspec)  &
               +paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3))*zaeg(kspec,3)                         &
               +paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))*zaeg(kspec,4)                         &
               +paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))*zaeg(kspec,5) )                       &
             / MAX( zaeods(j1),zepopd)
          zaeb0 (j1) = z1d8*(4.0_dp+zzg)/(1.0_dp+zzg)
        ENDDO
        !$acc end parallel

      ELSEIF ((.NOT. lrad_dust) .and. (lrad_seas) .and. (.NOT. lrad_aero)) THEN

        !$acc parallel
        !$acc loop gang vector
        DO j1 = ki1sc, ki1ec
          zaeoda(j1) = paeq1(j1,j3)*zaea(kspec,1)       &
                 + tau_abs_seas(j1,j3,kspec)            &
                 + paeq3(j1,j3)*zaea(kspec,3)           &
                 + paeq4(j1,j3)*zaea(kspec,4)           &
                 + paeq5(j1,j3)*zaea(kspec,5)

          zaeods(j1) =  paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1))             &
                  + tau_scat_seas(j1,j3,kspec)*(1.0_dp-(asym_seas(j1,j3,kspec)**2))  &
                  + paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3))                 &
                  + paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))                 &
                  + paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))

          zzg=   ( paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1))*zaeg(kspec,1)                  &
           +tau_scat_seas(j1,j3,kspec)*(1.0_dp-(asym_seas(j1,j3,kspec)**2))*asym_seas(j1,j3,kspec) &
           +paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3))*zaeg(kspec,3)                         &
           +paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))*zaeg(kspec,4)                         &
           +paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))*zaeg(kspec,5) )                       &
           / MAX( zaeods(j1),zepopd)
          zaeb0 (j1) = z1d8*(4.0_dp+zzg)/(1.0_dp+zzg)
        ENDDO
        !$acc end parallel

      ELSEIF ((lrad_dust) .AND. (lrad_seas) .AND. (lrad_aero)) THEN

        !$acc parallel
        !$acc loop gang vector
        DO j1 = ki1sc, ki1ec
          zaeoda(j1) = paeq1(j1,j3)*zaea(kspec,1) + tau_abs_dust(j1,j3,kspec) &
                + tau_abs_seas(j1,j3,kspec)                                   &
                + paeq3(j1,j3)*zaea(kspec,3) + tau_abs_aero(j1,j3,kspec)      &
                + paeq4(j1,j3)*zaea(kspec,4)                                  &
                + paeq5(j1,j3)*zaea(kspec,5)

          zaeods(j1) = paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1)) +          &
                  tau_scat_dust(j1,j3,kspec)*(1.0_dp-(asym_dust(j1,j3,kspec)**2))   &
                 +tau_scat_seas(j1,j3,kspec)*(1.0_dp-(asym_seas(j1,j3,kspec)**2))   &
                 +paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3)) +               &
                  tau_scat_aero(j1,j3,kspec)*(1.0_dp-(asym_aero(j1,j3,kspec)**2))   &
                 +paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))                 &
                 +paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))

          zzg= ( paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1))*zaeg(kspec,1) +                          &
               tau_scat_dust(j1,j3,kspec)*(1.0_dp-(asym_dust(j1,j3,kspec)**2))*asym_dust(j1,j3,kspec)       &
              +tau_scat_seas(j1,j3,kspec)*(1.0_dp-(asym_seas(j1,j3,kspec)**2))*asym_seas(j1,j3,kspec)       &
              +paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3))*zaeg(kspec,3)  +                           &
               tau_scat_aero(j1,j3,kspec)*(1.0_dp-(asym_aero(j1,j3,kspec)**2.0_dp))*asym_aero(j1,j3,kspec)  &
              +paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))*zaeg(kspec,4)                              &
              +paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))*zaeg(kspec,5) )                            &
             / MAX( zaeods(j1),zepopd)
          zaeb0 (j1) = z1d8*(4.0_dp+zzg)/(1.0_dp+zzg)
        ENDDO
        !$acc end parallel

      ELSEIF ((.NOT. lrad_dust) .AND. (lrad_seas) .AND. (lrad_aero)) THEN

        !$acc parallel
        !$acc loop gang vector
        DO j1 = ki1sc, ki1ec
          zaeoda(j1) = paeq1(j1,j3)*zaea(kspec,1)                         &
                + tau_abs_seas(j1,j3,kspec)                               &
                + paeq3(j1,j3)*zaea(kspec,3) + tau_abs_aero(j1,j3,kspec)  &
                + paeq4(j1,j3)*zaea(kspec,4)                              &
                + paeq5(j1,j3)*zaea(kspec,5)

          zaeods(j1) = paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1))                &
                 +tau_scat_seas(j1,j3,kspec)*(1.0_dp-(asym_seas(j1,j3,kspec)**2.0_dp))  &
                 +paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3)) +                   &
                  tau_scat_aero(j1,j3,kspec)*(1.0_dp-(asym_aero(j1,j3,kspec)**2.0_dp))  &
                 +paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))                     &
                 +paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))

          zzg=( paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1))*zaeg(kspec,1)                            &
             +tau_scat_seas(j1,j3,kspec)*(1.0_dp-(asym_seas(j1,j3,kspec)**2.0_dp))*asym_seas(j1,j3,kspec)  &
             +paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3))*zaeg(kspec,3) +                            &
              tau_scat_aero(j1,j3,kspec)*(1.0_dp-(asym_aero(j1,j3,kspec)**2.0_dp))*asym_aero(j1,j3,kspec)  &
             +paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))*zaeg(kspec,4)                              &
             +paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))*zaeg(kspec,5) )                            &
                 / MAX( zaeods(j1),zepopd)
          zaeb0 (j1) = z1d8*(4.0_dp+zzg)/(1.0_dp+zzg)
        ENDDO
        !$acc end parallel

      ENDIF
    ENDIF
#endif

#ifndef _OPENACC
    IF (ldebug .AND. jb==jindex) THEN
      print *,' zaeoda  (j1b,jb) :',zaeoda  (j1b    )
      print *,' zaeods  (j1b,jb) :',zaeods  (j1b    )
      print *,' zaeb0   (j1b,jb) :',zaeb0   (j1b    )
    ENDIF
#endif
 
    ! Combined effects
 
    ! a) cloud free part of each layer

    !$acc parallel
    !$acc loop gang vector
    DO j1 = ki1sc, ki1ec
      podaf(j1,j3) = MAX(zaeoda(j1), zepopd)   
      podsf(j1,j3) = MAX(zaeods(j1), zepopd)
      pbsff(j1,j3) =     zaeb0 (j1)
    ENDDO
    !$acc end parallel

#ifndef _OPENACC
    IF (ldebug .AND. jb==jindex) THEN
      print *,' podaf   (j1b,jb,j3) :',podaf   (j1b    ,j3)
      print *,' podsf   (j1b,jb,j3) :',podsf   (j1b    ,j3)
      print *,' pbsff   (j1b,jb,j3) :',pbsff   (j1b    ,j3)
    ENDIF
#endif
 
    ! b) cloudy part of each layer

    !$acc parallel
    !$acc loop gang vector
    DO j1 = ki1sc, ki1ec
      podac(j1,j3) = MAX( zlwoda (j1) + ziwoda (j1) + zaeoda(j1), zepopd)
      podsc(j1,j3) = MAX( zlwods (j1) + ziwods (j1) + zaeods(j1), zepopd)
      podsc(j1,j3) = MIN( podsc(j1,j3), (1.0_wp-zepssa) * (podac(j1,j3)+podsc(j1,j3)))
      pbsfc(j1,j3) = (  zlwb0 (j1)*zlwods (j1)                      &
                      + ziwb0 (j1)*ziwods (j1)                      &
                      + zaeb0 (j1)*zaeods (j1) ) / podsc(j1,j3)
    ENDDO
    !$acc end parallel

#ifndef _OPENACC
    IF (ldebug .AND. jb==jindex) THEN
      print *,' podac   (j1b,jb,j3) :',podac   (j1b    ,j3)
      print *,' podsc   (j1b,jb,j3) :',podsc   (j1b    ,j3)
      print *,' pbsfc   (j1b,jb,j3) :',pbsfc   (j1b    ,j3)
    ENDIF 
#endif
 
  ! End of vertical loop
  ! --------------------
  ENDDO       
 
!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

  !$acc end data

END SUBROUTINE opt_th

!==============================================================================
!==============================================================================
!+ Module procedure in "Radiation" 
!------------------------------------------------------------------------------

SUBROUTINE opt_so(prholwc  ,pdulwc  ,prhoiwc  ,pduiwc  ,                    &
                  paeq1    ,paeq2   ,paeq3    ,paeq4   , paeq5   ,          &
                  pdp      ,papra   ,psmu0    ,pqsmu0  ,                    &
                  ki1sd    ,ki1ed   ,ki3sd    ,ki3ed   ,                    &
                  kspec    ,ki1sc   ,ki1ec    ,ki3sc   ,ki3ec    ,          &
                  ldebug   ,jindex  ,                                       &
                  podac    ,podaf   ,podsc    ,podsf   , pbsfc   ,pbsff   , &
                  pusfc    ,pusff                                           )

!------------------------------------------------------------------------------
!
! Description:
!
!   The module procedure opt_so calculates the optical properties of the 
!   non-gaseous constituents for one spectral interval in the solar part 
!   of the spectrum.
!   Purpose is the calculation of (absorption and scattering) optical
!   depth and backward scattered fraction for diffuse and upward scattered 
!   fraction for direct solar radiation (excluding the contribution by 
!   gaseous constituents).
!
! Method:
!
! - determination of optical properies (i.e. extinction coefficient,
!   single scattering albedo and asymetry factor of the phase function)
!   using approximate relations for Rayleigh effects, cloud water,
!   cloud ice and a combination of five type of aerosols
!
! - calculation of optical depth (scattering and absorption), back-
!   scattered and upscattered fraction of radiation suitable for use
!   in an implicit delta-two-stream scheme
!     
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

! Input data
! ----------
  INTEGER, INTENT (IN) ::  &
     ki1sd,       & ! start index for first  array dimension
     ki1ed,       & ! end   index for first  array dimension
     ki3sd,       & ! start index for third  array dimension
     ki3ed,       & ! end   index for third  array dimension

   ! and the same for the computations
     ki1sc,       & ! start index for first  array computation
     ki1ec,       & ! end   index for first  array computation
     ki3sc,       & ! start index for third  array computation
     ki3ec,       & ! end   index for third  array computation
     kspec,       & ! selected spectral interval
     jindex         ! index of j-loop

  LOGICAL                 , INTENT (IN) ::  &
     ldebug         ! debug control switch       

  REAL    (KIND=dp)       , INTENT (IN) ::  &
           !  Liquid and ice water density and content within for the cloudy
           !  part of each layer
     prholwc(ki1sd:ki1ed,            ki3sd:ki3ed), &
     prhoiwc(ki1sd:ki1ed,            ki3sd:ki3ed), &
     pdulwc (ki1sd:ki1ed,            ki3sd:ki3ed), &
     pduiwc (ki1sd:ki1ed,            ki3sd:ki3ed), &

           !  Aerosole contents (optical depths at 0.55 micrometer) for 5 types
     paeq1  (ki1sd:ki1ed,            ki3sd:ki3ed), &
     paeq2  (ki1sd:ki1ed,            ki3sd:ki3ed), &
     paeq3  (ki1sd:ki1ed,            ki3sd:ki3ed), &
     paeq4  (ki1sd:ki1ed,            ki3sd:ki3ed), &
     paeq5  (ki1sd:ki1ed,            ki3sd:ki3ed), &

           !  pressure thickness of layers 
     pdp    (ki1sd:ki1ed,            ki3sd:ki3ed), &
 
           !  zenith angle and it's inverse
     psmu0  (ki1sd:ki1ed)                        , &
     pqsmu0 (ki1sd:ki1ed)             

  REAL    (KIND=dp)       , INTENT (INOUT) ::  &
           ! mean layer pressure (TOA on input)
     papra  (ki1sd:ki1ed)                            
  
! Output data
! -----------
  REAL    (KIND=dp)       , INTENT (OUT) ::  &
     podac (ki1sd:ki1ed,            ki3sd:ki3ed), & ! absorption optical depth
     podaf (ki1sd:ki1ed,            ki3sd:ki3ed), & ! in cloudy and free part
     podsc (ki1sd:ki1ed,            ki3sd:ki3ed), & ! scattering optical depth
     podsf (ki1sd:ki1ed,            ki3sd:ki3ed), & ! in cloudy and free part
     pbsfc (ki1sd:ki1ed,            ki3sd:ki3ed), & ! backscattering fraction 
     pbsff (ki1sd:ki1ed,            ki3sd:ki3ed), & ! in cloudy and free part
     pusfc (ki1sd:ki1ed,            ki3sd:ki3ed), & ! upward scattered fraction
     pusff (ki1sd:ki1ed,            ki3sd:ki3ed)    ! in cloudy and free part

! Local parameters: 
! ----------------
  REAL    (KIND=dp)       , PARAMETER ::  &
     z1dg   = 1.0_dp/9.80665_dp, & ! 1./g
     z1d8   = 0.125_dp         , & ! 1./8 
     zepopd = 1.0E-6_dp        , & ! Security constant for optical depth
     zepssa = 1.0E-6_dp            ! Security constant for single scattering albedo

! Local scalars:
! -------------
  INTEGER  ::  &
    j1,j2,j3,              & ! loop indices over spatial dimensions
    ja      ,              & ! local loop index
    isp                      ! (=kspec, but shorter for notation purposes)

  REAL    (KIND=dp)        ::  &
    ! individual optical properties of liquid water and ice 
    z_lwe, z_iwe,          & ! extinction coefficient
    z_lww, z_iww,          & ! single scattering coefficient
    z_lwg, z_iwg,          & ! asymetry factor
    z_lwf, z_iwf,          & ! forward scattered fraction 
    zzg
 
!- End of header
!==============================================================================
 
!------------------------------------------------------------------------------
! Start GPU data region
!------------------------------------------------------------------------------

  !$acc data &
  !---- Argument arrays - intent(in)
  !$acc present ( prholwc,pdulwc,prhoiwc,pduiwc,paeq1,paeq2,paeq3,paeq4 ) &
  !$acc present ( paeq5,pdp,psmu0,pqsmu0                                ) &
  !---- Argument arrays - intent(inout)
  !$acc present ( papra                                                 ) &
  !---- Argument arrays - intent(out)
  !$acc present ( podac,podaf,podsc,podsf,pbsfc,pbsff,pusfc,pusff       ) &
  !---- Data module arrays
  !$acc present ( zaeg,zaef,zlwg,zlww,zlwe,zlwemn,zlwemx,ziwg,ziww,ziwe ) &
  !$acc present ( ziwemn,ziwemx,zaea,zaes                               ) &
  !---- Local automatic arrays
  !$acc present ( zlwoda,zlwods,zlwb0,zlwb,ziwoda,ziwods,ziwb0,ziwb     ) &
  !$acc present ( zaeoda,zaeods,zaeb0,zaeb,zraods                       )

!------------------------------------------------------------------------------
! Begin Subroutine opt_so              
!------------------------------------------------------------------------------

  isp      = kspec
#ifndef _OPENACC
  IF (ldebug .AND. jb==jindex) THEN
      print *,' **** opt-so   start **********************'
      print *,' **** debug point : ',j1b,jb
      print *,' **** interval    : ',isp    
  ENDIF 
#endif
 
  !$acc parallel
  !$acc loop gang vector
  DO ja=1,5
    zaef(isp,ja)  = zaeg(isp,ja)**2 ! forward sc.fraction f.aerosols
  ENDDO
  !$acc end parallel

#ifndef _OPENACC
  IF (ldebug) THEN
     DO ja=1,5
       print *,'ja, zaef(isp,ja): ',ja,zaef(isp,ja)
     ENDDO 
  ENDIF 
#endif

  IF (ldebug) print *,' In opt-so   vor vertical loop'
 
! Vertical loop
! ------------ 

  DO j3=ki3sc,ki3ec
 
    if (ldebug) print *,' In opt-so   j3 = ',j3  

    !     Optical properties of liquid water and ice as function of the specific
    !     liquid water and ice content
 
#if defined TWOMOM_SB && defined COSMO_ART
    IF(iradpar_cloud == 1) THEN
#endif

      !$acc parallel
      !$acc loop gang vector
      DO j1=ki1sc,ki1ec
 
        ! liquid water effects
 
        z_lwg = zlwg(1,isp) + zlwg(2,isp)*prholwc(j1,j3)
        z_lwg = MAX(0.0_dp,MIN(1.0_dp,z_lwg))
        z_lwf = z_lwg*z_lwg 
        z_lww = zlww(1,isp) + zlww(2,isp)*prholwc(j1,j3)
        z_lww = MAX(zepssa,MIN(1.0_dp-zepssa,z_lww))
        z_lwe = z1dg * (zlwe(1,isp) + zlwe(2,isp)/ &
                (zlwe(3,isp)*prholwc(j1,j3)+zlwe(4,isp)))
        z_lwe = MAX(zlwemn(isp),MIN(zlwemx(isp),z_lwe))
 
        zlwoda(j1) = z_lwe*pdulwc(j1,j3)*(1.0_dp-z_lww)
        zlwods(j1) = z_lwe*pdulwc(j1,j3)*        z_lww *(1.0_dp-z_lwf)
        zlwb0 (j1) = z1d8*(4.0_dp+z_lwg)/(1.0_dp+z_lwg) 
        zlwb  (j1) = 0.5_dp-0.75_dp*psmu0(j1)*z_lwg/(1.0_dp+z_lwg)
 
        ! ice water effects
 
        z_iwg = ziwg(1,isp) + ziwg(2,isp)* LOG(prhoiwc(j1,j3))
        z_iwg = MAX(0.0_dp,MIN(1.0_dp,z_iwg))
        z_iwf = z_iwg*z_iwg
        z_iww = ziww(1,isp) + ziww(2,isp)* LOG(prhoiwc(j1,j3))
        z_iww = MAX(zepssa,MIN(1.0_wp-zepssa,z_iww))
        z_iwe = z1dg * (ziwe(1,isp) + ziwe(2,isp)/ &
                   (ziwe(3,isp)*prhoiwc(j1,j3)+ziwe(4,isp)))
        z_iwe = MAX(ziwemn(isp),MIN(ziwemx(isp),z_iwe  ))
 
        ziwoda(j1) = z_iwe*pduiwc(j1,j3)*(1.0_dp-z_iww)
        ziwods(j1) = z_iwe*pduiwc(j1,j3)*        z_iww *(1.0_dp-z_iwf)
        ziwb0 (j1) = z1d8*(4.0_dp+z_iwg)/(1.0_dp+z_iwg)
        ziwb  (j1) = 0.5_dp-0.75_dp*psmu0(j1)*z_iwg/(1.0_dp+z_iwg)
      END DO
      !$acc end parallel

#ifndef _OPENACC
      IF (ldebug .AND. jb==jindex) THEN
        print *,' prholwc (j1b,jb) :',prholwc (j1b    ,j3)   
        print *,' pdulwc  (j1b,jb) :',pdulwc  (j1b    ,j3)
        print *,' zlwoda  (j1b,jb) :',zlwoda  (j1b    )
        print *,' zlwods  (j1b,jb) :',zlwods  (j1b    )
        print *,' zlwb0   (j1b,jb) :',zlwb0   (j1b    )
        print *,' zlwb    (j1b,jb) :',zlwb    (j1b    )
        print *,' z_lwg            :',z_lwg                
        print *,' z_lwf            :',z_lwf                
        print *,' z_lww            :',z_lww                
        print *,' z_lwe            :',z_lwe                
        print *,' prhoiwc (j1b,jb) :',prhoiwc (j1b    ,j3)
        print *,' pduiwc  (j1b,jb) :',pduiwc  (j1b    ,j3)
        print *,' ziwoda  (j1b,jb) :',ziwoda  (j1b    )
        print *,' ziwods  (j1b,jb) :',ziwods  (j1b    )
        print *,' ziwb0   (j1b,jb) :',ziwb0   (j1b    )
        print *,' ziwb    (j1b,jb) :',ziwb    (j1b    )
      ENDIF 
#endif

#if defined TWOMOM_SB && defined COSMO_ART
!US Attention again: 2nd index removed
    ELSE
      !$acc parallel
      !$acc loop gang vector
      DO j1=ki1sc,ki1ec

        ! Use cloud optical properties based on effective radii
        zlwoda(j1) = zlwoda_so_prefac(j1,j3,kspec)*pdulwc(j1,j3)
        zlwods(j1) = zlwods_so_prefac(j1,j3,kspec)*pdulwc(j1,j3)
        zlwb0 (j1) = zlwb0_so(j1,j3,kspec)
        zlwb  (j1) = zlwb_so_prefac(j1,j3,kspec)*psmu0(j1)

        ziwoda(j1) = ziwoda_so_prefac(j1,j3,kspec)*pduiwc(j1,j3)
        ziwods(j1) = ziwods_so_prefac(j1,j3,kspec)*pduiwc(j1,j3)
        ziwb0 (j1) = ziwb0_so(j1,j3,kspec)
        ziwb  (j1) = ziwb_so_prefac(j1,j3,kspec)*psmu0(j1)

      END DO
      !$acc end parallel
    ENDIF
#endif

    ! Optical properties of five aerosol types combined
 
#ifdef COSMOART
    IF(.NOT. l_cosmo_art) THEN
#endif

      !$acc parallel
      !$acc loop gang vector
      DO j1=ki1sc,ki1ec

        zaeoda(j1)=                                          &
          paeq1(j1,j3)*zaea(isp,1)+paeq2(j1,j3)*zaea(isp,2)+ &
          paeq3(j1,j3)*zaea(isp,3)+paeq4(j1,j3)*zaea(isp,4)+ &
          paeq5(j1,j3)*zaea(isp,5)

        zaeods(j1)=                                            &
            ( paeq1(j1,j3)*zaes(isp,1)*(1.0_dp-zaef(isp,1)) )  &
           +( paeq2(j1,j3)*zaes(isp,2)*(1.0_dp-zaef(isp,2)) )  &
           +( paeq3(j1,j3)*zaes(isp,3)*(1.0_dp-zaef(isp,3)) )  &
           +( paeq4(j1,j3)*zaes(isp,4)*(1.0_dp-zaef(isp,4)) )  &
           +( paeq5(j1,j3)*zaes(isp,5)*(1.0_dp-zaef(isp,5)) )
 
        zzg=(                                                           &
           (paeq1(j1,j3)*zaes(isp,1)*(1.0_dp-zaef(isp,1)))*zaeg(isp,1)  &
          +(paeq2(j1,j3)*zaes(isp,2)*(1.0_dp-zaef(isp,2)))*zaeg(isp,2)  &
          +(paeq3(j1,j3)*zaes(isp,3)*(1.0_dp-zaef(isp,3)))*zaeg(isp,3)  &
          +(paeq4(j1,j3)*zaes(isp,4)*(1.0_dp-zaef(isp,4)))*zaeg(isp,4)  &
          +(paeq5(j1,j3)*zaes(isp,5)*(1.0_dp-zaef(isp,5)))*zaeg(isp,5)) &
          / MAX( zaeods(j1),zepopd)

        zaeb0(j1) = z1d8*(4.0_dp+zzg)/(1.0_dp+zzg)
        zaeb (j1) = 0.5_dp-0.75_dp*psmu0(j1)*zzg/(1.0_dp+zzg)

      ENDDO
      !$acc end parallel

#ifdef COSMOART
    ELSE   ! l_cosmo_art
!US Attention again: 2nd index removed
      IF ((.NOT.lrad_dust) .AND. (.NOT. lrad_seas) .AND. (.NOT.lrad_aero)) THEN
        !$acc parallel
        !$acc loop gang vector
        DO j1=ki1sc,ki1ec

          zaeoda(j1)=                                          &
            paeq1(j1,j3)*zaea(isp,1)+paeq2(j1,j3)*zaea(isp,2)+ &
            paeq3(j1,j3)*zaea(isp,3)+paeq4(j1,j3)*zaea(isp,4)+ &
            paeq5(j1,j3)*zaea(isp,5)

          zaeods(j1)=                                            &
              ( paeq1(j1,j3)*zaes(isp,1)*(1.0_dp-zaef(isp,1)) )  &
             +( paeq2(j1,j3)*zaes(isp,2)*(1.0_dp-zaef(isp,2)) )  &
             +( paeq3(j1,j3)*zaes(isp,3)*(1.0_dp-zaef(isp,3)) )  &
             +( paeq4(j1,j3)*zaes(isp,4)*(1.0_dp-zaef(isp,4)) )  &
             +( paeq5(j1,j3)*zaes(isp,5)*(1.0_dp-zaef(isp,5)) )

          zzg=(                                                           &
             (paeq1(j1,j3)*zaes(isp,1)*(1.0_dp-zaef(isp,1)))*zaeg(isp,1)  &
            +(paeq2(j1,j3)*zaes(isp,2)*(1.0_dp-zaef(isp,2)))*zaeg(isp,2)  &
            +(paeq3(j1,j3)*zaes(isp,3)*(1.0_dp-zaef(isp,3)))*zaeg(isp,3)  &
            +(paeq4(j1,j3)*zaes(isp,4)*(1.0_dp-zaef(isp,4)))*zaeg(isp,4)  &
            +(paeq5(j1,j3)*zaes(isp,5)*(1.0_dp-zaef(isp,5)))*zaeg(isp,5)) &
            / MAX( zaeods(j1),zepopd)

          zaeb0(j1) = z1d8*(4.0_dp+zzg)/(1.0_dp+zzg)
          zaeb (j1) = 0.5_dp-0.75_dp*psmu0(j1)*zzg/(1.0_dp+zzg)

        ENDDO
        !$acc end parallel

      ELSEIF ((lrad_dust) .AND. (.NOT.lrad_seas) .AND. (.NOT. lrad_aero)) THEN    

        !$acc parallel
        !$acc loop gang vector
        DO j1 = ki1sc, ki1ec
          zaeoda(j1) =   paeq1(j1,j3)*zaea(kspec,1) + tau_abs_dust(j1,j3,kspec) &
                  + paeq2(j1,j3)*zaea(kspec,2)                                  &
                  + paeq3(j1,j3)*zaea(kspec,3)                                  &
                  + paeq4(j1,j3)*zaea(kspec,4)                                  &
                  + paeq5(j1,j3)*zaea(kspec,5)

          zaeods(j1) = paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1))+         &
                  tau_scat_dust(j1,j3,kspec)*(1.0_dp-(asym_dust(j1,j3,kspec)**2)) &
                + paeq2(j1,j3)*zaes(kspec,2)*(1.0_dp-zaef(kspec,2))               &
                + paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3))               &
                + paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))               &
                + paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))

          zzg=( paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1))*zaeg(kspec,1) +                     &
                tau_scat_dust(j1,j3,kspec)*(1.0_dp-(asym_dust(j1,j3,kspec)**2))*asym_dust(j1,j3,kspec)&
               +paeq2(j1,j3)*zaes(kspec,2)*(1.0_dp-zaef(kspec,2))*zaeg(kspec,2)                       &
               +paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3))*zaeg(kspec,3)                       &
               +paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))*zaeg(kspec,4)                       &
               +paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))*zaeg(kspec,5) )                     &
             / MAX( zaeods(j1),zepopd)
          zaeb0 (j1) = z1d8*(4.0_dp+zzg)/(1.0_dp+zzg)
          zaeb  (j1) = 0.5_dp-0.75_dp*psmu0(j1)*zzg/(1.0_dp+zzg)
        ENDDO
        !$acc end parallel

      ELSEIF ((lrad_dust) .AND. (.NOT. lrad_seas) .AND. (lrad_aero)) THEN

        !$acc parallel
        !$acc loop gang vector
        DO j1 = ki1sc, ki1ec
          zaeoda(j1) = paeq1(j1,j3)*zaea(kspec,1) + tau_abs_dust(j1,j3,kspec) &
                + paeq2(j1,j3)*zaea(kspec,2)                                  &
                + paeq3(j1,j3)*zaea(kspec,3) + tau_abs_aero(j1,j3,kspec)      &
                + paeq4(j1,j3)*zaea(kspec,4)                                  &
                + paeq5(j1,j3)*zaea(kspec,5)

          zaeods(j1) = paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1)) +          &
                  tau_scat_dust(j1,j3,kspec)*(1.0_dp-(asym_dust(j1,j3,kspec)**2))   &
                + paeq2(j1,j3)*zaes(kspec,2)*(1.0_dp-zaef(kspec,2))                 &
                + paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3)) +               &
                  tau_scat_aero(j1,j3,kspec)*(1.0_dp-(asym_aero(j1,j3,kspec)**2))   &
                + paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))                 &
                + paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))

          zzg= (paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1))*zaeg(kspec,1) +                               &
                  tau_scat_dust(j1,j3,kspec)*(1.0_dp-(asym_dust(j1,j3,kspec)**2))*asym_dust(j1,j3,kspec)        &
                + paeq2(j1,j3)*zaes(kspec,2)*(1.0_dp-zaef(kspec,2))*zaeg(kspec,2)                               &
                + paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3))*zaeg(kspec,3)  +                            &
                  tau_scat_aero(j1,j3,kspec)*(1.0_dp-(asym_aero(j1,j3,kspec)**2.0_dp))*asym_aero(j1,j3,kspec)   &
                + paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))*zaeg(kspec,4)                               &
                + paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))*zaeg(kspec,5) )                             &
                / MAX( zaeods(j1),zepopd)
          zaeb0 (j1) = z1d8*(4.0_dp+zzg)/(1.0_dp+zzg)
          zaeb  (j1) = 0.5_dp-0.75_dp*psmu0(j1)*zzg/(1.0_dp+zzg)
        ENDDO
        !$acc end parallel

      ELSEIF ((.NOT. lrad_dust) .AND. (.NOT. lrad_seas) .AND. (lrad_aero)) THEN

        !$acc parallel
        !$acc loop gang vector
        DO j1 = ki1sc, ki1ec
          zaeoda(j1) = paeq1(j1,j3)*zaea(kspec,1)                               &
                     + paeq2(j1,j3)*zaea(kspec,2)                               &
                     + paeq3(j1,j3)*zaea(kspec,3) + tau_abs_aero(j1,j3,kspec)   &
                     + paeq4(j1,j3)*zaea(kspec,4)                               &
                     + paeq5(j1,j3)*zaea(kspec,5)

          zaeods(j1) = paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1))                     &
                      +paeq2(j1,j3)*zaes(kspec,2)*(1.0_dp-zaef(kspec,2))                     &
                      +paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3))+                    &
                       tau_scat_aero(j1,j3,kspec)*(1.0_dp-(asym_aero(j1,j3,kspec)**2.0_dp))  &
                      +paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))                     &
                      +paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))

          zzg=( paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1))*zaeg(kspec,1)                              &
               +paeq2(j1,j3)*zaes(kspec,2)*(1.0_dp-zaef(kspec,2))*zaeg(kspec,2)                              &
               +paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3))*zaeg(kspec,3)  +                           &
                tau_scat_aero(j1,j3,kspec)*(1.0_dp-(asym_aero(j1,j3,kspec)**2.0_dp))*asym_aero(j1,j3,kspec)  &
               +paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))*zaeg(kspec,4)                              &
               +paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))*zaeg(kspec,5) )                            &
                   / MAX( zaeods(j1),zepopd)
          zaeb0 (j1) = z1d8*(4.0_dp+zzg)/(1.0_dp+zzg)
          zaeb  (j1) = 0.5_dp-0.75_dp*psmu0(j1)*zzg/(1.0_dp+zzg)
        ENDDO
        !$acc end parallel
   
      ELSEIF ((lrad_dust) .AND. (lrad_seas) .AND. (.NOT. lrad_aero)) THEN   

        !$acc parallel
        !$acc loop gang vector
        DO j1 = ki1sc, ki1ec
          zaeoda(j1) = paeq1(j1,j3)*zaea(kspec,1) + tau_abs_dust(j1,j3,kspec)   &
                     + tau_abs_seas(j1,j3,kspec)                                &
                     + paeq3(j1,j3)*zaea(kspec,3)                               &
                     + paeq4(j1,j3)*zaea(kspec,4)                               &
                     + paeq5(j1,j3)*zaea(kspec,5)

          zaeods(j1) = paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1))+             &
                  tau_scat_dust(j1,j3,kspec)*(1.0_dp-(asym_dust(j1,j3,kspec)**2))     &
                + tau_scat_seas(j1,j3,kspec)*(1.0_dp-(asym_seas(j1,j3,kspec)**2))     &
                + paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3))                   &
                + paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))                   &
                + paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))

          zzg=( paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1))*zaeg(kspec,1) +                       &
                tau_scat_dust(j1,j3,kspec)*(1.0_dp-(asym_dust(j1,j3,kspec)**2))*asym_dust(j1,j3,kspec)  &
               +tau_scat_seas(j1,j3,kspec)*(1.0_dp-(asym_seas(j1,j3,kspec)**2))*asym_seas(j1,j3,kspec)  &
               +paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3))*zaeg(kspec,3)                         &
               +paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))*zaeg(kspec,4)                         &
               +paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))*zaeg(kspec,5) )                       &
             / MAX( zaeods(j1),zepopd)
          zaeb0 (j1) = z1d8*(4.0_dp+zzg)/(1.0_dp+zzg)
          zaeb  (j1) = 0.5_dp-0.75_dp*psmu0(j1)*zzg/(1.0_dp+zzg)
        ENDDO
        !$acc end parallel

      ELSEIF ((.NOT.lrad_dust) .AND. (lrad_seas) .AND. (.NOT. lrad_aero)) THEN

        !$acc parallel
        !$acc loop gang vector
        DO j1 = ki1sc, ki1ec
          zaeoda(j1) = paeq1(j1,j3)*zaea(kspec,1)   &
                     + tau_abs_seas(j1,j3,kspec)    &
                     + paeq3(j1,j3)*zaea(kspec,3)   &
                     + paeq4(j1,j3)*zaea(kspec,4)   &
                     + paeq5(j1,j3)*zaea(kspec,5)

          zaeods(j1) =  paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1))            &
                  + tau_scat_seas(j1,j3,kspec)*(1.0_dp-(asym_seas(j1,j3,kspec)**2))  &
                  + paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3))                &
                  + paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))                &
                  + paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))

          zzg=   ( paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1))*zaeg(kspec,1)                          &
                  +tau_scat_seas(j1,j3,kspec)*(1.0_dp-(asym_seas(j1,j3,kspec)**2))*asym_seas(j1,j3,kspec)   &
                  +paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3))*zaeg(kspec,3)                          &
                  +paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))*zaeg(kspec,4)                          &
                  +paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))*zaeg(kspec,5) )                        &
                  / MAX( zaeods(j1),zepopd)
          zaeb0 (j1) = z1d8*(4.0_dp+zzg)/(1.0_dp+zzg)
          zaeb  (j1) = 0.5_dp-0.75_dp*psmu0(j1)*zzg/(1.0_dp+zzg)
        ENDDO
        !$acc end parallel

      ELSEIF ((lrad_dust) .AND. (lrad_seas) .AND. (lrad_aero)) THEN

        !$acc parallel
        !$acc loop gang vector
        DO j1 = ki1sc, ki1ec
          zaeoda(j1) = paeq1(j1,j3)*zaea(kspec,1) + tau_abs_dust(j1,j3,kspec)    &
                + tau_abs_seas(j1,j3,kspec)                                      &
                + paeq3(j1,j3)*zaea(kspec,3) + tau_abs_aero(j1,j3,kspec)         &
                + paeq4(j1,j3)*zaea(kspec,4)                                     &
                + paeq5(j1,j3)*zaea(kspec,5)

          zaeods(j1) = paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1)) +           &
                  tau_scat_dust(j1,j3,kspec)*(1.0_dp-(asym_dust(j1,j3,kspec)**2))    &
                 +tau_scat_seas(j1,j3,kspec)*(1.0_dp-(asym_seas(j1,j3,kspec)**2))    &
                 +paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3)) +                &
                  tau_scat_aero(j1,j3,kspec)*(1.0_dp-(asym_aero(j1,j3,kspec)**2))    &
                 +paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))                  &
                 +paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))

          zzg= ( paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1))*zaeg(kspec,1) +                             &
                  tau_scat_dust(j1,j3,kspec)*(1.0_dp-(asym_dust(j1,j3,kspec)**2))*asym_dust(j1,j3,kspec)       &
                 +tau_scat_seas(j1,j3,kspec)*(1.0_dp-(asym_seas(j1,j3,kspec)**2))*asym_seas(j1,j3,kspec)       &
                 +paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3))*zaeg(kspec,3)  +                           &
                  tau_scat_aero(j1,j3,kspec)*(1.0_dp-(asym_aero(j1,j3,kspec)**2.0_dp))*asym_aero(j1,j3,kspec)  &
                 +paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))*zaeg(kspec,4)                              &
                 +paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))*zaeg(kspec,5) )                            &
                / MAX( zaeods(j1),zepopd)
          zaeb0 (j1) = z1d8*(4.0_dp+zzg)/(1.0_dp+zzg)
          zaeb  (j1) = 0.5_dp-0.75_dp*psmu0(j1)*zzg/(1.0_dp+zzg)
        ENDDO
        !$acc end parallel

      ELSEIF ((.NOT. lrad_dust) .AND. (lrad_seas) .AND. (lrad_aero)) THEN

        !$acc parallel
        !$acc loop gang vector
        DO j1 = ki1sc, ki1ec
          zaeoda(j1) = paeq1(j1,j3)*zaea(kspec,1)                        &
                + tau_abs_seas(j1,j3,kspec)                              &
                + paeq3(j1,j3)*zaea(kspec,3) + tau_abs_aero(j1,j3,kspec) &
                + paeq4(j1,j3)*zaea(kspec,4)                             &
                + paeq5(j1,j3)*zaea(kspec,5)

          zaeods(j1) = paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1))                 &
                 +tau_scat_seas(j1,j3,kspec)*(1.0_dp-(asym_seas(j1,j3,kspec)**2.0_dp))   &
                 +paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3)) +                    &
                  tau_scat_aero(j1,j3,kspec)*(1.0_dp-(asym_aero(j1,j3,kspec)**2.0_dp))   &
                 +paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))                      &
                 +paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))

          zzg=( paeq1(j1,j3)*zaes(kspec,1)*(1.0_dp-zaef(kspec,1))*zaeg(kspec,1)                              &
               +tau_scat_seas(j1,j3,kspec)*(1.0_dp-(asym_seas(j1,j3,kspec)**2.0_dp))*asym_seas(j1,j3,kspec)  &
               +paeq3(j1,j3)*zaes(kspec,3)*(1.0_dp-zaef(kspec,3))*zaeg(kspec,3) +                            &
                tau_scat_aero(j1,j3,kspec)*(1.0_dp-(asym_aero(j1,j3,kspec)**2.0_dp))*asym_aero(j1,j3,kspec)  &
               +paeq4(j1,j3)*zaes(kspec,4)*(1.0_dp-zaef(kspec,4))*zaeg(kspec,4)                              &
               +paeq5(j1,j3)*zaes(kspec,5)*(1.0_dp-zaef(kspec,5))*zaeg(kspec,5) )                            &
                   / MAX( zaeods(j1),zepopd)
          zaeb0 (j1) = z1d8*(4.0_dp+zzg)/(1.0_dp+zzg)
          zaeb  (j1) = 0.5_dp-0.75_dp*psmu0(j1)*zzg/(1.0_dp+zzg)
        ENDDO
        !$acc end parallel
   
      ENDIF
    ENDIF
#endif

#ifndef _OPENACC
    IF (ldebug .AND. jb==jindex) THEN
      print *,' zaeoda  (j1b,jb) :',zaeoda  (j1b)
      print *,' zaeods  (j1b,jb) :',zaeods  (j1b)
      print *,' zaeb0   (j1b,jb) :',zaeb0   (j1b)
      print *,' zaeb    (j1b,jb) :',zaeb    (j1b)
    ENDIF   
#endif

    ! Optical thickness for Rayleigh scattering

    !$acc parallel
    !$acc loop gang vector
    DO j1 = ki1sc, ki1ec
      zraods(j1) = zrsc(isp) * pdp(j1,j3) / (1.0_dp+ zrsc(isp) * &
                  (papra(j1) + 0.5_dp* pdp(j1,j3) * pqsmu0(j1)))
      papra (j1) = papra(j1) + pdp(j1,j3) * pqsmu0(j1)
    ENDDO
    !$acc end parallel

#ifndef _OPENACC
    IF (ldebug .AND. jb==jindex) THEN
      !  print *,' Rayleigh coefficient:',zrsc(isp)
      !  print *,' Papra  (j1b,jb)    :',papra   (j1b)
      !  print *,' Pqsmu0 (j1b,jb)    :',pqsmu0  (j1b)
      !  print *,' Pdp    (j1b,jb)    :',pdp     (j1b,j3)
         print *,' zraods (j1b,jb)    :',zraods  (j1b)
    ENDIF   
#endif

    !-----------------------------------------------------------------------
 
    ! Linear combination of individual contributions 
 

    ! a) cloud free part of layer

    !$acc parallel
    !$acc loop gang vector
    DO j1 = ki1sc, ki1ec
      podaf(j1,j3) = MAX ( zaeoda(j1)             , zepopd)
      podsf(j1,j3) = MAX ( zaeods(j1) + zraods(j1), zepopd)
      podsf(j1,j3) = MIN ( podsf (j1,j3),                             &
                 (1.0_dp-zepssa) * (podaf(j1,j3) + podsf(j1,j3)) )
      pbsff(j1,j3) = (zaeb0(j1) * zaeods(j1)                          &
                       + 0.5_dp * zraods(j1)) / podsf(j1,j3)
      pusff(j1,j3) = (zaeb (j1) * zaeods(j1)                          &
                       + 0.5_dp * zraods(j1)) / podsf(j1,j3)
    ENDDO
    !$acc end parallel
 
!-------------------------------------------------------------------------------
    ! b) cloudy part of layer

    !$acc parallel
    !$acc loop gang vector
    DO j1 = ki1sc, ki1ec
      podac(j1,j3) = zlwoda(j1) + ziwoda(j1) + zaeoda(j1)
      podsc(j1,j3) = zlwods(j1) + ziwods(j1) + zaeods(j1) + zraods(j1)
      podac(j1,j3) = MAX( podac(j1,j3), zepopd)
      podsc(j1,j3) = MAX( podsc(j1,j3), zepopd)
      podsc(j1,j3) = MIN( podsc(j1,j3),                                &
                (1.0_dp-zepssa) * (podac(j1,j3) + podsc(j1,j3)))

      pbsfc(j1,j3)= (zlwb0(j1) * zlwods(j1)                          &
                   + ziwb0(j1) * ziwods(j1) + zaeb0(j1) * zaeods(j1) &
                   +    0.5_dp * zraods(j1)) / podsc(j1,j3)
      pusfc(j1,j3)= (zlwb (j1) * zlwods(j1)                          &
                   + ziwb (j1) * ziwods(j1) + zaeb (j1) * zaeods(j1) &
                   +    0.5_dp * zraods(j1)) / podsc(j1,j3)
    ENDDO
    !$acc end parallel
 
  ! End of vertical loop
  END DO  
 
#ifndef _OPENACC
  IF (ldebug .AND. jb==jindex) THEN
      print *,' podaf   (j1b,jb,1) :',podaf   (j1b,1)
      print *,' podsf   (j1b,jb,1) :',podsf   (j1b,1)
      print *,' pbsff   (j1b,jb,1) :',pbsff   (j1b,1)
      print *,' pusff   (j1b,jb,1) :',pusff   (j1b,1)
      print *,' podac   (j1b,jb,1) :',podac   (j1b,1)
      print *,' podsc   (j1b,jb,1) :',podsc   (j1b,1)
      print *,' pbsfc   (j1b,jb,1) :',pbsfc   (j1b,1)
      print *,' pusfc   (j1b,jb,1) :',pusfc   (j1b,1)
      print *,' -------------------------------------------------'
  ENDIF     
#endif
 
!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

  !$acc end data

END SUBROUTINE opt_so

!==============================================================================
!==============================================================================
!+ Module procedure in "Radiation" 
!------------------------------------------------------------------------------

SUBROUTINE radiation_rg_wkarr_alloc(ki1ed,ki3ed,istat)

  INTEGER, INTENT(IN)  :: ki1ed,ki3ed
  INTEGER, INTENT(OUT) :: istat
  INTEGER              :: ki1sd, ki3sd

  istat = 0
  ki1sd = 1
  ki3sd = 1

  !fesft_sp
  IF (wp == sp) THEN
  ALLOCATE( pti_dp       (ki1sd:ki1ed,ki3sd:ki3ed+1), STAT=istat)
  ALLOCATE( pdp_dp       (ki1sd:ki1ed,ki3sd:ki3ed)  , STAT=istat)
  ALLOCATE( pclc_dp      (ki1sd:ki1ed,ki3sd:ki3ed)  , STAT=istat)
  ALLOCATE( pwv_dp       (ki1sd:ki1ed,ki3sd:ki3ed)  , STAT=istat)
  ALLOCATE( psw_dp       (ki1sd:ki1ed,ki3sd:ki3ed)  , STAT=istat)
  ALLOCATE( pqlwc_dp     (ki1sd:ki1ed,ki3sd:ki3ed)  , STAT=istat)
  ALLOCATE( pqiwc_dp     (ki1sd:ki1ed,ki3sd:ki3ed)  , STAT=istat)
  ALLOCATE( pduco2_dp    (ki1sd:ki1ed,ki3sd:ki3ed)  , STAT=istat)
  ALLOCATE( pduo3_dp     (ki1sd:ki1ed,ki3sd:ki3ed)  , STAT=istat)
  ALLOCATE( paeq1_dp     (ki1sd:ki1ed,ki3sd:ki3ed)  , STAT=istat)
  ALLOCATE( paeq2_dp     (ki1sd:ki1ed,ki3sd:ki3ed)  , STAT=istat)
  ALLOCATE( paeq3_dp     (ki1sd:ki1ed,ki3sd:ki3ed)  , STAT=istat)
  ALLOCATE( paeq4_dp     (ki1sd:ki1ed,ki3sd:ki3ed)  , STAT=istat)
  ALLOCATE( paeq5_dp     (ki1sd:ki1ed,ki3sd:ki3ed)  , STAT=istat)
  ALLOCATE( psmu0_dp     (ki1sd:ki1ed)              , STAT=istat)
  ALLOCATE( palth_dp     (ki1sd:ki1ed)              , STAT=istat)
  ALLOCATE( palso_dp     (ki1sd:ki1ed)              , STAT=istat)
  ALLOCATE( pskyview_dp  (ki1sd:ki1ed)              , STAT=istat)
  ALLOCATE( pfcor_dp     (ki1sd:ki1ed)              , STAT=istat)
  ALLOCATE( papre_dp     (ki1sd:ki1ed)              , STAT=istat)
  ALLOCATE( pflt_dp      (ki1sd:ki1ed,ki3sd:ki3ed+1), STAT=istat)
  ALLOCATE( pfls_dp      (ki1sd:ki1ed,ki3sd:ki3ed+1), STAT=istat)
  ALLOCATE( pflt_s_dp    (ki1sd:ki1ed)              , STAT=istat)
  ALLOCATE( pfls_s_dp    (ki1sd:ki1ed)              , STAT=istat)
  ALLOCATE( pflsdir_dp   (ki1sd:ki1ed,ki3sd:ki3ed+1), STAT=istat)
  ALLOCATE( pfltd_dp     (ki1sd:ki1ed)              , STAT=istat)
  ALLOCATE( pfltu_dp     (ki1sd:ki1ed)              , STAT=istat)
  ALLOCATE( pflsd_dp     (ki1sd:ki1ed)              , STAT=istat)
  ALLOCATE( pflsu_dp     (ki1sd:ki1ed)              , STAT=istat)
  ALLOCATE( pflsp_dp     (ki1sd:ki1ed)              , STAT=istat)
  ALLOCATE( pflpar_dp    (ki1sd:ki1ed)              , STAT=istat)
  ALLOCATE( pflsu_par_dp (ki1sd:ki1ed)              , STAT=istat)
  ALLOCATE( pflsd_par_dp (ki1sd:ki1ed)              , STAT=istat)
  ALLOCATE( pflsp_par_dp (ki1sd:ki1ed)              , STAT=istat)
  !$acc enter data &
  !$acc create ( pti_dp,pdp_dp,pclc_dp,pwv_dp,psw_dp,pqlwc_dp   ) &
  !$acc create ( pqiwc_dp,pduco2_dp,pduo3_dp,paeq1_dp,paeq2_dp  ) &
  !$acc create ( paeq3_dp,paeq4_dp,paeq5_dp,psmu0_dp,palth_dp   ) &
  !$acc create ( palso_dp,pskyview_dp,pfcor_dp,papre_dp,pflt_dp ) &
  !$acc create ( pfls_dp,pflt_s_dp,pfls_s_dp,pflsdir_dp         ) &
  !$acc create ( pfltd_dp,pfltu_dp,pflsd_dp,pflsu_dp,pflsp_dp   ) &
  !$acc create ( pflpar_dp,pflsu_par_dp,pflsd_par_dp            ) &
  !$acc create ( pflsp_par_dp                                   )
  ENDIF

  !fesft_dp
  ALLOCATE( zketyp   (jpther)                    , STAT=istat )
  ALLOCATE( ztetyp   (jpther)                    , STAT=istat )
  ALLOCATE( zflux    (ki1sd:ki1ed,ki3sd:ki3ed+1) , STAT=istat )
  ALLOCATE( zflux_c  (ki1sd:ki1ed,ki3sd:ki3ed+1) , STAT=istat )
  ALLOCATE( zfluxi   (ki1sd:ki1ed,ki3sd:ki3ed+1) , STAT=istat )
  ALLOCATE( zfluxu   (ki1sd:ki1ed,ki3sd:ki3ed+1) , STAT=istat )
  ALLOCATE( zfluxu_c (ki1sd:ki1ed,ki3sd:ki3ed+1) , STAT=istat )
  ALLOCATE( zfluxui  (ki1sd:ki1ed,ki3sd:ki3ed+1) , STAT=istat )
  ALLOCATE( zfluxd   (ki1sd:ki1ed,ki3sd:ki3ed+1) , STAT=istat )
  ALLOCATE( zfluxd_c (ki1sd:ki1ed,ki3sd:ki3ed+1) , STAT=istat )
  ALLOCATE( zfluxdi  (ki1sd:ki1ed,ki3sd:ki3ed+1) , STAT=istat )
  ALLOCATE( zfgas    (ki1sd:ki1ed,ki3sd:ki3ed+1) , STAT=istat )
  ALLOCATE( zfgasu   (ki1sd:ki1ed,ki3sd:ki3ed+1) , STAT=istat )
  ALLOCATE( zfgasd   (ki1sd:ki1ed,ki3sd:ki3ed+1) , STAT=istat )
  ALLOCATE( pbbr     (ki1sd:ki1ed,ki3sd:ki3ed+1) , STAT=istat )
  ALLOCATE( pflpt    (ki1sd:ki1ed)               , STAT=istat )
  ALLOCATE( palp     (ki1sd:ki1ed)               , STAT=istat )
  ALLOCATE( pqsmu0   (ki1sd:ki1ed)               , STAT=istat )
  ALLOCATE( palogt   (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( palogp   (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( papra    (ki1sd:ki1ed)               , STAT=istat )
  ALLOCATE( pduh2oc  (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( pduh2of  (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( pdulwc   (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( pduiwc   (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( prholwc  (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( prhoiwc  (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( zduetpc  (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( zduetpf  (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( zapre    (ki1sd:ki1ed)               , STAT=istat )
  ALLOCATE( ztm      (ki1sd:ki1ed)               , STAT=istat )
  ALLOCATE( zzwv     (ki1sd:ki1ed)               , STAT=istat )
  ALLOCATE( zcpo     (ki1sd:ki1ed)               , STAT=istat )
  ALLOCATE( zcpn     (ki1sd:ki1ed)               , STAT=istat )
  ALLOCATE( zcmo     (ki1sd:ki1ed)               , STAT=istat )
  ALLOCATE( zcmn     (ki1sd:ki1ed)               , STAT=istat )
  ALLOCATE( podac    (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( podaf    (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( podsc    (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( podsf    (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( pbsfc    (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( pbsff    (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( pusfc    (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( pusff    (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( pca1     (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( pcb1     (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( pcc1     (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( pcd1     (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( pca2     (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( pcb2     (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( pcc2     (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( pcd2     (ki1sd:ki1ed,ki3sd:ki3ed)   , STAT=istat )
  ALLOCATE( pflfd    (ki1sd:ki1ed,ki3sd:ki3ed+1) , STAT=istat )
  ALLOCATE( pflfu    (ki1sd:ki1ed,ki3sd:ki3ed+1) , STAT=istat )
  ALLOCATE( pflfp    (ki1sd:ki1ed,ki3sd:ki3ed+1) , STAT=istat )
  ALLOCATE( pflcd    (ki1sd:ki1ed,ki3sd:ki3ed+1) , STAT=istat )
  ALLOCATE( pflcu    (ki1sd:ki1ed,ki3sd:ki3ed+1) , STAT=istat )
  ALLOCATE( pflcp    (ki1sd:ki1ed,ki3sd:ki3ed+1) , STAT=istat )
  !$acc enter data &
  !$acc create  ( zketyp,ztetyp,zflux,zflux_c,zfluxi,zfluxu     ) &
  !$acc create  ( zfluxu_c,zfluxui,zfluxd,zfluxd_c,zfluxdi      ) &
  !$acc create  ( zfgas,zfgasu,zfgasd,pbbr,pflpt,palp,pqsmu0    ) &
  !$acc create  ( palogt,palogp,papra,pduh2oc,pduh2of,pdulwc    ) &
  !$acc create  ( pduiwc,prholwc,prhoiwc,zduetpc,zduetpf,zapre  ) &
  !$acc create  ( ztm,zzwv,zcpo,zcpn,zcmo,zcmn,podac,podaf      ) &
  !$acc create  ( podsc,podsf,pbsfc,pbsff,pusfc,pusff,pca1,pcb1 ) &
  !$acc create  ( pcc1,pcd1,pca2,pcb2,pcc2,pcd2,pflfd,pflfu     ) &
  !$acc create  ( pflfp,pflcd,pflcu,pflcp                       )

  !inv_th/so
  ALLOCATE( pa1c (ki1sd:ki1ed)            , STAT=istat)
  ALLOCATE( pa1f (ki1sd:ki1ed)            , STAT=istat)
  ALLOCATE( pa2c (ki1sd:ki1ed)            , STAT=istat)
  ALLOCATE( pa2f (ki1sd:ki1ed)            , STAT=istat)
  ALLOCATE( pa3c (ki1sd:ki1ed)            , STAT=istat)
  ALLOCATE( pa3f (ki1sd:ki1ed)            , STAT=istat)
  ALLOCATE( pa4c (ki1sd:ki1ed)            , STAT=istat)
  ALLOCATE( pa4f (ki1sd:ki1ed)            , STAT=istat)
  ALLOCATE( pa5c (ki1sd:ki1ed)            , STAT=istat)
  ALLOCATE( pa5f (ki1sd:ki1ed)            , STAT=istat)
  ALLOCATE( ztu1 (ki1sd:ki1ed,ki3sd:ki3ed), STAT=istat)
  ALLOCATE( ztu2 (ki1sd:ki1ed,ki3sd:ki3ed), STAT=istat)
  ALLOCATE( ztu3 (ki1sd:ki1ed,ki3sd:ki3ed), STAT=istat)
  ALLOCATE( ztu4 (ki1sd:ki1ed,ki3sd:ki3ed), STAT=istat)
  ALLOCATE( ztu5 (ki1sd:ki1ed,ki3sd:ki3ed), STAT=istat)
  ALLOCATE( ztu6 (ki1sd:ki1ed,ki3sd:ki3ed), STAT=istat)
  ALLOCATE( ztu7 (ki1sd:ki1ed,ki3sd:ki3ed), STAT=istat)
  ALLOCATE( ztu8 (ki1sd:ki1ed,ki3sd:ki3ed), STAT=istat)
  ALLOCATE( ztu9 (ki1sd:ki1ed,ki3sd:ki3ed), STAT=istat)
  !$acc enter data &
  !$acc create  ( pa1c,pa1f,pa2c,pa2f,pa3c,pa4f,pa4c,pa5f,pa5c,pa3f   ) &
  !$acc create  ( ztu1,ztu2,ztu3,ztu4,ztu5,ztu6,ztu7,ztu8,ztu9        )
 
  !opt_th/so
  ALLOCATE( zlwoda (ki1sd:ki1ed), STAT=istat )
  ALLOCATE( zlwods (ki1sd:ki1ed), STAT=istat )
  ALLOCATE( zlwb0  (ki1sd:ki1ed), STAT=istat )
  ALLOCATE( zlwb   (ki1sd:ki1ed), STAT=istat )
  ALLOCATE( ziwoda (ki1sd:ki1ed), STAT=istat )
  ALLOCATE( ziwods (ki1sd:ki1ed), STAT=istat )
  ALLOCATE( ziwb0  (ki1sd:ki1ed), STAT=istat )
  ALLOCATE( ziwb   (ki1sd:ki1ed), STAT=istat )
  ALLOCATE( zaeoda (ki1sd:ki1ed), STAT=istat )
  ALLOCATE( zaeods (ki1sd:ki1ed), STAT=istat )
  ALLOCATE( zaeb0  (ki1sd:ki1ed), STAT=istat )
  ALLOCATE( zaeb   (ki1sd:ki1ed), STAT=istat )
  ALLOCATE( zraods (ki1sd:ki1ed), STAT=istat )
  !$acc enter data &
  !$acc create ( zlwoda,zlwods,zlwb0,zlwb,ziwoda,ziwods,ziwb0,ziwb     ) &
  !$acc create ( zaeoda,zaeods,zaeb0,zaeb,zraods                       )

#ifdef COSMOART
    ALLOCATE( asym_dust     (ki1sd:ki1ed,ki3sd:ki3ed,8), STAT=istat)
    ALLOCATE( tau_abs_dust  (ki1sd:ki1ed,ki3sd:ki3ed,8), STAT=istat)
    ALLOCATE( tau_scat_dust (ki1sd:ki1ed,ki3sd:ki3ed,8), STAT=istat)
    asym_dust      = 0.0_wp
    tau_abs_dust   = 0.0_wp
    tau_scat_dust  = 0.0_wp
    ALLOCATE( asym_seas     (ki1sd:ki1ed,ki3sd:ki3ed,8), STAT=istat)
    ALLOCATE( tau_abs_seas  (ki1sd:ki1ed,ki3sd:ki3ed,8), STAT=istat)
    ALLOCATE( tau_scat_seas (ki1sd:ki1ed,ki3sd:ki3ed,8), STAT=istat)
    asym_seas      = 0.0_wp
    tau_abs_seas   = 0.0_wp
    tau_scat_seas  = 0.0_wp
    ALLOCATE( asym_aero     (ki1sd:ki1ed,ki3sd:ki3ed,8), STAT=istat)
    ALLOCATE( tau_abs_aero  (ki1sd:ki1ed,ki3sd:ki3ed,8), STAT=istat)
    ALLOCATE( tau_scat_aero (ki1sd:ki1ed,ki3sd:ki3ed,8), STAT=istat)
    asym_aero      = 0.0_wp
    tau_abs_aero   = 0.0_wp
    tau_scat_aero  = 0.0_wp
#endif

END SUBROUTINE radiation_rg_wkarr_alloc

!==============================================================================
!==============================================================================
!+ Module procedure in "Radiation" 
!------------------------------------------------------------------------------

SUBROUTINE radiation_rg_wkarr_dealloc(istat)

  INTEGER, INTENT(OUT) :: istat

  istat = 0

  !fesft_sp
  IF (wp == sp) THEN
  !$acc exit data &
  !$acc delete ( pti_dp,pdp_dp,pclc_dp,pwv_dp,psw_dp,pqlwc_dp   ) &
  !$acc delete ( pqiwc_dp,pduco2_dp,pduo3_dp,paeq1_dp,paeq2_dp  ) &
  !$acc delete ( paeq3_dp,paeq4_dp,paeq5_dp,psmu0_dp,palth_dp   ) &
  !$acc delete ( palso_dp,pskyview_dp,pfcor_dp,papre_dp,pflt_dp ) &
  !$acc delete ( pfls_dp,pflt_s_dp,pfls_s_dp,pflsdir_dp         ) &
  !$acc delete ( pfltd_dp,pfltu_dp,pflsd_dp,pflsu_dp,pflsp_dp   ) &
  !$acc delete ( pflpar_dp,pflsu_par_dp,pflsd_par_dp            ) &
  !$acc delete ( pflsp_par_dp                                   )
  DEALLOCATE( pti_dp       , STAT=istat)
  DEALLOCATE( pdp_dp       , STAT=istat)
  DEALLOCATE( pclc_dp      , STAT=istat)
  DEALLOCATE( pwv_dp       , STAT=istat)
  DEALLOCATE( psw_dp       , STAT=istat)
  DEALLOCATE( pqlwc_dp     , STAT=istat)
  DEALLOCATE( pqiwc_dp     , STAT=istat)
  DEALLOCATE( pduco2_dp    , STAT=istat)
  DEALLOCATE( pduo3_dp     , STAT=istat)
  DEALLOCATE( paeq1_dp     , STAT=istat)
  DEALLOCATE( paeq2_dp     , STAT=istat)
  DEALLOCATE( paeq3_dp     , STAT=istat)
  DEALLOCATE( paeq4_dp     , STAT=istat)
  DEALLOCATE( paeq5_dp     , STAT=istat)
  DEALLOCATE( psmu0_dp     , STAT=istat)
  DEALLOCATE( palth_dp     , STAT=istat)
  DEALLOCATE( palso_dp     , STAT=istat)
  DEALLOCATE( pskyview_dp  , STAT=istat)
  DEALLOCATE( pfcor_dp     , STAT=istat)
  DEALLOCATE( papre_dp     , STAT=istat)
  DEALLOCATE( pflt_dp      , STAT=istat)
  DEALLOCATE( pfls_dp      , STAT=istat)
  DEALLOCATE( pflt_s_dp    , STAT=istat)
  DEALLOCATE( pfls_s_dp    , STAT=istat)
  DEALLOCATE( pflsdir_dp   , STAT=istat)
  DEALLOCATE( pfltd_dp     , STAT=istat)
  DEALLOCATE( pfltu_dp     , STAT=istat)
  DEALLOCATE( pflsd_dp     , STAT=istat)
  DEALLOCATE( pflsu_dp     , STAT=istat)
  DEALLOCATE( pflsp_dp     , STAT=istat)
  DEALLOCATE( pflpar_dp    , STAT=istat)
  DEALLOCATE( pflsu_par_dp , STAT=istat)
  DEALLOCATE( pflsd_par_dp , STAT=istat)
  DEALLOCATE( pflsp_par_dp , STAT=istat)
  ENDIF
    
  !fesft
  !$acc exit data &
  !$acc delete ( zketyp,ztetyp,zflux,zflux_c,zfluxi,zfluxu     ) &
  !$acc delete ( zfluxu_c,zfluxui,zfluxd,zfluxd_c,zfluxdi      ) &
  !$acc delete ( zfgas,zfgasu,zfgasd,pbbr,pflpt,palp,pqsmu0    ) &
  !$acc delete ( palogt,palogp,papra,pduh2oc,pduh2of,pdulwc    ) &
  !$acc delete ( pduiwc,prholwc,prhoiwc,zduetpc,zduetpf,zapre  ) &
  !$acc delete ( ztm,zzwv,zcpo,zcpn,zcmo,zcmn,podac,podaf,podsc) &
  !$acc delete ( podsf,pbsfc,pbsff,pusfc,pusff,pca1,pcb1,pcc1  ) &
  !$acc delete ( pcd1,pca2,pcb2,pcc2,pcd2,pflfd,pflfu,pflfp    ) &
  !$acc delete ( pflcd,pflcu,pflcp                             )
  DEALLOCATE( zketyp   , STAT=istat )
  DEALLOCATE( ztetyp   , STAT=istat )
  DEALLOCATE( zflux    , STAT=istat )
  DEALLOCATE( zflux_c  , STAT=istat )
  DEALLOCATE( zfluxi   , STAT=istat )
  DEALLOCATE( zfluxu   , STAT=istat )
  DEALLOCATE( zfluxu_c , STAT=istat )
  DEALLOCATE( zfluxui  , STAT=istat )
  DEALLOCATE( zfluxd   , STAT=istat )
  DEALLOCATE( zfluxd_c , STAT=istat )
  DEALLOCATE( zfluxdi  , STAT=istat )
  DEALLOCATE( zfgas    , STAT=istat )
  DEALLOCATE( zfgasu   , STAT=istat )
  DEALLOCATE( zfgasd   , STAT=istat )
  DEALLOCATE( pbbr     , STAT=istat )
  DEALLOCATE( pflpt    , STAT=istat )
  DEALLOCATE( palp     , STAT=istat )
  DEALLOCATE( pqsmu0   , STAT=istat )
  DEALLOCATE( palogt   , STAT=istat )
  DEALLOCATE( palogp   , STAT=istat )
  DEALLOCATE( papra    , STAT=istat )
  DEALLOCATE( pduh2oc  , STAT=istat )
  DEALLOCATE( pduh2of  , STAT=istat )
  DEALLOCATE( pdulwc   , STAT=istat )
  DEALLOCATE( pduiwc   , STAT=istat )
  DEALLOCATE( prholwc  , STAT=istat )
  DEALLOCATE( prhoiwc  , STAT=istat )
  DEALLOCATE( zduetpc  , STAT=istat )
  DEALLOCATE( zduetpf  , STAT=istat )
  DEALLOCATE( zapre    , STAT=istat )
  DEALLOCATE( ztm      , STAT=istat )
  DEALLOCATE( zzwv     , STAT=istat )
  DEALLOCATE( zcpo     , STAT=istat )
  DEALLOCATE( zcpn     , STAT=istat )
  DEALLOCATE( zcmo     , STAT=istat )
  DEALLOCATE( zcmn     , STAT=istat )
  DEALLOCATE( podac    , STAT=istat )
  DEALLOCATE( podaf    , STAT=istat )
  DEALLOCATE( podsc    , STAT=istat )
  DEALLOCATE( podsf    , STAT=istat )
  DEALLOCATE( pbsfc    , STAT=istat )
  DEALLOCATE( pbsff    , STAT=istat )
  DEALLOCATE( pusfc    , STAT=istat )
  DEALLOCATE( pusff    , STAT=istat )
  DEALLOCATE( pca1     , STAT=istat )
  DEALLOCATE( pcb1     , STAT=istat )
  DEALLOCATE( pcc1     , STAT=istat )
  DEALLOCATE( pcd1     , STAT=istat )
  DEALLOCATE( pca2     , STAT=istat )
  DEALLOCATE( pcb2     , STAT=istat )
  DEALLOCATE( pcc2     , STAT=istat )
  DEALLOCATE( pcd2     , STAT=istat )
  DEALLOCATE( pflfd    , STAT=istat )
  DEALLOCATE( pflfu    , STAT=istat )
  DEALLOCATE( pflfp    , STAT=istat )
  DEALLOCATE( pflcd    , STAT=istat )
  DEALLOCATE( pflcu    , STAT=istat )
  DEALLOCATE( pflcp    , STAT=istat )

  !inv_th/so
  !$acc exit data &
  !$acc delete ( pa1c,pa1f,pa2c,pa2f,pa3c,pa4f,pa4c,pa5f,pa5c,pa3f   ) &
  !$acc delete ( ztu1,ztu2,ztu3,ztu4,ztu5,ztu6,ztu7,ztu8,ztu9        )
  DEALLOCATE( pa1c, STAT=istat)
  DEALLOCATE( pa1f, STAT=istat)
  DEALLOCATE( pa2c, STAT=istat)
  DEALLOCATE( pa2f, STAT=istat)
  DEALLOCATE( pa3c, STAT=istat)
  DEALLOCATE( pa3f, STAT=istat)
  DEALLOCATE( pa4c, STAT=istat)
  DEALLOCATE( pa4f, STAT=istat)
  DEALLOCATE( pa5c, STAT=istat)
  DEALLOCATE( pa5f, STAT=istat)
  DEALLOCATE( ztu1, STAT=istat)
  DEALLOCATE( ztu2, STAT=istat)
  DEALLOCATE( ztu3, STAT=istat)
  DEALLOCATE( ztu4, STAT=istat)
  DEALLOCATE( ztu5, STAT=istat)
  DEALLOCATE( ztu6, STAT=istat)
  DEALLOCATE( ztu7, STAT=istat)
  DEALLOCATE( ztu8, STAT=istat)
  DEALLOCATE( ztu9, STAT=istat)

  !opt_th/so
  !$acc exit data &
  !$acc delete ( zlwoda,zlwods,zlwb0,zlwb,ziwoda,ziwods,ziwb0,ziwb     ) &
  !$acc delete ( zaeoda,zaeods,zaeb0,zaeb,zraods                       )
  DEALLOCATE( zlwoda, STAT=istat )
  DEALLOCATE( zlwods, STAT=istat )
  DEALLOCATE( zlwb0 , STAT=istat )
  DEALLOCATE( zlwb  , STAT=istat )
  DEALLOCATE( ziwoda, STAT=istat )
  DEALLOCATE( ziwods, STAT=istat )
  DEALLOCATE( ziwb0 , STAT=istat )
  DEALLOCATE( ziwb  , STAT=istat )
  DEALLOCATE( zaeoda, STAT=istat )
  DEALLOCATE( zaeods, STAT=istat )
  DEALLOCATE( zaeb0 , STAT=istat )
  DEALLOCATE( zaeb  , STAT=istat )
  DEALLOCATE( zraods, STAT=istat )

#ifdef COSMOART
  DEALLOCATE( asym_dust,     STAT=istat)
  DEALLOCATE( tau_abs_dust,  STAT=istat)
  DEALLOCATE( tau_scat_dust, STAT=istat)
  DEALLOCATE( asym_seas,     STAT=istat)
  DEALLOCATE( tau_abs_seas,  STAT=istat)
  DEALLOCATE( tau_scat_seas, STAT=istat)
  DEALLOCATE( asym_aero,     STAT=istat)
  DEALLOCATE( tau_abs_aero,  STAT=istat)
  DEALLOCATE( tau_scat_aero, STAT=istat)
#endif

END SUBROUTINE radiation_rg_wkarr_dealloc

!==============================================================================

END MODULE radiation_rg

!==============================================================================
