!+ Utility module for post-processing utility routines
!------------------------------------------------------------------------------
 
MODULE pp_utilities
 
!------------------------------------------------------------------------------
!
! Description:
!   This module provides service utilities for post-processing.
!     - no routine uses other modules, except the declarations for the
!       KIND-type parameter; the data access is by parameter list only
!     - no routine allocates dynamic memory except for automatic arrays
!     - no derived data types are used
!
!   Routines (module procedures) currently contained:
!
!     - cal_conv_ind
!       Performs a parcel ascent and computes CAPE, CIN and Showalter Index and
!       SWISS indexes, the Lifting Condensation Level and the Level of Free
!       Convection
!
!     - calc_bulk_richardson
!       Compute the bulk Richardson Number at each grid cell
!
!     - calc_pbl_brn
!       Compute PBL height based on bulk Richardson Number
!
!     - calc_ceiling
!       Computes height above MSL, for which cloud coverage > 4/8
!
!     - calpmsl
!       Computes the mean sea-level pressusre
!
!     - calrelhum
!       Computes the relative humidity over water
!
!     - calomega
!       Computes the p-system vertical velocity (p-dot)
!
!     - calprsum
!       Computes the total amount of precipitation.
!
!     - calsnowlmt
!       Calculates height of the snowfall limit.
!
!     - caltopdc
!       Computes the top index of dry convection
!
!     - calhzero
!       Computes the height of the 0-Celsius isotherm
!
!     - calcldepth
!       Computes a normalized cloud depth for TV-presentation
!
!     - calclmod
!       Computes a modified total cloud cover for TV-presentation
!
!     - caliq    
!       Computes the vertically integrated mass of a humidity variable
!       with concentration q
!
!     - calztd
!       Computes the vertically integrated atmospheric dry, wet
!       and total refractivity
!
!     - potential_vorticity_rho
!       Computes the potential vorticity * rho
!
!     - radar_lm_ray
!       for computing grid point values of
!       synthetic radar reflectiviy for the 1mom microphysics schemes
!
!     - lightning_potential_index
!       for computing the LPI after Lynn et al. (2010) with the possibility of
!       an additional filtering of stable conditions to inhibit spurious flash
!       signals from deep orographic wave clouds. The LPI is based on the
!       coexistence of updrafts, supercooled water, graupel/hail and other ice
!       species (cloud ice, snow) in the same grid box in a temperature range 
!       of 0 to -20 deg C.
!
!     - lpi_spatial_filter
!       spatial filtering of the LPI:
!       - set LPI to 0 if the majority of vertical columns in a 10x10 km^2 
!         neighbourhood have MAX(W) < 0.5 m/s.
!       - set LPI to 0 if the smoothed (20x20 km^2) (or optionally unsmoothed)
!         total buoyancy between -50 and -550 hPa above the surface is < -1500 J/kg.
!
!     - max_in_pressure_layers
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49 69 8062 2739
!  fax :   +49 69 8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.32       1999/08/24 Guenther Doms
!  Initial release
! 1.33       1999/10/14 Christina Koepken
!  New routine  'caliq' for calculating the vertical integral of a humidity 
!  variable q (specific humidity or total water).
! 2.11       2001/09/28 Ulrich Schaettler
!  Limit the minimum value of the relative humidity to 0.01 (was 0 before)
! 2.14       2002/02/15 Ulrich Schaettler
!  Bug corrections in the SR calpmsl
!  new SR calsnowlmt for height of snow fall limit as a diagnostic quantity
! 3.3        2003/04/22 Guenther Doms
!  Bug corrction in calclmod.
! 3.5        2003/09/02 Ulrich Schaettler
!  New SR calztd to compute vertically integrated atmospheric dry, wet
!  and total refractivity  (by Guergana Guerova)
! 3.16       2005/07/22 Axel Seifert
!  Introduced new SR radar_lm_ray and gamma_fct for computing radar images
! 3.21       2006/12/04 Daniel Leuenberger, Jochen Foerstner, Axel Seifert
!  Introduction of routine cal_conv_ind, which calculates convective indices
!  Changed some constants in subroutine radar_lm_ray
! V3_23        2007/03/30 Michael Baldauf
!  Added new SR calc_ceiling to compute height above MSL, for which cloud
!  coverage is bigger than 4/8
! V4_1         2007/12/04 Michael Baldauf
!  Added new SR calc_sdi for computation of supercell detection indices
! V4_4         2008/07/16 Ulrich Schaettler
!  Vectorized version of SR cal_conv_ind
!  Changed order of SR arguments for some routines (first dimensions, then fields)
! V4_5         2008/09/10 Ulrich Schaettler
!  Bug correction (wrong dimesion of hhl) in SR calc_ceiling
! V4_8         2009/02/16 Guenther Zaengl, Oliver Fuhrer
!  Adapt SR calclmod for full consistency with new reference atmosphere 
!  implementation
!  Added computation of SWISSxx-Indices to cal_conv_ind
!  Added new SR calc_bulk_richardson and calc_pbl_brn (Daniel Leuenberger)
!  Adjusted interface of vert_avg
! V4_9         2009/07/16 Ulrich Schaettler, Christian Bollmann
!  Implemented a vectorized subroutine version of vert_avg
!  Corrected PRESENT checks in SR ascent (MCH)
! V4_11        2009/11/30 Oliver Fuhrer
!  Initialization of lcomp for lzcalc_si in SR cal_conv_ind
! V4_12        2010/05/11 Michael Baldauf, Oli Fuhrer, Ulrich Schaettler
!  New Subroutine potential_vorticity_rho
!  bug fix in scanning for PBL height in SR calc_pbl_brn
!  Adapted interface for cal_conv_ind according to the needs of FieldExtra
!  Moved SR calc_sdi to module src_output
!  Bug fix in SR calc_bulk_richardson for staggering of winds
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_15        2010/11/19 Ulrich Schaettler
!  Eliminated tgrlat, acrlat from interface to SR potential_vorticity_rho
! V4_18        2011/05/26 Guy deMorsier, Daniel Leuenberger
!  Modifications in SR calsnowlmt and SR calhzero: the search for snowfall and
!  0-degree limit now starts from the top of the model (was from bottom before) (GM)
!  Bug fix in SR cal_conv_ind for computation of cape_3km (a swiss product) (DL)
! V4_21        2011/12/06 Ulrich Blahak, Daniel Leuenberger
!  Fixed allocation/deallocation of zsli in SR cal_conv_ind;
!   also, zsli has to be computed in case of derivation of swiss12-Index
!   even when optional input parameter sli is not present.
!  Updates to the outdated radar_lm_ray for current hydci_pp
!   and hydci_pp_gr (provided by Axel Seifert).
!  Bugfix radar_lm_ray: A factor involving N0_snow ( n0s**(1-p_s) ) was doubled
!   and lead to underestimation of dbz_snow of 40 - 50 dB. This has been corrected.
!  Added reflectivity of cloud ice to radar_lm_ray, based on a
!   monodisperse size distribution with D_i = 200 microns.
!  Added reflectivity of cloud droplets to radar_lm_ray, based on a
!   monodisperse size distribution with D_c = 20 microns.
!  SR calhzero: if no 0-degree isotherm is found, hzerocl is set to missing 
!  value (-999.0) (Daniel Leuenberger)
! V4_24        2012/06/22 Ulrich Schaettler
!  Comment some WRITE statements with idebug in SR radar_lm_ray
! V4_25        2012/09/28 Ulrich Blahak
!  Introduced new subroutine radar_sb_ray for computing radar reflectivities 
!  when using the Seifert-Beheng 2-moment scheme
! V5_1         2014-11-28 Oliver Fuhrer, Ulrich Blahak
!  Replaced ireals by wp (working precision) (OF)
!  Cleaning up the code (unused variables)
!  Added optional argument z_radar_dp to subroutines radar_lm_ray and radar_sb_ray
!     (necessary for radar forward operator) (UB)
!  Optimized gamma_fct() to enable inlining (UB)
! V5_3         2015-10-09 Ulrich Blahak
!  Added routines to calculate the lightning potential index (UB)
! V5_3a        2015-11-24 Jochen Foerstner
!  Added new SR max_in_pressure_layers
! V5_4         2016-03-10 Ulrich Blahak
!  Modified error treatment in LPI computations, to avoid
!   model crashes due to problematic thermodynamic profiles
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
!
USE data_parameters , ONLY :   &
  wp, sp, dp, & ! KIND-type parameters for real variables
  iintegers     ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

IMPLICIT NONE

!------------------------------------------------------------------------------

CONTAINS

!==============================================================================
!==============================================================================

SUBROUTINE cal_conv_ind (te, qve, ue, ve, hsurf, prs_surf, prs, hhl,    &
     idim, jdim, kdim,                                                  &
     b1, b2w, b3, b4w, lh_v, cp_d, r_d, rdv, rvd_m_o, o_m_rdv, g,       &
     missing_value, idebug, lwarn, ierror, yerrmsg,                     &
     cape_mu, cin_mu, cape_ml, cape_3km, cin_ml, si, sli,               &
     swiss00, swiss12, lcl_ml, lfc_ml, lmissing_mask, idiagunit)

!US Arguments: prs_surf, lmissing_mask not used below!
!------------------------------------------------------------------------------
!
! Description:
!  Computation of Convective Available Potential Energy CAPE, 
!  Convective Inhibition CIN, Showalter Index and surface lifted index 
!  based on parcel theory.
! 
! Input:  
!         - Temperature, specific humidity and pressure of environment
!         - several thermodynamical constants
!         - missing value (output wherever a value is not defined)
!         - number of the calling PE
!
! Output: 
!         Several convective indices are computed and output dependent on the 
!         optional arguments of the routine:
!         - cape_mu/cin_mu: CAPE/CIN based on the most unstable parcel in the 
!           (lowest 300hPa of the) PBL
!         - cape_ml/cin_ml: CAPE/CIN based on a parcel with thermodynamical 
!           properties !           of the lowest mean layer in the PBL (50hPa)
!         - Showalter Index (si)
!         - Surface Lifted index (sli)
!         - Swiss Index 00Z (swiss00)
!         - Swiss Index 12Z (swiss12)
!         - lcl_ml/lfc_ml: Lifting Condensation Level/Level of Free Convection
!            based on a parcel with thermodynamical properties of the lowest
!            mean layer in the PBL (50hPa)(ABOVE GROUND level)
!         - cape_3km: CAPE based on a parcel with thermodynamical
!           properties of the lowest mean layer in the PBL (50hPa)with an
!           ascent until 3 km above ground
!      
! Motivation: 
!  Current parameter CAPE_CON is calculated in LM in the framework of the 
!  convective parametrisation scheme. Therefore this parameter is only available
!  at those gridpoints, where the scheme is called, but not continuously on the 
!  whole domain. This subroutine, on the other hand, provides continuous fields. 
!
! Method:
!  A dry/moist parcel ascent is performed following classic parcel theory.
!  Moist adiabatic ascent is calculated iteratively with an appropriate scheme.
!  Based on the temperature and moisture of the ascending parcel, CAPE and CIN
!  are computed, closely following the recommendations of Doswell and Rasmussen 
!  (1994), including a virtual temperature correction and searching for the 
!  most unstable parcel in the lower troposphere. Additionally, a mixed layer 
!  CAPE as well as the traditional Showalter Index and the surface lifted 
!  index are computed as further variables. 
!
!  References used during development: 
!  - C. A. Doswell and Rasmussen, E. N.: The Effect of Neglecting the 
!    Virtual Temperature Correction on CAPE Calculations. 
!    Weather and Forecasting, 9, 625-629.
!
!  - K. A. Emanuel (1994): Atmospheric Convection. Oxford University Press.
!
!  - H. Huntrieser et al. (1997): Comparison of Traditional and Newly Developed 
!    Thunderstorm Indices for Switzerland. Weather and Forecasting, 12, 
!    108-125.
!
!  - D. Bolton (1980): The Computation of Equivalent Potential Temperature. 
!    Monthly Weather Review, 108, 1046-1053
!
!  - Davies, J.M.,2002: On low-level thermodynamic parameters
!    associated with tornadic and nontornadic supercells.
!    Preprints, 21st Conf. On Severe Local Storms, San Antonio, Amer. Meteor. Soc.
!    http://members.cox.net/jondavies1/LLthermo.PDF

! Changes by:
!  Marco Stoll (sto) and Daniel Leuenberger (led), MeteoSwiss 24.7.2006
!  Daniel Leuenberg (led) and Anne Roches (roa), MeteoSwiss 12.4.2008
!
!  Vectorized Version by Uli Schaettler
!
!  Daniel Leuenberger and Anne Roches: SWISS00, SWISS12, cape_3km, heights
!  of the lifting condensation level and the level of free convection added, 6.2008

!------------------------------------------------------------------------------
  
! Input data
!----------- 
INTEGER (KIND=iintegers), INTENT (IN) ::  &
  idim ,                  & ! array dimension in zonal direction
  jdim ,                  & ! array dimension meridional direction
  kdim                      ! array dimension in vertical direction

REAL    (KIND=wp),        INTENT (IN) ::  &
  te  (idim,jdim,kdim),   & ! environment temperature
  qve (idim,jdim,kdim),   & ! environment humidity
  ue  (idim,jdim,kdim),   & ! environment zonal wind speed
  ve  (idim,jdim,kdim),   & ! environment meridional wind speed
  hsurf(idim,jdim),       & ! height of surface topography
  prs_surf(idim,jdim),    & ! surface pressure
  prs (idim,jdim,kdim),   & ! full level pressure
  hhl (idim,jdim,kdim+1)    ! height of half levels


! Physical constants
!-------------------
REAL    (KIND=wp),        INTENT (IN) ::  &
  cp_d,                   & ! specific heat capacity of dry air
  lh_v,                   & ! latent heat of vaporization
  r_d,                    & ! gas constant of dry air
  rdv,                    & ! R_d/R_v
  rvd_m_o,                & ! (R_v/R_d)-1
  o_m_rdv,                & ! 1. - (R_d/R_v)
  g,                      & ! gravity acceleration
  b1,b2w,b3,b4w,          & ! parameters for computing the saturation vapour pressure
  missing_value             ! Missing value for CIN (if no LFC/CAPE was found),
                            ! SI and SLI (if no LFC/CAPE was found)

! Auxiliaries
!-------------------
INTEGER    (KIND=iintegers),    INTENT (IN) ::  &
  idebug                    ! for controlling debug output

INTEGER    (KIND=iintegers),    INTENT(OUT) ::  &
  ierror                    ! for error code setting

CHARACTER (LEN=*),              INTENT(OUT) ::  &
  yerrmsg                   ! for error message

LOGICAL,                        INTENT(OUT) ::  &
  lwarn                     ! to indicate, if warnings have been issued

! Output data
!------------ 
REAL (KIND=wp),     INTENT (OUT), OPTIONAL :: &
  cape_mu  (idim,jdim),   & ! most unstable Convective Available Potential 
                            !                  Energy, CAPE_MU
  cin_mu   (idim,jdim),   & ! (most unstable) Convective INhibition, CIN_MU
  cape_ml  (idim,jdim),   & ! mixed layer CAPE_ML
  cape_3km (idim,jdim),   & ! 3km CAPE based on the mixed layer method (ascent until 3km ABOVE GROUND)
  cin_ml   (idim,jdim),   & ! mixed layer CIN_ML
  lcl_ml   (idim,jdim),   & ! mixed layer Lifting Condensation Level ABOVE GROUND level
  lfc_ml   (idim,jdim),   & ! mixed layer Level of Free Convection ABOVE GROUND level
  si       (idim,jdim),   & ! Showalter Index, SI
  sli      (idim,jdim),   & ! surface lifed index
  swiss00  (idim,jdim),   & ! SWISS00 index
  swiss12  (idim,jdim)      ! SWISS12 index

LOGICAL,                 INTENT (IN),  OPTIONAL :: &
  lmissing_mask(idim,jdim)

INTEGER(KIND=iintegers), INTENT (IN),  OPTIONAL ::  &
  idiagunit                 ! for printing diagnostic messages to a file
           ! NOTE: if used, this file must be open before calling this routine
                                     

! Local scalars and automatic arrays
!-----------------------------------
REAL    (KIND=wp)     ::       &      
  acape    (idim,jdim), & ! CAPE returned from subroutine ascent 
                          ! (the latter performs the parcel ascent)
  acin     (idim,jdim), & ! CIN returned from subroutine ascent
  theta_sum(idim,jdim), & ! help variable to calculate mixed layer average 
                          !    potential temperature
  qve_sum  (idim,jdim), & ! help variable to calculate mixed layer average 
                          !    specific humidity
  q_start  (idim,jdim), & ! parcel initial specific humidity when using mixed 
                          !    layer method
  t_start  (idim,jdim), & ! parcel initial temperature when using mixed 
                          !    layer method

     mup_lay_thck,  & ! thickness of layer in which most unstable parcel 
                      ! is searched for
     theta,         & ! potential temperature
     p0,            & ! reference pressure for calculation of potential 
                      !    temperature (1000hPa)
     ml_depth,      & ! mixed layer depth
     blt1,blt2,blt3,& ! constants used in calculation of dew point
     sistartprs,    & ! start pressure for SI calculation, per definition 850hPa
     sistopprs,     & ! upper limit pressure for SI calculation, per definition 500hPa
     e,             & ! water vapor pressure
     td,            & ! dew point temperature
     hl_l,hl_u,     & ! height of full levels in order to find the closest model level
     vh6000,        & ! norm of the horiz. wind vectors at approx. 6000m
     vh3000,        & ! norm of the horiz. wind vectors at approx. 3000m
     vhfirstlev       ! norm of the horiz. wind vectors at the first model level
  
REAL    (KIND=wp),     ALLOCATABLE ::       &
     zsi(:,:),      & ! local array for si
     zsli(:,:)        ! local array for sli 

INTEGER    (KIND=iintegers) :: &     
  i, j, k2, k1, k,       & ! Indices of input/output fields
  kstart(idim,jdim),     & ! Index of start level used during search of most unstable 
                           !        parcel
  k_ml(idim,jdim),       & ! Index for calculation of mixed layer averages 
                           ! (potential temperature, moisture)
  k_p_mean(idim,jdim),   & ! Model level approximately corresponding to mixed layer 
                           !        mean pressure
  ksi(idim,jdim),        & ! Index for SI calculation
  klcl(idim,jdim),       & ! Indices for Lifting Condensation Level LCL,
  klfc(idim,jdim),       & ! Level of Free Convection LFC and
  elprint (idim,jdim),   & ! Equilibrium Level EL, used for printing statements (debugging)
  k600(idim,jdim),       & ! k index corresponding to the model level closest to 600 hPa
  k650(idim,jdim),       & ! k index corresponding to the model level closest to 650 hPa
  k3000(idim,jdim),      & ! k index corresponding to the model level closest to 3000 m
  k6000(idim,jdim)         ! k index corresponding to the model level closest to 6000 m
 
LOGICAL        :: &
  lzcalc_mu,       & ! calculate most unstable cape/cin
  lzcalc_ml,       & ! calculate mean layer cape/cin
  lzcalc_fc,       & ! calculate level of free convection
  lzcalc_3km,      & ! calculate 3km cape
  lzcalc_si,       & ! showalter index
  lzcalc_sli,      & ! surface lifted index
  lzcalc_swiss00,  & ! SWISS00 index
  lzcalc_swiss12,  & ! SWISS12 index
  lzprint

LOGICAL                     :: &
  lcomp(idim,jdim)         ! to indicate, whether ascent has to be done for 
                           ! MU calculations

CHARACTER (LEN=18)       :: &
  sn                               ! Debugging message
  
CHARACTER (LEN=44)       :: &
  msg_mup,                  &      ! Some further debugging messages
  msg_si
  
! The following parameters are help values for the iterative calculation 
! of the parcel temperature during the moist adiabatic ascent
REAL    (KIND=wp)             :: r1,r2,esat,tguess1,tguess2,thetae1,thetae2
REAL    (KIND=wp), PARAMETER  :: eps_conv=0.03_wp

! Statement functions
! -------------------
REAL    (KIND=wp)        ::    &
  ztx,zpx,zgex,zrx,ze, &
  fesatw,              & ! Function for equilibrium vapour pressure over water
  fqvs,                & ! Function for saturation specific humidity
  fthetae,             & ! Function for equivalent potential temperature
  ftd                    ! Function for dewpoint

! Calculation of saturation specific humidity and equivalent potential 
! temperature, based upon Bolton's approximations (1980), see also Emanuel 
! (1994.)
! Sat.vap. pressure over water
fesatw(ztx)     = b1*EXP( b2w*(ztx-b3)/(ztx-b4w) )       

! Specific humidity at saturation
fqvs(zgex,zpx)  = rdv*zgex/( zpx - o_m_rdv*zgex )        
  
! Equivalent potential temperature to hold constant during ascent
fthetae(ztx,zpx,zrx) = (p0/zpx)**(r_d/cp_d)*ztx*exp(lh_v/cp_d*zrx/ztx)  
  
! Dew point temperature (based on Magnus-Tetens formula)
ftd (ze)    = ( b3*b2w -log(ze/b1)*b4w ) / ( b2w -log(ze/b1) )

!------------------------------------------------------------------------------
! Begin Subroutine cal_conv_ind
!------------------------------------------------------------------------------
  
! Initialize variables
!---------------------

  ierror         = 0_iintegers
  yerrmsg        = '          '
  lwarn          = .FALSE.

  lzcalc_mu      = PRESENT(cape_mu) .OR. PRESENT(cin_mu)
  lzcalc_ml      = PRESENT(cape_ml) .OR. PRESENT(cin_ml) .OR. &
                   PRESENT(lcl_ml)  .OR. PRESENT(lfc_ml) .OR. PRESENT(cape_3km)
  lzcalc_fc      = PRESENT(lcl_ml)  .OR. PRESENT(lfc_ml)
  lzcalc_3km     = PRESENT(cape_3km)
  lzcalc_si      = PRESENT(si)      .OR. PRESENT(swiss00)
  lzcalc_sli     = PRESENT(sli)     .OR. PRESENT(swiss12)
  lzcalc_swiss00 = PRESENT(swiss00)
  lzcalc_swiss12 = PRESENT(swiss12)
  lzprint        = PRESENT(idiagunit)

  IF (lzcalc_mu) THEN
    cape_mu(:,:)  = 0.0_wp 
    cin_mu (:,:)  = missing_value     
    ! Set to missing value, if no CAPE is present then CIN not defined!
  ENDIF

  acape(:,:)     = 0.0_wp
  acin (:,:)     = missing_value

  IF (lzcalc_ml) THEN
     cape_ml(:,:)  = 0.0_wp
     cin_ml(:,:)   = missing_value
  ENDIF
  IF (lzcalc_3km) THEN
     cape_3km(:,:) = 0.0_wp
  ENDIF
  IF (lzcalc_si) THEN
     ALLOCATE(zsi(idim,jdim))
     zsi(:,:)      = missing_value
  ENDIF
  IF (lzcalc_sli .OR. lzcalc_swiss12) THEN
     ALLOCATE(zsli(idim,jdim))
     zsli(:,:)     = missing_value
  END IF

  ! Reference pressure for potential temperature calculation
  p0             = 100000._wp    

! ! Mixed layer CAPE, real variables:
! q_start(:,:)   = 0.0_wp    
! t_start(:,:)   = 0.0_wp    

  ! Depth of mixed surface layer: 50hPa following Huntrieser, 1997.
  ! Other frequently used value is 100hPa.
  ml_depth       = 5000._wp      
  
  ! Level variables needed for profile print (1-D standalone version) and 
  ! debugging statements:
  klcl(:,:)   = 0_iintegers
  klfc(:,:)   = 0_iintegers
  elprint (:,:)   = 0_iintegers
  
  ! Pressure levels limiting the ascent during the calculation of Showalter 
  ! Index
  sistartprs = 85000.0_wp
  sistopprs  = 50000.0_wp
  
  ! Constants applied in Bolton's formula for dewpoint etc. calculation
  blt1 = 243.5_wp
  blt2 = 17.67_wp
  blt3 = 6.112_wp
  
  ! Thickness of the layer within which the most unstable parcel 
  ! is searched for -> smaller values of this parameter result in
  ! less computational cost, adapt if the code is running very slowly!
  mup_lay_thck = 30000._wp
  
  ! Print messages, use for debugging
  sn      = 'cal_conv_ind: '
  msg_mup = '     most unstable parcel updated: mucape = '
  msg_si  = '     start ascent for calculation of SI'
    
  
!------------------------------------------------------------------------------
! Start computation. Overview:
! 
! First calculate CAPE_MU/CIN_MU, then CAPE_ML/CIN_ML, then SI, then SLI
! at every single gridpoint. i and j loops go horizontally 
! through the whole domain, parcelloop varies the start level
! of the parcel in most unstable calculation method.
! One single ascent of a parcel is calculated within the 
! SUBROUTINE ascent further below. 
! Here we make sure that "ascent" is being called with the 
! correct initial temperature, moisture, model level and 
! gridpoint indices. 
!------------------------------------------------------------------------------
  
!US  ! horizontal loops
!US  jloop:DO j = 1, jdim
!US    iloop:DO i = 1, idim 
      
!------------------------------------------------------------------------------
! Section 1: Organize the calculation of most unstable CAPE
!
! The initial model level, from where the parcel starts ascending is varied 
! until the pressure at this level is lower than the threshold mup_lay_thck, 
! defined above.
!------------------------------------------------------------------------------
        
  IF (lzcalc_mu) THEN

    IF (idebug > 50) &
      PRINT *, sn,'entering calculations for CAPE/CIN (most unstable parcel)'

    parcelloop:  DO k1 = kdim, 1, -1

      DO j = 1, jdim
        DO i = 1, idim
          kstart(i,j) = k1
          lcomp (i,j) = (prs(i,j,k1) > (prs(i,j,kdim)-mup_lay_thck))
        ENDDO
      ENDDO

      IF (idebug > 50) WRITE(*,'(A,A,I2)') &
        sn,'   type of cape calculation: most unstable, kstart = ',kstart(1,1)
              
      ! Take temperature and moisture of environment profile at current 
      ! level as initial values for the ascending parcel, call "ascent"
      ! to perform the dry/moist adiabatic parcel ascent and get back
      ! the calculated CAPE/CIN
      IF (idebug > 500) THEN
        CALL ascent   (idim, jdim, kstart, te(:,:,k1), qve(:,:,k1),      &
                       lcomp(:,:), acape=acape, acin=acin)
      ELSE
        CALL ascent   (idim, jdim, kstart, te(:,:,k1), qve(:,:,k1),      &
                       lcomp(:,:), acape=acape, acin=acin, ael=elprint,  &
                       alfc=klfc, alcl=klcl)
      ENDIF
          
      ! current model level kstart is greater than the CAPE computed 
      ! before. If yes then save the CAPE from this ascent in help 
      ! variable mucape as the most unstable one, that has been found
      ! so far. 
      DO j = 1, jdim
        DO i = 1, idim
          IF (lcomp(i,j)) THEN

            ! debugging print statement 
            IF (idebug > 500) WRITE(*,'(A,A,I2,A,F6.1,A,F6.1)') sn,       &
             ' kstart = ',kstart(i,j),' acape = ',acape(i,j),' acin = ',acin(i,j)
      
            ! Check if CAPE returned from the parcel ascent starting at 

            IF ( acape(i,j) > cape_mu(i,j) ) THEN
              cape_mu(i,j)  = acape(i,j)
              cin_mu (i,j)  = acin (i,j)
         
              ! debugging print statement
              IF (idebug > 500) THEN 
                WRITE(*,'(A,A,F6.1,A,F6.1,A,F8.1)') &
                     sn, msg_mup, cape_mu(i,j),' mucin = ',cin_mu(i,j),&
                     ' startprs =',prs(i,j,kstart(i,j))/100.0_wp
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      
    ENDDO parcelloop

    ! Top of layer for most unstable parcel search has been reached, 

    ! write out the results
    IF (idebug > 500) THEN
      DO j = 1, jdim
        DO i = 1, idim
          WRITE (*,'(A,A,2I3,A,I2,A,F6.1,A,I2,A,I2,A,I2)')            &
             sn,'     i,j = ', i,j, '  kstart = ', kstart,            &
                '  mucape = ', cape_mu(i,j), ' LCL = ', klcl(i,j),&
                ' LFC = ', klfc(i,j), ' EL = ', elprint(i,j)
        ENDDO
      ENDDO
    ENDIF
        
  ! end most unstable method
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Organize calculation of mixed layer CAPE
! 
! A well mixed near surface layer is assumed (its depth is specified with 
! parameter ml_depth) Potential temperature and specific humidity are constant
! in this layer, they are calculated as arithmetical means of the corresponding
! variables of the environment (model) profile. The parcel starts from a level 
! approximately in the middle of this well mixed layer, with the average spec. 
! humidity and potential temperature as start values. 
!
! Meaning of the CAPE_3KM
! -----------------------
! When the CAPE_3KM is high, it tends to promote a stretching of the air columns
! and thus promote the development of tornadoes when associated with a high
! vertical vorticity near the ground and with a high humidity. Significant
! stretching when CAPE_3KM > 150 J/kg. (see Davies, 2002)
!
!------------------------------------------------------------------------------
        
  IF (lzcalc_ml) THEN

    IF (idebug > 50) &
      PRINT *, sn,'entering calculations for CAPE/CIN (mean layer parcel)'

    ! debugging print statement
    IF (idebug > 50) WRITE(*,'(A,A,F6.0)') &
       sn,'   type of cape calculation: mixed layer, ml_depth = ',ml_depth
            
    ! reinitialize the help variables first
    k_ml    (:,:)  = kdim  ! index used to step through the well mixed layer
    k_p_mean(:,:)  = kdim  ! index of model level corresponding to average 
                           ! mixed layer pressure

    qve_sum  (:,:) = 0.0_wp ! sum of specific humidities in well mixed layer
    theta_sum(:,:) = 0.0_wp ! sum of potential temperatures in well mixed layer
            
    ! now calculate the mixed layer average potential temperature and 
    ! specific humidity
    mixedlayerloop: DO k2 = kdim, 1, -1
      DO j = 1, jdim
        DO i = 1, idim
          IF (prs(i,j,k2) > (prs(i,j,kdim) - ml_depth)) THEN
            qve_sum(i,j)   = qve_sum(i,j) + qve(i,j,k2)
            theta          = te(i,j,k2)*(p0/prs(i,j,k2))**(r_d/cp_d)
            theta_sum(i,j) = theta_sum(i,j) + theta
               
            ! Find the level, where pressure approximately corresponds to the 
            ! average pressure of the well mixed layer. Simply assume a threshold
            ! of ml_depth/2 as average pressure in the layer, if this threshold 
            ! is surpassed the level with approximate mean pressure is found
            IF (prs(i,j,k2) > prs(i,j,kdim) - ml_depth/2.0_wp) THEN
              k_p_mean(i,j) = k2
            ENDIF

            k_ml(i,j) = k2 - 1
          ENDIF
        ENDDO
      ENDDO     
    ENDDO mixedlayerloop
        
    ! Calculate the start values for the parcel ascent, 
    ! temperature at level with mean pressure has to be 
    ! calculated from potential temperature:
    DO j = 1, jdim
      DO i = 1, idim
        q_start(i,j) =  qve_sum  (i,j) / (kdim-k_ml(i,j))
        t_start(i,j) = (theta_sum(i,j) / (kdim-k_ml(i,j)))    &
                       * (prs(i,j,k_p_mean(i,j))/p0)**(r_d/cp_d)
        lcomp(i,j)   = .TRUE.
      ENDDO
    ENDDO     
        
    ! do the ascent
    IF (idebug > 500) THEN
      IF (lzcalc_3km) THEN
        CALL ascent (idim, jdim, k_p_mean, t_start, q_start,               &
                     lcomp, acape=cape_ml, acin=cin_ml, acape3km=cape_3km, &
                     ael=elprint, alfc=klfc, alcl=klcl)
      ELSE
        CALL ascent (idim, jdim, k_p_mean, t_start, q_start,               &
                     lcomp, acape=cape_ml, acin=cin_ml,                    &
                     ael=elprint, alfc=klfc, alcl=klcl)
      ENDIF
    ELSE
      IF (lzcalc_3km) THEN
        CALL ascent (idim, jdim, k_p_mean, t_start, q_start,               &
                     lcomp, acape=cape_ml, acin=cin_ml, acape3km=cape_3km, &
                     alfc=klfc, alcl=klcl)
      ELSE
        CALL ascent (idim, jdim, k_p_mean, t_start, q_start,               &
                     lcomp, acape=cape_ml, acin=cin_ml,                    &
                     alfc=klfc, alcl=klcl)
      ENDIF
    ENDIF
        
    ! Compute LCL and LFC using klcl and klfc

    IF (lzcalc_fc) THEN
      DO j= 1, jdim
        DO i= 1, idim

          ! Lifting condensation level computations ABOVE GROUND level

          IF ((klcl(i,j).LE.0 ).OR.(klcl(i,j).GT.kdim))THEN ! if index not defined
             lcl_ml(i,j)=missing_value                              ! lcl=missing_val
          ELSE
             lcl_ml (i,j) = 0.5_wp * ( hhl(i,j,klcl(i,j)) + hhl(i,j,klcl(i,j)+1) )-hsurf(i,j)
          ENDIF

          ! Level of free condensation ABOVE GROUND level

          IF ((klfc(i,j).LE.0).OR.(klfc(i,j).GT.kdim))THEN ! if index not defined
             lfc_ml(i,j)=missing_value                             ! lfc=missing_value
          ELSE
             lfc_ml (i,j) = 0.5_wp * ( hhl(i,j,klfc(i,j)) + hhl(i,j,klfc(i,j)+1) )-hsurf(i,j)
          ENDIF

        ENDDO
      ENDDO
    ENDIF

    ! debugging print statement
    IF (idebug > 500) THEN
      IF (lzcalc_3km) THEN
        DO j = 1, jdim
          DO i = 1, idim
            WRITE (*,'(A,A,2I3,A,F6.1,A,F6.1,A,I2,A,I2,A,I2)')          &
               sn,'     i,j = ', i,j, '  cape_ml = ', cape_ml(i,j),     &
               ' cin_ml = ', cin_ml(i,j), 'cape _3km=', cape_3km(i,j),  &
               ' LCL = ', klcl(i,j),' LFC = ', klfc(i,j), ' EL = ', elprint(i,j)
          ENDDO
        ENDDO
      ELSE
        DO j = 1, jdim
          DO i = 1, idim
            WRITE (*,'(A,A,2I3,A,F6.1,A,F6.1,A,I2,A,I2,A,I2)')          &
               sn,'     i,j = ', i,j, '  cape_ml = ', cape_ml(i,j),     &
               ' cin_ml = ', cin_ml(i,j),                               &
               ' LCL = ', klcl(i,j),' LFC = ', klfc(i,j), ' EL = ', elprint(i,j)
          ENDDO
        ENDDO
      ENDIF
    ENDIF
        
  ! end mixed layer method
  ENDIF
        
!------------------------------------------------------------------------------
! Section 3: Showalter Index SI calculation
!
! Definition: Tp - Te at 500hPa, where Tp is temperature of parcel, ascending 
! from start level 850hPa and Te is environment temperature.
! Implementation here is done straightforward based on this definition.
!------------------------------------------------------------------------------

  IF (lzcalc_si) THEN

    IF (idebug > 50)                                       &
         PRINT *, sn,'entering calculations for showalter index'

    IF (idebug > 50) WRITE(*,'(A,A)') sn, msg_si
           
    ! Search a model level where pressure is approximately 850hPa
    ksi  (:,:) = kdim
    lcomp(:,:) = .TRUE.    ! fix by Oli Fuhrer: must be initialized

    siloop: DO k2 = kdim, 1, -1
      DO j = 1, jdim
        DO i = 1, idim
          IF (prs(i,j,k2) >=  sistartprs) THEN
            ksi(i,j) = k2           ! but this 1 has been added below:::  ksi(i,j)-1
          ENDIF
        ENDDO
      ENDDO
    ENDDO siloop

    DO j = 1, jdim
      DO i = 1, idim
        t_start(i,j) = te  (i,j,ksi(i,j))
        q_start(i,j) = qve (i,j,ksi(i,j))
      ENDDO
    ENDDO

    IF (idebug > 500) THEN
      DO j = 1, jdim
        DO i = 1, idim
          WRITE(*,'(A,A,I2)') sn,'     850hPa level found at k = ',ksi(i,j)-1
        ENDDO
      ENDDO
    ENDIF
           
    ! Approximate 850hPa level is reached, or lowest model level pressure
    ! is below 850hPa - in the latter case SI is set to a missing value
    CALL ascent   (idim, jdim, ksi, t_start, q_start, lcomp, asi=zsi)

    DO j = 1, jdim
      DO i = 1, idim
        IF ( prs(i,j,kdim) < sistartprs ) THEN 
          zsi(i,j) = missing_value
        ENDIF
      ENDDO
    ENDDO
           
    IF (idebug > 500) THEN
      DO j = 1, jdim
        DO i = 1, idim
          WRITE (*,'(A,A,F8.2)') sn,' SI calculation done, si =',zsi(i,j)
        ENDDO
      ENDDO
    ENDIF

  ! end SI calculation
  ENDIF
        
!------------------------------------------------------------------------------
! Section 4: surface lifed index SLI calculation
!
! Definition: Tp - Te at 500hPa, where Tp is temperature of parcel, ascending 
! from lowest level and Te is environment temperature.
! Implementation here is done straightforward based on this definition.
!------------------------------------------------------------------------------

  IF (lzcalc_sli .OR. lzcalc_swiss12) THEN

    IF (idebug > 50)                                                 &
         PRINT *, sn,'entering calculations for surface lifted index'
           
    IF (idebug > 50) WRITE (*,'(A,A)') sn,msg_si
    
    ! Calculate Tp - Te at 500hPa
    kstart (:,:) = kdim
    lcomp  (:,:) = .TRUE.

    CALL ascent   (idim, jdim, kstart, te(:,:,kdim), qve(:,:,kdim), lcomp, asi=zsli)
           
    IF (idebug > 500) THEN
      DO j = 1, jdim
        DO i = 1, idim
          WRITE (*,'(A,A,F8.2)') sn,' SI calculation done, sli =',zsli(i,j)
        ENDDO
      ENDDO
    ENDIF

  ! end SLI calculation
  ENDIF

    !------------------------------------------------------------------------------
    ! Section 5: SWISS indices
    !
    ! Definition: SWISS00 and SWISS12 are two statistically based indices developped
    !             using data originating from the meteorological station of
    !             Payerne in Switzerland. SWISS00 has been developed with observations at
    !             00:00UTC whereas SWISS12 has been develped with observations at 12:00UTC.
    !------------------------------------------------------------------------------

    IF (lzcalc_swiss00 .OR. lzcalc_swiss12) THEN

       !Initialisation of vertical indexes in order to treat correctly the points
       !located below 600hPa, below 650hPa, above 3000m and above 6000m resp.
       !(missing_value for SWISS indices)
       k600(:,:)  = -1
       k650(:,:)  = -1
       k3000(:,:) = -1
       k6000(:,:) = -1

       DO k = kdim, 2, -1
          DO j = 1,jdim
             DO i = 1,idim

                !Find the k index coressponding to the model level closest to 600 hPa
                IF ((prs (i,j,k) >= 60000._wp).AND.(prs (i,j,k-1) <= 60000._wp))THEN
                   IF (abs(prs (i,j,k)- 60000._wp) <= abs(prs(i,j,k-1)- 60000._wp))THEN
                      k600(i,j) = k
                   ELSE
                      k600(i,j) = k-1
                   ENDIF
                ENDIF

                !Find the k index coressponding to the model level closest to 650 hPa
                IF ((prs (i,j,k) >= 65000._wp).AND.(prs (i,j,k-1) <= 65000._wp))THEN
                   IF (abs(prs (i,j,k)- 65000._wp) <= abs(prs(i,j,k-1)- 65000._wp))THEN
                      k650(i,j) = k
                   ELSE
                      k650(i,j) = k-1
                   ENDIF
                ENDIF

                hl_l = 0.5_wp * (hhl(i,j,k  ) + hhl(i,j,k+1))
                hl_u = 0.5_wp * (hhl(i,j,k-1) + hhl(i,j,k  ))
                !Find the k index corresponding to the model level closest to 3000 m
                IF ((hl_l  <=  3000._wp) .AND. (hl_u  >=  3000._wp)) THEN
                   IF (abs(hl_l - 3000._wp) <= abs(hl_u - 3000._wp)) THEN
                      k3000(i,j) = k
                   ELSE
                      k3000(i,j) = k-1
                   ENDIF
                ENDIF

                !Find the k index corresponding to the model level closest to 6000 m
                IF ((hl_l  <=  6000._wp) .AND. (hl_u  >=  6000._wp)) THEN
                   IF (abs(hl_l - 6000._wp) <= abs(hl_u - 6000._wp)) THEN
                      k6000(i,j) = k
                   ELSE
                      k6000(i,j) = k-1
                   ENDIF
                ENDIF

             ENDDO
          ENDDO
       ENDDO
    ENDIF

    !------------------------------------------------------------------------------
    ! Section 5a: SWISS00 index
    !
    !             Its components are:
    !               - the showalter index (see above)
    !               - the wind shear between 3000m and 6000m
    !               - the dew point depression at 600hPa
    !
    !             A SWISS00 value less than 5.1 means "likely thunderstorms".
    !------------------------------------------------------------------------------

    IF (lzcalc_swiss00) THEN

       IF (idebug > 50) &
            print *, sn,'entering calculations for swiss00 index'

       DO j = 1,jdim
          DO i = 1,idim

             IF ( (k600(i,j)  <  0) .OR. (k3000(i,j)  <  0) .OR. (k6000(i,j)  <  0)  &
                  .OR.(zsi (i,j) == missing_value) ) THEN

                swiss00(i,j) = missing_value

             ELSE

                !Compute vapor pressure at approximately 600 hPa using index k600

                e = prs(i,j,k600(i,j))*qve(i,j,k600(i,j))/(qve(i,j,k600(i,j))+rdv*(1-qve(i,j,k600(i,j))))

                IF (e <= 0.0_wp) THEN

                   swiss00(i,j) = missing_value

                ELSE

                   !Compute dew point temperature at approximately 600 hPa
                   !Obtained using this identity: e = fesatw(td) at approximately 600 hPa
                   !It means that the vapor pressure is equal to the saturation vapor
                   !pressure at the dew point temperature

                   td = ftd (e)

                   !Compute norm of horizontal wind vectors at approximately 6000 m
                   !using k6000 index

                   vh6000 = sqrt ( ue (i,j,k6000(i,j)) **2 + ve (i,j,k6000(i,j)) **2)

                   !Compute norm of horizontal wind vectors at approximately 3000 m
                   !using k3000 index

                   vh3000 = sqrt ( ue (i,j,k3000(i,j)) **2 + ve (i,j,k3000(i,j)) **2)


                   !Compute SWISS00 Index with the above calculated parameters and the
                   !showalter index

                   swiss00(i,j) = zsi(i,j) + 0.4_wp*(vh6000-vh3000) + 0.1_wp*MAX(0.0_wp,(te(i,j,k600(i,j))-td))

                END IF

             ENDIF

             IF (idebug > 50) write(*,'(a,a,F8.2)') sn,' swiss00 calculation done, swiss00 =',swiss00(i,j)
          ENDDO
       ENDDO


       !end swiss00 calculation

    END IF

    !------------------------------------------------------------------------------
    ! Section 5b: SWISS12 index
    !
    !             Its components are:
    !               - the surface lifted index (see above)
    !               - the wind shear between the ground and 3000m
    !               - the dew point depression at 650hPa
    !
    !             The concept of "wind at the ground" is ambiguous: the wind at the
    !             surface is equal to zero but we could also consider the wind at
    !             10m for example as "ground value". We use the wind at the first
    !             model level as an approximation for the surface wind.
    !
    !             A SWISS12 value less than 0.6 means "likely thunderstorms".
    !------------------------------------------------------------------------------

    IF (lzcalc_swiss12) THEN

       IF (idebug > 50) &
            print *, sn,'entering calculations for swiss12 index'

       DO j = 1,jdim
          DO i=1,idim

             IF ( (k650(i,j)  <  0) .OR. (k3000(i,j)  <  0) ) THEN
                swiss12(i,j) = missing_value
             ELSE


                !Compute vapor pressure at approximately 650 hPa using index k650

                e = prs(i,j,k650(i,j))*qve(i,j,k650(i,j))/(qve(i,j,k650(i,j))+rdv*(1 - qve(i,j,k650(i,j))))

                !Compute dew point temperature at approximately 650 hPa
                !Obtained using this identity: e = fesatw(td) at approximately 650 hPa
                !It means that the vapor pressure is equal to the saturation vapor
                !pressure at the dew point temperature

                IF (e <= 0.0_wp) THEN
                   swiss12(i,j) = missing_value
                ELSE

                   td = ftd (e)

                   !Compute norm of horizontal wind vectors  at approximately 3000 m

                   vh3000 = sqrt ( ue (i,j,k3000(i,j)) **2 + ve (i,j,k3000(i,j)) **2)

                   !Compute norm of horizontal wind vectors  at the first model level
                   ! WARNING: Currently the first model level is at 10m a.s.l. but this
                   !          height could change!!!
                   vhfirstlev = sqrt ( ue (i,j,kdim) **2 + ve (i,j,kdim) **2)

                   !Compute SWISS12 Index with the above calculated parameters and the
                   !surface lifted index

                   swiss12(i,j) = zsli(i,j) - 0.3_wp*(vh3000-vhfirstlev) + 0.3_wp*(MAX(0.0_wp,(te(i,j,k650(i,j))-td)))

                END IF

             ENDIF


             IF (idebug > 50) write(*,'(a,a,F8.2)') sn,' swiss12 calculation done, swiss12 =',swiss12(i,j)
             !end swiss12 calculation
          ENDDO
       ENDDO

    END IF

    ! store local arrays of si and sli in optional arguments, if required
    IF (PRESENT(si)) THEN
       si = zsi
    END IF
    IF (ALLOCATED(zsi)) DEALLOCATE(zsi)

    IF (PRESENT(sli)) THEN
       sli = zsli
    END IF
    IF (ALLOCATED(zsli)) DEALLOCATE(zsli)

!==============================================================================

CONTAINS

  SUBROUTINE ascent ( idim, jdim, kstart, tp_start, qvp_start, lcomp,   & ! in
                   acape, acin, acape3km, ael, alfc, alcl, asi) ! out
!US removed tp_print, because it is not used and the implementation is buggy
!US                acape, acin, acape3km, ael, alfc, alcl, asi, tp_print) ! out
  
  !------------------------------------------------------------------------------
  !
  ! Description:
  !   Vectorized version of Subroutine ascent: computations are done
  !   concurrently for every grid point. Several debug prints have to be
  !   eliminated (only commented out), because they cannot be vectorized.
  !   (Ulrich Schaettler)
  !   A single parcel ascent is performed, based on the given start 
  !   values kstart (level), tp_start (initial parcel temperature) and
  !   qvp_start (initial parcel specific humidity). 
  !   CAPE, CIN, LCL, LFC, EL and SI are returned as optional arguments.
  !
  !------------------------------------------------------------------------------
      
  INTEGER, INTENT (IN) :: idim, jdim
  INTEGER, INTENT (IN) :: kstart(idim,jdim)   ! gridpoint and startlevel indices
      
  REAL (KIND=wp),     INTENT(IN)  ::       &
    tp_start (idim,jdim),    & ! initial temperature
    qvp_start(idim,jdim)       ! initial specific humidity
      
  LOGICAL,            INTENT(IN)  ::       &
    lcomp(idim,jdim)           ! to indicate, whether computations have to be done

  REAL (KIND=wp),     INTENT(OUT), OPTIONAL ::  &
    acape(idim,jdim),        & ! CAPE
    acin (idim,jdim),        & ! CIN             to be returned
    acape3km(idim,jdim),     & ! CAPE 3km
    asi  (idim,jdim)           ! SI
!US    tp_print(kdim)             ! for debug prints, but is not used
      
  INTEGER, INTENT(OUT), OPTIONAL ::        &
    ael (idim,jdim),        & ! level indices of EL
    alfc(idim,jdim),        & !                  LFC
    alcl(idim,jdim)           !              and LCL
      
  INTEGER ::                &
    i,j, idone,             & !
    k,                      & !
    k3000(idim,jdim),       & ! index k corresponding to 3000m
    ellev(idim,jdim),       & !
    lfclev(idim,jdim),      & !
    lcllev(idim,jdim),      & ! local EL, LFC and LCL indices
    lfcfound(idim,jdim),    & ! flag indicating if a LFC has already been found
                              ! below, in cases where several EL and LFC's occur
    icount  (idim,jdim)       ! counter for the iterative process
      
  REAL (KIND=wp)     ::                    & 
    ttt      (idim,jdim),       & ! for the iterative process
    thp_start(idim,jdim),       & ! potential temperature of parcel at start
    cape     (idim,jdim),       & ! CAPE computed in this parcel ascent
    cape3km  (idim,jdim),       & ! CAPE 3km computed in this parcel ascent
    cin      (idim,jdim),       & ! CIN computed in this parcel ascent
    cin_help (idim,jdim),       & ! help variable, the CIN above the LFC
    si       (idim,jdim),       & ! SI computed in this parcel ascent
    buo      (idim,jdim),       &   ! parcel buoyancy at level k
    tp       (idim,jdim,kdim),  & ! temperature profile of ascending parcel
    qvp      (idim,jdim,kdim),  & ! specific moisture profile of ascending parcel
    rp       (idim,jdim,kdim),  & ! mixing ratio profile of ascending parcel
    tvp,            &   ! virtual temperature of parcel at level k
    tve,            &   ! virtual temperature of environment at level k
    buo_belo,       &   ! parcel buoyancy of level k+1 below
    esatp,          &   ! saturation vapour pressure at level k
    qvsp,           &   ! saturation specific humidity at level k
    thp(idim,jdim), &   ! 1st guess theta_e of parcel for iterative 
                        ! calculation of moist adiabatic ascent
    cc_comp,        &   ! parameter used to weight CAPE and CIN amount
                        ! in LFC determination
    hl_l,           &   ! lowest interval level
    hl_u                ! highest interval level

  LOGICAL            :: &
    ldone    (idim,jdim),       & ! to indicate already processed points
    lzcompute(idim,jdim)          ! to indicate points for computation

  CHARACTER (len=20) :: &
    sn
      
  !------------------------------------------------------------------------------
  ! Begin Subroutine ascent
  !------------------------------------------------------------------------------
  
  ! Initialization
  
    sn = 'ascent:  '
    lcllev  (:,:)   = 0_iintegers
    lfclev  (:,:)   = 0_iintegers
    ellev   (:,:)   = 0_iintegers
    lfcfound(:,:)   = 0_iintegers
    cape    (:,:)   = 0.0_wp
    cin     (:,:)   = 0.0_wp
    cape3km (:,:)   = 0.0_wp
    cin_help(:,:)   = 0.0_wp
    si      (:,:)   = missing_value                   

    tp (:,:,:)      = 0.0_wp
    rp (:,:,:)      = 0.0_wp
    qvp(:,:,:)      = 0.0_wp               
  
    buo(:,:)   = 0.0_wp
    buo_belo   = 0.0_wp

    IF (PRESENT(acape3km)) THEN
       k3000(:,:)= -1_iintegers
       ! computation of the index corresponding to 3000m for the CAPE 3km computation
       DO k = kdim, 2, -1
          DO j = 1,jdim
             DO i = 1,idim
                hl_l = 0.5_wp * (hhl(i,j,k  ) + hhl(i,j,k+1)) - hsurf(i,j)
                hl_u = 0.5_wp * (hhl(i,j,k-1) + hhl(i,j,k  )) - hsurf(i,j)
                !Find the k index corresponding to the model level closest to 3000 m
                IF ((hl_l  <=  3000._wp) .AND. (hl_u  >=  3000._wp)) THEN
                   IF (abs(hl_l - 3000._wp) <= abs(hl_u - 3000._wp)) THEN
                      k3000(i,j) = k
                   ELSE
                      k3000(i,j) = k-1
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    ! this parameter helps to find the LFC above a capping inversion in cases, 
    ! where a LFC already was found in an unstable layer in the convective 
    ! boundary layer below. 
    cc_comp    = 2.0_wp                      
      
    ! potential temperature of initial parcel
    DO j = 1, jdim
      DO i = 1, idim
        thp_start(i,j)       = tp_start (i,j) * (p0/prs(i,j,kstart(i,j)))**(r_d/cp_d)  
        rp (i,j,kstart(i,j)) = qvp_start(i,j) / (1.0_wp-qvp_start(i,j))
      ENDDO
    ENDDO
  
    ! Loop over all model levels above kstart
    ! vectorization: this loop should start from kstart, but this now depends on i,j
    !                therefore it is done for all k, and inside the i-,j-loops
    !                it is checked, whether computations have to be done

    kloop: DO k = kdim, 1, -1

      DO j = 1, jdim
        DO i = 1, idim
          ! initialize the counter used within the iterations
          icount(i,j) = 0_iintegers    ! iterations counter
          ldone (i,j) = .FALSE.
        ENDDO
      ENDDO

      DO j = 1, jdim
        DO i = 1, idim
         
          computation_IF: IF ( (k <= kstart(i,j)) .AND. (lcomp(i,j)) ) THEN

            ! Dry ascent if below cloud base, assume first level is not saturated 
            ! (first approximation)
            IF (k > lcllev(i,j)) THEN
              tp (i,j,k) = thp_start(i,j)*(prs(i,j,k)/p0)**(r_d/cp_d)  ! dry adiabatic process
              qvp(i,j,k) = qvp_start(i,j)                              ! spec humidity conserved
              rp (i,j,k) = rp(i,j,kstart(i,j))                         ! mixing ratio conserved
            
              ! Calculate parcel saturation vapour pressure and saturation 
              ! specific humidity
              esatp = fesatw(tp(i,j,k))
              qvsp  = fqvs(esatp,prs(i,j,k))
            
              ! Check whether parcel is saturated or not and 
              ! no LCL was already found below
              IF ( (qvp(i,j,k) >= qvsp) .AND. (lcllev(i,j) == 0) ) THEN  
                lcllev(i,j) = k                                    ! LCL is reached
!               IF (idebug > 500) WRITE(*,'(A,A15,I2,A,F8.1)') &
!                     sn,'LCL at k = ',lcllev(i,j),' p = ',prs(i,j,lcllev(i,j))/100.0_wp
              ENDIF
            ENDIF
         
            ! Moist ascent above LCL, first calculate an approximate thetae to hold 
            ! constant during the remaining ascent
            IF ( k == lcllev(i,j) ) THEN
              thp(i,j) = fthetae( tp(i,j,k),prs(i,j,k),rp(i,j,k) )
            ENDIF
         
            ! Moist adiabatic process: the parcel temperature during this part of 
            ! the ascent is calculated iteratively using the iterative newton
            ! scheme, assuming the equivalent potential temperature of the parcel 
            ! at the LCL (thp) is held constant. The scheme converges usually within
            ! few (less than 10) iterations, its accuracy can be tuned with the 
            ! parameter "eps_conv", a value of 0.03 is tested and recommended. 
            
            IF ( k <= lcllev(i,j) ) THEN                                
              ! The scheme uses a first guess temperature, which is the parcel 
              ! temperature at the level below. If it happens that the initial 
              ! parcel is already saturated, the environmental temperature 
              ! is taken as first guess instead (tfg). 
              IF (  k == kstart(i,j) ) THEN
                ttt(i,j) = te(i,j,kstart(i,j))            
!               IF (idebug > 500) WRITE(*,'(A,A)') sn,'initial parcel is already saturated'
              ELSE
                ttt(i,j) = tp(i,j,k+1)
              END IF
            ENDIF
            
          ENDIF computation_IF
        ENDDO
      ENDDO

      DO j = 1, jdim
        DO i = 1, idim
          IF  ( (k <= kstart(i,j)) .AND. (lcomp(i,j)) .AND. (k <= lcllev(i,j)) ) THEN
            lzcompute(i,j) = .TRUE.
          ELSE
            lzcompute(i,j) = .FALSE.
          ENDIF
        ENDDO
      ENDDO

      idone = 0

      DO WHILE (idone < idim*jdim)

        idone = 0

        DO j = 1, jdim
          DO i = 1, idim

            computation_IF2: IF (lzcompute(i,j) .AND. (.NOT. ldone(i,j)) ) THEN

              ! Calculate iteratively parcel temperature from 
              ! thp, prs and 1st guess ttt
              tguess1    = ttt(i,j)
              esat       = fesatw(tguess1)
              r1         = rdv*esat/(prs(i,j,k)-esat)
              tguess2    = ttt(i,j) - 1.0_wp
              esat       = fesatw(tguess2)
              r2         = rdv*esat/(prs(i,j,k)-esat)
              thetae1    = fthetae(tguess1,prs(i,j,k),r1)
              thetae2    = fthetae(tguess2,prs(i,j,k),r2)
              ttt(i,j)   = ttt(i,j)+(thetae1-thp(i,j))/(thetae2-thetae1)
              tp(i,j,k)  = ttt(i,j)
              icount(i,j)= icount(i,j) + 1   
               
              ! Stop the calculation if the error of the parcel equivalent 
              ! potential temperature is less than the threshold defined above 
              ! (eps_conv) or in case more than 20 iterations were done
              IF ((ABS(thetae1-thp(i,j)) < eps_conv) .OR. (icount(i,j) >  20)) THEN
                ldone(i,j) = .TRUE.
              ENDIF

            ELSE   ! computation_IF2
              idone      = idone + 1
            ENDIF computation_IF2

          ENDDO
        ENDDO
      ENDDO

      DO j = 1, jdim
!CDIR NOMOVE
        DO i = 1, idim
         
          computation_IF3: IF ( (k <= kstart(i,j)) .AND. (lcomp(i,j)) ) THEN

            IF ( k <= lcllev(i,j) ) THEN                                
              ! update specific humidity of the saturated parcel for new temperature
              esatp  = fesatw(tp(i,j,k))
              qvp(i,j,k) = fqvs(esatp,prs(i,j,k))
            ENDIF
         
            ! Calculate virtual temperatures of parcel and environment
            tvp    = tp(i,j,k)  * (1.0_wp + rvd_m_o*qvp(i,j,k)/(1.0_wp - qvp(i,j,k)) )  
            tve    = te(i,j,k)  * (1.0_wp + rvd_m_o*qve(i,j,k)/(1.0_wp - qve(i,j,k)) ) 
         
            ! Calculate the buoyancy of the parcel at current level k, 
            ! save buoyancy from level k+1 below (buo_belo) to check if LFC or EL have been passed
            buo_belo = buo(i,j)
            buo(i,j) = tvp - tve

            ! Check for level of free convection (LFC) and set flag accordingly. 
            ! Basic LFC condition is that parcel buoyancy changes from negative to 
            ! positive (comparison of buo with buo_belo). Tests showed that very 
            ! often the LFC is already found within the boundary layer below even if 
            ! significant capping inversions are present above (and since CIN is only
            ! defined below the LFC no CIN was accumulated in these cases.)
            ! To handle these situations in a meteorologically meaningful way an 
            ! additional flag "lfcfound" was introduced which is initially zero but 
            ! set to 1 if a second LFC was found, under the condition that the CIN 
            ! within the capping inversion is greater than the CAPE in the convective
            ! boundary layer below times the factor cc_comp (cc_comp = 1 - 2 
            ! recommended.)
            ! Help variable CIN_HELP saves all contributions to the total cin above 
            ! the LFC and has to be subtracted at the end from the final CIN in order
            ! to get the CIN only below the LFC (this is necessary since we do not 
            ! know yet where exactly we will find an LFC when performing the ascent
            ! from bottom to top in a stepwise manner.)

            ! Find the first LFC
            IF ( (buo(i,j) > 0.0_wp) .AND. (buo_belo <= 0.0_wp)            &
                                    .AND. (lfcfound(i,j) == 0_iintegers) ) THEN
            
              ! Check whether it is an LFC at one of the lowest model levels 
              ! (indicated by CAPE=0)
              IF ( (cape(i,j) > 0.0_wp) .AND. (lfcfound(i,j) == 0_iintegers) ) THEN
                ! Check if there is a major capping inversion below, defined as 
                ! having CIN with an absolute value larger than the CAPE accumulated
                ! below times some arbitrary factor cc_comp - if this is the case the
                ! LFC index "lfclev" is updated to the current level k and 
                ! "lfcfound"-flag is now set to 1 assuming that we have found the 
                ! level of free convection finally. 
                IF ( cc_comp * ABS(cin_help(i,j)) > cape(i,j) ) THEN
                  lfclev(i,j)   = k
                  cape  (i,j)   = 0.0_wp
                  cape3km (i,j) = 0.0_wp
                  cin_help(i,j) = 0.0_wp
                  lfcfound(i,j) = 1_iintegers
!                 IF (idebug > 500) WRITE(*,'(A,A18,I2,A,F8.1,A,F8.2,A,I2)') sn,&
!                       'update LFC at k = ',lfclev(i,j),' p = ',prs(i,j,lfclev(i,j))/100.0_wp,&
!                       '  cinhelp =',cin_help(i,j),'  lfcfound =',lfcfound(i,j)
                ENDIF
              ELSE
                ! the LFC found is near the surface, set the LFC index to the current
                ! level k (lfclev) but do not set the flag "lfcfound" to zero to 
                ! indicate that a further LFC may be present above the boundary layer
                ! and an eventual capping inversion. Reset the CIN_HELP to zero to 
                ! store the contribution of CIN above this LFC.
                lfclev(i,j)   = k
                cin_help(i,j) = 0.0_wp
!               IF (idebug > 500) WRITE(*,'(A,A10,I2,A,F8.1,A,F8.2,A,I2)') sn,     &
!                   'LFC at k = ',lfclev(i,j),' p = ',prs(i,j,lfclev(i,j))/100.0_wp,&
!                   '  cinhelp =',cin_help(i,j),'  lfcfound =',lfcfound(i,j)
              ENDIF
            ENDIF
         
            ! Check whether an EL has been passed (bouyancy changing from positive to
            ! negative) under the condition that a LFC has been found below (profile 
            ! is absolutely stable otherwise and EL is not defined.)
            IF ( (buo(i,j) < 0_wp) .AND. (buo_belo >= 0_wp) .AND. (lfclev(i,j) /= 0) ) THEN
              ellev(i,j) = k
!             IF (idebug > 500) WRITE(*,'(A,A10,I2,A,F8.1)') sn,'EL  at k = ', &
!                                        ellev(i,j),   ' p = ',prs(i,j,ellev(i,j))
            ENDIF
         
            ! Accumulation of CAPE and CIN according to definition given in Doswell 
            ! and Rasmussen (1994), 
            IF ( (buo(i,j) >= 0.0_wp) .AND. (k <= lfclev(i,j)) ) THEN   
              cape(i,j)     = cape(i,j)     + (buo(i,j)/tve)*g*(hhl(i,j,k) - hhl(i,j,k+1))
              IF (PRESENT(acape3km)) THEN
                 IF ( k3000(i,j) > 0 .AND. k >= k3000(i,j) ) THEN
                    cape3km(i,j) = cape3km(i,j) + (buo(i,j)/tve)*g*(hhl(i,j,k) - hhl(i,j,k+1))
                 ENDIF
              END IF
            ELSEIF ( (buo(i,j) < 0.0_wp) .AND. (k < kstart(i,j)) ) THEN  
              cin(i,j)      = cin(i,j)      + (buo(i,j)/tve)*g*(hhl(i,j,k) - hhl(i,j,k+1))
              cin_help(i,j) = cin_help(i,j) + (buo(i,j)/tve)*g*(hhl(i,j,k) - hhl(i,j,k+1))
            ENDIF
         
            ! debugging statement
!           IF (idebug > 500) WRITE (*,'(A,A,I2,A,F8.2,A,F8.2,A,F7.1,A,F8.4)')     &
!                sn,' k = ', k, '   cape = ', cape(i,j),'   cin  =  ', cin(i,j), '   p =  ', &
!                    prs(i,j,k)/100.0_wp, '   buo = ', buo(i,j)
       
            ! If 500hPa level approximately reached, save parcel temperature for 
            ! calculation of Showalter Index. Do not check at lowest level of parcel
            ! ascent since pressure at level below is not defined if kstart=kdim 
            ! (parcel starting at lowest model level.)
            IF (present(asi) .AND. (k < kdim)) THEN
              ! Assume 500hPa level is approximately reached if model level pressure
              ! changes from >500hPa to <=500hPa
              IF ( (prs(i,j,k) <= sistopprs) .AND. (prs(i,j,k+1) > sistopprs) ) THEN
                asi(i,j) = te(i,j,k)-tp(i,j,k)
              ENDIF
            ENDIF

          ENDIF computation_IF3
        ENDDO
      ENDDO
    ENDDO  kloop       ! End k-loop over levels
      
    ! Subtract the CIN above the LFC from the total accumulated CIN to 
    ! get only contriubtions from below the LFC as the definition demands.
    DO j = 1, jdim
      DO i = 1, idim
        computation_IF4: IF (lcomp(i,j)) THEN
          ! make CIN positive
          cin(i,j) = ABS (cin(i,j) - cin_help(i,j))
      
          ! set the CIN to missing value if no LFC was found or no CAPE exists
          IF ( (lfclev(i,j) == 0_iintegers) .OR. (cape(i,j) == 0.0_wp)  )  THEN
            cin(i,j) = missing_value 
          ENDIF
        ENDIF computation_IF4
      ENDDO
    ENDDO
      
    ! check which options were given and assign the arguments to be returned
    DO j = 1, jdim
      DO i = 1, idim
        computation_IF5: IF (lcomp(i,j)) THEN
!US       IF (present(tp_print)) tp_print = tp
          IF (present(acape))    acape(i,j)    = cape  (i,j)
          IF (present(acin))     acin (i,j)    = cin   (i,j)
          IF (present(acape3km)) acape3km(i,j) = cape3km (i,j)
          IF (present(ael))      ael  (i,j)    = ellev (i,j)
          IF (present(alfc))     alfc (i,j)    = lfclev(i,j)
          IF (present(alcl))     alcl (i,j)    = lcllev(i,j)
        ENDIF computation_IF5
      ENDDO
    ENDDO
  
  !------------------------------------------------------------
  ! End Subroutine for parcel ascent
  !------------------------------------------------------------

  END SUBROUTINE ascent  

!==============================================================================

END SUBROUTINE cal_conv_ind

!==============================================================================
!==============================================================================

SUBROUTINE calc_bulk_richardson(brn,te,qve,ue,ve,prs,hsurf,prs_surf,t_surf,qv_surf,hhl,&
                                  idim,jdim,kdim,cp_d,r_d,rvd_m_o,g)
  !-------------------------------------------------------------------------------
  ! Description:
  !
  ! The bulk Richardson Number (BRN) is a dimensionless number defined as the bulk 
  ! ratio between the buoyant consumption term and the mechanical production term
  ! in the TKE equation (Stull, 1988).
  !
  ! BRN(z)=g* (z-topo)*(theta_v(z)-theta_v_ground)/(theta_v_mean * (u**2+v**2))
  ! 
  ! theta_v_mean is the virtual potential temperature in the layer comprised
  ! between the ground and the height z.
  !
  ! We here consider u=0 m/s and v=0 m/s at the ground (some authors consider
  ! u and v at 2meters). Some others consider theta_v_ground instead of 
  ! theta_v_mean.
  !
  ! References:
  ! ----------- 
  ! - Stull R.,B., 1988: An Introduction to Boundary Layer Meteorology. 
  !   Kluwer Academic Publishers. ISBN 90-277-2768-6
  !
  ! - Sorensen J.H., Rasmussen A., Svensmark H.,1996: Forecast of Atmospheric Boundary-Layer 
  !   Height Utilised for ETEX Real-Time Dispersion Modelling. 
  !   Phys. Chem. Earth, Vol. 21, No. 5-6, pp. 435-439, 1997 Elsevier Science Ltd.
  !
  ! - Jacobson M.Z., 1999: Fundamentals of atmospheric modeling.
  !   Cambridge University Press. ISBN 0-521-63143-2
  !-------------------------------------------------------------------------------

  ! Input data
  !----------- 
  INTEGER (KIND=iintegers), INTENT (IN) ::  &
       idim ,                  & ! array dimension in zonal direction
       jdim ,                  & ! array dimension meridional direction
       kdim                      ! array dimension in vertical direction
  
  REAL    (KIND=wp),        INTENT (IN) ::  &
       te  (idim,jdim,kdim),   & ! environment temperature
       qve (idim,jdim,kdim),   & ! environment humidity
       ue  (idim,jdim,kdim),   & ! environment zonal wind 
       ve  (idim,jdim,kdim),   & ! environment meridional wind
       prs (idim,jdim,kdim),   & ! full level pressure
       hsurf(idim,jdim),       & ! topography
       prs_surf(idim,jdim) ,   & ! surface pressure
       t_surf(idim,jdim)   ,   & ! surface temperature
       qv_surf(idim,jdim)  ,   & ! specific water vapor content on the surface   (kg/kg)
                                 ! in fact the 2m values are used here (on purpose!)
       hhl (idim,jdim,kdim+1)    ! height of half levels
  
  
  ! Physical constants
  !-------------------
  REAL    (KIND=wp),        INTENT (IN) ::  &
       cp_d,                   & ! specific heat of dry air at constant pressure
       r_d,                    & ! gas constant for dry air
       rvd_m_o,                & ! r_v/r_d - 1
       g                         ! acceleration due to gravity
  
  
  ! Output data
  !------------ 
  REAL (KIND=wp),     INTENT (OUT)    :: &
       brn(idim,jdim,kdim)      ! Bulk Richardson Number [-]
  
  
  ! Local scalars and automatic arrays
  !-----------------------------------
  REAL    (KIND=wp)     ::       &      
       u_avg(idim,jdim,kdim),   & ! averaged zonal wind at the center of the cells
       v_avg(idim,jdim,kdim),   & ! averaged meridional wind at the center of the cells
       hfl(idim,jdim,kdim),     & ! height of full levels
       theta_v(kdim),           & ! virtual potential temperature on each level
       pasl,                    & ! sea level pressure (ref)
       theta_s,                 & ! potential temperature at the surface 
       theta_v_s,               & ! virtual potential temperature at the surface 
       theta,                   & ! potential temperature on each level
       cumul,                   & ! cumul of virtual potential temperature
       theta_v_cum,             & ! cumulated virtual potential temperature
       theta_v_avg                ! mean virtual potential temperature between k and kdim
  
  INTEGER    (KIND=iintegers) :: &     
       i, j, k                   ! Indices of input/output fields
  
  !-----------------------------------------------------------------------------
  
  ! Reference pressure for potential temperature calculation
  pasl             = 100000._wp    
  
  ! Compute height on full LM COARSE levels
  DO k = 1, kdim
     hfl(:,:,k) = 0.5_wp * hhl(:,:,k+1) + 0.5_wp * hhl(:,:,k)
  ENDDO
  
  
  !-------------------------------------------------------------------------
  ! Compute averaged wind velocities at the center of the cells
  !-------------------------------------------------------------------------
  ! zonal wind
  
  ! NOTE: staggering has been changed in Version 4.12
  !       the old code did not de-stagger at the correct position
  DO k=1, kdim
    DO j=1, jdim
      DO i=1, idim
        u_avg(i,j,k) = 0.5_wp*(ue(MAX(1,i-1),j,k) + ue(i,j,k))
        v_avg(i,j,k) = 0.5_wp*(ve(i,MAX(1,j-1),k) + ve(i,j,k))
      ENDDO
    ENDDO
  ENDDO
 
  !Grid points loop
  
  DO j=1, jdim   !y loop
     DO i=1, idim  !x loop
        
        !-----------------------------------------------------------------------
        ! Calculation of all components of the bulk richardson number and 
        ! of the bulk Richardson number on each full level
        !-----------------------------------------------------------------------
        
        ! Potential temperature at the surface 
        theta_s = t_surf(i,j) * (pasl / prs_surf(i,j))**(r_d/cp_d)
        
        ! Virtual potential temperature at the surface           
        theta_v_s = theta_s * (1._wp + rvd_m_o * qv_surf(i,j))
        
        ! Initial cumulated virtual potential temperature
        cumul = 0.0_wp
        
        DO k = kdim, 1, -1  ! Level loop 
           
           ! Virtual potential temperature at each level
           theta = te(i,j,k) * (pasl / prs(i,j,k))**(r_d/cp_d)
           theta_v (k)= theta * (1 + rvd_m_o * qve(i,j,k)) 
           
           ! Mean virtual potential temperature between k and kdim
           theta_v_cum  = cumul + theta_v(k)
           theta_v_avg  = theta_v_cum / (kdim - k + 1)
           
           !Bulk Richardson number
           !**********************
           
           ! Limitation of the wind: 1cm/s minimum for each component
           
           IF((u_avg(i,j,k)**2 + v_avg(i,j,k)**2) .LE. 2.0E-4_wp)THEN
              brn(i,j,k) = g *(hfl(i,j,k)-hsurf(i,j))*(theta_v(k)-theta_v_s)/ &
                   (theta_v_avg*2.0E-4_wp)
              
           ELSE
              brn(i,j,k) = g *(hfl(i,j,k)-hsurf(i,j))*(theta_v(k)-theta_v_s)/ &
                   (theta_v_avg*(u_avg(i,j,k)**2+v_avg(i,j,k)**2))
           ENDIF
           
           
           !Update cumulated virtual potential temperature
           cumul = theta_v_cum
        ENDDO ! Level loop 
     ENDDO ! x
  ENDDO ! y
  
END  SUBROUTINE calc_bulk_richardson
!===============================================================================

SUBROUTINE calc_pbl_brn(te, qve, prs, hhl, hsurf,brn, idim, jdim, kdim, cp_d, r_d, &
                          rvd_m_o, missing_value, hpbl)
  !-------------------------------------------------------------------------------
  ! Description:
  !
  ! The PBL height is the estimated planetary boundary layer height in meters 
  ! above the ground. There are several methods to estimate it. In general they
  ! perform well for convective situations. The stable situations are more problematic.
  ! 
  ! We choose the bulk richardson method, which is a widely used and robust method. 
  ! The height of the PBL is given by the height at which the bulk Richardson
  ! number reaches a prescribed value, the critical Richardson number.
  ! In the literature, one finds critical values between say 0.2 and 1.
  ! We here consider a value of 0.33 under stable conditions (Wetzel, 1982)and
  ! of 0.22 under convective conditions (Vogelezang and Holtslag, 1996).
  !
  ! References:
  ! -----------
  ! - Stull R.,B., 1988: An Introduction to Boundary Layer Meteorology.
  !   Kluwer Academic Publishers. ISBN 90-277-2768-6
  !
  ! - Sorensen J.H., Rasmussen A., Svensmark H.,1996: Forecast of Atmospheric Boundary-Layer 
  !   Height Utilised for ETEX Real-Time Dispersion Modelling. 
  !   Phys. Chem. Earth, Vol. 21, No. 5-6, pp. 435-439, 1997 Elsevier Science Ltd.
  !
  ! - Seibert P. and all, 1999: Review and intercomparison of operational methods for the 
  !   determination of the mixing height.
  !   Atmospheric Environment 34 (2000) 1001-1027, 2000 Elsevier Science Ltd.  
  !
  ! - Vogelezang D.H.P., and Holtslag A.A.M., 1996: Evaluation and model impacts
  !   of alternative boundary-layer height formulations.
  !   Boundary-Layer Meteorol. 81, 245-269.
  !
  ! - Wetzel P.J., 1982: Toward parametrization of the stable boundary layer.
  !   J. Appl. Meteorol. 21, 7-13.
  !-------------------------------------------------------------------------------

  ! Input data
  !----------- 
  INTEGER (KIND=iintegers), INTENT (IN) ::  &
       idim ,                  & ! array dimension in zonal direction
       jdim ,                  & ! array dimension meridional direction
       kdim                      ! array dimension in vertical direction
  
  REAL    (KIND=wp),        INTENT (IN) ::  &
       te  (idim,jdim,kdim),   & ! environment temperature
       qve (idim,jdim,kdim),   & ! environment humidity
       prs (idim,jdim,kdim),   & ! full level pressure
       hhl (idim,jdim,kdim+1), & ! height of half levels
       hsurf(idim,jdim),       & ! topo
       brn (idim,jdim,kdim)      ! Bulk Richardson Number [-]
  
  ! Physical constants
  !-------------------
  REAL    (KIND=wp),        INTENT (IN) ::  &
       cp_d,                   & ! specific heat of dry air at constant pressure
       r_d,                    & ! gas constant for dry air
       rvd_m_o,                & ! r_v/r_d - 1
       missing_value             ! Missing value for CIN (if no LFC/CAPE was found),
  
  
  ! Output data
  !------------ 
  REAL (KIND=wp),     INTENT (OUT)    :: &
       hpbl(idim,jdim)            ! PBL height [m above ground]
  
  
  ! Local scalars and automatic arrays
  !-----------------------------------
  REAL    (KIND=wp)                   ::  &      
       hfl(idim,jdim,kdim),                 & ! height of full levels
       zh(idim,jdim),zh2(idim,jdim),        & ! working arrays for the linear regression
       ! zh=sum(hfl); zh2=sum(hfl**2)
       zttav(idim,jdim),zhttav(idim,jdim),  & ! working arrays for the linear regression
       ! zttav=sum(theta_v);zhttav=sum(hfl*theta_v)
       beta(idim,jdim),                     & ! coefficient of linear regression
       theta_v(kdim),                       & ! virtual potential temperature on 4 levels
       theta,                               & ! potential temperature on each level
       pasl,                                & ! sea level pressure (ref) 
       brn_cr                                 ! bulk richardson critical value
  
  
  INTEGER    (KIND=iintegers)       ::   &
       i,j,k,                               & ! working indexes
       k_pbl(idim,jdim)                       ! k-index corresponding to the top of the PBL
  
  !------------------------------------------------------------------------
  ! Initialisations
  !------------------------------------------------------------------------
  ! Sums for linear regression
  zh(:,:) = 0.0_wp
  zttav(:,:) = 0.0_wp
  zh2(:,:)= 0.0_wp
  zhttav(:,:)= 0.0_wp
  ! Index of the top of the PBL
  k_pbl(:,:)= -1
  ! ref pressure
  pasl = 100000.0_wp
  
  ! Compute height on full LM COARSE levels
  DO k = 1, kdim
     hfl(:,:,k) = 0.5_wp * hhl(:,:,k+1) + 0.5_wp * hhl(:,:,k)
  ENDDO
  
  
  ! Determine if stable or unstable (neutral) conditions
  !------------------------------------------------------------------------
  ! calculation of the linear regression coefficient in order to approximate
  ! the temperature gradient in the first 4 layers and evaluate the current
  ! conditions.
  ! if delta_theta_v/delta_z > 0: stable (beta >0)
  ! if delta_theta_v/delta_z < 0: unstable (or neutral) (beta <=0)
  
  DO j=1,jdim
     DO i=1,idim
        
        DO k=kdim-3,kdim
           theta = te(i,j,k) * (pasl / prs(i,j,k))**(r_d/cp_d) ! potential temp.
           theta_v (k)= theta * (1 + rvd_m_o * qve(i,j,k))     ! virtual potent. temp.
           
           zh (i,j) = zh (i,j) + hfl(i,j,k)                   ! sum of the heights        
           zh2(i,j) = zh2(i,j) + hfl(i,j,k)* hfl(i,j,k)       ! sum of the (heights)**2   
           zttav (i,j) = zttav (i,j) + theta_v(k)             ! sum of the theta_v
           zhttav(i,j) = zhttav(i,j) + hfl(i,j,k)* theta_v(k) ! sum of the (height*theta_v)
        ENDDO
        
        ! coefficient of linear regression
        beta (i,j)=(zhttav(i,j)-zh(i,j)*zttav(i,j)/4.0_wp)/    &
             (zh2(i,j)-zh(i,j)*zh(i,j)/4.0_wp)  
        
        ! set the critical bulk richardson number according to the current
        ! conditions 
        
        IF(beta (i,j) .GT. 0.0_wp) THEN     ! Stable conditions 
           brn_cr = 0.33_wp                     ! Wetzel value    
        ELSE                                    ! Unstable or neutral conditions
           brn_cr = 0.22_wp                     ! Vogelezang value   
        ENDIF
        
        !--------------------------------------------------------------------------
        ! Find k index corresponding to the upper limit of the BL and compute
        ! PBL height
        ! -------------------------------------------------------------------------
        
        ! Scan from the bottom upwards and find
        ! the first level where BRN > brn_cr
        
        DO k= kdim, 1, -1
           IF(brn(i,j,k) > brn_cr) THEN
              k_pbl(i,j) = k
              EXIT
           ENDIF
        ENDDO
        
        ! Compute the corresponding height. If no level found
        ! where BRN < brn_cr, assign a missing_value for the PBL
        ! height.
        
        IF (k_pbl(i,j).LE.0) THEN                  ! no level found
           hpbl(i,j) = missing_value
        ELSE                                       ! algorithm ok
           hpbl(i,j) =  hfl(i,j,k_pbl(i,j))-hsurf(i,j)
        ENDIF
        
     ENDDO !x loop
  ENDDO    !y loop
  
END  SUBROUTINE calc_pbl_brn

!==============================================================================
!==============================================================================

SUBROUTINE calc_ceiling( z_cld_ceiling_hgt, clc_sgs, hhl, ie, je, ke )

!------------------------------------------------------------------------------
!
! Description:
!   Ceiling = height above MSL, for which cloud coverage > 4/8
! Method:
!
!------------------------------------------------------------------------------

! Parameterlist
INTEGER (KIND=iintegers), INTENT(IN) :: &
  ie, je, ke                       ! field dimensions

REAL (KIND=wp),     INTENT(IN)       :: &
  clc_sgs (ie,je,ke),            & ! subgrid scale cloud cover
  hhl     (ie,je,ke+1)             ! height of half levels

REAL (KIND=wp),     INTENT(OUT)      :: &
  z_cld_ceiling_hgt(1:ie,1:je)     ! Cloud ceiling height above MSL (in m)

! Local variables
LOGICAL                 :: z_cld_base_found(1:ie, 1:je)
INTEGER(KIND=iintegers) :: i, j, k

!------------------------------------------------------------------------------

  z_cld_base_found(:,:)  = .FALSE.
  z_cld_ceiling_hgt(:,:) = hhl(:,:,1)  ! arbitrary default value

  DO k = ke, 1, -1
    DO j = 1, je
      DO i = 1, ie

        IF ( .NOT.(z_cld_base_found(i,j)) .AND.   &
             ( clc_sgs(i,j,k) > 0.5_wp ) )     THEN
          z_cld_ceiling_hgt(i,j) = 0.5_wp * ( hhl(i,j,k) + hhl(i,j,k+1) )
          z_cld_base_found(i,j) = .TRUE.
        ENDIF

      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE calc_ceiling

!==============================================================================
!==============================================================================

SUBROUTINE calpmsl ( pmsl, ps, t, rho0, dp0, hsurf, ie, je, g, r_d )

!------------------------------------------------------------------------------
!
! Description:
!   This routine computes the mean sea-level pressure by extrapolating the
!   surface pressure at height z = hsurf hydrostatically to z = 0 whereever
!   the surface height is greater than 0.01 m. For hsurf < 0.01 m, the
!   mean sea-level pressure is simply set to ps.
!   
!   The fields are passed twodimensional, that means the calling procedure has
!   to choose the proper level and time level.
!
! Method:
!   Formula to be described in the documentation
!
!------------------------------------------------------------------------------

! Parameter list:

INTEGER (KIND=iintegers), INTENT (IN)    ::    &
  ie, je              ! horizontal dimensions of the fields

REAL (KIND=wp),     INTENT (IN)          ::    &
  ps   (ie,je)    , & ! suface pressure      
  t    (ie,je)    , & ! temperature at lowest model level
  dp0  (ie,je)    , & ! full level pressure thickness of lowest layer
  rho0 (ie,je)    , & ! reference density at lowest full model level
  hsurf(ie,je)        ! geometrical height of surface

REAL (KIND=wp),     INTENT (OUT)         ::    &
  pmsl (ie,je)        ! mean sea-level pressure

REAL (KIND=wp),     INTENT (IN)          ::    &
  g,                & ! gravity acceleration
  r_d                 ! gas constant for dry air

! Local Variables

REAL (KIND=wp)        ::    &
  zlapse,           & ! lapse rate ( 6.5K / km )
  zrg,              & ! 1/g                      
  zr3,              & ! 1/3                      
  zdz,              & ! height difference: lowest level - surface 
  ztstar,           & ! temperature at surface (extrapolated from lowest level)
  ztmsl,            & ! temperature at mean sea level (extrapolated from tstar)
  zalph,            & ! zlapse*r_d/g
  zprt,             & ! 
  zprtal              !

INTEGER (KIND=iintegers)    ::    &
  i,j                 ! Loop indicees

!------------------------------------------------------------------------------

! Begin subroutine calpmsl

! Initializations
  zlapse = 0.0065_wp
  zrg    = 1.0_wp/g
  zr3    = 1.0_wp/3.0_wp

! Set mean sea-level pressure to ps
  pmsl (:,:) = ps (:,:)

  DO j = 1, je
    DO i = 1, ie

    ! Calculate surface temperature "tstar" and mean sea-level 
    ! temperature "tmsl" and correct lapse rate "zalph"
      zdz    = 0.5_wp*dp0(i,j)/ ( g*rho0(i,j) )
      ztstar = t(i,j) + zlapse*zdz
      IF( ztstar < 255.0_wp ) THEN
          ztstar = 0.5_wp*( 255.0_wp + ztstar )
      ENDIF
      ztmsl = ztstar + zlapse*hsurf(i,j)
      zalph = r_d*zlapse*zrg
      IF( ztstar > 290.5_wp ) THEN
          ztstar = 0.5_wp*( 290.5_wp + ztstar )
          ztmsl  = ztstar
          zalph  = 0.0_wp
      ENDIF
      IF( ztstar <= 290.5_wp .AND. ztmsl > 290.5_wp) THEN
          ztmsl  = 290.5_wp 
          zalph  = r_d*(ztmsl-ztstar)*zrg/MAX(hsurf(i,j), 0.01_wp)
      ENDIF

      ! Calculate mean sea level pressure for elevated terrain
      zprt   = g*hsurf(i,j)/(r_d*ztstar)
      zprtal = zprt*zalph 
      pmsl(i,j) = ps(i,j)*EXP &
                 ( zprt*(1.0_wp - zprtal*(0.5_wp - zr3*zprtal)) )
    ENDDO
  ENDDO

END SUBROUTINE calpmsl

!==============================================================================
!==============================================================================

SUBROUTINE calrelhum ( relhum, t, pp, p0, qv, ie, je, ke, &
                       b1, b2w, b3, b4w, rdv, o_m_rdv ) 

!------------------------------------------------------------------------------
!
! Description:
!   This routine computes the relative humidity (over water) in the domain
!   (ie,je,ke)
!   
!   The fields are passed threedimensional, that means the calling procedure 
!   has to choose the proper time level.
!
! Method:
!   Formula described in the documentation (Part I, Sec. 7.2).
!
!------------------------------------------------------------------------------

! Parameter list:

INTEGER (KIND=iintegers), INTENT (IN)    ::    &
  ie, je, ke          ! dimensions of the fields

REAL (KIND=wp),     INTENT (IN)          ::    &
  t    (ie,je,ke) , & ! temperature at full levels
  pp   (ie,je,ke) , & ! perturbation pressure at full levels
  p0   (ie,je,ke) , & ! reference pressure at full levels  
  qv   (ie,je,ke)     ! specific humidity

REAL (KIND=wp),     INTENT (OUT)       ::    &
  relhum (ie,je,ke)   ! relative humidity (resp. water saturation)

REAL (KIND=wp),     INTENT (IN)          ::    &
  b1, b2w, b3, b4w,  & ! constants to calculate the vapour pressure
  rdv,               & ! r_d/r_v
  o_m_rdv              ! 1 - r_d/r_v

!local variables
INTEGER (KIND=iintegers) ::  i, j, k    
REAL    (KIND=wp)        ::  zpvs, zp, zqvs
!------------------------------------------------------------------------------

DO k = 1, ke
  DO j = 1, je
    DO i = 1, ie
      zpvs = b1*EXP( b2w*(t(i,j,k)-b3) / (t(i,j,k)-b4w) )
      zp   = p0(i,j,k) + pp(i,j,k)
      zqvs = rdv*zpvs / (zp - o_m_rdv*zpvs)
      ! Set minimum value of relhum to 0.01 (was 0 before)
      relhum(i,j,k) = MAX( 0.01_wp, qv(i,j,k)/zqvs * 100_wp )
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE calrelhum

!==============================================================================
!==============================================================================

SUBROUTINE calomega ( omega, pp2, pp1, pptens, w, rho0, ie, je, ke, dt, g )

!------------------------------------------------------------------------------
!
! Description:
!   This routine computes the vertical velocity p-dot in a pressure coordinate
!   system
!   
!   The fields are passed threedimensional, that means the calling procedure 
!   has to choose the proper time level.
!
! Method:
!   Evaluation of p-dot
!
!------------------------------------------------------------------------------

! Parameter list:

INTEGER (KIND=iintegers), INTENT (IN)    ::    &
  ie, je, ke          ! dimensions of the fields

REAL (KIND=wp),     INTENT (IN)          ::    &
  pp2   (ie,je,ke) , & ! perturbation pressure at full levels for time n+1
  pp1   (ie,je,ke) , & ! perturbation pressure at full levels for time n  
  pptens(ie,je,ke),  & ! advective perturbation pressure tendency
  rho0  (ie,je,ke) , & ! reference density at full levels  
  w     (ie,je,ke+1)   ! vertical velocity

REAL (KIND=wp),     INTENT (OUT)       ::    &
  omega  (ie,je,ke)   ! p-dot

REAL (KIND=wp),     INTENT (IN)          ::    &
  dt,               & ! long time step      
  g                   ! gravity acceleration

!local variables
INTEGER (KIND=iintegers) ::  i, j, k    
REAL    (KIND=wp)        ::  zdtr
!------------------------------------------------------------------------------

zdtr = 1.0_wp/dt

DO k = 1, ke
  DO j = 1, je
    DO i = 1, ie
      omega(i,j,k) = ( pp2(i,j,k) - pp1(i,j,k) )*zdtr - pptens(i,j,k) &
                      -0.5_wp*g*rho0(i,j,k)*( w(i,j,k) + w(i,j,k+1) )
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE calomega

!==============================================================================
!==============================================================================

SUBROUTINE calprsum ( tot_prec, rrsn, rssn, rrkn, rskn, ie, je )

!------------------------------------------------------------------------------
!
! Description:
!   This routine computes the total amount of precipitation
!
! Method:
!   Summation of the precipitation amouts due to grid-scale as well as
!   convective rain and snow.
!
!------------------------------------------------------------------------------

! Parameter list:

INTEGER (KIND=iintegers), INTENT (IN)    ::    &
  ie, je              ! dimensions of the fields

REAL (KIND=wp),     INTENT (OUT)         ::    &
  tot_prec(ie,je)     ! total amount of precipitation

REAL (KIND=wp),     INTENT (IN)          ::    &
  rrsn (ie,je),     & !
  rssn (ie,je),     & !
  rrkn (ie,je),     & !
  rskn (ie,je)        !

!------------------------------------------------------------------------------

! Begin subroutine calprsum


  tot_prec(1:ie,1:je)  =   rrsn(1:ie,1:je) + rssn(1:ie,1:je) &
                         + rrkn(1:ie,1:je) + rskn(1:ie,1:je)

END SUBROUTINE calprsum

!==============================================================================
!==============================================================================

SUBROUTINE calsnowlmt ( snowlmt, t, pp, p0, qv, hhl, hhlr, ie, je, ke, t0, wbl)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine calculates height of the snowfall limit (snowlmt).
!
! Method:
!   In a first step the wet bulb temperature is derived from pp, t and qv.
!   In a second step the snowfall limit is evaluated from 8000m down to the
!   the lowest model level (ke) and linearly interpolated to the height where
!   the wet bulb temperature is >= wbl (=+1.3C after P. Haechler, MeteoSwiss).
!   A flag (-999) is set to indicate that no snowlmt was found.
!
!------------------------------------------------------------------------------
! Parameter list:

INTEGER (KIND=iintegers), INTENT (IN)    ::    &
  ie, je, ke        ! array dimensions of the I/O-fields

REAL (KIND=wp),     INTENT (OUT)   ::  &
  snowlmt(ie,je)    ! height of the snowfall limit in m above sea level

REAL (KIND=wp),     INTENT (IN)    ::  &
  t  (ie,je,ke) , & ! temperature
  pp (ie,je,ke) , & ! perturbation pressure at full levels
  p0 (ie,je,ke) , & ! reference pressure at full levels
  qv (ie,je,ke) , & ! specific humidity
  hhl(ie,je,ke+1),& ! height of model half levels
  hhlr(ke+1)        ! height of model half levels resp. sea level

REAL (KIND=wp),     INTENT (IN)    ::  &
  t0,             & ! freezing point temperature
  wbl               ! (empirical) wet bulb temperature at snowfall limit (1.3C)
!------------------------------------------------------------------------------
! Local variables

INTEGER (KIND=iintegers) ::    &
  i,j,k,                       & ! Loop indices
  ktopmin

LOGICAL                  ::    &
  lfound(ie,je)     ! Logical flag : =.TRUE when wet bulb temp corresponding to
                    !                  parameter "wbl" is found

REAL (KIND=wp)           ::    &
  za = 0.78588481_wp,      & ! local storage
  zb = 7.567_wp,           &
  zc = 2066.92605_wp,      &
  zd = 33.45_wp,           &
  ze = 0.622_wp,           &
  zf = 0.378_wp,           &
  zg = 0.5_wp,             &
  zh = 0.6_wp,             &
  zi = 700._wp,            &
  zl = 0.1_wp,             &
  zm = 6400._wp,           &
  zn = 11.564_wp,          &
  zo = 1742._wp,           &
  td,tl,tp,ppp,                &
  deltat,zt,zp,                &
  ep,const,                    &
  zh_bot, zh_top,              &
  zdt

REAL (KIND=wp)           ::    &
  wetblb(ie,je,ke)

!------------------------------------------------------------------------------

! Begin subroutine calsnowlmt

! Set the uppermost model level for the occurence of a wet bulb temperature (wbl)
! to about 8000m above surface
ktopmin = 2
DO k = ke+1, 1, -1
  IF ( hhlr(k) < 8000.0_wp ) THEN
    ktopmin = k
  ENDIF
ENDDO

! Initialize the definition mask and the output array snowlmt
lfound (:,:) = .FALSE.
snowlmt(:,:) = -999.0_wp

DO k = ktopmin, ke
  DO j = 1, je
    DO i = 1, ie
      zp     = (p0(i,j,k) + pp(i,j,k))/100._wp     ! in hPa
      ep     = MAX(1.0E-10_wp,qv(i,j,k))*zp /      &
               (ze + zf*MAX(1.0E-10_wp,qv(i,j,k)))
      ep     = MAX(ep,1.0E-10_wp)
      CONST  = LOG10(ep) - za
      td     = (zd*CONST-zc) / (CONST-zb)              ! in Kelvin
      ! Wet bulb temperature after Egger/Joss
      tl     = (t(i,j,k) - t0) *10._wp
      tp     = (td-t0) *10._wp
      ppp    = zp * 10._wp
      deltat = tl-tp
      zt     = tp + zg*deltat*(zh-tp/zi)
      wetblb(i,j,k) = zl * ( tp +                      & ! in Celsius
                      (deltat / (1._wp + zm*EXP(zn*zt/(zo+zt))/ppp)))

      IF ( wetblb(i,j,k) >= wbl ) THEN
      ! definition of snowlmt can be made in this column
        lfound (i,j) = .TRUE.
      ENDIF
    ENDDO
  ENDDO
ENDDO

DO k = ktopmin+1, ke
  DO j = 1, je
    DO i = 1, ie
      IF ( lfound(i,j) .AND. wetblb(i,j,k) >= wbl ) THEN
      ! definition of snowlmt is now made once
        lfound (i,j) = .FALSE.
        zh_bot       = 0.5_wp * ( hhl(i,j,k) + hhl(i,j,k+1) )
        zh_top       = 0.5_wp * ( hhl(i,j,k) + hhl(i,j,k-1) )
        zdt          = ( wbl - wetblb(i,j,k) ) /                   &
                       ( wetblb(i,j,k-1) - wetblb(i,j,k) )
        snowlmt(i,j) = zh_bot + (zh_top-zh_bot)*zdt
      ENDIF
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE calsnowlmt

!==============================================================================
!==============================================================================

SUBROUTINE caltopdc (top_dc, t, p0, pp, qv, hhl, hhlr, ie, je, ke,     &
                    b1, b2w, b3, b4w, rdrd, emrdrd, rddrm1, cpdr, g )  

!------------------------------------------------------------------------------
! Description:
!   This subroutine calculates the top layer-index of dry convection with
!   roots at the ground.
!
! Method:
!   An air parcel with environmental properties (plus a small excess
!   temperature) is lifted dry-adiabatically from the lowest model level.
!   The top height of dry convection is then defined to be the (full) 
!   level where the parcel either looses buoyancy or becomes saturated.
!   During the parcel's ascent, an instantaneous pressure adjustment is
!   assumed.
!
!------------------------------------------------------------------------------

! Parameter list:

INTEGER (KIND=iintegers), INTENT (IN)    ::    &
  ie, je, ke        ! dimensions of the fields

REAL (KIND=wp),     INTENT (OUT)         ::    &
  top_dc (ie,je)    ! top of dry convection (Height in m above sea level)

REAL (KIND=wp),     INTENT (IN)          ::    &
  t  (ie,je,ke),  & ! temperature 
  p0 (ie,je,ke),  & ! base-state pressure
  pp (ie,je,ke),  & ! pressure perturbation
  qv (ie,je,ke),  & ! specific humidity    
  hhl(ie,je,ke+1),& ! height of model half levels
  hhlr(ke+1)        ! height of model half levels resp. sea level

REAL    (KIND=wp),        INTENT (IN)    ::  &
  b1, b2w, b3, b4w, rdrd, emrdrd, rddrm1, cpdr, g  

! Local variables

INTEGER (KIND=iintegers) ::    &
  i,j,k,          & ! Loop indicees
  ktopmin,        & ! Minimum top level index of dry convection
  ilab  (ie,je)     ! Index for of buoyancy at a level

REAL (KIND=wp)           ::    &
  fgew, fgqd,     & ! Name of statement functions
  zt,zp,zge ,     & ! Dummy variables (for temperature, pressure and vapour
                    ! pressure in statement functions fgew (for saturation
                    ! vapour pressure) and fgqd (for sat. spec. humidity)
  ztve      ,     & ! Environmental virtual temperature
  zpe       ,     & ! Environmental (and parcel) pressure
  zdz       ,     & ! Height difference of model full levels (>0)
  zbuoy     ,     & ! Parcel buoyancy
  zqsat     ,     & ! Saturation specific humidity of parcel
  zcond     ,     & ! Condensation within parcel
  ztp(ie,je),     & ! Parcel temperature
  zqp(ie,je)        ! Parcel specific humidity

!------------------------------------------------------------------------------
! STATEMENT FUNCTIONS
  fgew(zt)       = b1*EXP( b2w*(zt-b3)/(zt-b4w) )
  fgqd(zge,zp)   = rdrd*zge/( zp - emrdrd*zge )

! Begin subroutine caltopdc
 
! Initialize the output array with zeros
  top_dc(:,:) = 0.0_wp

! Set the maximum top index of dry convection to about 3000 m above surface 
  ktopmin = 0
  DO k = ke+1, 1, -1
    IF ( hhlr(k) < 3000.0_wp ) THEN
      ktopmin = k
    ENDIF
  ENDDO

! Set the initial values for the rising thermal (starting at level ke)
  DO j = 1, je
    DO i = 1, ie
      ilab  (i,j) = 1 
      ztp   (i,j) =  t(i,j,ke) + 0.2_wp
      zqp   (i,j) = qv(i,j,ke)
    ENDDO
  ENDDO

! Lift the parcel adiabatically and check for buoyancy.
! IF the parcel is buoyant and the lifting condensation level
! is not reached, the ascent is continued.
  DO k = ke-1, ktopmin, -1
    DO j = 1, je
      DO i = 1, ie
        ztve = t(i,j,k)*( 1.0_wp + rddrm1*qv(i,j,k) )
        zpe  = p0(i,j,k) + pp(i,j,k)
        zdz  = 0.5_wp*( hhl(i,j,k) - hhl(i,j,k+2) )
        ztp(i,j) = ztp(i,j) - cpdr*g*zdz
        IF ( ilab(i,j) == 1 ) THEN
          zbuoy  = ztp(i,j)*( 1.0_wp + rddrm1*zqp(i,j) ) - ztve
          zqsat  = fgqd( fgew(ztp(i,j)), zpe )
          zcond  = zqp(i,j) - zqsat
          IF ( zcond < 0.0_wp .AND. zbuoy > 0.0_wp) THEN
             top_dc(i,j) =  hhl(i,j,k)
          ELSE
             ilab(i,j) = 0
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE caltopdc

!==============================================================================
!==============================================================================

SUBROUTINE calhzero ( hzerocl, t, hhl, hhlr, ie, je, ke, t0 )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine calculates the height of the 0 Celsius isotherm (hzerocl)
!   from the height of 8km downards to where t >= t0.
!
! Method:
!   The isotherm is evaluated from 8000m down to the lowest model level (ke)
!   and is set to the height where the melting temperature t0 (0 Celsius) is
!   reached for the first time.
!   In the case where an to isotherm exists, a linear interpolation of the
!   height is done where the temperature changes from t<t0 to t>=to in
!   the layer below.
!   If no t0 isotherm was found (either since the t0 isotherm is higher than 8km
!   or below model topography) hzerocl is set to missing value.
!
!------------------------------------------------------------------------------
! Parameter list:

INTEGER (KIND=iintegers), INTENT (IN)    ::    &
  ie, je, ke        ! array dimensions of the I/O-fields

REAL (KIND=wp),     INTENT (OUT)   ::  &
  hzerocl(ie,je)    ! height of the 0 Celsius isotherm in m above sea level

REAL (KIND=wp),     INTENT (IN)    ::  &
  t  (ie,je,ke),  & ! temperature
  hhl(ie,je,ke+1),& ! height of model half levels
  hhlr(ke+1)        ! height of model half levels resp. sea level

REAL (KIND=wp),     INTENT (IN)    ::  &
  t0                ! freezing point temperature
!------------------------------------------------------------------------------
! Local variables

INTEGER (KIND=iintegers) ::    &
  i,j,k,          & ! Loop indicees
  ktopmin           ! Minimum top level index of dry convection

LOGICAL                  ::    &
  lfound(ie,je)     ! Logical flag : =.TRUE when 0 Celsius level is found

REAL (KIND=wp)           ::    &
  zh_bot, zh_top, zdt  ! local storage

!------------------------------------------------------------------------------

! Begin subroutine calhzero

! Set the uppermost model level for the occurence of a 0 Celius temperature to
! about 8000 m above surface
  ktopmin = 2
  DO k = ke+1, 1, -1
    IF ( hhlr(k) < 8000.0_wp ) THEN
      ktopmin = k
    ENDIF
  ENDDO

! Initialize the output array hzerocl with a flag (-999) to indicate
! an undefined freezing height
  hzerocl(:,:) = -999.0_wp
  lfound (:,:) = .FALSE.

  DO k = ktopmin, ke
    DO j = 1, je
      DO i = 1, ie
        IF ( .NOT. lfound(i,j) .AND. t(i,j,k) >= t0 ) THEN
          ! definition of hzerocl is now made once
          lfound(i,j)  = .TRUE.
          zh_bot       = 0.5_wp * ( hhl(i,j,k) + hhl(i,j,k+1) )
          zh_top       = 0.5_wp * ( hhl(i,j,k) + hhl(i,j,k-1) )
          zdt          = ( t0 - t(i,j,k) ) / ( t(i,j,k-1) - t(i,j,k) )
          hzerocl(i,j) = zh_bot + (zh_top-zh_bot)*zdt
        ENDIF
      ENDDO
    ENDDO
  ENDDO

! This would set the 0-degree isotherm to about 8000 meters, which is 
! very unlikely
! WHERE ( .NOT. lfound(:,:) )
!   hzerocl(:,:) = hhl(:,:,ktopmin)
! ENDWHERE

END SUBROUTINE calhzero

!==============================================================================
!==============================================================================

SUBROUTINE calcldepth ( cldepth, clc_sgs, clc_con, dp0, ie, je, ke )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine calculates the normalized cloud depth 'cldepth'
!   as a modified cloud parameter for TV presentation.
!
! Method:
!   Normalize the vertical integal of cloud cover by 700 hPa.
!
!------------------------------------------------------------------------------
! Parameter list:

INTEGER (KIND=iintegers), INTENT (IN)    ::    &
  ie, je, ke        ! array dimensions of the I/O-fields

REAL (KIND=wp),     INTENT (OUT)   ::  &
  cldepth(ie,je)    ! normalized cloud depth

REAL (KIND=wp),     INTENT (IN)    ::  &
  clc_sgs(ie,je,ke),  & ! subgrid-scale cloud cover in model layers
  clc_con(ie,je,ke),  & ! convective    cloud cover in model layers
  dp0(ie,je,ke)         ! pressure thickness of model layers

!------------------------------------------------------------------------------
! Local variables

INTEGER (KIND=iintegers) ::    &
  k                 ! Index for vertical loop

!------------------------------------------------------------------------------

! Initialise the output array with 0
 
  cldepth  (:,:) = 0.0_wp
 
! Compute the cloud depth and scale it between 0 and 1

  DO k = 1, ke
    cldepth (:,:) = cldepth (:,:) + ( clc_sgs(:,:,k) &
             +  clc_con(:,:,k)*(1.0_wp -  clc_sgs(:,:,k)) ) * dp0(:,:,k)
  ENDDO

  cldepth (:,:) = MIN (1.0_wp, cldepth (:,:)/700.0E2_wp)

END SUBROUTINE calcldepth

!==============================================================================
!==============================================================================

SUBROUTINE calclmod ( clct_mod, clc_sgs, clc_con, p0hl, pi, ie, je, ke )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine calculates the total cloud cover 'clt_mod'
!   as a modified cloud parameter for TV presentation.
!
! Method:
!   Computation of the total cloud cover by minimum overlapping and
!   multiply it by an empirical 'red'-factor
!
!------------------------------------------------------------------------------
! Parameter list:

INTEGER (KIND=iintegers), INTENT (IN)    ::    &
  ie, je, ke        ! array dimensions of the I/O-fields

REAL (KIND=wp),     INTENT (OUT)   ::  &
  clct_mod(ie,je)   ! normalized cloud depth

REAL (KIND=wp),     INTENT (IN)    ::  &
  clc_sgs(ie,je,ke),  & ! subgrid-scale cloud cover in model layers
  clc_con(ie,je,ke),  & ! convective    cloud cover in model layers
  p0hl (ie,je,ke+1),  & ! reference pressure at half levels
  pi                    ! circle constant

!------------------------------------------------------------------------------
! Local variables

INTEGER (KIND=iintegers) ::    &
  mclbas(ie,je),  & ! base index of significant cloudiness (>= 0.5)
  mcltop(ie,je),  & ! top  Index of signigicant cloudiness (> =0.5)
  i,j,k             ! loop indexes

REAL (KIND=wp)           ::    &
  zclc(ie,je,ke ), & ! total cloud cover in a layer
  zpclbas(ie,je) , & ! pressure at base of significant cloudiness (>= 0.5)
  zpcltop(ie,je) , & ! pressure at top  of significant cloudiness (>= 0.5)
  zcl_bkn, zpclbas_min, zpclbas_max, zred  ! local storage

!------------------------------------------------------------------------------

! Initialise the output array with 0
  clct_mod (:,:) = 0.0_wp

! Compute the pressure at half levels (layer interfaces) and cloud cover
  DO k = 1, ke
    zclc(:,:,k) = clc_sgs(:,:,k)  &
               +  clc_con(:,:,k)*(1.0_wp -  clc_sgs(:,:,k))
  ENDDO

! Determine the j3-indices of the base and top of significant
! cloudiness (clc >= 0.5)

  zcl_bkn = 0.5_wp

  mclbas(:,:) = 1
  mcltop(:,:) = ke

  DO k = 2, ke
    DO j = 1, je
      DO i = 1, ie
        IF (zclc(i,j,k) >= zcl_bkn) THEN
          mclbas(i,j) = MAX (k + 1, mclbas(i,j) )
          mcltop(i,j) = MIN (k    , mcltop(i,j) )
        ENDIF
      ENDDO
    ENDDO
  ENDDO

! Interpolate the cloud top to the pressure where the cloud cover
! is less than "zcl_bkn"; interpolate the cloud base to the pressure
! where the cloud cover is greater than "zcl_bkn".
! Derive the cloud depth

  DO j = 1, je
    DO i = 1, ie
      IF (mclbas(i,j) == 1) THEN ! No cloud at this grid point
        zpclbas(i,j) = 0.0_wp
        zpcltop(i,j) = 0.0_wp
      ELSE                       ! Interpolate base and top pressure
        zpcltop(i,j) = (zclc(i,j,mcltop(i,j)) - zcl_bkn)     &
                / MAX (0.001_wp, zclc(i,j,mcltop(i,j))   &
                -  zclc(i,j,mcltop(i,j)-1))                  &
                * (p0hl (i,j,mcltop(i,j)-1) - p0hl (i,j,mcltop(i,j)) )  &
                + p0hl (i,j,mcltop(i,j))

        IF (mclbas(i,j) >= ke) THEN
            zpclbas (i,j) = p0hl (i,j,ke+1)
        ELSE
            zpclbas (i,j) = (zclc(i,j,mclbas(i,j)-1) - zcl_bkn)    &
                 /  MAX (0.001_wp, zclc(i,j,mclbas(i,j)-1)     &
                 -  zclc(i,j,mclbas(i,j)))                         &
                 *  (p0hl (i,j,mclbas(i,j)+1) - p0hl (i,j,mclbas(i,j)) ) &
                 +  p0hl (i,j,mclbas(i,j))
        ENDIF
      ENDIF
    ENDDO
  ENDDO

! Compute the modified total cloud cover; do not take the high clouds
! into account if they are the only clouds present at this grid
! point; the computation of the cloud cover uses maximum overlapping

  DO k = 1, ke
    clct_mod(:,:) = MAX (clct_mod(:,:), zclc(:,:,k))
  ENDDO

  zpclbas_min = 200.0E2_wp
  zpclbas_max = 600.0E2_wp

  DO j = 1, je
    DO i = 1, ie
      zred  = 1.0_wp
      IF (zpclbas(i,j) < zpclbas_min) THEN
        zred  = 0.0_wp
      ELSE IF (zpclbas(i,j) < zpclbas_max) THEN
        zred  = MAX (0.0_wp, COS(0.5_wp*pi/(zpclbas_min-zpclbas_max) &
                * (zpclbas(i,j) - zpclbas_max)) )
      ENDIF
      clct_mod(i,j) = zred*clct_mod(i,j)
    ENDDO
  ENDDO

END SUBROUTINE calclmod

!==============================================================================
!==============================================================================

SUBROUTINE caliq ( iq, rho, hhl, q, ie, je, ke )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine calculates the vertically integrated mass of a humidity
!   variable q (concentration)
!
! Method:
!   Vertical integral of 'rho*q' using time level nnow on input for q 
!
!------------------------------------------------------------------------------
! Parameter list:

INTEGER (KIND=iintegers), INTENT (IN)    ::    &
  ie, je, ke        ! array dimensions of the I/O-fields

REAL (KIND=wp),     INTENT (OUT)   ::  &
  iq  (ie,je)       ! vertically integrated mass for concentration q

REAL (KIND=wp),     INTENT (IN)    ::  &
  rho (ie,je,ke)  ,  & ! air density
  hhl (ie,je,ke+1),  & ! height of model half levels
  q   (ie,je,ke)       ! humidity mass concentration at time-level nnow

!------------------------------------------------------------------------------
! Local variables

INTEGER (KIND=iintegers) ::    &
  k              ! loop indexes

!------------------------------------------------------------------------------
! Begin of subroutine
!------------------------------------------------------------------------------

! integrate from top to bottom
    iq(:,:) = 0.0_wp

    DO k = 1, ke
       iq(:,:) = iq(:,:) + rho(:,:,k)*(hhl(:,:,k)-hhl(:,:,k+1)) * q(:,:,k)
    ENDDO

END SUBROUTINE caliq   

!==============================================================================
!==============================================================================

SUBROUTINE calztd ( ztd, zwd, zhd, rho, hhl, q, ps, t, hsurf, phi,        &
                     pi, ie, je, ke )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine calculates the vertically integrated dry, wet and total
!   refractivity as equivalent to zenith hydrostatic (ZHD), wet (ZWD) and total
!   (ZTD) delay in GPS
!
! Method:
!   The total, dry and wet delay are calculated following Bevis et al.
!   (1992, J. Geophys. Res. 97, p. 15 787 - 15 801)
!
!------------------------------------------------------------------------------
! Parameter list:

INTEGER (KIND=iintegers), INTENT (IN)    ::    &
  ie, je, ke        ! array dimensions of the I/O-fields

REAL (KIND=wp),     INTENT (IN)    ::  &
  ps  (ie,je)     , & ! surface pressure [Pa]
  phi (ie,je)     , & ! geographical latitude [deg N]
  t   (ie,je)     , & ! temperature at lowest model level [K]
  hsurf(ie,je)    , & ! geometrical height of surface [m]
  rho (ie,je,ke)  , & ! air density
  hhl (ie,je,ke+1), & ! height of model half levels
  q   (ie,je,ke)  , & ! humidity mass concentration at time-level nnow
  pi                  ! circle constant

REAL (KIND=wp),     INTENT (OUT)   ::  &
  zhd (ie,je)      , & ! integrated dry atmospheric refractivity
  zwd (ie,je)      , & ! integrated wet atmospheric refractivity
  ztd (ie,je)          ! integrated total atmospheric refractivity

!------------------------------------------------------------------------------
! Local variables

INTEGER (KIND=iintegers) ::    &
  i, j, k              ! loop indexes

REAL (KIND=wp)        ::    &
  ef              , & ! gravitational acceleration factor
  ak              , & ! conversion parameter
  Tm              , & ! weighted mean atmospheric temperature [K]
  iq  (ie,je)     , & ! vertically integrated mass for concentration q
  cc1, cc2, cc3       !local storage

!------------------------------------------------------------------------------
! Begin of subroutine
!------------------------------------------------------------------------------

  ! integrate from top to bottom
  iq(:,:) = 0.0_wp
  DO k = 1, ke
    iq(:,:) = iq(:,:) + rho(:,:,k)*(hhl(:,:,k)-hhl(:,:,k+1)) * q(:,:,k)
  ENDDO

  DO j = 1, je
    DO i = 1, ie
      cc1     = 0.0022768_wp
      cc2     = 0.00266_wp
      cc3     = 0.00028_wp

      ef      = 1.0_wp                                                &
                - cc2 * COS  (2.0_wp * phi(i,j) * pi / 180.0_wp) &
                - cc3 * hsurf(i,j) / 1000.0_wp

      Tm      = 70.2_wp + 0.72_wp * t(i,j)
      ak      = (4.6151_wp * (3.776E5_wp/Tm + 22_wp))*1E-6_wp

      zhd(i,j) = ( cc1 * ps(i,j) / 100.0_wp ) / ef
      zwd(i,j) = iq(i,j) * ak
      ztd(i,j) = zhd(i,j) + zwd(i,j)
    ENDDO
  ENDDO

END SUBROUTINE calztd

!==============================================================================
!==============================================================================

SUBROUTINE radar_lm_ray( ix, iy, iz, pi, rho_w, rho_ice, K_w, K_ice, t0_melt,&
                         klv850, myproc, igscp, idebug,                      &
                         T, q_cloud, q_rain, q_ice, q_snow, q_grau,          &
                         mue_rain_c_in, z_radar, z_radar_850, z_radar_cmax,  &
                         z_radar_dp)

!------------------------------------------------------------------------------
!
! Description:  Calculation of grid point values for effective radar
!               reflectivity factor Z in dBZ.
!
! Method:       Rayleigh-Approximation for the Back-Scattering
!               (no attenuation!), Debye-Approximation for the
!               effective refractive index of two-component ice-air-
!               mixture particles (dry ice, snow, and graupel).
!               Melting particles by substitution of the ice substance
!               with water when T_air > 273.16 K. For the applied
!               Rayleigh scattering, this is equal to assuming instantaneous
!               melting, because only the square of the total water mass of the particle
!               counts ( = the total square of its dipole moment).
!               
! Inputs:       ix, iy, iz   : field dimensions
!               pi           : Pi
!               rho_w        : bulk density of pure water [kg/m**3]
!               rho_ice      : bulk density of pure ice   [kg/m**3]
!               K_w          : dielectric constant of water
!               K_ice        : dielectric constant of ice
!               t0_melt      : melting temperature of ice
!               klv850       : model k-index of the height level next to 850 hPa
!               myproc       : index of this PE in the processor grid
!               igscp        : itype of grid scale precip (1,2,3,4)
!               T            : temperature field        [K]
!               q_cloud      : cloud water mass density [kg/m**3] 
!               q_rain       : rain water mass density  [kg/m**3] 
!               q_ice        : cloud ice mass density   [kg/m**3] 
!               q_snow       : snow mass density        [kg/m**3] 
!   OPTIONAL:   q_grau       : graupel mass density     [kg/m**3]
!   OPTIONAL:   mue_rain_c_in: assumed constant value of mue for rain for the DSD model
!                                N(D) = N_0 D^mue exp(-lambda*D)    (D = diameter)
!                              DEFAULT: depends on igsp!
!
! Outputs:      z_radar      : 3D field of Z                         [dBZ]
!               z_radar_850  : Z in the height level next to 850 hPa [dBZ]
!               z_radar_cmax : 2D field of column max Z              [dBz]  
! 
!
!------------------------------------------------------------------------------

! Parameterlist
!--------------

INTEGER(KIND=iintegers), INTENT(IN)       :: ix, iy, iz

REAL(KIND=wp)          , INTENT(IN)       :: pi, rho_w, rho_ice, K_w, K_ice, t0_melt

INTEGER(KIND=iintegers), INTENT(IN)       :: klv850, myproc, igscp, idebug

REAL(KIND=wp)          , INTENT(IN)       :: T(ix,iy,iz),          &
                                             q_cloud(ix,iy,iz),    &
                                             q_rain(ix,iy,iz),     &
                                             q_ice(ix,iy,iz),      &
                                             q_snow(ix,iy,iz)

REAL(KIND=wp),     INTENT(IN),  OPTIONAL  :: q_grau(ix,iy,iz), &
                                             mue_rain_c_in

REAL(KIND=wp),     INTENT(OUT), OPTIONAL  :: z_radar(ix,iy,iz)

REAL(KIND=wp),     INTENT(OUT), OPTIONAL  :: z_radar_850(ix,iy),   &
                                             z_radar_cmax(ix,iy)

REAL(KIND=dp),     INTENT(OUT), OPTIONAL  :: z_radar_dp(ix,iy,iz)

! Local Variables
!----------------

REAL(KIND=wp),     PARAMETER ::  &
  zxstar = 2.60e-10_wp,      & ! separating mass between cloud and rain
  eps    = 1.00E-15_wp

REAL(KIND=wp)     ::  z_local(ix,iy,iz)
REAL(KIND=wp)     ::  z_cloud(ix,iy,iz), z_rain(ix,iy,iz), z_snow(ix,iy,iz), &
                      z_ice(ix,iy,iz), z_grau(ix,iy,iz)

INTEGER(KIND=iintegers)       :: i, j, k
INTEGER(KIND=iintegers), SAVE :: firstcall = 0


REAL   (KIND=wp)              :: z_ice_dry, mom_fak

REAL (KIND=wp)                :: mma(10), mmb(10)
REAL (KIND=wp)                :: ztc, m2s, m3s, alf, bet, hlp, zn0s, D_i_mono, D_c_mono
INTEGER(KIND=iintegers)       :: nn

REAL   (KIND=wp),     SAVE    :: z_r, z_g, z_s, p_r, p_s, p_g
REAL   (KIND=wp),     SAVE    :: nor, nog, amg, bmg, nos, ams, bms, &
                                 ami, bmi, x_i_mono, x_c_mono, mue_rain_c

INTEGER (KIND=iintegers), PARAMETER :: isnow_n0temp   = 2

REAL   (KIND=wp),     PARAMETER ::       &
    zn0s1 = 13.5_wp * 5.65E5_wp, & ! parameter in N0S(T)
    zn0s2 = -0.107_wp                  ! parameter in N0S(T), Field et al

!------------------------------------------------------------------------------

! Coeffs for moment relation based on 2nd moment (Field 2005)
  mma = (/   5.065339_wp, -0.062659_wp, -3.032362_wp, 0.029469_wp, -0.000285_wp, &
             0.312550_wp,  0.000204_wp,  0.003199_wp, 0.000000_wp, -0.015952_wp /)
  mmb = (/   0.476221_wp, -0.015896_wp,  0.165977_wp, 0.007468_wp, -0.000141_wp, &
             0.060366_wp,  0.000079_wp,  0.000594_wp, 0.000000_wp, -0.003577_wp /)

!------------------------------------------------------------------------------

  IF (idebug > 5) THEN
    PRINT *, 'radar_lm_ray: igscp = ',igscp
  ENDIF

  z_cloud = 0.0_wp
  z_rain = 0.0_wp
  z_ice = 0.0_wp
  z_snow = 0.0_wp
  z_grau = 0.0_wp
  z_local = 0.0_wp

  z_ice_dry = (rho_w/rho_ice)**2 * K_ice/K_w
  mom_fak = (6.0_wp / (pi * rho_w))**2


  IF (firstcall /= 1) THEN
    IF (igscp < 3) THEN
      D_c_mono = 2.0E-5_wp
      x_c_mono = pi * rho_w / 6.0_wp * D_c_mono**3 
      D_i_mono = 2.0E-4_wp
      ami = 130.0_wp
      bmi = 3.0_wp
      x_i_mono = ami * D_i_mono**bmi
      ams = 0.069_wp
      bms = 2.0_wp
      nor = 8.0E6_wp
      nos = 8.0E5_wp
      p_r = 7.0_wp / 4.0_wp
      p_s = (2.0_wp*bms+1.0_wp)/(bms+1.0_wp)
      z_r = nor*gamma_fct(7.0_wp) * (pi*rho_w*nor)**(-p_r)
      z_s = mom_fak*ams**2 * nos*gamma_fct(2.0_wp*bms+1.0_wp) * (ams*nos*gamma_fct(bms+1.0_wp))**(-p_s)
      IF (idebug > 5) THEN
        WRITE (*, *) "radar_lm_ray: cloud ice scheme (using rain and snow)"
        WRITE (*,'(A,F10.3)') '     p_r = ', p_r
        WRITE (*,'(A,F10.3)') '     z_r = ', z_r
        WRITE (*,'(A,F10.3)') '     p_s = ', p_s
        WRITE (*,'(A,F10.3)') '     z_s = ', z_s
      ENDIF
      firstcall = 1
    ELSEIF (igscp == 3) THEN
      IF (.NOT. PRESENT(mue_rain_c_in) ) THEN
        mue_rain_c = 0.0_wp
      ELSE
        mue_rain_c = mue_rain_c_in
      END IF
      D_c_mono = 2.0E-5_wp
      x_c_mono = pi * rho_w / 6.0_wp * D_c_mono**3 
      D_i_mono = 2.0E-4_wp
      ami = 130.0_wp
      bmi = 3.0_wp
      x_i_mono = ami * D_i_mono**bmi
      ams = 0.069_wp
      bms = 2.0_wp
!!$ UB: moment relation of Ulbrich consistent to hydci_pp (the "*0.1" is a "needed" bug):
      nor = 8.0e6_wp * EXP(3.2_wp*mue_rain_c) * (0.01_wp)**(-mue_rain_c) * 0.1_wp
!!$      nos = 8.E5_wp  ! now computed below!
      p_r = (7.0_wp+mue_rain_c) / (4.0_wp+mue_rain_c)
      p_s = (2.0_wp*bms+1.0_wp)/(bms+1.0_wp)
      z_r = nor*gamma_fct(7.0_wp+mue_rain_c) * (pi*rho_w*nor)**(-p_r)
!!$ UB: the factor nos/nos**(-p_s) is a leftover from times when we had
!!$     constant n0s. However, this factor comes into play nowadays below,
!!$     so we have to cancel it here! Effect: Z_snow increases by 40 - 50 dB,
!!!     from much too low values to reasonable ones!!! 
!!$      z_s = mom_fak*ams**2 * nos*gamma_fct(2*bms+1) * (ams*nos*gamma_fct(bms+1))**(-p_s)
      z_s = mom_fak*ams**2 * gamma_fct(2.0_wp*bms+1.0_wp) * (ams*gamma_fct(bms+1.0_wp))**(-p_s)
      IF (idebug > 5) THEN
        WRITE (*, *) "radar_lm_ray: cloud ice scheme (using rain and snow)"
        WRITE (*,'(A,F10.3)') '     p_r = ', p_r
        WRITE (*,'(A,F10.3)') '     z_r = ', z_r
        WRITE (*,'(A,F10.3)') '     p_s = ', p_s
        WRITE (*,'(A,F10.3)') '     z_s = ', z_s
      ENDIF
      firstcall = 1
    ELSEIF (igscp == 4) THEN
      IF (.NOT. PRESENT(mue_rain_c_in) ) THEN
        mue_rain_c = 0.5_wp
      ELSE
        mue_rain_c = mue_rain_c_in
      END IF
      D_c_mono = 2.0E-5_wp
      x_c_mono = pi * rho_w / 6.0_wp * D_c_mono**3 
      D_i_mono = 2.0E-4_wp
      ami = 130.0_wp
      bmi = 3.0_wp
      x_i_mono = ami * D_i_mono**bmi
      ams = 0.038_wp
      bms = 2.0_wp
      amg = 169.6_wp
      bmg = 3.1_wp
!!$ UB: moment relation of Ulbrich consistent to hydci_pp_gr (the "*0.1" is a "needed" bug):
      nor = 8.0E6_wp * EXP(3.2_wp*mue_rain_c) * (0.01_wp)**(-mue_rain_c) * 0.1_wp
!!$      nos = 8.E5_wp  ! now computed below!
      nog = 4.E6_wp
      p_r = (7.0_wp+mue_rain_c) / (4.0_wp+mue_rain_c)
      p_s = (2.0_wp*bms+1.0_wp)/(bms+1.0_wp)
      p_g = (2.0_wp*bmg+1.0_wp)/(bmg+1.0_wp)
      z_r = nor*gamma_fct(7.0_wp+mue_rain_c) * &
            (pi*rho_w*nor*gamma_fct(4.0_wp+mue_rain_c)/6.0_wp)**(-p_r)
!!$ UB: the factor nos/nos**(-p_s) is a leftover from times when we had
!!$     constant n0s. However, this factor comes into play nowadays below,
!!$     so we have to cancel it here! Effect: Z_snow increases by 40 - 50 dB,
!!!     from much too low values to reasonable ones!!! 
!!$      z_s = mom_fak*ams**2 * nos*gamma_fct(2*bms+1.0_wp) *       &
!!$                                (ams*nos*gamma_fct(bms+1.))**(-p_s)
      z_s = mom_fak*ams**2 * gamma_fct(2.0_wp*bms+1.0_wp) *       &
                                (ams*gamma_fct(bms+1.0_wp))**(-p_s)
      z_g = mom_fak*amg**2 * nog*gamma_fct(2.0_wp*bmg+1.0_wp) *       &
                                (amg*nog*gamma_fct(bmg+1.0_wp))**(-p_g)
      IF (idebug > 5) THEN
        WRITE (*, *) "radar_lm_ray: graupel scheme (using rain, snow, graupel)"
        WRITE (*,'(A,F10.3)') '     p_r = ', p_r
        WRITE (*,'(A,F10.3)') '     z_r = ', z_r
        WRITE (*,'(A,F10.3)') '     p_s = ', p_s
        WRITE (*,'(A,F10.3)') '     z_s = ', z_s
        WRITE (*,'(A,F10.3)') '     p_g = ', p_g
        WRITE (*,'(A,F10.3)') '     z_g = ', z_g
      ENDIF
      firstcall = 1
    ELSE
      WRITE (*,*) 'ERROR in radar_lm_ray(), pp_utilities.f90: itype_gscp = ',igscp,' not implemented!'
      STOP
    ENDIF
  ENDIF

  DO k = 1, iz
    DO j = 1, iy
      DO i = 1, ix

        ! .. cloud droplets (monodisperse size distribution):
        z_cloud(i,j,k) = mom_fak * q_cloud(i,j,k) * x_c_mono
        z_local(i,j,k) = z_cloud(i,j,k)

        ! .. rain
        z_rain(i,j,k) = z_r * q_rain(i,j,k)**p_r
        z_local(i,j,k) = z_local(i,j,k) + z_rain(i,j,k)

        ! .. cloud ice (monodisperse size distribution):
        IF (T(i,j,k) < 273.15_wp) THEN
          z_ice(i,j,k) = z_ice_dry * mom_fak * q_ice(i,j,k) * x_i_mono
        ELSE
          z_ice(i,j,k) = mom_fak * q_ice(i,j,k) * x_i_mono
        ENDIF
        z_local(i,j,k) = z_local(i,j,k) + z_ice(i,j,k)

        ! .. snow
        IF (q_snow(i,j,k) > 1.0E-10_wp) THEN

          IF (isnow_n0temp == 1) THEN
            ! Calculate n0s using the temperature-dependent
            ! formula of Field et al. (2005)
            ztc = T(i,j,k) - t0_melt
            ztc = MAX(MIN(ztc,0.0_wp),-40.0_wp)
            zn0s = zn0s1*EXP(zn0s2*ztc)
            zn0s = MIN(zn0s,1e9_wp)
            zn0s = MAX(zn0s,1e6_wp)
          ELSEIF (isnow_n0temp == 2) THEN
            ! Calculate n0s using the temperature-dependent moment
            ! relations of Field et al. (2005)
            ztc = T(i,j,k) - t0_melt
            ztc = MAX(MIN(ztc,0.0_wp),-40.0_wp)
            nn  = 3
            hlp = mma(1)      +mma(2)*ztc      +mma(3)*nn       +mma(4)*ztc*nn+mma(5)*ztc**2 &
              + mma(6)*nn**2+mma(7)*ztc**2*nn+mma(8)*ztc*nn**2+mma(9)*ztc**3+mma(10)*nn**3
            alf = 10.0_wp**hlp
            bet = mmb(1)      +mmb(2)*ztc      +mmb(3)*nn       +mmb(4)*ztc*nn+mmb(5)*ztc**2 &
              + mmb(6)*nn**2+mmb(7)*ztc**2*nn+mmb(8)*ztc*nn**2+mmb(9)*ztc**3+mmb(10)*nn**3
!!$ UB: caution! Here is the exponent bms=2.0 hardwired! not ideal!
            m2s = q_snow(i,j,k) / ams
            m3s = alf*EXP(bet*LOG(m2s))
            hlp  = zn0s1*EXP(zn0s2*ztc)
!!$ UB: the 13.5 is actually 3^3 / gamma(3) ...
            zn0s = 13.50_wp * m2s*(m2s/m3s)**3
            zn0s = MAX(zn0s,0.5_wp*hlp)
            zn0s = MIN(zn0s,1e2_wp*hlp)
            zn0s = MIN(zn0s,1e9_wp)
            zn0s = MAX(zn0s,8e5_wp)
          ELSE
            ! Old constant n0s
            zn0s = 8.0e5_wp
          ENDIF
          hlp = z_s * zn0s**(1.0_wp-p_s)

          IF (T(i,j,k) < 273.15_wp) THEN
            z_snow(i,j,k) = hlp * q_snow(i,j,k)**p_s * z_ice_dry
          ELSE
            z_snow(i,j,k) = hlp * q_snow(i,j,k)**p_s
          ENDIF

          z_local(i,j,k) = z_local(i,j,k) + z_snow(i,j,k)

        END IF

        ! .. graupel
        IF (igscp >= 4) THEN
          IF ( PRESENT(q_grau) ) THEN
            IF (T(i,j,k) < 273.15_wp) THEN
              z_grau(i,j,k) = z_g * q_grau(i,j,k)**p_g * z_ice_dry
            ELSE
              z_grau(i,j,k) = z_g * q_grau(i,j,k)**p_g
            ENDIF
            z_local(i,j,k) = z_local(i,j,k) + z_grau(i,j,k)
          ENDIF
        ENDIF

      ENDDO
    ENDDO
  ENDDO

  z_local = 10.0_wp / LOG(10.0_wp) * LOG(z_local * 1.E18_wp + eps)

  IF (idebug > 5) THEN
    z_cloud= 10.0_wp / LOG(10.0_wp) * LOG(z_cloud* 1.E18_wp + eps)
    z_rain = 10.0_wp / LOG(10.0_wp) * LOG(z_rain * 1.E18_wp + eps)
    z_ice  = 10.0_wp / LOG(10.0_wp) * LOG(z_ice  * 1.E18_wp + eps)
    z_snow = 10.0_wp / LOG(10.0_wp) * LOG(z_snow * 1.E18_wp + eps)
    z_grau = 10.0_wp / LOG(10.0_wp) * LOG(z_grau * 1.E18_wp + eps)

    WRITE (*,'(A,f10.1)') "     MAX dBZ total = ", MAXVAL(z_local)
    WRITE (*,'(A,f10.1)') "     MAX dBZ cloud = ", MAXVAL(z_cloud)
    WRITE (*,'(A,f10.1)') "     MAX dBZ rain  = ", MAXVAL(z_rain)
    WRITE (*,'(A,f10.1)') "     MAX dBZ ice   = ", MAXVAL(z_ice)
    WRITE (*,'(A,f10.1)') "     MAX dBZ snow  = ", MAXVAL(z_snow)
    WRITE (*,'(A,f10.1)') "     MAX dBZ grau  = ", MAXVAL(z_grau)
  ENDIF

  IF ( PRESENT(z_radar) ) THEN
    DO k = 1, iz
      DO j = 1, iy
        DO i = 1, ix
          z_radar(i,j,k) = z_local(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF ( PRESENT(z_radar_dp) ) THEN
    DO k = 1, iz
      DO j = 1, iy
        DO i = 1, ix
          z_radar_dp(i,j,k) = REAL(z_local(i,j,k),kind=dp)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF ( PRESENT(z_radar_850) ) THEN
    DO j = 1, iy
      DO i = 1, ix
        z_radar_850(i,j) = z_local(i,j,klv850)
      ENDDO
    ENDDO
  ENDIF

  IF ( PRESENT(z_radar_cmax) ) THEN
    DO j = 1, iy
      DO i = 1, ix
        z_radar_cmax(i,j) = MAXVAL(z_local(i,j,:))
      ENDDO
    ENDDO
  ENDIF

END SUBROUTINE radar_lm_ray

!==============================================================================
!==============================================================================

FUNCTION gamma_fct_orig(x)

!------------------------------------------------------------------------------
!
! Description:
!  Gamma-function from Numerical Recipes (F77)
! Method:
!
!------------------------------------------------------------------------------

REAL   (KIND=wp)     :: gamma_fct_orig

REAL   (KIND=wp)     :: cof(6) = (/76.18009173_wp, -86.50532033_wp, &
                                   24.01409822_wp, -1.231739516_wp, &
                                   0.120858003E-2_wp, -0.536382E-5_wp/)
REAL   (KIND=wp)     :: stp=2.50662827465_wp,                           &
                        x, xx, tmp, ser, gamma

INTEGER(KIND=iintegers) ::  j

    xx  = x  - 1.0_wp
    tmp = xx + 5.5_wp
    tmp = (xx + 0.5_wp) * LOG(tmp) - tmp
    ser = 1.0_wp
    DO j = 1, 6
       xx  = xx  + 1.0_wp
       ser = ser + cof(j) / xx
    ENDDO
    gamma = tmp + LOG(stp*ser)
    gamma = EXP(gamma)

    gamma_fct_orig = gamma

END FUNCTION gamma_fct_orig

!==============================================================================
!==============================================================================

FUNCTION gamma_fct(x)
  !*******************************************************************************
  !                                                                              *
  !       Gammafunktion aus Numerical Recipes (F77)                              *
  !       (umformuliert wg. inlining und Vektorisierung)
  !*******************************************************************************
  IMPLICIT NONE

  REAL(KIND=wp) :: gamma_fct

  REAL(KIND=wp), INTENT(in) :: x

  REAL(KIND=wp) :: tmp, p

  REAL(KIND=wp), PARAMETER :: c1 =  76.18009173_wp
  REAL(KIND=wp), PARAMETER :: c2 = -86.50532033_wp
  REAL(KIND=wp), PARAMETER :: c3 =  24.01409822_wp
  REAL(KIND=wp), PARAMETER :: c4 = -1.231739516_wp
  REAL(KIND=wp), PARAMETER :: c5 =  0.120858003e-2_wp
  REAL(KIND=wp), PARAMETER :: c6 = -0.536382e-5_wp
  REAL(KIND=wp), PARAMETER :: stp = 2.50662827465_wp

  tmp = x + 4.5_wp;
  p = stp * (1.0_wp + c1/x + c2/(x+1.0_wp) + c3/(x+2.0_wp) + c4/(x+3.0_wp) + c5/(x+4.0_wp) + c6/(x+5.0_wp))
  gamma_fct = p * EXP( (x-0.5_wp) * LOG(tmp) - tmp )

END FUNCTION gamma_fct

!==============================================================================
!==============================================================================

SUBROUTINE potential_vorticity_rho( ie, je, ke, eddlon, eddlat, r_earth,     &
   fc, fccos, sqrtg_r_s, dzeta_dlam, dzeta_dphi, curl1, curl2, curl3,        &
   lmetr, Theta, u, v, w, pot_vort_rho )

!------------------------------------------------------------------------------
!
! Description:
!   calculate  rho * PV  =  del Theta * ( curl v + 2 Omega )
!
! this subroutine is for external use outside of the COSMO-model
! (e.g. in field_extra, ...) 
! (for use in the COSMO-model: see subr. 'potential_vorticity_rho')
!
!------------------------------------------------------------------------------

INTEGER (KIND=iintegers), INTENT(IN)  ::  &
  ie, je, ke

REAL (KIND=wp),     INTENT(IN)  ::  &
  eddlon,     & ! 1 / dlon, dlon            (in degrees)
  eddlat,     & ! 1 / dlat, dlat            (in degrees)
  r_earth       ! mean radius of the earth  (m)

REAL (KIND=wp),     INTENT(IN)  ::  &
  fc   (ie,je),         & ! coriolis-parameter (radial comp.)         ( 1/s )
  fccos(ie,je)            ! coriolis-parameter (zonal  comp.)         ( 1/s )

REAL (KIND=wp),     INTENT(IN)  ::  &
  sqrtg_r_s (ie,je,ke), & ! 1 / sqrt(G)       (at the scalar position)
  dzeta_dlam(ie,je,ke), & ! d zeta / d lambda (at the scalar position)
  dzeta_dphi(ie,je,ke), & ! d zeta / d phi    (at the scalar position)
  curl1     (ie,je,ke), & ! contravariant components of curl v
  curl2     (ie,je,ke), & ! 
  curl3     (ie,je,ke)    ! 

LOGICAL, INTENT(IN) ::    &
  lmetr         ! with metrics of the spherical earth

REAL (KIND=wp),     INTENT(IN)  ::  &
  Theta(ie,je,ke),      & ! potential temperature  (K)
  u    (ie,je,ke),      & ! meridional velocity    (m/s)
  v    (ie,je,ke)         ! zonal velocity         (m/s)

REAL (KIND=wp),     INTENT(IN)  ::  &
  w    (ie,je,ke+1)       ! vertical velocity      (m/s)

REAL (KIND=wp),     INTENT(OUT) ::  &
  pot_vort_rho(ie,je,ke)  ! potential vorticity * rho    (K s-1 m-1 )

! Local variables
REAL (KIND=wp)     :: dTheta_dx1, dTheta_dx2, dTheta_dx3
REAL (KIND=wp)     :: cor2, cor3  ! physical components of 2 * earth angular velocity vector
REAL (KIND=wp)     :: r_earth_inv

INTEGER (KIND=iintegers) :: i, j, k

!------------------------------------------------------------------------------

  r_earth_inv = 1.0_wp / r_earth

  pot_vort_rho(:,:,:) = 0.0_wp

  DO k=2, ke-1
    DO j=2, je-1
      DO i=2, ie-1

        dTheta_dx1 = ( Theta(i+1,j,k) - Theta(i-1,j,k) ) * 0.5_wp * eddlon   &
          &        + ( Theta(i,j,k+1) - Theta(i,j,k-1) ) * 0.5_wp * dzeta_dlam(i,j,k)

        dTheta_dx2 = ( Theta(i,j+1,k) - Theta(i,j-1,k) ) * 0.5_wp * eddlat   &
          &        + ( Theta(i,j,k+1) - Theta(i,j,k-1) ) * 0.5_wp * dzeta_dphi(i,j,k)

        dTheta_dx3 = ( Theta(i,j,k+1) - Theta(i,j,k-1) ) * 0.5_wp * (-sqrtg_r_s(i,j,k))

        cor2 = fccos(i,j) * r_earth_inv
        cor3 = 0.25_wp * ( fc(i,j) + fc(i-1,j) + fc(i,j-1) + fc(i-1,j-1) )

        pot_vort_rho(i,j,k) = dTheta_dx1 * ( curl1(i,j,k)        )  &
          &                 + dTheta_dx2 * ( curl2(i,j,k) + cor2 )  &
          &                 + dTheta_dx3 * ( curl3(i,j,k) + cor3 )

      END DO
    END DO
  END DO

  ! Top layer
  k=1
  DO j=2, je-1
    DO i=2, ie-1

      dTheta_dx1 = ( Theta(i+1,j,k) - Theta(i-1,j,k) ) * 0.5_wp * eddlon   &
        &        + ( Theta(i,j,k+1) - Theta(i,j,k) ) * dzeta_dlam(i,j,k)

      dTheta_dx2 = ( Theta(i,j+1,k) - Theta(i,j-1,k) ) * 0.5_wp * eddlat   &
        &        + ( Theta(i,j,k+1) - Theta(i,j,k) ) * dzeta_dphi(i,j,k)

      dTheta_dx3 = ( Theta(i,j,k+1) - Theta(i,j,k) ) * (-sqrtg_r_s(i,j,k))

      cor2 = fccos(i,j) * r_earth_inv
      cor3 = 0.25_wp * ( fc(i,j) + fc(i-1,j) + fc(i,j-1) + fc(i-1,j-1) )

      pot_vort_rho(i,j,k) = dTheta_dx1 * ( curl1(i,j,k)        )  &
        &                 + dTheta_dx2 * ( curl2(i,j,k) + cor2 )  &
        &                 + dTheta_dx3 * ( curl3(i,j,k) + cor3 )

    END DO
  END DO

  ! Bottom layer
  k=ke
  DO j=2, je-1
    DO i=2, ie-1

      dTheta_dx1 = ( Theta(i+1,j,k) - Theta(i-1,j,k) ) * 0.5_wp * eddlon   &
        &        + ( Theta(i,j,k) - Theta(i,j,k-1) ) * dzeta_dlam(i,j,k)

      dTheta_dx2 = ( Theta(i,j+1,k) - Theta(i,j-1,k) ) * 0.5_wp * eddlat   &
        &        + ( Theta(i,j,k) - Theta(i,j,k-1) ) * dzeta_dphi(i,j,k)

      dTheta_dx3 = ( Theta(i,j,k) - Theta(i,j,k-1) ) * (-sqrtg_r_s(i,j,k))

      cor2 = fccos(i,j) * r_earth_inv
      cor3 = 0.25_wp * ( fc(i,j) + fc(i-1,j) + fc(i,j-1) + fc(i-1,j-1) )

      pot_vort_rho(i,j,k) = dTheta_dx1 * ( curl1(i,j,k)        )  &
        &                 + dTheta_dx2 * ( curl2(i,j,k) + cor2 )  &
        &                 + dTheta_dx3 * ( curl3(i,j,k) + cor3 )

    END DO
  END DO

END SUBROUTINE potential_vorticity_rho

!==============================================================================
#ifdef TWOMOM_SB
!==============================================================================

SUBROUTINE radar_sb_ray( ix, iy, iz, pi, klv850, myproc, T,                  &
                         q_cloud, q_rain, q_ice, q_snow, q_graupel, q_hail,  &
                         n_cloud, n_rain, n_ice, n_snow, n_graupel, n_hail,  &
                         z_radar, z_radar_850, z_radar_cmax, z_radar_dp )

!------------------------------------------------------------------------------
!
! Description:
!
! Method:
!
!------------------------------------------------------------------------------

USE wolken_konstanten, ONLY:  &
     ice_typ, K_w, K_i, T_3, q_krit_aus, rho_ice, rho_w, &
     moment_gamma, rain, graupel, snow, ice, hail

! Parameterlist
!--------------

INTEGER(KIND=iintegers), INTENT(IN)       :: ix, iy, iz

REAL(KIND=wp)          , INTENT(IN)       :: pi
INTEGER(KIND=iintegers), INTENT(IN)       :: klv850, myproc

REAL(KIND=wp)          , INTENT(IN)       :: T(ix,iy,iz),          &
                                             q_cloud(ix,iy,iz),    &
                                             q_rain(ix,iy,iz),     &
                                             q_ice(ix,iy,iz),      &
                                             q_snow(ix,iy,iz),     &
                                             q_graupel(ix,iy,iz),  &
                                             q_hail(ix,iy,iz),     &
                                             n_cloud(ix,iy,iz),    &
                                             n_rain(ix,iy,iz),     &
                                             n_ice(ix,iy,iz),      &
                                             n_snow(ix,iy,iz),     &
                                             n_graupel(ix,iy,iz),  &
                                             n_hail(ix,iy,iz)

REAL(KIND=wp),     INTENT(OUT), OPTIONAL  :: z_radar(ix,iy,iz)

REAL(KIND=wp),     INTENT(OUT), OPTIONAL  :: z_radar_850(ix,iy),   &
                                             z_radar_cmax(ix,iy)

REAL(KIND=dp),     INTENT(OUT), OPTIONAL  :: z_radar_dp(ix,iy,iz)

! Local Variables
!----------------

REAL(KIND=wp),     PARAMETER ::  &
  eps     = 1.00E-15_wp

REAL(KIND=wp)     ::  z_local(ix,iy,iz)

INTEGER(KIND=iintegers)       :: i, j, k
INTEGER(KIND=iintegers), SAVE :: firstcall = 0

REAL   (KIND=wp)              :: T_a, q_c, q_l,                  &
                                 q_r, n_r, x_r,                  &
                                 q_g, n_g, x_g,                  &
                                 q_h, n_h, x_h,                  &
                                 q_s, n_s, x_s,                  &
                                 q_i, n_i, x_i,                  &
                                 z_fak_ice_dry, z_fak_ice_wet,   &
                                 mom_fak, p, q, q_krit_radar

REAL   (KIND=wp),     SAVE    :: z_fak_r, z_fak_g, z_fak_h, z_fak_s, z_fak_i

  z_fak_ice_dry = (rho_w/rho_ice)**2 * K_i/K_w
  z_fak_ice_wet = 1.0_wp

  IF (firstcall /= 1) THEN
    mom_fak = (6.0_wp / (pi * rho_w))**2
    z_fak_r = moment_gamma(rain,2)    * mom_fak
    z_fak_g = moment_gamma(graupel,2) * mom_fak
    z_fak_h = moment_gamma(hail,2)    * mom_fak
    z_fak_s = moment_gamma(snow,2)    * mom_fak
    z_fak_i = moment_gamma(ice,2)     * mom_fak
    firstcall = 1
    IF (myproc == 0) THEN
      WRITE (6, *) "sb_radar_ray:"
      WRITE (6,'(A,D10.3)') "     z_fak_r = ",z_fak_r
      WRITE (6,'(A,D10.3)') "     z_fak_g = ",z_fak_g
      WRITE (6,'(A,D10.3)') "     z_fak_h = ",z_fak_h
      WRITE (6,'(A,D10.3)') "     z_fak_s = ",z_fak_s
      WRITE (6,'(A,D10.3)') "     z_fak_i = ",z_fak_i
      WRITE (6,'(A,D10.3)') "     fak_ice = ",z_fak_ice_dry
    ENDIF
  ENDIF

  q_krit_radar = q_krit_aus

  z_local = 0.0_wp
  DO k= 1, iz
    DO j = 1, iy
      DO i = 1, ix
        T_a = T(i,j,k)
        q_r = q_rain(i,j,k)
        n_r = n_rain(i,j,k)
        q_g = q_graupel(i,j,k)
        n_g = n_graupel(i,j,k)
        q_h = q_hail(i,j,k)
        n_h = n_hail(i,j,k)
        q_s = q_snow(i,j,k)
        n_s = n_snow(i,j,k)
        q_i = q_ice(i,j,k)
        n_i = n_ice(i,j,k)
        q_l = q_cloud(i,j,k) + q_r
        x_r = MIN( MAX(q_r/(n_r+eps),rain%x_min),rain%x_max )
        x_g = MIN( MAX(q_g/(n_g+eps),graupel%x_min),graupel%x_max )
        x_h = MIN( MAX(q_h/(n_h+eps),hail%x_min),hail%x_max )
        x_s = MIN( MAX(q_s/(n_s+eps),snow%x_min),snow%x_max )
        x_i = MIN( MAX(q_i/(n_i+eps),ice%x_min),ice%x_max )
        IF (q_s < q_krit_radar) x_s = 0.0_wp
        IF (q_r < q_krit_radar) x_r = 0.0_wp
        IF (q_g < q_krit_radar) x_g = 0.0_wp
        IF (q_h < q_krit_radar) x_h = 0.0_wp
        IF (q_i < q_krit_radar) x_i = 0.0_wp

        z_local(i,j,k) = z_fak_r * q_r * x_r
        IF (T_a < T_3) THEN
          z_local(i,j,k) = z_local(i,j,k)                            &
            &       + z_fak_g * q_g * x_g * z_fak_ice_dry
          z_local(i,j,k) = z_local(i,j,k)                            &
            &       + z_fak_h * q_h * x_h * z_fak_ice_dry
          z_local(i,j,k) = z_local(i,j,k)                            &
            &       + z_fak_s * q_s * x_s * z_fak_ice_dry
          z_local(i,j,k) = z_local(i,j,k)                            &
            &       + z_fak_i * q_i * x_i * z_fak_ice_dry
        ELSE
          z_local(i,j,k) = z_local(i,j,k)                            &
            &       + z_fak_g * q_g * x_g
          z_local(i,j,k) = z_local(i,j,k)                            &
            &       + z_fak_h * q_h * x_h
          z_local(i,j,k) = z_local(i,j,k)                            &
            &       + z_fak_s * q_s * x_s
          z_local(i,j,k) = z_local(i,j,k)                            &
            &       + z_fak_i * q_i * x_i
        ENDIF

      END DO
    END DO
  END DO

  z_local = 10.0_wp / LOG(10.0_wp) * LOG(z_local * 1.E18_wp + eps)
  !WRITE (6,'(2(A,F10.1))') "  MAX dBZ = ",MAXVAL(z_local),"  MIN dBZ = ",MINVAL(z_local)

  IF ( PRESENT(z_radar) ) THEN
    DO k = 1, iz
      DO j = 1, iy
        DO i = 1, ix
          z_radar(i,j,k) = z_local(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF ( PRESENT(z_radar_dp) ) THEN
    DO k = 1, iz
      DO j = 1, iy
        DO i = 1, ix
          z_radar_dp(i,j,k) = REAL(z_local(i,j,k),kind=dp)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF ( PRESENT(z_radar_850) ) THEN
    DO j = 1, iy
      DO i = 1, ix
        z_radar_850(i,j) = z_local(i,j,klv850)
      ENDDO
    ENDDO
  ENDIF

  IF ( PRESENT(z_radar_cmax) ) THEN
    DO j = 1, iy
      DO i = 1, ix
        z_radar_cmax(i,j) = MAXVAL(z_local(i,j,:))
      ENDDO
    ENDDO
  ENDIF

END SUBROUTINE radar_sb_ray
#endif

!===========================================================================================
!===========================================================================================


SUBROUTINE lightning_potential_index( lpi, buo, wthresh, ie, je, ke, hhl, w, t, p, &
     qv, qc, qr, qi, qs, qg, t0_melt, r_d, r_v, cp_d, cp_v, lh_v, c_l, p0ref, b1, b2w, b3, b4w)

  !===========================================================================================
  !-------------------------------------------------------------------------------------------
  !
  !  DESCRIPTION:
  !
  !  Subroutine for computing the Lightning Potential Index (LPI) after
  !  Lynn et al. (2010) with the possibility of an additional
  !  filtering of stable conditions to inhibit spurious flash signals from deep orographic
  !  wave clouds. The LPI is a 2D field whose values exhibit some correlation to
  !  the flash rate and which is based on the coexistence of updrafts, supercooled water,
  !  graupel/hail and other ice species (cloud ice, snow) in the same grid box
  !  in a temperature range of 0 to -20 deg C. It is basically an integral over
  !  the square of the grid scale vertical speed over the height interval between
  !  0 to -20 deg C, weighted by a function depending on the coexistence of graupel/hail
  !  with other ice species.
  !
  !  There is another subroutine lpi_spatial_filter() to apply a filtering based on
  !  two neighbourhood criteria (column max. vertical speed and total buoyancy
  !  as a stability measure, see comment below) to the LPI. This routine has to be called WITH the
  !  LPI field for the total model domain, not the processor subdomain. Comments see
  !  below in the header of this subroutine.
  !
  !  The present subroutine computes also the total buoyancy (2D-field "buo") as a suitable column-based stability
  !  measure for filtering the LPI with regard to graupel regions of intense orographic wave
  !  related clouds, where lightning activity does not develop. It is basically the vertical integral
  !  over the (negative and positive) buoyancy between ps-50 hPa and ps-550 hPa:
  !
  !      buo = r_d * \int_{ps-50 hPa}^{ps-550 hPa} (Tv_parcel - Tv_surroundings) d ln(p)
  !
  !    Tv_parcel, the virtual temperature of the parcel, is computed from a pseudoadiabatic
  !    ascent of a mixed-layer parcel (average pot. T and qv over lowest 100 hPa) starting
  !    from 50 hPa above ground level. During this ascent, the reversible equivalent potential temperature
  !    (all water substance remains in the parcel) is assumed to be conserved (Emanuel, 1994).
  !
  !  We use buo instead of CAPE, because CAPE is by definiton 0 in the column of an upright
  !  convective updraft core and only counts the positive buoyancy parts of a parcel ascent.
  !  It is only a measure of instability, not stability. buo also takes into account the negative
  !  buoyancy contributions and can better distinguish between a convective updraft (buo ~ 0 J/kg)
  !  and orographic clouds in stable stratification (buo << 0 J/kg).
  !
  !  The actual filter step based on buo is not done in the present subroutine (elim_stable_conditions=.false. below), but,
  !  after spatial smoothing (20 x 20 km^2) of this field, in subroutine lpi_spatial_filter().
  !  A threshold of -1500 J/kg has been determined by U. Blahak based on COSMO-DE summer season data of 2014.
  !  If the smoothed buo falls below this threshold, LPI is set to 0.
  !
  !  NOTE: The computation of buoyancy is relatively expensive. Besides this, the pure function calls are
  !        expensive by themselves, even if all the code in the function body is commented out!
  !
  !  Author: Ulrich Blahak, DWD, 11.6.2015
  !
  !  Literature:
  !
  !       - B. Lynn, Y. Yair, 2010: Prediction of lightning flash density with the WRF model,
  !           Adv. Geosci., 23, 11-16
  !
  !       - Y. Yair et al., 2010: Predicting the potential for lightning activity in Mediterranean storms
  !           based on the Weather Research and Forecasting (WRF) model dynamic and microphysical fields,
  !           JGR, 115, D04205, doi:10.1029/2008JD010868
  !
  !       - K. Emanuel, 1994: Atmospheric Convection, Oxford University Press
  !
  !-------------------------------------------------------------------------------------------
  !===========================================================================================


  INTEGER(kind=iintegers), INTENT(in) :: ie, je, ke

  REAL (KIND=wp), INTENT(IN)  ::  wthresh

  REAL (KIND=wp), DIMENSION(ie,je,ke), INTENT(IN)  ::  &
       t, p, qv, qc, qr, qi, qs, qg
  REAL (KIND=wp), DIMENSION(ie,je,ke+1), INTENT(IN)  ::  &
       hhl, w

  REAL (KIND=wp), DIMENSION(ie,je), INTENT(OUT)  ::  lpi, buo

  REAL(kind=wp), INTENT(in) :: &
       t0_melt, &
       r_d, &
       r_v, &
       cp_d, &
       cp_v, &
       lh_v, &
       c_l, &
       p0ref, &
       b1, &
       b2w, &
       b3, &
       b4w

  REAL(kind=wp) :: &
       p0ref_recip, rdv, o_m_rdv, rvd_m_o, b234w, g, pthick, kappa_dry, esat, lpiincr
  REAL (KIND=wp)  ::  &
       zml, dz, wml, ql, qf, qtot(ie,ke), epsw
  REAL (KIND=wp), DIMENSION(ie,je)  ::  &
       denom
  REAL (KIND=wp), DIMENSION(ie)  ::  &
       t0, theta0, teq, p0, ps, qv0, phkn, thkn, kappa, tparcel, pmin, tdiff
  LOGICAL, DIMENSION(ie,je) :: lstable
  INTEGER(kind=iintegers) :: i, j, k, counter, error, kbot(ie), ktop(ie)

  ! .. Flag for eliminating LPI when the vertically integrated (negative) buoyancy
  !     is below a certain threshold. Note that this flag decides if the
  !     buoyancy criterion is applied within this subroutine for each column
  !     independently. There is another place in the code (subroutine "lpi_spatial_filter")
  !     where this can be done after some horizontal filtering of the buoyancy,
  !     which makes more sense physically.
  LOGICAL,       PARAMETER :: elim_stable_conditions = .FALSE.
  REAL(kind=wp), PARAMETER :: buo_thresh = -1500.0_wp  ! buoyancy threshold for elim_stable_conditions=.true.

  !====================================================================================
  ! .. Some derived constants which are used below:

  p0ref_recip = 1.0 / p0ref
  rdv         = r_d / r_v
  o_m_rdv     = 1.0_wp - rdv 
  rvd_m_o     = r_v / r_d - 1.0_wp
  b234w       = b2w * (b3 - b4w)
  g           = 9.816_wp
  pthick      = 100e2_wp   ! Pa
  kappa_dry   = r_d / cp_d

  !====================================================================================
  ! .. initializations:

  buo = 0.0_wp
  lpi = 0.0_wp

  !====================================================================================
  ! .. preparations:

  !====================================================================================
  ! .. Total buoyancy to determine stable conditions, where the LPI should be 0
  !     (e.g., in graupel regions of strong orographic waves)
  !    The total buoyancy is defined as the vertical integral over
  !
  !      buo = r_d * \int_{ps-50 hPa}^{ps-550 hPa} (Tv_parcel - Tv_surroundings) d ln(p)
  !
  !    Tv_parcel, the virtual temperature of the parcel, is computed from a pseudoadiabatic
  !    ascent of a mixed-layer parcel (average pot. T and qv over lowest 100 hPa) starting
  !    from 50 hPa above ground level. During this ascent, the reversible equivalent potential temperature
  !    (all water substance remains in the parcel) is assumed to be conserved (Emanuel, 1994).
  !

  ! 1) .. starting values for the virtual parcel, assuming no excess temperature:
  DO j = 1, je
    ! .. determine mixed layer mean temperature and qv (average potential temperature resp. qv):
    theta0 (:) = 0.0_wp
    qv0    (:) = 0.0_wp
    ps     (:) = p(:,j,ke) * EXP(g*0.5_wp*(hhl(:,j,ke)-hhl(:,j,ke+1))/(r_d*t(:,j,ke)))
    p0     (:) = ps - 0.5_wp*pthick ! Parcel starting pressure
    ktop   (:) = ke + 1
    qtot (:,:) = qv(:,j,:) + qc(:,j,:) + qr(:,j,:) + qi(:,j,:) + qs(:,j,:) + qg(:,j,:)
    ! .. Integral mean over p and qv, trapezoidal rule and konstant extrapolation to ps resp. ps - pthick:
    DO k = ke, 1, -1
      DO i = 1, ie
        IF (p(i,j,k) >= ps(i)-pthick) THEN
          ktop(i) = k
        END IF
      END DO
    END DO
    ! .. interior points of the p-integral:
    DO k = ke-1, MINVAL(ktop)+1, -1
      DO i = 1, ie
        IF (k > ktop(i)) THEN
          theta0(i) = theta0(i) + 0.5_wp * t(i,j,k)*EXP(kappa_dry*LOG(p0ref/p(i,j,k))) * &
               (p(i,j,k+1)-p(i,j,k-1))
          qv0(i)    = qv0(i) + 0.5_wp * qtot(i,k) * (p(i,j,k+1)-p(i,j,k-1))
        END IF
      END DO
    END DO
    ! .. boundary- and extrapolated points, division by pthick
    DO i = 1, ie
      theta0(i) = ( theta0(i) + &
           0.5_wp * t(i,j,ke)*EXP(kappa_dry*LOG(p0ref/p(i,j,ke))) * (2.0_wp*ps(i)-p(i,j,ke)-p(i,j,ke-1)) + &
           0.5_wp * t(i,j,ktop(i)) * (p(i,j,ktop(i)+1)+p(i,j,ktop(i))-2.0_wp*(ps(i)-pthick)) ) / pthick
      qv0(i) = ( qv0(i) + &
           0.5_wp * qtot(i,ke) * (2.0_wp*ps(i)-p(i,j,ke)-p(i,j,ke-1)) + &
           0.5_wp * qtot(i,ktop(i)) * (p(i,j,ktop(i)+1)+p(i,j,ktop(i))-2.0_wp*(ps(i)-pthick)) ) / pthick
    END DO
    kappa = (r_d*(1.0_wp-qv0) + r_v*qv0) / (cp_d*(1.0_wp-qv0) + cp_v*qv0)
    t0 = theta0 * EXP(kappa*LOG(p0/p0ref))

    ! 2) .. determine the LCL of the virtual parcel:
    CALL hkn_vec(ie, t0, p0, qv0, phkn, thkn, counter, error)

    ! 3) .. minimum pressure for computation of buoyancy (upper bound of integration):
    pmin = MAX(p0-500e2_wp, 250e2_wp)

    ! 4) .. compute buoyancy integral for the virtual parcel (trapezoidal rule)
    !       between p0 and p0 - 500 hPa:

    !    4a) .. vertical index bounds for integration:
    ktop = ke
    kbot = ke
    DO k = ke, 1, -1
      DO i = 1, ie
        IF (p(i,j,k) > p0(i)) THEN
          kbot(i) = k - 1
        END IF
        IF (p(i,j,k) > pmin(i) ) THEN
          ktop(i) = k
        END IF
      END DO
    END DO
    ktop = MIN(ktop, kbot)

    !    4b) .. Reversible equivalent potential temperature at the LCL:
    !            (This temperature is conserved during the moist ascent above LCL)
    CALL Oequiv_vec(ie, thkn, phkn, qv0, teq, error)

    ! 5) Vertical integration of total buoyancy:
    DO k = ke, 1, -1
      IF (k >= MINVAL(ktop)) THEN
        ! .. Temperature of virt. parcel as function of equiv. pot. temp. and pressure:
        CALL Tequiv_vec(ie, teq, p(:,j,k), qv0, tparcel, counter, error)
        IF (error == 0) THEN 
!dir$ safe_address
          DO i = 1, ie
            IF (p(i,j,k) >= phkn(i)) THEN
              ! .. reset T_parcel to dry-adiabatic value below LCL and compute T-difference:
              tparcel(i) = t0(i) * EXP(kappa(i)*LOG(p(i,j,k)/p0(i)))
              tdiff(i) = tparcel(i)*(1.0_wp+rvd_m_o*qv0(i)) - t(i,j,k)*(1.0_wp+rvd_m_o*qv(i,j,k))
            ELSE
              ! .. T-difference with T_parcel from equiv. pot. temp. above LCL:
              esat  = esat_w(t(i,j,k))
              tdiff(i) = tparcel(i)*(1.0_wp+rvd_m_o*spezfeuchte(esat,p(i,j,k))) - t(i,j,k)*(1.0_wp+rvd_m_o*qv(i,j,k))
            END IF
          END DO
!dir$ safe_address
          DO i = 1, ie
            ! Integration of the tdiff-points (p0,tdiff=0), (p(kbot),tdiff(kbot)), ... , (p(ktop),tdiff(ktop)), (pmin,tdiff(ktop))
            !  with trapezoidal rule:
            IF (k > ktop(i) .AND. k < kbot(i)) THEN
              ! full intervals:
              buo(i,j) = buo(i,j) + r_d * 0.5_wp * tdiff(i) / p(i,j,k) * ( p(i,j,k+1) - p(i,j,k-1) )
            ELSE IF (k == ktop(i)) THEN
              ! upper margin interval:
              buo(i,j) = buo(i,j) + r_d * 0.5_wp * tdiff(i) / p(i,j,k) * ( p(i,j,k+1) - p(i,j,k  ) )
              ! constant extrapolation of tdiff until pmin:
              buo(i,j) = buo(i,j) + r_d *          tdiff(i) / p(i,j,k) * ( p(i,j,k  ) - pmin(i)    )
            ELSE IF (k == kbot(i)) THEN
              ! lowest full interval (extra code because of p0 as integration limit):
              buo(i,j) = buo(i,j) + r_d * 0.5_wp * tdiff(i) / p(i,j,k) * ( p0(i)      - p(i,j,k-1) )
              ! The contribution of the lower margin interval is 0 because tdiff0 is 0. Therefore we can leave it out!
              !  buo(i,j) = buo(i,j) + r_d * 0.5_wp * tdiff0 / p0(i)   * ( p0(i)   - p(i,j,k  ) )
            END IF
          END DO
        END IF
      END IF
    END DO
  END DO

  ! .. decide if the filtering of stable points is done here or later as part
  !     of the spatial filtering. The latter has the advantage that the
  !     buoyancy itself can itself be smoothed in space.
  IF (elim_stable_conditions) THEN

    ! .. Assume stable conditions if the total buoyancy is below -1500 J/kg.
    !     We choose a negative threshold because near/within the pseudoadiabatic cores of
    !     grid scale convection, the CAPE can be slightly negative:
    DO j = 1, je
      DO i = 1, ie
        lstable (i,j) = buo(i,j) < buo_thresh
      END DO
    END DO

  ELSE
    lstable = .FALSE.
  END IF

  !====================================================================================
  ! .. Computation of LPI:

  denom = 0.0_wp
  DO k = 1, ke
!dir$ safe_address
    DO j = 1, je
      DO i = 1, ie

        zml = 0.5_wp * (hhl(i,j,k) + hhl(i,j,k+1))
        dz  =           hhl(i,j,k) - hhl(i,j,k+1)
        wml = 0.5_wp * (w  (i,j,k) + w  (i,j,k+1))
        ql = qc(i,j,k) + qr(i,j,k)
        qf = qg(i,j,k) * &
             ( SQRT(qi(i,j,k)*qg(i,j,k))/MAX(qi(i,j,k)+qg(i,j,k), 1e-20_wp) + &
             SQRT(qs(i,j,k)*qg(i,j,k))/MAX(qs(i,j,k)+qg(i,j,k), 1e-20_wp) ) 
        epsw = 2.0_wp * SQRT(ql*qf) / MAX(ql+qf, 1e-20_wp)

        IF (wml >= wthresh) THEN
          lpiincr = wml*wml * epsw * dz
        ELSE
          lpiincr = 0.0_wp
        END IF
        IF (t(i,j,k) <= 273.16_wp .AND. t(i,j,k) >= 253.16_wp) THEN
          lpi(i,j) = lpi(i,j) + lpiincr
          denom(i,j) = denom(i,j) + dz
        END IF

      END DO
    END DO
  END DO

  DO j = 1, je
    DO i = 1, ie
      IF (denom(i,j) >= 1e-20_wp .AND. .NOT.lstable(i,j)) THEN
        lpi(i,j) = lpi(i,j) / denom(i,j)
      ELSE
        lpi(i,j) = 0.0_wp
      END IF
    END DO
  END DO

  !====================================================================================
  !====================================================================================

CONTAINS

  FUNCTION spezfeuchte_vec(e,p,n) RESULT (qv)
    INTEGER(kind=iintegers), INTENT(in) :: n
    REAL(kind=wp),     INTENT(in)  :: e(n), p(n)
    REAL(kind=wp)                  :: qv(n), tmpe(n)

    ! spezfeuchte(e,p)
    ! Spec. humidity in kg/kg as function of vapor pressure e in Pa and
    ! total pressure p in Pa.

    tmpe = MIN(e, p)
    qv = rdv * e / (p-o_m_rdv*e)
  END FUNCTION spezfeuchte_vec

  FUNCTION spezfeuchte(e,p) RESULT (qv)
    REAL(kind=wp),     INTENT(in)  :: e, p
    REAL(kind=wp)                  :: qv, tmpe
    tmpe = MIN(e, p)   
    qv = rdv * e / (p-o_m_rdv*e)
  END FUNCTION spezfeuchte

  FUNCTION esat_w_vec(temp,n) RESULT(es)
    INTEGER(kind=iintegers), INTENT(in) :: n
    REAL(kind=wp),     INTENT(in)  :: temp(n)
    REAL(kind=wp)                  :: es(n)
    es = b1 * EXP( b2w * (temp-b3)/(temp-b4w) )
  END FUNCTION esat_w_vec

  FUNCTION esat_w(temp) RESULT(es)
    REAL(kind=wp),     INTENT(in)  :: temp
    REAL(kind=wp)                  :: es
    es = b1 * EXP( b2w * (temp-b3)/(temp-b4w) )
  END FUNCTION esat_w

  ! Derivative of saturation vapor pressure w.r.t. plain water surface (scalar version):
  FUNCTION desat_w_dT_vec(temp,n) RESULT(des)
    INTEGER(kind=iintegers), INTENT(in) :: n
    REAL(kind=wp),     INTENT(in) :: temp(n)
    REAL(kind=wp)                 :: des(n)
    des = b234w/(temp-b4w)**2 * b1 * EXP( b2w * (temp-b3)/(temp-b4w) )
  END FUNCTION desat_w_dT_vec

  FUNCTION desat_w_dT(temp) RESULT(des)
    REAL(kind=wp),     INTENT(in) :: temp
    REAL(kind=wp)                 :: des
    des = b234w/(temp-b4w)**2 * b1 * EXP( b2w * (temp-b3)/(temp-b4w) )
  END FUNCTION desat_w_dT

  ! The height of the lifting condensation level of a parcel as function of
  ! its intial temperature, pressure and moisture:

  SUBROUTINE hkn_vec(n, temp0, p0, qv0, phkn, thkn, counter, error)

    INTEGER(kind=iintegers), INTENT(in) :: n
    REAL(kind=wp), INTENT(in),  DIMENSION(n) :: temp0, p0, qv0 
    REAL(kind=wp), INTENT(out), DIMENSION(n) :: phkn, thkn
    INTEGER(kind=iintegers), INTENT(out)          :: counter, error

    REAL(kind=wp), DIMENSION(n) :: t_i, p_i, kappa, konst, konst2, f, df, pdiff, term1
    REAL(kind=wp) :: epsilon
    INTEGER(kind=iintegers) :: maxiter, i

    error = 0

    kappa = (r_d*(1.0-qv0)+r_v*qv0) / (cp_d*(1.0-qv0)+cp_v*qv0)

    epsilon = 1e-2_wp
    maxiter = 20
    counter = 0
    p_i = p0
    pdiff = 2.0_wp*epsilon
    konst = temp0/p0**kappa
    konst2 = o_m_rdv*qv0 + rdv

    DO WHILE (MAXVAL(ABS(pdiff)) >= epsilon .AND. counter < maxiter)
      counter = counter + 1
      DO i=1, n
        IF (ABS(pdiff(i)) >= epsilon) THEN
          t_i(i) = konst(i) * EXP(kappa(i)*LOG(p_i(i)))
          term1(i) = konst2(i) * esat_w(t_i(i))
          f(i) = qv0(i) - term1(i) / p_i(i)
          df(i) = (term1(i) - konst2(i) * desat_w_dT(t_i(i)) * kappa(i) * t_i(i)) / (p_i(i)*p_i(i))
          pdiff(i) = - f(i) / df(i)
          p_i(i) = p_i(i) + pdiff(i)
        END IF
      END DO
    END DO

    DO i=1, n
      IF (ABS(pdiff(i)) >= epsilon) THEN
        WRITE (*,*) 'WARNING: Computation of LCL with Newton-method did not converge (LPI)!'
        phkn(i) = -999.99_wp
        thkn(i) = -999.99_wp
        error = 1
      ELSE
        phkn(i) = p_i(i)
        thkn(i) = konst(i) * EXP(kappa(i)*LOG(phkn(i)))
      END IF
    END DO

  END SUBROUTINE hkn_vec

  SUBROUTINE Oequiv_vec(n, temp, druck, qgesamt, teq, error)

    INTEGER(kind=iintegers), INTENT(in) :: n
    REAL(kind=wp), INTENT(in),  DIMENSION(n) :: temp, druck, qgesamt
    REAL(kind=wp), INTENT(out), DIMENSION(n) :: teq
    INTEGER(kind=iintegers), INTENT(out) :: error

    REAL(kind=wp), DIMENSION(n) :: rgesamt, e, r, O, tmperror, tmptemp, tmpqgesamt, pdiff

    ! teq = Oequiv(temp, pressure, qtotal)
    ! Equivalent potential temperature in K after Emanuel (1994), S. 120. 
    ! Temperature in K, pressure in Pa, qtotal in kg/kg (Spec. humidity
    ! of the total water in the air parcel, does not change during
    ! an adiabatic ascent.
    !
    ! Oequiv is a conserved quantity under reversible adiabatic state changes
    ! of a saturated air parcel. All the condensed water remains in the parcel
    ! and accumulates during an ascent. Additional latent heat of freezing
    ! below the freezing level is not taken into account.
    !
    ! This formulation can only be used above the condensation level in a
    ! saturated air parcel, because the moisture is fixed at its saturation value. 
    !
    ! Vectorized version, because input- and output fields are vectors.


    WHERE (temp < 100.0_wp .OR. temp > 400.0_wp .OR. qgesamt >= 0.2_wp .OR. qgesamt < 0.0_wp)
      tmperror = 1
      tmptemp = MIN(MAX(temp, 100.0_wp), 400.0_wp)
      tmpqgesamt = MIN(MAX(qgesamt, 0.0_wp), 0.2_wp)
    ELSEWHERE
      tmperror = 0
      tmptemp = temp
      tmpqgesamt = qgesamt
    END WHERE
   
    error = MAXVAL(tmperror)

    rgesamt = tmpqgesamt / (1.0_wp-tmpqgesamt)
    e = MIN(esat_w_vec(temp, n), MAX(druck,1e2_wp)-0.1e2_wp)
    pdiff = MAX(druck,1e2_wp) - e
    r = spezfeuchte_vec(e,druck, n)      ! spec. saturation moisture
    r = r / (1.0_wp-r)                   ! Saturation mixing ratio
    O = tmptemp * EXP(r_d/(cp_d+c_l*rgesamt)*LOG(p0ref/pdiff))
    teq = O * EXP((lh_v*r)/(tmptemp*(cp_d+c_l*rgesamt)))

  END SUBROUTINE Oequiv_vec

  SUBROUTINE Oequiv(temp, druck, qgesamt, teq, error)

    REAL(kind=wp), INTENT(in)  :: temp, druck, qgesamt
    REAL(kind=wp), INTENT(out) :: teq
    INTEGER(kind=iintegers), INTENT(out) :: error

    REAL(kind=wp) :: rgesamt, e, r, O, tmptemp, tmpqgesamt, pdiff

    IF (temp < 100.0_wp .OR. temp > 400.0_wp .OR. qgesamt >= 0.2_wp .OR. qgesamt < 0.0_wp) THEN
      error = 1
      tmptemp = MIN(MAX(temp, 100.0_wp), 400.0_wp)
      tmpqgesamt = MIN(MAX(qgesamt, 0.0_wp), 0.2_wp)
    ELSE
      error = 0
      tmptemp = temp
      tmpqgesamt = qgesamt
    END IF

    rgesamt = tmpqgesamt / (1.0_wp-tmpqgesamt)
    e = MIN(esat_w(tmptemp), MAX(druck,1e2_wp)-0.1e2_wp)
    pdiff = MAX(druck,1e2_wp) - e
    r = spezfeuchte(e,druck)       ! spec. saturation moisture
    r = r / (1.0_wp-r)             ! Saturation mixing ratio
    O = tmptemp * EXP(r_d/(cp_d+c_l*rgesamt)*LOG(p0ref/pdiff))
    teq = O * EXP((lh_v*r)/(tmptemp*(cp_d+c_l*rgesamt)))

  END SUBROUTINE Oequiv

  SUBROUTINE Tequiv_vec(n, teq, druck, qgesamt, temp, counter, error)

    INTEGER(kind=iintegers), INTENT(in) :: n
    REAL(kind=wp), INTENT(in), DIMENSION(n)  :: teq, druck, qgesamt
    REAL(kind=wp), INTENT(out), DIMENSION(n) :: temp
    INTEGER(kind=iintegers), INTENT(out) :: counter, error

    REAL(kind=wp), DIMENSION(n) :: tneu, f, df, tdiff, Oe, OepdT
    REAL(kind=wp)               :: dT_newton, epsilon
    INTEGER(kind=iintegers) :: maxiter, i, zierr1, zierr2

    ! temp = Tequiv (teq, pressure, qtotal)
    ! Tempertature temp in K as function of the equivalent potential temperature teq in K and
    ! pressure in Pa and the total specific water content of the 
    ! air parcel under the assumption that the parcel is saturated and no condensed
    ! water leaves the parcel.
    !
    ! This formulation can only be used above the condensation level in a
    ! saturated air parcel, because the moisture is fixed at its saturation value. 
    !
    ! Vectorized version, because input- and output fields are vectors.

    ! Iterative solution. We make sure that the converged values of temp do
    ! not depend on the number of iterations necessary for other vector elements:

    error = 0

    epsilon = 1e-2_wp
    maxiter = 30
    dT_newton = 0.5_wp  ! K

    tdiff = 2.0_wp*epsilon
    tneu = teq * 0.9_wp
    counter = 0

    iteration: DO WHILE (MAXVAL(ABS(tdiff)) >= epsilon .AND. counter < maxiter)
      counter = counter + 1
      DO i=1, n
        IF (ABS(tdiff(i)) >= epsilon) THEN
          CALL Oequiv(tneu(i), druck(i), qgesamt(i), Oe(i), zierr1)
          CALL Oequiv(tneu(i)+dT_newton, druck(i), qgesamt(i), OepdT(i), zierr2)
          IF (zierr1 /= 0 .OR. zierr2 /= 0) THEN
            error = 2
            exit iteration
          END IF
          f(i) = Oe(i) - teq(i)
          df(i)  = (OepdT(i) - Oe(i)) / dT_newton
          tdiff(i) = - f(i) / df(i)
          tneu(i) = tneu(i) + tdiff(i)
        END IF
      END DO
    END DO iteration

    IF (error == 2) THEN
      WRITE (*,'(a)') 'WARNING: Computation of T from Tequiv with Newton-method' // &
           ' did not converge (LPI)! Unrealistic atmospheric state encountered'
      RETURN
    END IF

    DO i=1, n
      IF (ABS(tdiff(i)) >= epsilon) THEN
        WRITE (*,'(a,2(f12.1))') 'WARNING: Computation of T from Tequiv with Newton-method' // &
             ' did not converge (LPI)! min/max p = ', MINVAL(druck), MAXVAL(druck)
        temp(i) = -999.99_wp
        error = 1
      ELSE
        temp(i) = tneu(i)
      END IF
    END DO

  END SUBROUTINE Tequiv_vec

  SUBROUTINE Tequiv(teq, druck, qgesamt, temp, counter, error)

    REAL(kind=wp), INTENT(in)  :: teq, druck, qgesamt
    REAL(kind=wp), INTENT(out) :: temp
    INTEGER(kind=iintegers), INTENT(out) :: counter, error

    REAL(kind=wp) :: tneu, epsilon, f, df, tdiff, Oe, OepdT, dT_newton
    INTEGER(kind=iintegers) :: maxiter, zierr1, zierr2

    ! temp = Tequiv (teq, pressure, qtotal)
    ! Tempertature temp in K as function of the equivalent potential temperature teq in K and
    ! pressure in Pa and the total specific water content of the 
    ! air parcel under the assumption that the parcel is saturated and no condensed
    ! water leaves the parcel.
    !
    ! This formulation can only be used above the condensation level in a
    ! saturated air parcel, because the moisture is fixed at its saturation value. 

    error = 0

    epsilon = 1e-2_wp
    maxiter = 30
    dT_newton = 0.5_wp  ! K

    tdiff = 2.0_wp*epsilon
    tneu = teq * 0.9_wp
    counter = 0

    iteration: DO WHILE (ABS(tdiff) >= epsilon .AND. counter < maxiter)
      counter = counter + 1
      CALL Oequiv(tneu, druck, qgesamt, Oe, zierr1)
      CALL Oequiv(tneu+dT_newton, druck, qgesamt, OepdT, zierr2)
      IF (zierr1 /= 0 .OR. zierr2 /= 0) THEN
        error = 2
        EXIT iteration
      END IF
      f = Oe - teq
      df  = (OepdT - Oe) / dT_newton
      tdiff = - f / df
      tneu = tneu + tdiff
    END DO iteration

    IF (error == 2) THEN
      WRITE (*,'(a)') 'WARNING: Computation of T from Tequiv with Newton-method' // &
           ' did not converge (LPI)! Unrealistic atmospheric state encountered'
      RETURN
    END IF

    IF (abs(tdiff) >= epsilon) THEN
      WRITE (*,'(a,f12.1)') 'WARNING: Computation of T from Tequiv with Newton-method did not converge (LPI)! p = ', druck
      temp = -999.99_wp
      error = 1
    ELSE
      temp = tneu
    END IF

  END SUBROUTINE Tequiv

END SUBROUTINE lightning_potential_index

!====================================================================================

SUBROUTINE lpi_spatial_filter ( lpi_full, wmax, wthresh, buo, ie_tot, je_tot, nbl, &
                                dlon, dlat, startlat_tot, r_earth )

  !====================================================================================
  !-----------------------------------------------------------------------------------
  !
  ! Spatial filter for the lightning potential index:
  !
  !  - Majority of grid columns in a neighbourhood (+-neigh_wmax) must have wmax >= wthresh.
  !    A good value for COSMO-DE is wthresh = 1.1 m/s, as determined by U. Blahak
  !    based on summer season data from 2014.
  !
  !  - OPTIONAL by internal switch "elim_stable_conditions" below:
  !    The total buoyancy (between ps-50hPa and ps-550 hPa) must be >= -1500 J/Kg.
  !    This threshold has been determined by U. Blahak based on summer season data from 2014
  !    of COSMO-DE and lightning observations.
  !    Further, the input buoyancy field can optionally be smoothed over a certain
  !    neighbourhood region (+-neigh_buo) to obtain spatially more stable results.
  !    (internal switch "filter_buoyancy" below)
  !
  ! This routine has to be called with the LPI- and total buoyancy field
  !  for the total model domain, not the processor subdomain!
  !
  ! NOTE: - Filtering does not work correctly in case of periodic BCs!
  !       - If the interior domain is smaller than the filtering/smoothing
  !         neighbourhood(s), the neighbourhoods are reduced accordingly.
  !
  ! Author: Ulrich Blahak, DWD, 11.6.2015
  !
  ! Literature: The wmax-based filtering has been proposed by the original
  !             author(s) of the LPI. The total buoyancy criterion has been developed
  !             by U. Blahak.
  !
  !-----------------------------------------------------------------------------------
  !====================================================================================
  
  ! Parameter list:

  REAL    (KIND=sp),        INTENT (INOUT) :: lpi_full(ie_tot,je_tot)
  INTEGER (KIND=iintegers), INTENT (IN)    :: ie_tot, je_tot, nbl
  REAL    (KIND=wp),        INTENT (IN)    :: wthresh, dlon, dlat, startlat_tot, r_earth
  REAL    (KIND=wp),        INTENT (IN)    :: wmax(ie_tot,je_tot), buo(ie_tot,je_tot)

  ! Local Variables
  LOGICAL,           PARAMETER   :: elim_stable_conditions = .TRUE.
  LOGICAL,           PARAMETER   :: filter_buoyancy        = .TRUE.
  REAL    (kind=wp), PARAMETER   :: buo_thresh = -1500.0_wp ! buoyancy threshold for elim_stable_conditions=.true.
  REAL    (KIND=wp), PARAMETER   :: neigh_wmax = 5.0e3_wp   ! neighbourhood = point +- neigh [m]
  REAL    (KIND=wp), PARAMETER   :: neigh_buo  = 10.0e3_wp  ! neighbourhood = point +- neigh [m]
  REAL    (KIND=wp), PARAMETER   :: pi = 4.0_wp*atan(1.0_wp)
  REAL    (KIND=wp)              :: dx, dy, buo_sum(ie_tot, je_tot)
  INTEGER (KIND=iintegers)       :: i, j, im, jm, iu, ju, nyes(ie_tot, je_tot), neigh_i, neigh_j, neigh_ij

  dx = dlon * r_earth * pi / 180.0_wp * COS((startlat_tot+(je_tot/2)*dlat)*pi/180.0_wp)
  dy = dlat * r_earth * pi / 180.0_wp

  !=======================================================================
  ! Filtering with wmax-criterion, but limit neighourhood to domain size:
  !  (Note: periodic BCs not implemented so far!)

  neigh_i = MIN( NINT(neigh_wmax/dx), (ie_tot-2*nbl-1)/2 )
  neigh_j = MIN( NINT(neigh_wmax/dy), (je_tot-2*nbl-1)/2 )

  neigh_ij = (2*neigh_i+1) * (2*neigh_j+1)

  ! Sum up the number of points in the neighbourhood where wmax exceeds wthresh:
  nyes = 0
  DO j = -neigh_j, neigh_j
    DO i = -neigh_i, neigh_i
      DO jm = 1+nbl, je_tot-nbl
        DO im = 1+nbl, ie_tot-nbl
          
          ! actual working point:
          iu = im + i
          ju = jm + j
          ! points in the outer halo are replaced by inside mirror points:
          !   (does not make sense for periodic BCs!)
          iu = iu + 2*ABS(MIN(iu-nbl,0)) - 2*ABS(MIN(ie_tot-nbl-iu,0))
          ju = ju + 2*ABS(MIN(ju-nbl,0)) - 2*ABS(MIN(je_tot-nbl-ju,0))
          
          IF (wmax(iu,ju) >= wthresh) THEN
            nyes(im,jm) = nyes(im,jm) + 1
          END IF
          
        END DO
      END DO
    END DO
  END DO
  
  DO j = 1, je_tot
    DO i = 1, ie_tot
      IF (REAL(nyes(i,j),wp) < 0.5_wp*neigh_ij ) THEN
        ! NOTE: in the outer halo nyes is 0, so all lpi values are set to 0 there!
        lpi_full(i,j) = REAL(0.0, KIND=sp)
      END IF
    END DO
  END DO
  

  !==============================================================================
  ! Filtering with buoyancy-criterion, but limit smoothing region to domain size:
  !  (Note: periodic BCs not implemented so far!)

  IF (elim_stable_conditions) THEN

    IF (filter_buoyancy) THEN
      neigh_i = MIN( NINT(neigh_buo/dx), (ie_tot-2*nbl-1)/2 )
      neigh_j = MIN( NINT(neigh_buo/dy), (je_tot-2*nbl-1)/2 )
    ELSE
      ! no neighbourhood smoothing:
      neigh_i = 0
      neigh_j = 0
    END IF

    neigh_ij = (2*neigh_i+1) * (2*neigh_j+1)

    ! Sum up the number of points in the neighbourhood where wmax exceeds wthresh:
    buo_sum = 0.0_wp
    DO j = -neigh_j, neigh_j
      DO i = -neigh_i, neigh_i
        DO jm = 1+nbl, je_tot-nbl
          DO im = 1+nbl, ie_tot-nbl

            ! actual working point:
            iu = im + i
            ju = jm + j
            ! points in the outer halo are replaced by inside mirror points:
            iu = iu + 2*ABS(MIN(iu-nbl,0)) - 2*ABS(MIN(ie_tot-nbl-iu,0))
            ju = ju + 2*ABS(MIN(ju-nbl,0)) - 2*ABS(MIN(je_tot-nbl-ju,0))

            buo_sum(im,jm) = buo_sum(im,jm) + buo(iu,ju)

          END DO
        END DO
      END DO
    END DO

    DO j = 1+nbl, je_tot-nbl
      DO i = 1+nbl, ie_tot-nbl
        IF (buo_sum(i,j) < buo_thresh*neigh_ij ) THEN
          ! NOTE: in the outer halo lpi is not modified, because it
          !       has already been set to 0 above by the wmax criterion.
          lpi_full(i,j) = REAL(0.0, KIND=sp)
        END IF
      END DO
    END DO

  END IF

END SUBROUTINE lpi_spatial_filter

!===========================================================================================
!===========================================================================================

SUBROUTINE max_in_pressure_layers( field, pp, p0, ie, je, ke,        &
                                   player_max, pres_bot, pres_top )

INTEGER (KIND=iintegers), INTENT(IN)  ::  &
  ie, je, ke

REAL (KIND=wp), INTENT(IN)  ::  &
  field (ie,je,ke),    &
  pp    (ie,je,ke),    &
  p0    (ie,je,ke)

REAL (KIND=wp), INTENT(IN), OPTIONAL  ::  &
  pres_bot,              &
  pres_top

REAL (KIND=wp), INTENT(out)  ::  &
  player_max(ie,je)

REAL (KIND=wp) ::    &
  pressure(ie,je,ke),    &
  masked_field(ie,je,ke)

REAL (KIND=wp) ::    &
  maskval

!------------------------------------------------------------------------------

  pressure(:,:,:) = p0(:,:,:) + pp(:,:,:)

! JF:   maskval = -1.0E20_wp
  maskval = 0.0_wp

  masked_field(:,:,:) = field(:,:,:)

  IF ( PRESENT( pres_bot ) .AND. pres_bot > 0._wp ) &
    WHERE ( pressure(:,:,:) > pres_bot )  masked_field(:,:,:) = maskval
  IF ( PRESENT( pres_top ) .AND. pres_top > 0._wp ) &
    WHERE ( pressure(:,:,:) <= pres_top ) masked_field(:,:,:) = maskval

  player_max(:,:) = MAXVAL( masked_field(:,:,:), 3 )

END SUBROUTINE max_in_pressure_layers

!===========================================================================================

END MODULE pp_utilities
