!+ Data module for mathematical, physical and parametrizational constants
!-------------------------------------------------------------------------------

MODULE data_constants

!-------------------------------------------------------------------------------
!
! Description:
!  This module contains mathematical and physical constants as well as
!  variables that are used in several parameterization schemes.
!  These variables are constant through a model run, but some of them
!  can be determined via Namelist and are "tuning constants".
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8236 1493
!  email:  uschaettler@dwd.d400.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Ulrich Schaettler
!  Initial release
! 1.20       1999/01/07 Guenther Doms
!  Renaming of some global variables.
! 1.30       1999/06/24 Matthias Raschendorfer
!  Declaration of 5 new constants (lhocp, rcpv, rcpl, con_m, con_h)
! 2.18       2002/07/16 Reinhold Schrodin
!  Eliminated variable rhde (will be replaced by cf_snow from data_soil)
! 3.16       2005/07/22 Ulrich Schaettler
!  Added dielectric constants for water and ice; added density of ice
! 3.21       2006/12/04 Ulrich Schaettler
!  Added new Namelist parameters for group /TUNING/: clc_diag, q_crit, qi0, qc0
! V3_23        2007/03/30 Matthias Raschendorfer
!  Moving 'clc_diag', 'q_crit' and 'akt' to MODULE data_turbulence.
! V4_12        2010/05/11 Ulrich Schaettler
!  Renamed t0 to t0_melt because of conflicting names
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_20        2011/08/31 Matthias Raschendorfer
!  Introduction of the reference pressure p0ref used in the Exner function
! V4_27        2013/03/19 Michael Baldauf
!  Moved SR set_constants from src_setup to this module, so that it can
!   also be used easily by other programs
! V5_1         2014-11-28 Michael Baldauf, Oliver Fuhrer
!  Introduction of a new variable p0ref_recip (for 1/p0ref)
!  Replaced ireals by wp (working precision) (OF)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!
! Declarations:
!
! Modules used:
USE data_parameters , ONLY :   &
  wp,        & ! KIND-type parameters for real variables
  iintegers    ! KIND-type parameter for standard integer variables

!=======================================================================

IMPLICIT NONE

!=======================================================================
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
! Global (i.e. public) Declarations:

! 1. mathematical constants
! -------------------------

  REAL  (KIND=wp)                  ::           &
    pi              ! circle constant

! 2. physical constants and related variables
! -------------------------------------------

  REAL  (KIND=wp)                  ::           &
    t0_melt,      & ! melting temperature of ice
    r_d,          & ! gas constant for dry air
    r_v,          & ! gas constant for water vapor
    rdv,          & ! r_d / r_v
    o_m_rdv,      & ! 1 - r_d/r_v
    rvd_m_o,      & ! r_v/r_d - 1
    cp_d,         & ! specific heat of dry air at constant pressure
    cpdr,         & ! 1 / cp_d
    rdocp,        & ! r_d / cp_d
    gamma,        & ! 1 / (1 - rdocp)   ( = cp_d/cv_d)
    lh_v,         & ! latent heat of vapourization
    lh_f,         & ! latent heat of fusion
    lh_s,         & ! latent heat of sublimation
    lhocp,        & ! lh_v / cp_d
    rcpv,         & ! cp_d / cp_v - 1
    rcpl,         & ! cp_d / cp_l - 1
    con_m,        & ! kinematic viscosity (m2/s)
    con_h,        & ! scalar conductivity (m2/s)
    g,            & ! acceleration due to gravity
    gq,           & ! g * g
    gh,           & ! g / 2
    gr,           & ! 1 / g
    r_earth,      & ! mean radius of the earth (m)
    day_len,      & ! mean length of the day (s)
    rho_w,        & ! density of liquid water (kg/m^3)
    rho_ice,      & ! density of ice          (kg/m^3)
    K_w,          & ! dielectric constant for water
    K_ice,        & ! dielectric constant for ice
    sigma,        & ! Boltzmann-constant
    solc            ! solar constant

! 3. constants for parametrizations
! ---------------------------------

  REAL  (KIND=wp)                  ::           &
    p0ref,        & ! reference pressure for Exner-function (Pa)
    p0ref_recip,  & ! reciprocal of p0ref                   (Pa^-1)
    b1,           & ! variables for computing the saturation vapour pressure
    b2w,          & ! over water (w) and ice (i) according to Teten's formula
    b2i,          & !               -- " --
    b3,           & !               -- " --
    b4w,          & !               -- " --
    b4i,          & !               -- " --
    b234w,        & ! b2w * (b3 - b4w)
    uc1,          & ! variable for computing the rate of cloud cover in 
    uc2,          & ! the unsaturated case
    ucl,          & !               -- " --
    aks2,         & ! variable for horizontal diffusion of second order
    aks4            ! variable for horizontal diffusion of fourth order

! 4. tuning constants for radiation, cloud physics, turbulence
! ------------------------------------------------------------

  REAL  (KIND=wp)                  ::           &
    qi0,          & ! cloud ice threshold for autoconversion
    qc0             ! cloud water threshold for autoconversion


! 5. Precision-dependent security parameters (epsilons)
! ------------------------------------------------------

  REAL (KIND=wp), PARAMETER  ::                                      &

    repsilon = 1.0E8_wp*TINY(1.0_wp),                                &

      ! Very small number near zero.
      ! To be used mainly to avoid division by zero, e.g.
      ! eps_div = repsilon ; x = y / MAX(z,eps_div) ! for z > 0.
      ! Note that the factor 1.0E8 has been chosen rather
      ! arbitrarily to get some distance to zero.

    rprecision = 10.0_wp**(-PRECISION(1.0_wp))

      ! Precision of 1.0 in additions/subtractions.
      ! To be used for instance to check equality of reals, e.g.
      ! eps_cmpr = rprecision ; IF (ABS(a-b) < eps_cmpr) equal=.true.,
      ! or to increase the magnitude of an epsilon only in SP, e.g.
      ! epsilon = MAX(1.0E-8_wp,rprecision).

      !
      ! Approximate magnitudes:
      !
      !       | repsilon | rprecision
      !   ----+----------+------------
      !    SP |  1.0E-30 |   1.0E-7
      !   ----+----------+------------
      !    DP | 1.0E-300 +  1.0E-16
      !

!=======================================================================

CONTAINS

!==============================================================================
!==============================================================================
!+ Subroutine for initialization of constants used throughout the program
!------------------------------------------------------------------------------

SUBROUTINE set_constants

!------------------------------------------------------------------------------
!
! Description:
!   This routine initializes physical and mathematical constants used 
!   in the LM
!
! Method:
!   Arithmetical statements
!
!------------------------------------------------------------------------------

!- End of header
!==============================================================================

! Variables from data_constants

! mathematical constants
! ----------------------

  pi       = 4.0_wp * ATAN (1.0_wp)

! physical constants and related variables
! ----------------------------------------

  t0_melt  =   273.15_wp
  r_d      =   287.05_wp
  r_v      =   461.51_wp
  rdv      =         r_d / r_v
  o_m_rdv  =         1.0_wp - rdv 
  rvd_m_o  =         r_v/r_d - 1.0_wp
  cp_d     =  1005.0_wp
  cpdr     =         1.0_wp / cp_d
  rdocp    =         r_d / cp_d
  gamma    =         1.0_wp / (1.0_wp - rdocp)
  lh_v     =     2.501E6_wp
  lh_f     =     0.334E6_wp
  lh_s     =     2.835E6_wp
  lhocp    =         lh_v / cp_d
  rcpv     =  0.8396_wp
  rcpl     =  3.1733_wp
  con_m    =     1.50E-5_wp
  con_h    =     2.20E-5_wp
  g        =     9.80665_wp
  gq       =         g * g
  gh       =         g * 0.5_wp
  gr       =         1.0_wp / g
#ifndef MESSY
! NOTE: Unfortunately the radius of the Earth is not consistent between
!       ECHAM5 and COSMO. It needs to be tested, how the models respond
!       to a change of this parameter ... 
  r_earth  =  6371.229E3_wp
#else
  r_earth  =  6371000.0_wp
#endif
  day_len  = 86164.09054_wp
  rho_w    =  1000.0_wp
  rho_ice  =   900.0_wp
  K_w      =     0.930_wp
  K_ice    =     0.176_wp
  sigma    =     5.6697E-8_wp
  solc     =  1368.0_wp

! constants for parametrizations
! ------------------------------

  p0ref    =   1.0E5_wp
  p0ref_recip = 1.0_wp / p0ref
  b1       =   610.78_wp
  b2w      =    17.2693882_wp
  b2i      =    21.8745584_wp
  b3       =   273.16_wp
  b4w      =    35.86_wp
  b4i      =     7.66_wp
  b234w    =      b2w * (b3 - b4w)
  uc1      =     0.8_wp
  uc2      =  SQRT(3.0_wp)
  ucl      =     1.00_wp

END SUBROUTINE set_constants

END MODULE data_constants
