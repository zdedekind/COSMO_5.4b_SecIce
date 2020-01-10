!+ Source Module providing routines for numerical schemes used in "dynamics".
!------------------------------------------------------------------------------

MODULE numeric_utilities

!------------------------------------------------------------------------------
!
! Description:
!   This module provides some service utilities related to various numerical
!   schemes used in the dynamics.
!     - no routine uses other modules, except the declarations for the
!       KIND-type parameter; the data access is by parameter list only
!     - no routine allocates dynamic memory; work space needed is
!       provided via the parameter list
!     - no derived data types are used
!   The routines are written in plug-compatible format. Purists may 
!   identify the input/output variables to start with the letter 'p'.
!
!   Routines (module procedures) currently contained:
!     - hadv_cd2
!     - lap_2
!     - lap_4
!     - lap_4a
!     - lap_4am
!     - lap_4aml
!     - hadv_pd
!     - vadv_pd
!     - hadv_pd_2tl
!     - vadv_pd_2tl
!
!   Routines contained for the Semi-Lagrange advection for the
!   moisture variables
!     - backtraj_trilin_dt1_3tl
!     - backtraj_trilin_dt2_3tl
!     - interpol_sl_trilin
!     - interpol_sl_tricubic
!
!   Routines necessary for computing the supercell detection indices
!     - mean_over_box
!     - mean_cov_over_box
!     - vert_avg
!
!   Routines for metrics
!     - calc_sqrtg_r
!     - metric_coeffs
!     - curl
!     - calc_Theta_Tppp
!
! Current Code Owner: DWD, Michael Baldauf
!  phone:  +49  69  8062 2729
!  fax:    +49  69  8062 3721
!  email:  guenther.doms@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.32       1999/08/24 Guenther Doms
!  Initial release
! 1.34       1999/12/10 Ulrich Schaettler
!  Corrections for CALLs to intrinsic functions
! 1.37       2000/03/24 Guenther Doms
!  New subroutine 'lap_4a' added
! 2.9        2001/07/16 Guenther Doms
!  New subroutine 'lap_4am' for monotonic diffusion with the Xue limiter added.
!  New subroutine 'lap_4aml' for flux-limited diffusion including an 
!  orographic flux limiter added.
! 2.10       2001/07/24 Guenther Doms
!  Subroutine 'gauss' removed as it is no longer required.
!  Corrections for CALLs to intrinsic functions in routine lap_4aml.
! 2.11       2001/09/28 Ulrich Schaettler
!  Correction of a bug in lap_4aml
! 2.14       2002/02/15 Guenther Doms
!  Limitation of the vertical velocity in case of Courant numbers larger
!  than 1 in the routine 'vadv_pd' for monotonic explicit vertical advection.
! 2.18       2002/07/16 Almut Gassmann
!  New routines for horizontal and vertical advection used in the new
!  2 time level scheme.
! 3.7        2004/02/18 Jochen Foerstner + Michael Baldauf
!  New routines and functions for positive definite horizontal and vertical
!  advection and for normal advection for the Runge-Kutta scheme.
!  New routines for the semi-Lagrange advection of the moisture variables
!  Replaced cphi by crlat
! 3.9        2004/04/22 Michael Baldauf
!  Bug corrections in backtraj_trilin_dt1_3tl, backtraj_trilin_dt2_3tl
! 3.13       2004/12/03 Michael Baldauf
!  Correction of lower boundary in interpol_SL_trilin
! 3.14       2005/01/25 Jochen Foerstner
!  Changes in the vertical advection routine for Runge-Kutta scheme
! 3.15       2005/03/03 Michael Baldauf
!  Eliminated Subroutine adv_sl_trilin_dt1_3tl (not longer needed)
! 3.16       2005/07/22 Ulrich Schaettler, Michael Baldauf
!  Moved several routines to new utility-module numeric_utilities_rk
!  Added new subroutine interpol_sl_tricubic
!  Bug correction for computing the backward trajectory for SL advection
! 3.18       2006/03/03 Michael Baldauf
!  Corrected a warning in SL interpolation
! 3.21       2006/12/04 Ulrich Schaettler / Jochen Foerstner
!  Changed interfaces for Laplace-Operator subroutines to introduce
!  hd_mask*dcoeff
! V3_23        2007/03/30 Simone Campagna, Davide Cesari
!  Adapted some loop boundaries in SR lap_2
! V4_1         2007/12/04 Michael Baldauf
!  New routines for calculating meanvalues over boxes
! V4_4         2008/07/16 Ulrich Schaettler
!  Changed order of SR arguments for some routines (first dimensions, then fields)
!  Vectorization of computing area and volume integrals (by JO Beismann, NEC)
! V4_7         2008/12/12 Ulrich Schaettler
!  Remove a WRITE statement from SR interpol_sl_tricubic for vectorization
! V4_8         2009/02/16 Ulrich Schaettler
!  Eliminated use of data_modelconfig, data_fields
! V4_9         2009/07/16 Ulrich Schaettler, Christian Bollmann
!  Implemented 3D versions of routines lap_2, lap_4am, lap_4aml to save CALLs
!  Inserted Compiler Directives
!  Inserted SR calc_sqrtg_r to calculate reciprocal square root of G
! V4_12        2010/05/11 Michael Baldauf
!  Introduced new subroutines metric_coeffs, curl, calc_Theta_Tppp
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_15        2010/11/19 Michael Baldauf
!  new subroutines solve_5banddiag and the 'Numerical Recipes'-routines
!  bandec and banbks in 'vectorized' versions
! V4_18        2011/05/26 Michael Baldauf
!  Added optional "selective filling diffusion" for Semi-Lagrangian advection
!  New SR: remove_negative_values, diffusion, diffusion_x, diffusion_y, diffusion_z
!  Changed interface to SR calc_Theta_Tppp: only use full temperature (not tp+t0)
!   (Jean-Marie Bettems)
! V4_20        2011/08/31 Ulrich Schaettler
!  Corrected dimensions of acrlat, tgrlat in SR curl (tgrlat now 2 dimensions)
!  Print debug output only in case of ldebug_dyn
! V4_23        2012/05/10 Michael Baldauf 
!  Shifted SR calc_sqrtg_r, metric_coeffs to new module grid_metrics_utilities
!  Use fields sqrtg_r_* from grid_metrics_utilities
! V4_27        2013/03/19 Michael Baldauf
!  Modified some error messages for better understanding
! V4_30        2013/11/08 Ulrich Schaettler
!  Removed old CRAY compiler directives, which are not recognized any more
! V5_1         2014-11-28 Michael Baldauf, Oliver Fuhrer
!  Introduction of a new variable p0ref_recip (for 1/p0ref)
!  SR. curl: 
!   => improved weighted vertical interpolation for some z-derivatives.
!   => proper horizontal interpolation for the spherical metric correction
!   => one lateral boundary line is set to 0 (therefore, no uninitialised field elements)
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
!
USE data_parameters, ONLY :   &
    wp        ,& ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

USE data_constants, ONLY:     &
    repsilon     ! precision of 1.0 in current floating point format

!------------------------------------------------------------------------------
 
IMPLICIT NONE

!------------------------------------------------------------------------------

REAL (KIND=wp), PARAMETER :: &
  eps_div = repsilon ! small, precision-dependent value to be used in divisions
                     ! to avoid division by zero, e.g. a/MAX(b,eps_div)

!------------------------------------------------------------------------------

CONTAINS

!==============================================================================
!==============================================================================

SUBROUTINE curl    ( ie, je, ke, eddlon, eddlat, r_earth,    &
  acrlat, tgrlat, sqrtg_r_s,                                 &
  dzeta_dlam, dzeta_dphi, lmetr, wgtfac,                     &
  v1, v2, v3,                                                &
  l_phys_comp, curl1, curl2, curl3 )

!------------------------------------------------------------------------------
!
! Description:
!   calculate the contravariant (physical or 'non-physical') components of
!   curl v in spherical base vectors for given contravariant,
!   physical components of a vector field v (in spherical base vectors)
!
! Main input-variables:
!   the physical, contravariant components v1, v2, v3 of v
!   (typically: v1=u, v2=v, v3=w)
!   grid position: v1 at u-pos., v2 at v-pos., v3 at w-pos.
!   (other input variables, see below)
!
! Output:
!   the 3 contravariant components of curl v
!   grid position: curl1, curl2, curl3 at the scalar position
!   if l_phys_comp=.TRUE (default): physical components are produced
!   if l_phys_comp=.FALSE         : contravariant components are produced
!
! Method:
!   centered finite differences; no calculation at the lateral boundaries
!
!------------------------------------------------------------------------------

INTEGER (KIND=iintegers), INTENT(IN)  ::  &
  ie, je, ke
  
REAL (KIND=wp),     INTENT(IN)  ::  &
  eddlon,     & ! 1 / dlon, dlon            (in degrees)
  eddlat,     & ! 1 / dlat, dlat            (in degrees)
  r_earth       ! mean radius of the earth  (m)

REAL (KIND=wp),     INTENT(IN)  ::  &
  acrlat(je,2)  ! = 1 / ( r_earth * cos(lat) ), lat in rotated coord.

REAL (KIND=wp),     INTENT(IN)  ::  &
  tgrlat(je,2)  ! = tan(lat), lat in rotated coord.

REAL (KIND=wp),     INTENT(IN)  ::  &
  sqrtg_r_s (ie,je,ke), & ! 1 / sqrt(G)       (at the scalar position)
  dzeta_dlam(ie,je,ke), & ! d zeta / d lambda (at the scalar position)
  dzeta_dphi(ie,je,ke)    ! d zeta / d phi    (at the scalar position)

REAL (KIND=wp),     INTENT(IN)  ::  &
  wgtfac(ie,je,ke+1)    ! weighting factor for vertical interpolation ( 1 )

LOGICAL, INTENT(IN) ::  lmetr  ! with metrics of the spherical earth

REAL (KIND=wp),     INTENT(IN)  ::  &
  v1(ie,je,ke),         & ! 
  v2(ie,je,ke),         & !
  v3(ie,je,ke+1)          !

LOGICAL, INTENT(IN)             :: l_phys_comp

REAL (KIND=wp),     INTENT(OUT) ::  &
  curl1(ie,je,ke),      & !
  curl2(ie,je,ke),      & !
  curl3(ie,je,ke)         !

! Local Variables:

REAL (KIND=wp)     :: r_earth_inv
REAL (KIND=wp)     :: dv1_dlam, dv1_dphi, dv1_dzeta
REAL (KIND=wp)     :: dv2_dlam, dv2_dphi, dv2_dzeta
REAL (KIND=wp)     :: dv3_dlam, dv3_dphi, dv3_dzeta

REAL (KIND=wp)     :: v1_at_s, v2_at_s

INTEGER  (KIND=iintegers) :: i, j, k

!------------------------------------------------------------------------------

  r_earth_inv = 1.0_wp / r_earth

  ! first calculate the physical components 
  ! (without influence of spherical shape; see below)
  DO k=2, ke-1
    DO j=2, je-1
      DO i=2, ie-1

        ! derivatives at the scalar position:
        dv1_dlam  =   ( v1(i,j,k)   - v1(i-1,j,k) ) * eddlon
        dv1_dphi  = ( ( v1(i,j+1,k) + v1(i-1,j+1,k) )                   &
          &         - ( v1(i,j-1,k) + v1(i-1,j-1,k) ) ) * 0.25_wp * eddlat
        dv1_dzeta = 0.5_wp * (                                                 &
          &               wgtfac(i,j,k+1)  * ( v1(i,j,k+1) + v1(i-1,j,k+1) )   &
          &   + (1.0_wp - wgtfac(i,j,k+1)                                      &
          &             - wgtfac(i,j,k)  ) * ( v1(i,j,k)   + v1(i-1,j,k)   )   &
          &   - (1.0_wp - wgtfac(i,j,k)  ) * ( v1(i,j,k-1) + v1(i-1,j,k-1) ) )

        dv2_dlam  = ( ( v2(i+1,j,k) + v2(i+1,j-1,k) )                   &
          &         - ( v2(i-1,j,k) + v2(i-1,j-1,k) ) ) * 0.25_wp *eddlon
        dv2_dphi  =   ( v2(i,j,k)   - v2(i,j-1,k) ) * eddlat
        dv2_dzeta = 0.5_wp * (                                                 &
          &               wgtfac(i,j,k+1)  * ( v2(i,j,k+1) + v2(i,j-1,k+1) )   &
          &   + (1.0_wp - wgtfac(i,j,k+1)                                      &
          &             - wgtfac(i,j,k)  ) * ( v2(i,j,k)   + v2(i,j-1,k)   )   &
          &   - (1.0_wp - wgtfac(i,j,k)  ) * ( v2(i,j,k-1) + v2(i,j-1,k-1) ) )

        dv3_dlam  = ( ( v3(i+1,j,k+1) + v3(i+1,j,k) )                   &
          &         - ( v3(i-1,j,k+1) + v3(i-1,j,k) ) ) * 0.25_wp * eddlon
        dv3_dphi  = ( ( v3(i,j+1,k+1) + v3(i,j+1,k) )                   &
          &         - ( v3(i,j-1,k+1) + v3(i,j-1,k) ) ) * 0.25_wp * eddlat
        dv3_dzeta =   ( v3(i,j,k+1)   - v3(i,j,k) )

        curl1(i,j,k) = r_earth_inv * ( dv3_dphi + dzeta_dphi(i,j,k)*dv3_dzeta ) &
          &          + sqrtg_r_s(i,j,k) * dv2_dzeta

        curl2(i,j,k) = -sqrtg_r_s(i,j,k) * dv1_dzeta                            &
          &          - acrlat(j,1) * ( dv3_dlam + dzeta_dlam(i,j,k)*dv3_dzeta )

        curl3(i,j,k) = acrlat(j,1) * ( dv2_dlam + dzeta_dlam(i,j,k)*dv2_dzeta ) &
          &          - r_earth_inv * ( dv1_dphi + dzeta_dphi(i,j,k)*dv1_dzeta )

      END DO
    END DO
  END DO

  ! top level
  k=1
  DO j=2, je-1
    DO i=2, ie-1

        ! derivatives at the scalar position:
        dv1_dlam  =   ( v1(i,j,k)   - v1(i-1,j,k) ) * eddlon
        dv1_dphi  = ( ( v1(i,j+1,k) + v1(i-1,j+1,k) )                          &
          &         - ( v1(i,j-1,k) + v1(i-1,j-1,k) ) ) * 0.25_wp * eddlat
        dv1_dzeta = ( ( v1(i,j,k+1) + v1(i-1,j,k+1) )                          &
          &         - ( v1(i,j,k)   + v1(i-1,j,k  ) ) ) * 0.5_wp

        dv2_dlam  = ( ( v2(i+1,j,k) + v2(i+1,j-1,k) )                          &
          &         - ( v2(i-1,j,k) + v2(i-1,j-1,k) ) ) * 0.25_wp * eddlon
        dv2_dphi  =   ( v2(i,j,k)   - v2(i,j-1,k) ) * eddlat
        dv2_dzeta = ( ( v2(i,j,k+1) + v2(i,j-1,k+1) )                          &
          &         - ( v2(i,j,k)   + v2(i,j-1,k) ) ) * 0.5_wp

        dv3_dlam  = ( ( v3(i+1,j,k+1) + v3(i+1,j,k) )                          &
          &         - ( v3(i-1,j,k+1) + v3(i-1,j,k) ) ) * 0.25_wp * eddlon
        dv3_dphi  = ( ( v3(i,j+1,k+1) + v3(i,j+1,k) )                          &
          &         - ( v3(i,j-1,k+1) + v3(i,j-1,k) ) ) * 0.25_wp * eddlat
        dv3_dzeta =   ( v3(i,j,k+1)   - v3(i,j,k) )

        curl1(i,j,k) = r_earth_inv * ( dv3_dphi + dzeta_dphi(i,j,k)*dv3_dzeta ) &
          &          + sqrtg_r_s(i,j,k) * dv2_dzeta

        curl2(i,j,k) = -sqrtg_r_s(i,j,k) * dv1_dzeta                            &
          &          - acrlat(j,1) * ( dv3_dlam + dzeta_dlam(i,j,k)*dv3_dzeta )

        curl3(i,j,k) = acrlat(j,1) * ( dv2_dlam + dzeta_dlam(i,j,k)*dv2_dzeta ) &
          &          - r_earth_inv * ( dv1_dphi + dzeta_dphi(i,j,k)*dv1_dzeta )

    END DO
  END DO

  ! bottom level
  k = ke
  DO j=2, je-1
    DO i=2, ie-1

        ! derivatives at the scalar position:
        dv1_dlam  =   ( v1(i,j,k)   - v1(i-1,j,k) ) * eddlon
        dv1_dphi  = ( ( v1(i,j+1,k) + v1(i-1,j+1,k) )                    &
          &         - ( v1(i,j-1,k) + v1(i-1,j-1,k) ) ) * 0.25_wp * eddlat
        dv1_dzeta = ( ( v1(i,j,k)   + v1(i-1,j,k  ) )                    &
          &         - ( v1(i,j,k-1) + v1(i-1,j,k-1) ) ) * 0.5_wp

        dv2_dlam  = ( ( v2(i+1,j,k) + v2(i+1,j-1,k) )                    &
          &         - ( v2(i-1,j,k) + v2(i-1,j-1,k) ) ) * 0.25_wp *eddlon
        dv2_dphi  =   ( v2(i,j,k)   - v2(i,j-1,k) ) * eddlat
        dv2_dzeta = ( ( v2(i,j,k)   + v2(i,j-1,k) )                      &
          &         - ( v2(i,j,k-1) + v2(i,j-1,k-1) ) ) * 0.5_wp

        dv3_dlam  = ( ( v3(i+1,j,k+1) + v3(i+1,j,k) )                    &
          &         - ( v3(i-1,j,k+1) + v3(i-1,j,k) ) ) * 0.25_wp * eddlon
        dv3_dphi  = ( ( v3(i,j+1,k+1) + v3(i,j+1,k) )                    &
          &         - ( v3(i,j-1,k+1) + v3(i,j-1,k) ) ) * 0.25_wp * eddlat
        dv3_dzeta =   ( v3(i,j,k+1)   - v3(i,j,k) )

        curl1(i,j,k) = r_earth_inv * ( dv3_dphi + dzeta_dphi(i,j,k)*dv3_dzeta ) &
          &          + sqrtg_r_s(i,j,k) * dv2_dzeta

        curl2(i,j,k) = -sqrtg_r_s(i,j,k) * dv1_dzeta                            &
          &          - acrlat(j,1) * ( dv3_dlam + dzeta_dlam(i,j,k)*dv3_dzeta )

        curl3(i,j,k) = acrlat(j,1) * ( dv2_dlam + dzeta_dlam(i,j,k)*dv2_dzeta ) &
          &          - r_earth_inv * ( dv1_dphi + dzeta_dphi(i,j,k)*dv1_dzeta )

    END DO
  END DO

  IF ( lmetr ) THEN
    ! Corections due to the spherical shape of the earth
    DO k=1, ke
      DO j=2, je-1
        DO i=2, ie-1

          v1_at_s = 0.5_wp * ( v1(i,j,k) + v1(i-1,j,  k) )
          v2_at_s = 0.5_wp * ( v2(i,j,k) + v2(i,  j-1,k) )

          curl1(i,j,k) = curl1(i,j,k) - r_earth_inv * v2_at_s
          curl2(i,j,k) = curl2(i,j,k) + r_earth_inv * v1_at_s
          ! use tgrlat(j,1) for mass grid point and u-latitudes
          curl3(i,j,k) = curl3(i,j,k) + r_earth_inv * tgrlat(j,1) * v1_at_s

        END DO
      END DO
    END DO
  END IF

  IF ( .NOT. l_phys_comp ) THEN
    ! calculate the contravariant components from the physical components

    DO k=1, ke
      DO j=2, je-1
        DO i=2, ie-1

          curl1(i,j,k) = acrlat(j,1) * curl1(i,j,k)
          curl2(i,j,k) = r_earth_inv * curl2(i,j,k)
          ! curl3: contravar. and physical components are identical

        END DO
      END DO
    END DO

  END IF

  ! set boundary values to 0    (at interior processor boundaries 
  !   a possible exchange must be done elsewhere!)
  curl1(1, :, :) = 0.0_wp
  curl1(ie,:, :) = 0.0_wp
  curl1(:, 1, :) = 0.0_wp
  curl1(:, je,:) = 0.0_wp

  curl2(1, :, :) = 0.0_wp
  curl2(ie,:, :) = 0.0_wp
  curl2(:, 1, :) = 0.0_wp
  curl2(:, je,:) = 0.0_wp

  curl3(1 ,:, :) = 0.0_wp
  curl3(ie,:, :) = 0.0_wp
  curl3(:, 1, :) = 0.0_wp
  curl3(:, je,:) = 0.0_wp

END SUBROUTINE curl

!==============================================================================
!==============================================================================

SUBROUTINE calc_Theta_Tppp( t, pp, p0, ie, je, ke, r_d, cp_d, Theta )

!------------------------------------------------------------------------------
!
! Description:
!   calculate potential temperature Theta from T' and p'
!------------------------------------------------------------------------------

USE data_constants, ONLY: p0ref_recip

INTEGER (KIND=iintegers), INTENT(IN)    ::  &
  ie, je, ke

REAL (KIND=wp),     INTENT(IN)          ::  &
  t  (ie,je,ke),   & ! full temperature
  pp (ie,je,ke),   & ! pressure deviation
  p0 (ie,je,ke)      ! reference pressure

REAL (KIND=wp),     INTENT(IN)          ::  &
  r_d,          & ! gas constant for dry air
  cp_d            ! specific heat of dry air at constant pressure

REAL (KIND=wp),     INTENT(OUT)         ::  &
  Theta(ie,je,ke)

INTEGER  (KIND=iintegers) :: i, j, k
REAL     (KIND=wp)        :: rovcp

!------------------------------------------------------------------------------

  !PRINT*, "Subr. [calc_Theta_1] ..."

  rovcp = r_d/cp_d

  DO k=1, ke
    DO j=1, je
      DO i=1, ie

        Theta(i,j,k) = t(i,j,k) * EXP(-rovcp*LOG(               &
                          (p0(i,j,k)+pp(i,j,k)) *p0ref_recip) )

      END DO
    END DO
  END DO

END SUBROUTINE calc_Theta_Tppp

!==============================================================================
!==============================================================================

SUBROUTINE hadv_cd2 &
   ( s, sten, umw, vmw, dp0r, rdx, rdy, ie, je, istart, iend, jstart, jend)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine caculates the tendency of a prognostic variable
!   's' due to horizontal advection in a vertical layer k.
!   The advective tendency is added on the in/output array 'sten' for the
!   global tendency.
!   To include the metrics of a rotated lat/lon grid, the u-velocity on
!   input is set to the mass-weighted velocity umw = u*dp and the v-velocity
!   on input is set to v*dp*cos(phi). The grid-spacing in the
!   phi direction  has to be modified from a*dphi to a*dphi*cos(phi).
!
! Method:
!   The 2nd order centered difference Leapfrog scheme is used. The
!   discretization is described in Section 3.3 of the Scientific
!   Documentation.
!
!------------------------------------------------------------------------------
!
! Declarations:

! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER (KIND=iintegers),   INTENT (IN) ::  &
  ie, je,       & ! x- and y-dimension of the input/output arrays
  istart, iend, & ! computational start and end indices in i-direction
                  !(west-east for lat/lon grid)
  jstart, jend    ! computational start and end indices in j-direction
                  ! (south-north for lat/lon grid)

! Array arguments with intent(in):

REAL    (KIND=wp)       ,   INTENT(IN)  ::  &
  s(ie,je),     & ! scalar variable to be transported at current time-level
  umw(ie,je),   & ! mass-weighted zonal wind velocity at current time-level
  vmw(ie,je),   & ! mass-weighted meridional wind velocity (* cos(phi))
                  ! at current time-level
  dp0r(ie,je),  & ! reciprocal base state pressure thickness of layer k
  rdx(je) ,     & ! 0.5*reciprocal x-grid spacing (1/( a cos(phi)dlam )
  rdy(je)         ! 0.5*reciprocal y-grid spacing (1/( a cos(phi)dphi )

! Array arguments with intent(inout):

REAL    (KIND=wp)       ,   INTENT(INOUT) ::  &
  sten(ie,je)     ! tendency array of scalar variable s

!------------------------------------------------------------------------------
! Local variables

INTEGER (KIND=iintegers)    ::  &
  i, j            ! loop indices

REAL (KIND=wp)              ::  &
  z1, z2, z3, z4, zadv   ! intermediate storage

!------------------------------------------------------------------------------
! Begin Subroutine hadv_cd2
!------------------------------------------------------------------------------
  DO    j = jstart, jend
    DO  i = istart, iend
      z1        = umw(i  ,j  )  *( s(i+1,j) - s(i  ,j) )
      z2        = umw(i-1,j  )  *( s(i  ,j) - s(i-1,j) )
      z3        = vmw(i  ,j  )  *( s(i,j+1) - s(i,j  ) )
      z4        = vmw(i  ,j-1)  *( s(i,j  ) - s(i,j-1) )
      zadv      =  -  rdx(j)*(z1+z2) - rdy(j)*(z3+z4)
      sten(i,j) = sten(i,j) + zadv * dp0r(i,j)
    ENDDO
  ENDDO

END SUBROUTINE hadv_cd2

!==============================================================================
!==============================================================================

SUBROUTINE lap_2 (s, lap, crlato, crlatu, ie, je, ke)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine caculates the Laplacian (lap) of an input array (s),
!   normalized by the grid spacing.
!
! Method:
!   The 2nd order centered differences. 
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------
! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER (KIND=iintegers),   INTENT (IN) ::  &
  ie, je, ke      ! x-, y- and z-dimension of the input/output arrays

! Array arguments with intent(in):

REAL    (KIND=wp)       ,   INTENT(IN)  ::  &
  s(ie,je,ke),  & ! scalar variable to be diffused at current time-level
  crlato(je),   & ! cosine ratio at j+1/2 for metrics
  crlatu(je)      ! cosine ratio at j-1/2 for metrics

! Array arguments with intent(out):

REAL    (KIND=wp)       ,   INTENT(OUT) ::  &
  lap  (ie,je,ke) ! Laplacian of the input array

!------------------------------------------------------------------------------
! Local variables

INTEGER (KIND=iintegers)    ::  &
  i, j, k         ! loop indices

!------------------------------------------------------------------------------

  DO k = 1, ke
!CDIR OUTERUNROLL=4
    DO   j  = 2 , je-1
!CDIR ON_ADB(lap)
!CDIR ON_ADB(s)
      DO i  = 2 , ie-1
        lap (i,j,k)   = s (i+1,j,k) + s (i-1,j,k) - 2.0_wp*s(i,j,k) &
                        +  crlato(j)*(s(i,j+1,k) - s(i,j  ,k))      &
                        -  crlatu(j)*(s(i,j  ,k) - s(i,j-1,k))
      ENDDO
    ENDDO
  
    ! set free slip boundary data
    ! Take care of the edges, which are not set
!CDIR ON_ADB(lap)
    lap (2:ie-1,1 ,k) = lap (2:ie-1,2   ,k)
!CDIR ON_ADB(lap)
    lap (2:ie-1,je,k) = lap (2:ie-1,je-1,k)
!CDIR ON_ADB(lap)
    lap (1,2:je-1 ,k) = lap (2   ,2:je-1,k)
!CDIR ON_ADB(lap)
    lap (ie,2:je-1,k) = lap (ie-1,2:je-1,k)
    lap (1,1,k)       = lap (2,2,k)
    lap (1,je,k)      = lap (2,je-1,k)
    lap (ie,1,k)      = lap (ie-1,2,k)
    lap (ie,je,k)     = lap (ie-1,je-1,k)
  ENDDO

END SUBROUTINE lap_2

!==============================================================================
!==============================================================================

SUBROUTINE lap_4 &
   (s, sten, work, crlato, crlatu, dcoeff, ie, je, istart, iend, jstart, jend)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine caculates the linear forth order horizontal diffusion
!   operator - the Laplacian of the Laplacian - of an input array 's' on
!   a layer/level and adds it - multiplied by a diffusion coefficient
!   'dcoeff' - to the tendency array 'sten' .
!
! Method:
!   The 2nd order centered differences. See Section 5.2 of the
!   Scientific Documentation.
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------
! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER (KIND=iintegers),   INTENT (IN) ::  &
  ie, je,       & ! x- and y-dimension of the input/output arrays
  istart, iend, & ! computational start and end indices in i-direction
                  !(west-east for lat/lon grid)
  jstart, jend    ! computational start and end indices in j-direction
                  ! (south-north for lat/lon grid)

! Array arguments with intent(in):

REAL    (KIND=wp)       ,   INTENT(IN)  ::  &
  s(ie,je) ,    & ! scalar variable to be diffused at current time-level
  crlato(je),   & ! cosine ratio at j+1/2 for metrics
  crlatu(je),   & ! cosine ratio at j-1/2 for metrics
  dcoeff          ! 4th-order diffusion coefficient

! Array arguments with intent(out):

REAL    (KIND=wp)       ,   INTENT(OUT) ::  &
  work (ie,je)    ! work-array provided by the calling routine

! Array arguments with intent(inout):

REAL    (KIND=wp)       ,   INTENT(INOUT) ::  &
  sten (ie,je)    ! tendency array of scalar variable s where the diffusive
                  ! tendency is added on

!------------------------------------------------------------------------------
! Local variables

INTEGER (KIND=iintegers)    ::  &
  i, j            ! loop indices

!------------------------------------------------------------------------------

  DO   j  = jstart-1 , jend+1
    DO i  = istart-1 , iend+1
      work (i,j)     = s (i+1,j) + s (i-1,j) - 2.0_wp*s(i,j) &
                       +  crlato(j)*(s(i,j+1) - s(i,j  ))    &
                       -  crlatu(j)*(s(i,j  ) - s(i,j-1))
    ENDDO
  ENDDO
  DO    j = jstart , jend
    DO  i = istart , iend
      sten  (i,j) =  sten(i,j) - dcoeff  &
                     * ( work(i+1,j) + work(i-1,j) - 2.0_wp*work(i,j) &
                       + crlato(j)*(work(i,j+1) - work(i,j  ))        &
                       - crlatu(j)*(work(i,j  ) - work(i,j-1)) )
    ENDDO
  ENDDO

END SUBROUTINE lap_4

!==============================================================================
!==============================================================================

SUBROUTINE lap_4a &
   ( sup, s, work, crlato, crlatu, dcoeff, ie, je, istart, iend, jstart, jend)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine caculates the linear forth order horizontal diffusion
!   operator - the Laplacian of the Laplacian - of an input array 's' on
!   a layer/level and updates - multiplied by a diffusion coefficient
!   'dcoeff' - the corresponing field 'sup' accordingly.
!
! Method:
!   The 2nd order centered differences. See Section 5.2 of the
!   Scientific Documentation.
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------
! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER (KIND=iintegers),   INTENT (IN) ::  &
  ie, je,       & ! x- and y-dimension of the input/output arrays
  istart, iend, & ! computational start and end indices in i-direction
                  !(west-east for lat/lon grid)
  jstart, jend    ! computational start and end indices in j-direction
                  ! (south-north for lat/lon grid)

! Array arguments with intent(in):

REAL    (KIND=wp)       ,   INTENT(IN)  ::  &
  s(ie,je) ,    & ! scalar variable to be diffused at current time-level
  crlato(je),   & ! cosine ratio at j+1/2 for metrics
  crlatu(je),   & ! cosine ratio at j-1/2 for metrics
  dcoeff(ie,je)   ! 4th-order diffusion coefficient times the timestep
                  ! for the 3D masked array

! Array arguments with intent(inout):

REAL    (KIND=wp)       ,   INTENT(INOUT) ::  &
  work (ie,je), & ! work-array provided by the calling routine
  sup  (ie,je)    ! input field to be updated by diffusion of variable 's'

!------------------------------------------------------------------------------
! Local variables

INTEGER (KIND=iintegers)    ::  &
  i, j            ! loop indices

!------------------------------------------------------------------------------

  DO   j  = jstart-1 , jend+1
    DO i  = istart-1 , iend+1
      work (i,j)     = s (i+1,j) + s (i-1,j) - 2.0_wp*s(i,j)  &
                       +  crlato(j)*(s(i,j+1) - s(i,j  ))     &
                       -  crlatu(j)*(s(i,j  ) - s(i,j-1))
    ENDDO
  ENDDO
  DO    j = jstart , jend
    DO  i = istart , iend
      sup   (i,j) =  sup (i,j) - dcoeff (i,j)                          &
                     * ( work(i+1,j) + work(i-1,j) - 2.0_wp*work(i,j)  &
                       + crlato(j)*(work(i,j+1) - work(i,j  ))         &
                       - crlatu(j)*(work(i,j  ) - work(i,j-1)) )
    ENDDO
  ENDDO

END SUBROUTINE lap_4a

!==============================================================================
!==============================================================================

SUBROUTINE lap_4am &
   ( sup, s, work, crlato, crlatu, dcoeff, ie, je, ke,                   &
     istart, iend, jstart, jend, kstart, kend )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine caculates the linear forth order horizontal diffusion
!   operator - the Laplacian of the Laplacian - of an input array 's' on
!   a layer/level and updates - multiplied by a diffusion coefficient
!   'dcoeff' - the corresponing field 'sup' accordingly.
!   A simple monotonic flux limiter is applied. 
!
! Method:
!   The 2nd order centered differences. See Section 5.2 of the
!   Scientific Documentation. Flux limiter according to Ming Xue
!   (Mon. Wea. Rev., Vol. 128, 2000).
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------
! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER (KIND=iintegers),   INTENT (IN) ::  &
  ie, je, ke,   & ! x-, y- and z-dimension of the input/output arrays
  istart, iend, & ! computational start and end indices in i-direction
                  !(west-east for lat/lon grid)
  jstart, jend, & ! computational start and end indices in j-direction
                  ! (south-north for lat/lon grid)
  kstart, kend    ! computational start and end indices in k-direction

! Array arguments with intent(in):

REAL    (KIND=wp)       ,   INTENT(IN)  ::  &
  s(ie,je,ke),  & ! scalar variable to be diffused at current time-level
  crlato(je),   & ! cosine ratio at j+1/2 for metrics
  crlatu(je),   & ! cosine ratio at j-1/2 for metrics
  dcoeff(ie,je,ke)! 4th-order diffusion coefficient times the timestep
                  ! for the 3D masked array

! Array arguments with intent(inout):

REAL    (KIND=wp)       ,   INTENT(INOUT) ::  &
  work (ie,je,ke), & ! work-array provided by the calling routine
  sup  (ie,je,ke)    ! input field to be updated by diffusion of variable 's'

!------------------------------------------------------------------------------


! Local variables
REAL    (KIND=wp)       ::  &
  flx  (ie,je), & ! longitudinal diffusive flux
  fly  (ie,je)    ! latitudinal  diffusive flux

INTEGER (KIND=iintegers)    ::  &
  i, j, k         ! loop indices

!------------------------------------------------------------------------------

  DO k = kstart, kend

!CDIR OUTERUNROLL=4
    DO   j  = jstart-1 , jend+1
!CDIR ON_ADB(work)
!CDIR ON_ADB(s)
      DO i  = istart-1 , iend+1
        work (i,j,k)   = s (i+1,j,k) + s (i-1,j,k) - 2.0_wp*s(i,j,k)  &
                         +  crlato(j)*(s(i,j+1,k) - s(i,j  ,k))       &
                         -  crlatu(j)*(s(i,j  ,k) - s(i,j-1,k))
      ENDDO
    ENDDO

    ! Calculate the (3rd order) diffusive fluxes and apply simple limiter
    DO   j  = jstart-1 , jend
!CDIR ON_ADB(work)
!CDIR ON_ADB(s)
!CDIR ON_ADB(flx)
!CDIR ON_ADB(fly)
      DO i  = istart-1 , iend
        flx(i,j) = work(i+1,j,k) - work(i,j,k)
        fly(i,j) = crlato(j)*(work(i,j+1,k) - work(i,j,k))
        IF ( flx(i,j)*(s(i+1,j,k)-s(i,j,k)) > 0.0_wp ) flx(i,j) = 0.0_wp
        IF ( fly(i,j)*(s(i,j+1,k)-s(i,j,k)) > 0.0_wp ) fly(i,j) = 0.0_wp
      ENDDO
    ENDDO

!CDIR OUTERUNROLL=4
    DO    j = jstart , jend
!CDIR ON_ADB(flx)
!CDIR ON_ADB(fly)
      DO  i = istart , iend
        sup   (i,j,k) =  sup (i,j,k) - dcoeff (i,j,k)     &
                         * ( flx(i,j) - flx(i-1,j)        &
                           + fly(i,j) - fly(i,j-1) )
      ENDDO
    ENDDO

  END DO

END SUBROUTINE lap_4am

!==============================================================================
!==============================================================================

SUBROUTINE lap_4aml &
   ( sup, s, lap, fac_hdx, fac_hdy, crlato, crlatu, dcoeff,  &
     ie, je, ke, istart, iend, jstart, jend )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine caculates the linear forth order horizontal diffusion
!   operator - the Laplacian of the Laplacian - of an input array 's' on
!   a layer/level and updates - multiplied by a diffusion coefficient
!   'dcoeff' - the corresponing field 'sup' accordingly.
!   A simple monotonic flux limiter is applied. If the height difference
!   of neighbouring grid points exceeds a certain limit, the diffusion flux
!   is set to zero in the corresponding direction.
!
! Method:
!   The 2nd order centered differences. See Section 5.2 of the
!   Scientific Documentation. Flux limiter according to Ming Xue
!   (Mon. Wea. Rev., Vol. 128, 2000).
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------
! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER (KIND=iintegers),   INTENT (IN) ::  &
  ie, je, ke,   & ! x-, y- and z-dimension of the input/output arrays
  istart, iend, & ! computational start and end indices in i-direction
                  !(west-east for lat/lon grid)
  jstart, jend    ! computational start and end indices in j-direction
                  ! (south-north for lat/lon grid)

! Array arguments with intent(in):

REAL    (KIND=wp)       ,   INTENT(IN)  ::  &
  s  (ie,je,ke),& ! scalar variable to be diffused at current time-level
  lap(ie,je,ke),& ! Laplacian of the scalar input-variable
  crlato(je),   & ! cosine ratio at j+1/2 for metrics
  crlatu(je),   & ! cosine ratio at j-1/2 for metrics
  dcoeff(ie,je,ke)! 4th-order diffusion coefficient times the timestep
                  ! for the 3D masked array

REAL    (KIND=wp)       ,   INTENT(IN)  ::  &
  fac_hdx(ie,je,ke),& ! Limiter factor for reducing x-diff-flux 
  fac_hdy(ie,je,ke)   ! Limiter factor for reducing y-diff-flux 

! Array arguments with intent(inout):

REAL    (KIND=wp)       ,   INTENT(INOUT) ::  &
  sup  (ie,je,ke) ! input field to be updated by diffusion of variable 's'

!------------------------------------------------------------------------------

! Local variables
REAL    (KIND=wp)       ::  &
  flx  (ie,je), & ! longitudinal diffusive flux
  fly  (ie,je), & ! latitudinal  diffusive flux
  rxp  (ie,je), & ! limiter function for positive fluxes
  rxm  (ie,je)    ! limiter function for negative fluxes

REAL    (KIND=wp)       ::  &
  fmin, fmax, fluxin, fluxou

INTEGER (KIND=iintegers)    ::  &
  i, j, k         ! loop indices

!------------------------------------------------------------------------------

  DO k = 1, ke

    ! Calculate the (3rd order) diffusive fluxes 
!CDIR OUTERUNROLL=4
    DO   j  = 1 , je-1
!CDIR ON_ADB(lap)
!CDIR ON_ADB(flx)
!CDIR ON_ADB(fly)
      DO i  = 1 , ie-1
        flx(i,j) = dcoeff(i,j,k)*( lap(i+1,j,k) - lap(i,j,k) ) * fac_hdx(i,j,k)
        fly(i,j) = dcoeff(i,j,k)*( lap(i,j+1,k) - lap(i,j,k) ) * fac_hdy(i,j,k)*crlato(j)
      ENDDO
    ENDDO
    !Calculate limiter
!CDIR OUTERUNROLL=4
    DO   j  = jstart-1 , jend+1
!CDIR ON_ADB(s)
!CDIR ON_ADB(flx)
!CDIR ON_ADB(fly)
!CDIR ON_ADB(rxp)
!CDIR ON_ADB(rxm)
      DO i  = istart-1 , iend+1
        fmax   =   MAX( s(i-1,j,k), s(i,j,k), s(i+1,j,k), s(i,j-1,k), s(i,j+1,k) )
        fmin   =   MIN( s(i-1,j,k), s(i,j,k), s(i+1,j,k), s(i,j-1,k), s(i,j+1,k) )
        fluxin =   MAX(0.0_wp,flx(i-1,j)) - MIN(0.0_wp,flx(i,j)) &
                 + MAX(0.0_wp,fly(i,j-1)) - MIN(0.0_wp,fly(i,j))
        fluxou = - MIN(0.0_wp,flx(i-1,j)) + MAX(0.0_wp,flx(i,j))  &
                 - MIN(0.0_wp,fly(i,j-1)) + MAX(0.0_wp,fly(i,j))
        !IF (fluxin > 0.0 ) THEN
        !rxp (i,j) = ABS(fmax-s(i,j,k))/fluxin
        !ENDIF
        !IF (fluxou > 0.0 ) THEN
        !rxm (i,j) = ABS(fmin-s(i,j,k))/fluxou
        !ENDIF
        rxp (i,j) = ABS(fmax-s(i,j,k))/(fluxin + eps_div)
        rxm (i,j) = ABS(fmin-s(i,j,k))/(fluxou + eps_div)
      ENDDO
    ENDDO

    DO j = jstart-1, jend
!CDIR ON_ADB(flx)
!CDIR ON_ADB(fly)
!CDIR ON_ADB(rxp)
!CDIR ON_ADB(rxm)
      DO i = istart-1, iend
        IF(flx(i,j) >=  0.0_wp) THEN
           flx(i,j) = MIN(1.0_wp, rxp(i+1,j),rxm(i,j) )*flx(i,j)
        ELSE
           flx(i,j) = MIN(1.0_wp, rxp(i,j),rxm(i+1,j) )*flx(i,j)
        ENDIF
        IF(fly(i,j) >=  0.0_wp) THEN
           fly(i,j) = MIN(1.0_wp, rxp(i,j+1),rxm(i,j) )*fly(i,j)
        ELSE
           fly(i,j) = MIN(1.0_wp, rxp(i,j),rxm(i,j+1) )*fly(i,j)
        ENDIF
      ENDDO
    ENDDO

!CDIR OUTERUNROLL=4
    DO    j = jstart , jend
!CDIR ON_ADB(flx)
!CDIR ON_ADB(fly)
      DO  i = istart , iend
        sup   (i,j,k) =  sup (i,j,k) -  ( flx(i,j) - flx(i-1,j) &
                                        + fly(i,j) - fly(i,j-1) )
      ENDDO
    ENDDO

  ENDDO

END SUBROUTINE lap_4aml

!==============================================================================
!==============================================================================

SUBROUTINE hadv_pd &
   ( u, v, s, su, rdx, rdy, dt, ie, je, istart, iend, jstart, jend)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine updates a scalar prognostic
!   variable 's' due to horizontal advection in a vertical layer k.
!   The updated field is written on the output array 'su'.
!   To include the metrics of a rotated lat/lon grid the v-velocity on
!   input must be set to v*cos(phi) and the grid-spacing in the phi direction
!   has to be modified from a*dphi to a*dphi*cos(phi).
!
! Method:
!   Currently, the 2nd order vanLeer advection scheme is combined with the
!   Lin and Rood (MWR, 1996 (vol.124, no.9)) method to do the horizontal
!   advection fully 2D.
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------
! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER (KIND=iintegers),   INTENT (IN) ::  &
  ie, je,       & ! x- and y-dimension of the input/output arrays
  istart, iend, & ! computational start and end indices in i-direction
                  !(west-east for lat/lon grid)
  jstart, jend    ! computational start and end indices in j-direction
                  ! (south-north for lat/lon grid)

REAL    (KIND=wp)       ,   INTENT(IN)  ::  &
  dt              ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp)       ,   INTENT(IN)  ::  &
  u(ie,je),     & ! zonal wind velocity at current time-level
  v(ie,je),     & ! meridional wind velocity (* cos(phi)) at current time-level
  s(ie,je),     & ! scalar variable to be transported at current time-level
  rdx(je) ,     & ! reciprocal x-grid spacing (1/( a cos(phi)dlam )
  rdy(je)         ! reciprocal y-grid spacing (1/( a cos(phi)dphi )

! Array arguments with intent(out):

REAL    (KIND=wp)       ,   INTENT(INOUT) ::  &
  ! this variable had INTENT(OUT) before, but for the prognostic advection
  ! of cloud-ice, where hadv_pd is called twice, we need INOUT to save
  ! the values at the boundary zone of the total domain (are set in the
  ! calling routine)
  su(ie,je)       ! scalar variabel s updated due to horizontal advection

!------------------------------------------------------------------------------
! Local variables

INTEGER (KIND=iintegers)    ::  &
  i, j            ! loop indices

REAL (KIND=wp)              ::  &
  zdtbydx,      & ! const*dt/dx
  zdtbydy,      & ! const*dt/dy
  zdql, zdqr,   & ! Parameters related to
  zr,           & ! monotonic slopes
  zcrv, zcru,   & ! Courant numbers for intermediate updating
  zcrx, zcry,   & ! Courant numbers for final fluxes
  zsten           ! tendency due to horizontal advection

! Local automatic arrays
REAL (KIND=wp)              ::  &
  zfx(ie,je),   & ! Flux in x-direction at cell boundaries
  zfy(ie,je),   & ! Flux in y-direction at cell boundaries
  zsy(ie,je),   & ! intermediate update of s due to y-advection
  zsx(ie,je),   & ! intermediate update of s due to x-advection
  zdx(ie,je),   & ! monotonic slope in x-direction
  zdy(ie,je)      ! monotonic slope in y-direction

! End of header
!=======================================================================

!----------------------------------------------------------------------------
! Section 1: Compute y-adv of s with upwind for fully 2D advection
!            (Lin & Rood, MWR 1996)
!----------------------------------------------------------------------------

  DO j = jstart-1, jend+1
    zdtbydy = 0.125_wp*dt*rdy(j)
    DO i = istart-2, iend+2
      zcrv = ( v(i,j) + v(i,j-1) ) * zdtbydy
      zsy(i,j) = s(i,j) - zcrv*(s(i,j+1)-s(i,j-1)) &
                 + ABS(zcrv) * (s(i,j+1)+s(i,j-1)-2._wp*s(i,j))
    ENDDO
  ENDDO
  DO j = jstart-2, jend+2
    zdtbydx = 0.125_wp*dt*rdx(j)
    DO i = istart-1, iend+1
      zcru = ( u(i,j) + u(i-1,j) ) * zdtbydx
      zsx(i,j) = s(i,j) - zcru*(s(i+1,j)-s(i-1,j)) &
                 + ABS(zcru) * (s(i+1,j)+s(i-1,j)-2._wp*s(i,j))
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section 2: Compute monotonic slopes and form fluxes in x- and in
!            y- direction
!------------------------------------------------------------------------------

  ! Compute monotonic slopes
  DO j = jstart-1, jend+1
    DO i = istart-1, iend+1
      zdql     = zsy(i  ,j) - zsy(i-1,j)
      zdqr     = zsy(i+1,j) - zsy(i  ,j)
      zr       = zdqr / (zdql + eps_div)
      zdx(i,j) = zdql * MAX(0.0_wp, MIN(2.0_wp, 2.0_wp*zr, 0.5_wp*(1.0_wp+zr)))
      zdql     = zsx(i,j  ) - zsx(i,j-1)
      zdqr     = zsx(i,j+1) - zsx(i,j  )
      zr       = zdqr / (zdql + eps_div)
      zdy(i,j) = zdql * MAX(0.0_wp, MIN(2.0_wp, 2.0_wp*zr, 0.5_wp*(1.0_wp+zr)))
    ENDDO
  ENDDO

  ! Compute courant number and velo-flux at zone edges (u-v points)
  DO j = jstart-1, jend
    zdtbydx = 0.5_wp*dt*rdx(j)
    zdtbydy = 0.5_wp*dt*rdy(j)
    DO i = istart-1, iend
      zcrx     =   zdtbydx*u(i,j)
      zfx(i,j) =   MAX(u(i,j),0.0_wp)*(zsy(i  ,j) + (0.5_wp-zcrx)*zdx(i  ,j))&
                 + MIN(u(i,j),0.0_wp)*(zsy(i+1,j) - (0.5_wp+zcrx)*zdx(i+1,j))
      zcry     =   zdtbydy*v(i,j)
      zfy(i,j) =   MAX(v(i,j),0.0_wp)*(zsx(i,j  ) + (0.5_wp-zcry)*zdy(i,j  ))&
                 + MIN(v(i,j),0.0_wp)*(zsx(i,j+1) - (0.5_wp+zcry)*zdy(i,j+1))
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section 3: Update the zones
!------------------------------------------------------------------------------
  DO j = jstart, jend
    DO i = istart, iend
      zsten   =  - rdx(j) * ( (zfx(i,j) - zfx(i-1,j)) &
                         - s(i,j)*(u(i,j)-u(i-1,j)) ) &
                 - rdy(j) * ( (zfy(i,j) - zfy(i,j-1)) &
                         - s(i,j)*(v(i,j)-v(i,j-1)) )
      su(i,j) = MAX ( 0.0_wp, s(i,j) + zsten*dt )
    ENDDO
  ENDDO

END SUBROUTINE hadv_pd

!==============================================================================
!==============================================================================

SUBROUTINE vadv_pd &
   ( wc, dp0, st, dt, ie, je, ke, ke1, istart, iend, jstart, jend)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine updates a scalar prognostic
!   variable 's' due to vertical advection. The routine uses the
!   directional split method, i.e. vertical advection operates on
!   the updated field s = sold + st*dt on input.
!   The field is calculated at the end and overwrites the
!   incomming field s.
!
! Method:
!   Currently, the 2nd order vanLeer advection scheme is used.
!
! Declarations:
!
!------------------------------------------------------------------------------

! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER (KIND=iintegers),   INTENT (IN) ::  &
  ie, je,       & ! x- and y-dimension of the input/output arrays
  ke, ke1,      & ! z-dimension of the input/output arrays
  istart, iend, & ! computational start and end indices in i-direction
                  !(west-east for lat/lon grid)
  jstart, jend    ! computational start and end indices in j-direction
                  ! (south-north for lat/lon grid)

REAL    (KIND=wp)       ,   INTENT(IN)  ::  &
  dt              ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp)       ,   INTENT(IN)  ::  &
  wc (ie,je,ke1), & ! mass-weighted contravariant vertical velocity
                    ! (at half levels)
  dp0(ie,je,ke)     ! base-state-pressure thickness

! Array arguments with intent(out):

REAL    (KIND=wp)       ,   INTENT(INOUT) ::  &
  st(ie,je,ke)      ! scalar field on input (updated due to horizontal
                    ! advection) to be transportet vertically;
                    ! updated due to vertical advection on output.

!------------------------------------------------------------------------------
! Local variables

INTEGER (KIND=iintegers)    ::  &
  i, j, k         ! loop indices

REAL (KIND=wp)              ::  &
  zdqtop, zdqbot,  & ! Parameters related to
  zr,              & ! monotonic slopes
  zsnew,           & ! New s-value after vertical advection
  zcrzt, zcrzb,    & ! Courant numbers for intermediate updating
  zwccfl,          & ! max. cfl
  zeddt              ! 1.0 / dt

! Local automatic arrays
REAL (KIND=wp)              ::  &
  zfz (ie,je,ke1), & ! Flux in zeta-direction at cell boundary (w-point)
  zdzm(ie,je),     & ! monotonic slopes for upper and
  zdzp(ie,je),     & ! lower grid point (resp. the w-point)
  zwc(ie,je,ke1)     ! mass-weighted contravariant vertical velocity

LOGICAL lexist(ke)

! End of header
!=======================================================================

zeddt = 1.0_wp / dt

DO k = 1, ke
  lexist(k) = MAXVAL(st(:,:,k)) > 0.0_wp
ENDDO


!----------------------------------------------------------------------------
! Section 1: Compute monotonic slopes and form fluxes in zeta - direction
!----------------------------------------------------------------------------

  ! limit vertical velocity to maximum value (courant number 1)
  zwc(:,:,  1) = wc(:,:,  1)
  zwc(:,:,ke1) = wc(:,:,ke1)
  DO k = 2, ke
    DO j = jstart, jend
      DO i = istart, iend
         zwccfl     = 0.5_wp*(dp0(i,j,k-1)+dp0(i,j,k)) * zeddt
         zwc(i,j,k) = MIN(ABS(wc(i,j,k)),zwccfl)*SIGN(1.0_wp,wc(i,j,k))
      ENDDO
    ENDDO
  ENDDO

  zdzm(:,:) = 0.0_wp
  zfz (:,:,:) = 0.0_wp

  IF (lexist(1).OR.lexist(2)) THEN
    DO j = jstart, jend
      DO i = istart, iend
        zdzm(i,j ) = st(i,j,2 ) - st(i,j,1   )
      ENDDO
    ENDDO
  ENDIF
  ! Compute courant number and velo-flux at zone edges (w-points)
  DO k = 2, ke
    zdzp(:,:) = 0.0_wp
    IF ( k == ke ) THEN
      IF (lexist(ke).OR.lexist(ke-1)) THEN
        DO j = jstart, jend
          DO i = istart, iend
            zdzp(i,j) = st(i,j,ke) - st(i,j,ke-1)
          ENDDO
        ENDDO
      ENDIF
    ELSE ! Compute monotonic slopes
      IF (lexist(k).OR.lexist(k-1).OR.lexist(k+1)) THEN
      DO j = jstart, jend
        DO i = istart, iend
          zdqtop   = st(i,j,k  ) - st(i,j,k-1)
          zdqbot   = st(i,j,k+1) - st(i,j,k  )
          zr       = zdqbot / (zdqtop + eps_div)
          zdzp(i,j) = zdqtop * MAX(0.0_wp, MIN(2.0_wp, 2.0_wp*zr, 0.5_wp*(1.0_wp+zr)))
        ENDDO
      ENDDO
      ENDIF
    ENDIF
    ! Compute fluxes
    DO j = jstart, jend
      DO i = istart, iend
      zcrzt   =   0.5_wp*dt*zwc(i,j,k)/dp0(i,j,k-1)
      zcrzb   =   0.5_wp*dt*zwc(i,j,k)/dp0(i,j,k  )
      zfz(i,j,k) =                                                           &
             MAX(zwc(i,j,k),0.0_wp)*(st(i,j,k-1) + (0.5_wp-zcrzt)*zdzm(i,j)) &
           + MIN(zwc(i,j,k),0.0_wp)*(st(i,j,k  ) - (0.5_wp+zcrzb)*zdzp(i,j))
      zdzm(i,j) = zdzp(i,j)
      ENDDO
    ENDDO
  ENDDO

!----------------------------------------------------------------------------
! Section 2: Update the zones
!----------------------------------------------------------------------------
  DO k = 1, ke
  DO j = jstart, jend
    DO i = istart, iend
      zsnew     =  st(i,j,k) - dt/dp0(i,j,k)  &
                      * ( (zfz(i,j,k+1) - zfz(i,j,k)) &
                          - st(i,j,k)*(zwc(i,j,k+1)-zwc(i,j,k)) )
      st(i,j,k) = MAX ( 0.0_wp, zsnew )
    ENDDO
  ENDDO
  ENDDO

END SUBROUTINE vadv_pd

!==============================================================================
!==============================================================================

SUBROUTINE hadv_pd_2tl &
   ( u, v, s, su, rdx, rdy, dp0r, dt, ie, je, istart, iend, jstart, jend)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine updates a scalar prognostic
!   variable 's' due to horizontal advection in a vertical layer k.
!   The updated field is written on the output array 'su'.
!   To include the metrics of a rotated lat/lon grid the v-velocity on
!   input must be set to v*cos(phi) and the grid-spacing in the phi direction
!   has to be modified from a*dphi to a*dphi*cos(phi).
!   Mass weighted velocities are required.
!
! Method:
!   Currently, the 2nd order vanLeer advection scheme is combined with the
!   Lin and Rood (MWR, 1996 (vol.124, no.9)) method to do the horizontal
!   advection fully 2D.
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------
! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER (KIND=iintegers),   INTENT (IN) ::  &
  ie, je,       & ! x- and y-dimension of the input/output arrays
  istart, iend, & ! computational start and end indices in i-direction
                  !(west-east for lat/lon grid)
  jstart, jend    ! computational start and end indices in j-direction
                  ! (south-north for lat/lon grid)

REAL    (KIND=wp)       ,   INTENT(IN)  ::  &
  dt              ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp)       ,   INTENT(IN)  ::  &
  u(ie,je),     & ! zonal wind velocity at current time-level
  v(ie,je),     & ! meridional wind velocity (* cos(phi)) at current time-level
  s(ie,je),     & ! scalar variable to be transported at current time-level
  rdx(je) ,     & ! reciprocal x-grid spacing (1/( a cos(phi)dlam )
  rdy(je) ,     & ! reciprocal y-grid spacing (1/( a cos(phi)dphi )
  dp0r(ie,je)     ! reciprocal base state pressure thickness

! Array arguments with intent(out):

REAL    (KIND=wp)       ,   INTENT(OUT) ::  &
  su(ie,je)       ! scalar variabel s updated due to horizontal advection

!------------------------------------------------------------------------------
! Local variables

INTEGER (KIND=iintegers)    ::  &
  i, j            ! loop indices

REAL (KIND=wp)              ::  &
  zdtbydx,      & ! const*dt/dx
  zdtbydy,      & ! const*dt/dy
  zdql, zdqr,   & ! Parameters related to
  zr,           & ! monotonic slopes
  zcrv, zcru,   & ! Courant numbers for intermediate updating
  zcrx, zcry,   & ! Courant numbers for final fluxes
  zsten           ! tendency due to horizontal advection

! Local automatic arrays
REAL (KIND=wp)              ::  &
  zfx(ie,je),   & ! Flux in x-direction at cell boundaries
  zfy(ie,je),   & ! Flux in y-direction at cell boundaries
  zsy(ie,je),   & ! intermediate update of s due to y-advection
  zsx(ie,je),   & ! intermediate update of s due to x-advection
  zdx(ie,je),   & ! monotonic slope in x-direction
  zdy(ie,je)      ! monotonic slope in y-direction

! End of header
!=======================================================================

!----------------------------------------------------------------------------
! Section 1: Compute y-adv of s with upwind for fully 2D advection
!            (Lin & Rood, MWR 1996)
!----------------------------------------------------------------------------

  DO j = jstart-1, jend+1
    zdtbydy = 0.125_wp*dt*rdy(j)
    DO i = istart-2, iend+2
      zcrv = ( v(i,j) + v(i,j-1) ) * zdtbydy * dp0r(i,j)
      zsy(i,j) = s(i,j) - zcrv*(s(i,j+1)-s(i,j-1)) &
                 + ABS(zcrv) * (s(i,j+1)+s(i,j-1)-2._wp*s(i,j))
    ENDDO
  ENDDO
  DO j = jstart-2, jend+2
    zdtbydx = 0.125_wp*dt*rdx(j)
    DO i = istart-1, iend+1
      zcru = ( u(i,j) + u(i-1,j) ) * zdtbydx * dp0r(i,j)
      zsx(i,j) = s(i,j) - zcru*(s(i+1,j)-s(i-1,j)) &
                 + ABS(zcru) * (s(i+1,j)+s(i-1,j)-2._wp*s(i,j))
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section 2: Compute monotonic slopes and form fluxes in x- and in
!            y- direction
!------------------------------------------------------------------------------

  ! Compute monotonic slopes
  DO j = jstart-1, jend+1
    DO i = istart-1, iend+1
      zdql     = zsy(i  ,j) - zsy(i-1,j)
      zdqr     = zsy(i+1,j) - zsy(i  ,j)
      zr       = zdqr / (zdql + eps_div)
      zdx(i,j) = zdql * MAX(0.0_wp, MIN(2.0_wp, 2.0_wp*zr, 0.5_wp*(1.0_wp+zr)))
      zdql     = zsx(i,j  ) - zsx(i,j-1)
      zdqr     = zsx(i,j+1) - zsx(i,j  )
      zr       = zdqr / (zdql + eps_div)
      zdy(i,j) = zdql * MAX(0.0_wp, MIN(2.0_wp, 2.0_wp*zr, 0.5_wp*(1.0_wp+zr)))
    ENDDO
  ENDDO

  ! Compute courant number and velo-flux at zone edges (u-v points)
  DO j = jstart-1, jend
    zdtbydx = 0.5_wp*dt*rdx(j)
    zdtbydy = 0.5_wp*dt*rdy(j)
    DO i = istart-1, iend
      zcrx     =   zdtbydx*u(i,j)
      zfx(i,j) =   MAX(u(i,j),0.0_wp)*(zsy(i  ,j) + &
                      (0.5_wp-zcrx*dp0r(i  ,j))*zdx(i  ,j))&
                 + MIN(u(i,j),0.0_wp)*(zsy(i+1,j) - &
                      (0.5_wp+zcrx*dp0r(i+1,j))*zdx(i+1,j))
      zcry     =   zdtbydy*v(i,j)
      zfy(i,j) =   MAX(v(i,j),0.0_wp)*(zsx(i,j  ) + &
                      (0.5_wp-zcry*dp0r(i,j  ))*zdy(i,j  ))&
                 + MIN(v(i,j),0.0_wp)*(zsx(i,j+1) - &
                      (0.5_wp+zcry*dp0r(i,j+1))*zdy(i,j+1))
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section 3: Update the zones
!------------------------------------------------------------------------------
  DO j = jstart, jend
    DO i = istart, iend
      zsten   =  - rdx(j) * ( (zfx(i,j) - zfx(i-1,j)) &
                         !- s(i,j)*(u(i,j)-u(i-1,j))  &
                            ) &
                 - rdy(j) * ( (zfy(i,j) - zfy(i,j-1)) &
                         !- s(i,j)*(v(i,j)-v(i,j-1))  &
                            )
      su(i,j) = MAX ( 0.0_wp, s(i,j) + zsten*dt*dp0r(i,j) )
    ENDDO
  ENDDO

END SUBROUTINE hadv_pd_2tl

!==============================================================================
!==============================================================================

SUBROUTINE vadv_pd_2tl &
   ( wc, dp0, st, dt, ie, je, ke, ke1, istart, iend, jstart, jend)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine updates a scalar prognostic
!   variable 's' due to vertical advection. The routine uses the
!   directional split method, i.e. vertical advection operates on
!   the updated field s = sold + st*dt on input.
!   The field is calculated at the end and overwrites the
!   incomming field s.
!
! Method:
!   Currently, the 2nd order vanLeer advection scheme is used.
!
! Declarations:
!
!------------------------------------------------------------------------------

! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER (KIND=iintegers),   INTENT (IN) ::  &
  ie, je,       & ! x- and y-dimension of the input/output arrays
  ke, ke1,      & ! z-dimension of the input/output arrays
  istart, iend, & ! computational start and end indices in i-direction
                  !(west-east for lat/lon grid)
  jstart, jend    ! computational start and end indices in j-direction
                  ! (south-north for lat/lon grid)

REAL    (KIND=wp)       ,   INTENT(IN)  ::  &
  dt              ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp)       ,   INTENT(IN)  ::  &
  wc (ie,je,ke1), & ! mass-weighted contravariant vertical velocity
                    ! (at half levels)
  dp0(ie,je,ke)     ! base-state-pressure thickness

! Array arguments with intent(out):

REAL    (KIND=wp)       ,   INTENT(INOUT) ::  &
  st(ie,je,ke)      ! scalar field on input (updated due to horizontal
                    ! advection) to be transportet vertically;
                    ! updated due to vertical advection on output.

!------------------------------------------------------------------------------
! Local variables

INTEGER (KIND=iintegers)    ::  &
  i, j, k         ! loop indices

REAL (KIND=wp)              ::  &
  zdqtop, zdqbot,  & ! Parameters related to
  zr,              & ! monotonic slopes
  zsnew,           & ! New s-value after vertical advection
  zcrzt, zcrzb,    & ! Courant numbers for intermediate updating
  zwccfl,          & ! max. cfl
  zeddt              ! 1.0 / dt

! Local automatic arrays
REAL (KIND=wp)              ::  &
  zfz (ie,je,ke1), & ! Flux in zeta-direction at cell boundary (w-point)
  zdzm(ie,je),     & ! monotonic slopes for upper and
  zdzp(ie,je),     & ! lower grid point (resp. the w-point)
  zwc(ie,je,ke1)     ! mass-weighted contravariant vertical velocity

LOGICAL lexist(ke)

! End of header
!=======================================================================

zeddt = 1.0_wp / dt

DO k = 1, ke
!  lexist(k) = MAXVAL(st(:,:,k)) > 0.0_wp
  lexist(k) = .TRUE.
ENDDO


!----------------------------------------------------------------------------
! Section 1: Compute monotonic slopes and form fluxes in zeta - direction
!----------------------------------------------------------------------------

  ! limit vertical velocity to maximum value (courant number 1)
  zwc(:,:,  1) = wc(:,:,  1)
  zwc(:,:,ke1) = wc(:,:,ke1)
  DO k = 2, ke
    DO j = jstart, jend
      DO i = istart, iend
         zwccfl     = 0.5_wp*(dp0(i,j,k-1)+dp0(i,j,k)) * zeddt
         zwc(i,j,k) = MIN(ABS(wc(i,j,k)),zwccfl)*SIGN(1.0_wp,wc(i,j,k))
      ENDDO
    ENDDO
  ENDDO

  zdzm(:,:) = 0.0_wp
  zfz (:,:,:) = 0.0_wp

  IF (lexist(1).OR.lexist(2)) THEN
    DO j = jstart, jend
      DO i = istart, iend
        zdzm(i,j ) = st(i,j,2 ) - st(i,j,1   )
      ENDDO
    ENDDO
  ENDIF
  ! Compute courant number and velo-flux at zone edges (w-points)
  DO k = 2, ke
    zdzp(:,:) = 0.0_wp
    IF ( k == ke ) THEN
      IF (lexist(ke).OR.lexist(ke-1)) THEN
        DO j = jstart, jend
          DO i = istart, iend
            zdzp(i,j) = st(i,j,ke) - st(i,j,ke-1)
          ENDDO
        ENDDO
      ENDIF
    ELSE ! Compute monotonic slopes
      IF (lexist(k).OR.lexist(k-1).OR.lexist(k+1)) THEN
      DO j = jstart, jend
        DO i = istart, iend
          zdqtop   = st(i,j,k  ) - st(i,j,k-1)
          zdqbot   = st(i,j,k+1) - st(i,j,k  )
          zr       = zdqbot / (zdqtop + eps_div)
          zdzp(i,j) = zdqtop * MAX(0.0_wp, MIN(2.0_wp, 2.0_wp*zr, 0.5_wp*(1.0_wp+zr)))
        ENDDO
      ENDDO
      ENDIF
    ENDIF
    ! Compute fluxes
    DO j = jstart, jend
      DO i = istart, iend
      zcrzt   =   0.5_wp*dt*zwc(i,j,k)/dp0(i,j,k-1)
      zcrzb   =   0.5_wp*dt*zwc(i,j,k)/dp0(i,j,k  )
      zfz(i,j,k) =                                                           &
             MAX(zwc(i,j,k),0.0_wp)*(st(i,j,k-1) + (0.5_wp-zcrzt)*zdzm(i,j)) &
           + MIN(zwc(i,j,k),0.0_wp)*(st(i,j,k  ) - (0.5_wp+zcrzb)*zdzp(i,j))
      zdzm(i,j) = zdzp(i,j)
      ENDDO
    ENDDO
  ENDDO

!----------------------------------------------------------------------------
! Section 2: Update the zones
!----------------------------------------------------------------------------
  DO k = 1, ke
  DO j = jstart, jend
    DO i = istart, iend
      zsnew     =  st(i,j,k) - dt/dp0(i,j,k) * ( &
                         (zfz(i,j,k+1) - zfz(i,j,k)) &
                          !- st(i,j,k)*(zwc(i,j,k+1)-zwc(i,j,k)) &
                         )
      st(i,j,k) = MAX(0.0_wp,zsnew)
    ENDDO
  ENDDO
  ENDDO

END SUBROUTINE vadv_pd_2tl

!==============================================================================
!==============================================================================

SUBROUTINE backtraj_trilin_dt1_3tl                                        &
                     (ut, vt, wt, eddlon, eddlat, eddzeta,                &
                      dt, ie, je, ke, istart, iend, jstart, jend, idx, wght )

!--------------------------------------------------------------------------
!
! Description:
!   This subroutine calculates the backtrajectory due to the transformed
!   velocity field (ut, vt, wt) in 1. order dt
!   (designed originally for 3 time level schemes from n-1 to n+1, but it
!   can be used also in 2tl schemes with dt/2).
!   To include the metrics of a rotated lat/lon grid, the velocity on
!   input has to be set to the transformed velocities:
!   ut = u / a / cos(phi), vt = v / a, wt = d zeta/ dt
!   (so the transport takes place completely in the transformed space
!   (lambda, phi, zeta) )
!   This routine has to be called only once for all the subsequent
!   semi-Lagrange-interpolation-routines.
!
! Output: fields
!   'idx' for index shifts and
!   'wght' for the interpolation weights.
!   They are ONLY useable for the subsequent 'trilinear' and
!   'tricubic' interpolation (not for 'triquadratic').
!
! Method:
!   Remark: wt has ke+1 vertical levels (compared to ke levels of 'zqiwt'
!   in 'src_leapfrog.f90'!)
!
!------------------------------------------------------------------------------

! Declarations:

INTEGER (KIND=iintegers),   INTENT (IN) ::  &
  ie, je, ke,   & ! x-, y- and z-dimension of the input/output arrays
  istart, iend, & ! computational start and end indices in i-direction
                  !(west-east for lat/lon grid)
  jstart, jend    ! computational start and end indices in j-direction
                  ! (south-north for lat/lon grid)

! Array arguments with intent(in):
REAL    (KIND=wp)       ,   INTENT(IN)  ::  &
         ut(ie,je,ke),   & ! transf. zonal wind veloc. at current time-level
         vt(ie,je,ke),   & ! transf. merid. wind veloc. at current time-level
         eddlon,         & ! reciprocal lambda-grid spacing
         eddlat,         & ! reciprocal phi   -grid spacing 
         eddzeta,        & ! reciprocal zeta  -grid spacing 
         dt                ! time step

REAL    (KIND=wp)       ,   INTENT(INOUT)  ::  &
         wt(ie,je,ke+1)    ! transf. vertical wind veloc. at current time-level

INTEGER (KIND=iintegers),   INTENT(OUT) :: idx (:,:,:,:)
REAL    (KIND=wp)       ,   INTENT(OUT) :: wght(:,:,:,:)

! Local variables
!----------------

INTEGER (KIND=iintegers)    :: i, j, k          ! loop indices

REAL (KIND=wp)              :: ax, ay, az       ! backtraj. in index space

!------------------------------------------------------------------------------

  wt(:,:,ke+1) = 0.0_wp   ! Boundary condition at the surface

  DO k=1, ke
     DO j=jstart, jend
        DO i=istart, iend

          ! calculation of the backward trajectory in 1. order dt
          ax =  dt * eddlon  * ( ut(i,j,k) + ut(i-1,j,  k)   ) 
          ay =  dt * eddlat  * ( vt(i,j,k) + vt(i,  j-1,k)   )
          az =  dt * eddzeta * ( wt(i,j,k) + wt(i,  j,  k+1) )

          idx(i,j,k,1) = FLOOR( ax )
          idx(i,j,k,2) = FLOOR( ay )
          idx(i,j,k,3) = FLOOR( az )

          wght(i,j,k,1) = ax - REAL( idx(i,j,k,1), wp)
          wght(i,j,k,2) = ay - REAL( idx(i,j,k,2), wp)
          wght(i,j,k,3) = az - REAL( idx(i,j,k,3), wp)

        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE backtraj_trilin_dt1_3tl

!==============================================================================
!==============================================================================

SUBROUTINE backtraj_trilin_dt2_3tl                                          &
                (ut, vt, wt, eddlon, eddlat, eddzeta, dt, ie, je, ke,       &
                 istart, iend, jstart, jend, idx, wght )

!--------------------------------------------------------------------------
!
! Description:
!   This subroutine calculates the backtrajectory due to the transformed
!   velocity field (ut, vt, wt) in 2. order dt
!   (designed originally for 3 time level schemes from n-1 to n+1, but it
!   can be used also in 2tl schemes with dt/2).
!   To include the metrics of a rotated lat/lon grid, the velocity on
!   input has to be set to the transformed velocities:
!   ut = u / a / cos(phi), vt = v / a, wt = d zeta/ dt
!   (so the transport takes place completely in the transformed space
!   (lambda, phi, zeta) )
!   This routine has to be called only once for all the subsequent
!   semi-Lagrange-interpolation-routines.
!
! Output: fields
!   'idx' for index shifts and
!   'wght' for the interpolation weights.
!   They are ONLY useable for the subsequent 'trilinear' and
!   'tricubic' interpolation (not for 'triquadratic').
!
! Method:
!   Remark: wt has ke+1 vertical levels (compared to ke levels of 'zqiwt'
!   in 'src_leapfrog.f90'!)
!
!------------------------------------------------------------------------------

! Declarations:

INTEGER (KIND=iintegers),   INTENT (IN) ::  &
  ie, je, ke,   & ! x-, y- and z-dimension of the input/output arrays
  istart, iend, & ! computational start and end indices in i-direction
                  !(west-east for lat/lon grid)
  jstart, jend    ! computational start and end indices in j-direction
                  ! (south-north for lat/lon grid)

REAL    (KIND=wp)       ,   INTENT(IN)  ::  &
  ut(ie,je,ke),   & ! transf. zonal wind velocity at current time-level
  vt(ie,je,ke),   & ! transf. merid. wind velocity at current time-level
  eddlon,         & ! reciprocal lambda-grid spacing
  eddlat,         & ! reciprocal phi   -grid spacing 
  eddzeta,        & ! reciprocal zeta  -grid spacing 
  dt                ! time step

REAL    (KIND=wp)       ,   INTENT(INOUT)  ::  &
  wt(ie,je,ke+1)    ! transf. vertical wind veloc. at current time-level

! fields describing the backtrajectory
INTEGER (KIND=iintegers), INTENT(OUT) :: idx (:,:,:,:) ! index shift
REAL    (KIND=wp)       , INTENT(OUT) :: wght(:,:,:,:) ! interpolation weights

! Local variables
!----------------

!US erst mal auskommentiert:
!US integer, parameter :: bndy = 0  ! depth of a boundary layer with 
!US                                 ! undefined velocity values

INTEGER (KIND=iintegers)    ::  &
  i, j, k,        &  ! loop indices
  p, q, r,        &  ! index shifts
  ip0, ip1, jq0, jq1, kr0, kr1    ! auxiliary index shifts      

REAL (KIND=wp)              ::  &
  wx, wy, wz, wx1, wy1, wz1,    &  ! interpolation weigths
  ax, ay, az,                   &  ! intermediate storage
  u_intpl, v_intpl, w_intpl        ! interpol. veloc. for 1. iteration

!------------------------------------------------------------------------------

  wt(:,:,ke+1) = 0.0_wp   ! Boundary condition at the surface

  DO k=1, ke
     DO j=jstart, jend
        DO i=istart, iend

           ! calculation of the backward trajectory in 1. order dt (over 1*dt)
           ax =  dt * eddlon  * 0.5_wp * ( ut(i,j,k) + ut(i-1,j,  k)   ) 
           ay =  dt * eddlat  * 0.5_wp * ( vt(i,j,k) + vt(i,  j-1,k)   )
           az =  dt * eddzeta * 0.5_wp * ( wt(i,j,k) + wt(i,  j,  k+1) )

!US erst mal auskommentiert: siehe oben
!US          IF  (   ( i-FLOOR( ax+0.5_wp ) -1 >= 1 +bndy ) .AND.       &
!US               &  ( i-FLOOR( ax )           <= ie-bndy ) .AND.       &
!US               &  ( j-FLOOR( ay+0.5_wp ) -1 >= 1 +bndy ) .AND.       &
!US               &  ( j-FLOOR( ay )           <= je-bndy ) .AND.       &

           IF  (   ( i-FLOOR( ax+0.5_wp ) -1 >= 1 ) .AND.       &
                &  ( i-FLOOR( ax )           <= ie) .AND.       &
                &  ( j-FLOOR( ay+0.5_wp ) -1 >= 1 ) .AND.       &
                &  ( j-FLOOR( ay )           <= je) .AND.       &
                &  ( k-FLOOR( az+0.5_wp ) -1 >= 1 ) .AND.       &
                &  ( k-FLOOR( az )           <= ke) )     THEN
              
              ! backtraj. in 1. order lies in the domain

              ! 1. iteration = backward traj. in 2. order dt

              ! u - component

              p = FLOOR( ax + 0.5_wp )
              q = FLOOR( ay )
              r = FLOOR( az )

              ip0 = i - p
              ip1 = ip0 - 1
              jq0 = j - q
              jq1 = jq0 - 1
              kr0 = k - r
              kr1 = kr0 - 1

              wx = ax + 0.5_wp - REAL( p, wp)
              wy = ay          - REAL( q, wp)
              wz = az          - REAL( r, wp)

              wx1 = 1.0_wp - wx
              wy1 = 1.0_wp - wy
              wz1 = 1.0_wp - wz

              ! trilinear interpolation of the velocity field u
              u_intpl =                                                &
                   &   wx  * (  wy  * (  wz  * ut(ip1, jq1, kr1 )      &
                   &                   + wz1 * ut(ip1, jq1, kr0 ) )    &
                   &          + wy1 * (  wz  * ut(ip1, jq0, kr1 )      &
                   &                   + wz1 * ut(ip1, jq0, kr0 ) ) )  &
                   & + wx1 * (  wy  * (  wz  * ut(ip0, jq1, kr1 )      & 
                   &                   + wz1 * ut(ip0, jq1, kr0 ) )    & 
                   &          + wy1 * (  wz  * ut(ip0, jq0, kr1 )      &
                   &                   + wz1 * ut(ip0, jq0, kr0 ) ) )

              ! v - component

              p = FLOOR( ax )
              q = FLOOR( ay + 0.5_wp )
             !r = FLOOR( az )            ! not necessary -> optimization

              ip0 = i - p
              ip1 = ip0 - 1
              jq0 = j - q
              jq1 = jq0 - 1
             !kr0 = k - r                ! not nec. -> opt.
             !kr1 = kr0 - 1              ! not nec. -> opt.

              wx = ax          - REAL( p, wp)
              wy = ay + 0.5_wp - REAL( q, wp)
             !wz = az          - REAL( r, wp)  ! not nec. -> opt.

              wx1 = 1.0_wp - wx
              wy1 = 1.0_wp - wy
             !wz1 = 1.0_wp - wz             ! not nec. -> opt.

              ! trilinear interpolation of the velocity field v
              v_intpl =                                                &
                   &   wx  * (  wy  * (  wz  * vt(ip1, jq1, kr1 )      &
                   &                   + wz1 * vt(ip1, jq1, kr0 ) )    &
                   &          + wy1 * (  wz  * vt(ip1, jq0, kr1 )      &
                   &                   + wz1 * vt(ip1, jq0, kr0 ) ) )  &
                   & + wx1 * (  wy  * (  wz  * vt(ip0, jq1, kr1 )      &
                   &                   + wz1 * vt(ip0, jq1, kr0 ) )    & 
                   &          + wy1 * (  wz  * vt(ip0, jq0, kr1 )      &
                   &                   + wz1 * vt(ip0, jq0, kr0 ) ) )

              ! w - component

             !p = FLOOR( ax )            ! not nec. -> opt.
              q = FLOOR( ay )
              r = FLOOR( az - 0.5_wp )

             !ip0 = i - p                ! not nec. -> opt.
             !ip1 = ip0 - 1              ! not nec. -> opt.
              jq0 = j - q
              jq1 = jq0 - 1
              kr0 = k - r
              kr1 = kr0 - 1

             !wx = ax          - REAL( p, wp)  ! not nec. -> opt.
              wy = ay          - REAL( q, wp)
              wz = az - 0.5_wp - REAL( r, wp)

             !wx1 = 1.0_wp - wx             ! not nec. -> opt.
              wy1 = 1.0_wp - wy
              wz1 = 1.0_wp - wz

              ! trilinear interpolation of the velocity field w
              w_intpl =                                                &
                   &   wx  * (  wy  * (  wz  * wt(ip1, jq1, kr1 )      &
                   &                   + wz1 * wt(ip1, jq1, kr0 ) )    &
                   &          + wy1 * (  wz  * wt(ip1, jq0, kr1 )      &
                   &                   + wz1 * wt(ip1, jq0, kr0 ) ) )  &
                   & + wx1 * (  wy  * (  wz  * wt(ip0, jq1, kr1 )      &
                   &                   + wz1 * wt(ip0, jq1, kr0 ) )    &
                   &          + wy1 * (  wz  * wt(ip0, jq0, kr1 )      &
                   &                   + wz1 * wt(ip0, jq0, kr0 ) ) )

              ! calc. of 2. order dt

              ax = u_intpl * 2.0_wp * dt * eddlon
              ay = v_intpl * 2.0_wp * dt * eddlat
              az = w_intpl * 2.0_wp * dt * eddzeta

              idx(i,j,k,1) = FLOOR( ax )
              idx(i,j,k,2) = FLOOR( ay )
              idx(i,j,k,3) = FLOOR( az )

              wght(i,j,k,1) = ax - REAL( idx(i,j,k,1), wp)
              wght(i,j,k,2) = ay - REAL( idx(i,j,k,2), wp)
              wght(i,j,k,3) = az - REAL( idx(i,j,k,3), wp)

           ELSE 
              ! backtrajectory only in 1. order dt (over 2*dt)
              idx(i,j,k,1) = FLOOR( 2.0_wp * ax )
              idx(i,j,k,2) = FLOOR( 2.0_wp * ay )
              idx(i,j,k,3) = FLOOR( 2.0_wp * az )

              wght(i,j,k,1) = 2.0_wp * ax - REAL( idx(i,j,k,1), wp )
              wght(i,j,k,2) = 2.0_wp * ay - REAL( idx(i,j,k,2), wp )
              wght(i,j,k,3) = 2.0_wp * az - REAL( idx(i,j,k,3), wp )
           END IF

        END DO
     END DO
  END DO

END SUBROUTINE backtraj_trilin_dt2_3tl

!==============================================================================
!==============================================================================

SUBROUTINE  postprocess_backtraj       &
  (idx, weight, ie, je, ke, istart, iend, jstart, jend )

!--------------------------------------------------------------------------
!
! Description:
!   called after backtraj_... this routine checks 
!   - horizontal condition  | CFL_x,y | < 2  (2 for tricubic interpolation;
!     this is needed for reproducibility; the SL-method itself is 
!     of cause not CFL-limited)
!   - prepares backtrajectory for treatment of the upper and lower BC 'value=0'
!     of the 'halo-version' of interpol_SL_tricubic
!
!--------------------------------------------------------------------------

! Declarations:

USE data_parallel,            ONLY : &
    my_cart_id      ! rank of this subdomain in the cartesian communicator
USE environment,              ONLY :  &
    model_abort
USE data_runcontrol, ONLY:  idbg_level, ldebug_dyn

IMPLICIT NONE


INTEGER (KIND=iintegers),   INTENT (IN) ::  &
  ie, je, ke,   & ! x-, y- and z-dimension of the input/output arrays
  istart, iend, & ! computational start and end indices in i-direction
  jstart, jend    ! computational start and end indices in j-direction

INTEGER (KIND=iintegers),   INTENT(INOUT) :: idx(:,:,:,:)
REAL    (KIND=wp)       ,   INTENT(IN)    :: weight  (:,:,:,:)

! Local variables
!----------------

REAL    (KIND=wp),     ALLOCATABLE :: s_halo(:,:,:)

INTEGER (KIND=iintegers)    ::  &
  i, j, k              ! loop indices

INTEGER (KIND=iintegers)    ::  &
  i0, j0, k0   ! intermediate storage of index shifts

INTEGER (KIND=iintegers)    ::  &
  ie_m1, je_m1, ke_p1   ! pre-calculated index limits

LOGICAL        :: l_traj_error=.FALSE.
INTEGER        :: istat
CHARACTER(100) :: yzerrmsg


  IF ( ldebug_dyn .AND. ( idbg_level >= 10 ) .AND. ( my_cart_id == 0 ) ) THEN
    WRITE(*,*) "[Subr. postprocess_backtraj ...]"
  END IF

  ! possible adaptation of the backward trajectory

  ie_m1 = ie-1
  je_m1 = je-1
  ke_p1 = ke+1

  DO k=1, ke
    DO j=jstart, jend
      DO i=istart, iend

        i0  = i - idx(i,j,k,1)

        IF ( ( i0 >= 3 )  .AND. ( i0 <= ie_m1 ) ) THEN
          ! if ( i0-2 >= 1 and i0+1 <= ie ) then
          ! backward trajectory lies inside of the domain; OK
        ELSE 
          ! Courant number limitation |C|<2 due to reproducibility constraint is violated
          ! solution: you have to increase 'nboundlines'
          l_traj_error=.TRUE.  
        END IF

        j0  = j - idx(i,j,k,2)

        IF ( ( j0 >= 3 )  .AND. ( j0 <= je_m1 ) ) THEN
          ! if ( j0-2 >= 1 and j0+1 <= je ) then
          ! backward trajectory lies inside of the domain; OK
        ELSE 
          ! Courant number limitation |C|<2 due to reproducibility constraint is violated
          ! solution: you have to increase 'nboundlines'
          l_traj_error=.TRUE.  
        END IF

        k0  = k - idx(i,j,k,3)

        IF ( ( k0 >= 1 ) .AND. ( k0 <= ke_p1 ) ) THEN
          ! if ( k0-2 >= -1  and  k0+1 <= ke+2 ) then
          ! backward trajectory lies inside of the domain; OK
          ! a possible boundary treatment is done via a 'halo'
        ELSE 
          ! backward trajectory lies outside of the domain 
          ! --> value=0 boundary condition
          ! --> set the backward trajectory above of the model domain
          ! (--> later on there will be simply 0 interpolated)
          ! (possible extension: use 'weight' to better decide, if 
          !  trajectory starts inside of the domain)
          idx(i,j,k,3) = k + 3   ! --> k0 = -3
        END IF

      END DO
    END DO
  END DO

  IF (l_traj_error) THEN
    IF ( my_cart_id == 0 ) THEN
      WRITE(*,*) ""
      WRITE(*,*) "* max CFL for reproducibility is violated."
      WRITE(*,*) "* You may have to increase nboundlines for a proper working SL backward trajectory."
      WRITE(*,*) "* But most probably an instability occured and the error must be found elsewhere."
    END IF

    yzerrmsg = "max CFL for reproducibility is violated (see above)."
    CALL model_abort (my_cart_id, 100, yzerrmsg, 'postprocess_backtraj')
  ENDIF

END SUBROUTINE postprocess_backtraj

!==============================================================================
!==============================================================================

SUBROUTINE interpol_sl_trilin                                             &
            (s_old, s_new, idx, w, dt, ie, je, ke, istart, iend, jstart, jend )

!--------------------------------------------------------------------------
!
! Description:
!   This subroutine interpolates an input field s_old due to the
!   index shifts idx and interpolation weights w to the output field s_new.
!   idx and w are calculated either by the routine
!       backtraj_trilin_dt1_3tl    or
!       backtraj_trilin_dt2_3tl
!
! Input fields:
!   s_old : s at timelevel n-1 (or n)
!   idx, w: index shifts and interpolation weights for backtrajectory
!
! Ouput fields:
!   s_new : s at timelevel n+1
!
! Method:
!   trilinear (=3-dim. linear) Interpolation
!
!--------------------------------------------------------------------------

! Declarations:

INTEGER (KIND=iintegers),   INTENT (IN) ::  &
         ie, je, ke,   & ! x-, y- and z-dimension of the input/output arrays
         istart, iend, & ! computational start and end indices in i-direction
                         !(west-east for lat/lon grid)
         jstart, jend    ! computational start and end indices in j-direction
                         ! (south-north for lat/lon grid)

INTEGER (KIND=iintegers),   INTENT(IN) :: idx(:,:,:,:)
REAL    (KIND=wp)       ,   INTENT(IN) :: w  (:,:,:,:)

REAL    (KIND=wp)       ,   INTENT(IN)  ::  &
  s_old(ie,je,ke), & ! scalar var. to be transported at time-level n-1
  dt                 ! time step

REAL    (KIND=wp)       ,   INTENT(INOUT) ::  &
  s_new(ie,je,ke)    ! scalar variable s at time-level n+1

! Local variables
!----------------

INTEGER (KIND=iintegers)    ::  &
  i, j, k              ! loop indices

INTEGER (KIND=iintegers)    ::  &
  ip0, ip1, jq0, jq1, kr0, kr1   ! intermediate storage of index shifts

REAL (KIND=wp)              ::  &
  ax1, ay1, az1        ! intermediate storage

!------------------------------------------------------------------------------

  DO k=2, ke
    DO j=jstart, jend
      DO i=istart, iend

        ip0 = i-idx(i,j,k,1)
        ip1 = ip0-1
        jq0 = j-idx(i,j,k,2)
        jq1 = jq0-1
        kr0 = k-idx(i,j,k,3)
        kr1 = kr0-1

        IF  (  ( ip1 >= 1 ) .AND. ( ip0 <= ie ) .AND.   &
          &    ( jq1 >= 1 ) .AND. ( jq0 <= je ) .AND.   &
          &    ( kr1 >= 1 )                   ) THEN

          IF ( kr0 <= ke ) THEN
            ! back trajectory starts inside of the domain

            ax1 = 1.0_wp - w(i,j,k,1)
            ay1 = 1.0_wp - w(i,j,k,2)
            az1 = 1.0_wp - w(i,j,k,3)

            ! trilinear interpolation of the transported field
            s_new(i,j,k) =                                                   &
              w(i,j,k,1) *( w(i,j,k,2) *( w(i,j,k,3)* s_old(ip1,jq1,kr1)     &
                                             + az1  * s_old(ip1,jq1,kr0) )   &
                                 + ay1 *( w(i,j,k,3)* s_old(ip1,jq0,kr1)     &
                                             + az1  * s_old(ip1,jq0,kr0) ) ) &
               + ax1     *( w(i,j,k,2) *( w(i,j,k,3)* s_old(ip0,jq1,kr1)     &
                                             + az1  * s_old(ip0,jq1,kr0) )   &
                                 + ay1 *( w(i,j,k,3)* s_old(ip0,jq0,kr1)     &
                                              + az1 * s_old(ip0,jq0,kr0) ) )

            ! Tendency ds/dt:
            ! sten(i,j,k) = ( s_new(i,j,k) - s_old(i,j,k) ) / 2.0 / dt

          ELSE IF ( kr0 == ke+1 ) THEN
            !  back trajectory starts in the bottom layer

            ax1 = 1.0_wp - w(i,j,k,1)
            ay1 = 1.0_wp - w(i,j,k,2)

            ! trilinear interpolation of the transported field
            s_new(i,j,k) =                                       &
                w(i,j,k,1) *( w(i,j,k,2) * s_old(ip1,jq1,kr1 )   &
                                   + ay1 * s_old(ip1,jq0,kr1 ) ) &
                 + ax1     *( w(i,j,k,2) * s_old(ip0,jq1,kr1 )   &
                                   + ay1 * s_old(ip0,jq0,kr1 ) )

            ! Tendency ds/dt:
            ! sten(i,j,k) = ( s_new(i,j,k) - s_old(i,j,k) ) / 2.0 / dt

          ELSE
            ! back trajectory starts outside of the domain
            ! b.c.: no transport of s from outside

            s_new(i,j,k) = 0.0_wp

            !sten(i,j,k) = -s_old(i,j,k)    ! ???
          END IF
        ELSE
          ! back trajectory starts outside of the domain
          ! b.c.: no transport of s from outside

          s_new(i,j,k) = 0.0_wp

          !sten(i,j,k) = -s_old(i,j,k)    ! ???
        ENDIF

      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE interpol_sl_trilin

!==============================================================================
!==============================================================================

SUBROUTINE interpol_sl_tricubic     &
  (s_old, s_new, idx, weight, dt, ie, je, ke, istart, iend, jstart, jend )

!--------------------------------------------------------------------------
!
! Description:
!   This subroutine interpolates an input field s_old due to the
!   index shifts 'idx' and interpolation 'weights' to the output field s_new.
!   'idx' and 'weight' are calculated either by the routine
!       backtraj_trilin_dt1_3tl    or
!       backtraj_trilin_dt2_3tl
!
! Input fields:
!   s_old : s at timelevel n-1 (or n)
!   idx, weight: index shifts and interpolation weights for backtrajectory
!
! Ouput fields:
!   s_new : s at timelevel n+1
!
! Method:
!   tricubic (=3-dim. cubic) Interpolation
!
!   a clipping of over/undershoots is possible 
!   (l_avoid_under_overshoots=.TRUE.), but not recommended
!
!--------------------------------------------------------------------------

! Declarations:

USE data_parallel,            ONLY : &
  my_cart_id      ! rank of this subdomain in the cartesian communicator

USE environment,              ONLY :  &
  model_abort

USE data_runcontrol, ONLY:  idbg_level, ldebug_dyn

IMPLICIT NONE


INTEGER (KIND=iintegers),   INTENT (IN) ::  &
  ie, je, ke,   & ! x-, y- and z-dimension of the input/output arrays
  istart, iend, & ! computational start and end indices in i-direction
                                ! (west-east for lat/lon grid)
  jstart, jend    ! computational start and end indices in j-direction
! (south-north for lat/lon grid)

INTEGER (KIND=iintegers),   INTENT(IN) :: idx   (:,:,:,:)
REAL    (KIND=wp)       ,   INTENT(IN) :: weight(:,:,:,:)

REAL    (KIND=wp)       ,   INTENT(IN)  ::  &
  s_old(ie,je,ke), & ! scalar var. to be transported at time-level n-1
  dt                 ! time step

REAL    (KIND=wp)       ,   INTENT(INOUT) ::  &
  s_new(ie,je,ke)    ! scalar variable s at time-level n+1

! Local variables
!----------------

REAL    (KIND=wp),     ALLOCATABLE :: s_halo(:,:,:)

INTEGER (KIND=iintegers)    ::  &
  i, j, k              ! loop indices

INTEGER (KIND=iintegers)    ::  &
                                ! intermediate storage of index shifts
  im2, im1, i0, ip1, jm2, jm1, j0, jp1, km2, km1, k0, kp1   

REAL (KIND=wp)              ::  &
  wx, wy, wz

REAL (KIND=wp)              ::  &
  mini, maxi

REAL (KIND=wp)       :: ap1,  a0,  am1,  am2

REAL (KIND=wp)       :: ap1x, a0x, am1x, am2x
REAL (KIND=wp)       :: ap1y, a0y, am1y, am2y
REAL (KIND=wp)       :: ap1z, a0z, am1z, am2z
REAL (KIND=wp)       :: x

REAL (KIND=wp)       :: hilf1, hilf2, hilf3, hilf4

INTEGER              :: istat

LOGICAL, PARAMETER   :: l_avoid_under_overshoots = .FALSE.

CHARACTER(80)        :: yzerrmsg

! Statement functions
!--------------------

ap1(x) =  1.0_wp/6.0_wp              * x * (x+1.0_wp) * (x+2.0_wp)
a0 (x) =   -0.5_wp      * (x-1.0_wp)     * (x+1.0_wp) * (x+2.0_wp)
am1(x) =    0.5_wp      * (x-1.0_wp) * x              * (x+2.0_wp)
am2(x) = -1.0_wp/6.0_wp * (x-1.0_wp) * x * (x+1.0_wp)

!------------------------------------------------------------------------------

IF ( ldebug_dyn .AND. ( idbg_level >= 10 ) .AND. ( my_cart_id==0 ) ) THEN
  WRITE(*,*) "[Subr. interpol_sl_tricubic ...]"
END IF

ALLOCATE( s_halo( 1:ie, 1:je, -5:ke+2 ), STAT=istat )
IF ( istat /= 0 ) THEN
  yzerrmsg = "allocation of s_halo failed"
  CALL model_abort (my_cart_id, 100, yzerrmsg, 'interpol_sl_tricubic')
END IF


!CALL  FTRACE_REGION_BEGIN("SL-Intp_init_halo")

DO k=1, ke
  DO j=1, je
    DO i=1, ie
      s_halo(i,j,k) = s_old(i,j,k)
    END DO
  END DO
END DO
! alternativ:   s_halo(1:ie,1:je,1:ke) = s_old(1:ie,1:je,1:ke)


DO j=1, je
  DO i=1, ie
    ! halo values for a possible boundary treatment
    s_halo(i,j, 0)   = s_old(i,j,1)
    s_halo(i,j,-1)   = s_old(i,j,1)

    s_halo(i,j,ke+1) = s_old(i,j,ke)
    s_halo(i,j,ke+2) = s_old(i,j,ke)

    ! artificial halo values for the boundary condition 'value=0'
    s_halo(i,j,-2) = 0.0_wp
    s_halo(i,j,-3) = 0.0_wp
    s_halo(i,j,-4) = 0.0_wp
    s_halo(i,j,-5) = 0.0_wp

  END DO
END DO


DO k=1, ke
  DO j=jstart, jend
    DO i=istart, iend

      i0  = i - idx(i,j,k,1)
      im1 = i0 - 1
      im2 = i0 - 2
      ip1 = i0 + 1
      j0  = j - idx(i,j,k,2)
      jm1 = j0 - 1
      jm2 = j0 - 2
      jp1 = j0 + 1
      k0  = k - idx(i,j,k,3)
      km1 = k0 - 1
      km2 = k0 - 2
      kp1 = k0 + 1

      ! back trajectory always starts inside of the domain (extended by the halo)

      wx = - weight(i,j,k,1)
      wy = - weight(i,j,k,2)
      wz = - weight(i,j,k,3)

      ap1z = ap1(wz)
      a0z  = a0 (wz)
      am1z = am1(wz)
      am2z = am2(wz)

      !ap1y = ap1(wy)
      !a0y  = a0 (wy)
      !am1y = am1(wy)
      !am2y = am2(wy)

      !ap1x = ap1(wx)
      !a0x  = a0 (wx)
      !am1x = am1(wx)
      !am2x = am2(wx)


      ! tricubic interpolation of the transported field
      s_new(i,j,k) =                                                 &
        &   ap1(wx) * ( ap1(wy) * ( ap1z * s_halo(ip1, jp1, kp1 )    &
        &                         + a0z  * s_halo(ip1, jp1, k0  )    &
        &                         + am1z * s_halo(ip1, jp1, km1 )    &
        &                         + am2z * s_halo(ip1, jp1, km2 ) )  &
        &             + a0 (wy) * ( ap1z * s_halo(ip1, j0 , kp1 )    &
        &                         + a0z  * s_halo(ip1, j0 , k0  )    &
        &                         + am1z * s_halo(ip1, j0 , km1 )    &
        &                         + am2z * s_halo(ip1, j0 , km2 ) )  &
        &             + am1(wy) * ( ap1z * s_halo(ip1, jm1, kp1 )    &
        &                         + a0z  * s_halo(ip1, jm1, k0  )    &
        &                         + am1z * s_halo(ip1, jm1, km1 )    &
        &                         + am2z * s_halo(ip1, jm1, km2 ) )  &
        &             + am2(wy) * ( ap1z * s_halo(ip1, jm2, kp1 )    &
        &                         + a0z  * s_halo(ip1, jm2, k0  )    &
        &                         + am1z * s_halo(ip1, jm2, km1 )    &
        &                         + am2z * s_halo(ip1, jm2, km2 ) ) )

      s_new(i,j,k) = s_new(i,j,k) +                                  &
        &   a0 (wx) * ( ap1(wy) * ( ap1z * s_halo(i0 , jp1, kp1 )    &
        &                         + a0z  * s_halo(i0 , jp1, k0  )    &
        &                         + am1z * s_halo(i0 , jp1, km1 )    &
        &                         + am2z * s_halo(i0 , jp1, km2 ) )  &
        &             + a0 (wy) * ( ap1z * s_halo(i0 , j0 , kp1 )    &
        &                         + a0z  * s_halo(i0 , j0 , k0  )    &
        &                         + am1z * s_halo(i0 , j0 , km1 )    &
        &                         + am2z * s_halo(i0 , j0 , km2 ) )  &
        &             + am1(wy) * ( ap1z * s_halo(i0 , jm1, kp1 )    &
        &                         + a0z  * s_halo(i0 , jm1, k0  )    &
        &                         + am1z * s_halo(i0 , jm1, km1 )    &
        &                         + am2z * s_halo(i0 , jm1, km2 ) )  &
        &             + am2(wy) * ( ap1z * s_halo(i0 , jm2, kp1 )    &
        &                         + a0z  * s_halo(i0 , jm2, k0  )    &
        &                         + am1z * s_halo(i0 , jm2, km1 )    &
        &                         + am2z * s_halo(i0 , jm2, km2 ) ) )

      s_new(i,j,k) = s_new(i,j,k) +                                  &
        &   am1(wx) * ( ap1(wy) * ( ap1z * s_halo(im1, jp1, kp1 )    &
        &                         + a0z  * s_halo(im1, jp1, k0  )    &
        &                         + am1z * s_halo(im1, jp1, km1 )    &
        &                         + am2z * s_halo(im1, jp1, km2 ) )  &
        &             + a0 (wy) * ( ap1z * s_halo(im1, j0 , kp1 )    &
        &                         + a0z  * s_halo(im1, j0 , k0  )    &
        &                         + am1z * s_halo(im1, j0 , km1 )    &
        &                         + am2z * s_halo(im1, j0 , km2 ) )  &
        &             + am1(wy) * ( ap1z * s_halo(im1, jm1, kp1 )    &
        &                         + a0z  * s_halo(im1, jm1, k0  )    &
        &                         + am1z * s_halo(im1, jm1, km1 )    &
        &                         + am2z * s_halo(im1, jm1, km2 ) )  &
        &             + am2(wy) * ( ap1z * s_halo(im1, jm2, kp1 )    &
        &                         + a0z  * s_halo(im1, jm2, k0  )    &
        &                         + am1z * s_halo(im1, jm2, km1 )    &
        &                         + am2z * s_halo(im1, jm2, km2 ) ) )

      s_new(i,j,k) = s_new(i,j,k) +                                  &
        &   am2(wx) * ( ap1(wy) * ( ap1z * s_halo(im2, jp1, kp1 )    &
        &                         + a0z  * s_halo(im2, jp1, k0  )    &
        &                         + am1z * s_halo(im2, jp1, km1 )    &
        &                         + am2z * s_halo(im2, jp1, km2 ) )  &
        &             + a0 (wy) * ( ap1z * s_halo(im2, j0 , kp1 )    &
        &                         + a0z  * s_halo(im2, j0 , k0  )    &
        &                         + am1z * s_halo(im2, j0 , km1 )    &
        &                         + am2z * s_halo(im2, j0 , km2 ) )  &
        &             + am1(wy) * ( ap1z * s_halo(im2, jm1, kp1 )    &
        &                         + a0z  * s_halo(im2, jm1, k0  )    &
        &                         + am1z * s_halo(im2, jm1, km1 )    &
        &                         + am2z * s_halo(im2, jm1, km2 ) )  &
        &             + am2(wy) * ( ap1z * s_halo(im2, jm2, kp1 )    &
        &                         + a0z  * s_halo(im2, jm2, k0  )    &
        &                         + am1z * s_halo(im2, jm2, km1 )    &
        &                         + am2z * s_halo(im2, jm2, km2 ) ) )

      ! Tendency ds/dt:
      ! sten(i,j,k) = ( s_new(i,j,k) - s_halo(i,j,k) ) / 2.0 / dt

    ENDDO
  ENDDO
ENDDO

!CALL  FTRACE_REGION_END("SL-Interpolation")


IF ( l_avoid_under_overshoots ) THEN
  ! Clipping; avoids over-/undershootings

  DO k=1, ke
    DO j=jstart, jend
      DO i=istart, iend

        maxi = MAX ( s_old( i0 , j0 , k0  ), &
          &          s_old( i0 , j0 , km1 ), &
          &          s_old( i0 , jm1, k0  ), &
          &          s_old( i0 , jm1, km1 ), &
          &          s_old( im1, j0 , k0  ), &
          &          s_old( im1, j0 , km1 ), &
          &          s_old( im1, jm1, k0  ), &
          &          s_old( im1, jm1, km1 ) )

        mini = MIN ( s_old( i0 , j0 , k0  ), &
          &          s_old( i0 , j0 , km1 ), &
          &          s_old( i0 , jm1, k0  ), &
          &          s_old( i0 , jm1, km1 ), &
          &          s_old( im1, j0 , k0  ), &
          &          s_old( im1, j0 , km1 ), &
          &          s_old( im1, jm1, k0  ), &
          &          s_old( im1, jm1, km1 ) )

        IF ( s_new(i,j,k) < mini ) THEN
          s_new(i,j,k) = mini
        END IF

        IF ( s_new(i,j,k) > maxi ) THEN
          s_new(i,j,k) = maxi
        END IF
      ENDDO
    ENDDO
  ENDDO

END IF

DEALLOCATE ( s_halo )

END SUBROUTINE interpol_sl_tricubic


!==============================================================================
!==============================================================================


SUBROUTINE remove_negative_values &
  ( field, ie, je, ke, istart, iend, jstart, jend, &
  i_clipping_type, diffus_type )

!--------------------------------------------------------------------------
!
! Description:
!   This subroutine removes negative values from 'field'
!
! Input fields:
!   field(:,:,:)
!
! Ouput fields:
!   field(:,:,:)
!
! Method:
!   several methods are available:
!   i_clipping_type=0 : nothing is done here
!   i_clipping_type=1 : simple clipping of negative values 
!                      (in general this is a tremendous mass source!)
!   i_clipping_type=3 : 'selective filling diffusion' + clipping afterwards
!   i_clipping_type=4 : multiplicative filling (Rood, 1987)
!
!   diffus_type : sequence of 1D-diffusion ("xyz", "zyx" )
!
!--------------------------------------------------------------------------

! Declarations:

USE data_parallel,            ONLY : &
  my_cart_id      ! rank of this subdomain in the cartesian communicator

USE data_runcontrol,          ONLY:  idbg_level, ldebug_dyn

USE numeric_utilities_rk,     ONLY :  &
  multiplicative_filling,   & !
  multiplicative_filling_DDI  !

USE environment,              ONLY :  &
  model_abort

IMPLICIT NONE


INTEGER (KIND=iintegers),   INTENT (IN) ::  &
  ie, je, ke,   & ! x-, y- and z-dimension of the input/output arrays
  istart, iend, & ! computational start and end indices in i-direction
  jstart, jend    ! computational start and end indices in j-direction

REAL    (KIND=wp)       ,   INTENT(INOUT) ::  &
  field( 1:ie, 1:je, 1:ke)   ! scalar variable

INTEGER, INTENT(IN)       :: i_clipping_type
CHARACTER(*), INTENT(IN)  :: diffus_type   
       ! Type of artificial diffusion: ("3D", "xyz", "zyx")


! Local variables
!----------------

INTEGER (KIND=iintegers)    ::  &
  i, j, k              ! loop indices

LOGICAL, ALLOCATABLE        ::  clip_flag(:,:,:)

INTEGER (KIND=iintegers)    :: istat
INTEGER (KIND=iintegers)    :: count_clip_flags

LOGICAL        :: l_reproduce = .TRUE.
CHARACTER(80)  :: yzerrmsg

REAL (KIND=wp)     :: eps_diff_limit = 0.0_wp
  ! perform 'selective filling diffusion' only if value < eps_diff_limit
  ! recommendation: eps_diff_limit = 0 or slightly negative


IF ( ldebug_dyn .AND. ( idbg_level >= 10 ) .AND. ( my_cart_id==0 ) ) THEN
  WRITE(*,*) "[Subr. remove_negative_values ...]"
  WRITE(*,*) "  i_clipping_type=", i_clipping_type
END IF


SELECT CASE ( i_clipping_type )
CASE ( 0 )
  ! nothing to do

CASE ( 1 )

  ! simple clipping of negative values.
  ! This is in general a tremendous mass source!

  WHERE ( field(:,:,:) < 0.0_wp)
    field(:,:,:) = 0.0_wp
  END WHERE

CASE ( 3 )

  ! Selective Filling Diffusion (SFD)

  ALLOCATE( clip_flag(1:ie, 1:je, 1:ke), stat=istat)
  IF ( istat /= 0 ) THEN
    yzerrmsg = "allocation of clip_flag failed"
    CALL model_abort (my_cart_id, 100, yzerrmsg, 'remove_negative_values')
  END IF

  clip_flag(:,:,:) = .FALSE.

  DO k=1, ke
    DO j=jstart-2, jend+2
      DO i=istart-2, iend+2
        ! mark cells, where diffusion will take place
        ! (2 more rows and columns --> no exchange necessary between 
        !  the different diffusion directions)

        IF ( field(i,j,k) < eps_diff_limit ) THEN
          clip_flag(i  ,j  ,k  ) = .TRUE.
        END IF

      ENDDO
    ENDDO
  ENDDO

  IF ( diffus_type == "3D" ) THEN
    CALL diffusion( field, clip_flag, ie, je, ke )

  ELSE IF ( diffus_type == "xyz" ) THEN
    CALL diffusion_x( field, clip_flag, ie, je, ke )
    CALL diffusion_y( field, clip_flag, ie, je, ke )
    CALL diffusion_z( field, clip_flag, ie, je, ke )

  ELSE IF ( diffus_type == "zyx" ) THEN
    CALL diffusion_z( field, clip_flag, ie, je, ke )
    CALL diffusion_y( field, clip_flag, ie, je, ke )
    CALL diffusion_x( field, clip_flag, ie, je, ke )

  ELSE
    yzerrmsg = "false value in diffus_type"
    CALL model_abort (my_cart_id, 100, yzerrmsg, 'remove_negative_values')
  END IF

  ! .. and clipping:
  WHERE ( field(:,:,:) < 0.0_wp ) 
    field(:,:,:) = 0.0_wp
  END WHERE

  !IF ( .TRUE. ) THEN
  !  ! some statistics:
  !  count_clip_flags = 0
  !  DO k=1, ke
  !    DO j=jstart-1, jend+1
  !      DO i=istart-1, iend+1
  !        IF ( clip_flag(i,j,k) ) THEN
  !          count_clip_flags = count_clip_flags + 1
  !        END IF
  !      END DO
  !    END DO
  !  END DO

  !  WRITE(*,'(A,I2,A,I9,A,F9.6,A)') "[remove_negative_values] (mycart_id=", &
  !    my_cart_id, ") count_clip_flags = ", count_clip_flags,  &
  !    " ~ ", DBLE(count_clip_flags) / DBLE( (iend-istart+3)*(jend-jstart+3)*ke ) * 100.0, "% GPe"
  !END IF

  DEALLOCATE( clip_flag )

CASE (  4 )

  ! multiplicative filling (Rood 1987)

  IF ( l_reproduce ) THEN
    CALL multiplicative_filling_DDI( field(:,:,:), 6 )
  ELSE
    CALL multiplicative_filling( field(:,:,:), 6 )
  END IF

CASE DEFAULT

  yzerrmsg = "false value for i_clipping_type"
  CALL model_abort (my_cart_id, 100, yzerrmsg, 'remove_negative_values')

END SELECT


END SUBROUTINE remove_negative_values


!==============================================================================
!==============================================================================


SUBROUTINE diffusion( phi, clip_flag, ie, je, ke ) 

! 3D-Diffusion of a field phi(:,:,:)
! with an (artificial) diffusion coefficient alpha.
! phi is considered as a density ( [x/m^3] ).
! The diffusion is only performed on points with clip_flag(i,j,k) == .TRUE.
!   

USE data_constants, ONLY: r_earth, pi
USE data_runcontrol, ONLY:  idbg_level, ldebug_dyn
USE data_fields,    ONLY: crlat
USE grid_metrics_utilities,    ONLY: sqrtg_r_s
USE data_modelconfig, ONLY: dlon, dlat
USE data_parallel,            ONLY : &
    my_cart_id      ! rank of this subdomain in the cartesian communicator

INTEGER (KIND=iintegers),   INTENT (IN) ::  &
  ie, je, ke     ! x-, y- and z-dimension of the input/output arrays

REAL (KIND=wp),     INTENT(inout) :: phi(1:ie, 1:je, 1:ke)
LOGICAL, INTENT(in)            :: clip_flag(1:ie, 1:je, 1:ke)

REAL (KIND=wp),     ALLOCATABLE :: fx(:,:,:), fy(:,:,:), fz(:,:,:)

INTEGER            :: i, j, k
REAL (KIND=wp)     :: dx, dy, dz, vol
REAL (KIND=wp)     :: alpha      ! dim.less diffusion coefficient
REAL (KIND=wp),     PARAMETER        :: dzeta = 1.0_wp
INTEGER  :: istat

IF ( ldebug_dyn .AND. ( idbg_level >= 10 ) .AND. ( my_cart_id==0 ) ) THEN
  WRITE(*,*) "[Subr. diffusion] ..."
END IF

ALLOCATE( fx(1:ie+1, 1:je, 1:ke), stat=istat)
ALLOCATE( fy(1:ie, 1:je+1, 1:ke), stat=istat)
ALLOCATE( fz(1:ie, 1:je, 1:ke+1), stat=istat)

! 1.) Flux calculation where necessary

alpha = 0.25_wp / 3.0_wp

fx(:,:,:) = 0.0_wp
fy(:,:,:) = 0.0_wp
fz(:,:,:) = 0.0_wp

dy = dlat * R_earth * pi/180.0_wp

DO k=1, ke-1
  !DO j=jstart-1, jend+1
  DO j=1, je-1
    dx = dlon * R_earth * pi/180.0_wp * crlat(j,1)

    DO i=1, ie-1

        dz = dzeta / sqrtg_r_s(i,j,k)
        ! look only to the cells on the right, below and behind

        IF ( clip_flag(i,j,k) .OR. clip_flag(i+1,j,k) ) THEN
          ! diffusion flux * dt
          ! Index convention for flux i <--> i+1/2
          ! (artificial diffus.coeff: K = alpha * dx * dx / dt )
          fx(i,j,k) = - alpha * dx * (dy*dz) * ( phi(i+1,j,k) - phi(i,j,k) )
        END IF

        IF ( clip_flag(i,j,k) .OR. clip_flag(i,j+1,k) ) THEN
          ! diffusion flux * dt
          ! Index convention for flux j <--> j+1/2
          ! (artificial diffus.coeff: K = alpha * dy * dy / dt )
          fy(i,j,k) = - alpha * dy * (dx*dz) * ( phi(i,j+1,k) - phi(i,j,k) )
        END IF

        IF ( clip_flag(i,j,k) .OR. clip_flag(i,j,k+1) ) THEN
          ! diffusion flux * dt
          ! Index convention for flux k <--> k+1/2
          ! (artificial diffus.coeff: K = alpha * dz * dz / dt )
          fz(i,j,k) = - alpha * dz * (dx*dy) * ( phi(i,j,k+1) - phi(i,j,k) )
        END IF

    END DO
  END DO
END DO


! 2.) calculate divergence

DO k=2, ke-1
  DO j=2, je-1
    dx = dlon * R_earth * pi/180.0_wp * crlat(j,1)
    DO i=2, ie-1

      IF (       clip_flag(i  ,j  ,k  )   &
        &   .OR. clip_flag(i+1,j  ,k  )   &
        &   .OR. clip_flag(i-1,j  ,k  )   &
        &   .OR. clip_flag(i  ,j+1,k  )   &
        &   .OR. clip_flag(i  ,j-1,k  )   &
        &   .OR. clip_flag(i  ,j  ,k+1)   &
        &   .OR. clip_flag(i  ,j  ,k-1) ) THEN

        dz = dzeta / sqrtg_r_s(i,j,k)

        vol = dx * dy * dz

        phi(i,j,k) = phi(i,j,k) -       &
          ( fx(i,j,k) - fx(i-1,j,k)       &
          + fy(i,j,k) - fy(i,j-1,k)       &
          + fz(i,j,k) - fz(i,j,k-1) ) / vol
        ! this construction guarantees integral conservation of phi

      END IF

    END DO
  END DO
END DO

DEALLOCATE( fx, fy, fz )

END SUBROUTINE diffusion


SUBROUTINE diffusion_x( phi, clip_flag, ie, je, ke ) 

! Diffusion in x-direction of a field phi(:,:,:)
! with an (artificial) diffusion coefficient alpha.
! phi is considered as a density ( [x/m^3] ).
! The diffusion is only performed on points with clip_flag(i,j,k) == .TRUE.
!   

USE data_constants,   ONLY: r_earth, pi
USE data_runcontrol,  ONLY: idbg_level, ldebug_dyn
USE data_fields,      ONLY: crlat
USE grid_metrics_utilities,    ONLY: sqrtg_r_s, sqrtg_r_u
USE data_modelconfig, ONLY: dlon, dlat
USE data_parallel,    ONLY :    &
    my_cart_id      ! rank of this subdomain in the cartesian communicator

INTEGER (KIND=iintegers),   INTENT (IN) ::  &
  ie, je, ke    ! x-, y- and z-dimension of the input/output arrays

REAL (KIND=wp),  INTENT(inout) :: phi(1:ie, 1:je, 1:ke)
LOGICAL, INTENT(in)            :: clip_flag(1:ie, 1:je, 1:ke)

REAL (KIND=wp),     ALLOCATABLE :: fx(:,:,:)

INTEGER            :: i, j, k
REAL (KIND=wp)     :: dx, dy, dz, dx_aequ, dx_aequ_dy, dxdy, vol
REAL (KIND=wp)     :: alpha      ! dim.less diffusion coefficient
REAL (KIND=wp),     PARAMETER    :: dzeta = 1.0_wp
INTEGER  :: istat

IF ( ldebug_dyn .AND. ( idbg_level >= 10 ) .AND. ( my_cart_id==0 ) ) THEN
  WRITE(*,*) "[Subr. diffusion_x] ..."
END IF


ALLOCATE( fx(1:ie, 1:je, 1:ke), stat=istat)

! 1.) Flux calculation where necessary

alpha = 0.25_wp   ! dim.less diffusion coefficient

fx(:,:,:) = 0.0_wp

dy         = dlat * R_earth * pi/180.0_wp
dx_aequ    = dlon * R_earth * pi/180.0_wp
dx_aequ_dy = dx_aequ * dy 

DO k=1, ke
  DO j=1, je
    dxdy = dx_aequ_dy * crlat(j,1)

    DO i=1, ie-1

      IF ( clip_flag(i,j,k) .OR. clip_flag(i+1,j,k) ) THEN

        dz = dzeta / sqrtg_r_u(i,j,k)

        ! diffusion flux * dt
        ! (Index convention for flux:  i <--> i+1/2)
        ! (artificial diffus.coeff: K = alpha * dx * dx / dt )
        fx(i,j,k) = - alpha * dxdy * dz * ( phi(i+1,j,k) - phi(i,j,k) )

      END IF

    END DO
  END DO
END DO

! 2.) calculate divergence

DO k=1, ke
  DO j=1, je
    dxdy = dx_aequ_dy * crlat(j,1)
    DO i=2, ie-1

      IF ( clip_flag(i-1,j,k) .OR. clip_flag(i,j,k) .OR. clip_flag(i+1,j,k) ) THEN
        dz = dzeta / sqrtg_r_s(i,j,k)
        vol = dxdy * dz
        phi(i,j,k) = phi(i,j,k) - ( fx(i,j,k) - fx(i-1,j,k) ) / vol
        ! this construction guarantees integral conservation of phi
      END IF

    END DO
  END DO
END DO

DEALLOCATE( fx )

END SUBROUTINE diffusion_x



SUBROUTINE diffusion_y( phi, clip_flag, ie, je, ke ) 

! 1D-Diffusion of a field phi(:,:,:)
! with an (artificial) diffusion coefficient.
! phi is considered as a density ( [x/m^3] ).
! The diffusion is only performed on points with clip_flag(i,j,k) == .TRUE.
!   

USE data_constants,   ONLY: r_earth, pi
USE data_runcontrol,  ONLY:  idbg_level, ldebug_dyn
USE data_fields,      ONLY: crlat
USE grid_metrics_utilities,    ONLY: sqrtg_r_s, sqrtg_r_v
USE data_modelconfig, ONLY: dlon, dlat
USE data_parallel,    ONLY:   &
    my_cart_id      ! rank of this subdomain in the cartesian communicator

INTEGER (KIND=iintegers),   INTENT (IN) ::  &
  ie, je, ke    ! x-, y- and z-dimension of the input/output arrays

REAL (KIND=wp),     INTENT(inout) :: phi(1:ie, 1:je, 1:ke)
LOGICAL, INTENT(in)            :: clip_flag(1:ie, 1:je, 1:ke)

REAL (KIND=wp),     ALLOCATABLE :: fy(:,:,:)

INTEGER            :: i, j, k
REAL (KIND=wp)     :: dx, dy, dz, dx_aequ, dx_aequ_dy, dxdy, vol
REAL (KIND=wp)     :: alpha      ! dim.less diffusion coefficient
REAL (KIND=wp),     PARAMETER        :: dzeta = 1.0_wp
INTEGER  :: istat

IF ( ldebug_dyn .AND. ( idbg_level >= 10 ) .AND. ( my_cart_id==0 ) ) THEN
  WRITE(*,*) "[Subr. diffusion_y] ..."
END IF

ALLOCATE( fy(1:ie, 1:je, 1:ke), stat=istat)

! 1.) Flux calculation where necessary

alpha = 0.25_wp

fy(:,:,:) = 0.0_wp

dy         = dlat * R_earth * pi/180.0_wp
dx_aequ    = dlon * R_earth * pi/180.0_wp
dx_aequ_dy = dx_aequ * dy

DO k=1, ke
  DO j=1, je-1
    dxdy = dx_aequ_dy * crlat(j,1)

    DO i=1,ie

      IF ( clip_flag(i,j,k) .OR. clip_flag(i,j+1,k) ) THEN

        dz = dzeta / sqrtg_r_v(i,j,k)

          ! diffusion flux * dt
          ! Index convention for flux j <--> j+1/2
          ! (artificial diffus.coeff: K = alpha * dy * dy / dt )
          fy(i,j,k) = - alpha * dxdy *dz * ( phi(i,j+1,k) - phi(i,j,k) )

      END IF

    END DO
  END DO
END DO

! 2.) calculate divergence

DO k=1, ke
  DO j=2, je-1
    dxdy = dx_aequ_dy * crlat(j,1)
    DO i=1, ie

      IF ( clip_flag(i,j-1,k) .OR. clip_flag(i,j,k) .OR. clip_flag(i,j+1,k) ) THEN
        dz = dzeta / sqrtg_r_s(i,j,k)
        vol = dxdy * dz
        phi(i,j,k) = phi(i,j,k) - ( fy(i,j,k) - fy(i,j-1,k) ) / vol
        ! this construction guarantees integral conservation of phi
      END IF

    END DO
  END DO
END DO

DEALLOCATE( fy )

END SUBROUTINE diffusion_y


SUBROUTINE diffusion_z( phi, clip_flag, ie, je, ke ) 

! 1D-Diffusion of a field phi(:,:,:)
! with an (artificial) diffusion coefficient alpha.
! phi is considered as a density ( [x/m^3] ).
! The diffusion is only performed on points with clip_flag(i,j,k) == .TRUE.
!   

USE data_constants,   ONLY: r_earth, pi
USE data_runcontrol,  ONLY: idbg_level, ldebug_dyn
USE data_fields,      ONLY: crlat
USE grid_metrics_utilities,    ONLY: sqrtg_r_s, sqrtg_r_w
USE data_modelconfig, ONLY: dlon, dlat
USE data_parallel,    ONLY: &
    my_cart_id      ! rank of this subdomain in the cartesian communicator


INTEGER (KIND=iintegers),   INTENT (IN) ::  &
  ie, je, ke    ! x-, y- and z-dimension of the input/output arrays

REAL (KIND=wp),     INTENT(inout) :: phi(1:ie, 1:je, 1:ke)
LOGICAL, INTENT(in)            :: clip_flag(1:ie, 1:je, 1:ke)

REAL (KIND=wp),     ALLOCATABLE :: fz(:,:,:)

INTEGER            :: i, j, k
REAL (KIND=wp)     :: dx, dy, dz, dx_aequ, dx_aequ_dy, dxdy, vol
REAL (KIND=wp)     :: alpha      ! dim.less diffusion coefficient
REAL (KIND=wp),     PARAMETER        :: dzeta = 1.0_wp
INTEGER  :: istat

IF ( ldebug_dyn .AND. ( idbg_level >= 10 ) .AND. ( my_cart_id==0 ) ) THEN
  WRITE(*,*) "[Subr. diffusion_z] ..."
END IF

ALLOCATE( fz(1:ie, 1:je, 1:ke), stat=istat)

! 1.) Flux calculation where necessary

alpha = 0.25_wp

fz(:,:,:) = 0.0_wp

dy         = dlat * R_earth * pi/180.0_wp
dx_aequ    = dlon * R_earth * pi/180.0_wp
dx_aequ_dy = dx_aequ * dy

DO k=1, ke-1
  DO j=1, je
    dxdy = dx_aequ_dy * crlat(j,1)

    DO i=1, ie

      IF ( clip_flag(i,j,k) .OR. clip_flag(i,j,k+1) ) THEN

        dz = dzeta / sqrtg_r_w(i,j,k)

        ! diffusion flux * dt
        ! Index convention for flux k <--> k+1/2
        ! (artificial diffus.coeff: K = alpha * dz * dz / dt )
        fz(i,j,k) = - alpha * dz * dxdy * ( phi(i,j,k+1) - phi(i,j,k) )
      END IF

    END DO
  END DO
END DO

! 2.) calculate divergence

DO k=2, ke-1
  DO j=1, je
    dxdy = dx_aequ_dy * crlat(j,1)
    DO i=1, ie

      IF ( clip_flag(i,j,k-1) .OR. clip_flag(i,j,k) .OR. clip_flag(i,j,k+1) ) THEN
        dz = dzeta / sqrtg_r_s(i,j,k)
        vol = dxdy * dz
        phi(i,j,k) = phi(i,j,k) - ( fz(i,j,k) - fz(i,j,k-1) ) / vol
        ! this construction guarantees integral conservation of phi
      END IF

    END DO
  END DO
END DO

!  k=1
!  DO j=1, je
!    dxdy = dx_aequ_dy * crlat(j,1)
!    DO i=1, ie

!      IF ( clip_flag(i,j,k) .OR. clip_flag(i,j,k+1) ) THEN
!        dz = dzeta / sqrtg_r_s(i,j,k)
!        vol = dxdy * dz
!        phi(i,j,k) = phi(i,j,k) - ( fz(i,j,k) - 0.0_wp ) / vol
!        ! this construction guarantees integral conservation of phi
!      END IF

!    END DO
!  END DO

  k=ke
  DO j=1, je
    dxdy = dx_aequ_dy * crlat(j,1)
    DO i=1, ie

      IF ( clip_flag(i,j,k-1) .OR. clip_flag(i,j,k) ) THEN
        dz = dzeta / sqrtg_r_s(i,j,k)
        vol = dxdy * dz
        phi(i,j,k) = phi(i,j,k) - ( 0.0_wp - fz(i,j,k-1) ) / vol
        ! this construction guarantees integral conservation of phi
      END IF

    END DO
  END DO

DEALLOCATE( fz )

END SUBROUTINE diffusion_z


!==============================================================================
!==============================================================================


SUBROUTINE bandec( a, n, m1, m2, al, indx, d )

!------------------------------------------------------------------------------
!
! Description:
!
! performs an LU decomposition of a banddiagonal matrix a.
! Fortran90-Routine analogous to 'Numerical Recipes in C, 2. ed.'
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER (KIND=iintegers), INTENT(IN)  :: n, m1, m2
  REAL (KIND=wp),         INTENT(INOUT) :: a(1:n, 1:m1+m2+1)
  REAL (KIND=wp),         INTENT(OUT)   :: al(1:n, 1:m1), d
  INTEGER (KIND=iintegers), INTENT(OUT) :: indx(1:n) 

  REAL (KIND=wp),     PARAMETER :: TINY=1.E-20_wp
  INTEGER (KIND=iintegers)   :: k2,k3,k,nrbd,l,mm
  REAL (KIND=wp)             :: dum

  mm=m1+m2+1
  l=m1

  DO k2=1,m1
    DO  nrbd=m1+2-k2,mm
      a(k2,nrbd-l)=a(k2,nrbd)
    END DO
    l=l-1
    DO  nrbd=mm-l,mm
      a(k2,nrbd)=0._wp
    END DO
  END DO

  d=1.0_wp
  l=m1

  DO k=1,n

    dum=a(k,1)
    k2=k
    IF ( l<n ) l=l+1
    DO k3=k+1,l
      IF ( ABS( a(k3,1) ) > ABS(dum) ) THEN
        dum=a(k3,1)
        k2=k3
      ENDIF
    END DO
    indx(k)=k2
    IF ( dum == 0.0_wp ) a(k,1)=TINY
    IF ( k2 .NE. k )THEN
      d = -d
      DO nrbd=1,mm
        dum=a(k,nrbd)
        a(k,nrbd)=a(k2,nrbd)
        a(k2,nrbd)=dum
      END DO
    ENDIF
    DO k2=k+1,l
      dum=a(k2,1)/a(k,1)
      al(k,k2-k)=dum
      DO nrbd=2,mm
        a(k2,nrbd-1)=a(k2,nrbd)-dum*a(k,nrbd)
      END DO
      a(k2,mm)=0._wp
    END DO
  END DO
  RETURN

END SUBROUTINE bandec


SUBROUTINE banbks( a, n, m1, m2, al, indx, b)

!------------------------------------------------------------------------------
!
! Description:
!
! linear equation solver for general banddiagonal systems
! Fortran90-Routine analogous to 'Numerical Recipes in C, 2. ed.'
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER (KIND=iintegers), INTENT(IN)    :: n, m1, m2, indx(n)
  REAL (KIND=wp),           INTENT(IN)    :: a(n,m1+m2+1), al(n,m1)
  REAL (KIND=wp),           INTENT(INOUT) :: b(n)

  INTEGER (KIND=iintegers) :: k,k2,l,mm
  REAL (KIND=wp)       :: dum

  mm=m1+m2+1
  l=m1
  DO  k=1,n
    k2=indx(k)
    IF ( k2 .NE. k ) THEN
      dum=b(k)
      b(k)=b(k2)
      b(k2)=dum
    ENDIF
    IF ( l < n ) l=l+1
    DO  k2=k+1,l
      b(k2)=b(k2)-al(k,k2-k)*b(k)
    END DO
  END DO
  l=1
  DO  k2=n,1,-1
    dum=b(k2)
    DO  k=2,l
      dum=dum-a(k2,k)*b(k+k2-1)
    END DO
    b(k2)=dum/a(k2,1)
    IF ( l < mm ) l=l+1
  END DO

  RETURN

END SUBROUTINE banbks



SUBROUTINE solve_5banddiag( lgs, lgs_rhs, sol, &
  ie, je, istart, iend, jstart, jend, kend ) 

  !--------------------------------------------------------------------------
  !
  ! Description:
  !   This subroutine solves ie * je vertical implicit schemes which lead to 
  !   linear equation systems with 5-band-diagonal-matrices (e.g. implicit 
  !   vertical advection 3. order)
  !
  ! Input fields:
  !   The 5 diagonals  lgs(:,:,:,1:5)  and the right hand side  lgs_rhs
  !   are read as  ie*je*kend  (kend=ke or ke+1) fields. 
  !
  ! Ouput fields:
  !   the ie*je solution vectors  sol(ie, je, kend)
  !
  ! Method:
  !   This is a driver routine for the 'Numerical Recipes'-Routines
  !   'bandec' and 'banbks' for general banddiagonal-systems
  !
  !--------------------------------------------------------------------------

  ! Declarations:

  IMPLICIT NONE

  INTEGER (KIND=iintegers), INTENT(IN) :: ie, je
  INTEGER (KIND=iintegers), INTENT(IN) :: istart, iend, jstart, jend
  INTEGER (KIND=iintegers), INTENT(IN) :: kend

  REAL (KIND=wp),     INTENT(IN) ::  lgs( 1:ie, 1:je, 1:kend, 1:5 )
  REAL (KIND=wp),     INTENT(IN) ::  lgs_rhs( 1:ie, 1:je, 1:kend )

  REAL (KIND=wp),       INTENT(OUT) :: sol( 1:ie, 1:je, 1:kend)

  ! Local variables
  !----------------

  INTEGER (KIND=iintegers)    ::  &
    i, j, k              ! loop indices

  REAL (KIND=wp),           ALLOCATABLE :: q(:,:)
  REAL (KIND=wp),           ALLOCATABLE :: ql(:,:)
  REAL (KIND=wp),           ALLOCATABLE :: x(:)

  INTEGER (KIND=iintegers), ALLOCATABLE :: indx(:)
  REAL (KIND=wp)         :: vorzchn

  ! -------- Allokieren ---------------------------

  ALLOCATE( q(1:kend, 1:5) )
  ALLOCATE( ql(1:kend, 1:2) )
  ALLOCATE( x(1:kend) )

  ALLOCATE( indx(1:kend) )

  DO j=jstart, jend
    DO i=istart, iend

      ! --- rearrange the 5-banddiagonals in field q: ---
      DO k=1, kend
        q(k,1) = lgs(i,j,k,1)
        q(k,2) = lgs(i,j,k,2)
        q(k,3) = lgs(i,j,k,3)
        q(k,4) = lgs(i,j,k,4)
        q(k,5) = lgs(i,j,k,5)
        x(k)   = lgs_rhs(i,j,k) ! necessary: 'banbks' overwrites the rhs
      END DO

      ! --- call of 'Numerical Recipes'-routines ---

      CALL bandec( q, kend, 2_iintegers, 2_iintegers, ql, indx, vorzchn);
      CALL banbks( q, kend, 2_iintegers, 2_iintegers, ql, indx, x );

      DO k=1, kend
        sol(i,j,k) = x(k)
      END DO

    END DO
  END DO

  DEALLOCATE( q, ql, indx, x)

END SUBROUTINE solve_5banddiag


!==============================================================================
!==============================================================================


SUBROUTINE bandec_vec( a, n, m1, m2, al, indx, d )

!------------------------------------------------------------------------------
!
! Description:
!
! performs an LU decomposition of a banddiagonal matrix a.
! Fortran90-Routine analogous to 'Numerical Recipes in C, 2. ed.'
! optimized for vector computers
!
!------------------------------------------------------------------------------

  USE data_modelconfig, ONLY:  &
     ie, je,                   &
     istart, iend,             &
     jstart, jend

  IMPLICIT NONE

  INTEGER (KIND=iintegers), INTENT(IN)    :: n, m1, m2
  REAL (KIND=wp),           INTENT(INOUT) :: a(1:ie, 1:je, 1:n, 1:m1+m2+1)
  REAL (KIND=wp),           INTENT(OUT)   :: al(1:ie, 1:je, 1:n, 1:m1)
  REAL (KIND=wp),           INTENT(OUT)   :: d(1:ie, 1:je)
  INTEGER (KIND=iintegers), INTENT(OUT)   :: indx(1:ie, 1:je, 1:n) 

  REAL (KIND=wp),     PARAMETER :: TINY=1.E-20_wp
  INTEGER (KIND=iintegers)   :: k,k2,k3,l,mm
  INTEGER (KIND=iintegers)   :: nrbd
  INTEGER (KIND=iintegers)   :: i, j
  INTEGER (KIND=iintegers), ALLOCATABLE ::  k2f(:,:)
  REAL (KIND=wp),           ALLOCATABLE ::  dum(:,:)

  ALLOCATE( k2f( 1:ie, 1:je ) )
  ALLOCATE( dum( 1:ie, 1:je ) )

  mm=m1+m2+1
  l=m1

  DO k2=1,m1
    DO  nrbd=m1+2-k2,mm
      DO j=jstart, jend
        DO i=istart, iend
          a(i,j,k2,nrbd-l)=a(i,j,k2,nrbd)
        END DO
      END DO
    END DO
    l=l-1
    DO  nrbd=mm-l,mm
      DO j=jstart, jend
        DO i=istart, iend
          a(i,j,k2,nrbd)=0._wp
        END DO
      END DO
    END DO
  END DO

  d(:,:)=1.0_wp
  l=m1

  DO k=1,n

    DO j=jstart, jend
      DO i=istart, iend
        dum(i,j)=a(i,j,k,1)
        k2f(i,j)=k
      END DO
    END DO
    IF ( l<n ) l=l+1
    DO k3=k+1,l
      DO j=jstart, jend
        DO i=istart, iend
          IF ( ABS( a(i,j,k3,1) ) > ABS( dum(i,j) ) ) THEN
            dum(i,j)=a(i,j,k3,1)
            k2f(i,j)=k3
          ENDIF
        END DO
      END DO
    END DO
    DO j=jstart, jend
      DO i=istart, iend
        indx(i,j,k)=k2f(i,j)
        IF ( dum(i,j) == 0.0_wp ) a(i,j,k,1)=TINY
      END DO
    END DO

    DO j=jstart, jend
      DO i=istart, iend

        IF ( k2f(i,j) .NE. k )THEN
          d(i,j) = -d(i,j)
          DO nrbd=1,mm
            dum(i,j)=a(i,j,k,nrbd)
            a(i,j,k,nrbd)=a(i,j,k2f(i,j),nrbd)
            a(i,j,k2f(i,j),nrbd)=dum(i,j)
          END DO
        ENDIF
      END DO
    END DO

    DO k2=k+1,l
      DO j=jstart, jend
        DO i=istart, iend

          dum(i,j)=a(i,j,k2,1)/a(i,j,k,1)
          al(i,j,k,k2-k)=dum(i,j)
          DO nrbd=2,mm
            a(i,j,k2,nrbd-1)=a(i,j,k2,nrbd)-dum(i,j)*a(i,j,k,nrbd)
          END DO
          a(i,j,k2,mm)=0._wp
        END DO
      END DO
    END DO

  END DO

  DEALLOCATE( k2f , dum )

  RETURN

END SUBROUTINE bandec_vec



SUBROUTINE banbks_vec( a, n, m1, m2, al, indx, b)

!------------------------------------------------------------------------------
!
! Description:
!
! linear equation solver for general banddiagonal systems
! Fortran90-Routine analogous to 'Numerical Recipes in C, 2. ed.'
! optimized for vector computers
!
!------------------------------------------------------------------------------

  USE data_modelconfig, ONLY:    &
    ie, je,      &
    istart, iend,  &
    jstart, jend

  IMPLICIT NONE

  INTEGER (KIND=iintegers), INTENT(IN)    :: n, m1, m2
  INTEGER (KIND=iintegers), INTENT(IN)    :: indx(1:ie, 1:je, n)
  REAL (KIND=wp),           INTENT(IN)    :: a(1:ie, 1:je, n, m1+m2+1)
  REAL (KIND=wp),           INTENT(IN)    :: al(1:ie, 1:je, n, m1)
  REAL (KIND=wp),           INTENT(INOUT) :: b(1:ie, 1:je, n)

  INTEGER (KIND=iintegers) :: k2,k,l,mm
  INTEGER (KIND=iintegers) :: i,j 
  REAL (KIND=wp),     ALLOCATABLE   :: dum(:,:)

  ALLOCATE( dum(1:ie, 1:je) )

  mm=m1+m2+1
  l=m1

  DO  k=1,n

    DO j=jstart, jend
      DO i=istart, iend

        k2=indx(i,j,k)
        IF ( k2 .NE. k ) THEN
          dum(i,j)=b(i,j,k)
          b(i,j,k)=b(i,j,k2)
          b(i,j,k2)=dum(i,j)
        ENDIF
      END DO
    END DO

    IF ( l < n ) l=l+1

    DO  k2=k+1,l
      DO j=jstart, jend
        DO i=istart, iend
          b(i,j,k2)=b(i,j,k2)-al(i,j,k,k2-k)*b(i,j,k)
        END DO
      END DO
    END DO

  END DO

  l=1

  DO  k2=n,1,-1
    DO j=jstart, jend
      DO i=istart, iend
        dum(i,j)=b(i,j,k2)
      END DO
    END DO
    DO  k=2,l
      DO j=jstart, jend
        DO i=istart, iend
          dum(i,j)=dum(i,j)-a(i,j,k2,k)*b(i,j,k+k2-1)
        END DO
      END DO
    END DO
    DO j=jstart, jend
      DO i=istart, iend
        b(i,j,k2)=dum(i,j)/a(i,j,k2,1)
      END DO
    END DO

    IF ( l < mm ) l=l+1
  END DO

  DEALLOCATE( dum )

  RETURN

END SUBROUTINE banbks_vec



SUBROUTINE solve_5banddiag_vec( lgs, lgs_rhs, sol,  &
  ie, je, istart, iend, jstart, jend, kend ) 

  !--------------------------------------------------------------------------
  !
  ! Description:
  !   This subroutine solves ie * je vertical implicit schemes which lead to 
  !   linear equation systems with 5-band-diagonal-matrices (e.g. implicit 
  !   vertical advection 3. order)
  !
  ! Input fields:
  !   The 5 diagonals  lgs(:,:,:,1:5)  and the right hand side  lgs_rhs(:,:,:)
  !   are read as  ie*je*kend fields (kend=ke or ke+1).
  !   ATTENTION: lgs(:,:,:,:) is modified by this routine !
  !
  ! Ouput fields:
  !   the ie*je solution vectors  sol(ie, je, kend)
  !
  ! Method:
  !   This is a driver routine for the 'Numerical Recipes'-Routines
  !   'bandec' and 'banbks' for general banddiagonal-systems
  !   Optimized for Vector computers (should be fast on scalar machines, too)
  !--------------------------------------------------------------------------

  ! Declarations:

  IMPLICIT NONE

  INTEGER (KIND=iintegers), INTENT(IN) :: ie, je
  INTEGER (KIND=iintegers), INTENT(IN) :: istart, iend, jstart, jend
  INTEGER (KIND=iintegers), INTENT(IN) :: kend

  REAL (KIND=wp),     INTENT(INOUT)  :: lgs( 1:ie, 1:je, 1:kend, 1:5 )
  REAL (KIND=wp),     INTENT(IN)     :: lgs_rhs( 1:ie, 1:je, 1:kend )

  REAL (KIND=wp),     INTENT(OUT) :: sol( 1:ie, 1:je, 1:kend)

  ! Local variables
  !----------------

  INTEGER (KIND=iintegers)    ::  &
    i, j, k              ! loop indices

  REAL (KIND=wp),           ALLOCATABLE :: ql(:,:,:,:)

  INTEGER (KIND=iintegers), ALLOCATABLE :: indx(:,:,:)
  REAL (KIND=wp)          , ALLOCATABLE :: vorzchn(:,:)

  ! -------- Allokieren ---------------------------

  ALLOCATE( ql(1:ie, 1:je, 1:kend, 1:2) )
  ALLOCATE( vorzchn(1:ie, 1:je) )
  ALLOCATE( indx(1:ie, 1:je, 1:kend) )

  DO k=1, kend
    DO j=jstart, jend
      DO i=istart, iend

        ! necessary, because 'banbks' overwrites the 'right hand side':
        sol(i,j,k) = lgs_rhs(i,j,k) 

      END DO
    END DO
  END DO

  ! --- call of 'Numerical Recipes'-routines ---

  CALL bandec_vec( lgs, kend, 2_iintegers, 2_iintegers, ql, indx, vorzchn);
  CALL banbks_vec( lgs, kend, 2_iintegers, 2_iintegers, ql, indx, sol );

  DEALLOCATE( ql, indx, vorzchn)

END SUBROUTINE solve_5banddiag_vec

!==============================================================================
!==============================================================================

SUBROUTINE mean_over_box (f, f_mean_2d, k_center, d_idx_x, d_idx_y, d_idx_z, &
                          crlat, sqrtg_r_s,                                  &
                          ie, je, ke, istart, iend, jstart, jend)

!------------------------------------------------------------------------------
!
! Description:
!   calculate the mean values of f over gliding boxes (in the terrain
!   following system)
!       [ i_center-d_idx_x, i_center+d_idx_x ]
!     * [ j_center-d_idx_y, j_center+d_idx_y ]
!     * [ k_center-d_idx_z, k_center+d_idx_z]
!   for all center points (i_center, j_center, k_center).
!
!   please obey:
!   1.) the field  sqrtg_r_s  is needed (ok for RK-dynamics, not for leapfrog)
!   2.) this routine does not carry out a security check,
!       if   d_idx_x, ... <= boundlines
!   3.) the field f must be properly exchanged before entering this routine
!
! Method:
!
!------------------------------------------------------------------------------

! Declarations:

INTEGER (KIND=iintegers), INTENT(IN)     ::     &
  k_center,     & ! vertical center for computing mean values
  d_idx_x,      & ! extent in x-, y- and z-direction
  d_idx_y,      & !
  d_idx_z         !

INTEGER (KIND=iintegers), INTENT(IN)     ::     &
  ie,           & ! number of grid points in zonal direction
  je,           & ! number of grid points in meridional direction
  ke,           & ! number of grid points in vertical direction
  istart,       & ! start index for the forecast of w, t, qv, qc and pp
  iend,         & ! end index for the forecast of w, t, qv, qc and pp
  jstart,       & ! start index for the forecast of w, t, qv, qc and pp
  jend            ! end index for the forecast of w, t, qv, qc and pp

REAL    (KIND=wp),        INTENT(IN)     ::     &
  f(ie,je,ke)     ! field, for which means are computed

REAL    (KIND=wp),        INTENT(IN)     ::     &
  sqrtg_r_s(ie,je,ke),    & ! cosine of transformed latitude
  crlat    (je,2)

REAL    (KIND=wp),        INTENT(OUT)    ::     &
  f_mean_2d(ie,je)          ! Field with the results

! Local Variables
!REAL    (KIND=wp)        :: z_sum1_f,      z_sum_f
!REAL    (KIND=wp)        :: z_sum1_vol,  z_sum_vol
INTEGER (KIND=iintegers) :: i, j
INTEGER (KIND=iintegers) :: ibox, jbox, kbox  ! center of box
REAL    (KIND=wp)        :: z_sum2_f(ie,je), z_sum2_vol(ie,je)

!--------------------------------------------------------------------------

 z_sum2_f(:,:)   = 0.0_wp
 z_sum2_vol(:,:) = 0.0_wp
 f_mean_2d(:,:)  = 0.0_wp      ! eigentlich unnoetig

  DO jbox = -d_idx_y, d_idx_y
    DO kbox = -d_idx_z, d_idx_z
!CDIR OUTERUNROLL=3
      DO ibox = -d_idx_x, d_idx_x

        DO j = jstart, jend
!CDIR ON_ADB(z_sum2_f)
!CDIR ON_ADB(z_sum2_vol)
          DO i = istart, iend
            z_sum2_f(i,j)   = z_sum2_f(i,j)                             &
                + (f(i+ibox,j+jbox,k_center+kbox)                       &
                        / sqrtg_r_s(i+ibox,j+jbox,k_center+kbox))       &
                   * crlat(j+jbox, 1)
            z_sum2_vol(i,j) = z_sum2_vol(i,j) +                         &
             (1.0_wp / sqrtg_r_s(i+ibox,j+jbox,k_center+kbox))      &
                   * crlat(j+jbox, 1)
          ENDDO
        ENDDO

      ENDDO
    ENDDO
  ENDDO

  DO j = jstart, jend
!CDIR ON_ADB(z_sum2_f)
!CDIR ON_ADB(z_sum2_vol)
    DO i = istart, iend
      f_mean_2d(i,j) = z_sum2_f(i,j) / z_sum2_vol(i,j)
    ENDDO
  ENDDO

END SUBROUTINE mean_over_box

!==============================================================================
!==============================================================================

SUBROUTINE mean_cov_over_box (f1, f1_mean_2d, f2, f2_mean_2d, corr_2d, &
                              k_center, d_idx_x, d_idx_y, d_idx_z,     &
                              crlat, sqrtg_r_s,                        &
                              ie, je, ke, istart, iend, jstart, jend)

!------------------------------------------------------------------------------
!
! Description:
!
!   calculate the correlation of < f1' * f2' > over gliding boxes
!   (in terrain following system)
!       [ i_center-d_idx_x, i_center+d_idx_x ]
!     * [ j_center-d_idx_y, j_center+d_idx_y ]
!     * [ k_center-d_idx_z, k_center+d_idx_z]
!   for all center points (i_center, j_center, k_center)
!
!   please obey:
!   1.) the field  sqrtg_r_s  is needed (ok for RK-dynamics, not for leapfrog)
!   2.) this routine does not carry out a security check,
!       if   d_idx_x, ... <= boundlines
!   3.) the field f must be properly exchanged before entering this routine
!
! Method:
!
!------------------------------------------------------------------------------

! Declarations:

INTEGER (KIND=iintegers), INTENT(IN)     ::     &
  k_center,     & ! vertical center for computing mean values
  d_idx_x,      & ! extent in x-, y- and z-direction
  d_idx_y,      & !
  d_idx_z         !

INTEGER (KIND=iintegers), INTENT(IN)     ::     &
  ie,           & ! number of grid points in zonal direction
  je,           & ! number of grid points in meridional direction
  ke,           & ! number of grid points in vertical direction
  istart,       & ! start index for the forecast of w, t, qv, qc and pp
  iend,         & ! end index for the forecast of w, t, qv, qc and pp
  jstart,       & ! start index for the forecast of w, t, qv, qc and pp
  jend            ! end index for the forecast of w, t, qv, qc and pp

REAL    (KIND=wp),        INTENT(IN)     ::     &
  f1(ie,je,ke), & ! field, for which means are computed
  f2(ie,je,ke)    !

REAL    (KIND=wp),        INTENT(IN)     ::     &
  f1_mean_2d(ie,je), & !
  f2_mean_2d(ie,je)    !

REAL    (KIND=wp),        INTENT(IN)     ::     &
  sqrtg_r_s(ie,je,ke),    & ! cosine of transformed latitude
  crlat    (je,2)

REAL    (KIND=wp),        INTENT(OUT)    ::     &
  corr_2d(ie,je)            ! Field with the results

! Local Variables
!REAL    (KIND=wp)        :: z_sum1_ff,    z_sum_ff
!REAL    (KIND=wp)        :: z_sum1_vol,  z_sum_vol
REAL    (KIND=wp)        :: z_sum2_ff(ie,je), z_sum2_vol(ie,je)
REAL    (KIND=wp)        :: z_corr
INTEGER (KIND=iintegers) :: i, j, ibox, jbox, kbox

!------------------------------------------------------------------------------

  z_sum2_ff (:,:) = 0.0_wp
  z_sum2_vol(:,:) = 0.0_wp
  corr_2d   (:,:) = 0.0_wp

  DO jbox = -d_idx_y, d_idx_y
    DO kbox = -d_idx_z, d_idx_z
!CDIR OUTERUNROLL=3
      DO ibox = -d_idx_x, d_idx_x

        DO j = jstart, jend
!CDIR ON_ADB(z_sum2_ff)
!CDIR ON_ADB(z_sum2_vol)
          DO i = istart, iend
            z_corr = ( f1(i+ibox,j+jbox,k_center+kbox) - f1_mean_2d(i,j) )  &
                   * ( f2(i+ibox,j+jbox,k_center+kbox) - f2_mean_2d(i,j) )
            z_sum2_ff(i,j) = z_sum2_ff(i,j) + &
                 z_corr / sqrtg_r_s(i+ibox, j+jbox, k_center+kbox) * crlat(j+jbox,1)
            z_sum2_vol(i,j) = z_sum2_vol(i,j) + &
               1.0_wp / sqrtg_r_s(i+ibox,j+jbox,k_center+kbox) * crlat(j+jbox,1)
          ENDDO
        ENDDO

      ENDDO
    ENDDO
  ENDDO

  DO j = jstart, jend
!CDIR ON_ADB(z_sum2_ff)
!CDIR ON_ADB(z_sum2_vol)
    DO i = istart, iend
      corr_2d(i,j) = z_sum2_ff(i,j) / z_sum2_vol(i,j)
    ENDDO
  ENDDO

END SUBROUTINE mean_cov_over_box

!==============================================================================
!==============================================================================

SUBROUTINE vert_avg (vertavg, f, sqrtg_r_s, ie, je, ke, istart, iend,      &
                     jstart, jend, k_center, d_idx_z )

!--------------------------------------------------------------------------
!
! Description:
!   calculate the average of f over the vertical column
!     [ k_center-d_idx_z, k_center+d_idx_z]
!   at the horizontal position (i,j)
!
!--------------------------------------------------------------------------

! Declarations:

!USE data_modelconfig, ONLY :   &
!  ie,           & ! number of grid points in zonal direction
!  je,           & ! number of grid points in meridional direction
!  ke              ! number of grid points in vertical direction
!
!USE data_fields, ONLY: sqrtg_r_s

INTEGER (KIND=iintegers), INTENT(IN)  ::       &
  ie, je, ke,            & ! dimensions of fields
  istart, iend, jstart, jend, k_center, d_idx_z  !

REAL    (KIND=wp),        INTENT(IN)  ::       &
  f         (ie,je,ke),  & !
  sqrtg_r_s (ie,je,ke)

REAL    (KIND=wp),        INTENT(OUT) ::       &
  vertavg (ie,je)

INTEGER (KIND=iintegers) :: i,j,k
REAL    (KIND=wp)        :: z_sum, z_vol

!--------------------------------------------------------------------------

DO j=jstart,jend
  DO i=istart,iend
    z_sum = 0.0_wp
    z_vol = 0.0_wp
    DO k = k_center - d_idx_z, k_center + d_idx_z
      z_sum = z_sum + sqrtg_r_s(i,j,k) * f(i,j,k)
      z_vol = z_vol + sqrtg_r_s(i,j,k)
    END DO

    vertavg(i,j) = z_sum / z_vol
  END DO
END DO

END SUBROUTINE vert_avg

!==============================================================================

END MODULE numeric_utilities
