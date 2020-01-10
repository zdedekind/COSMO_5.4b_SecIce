!+ Module providing routines for numerical schemes used in Runge-Kutta scheme
!------------------------------------------------------------------------------

MODULE numeric_utilities_rk

!------------------------------------------------------------------------------
!
! Description:
!   This module provides some service utilities related to various numerical
!   schemes used in the Runge-Kutta scheme.
!     - no routine uses other modules, except the declarations for the
!       KIND-type parameter; the data access is by parameter list only
!     - no routine allocates dynamic memory; work space needed is
!       provided via the parameter list
!     - no derived data types are used
!   The routines are written in plug-compatible format. Purists may 
!   identify the input/output variables to start with the letter 'p'.
!
!   Routines contained for the Runge-Kutta schemes:
!
!     - udsdx                    udsdx_* interfaces (optimization for MeteoCH)
!     - udsdx_up5_xy
!     - udsdx_up3_z
!     - udsdx_up1_z
!
!     - udsdx_up1                these operators are in advection form:
!     - udsdx_cd2
!     - udsdx_up3                U * ds / dx
!     - udsdx_cd4
!     - udsdx_up5
!     - udsdx_cd6
!
!     - dfdx_up1                 these operators are in flux form
!     - dfdx_cd2
!     - dfdx_up3
!     - dfdx_cd4
!     - dfdx_up5
!     - dfdx_cd6
!     - flux_weno3
!     - flux_weno3_priv
!
!     - xadv_rk_cri_ppm
!     - yadv_rk_cri_ppm
!     - zadv_rk_cri_ppm
!     - xadv_rk_ppm
!     - yadv_rk_ppm
!     - zadv_rk_ppm
!     - xadv_pd_rk_cri_vanleer
!     - yadv_pd_rk_cri_vanleer
!     - zadv_pd_rk_cri_vanleer
!     - xadv_pd_rk_vanleer
!     - yadv_pd_rk_vanleer
!     - zadv_pd_rk_vanleer
!     - init_bott_coeffs
!     - xadv_pd_rk_cri_bott
!     - yadv_pd_rk_cri_bott
!     - zadv_pd_rk_cri_bott
!     - xadv_pd_rk_bott
!     - yadv_pd_rk_bott
!     - zadv_pd_rk_bott
!     - ufrac_crint_rk
!     - vfrac_crint_rk
!     - wcfrac_crint_rk
!
!     - clipping
!     - clipping2
!     - clipping_DDI
!     - integral_3d_wg
!     - integral_3d_wg_DDI
!     - multiplicative_filling
!     - multiplicative_filling_DDI
!     - sum_DDI
!
! Current Code Owner: DWD, Michael Baldauf
!  phone:  +49  69  8062 2733
!  fax:    +49  69  8062 3721
!  email:  Michael.Baldauf@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 3.16       2005/07/22 Jochen Foerstner
!  Initial release
! 3.18       2006/03/03 Jochen Foerstner
!  Update of the scheme; Optimizations for udsdx-interfaces
! 3.21       2006/12/04 Michael Baldauf, Jochen Foerstner
!  New subroutines clipping_2, integral_3d_wg, multiplicative_filling
! V4_1         2007/12/04 Ulrich Schaettler, Michael Baldauf
!  Global exchange of variable icr in SR wcfrac_crint_rk, to get reproducible results
!  New subroutines for reproducibility of results for SemiLagrange advection
! V4_4         2008/07/16 Michael Baldauf, Lucio Torrisi, Oliver Fuhrer
!  Consideration of the sign of dt for DFI (Lucio Torrisi)
!  Bugfix in the multiplicative filling (Oliver Fuhrer)
! V4_5         2008/09/10 NEC staff
!  Modify initial Bott coefficients, to avoid divisions during the forecast
! V4_8         2009/02/16 Oliver Fuhrer
!  Use global values only if num_compute greater 1
! V4_9         2009/07/16 Ulrich Schaettler, Christian Bollmann
!  Implemented index lists for bott-cri advection
!  Implemented 3D versions of several advection operators
!  Use gather_values (in Function sum_DDI) only if num_compute > 1 (Michael Baldauf)
! V4_10        2009/09/11 Christian Bollmann
!  Added compiler directive to use option _on_adb for NEC
! V4_12        2010/05/11 Stephan Pfahl (Ulrich Schaettler)
!  Bug correction in SR zadv_pd_rk_cri_bott for 4th order option
!  Bug correction in index lists for zadv_pd_rk_cri_bott 2nd order option
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_18        2011/05/26 Michael Baldauf
!  New NL switch y_scalar_advect replaces lsl_adv_qx, yef_adv_qx
! V4_20        2011/08/31 Ulrich Schaettler
!  Print debug output only in case of ldebug_dyn
! V4_23        2012/05/10 Michael Baldauf, Oliver Fuhrer
!  Allow new variants of the Bott-Advection schemes:
!    BOTT2_STRANG_B, BOTT4_STRANG_B: Strang Splitting only at the bottom (lowest 5 levels)
!    BOTT2_XYZYX, BOTT4_XYZYX: modified sequence compared to the current Strang splitting
!  Use field sqrtg_r_s from new module grid_metrics_utilities
!  Eliminated IN-argument of a field s to advection operators and used OUT-argument
!  of that field now as INOUT. (Oliver Fuhrer)
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer
!  Renamed variable ief_adv_qx to ief_adv_tracer and changed some comments
! V4_27        2013/03/19 Astrid Kerkweg
!  MESSy interface introduce: define eps_div=1.E-25 in Bott-routines
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
! V5_3         2015-10-09 Hans-Juergen Panitz, Ulrich Schaettler
!  Replace KIND=8 for INTEGER variables by generic definition i8 
!     from kind_parameters (HJP)
!  Removed all iinteger declarations (US)
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
USE kind_parameters , ONLY :   &
  wp        ,& ! KIND-type parameters for real variables
  i8           ! KIND-type parameter for integer*8 variables

USE data_constants,   ONLY :   &
  repsilon     ! precision of 1.0 in current floating point format

USE parallel_utilities, ONLY: global_values

!------------------------------------------------------------------------------ 
!
! Declarations:    in the module procedures
!
!------------------------------------------------------------------------------
 
IMPLICIT NONE

REAL (KIND=wp),  SAVE ::   abc(0:4,0:4,-2:2)
INTEGER,         SAVE ::   ief_adv_trcr = 0 ! to check only once,
                                            ! if Euler_forward schemes are chosen
! epsilons (security constants)
REAL (KIND=wp), PARAMETER :: &
  eps_div = repsilon , & ! avoid division by zero
#ifndef MESSY
  eps_adv = 1.0E-15_wp   ! general epsilon in advection
#else
  eps_adv = 1.0E-25_wp   ! 1.0E-15p is too large for small tracers as Radon
#endif

!------------------------------------------------------------------------------
!
! Module procedures of "numeric_utilities_rk"
!
!------------------------------------------------------------------------------

CONTAINS

!==============================================================================
!==============================================================================

SUBROUTINE xadv_rk_cri_ppm( u_frac, icr, rdx, s, dt, ie, je, &
                            istart, iend, jstart, jend, sqrtg_r )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine updates a scalar prognostic
!   variable 's' due to x-advection. The routine uses the
!   directional split method.
!
! Method:
!   Currently, the PPM advection scheme (Carpenter 1990) is used.
!   Courant number independent formulation  
!
! Declarations:
!
!------------------------------------------------------------------------------

! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER          ,   INTENT (IN) ::  &
  ie, je,         & ! dimensions of the input/output arrays
  istart, iend,   & ! computational start and end indices in i-direction
                    ! (west-east for lat/lon grid)
  jstart, jend      ! computational start and end indices in j-direction
                    ! (south-north for lat/lon grid)

REAL    (KIND=wp),   INTENT(IN)  ::  &
  dt                ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp),   INTENT(IN)  ::  &
  sqrtg_r(ie,je), & ! 1 / square root of G
  rdx(je),        & ! reciprocal x-grid spacing (1/( a cos(phi)dlam )
  u_frac(ie,je)     ! fractional u-velocity

INTEGER          ,   INTENT(IN)  ::  &
  icr(ie,je)        ! integer courant numbers

! Array arguments with intent(out):

REAL    (KIND=wp),   INTENT(INOUT) ::  &
  s(ie,je)          ! scalar field updated due to x-advection

!------------------------------------------------------------------------------
! Local variables

INTEGER        :: &
  i, j, l,        & ! loop indices
  ziif

REAL (KIND=wp) :: &
  zsgdx,          & !
  zcrn              ! Courant number for intermediate updating

! Local automatic arrays
REAL (KIND=wp) :: &
  sze(ie,je),     & ! scalar field to be transported at zone edges
  zfx(ie,je)        ! Flux in x-direction at cell boundary (u-point)

! End of header
!==============================================================================

  !----------------------------------------------------------------------------
  ! Section 1: Compute monotonic slopes and form fluxes in x - direction
  !----------------------------------------------------------------------------

  zfx (:,:) = 0.0_wp

  DO j = jstart, jend
    DO i = 2, ie-2
      sze(i,j) = ( 7.0_wp*( s(i+1,j)+s(i  ,j) )                  &
                            - ( s(i+2,j)+s(i-1,j) ) ) / 12.0_wp
    ENDDO
  ENDDO
  
  ! Compute courant number and velo-flux at zone edges (u-points)
  DO j = jstart, jend
    DO i = istart-1, iend
      ! Compute integer fluxes
      IF ( icr(i,j) >= 1 ) THEN
        zsgdx = 1.0_wp / rdx(j) / sqrtg_r(i,j)
        DO l = 1, icr(i,j)
          ziif = i-l
          zfx(i,j) = zfx(i,j) + s(ziif+1,j) * zsgdx
        ENDDO 
      ELSE IF ( icr(i,j) <= -1 ) THEN
        zsgdx = 1.0_wp / rdx(j) / sqrtg_r(i+1,j)
        DO l = -1, icr(i,j), -1
          ziif = i-l
          zfx(i,j) = zfx(i,j) - s(ziif  ,j) * zsgdx
        ENDDO
      ELSE
        ziif = i
      END IF
      ! Compute and add (fractional) fluxes
      zcrn = dt*u_frac(i,j)*rdx(j)
      IF ( zcrn > 0.0_wp ) THEN
        zcrn = zcrn * sqrtg_r(i,j)
        zfx(i,j) = zfx(i,j) + dt*u_frac(i,j) *  &
          ( sze(ziif,j) - zcrn*( sze(ziif  ,j)-s(ziif  ,j) )           &
                        - zcrn*(1.0_wp-zcrn)*                      &
                               ( sze(ziif-1,j)-2.0_wp*s(ziif  ,j)  &
                               + sze(ziif  ,j) ) )
      ELSE IF ( zcrn < 0.0_wp ) THEN
        zcrn = zcrn * sqrtg_r(i+1,j)
        zfx(i,j) = zfx(i,j) + dt*u_frac(i,j) *  &
          ( sze(ziif,j) + zcrn*( sze(ziif  ,j)-s(ziif+1,j) )           &
                        + zcrn*(1.0_wp+zcrn)*                      &
                               ( sze(ziif  ,j)-2.0_wp*s(ziif+1,j)  &
                               + sze(ziif+1,j) ) )
      END IF
    ENDDO
  ENDDO

  !----------------------------------------------------------------------------
  ! Section 2: Update the scalar field
  !----------------------------------------------------------------------------
  DO j = jstart, jend
    DO i = istart, iend
      s(i,j) = s(i,j) - ( zfx(i,j) - zfx(i-1,j) ) * rdx(j)*sqrtg_r(i,j)
    ENDDO
  ENDDO

END SUBROUTINE xadv_rk_cri_ppm

!==============================================================================
!============================================================================== 

SUBROUTINE yadv_rk_cri_ppm( v_frac, icr, rdy, s, dt, ie, je, &
                            istart, iend, jstart, jend, sqrtg_r )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine updates a scalar prognostic
!   variable 's' due to y-advection. The routine uses the
!   directional split method.
!
! Method:
!   Currently, the PPM advection scheme (Carpenter 1990) is used.
!   Courant number independent formulation  
!
! Declarations:
!
!------------------------------------------------------------------------------

! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER          ,   INTENT (IN) ::  &
  ie, je,         & ! dimensions of the input/output arrays
  istart, iend,   & ! computational start and end indices in i-direction
                    ! (west-east for lat/lon grid)
  jstart, jend      ! computational start and end indices in j-direction
                    ! (south-north for lat/lon grid)

REAL    (KIND=wp),   INTENT(IN)  ::  &
  dt                ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp),   INTENT(IN)  ::  &
  sqrtg_r(ie,je), & ! 1 / square root of G
  rdy(je,2),      & ! reciprocal y-grid spacing (1/( a cos(phi)dphi )
  v_frac(ie,je)     ! fractional v-velocity

INTEGER          ,   INTENT(IN)  ::  &
  icr(ie,je)        ! integer courant numbers

! Array arguments with intent(out):

REAL    (KIND=wp),   INTENT(INOUT) ::  &
  s(ie,je)          ! scalar field updated due to y-advection

!------------------------------------------------------------------------------
! Local variables

INTEGER        :: &
  i, j, l,        & ! loop indices
  zjif

REAL (KIND=wp) :: &
  zsgdy,          & !
  zcrn              ! Courant number for intermediate updating

! Local automatic arrays
REAL (KIND=wp) :: &
  sze(ie,je),     & ! scalar field to be transported at zone edges
  zfy(ie,je)        ! Flux in y-direction at cell boundary (v-point)

! End of header
!==============================================================================

  !----------------------------------------------------------------------------
  ! Section 1: Compute monotonic slopes and form fluxes in y - direction
  !----------------------------------------------------------------------------

  zfy (:,:) = 0.0_wp

  DO j = 2, je-2
    DO i = istart, iend
      sze(i,j) = ( 7.0_wp*( s(i,j+1)+s(i,j  ) )                  &
                            - ( s(i,j+2)+s(i,j-1) ) ) / 12.0_wp
    ENDDO
  ENDDO
  
  ! Compute courant number and velo-flux at zone edges (v-points)
  DO j = jstart-1, jend
    DO i = istart, iend
      ! Compute integer fluxes
      IF ( icr(i,j) >= 1 ) THEN
        zsgdy = 1.0_wp / rdy(j,2) / sqrtg_r(i,j)
        DO l = 1, icr(i,j) 
          zjif = j-l
          zfy(i,j) = zfy(i,j) + s(i,zjif+1) * zsgdy
        ENDDO 
      ELSE IF ( icr(i,j) <= -1 ) THEN
        zsgdy = 1.0_wp / rdy(j,2) / sqrtg_r(i,j+1)
        DO l = -1, icr(i,j), -1
          zjif = j-l
          zfy(i,j) = zfy(i,j) - s(i,zjif  ) * zsgdy
        ENDDO
      ELSE
        zjif = j
      END IF
      ! Compute and add (fractional) fluxes
      zcrn = dt*v_frac(i,j)*rdy(j,2)
      IF ( zcrn > 0.0_wp ) THEN
        zcrn = zcrn * sqrtg_r(i,j)
        zfy(i,j) = zfy(i,j) + dt*v_frac(i,j) *  &
          ( sze(i,zjif) - zcrn*( sze(i,zjif  )-s(i,zjif  ) )           &
                        - zcrn*(1.0_wp-zcrn)*                      &
                               ( sze(i,zjif-1)-2.0_wp*s(i,zjif  )  &
                               + sze(i,zjif  ) ) )
      ELSE IF ( zcrn < 0.0_wp ) THEN
        zcrn = zcrn * sqrtg_r(i,j+1)
        zfy(i,j) = zfy(i,j) + dt*v_frac(i,j) *  &
          ( sze(i,zjif) + zcrn*( sze(i,zjif  )-s(i,zjif+1) )           &
                        + zcrn*(1.0_wp+zcrn)*                      &
                               ( sze(i,zjif  )-2.0_wp*s(i,zjif+1)  &
                               + sze(i,zjif+1) ) )
      END IF
    ENDDO
  ENDDO

  !----------------------------------------------------------------------------
  ! Section 2: Update the scalar field
  !----------------------------------------------------------------------------
  DO j = jstart, jend
    DO i = istart, iend
      s(i,j) = s(i,j) - ( zfy(i,j) - zfy(i,j-1) ) * rdy(j,1)*sqrtg_r(i,j)
    ENDDO
  ENDDO

END SUBROUTINE yadv_rk_cri_ppm

!==============================================================================
!============================================================================== 

SUBROUTINE zadv_rk_cri_ppm( wc_frac, icr, s, dt, ie, je, ke, ke1, &
                            istart, iend, jstart, jend, sqrtg_r )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine updates a scalar prognostic
!   variable 's' due to vertical advection. The routine uses the
!   directional split method.
!
! Method:
!   Currently, the PPM advection scheme (Carpenter 1990) is used.
!   Courant number independent formulation  
!
! Declarations:
!
!------------------------------------------------------------------------------

! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER          ,   INTENT (IN) ::  &
  ie, je,            & ! x- and y-dimension of the input/output arrays
  ke, ke1,           & ! z-dimension of the input/output arrays
  istart, iend,      & ! computational start and end indices in i-direction
                       ! (west-east for lat/lon grid)
  jstart, jend         ! computational start and end indices in j-direction
                       ! (south-north for lat/lon grid)

REAL    (KIND=wp),   INTENT(IN)  ::  &
  dt                   ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp),   INTENT(IN)  ::  &
  sqrtg_r(ie,je,ke), & ! 1 / square root of G
  wc_frac(ie,je,ke1)   ! fractional contravariant vertical velocity
                       ! (at half levels)

INTEGER          ,   INTENT(IN)  ::  &
  icr(ie,je,ke1)       ! integer courant numbers (at half levels)

! Array arguments with intent(inout):

REAL    (KIND=wp),   INTENT(INOUT) ::  &
  s(ie,je,ke)          ! scalar field updated due to z-advection

!------------------------------------------------------------------------------
! Local variables

INTEGER        :: &
  i, j, k, l,     & ! loop indices
  zkif

REAL (KIND=wp) :: &
  zsgdz,          & !
  zcrn, zcrn_sgr    ! Courant number for intermediate updating

! Local automatic arrays
REAL (KIND=wp) :: &
  sze(ie,je,ke1), & ! scalar field to be transported at zone edges
  zfz(ie,je,ke1)    ! Flux in zeta-direction at cell boundary (w-point)

! End of header
!==============================================================================

  !----------------------------------------------------------------------------
  ! Section 1: Compute monotonic slopes and form fluxes in zeta - direction
  !----------------------------------------------------------------------------

  zfz (:,:,:) = 0.0_wp

  DO j = jstart, jend
    DO i = istart, iend
      sze(i,j,2)  = 0.5_wp*( s(i,j,2)+s(i,j,1) )
      sze(i,j,ke) = 0.5_wp*( s(i,j,ke)+s(i,j,ke-1) )
    ENDDO
  ENDDO
  DO k = 3, ke-1
    DO j = jstart, jend
      DO i = istart, iend
        sze(i,j,k) = ( 7.0_wp*( s(i,j,k  )+s(i,j,k-1) )                  &
                                - ( s(i,j,k+1)+s(i,j,k-2) ) ) / 12.0_wp
      ENDDO
    ENDDO
  ENDDO
  
  ! Compute courant number and velo-flux at zone edges (w-points)
  DO k = 2, ke
    DO j = jstart, jend
      DO i = istart, iend
        ! Compute integer fluxes
        IF ( icr(i,j,k) >= 1 ) THEN
          zsgdz = 1.0_wp / sqrtg_r(i,j,k-1)
          DO l = 1, icr(i,j,k) 
            zkif = k-l
            zfz(i,j,k) = zfz(i,j,k) + s(i,j,zkif  ) * zsgdz
          ENDDO 
        ELSE IF ( icr(i,j,k) <= -1 ) THEN
          zsgdz = 1.0_wp / sqrtg_r(i,j,k)
          DO l = -1, icr(i,j,k), -1
            zkif = k-l
            zfz(i,j,k) = zfz(i,j,k) - s(i,j,zkif-1) * zsgdz
          ENDDO
        ELSE
          zkif = k
        END IF
        ! Compute and add (fractional) fluxes
        zcrn = dt*wc_frac(i,j,k)
        IF ( zcrn > 0.0_wp ) THEN
          zcrn_sgr = zcrn * sqrtg_r(i,j,k-1)
          IF ( zkif-1 >= 2 ) THEN
            zfz(i,j,k) = zfz(i,j,k) + zcrn *  &
              ( sze(i,j,zkif) - zcrn_sgr*( sze(i,j,zkif  )-s(i,j,zkif-1) )   &
                              - zcrn_sgr*(1.0_wp-zcrn_sgr)*              &
                                 ( sze(i,j,zkif-1)-2.0_wp*s(i,j,zkif-1)  &
                                 + sze(i,j,zkif  ) ) )
          ELSE
            ! use upwind 1st order at upper boundary
            zfz(i,j,k) = zfz(i,j,k) + zcrn * s(i,j,zkif-1)
          END IF
        ELSE IF ( zcrn < 0.0_wp ) THEN
          zcrn_sgr = zcrn * sqrtg_r(i,j,k)
          IF ( zkif+1 <= ke ) THEN
            zfz(i,j,k) = zfz(i,j,k) + zcrn *  &
              ( sze(i,j,zkif) + zcrn_sgr*( sze(i,j,zkif  )-s(i,j,zkif  ) )   &
                              + zcrn_sgr*(1.0_wp+zcrn_sgr)*              &
                                 ( sze(i,j,zkif  )-2.0_wp*s(i,j,zkif  )  &
                                 + sze(i,j,zkif+1) ) )
          ELSE
            ! use upwind 1st order at lower boundary
            zfz(i,j,k) = zfz(i,j,k) + zcrn * s(i,j,zkif  )
          END IF
        END IF
      ENDDO
    ENDDO
  ENDDO

  !----------------------------------------------------------------------------
  ! Section 2: Update the scalar field
  !----------------------------------------------------------------------------
  DO k = 1, ke
    DO j = jstart, jend
      DO i = istart, iend
        s(i,j,k) = s(i,j,k) - ( zfz(i,j,k+1) - zfz(i,j,k) ) * sqrtg_r(i,j,k)
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE zadv_rk_cri_ppm

!==============================================================================
!============================================================================== 

SUBROUTINE xadv_rk_ppm( u, rdx, s, dt, ie, je, ke,           &
                        istart, iend, jstart, jend, kstart, kend, sqrtg_r )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine updates a scalar prognostic
!   variable 's' due to x-advection. The routine uses the
!   directional split method.
!
! Method:
!   Currently, the PPM advection scheme (Carpenter 1990) is used.
!
! Declarations:
!
!------------------------------------------------------------------------------

! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER          ,   INTENT (IN) ::  &
  ie, je, ke,     & ! dimensions of the input/output arrays
  istart, iend,   & ! computational start and end indices in i-direction
                    ! (west-east for lat/lon grid)
  jstart, jend,   & ! computational start and end indices in j-direction
                    ! (south-north for lat/lon grid)
  kstart, kend

REAL    (KIND=wp),   INTENT(IN)  ::  &
  dt                ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp),   INTENT(IN)  ::  &
  sqrtg_r(ie,je,kstart:kend), & ! 1 / square root of G
  rdx       (je),             & ! reciprocal x-grid spacing (1/( a cos(phi)dlam )
  u      (ie,je,kstart:kend)    ! u-velocity

! Array arguments with intent(out):

REAL    (KIND=wp),   INTENT(INOUT) ::  &
  s     (ie,je,kstart:kend)     ! scalar field updated due to x-advection

!------------------------------------------------------------------------------
! Local variables

INTEGER        :: &
  i, j, k, l        ! loop indices

REAL (KIND=wp) :: &
  zsgdx,          & !
  zcrn              ! Courant number for intermediate updating

! Local automatic arrays
REAL (KIND=wp) :: &
  sze(ie,je),     & ! scalar field to be transported at zone edges
  zfx(ie,je)        ! Flux in x-direction at cell boundary (u-point)

! End of header
!==============================================================================

  !----------------------------------------------------------------------------
  ! Section 1: Compute monotonic slopes and form fluxes in x - direction
  !----------------------------------------------------------------------------

  DO k = kstart, kend

    zfx (:,:) = 0.0_wp

    DO j = jstart, jend
      DO i = 2, ie-2
        sze(i,j) = ( 7.0_wp*( s(i+1,j,k)+s(i  ,j,k) )                  &
                              - ( s(i+2,j,k)+s(i-1,j,k) ) ) / 12.0_wp
      ENDDO
    ENDDO
  
    ! Compute courant number and velo-flux at zone edges (u-points)
    DO j = jstart, jend
      DO i = istart-1, iend
        ! Compute fluxes
        zcrn = dt*u(i,j,k)*rdx(j)
        IF ( zcrn > 0.0_wp ) THEN
          zcrn = zcrn * sqrtg_r(i,j,k)
          zfx(i,j) = dt*u(i,j,k) *  &
            ( sze(i,j) - zcrn*( sze(i  ,j)-s(i  ,j,k) )           &
                        - zcrn*(1.0_wp-zcrn)*                 &
                               ( sze(i-1,j)-2.0_wp*s(i  ,j,k) &
                               + sze(i  ,j) ) )
        ELSE IF ( zcrn < 0.0_wp ) THEN
          zcrn = zcrn * sqrtg_r(i+1,j,k)
          zfx(i,j) = dt*u(i,j,k) *  &
            ( sze(i,j) + zcrn*( sze(i  ,j)-s(i+1,j,k) )           &
                        + zcrn*(1.0_wp+zcrn)*                 &
                               ( sze(i  ,j)-2.0_wp*s(i+1,j,k) &
                               + sze(i+1,j) ) )
        END IF
      ENDDO
    ENDDO

  !----------------------------------------------------------------------------
  ! Section 2: Update the scalar field
  !----------------------------------------------------------------------------
    DO j = jstart, jend
      DO i = istart, iend
        s(i,j,k) = s(i,j,k) - ( zfx(i,j) - zfx(i-1,j) ) * rdx(j)*sqrtg_r(i,j,k)
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE xadv_rk_ppm

!==============================================================================
!============================================================================== 

SUBROUTINE yadv_rk_ppm( v, rdy, s, dt, ie, je, ke,             &
                        istart, iend, jstart, jend, kstart, kend, sqrtg_r )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine updates a scalar prognostic
!   variable 's' due to y-advection. The routine uses the
!   directional split method.
!
! Method:
!   Currently, the PPM advection scheme (Carpenter 1990) is used.
!
! Declarations:
!
!------------------------------------------------------------------------------

! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER          ,   INTENT (IN) ::  &
  ie, je, ke,     & ! dimensions of the input/output arrays
  istart, iend,   & ! computational start and end indices in i-direction
                    ! (west-east for lat/lon grid)
  jstart, jend,   & ! computational start and end indices in j-direction
                    ! (south-north for lat/lon grid)
  kstart, kend

REAL    (KIND=wp),   INTENT(IN)  ::  &
  dt                ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp),   INTENT(IN)  ::  &
  sqrtg_r(ie,je,kstart:kend), & ! 1 / square root of G
  rdy       (je, 2)         , & ! reciprocal y-grid spacing (1/( a cos(phi)dphi )
  v      (ie,je,kstart:kend)    ! v-velocity

! Array arguments with intent(out):

REAL    (KIND=wp),   INTENT(INOUT) ::  &
  s     (ie,je,kstart:kend)     ! scalar field updated due to y-advection

!------------------------------------------------------------------------------
! Local variables

INTEGER        :: &
  i, j, k, l        ! loop indices

REAL (KIND=wp) :: &
  zsgdy,          & !
  zcrn              ! Courant number for intermediate updating

! Local automatic arrays
REAL (KIND=wp) :: &
  sze(ie,je),     & ! scalar field to be transported at zone edges
  zfy(ie,je)        ! Flux in y-direction at cell boundary (v-point)

! End of header
!==============================================================================

  !----------------------------------------------------------------------------
  ! Section 1: Compute monotonic slopes and form fluxes in y - direction
  !----------------------------------------------------------------------------

  DO k = kstart, kend
    zfy (:,:) = 0.0_wp

    DO j = 2, je-2
      DO i = istart, iend
        sze(i,j) = ( 7.0_wp*( s(i,j+1,k)+s(i,j  ,k) )                  &
                              - ( s(i,j+2,k)+s(i,j-1,k) ) ) / 12.0_wp
      ENDDO
    ENDDO
  
    ! Compute courant number and velo-flux at zone edges (v-points)
    DO j = jstart-1, jend
      DO i = istart, iend
        ! Compute fluxes
        zcrn = dt*v(i,j,k)*rdy(j,2)
        IF ( zcrn > 0.0_wp ) THEN
          zcrn = zcrn * sqrtg_r(i,j,k)
          zfy(i,j) = dt*v(i,j,k) *  &
            ( sze(i,j) - zcrn*( sze(i,j  )-s(i,j  ,k) )              &
                          - zcrn*(1.0_wp-zcrn)*                  &
                                 ( sze(i,j-1)-2.0_wp*s(i,j  ,k)  &
                                 + sze(i,j  ) ) )
        ELSE IF ( zcrn < 0.0_wp ) THEN
          zcrn = zcrn * sqrtg_r(i,j+1,k)
          zfy(i,j) = dt*v(i,j,k) *  &
            ( sze(i,j) + zcrn*( sze(i,j  )-s(i,j+1,k) )              &
                          + zcrn*(1.0_wp+zcrn)*                  &
                                 ( sze(i,j  )-2.0_wp*s(i,j+1,k)  &
                                 + sze(i,j+1) ) )
        END IF
      ENDDO
    ENDDO

  !----------------------------------------------------------------------------
  ! Section 2: Update the scalar field
  !----------------------------------------------------------------------------
    DO j = jstart, jend
      DO i = istart, iend
        s(i,j,k) = s(i,j,k) - ( zfy(i,j) - zfy(i,j-1) ) * rdy(j,1)*sqrtg_r(i,j,k)
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE yadv_rk_ppm

!==============================================================================
!============================================================================== 

SUBROUTINE zadv_rk_ppm( wc, s, dt, ie, je, ke, ke1,        &
                        istart, iend, jstart, jend, sqrtg_r )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine updates a scalar prognostic
!   variable 's' due to vertical advection. The routine uses the
!   directional split method.
!
! Method:
!   Currently, the PPM advection scheme (Carpenter 1990) is used.
!
! Declarations:
!
!------------------------------------------------------------------------------

! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER           ,   INTENT (IN) ::  &
  ie, je,            & ! x- and y-dimension of the input/output arrays
  ke, ke1,           & ! z-dimension of the input/output arrays
  istart, iend,      & ! computational start and end indices in i-direction
                       ! (west-east for lat/lon grid)
  jstart, jend         ! computational start and end indices in j-direction
                       ! (south-north for lat/lon grid)

REAL    (KIND=wp),    INTENT(IN)  ::  &
  dt                   ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp),    INTENT(IN)  ::  &
  sqrtg_r(ie,je,ke), & ! 1 / square root of G
  wc(ie,je,ke1)        ! contravariant vertical velocity (at half levels)

! Array arguments with intent(inout):

REAL    (KIND=wp),    INTENT(INOUT) ::  &
  s(ie,je,ke)          ! scalar field updated due to z-advection

!------------------------------------------------------------------------------
! Local variables

INTEGER        :: &
  i, j, k, l        ! loop indices

REAL (KIND=wp) :: &
  zsgdz,          & !
  zcrn, zcrn_sgr    ! Courant number for intermediate updating

! Local automatic arrays
REAL (KIND=wp) :: &
  sze(ie,je,ke1), & ! scalar field to be transported at zone edges
  zfz(ie,je,ke1)    ! Flux in zeta-direction at cell boundary (w-point)

! End of header
!==============================================================================

  !----------------------------------------------------------------------------
  ! Section 1: Compute monotonic slopes and form fluxes in zeta - direction
  !----------------------------------------------------------------------------

  zfz (:,:,:) = 0.0_wp

  DO j = jstart, jend
    DO i = istart, iend
      sze(i,j,2)  = 0.5_wp*( s(i,j,2)+s(i,j,1) )
      sze(i,j,ke) = 0.5_wp*( s(i,j,ke)+s(i,j,ke-1) )
    ENDDO
  ENDDO
  DO k = 3, ke-1
    DO j = jstart, jend
      DO i = istart, iend
        sze(i,j,k) = ( 7.0_wp*( s(i,j,k  )+s(i,j,k-1) )                  &
                                - ( s(i,j,k+1)+s(i,j,k-2) ) ) / 12.0_wp
      ENDDO
    ENDDO
  ENDDO
  
  ! Compute courant number and velo-flux at zone edges (w-points)
  DO k = 2, ke
    DO j = jstart, jend
      DO i = istart, iend
        ! Compute fluxes
        zcrn = dt*wc(i,j,k)
        IF ( zcrn > 0.0_wp ) THEN
          zcrn_sgr = zcrn * sqrtg_r(i,j,k-1)
          IF ( k-1 >= 2 ) THEN
            zfz(i,j,k) = zcrn *  &
              ( sze(i,j,k) - zcrn_sgr*( sze(i,j,k  )-s(i,j,k-1) )      &
                              - zcrn_sgr*(1.0_wp-zcrn_sgr)*        &
                                 ( sze(i,j,k-1)-2.0_wp*s(i,j,k-1)  &
                                 + sze(i,j,k  ) ) )
          ELSE
            ! use upwind 1st order at upper boundary
            zfz(i,j,k) = zcrn * s(i,j,k-1)
          END IF
        ELSE IF ( zcrn < 0.0_wp ) THEN
          zcrn_sgr = zcrn * sqrtg_r(i,j,k)
          IF ( k+1 <= ke ) THEN
            zfz(i,j,k) = zcrn *  &
              ( sze(i,j,k) + zcrn_sgr*( sze(i,j,k  )-s(i,j,k  ) )      &
                              + zcrn_sgr*(1.0_wp+zcrn_sgr)*        &
                                 ( sze(i,j,k  )-2.0_wp*s(i,j,k  )  &
                                 + sze(i,j,k+1) ) )
          ELSE
            ! use upwind 1st order at lower boundary
            zfz(i,j,k) = zcrn * s(i,j,k  )
          END IF
        END IF
      ENDDO
    ENDDO
  ENDDO

  !----------------------------------------------------------------------------
  ! Section 2: Update the scalar field
  !----------------------------------------------------------------------------
  DO k = 1, ke
    DO j = jstart, jend
      DO i = istart, iend
        s(i,j,k) = s(i,j,k) - ( zfz(i,j,k+1) - zfz(i,j,k) ) * sqrtg_r(i,j,k)
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE zadv_rk_ppm

!==============================================================================
!============================================================================== 

SUBROUTINE xadv_pd_rk_cri_vanleer( u_frac, icr, rdx, s, dt, ie, je, &
                                   istart, iend, jstart, jend, sqrtg_r )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine updates a scalar prognostic
!   variable 's' due to x-advection. The routine uses the
!   directional split method.
!
! Method:
!   Currently, the 2nd order vanLeer advection scheme is used.
!   Courant number independent formulation  
!
! Declarations:
!
!------------------------------------------------------------------------------

! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER          ,   INTENT (IN) ::  &
  ie, je,         & ! dimensions of the input/output arrays
  istart, iend,   & ! computational start and end indices in i-direction
                    ! (west-east for lat/lon grid)
  jstart, jend      ! computational start and end indices in j-direction
                    ! (south-north for lat/lon grid)

REAL    (KIND=wp),   INTENT(IN)  ::  &
  dt                ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp),   INTENT(IN)  ::  &
  sqrtg_r(ie,je), & ! 1 / square root of G
  rdx(je),        & ! reciprocal x-grid spacing (1/( a cos(phi)dlam )
  u_frac(ie,je)     ! fractional u-velocity

INTEGER          ,   INTENT(IN)  ::  &
  icr(ie,je)        ! integer courant numbers

! Array arguments with intent(out):

REAL    (KIND=wp),   INTENT(INOUT) ::  &
  s(ie,je)          ! scalar field updated due to x-advection

!------------------------------------------------------------------------------
! Local variables

INTEGER        :: &
  i, j, l,        & ! loop indices
  ziif

REAL (KIND=wp) :: &
  zdql, zdqr, zr, & ! Parameters related to monotonic slopes
  zsgdx,          & !
  zcrn              ! Courant number for intermediate updating

! Local automatic arrays
REAL (KIND=wp) :: &
  zfx (ie,je),    & ! Flux in x-direction at cell boundary (u-point)
  zdxm(je),       & ! monotonic slopes for left and
  zdxp(je)          ! right grid point

! End of header
!==============================================================================

  !----------------------------------------------------------------------------
  ! Section 1: Compute monotonic slopes and form fluxes in x - direction
  !----------------------------------------------------------------------------

  zdxm(:)   = 0.0_wp
  zfx (:,:) = 0.0_wp

  i = istart-1
  DO j = jstart, jend
    ! Compute integer fluxes
    IF ( icr(i,j) >= 1 ) THEN
      ziif = i-icr(i,j)
    ELSE IF ( icr(i,j) <= -1 ) THEN
      ziif = i-icr(i,j)
    ELSE
      ziif = i
    END IF
    ! Compute monotonic slopes
    zdql = s(ziif  ,j) - s(ziif-1,j)
    zdqr = s(ziif+1,j) - s(ziif  ,j)
    zr   = zdqr / (zdql + eps_div)
    zdxm(j) = zdql * MAX(0.0_wp,   &
         MIN(2.0_wp, 2.0_wp*zr, 0.5_wp*(1.0_wp+zr)))
  ENDDO
  ! Compute courant number and velo-flux at zone edges (u-points)
  DO i = istart-1, iend
    zdxp(:) = 0.0_wp
    DO j = jstart, jend
      ! Compute integer fluxes
      IF ( icr(i,j) >= 1 ) THEN
        zsgdx = 1.0_wp / rdx(j) / sqrtg_r(i,j)
        DO l = 1, icr(i,j)
          ziif = i-l
          zfx(i,j) = zfx(i,j) + s(ziif+1,j) * zsgdx
        ENDDO 
      ELSE IF ( icr(i,j) <= -1 ) THEN
        zsgdx = 1.0_wp / rdx(j) / sqrtg_r(i+1,j)
        DO l = -1, icr(i,j), -1
          ziif = i-l
          zfx(i,j) = zfx(i,j) - s(ziif  ,j) * zsgdx
        ENDDO
      ELSE
        ziif = i
      END IF
      ! Compute monotonic slopes
      zdql = s(ziif+1,j) - s(ziif  ,j)
      zdqr = s(ziif+2,j) - s(ziif+1,j)
      zr   = zdqr / (zdql + eps_div)
      zdxp(j) = zdql * MAX(0.0_wp,   &
           MIN(2.0_wp, 2.0_wp*zr, 0.5_wp*(1.0_wp+zr)))
      ! Compute and add (fractional) fluxes
      zcrn = 0.5_wp*dt*u_frac(i,j)*rdx(j)
      zfx(i,j) = zfx(i,j) + dt * MAX( u_frac(i,j), 0.0_wp ) *                   &
                   ( s(ziif  ,j) + ( 0.5_wp - zcrn*sqrtg_r(i  ,j) ) * zdxm(j) ) &
                          + dt * MIN( u_frac(i,j), 0.0_wp ) *                   &
                   ( s(ziif+1,j) - ( 0.5_wp + zcrn*sqrtg_r(i+1,j) ) * zdxp(j) )
      zdxm(j) = zdxp(j)
    ENDDO
  ENDDO

  !----------------------------------------------------------------------------
  ! Section 2: Update the scalar field
  !----------------------------------------------------------------------------
  DO j = jstart, jend
    DO i = istart, iend
      s(i,j) = s(i,j) - ( zfx(i,j) - zfx(i-1,j) ) * rdx(j)*sqrtg_r(i,j)
    ENDDO
  ENDDO

END SUBROUTINE xadv_pd_rk_cri_vanleer

!==============================================================================
!============================================================================== 

SUBROUTINE yadv_pd_rk_cri_vanleer( v_frac, icr, rdy, s, dt, ie, je, &
                                   istart, iend, jstart, jend, sqrtg_r )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine updates a scalar prognostic
!   variable 's' due to y-advection. The routine uses the
!   directional split method.
!
! Method:
!   Currently, the 2nd order vanLeer advection scheme is used.
!   Courant number independent formulation  
!
! Declarations:
!
!------------------------------------------------------------------------------

! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER          ,   INTENT (IN) ::  &
  ie, je,         & ! dimensions of the input/output arrays
  istart, iend,   & ! computational start and end indices in i-direction
                    ! (west-east for lat/lon grid)
  jstart, jend      ! computational start and end indices in j-direction
                    ! (south-north for lat/lon grid)

REAL    (KIND=wp),   INTENT(IN)  ::  &
  dt                ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp),   INTENT(IN)  ::  &
  sqrtg_r(ie,je), & ! 1 / square root of G
  rdy(je,2),      & ! reciprocal y-grid spacing (1/( a cos(phi)dphi )
  v_frac(ie,je)     ! fractional v-velocity

INTEGER          ,   INTENT(IN)  ::  &
  icr(ie,je)        ! integer courant numbers

! Array arguments with intent(out):

REAL    (KIND=wp),   INTENT(INOUT) ::  &
  s(ie,je)          ! scalar field updated due to y-advection

!------------------------------------------------------------------------------
! Local variables

INTEGER        :: &
  i, j, l,        & ! loop indices
  zjif

REAL (KIND=wp) :: &
  zdql, zdqr, zr, & ! Parameters related to monotonic slopes
  zsgdy,          & !
  zcrn              ! Courant number for intermediate updating

! Local automatic arrays
REAL (KIND=wp) :: &
  zfy (ie,je),    & ! Flux in y-direction at cell boundary (v-point)
  zdym(ie),       & ! monotonic slopes for left and
  zdyp(ie)          ! right grid point

! End of header
!==============================================================================

  !----------------------------------------------------------------------------
  ! Section 1: Compute monotonic slopes and form fluxes in y - direction
  !----------------------------------------------------------------------------

  zdym(:)   = 0.0_wp
  zfy (:,:) = 0.0_wp

  j = jstart-1
  DO i = istart, iend
    ! Compute integer fluxes
    IF ( icr(i,j) >= 1 ) THEN
      zjif = j-icr(i,j)
    ELSE IF ( icr(i,j) <= -1 ) THEN
      zjif = j-icr(i,j)
    ELSE
      zjif = j
    END IF
    ! Compute monotonic slopes
    zdql = s(i,zjif  ) - s(i,zjif-1)
    zdqr = s(i,zjif+1) - s(i,zjif  )
    zr   = zdqr / (zdql + eps_div)
    zdym(i) = zdql * MAX(0.0_wp,   &
         MIN(2.0_wp, 2.0_wp*zr, 0.5_wp*(1.0_wp+zr)))
  ENDDO
  ! Compute courant number and velo-flux at zone edges (v-points)
  DO j = jstart-1, jend
    zdyp(:) = 0.0_wp
    DO i = istart, iend
      ! Compute integer fluxes
      IF ( icr(i,j) >= 1 ) THEN
        zsgdy = 1.0_wp / rdy(j,2) / sqrtg_r(i,j)
        DO l = 1, icr(i,j) 
          zjif = j-l
          zfy(i,j) = zfy(i,j) + s(i,zjif+1) * zsgdy
        ENDDO 
      ELSE IF ( icr(i,j) <= -1 ) THEN
        zsgdy = 1.0_wp / rdy(j,2) / sqrtg_r(i,j+1)
        DO l = -1, icr(i,j), -1
          zjif = j-l
          zfy(i,j) = zfy(i,j) - s(i,zjif  ) * zsgdy
        ENDDO
      ELSE
        zjif = j
      END IF
      ! Compute monotonic slopes
      zdql = s(i,zjif+1) - s(i,zjif  )
      zdqr = s(i,zjif+2) - s(i,zjif+1)
      zr   = zdqr / (zdql + eps_div)
      zdyp(i) = zdql * MAX(0.0_wp,   &
           MIN(2.0_wp, 2.0_wp*zr, 0.5_wp*(1.0_wp+zr)))
      ! Compute and add (fractional) fluxes
      zcrn = 0.5_wp*dt*v_frac(i,j)*rdy(j,2)
      zfy(i,j) = zfy(i,j) + dt * MAX( v_frac(i,j), 0.0_wp ) *                   &
                   ( s(i,zjif  ) + ( 0.5_wp - zcrn*sqrtg_r(i,j  ) ) * zdym(i) ) &
                          + dt * MIN( v_frac(i,j), 0.0_wp ) *                   &
                   ( s(i,zjif+1) - ( 0.5_wp + zcrn*sqrtg_r(i,j+1) ) * zdyp(i) )
      zdym(i) = zdyp(i)
    ENDDO
  ENDDO

  !----------------------------------------------------------------------------
  ! Section 2: Update the scalar field
  !----------------------------------------------------------------------------
  DO j = jstart, jend
    DO i = istart, iend
      s(i,j) = s(i,j) - ( zfy(i,j) - zfy(i,j-1) ) * rdy(j,1)*sqrtg_r(i,j)
    ENDDO
  ENDDO

END SUBROUTINE yadv_pd_rk_cri_vanleer

!==============================================================================
!============================================================================== 

SUBROUTINE zadv_pd_rk_cri_vanleer( wc_frac, icr, s, dt, ie, je, ke, ke1, &
                                   istart, iend, jstart, jend, sqrtg_r )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine updates a scalar prognostic
!   variable 's' due to vertical advection. The routine uses the
!   directional split method.
!
! Method:
!   Currently, the 2nd order vanLeer advection scheme is used.
!   Courant number independent formulation  
!
! Declarations:
!
!------------------------------------------------------------------------------

! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER          ,   INTENT (IN) ::  &
  ie, je,            & ! x- and y-dimension of the input/output arrays
  ke, ke1,           & ! z-dimension of the input/output arrays
  istart, iend,      & ! computational start and end indices in i-direction
                       ! (west-east for lat/lon grid)
  jstart, jend         ! computational start and end indices in j-direction
                       ! (south-north for lat/lon grid)

REAL    (KIND=wp),   INTENT(IN)  ::  &
  dt                   ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp),   INTENT(IN)  ::  &
  sqrtg_r(ie,je,ke), & ! 1 / square root of G
  wc_frac(ie,je,ke1)   ! fractional contravariant vertical velocity
                       ! (at half levels)

INTEGER          ,   INTENT(IN)  ::  &
  icr(ie,je,ke1)       ! integer courant numbers (at half levels)

! Array arguments with intent(inout):

REAL    (KIND=wp),   INTENT(INOUT) ::  &
  s(ie,je,ke)          ! scalar field updated due to z-advection

!------------------------------------------------------------------------------
! Local variables

INTEGER        :: &
  i, j, k, l,     & ! loop indices
  zkif

REAL (KIND=wp) :: &
  zdqtop, zdqbot, & ! Parameters related to
  zr,             & ! monotonic slopes
  zsgdz,          & !
  zcrn              ! Courant number for intermediate updating

! Local automatic arrays
REAL (KIND=wp) :: &
  zfz (ie,je,ke1),& ! Flux in zeta-direction at cell boundary (w-point)
  zdzm(ie,je),    & ! monotonic slopes for upper and
  zdzp(ie,je)       ! lower grid point (resp. the w-point)

! End of header
!==============================================================================

  !----------------------------------------------------------------------------
  ! Section 1: Compute monotonic slopes and form fluxes in zeta - direction
  !----------------------------------------------------------------------------

  zdzm(:,:) = 0.0_wp
  zfz (:,:,:) = 0.0_wp

  DO j = jstart, jend
    DO i = istart, iend
      zdzm(i,j) = s(i,j,2) - s(i,j,1)
    ENDDO
  ENDDO
  ! Compute courant number and velo-flux at zone edges (w-points)
  DO k = 2, ke
    zdzp(:,:) = 0.0_wp
    IF ( k == ke ) THEN
      DO j = jstart, jend
        DO i = istart, iend
          ! Compute integer fluxes
          IF ( icr(i,j,k) >= 1 ) THEN
            zsgdz = 1.0_wp / sqrtg_r(i,j,k-1)
            DO l = 1, icr(i,j,k) 
              zkif = k-l
              zfz(i,j,k) = zfz(i,j,k) + s(i,j,zkif  ) * zsgdz
            ENDDO 
          ELSE
            zkif = k
          END IF
          zdzp(i,j) = s(i,j,zkif  ) - s(i,j,zkif-1)
          ! Compute fluxes and add (fractional) fluxes
          zcrn = 0.5_wp*dt*wc_frac(i,j,k)
          zfz(i,j,k) = zfz(i,j,k) + dt * MAX( wc_frac(i,j,k), 0.0_wp ) *        &
             ( s(i,j,zkif-1) + ( 0.5_wp - zcrn*sqrtg_r(i,j,k-1) ) * zdzm(i,j) ) &
                                  + dt * MIN( wc_frac(i,j,k), 0.0_wp ) *        &
             ( s(i,j,zkif  ) - ( 0.5_wp + zcrn*sqrtg_r(i,j,k  ) ) * zdzp(i,j) )
        ENDDO
      ENDDO
    ELSE
      DO j = jstart, jend
        DO i = istart, iend
          ! Compute integer fluxes
          IF ( icr(i,j,k) >= 1 ) THEN
            zsgdz = 1.0_wp / sqrtg_r(i,j,k-1)
            DO l = 1, icr(i,j,k) 
              zkif = k-l
              zfz(i,j,k) = zfz(i,j,k) + s(i,j,zkif  ) * zsgdz
            ENDDO 
          ELSE IF ( icr(i,j,k) <= -1 ) THEN
            zsgdz = 1.0_wp / sqrtg_r(i,j,k)
            DO l = -1, icr(i,j,k), -1
              zkif = k-l
              zfz(i,j,k) = zfz(i,j,k) - s(i,j,zkif-1) * zsgdz
            ENDDO
          ELSE
            zkif = k
          END IF
          ! Compute monotonic slopes
          zdqtop   = s(i,j,zkif  ) - s(i,j,zkif-1)
          zdqbot   = s(i,j,zkif+1) - s(i,j,zkif  )
          zr       = zdqbot / (zdqtop + eps_div)
          zdzp(i,j) = zdqtop * MAX(0.0_wp,   &
               MIN(2.0_wp, 2.0_wp*zr, 0.5_wp*(1.0_wp+zr)))
          ! Compute and add (fractional) fluxes
          zcrn = 0.5_wp*dt*wc_frac(i,j,k)
          zfz(i,j,k) = zfz(i,j,k) + dt * MAX( wc_frac(i,j,k), 0.0_wp ) *        &
             ( s(i,j,zkif-1) + ( 0.5_wp - zcrn*sqrtg_r(i,j,k-1) ) * zdzm(i,j) ) &
                                  + dt * MIN( wc_frac(i,j,k), 0.0_wp ) *        &
             ( s(i,j,zkif  ) - ( 0.5_wp + zcrn*sqrtg_r(i,j,k  ) ) * zdzp(i,j) )
          zdzm(i,j) = zdzp(i,j)
        ENDDO
      ENDDO
    ENDIF
  ENDDO

  !----------------------------------------------------------------------------
  ! Section 2: Update the scalar field
  !----------------------------------------------------------------------------
  DO k = 1, ke
    DO j = jstart, jend
      DO i = istart, iend
        s(i,j,k) = s(i,j,k) - ( zfz(i,j,k+1) - zfz(i,j,k) ) * sqrtg_r(i,j,k)
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE zadv_pd_rk_cri_vanleer

!==============================================================================
!============================================================================== 

SUBROUTINE xadv_pd_rk_vanleer( u, rdx, s, dt, ie, je, ke,            &
                               istart, iend, jstart, jend, kstart, kend, sqrtg_r )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine updates a scalar prognostic
!   variable 's' due to x-advection. The routine uses the
!   directional split method.
!
! Method:
!   Currently, the 2nd order vanLeer advection scheme is used.
!
! Declarations:
!
!------------------------------------------------------------------------------

! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER          ,   INTENT (IN) ::  &
  ie, je, ke,     & ! dimensions of the input/output arrays
  istart, iend,   & ! computational start and end indices in i-direction
                    ! (west-east for lat/lon grid)
  jstart, jend,   & ! computational start and end indices in j-direction
                    ! (south-north for lat/lon grid)
  kstart, kend

REAL    (KIND=wp),   INTENT(IN)  ::  &
  dt                ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp),   INTENT(IN)  ::  &
  sqrtg_r(ie,je,kstart:kend), & ! 1 / square root of G
  rdx       (je)            , & ! reciprocal x-grid spacing (1/( a cos(phi)dlam )
  u      (ie,je,kstart:kend)    ! u-velocity

! Array arguments with intent(out):

REAL    (KIND=wp),   INTENT(INOUT) ::  &
  s      (ie,je,kstart:kend)    ! scalar field updated due to x-advection

!------------------------------------------------------------------------------
! Local variables

INTEGER        :: &
  i, j, k, l        ! loop indices

REAL (KIND=wp) :: &
  zdql, zdqr, zr, & ! Parameters related to monotonic slopes
  zsgdx,          & !
  zcrn              ! Courant number for intermediate updating

! Local automatic arrays
REAL (KIND=wp) :: &
  zfx (ie,je),    & ! Flux in x-direction at cell boundary (u-point)
  zdxm(je),       & ! monotonic slopes for left and
  zdxp(je)          ! right grid point

! End of header
!==============================================================================

  !----------------------------------------------------------------------------
  ! Section 1: Compute monotonic slopes and form fluxes in x - direction
  !----------------------------------------------------------------------------

  DO k = kstart, kend
    zdxm(:)   = 0.0_wp
    zfx (:,:) = 0.0_wp

    i = istart-1
    DO j = jstart, jend
      ! Compute monotonic slopes
      zdql = s(i  ,j,k) - s(i-1,j,k)
      zdqr = s(i+1,j,k) - s(i  ,j,k)
      zr   = zdqr / (zdql + eps_div)
      zdxm(j) = zdql * MAX(0.0_wp,   &
           MIN(2.0_wp, 2.0_wp*zr, 0.5_wp*(1.0_wp+zr)))
    ENDDO
    ! Compute courant number and velo-flux at zone edges (u-points)
    DO i = istart-1, iend
      zdxp(:) = 0.0_wp
      DO j = jstart, jend
        ! Compute monotonic slopes
        zdql = s(i+1,j,k) - s(i  ,j,k)
        zdqr = s(i+2,j,k) - s(i+1,j,k)
        zr   = zdqr / (zdql + eps_div)
        zdxp(j) = zdql * MAX(0.0_wp,   &
             MIN(2.0_wp, 2.0_wp*zr, 0.5_wp*(1.0_wp+zr)))
        ! Compute fluxes
        zcrn = 0.5_wp*dt*u(i,j,k)*rdx(j)
        zfx(i,j) = dt * MAX( u(i,j,k), 0.0_wp ) *                                 &
                    ( s(i  ,j,k) + ( 0.5_wp - zcrn*sqrtg_r(i  ,j,k) ) * zdxm(j))  &
                 + dt * MIN( u(i,j,k), 0.0_wp ) *                                 &
                    ( s(i+1,j,k) - ( 0.5_wp + zcrn*sqrtg_r(i+1,j,k) ) * zdxp(j))
        zdxm(j) = zdxp(j)
      ENDDO
    ENDDO

  !----------------------------------------------------------------------------
  ! Section 2: Update the scalar field
  !----------------------------------------------------------------------------
    DO j = jstart, jend
      DO i = istart, iend
        s(i,j,k) = s(i,j,k) - ( zfx(i,j) - zfx(i-1,j) ) * rdx(j)*sqrtg_r(i,j,k)
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE xadv_pd_rk_vanleer

!==============================================================================
!============================================================================== 

SUBROUTINE yadv_pd_rk_vanleer( v, rdy, s, dt, ie, je, ke,            &
                               istart, iend, jstart, jend, kstart, kend, sqrtg_r )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine updates a scalar prognostic
!   variable 's' due to y-advection. The routine uses the
!   directional split method.
!
! Method:
!   Currently, the 2nd order vanLeer advection scheme is used.
!
! Declarations:
!
!------------------------------------------------------------------------------

! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER          ,   INTENT (IN) ::  &
  ie, je, ke,     & ! dimensions of the input/output arrays
  istart, iend,   & ! computational start and end indices in i-direction
                    ! (west-east for lat/lon grid)
  jstart, jend,   & ! computational start and end indices in j-direction
                    ! (south-north for lat/lon grid)
  kstart, kend

REAL    (KIND=wp),   INTENT(IN)  ::  &
  dt                ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp),   INTENT(IN)  ::  &
  sqrtg_r(ie,je,kstart:kend), & ! 1 / square root of G
  rdy       (je,2),           & ! reciprocal y-grid spacing (1/( a cos(phi)dphi )
  v      (ie,je,kstart:kend)    ! v-velocity

! Array arguments with intent(out):

REAL    (KIND=wp),   INTENT(INOUT) ::  &
  s      (ie,je,kstart:kend)    ! scalar field updated due to y-advection

!------------------------------------------------------------------------------
! Local variables

INTEGER        :: &
  i, j, k, l        ! loop indices

REAL (KIND=wp) :: &
  zdql, zdqr, zr, & ! Parameters related to monotonic slopes
  zsgdy,          & !
  zcrn              ! Courant number for intermediate updating

! Local automatic arrays
REAL (KIND=wp) :: &
  zfy (ie,je),    & ! Flux in y-direction at cell boundary (v-point)
  zdym(ie),       & ! monotonic slopes for left and
  zdyp(ie)          ! right grid point

! End of header
!==============================================================================

  !----------------------------------------------------------------------------
  ! Section 1: Compute monotonic slopes and form fluxes in y - direction
  !----------------------------------------------------------------------------

  DO k = kstart, kend
    zdym(:)   = 0.0_wp
    zfy (:,:) = 0.0_wp

    j = jstart-1
    DO i = istart, iend
      ! Compute monotonic slopes
      zdql = s(i,j  ,k) - s(i,j-1,k)
      zdqr = s(i,j+1,k) - s(i,j  ,k)
      zr   = zdqr / (zdql + eps_div)
      zdym(i) = zdql * MAX(0.0_wp,   &
           MIN(2.0_wp, 2.0_wp*zr, 0.5_wp*(1.0_wp+zr)))
    ENDDO
    ! Compute courant number and velo-flux at zone edges (v-points)
    DO j = jstart-1, jend
      zdyp(:) = 0.0_wp
      DO i = istart, iend
        ! Compute monotonic slopes
        zdql = s(i,j+1,k) - s(i,j  ,k)
        zdqr = s(i,j+2,k) - s(i,j+1,k)
        zr   = zdqr / (zdql + eps_div)
        zdyp(i) = zdql * MAX(0.0_wp,   &
             MIN(2.0_wp, 2.0_wp*zr, 0.5_wp*(1.0_wp+zr)))
        ! Compute fluxes
        zcrn = 0.5_wp*dt*v(i,j,k)*rdy(j,2)
        zfy(i,j) = dt * MAX( v(i,j,k), 0.0_wp ) *                                 &
                    ( s(i,j  ,k) + ( 0.5_wp - zcrn*sqrtg_r(i,j  ,k) ) * zdym(i))  &
                 + dt * MIN( v(i,j,k), 0.0_wp ) *                                 &
                    ( s(i,j+1,k) - ( 0.5_wp + zcrn*sqrtg_r(i,j+1,k) ) * zdyp(i))
        zdym(i) = zdyp(i)
      ENDDO
    ENDDO

  !----------------------------------------------------------------------------
  ! Section 2: Update the scalar field
  !----------------------------------------------------------------------------
    DO j = jstart, jend
      DO i = istart, iend
        s(i,j,k) = s(i,j,k) - ( zfy(i,j) - zfy(i,j-1) ) * rdy(j,1)*sqrtg_r(i,j,k)
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE yadv_pd_rk_vanleer

!==============================================================================
!============================================================================== 

SUBROUTINE zadv_pd_rk_vanleer( wc, s, dt, ie, je, ke, ke1,        &
                               istart, iend, jstart, jend, sqrtg_r )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine updates a scalar prognostic
!   variable 's' due to vertical advection. The routine uses the
!   directional split method.
!
! Method:
!   Currently, the 2nd order vanLeer advection scheme is used.
!
! Declarations:
!
!------------------------------------------------------------------------------

! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER           ,   INTENT (IN) ::  &
  ie, je,            & ! x- and y-dimension of the input/output arrays
  ke, ke1,           & ! z-dimension of the input/output arrays
  istart, iend,      & ! computational start and end indices in i-direction
                       ! (west-east for lat/lon grid)
  jstart, jend         ! computational start and end indices in j-direction
                       ! (south-north for lat/lon grid)

REAL    (KIND=wp),    INTENT(IN)  ::  &
  dt                   ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp),    INTENT(IN)  ::  &
  sqrtg_r(ie,je,ke), & ! 1 / square root of G
  wc(ie,je,ke1)        ! contravariant vertical velocity (at half levels)

! Array arguments with intent(inout):

REAL    (KIND=wp),    INTENT(INOUT) ::  &
  s(ie,je,ke)          ! scalar field updated due to z-advection

!------------------------------------------------------------------------------
! Local variables

INTEGER        :: &
  i, j, k, l        ! loop indices

REAL (KIND=wp) :: &
  zdqtop, zdqbot, & ! Parameters related to
  zr,             & ! monotonic slopes
  zsgdz,          & !
  zcrn              ! Courant number for intermediate updating

! Local automatic arrays
REAL (KIND=wp) :: &
  zfz (ie,je,ke1),& ! Flux in zeta-direction at cell boundary (w-point)
  zdzm(ie,je),    & ! monotonic slopes for upper and
  zdzp(ie,je)       ! lower grid point (resp. the w-point)

! End of header
!==============================================================================

  !----------------------------------------------------------------------------
  ! Section 1: Compute monotonic slopes and form fluxes in zeta - direction
  !----------------------------------------------------------------------------

  zdzm(:,:) = 0.0_wp
  zfz (:,:,:) = 0.0_wp

  DO j = jstart, jend
    DO i = istart, iend
      zdzm(i,j) = s(i,j,2) - s(i,j,1)
    ENDDO
  ENDDO
  ! Compute courant number and velo-flux at zone edges (w-points)
  DO k = 2, ke
    zdzp(:,:) = 0.0_wp
    IF ( k == ke ) THEN
      DO j = jstart, jend
        DO i = istart, iend
          zdzp(i,j) = s(i,j,k  ) - s(i,j,k-1)
          ! Compute fluxes
          zcrn = 0.5_wp*dt*wc(i,j,k)
          zfz(i,j,k) = dt * MAX( wc(i,j,k), 0.0_wp ) *                        &
             ( s(i,j,k-1) + ( 0.5_wp - zcrn*sqrtg_r(i,j,k-1) ) * zdzm(i,j) )  &
                     + dt * MIN( wc(i,j,k), 0.0_wp ) *                        &
             ( s(i,j,k  ) - ( 0.5_wp + zcrn*sqrtg_r(i,j,k  ) ) * zdzp(i,j) )
        ENDDO
      ENDDO
    ELSE
      DO j = jstart, jend
        DO i = istart, iend
          ! Compute monotonic slopes
          zdqtop   = s(i,j,k  ) - s(i,j,k-1)
          zdqbot   = s(i,j,k+1) - s(i,j,k  )
          zr       = zdqbot / (zdqtop + eps_div)
          zdzp(i,j) = zdqtop * MAX(0.0_wp,   &
               MIN(2.0_wp, 2.0_wp*zr, 0.5_wp*(1.0_wp+zr)))
          ! Compute fluxes
          zcrn = 0.5_wp*dt*wc(i,j,k)
          zfz(i,j,k) = dt * MAX( wc(i,j,k), 0.0_wp ) *                        &
             ( s(i,j,k-1) + ( 0.5_wp - zcrn*sqrtg_r(i,j,k-1) ) * zdzm(i,j) )  &
                     + dt * MIN( wc(i,j,k), 0.0_wp ) *                        &
             ( s(i,j,k  ) - ( 0.5_wp + zcrn*sqrtg_r(i,j,k  ) ) * zdzp(i,j) )
          zdzm(i,j) = zdzp(i,j)
        ENDDO
      ENDDO
    ENDIF
  ENDDO

  !----------------------------------------------------------------------------
  ! Section 2: Update the scalar field
  !----------------------------------------------------------------------------
  DO k = 1, ke
    DO j = jstart, jend
      DO i = istart, iend
        s(i,j,k) = s(i,j,k) - ( zfz(i,j,k+1) - zfz(i,j,k) ) * sqrtg_r(i,j,k)
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE zadv_pd_rk_vanleer

!==============================================================================
!==============================================================================

SUBROUTINE init_bott_coeffs

  abc(:,:,:) = 0.0_wp

! To avoid too many divisions during the forecast, the initial coefficients
! are modified

  abc(2,0,-1:1) = (/ 1.0_wp, -26.0_wp, 1.0_wp /)
! abc(2,0,-1:1) = -1.0_wp/24.0_wp * abc(2,0,-1:1)  original
  abc(2,0,-1:1) = -1.0_wp/48.0_wp * abc(2,0,-1:1)

  abc(2,1,-1:1) = (/ -1.0_wp, 0.0_wp, 1.0_wp /)
! abc(2,1,-1:1) =  1.0_wp/ 2.0_wp * abc(2,1,-1:1)  original
  abc(2,1,-1:1) =  1.0_wp/16.0_wp * abc(2,1,-1:1)

  abc(2,2,-1:1) = (/ 1.0_wp, -2.0_wp, 1.0_wp /)
! abc(2,2,-1:1) =  1.0_wp/ 2.0_wp * abc(2,2,-1:1)  original
  abc(2,2,-1:1) =  1.0_wp/48.0_wp * abc(2,2,-1:1)


  abc(4,0,-2:2) = (/ 9.0_wp, -116.0_wp, 2134.0_wp, -116.0_wp, 9.0_wp /)
! abc(4,0,-2:2) =  1.0_wp/1920.0_wp * abc(4,0,-2:2) original
  abc(4,0,-2:2) =  1.0_wp/3840.0_wp * abc(4,0,-2:2)

  abc(4,1,-2:2) = (/ 5.0_wp, -34.0_wp, 0.0_wp, 34.0_wp, -5.0_wp /)
! abc(4,1,-2:2) =  1.0_wp/ 48.0_wp  * abc(4,1,-2:2) original
  abc(4,1,-2:2) =  1.0_wp/384.0_wp  * abc(4,1,-2:2)

  abc(4,2,-2:2) = (/ -3.0_wp, 36.0_wp, -66.0_wp, 36.0_wp, -3.0_wp /)
! abc(4,2,-2:2) =  1.0_wp/  48.0_wp * abc(4,2,-2:2) original
  abc(4,2,-2:2) =  1.0_wp/1152.0_wp * abc(4,2,-2:2)

  abc(4,3,-2:2) = (/ -1.0_wp, 2.0_wp, 0.0_wp, -2.0_wp, 1.0_wp /)
! abc(4,3,-2:2) =  1.0_wp/ 12.0_wp  * abc(4,3,-2:2) original
  abc(4,3,-2:2) =  1.0_wp/768.0_wp  * abc(4,3,-2:2)

  abc(4,4,-2:2) = (/ 1.0_wp, -4.0_wp, 6.0_wp, -4.0_wp, 1.0_wp /)
! abc(4,4,-2:2) =  1.0_wp/  12.0_wp * abc(4,4,-2:2) original
  abc(4,4,-2:2) =  1.0_wp/1920.0_wp * abc(4,4,-2:2)

END SUBROUTINE init_bott_coeffs

!==============================================================================
!==============================================================================

SUBROUTINE xadv_pd_rk_cri_bott( u_frac, icr, rdx, s, dt, ie, je, &
                                istart, iend, jstart, jend, sqrtg_r )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine updates a scalar prognostic
!   variable 's' due to x-advection. The routine uses the
!   directional split method.
!
! Method:
!   Currently, the Bott advection scheme (Bott 1989) is used.
!   Courant number independent formulation  
!
!------------------------------------------------------------------------------
  
!------------------------------------------------------------------------------
! Declarations:
!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &
  y_scalar_advect  ! type of scalar advection scheme (for qx, TKE and other tracers)

! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER          ,   INTENT (IN) ::  &
  ie, je,         & ! dimensions of the input/output arrays
  istart, iend,   & ! computational start and end indices in i-direction
                    ! (west-east for lat/lon grid)
  jstart, jend      ! computational start and end indices in j-direction
                    ! (south-north for lat/lon grid)

REAL    (KIND=wp),   INTENT(IN)  ::  &
  dt                ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp),   INTENT(IN)  ::  &
  sqrtg_r(ie,je), & ! 1 / square root of G
  rdx(je),        & ! reciprocal x-grid spacing (1/( a cos(phi)dlam )
  u_frac(ie,je)     ! fractional u-velocity

INTEGER          ,   INTENT(IN)  ::  &
  icr(ie,je)        ! integer courant numbers

! Array arguments with intent(out):

REAL    (KIND=wp),   INTENT(INOUT) ::  &
  s(ie,je)          ! scalar field updated due to x-advection

!------------------------------------------------------------------------------
! Local variables

INTEGER        ::    &
  i, j, l,           & ! loop indices
  ziif

REAL (KIND=wp) ::    &
  zdx,               & ! dx
  zsg,               & !
  zcrn                 ! Courant number for intermediate updating

! Local automatic arrays
REAL (KIND=wp) ::    &
  zalpha(ie,je,0:4), & !
  zfcr1(ie,je),      & !
  zfi_p(ie,je),      & !
  zfi_m(ie,je),      & !  
  zfx(ie,je)           ! Flux in x-direction at cell boundary (u-point)

! End of header
!==============================================================================

  !----------------------------------------------------------------------------
  ! Section 1: Compute monotonic slopes and form fluxes in x - direction
  !----------------------------------------------------------------------------

  zfx (:,:) = 0.0_wp

  SELECT CASE( TRIM(y_scalar_advect) )
    
  CASE( "BOTT2", "BOTT2_STRANG", "BOTT2_XYZYX", "BOTT2_STRANG_B" )

    ! Compute polynomials and fluxes for weighting (Cr=1)
    DO j = jstart, jend
      DO i = 2, ie-1
        zalpha(i,j,0) = abc(2,0,-1) * s(i-1,j)    &
                      + abc(2,0, 0) * s(i  ,j)    &
                      + abc(2,0, 1) * s(i+1,j)

        zalpha(i,j,1) = abc(2,1,-1) * s(i-1,j)    &
                      + abc(2,1, 1) * s(i+1,j)

        zalpha(i,j,2) = abc(2,2,-1) * s(i-1,j)    &
                      + abc(2,2, 0) * s(i  ,j)    &
                      + abc(2,2, 1) * s(i+1,j)
      ENDDO
    ENDDO
    
    ! Compute courant number and velo-flux at zone edges (u-points)
    DO j = jstart, jend
      DO i = istart-2, iend+1
        ! Compute integer fluxes
        IF ( icr(i,j) >= 1 ) THEN
          zsg = 1.0_wp / sqrtg_r(i,j)
          DO l = 1, icr(i,j)
            ziif = i-l
            zfx(i,j) = zfx(i,j) + s(ziif+1,j) * zsg
          ENDDO 
        ELSE IF ( icr(i,j) <= -1 ) THEN
          zsg = 1.0_wp / sqrtg_r(i+1,j)
          DO l = -1, icr(i,j), -1
            ziif = i-l
            zfx(i,j) = zfx(i,j) - s(ziif  ,j) * zsg
          ENDDO
        ELSE
          ziif = i
        END IF
        ! Compute integral fluxes
        zfi_p(i,j) = 0.0_wp
        zfi_m(i,j) = 0.0_wp
        zcrn = dt*u_frac(i,j)*rdx(j)
        IF ( zcrn > 0.0_wp ) THEN
          zcrn = 1.0_wp - 2.0_wp * zcrn * sqrtg_r(i,j)
          zfi_p(i,j) = zalpha(ziif  ,j,0) * ( 1.0_wp - zcrn )     &
                     + zalpha(ziif  ,j,1) * ( 1.0_wp - zcrn**2 )  &
                     + zalpha(ziif  ,j,2) * ( 1.0_wp - zcrn**3 )
          zfi_p(i,j) = MAX( 0.0_wp, zfi_p(i,j) )
        ELSE IF ( zcrn < 0.0_wp ) THEN
          zcrn = 1.0_wp + 2.0_wp * zcrn * sqrtg_r(i+1,j)
          zfi_m(i,j) = zalpha(ziif+1,j,0) * ( 1.0_wp - zcrn )     &
                     - zalpha(ziif+1,j,1) * ( 1.0_wp - zcrn**2 )  &
                     + zalpha(ziif+1,j,2) * ( 1.0_wp - zcrn**3 )
          zfi_m(i,j) = MAX( 0.0_wp, zfi_m(i,j) )
        END IF
        zfcr1(i,j) = ( zalpha(ziif,j,0) + zalpha(ziif,j,2) )  &
                   * 2.0_wp
      ENDDO
    ENDDO

  CASE( "BOTT4", "BOTT4_STRANG", "BOTT4_XYZYX", "BOTT4_STRANG_B" )
    
    ! Compute polynomials and fluxes for weighting (Cr=1)
    DO j = jstart, jend
      DO i = 3, ie-2
        zalpha(i,j,0) = abc(4,0,-2) * s(i-2,j)    &
                      + abc(4,0,-1) * s(i-1,j)    &
                      + abc(4,0, 0) * s(i  ,j)    &
                      + abc(4,0, 1) * s(i+1,j)    &
                      + abc(4,0, 2) * s(i+2,j)

        zalpha(i,j,1) = abc(4,1,-2) * s(i-2,j)    &
                      + abc(4,1,-1) * s(i-1,j)    &
                      + abc(4,1, 1) * s(i+1,j)    &
                      + abc(4,1, 2) * s(i+2,j)

        zalpha(i,j,2) = abc(4,2,-2) * s(i-2,j)    &
                      + abc(4,2,-1) * s(i-1,j)    &
                      + abc(4,2, 0) * s(i  ,j)    &
                      + abc(4,2, 1) * s(i+1,j)    &
                      + abc(4,2, 2) * s(i+2,j)

        zalpha(i,j,3) = abc(4,3,-2) * s(i-2,j)    &
                      + abc(4,3,-1) * s(i-1,j)    &
                      + abc(4,3, 1) * s(i+1,j)    &
                      + abc(4,3, 2) * s(i+2,j)

        zalpha(i,j,4) = abc(4,4,-2) * s(i-2,j)    &
                      + abc(4,4,-1) * s(i-1,j)    &
                      + abc(4,4, 0) * s(i  ,j)    &
                      + abc(4,4, 1) * s(i+1,j)    &
                      + abc(4,4, 2) * s(i+2,j)
      ENDDO
    ENDDO
    
    ! Compute courant number and velo-flux at zone edges (u-points)
    DO j = jstart, jend
      DO i = istart-2, iend+1
        ! Compute integer fluxes
        IF ( icr(i,j) >= 1 ) THEN
          zsg = 1.0_wp / sqrtg_r(i,j)
          DO l = 1, icr(i,j)
            ziif = i-l
            zfx(i,j) = zfx(i,j) + s(ziif+1,j) * zsg
          ENDDO 
        ELSE IF ( icr(i,j) <= -1 ) THEN
          zsg = 1.0_wp / sqrtg_r(i+1,j)
          DO l = -1, icr(i,j), -1
            ziif = i-l
            zfx(i,j) = zfx(i,j) - s(ziif  ,j) * zsg
          ENDDO
        ELSE
          ziif = i
        END IF
        ! Compute integral fluxes
        zfi_p(i,j) = 0.0_wp
        zfi_m(i,j) = 0.0_wp
        zcrn = dt*u_frac(i,j)*rdx(j)
        IF ( zcrn > 0.0_wp ) THEN
          zcrn = 1.0_wp - 2.0_wp * zcrn * sqrtg_r(i,j)
          zfi_p(i,j) = zalpha(ziif  ,j,0) * ( 1.0_wp - zcrn )     &
                     + zalpha(ziif  ,j,1) * ( 1.0_wp - zcrn**2 )  &
                     + zalpha(ziif  ,j,2) * ( 1.0_wp - zcrn**3 )  &
                     + zalpha(ziif  ,j,3) * ( 1.0_wp - zcrn**4 )  &
                     + zalpha(ziif  ,j,4) * ( 1.0_wp - zcrn**5 )
          zfi_p(i,j) = MAX( 0.0_wp, zfi_p(i,j) )
        ELSE IF ( zcrn < 0.0_wp ) THEN
          zcrn = 1.0_wp + 2.0_wp * zcrn * sqrtg_r(i+1,j)
          zfi_m(i,j) = zalpha(ziif+1,j,0) * ( 1.0_wp - zcrn )     &
                     - zalpha(ziif+1,j,1) * ( 1.0_wp - zcrn**2 )  &
                     + zalpha(ziif+1,j,2) * ( 1.0_wp - zcrn**3 )  &
                     - zalpha(ziif+1,j,3) * ( 1.0_wp - zcrn**4 )  &
                     + zalpha(ziif+1,j,4) * ( 1.0_wp - zcrn**5 )
          zfi_m(i,j) = MAX( 0.0_wp, zfi_m(i,j) )
        END IF
        zfcr1(i,j) = ( zalpha(ziif,j,0) + zalpha(ziif,j,2) + zalpha(ziif,j,4) )  &
                   * 2.0_wp
      ENDDO
    ENDDO

  END SELECT
  
  DO j = jstart, jend
    DO i = istart-1, iend+1
      ziif = i-icr(i,j)
      ! Compute weighting factor
      zfcr1(i,j) = MAX( zfcr1(i,j), zfi_p(i,j) + zfi_m(i-1,j) + eps_adv ) 
      zfcr1(i,j) = s(ziif,j) / zfcr1(i,j) / sqrtg_r(i,j)
    ENDDO
  ENDDO

  DO j = jstart, jend
    zdx = 1.0_wp / rdx(j)
    DO i = istart-1, iend
      ! Compute and add weighted (fractional) fluxes
      zfx(i,j) = zdx * ( zfx(i,j) +  &
                         zfi_p(i,j)*zfcr1(i,j) - zfi_m(i,j)*zfcr1(i+1,j) )
    ENDDO
  ENDDO
  
  !----------------------------------------------------------------------------
  ! Section 2: Update the scalar field
  !----------------------------------------------------------------------------
  DO j = jstart, jend
    DO i = istart, iend
      s(i,j) = s(i,j) - ( zfx(i,j) - zfx(i-1,j) ) * rdx(j)*sqrtg_r(i,j)
    ENDDO
  ENDDO

END SUBROUTINE xadv_pd_rk_cri_bott

!==============================================================================
!==============================================================================

SUBROUTINE yadv_pd_rk_cri_bott( v_frac, icr, rdy, s, dt, ie, je, &
                                istart, iend, jstart, jend, sqrtg_r )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine updates a scalar prognostic
!   variable 's' due to y-advection. The routine uses the
!   directional split method.
!
! Method:
!   Currently, the Bott advection scheme (Bott 1989) is used.
!   Courant number independent formulation  
!
!------------------------------------------------------------------------------
  
!------------------------------------------------------------------------------
! Declarations:
!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &
  y_scalar_advect  ! type of scalar advection scheme (for qx, TKE and other tracers)

! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER          ,   INTENT (IN) ::  &
  ie, je,         & ! dimensions of the input/output arrays
  istart, iend,   & ! computational start and end indices in i-direction
                    ! (west-east for lat/lon grid)
  jstart, jend      ! computational start and end indices in j-direction
                    ! (south-north for lat/lon grid)

REAL    (KIND=wp),   INTENT(IN)  ::  &
  dt                ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp),   INTENT(IN)  ::  &
  sqrtg_r(ie,je), & ! 1 / square root of G
  rdy(je,2),      & ! reciprocal y-grid spacing (1/( a cos(phi)dphi )
  v_frac(ie,je)     ! fractional v-velocity

INTEGER          ,   INTENT(IN)  ::  &
  icr(ie,je)        ! integer courant numbers

! Array arguments with intent(out):

REAL    (KIND=wp),   INTENT(INOUT) ::  &
  s(ie,je)          ! scalar field updated due to y-advection

!------------------------------------------------------------------------------
! Local variables

INTEGER          ::  &
  i, j, l,           & ! loop indices
  zjif

REAL (KIND=wp)   ::  &
  zdy,               & ! dy
  zsg,               & !
  zcrn                 ! Courant number for intermediate updating

! Local automatic arrays
REAL (KIND=wp)   ::  &
  zalpha(ie,je,0:4), & !
  zfcr1(ie,je),      & !
  zfi_p(ie,je),      & !
  zfi_m(ie,je),      & !  
  zfy(ie,je)        ! Flux in y-direction at cell boundary (v-point)

! End of header
!==============================================================================

  !----------------------------------------------------------------------------
  ! Section 1: Compute monotonic slopes and form fluxes in y - direction
  !----------------------------------------------------------------------------

  zfy (:,:) = 0.0_wp

  SELECT CASE( TRIM(y_scalar_advect) )
    
  CASE( "BOTT2", "BOTT2_STRANG", "BOTT2_XYZYX", "BOTT2_STRANG_B" )
    
    ! Compute polynomials and fluxes for weighting (Cr=1)
    DO j = 2, je-1
      DO i = istart, iend
        zalpha(i,j,0) = abc(2,0,-1) * s(i,j-1)    &
                      + abc(2,0, 0) * s(i,j  )    &
                      + abc(2,0, 1) * s(i,j+1)
        ! unncecessary because of changed Bott coefficients
        ! zalpha(i,j,0) = zalpha(i,j,0) / 2.0_wp

        zalpha(i,j,1) = abc(2,1,-1) * s(i,j-1)    &
                      + abc(2,1, 1) * s(i,j+1)
        ! unncecessary because of changed Bott coefficients
        ! zalpha(i,j,1) = zalpha(i,j,1) / 8.0_wp

        zalpha(i,j,2) = abc(2,2,-1) * s(i,j-1)    &
                      + abc(2,2, 0) * s(i,j  )    &
                      + abc(2,2, 1) * s(i,j+1)
        ! unncecessary because of changed Bott coefficients
        ! zalpha(i,j,2) = zalpha(i,j,2) / 24.0_wp
      ENDDO
    ENDDO
    
    ! Compute courant number and velo-flux at zone edges (v-points)
    DO j = jstart-2, jend+1
      DO i = istart, iend
        ! Compute integer fluxes
        IF ( icr(i,j) >= 1 ) THEN
          zsg = 1.0_wp / sqrtg_r(i,j)
          DO l = 1, icr(i,j) 
            zjif = j-l
            zfy(i,j) = zfy(i,j) + s(i,zjif+1) * zsg
          ENDDO 
        ELSE IF ( icr(i,j) <= -1 ) THEN
          zsg = 1.0_wp / sqrtg_r(i,j+1)
          DO l = -1, icr(i,j), -1
            zjif = j-l
            zfy(i,j) = zfy(i,j) - s(i,zjif  ) * zsg
          ENDDO
        ELSE
          zjif = j
        END IF
        ! Compute integral fluxes
        zfi_p(i,j) = 0.0_wp
        zfi_m(i,j) = 0.0_wp
        zcrn = dt*v_frac(i,j)*rdy(j,2)
        IF ( zcrn > 0.0_wp ) THEN
          zcrn = 1.0_wp - 2.0_wp * zcrn * sqrtg_r(i,j)
          zfi_p(i,j) = zalpha(i,zjif  ,0) * ( 1.0_wp - zcrn )     &
                     + zalpha(i,zjif  ,1) * ( 1.0_wp - zcrn**2 )  &
                     + zalpha(i,zjif  ,2) * ( 1.0_wp - zcrn**3 )
          zfi_p(i,j) = MAX( 0.0_wp, zfi_p(i,j) )
        ELSE IF ( zcrn < 0.0_wp ) THEN
          zcrn = 1.0_wp + 2.0_wp * zcrn * sqrtg_r(i,j+1)
          zfi_m(i,j) = zalpha(i,zjif+1,0) * ( 1.0_wp - zcrn )     &
                     - zalpha(i,zjif+1,1) * ( 1.0_wp - zcrn**2 )  &
                     + zalpha(i,zjif+1,2) * ( 1.0_wp - zcrn**3 )
          zfi_m(i,j) = MAX( 0.0_wp, zfi_m(i,j) )
        END IF
        zfcr1(i,j) = ( zalpha(i,zjif,0) + zalpha(i,zjif,2) )  &
                   * 2.0_wp
      ENDDO
    ENDDO

  CASE( "BOTT4", "BOTT4_STRANG", "BOTT4_XYZYX", "BOTT4_STRANG_B" )

    ! Compute polynomials and fluxes for weighting (Cr=1)
    DO j = 3, je-2
      DO i = istart, iend
        zalpha(i,j,0) = abc(4,0,-2) * s(i,j-2)    &
                      + abc(4,0,-1) * s(i,j-1)    &
                      + abc(4,0, 0) * s(i,j  )    &
                      + abc(4,0, 1) * s(i,j+1)    &
                      + abc(4,0, 2) * s(i,j+2)
        ! unncecessary because of changed Bott coefficients
        ! zalpha(i,j,0) = zalpha(i,j,0) / 2.0_wp

        zalpha(i,j,1) = abc(4,1,-2) * s(i,j-2)    &
                      + abc(4,1,-1) * s(i,j-1)    &
                      + abc(4,1, 1) * s(i,j+1)    &
                      + abc(4,1, 2) * s(i,j+2)
        ! unncecessary because of changed Bott coefficients
        ! zalpha(i,j,1) = zalpha(i,j,1) / 8.0_wp

        zalpha(i,j,2) = abc(4,2,-2) * s(i,j-2)    &
                      + abc(4,2,-1) * s(i,j-1)    &
                      + abc(4,2, 0) * s(i,j  )    &
                      + abc(4,2, 1) * s(i,j+1)    &
                      + abc(4,2, 2) * s(i,j+2)
        ! unncecessary because of changed Bott coefficients
        ! zalpha(i,j,2) = zalpha(i,j,2) / 24.0_wp

        zalpha(i,j,3) = abc(4,3,-2) * s(i,j-2)    &
                      + abc(4,3,-1) * s(i,j-1)    &
                      + abc(4,3, 1) * s(i,j+1)    &
                      + abc(4,3, 2) * s(i,j+2)
        ! unncecessary because of changed Bott coefficients
        ! zalpha(i,j,3) = zalpha(i,j,3) / 64.0_wp

        zalpha(i,j,4) = abc(4,4,-2) * s(i,j-2)    &
                      + abc(4,4,-1) * s(i,j-1)    &
                      + abc(4,4, 0) * s(i,j  )    &
                      + abc(4,4, 1) * s(i,j+1)    &
                      + abc(4,4, 2) * s(i,j+2)
        ! unncecessary because of changed Bott coefficients
        ! zalpha(i,j,4) = zalpha(i,j,4) / 160.0_wp
      ENDDO
    ENDDO
    
    ! Compute courant number and velo-flux at zone edges (v-points)
    DO j = jstart-2, jend+1
      DO i = istart, iend
        ! Compute integer fluxes
        IF ( icr(i,j) >= 1 ) THEN
          zsg = 1.0_wp / sqrtg_r(i,j)
          DO l = 1, icr(i,j) 
            zjif = j-l
            zfy(i,j) = zfy(i,j) + s(i,zjif+1) * zsg
          ENDDO 
        ELSE IF ( icr(i,j) <= -1 ) THEN
          zsg = 1.0_wp / sqrtg_r(i,j+1)
          DO l = -1, icr(i,j), -1
            zjif = j-l
            zfy(i,j) = zfy(i,j) - s(i,zjif  ) * zsg
          ENDDO
        ELSE
          zjif = j
        END IF
        ! Compute integral fluxes
        zfi_p(i,j) = 0.0_wp
        zfi_m(i,j) = 0.0_wp
        zcrn = dt*v_frac(i,j)*rdy(j,2)
        IF ( zcrn > 0.0_wp ) THEN
          zcrn = 1.0_wp - 2.0_wp * zcrn * sqrtg_r(i,j)
          zfi_p(i,j) = zalpha(i,zjif  ,0) * ( 1.0_wp - zcrn )     &
                     + zalpha(i,zjif  ,1) * ( 1.0_wp - zcrn**2 )  &
                     + zalpha(i,zjif  ,2) * ( 1.0_wp - zcrn**3 )  &
                     + zalpha(i,zjif  ,3) * ( 1.0_wp - zcrn**4 )  &
                     + zalpha(i,zjif  ,4) * ( 1.0_wp - zcrn**5 )
          zfi_p(i,j) = MAX( 0.0_wp, zfi_p(i,j) )
        ELSE IF ( zcrn < 0.0_wp ) THEN
          zcrn = 1.0_wp + 2.0_wp * zcrn * sqrtg_r(i,j+1)
          zfi_m(i,j) = zalpha(i,zjif+1,0) * ( 1.0_wp - zcrn )     &
                     - zalpha(i,zjif+1,1) * ( 1.0_wp - zcrn**2 )  &
                     + zalpha(i,zjif+1,2) * ( 1.0_wp - zcrn**3 )  &
                     - zalpha(i,zjif+1,3) * ( 1.0_wp - zcrn**4 )  &
                     + zalpha(i,zjif+1,4) * ( 1.0_wp - zcrn**5 )
          zfi_m(i,j) = MAX( 0.0_wp, zfi_m(i,j) )
        END IF
        zfcr1(i,j) = ( zalpha(i,zjif,0) + zalpha(i,zjif,2) + zalpha(i,zjif,4) )  &
                   * 2.0_wp
      ENDDO
    ENDDO

  END SELECT
      
  DO j = jstart-1, jend+1
    DO i = istart, iend
      zjif = j-icr(i,j)
      ! Compute weighting factor
      zfcr1(i,j) = MAX( zfcr1(i,j), zfi_p(i,j) + zfi_m(i,j-1) + eps_adv )
      zfcr1(i,j) = s(i,zjif) / zfcr1(i,j) / sqrtg_r(i,j)
    ENDDO
  ENDDO

  DO j = jstart-1, jend
    zdy = 1.0_wp / rdy(j,2)
    DO i = istart, iend
      ! Compute and add weighted (fractional) fluxes
      zfy(i,j) = zdy * ( zfy(i,j) +  &
                         zfi_p(i,j)*zfcr1(i,j) - zfi_m(i,j)*zfcr1(i,j+1) )
    ENDDO
  ENDDO
  
  !----------------------------------------------------------------------------
  ! Section 2: Update the scalar field
  !----------------------------------------------------------------------------
  DO j = jstart, jend
    DO i = istart, iend
      s(i,j) = s(i,j) - ( zfy(i,j) - zfy(i,j-1) ) * rdy(j,1)*sqrtg_r(i,j)
    ENDDO
  ENDDO

END SUBROUTINE yadv_pd_rk_cri_bott

!==============================================================================
!============================================================================== 

SUBROUTINE zadv_pd_rk_cri_bott( wc_frac, icr, s, dt, ie, je, ke, ke1, &
                                istart, iend, jstart, jend, sqrtg_r )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine updates a scalar prognostic
!   variable 's' due to vertical advection. The routine uses the
!   directional split method.
!
! Method:
!   Currently, the Bott advection scheme (Bott 1989) is used.
!   Courant number independent formulation  
!
!------------------------------------------------------------------------------
  
!------------------------------------------------------------------------------
! Declarations:
!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &
  y_scalar_advect  ! type of scalar advection scheme (for qx, TKE and other tracers)

! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER          ,   INTENT (IN) ::  &
  ie, je,            & ! x- and y-dimension of the input/output arrays
  ke, ke1,           & ! z-dimension of the input/output arrays
  istart, iend,      & ! computational start and end indices in i-direction
                       ! (west-east for lat/lon grid)
  jstart, jend         ! computational start and end indices in j-direction
                       ! (south-north for lat/lon grid)

REAL    (KIND=wp),   INTENT(IN)  ::  &
  dt                   ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp),   INTENT(IN)  ::  &
  sqrtg_r(ie,je,ke), & ! 1 / square root of G
  wc_frac(ie,je,ke1)   ! fractional contravariant vertical velocity
                       ! (at half levels)

INTEGER          ,   INTENT(IN)  ::  &
  icr(ie,je,ke1)       ! integer courant numbers (at half levels)

! Array arguments with intent(inout):

REAL    (KIND=wp),   INTENT(INOUT) ::  &
  s(ie,je,ke)          ! scalar field updated due to z-advection

!------------------------------------------------------------------------------
! Local variables

INTEGER              ::  &
  i, j, k, l,            & ! loop indices
  zkif

INTEGER              ::  &   !NEC_CB Indexgeneration to vectorize
  icrmax, ii, jj,        & !
  ic1(ke1), ic2(ke1),    & !
  index1(ie*je,ke1,2),   & !
  index2(ie*je,ke1,2),   & !
  zkif_v(ie,je,ke1)

REAL (KIND=wp)       ::  &
  zsg,                   & !
  zcrn                     ! Courant number for intermediate updating

! Local automatic arrays
REAL (KIND=wp)       ::  &
  zalpha(ie,je,ke1,0:4), & !
  zfcr1(ie,je,ke+2),     & !
  zfcr2a(ie*je,ke+2),    & ! NEC_CB Backup for zfcr1
  zfcr2b(ie*je,ke+2),    & ! NEC_CB Backup for zfcr1
  zfi_p(ie,je,ke1),      & !
  zfi_m(ie,je,ke1),      & !  
  zfz(ie,je,ke1)           ! Flux in zeta-direction at cell boundary (w-point)

!NEC_CB For statement function
REAL (KIND=wp)       ::  &
  x1, x2, x3, x4, x5, fzfcr1

! End of header
!==============================================================================

! Statement functions
! -------------------
fzfcr1(x1,x2,x3,x4,x5)=x5 / (MAX( x1, x2 + x3 + eps_adv ) * x4 )
!==============================================================================

  ! check y_scalar_advect only once
  IF (ief_adv_trcr == 0) THEN
    SELECT CASE( TRIM(y_scalar_advect) )
    CASE( "BOTT2", "BOTT2_STRANG", "BOTT2_XYZYX", "BOTT2_STRANG_B" )
      ief_adv_trcr = 2
    CASE( "BOTT4", "BOTT4_STRANG", "BOTT4_XYZYX", "BOTT4_STRANG_B" )
      ief_adv_trcr = 4
    END SELECT
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 1: Compute monotonic slopes and form fluxes in zeta - direction
  !----------------------------------------------------------------------------

  zfz (:,:,:)     = 0.0_wp
  zfcr1(:,:,1)    = 0.0_wp
  zfcr1(:,:,ke+2) = 0.0_wp

  ! Compute polynomials and fluxes for weighting (Cr=1)
  DO j = jstart, jend
    DO i = istart, iend
      ! k = 1
      zalpha(i,j,1,0) = 0.0_wp
      zalpha(i,j,1,1) = 0.0_wp
      zalpha(i,j,1,2) = 0.0_wp
      zalpha(i,j,1,3) = 0.0_wp
      zalpha(i,j,1,4) = 0.0_wp

      ! k = 2
      zalpha(i,j,2,0) = s(i,j,1) / 2.0_wp
      zalpha(i,j,2,1) = 0.0_wp
      zalpha(i,j,2,2) = 0.0_wp
      zalpha(i,j,2,3) = 0.0_wp
      zalpha(i,j,2,4) = 0.0_wp

      ! k = 3
      zalpha(i,j,3,0) = abc(2,0,-1) * s(i,j,1)        &
                      + abc(2,0, 0) * s(i,j,2)        &
                      + abc(2,0, 1) * s(i,j,3)

      zalpha(i,j,3,1) = abc(2,1,-1) * s(i,j,1)        &
                      + abc(2,1, 1) * s(i,j,3)

      zalpha(i,j,3,2) = abc(2,2,-1) * s(i,j,1)        &
                      + abc(2,2, 0) * s(i,j,2)        &
                      + abc(2,2, 1) * s(i,j,3)

      zalpha(i,j,3,3) = 0.0_wp
      zalpha(i,j,3,4) = 0.0_wp

      ! k = ke
      zalpha(i,j,ke,0) = abc(2,0,-1) * s(i,j,ke-2)    &
                       + abc(2,0, 0) * s(i,j,ke-1)    &
                       + abc(2,0, 1) * s(i,j,ke  )

      zalpha(i,j,ke,1) = abc(2,1,-1) * s(i,j,ke-2)    &
                       + abc(2,1, 1) * s(i,j,ke  )

      zalpha(i,j,ke,2) = abc(2,2,-1) * s(i,j,ke-2)    &
                       + abc(2,2, 0) * s(i,j,ke-1)    &
                       + abc(2,2, 1) * s(i,j,ke  )

      zalpha(i,j,ke,3) = 0.0_wp
      zalpha(i,j,ke,4) = 0.0_wp

      ! k = ke1
      zalpha(i,j,ke1,0) = s(i,j,ke) / 2.0_wp
      zalpha(i,j,ke1,1) = 0.0_wp
      zalpha(i,j,ke1,2) = 0.0_wp
      zalpha(i,j,ke1,3) = 0.0_wp
      zalpha(i,j,ke1,4) = 0.0_wp
    ENDDO
  ENDDO
  
  SELECT CASE( ief_adv_trcr )
    
  CASE( 2 )
    
    ! Compute polynomials and fluxes for weighting (Cr=1)
    DO k = 4, ke-1
      DO j = jstart, jend
        DO i = istart, iend
          zalpha(i,j,k,0) = abc(2,0,-1) * s(i,j,k-2)    &
                          + abc(2,0, 0) * s(i,j,k-1)    &
                          + abc(2,0, 1) * s(i,j,k  )

          zalpha(i,j,k,1) = abc(2,1,-1) * s(i,j,k-2)    &
                          + abc(2,1, 1) * s(i,j,k  )

          zalpha(i,j,k,2) = abc(2,2,-1) * s(i,j,k-2)    &
                          + abc(2,2, 0) * s(i,j,k-1)    &
                          + abc(2,2, 1) * s(i,j,k  )
        ENDDO
      ENDDO
    ENDDO

    DO k = 1, ke1
      icrmax=maxval(abs(icr(istart:iend,jstart:jend,k)))
      ii=0
      jj=0
      IF ( icrmax>0 ) THEN
        DO j = jstart, jend
          DO i = istart, iend
            IF ( icr(i,j,k) >= 1 ) THEN
              ii=ii+1
              index1(ii,k,1)=i
              index1(ii,k,2)=j
            ELSE IF ( icr(i,j,k) <= -1 ) THEN
              jj=jj+1
              index2(jj,k,1)=i
              index2(jj,k,2)=j
            END IF
          END DO
        END DO
      END IF
      ic1(k)=ii
      ic2(k)=jj

      ! Compute courant number and velo-flux at zone edges (w-points)

      ! Compute integer fluxes
      DO j = jstart, jend
        DO i = istart, iend
            zkif_v(i,j,k) = k
        ENDDO
      ENDDO

      DO j = jstart, jend
        DO i = istart, iend
          ! Compute integral fluxes
          zfi_p(i,j,k) = 0.0_wp
          zfi_m(i,j,k) = 0.0_wp
          zcrn = dt*wc_frac(i,j,k)
          IF ( zcrn > 0.0_wp ) THEN
            IF ( k >= 2 ) THEN
              zcrn = zcrn * sqrtg_r(i,j,k-1)
              zcrn = 1.0_wp - 2.0_wp * zcrn
              zfi_p(i,j,k) = zalpha(i,j,k  ,0)  &
                           * ( 1.0_wp - zcrn )      &
                           + zalpha(i,j,k  ,1)  &
                           * ( 1.0_wp - zcrn**2 )   &
                           + zalpha(i,j,k  ,2)  &
                           * ( 1.0_wp - zcrn**3 )
              zfi_p(i,j,k) = MAX( 0.0_wp, zfi_p(i,j,k) )
            END IF
          ELSE IF ( zcrn < 0.0_wp ) THEN
            IF ( k+1 <= ke1 ) THEN
              zcrn = zcrn * sqrtg_r(i,j,k)
              zcrn = 1.0_wp + 2.0_wp * zcrn
              zfi_m(i,j,k) = zalpha(i,j,k+1,0)  &
                           * ( 1.0_wp - zcrn )      &
                           - zalpha(i,j,k+1,1)  &
                           * ( 1.0_wp - zcrn**2 )   &
                           + zalpha(i,j,k+1,2)  &
                           * ( 1.0_wp - zcrn**3 )
              zfi_m(i,j,k) = MAX( 0.0_wp, zfi_m(i,j,k) )
            END IF
          END IF
          zfcr1(i,j,k) = ( zalpha(i,j,k,0)                  &
                         + zalpha(i,j,k,2) )  * 2.0_wp
        ENDDO
      ENDDO

      DO ii=1,ic1(k)
        i=index1(ii,k,1)
        j=index1(ii,k,2)

        zsg = 1.0_wp / sqrtg_r(i,j,k-1)
        DO l = 1, icr(i,j,k)
          zkif_v(i,j,k) = k-l
          zfz(i,j,k) = zfz(i,j,k) + s(i,j,zkif_v(i,j,k)  ) * zsg
        ENDDO
          ! Compute integral fluxes
          zfi_p(i,j,k) = 0.0_wp
          zfi_m(i,j,k) = 0.0_wp
          zcrn = dt*wc_frac(i,j,k)
            IF ( zkif_v(i,j,k) >= 2 ) THEN
              zcrn = zcrn * sqrtg_r(i,j,k-1)
              zcrn = 1.0_wp - 2.0_wp * zcrn
              zfi_p(i,j,k) = zalpha(i,j,zkif_v(i,j,k)  ,0)  &
                           * ( 1.0_wp - zcrn )      &
                           + zalpha(i,j,zkif_v(i,j,k)  ,1)  &
                           * ( 1.0_wp - zcrn**2 )   &
                           + zalpha(i,j,zkif_v(i,j,k)  ,2)  &
                           * ( 1.0_wp - zcrn**3 )
              zfi_p(i,j,k) = MAX( 0.0_wp, zfi_p(i,j,k) )
            END IF
          zfcr1(i,j,k) = ( zalpha(i,j,zkif_v(i,j,k),0)                  &
                         + zalpha(i,j,zkif_v(i,j,k),2) )  * 2.0_wp
          zfcr2a(ii,k)=zfcr1(i,j,k)
      END DO

      DO ii=1,ic2(k)
        i=index2(ii,k,1)
        j=index2(ii,k,2)

        zsg = 1.0_wp / sqrtg_r(i,j,k)
        DO l = -1, icr(i,j,k), -1
          zkif_v(i,j,k) = k-l
          zfz(i,j,k) = zfz(i,j,k) - s(i,j,zkif_v(i,j,k)-1) * zsg
        ENDDO

        ! Compute integral fluxes
        zfi_p(i,j,k) = 0.0_wp
        zfi_m(i,j,k) = 0.0_wp
        zcrn = dt*wc_frac(i,j,k)
        IF ( zkif_v(i,j,k)+1 <= ke1 ) THEN
          zcrn = zcrn * sqrtg_r(i,j,k)
          zcrn = 1.0_wp + 2.0_wp * zcrn
          zfi_m(i,j,k) = zalpha(i,j,zkif_v(i,j,k)+1,0)  &
                           * ( 1.0_wp - zcrn )      &
                       - zalpha(i,j,zkif_v(i,j,k)+1,1)  &
                           * ( 1.0_wp - zcrn**2 )   &
                       + zalpha(i,j,zkif_v(i,j,k)+1,2)  &
                           * ( 1.0_wp - zcrn**3 )
          zfi_m(i,j,k) = MAX( 0.0_wp, zfi_m(i,j,k) )
        END IF
        zfcr1(i,j,k) = ( zalpha(i,j,zkif_v(i,j,k),0)                  &
                       + zalpha(i,j,zkif_v(i,j,k),2) )  * 2.0_wp
        zfcr2b(ii,k)=zfcr1(i,j,k)
      ENDDO

    ENDDO !k

    ! the following part is only for the optimized version for ief_adv_trcr=2,
    ! not for ief_adv_trcr=4!
    !NEC_CB First do the major portion directly
    DO k = 2, ke1
      DO j = jstart, jend
        DO i = istart, iend
          ! Compute weighting factor
          !NEC_CB        zfcr1(i,j,k) = MAX( zfcr1(i,j,k), zfi_p(i,j,k) + zfi_m(i,j,k-1) + eps_adv )
          !NEC_CB        zfcr1(i,j,k) = s(i,j,k-1) / zfcr1(i,j,k) / sqrtg_r(i,j,k-1)
          zfcr1(i,j,k)=fzfcr1(zfcr1(i,j,k),zfi_p(i,j,k),zfi_m(i,j,k-1),sqrtg_r(i,j,k-1),s(i,j,k-1))
        ENDDO
      ENDDO
    ENDDO

    !NEC_CB Second do the points with icr<>0 (reading from zfcr2!)
    DO k = 2, ke1
!CDIR NODEP,VOVERTAKE,VOB
      DO ii=1,ic1(k)
        i=index1(ii,k,1)
        j=index1(ii,k,2)
        zkif = k-icr(i,j,k)
        ! Compute weighting factor
        !NEC_CB        zfcr1(i,j,k) = MAX( zfcr2a(ii,k), zfi_p(i,j,k) + zfi_m(i,j,k-1) + eps_adv )
        !NEC_CB        zfcr1(i,j,k) = s(i,j,zkif-1) / zfcr1(i,j,k) / sqrtg_r(i,j,k-1)
        zfcr1(i,j,k)=fzfcr1(zfcr2a(ii,k),zfi_p(i,j,k),zfi_m(i,j,k-1),sqrtg_r(i,j,k-1),s(i,j,zkif-1))
      END DO
!CDIR NODEP,VOVERTAKE,VOB
      DO ii=1,ic2(k)
        i=index2(ii,k,1)
        j=index2(ii,k,2)
        zkif = k-icr(i,j,k)
        ! Compute weighting factor
        !NEC_CB        zfcr1(i,j,k) = MAX( zfcr2b(ii,k), zfi_p(i,j,k) + zfi_m(i,j,k-1) + eps_adv )
        !NEC_CB        zfcr1(i,j,k) = s(i,j,zkif-1) / zfcr1(i,j,k) / sqrtg_r(i,j,k-1)
        zfcr1(i,j,k)=fzfcr1(zfcr2b(ii,k),zfi_p(i,j,k),zfi_m(i,j,k-1),sqrtg_r(i,j,k-1),s(i,j,zkif-1))
      END DO
    END DO

  CASE( 4 )
    
    ! Compute polynomials and fluxes for weighting (Cr=1)
    DO k = 4, ke-1
      DO j = jstart, jend
        DO i = istart, iend
          zalpha(i,j,k,0) = abc(4,0,-2) * s(i,j,k-3)    &
                          + abc(4,0,-1) * s(i,j,k-2)    &
                          + abc(4,0, 0) * s(i,j,k-1)    &
                          + abc(4,0, 1) * s(i,j,k  )    &
                          + abc(4,0, 2) * s(i,j,k+1)

          zalpha(i,j,k,1) = abc(4,1,-2) * s(i,j,k-3)    &
                          + abc(4,1,-1) * s(i,j,k-2)    &
                          + abc(4,1, 1) * s(i,j,k  )    &
                          + abc(4,1, 2) * s(i,j,k+1)

          zalpha(i,j,k,2) = abc(4,2,-2) * s(i,j,k-3)    &
                          + abc(4,2,-1) * s(i,j,k-2)    &
                          + abc(4,2, 0) * s(i,j,k-1)    &
                          + abc(4,2, 1) * s(i,j,k  )    &
                          + abc(4,2, 2) * s(i,j,k+1)

          zalpha(i,j,k,3) = abc(4,3,-2) * s(i,j,k-3)    &
                          + abc(4,3,-1) * s(i,j,k-2)    &
                          + abc(4,3, 1) * s(i,j,k  )    &
                          + abc(4,3, 2) * s(i,j,k+1)

          zalpha(i,j,k,4) = abc(4,4,-2) * s(i,j,k-3)    &
                          + abc(4,4,-1) * s(i,j,k-2)    &
                          + abc(4,4, 0) * s(i,j,k-1)    &
                          + abc(4,4, 1) * s(i,j,k  )    &
                          + abc(4,4, 2) * s(i,j,k+1)
        ENDDO
      ENDDO
    ENDDO
    
    ! Compute courant number and velo-flux at zone edges (w-points)
    DO k = 1, ke1
      DO j = jstart, jend

        DO i = istart, iend
          ! Compute integer fluxes
          IF ( icr(i,j,k) >= 1 ) THEN
            zsg = 1.0_wp / sqrtg_r(i,j,k-1)
            DO l = 1, icr(i,j,k)
              zkif_v(i,j,k) = k-l
              zfz(i,j,k)    = zfz(i,j,k) + s(i,j,zkif_v(i,j,k)  ) * zsg
            ENDDO
          ELSE IF ( icr(i,j,k) <= -1 ) THEN
            zsg = 1.0_wp / sqrtg_r(i,j,k)
            DO l = -1, icr(i,j,k), -1
              zkif_v(i,j,k) = k-l
              zfz(i,j,k)    = zfz(i,j,k) - s(i,j,zkif_v(i,j,k)-1) * zsg
            ENDDO
          ELSE
            zkif_v(i,j,k) = k
          END IF
        ENDDO

        DO i = istart, iend

          ! Compute integral fluxes
          zfi_p(i,j,k) = 0.0_wp
          zfi_m(i,j,k) = 0.0_wp
          zcrn = dt*wc_frac(i,j,k)
          IF ( zcrn > 0.0_wp ) THEN
            IF ( zkif_v(i,j,k) >= 2 ) THEN
              zcrn = zcrn * sqrtg_r(i,j,k-1)
              zcrn = 1.0_wp - 2.0_wp * zcrn
              zfi_p(i,j,k) = zalpha(i,j,zkif_v(i,j,k)  ,0)  &
                           * ( 1.0_wp - zcrn )          &
                           + zalpha(i,j,zkif_v(i,j,k)  ,1)  &
                           * ( 1.0_wp - zcrn**2 )       &
                           + zalpha(i,j,zkif_v(i,j,k)  ,2)  &
                           * ( 1.0_wp - zcrn**3 )       &
                           + zalpha(i,j,zkif_v(i,j,k)  ,3)  &
                           * ( 1.0_wp - zcrn**4 )       &
                           + zalpha(i,j,zkif_v(i,j,k)  ,4)  &
                           * ( 1.0_wp - zcrn**5 )
              zfi_p(i,j,k) = MAX( 0.0_wp, zfi_p(i,j,k) )
            END IF
          ELSE IF ( zcrn < 0.0_wp ) THEN
            IF ( zkif_v(i,j,k)+1 <= ke1 ) THEN
              zcrn = zcrn * sqrtg_r(i,j,k)
              zcrn = 1.0_wp + 2.0_wp * zcrn
              zfi_m(i,j,k) = zalpha(i,j,zkif_v(i,j,k)+1,0)  &
                           * ( 1.0_wp - zcrn )          &
                           - zalpha(i,j,zkif_v(i,j,k)+1,1)  &
                           * ( 1.0_wp - zcrn**2 )       &
                           + zalpha(i,j,zkif_v(i,j,k)+1,2)  &
                           * ( 1.0_wp - zcrn**3 )       &
                           - zalpha(i,j,zkif_v(i,j,k)+1,3)  &
                           * ( 1.0_wp - zcrn**4 )       &
                           + zalpha(i,j,zkif_v(i,j,k)+1,4)  &
                           * ( 1.0_wp - zcrn**5 )
              zfi_m(i,j,k) = MAX( 0.0_wp, zfi_m(i,j,k) )
            END IF
          END IF
          zfcr1(i,j,k) = ( zalpha(i,j,zkif_v(i,j,k),0) +    &
                           zalpha(i,j,zkif_v(i,j,k),2) +    &
                           zalpha(i,j,zkif_v(i,j,k),4) ) * 2.0_wp
        ENDDO
      ENDDO
    ENDDO

    DO k = 2, ke1
      DO j = jstart, jend
        DO i = istart, iend
          zkif = k-icr(i,j,k)
          ! Compute weighting factor
          zfcr1(i,j,k)=fzfcr1(zfcr1(i,j,k),zfi_p(i,j,k),zfi_m(i,j,k-1),        &
            sqrtg_r(i,j,k-1),s(i,j,zkif-1))
        ENDDO
      ENDDO
    ENDDO

  END SELECT

  DO k = 2, ke1
    DO j = jstart, jend
      DO i = istart, iend
        ! Compute and add weighted (fractional) fluxes
        zfz(i,j,k) = zfz(i,j,k) +  &
                     zfi_p(i,j,k)*zfcr1(i,j,k) - zfi_m(i,j,k)*zfcr1(i,j,k+1)
      ENDDO
    ENDDO
  ENDDO

  !----------------------------------------------------------------------------
  ! Section 2: Update the scalar field
  !----------------------------------------------------------------------------

  DO k = 1, ke
    DO j = jstart, jend
      DO i = istart, iend
        s(i,j,k) = s(i,j,k) - ( zfz(i,j,k+1) - zfz(i,j,k) ) * sqrtg_r(i,j,k)
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE zadv_pd_rk_cri_bott

!==============================================================================
!==============================================================================

SUBROUTINE xadv_pd_rk_bott( u, rdx, s, dt, ie, je, ke,           &
                            istart, iend, jstart, jend, kstart, kend, sqrtg_r )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine updates a scalar prognostic
!   variable 's' due to x-advection. The routine uses the
!   directional split method.
!
! Method:
!   Currently, the Bott advection scheme (Bott 1989) is used.
!
!------------------------------------------------------------------------------
  
!------------------------------------------------------------------------------
! Declarations:
!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &
  y_scalar_advect  ! type of scalar advection scheme (for qv, qc, qi, qr, qs, qg, TKE)

! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER          ,   INTENT (IN) ::  &
  ie, je, ke,     & ! dimensions of the input/output arrays
  istart, iend,   & ! computational start and end indices in i-direction
                    ! (west-east for lat/lon grid)
  jstart, jend,   & ! computational start and end indices in j-direction
                    ! (south-north for lat/lon grid)
  kstart, kend      !

REAL    (KIND=wp),   INTENT(IN)  ::  &
  dt                ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp),   INTENT(IN)  ::  &
  sqrtg_r(ie,je,kstart:kend), & ! 1 / square root of G
  rdx       (je)            , & ! reciprocal x-grid spacing (1/( a cos(phi)dlam )
  u      (ie,je,kstart:kend)    ! u-velocity

! Array arguments with intent(out):

REAL    (KIND=wp),   INTENT(INOUT) ::  &
  s      (ie,je,kstart:kend)    ! scalar field updated due to x-advection

!------------------------------------------------------------------------------
! Local variables

INTEGER          ::  &
  i, j, k              ! loop indices

REAL (KIND=wp)   ::  &
  zdx,               & ! dx
  zsg,               & !
  zcrn                 ! Courant number for intermediate updating

! Local automatic arrays
REAL (KIND=wp)               ::  &
  zalpha(ie,je,kstart:kend,0:4), & !
  zfcr1 (ie,je,kstart:kend),     & !
  zfi_p (ie,je,kstart:kend),     & !
  zfi_m (ie,je,kstart:kend),     & !  
  zfx   (ie,je,kstart:kend)        ! Flux in x-direction at cell boundary (u-point)

! End of header
!==============================================================================

  ! check y_scalar_advect only once
  IF (ief_adv_trcr == 0) THEN
    SELECT CASE( TRIM(y_scalar_advect) )
    CASE( "BOTT2", "BOTT2_STRANG", "BOTT2_XYZYX", "BOTT2_STRANG_B" )
      ief_adv_trcr = 2
    CASE( "BOTT4", "BOTT4_STRANG", "BOTT4_XYZYX", "BOTT4_STRANG_B" )
      ief_adv_trcr = 4
    END SELECT
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 1: Compute monotonic slopes and form fluxes in x - direction
  !----------------------------------------------------------------------------

  !NEC: zfx (:,:,:) = 0.0_wp

  SELECT CASE( ief_adv_trcr )
    
  CASE( 2 )
    
    ! Compute polynomials and fluxes for weighting (Cr=1)
    DO k = kstart, kend
      DO j = jstart, jend
        DO i = 2, ie-1
          zalpha(i,j,k,0) = abc(2,0,-1) * s(i-1,j,k)    &
                          + abc(2,0, 0) * s(i  ,j,k)    &
                          + abc(2,0, 1) * s(i+1,j,k)

          zalpha(i,j,k,1) = abc(2,1,-1) * s(i-1,j,k)    &
                          + abc(2,1, 1) * s(i+1,j,k)

          zalpha(i,j,k,2) = abc(2,2,-1) * s(i-1,j,k)    &
                          + abc(2,2, 0) * s(i  ,j,k)    &
                          + abc(2,2, 1) * s(i+1,j,k)
        ENDDO
      ENDDO
    ENDDO
    
    ! Compute courant number and velo-flux at zone edges (u-points)
    DO k = kstart, kend
      DO j = jstart, jend
        DO i = istart-2, iend+1
          ! Compute integral fluxes
          zfi_p(i,j,k) = 0.0_wp
          zfi_m(i,j,k) = 0.0_wp
          zcrn = dt*u(i,j,k)*rdx(j)
          IF ( zcrn > 0.0_wp ) THEN
            zcrn = 1.0_wp - 2.0_wp * zcrn * sqrtg_r(i,j,k)
            zfi_p(i,j,k) = zalpha(i  ,j,k,0) * ( 1.0_wp - zcrn )     &
                         + zalpha(i  ,j,k,1) * ( 1.0_wp - zcrn**2 )  &
                         + zalpha(i  ,j,k,2) * ( 1.0_wp - zcrn**3 )
            zfi_p(i,j,k) = MAX( 0.0_wp, zfi_p(i,j,k) )
          ELSE IF ( zcrn < 0.0_wp ) THEN
            zcrn = 1.0_wp + 2.0_wp * zcrn * sqrtg_r(i+1,j,k)
            zfi_m(i,j,k) = zalpha(i+1,j,k,0) * ( 1.0_wp - zcrn )     &
                         - zalpha(i+1,j,k,1) * ( 1.0_wp - zcrn**2 )  &
                         + zalpha(i+1,j,k,2) * ( 1.0_wp - zcrn**3 )
            zfi_m(i,j,k) = MAX( 0.0_wp, zfi_m(i,j,k) )
          END IF
          zfcr1(i,j,k) = ( zalpha(i,j,k,0) + zalpha(i,j,k,2) )  &
                       * 2.0_wp
        ENDDO
      ENDDO
    ENDDO

  CASE( 4 )
    
    ! Compute polynomials and fluxes for weighting (Cr=1)
    DO k = kstart, kend
      DO j = jstart, jend
        DO i = 3, ie-2
          zalpha(i,j,k,0) = abc(4,0,-2) * s(i-2,j,k)    &
                          + abc(4,0,-1) * s(i-1,j,k)    &
                          + abc(4,0, 0) * s(i  ,j,k)    &
                          + abc(4,0, 1) * s(i+1,j,k)    &
                          + abc(4,0, 2) * s(i+2,j,k)

          zalpha(i,j,k,1) = abc(4,1,-2) * s(i-2,j,k)    &
                          + abc(4,1,-1) * s(i-1,j,k)    &
                          + abc(4,1, 1) * s(i+1,j,k)    &
                          + abc(4,1, 2) * s(i+2,j,k)

          zalpha(i,j,k,2) = abc(4,2,-2) * s(i-2,j,k)    &
                          + abc(4,2,-1) * s(i-1,j,k)    &
                          + abc(4,2, 0) * s(i  ,j,k)    &
                          + abc(4,2, 1) * s(i+1,j,k)    &
                          + abc(4,2, 2) * s(i+2,j,k)

          zalpha(i,j,k,3) = abc(4,3,-2) * s(i-2,j,k)    &
                          + abc(4,3,-1) * s(i-1,j,k)    &
                          + abc(4,3, 1) * s(i+1,j,k)    &
                          + abc(4,3, 2) * s(i+2,j,k)

          zalpha(i,j,k,4) = abc(4,4,-2) * s(i-2,j,k)    &
                          + abc(4,4,-1) * s(i-1,j,k)    &
                          + abc(4,4, 0) * s(i  ,j,k)    &
                          + abc(4,4, 1) * s(i+1,j,k)    &
                          + abc(4,4, 2) * s(i+2,j,k)
        ENDDO
      ENDDO
    ENDDO
    
    ! Compute courant number and velo-flux at zone edges (u-points)
    DO k = kstart, kend
      DO j = jstart, jend
        DO i = istart-2, iend+1
          ! Compute integral fluxes
          zfi_p(i,j,k) = 0.0_wp
          zfi_m(i,j,k) = 0.0_wp
          zcrn = dt*u(i,j,k)*rdx(j)
          IF ( zcrn > 0.0_wp ) THEN
            zcrn = 1.0_wp - 2.0_wp * zcrn * sqrtg_r(i,j,k)
            zfi_p(i,j,k) = zalpha(i  ,j,k,0) * ( 1.0_wp - zcrn )     &
                         + zalpha(i  ,j,k,1) * ( 1.0_wp - zcrn**2 )  &
                         + zalpha(i  ,j,k,2) * ( 1.0_wp - zcrn**3 )  &
                         + zalpha(i  ,j,k,3) * ( 1.0_wp - zcrn**4 )  &
                         + zalpha(i  ,j,k,4) * ( 1.0_wp - zcrn**5 )
            zfi_p(i,j,k) = MAX( 0.0_wp, zfi_p(i,j,k) )
          ELSE IF ( zcrn < 0.0_wp ) THEN
            zcrn = 1.0_wp + 2.0_wp * zcrn * sqrtg_r(i+1,j,k)
            zfi_m(i,j,k) = zalpha(i+1,j,k,0) * ( 1.0_wp - zcrn )     &
                         - zalpha(i+1,j,k,1) * ( 1.0_wp - zcrn**2 )  &
                         + zalpha(i+1,j,k,2) * ( 1.0_wp - zcrn**3 )  &
                         - zalpha(i+1,j,k,3) * ( 1.0_wp - zcrn**4 )  &
                         + zalpha(i+1,j,k,4) * ( 1.0_wp - zcrn**5 )
            zfi_m(i,j,k) = MAX( 0.0_wp, zfi_m(i,j,k) )
          END IF
          zfcr1(i,j,k) = ( zalpha(i,j,k,0) + zalpha(i,j,k,2) + zalpha(i,j,k,4) )  &
                       * 2.0_wp
        ENDDO
      ENDDO
    ENDDO

  END SELECT
  
  DO k = kstart, kend
    DO j = jstart, jend
      DO i = istart-1, iend+1
        ! Compute weighting factor
        zfcr1(i,j,k) = MAX( zfcr1(i,j,k), zfi_p(i,j,k) + zfi_m(i-1,j,k) + eps_adv ) 
        zfcr1(i,j,k) = s(i,j,k) / zfcr1(i,j,k) / sqrtg_r(i,j,k)
      ENDDO
    ENDDO
  ENDDO

  DO k = kstart, kend
    DO j = jstart, jend
      zdx = 1.0_wp / rdx(j)
      DO i = istart-1, iend
        ! Compute weighted fluxes
        zfx(i,j,k) = zdx * ( zfi_p(i,j,k)*zfcr1(i,j,k) - zfi_m(i,j,k)*zfcr1(i+1,j,k) )
      ENDDO
    ENDDO
  ENDDO
  
  !----------------------------------------------------------------------------
  ! Section 2: Update the scalar field
  !----------------------------------------------------------------------------
  DO k = kstart, kend
    DO j = jstart, jend
      DO i = istart, iend
        s(i,j,k) = s(i,j,k) - ( zfx(i,j,k) - zfx(i-1,j,k) ) * rdx(j)*sqrtg_r(i,j,k)
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE xadv_pd_rk_bott

!==============================================================================
!==============================================================================

SUBROUTINE yadv_pd_rk_bott( v, rdy, s, dt, ie, je, ke,            &
                            istart, iend, jstart, jend, kstart, kend, sqrtg_r )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine updates a scalar prognostic
!   variable 's' due to y-advection. The routine uses the
!   directional split method.
!
! Method:
!   Currently, the Bott advection scheme (Bott 1989) is used.
!
!------------------------------------------------------------------------------
  
!------------------------------------------------------------------------------
! Declarations:
!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &
  y_scalar_advect  ! type of scalar advection scheme (for qv, qc, qi, qr, qs, qg, TKE)

! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER          ,   INTENT (IN) ::  &
  ie, je, ke,     & ! dimensions of the input/output arrays
  istart, iend,   & ! computational start and end indices in i-direction
                    ! (west-east for lat/lon grid)
  jstart, jend,   & ! computational start and end indices in j-direction
                    ! (south-north for lat/lon grid)
  kstart, kend

REAL    (KIND=wp),   INTENT(IN)  ::  &
  dt                ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp),   INTENT(IN)  ::  &
  sqrtg_r(ie,je,kstart:kend), & ! 1 / square root of G
  rdy       (je,2)          , & ! reciprocal y-grid spacing (1/( a cos(phi)dphi )
  v      (ie,je,kstart:kend)    ! v-velocity

! Array arguments with intent(out):

REAL    (KIND=wp),   INTENT(INOUT) ::  &
  s      (ie,je,kstart:kend)    ! scalar field updated due to y-advection

!------------------------------------------------------------------------------
! Local variables

INTEGER          ::  &
  i, j, k              ! loop indices

REAL (KIND=wp)   ::  &
  zdy,               & ! dy
  zsg,               & !
  zcrn                 ! Courant number for intermediate updating

! Local automatic arrays
REAL (KIND=wp)               ::  &
  zalpha(ie,je,kstart:kend,0:4), & !
  zfcr1 (ie,je,kstart:kend),     & !
  zfi_p (ie,je,kstart:kend),     & !
  zfi_m (ie,je,kstart:kend),     & !  
  zfy   (ie,je,kstart:kend)        ! Flux in y-direction at cell boundary (v-point)

! End of header
!==============================================================================

  ! check y_scalar_advect only once
  IF (ief_adv_trcr == 0) THEN
    SELECT CASE( TRIM(y_scalar_advect) )
    CASE( "BOTT2", "BOTT2_STRANG", "BOTT2_XYZYX", "BOTT2_STRANG_B" )
      ief_adv_trcr = 2
    CASE( "BOTT4", "BOTT4_STRANG", "BOTT4_XYZYX", "BOTT4_STRANG_B" )
      ief_adv_trcr = 4
    END SELECT
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 1: Compute monotonic slopes and form fluxes in y - direction
  !----------------------------------------------------------------------------

  !NEC: zfy (:,:,:) = 0.0_wp

  SELECT CASE( ief_adv_trcr )
    
  CASE( 2 )
    
    ! Compute polynomials and fluxes for weighting (Cr=1)
    DO k = kstart, kend
      DO j = 2, je-1
        DO i = istart, iend
          zalpha(i,j,k,0) = abc(2,0,-1) * s(i,j-1,k)    &
                          + abc(2,0, 0) * s(i,j  ,k)    &
                          + abc(2,0, 1) * s(i,j+1,k)

          zalpha(i,j,k,1) = abc(2,1,-1) * s(i,j-1,k)    &
                          + abc(2,1, 1) * s(i,j+1,k)

          zalpha(i,j,k,2) = abc(2,2,-1) * s(i,j-1,k)    &
                          + abc(2,2, 0) * s(i,j  ,k)    &
                          + abc(2,2, 1) * s(i,j+1,k)
        ENDDO
      ENDDO
    ENDDO
    
    ! Compute courant number and velo-flux at zone edges (v-points)
    DO k = kstart, kend
      DO j = jstart-2, jend+1
        DO i = istart, iend
          ! Compute integral fluxes
          zfi_p(i,j,k) = 0.0_wp
          zfi_m(i,j,k) = 0.0_wp
          zcrn = dt*v(i,j,k)*rdy(j,2)
          IF ( zcrn > 0.0_wp ) THEN
            zcrn = 1.0_wp - 2.0_wp * zcrn * sqrtg_r(i,j,k)
            zfi_p(i,j,k) = zalpha(i,j  ,k,0) * ( 1.0_wp - zcrn )     &
                         + zalpha(i,j  ,k,1) * ( 1.0_wp - zcrn**2 )  &
                         + zalpha(i,j  ,k,2) * ( 1.0_wp - zcrn**3 )
            zfi_p(i,j,k) = MAX( 0.0_wp, zfi_p(i,j,k) )
          ELSE IF ( zcrn < 0.0_wp ) THEN
            zcrn = 1.0_wp + 2.0_wp * zcrn * sqrtg_r(i,j+1,k)
            zfi_m(i,j,k) = zalpha(i,j+1,k,0) * ( 1.0_wp - zcrn )     &
                         - zalpha(i,j+1,k,1) * ( 1.0_wp - zcrn**2 )  &
                         + zalpha(i,j+1,k,2) * ( 1.0_wp - zcrn**3 )
            zfi_m(i,j,k) = MAX( 0.0_wp, zfi_m(i,j,k) )
          END IF
          zfcr1(i,j,k) = ( zalpha(i,j,k,0) + zalpha(i,j,k,2) )  &
                       * 2.0_wp
        ENDDO
      ENDDO
    ENDDO

  CASE( 4 )
    
    ! Compute polynomials and fluxes for weighting (Cr=1)
    DO k = kstart, kend
      DO j = 3, je-2
        DO i = istart, iend
          zalpha(i,j,k,0) = abc(4,0,-2) * s(i,j-2,k)    &
                          + abc(4,0,-1) * s(i,j-1,k)    &
                          + abc(4,0, 0) * s(i,j  ,k)    &
                          + abc(4,0, 1) * s(i,j+1,k)    &
                          + abc(4,0, 2) * s(i,j+2,k)

          zalpha(i,j,k,1) = abc(4,1,-2) * s(i,j-2,k)    &
                          + abc(4,1,-1) * s(i,j-1,k)    &
                          + abc(4,1, 1) * s(i,j+1,k)    &
                          + abc(4,1, 2) * s(i,j+2,k)

          zalpha(i,j,k,2) = abc(4,2,-2) * s(i,j-2,k)    &
                          + abc(4,2,-1) * s(i,j-1,k)    &
                          + abc(4,2, 0) * s(i,j  ,k)    &
                          + abc(4,2, 1) * s(i,j+1,k)    &
                          + abc(4,2, 2) * s(i,j+2,k)

          zalpha(i,j,k,3) = abc(4,3,-2) * s(i,j-2,k)    &
                          + abc(4,3,-1) * s(i,j-1,k)    &
                          + abc(4,3, 1) * s(i,j+1,k)    &
                          + abc(4,3, 2) * s(i,j+2,k)

          zalpha(i,j,k,4) = abc(4,4,-2) * s(i,j-2,k)    &
                          + abc(4,4,-1) * s(i,j-1,k)    &
                          + abc(4,4, 0) * s(i,j  ,k)    &
                          + abc(4,4, 1) * s(i,j+1,k)    &
                          + abc(4,4, 2) * s(i,j+2,k)
        ENDDO
      ENDDO
    ENDDO
    
    ! Compute courant number and velo-flux at zone edges (v-points)
    DO k = kstart, kend
      DO j = jstart-2, jend+1
        DO i = istart, iend
          ! Compute integral fluxes
          zfi_p(i,j,k) = 0.0_wp
          zfi_m(i,j,k) = 0.0_wp
          zcrn = dt*v(i,j,k)*rdy(j,2)
          IF ( zcrn > 0.0_wp ) THEN
            zcrn = 1.0_wp - 2.0_wp * zcrn * sqrtg_r(i,j,k)
            zfi_p(i,j,k) = zalpha(i,j  ,k,0) * ( 1.0_wp - zcrn )     &
                         + zalpha(i,j  ,k,1) * ( 1.0_wp - zcrn**2 )  &
                         + zalpha(i,j  ,k,2) * ( 1.0_wp - zcrn**3 )  &
                         + zalpha(i,j  ,k,3) * ( 1.0_wp - zcrn**4 )  &
                         + zalpha(i,j  ,k,4) * ( 1.0_wp - zcrn**5 )
            zfi_p(i,j,k) = MAX( 0.0_wp, zfi_p(i,j,k) )
          ELSE IF ( zcrn < 0.0_wp ) THEN
            zcrn = 1.0_wp + 2.0_wp * zcrn * sqrtg_r(i,j+1,k)
            zfi_m(i,j,k) = zalpha(i,j+1,k,0) * ( 1.0_wp - zcrn )     &
                         - zalpha(i,j+1,k,1) * ( 1.0_wp - zcrn**2 )  &
                         + zalpha(i,j+1,k,2) * ( 1.0_wp - zcrn**3 )  &
                         - zalpha(i,j+1,k,3) * ( 1.0_wp - zcrn**4 )  &
                         + zalpha(i,j+1,k,4) * ( 1.0_wp - zcrn**5 )
            zfi_m(i,j,k) = MAX( 0.0_wp, zfi_m(i,j,k) )
          END IF
          zfcr1(i,j,k) = ( zalpha(i,j,k,0) + zalpha(i,j,k,2) + zalpha(i,j,k,4) )  &
                       * 2.0_wp
        ENDDO
      ENDDO
    ENDDO

  END SELECT
      
  DO k = kstart, kend
    DO j = jstart-1, jend+1
      DO i = istart, iend
        ! Compute weighting factor
        zfcr1(i,j,k) = MAX( zfcr1(i,j,k), zfi_p(i,j,k) + zfi_m(i,j-1,k) + eps_adv )
        zfcr1(i,j,k) = s(i,j,k) / (zfcr1(i,j,k) * sqrtg_r(i,j,k))
      ENDDO
    ENDDO
  ENDDO

  DO k = kstart, kend
    DO j = jstart-1, jend
      zdy = 1.0_wp / rdy(j,2)
      DO i = istart, iend
        ! Compute weighted fluxes
        zfy(i,j,k) = zdy * ( zfi_p(i,j,k)*zfcr1(i,j,k) - zfi_m(i,j,k)*zfcr1(i,j+1,k) )
      ENDDO
    ENDDO
  ENDDO
  
  !----------------------------------------------------------------------------
  ! Section 2: Update the scalar field
  !----------------------------------------------------------------------------
  DO k = kstart, kend
    DO j = jstart, jend
      DO i = istart, iend
        s(i,j,k) = s(i,j,k) - ( zfy(i,j,k) - zfy(i,j-1,k) ) * rdy(j,1)*sqrtg_r(i,j,k)
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE yadv_pd_rk_bott

!==============================================================================
!============================================================================== 

SUBROUTINE zadv_pd_rk_bott( wc, s, dt, ie, je, ke, ke1,        &
                            istart, iend, jstart, jend, sqrtg_r )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine updates a scalar prognostic
!   variable 's' due to vertical advection. The routine uses the
!   directional split method.
!
! Method:
!   Currently, the Bott advection scheme (Bott 1989) is used.
!
!------------------------------------------------------------------------------
  
!------------------------------------------------------------------------------
! Declarations:
!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &
  y_scalar_advect  ! type of scalar advection scheme (for qv, qc, qi, qr, qs, qg, TKE)

! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER          ,   INTENT (IN) ::  &
  ie, je,            & ! x- and y-dimension of the input/output arrays
  ke, ke1,           & ! z-dimension of the input/output arrays
  istart, iend,      & ! computational start and end indices in i-direction
                       ! (west-east for lat/lon grid)
  jstart, jend         ! computational start and end indices in j-direction
                       ! (south-north for lat/lon grid)

REAL    (KIND=wp),   INTENT(IN)  ::  &
  dt                   ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp),   INTENT(IN)  ::  &
  sqrtg_r(ie,je,ke), & ! 1 / square root of G
  wc(ie,je,ke1)        ! contravariant vertical velocity (at half levels)

! Array arguments with intent(inout):

REAL    (KIND=wp),   INTENT(INOUT) ::  &
  s(ie,je,ke)          ! scalar field updated due to z-advection

!------------------------------------------------------------------------------
! Local variables

INTEGER              ::  &
  i, j, k                  ! loop indices

REAL (KIND=wp)       ::  &
  zsg,                   & !
  zcrn                     ! Courant number for intermediate updating

! Local automatic arrays
REAL (KIND=wp)       ::  &
  zalpha(ie,je,ke1,0:4), & !
  zfcr1(ie,je,ke+2),     & !
  zfi_p(ie,je,ke1),      & !
  zfi_m(ie,je,ke1),      & !  
  zfz(ie,je,ke1)           ! Flux in zeta-direction at cell boundary (w-point)

! End of header
!==============================================================================

  ! check y_scalar_advect only once
  IF (ief_adv_trcr == 0) THEN
    SELECT CASE( TRIM(y_scalar_advect) )
    CASE( "BOTT2", "BOTT2_STRANG", "BOTT2_XYZYX", "BOTT2_STRANG_B" )
      ief_adv_trcr = 2
    CASE( "BOTT4", "BOTT4_STRANG", "BOTT4_XYZYX", "BOTT4_STRANG_B" )
      ief_adv_trcr = 4
    END SELECT
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 1: Compute monotonic slopes and form fluxes in zeta - direction
  !----------------------------------------------------------------------------

  zfz (:,:,:)     = 0.0_wp
  zfcr1(:,:,1)    = 0.0_wp
  zfcr1(:,:,ke+2) = 0.0_wp

  ! Compute polynomials and fluxes for weighting (Cr=1)
  DO j = jstart, jend
    DO i = istart, iend
      ! k = 1
      zalpha(i,j,1,0) = 0.0_wp
      zalpha(i,j,1,1) = 0.0_wp
      zalpha(i,j,1,2) = 0.0_wp
      zalpha(i,j,1,3) = 0.0_wp
      zalpha(i,j,1,4) = 0.0_wp
      ! k = 2
      zalpha(i,j,2,0) = s(i,j,1) / 2.0_wp
      zalpha(i,j,2,1) = 0.0_wp
      zalpha(i,j,2,2) = 0.0_wp
      zalpha(i,j,2,3) = 0.0_wp
      zalpha(i,j,2,4) = 0.0_wp
      ! k = 3
      zalpha(i,j,3,0) = abc(2,0,-1) * s(i,j,1)        &
                      + abc(2,0, 0) * s(i,j,2)        &
                      + abc(2,0, 1) * s(i,j,3)

      zalpha(i,j,3,1) = abc(2,1,-1) * s(i,j,1)        &
                      + abc(2,1, 1) * s(i,j,3)

      zalpha(i,j,3,2) = abc(2,2,-1) * s(i,j,1)        &
                      + abc(2,2, 0) * s(i,j,2)        &
                      + abc(2,2, 1) * s(i,j,3)

      zalpha(i,j,3,3) = 0.0_wp
      zalpha(i,j,3,4) = 0.0_wp
      ! k = ke
      zalpha(i,j,ke,0) = abc(2,0,-1) * s(i,j,ke-2)    &
                       + abc(2,0, 0) * s(i,j,ke-1)    &
                       + abc(2,0, 1) * s(i,j,ke  )

      zalpha(i,j,ke,1) = abc(2,1,-1) * s(i,j,ke-2)    &
                       + abc(2,1, 1) * s(i,j,ke  )

      zalpha(i,j,ke,2) = abc(2,2,-1) * s(i,j,ke-2)    &
                       + abc(2,2, 0) * s(i,j,ke-1)    &
                       + abc(2,2, 1) * s(i,j,ke  )

      zalpha(i,j,ke,3) = 0.0_wp
      zalpha(i,j,ke,4) = 0.0_wp
      ! k = ke1
      zalpha(i,j,ke1,0) = s(i,j,ke) / 2.0_wp
      zalpha(i,j,ke1,1) = 0.0_wp
      zalpha(i,j,ke1,2) = 0.0_wp
      zalpha(i,j,ke1,3) = 0.0_wp
      zalpha(i,j,ke1,4) = 0.0_wp
    ENDDO
  ENDDO
  
  SELECT CASE( ief_adv_trcr )
    
  CASE( 2 )
    
    ! Compute polynomials and fluxes for weighting (Cr=1)
    DO k = 4, ke-1
      DO j = jstart, jend
        DO i = istart, iend
          zalpha(i,j,k,0) = abc(2,0,-1) * s(i,j,k-2)    &
                          + abc(2,0, 0) * s(i,j,k-1)    &
                          + abc(2,0, 1) * s(i,j,k  )

          zalpha(i,j,k,1) = abc(2,1,-1) * s(i,j,k-2)    &
                          + abc(2,1, 1) * s(i,j,k  )

          zalpha(i,j,k,2) = abc(2,2,-1) * s(i,j,k-2)    &
                          + abc(2,2, 0) * s(i,j,k-1)    &
                          + abc(2,2, 1) * s(i,j,k  )
        ENDDO
      ENDDO
    ENDDO
    
    ! Compute courant number and velo-flux at zone edges (w-points)
    DO k = 1, ke1
      DO j = jstart, jend
        DO i = istart, iend
          ! Compute integral fluxes
          zfi_p(i,j,k) = 0.0_wp
          zfi_m(i,j,k) = 0.0_wp
          zcrn = dt*wc(i,j,k)
          IF ( zcrn > 0.0_wp ) THEN
            IF ( k >= 2 ) THEN
              zcrn = zcrn * sqrtg_r(i,j,k-1)
              zcrn = 1.0_wp - 2.0_wp * zcrn
              zfi_p(i,j,k) = zalpha(i,j,k  ,0) * ( 1.0_wp - zcrn )     &
                           + zalpha(i,j,k  ,1) * ( 1.0_wp - zcrn**2 )  &
                           + zalpha(i,j,k  ,2) * ( 1.0_wp - zcrn**3 )
              zfi_p(i,j,k) = MAX( 0.0_wp, zfi_p(i,j,k) )
            END IF
          ELSE IF ( zcrn < 0.0_wp ) THEN
            IF ( k+1 <= ke1 ) THEN
              zcrn = zcrn * sqrtg_r(i,j,k)
              zcrn = 1.0_wp + 2.0_wp * zcrn
              zfi_m(i,j,k) = zalpha(i,j,k+1,0) * ( 1.0_wp - zcrn )     &
                           - zalpha(i,j,k+1,1) * ( 1.0_wp - zcrn**2 )  &
                           + zalpha(i,j,k+1,2) * ( 1.0_wp - zcrn**3 )
              zfi_m(i,j,k) = MAX( 0.0_wp, zfi_m(i,j,k) )
            END IF
          END IF
          zfcr1(i,j,k) = ( zalpha(i,j,k,0) + zalpha(i,j,k,2) )  &
                       * 2.0_wp
        ENDDO
      ENDDO
    ENDDO

  CASE( 4 )
    
    ! Compute polynomials and fluxes for weighting (Cr=1)
    DO k = 4, ke-1
      DO j = jstart, jend
        DO i = istart, iend
          zalpha(i,j,k,0) = abc(4,0,-2) * s(i,j,k-3)    &
                          + abc(4,0,-1) * s(i,j,k-2)    &
                          + abc(4,0, 0) * s(i,j,k-1)    &
                          + abc(4,0, 1) * s(i,j,k  )    &
                          + abc(4,0, 2) * s(i,j,k+1)

          zalpha(i,j,k,1) = abc(4,1,-2) * s(i,j,k-3)    &
                          + abc(4,1,-1) * s(i,j,k-2)    &
                          + abc(4,1, 1) * s(i,j,k  )    &
                          + abc(4,1, 2) * s(i,j,k+1)

          zalpha(i,j,k,2) = abc(4,2,-2) * s(i,j,k-3)    &
                          + abc(4,2,-1) * s(i,j,k-2)    &
                          + abc(4,2, 0) * s(i,j,k-1)    &
                          + abc(4,2, 1) * s(i,j,k  )    &
                          + abc(4,2, 2) * s(i,j,k+1)

          zalpha(i,j,k,3) = abc(4,3,-2) * s(i,j,k-3)    &
                          + abc(4,3,-1) * s(i,j,k-2)    &
                          + abc(4,3, 1) * s(i,j,k  )    &
                          + abc(4,3, 2) * s(i,j,k+1)

          zalpha(i,j,k,4) = abc(4,4,-2) * s(i,j,k-3)    &
                          + abc(4,4,-1) * s(i,j,k-2)    &
                          + abc(4,4, 0) * s(i,j,k-1)    &
                          + abc(4,4, 1) * s(i,j,k  )    &
                          + abc(4,4, 2) * s(i,j,k+1)
        ENDDO
      ENDDO
    ENDDO
    
    ! Compute courant number and velo-flux at zone edges (w-points)
    DO k = 1, ke1
      DO j = jstart, jend
        DO i = istart, iend
          ! Compute integral fluxes
          zfi_p(i,j,k) = 0.0_wp
          zfi_m(i,j,k) = 0.0_wp
          zcrn = dt*wc(i,j,k)
          IF ( zcrn > 0.0_wp ) THEN
            IF ( k >= 2 ) THEN
              zcrn = zcrn * sqrtg_r(i,j,k-1)
              zcrn = 1.0_wp - 2.0_wp * zcrn
              zfi_p(i,j,k) = zalpha(i,j,k  ,0) * ( 1.0_wp - zcrn )     &
                           + zalpha(i,j,k  ,1) * ( 1.0_wp - zcrn**2 )  &
                           + zalpha(i,j,k  ,2) * ( 1.0_wp - zcrn**3 )  &
                           + zalpha(i,j,k  ,3) * ( 1.0_wp - zcrn**4 )  &
                           + zalpha(i,j,k  ,4) * ( 1.0_wp - zcrn**5 )
              zfi_p(i,j,k) = MAX( 0.0_wp, zfi_p(i,j,k) )
            END IF
          ELSE IF ( zcrn < 0.0_wp ) THEN
            IF ( k+1 <= ke1 ) THEN
              zcrn = zcrn * sqrtg_r(i,j,k)
              zcrn = 1.0_wp + 2.0_wp * zcrn
              zfi_m(i,j,k) = zalpha(i,j,k+1,0) * ( 1.0_wp - zcrn )     &
                           - zalpha(i,j,k+1,1) * ( 1.0_wp - zcrn**2 )  &
                           + zalpha(i,j,k+1,2) * ( 1.0_wp - zcrn**3 )  &
                           - zalpha(i,j,k+1,3) * ( 1.0_wp - zcrn**4 )  &
                           + zalpha(i,j,k+1,4) * ( 1.0_wp - zcrn**5 )
              zfi_m(i,j,k) = MAX( 0.0_wp, zfi_m(i,j,k) )
            END IF
          END IF
          zfcr1(i,j,k) = ( zalpha(i,j,k,0) + zalpha(i,j,k,2) + zalpha(i,j,k,4) )  &
                       * 2.0_wp
        ENDDO
      ENDDO
    ENDDO

  END SELECT
      
  DO k = 2, ke1
    DO j = jstart, jend
      DO i = istart, iend
        ! Compute weighting factor
        zfcr1(i,j,k) = MAX( zfcr1(i,j,k), zfi_p(i,j,k) + zfi_m(i,j,k-1) + eps_adv )
        zfcr1(i,j,k) = s(i,j,k-1) / zfcr1(i,j,k) / sqrtg_r(i,j,k-1)
      ENDDO
    ENDDO
  ENDDO

  DO k = 2, ke1
    DO j = jstart, jend
      DO i = istart, iend
        ! Compute weighted fluxes
        zfz(i,j,k) = zfi_p(i,j,k)*zfcr1(i,j,k) - zfi_m(i,j,k)*zfcr1(i,j,k+1)
      ENDDO
    ENDDO
  ENDDO
  
  !----------------------------------------------------------------------------
  ! Section 2: Update the scalar field
  !----------------------------------------------------------------------------
  DO k = 1, ke
    DO j = jstart, jend
      DO i = istart, iend
        s(i,j,k) = s(i,j,k) - ( zfz(i,j,k+1) - zfz(i,j,k) ) * sqrtg_r(i,j,k)
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE zadv_pd_rk_bott

!==============================================================================
!============================================================================== 

SUBROUTINE ufrac_crint_rk( u, u_frac, icr, rdx, dt, ie, je, ke,       &
                           istart, iend, jstart, jend, kstart, kend,  &
                           intcr_max, sqrtg_r, lintcr_ne_zero )

!------------------------------------------------------------------------------
!
! Description:
!
!------------------------------------------------------------------------------

! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER          ,   INTENT (IN) ::  &
  intcr_max,      & ! max. allowed integer courant number
  ie, je, ke,     & ! dimensions of the input/output arrays
  istart, iend,   & ! computational start and end indices in i-direction
                    ! (west-east for lat/lon grid)
  jstart, jend,   & ! computational start and end indices in j-direction
                    ! (south-north for lat/lon grid)
  kstart, kend      !

REAL    (KIND=wp),   INTENT(IN)  ::  &
  dt                ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp),   INTENT(IN)  ::  &
  sqrtg_r(ie,je,kstart:kend), & ! 1 / square root of G
  rdx       (je),    & ! reciprocal x-grid spacing (1/( a cos(phi)dlam )
  u      (ie,je,kstart:kend)    ! u-velocity

! Scalar arguments with intent(out):

LOGICAL,             INTENT (OUT) ::  &
  lintcr_ne_zero

! Array arguments with intent(out):

REAL    (KIND=wp),   INTENT(OUT)  ::  &
  u_frac(ie,je,kstart:kend)     ! fractional u-velocity

INTEGER          ,   INTENT(OUT)  ::  &
  icr(ie,je,kstart:kend)        ! integer courant numbers

!------------------------------------------------------------------------------
! Local variables

INTEGER        :: &
  i, j, k           ! loop indices

REAL (KIND=wp) :: &
  zu_int,         & !
  zdt_o_dx,       & ! dt / dx
  zdx_o_dt          ! dx / dt

! End of header
!==============================================================================

  lintcr_ne_zero = .FALSE.
  
  ! Compute integer Courant numbers and (fractional) transport velocites
DO k= kstart, kend
  DO j = jstart, jend
    zdt_o_dx = dt*rdx(j)
    zdx_o_dt = 1.0_wp / zdt_o_dx
    DO i = istart-2, iend+1
      IF ( u(i,j,k) > 0.0_wp ) THEN
        icr(i,j,k) = INT( u(i,j,k)*zdt_o_dx*sqrtg_r(i,j,k) )
      ELSE IF ( u(i,j,k) < 0.0_wp ) THEN
        icr(i,j,k) = INT( u(i,j,k)*zdt_o_dx*sqrtg_r(i+1,j,k) )
      ELSE
        icr(i,j,k) = 0
      END IF
      IF ( icr(i,j,k) > 0) THEN
        icr(i,j,k) = MIN( icr(i,j,k), intcr_max )
        zu_int   = REAL( icr(i,j,k), wp) * zdx_o_dt/sqrtg_r(i,j,k)
        u_frac(i,j,k) = MOD( u(i,j,k), zu_int )
        lintcr_ne_zero = .TRUE.
      ELSE IF ( icr(i,j,k) < 0) THEN
        icr(i,j,k) = MAX( icr(i,j,k), -intcr_max )
        zu_int   = REAL( icr(i,j,k), wp) * zdx_o_dt/sqrtg_r(i+1,j,k)
        u_frac(i,j,k) = MOD( u(i,j,k), zu_int )
        lintcr_ne_zero = .TRUE.
      ELSE
        u_frac(i,j,k) = u(i,j,k)
      END IF
    ENDDO
  ENDDO
ENDDO

  IF ( intcr_max == 0) lintcr_ne_zero = .FALSE.

END SUBROUTINE ufrac_crint_rk

!==============================================================================
!============================================================================== 

SUBROUTINE vfrac_crint_rk( v, v_frac, icr, rdy, dt, ie, je, ke,       &
                           istart, iend, jstart, jend, kstart, kend,  &
                           intcr_max, sqrtg_r, lintcr_ne_zero )

!------------------------------------------------------------------------------
!
! Description:
!
!------------------------------------------------------------------------------

! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER          ,   INTENT (IN) ::  &
  intcr_max,      & ! max. allowed integer courant number
  ie, je, ke,     & ! dimensions of the input/output arrays
  istart, iend,   & ! computational start and end indices in i-direction
                    ! (west-east for lat/lon grid)
  jstart, jend,   & ! computational start and end indices in j-direction
                    ! (south-north for lat/lon grid)
  kstart, kend

REAL    (KIND=wp),   INTENT(IN)  ::  &
  dt                ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp),   INTENT(IN)  ::  &
  sqrtg_r(ie,je,kstart:kend), & ! 1 / square root of G
  rdy(je,2),                  & ! reciprocal y-grid spacing (1/( a cos(phi)dphi )
  v      (ie,je,kstart:kend)    ! v-velocity

! Scalar arguments with intent(out):

LOGICAL,             INTENT (OUT) ::  &
  lintcr_ne_zero

! Array arguments with intent(out):

REAL    (KIND=wp),   INTENT(OUT)  ::  &
  v_frac(ie,je,kstart:kend)     ! fractional v-velocity

INTEGER          ,   INTENT(OUT)  ::  &
  icr(ie,je,kstart:kend)        ! integer courant numbers

!------------------------------------------------------------------------------
! Local variables

INTEGER        :: &
  i, j, k           ! loop indices

REAL (KIND=wp) :: &
  zv_int,         & !
  zdt_o_dy,       & ! dt / dy
  zdy_o_dt          ! dy / dt

! End of header
!==============================================================================

  lintcr_ne_zero = .FALSE.
  
  ! Compute integer Courant numbers and (fractional) transport velocites
DO k = kstart, kend
  DO j = jstart-2, jend+1
    zdt_o_dy = dt*rdy(j,2)
    zdy_o_dt = 1.0_wp / zdt_o_dy
    DO i = istart, iend
      IF ( v(i,j,k) > 0.0_wp ) THEN
        icr(i,j,k) = INT( v(i,j,k)*zdt_o_dy*sqrtg_r(i,j,k) )
      ELSE IF ( v(i,j,k) < 0.0_wp ) THEN
        icr(i,j,k) = INT( v(i,j,k)*zdt_o_dy*sqrtg_r(i,j+1,k) )
      ELSE
        icr(i,j,k) = 0
      END IF
      IF ( icr(i,j,k) > 0) THEN
        icr(i,j,k) = MIN( icr(i,j,k), intcr_max )
        zv_int   = REAL( icr(i,j,k), wp) * zdy_o_dt/sqrtg_r(i,j,k)
        v_frac(i,j,k) = MOD( v(i,j,k), zv_int )
        lintcr_ne_zero = .TRUE.
      ELSE IF ( icr(i,j,k) < 0) THEN
        icr(i,j,k) = MAX( icr(i,j,k), -intcr_max )
        zv_int   = REAL( icr(i,j,k), wp) * zdy_o_dt/sqrtg_r(i,j+1,k)
        v_frac(i,j,k) = MOD( v(i,j,k), zv_int )
        lintcr_ne_zero = .TRUE.
      ELSE
        v_frac(i,j,k) = v(i,j,k)
      END IF
    ENDDO
  ENDDO
ENDDO

  IF ( intcr_max == 0) lintcr_ne_zero = .FALSE.

END SUBROUTINE vfrac_crint_rk

!==============================================================================
!============================================================================== 

!option! -pvctl _on_adb
SUBROUTINE wcfrac_crint_rk( wc, wc_frac, icr, dt, ie, je, ke, ke1,        &
                            istart, iend, jstart, jend, sqrtg_r,          &
                            lintcr_ne_zero, ivl_off_opt, num_compute,     &
                            icomm_cart, imp_integers)

!------------------------------------------------------------------------------
!
! Description:
!
!------------------------------------------------------------------------------

! Subroutine arguments

! Scalar arguments with intent(in):

INTEGER          ,   INTENT (IN) ::  &
  ie, je,            & ! x- and y-dimension of the input/output arrays
  ke, ke1,           & ! z-dimension of the input/output arrays
  istart, iend,      & ! computational start and end indices in i-direction
                       ! (west-east for lat/lon grid)
  jstart, jend,      & ! computational start and end indices in j-direction
                       ! (south-north for lat/lon grid)
  num_compute,       & ! number of compute PEs
  icomm_cart,        & ! communicator for the virtual cartesian topology
  imp_integers         ! determines the correct INTEGER type used in the

REAL    (KIND=wp),   INTENT(IN)  ::  &
  dt                   ! timestep

! Array arguments with intent(in):

REAL    (KIND=wp),   INTENT(IN)  ::  &
  sqrtg_r(ie,je,ke), & ! 1 / square root of G
  wc(ie,je,ke1)        ! fractional contravariant vertical velocity

! Scalar arguments with intent(out):

LOGICAL,             INTENT (OUT) ::  &
  lintcr_ne_zero

! Array arguments with intent(out):

REAL    (KIND=wp),   INTENT(OUT)  ::  &
  wc_frac(ie,je,ke1)   ! fractional contravariant vertical velocity

INTEGER          ,   INTENT(OUT)  ::  &
  icr(ie,je,ke1)       ! integer courant numbers

INTEGER          ,   INTENT(IN), OPTIONAL ::  &
  ivl_off_opt          ! offset differs between PPM and van Leer scheme
  
!------------------------------------------------------------------------------
! Local variables

INTEGER          ::  &
  i, j, k,           & ! loop indices
  ivl_off,           & ! offset differs between PPM and van Leer scheme
  kmin, kmax, izmax, izerror

REAL (KIND=wp)  ::   &
  zwc_int, zwc_tmp,  & !
  zrdt, zdz_o_dt       ! 1.0 / dt

CHARACTER (LEN=80) ::  yzerrmsg

! End of header
!==============================================================================

  izerror  = 0
  yzerrmsg = ''

  IF ( PRESENT(ivl_off_opt) ) THEN
    ivl_off = ivl_off_opt
  ELSE
    ivl_off = 1
  END IF
  
  zrdt = 1.0_wp / dt

  ! Compute integer Courant numbers and limit
  ! fractional contravariant vertical transport velocity
  ! to Cr=1.0 at upper and lower boundary
  DO k = 1, ke1
    kmin = MAX( k-1, 1 )
    kmax = MIN( k, ke )
    DO j = jstart, jend
      DO i = istart, iend
        IF ( wc(i,j,k) > 0.0_wp ) THEN
          icr(i,j,k) = INT( wc(i,j,k)*dt*sqrtg_r(i,j,kmin) )
        ELSE IF ( wc(i,j,k) < 0.0_wp ) THEN
          icr(i,j,k) = INT( wc(i,j,k)*dt*sqrtg_r(i,j,kmax) )
        ELSE
          icr(i,j,k) = 0
        END IF
        IF ( icr(i,j,k) > 0) THEN
          zdz_o_dt = zrdt/sqrtg_r(i,j,kmin)
          icr(i,j,k) = MAX( 0, MIN( k-2, icr(i,j,k) ) )
          IF ( icr(i,j,k) > 0) THEN
            zwc_int    = REAL( icr(i,j,k), wp) * zdz_o_dt
            zwc_tmp    = MOD( wc(i,j,k), zwc_int )
            wc_frac(i,j,k) = MIN( zwc_tmp,  zdz_o_dt )
          ELSE
            wc_frac(i,j,k) = zdz_o_dt
          ENDIF
        ELSE IF ( icr(i,j,k) < 0) THEN
          zdz_o_dt = zrdt/sqrtg_r(i,j,kmax)
          icr(i,j,k) = MIN( 0, MAX( k-ke+ivl_off, icr(i,j,k) ) )
          IF ( icr(i,j,k) < 0) THEN
            zwc_int    = REAL( icr(i,j,k), wp) * zdz_o_dt
            zwc_tmp    = MOD( wc(i,j,k), zwc_int )
            wc_frac(i,j,k) = MAX( zwc_tmp, -zdz_o_dt )
          ELSE
            wc_frac(i,j,k) = -zdz_o_dt
          ENDIF
        ELSE
          wc_frac(i,j,k) = wc(i,j,k)
        END IF
      ENDDO
    ENDDO
  ENDDO

  izmax = MAXVAL ( ABS(icr(istart:iend,jstart:jend,2:ke1) ) )
!NEC_CB Why? Doing cri_bott makes no difference if icr(i,j,k)=0.
  IF (num_compute > 1) THEN
    CALL global_values( izmax, 1, 'MAX', imp_integers, icomm_cart, -1, yzerrmsg, izerror )
  ENDIF

  IF (izmax == 0) THEN
    lintcr_ne_zero = .FALSE.
  ELSE
    lintcr_ne_zero = .TRUE.
  END IF

END SUBROUTINE wcfrac_crint_rk

!==============================================================================

FUNCTION udsdx( icase, iorder, s,             &
                ie, je, ke, i, j, k, im, ip,  &
                veloi_x, rdxy_x,              &
                veloi_y, rdxy_y, ssign )

!------------------------------------------------------------------------------

REAL   (KIND=wp)     :: udsdx

INTEGER                , INTENT(IN) :: icase  ! 0: x-advection + y-advection 
                                              ! 1: x-advection
                                              ! 2: y-advection
                                              ! 3: z-advection
INTEGER         , INTENT(IN) :: iorder
INTEGER         , INTENT(IN) :: ie, je, ke
REAL   (KIND=wp), INTENT(IN) :: s(ie,je,ke)
INTEGER         , INTENT(IN) :: im, ip
INTEGER         , INTENT(IN) :: i, j, k
REAL   (KIND=wp), INTENT(IN) :: veloi_x, rdxy_x
REAL   (KIND=wp), INTENT(IN) :: veloi_y, rdxy_y
REAL   (KIND=wp), INTENT(IN) :: ssign

! Local variables
!----------------

REAL   (KIND=wp)     :: s0, sim, sip, sm1, sm2, sm3, sp1, sp2, sp3
REAL   (KIND=wp)     :: sim_x, sip_x, sm1_x, sm2_x, sm3_x, sp1_x, sp2_x, sp3_x
REAL   (KIND=wp)     :: sim_y, sip_y, sm1_y, sm2_y, sm3_y, sp1_y, sp2_y, sp3_y

! End of header
!==============================================================================


  IF ( icase == 0 ) THEN
    !
    ! x-advection + y-advection
    !
    ! replace SELECT by IF...ELSEIF  (seems to be better on the NEC)
    IF     (iorder==1) THEN
      s0    = s(i,   j,   k)
      sim_x = s(i+im,j,   k)
      sip_x = s(i+ip,j,   k)
      sim_y = s(i,   j+im,k)
      sip_y = s(i,   j+ip,k)
      
      ! 1st order velo*ds/dx operator
      udsdx       = rdxy_x * 0.5_wp * (                              &
                veloi_x * ( - sim_x + sip_x )                            &
         + ssign * ABS(veloi_x) * ( - sim_x + 2._wp*s0 - sip_x ) )   &
                 +  rdxy_y * 0.5_wp * (                              &
                veloi_y * ( - sim_y + sip_y )                            &
         + ssign * ABS(veloi_y) * ( - sim_y + 2._wp*s0 - sip_y ) )
      
    ELSEIF (iorder==2) THEN
      sim_x = s(i+im,j,   k)
      sip_x = s(i+ip,j,   k)
      sim_y = s(i,   j+im,k)
      sip_y = s(i,   j+ip,k)
      
      ! 2nd order velo*ds/dx operator
      udsdx       = rdxy_x * 0.5_wp * veloi_x * ( - sim_x + sip_x )  &
                  + rdxy_y * 0.5_wp * veloi_y * ( - sim_y + sip_y )
      
    ELSEIF (iorder==3) THEN
      s0    = s(i,   j,   k)
      sm1_x = s(i-1, j,   k)
      sm2_x = s(i-2, j,   k)
      sp1_x = s(i+1, j,   k)
      sp2_x = s(i+2, j,   k)
      sm1_y = s(i,   j-1, k)
      sm2_y = s(i,   j-2, k)
      sp1_y = s(i,   j+1, k)
      sp2_y = s(i,   j+2, k)
      
      ! 3rd order velo*ds/dx operator
      udsdx       = rdxy_x * 1._wp/12._wp * (                    &
                veloi_x * ( sm2_x - 8._wp*(sm1_x-sp1_x) - sp2_x )    &
         + ssign * ABS(veloi_x) * ( sm2_x - 4._wp*(sm1_x+sp1_x)      &
                                  + 6._wp*s0 + sp2_x) )              &
                  + rdxy_y * 1._wp/12._wp * (                    &
                veloi_y * ( sm2_y - 8._wp*(sm1_y-sp1_y) - sp2_y )    &
         + ssign * ABS(veloi_y) * ( sm2_y - 4._wp*(sm1_y+sp1_y)      &
                                  + 6._wp*s0 + sp2_y) )
      
    ELSEIF (iorder==4) THEN
      sm1_x = s(i-1, j,   k)
      sm2_x = s(i-2, j,   k)
      sp1_x = s(i+1, j,   k)
      sp2_x = s(i+2, j,   k)
      sm1_y = s(i,   j-1, k)
      sm2_y = s(i,   j-2, k)
      sp1_y = s(i,   j+1, k)
      sp2_y = s(i,   j+2, k)
      
      ! 4th order velo*ds/dx operator
      udsdx       = rdxy_x * 1._wp/12._wp *                      &
                veloi_x * ( sm2_x - 8._wp*(sm1_x-sp1_x) - sp2_x )    &
                  + rdxy_y * 1._wp/12._wp *                      &
                veloi_y * ( sm2_y - 8._wp*(sm1_y-sp1_y) - sp2_y )
      
    ELSEIF (iorder==5) THEN
      s0    = s(i,   j,   k)
      sm1_x = s(i-1, j,   k)
      sm2_x = s(i-2, j,   k)
      sm3_x = s(i-3, j,   k)
      sp1_x = s(i+1, j,   k)
      sp2_x = s(i+2, j,   k)
      sp3_x = s(i+3, j,   k)
      sm1_y = s(i,   j-1, k)
      sm2_y = s(i,   j-2, k)
      sm3_y = s(i,   j-3, k)
      sp1_y = s(i,   j+1, k)
      sp2_y = s(i,   j+2, k)
      sp3_y = s(i,   j+3, k)
      
      ! 5th order velo*ds/dx operator
      udsdx       = rdxy_x * 1._wp/60._wp * (                    &
                veloi_x * ( - sm3_x + 9._wp*(sm2_x-sp2_x)            &
                                   - 45._wp*(sm1_x-sp1_x) + sp3_x )  &
         + ssign * ABS(veloi_x) * ( - sm3_x + 6._wp*(sm2_x+sp2_x)    &
                                   - 15._wp*(sm1_x+sp1_x)            &
                                   + 20._wp*s0 - sp3_x ) )           &
                  + rdxy_y * 1._wp/60._wp * (                    &
                veloi_y * ( - sm3_y + 9._wp*(sm2_y-sp2_y)            &
                                   - 45._wp*(sm1_y-sp1_y) + sp3_y )  &
         + ssign * ABS(veloi_y) * ( - sm3_y + 6._wp*(sm2_y+sp2_y)    &
                                   - 15._wp*(sm1_y+sp1_y)            &
                                   + 20._wp*s0 - sp3_y ) )
      
    ELSEIF (iorder==6) THEN
      sm1_x = s(i-1, j,   k)
      sm2_x = s(i-2, j,   k)
      sm3_x = s(i-3, j,   k)
      sp1_x = s(i+1, j,   k)
      sp2_x = s(i+2, j,   k)
      sp3_x = s(i+3, j,   k)
      sm1_y = s(i,   j-1, k)
      sm2_y = s(i,   j-2, k)
      sm3_y = s(i,   j-3, k)
      sp1_y = s(i,   j+1, k)
      sp2_y = s(i,   j+2, k)
      sp3_y = s(i,   j+3, k)
      
      ! 6th order velo*ds/dx operator
      udsdx       = rdxy_x * 1._wp/60._wp *                      &
                veloi_x * ( - sm3_x + 9._wp*(sm2_x-sp2_x)            &
                                   - 45._wp*(sm1_x-sp1_x) + sp3_x )  &
                  + rdxy_y * 1._wp/60._wp *                      &
                veloi_y * ( - sm3_y + 9._wp*(sm2_y-sp2_y)            &
                                   - 45._wp*(sm1_y-sp1_y) + sp3_y )
      
    ELSE
      udsdx=0.0_wp !NEC_CB To make the compiler see a result in any case and avoid him assuming dependencies
    ENDIF

  ELSE
    
    SELECT CASE(iorder)

    CASE(1)
      SELECT CASE(icase)
      CASE(1)
        ! x-advection
        s0  = s(i,   j,k)
        sim = s(i+im,j,k)
        sip = s(i+ip,j,k)
      CASE(2)
        ! y-advection
        s0  = s(i,j   ,k)
        sim = s(i,j+im,k)
        sip = s(i,j+ip,k)
      CASE(3)
        ! z-advection
        s0  = s(i,j,k   )
        sim = s(i,j,k+im)
        sip = s(i,j,k+ip)
      END SELECT

      ! 1st order velo*ds/dx operator
      udsdx       = rdxy_x * 0.5_wp * (                              &
                veloi_x * ( - sim + sip )                                &
         + ssign * ABS(veloi_x) * ( - sim + 2._wp*s0 - sip ) )
      
    CASE(2)
      SELECT CASE(icase)
      CASE(1)
        ! x-advection
        sim = s(i+im,j,k)
        sip = s(i+ip,j,k)
      CASE(2)
        ! y-advection
        sim = s(i,j+im,k)
        sip = s(i,j+ip,k)
      CASE(3)
        ! z-advection
        sim = s(i,j,k+im)
        sip = s(i,j,k+ip)
      END SELECT

      ! 2nd order velo*ds/dx operator
      udsdx       = rdxy_x * 0.5_wp * veloi_x * ( - sim + sip )
      
    CASE(3)
      SELECT CASE(icase)
      CASE(1)
        ! x-advection
        s0  = s(i,   j,k)
        sm1 = s(i-1, j,k)
        sm2 = s(i-1, j,k)
        sp1 = s(i+1, j,k)
        sp2 = s(i+2, j,k)
      CASE(2)
        ! y-advection
        s0  = s(i,j   ,k)
        sm1 = s(i,j-1, k)
        sm2 = s(i,j-2, k)
        sp1 = s(i,j+1, k)
        sp2 = s(i,j+2, k)
      CASE(3)
        ! z-advection
        s0  = s(i,j,k   )
        sm1 = s(i,j,k-1 )
        sm2 = s(i,j,k-2 )
        sp1 = s(i,j,k+1 )
        sp2 = s(i,j,k+2 )
      END SELECT

      ! 3rd order velo*ds/dx operator
      udsdx       = rdxy_x * 1._wp/12._wp * (                    &
                veloi_x * ( sm2 - 8._wp*(sm1-sp1) - sp2 )            &
         + ssign * ABS(veloi_x) * ( sm2 - 4._wp*(sm1+sp1)            &
                                + 6._wp*s0 + sp2) )
      
    CASE(4)
      SELECT CASE(icase)
      CASE(1)
        ! x-advection
        sm1 = s(i-1, j,k)
        sm2 = s(i-1, j,k)
        sp1 = s(i+1, j,k)
        sp2 = s(i+2, j,k)
      CASE(2)
        ! y-advection
        sm1 = s(i,j-1, k)
        sm2 = s(i,j-2, k)
        sp1 = s(i,j+1, k)
        sp2 = s(i,j+2, k)
      CASE(3)
        ! z-advection
        sm1 = s(i,j,k-1 )
        sm2 = s(i,j,k-2 )
        sp1 = s(i,j,k+1 )
        sp2 = s(i,j,k+2 )
      END SELECT

      ! 4th order velo*ds/dx operator
      udsdx       = rdxy_x * 1._wp/12._wp *                      &
                veloi_x * ( sm2 - 8._wp*(sm1-sp1) - sp2 )
      
    CASE(5)
      SELECT CASE(icase)
      CASE(1)
        ! x-advection
        s0  = s(i,   j,k)
        sm1 = s(i-1, j,k)
        sm2 = s(i-1, j,k)
        sm3 = s(i-3, j,k)
        sp1 = s(i+1, j,k)
        sp2 = s(i+2, j,k)
        sp3 = s(i+3, j,k)
      CASE(2)
        ! y-advection
        s0  = s(i,j   ,k)
        sm1 = s(i,j-1, k)
        sm2 = s(i,j-2, k)
        sm3 = s(i,j-3, k)
        sp1 = s(i,j+1, k)
        sp2 = s(i,j+2, k)
        sp3 = s(i,j+3, k)
      CASE(3)
        ! z-advection
        s0  = s(i,j,k   )
        sm1 = s(i,j,k-1 )
        sm2 = s(i,j,k-2 )
        sm3 = s(i,j,k-3 )
        sp1 = s(i,j,k+1 )
        sp2 = s(i,j,k+2 )
        sp3 = s(i,j,k+3 )
      END SELECT

      ! 5th order velo*ds/dx operator
      udsdx       = rdxy_x * 1._wp/60._wp * (                    &
                veloi_x * ( - sm3 + 9._wp*(sm2-sp2)                  &
                                 - 45._wp*(sm1-sp1) + sp3 )          &
         + ssign * ABS(veloi_x) * ( - sm3 + 6._wp*(sm2+sp2)          &
                                 - 15._wp*(sm1+sp1)                  &
                                 + 20._wp*s0 - sp3 ) )
      
    CASE(6)
      SELECT CASE(icase)
      CASE(1)
        ! x-advection
        sm1 = s(i-1, j,k)
        sm2 = s(i-1, j,k)
        sm3 = s(i-3, j,k)
        sp1 = s(i+1, j,k)
        sp2 = s(i+2, j,k)
        sp3 = s(i+3, j,k)
      CASE(2)
        ! y-advection
        sm1 = s(i,j-1, k)
        sm2 = s(i,j-2, k)
        sm3 = s(i,j-3, k)
        sp1 = s(i,j+1, k)
        sp2 = s(i,j+2, k)
        sp3 = s(i,j+3, k)
      CASE(3)
        ! z-advection
        sm1 = s(i,j,k-1 )
        sm2 = s(i,j,k-2 )
        sm3 = s(i,j,k-3 )
        sp1 = s(i,j,k+1 )
        sp2 = s(i,j,k+2 )
        sp3 = s(i,j,k+3 )
      END SELECT

      ! 6th order velo*ds/dx operator
      udsdx       = rdxy_x * 1._wp/60._wp *                      &
                veloi_x * ( - sm3 + 9._wp*(sm2-sp2)                  &
                                 - 45._wp*(sm1-sp1) + sp3 )
      
    END SELECT

  END IF

END FUNCTION udsdx

!==============================================================================

!==============================================================================

FUNCTION udsdx_up5_xy( s,                            &
                       ie, je, ke, i, j, k,          &
                       veloi_x, rdxy_x,              &
                       veloi_y, rdxy_y, ssign )

!------------------------------------------------------------------------------

REAL   (KIND=wp)     :: udsdx_up5_xy

INTEGER         , INTENT(IN) :: ie, je, ke
REAL   (KIND=wp), INTENT(IN) :: s(ie,je,ke)
INTEGER         , INTENT(IN) :: i, j, k
REAL   (KIND=wp), INTENT(IN) :: veloi_x, rdxy_x
REAL   (KIND=wp), INTENT(IN) :: veloi_y, rdxy_y
REAL   (KIND=wp), INTENT(IN) :: ssign

! Local variables
!----------------

REAL   (KIND=wp)     :: s0
REAL   (KIND=wp)     :: sm1_x, sm2_x, sm3_x, sp1_x, sp2_x, sp3_x
REAL   (KIND=wp)     :: sm1_y, sm2_y, sm3_y, sp1_y, sp2_y, sp3_y

! End of header
!==============================================================================

  ! x-advection + y-advection
  s0    = s(i,   j,   k)
  sm1_x = s(i-1, j,   k)
  sm2_x = s(i-2, j,   k)
  sm3_x = s(i-3, j,   k)
  sp1_x = s(i+1, j,   k)
  sp2_x = s(i+2, j,   k)
  sp3_x = s(i+3, j,   k)
  sm1_y = s(i,   j-1, k)
  sm2_y = s(i,   j-2, k)
  sm3_y = s(i,   j-3, k)
  sp1_y = s(i,   j+1, k)
  sp2_y = s(i,   j+2, k)
  sp3_y = s(i,   j+3, k)
  
  ! 5th order velo*ds/dx operator
  udsdx_up5_xy       = rdxy_x * 1._wp/60._wp * (             &
            veloi_x * ( - sm3_x + 9._wp*(sm2_x-sp2_x)            &
                               - 45._wp*(sm1_x-sp1_x) + sp3_x )  &
     + ssign * ABS(veloi_x) * ( - sm3_x + 6._wp*(sm2_x+sp2_x)    &
                               - 15._wp*(sm1_x+sp1_x)            &
                               + 20._wp*s0 - sp3_x ) )           &
                     + rdxy_y * 1._wp/60._wp * (             &
            veloi_y * ( - sm3_y + 9._wp*(sm2_y-sp2_y)            &
                               - 45._wp*(sm1_y-sp1_y) + sp3_y )  &
     + ssign * ABS(veloi_y) * ( - sm3_y + 6._wp*(sm2_y+sp2_y)    &
                               - 15._wp*(sm1_y+sp1_y)            &
                               + 20._wp*s0 - sp3_y ) )
  
END FUNCTION udsdx_up5_xy

!==============================================================================

!==============================================================================

FUNCTION udsdx_up3_z( s,                            &
                      ie, je, ke, i, j, k,          &
                      veloi, rdz, ssign )

!------------------------------------------------------------------------------

REAL   (KIND=wp)     :: udsdx_up3_z

INTEGER         , INTENT(IN) :: ie, je, ke
REAL   (KIND=wp), INTENT(IN) :: s(ie,je,ke)
INTEGER         , INTENT(IN) :: i, j, k
REAL   (KIND=wp), INTENT(IN) :: veloi, rdz, ssign

! Local variables
!----------------

REAL   (KIND=wp)     :: s0, sm1, sm2, sp1, sp2

! End of header
!==============================================================================


  ! z-advection
  s0  = s(i,j,k   )
  sm1 = s(i,j,k-1 )
  sm2 = s(i,j,k-2 )
  sp1 = s(i,j,k+1 )
  sp2 = s(i,j,k+2 )
  
  ! 3rd order velo*ds/dx operator
  udsdx_up3_z       = rdz * 1._wp/12._wp * (              &
            veloi * ( sm2 - 8._wp*(sm1-sp1) - sp2 )           &
     + ssign * ABS(veloi) * ( sm2 - 4._wp*(sm1+sp1)           &
                            + 6._wp*s0 + sp2) )

END FUNCTION udsdx_up3_z

!==============================================================================

!==============================================================================

FUNCTION udsdx_up1_z( s,                            &
                      ie, je, ke, i, j, k,          &
                      im, ip,                       &
                      veloi, rdz, ssign )

!------------------------------------------------------------------------------

REAL   (KIND=wp)     :: udsdx_up1_z

INTEGER         , INTENT(IN) :: ie, je, ke
REAL   (KIND=wp), INTENT(IN) :: s(ie,je,ke)
INTEGER         , INTENT(IN) :: im, ip
INTEGER         , INTENT(IN) :: i, j, k
REAL   (KIND=wp), INTENT(IN) :: veloi, rdz, ssign

! Local variables
!----------------

REAL   (KIND=wp)     :: s0, sim, sip

! End of header
!==============================================================================


  ! z-advection
  s0  = s(i,j,k   )
  sim = s(i,j,k+im)
  sip = s(i,j,k+ip)

  ! 1st order velo*ds/dx operator
  udsdx_up1_z       = rdz * 0.5_wp * (                        &
            veloi * ( - sim + sip )                               &
     + ssign * ABS(veloi) * ( - sim + 2._wp*s0 - sip ) )
  
END FUNCTION udsdx_up1_z

!==============================================================================

!==============================================================================

FUNCTION udsdx_up1(s,im,ip,veloi,rdxy, ssign)

!------------------------------------------------------------------------------

REAL   (KIND=wp), PARAMETER  :: c1 = 0.5_wp

REAL   (KIND=wp)     :: udsdx_up1

INTEGER         , INTENT(IN) ::   im, ip
REAL   (KIND=wp), INTENT(IN) ::   s(im:ip)
REAL   (KIND=wp), INTENT(IN) ::   veloi, rdxy, ssign

!------------------------------------------------------------------------------

  ! 1st order velo*ds/dx operator
  udsdx_up1 =  rdxy * c1 * (                                       &
              veloi * ( - s(im) + s(ip) )                          &
       + ssign * ABS(veloi) * ( - s(im) + 2._wp*s(0) - s(ip) ) &
       )

END FUNCTION udsdx_up1

!==============================================================================
!==============================================================================

FUNCTION udsdx_cd2(s,im,ip,veloi,rdxy)

!------------------------------------------------------------------------------

REAL   (KIND=wp), PARAMETER  :: c1 = 0.5_wp

REAL   (KIND=wp)     :: udsdx_cd2

INTEGER         , INTENT(IN) ::  im, ip
REAL   (KIND=wp), INTENT(IN) ::  s(im:ip)
REAL   (KIND=wp), INTENT(IN) ::  veloi, rdxy

!------------------------------------------------------------------------------

  ! 2nd order velo*ds/dx operator
  udsdx_cd2 =  rdxy * c1 * veloi * ( - s(im) + s(ip) )

END FUNCTION udsdx_cd2

!==============================================================================
!==============================================================================

FUNCTION udsdx_up3(s,im,ip,veloi,rdxy, ssign)

!------------------------------------------------------------------------------

REAL   (KIND=wp), PARAMETER  :: c1 = 1._wp/12._wp

REAL   (KIND=wp)     :: udsdx_up3

INTEGER         , INTENT(IN) ::   im, ip
REAL   (KIND=wp), INTENT(IN) ::   s(im:ip)
REAL   (KIND=wp), INTENT(IN) ::   veloi, rdxy, ssign

!------------------------------------------------------------------------------

  ! 3rd order velo*ds/dx operator
  udsdx_up3 =  rdxy * c1 * (                                              &
                     veloi  * ( s(-2) - 8._wp*(s(-1)-s(1)) - s(2) )   &
       + ssign * ABS(veloi) *                                             &
             ( s(-2) - 4._wp*(s(-1)+s(1)) + 6._wp*s(0) + s(2) )   &
       )

END FUNCTION udsdx_up3

!==============================================================================
!==============================================================================

FUNCTION udsdx_cd4(s,im,ip,veloi,rdxy)

!------------------------------------------------------------------------------

REAL   (KIND=wp), PARAMETER  :: c1 = 1._wp/12._wp

REAL   (KIND=wp)     :: udsdx_cd4

INTEGER         , INTENT(IN) ::  im, ip
REAL   (KIND=wp), INTENT(IN) ::  s(im:ip)
REAL   (KIND=wp), INTENT(IN) ::  veloi, rdxy

!------------------------------------------------------------------------------

  ! 4th order velo*ds/dx operator
  udsdx_cd4 =  rdxy * c1 * veloi * ( s(-2) - 8._wp*(s(-1)-s(1)) - s(2) )

END FUNCTION udsdx_cd4

!==============================================================================
!==============================================================================

FUNCTION udsdx_up5(s,im,ip,veloi,rdxy, ssign)

!------------------------------------------------------------------------------

REAL   (KIND=wp), PARAMETER  :: c1 = 1._wp/60._wp

REAL   (KIND=wp)     :: udsdx_up5

INTEGER         , INTENT(IN) ::  im, ip
REAL   (KIND=wp), INTENT(IN) ::  s(im:ip)
REAL   (KIND=wp), INTENT(IN) ::  veloi, rdxy, ssign

!------------------------------------------------------------------------------

    ! 5th order velo*ds/dx operator
    udsdx_up5 =  rdxy * c1 * (                                              &
                   veloi * ( - s(-3) + 9._wp*(s(-2)-s(2))               &
                             - 45._wp*(s(-1)-s(1)) + s(3) )             &
         +    ssign * ABS(veloi) * ( - s(-3) + 6._wp*(s(-2)+s(2))       &
                  - 15._wp*(s(-1)+s(1)) + 20._wp*s(0) - s(3))       &
         )

END FUNCTION udsdx_up5

!==============================================================================
!==============================================================================

FUNCTION udsdx_cd6(s,im,ip,veloi,rdxy)

!------------------------------------------------------------------------------

REAL   (KIND=wp), PARAMETER  :: c1 = 1._wp/60._wp

REAL   (KIND=wp)     :: udsdx_cd6

INTEGER         , INTENT(IN) ::  im, ip
REAL   (KIND=wp), INTENT(IN) ::  s(im:ip)
REAL   (KIND=wp), INTENT(IN) ::  veloi, rdxy

!------------------------------------------------------------------------------

  ! 6th order velo*ds/dx operator
  udsdx_cd6 =  rdxy * c1 * veloi * ( - s(-3) + 9._wp*(s(-2)-s(2)) &
                                            - 45._wp*(s(-1)-s(1)) + s(3) )

END FUNCTION udsdx_cd6

!==============================================================================
!==============================================================================

FUNCTION dfdx_up1(s,im,ip,veloim1,veloi,rdxy)

!------------------------------------------------------------------------------

REAL   (KIND=wp), PARAMETER  :: c1 = 0.5_wp

REAL   (KIND=wp)     :: dfdx_up1
INTEGER         , INTENT(IN) ::  im, ip
REAL   (KIND=wp), INTENT(IN) ::  s(im:ip)
REAL   (KIND=wp), INTENT(IN) ::  veloi, veloim1, rdxy

!------------------------------------------------------------------------------

  ! 1st order df/dx operator
  dfdx_up1 =  rdxy * ( &
       c1 * (    veloi * ( s(0) + s(1) )    &
       +    ABS(veloi) * ( s(0) - s(1) )    &
       +       veloim1 * ( - s(-1) - s(0) ) &
       +  ABS(veloim1) * ( - s(-1) + s(0) ) &
       ) &
! JF:         - s(0) * ( veloi - veloim1 )  &
         )

END FUNCTION dfdx_up1

!==============================================================================
!==============================================================================

FUNCTION dfdx_cd2(s,im,ip,veloim1,veloi,rdxy)

!------------------------------------------------------------------------------

REAL   (KIND=wp), PARAMETER  :: c1 = 0.5_wp

REAL   (KIND=wp)     :: dfdx_cd2

INTEGER         , INTENT(IN) ::  im, ip
REAL   (KIND=wp), INTENT(IN) ::  s(im:ip)
REAL   (KIND=wp), INTENT(IN) ::  veloi, veloim1, rdxy

!------------------------------------------------------------------------------

  ! 2nd order df/dx operator
  dfdx_cd2 =  rdxy * ( &
       c1 * (    veloi * ( s(0) + s(1) )    &
       +       veloim1 * ( - s(-1) - s(0) ) &
       ) &
! JF:          - s(0) * ( veloi - veloim1 )  &
       )

END FUNCTION dfdx_cd2

!==============================================================================
!==============================================================================

FUNCTION dfdx_up3(s,im,ip,veloim1,veloi,rdxy)

!------------------------------------------------------------------------------

REAL(KIND=wp), PARAMETER  :: c1 = 1._wp/12._wp

REAL(KIND=wp)     :: dfdx_up3

INTEGER      , INTENT(IN) ::  im, ip
REAL(KIND=wp), INTENT(IN) ::  s(im:ip)
REAL(KIND=wp), INTENT(IN) ::  veloi, veloim1, rdxy

!------------------------------------------------------------------------------

  ! 3rd order df/dx operator
  dfdx_up3 =  rdxy * ( &
       c1 * (    veloi * ( - s(-1) + 7._wp*(s(0)+s(1)) - s(2) )  &
       +    ABS(veloi) * ( - s(-1) + 3._wp*(s(0)-s(1)) + s(2) )  &
       +       veloim1 * (   s(-2) - 7._wp*(s(-1)+s(0)) + s(1) ) &
       +  ABS(veloim1) * (   s(-2) - 3._wp*(s(-1)-s(0)) - s(1) ) &
       ) &
! JF:          - s(0) * ( veloi - veloim1 )  &
       )

END FUNCTION dfdx_up3

!==============================================================================
!==============================================================================

FUNCTION dfdx_cd4(s,im,ip,veloim1,veloi,rdxy)

!------------------------------------------------------------------------------

REAL   (KIND=wp), PARAMETER  :: c1 = 1._wp/12._wp

REAL   (KIND=wp)     :: dfdx_cd4

INTEGER         , INTENT(IN) ::  im, ip
REAL   (KIND=wp), INTENT(IN) ::  s(im:ip)
REAL   (KIND=wp), INTENT(IN) ::  veloi, veloim1, rdxy

!------------------------------------------------------------------------------

  ! 4th order df/dx operator
  dfdx_cd4 =  rdxy * ( &
       c1 * (    veloi * ( - s(-1) + 7._wp*(s(0)+s(1))  - s(2) ) &
       +       veloim1 * (   s(-2) - 7._wp*(s(-1)+s(0)) + s(1) ) &
       ) &
! JF:          - s(0) * ( veloi - veloim1 )  &
       )

END FUNCTION dfdx_cd4

!==============================================================================
!==============================================================================

FUNCTION dfdx_up5(s,im,ip,veloim1,veloi,rdxy)

!------------------------------------------------------------------------------

REAL   (KIND=wp), PARAMETER  :: c1 = 1._wp/60._wp

REAL   (KIND=wp)     :: dfdx_up5

INTEGER         , INTENT(IN) ::  im, ip
REAL   (KIND=wp), INTENT(IN) ::  s(im:ip)
REAL   (KIND=wp), INTENT(IN) ::  veloi, veloim1, rdxy

!------------------------------------------------------------------------------

  ! 5th order df/dx operator
  dfdx_up5 =  rdxy * ( &
       c1 * (    veloi * (   s(-2) + 37._wp*(s(0)+s(1))  -           &
                                         8._wp*(s(-1)+s(2)) + s(3) ) &
       +    ABS(veloi) * (   s(-2) + 10._wp*(s(0)-s(1))  -           &
                                         5._wp*(s(-1)-s(2)) - s(3) ) &
       +       veloim1 * ( - s(-3) - 37._wp*(s(-1)+s(0)) +           &
                                         8._wp*(s(-2)+s(1)) - s(2) ) &
       +  ABS(veloim1) * ( - s(-3) - 10._wp*(s(-1)-s(0)) +           &
                                         5._wp*(s(-2)-s(1)) + s(2) ) &
       ) &
! JF:          - s(0) * ( veloi - veloim1 )  &
       )

END FUNCTION dfdx_up5

!==============================================================================
!==============================================================================

FUNCTION dfdx_cd6(s,im,ip,veloim1,veloi,rdxy)

!------------------------------------------------------------------------------

REAL   (KIND=wp), PARAMETER  :: c1 = 1._wp/60._wp

REAL   (KIND=wp)     :: dfdx_cd6
INTEGER         , INTENT(IN) ::  im, ip
REAL   (KIND=wp), INTENT(IN) ::  s(im:ip)
REAL   (KIND=wp), INTENT(IN) ::  veloi, veloim1, rdxy

!------------------------------------------------------------------------------

  ! 6th order df/dx operator
  dfdx_cd6 =  rdxy * ( &
    c1 * (    veloi * (   s(-2) + 37._wp*(s(0)+s(1))  -           &
                                      8._wp*(s(-1)+s(2)) + s(3) ) &
    +       veloim1 * ( - s(-3) - 37._wp*(s(-1)+s(0)) +           &
                                      8._wp*(s(-2)+s(1)) - s(2) ) &
    ) &
! JF:          - s(0) * ( veloi - veloim1 )  &
    )

END FUNCTION dfdx_cd6

!==============================================================================
!==============================================================================
  
FUNCTION flux_weno3(s,im,ip,veloi)

!------------------------------------------------------------------------------

REAL   (KIND=wp)     :: flux_weno3

INTEGER         , INTENT(IN) ::  im, ip
REAL   (KIND=wp), INTENT(IN) ::  s(im:ip)
REAL   (KIND=wp), INTENT(IN) ::  veloi

REAL   (KIND=wp)             ::  phi(5)

!------------------------------------------------------------------------------

  IF (veloi >= 0.0_wp) THEN
    phi(1) = s(-2)
    phi(2) = s(-1)
    phi(3) = s(0)
    phi(4) = s(1)
    phi(5) = s(2)
  ELSE
    phi(1) = s(3)
    phi(2) = s(2)
    phi(3) = s(1)
    phi(4) = s(0)
    phi(5) = s(-1)
  END IF
    
  flux_weno3 = flux_weno3_priv( veloi,phi )

END FUNCTION flux_weno3

!==============================================================================
!==============================================================================

FUNCTION flux_weno3_priv( c_adv, phi )

!------------------------------------------------------------------------------

REAL(KIND=wp)     :: flux_weno3_priv

REAL   (KIND=wp),     PARAMETER :: eps = 1.E-6_wp

! Eingabeparameter
! c_adv : Advektionsgeschwindigkeit
! phi: Differenzenstern der advehierten Groesse 

REAL   (KIND=wp),     INTENT(IN) :: c_adv
REAL   (KIND=wp),     INTENT(IN) :: phi(5)


REAL   (KIND=wp)                 :: v(5)
REAL   (KIND=wp)                 :: smeas(3), omega(3), q(3)
REAL   (KIND=wp)                 :: alpha(0:3)

!------------------------------------------------------------------------------

  v = c_adv * phi

  smeas(1) = 13._wp * ( v(1) - 2._wp*v(2) + v(3) )**2 &
            + 3._wp * ( v(1) - 4._wp*v(2) + 3._wp*v(3) )**2
  smeas(2) = 13._wp * ( v(2) - 2._wp*v(3) + v(4) )**2 &
            + 3._wp * ( v(2) - v(4) )**2
  smeas(3) = 13._wp * ( v(3) - 2._wp*v(4) + v(5) )**2 &
            + 3._wp * ( 3._wp*v(3) - 4._wp*v(4) + v(5) )**2

  alpha(1) = 1._wp / ( smeas(1) + eps )**2
  alpha(2) = 6._wp / ( smeas(2) + eps )**2
  alpha(3) = 3._wp / ( smeas(3) + eps )**2
  alpha(0) = 1._wp / ( alpha(1) + alpha(2) + alpha(3) )

  omega(1) = alpha(0) * alpha(1)
  omega(2) = alpha(0) * alpha(2)
  omega(3) = alpha(0) * alpha(3)

  q(1) =  1._wp/3._wp*v(1) - 7._wp/6._wp*v(2) + 11._wp/6._wp*v(3)
  q(2) = -1._wp/6._wp*v(2) + 5._wp/6._wp*v(3) + 1._wp/3._wp*v(4)
  q(3) =  1._wp/3._wp*v(3) + 5._wp/6._wp*v(4) - 1._wp/6._wp*v(5)

  !
  ! gewichteten Fluss des ENO-Verfahrens berechnen:
  !
  flux_weno3_priv = omega(1)*q(1) + omega(2)*q(2) + omega(3)*q(3)

END FUNCTION flux_weno3_priv

!==============================================================================
!==============================================================================

SUBROUTINE clipping( s, ie, je, ke, eps )

!--------------------------------------------------------------------------
!
! Description:
!   Clipping; reduction of undershootings
!
! Method:
!   negative values are set to 0
!   
!--------------------------------------------------------------------------

! Declarations:

INTEGER          ,   INTENT(IN) ::  &
  ie, je, ke    ! x-, y- and z-dimension of the input/output arrays

REAL    (KIND=wp),   INTENT(IN), OPTIONAL ::  &
  eps           ! security parameter to avoid division by zero

REAL    (KIND=wp),   INTENT(INOUT) ::  &
  s(ie,je,ke)   ! scalar var. to be transported at time-level n-1

! Local variables
!----------------

INTEGER         ::  &
  i, j, k              ! loop indices

  IF ( PRESENT(eps) ) THEN
    WHERE ( s(:,:,:) < eps )
      s = eps
    END WHERE
  ELSE
    WHERE ( s(:,:,:) < 0.0_wp )
      s = 0.0_wp
    END WHERE
  END IF

END SUBROUTINE clipping

!==============================================================================
!==============================================================================

SUBROUTINE multiplicative_filling( feld, clipping_typ )

!------------------------------------------------------------------------------
!
! Description:
!
! Clipping of negative values + global conservation
! (=driver routine for Subr. 'integral_3d' and 'clipping_2')
! Lit.: Rood (1987)
!
!------------------------------------------------------------------------------

USE environment,        ONLY: model_abort
USE parallel_utilities, ONLY: global_values

USE data_parallel,      ONLY :  &
    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed
    my_cart_id,      & !
    icomm_cart,      & ! communicator for the virtual cartesian topology
    num_compute,     & ! number of compute PEs 
    imp_reals          ! determines the correct REAL type used in the model

USE data_modelconfig,   ONLY :   &
    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ke              ! number of grid points in vertical direction

USE data_runcontrol,   ONLY :   &
    idbg_level, ldebug_dyn

REAL (KIND=wp),     INTENT(inout) :: feld(:,:,:)
INTEGER, INTENT(in)               :: clipping_typ

REAL (KIND=wp)     :: z_integ, factor, clip_limit
REAL (KIND=wp)     :: z_neg_integ    ! has a positive sign
CHARACTER (LEN=80) :: yzerrmsg
INTEGER            :: izerror

! End of header
!==============================================================================

  IF ( ldebug_dyn .AND. ( idbg_level >= 10 ) .AND. ( my_cart_id==0 ) ) THEN
    WRITE(*,*) "[Subr. multiplicative_filling ...]"
  END IF

  z_integ = integral_3d_wg ( feld, nboundlines+1, ie-nboundlines,   &
                                   nboundlines+1, je-nboundlines,   &
                                               1,             ke)

  IF ( clipping_typ == 6 ) THEN
    ! Clipping to 0
    CALL  clipping_2( feld, 0.0_wp, &
                      nboundlines+1, ie-nboundlines,   &
                      nboundlines+1, je-nboundlines,   &
                                  1,           ke  , z_neg_integ)
  ELSE IF ( clipping_typ == 7 ) THEN
    ! Clipping of slightly positive values, too
    ! Proposal of Rood (1987)
    clip_limit = MINVAL( feld )
    IF (num_compute > 1) THEN
      CALL global_values( clip_limit, 1, 'MIN', imp_reals, icomm_cart ,-1, yzerrmsg, izerror )
    ENDIF

    CALL  clipping_2( feld, -clip_limit/2,               &
                      nboundlines+1, ie-nboundlines,     &
                      nboundlines+1, je-nboundlines,     &
                                  1,           ke  , z_neg_integ)
  ELSE
    yzerrmsg = "false value in clipping_typ!"
    CALL model_abort( my_cart_id,  1, yzerrmsg, "multiplicative_filling")
  END IF

  IF (num_compute > 1) THEN
    CALL global_values( z_integ,     1, 'SUM', imp_reals, icomm_cart ,-1, yzerrmsg, izerror )
    CALL global_values( z_neg_integ, 1, 'SUM', imp_reals, icomm_cart ,-1, yzerrmsg, izerror )
  ENDIF

  IF ( z_integ >= 0.0_wp ) THEN

    factor = z_integ + z_neg_integ
    IF ( factor > 1.0e-20_wp ) THEN
      factor = z_integ / factor
    ELSE
      factor = 1.0_wp
    END IF

    !IF ( my_cart_id == 0 ) THEN
    !  WRITE(*,'(A,3E13.5)') " neg. values factor= ", factor, z_integ, z_neg_integ
    !END IF

    ! --- multiplicative reduction of values ---
    feld(:,:,:) = factor * feld(:,:,:)

  ELSE
    ! pathological case; nearest to global conservation is:
    feld(:,:,:) = 0.0_wp

  END IF

END SUBROUTINE multiplicative_filling

!==============================================================================
!==============================================================================


SUBROUTINE multiplicative_filling_DDI( feld, clipping_typ )

!------------------------------------------------------------------------------
!
! Description:
!
! Clipping of negative values + global conservation
! (=driver routine for Subr. 'integral_3d' and 'clipping_2')
! Lit.: Rood (1987)
! this is done in a domain decomposition invariant (DDI) manner
!
!------------------------------------------------------------------------------

USE environment,        ONLY: model_abort
USE parallel_utilities, ONLY: global_values

USE data_parallel,      ONLY :  &
  & nboundlines,     & ! number of boundary lines of the domain for which
  !                      no forecast is computed
  & my_cart_id,      & !
  & icomm_cart,      & ! communicator for the virtual cartesian topology
  & imp_reals          ! determines the correct REAL type used in the model

USE data_modelconfig,   ONLY :   &
  & ie,           & ! number of grid points in zonal direction
  & je,           & ! number of grid points in meridional direction
  & ke              ! number of grid points in vertical direction

USE data_runcontrol,   ONLY :   &
    idbg_level, ldebug_dyn

REAL (KIND=wp),     INTENT(inout) :: feld(:,:,:)
INTEGER, INTENT(in)               :: clipping_typ

REAL (KIND=wp)     :: z_integ_glob, factor, clip_limit
REAL (KIND=wp)     :: z_neg_integ_glob   ! has a positive sign
CHARACTER (LEN=80) :: yzerrmsg
INTEGER            :: izerror

! End of header
!==============================================================================

  IF ( ldebug_dyn .AND. ( idbg_level >= 10 ) .AND. ( my_cart_id==0 ) ) THEN
    WRITE(*,*) "[Subr. multiplicative_filling_DDI ...]"
  END IF

  z_integ_glob  = integral_3d_wg_DDI (feld, nboundlines+1, ie-nboundlines,   &
                                            nboundlines+1, je-nboundlines,   &
                                                        1,             ke)

  IF ( clipping_typ == 6 ) THEN
    ! Clipping to 0
    CALL  clipping_DDI( feld, 0.0_wp,                &
                        nboundlines+1, ie-nboundlines,   &
                        nboundlines+1, je-nboundlines,   &
                                    1,             ke, z_neg_integ_glob)
  ELSE IF ( clipping_typ == 7 ) THEN
    ! Clipping of slightly positive values, too
    ! Proposal of Rood (1987)
    clip_limit = MINVAL( feld )
    CALL global_values( clip_limit, 1, 'MIN', imp_reals, icomm_cart ,-1, yzerrmsg, izerror )

    CALL  clipping_DDI( feld, -clip_limit/2,               &
                        nboundlines+1, ie-nboundlines,     &
                        nboundlines+1, je-nboundlines,     &
                                    1,             ke, z_neg_integ_glob)
  ELSE
    yzerrmsg = "false value in clipping_typ!"
    CALL model_abort( my_cart_id,  1, yzerrmsg, "multiplicative_filling_DDI")
  END IF

  !not longer required:
  !CALL global_values( z_integ,     1, 'SUM', imp_reals, icomm_cart ,-1, yzerrmsg, izerror )
  !CALL global_values( z_neg_integ, 1, 'SUM', imp_reals, icomm_cart ,-1, yzerrmsg, izerror )

  IF ( z_integ_glob >= 0.0_wp ) THEN

    factor = z_integ_glob + z_neg_integ_glob
    IF ( factor > 1.0e-20_wp ) THEN
      factor = z_integ_glob / factor
    ELSE
      factor = 1.0_wp
    END IF

    ! --- multiplicative reduction of values ---
    feld(:,:,:) = factor * feld(:,:,:)
  ELSE
    ! pathological case; nearest to global conservation is:
    feld(:,:,:) = 0.0_wp

  END IF

END SUBROUTINE multiplicative_filling_DDI

!==============================================================================
!==============================================================================

REAL (KIND=wp)     FUNCTION integral_3d_wg( feld, istart, iend, jstart, jend, kstart, kend )

!------------------------------------------------------------------------------
!
! Description:
!
! calculate the integral of 'feld' in the domain of the processor
! here: sqrt(g) = measure of volume of the grid box
! is already calculated
!
!------------------------------------------------------------------------------

USE data_fields,      ONLY: crlat
USE grid_metrics_utilities, ONLY: sqrtg_r_s
USE data_constants,   ONLY: r_earth, pi
USE data_modelconfig, ONLY :   &
  & dlon,         & ! grid point distance in zonal direction (in degrees)
  & dlat            ! grid point distance in meridional direction (in degrees)

REAL (KIND=wp), INTENT(in) :: feld(:,:,:)
INTEGER       , INTENT(in) :: istart, iend, jstart, jend, kstart, kend
REAL (KIND=wp), PARAMETER  :: dzeta = 1.0_wp
REAL (KIND=wp)             :: h1, h2
INTEGER                    :: i, j, k

! End of header
!==============================================================================

  integral_3d_wg  = 0.0_wp

  DO j = jstart, jend
    h1 = 0.0_wp
      DO k = kstart, kend
      h2 = 0.0_wp
      DO i = istart, iend
        h2 = h2 + feld(i,j,k) / sqrtg_r_s(i,j,k)
      END DO
      h1 = h1 + h2
    END DO
    integral_3d_wg = integral_3d_wg  + h1 * crlat(j,1)
  END DO

  integral_3d_wg = integral_3d_wg         &
            * r_earth**2 * (pi/180.0_wp)**2 * dlon * dlat * dzeta

END FUNCTION integral_3d_wg

!==============================================================================
!==============================================================================

REAL (KIND=wp)     FUNCTION integral_3d_wg_DDI( feld, istart, iend,        &
                          jstart, jend, kstart, kend )

!------------------------------------------------------------------------------
!
! Description:
! calculate the integral of 'feld' in the domain of the processor
! The result is the global, domain decomposition invariant (DDI) integral
! here: sqrt(g) = measure of volume of the grid box
! is already calculated
!
!------------------------------------------------------------------------------

USE data_parallel,      ONLY : &
  my_cart_id

USE data_fields,      ONLY: crlat
USE grid_metrics_utilities, ONLY: sqrtg_r_s
USE data_constants,   ONLY: r_earth, pi
USE data_modelconfig, ONLY :   &
  & dlon,         & ! grid point distance in zonal direction (in degrees)
  & dlat            ! grid point distance in meridional direction (in degrees)

REAL    (KIND=wp), INTENT(in) :: feld(:,:,:)
INTEGER          , INTENT(in) :: istart, iend, jstart, jend, kstart, kend

REAL    (KIND=wp), PARAMETER  :: dzeta = 1.0_wp
INTEGER                              :: i, j, k
INTEGER                              :: istat
REAL    (KIND=wp), ALLOCATABLE :: vert_integral(:,:)

REAL    (KIND=wp)  :: sum2
CHARACTER (LEN=80) :: yzerrmsg
INTEGER            :: izerror

! End of header
!==============================================================================

  ALLOCATE(  vert_integral( istart:iend, jstart:jend ), STAT=istat )

  vert_integral(:,:) = 0.0_wp

  DO j = jstart, jend
    DO i = istart, iend

      DO k = kstart, kend
        vert_integral(i,j) = vert_integral(i,j)             &
          + feld(i,j,k) / sqrtg_r_s(i,j,k) * crlat(j,1)
      END DO

    END DO
  END DO

  integral_3d_wg_DDI = sum_DDI( vert_integral )

  integral_3d_wg_DDI = integral_3d_wg_DDI                &
    * r_earth**2 * (pi/180.0_wp)**2 * dlon * dlat * dzeta

  DEALLOCATE(  vert_integral, STAT=istat )

END FUNCTION integral_3d_wg_DDI

!==============================================================================
!==============================================================================

SUBROUTINE clipping_2 (feld, clip_limit, istart, iend, jstart, jend,    &
                  kstart, kend, int_neg_values )

!------------------------------------------------------------------------------
!
! Description:
! Clipping of negative values and also of small positive values
! (to balance approximately mass conservation (Rood, 1987) in 'feld'
! The integral of the clipped mass is delivered in 'int_neg_values'
!
!------------------------------------------------------------------------------

USE data_fields,      ONLY: crlat
USE grid_metrics_utilities, ONLY: sqrtg_r_s
USE data_constants,   ONLY: r_earth, pi
USE data_modelconfig, ONLY :   &
  & dlon,         & ! grid point distance in zonal direction (in degrees)
  & dlat            ! grid point distance in meridional direction (in degrees)

REAL (KIND=wp), INTENT(inout) :: feld(:,:,:)
INTEGER       , INTENT(in)    :: istart, iend, jstart, jend, kstart, kend
REAL (KIND=wp), INTENT(in)    :: clip_limit
REAL (KIND=wp), INTENT(out)   :: int_neg_values

REAL (KIND=wp), PARAMETER     :: dzeta = 1.0_wp
INTEGER                       :: i, j, k

! End of header
!==============================================================================

  int_neg_values = 0.0_wp

  DO k = kstart, kend
    DO j = jstart, jend
      DO i = istart, iend

        IF ( feld(i,j,k) < clip_limit )  THEN
          int_neg_values = int_neg_values - feld(i,j,k) * crlat(j,1)/ sqrtg_r_s(i,j,k)
          feld(i,j,k) = 0.0_wp
        END IF

      END DO
    END DO
  END DO

  int_neg_values = int_neg_values        &
           * r_earth**2 * (pi/180.0_wp)**2 * dlon * dlat * dzeta

END SUBROUTINE  clipping_2

!==============================================================================
!==============================================================================

SUBROUTINE clipping_DDI (feld, clip_limit, istart, iend, jstart, jend, &
                         kstart, kend, int_neg_values_glob )

!------------------------------------------------------------------------------
!
! Description:
! Clipping of negative values and also of small positive values
! (to balance approximately mass conservation (Rood, 1987) in 'feld'
! The DDI-integral of the clipped mass is delivered in 'int_neg_values_glob'
!
!------------------------------------------------------------------------------

USE data_fields,      ONLY: crlat
USE grid_metrics_utilities, ONLY: sqrtg_r_s
USE data_constants,   ONLY: r_earth, pi
USE data_modelconfig, ONLY :   &
  & dlon,         & ! grid point distance in zonal direction (in degrees)
  & dlat            ! grid point distance in meridional direction (in degrees)

REAL    (KIND=wp), INTENT(inout) :: feld(:,:,:)
INTEGER          , INTENT(in)    :: istart, iend, jstart, jend, kstart, kend
REAL    (KIND=wp), INTENT(in)    :: clip_limit
REAL    (KIND=wp), INTENT(out)   :: int_neg_values_glob

REAL    (KIND=wp), ALLOCATABLE   :: vert_integral(:,:)
REAL    (KIND=wp), PARAMETER     :: dzeta = 1.0_wp
INTEGER                          :: i, j, k
INTEGER                          :: istat

! End of header
!==============================================================================

  ALLOCATE(  vert_integral( istart:iend, jstart:jend ), STAT=istat )

  vert_integral(:,:) = 0.0_wp

  DO j = jstart, jend
    DO i = istart, iend

      DO k = kstart, kend

        IF ( feld(i,j,k) < clip_limit )  THEN
          vert_integral(i,j) = vert_integral(i,j)               &
                    - feld(i,j,k) * crlat(j,1)/ sqrtg_r_s(i,j,k)
          feld(i,j,k) = 0.0_wp
        END IF

      END DO

    END DO
  END DO

  int_neg_values_glob = sum_DDI( vert_integral )

  int_neg_values_glob = int_neg_values_glob  &
           * r_earth**2 * (pi/180.0_wp)**2 * dlon * dlat * dzeta

  DEALLOCATE(  vert_integral, STAT=istat )

END SUBROUTINE  clipping_DDI

!==============================================================================
!==============================================================================

REAL (KIND=wp) FUNCTION  sum_DDI( field_2d )

!------------------------------------------------------------------------------
!
! Description:
! DDI-summation of the 2-dim. field 'field_2D'
! i.e. the result is independent (or invariant)
! of the domain decomposition
!
!------------------------------------------------------------------------------

USE  data_modelconfig,  ONLY:    &
                      ie_tot, je_tot, ie, je, istart, iend, jstart, jend
USE data_parallel,      ONLY:  &
  imp_reals, imp_byte, imp_integers, & !
  icomm_cart,      & !
  nboundlines,     & ! number of boundary lines of the domain for which
                     ! no forecast is computed
  num_compute        ! number of compute PEs
USE parallel_utilities, ONLY:  &
  global_values,  & !
  gather_values

REAL    (KIND=wp), INTENT(IN)    :: field_2D(istart:iend, jstart:jend)

REAL    (KIND=wp)         :: field_max, field_min, df
INTEGER                   :: z_nmbr_gridpoints
INTEGER                   :: i, j
INTEGER (KIND=i8)         :: i_quant_range, quant_sum, h
INTEGER (KIND=i8)         :: iConst_2_32, iConst_2_31, iConst_2_16
INTEGER                   :: i_low, i_high
INTEGER, ALLOCATABLE      :: i_low_vec(:), i_high_vec(:)
LOGICAL                   :: is_8Byte_int_compil_implemented
LOGICAL                   :: is_8Byte_int_MPI_implemented
CHARACTER (LEN=80)        :: yzerrmsg
INTEGER                   :: izerror
INTEGER                   :: istat
REAL    (KIND=wp)         :: eps

! End of header
!==============================================================================

  eps = 1.0E-30_wp

  is_8Byte_int_MPI_implemented = .FALSE.
  iConst_2_16 = 256 * 256  ! = 2**16

  ! summation only in the interior of the domain:
  z_nmbr_gridpoints = ( ie_tot - 2 * nboundlines ) &
                    * ( je_tot - 2 * nboundlines )


  ! --- is 8 Byte-integer arithmetic implemented? ----
  IF ( HUGE( i_quant_range ) / ( iConst_2_16 * 256 ) > 1 ) THEN
    is_8Byte_int_compil_implemented = .TRUE.
    iConst_2_32 = iConst_2_16 * iConst_2_16     ! = 2**32
    iConst_2_31 = iConst_2_32 / 2               ! = 2**31
  ELSE
    is_8Byte_int_compil_implemented = .FALSE.
    PRINT*, "WARNING in sum_DDI: accuracy probably to low!"

  END IF
  ! --- determine the range of values in field_2D ----
  field_max = MAXVAL( field_2D )
  field_min = MINVAL( field_2D )

  IF (num_compute > 1) THEN
    CALL global_values( field_max, 1, 'MAX', imp_reals, icomm_cart ,-1, yzerrmsg, izerror )
    CALL global_values( field_min, 1, 'MIN', imp_reals, icomm_cart ,-1, yzerrmsg, izerror )
  ENDIF

  ! --- number of quantisation steps ---
  IF ( is_8Byte_int_compil_implemented  ) THEN
    i_quant_range = ( HUGE( i_quant_range ) / 2 ) / z_nmbr_gridpoints
  ELSE
    ! alternative solution not yet implemented;
    ! the following is probably too inaccurate:
    i_quant_range = ( HUGE( i_quant_range ) / 2 ) / z_nmbr_gridpoints
  END IF

  ! --- Quantisation of the field-elements and summation ---
  df = ( field_max - field_min ) / i_quant_range

  IF ( df < eps ) THEN
    df = eps
  END IF

  quant_sum = 0

  DO j = jstart, jend
    DO i = istart, iend
      h = INT ( (field_2D(i,j) - field_min ) / df + 0.5_wp, i8 )
      quant_sum = quant_sum + h
    END DO
  END DO

  IF ( is_8Byte_int_MPI_implemented  ) THEN
    CALL global_values( quant_sum, 1, 'SUM', 8*imp_byte, icomm_cart ,-1, &
                          yzerrmsg, izerror )
  ELSE

    ALLOCATE( i_low_vec (1:num_compute), STAT=istat )
    ALLOCATE( i_high_vec(1:num_compute), STAT=istat )

    i_high = INT(quant_sum / iConst_2_31)
    i_low  = INT(quant_sum - i_high * iConst_2_31)

    IF (num_compute > 1) THEN
      CALL gather_values ( i_low,  i_low_vec,  1, num_compute , imp_integers, -1,  &
                                icomm_cart, yzerrmsg, izerror)
      CALL gather_values ( i_high, i_high_vec, 1, num_compute , imp_integers, -1,  &
                                icomm_cart, yzerrmsg, izerror)
    ELSE
      ! single processor solution
      i_high_vec(1) = i_high
      i_low_vec (1) = i_low
    ENDIF

    quant_sum = 0
    DO i=1, num_compute
       quant_sum = quant_sum + ( i_high_vec(i) * iConst_2_31 + i_low_vec(i) )
    END DO

    DEALLOCATE( i_low_vec,  STAT=istat)
    DEALLOCATE( i_high_vec, STAT=istat)

  END IF

  ! --- backtransformation into real range ---

  sum_DDI = quant_sum * df + field_min * z_nmbr_gridpoints

END FUNCTION sum_DDI

!==============================================================================

END MODULE numeric_utilities_rk
