!+ Module for metric utilities
!------------------------------------------------------------------------------

MODULE grid_metrics_utilities

!------------------------------------------------------------------------------
!
! Description:
!   This module generates all the data describing the grid and the 
!   coordinate transformation.
!   The module contains only subroutines which don't need any 
!   other module (only exception: data_parameters).
!   These routines are called by 'init_grid_metrics' (possible
!   boundary exchange, periodic BC's or 2D-simulations
!   are also treated there)
!   At the end of the simulation, a finalization routine is called.
!
! Method:
!
! Current Code Owner: DWD, Michael Baldauf
!  phone:  +49  69  8062 2733
!  fax:    +49  69  8236 1493
!  email:  michael.baldauf@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! V4_23        2012/05/10 Michael Baldauf
!  Initial version: only some subroutines are included now, not to 
!  shuffle too much code around in one version change
! V5_1         2014-11-28 Ulrich Blahak
!  Removed l2dim-treatment, because this is done now within the call to
!  exchg_boundaries, which must be done in any case, also for sequential 
!  applications to treat periodic boundaries correct. 
!  Therefore removed conditional compilation for MPI
!  Replaced ireals by wp (working precision) (OF)
! V5_2         2015-05-21 Ulrich Schaettler
!  Implemented finalization routine final_grid_metrics, to deallocate fields
!  Reorganized numbering of MPI datatypes ID
! V5_3         2015-10-09 Michael Baldauf
!  New subroutine calc_dzeta_dxx_at_y, deallocate_dzeta_dxx_at_y
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

USE data_parameters,  ONLY:  &
  wp,        & ! KIND-type parameter for real variables
  iintegers    ! KIND-type parameter for standard integer variables

USE data_runcontrol,  ONLY:  &
  l2tls,            & ! forecast with 2-TL integration scheme
  l_dzeta_d_needed, &
  l2dim,            &
  lperi_x,          & ! if lartif_data=.TRUE.: periodic boundary conditions in x-dir.
                              ! or with Davies conditions (.FALSE.)
  lperi_y             ! if lartif_data=.TRUE.: periodic boundary conditions in y-dir.
                      ! or with Davies conditions (.FALSE.)

USE data_modelconfig, ONLY:  &
  ie,           & ! number of grid points in zonal direction
  je,           & ! number of grid points in meridional direction
  ke,           & ! number of grid points in vertical direction
  ke1,          & ! KE+1
  jstart,       & ! start index for the forecast of w, t, qv, qc and pp
  jend,         & ! end index for the forecast of w, t, qv, qc and pp
  jstartu,      & ! start index for the forecast of u
  jendu,        & ! end index for the forecast of u
  jstartv,      & ! start index for the forecast of v
  jendv,        & ! end index for the forecast of v
  jstartpar,    & ! start index for computations in the parallel program
  jendpar,      & ! end index for computations in the parallel program
  dlon,         & ! grid point distance in zonal direction (in degrees)
  dlat,         & ! grid point distance in meridional direction (in degrees)
  eddlon,       &
  eddlat,       &
  dt

USE data_fields,      ONLY:  &
  hhl               ! geometical height of half model levels        ( m )

USE data_parallel,    ONLY:  &
  my_cart_id,       & ! rank of this subdomain in the cartesian communicator
  my_cart_neigh,    & ! neighbors
  num_compute,      & ! number of compute PEs
  nboundlines,      & ! number of boundary lines of the domain for which
                      ! no forecast is computed = overlapping boundary
                      ! lines of the subdomains
  icomm_cart,       & ! communicator that belongs to the cartesian grid
  imp_reals,        & ! determines the correct REAL type used in the model
                      ! for MPI
  ncomm_type,       & ! type of communication
  ldatatypes,       & ! if .TRUE.: use MPI-Datatypes for some communications
  sendbuf,          & ! sending buffer for boundary exchange:
                      ! 1-4 are used for sending, 5-8 are used for receiving
  isendbuflen         ! length of one column of sendbuf

USE environment,      ONLY:  exchg_boundaries

!==============================================================================

IMPLICIT NONE

!==============================================================================

  REAL  (KIND=wp),     ALLOCATABLE  ::  &
    sqrtg_r_s (:,:,:),   & ! 1 / square root of G at scalar points       (1/m)
    sqrtg_r_u (:,:,:),   & ! 1 / square root of G at u points            (1/m)
    sqrtg_r_v (:,:,:),   & ! 1 / square root of G at v points            (1/m)
    sqrtg_r_w (:,:,:),   & ! 1 / square root of G at w points            (1/m)
    dzeta_dlam(:,:,:),   & ! d zeta / d lambda (for constant phi,    z)
                           ! at the scalar position                      ( 1 )
    dzeta_dphi(:,:,:)      ! d zeta / d phi    (for constant lambda, z)
                           ! at the scalar position                      ( 1 )

  REAL (KIND=wp),      ALLOCATABLE  ::  &
    dzeta_dlam_at_u (:,:,:), & ! d zeta / d lambda (for constant phi,z) at u-position
    dzeta_dlam_at_w (:,:,:), & ! d zeta / d lambda (for constant phi,z) at w-position
    dzeta_dlam_at_uw(:,:,:), & ! d zeta / d lambda (for constant phi,z) at uw-position
    dzeta_dphi_at_v (:,:,:), & ! d zeta / d phi    (for constant lambda,z) at v-position
    dzeta_dphi_at_w (:,:,:), & ! d zeta / d phi    (for constant lambda,z) at w-position
    dzeta_dphi_at_vw(:,:,:)    ! d zeta / d phi    (for constant lambda,z) at vw-position

  REAL  (KIND=wp),     ALLOCATABLE  ::  &
    wgtfac   (:,:,:),    & ! weighting factor for vertical interpolation ( 1 )
    wgtfac_u (:,:,:),    & ! weighting factor for vertical interpolation ( 1 )
    wgtfac_v (:,:,:),    & ! weighting factor for vertical interpolation ( 1 )
    wgtfacq  (:,:,:),    & ! weighting factor for vertical interpolation ( 1 )
    wgtfacq_u(:,:,:),    & ! weighting factor for vertical interpolation ( 1 )
    wgtfacq_v(:,:,:),    & ! weighting factor for vertical interpolation ( 1 )
    wgt_one_sided(:,:,:)   ! weighting factor for vertical interpolation ( 1 )

!==============================================================================

CONTAINS

!==============================================================================
!==============================================================================

SUBROUTINE init_grid_metrics (ierrorstat)

!------------------------------------------------------------------------------
!
! Description:
!   call of the routines which compute the necessary metric variables
!
!------------------------------------------------------------------------------

  INTEGER (KIND=iintegers), INTENT(OUT) :: ierrorstat

  INTEGER :: istat

  INTEGER (KIND=iintegers) :: kzdims(24)
  INTEGER (KIND=iintegers) :: izerror
  CHARACTER(100)           :: yzerrmsg

!------------------------------------------------------------------------------

  ierrorstat = 0

  istat = 0

  IF ( my_cart_id == 0 ) WRITE(*,*) "Subr.[init_grid_metrics] ..."

! Allocation of fields: are deallocated in final_grid_metrics
  IF (.NOT. ALLOCATED(sqrtg_r_s)) ALLOCATE ( sqrtg_r_s(ie,je,ke), STAT=istat )
  sqrtg_r_s=0.0_wp

  IF (.NOT. ALLOCATED(sqrtg_r_u)) ALLOCATE ( sqrtg_r_u(ie,je,ke), STAT=istat )
  sqrtg_r_u=0.0_wp

  IF (.NOT. ALLOCATED(sqrtg_r_v)) ALLOCATE ( sqrtg_r_v(ie,je,ke), STAT=istat )
  sqrtg_r_v=0.0_wp

  IF (.NOT. ALLOCATED(sqrtg_r_w)) ALLOCATE ( sqrtg_r_w(ie,je,ke1),STAT=istat )
  sqrtg_r_w=0.0_wp

  IF ( l_dzeta_d_needed ) THEN
    IF (.NOT. ALLOCATED(dzeta_dlam)) ALLOCATE ( dzeta_dlam(ie,je,ke), STAT=istat )
    dzeta_dlam=0.0_wp

    IF (.NOT. ALLOCATED(dzeta_dphi)) ALLOCATE ( dzeta_dphi(ie,je,ke), STAT=istat )
    dzeta_dphi=0.0_wp

  END IF

  IF (.NOT. ALLOCATED(wgtfac)) ALLOCATE ( wgtfac    (ie,je,ke1), STAT=istat )
  wgtfac    = 0.0_wp

  IF (.NOT. ALLOCATED(wgtfac_u)) ALLOCATE ( wgtfac_u  (ie,je,ke1), STAT=istat )
  wgtfac_u  = 0.0_wp

  IF (.NOT. ALLOCATED(wgtfac_v)) ALLOCATE ( wgtfac_v  (ie,je,ke1), STAT=istat )
  wgtfac_v  = 0.0_wp

  IF (.NOT. ALLOCATED(wgtfacq)) ALLOCATE ( wgtfacq   (ie,je,3),   STAT=istat )
  wgtfacq   = 0.0_wp

  IF (.NOT. ALLOCATED(wgtfacq_u)) ALLOCATE ( wgtfacq_u (ie,je,3),   STAT=istat )
  wgtfacq_u = 0.0_wp

  IF (.NOT. ALLOCATED(wgtfacq_v)) ALLOCATE ( wgtfacq_v (ie,je,3),   STAT=istat )
  wgtfacq_v = 0.0_wp

!------------------------------------------------------------------------------

  IF ( l2tls ) THEN
    CALL weighting_factors_full2half( hhl,                                &
              wgtfac, wgtfac_u, wgtfac_v, wgtfacq, wgtfacq_u, wgtfacq_v,  &
              ie, je, ke )
  END IF

  kzdims(1:24) =                           &
       (/ ke+1, ke+1, ke+1,   3,   3,      &
             3,    0,    0,   0,   0,      &
             0,    0,    0,   0,   0,      &
             0,    0,    0,   0,   0,      &
             0,    0,    0,   0  /)

  CALL exchg_boundaries                                                        &
       ( 9, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,  &
        kzdims, jstartpar, jendpar, nboundlines, nboundlines, my_cart_neigh,   &
        lperi_x, lperi_y, l2dim,                                               &
        0, ldatatypes, ncomm_type, izerror, yzerrmsg,                          &
        wgtfac   (:,:,:),     &
        wgtfac_u (:,:,:),     &
        wgtfac_v (:,:,:),     &
        wgtfacq  (:,:,:),     &
        wgtfacq_u(:,:,:),     &
        wgtfacq_v(:,:,:) )

  IF ( izerror /= 0_iintegers ) THEN
    ierrorstat = 20005
  END IF

  ! --- vertical weighting factors (2) ------------------------------

  IF (.NOT. ALLOCATED(wgt_one_sided)) ALLOCATE( wgt_one_sided (ie,je,3), STAT=istat )

  CALL weighting_factors_one_sided( hhl, wgt_one_sided, ie, je, ke )


  ! --- metric coeeficients of the terrain following coordinate -----

  ! initialize sqrt(G) for all time level schemes
  CALL calc_sqrtg_r (hhl, sqrtg_r_s, sqrtg_r_u, sqrtg_r_v, sqrtg_r_w, &
                     ie, je, ke )

  kzdims(1:24) =                         &
       (/ ke , ke ,  ke, ke+1,   0,      &
            0,   0,   0,    0,   0,      &
            0,   0,   0,    0,   0,      &
            0,   0,   0,    0,   0,      &
            0,   0,   0,    0  /)

  CALL exchg_boundaries                                                          &
       (10, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,    &
        kzdims, jstartpar, jendpar, nboundlines, nboundlines, my_cart_neigh,     &
        lperi_x, lperi_y, l2dim,                                                 &
        0, ldatatypes, ncomm_type, izerror, yzerrmsg,                            &
        sqrtg_r_s(:,:,:),     &
        sqrtg_r_u(:,:,:),     &
        sqrtg_r_v(:,:,:),     &
        sqrtg_r_w(:,:,:) )

  IF ( izerror /= 0_iintegers ) THEN
    ierrorstat = 20005
  END IF

  IF ( l_dzeta_d_needed ) THEN

    CALL metric_coeffs ( hhl, sqrtg_r_s, eddlon, eddlat,       &
      dzeta_dlam, dzeta_dphi, ie, je, ke )


    kzdims(1:24) =                        &
         (/ ke , ke ,   0,   0,   0,      &
              0,   0,   0,   0,   0,      &
              0,   0,   0,   0,   0,      &
              0,   0,   0,   0,   0,      &
              0,   0,   0,   0  /)

    CALL exchg_boundaries                                                          &
         (11, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,    &
          kzdims, jstartpar, jendpar, nboundlines, nboundlines, my_cart_neigh,     &
          lperi_x, lperi_y, l2dim,                                                 &
          0, ldatatypes, ncomm_type, izerror, yzerrmsg,                            &
          dzeta_dlam(:,:,:),     &
          dzeta_dphi(:,:,:)    )

    IF ( izerror /= 0_iintegers ) THEN
      ierrorstat = 20005
    END IF

  END IF

END SUBROUTINE init_grid_metrics

!==============================================================================
!==============================================================================

SUBROUTINE calc_sqrtg_r (hhl, sqrtg_r_s, sqrtg_r_u, sqrtg_r_v, sqrtg_r_w, &
                         ie, je, ke)

!------------------------------------------------------------------------------
!
! Description:
!   calculate the reciprocal square root of G (Dim.: [1/m]) at several grid positions
!   e.g. as a preliminary for Subr. 'integral_3D'.
!   This routine is only needed, if the 3-timelevel scheme is used
!   (because sqrtg_r_s, ... will not be calculated there)
!   A good place for calling this routine is e.g.
!   in Subr. 'init_dynamics' (file oranize_dynamics.f90)
!
! Method:
!
!------------------------------------------------------------------------------

  INTEGER (KIND=iintegers), INTENT(IN)  ::         &
    ie, je, ke              ! dimensions

  REAL (KIND = wp),         INTENT(IN)  ::         &
    hhl       (ie,je,ke+1)  ! height of half levels

  REAL (KIND = wp),         INTENT(OUT) ::         &
    sqrtg_r_s (ie,je,ke), & ! for scalar grid points
    sqrtg_r_u (ie,je,ke), & ! for u-grid points
    sqrtg_r_v (ie,je,ke), & ! for v-grid points
    sqrtg_r_w (ie,je,ke+1)  ! for w-grid points

  INTEGER (KIND=iintegers) :: i, j, k

!------------------------------------------------------------------------------

  ! for scalar , u- and v- points
  DO  k = 1, ke

    DO j = 1, je
      DO i = 1, ie
        sqrtg_r_s(i,j,k) = 1.0_wp / ( hhl(i,j,k)-hhl(i,j,k+1) )
      ENDDO
    ENDDO
    ! (--> no exchange necessary; only for periodic b.c.)

    DO j = 1, je
      DO i = 1, ie-1
        sqrtg_r_u(i,j,k) = 2.0_wp / ( hhl(i,j,k  )+hhl(i+1,j,k  )    &
                                         -hhl(i,j,k+1)-hhl(i+1,j,k+1) )
      ENDDO
      ! approx. for lateral boundaries:
      sqrtg_r_u(ie,j,k) = 1.0_wp / ( hhl(ie,j,k  )-hhl(ie,j,k+1) ) 
    ENDDO

    DO j = 1, je-1
      DO i = 1, ie
        sqrtg_r_v(i,j,k) = 2.0_wp / ( hhl(i,j,k  )+hhl(i,j+1,k  )    &
                                         -hhl(i,j,k+1)-hhl(i,j+1,k+1) )
      ENDDO
    ENDDO
    ! approx. for lateral boundaries:
    DO i = 1, ie
      sqrtg_r_v(i,je,k) = 1.0_wp / ( hhl(i,je,k)-hhl(i,je,k+1) )
    ENDDO

  ENDDO

  ! and for the w-points
  DO  k = 2, ke
    DO j = 1, je
      DO i = 1, ie
        ! version A:
         sqrtg_r_w(i,j,k) = 2.0_wp / ( hhl(i,j,k-1)-hhl(i,j,k+1) )

        ! version B:
        !sqrtg_r_w(i,j,k) =    wgtfac(i,j,k)   / ( hhl(i,j,k)   - hhl(i,j,k+1) )  &
        !      + (1.0_wp     - wgtfac(i,j,k) ) / ( hhl(i,j,k-1) - hhl(i,j,k)   )
      ENDDO
    ENDDO
  ENDDO

  ! one-sided differences in k=1 and k=ke+1
  DO j = 1, je
    DO i = 1, ie
      sqrtg_r_w(i,j,1)    = 1.0_wp / ( hhl(i,j,1 )-hhl(i,j,2   ) )
      sqrtg_r_w(i,j,ke+1) = 1.0_wp / ( hhl(i,j,ke)-hhl(i,j,ke+1) )
    ENDDO
  ENDDO

END SUBROUTINE calc_sqrtg_r

!==============================================================================
!==============================================================================

SUBROUTINE metric_coeffs( hhl, sqrtg_r_s, eddlon, eddlat,       &
                          dzeta_dlam, dzeta_dphi, ie, je, ke )

!------------------------------------------------------------------------------
!
! Description:
!   calculate the metric coefficients (see Doms, Schaettler (2002) section 3):
!     d zeta / d lambda (at constant phi,    z) = - dzeta/dz * dz/dlambda
!     d zeta / d phi    (at constant lambda, z) = - dzeta/dz * dz/dphi
!   at the scalar grid positions.
!   (remarks:
!     dzeta/dz = 1/( dz/dzeta ) = -1/sqrt(G) (see subroutine calc_sqrtg_r)
!     dz/dzeta = -sqrt(G)  )
!   This subroutine is closely related to calc_sqrtg_r and should be
!   called after that.
!
! Method:
!   centered differences using hhl for calculating dz/d...
!
!------------------------------------------------------------------------------

  INTEGER (KIND=iintegers), INTENT(IN)  ::         &
    ie, je, ke              ! dimensions

  REAL (KIND = wp),         INTENT(IN)  ::         &
    hhl       (ie,je,ke+1)  ! height of half levels

  REAL (KIND = wp),         INTENT(IN) ::          &
    sqrtg_r_s (ie,je,ke)    ! = 1 / sqrt(G)  (at scalar grid point position)

  REAL (KIND = wp),         INTENT(IN) ::          &
    eddlon,   & ! 1 / dlon, dlon in degrees (!)    &
    eddlat      ! 1 / dlat, dlat in degrees (!)

  REAL (KIND=wp),     DIMENSION(1:ie,1:je,1:ke), INTENT(OUT) ::   &
    dzeta_dlam, dzeta_dphi

  INTEGER (KIND=iintegers) :: i, j, k

!------------------------------------------------------------------------------

  DO k = 1, ke
    DO j = 1, je
      DO  i = 2, ie-1

        dzeta_dlam(i,j,k) = sqrtg_r_s(i,j,k) *           &
          ( ( hhl(i+1,j,k  ) - hhl(i-1,j,k  ) )      &
          + ( hhl(i+1,j,k+1) - hhl(i-1,j,k+1) ) )    &
          * 0.25_wp * eddlon

      END DO
    END DO
  END DO
  ! approximation for lateral boundaries:
  dzeta_dlam(1, :,:) = 0.0_wp
  dzeta_dlam(ie,:,:) = 0.0_wp


  DO k = 1, ke
    DO j = 2, je-1
      DO  i = 1, ie

        dzeta_dphi(i,j,k) = sqrtg_r_s(i,j,k) *           &
          ( ( hhl(i,j+1,k  ) - hhl(i,j-1,k  ) )      &
          + ( hhl(i,j+1,k+1) - hhl(i,j-1,k+1) ) )    &
          * 0.25_wp * eddlat

      END DO
    END DO
  END DO
  ! approximation for lateral boundaries:
  dzeta_dphi(:,1 ,:) = 0.0_wp
  dzeta_dphi(:,je,:) = 0.0_wp

END SUBROUTINE metric_coeffs

!==============================================================================
!==============================================================================

SUBROUTINE calc_dzeta_dxx_at_y( ierrorstat )

!------------------------------------------------------------------------------
!
! Description:
! analogous to subr. 'metric_coeffs'
! allocate and calculate the following metric coefficient fields:
!   dzeta_dlam_at_u, dzeta_dlam_at_w, dzeta_dlam_at_uw,
!   dzeta_dphi_at_v, dzeta_dphi_at_w, dzeta_dphi_at_vw
! i.e. at different grid positions.
! (This subroutine is only needed for 3D-diffusion.)
!
! Method:
!   centered differences using hhl for calculating dz/d...
!
!------------------------------------------------------------------------------

  INTEGER (KIND=iintegers), INTENT(OUT) :: ierrorstat

  INTEGER :: i, j, k

  INTEGER (KIND=iintegers) :: kzdims(24)
  INTEGER (KIND=iintegers) :: izerror
  CHARACTER(100)           :: yzerrmsg

  ierrorstat = 0

  ALLOCATE( dzeta_dlam_at_u ( ie, je, ke   ) )
  ALLOCATE( dzeta_dlam_at_w ( ie, je, ke+1 ) )
  ALLOCATE( dzeta_dlam_at_uw( ie, je, ke+1 ) )

  ALLOCATE( dzeta_dphi_at_v ( ie, je, ke   ) )
  ALLOCATE( dzeta_dphi_at_w ( ie, je, ke+1 ) )
  ALLOCATE( dzeta_dphi_at_vw( ie, je, ke+1 ) )

  DO k = 1, ke+1
    DO j = 1, je
      DO  i = 1, ie-1

        dzeta_dlam_at_uw(i,j,k) =                                       &
           &    0.5_wp * ( sqrtg_r_w(i,j,k) + sqrtg_r_w(i+1,j,k) )      &
           &           * ( hhl(i+1,j,k) - hhl(i,j,k ) ) * eddlon

      END DO
    END DO
  END DO
  ! approximation for lateral boundaries:
  dzeta_dlam_at_uw(ie,:,:) = 0.0_wp

  DO k = 1, ke+1
    DO j = 1, je-1
      DO  i = 1, ie

        dzeta_dphi_at_vw(i,j,k) =                                       &
           &    0.5_wp * ( sqrtg_r_w(i,j,k) + sqrtg_r_w(i,j+1,k) )      &
           &           * ( hhl(i,j+1,k) - hhl(i,j,k ) ) * eddlat

      END DO
    END DO
  END DO
  ! approximation for lateral boundaries:
  dzeta_dphi_at_vw(:,je,:) = 0.0_wp

  DO k = 1, ke+1
    DO j = 1, je
      DO  i = 2, ie-1

        dzeta_dlam_at_w(i,j,k) = sqrtg_r_w(i,j,k)         &
           &        * 0.5_wp * ( hhl(i+1,j,k) - hhl(i-1,j,k ) ) * eddlon

      END DO
    END DO
  END DO
  ! approximation for lateral boundaries:
  dzeta_dlam_at_w(  1,:,:) = 0.0_wp
  dzeta_dlam_at_w( ie,:,:) = 0.0_wp

  DO k = 1, ke+1
    DO j = 2, je-1
      DO  i = 1, ie

        dzeta_dphi_at_w(i,j,k) = sqrtg_r_w(i,j,k)         &
           &      * 0.5_wp * ( hhl(i,j+1,k) - hhl(i,j-1,k ) ) * eddlat

      END DO
    END DO
  END DO
  ! approximation for lateral boundaries:
  dzeta_dphi_at_w( :, 1,: ) = 0.0_wp
  dzeta_dphi_at_w( :,je,: ) = 0.0_wp

  DO k = 1, ke
    DO j = 1, je
      DO  i = 1, ie-1
        dzeta_dlam_at_u(i,j,k) = 0.5_wp * ( dzeta_dlam(i,j,k) + dzeta_dlam(i+1,j,k) )
      END DO
    END DO
  END DO
  dzeta_dlam_at_u(ie,:,:) = 0.0_wp

  DO k = 1, ke
    DO j = 1, je-1
      DO  i = 1, ie
        dzeta_dphi_at_v(i,j,k) = 0.5_wp * ( dzeta_dphi(i,j,k) + dzeta_dphi(i,j+1,k) )
      END DO
    END DO
  END DO
  dzeta_dphi_at_v(:,je,:) = 0.0_wp

  kzdims(1:24) =                           &
       (/ ke+1, ke+1, ke, ke+1, ke+1,      &
            ke,    0,    0,   0,   0,      &
             0,    0,    0,   0,   0,      &
             0,    0,    0,   0,   0,      &
             0,    0,    0,   0  /)

  ! nicht optimierte Version (icase=0 <--> ldatatypes=.FALSE.)
  CALL exchg_boundaries                                                &
       (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,  &
        kzdims, jstartpar, jendpar, nboundlines, nboundlines, my_cart_neigh,     &
        lperi_x, lperi_y, l2dim, &
        0, .FALSE., ncomm_type, izerror, yzerrmsg,       &
        dzeta_dlam_at_uw(:,:,:),     &
        dzeta_dlam_at_w (:,:,:),     &
        dzeta_dlam_at_u (:,:,:),     &
        dzeta_dphi_at_vw(:,:,:),     &
        dzeta_dphi_at_w (:,:,:),     &
        dzeta_dphi_at_v (:,:,:) )

  IF ( izerror /= 0_iintegers ) THEN
    ierrorstat = 20005
  END IF

END SUBROUTINE calc_dzeta_dxx_at_y

!==============================================================================
!==============================================================================

SUBROUTINE deallocate_dzeta_dxx_at_y

  ! deallocate the following metric coefficient fields:
  ! dzeta_dlam_at_u, dzeta_dlam_at_w, dzeta_dlam_at_uw,
  ! dzeta_dphi_at_v, dzeta_dphi_at_w, dzeta_dphi_at_vw
  ! (this is only needed for 3D-diffusion)

  DEALLOCATE( dzeta_dlam_at_u  )
  DEALLOCATE( dzeta_dlam_at_w  )
  DEALLOCATE( dzeta_dlam_at_uw )

  DEALLOCATE( dzeta_dphi_at_v  )
  DEALLOCATE( dzeta_dphi_at_w  )
  DEALLOCATE( dzeta_dphi_at_vw )

END SUBROUTINE deallocate_dzeta_dxx_at_y

!==============================================================================
!==============================================================================

SUBROUTINE weighting_factors_full2half( hhl,                        &
        wgtfac, wgtfac_u, wgtfac_v, wgtfacq, wgtfacq_u, wgtfacq_v,  &
        ie, je, ke )

!------------------------------------------------------------------------------
!
! Description:
! Weighting factors for linear interpolations from full levels 
! to half levels (for metric terms) are calculated. All fields are 
! set on the full (sub)-domain. In parallel programs the neighboring
! values have to be exchanged after calling this routine.
!
! Input:
!   hhl
!
! Output:
!   wgtfac, wgtfac_u, wgtfac_v, wgtfacq, wgtfacq_u, wgtfacq_v
!
!------------------------------------------------------------------------------

  INTEGER (KIND=iintegers), INTENT(IN) :: ie, je, ke

  REAL (KIND=wp),     INTENT(IN)  :: hhl       (ie,je,ke+1)

  REAL (KIND=wp),     INTENT(OUT) :: wgtfac    (ie,je,ke+1)
  REAL (KIND=wp),     INTENT(OUT) :: wgtfac_u  (ie,je,ke+1)
  REAL (KIND=wp),     INTENT(OUT) :: wgtfac_v  (ie,je,ke+1)
  REAL (KIND=wp),     INTENT(OUT) :: wgtfacq   (ie,je,   3)
  REAL (KIND=wp),     INTENT(OUT) :: wgtfacq_u (ie,je,   3)
  REAL (KIND=wp),     INTENT(OUT) :: wgtfacq_v (ie,je,   3)

  INTEGER            :: i, j, k
  REAL (KIND=wp)     :: z1, z2, z3

!------------------------------------------------------------------------------

  DO  k = 2, ke
    DO  j = 1, je
      !CDIR ON_ADB(hhl)
      DO  i = 1, ie
        wgtfac(i,j,k) = ( hhl(i,j,k-1) - hhl(i,j,k  ) )  &
          &           / ( hhl(i,j,k-1) - hhl(i,j,k+1) )
      ENDDO
    ENDDO
  ENDDO

  DO  j = 1, je
    !CDIR ON_ADB(hhl)
    !CDIR ON_ADB(wgtfac)
    DO  i = 1, ie
      wgtfac(i,j,1)    = ( hhl(i,j,2) - hhl(i,j,1)   )     &
        &              / ( hhl(i,j,3) - hhl(i,j,1)   )
      wgtfac(i,j,ke+1) = ( hhl(i,j,ke  ) - hhl(i,j,ke+1) )  &
        &              / ( hhl(i,j,ke-1) - hhl(i,j,ke+1) )
    ENDDO
  ENDDO

  DO  k = 1, ke+1
    DO  j = 1, je
      !CDIR ON_ADB(wgtfac)
      DO  i = 1, ie-1
        wgtfac_u(i,j,k) = 0.5_wp*(wgtfac(i+1,j,k)+wgtfac(i,j,k))
      ENDDO
    ENDDO
  ENDDO
  wgtfac_u(ie,:,:) = 0.5_wp  ! approximation for lateral boundaries

  DO  k = 1, ke+1
    DO  j = 1, je-1
      !CDIR ON_ADB(wgtfac)
      DO  i = 1, ie
        wgtfac_v(i,j,k) = 0.5_wp*(wgtfac(i,j+1,k)+wgtfac(i,j,k))
      ENDDO
    ENDDO
  ENDDO
  wgtfac_v(:,je,:) = 0.5_wp  ! approximation for lateral boundaries

  ! Factors for quadratic extrapolation of surface pressure 
  ! (used if ldyn_bbc = .FALSE.)
  DO  j = 1, je
    !CDIR ON_ADB(wgtfacq)
    DO  i = 1, ie
      z1 = 0.5_wp * ( hhl(i,j,ke  ) - hhl(i,j,ke+1 ) )
      z2 = 0.5_wp * ( hhl(i,j,ke  ) + hhl(i,j,ke-1) ) - hhl(i,j,ke+1)
      z3 = 0.5_wp * ( hhl(i,j,ke-1) + hhl(i,j,ke-2) ) - hhl(i,j,ke+1)
      wgtfacq(i,j,3) = z1*z2/(z2-z3)/(z1-z3)
      wgtfacq(i,j,2) = (z1-wgtfacq(i,j,3)*(z1-z3))/(z1-z2)
      wgtfacq(i,j,1) = 1._wp - (wgtfacq(i,j,2) + wgtfacq(i,j,3))
    ENDDO
  ENDDO

  DO k = 1, 3
    DO  j = 1, je
      !CDIR ON_ADB(wgtfacq)
      DO  i = 1, ie-1
        wgtfacq_u(i,j,k) = 0.5_wp*(wgtfacq(i+1,j,k)+wgtfacq(i,j,k))
      ENDDO
    ENDDO
  ENDDO
  wgtfacq_u(ie,:,:) = 0.5_wp   ! approximation for lateral boundaries

  DO k = 1, 3
    DO  j = 1, je-1
      !CDIR ON_ADB(wgtfacq)
      DO  i = 1, ie
        wgtfacq_v(i,j,k) = 0.5_wp*(wgtfacq(i,j+1,k)+wgtfacq(i,j,k))
      ENDDO
    ENDDO
  ENDDO
  wgtfacq_v(:,je,:) = 0.5_wp   ! approximation for lateral boundaries

END SUBROUTINE weighting_factors_full2half

!==============================================================================
!==============================================================================

SUBROUTINE weighting_factors_one_sided ( hhl, wgt_one_sided, &
                 ie, je, ke )

!------------------------------------------------------------------------------
!
! Description:
! calculate weighting factors for vertical 2nd order one sided 
! finite difference formula
! (valid for scalar variables like pressure p)
! dp/dz = wgt_one_sided(1) * p(ke)
!       + wgt_one_sided(2) * p(ke-1)
!       + wgt_one_sided(3) * p(ke-2)
!------------------------------------------------------------------------------

  INTEGER (KIND=iintegers), INTENT(IN)  :: ie, je, ke
  REAL    (KIND=wp),        INTENT(IN)  :: hhl          (ie,je,ke+1)

  REAL    (KIND=wp),        INTENT(OUT) :: wgt_one_sided(ie,je,   3)

  REAL    (KIND=wp)     :: delta1, delta2

  INTEGER :: i, j

!------------------------------------------------------------------------------

  DO j=1, je
    DO i=1, ie

      ! distance between main levels at ke-2 and ke:
      delta2 = 0.5_wp * ( ( hhl(i,j,ke-2) + hhl(i,j,ke-1) )   &
        &               - ( hhl(i,j,ke  ) + hhl(i,j,ke+1) ) )

      ! distance between main levels at ke-1 and ke:
      delta1 = 0.5_wp * ( ( hhl(i,j,ke-1) + hhl(i,j,ke)   )   &
        &               - ( hhl(i,j,ke  ) + hhl(i,j,ke+1) ) )

      ! quadratic weights:
      wgt_one_sided(i,j,3) = delta1 / delta2 / ( delta1 - delta2)
      wgt_one_sided(i,j,2) = ( 1.0_wp - wgt_one_sided(i,j,3) * delta2) / delta1
      wgt_one_sided(i,j,1) = -( wgt_one_sided(i,j,2) + wgt_one_sided(i,j,3) )

      ! linear weights:
      !wgt_one_sided(i,j,3) =  0.0_wp
      !wgt_one_sided(i,j,2) =  1.0_wp / delta1
      !wgt_one_sided(i,j,1) = -1.0_wp / delta1

    END DO
  END DO

  ! i=4
  ! j=4
  ! WRITE(*,'(A,3I4,3E14.5)') "wgt ", my_cart_id, i, j,          &
  !   & wgt_one_sided(i, j,1), wgt_one_sided(i, j,2), wgt_one_sided(i, j,3)

END SUBROUTINE weighting_factors_one_sided

!==============================================================================
!==============================================================================

SUBROUTINE final_grid_metrics (ierrorstat)

!------------------------------------------------------------------------------
!
! Description:
!   Deallocate all arrays
!
!------------------------------------------------------------------------------

  INTEGER (KIND=iintegers), INTENT(OUT) :: ierrorstat

!------------------------------------------------------------------------------

  ierrorstat = 0

  DEALLOCATE (sqrtg_r_s, sqrtg_r_u , sqrtg_r_v , sqrtg_r_w ,               &
              STAT=ierrorstat)

  IF ( l_dzeta_d_needed ) THEN
    DEALLOCATE (dzeta_dlam, dzeta_dphi, STAT=ierrorstat)
  END IF

  DEALLOCATE (wgtfac   , wgtfac_u , wgtfac_v , wgtfacq  ,                  &
              wgtfacq_u, wgtfacq_v, wgt_one_sided,   STAT=ierrorstat)

END SUBROUTINE final_grid_metrics

!==============================================================================

END MODULE grid_metrics_utilities
