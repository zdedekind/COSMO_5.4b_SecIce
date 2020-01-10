!+ Source module for organizing the spectral nudging
!-------------------------------------------------------------------------------

MODULE src_spectral_nudging

!-------------------------------------------------------------------------------
!
! Description:
! This module contains als subroutines for nudging the large scales of the 
! driving boundary data. 
! 
! Current Code Owner: GKSS, Burkhardt Rockel
!  phone:   +49 4152 87 2008
!  fax:     +49 4152 87 2020
!  email:   rockel@gkss.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 3.18       2006/03/03 Burkhardt Rockel (after H. Langenberg, F. Feser, 
!  Initial release                        and B. Mueller)
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_23        2012/05/10 Burkhardt Rockel
!  several changes for better parallelisation
!  only 4D prognostic fields allowed
! V4_28        2013/07/12 Ulrich Schaettler
!  Use variables for vertical grid from module vgrid_refatm_utils
! V4_29        2013/10/04 Ulrich Schaettler
!  For the COSMO-Model only use vcoord from vgrid_refatm_utils
! V5_1         2014-11-28 Oliver Fuhrer
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

USE data_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

USE data_parallel,      ONLY :  &
    lasync_io,      & ! if .TRUE.: the model runs with extra PEs for
                      ! asynchronous IO
    my_cart_id,     & ! rank of this subdomain in the cartesian communicator
    my_cart_neigh,  & ! neighbors
    isubpos,        & ! positions of the subdomains in the total domain. Given
                      ! are the i- and the j-indices of the lower left and the
                      ! upper right grid point in the order
                      !                  i_ll, j_ll, i_ur, j_ur.
                      ! Only the interior of the domains are considered, not
                      ! the boundary lines.
    num_compute,    & ! number of compute PEs
    nboundlines,    & ! number of boundary lines of the domain for which
                      ! no forecast is computed = overlapping boundary
                      ! lines of the subdomains
    icomm_cart,     & ! communicator that belongs to the cartesian grid 
    imp_reals,      & ! determines the correct REAL type used in the model
                      ! for MPI
    imp_grib,       & ! determines the REAL type for the GRIB library
    imp_integers,   & ! determines the correct INTEGER type used in the model
                      ! for MPI
    imp_logical       ! determines the correct LOGICAL   type used in the
                      ! model for MPI

!------------------------------------------------------------------------------

USE data_io,          ONLY:                              &
    yvarbd,         & ! list of variables for Output
    nyvar_b,        & ! number of variables for Output
    num_gribtabs,   & ! number of GRIB tables used in LM variable table
    var               ! array for LM variable table

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY : &
    ie_tot,       & ! total number of grid points in zonal direction
    je_tot,       & ! total number of grid points in meridional direction
    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction
    ke1             ! ke+1

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
    nlastbound,   & ! time step of the last boundary update
    nincbound,    & ! time step increment of boundary update
    nbd1,         & ! indices for permutation of the 
    nbd2,         & ! two boundary time levels
    nnew,         & ! corresponds to ntstep + 1
    yvarsn,       & ! list of fields for spectral nudging
    nyvar_sn,     & ! number of spectral variables
    list_sn,      & ! variable list for spectral nudging
    isc_sn,       & ! spectral nudging in i-direction
    jsc_sn,       & ! spectral nudging in j-direction
    pp_sn,        & ! lowest pressure level for spectral nudging
    alpha_sn        ! amplification factor for spectral nudging (0.<= alpha_sn<=1.)

!------------------------------------------------------------------------------

USE environment , ONLY :   &
    model_abort     ! aborts the program in case of errors
    
USE parallel_utilities, ONLY : &
    gather_field,     & ! gathers data from several processors to processor 0
    distribute_field    ! distributes data from a processor to all other processors

USE vgrid_refatm_utils, ONLY : &
    vcoord              ! parameters for the vertical grid

!==============================================================================

  IMPLICIT NONE

! include statements
INCLUDE "mpif.h"

!==============================================================================

  INTEGER :: &
    k_sn            ! lowest full level for spectral nudging 

  REAL (KIND=wp),     TARGET, ALLOCATABLE ::  &
      var_grid (:,:,:),        & ! grid array of actual variable
      var_grid_bd1(:,:,:),     & ! grid array of boundary variable
      var_grid_bd2(:,:,:),     & ! grid array of boundary variable
      var_spec (:,:),          & ! spectral array of actual variable
      var_spec_bd(:,:)           ! spectral array of boundary variable
     
  REAL (KIND=wp),     TARGET, ALLOCATABLE ::  &
      aa(:,:),             & !\   Fourier coefficients of the actual array
      ab(:,:),             & ! \  aa = cos cos        ab = cos sin
      ba(:,:),             & ! /  ba = sin cos        bb = sin sin
      bb(:,:),             & !/
      aa_bd(:,:),          & !\   Fourier coefficients of the boundary array
      ab_bd(:,:),          & ! \  aa = cos cos        ab = cos sin
      ba_bd(:,:),          & ! /  ba = sin cos        bb = sin sin
      bb_bd(:,:),          & !/
      aa_n(:,:),           & !\   Fourier coefficients of the nudged array
      ab_n(:,:),           & ! \  aa = cos cos        ab = cos sin
      ba_n(:,:),           & ! /  ba = sin cos        bb = sin sin
      bb_n(:,:)              !/

  REAL (KIND=wp),     TARGET, ALLOCATABLE ::  &
       wgtvd_sn(:)           ! vertical dependent weights for the nudging

  REAL (KIND=wp),     TARGET, ALLOCATABLE ::   &
       asem(:,:),         & !  
       bsem(:,:),         & !
       asem_n(:,:),       & !
       bsem_n(:,:),       & !
       suma(:,:),         & !
       sumb(:,:),         & !
       sum(:,:),          & !
       suma_n(:,:),       & !
       sumb_n(:,:),       & !
       sum_n(:,:)           !

  REAL (KIND=wp),     TARGET, ALLOCATABLE ::       &
      SinArg(:,:,:),      & ! working array for Sinus function in spectral decomposition
      COSArg(:,:,:),      & ! working array for Cosinus function in spectral decomposition
      SinArgDIV(:,:,:),   & ! SinArg divided by ie_tot and je_tot
      CosArgDIV(:,:,:)      ! CosArg divided by ie_tot and je_tot

!==============================================================================
! Module Procedures
!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure for intializing spectral nudging
!------------------------------------------------------------------------------

SUBROUTINE init_spectral_nudging

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure initializes values necessary for computing
!   the spectral nudging.
!
! Method:
!
!------------------------------------------------------------------------------

! Local scalars:
! -------------

  REAL (KIND=wp)           ::   &   
       twopi,      & ! 2. * Pi
       arg           ! help variable

  INTEGER (KIND=iintegers) ::   &
    k,             & ! loop variable
    ms,            & ! loop variable
    ns,            & ! loop variable
    ir,            & ! loop variable
    jr,            & ! loop variable
    istat,         & ! status variable
    isc_jsc_max,   & ! maximum of isc_sn and jsc_sn
    ie_je_max,     & ! maximum of ie_tot and je_tot
    iz1, iz2, iz3, n, i

  LOGICAL                  ::   &
    lzcontained
    
!----------- End of header ----------------------------------------------------
 
!------------------------------------------------------------------------------
! Begin Subroutine init_spectral_nudging
!------------------------------------------------------------------------------

!----------------------------------------------------------------------------
! Section 1: Initialize the description list list_sn
!----------------------------------------------------------------------------

  ! Check, whether all variables in yvar_sn are also contained in yvar_bd
    ! nudging variables must also be boundary variables
    DO i = 1, nyvar_sn
      lzcontained = .FALSE.
      bd_loop: DO n = 1, nyvar_b
        IF ( TRIM(yvarsn(i)) == TRIM(yvarbd(n)) ) THEN
          lzcontained = .TRUE.
          EXIT bd_loop
        ENDIF
      ENDDO bd_loop
      IF (.NOT. lzcontained) THEN
        IF (my_cart_id == 0) THEN
          PRINT *,' ERROR  *** yvarsn = ', TRIM(yvarsn(i)),             &
                  ' is no boundary variable ***'
        ENDIF
        CALL model_abort(my_cart_id, 2029,                    &
               'Wrong variable in list for spectral nudging', &
               'init_spectral_nudging')
      ENDIF
    ENDDO

  ! then initialize description list
  ALLOCATE(list_sn(nyvar_sn), STAT=istat)

  loop_list_snvar: DO n = 1, nyvar_sn
    list_sn(n)%name = yvarsn(n)
    DO iz3 = 1, num_gribtabs
      DO iz2 = 0, 255
        DO iz1 = 1,4
          IF (var(iz1,iz2,iz3)%name(1:LEN_TRIM(var(iz1,iz2,iz3)%name)) ==  &
                   list_sn (n)%name(1:LEN_TRIM(     list_sn (n)%name)) ) THEN
            ! Set location in variable table
            list_sn(n)%iloc1 = iz1
            list_sn(n)%iloc2 = iz2
            list_sn(n)%iloc3 = iz3

            ! Set dimension of this variable (1, ke or ke+1)
            SELECT CASE (var(iz1,iz2,iz3)%rank)
            CASE (4)
              list_sn(n)%idimvert = UBOUND (var(iz1,iz2,iz3)%p4,3)
            CASE (3)
              SELECT CASE (var(iz1,iz2,iz3)%name)
              CASE ('PS        ','T_G       ','QV_S      ','W_SNOW    ',  &
                    'T_S       ','T_M       ','W_G1      ','W_G2      ',  &
                    'W_G3      ','W_I       ','T_SNOW    ','RHO_SNOW  ',  &
                    'H_SNOW    ','H_ICE     ','T_ICE     ','T_MNW_LK  ',  &
                    'T_WML_LK  ','T_BOT_LK  ','T_B1_LK   ','C_T_LK    ',  &
                    'H_ML_LK   ','H_B1_LK   ')
                ! these are 2D variables with rank=3 because of time dependency
                list_sn(n)%idimvert = 1
              CASE DEFAULT
                ! these are real 3D variables
                list_sn(n)%idimvert = UBOUND (var(iz1,iz2,iz3)%p3,3)
              END SELECT
            CASE (2)
              list_sn(n)%idimvert = 1
            END SELECT

            CYCLE loop_list_snvar
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO loop_list_snvar

!----------------------------------------------------------------------------
! Section 2: allocate memory and calculate help variables
!----------------------------------------------------------------------------

  twopi = 8.0_wp*ATAN(1.0_wp) 
  isc_jsc_max = MAX(isc_sn, jsc_sn)
  ie_je_max   = MAX(ie_tot, je_tot)

  ! Allocate data fields fo spectral nudging
  ALLOCATE (var_spec (ie_tot, je_tot),           STAT=istat)
  ALLOCATE (var_spec_bd(ie_tot, je_tot),         STAT=istat)

  ALLOCATE (aa(isc_sn, jsc_sn),                  STAT=istat)
  ALLOCATE (ab(isc_sn, jsc_sn),                  STAT=istat)
  ALLOCATE (ba(isc_sn, jsc_sn),                  STAT=istat)
  ALLOCATE (bb(isc_sn, jsc_sn),                  STAT=istat)
  ALLOCATE (aa_bd(isc_sn, jsc_sn),               STAT=istat)
  ALLOCATE (ab_bd(isc_sn, jsc_sn),               STAT=istat)
  ALLOCATE (ba_bd(isc_sn, jsc_sn),               STAT=istat)
  ALLOCATE (bb_bd(isc_sn, jsc_sn),               STAT=istat)
  ALLOCATE (aa_n(isc_sn, jsc_sn),                STAT=istat)
  ALLOCATE (ab_n(isc_sn, jsc_sn),                STAT=istat)
  ALLOCATE (ba_n(isc_sn, jsc_sn),                STAT=istat)
  ALLOCATE (bb_n(isc_sn, jsc_sn),                STAT=istat)

  ALLOCATE (wgtvd_sn(ke1),                       STAT=istat)

  ALLOCATE (asem(isc_sn, je_tot),                STAT=istat)
  ALLOCATE (bsem(isc_sn, je_tot),                STAT=istat)
  ALLOCATE (asem_n(isc_sn, je_tot),              STAT=istat)
  ALLOCATE (bsem_n(isc_sn, je_tot),              STAT=istat)
  ALLOCATE (suma(isc_sn, je_tot),                STAT=istat)
  ALLOCATE (sumb(isc_sn, je_tot),                STAT=istat)
  ALLOCATE (sum(ie_tot, je_tot),                 STAT=istat)
  ALLOCATE (suma_n(isc_sn, je_tot),              STAT=istat)
  ALLOCATE (sumb_n(isc_sn, je_tot),              STAT=istat)
  ALLOCATE (sum_n(ie_tot, je_tot),               STAT=istat)

  ALLOCATE (SinArg(isc_jsc_max, ie_je_max, 2),      STAT=istat)
  ALLOCATE (COSArg(isc_jsc_max, ie_je_max, 2),      STAT=istat)
  ALLOCATE (SinArgDIV(isc_jsc_max, ie_je_max, 2),   STAT=istat)
  ALLOCATE (CosArgDIV(isc_jsc_max, ie_je_max, 2),   STAT=istat)

  IF(istat /= 0) THEN
    CALL model_abort(my_cart_id, 2029,                                    &
                 'Error allocating fields for spectral nudging', 'init_spectral_nudging')
  ENDIF

  ! calculate some Cosinus and Sinus dependent functions
  DO ms = 1,isc_sn
    DO ir = 1,ie_tot
      arg =twopi*REAL(ms-1,wp)*REAL(ir,wp)/REAL(ie_tot,wp)
      SinArg(ms,ir,1)=SIN(arg)
      CosArg(ms,ir,1)=COS(arg)
      SinArgDiv(ms,ir,1)=SinArg(ms,ir,1)/REAL(ie_tot,wp)
      CosArgDiv(ms,ir,1)=CosArg(ms,ir,1)/REAL(ie_tot,wp)
    ENDDO
  ENDDO

  DO ns = 1,jsc_sn
    DO jr = 1,je_tot
      arg =twopi*REAL(ns-1,wp)*REAL(jr,wp)/REAL(je_tot,wp)
      SinArg(ns,jr,2)=SIN(arg)
      CosArg(ns,jr,2)=COS(arg)
      SinArgDiv(ns,jr,2)=SinArg(ns,jr,2)/REAL(je_tot,wp)
      CosArgDiv(ns,jr,2)=CosArg(ns,jr,2)/REAL(je_tot,wp)
    ENDDO
  ENDDO

  ! Determine the lowest vertical level for spectral nudging
  k_sn = ke
  DO k = 1 , ke
! for reference pressure   p0sl = 10E5_wp
    IF (      (vcoord%sigm_coord(k)  <= pp_sn/1000.0_wp)   &
        .AND. (vcoord%sigm_coord(k+1) > pp_sn/1000.0_wp)) THEN
      k_sn = k
      EXIT
    ENDIF
  ENDDO

  ! Amplification for vertical nudging
  DO k = 1, k_sn
    wgtvd_sn(k) = ( 1.0_wp -                                               &
     MIN ( ((vcoord%sigm_coord(k)+vcoord%sigm_coord(k+1))*0.5E5_wp),       &
                     pp_sn*100.0_wp )/                                     &
                                           (pp_sn*100.0_wp) )**2*alpha_sn
  ENDDO

  ALLOCATE (var_grid    (ie_tot, je_tot, k_sn),   STAT=istat)
  ALLOCATE (var_grid_bd1(ie_tot, je_tot, k_sn),   STAT=istat)
  ALLOCATE (var_grid_bd2(ie_tot, je_tot, k_sn),   STAT=istat)

!------------------------------------------------------------------------------
! End of module procedure init_spectral_nudging
!------------------------------------------------------------------------------

END SUBROUTINE init_spectral_nudging

!==============================================================================
!==============================================================================
!+ Module procedure for spectral nudging
!------------------------------------------------------------------------------

SUBROUTINE spectral_nudging 

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure performs the spectral nudging of the boundary data
!
! Method:
!
!------------------------------------------------------------------------------

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
      i,                  & ! loop variable
      ii,                 & ! loop variable
      iik,                & ! loop variable
      j,                  & ! loop variable
      k,                  & ! loop variable
      kweight,            & ! height index for weighting coefficient
      kwgt(k_sn*nyvar_sn),& ! height index for weighting coefficient
      i1,                 & ! location in variable table
      i2,                 & ! location in variable table
      i3,                 & ! location in variable table
      kend_sn,            & ! lowest layer for spectral nudging
      ierror                ! error status
      
  REAL (KIND=wp)     ::   &
      zdtibd,             & ! weights the two boundary data arrays 
                            ! linearily in time 
      zemzdr                ! 1. - zdtibd

  CHARACTER (LEN= 25)        ::  &
    yroutine   ! name of this routine for error handling

!----------- End of header ----------------------------------------------------
 
!------------------------------------------------------------------------------
! Begin Subroutine spectral_nudging
!------------------------------------------------------------------------------

  ! Preset values
  yroutine = 'spectral_nudging'
  ierror = 0
  aa_n = 0.0_wp
  ab_n = 0.0_wp
  ba_n = 0.0_wp
  bb_n = 0.0_wp

  ! Coefficients for boundary interpolation linear in time
  zdtibd = REAL (ntstep+1-nlastbound, wp) / REAL (nincbound, wp)
  zemzdr = 1.0_wp - zdtibd

  kend_sn = k_sn

  IF (k_sn*nyvar_sn <= num_compute) THEN
  
  ! loop over spectral variables
    DO ii=1,nyvar_sn

    ! Location of this variable in the variable table
    i1 = list_sn(ii)%iloc1
    i2 = list_sn(ii)%iloc2
    i3 = list_sn(ii)%iloc3
 
    ! Write the actual data field and the boundary fields on extra fields

        DO  k = 1, kend_sn
        iik = k + (ii-1)*kend_sn
        kwgt(iik) = k
          CALL gather_field (var(i1,i2,i3)%p4(:,:,k,nnew), ie, je,        &
                           var_grid(:,:,1), ie_tot, je_tot, iik-1, ierror)
          CALL gather_field (var(i1,i2,i3)%p4_bd(:,:,k,nbd1), ie, je,     &
                           var_grid_bd1(:,:,1), ie_tot, je_tot, iik-1, ierror)
          CALL gather_field (var(i1,i2,i3)%p4_bd(:,:,k,nbd2), ie, je,     &
                           var_grid_bd2(:,:,1), ie_tot, je_tot, iik-1, ierror)
        ENDDO

    ENDDO
   
    IF (my_cart_id <  kend_sn*nyvar_sn) THEN

      kweight = kwgt(my_cart_id+1)
      ! Loop over all vertical levels that will be nudged
!      DO  k = 1, kend_sn

        ! Write actual grid field and boundary field (interpolated in time) 
        ! on 2D-fields
        DO j = 1, je_tot
          DO i = 1, ie_tot
            var_spec(i,j)    = var_grid(i,j,1)
            var_spec_bd(i,j) = zemzdr*var_grid_bd1(i,j,1) +      &
                               zdtibd*var_grid_bd2(i,j,1)
          ENDDO
        ENDDO
 
 
        ! Calculate Fourier coefficients for actual data for large scales
        CALL spectral_decomposition (var_spec(1:ie_tot,1:je_tot), aa(1:isc_sn, 1:jsc_sn), &
           ab(1:isc_sn, 1:jsc_sn), ba(1:isc_sn, 1:jsc_sn), bb(1:isc_sn, 1:jsc_sn))

        ! Calculate Fourier coefficients for boundary data for large scales
        CALL spectral_decomposition (var_spec_bd(1:ie_tot,1:je_tot), aa_bd(1:isc_sn, 1:jsc_sn), &
        ab_bd(1:isc_sn, 1:jsc_sn), ba_bd(1:isc_sn, 1:jsc_sn), bb_bd(1:isc_sn, 1:jsc_sn))

        ! Nudging of large scales only: j up to 3/2PI, i up to 5/2PI
        DO  i = 1,isc_sn
          DO  j = 1,jsc_sn

!            aa_n(i,j) = aa(i,j) - (wgtvd_sn(k))*( aa(i,j) - aa_bd(i,j) )
!            ab_n(i,j) = ab(i,j) - (wgtvd_sn(k))*( ab(i,j) - ab_bd(i,j) )
!            ba_n(i,j) = ba(i,j) - (wgtvd_sn(k))*( ba(i,j) - ba_bd(i,j) )
!            bb_n(i,j) = bb(i,j) - (wgtvd_sn(k))*( bb(i,j) - bb_bd(i,j) )
!            aa_n(i,j) = aa(i,j) - (wgtvd_sn(my_cart_id+1))*( aa(i,j) - aa_bd(i,j) )
!            ab_n(i,j) = ab(i,j) - (wgtvd_sn(my_cart_id+1))*( ab(i,j) - ab_bd(i,j) )
!            ba_n(i,j) = ba(i,j) - (wgtvd_sn(my_cart_id+1))*( ba(i,j) - ba_bd(i,j) )
!            bb_n(i,j) = bb(i,j) - (wgtvd_sn(my_cart_id+1))*( bb(i,j) - bb_bd(i,j) )
            aa_n(i,j) = aa(i,j) - (wgtvd_sn(kweight))*( aa(i,j) - aa_bd(i,j) )
            ab_n(i,j) = ab(i,j) - (wgtvd_sn(kweight))*( ab(i,j) - ab_bd(i,j) )
            ba_n(i,j) = ba(i,j) - (wgtvd_sn(kweight))*( ba(i,j) - ba_bd(i,j) )
            bb_n(i,j) = bb(i,j) - (wgtvd_sn(kweight))*( bb(i,j) - bb_bd(i,j) )

          ENDDO
        ENDDO

       ! Composition of spectral coefficients back to grid field
        CALL spectral_composition                                          &
              (aa(1:isc_sn, 1:jsc_sn), ab(1:isc_sn, 1:jsc_sn), ba(1:isc_sn, 1:jsc_sn), &
              bb(1:isc_sn, 1:jsc_sn), aa_n(1:isc_sn, 1:jsc_sn), ab_n(1:isc_sn, 1:jsc_sn), &
              ba_n(1:isc_sn, 1:jsc_sn), bb_n(1:isc_sn, 1:jsc_sn), var_spec(1:ie_tot,1:je_tot))

        ! Write nudged 2D-field to 3D-field
        DO j = 1, je_tot
          DO i = 1, ie_tot
            var_grid(i,j,1) = var_spec (i,j)
          ENDDO
        ENDDO 

      ENDIF

    DO ii=1,nyvar_sn

    ! Location of this variable in the variable table
      i1 = list_sn(ii)%iloc1
      i2 = list_sn(ii)%iloc2
      i3 = list_sn(ii)%iloc3
    
        DO  k = 1, kend_sn
        iik = k + (ii-1)*kend_sn
        CALL distribute_field (var_grid(1:ie_tot,1:je_tot,1), ie_tot, je_tot,         &
                               var(i1,i2,i3)%p4(:,:,k,nnew), ie, je, iik-1, ierror)
        ENDDO

    ENDDO



  ELSE IF (k_sn <= num_compute) THEN
!  IF (k_sn <= num_compute) THEN

  ! loop over spectral variables
  DO ii=1,nyvar_sn

    ! Location of this variable in the variable table
    i1 = list_sn(ii)%iloc1
    i2 = list_sn(ii)%iloc2
    i3 = list_sn(ii)%iloc3
 
    ! Write the actual data field and the boundary fields on extra fields

         DO  k = 1, kend_sn
          CALL gather_field (var(i1,i2,i3)%p4(:,:,k,nnew), ie, je,        &
                             var_grid(:,:,1), ie_tot, je_tot, k-1, ierror)
          CALL gather_field (var(i1,i2,i3)%p4_bd(:,:,k,nbd1), ie, je,     &
                             var_grid_bd1(:,:,1), ie_tot, je_tot, k-1, ierror)
          CALL gather_field (var(i1,i2,i3)%p4_bd(:,:,k,nbd2), ie, je,     &
                             var_grid_bd2(:,:,1), ie_tot, je_tot, k-1, ierror)
        ENDDO

   

    IF (my_cart_id <  kend_sn) THEN

      ! Loop over all vertical levels that will be nudged
!      DO  k = 1, kend_sn

        ! Write actual grid field and boundary field (interpolated in time) 
        ! on 2D-fields
        DO j = 1, je_tot
          DO i = 1, ie_tot
            var_spec(i,j)    = var_grid(i,j,1)
            var_spec_bd(i,j) = zemzdr*var_grid_bd1(i,j,1) +      &
                               zdtibd*var_grid_bd2(i,j,1)
          ENDDO
        ENDDO
 
 
        ! Calculate Fourier coefficients for actual data for large scales
        CALL spectral_decomposition (var_spec(1:ie_tot,1:je_tot), aa(1:isc_sn, 1:jsc_sn), &
           ab(1:isc_sn, 1:jsc_sn), ba(1:isc_sn, 1:jsc_sn), bb(1:isc_sn, 1:jsc_sn))

        ! Calculate Fourier coefficients for boundary data for large scales
        CALL spectral_decomposition (var_spec_bd(1:ie_tot,1:je_tot), aa_bd(1:isc_sn, 1:jsc_sn), &
        ab_bd(1:isc_sn, 1:jsc_sn), ba_bd(1:isc_sn, 1:jsc_sn), bb_bd(1:isc_sn, 1:jsc_sn))

        ! Nudging of large scales only: j up to 3/2PI, i up to 5/2PI
        DO  i = 1,isc_sn
          DO  j = 1,jsc_sn

            aa_n(i,j) = aa(i,j) - (wgtvd_sn(my_cart_id+1))*( aa(i,j) - aa_bd(i,j) )
            ab_n(i,j) = ab(i,j) - (wgtvd_sn(my_cart_id+1))*( ab(i,j) - ab_bd(i,j) )
            ba_n(i,j) = ba(i,j) - (wgtvd_sn(my_cart_id+1))*( ba(i,j) - ba_bd(i,j) )
            bb_n(i,j) = bb(i,j) - (wgtvd_sn(my_cart_id+1))*( bb(i,j) - bb_bd(i,j) )

          ENDDO
        ENDDO

       ! Composition of spectral coefficients back to grid field
        CALL spectral_composition                                          &
              (aa(1:isc_sn, 1:jsc_sn), ab(1:isc_sn, 1:jsc_sn), ba(1:isc_sn, 1:jsc_sn), &
              bb(1:isc_sn, 1:jsc_sn), aa_n(1:isc_sn, 1:jsc_sn), ab_n(1:isc_sn, 1:jsc_sn), &
              ba_n(1:isc_sn, 1:jsc_sn), bb_n(1:isc_sn, 1:jsc_sn), var_spec(1:ie_tot,1:je_tot))

        ! Write nudged 2D-field to 3D-field
        DO j = 1, je_tot
          DO i = 1, ie_tot
            var_grid(i,j,1) = var_spec (i,j)
          ENDDO
        ENDDO 

!      ENDDO    ! End of loop over vertical levels
      
    ENDIF

    ! Write the calculated extra field to the main field

      DO  k = 1, kend_sn
        CALL distribute_field (var_grid(1:ie_tot,1:je_tot,1), ie_tot, je_tot,         &
                               var(i1,i2,i3)%p4(:,:,k,nnew), ie, je, k-1, ierror)
      ENDDO

    ENDDO  ! variable_loop


 ELSE IF (nyvar_sn <= num_compute) THEN  !+++++++++++++++++++++++++++++++++++

  ! loop over spectral variables
    DO ii=1,nyvar_sn

    ! Location of this variable in the variable table
      i1 = list_sn(ii)%iloc1
      i2 = list_sn(ii)%iloc2
      i3 = list_sn(ii)%iloc3
 
    ! Write the actual data field and the boundary fields on extra fields
        DO  k = 1, kend_sn
          CALL gather_field (var(i1,i2,i3)%p4(:,:,k,nnew), ie, je,        &
                             var_grid(:,:,k), ie_tot, je_tot, ii-1, ierror)
          CALL gather_field (var(i1,i2,i3)%p4_bd(:,:,k,nbd1), ie, je,     &
                             var_grid_bd1(:,:,k), ie_tot, je_tot, ii-1, ierror)
          CALL gather_field (var(i1,i2,i3)%p4_bd(:,:,k,nbd2), ie, je,     &
                             var_grid_bd2(:,:,k), ie_tot, je_tot, ii-1, ierror)
        ENDDO

    ENDDO  ! variable_loop

    IF (my_cart_id < nyvar_sn) THEN

      ! Loop over all vertical levels that will be nudged
      DO  k = 1, kend_sn

        ! Write actual grid field and boundary field (interpolated in time) 
        ! on 2D-fields
        DO j = 1, je_tot
          DO i = 1, ie_tot
            var_spec(i,j)    = var_grid(i,j,k)
            var_spec_bd(i,j) = zemzdr*var_grid_bd1(i,j,k) +      &
                               zdtibd*var_grid_bd2(i,j,k)
          ENDDO
        ENDDO
    
        ! Calculate Fourier coefficients for actual data for large scales
        CALL spectral_decomposition (var_spec(1:ie_tot,1:je_tot), aa(1:isc_sn, 1:jsc_sn), &
           ab(1:isc_sn, 1:jsc_sn), ba(1:isc_sn, 1:jsc_sn), bb(1:isc_sn, 1:jsc_sn))

        ! Calculate Fourier coefficients for boundary data for large scales
        CALL spectral_decomposition (var_spec_bd(1:ie_tot,1:je_tot), aa_bd(1:isc_sn, 1:jsc_sn), &
        ab_bd(1:isc_sn, 1:jsc_sn), ba_bd(1:isc_sn, 1:jsc_sn), bb_bd(1:isc_sn, 1:jsc_sn))

        ! Nudging of large scales only: j up to 3/2PI, i up to 5/2PI
        DO  i = 1,isc_sn
          DO  j = 1,jsc_sn

            aa_n(i,j) = aa(i,j) - (wgtvd_sn(k))*( aa(i,j) - aa_bd(i,j) )
            ab_n(i,j) = ab(i,j) - (wgtvd_sn(k))*( ab(i,j) - ab_bd(i,j) )
            ba_n(i,j) = ba(i,j) - (wgtvd_sn(k))*( ba(i,j) - ba_bd(i,j) )
            bb_n(i,j) = bb(i,j) - (wgtvd_sn(k))*( bb(i,j) - bb_bd(i,j) )

          ENDDO
        ENDDO

        ! Composition of spectral coefficients back to grid field
        CALL spectral_composition                                          &
              (aa(1:isc_sn, 1:jsc_sn), ab(1:isc_sn, 1:jsc_sn), ba(1:isc_sn, 1:jsc_sn), &
              bb(1:isc_sn, 1:jsc_sn), aa_n(1:isc_sn, 1:jsc_sn), ab_n(1:isc_sn, 1:jsc_sn), &
              ba_n(1:isc_sn, 1:jsc_sn), bb_n(1:isc_sn, 1:jsc_sn), var_spec(1:ie_tot,1:je_tot))

        ! Write nudged 2D-field to 3D-field
        DO j = 1, je_tot
          DO i = 1, ie_tot
            var_grid(i,j,k) = var_spec (i,j)
          ENDDO
        ENDDO 

      ENDDO    ! End of loop over vertical levels

    ENDIF

    ! Write the calculated extra field to the main field
   ! loop over spectral variables
    DO ii=1,nyvar_sn

    ! Location of this variable in the variable table
      i1 = list_sn(ii)%iloc1
      i2 = list_sn(ii)%iloc2
      i3 = list_sn(ii)%iloc3
        DO  k = 1, kend_sn
          CALL distribute_field (var_grid(1:ie_tot,1:je_tot,k), ie_tot, je_tot,         &
                                 var(i1,i2,i3)%p4(:,:,k,nnew), ie, je, ii-1, ierror)
        ENDDO

    ENDDO  ! variable_loop


  ELSE IF (num_compute > 1) THEN   ! +++++++++++++++++++++++++++++++++++++
  
    ! loop over spectral variables
    DO ii=1,nyvar_sn

    ! Location of this variable in the variable table
      i1 = list_sn(ii)%iloc1
      i2 = list_sn(ii)%iloc2
      i3 = list_sn(ii)%iloc3
 
    ! Write the actual data field and the boundary fields on extra fields
      kend_sn = k_sn
      DO  k = 1, kend_sn
        CALL gather_field (var(i1,i2,i3)%p4(:,:,k,nnew), ie, je,        &
                           var_grid(:,:,k), ie_tot, je_tot, 0, ierror)
        CALL gather_field (var(i1,i2,i3)%p4_bd(:,:,k,nbd1), ie, je,     &
                           var_grid_bd1(:,:,k), ie_tot, je_tot, 0, ierror)
        CALL gather_field (var(i1,i2,i3)%p4_bd(:,:,k,nbd2), ie, je,     &
                           var_grid_bd2(:,:,k), ie_tot, je_tot, 0, ierror)
      ENDDO

      IF (my_cart_id == 0) THEN

      ! Loop over all vertical levels that will be nudged
        DO  k = 1, kend_sn

        ! Write actual grid field and boundary field (interpolated in time) 
        ! on 2D-fields
          DO j = 1, je_tot
            DO i = 1, ie_tot
              var_spec(i,j)    = var_grid(i,j,k)
              var_spec_bd(i,j) = zemzdr*var_grid_bd1(i,j,k) +      &
                                 zdtibd*var_grid_bd2(i,j,k)
            ENDDO
          ENDDO
    
        ! Calculate Fourier coefficients for actual data for large scales
        CALL spectral_decomposition (var_spec(1:ie_tot,1:je_tot), aa(1:isc_sn, 1:jsc_sn), &
           ab(1:isc_sn, 1:jsc_sn), ba(1:isc_sn, 1:jsc_sn), bb(1:isc_sn, 1:jsc_sn))

        ! Calculate Fourier coefficients for boundary data for large scales
        CALL spectral_decomposition (var_spec_bd(1:ie_tot,1:je_tot), aa_bd(1:isc_sn, 1:jsc_sn), &
        ab_bd(1:isc_sn, 1:jsc_sn), ba_bd(1:isc_sn, 1:jsc_sn), bb_bd(1:isc_sn, 1:jsc_sn))

        ! Nudging of large scales only: j up to 3/2PI, i up to 5/2PI
          DO  i = 1,isc_sn
            DO  j = 1,jsc_sn

              aa_n(i,j) = aa(i,j) - (wgtvd_sn(k))*( aa(i,j) - aa_bd(i,j) )
              ab_n(i,j) = ab(i,j) - (wgtvd_sn(k))*( ab(i,j) - ab_bd(i,j) )
              ba_n(i,j) = ba(i,j) - (wgtvd_sn(k))*( ba(i,j) - ba_bd(i,j) )
              bb_n(i,j) = bb(i,j) - (wgtvd_sn(k))*( bb(i,j) - bb_bd(i,j) )

            ENDDO
          ENDDO

       ! Composition of spectral coefficients back to grid field
        CALL spectral_composition                                          &
              (aa(1:isc_sn, 1:jsc_sn), ab(1:isc_sn, 1:jsc_sn), ba(1:isc_sn, 1:jsc_sn), &
              bb(1:isc_sn, 1:jsc_sn), aa_n(1:isc_sn, 1:jsc_sn), ab_n(1:isc_sn, 1:jsc_sn), &
              ba_n(1:isc_sn, 1:jsc_sn), bb_n(1:isc_sn, 1:jsc_sn), var_spec(1:ie_tot,1:je_tot))

        ! Write nudged 2D-field to 3D-field
          DO j = 1, je_tot
            DO i = 1, ie_tot
              var_grid(i,j,k) = var_spec (i,j)
            ENDDO
          ENDDO 

        ENDDO    ! End of loop over vertical levels

      ENDIF

    ! Write the calculated extra field to the main field
        DO  k = 1, kend_sn
          CALL distribute_field (var_grid(:,:,k), ie_tot, je_tot,         &
                                 var(i1,i2,i3)%p4(:,:,k,nnew), ie, je, 0, ierror)
        ENDDO

    ENDDO  ! variable_loop

  
  ELSE  ! seriell +++++++++++++++++++++++++++++++++++++++
  
    ! loop over spectral variables
    DO ii=1,nyvar_sn

    ! Location of this variable in the variable table
      i1 = list_sn(ii)%iloc1
      i2 = list_sn(ii)%iloc2
      i3 = list_sn(ii)%iloc3
 
    ! Write the actual data field and the boundary fields on extra fields
    
      kend_sn = k_sn
      DO  k = 1, kend_sn
        var_grid(:,:,k)     = var(i1,i2,i3)%p4(:,:,k,nnew)
        var_grid_bd1(:,:,k) = var(i1,i2,i3)%p4_bd(:,:,k,nbd1)
        var_grid_bd2(:,:,k) = var(i1,i2,i3)%p4_bd(:,:,k,nbd2)
      ENDDO
      

      ! Loop over all vertical levels that will be nudged
      DO  k = 1, kend_sn

        ! Write actual grid field and boundary field (interpolated in time) 
        ! on 2D-fields
        DO j = 1, je_tot
          DO i = 1, ie_tot
            var_spec(i,j)    = var_grid(i,j,k)
            var_spec_bd(i,j) = zemzdr*var_grid_bd1(i,j,k) +      &
                               zdtibd*var_grid_bd2(i,j,k)
          ENDDO
        ENDDO
    
        ! Calculate Fourier coefficients for actual data for large scales
        CALL spectral_decomposition (var_spec(1:ie_tot,1:je_tot), aa(1:isc_sn, 1:jsc_sn), &
           ab(1:isc_sn, 1:jsc_sn), ba(1:isc_sn, 1:jsc_sn), bb(1:isc_sn, 1:jsc_sn))

        ! Calculate Fourier coefficients for boundary data for large scales
        CALL spectral_decomposition (var_spec_bd(1:ie_tot,1:je_tot), aa_bd(1:isc_sn, 1:jsc_sn), &
        ab_bd(1:isc_sn, 1:jsc_sn), ba_bd(1:isc_sn, 1:jsc_sn), bb_bd(1:isc_sn, 1:jsc_sn))

        ! Nudging of large scales only: j up to 3/2PI, i up to 5/2PI
        DO  i = 1,isc_sn
          DO  j = 1,jsc_sn

            aa_n(i,j) = aa(i,j) - (wgtvd_sn(k))*( aa(i,j) - aa_bd(i,j) )
            ab_n(i,j) = ab(i,j) - (wgtvd_sn(k))*( ab(i,j) - ab_bd(i,j) )
            ba_n(i,j) = ba(i,j) - (wgtvd_sn(k))*( ba(i,j) - ba_bd(i,j) )
            bb_n(i,j) = bb(i,j) - (wgtvd_sn(k))*( bb(i,j) - bb_bd(i,j) )

          ENDDO
        ENDDO

       ! Composition of spectral coefficients back to grid field
        CALL spectral_composition                                          &
              (aa(1:isc_sn, 1:jsc_sn), ab(1:isc_sn, 1:jsc_sn), ba(1:isc_sn, 1:jsc_sn), &
              bb(1:isc_sn, 1:jsc_sn), aa_n(1:isc_sn, 1:jsc_sn), ab_n(1:isc_sn, 1:jsc_sn), &
              ba_n(1:isc_sn, 1:jsc_sn), bb_n(1:isc_sn, 1:jsc_sn), var_spec(1:ie_tot,1:je_tot))

        ! Write nudged 2D-field to 3D-field
        DO j = 1, je_tot
          DO i = 1, ie_tot
            var_grid(i,j,k) = var_spec (i,j)
          ENDDO
        ENDDO 

      ENDDO    ! End of loop over vertical levels
    
        DO  k = 1, kend_sn
          var(i1,i2,i3)%p4(:,:,k,nnew) = var_grid(:,:,k)
        ENDDO

    ENDDO  ! variable_loop

  
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure spectral_nudging
!-------------------------------------------------------------------------------

END SUBROUTINE spectral_nudging

!==============================================================================
!==============================================================================
!+ Module procedure spectral decomposition
!-------------------------------------------------------------------------------

SUBROUTINE spectral_decomposition (array, aafin, abfin, bafin, bbfin)

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure decomposes the grid array into spectral coefficients
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
! Array arguments with intent(in):
  REAL (KIND=wp),     INTENT(IN) :: &
       array(:,:)      ! array for which the Fourier coefficients 
                       ! will be calculated 

! Array arguments with intent(out):
  REAL (KIND=wp),     INTENT(OUT) :: &
       aafin(:,:),   & !\   Fourier coefficients of the actual array
       abfin(:,:),   & ! \  aa = cos cos         ab = cos sin
       bafin(:,:),   & ! /  ba = sin cos         bb = sin sin
       bbfin(:,:)      !/
       
!------------------------------------------------------------------------------
!
! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    n,             & ! loop variable
    ms,            & ! loop variable
    ns,            & ! loop variable
    ir,            & ! loop variable
    jr               ! loop variable

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine spectral_decomposition
!------------------------------------------------------------------------------

  ! Preset some variables
  aafin = 0.0_wp
  bafin = 0.0_wp
  abfin = 0.0_wp
  bbfin = 0.0_wp
  asem = 0.0_wp
  bsem = 0.0_wp

  !First decomposition in i-direction: calculate asem(is,j), bsem(is,j)
  DO n = 1,je_tot
    DO ms = 1,isc_sn
      DO ir = 1,ie_tot
        asem(ms,n)=asem(ms,n)+ array(ir,n)*CosArgDiv(ms,ir,1)
        bsem(ms,n)=bsem(ms,n)+ array(ir,n)*SinArgDiv(ms,ir,1)
      ENDDO
    ENDDO
  ENDDO

  !Second decomposition in j-direction: xsem coefficients as function values
  DO ms = 1,isc_sn
    DO ns = 1,jsc_sn
      DO jr = 1,je_tot
        aafin(ms,ns)=aafin(ms,ns)+ asem(ms,jr)*CosArgDiv(ns,jr,2)
        bafin(ms,ns)=bafin(ms,ns)+ asem(ms,jr)*SinArgDiv(ns,jr,2)
        abfin(ms,ns)=abfin(ms,ns)+ bsem(ms,jr)*CosArgDiv(ns,jr,2)
        bbfin(ms,ns)=bbfin(ms,ns)+ bsem(ms,jr)*SinArgDiv(ns,jr,2)
      ENDDO
    ENDDO
  ENDDO

!-------------------------------------------------------------------------------
! End of module procedure spectral_decomposition
!-------------------------------------------------------------------------------

END SUBROUTINE spectral_decomposition

!==============================================================================
!==============================================================================
!+ Module procedure spectral_composition
!-------------------------------------------------------------------------------

SUBROUTINE spectral_composition (aafin,   abfin,   bafin,   bbfin,   &
                                 aafin_n, abfin_n, bafin_n, bbfin_n, &
                                 array) 

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure converts spectral coefficients back to grid field
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
! Array arguments with intent(in):
  REAL (KIND=wp),     INTENT(IN) :: &
       aafin(:,:),       & !\   Fourier coefficients of the actual array
       abfin(:,:),       & ! \  aa = cos cos         ab = cos sin
       bafin(:,:),       & ! /  ba = sin cos         bb = sin sin
       bbfin(:,:),       & !/
       aafin_n(:,:),     & !\   Fourier coefficients of the nudged array
       abfin_n(:,:),     & ! \  aa = cos cos         ab = cos sin
       bafin_n(:,:),     & ! /  ba = sin cos         bb = sin sin
       bbfin_n(:,:)        !/

! Array arguments with intent(inout):
  REAL (KIND=wp),     INTENT(INOUT) :: &
       array(:,:)          ! grid point array
       
!------------------------------------------------------------------------------
!
! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
       n,                    & ! loop variable
       m,                    & ! loop variable
       ms,                   & ! loop variable
       ns                      ! loop variable

!------------ End of header ---------------------------------------------------



!------------------------------------------------------------------------------
! Begin Subroutine spectral_composition
!------------------------------------------------------------------------------

  ! Presettings
  sum  = 0.0_wp
  sum_n= 0.0_wp
  suma = 0.0_wp
  sumb = 0.0_wp
  asem = 0.0_wp
  bsem = 0.0_wp
  suma_n = 0.0_wp
  sumb_n = 0.0_wp
  asem_n = 0.0_wp
  bsem_n = 0.0_wp

  ! spectral composition, first dimension is j
  DO ms = 1,isc_sn
    DO n  = 1,je_tot
      DO ns = 2,jsc_sn
        suma(ms,n)=suma(ms,n) + aafin(ms,ns)*CosArg(ns,n,2)       &
                    + bafin(ms,ns)*SinArg(ns,n,2)
        sumb(ms,n)=sumb(ms,n) + abfin(ms,ns)*CosArg(ns,n,2)       &
                    + bbfin(ms,ns)*SinArg(ns,n,2)
        suma_n(ms,n)=suma_n(ms,n) + aafin_n(ms,ns)*CosArg(ns,n,2) &
                    + bafin_n(ms,ns)*SinArg(ns,n,2)
        sumb_n(ms,n)=sumb_n(ms,n) + abfin_n(ms,ns)*CosArg(ns,n,2) &
                    + bbfin_n(ms,ns)*SinArg(ns,n,2)
      ENDDO
      asem(ms,n) = aafin(ms,1) + 2.0_wp*suma(ms,n)
      bsem(ms,n) = abfin(ms,1) + 2.0_wp*sumb(ms,n)
      asem_n(ms,n) = aafin_n(ms,1) + 2.0_wp*suma_n(ms,n)
      bsem_n(ms,n) = abfin_n(ms,1) + 2.0_wp*sumb_n(ms,n)
    ENDDO
  ENDDO

  ! Spectral composition, second dimension is i
  DO ms = 2,isc_sn
    DO n = 1,je_tot
      DO m = 1,ie_tot
        sum(m,n)   = sum(m,n)   + asem(ms,n)  *CosArg(ms,m,1)  &
                                + bsem(ms,n)  *SinArg(ms,m,1)
        sum_n(m,n) = sum_n(m,n) + asem_n(ms,n)*CosArg(ms,m,1)  &
                                + bsem_n(ms,n)*SinArg(ms,m,1)
      ENDDO
    ENDDO
  ENDDO
  DO n = 1,je_tot
    DO m = 1,ie_tot
     ! replace the original values in the grid point array by the nudged ones
      array(m,n) = array(m,n)  - asem(1,n)   - 2.0_wp*sum(m,n)    &
                               + asem_n(1,n) + 2.0_wp*sum_n(m,n)
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! End of module procedure spectral_composition
!------------------------------------------------------------------------------

END SUBROUTINE spectral_composition

!==============================================================================

END MODULE src_spectral_nudging
