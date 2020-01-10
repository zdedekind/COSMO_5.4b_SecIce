!+ External Procedure for organizing the semi-implicit scheme
!------------------------------------------------------------------------------

SUBROUTINE semi_implicit (kzflat)

!------------------------------------------------------------------------------
!
! Description:
!   This routine is the organization routine for the fully 3D semi-implicit 
!   scheme. All other routines necessary for this scheme are contained 
!   herein. The semi-implicit scheme for the LM was developed by 
!   Steve Thomas (University of Montreal, now NCAR).
!
! Method:
!   Allocate space for basic operators, compute right-hand
!   sides of the governing equations and wave equation.
!   Iteratively solve the wave equation and back solve to
!   update prognostic variables to time level nnew.
!
! Current Code Owner: Steve Thomas
!  phone:  1 303 497 8194
!  fax:    1 303 497 8171
!  email:  thomas@ncar.ucar.edu
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.34       1999/12/10 Steve Thomas
!  Initial release
! 1.39       2000/05/03 Ulrich Schaettler
!  Adapted call to timing routine get_timing.
!  Changed names for variables concerned to latitude or longitude.
!  Introduced switch lw_freeslip for boundary treatment of w.
! 2.8        2001/07/06 Ulrich Schaettler
!  Corrected bug in time-measuring
! 3.5        2003/09/02 Ulrich Schaettler
!  Adaptations for the interface of exchg_boundaries
! 3.7        2004/02/18 Ulrich Schaettler
!  Renamed cphi (crlat), acphir (acrlat)
! 3.18       2006/03/03 Ulrich Schaettler
!  Changed treatment of ASCII files to introduce restart possibility
! 3.21       2006/12/04 Ulrich Schaettler
!  Adapted output for YUSOLVER
! V4_4         2008/07/16 Ulrich Schaettler
!  Adapted interface of get_timings
! V4_9         2009/07/16 Ulrich Schaettler
!  Check calls to get_timing with IF ltime
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_17        2011/02/24 Ulrich Blahak
!  Adapted interface of exchg_boundaries; corrected kzdims(1:20) -> kzdims(1:24);
!  eliminated my_peri_neigh
! V4_23        2012/05/10 Oliver Fuhrer
!  Eliminated obsolete computed GOTO in SR fmgres
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer, Hans-Juergen Panitz
!  Removed qx, qx-tens variables from declaration list
!  Introduced nexch_tag for MPI boundary exchange tag to replace ntstep (HJP)
! V4_28        2013/07/12 Ulrich Schaettler
!  Pass kflat as parameter to avoid dependency on vgrid_refatm_utils
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
! V5_2         2015-05-21 Ulrich Schaettler
!  Reorganized numbering of MPI datatypes ID
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

USE data_parameters,         ONLY:  wp, iintegers

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------
    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction
    ke1,          & ! ke+1

! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from
!    the other ones because of the use of the staggered Arakawa-B-grid.
!
!   zonal direction
    istart,       & ! start index for the forecast of w, t, qv, qc and pp
    iend,         & ! end index for the forecast of w, t, qv, qc and pp
    istartu,      & ! start index for the forecast of u
    iendu,        & ! end index for the forecast of u
    istartv,      & ! start index for the forecast of v
    iendv,        & ! end index for the forecast of v
    istartpar,    & ! start index for computations in the parallel program
    iendpar,      & ! end index for computations in the parallel program

!   meridional direction
    jstart,       & ! start index for the forecast of w, t, qv, qc and pp
    jend,         & ! end index for the forecast of w, t, qv, qc and pp
    jstartu,      & ! start index for the forecast of u
    jendu,        & ! start index for the forecast of u
    jstartv,      & ! start index for the forecast of v
    jendv,        & ! end index for the forecast of v
    jstartpar,    & ! start index for computations in the parallel program
    jendpar,      & ! end index for computations in the parallel program

! 4. constants for the horizontal rotated grid and related variables
! ------------------------------------------------------------------

    eddlon,       & ! 1 / dlon
    eddlat,       & ! 1 / dlat
    edadlat,      & ! 1 / (radius of the earth * dlat)

! 5. variables for the time discretization and related variables
! --------------------------------------------------------------

    dt,           & ! long time-step
    dt2,          & ! 2 * dt
    betasw          ! beta-variable for treatment of soundwaves

! end of data_modelconfig

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------
    rho0       ,    & ! reference density at the full model levels    (kg/m3)
    dp0        ,    & ! reference pressure thickness of layers        ( Pa)
    p0         ,    & ! reference pressure at main levels             ( Pa)
    hhl        ,    & ! geometical height of half model levels        ( m )

! 2. external parameter fields                                        (unit)
! ----------------------------
    crlat      ,    & ! cosine of transformed latitude
    acrlat     ,    & ! 1 / ( crlat * radius of the earth )           ( 1/m )

! 3. prognostic variables                                             (unit)
! -----------------------
    u          ,    & ! zonal wind speed                              ( m/s )
    v          ,    & ! meridional wind speed                         ( m/s )
    w          ,    & ! vertical wind speed (defined on half levels)  ( m/s )
    t          ,    & ! temperature                                   (  k  )
    pp         ,    & ! deviation from the reference pressure         ( pa  )

! 4. tendency fields for the prognostic variables                     (unit )
! -----------------------------------------------
!    timely deviation  by diabatic and adiabatic processes
!    without sound-wave terms
    utens        ,  & ! u-tendency without sound-wave terms           ( m/s2)
    vtens        ,  & ! v-tendency without sound-wave terms           ( m/s2)
    wtens        ,  & ! w-tendency without sound-wave terms           ( m/s2
                      ! (defined on half levels )
    ttens        ,  & ! t-tendency without sound-wave terms           ( K/s )
    pptens       ,  & ! pp-tendency without sound-wave terms          (Pa/s )

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------
    rho               ! density of moist air

! end of data_fields

!------------------------------------------------------------------------------

USE data_parallel,      ONLY :  &
    num_compute,     & ! number of compute PEs
    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains
    ldatatypes,      & ! if .TRUE.: use MPI-Datatypes for some communications
    ltime_barrier,   & ! if .TRUE.: use additional barriers for determining the
                       ! load-imbalance
    ncomm_type,      & ! type of communication
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    my_cart_neigh,   & ! neighbors of this subdomain in the cartesian grid
    icomm_cart,      & ! communicator for the virtual cartesian topology
    iexch_req,       & ! stores the sends requests for the neighbor-exchange
                       ! that can be used by MPI_WAIT to identify the send
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    nexch_tag,       & ! tag to be used for MPI boundary exchange
                       !  (in calls to exchg_boundaries)
    sendbuf,         & ! sending buffer for boundary exchange:
                       ! 1-4 are used for sending, 5-8 are used for receiving
    isendbuflen        ! length of one column of sendbuf

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
    nold,         & ! corresponds to ntstep - 1
    nnow,         & ! corresponds to ntstep
    nnew,         & ! corresponds to ntstep + 1

! 4. controlling the dynamics
! ---------------------------
    nusolver,     & ! unit number for file YUSOLVER
    yusolver,     & ! file name for file YUSOLVER
    ikrylow_si,   & ! dimension of the Krylow space used in the elliptic
                    ! solver for the semi-implicit scheme
    maxit_si,     & ! maximum number of iterations for the elliptic solver
    iprint_si,    & ! to control whether statistics of the solver are printed
    eps_si,       & ! precision limit for the elliptic solver

! 7. additional control variables
! -------------------------------
    ltime,        & ! detailled timings of the program are given
    lw_freeslip,  & ! if .TRUE.: with free slip lateral boundary condition and
                    ! if .FALSE. specified lateral boundary values for w
    lperi_x,      & ! if lartif_data=.TRUE.: periodic boundary conditions in x-dir.
                    ! or with Davies conditions (.FALSE.)
    lperi_y,      & ! if lartif_data=.TRUE.: periodic boundary conditions in y-dir.
                    ! or with Davies conditions (.FALSE.)
    l2dim           ! 2 dimensional runs

!-------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 2. physical constants and related variables
! -------------------------------------------
    r_d,          & ! gas constant for dry air
    cpdr,         & ! 1 / cp_d
    gamma,        & ! 1 / (1 - rdocp)   ( = cp_d/cv_d)
    g,            & ! acceleration due to gravity
    gq,           & ! g * g
    rdocp           ! r_d / cp_d

!-------------------------------------------------------------------------------

USE parallel_utilities, ONLY :  global_values
USE environment       , ONLY :  exchg_boundaries, comm_barrier
USE time_utilities    , ONLY :  get_timings, i_dyn_computations,            &
                                i_semi_implicit, i_semi_impl_comm,          &
                                i_semi_impl_barrier

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Subroutine arguments:
! ---------------------
  INTEGER (KIND=iintegers), INTENT (IN) ::  &
    kzflat                ! level-index where half-levels bcome flat

! Local (automatic) arrays:
! -------------------------
  REAL (KIND = wp)                  :: &
    zqu (ie,je,ke),    & ! u right-hand side
    zqv (ie,je,ke),    & ! v right-hand side
    zqw (ie,je,ke),    & ! w right-hand side
    zqp (ie,je,ke),    & ! p right-hand side
    zqt (ie,je,ke),    & ! T right-hand side
    zqws(ie,je,ke),    & ! auxilary w right-hand side
    zqps(ie,je,ke),    & ! wave equation right-hand side
    za  (ie,je,ke),    & ! lower diagonal of N matrix
    zb  (ie,je,ke),    & ! main  diagonal of N matrix
    zc  (ie,je,ke),    & ! upper diagonal of N matrix
    za1 (ie,je,ke),    & ! lower diagonal of Dz1
    zb1 (ie,je,ke),    & ! main  diagonal of Dz1
    za2 (ie,je,ke),    & ! lower diagonal of Dz2
    zb2 (ie,je,ke),    & ! main  diagonal of Dz2
    za3 (ie,je,ke),    & ! lower diagonal of Dz2.Ni.Dz1 matrix
    zb3 (ie,je,ke),    & ! main  diagonal of Dz2.Ni.Dz1 matrix
    zc3 (ie,je,ke)       ! upper diagonal of Dz2.Ni.Dz1 matrix

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    im,                  & ! Number of Krylov vectors
    im1,                 & ! im+1
    kzdims(24),          & ! Loop index in vertical direction 
    izerror,             & ! error status variable
    niostat                ! error status for IO of Ascii file

! End of header
!==============================================================================
 
!------------------------------------------------------------------------------
! Begin Subroutine semi_implicit   
!------------------------------------------------------------------------------

  kzdims(:) = 0_iintegers
  izerror   = 0_iintegers
  IF (ltime) CALL get_timings (i_dyn_computations, ntstep, dt, izerror)

!------------------------------------------------------------------------------
! Section 1: Allocate memory for the right-hand sides and set to zero
!------------------------------------------------------------------------------

  zqu (:,:,:) = 0.0_wp 
  zqv (:,:,:) = 0.0_wp 
  zqw (:,:,:) = 0.0_wp 
  zqp (:,:,:) = 0.0_wp 
  zqt (:,:,:) = 0.0_wp 
  zqws(:,:,:) = 0.0_wp 
  zqps(:,:,:) = 0.0_wp 
  za  (:,:,:) = 0.0_wp 
  zb  (:,:,:) = 0.0_wp 
  zc  (:,:,:) = 0.0_wp 
  za1 (:,:,:) = 0.0_wp 
  zb1 (:,:,:) = 0.0_wp 
  za2 (:,:,:) = 0.0_wp 
  zb2 (:,:,:) = 0.0_wp 
  za3 (:,:,:) = 0.0_wp 
  zb3 (:,:,:) = 0.0_wp 
  zc3 (:,:,:) = 0.0_wp 

!------------------------------------------------------------------------------
! Section 2: Initialize vertical operators and preconditioner
!------------------------------------------------------------------------------

  IF ( (iprint_si > 0) .AND. (my_cart_id == 0) ) THEN
    OPEN(nusolver, FILE=yusolver, FORM='FORMATTED', STATUS='OLD',     &
                   POSITION='APPEND', IOSTAT=niostat)
  ENDIF

  CALL init_implicit
  CALL init_precond

!------------------------------------------------------------------------------
! Section 3: Compute the right-hand sides
!------------------------------------------------------------------------------

  CALL rhs

!------------------------------------------------------------------------------
! Section 4: Compute the right-hand side of the wave equation
!------------------------------------------------------------------------------

  CALL wave_rhs

!------------------------------------------------------------------------------
! Section 5: Solve the wave equation for the pressure tendency
!------------------------------------------------------------------------------

  CALL wave_solve

!------------------------------------------------------------------------------
! Section 6: Back solve for tendencies and update dynamics fields
!------------------------------------------------------------------------------

  CALL back_solve

  IF ( (iprint_si > 0) .AND. (my_cart_id == 0) ) THEN
    CLOSE (nusolver, STATUS='KEEP')
  ENDIF

  IF (ltime) CALL get_timings (i_semi_implicit, ntstep, dt, izerror)

!==============================================================================

! Internal procedures for the semi-implicit scheme

CONTAINS

!==============================================================================
!+ Internal procedure in "semi_implicit" for initialization
!------------------------------------------------------------------------------

SUBROUTINE init_implicit

!------------------------------------------------------------------------------
!
! Description:
!   This internal procedure of "semi_implicit" computes matrix elements 
!   of the vertical operators N, Dz1, and Dz2.
!
! Method:
!   The 3D semi-implicit scheme employs a linearisation about the
!   the mid-point (time level n) values of p, rho, and T. The vertical
!   coordinate is also based on the base state pressure p0, rho0, dp0.
!   Therefore the matrices N, Dz1 and Dz2 have coefficients which vary
!   with the time step and also with horizontal coordinates.
!
!------------------------------------------------------------------------------

! Local scalars:
! -----------------------------
  INTEGER (KIND=iintegers) ::  &
    i,  j,  k              !  Loop indices in longitudinal, latitudinal and
                           !  vertical direction 

  REAL    (KIND=wp   )     ::  &
    zbt,                 & !
    zpn,                 & !
    zpfac, zcp,          & !
    zrofac, zcw1,        & !
    zbk1, zbk2,          & !
    zn02                   !
   
! Local (automatic) arrays:
! -----------------------------
! REAL    (KIND=wp   )     ::  &

!
! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine init_implicit
!------------------------------------------------------------------------------
 
zbt = (1.0_wp + betasw) * dt

! Initialize dimensions of GMRES work vectors
im     = ikrylow_si
im1    = im + 1

! Compute elements of the vertical derivative matrix Dz2

DO k = 1, ke
  DO j = jstart, jend
    DO  i = istart, iend
      zpn        = p0(i,j,k) + pp(i,j,k,nnow)
      zpfac      = 0.5_wp * dp0(i,j,k) / ( gamma * zpn )
      zcp        = - g * rho0(i,j,k) / dp0(i,j,k)
      zb2(i,j,k) = zcp * ( 1.0_wp + zpfac )
      za2(i,j,k) = zcp * ( 1.0_wp - zpfac )
    ENDDO
  ENDDO
ENDDO

! Compute elements of the vertical derivative matrix Dz1

DO k = 2, ke
  DO j = jstart, jend
    DO  i = istart, iend
      zrofac = ( rho0(i,j,k)*dp0(i,j,k-1) + rho0(i,j,k-1)*dp0(i,j,k) ) /    &
               ( rho (i,j,k)*dp0(i,j,k-1) + rho (i,j,k-1)*dp0(i,j,k) )
      zcw1   = - 2.0_wp * g * zrofac / ( dp0(i,j,k-1) + dp0(i,j,k) )
      zpfac  = 0.5_wp * dp0(i,j,k-1) / ( gamma * p0(i,j,k) )
      zb1(i,j,k) = zcw1 * ( 1 - zpfac )
      zpfac  = 0.5_wp * dp0(i,j,k)   / ( gamma * p0(i,j,k-1) )
      za1(i,j,k) = zcw1 * ( 1 + zpfac )
    ENDDO
  ENDDO
ENDDO

! Compute elements of the Brunt-Vaisalla matrix N

DO k = 2, ke
  DO j = jstart, jend
    DO  i = istart, iend
      zrofac = ( rho0(i,j,k)*dp0(i,j,k-1) + rho0(i,j,k-1)*dp0(i,j,k) ) /    &
               ( rho (i,j,k)*dp0(i,j,k-1) + rho (i,j,k-1)*dp0(i,j,k) )
      zcw1   = zbt * zbt * 0.5_wp * zrofac / ( dp0(i,j,k-1) + dp0(i,j,k) )
      zn02   = gq * rdocp * rho0(i,j,k) / p0(i,j,k)
      zc(i,j,k) = zcw1 * dp0(i,j,k-1) * zn02
      zn02   = gq * rdocp * rho0(i,j,k-1) / p0(i,j,k-1)
      za(i,j,k) = zcw1 * dp0(i,j,k)   * zn02
      zb(i,j,k) = 1.0_wp + za(i,j,k) + zc(i,j,k)
!     zrofac = ( rho0(i,j,k) / rho (i,j,k) )
!     zn02   = gq * rdocp * rho0(i,j,k-1) / p0(i,j,k-1)
!     zb(i,j,k) = 1.0_wp/(1.0_wp + zn02 * zbt * zbt * zrofac )
!     zb(i,j,k) = 1.0_wp/(1.0_wp + zn02 * zbt * zbt )
    ENDDO
  ENDDO
ENDDO

! Compute tridiagonal elements of Dz2.Ni.Dz1 (mass-lumped N)

DO k = 2, ke-1
  DO j = jstart, jend
    DO  i = istart, iend
      zbk1       = 1.0_wp/( 1.0_wp + 2.0_wp*( za(i,j,k)   + zc(i,j,k) ) )
      zbk2       = 1.0_wp/( 1.0_wp + 2.0_wp*( za(i,j,k+1) + zc(i,j,k+1) ) )
      za3(i,j,k) = za2(i,j,k) * zbk1 * za1(i,j,k)
      zc3(i,j,k) = zb2(i,j,k) * zbk2 * zb1(i,j,k+1)
      zb3(i,j,k) = - za2(i,j,k) * zbk1 * zb1(i,j,k) -                      &
                     zb2(i,j,k) * zbk2 * za1(i,j,k+1)
    ENDDO
  ENDDO
ENDDO

DO j = jstart, jend
  DO  i = istart, iend
    zbk2       = 1.0_wp/( 1.0_wp + 2.0_wp*( za(i,j,2) + zc(i,j,2) ) )
    zc3(i,j,1) = zb2(i,j,1) * zbk2 * zb1(i,j,2)
    zb3(i,j,1) = - zb2(i,j,1) * zbk2 * za1(i,j,2)
  ENDDO
ENDDO

DO j = jstart, jend
  DO  i = istart, iend
    zbk1        = 1.0_wp/( 1.0_wp + 2.0_wp*( za(i,j,ke) + zc(i,j,ke) ) )
    za3(i,j,ke) = za2(i,j,ke) * zbk1 * za1(i,j,ke)
    zb3(i,j,ke) = - za2(i,j,ke) * zbk1 * zb1(i,j,ke)
  ENDDO
ENDDO

!------------------------------------------------------------------------------
!  End of internal procedure init_implicit
!------------------------------------------------------------------------------

END SUBROUTINE init_implicit

!==============================================================================
!+ Internal procedure in "semi_implicit" to initialize preconditioner
!------------------------------------------------------------------------------

SUBROUTINE init_precond

!------------------------------------------------------------------------------
!
! Description:
!   This internal procedure of "semi_implicit" initializes
!   the vertical line Jacobi preconditioner.
!
! Method:
!   Compute the diagonal coefficients of the horizontal 
!   elliptic operator and then add to the vertical.
!   Compute the LU decomposition of the resulting
!   tridiagonal preconditioning matrix P.
!
!------------------------------------------------------------------------------

! Local scalars:
! -----------------------------
  INTEGER (KIND=iintegers) ::  &
    i,  j,  k,           & !  Loop indices in longitudinal, latitudinal and
                           !  vertical direction 
    ilowu, ilowv,        & !  start- and end-indices for subdomains with
    jlowu, jlowv,        & !  or without neighbours
    iup,   jup             !  

  REAL    (KIND=wp   )     ::  &
    zbt,                 & !
    zfy,                 & !
    zrhoqx, zrhoqy,      & !
    zpn,                 & !
    zpfac,               & !
    zdenom

   
! Local (automatic) arrays:
! -----------------------------
  REAL    (KIND=wp   )     ::  &
    zfx     (je   )     ,      & !
    zfydn   (je   )     ,      & !
    zfyds   (je   )     ,      & !
    zrhoqx_i(ie,je,ke)  ,      & !
    zrhoqy_i(ie,je,ke)  ,      & !
    zdi     (ie,je,ke)           !
!
! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine init_precond
!------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------
!  Section 1: Setup of subdomain indices 
!------------------------------------------------------------------------------

zbt = (1.0_wp + betasw) * dt

! Compute local subdomain start and end indices

IF (my_cart_neigh(1) == -1) THEN ! west
  ilowu = istartu
  ilowv = istartv
ELSE
  ilowu = istartu-1
  ilowv = istartv-1
ENDIF

IF (my_cart_neigh(2) == -1) THEN ! north
  jup = jend
ELSE
  jup = jend+1
ENDIF

IF (my_cart_neigh(3) == -1) THEN ! east
  iup = iend
ELSE
  iup = iend+1
ENDIF

IF (my_cart_neigh(4) == -1) THEN ! south
  jlowu = jstartu
  jlowv = jstartv
ELSE
  jlowu = jstartu-1
  jlowv = jstartv-1
ENDIF

!------------------------------------------------------------------------------
!  Section 2: Some preparations to avoid division
!------------------------------------------------------------------------------

zfy = edadlat*edadlat
DO j = jstart-1, jend+1
   zfx(j)   = eddlon*acrlat(j,1)
   zfx(j)   = zfx(j)*zfx(j)
   zfydn(j) = crlat(j  ,2)/crlat(j,1)
   zfyds(j) = crlat(j-1,2)/crlat(j,1)
ENDDO          

! Avoid division by zrhoqx and zrhoqy 

DO  k = 1, ke
  DO  j = jlowu, jendu
    DO  i = ilowu, iendu
      zrhoqx          = 0.5_wp*( rho(i,j,k) + rho(i+1,j,k) )
      zrhoqx_i(i,j,k) = zfx(j) / zrhoqx
    ENDDO
    IF ( j <= jendv ) THEN
      DO  i = ilowv, iendv
        zrhoqy          = 0.5_wp*( rho(i,j,k) + rho(i,j+1,k) )
        zrhoqy_i(i,j,k) = zfy / zrhoqy
      ENDDO
    ENDIF
  ENDDO
ENDDO

! Compute diagonal of horizontal -div.1/rho.grad

! South-West corner i=istart, j=jstart

IF ( (my_cart_neigh(1) == -1) .AND.                                     &
     (my_cart_neigh(4) == -1) ) THEN 
DO  k = 1, ke
  zdi(istart,jstart,k) =                                                &
    zrhoqx_i(istart,jstart,k) + zrhoqy_i(istart,jstart,k)*zfydn(jstart)
ENDDO
ENDIF

! South boundary j=jstart

IF (my_cart_neigh(4) == -1) THEN ! south
DO  k = 1, ke
  DO  i = ilowu+1, iup-1
    zdi(i,jstart,k) =                                                   &
      zrhoqx_i(i,jstart,k) + zrhoqx_i(i-1,jstart,k) +                   &
      zrhoqy_i(i,jstart,k)*zfydn(jstart)
  ENDDO
ENDDO
ENDIF

! South-East corner i=iend, j=jstart

IF ( (my_cart_neigh(3) == -1) .AND.                                     &
     (my_cart_neigh(4) == -1) ) THEN 
DO  k = 1, ke
  zdi(iend,jstart,k) =                                                  &
    zrhoqx_i(iend-1,jstart,k) + zrhoqy_i(iend,jstart,k)*zfydn(jstart)
ENDDO
ENDIF

! West boundary i=istart

IF (my_cart_neigh(1) == -1) THEN ! west
DO  k = 1, ke
  DO  j = jlowv+1, jup-1
    zdi(istart,j,k) =                                                   &
      zrhoqx_i(istart,j,k) +                                            &
      zrhoqy_i(istart,j,k)*zfydn(j) + zrhoqy_i(istart,j-1,k)*zfyds(j)
  ENDDO
ENDDO
ENDIF

! Internal Points

DO  k = 1, ke
  DO  j = jlowv+1, jup-1
    DO  i = ilowu+1, iup-1
      zdi(i,j,k) =                                                      &
        zrhoqx_i(i,j,k) + zrhoqx_i(i-1,j,k) +                           &
        zrhoqy_i(i,j,k)*zfydn(j) + zrhoqy_i(i,j-1,k)*zfyds(j)
    ENDDO
  ENDDO
ENDDO

! East boundary i=iend

IF (my_cart_neigh(3) == -1) THEN ! east 
DO  k = 1, ke
  DO  j = jlowv+1, jup-1
    zdi(iend,j,k) =                                                     &
      zrhoqx_i(iend-1,j,k) +                                            &
      zrhoqy_i(iend,j,k)*zfydn(j) + zrhoqy_i(iend,j-1,k)*zfyds(j)
  ENDDO
ENDDO
ENDIF

! North-West corner i=istart, j=jend

IF ( (my_cart_neigh(1) == -1) .AND.                                     &
     (my_cart_neigh(2) == -1) ) THEN 
DO  k = 1, ke
  zdi(istart,jend,k) =                                                  &
    zrhoqx_i(istart,jend,k) + zrhoqy_i(istart,jend-1,k)*zfyds(jend)
ENDDO
ENDIF

! North boundary j=jend

IF (my_cart_neigh(2) == -1) THEN ! north
DO  k = 1, ke
  DO  i = ilowu+1, iup-1
    zdi(i,jend,k) =                                                     &
      zrhoqx_i(i,jend,k) + zrhoqx_i(i-1,jend,k) +                       &
      zrhoqy_i(i,jend-1,k)*zfyds(jend)
  ENDDO
ENDDO
ENDIF

! North-East corner i=iend, j=jend

IF ( (my_cart_neigh(3) == -1) .AND.                                     &
     (my_cart_neigh(2) == -1) ) THEN 
DO  k = 1, ke
  zdi(iend,jend,k) =                                                    &
    zrhoqx_i(iend-1,jend,k) + zrhoqy_i(iend-1,jend-1,k)*zfyds(jend)
ENDDO
ENDIF

! Compute diagonal of  -div.1/rho.grad - Dz2.Ni.Dz1.d(pp)

DO j = jstart, jend
  DO  i = istart, iend
    zdi(i,j,1) = zdi(i,j,1) - zb3(i,j,1)
  ENDDO
ENDDO

DO k = 2, ke-1
  DO j = jstart, jend
    DO  i = istart, iend
      zdi(i,j,k) = zdi(i,j,k) - zb3(i,j,k)
    ENDDO
  ENDDO
ENDDO

DO j = jstart, jend
  DO  i = istart, iend
    zdi(i,j,ke) = zdi(i,j,ke) - zb3(i,j,ke)
  ENDDO
ENDDO

! P = - [1/gamma.pn + dt^2.( -div.(1/rho).grad - Dz2.Ni.Dz1 ) ]
! P = -  1/gamma.pn - dt^2.( -div.(1/rho).grad - Dz2.Ni.Dz1 ) 

DO k = 1, ke
  DO j = jstart, jend
    DO  i = istart, iend
      zpn         = p0(i,j,k) + pp(i,j,k,nnow)
      zpfac       = 1.0_wp / ( gamma * zpn )
      zb3(i,j,k)  = - zpfac - zbt * zbt * zdi(i,j,k)
      zc3(i,j,k)  = zbt * zbt * zc3(i,j,k)
      za3(i,j,k)  = zbt * zbt * za3(i,j,k)
    ENDDO
  ENDDO
ENDDO

! LU Decomposition P = L.U

DO j = jstart, jend
  DO  i = istart, iend
    zdenom     = zb3(i,j,1)
    zb3(i,j,1) = 1.0_wp / zdenom
    zc3(i,j,1) = - zc3(i,j,1) * zb3(i,j,1)
  ENDDO
ENDDO

DO k = 2, ke
  DO j = jstart, jend
    DO  i = istart, iend
      zdenom     = za3(i,j,k) * zc3(i,j,k-1) + zb3(i,j,k)
      zb3(i,j,k) = 1.0_wp / zdenom
      zc3(i,j,k) = - zc3(i,j,k) * zb3(i,j,k)
    ENDDO
  ENDDO
ENDDO


!------------------------------------------------------------------------------
!  End of internal procedure init_precond
!------------------------------------------------------------------------------

END SUBROUTINE init_precond

!==============================================================================
!+ Internal procedure in "semi_implicit" for computing rhs of prog. equations
!------------------------------------------------------------------------------

SUBROUTINE rhs

!------------------------------------------------------------------------------
!
! Description:
!   This internal procedure of "semi_implicit" calculates the 
!   right-hand sides of the prognostic equations at time level n-1
!   for the semi-implicit scheme (zqu, zqv, zqw, zqp, zqt).
!
! Method:
!   A three time level semi-implicit scheme is implemented with
!   a leapfrog interval 2*dt. Slow tendencies of the prognostic variables
!   due to advection and diffusion (which have been calculated before and
!   are stored in the corresponding tendency fields) are evaluated at 
!   the midpoint time level (nnow).
!
!------------------------------------------------------------------------------

! Local scalars:
! -----------------------------
  INTEGER (KIND=iintegers) ::  &
    i,  j,  k,           & !  Loop indices in longitudinal, latitudinal and
                           !  vertical direction 
    kp1, km1,            & !  k+1, k-1
    ilowu, ilowv,        & !  start- and end-indices for subdomains with
    jlowu, jlowv,        & !  or without neighbours
    iup,   jup             !  

  REAL    (KIND=wp   )     ::  &
    zfy,                 & !
    zrhoqx, zrhoqy,      & !
    zpdx, zpdy,          & !
    zdppz,               & !
    zdpq1, zdpq2, zdpq3, & !
    zppqx, zppqy,        & !
    zpgradx, zpgrady,    & !
    ztp, zrofac,         & !
    zbuoy, zfact,        & !
    zcw1, zpn,           & !
    zcp1, zcp2, zct1,    & !
    zpfac, ztfac           !

   
! Local (automatic) arrays:
! -----------------------------
  REAL    (KIND=wp   )     ::  &
    zpb     (ie,je,ke)  ,      & !
    zppqz   (ie,je)     ,      & !
    ztb     (ie,je,ke)  ,      & !
    zfx     (je   )     ,      & !
    zfydn   (je   )     ,      & !
    zfyds   (je   )     ,      & !
    zpi     (ie,je,ke)  ,      & !
    zp0u    (ie,je,ke)  ,      & !
    zp0v    (ie,je,ke)  ,      & !
    zrhoqx_i(ie,je,ke)  ,      & !
    zrhoqy_i(ie,je,ke)           !

! For error handling
! ------------------
  INTEGER (KIND=iintegers) ::  &
    izerror

  CHARACTER (LEN=80)       ::  &
    yzerrmsg

!
! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine rhs
!------------------------------------------------------------------------------
 
izerror  = 0
yzerrmsg = '   '

!------------------------------------------------------------------------------
!  Section 1: Setup of subdomain indices to compute pressure gradients
!------------------------------------------------------------------------------

! zpi has to be set to 0 (at least at the boundary)
zpi(:,:,:) = 0.0_wp 

! Communication of 1 boundary line for utens and vtens
IF (ltime) THEN
  CALL get_timings (i_semi_implicit, ntstep, dt, izerror)
  IF (ltime_barrier) THEN
    CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
    CALL get_timings (i_semi_impl_barrier, ntstep, dt, izerror)
  ENDIF
ENDIF

kzdims(1:24)=(/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
CALL exchg_boundaries                                                   &
   (15  ,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,    &
    ie, je, kzdims, jstartpar, jendpar, 1, nboundlines,                 &
    my_cart_neigh, lperi_x, lperi_y, l2dim,                             &
    10000+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,         &
    utens(:,:,:), vtens(:,:,:))
IF (ltime) CALL get_timings (i_semi_impl_comm, ntstep, dt, izerror)

! Compute local subdomain start and end indices

IF (my_cart_neigh(1) == -1) THEN ! west
  ilowu = istartu
  ilowv = istartv
ELSE
  ilowu = istartu-1
  ilowv = istartv-1
ENDIF

IF (my_cart_neigh(2) == -1) THEN ! north
  jup = jend
ELSE
  jup = jend+1
ENDIF

IF (my_cart_neigh(3) == -1) THEN ! east
  iup = iend
ELSE
  iup = iend+1
ENDIF

IF (my_cart_neigh(4) == -1) THEN ! south
  jlowu = jstartu
  jlowv = jstartv
ELSE
  jlowu = jstartu-1
  jlowv = jstartv-1
ENDIF

!------------------------------------------------------------------------------
!  Section 2: Some preparations to avoid division
!------------------------------------------------------------------------------

zfy = edadlat
DO j = jstart-1, jend+1
   zfx(j)   = eddlon*acrlat(j,1)
   zfydn(j) = eddlat*acrlat(j,1)*crlat(j  ,2)
   zfyds(j) = eddlat*acrlat(j,1)*crlat(j-1,2)
ENDDO          

! Avoid division by zrhoqx and zrhoqy 

DO  k = 1, ke
  DO  j = jlowu, jendu
    DO  i = ilowu, iendu
      zrhoqx          = 0.5_wp*( rho(i,j,k) + rho(i+1,j,k) )
      zrhoqx_i(i,j,k) = zfx(j) / zrhoqx
    ENDDO
    IF ( j <= jendv ) THEN
      DO  i = ilowv, iendv
        zrhoqy          = 0.5_wp*( rho(i,j,k) + rho(i,j+1,k) )
        zrhoqy_i(i,j,k) = zfy / zrhoqy
      ENDDO
    ENDIF
  ENDDO
ENDDO

DO  k = 1, ke
  kp1 = MIN(k+1, ke)
  km1 = MAX(k-1, 1 )
  DO  j = jstart-1, jend
    DO  i = istart-1, iend
      zdpq1   =   p0(i  ,j  ,kp1) - p0(i  ,j  ,km1)
      zdpq2   =   p0(i+1,j  ,kp1) - p0(i+1,j,  km1)
      zdpq3   =   p0(i  ,j+1,kp1) - p0(i  ,j+1,km1)
      zp0u(i,j,k) =  (p0(i+1,j,k) - p0(i  ,j,k)) / (zdpq2 + zdpq1)
      zp0v(i,j,k) =  (p0(i,j+1,k) - p0(i  ,j,k)) / (zdpq3 + zdpq1)
    ENDDO
  ENDDO
ENDDO

!------------------------------------------------------------------------------
!  Section 3: Compute the right-hand sides zqu, zqv, zqw, zqp, zqt
!------------------------------------------------------------------------------

! Save pressure perturbation pp(nold)

DO  k = 1, ke
  zpi(istart-1:iend+1,jstart-1:jend+1,k) =                              &
                pp(istart-1:iend+1,jstart-1:jend+1,k,nold)
ENDDO         

! Lateral boundary conditions on zqu and zqv

IF (my_cart_neigh(1) == -1) THEN  ! west
  DO k = 1, ke
    DO j = jstartu, jendu
      zqu(istartu-1,j,k) =                                              &
        ( u(istartu-1,j,k,nnew) - u(istartu-1,j,k,nold) ) / dt2
    ENDDO
  ENDDO
END IF

IF (my_cart_neigh(3) == -1) THEN  ! east
  DO k = 1, ke
    DO j = jstartu, jendu
      zqu(iendu+1,j,k) =                                                &
        ( u(iendu+1,j,k,nnew) - u(iendu+1,j,k,nold) ) / dt2
    ENDDO
  ENDDO
END IF

IF (my_cart_neigh(4) == -1) THEN  ! south
  DO k = 1, ke
    DO i = istartv, iendv
      zqv(i,jstartv-1,k) =                                              &
        ( v(i,jstartv-1,k,nnew) - v(i,jstartv-1,k,nold) ) / dt2
    ENDDO
  ENDDO
END IF

IF (my_cart_neigh(2) == -1) THEN  ! north
  DO k = 1, ke
    DO i = istartv, iendv
      zqv(i,jendv+1,k) =                                                &
        ( v(i,jendv+1,k,nnew) - v(i,jendv+1,k,nold) ) / dt2
    ENDDO
  ENDDO
END IF

! Compute zqu and zqv (horizontal momentum)

! No communication since zpi available in halo

DO  k = 1, ke

  kp1   = MIN ( ke, k+1 )
  km1   = MAX ( 1 , k-1 )

  IF ( k < kzflat ) THEN
    DO  j = jlowu, jendu
      DO  i = ilowu, iendu   
        zpgradx =  (zpi(i+1,j,k) - zpi(i,j,k)) * zrhoqx_i(i,j,k)
        zqu(i,j,k) = utens(i,j,k) - zpgradx
      ENDDO 
      IF ( j <= jendv ) THEN
        DO  i = ilowv, iendv   
          zpgrady =  (zpi(i,j+1,k) - zpi(i,j,k)) * zrhoqy_i(i,j,k)
          zqv(i,j,k) = vtens(i,j,k) - zpgrady
        ENDDO 
      ENDIF
    ENDDO   
  ELSE 
    zppqz(istart-1:iend+1,jstart-1:jend+1) =                  &
                  zpi(istart-1:iend+1,jstart-1:jend+1,kp1)    &
                 -zpi(istart-1:iend+1,jstart-1:jend+1,km1)

    DO  j = jlowu, jendu
      DO  i = ilowu, iendu   
        zppqx   =  zppqz(i+1,j) + zppqz(i,j)
        zdppz   =  zppqx * zp0u(i,j,k)
        zpdx    =  zpi(i+1,j,k) - zpi(i,j,k)
        zpgradx =  zpdx - zdppz
        zqu(i,j,k) = utens(i,j,k) - zpgradx * zrhoqx_i(i,j,k)
      ENDDO 
      IF ( j <= jendv ) THEN
        DO  i = ilowv, iendv   
          zppqy   =  zppqz(i,j+1) + zppqz(i,j)
          zdppz   =  zppqy * zp0v(i,j,k)
          zpdy    =  zpi(i,j+1,k) - zpi(i,j,k)
          zpgrady =  zpdy - zdppz
          zqv(i,j,k) = vtens(i,j,k) - zpgrady * zrhoqy_i(i,j,k)
        ENDDO 
      ENDIF
    ENDDO   
  ENDIF

ENDDO

! Compute zqw (vertical momentum)

DO k = 1, ke
  DO j = jstart, jend
    DO  i = istart, iend
      zpn         = p0(i,j,k) + pp(i,j,k,nnow)
      zpb (i,j,k) = 0.5_wp * pp(i,j,k,nold) / p0(i,j,k)
      ztp         = t(i,j,k,nold) - p0(i,j,k) / ( r_d * rho0(i,j,k) )
      ztb (i,j,k) = 0.5_wp * ztp / t(i,j,k,nnow) * zpn / p0(i,j,k)
    ENDDO
  ENDDO
ENDDO

DO k = 2, ke
  DO j = jstart, jend
    DO  i = istart, iend
      zrofac = ( rho0(i,j,k)*dp0(i,j,k-1) + rho0(i,j,k-1)*dp0(i,j,k) ) /    &
               ( rho (i,j,k)*dp0(i,j,k-1) + rho (i,j,k-1)*dp0(i,j,k) )
      zcw1   = 2.0_wp * g * zrofac / ( dp0(i,j,k-1) + dp0(i,j,k) )
      zbuoy  = zcw1 * (                                                     &
              ( ztb(i,j,k)   - zpb(i,j,k)   ) * dp0(i,j,k-1) +              &
              ( ztb(i,j,k-1) - zpb(i,j,k-1) ) * dp0(i,j,k)   )
      zqw(i,j,k) = wtens(i,j,k) + zcw1*( zpi(i,j,k) - zpi(i,j,k-1) ) + zbuoy
    ENDDO
  ENDDO
ENDDO

! Compute horizontal divergence Dh( u(nold), v(nold) )

! No communication since u(nold) and v(nold) available in halo

DO k = 1, ke
  kp1 = MIN( ke, k+1 )
  km1 = MAX( 1 , k-1 )
  zfact = 1.0_wp
  IF(k==1 .OR. k==ke) zfact=0.5_wp
  IF ( k < kzflat ) THEN
    DO j = jstart, jend
    DO i = istart, iend
      zpi(i,j,k) = zfx(j)*( u(i,j,k,nold) - u(i-1,j,k,nold) ) &
                  + zfydn(j)*v(i,j,k,nold) - zfyds(j)*v(i,j-1,k,nold)
    ENDDO
    ENDDO
  ELSE
    DO j = jstart, jend
    DO i = istart, iend
      zpi(i,j,k) =                                                         &
       zfx(j)*( u(i,j,k,nold) - u(i-1,j,k,nold)                            &
        - zfact * (u(i-1,j,kp1,nold) - u(i-1,j,km1,nold))* zp0u(i-1,j,k)   &
        - zfact * (u(i  ,j,kp1,nold) - u(i  ,j,km1,nold))* zp0u(i  ,j,k) ) &
        + zfydn(j)*(  v(i,j  ,k,nold)                                      &
        - zfact * (v(i,j,kp1,nold)   - v(i,j,km1,nold))  * zp0v(i,j,k) )   &
        + zfyds(j)*(- v(i,j-1,k,nold)                                      &
        - zfact * (v(i,j-1,kp1,nold) - v(i,j-1,km1,nold)) * zp0v(i,j-1,k) )
    ENDDO
    ENDDO
  ENDIF
ENDDO      

! Compute zqp (pressure)

DO k = 1, ke
  DO j = jstart, jend
    DO  i = istart, iend
      zpn        = p0(i,j,k) + pp(i,j,k,nnow)
      zpfac      = gamma * zpn 
      zcp1       = zpfac * g * rho0(i,j,k) / dp0(i,j,k)
      zcp2       = 0.5_wp * g * rho0(i,j,k)
      zqp(i,j,k) = pptens(i,j,k) - zpfac * zpi(i,j,k) +                   &
                   zcp1 * ( w(i,j,k+1,nold) - w(i,j,k,nold) ) +           &
                   zcp2 * ( w(i,j,k+1,nold) + w(i,j,k,nold) )
    ENDDO
  ENDDO
ENDDO

! Compute zqt (temperature)

DO k = 1, ke
  DO j = jstart, jend
    DO  i = istart, iend
      ztfac      = gamma * cpdr * r_d * t(i,j,k,nnow)
      zct1       = ztfac * g * rho0(i,j,k) / dp0(i,j,k)
      zqt(i,j,k) = ttens(i,j,k) - ztfac * zpi(i,j,k) +                    &
                   zct1 * ( w(i,j,k+1,nold) - w(i,j,k,nold) ) 
    ENDDO
  ENDDO
ENDDO

!------------------------------------------------------------------------------
!  End of internal procedure rhs
!------------------------------------------------------------------------------

END SUBROUTINE rhs

!==============================================================================
!+ Internal procedure in "semi_implicit" for computing rhs of wave equation
!------------------------------------------------------------------------------

SUBROUTINE wave_rhs

!------------------------------------------------------------------------------
!
! Description:
!   This internal procedure of "semi_implicit" calculates the 
!   right-hand side of the implicit wave equation.
!
! Method:
!   Elimination of the buoyancy and divergence from the prognostic
!   equations leads to a fully 3D semi-implicit wave equation for
!   the pressure perturbation pp tendency. The elimination steps
!   applied to the right-hand sides (zqu, zqv, zqw, zqp, zqt)
!   result in three auxilary equations with rhs (zqws, zqpp, zqps).
!   zqps is the rhs of the wave equation. zqpp and zqps share the
!   the same storage. Top (k=1) and bottom (k=ke+1) boundary conditions 
!   are imposed on zqws.
!
!------------------------------------------------------------------------------

! Local scalars:
! -----------------------------
  INTEGER (KIND=iintegers) ::  &
    i,  j,  k,           & !  Loop indices in longitudinal, latitudinal and
                           !  vertical direction 
    kp1, km1               !  k+1, k-1

  REAL    (KIND=wp   )     ::  &
    zbt,                 & !
    zdpq1, zdpq2, zdpq3, & !
           ztp, zrofac,  & !
    zbuoy, zfact,        & !
    zcw1,  zpn,          & !
    zpfac,               & !
    zjpvn, zjpvs,        & !
    zjlur, zjlul, zdenom   !

   
! Local (automatic) arrays:
! -----------------------------
  REAL    (KIND=wp   )     ::  &
    zpb     (ie,je,ke)  ,      & !
    ztb     (ie,je,ke)  ,      & !
    zfx     (je   )     ,      & !
    zfydn   (je   )     ,      & !
    zfyds   (je   )     ,      & !
    zpi     (ie,je,ke)  ,      & !
    zp0u    (ie,je,ke)  ,      & !
    zp0v    (ie,je,ke)  ,      & !
    zd      (ie,je,ke1)          !
!
! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine wave_rhs
!------------------------------------------------------------------------------
 
! zpi has to be set to 0 (at least at the boundary)
zpi(:,:,:) = 0.0_wp 

zbt = (1.0_wp + betasw) * dt

DO j = jstart-1, jend+1
   zfx(j)   = eddlon*acrlat(j,1)
   zfydn(j) = eddlat*acrlat(j,1)*crlat(j  ,2)
   zfyds(j) = eddlat*acrlat(j,1)*crlat(j-1,2)
ENDDO          

DO  k = 1, ke
  kp1 = MIN(k+1, ke)
  km1 = MAX(k-1, 1 )
  DO  j = jstart-1, jend
    DO  i = istart-1, iend
      zdpq1   =   p0(i  ,j  ,kp1) - p0(i  ,j  ,km1)
      zdpq2   =   p0(i+1,j  ,kp1) - p0(i+1,j,  km1)
      zdpq3   =   p0(i  ,j+1,kp1) - p0(i  ,j+1,km1)
      zp0u(i,j,k) =  (p0(i+1,j,k) - p0(i  ,j,k)) / (zdpq2 + zdpq1)
      zp0v(i,j,k) =  (p0(i,j+1,k) - p0(i  ,j,k)) / (zdpq3 + zdpq1)
    ENDDO
  ENDDO
ENDDO

! Top and bottom boundary conditions on zqws

DO j = jstart, jend
  DO  i = istart, iend
    zjpvn = ( hhl(i  ,j+1,ke1) - hhl(i  ,j  ,ke1) ) * zqv(i,j  ,ke)
    zjpvs = ( hhl(i  ,j  ,ke1) - hhl(i  ,j-1,ke1) ) * zqv(i,j-1,ke)
    zjlur = ( hhl(i+1,j  ,ke1) - hhl(i  ,j  ,ke1) ) * zqu(i  ,j,ke)
    zjlul = ( hhl(i  ,j  ,ke1) - hhl(i-1,j  ,ke1) ) * zqu(i-1,j,ke)
    zqws(i,j,ke1) = 0.5_wp*(  zfx(j)*( zjlur + zjlul ) +                        &
                           zfydn(j)*zjpvn + zfyds(j)*zjpvs )
    zqws(i,j,1)   = 0.0_wp
    zd(i,j,ke1)   = zqws(i,j,ke1)
    zd(i,j,1)     = zqws(i,j,1)
  ENDDO
ENDDO

! Compute zqws 

DO k = 1, ke
  DO j = jstart, jend
    DO  i = istart, iend
      zpn         = p0(i,j,k) + pp(i,j,k,nnow)
      zpb (i,j,k) = 0.5_wp * rdocp * zqp(i,j,k) / p0(i,j,k)
      ztp         = zqt(i,j,k)
      ztb (i,j,k) = 0.5_wp * ztp / t(i,j,k,nnow) * zpn / p0(i,j,k)
    ENDDO
  ENDDO
ENDDO

DO k = 2, ke
  DO j = jstart, jend
    DO  i = istart, iend
      zrofac = ( rho0(i,j,k)*dp0(i,j,k-1) + rho0(i,j,k-1)*dp0(i,j,k) ) /     &
               ( rho (i,j,k)*dp0(i,j,k-1) + rho (i,j,k-1)*dp0(i,j,k) )
      zcw1   = 2.0_wp * g * zrofac / ( dp0(i,j,k-1) + dp0(i,j,k) )
      zbuoy  = zcw1 * (                                                      &
              ( ztb(i,j,k)   - zpb(i,j,k)   ) * dp0(i,j,k-1) +               &
              ( ztb(i,j,k-1) - zpb(i,j,k-1) ) * dp0(i,j,k)   )
      zqws(i,j,k) = zqw(i,j,k) + zbt * zbuoy
    ENDDO
  ENDDO
ENDDO

! Compute horizontal divergence Dh( zqu, zqv )

! No communication since zqu and zqv available in halo

DO k = 1, ke
  kp1 = MIN( ke, k+1 )
  km1 = MAX( 1 , k-1 )
  zfact = 1.0_wp
  IF(k==1 .OR. k==ke) zfact=0.5_wp
  IF ( k < kzflat ) THEN
    DO j = jstart, jend
    DO i = istart, iend
      zpi(i,j,k) = zfx(j) * ( zqu(i,j,k) - zqu(i-1,j,k) ) &
                  + zfydn(j) * zqv(i,j,k) - zfyds(j) * zqv(i,j-1,k)
    ENDDO
    ENDDO
  ELSE
    DO j = jstart, jend
    DO i = istart, iend
      zpi(i,j,k) =                                                   &
       zfx(j)*( zqu(i,j,k) - zqu(i-1,j,k)                            &
        - zfact * (zqu(i-1,j,kp1) - zqu(i-1,j,km1))* zp0u(i-1,j,k)   &
        - zfact * (zqu(i  ,j,kp1) - zqu(i  ,j,km1))* zp0u(i  ,j,k) ) &
        + zfydn(j)*(  zqv(i,j ,k)                                    &
        - zfact * (zqv(i,j,kp1) - zqv(i,j,km1)) * zp0v(i,j,k) )      &
        + zfyds(j)*(- zqv(i,j-1,k)                                   &
        - zfact * (zqv(i,j-1,kp1) - zqv(i,j-1,km1)) * zp0v(i,j-1,k) )
    ENDDO
    ENDDO
  ENDIF
ENDDO      

! Compute zqpp (store in zqps)

DO k = 1, ke
  DO j = jstart, jend
    DO  i = istart, iend
      zpn          = p0(i,j,k) + pp(i,j,k,nnow)
      zpfac        = 1.0_wp / ( gamma * zpn )
      zqps(i,j,k)  = zpfac * zqp(i,j,k) - zbt * zpi(i,j,k)
    ENDDO
  ENDDO
ENDDO

! LU Decomposition N = L.U

DO j = jstart, jend
  DO  i = istart, iend
    zdenom    = zb(i,j,2)
    zb(i,j,2) = 1.0_wp / zdenom
    zc(i,j,2) = - zc(i,j,2) * zb(i,j,2)
  ENDDO
ENDDO

DO k = 3, ke
  DO j = jstart, jend
    DO  i = istart, iend
      zdenom    = za(i,j,k) * zc(i,j,k-1) + zb(i,j,k)
      zb(i,j,k) = 1.0_wp / zdenom
      zc(i,j,k) = - zc(i,j,k) * zb(i,j,k)
    ENDDO
  ENDDO
ENDDO

! Tridiagonal solve L.U.zd = zqws

DO k = 2, ke
  DO j = jstart, jend
    DO  i = istart, iend
      zd(i,j,k) = ( zqws(i,j,k) - za(i,j,k) * zd(i,j,k-1) ) * zb(i,j,k)
!     zd(i,j,k) = zqws(i,j,k) * zb(i,j,k)
    ENDDO
  ENDDO
ENDDO

DO k = ke, 2, -1
  DO j = jstart, jend
    DO  i = istart, iend
      zd(i,j,k) = zc(i,j,k) * zd(i,j,k+1) + zd(i,j,k)
    ENDDO
  ENDDO
ENDDO

! Compute zqps = zqpp - dt.dz2.Ni.zqws

DO k = 1, ke
  DO j = jstart, jend
    DO  i = istart, iend
      zqps(i,j,k) = zqps(i,j,k) - zbt * (                             &
                    zb2(i,j,k)*zd(i,j,k+1) - za2(i,j,k)*zd(i,j,k) )
    ENDDO
  ENDDO
ENDDO


!------------------------------------------------------------------------------
!  End of internal procedure wave_rhs
!------------------------------------------------------------------------------

END SUBROUTINE wave_rhs

!==============================================================================
!+ Internal procedure in "semi_implicit" for wave equation solver
!------------------------------------------------------------------------------

SUBROUTINE wave_solve

!------------------------------------------------------------------------------
!
! Description:
!   This internal procedure of "semi_implicit" serves as the
!   driver to compute the pressure perturbation tendency.
!
! Method:
!   The GMRES iterative Krylov solver is invoked with a
!   reverse communication strategy. When the solver requires
!   application of the elliptic operator (matvec) or the
!   preconditioner (precond), the GMRES subroutine saves
!   the current state and exits to wave_solve to invoke
!   the operation according to `icode'. Since GMRES 
!   works with 1D arrays we must copy between 3D halo'ed
!   arrays (ie,je,ke) and non-halo'ed (nn).
!
!------------------------------------------------------------------------------

! Local scalars:
! -----------------------------
  INTEGER (KIND=iintegers) ::  &
    i,  j,  k,           & !  Loop indices in longitudinal, latitudinal and
                           !  vertical direction
    kp1, km1,            & !  k+1, k-1
    ilowu, ilowv,        & !  start- and end-indices for subdomains with
    jlowu, jlowv,        & !  or without neighbours
    iup,   jup,          & !
    ik,                  & ! GMRES iteration
    icode                  ! precond (icode=1), matvec (icode>=2)

  REAL    (KIND=wp   )     ::  &
    zfy,                 & !
    zrhoqx, zrhoqy,      & !
    zdpq1, zdpq2, zdpq3    !

! Local (automatic) arrays:
! -----------------------------
  REAL    (KIND=wp   )     ::  &
    zpi     (ie,je,ke) ,       & ! 
    zfx     (je   )     ,      & !
    zfydn   (je   )    ,       & !
    zfyds   (je   )    ,       & !
    zp0u    (ie,je,ke) ,       & !
    zp0v    (ie,je,ke) ,       & !
    zrhoqx_i(ie,je,ke) ,       & !
    zrhoqy_i(ie,je,ke) ,       & !
    wk1     (ie,je,ke) ,       & !
    wk2     (ie,je,ke) ,       & !
    vv      (ie,je,ke,im1),    & !
    w       (ie,je,ke,im)        !
!
! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine wave_solve
!------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------
!  Section 1: Setup of subdomain indices to compute pressure gradients
!------------------------------------------------------------------------------

! Write headline in file YUSOLVER
IF ( (iprint_si > 0) .AND. (my_cart_id == 0) ) THEN
  WRITE (nusolver, '(A,I7)') '   Convergence in step ', ntstep
ENDIF

! Compute local subdomain start and end indices

IF (my_cart_neigh(1) == -1) THEN ! west
  ilowu = istartu
  ilowv = istartv
ELSE
  ilowu = istartu-1
  ilowv = istartv-1
ENDIF

IF (my_cart_neigh(2) == -1) THEN ! north
  jup = jend
ELSE
  jup = jend+1
ENDIF

IF (my_cart_neigh(3) == -1) THEN ! east
  iup = iend
ELSE
  iup = iend+1
ENDIF

IF (my_cart_neigh(4) == -1) THEN ! south
  jlowu = jstartu
  jlowv = jstartv
ELSE
  jlowu = jstartu-1
  jlowv = jstartv-1
ENDIF

!------------------------------------------------------------------------------
!  Section 2: Some preparations to avoid division
!------------------------------------------------------------------------------

zfy = edadlat
DO j = jstart-1, jend+1
   zfx(j)   = eddlon*acrlat(j,1)
   zfydn(j) = eddlat*acrlat(j,1)*crlat(j  ,2)
   zfyds(j) = eddlat*acrlat(j,1)*crlat(j-1,2)
ENDDO          

! Avoid division by zrhoqx and zrhoqy 

DO  k = 1, ke
  DO  j = jlowu, jendu
    DO  i = ilowu, iendu
      zrhoqx          = 0.5_wp*( rho(i,j,k) + rho(i+1,j,k) )
      zrhoqx_i(i,j,k) = zfx(j) / zrhoqx
    ENDDO
    IF ( j <= jendv ) THEN
      DO  i = ilowv, iendv
        zrhoqy          = 0.5_wp*( rho(i,j,k) + rho(i,j+1,k) )
        zrhoqy_i(i,j,k) = zfy / zrhoqy
      ENDDO
    ENDIF
  ENDDO
ENDDO

DO  k = 1, ke
  kp1 = MIN(k+1, ke)
  km1 = MAX(k-1, 1 )
  DO  j = jstart-1, jend
    DO  i = istart-1, iend
      zdpq1   =   p0(i  ,j  ,kp1) - p0(i  ,j  ,km1)
      zdpq2   =   p0(i+1,j  ,kp1) - p0(i+1,j,  km1)
      zdpq3   =   p0(i  ,j+1,kp1) - p0(i  ,j+1,km1)
      zp0u(i,j,k) =  (p0(i+1,j,k) - p0(i  ,j,k)) / (zdpq2 + zdpq1)
      zp0v(i,j,k) =  (p0(i,j+1,k) - p0(i  ,j,k)) / (zdpq3 + zdpq1)
    ENDDO
  ENDDO
ENDDO

icode = 0

! Initial guess for the the pressure tendency 
zpi(:,:,:) = 0.0_wp 
wk1(:,:,:) = 0.0_wp 
wk2(:,:,:) = 0.0_wp 

1 CONTINUE

  CALL fgmres(ie, je, ke, im, zqps, zpi, ik, vv, w, wk1, wk2, eps_si,   &
              maxit_si, icode)

IF ( icode == 1 ) THEN

  CALL precond ( wk2, wk1, ie, je, ke )

ELSEIF ( icode >= 2 ) THEN

  CALL matvec ( wk2, wk1, zfx, zfydn, zfyds, zp0u, zp0v, &
                zrhoqx_i, zrhoqy_i, ie, je, ke )

ENDIF

IF ( (icode == 1) .OR. (icode >= 2) ) GOTO 1

! Compute the new pressure pp(nnew)

DO  k = 1, ke
  pp(istart:iend,jstart:jend,k,nnew) =                     &
    pp(istart:iend,jstart:jend,k,nold) +                   &
    zpi(istart:iend,jstart:jend,k) * dt2
!   zpi(istart:iend,jstart:jend,k)
ENDDO         

!------------------------------------------------------------------------------
!  End of internal procedure wave_solve
!------------------------------------------------------------------------------

END SUBROUTINE wave_solve

!==============================================================================
!+ Internal procedure in "semi_implicit" for updating prognostic variables
!------------------------------------------------------------------------------

SUBROUTINE back_solve

!------------------------------------------------------------------------------
!
! Description:
!   This internal procedure of "semi_implicit" uses the pressure
!   perturbation tendency obtained from the GMRES solver to
!   update the prognostic variables to time level n+1.
!
! Method:
!   Once the pressure tendency has been computed in "wave_solve"
!   the implicit equations for the remaining prognostic variables
!   can be solved. A tridiagonal solver is required in order to 
!   invert the Brunt-Vaisalla matrix N in the computation of the 
!   vertical velocity tendency.
!
!------------------------------------------------------------------------------

! Local scalars:
! -----------------------------
  INTEGER (KIND=iintegers) ::  &
    i,  j,  k,           & !  Loop indices in longitudinal, latitudinal and
                           !  vertical direction 
    kp1, km1,            & !  k+1, k-1
    ilowu, ilowv,        & !  start- and end-indices for subdomains with
    jlowu, jlowv,        & !  or without neighbours
    iup,   jup             !  

  REAL    (KIND=wp   )     ::  &
    zbt,                 & !
    zfy,                 & !
    zrhoqx, zrhoqy,      & !
    zpdx, zpdy,          & !
    zdppz,               & !
    zdpq1, zdpq2, zdpq3, & !
    zppqx, zppqy,        & !
    zpgradx, zpgrady,    & !
    zct1, zfact,         & !
           ztfac,        & !
    zjpvn, zjpvs,        & !
    zjlur, zjlul           !

   
! Local (automatic) arrays:
! -----------------------------
  REAL    (KIND=wp   )     ::  &
    zppqz   (ie,je)     ,      & !
    zfx     (je   )     ,      & !
    zfydn   (je   )     ,      & !
    zfyds   (je   )     ,      & !
    zpi     (ie,je,ke)  ,      & !
    zp0u    (ie,je,ke)  ,      & !
    zp0v    (ie,je,ke)  ,      & !
    zrhoqx_i(ie,je,ke)  ,      & !
    zrhoqy_i(ie,je,ke)  ,      & !
    zdu     (ie,je,ke)  ,      & !
    zdv     (ie,je,ke)  ,      & !
    zdw     (ie,je,ke1) ,      & !
    zd      (ie,je,ke1) ,      & !
    zdt     (ie,je,ke)  ,      & !
    zu      (ie,je,ke)  ,      & !
    zv      (ie,je,ke)  ,      & !
    zw      (ie,je,ke1) ,      & !
    zt      (ie,je,ke)           !

! For error handling
! ------------------
  INTEGER (KIND=iintegers) ::  &
    izerror

  CHARACTER (LEN=80)       ::  &
    yzerrmsg
!
! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine back_solve
!------------------------------------------------------------------------------
 
izerror  = 0
yzerrmsg = '   '

!------------------------------------------------------------------------------
!  Section 1: Setup of subdomain indices to compute pressure gradients
!------------------------------------------------------------------------------

! zpi has to be set to 0 (at least at the boundary)
zpi(:,:,:) = 0.0_wp 

zu(:,:,:)  = u(:,:,:,nold)
zv(:,:,:)  = v(:,:,:,nold)
zw(:,:,:)  = w(:,:,:,nold)
zt(:,:,:)  = t(:,:,:,nold)

zbt = (1.0_wp + betasw) * dt

! Compute local subdomain start and end indices

IF (my_cart_neigh(1) == -1) THEN ! west
  ilowu = istartu
  ilowv = istartv
ELSE
  ilowu = istartu-1
  ilowv = istartv-1
ENDIF

IF (my_cart_neigh(2) == -1) THEN ! north
  jup = jend
ELSE
  jup = jend+1
ENDIF

IF (my_cart_neigh(3) == -1) THEN ! east
  iup = iend
ELSE
  iup = iend+1
ENDIF

IF (my_cart_neigh(4) == -1) THEN ! south
  jlowu = jstartu
  jlowv = jstartv
ELSE
  jlowu = jstartu-1
  jlowv = jstartv-1
ENDIF

!-------------------------------------------------------------------------------
!  Section 2: Some preparations to avoid division
!-------------------------------------------------------------------------------

zfy = edadlat
DO j = jstart-1, jend+1
   zfx(j)   = eddlon*acrlat(j,1)
   zfydn(j) = eddlat*acrlat(j,1)*crlat(j  ,2)
   zfyds(j) = eddlat*acrlat(j,1)*crlat(j-1,2)
ENDDO          

! Avoid division by zrhoqx and zrhoqy 

DO  k = 1, ke
  DO  j = jlowu, jendu
    DO  i = ilowu, iendu
      zrhoqx          = 0.5_wp*( rho(i,j,k) + rho(i+1,j,k) )
      zrhoqx_i(i,j,k) = zfx(j) / zrhoqx
    ENDDO
    IF ( j <= jendv ) THEN
      DO  i = ilowv, iendv
        zrhoqy          = 0.5_wp*( rho(i,j,k) + rho(i,j+1,k) )
        zrhoqy_i(i,j,k) = zfy / zrhoqy
      ENDDO
    ENDIF
  ENDDO
ENDDO

DO  k = 1, ke
  kp1 = MIN(k+1, ke)
  km1 = MAX(k-1, 1 )
  DO  j = jstart-1, jend
    DO  i = istart-1, iend
      zdpq1   =   p0(i  ,j  ,kp1) - p0(i  ,j  ,km1)
      zdpq2   =   p0(i+1,j  ,kp1) - p0(i+1,j,  km1)
      zdpq3   =   p0(i  ,j+1,kp1) - p0(i  ,j+1,km1)
      zp0u(i,j,k) =  (p0(i+1,j,k) - p0(i  ,j,k)) / (zdpq2 + zdpq1)
      zp0v(i,j,k) =  (p0(i,j+1,k) - p0(i  ,j,k)) / (zdpq3 + zdpq1)
    ENDDO
  ENDDO
ENDDO

! Save pressure perturbation tendency ( pp(nnew) - pp(nold) ) / dt2

DO  k = 1, ke
  zpi(istart-1:iend+1,jstart-1:jend+1,k) =                              &
  ( pp(istart-1:iend+1,jstart-1:jend+1,k,nnew) -                        &
    pp(istart-1:iend+1,jstart-1:jend+1,k,nold) ) / dt2
ENDDO         

! Lateral boundary conditions on velocity tendencies zdu and zdv

IF (my_cart_neigh(1) == -1) THEN  ! west
  DO k = 1, ke
    DO j = jstartu, jendu
      zdu(istartu-1,j,k) = zqu(istartu-1,j,k)
      zu (istartu-1,j,k) = u(istartu-1,j,k,nnew)
    ENDDO
  ENDDO
END IF

IF (my_cart_neigh(3) == -1) THEN  ! east
  DO k = 1, ke
    DO j = jstartu, jendu
      zdu(iendu+1,j,k) = zqu(iendu+1,j,k)
      zu (iendu+1,j,k) = u(iendu+1,j,k,nnew)
    ENDDO
  ENDDO
END IF

IF (my_cart_neigh(4) == -1) THEN  ! south
  DO k = 1, ke
    DO i = istartv, iendv
      zdv(i,jstartv-1,k) = zqv(i,jstartv-1,k)
      zv (i,jstartv-1,k) = v(i,jstartv-1,k,nnew)
    ENDDO
  ENDDO
END IF

IF (my_cart_neigh(2) == -1) THEN  ! north
  DO k = 1, ke
    DO i = istartv, iendv
      zdv(i,jendv+1,k) = zqv(i,jendv+1,k)
      zv (i,jendv+1,k) = v(i,jendv+1,k,nnew)
    ENDDO
  ENDDO
END IF

! Compute horizontal pressure gradients zdu and zdv 
IF (ltime) THEN
  CALL get_timings (i_semi_implicit, ntstep, dt, izerror)
  IF (ltime_barrier) THEN
    CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
    CALL get_timings (i_semi_impl_barrier, ntstep, dt, izerror)
  ENDIF
ENDIF

kzdims(1:24)=(/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
CALL exchg_boundaries                                                    &
   (19,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,       &
    ie, je, kzdims, jstartpar, jendpar, 1, nboundlines,                  &
    my_cart_neigh, lperi_x, lperi_y, l2dim,                              &
    10000+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,          &
    zpi(:,:,:) )
IF (ltime) CALL get_timings (i_semi_impl_comm, ntstep, dt, izerror)

! Compute tendencies zdu and zdv. Compute u(nnew) and v(nnew). 

DO  k = 1, ke

  kp1   = MIN ( ke, k+1 )
  km1   = MAX ( 1 , k-1 )

  IF ( k < kzflat ) THEN
    DO  j = jlowu, jendu
      DO  i = ilowu, iendu   
        zpgradx =  (zpi(i+1,j,k) - zpi(i,j,k)) * zrhoqx_i(i,j,k)
        zdu(i,j,k) = zqu(i,j,k) - zbt * zpgradx
        zu(i,j,k)  = u(i,j,k,nold) + dt2*zdu(i,j,k)
      ENDDO 
      IF ( j <= jendv ) THEN
        DO  i = ilowv, iendv   
          zpgrady =  (zpi(i,j+1,k) - zpi(i,j,k)) * zrhoqy_i(i,j,k)
          zdv(i,j,k) = zqv(i,j,k) - zbt * zpgrady
          zv(i,j,k)  = v(i,j,k,nold) + dt2*zdv(i,j,k)
        ENDDO 
      ENDIF
    ENDDO   
  ELSE 
    zppqz(istart-1:iend+1,jstart-1:jend+1) =                  &
                  zpi(istart-1:iend+1,jstart-1:jend+1,kp1)    &
                 -zpi(istart-1:iend+1,jstart-1:jend+1,km1)

    DO  j = jlowu, jendu
      DO  i = ilowu, iendu   
        zppqx   =  zppqz(i+1,j) + zppqz(i,j)
        zdppz   =  zppqx * zp0u(i,j,k)
        zpdx    =  zpi(i+1,j,k) - zpi(i,j,k)
        zpgradx =  zpdx - zdppz
        zdu(i,j,k) = zqu(i,j,k) - zbt * zpgradx*zrhoqx_i(i,j,k)
        zu(i,j,k)  = u(i,j,k,nold) + dt2*zdu(i,j,k)
      ENDDO 
      IF ( j <= jendv ) THEN
        DO  i = ilowv, iendv   
          zppqy   =  zppqz(i,j+1) + zppqz(i,j)
          zdppz   =  zppqy * zp0v(i,j,k)
          zpdy    =  zpi(i,j+1,k) - zpi(i,j,k)
          zpgrady =  zpdy - zdppz
          zdv(i,j,k) = zqv(i,j,k) - zbt * zpgrady*zrhoqy_i(i,j,k)
          zv(i,j,k)  = v(i,j,k,nold) + dt2*zdv(i,j,k)
        ENDDO 
      ENDIF
    ENDDO   
  ENDIF

ENDDO

! Top and bottom boundary conditions 

DO j = jstart, jend
  DO  i = istart, iend
    zjpvn = ( hhl(i  ,j+1,ke1) - hhl(i  ,j  ,ke1) ) * zv(i,j  ,ke)
    zjpvs = ( hhl(i  ,j  ,ke1) - hhl(i  ,j-1,ke1) ) * zv(i,j-1,ke)
    zjlur = ( hhl(i+1,j  ,ke1) - hhl(i  ,j  ,ke1) ) * zu(i  ,j,ke)
    zjlul = ( hhl(i  ,j  ,ke1) - hhl(i-1,j  ,ke1) ) * zu(i-1,j,ke)
    zdw(i,j,ke1) = 0.5_wp*(  zfx(j)*( zjlur + zjlul ) +                        &
                          zfydn(j)*zjpvn + zfyds(j)*zjpvs )
    zdw(i,j,1)   = 0.0_wp
    zd(i,j,ke1)  = ( zdw(i,j,ke1) - w(i,j,ke1,nold) ) / dt2
    zdw(i,j,ke1) =  zd(i,j,ke1)
    zd(i,j,1)    = 0.0_wp
    zw(i,j,ke1)  = zdw(i,j,ke1)
    zw(i,j,1)    = zdw(i,j,1)
  ENDDO
ENDDO

! Compute zdw = qws - dt * Dz1.d(pp)

DO k = 2, ke
  DO j = jstart, jend
    DO  i = istart, iend
      zdw(i,j,k) = zqws(i,j,k) -                                            &
                   zbt * ( zb1(i,j,k)*zpi(i,j,k) - za1(i,j,k)*zpi(i,j,k-1) )
    ENDDO
  ENDDO
ENDDO

! Tridiagonal solve L.U.zd = zdw for vertical velocity tendency

DO k = 2, ke
  DO j = jstart, jend
    DO  i = istart, iend
      zd(i,j,k) = ( zdw(i,j,k) - za(i,j,k) * zd(i,j,k-1) ) * zb(i,j,k)
!     zd(i,j,k) = zdw(i,j,k) * zb(i,j,k)
    ENDDO
  ENDDO
ENDDO

DO k = ke, 2, -1
  DO j = jstart, jend
    DO  i = istart, iend
      zd(i,j,k) = zc(i,j,k) * zd(i,j,k+1) + zd(i,j,k)
    ENDDO
  ENDDO
ENDDO

! Compute horizontal divergence Dh( zdu, zdv )

! No communication since zdu and zdv available in halo

DO k = 1, ke
  kp1 = MIN( ke, k+1 )
  km1 = MAX( 1 , k-1 )
  zfact = 1.0_wp
  IF(k==1 .OR. k==ke) zfact=0.5_wp
  IF ( k < kzflat ) THEN
    DO j = jstart, jend
    DO i = istart, iend
      zpi(i,j,k) = zfx(j)*( zdu(i,j,k) - zdu(i-1,j,k) ) &
                  + zfydn(j)*zdv(i,j,k) - zfyds(j)*zdv(i,j-1,k)
    ENDDO
    ENDDO
  ELSE
    DO j = jstart, jend
    DO i = istart, iend
      zpi(i,j,k) =                                                         &
       zfx(j)*( zdu(i,j,k) - zdu(i-1,j,k)                                  &
        - zfact * (zdu(i-1,j,kp1) - zdu(i-1,j,km1))* zp0u(i-1,j,k)         &
        - zfact * (zdu(i  ,j,kp1) - zdu(i  ,j,km1))* zp0u(i  ,j,k) )       &
        + zfydn(j)*(  zdv(i,j  ,k)                                         &
        - zfact * (zdv(i,j,kp1)   - zdv(i,j,km1))  * zp0v(i,j,k) )         &
        + zfyds(j)*(- zdv(i,j-1,k)                                         &
        - zfact * (zdv(i,j-1,kp1) - zdv(i,j-1,km1)) * zp0v(i,j-1,k) )
    ENDDO
    ENDDO
  ENDIF
ENDDO      

! Compute w(nnew)

DO k = 1, ke1
  DO j = jstart, jend
    DO  i = istart, iend
      zw(i,j,k) = w(i,j,k,nold) + dt2*zd(i,j,k)
    ENDDO
  ENDDO
ENDDO 

! Compute zdt (temperature)

DO k = 1, ke
  DO j = jstart, jend
    DO  i = istart, iend
      ztfac      = gamma * cpdr * r_d * t(i,j,k,nnow)
      zct1       = ztfac * g * rho0(i,j,k) / dp0(i,j,k)
      zdt(i,j,k) = zqt(i,j,k) - zbt * ztfac * zpi(i,j,k) +                 &
                   zbt * zct1 * ( zd(i,j,k+1) - zd(i,j,k) )
    ENDDO
  ENDDO
ENDDO

! Compute tw(nnew)

DO k = 1, ke
  DO j = jstart, jend
    DO  i = istart, iend
      zt(i,j,k) = t(i,j,k,nold) + dt2 * zdt(i,j,k)
    ENDDO
  ENDDO
ENDDO 

! u(:,:,:,nnew) = zu
! v(:,:,:,nnew) = zv
! w(:,:,:,nnew) = zw
! t(:,:,:,nnew) = zt

DO  k = 1, ke
  DO    j = jstartu  , jendu
    DO  i = istartu  , iendu
        u (i,j,k,nnew)  = zu (i,j,k)
    ENDDO
  ENDDO
  DO    j = jstartv  , jendv
    DO  i = istartv  , iendv
        v (i,j,k,nnew)  = zv (i,j,k)
    ENDDO
  ENDDO
  DO    j = jstart, jend
    DO  i = istart, iend
        w (i,j,k,nnew)  = zw (i,j,k)
        t (i,j,k,nnew)  = zt (i,j,k)
    ENDDO
  ENDDO
!  DO    j = jlowv, jup
!    DO  i = ilowu, iup
!        pp(i,j,k,nnew)  = pp(i,j,k)
!    ENDDO
!  ENDDO
ENDDO
w(istart:iend,jstart:jend,ke1,nnew ) = zw(istart:iend,jstart:jend,ke1)

!  Set free-slip lateral boundary conditions on w or
!  periodic boundary conditions on all variables if required

IF ( lw_freeslip ) THEN
  IF (my_cart_neigh(1) == -1) THEN ! west
    DO k = 1, ke1
      DO j = jstart, jend
        DO i = 1, nboundlines
          w(i,j,k,nnew) = w(istart,j,k,nnew)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  IF (my_cart_neigh(3) == -1) THEN ! east
    DO k = 1, ke1
      DO j = jstart, jend
        DO i = ie-nboundlines+1, ie
          w(i,j,k,nnew) = w(iend  ,j,k,nnew)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  IF (my_cart_neigh(4) == -1) THEN ! south
    DO k = 1, ke1
      DO j = 1, nboundlines
        DO i = 1,ie
          w(i,j,k,nnew) = w(i,jstart,k,nnew)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  IF (my_cart_neigh(2) == -1) THEN ! north
    DO k = 1, ke1
      DO j = je-nboundlines+1, je
        DO i = 1,ie
          w(i,j,k,nnew) = w(i,jend  ,k,nnew)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
ENDIF

!------------------------------------------------------------------------------
!  End of internal procedure back_solve
!------------------------------------------------------------------------------

END SUBROUTINE back_solve

!==============================================================================
!+ Internal procedure in "semi_implicit" for the elliptic solver
!------------------------------------------------------------------------------

SUBROUTINE fgmres (ie, je, ke, im, rhs, sol, i, vv, w, wk1, wk2,            &
                   eps, maxits, icode)

!------------------------------------------------------------------------------
!
! Description:
! flexible GMRES routine which allows a variable preconditioner. 
! Implemented with a reverse communication protocol for flexibility -
! explicit (exact) residual norms for restarts  
!
! written by Y. Saad, modified by A. Malevsky, version February 1, 1995
! (Reverse Communication Implementation)
!
!   USAGE: (see also comments for icode below). CGMRES
!   should be put in a loop and the loop should be active for as
!   long as icode is not equal to 0. On return fgmres will
!      1) either be requesting the new preconditioned vector applied
!         to wk1 in case icode.eq.1 (result should be put in wk2) 
!      2) or be requesting the product of A applied to the vector wk1
!         in case icode.eq.2 (result should be put in wk2) 
!      3) or be terminated in case icode .eq. 0. 
!   on entry always set icode = 0. So icode should be set back to zero
!   upon convergence.
!
!   Here is a typical way of running fgmres: 
!
!        icode = 0
!   1    continue
!        call fgmres (n,im,rhs,sol,i,vv,w,wk1,wk2,eps,maxits,icode)
!
!        if (icode .eq. 1) then
!           call  precon(n, wk1, wk2)    <--- user's variable preconditioning
!           goto 1
!        else if (icode .ge. 2) then
!           call  matvec (n,wk1, wk2)    <--- user's matrix vector product. 
!           goto 1
!        else 
!           ----- done ---- 
!           .........
!
!-------- INPUT --------
!
! ie,je,ke == integer. the dimension of the problem
!
! im    == size of Krylov subspace:  should not exceed 50 in this
!          version (can be reset in code. looking at comment below)
!
! rhs   == vector of length n containing the right hand side
!
! sol   == initial guess on input, approximate solution on output
!
! vv    == work space of size n x (im+1)
!
! w     == work space of length n x im 
!
! wk1,
! wk2,  == two work vectors of length n each used for the reverse
!          communication protocole. When on return (icode .ne. 1)
!          the user should call fgmres again with wk2 = precon * wk1
!          and icode untouched. When icode.eq.1 then it means that
!          convergence has taken place.
!          
! eps   == tolerance for stopping criterion. process is stopped
!          as soon as ( ||.|| is the euclidean norm):
!          || current residual||/||initial residual|| <= eps
!
! maxits== maximum number of iterations allowed
!
! icode = integer. indicator for the reverse communication protocole.
!         ON ENTRY : icode should be set to icode = 0.
!         ON RETURN: 
!       * icode .eq. 1 value means that fgmres has not finished
!         and that it is requesting a preconditioned vector before
!         continuing. The user must compute M**(-1) wk1, where M is
!         the preconditioing  matrix (may vary at each call) and wk1 is
!         the vector as provided by fgmres upun return, and put the 
!         result in wk2. Then fgmres must be called again without
!         changing any other argument. 
!       * icode .eq. 2 value means that fgmres has not finished
!         and that it is requesting a matrix vector product before
!         continuing. The user must compute  A * wk1, where A is the
!         coefficient  matrix and wk1 is the vector provided by 
!         upon return. The result of the operation is to be put in
!         the vector wk2. Then fgmres must be called again without
!         changing any other argument. 
!       * icode .eq. 0 means that fgmres has finished and sol contains 
!         the approximate solution.
!         comment: typically fgmres must be implemented in a loop
!         with fgmres being called as long icode is returned with 
!         a value .ne. 0. 
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------
  INTEGER (KIND=iintegers), INTENT (IN)    ::  &
    ie, je, ke,   & ! the dimension of the problem
    im,           & ! size of Krylov subspace
    maxits          ! maximum number of iterations allowed

  INTEGER (KIND=iintegers), INTENT (INOUT) ::  &
    icode           ! indicator for the reverse communication protocole

  REAL    (KIND=wp   ),     INTENT (IN)    ::  &
    rhs(ie,je,ke),& ! vector of length n containing the right hand side
    wk2(ie,je,ke),& ! work vector of length n
    eps             ! tolerance for stopping criterion

  REAL    (KIND=wp   ),     INTENT (INOUT) ::  &
    w(ie,je,ke,im),   & ! work space of length n x im
    wk1(ie,je,ke),    & ! work vector of length n
    vv(ie,je,ke,im+1),& ! work space of size n x (im+1)
    sol(ie,je,ke)       ! initial guess on input, approximate solution on output


! Local variables:
! ----------------

  INTEGER (KIND=iintegers) ::  &
    i,j,k,i1,k1,its,ii,jj, i2, j2, k2, izerror

  REAL    (KIND=wp   )     ::  &
    hh(51  ,50), c(50), s(50), rs(51), t, ro, eps1, r0, gam,    &
    epsmac=1.0E-6_wp

  CHARACTER (LEN=80)       ::  yzerrmsg

SAVE
!
! End of header
!==============================================================================

      izerror  = 0
      yzerrmsg = '   '

!------------------------------------------------------------------------------
! Begin Subroutine fgmres
!------------------------------------------------------------------------------

!     computed goto 
!     GOTO (100,200,300,11) icode +1   ! obsolete Fortran Feature
      IF (icode == 0) GOTO 100
      IF (icode == 1) GOTO 200
      IF (icode == 2) GOTO 300
      IF (icode == 3) GOTO 11

 100  CONTINUE

      its = 0

!     outer loop starts here..

!--------------compute initial residual vector --------------

      DO k2 = 1, ke
        DO j2 = jstart, jend
          DO i2 = istart, iend
            wk1(i2,j2,k2) = sol(i2,j2,k2)
          ENDDO
        ENDDO
      ENDDO

      icode = 3
      GOTO 9991    ! return

 11   CONTINUE

      DO k2 = 1, ke
        DO j2 = jstart, jend
          DO i2 = istart, iend
            vv(i2,j2,k2,1) = rhs(i2,j2,k2) - wk2(i2,j2,k2)
          ENDDO
        ENDDO
      ENDDO

 20   CONTINUE

      ! compute the parallel dot-product vv(:,1)*vv(:,1))
      ro = 0.0_wp
      DO k2 = 1, ke
        DO j2 = jstart, jend
          DO i2 = istart, iend
            ro = ro + vv(i2,j2,k2,1) * vv(i2,j2,k2,1)
          ENDDO
        ENDDO
      ENDDO
      IF (num_compute > 1) THEN
        IF (ltime) THEN
          CALL get_timings (i_semi_implicit, ntstep, dt, izerror)
          CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
          CALL get_timings (i_semi_impl_barrier, ntstep, dt, izerror)
        ENDIF
        CALL global_values (ro, 1, 'SUM', imp_reals, icomm_cart, -1,   &
                            yzerrmsg, izerror)
        IF (ltime) CALL get_timings (i_semi_impl_comm, ntstep, dt, izerror)
      ENDIF
      ro = SQRT(ro)
      IF (ro == 0.0_wp) GOTO 999 

      t = 1.0_wp/ ro
      DO k2 = 1, ke
        DO j2 = jstart, jend
          DO i2 = istart, iend
            vv(i2,j2,k2,1) = vv(i2,j2,k2,1) * t
          ENDDO
        ENDDO
      ENDDO

      IF (its .eq. 0) eps1=eps
      IF (its .eq. 0) r0 = ro
      IF ( (iprint_si > 0) .AND. (my_cart_id == 0) ) THEN
        IF (MOD(its,10) == 0) THEN
          WRITE (nusolver, '(A,I4,A,F20.15)')                         &
                 '   -- fgmres its = ', its, ' res. norm = ', ro
        ENDIF
      ELSEIF ( (iprint_si < 0) .AND. (my_cart_id == 0) ) THEN
        PRINT *, '   -- fgmres its = ', its, ' res. norm = ', ro
      ENDIF

!     initialize 1-st term  of rhs of hessenberg system..
      rs(1) = ro
      i = 0
 4    i=i+1
      its = its + 1
      i1 = i + 1
      DO k2 = 1, ke
        DO j2 = jstart, jend
          DO i2 = istart, iend
            wk1 (i2,j2,k2) = vv(i2,j2,k2,i)
          ENDDO
        ENDDO
      ENDDO

      icode = 1
      GOTO 9991   ! return

 200  CONTINUE

      DO k2 = 1, ke
        DO j2 = jstart, jend
          DO i2 = istart, iend
            w(i2,j2,k2,i) = wk2(i2,j2,k2)
          ENDDO
        ENDDO
      ENDDO

!     call matvec operation
      icode = 2
      DO k2 = 1, ke
        DO j2 = jstart, jend
          DO i2 = istart, iend
            wk1(i2,j2,k2) = wk2(i2,j2,k2)
          ENDDO
        ENDDO
      ENDDO

      GOTO 9991   ! return


 300  CONTINUE

!     first call to ope corresponds to intialization goto back to 11.
      DO k2 = 1, ke
        DO j2 = jstart, jend
          DO i2 = istart, iend
            vv(i2,j2,k2,i1) = wk2(i2,j2,k2)
          ENDDO
        ENDDO
      ENDDO

!     modified gram - schmidt...
      DO j=1, i
        ! compute the parallel dot-product vv(:,j)*vv(:,i1))
        ro = 0.0_wp
        DO k2 = 1, ke
          DO j2 = jstart, jend
            DO i2 = istart, iend
              ro = ro + vv(i2,j2,k2,j) * vv(i2,j2,k2,i1)
            ENDDO
          ENDDO
        ENDDO
        IF (num_compute > 1) THEN
          IF (ltime) THEN
            CALL get_timings (i_semi_implicit, ntstep, dt, izerror)
            CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
            CALL get_timings (i_semi_impl_barrier, ntstep, dt, izerror)
          ENDIF
          CALL global_values (ro, 1, 'SUM', imp_reals, icomm_cart, -1,   &
                              yzerrmsg, izerror)
          IF (ltime) CALL get_timings (i_semi_impl_comm, ntstep, dt, izerror)
        ENDIF
        t = ro
        hh(j,i) = t

        DO k2 = 1, ke
          DO j2 = jstart, jend
            DO i2 = istart, iend
              vv(i2,j2,k2,i1) = vv(i2,j2,k2,i1) - t * vv(i2,j2,k2,j)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      ! compute the parallel dot-product vv(:,i1)*vv(:,i1))
      ro = 0.0_wp
      DO k2 = 1, ke
        DO j2 = jstart, jend
          DO i2 = istart, iend
            ro = ro + vv(i2,j2,k2,i1) * vv(i2,j2,k2,i1)
          ENDDO
        ENDDO
      ENDDO
      IF (num_compute > 1) THEN
        IF (ltime) THEN
          CALL get_timings (i_semi_implicit, ntstep, dt, izerror)
          CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
          CALL get_timings (i_semi_impl_barrier, ntstep, dt, izerror)
        ENDIF
        CALL global_values (ro, 1, 'SUM', imp_reals, icomm_cart, -1,   &
                            yzerrmsg, izerror)
        IF (ltime) CALL get_timings (i_semi_impl_comm, ntstep, dt, izerror)
      ENDIF
      t = SQRT (ro)
      hh(i1,i) = t

      IF (t == 0.0_wp) GOTO 58
      t = 1.0_wp / t
      DO k2 = 1, ke
        DO j2 = jstart, jend
          DO i2 = istart, iend
            vv(i2,j2,k2,i1) = vv(i2,j2,k2,i1) * t
          ENDDO
        ENDDO
      ENDDO

!     done with modified gram schimdt and arnoldi step. 
!     now update factorization of hh

 58   IF (i == 1) GOTO 121

!     perform previous transformations on i-th column of h
      DO k=2,i
         k1 = k-1
         t = hh(k1,i)
         hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
         hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
      ENDDO

 121  gam = SQRT(hh(i,i)**2 + hh(i1,i)**2)
      IF (gam == 0.0_wp) gam = epsmac

!     determine next plane rotation
      c(i) = hh(i,i)/gam
      s(i) = hh(i1,i)/gam
      rs(i1) = -s(i)*rs(i)
      rs(i) =  c(i)*rs(i)

!     determine residual norm and test for convergence
      hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
      ro = ABS(rs(i1))

      IF ( (iprint_si > 0) .AND. (my_cart_id == 0) ) THEN
        IF (MOD(its,10) == 0) THEN
          WRITE (nusolver, '(A,I4,A,F20.15)')                         &
               '   -- fgmres its = ', its, ' res. norm = ', ro
        ENDIF
      ELSEIF ( (iprint_si < 0) .AND. (my_cart_id == 0) ) THEN
        PRINT *, '   -- fgmres its = ', its, ' res. norm = ', ro
      ENDIF
      IF (i < im .AND. (ro > eps1))  GOTO 4

!     now compute solution. first solve upper triangular system.
      rs(i) = rs(i)/hh(i,i)
      DO ii=2,i
         k=i-ii+1
         k1 = k+1
         t=rs(k)
         DO j=k1,i
            t = t-hh(k,j)*rs(j)
         ENDDO
         rs(k) = t/hh(k,k)
      ENDDO

!     done with back substitution..
!     now form linear combination to get solution
      DO j=1, i
        t = rs(j)
        DO k2 = 1, ke
          DO j2 = jstart, jend
            DO i2 = istart, iend
              sol(i2,j2,k2) = sol(i2,j2,k2) + t * w(i2,j2,k2,j)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!     test for return 
      IF (ro <= eps1 .OR. its >= maxits) GOTO 999

!     else compute residual vector and continue..
      DO j=1,i
         jj = i1-j+1
         rs(jj-1) = -s(jj-1)*rs(jj)
         rs(jj) = c(jj-1)*rs(jj)
      ENDDO

      DO j=1,i1
        t = rs(j)
        if (j == 1)  t = t-1.0_wp
        DO k2 = 1, ke
          DO j2 = jstart, jend
            DO i2 = istart, iend
              vv(i2,j2,k2,1) = vv(i2,j2,k2,1) + t * vv(i2,j2,k2,j)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!     restart outer loop.

      GOTO 20
 999  icode = 0

 9991 CONTINUE

END SUBROUTINE fgmres

!==============================================================================
!+ Internal procedure in "semi_implicit" for preconditioner
!------------------------------------------------------------------------------

SUBROUTINE precond ( zd, zw, ie, je, ke )

!------------------------------------------------------------------------------
!
! Description:
!   This internal procedure of "semi_implicit" applies
!   the vertical line Jacobi preconditioner.
!
! Method:
!   Perform one iteration of a Jacobi line relaxation
!   in the vertical direction. 
!
!------------------------------------------------------------------------------
! Subroutine arguments:
! --------------------
  INTEGER (KIND=iintegers), INTENT (IN   ) ::  &
    ie, je, ke            ! Array dimensions

  REAL    (KIND=wp   ),     INTENT (IN   ) ::  &
    zw(ie,je,ke)          ! Right-hand side

  REAL    (KIND=wp   ),     INTENT (  OUT) ::  &
    zd(ie,je,ke)          ! Solution

! Local scalars:
! -----------------------------
  INTEGER (KIND=iintegers) ::  &
    i,  j,  k              !  Loop indices in longitudinal, latitudinal and

! Local (automatic) arrays:
! -----------------------------
! REAL    (KIND=wp   )     ::  &
!
! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine precond
!------------------------------------------------------------------------------
 
! Tridiagonal solve L.U.zd = zw

DO j = jstart, jend
  DO  i = istart, iend
    zd(i,j,1) = zw(i,j,1) * zb3(i,j,1)
  ENDDO
ENDDO

DO k = 2, ke
  DO j = jstart, jend
    DO  i = istart, iend
      zd(i,j,k) = ( zw(i,j,k) - za3(i,j,k) * zd(i,j,k-1) ) * zb3(i,j,k)
    ENDDO
  ENDDO
ENDDO

DO k = ke-1, 1, -1
  DO j = jstart, jend
    DO  i = istart, iend
      zd(i,j,k) = zc3(i,j,k) * zd(i,j,k+1) + zd(i,j,k)
    ENDDO
  ENDDO
ENDDO

!------------------------------------------------------------------------------
!  End of internal procedure precond
!------------------------------------------------------------------------------

END SUBROUTINE precond

!==============================================================================
!+ Internal procedure in "semi_implicit" for matrix-vector product
!------------------------------------------------------------------------------

SUBROUTINE matvec ( zp, zpi, zfx, zfydn, zfyds, zp0u, zp0v, &
                    zrhoqx_i, zrhoqy_i, ie, je, ke )

!------------------------------------------------------------------------------
!
! Description:
!   This internal procedure of "semi_implicit" applies
!   the elliptic operator to the pressure perturbation
!   tendency in the form a matrix-vector product.     
!
! Method:
!   Application of the elliptic operator to the current solution
!   vector (wk1) is required during iteration of the GMRES solver.
!   The GMRES subroutine employs a reverse communication protocol:
!   requiring the application (LM) to invoke matvec and return
!   the matrix-vector product to GMRES in the vector (wk2) . 
!
!------------------------------------------------------------------------------
! Subroutine arguments:
! --------------------
  INTEGER (KIND=iintegers), INTENT (IN   ) ::  &
    ie, je, ke             !  Array dimensions

  REAL    (KIND=wp   ),     INTENT (INOUT) ::  &
    zpi     (ie,je,ke)  ,      & ! Pressure Tendency Field
    zfx     (je   )     ,      & !
    zfydn   (je   )     ,      & !
    zfyds   (je   )     ,      & !
    zp0u    (ie,je,ke)  ,      & !
    zp0v    (ie,je,ke)  ,      & !
    zrhoqx_i(ie,je,ke)  ,      & !
    zrhoqy_i(ie,je,ke)           !

  REAL    (KIND=wp   ),     INTENT (  OUT) ::  &
    zp (ie,je,ke)          !

! Local scalars:
! -----------------------------
  INTEGER (KIND=iintegers) ::  &
    i,  j,  k,           & !  Loop indices in longitudinal, latitudinal and
                           !  vertical direction 
    kp1, km1,            & !  k+1, k-1
    ilowu, ilowv,        & !  start- and end-indices for subdomains with
    jlowu, jlowv,        & !  or without neighbours
    iup,   jup             !  

  REAL    (KIND=wp   )     ::  &
    zbt,                 & !
    zpdx, zpdy,          & !
    zdppz,               & !
    zppqx, zppqy,        & !
    zpgradx, zpgrady,    & !
    zfact,               & !
    zpn,                 & !
    zpfac,               & !
    zjpvn, zjpvs,        & !
    zjlur, zjlul           !

   
! Local (automatic) arrays:
! -----------------------------
  REAL    (KIND=wp   )     ::  &
    zppqz   (ie,je)     ,      & !
    zdu     (ie,je,ke)  ,      & !
    zdv     (ie,je,ke)  ,      & !
    zdw     (ie,je,ke1) ,      & !
    zd      (ie,je,ke1)          !

! For error handling
! ------------------
  INTEGER (KIND=iintegers) ::  izerror

  CHARACTER (LEN=80)       ::  yzerrmsg

!
! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine matvec
!------------------------------------------------------------------------------
 
izerror  = 0
yzerrmsg = '   '

!------------------------------------------------------------------------------
!  Section 1: Setup of subdomain indices to compute pressure gradients
!------------------------------------------------------------------------------

zbt = (1.0_wp + betasw) * dt

! Compute local subdomain start and end indices

IF (my_cart_neigh(1) == -1) THEN ! west
  ilowu = istartu
  ilowv = istartv
ELSE
  ilowu = istartu-1
  ilowv = istartv-1
ENDIF

IF (my_cart_neigh(2) == -1) THEN ! north
  jup = jend
ELSE
  jup = jend+1
ENDIF

IF (my_cart_neigh(3) == -1) THEN ! east
  iup = iend
ELSE
  iup = iend+1
ENDIF

IF (my_cart_neigh(4) == -1) THEN ! south
  jlowu = jstartu
  jlowv = jstartv
ELSE
  jlowu = jstartu-1
  jlowv = jstartv-1
ENDIF

! Problem is that zp must be computed this way
  zp(:,:,:) = zpi(:,:,:)

! Lateral boundary conditions (Neumann) on pressure gradients zdu and zdv

IF (my_cart_neigh(1) == -1) THEN  ! west
  DO k = 1, ke
    DO j = jstartu, jendu
      zdu(istartu-1,j,k) = 0.0_wp
    ENDDO
  ENDDO
END IF

IF (my_cart_neigh(3) == -1) THEN  ! east
  DO k = 1, ke
    DO j = jstartu, jendu
      zdu(iendu+1,j,k) = 0.0_wp
    ENDDO
  ENDDO
END IF

IF (my_cart_neigh(4) == -1) THEN  ! south
  DO k = 1, ke
    DO i = istartv, iendv
      zdv(i,jstartv-1,k) = 0.0_wp
    ENDDO
  ENDDO
END IF

IF (my_cart_neigh(2) == -1) THEN  ! north
  DO k = 1, ke
    DO i = istartv, iendv
      zdv(i,jendv+1,k) = 0.0_wp
    ENDDO
  ENDDO
END IF

! Compute horizontal pressure gradients zdu and zdv

IF (ltime) THEN
  CALL get_timings (i_semi_implicit, ntstep, dt, izerror)
  IF (ltime_barrier) THEN
    CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
    CALL get_timings (i_semi_impl_barrier, ntstep, dt, izerror)
  ENDIF
ENDIF

kzdims(1:24)=(/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
CALL exchg_boundaries                                                  &
   (19,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
    ie, je, kzdims, jstartpar, jendpar, 1, nboundlines,                &
    my_cart_neigh, lperi_x, lperi_y, l2dim,                            &
    10000+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,        &
    zpi(:,:,:) )
IF (ltime) CALL get_timings (i_semi_impl_comm, ntstep, dt, izerror)

DO  k = 1, ke

  kp1   = MIN ( ke, k+1 )
  km1   = MAX ( 1 , k-1 )

  IF ( k < kzflat ) THEN
    DO  j = jlowu, jendu
      DO  i = ilowu, iendu   
        zpgradx =  (zpi(i+1,j,k) - zpi(i,j,k)) * zrhoqx_i(i,j,k)
        zdu(i,j,k) = zpgradx
      ENDDO 
      IF ( j <= jendv ) THEN
        DO  i = ilowv, iendv   
          zpgrady  =  (zpi(i,j+1,k) - zpi(i,j,k)) * zrhoqy_i(i,j,k)
          zdv(i,j,k) = zpgrady
        ENDDO   
      ENDIF   
    ENDDO   
  ELSE
    zppqz(istart-1:iend+1,jstart-1:jend+1) =                  &
                  zpi(istart-1:iend+1,jstart-1:jend+1,kp1)    &
                 -zpi(istart-1:iend+1,jstart-1:jend+1,km1)

    DO  j = jlowu, jendu
      DO  i = ilowu, iendu   
        zppqx   =  zppqz(i+1,j) + zppqz(i,j)
        zdppz   =  zppqx * zp0u(i,j,k)
        zpdx    =  zpi(i+1,j,k) - zpi(i,j,k)
        zpgradx =  zpdx - zdppz
        zdu(i,j,k) = zpgradx*zrhoqx_i(i,j,k)
      ENDDO 
      IF ( j <= jendv ) THEN
        DO  i = ilowv, iendv   
          zppqy   =  zppqz(i,j+1) + zppqz(i,j)
          zdppz   =  zppqy * zp0v(i,j,k)
          zpdy    =  zpi(i,j+1,k) - zpi(i,j,k)
          zpgrady =  zpdy - zdppz
          zdv(i,j,k) = zpgrady*zrhoqy_i(i,j,k)
        ENDDO 
      ENDIF
    ENDDO   
  ENDIF

ENDDO

! Compute horizontal divergence Dh( zdu, zdv )

DO k = 1, ke
  kp1 = MIN( ke, k+1 )
  km1 = MAX( 1 , k-1 )
  zfact = 1.0_wp
  IF(k==1 .OR. k==ke) zfact=0.5_wp
  IF ( k < kzflat ) THEN
    DO j = jstart, jend
    DO i = istart, iend
      zpi(i,j,k) = zfx(j)*( zdu(i,j,k) - zdu(i-1,j,k) ) &
                  + zfydn(j)*zdv(i,j,k) - zfyds(j)*zdv(i,j-1,k)
    ENDDO
    ENDDO
  ELSE
    DO j = jstart, jend
    DO i = istart, iend
      zpi(i,j,k) =                                                         &
       zfx(j)*( zdu(i,j,k) - zdu(i-1,j,k)                                  &
        - zfact * (zdu(i-1,j,kp1) - zdu(i-1,j,km1))* zp0u(i-1,j,k)         &
        - zfact * (zdu(i  ,j,kp1) - zdu(i  ,j,km1))* zp0u(i  ,j,k) )       &
        + zfydn(j)*(  zdv(i,j  ,k)                                         &
        - zfact * (zdv(i,j,kp1)   - zdv(i,j,km1))  * zp0v(i,j,k) )         &
        + zfyds(j)*(- zdv(i,j-1,k)                                         &
        - zfact * (zdv(i,j-1,kp1) - zdv(i,j-1,km1)) * zp0v(i,j-1,k) )
    ENDDO
    ENDDO
  ENDIF
ENDDO

! Top and bottom boundary conditions

DO j = jstart, jend
  DO  i = istart, iend
    zjpvn = ( hhl(i  ,j+1,ke1) - hhl(i  ,j  ,ke1) ) * zdv(i,j  ,ke)
    zjpvs = ( hhl(i  ,j  ,ke1) - hhl(i  ,j-1,ke1) ) * zdv(i,j-1,ke)
    zjlur = ( hhl(i+1,j  ,ke1) - hhl(i  ,j  ,ke1) ) * zdu(i  ,j,ke)
    zjlul = ( hhl(i  ,j  ,ke1) - hhl(i-1,j  ,ke1) ) * zdu(i-1,j,ke)
    zdw(i,j,ke1) = 0.5_wp*(  zfx(j)*( zjlur + zjlul ) +                        &
                          zfydn(j)*zjpvn + zfyds(j)*zjpvs )
    zdw(i,j,1)   = 0.0_wp
    zd(i,j,ke1)  = zdw(i,j,ke1)
    zd(i,j,1)    = zdw(i,j,1)
  ENDDO
ENDDO

! Compute zdw = Dz1.d(pp)

DO k = 2, ke
  DO j = jstart, jend
    DO  i = istart, iend
      zdw(i,j,k) = zb1(i,j,k)*zp(i,j,k) - za1(i,j,k)*zp(i,j,k-1)
    ENDDO
  ENDDO
ENDDO

! Tridiagonal solve L.U.zd = zdw

DO k = 2, ke
  DO j = jstart, jend
    DO  i = istart, iend
      zd(i,j,k) = ( zdw(i,j,k) - za(i,j,k) * zd(i,j,k-1) ) * zb(i,j,k)
    ENDDO
  ENDDO
ENDDO

DO k = ke, 2, -1
  DO j = jstart, jend
    DO  i = istart, iend
      zd(i,j,k) = zc(i,j,k) * zd(i,j,k+1) + zd(i,j,k)
    ENDDO
  ENDDO
ENDDO

! zp = [ 1/gamma.pn - dt^2.div.1/rho.grad - dt^2.Dz2.Ni.Dz1 ].d(pp)

DO k = 1, ke
  DO j = jstart, jend
    DO  i = istart, iend
      zpn        = p0(i,j,k) + pp(i,j,k,nnow)
      zpfac      = 1.0_wp / ( gamma * zpn )
      zp(i,j,k)  = zpfac * zp(i,j,k) - zbt * zbt * zpi(i,j,k)
      zp(i,j,k)  = zp(i,j,k) - zbt * zbt * (                               &
                   zb2(i,j,k)*zd(i,j,k+1) - za2(i,j,k)*zd(i,j,k) )
    ENDDO
  ENDDO
ENDDO

!------------------------------------------------------------------------------
!  End of internal procedure matvec
!------------------------------------------------------------------------------

END SUBROUTINE matvec

!==============================================================================

!------------------------------------------------------------------------------
! End of external procedure semi_implicit
!------------------------------------------------------------------------------

END SUBROUTINE semi_implicit
