!+ Source module for handling the boundary conditions
!------------------------------------------------------------------------------

MODULE src_lbc

!------------------------------------------------------------------------------
!
! Description:
!   This module contains subroutines to handle the boundary conditions
!   in the model. Each time a stencil is computed, one should check
!   if something has to be done at the boundaries before applying the
!   next stencil in order to respect the boundary condition type
!   specified either implicitely (as it is the case currently
!   for T, U, V, W, PP and other variables) or explicitly
!   (as done for the tracers). The subroutines available in this
!   modules come to hand for performing such operations.
!
! Current Code Owner: MeteoSwiss, Oliver Fuhrer
!  phone:  +41 58 460 9359
!  fax:    +41 58 460 9278
!  email:  oliver.fuhrer@meteoswiss.ch
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V5_4b        2016-07-12 Oliver Fuhrer, Carlos Osuna, Anne Roches, Pascal Spoerri
!  Initial release
!  Renamed lbc_mirror to lbc_zerograd. Replaced calls to STOP with model_abort
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!
!==============================================================================
!
! Declarations:
!
! Modules used:

USE kind_parameters,    ONLY :  &
    wp                ! KIND-type parameters for real variables

USE data_modelconfig,   ONLY :  &
    ie         ,    & ! number of grid points in zonal direction
    je         ,    & ! number of grid points in meridional direction
    ke                ! number of grid points in vertical direction

USE data_parallel,      ONLY :  &
   nboundlines,     & ! number of boundary lines of the domain for which
                      ! no forecast is computed = overlapping boundary
    my_cart_id,     & ! rank of this subdomain in the cartesian communicator
    my_cart_neigh     ! neighbors of this subdomain in the cartesian grid

USE data_runcontrol,    ONLY :  &
    nlastbound_orig => nlastbound, & ! time step of the last boundary update
    nincbound_orig  => nincbound,  & ! time step increment of boundary update
    ntstep
 
USE environment,        ONLY :  &
    model_abort       ! process stops the whole parallel program

!==============================================================================

IMPLICIT NONE

!==============================================================================

! default private
PRIVATE

!------------------------------------------------------------------------------
! Declarations
!------------------------------------------------------------------------------

! Parameters
!-----------

! cardinal directions (same order chosen as for my_cart_neigh() indexing)
INTEGER, PARAMETER :: &
  WestBoundary     = 1, &
  NorthBoundary    = 2, &
  EastBoundary     = 3, &
  SouthBoundary    = 4, &
  SouthWestCorner  = 5, &
  SouthEastCorner  = 6, &
  NorthWestCorner  = 7, &
  NorthEastCorner  = 8

! number of boundaries
INTEGER, PARAMETER :: &
  NBoundaries = 4

! array containing all boundaries
INTEGER, PARAMETER :: &
  AllBoundaries(NBoundaries) =         &
    (/ WestBoundary, SouthBoundary, EastBoundary, NorthBoundary /)

! type of field, that will determine behaviour of some algorithms applied for boundary conditions
INTEGER, PARAMETER :: &
  BCFieldType_Scalar  = 1, &
  BCFieldType_VectorI = 2, &
  BCFieldType_VectorJ = 3

! Public methods
!---------------

! exported constants
PUBLIC :: BCFieldType_Scalar, BCFieldType_VectorI,         &
          BCFieldType_VectorJ

! exported methods
PUBLIC :: lbc_value, lbc_copy, lbc_zerograd, lbc_interpolate, &
          lbc_tendency, lbc_compute_tendency

!==============================================================================
! Module Procedures in src_lbc
!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure lbc_value to copy a value for a 3D var
!------------------------------------------------------------------------------

SUBROUTINE lbc_value(field, value, id, jd, kd, ist, ien, jst, jen, kst, ken,       &
  field_type, nlines, doEW, doNS, doCorners, mask, lacc)

!------------------------------------------------------------------------------
! Description:
!  Implementation of the copy value BC. A single value is copied to 
!    all boundary points.
!
! Arguments:
!  field               , field where to apply bc
!  lacc                , flag to run on host CPU or GPU
!  id, jd, kd          , dimensions of field
!  ist, ien            , start/end points of the inner domain
!  jst, jen
!  kst, ken
!  field_type          , type of field (scalar or vector component)
!
! Optional arguments:
!  nlines              , number of lines at the boundary to apply bc
!  doEW                , if TRUE apply bc to WestEast boundaries
!  doNS                , if TRUE apply bc to NorthSouth boundaries
!  doCorners           , if TRUE apply bc to corners
!  mask                , mask that determines for every point whether bc are applied or not
!  lacc                , if TRUE apply bc on GPU
!
!------------------------------------------------------------------------------

  ! Subroutine arguments:
  ! --------------------
  INTEGER,       INTENT(IN)           :: id, jd, kd
  REAL(KIND=wp), INTENT(INOUT)        :: field(id,jd,kd)
  REAL(KIND=wp), INTENT(IN)           :: value
  INTEGER,       INTENT(IN)           :: ist, ien, jst, jen, kst, ken
  INTEGER,       INTENT(IN), OPTIONAL :: field_type
  INTEGER,       INTENT(IN), OPTIONAL :: nlines
  LOGICAL,       INTENT(IN), OPTIONAL :: doEW, doNS, doCorners
  LOGICAL,                  INTENT(IN), OPTIONAL :: mask(id,jd,kd)
  LOGICAL,                  INTENT(IN), OPTIONAL :: lacc

  ! Local variables
  INTEGER  :: i, j, k, idir
  INTEGER  :: ibc_start, ibc_end, jbc_start, jbc_end
  LOGICAL  :: wbc, ebc, nbc, sbc
  INTEGER  :: znlines, zfield_type
  LOGICAL  :: zdoEW, zdoNS, zdoCorners
  LOGICAL  :: lzacc
  
!------------------------------------------------------------------------------
! Begin module procedure lbc_value
!------------------------------------------------------------------------------

  CALL check_setup(zfield_type, znlines, zdoEW, zdoNS, zdoCorners, lzacc, &
    wbc, ebc, nbc, sbc, field_type, nlines, doEW, doNS, doCorners, lacc )

  !NOacc data present(field, mask) if (lzacc)

  ! loop over sides of subdomain
  DO idir = 1, NBoundaries
    IF ( apply_bc( AllBoundaries(idir), wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners ) ) THEN

      ! get index range of boundary grid points
      CALL get_bc_ranges( AllBoundaries(idir), znlines, ist, ien, jst, jen,    &
                          wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners,        &
                          ibc_start, ibc_end, jbc_start, jbc_end )

      ! apply value to the boundary grid points
      IF ( PRESENT(mask) ) THEN
        !NOacc kernels if (lzacc)
        DO k = kst, ken
          DO j = jbc_start, jbc_end
            DO i = ibc_start, ibc_end
              IF ( mask(i,j,k) ) THEN
                field(i,j,k) = value
              END IF
            END DO
          END DO
        END DO
        !NOacc end kernels
      ELSE
        !NOacc kernels if (lzacc) 
        DO k = kst, ken
          DO j = jbc_start, jbc_end
            DO i = ibc_start, ibc_end
              field(i,j,k) = value
            END DO
          END DO
        END DO
        !NOacc end kernels
      END IF

    END IF
  END DO
  !NOacc end data

!------------------------------------------------------------------------------
! End of module procedure lbc_value
!------------------------------------------------------------------------------

END SUBROUTINE lbc_value


!==============================================================================
!+ Module procedure lbc_copy to copy a field for a 3D var
!------------------------------------------------------------------------------

SUBROUTINE lbc_copy( field, src, id, jd, kd, ist, ien, jst, jen, kst, ken,  &
  field_type, nlines, doEW, doNS, doCorners, mask, lacc )

!------------------------------------------------------------------------------
! Description:
!  Implementation of the copy field BC. A source field is copied to
!    all boundary points.
!
! Arguments:
!  field               , field where to apply bc
!  src                 , source field for boundaries values
!  id, jd, kd          , dimensions of field
!  ist, ien            , start/end points of the inner domain
!  jst, jen
!  kst, ken
!
! Optional arguments:
!  field_type          , type of field (scalar or vector component)
!  nlines              , number of lines at the boundary to apply bc
!  doEW                , if TRUE apply bc to WestEast boundaries
!  doNS                , if TRUE apply bc to NorthSouth boundaries
!  doCorners           , if TRUE apply bc to corners
!  mask                , mask that determines for every point whether bc are applied or not
!  lacc                , flag to run on host CPU or GPU
!
!------------------------------------------------------------------------------

  ! Subroutine arguments:
  ! --------------------
  INTEGER,       INTENT(IN)           :: id, jd, kd
  REAL(KIND=wp), INTENT(INOUT)        :: field(id,jd,kd)
  REAL(KIND=wp), INTENT(IN)           :: src(id,jd,kd)
  INTEGER,       INTENT(IN)           :: ist, ien, jst, jen, kst, ken
  INTEGER,       INTENT(IN), OPTIONAL :: field_type
  INTEGER,       INTENT(IN), OPTIONAL :: nlines
  LOGICAL,       INTENT(IN), OPTIONAL :: doEW, doNS, doCorners
  LOGICAL,       INTENT(IN), OPTIONAL :: mask(id,jd,kd)
  LOGICAL,       INTENT(IN), OPTIONAL :: lacc

  ! Local variables
  INTEGER  :: i, j, k, idir
  INTEGER  :: ibc_start, ibc_end, jbc_start, jbc_end
  LOGICAL  :: wbc, ebc, nbc, sbc
  INTEGER  :: znlines, zfield_type
  LOGICAL  :: zdoEW, zdoNS, zdoCorners
  LOGICAL  :: lzacc

!------------------------------------------------------------------------------
! Begin module procedure lbc_copy
!------------------------------------------------------------------------------

  CALL check_setup(zfield_type, znlines, zdoEW, zdoNS, zdoCorners, lzacc, &
    wbc, ebc, nbc, sbc, field_type, nlines, doEW, doNS, doCorners, lacc )

  !NOacc data present(field, src, mask) if (lzacc)

  ! loop over sides of subdomain
  DO idir = 1, NBoundaries
    IF ( apply_bc( AllBoundaries(idir), wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners ) ) THEN

      ! get index range of boundary grid points
      CALL get_bc_ranges( AllBoundaries(idir), znlines, ist, ien, jst, jen,   &
             wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners,                    &
             ibc_start, ibc_end, jbc_start, jbc_end )

      ! apply field to the boundary grid points
      IF ( PRESENT(mask) ) THEN
        !NOacc kernels if (lzacc)
        DO k = kst, ken
          DO j = jbc_start, jbc_end
            DO i = ibc_start, ibc_end
              IF ( mask(i,j,k) ) THEN
                field(i,j,k) = src(i,j,k)
              END IF
            END DO
          END DO
        END DO
        !NOacc end kernels
      ELSE
        !NOacc kernels if (lzacc)
        DO k = kst, ken
          DO j = jbc_start, jbc_end
            DO i = ibc_start, ibc_end
              field(i,j,k) = src(i,j,k)
            END DO
          END DO
        END DO
        !NOacc end kernels
      END IF

    END IF
  END DO
  !NOacc end data

!------------------------------------------------------------------------------
! End of module procedure lbc_value
!------------------------------------------------------------------------------

END SUBROUTINE lbc_copy


!==============================================================================
!+ Module procedure lbc_zerograd  to apply a 0gradient BC for a 3D var
!------------------------------------------------------------------------------

SUBROUTINE lbc_zerograd( field, id, jd, kd, ist, ien, jst, jen, kst, ken, &
  field_type, nlines, doEW, doNS, doCorners, mask, lacc )

!------------------------------------------------------------------------------
! Description:
!  Implementation of the zero-gradient BC. Points in the boundary zone
!    are copied from corresponding points in the inner domain in order
!    to impose a zero-gradient boundary condition.
!
! Arguments:
!  field               , field where to apply bc
!  id, jd, kd          , dimensions of field
!  ist, ien            , start/end points of the inner domain
!  jst, jen
!  kst, ken
!
! Optional arguments:
!  field_type          , type of field (scalar or vector component)
!  nlines              , number of lines at the boundary to apply bc
!  doEW                , if TRUE apply bc to WestEast boundaries
!  doNS                , if TRUE apply bc to NorthSouth boundaries
!  doCorners           , if TRUE apply bc to corners
!  mask                , mask that determines for every point whether bc are applied or not
!  lacc                , flag to run on host CPU or GPU
!
!------------------------------------------------------------------------------

  ! Subroutine arguments:
  ! --------------------
  INTEGER,       INTENT(IN)           :: id, jd, kd
  REAL(KIND=wp), INTENT(INOUT)        :: field(id,jd,kd)
  INTEGER,       INTENT(IN)           :: ist, ien, jst, jen, kst, ken
  INTEGER,       INTENT(IN), OPTIONAL :: field_type
  INTEGER,       INTENT(IN), OPTIONAL :: nlines
  LOGICAL,       INTENT(IN), OPTIONAL :: doEW, doNS, doCorners
  LOGICAL,       INTENT(IN), OPTIONAL :: mask(id, jd, kd)
  LOGICAL,       INTENT(IN), OPTIONAL :: lacc

  ! Local variables
  INTEGER  :: i, j, k
  INTEGER  :: ibc_start, ibc_end, jbc_start, jbc_end
  LOGICAL  :: wbc, ebc, nbc, sbc
  INTEGER  :: znlines, zfield_type
  LOGICAL  :: zdoEW, zdoNS, zdoCorners
  LOGICAL  :: lzacc

!------------------------------------------------------------------------------
! Begin module procedure lbc_zerograd
!------------------------------------------------------------------------------

  CALL check_setup(zfield_type, znlines, zdoEW, zdoNS, zdoCorners, lzacc, &
    wbc, ebc, nbc, sbc, field_type, nlines, doEW, doNS, doCorners, lacc )

  !NOacc data present(field, mask) if (lzacc)

  ! West boundary
  ! ---------------
  IF ( apply_bc( WestBoundary, wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners ) ) THEN

    ! get index range of boundary grid points
    CALL get_bc_ranges( WestBoundary, znlines, ist, ien, jst, jen,   &
           wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners,             &
           ibc_start, ibc_end, jbc_start, jbc_end )

    ! apply zero gradient
    IF ( PRESENT(mask) ) THEN
      !NOacc kernels if (lzacc)
      DO k = kst, ken
        DO j = jbc_start, jbc_end
          DO i = ibc_start, ibc_end
            IF ( mask(i,j,k) ) THEN
              field(i,j,k) = field(ist,j,k)
            END IF
          END DO
        END DO
      END DO
      !NOacc end kernels
    ELSE
      !NOacc kernels if (lzacc)
      DO k = kst, ken
        DO j = jbc_start, jbc_end
          DO i = ibc_start, ibc_end
            field(i,j,k) = field(ist,j,k)
          END DO
        END DO
      END DO
      !NOacc end kernels
    END IF
    
  END IF


  ! East boundary
  ! ---------------
  IF ( apply_bc( EastBoundary, wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners ) ) THEN

    ! get index range of boundary grid points
    CALL get_bc_ranges( EastBoundary, znlines, ist, ien, jst, jen,   &
           wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners,             &
           ibc_start, ibc_end, jbc_start, jbc_end )

    ! apply zero gradient
    IF ( PRESENT(mask) ) THEN
      !NOacc kernels if (lzacc)
      DO k = kst, ken
        DO j = jbc_start, jbc_end
          DO i = ibc_start, ibc_end
            IF ( mask(i,j,k) ) THEN
              field(i,j,k) = field(ien,j,k)
            END IF
          END DO
        END DO
      END DO
      !NOacc end kernels
    ELSE
      !NOacc kernels if (lzacc)
      DO k = kst, ken
        DO j = jbc_start, jbc_end
          DO i = ibc_start, ibc_end
            field(i,j,k) = field(ien,j,k)
          END DO
        END DO
      END DO
      !NOacc end kernels
    END IF

  END IF


  ! North boundary
  ! ---------------
  IF ( apply_bc( NorthBoundary, wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners ) ) THEN

    ! get index range of boundary grid points
    CALL get_bc_ranges( NorthBoundary, znlines, ist, ien, jst, jen,   &
           wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners,              &
           ibc_start, ibc_end, jbc_start, jbc_end )

    ! apply zero gradient
    IF ( PRESENT(mask) ) THEN
      !NOacc kernels  if (lzacc)
      DO k = kst, ken
        DO j = jbc_start, jbc_end
          DO i = ibc_start, ibc_end
            IF ( mask(i,j,k) ) THEN
              field(i,j,k) = field(i,jen,k)
            END IF
          END DO
        END DO
      END DO
      !NOacc end kernels
    ELSE
       !NOacc kernels  if (lzacc)
      DO k = kst, ken
        DO j = jbc_start, jbc_end
          DO i = ibc_start, ibc_end
            field(i,j,k) = field(i,jen,k)
          END DO
        END DO
      END DO
      !NOacc end kernels
    END IF

  END IF


  ! South boundary
  ! ---------------
  IF ( apply_bc( SouthBoundary, wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners ) ) THEN

    ! get index range of boundary grid points
    CALL get_bc_ranges( SouthBoundary, znlines, ist, ien, jst, jen,   &
           wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners,              &
           ibc_start, ibc_end, jbc_start, jbc_end )

    ! apply zero gradient
    IF ( PRESENT(mask) ) THEN
      !NOacc kernels if (lzacc)
      DO k = kst, ken
        DO j = jbc_start, jbc_end
          DO i = ibc_start, ibc_end
            IF ( mask(i,j,k) ) THEN
              field(i,j,k) = field(i,jst,k)
            END IF
          END DO
        END DO
      END DO
      !NOacc end kernels
    ELSE
      !NOacc kernels if (lzacc)
      DO k = kst, ken
        DO j = jbc_start, jbc_end
          DO i = ibc_start, ibc_end
            field(i,j,k) = field(i,jst,k)
          END DO
        END DO
      END DO
      !NOacc end kernels
    END IF

  END IF


  ! South-West corner
  ! -----------------
  IF ( apply_bc( SouthWestCorner, wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners ) ) THEN

    ! get index range of boundary grid points
    CALL get_bc_ranges( SouthWestCorner, znlines, ist, ien, jst, jen,   &
           wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners,                &
           ibc_start, ibc_end, jbc_start, jbc_end )

    ! apply zero gradient
    IF ( PRESENT(mask) ) THEN
      !NOacc kernels if (lzacc)
      DO k = kst, ken
        DO j = jbc_start, jbc_end
          DO i = ibc_start, ibc_end
            IF ( mask(i,j,k) ) THEN
              field(i,j,k) = field(ist,jst,k)
            END IF
          END DO
        END DO
      END DO
      !NOacc end kernels
    ELSE
      !NOacc kernels if (lzacc)
      DO k = kst, ken
        DO j = jbc_start, jbc_end
          DO i = ibc_start, ibc_end
            field(i,j,k) = field(ist,jst,k)
          END DO
        END DO
      END DO
      !NOacc end kernels
    END IF

  END IF


  ! South-East corner
  ! -----------------
  IF ( apply_bc( SouthEastCorner, wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners ) ) THEN

    ! get index range of boundary grid points
    CALL get_bc_ranges( SouthEastCorner, znlines, ist, ien, jst, jen,   &
           wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners,                &
           ibc_start, ibc_end, jbc_start, jbc_end )

    ! apply zero gradient
    IF ( PRESENT(mask) ) THEN
     !NOacc kernels if (lzacc)
      DO k = kst, ken
        DO j = jbc_start, jbc_end
          DO i = ibc_start, ibc_end
            IF ( mask(i,j,k) ) THEN
              field(i,j,k) = field(ien,jst,k)
            END IF
          END DO
        END DO
      END DO
      !NOacc end kernels
    ELSE
      !NOacc kernels if (lzacc)
      DO k = kst, ken
        DO j = jbc_start, jbc_end
          DO i = ibc_start, ibc_end
            field(i,j,k) = field(ien,jst,k)
          END DO
        END DO
      END DO
      !NOacc end kernels
    END IF

  END IF


  ! North-West corner
  ! -----------------
  IF ( apply_bc( NorthWestCorner, wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners ) ) THEN

    ! get index range of boundary grid points
    CALL get_bc_ranges( NorthWestCorner, znlines, ist, ien, jst, jen,   &
           wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners,                &
           ibc_start, ibc_end, jbc_start, jbc_end )

    ! apply zero gradient
    IF ( PRESENT(mask) ) THEN
      !NOacc kernels if (lzacc)
      DO k = kst, ken
        DO j = jbc_start, jbc_end
          DO i = ibc_start, ibc_end
            IF ( mask(i,j,k) ) THEN
              field(i,j,k) = field(ist,jen,k)
            END IF
          END DO
        END DO
      END DO
      !NOacc end kernels
    ELSE
     !NOacc kernels if (lzacc)
      DO k = kst, ken
        DO j = jbc_start, jbc_end
          DO i = ibc_start, ibc_end
            field(i,j,k) = field(ist,jen,k)
          END DO
        END DO
      END DO
      !NOacc end kernels
    END IF

  END IF


  ! North-East corner
  ! -----------------
  IF ( apply_bc( NorthEastCorner, wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners ) ) THEN

    ! get index range of boundary grid points
    CALL get_bc_ranges( NorthEastCorner, znlines, ist, ien, jst, jen,    &
           wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners,                 &
           ibc_start, ibc_end, jbc_start, jbc_end )

    ! apply zero gradient
    IF ( PRESENT(mask) ) THEN
      !NOacc kernels if (lzacc)
      DO k = kst, ken
        DO j = jbc_start, jbc_end
          DO i = ibc_start, ibc_end
            IF ( mask(i,j,k) ) THEN
              field(i,j,k) = field(ien,jen,k)
            END IF
          END DO
        END DO
      END DO
      !NOacc end kernels
    ELSE
      !NOacc kernels if (lzacc)
      DO k = kst, ken
        DO j = jbc_start, jbc_end
          DO i = ibc_start, ibc_end
            field(i,j,k) = field(ien,jen,k)
          END DO
        END DO
      END DO
      !NOacc end kernels
    END IF

  END IF

  !NOacc end data

!------------------------------------------------------------------------------
! End of module procedure lbc_zerograd
!------------------------------------------------------------------------------

END SUBROUTINE lbc_zerograd


!==============================================================================
!+ Module procedure lbc_interpolate to interpolate bd fields for 3D var 
!------------------------------------------------------------------------------

SUBROUTINE lbc_interpolate( field, bd1, bd2, id, jd, kd,     &
  ist, ien, jst, jen, kst, ken,                              &
  field_type, nlines, doEW, doNS, doCorners, nlastbound,     &
  nincbound, mask, lacc)

!------------------------------------------------------------------------------
! Description:
!  Implementation of the "interpolate" BC. Points in the boundary zone
!    are computed as a linear combination of bd1 and bd2.
!
! Arguments:
!  field               , field where to apply bc
!  bd1                 , first boundary field used to interpolate boundary values
!  bd2                 , second boundary field used to interpolate boundary values
!  id, jd, kd          , dimensions of field
!  ist, ien            , start/end points of the inner domain
!  jst, jen
!  kst, ken
!
! Optional Arguments:
!  field_type          , type of field (scalar or vector component)
!  nlines              , number of lines at the boundary to apply bc
!  doEW                , if TRUE apply bc to WestEast boundaries
!  doNS                , if TRUE apply bc to NorthSouth boundaries
!  doCorners           , if TRUE apply bc to corners
!  nlastbound          , # last time step when boundary fields were updated
!  nincbound           , # time steps between two boundary fields
!  mask                , mask that determines for every point whether bc are applied or not
!  lacc                , flags to run host CPU or GPU
!
!------------------------------------------------------------------------------

  ! Subroutine arguments:
  ! --------------------
  INTEGER,       INTENT(IN)           :: id, jd, kd
  REAL(KIND=wp), INTENT(INOUT)        :: field(id,jd,kd)
  REAL(KIND=wp), INTENT(IN)           :: bd1(id,jd,kd), bd2(id,jd,kd)
  INTEGER,       INTENT(IN)           :: ist, ien, jst, jen, kst, ken
  INTEGER,       INTENT(IN), OPTIONAL :: field_type
  INTEGER,       INTENT(IN), OPTIONAL :: nlines
  LOGICAL,       INTENT(IN), OPTIONAL :: doEW, doNS, doCorners
  INTEGER,       INTENT(IN), OPTIONAL :: nlastbound
  INTEGER,       INTENT(IN), OPTIONAL :: nincbound
  LOGICAL,       INTENT(IN), OPTIONAL :: mask(id,jd,kd)
  LOGICAL,       INTENT(IN), OPTIONAL :: lacc

  ! Local variables
  INTEGER       :: i, j, k, idir
  REAL(KIND=wp) :: z1, z2
  INTEGER       :: ibc_start, ibc_end, jbc_start, jbc_end
  LOGICAL       :: wbc, ebc, nbc, sbc
  INTEGER       :: znlines, zfield_type
  LOGICAL       :: zdoEW, zdoNS, zdoCorners
  INTEGER       :: znlastbound, znincbound
  LOGICAL       :: lzacc

!------------------------------------------------------------------------------
! Begin module procedure lbc_interpolate
!------------------------------------------------------------------------------

  CALL check_setup(zfield_type, znlines, zdoEW, zdoNS, zdoCorners, lzacc, &
    wbc, ebc, nbc, sbc, field_type, nlines, doEW, doNS, doCorners, lacc )

  IF ( PRESENT(nlastbound) ) THEN
    znlastbound = nlastbound
  ELSE
    znlastbound = nlastbound_orig
  END IF
  IF ( PRESENT(nincbound) ) THEN
    znincbound = nincbound
  ELSE
    znincbound = nincbound_orig
  END IF

  z2 = REAL( ntstep + 1 - znlastbound, wp) / REAL(znincbound, wp)
  z2 = MIN( 1.0_wp, z2)
  z1 = 1.0_wp - z2

  !NOacc data present(field, bd1, bd2, mask) if (lzacc)

  ! loop over sides of subdomain
  DO idir = 1, NBoundaries
    IF ( apply_bc( AllBoundaries(idir), wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners ) ) THEN

      ! get index range of boundary grid points
      CALL get_bc_ranges( AllBoundaries(idir), znlines, ist, ien, jst, jen,  &
             wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners,                   &
             ibc_start, ibc_end, jbc_start, jbc_end )

      ! apply linear combination
      IF ( PRESENT(mask) ) THEN
        !NOacc kernels if (lzacc)
        DO k = kst, ken
          DO j = jbc_start, jbc_end
            DO i = ibc_start, ibc_end
              IF ( mask(i,j,k) ) THEN
                field(i,j,k) = z1 * bd1(i,j,k) + z2 * bd2(i,j,k)
              END IF
            END DO
          END DO
        END DO
        !NOacc end kernels
      ELSE
        !NOacc kernels if (lzacc)
        DO k = kst, ken
          DO j = jbc_start, jbc_end
            DO i = ibc_start, ibc_end
              field(i,j,k) = z1 * bd1(i,j,k) + z2 * bd2(i,j,k)
            END DO
          END DO
        END DO
        !NOacc end kernels
      END IF

    END IF
  END DO
  !NOacc end data

!------------------------------------------------------------------------------
! End of module procedure lbc_interpolate
!------------------------------------------------------------------------------

END SUBROUTINE lbc_interpolate


!==============================================================================
!+ Module procedure lbc_tendency to derive BC from tendency for 3D var 
!------------------------------------------------------------------------------

SUBROUTINE lbc_tendency( field, tend, dt, id, jd, kd, ist, ien, jst, jen, kst, ken,  &
  field_type, nlines, doEW, doNS, doCorners, mask, lacc )

!------------------------------------------------------------------------------
! Description:
!  Compute tendency from the difference of two fields
!
! Arguments:
!  field               , field where to apply bc
!  tend                , tendency for boundaries values
!  dt                  , time increment (for tendency)
!  id, jd, kd          , dimensions of field
!  ist, ien            , start/end points of the inner domain
!  jst, jen
!  kst, ken
!
! Optional Arguments:
!  field_type          , type of field (scalar or vector component)
!  nlines              , number of lines at the boundary to apply bc
!  doEW                , if TRUE apply bc to WestEast boundaries
!  doNS                , if TRUE apply bc to NorthSouth boundaries
!  doCorners           , if TRUE apply bc to corners
!  mask                , mask that determines for every point whether bc are applied or not
!  lacc                , flag to run on host CPU or GPU
!
!------------------------------------------------------------------------------

  ! Subroutine arguments:
  ! --------------------
  INTEGER,       INTENT(IN)           :: id, jd, kd
  REAL(KIND=wp), INTENT(INOUT)        :: field(id,jd,kd)
  REAL(KIND=wp), INTENT(IN)           :: tend(id,jd,kd)
  REAL(KIND=wp), INTENT(IN)           :: dt
  INTEGER,       INTENT(IN)           :: ist, ien, jst, jen, kst, ken
  INTEGER,       INTENT(IN), OPTIONAL :: field_type
  INTEGER,       INTENT(IN), OPTIONAL :: nlines
  LOGICAL,       INTENT(IN), OPTIONAL :: doEW, doNS, doCorners
  LOGICAL,       INTENT(IN), OPTIONAL :: lacc
  LOGICAL,       INTENT(IN), OPTIONAL :: mask(id,jd,kd)

  ! Local variables
  INTEGER  :: i, j, k, idir
  INTEGER  :: ibc_start, ibc_end, jbc_start, jbc_end
  LOGICAL  :: wbc, ebc, nbc, sbc
  INTEGER  :: znlines, zfield_type
  LOGICAL  :: zdoEW, zdoNS, zdoCorners
  LOGICAL  :: lzacc

!------------------------------------------------------------------------------
! Begin module procedure lbc_tendency
!------------------------------------------------------------------------------

  CALL check_setup( zfield_type, znlines, zdoEW, zdoNS, zdoCorners, lzacc, &
    wbc, ebc, nbc, sbc, field_type, nlines, doEW, doNS, doCorners, lacc )

  !NOacc data present(field, tend, mask) if (lzacc)

  ! loop over sides of subdomain
  DO idir = 1, NBoundaries
    IF ( apply_bc( AllBoundaries(idir), wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners ) ) THEN

      ! get index range of boundary grid points
      CALL get_bc_ranges( AllBoundaries(idir), znlines, ist, ien, jst, jen,  &
             wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners,                   &
             ibc_start, ibc_end, jbc_start, jbc_end )

      ! derive BC using the tendency field
      IF ( PRESENT(mask) ) THEN
        !NOacc kernels if (lzacc)
        DO k = kst, ken
          DO j = jbc_start, jbc_end
            DO i = ibc_start, ibc_end
              IF ( mask(i,j,k) ) THEN
                field(i,j,k) = field(i,j,k) + tend(i,j,k) * dt
              END IF
            END DO
          END DO
        END DO
        !NOacc end kernels
      ELSE
        !NOacc kernels if (lzacc)
        DO k = kst, ken
          DO j = jbc_start, jbc_end
            DO i = ibc_start, ibc_end
              field(i,j,k) = field(i,j,k) + tend(i,j,k) * dt
            END DO
          END DO
        END DO
        !NOacc end kernels
      END IF

    END IF
  END DO

  !NOacc end data

!------------------------------------------------------------------------------
! End of module procedure lbc_tendency
!------------------------------------------------------------------------------

END SUBROUTINE lbc_tendency


!==============================================================================
!+ Module procedure lbc_compute_tendency to compute tendency for tendency BC
!------------------------------------------------------------------------------

SUBROUTINE lbc_compute_tendency( tend, bd1, bd2, dt, id, jd, kd,   &
  ist, ien, jst, jen, kst, ken, field_type, nlines, doEW, doNS, doCorners, lacc )

!------------------------------------------------------------------------------
! Description:
!  Compute a 3d boundary tendency field from two boundary fields.
!
! Arguments:
!  tend                , tendency of boundary field
!  bd1                 , first boundary field used to compute tendency
!  bd2                 , second boundary field used to compute tendency
!  dt                  , time step between two boundary fields
!  id, jd, kd          , dimensions of field
!  ist, ien            , end points of the inner domain
!  jst, jen
!  kst, ken
!
! Optional arguments:
!  field_type          , type of field (scalar or vector component)
!  nlines              , number of lines at the boundary to apply bc
!  doEW                , if TRUE apply bc to WestEast boundaries
!  doNS                , if TRUE apply bc to NorthSouth boundaries
!  doCorners           , if TRUE apply bc to corners
!  lacc                , flag to run on host CPU or GPU
!
!------------------------------------------------------------------------------

  ! arguments
  INTEGER,       INTENT(IN)           :: id, jd, kd
  REAL(KIND=wp), INTENT(INOUT)        :: tend(id,jd,kd)
  REAL(KIND=wp), INTENT(IN)           :: bd1(id,jd,kd), bd2(id,jd,kd)
  REAL(KIND=wp), INTENT(IN)           :: dt
  INTEGER,       INTENT(IN)           :: ist, ien, jst, jen, kst, ken
  INTEGER,       INTENT(IN), OPTIONAL :: field_type
  INTEGER,       INTENT(IN), OPTIONAL :: nlines
  LOGICAL,       INTENT(IN), OPTIONAL :: doEW, doNS, doCorners
  LOGICAL,       INTENT(IN), OPTIONAL :: lacc

  ! local variables
  INTEGER  :: znlines, zfield_type
  LOGICAL  :: wbc, ebc, nbc, sbc
  LOGICAL  :: zdoEW, zdoNS, zdoCorners
  INTEGER  :: ibc_start, ibc_end, jbc_start, jbc_end
  INTEGER  :: i, j, k, idir
  LOGICAL  :: lzacc

!------------------------------------------------------------------------------
! Begin module procedure lbc_compute_tendency
!------------------------------------------------------------------------------

  CALL check_setup(zfield_type, znlines, zdoEW, zdoNS, zdoCorners, lzacc, &
    wbc, ebc, nbc, sbc, field_type, nlines, doEW, doNS, doCorners, lacc )

  !NOacc data present(tend, bd1, bd2) if (lzacc)

  ! loop over sides of subdomain
  DO idir = 1, NBoundaries
    IF ( apply_bc( AllBoundaries(idir), wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners ) ) THEN

      CALL get_bc_ranges( AllBoundaries(idir), znlines, ist, ien, jst, jen,  &
             wbc, ebc, nbc, sbc, zdoEW, zdoNS, zdoCorners,                   &
             ibc_start, ibc_end, jbc_start, jbc_end )

      !NOacc kernels if (lzacc)
      DO k = kst, ken
        DO j = jbc_start, jbc_end
          DO i = ibc_start, ibc_end
            tend(i,j,k) =  (bd1(i,j,k) - bd2(i,j,k)) / dt
          END DO
        END DO
      END DO
      !NOacc end kernels

     END IF
  END DO
  !NOacc end data

!------------------------------------------------------------------------------
! End of module procedure lbc_compute_tendency
!------------------------------------------------------------------------------

END SUBROUTINE lbc_compute_tendency


!==============================================================================
!+ Module procedure check_setup to setup some default values
!------------------------------------------------------------------------------

SUBROUTINE check_setup(zfield_type, znlines, zdoEW, zdoNS, zdoCorners, lzacc, &
  wbc, ebc, nbc, sbc, field_type, nlines, doEW, doNS, doCorners, lacc )

!------------------------------------------------------------------------------
! Description:
!  General backend for applying boundary conditions both at u-, v-, and mass-point 
!  for a 3d field.
! 
! Arguments (output): 
!  zfield_type , type of field (scalar or vector component)
!  znlines     , number of lines at the boundary to be exchanged
!  zdoEW       , if TRUE apply bc to WestEast boundaries (default TRUE)
!  zdoNS       , if TRUE apply bc to NorthSouth boundaries (default TRUE)
!  zdoCorners  , if TRUE apply bc to corners (default TRUE if doEW=doNS=TRUE, FALSE if only one is TRUE)
!  lzacc       , flag to run on host CPU or GPU
!
! Optional Arguments:
!  field_type  , type of field (scalar or vector component)
!  nlines      , number of lines at the boundary to be exchanged
!  doEW        , if TRUE apply bc to WestEast boundaries (default TRUE)
!  doNS        , if TRUE apply bc to NorthSouth boundaries (default TRUE)
!  doCorners   , if TRUE apply bc to corners (default TRUE if doEW=doNS=TRUE, FALSE if only one is TRUE)
!  lacc        , flag to run on host CPU or GPU
!
!------------------------------------------------------------------------------

  ! Subroutine arguments:
  ! --------------------
  INTEGER, INTENT(OUT)           :: zfield_type
  INTEGER, INTENT(OUT)           :: znlines
  LOGICAL, INTENT(OUT)           :: zdoEW, zdoNS, zdoCorners
  LOGICAL, INTENT(OUT)           :: lzacc
  LOGICAL, INTENT(OUT)           :: wbc, ebc, nbc, sbc
  INTEGER, INTENT(IN),  OPTIONAL :: field_type
  INTEGER, INTENT(IN),  OPTIONAL :: nlines
  LOGICAL, INTENT(IN),  OPTIONAL :: doEW, doNS, doCorners
  LOGICAL, INTENT(IN),  OPTIONAL :: lacc

!------------------------------------------------------------------------------
! Begin module procedure check_setup
!------------------------------------------------------------------------------

  ! determine field type (default is scalar)
  IF ( PRESENT(field_type) ) THEN
    zfield_type = field_type
  ELSE
    zfield_type = BCFieldType_Scalar
  END IF
  IF (.NOT. ANY( zfield_type == (/  BCFieldType_Scalar, BCFieldType_VectorI,  &
                                    BCFieldType_VectorJ /) ) ) THEN
    CALL model_abort( my_cart_id, 42010, 'LBC: Unknown BCType encountered','check_setup' )
  END IF

  ! determine number of lines to apply BC
  ! (default is nboundlines + 1 since staggered fields have one boundary line more)
  IF ( PRESENT(nlines) ) THEN
    znlines = nlines
  ELSE
    znlines = nboundlines + 1
  END IF
  ! NOTE: updating nboundlines+1 lines is possible, since staggered fields
  !       at the u- or v-point may have 4 boundary lines
  IF ( znlines < 0 .OR. znlines > nboundlines + 1 ) THEN
    CALL model_abort( my_cart_id, 42011, 'LBC: Invalid number of boundary lines (nlines)','check_setup' )
  END IF

  ! determine flags where to apply BC
  IF ( PRESENT(doEW) ) THEN
    zdoEW = doEW
  ELSE
    zdoEW = .TRUE.
  END IF
  IF ( PRESENT(doNS) ) THEN
    zdoNS = doNS
  ELSE
    zdoNS = .TRUE.
  END IF
  IF ( PRESENT(doCorners) ) THEN
    zdoCorners = doCorners
  ELSE
    ! NOTE: default value for doCorners is TRUE if we are applying BC to all sides, and
    ! false in case the user has request on EW or NS via the doEW/doNS switches
    zdoCorners = .FALSE.
    IF (zdoEW .AND. zdoNS) THEN
      zdoCorners = .TRUE.
    END IF
  END IF
  IF ( zdoCorners .AND. (.NOT. zdoEW .AND. .NOT. zdoNS) ) THEN
    CALL model_abort( my_cart_id, 42012, 'LBC: Can only do corners if either doEW or doNS is .TRUE.', 'check_setup' )
  END IF
  IF ( .NOT. zdoCorners .AND. (zdoEW .AND. zdoNS) ) THEN
    CALL model_abort( my_cart_id, 42013, 'LBC: Can not switch off corners if both doEW and doNS are .TRUE.','check_setup' )
  END IF

  ! set acc flag (to run on gpu or not)
  IF ( PRESENT(lacc) ) THEN
     lzacc=lacc
  ELSE
     lzacc=.FALSE.
  END IF

  ! get position of this PE
  wbc = PE_at_western_boundary() 
  ebc = PE_at_eastern_boundary() 
  nbc = PE_at_northern_boundary()
  sbc = PE_at_southern_boundary()

!------------------------------------------------------------------------------
! End of module procedure check_setup
!------------------------------------------------------------------------------

END SUBROUTINE check_setup


!==============================================================================
!+ Module function apply_bc to check whether to apply a BC
!------------------------------------------------------------------------------

LOGICAL FUNCTION apply_bc( bc_dir, wbc, ebc, nbc, sbc, doEW, doNS, doCorners )

!------------------------------------------------------------------------------
! Description:
!  Function to determine whether we need to apply a BC or not.
!
! Arguments:
!  bc_dir              , indicates which side of the subdomain of this PE we are considering
!  wbc, ebc, nbc, sbc  , flags to enable/disable bc at each boundary (W,E,N,S)
!  doEW                , if TRUE apply bc to WestEast boundaries
!  doNS                , if TRUE apply bc to NorthSouth boundaries
!  doCorners           , if TRUE apply bc to corners
!
! Optional Arguments:
!  mask                , mask that determines for every point whether bc are applied or not
!
!------------------------------------------------------------------------------

  ! arguments
  INTEGER, INTENT(IN) :: bc_dir
  LOGICAL, INTENT(IN) :: wbc, ebc, nbc, sbc
  LOGICAL, INTENT(IN) :: doEW, doNS, doCorners

  ! default to false (don't apply BC)
  apply_bc = .FALSE.

  ! check which sub-domain border we are being queried for
  ! and return FALSE if we are not sitting at the same global domain border
  ! or the user has specifically asked not to apply BCs (via doEW or doNS)
  SELECT CASE(bc_dir)
    CASE(WestBoundary)       ! Western BC requested
      IF (.NOT. wbc) RETURN  ! -> no need for Western BC if we are not at Western border of global domain
      IF (.NOT. doEW) RETURN ! -> user did not request East/West BCs
    CASE(EastBoundary)       ! Eastern BC requested
      IF (.NOT. ebc) RETURN  ! -> no need for Eastern BC if we are not at Eastern border of global domain
      IF (.NOT. doEW) RETURN ! -> user did not request East/West BCs
    CASE(NorthBoundary)      ! Northern BC requested
      IF (.NOT. nbc) RETURN  ! -> no need for Northern BC if we are not at Northern border of global domain
      IF (.NOT. doNS) RETURN ! -> user did not request North/South BCs
    CASE(SouthBoundary)      ! Southern BC requested
      IF (.NOT. sbc) RETURN  ! -> no need for Southern BC if we are not at Southern border of global domain
      IF (.NOT. doNS) RETURN ! -> user did not request North/South BCs
    CASE (SouthWestCorner)
      IF (.NOT. doCorners) RETURN ! -> user did not request corners
      IF (.NOT. sbc) RETURN  ! -> no need for Southern corner if we are not at Southern border of global domain
      IF (.NOT. wbc) RETURN  ! -> no need for Western corner if we are not at Western border of global domain
    CASE (SouthEastCorner)
      IF (.NOT. doCorners) RETURN ! -> user did not request corners
      IF (.NOT. sbc) RETURN  ! -> no need for Southern corner if we are not at Southern border of global domain
      IF (.NOT. ebc) RETURN  ! -> no need for Eastern corner if we are not at Eastern border of global domain
    CASE (NorthWestCorner)
      IF (.NOT. doCorners) RETURN ! -> user did not request corners
      IF (.NOT. nbc) RETURN  ! -> no need for Northern corner if we are not at Northern border of global domain
      IF (.NOT. wbc) RETURN  ! -> no need for Western corner if we are not at Western border of global domain
    CASE (NorthEastCorner)
      IF (.NOT. doCorners) RETURN ! -> user did not request corners
      IF (.NOT. nbc) RETURN  ! -> no need for Northern corner if we are not at Northern border of global domain
      IF (.NOT. ebc) RETURN  ! -> no need for Eastern corner if we are not at Eastern border of global domain
    CASE DEFAULT
    CALL model_abort( my_cart_id, 42020, 'LBC: Invalid option in apply_bc encountered', 'apply_bc' ) 
  END SELECT

  ! we are actually sitting at the global domain border equal to
  ! the BC direction we are being queried for and should apply
  ! BCs on this side, so return TRUE
  apply_bc = .TRUE.

END FUNCTION apply_bc


!==============================================================================
!+ Module function get_bc_ranges to compute indices
!------------------------------------------------------------------------------

SUBROUTINE get_bc_ranges( bc_dir, nlines, ist, ien, jst, jen, wbc, ebc, nbc, sbc,  &
  doEW, doNS, doCorners, ibc_start, ibc_end, jbc_start, jbc_end )

!------------------------------------------------------------------------------
! Description:
!  Computes boundary range for a 2d block.
!
! Arguments:
!  bc_dir              , indicates which side of the subdomain of this PE we are considering
!  nlines              , number of lines at the boundary to apply bc
!  ist, ien            , start/end points of the inner domain
!  jst, jen
!  wbc, ebc, nbc, sbc  , flags to enable/disable bc at each boundary (W,E,N,S)
!  doEW                , if TRUE apply bc to WestEast boundaries
!  doNS                , if TRUE apply bc to NorthSouth boundaries
!  doCorners           , if TRUE apply bc to corners
!  ibc_start           , output values of boundary range
!  ibc_end
!  jbc_start
!  jbc_end 
!
! Note:
!  It is assumed that as a pre-condition for calling get_bc_ranges() the user has
!  called apply_bc() and that this function has returned TRUE. It can thus be assumed
!  that...   bc_dir == WestBoundary    requires    wbc = TRUE .AND. doEW = TRUE
!            bc_dir == EastBoundary    requires    ebc = TRUE .AND. doEW = TRUE
!            bc_dir == NorthBoundary   requires    nbc = TRUE .AND. doNS = TRUE
!            bc_dir == SouthBoundary   requires    sbc = TRUE .AND. doNS = TRUE
!
!------------------------------------------------------------------------------

  ! arguments
  INTEGER, INTENT(IN)    :: bc_dir
  INTEGER, INTENT(INOUT) :: ibc_start, ibc_end, jbc_start, jbc_end
  INTEGER, INTENT(IN)    :: nlines
  INTEGER, INTENT(IN)    :: ist, ien, jst, jen
  LOGICAL, INTENT(IN)    :: wbc, ebc, nbc, sbc
  LOGICAL, INTENT(IN)    :: doEW, doNS, doCorners

  ! this is a request for Eastern/Western BCs
  IF ( bc_dir == WestBoundary .OR.  bc_dir == EastBoundary ) THEN

    ! NOTE: corners are handled by default in the NS-direction, so here they
    !       are only added in the special cases where they have not already been
    !       consideren in the NS-direction

    ! no corners by default
    jbc_start = jst
    jbc_end = jen
    ! add corners if required
    IF ( doCorners ) THEN
      ! check if Southern corner has already been done in NS BCs, otherwise add
      IF ( .NOT. (sbc .AND. doNS) ) THEN
        jbc_start = jst - nlines
      END IF
      ! check if Northern corner has already been done in NS BCs, otherwise add
      IF ( .NOT. (nbc .AND. doNS) )  THEN
        jbc_end = jen + nlines
      END IF
    END IF
    ! setup EW indices to include nlines in the boundary zone
    IF ( bc_dir == WestBoundary ) THEN
      ibc_start = ist - nlines
      ibc_end = ist - 1
    ELSE IF (bc_dir == EastBoundary ) THEN
      ibc_start = ien + 1
      ibc_end = ien + nlines
    END IF

  ! this is a request for Northern/Southern BCs
  ELSE IF ( bc_dir == NorthBoundary .OR.  bc_dir == SouthBoundary ) THEN

    ! no corners by default
    ibc_start = ist
    ibc_end = ien
    ! add corners if required
    IF ( doCorners ) THEN
      ibc_start = ist - nlines
      ibc_end   = ien + nlines
    END IF

    ! setup NS indices to include nlines in the boundary zone
    IF ( bc_dir == NorthBoundary ) THEN
      jbc_start = jen + 1
      jbc_end = jen + nlines
    ELSE IF ( bc_dir == SouthBoundary ) THEN
      jbc_start = jst - nlines
      jbc_end = jst - 1
    END IF

  ! this is a request for South-West Corner
  ELSE IF ( bc_dir == SouthWestCorner ) THEN

    ibc_start = ist
    ibc_end   = ist - 1
    IF ( doCorners ) THEN
      ibc_start = ist - nlines
    END IF
    jbc_start = jst
    jbc_end   = jst - 1
    IF ( doCorners ) THEN
      jbc_start = jst - nlines
    END IF

  ! this is a request for South-East Corner
  ELSE IF ( bc_dir == SouthEastCorner ) THEN

    ibc_start = ien + 1
    ibc_end   = ien
    IF ( doCorners ) THEN
      ibc_end = ien + nlines
    END IF
    jbc_start = jst
    jbc_end   = jst - 1
    IF ( doCorners ) THEN
      jbc_start = jst - nlines
    END IF

  ! this is a request for Norht-West Corner
  ELSE IF ( bc_dir == NorthWestCorner ) THEN

    ibc_start = ist
    ibc_end   = ist - 1
    IF ( doCorners ) THEN
      ibc_start = ist - nlines
    END IF
    jbc_start = jen + 1
    jbc_end   = jen
    IF ( doCorners ) THEN
      jbc_end = jen + nlines
    END IF

  ! this is a request for North-East Corner
  ELSE IF ( bc_dir == NorthEastCorner ) THEN

    ibc_start = ien + 1
    ibc_end   = ien
    IF ( doCorners ) THEN
      ibc_end = ien + nlines
    END IF
    jbc_start = jen + 1
    jbc_end   = jen
    IF ( doCorners ) THEN
      jbc_end = jen + nlines
    END IF

  ELSE
    CALL model_abort( my_cart_id, 42030, 'LBC: Illegal condition', 'get_bc_ranges' )
  END IF

  ! safety bounds for indices
  ibc_start = MIN( ie, MAX( 1, ibc_start ) )
  ibc_end   = MIN( ie, MAX( 1, ibc_end   ) )
  jbc_start = MIN( je, MAX( 1, jbc_start ) )
  jbc_end   = MIN( je, MAX( 1, jbc_end   ) )

END SUBROUTINE get_bc_ranges


!==============================================================================
!+ Module functions to check whether we are at global border
!------------------------------------------------------------------------------

logical FUNCTION PE_at_western_boundary()
  IF ( my_cart_neigh(WestBoundary) == -1 ) THEN
    PE_at_western_boundary = .TRUE.
  ELSE
    PE_at_western_boundary = .FALSE.
  END IF

END FUNCTION PE_at_western_boundary


logical FUNCTION PE_at_southern_boundary()
  IF ( my_cart_neigh(SouthBoundary) == -1 ) THEN
    PE_at_southern_boundary = .TRUE.
  ELSE
    PE_at_southern_boundary = .FALSE.
  END IF

END FUNCTION PE_at_southern_boundary


logical FUNCTION PE_at_eastern_boundary()
  IF ( my_cart_neigh(EastBoundary) == -1 ) THEN
    PE_at_eastern_boundary = .TRUE.
  ELSE
    PE_at_eastern_boundary = .FALSE.
  END IF

END FUNCTION PE_at_eastern_boundary


logical FUNCTION PE_at_northern_boundary()
  IF ( my_cart_neigh(NorthBoundary) == -1 ) THEN
    PE_at_northern_boundary = .TRUE.
  ELSE
    PE_at_northern_boundary = .FALSE.
  END IF

END FUNCTION PE_at_northern_boundary

END MODULE src_lbc
