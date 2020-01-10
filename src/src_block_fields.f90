!+ Module src_block_fields
!------------------------------------------------------------------------------

MODULE src_block_fields

!------------------------------------------------------------------------------
!
! Description : 
!  Contains data structures and functions to automatically copy fields from ijk
!  format to block format (used in the common COSMO-ICON physics) and vice-versa
!
!  Correspondence between ijk and block fields is stored in the blockFieldTable 
!  data structure using pointers.
!
!  For each physics package lists of fields to be copied to block and back are 
!  stored in a structure of type CopylistStruct which contains index lists
!  referring to the blockFieldTable
!
!  The function requestCopy is used to set the copy flag at a given time step.
!
!  The routines copy_to_block and copy_from_block are used to copy all elements 
!  of a list to block, respectively from block to ijk. 
!  This should follow a call to request copy
!
!  Note : All copy routines use assumed shape array arguments. 
!         Arrays with a 3rd dimension of either 0:n-1 or 1:n would
!         be both treated in the same manner as an array with a 3rd dimension n 
!         (with indices from 1 to n).
! 
!  Convention: Private routine use camelCase, Public routine under_score
!
! Routines contained:
!
! - init_block_fields
!      Initialize blockFieldTable registry , set state
!
! - finalize_block_fields
!      Finalize, clean up data structures
!
! - registerBlockField3dRealTimelevel
! - registerBlockField3dReal
! - registerBlockField2dRealTimelevel
! - registerBlockField2dReal
! - registerBlockField2dInt
! - registerBlockField2dLogical
!      procedures to register field and corresponding block field into 
!      data structure
!
! - checkRegisterBlockFieldPointer
!      Helper function checkRegisterBlockField for register_block_field
!      Checks name, allocation and state returns errorcode and message
!
! - checkRegisterBlockFieldName
!      Helper function checkRegisterBlockFieldName for register_block_field
!      Check name is unique and status is appropriate
!
! - init_copy_list
!      Initialize CopyListStruct structure. 
!      This needs to be called before register_copy
!
! - registerCopy2dReal
! - registerCopy2dInt
! - registerCopy2dLogical
! - registerCopy3dReal
!      Procedures "register_copyX" to register required copies 
!
! - reset_copy_flags
!      Procedure to reset all copy flags in blockFieldTable
!
! - request_copy
!      Set copy to block, and back to ijk flags in blockFieldTable for given 
!      indices copy list structure CLS
!
! - setCopyToBlockFlag
! - setCopyFromBlockFlag
!      Helper function for request_copy 
!      Set copy flag for a list of indices
!
! - copy_to_block
! - copyToBlock3d
! - copyToBlock2d
! - copy_from_block
! - copyFromBlock3d
! - copyFromBlock2d
!      Procedures to apply requested copies
!
! - finalize_copy
!      Checks that all requested copied have been done, i.e. all lcopyToBlock
!       and lcopyFromBlock are FALSE
!
! - checkCopyBlockState
!      Check that current state machine is between ValidMin and Max
!
! - clearBlockFieldTable
!      Clear BlockFieldTable fields, nulify pointers
!
! Functions contained:
!
! - nextFreeBlockId
!      Scan blockFieldTable, return first unused id
!
! - getCopyType
!      Get copy type index based on ksize
!      Store ksize information in global array copyTypeKSize
!      Assumes init_block_struct was called once
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
!  Adaptations to MCH version
! V5_4         2016-03-10 Xavier Lapillonne
!  Deactivate OpenACC statements for the moment being
!  Restructured module as a preparation for future optimization
!  (multi-copy)
! V5_4a        2016-05-10 Xavier Lapillonne
!  Added support for integer and logical
!  Renamed internal type specific routines with trailing 
!  Real/Int/Logical
! V5_4b        2016-07-12 Ulrich Schaettler
!  Increased nksizeMax to 8 for implementation of blocked convection
!  Added kind of function to error message in registerCopy routines
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Modules used:

USE kind_parameters,       ONLY : wp

USE data_parallel,         ONLY : my_cart_id

USE environment,           ONLY : model_abort

!------------------------------------------------------------------------------

IMPLICIT NONE

!------------------------------------------------------------------------------

PRIVATE

!==============================================================================

!------------------------------------------------------------------------------
! Global (i.e. public) Declarations:

! Types
PUBLIC :: CopylistStruct

! Variables
PUBLIC :: copyToBlockF, copyFromBlockF, mind_ilon, mind_jlat

! Routines
PUBLIC :: init_block_fields, finalize_block_fields, init_copy_list,           &
          register_block_field, register_copy, reset_copy_flags,              &
          request_copy, copy_to_block, copy_from_block,                       &
          finalize_copy

!==============================================================================
! Parameters

INTEGER, PARAMETER :: blockFieldTableSize = 300 !max number of block fields
INTEGER, PARAMETER :: nksizeMax = 8             !max number of possible different
                                                !k dimensions

! Flag used to copy to block and back from block to ijk
INTEGER, PARAMETER :: copyToBlockF   = 1
INTEGER, PARAMETER :: copyFromBlockF = 2
INTEGER, PARAMETER :: ncopyFlags     = 2   ! Should exaclty correspond to largest
                                           ! icopyFlag

!Copy type information
! icopyTypeX : define field type
INTEGER, PARAMETER :: icopyTypeLogical2D = 1
INTEGER, PARAMETER :: icopyTypeInt2D     = 2
INTEGER, PARAMETER :: icopyTypeReal2D    = 3
INTEGER, PARAMETER :: icopyTypeReal3D    = 4
INTEGER, PARAMETER :: ncopyType          = 4  ! Should exaclty correspond to largest
                                              ! icopyTypeX value

! Possible states of src_block_fields module
! (Used to ensure calling order consistency)
INTEGER, PARAMETER :: &
  iStateStart    = -1, &    ! before initialization
  iStateInit     =  0, &    ! after initialization
  iStateRegField =  1, &    ! after field registration
  iStateRegCopy  =  2, &    ! after copy registration
  iStateCopy     =  3       ! request and copy

! Variable for debugging
LOGICAL, PARAMETER :: ldebug = .FALSE.

!==============================================================================
! Type definitions

! Structure holding pointers to field and block fields
TYPE block_struct

  ! Field information
  LOGICAL            :: lused          ! is this element attributed
  CHARACTER (LEN=10) :: yname          ! field name
  INTEGER            :: rank
  LOGICAL            :: timelevel      ! is there a time level
  INTEGER, POINTER   :: ntime          ! time level to copy
  INTEGER            :: ksize          ! 0 for 2d field

  ! Pointer to ijk-field
  REAL(KIND=wp), POINTER :: pfld3dRtime(:,:,:,:)
  REAL(KIND=wp), POINTER :: pfld3dR    (:,:,:)
  REAL(KIND=wp), POINTER :: pfld2dRtime(:,:,:)
  REAL(KIND=wp), POINTER :: pfld2dR    (:,:)
  INTEGER, POINTER       :: pfld2dI   (:,:)
  LOGICAL, POINTER       :: pfld2dL   (:,:)


  ! Pointer to block-field
  REAL(KIND=wp), POINTER :: pfld3dR_b  (:,:)
  REAL(KIND=wp), POINTER :: pfld2dR_b  (:)
  INTEGER, POINTER       :: pfld2dI_b (:)
  LOGICAL, POINTER       :: pfld2dL_b (:)

  ! Copy flags
  INTEGER :: copyType                ! used to pack similar copies together
  LOGICAL :: lcopyFlag(ncopyFlags)   ! .false. if field does not need to be copied


END TYPE block_struct

! Structure containing in (copytoblock) and out (copyback) variables lists
TYPE CopylistStruct
  INTEGER :: IN(blockFieldTableSize)     ! list of fields to be copied to block
  INTEGER :: ncopyin                     ! number of register copy in fields
  INTEGER :: OUT(blockFieldTableSize)    ! list of fields to be copied from block to ijk
  INTEGER :: ncopyout                    ! number of register copy out fields
END TYPE CopylistStruct

!==============================================================================
! Module variables

INTEGER :: nBlockFields = 0 ! actual number of registered block fields

INTEGER :: copyTypeNKSize(ncopyType)               ! number of different ksize for a given copytype
                                           ! (used to group field of same size in the same copy)
INTEGER :: copyTypeKSize(ncopyType,nksizeMax)      ! list of all different ksize of this copytype

! Current state
INTEGER :: icopyBlockState = iStateStart

! Table containing ijk and corresponding block field pointers
TYPE(block_struct) :: blockFieldTable(blockFieldTableSize)

! Fields for index-reorganization
INTEGER,  ALLOCATABLE ::    &
  mind_ilon(:,:),  & ! keeps the i (longitude) index
  mind_jlat(:,:)     ! keeps the j (latitude)  index
                     ! for every grid point in the (nproma,nblock) ordering

!==============================================================================
! Interface definition

INTERFACE register_copy
  MODULE PROCEDURE         &
    registerCopy3dReal,        &
    registerCopy2dReal,    &
    registerCopy2dInt,     &
    registerCopy2dLogical
END INTERFACE

INTERFACE register_block_field
  MODULE PROCEDURE                  &
    registerBlockField3dRealTimelevel,  &
    registerBlockField3dReal,           &
    registerBlockField2dRealTimelevel,  &
    registerBlockField2dReal,           &
    registerBlockField2dInt,        &
    registerBlockField2dLogical
END INTERFACE

!==============================================================================
! Module procedures in src_block_fields
!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure "init_block_fields" in "src_block_fields"
!------------------------------------------------------------------------------

SUBROUTINE init_block_fields(isc, iec, jsc, jec, nproma, nlastproma, nblock, &
                             ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   Initialize blockFieldTable registry , set state
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  INTEGER, INTENT(IN) :: isc, iec, jsc, jec   !  compute indices in lat-long directions
  INTEGER, INTENT(IN) :: nproma, nlastproma, nblock ! indices in block format
  INTEGER, INTENT(OUT) :: ierror
  CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg

! Locals:
! --------------------
  INTEGER :: i, j, ip, ipend, ib, nij
  INTEGER :: izerr


!------------------------------------------------------------------------------
! Begin Subroutine init_block_fields
!------------------------------------------------------------------------------

  ierror = 0
  yerrmsg = ''

  ! Check and set state
  ! Allowed state :  iStateStart
  CALL checkCopyBlockState(iStateStart, iStateStart, ierror, yerrmsg)
  IF (ierror > 0) RETURN

  ! Set state flag to init value
  icopyBlockState = iStateInit

  ! Init block field table
  CALL clearBlockFieldTable()

  ! Init copyTypes
  copyTypeNKSize(:)  = 0
  copyTypeKSize(:,:) = -1

  ! Allocate and initialize mind_ilon, mind_jlat
  ! compute the index-arrays

  ! Check that nproma and nblock have consistent sizes
  nij = (iec - isc + 1)*(jec - jsc + 1)
  IF ( (nproma < 1) .OR. (nproma > nij) .OR. &
       (nlastproma < 1) .OR. (nlastproma > nproma) .OR. &
       (nblock < 1) .OR. (nblock > nij) ) THEN
    ierror = 1034
    WRITE(yerrmsg, "(A,4I8)") 'ERROR: Wrong values for nproma, nlastproma, nblock, nij', &
          nproma, nlastproma, nblock, nij
    RETURN
  ENDIF

  ! Allocation
  ALLOCATE (mind_ilon(nproma, nblock), mind_jlat(nproma, nblock), STAT=izerr)
  IF (izerr>0) THEN
    ierror = 10135
    WRITE(yerrmsg, "(A,I8)") 'ERROR: Index arrays allocation failed , status=', izerr
    RETURN
  END IF

!$acc enter data &
!$acc create (mind_ilon, mind_jlat)

  ! Set index
  i = isc - 1
  j = jsc

  ipend = nproma       ! compute domain for all but last block is 1:nproma
  DO ib = 1, nblock
    IF (ib == nblock) ipend = nlastproma          ! last block 1:nlastproma
    DO ip = 1, ipend
      i = i + 1
      IF (i > iec) THEN
        j = j + 1
        i = isc
      ENDIF
      IF (j > jec) THEN
        ! Error: this must not happen, because
        ! nproma*nblock <= (iec-isc+1)*(jec-jsc+1)
        ierror = 10136
        WRITE(yerrmsg, "(A,I8,I8,I8)") ' ERROR: inconsistent sizes nproma*nblock, nlastproma, nij :', &
             nproma*nblock, nlastproma, nij
        RETURN
      ENDIF

      mind_ilon (ip,ib) = i
      mind_jlat (ip,ib) = j
    ENDDO
  ENDDO

!$acc update device (mind_ilon, mind_jlat)

!------------------------------------------------------------------------------
! End of module procedure init_block_fields
!------------------------------------------------------------------------------

END SUBROUTINE init_block_fields

!==============================================================================
!==============================================================================
!+ Module procedure "init_block_fields" in "src_block_fields"
!------------------------------------------------------------------------------

SUBROUTINE finalize_block_fields(ierror)

!------------------------------------------------------------------------------
!
! Description:
!   Finalize, clean up data structures
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  INTEGER, INTENT(OUT) :: ierror

!----------- End of header ----------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine finalize_block_fields
!------------------------------------------------------------------------------

  ierror = 0

  ! Reset state (callable from any state)
  icopyBlockState = iStateStart

  ! Clean up
  CALL clearBlockFieldTable()

  ! Init copyTypes
  copyTypeNKSize(:)   = 0
  copyTypeKSize(:,:) = 0

  ! Deallocate index-arrays

  !$acc exit data &
  !$acc delete (mind_ilon, mind_jlat)

  DEALLOCATE(mind_ilon, mind_jlat, STAT=ierror)

!------------------------------------------------------------------------------
! End of module procedure finalize_block_fields
!------------------------------------------------------------------------------

END SUBROUTINE finalize_block_fields

!==============================================================================
!==============================================================================
!+ register_block_field : 3d + timelevel field
!------------------------------------------------------------------------------    

SUBROUTINE registerBlockField3dRealTimelevel(yname, field, field_b, ntime)

!------------------------------------------------------------------------------
!
! Description:
!  Module procedures "register_block_fieldX" in "src_block_fields" to
!  register field and corresponding block field into data structure
!   The following routines store pointers to field and corresponding block 
!   field in table blockFieldTable
!   Fields with a time level dimension requires to be given the appropiate
!   level. The time level variables given in argument needs to have a scope 
!   larger than all calls to copy to block as a pointer to it is used.
!   It is an error to register fields which are not allocated
!   This routine may abort program
!
! Subroutines :
!   registerBlockField3dRealTimelevel
!   registerBlockField3dReal
!   registerBlockField2dRealTimelevel
!   registerBlockField2dReal
!   registerBlockField2dInt
!   registerBlockField2dLogical
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------

  CHARACTER (LEN=*), INTENT(IN) :: yname
  REAL(KIND=wp), TARGET, INTENT(IN) :: field(:,:,:,:)
  REAL(KIND=wp), TARGET, INTENT(IN) :: field_b(:,:)
  INTEGER, TARGET, INTENT(IN) :: ntime

! Local parameters:
! ----------------
  INTEGER :: id, ksize
  INTEGER :: ierror
  CHARACTER(LEN=100) :: yerrmsg

!----------- End of header ----------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine registerBlockField3dRealTimelevel
!------------------------------------------------------------------------------

  ierror = 0
  yerrmsg = ''

  ! Check yname is unique and status is correct
  CALL checkRegisterBlockFieldName(yname, ierror, yerrmsg)
  IF (ierror > 0)  &
    CALL model_abort(my_cart_id, ierror, yerrmsg, 'checkRegisterBlockFieldName (3dRealTimelevel)')
  icopyBlockState = iStateRegField

  id = nextFreeBlockId()
  blockFieldTable(id)%lused = .TRUE.

  ksize = SIZE(field, 3)
  blockFieldTable(id)%yname      = TRIM(yname)
  blockFieldTable(id)%rank       = 3
  blockFieldTable(id)%timelevel  = .TRUE.
  blockFieldTable(id)%ntime      => ntime
  blockFieldTable(id)%ksize      = ksize
  blockFieldTable(id)%copyType   = icopyTypeReal3D
  blockFieldTable(id)%pfld3dRtime => field
  blockFieldTable(id)%pfld3dR_b   => field_b

  CALL setCopyTypeKSize(icopyTypeReal3D, ksize)

  ! Check that pointer are associated, and dimension consistent
  ! This would trigger an error in case field or field_b
  ! are not allocated
  CALL checkRegisterBlockFieldPointer(yname, ASSOCIATED(blockFieldTable(id)%pfld3dRtime), &
       ASSOCIATED(blockFieldTable(id)%pfld3dR_b), SIZE(field, 4), ntime, ksize,           &
       SIZE(field_b, 2), ierror, yerrmsg)
  IF (ierror > 0)  &
    CALL model_abort(my_cart_id, ierror, yerrmsg, 'checkRegisterBlockFieldPointer (3dRealTimeLevel)')

!------------------------------------------------------------------------------
! End of module procedure registerBlockField3dRealTimelevel
!------------------------------------------------------------------------------

END SUBROUTINE registerBlockField3dRealTimelevel

!==============================================================================
!==============================================================================
! register_block_field : 3d field
!------------------------------------------------------------------------------

SUBROUTINE registerBlockField3dReal(yname, field, field_b)

!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  CHARACTER (LEN=*), INTENT(IN) :: yname
  REAL(KIND=wp), TARGET, INTENT(IN) :: field(:,:,:)
  REAL(KIND=wp), TARGET, INTENT(IN) :: field_b(:,:)

! Local parameters:
! ----------------
  INTEGER :: id, ksize
  INTEGER :: ierror
  CHARACTER(LEN=64) :: yerrmsg

!----------- End of header ----------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine registerBlockField3dReal
!------------------------------------------------------------------------------

  ierror = 0
  yerrmsg = ''

  ! Check yname is unique and status is correct
  CALL checkRegisterBlockFieldName(yname, ierror, yerrmsg)
  IF (ierror > 0)  &
    CALL model_abort(my_cart_id, ierror, yerrmsg, 'checkRegisterBlockFieldName (3dReallevel)')
  icopyBlockState = iStateRegField


  id = nextFreeBlockId()
  blockFieldTable(id)%lused = .TRUE.

  ksize = SIZE(field, 3)

  blockFieldTable(id)%yname      = TRIM(yname)
  blockFieldTable(id)%rank       = 3
  blockFieldTable(id)%timelevel  = .FALSE.
  blockFieldTable(id)%ksize      = ksize
  blockFieldTable(id)%copyType   = icopyTypeReal3D
  blockFieldTable(id)%pfld3dR     => field
  blockFieldTable(id)%pfld3dR_b   => field_b

  CALL setCopyTypeKSize(icopyTypeReal3D, ksize)

  ! Check that pointer are associated, and dimension consistent
  ! This would trigger an error in case field or field_b
  ! are not allocated
  CALL checkRegisterBlockFieldPointer(yname, ASSOCIATED(blockFieldTable(id)%pfld3dR),     &
         ASSOCIATED(blockFieldTable(id)%pfld3dR_b), 0, 0, ksize, SIZE(field_b, 2), ierror,&
         yerrmsg)
  IF (ierror > 0)  &
    CALL model_abort(my_cart_id, ierror, yerrmsg, 'checkRegisterBlockFieldPointer (3dReal)')

!------------------------------------------------------------------------------
! End of module procedure registerBlockField3dReal
!------------------------------------------------------------------------------

END SUBROUTINE registerBlockField3dReal

!==============================================================================
!==============================================================================
! register_block_field : 2d + timelevel field
!------------------------------------------------------------------------------

SUBROUTINE registerBlockField2dRealTimelevel(yname, field, field_b, ntime)

!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  CHARACTER (LEN=*), INTENT(IN) :: yname
  REAL(KIND=wp), TARGET, INTENT(IN) :: field(:,:,:)
  REAL(KIND=wp), TARGET, INTENT(IN) :: field_b(:)
  INTEGER, TARGET :: ntime

! Local parameters:
! ----------------
  INTEGER :: id
  INTEGER :: ierror
  CHARACTER(LEN=64) :: yerrmsg

!----------- End of header ----------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine registerBlockField2dRealTimelevel
!------------------------------------------------------------------------------

  ierror = 0
  yerrmsg = ''

  ! Check yname is unique and status is correct
  CALL checkRegisterBlockFieldName(yname, ierror, yerrmsg)
  IF (ierror > 0)   &
    CALL model_abort(my_cart_id, ierror, yerrmsg, 'checkRegisterBlockFieldName (2dRealTimeLevel)')
  icopyBlockState = iStateRegField

  id = nextFreeBlockId()
  blockFieldTable(id)%lused = .TRUE.

  blockFieldTable(id)%yname      = TRIM(yname)
  blockFieldTable(id)%rank       = 2
  blockFieldTable(id)%timelevel  = .TRUE.
  blockFieldTable(id)%ntime      => ntime
  blockFieldTable(id)%ksize      = 0
  blockFieldTable(id)%copyType   = icopyTypeReal2D
  blockFieldTable(id)%pfld2dRtime => field
  blockFieldTable(id)%pfld2dR_b   => field_b

  CALL setCopyTypeKSize(icopyTypeReal2D, 0)

  ! Check that pointer are associated, and dimension consistent
  ! This would trigger an error in case field or field_b
  ! are not allocated
  CALL checkRegisterBlockFieldPointer(yname, ASSOCIATED(blockFieldTable(id)%pfld2dRtime), &
       ASSOCIATED(blockFieldTable(id)%pfld2dR_b), SIZE(field, 3), ntime, 0, 0, ierror,    &
       yerrmsg)
  IF (ierror > 0)  &
       CALL model_abort(my_cart_id, ierror, yerrmsg, 'checkRegisterBlockFieldPointer (2dRealTimeLevel)')

!------------------------------------------------------------------------------
! End of module procedure registerBlockField2dRealTimelevel
!------------------------------------------------------------------------------

END SUBROUTINE registerBlockField2dRealTimelevel

!==============================================================================
!==============================================================================
! register_block_field : 2d
!------------------------------------------------------------------------------

SUBROUTINE registerBlockField2dReal(yname, field, field_b)

!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  CHARACTER (LEN=*), INTENT(IN) :: yname
  REAL(KIND=wp), TARGET, INTENT(IN) :: field(:,:)
  REAL(KIND=wp), TARGET, INTENT(IN) :: field_b(:)

! Local parameters:
! ----------------
  INTEGER :: id
  INTEGER :: ierror
  CHARACTER(LEN=64) :: yerrmsg


!----------- End of header ----------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine registerBlockField2dReal
!------------------------------------------------------------------------------

  ierror = 0
  yerrmsg = ''

  ! Check yname is unique and status is correct
  CALL checkRegisterBlockFieldName(yname, ierror, yerrmsg)
  IF (ierror > 0)  &
    CALL model_abort(my_cart_id, ierror, yerrmsg, 'checkRegisterBlockFieldName (2dReal)')
  icopyBlockState = iStateRegField

  id = nextFreeBlockId()
  blockFieldTable(id)%lused = .TRUE.

  blockFieldTable(id)%yname      = yname
  blockFieldTable(id)%rank       = 2
  blockFieldTable(id)%timelevel  = .FALSE.
  blockFieldTable(id)%ksize      = 0
  blockFieldTable(id)%copyType   = icopyTypeReal2D
  blockFieldTable(id)%pfld2dR     => field
  blockFieldTable(id)%pfld2dR_b   => field_b

  CALL setCopyTypeKSize(icopyTypeReal2D, 0)

  ! Check that pointer are associated, and dimension consistent
  ! This would trigger an error in case field or field_b
  ! are not allocated
  CALL checkRegisterBlockFieldPointer(yname, ASSOCIATED(blockFieldTable(id)%pfld2dR), &
       ASSOCIATED(blockFieldTable(id)%pfld2dR_b), 0, 0, 0, 0, ierror, yerrmsg)
  IF (ierror > 0)  &
       CALL model_abort(my_cart_id, ierror, yerrmsg, 'checkRegisterBlockFieldPointer (2dReal)')

!------------------------------------------------------------------------------
! End of module procedure registerBlockField2dReal
!------------------------------------------------------------------------------

END SUBROUTINE registerBlockField2dReal

!==============================================================================
!==============================================================================
! register_block_field : 2d integer
!------------------------------------------------------------------------------

SUBROUTINE registerBlockField2dInt(yname, field, field_b)

!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  CHARACTER (LEN=*), INTENT(IN) :: yname
  INTEGER, TARGET, INTENT(IN) :: field(:,:)
  INTEGER, TARGET, INTENT(IN) :: field_b(:)

! Local parameters:
! ----------------
  INTEGER :: id
  INTEGER :: ierror
  CHARACTER(LEN=64) :: yerrmsg


!----------- End of header ----------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine registerBlockField2dInt
!------------------------------------------------------------------------------

  ierror = 0
  yerrmsg = ''

  ! Check yname is unique and status is correct
  CALL checkRegisterBlockFieldName(yname, ierror, yerrmsg)
  IF (ierror > 0)  &
    CALL model_abort(my_cart_id, ierror, yerrmsg, 'checkRegisterBlockFieldName (2dInt)')
  icopyBlockState = iStateRegField

  id = nextFreeBlockId()
  blockFieldTable(id)%lused = .TRUE.

  blockFieldTable(id)%yname      = yname
  blockFieldTable(id)%rank       = 2
  blockFieldTable(id)%timelevel  = .FALSE.
  blockFieldTable(id)%ksize      = 0
  blockFieldTable(id)%copyType   = icopyTypeInt2D
  blockFieldTable(id)%pfld2dI    => field
  blockFieldTable(id)%pfld2dI_b  => field_b

  CALL setCopyTypeKSize(icopyTypeInt2D, 0)

  ! Check that pointer are associated, and dimension consistent
  ! This would trigger an error in case field or field_b
  ! are not allocated
  CALL checkRegisterBlockFieldPointer(yname, ASSOCIATED(blockFieldTable(id)%pfld2dI), &
       ASSOCIATED(blockFieldTable(id)%pfld2dI_b), 0, 0, 0, 0, ierror, yerrmsg)
  IF (ierror > 0)  &
       CALL model_abort(my_cart_id, ierror, yerrmsg, 'checkRegisterBlockFieldPointer (2dInt)')

!------------------------------------------------------------------------------
! End of module procedure registerBlockField2dInt
!------------------------------------------------------------------------------

END SUBROUTINE registerBlockField2dInt

!==============================================================================
!==============================================================================
! register_block_field : 2d logical
!------------------------------------------------------------------------------

SUBROUTINE registerBlockField2dLogical(yname, field, field_b)

!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  CHARACTER (LEN=*), INTENT(IN) :: yname
  LOGICAL, TARGET, INTENT(IN) :: field(:,:)
  LOGICAL, TARGET, INTENT(IN) :: field_b(:)

! Local parameters:
! ----------------
  INTEGER :: id
  INTEGER :: ierror
  CHARACTER(LEN=64) :: yerrmsg


!----------- End of header ----------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine registerBlockField2dLogical
!------------------------------------------------------------------------------

  ierror = 0
  yerrmsg = ''

  ! Check yname is unique and status is correct
  CALL checkRegisterBlockFieldName(yname, ierror, yerrmsg)
  IF (ierror > 0)  &
    CALL model_abort(my_cart_id, ierror, yerrmsg, 'checkRegisterBlockFieldName (2dLogical)')
  icopyBlockState = iStateRegField

  id = nextFreeBlockId()
  blockFieldTable(id)%lused = .TRUE.

  blockFieldTable(id)%yname      = yname
  blockFieldTable(id)%rank       = 2
  blockFieldTable(id)%timelevel  = .FALSE.
  blockFieldTable(id)%ksize      = 0
  blockFieldTable(id)%copyType   = icopyTypeLogical2D
  blockFieldTable(id)%pfld2dL     => field
  blockFieldTable(id)%pfld2dL_b   => field_b

  CALL setCopyTypeKSize(icopyTypeLogical2D, 0)

  ! Check that pointer are associated, and dimension consistent
  ! This would trigger an error in case field or field_b
  ! are not allocated
  CALL checkRegisterBlockFieldPointer(yname, ASSOCIATED(blockFieldTable(id)%pfld2dL), &
       ASSOCIATED(blockFieldTable(id)%pfld2dL_b), 0, 0, 0, 0, ierror, yerrmsg)
  IF (ierror > 0)  &
       CALL model_abort(my_cart_id, ierror, yerrmsg, 'checkRegisterBlockFieldPointer (2dLogical)')

!------------------------------------------------------------------------------
! End of module procedure registerBlockField2dLogical
!------------------------------------------------------------------------------

END SUBROUTINE registerBlockField2dLogical

!==============================================================================
!==============================================================================
! Helper function checkRegisterBlockField for register_block_field
!------------------------------------------------------------------------------

SUBROUTINE checkRegisterBlockFieldPointer(yname, lallocField, lallocfield_b,   &
                          nsizetime, ntime, ksize, ksize_b, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   Checks yname, allocation and state
!   Returns errorcode and message
!
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  CHARACTER (LEN=*), INTENT(IN) :: yname
  LOGICAL, INTENT(IN) :: lallocField, lallocField_b
  INTEGER, INTENT(IN) :: nsizetime, ntime, ksize, ksize_b
  INTEGER, INTENT(OUT) :: ierror
  CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg

! Local parameters:
! ----------------
  INTEGER :: id

  ierror = 0
  yerrmsg = ''
  id = 1

!------------------------------------------------------------------------------
! Begin Subroutine registerBlockFieldPointer
!------------------------------------------------------------------------------

  ! Check allocation
  IF (.NOT. lallocField) THEN
    ierror = 10102
    WRITE(yerrmsg, "(A,A)") 'ERROR: Field not allocated/associated : ', yname
    RETURN
  END IF

  IF (.NOT. lallocField_b) THEN
    ierror = 10103
    WRITE(yerrmsg, "(A,A)") 'ERROR: Block Field not allocated/associated : ', yname
    RETURN
  END IF

  ! Check time dimension
  IF (nsizetime > 0) THEN
    IF (ntime > nsizetime) THEN
      ierror = 10104
      WRITE(yerrmsg, "(A,A)") 'ERROR: Time dimension larger than actual size for field : ',  yname
      RETURN
    END IF
  END IF

  ! Check k dimension
  IF (ksize /= ksize_b) THEN
    ierror = 10105
    WRITE(yerrmsg, "(A,A,I8,I8)") 'ERROR: inconsistent k dimension field : ',  yname, ksize, ksize_b
    RETURN
  END IF

  IF (ldebug.AND.(my_cart_id==0)) WRITE(*,"(A,A16,I8,I8)") 'Register block field, k, ntz : ', TRIM(yname), ksize, nsizetime

END SUBROUTINE checkRegisterBlockFieldPointer

!==============================================================================
!==============================================================================
! Helper function checkRegisterBlockFieldName for register_block_field
! Check yname is unique and status is appropriate
!------------------------------------------------------------------------------

SUBROUTINE checkRegisterBlockFieldName(yname, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   Checks yname, allocation and state
!   Returns errorcode and message. 
!   Returns 0 if fieldname is not alread used
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  CHARACTER (LEN=*), INTENT(IN) :: yname
  INTEGER, INTENT(OUT) :: ierror
  CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg

! Local parameters:
! ----------------
  INTEGER :: id

  ierror = 0
  yerrmsg = ''
  id = 1

!------------------------------------------------------------------------------
! Begin Subroutine checkRegisterBlockFieldName
!------------------------------------------------------------------------------

  ! Check State
  CALL checkCopyBlockState(iStateInit, iStateRegField, ierror, yerrmsg)
  IF (ierror > 0) RETURN

  ! Check yname is not already used
  DO WHILE( blockFieldTable(id)%lused )
    IF (TRIM(blockFieldTable(id)%yname) == TRIM(yname)) THEN
      ierror = 10101
      WRITE(yerrmsg, "(A,A)") 'ERROR: Field name already used : ', yname
      RETURN
    END IF
    id = id + 1
    IF (id > blockFieldTableSize) EXIT
  END DO

END SUBROUTINE checkRegisterBlockFieldName

!==============================================================================
!==============================================================================
! Helper function nextFreeBlockId for register_block_field
!------------------------------------------------------------------------------

FUNCTION nextFreeBlockId()

!------------------------------------------------------------------------------
!
! Description:
!   Scan blockFieldTable, return first unused id
!   Side effect: increment global counter
!
!------------------------------------------------------------------------------

  IMPLICIT NONE
  
! Return value
  INTEGER :: nextFreeBlockId
  
! Locals
  INTEGER :: id, ierror
  CHARACTER(LEN=64) :: yerrmsg

!------------------------------------------------------------------------------
! Begin function nextFreeBlockId
!------------------------------------------------------------------------------

  ! initialize ID
  id = 1

  ! find free position in block-field table
  DO WHILE( blockFieldTable(id)%lused )
     id = id + 1
     IF (id > blockFieldTableSize) EXIT
  END DO

  ! error if id > blockFieldTableSize
  IF (id > blockFieldTableSize) THEN
    ierror = 10108
    yerrmsg = 'No slot available. Need to increase blockFieldTableSize'
    CALL model_abort(my_cart_id, ierror, yerrmsg, 'nextFreeBlockId')
  END IF

  nBlockFields = id    ! increment global counter
  nextFreeBlockId = id ! return nextFreeBlockId

END FUNCTION nextFreeBlockId

!==============================================================================
!==============================================================================
! Helper routine setCopyTypeKSize for register_block_field
!------------------------------------------------------------------------------

SUBROUTINE setCopyTypeKSize(icopyType, ksize)

!------------------------------------------------------------------------------
!
! Description:
!   Increase copyTypeNKSize if this ksize do not already exist for this type
!   do nothing otherwise
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  INTEGER, INTENT(IN) :: icopyType, ksize

! Local parameters:
! ----------------
  INTEGER :: ik, nksize
  LOGICAL :: lexistKSize
  INTEGER :: ierror
  CHARACTER(LEN=64) :: yerrmsg

  lexistKSize = .FALSE.

  ! Check if this copyTypeKSize already exists
  DO ik = 1, copyTypeNKSize(icopyType)
    IF ( copyTypeKSize(icopytype, ik) == ksize )  &
          lexistKSize = .TRUE.
  END DO

  ! Increment NcopyTypeksize if new ksize
  IF (.NOT.lexistKSize) THEN
    IF ( copyTypeNKSize(icopyType) < nksizeMax ) THEN
      nksize = copyTypeNKSize(icopyType) + 1
      copyTypeNKSize(icopyType) = nksize
      copyTypeKSize(icopytype, nksize) = ksize
    ELSE
      ierror = 10109
      yerrmsg = 'ERROR: Need to increase nksizeMax'
      CALL model_abort(my_cart_id, ierror, yerrmsg, 'setCopyTypeKSize')
    END IF
  END IF

END SUBROUTINE setCopyTypeKSize

!==============================================================================
!==============================================================================
!+ Module procedure "init_copy_list" in "src_block_fields"
!------------------------------------------------------------------------------

SUBROUTINE init_copy_list(CList)

!------------------------------------------------------------------------------
!
! Description:
!   Initialize CopyListStruct structure. This needs to be called
!   before register_copy
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------

  TYPE(CopyListStruct), INTENT(INOUT) :: CList

!------------------------------------------------------------------------------
! Begin function init_copy_list
!------------------------------------------------------------------------------

  Clist%IN(:)      = -1
  Clist%ncopyin    =  0

  Clist%OUT(:)     = -1
  Clist%ncopyout   =  0

END SUBROUTINE init_copy_list

!==============================================================================
!==============================================================================
!+ Module procedures "register_copyX" in "src_block_fields" to
!  register required copies
!------------------------------------------------------------------------------
SUBROUTINE registerCopy2dReal(Field_b, CLS, copyFlag)

!------------------------------------------------------------------------------
!
! Description:
!
!   Find field id in blockFieldTable and add it to the
!   copy to block (IN) or copy from block (OUT) list depending
!   on the CopyFlag argument
!   3d fields are added to list3d depending on copyType
!   2d fields are added to list2d
!   Corresponding ncopy variable is incremented
!   Assumes CL%List is once initialized to -1
!   It is an error to register copy for fields which are not registered
!   This routine may abort model run
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  REAL(KIND=wp), INTENT(IN), TARGET :: Field_b(:)
  TYPE(CopylistStruct), INTENT(INOUT), TARGET :: CLS
  INTEGER, INTENT(IN) :: copyFlag

! Local parameters:
! ----------------

  INTEGER, POINTER :: CL(:), ncopy
  INTEGER :: id
  LOGICAL :: lsuccess
  INTEGER :: ierror
  CHARACTER(LEN=64) :: yerrmsg
  CHARACTER(LEN=64) :: yroutine
  CHARACTER(LEN=32) :: yinfo

!------------------------------------------------------------------------------
! Begin subroutine registerCopy2dReal
!------------------------------------------------------------------------------

  yroutine='register_copy (2dReal)'

  ! Check State : allowed iStateRegField or iStateRegCopy
  CALL checkCopyBlockState(iStateRegField, iStateRegCopy, ierror, yerrmsg)
  IF (ierror > 0) CALL model_abort(my_cart_id, ierror, yerrmsg, yroutine)
  icopyBlockState = iStateRegCopy

  ! Set copy list (CL) and ncopy pointers to IN or OUT depending on the copy flag
  CALL setCopyListPointer(CLS, copyFlag, CL, ncopy, yinfo, ierror, yerrmsg)
  IF (ierror > 0) CALL model_abort(my_cart_id, ierror, yerrmsg, yroutine)

  ! Scan throught blockFieldTable
  id = 1
  lsuccess = .FALSE.
  DO WHILE( (.NOT. lsuccess) .AND. (id <= nBlockFields) )
    IF ( ASSOCIATED(blockFieldTable(id)%pfld2dR_b, Field_b) ) THEN
      ncopy = ncopy + 1                          ! increment counter
      IF (ncopy > nBlockFields) THEN
        WRITE(yerrmsg, "(A,I8,I8)") 'ERROR: ncopy> nBlockFields', &
             ncopy, nBlockFields
        ierror = 10111
        CALL model_abort(my_cart_id, ierror, yerrmsg, yroutine)
      END IF
      CL(ncopy) = id
      lsuccess = .TRUE.
      IF (ldebug.AND.(my_cart_id == 0)) WRITE(*,"(A,A)") TRIM(yinfo), &
           TRIM(blockFieldTable(id)%yname)
    END IF
    id = id + 1
  END DO

  IF ( .NOT. lsuccess ) THEN ! no corresponding field
    WRITE(yerrmsg, "(A)") 'ERROR: no corresponding field in  blockFieldTable (2dReal)'
    ierror = 10112
    CALL model_abort(my_cart_id, ierror, yerrmsg, yroutine)
  END IF

END SUBROUTINE registerCopy2dReal


!==============================================================================

SUBROUTINE registerCopy2dInt(Field_b, CLS, copyFlag)

!------------------------------------------------------------------------------
! Begin subroutine registerCopy2dInt
!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  INTEGER, INTENT(IN), TARGET :: Field_b(:)
  TYPE(CopylistStruct), INTENT(INOUT), TARGET :: CLS
  INTEGER, INTENT(IN) :: copyFlag

! Local parameters:
! ----------------

  INTEGER, POINTER :: CL(:), ncopy
  INTEGER :: id
  LOGICAL :: lsuccess
  INTEGER :: ierror
  CHARACTER(LEN=64) :: yerrmsg
  CHARACTER(LEN=64) :: yroutine
  CHARACTER(LEN=32) :: yinfo

!------------------------------------------------------------------------------
! Begin subroutine registerCopy2dInt
!------------------------------------------------------------------------------

  yroutine='register_copy (2dInt)'

  ! Check State : allowed iStateRegField or iStateRegCopy
  CALL checkCopyBlockState(iStateRegField, iStateRegCopy, ierror, yerrmsg)
  IF (ierror > 0) CALL model_abort(my_cart_id, ierror, yerrmsg, yroutine)
  icopyBlockState = iStateRegCopy

  ! Set copy list (CL) and ncopy pointers to IN or OUT depending on the copy flag
  CALL setCopyListPointer(CLS, copyFlag, CL, ncopy, yinfo, ierror, yerrmsg)
  IF (ierror > 0) CALL model_abort(my_cart_id, ierror, yerrmsg, yroutine)

  ! Scan throught blockFieldTable
  id = 1
  lsuccess = .FALSE.
  DO WHILE( (.NOT. lsuccess) .AND. (id <= nBlockFields) )
    IF ( ASSOCIATED(blockFieldTable(id)%pfld2dI_b, Field_b) ) THEN
      ncopy = ncopy + 1                          ! increment counter
      IF (ncopy > nBlockFields) THEN
        WRITE(yerrmsg, "(A,I8,I8)") 'ERROR: ncopy> nBlockFields', &
             ncopy, nBlockFields
        ierror = 10111
        CALL model_abort(my_cart_id, ierror, yerrmsg, yroutine)
      END IF
      CL(ncopy) = id
      lsuccess = .TRUE.
      IF (ldebug.AND.(my_cart_id == 0)) WRITE(*,"(A,A)") TRIM(yinfo), &
           TRIM(blockFieldTable(id)%yname)
    END IF
    id = id + 1
  END DO

  IF ( .NOT. lsuccess ) THEN ! no corresponding field
    WRITE(yerrmsg, "(A)") 'ERROR: no corresponding field in  blockFieldTable (2dInt)'
    ierror = 10112
    CALL model_abort(my_cart_id, ierror, yerrmsg, yroutine)
  END IF

END SUBROUTINE registerCopy2dInt

!==============================================================================
SUBROUTINE registerCopy2dLogical(Field_b, CLS, copyFlag)

!------------------------------------------------------------------------------
! Begin subroutine registerCopy2dInt
!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  LOGICAL, INTENT(IN), TARGET :: Field_b(:)
  TYPE(CopylistStruct), INTENT(INOUT), TARGET :: CLS
  INTEGER, INTENT(IN) :: copyFlag

! Local parameters:
! ----------------

  INTEGER, POINTER :: CL(:), ncopy
  INTEGER :: id
  LOGICAL :: lsuccess
  INTEGER :: ierror
  CHARACTER(LEN=64) :: yerrmsg
  CHARACTER(LEN=64) :: yroutine
  CHARACTER(LEN=32) :: yinfo

!------------------------------------------------------------------------------
! Begin subroutine registerCopy2dLogical
!------------------------------------------------------------------------------

  yroutine='register_copy (2dLogical)'

  ! Check State : allowed iStateRegField or iStateRegCopy
  CALL checkCopyBlockState(iStateRegField, iStateRegCopy, ierror, yerrmsg)
  IF (ierror > 0) CALL model_abort(my_cart_id, ierror, yerrmsg, yroutine)
  icopyBlockState = iStateRegCopy

  ! Set copy list (CL) and ncopy pointers to IN or OUT depending on the copy flag
  CALL setCopyListPointer(CLS, copyFlag, CL, ncopy, yinfo, ierror, yerrmsg)
  IF (ierror > 0) CALL model_abort(my_cart_id, ierror, yerrmsg, yroutine)

  ! Scan throught blockFieldTable
  id = 1
  lsuccess = .FALSE.
  DO WHILE( (.NOT. lsuccess) .AND. (id <= nBlockFields) )
    IF ( ASSOCIATED(blockFieldTable(id)%pfld2dL_b, Field_b) ) THEN
      ncopy = ncopy + 1                          ! increment counter
      IF (ncopy > nBlockFields) THEN
        WRITE(yerrmsg, "(A,I8,I8)") 'ERROR: ncopy> nBlockFields', &
             ncopy, nBlockFields
        ierror = 10111
        CALL model_abort(my_cart_id, ierror, yerrmsg, yroutine)
      END IF
      CL(ncopy) = id
      lsuccess = .TRUE.
      IF (ldebug.AND.(my_cart_id == 0)) WRITE(*,"(A,A)") TRIM(yinfo), &
           TRIM(blockFieldTable(id)%yname)
    END IF
    id = id + 1
  END DO

  IF ( .NOT. lsuccess ) THEN ! no corresponding field
    WRITE(yerrmsg, "(A)") 'ERROR: no corresponding field in  blockFieldTable (2dLogical)'
    ierror = 10112
    CALL model_abort(my_cart_id, ierror, yerrmsg, yroutine)
  END IF

END SUBROUTINE registerCopy2dLogical
 
!==============================================================================

SUBROUTINE registerCopy3dReal(Field_b, CLS, CopyFlag)

!------------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------------
! Begin subroutine registerCopy3dReal
!------------------------------------------------------------------------------

  ! Subroutine arguments:
  ! --------------------
  REAL(KIND=wp), INTENT(IN), TARGET :: Field_b(:,:)
  TYPE(CopylistStruct), INTENT(INOUT), TARGET :: CLS
  INTEGER, INTENT(IN) :: copyFlag

  INTEGER, POINTER :: CL(:), ncopy
  INTEGER :: id
  LOGICAL :: lsuccess
  INTEGER :: ierror
  CHARACTER(LEN=64) :: yerrmsg
  CHARACTER(LEN=64) :: yroutine
  CHARACTER(LEN=32) :: yinfo

  yroutine='register_copy (3dReal)'

  ! Check State : allowed iStateRegField or iStateRegCopy
  CALL checkCopyBlockState(iStateRegField, iStateRegCopy, ierror, yerrmsg)
  IF (ierror > 0) CALL model_abort(my_cart_id, ierror, yerrmsg, yroutine)
  icopyBlockState = iStateRegCopy

  ! Set copy list (CL) and ncopy pointers to IN or OUT depending on the copy flag
  CALL setCopyListPointer(CLS, copyFlag, CL, ncopy, yinfo, ierror, yerrmsg)
  IF (ierror > 0) CALL model_abort(my_cart_id, ierror, yerrmsg, yroutine)

  ! Scan throught blockFieldTable
  id = 1
  lsuccess = .FALSE.
  DO WHILE( (.NOT. lsuccess) .AND. (id <= nBlockFields) )
    IF ( ASSOCIATED(blockFieldTable(id)%pfld3dR_b, Field_b) ) THEN
      ncopy = ncopy + 1
      IF (ncopy >  nBlockFields )  THEN
        WRITE(yerrmsg, "(A,I8,I8)") 'ERROR: ncopy> nBlockFields', &
             ncopy, nBlockFields
        ierror = 10133
        CALL model_abort(my_cart_id, ierror, yerrmsg, yroutine)
      END IF
      CL(ncopy) = id
      lsuccess = .TRUE.
      IF (ldebug.AND.(my_cart_id == 0)) WRITE(*,"(A,A)") TRIM(yinfo), &
           TRIM(blockFieldTable(id)%yname)
    END IF
    id = id + 1
  END DO

  IF ( .NOT. lsuccess ) THEN ! no corresponding field
    WRITE(yerrmsg, "(A)") 'ERROR: no corresponding field in blockFieldTable (3dReal)'
    ierror = 10115
    CALL model_abort(my_cart_id, ierror, yerrmsg, yroutine)
  END IF

END SUBROUTINE registerCopy3dReal

!==============================================================================
!==============================================================================
!+ Helper procedure setCopyListPointer for registerCopy in "src_block_fields"
!------------------------------------------------------------------------------

SUBROUTINE setCopyListPointer(CLS, copyFlag, CL, ncopy, yinfo, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!    Set copy list (CL) and ncopy pointers to IN or OUT component of the copy 
!    structure depending on the copy flag
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  TYPE(CopylistStruct), INTENT(IN), TARGET :: CLS
  INTEGER, INTENT(IN)             :: copyFlag
  INTEGER, POINTER, INTENT(INOUT) :: CL(:), ncopy
  CHARACTER(LEN=*), INTENT(OUT)   :: yinfo
  INTEGER, INTENT(OUT)            :: ierror
  CHARACTER(LEN=*), INTENT(OUT)   :: yerrmsg

!------------------------------------------------------------------------------
! Begin subroutine setCopyListPointer
!------------------------------------------------------------------------------
  
   IF ( copyFlag == copyToBlockF ) THEN
    CL => CLS%IN
    ncopy => CLS%ncopyin
    yinfo = 'Register copy TO block : '   ! for debugging
  ELSE IF ( copyFlag == copyFromBlockF) THEN
    CL => CLS%OUT
    ncopy => CLS%ncopyout
    yinfo = 'Register copy FROM block : ' ! for debugging
  ELSE
    WRITE(yerrmsg, "(A,I8)") 'ERROR: Invalid copyFlag ', copyFlag
    ierror = 10110
  END IF
 

END SUBROUTINE setCopyListPointer

!==============================================================================
!==============================================================================
!+ Module procedure "reset_copy_flags" in "src_block_fields" to
!  reset all copy flag in  blockFieldTable
!------------------------------------------------------------------------------

SUBROUTINE reset_copy_flags

!------------------------------------------------------------------------------
!
! Description:
!   Set all copy flag in the blockfield table to False
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------------
! Begin subroutine reset_copy_flags
!------------------------------------------------------------------------------

  blockFieldTable(:)%lcopyFlag(copyToBlockF) = .FALSE.
  blockFieldTable(:)%lcopyFlag(copyFromBlockF) = .FALSE.

END SUBROUTINE reset_copy_flags

!==============================================================================
!==============================================================================
!+ Module procedure "request_copy" in "src_block_fields" to
!  request copy
!------------------------------------------------------------------------------

SUBROUTINE request_copy(CLS, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   Set copy to block, and back to ijk flags in blockFieldTable for given
!   indices copy list structure CLS
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  TYPE(CopylistStruct), INTENT(IN) :: CLS
  INTEGER, INTENT(OUT) :: ierror
  CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg

! Local parameters:
! ----------------
  INTEGER :: i

!------------------------------------------------------------------------------
! Begin subroutine request_copy
!------------------------------------------------------------------------------

  ierror = 0
  yerrmsg = ''

  ! Check State : allowed 2 or 3
  CALL checkCopyBlockState(iStateRegCopy, iStateCopy, ierror, yerrmsg)
  IF (ierror > 0) RETURN
  icopyBlockState = iStateCopy

  ! set copytoblock flag
  CALL setCopyFlag(CLS%IN, CLS%ncopyin, copyToBlockF, ierror, yerrmsg)
  IF (ierror > 0) RETURN

  ! set copyback flag
  CALL setCopyFlag(CLS%OUT, CLS%ncopyout, copyFromBlockF, ierror, yerrmsg)
  IF (ierror > 0) RETURN

END SUBROUTINE request_copy

!==============================================================================
!==============================================================================
! Helper function for request_copy
! Set copy flag for a list of indices
!------------------------------------------------------------------------------

SUBROUTINE setCopyFlag(clist, ncopy, copyFlag, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   Set copy flags "copyFlag" to True for a list of indices
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  INTEGER, INTENT(IN) :: clist(:)
  INTEGER, INTENT(IN) :: ncopy
  INTEGER, INTENT(IN) :: copyFlag
  INTEGER, INTENT(OUT) :: ierror
  CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg

! Local parameters:
! ----------------
  INTEGER :: i, id

!------------------------------------------------------------------------------
! Begin subroutine setCopyToBlockFlag
!------------------------------------------------------------------------------

  ierror = 0
  yerrmsg = ''

  IF ( ncopy > nBlockFields ) THEN
    ierror = 10116
    yerrmsg = 'ERROR: ncopy > nBlockFields'
    RETURN
  END IF

  IF ( copyFlag > ncopyFlags ) THEN
    ierror = 10116
    yerrmsg = 'ERROR:  copyFlag > ncopyFlag'
    RETURN
  END IF

  DO i = 1, ncopy
    id = clist(i)
    IF ( id > nBlockFields ) THEN
      ierror = 10117
      yerrmsg = 'ERROR: id > nBlockFields'
      RETURN
    END IF
    blockFieldTable(id)%lcopyFlag(copyFlag) = .TRUE.
  END DO

END SUBROUTINE setCopyFlag

!==============================================================================
!==============================================================================
!+ Module procedures "copy_to_block" in "src_block_fields" to apply
!  requested copies
!------------------------------------------------------------------------------

SUBROUTINE copy_to_block(CLS, ipend, ib, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   Apply copy to block for given list structure CLS
!   Copy are applied for each type and ksize
!   Active copies, i.e. with copyflag still true, are selected
!   using the getActiveCopyList routine
!   This routine shoud be called after a call to requestcopy
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  TYPE(CopylistStruct), INTENT(INOUT) :: CLS
  INTEGER, INTENT(IN) :: ipend, ib  ! block end, current block
  INTEGER, INTENT(OUT) :: ierror
  CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg

! Locals :
! ----------------
  INTEGER :: icpt, iks      ! loop indices
  INTEGER :: ksize
  INTEGER :: CopyList(nBlockFields) ! list of active copy
  INTEGER :: ncopy                  ! number of active copy

!----------- End of header ----------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine  copy_to_block
!------------------------------------------------------------------------------

  ierror = 0
  yerrmsg = ''

  ! Check Status
  CALL checkCopyBlockState(iStateCopy, iStateCopy, ierror, yerrmsg)
  IF (ierror > 0) RETURN


  ! Loop over all copy type
  DO icpt = 1, ncopyType
    ! Loop over different ksize
    DO iks = 1, copyTypeNKSize(icpt)

      ksize = copyTypeKSize(icpt, iks)

      ! Get list of requested active copy of this type and ksize
      CALL getActiveCopyList(CLS%IN, CLS%ncopyin, copyToBlockF, icpt, ksize, &
           CopyList, ncopy, ierror, yerrmsg)
      IF (ierror > 0) RETURN

      ! Apply
      CALL doCopy(Copylist, ncopy, copyToBlockF, icpt, ksize, ipend, ib, &
           ierror, yerrmsg)
      IF (ierror > 0) RETURN

      ! Set copyToblock flags to .FALSE.
      CALL resetCopyFlag(Copylist, ncopy, copyToBlockF, ierror, yerrmsg)
      IF (ierror > 0) RETURN

    END DO
  END DO

!------------------------------------------------------------------------------
! End of module procedure copy_to_block
!------------------------------------------------------------------------------

END SUBROUTINE copy_to_block

!==============================================================================
!==============================================================================
!+ Module procedures "copy_from_block" in "src_block_fields" to apply
!  requested copies
!------------------------------------------------------------------------------

SUBROUTINE copy_from_block(CLS, ipend, ib, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   Apply copy from block for given list structure CLS
!   Copy are applied for each type and ksize
!   Active copies, i.e. with copyflag still true, are selected
!   using the getActiveCopyList routine
!   This routine shoud be called after a call to requestcopy
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  TYPE(CopylistStruct), INTENT(INOUT) :: CLS
  INTEGER, INTENT(IN) :: ipend, ib  ! block end, current block
  INTEGER, INTENT(OUT) :: ierror
  CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg

! Locals :
! ----------------
  INTEGER :: icpt, iks      ! loop indices
  INTEGER :: ksize
  INTEGER :: CopyList(nBlockFields) ! active copy list
  INTEGER :: ncopy                  ! number of active copy

!----------- End of header ----------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine  copy_from_block
!------------------------------------------------------------------------------

  ierror = 0
  yerrmsg = ''

  ! Check Status
  CALL checkCopyBlockState(iStateCopy, iStateCopy, ierror, yerrmsg)
  IF (ierror > 0) RETURN

  ! Loop over all copy type
  DO icpt = 1, ncopyType
    ! Loop over different ksize
    DO iks = 1, copyTypeNKSize(icpt)

      ksize = copyTypeKSize(icpt, iks)

      ! Get list of requested active copy of this type and ksize
      CALL getActiveCopyList(CLS%OUT, CLS%ncopyout, copyFromBlockF, icpt, ksize, &
           CopyList, ncopy, ierror, yerrmsg)
      IF (ierror > 0) RETURN

      ! Apply
      CALL doCopy(Copylist, ncopy, copyFromBlockF, icpt, ksize, ipend, ib, &
           ierror, yerrmsg)
      IF (ierror > 0) RETURN

      ! Set copyToblock flags to .FALSE.
      CALL resetCopyFlag(Copylist, ncopy, copyFromBlockF, ierror, yerrmsg)
      IF (ierror > 0) RETURN

    END DO
  END DO

!------------------------------------------------------------------------------
! End of module procedure copy_from_block
!------------------------------------------------------------------------------

END SUBROUTINE copy_from_block
!==============================================================================


!------------------------------------------------------------------------------
! Helper function for copy_to_block
! Set copy flag to .FALSE. for a list of indices
!------------------------------------------------------------------------------

SUBROUTINE resetCopyFlag(clist, ncopy, copyFlag, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   Set copy flags "copyFlag" to False for a list of indices
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  INTEGER, INTENT(IN) :: clist(:)
  INTEGER, INTENT(IN) :: ncopy
  INTEGER, INTENT(IN) :: copyFlag
  INTEGER, INTENT(OUT) :: ierror
  CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg

! Local parameters:
! ----------------
  INTEGER :: i, id

!------------------------------------------------------------------------------
! Begin subroutine setCopyToBlockFlag
!------------------------------------------------------------------------------

  ierror = 0
  yerrmsg = ''

  ! This routine is called very offen.
  IF ( ncopy > nBlockFields ) THEN
    ierror = 10122
    yerrmsg = 'ERROR: ncopy > nBlockFields'
    RETURN
  END IF

  IF ( copyFlag > ncopyFlags ) THEN
    ierror = 10122
    yerrmsg = 'ERROR:  copyFlag > ncopyFlag'
    RETURN
  END IF

  DO i = 1, ncopy
    id = clist(i)
    IF ( id > nBlockFields ) THEN
      ierror = 10117
      yerrmsg = 'ERROR: id > nBlockFields'
      RETURN
    END IF
    blockFieldTable(id)%lcopyFlag(copyFlag) = .FALSE.
  END DO

END SUBROUTINE resetCopyFlag

!==============================================================================
!==============================================================================
!+ Module procedures "doCopy" in "src_block_fields" to apply
!  all copies for a given copy list
!------------------------------------------------------------------------------

SUBROUTINE doCopy(CopyList, ncopy, copyFlag, copytype, ksize, ipend, ib, &
                  ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   Apply copy given list.
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  INTEGER, INTENT(IN) :: CopyList(:)         !list of field to be copeid
  INTEGER, INTENT(IN) :: ncopy, copyFlag     !copy size and Flag (to or from)
  INTEGER, INTENT(IN) :: copytype, ksize     !copy type
  INTEGER, INTENT(IN) :: ipend, ib ! block end, current block
  INTEGER, INTENT(OUT):: ierror
  CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg

! Locals :
! ----------------
  INTEGER :: i, id, ntime
! Helper pointers to fields
  REAL(KIND=wp), POINTER :: pfld3dR(:,:,:), pfld3dR_b(:,:)
  REAL(KIND=wp), POINTER :: pfld2dR(:,:), pfld2dR_b(:)
  INTEGER, POINTER       :: pfld2dI(:,:), pfld2dI_b(:)
  LOGICAL, POINTER       :: pfld2dL(:,:), pfld2dL_b(:)
  CHARACTER(LEN=32) :: copyInfo

!----------- End of header ----------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine  doCopy
!------------------------------------------------------------------------------

  ierror = 0
  yerrmsg = ''

  ! Nullify helper pointer
  NULLIFY(pfld3dR,pfld3dR_b,pfld2dR,pfld2dR_b,pfld2dI,pfld2dI_b,pfld2dL,pfld2dL_b)

  ! Set print info
  IF ( ldebug ) THEN
    IF (copyFlag == copyToBlockF) THEN
      copyInfo = 'Copy to block :'
    ELSE
      copyInfo = 'Copy from block :'
    END IF
  END IF

  ! Check Status
  CALL checkCopyBlockState(iStateCopy, iStateCopy, ierror, yerrmsg)
  IF (ierror > 0) RETURN

  DO i = 1, ncopy
    id = CopyList(i)
    IF ( ldebug .AND. (ib == 1).AND.(my_cart_id == 0) )  &
      PRINT *, TRIM(copyInfo), TRIM(blockFieldTable(id)%yname)

    ! Handle all possible copy types
    SELECT CASE (copytype)

    CASE (icopyTypeReal2D) ! 2d real fields
      IF ( blockFieldTable(id)%timelevel ) THEN
        ntime = blockFieldTable(id)%ntime
        pfld2dR => blockFieldTable(id)%pfld2dRtime(:,:,ntime)
      ELSE
        pfld2dR => blockFieldTable(id)%pfld2dR
      END IF
      pfld2dR_b => blockFieldTable(id)%pfld2dR_b

      ! Check pointer association
      IF (.NOT.( ASSOCIATED(pfld2dR) .AND. ASSOCIATED(pfld2dR_b))) THEN
        ierror = 10120
        yerrmsg = 'ERROR: field2d/_b not associated'
        RETURN
      END IF

      IF (copyFlag == copyToBlockF) THEN
        CALL copyToBlock2d(pfld2dR, pfld2dR_b, ipend, ib)
      ELSE IF (copyFlag == copyFromBlockF) THEN
        CALL copyFromBlock2d(pfld2dR_b, pfld2dR, ipend, ib)
      ELSE 
        ierror = 10123
        yerrmsg = 'ERROR: invalid copyFlag'
        RETURN
      END IF

      NULLIFY(pfld2dR,pfld2dR_b)

    CASE (icopyTypeInt2D) ! 2d integer fields

      pfld2dI => blockFieldTable(id)%pfld2dI
      pfld2dI_b => blockFieldTable(id)%pfld2dI_b

      ! Check pointer association
      IF (.NOT.( ASSOCIATED(pfld2dI) .AND. ASSOCIATED(pfld2dI_b))) THEN
        ierror = 10126
        yerrmsg = 'ERROR: field2dI/_b not associated'
        RETURN
      END IF

      IF (copyFlag == copyToBlockF) THEN
        CALL copyToBlock2dI(pfld2dI, pfld2dI_b, ipend, ib)
      ELSE IF (copyFlag == copyFromBlockF) THEN
        CALL copyFromBlock2dI(pfld2dI_b, pfld2dI, ipend, ib)
      ELSE 
        ierror = 10127
        yerrmsg = 'ERROR: invalid copyFlag'
        RETURN
      END IF

      NULLIFY(pfld2dI,pfld2dI_b)

    CASE (icopyTypeLogical2D) ! 2d logical fields

      pfld2dL => blockFieldTable(id)%pfld2dL
      pfld2dL_b => blockFieldTable(id)%pfld2dL_b

      ! Check pointer association
      IF (.NOT.( ASSOCIATED(pfld2dL) .AND. ASSOCIATED(pfld2dL_b))) THEN
        ierror = 10128
        yerrmsg = 'ERROR: field2dL/_b not associated'
        RETURN
      END IF

      IF (copyFlag == copyToBlockF) THEN
        CALL copyToBlock2dL(pfld2dL, pfld2dL_b, ipend, ib)
      ELSE IF (copyFlag == copyFromBlockF) THEN
        CALL copyFromBlock2dL(pfld2dL_b, pfld2dL, ipend, ib)
      ELSE 
        ierror = 10129
        yerrmsg = 'ERROR: invalid copyFlag'
        RETURN
      END IF

      NULLIFY(pfld2dL,pfld2dL_b)

    CASE (icopyTypeReal3D) ! 3d real fields

      ! retrieve correct point to source field
      IF ( blockFieldTable(id)%timelevel ) THEN
        ntime = blockFieldTable(id)%ntime
        pfld3dR => blockFieldTable(id)%pfld3dRtime(:,:,:,ntime)
      ELSE
        pfld3dR => blockFieldTable(id)%pfld3dR
      END IF
      pfld3dR_b => blockFieldTable(id)%pfld3dR_b

      ! Check pointer association
      IF (.NOT.( ASSOCIATED(pfld3dR) .AND. ASSOCIATED(pfld3dR_b))) THEN
        ierror = 10124
        yerrmsg = 'ERROR: field3d/_b not associated'
        RETURN
      END IF

      IF (copyFlag == copyToBlockF) THEN
        CALL copyToBlock3d(pfld3dR, pfld3dR_b, ipend, ksize, ib)
      ELSE IF (copyFlag == copyFromBlockF) THEN
        CALL copyFromBlock3d(pfld3dR_b, pfld3dR, ipend, ksize, ib)
      ELSE
        ierror = 10125
        yerrmsg = 'ERROR: invalid copyFlag'
      END IF
      
      NULLIFY(pfld3dR,pfld3dR_b)

    END SELECT
  END DO

!------------------------------------------------------------------------------
! End of module procedure doCopy
!------------------------------------------------------------------------------

END SUBROUTINE doCopy

!==============================================================================
!==============================================================================
!+ Module procedures "copyToBlockX" in "src_block_fields"
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!
! Description:
!   Actual implementation of copy to block for 2d and 3d fields
!
!------------------------------------------------------------------------------

!==============================================================================
! copytoblock  : 3d fields
!-----------------------------------------------------------------------------

SUBROUTINE copyToBlock3d(fld, fld_b, ipend, ksize, ib)

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  REAL(KIND=wp), INTENT(IN) :: fld(:,:,:)
  REAL(KIND=wp), INTENT(OUT) :: fld_b(:,:)
  INTEGER, INTENT(IN) :: ipend, ksize, ib

! Local parameters:
! ----------------
  INTEGER :: ip, i, j, k

!------------------------------------------------------------------------------
! Begin Subroutine  copyToBlock3d
!------------------------------------------------------------------------------

  !NOacc data present(fld, fld_b, mind_ilon, mind_jlat)
  !NOacc parallel
  DO k = 1, ksize
    !NOacc loop gang vector
    DO ip = 1, ipend
      i = mind_ilon(ip,ib)
      j = mind_jlat(ip,ib)
      fld_b(ip,k) = fld(i,j,k)
    END DO
  END DO
  !NOacc end parallel
  !NOacc end data
  
END SUBROUTINE copyToBlock3d

!==============================================================================
!==============================================================================
! copytoblock  : 2d fields
!-----------------------------------------------------------------------------

SUBROUTINE copyToBlock2d(fld, fld_b, ipend, ib)

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  REAL(KIND=wp), INTENT(IN) :: fld(:,:)
  REAL(KIND=wp), INTENT(OUT) :: fld_b(:)
  INTEGER, INTENT(IN) :: ipend, ib

! Local parameters:
! ----------------
  INTEGER :: ip, i, j

!------------------------------------------------------------------------------
! Begin Subroutine  copyToBlock2d
!------------------------------------------------------------------------------

  !NOacc data present(fld, fld_b, mind_ilon, mind_jlat)
  !NOacc parallel
  !NOacc loop gang vector
  DO ip = 1, ipend
    i = mind_ilon(ip,ib)
    j = mind_jlat(ip,ib)
    fld_b(ip) = fld(i,j)
  END DO
  !NOacc end parallel
  !NOacc end data

END SUBROUTINE copyToBlock2d

!==============================================================================
! copytoblock  : 2d fields (integer)
!-----------------------------------------------------------------------------

SUBROUTINE copyToBlock2dI(fld, fld_b, ipend, ib)

! Subroutine arguments:
! --------------------
  INTEGER, INTENT(IN)  :: fld(:,:)
  INTEGER, INTENT(OUT) :: fld_b(:)
  INTEGER, INTENT(IN) :: ipend, ib

! Local parameters:
! ---------------- 
  INTEGER :: ip, i, j

!------------------------------------------------------------------------------
! Begin Subroutine  copyToBlock2dI
!------------------------------------------------------------------------------

  !NOacc data present(fld, fld_b, mind_ilon, mind_jlat)
  !NOacc parallel
  !NOacc loop gang vector
  DO ip = 1, ipend
     i = mind_ilon(ip,ib)
     j = mind_jlat(ip,ib)
     fld_b(ip) = fld(i,j)
  END DO
  !NOacc end parallel
  !NOacc end data

END SUBROUTINE copyToBlock2dI

!==============================================================================
!==============================================================================
! copytoblock  : 2d fields (logical)
!-----------------------------------------------------------------------------

SUBROUTINE copyToBlock2dL(fld, fld_b, ipend, ib)

! Subroutine arguments:
! --------------------
  LOGICAL, INTENT(IN)  :: fld(:,:)
  LOGICAL, INTENT(OUT) :: fld_b(:)
  INTEGER, INTENT(IN) :: ipend, ib

! Local parameters:
! ---------------- 
  INTEGER :: ip, i, j

!------------------------------------------------------------------------------
! Begin Subroutine  copyToBlock2dL
!------------------------------------------------------------------------------

  !NOacc data present(fld, fld_b, mind_ilon, mind_jlat)
  !NOacc parallel
  !NOacc loop gang vector
  DO ip = 1, ipend
     i = mind_ilon(ip,ib)
     j = mind_jlat(ip,ib)
     fld_b(ip) = fld(i,j)
  END DO
  !NOacc end parallel
  !NOacc end data

END SUBROUTINE copyToBlock2dL


!==============================================================================
!==============================================================================
!+ Module procedures "copyFromBlockX" in "src_block_fields" 
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!
! Description:
!   Actual implementation of copy from block for 2d and 3d fields
!
!------------------------------------------------------------------------------

!==============================================================================
! copyFromBlock3d  : 3d fields
!------------------------------------------------------------------------------

SUBROUTINE copyFromBlock3d(fld_b, fld, ipend, ksize, ib)

!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  REAL(KIND=wp), INTENT(IN) :: fld_b(:,:)                ! src : ijk
  REAL(KIND=wp), INTENT(INOUT) :: fld(:,:,:)                ! intent inout
  INTEGER, INTENT(IN) :: ipend, ksize, ib

! Local parameters:
! ----------------
  INTEGER :: ip, i, j, k

!------------------------------------------------------------------------------
! Begin Subroutine  copyFromBlock3d
!------------------------------------------------------------------------------

  !NOacc data present(fld_b, fld, mind_ilon, mind_jlat)
  !NOacc parallel
  DO k = 1, ksize
     !NOacc loop gang vector
    DO ip = 1, ipend
      i = mind_ilon(ip,ib)
      j = mind_jlat(ip,ib)
      fld(i,j,k) = fld_b(ip,k)
    END DO
  END DO
  !NOacc end parallel
  !NOacc end data

END SUBROUTINE copyFromBlock3d

!==============================================================================
! copytoblock  : 2d fields
!------------------------------------------------------------------------------

SUBROUTINE copyFromBlock2d(fld_b, fld, ipend, ib)

!-----------------------------------------------------------------------------

  IMPLICIT NONE


! Subroutine arguments:
! --------------------
  REAL(KIND=wp), INTENT(IN) :: fld_b(:)
  REAL(KIND=wp), INTENT(INOUT) :: fld(:,:)
  INTEGER, INTENT(IN) :: ipend, ib

! Local parameters:
! ----------------
  INTEGER :: ip, i, j

!------------------------------------------------------------------------------
! Begin Subroutine  copyFromBlock2d
!------------------------------------------------------------------------------

  !NOacc data present(fld_b, fld, mind_ilon, mind_jlat)
  !NOacc parallel
  !NOacc loop gang vector
  DO ip = 1, ipend
    i = mind_ilon(ip,ib)
    j = mind_jlat(ip,ib)
    fld(i,j) = fld_b(ip)
  END DO
  !NOacc end parallel
  !NOacc end data

END SUBROUTINE copyFromBlock2d

!==============================================================================
! copytoblock  : 2d fields (integer)
!------------------------------------------------------------------------------

SUBROUTINE copyFromBlock2dI(fld_b, fld, ipend, ib)

!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------
  INTEGER, INTENT(IN)    :: fld_b(:)
  INTEGER, INTENT(INOUT) :: fld(:,:)
  INTEGER, INTENT(IN) :: ipend, ib

! Local parameters:
! ---------------- 
  INTEGER :: ip, i, j

!------------------------------------------------------------------------------
! Begin Subroutine  copyFromBlock2dI
!------------------------------------------------------------------------------

  !NOacc data present(fld_b, fld, mind_ilon, mind_jlat)
  !NOacc parallel
  !NOacc loop gang vector
  DO ip = 1, ipend
     i = mind_ilon(ip,ib)
     j = mind_jlat(ip,ib)
     fld(i,j) = fld_b(ip)
  END DO
  !NOacc end parallel
  !NOacc end data

END SUBROUTINE copyFromBlock2dI

!==============================================================================
!==============================================================================
! copytoblock  : 2d fields (logical)
!------------------------------------------------------------------------------

SUBROUTINE copyFromBlock2dL(fld_b, fld, ipend, ib)

!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------
  LOGICAL, INTENT(IN)    :: fld_b(:)
  LOGICAL, INTENT(INOUT) :: fld(:,:)
  INTEGER, INTENT(IN) :: ipend, ib

! Local parameters:
! ---------------- 
  INTEGER :: ip, i, j

!------------------------------------------------------------------------------
! Begin Subroutine  copyFromBlock2dL
!------------------------------------------------------------------------------

  !NOacc data present(fld_b, fld, mind_ilon, mind_jlat)
  !NOacc parallel
  !NOacc loop gang vector
  DO ip = 1, ipend
     i = mind_ilon(ip,ib)
     j = mind_jlat(ip,ib)
     fld(i,j) = fld_b(ip)
  END DO
  !NOacc end parallel
  !NOacc end data

END SUBROUTINE copyFromBlock2dL

!==============================================================================
!==============================================================================
! Helper function for copy to/from block
! Returns a list all active copy in a given copy list 
!------------------------------------------------------------------------------

SUBROUTINE getActiveCopyList(clist, ncopy, copyFlag, copyType, ksize, &
     activeClist, nactiveCopy, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   Returns a list of all active copy in a given copy list (i.e. with copy flag=T)
!   with same copy type and ksize
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  INTEGER, INTENT(IN) :: clist(:)
  INTEGER, INTENT(IN) :: ncopy
  INTEGER, INTENT(IN) :: copyFlag, copyType, ksize
  INTEGER, INTENT(OUT) :: activeClist(:)
  INTEGER, INTENT(OUT) :: nactiveCopy
  INTEGER, INTENT(OUT) :: ierror
  CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg

! Local parameters:
! ----------------
  INTEGER :: i, id
  LOGICAL, POINTER :: lcopyFlag

!------------------------------------------------------------------------------
! Begin subroutine setCopyFromBlockFlag
!------------------------------------------------------------------------------

  ierror = 0
  yerrmsg = ''

  nactiveCopy = 0

  IF ( ncopy > nBlockFields ) THEN
    ierror = 10118
    yerrmsg = 'ERROR: ncopy > nBlockFields'
    RETURN
  END IF

  DO i = 1, ncopy
    id = clist(i)
    IF ( blockFieldTable(id)%lcopyFlag(copyFlag) &
         .AND. blockFieldTable(id)%copyType == copyType &
         .AND. blockFieldTable(id)%ksize == ksize ) THEN
      nactiveCopy = nactiveCopy + 1
      activeClist(nactiveCopy) = id
    END IF
  END DO

END SUBROUTINE getActiveCopyList

!==============================================================================
!==============================================================================
!+ Module procedures "finalize_copy" in "src_block_fields"
!------------------------------------------------------------------------------

SUBROUTINE finalize_copy(ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   Checks that all requested copied have beed done, i.e. all lcopyToBlock
!   and lcopyFromBlock are FALSE
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  INTEGER, INTENT(OUT) :: ierror
  CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg

! Local parameters:
! ----------------
  INTEGER :: i

!------------------------------------------------------------------------------
! Begin Subroutine  finalize_copy
!------------------------------------------------------------------------------

  ierror = 0
  yerrmsg = ''

  DO i = 1, blockFieldTableSize
    ! One of this error means that there is a missing copy_to_block or
    ! copy_from_block
    IF (blockFieldTable(i)%lcopyFlag(copyToBlockF)) THEN
      ierror = 10130
      WRITE(yerrmsg, "(A,A)") 'ERROR: Field was not copied to block : ',  &
            TRIM(blockFieldTable(i)%yname)
      RETURN
    END IF
    IF (blockFieldTable(i)%lcopyFlag(copyFromBlockF)) THEN
      ierror = 10131
      WRITE(yerrmsg, "(A,A)") 'ERROR: Field was not copied from block : ',  &
            TRIM(blockFieldTable(i)%yname)
      RETURN
    END IF
  END DO

END SUBROUTINE finalize_copy

!==============================================================================
!==============================================================================
! Helper function for src_block_fields 
!------------------------------------------------------------------------------

SUBROUTINE checkCopyBlockState(ValidMin, ValidMax, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!    Check that current state machine is between ValidMin and Max
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  INTEGER, INTENT(IN) :: ValidMin, ValidMax
  INTEGER, INTENT(OUT) :: ierror
  CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg

!------------------------------------------------------------------------------
! Begin Subroutine  checkCopyBlockState
!------------------------------------------------------------------------------

  ierror = 0
  yerrmsg = ''

  IF ( .NOT.( (icopyBlockState >= ValidMin) .OR. (icopyBlockState <= ValidMax) ) ) THEN
    ierror = 10107
    WRITE(yerrmsg, "(A,I8)") 'ERROR: Not valid state :', icopyBlockState
    RETURN
  END IF

END SUBROUTINE checkCopyBlockState

!==============================================================================
!==============================================================================
! Helper Function for init and finalize
!------------------------------------------------------------------------------

SUBROUTINE clearBlockFieldTable

!------------------------------------------------------------------------------
!
! Description:
!  Clear BlockFieldTable fields, nulify pointers
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

! Local parameters:
! ----------------
  INTEGER :: icount

!------------------------------------------------------------------------------
! Begin Subroutine  clearBlockFieldTable
!------------------------------------------------------------------------------

  ! set all flag to unused and pointer to fields
  blockFieldTable(:)%lused        = .FALSE.
  blockFieldTable(:)%yname        = ''
  blockFieldTable(:)%rank         = -1
  blockFieldTable(:)%timelevel    = .FALSE.
  blockFieldTable(:)%ksize        = -1
  DO icount = 1, blockFieldTableSize
    NULLIFY( blockFieldTable(icount)%ntime )
    NULLIFY( blockFieldTable(icount)%pfld3dRtime )
    NULLIFY( blockFieldTable(icount)%pfld3dR )
    NULLIFY( blockFieldTable(icount)%pfld2dRtime )
    NULLIFY( blockFieldTable(icount)%pfld2dR )
    NULLIFY( blockFieldTable(icount)%pfld2dI )
    NULLIFY( blockFieldTable(icount)%pfld2dL )
    NULLIFY( blockFieldTable(icount)%pfld3dR_b )
    NULLIFY( blockFieldTable(icount)%pfld2dR_b )
    NULLIFY( blockFieldTable(icount)%pfld2dI_b )
    NULLIFY( blockFieldTable(icount)%pfld2dL_b )
  END DO
  blockFieldTable(:)%lcopyFlag(copyToBlockF) = .FALSE.
  blockFieldTable(:)%lcopyFlag(copyFromBlockF) = .FALSE.
  blockFieldTable(:)%copyType = -1

  ! Reset number of block fields
  nBlockFields = 0

END SUBROUTINE clearBlockFieldTable

!==============================================================================

END MODULE src_block_fields
