!+ Data structure and constants to put/get arbitrary metadata to/from a buffer
!------------------------------------------------------------------------------

MODULE data_tracer_metadata

!------------------------------------------------------------------------------
! Description:
!   The module "data_tracer_metadata" contains the global declarations for the 
!   module "src_tracer_metadata" which allows the user to store arbitrary (i.e.
!   data of type integer, real, and string as scalar
!   or vector) in a buffer. This is convenient to provide metadata
!   to an field or I/O where the requirements of the user are not known
!   exactly in advance or might change rapidly.
!------------------------------------------------------------------------------
!
! Current Code Owner: MeteoSwiss, Oliver Fuhrer
!  phone: +41 58 460 9359
!  fax:   +41 58 460 9278
!  email: oliver.fuhrer@meteoswiss.ch
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer
!  Initial release
! V4_26        2012/12/06 Anne Roches
!  Addition of the pointer support in the metadata. This function is required
!  for handling associated fields (surface field, emissions, ...) gracefully.
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
! V5_4         2016-03-10 Oliver Fuhrer
!  Updated code owner information
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Declarations:
!
! Modules used:
!
USE iso_c_binding, ONLY : c_ptr

USE data_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================

!------------------------------------------------------------------------------
! Local Declarations:

! Local Parameters:

INTEGER (KIND=iintegers), PARAMETER ::    &
  I_MAX_KEYLEN          =   32 , & ! maximum string length for key
  I_MAX_STRLEN          =   80 , & ! maximum string length for storage
  I_STD_NUMBUF          =  128 , & ! maximum number of buffers
  I_STD_NUMKEY          =  128 , & ! default number of keys
  I_STD_BUFLEN          = 4096     ! default storage size in Bytes
INTEGER (KIND=iintegers), PARAMETER :: I_PTR_MAXDIM = 7    ! maximum number of dimensions for pointer storage

! Definition for data types
INTEGER (KIND=iintegers), PARAMETER ::    &
  I_DATA_TYPE_UNDEF     =    0 , & !
  I_DATA_TYPE_INTEGER   =    1 , & !
  I_DATA_TYPE_REAL      =    2 , & !
  I_DATA_TYPE_DOUBLE    =    3 , & !
  I_DATA_TYPE_STRING    =    4 , & !
  I_DATA_TYPE_LOGICAL   =    5 , & !
  I_DATA_TYPE_POINTER   =    6 , & !
  I_NUM_DATA_TYPES      =    6     !

INTEGER (KIND=iintegers), PARAMETER ::    &
  I_STATUS_READ         =    1 , & ! can be read
  I_STATUS_WRITE        =    2 , & ! can be written
  I_STATUS_DELETE       =    4 , & ! can be deleted

  I_STATUS_DEFAULT      =        & ! default (undef, read, write, delete)
       I_STATUS_READ    +        &
       I_STATUS_WRITE   +        &
       I_STATUS_DELETE         , & !

  I_STATUS_MASK         =        & ! mask for valid attribute bits
       I_STATUS_READ    +        &
       I_STATUS_WRITE   +        &
       I_STATUS_DELETE             !


! special buffer number to write to all buffers
INTEGER (KIND=iintegers), PARAMETER ::    &
  I_ALL_BUF             = -98765, & !
  I_GETBYNAME           = -98765, & !
  I_PUTBYNAME           = -98765

! data type short string (for pretty print)
CHARACTER(LEN=12), PARAMETER :: &  
   Y_DATA_TYPE(0:I_NUM_DATA_TYPES) =                  &
      (/'UNDEF       ','INTEGER     ','REAL        ', &
        'DOUBLE      ','STRING      ','LOGICAL     ', &
        'POINTER     '/)


! Local Type Definitions:

TYPE t_fld_pointer
  TYPE(c_ptr) :: fld_cptr
  INTEGER (KIND=iintegers) :: fld_rank
  INTEGER (KIND=iintegers) :: fld_shape(I_PTR_MAXDIM)
END TYPE t_fld_pointer

TYPE t_metadata
  LOGICAL                                        :: lisready = .FALSE. ! ready for use?
  INTEGER (KIND=iintegers)                       :: imaxkeys           ! maximum number of elements
  INTEGER (KIND=iintegers)                       :: inumkey            ! number of elements currently stored
  INTEGER (KIND=iintegers)                       :: iunique            ! counter for unique key indexing
  INTEGER (KIND=iintegers), ALLOCATABLE          :: iidx(:)            ! unique index of element
  CHARACTER(LEN=I_MAX_KEYLEN), ALLOCATABLE       :: ykey(:)            ! name of element
  INTEGER (KIND=iintegers), ALLOCATABLE          :: iattr(:)           ! attribute bits
  INTEGER (KIND=iintegers), ALLOCATABLE          :: itype(:)           ! type of element (of I_DATA_TYPE_*)
  INTEGER (KIND=iintegers), ALLOCATABLE          :: ipos(:)            ! position of elements in ybuf/ydefault array
  INTEGER (KIND=iintegers), ALLOCATABLE          :: ilen(:)            ! length of elements in ybuf/ydefault array
  INTEGER (KIND=iintegers)                       :: imaxbuf            ! maximum number of buffers
  INTEGER (KIND=iintegers)                       :: inumbuf            ! number of buffers in use
  INTEGER (KIND=iintegers)                       :: imaxbuflen         ! size of ybuf/ydefault
  CHARACTER(LEN=1), ALLOCATABLE                  :: ydefault(:)        ! default value storage
  CHARACTER(LEN=1), ALLOCATABLE                  :: ybuf(:,:)          ! data storage
END TYPE t_metadata

! Local Arrays:

INTEGER (KIND=iintegers) :: I_DATA_TYPE_SIZE(I_NUM_DATA_TYPES)         ! data type size in bytes

!------------------------------------------------------------------------------

END MODULE data_tracer_metadata
