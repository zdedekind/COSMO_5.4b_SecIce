!+ Data structure and methods to set/get arbitrary metadata to/from a buffer
!------------------------------------------------------------------------------

MODULE src_tracer_metadata

!------------------------------------------------------------------------------
! Description:
!   The module "src_tracer_metadata" allows the user to store arbitrary (i.e.
!   data of type integer, real, and string as scalar
!   or vector) in a buffer. This is convenient to provide metadata
!   to an field or I/O where the requirements of the user are not known
!   exactly in advance or might change rapidly.
!
!   The data is typecasted into a character stream and stored in a
!   character buffer of predefined dimensions. The individual data
!   elements can be retrieved either by name or by index in the meta-
!   data buffer.
!------------------------------------------------------------------------------
!
! Current Code Owner: MeteoSwiss, Oliver Fuhrer
!  phone:  +41 58 460 9359
!  fax:    +41 58 460 9278
!  email:  oliver.fuhrer@meteoswiss.ch
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V4_25        2012/09/28 Oliver Fuhrer
!  Initial release
! V4_26        2012/12/06 Anne Roches
!  Addition of the pointer support in the metadata. This function is required
!   for handling associated fields (surface field, emissions, ...) gracefully.
!  Cleanup and cosmetics (removal of unnecessary calls, correction in
!   comments, update of headers, ...)
! V4_27        2013/03/19 Ulrich Blahak
!  Corrected syntax (complaints by Intel compiler)
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
USE iso_c_binding, ONLY : C_PTR, C_NULL_PTR, C_LOC, C_F_POINTER, C_ASSOCIATED

USE data_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    sp,        & ! KIND-type parameter for real variables (single precision)
    dp,        & ! KIND-type parameter for real variables (double precision)
    iintegers    ! KIND-type parameter for standard integer variables

USE data_tracer_metadata

!------------------------------------------------------------------------------

IMPLICIT NONE

!------------------------------------------------------------------------------

PRIVATE

!==============================================================================

!------------------------------------------------------------------------------
! Global (i.e. public) Declarations:

PUBLIC :: metadata_init, metadata_modify, metadata_finish, metadata_print, &
          metadata_define, metadata_set, metadata_get, metadata_delete, &
          metadata_error

!------------------------------------------------------------------------------
! Local (i.e. private) Declarations:

! Local Scalars:

LOGICAL :: lfirst_call = .TRUE.      ! used to initialize globals

INTEGER :: imodule_errstat           ! error return value for module
CHARACTER(LEN=256) :: ymodule_errstr ! error description

! Operator definitions:

INTERFACE metadata_define
  MODULE PROCEDURE meta_define_real4   , meta_define_real8,                   &
                   meta_define_integer , meta_define_string,                  &
                   meta_define_logical , meta_define_pointer2,                &
                   meta_define_pointer3, meta_define_pointer4
END INTERFACE

INTERFACE metadata_delete
  MODULE PROCEDURE meta_delete_byname, meta_delete_byindex
END INTERFACE

INTERFACE metadata_set
  MODULE PROCEDURE meta_set_real4_byname       , meta_set_real8_byname,       &
                   meta_set_integer_byname     , meta_set_string_byname,      &
                   meta_set_logical_byname     , meta_set_pointer2_byname,    &
                   meta_set_pointer3_byname    , meta_set_pointer4_byname,    &
                   meta_set_real4_arr_byname   , meta_set_real8_arr_byname,   &
                   meta_set_integer_arr_byname , meta_set_string_arr_byname,  &
                   meta_set_logical_arr_byname ,                              &
                   meta_set_real4_byindex      , meta_set_real8_byindex,      &
                   meta_set_integer_byindex    , meta_set_string_byindex,     &
                   meta_set_logical_byindex    , meta_set_pointer2_byindex,   &
                   meta_set_pointer3_byindex   , meta_set_pointer4_byindex,   &
                   meta_set_real4_arr_byindex  , meta_set_real8_arr_byindex,  &
                   meta_set_integer_arr_byindex, meta_set_string_arr_byindex, &
                   meta_set_logical_arr_byindex
END INTERFACE

INTERFACE metadata_get
  MODULE PROCEDURE meta_get_real4_byname       , meta_get_real8_byname,       &
                   meta_get_integer_byname     , meta_get_string_byname,      &
                   meta_get_logical_byname     , meta_get_pointer2_byname,    &
                   meta_get_pointer3_byname    , meta_get_pointer4_byname,    &
                   meta_get_real4_arr_byname   , meta_get_real8_arr_byname,   &
                   meta_get_integer_arr_byname , meta_get_string_arr_byname,  & 
                   meta_get_logical_arr_byname ,                              &
                   meta_get_real4_byindex      , meta_get_real8_byindex,      &
                   meta_get_integer_byindex    , meta_get_string_byindex,     &
                   meta_get_logical_byindex    , meta_get_pointer2_byindex,   &
                   meta_get_pointer3_byindex   , meta_get_pointer4_byindex,   &
                   meta_get_real4_arr_byindex  , meta_get_real8_arr_byindex,  &
                   meta_get_integer_arr_byindex, meta_get_string_arr_byindex, &
                   meta_get_logical_arr_byindex
END INTERFACE

!------------------------------------------------------------------------------

!==============================================================================
! Module procedures in "src_tracer_metadata"
!==============================================================================

CONTAINS

!==============================================================================
!+ initialize a new metadata storage container
!------------------------------------------------------------------------------

FUNCTION metadata_error(yerrmsg)
!------------------------------------------------------------------------------
!
! Description:
!   This public function returns the error status of the metadata module and a
!   string describing the nature of the error. All internal error handling is
!   done via setting the internal module error status and the caller of any
!   module procedure can query the success/failure of his call by using this
!   function.
!
! Method:
!
!==============================================================================

IMPLICIT NONE

! Argument list
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: yerrmsg

! Return value
INTEGER :: metadata_error

!==============================================================================

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin FUNCTION metadata_error
!------------------------------------------------------------------------------

if (imodule_errstat /= 0_iintegers) then
  IF (PRESENT(yerrmsg)) yerrmsg = trim(ymodule_errstr)
  metadata_error = imodule_errstat
else
  IF (PRESENT(yerrmsg)) yerrmsg = ''
  metadata_error = 0_iintegers
end if

!------------------------------------------------------------------------------
!  End of the Function
!------------------------------------------------------------------------------

END FUNCTION metadata_error


!==============================================================================
!+ initialize a new metadata storage container
!------------------------------------------------------------------------------

SUBROUTINE metadata_init (tm, imaxbuf, imaxkey, ibuflen, tsrc)
!------------------------------------------------------------------------------
!
! Description:
!   This public procedure initializes a new metadata storage container (tm).
!   To this end it allocates the necessary memory and initializes the fields
!   to their default values. Optionally, the caller can specify the number of
!   metadata elements that can be stored in the container and the size of the
!   metadata storage buffer. If required, another metadata storage (tsrc) can be
!   supplied as an optional argument and the new metadata storage is initialized
!   by copying all content of tsrc to the new metadata storage.
!
! Error numbers: 1-49
!
! Method:
!   If module called for the first time: global init
!   Check arguments
!   Allocate memory
!   Initialize to default values or by copying from an existing metadata storage
!
!==============================================================================

IMPLICIT NONE

! Argument list
TYPE(t_metadata)       , INTENT(INOUT) :: tm      ! metadata storage to initialize
INTEGER         , OPTIONAL, INTENT(IN) :: imaxbuf ! maximum number of buffers to store
INTEGER         , OPTIONAL, INTENT(IN) :: imaxkey ! maximum number of elements to store
INTEGER         , OPTIONAL, INTENT(IN) :: ibuflen ! size of storage buffer
TYPE(t_metadata), OPTIONAL, INTENT(IN) :: tsrc    ! copy values from here

! Local variables
! initialization needed for ifort-compiler, otherwise crash in TRANSFER() below!
INTEGER(KIND=iintegers)  :: iz = 0_iintegers
REAL(KIND=sp)            :: z4 = 0.0_sp
REAL(KIND=dp)            :: z8 = 0.0_dp
CHARACTER(LEN=1)         :: yz = ' ', ybuf(64)
LOGICAL                  :: lz = .TRUE.
INTEGER :: iznumbuf, izmaxbuf, izmaxkey, iznumkey, izbuflen, izerr, izstat
TYPE(t_fld_pointer) :: zp

!==============================================================================

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE metadata_init
!------------------------------------------------------------------------------

imodule_errstat = 0_iintegers

!------------------------------------------------------------------------------
!  Section 1: Module initialization
!------------------------------------------------------------------------------

IF (lfirst_call) THEN
  ! tm%ybuf not allocated yet! Use local character array ybuf instead!
  I_DATA_TYPE_SIZE(I_DATA_TYPE_INTEGER) = SIZE(TRANSFER(iz,ybuf))
  I_DATA_TYPE_SIZE(I_DATA_TYPE_REAL   ) = SIZE(TRANSFER(z4,ybuf))
  I_DATA_TYPE_SIZE(I_DATA_TYPE_DOUBLE ) = SIZE(TRANSFER(z8,ybuf))
  I_DATA_TYPE_SIZE(I_DATA_TYPE_STRING ) = SIZE(TRANSFER(yz,ybuf))
  I_DATA_TYPE_SIZE(I_DATA_TYPE_LOGICAL) = SIZE(TRANSFER(lz,ybuf))
  I_DATA_TYPE_SIZE(I_DATA_TYPE_POINTER) = SIZE(TRANSFER(zp,ybuf))
  lfirst_call = .FALSE.
END IF

!------------------------------------------------------------------------------
!  Section 2: Check arguments
!------------------------------------------------------------------------------

! check optional arguments
IF (PRESENT(imaxbuf)) THEN
  izmaxbuf = imaxbuf
ELSE
  IF (PRESENT(tsrc)) THEN
    izmaxbuf = tsrc%imaxbuf
  ELSE
    izmaxbuf = I_STD_NUMBUF
  END IF
END IF
IF (izmaxbuf < 1_iintegers) THEN
  imodule_errstat = 5_iintegers
  ymodule_errstr = 'metadata_init: number of buffers must be >= 1'
  RETURN
END IF
IF (PRESENT(imaxkey)) THEN
  izmaxkey = imaxkey
ELSE
  IF (PRESENT(tsrc)) THEN
    izmaxkey = tsrc%imaxkeys
  ELSE
    izmaxkey = I_STD_NUMKEY
  END IF
END IF
IF (izmaxkey < 1_iintegers) THEN
  imodule_errstat = 10_iintegers
  ymodule_errstr = 'metadata_init: number of key must be >= 1'
  RETURN
END IF
IF (PRESENT(ibuflen)) THEN
  izbuflen = ibuflen
ELSE
  IF (PRESENT(tsrc)) THEN
    izbuflen = tsrc%imaxbuflen
  ELSE
    izbuflen = I_STD_BUFLEN
  END IF
END IF
IF (izbuflen < 1_iintegers) THEN
  imodule_errstat = 15_iintegers
  ymodule_errstr = 'metadata_init: buffer size must be >= 1'
  RETURN
END IF

! check if already allocated
IF (tm%lisready) THEN
  imodule_errstat = 20_iintegers
  ymodule_errstr = 'metadata_init: metadata storage already initialized'
  RETURN
END IF

! check if tsrc is present and ok
IF (PRESENT(tsrc)) THEN
  IF (.NOT.tsrc%lisready) THEN
    imodule_errstat = 24_iintegers
    ymodule_errstr = 'metadata_init: source metadata storage not initialized'
    RETURN
  END IF
  IF (izmaxbuf < tsrc%inumbuf) THEN
    imodule_errstat = 26_iintegers
    ymodule_errstr = 'metadata_init: not enough buffers specified to copy from source'
    RETURN
  END IF
  IF (izmaxkey < tsrc%inumkey) THEN
    imodule_errstat = 28_iintegers
    ymodule_errstr = 'metadata_init: not enough keys specified to copy from source'
    RETURN
  END IF
  IF (tsrc%inumkey > 0_iintegers) THEN
    IF (izbuflen < tsrc%ipos(tsrc%inumkey)+tsrc%ilen(tsrc%inumkey)) THEN
      imodule_errstat = 30_iintegers
      ymodule_errstr = 'metadata_init: storage buffer not large enough to copy from source'
      RETURN
    END IF
  END IF
END IF

!------------------------------------------------------------------------------
!  Section 3: Allocate memory
!------------------------------------------------------------------------------

izerr = 0_iintegers
ALLOCATE(tm%iidx    (izmaxkey),stat=izstat); izerr = izerr + izstat
ALLOCATE(tm%ykey    (izmaxkey),stat=izstat); izerr = izerr + izstat
ALLOCATE(tm%iattr   (izmaxkey),stat=izstat); izerr = izerr + izstat
ALLOCATE(tm%itype   (izmaxkey),stat=izstat); izerr = izerr + izstat
ALLOCATE(tm%ipos    (izmaxkey),stat=izstat); izerr = izerr + izstat
ALLOCATE(tm%ilen    (izmaxkey),stat=izstat); izerr = izerr + izstat
ALLOCATE(tm%ydefault(izbuflen),stat=izstat); izerr = izerr + izstat
ALLOCATE(tm%ybuf    (izbuflen,izmaxbuf),stat=izstat); izerr = izerr + izstat
IF (izerr /= 0_iintegers) THEN
  imodule_errstat = 40_iintegers
  ymodule_errstr = 'metadata_get: problem allocating memory'
  RETURN
END IF

!------------------------------------------------------------------------------
!  Section 4: Initialize metadata storage
!------------------------------------------------------------------------------

! init dimensional parameters
tm%lisready   = .TRUE.
tm%imaxbuf    = izmaxbuf
tm%imaxkeys   = izmaxkey
tm%imaxbuflen = izbuflen

IF (PRESENT(tsrc)) THEN
  ! copy default values for metadata storage
  iznumbuf                 = tsrc%inumbuf
  iznumkey                 = tsrc%inumkey
  izbuflen                 = tsrc%ipos(iznumkey)+tsrc%ilen(iznumkey)
  tm%inumbuf               = iznumbuf
  tm%inumkey               = iznumkey
  tm%iunique               = tsrc%iunique
  tm%iidx(1:iznumkey)      = tsrc%iidx(1:iznumkey)
  tm%ykey(1:iznumkey)      = tsrc%ykey(1:iznumkey)
  tm%iattr(1:iznumkey)     = tsrc%iattr(1:iznumkey)
  tm%itype(1:iznumkey)     = tsrc%itype(1:iznumkey)
  tm%ipos(1:iznumkey)      = tsrc%ipos(1:iznumkey)
  tm%ilen(1:iznumkey)      = tsrc%ilen(1:iznumkey)
  tm%ydefault(1:izbuflen)  = tsrc%ydefault(1:izbuflen)
  tm%ybuf(1:izbuflen,1:iznumbuf) = tsrc%ybuf(1:izbuflen,1:iznumbuf)
ELSE
  ! setup default values for metadata storage
  tm%inumbuf     = 0_iintegers   ! initially no buffers active
  tm%inumkey     = 0_iintegers   ! initially no keys stored
  tm%iunique     = 42_iintegers
  tm%iidx(:)     = 0_iintegers
  tm%ykey(:)     = REPEAT(' ',I_MAX_KEYLEN)
  tm%iattr(:)    = I_STATUS_DEFAULT
  tm%itype(:)    = I_DATA_TYPE_UNDEF
  tm%ipos(:)     = 0_iintegers
  tm%ilen(:)     = 0_iintegers
  tm%ydefault(:) = ' '
  tm%ybuf(:,:)   = ' '
END IF

!------------------------------------------------------------------------------
!  End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE metadata_init


!==============================================================================
!+ initialize a new metadata storage container
!------------------------------------------------------------------------------

SUBROUTINE metadata_modify (tm, inumbuf)
!------------------------------------------------------------------------------
!
! Description:
!   This public procedure modifies the properties of an existing metadata
!   storage container (tm).
!
! Error numbers: 50-99
!
! Method:
!   Check arguments
!   Make modifications
!
!==============================================================================

IMPLICIT NONE

! Argument list
TYPE(t_metadata), INTENT(INOUT) :: tm         ! metadata storage to initialize
INTEGER, OPTIONAL, INTENT(IN) :: inumbuf      ! number of active buffers to store

! Local variables
INTEGER :: iz
INTEGER :: iznumbuf

!==============================================================================

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE metadata_modify
!------------------------------------------------------------------------------

imodule_errstat = 0_iintegers

!------------------------------------------------------------------------------
!  Section 1: Check arguments
!------------------------------------------------------------------------------

! check if ready
IF (.NOT. tm%lisready) THEN
  imodule_errstat = 52_iintegers
  ymodule_errstr = 'metadata_modify: metadata storage not initialized'
  RETURN
END IF

! check optional arguments
IF (PRESENT(inumbuf)) THEN
  IF (inumbuf < 1_iintegers) THEN
    imodule_errstat = 55_iintegers
    ymodule_errstr = 'metadata_modify: number of active buffers must be >= 1'
    RETURN
  END IF
  IF (inumbuf > tm%imaxbuf) THEN
    imodule_errstat = 57_iintegers
    ymodule_errstr = 'metadata_modify: maximum number of active buffers exceeded'
    RETURN
  END IF
END IF

!------------------------------------------------------------------------------
!  Section 2: Modify number of active buffers
!------------------------------------------------------------------------------

IF (PRESENT(inumbuf)) THEN
  ! if number of buffers is being increased, initialize with default values
  iznumbuf = tm%inumbuf
  IF (iznumbuf < inumbuf) THEN
    DO iz = iznumbuf+1, inumbuf
      tm%ybuf(:,iz) = tm%ydefault(:)
    ENDDO
  ENDIF
  tm%inumbuf = inumbuf
END IF

!------------------------------------------------------------------------------
!  End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE metadata_modify


!==============================================================================
!+ Module procedure to dispose of metadata storage
!------------------------------------------------------------------------------

SUBROUTINE metadata_finish(tm)
!------------------------------------------------------------------------------
!
! Description:
!   This public procedure diposes of a metadata storage container. To this
!   end it deallocates the memory and resets internal variables to mark the
!   metadata container as unused.
!
! Error numbers: 100-199
!
! Method:
!   Check arguments
!   Deallocate memory
!   Reset metadata storage
!
!==============================================================================

IMPLICIT NONE

! Argument list
TYPE(t_metadata), INTENT(INOUT) :: tm

! Local variables
INTEGER :: izerr, izstat

!==============================================================================

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE metadata_finish
!------------------------------------------------------------------------------

imodule_errstat = 0_iintegers

!------------------------------------------------------------------------------
!  Section 1: Check arguments
!------------------------------------------------------------------------------

! check if already allocated
IF (.NOT. tm%lisready) THEN
  imodule_errstat = 100_iintegers
  ymodule_errstr = 'metadata_finish: metadata storage is not initialized'
  RETURN
END IF

!------------------------------------------------------------------------------
!  Section 2: Deallocate memory
!------------------------------------------------------------------------------

izerr = 0_iintegers
DEALLOCATE(tm%iidx    ,stat=izstat); izerr = izerr + izstat
DEALLOCATE(tm%ykey    ,stat=izstat); izerr = izerr + izstat
DEALLOCATE(tm%iattr   ,stat=izstat); izerr = izerr + izstat
DEALLOCATE(tm%itype   ,stat=izstat); izerr = izerr + izstat
DEALLOCATE(tm%ipos    ,stat=izstat); izerr = izerr + izstat
DEALLOCATE(tm%ilen    ,stat=izstat); izerr = izerr + izstat
DEALLOCATE(tm%ydefault,stat=izstat); izerr = izerr + izstat
DEALLOCATE(tm%ybuf    ,stat=izstat); izerr = izerr + izstat
IF (izerr /= 0_iintegers) THEN
  imodule_errstat = 110_iintegers
  ymodule_errstr = 'metadata_finish: problem deallocating memory'
  RETURN
END IF

!------------------------------------------------------------------------------
!  Section 3: Reset metadata storage
!------------------------------------------------------------------------------

tm%lisready = .FALSE.
tm%inumkey  = 0_iintegers

!------------------------------------------------------------------------------
!  End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE metadata_finish


!==============================================================================
!+ Module procedures to delete an item from metadata storage
!------------------------------------------------------------------------------

SUBROUTINE meta_delete_byname(tm,yname)
!------------------------------------------------------------------------------
!
! Description:
!   This public procedure deletes a metadata_item in the supplied metadata
!   storage (tm)
!
! Error numbers: 300-349
!
! Method:
!   Find metadata item in storage buffer
!   Shift data in buffer to delete item
!
!==============================================================================

IMPLICIT NONE

! Argument list
TYPE(t_metadata), INTENT(INOUT) :: tm
CHARACTER(LEN=*), INTENT(IN) :: yname

! Local variables
INTEGER :: izcurkey, iznumkey, izcurpos, iznextpos, izendpos
CHARACTER(LEN=I_MAX_KEYLEN+7) :: yzwhere

!==============================================================================

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE meta_delete_byname
!------------------------------------------------------------------------------

imodule_errstat = 0_iintegers
yzwhere = ' (key=' // TRIM(yname) // ')'

!------------------------------------------------------------------------------
!  Section 1: Find metadata item in storage buffer
!------------------------------------------------------------------------------

! check if metadata storage is ready for use
IF (.NOT. tm%lisready) THEN
  imodule_errstat = 300_iintegers
  ymodule_errstr = 'metadata_delete: metadata storage is not initialized' // yzwhere
  RETURN
END IF

! check string length
IF (LEN_TRIM(yname) > I_MAX_KEYLEN) THEN
  imodule_errstat = 310_iintegers
  ymodule_errstr = 'metadata_delete: key string too long' // yzwhere
  RETURN
END IF

! find key by name
DO izcurkey = 1, tm%inumkey
  IF (TRIM(tm%ykey(izcurkey)) == TRIM(yname)) EXIT
END DO

! check if key is present
IF (izcurkey < 1_iintegers .OR. izcurkey > tm%inumkey) THEN
  imodule_errstat = 320_iintegers
  ymodule_errstr = 'metadata_delete: key not found or out of range' // yzwhere
  RETURN
END IF

! check if delete protected
IF (IAND(tm%iattr(izcurkey),I_STATUS_DELETE) == 0_iintegers) THEN
  ! key exists, but is delete protected -> error
  imodule_errstat = 330_iintegers
  ymodule_errstr = 'metadata_delete: key is delete protected' // yzwhere
  RETURN
END IF

!------------------------------------------------------------------------------
!  Section 2: Shift trailing elements to delete item
!------------------------------------------------------------------------------

! check if last entry
IF (izcurkey == tm%inumkey) THEN
  tm%inumkey = tm%inumkey - 1_iintegers
  RETURN
END IF

! compute amount of shifting required
iznumkey  = tm%inumkey
izcurpos  = tm%ipos(izcurkey)
iznextpos = tm%ipos(izcurkey+1)
izendpos  = tm%ipos(tm%inumkey) + tm%ilen(tm%inumkey)

! shift trailing elements forward
tm%inumkey                       = tm%inumkey - 1_iintegers
tm%iidx    (izcurkey:iznumkey-1) = tm%iidx    (izcurkey+1:iznumkey)
tm%ykey    (izcurkey:iznumkey-1) = tm%ykey    (izcurkey+1:iznumkey)
tm%iattr   (izcurkey:iznumkey-1) = tm%iattr   (izcurkey+1:iznumkey)
tm%itype   (izcurkey:iznumkey-1) = tm%itype   (izcurkey+1:iznumkey)
tm%ipos    (izcurkey:iznumkey-1) = tm%ipos    (izcurkey+1:iznumkey)-tm%ilen(izcurkey)
tm%ilen    (izcurkey:iznumkey-1) = tm%ilen    (izcurkey+1:iznumkey)
tm%ydefault(izcurpos:izcurpos+izendpos-iznextpos) = tm%ydefault(iznextpos:izendpos)
tm%ybuf    (izcurpos:izcurpos+izendpos-iznextpos,1:tm%inumbuf) &
  = tm%ybuf(iznextpos:izendpos,1:tm%inumbuf)

!------------------------------------------------------------------------------
!  End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE meta_delete_byname


SUBROUTINE meta_delete_byindex(tm,iidx)
!------------------------------------------------------------------------------
!
! Description:
!   This public procedure deletes a metadata_item in the supplied metadata
!   storage (tm)
!
! Error numbers: 350-399
!
! Method:
!   Find metadata item in storage buffer
!   Shift data in buffer to delete item
!
!==============================================================================

IMPLICIT NONE

! Argument list
TYPE(t_metadata), INTENT(INOUT) :: tm
INTEGER, INTENT(IN) :: iidx

! Local variables
INTEGER :: izcurkey, iznumkey, izcurpos, iznextpos, izendpos
CHARACTER(LEN=I_MAX_KEYLEN+7) :: yzwhere

!==============================================================================

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE meta_delete_byindex
!------------------------------------------------------------------------------

imodule_errstat = 0_iintegers
write(yzwhere,*) iidx
yzwhere = ' (key=' // TRIM(ADJUSTL(yzwhere)) // ')'

!------------------------------------------------------------------------------
!  Section 1: Find metadata item in storage buffer
!------------------------------------------------------------------------------

! check if metadata storage is ready for use
IF (.NOT. tm%lisready) THEN
  imodule_errstat = 350_iintegers
  ymodule_errstr = 'metadata_delete: metadata storage is not initialized' // yzwhere
  RETURN
END IF

! find key by index
DO izcurkey = 1, tm%inumkey
  IF (tm%iidx(izcurkey) == iidx) EXIT
END DO

! check if key is present
IF (izcurkey < 1 .OR. izcurkey > tm%inumkey) THEN
  imodule_errstat = 370_iintegers
  ymodule_errstr = 'metadata_delete: key not found or out of range' // yzwhere
  RETURN
END IF

! check if delete protected
IF (IAND(tm%iattr(izcurkey),I_STATUS_DELETE) == 0_iintegers) THEN
  ! key exists, but is delete protected -> error
  imodule_errstat = 380_iintegers
  ymodule_errstr = 'metadata_delete: key is delete protected' // yzwhere
  RETURN
END IF

!------------------------------------------------------------------------------
!  Section 2: Shift trailing elements to delete item
!------------------------------------------------------------------------------

! check if last entry
IF (izcurkey == tm%inumkey) THEN
  tm%inumkey = tm%inumkey - 1
  RETURN
END IF

! shift trailing elements forward
iznumkey  = tm%inumkey
izcurpos  = tm%ipos(izcurkey)
iznextpos = tm%ipos(izcurkey+1)
izendpos  = tm%ipos(tm%inumkey) + tm%ilen(tm%inumkey)
tm%inumkey                       = tm%inumkey - 1_iintegers
tm%iidx    (izcurkey:iznumkey-1) = tm%iidx    (izcurkey+1:iznumkey)
tm%ykey    (izcurkey:iznumkey-1) = tm%ykey    (izcurkey+1:iznumkey)
tm%iattr   (izcurkey:iznumkey-1) = tm%iattr   (izcurkey+1:iznumkey)
tm%itype   (izcurkey:iznumkey-1) = tm%itype   (izcurkey+1:iznumkey)
tm%ipos    (izcurkey:iznumkey-1) = tm%ipos    (izcurkey+1:iznumkey)-tm%ilen(izcurkey)
tm%ilen    (izcurkey:iznumkey-1) = tm%ilen    (izcurkey+1:iznumkey)
tm%ydefault(izcurpos:izcurpos+izendpos-iznextpos) = tm%ydefault(iznextpos:izendpos)
tm%ybuf    (izcurpos:izcurpos+izendpos-iznextpos,1:tm%inumbuf) &
  = tm%ybuf(iznextpos:izendpos,1:tm%inumbuf)

!------------------------------------------------------------------------------
!  End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE meta_delete_byindex


!==============================================================================
!+ Module procedure to print content of metadata storage to a I/O unit
!------------------------------------------------------------------------------

SUBROUTINE metadata_print(tm,ibuf,iounit,iindent)
!------------------------------------------------------------------------------
!
! Description:
!   This public procedure prints the contents of a metadata storage buffer
!   to a user specified I/O unit in pretty print format.
!
! Error numbers: 400-499
!
! Method:
!
!==============================================================================

IMPLICIT NONE

! Argument list
TYPE(t_metadata) , INTENT(INOUT) :: tm
INTEGER, OPTIONAL, INTENT(IN)    :: ibuf
INTEGER, OPTIONAL, INTENT(IN)    :: iounit
INTEGER, OPTIONAL, INTENT(IN)    :: iindent

! Local variables
INTEGER                     :: iz
REAL(KIND=sp)               :: zr4
REAL(KIND=dp)               :: zr8
LOGICAL                     :: lz
CHARACTER(LEN=I_MAX_STRLEN) :: yy
CHARACTER(LEN=12)           :: yindentf
CHARACTER(LEN=18)           :: ykeystr
CHARACTER(LEN=12)           :: ytypestr
CHARACTER(LEN=8)            :: yattr
CHARACTER(LEN=12)           :: yvalstr1
INTEGER :: izcurkey, izunit
INTEGER :: izerr, izstat, izindent, izbuf, izbufstart, izbufend

!==============================================================================

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE metadata_print
!------------------------------------------------------------------------------

imodule_errstat = 0_iintegers

! check if metadata storage is ready for use
IF (.NOT. tm%lisready) THEN
  imodule_errstat = 400_iintegers
  ymodule_errstr = 'metadata_print: metadata storage is not initialized'
  RETURN
END IF

! setup buffer entries to print
IF (PRESENT(ibuf)) THEN
  izbufstart = ibuf
  izbufend = ibuf
ELSE
  izbufstart = 1_iintegers
  izbufend = tm%inumbuf
END IF

! check if requested buffer indices are ok
IF (izbufstart < 1_iintegers .OR. izbufend > tm%inumbuf) THEN
  imodule_errstat = 410_iintegers
  ymodule_errstr = 'metadata_print: requested buffer index is out of range'
  RETURN
END IF

! check if output unit is present, otherwise use STDOUT
IF (PRESENT(iounit)) THEN
  izunit = iounit
ELSE
  izunit = 6_iintegers
END IF

! check if indent is present, other use default
IF (PRESENT(iindent)) THEN
  izindent = iindent
ELSE
  izindent = 6_iintegers
END IF

! prepare indent format string
WRITE(yindentf,'(i5)') izindent
yindentf = '('//trim(adjustl(yindentf))//'x)'

DO izbuf = izbufstart, izbufend

  ! construct key string including buffer number
  WRITE(ykeystr,*) izbuf
  ykeystr = ADJUSTL(TRIM(ykeystr))
  ykeystr = '#' // TRIM(ykeystr) // REPEAT('.',7-LEN_TRIM(ykeystr)) // 'KEY.......'
  
  ! print header
  WRITE(izunit,yindentf,advance="no")
  WRITE(izunit,'(1x,a18,2x,a30,2x,a12,2x,a8)') &
    ykeystr, &
    '.............VALUE............', &
    '....TYPE....'!!! , '..ATTR..'
  
  ! loop over all keys
  DO izcurkey = 1, tm%inumkey
    IF (LEN_TRIM(tm%ykey(izcurkey))>18_iintegers) THEN
      ykeystr = tm%ykey(izcurkey)(1:15)//'...'
    ELSE
      ykeystr = ADJUSTL(tm%ykey(izcurkey))
    END IF
    ytypestr = ADJUSTL(Y_DATA_TYPE(tm%itype(izcurkey)))
    yattr = ''
!!!    IF (IAND(tm%iattr(izcurkey),I_STATUS_DEFINED) /= 0) THEN
!!!      yattr = TRIM(yattr) // 's'
!!!    ELSE
!!!      yattr = TRIM(yattr) // '-'
!!!    END IF
!!!    IF (IAND(tm%iattr(izcurkey),I_STATUS_READ) /= 0) THEN
!!!      yattr = TRIM(yattr) // 'r'
!!!    ELSE
!!!      yattr = TRIM(yattr) // '-'
!!!    END IF
!!!    IF (IAND(tm%iattr(izcurkey),I_STATUS_WRITE) /= 0) THEN
!!!      yattr = TRIM(yattr) // 'w'
!!!    ELSE
!!!      yattr = TRIM(yattr) // '-'
!!!    END IF
!!!    IF (IAND(tm%iattr(izcurkey),I_STATUS_DELETE) /= 0) THEN
!!!      yattr = TRIM(yattr) // 'd'
!!!    ELSE
!!!      yattr = TRIM(yattr) // '-'
!!!    END IF
    yattr = ADJUSTL(yattr)
    WRITE(izunit,yindentf,advance="no")
    WRITE(izunit,'(1x,a18,2x)',advance="no") ykeystr
    SELECT CASE (tm%itype(izcurkey))
      CASE (I_DATA_TYPE_UNDEF)
      CASE (I_DATA_TYPE_INTEGER)
        CALL metadata_get(tm,izbuf,tm%ykey(izcurkey),iz)
        IF (imodule_errstat /= 0_iintegers) RETURN
        WRITE(yvalstr1,'(i10)') iz
        WRITE(izunit,'(a30)',advance="no") TRIM(yvalstr1)
      CASE (I_DATA_TYPE_REAL)
        CALL metadata_get(tm,izbuf,tm%ykey(izcurkey),zr4)
        IF (imodule_errstat /= 0_iintegers) RETURN
        WRITE(yvalstr1,'(e12.4)') zr4
        WRITE(izunit,'(a30)',advance="no") TRIM(yvalstr1)
      CASE (I_DATA_TYPE_DOUBLE)
        CALL metadata_get(tm,izbuf,tm%ykey(izcurkey),zr8)
        IF (imodule_errstat /= 0_iintegers) RETURN
        WRITE(yvalstr1,'(e12.4)') zr8
        WRITE(izunit,'(a30)',advance="no") TRIM(yvalstr1)
      CASE (I_DATA_TYPE_STRING)
        CALL metadata_get(tm,izbuf,tm%ykey(izcurkey),yy)
        IF (imodule_errstat /= 0_iintegers) RETURN
        WRITE(izunit,'(a30)',advance="no") "'"//TRIM(yy)//"'"
      CASE (I_DATA_TYPE_LOGICAL)
        CALL metadata_get(tm,izbuf,tm%ykey(izcurkey),lz)
        IF (imodule_errstat /= 0_iintegers) RETURN
        IF (lz) THEN
          yvalstr1 = 'TRUE'
        ELSE
          yvalstr1 = 'FALSE'
        END IF
        WRITE(izunit,'(a30)',advance="no") TRIM(yvalstr1)
      CASE (I_DATA_TYPE_POINTER)
        yvalstr1 = ''
        WRITE(izunit,'(a30)',advance="no") TRIM(yvalstr1)
    END SELECT
    WRITE(izunit,'(2x,a12,2x,a8)',advance="no") ytypestr,yattr
    WRITE(izunit,'(a)') ''
  END DO
  
  ! print footer
  WRITE(izunit,yindentf,advance="no")
  WRITE(izunit,'(a79)') ' '//REPEAT(' ',78)

END DO ! loop over buffer indices

!------------------------------------------------------------------------------
!  End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE metadata_print


!==============================================================================
!+ Module procedures to define metadata item
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!
! Description:
!   The following procedures overload the procedure metadata_define which is
!   called upon a request to define a data item in the metadata storage. For
!   each data type a separate procedure  has to be created. The procedures
!   all have exactly the same structure with the exception of type of the ydefault
!   argument as well as some local type-dependent parameters. The procedures rely
!   heavily on the meta_define functions which do most of the work.
!
! Subroutines:
!
!   meta_define_integer
!   meta_define_real4
!   meta_define_real8
!   meta_define_string
!   meta_define_logical
!   meta_define_pointer2
!   meta_define_pointer3
!   meta_define_pointer4
!
! Error numbers: 1000-1099
!
! Method:
!   Check define request
!   Define metadata entry
!
!==============================================================================

SUBROUTINE meta_define_integer(tm,yname,ydefault,iidx,lprotect)
  IMPLICIT NONE
  TYPE(t_metadata) , INTENT(INOUT) :: tm
  CHARACTER(LEN=*) , INTENT(IN)    :: yname
  INTEGER          , INTENT(IN)    :: ydefault
  INTEGER, OPTIONAL, INTENT(OUT)   :: iidx
  LOGICAL, OPTIONAL, INTENT(IN)    :: lprotect
  INTEGER :: izkey, izsize, izcurpos
  izsize = SIZE(TRANSFER(ydefault, tm%ybuf))
  izcurpos = meta_define_check(tm, yname, I_DATA_TYPE_INTEGER, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_define_store(tm, izkey, yname, I_DATA_TYPE_INTEGER, &
         izcurpos, izsize, TRANSFER(ydefault, tm%ybuf), lprotect=lprotect)
  IF (PRESENT(iidx)) THEN
    iidx = tm%iidx(izkey)
  END IF
END SUBROUTINE meta_define_integer

SUBROUTINE meta_define_real4(tm,yname,ydefault,iidx,lprotect)
  IMPLICIT NONE
  TYPE(t_metadata) , INTENT(INOUT) :: tm
  CHARACTER(LEN=*) , INTENT(IN)    :: yname
  REAL(KIND=sp)    , INTENT(IN)    :: ydefault
  INTEGER, OPTIONAL, INTENT(OUT)   :: iidx
  LOGICAL, OPTIONAL, INTENT(IN)    :: lprotect
  INTEGER :: izkey, izsize, izcurpos
  izsize = SIZE(TRANSFER(ydefault, tm%ybuf))
  izcurpos = meta_define_check(tm, yname, I_DATA_TYPE_REAL, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_define_store(tm, izkey, yname, I_DATA_TYPE_REAL, &
         izcurpos, izsize, TRANSFER(ydefault, tm%ybuf), lprotect=lprotect)
  IF (PRESENT(iidx)) THEN
    iidx = tm%iidx(izkey)
  END IF
END SUBROUTINE meta_define_real4

SUBROUTINE meta_define_real8(tm,yname,ydefault,iidx,lprotect)
  IMPLICIT NONE
  TYPE(t_metadata) , INTENT(INOUT) :: tm
  CHARACTER(LEN=*) , INTENT(IN)    :: yname
  REAL(KIND=dp)    , INTENT(IN)    :: ydefault
  INTEGER, OPTIONAL, INTENT(OUT)   :: iidx
  LOGICAL, OPTIONAL, INTENT(IN)    :: lprotect
  INTEGER :: izkey, izsize, izcurpos
  izsize = SIZE(TRANSFER(ydefault, tm%ybuf))
  izcurpos = meta_define_check(tm, yname, I_DATA_TYPE_DOUBLE, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_define_store(tm, izkey, yname, I_DATA_TYPE_DOUBLE, &
         izcurpos, izsize, TRANSFER(ydefault, tm%ybuf), lprotect=lprotect)
  IF (PRESENT(iidx)) THEN
    iidx = tm%iidx(izkey)
  END IF
END SUBROUTINE meta_define_real8

SUBROUTINE meta_define_string(tm,yname,ydefault,iidx,lprotect)
  IMPLICIT NONE
  TYPE(t_metadata) , INTENT(INOUT) :: tm
  CHARACTER(LEN=*) , INTENT(IN)    :: yname
  CHARACTER(LEN=*) , INTENT(IN)    :: ydefault
  CHARACTER(LEN=I_MAX_STRLEN)      :: yzstr
  INTEGER, OPTIONAL, INTENT(OUT)   :: iidx
  LOGICAL, OPTIONAL, INTENT(IN)    :: lprotect
  INTEGER :: izkey, izsize, izcurpos
  izsize = SIZE(TRANSFER(TRIM(ydefault), tm%ybuf))
  izcurpos = meta_define_check(tm, yname, I_DATA_TYPE_STRING, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  yzstr = TRIM(ydefault)
  CALL meta_define_store(tm, izkey, yname, I_DATA_TYPE_STRING, &
         izcurpos, izsize, TRANSFER(yzstr, tm%ybuf), lprotect=lprotect)
  IF (PRESENT(iidx)) THEN
    iidx = tm%iidx(izkey)
  END IF
END SUBROUTINE meta_define_string

SUBROUTINE meta_define_logical(tm,yname,ydefault,iidx,lprotect)
  IMPLICIT NONE
  TYPE(t_metadata) , INTENT(INOUT) :: tm
  CHARACTER(LEN=*) , INTENT(IN)    :: yname
  LOGICAL          , INTENT(IN)    :: ydefault
  INTEGER, OPTIONAL, INTENT(OUT)   :: iidx
  LOGICAL, OPTIONAL, INTENT(IN)    :: lprotect
  INTEGER :: izkey, izsize, izcurpos
  izsize = SIZE(TRANSFER(ydefault, tm%ybuf))
  izcurpos = meta_define_check(tm, yname, I_DATA_TYPE_LOGICAL, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_define_store(tm, izkey, yname, I_DATA_TYPE_LOGICAL, &
         izcurpos, izsize, TRANSFER(ydefault, tm%ybuf), lprotect=lprotect)
  IF (PRESENT(iidx)) THEN
    iidx = tm%iidx(izkey)
  END IF
END SUBROUTINE meta_define_logical

SUBROUTINE meta_define_pointer2(tm,yname,ydefault,iidx,lprotect)
  IMPLICIT NONE
  TYPE(t_metadata)      , INTENT(INOUT) :: tm
  CHARACTER(LEN=*)      , INTENT(IN)    :: yname
  REAL(KIND=wp), POINTER, INTENT(IN)    :: ydefault(:,:)
  INTEGER, OPTIONAL     , INTENT(OUT)   :: iidx
  LOGICAL, OPTIONAL     , INTENT(IN)    :: lprotect
  INTEGER :: izkey, izsize, izcurpos
  TYPE(t_fld_pointer) :: zp
  imodule_errstat = 0_iintegers
  IF (ASSOCIATED(ydefault)) THEN
    IF (I_PTR_MAXDIM < 2) THEN
      imodule_errstat = 1000_iintegers
      ymodule_errstr  = 'meta_define_pointer: maximal number of dimensions is not sufficient (<2)'
      RETURN
    ENDIF
    IF (ANY(lbound(ydefault) /= 1_iintegers)) THEN
      imodule_errstat = 1001_iintegers
      ymodule_errstr  = 'meta_define_pointer: some pointer lower bound are not = 0'
      RETURN
    ENDIF
    zp%fld_cptr = C_LOC(ydefault(1,1))
    zp%fld_rank = 2_iintegers
    zp%fld_shape(:) = -1_iintegers
    zp%fld_shape(1:2) = ubound(ydefault) - lbound(ydefault) + 1_iintegers
  ELSE
    zp%fld_cptr = C_NULL_PTR
  END IF
  izsize = SIZE(TRANSFER(zp, tm%ybuf))
  izcurpos = meta_define_check(tm, yname, I_DATA_TYPE_POINTER, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_define_store(tm, izkey, yname, I_DATA_TYPE_POINTER, &
         izcurpos, izsize, TRANSFER(zp, tm%ybuf), lprotect=lprotect)
  IF (PRESENT(iidx)) THEN
    iidx = tm%iidx(izkey)
  END IF
END SUBROUTINE meta_define_pointer2

SUBROUTINE meta_define_pointer3(tm,yname,ydefault,iidx,lprotect)
  IMPLICIT NONE
  TYPE(t_metadata)      , INTENT(INOUT) :: tm
  CHARACTER(LEN=*)      , INTENT(IN)    :: yname
  REAL(KIND=wp), POINTER, INTENT(IN)    :: ydefault(:,:,:)
  INTEGER, OPTIONAL     , INTENT(OUT)   :: iidx
  LOGICAL, OPTIONAL     , INTENT(IN)    :: lprotect
  INTEGER :: izkey, izsize, izcurpos
  TYPE(t_fld_pointer) :: zp
  imodule_errstat = 0_iintegers
  IF (ASSOCIATED(ydefault)) THEN
    IF (I_PTR_MAXDIM < 3) THEN
      imodule_errstat = 1005_iintegers
      ymodule_errstr  = 'meta_define_pointer: maximal number of dimensions is not sufficient (<3)'
      RETURN
    ENDIF
    IF (ANY(lbound(ydefault) /= 1_iintegers)) THEN
      imodule_errstat = 1006_iintegers
      ymodule_errstr  = 'meta_define_pointer: some pointer lower bound are not = 0'
      RETURN
    ENDIF
    zp%fld_cptr = C_LOC(ydefault(1,1,1))
    zp%fld_rank = 3_iintegers
    zp%fld_shape(:) = -1_iintegers
    zp%fld_shape(1:3) = ubound(ydefault) - lbound(ydefault) + 1_iintegers
  ELSE
    zp%fld_cptr = C_NULL_PTR
  END IF
  izsize = SIZE(TRANSFER(zp, tm%ybuf))
  izcurpos = meta_define_check(tm, yname, I_DATA_TYPE_POINTER, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_define_store(tm, izkey, yname, I_DATA_TYPE_POINTER, &
         izcurpos, izsize, TRANSFER(zp, tm%ybuf), lprotect=lprotect)
  IF (PRESENT(iidx)) THEN
    iidx = tm%iidx(izkey)
  END IF
END SUBROUTINE meta_define_pointer3

SUBROUTINE meta_define_pointer4(tm,yname,ydefault,iidx,lprotect)
  IMPLICIT NONE
  TYPE(t_metadata)      , INTENT(INOUT) :: tm
  CHARACTER(LEN=*)      , INTENT(IN)    :: yname
  REAL(KIND=wp), POINTER, INTENT(IN)    :: ydefault(:,:,:,:)
  INTEGER, OPTIONAL     , INTENT(OUT)   :: iidx
  LOGICAL, OPTIONAL     , INTENT(IN)    :: lprotect
  INTEGER :: izkey, izsize, izcurpos
  TYPE(t_fld_pointer) :: zp
  imodule_errstat = 0_iintegers
  IF (ASSOCIATED(ydefault)) THEN
    IF (I_PTR_MAXDIM < 4) THEN
      imodule_errstat = 1010_iintegers
      ymodule_errstr  = 'meta_define_pointer: maximal number of dimensions is not sufficient (<4)'
      RETURN
    ENDIF
    IF (ANY(lbound(ydefault) /= 1_iintegers)) THEN
      imodule_errstat = 1011_iintegers
      ymodule_errstr  = 'meta_define_pointer: some pointer lower bound are not = 0'
      RETURN
    ENDIF
    zp%fld_cptr = C_LOC(ydefault(1,1,1,1))
    zp%fld_rank = 4_iintegers
    zp%fld_shape(:) = -1_iintegers
    zp%fld_shape(1:4) = ubound(ydefault) - lbound(ydefault) + 1_iintegers
  ELSE
    zp%fld_cptr = C_NULL_PTR
  END IF
  izsize = SIZE(TRANSFER(zp, tm%ybuf))
  izcurpos = meta_define_check(tm, yname, I_DATA_TYPE_POINTER, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_define_store(tm, izkey, yname, I_DATA_TYPE_POINTER, &
         izcurpos, izsize, TRANSFER(zp, tm%ybuf), lprotect=lprotect)
  IF (PRESENT(iidx)) THEN
    iidx = tm%iidx(izkey)
  END IF
END SUBROUTINE meta_define_pointer4

!==============================================================================
!+ Module function to check data before storing it into metadata storage
!------------------------------------------------------------------------------

FUNCTION meta_define_check(tm,yname,itype,isize,ikey)
!------------------------------------------------------------------------------
!
! Description:
!   This private function checks a user request to store a metadata item
!   in a metadata storage. It executes several consistency checks and if
!   successfull returns the position in the buffer where to store the data
!   item.
!
! Error numbers: 1100-1199
!
! Method:
!   Check store request
!
!==============================================================================

IMPLICIT NONE

! Argument list
TYPE(t_metadata), INTENT(IN)    :: tm
CHARACTER(LEN=*), INTENT(IN)    :: yname
INTEGER         , INTENT(IN)    :: itype
INTEGER         , INTENT(INOUT) :: isize
INTEGER         , INTENT(OUT)   :: ikey

! Return value
INTEGER :: meta_define_check

! Local variables
INTEGER :: iz, izcurpos, izcurkey
CHARACTER(LEN=I_MAX_KEYLEN+7) :: yzwhere

!==============================================================================

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE meta_define_check
!------------------------------------------------------------------------------

imodule_errstat = 0_iintegers
yzwhere = ' (key=' // TRIM(yname) // ')'

!------------------------------------------------------------------------------
!  Section 1: Check store request
!------------------------------------------------------------------------------

! check if metadata storage is ready for use
IF (.NOT. tm%lisready) THEN
  imodule_errstat = 1100_iintegers
  ymodule_errstr = 'metadata_define: metadata storage is not initialized' // yzwhere
  RETURN
END IF

! check key string length
IF (LEN_TRIM(yname) > I_MAX_KEYLEN) THEN
  imodule_errstat = 1110_iintegers
  ymodule_errstr = 'metadata_define: key string too long' // yzwhere
  RETURN
END IF

! check for empty arrays
IF (isize == 0_iintegers) THEN
  imodule_errstat = 1130_iintegers
  ymodule_errstr = 'metadata_define: cannot store zero sized arrays' // yzwhere
  RETURN
END IF

! check data size
IF (itype /= I_DATA_TYPE_STRING) THEN
  IF (isize /= I_DATA_TYPE_SIZE(itype)) THEN
    imodule_errstat = 1140_iintegers
    ymodule_errstr = 'metadata_define: size mismatch' // yzwhere
    RETURN
  END IF
ELSE
  IF (MOD(isize,I_DATA_TYPE_SIZE(itype)) /= 0_iintegers) THEN
    imodule_errstat = 1140_iintegers
    ymodule_errstr = 'metadata_define: size mismatch' // yzwhere
    RETURN
  END IF
  IF (isize > I_MAX_STRLEN*I_DATA_TYPE_SIZE(itype)) THEN
    imodule_errstat = 1142_iintegers
    ymodule_errstr = 'metadata_define: string length exceeds maximum (I_MAX_STRLEN)' &
      // yzwhere
    RETURN
  END IF
  isize = I_MAX_STRLEN*I_DATA_TYPE_SIZE(itype)
END IF

! check if key is already taken
DO iz = 1, tm%inumkey
  IF (TRIM(tm%ykey(iz)) == TRIM(yname)) THEN
    ! key exists, don't overwrite -> error
    imodule_errstat = 1155_iintegers
    ymodule_errstr = 'metadata_define: key already in use' // yzwhere
    RETURN
  END IF
END DO

! check if keys left
izcurkey = tm%inumkey
IF (izcurkey+1 > tm%imaxkeys) THEN
  imodule_errstat = 1120_iintegers
  ymodule_errstr = 'metadata_define: maximum number of keys exceeded' // yzwhere
  RETURN
END IF

! check if enough space left
izcurpos = 1_iintegers
IF (izcurkey > 0_iintegers) izcurpos = tm%ipos(izcurkey) + tm%ilen(izcurkey)
IF (isize+izcurpos > tm%imaxbuflen) THEN
  imodule_errstat = 1160
  ymodule_errstr = 'metadata_define: storage buffer size exceeded' // yzwhere
  RETURN
END IF

! append to end of metadata storage
iz = izcurkey + 1_iintegers

! ok, return key number and buffer position
ikey = iz
meta_define_check = izcurpos

!------------------------------------------------------------------------------
!  End of the Subroutine
!------------------------------------------------------------------------------

END FUNCTION meta_define_check

!==============================================================================
!+ Module procedure to store a data item in a metadata storage
!------------------------------------------------------------------------------

SUBROUTINE meta_define_store(tm,ikey,yname,itype,ipos,isize,ydata,lprotect)
!------------------------------------------------------------------------------
!
! Description:
!   This private procedure stores a metadata item which has been checked for
!   consistency into a metadata storage. It requires the storage position and
!   data item size.
!
! Error numbers: 1300-1399
!
! Method:
!   Store data in buffer
!
!==============================================================================

IMPLICIT NONE

! Argument list
TYPE(t_metadata) , INTENT(INOUT) :: tm
INTEGER          , INTENT(IN)    :: ikey
CHARACTER(LEN=*) , INTENT(IN)    :: yname
INTEGER          , INTENT(IN)    :: itype
INTEGER          , INTENT(IN)    :: ipos
INTEGER          , INTENT(IN)    :: isize
CHARACTER(LEN=1) , INTENT(IN)    :: ydata(:)
LOGICAL, OPTIONAL, INTENT(IN)    :: lprotect

! Local variables
INTEGER :: izattr, iz

!==============================================================================

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE meta_define_store
!------------------------------------------------------------------------------

imodule_errstat = 0_iintegers

! check optional arguments
IF (PRESENT(lprotect)) THEN
  IF (lprotect) THEN
    izattr = I_STATUS_DEFAULT - I_STATUS_DELETE
  ELSE
    izattr = I_STATUS_DEFAULT
  END IF
ELSE
  izattr = I_STATUS_DEFAULT
END IF

!------------------------------------------------------------------------------
!  Section 1: Store data in buffer
!------------------------------------------------------------------------------

tm%inumkey        = ikey
tm%iidx(ikey)     = tm%iunique; tm%iunique = tm%iunique + 1_iintegers
tm%ykey(ikey)     = TRIM(yname)
tm%iattr(ikey)    = izattr
tm%itype(ikey)    = itype
tm%ipos(ikey)     = ipos
tm%ilen(ikey)     = isize
tm%ydefault(ipos:ipos+isize-1)  = ydata
DO iz = 1, tm%inumbuf
  tm%ybuf(ipos:ipos+isize-1,iz) = ydata
END DO

!------------------------------------------------------------------------------
!  End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE meta_define_store


!==============================================================================
!+ Module procedures to store data into metadata storage
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!
! Description:
!   The following procedures overload the procedure metadata_set which is
!   called upon a request to store a data item in the metadata storage. For
!   each data type as well as for scalar/vector data a separate procedure
!   has to be created. The procedures all have exactly the same structure with
!   the exception of dimension and type of the ydata argument as well as
!   some local type-dependent parameters. The procedures rely heavily on the
!   meta_set_check and meta_set_store functions which do most of the work.
!
! Subroutines:
!
!   meta_set_integer_byname    meta_set_integer_arr_byname
!   meta_set_real4_byname      meta_set_real4_arr_byname
!   meta_set_real8_byname      meta_set_real8_arr_byname
!   meta_set_string_byname     meta_set_string_arr_byname
!   meta_set_logical_byname    meta_set_logical_arr_byname
!   meta_set_pointer2_byname  
!   meta_set_pointer3_byname
!   meta_set_pointer4_byname
!   meta_set_integer_byindex   meta_set_integer_arr_byindex
!   meta_set_real4_byindex     meta_set_real4_arr_byindex
!   meta_set_real8_byindex     meta_set_real8_arr_byindex
!   meta_set_string_byindex    meta_set_string_arr_byindex
!   meta_set_logical_byindex   meta_set_logical_arr_byindex
!   meta_set_pointer2_byindex
!   meta_set_pointer3_byindex
!   meta_set_pointer4_byindex
!
! Error numbers: 1200-1299
!
! Method:
!   Check set request
!   Put data into metadata storage
!
!==============================================================================

SUBROUTINE meta_set_integer_byname(tm,ibuf,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: ibuf
  CHARACTER(LEN=*), INTENT(IN)    :: yname
  INTEGER         , INTENT(IN)    :: ydata
  INTEGER :: iidx, izkey, izsize, izcurpos
  iidx = I_PUTBYNAME
  izsize = SIZE(TRANSFER(ydata, tm%ybuf))
  izcurpos = meta_set_check(tm, ibuf, yname, iidx, I_DATA_TYPE_INTEGER, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_set_store(tm, ibuf, izkey, izcurpos, izsize, TRANSFER(ydata, tm%ybuf))
END SUBROUTINE meta_set_integer_byname

SUBROUTINE meta_set_real4_byname(tm,ibuf,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata)   , INTENT(INOUT) :: tm
  INTEGER            , INTENT(IN)    :: ibuf
  CHARACTER(LEN=*)   , INTENT(IN)    :: yname
  REAL(KIND=sp)      , INTENT(IN)    :: ydata
  INTEGER :: iidx, izkey, izsize, izcurpos
  iidx = I_PUTBYNAME
  izsize = SIZE(TRANSFER(ydata, tm%ybuf))
  izcurpos = meta_set_check(tm, ibuf, yname, iidx, I_DATA_TYPE_REAL, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_set_store(tm, ibuf, izkey, izcurpos, izsize, TRANSFER(ydata, tm%ybuf))
END SUBROUTINE meta_set_real4_byname

SUBROUTINE meta_set_real8_byname(tm,ibuf,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: ibuf
  CHARACTER(LEN=*), INTENT(IN)    :: yname
  REAL(KIND=dp)   , INTENT(IN)    :: ydata
  INTEGER :: iidx, izkey, izsize, izcurpos
  iidx = I_PUTBYNAME
  izsize = SIZE(TRANSFER(ydata, tm%ybuf))
  izcurpos = meta_set_check(tm, ibuf, yname, iidx, I_DATA_TYPE_DOUBLE, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_set_store(tm, ibuf, izkey, izcurpos, izsize, TRANSFER(ydata, tm%ybuf))
END SUBROUTINE meta_set_real8_byname

SUBROUTINE meta_set_string_byname(tm,ibuf,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: ibuf
  CHARACTER(LEN=*), INTENT(IN)    :: yname
  CHARACTER(LEN=*), INTENT(IN)    :: ydata
  CHARACTER(LEN=I_MAX_STRLEN)     :: yzstr
  INTEGER :: iidx, izkey, izsize, izcurpos
  iidx = I_PUTBYNAME
  izsize = SIZE(TRANSFER(TRIM(ydata), tm%ybuf))
  izcurpos = meta_set_check(tm, ibuf, yname, iidx, I_DATA_TYPE_STRING, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  yzstr = TRIM(ydata)
  CALL meta_set_store(tm, ibuf, izkey, izcurpos, izsize, TRANSFER(yzstr, tm%ybuf))
END SUBROUTINE meta_set_string_byname

SUBROUTINE meta_set_logical_byname(tm,ibuf,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: ibuf
  CHARACTER(LEN=*), INTENT(IN)    :: yname
  LOGICAL         , INTENT(IN)    :: ydata
  INTEGER :: iidx, izkey, izsize, izcurpos
  iidx = I_PUTBYNAME
  izsize = SIZE(TRANSFER(ydata, tm%ybuf))
  izcurpos = meta_set_check(tm, ibuf, yname, iidx, I_DATA_TYPE_LOGICAL, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_set_store(tm, ibuf, izkey, izcurpos, izsize, TRANSFER(ydata, tm%ybuf))
END SUBROUTINE meta_set_logical_byname

SUBROUTINE meta_set_pointer2_byname(tm,ibuf,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata)      , INTENT(INOUT) :: tm
  INTEGER               , INTENT(IN)    :: ibuf
  CHARACTER(LEN=*)      , INTENT(IN)    :: yname
  REAL(KIND=wp), POINTER, INTENT(IN)    :: ydata(:,:)
  INTEGER :: iidx, izkey, izsize, izcurpos
  TYPE(t_fld_pointer) :: zp
  imodule_errstat = 0_iintegers
  IF (I_PTR_MAXDIM < 2) THEN
    imodule_errstat = 1200
    ymodule_errstr  = 'meta_set_pointer_byname: maximal number of dimensions is not sufficient (<2)'
    RETURN
  ENDIF
  IF (ASSOCIATED(ydata)) THEN
    IF (ANY(lbound(ydata) /= 1_iintegers)) THEN
      imodule_errstat = 1201
      ymodule_errstr  = 'meta_set_pointer_byname:  some pointer lower bound are not = 0'
      RETURN
    ENDIF
    zp%fld_cptr = C_LOC(ydata(1,1))
    zp%fld_rank = 2_iintegers
    zp%fld_shape(:) = -1_iintegers
    zp%fld_shape(1:2) = ubound(ydata) - lbound(ydata) + 1_iintegers
  ELSE
    zp%fld_cptr = C_NULL_PTR
  END IF
  iidx = I_PUTBYNAME
  izsize = SIZE(TRANSFER(zp, tm%ybuf))
  izcurpos = meta_set_check(tm, ibuf, yname, iidx, I_DATA_TYPE_POINTER, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_set_store(tm, ibuf, izkey, izcurpos, izsize, TRANSFER(zp, tm%ybuf))
END SUBROUTINE meta_set_pointer2_byname

SUBROUTINE meta_set_pointer3_byname(tm,ibuf,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata)      , INTENT(INOUT) :: tm
  INTEGER               , INTENT(IN)    :: ibuf
  CHARACTER(LEN=*)      , INTENT(IN)    :: yname
  REAL(KIND=wp), POINTER, INTENT(IN)    :: ydata(:,:,:)
  INTEGER :: iidx, izkey, izsize, izcurpos
  TYPE(t_fld_pointer) :: zp
  imodule_errstat = 0_iintegers
  IF (I_PTR_MAXDIM < 3) THEN
    imodule_errstat = 1205
    ymodule_errstr  = 'meta_set_pointer_byname: maximal number of dimensions is not sufficient (<3)'
    RETURN
  ENDIF
  IF (ASSOCIATED(ydata)) THEN
    IF (ANY(lbound(ydata) /= 1_iintegers)) THEN
      imodule_errstat = 1206
      ymodule_errstr  = 'meta_set_pointer_byname:  some pointer lower bound are not = 0'
      RETURN
    ENDIF
    zp%fld_cptr = C_LOC(ydata(1,1,1))
    zp%fld_rank = 3_iintegers
    zp%fld_shape(:) = -1_iintegers
    zp%fld_shape(1:3) = ubound(ydata) - lbound(ydata) + 1_iintegers
  ELSE
    zp%fld_cptr = C_NULL_PTR
  END IF
  iidx = I_PUTBYNAME
  izsize = SIZE(TRANSFER(zp, tm%ybuf))
  izcurpos = meta_set_check(tm, ibuf, yname, iidx, I_DATA_TYPE_POINTER, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_set_store(tm, ibuf, izkey, izcurpos, izsize, TRANSFER(zp, tm%ybuf))
END SUBROUTINE meta_set_pointer3_byname

SUBROUTINE meta_set_pointer4_byname(tm,ibuf,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata)      , INTENT(INOUT) :: tm
  INTEGER               , INTENT(IN)    :: ibuf
  CHARACTER(LEN=*)      , INTENT(IN)    :: yname
  REAL(KIND=wp), POINTER, INTENT(IN)    :: ydata(:,:,:,:)
  INTEGER :: iidx, izkey, izsize, izcurpos
  TYPE(t_fld_pointer) :: zp
  imodule_errstat = 0_iintegers
  IF (I_PTR_MAXDIM < 4) THEN
    imodule_errstat = 1210
    ymodule_errstr  = 'meta_set_pointer_byname: maximal number of dimensions is not sufficient (<4)'
    RETURN
  ENDIF
  IF (ASSOCIATED(ydata)) THEN
    IF (ANY(lbound(ydata) /= 1_iintegers)) THEN
      imodule_errstat = 1211
      ymodule_errstr  = 'meta_set_pointer_byname:  some pointer lower bound are not = 0'
      RETURN
    ENDIF
    zp%fld_cptr = C_LOC(ydata(1,1,1,1))
    zp%fld_rank = 4_iintegers
    zp%fld_shape(:) = -1_iintegers
    zp%fld_shape(1:4) = ubound(ydata) - lbound(ydata) + 1_iintegers
  ELSE
    zp%fld_cptr = C_NULL_PTR
  END IF
  iidx = I_PUTBYNAME
  izsize = SIZE(TRANSFER(zp, tm%ybuf))
  izcurpos = meta_set_check(tm, ibuf, yname, iidx, I_DATA_TYPE_POINTER, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_set_store(tm, ibuf, izkey, izcurpos, izsize, TRANSFER(zp, tm%ybuf))
END SUBROUTINE meta_set_pointer4_byname

SUBROUTINE meta_set_integer_arr_byname(tm,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  CHARACTER(LEN=*), INTENT(IN)    :: yname
  INTEGER         , INTENT(IN)    :: ydata(:)
  INTEGER :: iidx, izkey, izsize, izcurpos
  iidx = I_PUTBYNAME
  izsize = SIZE(TRANSFER(ydata, tm%ybuf))
  izcurpos = meta_set_check(tm, I_ALL_BUF, yname, iidx, I_DATA_TYPE_INTEGER, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_set_store(tm, I_ALL_BUF, izkey, izcurpos, izsize, TRANSFER(ydata, tm%ybuf))
END SUBROUTINE meta_set_integer_arr_byname

SUBROUTINE meta_set_real4_arr_byname(tm,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  CHARACTER(LEN=*), INTENT(IN)    :: yname
  REAL(KIND=sp)   , INTENT(IN)    :: ydata(:)
  INTEGER :: iidx, izkey, izsize, izcurpos
  iidx = I_PUTBYNAME
  izsize = SIZE(TRANSFER(ydata, tm%ybuf))
  izcurpos = meta_set_check(tm, I_ALL_BUF, yname, iidx, I_DATA_TYPE_REAL, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_set_store(tm, I_ALL_BUF, izkey, izcurpos, izsize, TRANSFER(ydata, tm%ybuf))
END SUBROUTINE meta_set_real4_arr_byname

SUBROUTINE meta_set_real8_arr_byname(tm,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  CHARACTER(LEN=*), INTENT(IN)    :: yname
  REAL(KIND=dp)   , INTENT(IN)    :: ydata(:)
  INTEGER :: iidx, izkey, izsize, izcurpos
  iidx = I_PUTBYNAME
  izsize = SIZE(TRANSFER(ydata, tm%ybuf))
  izcurpos = meta_set_check(tm, I_ALL_BUF, yname, iidx, I_DATA_TYPE_DOUBLE, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_set_store(tm, I_ALL_BUF, izkey, izcurpos, izsize, TRANSFER(ydata, tm%ybuf))
END SUBROUTINE meta_set_real8_arr_byname

SUBROUTINE meta_set_string_arr_byname(tm,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  CHARACTER(LEN=*), INTENT(IN) :: yname
  CHARACTER(LEN=*), INTENT(IN) :: ydata(:)
  CHARACTER(LEN=I_MAX_STRLEN) :: yzstr(tm%inumbuf)
  INTEGER :: iidx, izkey, izsize, izcurpos, iz
  iidx = I_PUTBYNAME
  izsize = SIZE(TRANSFER(ydata, tm%ybuf))
  izcurpos = meta_set_check(tm, I_ALL_BUF, yname, iidx, I_DATA_TYPE_STRING, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  DO iz = 1, tm%inumbuf
    yzstr(iz) = TRIM(ydata(iz))
  END DO
  CALL meta_set_store(tm, I_ALL_BUF, izkey, izcurpos, izsize, TRANSFER(yzstr, tm%ybuf))
END SUBROUTINE meta_set_string_arr_byname

SUBROUTINE meta_set_logical_arr_byname(tm,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  CHARACTER(LEN=*), INTENT(IN)    :: yname
  LOGICAL         , INTENT(IN)    :: ydata(:)
  INTEGER :: iidx, izkey, izsize, izcurpos
  iidx = I_PUTBYNAME
  izsize = SIZE(TRANSFER(ydata, tm%ybuf))
  izcurpos = meta_set_check(tm, I_ALL_BUF, yname, iidx, I_DATA_TYPE_LOGICAL, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_set_store(tm, I_ALL_BUF, izkey, izcurpos, izsize, TRANSFER(ydata, tm%ybuf))
END SUBROUTINE meta_set_logical_arr_byname

SUBROUTINE meta_set_integer_byindex(tm,ibuf,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: ibuf
  INTEGER         , INTENT(IN)    :: iidx
  INTEGER         , INTENT(IN)    :: ydata
  INTEGER :: izkey, izsize, izcurpos
  izsize = SIZE(TRANSFER(ydata, tm%ybuf))
  izcurpos = meta_set_check(tm, ibuf, '', iidx, I_DATA_TYPE_INTEGER, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_set_store(tm, ibuf, izkey, izcurpos, izsize, TRANSFER(ydata, tm%ybuf))
END SUBROUTINE meta_set_integer_byindex

SUBROUTINE meta_set_real4_byindex(tm,ibuf,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: ibuf
  INTEGER         , INTENT(IN)    :: iidx
  REAL(KIND=sp)   , INTENT(IN)    :: ydata
  INTEGER :: izkey, izsize, izcurpos
  izsize = SIZE(TRANSFER(ydata, tm%ybuf))
  izcurpos = meta_set_check(tm, ibuf, '', iidx, I_DATA_TYPE_REAL, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_set_store(tm, ibuf, izkey, izcurpos, izsize, TRANSFER(ydata, tm%ybuf))
END SUBROUTINE meta_set_real4_byindex

SUBROUTINE meta_set_real8_byindex(tm,ibuf,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: ibuf
  INTEGER         , INTENT(IN)    :: iidx
  REAL(KIND=dp)   , INTENT(IN)    :: ydata
  INTEGER :: izkey, izsize, izcurpos
  izsize = SIZE(TRANSFER(ydata, tm%ybuf))
  izcurpos = meta_set_check(tm, ibuf, '', iidx, I_DATA_TYPE_DOUBLE, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_set_store(tm, ibuf, izkey, izcurpos, izsize, TRANSFER(ydata, tm%ybuf))
END SUBROUTINE meta_set_real8_byindex

SUBROUTINE meta_set_string_byindex(tm,ibuf,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: ibuf
  INTEGER         , INTENT(IN)    :: iidx
  CHARACTER(LEN=*), INTENT(IN)    :: ydata
  CHARACTER(LEN=I_MAX_STRLEN)     :: yzstr
  INTEGER :: izkey, izsize, izcurpos
  izsize = SIZE(TRANSFER(TRIM(ydata), tm%ybuf))
  izcurpos = meta_set_check(tm, ibuf, '', iidx, I_DATA_TYPE_STRING, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  yzstr = TRIM(ydata)
  CALL meta_set_store(tm, ibuf, izkey, izcurpos, izsize, TRANSFER(yzstr, tm%ybuf))
END SUBROUTINE meta_set_string_byindex

SUBROUTINE meta_set_logical_byindex(tm,ibuf,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: ibuf
  INTEGER         , INTENT(IN)    :: iidx
  LOGICAL         , INTENT(IN)    :: ydata
  INTEGER :: izkey, izsize, izcurpos
  izsize = SIZE(TRANSFER(ydata, tm%ybuf))
  izcurpos = meta_set_check(tm, ibuf, '', iidx, I_DATA_TYPE_LOGICAL, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_set_store(tm, ibuf, izkey, izcurpos, izsize, TRANSFER(ydata, tm%ybuf))
END SUBROUTINE meta_set_logical_byindex

SUBROUTINE meta_set_pointer2_byindex(tm,ibuf,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata)      , INTENT(INOUT) :: tm
  INTEGER               , INTENT(IN)    :: ibuf
  INTEGER               , INTENT(IN)    :: iidx
  REAL(KIND=wp), POINTER, INTENT(IN)    :: ydata(:,:)
  INTEGER :: izkey, izsize, izcurpos
  TYPE(t_fld_pointer) :: zp
  imodule_errstat = 0_iintegers
  IF (I_PTR_MAXDIM < 2) THEN
    imodule_errstat = 1220
    ymodule_errstr  = 'meta_set_pointer_byindex: maximal number of dimensions is not sufficient (<2)'
    RETURN
  ENDIF
  IF (ASSOCIATED(ydata)) THEN
    IF (ANY(lbound(ydata) /= 1_iintegers)) THEN
      imodule_errstat = 1221
      ymodule_errstr  = 'meta_set_pointer_byindex:  some pointer lower bound are not = 0'
      RETURN
    ENDIF
    zp%fld_cptr = C_LOC(ydata(1,1))
    zp%fld_rank = 2_iintegers
    zp%fld_shape(:) = -1_iintegers
    zp%fld_shape(1:2) = ubound(ydata) - lbound(ydata) + 1_iintegers
  ELSE
    zp%fld_cptr = C_NULL_PTR
  END IF
  izsize = SIZE(TRANSFER(zp, tm%ybuf))
  izcurpos = meta_set_check(tm, ibuf, '', iidx, I_DATA_TYPE_POINTER, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_set_store(tm, ibuf, izkey, izcurpos, izsize, TRANSFER(zp, tm%ybuf))
END SUBROUTINE meta_set_pointer2_byindex

SUBROUTINE meta_set_pointer3_byindex(tm,ibuf,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata)      , INTENT(INOUT) :: tm
  INTEGER               , INTENT(IN)    :: ibuf
  INTEGER               , INTENT(IN)    :: iidx
  REAL(KIND=wp), POINTER, INTENT(IN)    :: ydata(:,:,:)
  INTEGER :: izkey, izsize, izcurpos
  TYPE(t_fld_pointer) :: zp
  imodule_errstat = 0_iintegers
  IF (I_PTR_MAXDIM < 3) THEN
    imodule_errstat = 1225
    ymodule_errstr  = 'meta_set_pointer_byindex: maximal number of dimensions is not sufficient (<3)'
    RETURN
  ENDIF
  IF (ASSOCIATED(ydata)) THEN
    IF (ANY(lbound(ydata) /= 1_iintegers)) THEN
      imodule_errstat = 1226
      ymodule_errstr  = 'meta_set_pointer_byindex:  some pointer lower bound are not = 0'
      RETURN
    ENDIF
    zp%fld_cptr = C_LOC(ydata(1,1,1))
    zp%fld_rank = 3_iintegers
    zp%fld_shape(:) = -1_iintegers
    zp%fld_shape(1:3) = ubound(ydata) - lbound(ydata) + 1_iintegers
  ELSE
    zp%fld_cptr = C_NULL_PTR
  END IF
  izsize = SIZE(TRANSFER(zp, tm%ybuf))
  izcurpos = meta_set_check(tm, ibuf, '', iidx, I_DATA_TYPE_POINTER, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_set_store(tm, ibuf, izkey, izcurpos, izsize, TRANSFER(zp, tm%ybuf))
END SUBROUTINE meta_set_pointer3_byindex

SUBROUTINE meta_set_pointer4_byindex(tm,ibuf,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata)      , INTENT(INOUT) :: tm
  INTEGER               , INTENT(IN)    :: ibuf
  INTEGER               , INTENT(IN)    :: iidx
  REAL(KIND=wp), POINTER, INTENT(IN)    :: ydata(:,:,:,:)
  INTEGER :: izkey, izsize, izcurpos
  TYPE(t_fld_pointer) :: zp
  imodule_errstat = 0_iintegers
  IF (I_PTR_MAXDIM < 4) THEN
    imodule_errstat = 1230
    ymodule_errstr  = 'meta_set_pointer_byindex: maximal number of dimensions is not sufficient (<4)'
    RETURN
  ENDIF
  IF (ASSOCIATED(ydata)) THEN
    IF (ANY(lbound(ydata) /= 1_iintegers)) THEN
      imodule_errstat = 1231
      ymodule_errstr  = 'meta_set_pointer_byindex:  some pointer lower bound are not = 0'
      RETURN
    ENDIF
    zp%fld_cptr = C_LOC(ydata(1,1,1,1))
    zp%fld_rank = 4_iintegers
    zp%fld_shape(:) = -1_iintegers
    zp%fld_shape(1:4) = ubound(ydata) - lbound(ydata) + 1_iintegers
  ELSE
    zp%fld_cptr = C_NULL_PTR
  END IF
  izsize = SIZE(TRANSFER(zp, tm%ybuf))
  izcurpos = meta_set_check(tm, ibuf, '', iidx, I_DATA_TYPE_POINTER, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_set_store(tm, ibuf, izkey, izcurpos, izsize, TRANSFER(zp, tm%ybuf))
END SUBROUTINE meta_set_pointer4_byindex

SUBROUTINE meta_set_integer_arr_byindex(tm,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: iidx
  INTEGER         , INTENT(IN)    :: ydata(:)
  INTEGER :: izkey, izsize, izcurpos
  izsize = SIZE(TRANSFER(ydata, tm%ybuf))
  izcurpos = meta_set_check(tm, I_ALL_BUF, '', iidx, I_DATA_TYPE_INTEGER, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_set_store(tm, I_ALL_BUF, izkey, izcurpos, izsize, TRANSFER(ydata, tm%ybuf))
END SUBROUTINE meta_set_integer_arr_byindex

SUBROUTINE meta_set_real4_arr_byindex(tm,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: iidx
  REAL(KIND=sp)   , INTENT(IN)    :: ydata(:)
  INTEGER :: izkey, izsize, izcurpos
  izsize = SIZE(TRANSFER(ydata, tm%ybuf))
  izcurpos = meta_set_check(tm, I_ALL_BUF, '', iidx, I_DATA_TYPE_REAL, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_set_store(tm, I_ALL_BUF, izkey, izcurpos, izsize, TRANSFER(ydata, tm%ybuf))
END SUBROUTINE meta_set_real4_arr_byindex

SUBROUTINE meta_set_real8_arr_byindex(tm,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: iidx
  REAL(KIND=dp)   , INTENT(IN)    :: ydata(:)
  INTEGER :: izkey, izsize, izcurpos
  izsize = SIZE(TRANSFER(ydata, tm%ybuf))
  izcurpos = meta_set_check(tm, I_ALL_BUF, '', iidx, I_DATA_TYPE_DOUBLE, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_set_store(tm, I_ALL_BUF, izkey, izcurpos, izsize, TRANSFER(ydata, tm%ybuf))
END SUBROUTINE meta_set_real8_arr_byindex

SUBROUTINE meta_set_string_arr_byindex(tm,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: iidx
  CHARACTER(LEN=*), INTENT(IN)    :: ydata(:)
  CHARACTER(LEN=I_MAX_STRLEN)     :: yzstr(tm%inumbuf)
  INTEGER :: izkey, izsize, izcurpos, iz
  izsize = SIZE(TRANSFER(ydata, tm%ybuf))
  izcurpos = meta_set_check(tm, I_ALL_BUF, '', iidx, I_DATA_TYPE_STRING, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  DO iz = 1, tm%inumbuf
    yzstr(iz) = TRIM(ydata(iz))
  END DO
  CALL meta_set_store(tm, I_ALL_BUF, izkey, izcurpos, izsize, TRANSFER(yzstr, tm%ybuf))
END SUBROUTINE meta_set_string_arr_byindex

SUBROUTINE meta_set_logical_arr_byindex(tm,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: iidx
  LOGICAL         , INTENT(IN)    :: ydata(:)
  INTEGER :: izkey, izsize, izcurpos
  izsize = SIZE(TRANSFER(ydata, tm%ybuf))
  izcurpos = meta_set_check(tm, I_ALL_BUF, '', iidx, I_DATA_TYPE_LOGICAL, izsize, izkey)
  IF (imodule_errstat /= 0_iintegers) RETURN
  CALL meta_set_store(tm, I_ALL_BUF, izkey, izcurpos, izsize, TRANSFER(ydata, tm%ybuf))
END SUBROUTINE meta_set_logical_arr_byindex


!==============================================================================
!+ Module function to check data before storing it into metadata storage
!------------------------------------------------------------------------------

FUNCTION meta_set_check(tm,ibuf,yname,iidx,itype,isize,ikey)
!------------------------------------------------------------------------------
!
! Description:
!   This private function checks a user request to store a metadata item
!   in a metadata storage. It executes several consistency checks and if
!   successfull returns the position in the buffer where to store the data
!   item.
!
! Error numbers: 600-699
!
! Method:
!   Check store request
!
!==============================================================================

IMPLICIT NONE

! Argument list
TYPE(t_metadata), INTENT(IN)    :: tm
INTEGER         , INTENT(IN)    :: ibuf
CHARACTER(LEN=*), INTENT(IN)    :: yname
INTEGER         , INTENT(IN)    :: iidx
INTEGER         , INTENT(IN)    :: itype
INTEGER         , INTENT(INOUT) :: isize
INTEGER         , INTENT(OUT)   :: ikey

! Return value
INTEGER :: meta_set_check

! Local variables
INTEGER :: izcurpos, izcurkey, izbufcount
CHARACTER(LEN=I_MAX_KEYLEN+7) :: yzwhere

!==============================================================================

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin FUNCTION meta_set_check
!------------------------------------------------------------------------------

imodule_errstat = 0_iintegers
IF (iidx == I_PUTBYNAME) THEN
  yzwhere = ' (key=' // TRIM(yname) // ')'
ELSE
  WRITE(yzwhere,*) iidx
  yzwhere = ' (key=' // TRIM(ADJUSTL(yzwhere)) // ')'
END IF

!------------------------------------------------------------------------------
!  Section 1: Check store request
!------------------------------------------------------------------------------

! check if metadata storage is ready for use
IF (.NOT. tm%lisready) THEN
  imodule_errstat = 600_iintegers
  ymodule_errstr = 'metadata_set: metadata storage is not initialized' // yzwhere
  RETURN
END IF

! check for empty arrays
IF (isize == 0_iintegers) THEN
  imodule_errstat = 630_iintegers
  ymodule_errstr = 'metadata_set: cannot store zero sized arrays' // yzwhere
  RETURN
END IF

! check buffer number
IF (ibuf == I_ALL_BUF) THEN
  izbufcount = tm%inumbuf
ELSE
  izbufcount = 1_iintegers
  IF (ibuf <= 0_iintegers .OR. ibuf > tm%inumbuf) THEN
    imodule_errstat = 652_iintegers
    ymodule_errstr = 'metadata_set: requested buffer out of bounds' // yzwhere
    RETURN
  END IF
END IF

! check data size (strings must simply be a multiple of character size)
IF (itype /= I_DATA_TYPE_STRING) THEN
  IF (isize /= izbufcount*I_DATA_TYPE_SIZE(itype)) THEN
    imodule_errstat = 640_iintegers
    ymodule_errstr = 'metadata_set: size mismatch' // yzwhere
    RETURN
  END IF
  isize = isize/izbufcount
ELSE
  IF (MOD(isize,I_DATA_TYPE_SIZE(itype)) /= 0_iintegers) THEN
    imodule_errstat = 641_iintegers
    ymodule_errstr = 'metadata_set: size mismatch' // yzwhere
    RETURN
  END IF
  IF (isize > izbufcount*I_MAX_STRLEN*I_DATA_TYPE_SIZE(itype)) THEN
    imodule_errstat = 642_iintegers
    ymodule_errstr = 'metadata_set: string length exceeds maximum (I_MAX_STRLEN)' &
       // yzwhere
    RETURN
  END IF
  isize = I_MAX_STRLEN*I_DATA_TYPE_SIZE(itype)
END IF

! did user specify key name or index
IF (iidx == I_PUTBYNAME) THEN

  ! check key string length
  IF (LEN_TRIM(yname) > I_MAX_KEYLEN) THEN
    imodule_errstat = 610_iintegers
    ymodule_errstr = 'metadata_set: key string too long' // yzwhere
    RETURN
  END IF
  
  ! find key by name
  DO izcurkey = 1, tm%inumkey
    IF (TRIM(tm%ykey(izcurkey)) == TRIM(yname)) EXIT
  END DO

ELSE

  ! find key by index
  DO izcurkey = 1, tm%inumkey
    IF (tm%iidx(izcurkey) == iidx) EXIT
  END DO

END IF

! if key not found, abort
IF (izcurkey < 1_iintegers .OR. izcurkey > tm%inumkey) THEN
  imodule_errstat = 620_iintegers
  ymodule_errstr = 'metadata_set: key not found or out of range' // yzwhere
  RETURN
END IF

! key exists, extract position
izcurpos = tm%ipos(izcurkey)

! check type
IF (itype /= tm%itype(izcurkey)) THEN
  ! type of existing entry does not match request
  imodule_errstat = 648_iintegers
  ymodule_errstr = 'metadata_set: type mismatch with existing entry' // yzwhere
  RETURN
END IF

! check size
IF (isize /= tm%ilen(izcurkey)) THEN
  ! size of existing entry does not match request
  imodule_errstat = 650_iintegers
  ymodule_errstr = 'metadata_set: size mismatch with existing entry' // yzwhere
  RETURN
END IF

! check if write protected (or undefined)
IF (IAND(tm%iattr(izcurkey),I_STATUS_WRITE) == 0_iintegers) THEN
  ! key exists, but is write protected -> error
  imodule_errstat = 655_iintegers
  ymodule_errstr = 'metadata_set: key is write protected (and defined)' // yzwhere
  RETURN
END IF

! ok, return key number and buffer position
ikey = izcurkey
meta_set_check = izcurpos

!------------------------------------------------------------------------------
!  End of the Function
!------------------------------------------------------------------------------

END FUNCTION meta_set_check

!==============================================================================
!+ Module procedure to store a data item in a metadata storage
!------------------------------------------------------------------------------

SUBROUTINE meta_set_store(tm,ibuf,ikey,ipos,isize,ydata)
!------------------------------------------------------------------------------
!
! Description:
!   This private procedure stores a metadata item which has been checked for
!   consistency into a metadata storage. It requires the storage position and
!   data item size.
!
! Error numbers: 700-799
!
! Method:
!   Store data in buffer
!
!==============================================================================

  IMPLICIT NONE

! Argument list
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER,          INTENT(IN)    :: ibuf
  INTEGER,          INTENT(IN)    :: ikey
  INTEGER,          INTENT(IN)    :: ipos
  INTEGER,          INTENT(IN)    :: isize
  CHARACTER(LEN=1), INTENT(IN)    :: ydata(:)

!==============================================================================

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE meta_set_store
!------------------------------------------------------------------------------

imodule_errstat = 0_iintegers

!------------------------------------------------------------------------------
!  Section 1: Store data in buffer
!------------------------------------------------------------------------------

IF (ibuf == I_ALL_BUF) THEN
  tm%ybuf(ipos:ipos+isize-1,1:tm%inumbuf) = RESHAPE(ydata,(/isize,tm%inumbuf/))
ELSE
  tm%ybuf(ipos:ipos+isize-1,ibuf) = ydata
END IF

!------------------------------------------------------------------------------
!  Section 2: Update defined status
!------------------------------------------------------------------------------

!!!! tm%iattr(ikey) = IOR(tm%iattr(ikey),I_STATUS_DEFINED)

!------------------------------------------------------------------------------
!  End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE meta_set_store


!==============================================================================
!+ Module procedures to get data from metadata storage
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!
! Description:
!   The following procedures overload the procedure metadata_get which is
!   called upon a request to retrieve a data item from the metadata storage. For
!   each data type as well as for scalar/vector data a separate procedure
!   has to be created. The procedures all have exactly the same structure with
!   the exception of dimension and type of the ydata argument as well as
!   some local type-dependent parameters. The procedures rely heavily on the
!   meta_get_check function which do most of the work.
!
! Subroutines:
!
!  meta_get_integer_byname   meta_get_integer_arr_byname
!  meta_get_real4_byname     meta_get_real4_arr_byname
!  meta_get_real8_byname     meta_get_real8_arr_byname
!  meta_get_string_byname    meta_get_string_arr_byname
!  meta_get_logical_byname   meta_get_logical_arr_byname
!  meta_get_pointer2_byname
!  meta_get_pointer3_byname
!  meta_get_pointer4_byname
!  meta_get_integer_byindex  meta_get_integer_arr_byindex
!  meta_get_real4_byindex    meta_get_real4_arr_byindex
!  meta_get_real8_byindex    meta_get_real8_arr_byindex
!  meta_get_string_byindex   meta_get_string_arr_byindex
!  meta_get_logical_byindex  meta_get_logical_arr_byindex
!  meta_get_pointer2_byindex
!  meta_get_pointer3_byindex
!  meta_get_pointer4_byindex
!
! Error numbers: 800-899
!
! Method:
!   Check retrieve request
!   Put data into metadata storage
!
!==============================================================================

SUBROUTINE meta_get_integer_byname(tm,ibuf,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: ibuf
  CHARACTER(LEN=*), INTENT(IN)    :: yname
  INTEGER         , INTENT(OUT)   :: ydata
  INTEGER :: izkey, izcurpos, izcurlen
  izkey = I_GETBYNAME
  izkey = meta_get_check(tm, ibuf, yname, izkey, I_DATA_TYPE_INTEGER, &
               SIZE(TRANSFER(ydata, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  ydata = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,ibuf),ydata)
END SUBROUTINE meta_get_integer_byname

SUBROUTINE meta_get_real4_byname(tm,ibuf,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: ibuf
  CHARACTER(LEN=*), INTENT(IN)    :: yname
  REAL(KIND=sp)   , INTENT(OUT)   :: ydata
  INTEGER :: izkey, izcurpos, izcurlen
  izkey = I_GETBYNAME
  izkey = meta_get_check(tm, ibuf, yname, izkey, I_DATA_TYPE_REAL, &
               SIZE(TRANSFER(ydata, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  ydata = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,ibuf),ydata)
END SUBROUTINE meta_get_real4_byname

SUBROUTINE meta_get_real8_byname(tm,ibuf,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: ibuf
  CHARACTER(LEN=*), INTENT(IN)    :: yname
  REAL(KIND=dp)   , INTENT(OUT)   :: ydata
  INTEGER :: izkey, izcurpos, izcurlen
  izkey = I_GETBYNAME
  izkey = meta_get_check(tm, ibuf, yname, izkey, I_DATA_TYPE_DOUBLE, &
               SIZE(TRANSFER(ydata, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  ydata = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,ibuf),ydata)
END SUBROUTINE meta_get_real8_byname

SUBROUTINE meta_get_string_byname(tm,ibuf,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: ibuf
  CHARACTER(LEN=*), INTENT(IN)    :: yname
  CHARACTER(LEN=*), INTENT(OUT)   :: ydata
  CHARACTER(LEN=I_MAX_STRLEN)     :: yzstr
  INTEGER :: izkey, izcurpos, izcurlen
  izkey = I_GETBYNAME
  izkey = meta_get_check(tm, ibuf, yname, izkey, I_DATA_TYPE_STRING, &
               SIZE(TRANSFER(ydata, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  yzstr = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,ibuf),yzstr(1:izcurlen))
  ydata = TRIM(yzstr)
END SUBROUTINE meta_get_string_byname

SUBROUTINE meta_get_logical_byname(tm,ibuf,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: ibuf
  CHARACTER(LEN=*), INTENT(IN)    :: yname
  LOGICAL         , INTENT(OUT)   :: ydata
  INTEGER :: izkey, izcurpos, izcurlen
  izkey = I_GETBYNAME
  izkey = meta_get_check(tm, ibuf, yname, izkey, I_DATA_TYPE_LOGICAL, &
               SIZE(TRANSFER(ydata, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  ydata = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,ibuf),ydata)
END SUBROUTINE meta_get_logical_byname

SUBROUTINE meta_get_pointer2_byname(tm,ibuf,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata)      , INTENT(INOUT) :: tm
  INTEGER               , INTENT(IN)    :: ibuf
  CHARACTER(LEN=*)      , INTENT(IN)    :: yname
  REAL(KIND=wp), POINTER, INTENT(OUT)   :: ydata(:,:)
  INTEGER :: izkey, izcurpos, izcurlen
  TYPE(t_fld_pointer) :: zp
  imodule_errstat = 0_iintegers
  izkey = I_GETBYNAME
  izkey = meta_get_check(tm, ibuf, yname, izkey, I_DATA_TYPE_POINTER, &
               SIZE(TRANSFER(zp, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  zp = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,ibuf),zp)
  IF (C_ASSOCIATED(zp%fld_cptr)) THEN
    IF (zp%fld_rank /= 2) THEN
      imodule_errstat = 800_iintegers
      ymodule_errstr  = 'meta_get_pointer_byname: wrong pointer rank (/=2)'
      RETURN
    ENDIF
    IF (ANY(zp%fld_shape(1:2) <= 0_iintegers)) THEN
      imodule_errstat = 801_iintegers
      ymodule_errstr  = 'meta_get_pointer_byname: one or several wrong dimension(s) for pointer (<=0)'
      RETURN
    ENDIF
    CALL C_F_POINTER(zp%fld_cptr, ydata, zp%fld_shape(1:2))
  ELSE
    NULLIFY(ydata)
  END IF
END SUBROUTINE meta_get_pointer2_byname

SUBROUTINE meta_get_pointer3_byname(tm,ibuf,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata)      , INTENT(INOUT) :: tm
  INTEGER               , INTENT(IN)    :: ibuf
  CHARACTER(LEN=*)      , INTENT(IN)    :: yname
  REAL(KIND=wp), POINTER, INTENT(OUT)   :: ydata(:,:,:)
  INTEGER :: izkey, izcurpos, izcurlen
  TYPE(t_fld_pointer) :: zp
  imodule_errstat = 0_iintegers
  izkey = I_GETBYNAME
  izkey = meta_get_check(tm, ibuf, yname, izkey, I_DATA_TYPE_POINTER, &
               SIZE(TRANSFER(zp, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  zp = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,ibuf),zp)
  IF (C_ASSOCIATED(zp%fld_cptr)) THEN
    IF (zp%fld_rank /= 3) THEN
      imodule_errstat = 805_iintegers
      ymodule_errstr  = 'meta_get_pointer_byname: wrong pointer rank (/=3)'
      RETURN
    ENDIF
    IF (ANY(zp%fld_shape(1:3) <= 0_iintegers)) THEN
      imodule_errstat = 806_iintegers
      ymodule_errstr  = 'meta_get_pointer_byname: one or several wrong dimension(s) for pointer (<=0)'
      RETURN
    ENDIF
    CALL C_F_POINTER(zp%fld_cptr, ydata, zp%fld_shape(1:3))
  ELSE
    NULLIFY(ydata)
  END IF
END SUBROUTINE meta_get_pointer3_byname

SUBROUTINE meta_get_pointer4_byname(tm,ibuf,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata)      , INTENT(INOUT) :: tm
  INTEGER               , INTENT(IN)    :: ibuf
  CHARACTER(LEN=*)      , INTENT(IN)    :: yname
  REAL(KIND=wp), POINTER, INTENT(OUT)   :: ydata(:,:,:,:)
  INTEGER :: izkey, izcurpos, izcurlen
  TYPE(t_fld_pointer) :: zp
  imodule_errstat = 0_iintegers
  izkey = I_GETBYNAME
  izkey = meta_get_check(tm, ibuf, yname, izkey, I_DATA_TYPE_POINTER, &
               SIZE(TRANSFER(zp, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  zp = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,ibuf),zp)
  IF (C_ASSOCIATED(zp%fld_cptr)) THEN
    IF (zp%fld_rank /= 4) THEN
      imodule_errstat = 810_iintegers
      ymodule_errstr  = 'meta_get_pointer_byname: wrong pointer rank (/=4)'
      RETURN
    ENDIF
    IF (ANY(zp%fld_shape(1:4) <= 0_iintegers)) THEN
      imodule_errstat = 811_iintegers
      ymodule_errstr  = 'meta_get_pointer_byname: one or several wrong dimension(s) for pointer (<=0)'
      RETURN
    ENDIF
    CALL C_F_POINTER(zp%fld_cptr, ydata, zp%fld_shape(1:4))
  ELSE
    NULLIFY(ydata)
  END IF
END SUBROUTINE meta_get_pointer4_byname

SUBROUTINE meta_get_integer_arr_byname(tm,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  CHARACTER(LEN=*), INTENT(IN)    :: yname
  INTEGER         , INTENT(OUT)   :: ydata(:)
  INTEGER :: izkey, izcurpos, izcurlen
  izkey = I_GETBYNAME
  izkey = meta_get_check(tm, I_ALL_BUF, yname, izkey, I_DATA_TYPE_INTEGER, &
               SIZE(TRANSFER(ydata, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  ydata(:) = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,1:tm%inumbuf),ydata)
END SUBROUTINE meta_get_integer_arr_byname

SUBROUTINE meta_get_real4_arr_byname(tm,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  CHARACTER(LEN=*), INTENT(IN)    :: yname
  REAL(KIND=sp)   , INTENT(OUT)   :: ydata(:)
  INTEGER :: izkey, izcurpos, izcurlen
  izkey = I_GETBYNAME
  izkey = meta_get_check(tm, I_ALL_BUF, yname, izkey, I_DATA_TYPE_REAL, &
               SIZE(TRANSFER(ydata, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  ydata(:) = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,1:tm%inumbuf),ydata)
END SUBROUTINE meta_get_real4_arr_byname

SUBROUTINE meta_get_real8_arr_byname(tm,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  CHARACTER(LEN=*), INTENT(IN)    :: yname
  REAL(KIND=dp)   , INTENT(OUT)   :: ydata(:)
  INTEGER :: izkey, izcurpos, izcurlen
  izkey = I_GETBYNAME
  izkey = meta_get_check(tm, I_ALL_BUF, yname, izkey, I_DATA_TYPE_DOUBLE, &
               SIZE(TRANSFER(ydata, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  ydata(:) = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,1:tm%inumbuf),ydata)
END SUBROUTINE meta_get_real8_arr_byname

SUBROUTINE meta_get_string_arr_byname(tm,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  CHARACTER(LEN=*), INTENT(IN) :: yname
  CHARACTER(LEN=*), INTENT(OUT) :: ydata(:)
  CHARACTER(LEN=I_MAX_STRLEN) :: yzstr(tm%inumbuf)
  INTEGER :: izkey, izcurpos, izcurlen, iz
  izkey = I_GETBYNAME
  izkey = meta_get_check(tm, I_ALL_BUF, yname, izkey, I_DATA_TYPE_STRING, &
               SIZE(TRANSFER(ydata, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  yzstr(:) = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,1:tm%inumbuf),yzstr)
  DO iz = 1, tm%inumbuf
    ydata(iz) = TRIM(yzstr(iz))
  END DO
END SUBROUTINE meta_get_string_arr_byname

SUBROUTINE meta_get_logical_arr_byname(tm,yname,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  CHARACTER(LEN=*), INTENT(IN)    :: yname
  LOGICAL         , INTENT(OUT)   :: ydata(:)
  INTEGER :: izkey, izcurpos, izcurlen
  izkey = I_GETBYNAME
  izkey = meta_get_check(tm, I_ALL_BUF, yname, izkey, I_DATA_TYPE_LOGICAL, &
               SIZE(TRANSFER(ydata, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  ydata(:) = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,1:tm%inumbuf),ydata)
END SUBROUTINE meta_get_logical_arr_byname

SUBROUTINE meta_get_integer_byindex(tm,ibuf,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: ibuf
  INTEGER         , INTENT(IN)    :: iidx
  INTEGER         , INTENT(OUT)   :: ydata
  INTEGER :: izkey, izcurpos, izcurlen
  izkey = iidx
  izkey = meta_get_check(tm, ibuf, '', izkey, I_DATA_TYPE_INTEGER, &
               SIZE(TRANSFER(ydata, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  ydata = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,ibuf),ydata)
END SUBROUTINE meta_get_integer_byindex

SUBROUTINE meta_get_real4_byindex(tm,ibuf,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: ibuf
  INTEGER         , INTENT(IN)    :: iidx
  REAL(KIND=sp)   , INTENT(OUT)   :: ydata
  INTEGER :: izkey, izcurpos, izcurlen
  izkey = iidx
  izkey = meta_get_check(tm, ibuf, '', izkey, I_DATA_TYPE_REAL, &
               SIZE(TRANSFER(ydata, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  ydata = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,ibuf),ydata)
END SUBROUTINE meta_get_real4_byindex

SUBROUTINE meta_get_real8_byindex(tm,ibuf,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: ibuf
  INTEGER         , INTENT(IN)    :: iidx
  REAL(KIND=dp)   , INTENT(OUT)   :: ydata
  INTEGER :: izkey, izcurpos, izcurlen
  izkey = iidx
  izkey = meta_get_check(tm, ibuf, '', izkey, I_DATA_TYPE_DOUBLE, &
               SIZE(TRANSFER(ydata, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  ydata = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,ibuf),ydata)
END SUBROUTINE meta_get_real8_byindex

SUBROUTINE meta_get_string_byindex(tm,ibuf,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: ibuf
  INTEGER         , INTENT(IN)    :: iidx
  CHARACTER(LEN=*), INTENT(OUT)   :: ydata
  CHARACTER(LEN=I_MAX_STRLEN)     :: yzstr
  INTEGER :: izkey, izcurpos, izcurlen
  izkey = iidx
  izkey = meta_get_check(tm, ibuf, '', izkey, I_DATA_TYPE_STRING, &
               SIZE(TRANSFER(ydata, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  yzstr = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,ibuf),yzstr(1:izcurlen))
  ydata = TRIM(yzstr)
END SUBROUTINE meta_get_string_byindex

SUBROUTINE meta_get_logical_byindex(tm,ibuf,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: ibuf
  INTEGER         , INTENT(IN)    :: iidx
  lOGICAL         , INTENT(OUT)   :: ydata
  INTEGER :: izkey, izcurpos, izcurlen
  izkey = iidx
  izkey = meta_get_check(tm, ibuf, '', izkey, I_DATA_TYPE_LOGICAL, &
               SIZE(TRANSFER(ydata, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  ydata = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,ibuf),ydata)
END SUBROUTINE meta_get_logical_byindex

SUBROUTINE meta_get_pointer2_byindex(tm,ibuf,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata)      , INTENT(INOUT) :: tm
  INTEGER               , INTENT(IN)    :: ibuf
  INTEGER               , INTENT(IN)    :: iidx
  REAL(KIND=wp), POINTER, INTENT(OUT)   :: ydata(:,:)
  INTEGER :: izkey, izcurpos, izcurlen
  TYPE(t_fld_pointer) :: zp
  imodule_errstat = 0_iintegers
  izkey = iidx
  izkey = meta_get_check(tm, ibuf, '', izkey, I_DATA_TYPE_POINTER, &
               SIZE(TRANSFER(zp, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  zp = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,ibuf),zp)
  IF (C_ASSOCIATED(zp%fld_cptr)) THEN
    IF (zp%fld_rank /= 2) THEN
      imodule_errstat = 820_iintegers
      ymodule_errstr  = 'meta_get_pointer_byindex: wrong pointer rank (/=2)'
      RETURN
    ENDIF
    IF (ANY(zp%fld_shape(1:2) <= 0_iintegers)) THEN
      imodule_errstat = 821_iintegers
      ymodule_errstr  = 'meta_get_pointer_byname: one or several wrong dimension(s) for pointer (<=0)'
      RETURN
    ENDIF
    CALL C_F_POINTER(zp%fld_cptr, ydata, zp%fld_shape(1:2))
  ELSE
    NULLIFY(ydata)
  END IF
END SUBROUTINE meta_get_pointer2_byindex

SUBROUTINE meta_get_pointer3_byindex(tm,ibuf,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata)      , INTENT(INOUT) :: tm
  INTEGER               , INTENT(IN)    :: ibuf
  INTEGER               , INTENT(IN)    :: iidx
  REAL(KIND=wp), POINTER, INTENT(OUT)   :: ydata(:,:,:)
  INTEGER :: izkey, izcurpos, izcurlen
  TYPE(t_fld_pointer) :: zp
  imodule_errstat = 0_iintegers
  izkey = iidx
  izkey = meta_get_check(tm, ibuf, '', izkey, I_DATA_TYPE_POINTER, &
               SIZE(TRANSFER(zp, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  zp = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,ibuf),zp)
  IF (C_ASSOCIATED(zp%fld_cptr)) THEN
    IF (zp%fld_rank /= 3) THEN
      imodule_errstat = 825_iintegers
      ymodule_errstr  = 'meta_get_pointer_byindex: wrong pointer rank (/=3)'
      RETURN
    ENDIF
    IF (ANY(zp%fld_shape(1:3) <= 0_iintegers)) THEN
      imodule_errstat = 826_iintegers
      ymodule_errstr  = 'meta_get_pointer_byname: one or several wrong dimension(s) for pointer (<=0)'
      RETURN
    ENDIF
    CALL C_F_POINTER(zp%fld_cptr, ydata, zp%fld_shape(1:3))
  ELSE
    NULLIFY(ydata)
  END IF
END SUBROUTINE meta_get_pointer3_byindex

SUBROUTINE meta_get_pointer4_byindex(tm,ibuf,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata)      , INTENT(INOUT) :: tm
  INTEGER               , INTENT(IN)    :: ibuf
  INTEGER               , INTENT(IN)    :: iidx
  REAL(KIND=wp), POINTER, INTENT(OUT)   :: ydata(:,:,:,:)
  INTEGER :: izkey, izcurpos, izcurlen
  TYPE(t_fld_pointer) :: zp
  imodule_errstat = 0_iintegers
  izkey = iidx
  izkey = meta_get_check(tm, ibuf, '', izkey, I_DATA_TYPE_POINTER, &
               SIZE(TRANSFER(zp, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  zp = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,ibuf),zp)
  IF (C_ASSOCIATED(zp%fld_cptr)) THEN
    IF (zp%fld_rank /= 4) THEN
      imodule_errstat = 830_iintegers
      ymodule_errstr  = 'meta_get_pointer_byindex: wrong pointer rank (/=4)'
      RETURN
    ENDIF
    IF (ANY(zp%fld_shape(1:4) <= 0_iintegers)) THEN
      imodule_errstat = 831_iintegers
      ymodule_errstr  = 'meta_get_pointer_byname: one or several wrong dimension(s) for pointer (<=0)'
      RETURN
    ENDIF
    CALL C_F_POINTER(zp%fld_cptr, ydata, zp%fld_shape(1:4))
  ELSE
    NULLIFY(ydata)
  END IF
END SUBROUTINE meta_get_pointer4_byindex

SUBROUTINE meta_get_integer_arr_byindex(tm,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: iidx
  INTEGER         , INTENT(OUT)   :: ydata(:)
  INTEGER :: izkey, izcurpos, izcurlen
  izkey = iidx
  izkey = meta_get_check(tm, I_ALL_BUF, '', izkey, I_DATA_TYPE_INTEGER, &
               SIZE(TRANSFER(ydata, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  ydata(:) = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,1:tm%inumbuf),ydata)
END SUBROUTINE meta_get_integer_arr_byindex

SUBROUTINE meta_get_real4_arr_byindex(tm,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: iidx
  REAL(KIND=sp)   , INTENT(OUT)   :: ydata(:)
  INTEGER :: izkey, izcurpos, izcurlen
  izkey = iidx
  izkey = meta_get_check(tm, I_ALL_BUF, '', izkey, I_DATA_TYPE_REAL, &
               SIZE(TRANSFER(ydata, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  ydata(:) = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,1:tm%inumbuf),ydata)
END SUBROUTINE meta_get_real4_arr_byindex

SUBROUTINE meta_get_real8_arr_byindex(tm,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: iidx
  REAL(KIND=dp)   , INTENT(OUT)   :: ydata(:)
  INTEGER :: izkey, izcurpos, izcurlen
  izkey = iidx
  izkey = meta_get_check(tm, I_ALL_BUF, '', izkey, I_DATA_TYPE_DOUBLE, &
               SIZE(TRANSFER(ydata, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  ydata(:) = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,1:tm%inumbuf),ydata)
END SUBROUTINE meta_get_real8_arr_byindex

SUBROUTINE meta_get_string_arr_byindex(tm,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: iidx
  CHARACTER(LEN=*), INTENT(OUT)   :: ydata(:)
  CHARACTER(LEN=I_MAX_STRLEN)     :: yzstr(tm%inumbuf)
  INTEGER :: izkey, izcurpos, izcurlen, iz
  izkey = iidx
  izkey = meta_get_check(tm, I_ALL_BUF, '', izkey, I_DATA_TYPE_STRING, &
               SIZE(TRANSFER(ydata, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  yzstr(:) = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,1:tm%inumbuf),yzstr)
  DO iz = 1, tm%inumbuf
    ydata(iz) = TRIM(yzstr(iz))
  END DO
END SUBROUTINE meta_get_string_arr_byindex

SUBROUTINE meta_get_logical_arr_byindex(tm,iidx,ydata)
  IMPLICIT NONE
  TYPE(t_metadata), INTENT(INOUT) :: tm
  INTEGER         , INTENT(IN)    :: iidx
  LOGICAL         , INTENT(OUT)   :: ydata(:)
  INTEGER :: izkey, izcurpos, izcurlen
  izkey = iidx
  izkey = meta_get_check(tm, I_ALL_BUF, '', izkey, I_DATA_TYPE_LOGICAL, &
               SIZE(TRANSFER(ydata, tm%ybuf)))
  IF (imodule_errstat /= 0_iintegers) RETURN
  izcurpos = tm%ipos(izkey)
  izcurlen = tm%ilen(izkey)
  ydata(:) = TRANSFER(tm%ybuf(izcurpos:izcurpos+izcurlen-1,1:tm%inumbuf),ydata)
END SUBROUTINE meta_get_logical_arr_byindex


!==============================================================================
!+ Module function to check request before getting it from metadata storage
!------------------------------------------------------------------------------

FUNCTION meta_get_check(tm,ibuf,yname,iidx,itype,isize)
!------------------------------------------------------------------------------
!
! Description:
!   This private fuction checks a user request to retrieve a metadata item
!   in a metadata storage. It executes several consistency checks and if
!   successfull returns the position of the metadata item in the storage buffer.
!
! Error numbers: 900-999
!
! Method:
!   Check retrieve request
!
!==============================================================================

IMPLICIT NONE

! Argument list
TYPE(t_metadata), INTENT(IN) :: tm
INTEGER         , INTENT(IN) :: ibuf
CHARACTER(LEN=*), INTENT(IN) :: yname
INTEGER         , INTENT(IN) :: iidx
INTEGER         , INTENT(IN) :: itype
INTEGER         , INTENT(IN) :: isize

! Return value
INTEGER :: meta_get_check

! Local variables
INTEGER :: izcurkey, izbufcount
CHARACTER(LEN=I_MAX_KEYLEN+7) :: yzwhere

!==============================================================================

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin FUNCTION meta_get_check
!------------------------------------------------------------------------------

imodule_errstat = 0_iintegers
IF (iidx == I_GETBYNAME) THEN
  yzwhere = ' (key=' // TRIM(yname) // ')'
ELSE
  write(yzwhere,*) iidx
  yzwhere = ' (key=' // TRIM(ADJUSTL(yzwhere)) // ')'
END IF

!------------------------------------------------------------------------------
!  Section 1: Check retrieve request
!------------------------------------------------------------------------------

! check if metadata storage is ready for use
IF (.NOT. tm%lisready) THEN
  imodule_errstat = 900_iintegers
  ymodule_errstr = 'metadata_get: metadata storage not initialized' // yzwhere
  RETURN
END IF

! are we looking for by name
IF (iidx == I_GETBYNAME) THEN

  ! check key string length
  IF (LEN_TRIM(yname) > I_MAX_KEYLEN) THEN
    imodule_errstat = 910_iintegers
    ymodule_errstr = 'metadata_get: key string too long' // yzwhere
    RETURN
  END IF
  
  ! find key by name
  DO izcurkey = 1, tm%inumkey
    IF (TRIM(tm%ykey(izcurkey)) == TRIM(yname)) EXIT
  END DO

ELSE

  ! find key by index
  DO izcurkey = 1, tm%inumkey
    IF (tm%iidx(izcurkey) == iidx) EXIT
  END DO
  
END IF

! check if key is found and valid
IF (izcurkey < 1_iintegers .OR. izcurkey > tm%inumkey) THEN
  imodule_errstat = 920_iintegers
  ymodule_errstr = 'metadata_get: specified key not found or out of range' // yzwhere
  RETURN
END IF

! check if type is ok
IF (tm%itype(izcurkey) /= itype) THEN
  imodule_errstat = 930_iintegers
  ymodule_errstr = 'metadata_get: type mismatch' // yzwhere
  RETURN
END IF

! check buffer number
IF (ibuf == I_ALL_BUF) THEN
  izbufcount = tm%inumbuf
ELSE
  izbufcount = 1_iintegers
  IF (ibuf <= 0_iintegers .OR. ibuf > tm%inumbuf) THEN
    imodule_errstat = 952_iintegers
    ymodule_errstr = 'metadata_get: requested buffer out of bounds' // yzwhere
    RETURN
  END IF
END IF

! check if size is ok
IF (itype /= I_DATA_TYPE_STRING) THEN
  IF (tm%ilen(izcurkey)*izbufcount /= isize) THEN
    imodule_errstat = 960_iintegers
    ymodule_errstr = 'metadata_get: size mismatch' // yzwhere
    RETURN
  END IF
!ELSE
!  IF (tm%ilen(izcurkey) > isize) THEN
!    imodule_errstat = 950_iintegers
!    ymodule_errstr = 'metadata_get: size mismatch' // yzwhere
!    RETURN
!  END IF
END IF

! check if read protected (or undefined)
IF (IAND(tm%iattr(izcurkey),I_STATUS_READ) == 0_iintegers) THEN
  ! key exists, but is read protected -> error
  imodule_errstat = 970_iintegers
  ymodule_errstr = 'metadata_get: key is read protected' // yzwhere
  RETURN
END IF

! ok, return element index
meta_get_check = izcurkey

!------------------------------------------------------------------------------
!  End of the Function
!------------------------------------------------------------------------------

END FUNCTION meta_get_check

END MODULE src_tracer_metadata

