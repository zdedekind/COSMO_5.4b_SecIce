!==============================================================================
!+ Module for parallel I/O in meteorological codes
!------------------------------------------------------------------------------

MODULE mpe_io2

!------------------------------------------------------------------------------
!
! Description:
!   This module provides routines and data structures for 
!   conducting parallel I/O in several modes and on different
!   devices (DWD file I/O, DWD database I/O). The design is
!   kept in layers to support the inclusion of new devices.
!
! Note:
!   Some constants are contained in the module interface, which
!   may be redefined for performance/debugging reasons.
!   These are marked by "!<<"
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date          Name
! 1.1        1999/01/21    Heinrich Bockhorst (PALLAS)
!  Initial release
! 1.2        1999/02/26    H. J. Plum / H. Bockhorst (PALLAS)
!  Integration of namelist variables database/idbg_level, database/iretry;
!  stdout tidied up
! 1.3        2004/10/23    Ulrich Schaettler
!  Optimizations; get rid of the database, which was never used
! 1.4        2011/07/21    F. Prill
!  Extended previous work from U. Schaettler: non-blocking
!  communication; buffered pre-fetching of data; changed interfaces
! 1.5        2012/10/25    Ulrich Schaettler, Burkhardt Rockel
!  SR mpe_io_read: Move MPI_Wait in loop over PEs outside the "IF (ilength > 0)"
!   (otherwise the former message could be damaged) (US)
!  Include missing #ifdef GRIBDWD (BR)
! 1.6          2013/04/26  Ulrich Schaettler
!  Implemented grib-api calls in mpe_io_open, mpe_io_read, mpe_io_write, mpe_io_close
!  Re-interpretation of the 3 characters in ymode (only first one used up to now)
!    1st:  'r', 'w' for reading / writing
!    2nd:  for library used:  'g','1','2' for DWD griblib, grib_api/Grib1, grib_api/Grib2
!    3rd:  for speed of writing: 's' or ' ' (slow/sorted), 'f' (fast/unsorted)
!  Inlined SR buffer_list_probe_head into SR buffer_list_write_action
!  SR buffer_list_finalize_write: Removed the last call to SR buffer_list_write_action 
!    (because of problems on IBM pwr7 when turning on array bound checking)
! 1.7          2014/03/26  Ulrich Schaettler, Helmut Frank, Florian Prill
!  Replaced MPI_COMM_WORLD by icomm_ori in SR mpe_io_init, mpe_io_shutdown
!  Added possibility to append GRIB files with grib_api by choosing 'a' 
!    for 1st character in ymode
!  Implemented FIFO queue for incoming messages and limited receive requests
!    by a parameter MAX_RECV_REQUESTS
! 1.8          2015/04/30  Ulrich Schaettler
!  Only allocate write_buffer in IO-PE for writing, not for reading
!
!==============================================================================

! still use include mpif.h, for running sequentially
!USE MPI

#ifdef GRIBAPI
! grib_api interface
USE grib_api
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================

INCLUDE "mpif.h"

!==============================================================================

! public interfaces of mpe_io2

PUBLIC  &
  mpe_io_init                       ,   & ! initialization
  mpe_io_open                       ,   & ! opens a file
  mpe_io_read                       ,   & ! read from file
  mpe_io_write                      ,   & ! writes to file
  mpe_io_complete                   ,   & ! complete a (collective) mpe_io_write
  mpe_io_close                      ,   & ! close file on I/O PEs
  mpe_io_node                       ,   & ! program for "dumb" io node
  mpe_io_shutdown                   ,   & ! send shutdown signal to io_nodes
  mpe_io_ready                      ,   & ! actions at the end of a timestep
  mpe_is_compute                    ,   & ! identifies compute PEs
  mpe_is_io                               ! identifies io      PEs

! everything else is private

PRIVATE

!------------------------------------------------------------------------------

! integer constant determining MPI data type:
INTEGER, PARAMETER :: MPI_DATATYPE = MPI_INTEGER

! KIND -parameters for declaring the variables that are passed to the
! grib library:  values for the DWDLIB
INTEGER, PARAMETER :: gribc_kind = KIND(1)
INTEGER, PARAMETER :: gribf_kind = KIND(1)

! integer precision necessary for grib_api in interfaces where length of
! message in bytes is involved
INTEGER, PARAMETER ::                                         &
#ifdef GRIBAPI
       gapi_int  = kindOfSize              ! should be INTEGER *8 where necessary
#else
       gapi_int  = SELECTED_INT_KIND (8)   ! INTEGER *4
#endif


!------------------------------------------------------------------------------

! Communicator and process identification

INTEGER :: io_config            ! identify io configuration

INTEGER :: MPE_COMM_IO,       & ! communicator for IO processors
  MPE_COMM_COMPUTE,           & ! communicator for compute processors
  MPE_COMM_CO2IO,             & ! communicator between compute and IO procs
  MPE_COMM_CLOSEFILE,         & ! MPI communicator for close-file messages
  MPE_COMM_CONTEXT,           & ! MPI communicators containing only pairs of IO/compute PEs
  MPE_COMM_WRITE,             & ! MPI communicator for grib write
  MPE_COMM_READ,              & ! MPI communicator for grib read
  !
  me_ori,                     & ! rank in original communicator
  num_io,                     & ! number of processors in IO communicator
  me_io,                      & ! rank of a processor in IO communicator
  num_compute,                & ! number of compute processors
  me_compute                    ! rank of a processor in compute communicator

INTEGER, ALLOCATABLE ::  &
  iclose_comm_compute_ranks(:), & ! process rank translation mpe_comm_compute -> MPE_COMM_CLOSEFILE
  iclose_comm_io_ranks(:)         ! process rank translation mpe_comm_io      -> MPE_COMM_CLOSEFILE

!------------------------------------------------------------------------------

! special tags used inside the mpe_io2 system
INTEGER, PARAMETER :: TAG_OPEN      =  6666   ! tag for "end-message" of write 

! flags (stored at position MPE_FLAG of field messages)
INTEGER, PARAMETER :: MPE_FLAG            =  0 ! (= array index)
INTEGER, PARAMETER :: MPE_FLAG_WRITE_END  =  -1
INTEGER, PARAMETER :: MPE_FLAG_READ_END   =  -2

! tags identifying the type of a context message
INTEGER, PARAMETER :: MSG_ID_DATNAME               =  0
INTEGER, PARAMETER :: MSG_ID_STOP                  = -1
INTEGER, PARAMETER :: MSG_ID_WRITE_READY_FILE      = -2
INTEGER, PARAMETER :: MSG_ID_DONT_WRITE_READY_FILE = -3

! Field description 
INTEGER :: field_snd          ! index of the current field to be sent
INTEGER :: field_rcv          ! index of the current field to be read

LOGICAL :: &
  is_compute_pe, is_io_pe

!--- Buffer environment and management ----------------------------------------

INTEGER ::   &
  BUFF_length                   ! value for dwdlib, corresponds to lfd (Dimension for iblock)

INTEGER(KIND=gribf_kind), ALLOCATABLE ::   &
  irecv_buffer(:,:),          & ! buffer for data that is received from compute 
                                ! processors and written to disk
  isend_buffer(:,:),          & ! buffer for data that is send from compute 
                                ! processors to IO processors
  iclose_msg(:,:),            & ! received close-file messages
  iread_recv_buffer(:)          ! buffer for receiving fields from IO PE in read mode

!<< number of maximal possible send requests
!     this may cost memory, but too small values slow down
!     the contributing compute PEs.
INTEGER, PARAMETER                    ::                   &
  MAX_SEND_REQUESTS = 5

INTEGER, PARAMETER                    ::                   &
  MAX_RECV_REQUESTS = 5  

INTEGER                               ::                   &
  isend_current_request,      & ! to identify the current send request
  irecv_current_request,      & ! current request
  iread_recv_request            ! request for IRECV in read mode

INTEGER, ALLOCATABLE                  ::                   &
  irecv_requests(:),          & ! for storing irecv requests
  irecv_isource(:),           & ! source rank for irecv request
  irecv_status(:,:),          & ! for storing the status of the irecv's
  isend_requests(:),          & ! for storing at least two different send 
                                ! requests in the compute processors
  isend_status(:),            & ! for storing the status of the sends
  iclose_rcv_request(:),      & ! identifies received close-file file messages
  iclose_rcv_status(:,:)

!--- Device specific stuff: ---------------------------------------------------

CHARACTER(LEN=1)     :: current_mode            ! 'r'/'w'/'a' : read/write/append mode
CHARACTER(LEN=1)     :: current_speed           ! 's' for slow/sorted
                                                ! 'f' for fast/unsorted
CHARACTER(LEN=4)     :: current_lib             ! for grib-library used
LOGICAL              :: lprefetching_started    ! Flag. True, if prefetching process under way
LOGICAL              :: lprefetching_available  ! Flag. True, if prefetched data is up to date

INTEGER, PARAMETER   :: MAX_NAME_LEN   = 256    !<< maximal length for file names (is 250 in the model)
INTEGER, PARAMETER   :: ISAFETY_MARGIN =   4    !<< no. of additional ints sent/received

!--- data type for context messages (shutdown, ready file, open, ...) ---------

TYPE mpe_context_msg
  INTEGER ::  &
    id,       &                 ! message type: shutdown, ready mode, ...
    length,   &                 ! length of filename
    iunit,    &                 ! unit number for file    
    itag                        ! additional key value (for counting open/close)
  CHARACTER (LEN=MAX_NAME_LEN) :: &
    yfilename                   ! file name
  CHARACTER (LEN=3) ::  &
    ymode                       ! file mode
  CHARACTER (LEN=4) ::  &
    ylib                        ! data format (or grib library used)
END TYPE mpe_context_msg

!USUS  INTEGER, PARAMETER ::  SIZE_MPE_CONTEXT_MSG = MAX_NAME_LEN/4 + 6 + ISAFETY_MARGIN ! (unit: INTEGERs)
INTEGER, PARAMETER ::  SIZE_MPE_CONTEXT_MSG = MAX_NAME_LEN + 6 + ISAFETY_MARGIN ! (unit: INTEGERs)

!--- ready message handling (compute PE #0 -> IO PE #0) -----------------------

! variable counting opened files, s.t. it is possible to check if
! the last closed file is the one that was meant by the call to
! "mpe_io_ready":
INTEGER                               ::                   &
  icurrent_file

! Flag. If .TRUE., then compute PE #0 has already sent a notification
! message to IO PE #0, whether to write a ready file or not.
LOGICAL                               ::  lready_signal_sent

! buffer for context message received by IO PE #0 from compute PE #0
INTEGER            ::  ireadymsg(SIZE_MPE_CONTEXT_MSG)

!--- Linked list buffer for fields received from compute PEs ------------------

TYPE linked_list_type
  INTEGER(KIND=gribf_kind), POINTER  :: field(:)
  INTEGER                            :: field_tag      ! sorting key = MPI tag = field no.
  INTEGER                            :: field_len      ! unit is in bytes!
  INTEGER                            :: istorage       ! if >=0 field lies in preallocated storage
  TYPE(linked_list_type), POINTER    :: next
END TYPE linked_list_type

TYPE buffer_type
  TYPE (linked_list_type), POINTER   :: head           ! buffer list head
  INTEGER                            :: length         ! no. of entries in buffer
  INTEGER                            :: nudat          ! file descriptor
  INTEGER                            :: next_work_tag  ! next list tag for write-out
  INTEGER                            :: next_work_tag0 ! list tag of first compute PE
  LOGICAL, POINTER                   :: completed(:)   ! ith entry TRUE, if write-end-message recv'd
END TYPE buffer_type

TYPE (buffer_type) :: write_buffer  ! dynamically growing field buffer (on IO PEs)
TYPE (buffer_type) :: read_buffer   ! dynamically growing field buffer (on compute PEs)

!--- pre-allocated pool of fields for faster memory handling ------------------

!<< size of pre-allocated field buffer, set to zero to disable pre-allocation
INTEGER, PARAMETER :: MAX_STORAGE_PREALLOC = 1

! field buffer, pre-allocated for both, read and write data
INTEGER(KIND=gribf_kind), SAVE, TARGET, ALLOCATABLE  ::  field_storage(:,:) ! (zero-based array)
LOGICAL,                  SAVE, ALLOCATABLE          ::  loccupied(:)       ! (zero-based array)

!--- FIFO queue of receive requests to be launched ----------------------------

! Since we buffer only MAX_RECV_REQUESTS different IRECV's at once, an
! additional data structure is required to keep the rest of not-yet-launched
! receive requests:

TYPE recv_request_type
   INTEGER                           :: count          ! number of elements in receive buffer
   INTEGER                           :: source         ! rank of source
   INTEGER                           :: tag            ! message tag
END TYPE recv_request_type

TYPE recv_fifo_type
   TYPE(recv_request_type), ALLOCATABLE :: contents(:)
   INTEGER                           :: front          ! front element index in queue
   INTEGER                           :: count          ! number of elements in queue
END TYPE recv_fifo_type

TYPE (recv_fifo_type)        :: recv_fifo              ! FIFO queue of receive requests to be launched


!------- status parameters for error and diagnosis ----------------------------

INTEGER                      :: MPE_DBG_LEVEL = 2          !<< output verbosity level

  ! preprocessor macro for debugging print-out:
#define DEBUG_WRITE(lvl) IF (MPE_DBG_LEVEL>lvl) \
 WRITE (*,*) "proc [",me_compute,",",me_io,"] :  ",

!==============================================================================

CONTAINS

!==============================================================================
!==============================================================================
!+ Initialization of the io module (Constructor)
!------------------------------------------------------------------------------

SUBROUTINE mpe_io_init (icomm_ori, icompute_pes, io_pes, ibuf_len, &
                        icomm_compute, idebug_lvl, ierror)

!------------------------------------------------------------------------------
!
! Definitions  : compute PE : PE that does computations and calls
!                             mpe_io_write to store arrays 
!
!                io      PE : PE that receives data from compute 
!                             PEs and dumps it directly to a 
!                             device
!                             PE can be io- and compute PE 
!
!
! Description  : Build the communicator for this 
!                io system and a new communicator for compute-PEs.
!                Initialize the chosen io devices.     
!                          
!                This subroutine has to be called after MPI_INIT! 
!
! Method       : 
!
!                There are 2 IO Configurations:
!
!                1. Every io PE is a compute PE ( traditional method )
!                2. Every io PE is not a compute PE ( separate io PEs )
!
!                The two configurations are determined by the number 
!                of compute PEs ( icompute_pes ), the number 
!                of io PEs ( io_pes ) together with the number 
!                of PEs inside icomm_ori (ori_pes).
!                
!                1: icompute_pes          = ori_pes
!                2: icompute_pes + io_pes = ori_pes
!
!                Currently io_pes <= icompute_pes is necessary.
!
!                Communication is made more safe by using different
!                communicators for different (sub)systems. 
! 
!                Communicator:       Purpose:
!                --------------------------------------------
!
!                MPI_COMM_COMPUTE    communicator for computation PEs
!                MPI_COMM_IO         communicator for io PEs
!                MPI_COMM_CO2IO      communicator for communication
!                                    between io and compute PEs
!
!                In case of separate io and compute group 
!                MPI_COMM_CO2IO is an intercommunicator.
!
!                MPE_COMM_CONTEXT    communicator between a compute PE and
!                                    a pure IO PE, for messages that each
!                                    IO PE must receive only once.
!                MPE_COMM_CLOSEFILE  communicator with IO #0 as rank 0 proc.
!                                    The task of this communicator is
!                                    to collect close-file messages from
!                                    other IO PEs.
!
! error treatment:
!
!             -  MPI is not initialized 
!             -  Use of a nonsense configuration                            
!             -  Device failed ( no Database response )    
!             -  allocation error   
!
!------------------------------------------------------------------------------

! Parameters
INTEGER, INTENT(IN)  :: &
  icomm_ori       ,     & ! original communicator
  icompute_pes    ,     & ! number of PEs for computation
  io_pes          ,     & ! number of PEs for io
  ibuf_len        ,     & ! length of integer buffer for grib fields
  idebug_lvl              ! debug level for output

INTEGER, INTENT(OUT) :: &
  icomm_compute           ! new communicator for compute PEs
INTEGER, INTENT(INOUT) :: &
  ierror                  ! error return value (will not be overwritten)

! Local parameters :
INTEGER              :: &
  num_ori         ,     & ! number of PEs in ori. communicator
  key             ,     & ! key for MPI_COMM_SPLIT
  color_io        ,     & ! color   MPI_COMM_SPLIT (I/O)
  color_comp      ,     & ! color   MPI_COMM_SPLIT (compute)
  ierr            ,     & ! error parameter
  r_leader        ,     & ! remote leader PE for intra comm..
  comm_local      ,     & ! local communicator
  color           ,     & ! color value for communicator split
  i               ,     & ! loop counter
  iclose_group    ,     & ! group ID for close-file communicator
  iori_group              ! group ID for original communicator
INTEGER, ALLOCATABLE :: isource_ranks(:)

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  icomm_compute = 0
  MPE_DBG_LEVEL = idebug_lvl

  ! copy dummy parameter to module variable
  num_io      = io_pes           
  num_compute = icompute_pes

  ! determine number of PEs
  CALL MPI_COMM_SIZE( icomm_ori, num_ori, ierr) ; CALL CHKERR(ierr, ierror, 1)
  CALL MPI_COMM_RANK( icomm_ori, me_ori , ierr) ; CALL CHKERR(ierr, ierror)

  ! correct nonsense configurations
  IF (num_io > num_ori) num_io = num_ori

  ! Currently, no. of IO procs <= no. of compute procs is necessary:
  IF (num_compute < num_io) THEN
    CALL CHKERR(ierrcode=1, iglobalerr=ierror)
  END IF

  IF (num_ori == 1) THEN   ! no sep io PEs
    num_compute = num_ori
    num_io      = 1
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Initialization of the io-configurations
!------------------------------------------------------------------------------

  IF (num_compute == num_ori) THEN    

    !--------------------------------------------------------------------------
    ! Section 2.1: First configuration: included IO PEs
    !--------------------------------------------------------------------------

    io_config = 1

    ! in "organize_data::input_ioctl", the value for "num_io" is set
    ! to "1", but we check this here for safety reasons:
    IF (num_io /= 1) THEN
      CALL CHKERR(ierrcode=1, iglobalerr=ierror)
    END IF

    IF( me_ori < num_io ) THEN 
      is_io_pe    = .TRUE.
      color_io    =  0
    ELSE
      is_io_pe    = .FALSE.
      color_io    =  MPI_UNDEFINED 
    ENDIF

    ! define I/O communicator
    key = 0
    CALL MPI_COMM_SPLIT( icomm_ori, color_io, key, MPE_COMM_IO, ierr)
    CALL CHKERR(ierr, ierror)

    IF ( is_io_pe ) THEN
      CALL MPI_COMM_RANK( MPE_COMM_IO, me_io, ierr)
      CALL CHKERR(ierr, ierror)
    ELSE
      me_io = -1
    ENDIF

    ! compute communicator is unchanged
    icomm_compute = icomm_ori
    is_compute_pe = .TRUE.

    !US ! "inter" communicator MPE_COMM_CO2IO is just icomm_compute

  ELSEIF ((num_compute + num_io) == num_ori) THEN 

    !--------------------------------------------------------------------------
    ! Section 2.2: Second configuration: separate IO PEs
    !--------------------------------------------------------------------------

    io_config = 2
    r_leader  = num_ori - num_io ! remote leader PE

    IF( me_ori >= r_leader ) THEN     ! last PEs = I/O PE
      is_io_pe      = .TRUE.
      color_io      =  0
      is_compute_pe = .FALSE.
      color_comp    =  MPI_UNDEFINED
    ELSE
      is_io_pe      = .FALSE.
      color_io      =  MPI_UNDEFINED
      is_compute_pe = .TRUE.
      color_comp    =  0 
    ENDIF

    ! define I/O communicator
    key = 0
    CALL MPI_COMM_SPLIT( icomm_ori, color_io, 0, MPE_COMM_IO, ierr)
    CALL CHKERR(ierr, ierror)

    IF ( is_io_pe ) THEN
      CALL MPI_COMM_RANK( MPE_COMM_IO, me_io, ierr)
      CALL CHKERR(ierr, ierror)
      me_compute = -1
    ELSE
      me_io = -1
    ENDIF

    ! PE 0,..., num_pe-io_pes-1 are compute PEs
    CALL MPI_COMM_SPLIT( icomm_ori, color_comp, key, icomm_compute, ierr ) 
    CALL CHKERR(ierr, ierror)

    ! io group and comm group are disjunct the will communicate
    ! via the true inter communicator MPE_COMM_CO2IO
    ! first we have to build the inter communicator
    IF (is_compute_pe) THEN
      comm_local = icomm_compute
    ELSE
      comm_local = MPE_COMM_IO
      r_leader   = 0
    ENDIF

  ELSE

    !--------------------------------------------------------------------------
    ! Section 2.3: Error in all other cases
    !--------------------------------------------------------------------------

    WRITE(*,*) " non implemented config: !"
    WRITE(*,*) " num_compute =            ", icompute_pes
    WRITE(*,*) " num_io      =            ", num_io
    WRITE(*,*) " num_ori     =            ", num_ori

    CALL CHKERR(ierrcode=5, iglobalerr=ierror)
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Duplicate communicator for internal use
!------------------------------------------------------------------------------

  IF( is_compute_pe )  THEN
    ! duplicate this communicator for internal use
    CALL MPI_COMM_DUP  (icomm_compute, MPE_COMM_COMPUTE, ierr) ; CALL CHKERR(ierr, ierror)
    CALL MPI_COMM_SIZE (icomm_compute, num_compute     , ierr) ; CALL CHKERR(ierr, ierror)
    CALL MPI_COMM_RANK (icomm_compute, me_compute      , ierr) ; CALL CHKERR(ierr, ierror)
  ENDIF

  ! now the values me_compute, me_io are determined and we can do the first
  ! debug print:
  IF (MPE_DBG_LEVEL > 5) THEN
    WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',            &
      'MPE_IO2 init: initialized IO configuration: ', io_config
  ENDIF

!------------------------------------------------------------------------------
! Section 4: Allocate memory for Buffer management
!------------------------------------------------------------------------------

  ! Note: The ibuf_len already contains a safety margin
  ! (required, e.g., for endian problem)
  BUFF_length = ibuf_len

  IF (MPE_DBG_LEVEL > 10) THEN
    WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',            &
      'MPE_IO2 init: buffer length = ', BUFF_length
  ENDIF

  ! --- Initialize max- and current-requests

  isend_current_request = 1     ! initial send request
  irecv_current_request = 1     ! initial receive request

  ! --- Allocate receive/send buffers
  ! note: buffers are 0-indexed, since position 0 contains a message flag!
  IF (is_io_pe) THEN
    ALLOCATE (irecv_buffer   (0:BUFF_length,MAX_RECV_REQUESTS) , STAT = ierr)
    CALL CHKERR(ierr, ierror)
    ALLOCATE (irecv_requests(               MAX_RECV_REQUESTS)  , STAT = ierr)
    CALL CHKERR(ierr, ierror)
    ALLOCATE (irecv_isource(                MAX_RECV_REQUESTS)  , STAT = ierr)
    CALL CHKERR(ierr, ierror)
    ALLOCATE (irecv_status(MPI_STATUS_SIZE, MAX_RECV_REQUESTS)  , STAT = ierr)
    CALL CHKERR(ierr, ierror)
    irecv_requests(:) = MPI_REQUEST_NULL
    irecv_isource(:)  = -1
    CALL init_recv_fifo(recv_fifo, 2*num_compute, ierror)
  ENDIF

  IF (is_compute_pe .OR. (me_io == 0)) THEN
    ALLOCATE (isend_buffer   (0:BUFF_length, MAX_SEND_REQUESTS), STAT = ierr)
    CALL CHKERR(ierr, ierror)
    ALLOCATE (isend_requests              (MAX_SEND_REQUESTS)  , STAT = ierr)
    CALL CHKERR(ierr, ierror)
    ALLOCATE (isend_status   (MPI_STATUS_SIZE)                 , STAT = ierr)
    CALL CHKERR(ierr, ierror)
    isend_requests(:) = MPI_REQUEST_NULL
  ENDIF
  ! Allocate buffers for "close-file" messages from IO PEs to IO #0
  IF (is_io_pe .AND. (me_io == 0)) THEN
    ALLOCATE (iclose_rcv_request(num_io)                       , STAT = ierr)
    CALL CHKERR(ierr, ierror)
    iclose_rcv_request(:) = MPI_REQUEST_NULL
    ALLOCATE (iclose_msg(1, num_io)                            , STAT = ierr)
    CALL CHKERR(ierr, ierror)
    ALLOCATE (iclose_rcv_status(MPI_STATUS_SIZE, num_io)       , STAT = ierr)
    CALL CHKERR(ierr, ierror)
  ENDIF

  IF (is_compute_pe) THEN
    ALLOCATE (iread_recv_buffer(0:BUFF_length)                 , STAT = ierr)
    CALL CHKERR(ierr, ierror)
    iread_recv_request = MPI_REQUEST_NULL
  END IF

  ! introduce MPI communicator for close-file messages (note the
  ! rank reordering)
  color = 0
  key   = 2
  IF (is_compute_pe) THEN
    key = 1
  END IF
  IF (is_io_pe .AND. (me_io == 0)) THEN
    key = 0 
  END IF
  ! ordering: [ IO#0 ; compute PEs ; IO PEs ]
  CALL MPI_COMM_SPLIT(icomm_ori, color, key, MPE_COMM_CLOSEFILE, ierr)
  CALL CHKERR(ierr, ierror)

  ! Reset counter. This variable counts opening of files, s.t. it is
  ! possible to check if the last closed file is the one that was
  ! meant by the call to "mpe_io_ready".
  icurrent_file  = 1

  ! create rank translation tables
  IF (MPE_DBG_LEVEL > 10) THEN
    WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',            &
      'MPE_IO2 init: Query MPI comm groups.'
  ENDIF

  CALL MPI_COMM_GROUP(MPE_COMM_CLOSEFILE, iclose_group,  ierr) ; CALL CHKERR(ierr, ierror)
  CALL MPI_COMM_GROUP(icomm_ori,          iori_group,    ierr) ; CALL CHKERR(ierr, ierror)

  IF (MPE_DBG_LEVEL > 10) THEN
    WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',            &
      'MPE_IO2 init: Create rank translation tables.'
  ENDIF

  ALLOCATE(isource_ranks(num_ori),                 stat=ierr)
  CALL CHKERR(ierr, ierror)
  ALLOCATE(iclose_comm_compute_ranks(num_compute), stat=ierr)
  CALL CHKERR(ierr, ierror)
  ALLOCATE(iclose_comm_io_ranks(num_io),           stat=ierr)
  CALL CHKERR(ierr, ierror)
  isource_ranks = (/ ( i, i=0,(num_ori-1)) /)
  ! process rank translation mpe_comm_compute -> MPE_COMM_CLOSEFILE
  CALL MPI_GROUP_TRANSLATE_RANKS(iori_group, num_compute, isource_ranks(:), iclose_group, &
    &                            iclose_comm_compute_ranks(1:num_compute), ierr)
  CALL CHKERR(ierr, ierror)
  ! process rank translation mpe_comm_io      -> MPE_COMM_CLOSEFILE
  IF (io_config == 2) THEN
    CALL MPI_GROUP_TRANSLATE_RANKS(iori_group, num_io,                                           &
      &                            isource_ranks((num_compute+1):(num_compute+num_io)),          &
      &                            iclose_group, iclose_comm_io_ranks(1:num_io), ierr)
    CALL CHKERR(ierr, ierror)
  ELSE
    CALL MPI_GROUP_TRANSLATE_RANKS(iori_group, num_io,                                           &
      &                            isource_ranks(1:num_io), iclose_group,                        &
      &                            iclose_comm_io_ranks(1:num_io), ierr)
    CALL CHKERR(ierr, ierror)
  END IF
  DEALLOCATE(isource_ranks)

  ! allocate context communicator pairs:
  color = MPI_UNDEFINED
  IF (is_compute_pe .AND. (me_compute < num_io)) THEN    
    color = me_compute
    key   = 0 ! 0: sender
  ELSE IF (is_io_pe) THEN
    color = me_io
    key   = 1 
  END IF
  CALL MPI_COMM_SPLIT(icomm_ori, color, key, MPE_COMM_CONTEXT, ierr) ; CALL CHKERR(ierr, ierror)

  ! Create MPI communicators for read and write data exchange
  ! (Simply duplicates of icomm_ori)
  CALL MPI_COMM_DUP(icomm_ori, MPE_COMM_READ , ierr) ; CALL CHKERR(ierr, ierror)
  CALL MPI_COMM_DUP(icomm_ori, MPE_COMM_WRITE, ierr) ; CALL CHKERR(ierr, ierror)

  ! debug print-out
  IF (MPE_DBG_LEVEL > 10) THEN
    WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',            &
      'MPE_IO2 init: MPI communicators:'              , &
       '  icomm_ori          = ', icomm_ori           , &
       ', icomm_compute      = ', icomm_compute       , &
       ', mpe_comm_io        = ', MPE_COMM_IO         , &
       ', mpe_comm_compute   = ', MPE_COMM_COMPUTE    , &
       ', mpe_comm_co2io     = ', MPE_COMM_CO2IO      , &
       ', MPE_COMM_CLOSEFILE = ', MPE_COMM_CLOSEFILE  , &
       ', MPE_COMM_CONTEXT   = ', MPE_COMM_CONTEXT         
  ENDIF

  ! Initialize flag which prevents prefetching loop from being
  ! started without need:
  lprefetching_available = .FALSE.
  lprefetching_started   = .FALSE.

  ! --- pre-allocate some memory for read/write buffers:
  ALLOCATE (loccupied    (             0:(MAX_STORAGE_PREALLOC-1)), STAT = ierr)
  CALL CHKERR(ierr, ierror)
  ALLOCATE (field_storage(BUFF_length, 0:(MAX_STORAGE_PREALLOC-1)), STAT = ierr)
  CALL CHKERR(ierr, ierror)
  loccupied(:) = .FALSE.

  ! --- initialize ready file flag
  lready_signal_sent = .TRUE. 

  IF (MPE_DBG_LEVEL > 5) THEN
    WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',            &
      'MPE_IO2 init finished '
  ENDIF

END SUBROUTINE mpe_io_init

!==============================================================================
!==============================================================================
!+ Opens Files only on IO PEs
!------------------------------------------------------------------------------

SUBROUTINE mpe_io_open (nudat, ydatname, ymode, lstop_io, ierror)

!------------------------------------------------------------------------------
!
! Description:
!  Called by all PEs (compute and IO). Besides handling the file
!  opening (IO processes), this subroutine exchanges status
!  information ("context messages"): file name, mode;
!  it issues requests for a mpe_io_ready()-event and asks the
!  contributing compute PEs for the first fields.
!------------------------------------------------------------------------------

! Parameters
INTEGER             , INTENT(INOUT) ::  &
  nudat           ! file descriptor
CHARACTER (LEN = *) , INTENT(IN)    ::  &
  ydatname   ,  & ! name of data file
  ymode           ! file mode
LOGICAL             , INTENT(OUT)   ::  &
  lstop_io        ! tells pure IO nodes to stop
INTEGER             , INTENT(INOUT) ::  &
  ierror          ! error return value (will not be overwritten)

! Local parameters
TYPE(mpe_context_msg)      ::   &
  current_context_msg             ! context message
INTEGER                    ::   &
  i,                            & ! loop counter
  dat_len,                      & ! length for datname and mode
  ierr,                         & ! local error parameter
  status(MPI_STATUS_SIZE),      &
  compute_pe, io_pe, field_tag, &
  source_pe
INTEGER :: imsg(SIZE_MPE_CONTEXT_MSG) ! compressed context message
LOGICAL :: lskip_context_msg, lready_file
TYPE(mpe_context_msg)      ::   &
  ready_context_msg ! context message containing ready file info

!------------------------------------------------------------------------------

IF (MPE_DBG_LEVEL > 3) THEN
  WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',            &
              'step into mpe_io_open, mode = ', ymode(1:1)
ENDIF

!------------------------------------------------------------------------------
! --- receive context messages

  lready_file = .FALSE.

  current_context_msg%ymode = ymode
  SELECT CASE (current_context_msg%ymode(2:2))
  CASE ('g',' ')
    current_context_msg%ylib = 'grb1'
  CASE ('x')
    current_context_msg%ylib = 'apix'
  CASE ('1')
    current_context_msg%ylib = 'api1'
  CASE ('2')
    current_context_msg%ylib = 'api2'
  END SELECT

  lskip_context_msg = (current_context_msg%ymode(1:1) == 'r' ) &
    &                 .AND.  is_compute_pe .AND. (me_compute /= 0)

  lstop_io       = .FALSE.

  IF (.NOT. lskip_context_msg) THEN

    dat_len   = LEN_TRIM(ydatname)
    current_context_msg%id                   = MSG_ID_DATNAME
    current_context_msg%length               = dat_len
    current_context_msg%yfilename(1:dat_len) = ydatname(1:dat_len)

    IOCONFIG2 : IF ( io_config == 2 ) THEN
      ! datname and mode have to be communicated to io nodes
      ! code datname and mode on dat_info

      ISCOMPUTE : IF ( is_compute_pe ) THEN

        ! send dat_info to IO PEs
        IF( me_compute < num_io ) THEN

          DEBUG_WRITE(3) &
            'send open message'

          CALL compress_mpe_context_msg(current_context_msg, imsg)
          CALL MPI_SEND(imsg, SIZE_MPE_CONTEXT_MSG, MPI_INTEGER, 1, TAG_OPEN, &
                        MPE_COMM_CONTEXT, ierr ) ; CALL CHKERR(ierr, ierror, 1)
        ENDIF

      ELSE

        OPENLOOP : DO

          CALL MPI_RECV (imsg, SIZE_MPE_CONTEXT_MSG, MPI_INTEGER, 0, TAG_OPEN, &
                         MPE_COMM_CONTEXT, status, ierr ) ; CALL CHKERR(ierr, ierror)
          ! decode datname (or input_order) and mode
          CALL decode_mpe_context_msg(imsg, current_context_msg)

          DEBUG_WRITE(3) &
            'received open message with dat_len = ', dat_len

          IF (current_context_msg%id == MSG_ID_DATNAME) THEN
            
            ! -- Case #1 : simple "open-file" message:

            EXIT OPENLOOP
          ELSE IF (current_context_msg%id == MSG_ID_STOP) THEN

            ! -- Case # 2: end message, i.e. compute PEs have
            ! finished the forecast and send a message (in
            ! mpe_io_shutdown) to stop the IO-PEs

            DEBUG_WRITE(10) &
              "Received shutdown signal."
            lstop_io = .TRUE.
            RETURN
          ELSE IF (current_context_msg%id == MSG_ID_WRITE_READY_FILE) THEN

            ! -- Case #3: ready-file message, ie. the SR
            ! "mpe_io_ready" has been called on compute PE #0, which
            ! notifies IO PE #0:

            IF (current_context_msg%itag <= icurrent_file) THEN

              ! Wait for closing files of other IO PEs:
              DEBUG_WRITE(3) &
                "Stopped open message loop, wait for closing files, tag = ", icurrent_file
              CALL MPI_WAITALL( num_io, iclose_rcv_request(:),  &
                                iclose_rcv_status(:,:), ierr ) ; CALL CHKERR(ierr, ierror)
              DEBUG_WRITE(3) &
                "... all messages arrived"

              ! create ready file:        
              CALL DEV_ready ( current_context_msg%iunit       , &
                &              current_context_msg%yfilename   , &
                &              current_context_msg%length )
              lready_file = .TRUE.
              
            ELSE
              
              ! the "ready" event corresponds to a file which has
              ! not yet been closed by this IO PE

              DEBUG_WRITE(3) &
                "ready message tag = ", current_context_msg%itag, ", ", &
                "current tag (IO)  = ", icurrent_file
    
            END IF
          ELSE IF (current_context_msg%id == MSG_ID_DONT_WRITE_READY_FILE) THEN

            ! -- Case #4: On compute PE #0 the SR "mpe_io_ready" was
            ! NOT invoked, ie. the recently closed file was an
            ! ordinary one (not the last one).

            DEBUG_WRITE(3) &
              "Receive DONT_WRITE_READY_FILE"

          ELSE
            CALL CHKERR(ierrcode=1, iglobalerr=ierror) ! Error: Unknown context message
          ENDIF
        END DO OPENLOOP

        DEBUG_WRITE(3) &
          ' ymode = ',current_context_msg%ymode,' in_device = file'
      ENDIF ISCOMPUTE

    ENDIF IOCONFIG2

  END IF

  current_mode = current_context_msg%ymode(1:1)
  SELECT CASE (current_context_msg%ymode(3:3))
  CASE ('s',' ')
    current_speed = 's'
  CASE ('b')
    current_speed = 'b'
  CASE ('f')
    current_speed = 'f'
  END SELECT
  current_lib  = current_context_msg%ylib (1:4)

  IF ((current_context_msg%ymode(1:1) /= 'r')  .AND.   &
      (current_context_msg%ymode(1:1) /= 'w')  .AND.   &
      (current_context_msg%ymode(1:1) /= 'a')) THEN
    ! Error: Access mode unknown
    CALL CHKERR(ierrcode=1, iglobalerr=ierror)
  END IF

  ! only the IO PEs actually open the file
  IF (is_io_pe) THEN
    DEBUG_WRITE(3) &
      "OPEN file ", current_context_msg%yfilename(1:current_context_msg%length)
    DEBUG_WRITE(3) &
      'file ymode_new = ', current_context_msg%ymode, num_io
    CALL DEV_Open(nudat, current_context_msg%yfilename, current_context_msg%length, &
                  current_context_msg%ymode, ierror)
  ENDIF

  !------------------------------------------------------------------------------
  ! --- initialize data structures for write mode

  IF (is_io_pe) THEN
    IF (current_context_msg%ymode(1:1) == 'w' .OR. current_context_msg%ymode(1:1) == 'a') THEN
      CALL buffer_list_init(nudat, write_buffer, ierror)
    ENDIF
  ENDIF

  IF (is_compute_pe .AND. ( (current_mode=='w') .OR. (current_mode=='a') )) THEN
    field_snd = me_compute - num_compute 
  END IF

  ! if in write mode, launch initial irecv's on IO PEs:
  IF (is_io_pe .AND. ( (current_mode=='w') .OR. (current_mode=='a') )) THEN
    DO compute_pe = 0, (num_compute-1)
       ! the following formula defines which fields are RECEIVED
       ! from a given compute PE:
       field_tag = compute_pe
       IF ( (MOD(field_tag, num_io) == me_io) .AND.   &
            (.NOT. (compute_pe == me_compute))) THEN
          ! Note: including the MPE_FLAG, we receive (BUFF_length+1) values
          CALL push_recv_fifo(recv_fifo, icount=(BUFF_length+1), isource=compute_pe, &
               &              itag=field_tag, ierror=ierror)
       END IF
    END DO
    DEBUG_WRITE(10) "initial irecvs"
    CALL launch_irecvs(recv_fifo, ierror)
  END IF

  !------------------------------------------------------------------------------
  ! --- initialize data structures for read mode

  IF (.NOT. lprefetching_started) THEN
    
    IF (is_compute_pe) THEN
      CALL buffer_list_init(nudat, read_buffer, ierror)
    ENDIF

    IF ( is_io_pe  .AND. (current_mode=='r') ) THEN
      field_snd = 0
    END IF

    ! if in read mode, launch initial irecv's on compute PEs:
    IF ( is_compute_pe  .AND.  (current_mode=='r') ) THEN

      IF (io_config == 2) THEN
        io_pe = num_compute
      ELSE
        io_pe = 0
      END IF

      IF (io_pe /= me_io) THEN
        field_rcv = me_compute
        IF (MPE_DBG_LEVEL > 10) THEN
          WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',                  &
            'Launch IRECV(tag= ', field_rcv, ');  iread_recv_request = ', iread_recv_request
        ENDIF
        ! Note: including the MPE_FLAG, we receive (BUFF_length+1) values
        CALL MPI_IRECV(iread_recv_buffer(:), (BUFF_length+1),            &
                       MPI_DATATYPE, io_pe, field_rcv, MPE_COMM_READ,    &
                       iread_recv_request, ierr) ; CALL CHKERR(ierr, ierror)
        IF (MPE_DBG_LEVEL > 10) THEN
          WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',                  &
            'Now launched IRECV(tag= ', field_rcv, ');  iread_recv_request = ', iread_recv_request
        ENDIF
        field_rcv = field_rcv + num_compute
      END IF
    END IF

  END IF

  !------------------------------------------------------------------------------
  ! --- ready message handling:

  IF (is_io_pe  .AND.  (me_io == 0)  .AND.  &
      ( (current_context_msg%ymode(1:1) == 'w') .OR.  &
        (current_context_msg%ymode(1:1) == 'a')) ) THEN

    DEBUG_WRITE(3) &
      "Issue close-file-requests on IO PEs, tag=", icurrent_file
    
    ! compute PE #0 receives "close file" messages from all IO PEs:
    DO i=1, num_io
      source_pe = iclose_comm_io_ranks(i)
      
      ! (From the ordering of MPE_COMM_CLOSEFILE, we know, that IO PE#0
      ! corresponds to i==1)
      IF ((i==1) .AND. (source_pe /= 0)) THEN
        CALL CHKERR(ierrcode=1, iglobalerr=ierror)
      END IF
      
      IF (source_pe /= me_io) THEN
        CALL MPI_IRECV (iclose_msg(:,i), 1, MPI_INTEGER, source_pe,                   &
          &             icurrent_file, MPE_COMM_CLOSEFILE, iclose_rcv_request(i), ierr )
        CALL CHKERR(ierr, ierror)
      END IF
    END DO
    icurrent_file = icurrent_file + 1
  END IF

  ! "release" IO PE #0 waiting for information whether to write a
  ! ready file or not:
  IF ( is_compute_pe  .AND.  (me_compute == 0)  .AND.   &
       (io_config == 2)                         .AND.   &
       .NOT. lready_signal_sent) THEN

    DEBUG_WRITE(3) &
      "Send DONT_WRITE_READY_FILE"

    ! Send a message to IO PE#0, indicating that NO ready file
    ! should be written.
    ready_context_msg%id     = MSG_ID_DONT_WRITE_READY_FILE
    ready_context_msg%length = 0 
    ready_context_msg%itag   = -1

    CALL compress_mpe_context_msg(ready_context_msg, imsg)

    DEBUG_WRITE(3) &
      "Sending message for DONT ready notification from compute PE #0", &
      ", tag = ", icurrent_file

    CALL MPI_SEND(imsg, SIZE_MPE_CONTEXT_MSG, MPI_INTEGER, 1,           &
                   TAG_OPEN, MPE_COMM_CONTEXT,  ierr) ; CALL CHKERR(ierr, ierror)
    icurrent_file = icurrent_file + 1
  
  END IF

  DEBUG_WRITE(3) &
    "step off mpe_io_open"

END SUBROUTINE mpe_io_open

!==============================================================================
!==============================================================================
! Actually opens file for reading/writing
!------------------------------------------------------------------------------

SUBROUTINE DEV_Open( nudat, datname, dat_len, filemode, ierror )

  ! Parameters
  INTEGER           , INTENT(INOUT)   :: nudat
  CHARACTER  (LEN=*), INTENT(IN)      :: datname
  INTEGER           , INTENT(IN)      :: dat_len
  CHARACTER         , INTENT(INOUT)   :: filemode
  ! error return value (will not be overwritten):
  INTEGER           , INTENT(INOUT)   :: ierror   

  ! Local parameters
  INTEGER                             :: new_len
  CHARACTER    (LEN= dat_len+3)       :: datname_new
  INTEGER      (KIND=gribc_kind)      :: &
    nudatc,     & ! local file descriptor
    ierrc         ! local file error
  CHARACTER    (LEN=3)                :: mode

  ! copy mode character to dummy array:
  mode = '   '
  mode(1:1) = filemode

  new_len     = dat_len
  datname_new = datname(1:dat_len)

  nudatc = nudat

  ! open (possibly multiple) files for writing:
  IF(current_mode == 'w' .OR. current_mode == 'a' ) THEN

    IF( num_io > 1 ) THEN

      IF( me_io < 10 ) THEN
        WRITE( datname_new(dat_len+1:dat_len+2),'(A1,I1)') '_',me_io
      ELSE IF( me_io < 100 ) THEN
        WRITE( datname_new(dat_len+1:dat_len+3),'(A1,I2)') '_',me_io
      ELSE IF( me_io < 1000 ) THEN
        WRITE( datname_new(dat_len+1:dat_len+4),'(A1,I3)') '_',me_io
      ENDIF

      new_len = dat_len+3

    ENDIF

    IF(current_mode == 'w' ) THEN
      DEBUG_WRITE(10) &
      "Open file ", datname_new(1:new_len), " for writing."
    ELSE IF(current_mode == 'a' ) THEN
      DEBUG_WRITE(10) &
      "Open file ", datname_new(1:new_len), " for appending."
    END IF

    SELECT CASE (current_lib)

#ifdef GRIBDWD
    CASE ('grb1')
      CALL copen(nudatc,datname_new(1:new_len), mode(1:3), ierrc)
#endif

#ifdef GRIBAPI
    CASE ('api1','api2')
      CALL grib_open_file(nudatc, datname_new(1:new_len), current_mode, ierrc)
      IF (ierrc /= GRIB_SUCCESS) THEN
        IF ( current_mode == 'w' ) THEN
          PRINT *, ' *** ERROR opening grib_api file for writing'
        ELSE IF ( current_mode == 'a' ) THEN
          PRINT *, ' *** ERROR opening grib_api file for appending'
        ENDIF
      ENDIF
#endif

    CASE DEFAULT
      PRINT *, ' *** ERROR: Wrong value for data format: ', current_lib
    END SELECT
    CALL CHKERR(ierrc, ierror, 1)

  ENDIF

  ! open input file for reading
  IF ( (current_mode == 'r')  .AND.  (me_io == 0) ) THEN

    SELECT CASE (current_lib)

#ifdef GRIBDWD
    CASE ('grb1')
      CALL copen(nudatc,datname_new(1:new_len), mode(1:3), ierrc)
#endif

#ifdef GRIBAPI
    CASE ('apix')
      CALL grib_open_file(nudatc, datname_new(1:new_len), 'r', ierrc)
      IF (ierrc /= GRIB_SUCCESS) THEN
        PRINT *, ' *** ERROR opening grib_api file for reading'
      ENDIF
#endif

    CASE DEFAULT
      PRINT *, ' *** ERROR: Wrong value for data format: ', current_lib
    END SELECT
    CALL CHKERR(ierrc, ierror, 1)

  ENDIF
  ! return file handle
  nudat = nudatc

END SUBROUTINE DEV_Open

!==============================================================================
!==============================================================================
!+ Closes IO device
!------------------------------------------------------------------------------

SUBROUTINE mpe_io_close (nudat, ierror)

!------------------------------------------------------------------------------
!
! Description:
!    IO PEs close their output/input files and inform the IO PE #0
!    by sending an MPI message.
!
!------------------------------------------------------------------------------

! Parameters
INTEGER, INTENT(IN)    :: nudat  ! file descriptor
INTEGER, INTENT(INOUT) :: ierror ! error return value (will not be overwritten)

!------------------------------------------------------------------------------

! Local parameters
INTEGER  (KIND= gribc_kind ) ::    &
  nudatc,     & ! file descriptor for C-routine      
  ierrc         ! error parameter for C-routine
INTEGER :: imsg, ierr

!------------------------------------------------------------------------------

  DEBUG_WRITE(10) &
    "into mpe_io_close"
  
  ! Error checks:
  ierrc  = 0

  IF (is_io_pe) THEN

    ! close current file
    nudatc = nudat
    IF (  (current_mode == 'w')  .OR.  (current_mode == 'a') .OR. &
         ((current_mode == 'r')  .AND. (me_io == 0)) ) THEN

      SELECT CASE (current_lib)

#ifdef GRIBDWD
      CASE ('grb1')
        CALL cclose(nudatc,'exi',ierrc)
#endif

#ifdef GRIBAPI
      CASE ('apix','api1','api2')
        CALL grib_close_file(nudatc)
#endif
      CASE DEFAULT
        PRINT *, ' *** ERROR: Wrong value for library used:  ', current_lib
      END SELECT
      CALL CHKERR(ierrc, ierror, 1)
    ENDIF

    ! send close-file message to IO PE #0
    IF (me_io /= 0) THEN

      IF (current_mode == 'w' .OR. current_mode == 'a') THEN
        ! (note: checking for write mode should be unnecessary)
        DEBUG_WRITE(3) &
          & "sending close-file message to IO PE #0, tag = ", icurrent_file

        CALL MPI_SEND (imsg, 1, MPI_INTEGER, 0, &
          &            icurrent_file, MPE_COMM_CLOSEFILE, ierr ) ; CALL CHKERR(ierr, ierror)

        icurrent_file = icurrent_file + 1
      END IF

    ENDIF

  END IF

  ! clear read buffer
  IF (is_compute_pe) THEN
    CALL buffer_list_destroy(read_buffer, ierror)
  END IF

  DEBUG_WRITE(10) &
    "Reset prefetching flags."
  lprefetching_available = .FALSE.
  lprefetching_started   = .FALSE.

  ! Reset ready file signal (indication whether compute PE #0 has
  ! sent a message to IO PE #0). This flag will be enabled again
  ! either during a call to "mpe_io_ready()" or, otherwise, during
  ! the next call to "mpe_io_open"
  IF (is_compute_pe .AND. (me_compute == 0) .AND. (io_config == 2) .AND. &
      ((current_mode == 'w') .OR. (current_mode == 'a')) ) THEN
    lready_signal_sent = .FALSE.
  END IF

END SUBROUTINE mpe_io_close

!==============================================================================
!==============================================================================
!+ Completes a mpe_io_write call
!------------------------------------------------------------------------------

SUBROUTINE mpe_io_complete(ierror)

!------------------------------------------------------------------------------
!
! Description:
!  This routine is called by all compute PEs (which might be IO at the same
!  time). The call is issued, when the last grib field has been ISENDed by
!  a compute PE. Since the IO PEs must wait for all their contributing PEs
!  they are informed by a write-end-message.
!
!------------------------------------------------------------------------------

! Parameters:
INTEGER, INTENT(INOUT) :: ierror ! error return value (will not be overwritten)
! Local parameters
INTEGER :: ierr, pe_target
INTEGER :: istatus_list(MPI_STATUS_SIZE, MAX_SEND_REQUESTS)

!------------------------------------------------------------------------------

  DEBUG_WRITE(20) &
    'jump into  complete'

  COMPUTE : IF (is_compute_pe) THEN

    ! Finish pending requests:

    ! Technically, the following WAITALL() statement is not
    ! required. However, it has some performance benefits, since it
    ! "releases" the IO PE #0 which is also a compute PE
    IF (io_config == 1) THEN
      CALL MPI_WAITALL(MAX_SEND_REQUESTS, isend_requests(:), &
        &              istatus_list, ierr) ; CALL CHKERR(ierr, ierror, 1)
    END IF

    field_snd  = field_snd + num_compute
    IF (io_config == 1) THEN
      pe_target  = MOD(me_compute, num_io)
    ELSEIF (io_config == 2) THEN
      pe_target  = MOD(me_compute, num_io) + num_compute
    ENDIF
    
    DEBUG_WRITE(20) &
      & "sending end message (tag=", field_snd, ") in mpe_io_complete to ", &
      & pe_target
    
    ! don't send to myself
    IF (pe_target == me_io) THEN
      write_buffer%completed(me_io) = .TRUE.
      pe_target = MPI_PROC_NULL
    ENDIF
    
    CALL MPI_SEND(MPE_FLAG_WRITE_END, 1, MPI_DATATYPE,      &
      &            pe_target, field_snd, MPE_COMM_WRITE,     &
      &            ierr) ; 
    CALL CHKERR(ierr, ierror)
    
  ENDIF COMPUTE

  IF (is_io_pe) THEN
    CALL buffer_list_finalize_write(write_buffer, ierror)
    CALL buffer_list_destroy(write_buffer, ierror)
  END IF

END SUBROUTINE mpe_io_complete

!==============================================================================
!==============================================================================
!+ checks (and signals) whether a postprocessing step has finished
!------------------------------------------------------------------------------

SUBROUTINE mpe_io_ready (file_name, ifn_len, nunit, ierror)

!------------------------------------------------------------------------------
!
! Description:
!  This subroutine creates a ready file (by IO proc #0) as soon as 
!  IO procs report closing the last file.
!
!------------------------------------------------------------------------------

! Parameters
INTEGER                , INTENT(IN)    :: ifn_len     ! length of file name
INTEGER                , INTENT(IN)    :: nunit       ! unit number for ready file
CHARACTER (LEN=ifn_len), INTENT(IN)    :: file_name   ! filename of the ready file
INTEGER                , INTENT(INOUT) :: ierror      ! error value (will not be overwritten)

! Local Parameters
INTEGER               :: ierr ! local error code
TYPE(mpe_context_msg) ::       &
  ready_context_msg               ! context message containing ready file info
INTEGER               ::       &
  imsg(SIZE_MPE_CONTEXT_MSG)      ! compressed context message
  
!------------------------------------------------------------------------------

  DEBUG_WRITE(3) &
    'jump into  ready file'

  ! "release" IO PE #0 waiting for information whether to write a
  ! ready file or not:
  IF (is_compute_pe .AND. (me_compute == 0)) THEN
    
    IF (io_config == 2) THEN
      ! Send a message to IO PE#0, indicating that a ready file
      ! should be written.

      ready_context_msg%id        = MSG_ID_WRITE_READY_FILE
      ready_context_msg%yfilename = file_name
      ready_context_msg%length    = ifn_len
      ready_context_msg%iunit     = nunit
      ready_context_msg%itag      = icurrent_file

      ! pack ready file name into integer array
      CALL compress_mpe_context_msg(ready_context_msg, imsg)

      ! send message (non-blocking)
      DEBUG_WRITE(3) &
        "Sending message for ready notification from compute PE #0", &
        ", tag = ", icurrent_file

      CALL MPI_SEND(imsg, SIZE_MPE_CONTEXT_MSG, MPI_INTEGER, 1, TAG_OPEN, &
        &            MPE_COMM_CONTEXT, ierr) ; CALL CHKERR(ierr, ierror)
      icurrent_file = icurrent_file + 1

      lready_signal_sent = .TRUE.

    ELSE

      ! for io_config==1: directly create ready file:        
      CALL DEV_ready ( nunit      , &
        &              file_name  , &
        &              ifn_len      )

    END IF

  END IF

  DEBUG_WRITE(3) &
    'jump off  ready file'

END SUBROUTINE mpe_io_ready

!==============================================================================
!==============================================================================
!+ Creates a "ready file", i.e. a small file indicating some event.
!------------------------------------------------------------------------------

SUBROUTINE DEV_ready ( iunit, yfilename, ifn_len )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine is called collectively by all IO processes. The first
!   IO process then creates a ready file.
!
!------------------------------------------------------------------------------

  ! Parameters
  INTEGER           , INTENT(IN)   :: iunit
  CHARACTER  (LEN=*), INTENT(IN)   :: yfilename
  INTEGER           , INTENT(IN)   :: ifn_len

  DEBUG_WRITE(3) &
    "jump into DEV_ready, me_io = ", me_io
  ! write ready file on first IO PE
  IF( is_io_pe .AND. (me_io == 0)) THEN
    
    DEBUG_WRITE(3) &
      "Write ready file ", yfilename(1:ifn_len)
    OPEN (iunit, file=yfilename(1:ifn_len), form='formatted')
    WRITE(iunit, '(A)') 'ready'
    CLOSE(iunit)

  END IF

END SUBROUTINE DEV_ready

!==============================================================================
!==============================================================================
!+ Shutdown io nodes
!------------------------------------------------------------------------------

SUBROUTINE mpe_io_shutdown(ierror)

!------------------------------------------------------------------------------
!
! Description:
!   Subroutine for clean-up of dynamic data structures.
!   In addition, if this process is a compute PE, this routine informs a
!   corresponding pure IO process (me_compute+num_compute), about the
!   shutdown by sending an appropriate context message. There, this ends the 
!   loop in mpe_io_node().
!
!------------------------------------------------------------------------------

! Parameters:
INTEGER, INTENT(INOUT)  :: ierror ! error return value (will not be overwritten)
! Local parameters
TYPE(mpe_context_msg) :: &
  context_msg                                          ! context message
INTEGER :: &
  imsg(SIZE_MPE_CONTEXT_MSG)                       , & ! compressed context message
  ierr                                                 ! local error code

!------------------------------------------------------------------------------

  DEBUG_WRITE(10) &
    "into mpe_io_shutdown"

  IF ((io_config == 2) .AND. is_compute_pe) THEN
      IF ( me_compute < num_io ) THEN
        DEBUG_WRITE(5) &
          "Preparing shutdown message"
        context_msg%id     = MSG_ID_STOP
        context_msg%length = 0
        CALL compress_mpe_context_msg(context_msg, imsg)

        DEBUG_WRITE(5) &
          "Sending shutdown message"
        CALL MPI_SEND (imsg, SIZE_MPE_CONTEXT_MSG, MPI_INTEGER, 1, &
          &            TAG_OPEN, MPE_COMM_CONTEXT, ierr) ; CALL CHKERR(ierr, ierror, 1)
      ENDIF
    ENDIF

  !! CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) ; CALL CHKERR(ierr, ierror, 1)
  !! but MPE_COMM_WRITE should be also icomm_ori=icomm_world
  CALL MPI_BARRIER(MPE_COMM_WRITE,ierr) ; CALL CHKERR(ierr, ierror, 1)

  IF(is_compute_pe .OR. (me_io==0)) THEN
    DEALLOCATE (isend_buffer, isend_requests, isend_status, STAT = ierr)
    CALL CHKERR(ierr, ierror)
  END IF
  IF(is_io_pe) THEN
    DEALLOCATE (irecv_buffer, irecv_requests, irecv_isource, irecv_status, STAT = ierr)
    CALL CHKERR(ierr, ierror)
    CALL finalize_recv_fifo(recv_fifo, ierror)
  END IF
  IF(is_io_pe .AND. (me_io == 0)) THEN
    DEALLOCATE (iclose_rcv_request, iclose_msg, iclose_rcv_status, STAT = ierr)
    CALL CHKERR(ierr, ierror)
  ENDIF

  IF (is_compute_pe) THEN
    DEALLOCATE (iread_recv_buffer, STAT = ierr)
    CALL CHKERR(ierr, ierror)
  END IF

  DEALLOCATE (iclose_comm_compute_ranks,          &
    &         iclose_comm_io_ranks,               &
    &         field_storage,                      &
    &         loccupied,                          &
    &         STAT = ierr)
  CALL CHKERR(ierr, ierror)

END SUBROUTINE mpe_io_shutdown

!==============================================================================
!==============================================================================
!+ Driver routine for the IO PEs
!------------------------------------------------------------------------------

SUBROUTINE mpe_io_node(ierror)

!------------------------------------------------------------------------------
!
! Description:
!  In this subroutine, pure IO PEs are cycling through receive, send, open,
!  close states. The loop is ended by an mpe_io_shutdown() call.
!
!------------------------------------------------------------------------------

! Parameters
INTEGER                      :: ierror ! error return value (will not be overwritten)
! Local parameters :
INTEGER                      :: nudat
INTEGER(KIND=gapi_int)       :: ilength
LOGICAL                      :: lstop_io
CHARACTER (LEN=1)            :: yname_dummy, ymode_dummy

!------------------------------------------------------------------------------

  lstop_io   = .FALSE.

  ENDLESSLOOP: DO

    DEBUG_WRITE(3) &
      'before open in mpe_io_node ', me_io
    CALL mpe_io_open (nudat, yname_dummy, ymode_dummy, lstop_io, ierror)
    DEBUG_WRITE(3) &
      'open on io node done ', me_io

    ! this flag is set after receiving an appropriate context
    ! message
    IF (lstop_io) EXIT ENDLESSLOOP

     IF (current_mode == 'w' .OR. current_mode == 'a') THEN
 
      ! Loop until no irecv requests left:
      DO
        IF (ALL(write_buffer%completed)) EXIT
        
        DEBUG_WRITE(30) &
          'before write on node'
        CALL mpe_io_write(ierror=ierror)
        DEBUG_WRITE(30) &
          'write on node done'
      END DO

      CALL buffer_list_finalize_write(write_buffer, ierror)
      CALL buffer_list_destroy(write_buffer, ierror)
    ENDIF

    IF( current_mode == 'r' ) THEN

      IF (.NOT. lprefetching_available) THEN

        ilength = -1_gapi_int
        ! loop until whole file has been read
        DO
          IF (ilength == 0) EXIT
          
          DEBUG_WRITE(30) &
            'before read on io-node'
          CALL mpe_io_read (nudat, ilength=ilength, ierror=ierror)
          DEBUG_WRITE(30) &
            'read on node done'
        END DO

      ELSE
        DEBUG_WRITE(10) &
          "Prefetching data already available!"
      END IF

    ENDIF

    DEBUG_WRITE(3) &
      'before close on io node', me_io
    CALL mpe_io_close (nudat, ierror)
    DEBUG_WRITE(3) &
      'close on io node done ', me_io

  ENDDO ENDLESSLOOP

  DEBUG_WRITE(3) &
    'io-node =', me_io, ' has been shut down'

  ! clean up
  CALL mpe_io_shutdown(ierror)

END SUBROUTINE mpe_io_node

!==============================================================================
!==============================================================================
!+ Read data on I/O PEs and send it to compute PEs
!------------------------------------------------------------------------------

SUBROUTINE mpe_io_read (nudat, ifield, ilength, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   The implementation uses internal, dynamically growing read buffers
!   stored on the compute PEs and filled with data from the IO PEs.
!   This subroutine serves both as a prefetching mechanism and as a direct
!   read-in:
!   Pre-fetching operation only reads and fills the internal
!   buffer. In normal operation mode, if there are entries
!   available in the internal buffer, these are simply fetched,
!   otherwise we have to wait for completion of the current
!   iread_recv_request.
!
! Note:
!   Sending only one field to each compute PE without acknowledgment might be
!   problematic, since this means that prefetching is rather slow.
!   Think about sending more than one field, but this does not only require
!   larger isend buffers but also a different handling of irecv cancellation.
!   The slow-down might be also less relevant for a realistic simulation with
!   more than 3 compute PEs.
!------------------------------------------------------------------------------

! Parameters
INTEGER, INTENT(IN)                                ::    &
  nudat          ! file descriptor                
INTEGER(KIND=gapi_int), INTENT(OUT)                ::    &
  ilength        ! number of Bytes read           
INTEGER , INTENT(INOUT)                            ::    &
  ierror         ! error return value (will not be overwritten)
INTEGER(KIND=gribf_kind), OPTIONAL, INTENT(INOUT)  ::    &
  ifield(:)       ! array for read data 

! Local Parameters
INTEGER  (KIND= gribc_kind ) ::  &
  nudatc, imaxlenc, ilengthc, ierrc       ! corresponding variables for C-routines
INTEGER        ::   &
      istatus(MPI_STATUS_SIZE), icomp_pes,                &
      io_pe, byte_length,                                 &
      istatus_list(MPI_STATUS_SIZE, MAX_SEND_REQUESTS),   &
      ifirst_zero_record,                                 &
      imaxlen, ilengthf,                                  &    ! buffer size in BYTES!!
      dat_length,                                         &    ! actual message size
      ierr                                                     ! local error code

INTEGER(KIND=gapi_int)       ::  &
      ilength_ga

LOGICAL        ::           &
      lprefetch,            &    ! Flag. True, if only buffer should be filled (no data returned)
      ldirect_copy               ! Flag. True, if local copy (ie. no MPI send/recv required)

!------------------------------------------------------------------------------

  IF (PRESENT(ifield)  .OR.  (is_io_pe .AND. (io_config == 2))) THEN
    lprefetch = .FALSE.
  ELSE

    lprefetch = .TRUE.
    IF (lprefetching_available) THEN
      ! immediately return, if read buffer is already prepared
      IF (MPE_DBG_LEVEL > 5) THEN
        WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',            &
            'prefetching: nothing to do.'
      ENDIF
      RETURN
    ELSE
      IF (MPE_DBG_LEVEL > 5) THEN
        WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',            &
            'starting prefetching'
      ENDIF
      lprefetching_started = .TRUE.
    END IF
    
  END IF

  imaxlen = datatype_to_byte(BUFF_length)    ! length of buffer in bytes

  IF (MPE_DBG_LEVEL > 5) THEN
    WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',            &
        'into mpe_io_read '
  ENDIF

  ilength    = -1_gapi_int
  ilengthf   = -1_gribf_kind
  ilengthc   = -1_gribc_kind
  ilength_ga = -1_gapi_int
  ifirst_zero_record = -1

  !------------------------------------------------------------------------------
  ! --- IO process: read data from file, send it to compute PEs

  IO0 : IF ((me_io == 0) .AND. .NOT. lprefetching_available) THEN
    DO icomp_pes = 0, (num_compute-1)
      ldirect_copy = (icomp_pes == 0) .AND. (io_config == 1)

      ! this wait must not be within the IF (ilength > 0) block, because the 
      ! former message could be damaged before it is send otherwise
      CALL MPI_WAIT  (isend_requests(isend_current_request), &
                      istatus, ierr) ; CALL CHKERR(ierr, ierror)

      IF (ilength /= 0) THEN
        
        ! get a field entry from input file:
!US     CALL MPI_WAIT  (isend_requests(isend_current_request), &
!US                     istatus, ierr) ; CALL CHKERR(ierr, ierror)
        imaxlenc  = INT (imaxlen, gribc_kind)
        nudatc    = INT (nudat  , gribc_kind)

        SELECT CASE (current_lib)
#ifdef GRIBDWD
        CASE ('grb1')
          ! call to DWD grib library: read in grib file
          CALL cuegin (nudatc, imaxlenc, isend_buffer(1:,isend_current_request), ilengthc, ierrc)
          ilengthf = INT (ilengthc, gribf_kind) ! ilength is given in bytes!
#endif
#ifdef GRIBAPI
        CASE ('apix')
          ! call to grib-api library: read in grib file
          ilength_ga  = INT (imaxlen, gapi_int)
          CALL grib_read_from_file (nudatc, isend_buffer(1:,isend_current_request), ilength_ga, ierrc)

          IF ( ierrc /= GRIB_SUCCESS ) THEN
            IF (ierrc /= GRIB_END_OF_FILE) THEN
              CALL grib_check(ierrc, 'read_from_file', '')
              PRINT *, ' *** ERROR:  in grib_read_from_file'
            ELSE
              ierrc      = 0_gribc_kind
              ilength_ga = 0_gapi_int
              ilengthc   = 0_gribc_kind
              ilengthf   = 0_gribf_kind
            ENDIF
          ELSE
            IF (ilength_ga < INT( HUGE (1_gribf_kind), gapi_int) ) THEN
              ilengthf = INT (ilength_ga, gribf_kind) ! ilength is given in bytes!
            ELSE
              ! now we are in trouble: the message is too big:
              PRINT *, ' *** ERROR in mpe_io_read:  Message length in bytes is too big for standard integer ***'
              RETURN
              ierror = 10 
            ENDIF
          ENDIF
#endif
        CASE DEFAULT
          PRINT *, ' *** ERROR: wrong value for used library: ', current_lib
        END SELECT
        CALL CHKERR(ierrc, ierror)

        ! now set also ilength
        ilength = INT (ilengthf, gapi_int)

        IF (MPE_DBG_LEVEL > 10) THEN
          WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',            &
              'grib read done, len = ', ilength, imaxlenc, ';  isend_current_request = ', isend_current_request
        ENDIF
      ENDIF    ! ilength /= 0

      ! mark field as READ_END if end of file has been reached.
      IF (ilengthf > 0) THEN
        isend_buffer(MPE_FLAG, isend_current_request) = ilengthf
      ELSE
        isend_buffer(MPE_FLAG, isend_current_request) = MPE_FLAG_READ_END
        IF (ifirst_zero_record == -1) &
          ifirst_zero_record = icomp_pes
      END IF

      dat_length = byte_to_datatype(ilengthf)        
      IF (.NOT. ldirect_copy) THEN
        IF (MPE_DBG_LEVEL > 10) THEN
          WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',                  &
              'send data(tag= ', field_snd, ') to me_compute = ', icomp_pes, &
              'len = ', ilengthf, dat_length, dat_length+ISAFETY_MARGIN+1
        ENDIF
        
        ! because of the little endian problem, send some additional bytes
        CALL MPI_ISEND (isend_buffer(:,isend_current_request),             &
                        dat_length+ISAFETY_MARGIN+1,                       &
                        MPI_DATATYPE, icomp_pes, field_snd, MPE_COMM_READ, & 
                        isend_requests(isend_current_request), ierr)
        CALL CHKERR(ierr, ierror)
        field_snd = field_snd + 1
      ELSEIF (ilengthf > 0) THEN
        ! if the target PE is the local process, just insert record
        ! into internal read_buffer
        IF (MPE_DBG_LEVEL > 10) THEN
          WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',                  &
              'keep data(tag= ', field_snd, ') in me_compute = ', icomp_pes, &
              'len = ', ilengthf
        ENDIF

        CALL buffer_list_insert(read_buffer,                                       &
                                isend_buffer(1:dat_length, isend_current_request), &
                                ilengthf, field_snd, ierror)
        field_snd = field_snd + 1
      ELSE
        field_snd = field_snd + 1
        read_buffer%completed = .TRUE.
      ENDIF

      isend_current_request = isend_current_request + 1
      IF (isend_current_request > MAX_SEND_REQUESTS) isend_current_request = 1
    ENDDO

  ENDIF IO0

  !------------------------------------------------------------------------------
  ! --- Compute process: wait for new input from IO PE

  IF_COMPUTE : IF (is_compute_pe) THEN

    ldirect_copy = (me_compute == 0) .AND. (io_config == 1)

    RECEIVE : IF ((.NOT. ASSOCIATED(read_buffer%head)  .OR.  lprefetch)      &
      &            .AND.  .NOT. ldirect_copy) THEN      

      IF (MPE_DBG_LEVEL > 10) THEN
        WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',                  &
             'waiting for a message:  ', iread_recv_request, ALLOCATED(iread_recv_buffer)
      ENDIF

      ! if no read data is available: wait until it arrives
      CALL MPI_WAIT(iread_recv_request, istatus, ierr)            ; CALL CHKERR(ierr, ierror)
      CALL MPI_GET_COUNT(istatus, MPI_DATATYPE, dat_length, ierr) ; CALL CHKERR(ierr, ierror)

      ! because of the little endian problems, some additional bytes were sent
      IF (MPE_DBG_LEVEL > 10) THEN
        WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',                  &
             dat_length , ' INTs received (tag= ', istatus(MPI_TAG) , ') on compute node ', &
             me_compute, ' with FLAG = ', iread_recv_buffer(MPE_FLAG)
      ENDIF

      ! insert record into internal read_buffer
      dat_length  = dat_length - ISAFETY_MARGIN - 1
      ilengthf    = dat_length
      IF (iread_recv_buffer(MPE_FLAG) /= MPE_FLAG_READ_END) THEN
        byte_length = iread_recv_buffer(MPE_FLAG)
      ELSE
        byte_length = 0
      END IF

      CALL buffer_list_insert(read_buffer, iread_recv_buffer(1:),  &
                              byte_length, istatus(MPI_TAG), ierror)

      ! launch new IRECV() call, if this was not the final message
      IF (iread_recv_buffer(MPE_FLAG) /= MPE_FLAG_READ_END) THEN

        IF (io_config == 2) THEN
          io_pe = num_compute
        ELSE
          io_pe = 0
        END IF

        IF (io_pe /= me_io) THEN
          IF (MPE_DBG_LEVEL > 10) THEN
            WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',                  &
                 'Launch IRECV(tag= ', field_rcv, ') on me_compute = ', me_compute, &
                 '; IO source: ', io_pe
          ENDIF
          ! Note: including the MPE_FLAG, we receive (BUFF_length+1) values
          CALL MPI_IRECV(iread_recv_buffer(:), (BUFF_length+1),            &
                         MPI_DATATYPE, io_pe, field_rcv, MPE_COMM_READ,    &
                         iread_recv_request, ierr) ; CALL CHKERR(ierr, ierror)
          field_rcv = field_rcv + num_compute
        END IF

      ELSE
        read_buffer%completed = .TRUE.
      END IF

    END IF RECEIVE

    !------------------------------------------------------------------------------
    ! --- Read mode (not prefetching): return data from buffer
    
    IF (.NOT. lprefetch .AND. PRESENT(ifield)) THEN      

      IF (MPE_DBG_LEVEL > 10) THEN
        WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',                  &
           'Retrieving field from buffer.'
      ENDIF
      
      ! There should be data available in the internal
      ! read_buffer, fetch this:
      IF (.NOT. ASSOCIATED(read_buffer%head)) THEN
        ilengthf = 0
      ELSE
        ilengthf   = byte_to_datatype(read_buffer%head%field_len)
        dat_length = ilengthf
        DEBUG_WRITE(10) &
          "dat_length(tag=", read_buffer%head%field_tag, ") = ", dat_length
        ifield(1:dat_length) = read_buffer%head%field(1:dat_length)
        CALL buffer_list_remove_head(read_buffer)
      END IF

    END IF
      
  ENDIF IF_COMPUTE

  !------------------------------------------------------------------------------
  ! --- Finalize procedure

  ! if we have encountered a zero record (which finishes the read
  ! processes) there are possibly IRECVs triggered which must be
  ! satisfied:
  EOF : IF ((me_io == 0) .AND. (ifirst_zero_record /= -1)) THEN
    IF (MPE_DBG_LEVEL > 10) THEN
      WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',                  &
         'Satisfying obsolete pending IRECVs'
    ENDIF
    DO icomp_pes = 0, (ifirst_zero_record-1)
      ldirect_copy = (icomp_pes == 0) .AND. (io_config == 1)        

      CALL MPI_WAIT  (isend_requests(isend_current_request), &
        &             istatus, ierr) ; CALL CHKERR(ierr, ierror)
      isend_buffer(MPE_FLAG, isend_current_request) = MPE_FLAG_READ_END
      
      IF (.NOT. ldirect_copy) THEN
        DEBUG_WRITE(5) &
          &   "send data(tag=", field_snd, ") to me_compute = ", icomp_pes
        ! because of the little endian problem, send some additional bytes
        CALL MPI_ISEND (isend_buffer(:,isend_current_request),             &
          &             0+ISAFETY_MARGIN+1,                                &
          &             MPI_DATATYPE, icomp_pes, field_snd, MPE_COMM_READ, &
          &             isend_requests(isend_current_request), ierr) ; CALL CHKERR(ierr, ierror)
      END IF
      field_snd = field_snd + 1
      
      isend_current_request = isend_current_request + 1
      IF (isend_current_request > MAX_SEND_REQUESTS) isend_current_request = 1
    END DO
  END IF EOF

  ! Technically, the following WAITALL() statement is not
  ! required. However, it should have some performance benefits
  IF (is_io_pe .AND. (io_config == 2)) THEN
    IF (MPE_DBG_LEVEL > 10) THEN
      WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',                  &
         'Waiting for pending ISENDs'
    ENDIF
    CALL MPI_WAITALL(MAX_SEND_REQUESTS, isend_requests(:), &
      &              istatus_list, ierr) ; CALL CHKERR(ierr, ierror)
    IF (MPE_DBG_LEVEL > 10) THEN
      WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',                  &
         '...done'
    ENDIF
  END IF

  IF (ilengthf == 0) THEN
    IF (MPE_DBG_LEVEL > 10) THEN
      WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',                  &
         'Prefetching process completed.'
    ENDIF
    lprefetching_available = .TRUE.
    lprefetching_started   = .FALSE.
  END IF

  ilength = INT (ilengthf, gapi_int)

END SUBROUTINE mpe_io_read

!==============================================================================
!==============================================================================
!+ Send data from compute PEs to I/O PEs and write it on I/O PEs
!------------------------------------------------------------------------------

SUBROUTINE mpe_io_write (field, ilen_send, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   mpe_io_write is called by all compute- and IO-processors. 
!   -  Compute PEs: call event during grib_write()
!   -  Pure IO PEs: call event during mpe_io_node() cycling
!
!   For each contributing compute PEs the IO process always has
!   a pending IRECV(...) call. If a MPI_TESTSOME() returns the
!   result, that this call has finally arrived, the message is
!   copied to a list buffer and a new IRECV(...) is launched,
!   as long as the last received message did not contain an
!   MPE_FLAG_WRITE_END.
!
!------------------------------------------------------------------------------

! Parameters
  
! ilen_send: number of Bytes to be written  (length of GRIB record in bytes)
INTEGER(KIND=gribf_kind), INTENT(IN), OPTIONAL    ::    field(:)    ! array to be written  
INTEGER(KIND=gapi_int)  , INTENT(IN), OPTIONAL    ::    ilen_send
INTEGER                 , INTENT(INOUT)           ::    ierror      ! for error handling

! Local parameters
INTEGER  :: &
            pe_target, tag_target,                    &     ! pe and tag for target PE
            i,                                        &
            ilen_buff,                                &     ! number of INTs to be written
            imsg_length, istatus(MPI_STATUS_SIZE),    &
            compute_pe, tag_source, ioutcount,        &
            ioutstatus(MPI_STATUS_SIZE, num_compute), &
            ioutindex(num_compute),                   &
            dat_msg_length,                           &
            ierr,                                     &     ! local error code
            req_index

INTEGER(KIND=gribf_kind)  :: ilen_sendf        ! length in bytes, but only 4-byte int
!------------------------------------------------------------------------------

  ! Initializations

  ! compute processors that have to send data get
  ! value for ilen from application program
  IF (PRESENT(ilen_send)) THEN
    IF (ilen_send > INT(HUGE(1_gribf_kind),gapi_int) ) THEN
      ! now we are in trouble: the message is too big:
      PRINT *, ' *** ERROR in mpe_io_write:  Message length in bytes is too big for standard integer ***'
      RETURN
    ELSE
      ilen_sendf   = INT(ilen_send, gribf_kind)
      ilen_buff    = byte_to_datatype(ilen_sendf)
    ENDIF
  ELSE
    ilen_buff    = 0
  END IF

  DEBUG_WRITE(30) &
    'into mpe_io_write'

  ! Send data to I/O PE
  IS_COMP: IF (is_compute_pe) THEN

    field_snd  = field_snd + num_compute
    IF (io_config == 1) THEN
      pe_target  = MOD(me_compute, num_io)
    ELSEIF (io_config == 2) THEN
      pe_target  = MOD(me_compute, num_io) + num_compute
    ENDIF
    tag_target = field_snd

    IF (pe_target /= me_io) THEN

      DEBUG_WRITE(10) &
        &   'Send field #',field_snd,' from ', me_compute, ' to ',   &
        &   pe_target, ' with ilen = ', ilen_buff, '; tag_target = ', tag_target

      CALL MPI_WAIT  (isend_requests(isend_current_request), &
        &             istatus, ierr) ; CALL CHKERR(ierr, ierror, 1)

      isend_buffer(MPE_FLAG, isend_current_request) = ilen_sendf

      IF (PRESENT(field)) THEN
        isend_buffer(1:ilen_buff, isend_current_request) = field(1:ilen_buff)
      ELSE
        CALL CHKERR(ierrcode=1, iglobalerr=ierror) ! Error: Missing data field
      END IF

      ! because of the little endian problem send some additional bytes
      dat_msg_length = ilen_buff + ISAFETY_MARGIN + 1
      CALL MPI_ISEND (isend_buffer(:,isend_current_request), dat_msg_length,           &
        &             MPI_DATATYPE, pe_target, tag_target, MPE_COMM_WRITE,             &
        &             isend_requests(isend_current_request), ierr) ; CALL CHKERR(ierr, ierror)

      isend_current_request = isend_current_request + 1
      IF (isend_current_request > MAX_SEND_REQUESTS) isend_current_request = 1

    ELSE

      ! directly push data into write buffer
      IF (PRESENT(field)) THEN
        CALL buffer_list_insert(write_buffer, field, ilen_sendf, tag_target, ierror)
      ELSE
        CALL CHKERR(ierrcode=1, iglobalerr=ierror) ! Error: Missing data field
      END IF
      CALL buffer_list_write_action(write_buffer, ierror)

    ENDIF

  ENDIF IS_COMP

  ! receive buffer
  IS_IO: IF (is_io_pe) THEN

    IF (current_speed == 'f' .OR. current_speed == 'b') THEN
      ! -- receive available fields from any source and put them into
      ! -- the field buffer
      
      ! test for any message in ireceive buffer
      CALL MPI_TESTSOME(MAX_RECV_REQUESTS, irecv_requests(:), &
        &               ioutcount, ioutindex(:), ioutstatus(:,:), ierr)
      CALL CHKERR(ierr, ierror)
      DEBUG_WRITE(30) &
        "MPI_TESTSOME: ", ioutcount, " completions."
      
      ! if no irecv requests completed, leave subroutine, otherwise
      ! check on all pending requests:
      DO_INCOMING : DO i=1, ioutcount
        compute_pe = ioutstatus(MPI_SOURCE, i)
        ! test if message has arrived only recently:
        IF_ARRIVEDRECENTLY : IF ( .NOT. write_buffer%completed(compute_pe)) THEN
  
          tag_source = ioutstatus(MPI_TAG, i)
          DEBUG_WRITE(10) &
            &   "Receive message from source=", ioutstatus(MPI_SOURCE,i), &
            &   "with tag=",tag_source
   
          ! determine message length
          CALL MPI_GET_COUNT (ioutstatus(:,i), MPI_DATATYPE, imsg_length, ierr)
          CALL CHKERR(ierr, ierror)
          ! because of the little endian problem we received some
          ! additional bytes (and, additionally, we have a status flag
          ! at position 0).
          DEBUG_WRITE(20) &
            &   "Message length: ", imsg_length, " byte(s)."
        
          ! check 0th array entry, it tells us, if we must issue a
          ! consecutive irecv:

          ! find irecv request corresponding to "compute_pe", throw an error if there is none.
          req_index = get_recv_request(recv_fifo, compute_pe, ierr)

          IF (irecv_buffer(MPE_FLAG,req_index) /= MPE_FLAG_WRITE_END) THEN

            ! insert new field into buffer
            CALL buffer_list_insert(write_buffer, irecv_buffer(1:,req_index), &
                                    irecv_buffer(MPE_FLAG,req_index), tag_source, ierror)
            CALL buffer_list_write_action(write_buffer, ierror)
  
            ! if this is flagged as a data field, then launch new irecv for next field
            ! put IRECV on the FIFO queue (and launch after DO_INCOMING)
            CALL push_recv_fifo(recv_fifo, icount=(BUFF_length+1), isource=compute_pe, &
                                itag=tag_source + num_compute, ierror=ierror)
          ELSE
            write_buffer%completed(compute_pe) = .TRUE.
          END IF
        END IF IF_ARRIVEDRECENTLY
      END DO DO_INCOMING

      ! the new MPI_IRevcs have to be launched here for all finished requests
      CALL launch_irecvs(recv_fifo, ierror)

    ELSEIF (current_speed == 's') THEN

      DO compute_pe = 0, num_compute-1

        IF ( (MOD(compute_pe, num_io) == me_io) .AND.   &
             (.NOT. (compute_pe == me_compute))) THEN

          ! find irecv request corresponding to "compute_pe", throw an error if there is none.
          req_index = get_recv_request(recv_fifo, compute_pe, ierr)

          ! wait for the message from compute_pe to arrive
          CALL MPI_WAIT  (irecv_requests(req_index), istatus, ierr)
          CALL CHKERR(ierr, ierror)

          tag_source = istatus(MPI_TAG)

          IF (MPE_DBG_LEVEL > 10) THEN
            WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',                  &
                 'Received message(tag= ', istatus(MPI_TAG), '; source= ',     &
                 istatus(MPI_SOURCE), '  from recv_request:  ', req_index,     &
                 irecv_buffer(MPE_FLAG, req_index)
          ENDIF

          ! check end of file
          IF (irecv_buffer(MPE_FLAG, req_index) /= MPE_FLAG_WRITE_END) THEN

            ! insert new field into buffer
            CALL buffer_list_insert(write_buffer, irecv_buffer(1:,req_index), &
                          irecv_buffer(MPE_FLAG,req_index), tag_source, ierror)
            CALL buffer_list_write_action(write_buffer, ierror)

            ! if this is flagged as a data field, then launch new irecv for next field
            ! put IRECV on the FIFO queue (and launch after IF)
            CALL push_recv_fifo(recv_fifo, icount=(BUFF_length+1), isource=compute_pe, &
                 &              itag=tag_source + num_compute, ierror=ierror)
          ELSE
            write_buffer%completed(compute_pe) = .TRUE.
          ENDIF

          ! the new MPI_IRevcs have to be launched here for all finished requests
          CALL launch_irecvs(recv_fifo, ierror)

          ! else case (this compute pe is also io-pe) is handled above in the IS_COMP part
        ENDIF

      ENDDO

    ENDIF   ! speed mode of writing

  ENDIF IS_IO

END SUBROUTINE mpe_io_write

!==============================================================================
!==============================================================================
!+ Inquiry function
!------------------------------------------------------------------------------

PURE LOGICAL FUNCTION mpe_is_compute()
  
  mpe_is_compute = is_compute_pe
  
END FUNCTION mpe_is_compute

!==============================================================================
!==============================================================================
!+ Inquiry function
!------------------------------------------------------------------------------

PURE LOGICAL FUNCTION mpe_is_io()
  
  mpe_is_io = is_io_pe
  
END FUNCTION mpe_is_io

!==============================================================================
!==============================================================================
!+ Utility function: convert context message to MPI sendable integer array.
!------------------------------------------------------------------------------

SUBROUTINE compress_mpe_context_msg(context_msg, imsg)

  ! Parameters
  TYPE(mpe_context_msg), INTENT(IN)   :: context_msg    ! context message
  INTEGER, INTENT(OUT)  :: imsg(SIZE_MPE_CONTEXT_MSG)   ! compressed context message

  imsg = TRANSFER(context_msg, imsg, SIZE_MPE_CONTEXT_MSG)

END SUBROUTINE compress_mpe_context_msg

!==============================================================================
!==============================================================================
!+ Utility function: decode context message from MPI sendable integer array.
!------------------------------------------------------------------------------

SUBROUTINE decode_mpe_context_msg(imsg, context_msg)

  ! Parameters
  INTEGER, INTENT(IN)  :: imsg(SIZE_MPE_CONTEXT_MSG)   ! compressed context message
  TYPE(mpe_context_msg), INTENT(OUT)   :: context_msg  ! context message

  context_msg = TRANSFER(imsg, context_msg)

END SUBROUTINE decode_mpe_context_msg

!==============================================================================
!==============================================================================
!+ Utility function: compute length of integer array required to store XX bytes
!------------------------------------------------------------------------------

FUNCTION byte_to_datatype(ibytes)

  ! Parameters
  INTEGER, INTENT(IN)  :: ibytes           ! no. of bytes
  INTEGER              :: byte_to_datatype ! result

  byte_to_datatype = (ibytes + gribf_kind - 1)/gribf_kind

END FUNCTION byte_to_datatype

!==============================================================================
!==============================================================================
!+ Utility function: compute length of byte array required to store XX integers
!------------------------------------------------------------------------------

FUNCTION datatype_to_byte(iints)

  ! Parameters
  INTEGER, INTENT(IN)  :: iints            ! no. of words
  INTEGER              :: datatype_to_byte ! result

  datatype_to_byte = iints*gribf_kind

END FUNCTION datatype_to_byte

!------------------------------------------------------------------------------
! Subroutines for FIFO queue of IRECV requests
!------------------------------------------------------------------------------

!==============================================================================
!==============================================================================
!+ Utility function: allocate data structure of size MAX_RECV_FIFO_SIZE
!------------------------------------------------------------------------------

SUBROUTINE init_recv_fifo(fifo, max_recv_fifo_size, ierror)
  TYPE(recv_fifo_type), INTENT(INOUT) :: fifo
  INTEGER             , INTENT(IN)    :: max_recv_fifo_size
  INTEGER             , INTENT(INOUT) :: ierror      ! for error handling
  ! local variables
  INTEGER :: ierr
  ! allocate data structure. Note that we have a zero-indexed array (to make
  ! modulo operations easier):
  ALLOCATE(fifo%contents(0:(max_recv_fifo_size-1)), STAT=ierr)
  CALL CHKERR(ierr, ierror)
  fifo%front = 0
  fifo%count = 0
END SUBROUTINE init_recv_fifo

!==============================================================================
!==============================================================================
!+ Utility function: deallocate FIFO queue data structure
!------------------------------------------------------------------------------
SUBROUTINE finalize_recv_fifo(fifo, ierror)
  TYPE(recv_fifo_type), INTENT(INOUT) :: fifo
  INTEGER             , INTENT(INOUT) :: ierror      ! for error handling
  ! local variables
  INTEGER :: ierr
  DEALLOCATE(fifo%contents, STAT=ierr)
  CALL CHKERR(ierr, ierror)
  fifo%front = 0
  fifo%count = 0
END SUBROUTINE finalize_recv_fifo

!==============================================================================
!==============================================================================
!+ Utility function: insert element into FIFO queue data structure
!------------------------------------------------------------------------------
SUBROUTINE push_recv_fifo(fifo, icount, isource, itag, ierror)
  TYPE(recv_fifo_type), INTENT(INOUT) :: fifo
  INTEGER,              INTENT(IN)    :: icount, isource, itag
  INTEGER             , INTENT(INOUT) :: ierror                 ! for error handling
  ! local variables
  INTEGER :: new_element, ierr

  IF (fifo%count >= SIZE(fifo%contents)) THEN
     WRITE (0,*) "PUSH on full queue."
     ierr = -1
     CALL CHKERR(ierr, ierror)
  END IF

  ! calculate index at which to put next element
  new_element = MOD( (fifo%front + fifo%count), SIZE(fifo%contents) )

  fifo%contents(new_element)%count  = icount
  fifo%contents(new_element)%source = isource
  fifo%contents(new_element)%tag    = itag
  fifo%count = fifo%count + 1
END SUBROUTINE push_recv_fifo

!==============================================================================
!==============================================================================
!+ Utility function: remove element from the front of the FIFO queue data structure
!------------------------------------------------------------------------------
SUBROUTINE pop_recv_fifo(fifo, icount, isource, itag, ierror)
  TYPE(recv_fifo_type), INTENT(INOUT) :: fifo
  INTEGER,              INTENT(OUT)   :: icount, isource, itag
  INTEGER             , INTENT(INOUT) :: ierror                 ! for error handling
  ! local variables
  INTEGER :: ierr

  IF (fifo%count <= 0) THEN
     WRITE (0,*) "POP on empty queue."
     ierr = -1
     CALL CHKERR(ierr, ierror)
  END IF
  icount  = fifo%contents(fifo%front)%count
  isource = fifo%contents(fifo%front)%source
  itag    = fifo%contents(fifo%front)%tag

  ! advance the index of the front, making sure it wraps around the array
  ! properly.
  fifo%front = fifo%front + 1
  fifo%front = MOD(fifo%front, SIZE(fifo%contents))
  fifo%count = fifo%count - 1
END SUBROUTINE pop_recv_fifo

!==============================================================================
!==============================================================================
!+ Utility function: return if FIFO queue is empty
!------------------------------------------------------------------------------
FUNCTION is_empty_fifo(fifo)
  LOGICAL :: is_empty_fifo
  TYPE(recv_fifo_type), INTENT(IN) :: fifo
  is_empty_fifo = (fifo%count == 0)
END FUNCTION is_empty_fifo

!==============================================================================
!==============================================================================
!+ Utility function: launch IRECV requests until FIFO queue is empty or
!+                   MAX_RECV_REQUESTS requests are pending.
!------------------------------------------------------------------------------
SUBROUTINE launch_irecvs(fifo, ierror)
  TYPE(recv_fifo_type), INTENT(INOUT) :: fifo
  INTEGER             , INTENT(INOUT) :: ierror                 ! for error handling
  ! local variables
  INTEGER :: i, icount, isource, itag, ierr
  DO i=1,MAX_RECV_REQUESTS
     IF (is_empty_fifo(fifo)) EXIT
     IF (irecv_requests(i) == MPI_REQUEST_NULL) THEN
        ! pop request data from local FIFO queue
        CALL pop_recv_fifo(fifo, icount, isource, itag, ierror)
        ! launch a corresponding IRECV call
        CALL MPI_IRECV(irecv_buffer(:,i), icount,                      &
                       MPI_DATATYPE, isource, itag, MPE_COMM_WRITE,    &
                       irecv_requests(i), ierr)
        CALL CHKERR(ierr, ierror)
        irecv_isource(i) = isource
        IF (MPE_DBG_LEVEL > 10) THEN
          WRITE (*,*) 'proc [',me_compute,',',me_io,'] :  ',                 &
               "Launching new IRECV: source = ", isource,                    &
               ", tag = ", itag, "in index:  ", i, irecv_requests(i), MPI_REQUEST_NULL
        ENDIF
     END IF
  END DO
END SUBROUTINE launch_irecvs


!==============================================================================
!==============================================================================
!+ Utility function: find irecv request corresponding to "compute_pe", throw
!                    an error if there is none.
!------------------------------------------------------------------------------
FUNCTION get_recv_request(fifo, compute_pe, ierror)
  INTEGER :: get_recv_request
  TYPE(recv_fifo_type), INTENT(INOUT) :: fifo
  INTEGER,              INTENT(IN)    :: compute_pe
  INTEGER             , INTENT(INOUT) :: ierror                 ! for error handling
  ! local variables
  INTEGER :: i, ierr

  get_recv_request = MAX_RECV_REQUESTS+1
  get_req: DO i=1,MAX_RECV_REQUESTS
     IF (irecv_isource(i) == compute_pe) THEN
       get_recv_request = i
       EXIT get_req
     ENDIF
  END DO get_req
  IF (get_recv_request > MAX_RECV_REQUESTS) THEN
     WRITE (*,*) 'proc [',me_compute,',',me_io,'] :  ',                  &
                 "Requested IRECV request that has never been launched in ", me_io, ' for ', &
                  fifo%contents(compute_pe)%source, fifo%contents(compute_pe)%tag
     ierr = -1
     CALL CHKERR(ierr, ierror)
  END IF
END FUNCTION get_recv_request

!------------------------------------------------------------------------------
! Subroutines for linked list of field buffers
!------------------------------------------------------------------------------

!==============================================================================
!==============================================================================
!+ Utility function: Initialize linked list
!------------------------------------------------------------------------------

SUBROUTINE buffer_list_init(nudat, list, ierror)

! Parameters
INTEGER                         , INTENT(IN)    :: nudat        ! file descriptor
TYPE (buffer_type)              , INTENT(INOUT) :: list
INTEGER                         , INTENT(INOUT) :: ierror
! Local parameters
INTEGER                                         :: &
  ierrstat, icompute_pe

  IF (MPE_DBG_LEVEL > 5) THEN
    WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ', 'start buffer_list_init:  '
  ENDIF

  list%head => NULL()
  list%length         = 0
  list%next_work_tag0 = me_io
  list%next_work_tag  = list%next_work_tag0
  list%nudat          = nudat

  ! allocate flag array for receive completion
  ! for write buffers we have >1 contributing PEs
  IF (is_io_pe) THEN
    ALLOCATE(list%completed(0:(num_compute-1)), stat=ierrstat)
    CALL CHKERR(ierrstat, ierror, 1)
    ! set all those flags to TRUE which send to other IO PEs
    FORALL ( icompute_pe=0:(num_compute-1) )
      list%completed(icompute_pe) = (MOD (icompute_pe, num_io) /= me_io)
    END FORALL
    IF (MPE_DBG_LEVEL > 5) THEN
      WRITE (*,*) 'io-proc [',me_io,']:  ', 'allocated list%completed '
    ENDIF
  ELSEIF (is_compute_pe) THEN
    ALLOCATE(list%completed(0:0), stat=ierrstat)
    CALL CHKERR(ierrstat, ierror, 1)
    list%completed = .FALSE.
    IF (MPE_DBG_LEVEL > 5) THEN
      WRITE (*,*) 'compute-proc [',me_compute,']:  ', 'allocated list%completed '
    ENDIF
  ELSE
    CALL CHKERR(ierrcode=1, iglobalerr=ierror, ierrkey=1)
  END IF

  DEBUG_WRITE(10) &
    "step off buffer_list_init"
  
END SUBROUTINE buffer_list_init
 
!==============================================================================
!==============================================================================
!+ Utility function: Destroy linked list
!------------------------------------------------------------------------------

SUBROUTINE buffer_list_destroy(list, ierror)

! Parameters
TYPE (buffer_type)              , INTENT(INOUT) :: list
INTEGER                         , INTENT(INOUT) :: ierror
! Local parameters
TYPE (linked_list_type), POINTER                :: tmp
INTEGER                                         :: total_size

!------------------------------------------------------------------------------
  
  IF (MPE_DBG_LEVEL > 5) THEN
    WRITE (*,*) 'proc [',me_compute,',',me_io,']:  ',            &
        'buffer_list_destroy:  ', ASSOCIATED(list%completed)
  ENDIF

  IF (ASSOCIATED(list%head)) THEN
    total_size = 0
    DEBUG_WRITE(10) "Tags left:"
    tmp => list%head
    DO
      IF (.NOT. ASSOCIATED(tmp)) EXIT
      total_size = total_size + tmp%field_len
      DEBUG_WRITE(10) &
        "tag = ", tmp%field_tag
      tmp => tmp%next
    END DO
    IF (total_size > 0) THEN
      DEBUG_WRITE(10) "Tags left, total size: ", total_size
      CALL CHKERR(ierrcode=1, iglobalerr=ierror, ierrkey=1) ! Error: tags left
    END IF
  ELSE
    DEBUG_WRITE(10) "No tags left."
  END IF

  ! --- clean up
  DEBUG_WRITE(20) &
    "Clean up list buffer."
  DO
    IF (.NOT. ASSOCIATED(list%head)) EXIT
    ! remove list head
    tmp => list%head%next
    IF (list%head%istorage >= 0) THEN
      NULLIFY(list%head%field)     ! field was pre-allocated
      loccupied(list%head%istorage) = .FALSE.
    ELSE
      ! field was dynamically allocated
      DEALLOCATE(list%head%field)
    END IF
    DEALLOCATE(list%head)
    list%head => tmp
  END DO
  list%length = 0
  DEALLOCATE(list%completed)

END SUBROUTINE buffer_list_destroy

!==============================================================================
!==============================================================================
!+ Utility function: Insert data field into buffer list, ordered by field no.
!------------------------------------------------------------------------------

SUBROUTINE buffer_list_insert(list, field, field_len, field_tag, ierror)

! Parameters
TYPE (buffer_type)                     , INTENT(INOUT) :: list
INTEGER(KIND=gribf_kind)               , INTENT(IN)    :: field(:)
INTEGER                                , INTENT(IN)    :: field_len
INTEGER                                , INTENT(IN)    :: field_tag
INTEGER                                , INTENT(INOUT) :: ierror

! Local parameters
TYPE (linked_list_type), POINTER                       :: tmp, nb_left, element
INTEGER                                                :: errstat, int_length

!------------------------------------------------------------------------------

  ! find the correct position in list:
  tmp     => list%head
  nb_left => NULL()
  DO
    IF (.NOT. ASSOCIATED(tmp))      EXIT
    IF (tmp%field_tag >= field_tag) EXIT

    nb_left => tmp
    tmp     => tmp%next
  END DO

  ! insert element into linked list
  IF (ASSOCIATED(nb_left)) THEN
    ALLOCATE(nb_left%next, stat=errstat)
    CALL CHKERR(errstat, ierror, 1)
    element => nb_left%next
  ELSE
    ALLOCATE(list%head, stat=errstat)
    CALL CHKERR(errstat, ierror, 1)
    element => list%head
  END IF
  list%length = list%length + 1

  ! fill new list element with data
  element%field_tag = field_tag
  element%field_len = field_len
  element%next      => tmp
  
  ! reserve (if possible) a free entry in pre-allocated storage
  element%istorage = buffer_list_get_free_storage()
  IF (element%istorage >= 0) THEN
    ! using pre-allocated storage entry "element%istorage"
    loccupied(element%istorage) = .TRUE.
    element%field => field_storage(:, element%istorage)
  ELSE
    ALLOCATE(element%field(SIZE(field)) , stat=errstat)
    CALL CHKERR(errstat, ierror)
  END IF
  ! fill in field content
  int_length = byte_to_datatype(field_len)
  element%field(1:int_length)  = field(1:int_length)

  DEBUG_WRITE(20) &
    "Insert ", field_tag, " into buffer list."

END SUBROUTINE buffer_list_insert

!==============================================================================
!==============================================================================
!+ Utility function: Write out grib data, if possible.
!------------------------------------------------------------------------------

SUBROUTINE buffer_list_write_action(list, ierror)

! Parameters
TYPE (buffer_type), INTENT(INOUT)  ::  list
INTEGER           , intent(INOUT)  ::  ierror

! Local parameters
INTEGER                            ::  &
  ihead_tag
INTEGER  (KIND= gribc_kind )       ::  &
  nudatc, iwrite_lenc, ierrc           ! corresponding variables for C-routines
INTEGER  (KIND= gapi_int)          ::  &
  iwrite_lenga

!------------------------------------------------------------------------------

! @note :  As grib output may also contain records in
! arbitrary order, we can replace the following part by an
! "immediate flush policy"

!------------------------------------------------------------------------------

  ! --- write out as many records as possible (i.e. as there are
  ! --- available in the buffer list)

  WRITEACTION : DO

!US CALL buffer_list_probe_head(list, ihead_tag)
      
  IF (.NOT. ASSOCIATED(list%head)) THEN
    ihead_tag = -1
  ELSE
    ihead_tag = list%head%field_tag
  END IF

    ! -- if the field at the buffer head is due for write-out, call
    ! -- the DWD grib library
    IF     (current_speed == 'b') THEN
      ! records are eventually buffered and written in a sorted manner
      ! (as if done by a sequential program)
      ! this can lead to memory problems because of exhausting buffering
      IF (ihead_tag /= list%next_work_tag) EXIT WRITEACTION
    ELSEIF (current_speed == 'f') THEN
      ! this does not buffer, therefore grib records could be written
      ! in different order
      IF (ihead_tag == -1) EXIT WRITEACTION
    ELSEIF (current_speed == 's') THEN
      ! this does also not buffer, but records are received in ordered manner
      ! and are therefore also written in ordered manner
      IF (ihead_tag == -1) EXIT WRITEACTION
    ENDIF
    
    IF (list%head%field_len > 0) THEN
      ! "translate" C parameter
      iwrite_lenc  = INT (list%head%field_len,  gribc_kind) ! given in bytes
      nudatc       = INT (list%nudat, gribc_kind)
      ! call to DWD grib library: write out to grib file

        SELECT CASE (current_lib)
#ifdef GRIBDWD
        CASE ('grb1')
          ! call to DWD grib library: write grib file
          CALL cuegex (nudatc, list%head%field, iwrite_lenc, ierrc)
#endif
#ifdef GRIBAPI
        CASE ('api1','api2')
          ! call to grib-api library: write grib file
          ! Transfer integer-array back to character array
          iwrite_lenga = INT (list%head%field_len,  gapi_int)   ! given in bytes
          CALL grib_write_bytes (nudatc, list%head%field, iwrite_lenga, ierrc)
#endif
        CASE DEFAULT
          PRINT *, ' *** ERROR: wrong value for used library: ', current_lib
        END SELECT

      CALL CHKERR(ierrc, ierror, 1)
      
      DEBUG_WRITE(20) 'wrote grib'
    ELSE
      DEBUG_WRITE(20) &
        "skipped zero-size message."
    ENDIF

    CALL buffer_list_remove_head(list)

    ! increment field counter
    list%next_work_tag  = list%next_work_tag + num_io
    IF ((list%next_work_tag - list%next_work_tag0 + me_io) >= num_compute) THEN
      list%next_work_tag0 = list%next_work_tag0 + num_compute
      list%next_work_tag  = list%next_work_tag0
    END IF
    DEBUG_WRITE(10) &
      "next tag for write-out: ", list%next_work_tag

  END DO WRITEACTION

END SUBROUTINE buffer_list_write_action

!==============================================================================
!==============================================================================
!+ Utility function: wait and receive pending messages, perform write action.
!------------------------------------------------------------------------------

SUBROUTINE buffer_list_finalize_write(list, ierror)

! Parameters
TYPE (buffer_type)              , INTENT(INOUT) :: list
INTEGER                         , INTENT(INOUT) :: ierror
! Local parameters
INTEGER                                         ::   &
     ilen_recv, imsg_length, ioutcount,              &
     icompute_pe, i, ierr, tag_source, req_index

INTEGER :: ioutstatus(MPI_STATUS_SIZE, num_compute), &
     ioutindex(num_compute)

!------------------------------------------------------------------------------

  DEBUG_WRITE(10) &
    "into buffer_list_finalize_write"

  DEBUG_WRITE(20) &
    "Waiting for pending messages.", list%completed

  ! Note: The following loop may impede the compute PEs after
  ! the write process, since these are waiting in
  ! mpe_io_open(). Wouldn't it be possible to postpone the finalize
  ! until after the open process?

  ! --- loop until all contributing compute PEs have acknowledged completion
  DO_FINALPENDING : DO
    IF (ALL(list%completed)) EXIT DO_FINALPENDING

    CALL MPI_TESTSOME(MAX_RECV_REQUESTS, irecv_requests(:), &
                      ioutcount, ioutindex(:), ioutstatus(:,:), ierr)
    CALL CHKERR(ierr, ierror, 1)

    ! otherwise check on all pending requests:
    DO_INCOMING : DO i=1, ioutcount
      icompute_pe = ioutstatus(MPI_SOURCE, i)
      ! test if message has arrived only recently:
      IF_ARRIVEDRECENTLY : IF ( .NOT. list%completed(icompute_pe)) THEN

        ! find irecv request corresponding to "compute_pe", throw an error if there is none.
        req_index = get_recv_request(recv_fifo, icompute_pe, ierr)

        tag_source = ioutstatus(MPI_TAG, i)
        DEBUG_WRITE(10) &
               "Receive message from source=", ioutstatus(MPI_SOURCE,i), &
               "with tag=",tag_source
 
        ! determine message length
        CALL MPI_GET_COUNT (ioutstatus(:,i), MPI_DATATYPE, imsg_length, ierr)
        CALL CHKERR(ierr, ierror)
        ! because of the little endian problem we received some
        ! additional bytes (and, additionally, we have a status flag
        ! at position 0).
        ilen_recv = imsg_length - ISAFETY_MARGIN - 1

        ! check 0th array entry, if we must issue a consecutive irecv:
        IF (irecv_buffer(MPE_FLAG, req_index) /= MPE_FLAG_WRITE_END) THEN
          ! insert new field into buffer
          CALL buffer_list_insert(list, irecv_buffer(1:,i), &
                                  irecv_buffer(MPE_FLAG, i), tag_source, ierror)
          CALL buffer_list_write_action(list, ierror)

          ! if this is flagged as a data field, then put IRECV on the FIFO queue
          ! and launch new IRECV after the DO_INCOMING
          CALL push_recv_fifo(recv_fifo, icount=(BUFF_length+1), isource=icompute_pe, &
                              itag=tag_source + num_compute, ierror=ierror)
          DEBUG_WRITE(10) "launch_irecvs: buffer_list_finalize_write"
        ELSE
          list%completed(icompute_pe) = .TRUE.
          ! and launch new IRECV after the DO_INCOMING
        END IF
      END IF IF_ARRIVEDRECENTLY
    END DO DO_INCOMING

    ! launch new irecvs for completed messages
    CALL launch_irecvs(recv_fifo, ierror)

  END DO DO_FINALPENDING

!print *, ' call buffer_list_write_action from buffer_list_finalize_write 2'
! CALL buffer_list_write_action(list, ierror) ! paranoia
 
END SUBROUTINE buffer_list_finalize_write


!==============================================================================
!==============================================================================
!+ Utility function: get field tag of head item from linked list
!------------------------------------------------------------------------------

SUBROUTINE buffer_list_probe_head(list, ifield_tag)

! Parameters
TYPE (buffer_type)              , INTENT(IN)    :: list
INTEGER                         , INTENT(OUT)   :: ifield_tag

  IF (.NOT. ASSOCIATED(list%head)) THEN
    ifield_tag = -1
  ELSE
    ifield_tag = list%head%field_tag
  END IF

END SUBROUTINE buffer_list_probe_head

!==============================================================================
!==============================================================================
!+ Utility function: pop head item from linked list
!------------------------------------------------------------------------------

SUBROUTINE buffer_list_remove_head(list)

! Parameters
TYPE (buffer_type)              , INTENT(INOUT) :: list

! Local parameters
TYPE (linked_list_type), POINTER                :: tmp

  DEBUG_WRITE(10) &
    "buffer list: remove head ", list%head%field_tag, ", total list length=", list%length

  IF (.NOT. ASSOCIATED(list%head))  RETURN 
  
  ! remove list head
  tmp => list%head%next
  IF (list%head%istorage >= 0) THEN
    NULLIFY(list%head%field)     ! field was pre-allocated
    loccupied(list%head%istorage) = .FALSE.
  ELSE
    ! field was dynamically allocated
    DEALLOCATE(list%head%field)
  END IF
  DEALLOCATE(list%head)
  list%head => tmp    
  list%length = list%length - 1

END SUBROUTINE buffer_list_remove_head

!==============================================================================
!==============================================================================
!+ Utility function: returns next unoccupied storage index
!------------------------------------------------------------------------------

FUNCTION buffer_list_get_free_storage()

! Parameters
INTEGER :: buffer_list_get_free_storage ! result
INTEGER :: i, index

INTEGER, SAVE :: iloop_index = 0        ! note: SAVE attribute

  buffer_list_get_free_storage = -1 ! i.e. no storage available
  LOOP : DO i = 1,MAX_STORAGE_PREALLOC

    index = MOD(i+iloop_index, MAX_STORAGE_PREALLOC)

    IF (.NOT. loccupied(index)) THEN
      iloop_index                  = index
      buffer_list_get_free_storage = index
      EXIT LOOP
    END IF
  END DO LOOP

END FUNCTION buffer_list_get_free_storage

!==============================================================================
!==============================================================================
!+ Checks and merges error codes
!------------------------------------------------------------------------------

SUBROUTINE CHKERR(ierrcode, iglobalerr, ierrkey)

!------------------------------------------------------------------------------
!
! Description:
!   The purpose of this subroutine is to merge error codes
!   returned by different, successively executed statements. The
!   final, merged error code may then be returned to the caller of
!   the subroutine.
!
!------------------------------------------------------------------------------

! Parameters
INTEGER, INTENT(IN)           ::  ierrcode   ! local error results
INTEGER, INTENT(INOUT)        ::  iglobalerr ! "global" error (ie. for whole SR)
INTEGER, INTENT(IN), OPTIONAL ::  ierrkey    ! error key value (for identification)
! Local parameters
LOGICAL, PARAMETER    ::  lmsg_enabled = .TRUE. !<< enable flag for detailed messages
INTEGER, SAVE         ::  icurr_errnum = 1

  IF (PRESENT(ierrkey)) icurr_errnum = ierrkey

  ! Rule that merges error result with previous error state.
  ! Note: one could also take the MAX value or anything similar.
  IF (ierrcode /= 0) &
    iglobalerr = ierrcode

  ! usually, MPI libraries handle errors to be "fatal" anyway...
  IF (lmsg_enabled .AND. (ierrcode /= 0)) THEN
    DEBUG_WRITE(0) &
      "MODULE mpe_io2: Error! Key value = ", icurr_errnum
  END IF

  icurr_errnum = icurr_errnum + 1
  
END SUBROUTINE CHKERR

!================================================================================

END MODULE mpe_io2
