!==============================================================================
!+ Module for (parallel/asynchronous) NetCDF I/O
!------------------------------------------------------------------------------

MODULE netcdf_io

#ifdef NETCDF
!------------------------------------------------------------------------------
!
! Description:
!   This module provides routines and data structures for
!   netcdf I/O. In addition to routines to write netcdf
!   headers and variables into a netcdf file, it also
!   provides functionality for asynchronous I/O, where
!   dedicates PE collect the output data and have access
!   to output netcdf files.
!
! Current Code Owner: MeteoSwiss, Carlos Osuna
!  phone: +41 58 460 9877
!  fax:    
!  email:  carlos.osuna@meteoswiss.ch
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V4_25        2012/09/28 Carlos Osuna
!  Initial release
! V4_26        2012/12/06 Burkhardt Rockel, Hans-Juergen Panitz, Ulrich Schaettler
!  Changes in NetCDF I/O: change some variables of dimension with lenth 1 to scalar
!  and reduce the number of dimension IDs
!  Correction of restart related to asynchronous NetCDF I/O
!  In case that the 'c' file is written asynchronically the file name
!   will always contain the date of the very initialisation of a simulation,
!   also in case of a restart job
!  Renamed argument variable "namelist" in subroutines to "outblock" (because
!   namelist is a Fortran Keyword)
! V4_27        2013-03-19 Ulrich Schaettler, Hans-Juergen Panitz
!  Use nmsgchan from data_satellites (US)
!  Bug fix in calculation of list_out_nsteps (HJP)
!  now%nextstep has to be increased by one in order to get the correct output
!    times for the calculation of time_bnds (HJP)
! V4_28        2013/07/12 Ulrich Schaettler
!  Use subroutines and variables for vertical grid and reference atmospheres 
!    from module vgrid_refatm_utils
! V4_29        2013/10/04 Ulrich Schaettler
!  For the COSMO-Model only use vcoord and refatm from vgrid_refatm_utils
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
! V5_4         2016-03-10 Burkhardt Rockel, Oliver Fuhrer
!  Bug fix: after reading delta_t and h_scal these variables have not been 
!    distributed to the other processors in SR distribute_values_asynio
!  Updated code owner information
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

USE netcdf,           ONLY :   &
  nf90_def_dim,            &
  nf90_def_var,            &
  nf90_enddef,             &
  nf90_redef,              &
  nf90_put_att,            &
  nf90_put_var,            &
  nf90_noerr,              &
  nf90_strerror,           &
  NF90_create,             &
  NF90_clobber,            &
#ifndef MESSY
#ifdef PNETCDF
  NF90_netcdf4,            &
#endif
#endif
  NF90_set_fill,           &
  NF90_nofill,             &
  NF90_open,               &
  NF90_write,              &
  NF90_close,              & 
  NF90_CHAR,               &
  NF90_DOUBLE,             &
  NF90_FLOAT,              &
  NF90_GLOBAL,             &
  NF90_UNLIMITED

!------------------------------------------------------------------------------

USE environment,      ONLY :   &
  model_abort      ! aborts the program in case of errors

!------------------------------------------------------------------------------

USE data_parameters , ONLY :   &
  wp,        & ! KIND-type parameters for real variables
  sp,        & ! KIND-type parameters for single precision variables
  iintegers, & ! KIND-type parameter for standard integer variables
  irealgrib, & ! KIND-type parameter for the real variables in the grib library
  iwlength     ! length of a integer word of the grib library in byte

!------------------------------------------------------------------------------

USE parallel_utilities,  ONLY  : &
  distribute_values,   &
  distribute_vartab

!------------------------------------------------------------------------------

USE data_io,  ONLY : &
   pp_nl,                 & ! structure for gribout namelist 
   root,                  & ! pointer to the root of gribout namelists
   ydate_ini,             & ! start of the forecast
   ymode_write,           & ! mode for opening the (write) output files
   idims_id_out,          & ! array for the IDs of the dimensions of netCDF 
                            ! formatted output
   lbdclim,               & ! boundary data in climate model     ! PIK  (D.Hauffe)
   yncglob_institution,   & ! originating center name
   yncglob_title,         & ! title string for the output
   yncglob_source,        & ! program name and version
   yncglob_project_id,    & ! identification of the project of simulation
   yncglob_experiment_id, & ! identification of the experiment of simulation
   yncglob_contact,       & ! contact e.g. email address
   yncglob_references,    & ! URL, report etc.
   ncglob_realization,    & ! number of the realization of the experiment
   nzmxid,                & ! maximum number of NetCDF variabe IDs
   var,                   & ! array for LM variable table
   undefncdf,             & ! value for "undefined" in the netcdf routines
   ngribout,              & ! number of GRIBOUT namelist groups 
   num_gribtabs,          & ! number of GRIB tables used in LM variable table
   lmmss                    ! 10/14 digits date format

!------------------------------------------------------------------------------

USE data_parallel,  ONLY : & 
  num_compute,             & ! number of compute PEs
  nc_asyn_io,              & ! number of asynchronous I/O PEs (netcdf)
  num_asynio_comm,         & ! number of asynchronous I/O communicators (netcdf)
  num_iope_percomm,        & ! number of asynchronous I/O PE per communicator (netcdf)
  my_cart_id,              & ! rank of this subdomain in the cartesian communicator
  icomm_asynio,            & ! communicator for the group of I/O PEs (netcdf).
  icomm_cart,              & ! communicator that belongs to the cartesian grid
  intercomm_asynio,        & ! intercommunicator between compute and I/O groups
  imp_integers,            & ! determines the correct INTEGER type used in the model
                             ! for MPI
  imp_character,           & ! determines the correct CHARACTER type used in the
                             ! model for MPI
  imp_reals,               & ! determines the correct REAL type used in the model
                             ! for MPI
  nproc,                   & ! total number of processors: nprocx * nprocy + <nprocIO>
  my_world_id,             & ! rank of this subdomain in the global communicator
  lasync_io                  ! if .TRUE.: the model runs with extra PEs for
                             ! asynchronous IO

!------------------------------------------------------------------------------

USE data_runcontrol,  ONLY : &
  ldebug_io,               & ! if .TRUE., debug output for I/O
  lprintdeb_all,           & ! whether all or only one task prints debug output
  idbg_level,              & ! to control the verbosity of debug output
  itype_calendar,          & ! for specifying the calendar used
  nstart,                  & ! first time step of the forecast
  nstop,                   & ! last time step of the forecast
  lmulti_layer,            & ! run multi-layer soil model
  ntstep,                  & ! actual time step
  nfinalstop,              & ! last time step of the total forecast
  luse_rttov,              & ! if rttov-library is used
  isynsat_stat,            & ! status of synthetic satellite images
  ltime,                   & ! detailed timings of the program are given
  nhori,                   & ! number of sectors for the horizont array by the topographic
                             ! correction of the radiation
  lradtopo,                & ! if .TRUE., calculate topographic correction of radiation
  lmulti_snow                ! run multi-layer snow model

!------------------------------------------------------------------------------

USE data_modelconfig,   ONLY : &
  dt,                          & ! long time-step
  ke_soil,                     & ! number of layers in the multi-layer soil model
  ke,                          & ! number of grid points in vertical direction
  ke1,                         & ! KE+1
  ke_snow,                     & ! number of layers in multi-layer snow model
  ie_tot,                      & ! number of grid points in zonal direction
  je_tot,                      & ! number of grid points in meridional direction total
  czmls,                       & ! depth of the main soil layers in m
  czhls,                       & ! depth of the half soil layers in m
  pollat,                      & ! latitude of the rotated north pole (in degrees, N>0)
  pollon,                      & ! longitude of the rotated north pole (in degrees, E>0)
  polgam,                      & ! angle between the north poles of the systems
  startlon_tot,                & ! transformed longitude of the lower left grid point
  startlat_tot,                & ! transformed latitude of the lower left grid point
  dlat,                        & ! grid point distance in meridional direction (in degrees)
  dlon                           ! grid point distance in zonal direction (in degrees)

!------------------------------------------------------------------------------

USE data_fields, ONLY :  &
    hhl           ! geometrical height of model half levels

!------------------------------------------------------------------------------

#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
USE data_satellites, ONLY :  &
    nmsgchan      ! for output of MSG-variables
#endif

!------------------------------------------------------------------------------

USE utilities,         ONLY :  &
  get_utc_date,                &
  rlarot2rla,                  &
  phirot2phi

!------------------------------------------------------------------------------

USE io_utilities,      ONLY :  &
  make_fn, open_file                    

!------------------------------------------------------------------------------

USE vgrid_refatm_utils, ONLY :  &
  refatm, vcoord, nfltvc, svc1, svc2

!==============================================================================

IMPLICIT NONE

!==============================================================================

! include statements
INCLUDE "mpif.h"

!==============================================================================

INTEGER:: me_asynio,        & ! rank of I/O PE within I/O communicator
          iogroup,          & ! MPI group where me_asynio belongs to
          total_n_outsteps    ! total number of steps where output is requested

LOGICAL:: is_asynio_pe,   &   ! logical flag to determine if current PE is an I/O PE
          is_compute_pe       ! logical flag to determine if current PE is a compute PE

INTEGER (KIND=iintegers), ALLOCATABLE    ::  &
    list_out_nsteps(:)        ! list of steps where output is requested

CHARACTER (LEN=14), ALLOCATABLE          :: &
    list_yakdat(:)            ! list of dates for every output step 

INTEGER (KIND=iintegers), ALLOCATABLE ::   & 
    isend_request(:)          ! array of mpi request to send data to I/O PE
INTEGER (KIND=iintegers) :: &
    isend_flush_request       ! mpi request to send flush signal to I/O PE

INTEGER (KIND=iintegers) :: &
    ivar_id(nzmxid)           ! array of indexes for netcdf variables

INTEGER (KIND=iintegers), PARAMETER :: &
     nc_orgmdata_length=3     ! size of metadata related to netcdf file 
INTEGER (KIND=iintegers), PARAMETER :: &
     nc_varmdata_length=3     ! size of metadata related to field


INTEGER (KIND=iintegers),ALLOCATABLE :: &
     isend_buffer(:,:)        ! buffer for asynchronous communications

INTEGER (KIND=iintegers) :: &
     real2int_x               ! ratio in bytes between a real and integer

INTEGER (KIND=iintegers) :: &
     isend_current_request    ! current position in the sending buffer

INTEGER (KIND=iintegers) :: &
     ngribouts                ! number of gribout sections

INTEGER (KIND=iintegers) :: &
     nbits_ngribout           !position of bit were to stamp the gribout
                              !index in every MPI_TAG

INTEGER :: &
     impi_tag_ub              ! maximum mpi tag value allowed

TYPE :: buffers_container     ! structure that keeps blocks of isend buffers and list of request idx
    INTEGER(KIND=iintegers), POINTER :: buffer(:,:)
    INTEGER(KIND=iintegers), POINTER :: request_idx(:)
    TYPE( buffers_container ), POINTER :: next
    TYPE( buffers_container ), POINTER :: previous
END TYPE buffers_container

TYPE( buffers_container ), POINTER ::   & ! pointer to first block of isend buffers
   isend_buffers_root
TYPE( buffers_container ), POINTER ::   & ! pointer to current block of isend buffers 
   isend_buffers_now

INTEGER (KIND=iintegers)  :: block_capacity  ! number of elements of arrays in every isend block 

INTEGER (KIND=iintegers) ::             &
   num_blocks,                          & ! current number of isend blocks allocated 
   max_num_blocks                         ! maximum number of isend blocks allocatable

INTEGER (KIND=iintegers) :: izdebug ! verbosity level

!==============================================================================

CONTAINS

!==============================================================================
 
SUBROUTINE shutdown_netcdfio_sendbuffers() 

!------------------------------------------------------------------------------
!
! Description:
!  The routine closes the buffers used by compute PE to allocate data that was
!  send to a I/O PE. It checks if all the asynchronous MPI communications are
!  finilized and deallocates the buffers.
!
!------------------------------------------------------------------------------


INTEGER (KIND=iintegers) :: i
INTEGER (KIND=iintegers) :: izstatus(MPI_STATUS_SIZE), izerror
INTEGER (KIND=iintegers) :: ireq, i_block

isend_buffers_now => isend_buffers_root
DO i_block=1, num_blocks
  IF( izdebug > 10) THEN
    PRINT *,'Deallocating block ',i_block,my_world_id,num_blocks
  ENDIF
  DO ireq = 1,block_capacity
    IF( isend_buffers_now%request_idx(ireq) /= MPI_REQUEST_NULL) THEN
      CALL MPI_WAIT(isend_buffers_now%request_idx(ireq),izstatus,izerror)
    ENDIF
  ENDDO

  DEALLOCATE (isend_buffers_now%buffer )

  DEALLOCATE (isend_buffers_now%request_idx )

  IF( i_block /= num_blocks ) THEN
    isend_buffers_now => isend_buffers_now%next

    IF( ASSOCIATED( isend_buffers_now%previous ) ) THEN
     DEALLOCATE (isend_buffers_now%previous )
     ENDIF
  ENDIF

ENDDO


IF( ASSOCIATED( isend_buffers_now ) ) THEN
  IF (izdebug > 10) THEN
    PRINT *,'deallocate last sending block structure'
  ENDIF
  DEALLOCATE( isend_buffers_now )
ENDIF

END SUBROUTINE shutdown_netcdfio_sendbuffers

SUBROUTINE cdf_io_init( nproc, myw_id, yerrmsg, ierror )

!------------------------------------------------------------------------------
!
! Description:
!  The routine initializes the environment for a netcdf I/O PE
!  It creates MPI groups and communicators based on namelist options
!  for the I/O PE
!
!------------------------------------------------------------------------------
 
  INTEGER (KIND=iintegers), INTENT(IN)   ::  &
    nproc,     &    ! total number of processors
    myw_id          ! ID of this PE in igroup_world

  INTEGER (KIND=iintegers),  INTENT(OUT) :: &
    ierror          ! error returned by the routine

  CHARACTER (LEN=50),       INTENT(OUT)  ::       &
    yerrmsg         ! for MPI error message

  INTEGER (KIND=iintegers) :: &
    izmplcode       ! internal error returned by calls to MPI functions 
  INTEGER (KIND=iintegers)     ::     &
    nziope0,      & ! rank in MPI_COMM_WORLD of first I/O PE
    color_io        ! color to split communicators                 
                       
  INTEGER (KIND=iintegers) :: i

  LOGICAL                  :: attr_found


  ierror     = 0
  izmplcode  = 0
  yerrmsg    = '    '
  me_asynio = -1
  impi_tag_ub = 0
  num_blocks = 0
  isend_current_request = 1

  ! Initialize, whether additional debug output shall be done
  IF (ldebug_io) THEN
    IF (lprintdeb_all) THEN
      izdebug = idbg_level
    ELSE
      IF (my_cart_id == 0) THEN
        izdebug = idbg_level
      ELSE
        izdebug = 0
      ENDIF
    ENDIF
  ELSE
    izdebug = 0
  ENDIF
  
  CALL MPI_COMM_GET_ATTR(MPI_COMM_WORLD, MPI_TAG_UB, impi_tag_ub, attr_found, izmplcode)
  IF( izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_COMM_GET_ATTR extracting MPI_TAG_UB failed'
    RETURN
  ENDIF

  IF( attr_found .EQV. .FALSE.) THEN
    ierror  = 3335
    yerrmsg = 'MPI_TAG_UB_ not found'
    RETURN
  ENDIF

  IF (nc_asyn_io > 0 ) THEN
    nziope0 = nproc - nc_asyn_io ! first PE of IO-group 

    IF( myw_id >= nziope0) THEN     ! last PEs = I/O PE
      is_asynio_pe = .TRUE.
      is_compute_pe = .FALSE.
    ELSE
      is_asynio_pe = .FALSE.
      is_compute_pe = .TRUE.
    ENDIF


    ! create all the I/O communicators specified by num_iope_percomm
    IF( myw_id >= nziope0) THEN
      iogroup = (myw_id - nziope0) / num_iope_percomm
    ELSE
      iogroup = MPI_UNDEFINED
    ENDIF

    CALL MPI_COMM_SPLIT( MPI_COMM_WORLD, iogroup, 0, icomm_asynio, &
                     izmplcode)
    IF (izmplcode /= 0) THEN
      ierror  = izmplcode
      yerrmsg = 'MPI_COMM_SPLIT with icomm_asynio'
      RETURN
    ENDIF

    ! create an intercommunicator between compute PE and I/O PE
    IF( myw_id >= nziope0 .OR. myw_id == 0) THEN     
      color_io    =  0
    ELSE
      color_io    =  MPI_UNDEFINED
    ENDIF

    CALL MPI_COMM_SPLIT( MPI_COMM_WORLD, color_io, 0, intercomm_asynio, &
                     izmplcode)

    IF (izmplcode /= 0) THEN
      ierror  = izmplcode
      yerrmsg = 'MPI_COMM_SPLIT with intercomm_asynio'
      RETURN
    ENDIF

  ENDIF

  ! Get rank of current PE within the I/O communicator where it belongs to
  IF( is_asynio_pe ) THEN

    CALL MPI_COMM_RANK ( icomm_asynio, me_asynio,izmplcode) 
    IF (izmplcode /= 0) THEN
      ierror  = izmplcode
      yerrmsg = 'WRONG comm getting rank'
      RETURN
    ENDIF

  ENDIF



  real2int_x = BIT_SIZE( irealgrib ) / BIT_SIZE( iintegers )

END SUBROUTINE cdf_io_init

SUBROUTINE  allocate_io_sendbuffer( yerrmsg, ierr ) 

!------------------------------------------------------------------------------
!
! Description:
!  The routine allocates the buffers used by compute PE to buffer
!  data that is sent asynchronously to I/O PE 
!
!------------------------------------------------------------------------------

CHARACTER (LEN= *),       INTENT(OUT)  ::       &
    yerrmsg

INTEGER (KIND=iintegers),INTENT(OUT) :: ierr

INTEGER (KIND=iintegers) :: izerror

INTEGER (KIND=iintegers), POINTER :: buffer_block(:,:) 

INTEGER (KIND=iintegers), POINTER :: requestidx_block(:)


! the block capacity for isend buffers will be 10 MB
block_capacity = (10*1024*1024/8 ) / ( ie_tot*je_tot )
! maximum allocatable isend buffers = 1 GB
max_num_blocks =  ( (1000*1024*1024/8 ) / (ie_tot*je_tot ) )/ block_capacity
 

IF( block_capacity < 1 .OR. max_num_blocks < 1 ) THEN
  yerrmsg = 'can not buffer data for isend operations with less than 1 block or blocks sizes < 1'
  ierr = -1
  RETURN
ENDIF

ALLOCATE( isend_buffers_root,STAT=izerror )
ALLOCATE( isend_buffers_now, STAT=izerror )
ALLOCATE( buffer_block( (ie_tot*je_tot)*real2int_x +  nc_varmdata_length ,  block_capacity ), STAT=izerror )
ALLOCATE( requestidx_block ( block_capacity ) , STAT=izerror )

num_blocks = num_blocks + 1
buffer_block(:,:) = 0
requestidx_block(:) = MPI_REQUEST_NULL

NULLIFY( isend_buffers_root%next )
NULLIFY( isend_buffers_root%previous )
isend_buffers_root%buffer => buffer_block
isend_buffers_root%request_idx => requestidx_block
  
isend_buffers_now => isend_buffers_root

isend_flush_request = MPI_REQUEST_NULL

END SUBROUTINE allocate_io_sendbuffer


SUBROUTINE receive_io_data(outblock,istep,yextension, yerrmsg, ierror ) 

!------------------------------------------------------------------------------
!
! Description:
!  The routine is called by I/O PE once for every netcdf file that a I/O 
!  communicator has to write. It prepares the file, open it, write netcdf headers
!  and receives all the fields from the compute PE, until it receives the 
!  end of communication message from all the compute PE associated to a given I/O PE. 
!  Every message received is unpacked and written to the corresponding variable in 
!  the netcdf file.
!
!------------------------------------------------------------------------------

TYPE(pp_nl),       INTENT(IN) :: outblock
CHARACTER (LEN=1), INTENT(IN) :: yextension
INTEGER, INTENT(IN) :: istep

INTEGER (KIND=iintegers),  INTENT(OUT) :: ierror

CHARACTER (LEN= *),       INTENT(OUT)  ::       &
    yerrmsg


CHARACTER (LEN=250) :: yname

INTEGER (KIND=iintegers) :: izerror

CHARACTER (LEN=80) :: errmsg
INTEGER (KIND=iintegers) :: nuedat

INTEGER (KIND=iintegers) :: current_out_step

INTEGER (KIND=iintegers) :: nl_index_sh

INTEGER (KIND=iintegers) :: RECBUFF_len,data_len

INTEGER (KIND=iintegers) , ALLOCATABLE  :: recbuff(:)

REAL (KIND=irealgrib), ALLOCATABLE  :: rec_data(:)

INTEGER (KIND=iintegers) :: num_recvs, num_src_pes,ilen_recv, &
        var_idx,k_lev,tag_source, tag_source_end
INTEGER :: istatus(MPI_STATUS_SIZE),ostatus(MPI_STATUS_SIZE)

CHARACTER (LEN=3) :: yzhead

CHARACTER (LEN= 25)      :: yroutine

INTEGER (KIND=iintegers)  :: dummint

INTEGER (KIND=iintegers) :: i1,i2,i3,nyvar

INTEGER (KIND=iintegers),ALLOCATABLE :: ilist(:,:)

CHARACTER (LEN=14) :: stepDate

REAL (KIND=irealgrib) :: min,max,sum
INTEGER (KIND=iintegers) :: i,j

INTEGER (KIND=iintegers) :: yext_code

INTEGER (KIND=iintegers) :: klevels

  current_out_step =0

  yroutine = 'receive_io_data'

  num_recvs = 0
  ierror = 0
  yerrmsg = '     '

  num_src_pes= num_compute/num_iope_percomm
  IF( me_asynio < MOD(num_compute , num_iope_percomm ) ) THEN
    num_src_pes = num_src_pes + 1
  ENDIF


  RECBUFF_len = ( outblock%ie_out_tot * outblock%je_out_tot ) * real2int_x + nc_varmdata_length
  data_len = ( outblock%ie_out_tot * outblock%je_out_tot )

  ALLOCATE(recbuff(RECBUFF_len) , STAT=izerror)
  IF (izerror /= 0) THEN
    ierror = 7
    yerrmsg = 'can not allocate receiving buffer'
    RETURN
  ENDIF

  ALLOCATE(rec_data(data_len),  STAT=izerror)
  IF (izerror /= 0) THEN
    ierror = 7
    yerrmsg = 'can not allocate receiving buffer'
    RETURN
  ENDIF


  IF( outblock%lanalysis .AND. yextension /='c' ) THEN
    yzhead = 'la'//outblock%ydomain 
  ELSE
    yzhead = 'lf'//outblock%ydomain
  ENDIF

  IF( istep == 0 ) THEN
! istep=0 is only related to the output of the 'c' file
! Since this file always gets the date of the initialisation of a simulation
!  stepDate is set to ydate_ini
! Otherwise, in case of a restart job, the 'c' file would get the restart date
!   stepDate = list_yakdat(1)
    stepDate = ydate_ini
    current_out_step = list_out_nsteps(1)
  ELSE
    current_out_step = list_out_nsteps(istep)
    stepDate = list_yakdat(istep)
  ENDIF 

  CALL make_fn ( yzhead, stepDate, ydate_ini, outblock%ytunit, yextension,    &
                 current_out_step , dt, outblock%lhour, itype_calendar,       &
                 outblock%ydir, yname, lmmss, izdebug, izerror)

  yname = yname(1:LEN_TRIM(yname)) // '.nc'
 
  ! open netcdf file for writing
  CALL open_file (nuedat, yname, ymode_write, 'ncdf',  &
                  icomm_asynio, me_asynio, -1,lasync_io, 2, errmsg, izerror)
  IF (izerror /= 0) THEN
    ierror = 3322
    yerrmsg = errmsg
    RETURN
  ENDIF

  ! write netcdf headers
  CALL write_nc_gdefs (nuedat, outblock, icomm_asynio, nc_asyn_io,       &
                         yextension, istep, errmsg, izerror)

  IF (izerror /= 0) THEN
    ierror = 3323
    yerrmsg = errmsg
    RETURN
  ENDIF

  ! write variable definitions
  IF (yextension == 'c' ) THEN


    CALL write_nc_vdefs (nuedat, outblock%nyvar_c, outblock%ilist_c, ivar_id,  &
                           outblock%luvmasspoint, icomm_asynio, nc_asyn_io,  &
                           yextension, errmsg, izerror)
  ELSEIF (yextension == ' ') THEN

    CALL write_nc_vdefs (nuedat, outblock%nyvar_m, outblock%ilist_ml, ivar_id,  &
                     outblock%luvmasspoint, icomm_asynio, nc_asyn_io,  &
                     yextension, errmsg, izerror)
  ELSEIF (yextension == 'p') THEN

    CALL write_nc_vdefs (nuedat, outblock%nyvar_p, outblock%ilist_pl, ivar_id,  &
                     outblock%luvmasspoint, icomm_asynio, nc_asyn_io,  &
                     yextension, errmsg, izerror)
  ELSEIF (yextension == 'z') THEN

    CALL write_nc_vdefs (nuedat, outblock%nyvar_z, outblock%ilist_zl, ivar_id,  &
                     outblock%luvmasspoint, icomm_asynio, nc_asyn_io,  &
                     yextension, errmsg, izerror)
  ELSEIF (yextension == 's') THEN

    CALL write_nc_vdefs (nuedat, outblock%nyvar_s, outblock%ilist_sl, ivar_id,  &
                     outblock%luvmasspoint, icomm_asynio, nc_asyn_io,  &
                     yextension, errmsg, izerror)
  ENDIF

  IF (izerror /= 0) THEN
    ierror = 3324
    yerrmsg = errmsg
    RETURN
  ENDIF


  SELECT CASE (yextension )
    CASE('c')
      yext_code = 0
    CASE(' ')
      yext_code = 1
    CASE('p')
      yext_code = 2
    CASE('z')
      yext_code = 3
    CASE('s')
      yext_code = 4
    CASE DEFAULT
      yerrmsg = 'Invalid yextension '//yextension
      ierror = 3334
      RETURN
  END SELECT

  IF(ierror /= 0) THEN
    RETURN
  ENDIF

  nl_index_sh = ISHFT(outblock%nl_index, iwlength*8/2 )
  tag_source = IOR( nl_index_sh,IOR( ISHFT( istep,2), yext_code ) )

  IF( tag_source < 0 .OR. tag_source >= impi_tag_ub ) THEN
    yerrmsg = 'Trying to receive MPI communication with invalid tag'
    ierror= 3333
    RETURN
  ENDIF

  ! loop over expected messages with full levels from compute PE
  DO WHILE( num_recvs < num_src_pes)

    CALL MPI_PROBE (MPI_ANY_SOURCE, tag_source , MPI_COMM_WORLD, istatus, izerror)

    IF( izerror /= 0) THEN
      ierror = 3325
      yerrmsg = 'MPI_PROBE failed'
      RETURN
    ENDIF

    ! get the message size
    CALL MPI_GET_COUNT(istatus, MPI_INTEGER, ilen_recv,izerror )
    IF( izerror /= 0) THEN
      ierror = 3326

      yerrmsg = 'MPI_GET_COUNT failed'
      RETURN
    ENDIF

    ! receive an end of communication message (with size=1)
    IF ( ilen_recv == 1 ) THEN
      IF( izdebug > 10) THEN
        PRINT *,'RECV ',my_world_id, ' @ ', istatus(MPI_SOURCE), ' T ',istatus(MPI_TAG), ' S ', 1
      ENDIF
      CALL MPI_RECV( dummint, 1, MPI_INTEGER, istatus(MPI_SOURCE),  &
         istatus(MPI_TAG), MPI_COMM_WORLD, ostatus, izerror) 
      IF( izerror /= 0 ) THEN
        ierror = 3327
        yerrmsg = 'MPI_RECV failed'
        RETURN
      ENDIF

      num_recvs = num_recvs + 1

    ! receive data from a level of an output variable
    ELSE IF( ilen_recv > 1) THEN


      IF( ilen_recv /= RECBUFF_len ) THEN
        ierror = 3328
        yerrmsg = 'Length of message differs from expected'
        RETURN
      ENDIF
      IF( izdebug > 10) THEN
        PRINT *,'RECV ',my_world_id, ' @ ', istatus(MPI_SOURCE), ' T ',istatus(MPI_TAG), ' S ', ilen_recv, RECBUFF_len
      ENDIF
      CALL MPI_RECV( recbuff, RECBUFF_len,MPI_INTEGER,istatus(MPI_SOURCE), &
          istatus(MPI_TAG), MPI_COMM_WORLD, ostatus,izerror)

      IF( izerror /= 0 ) THEN
        ierror = 3329
        yerrmsg = 'MPI_RECV failed'
        RETURN
      ENDIF

      ! deserialize
      rec_data = transfer( recbuff, rec_data, SIZE=data_len)

      ! extract metadata packed in the message
      var_idx = INT( recbuff( ( outblock%ie_out_tot * outblock%je_out_tot )*real2int_x  + 1 ) )
      k_lev = INT( recbuff( ( outblock%ie_out_tot * outblock%je_out_tot )*real2int_x +2 ) )
      klevels = INT( recbuff( ( outblock%ie_out_tot * outblock%je_out_tot )*real2int_x +3 ) )


    
      ! write the level into the netcdf output file
      IF( yextension == 'c') THEN
        CALL write_netcdf_asyn( nuedat, rec_data, &
            outblock%ie_out_tot, outblock%je_out_tot, &
            outblock%ilist_c(1, var_idx ),outblock%ilist_c(2, var_idx ), &
            outblock%ilist_c(3, var_idx ), var_idx, k_lev, klevels, izerror,errmsg)       
      ELSEIF (yextension == ' ') THEN
        CALL write_netcdf_asyn( nuedat, rec_data, &
            outblock%ie_out_tot, outblock%je_out_tot, &
            outblock%ilist_ml(1, var_idx ),outblock%ilist_ml(2, var_idx ), &
            outblock%ilist_ml(3, var_idx ), var_idx, k_lev, klevels, izerror,errmsg)
      ELSEIF (yextension == 'p') THEN
       
        CALL write_netcdf_asyn( nuedat, rec_data, &
                outblock%ie_out_tot, outblock%je_out_tot, &
                outblock%ilist_pl(1, var_idx ),outblock%ilist_pl(2, var_idx ), &
                outblock%ilist_pl(3, var_idx ), var_idx, k_lev, klevels, izerror,errmsg)
      ELSEIF (yextension == 'z') THEN
        CALL write_netcdf_asyn( nuedat, rec_data, &
                outblock%ie_out_tot, outblock%je_out_tot, &
                outblock%ilist_zl(1, var_idx ),outblock%ilist_zl(2, var_idx ), &
                outblock%ilist_zl(3, var_idx ), var_idx, k_lev, klevels, izerror,errmsg)
      ELSEIF (yextension == 's') THEN
        CALL write_netcdf_asyn( nuedat, rec_data, &
                outblock%ie_out_tot, outblock%je_out_tot, &
                outblock%ilist_sl(1, var_idx ),outblock%ilist_sl(2, var_idx ), &
                outblock%ilist_sl(3, var_idx ), var_idx, k_lev, klevels, izerror,errmsg)
      ENDIF

      IF ( izerror /= 0 ) THEN
PRINT *,'COSUNA c',var_idx, k_lev,klevels
        ierror = 3330
        yerrmsg = errmsg
        RETURN
      ENDIF

    ENDIF
  ENDDO


  ! close the file
  CALL close_netcdf_file (nuedat,errmsg, izerror )

  ! reset
  CALL cleanup_env()

  IF ( izerror /= 0 ) THEN
    ierror = 3331
    yerrmsg = errmsg
    RETURN
  ENDIF

  DEALLOCATE( recbuff )
  DEALLOCATE( rec_data )
 

END SUBROUTINE receive_io_data

SUBROUTINE cleanup_env()
  idims_id_out(:)= -1

END SUBROUTINE cleanup_env

SUBROUTINE start_ionode ( yerrmsg, ierror )
!------------------------------------------------------------------------------
!
! Description:
!  The routine is called by I/O PE after initialization. 
!  It loops over all the time steps where output is expected, and call routines
!  to listen for messages with data from the compute PE.
!
!------------------------------------------------------------------------------

INTEGER (KIND=iintegers), INTENT(OUT) ::   &
  ierror

CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
  yerrmsg

CHARACTER (LEN=80)        :: yzerrmsg
INTEGER (KIND=iintegers)  :: izerror

TYPE (pp_nl), POINTER                 :: now
CHARACTER  (LEN=260)         :: yname_dummy

INTEGER (KIND=iintegers) :: i, istep
LOGICAL                  :: gribout_active

INTEGER (KIND=iintegers)  :: file_idx,gribout_cnt

ivar_id(:) = 0_iintegers

ierror   = 0
yerrmsg  = '         '

CALL distribute_values_asynio(izerror, yzerrmsg)
IF ( izerror /= 0_iintegers ) THEN
  ierror = izerror
  yerrmsg = yzerrmsg
  RETURN
ENDIF

CALL compute_outsteps(izerror, yzerrmsg )
IF ( izerror /= 0_iintegers ) THEN
  ierror = izerror
  yerrmsg = yzerrmsg
  RETURN
ENDIF

gribout_cnt = 0

now => root

gribout_loop_c: DO
  IF ( now%lwrite_const .EQV. .TRUE.) THEN
    IF ( MOD(gribout_cnt, num_asynio_comm) == iogroup ) THEN 
      CALL receive_io_data( now,0,'c',yzerrmsg, izerror )
      IF( izerror /= 0) THEN
        ierror = izerror
        yerrmsg = yzerrmsg
        RETURN
      ENDIF

    ENDIF
    gribout_cnt = gribout_cnt + 1

  ENDIF


  IF (ASSOCIATED(now%next) ) THEN
    now => now%next
  ELSE
    EXIT gribout_loop_c
  ENDIF

ENDDO gribout_loop_c

DO istep=1,total_n_outsteps

  now => root
 
  gribout_loop: DO

    gribout_active = .FALSE.

    i = 1
    DO WHILE ( i <=  now%outsteps .AND.    &
           now%ngrib(i) <= list_out_nsteps(istep) )
      
      IF( now%ngrib(i) == list_out_nsteps(istep) ) THEN
        gribout_active = .TRUE.
      ENDIF
      i=i+1
    ENDDO

    IF( gribout_active ) THEN

      file_idx = istep*ngribout + now%nl_index

      IF( MOD(gribout_cnt, num_asynio_comm) == iogroup ) THEN

        IF (luse_rttov .AND. now%nyvar_s > 0 .AND. isynsat_stat == 0) THEN
          CALL receive_io_data(now,istep,'s',yzerrmsg, izerror )
          IF( izerror /= 0) THEN
            ierror = izerror
            yerrmsg = yzerrmsg
            RETURN
          ENDIF

        ENDIF 

        IF (now%nyvar_p > 0 ) THEN
          CALL receive_io_data( now,istep,'p',yzerrmsg, izerror )
          IF( izerror /= 0) THEN
            ierror = izerror
            yerrmsg = yzerrmsg
            RETURN
          ENDIF

        ENDIF

        IF( now%nyvar_z > 0 ) THEN
          CALL receive_io_data( now,istep,'z', yzerrmsg, izerror )
          IF( izerror /= 0) THEN
            ierror = izerror
            yerrmsg = yzerrmsg
            RETURN
          ENDIF

        ENDIF

        IF( now%nyvar_m > 0 ) THEN
          CALL receive_io_data(now, istep,' ', yzerrmsg, izerror)
          IF( izerror /= 0) THEN
            ierror = izerror
            yerrmsg = yzerrmsg
            RETURN
          ENDIF

          ! now%nextstep has to be increased by one in order to
          ! get the correct output times for the calculation of time_bnds
          now%nextstep = MIN( now%nextstep + 1, now%outsteps)
          PRINT *, ' start_ionode:  Next asynchronous output will be done in step:  ', &
                     now%ngrib(now%nextstep), now%nextstep

        ENDIF

      ENDIF

      gribout_cnt = gribout_cnt + 1
    ENDIF

    IF (ASSOCIATED(now%next) ) THEN
      now => now%next
    ELSE
      EXIT gribout_loop
    ENDIF

  ENDDO gribout_loop


ENDDO


END SUBROUTINE start_ionode

SUBROUTINE write_netcdf_asyn ( nuedat, buff, i_tot, j_tot, &
             i1,i2,i3,var_idx, k_lev, klevels, ierror, yerrmsg )

!------------------------------------------------------------------------------
!
! Description:
!   write a level of a field into a netcdf files.
!
!------------------------------------------------------------------------------
!
! Subroutine arguments
INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
  nuedat,     & ! internal file descriptor
  i_tot,     & ! i-dimension for the array data
  j_tot,     & ! j-dimension for the array data
  i1,i2,i3,  & ! indices of variable in vartab table
  var_idx,   & ! index of variable
  k_lev,     & ! k level of data
  klevels      ! total number of k levels of field

INTEGER (KIND=iintegers), INTENT(OUT) ::   &
  ierror

CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
  yerrmsg      ! string for error messages
  


! Array arguments with intent(in):
REAL (KIND=irealgrib),     INTENT(IN)    ::  &
  buff(i_tot*j_tot)   
               ! array to be written


!------------------------------------------------------------------------------

! Local scalars
INTEGER  (KIND=iintegers)                ::  &
  ilen,          & ! i_tot*j_tot
  implcode,      & ! Error variable for MPI
  izerr            ! another error variable

INTEGER  (KIND=iintegers)                ::  &
  istatus(MPI_STATUS_SIZE)                     ! for status in MPI calls

REAL(KIND=irealgrib)                     ::  &
  rbuf(i_tot*j_tot), & ! buffer for receiving the data
  reof                 ! for sending end-of-file message


!- End of header
!------------------------------------------------------------------------------

izerr    = 0_iintegers
ierror   = 0
reof     = -999999.0_irealgrib
yerrmsg  = '         '

  IF ( klevels == 1) THEN
    ! this is a 2D variable
    izerr = nf90_put_var (nuedat, ivar_id( var_idx)  , buff,              &
          start=(/ 1, 1, 1 /), count=(/ i_tot, j_tot, 1 /))
    IF (izerr /= nf90_noerr) THEN
      yerrmsg = 'Error writing netcdf 2D variable'
      ierror  = 1
      RETURN
    ENDIF
  ELSEIF ( klevels > 1 ) THEN
    ! this is a 2D-slice of a 3D variable
    izerr = nf90_put_var (nuedat, ivar_id( var_idx), buff,             &
            start=(/ 1, 1, k_lev  ,1 /), count=(/ i_tot, j_tot, 1, 1 /))
    IF (izerr /= nf90_noerr) THEN
      yerrmsg = 'Error writing netcdf 3D variable'
      ierror  = 2
      RETURN
    ENDIF
  ELSE
    yerrmsg = 'Wrong value on klevels ' 
    ierror = 3
    RETURN 
  ENDIF


END SUBROUTINE write_netcdf_asyn

SUBROUTINE shutdown_io()

!------------------------------------------------------------------------------
!
! Description:
!   clean and deallocate I/O PE
!
!------------------------------------------------------------------------------

DEALLOCATE( list_out_nsteps )
DEALLOCATE( list_yakdat )
DEALLOCATE( czhls )
!US DEALLOCATE( vcoord )
DEALLOCATE( hhl )

END SUBROUTINE shutdown_io

SUBROUTINE compute_outsteps(izerror, yzerrmsg )

!------------------------------------------------------------------------------
!
! Description:
!   This routine computes in advance the index of the time steps where
!   output data is required. It also computes information needed at this time
!   steps, like the date.
!
!------------------------------------------------------------------------------

INTEGER (KIND=iintegers), INTENT(INOUT)     :: izerror
CHARACTER (LEN=80), INTENT(INOUT)           :: yzerrmsg


! Local scalars
TYPE (pp_nl), POINTER                 :: now


INTEGER (KIND=iintegers) :: istat, &
     nextmin, minamongvars,thismin,i,tstep

LOGICAL     ::  keeplooking 


CHARACTER (LEN=14) :: ydate1_l ! variable to keep date at certain time step
CHARACTER (LEN=28) :: ydate2_l
INTEGER (KIND=iintegers) :: nactday_l
REAL (KIND=wp)     :: acthour_l
INTEGER (KIND=iintegers ) :: max_val

nextmin = -1
total_n_outsteps =0
keeplooking = .TRUE.


! compute the number of gribout sections
now => root
ngribouts=0
comp_ngrib_loop: DO

  ngribouts = ngribouts + 1

  IF (ASSOCIATED(now%next) ) THEN
    now => now%next
  ELSE
    EXIT comp_ngrib_loop
  ENDIF

ENDDO comp_ngrib_loop

max_val = ngribouts
nbits_ngribout = 0
DO WHILE (max_val /= 0)
  nbits_ngribout = nbits_ngribout + 1
  max_val = ISHFT(max_val, -1)
ENDDO

! compute total number of steps where output is required
DO WHILE( keeplooking  )
  keeplooking = .FALSE.
  now => root
  minamongvars=-1 
  pregribout_loop: DO

    thismin=-1
    outstep_loop: DO i=1,now%outsteps
      IF ( now%ngrib(i) > nextmin .AND. &
         ( thismin .EQ. -1 .OR. now%ngrib(i) < thismin )) THEN
        thismin = now%ngrib(i)
        EXIT outstep_loop
      ENDIF         

    ENDDO outstep_loop

    IF( thismin .NE. -1 .AND.                      &
            ( thismin < minamongvars .OR. minamongvars .EQ. -1 )) THEN

      minamongvars = thismin
    ENDIF

    IF (ASSOCIATED(now%next) ) THEN
      now => now%next
    ELSE
      EXIT pregribout_loop
    ENDIF

  ENDDO pregribout_loop
  IF( minamongvars .NE. -1 .AND.    &
     ( (minamongvars < nstop - nstart) .OR. (minamongvars==nfinalstop - nstart) ) ) THEN

     keeplooking = .TRUE. 
     nextmin = minamongvars
     total_n_outsteps = total_n_outsteps +1

  ENDIF

ENDDO

! allocate arrays to keep information of each output step
ALLOCATE( list_out_nsteps( total_n_outsteps) , STAT=istat)
IF( istat /= 0) THEN
  izerror = 33344
  yzerrmsg = 'ERROR allocating list of output time steps'  
ENDIF
ALLOCATE( list_yakdat( total_n_outsteps) , STAT=istat)
IF( istat /= 0) THEN
  izerror = 33344
  yzerrmsg = 'ERROR allocating list of times of output time steps' 
ENDIF

nextmin = -1
tstep=0
keeplooking = .TRUE.

! loop over time steps
DO WHILE( keeplooking  )
  keeplooking = .FALSE.
  now => root
  minamongvars=-1 

  !loop over gribouts 
  gribout_loop2: DO

    thismin=-1
    outstep_loop2: DO i=1,now%outsteps
      IF ( now%ngrib(i) > nextmin .AND. &
         ( thismin .EQ. -1 .OR. now%ngrib(i) < thismin )) THEN
        thismin = now%ngrib(i)
        EXIT outstep_loop2
      ENDIF         

    ENDDO outstep_loop2

    IF( thismin .NE. -1 .AND.                      &
            ( thismin < minamongvars .OR. minamongvars .EQ. -1 )) THEN

      minamongvars = thismin
    ENDIF

    IF (ASSOCIATED(now%next) ) THEN
      now => now%next
    ELSE
      EXIT gribout_loop2
    ENDIF

  ENDDO gribout_loop2

  ! get date for next file (gribout) that will be written 
  IF( minamongvars .NE. -1 .AND.   & 
    ((minamongvars < nstop - nstart) .OR. (minamongvars==nfinalstop - nstart)) ) THEN
     keeplooking = .TRUE. 
    nextmin = minamongvars
    tstep = tstep + 1

    CALL get_utc_date(minamongvars+nstart, ydate_ini, dt, itype_calendar, ydate1_l,    &
                      ydate2_l, nactday_l, acthour_l)

    list_yakdat( tstep ) = ydate1_l
    ! nstart has to added to minamongvars
    ! otherwise the calculation of variables "time" and "time_bnds" goes wrong
    ! list_out_nsteps( tstep ) = minamongvars
    list_out_nsteps( tstep ) = minamongvars + nstart
  ENDIF

ENDDO


END SUBROUTINE compute_outsteps


SUBROUTINE send_asyn_io ( data, irec_len, my_orgdata_asyn,my_vardata_asyn, &
     outblock,ie_tot_l,je_tot_l, lflush,yerrmsg, ierror  )

!------------------------------------------------------------------------------
!
! Description:
!   This routine will be called from a compute PE to send, asynchronously, 
!   a full record to a given I/O PE. It is also used to send an end of message
!   (one integer size message) to the associated I/O PE to indicate that
!   the communication of data related to the current netcdf file has finished.
!
!------------------------------------------------------------------------------

! Subroutine arguments

INTEGER (KIND=iintegers), INTENT(IN)     ::  &
  irec_len                  ! size of message

INTEGER (KIND=iintegers), INTENT(IN)     ::  &
  ie_tot_l, je_tot_l        ! domain size of record being sent. 


REAL (KIND=irealgrib),     INTENT(IN)    ::  &
  data(ie_tot_l*je_tot_l)   ! buffer that contains the data 


INTEGER (KIND=iintegers), INTENT(IN)     ::  &
  my_orgdata_asyn( nc_orgmdata_length  ),  & ! array with metadata of netcdf file
  my_vardata_asyn( nc_varmdata_length  )     ! array with metadata of variable being sent 

LOGICAL, INTENT(IN)         :: lflush        ! flag to indicate end of communication

TYPE(pp_nl),              INTENT(IN)     ::    &
  outblock        ! pointer to the namelist group

CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
  yerrmsg        ! string for error messages

INTEGER  (KIND=iintegers), INTENT(OUT)   ::  &
  ierror         ! error status


! Local variables

REAL (KIND=irealgrib)   ::  &
  buffer(ie_tot_l*je_tot_l + nc_varmdata_length ) ! local buffer to pack data + metadata

INTEGER (KIND=iintegers) :: nl_index_sh  ! shifted index of a gribout in namelists

INTEGER (KIND=iintegers) ::  &
   tag_target, pe_target  ! mpi tag and pe for the target of the current message

INTEGER (KIND=iintegers) :: i,record_size 

INTEGER (KIND=iintegers) :: dummy_buf

INTEGER (KIND=iintegers) :: izstatus(MPI_STATUS_SIZE), izerror

INTEGER (KIND=iintegers) :: file_idx ! index for output netcdf file

INTEGER (KIND=iintegers) :: grb_idx  ! index of gribout section

CHARACTER (LEN=12)       :: errortag_str

CHARACTER (LEN=12)       :: impi_tag_ub_str

LOGICAL :: complete

INTEGER (KIND=iintegers), POINTER :: buffer_block(:,:)
INTEGER (KIND=iintegers), POINTER :: requestidx_block(:)
TYPE( buffers_container ) , POINTER :: buffer_element
INTEGER (KIND=iintegers) :: ireq

dummy_buf = 0


file_idx = my_orgdata_asyn(1)*ngribout + outblock%nl_index

grb_idx = my_orgdata_asyn(3)

nl_index_sh = ISHFT(outblock%nl_index, iwlength*8/2 )
tag_target = IOR( nl_index_sh, IOR( ISHFT( my_orgdata_asyn(1),2), my_orgdata_asyn(2) ) )

IF( tag_target < 0 .OR. tag_target >= impi_tag_ub ) THEN
  PRINT *,'MPI_TAG_UB ERROR ',impi_tag_ub
  write(unit=errortag_str, fmt='(I12)') tag_target
  write(unit=impi_tag_ub_str, fmt='(I12)') impi_tag_ub
  yerrmsg = 'Trying to establish MPI_ISEND communication with invalid tag :'//errortag_str//'. Max: '//impi_tag_ub_str
  ierror= 3332
  RETURN
ENDIF
 
pe_target = MOD( my_world_id,num_iope_percomm) + num_compute + &
      MOD( grb_idx,num_asynio_comm) * num_iope_percomm

 
record_size = irec_len * real2int_x  +nc_varmdata_length 

! prepare data + metadata in the buffer and send it asynchronously
IF ( irec_len > 0 ) THEN

  IF( isend_current_request > block_capacity) THEN
    
    isend_current_request = 1
    CALL MPI_TEST( isend_buffers_root%request_idx( block_capacity), complete, izstatus, izerror  ) 
    IF ( complete .EQV. .TRUE. ) THEN
      IF(izdebug > 5) THEN
        PRINT *,'block number ',num_blocks,' was filled, but the oldest one was released.'
      ENDIF

      DO ireq = 1,block_capacity
        CALL MPI_WAIT(isend_buffers_now%request_idx(ireq),izstatus,izerror)
      ENDDO
 
      isend_buffers_root%previous => isend_buffers_now
      isend_buffers_now%next => isend_buffers_root
      isend_buffers_now => isend_buffers_root
      isend_buffers_root => isend_buffers_root%next
    ELSE IF ( num_blocks < max_num_blocks ) THEN

      IF( izdebug > 5) THEN
        PRINT *, 'block number ',num_blocks,' was filled. Allocating a new block' 
      ENDIF
      ALLOCATE( buffer_block( (ie_tot*je_tot)*real2int_x +  nc_varmdata_length ,  block_capacity ), STAT=izerror )
      ALLOCATE( requestidx_block ( block_capacity ) , STAT=izerror )
      ALLOCATE( buffer_element, STAT=izerror )

      num_blocks = num_blocks + 1

      buffer_block(:,:) = 0
      requestidx_block(:) = MPI_REQUEST_NULL

      buffer_element%buffer => buffer_block
      buffer_element%request_idx => requestidx_block
      NULLIFY( buffer_element%next )
      buffer_element%previous => isend_buffers_now
      isend_buffers_now%next => buffer_element

      isend_buffers_now => buffer_element

    ELSE
      ierror  = 3336
      write(unit=errortag_str, fmt='(I12)') max_num_blocks
      yerrmsg = 'Reached maximum number of buffering blocks :'//errortag_str
      RETURN
    ENDIF

  ENDIF

  IF( isend_buffers_now%request_idx(isend_current_request) /= MPI_REQUEST_NULL ) THEN
    CALL MPI_WAIT( isend_buffers_now%request_idx(isend_current_request),izstatus,izerror)
  ENDIF

  isend_buffers_now%buffer(1:ie_tot_l*je_tot_l * real2int_x,isend_current_request) = &
            transfer( data(1:ie_tot_l*je_tot_l), isend_buffers_now%buffer(1,isend_current_request), &
            SIZE = ie_tot_l*je_tot_l * real2int_x )

  DO i=1,nc_varmdata_length
    isend_buffers_now%buffer(ie_tot_l*je_tot_l * real2int_x +i,isend_current_request) = my_vardata_asyn(i)
  ENDDO

  IF( izdebug > 10) THEN
    PRINT *,'ISEND ',my_world_id, ' @ ',pe_target, ' T ' , tag_target ,' S ' , record_size
  ENDIF
  CALL MPI_ISEND( isend_buffers_now%buffer (1, isend_current_request) ,record_size, MPI_INTEGER, &
    pe_target, tag_target, MPI_COMM_WORLD, isend_buffers_now%request_idx(isend_current_request), izerror) 
  IF( izerror /= 0 ) THEN
    yerrmsg = 'Error in MPI_ISEND record'
    ierror = 7
    RETURN
  ENDIF


  isend_current_request = isend_current_request + 1

ENDIF

! if flush flag, send an end of communications message
IF (lflush) THEN
  IF( izdebug > 10) THEN
    PRINT *,'ISEND ',my_world_id, ' @ ',pe_target, ' T ' , tag_target ,' S 1'
  ENDIF
  CALL MPI_ISEND(dummy_buf, 1,MPI_INTEGER, pe_target, tag_target, MPI_COMM_WORLD, &
        isend_flush_request, izerror )
  IF( izerror /= 0 ) THEN
    yerrmsg = 'Error in MPI_ISEND terminate record'
    ierror = 8
    RETURN
  ENDIF
ENDIF


END SUBROUTINE send_asyn_io


SUBROUTINE write_nc_gdefs (ncid, outblock, icomm, npes, yextension,      &
                           step, yerrmsg, ierror)
 
!------------------------------------------------------------------------------
!
! Description:
!   This routine initializes global definitions for a NetCDF output.
!   It writes the latitude and longitudes values and the vertical coordinates.
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
  TYPE(pp_nl),              INTENT(IN)     ::    &
    outblock         ! pointer to the namelist group

! Scalar arguments with intent(in):
  INTEGER (KIND=iintegers),   INTENT(IN)   ::    &
    ncid            ! NetCDF file IDr

  INTEGER (KIND=iintegers),   INTENT(IN)   ::    &
   icomm,    & ! MPI communicator
   npes        ! number of PEs

  CHARACTER (LEN=1),          INTENT(IN)   ::    &
    yextension    ! indicates model variables (''), p-('p') or z-levels ('z')

  INTEGER (KIND=iintegers),   INTENT(IN)   ::    &
    step

  CHARACTER (LEN=*),          INTENT(OUT)  ::    &
    yerrmsg     ! string for error messages

! Scalar arguments with intent(out):
  INTEGER (KIND=iintegers),   INTENT(OUT)  ::    &
    ierror          ! error index

!------------------------------------------------------------------------------
!
! Local scalars:

INTEGER (KIND=iintegers), PARAMETER  ::          &
  nbnds=2,   & !
  ntime=1,   & !
  nheight=1, & !
  nsynmsg=32

INTEGER (KIND=iintegers)  :: i, j, timestep, isect

INTEGER (KIND=iintegers)  :: &
    jgridVarID,     & ! NetCDF ID for rotated_pole
    jlonVarID,      & ! NetCDF ID for longitude
    jlatVarID,      & ! NetCDF ID for latitude
    jrlonVarID,     & ! NetCDF ID for rotated longitude
    jrlatVarID,     & ! NetCDF ID for rotated latitude
    jslonuVarID,    & ! NetCDF ID for U-direction shifted longitude
    jslatuVarID,    & ! NetCDF ID for U-direction shifted latititude
    jslonvVarID,    & ! NetCDF ID for V-direction shifted longitude
    jslatvVarID,    & ! NetCDF ID for V-direction shifted latitude
    jsrlonVarID,    & ! NetCDF ID for shifted rotated longitude
    jsrlatVarID,    & ! NetCDF ID for shifted rotated latitude
    jvcVarID,       & ! NetCDF ID for the vertical component
    jsectVarID,     & ! NetCDF ID for nhori (number of sectors for the horizon)
    jh2VarID,       & ! NetCDF ID for the 2m height
    jh10VarID,      & ! NetCDF ID for the 10m height
    jhtoaVarID,     & ! NetCDF ID for the TOA height
    jwbt13VarID,    & ! NetCDF ID for the 1.3 Celsius wet bulb temperature
    jsynmsgVarID,   & ! NetCDF ID for the synthetic satellite data (only: SYNMSG)
    jmsgchanVarID,  & ! NetCDF ID for the synthetic satellite data (grouped by channel)
    jsoilVarID,     & ! NetCDF ID for the multi soil layer component
    jsoilbdsID,     & ! NetCDF ID for the multi soil layer bounds
    jsnowVarID,     & ! NetCDF ID for the multi snow layer component !_br 23.01.12
    jtimeID,        & ! NetCDF ID for the time
    jtbdsID           ! NetCDF ID for the time bounds
    
CHARACTER (LEN=40) :: ydate

CHARACTER (LEN= 1) :: &
  cvctype           ! character version of ivctype

CHARACTER (LEN=19) :: &
  creation_date     ! actual time when the data file is written

!Local arrays:
REAL (KIND=wp)     :: &
  zlongitude(outblock%ie_out_tot,outblock%je_out_tot),  & ! geographic longitudes
  zlatitude (outblock%ie_out_tot,outblock%je_out_tot),  & ! geographic latitudes
  zrotlon   (outblock%ie_out_tot),                      & ! rotated longitudes
  zrotlat   (outblock%je_out_tot),                      & ! rotated latitudes
  zslonu    (outblock%ie_out_tot,outblock%je_out_tot),  & ! U-direction shifted longitudes
  zslatu    (outblock%ie_out_tot,outblock%je_out_tot),  & ! U-direction shifted latitudes
  zslonv    (outblock%ie_out_tot,outblock%je_out_tot),  & ! V-direction shifted longitudes
  zslatv    (outblock%ie_out_tot,outblock%je_out_tot),  & ! V-direction shifted latitudes
  zsrlon    (outblock%ie_out_tot),                      & ! shifted rotated longitudes
  zsrlat    (outblock%je_out_tot)                         ! shifted rotated latitudes

REAL (KIND=wp)     :: &
  time(ntime), time_bnds(nbnds, ntime)
    
REAL (KIND=irealgrib) :: &
  zsoil_bnds(nbnds, ke_soil+1), &
  zczhls(0:ke_soil+1)
                      
REAL (KIND=wp),     ALLOCATABLE   :: &
  zvcoord(:)        ! vertical coordinate
    
REAL (KIND=wp)     :: &
     zheight_2m(nheight),     & !  2m height
     zheight_10m(nheight),    & ! 10m height
     zheight_toa(nheight),    & ! TOA height (actually model TOA)
     zwbtemp_13c(nheight)       ! 1.3 Celsius wet bulb temperature

#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
REAL (KIND=wp)     :: &
     zmsgchan_wave(nmsgchan)    ! synthetic satellite channel wavelenghts
#endif

INTEGER :: &
    k                           ! loop index
REAL (KIND=irealgrib) :: &
    zind_snow(ke_snow)          ! Indices of snow layers

INTEGER (KIND=iintegers) :: &
  zsect(nhori),   & !
  ztime_values(8)   ! holds the date and time information taken from 
                    ! internal subroutine date_and_time

INTEGER (KIND=iintegers) :: my_comm_id, implcode

INTEGER (KIND=iintegers) :: zvcoord_len          

! string variable to hold grid information
CHARACTER (LEN=200)     grid_mapping

 
!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

zvcoord_len = 0

IF (npes > 1) THEN
 ! Get id in communicator comm
  CALL MPI_COMM_RANK(icomm,my_comm_id,implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'MPI_COMM_RANK failed'
    ierror  = 2
    RETURN
  ENDIF
ELSE
  my_comm_id = 0
ENDIF


! The next block actually defines dimensions and variables
! It is executed by all the IO processors in the asynchronous netcdf
! or just by PE 0 in the sequential netcdf
IF( lasync_io .OR. my_comm_id == 0 ) THEN

! some pre-settings
  ierror = 0
  IF( lasync_io ) THEN
    IF( step == 0 ) THEN
      timestep = list_out_nsteps(1)
    ELSE
      timestep = list_out_nsteps( step )
    ENDIF
  ELSE
    timestep = ntstep
  ENDIF
 
  time = 0._wp
  time_bnds = 0._wp
  

! determine the creation date
  CALL DATE_AND_TIME(values=ztime_values)
  WRITE (creation_date,'(I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2)') &
            ztime_values(1),'-',ztime_values(2),'-',ztime_values(3),' ', &
            ztime_values(5),':',ztime_values(6),':',ztime_values(7)
  
! write the global attributes
!_br 23.01.12
  IF (TRIM(yncglob_title) /= '-') THEN
    ierror=nf90_put_att(ncid, NF90_GLOBAL, "title",    TRIM(yncglob_title))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
  ENDIF
  IF (TRIM(yncglob_institution) /= '-') THEN
    ierror=nf90_put_att(ncid, NF90_GLOBAL, "institution",               &
                                            TRIM(yncglob_institution))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
  ENDIF
  IF (TRIM(yncglob_source) /= '-') THEN
    ierror=nf90_put_att(ncid, NF90_GLOBAL, "source",                    &
                                            TRIM(yncglob_source))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
  ENDIF
  IF (TRIM(yncglob_project_id) /= '-') THEN
    ierror=nf90_put_att(ncid, NF90_GLOBAL, "project_id",                &
                                            TRIM(yncglob_project_id))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
  ENDIF
  IF (TRIM(yncglob_experiment_id) /= '-') THEN
    ierror=nf90_put_att(ncid, NF90_GLOBAL, "experiment_id",             &
                                            TRIM(yncglob_experiment_id))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
  ENDIF
  IF (ncglob_realization /= -999) THEN
    ierror=nf90_put_att(ncid, NF90_GLOBAL, "realization",ncglob_realization)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
  ENDIF
  ierror=nf90_put_att(ncid, NF90_GLOBAL, "Conventions","CF-1.4")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, NF90_GLOBAL, "conventionsURL", &
                         "http://www.cfconventions.org/")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  IF (TRIM(yncglob_contact) /= '-') THEN
    ierror=nf90_put_att(ncid, NF90_GLOBAL, "contact",                  &
                                           TRIM(yncglob_contact))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
  ENDIF
  IF (TRIM(yncglob_references) /= '-') THEN
    ierror=nf90_put_att(ncid, NF90_GLOBAL,"references",                &
                                           TRIM(yncglob_references))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
  ENDIF

  ierror=nf90_put_att(ncid, NF90_GLOBAL, "creation_date", creation_date)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

! define dimensions of time and time_bounds
  IF (lasync_io ) THEN
#ifdef PNETCDF
 ! unlimited dimension currently buggy with PNETCDF 
    ierror=nf90_def_dim(ncid,"time",  1 , idims_id_out(5))
#else
    ierror=nf90_def_dim(ncid,"time",  NF90_UNLIMITED , idims_id_out(5))
#endif
  ELSE
    ierror=nf90_def_dim(ncid,"time",  NF90_UNLIMITED , idims_id_out(5))
  ENDIF
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_def_var(ncid, "time", NF90_DOUBLE, (/ idims_id_out(5) /), jtimeID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_att(ncid, jtimeID, "standard_name", "time")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_att(ncid, jtimeID, "long_name", "time")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_def_dim(ncid,"bnds", nbnds, idims_id_out(6))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_def_var(ncid, "time_bnds", NF90_DOUBLE,       &
                (/ idims_id_out(6), idims_id_out(5) /), jtbdsID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_att(ncid, jtbdsID, "long_name", "time bounds")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

! determine the actual forecast time

  IF (outblock%lanalysis) THEN
    ! When nudging is active the output fields are treated as analyses
    ! and the reference time is the actual forecast time
    ydate = 'seconds since '//list_yakdat( step)(1:4)//'-'//list_yakdat(step)(5:6)//'-'//     &
             list_yakdat(step)(7:8)//' '//list_yakdat(step)(9:10)//':'//list_yakdat(step)(11:12)//':'  &
             //list_yakdat(step)(13:14)
  ELSE
    ydate = 'seconds since '// ydate_ini(1:4)//'-'//ydate_ini(5:6)//'-'  &
             //ydate_ini(7:8)//' '//ydate_ini(9:10)//':'//ydate_ini(11:12)//':'&
             //ydate_ini(13:14)
  ENDIF
  ierror = nf90_put_att (ncid, jtimeID, "units", ydate)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror = nf90_put_att (ncid, jtbdsID, "units", ydate)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  IF     (itype_calendar == 0) THEN
    ierror = nf90_put_att (ncid, jtimeID, "calendar", "proleptic_gregorian")
  ELSEIF (itype_calendar == 1) THEN
    ierror = nf90_put_att (ncid, jtimeID, "calendar", "360_day")
  ELSEIF (itype_calendar == 2) THEN
    ierror = nf90_put_att (ncid, jtimeID, "calendar", "365_day")
  ENDIF
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror = nf90_put_att (ncid, jtimeID, "bounds", "time_bnds")

! the time will be in seconds

  IF( outblock%lanalysis ) THEN
    time(1)        = 0.0_wp
  ELSE
    time(1)        = timestep*dt
  ENDIF

!
!HJP Begin Change relevant for restart
!
! IF(lbdclim .AND. (outblock%nextstep /= 1)) THEN
!    time_bnds(1,1) = time(1) - &
!      (outblock%ngrib(outblock%nextstep) - outblock%ngrib(outblock%nextstep-1))*dt
! ELSE
!    time_bnds(1,1) = 0._wp
! ENDIF
!

  time_bnds(1,1) = 0._wp

  IF (lbdclim) THEN

     IF (outblock%nextstep /= 1) THEN
       time_bnds(1,1) = time(1) - &
           (outblock%ngrib(outblock%nextstep) - outblock%ngrib(outblock%nextstep-1))*dt
     ELSE
       IF (outblock%ngrib( outblock%nextstep) .GT. 0 ) THEN
         time_bnds(1,1) = time(1) -                             &
           (outblock%ngrib( outblock%nextstep+1 )-outblock%ngrib( outblock%nextstep ))*dt
       ENDIF
     ENDIF

  ENDIF
!
!HJP End Change relevant for restart
!
  time_bnds(2,1) = time(1)

! set the values of rotated North Pole in the grid_mapping attribute
  grid_mapping = 'rotated_pole'

  ierror=nf90_def_var(ncid, "rotated_pole", NF90_CHAR, jgridVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_att(ncid, jgridVarID, "long_name",                 &
                            "coordinates of the rotated North Pole")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  
  ierror=nf90_put_att(ncid, jgridVarID, "grid_mapping_name",                &
                            "rotated_latitude_longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  
  ierror=nf90_put_att(ncid, jgridVarID, "grid_north_pole_latitude",         &
                                        REAL(pollat,sp))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  
  ierror=nf90_put_att(ncid, jgridVarID, "grid_north_pole_longitude",        &
                                        REAL(pollon,sp))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  IF (polgam /= 0._wp) THEN

   ierror=nf90_put_att(ncid, jgridVarID, "north_pole_grid_longitude",      &
                                        REAL(polgam,sp))

    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

  ENDIF

! determine rotated lat/lon values

  DO i = 1, outblock%ie_out_tot
    zrotlon(i) = startlon_tot + (outblock%i_out_start-1)*dlon + (i-1)*dlon
  ENDDO
  DO j = 1, outblock%je_out_tot
    zrotlat(j) = startlat_tot + (outblock%j_out_start-1)*dlon + (j-1)*dlat
  ENDDO

! define the attributes of longitude and latitude
  ierror=nf90_def_dim(ncid,"rlon",outblock%ie_out_tot, idims_id_out(1))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_def_dim(ncid,"rlat",outblock%je_out_tot, idims_id_out(2))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ! rotated longitude
  ierror=nf90_def_var(ncid, "rlon", NF90_FLOAT, (/ idims_id_out(1) /),     &
                                                jrlonVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_att(ncid, jrlonVarID, "standard_name", "grid_longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jrlonVarID, "long_name", "rotated longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jrlonVarID, "units", "degrees")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ! rotated latitude
  ierror=nf90_def_var(ncid, "rlat", NF90_FLOAT, (/ idims_id_out(2) /),      &
                                                jrlatVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jrlatVarID, "standard_name", "grid_latitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jrlatVarID, "long_name", "rotated latitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jrlatVarID, "units", "degrees")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  IF (.NOT. outblock%luvmasspoint .AND. yextension /= 'c') THEN  
    zsrlon = zrotlon + dlon*0.5_wp
    zsrlat = zrotlat + dlat*0.5_wp

    ierror=nf90_def_dim(ncid,"srlon",outblock%ie_out_tot, idims_id_out(9))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

    ierror=nf90_def_dim(ncid,"srlat",outblock%je_out_tot, idims_id_out(10))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

    ! shifted rotated longitude
    ierror=nf90_def_var(ncid, "srlon", NF90_FLOAT, (/ idims_id_out(9) /),     &
                                                 jsrlonVarID)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

    ierror=nf90_put_att(ncid, jsrlonVarID, "standard_name", "grid_longitude")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jsrlonVarID, "long_name",                       &
                                           "staggered rotated longitude")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jsrlonVarID, "units", "degrees")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

    ! shifted rotated latitude
    ierror=nf90_def_var(ncid, "srlat", NF90_FLOAT, (/ idims_id_out(10) /),    &
                                                 jsrlatVarID)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jsrlatVarID, "standard_name", "grid_latitude")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jsrlatVarID, "long_name",                       &
                                           "staggered rotated latitude")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jsrlatVarID, "units", "degrees")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF


    DO j = 1, outblock%je_out_tot
      DO i = 1, outblock%ie_out_tot
        zslonu(i,j) = rlarot2rla (zrotlat(j), zsrlon(i), pollat, pollon, polgam)
        zslatu(i,j) = phirot2phi (zrotlat(j), zsrlon(i), pollat, pollon, polgam)
        zslonv(i,j) = rlarot2rla (zsrlat(j), zrotlon(i), pollat, pollon, polgam)
        zslatv(i,j) = phirot2phi (zsrlat(j), zrotlon(i), pollat, pollon, polgam)
      ENDDO
    ENDDO

  ENDIF  ! luvmasspoint

  DO j = 1, outblock%je_out_tot
    DO i = 1, outblock%ie_out_tot
      zlongitude(i,j) = rlarot2rla (zrotlat(j), zrotlon(i), pollat, pollon, polgam)
      zlatitude(i,j)  = phirot2phi (zrotlat(j), zrotlon(i), pollat, pollon, polgam)
    ENDDO
  ENDDO

! set the values for rotated and true geographic longitudes and latitudes

! geographic longitude
  ierror=nf90_def_var(ncid, "lon", NF90_FLOAT, (/ idims_id_out(1),        &
                                               idims_id_out(2) /), jlonVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_att(ncid, jlonVarID, "standard_name", "longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jlonVarID, "long_name", "longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jlonVarID, "units", "degrees_east")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ! geographic latitude
  ierror=nf90_def_var(ncid, "lat", NF90_FLOAT, (/ idims_id_out(1),         &
                                               idims_id_out(2) /), jlatVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jlatVarID, "standard_name", "latitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jlatVarID, "long_name", "latitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jlatVarID, "units", "degrees_north")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF


  IF (.NOT. outblock%luvmasspoint .AND. yextension /= 'c') THEN  

    ! shifted geographic longitude or U wind component
    ierror=nf90_def_var(ncid, "slonu", NF90_FLOAT,                         &
                        (/ idims_id_out(9), idims_id_out(2) /), jslonuVarID)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

    ierror=nf90_put_att(ncid, jslonuVarID, "standard_name", "longitude")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jslonuVarID, "long_name",                    &
                                           "staggered U-wind longitude")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jslonuVarID, "units", "degrees_east")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

    ! shifted geographic latitude in U wind component
    ierror=nf90_def_var(ncid, "slatu", NF90_FLOAT,                         &
                       (/ idims_id_out(9), idims_id_out(2) /), jslatuVarID)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jslatuVarID, "standard_name", "latitude")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jslatuVarID, "long_name",                    &
                                           "staggered U-wind latitude")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jslatuVarID, "units", "degrees_north")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

    ! shifted geographic longitude or V wind component
    ierror=nf90_def_var(ncid, "slonv", NF90_FLOAT,                           &
                       (/ idims_id_out(1), idims_id_out(10) /), jslonvVarID)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

    ierror=nf90_put_att(ncid, jslonvVarID, "standard_name", "longitude")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jslonvVarID, "long_name",                      &
                                           "staggered V-wind longitude")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jslonvVarID, "units", "degrees_east")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

    ! shifted geographic latitude in V wind component
    ierror=nf90_def_var(ncid, "slatv", NF90_FLOAT,                           &
                       (/ idims_id_out(1), idims_id_out(10) /), jslatvVarID)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jslatvVarID, "standard_name", "latitude")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jslatvVarID, "long_name",                      &
                                           "staggered V-wind latitude")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jslatvVarID, "units", "degrees_north")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

  ENDIF ! luvmasspoint

! take care of different height co-ordinates and define the vertical 
! axis accordingly
  IF (yextension == 'p') THEN   ! pressure co-ordinate

    zvcoord_len = outblock%kepin
    ALLOCATE (zvcoord(outblock%kepin))
    zvcoord(1:outblock%kepin) = outblock%plev(1:outblock%kepin)
    
    ierror=nf90_def_dim(ncid,"pressure",outblock%kepin,idims_id_out(3))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_def_var(ncid, "pressure", NF90_FLOAT,                    &
                       (/ idims_id_out(3) /), jvcVarID)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jvcVarID, "long_name", "pressure")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jvcVarID, "units", "Pa")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jvcVarID, "positive", "down")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

  ELSEIF (yextension == 'z') THEN  ! geometric height co-ordinate
 
    zvcoord_len = outblock%kezin 
    ALLOCATE (zvcoord(outblock%kezin))
    zvcoord(1:outblock%kezin) = outblock%zlev(1:outblock%kezin)

    ierror=nf90_def_dim(ncid,"altitude",outblock%kezin,idims_id_out(3))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_def_var(ncid, "altitude", NF90_FLOAT,                    &
                       (/ idims_id_out(3) /), jvcVarID)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jvcVarID, "long_name",                     &
                                        "height above mean sea level")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jvcVarID, "units", "m")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jvcVarID, "positive", "up")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

!_br 25.11.08
  ELSE ! model coordinate

      zvcoord_len = ke1
      ALLOCATE (zvcoord(ke1))
      IF     (vcoord%ivctype == 1) THEN
        zvcoord(1:ke1) = vcoord%sigm_coord(1:ke1)
      ELSEIF ( ANY( vcoord%ivctype == (/2,3,4/) ) ) THEN
        zvcoord(1:ke1) = vcoord%vert_coord(1:ke1)
      ENDIF
      
      ierror=nf90_def_dim(ncid,"level",ke, idims_id_out(3))
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror=nf90_def_dim(ncid,"level1",ke1, idims_id_out(4))
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror=nf90_def_var(ncid, "vcoord", NF90_FLOAT,                     &
                         (/ idims_id_out(4) /), jvcVarID)
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      IF (vcoord%ivctype == 1) THEN
        ierror=nf90_put_att(ncid, jvcVarID, "long_name",                     &
                                           "Pressure based hybrid coordinate")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ELSE IF (vcoord%ivctype == 2) THEN
        ierror=nf90_put_att(ncid, jvcVarID, "long_name",                     &
                                    "Height-based hybrid Gal-Chen coordinate")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ELSE IF (vcoord%ivctype == 3) THEN
        ierror=nf90_put_att(ncid, jvcVarID, "long_name",                     &
                                                           "SLEVE coordinate")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ELSE IF (vcoord%ivctype == 4) THEN
        ierror=nf90_put_att(ncid, jvcVarID, "long_name",                     &
                                                           "SLEVE2 coordinate")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ENDIF

      ierror=nf90_put_att(ncid, jvcVarID, "units", "Pa")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror=nf90_put_att(ncid, jvcVarID, "ivctype", vcoord%ivctype)
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror=nf90_put_att(ncid, jvcVarID, "irefatm", refatm%irefatm)
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror=nf90_put_att(ncid, jvcVarID, "p0sl", refatm%p0sl)
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror=nf90_put_att(ncid, jvcVarID, "t0sl", refatm%t0sl)
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror=nf90_put_att(ncid, jvcVarID, "dt0lp", refatm%dt0lp)
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror=nf90_put_att(ncid, jvcVarID, "vcflat", vcoord%vcflat)
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF

      IF ( ANY( vcoord%ivctype == (/3,4/) ) ) THEN
        ierror=nf90_put_att(ncid, jvcVarID, "svc1", svc1)
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jvcVarID, "svc2", svc2)
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jvcVarID, "nfltvc", nfltvc)
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ENDIF
      IF (refatm%irefatm == 2) THEN
        ierror=nf90_put_att(ncid, jvcVarID, "delta_t", refatm%delta_t)
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jvcVarID, "h_scal", refatm%h_scal)
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ELSE IF (refatm%irefatm == 3) THEN
        ierror=nf90_put_att(ncid, jvcVarID, "bvref", refatm%bvref)
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ENDIF

      IF (lradtopo) THEN
         ierror=nf90_def_dim(ncid,"nhori", nhori, idims_id_out(11))
         IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
         ENDIF
         ierror=nf90_def_var(ncid, "nhori", NF90_FLOAT,                 &
              (/ idims_id_out(11) /), jsectVarID)
         IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
         ENDIF
         ierror=nf90_put_att(ncid, jsectVarID, "standard_name", "sectors")
         IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
         ENDIF
         ierror=nf90_put_att(ncid, jsectVarID, "long_name", "sectors of horizontal angles around grid cells")
         IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
         ENDIF
         ierror=nf90_put_att(ncid, jsectVarID, "units", "sector number starting from the North clockwise")
         IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
         ENDIF
      ENDIF

      IF (yextension /= 'c') THEN
        ! height_2m, height_10m, height_toa, and wbtemp_13c
        ! not needed in files holding the constant fields

        zheight_2m  = 2._wp
        zheight_10m = 10._wp
        zheight_toa = hhl(1,1,1)
        zwbtemp_13c = 1.3_wp

        ierror=nf90_def_var(ncid, "height_2m", NF90_FLOAT, jh2VarID)
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jh2VarID, "standard_name", "height")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jh2VarID, "long_name",                   &
                                            "height above the surface")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jh2VarID, "units", "m")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jh2VarID, "positive", "up")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_def_var(ncid, "height_10m", NF90_FLOAT, jh10VarID)
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jh10VarID, "standard_name", "height")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jh10VarID, "long_name",                  &
                                             "height above the surface")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jh10VarID, "units", "m")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jh10VarID, "positive", "up")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_def_var(ncid, "height_toa", NF90_FLOAT, jhtoaVarID)
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jhtoaVarID, "standard_name", "height")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jhtoaVarID, "long_name", "height of top of model")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jhtoaVarID, "units", "m")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jhtoaVarID, "positive", "up")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_def_var(ncid, "wbt_13c", NF90_FLOAT, jwbt13VarID)
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jwbt13VarID, "long_name", "wet bulb temperature")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jwbt13VarID, "units", "Celsius")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF

      ENDIF
!_br 25.11.08 end

  ENDIF

! define fields of variable multi layer soil model
  IF (lmulti_layer .AND. &
     yextension /= 'c' .AND. yextension /= 'p' .AND. yextension /= 'z') THEN
  

    ierror=nf90_def_dim(ncid, "soil", ke_soil, idims_id_out(7))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_def_dim(ncid, "soil1", ke_soil+1, idims_id_out(8))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF


    ierror=nf90_def_var(ncid, "soil1", NF90_FLOAT, (/ idims_id_out(8) /),  &
                                                   jsoilVarID)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
 
     ierror=nf90_put_att(ncid, jsoilVarID, "standard_name", "depth")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jsoilVarID, "long_name", "depth of soil layers")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jsoilVarID, "units", "m")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jsoilVarID, "positive", "down")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jsoilVarID, "bounds", "soil1_bnds")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

    ierror=nf90_def_var(ncid, "soil1_bnds", NF90_FLOAT,                    &
                       (/ idims_id_out(6), idims_id_out(8) /), jsoilbdsID)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jsoilbdsID, "long_name",                     &
                                          "boundaries of soil layers")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

    ierror=nf90_put_att(ncid, jsoilbdsID, "units", "m")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

    ! define multiple snow layers
    IF (lmulti_snow) THEN
       ierror=nf90_def_dim(ncid, "snow_layer", ke_snow, idims_id_out(14))
       IF (ierror /= NF90_NOERR) THEN
         yerrmsg = TRIM(NF90_strerror(ierror))
         RETURN
       ENDIF
      ierror=nf90_def_var(ncid, "snow_layer", NF90_FLOAT, (/ idims_id_out(14) /),  &
                                                     jsnowVarID)
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror=nf90_put_att(ncid, jsnowVarID, "long_name", "index of snow layers")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror=nf90_put_att(ncid, jsnowVarID, "units", "1")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror=nf90_put_att(ncid, jsnowVarID, "positive", "down")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ENDIF

  ENDIF  ! end defining fields of variable multi layer soil model

#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
  IF (luse_rttov .AND. yextension == 's') THEN
    ierror=nf90_def_dim(ncid, "nsynmsg", nsynmsg, idims_id_out(12))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

    ierror=nf90_def_dim(ncid, "msgchan", nmsgchan, idims_id_out(13))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

    zmsgchan_wave = (/ 1.6_wp, 6.2_wp, 7.3_wp, 3.8_wp, 8.7_wp, 10.8_wp, 12.0_wp, 9.7_wp /)
    ierror=nf90_def_var(ncid, "msgchan_wave", NF90_FLOAT, (/ idims_id_out(13) /),  &
                                                   jmsgchanVarID)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

    ierror=nf90_put_att(ncid, jmsgchanVarID, "standard_name", "channel_wave_length")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jmsgchanVarID, "long_name", "channel wave length")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jmsgchanVarID, "units", "micrometer")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
  ENDIF
#endif

! End of definition mode
!!! no more attribute definitions beyond this line !!!

  ! Enddef routine: reserve space on the header (100K) of the netcdf, 
  ! then if we need to add variables (after some data was dumped to the file) 
  ! or more info to the header, 
  ! it won't involve an expensive relocation of all the data present
  ! Nevertheless, for some reasons not clear to me, this parameters seems
  ! to be ignored when opening the file in parallel among all PE in a communicator
  ierror=nf90_enddef(ncid, h_minfree=102400, v_minfree=102400 )
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

ENDIF

IF (my_comm_id == 0) THEN

! write the values of longitude, latitude and vertical axis to file
  ierror=nf90_put_var(ncid, jlonVarID, zlongitude, start=(/ 1,1 /),&
     count=(/ outblock%ie_out_tot,outblock%je_out_tot /))

  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_var(ncid, jlatVarID, zlatitude, start=(/ 1,1 /),&
             count=(/ outblock%ie_out_tot,outblock%je_out_tot /))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_var(ncid, jrlonVarID, zrotlon, start=(/ 1 /),   &
             count=(/ outblock%ie_out_tot  /))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_var(ncid, jrlatVarID, zrotlat, start=(/ 1 /),   &
             count=(/ outblock%je_out_tot /))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  IF (.NOT. outblock%luvmasspoint .AND. yextension /= 'c') THEN

    ierror=nf90_put_var(ncid, jsrlonVarID, zsrlon, start=(/ 1 /),   &
             count=(/ outblock%ie_out_tot /))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_var(ncid, jsrlatVarID, zsrlat, start=(/ 1 /),   &
             count=(/ outblock%je_out_tot /))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_var(ncid, jslonuVarID, zslonu, start=(/ 1,1 /),   &
             count=(/ outblock%ie_out_tot,outblock%je_out_tot /))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_var(ncid, jslatuVarID, zslatu, start=(/ 1,1 /),   &
             count=(/ outblock%ie_out_tot,outblock%je_out_tot  /))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_var(ncid, jslonvVarID, zslonv, start=(/ 1,1 /),   &
             count=(/ outblock%ie_out_tot,outblock%je_out_tot /))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_var(ncid, jslatvVarID, zslatv, start=(/ 1,1 /),   &
             count=(/ outblock%ie_out_tot,outblock%je_out_tot /))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

  ENDIF


  ierror=nf90_put_var(ncid, jvcVarID, zvcoord, start=(/ 1 /),   &
             count=(/ zvcoord_len /))   
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  IF ( yextension /= 'c' .AND. yextension /= 'p' .AND. yextension /= 'z') THEN 
    ierror=nf90_put_var(ncid, jh2VarID, zheight_2m, start=(/ 1 /),   &
             count=(/ 1 /))   
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_var(ncid, jh10VarID, zheight_10m, start=(/ 1 /),   &
             count=(/ 1 /))   
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_var(ncid, jhtoaVarID, zheight_toa,  start=(/ 1 /),   &
             count=(/ 1 /))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_var(ncid, jwbt13VarID, zwbtemp_13c, start=(/ 1 /),   &
             count=(/ 1 /))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

    IF (lmulti_layer .AND. &
     yextension /= 'c' .AND. yextension /= 'p' .AND. yextension /= 'z') THEN
      ierror=nf90_put_var(ncid, jsoilVarID, czmls(1:ke_soil+1), start=(/ 1 /),   &
             count=(/ ke_soil+1 /))
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      zsoil_bnds(1, 1) = 0._irealgrib
      zsoil_bnds(1, 2:ke_soil+1) = REAL(czhls(1:ke_soil),irealgrib)
      zsoil_bnds(2, 1:ke_soil+1) = REAL(czhls(1:ke_soil+1),irealgrib)
      ierror=nf90_put_var(ncid, jsoilbdsID, zsoil_bnds, start=(/ 1,1 /),   &
             count=(/ nbnds, ke_soil+1 /))
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF

      IF (lmulti_snow) THEN
        DO k = 1, ke_snow
          zind_snow(k) = REAL(k,irealgrib)
        ENDDO
        ierror=nf90_put_var(ncid, jsnowVarID, zind_snow(1:ke_snow), start=(/ 1 /),   &
             count=(/ ke_snow /))
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ENDIF

    ENDIF
  ENDIF

  IF (lradtopo .AND. yextension /= 'p' .AND. yextension /= 'z') THEN
     DO isect=1,nhori
        zsect(isect)=isect
     ENDDO
     ierror=nf90_put_var(ncid, jsectVarID, zsect, start=(/ 1 /),   &
             count=(/ nhori /))
     IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
     ENDIF
  ENDIF

#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
  ! for the synthetic satellite images
  IF (luse_rttov .AND. yextension == 's') THEN
    ierror=nf90_put_var(ncid, jmsgchanVarID, zmsgchan_wave, start=(/ 1 /),   &
             count=(/ nmsgchan /))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
  ENDIF
#endif

! write the values of time and time_bnds to file
  ierror=nf90_put_var(ncid, jtimeID, time, start = (/ 1 /),  &
         count = (/ ntime  /) )   
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_var(ncid, jtbdsID, time_bnds, start = (/ 1,1 /),  &
         count = (/ nbnds, ntime  /))  
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

ENDIF ! my_comm_id == 0

IF(lasync_io .OR. my_comm_id == 0) THEN
  DEALLOCATE (zvcoord)
ENDIF

IF (npes > 1) THEN
  CALL MPI_BARRIER(icomm,implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'MPI_BARRIER failed'
    ierror  = 1
    RETURN
  ENDIF
ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine 
!------------------------------------------------------------------------------

END SUBROUTINE write_nc_gdefs

!==============================================================================
!+ Module procedure in "io_utilities" for writing a NetCDF file
!------------------------------------------------------------------------------

SUBROUTINE write_nc_vdefs (ncid, numlist, ilist, var_id, luvmasspoint,      &
                           icomm, npes, yextension, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   write NetCDF attributes for each output parameter
!
!------------------------------------------------------------------------------

! Scalar arguments with intent(in):
  INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
    ncid             ! NetCDF file ID

  INTEGER (KIND=iintegers), INTENT(IN)     ::  &
    numlist          ! number of variables for output

  INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
    icomm,    & ! MPI communicator
    npes        ! number of PEs

   LOGICAL                   , INTENT(IN) :: &
   luvmasspoint     ! interpolate horizontal winds to mass grid points

  CHARACTER (LEN=1),          INTENT(IN)   ::  &
    yextension    ! indicates model variables (''), p-('p') or z-levels ('z')

! Array arguments with intent(in):
  INTEGER (KIND=iintegers), INTENT(IN)     ::  &
    ilist(3,numlist) ! 

! Array arguments with intent(out):
  INTEGER (KIND=iintegers), INTENT(OUT)    ::  &
    var_id(nzmxid)      ! NetCDF-ID of each variable 

  CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
    yerrmsg     ! string for error messages

! Scalar arguments with intent(out):
  INTEGER (KIND=iintegers), INTENT(OUT)     ::  &
    ierror          ! error index

!------------------------------------------------------------------------------

! Local scalars:
  INTEGER (KIND=iintegers)  :: i1,i2,i3, iid, jid, n, timestep
  CHARACTER (LEN= 60)       :: ymethod

  INTEGER (KIND=iintegers)  :: my_comm_id, implcode

  REAL (KIND=wp)            :: zdtime

! string variable to hold grid information
CHARACTER (LEN=200)     grid_mapping

! k dim of 3d and 4d variables 
INTEGER (KIND=iintegers)  :: vardim3_p3,vardim3_p4

!
! Local arrays:
!
!- End of header
!==============================================================================

vardim3_p3=-1
vardim3_p4=-1

grid_mapping = 'rotated_pole'

IF (npes > 1) THEN
 ! Get id in communicator comm
  CALL MPI_COMM_RANK(icomm,my_comm_id,implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'MPI_COMM_RANK failed'
    ierror  = 2
    RETURN
  ENDIF
ELSE
  my_comm_id = 0
ENDIF

var_id(:) = 0

IF( lasync_io .OR. my_comm_id == 0) THEN
  ierror=0

! Re-Enter definition mode
  ierror = NF90_redef (ncid)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ! loop over all variables that should be written and loop over all variables
  ! in the LM variable table until equal elements are found
  DO n = 1, numlist
    ! indices of field in variable table
    i1 = ilist(1,n)
    i2 = ilist(2,n)
    i3 = ilist(3,n)

    ! select netCDF dimension IDs
    IF( ( (TRIM(var(i1,i2,i3)%name) == 'U')            .OR.      &
          (TRIM(var(i1,i2,i3)%name) == 'AUMFL_S') )    .AND.     &
        (.NOT. luvmasspoint) ) THEN
       iid = 9
       jid = 2
    ELSEIF( ( (TRIM(var(i1,i2,i3)%name) == 'V')        .OR.      &
          (TRIM(var(i1,i2,i3)%name) == 'AVMFL_S') )    .AND.     &
        (.NOT. luvmasspoint) ) THEN
       iid = 1
       jid = 10
    ELSE
       iid = 1
       jid = 2
    ENDIF

! set the dimensions of the variable regarding to its rank and other criterea

    IF(lasync_io ) THEN
      vardim3_p3 = var(i1,i2,i3)%idimvert
      vardim3_p4 = var(i1,i2,i3)%idimvert
    ELSE
      IF( .NOT. ASSOCIATED(var(i1,i2,i3)%p4) ) THEN
        vardim3_p4 = -1
      ELSE
        vardim3_p4 = UBOUND(var(i1,i2,i3)%p4,3)
      ENDIF
      IF( .NOT. ASSOCIATED(var(i1,i2,i3)%p3) ) THEN
        vardim3_p3 = -1
      ELSE
        vardim3_p3 = UBOUND(var(i1,i2,i3)%p3,3)
      ENDIF

    ENDIF

    SELECT CASE (var(i1,i2,i3)%rank)

    CASE(4)

      IF ( (yextension == 'p') .OR. (yextension == 'z') .OR.               &
                            (vardim3_p4 == ke) ) THEN
        ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT,  &
                          (/ idims_id_out(iid), idims_id_out(jid),         &
                             idims_id_out(3), idims_id_out(5) /), var_id(n))
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ELSEIF ( vardim3_p4 == ke1) THEN
        ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT,  &
                          (/ idims_id_out(iid), idims_id_out(jid),         &
                             idims_id_out(4), idims_id_out(5) /), var_id(n))
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ELSEIF ( vardim3_p4 == ke_soil+1) THEN
        ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT,  &
                          (/ idims_id_out(iid), idims_id_out(jid),         &
                             idims_id_out(8), idims_id_out(5) /), var_id(n))
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF

      ELSEIF (UBOUND(var(i1,i2,i3)%p4,3) == ke_snow) THEN
        ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT,  &
                          (/ idims_id_out(iid), idims_id_out(jid),         &
                             idims_id_out(14), idims_id_out(5) /), var_id(n))
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF

      ENDIF

    CASE(3)

      IF (yextension == 'p' .OR. yextension == 'z') THEN
        ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT,  &
                          (/ idims_id_out(iid), idims_id_out(jid),         &
                             idims_id_out(3), idims_id_out(5) /), var_id(n))
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ELSE
        IF ( vardim3_p3 == -1 ) THEN
          IF (var(i1,i2,i3)%levtyp == 109) THEN
            ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT,&
                            (/ idims_id_out(iid), idims_id_out(jid),         &
                               idims_id_out(4), idims_id_out(5) /), var_id(n))
            IF (ierror /= NF90_NOERR) THEN
              yerrmsg = TRIM(NF90_strerror(ierror))
              RETURN
            ENDIF
          ELSEIF (var(i1,i2,i3)%levtyp == 110) THEN
            ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT,&
                            (/ idims_id_out(iid), idims_id_out(jid),         &
                               idims_id_out(3), idims_id_out(5) /), var_id(n))
            IF (ierror /= NF90_NOERR) THEN
              yerrmsg = TRIM(NF90_strerror(ierror))
              RETURN
            ENDIF
          ELSE
            ierror = 1
            WRITE(yerrmsg,'(A,I4,A)') 'levtyp = ',var(i1,i2,i3)%levtyp,     &
                                      ' not implemented'
          ENDIF
        ELSEIF (  vardim3_p3 <= 3) THEN
          ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT, &
                            (/ idims_id_out(iid), idims_id_out(jid),        &
                               idims_id_out(5) /), var_id(n))
          IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
          ENDIF
        ELSEIF (  vardim3_p3 == ke) THEN
          ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT, &
                            (/ idims_id_out(iid), idims_id_out(jid),        &
                               idims_id_out(3), idims_id_out(5) /), var_id(n))
          IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
          ENDIF
        ELSEIF (  vardim3_p3 == ke1) THEN
          ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT, &
                            (/ idims_id_out(iid), idims_id_out(jid),        &
                               idims_id_out(4), idims_id_out(5) /), var_id(n))
          IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
          ENDIF
        ELSEIF ( vardim3_p3 == nhori) THEN
          ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT, &
               (/ idims_id_out(iid), idims_id_out(jid),        &
               idims_id_out(11), idims_id_out(5) /), var_id(n))
          IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
          ENDIF
        ELSE
          ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT, &
                            (/ idims_id_out(iid), idims_id_out(jid),        &
                               idims_id_out(3), idims_id_out(5) /), var_id(n))
          IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
          ENDIF
        ENDIF
      ENDIF
    
      ! for the synthetic satellite images
      IF (var(i1,i2,i3)%levtyp == 222) THEN
        IF (TRIM(var(i1,i2,i3)%name) == 'SYNMSG') THEN
            ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT,&
                            (/ idims_id_out(iid), idims_id_out(jid),         &
                               idims_id_out(12), idims_id_out(5) /), var_id(n))
            IF (ierror /= NF90_NOERR) THEN
              yerrmsg = TRIM(NF90_strerror(ierror))
              RETURN
            ENDIF
        ELSE
            ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT,&
                            (/ idims_id_out(iid), idims_id_out(jid),         &
                               idims_id_out(13), idims_id_out(5) /), var_id(n))
            IF (ierror /= NF90_NOERR) THEN
              yerrmsg = TRIM(NF90_strerror(ierror))
              RETURN
            ENDIF
        ENDIF
      ENDIF

    CASE(2)

      IF (var(i1,i2,i3)%levtyp == 105) THEN
      
        ! variables on a specific altitude
        IF (var(i1,i2,i3)%levbot == 2) THEN
          ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT, &
                            (/ idims_id_out(iid), idims_id_out(jid),        &
                               idims_id_out(5) /), var_id(n))
          IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
          ENDIF
        ELSE IF (var(i1,i2,i3)%levbot == 10) THEN
          ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT, &
                            (/ idims_id_out(iid), idims_id_out(jid),        &
                               idims_id_out(5) /), var_id(n))
          IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
          ENDIF
        ELSE
          ierror = 1
          yerrmsg = 'invalid "levbot" '
          RETURN
        ENDIF
       
      ELSE IF (var(i1,i2,i3)%levtyp == 8) THEN

        ! variables at TOA (actually top of model)
        ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT, &
                       (/ idims_id_out(iid), idims_id_out(jid),           &
                          idims_id_out(5) /), var_id(n))
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF

      ELSE IF (TRIM(var(i1,i2,i3)%name) == 'SNOWLMT') THEN

        ! Snow limit has coordinate variable wet bulb temp.
        ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT, &
                       (/ idims_id_out(iid), idims_id_out(jid),           &
                          idims_id_out(5) /), var_id(n))
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF

      ELSE
      
        ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT,  &
                          (/ idims_id_out(iid), idims_id_out(jid),         &
                             idims_id_out(5) /), var_id(n))
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        
      ENDIF
    
    END SELECT
   
! set the attributes of the variable
    IF (var(i1,i2,i3)%standard_name(1:1) /= '-') THEN
      ierror = nf90_put_att (ncid, var_id(n), "standard_name",               &
                                              TRIM(var(i1,i2,i3)%standard_name))
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ENDIF
    ierror = nf90_put_att (ncid, var_id(n), "long_name",                   &
                                            TRIM(var(i1,i2,i3)%long_name))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror = nf90_put_att (ncid, var_id(n), "units", TRIM(var(i1,i2,i3)%units))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror = nf90_put_att (ncid, var_id(n), "grid_mapping", grid_mapping)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    IF(  (TRIM(var(i1,i2,i3)%name) == 'HTOP_CON')     .OR.      &
         (TRIM(var(i1,i2,i3)%name) == 'HBAS_CON')     .OR.      &
         (TRIM(var(i1,i2,i3)%name) == 'HZEROCL')      .OR.      &
         (TRIM(var(i1,i2,i3)%name) == 'SNOWLMT')      .OR.      &
         (var(i1,i2,i3)%lsm /= ' ') ) THEN
      ierror = nf90_put_att (ncid, var_id(n), "_FillValue", undefncdf)
      IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
      ENDIF
    ENDIF

    IF  (TRIM(var(i1,i2,i3)%name) == 'HMO3') THEN
      ierror = nf90_put_att (ncid, var_id(n), "cell_methods", &
                          "mole_fraction_of_ozone_in_air: maximum")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ENDIF

    IF (TRIM(var(i1,i2,i3)%name) == 'SOILTYP') THEN
      ierror = nf90_put_att (ncid, var_id(n), "flag_values", &
                            (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 0 /))
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror = nf90_put_att (ncid, var_id(n), "flag_meanings", &
         "ice rock sand sandy_loam loam clay_loam clay peat sea_water sea_ice")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ENDIF

! set attribute for vertical coordinate
    IF (TRIM(var(i1,i2,i3)%name) == 'HHL') THEN
      ierror=nf90_put_att(ncid, var_id(n), "positive", "up")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ENDIF

    IF( ((TRIM(var(i1,i2,i3)%name) == 'U')            .OR.      &
        (TRIM(var(i1,i2,i3)%name) == 'AUMFL_S'))      .AND.     &
        (.NOT. luvmasspoint) ) THEN
      ierror = nf90_put_att (ncid, var_id(n), "coordinates", "slonu slatu")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ELSEIF( ((TRIM(var(i1,i2,i3)%name) == 'V')        .OR.      &
        (TRIM(var(i1,i2,i3)%name) == 'AVMFL_S'))      .AND.     &
        (.NOT. luvmasspoint) ) THEN
      ierror = nf90_put_att (ncid, var_id(n), "coordinates", "slonv slatv")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ELSEIF (var(i1,i2,i3)%levtyp == 105) THEN
      IF (var(i1,i2,i3)%levbot == 2) THEN
        ierror = nf90_put_att (ncid, var_id(n), "coordinates", "lon lat height_2m")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ENDIF
      IF (var(i1,i2,i3)%levbot == 10) THEN
        ierror = nf90_put_att (ncid, var_id(n), "coordinates", "lon lat height_10m")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ENDIF
    ELSE IF (var(i1,i2,i3)%levtyp == 8) THEN
      ierror = nf90_put_att (ncid, var_id(n), "coordinates", "lon lat height_toa")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ELSE IF (TRIM(var(i1,i2,i3)%name) == 'SNOWLMT') THEN
      ierror = nf90_put_att (ncid, var_id(n), "coordinates", "lon lat wbtemp_13c")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ELSE
      ierror = nf90_put_att (ncid, var_id(n), "coordinates", "lon lat")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ENDIF

! Time Range Indicator determines the cell_methods attribute
! time will be in seconds
    SELECT CASE(var(i1,i2,i3)%ntri)
    CASE(2)
      IF(var(i1,i2,i3)%name == 'TMIN_2M') THEN 
        ierror = nf90_put_att (ncid, var_id(n), "cell_methods", "time: minimum")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ELSE IF(var(i1,i2,i3)%name == 'TMAX_2M') THEN
        ierror = nf90_put_att (ncid, var_id(n), "cell_methods", "time: maximum")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ELSE IF(var(i1,i2,i3)%name == 'VMAX_10M') THEN
        ierror = nf90_put_att (ncid, var_id(n), "cell_methods", "time: maximum")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ELSE IF(var(i1,i2,i3)%name == 'VABSMX_10M') THEN
        ierror = nf90_put_att (ncid, var_id(n), "cell_methods", "time: maximum")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ELSE IF(var(i1,i2,i3)%name == 'VGUST_CON') THEN
        ierror = nf90_put_att (ncid, var_id(n), "cell_methods", "time: maximum")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ELSE IF(var(i1,i2,i3)%name == 'VGUST_DYN') THEN
        ierror = nf90_put_att (ncid, var_id(n), "cell_methods", "time: maximum")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ENDIF
    CASE(3)
      ierror = nf90_put_att (ncid, var_id(n), "cell_methods", "time: mean")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    CASE(4)
      ierror = nf90_put_att (ncid, var_id(n), "cell_methods", "time: sum")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    END SELECT 

! Set attribute for layer means of the soil layers
    IF (  vardim3_p4  /= -1) THEN
      IF ( (var(i1,i2,i3)%rank == 4) .AND.                                 &
           ( vardim3_p4  == ke_soil+1) ) THEN
        ierror = nf90_put_att (ncid, var_id(n), "cell_methods", "soil1: mean")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ENDIF
    ENDIF

  ENDDO     ! End of loop over all variables


! End of definition mode
!!! no more attribute definitions beyond this line !!!

  ierror=nf90_enddef(ncid)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

ENDIF

IF (npes > 1) THEN
  CALL distribute_values  (var_id, nzmxid, 0, imp_integers,  icomm, implcode)
  CALL MPI_BARRIER(icomm,implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'MPI_BARRIER failed'
    ierror  = 1
    RETURN
  ENDIF
ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine 
!------------------------------------------------------------------------------

END SUBROUTINE write_nc_vdefs


SUBROUTINE distribute_values_asynio(izerror, yzerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   This routine distributes initialized values, namelist properties, etc
!   from a compute PE (i.e. PE 0) to all the I/O PE, so that all the I/O PE
!   have the same information needed for their I/O tasks.
!
!------------------------------------------------------------------------------

INTEGER (KIND=iintegers), INTENT(INOUT)     :: izerror
CHARACTER (LEN=80), INTENT(INOUT)           :: yzerrmsg

! Local variables

! pointer fo gribout information from namelist
TYPE (pp_nl), POINTER                 :: now


INTEGER (KIND=iintegers)  :: ibuflen,rbuflen, istat,i,j,ierr,  &
       nziope0, realcnt,intcnt


! integer buffer to pack variables 
INTEGER (KIND=iintegers), ALLOCATABLE ::   intbuf  (:)

! real buffer to pack variables
REAL    (KIND=wp)       , ALLOCATABLE :: realbuf(:)

INTEGER (KIND=iintegers)              :: ar_des_size
INTEGER (KIND=iintegers)              :: vartab_bufsize, var_size
CHARACTER (LEN=1), ALLOCATABLE         :: vartab_buf(:)

  nziope0 = nproc - nc_asyn_io
  izerror = 0

  IF( nc_asyn_io > 0 .AND. &
      (my_world_id == 0 .OR. my_world_id >= nziope0 ) ) THEN

    ibuflen = 0

    IF( my_world_id == 0 ) THEN

      now=>root 

      gribout_cnt: DO

        ibuflen = ibuflen + 5

        IF (ASSOCIATED(now%next) ) THEN
          now => now%next
        ELSE
          EXIT gribout_cnt
        ENDIF

      ENDDO gribout_cnt

      ALLOCATE( intbuf(ibuflen),  STAT=istat )

      intcnt=0
      now=>root
      gribout_set: DO

        intbuf(intcnt+1) = now%nyvar_m
        intbuf(intcnt+2) = now%nyvar_p
        intbuf(intcnt+3) = now%nyvar_z
        intbuf(intcnt+4) = now%nyvar_s
        intbuf(intcnt+5) = now%nyvar_c
        intcnt = intcnt + 5

        IF (ASSOCIATED(now%next) ) THEN
          now => now%next
        ELSE
          EXIT gribout_set
        ENDIF

      ENDDO gribout_set

    ENDIF

    
    CALL MPI_BCAST(intcnt,1,imp_integers,0,intercomm_asynio,ierr)

    IF( my_world_id >= nziope0 ) THEN
      ALLOCATE ( intbuf(intcnt)   , STAT=istat )
      IF( istat /= 0) THEN
        izerror = 33340
        yzerrmsg = 'ERROR allocating buffer of integers for IO distribution'
        RETURN
      ENDIF
    ENDIF

    CALL distribute_values  (intbuf , intcnt, 0, imp_integers, intercomm_asynio, ierr)

    IF( my_world_id >= nziope0 ) THEN

      intcnt=0
      now=>root
      gribout_set2: DO

        now%nyvar_m = intbuf(intcnt+1)
        now%nyvar_p = intbuf(intcnt+2)
        now%nyvar_z = intbuf(intcnt+3)
        now%nyvar_s = intbuf(intcnt+4)
        now%nyvar_c = intbuf(intcnt+5)
        intcnt = intcnt + 5

        IF (ASSOCIATED(now%next) ) THEN
          now => now%next
        ELSE
          EXIT gribout_set2
        ENDIF

      ENDDO gribout_set2

    ENDIF
   
    DEALLOCATE(intbuf) 

    ibuflen  = 2    ! should be long enough for very long Namelists

    rbuflen = ke1+7+ke_soil+1

    ALLOCATE ( realbuf(rbuflen),  STAT=istat )
    IF( istat /= 0) THEN
        izerror = 33340
        yzerrmsg = 'ERROR allocating buffer of reals for IO distribution'
        RETURN
    ENDIF
 
    IF(my_world_id == 0) THEN

      now => root
      gribout_lpt: DO

        ibuflen = ibuflen + 3*now%nyvar_m + 3*now%nyvar_p + 3*now%nyvar_z + &
          3*now%nyvar_s + 3*now%nyvar_c + 1

        IF (ASSOCIATED(now%next) ) THEN
          now => now%next
        ELSE
          EXIT gribout_lpt
        ENDIF

      ENDDO gribout_lpt

      ALLOCATE ( intbuf(ibuflen)   , STAT=istat )
      IF( istat /= 0) THEN
        izerror = 33340
        yzerrmsg = 'ERROR allocating buffer of integers for IO distribution'
        RETURN
      ENDIF
 
      intbuf(1) = vcoord%ivctype
      intbuf(2) = refatm%irefatm

      now => root
      intcnt=2
      gribout_loop_snd: DO

        DO i=1,3
          DO j=1,now%nyvar_m
            intcnt = intcnt+1
            intbuf(intcnt) = now%ilist_ml(i,j)
          ENDDO
          DO j=1,now%nyvar_p
            intcnt = intcnt+1
            intbuf(intcnt) = now%ilist_pl(i,j)
          ENDDO
          DO j=1,now%nyvar_z
            intcnt = intcnt+1
            intbuf(intcnt) = now%ilist_zl(i,j)
          ENDDO
          DO j=1,now%nyvar_s
            intcnt = intcnt+1
            intbuf(intcnt) = now%ilist_sl(i,j)
          ENDDO
          DO j=1,now%nyvar_c
            intcnt = intcnt+1
            intbuf(intcnt) = now%ilist_c(i,j)
          ENDDO

        ENDDO

        intcnt = intcnt+1
        IF( now%lwrite_const ) THEN
          intbuf(intcnt) = 1
        ELSE
          intbuf(intcnt) = 0
        ENDIF


        IF (ASSOCIATED(now%next) ) THEN
          now => now%next
        ELSE
          EXIT gribout_loop_snd
        ENDIF

      ENDDO gribout_loop_snd


      realcnt = 0
      IF     (vcoord%ivctype == 1) THEN
        DO i=1,ke1
          realcnt = realcnt+1
          realbuf(realcnt) = vcoord%sigm_coord(i)
        ENDDO
      ELSEIF ( ANY( vcoord%ivctype == (/2,3,4/) ) ) THEN
        DO i=1,ke1
          realcnt = realcnt+1
          realbuf(realcnt) = vcoord%vert_coord(i)
        ENDDO
      ENDIF

      realcnt = realcnt+1
      realbuf(realcnt) = hhl(1,1,1)   

      realcnt = realcnt+1
      realbuf(realcnt) = refatm%p0sl

      realcnt = realcnt+1
      realbuf(realcnt) = refatm%t0sl
 
      realcnt = realcnt+1
      realbuf(realcnt) = refatm%dt0lp

      ! two variables of the reference atmosphere 2
      realcnt = realcnt+1
      realbuf(realcnt) = refatm%delta_t

      realcnt = realcnt+1
      realbuf(realcnt) = refatm%h_scal

      realcnt = realcnt+1
      realbuf(realcnt) = vcoord%vcflat 

      DO i=1,ke_soil + 1
        realcnt = realcnt+1
        realbuf( realcnt ) = czhls(i)
      ENDDO

    ENDIF

    CALL MPI_BCAST(intcnt,1,imp_integers,0,intercomm_asynio,ierr) 


    IF( my_world_id >= nziope0 ) THEN
      ALLOCATE ( intbuf(intcnt)   , STAT=istat )
      IF( istat /= 0) THEN
        izerror = 33340
        yzerrmsg = 'ERROR allocating buffer of integers for IO distribution'
        RETURN
      ENDIF
 
    ENDIF 


    CALL distribute_values  (intbuf , intcnt, 0, imp_integers, intercomm_asynio, ierr)
    IF( ierr /= 0) THEN
        izerror = 33343
        yzerrmsg = 'ERROR distributing values to IO sector'
        RETURN
    ENDIF
 

    CALL distribute_values  ( realbuf , rbuflen , 0, imp_reals, intercomm_asynio, ierr)
    IF( ierr /= 0) THEN
        izerror = 33343
        yzerrmsg = 'ERROR distributing values to IO sector'
        RETURN
    ENDIF
  
    IF( my_world_id >= nziope0 ) THEN

      vcoord%ivctype = intbuf(1)
      refatm%irefatm = intbuf(2)


      now => root

      intcnt=2
 
      gribout_loop_rec: DO

        DO i=1,3
          DO j=1,now%nyvar_m
            intcnt = intcnt+1
            now%ilist_ml(i,j) = intbuf(intcnt) 
          ENDDO
          DO j=1,now%nyvar_p
            intcnt = intcnt+1
            now%ilist_pl(i,j) = intbuf(intcnt) 
          ENDDO
          DO j=1,now%nyvar_z
            intcnt = intcnt+1
            now%ilist_zl(i,j) = intbuf(intcnt) 
          ENDDO
          DO j=1,now%nyvar_s
            intcnt = intcnt+1
            now%ilist_sl(i,j) = intbuf(intcnt) 
          ENDDO
          DO j=1,now%nyvar_c
            intcnt = intcnt+1
            now%ilist_c(i,j) = intbuf(intcnt) 
          ENDDO

        ENDDO

        intcnt = intcnt+1
        IF ( intbuf(intcnt) == 1 ) THEN
          now%lwrite_const = .TRUE. 
        ELSE
          now%lwrite_const = .FALSE. 
        ENDIF

        IF (ASSOCIATED(now%next) ) THEN
          now => now%next
        ELSE
          EXIT gribout_loop_rec
        ENDIF

      ENDDO gribout_loop_rec


!US      ALLOCATE ( vcoord (ke+1)      , STAT=istat )
      ALLOCATE (czhls(0:ke_soil+1) , STAT=istat )
      ALLOCATE (hhl (1,1,1), STAT=istat )
  

      realcnt = 0
      IF     (vcoord%ivctype == 1) THEN
        DO i=1,ke1
          realcnt = realcnt+1
          vcoord%sigm_coord(i) = realbuf(realcnt)
        ENDDO
      ELSEIF ( ANY( vcoord%ivctype == (/2,3,4/) ) ) THEN
        DO i=1,ke1
          realcnt = realcnt+1
          vcoord%vert_coord(i) = realbuf(realcnt)
        ENDDO
      ENDIF

      realcnt = realcnt+1
      hhl(1,1,1) =  realbuf(realcnt)

      realcnt = realcnt+1
      refatm%p0sl =  realbuf(realcnt)

      realcnt = realcnt+1
      refatm%t0sl =  realbuf(realcnt)

      realcnt = realcnt+1
      refatm%dt0lp =  realbuf(realcnt)

      ! two variables of the reference atmosphere 2
      realcnt = realcnt+1
      refatm%delta_t =  realbuf(realcnt)

      realcnt = realcnt+1
      refatm%h_scal =  realbuf(realcnt)

      realcnt = realcnt+1
      vcoord%vcflat =  realbuf(realcnt)

      DO i=1,ke_soil + 1
        realcnt = realcnt+1
        czhls(i) = realbuf( realcnt )
      ENDDO


    ENDIF

    DEALLOCATE( intbuf )

    DEALLOCATE( realbuf )


    ! the following block will distribute the vartab table to all the I/O PE
    IF( my_world_id == 0) THEN

      ! size of ar_des data structure
      ar_des_size = 944_iintegers + 1024_iintegers ! 944 is estimated size + 1024 extra margin
      vartab_bufsize = SIZE(var,1)*SIZE(var,2)*SIZE(var,3)*ar_des_size

      ! allocate memory for buffer
      ALLOCATE(vartab_buf(vartab_bufsize), STAT=istat)
      IF( istat /= 0) THEN
        izerror = 33340
        yzerrmsg = 'ERROR allocating buffer for vartab distribution'
        RETURN
      ENDIF
 
      ! transfer to buffer and check size
      var_size = SIZE(TRANSFER(var,vartab_buf))
      IF (var_size > vartab_bufsize) THEN
        yzerrmsg = 'ERROR: size of buffer to distribute vartab table < effective size'
        izerror = 33341
        RETURN
      ENDIF

      ! transfer into linear byte (i.e. character) buffer
      vartab_buf(1:var_size) = TRANSFER(var,vartab_buf)

    ENDIF

    CALL MPI_BCAST(var_size,1,imp_integers,0,intercomm_asynio,ierr) 
    IF( ierr /= 0) THEN
        izerror = 33343
        yzerrmsg = 'ERROR broadcast size of buffer'
        RETURN
    ENDIF
 

    IF( my_world_id >= nziope0 ) THEN
      ! allocate memory for buffer
      ALLOCATE(vartab_buf(var_size), STAT=istat)
      IF( istat /= 0) THEN
        izerror = 33340
        yzerrmsg = 'ERROR allocating buffer for vartab distribution'
        RETURN
      ENDIF
 
    ENDIF
  
    CALL distribute_vartab (vartab_buf, var_size, 0, imp_integers, intercomm_asynio, ierr)
    IF( ierr /= 0) THEN
        izerror = 33343
        yzerrmsg = 'ERROR distributing values to IO sector'
        RETURN
    ENDIF
 

    IF( my_world_id >= nziope0 ) THEN
      IF ( .NOT. ALLOCATED( var ) ) THEN
        ! Allocate the array var
        istat = 0
        ALLOCATE(var(4,0:255,num_gribtabs), STAT=istat)
        IF(istat /= 0) THEN
          yzerrmsg = 'Error allocating var from IO PE'
          izerror = 33340
          RETURN
        ENDIF
      ENDIF
      ! transfer back from linear buffer into vartab (this would be done on I/O PE)
      var = RESHAPE(TRANSFER(vartab_buf,var),shape(var))
    ENDIF

    ! free memory
    DEALLOCATE(vartab_buf)

  ENDIF

END SUBROUTINE distribute_values_asynio


SUBROUTINE close_netcdf_file (nudat,yerrmsg,ierror)

!------------------------------------------------------------------------------
!
! Description:
!   Close a netcdf file. 
!
!------------------------------------------------------------------------------
!
! Scalar arguments with intent(in):
INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
  nudat ! internal file descriptor

CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
  yerrmsg     ! string for error messages

INTEGER  (KIND=iintegers), INTENT(OUT)   ::  &
  ierror      ! error status

!------------------------------------------------------------------------------

! Local scalars:
INTEGER  (KIND=iintegers)           :: &
  implcode                               ! Error variable for MPI


! End of Header
!------------------------------------------------------------------------------

  implcode = 0
  yerrmsg = '     '

  IF( izdebug > 5) THEN
    PRINT *,'CLOSING netcdf FILE'
  ENDIF

  implcode = nf90_close(nudat)

  IF (implcode /= NF90_NOERR) THEN
    PRINT *, 'error closing file',TRIM(NF90_strerror(implcode))
    ierror  = 3
    yerrmsg = 'error closing NetCDF file'
    RETURN
  ENDIF


END SUBROUTINE close_netcdf_file

#endif

END MODULE netcdf_io 
