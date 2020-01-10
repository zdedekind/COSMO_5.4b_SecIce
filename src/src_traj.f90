!+ Module providing routines for online trajectory calculation
!------------------------------------------------------------------------------

MODULE src_traj

!------------------------------------------------------------------------------
!
! Description:
!   This module performs the iterative Euler time-step to calculate the new 
!   trajectory position.
!
! Current Code Owner: IAC ETH Zurich,Stephan Pfahl 
!  phone:  +41 44 632 93 65
!  email:  stephan.pfahl@env.ethz.ch
!
! History:
! Version      Date       Name
! ----------   ---------- ----
! V5_1         2014-11-28 Annette Miltenberger/Anne Roches
!  Initial release
! V5_2         2015-05-21 Ulrich Schaettler, Tobias Selz
!  Introduced some ifdef NETCDF directives to be able to compile without NetCDF
!  Corrected dimensions for MPI_WAITALL status variable in SR comm_traj_2buf
!   (MPI_STATUS_SIZE was missing)
!  Check for any trajectory in output (SR: output_traj), not only for the first one
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Modules used:
!
USE data_parameters ,   ONLY :   &
  intgribf,                & ! KIND-type parameter for fortran files in the grib library
  wp,                      & ! KIND-type parameters for real variables
  iintegers                  ! KIND-type parameter for standard integer variables

USE data_modelconfig,   ONLY :   &
  idt_qv,                  & ! index for QV in the tracer structure
  idt_qc,                  & ! index for QC in the tracer structure
  idt_qi,                  & ! index for QI in the tracer structure
  idt_qr,                  & ! index for QR in the tracer structure
  idt_qs,                  & ! index for QS in the tracer structure
  idt_qg,                  & ! index for QG in the tracer structure
  ie_tot,                  & ! number of grid points in zonal direction
  je_tot,                  & ! number of grid points in meridional direction
  ke_tot,                  & ! number of grid points in vertical direction
  startlon_tot,            & ! longitude at lower left corner of domain 
  startlat_tot,            & ! latitude at lower left corner of domain 
  dlon,                    & ! grid spacing in x in degree
  dlat,                    & ! grid spacing in y in degree
  eddlon,                  & ! 1 / dlon
  eddlat,                  & ! 1 / dlat
  pollon,                  & ! longitude of the rotated north pole (in degrees, E>0)
  pollat,                  & ! latitude  of the rotated north pole (in degrees, N>0)
  dt,                      & ! long time-step
  ie, je, ke                 ! number of grid points in zonal/meridio./vertic. direction

USE data_runcontrol,    ONLY :   &
  lmetr,                   & ! lartif_data=.TRUE.:  with metric terms
                             !            =.FALSE.: or without metric terms
  lcori_deep,              & ! if =.TRUE.: take cos(phi) coriolis terms into account
  lprintdeb_all,           & ! if debug should be done on all PEs or only on PE0
  idbg_level,              & ! control the verbosity of debug output
  nuspecif,                & ! unit number for file YUSPECIF
  ntstep,                  & ! actual time step
  hstart,                  & ! start of the forecast in full hours
  hstop,                   & ! end of the forecast in hours
  nstart,                  & ! first time step of the forecast
  nstop,                   & ! last time step of the forecast
  nnew,                    & ! time level after the dynamics 
  nnow,                    & ! time level before the dynamics
  itype_calendar,          & ! for specifying the calendar used
  lperi_x,                 & ! if lartif_data=.TRUE.: periodic boundary conditions in x-dir.
  lperi_y,                 & ! if lartif_data=.TRUE.: periodic boundary conditions in y-dir.
  l2dim                      ! 2dimensional model version

USE data_fields,        ONLY :   &
  fc,                      & ! coriolis-parameter                          ( 1/s )
  fccos,                   & ! horizontal coriolis-parameter               ( 1/s )
  acrlat,                  & ! 1 / ( crlat * radius of the earth )         ( 1/m )
  tgrlat,                  & ! tangens of transformed latitude               --
  hhl,                     & ! geometrical height of half levels           ( m   )
  u, v, w,                 & ! zonal/meridional/vertical wind speed        ( m/s )
  pp, p0,                  & ! deviation from and reference pressure       ( pa  )
  t                          ! temperature                                 ( K   ) 

USE data_parallel,      ONLY :   &
  isubpos,                 & ! positions of the subdomains in the total domain
  logbuf, intbuf,          & ! buffers for logical and integer
  realbuf, charbuf,        & ! buffers for real and character
  my_cart_id,              & ! rank of this subdomain in the cartesian communicator
  my_cart_pos,             & ! position of this subdomain in the cartesian grid
  icomm_cart,              & ! communicator for the virtual cartesian topology
  icomm_world,             & ! communicator that belongs to igroup_world
  nprocx,                  & ! number of processors in x-direction
  nprocy,                  & ! number of processors in y-direction
  nproc,                   & ! total number of processors
  nboundlines,             & ! number of overlapping boundary lines of the subdomains
  imp_reals, imp_integers, & ! datatypes for MPI
  imp_character,           & ! datatype for MPI
  num_compute                ! number of compute PEs

USE data_io,            ONLY :   &
  ydate_ini,               & ! start of the forecast
  ydir_restart_out,        & ! directory of restart data
  ytunit_restart,          & ! unit of timescale
  max_gribrep,             & ! maximum number of GRIB tables in LM variable table
  num_gribtabs,            & ! maximum number of GRIB tables in LM variable table
  var                        ! array for LM variable table

USE data_constants,     ONLY :   &
  r_earth,                 & ! mean radius of the earth (m)
  b1,b2w,b3,b4w,rdv,       & ! saturation constants
  o_m_rdv,r_d,cp_d,rvd_m_o   ! saturation constants

USE data_traj

USE data_tracer,        ONLY : T_ERR_NOTFOUND

USE grid_metrics_utilities,        ONLY :  &
  sqrtg_r_s,    & ! 1 / square root of G at scalar points       ( 1/m )
  dzeta_dlam,   & ! d zeta / d lambda (for constant phi,    z)
                  ! at the scalar position                      ( 1   )
  dzeta_dphi,   & ! d zeta / d phi    (for constant lambda, z)
                  ! at the scalar position                      ( 1   )
  wgtfac          ! weighting factor for vertical interpolation ( 1 )

USE environment,        ONLY :   &
  model_abort

USE io_utilities,       ONLY :   &
  make_fn, get_utc_date

USE parallel_utilities, ONLY :   &
  distribute_values

USE pp_utilities,       ONLY :   &
  calrelhum, potential_vorticity_rho

USE numeric_utilities,  ONLY :   &
  curl, calc_Theta_Tppp

USE meteo_utilities,    ONLY :   &
  calrho

USE src_tracer,         ONLY :   &
  trcr_get, trcr_errorstr

#ifdef NETCDF
USE netcdf,             ONLY:    &
   nf90_def_dim, nf90_def_var, nf90_enddef,                                   &                         
   nf90_put_att, nf90_put_var, nf90_inq_varid,                                &
   nf90_clobber, nf90_close, nf90_create, nf90_open,                          &
   nf90_set_fill,nf90_NOFILL, NF90_WRITE,                                     &
   nf90_64BIT_OFFSET, nf90_GLOBAL, nf90_FLOAT, nf90_strerror, nf90_NOERR
#endif

!------------------------------------------------------------------------------

IMPLICIT NONE

!------------------------------------------------------------------------------

INCLUDE "mpif.h"

!------------------------------------------------------------------------------

PRIVATE

! Local variables

INTEGER (KIND=iintegers)                                :: &
  my_cart_8neigh(8), &     ! 8 neighbor PEs of each compute PE 
  ntracevar_traj           ! number of traced variables along trajectories  [-]

INTEGER (KIND=iintegers)                                :: &
  TimeDimId                ! NetCDF ID for the time dimension

INTEGER (KIND=iintegers)                                :: &
  n_files_tot              ! number of output files


INTEGER(KIND=iintegers), ALLOCATABLE                    :: &
  ncID     (:),    &       ! ID of the netcdf files
  IdDimId  (:),    &       ! NetCDF ID for the trajectory IDs dimension         
  varLonId (:),    &       ! NetCDF ID for the longitudes
  varLatId (:),    &       ! NetCDF ID for the latitudes
  varZId   (:),    &       ! NetCDF ID for the height
  varTimeId(:),    &       ! NetCDF ID for the time
  outstep  (:)             ! number of output steps 

INTEGER(KIND=iintegers), ALLOCATABLE                    :: &
  varId  (:,:),    &       ! NetCDF ID for the traced variables
  id_file(:,:)             ! index of the 1st and last traj. in each file 

REAL (KIND=wp),     ALLOCATABLE                         :: &
  traj_startinfo(:,:)   ,& ! array containing the start position of the traj.
  traj_startinfo_time(:)   ! array containing the start times of the traj.

! Declare public entities

PUBLIC :: &
  organize_traj, organize_traj_restart

!==============================================================================
! Module procedures
!==============================================================================

CONTAINS

!==============================================================================
!+ organize trajectory computation
!------------------------------------------------------------------------------

SUBROUTINE organize_traj (yaction, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   This procedure is the driving routine for the online trajectory calculation
!   It is called first at the beginning of the program to read the namelist
!   trajctl and the startfile. It is then called at initialization to 
!   initialize the trajectories and their output. During the simulation,
!   it is then called at each time step in order to integrate the new
!   trajectory position and possibly update the value of the traced
!   variables. A last call is performed at the end of the simulation in order
!   to close the NetCDF output files and clean the memory. 
!
! Method:
!   Depending on the calling action, the namelists are read, the module is
!   initialized, the trajectories are computed and written out or a final
!   cleanup is performed.
!
!==============================================================================

! Parameter list:
CHARACTER (LEN= *),         INTENT(IN)            ::                          &
  yaction      ! action to be performed

INTEGER   (KIND=iintegers), INTENT(OUT)           ::                          &
  ierror       ! error status

CHARACTER (LEN= *),         INTENT(OUT)           ::                          &
  yerrmsg      ! error message


! Local variables: 
INTEGER (KIND=iintegers)                  ::  &
  i, k,              & ! loop indices
  izdebug,           & ! debug level 
  nuin,              & ! unit number for Namelist INPUT file
  iz_proc,           & ! PE index
  izerr,             & ! error code
  izerrstat            ! error status  

CHARACTER (LEN=10)                        ::  &
  yzinput              ! namelist INPUT file


INTEGER (KIND=iintegers)                  ::  &
  nztraj_buff   (8)    ! number of trajectories to be sent to the 8 neigh. PEs


REAL (KIND=wp),     ALLOCATABLE           ::  & 
  traj_buff(:,:)       ! sendbuffer for trajectory communication to neigh. PEs


!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------
!- Begin Subroutine organize_traj
!------------------------------------------------------------------------------

  izerrstat = 0_iintegers
  ierror    = 0_iintegers
  yerrmsg   = '    '
  izerr     = 0_iintegers


  ! Initalize, whether debug output shall be done
  IF (lprintdeb_all) THEN
    izdebug = idbg_level
  ELSE
    IF (my_cart_id == 0_iintegers) THEN
      izdebug = idbg_level
    ELSE
      izdebug = 0_iintegers
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
!- Section 1: Input of the Namelist and read of the start file
!------------------------------------------------------------------------------

  IF (yaction == 'input') THEN

    !--------------------------------------------------------------------------
    ! Section 1.1: Read namelist
    !--------------------------------------------------------------------------

    ! Open NAMELIST-INPUT file
    IF (my_cart_id == 0_iintegers) THEN
      IF (idbg_level > 0_iintegers) THEN
        PRINT *,'    INPUT OF THE NAMELISTS FOR TRAJECTORY SIMULATION'
      ENDIF
      yzinput = 'INPUT_TRAJ'
      nuin    =  1_iintegers
      OPEN(nuin, FILE=yzinput, FORM=  'FORMATTED', STATUS='UNKNOWN',          &
           IOSTAT=izerrstat)
      IF(izerrstat /= 0) THEN
        yerrmsg  = ' ERROR    *** Error while opening file INPUT_TRAJ *** '
        ierror   = 10000_iintegers
        RETURN
      ENDIF
    ENDIF

    ! Read NAMELIST group TRAJCTL
    CALL input_trajctl (izdebug, nuspecif, nuin, izerrstat)

    IF (izerrstat /= 0_iintegers) THEN
      yerrmsg = 'ERROR    *** while reading NAMELIST Group /TRAJCTL/ ***'
      ierror  = 10003_iintegers
      RETURN
    ENDIF


    IF (my_cart_id == 0_iintegers) THEN
      ! Close file for input of the NAMELIST
      CLOSE (nuin, STATUS='KEEP')
    ENDIF

    !--------------------------------------------------------------------------
    ! Section 1.2: Read start file for trajectories
    !--------------------------------------------------------------------------

    ! Open start file
     
    IF (my_cart_id == 0_iintegers) THEN  
      IF (izdebug > 0_iintegers) THEN
        PRINT *,'    INPUT OF STARTFILE FOR TRAJECTORY CALCULATION'
      ENDIF
      nuin    = 1_iintegers
      OPEN(nuin, FILE=startfile_traj, FORM=  'FORMATTED', STATUS='OLD',       &
           ACTION='READ', IOSTAT=izerrstat)
      IF(izerrstat /= 0) THEN
        yerrmsg  = ' ERROR    *** Error while opening file STARTFILE *** '
        ierror   = 10006_iintegers
        RETURN
      ENDIF
    ENDIF

    ! Read start file

    CALL read_startf_traj(izdebug, nuin, izerrstat, yerrmsg)
    IF (izerrstat /= 0_iintegers) THEN
      ierror = izerrstat
      RETURN
    ENDIF

    IF (my_cart_id == 0_iintegers) THEN
      ! Close file for input of the NAMELIST
      CLOSE (nuin, STATUS='KEEP')
    ENDIF


!------------------------------------------------------------------------------
!- Section 2: Initialize trajectory calculation
!------------------------------------------------------------------------------

  ELSEIF (yaction == 'init') THEN

    !--------------------------------------------------------------------------
    ! Section 2.1: Precomputations
    !--------------------------------------------------------------------------

    IF (my_cart_id == 0_iintegers) THEN
      IF (izdebug > 0_iintegers) THEN
        PRINT *,'*** INITIALIZATION OF TRAJECTORY CALCULATIONS ***'
      ENDIF
    ENDIF

    ! Allocate the variable holding the status of the trajectory
    ALLOCATE (istat_traj (num_start_pt_tot*num_start_t), STAT=izerrstat)
    IF (izerrstat /= 0_iintegers) THEN
      ierror = 10009_iintegers
      yerrmsg = 'allocation of istat_traj failed'
      RETURN
    ENDIF
    istat_traj(:) = idead

    ! Allocate the trajectory structure
    !   dim 1:  # of startpoints * # of start times (== # trajectories) 
    !   dim 2:  # of spatial coordinates (3)      
    ALLOCATE (traj (num_start_pt_tot*num_start_t,3), STAT=izerrstat)
    IF (izerrstat /= 0_iintegers) THEN
      ierror = 10018_iintegers
      yerrmsg = 'allocation of space for traj failed'
      RETURN
    ENDIF
    traj(:,:) = -999.0_wp

    ! Allocate and calculate geometrical height of full levels
    ALLOCATE (hfl (ie, je, ke), STAT=izerrstat)
    IF (izerrstat /= 0_iintegers) THEN
      ierror = 10012_iintegers
      yerrmsg = 'allocation of space for hl failed'
      RETURN
    ENDIF
    DO k = 1,ke
      hfl(:,:,k) = 0.5_wp * (hhl(:,:,k) + hhl(:,:,k+1))
    ENDDO

    ! Allocate the data structure holding the variables to trace along
    ! the trajectories
    ALLOCATE (traj_trace(ntracevar_traj), STAT=izerrstat)
    IF (izerrstat /= 0_iintegers) THEN
      ierror = 10015_iintegers
      yerrmsg = 'allocation of traj_trace failed'
      RETURN
    ENDIF

    ! Define the units, name, etc. of the variables to trace
    CALL define_traj_trace

    ! Define for each PE the 8 surrounding PEs (horizontal domaine decomp.)
    ! i.e. the 8 PEs a trajectory could potentially travel to
    IF (num_compute > 1_iintegers) THEN
      CALL define_8neigh_proc(my_cart_8neigh,ierror, yerrmsg)
      IF (ierror /= 0_iintegers) RETURN
    ENDIF


    ! Define the total number of trajectories
    ntraj = num_start_pt_tot * num_start_t


    !--------------------------------------------------------------------------
    ! Section 2.2: Define the trajectory spatial position and status at start
    !--------------------------------------------------------------------------

    IF ( nstart == 0_iintegers ) THEN

      ! Intialize the trajectory spatial positions using the starting positions 
      ! read from the startfile
      traj(1:num_start_pt_tot,1) = traj_startinfo(1:num_start_pt_tot,1)
      traj(1:num_start_pt_tot,2) = traj_startinfo(1:num_start_pt_tot,2)
      traj(1:num_start_pt_tot,3) = traj_startinfo(1:num_start_pt_tot,3)
      
      IF (num_start_t > 1_iintegers) THEN
        DO i = 1, num_start_t-1
          traj(num_start_pt_tot*i+1_iintegers:num_start_pt_tot*(i+1),:) = traj(1:num_start_pt_tot,:)
        ENDDO
      ENDIF
      
      IF (ANY(istart_mode_traj == (/1,2/))) THEN
        istat_traj(1:num_start_pt_tot) = nstart_traj_gp(1)
        IF (num_start_t > 1_iintegers) THEN
          DO i = 1, num_start_t-1
            istat_traj(num_start_pt_tot*i+1:num_start_pt_tot*(i+1)) = nstart_traj_gp(i+1)
          ENDDO
        ENDIF
      ELSEIF (istart_mode_traj == 3_iintegers) THEN
        istat_traj(1:num_start_pt_tot) = nstart_traj_iv(1:num_start_pt_tot)
      ENDIF

    !--------------------------------------------------------------------------
    ! Section 2.3: Read the trajectory spatial position at restart
    !--------------------------------------------------------------------------

    ELSE
      !Read restart file
      CALL organize_traj_restart('read')

    ENDIF
    
    !--------------------------------------------------------------------------
    ! Section 2.4: Open NetCDF and initialize NetCDF files
    !--------------------------------------------------------------------------
    
    ! All NetCDF output files are created, global attributes, variables and 
    ! their attributes are defined, and a first value for the variable
    ! time (-1) is written out
    CALL output_traj_init(izdebug, ierror, yerrmsg)
    IF (ierror /= 0_iintegers) THEN
      RETURN
    ENDIF

    !--------------------------------------------------------------------------
    ! Section 2.5: Final operations 
    !--------------------------------------------------------------------------
    
    ! Deallocate variable required for initialization
    DEALLOCATE (traj_startinfo, STAT=izerr)
    IF (izerr /= 0_iintegers) THEN
      ierror  = 10021_iintegers
      yerrmsg = 'Deallocation of traj_startinfo variable for trajectories failed.'
      RETURN
    ENDIF

!------------------------------------------------------------------------------
!- Section 3: Computation of trajectories
!------------------------------------------------------------------------------

  ELSEIF (yaction == 'compute') THEN 

    IF (my_cart_id == 0_iintegers) THEN
      IF (izdebug > 0_iintegers) THEN
        PRINT *,'      TRAJECTORY COMPUTATION'
      ENDIF
    ENDIF


    ! For the whole trajectory simulation period
    IF ((ntstep >= nstartfirst_traj) .AND. (ntstep <= nstop_traj)) THEN

      IF (num_compute > 1_iintegers) THEN
        ! Allocate space for trajectory buffer
        !  dim 1: 8 neighbor PEs for each compute PE
        !  dim 2: for each traj, save traj_idx and (x,y,z) -> 4 positions/traj
        ALLOCATE (traj_buff (ntraj*4, 8), STAT=izerrstat)  
        IF (izerrstat /= 0_iintegers) THEN
          ierror  = 10024_iintegers    
          yerrmsg = 'allocation of space for traj_buff failed'
          RETURN
        ENDIF
        traj_buff (:,:)= 0.0_wp
      ENDIF

      nztraj_buff(:) = 0_iintegers


      DO i = 1, ntraj


        ! Check whether computation of a trajectory should start
        IF (istat_traj(i) == ntstep) istat_traj(i) = ialive

        ! Perform trajectory computation if trajectory is on this PE
        ! and has to be computed
        IF (istat_traj(i) == ialive) THEN
          iz_proc = ipe_location(traj(i,1), traj(i,2), startlon_tot,          &
                                 startlat_tot, isubpos(my_cart_id,1),         &
                                 isubpos(my_cart_id,2), dlon, dlat, idxpe_n,  &
                                 idxpe_ne, idxpe_e,  idxpe_se, idxpe_s,       &
                                 idxpe_sw, idxpe_w, idxpe_nw, idxpe_this,     &
                                 nboundlines, ie, je)

          IF(iz_proc == idxpe_this) THEN

            ! Forward integration to compute the new position of the trajectory
            CALL euler_ts(i,izdebug) 
            
            IF (num_compute > 1_iintegers) THEN
              ! Check if trajectory is still within the PE domain and if not prepare 
              ! sending its position to the relevant PE
              CALL prep_comm_traj(i,nztraj_buff,traj_buff)

            ENDIF

          ENDIF 
        ENDIF
      ENDDO

      IF (num_compute > 1_iintegers) THEN
        ! each PE send nztraj_buff(j) trajectories to its neighbor j, and for 
        ! each trajectory the id and (x,y,z) info are sent
        nztraj_buff(:) = nztraj_buff(:) * 4_iintegers

        ! Communication with other PEs
        CALL do_comm_traj(nztraj_buff,traj_buff,my_cart_8neigh,ierror,yerrmsg)
        IF (ierror /= 0_iintegers) RETURN
      ENDIF
 
      ! In case of output timestep, interpolate fields to trajectory location,
      ! gather the information to process 0 and write it to netcdf-file
      IF ((MODULO(ntstep,ninc_out_traj) == 0_iintegers)) THEN
        CALL output_traj(ntstep,dt,ierror,yerrmsg)
        IF (ierror /= 0_iintegers) RETURN
      ENDIF

      IF (num_compute > 1_iintegers) THEN
        DEALLOCATE (traj_buff, STAT=izerr)
        IF (izerr /= 0_iintegers) THEN
          ierror  = izerr
          yerrmsg = 'Deallocation of traj_buff failed.'
          RETURN
        ENDIF
      ENDIF

    ENDIF

!------------------------------------------------------------------------------
!- Section 4: Finalize trajectory computation
!------------------------------------------------------------------------------

  ELSEIF (yaction == 'finalize') THEN

#ifdef NETCDF
    ! Close netcdf files
    IF (my_cart_id == 0_iintegers) THEN
      DO k=1,n_files_tot
        izerrstat = nf90_close(ncID(k))
        IF (izerrstat /= nf90_noerr) THEN
          ierror  = 10027_iintegers 
          yerrmsg = TRIM(nf90_strerror(izerrstat))
          RETURN
        ENDIF
      ENDDO
    ENDIF
#endif

    ! Deallocate variables for trajectory calculation
    DEALLOCATE(traj_trace,    STAT=izerrstat); izerr=izerr+izerrstat
    DEALLOCATE(traj,          STAT=izerrstat); izerr=izerr+izerrstat
    DEALLOCATE(istat_traj,    STAT=izerrstat); izerr=izerr+izerrstat
    DEALLOCATE(hfl,           STAT=izerrstat); izerr=izerr+izerrstat
    IF (izerr /= 0_iintegers) THEN
      ierror  = 10030_iintegers
      yerrmsg = 'Deallocation of variables for trajectories failed.'
      RETURN
    ENDIF
    IF (my_cart_id == 0_iintegers) THEN
      DEALLOCATE(outstep,     STAT=izerrstat); izerr=izerr+izerrstat
      DEALLOCATE(varId,       STAT=izerrstat); izerr=izerr+izerrstat
      DEALLOCATE(ncId,        STAT=izerrstat); izerr=izerr+izerrstat
      DEALLOCATE(IdDimId,     STAT=izerrstat); izerr=izerr+izerrstat
      DEALLOCATE(varLonId,    STAT=izerrstat); izerr=izerr+izerrstat
      DEALLOCATE(varLatId,    STAT=izerrstat); izerr=izerr+izerrstat
      DEALLOCATE(varZId,      STAT=izerrstat); izerr=izerr+izerrstat
      DEALLOCATE(varTimeId,   STAT=izerrstat); izerr=izerr+izerrstat
      DEALLOCATE(id_file,     STAT=izerrstat); izerr=izerr+izerrstat
      IF (izerr /= 0_iintegers) THEN
        ierror  = 10033_iintegers
        yerrmsg = 'Deallocation of NetCDF variables for trajectories failed.'
        RETURN
      ENDIF
    ENDIF

!------------------------------------------------------------------------------
!- Section X: All other actions are wrong
!------------------------------------------------------------------------------

  ELSE

    ierror  = 10036_iintegers
    yerrmsg = 'ERROR *** No valid action for organize_traj ***'

  ENDIF


!------------------------------------------------------------------------------
!- End of Subroutine organize_traj
!------------------------------------------------------------------------------

END SUBROUTINE organize_traj

!==============================================================================
!+ Module procedure for the definition of the elements of traj_trace
!------------------------------------------------------------------------------

SUBROUTINE define_traj_trace

!------------------------------------------------------------------------------
!
! Description: 
!   A check is performed that the variable required by the user to be traced 
!   along the trajectories is present in the setup_vartab table, that it is 
!   really a 3D spatial variable (x,y,z) and that it is an instantaneous 
!   variable. The name, level type, time range indicator, rank, definition 
!   status, units, short name and long name are defined and the pointers
!   are set.
!              
! Method:
!  Comparison against the setup_vartab table (using variable name) and filling
!  of the structure traj_trace for all traced variables using the information
!  available for the variable in the setup_vartab table.
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments

 
! Other local variables
INTEGER (KIND=iintegers)               :: &
  k,                 &     ! loop index
  iz1, iz2, iz3,     &     ! indices for location in vartab
  izerror                  ! error code

CHARACTER(LEN=80)                      :: &
  yzerror,           &     ! error message
  yzroutine                ! subroutine name

INTEGER (KIND=iintegers), PARAMETER    :: &
  notdef = -999_intgribf   ! initialization value         

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE define_traj_trace
!------------------------------------------------------------------------------

  izerror   = 0_iintegers
  yzroutine = 'define_traj_trace'
  yzerror   = '  '

!------------------------------------------------------------------------------
!- Section 1: Initialization
!------------------------------------------------------------------------------

  traj_trace(:)%name      = ' '
  traj_trace(:)%levtyp    = notdef
  traj_trace(:)%ntri      = notdef
  traj_trace(:)%rank      = notdef
  traj_trace(:)%idef_stat = notdef
  traj_trace(:)%istag     = 0_iintegers
  traj_trace(:)%units     = ' '
  traj_trace(:)%sdname    = ' '
  traj_trace(:)%lgname    = ' '

! Note: the pointer elements have been initialized in data_traj.f90 already

!------------------------------------------------------------------------------
!- Section 2: Assignment
!------------------------------------------------------------------------------ 

  DO k = 1, ntracevar_traj
    DO iz3 = 1,num_gribtabs
      DO iz2 = 0,255
        DO iz1 = 1,max_gribrep
          IF (tracevar_traj(k) == TRIM(var(iz1,iz2,iz3)%name)) THEN
            IF (var(iz1,iz2,iz3)%rank == 4_iintegers) THEN
              IF (var(iz1,iz2,iz3)%ntri == 0_iintegers) THEN
                IF (ANY(var(iz1,iz2,iz3)%levtyp == (/109,110/))) THEN
                  traj_trace(k)%name      = TRIM(tracevar_traj(k))
                  traj_trace(k)%levtyp    = var(iz1,iz2,iz3)%levtyp
                  traj_trace(k)%ntri      = var(iz1,iz2,iz3)%ntri
                  traj_trace(k)%rank      = var(iz1,iz2,iz3)%rank
                  traj_trace(k)%p4        => var(iz1,iz2,iz3)%p4
                  traj_trace(k)%idef_stat = var(iz1,iz2,iz3)%idef_stat
                  traj_trace(k)%units     = var(iz1,iz2,iz3)%units
                  traj_trace(k)%sdname    = var(iz1,iz2,iz3)%standard_name
                  traj_trace(k)%lgname    = var(iz1,iz2,iz3)%long_name
                ELSE
                  izerror = 10040_iintegers
                  yzerror = 'TRAJECTORY module: only variables of levtyp=109,110 can be traced'
                  CALL model_abort(my_cart_id, izerror, yzerror, yzroutine)
                ENDIF
              ELSE
                izerror = 10043_iintegers
                yzerror = 'TRAJECTORY module: only variables of ntri=0 can be traced'
                CALL model_abort(my_cart_id, izerror, yzerror, yzroutine)
              ENDIF
            ELSEIF (var(iz1,iz2,iz3)%rank == 3_iintegers) THEN
              IF (var(iz1,iz2,iz3)%ntri == 0_iintegers) THEN
                IF (ANY(var(iz1,iz2,iz3)%levtyp == (/109,110/))) THEN
                  traj_trace(k)%name      = TRIM(tracevar_traj(k))
                  traj_trace(k)%levtyp    = var(iz1,iz2,iz3)%levtyp
                  traj_trace(k)%ntri      = var(iz1,iz2,iz3)%ntri
                  traj_trace(k)%rank      = var(iz1,iz2,iz3)%rank
                  traj_trace(k)%p3        => var(iz1,iz2,iz3)%p3
                  traj_trace(k)%idef_stat = var(iz1,iz2,iz3)%idef_stat
                  traj_trace(k)%units     = var(iz1,iz2,iz3)%units
                  traj_trace(k)%sdname    = var(iz1,iz2,iz3)%standard_name
                  traj_trace(k)%lgname    = var(iz1,iz2,iz3)%long_name
                ELSE
                  izerror = 10046_iintegers
                  yzerror = 'TRAJECTORY module: only variables of levtyp=109,110 can be traced'
                  CALL model_abort(my_cart_id, izerror, yzerror, yzroutine)
                ENDIF
              ELSE
                izerror = 10049_iintegers
                yzerror = 'TRAJECTORY module: only variables of ntri=0 can be traced'
                CALL model_abort(my_cart_id, izerror, yzerror, yzroutine)
              ENDIF
            ELSE
              izerror = 10052_iintegers
              yzerror = 'TRAJECTORY module: only variables of rank=3,4 can be traced'
              CALL model_abort(my_cart_id, izerror, yzerror, yzroutine)
            ENDIF            
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    IF (traj_trace(k)%name == ' ') THEN
      izerror = 10055_iintegers
      yzerror = 'TRAJECTORY module: Variable '                                &
                //TRIM(tracevar_traj(k))//' not found in the GRIB table.'
      CALL model_abort(my_cart_id, izerror, yzerror, yzroutine)
    ENDIF
  ENDDO

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE define_traj_trace

!==============================================================================
!+ Module procedure for the input of NAMELIST trajctl
!------------------------------------------------------------------------------

SUBROUTINE input_trajctl (idebug, nuspecif, nuin, ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine organizes the input of the NAMELIST-group trajctl. 
!   The group trajctl contains variables for the organization of online 
!   trajectory calculation.
!
! Method:
!   All variables are initialized with default values and then read in from
!   the file INPUT_TRAJ.
! 
!   The input values are checked for errors and for
!   consistency. If wrong input values are detected the program prints
!   an error message. The program is not stopped in this routine but an
!   error code is returned to the calling routine that aborts the program after
!   reading in all other namelists.
!   In parallel mode, the variables are distributed to all nodes with the
!   environment-routine distribute_values.
!   Both, default and input values are written to the file YUSPECIF
!   (specification of the run).
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
  idebug,       & ! Debug level
  nuspecif,     & ! Unit number for protocolling the task
  nuin            ! Unit number for Namelist INPUT file

INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
  ierrstat        ! error status variable


! Variables for default values
INTEGER   (KIND=iintegers)   ::          & 
  ninc_out_traj_d,                       & ! time interval for trajectory output
  nstop_traj_d,                          & ! final time step of trajectory computation
  nstart_traj_d(nmax_traj_starttime),    & ! start time steps of the trajectories
  ncomb_start_traj_d(3),                 & ! triplet(start,stop,dt) for start time steps
  istart_mode_traj_d                       ! start mode of the trajectories

REAL      (KIND=wp)          ::          &
  hinc_out_traj_d,                       & ! time interval for trajectory output
  hstop_traj_d,                          & ! end of trajectory computation 
  hstart_traj_d(nmax_traj_starttime),    & ! start times of the trajectories 
  hcomb_start_traj_d(3)                    ! triplet (start,stop,dt) for start time


CHARACTER (LEN = 100)        ::          & 
  ydir_traj_d,                           & ! output directory for trajectories
  startfile_traj_d                         ! path to startfile


! Other local variables
INTEGER   (KIND=iintegers)   ::          &
  iznoval,                               & ! no value
  iz_err,                                & ! local error code
  i                                        ! loop index  
  
REAL      (KIND=wp)          ::          &     
  zrnoval                                  ! no value 

REAL(KIND=wp),     PARAMETER::           &
  zeps = 1.E-3_wp                      ! tolerance for real comparison


! Define the namelist group
NAMELIST /trajctl/ hstop_traj, nstop_traj, istart_mode_traj,                  &
                   hstart_traj, nstart_traj, hcomb_start_traj,                &
                   ncomb_start_traj, startfile_traj, ydir_traj,               &
                   hinc_out_traj, ninc_out_traj, tracevar_traj

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE input_trajctl
!------------------------------------------------------------------------------

  ierrstat           = 0_iintegers
  iz_err             = 0_iintegers
  charbuf(:)         = ''
  iznoval            = -999
  zrnoval            = -999.0_wp
  nstart_traj_gp(:)  = 0_iintegers
  tracevar_traj(:)   = 'novar'
  ntracevar_traj     = 0_iintegers

  IF (my_cart_id == 0_iintegers) THEN

!------------------------------------------------------------------------------
!- Section 1: Define defaults
!------------------------------------------------------------------------------

    hstop_traj_d           = 0.0_wp
    nstop_traj_d           = 0_iintegers
    istart_mode_traj_d     = 1_iintegers
    hstart_traj_d(:)       = zrnoval
    nstart_traj_d(:)       = iznoval
    hcomb_start_traj_d(:)  = zrnoval
    ncomb_start_traj_d(:)  = iznoval
    startfile_traj_d       = ' '
    ydir_traj_d            = ' '
    hinc_out_traj_d        = zrnoval
    ninc_out_traj_d        = iznoval

!------------------------------------------------------------------------------
!- Section 2: Initialize variables with defaults
!------------------------------------------------------------------------------

    hstop_traj             = hstop_traj_d
    nstop_traj             = nstop_traj_d
    istart_mode_traj       = istart_mode_traj_d
    hstart_traj  (:)       = hstart_traj_d  (:)
    nstart_traj  (:)       = nstart_traj_d  (:)
    hcomb_start_traj(:)    = hcomb_start_traj_d(:)
    ncomb_start_traj(:)    = ncomb_start_traj_d(:)
    startfile_traj         = startfile_traj_d
    ydir_traj              = ydir_traj_d
    hinc_out_traj          = hinc_out_traj_d
    ninc_out_traj          = ninc_out_traj_d

!------------------------------------------------------------------------------
!- Section 3: Input of the namelist values
!------------------------------------------------------------------------------

    READ (nuin, trajctl, IOSTAT=iz_err)
  ENDIF

  IF (nproc > 1) THEN
    ! distribute error status to all processors
    CALL distribute_values  (iz_err, 1, 0, imp_integers,  icomm_world, ierrstat)
  ENDIF

  IF (iz_err /= 0_iintegers) THEN
    ierrstat = -1_iintegers
    RETURN
  ENDIF

!-----------------------------------------------------------------------------
!- Section 4: Check values for errors and consistency
!------------------------------------------------------------------------------

  IF (my_cart_id == 0_iintegers) THEN

  !----------------------------------------------------------------------------
  !- Section 4.1: Check hstop_traj and nstop_traj
  !----------------------------------------------------------------------------

    ! nstop_traj has priority over hstop_traj if both are specified
    IF (hstop_traj /= hstop_traj_d .AND. nstop_traj == nstop_traj_d) THEN
      nstop_traj = NINT(hstop_traj * 3600.0_wp / dt)
      IF (ABS(nstop_traj - (hstop_traj * 3600.0_wp / dt)) >               &
                    zeps * (hstop_traj * 3600.0_wp / dt)) THEN
        PRINT *, ' WARNING  *** TRAJECTORY module *** '
        PRINT *, '          *** hstop_traj is not a multiple of dt.'
        PRINT *, '          *** nstop_traj is thus based on a rounded hstop_traj.'
      ENDIF
    ELSEIF (hstop_traj == hstop_traj_d .AND. nstop_traj /= nstop_traj_d) THEN
      hstop_traj = nstop_traj * dt / 3600.0_wp
    ELSEIF (hstop_traj /= hstop_traj_d .AND. nstop_traj /= nstop_traj_d) THEN
      PRINT *, ' WARNING  *** TRAJECTORY module *** '
      PRINT *, '          *** hstop_traj and nstop_traj are both specified.'
      PRINT *, '          *** hstop_traj will be ignored.'
      hstop_traj = nstop_traj * dt / 3600.0_wp
    ELSE
      PRINT *, ' WARNING  *** TRAJECTORY module *** '
      PRINT *, '          *** neither hstop_traj nor nstop_traj are specified.'
      PRINT *, '          *** hstop_traj/nstop_traj will be set to hstop/nstop.'
      hstop_traj = hstop
      nstop_traj = nstop   
    ENDIF      

    ! hstop_traj/nstop_traj must be within [hstart/nstart;hstop/nstop]
    IF (nstop_traj < nstart) THEN
      PRINT *, ' ERROR  *** TRAJECTORY module *** '
      PRINT *, '        h/nstop_traj is smaller than h/nstart.'
      ierrstat = 10060_iintegers
      RETURN
    ENDIF
    IF (nstop_traj > nstop) THEN
      PRINT *, ' WARNING  *** TRAJECTORY module *** '
      PRINT *, '          *** hstop_traj/nstop_traj is larger than hstop/nstop.'
      PRINT *, '          *** hstop_traj/nstop_traj are set to hstop/nstop.'
      hstop_traj = hstop
      nstop_traj = nstop
    ENDIF


 
  !----------------------------------------------------------------------------
  !- Section 4.2: Check istart_mode_traj and derive nstart_traj_gp
  !----------------------------------------------------------------------------

    ! First initialization mode: the start times are provided by the user 
    ! -------------------------------------------------------------------
    IF (istart_mode_traj == 1) THEN
      IF (ANY(hstart_traj(:) /= hstart_traj_d(:)) .AND.                     &
          ALL(nstart_traj(:) == nstart_traj_d(:))) THEN
        DO i = 1, nmax_traj_starttime
          IF (hstart_traj(i) /= zrnoval) THEN
            nstart_traj(i) = NINT(hstart_traj(i) * 3600.0_wp / dt)
            IF (ABS(nstart_traj(i) - (hstart_traj(i) * 3600.0_wp / dt)) &
                            > zeps * (hstart_traj(i) * 3600.0_wp / dt)) THEN
              PRINT *, ' WARNING  *** TRAJECTORY module *** '
              PRINT *, '          *** hstart_traj(',i,') is not a multiple of dt.'
              PRINT *, '          *** nstart_traj(',i,') is thus based on a  rounded hstart_traj(',i,').'
            ENDIF
          ELSE
            nstart_traj(i) = iznoval
          ENDIF
        ENDDO
      ELSEIF (ALL(hstart_traj(:) == hstart_traj_d(:)) .AND.                 &
              ANY(nstart_traj(:) /= nstart_traj_d(:))) THEN
        DO i = 1, nmax_traj_starttime
          IF (nstart_traj(i) /= iznoval) THEN
            hstart_traj(i) = nstart_traj(i) * dt / 3600.0_wp
          ELSE
            hstart_traj(i) = zrnoval
          ENDIF
        ENDDO
      ELSEIF (ANY(hstart_traj(:) /= hstart_traj_d(:)) .AND.                 &
              ANY(nstart_traj(:) /= nstart_traj_d(:))) THEN
        PRINT *, ' WARNING  *** TRAJECTORY module *** '
        PRINT *, '          *** hstart_traj and nstart_traj are both specified.'
        PRINT *, '          *** hstart_traj will be ignored.'
        DO i = 1, nmax_traj_starttime
          IF (nstart_traj(i) /= iznoval) THEN
            hstart_traj(i) = nstart_traj(i) * dt / 3600.0_wp
          ELSE
            hstart_traj(i) = zrnoval
          ENDIF
        ENDDO
      ELSEIF (ALL(hstart_traj(:) == hstart_traj_d(:)) .AND.                   &
              ALL(nstart_traj(:) == nstart_traj_d(:))) THEN
        PRINT *, ' ERROR    *** TRAJECTORY module *** '
        PRINT *, '          *** all hstart_traj and nstart_traj  are default values'
        ierrstat = 10063_iintegers
        RETURN
      ENDIF      
        
      IF (ANY(hcomb_start_traj(:) /= zrnoval) .OR.                            &
          ANY(ncomb_start_traj(:) /= iznoval)) THEN
        PRINT *, ' WARNING  *** TRAJECTORY module *** '  
        PRINT *, '          *** hcomb/ncomb_start_traj must not be set for istart_mode_traj=1'
        PRINT *, '          *** hcomb/ncomb_start_traj will be set to default values'
        hcomb_start_traj(:) = hcomb_start_traj_d(:)
        ncomb_start_traj(:) = ncomb_start_traj_d(:)
      ENDIF    
        
      num_start_t = COUNT(nstart_traj(:) /= iznoval)
      IF (num_start_t > nmax_traj_starttime) THEN
        PRINT *, ' ERROR  *** TRAJECTORY module *** '
        PRINT *, '        *** the number of hstart_traj/nstart_traj is larger than',nmax_traj_starttime,'.'
        ierrstat = 10066_iintegers
        RETURN
      ENDIF

      DO i = 1, num_start_t
        nstart_traj_gp(i) = nstart_traj(i)
        IF (nstart == 0_iintegers) THEN
          IF (nstart_traj_gp(i) <  nstart .OR. nstart_traj_gp(i) > nstop_traj) THEN
            PRINT *, ' ERROR  *** TRAJECTORY module *** '
            PRINT *, '        *** hstart_traj/nstart_traj(',i,') is not within [hstart/nstart;hstop_/nstop_traj]'
            ierrstat = 10069_iintegers
            RETURN
          ENDIF
        ELSE
          IF (nstart_traj_gp(i) > nstop_traj) THEN
            PRINT *, ' ERROR  *** TRAJECTORY module *** '
            PRINT *, '        *** hstart_traj/nstart_traj(',i,') is > hstop_/nstop_traj'
            ierrstat = 10070_iintegers
            RETURN
          ENDIF
        ENDIF
      ENDDO
      DO i = 1, num_start_t-1
        IF (nstart_traj_gp(i) >= nstart_traj_gp(i+1)) THEN
          PRINT *, ' ERROR    *** TRAJECTORY module *** '
          PRINT *, '          hstart_traj/nstart_traj(',i,') is >= hstart_traj/nstart_traj(',i+1,')'
          ierrstat = 10072_iintegers
          RETURN
        ENDIF
        IF (REAL((nstart_traj_gp(i+1)- nstart_traj_gp(i))*dt,wp) < 60.0_wp ) THEN
          PRINT *, ' ERROR    *** TRAJECTORY module *** '
          PRINT *, '          trajectories',i,'and ',i+1,'are started appart in less than a minute.'
          ierrstat = 10075_iintegers
          RETURN
        ENDIF
      ENDDO
      IF (REAL(nstart_traj_gp(num_start_t)*dt/60._wp,wp) > 999999.0_wp) THEN
        PRINT *, ' ERROR    *** TRAJECTORY module *** '
        PRINT *, '          hstart_traj/nstart_traj(',num_start_t,') is larger than 999999 min.'
        ierrstat = 10078_iintegers
        RETURN
      ENDIF

      nstartfirst_traj = nstart_traj_gp(1)
      num_start_t_tot  = num_start_t

    ! Second initialization mode: the start triplet is given by the user 
    ! -------------------------------------------------------------------
    ELSEIF (istart_mode_traj == 2) THEN
      IF (ALL(hcomb_start_traj(:) /= hcomb_start_traj_d(:)) .AND.                     &
          ALL(ncomb_start_traj(:) == ncomb_start_traj_d(:))) THEN
        ncomb_start_traj(:) = NINT(hcomb_start_traj(:) * 3600.0_wp / dt)
        IF (ANY(ABS(ncomb_start_traj(:) - (hcomb_start_traj(:) * 3600.0_wp / dt)) &
                                 > zeps * (hcomb_start_traj(:) * 3600.0_wp / dt))) THEN
          PRINT *, ' WARNING  *** TRAJECTORY module *** '
          PRINT *, '          *** one/some element(s) of hcomb_start_traj is/are not a multiple of dt.'
          PRINT *, '          *** ncomb_start_traj elements are thus based on rounded hcomb_start_traj elements.'
        ENDIF
      ELSEIF (ALL(hcomb_start_traj(:) == hcomb_start_traj_d(:)) .AND.                 &
              ALL(ncomb_start_traj(:) /= ncomb_start_traj_d(:))) THEN
        hcomb_start_traj(:) = ncomb_start_traj(:) * dt / 3600.0_wp
      ELSEIF (ALL(hcomb_start_traj(:) /= hcomb_start_traj_d(:)) .AND.                 &
              ALL(ncomb_start_traj(:) /= ncomb_start_traj_d(:))) THEN
        PRINT *, ' WARNING  *** TRAJECTORY module *** '
        PRINT *, '          *** hcomb_start_traj and ncomb_start_traj are both specified.'
        PRINT *, '          *** hcomb_start_traj will be ignored.'
        hcomb_start_traj(:) = ncomb_start_traj(:) * dt / 3600.0_wp
      ELSEIF (ALL(hcomb_start_traj(:) == hcomb_start_traj_d(:)) .AND.                 &
              ALL(ncomb_start_traj(:) == ncomb_start_traj_d(:))) THEN
        PRINT *, ' ERROR    *** TRAJECTORY module *** '
        PRINT *, '          *** both hcomb_start_traj and ncomb_start_traj  are default values'
        ierrstat = 10081_iintegers
        RETURN
      ELSE
        PRINT *, ' ERROR    *** TRAJECTORY module *** '
        PRINT *, '          *** something went wrong with the specification of hcomb_start_traj/ncomb_start_traj'
        ierrstat = 10084_iintegers
        RETURN
      ENDIF      
      IF (nstart == 0_iintegers) THEN
        IF (ncomb_start_traj(1) < nstart .OR. ncomb_start_traj(1) > nstop_traj) THEN
          PRINT *, ' ERROR    *** TRAJECTORY module *** '
          PRINT *, '          *** hcomb_/ncomb_start_traj(1) must be within [hstart/nstart;hstop_/nstop_traj]'
          ierrstat = 10087_iintegers
          RETURN
        ENDIF
        IF (ncomb_start_traj(2) < nstart .OR. ncomb_start_traj(2) > nstop_traj) THEN
          PRINT *, ' ERROR    *** TRAJECTORY module *** '
          PRINT *, '          *** hcomb_/ncomb_start_traj(2) must be within [hstart/nstart;hstop_/nstop_traj]'
          ierrstat = 10090_iintegers
          RETURN
        ENDIF
      ELSE
        IF (ncomb_start_traj(1) > nstop_traj) THEN
          PRINT *, ' ERROR    *** TRAJECTORY module *** '
          PRINT *, '          *** hcomb_/ncomb_start_traj(1) is > hstop_/nstop_traj'
          ierrstat = 10088_iintegers
          RETURN
        ENDIF
        IF (ncomb_start_traj(2) > nstop_traj) THEN
          PRINT *, ' ERROR    *** TRAJECTORY module *** '
          PRINT *, '          *** hcomb_/ncomb_start_traj(2) is > hstop_/nstop_traj'
          ierrstat = 10091_iintegers
          RETURN
        ENDIF
      ENDIF
      IF (ncomb_start_traj(3) < 1 .OR. ncomb_start_traj(3) >                  &
                               (ncomb_start_traj(2)-ncomb_start_traj(1))) THEN
        PRINT *, ' ERROR    *** TRAJECTORY module *** '
        PRINT *, '          *** hcomb_/ncomb_start_traj(3) must be within [dt/1;(h/ncomb_start_traj(2)-h/ncomb_start_traj(1))]'
        ierrstat = 10093_iintegers
        RETURN
      ENDIF
   
      IF (ANY(hstart_traj(:) /= zrnoval) .OR. ANY(nstart_traj(:) /= iznoval)) THEN
        PRINT *, ' WARNING  *** TRAJECTORY module *** '  
        PRINT *, '          *** hstart_/nstart_traj must not be set for istart_mode_traj=2'
        PRINT *, '          *** hstart_/nstart_traj will be set to default values'
        hstart_traj(:) = hstart_traj_d(:)
        nstart_traj(:) = nstart_traj_d(:)
      ENDIF 


      num_start_t = (ncomb_start_traj(2)-ncomb_start_traj(1))/ncomb_start_traj(3) + 1_iintegers
      IF (num_start_t > nmax_traj_starttime) THEN
        PRINT *, ' ERROR  *** TRAJECTORY module *** '
        PRINT *, '        *** the number of start times is larger than',nmax_traj_starttime,'.'
        ierrstat = 10095_iintegers
        RETURN
      ENDIF

      DO i = 1, num_start_t
        nstart_traj_gp(i) = ncomb_start_traj(1) + (i-1)*ncomb_start_traj(3)
      ENDDO
      DO i = 1, num_start_t-1
        IF (nstart_traj_gp(i) >= nstart_traj_gp(i+1)) THEN
          PRINT *, ' ERROR    *** TRAJECTORY module *** '
          PRINT *, '          hcomb/ncomb_start_traj specified causes start of traj(',i,') to be larger than start of traj(',i+1,')'
          ierrstat = 10098_iintegers
          RETURN
        ENDIF
        IF (REAL((nstart_traj_gp(i+1)- nstart_traj_gp(i))*dt,wp) < 60.0_wp ) THEN
          PRINT *, ' ERROR    *** TRAJECTORY module *** '
          PRINT *, '          trajectories',i,'and ',i+1,'are started appart in less than a minute.'
          ierrstat = 10101_iintegers
          RETURN
        ENDIF
      ENDDO
      IF (REAL(nstart_traj_gp(num_start_t)*dt/60._wp,wp) > 999999.0_wp) THEN
        PRINT *, ' ERROR    *** TRAJECTORY module *** '
        PRINT *, '          hcomb/ncomb_start_traj specified causes start of traj(',num_start_t,') to be larger than 999999 min.'
        ierrstat = 10104_iintegers
        RETURN
      ENDIF

      nstartfirst_traj = nstart_traj_gp(1)
      num_start_t_tot  = num_start_t

    ! Third initialization mode: the starting dates are provided in the startfile 
    ! ----------------------------------------------------------------------------
    ELSEIF (istart_mode_traj == 3_iintegers) THEN
      IF (ANY(hstart_traj(:) /= zrnoval) .OR. ANY(nstart_traj(:) /= iznoval)) THEN
        PRINT *, ' WARNING  *** TRAJECTORY module *** '  
        PRINT *, '          *** hstart_/nstart_traj must not be set for istart_mode_traj=3'
        PRINT *, '          *** hstart_/nstart_traj will be set to default values'
        hstart_traj(:) = hstart_traj_d(:)
        nstart_traj(:) = nstart_traj_d(:)
      ENDIF 

      IF (ANY(hcomb_start_traj(:) /= zrnoval) .OR.                            &
          ANY(ncomb_start_traj(:) /= iznoval)) THEN
        PRINT *, ' WARNING  *** TRAJECTORY module *** '  
        PRINT *, '          *** hcomb/ncomb_start_traj must not be set for istart_mode_traj=3'
        PRINT *, '          *** hcomb/ncomb_start_traj will be set to default values'
        hcomb_start_traj(:) = hcomb_start_traj_d(:)
        ncomb_start_traj(:) = ncomb_start_traj_d(:)
      ENDIF    

    ! Any other mode doesn't exist
    ! -------------------------------------------------------------------
    ELSE
      PRINT *, ' ERROR    *** TRAJECTORY module *** '
      PRINT *, '          *** istart_mode_traj must be equal to 1, 2 or 3'
      ierrstat = 10107_iintegers
      RETURN
    ENDIF

  !----------------------------------------------------------------------------
  !- Section 4.3: Check the start file
  !----------------------------------------------------------------------------

    IF (TRIM(startfile_traj) == ' ') THEN
      PRINT *, ' ERROR    *** TRAJECTORY module *** '
      PRINT *, '          *** a startfile (startfile_traj) must be provided'
      ierrstat = 10110_iintegers
      RETURN    
    ENDIF
 
  !----------------------------------------------------------------------------
  !- Section 4.4: Check output information: ydir, hinc_out_traj, ninc_out_traj
  !----------------------------------------------------------------------------

    IF (TRIM(ydir_traj) == ' ') THEN
      PRINT *, ' ERROR  *** TRAJECTORY module *** '
      PRINT *, '        *** an output directory (ydir_traj) must be provided'
      ierrstat = 10113_iintegers
      RETURN    
    ENDIF

    IF (hinc_out_traj /= hinc_out_traj_d .AND. ninc_out_traj == ninc_out_traj_d) THEN
      ninc_out_traj = NINT(hinc_out_traj * 3600.0_wp / dt)
      IF (ABS(ninc_out_traj-((hinc_out_traj * 60._wp) / dt)) > zeps *(hinc_out_traj * 60._wp) / dt) THEN
        PRINT *, ' WARNING  *** TRAJECTORY module *** '
        PRINT *, '          *** hinc_out_traj is not a multiple of dt.'
        PRINT *, '          *** ninc_out_traj will be based on a rounded hinc_out_traj.'
      ENDIF 
    ELSEIF (hinc_out_traj == hinc_out_traj_d .AND. ninc_out_traj /= ninc_out_traj_d) THEN
      hinc_out_traj = ninc_out_traj * dt / 3600.0_wp
    ELSEIF (hinc_out_traj /= hinc_out_traj_d .AND. ninc_out_traj /= ninc_out_traj_d) THEN
      PRINT *, ' WARNING  *** TRAJECTORY module *** '
      PRINT *, '          *** hinc_out_traj and ninc_out_traj are both specified.'
      PRINT *, '          *** hinc_out_traj will be ignored.'
      hinc_out_traj = ninc_out_traj * dt / 3600.0_wp
    ELSE
      PRINT *, ' ERROR    *** TRAJECTORY module *** '
      PRINT *, '          *** neither hinc_out_traj nor ninc_out_traj are specified.'
      ierrstat = 10116_iintegers
      RETURN 
    ENDIF      

    IF (ninc_out_traj < 1 .OR. ninc_out_traj > nstop_traj) THEN
      PRINT *, ' ERROR    *** TRAJECTORY module *** '
      PRINT *, '          *** h/ninc_out_traj is not within [dt/1;h/nstop_traj].'
      ierrstat = 10119_iintegers
      RETURN    
    ENDIF

  !----------------------------------------------------------------------------
  !- Section 4.5: Check the info about the traced variables and derive new info
  !----------------------------------------------------------------------------

    ntracevar_traj = COUNT(tracevar_traj(:) /= 'novar')

  ENDIF
!------------------------------------------------------------------------------
!- Section 5: Distribute variables to all nodes
!------------------------------------------------------------------------------


  IF (nproc > 1_iintegers) THEN
    
    IF (my_cart_id == 0_iintegers) THEN

      ! integer variables
      intbuf ( 1)  = ntracevar_traj
      intbuf ( 2)  = istart_mode_traj
      intbuf ( 3)  = nstop_traj
      intbuf ( 4)  = ninc_out_traj
      intbuf ( 5)  = ntracevar_traj
      intbuf ( 6)  = nstartfirst_traj 
      intbuf (7:6+nmax_traj_starttime) = nstart_traj_gp(:)
      intbuf (7+nmax_traj_starttime)   = num_start_t
      intbuf (8+nmax_traj_starttime:10+nmax_traj_starttime) = ncomb_start_traj(:)

      ! real variables
      realbuf(1)      = hinc_out_traj
      realbuf(2)      = hstop_traj
      realbuf(3:2+nmax_traj_starttime)  = hstart_traj(:)
      realbuf(2+nmax_traj_starttime+1:2+nmax_traj_starttime+3)= hcomb_start_traj(:)
      
      charbuf (1)                 = ydir_traj
      charbuf (2)                 = startfile_traj
      DO i=1,ntracevar_traj
        charbuf (2+i)(1:15)       = tracevar_traj(i)(1:15)
      ENDDO
    ENDIF

    CALL distribute_values (intbuf, 10+nmax_traj_starttime  , 0, imp_integers ,icomm_world,iz_err)
    CALL distribute_values (realbuf,2+nmax_traj_starttime+3, 0, imp_reals    ,icomm_world,iz_err)
    CALL distribute_values (charbuf,2+nmax_traced_var,       0, imp_character,icomm_world,iz_err)

    IF (my_cart_id /= 0_iintegers) THEN
      
      ntracevar_traj        = intbuf ( 1)
      istart_mode_traj      = intbuf ( 2)
      nstop_traj            = intbuf ( 3)
      ninc_out_traj         = intbuf ( 4)
      ntracevar_traj        = intbuf ( 5)
      nstartfirst_traj      = intbuf ( 6)
      nstart_traj_gp(:)     = intbuf ( 7:6+nmax_traj_starttime)
      num_start_t           = intbuf ( 7+nmax_traj_starttime)
      ncomb_start_traj(:)   = intbuf ( 8+nmax_traj_starttime:10+nmax_traj_starttime)

      hinc_out_traj         = realbuf(1)
      hstop_traj            = realbuf(2)
      hstart_traj  (:)      = realbuf(3:2+nmax_traj_starttime)
      hcomb_start_traj(:)   = realbuf(2+nmax_traj_starttime+1:2+nmax_traj_starttime+3)

      ydir_traj             = charbuf(1)
      startfile_traj        = charbuf(2)

      DO i=1,ntracevar_traj
        tracevar_traj(i)(1:15)    = charbuf(2+i)(1:15)
      ENDDO
    ENDIF 
  ENDIF

  IF (idebug > 0_iintegers .AND. my_cart_id == 0_iintegers) THEN
    PRINT *, '  *** TRAJECTORY module: ', ntracevar_traj, ' variables will be traced *** '
  ENDIF

!------------------------------------------------------------------------------
!- Section 6: Output of the namelist variables and their default values
!------------------------------------------------------------------------------

  IF (my_cart_id == 0_iintegers) THEN
    WRITE (nuspecif, '(A2)')  '  '
    WRITE (nuspecif, '(A24)') '0     NAMELIST:  trajctl'
    WRITE (nuspecif, '(A24)') '      -----------------'
    WRITE (nuspecif, '(A2)')  '  '
    WRITE (nuspecif, '(T7,A,T33,A,T51,A,T70,A)') 'Variable', 'Actual Value',  &
                                                 'Default Value', 'Format'

    WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                     &
                   'hstop_traj'     , hstop_traj     , hstop_traj_d     , ' R '
    WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                     &
                   'nstop_traj'     , nstop_traj     , nstop_traj_d     , ' I '

    WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                     &
                   'istart_mode_traj', istart_mode_traj, istart_mode_traj_d, ' I '

    WRITE (nuspecif, '(T8,A,T71,A3)') 'Start times of the trajectories:', ' R '
    WRITE (nuspecif, '(T10,A)') '  trajectory ID          start hour    start time step'  
    IF (istart_mode_traj == 1_iintegers) THEN  
      DO i = 1, nmax_traj_starttime
        IF (hstart_traj(i) /= zrnoval) THEN
          WRITE (nuspecif, '(T13,I12,A,F12.4, I12)') i, ':   ', hstart_traj(i),nstart_traj_gp(i)
        ENDIF
      ENDDO
    ELSEIF (istart_mode_traj == 2_iintegers) THEN  
      WRITE (nuspecif, '(T8,A,T71,A3)') 'hcomb_start_traj triplets:',       ' R '
      WRITE (nuspecif, '(T13,A,T33,F12.4,T52,F12.4)') 'begin ', hcomb_start_traj(1), &
                                                                hcomb_start_traj_d(1)
      WRITE (nuspecif, '(T13,A,T33,F12.4,T52,F12.4)') 'end   ', hcomb_start_traj(2), &
                                                                hcomb_start_traj_d(2)
      WRITE (nuspecif, '(T13,A,T33,F12.4,T52,F12.4)') 'incr  ', hcomb_start_traj(3), &
                                                                hcomb_start_traj_d(3)
      WRITE (nuspecif, '(T8,A,T71,A3)') 'ncomb_start_traj triplets:',       ' I '
      WRITE (nuspecif, '(T13,A,T33,I12,T52,I12.4)') 'begin ', ncomb_start_traj(1), &
                                                                ncomb_start_traj_d(1)
      WRITE (nuspecif, '(T13,A,T33,I12,T52,I12.4)') 'end   ', ncomb_start_traj(2), &
                                                                ncomb_start_traj_d(2)
      WRITE (nuspecif, '(T13,A,T33,I12,T52,I12.4)') 'incr  ', ncomb_start_traj(3), &
                                                                ncomb_start_traj_d(3)
    ENDIF
    WRITE (nuspecif, '(T8,A,T45,  A               )')                         &
                   'startfile_traj (C*300)   ',TRIM(startfile_traj)
    WRITE (nuspecif, '(T8,A,T45,  A               )')                         &
                   'ydir_traj (C*300)   ',TRIM(ydir_traj)
    WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                     & 
                   'hinc_out_traj ' , hinc_out_traj ,  hinc_out_traj_d  , ' R '
    WRITE (nuspecif, '(T8,A,T33,I12,T52,I12,T71,A3)')                     & 
                   'ninc_out_traj ' , ninc_out_traj ,  ninc_out_traj_d  , ' R '
    IF(ntracevar_traj >= 0) THEN
      WRITE(nuspecif, '(A)')  '       Traced variables along trajectories'
      IF (ntracevar_traj == 0_iintegers) THEN
        WRITE(nuspecif, '(A)')  '           none '  
      ELSE
        DO i=1,ntracevar_traj
          WRITE (nuspecif, '(T8,A,I3,A,T30,A10,T72,A)') 'tracevar_traj(',i,')', &
                                                   TRIM(tracevar_traj(i)),'C15'
        ENDDO
      ENDIF
    ENDIF
    WRITE (nuspecif, '(A2)') '  '
  ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE input_trajctl

!==============================================================================
!+ Module procedure for the read of the startfile
!------------------------------------------------------------------------------

SUBROUTINE read_startf_traj (idebug, nuin, ierrstat, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   A startfile containing the start positions of the trajectories has to be
!   provided. It is specified by the user via the namelist parameter 
!   startfile_traj. 
!   It is a text file and its format is strictly defined.
!   It contains:
!     - 3 header lines (no space at the beginning of the lines):
!        * 1st line looks like:"Reference Date 19870725_0000"
!        * 2nd line looks like:"time lon lat z" or "lon lat z"
!        * 3rd line looks like:"-------"
!     - several content lines that look like
!       " 1.00   -8.000  -2.500 9000.000" or
!       "-8.000  -2.500 9000.000"
!
!  For the content lines, the format is thus the following one:
!   - 1 empty position, 4 positions for the time, 1 empty position,
!     8 positions for the longitudes, 1 empty position, 7 positions for the
!     latitudes, 1 empty position, 8 positions for the height if time
!     is present in the file
!   - 1 empty position, 8 positions for the longitudes, 1 empty position, 
!     7 positions for the latitudes, 1 empty position, 8 positions for the 
!     height if time is not present
!
! The time is given in minutes and can have max. 2 decimal positions and
! can go up to 99 minutes only.
! The longitudes are rotated longitudes and can go from -180 to
! 180 degrees. They must be expressed in the same referential
! (pollon, pollat) than the COSMO grid. They can have max. 3 decimal
! positions.
! The latitudes are rotated latitudes and can go from -90 to 90
! degrees. They must be expressed in the same referential
! (pollon, pollat) than the COSMO grid. They can have max. 3 decimal
! positions.
! The height is above see level in meters. It can go up to
! 99999.999 meters. It can have up to 3 decimal positions.
!
! Make sure that the startfile respects all these format constraints.
! 
!  
! Method:
!   The 3 header lines are read and the 2nd is checked in order
!   to determine if time is present in the startfile or not.  
!   In an endless loop the rest of the startfile is read. 
!   Once there are no more records to read, the loop is stopped (exit).
!   The values read by PE0 are stored in traj_startinfo and communicated
!   to all PEs.
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
  idebug             ! debug level

INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
  nuin               ! unit number for startfile file

INTEGER   (KIND=iintegers),   INTENT (OUT)     ::        &
  ierrstat           ! error status variable

CHARACTER(LEN=*),             INTENT (OUT)     ::        &
  yerrmsg            ! error message

 
! Other local variables
INTEGER (KIND=iintegers)                       ::        &
  izerr, izstat,   & ! local error codes
  izcount,         & ! number of lines read from the file
  izout,           & ! number of start positions out of the domain
  i                  ! loop index

CHARACTER (LEN=100)                            ::        &
  yzdummy            ! dummy to check if "time" is present in startfile

CHARACTER (LEN=4)                              ::        &
  yzdummy2           ! dummy to check if "time" is present in startfile

REAL(KIND=wp),     ALLOCATABLE                 ::        &
  zrealbuf(:)        ! buffer to exchange the traj_startinfo values

INTEGER(KIND=iintegers), ALLOCATABLE           ::        &
  zintbuf(:)         ! buffer to exchange the nstart_traj_iv values and others

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE read_startf_traj
!------------------------------------------------------------------------------

  yerrmsg  = ''
  ierrstat = 0_iintegers
  izerr    = 0_iintegers
  izstat   = 0_iintegers
  izcount  = 0_iintegers
  izout    = 0_iintegers

!------------------------------------------------------------------------------
!- Section 1: Allocate structures to hold the 3 spatial coord. and time
!------------------------------------------------------------------------------


  ALLOCATE(traj_startinfo(nmax_traj_startpt,3),STAT=izstat)

  IF (izstat /= 0_iintegers) THEN
    ierrstat = 10130_iintegers
    yerrmsg  = 'Allocation ot traj_startinfo failed.'
    RETURN
  ENDIF     

  traj_startinfo(:,:) = -999.0_wp

  ALLOCATE(traj_startinfo_time(nmax_traj_startpt),STAT=izstat)

  IF (izstat /= 0_iintegers) THEN
    ierrstat = 10133_iintegers
    yerrmsg  = 'Allocation ot traj_startinfo_time failed.'
    RETURN
  ENDIF     

  traj_startinfo_time(:) = -999.0_wp

!-----------------------------------------------------------------------------
!- Section 2: Read the 3 spatial coord. and the start time (istart_mode_traj =3)
!-----------------------------------------------------------------------------
     
  IF ( my_cart_id == 0 ) THEN

    DO  ! endless loop until EOF

      izcount = izcount + 1_iintegers

      IF (izcount > nmax_traj_startpt+3_iintegers) THEN
        PRINT *, ' ERROR    *** TRAJECTORY module *** '
        PRINT *, '          *** the startfile contains more than ',nmax_traj_startpt,' lines'
        ierrstat = 10136_iintegers
        RETURN   
      ENDIF
         

      ! skip 1st and 3rd lines
      IF (ANY(izcount == (/1, 3/))) THEN
        READ(nuin, *, IOSTAT=izerr) 
        IF (izerr > 0_iintegers) THEN
          yerrmsg  = 'Cannot read 1st or 3rd header line of the startfile.'
          ierrstat = 10139_iintegers
          RETURN
        ELSEIF (izerr < 0_iintegers) THEN
          izcount = izcount - 1_iintegers
          EXIT ! EOF reached
        ENDIF
        
      ! test 2nd line for time
      ELSEIF (izcount == 2_iintegers) THEN
        READ(nuin,'(A)',IOSTAT=izerr) yzdummy
        IF (izerr > 0_iintegers) THEN
          yerrmsg  = 'Cannot read 2nd header line of the startfile.'
          ierrstat = 10142_iintegers
          RETURN
        ELSEIF (izerr < 0_iintegers) THEN
          izcount = izcount - 1_iintegers
          EXIT ! EOF reached
        ENDIF
        yzdummy2 = TRIM(ADJUSTL(yzdummy))

        IF (istart_mode_traj == 3 .AND. .NOT.(yzdummy2 == 'time')) THEN
          PRINT *, ' ERROR    *** TRAJECTORY module *** '
          PRINT *, '          *** the startfile (startfile_traj) must contain'
          PRINT *, '          *** a column with start times'
          ierrstat = 10145_iintegers
          RETURN    
        ENDIF

      ! read from the 4th line on
      ELSEIF (izcount > 3_iintegers) THEN
        IF (yzdummy2 == 'time') THEN
          READ (nuin,fmt=*,IOSTAT=izerr) traj_startinfo_time(izcount-3), traj_startinfo(izcount-3,1), &
                traj_startinfo(izcount-3,2), traj_startinfo(izcount-3,3)
          IF (izerr > 0_iintegers) THEN
            yerrmsg  = 'Cannot read content line of the startfile.'
            ierrstat = 10148_iintegers
            RETURN
          ELSEIF (izerr < 0_iintegers) THEN
            izcount = izcount - 4_iintegers
            EXIT ! EOF reached
          ENDIF
        ELSE
          READ(nuin,fmt=*,IOSTAT=izerr) traj_startinfo(izcount-3,1), &
               traj_startinfo(izcount-3,2), traj_startinfo(izcount-3,3)
          IF (izerr > 0_iintegers) THEN
            yerrmsg  = 'Cannot read content line of the startfile.'
            ierrstat = 10151_iintegers
            RETURN
          ELSEIF (izerr < 0_iintegers) THEN
            izcount = izcount - 4_iintegers
            EXIT ! EOF reached
          ENDIF
        ENDIF
      ENDIF      
    ENDDO  ! end of endless loop

    ! The file is empty (except header)
    IF (izcount < 1_iintegers) THEN
      yerrmsg  = 'Startfile is empty (up to the header).'
      ierrstat = 10154_iintegers
      RETURN
    ! The file is too long
    ELSEIF (izcount > nmax_traj_startpt) THEN
      yerrmsg  = 'Startfile contains too many start points.'
      ierrstat = 10157_iintegers
      RETURN
    ENDIF
  ENDIF
  
!-----------------------------------------------------------------------------
!- Section 3: Determine the number of start positions on all PEs
!-----------------------------------------------------------------------------

  IF (nproc > 1_iintegers) THEN
    ! distribute error status to all processors
    CALL distribute_values  (izstat, 1, 0, imp_integers,  icomm_world, izerr)
  ENDIF
  IF (izstat /= 0_iintegers) THEN
    yerrmsg  = 'Distributing the values of izstat failed.'
    ierrstat = 10160_iintegers
    RETURN
  ENDIF

  IF (nproc > 1_iintegers) THEN
    ! distribute izcount to all processors
    CALL distribute_values  (izcount, 1, 0, imp_integers,  icomm_world, izerr)
  ENDIF


  num_start_pt_tot = izcount
  num_start_pt(:)  = izcount


  ALLOCATE(zrealbuf(3*num_start_pt_tot), STAT=izerr)
  IF (izerr /= 0_iintegers) THEN
    yerrmsg  = 'Allocation of zrealbuf failed.'
    ierrstat = 10163_iintegers
    RETURN
  ENDIF 

  ALLOCATE(zintbuf(4+2*nmax_traj_starttime+nmax_traj_startpt), STAT=izerr)
  IF (izerr /= 0_iintegers) THEN
    yerrmsg  = 'Allocation of zintbuf failed.'
    ierrstat = 10164_iintegers
    RETURN
  ENDIF 

!------------------------------------------------------------------------------
!- Section 4: Perform basic checks
!------------------------------------------------------------------------------

  IF ( my_cart_id == 0 ) THEN

    ! Check that the start points are located in the domain
    DO i=1,num_start_pt_tot
      IF ((traj_startinfo(i,1) < startlon_tot + dlon*nboundlines)         .OR.&
          (traj_startinfo(i,1) > startlon_tot + dlon*(ie_tot-nboundlines)).OR.&
          (traj_startinfo(i,2) < startlat_tot + dlat*nboundlines)         .OR.&
          (traj_startinfo(i,2) > startlat_tot + dlat*(je_tot-nboundlines))) THEN
        izout = izout + 1_iintegers
        PRINT *, ' WARNING  *** TRAJECTORY module *** '
        PRINT *, '          *** Trajectory with start coordinates (',         &
                 traj_startinfo(i,1),',',traj_startinfo(i,2),                 &
                 ') is out of the domain and will not be computed.'
      ENDIF      
    ENDDO

    IF ( izout == num_start_pt_tot ) THEN
      yerrmsg  = 'All trajectory startpoints are out of the domain.'
      ierrstat = 10166_iintegers
      RETURN
    ENDIF

    ! Debug information
    IF (idebug > 0_iintegers) THEN
      PRINT *, '  *** TRAJECTORY module: ',num_start_pt_tot, ' startpoints read from startfile *** '
      PRINT *, '  *** TRAJECTORY module: min/max values read from startfile *** '
      PRINT *, '  ***            (lon/lat are in rotated syst. in degrees, z is in meters AMSL)'
      PRINT *, '                 min(lon)/max(lon):',MINVAL(traj_startinfo(1:num_start_pt_tot,1)),&
                                                     MAXVAL(traj_startinfo(1:num_start_pt_tot,1))
      PRINT *, '                 min(lat)/max(lat):',MINVAL(traj_startinfo(1:num_start_pt_tot,2)),&
                                                     MAXVAL(traj_startinfo(1:num_start_pt_tot,2))
      PRINT *, '                 min(z)/max(z):',    MINVAL(traj_startinfo(1:num_start_pt_tot,3)),&
                                                     MAXVAL(traj_startinfo(1:num_start_pt_tot,3))
      IF (istart_mode_traj == 3) THEN
        PRINT *, '  ***            (t is in hours after model start)'
        PRINT *, '                 min(t)/max(t):',  MINVAL(traj_startinfo_time(1:num_start_pt_tot)),&
                                                     MAXVAL(traj_startinfo_time(1:num_start_pt_tot))
      ENDIF
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
!- Section 5: Determine start times of trajectories (istart_mode_traj = 3)
!------------------------------------------------------------------------------ 
  
  IF (istart_mode_traj == 3_iintegers .AND. my_cart_id == 0_iintegers) THEN

    ! Convert start times in model time steps
    nstart_traj_iv(:) = NINT(traj_startinfo_time(:)*3600_wp/dt)

    ! Determine number of start times
    num_start_t       = 1_iintegers
    num_start_t_tot   = 1_iintegers
    num_start_pt(:)   = 0_iintegers
    nstart_traj_gp(1) = nstart_traj_iv(1)

    DO i = 1,num_start_pt_tot
      IF (nstart == 0_iintegers) THEN
        IF (nstart_traj_iv(i) < nstart .OR. nstart_traj_iv(i) > nstop_traj) THEN
          PRINT *, ' ERROR  *** TRAJECTORY module *** '
          PRINT *, '        *** start time(',i,') is not within [hstart/nstart;hstop_/nstop_traj]'
          ierrstat = 10169_iintegers
          RETURN
        ENDIF
      ELSE
        IF (nstart_traj_iv(i) > nstop_traj) THEN
          PRINT *, ' ERROR  *** TRAJECTORY module *** '
          PRINT *, '        *** start time(',i,') is > hstop_/nstop_traj'
          ierrstat = 10170_iintegers
          RETURN
        ENDIF
      ENDIF
      IF (nstart_traj_iv(i) * dt / 60._wp  > 999999.0_wp) THEN
        PRINT *, ' ERROR    *** TRAJECTORY module *** '
        PRINT *, '          start time(',i,') is larger than 999999 min.'
        ierrstat = 10172_iintegers
        RETURN
      ENDIF
    ENDDO

    DO i = 2,num_start_pt_tot
      IF (nstart_traj_iv(i) < nstart_traj_iv(i-1)) THEN
        PRINT *, ' ERROR  *** TRAJECTORY module *** '
        PRINT *, '        *** start time (',i,') is < start time (',i-1,')'
        ierrstat = 10175_iintegers
        RETURN
      ENDIF

      ! Check if the starttime is the same as the previous one and if not...
      IF (.NOT.(ANY(nstart_traj_gp(1:num_start_t_tot) == nstart_traj_iv(i)))) THEN
        IF (num_start_t_tot <= nmax_traj_starttime) THEN
          ! determine the number of start points belonging 
          ! to the previous starttime ("current group")
          num_start_pt(num_start_t_tot) = (i-1) -                             &
                   SUM(num_start_pt(1:MAX(1,num_start_t_tot - 1)))
          ! increase the group index
          num_start_t_tot = num_start_t_tot + 1_iintegers
          ! define the new starttime for the next group
          nstart_traj_gp(num_start_t_tot) = nstart_traj_iv(i)
        ELSE
          PRINT *,'ERROR  *** TRAJECTORY module *** '
          PRINT *,'       number of starting times is larger than ',nmax_traj_starttime,'.'
          ierrstat = 10178_iintegers 
          RETURN
        ENDIF
      ENDIF
    ENDDO

    num_start_pt(num_start_t_tot) = num_start_pt_tot-                         &
                                SUM(num_start_pt(1:num_start_t_tot-1_iintegers))

    nstartfirst_traj = MINVAL(nstart_traj_gp(1:num_start_t_tot))

  ENDIF

!------------------------------------------------------------------------------
!- Section 6: Distribute variables to all nodes
!------------------------------------------------------------------------------

  IF (nproc > 1_iintegers) THEN 
    IF (istart_mode_traj == 3_iintegers) THEN
      IF (my_cart_id == 0_iintegers) THEN
        zintbuf ( 1) = nstartfirst_traj
        zintbuf ( 2) = num_start_t
        zintbuf ( 3) = num_start_t_tot
        zintbuf ( 4) = num_start_pt_tot
        zintbuf (5:4+nmax_traj_starttime) = nstart_traj_gp(:)
        zintbuf (5+nmax_traj_starttime:4+2*nmax_traj_starttime) = num_start_pt(:)
        zintbuf (5+2*nmax_traj_starttime:                                     &
                4+2*nmax_traj_starttime+nmax_traj_startpt) = nstart_traj_iv(:)
      ENDIF
       
      CALL distribute_values(zintbuf,4+2*nmax_traj_starttime+nmax_traj_startpt,&
                             0, imp_integers, icomm_world, izerr)

      IF (my_cart_id /= 0_iintegers) THEN
        nstartfirst_traj   = zintbuf( 1)
        num_start_t        = zintbuf( 2)
        num_start_t_tot    = zintbuf( 3)
        num_start_pt_tot   = zintbuf( 4)
        nstart_traj_gp(:)  = zintbuf(5:4+nmax_traj_starttime)
        num_start_pt(:)    = zintbuf(5+nmax_traj_starttime:4+2*nmax_traj_starttime)
        nstart_traj_iv(:)  = zintbuf (5+2*nmax_traj_starttime:                &
                                     4+2*nmax_traj_starttime+nmax_traj_startpt)
      ENDIF

    ENDIF
    
    IF (my_cart_id == 0_iintegers) THEN
      zrealbuf ( 1                    :  num_start_pt_tot)  = traj_startinfo(1:num_start_pt_tot,1)
      zrealbuf ( 1+  num_start_pt_tot :2*num_start_pt_tot)  = traj_startinfo(1:num_start_pt_tot,2)
      zrealbuf ( 1+2*num_start_pt_tot :3*num_start_pt_tot)  = traj_startinfo(1:num_start_pt_tot,3)
    ENDIF

    CALL distribute_values(zrealbuf,3*num_start_pt_tot,0,imp_reals,icomm_world,izerr)

    IF (my_cart_id /= 0_iintegers) THEN
      traj_startinfo(1:num_start_pt_tot,1)    = zrealbuf ( 1                    :  num_start_pt_tot)
      traj_startinfo(1:num_start_pt_tot,2)    = zrealbuf ( 1+  num_start_pt_tot :2*num_start_pt_tot)
      traj_startinfo(1:num_start_pt_tot,3)    = zrealbuf ( 1+2*num_start_pt_tot :3*num_start_pt_tot)
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
!- Section 7: Finalization
!------------------------------------------------------------------------------

  DEALLOCATE(zrealbuf, STAT=izerr)
  IF (izerr /= 0_iintegers) THEN
    yerrmsg  = 'Deallocation of zrealbuf failed.'
    ierrstat = 10190_iintegers
    RETURN
  ENDIF  

  DEALLOCATE(zintbuf, STAT=izerr)
  IF (izerr /= 0_iintegers) THEN
    yerrmsg  = 'Deallocation of zintbuf failed.'
    ierrstat = 10191_iintegers
    RETURN
  ENDIF  

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE read_startf_traj

!==============================================================================
!+ Module procedure for the output initialization
!------------------------------------------------------------------------------

SUBROUTINE output_traj_init(izdebug, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   Output is initialized: the required variables are allocated, the needed
!   number of output files is computed, the files are created, the dimensions 
!   and global definitions are written out and the variables defined.
!   At least one output file per start time is written out. If there is
!   not enough space in one file to contain all trajectories that start
!   at this start time, additional files are created.
!   The name of the files is constructed in this way:
!     traj_tXXXXXX_pYYY.nc where:
!        XXXXXX: start time in minutes
!           YYY: number of the file for this start time
!
!   Only the NetCDF format is supported.
! 
! Method:
!   The trajectories are split over several NetCDF files, each of them having
!   a maximal size of 2 GB. The trajectories are split first according to
!   their start time and then further split according to their ID if needed. 
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
INTEGER   (KIND=iintegers), INTENT(IN)            ::                          &
  izdebug      ! debug level

INTEGER   (KIND=iintegers), INTENT(OUT)           ::                          &
  ierror       ! error status

CHARACTER (LEN= *),         INTENT(OUT)           ::                          &
  yerrmsg      ! error message
 
! Other local variables
INTEGER (KIND=iintegers)    :: &             
  n, k, i, j,                  & ! loop indices
  istat, izerr,                & ! local error codes
  nzfiles,                     & ! number of output files per start time
  nztraj,                      & ! number of traj. start points per output file
  nzout,                       & ! number of trajectory output steps 
                                 ! (for each start time)
  iztime_tmp                     ! tmp variable to hold time information (integer)

REAL(KIND=wp)               :: &
  rztime_tmp                     ! tmp variable to hold time information (real)  

CHARACTER (LEN=30)          :: &
  yzout_filename                 ! filename for NetCDF output

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE output_traj_init
!------------------------------------------------------------------------------

  nzout          = 0_iintegers
  j              = 0_iintegers
  izerr          = 0_iintegers
  istat          = 0_iintegers
  ierror         = 0_iintegers
  yerrmsg        = ' '
  yzout_filename = '  '



!------------------------------------------------------------------------------
!- Section 1: Precomputations
!------------------------------------------------------------------------------


#ifdef NETCDF
  IF (my_cart_id == 0_iintegers) THEN

    ! Allocate NetCDF variables

    ALLOCATE (outstep  (             nmax_out_per_start_t*num_start_t_tot)  ,STAT=istat)
    izerr=izerr+istat
    ALLOCATE (varId(ntracevar_traj+3,nmax_out_per_start_t*num_start_t_tot)  ,STAT=istat)
    izerr=izerr+istat
    ALLOCATE (ncID     (             nmax_out_per_start_t*num_start_t_tot)  ,STAT=istat)
    izerr=izerr+istat
    ALLOCATE (IdDimId  (             nmax_out_per_start_t*num_start_t_tot)  ,STAT=istat)
    izerr=izerr+istat
    ALLOCATE (varLonId (             nmax_out_per_start_t*num_start_t_tot)  ,STAT=istat)
    izerr=izerr+istat
    ALLOCATE (varLatId (             nmax_out_per_start_t*num_start_t_tot)  ,STAT=istat)
    izerr=izerr+istat
    ALLOCATE (varZId   (             nmax_out_per_start_t*num_start_t_tot)  ,STAT=istat)
    izerr=izerr+istat
    ALLOCATE (varTimeId(             nmax_out_per_start_t*num_start_t_tot)  ,STAT=istat)
    izerr=izerr+istat
    ALLOCATE (id_file  (             nmax_out_per_start_t*num_start_t_tot,2),STAT=istat)
    izerr=izerr+istat

    IF (izerr /= 0_iintegers) THEN
      ierror  = 10200_iintegers
      yerrmsg = 'Allocation of NetCDF variables for trajectories failed.'
      RETURN
    ENDIF

    outstep(:) = 0_iintegers

    DO n = 1, num_start_t_tot
      ! number of output steps
      nzout = (nstop_traj-nstart_traj_gp(n))/ninc_out_traj + 1_iintegers
 
      ! number of output files required per starting time
      ! it is [the number of output steps * number of starting positions * 
      ! (number of traced variables + lon/lat/z/time) * number of byes]/
      ! [2 GB -2MB for header]
      ! IMPORTANT NOTE: all terms that appear in the computation below
      !                 are of kind iintegers, which corresponds to the 
      !                 default representation of an integer. On most
      !                 systems, it is 4 bytes. Values ranging from
      !                 -2**31 to (2**31-1) can thus be represented.
      !                 However, multiplying the number of start points
      !                 by the number of output steps can produce
      !                 larger numbers. It is thus important to convert
      !                 all terms in a 8-bytes representation (since the
      !                 operation "format" is performed according to the
      !                 "format of the variables used in the computation). 

      nzfiles = MAX(1_iintegers, INT(INT(nzout,8) * INT(MAXVAL(num_start_pt(1:num_start_t_tot)),8)       &
               * INT((4_iintegers + ntracevar_traj),8)*INT(wp,8)          &
               /INT((max_nc_size - header_nc_size),8),4)+1_iintegers) 
      IF (nzfiles > nmax_out_per_start_t) THEN
        ierror  = 10203_iintegers    
        yerrmsg = 'Only up to 100 files are allowed per initialization time for the trajectory module.'
        RETURN  
      ENDIF

      ! if restart, set output position
      IF (nstart > 0_iintegers) THEN
        IF (ntstep >= nstart_traj_gp(n)) THEN
          outstep(n) = (ntstep-nstart_traj_gp(n))/ninc_out_traj
        ENDIF
      ENDIF   
       
      ! number of trajectory (starting point) per file
      nztraj = num_start_pt(n)/nzfiles     

      IF (izdebug > 0_iintegers) THEN
        PRINT*,' TRAJECTORY output '
        PRINT*,'   number of output steps:',nzout, ' number of output files: ',nzfiles, &
               '   number of traced variables: ',ntracevar_traj,                        &
               '   number of starting points: ' ,num_start_pt(n),                             &
               '   number of starting points per file: ', nztraj
      ENDIF

      DO k = 1, nzfiles
        ! Compute the index of the 1st and last trajectory in each file
        IF ( k == nzfiles ) THEN
          id_file(k+j,1)  = (k-1)*nztraj+1_iintegers+SUM(num_start_pt(1:n-1))
          id_file(k+j,2)  = num_start_pt(n)+SUM(num_start_pt(1:n-1))
          nztraj          = num_start_pt(n)-(k-1)*nztraj
        ELSE
          id_file(k+j,1) = (k-1)* nztraj + SUM(num_start_pt(1:n-1)) + 1_iintegers
          id_file(k+j,2) =    k * nztraj + SUM(num_start_pt(1:n-1))
        ENDIF

        !----------------------------------------------------------------------
        !- Section 2: Create output files in case of conventionnal model start
        !----------------------------------------------------------------------
        
        ! Construct NetCDF filenames
        WRITE (yzout_filename,"(A6,I6.6,A2,I3.3,A3)") 'traj_t',               &
               INT(nstart_traj_gp(n)*dt/60._wp,KIND=iintegers) ,'_p',k,'.nc'

        IF ( nstart == 0 ) THEN        
          
          IF (izdebug > 0_iintegers) THEN
            PRINT *,'Trajectory output written to ',TRIM(TRIM(ydir_traj)//'/'// &
                                                    TRIM(yzout_filename))
          ENDIF
          
          !--------------------------------------------------------------------
          !- Section 2.1: Create file and define global attributes and dimensions
          !--------------------------------------------------------------------
          
          ! Create NetCDF files
          istat = nf90_create(TRIM(TRIM(ydir_traj)//'/'//TRIM(yzout_filename)),&
                  IOR(nf90_CLOBBER,nf90_64bit_offset),ncID(k+j))
          IF (istat /= nf90_noerr) THEN
            ierror  = 10206_iintegers 
            yerrmsg = TRIM(nf90_strerror(istat))
            RETURN
          ENDIF
          istat = nf90_set_fill(ncID(k+j),nf90_NOFILL,izerr)
          IF (istat /= nf90_noerr) THEN
            ierror  = 10209_iintegers 
            yerrmsg = TRIM(nf90_strerror(istat))
            RETURN
          ENDIF
          
          ! Define global attributes
          READ(ydate_ini(1:4),'(i4)') iztime_tmp
          istat = nf90_put_att(ncID(k+j), nf90_GLOBAL, 'ref_year', iztime_tmp)
          IF (istat /= nf90_noerr) THEN
            ierror  = 10212_iintegers 
            yerrmsg = TRIM(nf90_strerror(istat))
            RETURN
          ENDIF
          READ(ydate_ini(5:6),'(i2)') iztime_tmp
          istat = nf90_put_att(ncID(k+j), nf90_GLOBAL, 'ref_month', iztime_tmp)
          IF (istat /= nf90_noerr) THEN
            ierror  = 10215_iintegers 
            yerrmsg = TRIM(nf90_strerror(istat))
            RETURN
          ENDIF
          READ(ydate_ini(7:8),'(i2)') iztime_tmp
          istat = nf90_put_att(ncID(k+j), nf90_GLOBAL, 'ref_day', iztime_tmp)
          IF (istat /= nf90_noerr) THEN
            ierror  = 10218_iintegers 
            yerrmsg = TRIM(nf90_strerror(istat))
            RETURN
          ENDIF
          READ(ydate_ini(9:10),'(i2)') iztime_tmp
          istat = nf90_put_att(ncID(k+j), nf90_GLOBAL, 'ref_hour', iztime_tmp)
          IF (istat /= nf90_noerr) THEN
            ierror  = 10221_iintegers 
            yerrmsg = TRIM(nf90_strerror(istat))
            RETURN
          ENDIF
          READ(ydate_ini(11:12),'(i2)') iztime_tmp
          istat = nf90_put_att(ncID(k+j), nf90_GLOBAL, 'ref_min', iztime_tmp)
          IF (istat /= nf90_noerr) THEN
            ierror  = 10224_iintegers 
            yerrmsg = TRIM(nf90_strerror(istat))
            RETURN
          ENDIF
          READ(ydate_ini(13:14),'(i2)') iztime_tmp
          istat = nf90_put_att(ncID(k+j), nf90_GLOBAL, 'ref_sec', iztime_tmp)
          IF (istat /= nf90_noerr) THEN
            ierror  = 10227_iintegers 
            yerrmsg = TRIM(nf90_strerror(istat))
            RETURN
          ENDIF
          rztime_tmp = REAL(nzout * ninc_out_traj - 1_iintegers,wp) * dt
          istat = nf90_put_att(ncID(k+j), nf90_GLOBAL, 'duration_in_sec', rztime_tmp)
          IF (istat /= nf90_noerr) THEN
            ierror  = 10230_iintegers 
            yerrmsg = TRIM(nf90_strerror(istat))
            RETURN
          ENDIF
          istat = nf90_put_att(ncID(k+j), nf90_GLOBAL, 'pollon', pollon)
          IF (istat /= nf90_noerr) THEN
            ierror  = 10233_iintegers 
            yerrmsg = TRIM(nf90_strerror(istat))
            RETURN
          ENDIF
          istat = nf90_put_att(ncID(k+j), nf90_GLOBAL, 'pollat', pollat)
          IF (istat /= nf90_noerr) THEN
            ierror  = 10236_iintegers 
            yerrmsg = TRIM(nf90_strerror(istat))
            RETURN
          ENDIF
          istat = nf90_put_att(ncID(k+j), nf90_GLOBAL, 'output_timestep_in_sec', ninc_out_traj*dt)
          IF (istat /= nf90_noerr) THEN
            ierror  = 10239_iintegers 
            yerrmsg = TRIM(nf90_strerror(istat))
            RETURN
          ENDIF
          
          ! Define dimensions
          istat = nf90_def_dim(ncID(k+j), 'id', nztraj, IdDimId(k+j))
          IF (istat /= nf90_noerr) THEN
            ierror  = 10242_iintegers 
            yerrmsg = TRIM(nf90_strerror(istat))
            RETURN
          ENDIF
          istat = nf90_def_dim(ncID(k+j), 'time', nzout, TimeDimId)
          IF (istat /= nf90_noerr) THEN
            ierror  = 10245_iintegers 
            yerrmsg = TRIM(nf90_strerror(istat))
            RETURN
          ENDIF
          
          !----------------------------------------------------------------------
          !- Section 2.2: Define variables and their attributes
          !----------------------------------------------------------------------
          
          ! Time
          istat = nf90_def_var(ncID(k+j), 'time', nf90_FLOAT, (/ TimeDimID /),  & 
                               varTimeId(k+j))
          IF (istat /= nf90_noerr) THEN
            ierror  = 10248_iintegers 
            yerrmsg = TRIM(nf90_strerror(istat))
            RETURN
          ENDIF
          istat=nf90_put_att(ncID(k+j), varTimeId(k+j), "standard_name", "time")
          IF (istat /= nf90_noerr) THEN
            ierror  = 10251_iintegers
            yerrmsg = TRIM(nf90_strerror(ierror))
            RETURN
          ENDIF
          istat=nf90_put_att(ncID(k+j), varTimeId(k+j), "long_name", "time")
          IF (istat /= nf90_noerr) THEN
            ierror  = 10254_iintegers
            yerrmsg = TRIM(nf90_strerror(ierror))
            RETURN
          ENDIF
          istat=nf90_put_att(ncID(k+j), varTimeId(k+j), "units", "seconds")
          IF (istat /= nf90_noerr) THEN
            ierror  = 10257_iintegers
            yerrmsg = TRIM(nf90_strerror(ierror))
            RETURN
          ENDIF
          
          
          ! Rotated longitudes
          istat = nf90_def_var(ncID(k+j), 'longitude', nf90_FLOAT,              &
                              (/ IdDimID(k+j), TimeDimID /), varLonId(k+j))
          IF (istat /= nf90_noerr) THEN
            ierror  = 10260_iintegers 
            yerrmsg = TRIM(nf90_strerror(istat))
            RETURN
          ENDIF
          istat=nf90_put_att(ncID(k+j), varLonId(k+j), "standard_name", "grid_longitude")
          IF (istat /= nf90_noerr) THEN
            ierror  = 10263_iintegers
            yerrmsg = TRIM(nf90_strerror(ierror))
            RETURN
          ENDIF
          istat=nf90_put_att(ncID(k+j), varLonId(k+j), "long_name", "rotated longitudes")
          IF (istat /= nf90_noerr) THEN
            ierror  = 10266_iintegers
            yerrmsg = TRIM(nf90_strerror(ierror))
            RETURN
          ENDIF
          istat=nf90_put_att(ncID(k+j), varLonId(k+j), "units", "degrees")
          IF (istat /= nf90_noerr) THEN
            ierror  = 10269_iintegers
            yerrmsg = TRIM(nf90_strerror(ierror))
            RETURN
          ENDIF
          
          
          ! Rotated latitudes
          istat = nf90_def_var(ncID(k+j), 'latitude', nf90_FLOAT,               &
                              (/ IdDimID(k+j), TimeDimID /), varLatId(k+j))
          IF (istat /= nf90_noerr) THEN
            ierror  = 10272_iintegers 
            yerrmsg = TRIM(nf90_strerror(istat))
            RETURN
          ENDIF
          istat=nf90_put_att(ncID(k+j), varLatId(k+j), "standard_name", "grid_latitude")
          IF (istat /= nf90_noerr) THEN
            ierror  = 10275_iintegers
            yerrmsg = TRIM(nf90_strerror(ierror))
            RETURN
          ENDIF
          istat=nf90_put_att(ncID(k+j), varLatId(k+j), "long_name", "rotated latitudes")
          IF (istat /= nf90_noerr) THEN
            ierror  = 10278_iintegers
            yerrmsg = TRIM(nf90_strerror(ierror))
            RETURN
          ENDIF
          istat=nf90_put_att(ncID(k+j), varLatId(k+j), "units", "degrees")
          IF (istat /= nf90_noerr) THEN
            ierror  = 10281_iintegers
            yerrmsg = TRIM(nf90_strerror(ierror))
            RETURN
          ENDIF
          
          
          ! Height 
          istat = nf90_def_var(ncID(k+j), 'z', nf90_FLOAT,                      &
                              (/ IdDimID(k+j), TimeDimID /), varZId(k+j))
          IF (istat /= nf90_noerr) THEN
            ierror  = 10284_iintegers 
            yerrmsg = TRIM(nf90_strerror(istat))
            RETURN
          ENDIF
          istat=nf90_put_att(ncID(k+j), varZId(k+j), "standard_name", "height")
          IF (istat /= nf90_noerr) THEN
            ierror  = 10287_iintegers
            yerrmsg = TRIM(nf90_strerror(ierror))
            RETURN
          ENDIF
          istat=nf90_put_att(ncID(k+j), varZId(k+j), "long_name", "height above mean sea level")
          IF (istat /= nf90_noerr) THEN
            ierror  = 10290_iintegers
            yerrmsg = TRIM(nf90_strerror(ierror))
            RETURN
          ENDIF
          istat=nf90_put_att(ncID(k+j), varZId(k+j), "units", "m AMSL")
          IF (istat /= nf90_noerr) THEN
            ierror  = 10293_iintegers
            yerrmsg = TRIM(nf90_strerror(ierror))
            RETURN
          ENDIF
          
          
          ! Traced variables
          DO i=1,ntracevar_traj
            istat = nf90_def_var(ncID(k+j), TRIM(traj_trace(i)%name), nf90_FLOAT, &
                                (/ IdDimID(k+j), TimeDimID /), varId(i,k+j))
            IF (istat /= nf90_noerr) THEN
              ierror  = 10302_iintegers 
              yerrmsg = TRIM(nf90_strerror(istat))
              RETURN
            ENDIF
            istat=nf90_put_att(ncID(k+j),varId(i,k+j),"standard_name",TRIM(traj_trace(i)%sdname))
            IF (istat /= nf90_noerr) THEN
              ierror  = 10305_iintegers 
              yerrmsg = TRIM(nf90_strerror(istat))
              RETURN
            ENDIF
            istat=nf90_put_att(ncID(k+j),varId(i,k+j),"long_name",TRIM(traj_trace(i)%lgname))
            IF (istat /= nf90_noerr) THEN
              ierror  = 10308_iintegers 
              yerrmsg = TRIM(nf90_strerror(istat))
              RETURN
            ENDIF
            istat=nf90_put_att(ncID(k+j),varId(i,k+j),"units",TRIM(traj_trace(i)%units))
            IF (istat /= nf90_noerr) THEN
              ierror  = 10311_iintegers 
              yerrmsg = TRIM(nf90_strerror(istat))
              RETURN
            ENDIF
          ENDDO
          
          
          ! End of definition section
          istat = nf90_enddef(ncID(k+j))
          IF (istat /= nf90_noerr) THEN
            ierror  = 10314_iintegers 
            yerrmsg = TRIM(nf90_strerror(istat))
            RETURN
          ENDIF
 
        !----------------------------------------------------------------------
        !- Section 3: Inquire variables from previous output in case of restart
        !----------------------------------------------------------------------

        ELSE

          !Check if trajectory files from former simulation are present
          istat = nf90_open(TRIM(TRIM(ydir_traj)//'/'//TRIM(yzout_filename)), &
                            IOR(NF90_WRITE,nf90_64bit_offset), ncID(k+j))
          IF (istat /= nf90_noerr) THEN
            ierror  = 10315_iintegers
            yerrmsg = 'TRAJECTORY module: Former trajectory output file missing.'
            RETURN
          ENDIF
          !Read variable IDs for t, lon, lat, z
          istat = nf90_inq_varid(ncID(k+j),"time",varTimeId(k+j))
          istat = nf90_inq_varid(ncID(k+j),"longitude",varLonId(k+j))
          istat = nf90_inq_varid(ncID(k+j),"latitude",varLatId(k+j))
          istat = nf90_inq_varid(ncID(k+j),"z" ,varZId(k+j))
          !Read variable IDs of traced variables
          DO i = 1,ntracevar_traj
            istat = nf90_inq_varid(ncID(k+j),TRIM(traj_trace(i)%name),varId(i,k+j))
          ENDDO

        ENDIF

      ENDDO
      j = j+nzfiles
    ENDDO
    n_files_tot = j

  ENDIF
#endif

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE output_traj_init

!==============================================================================
!+ euler time-step for new trajectory position
!------------------------------------------------------------------------------

SUBROUTINE euler_ts(tra_id,idebug)

!------------------------------------------------------------------------------
!
! Description:
!   New trajectory position is calculated based on the actual wind field.
!
! Method:
!   The new position is calculated with an Euler forward integration. The 
!   wind-fields are interpolated to this position and with this new velocities 
!   the trajectory position is calculated again. This procedure is performed
!   n times, as specified by numit (see below). The number of iteration is set
!   to 3 by default since tests have shown that 3 is enough to reach a good accuracy
!   in most cases. In some convective studies however, one might want to increase
!   this number in order to reach a good accuracy.
!
!   A parameter ljump prevents the "loss" of trajectories at the lower boundary, 
!   i.e.it prevents that the trajectories "enter" the ground. 
!   For any trajectory reaching the ground, its vertical position is reset to
!   hjump [m above ground] when ljump is set to TRUE (default, see below).
!   Tests have shown that in a normal case, about 40% of the trajectories
!   are "lost" at the lower boundary if this parameter is not set.
!
!   For more details, look at Miltenberger A.K. et al.: 
!   "An online trajectory module (version 1.0) for the nonhydrostatic numerical weather prediction model COSMO"
!   in Geosci. Model Dev., 6, 1989-2004, 2013
!   doi:10.5194/gmd-6-1989-2013
!------------------------------------------------------------------------------

! Declarations:

! Subroutine arguments:
! ---------------------
INTEGER (KIND=iintegers), INTENT(IN)     ::  &
  tra_id,                 & ! id of trajectory to be integrated  
  idebug                    ! debug level

! Local variables:
! ----------------

INTEGER (KIND=iintegers)                 ::  &
  j                         ! loop indices                     

REAL    (KIND=wp)                        ::  &
  u0, v0, w0,             & ! interpolated wind speed at initial trajectory
                            ! location                                           [m/s]
  u1, v1, w1,             & ! interpolated wind speed at new trajectory location [m/s]
  uu, vv, ww,             & ! wind speed for advection                           [m/s]
  hsurft,                 & ! surface elevation at new traj. location            [ m ]
  topt,                   & ! highest level at new traj. location                [ m ]
  tra_cur(3)                ! current trajectory position (before the Euler integration)          
                            !        (rotlon, rotlat, z) [deg, deg, m]  

REAL    (KIND=wp)                        ::  &
  x_trajpos,              & ! x-position of the trajectory in the grid index space [-] 
  y_trajpos,              & ! y-position of the trajectory in the grid index space [-]
  z_trajpos,              & ! z-position of the trajectory in the grid index space [-]
  zh_trajpos                ! z-position of the trajectory in the grid index space [-]


LOGICAL                                  ::  &
  lzleft_hori               ! has left horizontal total domain?

REAL    (KIND=wp),     PARAMETER         ::  &
  deltay = 1.112E5_wp,    & ! = (2*pi*earth_radius/360) [m/deg]
  pi     = 3.1415927_wp     !                           [-]


INTEGER(KIND=iintegers), PARAMETER       ::  &
  numit = 3_iintegers       ! Number of iterations          [-]
                            ! Tests have shown that in most cases, 3 iterations
                            ! are sufficient in most cases. For convective cases, higher
                            ! numbers of iterations might bring an added value.

LOGICAL, PARAMETER                       ::  &
  ljump = .TRUE.            ! This parameter resets the vertical position
                            ! of a trajectory to hjump [m above ground] each time 
                            ! a trajectory reaches the ground.

REAL(KIND=wp),     PARAMETER       ::  &
  hjump = 10.0_wp       ! height above ground at which trajectory is reset 
                            ! in case it reaches the ground and ljump=.TRUE.  [m] 

! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin subroutine euler_ts
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1 - Interpolate wind fields to actual position traj(tra_id,:)
!------------------------------------------------------------------------------

  ! Store current trajectory position
  tra_cur = traj(tra_id,:)

  ! Compute current position in the grid index space 
  !    (x_trajpos, y_trajpos, z_trajpos, zh_trajpos)
  CALL indpos(startlon_tot, startlat_tot, dlon, dlat, traj(tra_id,1),         &
              traj(tra_id,2), traj(tra_id,3), x_trajpos, y_trajpos,           &
              z_trajpos, zh_trajpos,debug=.true.)

  ! Interpolate wind speed at current position
  u0 = int3d(u(:,:,:,nnow), ie, je, ke  , x_trajpos, y_trajpos, z_trajpos , 1_iintegers)
  IF (l2dim) THEN
    v0 = 0.0_wp
  ELSE
    v0 = int3d(v(:,:,:,nnow), ie, je, ke  , x_trajpos, y_trajpos, z_trajpos , 2_iintegers)
  ENDIF
  w0 = int3d(w(:,:,:,nnow), ie, je, ke+1, x_trajpos, y_trajpos, zh_trajpos, 0_iintegers) 
  
!------------------------------------------------------------------------------
! Section 2 - Iterative calculation of new position
! -----------------------------------------------------------------------------

  DO j = 1, numit

    ! Compute new position in the grid index space
    CALL indpos(startlon_tot, startlat_tot, dlon, dlat, traj(tra_id,1),       &
                traj(tra_id,2), traj(tra_id,3), x_trajpos, y_trajpos,         &
                z_trajpos, zh_trajpos)

    ! Interpolate wind speed at new position    
    u1 = int3d(u(:,:,:,nnew), ie, je, ke  , x_trajpos, y_trajpos, z_trajpos , 1_iintegers)
    IF (l2dim) THEN
      v1 = 0.0_wp
    ELSE
      v1 = int3d(v(:,:,:,nnew), ie, je, ke  , x_trajpos, y_trajpos, z_trajpos , 2_iintegers)
    ENDIF
    w1 = int3d(w(:,:,:,nnew), ie, je, ke+1, x_trajpos, y_trajpos, zh_trajpos, 0_iintegers)

    ! Compute wind speed for advection (average of wind before and after
    ! integration)
    uu = 0.5_wp * (u0 + u1)
    vv = 0.5_wp * (v0 + v1)
    ww = 0.5_wp * (w0 + w1)

    ! Compute new position in physical space (proper integration)
    traj(tra_id,1) = tra_cur(1) + uu*dt/(deltay*cos(traj(tra_id,2)*pi/180.0_wp))
    traj(tra_id,2) = tra_cur(2) + vv*dt/deltay
    traj(tra_id,3) = tra_cur(3) + ww*dt

    ! Recompute new position in the grid index space
    CALL indpos(startlon_tot, startlat_tot, dlon, dlat, traj(tra_id,1),       &
                traj(tra_id,2), traj(tra_id,3), x_trajpos, y_trajpos,         &
                z_trajpos, zh_trajpos)

    ! Interpolate orography and domain top at new position 
    hsurft = int2d(hhl(:,:,:), ie, je, x_trajpos, y_trajpos, REAL(ke+1,wp))
    topt   = int2d(hhl(:,:,:), ie, je, x_trajpos, y_trajpos, REAL(1   ,wp))

    ! Check if trajectory has left the global computation domain
    lzleft_hori = lquit_hori_tot_domain(traj(tra_id,1), traj(tra_id,2),       & 
                                startlon_tot, startlat_tot, dlon, dlat,       &
                                nboundlines, ie_tot, je_tot)
    IF (lzleft_hori .OR. (traj(tra_id,3) > topt))  THEN
      IF (idebug > 50_iintegers) THEN
        PRINT *,'*** TRAJECTORY MODULE                                ***'
        PRINT *,'Trajectory ', tra_id,' has left the domain (',               &
                startlon_tot+dlon*nboundlines,':',                            &
                startlon_tot+dlon*(ie_tot-nboundlines),',',                   &
                startlat_tot+dlat*nboundlines,':',                            &
                startlat_tot+dlat*(je_tot-nboundlines),',', topt,')',         &
                ' at an horizontal border or at the top of the domain.'     
      ENDIF
      istat_traj(tra_id) = idead
      RETURN
    ELSEIF ((traj(tra_id,3) < hsurft)) THEN
      IF ( j /= numit ) THEN
        ! Reset the trajectory to the ground level as long as we are iterating
        traj(tra_id,3) = hsurft
      ELSE
        ! At last iteration either reset it to 10m above ground or kill it
        IF (ljump) THEN
          traj(tra_id,3) = hsurft+hjump
        ELSE
          IF (idebug > 50_iintegers) THEN
            PRINT *,'*** TRAJECTORY MODULE                                ***'
            PRINT *,'Trajectory ', tra_id,' has left the domain.',               &
                    ' It has reached the ground: ', hsurft,'.',                  &
                    ' The current position is: (',traj(tra_id,1),',',            &
                    traj(tra_id,2),',',traj(tra_id,3),').',                      &
                    ' The ground at the 4 surrounding horizontal points is at: ',&
                    hhl(FLOOR(x_trajpos)  ,FLOOR(y_trajpos)  ,ke+1),'//',        &
                    hhl(FLOOR(x_trajpos)+1,FLOOR(y_trajpos)  ,ke+1),'//',        &
                    hhl(FLOOR(x_trajpos)+1,FLOOR(y_trajpos)+1,ke+1),'//',        &
                    hhl(FLOOR(x_trajpos)  ,FLOOR(y_trajpos)+1,ke+1), ' meters.'
          ENDIF
          istat_traj(tra_id) = idead
        ENDIF
      ENDIF
    ENDIF

  ENDDO
  
!------------------------------------------------------------------------------
! End of subroutine euler_ts
!------------------------------------------------------------------------------

END SUBROUTINE euler_ts

!==============================================================================
!+ define 8 surrounding PEs for each PE
!------------------------------------------------------------------------------

SUBROUTINE define_8neigh_proc(my_cart_8neigh,ierror,yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   Define the 8 surroundings PEs for each PE.
!   A trajectory currently located on a given PE could travel to one of those
!   8 PEs after having been advected to its new position.
!
! Method:
!   The layout of the PEs in the cartesian communicator is not known from the
!   user side and might change depending on the MPI library implementation.
!   We thus cannot make any assumption on this layout and cannot address the
!   PE positions directly.
!   my_cart_pos gives the position (coordinates) of the PE in the cartesian
!   communicator and a call to MPI_CART_RANK provides the rank (ID) 
!   for the neighbors PEs.
!      
!------------------------------------------------------------------------------

! Declarations:

! Subroutine arguments:
! ---------------------
INTEGER   (KIND=iintegers), INTENT(OUT)                 ::  &
  my_cart_8neigh(:) ! 8 surrounding PEs  

INTEGER   (KIND=iintegers), INTENT(OUT)                 ::  &
  ierror            ! error status

CHARACTER (LEN= *),         INTENT(OUT)                 ::  &
  yerrmsg           ! error message

! Local variables:
!----------------
INTEGER   (KIND=iintegers)                              ::  & 
  i,           &    ! loop index
  izerror           ! local error code

INTEGER   (KIND=iintegers)                              ::  & 
  ij (2,8)          ! position on cartesian grid in i and j direction for
                    ! all neighboring PEs


! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin subroutine define_8neigh_proc
!------------------------------------------------------------------------------

  ierror            =  0_iintegers
  yerrmsg           =  '    '
  izerror           =  0_iintegers
  my_cart_8neigh(:) = -1_iintegers


  ! Compute indices for all neighbor PEs
  !-------------------------------------

  ! N neighbor PE
  ij (1,idxpe_n)  =  my_cart_pos(1)
  ij (2,idxpe_n)  =  my_cart_pos(2)+1
  ! NE neighbor PE
  ij (1,idxpe_ne) =  my_cart_pos(1)+1
  ij (2,idxpe_ne) =  my_cart_pos(2)+1
  ! E neighbor PE
  ij (1,idxpe_e)  =  my_cart_pos(1)+1
  ij (2,idxpe_e)  =  my_cart_pos(2)
  ! SE neighbor PE
  ij (1,idxpe_se) =  my_cart_pos(1)+1
  ij (2,idxpe_se) =  my_cart_pos(2)-1
  ! S neighbor PE
  ij (1,idxpe_s)  =  my_cart_pos(1)
  ij (2,idxpe_s)  =  my_cart_pos(2)-1
  ! SW neighbor PE
  ij (1,idxpe_sw) =  my_cart_pos(1)-1
  ij (2,idxpe_sw) =  my_cart_pos(2)-1
  ! W neighbor PE
  ij (1,idxpe_w)  =  my_cart_pos(1)-1
  ij (2,idxpe_w)  =  my_cart_pos(2)
  ! NW neighbor PE
  ij (1,idxpe_nw) =  my_cart_pos(1)-1
  ij (2,idxpe_nw) =  my_cart_pos(2)+1


  ! Determine the rank for all neighbor PEs
  !----------------------------------------

  DO i=1,8
    IF (ij(1,i) < 0_iintegers) THEN
      IF (lperi_x) THEN
        ij(1,i) = ij(1,i) + nprocx
      ELSE
        ij(1,i) = -1_iintegers
      ENDIF
    ENDIF
    IF (ij(1,i) > nprocx-1) THEN
      IF (lperi_x) THEN
        ij(1,i) = ij(1,i) - nprocx
      ELSE
        ij(1,i) = -1_iintegers
      ENDIF
    ENDIF
    IF (ij(2,i) < 0_iintegers) THEN
      IF (lperi_y) THEN
        ij(2,i) = ij(2,i) + nprocy
      ELSE
        ij(2,i) = -1_iintegers
      ENDIF
    ENDIF
    IF (ij(2,i) > nprocy-1) THEN
      IF (lperi_y) THEN
        ij(2,i) = ij(2,i) - nprocy 
      ELSE
        ij(2,i) = -1_iintegers
      ENDIF
    ENDIF
    IF (ij(1,i) < 0_iintegers .OR. ij(2,i) < 0_iintegers) THEN
      my_cart_8neigh(i) = -1
    ELSE
      CALL MPI_CART_RANK (icomm_cart, ij(:,i), my_cart_8neigh(i), izerror)
       IF (izerror /= 0_iintegers) THEN
         ierror  = 10330_iintegers
         yerrmsg = 'Error while determining the 8 neighbor PEs'
       ENDIF
    ENDIF
  ENDDO
 
!------------------------------------------------------------------------------
! End of subroutine define_8neigh_proc
!------------------------------------------------------------------------------

END SUBROUTINE define_8neigh_proc

!==============================================================================
!+ preparation of communication with other processes
!------------------------------------------------------------------------------

SUBROUTINE prep_comm_traj(traj_id, ntrajbuf, buf)

!------------------------------------------------------------------------------
!
! Description:
!   Prepare the buffer for sending the position and ID of the trajectories 
!   leaving the PE subdomain to the neighbor PEs.
!
! Method:
!   The position of the trajectory is compared against the subdomain limits
!   and in case one is exceeded the neighbor PE is identified.
!   The ID and the 3 components of the trajectory position are packed in
!   the buffer corresponding to the identifed neighbor PE.
!      
!------------------------------------------------------------------------------

! Declarations:

! Subroutine arguments:
! ---------------------
INTEGER (KIND=iintegers), INTENT(IN)                 ::  &
  traj_id                   ! id of trajectory to be integrated  

INTEGER (KIND=iintegers), INTENT(INOUT)              ::  &
  ntrajbuf(:)               ! number of trajectories to be send
 
REAL (KIND=wp),           INTENT(INOUT)              ::  &
  buf(:,:)                  ! positions of trajectories to be send
  
! Local variables:
! ----------------
INTEGER (KIND=iintegers)                             ::  &
  izpe_id,               &  ! index of the neighbor PE
  izshift(8),            &  ! position in the buffer
  i                         ! loop index


! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin subroutine prep_comm_traj
!------------------------------------------------------------------------------

  izshift(:) = 0_iintegers
  izpe_id    = 0_iintegers

  ! Identify to which neighbor PE the trajectory has to be sent
  IF ((istat_traj(traj_id) == ialive)) THEN 
    izpe_id = ipe_location(traj(traj_id,1), traj(traj_id,2), startlon_tot,    &
                           startlat_tot, isubpos(my_cart_id,1),               &
                           isubpos(my_cart_id,2), dlon, dlat, idxpe_n,        &
                           idxpe_ne, idxpe_e,  idxpe_se, idxpe_s,             &
                           idxpe_sw, idxpe_w, idxpe_nw, idxpe_this,           &
                           nboundlines, ie, je)

    ! Compute the number of traj to be sent to a given neighbor PE,
    ! compute the position in the corresponding buffer,
    ! write the traj ID, as well as the traj position (lon, lat,z)
    ! in the correct position of the corresponding buffer 
    DO i=1,8
      IF (i==izpe_id) THEN
        ntrajbuf(i)         = ntrajbuf(i) + 1
        izshift(i)          = (ntrajbuf(i)-1) * 4 + 1
        buf(izshift(i),i)   = REAL(traj_id,wp)
        buf(izshift(i)+1:izshift(i)+3,i) = traj(traj_id,:)
        ! deactivate trajectory leaving this PE
        !istat_traj(traj_id) == idead
      ENDIF
    ENDDO
  ENDIF

!------------------------------------------------------------------------------
! End of subroutine prep_comm_traj
!------------------------------------------------------------------------------

END SUBROUTINE prep_comm_traj

!==============================================================================
!+ Send position of trajectories leaving the subdomain to new processor
!------------------------------------------------------------------------------

SUBROUTINE do_comm_traj(ntrajbuf,buf,irank_8neigh,ierror,yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   Communication between processors.
!
! Method:
!   The communication of the trajectories is done "pairwise".
!   Two buffers (e.g. North and South)are communicated at the same time.
!   This should be faster than communicating in each direction one after the
!   other.
!   On the other hand, less memory needs to be allocated as compared to a 
!   solution where all communications would be done at the same time.
!   Here we allocate a buffer of size 2 and re-use it 4 times (instead
!   of having a buffer of size 8). This is relevant since it is allocated on all
!   processors with the number of trajectories times 4!
!   This solution should thus be a compromise between performance in terms of
!   communication speed and memory consumption.
!   This could easily be changed if desired.
!
!------------------------------------------------------------------------------

! Declarations:

! Subroutine arguments:
! ---------------------
INTEGER (KIND=iintegers), INTENT(IN)     ::  &
  ntrajbuf    (:),           & ! number of trajectories to be send
  irank_8neigh(:)              ! rank of the 8 neighbor PEs
 
REAL    (KIND=wp),        INTENT(IN)     ::  &
  buf(:,:)                     ! trajectory data to be sent

INTEGER (KIND=iintegers), INTENT(OUT)    ::  &
  ierror                       ! error index

CHARACTER (LEN= *),       INTENT(OUT)    ::  &
  yerrmsg                      ! error message


! Local variables:
! ----------------
INTEGER (KIND=iintegers)                 ::  &
  izerrstat,                 & ! local error code
  i, j, k,                   & ! loop indices
  iztraj_id,                 & ! id of the trajectory
  izpe(8),                   & ! indices for the 8 neighbor PEs 
  izshift                      ! position in the buffer

REAL (KIND=wp),     ALLOCATABLE          ::  &
  zrcv_buf(:,:)                ! receive buffer for trajectory data


! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin subroutine do_comm_traj
!------------------------------------------------------------------------------

  ierror    = 0_iintegers
  yerrmsg   = ' '
  izerrstat = 0_iintegers
  izshift   = 0_iintegers

!------------------------------------------------------------------------------
! Section 1: Preparation
!------------------------------------------------------------------------------

  ALLOCATE(zrcv_buf(4*ntraj,2),STAT=izerrstat)
  IF (izerrstat /= 0_iintegers) THEN
    ierror  = 10340_iintegers
    yerrmsg = 'Allocation of zrcv_buf failed.'
    RETURN
  ENDIF
  zrcv_buf(:,:) = -999.0_wp

!------------------------------------------------------------------------------
! Section 2: Communication
!------------------------------------------------------------------------------

  izpe(1) = idxpe_n
  izpe(2) = idxpe_s
  izpe(3) = idxpe_e
  izpe(4) = idxpe_w
  izpe(5) = idxpe_ne
  izpe(6) = idxpe_sw
  izpe(7) = idxpe_se
  izpe(8) = idxpe_nw

  DO j = 1,4
    ! Do the communication of the 2 buffers
    CALL comm_traj_2buf(ntrajbuf(:), buf(:,:), irank_8neigh(:),               &
                        zrcv_buf(:,:), izpe((2*j)-1:2*j), ierror, yerrmsg)
    IF (ierror /= 0_iintegers) THEN
      RETURN
    ENDIF
    DO i = 1,2
      ! Unpack the 2 receive buffers and update trajectories
      IF (MAXVAL(zrcv_buf(:,i)) > -999.0_wp) THEN
        DO k = 1,ntraj
          izshift = 4*(k-1)+1
          IF (zrcv_buf(izshift,i) /= -999.0_wp) THEN
            iztraj_id         = NINT(zrcv_buf(izshift,i),KIND=iintegers)
            traj(iztraj_id,:) = zrcv_buf(izshift+1:izshift+3,i)
            ! activate trajectory
            !istat_traj(traj_id) == ialive
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section 3: Cleanup
!------------------------------------------------------------------------------

  DEALLOCATE (zrcv_buf, STAT = izerrstat)
  IF (izerrstat /=0) THEN
    ierror  = 10343_iintegers
    yerrmsg = 'Deallocation of zrcv_buf failed.'
    RETURN
  ENDIF

!------------------------------------------------------------------------------
! End of subroutine do_comm_traj
!------------------------------------------------------------------------------

END SUBROUTINE do_comm_traj

!==============================================================================
!+ Non-blocking communication of 2 buffers
!------------------------------------------------------------------------------

SUBROUTINE comm_traj_2buf(ntrajbuf,buf,irank_neigh,rcv_buf,ipe,ierror,yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   Communication between processors.
!
! Method:
!   "Pair-wise" communication of the trajectory data.
!   A non-blocking strategy is used.
!   First the 2 non-blocking receive operations are pre-posted and then the
!   2 non-blocking send operations are posted. A waitall ensures that the
!   2 sends and the 2 receives are completed before the buffers are re-used
!   thus guarantying a safe communication.
!   
!------------------------------------------------------------------------------

! Declarations:

! Subroutine arguments:
! ---------------------
INTEGER (KIND=iintegers), INTENT(IN)          ::  &
  ntrajbuf   (:),        & ! number of trajectories to be sent
  irank_neigh(:),        & ! rank of the 2 neighbor PEs
  ipe        (:)           ! index for the 2 neighbor PEs

REAL    (KIND=wp),        INTENT(IN)          ::  &
  buf(:,:)                 ! trajectory data to be sent

REAL    (KIND=wp),        INTENT(INOUT)       ::  &
  rcv_buf(:,:)             ! receiving buffer for trajectory data

INTEGER (KIND=iintegers), INTENT(OUT)         ::  &
  ierror                   ! error index

CHARACTER (LEN= *),       INTENT(OUT)         ::  &
  yerrmsg                  ! error message

! Local variables:
! ----------------

INTEGER (KIND=iintegers)                      ::  &
  izerr,                 & ! local error index
  nzcom,                 & ! counter for the communications (4 in total)  
  izreq(4),              & ! request object for the 4 communications
  izst(MPI_STATUS_SIZE,4)  ! status of the 4 communications


! End of header
!==============================================================================

  ierror  = 0_iintegers
  yerrmsg = ' '

!------------------------------------------------------------------------------
! Begin subroutine comm_traj_2buf
!------------------------------------------------------------------------------

  ! Initialize variables
  izerr    = 0_iintegers
  nzcom    = 0_iintegers
  izreq(:) = 0_iintegers
  izst(:,:)= 0_iintegers


  ! Post receive from first neighbor PE (N, E, WE, SE for the 4 calls)
  IF (irank_neigh(ipe(1)) /= -1_iintegers) THEN
    nzcom = nzcom + 1_iintegers
    CALL MPI_IRECV(rcv_buf(:,1), 4*ntraj, imp_reals, irank_neigh(ipe(1)),&
                   202, icomm_cart, izreq(nzcom), izerr)  
    IF (izerr /= 0_iintegers) THEN
      ierror  = 10350_iintegers
      yerrmsg = '1st MPI_IRECV for trajectories failed.'
      RETURN
    ENDIF 
  ENDIF

  ! Post receive from second neighbor PE (S, W, SW, NW for the 4 calls)
  IF (irank_neigh(ipe(2)) /= -1_iintegers) THEN
    nzcom = nzcom + 1_iintegers
    CALL MPI_IRECV(rcv_buf(:,2), 4*ntraj, imp_reals, irank_neigh(ipe(2)),&
                   201, icomm_cart, izreq(nzcom), izerr)
    IF (izerr /= 0_iintegers) THEN
      ierror  = 10353_iintegers
      yerrmsg = '2nd MPI_IRECV for trajectories failed.'
      RETURN
    ENDIF 
  ENDIF

  ! Post send to first neighbor PE (N, E, WE, SE for the 4 calls)
  IF (irank_neigh(ipe(1)) /= -1_iintegers) THEN
    nzcom = nzcom + 1_iintegers
    CALL MPI_ISEND(buf(1:ntrajbuf(ipe(1)),ipe(1)), ntrajbuf(ipe(1)),         &
                   imp_reals, irank_neigh(ipe(1)), 201, icomm_cart,           &
                   izreq(nzcom), izerr)
    IF (izerr /= 0_iintegers) THEN
      ierror  = 10356_iintegers
      yerrmsg = '1st MPI_ISEND for trajectories failed.'
      RETURN
    ENDIF
  ENDIF

  ! Post send to second neighbor PE (S, W, SW, NW for the 4 calls)
  IF (irank_neigh(ipe(2)) /= -1_iintegers) THEN
    nzcom = nzcom + 1_iintegers
    CALL MPI_ISEND(buf(1:ntrajbuf(ipe(2)),ipe(2)), ntrajbuf(ipe(2)),         &
                   imp_reals, irank_neigh(ipe(2)), 202, icomm_cart,           &
                   izreq(nzcom), izerr)
    IF (izerr /= 0_iintegers) THEN
      ierror  = 10359_iintegers
      yerrmsg = '2nd MPI_ISEND for trajectories failed.'
      RETURN
    ENDIF
  ENDIF


  ! Wait until communications (2send and 2receive) for this PE are finished
  IF (nzcom > 0_iintegers) THEN
    CALL MPI_WAITALL(nzcom, izreq, izst, izerr)
    IF (izerr /= 0_iintegers) THEN
      ierror  = 10362_iintegers
      yerrmsg = 'MPI_WAITALL for trajectories failed.'
      RETURN
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! End of subroutine comm_traj_2buf
!------------------------------------------------------------------------------

END SUBROUTINE comm_traj_2buf


!==============================================================================
!+ interpolation of fields to new trajectory position
!------------------------------------------------------------------------------

SUBROUTINE intp_traj(nsendl,sendbufl)

!------------------------------------------------------------------------------
!
! Description:
!   Interpolates the traced variables to the actual trajectory location.
!
! Method:
!   A loop over all variables to trace is performed. According to their
!   leveltype and name, the staggering flag and vertical dimension are set. 
!   The variables that are not ready to be written out are computed.
!   In case the trajectory is located on the PE, its position in the index
!   space is computed and the values of the traced variables are interpolated
!   to the trajectory position.
!   
!------------------------------------------------------------------------------

! Declarations:

! Subroutine arguments:
! ---------------------
INTEGER (KIND=iintegers), INTENT(INOUT)       ::  &
  nsendl                       ! number of trajectories to be send

REAL    (KIND = wp),     INTENT(INOUT)        ::  &
  sendbufl(:)                  ! send buffer

! Local variables:
! ----------------
INTEGER (KIND=iintegers)                      ::  &
  izerror,           & ! local error code
  k, i,              & ! loop indices
  ind,               & ! position of the traj in the buffer
  iz_proc,           & ! PE index
  zkdim                ! local vertical dimension        

CHARACTER(LEN=255)                            ::  &
  yzerror,           & ! error message
  yzroutine            ! subroutine name

REAL    (KIND = wp)                           ::  &
  ind_x, ind_y, ind_z, ind_zh  ! position in the grid index space

REAL    (KIND=wp)                             ::  &
  zhelp1 (ie,je,ke)   , &      !
  zhelp2 (ie,je,ke)   , &      !
  zhelp3 (ie,je,ke)   , &      !
  zhelp4 (ie,je,ke)   , &      !
  zhelp2d(ie,je)      , &      ! help variables
  zqrs   (ie,je,ke)   , &      ! sum of precipitating species
  zrho   (ie,je,ke)   , &      ! air density 
  var_int(ie,je,ke+1)          ! variable to be traced

REAL    (KIND=wp)                             ::  &
  int_value                    ! values of the traced variables at the 
                               ! trajectory position (interpolated values)

! Tracer pointers
!----------------
REAL (KIND=wp),     POINTER                   ::  &
  qv  (:,:,:) => NULL(), &     ! QV at tlev=nnew
  qc  (:,:,:) => NULL(), &     ! QC at tlev=nnew
  qi  (:,:,:) => NULL(), &     ! QI at tlev=nnew
  qr  (:,:,:) => NULL(), &     ! QR at tlev=nnew
  qs  (:,:,:) => NULL(), &     ! QS at tlev=nnew
  qg  (:,:,:) => NULL()        ! QG at tlev=nnew


! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin subroutine intp_traj
!------------------------------------------------------------------------------


  izerror    = 0_iintegers
  yzerror    = '         '    
  yzroutine  = 'intp_traj'
  zkdim      = 0_iintegers

  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnew, ptr = qv)
  IF (izerror /= 0_iintegers) THEN
    yzerror = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerror, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nnew, ptr = qc)
  IF (izerror /= 0_iintegers) THEN
    yzerror = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerror, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qi, ptr_tlev = nnew, ptr = qi)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerror = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerror, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qr, ptr_tlev = nnew, ptr = qr)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerror = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerror, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qs, ptr_tlev = nnew, ptr = qs)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerror = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerror, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qg, ptr_tlev = nnew, ptr = qg)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerror = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerror, yzroutine)
  ENDIF

  DO k = 1,ntracevar_traj
    ind = 0_iintegers
    ! variable is ready to be printed out
    IF (traj_trace(k)%idef_stat == 1_iintegers) THEN
      IF (traj_trace(k)%rank == 4_iintegers) THEN
        IF (traj_trace(k)%levtyp == 109_iintegers) THEN
          traj_trace(k)%istag = 3_iintegers
          zkdim               = ke + 1_iintegers
          var_int(:,:,1:ke+1) = traj_trace(k)%p4(:,:,:,nnew)
        ELSE
          zkdim               = ke
          var_int(:,:,1:ke)   = traj_trace(k)%p4(:,:,:,nnew)
        ENDIF
        IF (TRIM(traj_trace(k)%name) == 'U') THEN
          traj_trace(k)%istag = 1_iintegers
        ELSEIF (TRIM(traj_trace(k)%name) == 'V') THEN
          traj_trace(k)%istag = 2_iintegers
        ENDIF
      ELSEIF (traj_trace(k)%rank == 3_iintegers) THEN
        IF (traj_trace(k)%levtyp == 109_iintegers) THEN
          traj_trace(k)%istag = 3_iintegers
          zkdim               = ke + 1_iintegers
          var_int(:,:,1:ke+1) = traj_trace(k)%p3(:,:,:)
        ELSE
          zkdim               = ke
          var_int(:,:,1:ke)   = traj_trace(k)%p3(:,:,:)
        ENDIF
      ELSE
        ! wrong but has been tested already
      ENDIF
    ! variable has to be computed
    ELSEIF (traj_trace(k)%idef_stat == 0_iintegers) THEN
      ! test name and according to name, compute and fill the elements of the data structure
      ! if they are different from the default values
      IF (TRIM(traj_trace(k)%name) == 'P') THEN
        var_int(:,:,1:ke) = (pp(:,:,:,nnew) + p0(:,:,:))
      ELSEIF (TRIM(traj_trace(k)%name) == 'RELHUM') THEN
        CALL calrelhum(var_int(:,:,1:ke), t(:,:,:,nnew), pp(:,:,:,nnew),      &
                       p0(:,:,:), qv(:,:,:), ie, je, ke,                      &
                       b1, b2w, b3, b4w, rdv, o_m_rdv )
      ELSEIF (TRIM(traj_trace(k)%name) == 'POT_VORTIC') THEN
        zqrs = qr
        IF (ASSOCIATED(qi)) THEN
          zqrs = zqrs + qi
        ENDIF
        IF ( ASSOCIATED(qs) ) THEN
          zqrs = zqrs + qs
        END IF
        IF ( ASSOCIATED(qg) ) THEN
          zqrs = zqrs + qg
        END IF
        CALL calrho(t(:,:,:,nnew), pp(:,:,:,nnew), qv, &
                    qc, zqrs, p0, zrho, ie, je, ke, r_d,     &
                    rvd_m_o )
        CALL calc_Theta_Tppp(t(:,:,:,nnew), pp(:,:,:,nnew), p0, ie, je, ke,   &
                             r_d, cp_d, zhelp4)
        CALL curl(ie, je, ke, eddlon, eddlat, r_earth, acrlat, tgrlat,        &
                  sqrtg_r_s, dzeta_dlam, dzeta_dphi, lmetr,                   &
                  u(:,:,:,nnew), v(:,:,:,nnew), w(:,:,:,nnew), wgtfac,        &
                  .FALSE., zhelp1, zhelp2, zhelp3)
        IF ( .NOT. lcori_deep ) THEN
          zhelp2d(:,:) = 0.0_wp
        ELSE
          zhelp2d(:,:) = fccos(:,:)
        ENDIF
        CALL potential_vorticity_rho(ie, je, ke, eddlon, eddlat, r_earth,     &
                   fc, zhelp2d, sqrtg_r_s, dzeta_dlam, dzeta_dphi,            &
                   zhelp1, zhelp2, zhelp3, lmetr, zhelp4,                     &
                   u(:,:,:,nnew), v(:,:,:,nnew), w(:,:,:,nnew), var_int(:,:,1:ke))
        var_int(:,:,1:ke) = var_int(:,:,1:ke) / zrho(:,:,1:ke)
      ELSE
        izerror = 10370_iintegers
        yzerror = ' ERROR  *** TRAJECTORY module: '//TRIM(traj_trace(k)%name)//&
                  ' cannot be traced (not implemented) *** '
        CALL model_abort(my_cart_id, izerror, yzerror, yzroutine)      
      ENDIF    
    ELSE
        izerror = 10373_iintegers
        yzerror = ' ERROR  *** TRAJECTORY module: variable status can only be 0 or 1 *** '
        CALL model_abort(my_cart_id, izerror, yzerror, yzroutine)            
    ENDIF

    DO i = 1,ntraj
      ! Check if trajectory is in my domain


      IF (istat_traj(i) == ialive) THEN
        iz_proc = ipe_location(traj(i,1), traj(i,2), startlon_tot,            &
                               startlat_tot, isubpos(my_cart_id,1),           &
                               isubpos(my_cart_id,2), dlon, dlat, idxpe_n,    &
                               idxpe_ne, idxpe_e,  idxpe_se, idxpe_s,         &
                               idxpe_sw, idxpe_w, idxpe_nw, idxpe_this,       &
                               nboundlines, ie, je)

        IF(iz_proc == idxpe_this) THEN
       
          ! Calculate position in the grid index space
          CALL indpos(startlon_tot,startlat_tot,dlon,dlat,traj(i,1),            &
                      traj(i,2),traj(i,3),ind_x,ind_y,ind_z,ind_zh)
          
          IF (traj_trace(k)%istag == 3_iintegers) THEN
            ! Interpolation to the trajectory position
            int_value = int3d(var_int,ie,je,zkdim,ind_x,ind_y,ind_zh, traj_trace(k)%istag)
          ELSE
            ! Interpolation to the trajectory position
            int_value = int3d(var_int,ie,je,zkdim,ind_x,ind_y,ind_z, traj_trace(k)%istag)
          ENDIF
          
          ! Prepare sendbuffers
          ind = ind+1_iintegers
          IF (k==1) THEN
            nsendl                = nsendl+4_iintegers+ntracevar_traj
            sendbufl(ind)         = REAL(i,wp)
            sendbufl(ind+1:ind+3) = traj(i,:)
            sendbufl(ind+4)       = int_value
          ELSE
            sendbufl(ind+3+k) = int_value
          ENDIF
          ind = ind+3_iintegers+ntracevar_traj     
        ENDIF   
      ENDIF 

    ENDDO

  ENDDO ! k-loop

!------------------------------------------------------------------------------
! End of subroutine intp_traj
!------------------------------------------------------------------------------

END SUBROUTINE intp_traj

!==============================================================================
!+ NetCDF output of the trajectories
!------------------------------------------------------------------------------

SUBROUTINE output_traj(nstep,dt,ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   Performs the output of the online trajectory data to the NetCDF files.
!
! Method:
!   First the traced variables are interpolated to the trajectory position
!   and the sendbuffer are prepared (in intp_traj).
!   Then all trajectories to be sent are gathered on PE0, which then writes
!   the data out to the NetCDF files. The files are open already (and kept
!   open during the whole trajectory simulation) so the data can be written
!   directly.
!   
!------------------------------------------------------------------------------

! Declarations:

! Subroutine arguments:
! ---------------------
 
REAL      (KIND=wp),        INTENT(IN)            :: &
  dt

INTEGER   (KIND=iintegers), INTENT(IN)            :: &
  nstep

INTEGER   (KIND=iintegers), INTENT(OUT)           :: &
  ierror  

CHARACTER (LEN=*),          INTENT(OUT)           :: &
  yerrmsg
 
! Local variables:
! ----------------
INTEGER   (KIND=iintegers)                        :: &
  izerr, izerrstat,      & ! local error codes 
  nsend,                 & ! number of trajectories to be sent
  i, k,                  & ! loop indices
  ind,                   & ! position of the trajectory in the tranew structure
  id_tra,                & ! ID of the trajectory
  istat                    ! NetCDF status variable
 
REAL      (KIND=wp)                               :: &
  time_now                 ! current time

REAL      (KIND=wp),        ALLOCATABLE           :: &
  tranew(:),             & ! receiving buffer
  sendbuf(:)               ! send buffer

INTEGER   (KIND=iintegers), ALLOCATABLE           :: &
  displs(:),             & ! displacement in receive buffer
  all_send(:)              ! number of datasets received form processes
 
REAL      (KIND=wp)                               :: &
  intp_fields(ntraj,ntracevar_traj) ! array of traced variables for output


! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin subroutine output_traj
!------------------------------------------------------------------------------

  yerrmsg   = ''
  ierror    = 0_iintegers
  izerr     = 0_iintegers
  izerrstat = 0_iintegers

!------------------------------------------------------------------------------
! Section 1: Prepare sending buffers and gather data to process 0
!------------------------------------------------------------------------------


  ! Allocate and initialize variables
  ALLOCATE (all_send(num_compute),STAT=izerrstat);izerr=izerr+izerrstat
  all_send(:) = 0_iintegers

  ALLOCATE (sendbuf(ntraj*(4+ntracevar_traj)),STAT=izerrstat);izerr=izerr+izerrstat
  sendbuf(:) = 0.0_wp

  ALLOCATE (displs(num_compute), STAT=izerrstat);izerr=izerr+izerrstat
  displs = 0_iintegers

  ALLOCATE (tranew(ntraj*(4+ntracevar_traj)), STAT=izerrstat);izerr=izerr+izerrstat
  tranew(:) = 0.0_wp

  IF (izerr /= 0_iintegers) THEN
    ierror  = 10400_iintegers
    yerrmsg = 'Allocation of variables in output_traj failed.'
    RETURN
  ENDIF

  nsend     = 0_iintegers


  ! Interpolate the traced variable values
  CALL intp_traj(nsend,sendbuf)


  IF (num_compute > 1_iintegers) THEN
    ! Gather number of trajectory datasets to be send from each process to 
    ! process 0
    CALL MPI_GATHER(nsend,1,imp_integers,all_send,1,imp_integers,0,        &
                    icomm_cart,izerr)
    IF (izerr /= 0_iintegers) THEN
      ierror  = 10403_iintegers
      yerrmsg = ' MPI_GATHER of nsend failed.'
    ENDIF
  ELSE
    all_send(1) = nsend
  ENDIF
 
  ! Prepare variables for GATHERV
  displs(1) = 0_iintegers
  DO i = 1,(num_compute-1)
    displs(i+1) = SUM(all_send(1:i))
  ENDDO
  IF (num_compute > 1_iintegers) THEN
    ! Gather all trajectory datasets to process 0
    CALL MPI_GATHERV(sendbuf(1:nsend),nsend,imp_reals,tranew,all_send,displs,&
                     imp_reals,0,icomm_cart,izerr)
    IF (izerr /= 0_iintegers) THEN
      ierror  = 10406_iintegers
      yerrmsg = ' MPI_GATHER of sendbuf failed.'
    ENDIF
  ELSE
    tranew(1:nsend) = sendbuf(1:nsend)
  ENDIF
                                      
  ! Update traj and fill intp_fields on process 0
  IF (my_cart_id == 0) THEN
    DO i=1,ntraj
      ind = (4+ntracevar_traj)*(i-1)+1_iintegers

      IF (tranew(ind) /= 0) THEN
        id_tra = INT(tranew(ind))
        traj(id_tra,:) = tranew(ind+1:ind+3)
        IF (ntracevar_traj > 0) THEN
          intp_fields(id_tra,:) = tranew(ind+4:ind+ntracevar_traj+3)
        ENDIF
      ENDIF
    ENDDO
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Write output to netcdf file
!------------------------------------------------------------------------------

#ifdef NETCDF
  IF (my_cart_id == 0_iintegers) THEN
    ! Determine writing position in NetCDF file and the time index
    time_now = REAL(nstep,wp)*dt
    
    ! Write data to NetCDF files
    DO k=1,n_files_tot
      ! Check whether trajectories in this file have already been initialized
      ! but check for ANY trajectory, not only for the first one
      IF (ANY(istat_traj(id_file(k,1):id_file(k,2)) == ialive))  THEN

        outstep(k) = outstep(k)+1

        istat = NF90_put_var(ncID(k),varTimeId(k),time_now,start=(/outstep(k)/))
        istat = NF90_put_var(ncID(k),varLonId(k),traj(id_file(k,1):id_file(k,2),1), &
                             start=(/1,outstep(k)/))
        istat = NF90_put_var(ncID(k),varLatId(k),traj(id_file(k,1):id_file(k,2),2), &
                             start=(/1,outstep(k)/))
        istat = NF90_put_var(ncID(k),varZId(k),traj(id_file(k,1):id_file(k,2),3),   &
                             start=(/1,outstep(k)/))
 
        IF (ntracevar_traj > 0) THEN
          DO i=1,ntracevar_traj
            istat = NF90_put_var(ncID(k),varId(i,k),    &
                    intp_fields(id_file(k,1):id_file(k,2),i),start=(/1,outstep(k)/))
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF
#endif

!------------------------------------------------------------------------------
! Section 3: Final actions
!------------------------------------------------------------------------------

  DEALLOCATE (all_send,STAT=izerrstat); izerr=izerr+izerrstat
  DEALLOCATE (sendbuf, STAT=izerrstat); izerr=izerr+izerrstat
  DEALLOCATE (displs,  STAT=izerrstat); izerr=izerr+izerrstat
  DEALLOCATE (tranew,  STAT=izerrstat); izerr=izerr+izerrstat
  IF (izerr /= 0_iintegers) THEN
    ierror  = 10409_iintegers
    yerrmsg = 'Deallocation of variables in output_traj failed.'
    RETURN
  ENDIF

!------------------------------------------------------------------------------
! End of subroutine output_traj
!------------------------------------------------------------------------------

END SUBROUTINE output_traj

!==============================================================================
!+ calculates position in grid index space from physical space
!------------------------------------------------------------------------------

SUBROUTINE indpos (minlon, minlat, deltalon, deltalat, lon_pos, lat_pos,      &
                   height_pos, x_pos, y_pos, z_pos, zh_pos, debug)

!------------------------------------------------------------------------------
!
! Description:
!   Calculate the trajectory position in grid index space from the physical 
!   trajectory position.
!
! Method:
!
! Naming conventions:
!   The physical space is described in terms of rotated longitudes, 
!   rotated latitudes and height above mean sea level. 
!   The grid is caracterized by the lon/lat of the lower left
!   corner (minlon,minlat) and by the grid spacing in both directions
!   (deltalon, deltalat). These trajectory positions in the physical space
!   are given by (lon_pos, lat_pos, height_pos).
!   
!   The grid space is considered as a cartesian space and positions in the 
!   grid space are denoted with (x_pos,y_pos,z_pos/zh_pos). x_pos and y_pos
!   correspond respectively to lon and lat in physical space.
!   z_pos is relative to the full levels and zh_pos to the half levels.
!   
!------------------------------------------------------------------------------

! Declarations:


! Subroutine arguments:
! ---------------------

REAL (KIND = wp),     INTENT(IN)   ::    &
  minlon, minlat,                  & !physical coordinates of the lower left 
                                     !corner of the domain                  [deg]
  deltalon, deltalat,              & !grid spacing                          [deg]
  lon_pos, lat_pos,                & !trajectory position along lon/lat     [deg]
  height_pos                         !trajectory position along vertic.     [ m ]

REAL (KIND = wp),     INTENT(OUT)  ::    &
  x_pos, y_pos,                    & !trajectory position in the grid index space [-] 
  z_pos, zh_pos                      !trajectory position in the grid index space [-] 

LOGICAL, OPTIONAL, INTENT(IN) :: debug

! Local variables:
! ----------------
INTEGER (KIND = iintegers) ::  &
  k         ! loop index   

REAL    (KIND = wp)        ::  &
  zt, zb    ! interpolated half-level or full-level
            ! at the horizontal location of the traj for k and k+1 resp.

INTEGER (KIND = iintegers), PARAMETER ::  &
  hmin = 10_iintegers   ! minimal vertical height above ground [m]

LOGICAL :: ldebug


! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin subroutine indpos
!------------------------------------------------------------------------------

  IF (PRESENT(debug)) THEN
    ldebug = debug
  ELSE
    ldebug = .FALSE.
  ENDIF

! Calculate the horizontal position in the grid index space
!------------------------------------------------------------------------------

  x_pos = (lon_pos-minlon)/deltalon - REAL(isubpos(my_cart_id,1)-nboundlines-2,wp)
  y_pos = (lat_pos-minlat)/deltalat - REAL(isubpos(my_cart_id,2)-nboundlines-2,wp)


! Calculate the vertical position in the grid index space
!------------------------------------------------------------------------------

  ! Position relative to the half-levels 
  zh_pos = 0.0_wp

  DO k = 1, ke_tot
    zt = int3d(hhl, ie, je, ke+1, x_pos, y_pos, REAL(k  , wp),0)
    zb = int3d(hhl, ie, je, ke+1, x_pos, y_pos, REAL(k+1, wp),0)
  
    IF ((height_pos <= zt) .AND. (height_pos > zb)) THEN 
      zh_pos = k + (height_pos-zt)/(zb-zt)
      IF (k==ke_tot) THEN
        !zt lowest level; zb surface
        IF (height_pos < zb + hmin) THEN
          zh_pos = ke_tot + 1.0_wp
        ELSE
          zh_pos = ke_tot + (height_pos - zb - hmin)/(zt - zb - hmin)
        ENDIF
      ENDIF
      EXIT
    ENDIF
  ENDDO
  
  ! Position relative to the full-levels 
  z_pos = 0.0_wp
  zt = int3d(hfl, ie, je, ke, x_pos, y_pos, REAL(ke_tot, wp),0)       ! lowest level
  IF (height_pos < zt) THEN
    zb = int3d(hhl, ie, je, ke+1, x_pos, y_pos, REAL(ke_tot+1, wp),0) ! surface
    IF (height_pos < zb + hmin) THEN
      z_pos = ke_tot + 1.0_wp
    ELSE
      z_pos = ke_tot + (height_pos - zb - hmin)/(zt - zb - hmin)
    ENDIF
  ELSE
    DO k = 1, ke_tot-1
      zt = int3d(hfl, ie, je, ke, x_pos, y_pos, REAL(k, wp),0)
      zb = int3d(hfl, ie, je, ke, x_pos, y_pos, REAL(k+1, wp),0)
      IF ((zt >= height_pos) .AND. (zb <= height_pos)) THEN
        z_pos = k + (height_pos - zt)/(zb - zt)
        EXIT
      ELSEIF (k == ke_tot - 1)THEN
        z_pos = ke_tot + 1.0_wp
      ENDIF
    ENDDO
  ENDIF
  
!------------------------------------------------------------------------------
! End of subroutine indpos
!------------------------------------------------------------------------------

END SUBROUTINE indpos

!==============================================================================
!+ interpolation of 3d fields to arbitrary location
!------------------------------------------------------------------------------

REAL (KIND=wp)     FUNCTION int3d(ar,n1,n2,n3,rid,rjd,rkd,istagger)

!------------------------------------------------------------------------------
!
! Description:
!   Interpolates 3d-array to an arbitrary location within the grid.
!
! Method:
!   
!------------------------------------------------------------------------------

! Declarations:

! Function arguments:
! ---------------------
 
INTEGER (KIND=iintegers),INTENT(IN)   ::  &
  n1, n2, n3,istagger     ! size of the 3d-field   

REAL    (KIND=wp),       INTENT(IN)   ::  &
  ar(:,:,:)               ! intput 3d-field

REAL    (KIND=wp),       INTENT(IN)   ::  &
  rid, rjd, rkd           ! location to which the field is interpolated

! Local variables:
! ----------------

INTEGER (KIND=iintegers)              ::  &
  iz,       & ! grid point West  to the point to interpolate
  jz,       & ! grid point South to the point to interpolate
  kz,       & ! grid point below    the point to interpolate
  izp1,     & ! grid point East  to the point to interpolate
  jzp1,     & ! grid point North to the point to interpolate
  kzp1        ! grid point above    the point to interpolate

REAL    (KIND=wp)                     ::  &
! weights for the surrounding grid points in the interpolation
  zfrac0i,  & ! weight for izp1 (grid point East  to the point to interpolate)
  zfrac0j,  & ! weight for jzp1 (grid point North to the point to interpolate)
  zfrac0k,  & ! weight for kzp1 (grid point above    the point to interpolate)
  zfrac1i,  & ! weight for iz   (grid point West  to the point to interpolate)
  zfrac1j,  & ! weight for jz   (grid point South to the point to interpolate)
  zfrac1k     ! weight for kz   (grid point below    the point to interpolate)

REAL    (KIND=wp)                     ::  &
  zri, zrj, zrk ! location to which the field is interpolated (with domain check)

REAL    (KIND=wp)                     ::  &
  zar_ip(2,2,2) ! value of the 3d-field at the 8 surrounding points


INTEGER (KIND=iintegers)              ::  &
  iz1, jz1
REAL    (KIND=wp)                     ::  &
  zex_fac


! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin function int3d
!------------------------------------------------------------------------------
 
  ! Secure that the point to which we should interpolate is within the 
  ! domain and handle staggering

  IF (istagger == 1_iintegers) THEN
    zri = MAX(1.0_wp,MIN(REAL(n1,wp),rid-0.5_wp))
  ELSE
    zri = MAX(1.0_wp,MIN(REAL(n1,wp),rid))
  ENDIF

  IF (istagger == 2_iintegers) THEN
    zrj = MAX(1.0_wp,MIN(REAL(n2,wp),rjd-0.5_wp))
  ELSE
    zrj = MAX(1.0_wp,MIN(REAL(n2,wp),rjd))
  ENDIF

  IF ((rkd > REAL((ke_tot),wp)).AND.(rkd <= REAL((ke_tot+1),wp))) THEN
    zrk = MAX(1.0_wp,REAL(rkd,wp))
  ELSE
    zrk = MAX(1.0_wp,MIN(REAL(n3,wp),rkd))
  ENDIF

  ! Check for interpolation in i
  IF (ABS(REAL(NINT(zri),wp)-zri) < 1.E-3_wp) THEN
    iz   = NINT(zri)
    izp1 = NINT(zri)
  ELSE
    iz   = MAX(MIN(INT(zri,KIND=iintegers),n1-1),1_iintegers)
    izp1 = iz+1
  ENDIF
  
  ! Check for interpolation in j
  IF (ABS(REAL(NINT(zrj),wp)-zrj) < 1.E-3_wp) THEN
    jz   = NINT(zrj)
    jzp1 = NINT(zrj)
  ELSE
    jz    = MAX(MIN(INT(zrj,KIND=iintegers),n2-1),1_iintegers)
    jzp1  = jz+1
  ENDIF
  
  ! Check for interpolation in k
  ! Needed because a trajectory can go below ke_tot even
  ! for a field that is not staggered vertically
  IF (ABS(REAL(NINT(zrk),wp)-zrk) < 1.E-3_wp) THEN
    kz   = NINT(zrk)
    kzp1 = NINT(zrk)
  ELSE
    kz   = MAX(MIN(INT(zrk,KIND=iintegers),ke_tot),1_iintegers)
    kzp1 = kz+1
  ENDIF

  ! Select necessary grid-points from ar and if necessary do extrapolation to surface
  IF (rkd > REAL(ke_tot,wp)) THEN
    IF ((zri-REAL(iz,wp)) >= 0.5_wp) THEN
      iz1 = iz+1
    ELSE
      iz1 = iz
    ENDIF
    IF ((zrj-REAL(jz,wp)) >= 0.5_wp) THEN
      jz1 = jz+1
    ELSE
      jz1 = jz
    ENDIF
     
    IF (istagger == 1_iintegers) THEN

      !Extrapolate field to surface
      zex_fac = (hhl(iz1,jz,ke_tot+1)-hhl(iz1,jz,ke_tot))/                         &
                (hhl(iz1,jz,ke_tot-1)-hhl(iz1,jz,ke_tot+1))
      zar_ip(1,1,1) = 0.5_wp*(ar(iz1-1,jz,ke)+ar(iz1,jz,ke)+                   &
                      zex_fac*(ar(iz1-1,jz,ke-1)-ar(iz1-1,jz,ke)+ar(iz1,jz,ke-1)-ar(iz1,jz,ke)))
        
      zex_fac = (hhl(iz1,jz+1,ke_tot+1)-hhl(iz1,jz+1,ke_tot))/                     &
                (hhl(iz1,jz1+1,ke_tot-1)-hhl(iz1,jz1+1,ke_tot+1))
      zar_ip(1,2,1) = 0.5_wp*(ar(iz1-1,jzp1,ke)+ar(iz1,jzp1,ke)+             &
                      zex_fac*(ar(iz1-1,jzp1,ke-1)-ar(iz1-1,jzp1,ke)+ar(iz1,jzp1,ke-1)-ar(iz1,jzp1,ke)))
        
      zex_fac = (hhl(iz1+1,jz,ke_tot+1)-hhl(iz1+1,jz,ke_tot))/                     &
                (hhl(iz1+1,jz,ke_tot-1)-hhl(iz1+1,jz,ke_tot+1))
      zar_ip(2,1,1) = 0.5_wp*(ar(iz1,jz,ke)+ar(iz1+1,jz,ke)+                   &
                      zex_fac*(ar(iz1,jz,ke-1)-ar(iz1,jz,ke)+ar(iz1+1-1,jz,ke-1)-ar(iz1+1-1,jz,ke)))
        
      zex_fac = (hhl(iz1+1,jz+1,ke_tot+1)-hhl(iz1+1,jz+1,ke_tot))/                 &
                (hhl(iz1+1,jz+1,ke_tot-1)-hhl(iz1+1,jz+1,ke_tot+1))
      zar_ip(2,2,1) = 0.5_wp*(ar(iz1,jzp1,ke)+ar(iz1+1,jzp1,ke)+             &
                      zex_fac*(ar(iz1,jzp1,ke-1)-ar(iz1,jzp1,ke)+ar(iz1+1,jzp1,ke-1)-ar(iz1+1,jzp1,ke)))

      zar_ip(1,1,2) = 0.5_wp*(ar(iz1-1,jz,ke_tot)+ar(iz1,jz,ke_tot)) 
      zar_ip(1,2,2) = 0.5_wp*(ar(iz1-1,jzp1,ke_tot)+ar(iz1,jzp1,ke_tot))
      zar_ip(2,1,2) = 0.5_wp*(ar(iz1,jz,ke_tot)+ar(iz1+1,jz,ke_tot))
      zar_ip(2,2,2) = 0.5_wp*(ar(iz1,jzp1,ke_tot)+ar(iz1+1,jzp1,ke_tot))

      zri = zri+0.5_wp
      iz = iz1

    ELSEIF (istagger == 2_iintegers) THEN

      !Extrapolate field to surface
      zex_fac = (hhl(iz,jz1,ke_tot+1)-hhl(iz,jz1,ke_tot))/                         &
                (hhl(iz,jz1,ke_tot-1)-hhl(iz,jz1,ke_tot+1))
      zar_ip(1,1,1) = 0.5_wp*(ar(iz,jz1-1,ke)+ar(iz,jz1,ke)+                   &
                      zex_fac*(ar(iz,jz1-1,ke-1)-ar(iz,jz1-1,ke)+ar(iz,jz1,ke-1)-ar(iz,jz1,ke)))
        
      zex_fac = (hhl(iz,jz1+1,ke_tot+1)-hhl(iz,jz1+1,ke_tot))/                     &
                (hhl(iz,jz1+1,ke_tot-1)-hhl(iz,jz1+1,ke_tot+1))
      zar_ip(1,2,1) = 0.5_wp*(ar(iz,jz1,ke)+ar(iz,jz1+1,ke)+                   &
                      zex_fac*(ar(iz,jz1,ke-1)-ar(iz,jz1+1,ke)+ar(iz,jz1,ke-1)-ar(iz,jz1+1,ke)))
        
      zex_fac = (hhl(iz+1,jz1,ke_tot+1)-hhl(iz+1,jz1,ke_tot))/                     &
                (hhl(iz+1,jz1,ke_tot-1)-hhl(iz+1,jz1,ke_tot+1))
      zar_ip(2,1,1) = 0.5_wp*(ar(izp1,jz1-1,ke)+ar(izp1,jz1,ke)+             &
                      zex_fac*(ar(izp1,jz1-1,ke-1)-ar(izp1,jz1-1,ke)+ar(izp1,jz1,ke-1)-ar(izp1,jz1,ke)))
        
      zex_fac = (hhl(iz+1,jz1+1,ke_tot+1)-hhl(iz+1,jz1+1,ke_tot))/                 &
                (hhl(iz+1,jz1+1,ke_tot-1)-hhl(iz+1,jz1+1,ke_tot+1))
      zar_ip(2,2,1) = 0.5_wp*(ar(izp1,jz1,ke)+ar(izp1,jz1+1,ke)+             &
                      zex_fac*(ar(izp1,jz1,ke-1)-ar(izp1,jz1,ke)+ar(izp1,jz1+1,ke-1)-ar(izp1,jz1+1,ke)))

      zar_ip(1,1,2) = 0.5_wp*(ar(iz,jz1-1,ke_tot)+ar(iz,jz1,ke_tot)) 
      zar_ip(1,2,2) = 0.5_wp*(ar(iz,jz1,ke_tot)+ar(iz,jz1+1,ke_tot))
      zar_ip(2,1,2) = 0.5_wp*(ar(izp1,jz1-1,ke_tot)+ar(izp1,jz1,ke_tot))
      zar_ip(2,2,2) = 0.5_wp*(ar(izp1,jz1,ke_tot)+ar(izp1,jz1+1,ke_tot))

      zrj = zrj + 0.5_wp
      jz  = jz1

    ELSE
      ! extrapolate field to the surface
      zex_fac = (hhl(iz,jz,ke_tot+1)-hhl(iz,jz,ke_tot))/                          &
                (hhl(iz,jz,ke_tot-1)-hhl(iz,jz,ke_tot+1))
      zar_ip(1,1,1) = ar(iz,jz,ke)+zex_fac*(ar(iz,jz,ke-1)-ar(iz,jz,ke))
        
      zex_fac = (hhl(iz,jz+1,ke_tot+1)-hhl(iz,jz+1,ke_tot))/                      &
                (hhl(iz,jz+1,ke_tot-1)-hhl(iz,jz+1,ke_tot+1))
      zar_ip(1,2,1) = ar(iz,jzp1,ke)+zex_fac*(ar(iz,jzp1,ke-1)-ar(iz,jzp1,ke))
       
      zex_fac = (hhl(iz+1,jz,ke_tot+1)-hhl(iz+1,jz,ke_tot))/                      &
                (hhl(iz+1,jz,ke_tot-1)-hhl(iz+1,jz,ke_tot+1))
      zar_ip(2,1,1) = ar(izp1,jz,ke)+zex_fac*(ar(izp1,jz,ke-1)-ar(izp1,jz,ke))
        
      zex_fac = (hhl(iz+1,jz+1,ke_tot+1)-hhl(iz+1,jz+1,ke_tot))/                  &
                (hhl(iz1+1,jz1+1,ke_tot-1)-hhl(iz1+1,jz1+1,ke_tot+1))
      zar_ip(2,2,1) = ar(izp1,jzp1,ke)+zex_fac*(ar(izp1,jzp1,ke-1)-ar(izp1,jzp1,ke))
        
      zar_ip(1,1,2) = ar(iz,jz,ke_tot) 
      zar_ip(1,2,2) = ar(iz,jzp1,ke_tot)
      zar_ip(2,1,2) = ar(izp1,jz,ke_tot)
      zar_ip(2,2,2) = ar(izp1,jzp1,ke_tot)
    ENDIF
  ELSE
    zar_ip(1,1,1) = ar(iz,jz,kz) 
    zar_ip(1,2,1) = ar(iz,jzp1,kz)
    zar_ip(2,1,1) = ar(izp1,jz,kz)
    zar_ip(2,2,1) = ar(izp1,jzp1,kz)
    zar_ip(1,1,2) = ar(iz,jz,kzp1) 
    zar_ip(1,2,2) = ar(iz,jzp1,kzp1)
    zar_ip(2,1,2) = ar(izp1,jz,kzp1)
    zar_ip(2,2,2) = ar(izp1,jzp1,kzp1)
  ENDIF
  
  ! Do interpolation.
  IF(kz == kzp1) THEN
    IF (l2dim) THEN
      zfrac0i = zri-REAL(iz,wp)
      zfrac1i = 1._wp-zfrac0i
      int3d = zar_ip(1,1,1)*zfrac1i + zar_ip(2,1,1)*zfrac0i
    ELSEIF((iz == izp1) .AND. (jz == jzp1)) THEN
      int3d = zar_ip(1,1,1)
    ELSE
      zfrac0i = zri-REAL(iz,wp)
      zfrac0j = zrj-REAL(jz,wp)
      zfrac1i = 1._wp-zfrac0i
      zfrac1j = 1._wp-zfrac0j
      int3d = zar_ip(1,1,1)*zfrac1i*zfrac1j + zar_ip(1,2,1)*zfrac1i*zfrac0j + &
              zar_ip(2,1,1)*zfrac0i*zfrac1j + zar_ip(2,2,1)*zfrac0i*zfrac0j
    ENDIF
  ELSE
    zfrac0k = zrk-REAL(kz,wp) 
    zfrac1k = 1._wp-zfrac0k
    IF(l2dim) THEN
      zfrac0i = zri-REAL(iz,wp)
      zfrac1i = 1._wp-zfrac0i
      int3d = zar_ip(1,1,1)*zfrac1i*zfrac1k +                                 &
              zar_ip(2,1,1)*zfrac0i*zfrac1k +                                 &
              zar_ip(1,1,2)*zfrac1i*zfrac0k +                                 &
              zar_ip(2,1,2)*zfrac0i*zfrac0k
    ELSEIF ((iz == izp1) .AND. (jz == jzp1)) THEN
      int3d = zar_ip(1,1,1)*zfrac1k + zar_ip(1,1,2)*zfrac0k
    ELSE
      zfrac0i = zri-REAL(iz,wp)
      zfrac0j = zrj-REAL(jz,wp)
      zfrac1i = 1._wp-zfrac0i
      zfrac1j = 1._wp-zfrac0j
      int3d = zar_ip(1,1,1)*zfrac1i*zfrac1j*zfrac1k +                         &
              zar_ip(1,2,1)*zfrac1i*zfrac0j*zfrac1k +                         &
              zar_ip(2,1,1)*zfrac0i*zfrac1j*zfrac1k +                         &
              zar_ip(2,2,1)*zfrac0i*zfrac0j*zfrac1k +                         &
              zar_ip(1,1,2)*zfrac1i*zfrac1j*zfrac0k +                         &
              zar_ip(1,2,2)*zfrac1i*zfrac0j*zfrac0k +                         &
              zar_ip(2,1,2)*zfrac0i*zfrac1j*zfrac0k +                         &
              zar_ip(2,2,2)*zfrac0i*zfrac0j*zfrac0k
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! End of function int3d
!------------------------------------------------------------------------------

END FUNCTION int3d

!==============================================================================
!+ interpolation of 2d fields to arbitrary location
!------------------------------------------------------------------------------

REAL (KIND=wp)     FUNCTION int2d(ar,n1,n2,rid,rjd,rkd)

!------------------------------------------------------------------------------
!
! Description:
!   Interpolates 2d-array to an arbitrary location within the grid.
!
! Method:
!   
!------------------------------------------------------------------------------

! Declarations:

! Function arguments:
! ---------------------

INTEGER (KIND = iintegers),INTENT(IN)        ::  &
  n1, n2             ! size of 2d-field

REAL (KIND = wp),    INTENT(IN)              ::  &
  ar(:,:,:)          ! input 2d-field

REAL (KIND=wp),    INTENT(IN)                              ::  &
  rid, rjd, rkd      ! location to which is interpolated

! Local variables:
! ----------------

INTEGER (KIND=iintegers)              ::  &
  iz,       & ! grid point West  to the point to interpolate
  jz,       & ! grid point South to the point to interpolate
  kz,       & ! point to which is interpolated in the vertical
  izp1,     & ! grid point East  to the point to interpolate
  jzp1        ! grid point North to the point to interpolate

REAL    (KIND=wp)                     ::  &
! weights for the surrounding grid points in the interpolation
  zfrac0i,  & ! weight for izp1 (grid point East  to the point to interpolate)
  zfrac0j,  & ! weight for jzp1 (grid point North to the point to interpolate)
  zfrac1i,  & ! weight for iz   (grid point West  to the point to interpolate)
  zfrac1j     ! weight for jz   (grid point South to the point to interpolate)

REAL    (KIND=wp)                     ::  &
  zri, zrj    ! location to which the field is interpolated (with domain check)



! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin function int2d
!------------------------------------------------------------------------------

  zri = MAX(1._wp, MIN(REAL(n1,wp),rid))
  zrj = MAX(1._wp, MIN(REAL(n2,wp),rjd))

  ! Check for interpolation in i
  IF (ABS(REAL(NINT(zri),wp)-zri) < 1.E-3_wp) THEN
    iz   = NINT(zri)
    izp1 = NINT(zri)
  ELSE
    iz   = MAX(MIN(INT(zri,KIND=iintegers),n1-1_iintegers), 1_iintegers)
    izp1 = iz + 1_iintegers
  ENDIF
  
  ! Check for interpolation in j
  IF (ABS(REAL(NINT(zrj),wp)-zrj) < 1.E-3_wp) THEN
    jz   = NINT(zrj)
    jzp1 = NINT(zrj)
  ELSE
    jz   = MAX(MIN(INT(zrj,KIND=iintegers),n2-1_iintegers), 1_iintegers)
    jzp1 = jz + 1_iintegers
  ENDIF

  kz = INT(rkd,KIND=iintegers)
  
  ! Do interpolation.
  IF (l2dim) THEN
    zfrac0i = zri-REAL(iz,wp)
    zfrac1i = 1._wp-zfrac0i
    int2d = ar(iz,jz,kz)*zfrac1i + ar(izp1,jz,kz)*zfrac0i
  ELSEIF ((iz == izp1) .AND. (jz == jzp1)) THEN
    int2d = ar(iz,jz,kz)
  ELSE
    zfrac0i = zri-REAL(iz,wp)
    zfrac0j = zrj-REAL(jz,wp)
    zfrac1i = 1._wp-zfrac0i
    zfrac1j = 1._wp-zfrac0j
    int2d = hhl(iz,jz,kz)  *zfrac1i*zfrac1j + hhl(iz,jzp1,kz)  *zfrac1i*zfrac0j + &
            hhl(izp1,jz,kz)*zfrac0i*zfrac1j + hhl(izp1,jzp1,kz)*zfrac0i*zfrac0j
  ENDIF
  
!------------------------------------------------------------------------------
! End of function intp2d
!------------------------------------------------------------------------------

END FUNCTION int2d

!==============================================================================
!+ determine on which PE the trajectory is
!------------------------------------------------------------------------------

PURE FUNCTION ipe_location(pos_x, pos_y, startlon_tot, startlat_tot,          &
                           isubpos_x, isubpos_y, dlon, dlat, idxpe_n,         &
                           idxpe_ne, idxpe_e, idxpe_se, idxpe_s, idxpe_sw,    &
                           idxpe_w, idxpe_nw, idxpe_this, nboundlines,        &
                           ie, je)

!------------------------------------------------------------------------------
!
! Description:
!   Determine if the trajectory is still on the PE or if it is on one
!   of the 8 neighbor PEs.
!
! Method:
!   Compare the trajectory position against the boundaries coordinatesa 
!   of the PE.
!   
!------------------------------------------------------------------------------

! Declarations:

! Function arguments:
! ---------------------

INTEGER (KIND=iintegers) :: ipe_location ! processor index
 
INTEGER (KIND=iintegers),INTENT(IN)   :: & 
  isubpos_x,     & ! x-position of this subdomain in the total domain
  isubpos_y,     & ! y-position of this subdomain in the total domain
  idxpe_this,    & ! index for this                    PE
  idxpe_n   ,    & ! index for the North      neighbor PE 
  idxpe_ne  ,    & ! index for the North-East neighbor PE 
  idxpe_e   ,    & ! index for the East       neighbor PE 
  idxpe_se  ,    & ! index for the South-East neighbor PE 
  idxpe_s   ,    & ! index for the South      neighbor PE 
  idxpe_sw  ,    & ! index for the South-West neighbor PE 
  idxpe_w   ,    & ! index for the West       neighbor PE 
  idxpe_nw  ,    & ! index for the North-West neighbor PE 
  nboundlines,   & ! number of overlapping boundary lines of the subdomains
  ie,            & ! number of grid points in zonal direction
  je               ! number of grid points in meridional direction

REAL (KIND=wp)          ,INTENT(IN)   ::  &
  pos_x     ,    & ! trajectory position in zonal direction
  pos_y     ,    & ! trajectory position in merid. direction
  startlon_tot,  & ! longitude of lower left corner in tot. domain
  startlat_tot,  & ! latitude  of lower left corner in tot. domain
  dlon,          & ! grid spacing in x in degree
  dlat             ! grid spacing in y in degree


! Local variables:
! ----------------


! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin function ipe_location
!------------------------------------------------------------------------------

  IF (pos_x < startlon_tot+(isubpos_x-1)*dlon) THEN
    IF (pos_y < startlat_tot+(isubpos_y-1)*dlat) THEN
      ipe_location = idxpe_sw
    ELSEIF (pos_y >= startlat_tot+(je+isubpos_y-1-2*nboundlines)*dlat) THEN
      ipe_location = idxpe_nw
    ELSE
      ipe_location = idxpe_w
    ENDIF
  ELSEIF (pos_x >= startlon_tot+(ie+isubpos_x-1-2*nboundlines)*dlon) THEN
    IF (pos_y < startlat_tot+(isubpos_y-1)*dlat) THEN
      ipe_location = idxpe_se
    ELSEIF (pos_y >= startlat_tot+(je+isubpos_y-1-2*nboundlines)*dlat) THEN
      ipe_location = idxpe_ne
    ELSE
      ipe_location = idxpe_e
    ENDIF
  ELSEIF ((pos_x >= startlon_tot+(isubpos_x-1)*dlon) .AND.        &
          (pos_x <  startlon_tot+(ie+isubpos_x-1-2*nboundlines)*dlon)) THEN
    IF (pos_y <  startlat_tot+(isubpos_y-1)*dlat) THEN
      ipe_location = idxpe_s
    ELSEIF (pos_y >= startlat_tot+(je+isubpos_y-1-2*nboundlines)*dlat) THEN
      ipe_location = idxpe_n
    ELSE
      ipe_location = idxpe_this
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! End of function ipe_location
!------------------------------------------------------------------------------

END FUNCTION ipe_location

!==============================================================================
!+ determine if trajectory has left the horizontal total domain
!------------------------------------------------------------------------------

PURE FUNCTION lquit_hori_tot_domain(pos_x, pos_y, startlon_tot, startlat_tot, &
                                    dlon, dlat, nboundlines, ie_tot, je_tot)

!------------------------------------------------------------------------------
!
! Description:
!   Determine if the trajectory has left the horizontal total domain.
!
! Method:
!   Compare the trajectory position against the boundaries coordinates
!   of the total domain.
!   
!------------------------------------------------------------------------------

! Declarations:

! Function arguments:
! ---------------------

LOGICAL    :: lquit_hori_tot_domain ! has left total domain?
 
INTEGER (KIND=iintegers),INTENT(IN)   ::  &
  ie_tot,        & ! number of grid points in zonal direction
  je_tot,        & ! number of grid points in meridional direction
  nboundlines      ! number of overlapping boundary lines of the subdomains
 

REAL (KIND=wp)          ,INTENT(IN)   ::  &
  pos_x     ,    & ! trajectory position in zonal direction
  pos_y     ,    & ! trajectory position in merid. direction
  startlon_tot,  & ! longitude of lower left corner in tot. domain
  startlat_tot,  & ! latitude  of lower left corner in tot. domain
  dlon,          & ! grid spacing in x in degree
  dlat             ! grid spacing in y in degree


! Local variables:
! ----------------


! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin function lquit_hori_tot_domain
!------------------------------------------------------------------------------

  IF ((pos_x < startlon_tot + dlon*nboundlines)          .OR.      &
      (pos_x > startlon_tot + dlon*(ie_tot-nboundlines)) .OR.      &
      (pos_y < startlat_tot + dlat*nboundlines)          .OR.      &
      (pos_y > startlat_tot + dlat*(je_tot-nboundlines))) THEN      
    lquit_hori_tot_domain = .TRUE.
  ELSE
    lquit_hori_tot_domain = .FALSE.
  ENDIF 

!------------------------------------------------------------------------------
! End of function lquit_hori_tot_domain
!------------------------------------------------------------------------------

END FUNCTION lquit_hori_tot_domain

!==============================================================================
!+ Module procedure for restart handling
!------------------------------------------------------------------------------

SUBROUTINE organize_traj_restart(yaction)

!------------------------------------------------------------------------------
!
! Description:
!   Read or write an extra restart file for the trajectory data.
!
! Method:
!
!      
!------------------------------------------------------------------------------

! Declarations:

! Subroutine arguments:
! ---------------------

CHARACTER (LEN= *),       INTENT(IN)            ::                            &
  yaction      ! action to be performed

! Local variables:
! ----------------

CHARACTER (LEN= 14)        :: yzdat1
CHARACTER (LEN= 28)        :: yzdat2
CHARACTER (LEN= 80)        :: yzerror, yzroutine
CHARACTER (LEN=250)        :: yfname
INTEGER   (KIND=iintegers) :: nzjulianday, izerror, izerr, nsend, i, ind,    &
                              iz_proc, ntraj_restart
REAL (KIND=wp)             :: zacthour

INTEGER (KIND=iintegers), ALLOCATABLE           ::  &
  displs(:),             & ! displacement in receive buffer
  all_send(:)              ! number of datasets received form 
REAL (KIND=wp)     , ALLOCATABLE                ::  &
  tranew(:),             & ! receiving buffer
  sendbuf(:)               ! send buffer
  
! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin subroutine organize_traj_restart
!------------------------------------------------------------------------------

  izerror   = 0_iintegers
  yzroutine = 'organize_traj_restart'
  yzerror   = '  '

!------------------------------------------------------------------------------
! Section 1: Precomputations
!------------------------------------------------------------------------------

  !Create restart filename for trajectories
  IF (my_cart_id == 0) THEN
    CALL get_utc_date(ntstep, ydate_ini, dt, itype_calendar, yzdat1,yzdat2, nzjulianday, zacthour)
    CALL make_fn ('lrt', yzdat1, ydate_ini, ytunit_restart,'o',ntstep, dt, .TRUE., itype_calendar, &
                   ydir_restart_out, yfname, .FALSE., idbg_level, izerr)
    IF (izerr/=0) THEN
      izerror = 10420_iintegers
      yzerror = 'TRAJECTORY module: Generation of restart file name failed.'
      CALL model_abort(my_cart_id, izerror, yzerror, yzroutine)
    ENDIF
    yfname=TRIM(yfname)
  ENDIF
  
!------------------------------------------------------------------------------
! Section 2: Read of restart
!------------------------------------------------------------------------------

  IF (yaction == 'read') THEN
    !Allocate sendbuffer
    ALLOCATE (sendbuf(ntraj*5),STAT=izerr)
    IF (izerr/=0) THEN
      izerror = 10423_iintegers
      yzerror = 'TRAJECTORY module: allocation of sendbuf failed.'
      CALL model_abort(my_cart_id,izerror,yzerror,yzroutine)
    ENDIF
    sendbuf(:) = 0.0_wp
    IF (my_cart_id == 0) THEN
      IF (idbg_level > 1) THEN
        PRINT *, 'OPEN: bina-file: ', yfname
      ENDIF
!ROA: here the function get_free_unit should be used instead of hard-coded 99 
      OPEN(unit=99,file=yfname,FORM='unformatted',action='READ',access='STREAM',iostat=izerr)
      IF (izerr/=0) THEN
        izerror = 10426_iintegers
        yzerror = 'TRAJECTORY module: Open of restart file failed.'
        CALL model_abort(my_cart_id, izerror, yzerror, yzroutine)
      ENDIF
      !Read number of entries
      READ(99) ntraj_restart
      IF (ntraj_restart > ntraj) THEN
         izerror = 10429_iintegers
         yzerror = 'TRAJECTORY module: ntraj_restart > ntraj.'
         CALL model_abort(my_cart_id, izerror,yzerror,yzroutine)
      ENDIF
      !Read positions
      sendbuf(:) = 0.0_wp
      READ(99) sendbuf(1:5*ntraj_restart)
      IF (idbg_level > 1) THEN
        PRINT *, 'CLOSING bina FILE'
      ENDIF
      CLOSE(99)
!AKMITODO: Error message if closing of file fails
    ENDIF
    !Distribute buffer
    IF (num_compute > 1_iintegers) THEN
      CALL MPI_Bcast(sendbuf,5*ntraj, imp_reals,0,icomm_cart,izerr)
      IF (izerr/=0) THEN
        izerror = 10432_iintegers
        yzerror = 'TRAJECTORY module: MPI_Bcast failed.'
        CALL model_abort(my_cart_id, izerror, yzerror, yzroutine)
      ENDIF
    ENDIF
    !Set all Trajectories to dead
    istat_traj(:) = -999_iintegers
    DO i=1,ntraj
      !Check if end of restart list is reached
      IF (NINT(sendbuf((i-1)*5+1)) == 0) EXIT
      !Check domain and copy
        iz_proc = ipe_location(sendbuf((i-1)*5+2), sendbuf((i-1)*5+3), startlon_tot, &
                                       startlat_tot, isubpos(my_cart_id,1),          &
                                       isubpos(my_cart_id,2), dlon, dlat, -1,        &
                                       -1, -1, -1, -1,                               &
                                       -1, -1, -1, my_cart_id,                       &
                                       nboundlines, ie, je)
      IF (iz_proc==my_cart_id) THEN
        traj(NINT(sendbuf((i-1)*5+1)),:) = sendbuf((i-1)*5+2:(i-1)*5+4)
      ELSE
        !Set some coordinates that are outside the subdomain
        traj(NINT(sendbuf((i-1)*5+1)),1) = -999.0_wp
        traj(NINT(sendbuf((i-1)*5+1)),2) = -999.0_wp
        traj(NINT(sendbuf((i-1)*5+1)),3) =    0.0_wp
      ENDIF
      istat_traj(NINT(sendbuf((i-1)*5+1))) = NINT(sendbuf(i*5))
    ENDDO

    ! Deallocate sendbuffer
    DEALLOCATE (sendbuf, STAT=izerr)
    IF (izerr /= 0_iintegers) THEN
      izerror  = 10435_iintegers
      yzerror = 'Deallocation of sendbuf for reading traj. restart file failed.'
      CALL model_abort(my_cart_id, izerror, yzerror, yzroutine)
    ENDIF

!------------------------------------------------------------------------------
! Section 3: Write of the restart
!------------------------------------------------------------------------------

  ELSEIF (yaction == 'write') THEN
    !Allocate fields for data exchange
    ALLOCATE (all_send(num_compute),STAT=izerr) 
    all_send(:) = 0_iintegers
    ALLOCATE (sendbuf(ntraj*5),STAT=izerr)
    sendbuf(:) = 0.0_wp
    ALLOCATE (displs(num_compute), STAT=izerr)
    displs = 0_iintegers
    ALLOCATE (tranew(ntraj*5), STAT=izerr)
    IF (izerr/=0) THEN
      izerror = 10438_iintegers
      yzerror = 'TRAJECTORY module: Allocation of tranew failed.'
      CALL model_abort(my_cart_id, izerror, yzerror, yzroutine)
    ENDIF
  
    tranew(:) = 0.0_wp
    
    !Fill sendbuffers
    ind = 1_iintegers
    
    DO i = 1,ntraj
      ! Check if trajectory is in my domain
      iz_proc = ipe_location(traj(i,1), traj(i,2), startlon_tot,          &
                            startlat_tot, isubpos(my_cart_id,1),          &
                            isubpos(my_cart_id,2), dlon, dlat, -1,        &
                            -1, -1, -1, -1,                               &
                            -1, -1, -1, my_cart_id,                       &
                            nboundlines, ie, je)
    
      IF(iz_proc == my_cart_id) THEN
        sendbuf(ind)         = REAL(i,wp)
        sendbuf(ind+1:ind+3) = traj(i,:)
        sendbuf(ind+4)       = REAL(istat_traj(i),wp)
        ind = ind+5_iintegers  
      ENDIF
    ENDDO
    
    nsend = ind-1_iintegers
    
    IF (num_compute > 1_iintegers) THEN
    ! Gather number of trajectory datasets to be send from each process to 
    ! process 0
      CALL MPI_GATHER(nsend,1,imp_integers,all_send,1,imp_integers,0,        &
                      icomm_cart,izerr)
      IF (izerr /= 0_iintegers) THEN
        izerror = 10441_iintegers
        yzerror = 'TRAJECTORY module: GATHER failed while writing restart file.'
        CALL model_abort(my_cart_id, izerror, yzerror, yzroutine)
      ENDIF
    ELSE
      all_send(1) = nsend
    ENDIF
 
    ! Prepare variables for GATHERV
    displs(1) = 0_iintegers
    DO i = 1,(num_compute-1)
      displs(i+1) = SUM(all_send(1:i))
    ENDDO
    IF (num_compute > 1_iintegers) THEN
      ! Gather all trajectory datasets to process 0
      CALL MPI_GATHERV(sendbuf(1:nsend),nsend,imp_reals,tranew,all_send,displs,&
                       imp_reals,0,icomm_cart,izerror)
      IF (izerr /= 0_iintegers) THEN
        izerror = 10444_iintegers
        yzerror = 'TRAJECTORY module: GATHER failed while writing restart file.'
        CALL model_abort(my_cart_id, izerror, yzerror, yzroutine)
      ENDIF
    ELSE
      tranew(1:nsend) = sendbuf(1:nsend)
    ENDIF

    !Open and write trajectory restart file
    IF (my_cart_id == 0) THEN
      !Determine how many trajectories have to be written into the restart file
      !Dead trajectories don't have to be saved.
      ntraj_restart=ntraj
      DO i=1,ntraj
        IF (tranew(1+5*(i-1)) == 0.0_wp) THEN
          ntraj_restart=i-1
          EXIT
        ENDIF
      ENDDO
      IF (idbg_level > 1) THEN
        PRINT *, 'OPEN: bina-file: ', yfname
      ENDIF
      OPEN(unit=99,file=yfname,FORM='unformatted',access='stream',iostat=izerr)
      IF (izerr /= 0) THEN
        izerror = 10447_iintegers
        yzerror = 'TRAJECTORY module: Error while opening restart file.'
        CALL model_abort(my_cart_id, izerror, yzerror, yzroutine)
      ENDIF
      !Write ntraj_restart
      WRITE(99,iostat=izerr) ntraj_restart
      !Write position and status of trajectories
      WRITE(99,iostat=izerr) tranew(1:5*ntraj_restart)
      IF (izerr /= 0) THEN
        izerror = 10450_iintegers
        yzerror = 'TRAJECTORY module: Error while writing restart file.'
        CALL model_abort(my_cart_id, izerror, yzerror, yzroutine)
      ENDIF
      !Close file
      IF (idbg_level > 1) THEN
        PRINT *, 'CLOSING bina FILE'
      ENDIF
      CLOSE(99)
    ENDIF
  
    !Deallocate space
    DEALLOCATE(all_send, STAT=izerr)
    DEALLOCATE(sendbuf,STAT=izerr)
    DEALLOCATE(displs, STAT=izerr)
    DEALLOCATE(tranew, STAT=izerr)

  ENDIF

!------------------------------------------------------------------------------
! End of subroutine organize_traj_restart
!------------------------------------------------------------------------------ 
 
END SUBROUTINE organize_traj_restart

!==============================================================================
!==============================================================================

END MODULE src_traj
