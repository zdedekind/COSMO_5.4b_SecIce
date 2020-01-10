!+ Data module for variables of the online trajectory calculation
!------------------------------------------------------------------------------

MODULE data_traj

!------------------------------------------------------------------------------
!
! Description:
!  This module contains variables and fields for running online trajectory 
!  calculations. The primary control switch ltraj is defined in 
!  data_runcontrol.f90. 
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
!
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
  iintegers, & ! KIND-type parameter for standard integer variables
  intgribf     ! KIND-type parameter for fortran files in the grib library

USE data_io, ONLY : &
  clen

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Parameters 
! -----------

INTEGER(KIND=iintegers), PARAMETER :: &
  nmax_traj_starttime  =   100_iintegers,   & ! maximum number of traject. start times  [-]
  nmax_traj_startpt    = 50000_iintegers,   & ! maximum number of traject. start points [-]
  nmax_traced_var      =    30_iintegers,   & ! maximum number of traced variables      [-] 
  nmax_out_per_start_t =   100_iintegers,   & ! maximum number of output file per       [-]
                                              ! start time
  max_nc_size          = 2**30-1+2**30_iintegers ,& ! maximum size of a trajectory NetCDF output file (about2GB)
  header_nc_size       = 2*2**20_iintegers    ! reserved space for NetCDF header (2MB)


INTEGER(KIND=iintegers), PARAMETER :: &
  idxpe_this = 0_iintegers,    & ! index for this                    PE
  idxpe_n    = 1_iintegers,    & ! index for the North      neighbor PE 
  idxpe_ne   = 2_iintegers,    & ! index for the North-East neighbor PE 
  idxpe_e    = 3_iintegers,    & ! index for the East       neighbor PE 
  idxpe_se   = 4_iintegers,    & ! index for the South-East neighbor PE 
  idxpe_s    = 5_iintegers,    & ! index for the South      neighbor PE 
  idxpe_sw   = 6_iintegers,    & ! index for the South-West neighbor PE 
  idxpe_w    = 7_iintegers,    & ! index for the West       neighbor PE 
  idxpe_nw   = 8_iintegers       ! index for the North-West neighbor PE 

INTEGER(KIND=iintegers), PARAMETER :: &
  idead      = -999_iintegers, & ! status of a trajectory that is not on this PE
  ialive     =   -1_iintegers    ! status of a trajectory that is running on this PE

! Namelist variables
!--------------------

INTEGER (KIND=iintegers)                 :: &
  nstop_traj,            & ! last time step of the trajectory computation  [-]
  ninc_out_traj,         & ! time interval for trajectory output in steps  [-]
  nstart_traj(nmax_traj_starttime), & ! start time steps of the traject.   [-]
  ncomb_start_traj(3),   & ! triplet(start,stop,dt) for start time steps   [-]
  istart_mode_traj         ! 1: start trajectories at user-specified times [-]
                           ! 2: start trajectories at predefined intervals

REAL      (KIND=wp)                       :: &
  hstop_traj,            & ! end of trajectory computation                 [h]
  hinc_out_traj,         & ! time interval for trajectory output           [h]
  hstart_traj(nmax_traj_starttime), & ! start times of the traject.        [h]
  hcomb_start_traj(3)      ! triplet(start,stop,dt) for start time         [h]

CHARACTER (LEN=100)                       :: &
  ydir_traj,             & ! output directory for trajectories             [-]
  startfile_traj           ! path to startfile                             [-]

CHARACTER (LEN=15)                        :: &
  tracevar_traj(nmax_traced_var)  ! variables traced along traject.   [depend on varia]


! Variables used in the trajectory module
!----------------------------------------

INTEGER (KIND=iintegers)                 :: &
  num_start_t                        ,& ! number of starting times for a same
                                        ! spatial starting position                      [-]
  num_start_t_tot                    ,& ! total number of starting times                 [-]
  num_start_pt(nmax_traj_starttime)  ,& ! number of traj. start points per starting time [-]
  num_start_pt_tot                      ! total number of start points                   [-]

INTEGER (KIND=iintegers)                 :: &
  ntraj                ! number of trajectories                            [-]

INTEGER (KIND=iintegers)                 :: &
  nstart_traj_gp(nmax_traj_starttime),& ! traj. group starttimes  (in model t step)[-]
  nstart_traj_iv(nmax_traj_startpt),  & ! indiv. traj. starttimes (in model t step)[-]
  nstartfirst_traj                      ! time step to start first trajectories    [-]

REAL (KIND=wp),     ALLOCATABLE, TARGET           :: &
  hfl  (:,:,:),    &   ! geometrical height of full levels                 [m]
  traj(:,:)            ! trajectory spatial information    
                       !     dim 1:  # of startpoints * # of init times (== # trajectories) 
                       !     dim 2:  # of spatial coordinates (3)                           
                       ! traj(:,1)/traj(:,2) are rotated lon/lat           [deg]
                       ! traj(:,3) is height above mean sea level          [m]


INTEGER (KIND=iintegers), ALLOCATABLE     :: &
  istat_traj(:)        ! indicates status of a trajectory              [-]
                       ! -1 : trajectory has started and not yet left 
                       !      the domain, should be computed ("alive")
                       !>=0 : trajectory has not yet started and will start
                       !      at ntstep=istat_traj           ("embryonic")
                       !-999: trajectory has left the domain ("dead")
 


TYPE trajtrace_type
  CHARACTER(LEN=clen)                  :: name
  INTEGER  (KIND=intgribf)             :: levtyp
  INTEGER  (KIND=intgribf)             :: ntri
  INTEGER  (KIND=intgribf)             :: rank
  REAL     (KIND=wp),        POINTER   :: p4   (:,:,:,:) => NULL()
  REAL     (KIND=wp),        POINTER   :: p3   (:,:,:)   => NULL()
  INTEGER  (KIND=iintegers)            :: idef_stat
  INTEGER  (KIND=iintegers)            :: istag
  CHARACTER (LEN=20)                   :: units
  CHARACTER (LEN=80)                   :: sdname
  CHARACTER (LEN=80)                   :: lgname
END TYPE trajtrace_type


TYPE(trajtrace_type), ALLOCATABLE :: traj_trace(:)



!==============================================================================

END MODULE data_traj
