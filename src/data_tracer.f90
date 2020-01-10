!+ Data structure for general treatment of tracers
!------------------------------------------------------------------------------

MODULE data_tracer

!------------------------------------------------------------------------------
!
! Description:
!  This module declares parameters for all tracers that have to reside in
!  the long term storage, i.e. that are used in more than one module.
!  Fields included are
!    - moisture related tracers (microphysics)
!    - chemical species (active tracers)
!    - pollens (quasi passive tracers)
!    - aerosols (active tracers)
!
!  All tracers are integrated in a data structure. They are allocated in the
!  setup of the model and deallocated in the cleanup at the end of the
!  program.
!
! Current Code Owner: MeteoSwiss, Oliver Fuhrer
!  phone:  +41 58 460 9359
!  fax:    +41 58 460 9278
!  email:  oliver.fuhrer@meteoswiss.ch
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer 
!  Initial release
! V4_26        2012/12/06 Anne Roches
!  Changes and technical adaptations to the tracer handling
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  Introduced MESSy interface
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Adaption to MESSy tracer interface
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
!
! Declarations:
!
! Modules used:
USE data_parameters,       ONLY:                                             &
       wp,        & ! KIND-type parameter for real variables
       iintegers    ! KIND-type parameter for standard integer variables

USE data_io,               ONLY:                                             &
       clen         ! length of character string for variable names

USE data_tracer_metadata,  ONLY:                                             &
       t_metadata   ! structure for metadata storage

#ifdef MESSY
USE messy_main_tracer, ONLY: I_GRIBPARAM, I_GRIBTAB, I_ADVECT,  I_HDIFF &
                           , I_VDIFF,     I_CONVECT, I_INITIAL, I_LBC   &
                           , I_RELAX,     I_DAMP
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================


! 1. Global parameters, dimensions
! ------------------------------------------------

INTEGER (KIND=iintegers), PARAMETER  :: &
  n_trcr_max    = 300_iintegers,    & ! maximal number of tracers
  n_mkey_max    = 50_iintegers,     & ! maximal number of metadata entries
  n_mbuf_max    = 1024_iintegers      ! maximal buffer size for metadata storage (Bytes)

INTEGER (KIND=iintegers), PARAMETER  :: &
  ilen_sn       = clen,             & ! length of name (character)
  ilen_ln       = 80_iintegers,     & ! length of parent, standard name and long name 
  ilen_en       = 255_iintegers,    & ! length of error string
  ilen_un       = 20_iintegers        ! length of unit name and level type


! 2. Enumerate parameters
! ------------------------------------------------

! tracer module status
INTEGER (KIND=iintegers), PARAMETER  :: &
  T_STAT_START   = 0_iintegers,    &  ! before init at start of executation
  T_STAT_DEFINE  = 1_iintegers,    &  ! after init, while defining tracers
  T_STAT_ALLOC   = 2_iintegers,    &  ! after memory allocation
  T_STAT_FINISH  = 3_iintegers        ! after memory deallocation

! missing index
INTEGER (KIND=iintegers), PARAMETER  :: &
  T_MISSING      = -999_iintegers     ! index value for missing locations

! error for tracer not found
INTEGER (KIND=iintegers), PARAMETER  :: &
  T_ERR_NOTFOUND = -888_iintegers     ! error value for not found error


! 3. Standard tracer meta-data
! ------------------------------------------------

! name of tracer (string)
CHARACTER (LEN=*), PARAMETER ::         &
  T_NAME_KEY = 'NAME'
INTEGER (KIND=iintegers) ::             &
  T_NAME_ID = -1
CHARACTER (LEN=*), PARAMETER ::         &
  T_NAME_DEFAULT = 'undefined'

! parent module of tracer (string)
CHARACTER (LEN=*), PARAMETER ::         &
  T_PARENT_KEY = 'PARENT'
INTEGER (KIND=iintegers) ::             &
  T_PARENT_ID = -2
CHARACTER (LEN=*), PARAMETER ::         &
  T_PARENT_DEFAULT = 'undefined'

! units of tracer (string)
CHARACTER (LEN=*), PARAMETER ::         &
  T_UNITS_KEY = 'UNITS'
INTEGER (KIND=iintegers) ::             &
  T_UNITS_ID = -3
CHARACTER (LEN=*), PARAMETER ::         &
  T_UNITS_DEFAULT = 'undefined'

! NetCDF standard name of tracer (string)
CHARACTER (LEN=*), PARAMETER ::         &
  T_NCSTDNAME_KEY = 'NC_STDNAME'
INTEGER (KIND=iintegers) ::             &
  T_NCSTDNAME_ID = -4
CHARACTER (LEN=*), PARAMETER ::         &
  T_NCSTDNAME_DEFAULT = 'undefined'

! NetCDF long name of tracer (string)
CHARACTER (LEN=*), PARAMETER ::         &
  T_NCLONGNAME_KEY = 'NC_LONGNAME'
INTEGER (KIND=iintegers) ::             &
  T_NCLONGNAME_ID = -5
CHARACTER (LEN=*), PARAMETER ::         &
  T_NCLONGNAME_DEFAULT = 'undefined'

! grib parameter number of tracer (integer)
CHARACTER (LEN=*), PARAMETER ::         &
  T_GRBPARAM_KEY = 'GRB_PARAM'
#ifndef MESSY
INTEGER (KIND=iintegers) ::             &
  T_GRBPARAM_ID = -6
#else
  INTEGER (KIND=iintegers), PARAMETER ::  T_GRBPARAM_ID = I_GRIBPARAM
#endif
INTEGER (KIND=iintegers) ::             &
  T_GRBPARAM_DEFAULT = -999_iintegers

! grib parameter number of tracer (integer)
CHARACTER (LEN=*), PARAMETER ::         &
  T_GRBTABLE_KEY = 'GRB_TABLE'
#ifndef MESSY
INTEGER (KIND=iintegers) ::             &
  T_GRBTABLE_ID = -7
#else
INTEGER (KIND=iintegers), PARAMETER :: T_GRBTABLE_ID = I_GRIBTAB
#endif
INTEGER (KIND=iintegers) ::             &
  T_GRBTABLE_DEFAULT = -999_iintegers

! advection (integer)
CHARACTER (LEN=*), PARAMETER ::         &
  T_ADV_KEY = 'ADV'
#ifndef MESSY
INTEGER (KIND=iintegers) ::             &
  T_ADV_ID = -8
#else
INTEGER (KIND=iintegers), PARAMETER :: T_ADV_ID = I_ADVECT
#endif
INTEGER (KIND=iintegers), PARAMETER  :: &
  T_ADV_OFF     = 0_iintegers,      & ! don't advect this tracer
  T_ADV_ON      = 1_iintegers,      & ! advect this tracer
  T_ADV_MAX     = 1_iintegers         ! maximum allowable parameter value
INTEGER (KIND=iintegers) ::             &
  T_ADV_DEFAULT = T_ADV_OFF

! horizontal hyperdiffusion (integer)
CHARACTER (LEN=*), PARAMETER ::         &
  T_DIFF_KEY = 'DIFF'
#ifndef MESSY
INTEGER (KIND=iintegers) ::             &
  T_DIFF_ID = -9
#else
INTEGER (KIND=iintegers), PARAMETER :: T_DIFF_ID = I_HDIFF
#endif
INTEGER (KIND=iintegers), PARAMETER  :: &
  T_DIFF_OFF    = 0_iintegers,      & ! don't diffuse this tracer
  T_DIFF_ON     = 1_iintegers,      & ! diffuse this tracer
  T_DIFF_MAX    = 1_iintegers         ! maximum allowable parameter value
INTEGER (KIND=iintegers) ::             &
  T_DIFF_DEFAULT = T_DIFF_OFF

! turbulent mixing (integer)
CHARACTER (LEN=*), PARAMETER ::         &
  T_TURB_KEY = 'TURB'
#ifndef MESSY
INTEGER (KIND=iintegers) ::             &
  T_TURB_ID = -10
#else
INTEGER (KIND=iintegers), PARAMETER :: T_TURB_ID = I_VDIFF
#endif
INTEGER (KIND=iintegers), PARAMETER  :: &
  T_TURB_OFF    = 0_iintegers,      & ! don't do turbulent mixing with this tracer
  T_TURB_1D     = 1_iintegers,      & ! do 1-dimensional (vertical) turbulent mixing with this tracer
  T_TURB_3D     = 2_iintegers,      & ! do 3-dimensional turbulent mixing with this tracer
  T_TURB_MAX    = 2_iintegers         ! maximum allowable parameter value
INTEGER (KIND=iintegers) ::             &
  T_TURB_DEFAULT = T_TURB_OFF

! passive transport by convection (integer)
CHARACTER (LEN=*), PARAMETER ::         &
  T_CONV_KEY = 'CONV'
#ifndef MESSY
INTEGER (KIND=iintegers) ::             &
  T_CONV_ID = -11
#else
INTEGER (KIND=iintegers), PARAMETER :: T_CONV_ID = I_CONVECT
#endif
INTEGER (KIND=iintegers), PARAMETER  :: &
  T_CONV_OFF    = 0_iintegers,      & ! no convective passive transport for this tracer
  T_CONV_ON     = 1_iintegers,      & ! convective passive transport for this tracer
  T_CONV_MAX    = 1_iintegers         ! maximum allowable parameter value
INTEGER (KIND=iintegers) ::             &
  T_CONV_DEFAULT = T_CONV_OFF

! type of initial condition (integer)
CHARACTER (LEN=*), PARAMETER ::         &
  T_INI_KEY = 'INI'
#ifndef MESSY
INTEGER (KIND=iintegers) ::             &
  T_INI_ID = -12
#else
INTEGER (KIND=iintegers), PARAMETER :: T_INI_ID = I_INITIAL
#endif
INTEGER (KIND=iintegers), PARAMETER  :: &
  T_INI_ZERO    = 0_iintegers,      & ! initialize to zero
  T_INI_FILE    = 1_iintegers,      & ! initialize tracer from analysis file(laf)
  T_INI_USER    = 2_iintegers,      & ! initialization let to the user
  T_INI_MAX     = 2_iintegers         ! maximum allowable parameter value
INTEGER (KIND=iintegers) ::             &
  T_INI_DEFAULT = T_INI_ZERO

! type of lateral boundary condition (integer)
CHARACTER (LEN=*), PARAMETER ::         &
  T_LBC_KEY = 'LBC'
#ifndef MESSY
INTEGER (KIND=iintegers) ::             &
  T_LBC_ID = -13
#else
INTEGER (KIND=iintegers), PARAMETER :: T_LBC_ID = I_LBC
#endif
INTEGER (KIND=iintegers), PARAMETER  :: &
  T_LBC_ZERO     = 0_iintegers,     & ! initialize lateral BC to zero
  T_LBC_FILE     = 1_iintegers,     & ! initialize tracer from boundary condition file (lbf)
  T_LBC_CST      = 2_iintegers,     & ! constant lateral boundary conditions
  T_LBC_ZEROGRAD = 3_iintegers,     & ! zero-gradient lateral boundary conditions
  T_LBC_USER     = 4_iintegers,     & ! lateral boundary conditions have to be defined by the user
  T_LBC_MAX      = 4_iintegers        ! maximum allowable parameter value
INTEGER (KIND=iintegers) ::             &
  T_LBC_DEFAULT = T_LBC_ZERO

! type of bottom boundary conditions (integer)
CHARACTER (LEN=*), PARAMETER ::         &
  T_BBC_KEY = 'BBC'
INTEGER (KIND=iintegers) ::             &
  T_BBC_ID = -14
INTEGER (KIND=iintegers), PARAMETER  :: &
  T_BBC_ZEROFLUX     = 0_iintegers, & ! zero flux at the bottom 
  T_BBC_ZEROVAL      = 1_iintegers, & ! zero value at the bottom
  T_BBC_SURF_VAL     = 2_iintegers, & ! bottom values provided in a surface field
  T_BBC_MAX          = 2_iintegers    ! maximum allowable parameter value
INTEGER (KIND=iintegers) ::             &
  T_BBC_DEFAULT = T_BBC_ZEROFLUX

! type of relaxation (integer)
CHARACTER (LEN=*), PARAMETER ::         &
  T_RELAX_KEY = 'RELAX'
#ifndef MESSY
INTEGER (KIND=iintegers) ::             &
  T_RELAX_ID = -15
#else
INTEGER (KIND=iintegers), PARAMETER :: T_RELAX_ID = I_RELAX
#endif
INTEGER (KIND=iintegers), PARAMETER  :: &
  T_RELAX_OFF        = 0_iintegers, & ! no relaxation
  T_RELAX_FULL       = 1_iintegers, & ! full relaxation (at the 4 boundaries)
  T_RELAX_INFLOW     = 2_iintegers, & ! relaxation at inflow boundaries only
  T_RELAX_MAX        = 2_iintegers    ! maximum allowable parameter value
INTEGER (KIND=iintegers) ::             &
  T_RELAX_DEFAULT = T_RELAX_OFF

! type of Rayleigh damping (integer)
CHARACTER (LEN=*), PARAMETER ::         &
  T_DAMP_KEY = 'DAMP'
#ifndef MESSY
INTEGER (KIND=iintegers) ::             &
  T_DAMP_ID = -16
#else
INTEGER (KIND=iintegers), PARAMETER :: T_DAMP_ID = I_DAMP
#endif
INTEGER (KIND=iintegers), PARAMETER  :: &
  T_DAMP_OFF         = 0_iintegers, & ! no Rayleigh damping
  T_DAMP_ON          = 1_iintegers, & ! Rayleigh damping
  T_DAMP_MAX         = 1_iintegers    ! maximum allowable parameter value
INTEGER (KIND=iintegers) ::             &
  T_DAMP_DEFAULT = T_DAMP_OFF

! type of clipping to apply (integer)
CHARACTER (LEN=*), PARAMETER ::         &
  T_CLP_KEY = 'CLP'
#ifndef MESSY
INTEGER (KIND=iintegers) ::             &
  T_CLP_ID = -17
#else
INTEGER (KIND=iintegers), PARAMETER :: T_CLP_ID = -99
#endif
INTEGER (KIND=iintegers), PARAMETER  :: &
  T_CLP_OFF     = 0_iintegers,      & ! no clipping
  T_CLP_ON      = 1_iintegers,      & ! clipping to ensure positive definiteness 
  T_CLP_MAX     = 1_iintegers         ! maximum allowable parameter value
INTEGER (KIND=iintegers) ::             &
  T_CLP_DEFAULT = T_CLP_OFF

! 3. Defined derived data type for the tracer
! ------------------------------------------------

TYPE trcr_type
  ! (internal) metadata
  INTEGER (KIND=iintegers)    :: igribrep       ! repetition of the same combination grib param-grib table (=iloc1)

  ! index for the variables holding the tracer data
  INTEGER (KIND=iintegers)    :: idx_trcr       ! index for the tracer data
  INTEGER (KIND=iintegers)    :: idx_bd         ! index for the tracer boundary data
END TYPE trcr_type

TYPE (trcr_type), SAVE        :: trcr(n_trcr_max) ! list of all tracers

TYPE (t_metadata), SAVE       :: mtrcr          ! tracer meta-data storage


! 4. Variable for tracer data holding
! -----------------------------------

#ifndef MESSY
REAL (KIND=wp),     TARGET, ALLOCATABLE :: &
#else
REAL (KIND=wp),     POINTER             :: &
#endif
  trcr_data      (:,:,:,:,:),          & ! tracer data storage
  trcr_data_bd   (:,:,:,:,:),          & ! tracer boundary data storage
  trcr_data_tens (:,:,:,:)               ! tracer tendency data storage

!==============================================================================

END MODULE data_tracer
