!+ Variables for the computation of synthetic satellite images
!------------------------------------------------------------------------------

MODULE data_satellites

!------------------------------------------------------------------------------
!
! Description:
!  This data module contains all data necessary for the computation of 
!  synthetic satellite images.
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 3.7        2004/02/18 Ulrich Schaettler
!  Initial release
! V4_9         2009/07/16 Ulrich Schaettler
!  Put some fields to longtime memory for calling RTTVI only once
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_18        2011/05/26 Robin Faulwetter
!  Additional data necessary for RTTOV9
! V4_26        2012/12/06 Andreas Messer
!  Modifications for using also RTTOV10 and for satellite observation processing
! V4_27        2013/03/19 Ulrich Schaettler
!  Introduced nmsgchan as global variable (was local in src_output, netcdf_io before)
! V4_28        2013/07/12 Ulrich Schaettler
!  Extensions to provide GRIB2 shortnames and additional meta data
!  Renamed nlist_chan to nchan_list (for consistent naming)
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
! V5_2         2015-05-21 Annika Schomburg
!  Two new variables to read NWCSAF SEVIRI cloud products from grib-files
!   (only in affect if lobsrad=.TRUE.)
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
!
USE data_parameters, ONLY : &
    wp,            & ! KIND-type parameters for real variables
    iintegers        ! KIND-type parameter for standard integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Global (i.e. public) Declarations:

! 1. Global parameters from the RTTOV-module MOD_CPARAM
!------------------------------------------------------

!------------------------------------------------------------------------------
! These are global parameters which users may want to edit to optimise
! for their application
!------------------------------------------------------------------------------

  INTEGER (KIND=iintegers), PARAMETER ::    &
     jppf   =  1,    & ! maximal number of profiles per RTTOV call
     jpch   = 10,    & ! maximal number of channels
     jpchus = 10,    & ! maximal number of channels computed/call
     jpnsat =  9,    & ! maximal number of sensors to be used
     jplev  = 43,    & ! number of pressure levels
     jpnav  =  4,    & ! number of profile variables
     jpnsav =  5,    & ! number of surface air variables
     jpnssv =  6       ! number of skin variables

  REAL    (KIND=wp)       , PARAMETER ::    &
     rcnv = 6.03504E5_wp,    & ! kg/kg--> ppmv ozone
     rcnw = 1.60771704E6_wp    ! kg/kg--> ppmv water vapour


! 2. Variables for the organisation of the satellite computations
!----------------------------------------------------------------

TYPE sat_check_type
  CHARACTER (LEN= 8)      :: ysatellite         ! platform, e.g. METEOSAT, MSG,
  INTEGER(KIND=iintegers) :: nsat_id            ! Satellite identification 
                                                !   within rttov
  INTEGER(KIND=iintegers) :: nsat_id_min        ! for range of satellite ids
  INTEGER(KIND=iintegers) :: nsat_id_max        ! for range of satellite ids
  CHARACTER (LEN=12)      :: ysensor            ! Name of sensor used
  INTEGER(KIND=iintegers) :: nrttov_id          ! sensor identification within 
                                                !   rttov
  INTEGER(KIND=iintegers) :: nchan_min          ! for range of channels
  INTEGER(KIND=iintegers) :: nchan_max          ! for range of channels  
  INTEGER(KIND=iintegers) :: num_chan_max       ! Max. Number of channels for 
                                                !   that sensor
END TYPE sat_check_type

TYPE sat_org_type
  CHARACTER (LEN= 8)      :: ysatellite         ! platform, e.g. METEOSAT, MSG,
  INTEGER(KIND=iintegers) :: nsatell_table_id   ! entry in rttov satellite table
  INTEGER(KIND=iintegers) :: nsat_id            ! Satellite identification
  INTEGER(KIND=iintegers) :: wmo_satid          ! WMO satellite identification
  CHARACTER (LEN=12)      :: ysensor            ! Name of sensor used
  INTEGER(KIND=iintegers) :: nsensor_table_id   ! entry in rttov sensor table
  REAL   (KIND=wp)        :: longitude          ! position of geost. satellite
  INTEGER(KIND=iintegers) :: num_chan           ! Number of channels used
  INTEGER(KIND=iintegers) :: nchan_list(jpch)   ! List of channels used
  CHARACTER (LEN=10)      :: ychan_name(jpch)   ! Names of channels used (IRx.y, WVx.y)
  REAL   (KIND=wp)        :: emissivity(jpch)   ! emissivities for all channels
  LOGICAL                 :: lclear_rad         ! for clear sky radiance 
  LOGICAL                 :: lcloud_rad         ! for cloudy sky radiance 
  LOGICAL                 :: lclear_tem         ! for clear sky temperature
  LOGICAL                 :: lcloud_tem         ! for cloudy sky temperature
  INTEGER(KIND=iintegers) :: ngrib_chan(4*jpch) ! list of channels for grib 
                                                !   output
  INTEGER(KIND=iintegers) :: ngrib_aees(4*jpch) ! list of additional element
                                                ! numbers for grib output
  INTEGER(KIND=iintegers) :: ndim3_field        ! 3. dimension for LM variables
END TYPE sat_org_type

TYPE (sat_org_type)  :: sat_compute(jpnsat)

! 3. Variables for the organisation of the Namelist Input
!--------------------------------------------------------

TYPE sat_input_type
  CHARACTER (LEN= 8)      :: ysatellite         ! platform, e.g. METEOSAT, MSG,
  INTEGER(KIND=iintegers) :: nsat_id            ! Satellite identification
  CHARACTER (LEN=12)      :: ysensor            ! Name of sensor used
  INTEGER(KIND=iintegers) :: num_chan           ! Number of channels used
  LOGICAL                 :: lclear_rad         ! for clear sky radiance 
  LOGICAL                 :: lcloud_rad         ! for cloudy sky radiance 
  LOGICAL                 :: lclear_tem         ! for clear sky temperature
  LOGICAL                 :: lcloud_tem         ! for cloudy sky temperature
END TYPE sat_input_type

LOGICAL                          ::           &
  lcon_clw        ! if .TRUE.: convective liquid water used in rttov

! 4. Additional control variables
!--------------------------------

INTEGER (KIND=iintegers)  ::                    &
   itype_rttov = 7, & ! Version of RTTOV model
   num_sensors,     & ! No. of sensors used
   nsat_next          ! for the organization of the computations

INTEGER(KIND=iintegers), ALLOCATABLE   :: nsat_steps(:)
         ! list of time steps for which satellite computations must be done


! the following are the RTTOV satellite and sensor tables 
! (Table 3 from the RTTOV documentation)

CHARACTER (LEN= 8)    :: yrttov_satell_table(  13)
CHARACTER (LEN=12)    :: yrttov_sensor_table(0:26)

LOGICAL                          ::           &
  lsynsat,          & ! to produce the synthetic satellite images
  lobsrad,          & ! to process satellite observations and compute radiances
  linterp             ! do interpolation of p, t, q to half-levels

! 5. Variables that are initialized by RTTVI
!-------------------------------------------

! These variables are the same for all instruments, but are initialized
! by the call to RTTVI, which is done only once for the whole program

INTEGER (KIND=iintegers) ::           &
  maxknchpf          ,  & ! maximum number of output radiances
  numchans   (jpnsat),  & ! Number of valid channels
  kiu1       (jpnsat)     ! for input-file unit number

REAL  (KIND=wp),     ALLOCATABLE ::               &
  o3_ref     (:)          ! default ozone values on the prescribed levels

! These variables can be different for every instrument and are initialized
! by the call to RTTVI, which is done only once for the whole program

REAL  (KIND=wp),     ALLOCATABLE ::               &
  ppres_d    (:,:),     & ! default pressure on the prescribed levels
  utmx       (:,:),     & ! maximum temperature
  utmn       (:,:),     & ! minimum temperature
  uqmx       (:,:),     & ! maximum humidity
  uqmn       (:,:),     & ! minimum humidity
  uomx       (:,:),     & ! maximum ozone
  uomn       (:,:)        ! minimum ozone

INTEGER (KIND=iintegers), ALLOCATABLE :: &
  ivch       (:,:)

! 6. Some constant variables
!---------------------------

REAL  (KIND=wp)     ::           &
  const_aheat, r_sat

! for output of MSG-variables
INTEGER (KIND=iintegers), PARAMETER  ::          &
  nmsgchan = 8


! 7. Data structures and Variables for RTTOV9 and higher
!-------------------------------------------------------

! Required for initialization of RTTOV9/10
! (function rttov_init of mo_rttov_ifc)

INTEGER(KIND=iintegers), ALLOCATABLE   :: instruments(:,:)
  ! for every sensor (2nd dimension):
  !    nsatell_table_id
  !    nsat_id
  !    nsensor_table_id

INTEGER(KIND=iintegers), ALLOCATABLE   :: channels(:,:)
  ! for every sensor (2nd dimension):
  !    wmo_satid
  !    nsatell_table_id

INTEGER(KIND=iintegers), ALLOCATABLE   :: n_chans(:)
  ! for every sensor the number of channels

LOGICAL,                 ALLOCATABLE   :: addclouds(:)
  ! for every sensor if it shall use ir cloud scattering

INTEGER(KIND=iintegers)                :: mchans
  ! maximum number of channels for one sensor (MAX (numchans))

REAL(Kind=wp), PARAMETER :: zenmax9  = 86.5_wp        ! deg
REAL(Kind=wp), PARAMETER :: zenmax10 = 75.0_wp        ! deg
  ! maximum satellite zenith angles

REAL(KIND=wp),     PARAMETER ::              &
  p_top = 9.9_wp,                 & ! Highest pressure level (in Pa),
                                    !     that is used for input to RTTOV
  t_top = 231.6_wp,               & ! Standard temperature  [K]
  w_top = 0.349555E-05_wp,        & ! Standard mixing ratio [Kg/Kg]
  q_top = w_top / (1._wp + w_top)   !
                                    !  all at highest pressure level

   ! Bits for extrapolation of RTTOV input profiles above model top:
INTEGER(KIND=iintegers), PARAMETER :: extrp_logp    = 1_iintegers !extrapolate log(p) linearly instead of p
INTEGER(KIND=iintegers), PARAMETER :: extrp_const   = 2_iintegers !extrapolate t, q with constant values
INTEGER(KIND=iintegers), PARAMETER :: extrp_lin     = 4_iintegers !extrapolate t, q linearly
INTEGER(KIND=iintegers), PARAMETER :: extrp_clim    = 8_iintegers !extrapolate t, q to climatological value
INTEGER(KIND=iintegers)            :: extrp_type    = extrp_const

INTEGER(KIND=iintegers) :: iwc2effdiam   = 4_iintegers
     ! 1: Scheme by Ou and Liou, 1995, Atmos. Res., 35, 127-138.
     ! 2: Scheme by Wyser et al. (see McFarquhar et al. (2003))
     ! 3: Scheme by Boudala et al., 2002, Int. J. Climatol., 22, 1267-1284.
     ! 4: Scheme by McFarquhar et al. (2003)
     !RF  Don't use 1-3, these values might cause floating point exceptions!!
INTEGER(KIND=iintegers) :: iceshape = 1_iintegers
     ! 1: hexagonal
     ! 2: ice aggregates

! 8. Namelist variables for reading SEVIRI NWCSAF cloud products
!---------------------------------------------------------------------------
LOGICAL             :: lread_ct              ! namelist switch: read cloud type?
CHARACTER (LEN=100) :: yclouddir             ! directory of cloud data

!==============================================================================

END MODULE data_satellites
