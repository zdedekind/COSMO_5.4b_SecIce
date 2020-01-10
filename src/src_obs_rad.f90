!+ Module for generating feedback files for satellite radiances.
!------------------------------------------------------------------------------

MODULE src_obs_rad

!------------------------------------------------------------------------------
!
! Description:
!   This module contains data and routines for calculating first guess for
!   satellite observations.
!
! Method:
!
! Current Code Owner: DWD, Robin Faulwetter
!  phone:  +49  69  8062 2746
!  fax:    +49  69  8062 3721
!  email:  robin.faulwetter@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! V4_26        2012/12/06 Andreas Messer
!  Initial Release
! V4_27        2013/03/19 Ulrich Schaettler
!  Introduced conditional compilation for Nudging to separate nudging parts
!  from SYNSAT parts.
! V4_28        2013/07/12 Ulrich Schaettler
!  Changed pointer assignment to standard assignment, because variable is
!   no pointer any more (SX accepted that gracefully, gfortran reports error)
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
! V4_30        2013/11/08 Ulrich Schaettler
!  Changed intent attribute of variable nlev to INOUT in SR prepare_rttov_input
!  Changed intent attribute of variable comm to INOUT in SR p_bcast_rad_set
! V5_1         2014-11-28 Ulrich Blahak, Oliver Fuhrer
!  Changed the format of some YUSPECIF entries for the CLM namelist tool. (UB)
!  Implemented F2003 IOMSG-mechanism for better namelist error messages. (UB)
!  Initialize nsensors with 0; otherwise it is not initialized, if lobsrad=False (OF)
!  Replaced ireals by wp (working precision) (OF)
! V5_2         2015-05-21 Annika Schomburg, Robin Faulwetter
!  Added routines cloud_read_prep, grb_open, grb_read for reading NWCSAF SEVIRI 
!    cloud type and cloud top height information mapped on COSMO grid (from GRIB-file)
!  Some changes to determine HRIT file names and edges internally, not by 
!    namelist variables, some changes to namelist TOVS_OBS.
!  Fixed MPI communication for asynchronous IO.
!  Added superobbing option for HRIT.
!  Improved calculation of input profiles for RTTOV.
! V5_3         2015-10-09 Robin Faulwetter
!  Corrected interpolation of input for RTTOV (linterp=T)
! V5_3a        2015-11-24 Ulrich Schaettler
!  Replaced a GRIBDWD by the correct GRIBAPI
! V5_4         2016-03-10 Robin Faulwetter
!  Bug fix in SR input_obs_satpp: had to move a part outside an IF-statement
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

#ifdef GRIBAPI
USE grib_api
#endif

USE iso_fortran_env, ONLY: stderr => error_unit, stdout => output_unit

USE data_parameters, ONLY: iintegers, wp, sp, intgribf, irealgrib, intgribc

USE data_constants, ONLY: B1, B2W, B3, B4W, rdv, o_m_rdv, b2i, b4i, &
                            t0_melt,      & ! melting temperature of ice
                            rdocp,        & ! r_d / cp_d
                            lh_v,         & ! latent heat of vapourization
                            cpdr,         & ! 1 / cp_d
                            uc1,          & ! variable for computing the rate of cloud cover in 
                            uc2,          & ! the unsaturated case
                            ucl


USE data_parallel, ONLY: &
  nproc,            & ! total number of processors
  num_compute,      & ! number of compute pe's
  nboundlines,      & ! no of overlapping boundary lines of the subdomains
  isubpos,          & ! positions of subdomains in total domain
  my_world_id,      & ! rank of this subdomain in the global    communicator
  my_cart_id,       & ! rank of this subdomain in the cartesian communicator
  icomm_world,      & ! global communicator
  icomm_cart,       & ! cartesian communicator
  imp_character,    & ! character type for MPI
  imp_reals,        & ! character type for MPI
  imp_integers,     & ! integer type for MPI
  imp_logical         ! logical type for MPI


USE data_modelconfig, ONLY: &
  ie,             & ! no of gridpoints in zonal direction (local)
  je,             & ! no of gridpoints in meridional direction (local)
  ke,             & ! no of gridpoints in vertical direction
  ie_tot,         & ! no of gridpoints in zonal direction
  je_tot,         & ! no of gridpoints in meridional direction
  ke_tot,         & ! no of gridpoints in vertical direction
  pollon,         & ! logitude of the rotated north pole
  pollat,         & ! latitude of the rotated north pole
  polgam,         & ! angle between the north poles of the system
  dlon,           & ! grid point distance in zonal direction
  dlat,           & ! grid point distance in meridonal direction
  startlat_tot,   & ! rot latitude lower left gridpoint 
  startlon_tot,   & ! rot longitude lower left gridpoint
  degrad,         & ! factor for transforming degree to rad
  raddeg,         & ! factor for transforming rad to degree
  dt,             & ! Length of Timestep
  idt_qv, idt_qc, idt_qi, idt_qs, idt_qg

#ifdef NUDGING
USE data_nudge_all, ONLY: &
  nolbc,            & ! no of grid rows at lateral boundaries where obs are neg.
  doromx            ! station heights
#endif

USE data_runcontrol, ONLY: &
  nstart,           & ! First timestep of forecast
  nstop,            & ! Last timestep of forecast
  ntstep,           & ! Current timestep of forecast
  nnow,             & ! field index for last timestep
  itype_calendar,   & ! for specifying the calendar used
  lseaice,          & ! 
  hstart,           & ! 
  hstop,            & ! 
  leps,             & ! 
  iepsmem,          & ! 
  nvers,            & ! 
  lprog_qi            ! prognostic cloud ice?                

USE data_io, ONLY: &
  ydate_ini,           & ! Start date of forecast
  yncglob_source,      & ! Start date of forecast
  yncglob_institution, & ! Start date of forecast
  root,                & !
  npds,                & ! Dimension for product definition section (pds)
  ngds,                & ! Dimension for grid description section (gds)
  nbms,                & ! Dimension for bit map section (bms)
  nbds,                & ! Dimension for binary data section
  ndims                  ! Dimension for idims

#ifdef NUDGING
USE data_obs_record, ONLY: &
  imdi,             & ! missing data indicator
  fdoro               ! scaling factor to vertical distances btw model
#endif

USE data_fields, ONLY: &
  p0,               & ! presure
  pp,               & ! presure
  ps,               & ! surface presure
  t,                & ! temperature
  t_g,              & ! surface temperature
  t_2m,             & ! 2m temperature
! qv,               & ! humidity
! qc,               & ! humidity
! qi,               & ! humidity
! qs,               & ! humidity
! qg,               & ! humidity
  qv_2m,            & ! 2m humidity
  u_10m,            & ! 10m wind u-component
  v_10m,            & ! 10m wind v-component
  rlat,             & ! rotated lon lat
  hsurf,            & ! surface height
  hhl,              & ! geommetricla height of half model levels
  soiltyp,          & ! soil type
  lseamask,         & ! land/sea mask
  llandmask,        & ! landpoint mask
  fr_land,          & ! geommetricla height of half model levels
  sun_el,           & ! sun
  h_ice,            & ! 
  clc_sgs,          & ! 
  clw_con,          & ! 
  clc_con,          & ! 
  rho                 ! 

USE data_satellites, ONLY: &
  extrp_logp,       &
  extrp_const,      &
  extrp_lin,        &
  extrp_clim,       &
  extrp_type,       &
  p_top,            & ! Top level pressure
  t_top,            & ! Top level temperature
  q_top,            & ! Top level humidity
  rcnw,             & ! kg/kg--> ppmv water vapour
  iceshape,         & !
  iwc2effdiam,      & !
  lcon_clw,         & !
  num_sensors,      & ! Number of sensors in use for synsat
  zenmax9,          & ! Maximum satellite zenith angle for RTTOV9
  zenmax10,         & ! Maximum satellite zenith angle for RTTOV10
  lread_ct,         & ! read NWCSAF cloud type classification?
  yclouddir,        & ! directory of cloud products
  linterp

#ifdef NUDGING
USE data_obs_lib_cosmo, ONLY: &
  rmdi,  &
  rmdich
#endif

USE src_tracer,       ONLY :  trcr_get, trcr_errorstr
USE data_tracer,      ONLY :  T_ERR_NOTFOUND

USE environment, ONLY: &
  model_abort,         & ! aborts the program
  get_free_unit,       & !
  comm_barrier

#ifdef NUDGING
USE src_obs_cdfin_util, ONLY: &
  obs_assign_gridpt   ,& ! Assign gridpoint to obs according lon/lat
  obs_assign_sort_node,& ! Assign node to gridpoint
  get_global_surface     ! get global surface
#endif

USE parallel_utilities, ONLY: &  
  distribute_field, & ! Distribute 2D field to all PEs  
  distribute_values,& ! Distribute values to all PEs
  gather_field,     & ! 
  gather_values,    & ! 
  global_values,    & ! 
  gatherv_values,   &
  scatterv_values

USE utilities, ONLY: &
  diff_minutes,     & ! Calculate time difference in minutes
  rla2rlarot,       & ! Convert lambda to rotated latlon
  phi2phirot,       & ! Convert phi to rotated latlon
  phirot2phi,       &
  rlarot2rla,       &
  get_utc_date        ! Determine date from timestep

#ifdef RTTOV10
USE mo_rad, ONLY: &
  t_rad_set,             & ! Type of an radiance set 
  t_radv,                & ! Type of a set of radiances
  t_obserr_par,          & ! Type for radiance obserrors
  rad_set,               & ! Array of radiance sets 
  read_tovs_obs_chan_nml,& ! Subroutine to read "TOVS_OBS_CHAN" Namelist
  read_satpp_feedbk     ,& ! Subroutine to read a satpp file
  link_rad              ,& ! Combine Information from Namelists with Satpp
  sat_id2name,           & ! Get the satellite name for a given bufr id
  destruct,              & ! Release memory  
  assignment (=),        & !
  USE_PASSIVE,           & ! Passive bit
  construct,             & ! constructor routines
  read_hrit,             & ! read hrit input files
  write_rtovp,           & ! write monRTOVP files
  WRM_OPEN,              & ! mode for write_rtovp
  WRM_WRITE,             & ! mode for write_rtovp
  WRM_CLOSE,             & ! mode for write_rtovp
  MP_H,                  & ! parameter to indicate that jacobians shall be written
  print_radv,            & ! print t_radv
  calc_obserr,           & ! calculate obserror based on namelist input
  superob_radv             ! superobbing routine

USE mo_rttov_ifc, ONLY:  &
  rttov_ifc_version,     & !
  rttov_fill_input,      & !
  rttov_set_input,       & !
  rttov_direct_ifc,      & !
  rttov_k_ifc,           & !
  rttov_ifc_errmsg,      & !
  rttov_print_profiles,  & !
  OUT_ASB,               & !
  OUT_CSB,               & !
  OUT_ASR,               & !
  OUT_CSR                  !
#endif

#ifdef NUDGING
USE mo_fdbk, ONLY: &
  t_fdbk,                & ! Type of feedback file
  setup_fdbk,            & ! setup a fdbk files variables
  create_fdbk,           & ! create the fdbk file
  close_fdbk,            & ! create the fdbk file
  open_fdbk_read,        & ! create the fdbk file
  open_fdbk_write,       & ! create the fdbk file
  add_verification         ! create the fdbk file

USE mo_fdbk_tables, ONLY: &
  vt_firstguess,          &
  vt_analysis,            &
  rc_ass,                 &
  OT_RAD,                 &
  VN_RAWBT,               &
  ST_SEA,                 &
  ST_ICE,                 &
  ST_LAND,                &
  ST_HIGHLAND,            &
  ST_MISMATCH,            &
  MS_LAND,                &
  MS_SEA,                 &
  MS_ICE,                 &
  MS_NO_ICE,              &
  MS_SNOW,                &
  MS_NO_SNOW,             &
  ST_ACCEPTED,            &
  ST_PASSIVE,             &
  ST_PAS_REJ,             &
  ST_DISMISS,             &
  FL_TIME,                &
  FL_AREA,                &
  FL_SURF,                &
  FL_OBSTYPE,             &
  FL_OPERATOR,            &
  FL_NONE,                &
  FL_DATASET,             &
  ST_ACTIVE,              &
  OF_RAD_CLEAR_SKY,       &
  OF_RAD_CLOUDY,          &   
  OF_BT_CLEAR_SKY 

USE mo_fdbk_cosmo, ONLY: &
  t_account,             &
  t_acc_header,          &
  t_acc_body,            &
  write_report

USE mo_netcdf_param  
#endif

IMPLICIT NONE

PRIVATE

#ifdef NUDGING
PUBLIC :: input_obs_radctl
PUBLIC :: input_obs_satpp
PUBLIC :: calc_obs_satpp
PUBLIC :: cloud_read_prep
#endif
PUBLIC :: t_sensor
PUBLIC :: nsensors
PUBLIC :: sensors
PUBLIC :: prepare_rttov_input
PUBLIC :: get_n_add

!==============================================================================
! MPI Include
!------------------------------------------------------------------------------
#include "mpif.h"

!==============================================================================
! Module Parameters
!------------------------------------------------------------------------------

  CHARACTER(LEN=*), PARAMETER :: &
    yurad = 'YURAD'                ! iradiance obs diagnostic file
  INTEGER, PARAMETER ::                &
    max_satpp_files         = 16    ! maximum number of satpp files
  INTEGER, PARAMETER ::                &
       max_hrit_files         = 300    ! maximum number of hrit files
  REAL (KIND=wp), PARAMETER :: &
    tmax = 399.0_wp,           & ! [K]
    tmin = 91.0_wp,            & ! [K]
    qmax = 0.372_wp,           & ! [kg/kg]
    qmin = TINY(0._wp)            ! [kg/kg]

#ifdef NUDGING
  ! The following array maps a check to its assigned state  
  INTEGER (KIND=iintegers), PARAMETER :: &
    map_flg2state(0:21) = (/ST_PASSIVE,  & ! FL_OBSTYPE
                            ST_DISMISS,  & ! FL_BLACKLIST
                            ST_PAS_REJ,  & ! FL_SUSP_LOCT
                            ST_DISMISS,  & ! FL_TIME
                            ST_DISMISS,  & ! FL_AREA
                            ST_PASSIVE,  & ! FL_HEIGHT
                            ST_PASSIVE,  & ! FL_SURF
                            ST_PASSIVE,  & ! FL_CLOUD
                            ST_PAS_REJ,  & ! FL_PRACTICE
                            ST_PAS_REJ,  & ! FL_DATASET    (3dvar TOVS convetion)
                            ST_PASSIVE,  & ! FL_REDUNDANT
                            ST_PAS_REJ,  & ! FL_FLIGHTTRACK
                            ST_DISMISS,  & ! FL_MERGE
                            ST_DISMISS,  & ! FL_THIN
                            ST_PAS_REJ,  & ! FL_RULE
                            ST_PAS_REJ,  & ! FL_OBS_ERR
                            ST_PASSIVE,  & ! FL_GROSS
                            ST_PASSIVE,  & ! FL_NO_BIASCOR
                            ST_PASSIVE,  & ! FL_FG
                            ST_PAS_REJ,  & ! FL_NO_OBS
                            ST_PAS_REJ,  & ! FL_OPERATOR
                            ST_DISMISS/)   ! FL_FG_LBC

  ! Category and subcategory
  INTEGER (KIND=iintegers), PARAMETER :: &
    fdbk_rad_cat    = 21, & ! value to use for category field in fdbk file
    fdbk_rad_subcat = -1    ! value to use fur sub_cat.. field in fdbk file


  INTEGER (KIND=iintegers), PARAMETER :: &
    fdbk_height_method = 2  ! Which method to use for calculating the height
                            ! of an obs:
                            ! 1 = max of temp jacobian
                            ! 2 = mean of jacobian weighted by temp/hum
                            ! 3 = like (2) but squared sensitivities
  REAL (KIND=wp), PARAMETER :: &
    fdbk_e_bg_t  = 0.5_wp, & ! Typical background temp error
    fdbk_e_bg_rh = 0.1_wp    ! Typical background hum error
#endif


  REAL   (KIND=irealgrib), ALLOCATABLE :: &
    ds  (:)          ! array for unpacked data

  INTEGER (KIND=intgribf)  ::  &
    idims (ndims) ,& ! array for all dimension reading a GRIB file
    ipds  (npds)  ,& ! product definition section for GRIB file
    igds  (ngds)  ,& ! grid description section for GRIB file
    ibms  (nbms)  ,& ! bit map section for GRIB file
    ibds  (nbds)  ,& ! binary data section for GRIB file
    iednr         ,& ! grib edition number
    lds           ,&
    lbm           ,&
    lfd

!==============================================================================
! Module Types
!------------------------------------------------------------------------------

  TYPE t_sensor
    INTEGER (KIND=iintegers) :: &
      id(3)                 ! RTTOV instrument triplet (platform,sat,instr)
    INTEGER (KIND=iintegers), POINTER :: &
      channels(:) => NULL() ! RTTOV channels for this instrument
    LOGICAL :: &
      addcloud              ! Use cloud information for this sensor
  END TYPE t_sensor

  TYPE t_time
    INTEGER (KIND=iintegers) :: &
      year,  &
      month, &
      day,   &
      hour,  &
      minute
  END TYPE t_time

!==============================================================================
! Module Variables
!------------------------------------------------------------------------------

  INTEGER  (KIND=iintegers), SAVE ::  &
    nsensors=0_iintegers              ! Number of sensors
  TYPE (t_sensor), SAVE, ALLOCATABLE, TARGET :: &
    sensors(:)                          ! Array of sensors
  CHARACTER(LEN=128), SAVE ::         &
    input_path = ''                     ! path to input files (from satpp)
  CHARACTER(LEN=128), SAVE ::         &
    satpp_files(max_satpp_files) = ''   ! satpp files to read in
  CHARACTER(LEN=128), SAVE ::         &
    input_path_hrit = ''                ! path to hrit input files
  CHARACTER(LEN=128), SAVE ::         &
    hrit_files(max_hrit_files) = ''     ! hrit files to read in
  REAL (KIND=wp), SAVE ::             &
    minlat =  -90._wp                   ! minimum latitude for data from hrit files
  REAL (KIND=wp), SAVE ::             &
    maxlat =   90._wp                   ! maximum latitude for data from hrit files
  REAL (KIND=wp), SAVE ::             &
    minlon = -360._wp                   ! minimum longitude for data from hrit files
  REAL (KIND=wp), SAVE ::             &
    maxlon =  360._wp                   ! maximum longitude for data from hrit files

#ifdef RTTOV10    
  INTEGER  (KIND=iintegers), SAVE ::  &
    nurad                               ! radiance obs information file

  INTEGER  (KIND=iintegers), SAVE, ALLOCATABLE :: &
    nobspe_gath(:,:),                 & ! number of obs gathered (pe,radv) (local)
    nobspe_tot(:,:),                  & ! number of total obs   (pe,radv) (global)
    fdbk_offset(:)                      ! offset within feedback fie (pe 0 only)
  INTEGER   (KIND=iintegers), SAVE :: &
    nprof_rttov= -1                     ! maximum number of profiles in RTTOV-call
  INTEGER   (KIND=iintegers), SAVE :: &
    monitor_prof = 0                    ! control of monRTOVP files
  INTEGER   (KIND=iintegers), SAVE :: &
    usesurftype = 1                     ! Which surfacetype to use: 0,1 = cosmo,satpp
  TYPE (t_radv), SAVE, ALLOCATABLE, TARGET :: &
    obs_rad(:)                          ! Array of radiances for each satellite dataset
  INTEGER (KIND=iintegers), SAVE   :: &
    nsatpp,                           & ! number of satpp files
    nso,                              & ! superobbing box size
    nbox                                ! number of superobbing boxes

#ifdef NUDGING
  TYPE (t_fdbk), SAVE, ALLOCATABLE :: &
    fdbk(:)                             ! array of fdbk files
  TYPE (t_time), SAVE ::              &
    reftime                             ! reference time used in varous places
#endif

  REAL (KIND=wp), SAVE ::             &
    zenmax = zenmax10                ! maximum satellite zenith angle
#elif defined(RTTOV9)
  REAL (KIND=wp), SAVE ::             &
    zenmax = zenmax9                 ! maximum satellite zenith angle
#endif

  INTEGER  (KIND=iintegers),SAVE ::   &
#ifdef RTTOV10
       rad_output   = ibset(0,OUT_ASB),&! radiance output
#endif
       thin_line    = 1,              & ! Line thinning for HRIT files
       thin_col     = 1,              & ! Line thinning for HRIT files
       hrit_superob = 1                 ! Superobbing box size
  LOGICAL, SAVE ::                    &
       calc_k    = .false.,           & ! Calculate k-matrices for radiances
       l_H       = .false.              ! write Jacobians in feedback files

  CHARACTER (LEN=4)  ::               &
       hrit_name    = ''                ! name of geostationary satellite
  CHARACTER (LEN=6), SAVE  ::         &
       hrit_chan(8) = ''                ! channel names
  REAL (KIND = wp), SAVE ::           &
       hrit_hh(3) = (/0.0,1.0,1.0/)     ! triple at which hours HRIT-files should be read

!==============================================================================
! Module Interfaces
!------------------------------------------------------------------------------

#ifdef NUDGING
#ifdef RTTOV10
  INTERFACE p_bcast
    MODULE PROCEDURE p_bcast_rad_set
    MODULE PROCEDURE p_bcast_obserr_par
  END INTERFACE
#endif  
#endif

!==============================================================================

CONTAINS

#ifdef NUDGING
#ifdef RTTOV10
!==============================================================================
!+ Procedure in "src_obs_sat" for distributing the rad container
!------------------------------------------------------------------------------

SUBROUTINE gather_radv (c_loc, start_loc, end_loc, &
                        c_tot, start_tot, end_tot, &
                        ireceiver, ierrstat,       &
                        numreceived)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine collects radiances in one place
!
!------------------------------------------------------------------------------

! Parameters
  INTEGER   (KIND=iintegers),   INTENT (IN)   :: &
    start_loc,    & ! First Record to distribute
    end_loc,      & ! Last Record        
    start_tot,    & ! First Record to receive
    ireceiver       ! Node which shall receive the data
  TYPE(t_radv), INTENT(INOUT), TARGET :: &
    c_loc     ! Input radiances
  TYPE(t_radv), INTENT(INOUT) :: &
    c_tot     ! Collected radiances
  INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
    ierrstat        ! error status variable
  INTEGER   (KIND=iintegers),   INTENT (INOUT), OPTIONAL :: &
    numreceived(num_compute)      ! number of records received  
  INTEGER   (KIND=iintegers),   INTENT (OUT)   :: &
    end_tot         ! Last Record received 

! Local variables
  INTEGER   (KIND=iintegers) :: &
    inumelems(num_compute),     &  ! Number of obs per node
    icnt,                       &  ! Loop counter
    idum
  CHARACTER(LEN=256) :: &
    ymsg            ! error message

!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE gather_radv
!-------------------------------------------------------------------------------

  inumelems = 0
  inumelems(my_cart_id+1) = end_loc - start_loc + 1
  
  idum = inumelems(my_cart_id+1)
  CALL gather_values(idum, inumelems(:), 1, &
                     num_compute, imp_integers, ireceiver, icomm_cart, &
                     ymsg, ierrstat)
  IF (ierrstat /= 0) CALL model_abort(my_cart_id,99903, 'gather_values', &
                                   'gather_radv', ierrstat)

  end_tot = start_tot + SUM(inumelems) - 1

  CALL loc_gatherv_integers   (c_loc%ntstep,     c_tot%ntstep    )
  CALL loc_gatherv_integers   (c_loc%i_box ,     c_tot%i_box     )
  CALL loc_gatherv_integers   (c_loc%j_box ,     c_tot%j_box     )
  CALL loc_gatherv_integers   (c_loc%i_reprt,    c_tot%i_reprt   )
  CALL loc_gatherv_reals2d    (c_loc%plevel,     c_tot%plevel    )
  CALL loc_gatherv_integers   (c_loc%obsnum,     c_tot%obsnum    )
  CALL loc_gatherv_integers   (c_loc%date,       c_tot%date      )
  CALL loc_gatherv_integers   (c_loc%time,       c_tot%time      )
  CALL loc_gatherv_integers   (c_loc%date_d,     c_tot%date_d    )
  CALL loc_gatherv_integers   (c_loc%time_d,     c_tot%time_d    )
  CALL loc_gatherv_reals      (c_loc%dlat,       c_tot%dlat      )
  CALL loc_gatherv_reals      (c_loc%dlon,       c_tot%dlon      )
  CALL loc_gatherv_integers   (c_loc%fov,        c_tot%fov       )
  CALL loc_gatherv_integers   (c_loc%scanl,      c_tot%scanl     )
  CALL loc_gatherv_reals      (c_loc%stzen,      c_tot%stzen     )
  CALL loc_gatherv_reals      (c_loc%stazi,      c_tot%stazi     )
  CALL loc_gatherv_reals      (c_loc%sunzen,     c_tot%sunzen    )
  CALL loc_gatherv_reals      (c_loc%sunazi,     c_tot%sunazi    )
  CALL loc_gatherv_reals2d    (c_loc%landfr,     c_tot%landfr    )
  CALL loc_gatherv_integers2d (c_loc%stype,      c_tot%stype     )
  CALL loc_gatherv_reals2d    (c_loc%shgt,       c_tot%shgt      )
  CALL loc_gatherv_integers   (c_loc%center,     c_tot%center    )
  CALL loc_gatherv_integers   (c_loc%subcenter,  c_tot%subcenter )
  CALL loc_gatherv_integers   (c_loc%mdlsfc,     c_tot%mdlsfc    )
  CALL loc_gatherv_integers2d (c_loc%state,      c_tot%state     )
  CALL loc_gatherv_integers2d (c_loc%flags,      c_tot%flags     )
  CALL loc_gatherv_logical2d  (c_loc%valid,      c_tot%valid     )
  CALL loc_gatherv_reals2d    (c_loc%bt_obs,     c_tot%bt_obs    )
  CALL loc_gatherv_reals2d    (c_loc%bt_fg,      c_tot%bt_fg     )
  CALL loc_gatherv_reals2d    (c_loc%bt_bcor,    c_tot%bt_bcor   )
  CALL loc_gatherv_reals2d    (c_loc%bcor_,      c_tot%bcor_     )
  CALL loc_gatherv_reals2d    (c_loc%bt_fg_cs,   c_tot%bt_fg_cs  )
  CALL loc_gatherv_reals2d    (c_loc%rad_fg,     c_tot%rad_fg    )
  CALL loc_gatherv_reals2d    (c_loc%rad_fg_cs,  c_tot%rad_fg_cs )
  CALL loc_gatherv_reals2d    (c_loc%sinfl,      c_tot%sinfl     )
  CALL loc_gatherv_reals2d    (c_loc%emiss,      c_tot%emiss     )
  CALL loc_gatherv_reals2d    (c_loc%emis_pc,    c_tot%emis_pc   )
  CALL loc_gatherv_reals2d    (c_loc%p,          c_tot%p         )
  CALL loc_gatherv_reals2d    (c_loc%hh,         c_tot%hh        )
  CALL loc_gatherv_reals2d    (c_loc%h_fg,       c_tot%h_fg      )
  CALL loc_gatherv_reals2d    (c_loc%t_fg,       c_tot%t_fg      )
  CALL loc_gatherv_reals2d    (c_loc%q_fg,       c_tot%q_fg      )
  CALL loc_gatherv_reals      (c_loc%t2m,        c_tot%t2m       )
  CALL loc_gatherv_reals      (c_loc%q2m,        c_tot%q2m       )
  CALL loc_gatherv_reals      (c_loc%ps_fg,      c_tot%ps_fg     )
  CALL loc_gatherv_reals      (c_loc%ts_fg,      c_tot%ts_fg     )
  CALL loc_gatherv_reals      (c_loc%u10_fg,     c_tot%u10_fg    )
  CALL loc_gatherv_reals      (c_loc%v10_fg,     c_tot%v10_fg    )
  CALL loc_gatherv_reals      (c_loc%v10_abs_fg, c_tot%v10_abs_fg)
  CALL loc_gatherv_reals      (c_loc%cld_top,    c_tot%cld_top   )
  CALL loc_gatherv_reals3d    (c_loc%cld_fg,     c_tot%cld_fg    )
  CALL loc_gatherv_reals3d    (c_loc%cld_frc,    c_tot%cld_frc   )
  CALL loc_gatherv_reals3d    (c_loc%H_t,        c_tot%H_t       )
  CALL loc_gatherv_reals3d    (c_loc%H_q,        c_tot%H_q       )
  CALL loc_gatherv_reals2d    (c_loc%H_ts,       c_tot%H_ts      )
  CALL loc_gatherv_reals2d    (c_loc%H_ps,       c_tot%H_ps      )

CONTAINS

  SUBROUTINE loc_gatherv_integers(loc,tot)
    INTEGER (KIND=iintegers), INTENT(IN)   , POINTER :: loc(:)
    INTEGER (KIND=iintegers), INTENT(INOUT), POINTER :: tot(:)
    INTEGER (KIND=iintegers), POINTER :: locp(:)
    INTEGER (KIND=iintegers), TARGET  :: locd(1)
    INTEGER (KIND=iintegers) :: s_loc, numtot,i
    IF (ASSOCIATED(loc)) THEN
      if (size(loc) > 0) then
        locp => loc(start_loc:end_loc)
      else
        locp => locd
      end if
      call gatherv_values(locp,inumelems(my_cart_id+1),inumelems(my_cart_id+1),num_compute,&
                          tot, numtot, size(tot), ireceiver, imp_integers, icomm_cart,&
                          ymsg, ierrstat, init=.false.)
      IF (ierrstat /= 0) CALL model_abort(my_cart_id,99904, 'MPI_GATHERV', &
           'gather_radv', ierrstat)
    END IF
  END SUBROUTINE loc_gatherv_integers


  SUBROUTINE loc_gatherv_integers2d (loc,tot)
    INTEGER (KIND=iintegers), INTENT(IN)   , POINTER :: loc(:,:)
    INTEGER (KIND=iintegers), INTENT(INOUT), POINTER :: tot(:,:)
    INTEGER (KIND=iintegers), POINTER :: locp(:,:)
    INTEGER (KIND=iintegers), TARGET  :: locd(1,1)
    INTEGER (KIND=iintegers) :: s_loc, numtot, i
    IF (ASSOCIATED(loc)) THEN
      if (size(loc) > 0) then
        locp => loc(:,start_loc:end_loc)
      else
        locp => locd
      end if
      call gatherv_values(locp,inumelems(my_cart_id+1),inumelems(my_cart_id+1),num_compute,&
           tot, numtot, size(tot,2), size(loc,1), ireceiver, imp_integers, icomm_cart,&
           ymsg, ierrstat, init=.false.)
      IF (ierrstat /= 0) CALL model_abort(my_cart_id,99904, 'MPI_GATHERV', &
           'gather_radv', ierrstat)
    END IF
  END SUBROUTINE loc_gatherv_integers2d

  SUBROUTINE loc_gatherv_reals (loc,tot)
    REAL (KIND=wp), INTENT(IN)   , POINTER :: loc(:)
    REAL (KIND=wp), INTENT(INOUT), POINTER :: tot(:)
    REAL (KIND=wp), POINTER :: locp(:)
    REAL (KIND=wp), TARGET  :: locd(1)
    INTEGER (KIND=iintegers) :: s_loc, numtot
    IF (ASSOCIATED(loc)) THEN
      if (size(loc) > 0) then
        locp => loc(start_loc:end_loc)
      else
        locp => locd
      end if
      call gatherv_values(locp,inumelems(my_cart_id+1),inumelems(my_cart_id+1),num_compute,&
           tot, numtot, size(tot), ireceiver, imp_reals, icomm_cart,&
           ymsg, ierrstat, init=.false.)
      IF (ierrstat /= 0) CALL model_abort(my_cart_id,99904, 'MPI_GATHERV', &
           'gather_radv', ierrstat)
    END IF
  END SUBROUTINE loc_gatherv_reals

  SUBROUTINE loc_gatherv_reals2d (loc,tot)
    REAL (KIND=wp), INTENT(IN)   , POINTER :: loc(:,:)
    REAL (KIND=wp), INTENT(INOUT), POINTER :: tot(:,:)
    REAL (KIND=wp), POINTER :: locp(:,:)
    REAL (KIND=wp), TARGET  :: locd(1,1)
    INTEGER (KIND=iintegers) :: s_loc, numtot, i
    IF (ASSOCIATED(loc)) THEN
      if (size(loc) > 0) then
        locp => loc(:,start_loc:end_loc)
      else
        locp => locd
      end if
      call gatherv_values(locp,inumelems(my_cart_id+1),inumelems(my_cart_id+1),num_compute,&
           tot, numtot, size(tot,2), SIZE(loc,1), ireceiver, imp_reals, icomm_cart,&
           ymsg, ierrstat, init=.false.)
      IF (ierrstat /= 0) CALL model_abort(my_cart_id,99904, 'MPI_GATHERV', &
           'gather_radv', ierrstat)
    END IF
  END SUBROUTINE loc_gatherv_reals2d

  SUBROUTINE loc_gatherv_reals3d (loc,tot)
    REAL (KIND=wp), INTENT(IN)   , POINTER :: loc(:,:,:)
    REAL (KIND=wp), INTENT(INOUT), POINTER :: tot(:,:,:)
    REAL (KIND=wp), POINTER :: locp(:,:,:)
    REAL (KIND=wp), TARGET  :: locd(1,1,1)
    INTEGER (KIND=iintegers) :: s_loc, numtot, i
    IF (ASSOCIATED(loc)) THEN
      if (size(loc) > 0) then
        locp => loc(:,:,start_loc:end_loc)
      else
        locp => locd
      end if
      call gatherv_values(locp,inumelems(my_cart_id+1),inumelems(my_cart_id+1),num_compute,&
           tot, numtot, size(tot,3), SIZE(loc,1), SIZE(loc,2), ireceiver, imp_reals, icomm_cart,&
           ymsg, ierrstat, init=.false.)
      IF (ierrstat /= 0) CALL model_abort(my_cart_id,99904, 'MPI_GATHERV', &
           'gather_radv', ierrstat)
    END IF
  END SUBROUTINE loc_gatherv_reals3d

  SUBROUTINE loc_gatherv_logical2d (loc,tot)
    LOGICAL, INTENT(IN)   , POINTER :: loc(:,:)
    LOGICAL, INTENT(INOUT), POINTER :: tot(:,:)
    LOGICAL, POINTER :: locp(:,:)
    LOGICAL, TARGET  :: locd(1,1)
    INTEGER (KIND=iintegers) :: s_loc, numtot, i
    IF (ASSOCIATED(loc)) THEN
      if (size(loc) > 0) then
        locp => loc(:,start_loc:end_loc)
      else
        locp => locd
      end if
      call gatherv_values(locp,inumelems(my_cart_id+1),inumelems(my_cart_id+1),num_compute,&
           tot, numtot, size(tot,2), SIZE(loc,1), ireceiver, imp_logical, icomm_cart,&
           ymsg, ierrstat, init=.false.)
      IF (ierrstat /= 0) CALL model_abort(my_cart_id,99904, 'MPI_GATHERV', &
           'gather_radv', ierrstat)
    END IF
  END SUBROUTINE loc_gatherv_logical2d

END SUBROUTINE gather_radv

!==============================================================================
!+ Procedure in "src_obs_sat" for broadcasting the t_rad_set type
!------------------------------------------------------------------------------

SUBROUTINE p_bcast_rad_set(buffer,source,comm)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine broadcasts a t_Rad_set structure to all PEs
!
!------------------------------------------------------------------------------
! Parameters
  TYPE(t_rad_set),   INTENT(INOUT) :: buffer(:)
  INTEGER,           INTENT(IN)    :: source
  INTEGER, OPTIONAL, INTENT(INOUT) :: comm

! Local variables
  INTEGER :: iradset, real_type
  INTEGER :: lcom, errorcode
!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE p_bcast_rad_set
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Section 1: Broadcast derived type 
!-------------------------------------------------------------------------------
 
  lcom = MPI_COMM_WORLD ; IF(PRESENT(comm)) lcom = comm

  CALL MPI_Bcast(buffer,SIZE(buffer) * SIZE(TRANSFER(buffer(1),(/' '/))), & 
                 MPI_BYTE, source,lcom,errorcode)

  IF (errorcode /= MPI_SUCCESS) THEN
    PRINT *, 'MPI ERROR in MPI_Bcast: ', errorcode
    STOP 'MPI ERROR'
  ENDIF

!-------------------------------------------------------------------------------
!- Section 2: Broadcast all allocated arrays 
!-------------------------------------------------------------------------------
  DO iradset=1,SIZE(buffer)
    IF (ASSOCIATED(buffer(iradset)%chan)) THEN
      IF (my_world_id /= source) &
        ALLOCATE(buffer(iradset)%chan(buffer(iradset)%n_chan))

      CALL MPI_Bcast(buffer(iradset)%chan, buffer(iradset)%n_chan, &
                     imp_integers, source, lcom,errorcode)
    ENDIF

    IF (ASSOCIATED(buffer(iradset)%ichan)) THEN
      IF (my_world_id /= source) &
        ALLOCATE(buffer(iradset)%ichan(buffer(iradset)%n_chan))

      CALL MPI_Bcast(buffer(iradset)%ichan, buffer(iradset)%n_chan, &
                     imp_integers, source, lcom,errorcode)
    ENDIF

    IF (ASSOCIATED(buffer(iradset)%band)) THEN
      IF (my_world_id /= source) &
        ALLOCATE(buffer(iradset)%band(buffer(iradset)%n_chan))

        CALL MPI_Bcast(buffer(iradset)%band, buffer(iradset)%n_chan, &
                      imp_integers, source, lcom,errorcode)
      ENDIF

    IF (ASSOCIATED(buffer(iradset)%flag)) THEN
      IF (my_world_id /= source) &
        ALLOCATE(buffer(iradset)%flag(buffer(iradset)%n_chan))

      CALL MPI_Bcast(buffer(iradset)%flag, buffer(iradset)%n_chan, &
                     imp_integers, source, lcom,errorcode)
      ENDIF

    IF (ASSOCIATED(buffer(iradset)%var)) THEN
      IF (my_world_id /= source) &
        ALLOCATE(buffer(iradset)%var(buffer(iradset)%n_chan))

      IF (KIND(buffer(1)%var) == 4) THEN
        real_type = MPI_REAL4
      ELSE
        real_type = MPI_DOUBLE_PRECISION
      ENDIF

      CALL MPI_Bcast(buffer(iradset)%var, buffer(iradset)%n_chan, &
                     real_type, source, lcom,errorcode)
    ENDIF
   
    IF (associated(buffer(iradset)%sensor_instr)) THEN
      IF (my_world_id /= source) &
        ALLOCATE(buffer(iradset)%sensor_instr(buffer(iradset)%n_sens))
      
      CALL MPI_Bcast(buffer(iradset)%sensor_instr, buffer(iradset)%n_sens, &
                     imp_integers, source, lcom,errorcode)
    ENDIF

    IF (associated(buffer(iradset)%oe_par)) THEN
      IF (my_world_id /= source) &
        ALLOCATE(buffer(iradset)%oe_par(buffer(iradset)%n_chan))

      call p_bcast(buffer(iradset)%oe_par, source, comm)
    ENDIF
  ENDDO

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE p_bcast_rad_set
!==============================================================================

!==============================================================================
!+ Procedure in "src_obs_sat" for broadcasting the t_obserr_par type
!------------------------------------------------------------------------------

SUBROUTINE p_bcast_obserr_par(buffer,source,comm)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine broadcasts a t_obserr_par structure to all PEs
!
!------------------------------------------------------------------------------
! Parameters
  TYPE(t_obserr_par),   INTENT(INOUT) :: buffer(:)
  INTEGER,              INTENT(IN)    :: source
  INTEGER, OPTIONAL,    INTENT(INOUT) :: comm

! Local variables
  INTEGER :: i, real_type
  INTEGER :: lcom, errorcode
!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE p_bcast_obserr_par
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Section 1: Broadcast derived type 
!-------------------------------------------------------------------------------

  lcom = MPI_COMM_WORLD ; IF(PRESENT(comm)) lcom = comm

  CALL MPI_Bcast(buffer,SIZE(buffer) * SIZE(TRANSFER(buffer(1),(/' '/))), &
                 MPI_BYTE, source,lcom,errorcode)

  IF (errorcode /= MPI_SUCCESS) THEN
    PRINT *, 'MPI ERROR in MPI_Bcast: ', errorcode
    STOP 'MPI ERROR'
  ENDIF

!-------------------------------------------------------------------------------
!- Section 2: Broadcast all allocated arrays 
!-------------------------------------------------------------------------------
  DO i=1,SIZE(buffer)
    IF (ASSOCIATED(buffer(i)%par)) THEN
      IF (my_world_id /= source) &
        ALLOCATE(buffer(i)%par(buffer(i)%n))

      CALL MPI_Bcast(buffer(i)%par, buffer(i)%n, &
                     imp_reals, source, lcom,errorcode)
    ENDIF
  END DO
END SUBROUTINE p_bcast_obserr_par

!==============================================================================
!+ Procedure in "src_obs_sat" for updating the state of an observable
!------------------------------------------------------------------------------
SUBROUTINE change_use_rpt_rad(radv,ichan,irpt,flag)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine updates the state of an report in a radv structure
!
!------------------------------------------------------------------------------

! Parameters
  TYPE(t_radv), INTENT(INOUT) :: &
    radv ! input radv structure
  INTEGER (KIND=iintegers), INTENT(IN) :: &
    ichan, & ! channel
    irpt,  & ! report index
    flag     ! the flag of the check triggered
! Local Variabes
  INTEGER (KIND=iintegers) :: &
    state   ! the new state
    
!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE change_use_rpt_rad
!-------------------------------------------------------------------------------
  state = map_flg2state(flag)

  IF (radv%state(ichan,irpt) < state) THEN
    radv%state(ichan,irpt) = state
  ENDIF

  radv%flags(ichan,irpt) = IBSET(radv%flags(ichan,irpt),flag)
!-------------------------------------------------------------------------------
!- End SUBROUTINE change_use_rpt_rad
!-------------------------------------------------------------------------------

END SUBROUTINE change_use_rpt_rad
!==============================================================================
#endif

!==============================================================================
!+ Procedure in "src_obs_sat" for the input of NAMELIST obs_radctl
!------------------------------------------------------------------------------

SUBROUTINE input_obs_radctl (nuspecif, nuin, ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine organizes the input of the NAMELIST file for satellite
!   observations processing. This file contains multiple NAMELIST groups 
!   TOVS_OBS_CHAN describing a single satellite dataset.
!
! Method:
!
!------------------------------------------------------------------------------

! Parameters
  INTEGER   (KIND=iintegers),   INTENT (IN)    :: &
    nuspecif,                                     & ! Unit number for protocolling the task
    nuin                                            ! Unit number for Namelist INPUT file

  INTEGER   (KIND=iintegers),   INTENT (INOUT) :: &
    ierrstat                                        ! error status variable
! Local variables
#ifdef RTTOV10
  TYPE      (t_rad_set) ::                        &
    tmp                                             ! Temporary buffer for a rad_set

  INTEGER   (KIND=iintegers) ::                   &
    iradset,                                      & ! counting variable for radiance set
    istat,                                        & ! result of calls to mo_rad subroutines
    i,                                            & ! counter variable
    izerr, ierr                                     ! error status variables

  ! Definition of the namelist groups
  NAMELIST /TOVS_OBS/ input_path, satpp_files, usesurftype,                  &
                      input_path_hrit, hrit_name, hrit_chan, hrit_hh,        &
                      calc_k, thin_line, thin_col, nprof_rttov, monitor_prof,&
                      hrit_superob, rad_output
#endif  
    
!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE input_obs_radctl
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Section 2: Read in namelist group TOVS_OBS and all following groups i
!-            TOVS_OBS_CHANs
!-------------------------------------------------------------------------------
#ifdef RTTOV10
IF (my_world_id == 0) THEN
  iradset  = 0
  istat    = 0

  read(nuin,nml=TOVS_OBS, iostat=istat)

  IF (istat /= 0) RETURN

  rad_output=ibset(rad_output,OUT_ASB)

  DO WHILE (istat == 0)
    call read_tovs_obs_chan_nml(nuin, tmp, istat, iradset + 1)

    IF (istat /= 0) EXIT

    IF (iradset + 1 > SIZE(rad_set)) THEN
      ierrstat = 99002
      PRINT *, ' ERROR    *** Maximum number of radiance sets reached. *** '
      EXIT
    ENDIF  

    iradset = iradset + 1
    rad_set(iradset) = tmp 
  ENDDO
ENDIF

!-------------------------------------------------------------------------------
!- Section 3: Evaluate parameters
!-------------------------------------------------------------------------------

i = LEN_TRIM(input_path)
IF (i > 0) THEN
  IF (input_path(i:i) /= '/') input_path = TRIM(input_path) // '/'
ENDIF
i = LEN_TRIM(input_path_hrit)
IF (i > 0) THEN
  IF (input_path_hrit(i:i) /= '/') input_path_hrit = TRIM(input_path_hrit) // '/'
ENDIF

!-------------------------------------------------------------------------------
!- Section 4: Distribute variables to all nodes
!-------------------------------------------------------------------------------

izerr = 0

IF (nproc > 1) THEN
  call distribute_values(input_path, len(input_path), 0, &
       imp_character, icomm_world,ierr)
  izerr = izerr + ierr
  call distribute_values(input_path_hrit, len(input_path_hrit), 0, &
       imp_character, icomm_world)
  izerr = izerr + ierr
  call distribute_values(usesurftype, 1, 0, imp_integers, icomm_world)
  izerr = izerr + ierr
  DO i = 1,size(satpp_files)
    call distribute_values(satpp_files(i), len(satpp_files(i)), 0, &
         imp_character, icomm_world)
    izerr = izerr + ierr
  ENDDO
  call distribute_values(input_path_hrit, len(input_path_hrit), 0, &
       imp_character, icomm_world)
  izerr = izerr + ierr
  CALL distribute_values(hrit_hh, SIZE(hrit_hh), 0, imp_reals, icomm_world)
  izerr = izerr + ierr
  CALL distribute_values(hrit_name,LEN(hrit_name),0,imp_character, icomm_world)
  izerr = izerr + ierr
  DO i= 1, SIZE(hrit_chan)
    CALL distribute_values(hrit_chan(i), LEN(hrit_chan(i)), 0, imp_character, icomm_world)
    izerr = izerr + ierr
  ENDDO
  call distribute_values(calc_k, 1, 0, imp_logical, icomm_world)
  izerr = izerr + ierr
  call distribute_values(thin_line, 1, 0, imp_integers, icomm_world)
  izerr = izerr + ierr
  call distribute_values(thin_col , 1, 0, imp_integers, icomm_world)
  izerr = izerr + ierr
  call distribute_values(hrit_superob, 1, 0, imp_integers, icomm_world)
  izerr = izerr + ierr
  call distribute_values(rad_output, 1, 0, imp_integers, icomm_world)
  izerr = izerr + ierr
  call distribute_values(nprof_rttov, 1, 0, imp_integers, icomm_world)
  izerr = izerr + ierr
  call distribute_values(monitor_prof, 1, 0, imp_integers, icomm_world)
  izerr = izerr + ierr
  IF (izerr /= 0) THEN
    PRINT*, ' distributing namelist information to nodes failed '
    ierrstat = -1
    RETURN
  ENDIF

  call p_bcast(rad_set, 0,icomm_world)

ENDIF

!-------------------------------------------------------------------------------
!- Section 4: Output of the namelist variables
!-------------------------------------------------------------------------------

IF (my_world_id == 0) THEN
  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(A28)') '0  NAMELIST tovs_obs'
  WRITE (nuspecif, '(A28)') '   -----------------'
  WRITE (nuspecif, '(A2)')  '  '

  WRITE (nuspecif, '(T7,A,T24,A,T79,A)')     'Variable', 'Value', 'Format'
  WRITE (nuspecif, '(T7,A,T24,A,T79,A)')     'input_path', trim(input_path), 'A'
  WRITE (nuspecif, '(T7,A,T24,A,T79,A)')     'input_path_hrit', trim(input_path_hrit), 'A'
  WRITE (nuspecif, '(T7,A,T24,L1,T79,A)')    'calc_k', calc_k, 'L'
  WRITE (nuspecif, '(T7,A,T24,I3,T79,A)')    'thin_line', thin_line, 'I'
  WRITE (nuspecif, '(T7,A,T24,I3,T79,A)')    'thin_col', thin_col, 'I'
  WRITE (nuspecif, '(T7,A,T24,A,T79,A)')     'hrit_name', hrit_name, 'A'
  DO i = 1,COUNT(hrit_chan /= '')
    WRITE (nuspecif, '(T7,A,T24,A,T79,A)')   'hrit channel', hrit_chan(i),'A'
  ENDDO
  WRITE (nuspecif, '(T7,A,T24,3F5.1,T79,A)') 'hrit_hh', hrit_hh(:), 'E*3'
  WRITE (nuspecif, '(T7,A,T24,I3,T79,A)')    'hrit_superob', hrit_superob, 'I'
  WRITE (nuspecif, '(T7,A,T24,I8,T79,A)')    'nprof_rttov', nprof_rttov, 'I'
  WRITE (nuspecif, '(T7,A,T24,I3,T79,A)')    'monitor_prof', monitor_prof, 'I'
  WRITE (nuspecif, '(T7,A,T24,I3,T79,A)')    'rad_output', rad_output, 'I'

  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(T7,A)')  'Input Files:'
  WRITE (nuspecif, '(A2)')  '  '
  DO i = 1,COUNT(satpp_files /= '')
    WRITE (nuspecif, '(T7,A,T24,A)') 'SATPP', trim(satpp_files(i))
  ENDDO
  
  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(T7,A)')  'TOVS_OBS_CHAN'
  WRITE (nuspecif, '(A2)')  '  '

  WRITE (nuspecif, '(T7,A,T24,A,T39,A)')  'Variable', 'Value', 'Format'
  DO i = 1,iradset
    WRITE (nuspecif, '(T7,A,T24,A,T39,A)')  'c_satellite', &
                                            sat_id2name(rad_set(i)%satid), &
                                            'A'
    WRITE (nuspecif, '(T7,A,T24,I8.2,T39,A)')  'grid', rad_set(i)%grid, 'I'
  ENDDO
ENDIF

l_H = calc_k .and. (IAND(monitor_prof, MP_H) > 0)

#endif

ierrstat = 0

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE input_obs_radctl
!==============================================================================

!==============================================================================
!+ Procedure in "src_obs_sat" for the input of satpp files
!------------------------------------------------------------------------------

SUBROUTINE input_obs_satpp (ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine organizes the input of the satpp files
!
!------------------------------------------------------------------------------

! Parameters
  INTEGER (KIND=iintegers), INTENT (INOUT) :: &
    ierrstat        ! error status variable
    
! Local variables
  INTEGER   (KIND=iintegers) ::               &
    nhrit,                                    & ! number of hrit files
    nsat                                        ! number of satellite datasets
#ifdef RTTOV10    
  CHARACTER(LEN=256) ::                       &
    ymsg,                                     & ! error message
    path                                        ! path to input file

  CHARACTER(LEN=64) :: yversion         
  CHARACTER(LEN=45) :: yfdbk_descript
  CHARACTER(LEN=12) :: yveri_ref_datim

  INTEGER   (KIND=iintegers) ::               &
    istat,                                    & ! status variable
    isat,                                     & ! couting variable for satpp files
    i, cc,                                    & ! loop index
    nrec,                                     & ! number of satpp records to read in
    irec,                                     & ! counting variable for record
    inobs,                                    & ! counting variable for obs
    ioffs,                                    & ! 
    ichan, imissing,                          & ! counting variable for channel
    iinstr,                                   & ! counting variable for channel
    icnt, istart,iend,                        & ! counting variable
    idiffm,                                   & ! time difference [minutes]
    imoyy, imomm, imodd, imohh,               & ! Start date of forecast               
    isthght,                                  & ! Station Height
    inobs_node_loc(num_compute),              & ! number of obs per node (local)
    inobs_dropped_area,                       & ! no of obs dropped for area 
    inobs_dropped_tstart,                     & ! no of obs dropped for too early 
    inobs_dropped_tend,                       & ! no of obs dropped for too late 
    inobs_dropped,                            & ! no of obs dropped
    iveri_ref_date,                           & ! verification reference date
    iveri_ref_time,                           & ! verification reference date
    iveri_start,                              & ! verification reference date
    iveri_end,                                & ! verification reference date
    domain(3),                                & ! 
    iveri_run_type,                           & ! 
    iveri_run_class,                          & ! 
    iveri_exp_id,                             & ! 
    iveri_forecast_time,                      & ! 
    iveri_epsmem,                             & ! 
    varid,                                    & ! 
    satstep,                                  & ! time step of satellite input
    n_add,                                    & ! number of additional levels
    nb_t,                                     & ! superobbing time box size
    nb2_t,                                    & ! superobbing half time box size
    mn_t,                                     & ! superobbing min time step
    mx_t,                                     & ! superobbing max time step
    id_t,                                     & ! superobbing step size of time step
    nb_c,                                     & ! superobbing column box size
    nb2_c,                                    & ! superobbing half column box size
    mn_c,                                     & ! superobbing min column
    mx_c,                                     & ! superobbing max column
    id_c,                                     & ! superobbing step size of column
    nb_l,                                     & ! superobbing line box size
    nb2_l,                                    & ! superobbing line half box size
    mn_l,                                     & ! superobbing min line
    mx_l,                                     & ! superobbing max line
    id_l,                                     & ! superobbing step size of line
    mso,                                      & ! superobbing box size (work)
    it,                                       & ! superobbing time step
    il,                                       & ! superobbing line
    ic,                                       & ! superobbing column
    ibox                                        ! superobbing box number
  INTEGER (KIND=iintegers), ALLOCATABLE ::    &
    ifcaststep(:),                            & ! timesteps assigned to obs
    iob_tot(:),                               & ! zonal idx of gridpt 
    job_tot(:),                               & ! meridional idx of gridpt
    iindex(:),                                & ! index vector of rads to use
    isort_index(:),                           & ! index vector of rads to use
    inode_index(:),                           & ! index vector of rads to use
    iline(:),                                 & ! scanline of FOV (for superobbing)
    icol (:),                                 & ! FOV number (for superobbing)
    istep(:),                                 & ! Time step
    iso(:),                                   & ! working variable for superobbing
    ibox_so(:)                                  ! Superobbing box ID
  INTEGER (KIND=iintegers), POINTER ::        &  
    idummy(:)                                   ! Array of channels
  REAL      (KIND=wp)     ::              &
    rlon,rlat,                                & ! Rotated Latitude/Longitude
    zio_tot,zjo_tot,                          & ! obs location in grid points
    rdummy                                      ! dummy variable
  REAL :: & 
    pole(2),lower_left(2),upper_right(2),     & ! 
    resolution(2)                               ! 
  LOGICAL :: &
    valid,                                    & ! result of satpp reading routine
    lsea,                                     & ! If the obs is above sea
    ldummy
  TYPE (t_radv), TARGET  ::                   &
    radv                                        ! Radiances from satpp
  TYPE (t_radv) ::                            &
    obs_rad_loc                                 ! Radiances for COSMO
  TYPE (t_rad_set), POINTER  ::               &
    rs              
  TYPE (t_radv), POINTER ::                   &
    r                                           ! Radiances from satpp
  TYPE (t_sensor), POINTER  ::                &
    sensor

  REAL(KIND=wp) ::                        & ! coordinates for reading HRIT satellite files
    endlat_tot,                               &
    endlon_tot,                               &
    shh                                         ! hour index

  !for reading NWCSAF SEVIRI cloud type
  CHARACTER (LEN=14)         :: yzactdate
  CHARACTER (LEN=28)         :: yzdum           ! actual date in the form   wd   dd.mm.yy  hh UTC

  REAL (KIND= wp)        ::               &
    rzdum

  INTEGER (KIND=iintegers)   ::               &
    izdum,                                    & ! dummy variable
    nhrit_chan,                               & ! number of channels
    nstep                                       ! number of hrit images

!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE input_obs_satpp
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Section 1: Initialize variables 
!-------------------------------------------------------------------------------

  if (hrit_name /= '') then
    ! determine HRIT-file list, i.e. which files to open
    nhrit_chan   = COUNT(hrit_chan /= '')

    cc=0
    IF (hrit_hh(3) > 0.0_wp) THEN
      nstep = NINT( (hrit_hh(2)-hrit_hh(1)) / hrit_hh(3) ) + 1
    ELSE
      nstep = 0
    END IF
    DO i = 1, nstep
      shh = hrit_hh(1) + (i-1) * hrit_hh(3)
      satstep = INT((shh*3600.0_wp)/dt,iintegers)
      CALL get_utc_date ( satstep, ydate_ini, dt, itype_calendar, yzactdate       &
           , yzdum, izdum,rzdum)
      DO ichan = 1,nhrit_chan
        cc = cc +1
        hrit_files(cc)='H-000-'//hrit_name//'__-'//hrit_name//'________-'//hrit_chan(ichan)//'___-000007___-'//  &
             yzactdate(1:12)//'-__'
        cc = cc +1
        hrit_files(cc)='H-000-'//hrit_name//'__-'//hrit_name//'________-'//hrit_chan(ichan)//'___-000008___-'//  &
             yzactdate(1:12)//'-__'
      ENDDO
      cc = cc +1
      hrit_files(cc)='H-000-'//hrit_name//'__-'//hrit_name//'________-_________-PRO______-'//  &
           yzactdate(1:12)//'-__'
    ENDDO
  end if

  nsatpp   = COUNT(satpp_files /= '')
  nsat     = nsatpp
  nhrit    = COUNT(hrit_files /= '')
  IF (nhrit > 0) nsat = nsat + 1
  nsensors = 0

  ALLOCATE (obs_rad(nsat),                &
            nobspe_gath(num_compute,nsat),&
            nobspe_tot(num_compute,nsat), &
            stat = istat)

  IF (my_cart_id == 0 .AND. istat == 0) THEN
    ALLOCATE(fdbk(nsat),        & 
             fdbk_offset(nsat), &
              stat = istat)
  ENDIF

  IF (istat /= 0) THEN
    ierrstat = 99101
    PRINT *, ' ERROR   *** Allocating memory failed. ***'
    RETURN
  ENDIF  
 
  READ (ydate_ini,'(I4,I2,I2,I2)') imoyy, imomm, imodd, imohh
 
  IF (my_cart_id == 0) THEN
    call get_free_unit(nurad)
    OPEN(nurad, FILE=yurad, FORM='FORMATTED', STATUS='UNKNOWN', &
                POSITION='APPEND', IOSTAT=istat)

    IF (istat /= 0) THEN
      ierrstat = 99102
      PRINT *, ' ERROR   *** OPENING OF FILE yurad FAILED. ***'
      RETURN
    ENDIF  
    
    WRITE (nurad,*) '  '
    WRITE (nurad,'(T2,A)') 'SATELLITE DATA Input'
    WRITE (nurad,'(T2,A)') '--------------------'
  ENDIF  
  
  ! compute edges in geographical coordinates needed for reading HRIT files
  endlat_tot= startlat_tot + dlat * je_tot
  endlon_tot= startlon_tot + dlon * ie_tot

  maxlat = MAX(phirot2phi(endlat_tot,startlon_tot,pollat,pollon, polgam), &
       phirot2phi(endlat_tot,endlon_tot,pollat,pollon, polgam))
  minlat = MIN(phirot2phi(startlat_tot, startlon_tot, pollat,pollon,polgam), &
       phirot2phi(startlat_tot, endlon_tot, pollat,pollon,polgam))

  maxlon = MAX(rlarot2rla(endlat_tot,endlon_tot,pollat,pollon, polgam),  &
       rlarot2rla(startlat_tot,endlon_tot,pollat,pollon, polgam))
  minlon = MIN(rlarot2rla(startlat_tot, startlon_tot, pollat,pollon,polgam), &
       rlarot2rla(endlat_tot, startlon_tot, pollat,pollon,polgam))

  DO isat=1,nsat
!-------------------------------------------------------------------------------
!- Section 2: Read in the sat file
!-------------------------------------------------------------------------------
    IF (isat <= nsatpp) THEN
      !-------------------------------------------------------------------------------
      !- Section 2a: Read in the satpp file
      !-------------------------------------------------------------------------------

      path = TRIM(input_path) // TRIM(satpp_files(isat))

      ! find out about the size of the satpp
      call read_satpp_feedbk(TRIM(path),radv,valid,lprint=.False., lread=.False.)

      IF (.NOT. valid) THEN
        ierrstat = 99103
        PRINT *, ' ERROR   *** Inquiring length of satpp failed. ***'
        RETURN
      ENDIF

      ! The satpp file will be read in parallel on all processors
      ! Therefore each processor has a part of the file in its memory
      istart = nint(radv%n_rec* my_cart_id   /(1.*num_compute))+1
      iend   = nint(radv%n_rec*(my_cart_id+1)/(1.*num_compute))
      nrec = istart - iend + 1
      call read_satpp_feedbk(TRIM(path),radv,valid, &
                             istart = istart,       &
                             iend   = iend,         &
                             lprint = .false.)

      IF (.NOT. valid) THEN
        ierrstat = 99104
        WRITE (stderr,*) ' ERROR   *** Reading satpp failed. ***'
        RETURN
      ENDIF

      path = radv % filename

    ELSE
      !-------------------------------------------------------------------------------
      !- Section 2a: Read in data from HRIT files
      !-------------------------------------------------------------------------------
      DO i = 1, nhrit
        hrit_files(i) = TRIM(input_path_hrit)//TRIM(hrit_files(i))
      END DO
      call read_hrit(hrit_files(1:nhrit), radv, lread=.false.,                  &
           minlat=real(minlat), maxlat=real(maxlat), minlon=real(minlon),       &
           maxlon=real(maxlon), zenmax=real(zenmax),                            &
           thin_line=thin_line, thin_col=thin_col,                              &
           lprint=(my_cart_id == 0), unit=nurad)

      ! The hrit files will be read in parallel on all processors
      ! Therefore each processor has a part of the file in its memory
      istart = nint(radv%n_rec* my_cart_id   /(1.*num_compute))+1
      iend   = nint(radv%n_rec*(my_cart_id+1)/(1.*num_compute))
      nrec = istart - iend + 1
      call read_hrit(hrit_files(1:nhrit), radv, lread=.true.,                       &
                     istart = nint(radv%n_rec* my_cart_id   /(1.*num_compute))+1,   &
                     iend   = nint(radv%n_rec*(my_cart_id+1)/(1.*num_compute))  ,   &
                     minlat=real(minlat), maxlat=real(maxlat), minlon=real(minlon), &
                     maxlon=real(maxlon), zenmax=real(zenmax), time_mode='n',       &
                     thin_line=thin_line, thin_col=thin_col,                        &
                     lprint=(my_cart_id == 0), unit=nurad)
      
      IF (radv%i%instr(1) == 21) THEN
        ! For SEVIRI the HRIT library yields channel numbers from 1-12. i.e. all 
        ! available channels. However, in RTTOV the first 3 and the last channel are
        ! not available. Thus 3 has to be subtracted from the original channel number 
        ! to obtain the RTTOV channel number.
        radv%i%chan(:) = radv%i%chan(:) - 3
        IF (any(radv%i%chan(:) < 1 .or. radv%i%chan(:) > 8)) THEN
          ierrstat = 99105
          PRINT *, ' ERROR   *** Invalid channel number for SEVIRI (read from HRIT). ***'
          RETURN
        END IF
      END IF

      path = '(HRIT)'

    END IF

    nrec = radv%n_rec
    CALL global_values(nrec,1,'SUM',imp_integers,icomm_cart,0,ymsg,istat)
    
    IF (my_cart_id == 0) THEN
      WRITE (nurad,*) '  '
      WRITE (nurad,'(T4,A,T20,A)') 'Radiance input', trim(path)
      WRITE (nurad,'(T6,A,T30,I15)') 'Number Of Records: ', nrec
      ! flush YURAD file
      CLOSE (nurad)
      OPEN(nurad, FILE=yurad, FORM='FORMATTED', STATUS='OLD', &
           POSITION='APPEND', IOSTAT=istat)      
    ENDIF    

!-------------------------------------------------------------------------------
!- Section 3a: Link Radiance data to namelist information 
!-------------------------------------------------------------------------------
    if (radv%n_rec > 0) CALL link_rad(radv,istat)

    IF (istat /= 0) THEN
      ierrstat = 99106
      WRITE (stderr,*) &
        ' ERROR   *** Associating TOVS_OBS_CHAN namelists with input files failed ***'
      RETURN
    ENDIF
  
    rs => radv%i
!-------------------------------------------------------------------------------
!- Section 3b: Prepare additional fields 
!-------------------------------------------------------------------------------
    ALLOCATE(radv%state(rs%n_chan,radv%n_rec), &
             radv%flags(rs%n_chan,radv%n_rec), &
!             radv%check(rs%n_chan,radv%n_rec), &
             stat = istat)

    IF (istat /= 0) THEN
      ierrstat = 99107
      PRINT *, ' ERROR   *** Allocating memory failed. ***'
      RETURN
    ENDIF  

    DO irec=1,radv%n_rec
      DO ichan=1,rs%n_chan
        radv%state(ichan,irec) = ST_ACTIVE
        radv%flags(ichan,irec) = 0

        ! set initial state
        IF (.NOT. radv%valid(ichan,irec)) THEN
          CALL change_use_rpt_rad(radv,ichan,irec,FL_DATASET)
        ELSE IF (BTEST(rs%flag(ichan), USE_PASSIVE)) THEN
          CALL change_use_rpt_rad(radv,ichan,irec,FL_OBSTYPE)
        ENDIF
      ENDDO
    ENDDO  
    
!-------------------------------------------------------------------------------
!- Section 3c: Evaluate obs location and time 
!-------------------------------------------------------------------------------
    ALLOCATE(ifcaststep(radv%n_rec), iob_tot(radv%n_rec),          &
             job_tot(radv%n_rec), iindex(radv%n_rec+1),            &
             isort_index(radv%n_rec+1), inode_index(radv%n_rec+1), &
             stat = istat)

    IF (istat /= 0) THEN
      ierrstat = 99108
      PRINT *, ' ERROR   *** Allocating memory failed. ***'
      RETURN
    ENDIF  

    ! The following is really weird, we should get rid of this
    ! slow function, somehow
    DO irec = 1,radv%n_rec
      CALL diff_minutes(imoyy,imomm,imodd,imohh,0,                 &
                        radv % date(irec)     / 10000     , &
                        MOD(radv % date(irec) / 100  ,100), &
                        MOD(radv % date(irec)        ,100), &
                        radv % time(irec)     / 10000     , &
                        MOD(radv % time(irec) / 100  ,100), &
                        itype_calendar, idiffm,             &
                        istat                                        )

      IF (istat /= 0) THEN
        ifcaststep(irec) = -1
      ELSE  
        ifcaststep(irec) = int(real(idiffm) * 60.0_wp / dt)
      ENDIF
    ENDDO
    
    DO irec = 1,radv%n_rec
      rlon = rla2rlarot(radv%dlat(irec), radv%dlon(irec), &
                        pollat, pollon, polgam)

      rlat = phi2phirot(radv%dlat(irec), radv%dlon(irec), &
                        pollat, pollon)

      zio_tot = 1._wp + (rlon - startlon_tot) / dlon
      zjo_tot = 1._wp + (rlat - startlat_tot) / dlat
     
      ! the following was stolen from func obs_assign_gridpt
      rdummy = MAX(nolbc, 1+nboundlines) + 1.E-8_wp

      IF (zio_tot > rdummy .AND. zio_tot < ie_tot + 1.0_wp - rdummy .AND. &
          zjo_tot > rdummy .AND. zjo_tot < je_tot + 1.0_wp - rdummy) THEN
          
        iob_tot(irec) = NINT(zio_tot)
        job_tot(irec) = NINT(zjo_tot)

        iob_tot(irec) = MAX(iob_tot(irec), 1 + nolbc)
        iob_tot(irec) = MIN(iob_tot(irec), ie_tot - nolbc)
        job_tot(irec) = MAX(job_tot(irec), 1 + nolbc)
        job_tot(irec) = MIN(job_tot(irec), je_tot - nolbc)
      ELSE
        iob_tot(irec) = 0
        job_tot(irec) = 0
      ENDIF  
    ENDDO  

!-------------------------------------------------------------------------------
!- Section 4: Filter by timestep and location 
!-------------------------------------------------------------------------------
    inobs_dropped_area   = 0
    inobs_dropped_tstart = 0
    inobs_dropped_tend   = 0
    inobs_dropped        = 0
      
    inobs = 0
    DO irec = 1,radv%n_rec
      IF (iob_tot(irec) <= 0 .OR.  job_tot(irec) <= 0) THEN
        DO ichan = 1,rs%n_chan
          CALL change_use_rpt_rad(radv,ichan,irec,FL_AREA)
        ENDDO
        inobs_dropped_area = inobs_dropped_area + 1
      END IF
      IF (ifcaststep(irec) <= nstart) THEN ! Do not assimilate obs before first forecast step
        DO ichan = 1,rs%n_chan
          CALL change_use_rpt_rad(radv,ichan,irec,FL_TIME)
        ENDDO
        inobs_dropped_tstart = inobs_dropped_tstart + 1
      ELSE IF (ifcaststep(irec) > nstop) THEN
        DO ichan = 1,rs%n_chan
          CALL change_use_rpt_rad(radv,ichan,irec,FL_TIME)
        ENDDO
        inobs_dropped_tend = inobs_dropped_tend + 1
      ENDIF

      IF (ANY(radv%state(:,irec) < ST_DISMISS)) THEN
        inobs = inobs + 1    
        iindex(inobs) = irec
      ELSE
        inobs_dropped = inobs_dropped + 1        
      ENDIF    
    ENDDO

    CALL global_values(inobs_dropped_area,1,'SUM',imp_integers,&
                       icomm_cart,0,ymsg,istat)
    CALL global_values(inobs_dropped_tstart,1,'SUM',imp_integers,&
                       icomm_cart,0,ymsg,istat)
    CALL global_values(inobs_dropped_tend,1,'SUM',imp_integers,&
                       icomm_cart,0,ymsg,istat)
    CALL global_values(inobs_dropped,1,'SUM',imp_integers,&
                       icomm_cart,0,ymsg,istat)

    IF (my_cart_id == 0) THEN
      WRITE (nurad,'(T6,A,T30,I15)') '# Obs in wrong area:', inobs_dropped_area
      WRITE (nurad,'(T6,A,T30,I15)') '# Obs too early:',     inobs_dropped_tstart
      WRITE (nurad,'(T6,A,T30,I15)') '# Obs too late:',      inobs_dropped_tend
      WRITE (nurad,'(T6,A,T30,I15)') '# Obs dropped:',       inobs_dropped
      WRITE (nurad,'(T6,A,T30,I15)') '# Obs retained:',      nrec - inobs_dropped
      ! flush YURAD file
      CLOSE (nurad)
      OPEN(nurad, FILE=yurad, FORM='FORMATTED', STATUS='OLD', &
           POSITION='APPEND', IOSTAT=istat)      
    ENDIF  

!-------------------------------------------------------------------------------
!- Section 5: Assign selected obs to nodes 
!-------------------------------------------------------------------------------
    IF (inobs > 0) THEN
      CALL obs_assign_sort_node(SIZE(iob_tot), inobs, iindex(1:inobs+1), &
                                iob_tot,job_tot, num_compute,       &
                                isubpos, nboundlines, my_cart_id,       &
                                inode_index(1:inobs+1),isort_index(1:inobs+1))
      DO icnt=1,num_compute
        inobs_node_loc(icnt) = COUNT(inode_index(1:inobs) == (icnt - 1))
      ENDDO
    ELSE
      inobs_node_loc(:) = 0
    ENDIF  

    IF (SUM(inobs_node_loc) /= inobs) THEN
      ierrstat = 99109
      WRITE (stderr,*) ' ERROR   *** Inconsistent obs numbers. ***'
      RETURN
    ENDIF

    nobspe_tot(:,isat) = inobs_node_loc(:)
    
    CALL global_values(nobspe_tot(:,isat),SIZE(nobspe_tot,1),'SUM',imp_integers,&
                       icomm_cart,0,ymsg,istat)
    IF (istat /= 0) CALL model_abort(my_cart_id,99901, 'global_values', &
                                     'input_obs_satpp', istat)

    CALL distribute_values(nobspe_tot(:,isat),SIZE(nobspe_tot,1), 0, imp_integers, &
                       icomm_cart,istat)
    IF (istat /= 0) CALL model_abort(my_cart_id,99902, 'distribute_values', &
                                     'input_obs_satpp', istat)
    
    IF (my_cart_id == 0) THEN
      WRITE (ymsg,*) '(T6,A,T30,',num_compute,'(1X,I8))'
      WRITE (nurad,'(T6,A,T15,A)') '    Node', 'Assigned obs'
      DO icnt=0,num_compute-1
        WRITE (nurad,'(T6,I8,T15,I12)') icnt,nobspe_tot(icnt+1,isat)
      END DO
      ! flush YURAD file
      CLOSE (nurad)
      OPEN(nurad, FILE=yurad, FORM='FORMATTED', STATUS='OLD', &
           POSITION='APPEND', IOSTAT=istat)      
    ENDIF  

!-------------------------------------------------------------------------------
!- Section 6: Create obs container with obs to be distributed from this node
!-            to other nodes (ordered by nodeid) 
!-------------------------------------------------------------------------------
    call construct(obs_rad_loc, inobs, istat, radv) 

    IF (istat /= 0) THEN
      ierrstat = 99110
      PRINT *, ' ERROR   *** Allocating memory failed. ***'
      RETURN
    ENDIF 

    obs_rad_loc%ntstep(:)  = ifcaststep(isort_index(1:inobs))
    obs_rad_loc%i_box  (:) = iob_tot(isort_index(1:inobs))
    obs_rad_loc%j_box  (:) = job_tot(isort_index(1:inobs))

    IF (ASSOCIATED(radv%obsnum   )) obs_rad_loc%obsnum   (:) = radv%obsnum   (  isort_index(1:inobs))
    IF (ASSOCIATED(radv%date     )) obs_rad_loc%date     (:) = radv%date     (  isort_index(1:inobs))
    IF (ASSOCIATED(radv%time     )) obs_rad_loc%time     (:) = radv%time     (  isort_index(1:inobs))
    IF (ASSOCIATED(radv%date_d   )) obs_rad_loc%date_d   (:) = radv%date_d   (  isort_index(1:inobs))
    IF (ASSOCIATED(radv%time_d   )) obs_rad_loc%time_d   (:) = radv%time_d   (  isort_index(1:inobs))
    IF (ASSOCIATED(radv%dlat     )) obs_rad_loc%dlat     (:) = radv%dlat     (  isort_index(1:inobs))
    IF (ASSOCIATED(radv%dlon     )) obs_rad_loc%dlon     (:) = radv%dlon     (  isort_index(1:inobs))
    IF (ASSOCIATED(radv%fov      )) obs_rad_loc%fov      (:) = radv%fov      (  isort_index(1:inobs))
    IF (ASSOCIATED(radv%scanl    )) obs_rad_loc%scanl    (:) = radv%scanl    (  isort_index(1:inobs))
    IF (ASSOCIATED(radv%stzen    )) obs_rad_loc%stzen    (:) = radv%stzen    (  isort_index(1:inobs))
    IF (ASSOCIATED(radv%stazi    )) obs_rad_loc%stazi    (:) = radv%stazi    (  isort_index(1:inobs))
    IF (ASSOCIATED(radv%sunzen   )) obs_rad_loc%sunzen   (:) = radv%sunzen   (  isort_index(1:inobs))
    IF (ASSOCIATED(radv%sunazi   )) obs_rad_loc%sunazi   (:) = radv%sunazi   (  isort_index(1:inobs))
    IF (ASSOCIATED(radv%landfr   )) obs_rad_loc%landfr (:,:) = radv%landfr   (:,isort_index(1:inobs))
    IF (ASSOCIATED(radv%stype    )) obs_rad_loc%stype  (:,:) = radv%stype    (:,isort_index(1:inobs))
    IF (ASSOCIATED(radv%shgt     )) obs_rad_loc%shgt   (:,:) = radv%shgt     (:,isort_index(1:inobs))
    IF (ASSOCIATED(radv%center   )) obs_rad_loc%center   (:) = radv%center   (  isort_index(1:inobs))
    IF (ASSOCIATED(radv%subcenter)) obs_rad_loc%subcenter(:) = radv%subcenter(  isort_index(1:inobs))
    IF (ASSOCIATED(radv%valid    )) obs_rad_loc%valid  (:,:) = radv%valid    (:,isort_index(1:inobs))
    IF (ASSOCIATED(radv%state    )) obs_rad_loc%state  (:,:) = radv%state    (:,isort_index(1:inobs))
    IF (ASSOCIATED(radv%flags    )) obs_rad_loc%flags  (:,:) = radv%flags    (:,isort_index(1:inobs))
    IF (ASSOCIATED(radv%bt_obs   )) obs_rad_loc%bt_obs (:,:) = radv%bt_obs   (:,isort_index(1:inobs))

!-------------------------------------------------------------------------------
!- Section 7: Merge Obs from all nodes at assigned node 
!-------------------------------------------------------------------------------
    call construct(obs_rad(isat), nobspe_tot(my_cart_id+1,isat), istat, obs_rad_loc)

    IF (istat /= 0) THEN
      ierrstat = 99111
      PRINT *, ' ERROR   *** Allocating memory failed. ***'
      RETURN
    ENDIF 
    
    ioffs = 1
    DO icnt = 1, num_compute
        call gather_radv (obs_rad_loc, ioffs, ioffs + inobs_node_loc(icnt) - 1, &
             obs_rad(isat), 1, irec, icnt-1,istat, nobspe_gath(:,isat))
      ioffs = ioffs + inobs_node_loc(icnt)
    ENDDO

!-------------------------------------------------------------------------------
!- Section 8: Clean up
!-------------------------------------------------------------------------------
!CDIR NOIEXPAND
    call destruct(obs_rad_loc)
    
!CDIR NOIEXPAND
    call destruct(radv)

    DEALLOCATE(ifcaststep, iob_tot, job_tot, iindex, isort_index, inode_index)
   
!-------------------------------------------------------------------------------
!- Section 9: Allocate  & Prepare Fields fg calculations
!-------------------------------------------------------------------------------
  
    r     => obs_rad(isat)
    rs    => r%i
    inobs =  r%n_rec
    call get_n_add(n_add)
    IF (my_cart_id == 0) THEN
      WRITE (nurad,*) '  '
      WRITE (nurad,'(T4,A,T44,I4)') 'Extrapolation type         :', extrp_type
      WRITE (nurad,'(T4,A,T44,I4)') 'Number of additional levels:', n_add
      ! flush YURAD file
      CLOSE (nurad)
      OPEN(nurad, FILE=yurad, FORM='FORMATTED', STATUS='OLD', &
           POSITION='APPEND', IOSTAT=istat)      
    END IF
    r%n_lev = ke + 1 + n_add

    ALLOCATE(r%mdlsfc  (             inobs),              &
             r%plevel  (  rs%n_chan ,inobs),              &
             r%bt_fg   (  rs%n_chan ,inobs),              &
             r%bt_bcor (  rs%n_chan ,inobs),              &
             r%bcor_   (  rs%n_chan ,inobs),              &
             r%p       (  r%n_lev   ,inobs),              &
             r%hh      (  r%n_lev   ,inobs),              &
             r%h_fg    (  rs%n_instr,inobs),              &
             r%t_fg    (  r%n_lev   ,inobs),              &
             r%q_fg    (  r%n_lev   ,inobs),              &
             r%t2m     (             inobs),              &
             r%q2m     (             inobs),              &
             r%ps_fg   (             inobs),              &
             r%ts_fg   (             inobs),              &
             r%u10_fg  (             inobs),              &
             r%v10_fg  (             inobs),              &
             r%cld_fg  (6,r%n_lev-1 ,inobs),              &
             r%cld_frc (6,r%n_lev-1 ,inobs),              &
             stat = istat)
    
    IF (istat==0 .and. btest(rad_output, OUT_CSB)) ALLOCATE(r%bt_fg_cs (rs%n_chan,inobs),stat=istat)
    IF (istat==0 .and. btest(rad_output, OUT_ASR)) ALLOCATE(r%rad_fg   (rs%n_chan,inobs),stat=istat)
    IF (istat==0 .and. btest(rad_output, OUT_CSR)) ALLOCATE(r%rad_fg_cs(rs%n_chan,inobs),stat=istat)

    IF (istat /= 0) THEN
      ierrstat = 99112
      PRINT *, ' ERROR   *** Allocating memory failed. ***'
      RETURN
    ENDIF 

    r%mdlsfc(:)    = 0
    r%bt_fg(:,:)   = rmdi
    r%bt_bcor(:,:) = rmdi
    r%bcor_(:,:)   = 0._wp
    IF (btest(rad_output, OUT_CSB)) r%bt_fg_cs (:,:)   = rmdi
    IF (btest(rad_output, OUT_ASR)) r%rad_fg   (:,:)   = rmdi  
    IF (btest(rad_output, OUT_CSR)) r%rad_fg_cs(:,:)   = rmdi
    r%p(:,:)       = rmdi
    r%hh(:,:)      = rmdi
    r%h_fg(:,:)    = rmdi
    r%t_fg(:,:)    = rmdi
    r%q_fg(:,:)    = rmdi
    r%t2m(:)       = rmdi
    r%q2m(:)       = rmdi
    r%ps_fg(:)     = rmdi
    r%ts_fg(:)     = rmdi
    r%u10_fg(:)    = rmdi
    r%v10_fg(:)    = rmdi
    r%plevel(:,:)  = rmdi
    r%sinfl(:,:)   = rmdi

    IF (l_H) THEN
      ALLOCATE(r%H_t (rs%n_chan, r%n_lev, inobs),         &
               r%H_q (rs%n_chan, r%n_lev, inobs),         &
               r%H_ts(rs%n_chan,          inobs),         &
               r%H_ps(rs%n_chan,          inobs),         &
               stat = istat)
      IF (istat /= 0) THEN
        ierrstat = 99113
        PRINT *, ' ERROR   *** Allocating memory failed. ***'
        RETURN
      ENDIF
    END IF

    IF (.NOT. ASSOCIATED(r%sunzen)) THEN
      ALLOCATE(r%sunzen(inobs),stat=istat)

      IF (istat /= 0) THEN
        ierrstat = 99114
        PRINT *, ' ERROR   *** Allocating memory failed. ***'
        RETURN
      ENDIF 

      r%sunzen(:) = rmdi
    END IF  
  ENDDO
  
!-------------------------------------------------------------------------------
!- Section 10: Compute RTTOV instrument information
!-------------------------------------------------------------------------------
  ! This may reserve more space than necessary if same instruments
  ! are input. Anyway, I dont see much benefit from saving some bytes
  ! to double looping over all the stuff again
    ALLOCATE (sensors(SUM(obs_rad(:)%i%n_instr)), stat=istat)

  IF (istat /= 0) THEN
    ierrstat = 99115
    PRINT *, ' ERROR   *** Allocating memory failed. ***'
    RETURN
  ENDIF 
  
  !-------------------------------------------------------------------------------
  !+ Section 10.1: Merge channel information for all input files 
  !-------------------------------------------------------------------------------
  
  DO isat=1,nsat
    
    r  => obs_rad(isat)
    rs => r%i
   
    DO iinstr = 1, rs % n_instr
      sensors(nsensors+1)%id = &
                   (/ rs%platform   , &
                      rs%rttov_satid, &
                      rs%instr(iinstr)/)
                   
      ! Search already defined sensors for sensor
      DO icnt=1,nsensors
        IF (ALL(sensors(icnt)%id(:) == sensors(nsensors+1)%id(:))) THEN
          rs%rttov_indx(iinstr) = icnt
          EXIT
        ENDIF
      ENDDO

      ! sensor not yet in use by any other radset, store it
      IF (rs%rttov_indx(iinstr) == -1) THEN
        nsensors = nsensors + 1
        rs%rttov_indx(iinstr) = nsensors
      ENDIF
      
      sensor => sensors(rs%rttov_indx(iinstr))

      ! check if the sensor entry misses some of the channels in this rad set
      IF (ASSOCIATED(sensor%channels)) THEN
        imissing = 0
        DO icnt=rs%o_ch_i(iinstr) + 1,rs%o_ch_i(iinstr) + rs%n_ch_i(iinstr)
          IF (.NOT. ANY(sensor%channels == rs%chan(icnt))) THEN
            imissing = imissing+1
          ENDIF
        ENDDO
      ELSE
        imissing = rs%n_ch_i(iinstr)
      ENDIF  

      IF (imissing > 0) THEN
        istart = MINVAL(rs%chan(rs%o_ch_i(iinstr) + 1:&
                                rs%o_ch_i(iinstr) + rs%n_ch_i(iinstr)))
        iend   = MAXVAL(rs%chan(rs%o_ch_i(iinstr) + 1:&
                                rs%o_ch_i(iinstr) + rs%n_ch_i(iinstr)))

        IF (ASSOCIATED(sensor%channels)) THEN
          ALLOCATE(idummy(imissing + SIZE(sensor%channels)))
          istart = MIN(istart,MINVAL(sensor%channels))
          iend   = MAX(iend,  MAXVAL(sensor%channels))
        ELSE
          ALLOCATE(idummy(imissing))
        ENDIF  

        ! merge existing sensor channels with radset's channels 
        ! sensors are automatically sorted  
        icnt = 0
        DO ichan=istart,iend
          IF (ASSOCIATED(sensor%channels)) THEN
            IF (ANY(ichan == sensor%channels)) THEN
              icnt = icnt+1
              idummy(icnt)= ichan
              CYCLE
            ENDIF
          ENDIF  
            
          IF (ANY(ichan == rs%chan(rs%o_ch_i(iinstr) + 1:                      &
                                   rs%o_ch_i(iinstr) + rs%n_ch_i(iinstr)))) THEN
              icnt = icnt+1
              idummy(icnt)= ichan
              CYCLE
          ENDIF    
        ENDDO

        ! interchange arrays
        IF (ASSOCIATED(sensor%channels)) DEALLOCATE(sensor%channels)
        sensor%channels => idummy
      ENDIF
    ENDDO
  ENDDO
 
  !-------------------------------------------------------------------------------
  !+ Section 10.2: map input file channels to rttov channel indices
  !-------------------------------------------------------------------------------
  ! The following is needed, as we only load these channels in rttov we
  ! really need. Therefore the channel indices of rttov will change e.g.
  ! input data channels 4-7 -> rttov channels 1-3
  DO isat=1,nsat
    r  => obs_rad(isat)
    rs => r%i
    
    ALLOCATE(rs%ichan(rs%n_chan), stat = istat)

    IF (istat /= 0) THEN
      ierrstat = 99116
      PRINT *, ' ERROR   *** Allocating memory failed. ***'
      RETURN
    ENDIF 
    
    DO iinstr = 1, rs % n_instr
      sensor => sensors(rs%rttov_indx(iinstr))
      DO icnt=rs%o_ch_i(iinstr) + 1,rs%o_ch_i(iinstr) + rs%n_ch_i(iinstr)
        rs%ichan(icnt) = &
          COUNT(sensor%channels(:) <= rs%chan(icnt))
      ENDDO
    ENDDO  

  ENDDO

  !> \todo allow user control of clouds in RTTOV
  WHERE (sensors(:)%id(3) == 21) ! SEVIRI with clouds
    sensors(:)%addcloud = .True.
  ELSEWHERE
    sensors(:)%addcloud = .False.
  END WHERE  
 
!-------------------------------------------------------------------------------
!- Section 11: Prepare superobbing
!-------------------------------------------------------------------------------
  IF (hrit_superob > 1 .and. nsat > nsatpp) THEN
    isat  =  nsat
    r     => obs_rad(isat)
    rs    => r%i

    IF (my_cart_id == 0) THEN
      inobs = sum(nobspe_tot(:,isat))
    ELSE
      inobs = 1
    END IF

    ALLOCATE(iline(inobs), icol(inobs), istep(inobs), ibox_so(inobs), stat=istat)
    IF (istat /= 0) THEN
      ierrstat = 99136
      PRINT *, ' ERROR   *** Allocating memory failed. ***'
      RETURN
    ENDIF

    CALL gatherv_values(r%ntstep(:),nobspe_tot(my_cart_id+1,isat),nobspe_tot(my_cart_id+1,isat),&
         num_compute,istep,izdum,inobs,0,imp_integers,icomm_cart,ymsg,ierrstat,init=.false.)
    CALL gatherv_values(r%scanl (:),nobspe_tot(my_cart_id+1,isat),nobspe_tot(my_cart_id+1,isat),&
         num_compute,iline,izdum,inobs,0,imp_integers,icomm_cart,ymsg,ierrstat,init=.false.)
    CALL gatherv_values(r%fov   (:),nobspe_tot(my_cart_id+1,isat),nobspe_tot(my_cart_id+1,isat),&
         num_compute,icol, izdum,inobs,0,imp_integers,icomm_cart,ymsg,ierrstat,init=.false.)
    IF (my_cart_id == 0) THEN

      nb_t  = 1
      nb_c  = hrit_superob
      nb_l  = hrit_superob
      nso   = nb_t * nb_c * nb_l
      nb2_t = nb_t / 2
      id_t  = 1
      mn_t  = minval(istep(:))
      mx_t  = maxval(istep(:))
      nb2_c = hrit_superob / 2
      mn_c  = minval(icol (:))
      mx_c  = maxval(icol (:))
      id_c  = thin_col
      nb2_l = hrit_superob / 2
      mn_l  = minval(iline(:))
      mx_l  = maxval(iline(:))
      id_l  = thin_line

      WRITE (nurad,*) '  '
      WRITE (nurad,*) '  Superobbing:'
      WRITE (nurad,*) '    time   box size: ',nb_t
      WRITE (nurad,*) '    line   box size: ',nb_l
      WRITE (nurad,*) '    column box size: ',nb_c
      WRITE (nurad,*) '    box size       : ',nso

      ALLOCATE(iso(nso), stat=istat)
      IF (istat /= 0) THEN
        ierrstat = 99137
        PRINT *, ' ERROR   *** Allocating memory failed. ***'
        RETURN
      ENDIF
      ibox_so(:)  = 0

      it    = mn_t + nb2_t * id_t
      il    = mn_l + nb2_l * id_l
      ic    = mn_c + nb2_c * id_c
      ibox  = 0

      DO WHILE (it <= mx_t - nb2_t*id_t .and. &
                il <= mx_l - nb2_l*id_l .and. &
                ic <= mx_c - nb2_c*id_c)
        mso = 0
        DO i = 1, inobs
          IF (ibox_so(i) == 0) THEN
            IF (istep(i) >= it - nb2_t*id_t .and. istep(i) <= it + nb2_t*id_t) THEN
              IF (iline(i) >= il - nb2_l*id_l .and. iline(i) <= il + nb2_l*id_l) THEN
                IF (icol (i) >= ic - nb2_c*id_c .and. icol (i) <= ic + nb2_c*id_c) THEN
                  mso = mso + 1
                  iso(mso) = i
                END IF
              END IF
            END IF
          END IF
        END DO
        IF (mso > nso) THEN
          ierrstat = 99137
          PRINT *, ' ERROR   *** Superobbing failed. *** '
          RETURN
        ELSEIF (mso == nso) THEN
          ibox = ibox + 1
          ibox_so(iso(1:nso)) = ibox
          DO i = 1, nso
            IF (istep(iso(i))==it .and. iline(iso(i))==il .and. icol(iso(i))==ic) ibox_so(iso(i)) = -ibox
          END DO
          ic = ic + nb_c*id_c
        ELSE
          ic = ic + id_c
        END IF
        IF (ic > mx_c - nb2_c*id_c) THEN
          il = il + id_l
          ic = mn_c + nb2_c*id_c
          IF (il > mx_l - nb2_l*id_l) THEN
            it = it + id_t
            DO WHILE (it <= mx_t - nb2_t*id_t .and. &
                      .not.any(istep(:) == it))
              it = it + id_t
            END DO
            il = mn_l + nb2_l*id_l
          END IF
        END IF
      END DO
      nbox = ibox

      DEALLOCATE(iso)
      WRITE (nurad,*) '    number of boxes: ',nbox
      CLOSE (nurad)
      OPEN(nurad, FILE=yurad, FORM='FORMATTED', STATUS='OLD', &
           POSITION='APPEND', IOSTAT=istat)

    END IF

    ! Misuse i_reprt as ibox_so, i.e. number of superobbing box
    ALLOCATE(r%i_reprt(r%n_rec), stat=istat)
    IF (istat /= 0) THEN
      ierrstat = 99139
      PRINT *, ' ERROR   *** Allocating memory failed. ***'
      RETURN
    ENDIF
    call scatterv_values(ibox_so, inobs, 0, nobspe_tot(:,isat), r%i_reprt, &
         r%n_rec, my_cart_id, num_compute, imp_integers, icomm_cart, ymsg, istat)

    IF (istat /= 0) THEN
      ierrstat = 99140
      PRINT *, ' ERROR   *** scatterv_values ibox_so failed. *** ', istat
      RETURN
    ENDIF

    DEALLOCATE(iline,icol,istep,ibox_so)

  END IF

!-------------------------------------------------------------------------------
!- Section 12: Prepare some globally used variables
!-------------------------------------------------------------------------------

  READ (ydate_ini,'(I4,I2,I2,I2)') reftime%year, reftime%month, reftime%day, &
                                   reftime%hour
  
  reftime%minute = 0
!-------------------------------------------------------------------------------
!- Section 13: Prepare fdbk file output (pe 0 only
!-------------------------------------------------------------------------------
  IF (.not.any(sensors(rs%rttov_indx(1:rs%n_instr))%addcloud)) THEN
    ! Avoid writing of cloud parameters if not necessary
    if (associated(obs_rad(isat)%cld_fg)) deallocate(obs_rad(isat)%cld_fg)
    if (associated(obs_rad(isat)%cld_frc)) deallocate(obs_rad(isat)%cld_frc)
  END IF
  IF (my_cart_id == 0) THEN
    fdbk_offset(:) = 0
    yversion = yncglob_source 
    IF (yversion(1:5) == 'COSMO') yversion = yversion(7:LEN_TRIM(yversion)-6)

    iveri_ref_date = reftime%year*10000 + reftime%month*100 + reftime%day
    iveri_ref_time = reftime%hour * 100 + reftime%minute 
    
    iveri_start = hstart * 60.0_wp
    iveri_end   = hstop  * 60.0_wp

    resolution  = (/dlat, dlon/)
    domain      = (/ie_tot, je_tot, ke_tot/)
    pole        = (/startlat_tot, startlon_tot/)
    lower_left  = pole
    upper_right = pole + resolution * (/je_tot-1,ie_tot-1/)

    WRITE (yveri_ref_datim,'(I8.8,I4.4)') iveri_ref_date, iveri_ref_time

    IF (nvers >= 0 .AND. nvers < 1048576) THEN
      iveri_exp_id    = MOD(nvers,16384)
      iveri_run_class = nvers/16384
    ELSE
      iveri_exp_id    = imdi
      iveri_run_class = imdi
    ENDIF

    IF (leps) THEN
      WRITE (yfdbk_descript,'(A,1X,I3)') 'ensemble forecast member', &
                                        iepsmem
      iveri_epsmem = iepsmem                                  
    ELSE
      yfdbk_descript = 'deterministic forecast'
      iveri_epsmem   = -1                                  
    ENDIF  

    IF ( iveri_run_class == rc_ass .OR. iveri_end <= 300) THEN
      iveri_run_type = vt_firstguess
    ELSE
      iveri_run_type = vt_analysis
    ENDIF  

    iveri_forecast_time = INT(hstop * 100)
    
    WRITE (nurad,*) '  '
    WRITE (nurad,'(T4,A,T14,A,T24,A,T32,A)') '     #obs' , '#obs fdbk','Superob', 'file'
    DO isat=1,nsat
      r     => obs_rad(isat)
      rs    => r%i
      IF (hrit_superob > 1 .and. isat > nsatpp) THEN
        inobs = nbox
      ELSE
        inobs = SUM(nobspe_tot(:,isat))
      END IF

      IF (inobs <= 0) CYCLE

      CALL setup_fdbk(fdbk(isat)%nc, latex = .FALSE.)

      DO iinstr=1,rs%n_instr
        WRITE (ymsg((iinstr-1)*3 + 1:),'(I2.2,A1)') rs%instr(iinstr), '-'
      ENDDO
      
      WRITE (path,'(A,A,I3.3,A,I2.2,A,A14,A)')          &
           TRIM(root%ydir),'/fof_RAD_', rs%satid,'_instr', rs%grid, &
           '_', ydate_ini(1:14) ,'.nc'
      WRITE (nurad,'(T4,I9,T14,I9,T24,L7,T32,A)') SUM(nobspe_tot(:,isat)),inobs, &
           (hrit_superob>1 .and. isat>nsatpp),trim(path)

      CALL create_fdbk(fdbk(isat),                          &
                       TRIM(path),                          &
                       'COSMO',                             &
                       yversion, yncglob_institution,       &
                       inobs, inobs * rs%n_chan,            &
                       iveri_ref_date,                      &
                       iveri_ref_time,                      &
                       iveri_start,                         &
                       iveri_end,                           &
                       resolution,                          &
                       domain,                              &
                       yfdbk_descript,                      &
                       yveri_ref_datim,                     &
                       pole        = pole,                  &
                       lower_left  = lower_left,            &
                       upper_right = upper_right,           &
                       opt         ='RAD'          )

      CALL add_verification(fdbk(isat),                     &
                            'COSMO',                        &
                            iveri_run_type,                 &
                            iveri_run_class,                &
                            yveri_ref_datim,                &
                            iveri_forecast_time,            &
                            resolution,                     &
                            domain,                         &
                            'first guess',                  &
                            iveri_epsmem,                   &
                            iveri_exp_id,                   &
                            varid                    )

      IF (btest(rad_output, OUT_CSB)) THEN
        CALL add_verification(fdbk(isat),                   &
                              'COSMO',                      &
                              iveri_run_type,               &
                              iveri_run_class,              &
                              yveri_ref_datim,              &
                              iveri_forecast_time,          &
                              resolution,                   &
                              domain,                       &
                              'first guess',                &
                              iveri_epsmem,                 &
                              iveri_exp_id,                 &
                              varid,                        &
                              operator_flag=OF_BT_CLEAR_SKY)
      END IF
      IF (btest(rad_output, OUT_ASR)) THEN
        CALL add_verification(fdbk(isat),                   &
                              'COSMO',                      &
                              iveri_run_type,               &
                              iveri_run_class,              &
                              yveri_ref_datim,              &
                              iveri_forecast_time,          &
                              resolution,                   &
                              domain,                       &
                              'first guess',                &
                              iveri_epsmem,                 &
                              iveri_exp_id,                 &
                              varid,                        &
                              operator_flag=OF_RAD_CLOUDY)
      END IF
      IF (btest(rad_output, OUT_CSR)) THEN
        CALL add_verification(fdbk(isat),                   &
                              'COSMO',                      &
                              iveri_run_type,               &
                              iveri_run_class,              &
                              yveri_ref_datim,              &
                              iveri_forecast_time,          &
                              resolution,                   &
                              domain,                       &
                              'first guess',                &
                              iveri_epsmem,                 &
                              iveri_exp_id,                 &
                              varid,                        &
                              operator_flag=OF_RAD_CLEAR_SKY)
      END IF

      WRITE (path,'(A,A,I3.3,A,I2.2,A,A12,A)')          &
             TRIM(root%ydir),'/monRTOVP_', rs%satid,'_instr', rs%grid, &
             '_', ydate_ini(1:12) ,'.nc'
      obs_rad(isat)%filename = trim(path)
      if (any(sensors(rs%rttov_indx(1:rs%n_instr))%addcloud)) then
        CALL write_rtovp((/obs_rad(isat)/), monitor_prof, mode=WRM_OPEN, status=istat, &
             name=trim(path), n_prof=inobs, cloud_mask=(/1, 3, 6/), date=ydate_ini,    &
             lH=l_H)
      else
        CALL write_rtovp((/obs_rad(isat)/), monitor_prof, mode=WRM_OPEN, status=istat, &
             name=trim(path), n_prof=inobs, date=ydate_ini,                            &
             lH=l_H)
      end if
      IF (istat /= 0) THEN
        ierrstat = 99117
        PRINT *, ' ERROR   *** Opening monRTOVP.nc failed. ***'
        PRINT *, ' ERROR   *** ',istat
        RETURN
      ENDIF
      
    ENDDO
  ENDIF

!-------------------------------------------------------------------------------
!- Section 12: Writeout some information
!-------------------------------------------------------------------------------
  IF (my_cart_id == 0) THEN
    WRITE (nurad,*) '  '
    WRITE (nurad,'(T4,A,T20,A,T77,A)') 'Sensor Triplet' , 'Channels', 'Addclouds'
    WRITE (nurad,*) '  '

    DO icnt=1,nsensors
      ymsg = ''

      DO ichan = 1,SIZE(sensors(icnt)%channels)
        IF (LEN_TRIM(ymsg) > 50) THEN
          WRITE (ymsg(LEN_TRIM(ymsg)+1:),*) '...'
          EXIT
        ENDIF
        
        IF (ichan == 1) THEN
          WRITE (ymsg(LEN_TRIM(ymsg)+1:),*) sensors(icnt)%channels(ichan)
        ELSEIF (ichan == SIZE(sensors(icnt)%channels)) THEN
          IF (sensors(icnt)%channels(ichan) - 1 == &
              sensors(icnt)%channels(ichan-1)      ) THEN
            WRITE (ymsg(LEN_TRIM(ymsg)+1:),*) &
                   '-',sensors(icnt)%channels(ichan)
          ENDIF         
        ELSE
          IF (sensors(icnt)%channels(ichan) - 1 /= &
              sensors(icnt)%channels(ichan-1)      ) THEN
            IF (ichan > 2) THEN
              IF (sensors(icnt)%channels(ichan-1) - 1 == &
                   sensors(icnt)%channels(ichan-2)        ) THEN
                WRITE (ymsg(LEN_TRIM(ymsg)+1:),*) &
                     '-',sensors(icnt)%channels(ichan-1)
              ENDIF
            END IF
            WRITE (ymsg(LEN_TRIM(ymsg)+1:),*) &
                  ', ', sensors(icnt)%channels(ichan)
          ENDIF
        ENDIF
      ENDDO

      WRITE (nurad,'(T4,3(1x,I2),T20,A,T77,L1)') sensors(icnt)%id(:), TRIM(ymsg), sensors(icnt)%addcloud
    ENDDO
    ! flush YURAD file
    CLOSE (nurad)
    OPEN(nurad, FILE=yurad, FORM='FORMATTED', STATUS='OLD', &
         POSITION='APPEND', IOSTAT=istat)      
  ENDIF
!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------
#endif
  ierrstat = 0

END SUBROUTINE input_obs_satpp

!==============================================================================
!+ Procedure in "src_obs_sat" for calculating the radiances
!------------------------------------------------------------------------------

SUBROUTINE calc_obs_satpp (ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine calculates the radiances
!
!------------------------------------------------------------------------------

! Parameters
  INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
    ierrstat        ! error status variable

! Local variables
#ifdef RTTOV10
  INTEGER   (KIND=iintegers) ::  &
    istat,                       & ! Status variable
    inobs,                       & ! Number of observables to calc fg for
    inobs_all,                   & ! counting variable for obs
    iradv,                       & ! Counter for radv (satpp input file)
    iobs,                        & ! Counter for obs
    ichan,                       & ! Counter for channels
    dummy,                       & ! surface type
    isurfsens,                   & ! grid instrument index
    icnt,iinstr,                 & ! Counting variables
    ipstart,ipend,               & ! Array slice offsets
    icstart,icend,               & ! Array slice offsets
    iprof,                       & ! profile number counter
    icalc,                       & ! calc. counter
    j,                           & ! loop index
    ncalcs,                      & ! number of calculations for rttov
    nlev,                        & ! number of levels for rttov calcs
    nchan,                       & ! number of channel for instrument
    ierr                           ! error status
  INTEGER   (KIND=iintegers), ALLOCATABLE :: &
    iob_loc(:),                  & ! array of zonal gridpoint indices (local)
    job_loc(:),                  & ! array of meridional gridpoint indices (local)
    profs(:),                    & ! array of profile indizes
    chans(:),                    & ! array of channel indizes
    idg(:),                      &
    ish(:),                      &
    wtype(:),                    &
    stype(:),                    &
    ibox_so_loc(:)                 ! Superobbing box ID
  REAL      (KIND=wp), ALLOCATABLE :: &
    emiss(:,:),                  &
    emiss_k(:,:),                &
    temp_k(:,:,:),               &
    humi_k(:,:,:),               &
    t2m_k (:,:),                 &
    q2m_k (:,:),                 &
    stemp_k (:,:),               &
    psurf_k (:,:)
  REAL      (KIND=sp), ALLOCATABLE :: &
    oe_var(:)
  REAL      , ALLOCATABLE, TARGET :: &
    report_veridata(:,:,:)     
  REAL      (KIND=wp), POINTER :: &
    h_ice_dummy(:,:,:),          & !
    qg_dummy   (:,:,:)
  CHARACTER(LEN=256) ::          &
    ymsg                           ! error message
  TYPE(t_radv), POINTER ::       &
    r  
  TYPE(t_rad_set), POINTER ::    &
    rs  
  TYPE(t_sensor), POINTER ::     &
    sensor  
  TYPE(t_radv), TARGET ::        &
    obs
  TYPE(t_account), ALLOCATABLE:: &
    reports(:)
  TYPE(t_acc_body), ALLOCATABLE, TARGET :: &
    report_body(:,:)
#endif    

! Tracer pointers
REAL (KIND=wp), POINTER ::   &
  qv  (:,:,:) => NULL(),         & ! QV at nnow
  qc  (:,:,:) => NULL(),         & ! QC at nnow
  qg  (:,:,:) => NULL(),         & ! QG at nnow
  qi  (:,:,:) => NULL(),         & ! QI at nnow
  qs  (:,:,:) => NULL()            ! QS at nnow

CHARACTER (LEN=50)  :: yerror       ! error message
CHARACTER (LEN=25)  :: yzroutine

!for reading NWCSAF SEVIRI cloud type
CHARACTER (LEN=14)         :: yzactdate
CHARACTER (LEN=28)         :: yzdum  ! actual date in the form   wd   dd.mm.yy  hh UTC

REAL (KIND= wp)        ::    &
     rzdum,                      & ! 
     ctype(ie,je),               & ! cloud type
     ctype_tot(ie_tot,je_tot),   &
     cth(ie,je),                 &
     cth_tot(ie_tot,je_tot)

INTEGER (KIND=iintegers)   ::    &
     izdum,                      & ! dummy variable
     npr_block,                  & ! Number of profiles in block
     jpstart,                    & ! First profile in block
     jpend,                      & ! Last profile in block
     kpstart,                    & ! First profile in block
     kpend,                      & ! Last profile in block
     jnobs,                      & ! Number of profiles in current block
     nveri                         ! Number of d_veri entries in feedback files

!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE calc_obs_satpp
!-------------------------------------------------------------------------------

  ierrstat = 0
  yerror   = '   '
  yzroutine= 'calc_obs_satpp'

  ! retrieve the required microphysics tracers
  CALL trcr_get(ierrstat, idt_qv, ptr_tlev = nnow, ptr=qv)
  IF (ierrstat /= 0) THEN
    yerror = trcr_errorstr(ierrstat)
    CALL model_abort(my_cart_id, ierrstat, yerror, yzroutine)
  ENDIF
  CALL trcr_get(ierrstat, idt_qc, ptr_tlev=nnow, ptr=qc)
  IF (ierrstat /= 0) THEN
    yerror = trcr_errorstr(ierrstat)
    CALL model_abort(my_cart_id, ierrstat, yerror, yzroutine)
  ENDIF
  CALL trcr_get(ierrstat, idt_qg, ptr_tlev=nnow, ptr=qg)
  IF (ierrstat /= 0 .AND. ierrstat /= T_ERR_NOTFOUND) THEN
    yerror = trcr_errorstr(ierrstat)
    CALL model_abort(my_cart_id, ierrstat, yerror, yzroutine)
  ENDIF
  CALL trcr_get(ierrstat, idt_qi, ptr_tlev=nnow, ptr=qi)
  IF (ierrstat /= 0 .AND. ierrstat /= T_ERR_NOTFOUND) THEN
    yerror = trcr_errorstr(ierrstat)
    CALL model_abort(my_cart_id, ierrstat, yerror, yzroutine)
  ENDIF
  CALL trcr_get(ierrstat, idt_qs, ptr_tlev=nnow, ptr=qs)
  IF (ierrstat /= 0 .AND. ierrstat /= T_ERR_NOTFOUND) THEN
    yerror = trcr_errorstr(ierrstat)
    CALL model_abort(my_cart_id, ierrstat, yerror, yzroutine)
  ENDIF

#ifdef RTTOV10
  IF (ALLOCATED(h_ice)) THEN
    h_ice_dummy => h_ice(:,:,:)
  ELSE
    ALLOCATE(h_ice_dummy(0,0,nnow))
  ENDIF
  IF (ASSOCIATED(qg)) THEN
    qg_dummy => qg(:,:,:)
  ELSE
    ALLOCATE(qg_dummy(0,0,0))
  ENDIF
 
  radv_loop: DO iradv=1,SIZE(obs_rad)
    
    r    => obs_rad(iradv)
    
    rs   => r%i
    nlev =  r%n_lev
    
!-------------------------------------------------------------------------------
!- Section 1: Prepare obs to be computede
!-------------------------------------------------------------------------------
    ipstart = COUNT((r%ntstep(:) - ntstep) < 0 ) + 1
    ipend   = COUNT((r%ntstep(:) - ntstep) <= 0)
    inobs = ipend - ipstart + 1

    !- Check, whether radiances are to be processed at this time step
    inobs_all = inobs
    CALL global_values(inobs_all,1,'SUM',imp_integers,icomm_cart,-1,ymsg,istat)
    IF (inobs_all <= 0) CYCLE radv_loop


    ! find the index of the grid instrument
    isurfsens = 0
    
    IF (rs%flag_instr >= 0) THEN
      DO icnt=1,rs%n_instr
        IF (rs%flag_instr == rs%instr(icnt)) isurfsens = icnt
      END DO  
    ELSE
      DO icnt=1,rs%n_instr
        IF (rs%grid == rs%instr(icnt)) isurfsens = icnt
      END DO  
    ENDIF  

    IF (isurfsens <= 0) THEN
      ierrstat = 99118
      PRINT *, ' ERROR   *** invalid flag_instr in namelist. ***'
      RETURN
    ENDIF    
      
    IF (inobs > 0) THEN
      IF (ANY(r%ntstep(ipstart:ipend) /= ntstep)) THEN
        ierrstat = 99119
        PRINT *, ' ERROR   *** bad order of profiles. ***'
        RETURN
      ENDIF

      ALLOCATE(iob_loc(inobs),                   &
               job_loc(inobs),                   &
               wtype(inobs),                     &
               stype(inobs),                     &
               idg(inobs),                       &
               ish(inobs),                       &
               stat = istat)

      IF (istat /= 0) THEN
          ierrstat = 99120
          PRINT *, ' ERROR   *** Allocating memory failed. ***'
          RETURN
      ENDIF    

      DO iobs=ipstart,ipend
        icnt = iobs - ipstart + 1
        
        iob_loc(icnt) = r%i_box(iobs) - isubpos(my_cart_id,1) + 1 +nboundlines
        job_loc(icnt) = r%j_box(iobs) - isubpos(my_cart_id,2) + 1 +nboundlines
        
        ! fall back to COSMO values if no satpp value available
        IF (r%sunzen(iobs) < rmdich) THEN
          r%sunzen(iobs) = sun_el(iob_loc(icnt),job_loc(icnt))
        ENDIF  
      ENDDO

      !<AS: Reading of NWCSAF SEVIRI cloud type
      IF (lread_ct) THEN
        CALL get_utc_date ( ntstep, ydate_ini, dt, itype_calendar, yzactdate       &
             , yzdum, izdum,rzdum)
        !     =================
        CALL cloud_read_prep(yzactdate,'compute',ctype,cth,istat)
        IF (istat /= 0) THEN
          PRINT *, ' ERROR *** while reading cloud type GRIB file *** '
          ierrstat = 9100
          RETURN
        ENDIF
      ELSE
        ctype = -999._wp
        cth   = -999._wp
      ENDIF !lread_ct
      !AS>


!-------------------------------------------------------------------------------
!- Section 2: Prepare RTTOV profile data from model statee
!-------------------------------------------------------------------------------


      CALL prepare_rttov_input(iob_loc, job_loc, nlev, nlev-ke-1,          &
                               pp(:,:,:,nnow),ps(:,:,nnow),                &
                               t(:,:,:,nnow),t_g(:,:,nnow),                &
                               qv(:,:,:),qc(:,:,:),                        &
                               qi(:,:,:),qs(:,:,:),                        &
                               qg_dummy(:,:,:),h_ice_dummy(:,:,nnow),      &
                               r%p   (:,ipstart:ipend),                    &
                               r%t_fg(:,ipstart:ipend),                    &
                               r%q_fg(:,ipstart:ipend),                    &
                               r%t2m(ipstart:ipend),                       &
                               r%q2m(ipstart:ipend),                       &
                               r%ps_fg(ipstart:ipend),                     &
                               r%h_fg(1,ipstart:ipend),                    &
                               r%u10_fg(ipstart:ipend),                    &
                               r%v10_fg(ipstart:ipend),                    &
                               r%ts_fg(ipstart:ipend),                     &
                               stype(:),                                   &
                               wtype(:),                                   &
                               cloud = r%cld_fg(:,:,ipstart:ipend),        &
                               cfrac = r%cld_frc(:,:,ipstart:ipend),       &
                               hh    = r%hh  (:,ipstart:ipend),            &
                               idg   = idg,                                &
                               ish   = ish,                                &
                               lcloud= any(sensors(rs%rttov_indx(          &
                                               1:rs%n_instr))%addcloud),   &
                               ierrstat = istat) 

      IF (istat /= 0) THEN
          ierrstat = 99121
          PRINT *, ' ERROR   *** prepare rttov_input failed. ***'
          RETURN
      ENDIF    
 
      ! short infor on formats:
      ! surftype satpp: sea,mixed,land (0,1,2)
      ! surftype rttov: land,sea,sea-ice (0,1,2)
      ! we use the satpp encoding for radv instances here!!!
      
      IF (usesurftype == 1 .and. iradv <= nsatpp) THEN
        DO iinstr=1,rs%n_instr
          IF (rs%flag_instr == -2) THEN
            r%h_fg(iinstr,ipstart:ipend) = MAXVAL(r%shgt(:,ipstart:ipend),1)
          ELSE
            r%h_fg(iinstr,ipstart:ipend) = r%shgt(isurfsens,ipstart:ipend)
          ENDIF
        END DO                   

        DO iobs=ipstart,ipend
          IF (rs%flag_instr == -2) THEN
            dummy = MAXVAL(r%stype(:,iobs))
          ELSE
            dummy = r%stype(isurfsens,iobs)
          ENDIF 

          IF (dummy == 2) THEN
           ! if satpp says land(2) force land only
            r%mdlsfc(iobs) = IBSET(r%mdlsfc(iobs),MS_LAND)
            stype(iobs) = 0
          ELSE 
            IF (dummy == 0 .OR. dummy == 1) THEN
              r%mdlsfc(iobs) = IBSET(r%mdlsfc(iobs),MS_SEA)
            ENDIF
            
            IF (dummy == 1) r%mdlsfc(iobs) = IBSET(r%mdlsfc(iobs),MS_LAND)

            IF (stype(iobs) == 2) THEN
              r%mdlsfc(iobs) = IBSET(r%mdlsfc(iobs),MS_ICE)
            ELSE  
              r%mdlsfc(iobs) = IBSET(r%mdlsfc(iobs),MS_NO_ICE)
            END IF
            ! satpp says land/mixed, what now?
          END IF  
          
          IF (rs%flag_instr == -3) THEN
            IF (ANY(r%stype(:,iobs) /= r%stype(1,iobs))) THEN
              DO ichan = 1, rs%n_chan
                CALL change_use_rpt_rad(r,ichan,iobs,FL_SURF)
              ENDDO  
            ENDIF
          ENDIF
        END DO  
      ELSE
        ! use surface from cosmo
        ! prepare_rttov_input generates kilometer heights as
        ! required by rttov, but h_fg in radv is in meters
        r%h_fg(1,ipstart:ipend) = r%h_fg(1,ipstart:ipend) * 1000.0

!CDIR NOLOOPCHG      
        DO iinstr=2,rs%n_instr
          r%h_fg(iinstr,ipstart:ipend) = r%h_fg(1,ipstart:ipend)
        END DO  

        WHERE (stype(:) == 0)
          r%mdlsfc(ipstart:ipend) = IBSET(r%mdlsfc(ipstart:ipend),MS_LAND)
        ELSEWHERE (stype(:) == 1)
          r%mdlsfc(ipstart:ipend) = IBSET(r%mdlsfc(ipstart:ipend),MS_SEA)
          r%mdlsfc(ipstart:ipend) = IBSET(r%mdlsfc(ipstart:ipend),MS_NO_ICE)
        ELSEWHERE (stype(:) == 2)
          r%mdlsfc(ipstart:ipend) = IBSET(r%mdlsfc(ipstart:ipend),MS_SEA)
          r%mdlsfc(ipstart:ipend) = IBSET(r%mdlsfc(ipstart:ipend),MS_ICE)
        END WHERE  

        DO iinstr=2,rs%n_instr
          WHERE (stype(1:inobs) == 1)
            r%stype(iinstr,1:inobs) = 0
          ELSEWHERE (stype(1:inobs) == 2)
            r%stype(iinstr,1:inobs) = 0
          END WHERE
        END DO
        DO iobs = ipstart, ipend
          r%landfr(:,iobs) = fr_land(r%i_box(iobs)-isubpos(my_cart_id,1)+1+nboundlines,&
                                     r%j_box(iobs)-isubpos(my_cart_id,2)+1+nboundlines)
          r%shgt  (:,iobs) = hsurf  (r%i_box(iobs)-isubpos(my_cart_id,1)+1+nboundlines,&
                                     r%j_box(iobs)-isubpos(my_cart_id,2)+1+nboundlines)
        END DO
      END IF

      IF (nprof_rttov > 0) THEN
        npr_block = nprof_rttov
      ELSE
        npr_block = inobs
      END IF

      jpstart = ipstart
      jpend   = ipstart - 1
      block_loop: DO WHILE(jpend < ipend)
        jpend = min(jpstart - 1 + npr_block, ipend)
        jnobs = jpend - jpstart + 1

        kpstart = jpstart - ipstart + 1
        kpend   = jpend   - ipstart + 1

        istat = rttov_fill_input(                          &
            r%p(1:nlev,jpstart:jpend),                     &
            r%t_fg(1:nlev,jpstart:jpend),                  &
            r%q_fg(1:nlev,jpstart:jpend),                  &
            r%t2m(jpstart:jpend),                          &
            r%q2m(jpstart:jpend),                          &
            r%ps_fg(jpstart:jpend),                        &
            r%h_fg(1,jpstart:jpend) * 0.001,               &
            r%u10_fg(jpstart:jpend),                       &
            r%v10_fg(jpstart:jpend),                       &
            r%ts_fg(jpstart:jpend),                        &
            stype(kpstart:kpend),                          &
            r%dlat (jpstart:jpend),                        &
            r%stzen(jpstart:jpend),                        &
            r%sunzen(jpstart:jpend),                       &
            satazim       = r%stazi(jpstart:jpend),        &
            cloud         = r%cld_fg(:,:,jpstart:jpend),   &
            cfrac         = r%cld_frc(:,:,jpstart:jpend),  &
            idg           = idg(kpstart:kpend),            &
            ish           = ish(kpstart:kpend),            &
            watertype     = wtype(kpstart:kpend),          &
            addsolar      = .FALSE.,                       &
            addrefrac     = .TRUE.,                        &
            addinterp     = .TRUE.,                        &
            ivect         = 1,                             &
            pe            = my_cart_id)
!           rttov9_compat = .FALSE.)

        IF (istat /= 0) THEN
          ierrstat = 99122
          PRINT *, ' ERROR   *** rttov_fill_input failed. ***'
          PRINT *, '         *** return code:',istat
          RETURN
        ENDIF    
      
!-------------------------------------------------------------------------------
!- Section 3: Compute First guess
!-------------------------------------------------------------------------------
        ALLOCATE(profs(jnobs * MAXVAL(rs%n_ch_i)),              &
                 chans(jnobs * MAXVAL(rs%n_ch_i)),              &
                 emiss(MAXVAL(rs%n_ch_i), jnobs),               &
                 stat = istat)

        IF (istat /= 0) THEN
          ierrstat = 99123
          PRINT *, ' ERROR   *** Allocating memory failed. ***'
          RETURN
        ENDIF    

        IF (calc_k) THEN
          ALLOCATE(emiss_k(MAXVAL(rs%n_ch_i), jnobs),            &
                   temp_k(nlev,MAXVAL(rs%n_ch_i),jpstart:jpend), &
                   humi_k(nlev,MAXVAL(rs%n_ch_i),jpstart:jpend), &
                   t2m_k (MAXVAL(rs%n_ch_i),jnobs),              &
                   q2m_k (MAXVAL(rs%n_ch_i),jnobs),              &
                   stemp_k (MAXVAL(rs%n_ch_i),jnobs),            &
                   psurf_k(     MAXVAL(rs%n_ch_i),jnobs),        &
                   stat = istat)
          IF (istat /= 0) THEN
            ierrstat = 99124
            PRINT *, ' ERROR   *** Allocating memory failed. *** ',nlev,MAXVAL(rs%n_ch_i),jnobs,jpstart,jpend
            RETURN
          ENDIF
        ENDIF

        instr_loop: DO iinstr = 1, rs%n_instr
          sensor => sensors(rs%rttov_indx(iinstr))
          nchan = rs%n_ch_i(iinstr)
           
          DO iobs = 1, jnobs
            icstart = 1 + (iobs - 1) * nchan
            icend   =     (iobs)     * nchan
            profs(icstart:icend) = iobs
            chans(icstart:icend) = rs%ichan(&
                                       1+rs%o_ch_i(iinstr):nchan+rs%o_ch_i(iinstr))
          ENDDO

          ncalcs = jnobs * nchan

          icstart = 1     + rs%o_ch_i(iinstr)
          icend   = nchan + rs%o_ch_i(iinstr)
            
          emiss(1:nchan,1:jnobs) = 0.0_wp

          istat = rttov_set_input(                            &
                 hsurf=r%h_fg(iinstr,jpstart:jpend)*0.001)

          IF (istat /= 0) THEN
            ierrstat = 99125
            PRINT *, ' ERROR   *** rttov_set_input failed. ***'
            RETURN
          ENDIF

          IF (calc_k) THEN
            emiss_k(1:nchan,1:jnobs) = 0.0_wp

            ! We use the k matrix rttov interface, as we have to calculate
            ! the most sensitive height as well (for LETKF)
            istat = rttov_k_ifc(                                             &
                    rs%rttov_indx(iinstr) + num_sensors,                     &
                    profs(1:ncalcs),                                         &
                    chans(1:ncalcs),                                         &
                    emiss(1:nchan,1:jnobs),                                  &
                    emiss_k(1:nchan,1:jnobs),                                &
                    temp_k (:,1:nchan,jpstart:jpend),                        &
                    humi_k (:,1:nchan,jpstart:jpend),                        &
                    t2m_k  (  1:nchan,:),                                    &
                    q2m_k  (  1:nchan,:),                                    &
                    stemp_k(  1:nchan,:),                                    &
                    t_b         = r%bt_fg(icstart:icend, jpstart:jpend),     &
                    t_b_clear   = r%bt_fg_cs (icstart:icend, jpstart:jpend), &
                    rad         = r%rad_fg   (icstart:icend, jpstart:jpend), &
                    radclear    = r%rad_fg_cs(icstart:icend, jpstart:jpend), &
                    psurf_k     = psurf_k(1:nchan,:),                        &
                    dealloc     = (jpend >= ipend),                          &
                    rad_out_flg = rad_output,                                &
                    pe          = my_cart_id,                                &
                    iprint      = -1)

            IF (istat == 0) THEN
              DO iobs=jpstart,jpend
                SELECT CASE(fdbk_height_method)
                CASE (1)
                  r%plevel(icstart:icend,iobs) = 100.0_wp * &
                       r%p(MAXLOC(ABS(temp_k(:,1:nchan,iobs)),1),iobs)
                CASE (2)
                  IF (ANY( SUM(temp_k(:,1:nchan,iobs) + &
                       humi_k(:,1:nchan,iobs) ,1) == 0)) THEN
                    PRINT *, ' WARNING   *** Normalization not possible. ***'
                  ELSE
                    r%plevel(icstart:icend,iobs) = 100.0_wp * &
                         SUM((ABS(temp_k(:,1:nchan,iobs) * fdbk_e_bg_t) +      &
                              ABS(humi_k(:,1:nchan,iobs) * fdbk_e_bg_rh)) *    &
                             SPREAD(r%p(:,iobs),2,nchan),1) /                 &
                         SUM((ABS(temp_k(:,1:nchan,iobs) * fdbk_e_bg_t) +      &
                              ABS(humi_k(:,1:nchan,iobs) * fdbk_e_bg_rh)), 1)
                  ENDIF
                CASE (3)
                  IF (ANY( SUM(temp_k(:,1:nchan,iobs) + &
                       humi_k(:,1:nchan,iobs) ,1) == 0)) THEN
                    PRINT *, ' WARNING   *** Normalization not possible. ***'
                  ELSE
                    r%plevel(icstart:icend,iobs) = 100.0_wp * &
                         SUM((SQUARE(temp_k(:,1:nchan,iobs) * fdbk_e_bg_t) +   &
                              SQUARE(humi_k(:,1:nchan,iobs) * fdbk_e_bg_rh)) * &
                             SPREAD(r%p(:,iobs),2,nchan),1) /                 &
                         SUM((SQUARE(temp_k(:,1:nchan,iobs) * fdbk_e_bg_t) +   &
                              SQUARE(humi_k(:,1:nchan,iobs) * fdbk_e_bg_rh)), 1)
                  ENDIF
                END SELECT

                ! If the highest level has the highest sensitivity
                ! we peak most likely outside the cosmo levels
                ! so set the plevel to a very small pressure
                WHERE (MAXLOC(ABS(temp_k(:,1:nchan,iobs)),1) == 1)
                  r%plevel(icstart:icend,iobs) = 0.000001
                END WHERE

                ! Calculate influence by surface
                WHERE (emiss(1:nchan,iobs-jpstart+1) > 0._wp .and. emiss(1:nchan,iobs-jpstart+1) <=1._wp)
                  r%sinfl(icstart:icend,iobs) = stemp_k(1:nchan,iobs-jpstart+1) / emiss(1:nchan,iobs-jpstart+1)
                ELSEWHERE
                  r%sinfl(icstart:icend,iobs) = rmdi
                END WHERE

                IF (l_H) THEN
                  r%H_t (icstart:icend,:,iobs) = transpose(temp_k (:,1:nchan,iobs))
                  r%H_q (icstart:icend,:,iobs) = transpose(humi_k (:,1:nchan,iobs))
                  r%H_ts(icstart:icend,  iobs) =           stemp_k(  1:nchan,iobs-jpstart+1)
                  r%H_ps(icstart:icend,  iobs) =           psurf_k(  1:nchan,iobs-jpstart+1)
                END IF

              ENDDO

            ELSE

              WRITE (stderr,*) ' WARNING   *** rttov_k_ifc failed: ',istat,' ***'
              DO iobs = jpstart,jpend
                DO ichan = icstart,icend
                  CALL change_use_rpt_rad(r,ichan,iobs,FL_OPERATOR)
                  r%bt_fg(ichan, iobs) = rmdi
                END DO
              END DO  
            ENDIF  

          ELSE ! calc_k

            istat = rttov_direct_ifc(                                        &
                    rs%rttov_indx(iinstr) + num_sensors,                     &
                    profs(1:ncalcs),                                         &
                    chans(1:ncalcs),                                         &
                    emiss(1:nchan,1:jnobs),                                  &
                    t_b         = r%bt_fg    (icstart:icend, jpstart:jpend), &
                    t_b_clear   = r%bt_fg_cs (icstart:icend, jpstart:jpend), &
                    rad         = r%rad_fg   (icstart:icend, jpstart:jpend), &
                    radclear    = r%rad_fg_cs(icstart:icend, jpstart:jpend), &
                    dealloc     = (jpend >= ipend),                          &
                    rad_out_flg = rad_output,                                &
                    pe          = my_cart_id,                                &
                    iprint      =-1)

            IF (istat /= 0) THEN
              WRITE (stderr,*) ' WARNING   *** rttov_direct_ifc failed: ',istat,' ***'
              DO iobs = jpstart,jpend
                DO ichan = icstart,icend
                  CALL change_use_rpt_rad(r,ichan,iobs,FL_OPERATOR)
                  r%bt_fg(ichan, iobs) = rmdi
                END DO
              END DO
            ENDIF !istat

          END IF ! (calc_k)

        ENDDO instr_loop

        DEALLOCATE(profs, chans, emiss)
        IF (calc_k) DEALLOCATE(emiss_k, temp_k, humi_k, t2m_k,q2m_k, stemp_k, psurf_k)

        jpstart = jpstart + npr_block
      END DO block_loop

      DEALLOCATE(wtype,stype,idg,ish,iob_loc,job_loc)

!-------------------------------------------------------------------------------
!- Section 4.1: Calculate Bias correction and bias corrected obs
!-------------------------------------------------------------------------------

      ! TODO: check if the following is correct
      ! WHERE ( r%state(:,ipstart:ipend) < ST_DISMISS)
      r%bt_bcor(:,ipstart:ipend) =  r%bt_obs(:,ipstart:ipend) - &
           r%bcor_(:,ipstart:ipend)
      ! END WHERE    

    ENDIF  ! inobs
    
!-------------------------------------------------------------------------------
!- Section 4.2: Collect Results at PE 0
!-------------------------------------------------------------------------------

    CALL global_values(inobs,1,'SUM',imp_integers,icomm_cart,0,ymsg,istat)

    IF (my_world_id /= 0) inobs = 0

    CALL construct(obs,inobs,istat,r)
    IF (istat == 0) ALLOCATE(obs%plevel(rs%n_chan,inobs), stat=istat)
   
    IF (istat == 0) THEN
      CALL gather_radv(r, ipstart, ipend, obs, 1, ncalcs, 0,istat)
      ierrstat = istat + ierrstat
      IF (lread_ct) THEN
        CALL gather_field(ctype(:,:),ie,je,ctype_tot(:,:),ie_tot,je_tot,0,istat)
        ierrstat = istat + ierrstat
        CALL gather_field(cth(:,:),ie,je,cth_tot(:,:),ie_tot,je_tot,0,istat)
        ierrstat = istat + ierrstat
        IF ( ierrstat /= 0 ) THEN
          ierrstat = 99126
          PRINT *, ' ERROR   *** failed to gather fields at PE 0 in calc_obs_satpp ***'
        ENDIF
      END IF
    ELSE
       ierrstat = 99127
       PRINT *, ' ERROR   *** failed to allocate memory ***'
    ENDIF

!-------------------------------------------------------------------------------
!- Section 5: Superobbing and append data to feed back file (PE 0 only)
!-------------------------------------------------------------------------------
    r  => obs
    rs => r%i
    IF (.not.any(sensors(rs%rttov_indx(1:rs%n_instr))%addcloud)) THEN
      ! Avoid writing of cloud parameters if not necessary
      IF (associated(r%cld_fg )) deallocate(r%cld_fg )
      IF (associated(r%cld_frc)) deallocate(r%cld_frc)
    END IF

    IF (my_cart_id == 0 .AND. inobs > 0 .AND. istat == 0) THEN

      IF (inobs /= ncalcs ) THEN
        ierrstat = 99128
        PRINT *, ' ERROR   *** obs number mismatch ***'
        RETURN
      ENDIF

  !-----------------------------------------------------------------------------
  !+ Section 5.1: Superobbing
  !-----------------------------------------------------------------------------
      IF (hrit_superob > 1 .and. iradv > nsatpp) THEN
        ALLOCATE(ibox_so_loc(inobs),stat = istat)
        IF (istat /= 0) THEN
          ierrstat = 99141
          PRINT *, ' ERROR   *** allocating memory failed. ***'
          RETURN
        ENDIF
        ibox_so_loc(1:inobs) = r%i_reprt(1:inobs)
        IF (lread_ct) THEN
          ALLOCATE(r%cld_top(inobs), stat=istat)
          IF (istat /= 0) THEN
            ierrstat = 99142
            PRINT *, ' ERROR   *** allocating memory failed. ***'
            RETURN
          ENDIF
          DO iobs = 1, inobs
            r%cld_top(iobs) = cth_tot(r%i_box(iobs),r%j_box(iobs))
          END DO
        END IF
        call superob_radv(r, ibox_so_loc, nso, rmdi, istat)
        IF (istat /= 0) THEN
          ierrstat = 99138
          PRINT *, ' ERROR   *** superobbing failed. ***'
          RETURN
        ENDIF
        inobs = r%n_rec
        DEALLOCATE(ibox_so_loc)
      END IF

  !-----------------------------------------------------------------------------
  !+ Section 5.2: Populate report structure with our values 
  !-----------------------------------------------------------------------------
      
      nveri = 0
      if (btest(rad_output, OUT_ASB)) nveri = nveri + 1
      if (btest(rad_output, OUT_CSB)) nveri = nveri + 1
      if (btest(rad_output, OUT_ASR)) nveri = nveri + 1
      if (btest(rad_output, OUT_CSR)) nveri = nveri + 1
      
      ALLOCATE(reports(inobs),                         &
               report_veridata(nveri,rs%n_chan,inobs), &
               report_body(rs%n_chan,inobs),           &
               stat = istat)

      IF (istat /= 0) THEN
        ierrstat = 99129
        PRINT *, ' ERROR   *** allocating memory failed. ***'
        RETURN
      ENDIF

      DO iobs=1,inobs   !loop over FOVs
        reports(iobs)%offset        = (fdbk_offset(iradv) + iobs-1) * rs%n_chan + 1 
        reports(iobs)%header%i_body = (fdbk_offset(iradv) + iobs-1) * rs%n_chan + 1

        ! Dont allocate each pointer on its own, this
        ! causes memory fragmentation and tons of allocs (slow on sx9)
        ! use a container alloc instead and just point to the memory
        reports(iobs)%body => report_body(:,iobs)


        CALL diff_minutes(reftime%year,reftime%month,reftime%day,         &
                          reftime%hour, reftime%minute,                   &
                          r%date(iobs)     / 10000     ,                  &
                          MOD(r%date(iobs) / 100  ,100),                  &
                          MOD(r%date(iobs)        ,100),                  &
                          r%time(iobs)     / 10000     ,                  &
                          MOD(r%time(iobs) / 100  ,100),                  &
                          itype_calendar, reports(iobs)%header%time,      &
                          istat                                   )

        IF (associated(r%date_d) .and. associated(r%time_d)) THEN
          CALL diff_minutes(reftime%year,reftime%month,reftime%day,          &
                            reftime%hour, reftime%minute,                    &
                            r%date_d(iobs)     / 10000     ,                 &
                            MOD(r%date_d(iobs) / 100  ,100),                 &
                            MOD(r%date_d(iobs)        ,100),                 &
                            r%time_d(iobs)     / 10000     ,                 &
                            MOD(r%time_d(iobs) / 100  ,100),                 &
                            itype_calendar, reports(iobs)%header%time_dbase, &
                            istat                                   )
        ELSE
          reports(iobs)%header%time_dbase = -999
        ENDIF
        
        allocate(oe_var(rs%n_chan))
        call calc_obserr(rs%oe_par(:), oe_var(:), r%stzen(iobs), r%fov(iobs), r%dlat(iobs))
        WHERE (oe_var(:) > 0._wp) 
          report_body(:,iobs)%e_o = SQRT(oe_var(:))
        ELSEWHERE
          report_body(:,iobs)%e_o = -999._wp
        END WHERE
        deallocate(oe_var)
      END DO  

      DO iinstr=1,rs%n_instr
        DO ichan=1+rs%o_ch_i(iinstr),rs%n_ch_i(iinstr)+ rs%o_ch_i(iinstr)
          report_body(ichan,:)%level     = rs%chan(ichan)
          report_body(ichan,:)%level_sig = rs%instr_wmo(iinstr)
        END DO
      END DO  
 
      report_body(:,:)%varno     = VN_RAWBT
      report_body(:,:)%obs       = r%bt_bcor(:,1:inobs)
      report_body(:,:)%bcor      = - r%bcor_  (:,1:inobs)
      report_body(:,:)%level_typ = 0
      report_body(:,:)%qual      = imdi
      report_body(:,:)%plevel    = r%plevel(:,1:inobs) !TODO

      nveri = 0
      IF (btest(rad_output,OUT_ASB)) THEN
        nveri = nveri + 1
        WHERE (r%state(:,1:inobs) < ST_DISMISS)
          report_veridata(nveri,:,:) = r%bt_fg(:,1:inobs)
        ELSEWHERE
          report_veridata(nveri,:,:) = REAL(rmdi)
        END WHERE
      END IF
      IF (btest(rad_output,OUT_CSB)) THEN
        nveri = nveri + 1
        WHERE (r%state(:,1:inobs) < ST_DISMISS)
          report_veridata(nveri,:,:) = r%bt_fg_cs(:,1:inobs)
        ELSEWHERE
          report_veridata(nveri,:,:) = REAL(rmdi)
        END WHERE
      END IF
      IF (btest(rad_output,OUT_ASR)) THEN
        nveri = nveri + 1
        WHERE (r%state(:,1:inobs) < ST_DISMISS)
          report_veridata(nveri,:,:) = r%rad_fg(:,1:inobs)
        ELSEWHERE
          report_veridata(nveri,:,:) = REAL(rmdi)
        END WHERE
      END IF
      IF (btest(rad_output,OUT_CSR)) THEN
        nveri = nveri + 1
        WHERE (r%state(:,1:inobs) < ST_DISMISS)
          report_veridata(nveri,:,:) = r%rad_fg_cs(:,1:inobs)
        ELSEWHERE
          report_veridata(nveri,:,:) = REAL(rmdi)
        END WHERE
      END IF

      !report veri_data needs to be filled AFTER report_veridata has been filled
      DO iobs=1,inobs   !loop over FOVs
        DO icnt = 1, rs%n_chan   ! loop over channels
          reports(iobs)%body(icnt)%veri_data => report_veridata(:,icnt,iobs)
        END DO
      END DO
        
      report_body(1:rs%n_chan,1:inobs)%state = r%state(1:rs%n_chan,1:inobs)
      report_body(1:rs%n_chan,1:inobs)%flags = r%flags(1:rs%n_chan,1:inobs)
      ! Set check variable
      WHERE (btest(report_body(1:rs%n_chan,1:inobs)%flags, FL_AREA))
        report_body(1:rs%n_chan,1:inobs)%check = FL_AREA
      ELSEWHERE (btest(report_body(1:rs%n_chan,1:inobs)%flags, FL_TIME))
        report_body(1:rs%n_chan,1:inobs)%check = FL_TIME
      ELSEWHERE (btest(report_body(1:rs%n_chan,1:inobs)%flags, FL_SURF))
        report_body(1:rs%n_chan,1:inobs)%check = FL_SURF
      ELSEWHERE (btest(report_body(1:rs%n_chan,1:inobs)%flags, FL_OPERATOR))
        report_body(1:rs%n_chan,1:inobs)%check = FL_OPERATOR
      ELSEWHERE
        report_body(1:rs%n_chan,1:inobs)%check = 0
      END WHERE
      ! Report headers
      reports(:)%header%r_state = MINVAL(r%state(:,1:inobs),1)
      reports(:)%header%r_flags = 0
      DO ichan=1,rs%n_chan
        reports(:)%header%r_flags = IOR(r%state(ichan,1:inobs),        &
                                        reports(:)%header%r_flags)
      END DO
      reports(:)%header%r_check       = FL_NONE
      reports(:)%len                  = rs%n_chan
      reports(:)%header%l_body        = rs%n_chan
      reports(:)%header%n_level       = rs%n_chan
      reports(:)%header%data_category = fdbk_rad_cat
      reports(:)%header%sub_category  = fdbk_rad_subcat
      reports(:)%header%obstype       = OT_RAD
      reports(:)%header%codetype      = 0  ! TODO: check this value
      reports(:)%header%ident         = rs%satid
      reports(:)%header%statid        = sat_id2name(rs%satid)
      reports(:)%header%lat           = r%dlat(1:inobs)
      reports(:)%header%lon           = r%dlon(1:inobs)
      reports(:)%header%time_nomi     = reports(:)%header%time
      reports(:)%header%time_dbase    = reports(:)%header%time_dbase
      reports(:)%header%z_station     = -999
      reports(:)%header%z_modsurf     = r%h_fg(isurfsens,1:inobs) ! TODO: instr. dep. height
      reports(:)%header%sta_corr      = 0
      reports(:)%header%index_x       = r%i_box(1:inobs)
      reports(:)%header%index_y       = r%j_box(1:inobs)
      reports(:)%header%dbkz          = -999

      ! TODO: superobbing not implemented for these quantities
      DO iobs=1,inobs
        IF (lread_ct) THEN
          reports(iobs)%header%ct_nwc = ctype_tot(r%i_box(iobs),r%j_box(iobs))
        ELSE
          reports(iobs)%header%ct_nwc = imdi
        END IF
        IF (associated(r%cld_top)) THEN
          reports(iobs)%header%ch_nwc = r%cld_top(iobs)
        ELSE IF (lread_ct) THEN
          reports(iobs)%header%ch_nwc = cth_tot(r%i_box(iobs),r%j_box(iobs))
        ELSE
          reports(iobs)%header%ch_nwc = REAL(rmdich)          
        END IF
      ENDDO

      DO iinstr=1,rs%n_instr
        IF (rs%grid == rs%instr(iinstr)) reports(:)%header%instype = rs%instr_wmo(iinstr)
      END DO

      reports(:)%header%phase         = r%fov(1:inobs)
      reports(:)%header%surftype      = 0

      ! if fl_surf flag ist set, then for all obs of one report
      WHERE (BTEST(r%flags(1,1:inobs),FL_SURF))
        reports(:)%header%surftype = IBSET(0,ST_MISMATCH)
      ELSEWHERE (r%h_fg(isurfsens,1:inobs) > 1000.0_wp)
        reports(:)%header%surftype = IBSET(0,ST_HIGHLAND)
      ELSEWHERE (BTEST(r%mdlsfc(1:inobs),MS_LAND))
        reports(:)%header%surftype = IBSET(0,ST_LAND)
      ELSEWHERE (BTEST(r%mdlsfc(1:inobs),MS_ICE))
        reports(:)%header%surftype = IBSET(0,ST_ICE)
      ELSEWHERE (BTEST(r%mdlsfc(1:inobs),MS_SEA))
        reports(:)%header%surftype = IBSET(0,ST_SEA)
      END WHERE
        
      reports(:)%header%source        = iradv
      reports(:)%header%subset        = r%fov(1:inobs)
      reports(:)%header%sat_zenit     = r%stzen(1:inobs)
      reports(:)%header%mdlsfc        = r%mdlsfc(1:inobs)
      reports(:)%header%sun_zenit     = r%sunzen(1:inobs)

      IF (ASSOCIATED(r%center)) THEN
        reports(:)%header%center = r%center(1:inobs)
      ELSE
        reports(:)%header%center = -1
      ENDIF  

      IF (ASSOCIATED(r%subcenter)) THEN
        reports(:)%header%sub_center = r%subcenter(1:inobs)
      ELSE
        reports(:)%header%sub_center = -1
      ENDIF  

      IF (ASSOCIATED(r%scanl)) THEN
        reports(:)%header%record = r%scanl(1:inobs) 
      ELSE 
        reports(:)%header%record = imdi 
      ENDIF  
      
      IF (ASSOCIATED(r%flg_prc)) THEN
        reports(:)%header%flg_1dvar     = r%flg_prc(1:inobs)
      ELSE  
        reports(:)%header%flg_1dvar     = imdi
      ENDIF

      IF (ASSOCIATED(r%flg_cld)) THEN
        reports(:)%header%flg_cld       = r%flg_cld(1:inobs)
      ELSE  
        reports(:)%header%flg_cld       = imdi
      ENDIF
      
      IF (ASSOCIATED(r%obsnum)) THEN
        reports(:)%header%obs_id        = r%obsnum(1:inobs)
      ELSE
        reports(:)%header%obs_id        = imdi
      ENDIF
     
  !-----------------------------------------------------------------------------
  !+ Section 5.3: Write data to file 
  !-----------------------------------------------------------------------------

      CALL write_report(fdbk(iradv), reports, inobs, inobs*rs%n_chan,    &
                    fdbk_offset(iradv)+1,fdbk_offset(iradv)*rs%n_chan+1, &
                    imdi, REAL(rmdich), istat,ymsg)

      r%p_unit = 'hpa'
      IF (associated(r%cld_top)) DEALLOCATE(r%cld_top) ! Do not write cld_top in monRTOVP, since it is not derived from the 
                                                       ! from the model profile here (it is the NWCSAF cloud top height)
      CALL write_rtovp((/r/), monitor_prof, mode=WRM_WRITE, status=istat, &
           name=trim(r%filename), offset=fdbk_offset(iradv), lH=l_H)
      IF (istat /= 0) THEN
        ierrstat = 99130
        PRINT *, ' ERROR   *** Writing monRTOVP.nc failed. ***'
        RETURN
      ENDIF

      DEALLOCATE(reports, report_veridata, report_body)
      fdbk_offset(iradv) = fdbk_offset(iradv) + inobs
      
      IF (istat /= NF_NOERR) THEN
        ierrstat = 99131
        WRITE (stderr,*) ' ERROR   *** write_report error ',istat,':',TRIM(ymsg),'***'
      ENDIF

  !-----------------------------------------------------------------------------
  !+ Section 5.4: Close file if done 
  !-----------------------------------------------------------------------------

      ! close file if all obs have been written
      !!!! TODO: does not work for superobbing:
      IF (fdbk_offset(iradv) == SUM(nobspe_tot(:,iradv)) .or. &
           hrit_superob > 1 .and. iradv > nsatpp .and. fdbk_offset(iradv) == nbox) THEN
        fdbk(iradv)% nc% error = nf_redef ( fdbk(iradv)% nc% ncid )
        IF (fdbk(iradv)% nc% error /= 0) THEN
          ierrstat = 99132
          write(stderr,*) 'ERROR   *** nf_redef failed ***'
          RETURN
        END IF

        fdbk(iradv)% n_hdr     = fdbk_offset(iradv)
        fdbk(iradv)% nc% error = nf_put_att_int ( fdbk(iradv)%nc% ncid, NF_GLOBAL, &
             'n_hdr',NF_INT, 1, fdbk(iradv)% n_hdr)
        fdbk(iradv)% n_body    = fdbk(iradv)% n_hdr * rs%n_chan
        fdbk(iradv)% nc% error = nf_put_att_int ( fdbk(iradv)%nc% ncid, NF_GLOBAL, &
             'n_body',NF_INT, 1, fdbk(iradv)% n_body )

        fdbk(iradv)% nc% error = nf_enddef ( fdbk(iradv)% nc% ncid )
        IF (fdbk(iradv)% nc% error /= 0) THEN
          ierrstat = 99133
          write(stderr,*) 'ERROR   *** nf_enddef failed ***'
          RETURN
        END IF

        CALL close_fdbk(fdbk(iradv))

        CALL write_rtovp((/r/), monitor_prof, mode=WRM_CLOSE, status=istat, &
             name=trim(r%filename))
        IF (istat /= 0) THEN
          ierrstat = 99134
          PRINT *, ' ERROR   *** Closing '//trim(r%filename)//' failed. ***'
          RETURN
        ENDIF
      ENDIF
    ENDIF

    CALL destruct(obs)

    IF (ierrstat /= 0) RETURN
  ENDDO radv_loop

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------
  IF (.NOT.ALLOCATED(h_ice)) DEALLOCATE(h_ice_dummy)
  IF (.NOT.ASSOCIATED  (qg)) DEALLOCATE(qg_dummy)
#endif

CONTAINS
  ELEMENTAL REAL(KIND=wp) FUNCTION square(val)
  REAL (KIND=wp),INTENT(IN) :: val
    square = val * val
  END FUNCTION square
END SUBROUTINE calc_obs_satpp

#endif

!==============================================================================
!+ Procedure in "src_obs_sat" for preventing supersaturation
!------------------------------------------------------------------------------
ELEMENTAL REAL(KIND=wp) FUNCTION q_sat(p, t)
  REAL(KIND=wp), intent(in) :: p, t
  REAL(KIND=wp) :: p_sat

  p_sat = B1 * exp(B2w * (t - B3) / (t - B4w))
  q_sat = Rdv * p_sat / max((p*100 - O_m_rdv * p_sat), 1.0_wp)

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END FUNCTION q_sat

  !==============================================================================
  !+ Procedure in "src_obs_sat" for calculating additional levels for rttov input
  !------------------------------------------------------------------------------
SUBROUTINE get_n_add(n_add)
  INTEGER(KIND=iintegers), INTENT(OUT) :: n_add
  
  REAL   (KIND=wp)    :: hhlm(2), pprof(2), avg_p(2), fact
  INTEGER(KIND=iintegers) :: i, j, k, iprof, istat
  CHARACTER(LEN=256)      ::  ymsg            ! error message
  
  n_add = 0
  IF (IAND(extrp_type, extrp_const+extrp_lin+extrp_clim) /= 0) THEN
    avg_p(1) = 0
    avg_p(2) = 1
    DO i = 1, ie
      DO j = 1, je
        IF (linterp) THEN
          DO k = 1, 2
            ! height of layers
            hhlm(k) = 0.5_wp*(hhl(i,j,k+1)+hhl(i,j,k))
          END DO
          fact = (hhlm(1)-hhl(i,j,1))/(hhlm(2)-hhlm(1))
          pprof(1) = LOG(p0(i,j,1)) - fact * (LOG(p0(i,j,2)) - LOG(p0(i,j,1)))
          fact = (hhl(i,j,2)-hhlm(2-1))/(hhlm(2)-hhlm(1))
          pprof(2) = LOG(p0(i,j,1)) + fact * (LOG(p0(i,j,2)) - LOG(p0(i,j,1)))
          avg_p(1) = avg_p(1) + pprof(1)
          avg_p(2) = avg_p(2) + pprof(2)
        ELSE
          avg_p(1) = avg_p(1) + LOG(p0(i,j,1))
          avg_p(2) = avg_p(2) + LOG(p0(i,j,2))
        END IF
      END DO
    ENDDO
    
    CALL global_values(avg_p,2,'SUM',imp_reals  ,icomm_cart,0,ymsg,istat)
    IF (istat /= 0) CALL model_abort(my_cart_id,99905, 'gather_values:'//trim(ymsg), &
         'get_n_add', istat)
    IF (my_cart_id == 0) THEN
      IF (IAND(extrp_type, extrp_logp) /= 0) THEN
        avg_p(:) = avg_p(:) / (ie_tot * je_tot)
        IF (avg_p(1) > log(p_top)) n_add = MAX(1,NINT(ABS((log(p_top)-avg_p(1))/(avg_p(1)-avg_p(2)))))
      ELSE
        avg_p(:) = exp(avg_p(:) / (ie_tot * je_tot))
        IF (avg_p(1) > p_top) n_add = MAX(1,NINT(ABS((p_top-avg_p(1))/(avg_p(1)-avg_p(2)))))
      END IF
    END IF
    CALL distribute_values(n_add, 1, 0, imp_integers, icomm_cart, istat)
    IF (istat /= 0) CALL model_abort(my_cart_id,99905, 'distribute_values', &
         'get_n_add', istat)
  ENDIF
  
END SUBROUTINE get_n_add

!==============================================================================
!+ Procedure in "src_obs_sat" for preparing rttov input fields
!------------------------------------------------------------------------------

SUBROUTINE prepare_rttov_input(iob_loc, job_loc, nlev, n_add,           &
                               pp, ps,t, t_g, qv, qc, qi, qs, qg, h_ice,&
                               pres, temp, humidity,                    &
                               t2m, hum2m, psurf, s_hgt, u10m, v10m,t_s,&
                               stype,wtype,                             &
                               cloud, cfrac, hh, idg, ish, lcloud,      &
                               ierrstat)


!------------------------------------------------------------------------------
!
! Description:
!   This subroutine prepares the profiles for rttov
!
!------------------------------------------------------------------------------

! Parameters
  INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
    iob_loc(:),                & ! zonal indices of profile
    job_loc(:)                   ! meridional indices of profile
  REAL   (KIND=wp),   INTENT (IN)   ::        &
    pp(:,:,:),                 & ! pressure state
    ps(:,:),                   & ! surface pressure
    t(:,:,:),                  & ! temperature state
    t_g(:,:),                  & ! surface temperature state
    qv(:,:,:),                 & ! humidity state
    qc(:,:,:),                 & ! cloud water
    qi(:,:,:),                 & ! cloud ice
    qs(:,:,:),                 & ! snow
    h_ice(:,:),                & ! 
    qg(:,:,:)                    ! graupel
  INTEGER   (KIND=iintegers),   INTENT (IN) :: &
    nlev,                      & ! number of model levels
    n_add                        ! number of additional levels for RTTOV
  INTEGER   (KIND=iintegers),   INTENT (OUT) :: &
    stype(:),                  & ! surface type
    wtype(:)                     ! water type of grounde
  LOGICAL,            INTENT (IN), OPTIONAL :: &
    lcloud                       ! whether cloud arrays shall be filled
  REAL   (KIND=wp),   INTENT (OUT) :: &
    pres(:,:),                 & ! presure profiles
    temp(:,:),                 & ! temperature profiles
    humidity(:,:),             & ! humidity profiles
    t2m(:),                    & ! 2m temperature
    hum2m(:),                  & ! 2m humidity
    psurf(:),                  & ! surface pressure
    s_hgt(:),                  & ! surface height
    u10m(:),                   & ! 10m u wind-component
    v10m(:),                   & ! 10m u wind-component
    t_s(:)                       ! surface temperature
  REAL   (KIND=wp),   INTENT (OUT), OPTIONAL :: &
    cloud(:,:,:),              & ! cloud water/ice content
    cfrac(:,:,:),              & ! cloud fraction
    hh(:,:)                      ! height of layers
  INTEGER   (KIND=iintegers),  INTENT (OUT), OPTIONAL :: &
    idg(:),                    & ! ice water cloud scheme
    ish(:)                       ! ice crystal shape
  INTEGER   (KIND=iintegers),   INTENT (OUT), OPTIONAL ::   &
    ierrstat                     ! error status
! Constants
  REAL   (KIND=wp), PARAMETER :: &
    psmax = 1099.9_wp        ! Max surface pressure [hPa]
! Local variables  
  INTEGER (KIND=iintegers) ::  &
    i,j,k,                     & ! indices for gridpoint, level
    iprof,                     & ! indice for profile
    nprof,                     & ! number of profiles
!   n_add,                     & ! additional profiles on model top
    izstat                       ! error status
  REAL (KIND=wp) ::        & !
    fgew   ,fgee   , fgqv   ,  & ! name of statement functions
    ztt    ,zzpv   , zzpa   ,  & ! dummy arguments of stat. func.
    zph, zsigma, zdthdz,       & !
    zsex, zqdw,                & !
    zf_ice,  zcs,              & !
    zt_ice1,zt_ice2,           & !
    zpio   ,zpiu   , zpim   ,  & !
    zthvo  ,zthvu  , zuc    ,  & !
    zclics ,zclws  ,           & !
    zclick ,zclwck , zclwk,    & !
    zclwcs ,                   & !
    zclc,                      & ! Cloud cover in each layer
    zsw,                       & ! Saturation water vapour mixing ratio over water
    zse,                       & ! Saturation water vapour mixing ratio over ice
    zciwc,                     & ! ice mixing ratio
    fact,                      & ! factor for interpolation
    hhlm(ke),                  & ! height at layer centers
    qih,                       & ! cloud ice (plus part of snow and graupel) in g/m^3
    qwh,                       & ! stratiform cloud water in g/m^3
    qwch                         !  convective cloud water

  LOGICAL ::                   & !
    lflag,                     & !
    lcld
    
  REAL    (KIND=wp   ), PARAMETER ::  &
       zepclc = 1.0E-8_wp, & ! avoids cloud cover =1.0 and = 0.0 
       zclwcm = 1.0E-9_wp, & ! avoids cloud water content = 0.0
       tune_qsqg_fac = 0.1,    & ! tuning constant: how much of snow / graupel is transferred to RTTOV
       log_p_top = log(p_top)    ! top of extrapolation above model top

  REAL (KIND=wp ), ALLOCATABLE :: &
       pprof(:,:),             & !pressure profiles at levels
       Tprof(:,:),             & !temperature profile at levels
       qprof(:,:)                !humidity profile at levels


  CHARACTER (LEN=80)     :: yerrmsg
  CHARACTER (LEN=20)     :: yroutine    ! name of the subroutine
    
!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE prepare_rttov_input
!-------------------------------------------------------------------------------

   ! statement function to calculate saturation vapour pressure over water
   fgew(ztt)       = b1 * EXP( b2w*(ztt - b3)/(ztt - b4w) ) ! ztt: temperature

   ! statement function to calculate saturation vapour pressure over ice
   fgee(ztt)       = b1 * EXP( b2i*(ztt - b3)/(ztt - b4i) ) ! ztt: temperature

   ! statement function to calculate specific humitdity  
   fgqv(zzpv,zzpa) = rdv*zzpv/(zzpa - o_m_rdv*zzpv)   ! zzpv: vapour pressure
   ! zzpa: total air pressure

!-------------------------------------------------------------------------------
!- Section 1: Initialize some variables
!-------------------------------------------------------------------------------
   lcld = present(cloud)
   if (present(lcloud)) lcld = (lcloud .and. lcld)

   nprof = SIZE(iob_loc)

   ierrstat = 0

   IF (linterp) THEN
     !---------------------------------------------------------------------------------
     ! - Section 1.1: Interpolate variables to half levels (RTTOV needs p,T,q on  
     !     levels, the cloud variables on layers). Pressure is interpolated linear in 
     !---------------------------------------------------------------------------------
     ALLOCATE(pprof(ke+1,nprof), STAT=izstat)
     ierrstat = ierrstat + izstat
     ALLOCATE(Tprof(ke+1,nprof), STAT=izstat)
     ierrstat = ierrstat + izstat
     ALLOCATE(qprof(ke+1,nprof), STAT=izstat)
     ierrstat = ierrstat + izstat
     IF (ierrstat /= 0) THEN
       yerrmsg =' ERROR  *** allocating space in prepare_rttov_input failed'
       CALL model_abort (my_cart_id,65011,yerrmsg,yroutine)
     ENDIF

     DO iprof = 1,nprof
       i = iob_loc(iprof)
       j = job_loc(iprof)
       DO k = 1,ke
         ! height of layers
         hhlm(k) = 0.5_wp*(hhl(i,j,k+1)+hhl(i,j,k))
       ENDDO
       !IF (.TRUE.) THEN  ! interpolate to half levels??
       DO k=1,ke
         IF (k==1) THEN   !top layer
           fact = (hhlm(k)-hhl(i,j,k))/(hhlm(k+1)-hhlm(k))
           pprof(k,iprof) = LOG(p0(i,j,k)+pp(i,j,k))- fact * (LOG(p0(i,j,k+1)+pp(i,j,k+1))-LOG(p0(i,j,k)+pp(i,j,k)))
           pprof(k,iprof) = EXP(pprof(k,iprof))
           Tprof(k,iprof) = t(i,j,k) - fact * (t(i,j,k+1)-t(i,j,k))
           qprof(k,iprof) = qv(i,j,k) - fact * (qv(i,j,k+1)-qv(i,j,k))
         ELSE   !all layers in between
           fact = (hhl(i,j,k)-hhlm(k-1))/(hhlm(k) - hhlm(k-1))
           pprof(k,iprof) = LOG(p0(i,j,k-1)+pp(i,j,k-1)) + fact * (LOG(p0(i,j,k) + pp(i,j,k)) - LOG(p0(i,j,k-1)+pp(i,j,k-1))) 
           pprof(k,iprof) = EXP(pprof(k,iprof))
           Tprof(k,iprof) = t(i,j,k-1) + fact * (t(i,j,k) -t(i,j,k-1))
           qprof(k,iprof) = qv(i,j,k-1) + fact * (qv(i,j,k) -qv(i,j,k-1))
         ENDIF
       ENDDO
       !ELSE !no interpolation to half levels
     ENDDO
   END IF

   IF (nlev < size(pres,1)) THEN
     ierrstat = 99135
     WRITE(6, '(" ERROR   *** too less levels available (nlev=",I3,", size(pres,1)=",I3,") ***")') nlev, size(pres,1)
     RETURN
   ENDIF

!-------------------------------------------------------------------------------
!- Section 2.1: Fill atmospheric fields
!-------------------------------------------------------------------------------
  DO k=1, ke
    DO iprof = 1,nprof
      i = iob_loc(iprof)
      j = job_loc(iprof)

      IF (linterp) THEN
        pres    (k+n_add,iprof) = pprof(k,iprof)/ 100.0_wp
        humidity(k+n_add,iprof) = qprof(k,iprof)
        temp    (k+n_add,iprof) = tprof(k,iprof)
      ELSE
        pres    (k+n_add,iprof) = (p0(i,j,k) + pp(i,j,k)) / 100.0_wp
        humidity(k+n_add,iprof) = qv(i,j,k)
        temp    (k+n_add,iprof) = t(i,j,k)
      END IF
    ENDDO  
  ENDDO

  IF (PRESENT(hh)) THEN
    DO iprof = 1,nprof
      i = iob_loc(iprof)
      j = job_loc(iprof)
      DO k = 1,ke+1
        ! height of layers
        hh(k,iprof) = hhl(i,j,k)
      ENDDO
    END DO
  ENDIF
  
  DO k = 1,n_add
    DO iprof = 1,nprof
      IF (IAND(extrp_type, extrp_logp) /= 0) THEN
        pres(k,iprof) = exp(log_p_top + &
             (k-1)/REAL(n_add) * (log(pres(n_add+1,iprof)*100) - log_p_top))/ 100._wp
      ELSE
        pres(k,iprof) = (p_top + &
             (k-1)/REAL(n_add) * (pres(n_add+1,iprof)*100 - p_top)) / 100._wp
      END IF
      IF (IAND(extrp_type, extrp_lin) /= 0) THEN
        temp(k,iprof) = temp(n_add + 1, iprof) + (n_add - k + 1) * &
                        (temp(n_add+1, iprof) - temp(n_add+2, iprof))
        humidity(k,iprof) = humidity(n_add + 1, iprof) + (n_add - k + 1) * &
                            (humidity(n_add+1, iprof) - humidity(n_add+2, iprof))
      ELSEIF (IAND(extrp_type, extrp_clim) /= 0) THEN
        temp(k,iprof) = temp(n_add + 1, iprof) + (n_add - k + 1) * &
                        (t_top - temp(n_add+1, iprof)) / REAL(n_add)
        
        humidity(k,iprof) = humidity(n_add + 1, iprof) + (n_add - k + 1) * &
                            (q_top - humidity(n_add+1, iprof)) / REAL(n_add)
      ELSEIF (IAND(extrp_type, extrp_const) /= 0) THEN
        temp(k,iprof)     = temp(n_add + 1, iprof)
        humidity(k,iprof) = humidity(n_add + 1, iprof)
      ENDIF

      temp(k,iprof) = MIN(tmax,MAX(tmin,temp(k,iprof)))
      humidity(k,iprof) = MAX(qmin,MIN(qmax,                               &
                                       q_sat(pres(k,iprof), temp(k,iprof)),&
                                       humidity(k,iprof)))
    ENDDO  
  ENDDO

!-------------------------------------------------------------------------------
!- Section 2.2: Fill Surface information
!-------------------------------------------------------------------------------

  lflag = ALL(SHAPE(h_ice) > 0)


  DO iprof = 1,nprof
    i = iob_loc(iprof)
    j = job_loc(iprof)
    
    psurf(iprof) = ps(i,j) / 100.0_wp
    psurf(iprof) = MIN(psurf(iprof), psmax)

    pres(nlev,iprof)     = psurf(iprof)
    humidity(nlev,iprof) = qv_2m(i,j)
    temp(nlev,iprof)     = t_2m(i,j)

    t2m(iprof)   = t_2m(i,j)
    hum2m(iprof) = MIN(qv_2m(i,j), q_sat(psurf(iprof), t2m(iprof)))

    s_hgt(iprof)         = hsurf(i,j) * 0.001_wp
    u10m(iprof)          = u_10m(i,j)
    v10m(iprof)          = v_10m(i,j)
    t_s(iprof)           = t_g(i,j)

    IF (soiltyp(i,j) == 9) THEN
      stype(iprof) = 1
      IF (lflag) THEN
        IF (h_ice(i,j) > 0.01_wp) stype(iprof) = 2
      ENDIF  
    ELSE IF (soiltyp(i,j) == 10) THEN   
      stype(iprof) = 2
    ELSE  
      stype(iprof) = 0
    ENDIF

    IF (lseaice .AND. stype(iprof) /= 0) THEN
      IF (lseamask(i,j)) THEN
        wtype(iprof) = 1
      ELSE
        wtype(iprof) = 0
      ENDIF
    ELSE  
      wtype(iprof) = 0
    ENDIF  
  ENDDO

!-------------------------------------------------------------------------------
!- Section 2.3: Convert Humidity
!-------------------------------------------------------------------------------
  humidity(1:nlev,:) = humidity(1:nlev,:) / (1.0_wp - humidity(1:nlev,:)) &
                       * rcnw
  hum2m(:) = hum2m(:) / (1.0_wp - hum2m(:)) * rcnw

!-------------------------------------------------------------------------------
!- Section 2.4: Compute Cloud Cover and Contents
!-------------------------------------------------------------------------------
  IF (lcld) THEN
!CDIR COLLAPSE  
    cloud(:,:,:) = 0._wp
!CDIR COLLAPSE  
    cfrac(:,:,:) = 0._wp

    idg  (:)     = iwc2effdiam
    ish  (:)     = iceshape

    lflag = ALL(SHAPE(qg) > 0)

    DO k=1, ke
      DO iprof = 1,nprof
        i = iob_loc(iprof)
        j = job_loc(iprof)

        !<ASb
        !  Here similar implementation as in src_radiation to obtain total 
        !  (grid-scale plus subgrid-scale) water- and ice content 

        ! Critical relative humidity as function of thermal stability
        zph      = p0(i,j,k) + pp(i,j,k)
        zsigma   = zph / ps(i,j)
        zdthdz   = 0.0
        ! specific humidity at saturation
        ! over water (zsw ) and ice (zse)
        zsw      = fgqv ( fgew(t(i,j,k)), zph)
        zse      = fgqv ( fgee(t(i,j,k)), zph)
        zt_ice1= t0_melt -  5.0_wp
        zt_ice2= t0_melt - 25.0_wp

        zsex      = zsw
        zqdw      = qv(i,j,k) + qc(i,j,k)
        IF (lprog_qi) THEN
          zf_ice      = 1.0_wp - MIN( 1.0_wp, MAX( 0.0_wp, &
               (t(i,j,k)-zt_ice2)/(zt_ice1-zt_ice2) ) )
          zqdw        = zqdw      + qi(i,j,k)
          zsex        = zsw * (1.0_wp - zf_ice) + zse*zf_ice
        ENDIF

        IF(k == ke) THEN
          zpio    = ( 1.0E-5 *( p0(i,j,k)+pp(i,j,k) ) )**rdocp
          zpiu    = ( 1.0E-5 * ps(i,j)  )**rdocp
          zpim    = 0.5*(zpio+zpiu)
          zthvo   = t  (i,j,k  )/zpio
          zthvu   = t_g(i,j    )/zpiu
          zdthdz  = zthvo - zthvu
        ELSE IF(zsigma.GT.0.95) THEN
          zpio    = ( 1.0E-5 *( p0(i,j,k  )+pp(i,j,k  ) ) )**rdocp
          zpiu    = ( 1.0E-5 *( p0(i,j,k+1)+pp(i,j,k+1) ) )**rdocp
          zpim    = 0.5*(zpio+zpiu)
          zthvo   = t(i,j,k  )/zpio
          zthvu   = t(i,j,k+1)/zpiu
          zdthdz  = zthvo - zthvu + (lh_v*cpdr/zpim)*(qv(i,j,k)-qv(i,j,k+1))
        ENDIF

        ! subgrid scale cloud cover as function of relative humidity
        zuc     = 0.95 - uc1*zsigma*(1.-zsigma)*(1.+uc2*(zsigma-0.5))
        zcs     = MAX ( 0.0_wp,                &
             MIN ( 1.0_wp, (zqdw/zsex-zuc)/(ucl-zuc) ) )**2

        ! Corrections and limitations
        IF ( (zsigma > 0.95_wp) .AND. (zdthdz < 0.0_wp) ) THEN
          zcs = 0.0_wp  ! no cloud cover in unstable stratification
        ENDIF
        IF ( qc(i,j,k) > 0.0_wp ) THEN  ! grid scale clouds
          IF ( llandmask(i,j) .AND. k < ke ) zcs = 1.0_wp
        ENDIF
        IF (lprog_qi) THEN
          IF (qi(i,j,k) > 1.0E-7_wp) THEN
            zcs = 1.0_wp ! grid scale clouds with cloud ice
          ENDIF
        ENDIF

        ! Maximum in-cloud water content:  1.0% of specific humidity at saturation
        !                                  except for convective clouds (fixed)
        ! Standard diagnosis

        IF (lprog_qi) THEN
          zclws  = 0.005*zsex
          zclwcs = zclws*(1.0_wp-zf_ice)
          zclics = zclws*zf_ice

          ! Check for grid-scale water or ice-clouds
          ! Now change in zclwcs only if qc(i,j,k,nzx) > 0.0
          IF ( qc(i,j,k) > 0.0_wp ) THEN      ! grid scale cloud water
            zclwcs = MAX( zclwcs, 0.5_wp*qc(i,j,k) )
          ENDIF
          ! Now change in zclics only if qi(i,j,k,nzx) > 1.0E-7
          IF ( qi(i,j,k) > 1.0E-7_wp ) THEN  ! grid scale cloud ice
            zclics = MAX( zclics, 0.5_wp*qi(i,j,k) )
          ENDIF

          ! Convective cloud water / ice content
          zclwk  = MAX( 2.0_wp*zclws, 0.0002_wp )
          zclwck = zclwk*(1.0_wp-zf_ice)
          zclick = zclwk*zf_ice

          ! Reduce the cloud cover of ice clouds in the upper troposphere
          ! for the diagnosis of clch and clct
          ! Changed by Axel: 
          IF ( (zclwcs <= 1.0E-10_wp) .AND. (zclics  > 0.0_wp) ) THEN
            zcs = zcs*MIN( 1._wp, MAX(0.0_wp, &
                  ( LOG(zclics)       - LOG(1.E-7_wp) )/ &
                  ( LOG(8.E-6_wp) - LOG(1.E-7_wp) )) )
          ENDIF

          zciwc = zclick*clc_con(i,j,k) + &
               zclics*zcs*(1.0_wp-clc_con(i,j,k))

        ELSE  !no prognostic qi

          zclwcs = 0.005*zsw
          zclwck = MAX( zclwcs, 0.0002_wp )
          IF ( qc(i,j,k) > 0.0 ) THEN  ! grid scale clouds
            zclwcs = MAX( zclwcs, 0.5*qc(i,j,k) )
          ENDIF
          ! set area average cloud ice content (in-cloud)
          zciwc = 0.0_wp

        ENDIF   !lprog_qi

        ! calculate combined cloud cover
        zclc   = zcs + clc_con(i,j,k)*( 1.0_wp - zcs )

        IF (zclc > 0._wp) THEN   ! if clouds present

          ! Restrictions for radiative calculations, conversion to [g/m^3]
          ! ------------------------------------
          zciwc = MAX( zclwcm, zciwc )
          zclwcs = MAX( zclwcm, zclwcs )
          zclwck = MAX( zclwcm, zclwck )

          IF (lflag) THEN  ! qg available
            ! also consider radiative effects of snow and graupel:
            qih =  (zciwc + tune_qsqg_fac * (qs(i,j,k)+qg(i,j,k))) * 1000._wp * rho(i,j,k)
          ELSE             ! qg not available
            qih =  (zciwc + tune_qsqg_fac * (qs(i,j,k)          )) * 1000._wp * rho(i,j,k)
          ENDIF

          ! More than 10g/m^3 ice/snow is physically not reasonable
          ! and might cause a crash of RTTOV
          qih = min(qih,10._wp)
          qwh = zclwcs * 1000._wp * rho(i,j,k)   ! stratiform cloud water
          qwch = zclwck * 1000._wp * rho(i,j,k)  ! convective cloud water
          zclc  = MAX( zepclc, MIN(1.0_wp-zepclc,zclc) )

          ! Fill RTTOV structures

#ifdef RTTOV10
          IF (rttov_ifc_version >= 10) THEN
            ! RTTOV10 can handle multiple cloud types at the same time
            cloud(6,k+n_add,iprof) = qih
            cloud(3,k+n_add,iprof) = qwch
            cloud(1,k+n_add,iprof) = qwh
            ! the cloud fraction must be specified once only - not type
            ! specific
            cfrac(1,k+n_add,iprof) = zclc
#else
          IF (.false.) THEN
#endif          
          ELSEIF ( (qwh) <= 0._wp .AND. qih > 0._wp) THEN
            !Cirrus clouds
            cloud(6,k+n_add,iprof) = qih
            cfrac(6,k+n_add,iprof) = zclc
          ELSEIF ( qwh <= 0._wp   .AND. clc_con(i,j,k) > 0._wp) THEN
            !convective clouds
            cloud(3,k+n_add,iprof) = qwch
            cfrac(3,k+n_add,iprof) = zclc
          ELSEIF ( (qwh + qih ) > 0._wp) THEN
            !stratiform clouds
            cloud(1,k+n_add,iprof) = qwh + qwch + qih
            cfrac(1,k+n_add,iprof) = zclc
          ENDIF
        ENDIF ! ! (zclc > 0._wp)
      ENDDO  !nprof
    ENDDO !ke
  ENDIF ! present(cloud)

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE prepare_rttov_input

!==============================================================================

! this ifdef holds for SR cloud_read_prep, grb_open, grb_read
#ifdef NUDGING
!==============================================================================
!+ Module procedure in "src_cloudana" preparing cloud top observations
!-------------------------------------------------------------------------------

SUBROUTINE cloud_read_prep(yactdate,yaction,ctype,cth,ierr)

!-------------------------------------------------------------------------------
!
! Description:
!   This procedure of the module "src_obs_rad" reads in the cloud type and cloud top height
!    (NWCSAF SEVIRI product mapped on COSMO grid)
! 
! Method:
!   1: Read the cloud data from binary file 
!   2: Distribute the fields with new cloud type to different PE's
!
!-------------------------------------------------------------------------------

! Subroutine / Function arguments:
!----------------------------------

CHARACTER (LEN=12), INTENT(IN)        ::  yactdate

CHARACTER (LEN= *), INTENT(IN)        ::  yaction      ! action to be performed

INTEGER (KIND=iintegers), INTENT(OUT) :: ierr

REAL (KIND=wp), INTENT(OUT)       ::    &
                         ctype(ie,je),      &  ! cloud type field
                         cth(ie,je)            ! cloud top height

! Local parameters, scalars, arrays :
!-------------------------------------------------------------------------------
! Local scalars:

CHARACTER (LEN=20)    :: yroutine    ! name of the subroutine

CHARACTER (LEN=80)    :: yerrmsg

CHARACTER (LEN=100)   :: yzfile_ct, yzfile_cth


!---------------
! Variables for input of observational data
! 
INTEGER (KIND=iintegers), SAVE             ::       &
     nugrib_ct, nugrib_cth    ! unit number for cloudana GRIB file

! 
INTEGER (KIND=iintegers)                   ::       &
       izlocstat     ,&   ! status variable
       istat         ,&     ! --- "---------
       ntimes        ,&   ! how often will the cloudanalysis be carried out
       ierr_p        ,&  ! error flag for communication
       mocc          ,&  ! current century
       moyy          ,&  ! current year
       momm          ,&  ! current month
       modd          ,&  ! current day
       mohh          ,&  ! current hour
       momi          ,&  ! current minute
       moss          ,&  ! current second
       iwlength      ,&
       ibits

REAL (KIND=wp)   ::   &     !global fields
          ctype_tot(ie_tot,je_tot),  &  
          cth_tot(ie_tot, je_tot)       

!- End of header
!===============================================================================

!-------------------------------------------------------------------------------
! Begin of subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Section 0: Initialization of some variables
!           - Get start date and time of this model run
!           - Open data file and read general header information
!           - Allocate space for temporal fields
!
!-------------------------------------------------------------------------------

ibits       = BIT_SIZE (1_intgribf)
iwlength    = INT (ibits, iintegers) / 8


  yroutine='cloud_read_prep'
  ierr = 0
  ierr_p = 0
  istat = 0

  ctype_tot(:,:)= rmdi
  cth_tot(:,:)  = rmdi


  !Some initializations
  lds = ie_tot*je_tot
  lbm = 0_iintegers
  lfd = lds * 2_iintegers /iwlength + 5000_iintegers
  READ(yactdate,'(7I2)') mocc,moyy, momm, modd, mohh, momi, moss 
        
  !Filenames
  yzfile_ct =  yclouddir(1:LEN_TRIM(yclouddir))//'/ct_nwcsaf_'//ydate_ini(1:12)//'.grb2'
  yzfile_cth = yclouddir(1:LEN_TRIM(yclouddir))//'/cth_nwcsaf_'//ydate_ini(1:12)//'.grb2'

  ! Open data file and allocations in first time step
  IF (yaction == 'input-satpp') THEN  !first call: open file, allocations
    IF (my_cart_id == 0) THEN  !Reading is only done by Processor 0
       ALLOCATE ( ds(lds), STAT=izlocstat )
       istat = istat + izlocstat
       IF (istat /= 0) THEN
          yerrmsg =' ERROR  *** allocating space for reading cloud product grib files failed'
          CALL model_abort (my_cart_id,6501,yerrmsg,yroutine)
       ENDIF
       !open cloud type file

       CALL grb_open (yzfile_ct(1:LEN_TRIM(yzfile_ct)),nugrib_ct,ierr)
       IF (ierr /= 0 .OR. nugrib_ct == 0) THEN
          PRINT*, 'WARNING: Error when opening cloud type data file.'
         ! set complete field ctype_tot to missing values
          ctype_tot(:,:)=rmdi
          PRINT*, 'Complete field ctype_tot has been set to missval.'
          RETURN
        ENDIF
       !open cloud top height file
        CALL grb_open (yzfile_cth(1:LEN_TRIM(yzfile_cth)),nugrib_cth,ierr)
        !===================
       IF (ierr /= 0 .OR. nugrib_cth == 0) THEN
          PRINT*, 'WARNING: Error when opening cloud top height data file.'
         ! set complete field ctype_tot to missing values
         ctype_tot(:,:)=rmdi
         PRINT*, 'Complete field cth_tot has been set to missval.'
         RETURN
       ENDIF
    ENDIF !my_cart id
  ELSE IF (yaction == 'compute') THEN 
    IF (my_cart_id == 0) THEN  !Reading is only done by Processor 0
       CALL grb_read(nugrib_ct,mocc,moyy,momm,modd,mohh,momi,ctype_tot,ierr)
        !======================
       IF (ierr /= 0 .OR. nugrib_ct == 0) PRINT*, 'WARNING: Error when reading cloud type data.'
       CALL grb_read(nugrib_cth,mocc,moyy,momm,modd,mohh,momi,cth_tot,ierr)
        !======================
       IF (ierr /= 0 .OR. nugrib_cth == 0) PRINT*, 'WARNING: Error when reading cloud top height data.'
        !forward information to all processors
     ENDIF !my_cart_id
     CALL distribute_values (ierr,1,0,imp_integers,icomm_cart,ierr_p)
     IF (num_compute > 1) THEN
        CALL distribute_field (ctype_tot(1:ie_tot,1:je_tot),ie_tot,je_tot,ctype(1:ie,1:je),ie,je,0,ierr_p)
            !======================
        CALL distribute_field (cth_tot(1:ie_tot,1:je_tot), ie_tot,je_tot,cth(1:ie,1:je),ie,je,0,ierr_p)
          !======================
      ELSE ! we are running on one PE
         ctype(:,:) = ctype_tot(:,:)
         cth(:,:)   = cth_tot(:,:)
      ENDIF
      
  ELSE IF (yaction == 'cleanup') THEN 
#ifdef GRIBAPI
  !--  CLEANUP --
       IF (my_cart_id == 0) THEN
          IF (nugrib_ct > 0 .AND. yaction == 'cleanup') THEN
             PRINT*, '+++ Close cloud type grib file +++'
             CALL grib_close_file (nugrib_ct)
              !============
          ENDIF
          IF (nugrib_cth > 0 .AND. yaction == 'cleanup') THEN   
             CALL grib_close_file (nugrib_cth)
              !============
          ENDIF
          DEALLOCATE(ds)
        ENDIF !my_cart_id
#endif
  ENDIF !yaction

!-------------------------------------------------------------------------------
! End of subroutine 
!-------------------------------------------------------------------------------
   
END SUBROUTINE cloud_read_prep

!===============================================================================
!+ Module procedure in "cloudana" reading cloud-top info from input file
!-------------------------------------------------------------------------------
 
SUBROUTINE grb_open (yfile,nugrib,ierr)

!-------------------------------------------------------------------------------
!
! Description:
!   This procedure opens the cloud data input file and reads the general header information
!  
! Method:
!   Call of GRIB I/O routines of DWDLIB.
!
! Input files:
!   GRIB file.
!   
!-------------------------------------------------------------------------------

! Subroutine scalars (intent: out):
!----------------------------------

CHARACTER (LEN=*), INTENT(IN)                   ::       &
     yfile

INTEGER (KIND=iintegers), INTENT(OUT)            ::       &
     nugrib        ,& ! unit number for GRIB file
     ierr             ! status/error code

INTEGER (KIND=intgribc)  ::  &
     nudatc,  ierrc ! variables for C routines

! Local scalars, arrays:
!-------------------------------------------------------------------------------

CHARACTER (LEN=14)     :: yroutine
CHARACTER (LEN=50)     :: yerrmsg

!- End of header
!===============================================================================

!-------------------------------------------------------------------------------
! Begin of subroutine
!-------------------------------------------------------------------------------

  yroutine = 'grb_open'
  ierr = 0
  ierrc = 0

  ! Open input file and read the header information
  !-------------------------------------------------------------------------------
  ! open input file
  nugrib=0_iintegers

  ! open GRIB file
  nudatc=INT(nugrib, intgribc)   


#ifdef GRIBAPI
  CALL grib_open_file(nudatc, yfile, 'r', ierrc)
#endif

  nugrib=INT(nudatc, iintegers)
  PRINT *,'opened input file: cloud data with unit number: ', nugrib
  ierr=INT(ierrc, iintegers)
  IF (ierr /= 0) THEN
    yerrmsg='Error in opening the satellite cloud product data file.'
    CALL model_abort(my_cart_id,6010,yerrmsg,yroutine)
  ELSE
    PRINT*, 'Opened ', yfile ,' successfully '
  ENDIF

!-------------------------------------------------------------------------------
! End of subroutine 
!-------------------------------------------------------------------------------

END SUBROUTINE grb_open

!===============================================================================
!+ Module procedure in "src_obs_rad" reading cloud products on COSMO grid
!-------------------------------------------------------------------------------

SUBROUTINE grb_read (nugrib,mocc,moyy,momm,modd,mohh,momin,datafield,ierr)

!-------------------------------------------------------------------------------
!
! Description:
!   This procedure is called by "cloud_obs_prep" and
!   reads data from the cloud type grib file.
!
! Method:
!   Simple read from sequential unformatted file.
!
! Input files:
!   Sequential unformatted file 
!
!-------------------------------------------------------------------------------
! 
! Subroutine arguments: 
!-------------------------------------------------------------------------------
! Subroutine scalars (intent: in):
!-----------------------------------

INTEGER (KIND=iintegers), INTENT(INOUT)          ::       &
     nugrib          ! unit number for input GRIB file

! Subroutine scalars (intent: out):
!----------------------------------

INTEGER (KIND=iintegers), INTENT(IN)            ::       &
     mocc         ,& ! date of requested observations : century
     moyy         ,& ! date of requested observations : year
     momm         ,& ! date of requested observations : month
     modd         ,& ! date of requested observations : day
     mohh         ,& ! date of requested observations : hour
     momin           ! date of requested observations : minute

INTEGER (KIND=iintegers), INTENT(OUT)           ::       &
     ierr             ! status/error code

! Subroutine arrays (intent: out):          
!---------------------------------

REAL (KIND=wp), INTENT(OUT)                   ::       &
     datafield(ie_tot,je_tot)        ! field for decoded data

! Local parameters, scalars, arrays :
!-------------------------------------------------------------------------------
! Local scalars :
!----------------

CHARACTER (LEN=14)     :: yroutine

INTEGER (KIND=iintegers)             ::       &
     ilen             !length of grib record

INTEGER (KIND=iintegers)         ::       &
     igrib       ,& ! grib handle
     inrpoints   ,& ! number of points in GRIB record
     icc         ,& ! date of observations : century
     iyy         ,& ! date of observations : year
     imm         ,& ! date of observations : month
     idd         ,& ! date of observations : day
     ihh         ,& ! date of observations : hour
     imin        ,& ! date of observations : minute
     iwlength    ,&
     ibits

! ! variable for GRIB-input
INTEGER (KIND=intgribc)  ::  &
    ierrc ! variables for C routines


CHARACTER (LEN=20) :: ydate
CHARACTER (LEN=20) :: ytime

!- End of header
!===============================================================================
#ifdef GRIBAPI
!-------------------------------------------------------------------------------
! Begin of subroutine
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Section 0: Initialization of some variables
!-------------------------------------------------------------------------------

  yroutine = 'grb_read'
  PRINT*, 'Beginning of routine grb_read at ntstep', ntstep
  ! 
  datafield = rmdi

  ibits     = BIT_SIZE (1_intgribf)
  iwlength  = INT (ibits, iintegers) / 8

  ierr = 5

  readloop: DO !until no more record

    ! read the record with grib api. Build grib handler igrib
    CALL grib_new_from_file(nugrib, igrib,ierr)

    !      

    IF (ierr /= GRIB_SUCCESS) THEN
       IF (ierr == GRIB_END_OF_FILE) THEN
       !end of file reached
          PRINT*,'end of cloud obs file reached and no cloud data found: ntstep: ', ntstep
        ENDIF
       CALL grib_release(igrib)
       RETURN
    ELSE
      PRINT*, 'read record successfully '
       !Get meta data from grib handle
       CALL grib_get_size(igrib, 'values', inrpoints, ierr)
       IF ( ierr /= 0 ) THEN
          PRINT*, 'ERROR while getting size of grib file '
       ENDIF
       IF (inrpoints > lds) THEN
          PRINT*, ' ERROR: size of message is too big for allocated field: ', inrpoints, ' > ', lds
       ENDIF
       ! Get metadata for checking the date and time
       CALL grib_get (igrib, 'dataDate', ydate, ierr)
       CALL grib_get (igrib, 'dataTime', ytime, ierr)
       PRINT*, 'ydate ', ydate, 'ytime ', ytime

        ! check for correct date and time
       read( ydate(1:2), '(i2)' ) icc
       read(ydate(3:4),'(i2)') iyy
       read(ydate(5:6),'(i2)') imm
       read(ydate(7:8),'(i2)') idd
       read(ytime(1:2),'(i2)') ihh
       read(ytime(3:4),'(i2)') imin
       !PRINT*, '+++READ GRIB File with time info :', icc, iyy, imm, idd, ihh, imin
       !PRINT*, '### looking for time stamp:   ', mocc, moyy, momm, modd, mohh, momin

       IF (icc==mocc .AND. iyy==moyy .AND. imm==momm .AND. idd==modd .AND. ihh==mohh .AND. &
           (ABS(imin-momin) < 5)) THEN
       ! get field from grib file
         CALL grib_get(igrib,'values',ds)
         datafield=RESHAPE(REAL(ds,wp),(/ie_tot,je_tot/))
         !PRINT*, 'Found correct timestamp in data: ', icc, iyy, imm, idd, ihh, imin
         ierr = 0
         CALL grib_release(igrib)
         RETURN
       ENDIF   ! matching time stamp
       CALL grib_release(igrib)
    ENDIF  !ierr

  ENDDO readloop
#endif
!------------------------------------------------------------------------------
! End of subroutine 
!------------------------------------------------------------------------------

END SUBROUTINE grb_read

! this ends ifdef NUDGING
#endif

!==============================================================================

!------------------------------------------------------------------------------
! End of module src_obs_rad
!------------------------------------------------------------------------------

END MODULE src_obs_rad
