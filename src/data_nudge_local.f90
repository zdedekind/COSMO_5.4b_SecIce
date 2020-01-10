!+ Data module for the variables only used to compute the local information
!-------------------------------------------------------------------------------

MODULE data_nudge_local

!-------------------------------------------------------------------------------
!
! Description:
!   This module contains those variables (scalars and arrays) which are merely
!   concerned with the computation of all the required 'local' information,
!   and which are not already contained in 'data_obs_qc_limits.f90'.
!   The 'local information' consists of observation increments, and further
!   parameters on the observations and their location used for the nudging.
!   This information will later be broadcast to other processors for the
!   spreading of the observational information.
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.10       1998/09/29 Christoph Schraff
!  Initial release
! 1.13       1998/10/22 Christoph Schraff
!  Removal of some temporal varialbes.
! 1.15       1998/11/02 Christoph Schraff
!  Global allocatable arrays moved from modules src_mult_local, src_sing_local.
! 1.19       1998/12/11 Christoph Schraff
!  Tables for quality control thresholds modified.
! 1.38       2000/04/06 Christoph Schraff
!  Radiosonde height thresholds reduced to values commonly used for US-network.
!  Radiosonde height error correlation matrix added. 'zvispro', 'zvifaco' added.
! 2.4        2001/01/29 Christoph Schraff
!  Variance / threshold tables for wind adjusted to values used at ECMWF.
!  Introduction of several additional factors to thresholds for selected obser-
!  vations following ECMWF, and to surface pressure thresholds used for the
!  spatial consistency check.
! 2.13       2002/01/18 Christoph Schraff
!  To allow for the 'hot' compiler option on IBM: Data statements replaced by
!  direct assignment, and parameter attribute added to corresponding variables.
! 3.3        2003/04/22 Maria Tomassini + Christoph Schraff
!  Administration arrays for GPS reports introduced.
! 3.12       2004/09/15 Christoph Schraff
!  Extension to (prepare to) include assimilation of satellite retrievals.
!  Introduction of several limits / errors for revised optional humidity
!  thresholds. Introduction of analysis layers and standard layers for revised
!  optional multi-level check.
! 3.18       2006/03/03 Christoph Schraff
!  Introduction of further factors for revised optional humidity thresholds.
!  Increased basic humidity quality control thresholds ( = new GME-OI values).
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_22        2012/01/31 Christoph Schraff
!  Introduction of SOR (Simulated Observation Record).
! V4_28        2013/07/12 Christoph Schraff
!  SOR moved to 'data_obs_record.f90'; quality control thresholds moved to new
!  module 'data_obs_qc_limits.f90'. Some other variables moved inside
!  'src_mult_local.f90'.
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!===============================================================================
!
! Modules used:
!
!-------------------------------------------------------------------------------

USE data_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!-------------------------------------------------------------------------------

IMPLICIT NONE

!===============================================================================

! Local Declarations:

!-------------------------------------------------------------------------------
! Global (i.e. public) Declarations:
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Section 1:  General parameters and other variables
!-------------------------------------------------------------------------------

!      1.1  Parameters
!      ---------------

  REAL (KIND=wp)           , PARAMETER  :: &
    ptropop = 22500.0_wp ,& ! pressure at tropopause (standard atmosphere)
    thdzt   =  0.0041_wp ,& ! mean vert. pot. temp. gradient in troposphere
    fsvcut  =     0.9_wp    ! fraction with which the adjusted scale S_adj,
                                ! rather than the scale S_nl as specified in the
                                ! namelist, is used for the determination of the
                                ! vertical cut-off in log(p) or theta(e)-units:
                                !   rcutoff = fcutoff * (     fsvcut * S_adj
                                !                        + (1-fsvcut)* S_nl )


!      1.2  Other arrays and variables
!      -------------------------------

  REAL (KIND=wp)           , ALLOCATABLE  :: &
    wtml     (:,:) ,& ! temporal nudging weights for multi-level reports
    wtsg     (:,:) ,& ! temporal nudging weights for single-level reports
    wtgp     (:,:) ,& ! temporal nudging weights for GPS reports
    wttv     (:,:) ,& ! temporal nudging weights for satellite retrievals
                      ! format of time administrator:
                      ! index 1: 2nd part of nudging time window for 1st report
                      ! index 2: 1st part of nudging time window for 2nd report
                      ! index 3: 1st part of nudging time window for 2nd report
    tmladm   (:,:) ,& ! administrator of multi-level ODR
    tsgadm   (:,:) ,& ! administrator of single-level ODR
    tgpadm   (:,:) ,& ! administrator of GPS ODR
    ttvadm   (:,:) ,& ! administrator of satellite retrievals
    zsdni    (:,:)    ! equidistant pts. used for non-isotropic horiz. correlat.

  INTEGER (KIND=iintegers) , ALLOCATABLE  :: &
                      ! format of integer administrator:
                      ! index 1: ODR index of first  (past)   obs report
                      ! index 2: ODR index of second (future) obs report
                      ! index 3: longitudinal grid point assigned to report
                      ! index 4: latitudinal  grid point assigned to report
                      ! index 5: non-zero ODR index of one of the 2 reports
                      ! index 6: = 1: linear temporal interpolation; = 0: else
                      ! index 7: index 'icps' of surface pressure incr. array
                      ! index 8: index 'icsu / icua' of single-lev. incr. array
    imladm   (:,:) ,& ! administrator of multi-level ODR
    isgadm   (:,:) ,& ! administrator of single-level ODR
    igpadm   (:,:) ,& ! administrator of GPS ODR
    itvadm   (:,:) ,& ! administrator of satellite retrievals
    ilvz     (:)      ! index of (surface) level of multi-level report for which
                      !   simulated height is written to 'smlbdy' in
                      !   'ps_obs_operator_mult' (not in 'mult_obs_operator_z')

  INTEGER (KIND=iintegers)  ::    &
!   nmasta         ,& ! total number of sorted multi-level aircraft reports
    nmlsta         ,& ! total number of sorted multi-level stations
    ngpsta         ,& ! total number of sorted GPS stations
    ntvsta         ,& ! total number of sorted satellite retrieval 'stations'
    kml300            ! index of full model level corresponding to about 300 hPa

  LOGICAL                   ::    &
    lqcall            ! .TRUE if threshold quality control for all observations
                      !          at current timestep

!-------------------------------------------------------------------------------

END MODULE data_nudge_local
