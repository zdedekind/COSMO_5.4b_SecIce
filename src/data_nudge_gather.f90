!+ Data module for the variables used generally for the nudging
!-------------------------------------------------------------------------------

MODULE data_nudge_gather

!-------------------------------------------------------------------------------
!
! Description:
!   This module contains the variables (scalars and arrays) which are both
!   concerned with the computation and broadcasting to other processors of
!   all the required 'local' information (i.e. observation increments, and
!   further parameters on the observations and their location) and the
!   spreading of the observation increments. These are
!    - some commonly used parameters (constants) and general variables
!    - the arrays containing the 'local' information to be broadcast / gathered
!    - the format of the observation increment array for multi-level reports
!    - arrays on local informations for printing for control
!    - variables defining the size of the arrays
!
!   The arrays are declared as allocatable arrays and are allocated in
!   procedure 'local_info_aux' and deallocated in 'gather_spread_aux'.
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
!  Variable 'acthr' moved to 'data_nudge_all'.
! 1.36       2000/02/24 Michael Buchhold
!  Introduction of ACAR aircraft reports
! 2.4        2001/01/29 Christoph Schraff
!  Addition of local information and geometrical fields used for the spatial 
!  consistency check of surface pressure observations.
! 2.13       2002/01/18 Christoph Schraff
!  To allow for the 'hot' compiler option on IBM: Data statements replaced by
!  direct assignment, and parameter attribute added to corresponding variables.
! 2.19       2002/10/24 Christoph Schraff
!  Variable 'zpblsu' introduced, variable 'zp0ops' removed.
! 3.18       2006/03/03 Christoph Schraff
!  Addition of local information used for spatial consistency check of IWV.
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_22        2012/01/31 Christoph Schraff
!  Introduction of 'kcdtyml' and 'maxmloi_tot' for more consistent array sizes.
!  ANSI checking: Avoid more than 39 continuation lines.
!  Variable 'kobtysu' introduced (for scatterometer), and also 'kcdtyua'
!  'kcdtysu' (for separate net increments from different sets of obs systems).
! V4_28        2013/07/12 Christoph Schraff
!  Variable 'zoips_b' introduced for LBC QC checks.
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
! Global (i.e. public) Declarations:
!-------------------------------------------------------------------------------

! 1. Parameters and general variables
! -----------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    ilstidg = 9              ,& ! char. length used for gathering the station ID
    ilstidp = 9              ,& ! char. length used for printing the station ID
                                ! Note: (ilstid >= ilstidg >= ilstidp) !!!
                                !        ilstid is defined in 'data_obs_record'
    mxispr  = 42                ! number of equidistant vertical points (used
                                ! for non-isotropic horizontal correlations)

  REAL (KIND = wp)         , PARAMETER  :: &
    zalllow = -400.0_wp  ,& ! smallest allowed height within model domain
    zallhig = 100000.0_wp,& ! largest  allowed height within model domain
    p0r     = 100000.0_wp,& ! reference pressure for potential temperature
    qst     = 4.0_wp     ,& ! quotient of the mean vertical potential temp-
                                ! erature gradient in the stratosphere to the
                                ! analogous gradient (thdzt) in the troposphere
    xezispr = 1.05_wp    ,& !   / parameters used to
    dezispr = -0.025_wp  ,& !  /  define the vertically
    sezispr = 8000.0_wp  ,& ! <   equidistant points
    xthispr = 242.5_wp   ,& !  \  (for non-isotropic
    dthispr = 2.5_wp     ,& !   \ horizontal correlations)
    vcutnit = 2499.9_wp  ,& ! horiz. correlations are non-isotropic if the
    vcutnip = 99.9_wp       ! \  vertical scales  'rdsprni' < 'zcutnit,p'

  INTEGER (KIND = iintegers)          :: &
    ktp            ,& ! lowermost purely horizontal model main level
    ktth              ! top model level with spreading along isentropic surfaces

  CHARACTER (LEN=ilstidp) ystid      ! obs. station identity to be printed


! 2. All the required local information on the observations and their
!    location which is broadcast to other processors for the spreading 
! --------------------------------------------------------------------

  REAL    (KIND = wp)       , ALLOCATABLE :: &
! local information on multi-level reports
    oiml     (:,:) ,& ! observation increment record OIR for multi-level reports
    zwtml  (:,:,:) ,& ! temporal weights
    zspobml(:,:,:) ,& ! spreading parameter at the base / top of the profile
    fcorlml(:,:,:) ,& ! reduction to vertical correl. at the profile's base /top
    zvcutml(:,:,:) ,& ! vertical cut-off at the base / top of the profile
    zrtdgml(:,:,:) ,& ! convertor from height to pressure units for vertical
                      ! correlations at the base / top of the profile
    zpobml (:,:,:) ,& ! pressure (p) at the base / top of the profile
    zemkpml(:,:,:) ,& ! EXP( -R/cp *p ) at the base / top of the profile
    zriflml(:,:,:) ,& ! upper estimate to horizontal radius of influence
    zstpml     (:) ,& ! spreading parameter at the tropopause & obs. location
    zsprml   (:,:) ,& ! spreading parameter at model levels & obs. location
    zlopml   (:,:) ,& ! log( pressure ) at model levels & obs. location
    zrtgkml  (:,:) ,& ! convertor from height to pressure units for vertical
                      ! correlations at model levels & obs. location
    znisml   (:,:) ,& ! parameter used for non-isotropic horizontal correlations
                      ! at model levels & obs. location
    znismq   (:,:) ,& ! parameter used for non-isotropic horizontal correlations
                      ! at vertically equidistant points & obs. location
! local information on upper-air single-level reports
    xoiua  (:,:,:) ,& ! observation increments
    zwtua  (:,:,:) ,& ! temporal weights
    zsprob     (:) ,& ! spreading parameter at the obs. point
    fcorlua(:,:,:) ,& ! reduction to vertical correl. below / above the obs.
    zvcutua(:,:,:) ,& ! vertical cut-off below / above the obs.
    zrtdgua    (:) ,& ! convertor height to pressure units for vertic. correlat.
    zpobua     (:) ,& ! pressure at the obs. point
    zzobua     (:) ,& ! height at the obs. point
    zriflua(:,:,:) ,& ! upper estimate to horizontal radius of influence
    zvidua     (:) ,& ! vertical interpolation distance of adjacent model levels
    zqualua(:,:,:) ,& ! quality weight factor
    zsprtp     (:) ,& ! spreading parameter at the tropopause & obs. location
    zsprua   (:,:) ,& ! spreading parameter at model levels & obs. location
    zlopua   (:,:) ,& ! log( pressure ) at model levels & obs. location
    znisua   (:,:) ,& ! parameter used for non-isotropic horizontal correlations
                      ! at model levels & obs. location
    znisuq   (:,:)    ! parameter used for non-isotropic horizontal correlations
                      ! at vertically equidistant points & obs. location

  REAL    (KIND = wp)       , ALLOCATABLE :: &
! local information on surface-level reports
    xoisu  (:,:,:) ,& ! observation increments
    zwtsu  (:,:,:) ,& ! temporal weights
    zsposu     (:) ,& ! spreading parameter at the obs. point
    zvcutsu  (:,:) ,& ! vertical cut-off above the obs.
    zrtdgsu    (:) ,& ! convertor height to pressure units for vertic. correlat.
    zpobsu     (:) ,& ! pressure at the obs. point
    zzobsu     (:) ,& ! height at the obs. point
    zriflsu  (:,:) ,& ! upper estimate to horizontal radius of influence
    zqualsu(:,:,:) ,& ! quality weight factor
    zsprsu   (:,:) ,& ! spreading parameter at model levels & obs. location
    zpblsu   (:,:) ,& ! potential temperature difference rel. to the obs. level
    znissu   (:,:) ,& ! parameter used for non-isotropic horizontal correlations
                      ! at model levels & obs. location
    znissq   (:,:) ,& ! parameter used for non-isotropic horizontal correlations
                      ! at vertically equidistant points & obs. location
! local information on 'surface' (i.e. on lowest model level) pressure data
    zoips    (:,:) ,& ! observation increments
    omykps   (:,:) ,& ! local part of nudging weight, i.e. exclud. horiz. weight
    r1ifps     (:) ,& ! horizontal correlation scale
    zmassps    (:) ,& ! total mass affected by the 'temperature correction'
    ztdpps     (:) ,& ! 'temperature / pressure' at obs. point
    qcimps   (:,:) ,& ! model pressure interpolated to obs. point
    qctmps   (:,:) ,& ! time of all observations
    qcqfps   (:,:) ,& ! weight factor related to obs. quality (representiveness)
    zoips_b  (:,:) ,& ! increments: observations minus lateral boundary fields
! local information on integrated water vapour increments
    zoiciv   (:,:) ,& ! observation increments
    zqcfiv   (:,:) ,& ! local part of nudging weight, i.e. exclud. horiz. weight
    zmodiv     (:) ,& ! model integrated water vapour IWV
    zsativ     (:) ,& ! IWV of saturated model temperature profile
    zactiv   (:,:)    ! observation time


  INTEGER (KIND = iintegers), ALLOCATABLE :: &
! local information on multi-level reports
    ioml       (:) ,& ! longitudinal index of station location on local domain
    joml       (:) ,& ! latitudinal  index of station location on local domain
    ioml_tot   (:) ,& ! longitudinal index of station location on global domain
    joml_tot   (:) ,& ! latitudinal  index of station location on global domain
    kviflml(:,:,:) ,& ! vertical range of possibly influenced model levels
    kkoiml   (:,:) ,& ! level index in the OIR (observation increment record)
                      ! for the lowest levels of the reports
    ksprml     (:) ,& ! lowest level with spreading along model levels
    kobtyml    (:) ,& ! CMA observation type
    kcdtyml    (:) ,& ! CMA observation code type
    mszlev   (:,:) ,& ! number of vertical levels with obs. incr. in the OIR
    mbotlv (:,:,:) ,& ! number of model levels above the surface without
                      ! observation increment
    mtoplv (:,:,:) ,& ! number of model levels at the top of the model without
                      ! observation increment
! local information on upper-air single-level reports
    ioua       (:) ,& ! longitudinal index of station location on local domain
    joua       (:) ,& ! latitudinal  index of station location on local domain
    ioua_tot   (:) ,& ! longitudinal index of station location on global domain
    joua_tot   (:) ,& ! latitudinal  index of station location on global domain
    kviflua(:,:,:) ,& ! vertical range of possibly influenced model levels
    kobtyua    (:) ,& ! CMA observation type
    kcdtyua    (:) ,& ! CMA observation code type
! local information on surface-level reports
    iosu       (:) ,& ! longitudinal index of station location on local domain
    josu       (:) ,& ! latitudinal  index of station location on local domain
    iosu_tot   (:) ,& ! longitudinal index of station location on global domain
    josu_tot   (:) ,& ! latitudinal  index of station location on global domain
    kviflsu  (:,:) ,& ! vertical range of possibly influenced model levels
    kobtysu    (:) ,& ! CMA observation type
    kcdtysu    (:) ,& ! CMA observation code type
! local information on 'surface' (i.e. on lowest model level) pressure data
    iops       (:) ,& ! longitudinal index of station location on local domain
    jops       (:) ,& ! latitudinal  index of station location on local domain
    iops_tot   (:) ,& ! longitudinal index of station location on global domain
    jops_tot   (:) ,& ! latitudinal  index of station location on global domain
    iqclps     (:) ,& ! indicator for current obs. to undergo spatial check now
    iqcfps     (:) ,& ! indicator for obs. to have passed latest threshold QC
    iqcnps     (:)    ! index of local administrator of (local) ODR

  INTEGER (KIND = iintegers), ALLOCATABLE :: &
! local information on integrated water vapour increments
    ioiv       (:) ,& ! longitudinal index of station location on local domain
    joiv       (:) ,& ! latitudinal  index of station location on local domain
    ioiv_tot   (:) ,& ! longitudinal index of station location on global domain
    joiv_tot   (:) ,& ! latitudinal  index of station location on global domain
    iqcliv     (:) ,& ! indicator for current obs. to undergo spatial check now
    iqcfiv     (:) ,& ! indicator for obs. to have passed latest threshold QC
    kiobiv     (:) ,& ! index of local administrator of (local) ODR
    kioiiv     (:) ,& ! index of local information array for multi-level data
    ktypiv     (:)    ! observation type


  LOGICAL                   , ALLOCATABLE :: &
! local information on multi-level reports
    ltiml    (:,:) ,& ! .TRUE if temporal linear interpolation is applied
! local information on upper-air single-level reports
    ltiua    (:,:) ,& ! .TRUE if temporal linear interpolation is applied
! local information on surface-level reports
    ltisu    (:,:) ,& ! .TRUE if temporal linear interpolation is applied
! local information on 'surface' (i.e. on lowest model level) pressure data
    ltips      (:) ,& ! .TRUE if temporal linear interpolation is applied
    lmlps      (:)    ! .TRUE if datum is derived from multi-level report

  CHARACTER (LEN=ilstidg)   , ALLOCATABLE :: &
    ystidml    (:) ,& ! station identity of multi-level station
    ystidua    (:) ,& ! station identity of upper-air single-level station
    ystidsu    (:) ,& ! station identity of surface-level station
    ystidps    (:) ,& ! station identity of 'surface' pressure station
    ystidiv    (:)    ! station identity of integrated water vapour station


! 3. Observation increment record for multi-level reports 'oiml'
! --------------------------------------------------------------

! Observation increment record for multi-level reports 'oiml'
  INTEGER (KIND = iintegers) , PARAMETER  :: &
    maxnoi = 15  ,& ! length of observation increment record
    noiu   =  1  ,& ! zonal wind observation increment
    noiv   =  2  ,& ! meridional wind observation increment
    noit   =  3  ,& ! temperature observation increment
    noiqd  =  4  ,& ! specific humidity observation increment
    noirh  =  5  ,& ! relative humidity observation increment
    noiz   =  6  ,& ! height (of obs. increment level)
    noith  =  7  ,& ! potential temperature
    noilp  =  8  ,& ! log( pressure )
    noiuqc =  9  ,& ! quality weight to wind obs. increment
    noitqc = 10  ,& ! quality weight to temperature obs. increment
    noiqqc = 11  ,& ! quality weight to humidity obs. increment
    noivcr = 12  ,& ! normalization factor to vertical interpolation distances
    noiulr = 13  ,& ! vertical interpolation   /  for wind
    noitlr = 14  ,& ! distance to the         <   for temperature
    noiqlr = 15     ! obs. increment level     \  for humidity


! 4. (Local) information gathered by 1 (or 2) nodes for printing for control
! --------------------------------------------------------------------------

  REAL    (KIND = wp)       , ALLOCATABLE :: &
    oyqc     (:,:)    ! on data rejected by the threshold quality control

  INTEGER (KIND = iintegers), ALLOCATABLE :: &
    myqc     (:,:)    ! on data rejected by the threshold quality control

  CHARACTER (LEN=ilstidp)   , ALLOCATABLE :: &
    yyqc       (:)    ! on data rejected by the threshold quality control

  REAL    (KIND = wp)       , ALLOCATABLE :: &
    oysu     (:,:)    ! on interpol. of surface-level data to lowest model level

  CHARACTER (LEN=ilstidp)   , ALLOCATABLE :: &
    yysu       (:)    ! on interpol. of surface-level data to lowest model level


! 5. Variables defining the size of the arrays containing the local information
! -----------------------------------------------------------------------------

  INTEGER (KIND = iintegers)          :: &
    maxoil         ,& ! max.  number of vertical levels in the OIR 'oiml'
    maxmloi_tot    ,& ! length of arrays with local info on multi-level reports
                      !   (except 'oiml'):
                      !   maxmloi_tot >= maxmlo + maxgpo/2 + maxtvo
    maxpso         ,& ! length of arrays with local info on sfc. press. reports
    maxivq         ,& ! length of arrays for IWV increm. used for spatial check
    maxqcp         ,& ! length of arrays for QC rejected data to be printed
    maxysu         ,& ! length of arrays for extrapolated printed surface reps.
    nmloit         ,& ! total number of active profiles of obs. incr. in the OIR
    nmltot         ,& ! total number of active multi-level stations
    nuatot         ,& ! total number of active upper-air single-level stations
    nsutot         ,& ! total number of active surface-level stations
    npstot         ,& ! total number of active surface-pressure stations
    nivtot         ,& ! total number of IWV reports used for spatial check
    ntotqc         ,& ! total number of rejected data to be printed per timestep
    ntotys         ,& ! total number of printed interpol. surface-level reports
    nuaex          ,& ! number of local upper-air obs not used due to ODR size
    ktopsu            ! uppermost model level possibly influenced by surface obs

  LOGICAL                           :: &
    lnissu         ,& ! non-isotrophic correlations for surface-level data
    lnisua            ! non-isotrophic correlat. for upper-air single-lev. data

! 6. Geometrics and variables used for spatial consistency check of pressure obs
! ------------------------------------------------------------------------------

  REAL    (KIND = wp)       , ALLOCATABLE :: &
    pyyd       (:) ,& ! latitudinal (merid.) distance on tangent cone projection
    pxxd2    (:,:) ,& ! square of zonal distance on tangent cone projection
    pxsalpa  (:,:)    ! factor used for distances on tangent cone projection

  INTEGER (KIND = iintegers), ALLOCATABLE :: &
    isrtpqc    (:) ,& ! (sorted) list of stations with 'surface' pressure data
    isrtvqc    (:)    ! (sorted) list of stations with IWV

!-------------------------------------------------------------------------------

END MODULE data_nudge_gather
