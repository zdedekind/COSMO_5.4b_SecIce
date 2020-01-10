!+ Data module for variables only used for the spreading and nudging equations
!-------------------------------------------------------------------------------

MODULE data_nudge_spread

!-------------------------------------------------------------------------------
!
! Description:
!   This module contains those variables (scalars and arrays) which are merely
!   concerned with the spreading of observation increments, computation of the
!   analysis increments, and updating of the prognostic variables by nudging,
!   i.e. after the processor has received all the required 'local' information
!   from all observations. These are
!    - the analysis increment fields
!    - the output of spreading procedures (sum of weights and weighted increm.)
!    - further fields used to compute analysis increments
!    - geometrical fields used for horizontal distances and wind correlations
!    - the horizontal fields used for the spreading of the obs. increments
!    - other fields, variables, and indices used for the spreading
!
!   The arrays are declared as allocatable arrays and are allocated after the
!   gathering of the 'local information' and deallocated after updating the
!   prognostic fields due to the nudging, except for some analysis increment
!   fields (section 1), which are allocated in the setup and deallocated in
!   the cleanup of the model.
!   The scalar variables are initialized after the gathering of the 'local
!   information' or during the spreading.
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
! 1.15       1998/11/02 Christoph Schraff
!  All global allocatable arrays moved from module 'src_mult_spread'.
! 2.4        2001/01/29 Christoph Schraff
!  'tairh' and 'taips' added, to differentiate humidity balance at nudging T.
!  Variables for '?cutof' compression added, for efficiency on vector processor.
! 3.18       2006/03/03 Christoph Schraff
!  Section 9 on observation types in weighted increment arrays 'zwi' introduced.
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_22        2012/01/31 Christoph Schraff
!  Fields 'psaigeo' and 'ztpgeo' added, for geostrophic pressure increments.
!  'zablpo' introduced, 'zablfc' re-defined for speed-up. 'k', 'lconai' removed.
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

! 1. Analysis increment fields (to be kept in the long-term storage)
! ------------------------------------------------------------------

  REAL    (KIND = wp)       , ALLOCATABLE :: &
    uanai  (:,:,:) ,& ! analysis increments of zonal wind
    vanai  (:,:,:) ,& ! analysis increments of meridional wind
    tanai  (:,:,:) ,& ! analysis increments of temperature
    qanai  (:,:,:) ,& ! analysis increments of specific humidity

    psanai   (:,:) ,& ! analysis increments of pressure at lowest model level
    psaigeo  (:,:) ,& ! analysis increments of pressure at lowest model level
                      !   from geostrophic balancing of wind analysis increments
    taips  (:,:,:)    ! analysis increment part from pressure nudging

! 2. Analysis increment fields
! ----------------------------

  REAL    (KIND = wp)       , ALLOCATABLE :: &
    zpai   (:,:,:) ,& ! analysis incr. of pressure, condens./evapor. not incl.
    zroi   (:,:,:) ,& ! analysis incr. of density , condens./evapor. not incl.
    ztwips (:,:,:) ,& ! temperature increments implied (physic'ly or statistic.)
                      ! by pressure analysis increments at the lowest model lev.
    ztpgeo (:,:,:) ,& ! temperature increments implied (physic'ly or statistic.)
                      ! by geostrophic (surface) pressure analysis increments

! 3. Output of spreading procedures
! ---------------------------------

    omy  (:,:,:,:) ,& ! sum of  spatial * temporal * quality 'spreading weights'
    om2  (:,:,:,:) ,& ! sum of  squares of 'spreading weights'
    zwi  (:,:,:,:) ,& ! sum of weighted (observation) increments (weights are
                      !        squares of 'spreading weights')

! 4. Further fields used to compute analysis increments
! -----------------------------------------------------

    faclbc (:,:,:) ,& ! reduction of nudging weights (near lateral boundaries)

! 5. Geometrical fields used for horizontal distances and wind correlations
! -------------------------------------------------------------------------

    yyd        (:) ,& ! latitudinal (merid.) distance on tangent cone projection
    yyd2       (:) ,& ! yyd **2
    xxd2     (:,:) ,& ! square of zonal distance on tangent cone projection
    c2alpa   (:,:) ,& !  / further factors used to compute
    scalpa   (:,:) ,& ! /  the 2-dim. horizontal wind correlations
    xcalpa   (:,:) ,& ! \  and the horizontal distances
    xsalpa   (:,:) ,& !  \ on the tangent cone projection

! 6. Further horizontal input fields for the spreading of obs. increments
! -----------------------------------------------------------------------

    zlop     (:,:) ,& ! log( pressure )
    zeklop   (:,:) ,& ! exp( R/cp * log(p) )
    zthvg    (:,:) ,& ! vertical gradient of potential temperature
    zpk      (:,:) ,& ! pressure
!   zqdmod   (:,:) ,& ! specific water vapour content
!   zewmod   (:,:) ,& ! vapour pressure
!   zewsat   (:,:) ,& ! saturation vapour pressure
    zspr   (:,:,:) ,& ! spreading parameter , param. def. non-isotropic weights
    zdds   (:,:,:) ,& ! scaled horizontal distance betw. obs. and target grid pt
    zcoruu   (:,:) ,& ! zonal  wind - zonal  wind correlation  \  (without
    zcoruv   (:,:) ,& ! zonal  wind - merid. wind correlation   \  EXP( -zdds )
    zcorvu   (:,:) ,& ! merid. wind - zonal  wind correlation   /  -term )
    zcorvv   (:,:) ,& ! merid. wind - merid. wind correlation  /
    rhscal   (:,:) ,& ! data density-dependent factor to ps horiz. correl. scale
    zablpo   (:,:) ,& ! =1 if grid pt. is (fully) above the ABL, =0 otherwise
    zablfc (:,:,:,:)  ! reduced weighting inside/above ABL for upper-air/sfc obs

! 7. Further fields used for or during the spreading
! --------------------------------------------------

  INTEGER (KIND = iintegers), ALLOCATABLE :: &
    klva   (:,:,:) ,& ! indices of upper  __\  obs. increment levels used for
    klvb   (:,:,:)    ! indices of lower    /  vertical interpolation to grid pt

  REAL    (KIND = wp)       , ALLOCATABLE :: &
    zoic1    (:,:) ,& ! spreaded obs. incr. from first multi-level report
    zoiv1    (:,:) ,& ! 2nd component of spreaded obs. incr. from first report
    omys1    (:,:)    ! weight from first report

  REAL    (KIND = wp)       , ALLOCATABLE :: &
    zsdnis   (:,:) ,& ! equidistant pts. used for non-isotropic horiz. correlat.
    zsdnid   (:,:) ,& ! distance between two adjacent 'equidistant' pts.
    gppkmi     (:)    ! convertor for zonal dist. from 'km' to rotated grid pts.

  LOGICAL                   , ALLOCATABLE :: &
    lcutof (:,:,:)    ! .TRUE if grid pt. is within area of influence of report

!   variables for cutof compression
  INTEGER (KIND=iintegers)  , ALLOCATABLE :: &
    icutof   (:,:) ,& ! x-coordinate of grid pt. within area of influence
    jcutof   (:,:) ,& ! x-coordinate of grid pt. within area of influence
    icutmp     (:) ,& ! x-coordinate of grid pt. within area of influence, tmp
    jcutmp     (:)    ! x-coordinate of grid pt. within area of influence, tmp

  INTEGER (KIND=iintegers) :: &
    ncutof     (2) ,& ! number of grid points within area of influence
    idimcut           ! dimension of 'icutof', 'jcutof', etc.

! 8. Indices and lists of indices
! -------------------------------

  INTEGER (KIND=iintegers)  , ALLOCATABLE :: &
    isortps    (:) ,& ! (sorted) list of stations with 'surface' pressure data
    isortsu    (:) ,& ! (sorted) list of stations with surface-level data
    isortua    (:) ,& ! (sorted) list of stations with upper-air single-lv. data
    isortml    (:)    ! (sorted) list of stations with multi-level data

  INTEGER (KIND=iintegers)          :: &
    io     ,jo     ,& ! local  indices of location of observation
    io_tot ,jo_tot ,& ! global indices of location of observation
    istaspr,jstaspr,& ! lower left corner of domain containing area of influence
    iendspr,jendspr,& ! upper right corner of domain containing area of infl.
    jrs_tot,jre_tot,& ! index range for convertor for zonal distances 'gppkmi'
    ista           ,& ! index of observing station
    itim           ,& ! (time) index over obs. at one observing station
    kml250            ! index of full model level corresponding to about 250 hPa

! 9. Observation types in weighted increment arrays 'zwi'
! -------------------------------------------------------

  INTEGER (KIND=iintegers)  , PARAMETER  :: &
    mxotyp   =  5  ,& ! number of observation types
    noiras   =  1  ,& ! radiosonde
    noisfc   =  2  ,& ! surface-level (SYNOP, SHIP, BUOY)
    noiair   =  3  ,& ! aircraft
    noisat   =  4  ,& ! satellite (SATOB, ATOVS, MSG)
    noigps   =  5     ! GPS

! 10. Other variables
! -------------------

  REAL    (KIND = wp)       , ALLOCATABLE :: &
    wvgeo      (:)    ! vertical weights to the geostrophic wind correction

  REAL    (KIND = wp)        :: &
    rtinfl         ,& ! horizontal correlation scale
    gppkmj            ! convertor for merid. dist. from 'km' to rotated grid pts

  INTEGER (KIND=iintegers)   :: &
    nhisto (15)       ! data density histogram for surface pressure obs
                      !  (used for output only)


!-------------------------------------------------------------------------------

END MODULE data_nudge_spread
