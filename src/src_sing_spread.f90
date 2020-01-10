!+ Source module for spreading single-level obs. increments for the nudging
!-------------------------------------------------------------------------------

MODULE src_sing_spread

!-------------------------------------------------------------------------------
!
! Description:
!   The module "src_sing_spread" performs the spreading of observational
!   information from all single-level data (upper-air single-level data,
!   surface-level (i.e. screen-level) data, and surface pressure data) from
!   the total model domain to the model grid points on the local sub-domain.
!
! Method:
!   This module contains the following procedures:
!    - ps_spreading      : spreading of the 'surface' pressure obs increments
!    - surf_upair_spread : spreading of single-level observation increments
!    - surf_org_spread   : organize spreading of screen-level obs increments
!    - upair_org_spread  : organize spreading of upper-air single-level obs inc.
!   Driving routine is "organize_nudging" of module "src_obs_use_org.f90".
!   'surf_upair_spread' is called in 'surf_org_spread' and 'upair_org_spread'.
!
!   This module also contains an elemental function, formerly statement funct.:
!    - rmod              : MOD function for positive reals
!
!   Note: This module does not contain any communication between PE's.
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.13       1998/10/22 Christoph Schraff
!  Initial release
! 1.19       1998/12/11 Christoph Schraff
!  ANSI violations removed.
! 1.20       1999/01/07 Guenhter Doms
!  Renaming of some global variables
! 1.31       1999/07/01 Christoph Schraff
!  Quantities related to MPI communicator 'icomm_world' replaced.
! 2.4        2001/01/29 Christoph Schraff
!  Reorganisation of loops to enhance efficiency on vector processors and to
!  render spreading loops more uniform.
! 2.5        2001/06/01 Christoph Schraff
!  Savety test at array allocations.
! 2.19       2002/10/24 Christoph Schraff
!  Bug correction at lateral spreading of surface pressure data.
!  Optional stability-dependent vertical weights for surface-level data.
! 3.6        2003/12/11 Christoph Schraff / NEC (/Ulrich Schaettler)
!  Introduction of SATOB reports.
!  Splitting of a Loop for better vectorization.
! 3.12       2004/09/15 Christoph Schraff
!  Minor bug correction related to the spreading in the lowest model layer.
! 3.18       2006/03/03 Christoph Schraff
!  New option for alternative weighting for multiple observations.
!  Option for separate weighting for different observation types.
!  Option to render horizontal correlation scale for surface pressure
!  dependent on data density. Flushing of output files.
! V3_23        2007/03/30 Ulrich Schaettler
!  Introduced ldump_ascii for flushing the ASCII files
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_22        2012/01/31 Christoph Schraff
!  - Refined option for computing net increments and weights separately for
!    different pre-specified sets of observing systems.
!  - Optimisations for decreased cost on NEC-SX9 (e.g to allow for vectorisation
!    in section 1 of subroutines ps_spreading, surf_org_spread, upair_org_spread
!    or in section 3 of upair_org_spread for horizontal correlation scales).
!  - Some variables moved from 'data_nudge_all' to 'data_obs_lib_cosmo'.
!  - Subroutine arguments lists extended by vertical model level 'k'.
! V4_28        2013/07/12 Christoph Schraff
!  Statement function 'rmod' replaced by elemental function. Improved comments.
! V5_1         2014-11-28 Christoph Schraff, Oliver Fuhrer
!  Processing of Mode-S aircraft obs introduced (CS)
!  Replaced ireals by wp (working precision) (OF)
! V5_2         2015-05-21 Ulrich Schaettler
!  Initialized lwobs in SR upair_org_spread, surf_org_spread to .FALSE.
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!===============================================================================
!
! Declarations:
!
! Modules used:
!
!-------------------------------------------------------------------------------

USE data_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!-------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------

    ie_tot,       & ! number of grid points in zonal direction
    je_tot,       & ! number of grid points in meridional direction
    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction

! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the 
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from 
!    the other ones because of the use of the staggered Arakawa-B-grid.
!    
!   zonal direction
    istart,       & ! start index for the forecast of w, t, qd, qw and pp
    iend,         & ! end index for the forecast of w, t, qd, qw and pp

!   meridional direction
    jstart,       & ! start index for the forecast of w, t, qd, qw and pp
    jend,         & ! end index for the forecast of w, t, qd, qw and pp

! 5. variables for the time discretization and related variables
! --------------------------------------------------------------

    dt,           & ! long time-step
    dt2,          & ! 2 * dt
    dtdeh           ! dt / 3600 seconds

! end of data_modelconfig

!-------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

    rdocp           ! r_d / cp_d

! end of data_constants

!-------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------

    p0         ,    & ! reference pressure at main levels             ( Pa)
    hhl        ,    & ! geometical height of half model levels        ( m )

! 3. prognostic variables                                             (unit)
! -----------------------

    t          ,    & ! temperature                                   (  k  )
    pp                ! deviation from the reference pressure         ( pa  )

! end of data_fields

!-------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------

    nstop,        & ! last time step of the forecast
    ntstep,       & ! actual time step
                    ! index for permutation of three time levels
    nnew,         & ! corresponds to ntstep + 1

! 7. additional control variables
! -------------------------------
    ldump_ascii,  & ! for flushing (close and re-open) the ASCII files
    l2tls           ! forecast with 2-TL integration scheme

! end of data_runcontrol 

!-------------------------------------------------------------------------------

USE data_nudge_all , ONLY :   &

! 1. parameters and related variables
! -----------------------------------

    lwonl        ,& ! .TRUE if grid pt. (ionl ,jonl ) lies in the local domain
    lwonl2       ,& ! .TRUE if grid pt. (ionl2,jonl2) lies in the local domain
    lcroot       ,& ! .TRUE if (my_cart_id == 0) (for print on std. output)
    lfirst       ,& ! .TRUE if 'organize_nudging' is called the first time
    acthr        ,& ! actual forecast hour

! 2. namelist variables controlling the data assimilation
! -------------------------------------------------------

    nwtyp        ,& ! 1   : if > 1 then compute net obs. increments for 'nwtyp'
                    !       different sets of observing systems separately
    niwtyp       ,& ! 1,0,0,..: number of obs or code types which belong to
                    !           the sets of observing systems
    iwtyp        ,& ! 0,0,0,..: obs or code types belonging to a set of obs
                    !           system, specified successively for each set
    kwtyp        ,& ! 1   : function for weights W for multiple observations
    nudgend      ,& ! 0   : end of nudging period in timesteps
    tconbox      ,& ! 6*dt: timestep [s] for computing analysis increments
                    !       (i.e. time box of constant analysis increments)
    luvgcor      ,& ! .t. : .t. ==> geostrophic wind correction applied
    mpsgcor      ,& ! 1   : mode to apply geostrophic pressure correction
    ntpscor      ,& ! 1   : switch for temperature correction when the 'surface'
                    !       pressure (i.e. p. at lowest model level) is nudged
    ptpstop      ,& ! 400.: pressure [hPa] at upper boundary of the temperature
                    !       correction for 'surface' pressure nudging
    gnudg        ,& ! 6,12,6,6*10^-4: nudging coefficients for TEMP / PILOT data
    gnudgar      ,& ! 6, 0,6,0*10^-4: nudging coeffic. for AIRCRAFT data [1/s]
    gnudgms      ,& ! 6, 0,6,0*10^-4: nudging coeffic. for Mode-S aircraft [1/s]
    gnudgsu      ,& ! 6,12,0,6*10^-4: nudging coef. for surface-level data [1/s]
    wtukrsa      ,& ! 3.0 : temporal radius of influence towards the past
                    !       relative to the obs. time .   for TEMP / PILOT
    wtukara      ,& ! 1.5 : temporal radius of influence towards the past
                    !       relative to the obs. time .   for AIRCRAFT data
    wtuksua         ! 1.5 : temporal radius of influence towards the past
                    !       relative to the obs. time.   for surface-level data

USE data_nudge_all , ONLY :   &

    msprpar      ,& ! 1   : switch specifying the surfaces along which observat.
                    !       increments of upper-air data are (primarily) spread 
    msprpsu      ,& ! 0   : switch specifying the surface along which surface-  
                    !       level data increments are primarily spreaded
    vcorls       ,& ! 2*.333,: square of the vertical correlation scale,
                    ! 2*.04    i.e. of the Gaussian vertical influence 'radius'
                    !          in log( pressure ) if (msprpar <= 1), or
                    !          in potential temperature if (msprpar == 2)
    wablua       ,& ! 4*  1. : factor to weights within the ABL
                    !          (atmospheric boundary layer, i.e. mixed layer)
    vcorlsu      ,& ! 2*.013,: square of the vertical correlation scale,
                    ! 2*.002   i.e. of the Gaussian vertical influence 'radius'
                    !          in log( pressure ) if (msprpsu == 1) or
                    !          in potential temperature if (msprpsu == 2)
    vpblsu       ,& ! 2*99. ,: Gaussian vertical influence 'radius' of potential
                    ! 2*99.    temperature differences between obs. level and
                    !          model level, to render vertical weights depend on
                    !          the stability (in the PBL) even if (msprpsu <= 1)
    wablsu       ,& ! 4*  1. : factor to weights above the ABL
                    !          (atmospheric boundary layer, i.e. mixed layer)
    rhinfl       ,& ! 0.,70.,: constant part of the 'correlation scale of the
                    ! 0.,0.    autoregressive horizontal correlation function'
                    !          (='COSAC') [km]
    rhvfac       ,& ! 1., 0.,: multiplication factor to the vertically varying
                    ! 2* .83   part of the 'COSAC' (as def. in 'data_nudge_all')
    rhtfac       ,& ! 1.3,   : scaling factor of the total 'COSAC' for the
                    ! 1.43,    beginning and end of the nudging period for 1 obs
                    ! 1.3,     relative to the 'COSAC' at the obs. time as given
                    ! 1.3      by 'rhinfl', 'rhvfac')
    cutofr       ,& ! 4* 3.5 : cut-off in 'COSAC' units of the horizontal
                    !          correlation function
    vcsni        ,& ! 4*2500.: square of Gaussian vertical influence 'radius'
                    !          in potential temperature (if msprpar <= 1) or
                    !          log( pressure ) (if msprpar == 2) on surfaces
                    !          along which obs. increments are spread laterally
    rhiflsu      ,& ! 2 *70.,: constant part of the 'correlation scale of the
                    !   100.,  autoregressive horiz. correlation function'
                    !    70.   (='COSAC') [km]  (at the obs. time)
    rhtfsu       ,& ! 1.,    : scaling factor of the total 'COSAC' for the
                    ! 1.43,    beginning and end of the nudging period for 1 obs
                    ! 1.,      relative to the 'COSAC' at the obs. time as given
                    ! 1.,      by 'rhiflsu')
    cutofsu      ,& ! 2., 3.5: cut-off in 'COSAC' units of the horizontal
                    ! 2., 2.   correlation function
    vcsnisu      ,& ! 2*2500.: square of Gaussian vertical influence 'radius'
                    ! 2*   9.  in potential temperature (if msprpsu <= 1) or
                    !          log( pressure ) (if msprpsu == 2) on surfaces
                    !          along which obs. increments are spread laterally 
    cnondiv      ,& ! .1  : constant part of the factor to the non-divergence
                    !       correction ('fanodicor') in 2-dim wind correlations
    fnondiv      ,& ! .8  : multiplication factor to the vertically varying part
                    !       of the factor to the non-divergence correction
                    !       (as defined in 'data_nudge_all')
    tnondiv      ,& ! 1.1 : temporal multiplication factor to the factor to the
                    !       non-divergence correction for the beginning and end
                    !       of the nudging period for 1 obs relative to that
                    !       given by 'cnondiv', 'fnondiv' for the obs. time
    rhfpsdd      ,& ! 1.0    : minimum scaling (reduction) factor of the total
                    !          'COSAC' for surface pressure dep. on data density
    ionl         ,& ! 167    : / grid point coordinates
    jonl         ,& ! 103    : \ for standard output on nudging
    ionl2        ,& ! 167    : / 2nd grid pt coordinates
    jonl2        ,& ! 103    : \ for other standard output on nudging

! 5. Miscellany
! -------------

    isetyp0      ,& ! index in 'niwtyp' which points to the first index with
                    ! observation type '0' in 'iwtyp' (denoting the set of obs
                    ! systems with contains all ('remaining') obs types that are
                    ! not specified explicitly in 'iwtyp')
    isetyp          ! defines for each observation or code type in 'iwtyp' which
                    ! set of observing systems (index of 'niwtyp') it belongs to

! end of data_nudge_all

!-------------------------------------------------------------------------------

USE data_obs_lib_cosmo , ONLY :   &

! 1. General parameters
! ---------------------

    c0         ,& ! standard real constant 0.0
    c1         ,& ! standard real constant 1.0
    c2         ,& ! standard real constant 2.0
    c05        ,& ! standard real constant 0.5
    c3600      ,& ! standard real constant 3600.0
    epsy       ,& ! = 1.E-8_wp : commonly used very small value > 0
    i1         ,& ! standard integer constant 1

! 3. pressure dependent scales and geometry of horizontal correlations
! --------------------------------------------------------------------

    ncolev     ,& ! number of levels in the correlation scale tables
    tabcolp    ,& ! LOG( tabcop )
    rhvsond    ,& ! upper-air wind horizontal correlation scales
                  ! (pressure dependent part)
    rhtsond    ,& ! upper-air temperature horiz. correlation scales
    rhqsond    ,& ! upper-air humidity horiz. correlation scales
    fnodivc    ,& ! factor to non-divergence correction in 2-d wind correlation

! 4. I/O device numbers for obs processing / nudging and file names
! -----------------------------------------------------------------

    yuprint    ,& ! file name for all the remaining information
    nupr       ,& ! unit number of file for all the remaining information

! 5. CMA observation type and code type numbers
! ---------------------------------------------

    nsynop     ,& ! SYNOP reports
    nairep     ,& ! AIREP reports (all aircraft reports)
    nsatob     ,& ! SATOB reports
    ndribu     ,& ! DRIBU reports
    ntemp      ,& ! TEMP  reports
    npilot     ,& ! PILOT reports
    nscatt     ,& ! SCATT reports
    nmodes        !   mode-s report

! end of data_obs_lib_cosmo

!-------------------------------------------------------------------------------

USE data_nudge_gather , ONLY :   &

! 1. Parameters and general variables
! -----------------------------------

    ilstidp      ,& ! char. length used for printing the station ID
    ystid        ,& ! obs. station identity to be printed
    qst          ,& ! quotient of the mean vertical potential temperature grad.
                    ! in the stratosphere to the analog. tropospheric gradient
    mxispr       ,& !     number of equidistant vertical points
    xezispr      ,& !   / parameters used to
    dezispr      ,& !  /  define the vertically
    sezispr      ,& ! <   equidistant points
    xthispr      ,& !  \  (for non-isotropic
    dthispr      ,& !   \ horizontal correlations)
    vcutnit      ,& ! horiz. correlations are non-isotropic if the
    vcutnip      ,& ! \  vertical scales  'rdsprni' < 'zcutnit,p'
    ktp          ,& ! lowermost purely horizontal model main level
    ktth            ! top model level with spreading along isentropic surfaces

USE data_nudge_gather , ONLY :   &

! 2. The required local information on the observations and their location
! ------------------------------------------------------------------------

! local information on upper-air single-level reports
    xoiua        ,& ! observation increments
    zwtua        ,& ! temporal weights
    zsprob       ,& ! spreading parameter at the obs. point
    fcorlua      ,& ! reduction to vertical correl. below / above the obs.
    zvcutua      ,& ! vertical cut-off below / above the obs.
    zrtdgua      ,& ! convertor height to pressure units for vertic. correlat.
    zpobua       ,& ! pressure at the obs. point
    zzobua       ,& ! height at the obs. point
    zriflua      ,& ! upper estimate to horizontal radius of influence
    zvidua       ,& ! vertical interpolation distance of adjacent model levels
    zqualua      ,& ! quality weight factor
    zsprtp       ,& ! spreading parameter at the tropopause & obs. location
    zsprua       ,& ! spreading parameter at model levels & obs. location
    zlopua       ,& ! log( pressure ) at model levels & obs. location
    znisua       ,& ! parameter used for non-isotropic horizontal correlations
                    ! at model levels & obs. location
    znisuq       ,& ! parameter used for non-isotropic horizontal correlations
                    ! at vertically equidistant points & obs. location
    ioua         ,& ! longitudinal index of station location on local domain
    joua         ,& ! latitudinal  index of station location on local domain
    ioua_tot     ,& ! longitudinal index of station location on global domain
    joua_tot     ,& ! latitudinal  index of station location on global domain
    kviflua      ,& ! vertical range of possibly influenced model levels
    kobtyua      ,& ! CMA observation type
    kcdtyua      ,& ! CMA observation code type
    ltiua        ,& ! .TRUE if temporal linear interpol. for upper-air sing.l.
    ystidua         ! station identity of upper-air single-level station

USE data_nudge_gather , ONLY :   &

! local information on surface-level reports
    xoisu        ,& ! observation increments
    zwtsu        ,& ! temporal weights
    zsposu       ,& ! spreading parameter at the obs. point
    zvcutsu      ,& ! vertical cut-off above the obs.
    zrtdgsu      ,& ! convertor height to pressure units for vertic. correlat.
    zpobsu       ,& ! pressure at the obs. point
    zzobsu       ,& ! height at the obs. point
    zriflsu      ,& ! upper estimate to horizontal radius of influence
    zqualsu      ,& ! quality weight factor
    zsprsu       ,& ! spreading parameter at model levels & obs. location
    zpblsu       ,& ! potential temperature difference rel. to the obs. level
    znissu       ,& ! parameter used for non-isotropic horizontal correlations
                    ! at model levels & obs. location
    znissq       ,& ! parameter used for non-isotropic horizontal correlations
                    ! at vertically equidistant points & obs. location
    iosu         ,& ! longitudinal index of station location on local domain
    josu         ,& ! latitudinal  index of station location on local domain
    iosu_tot     ,& ! longitudinal index of station location on global domain
    josu_tot     ,& ! latitudinal  index of station location on global domain
    kviflsu      ,& ! vertical range of possibly influenced model levels
    kobtysu      ,& ! CMA observation type
    kcdtysu      ,& ! CMA observation code type
    ltisu        ,& ! .TRUE if temporal linear interpol. for surface-level rep
    ystidsu      ,& ! station identity of surface-level station
! local information on 'surface' (i.e. on lowest model level) pressure data
    zoips        ,& ! observation increments
    omykps       ,& ! local part of nudging weight, i.e. exclud. horiz. weight
    r1ifps       ,& ! horizontal correlation scale
    zmassps      ,& ! total mass affected by the 'temperature correction'
    ztdpps       ,& ! 'temperature / pressure' at obs. point
    iops         ,& ! longitudinal index of station location on local domain
    jops         ,& ! latitudinal  index of station location on local domain
    iops_tot     ,& ! longitudinal index of station location on global domain
    jops_tot     ,& ! latitudinal  index of station location on global domain
    ltips        ,& ! .TRUE if temporal linear interpol. for surface pressure
    lmlps        ,& ! .TRUE if datum is derived from multi-level report
    ystidps      ,& ! station identity of 'surface' pressure station

! 5. Variables defining the size of the arrays containing the local information
! -----------------------------------------------------------------------------

    nuatot       ,& ! total number of active upper-air single-level stations
    nsutot       ,& ! total number of active surface-level stations
    npstot       ,& ! total number of active surface-pressure stations
    lnissu       ,& ! non-isotrophic correlations for surface-level data
    lnisua          ! non-isotrophic correlat. for upper-air single-lev. data

! end of data_nudge_gather

!-------------------------------------------------------------------------------

USE data_nudge_spread , ONLY :   &

! 1. Analysis increment fields (to be kept in the long-term storage)
! ------------------------------------------------------------------

    psanai       ,& ! analysis increments of pressure at lowest model level

! 3. Output of spreading procedures
! ---------------------------------

    omy          ,& ! sum of  spatial * temporal * quality 'spreading weights'
    om2          ,& ! sum of  squares of 'spreading weights'
    zwi          ,& ! sum of weighted (observation) increments (weights are
                    !        squares of 'spreading weights')

! 4. Further fields used to compute analysis increments
! -----------------------------------------------------

    faclbc       ,& ! reduction of nudging weights (near lateral boundaries)

! 5. Geometrical fields used for horizontal distances and wind correlations
! -------------------------------------------------------------------------

    yyd          ,& ! latitudinal (merid.) distance on tangent cone projection
    yyd2         ,& ! yyd **2
    xxd2         ,& ! square of zonal distance on tangent cone projection
    xsalpa       ,& ! factor used for horiz. distances on tangent cone project.

! 6. Further horizontal input fields for the spreading of obs. increments
! -----------------------------------------------------------------------

    zeklop       ,& ! exp( R/cp * log(p) )
    zspr         ,& ! spreading parameter , param. def. non-isotropic weights
    zdds         ,& ! scaled horizontal distance betw. obs. and target grid pt
    zcoruu       ,& ! zonal  wind - zonal  wind correlation  \  (without
    zcoruv       ,& ! zonal  wind - merid. wind correlation   \  EXP( -zdds )
    zcorvu       ,& ! merid. wind - zonal  wind correlation   /  -term )
    zcorvv       ,& ! merid. wind - merid. wind correlation  /
    zablfc       ,& ! reduced weighting inside/above ABL for upper-air/sfc obs
    rhscal       ,& ! data density-dependent factor to ps horiz. correl. scale

! 7. Further fields and variables used for or during the spreading
! ----------------------------------------------------------------

    zsdnis       ,& ! equidistant pts. used for non-isotropic horiz. correlat.
    zsdnid       ,& ! distance between two adjacent 'equidistant' pts.
    gppkmi       ,& ! convertor for zonal dist. from 'km' to rotated grid pts.
    gppkmj       ,& ! convertor for merid. dist. from 'km' to rotated grid pts
    rtinfl       ,& ! horizontal correlation scale
    lcutof       ,& ! .TRUE if grid pt. is within area of influence of report
    icutof       ,& ! x-coordinate of grid pt. within area of influence
    jcutof       ,& ! x-coordinate of grid pt. within area of influence
    ncutof          ! number of grid points within area of influence

USE data_nudge_spread , ONLY :   &

! 8. Indices and lists of indices
! -------------------------------

    isortps      ,& ! (sorted) list of stations with 'surface' pressure data
    isortsu      ,& ! (sorted) list of stations with surface-level data
    isortua      ,& ! (sorted) list of stations with upper-air single-level data
    io    ,jo    ,& ! local  indices of location of observation
    io_tot,jo_tot,& ! global indices of location of observation
    istaspr      ,& ! /  lower left corner of domain
    jstaspr      ,& ! \  containing area of influence
    iendspr      ,& ! /  upper right corner of domain
    jendspr      ,& ! \  containing area of influence
    jrs_tot      ,& ! /  index range for convertor
    jre_tot      ,& ! \  for zonal distances 'gppkmi'
    ista         ,& ! index of observing station

! 9. Observation types in weighted increment arrays 'zwi'
! -------------------------------------------------------

    mxotyp       ,& ! number of observation types
    noiras       ,& ! radiosonde
    noisfc       ,& ! surface-level (SYNOP, SHIP, BUOY)
    noiair       ,& ! aircraft
    noisat       ,& ! satellite (SATOB, ATOVS, MSG)
    noigps       ,& ! GPS

! 10. Other variables
! -------------------

    nhisto          ! data density histogram for surface pressure obs 

! end of data_nudge_spread

!-------------------------------------------------------------------------------

USE src_correl_cutoff,        ONLY :  &
      cutoff_wind_correl        ! areas of influence of observations


!-------------------------------------------------------------------------------
! Local declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE


!===============================================================================

!-------------------------------------------------------------------------------
! Public and Private Subroutines
!-------------------------------------------------------------------------------

CONTAINS

! INCLUDE "ps_spreading.incf"
! INCLUDE "surf_upair_spread.incf"
! INCLUDE "surf_org_spread.incf"
! INCLUDE "upair_org_spread.incf"

!-------------------------------------------------------------------------------
!+ Module procedure for spreading of the 'surface' pressure obs increments
!-------------------------------------------------------------------------------

SUBROUTINE ps_spreading

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure performs the weighting and spreading of the
!   'surface' pressure observation increments, and computes the analysis
!   increments for pressure at the lowest model level.
!
! Method:
!   The area of influence is determined using the tangent cone projection for
!   horizontal distances (cf. Bell et al, 1996). Spreading correction factors
!   to the observation increments are determined. They depend on the applied
!   'temperature correction', and are designed to avoid orographic footprints
!   in pressure and (possibly) temperature analysis increments on isobaric
!   surfaces.
!   Data density at grid points is accounted for by a weighted mean of weighted
!   observation increments according to Benjamin and Seaman (MWR, 1985).
!
! Written by        :  Christoph Schraff, DWD  (original version: 25.06.97)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments: None
! --------------------

! Local parameters:
! ----------------

  REAL    (KIND=wp   )     , PARAMETER  ::  &
    rcps    = 0.0_wp    ! constant part of autoregressive horizontal
                            ! structure function for surface pressure (ps) data

  LOGICAL                  , PARAMETER  ::  &
    lautorg = .TRUE.        ! .FALSE. ==> Cressman instead of autoregressive
                            !             horizontal correlation function

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    istaps           ,& ! index of observing station
    ihisto           ,& ! index for data density histogram
    i      , j       ,& ! loop indices in horizontal direction
    irange , jrange  ,& ! max. distance from the obs. within this domain
    idix   , jdix    ,& ! indices depending on the distance from the obs.
    kwtyp_ps         ,& ! function for weights W for multiple surf. pressure obs
    nstat               ! 

  REAL    (KIND=wp   )     ::  &
    romy             ,& ! 
    zdt              ,& ! timestep relevant for the nudging
    fmultw           ,& ! 0. or 1., depending on multiple observation weighting
!   zpsoi            ,& ! surface pressure observation increment
!   zpsoi1           ,& ! ps obs. incr. of 2nd obs. at current obs. station
!   omyk             ,& ! local part of nudging weight
!   omyk1            ,& ! local part of nudging weight of 2nd obs.
    zpor             ,& ! reciproke of pressure at obs. point
    zmasstr          ,& ! reciproke of total mass affected by the 'temp. corr.'
    zcutop           ,& ! cut-off radius for area of influence
    zcutopl          ,& ! enlarged cut-off radius if geostrophic wind correction
    zcutof           ,& ! cut-off in horiz. correlation scale units
    omysedge         ,& ! value of horiz. correlation at the selected cut-off
!   zpke             ,& ! model pressure at lowest model level
!   cfitop           ,& ! weighted conveyor of geopot. correl. to p- correlation
    zqmass           ,& ! mass affected by 'temp. corr.' devided by 'zmassps'
    zdist            ,& ! distance from the observation location
    zddsp            ,& ! 'zdist' scaled by horizontal correlation scale
    smyps            ,& ! sum of relatively weighted weights
    zgnups           ,& ! nudging coefficient for surface pressure data
    epsw             ,& ! scaled epsilon (very small value)
    pswinc1, pswinc2 ,& ! (partially) weighted obs. incr. of surface pressure
    pswinct, pswint2 ,& ! sum of the 2 partially weighted scalar obs. increments
    omykt            ,& ! sum of the local parts of nudging weights of 2 obs
    omykm            ,& ! omykt *omysc : sum of the 2 correlations (weight)
    omysc2           ,& ! omysc **2
    omyk2t           ,& ! omyk1 **2 + omyk2 **2
    omykh            ,& ! omykt *omysc *omysc
    r1iflpr          ,& ! reciproke of horizontal correlation scale
    r2iflpr          ,& ! r1iflpr **2
    zcutop2, zcutpl2 ,& ! zcutop  **2    ;   zcutopl **2
    zdcutopr         ,& ! c1 / (zcutopl - zcutop)
    zqmass2, zqmass3    ! zqmass  **2    ;   zqmass **3

  LOGICAL                  ::  &
    ledge               ! small gradients for weights near edge of influ. area


! Local (automatic) arrays:
! -------------------------

  LOGICAL                  ::  &
    lsprsv  (npstot)     ! temporal and horiz. domain of node influenced by obs.

  INTEGER (KIND=iintegers) ::  &
    ijspr   (4,npstot)   ! domain containing area of influence

  REAL    (KIND=wp   )     ::  &
    zrcutov (npstot,3)   ! horizontal cut-off radius [km]

  REAL    (KIND=wp   )     ::  &
    zpswi   (istart:iend,jstart:jend) ,& ! sum of weighted observat. increments
    omyps   (istart:iend,jstart:jend) ,& ! sum of spreading weights
    om2ps   (istart:iend,jstart:jend) ,& ! sum of square of spreading weights
    zdist2  (istart:iend,jstart:jend) ,& ! (distance from the obs. location) ^2
    omysc   (istart:iend,jstart:jend) ,& ! horizontal nudging weight (for 1 obs)
    zsprcor (istart:iend,jstart:jend)    ! spreading correction factor to obs.
                                         ! increment; scalor dep on data density

  REAL    (KIND=wp   )     ::  &
    zpke_v  (istart:iend,jstart:jend) ,& ! model pressure at lowest model level
    cfitop_v(istart:iend,jstart:jend)    ! weighted conveyor of geopot. correl. 
                                         ! to p- correlation

  LOGICAL                  ::  &
    larea   (istart:iend,jstart:jend) ,& ! grid pt lies within area of influence
    lareaex (istart:iend,jstart:jend)    ! grid pt lies outside area of influ-
                                         ! ence as defined by namelist cutoff
!
!------------ End of header ----------------------------------------------------


 
!-------------------------------------------------------------------------------
! Begin Subroutine ps_spreading
!-------------------------------------------------------------------------------
  
!-------------------------------------------------------------------------------
!  Section 1: Preliminaries
!-------------------------------------------------------------------------------

  zgnups = gnudg(2)
  IF (gnudg(2) < epsy) zgnups = gnudgsu(2)

  DO   j = jstart, jend
    DO i = istart, iend
      psanai(i,j) = c0
      zpswi (i,j) = c0
      omyps (i,j) = c0
      om2ps (i,j) = c0
      zdist2(i,j) = c0
      zsprcor(i,j) = c0
    ENDDO
  ENDDO

  ledge  = (luvgcor) .AND. (lautorg)
  cfitop_v(istart:iend,jstart:jend) = c1
  kwtyp_ps = kwtyp(nwtyp+2)
  IF (nwtyp == 1)  kwtyp_ps = kwtyp(1)
  fmultw = REAL ( kwtyp_ps - 1, wp )

  IF ((npstot == 0) .AND. ((lfirst) .OR. (rmod( ntstep*dt,c3600 ) < tconbox))  &
                    .AND. (lcroot))                                            &
    PRINT '(''hour '',F3.0,'' : no surface pressure data'')' , acthr

! check if this model sub-domain overlaps with the
! rectangular domains containing the areas of influence
! -----------------------------------------------------

  DO istaps = 1 , npstot
    ista = isortps(istaps)
    io   = iops (ista)
    jo   = jops (ista)
                     zrcutov (istaps,1)  =  cutofsu(2)
    IF (lmlps(ista)) zrcutov (istaps,1)  =  cutofr(2)
    zrcutov (istaps,2)  =  zrcutov(istaps,1) * r1ifps(ista)
    zrcutov (istaps,3)  =  zrcutov(istaps,2)
    IF (ledge)  zrcutov (istaps,3)  =  zrcutov(istaps,2)                       &
                                     + r1ifps(ista) *(c1 + c1/zrcutov(istaps,1))
    jrange  = INT( zrcutov(istaps,3) * gppkmj )
    jstaspr = MAX( INT( jo - jrange ,iintegers) , jrs_tot )
    jendspr = MIN( INT( jo + jrange ,iintegers) , jre_tot )
    irange  = INT( zrcutov(istaps,3) * MAX( gppkmi(jstaspr), gppkmi(jendspr) ) )
!   irange  = INT( zrcutov(istaps,3) * MAX( gppkmi(MIN( jstaspr,jend   ))      &
!                                         , gppkmi(MAX( jendspr,jstart )) ) )
    ijspr (3,istaps) = MAX( INT( jo - jrange ,iintegers) , jstart )
    ijspr (4,istaps) = MIN( INT( jo + jrange ,iintegers) , jend )
    ijspr (1,istaps) = MAX( INT( io - irange ,iintegers) , istart )
    ijspr (2,istaps) = MIN( INT( io + irange ,iintegers) , iend )
    lsprsv  (istaps) = (      (ijspr(3,istaps) <= jend)                        &
                        .AND. (ijspr(4,istaps) >= jstart)                      &
                        .AND. (ijspr(1,istaps) <= iend)                        &
                        .AND. (ijspr(2,istaps) >= istart) )
  ENDDO


!-------------------------------------------------------------------------------
!  Section 2: Weighting and spreading of observation increments
!-------------------------------------------------------------------------------

loop_over_stations:  DO istaps = 1 , npstot

  IF (lsprsv(istaps)) THEN

    ista = isortps(istaps)

! set scalar quantities

    io     = iops     (ista)
    jo     = jops     (ista)
!   io_tot = iops_tot (ista)
!   jo_tot = jops_tot (ista)
    ystid  = ystidps  (ista) (1:ilstidp)

    zcutof  = zrcutov(istaps,1)
    zcutop  = zrcutov(istaps,2)
    zcutopl = zrcutov(istaps,3)
    istaspr = ijspr(1,istaps)
    iendspr = ijspr(2,istaps)
    jstaspr = ijspr(3,istaps)
    jendspr = ijspr(4,istaps)

    pswinc1 = omykps(ista,1) * zoips(ista,1)
    pswinc2 = omykps(ista,2) * zoips(ista,2)
    omykt   = omykps(ista,1) + omykps(ista,2)
    pswinct = pswinc1 + pswinc2
    IF (.NOT. ltips(ista)) THEN
      pswinc1 = omykps(ista,1) * pswinc1
      pswinc2 = omykps(ista,2) * pswinc2
      omyk2t  = omykps(ista,1) *omykps(ista,1) + omykps(ista,2) *omykps(ista,2)
      pswint2 = pswinc1 + pswinc2
    ENDIF

    r1iflpr  = c1 / r1ifps(ista)
    IF (.NOT.lautorg)  r2iflpr = r1iflpr **2
    zcutop2  = zcutop  **2
    zcutpl2  = zcutopl **2
    IF (ledge) zdcutopr = c1 / (zcutopl - zcutop)
    omysedge = (c1 + zcutof) * EXP( -zcutof )
    zpor     = c1 /(zmassps(ista) + ptpstop)
    zmasstr  = c1 / zmassps(ista)

! get area of influence and horizontal distance betw. obs. and target grid pt.

    DO   j = jstaspr, jendspr
      DO i = istaspr, iendspr
        idix  = i - io + ie_tot
        jdix  = j - jo + je_tot
        zdist2 (i,j) = xxd2(idix,j) + yyd2(jdix) - yyd(jdix) *xsalpa(idix,j)
        larea  (i,j) = (zdist2(i,j) <= zcutpl2)
        lareaex(i,j) = (zdist2(i,j) >  zcutop2)
      ENDDO
    ENDDO

! get horizontal weights, and spreading correction factor to obs. increment
! (the latter is designed to avoid orographic footprints in pressure and
!  temperature analysis increments at upper levels.)

    DO   j = jstaspr, jendspr
      DO i = istaspr, iendspr
        IF (larea(i,j)) THEN
          zpke_v(i,j) = p0(i,j,ke) + pp(i,j,ke,nnew)
          zsprcor (i,j) = zpke_v(i,j) * zpor
          IF (ntpscor >= 3) THEN
            cfitop_v(i,j) = ztdpps(ista)* zpke_v(i,j) / t(i,j,ke,nnew)
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    DO   j = jstaspr, jendspr
      DO i = istaspr, iendspr
        IF (larea(i,j)) THEN
          IF (ntpscor >= 1) THEN
            zqmass  = (zpke_v(i,j) - ptpstop) * zmasstr
            IF (MOD(ntpscor,2) == 1) THEN
              zqmass2 = zqmass *zqmass
              zqmass3 = zqmass *zqmass2
              zsprcor (i,j) =  cfitop_v(i,j) * zqmass2                         &
                             * EXP( (c1-zqmass3)*0.125_wp / MAX(c1,zqmass3))
            ELSE
              zsprcor (i,j) =  cfitop_v(i,j) * zqmass * MAX(c1 , zqmass)       &
                             * (c05*(c1+zqmass)) **NINT(SIGN( c1, c1-zqmass ))
            ENDIF
          ENDIF
          zdist        = SQRT( zdist2(i,j) )
          zddsp        = zdist * r1iflpr / rhscal(i,j)
          omysc (i,j)  = (c1 + zddsp) *EXP( -zddsp )
!         omysc (i,j)  = rcps + (c1-rcps) *(c1+zddsp) *EXP( -zddsp )
! if the geostrophic wind correction is applied, avoid strong gradients at the
! edge of the area of influence
! --> wanted: linear decrease from omysc(at zcutop) to zero within lareaex 
!             (i.e. for zcutop <= zdist <= zcutopl)
!     problem: omysc(at zcutop) /= omysedge and is not known at zdist > zcutop
!              (unless rhscal=1 everywhere)
! --> assume: rhscal(at zcutop) = rhscal(i,j) (no abrupt change of data density)
          IF (lareaex(i,j)) omysc(i,j) =  (zcutopl - zdist)* zdcutopr          &
                                        * (c1 + zcutof /rhscal(i,j))           &
                                        * EXP( -zcutof /rhscal(i,j) )
!                                       * (c1-rcps) + rcps)
!         IF (lareaex(i,j)) omysc(i,j) = (zcutopl - zdist) * zdcutopr * omysedge
!         IF (ledge) omysc = omysc *(zcutopl - MAX( zdist ,zcutop )) *zdcutopr
          IF (.NOT.lautorg) omysc(i,j) = MAX( c0 , c1 - zdist2(i,j) *r2iflpr ) &
                                        /        ( c1 + zdist2(i,j) *r2iflpr )
        ENDIF
      ENDDO
    ENDDO

! Prepare the lateral spreading

! Execute lateral spreading of obs. increments, determine weights

    IF (ltips(ista)) THEN
      DO   j = jstaspr, jendspr
        DO i = istaspr, iendspr
          IF (larea(i,j)) THEN
            omykm  =  omysc(i,j) * omykt
            omykh  =  omysc(i,j) *(omykm + fmultw)
            omyps (i,j) = omyps(i,j) + omykm
            om2ps (i,j) = om2ps(i,j) + omykm *omykm
            zpswi (i,j) = zpswi(i,j) + omykh *pswinct *zsprcor(i,j)
          ENDIF
        ENDDO
      ENDDO
    ELSE
      DO   j = jstaspr, jendspr
        DO i = istaspr, iendspr
          IF (larea(i,j)) THEN
            omysc2 =  omysc(i,j) * omysc(i,j)
            omyps (i,j) = omyps(i,j) + omysc(i,j) *omykt
            om2ps (i,j) = om2ps(i,j) + omysc2     *omyk2t
            zpswi (i,j) = zpswi(i,j) + omysc2     *pswint2 *zsprcor(i,j)       &
                                     + omysc(i,j) *pswinct *zsprcor(i,j) *fmultw
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    IF (      (lwonl2) .AND. (jstaspr <= jonl2) .AND. (jendspr >= jonl2)       &
                       .AND. (istaspr <= ionl2) .AND. (iendspr >= ionl2)       &
        .AND. (    (lfirst)                                                    &
              .OR. (    (rmod( ntstep*dt,21600.0_wp ) < tconbox)           &
                  .AND. (REAL(ntstep,wp) >= REAL(MIN(nudgend,nstop),wp)-21599.9_wp/dt)))) THEN
      zdist = SQRT( zdist2(ionl2,jonl2) )
      WRITE( nupr,'(''ps_spread '',A ,'', weights'',2F5.2,F6.2                 &
                  &,'', (wght.) incr'',F8.2,2F7.1,'', dist'',F6.0,L2,F5.2)' )  &
             ystid, omykps(ista,1), omykps(ista,2), omyps(ionl2,jonl2)         &
                  , zoips (ista,1), zoips (ista,2), zpswi(ionl2,jonl2)         &
                  , zdist, larea(ionl2,jonl2), zsprcor(ionl2,jonl2)
    ENDIF
!   IF ((ntstep == nstop) .AND. (jstaspr <= jonl2) .AND. (jendspr >= jonl2)    &
!         .AND. (lwonl2) .AND. (istaspr <= ionl2) .AND. (iendspr >= ionl2)) THEN
!     PRINT *,'ps_spr ', ystid, ntstep, omyps(ionl2,jonl2), zpswi(ionl2,jonl2)
!     PRINT *,'ps_spr ', ystid, ntstep, omykps(ista,1), omykps(ista,2)
!     PRINT *,'ps_spr ', ystid, ntstep, zoips (ista,1), zoips (ista,2)
!     PRINT *,'ps_spr ', ystid, ntstep, io, jo, lmlps(ista), ltips(ista)       &
!                             , larea(ionl2,jonl2), zsprcor(ionl2,jonl2)
!     PRINT *
!   ENDIF

  ENDIF

ENDDO loop_over_stations


!-------------------------------------------------------------------------------
!  Section 3.0: Determine factor 'fs' used to scale the horizontal correlation 
!               scale depending on the sum of individual weights 'Sw' as a proxy
!               for data density
!               ( fs = 1                             for  Sw < 1  ,
!                 fs = k + (1-k)* 1/(1+(Sw-1)**2)    for  Sw > 1  ,
!                 fs = k + (1-k)* (1+Sw/r)*exp(-Sw/r)    for  Sw > 1  ,
!                 k being the minimum scale factor to the correlation scale )
!-------------------------------------------------------------------------------

! romy    =  2.0_wp
! rhfpsdd =  0.7_wp
! romy    =  4.0_wp
! rhfpsdd =  0.5_wp
  romy    =  3.0_wp
  IF (rhfpsdd < c1-epsy) THEN
    DO   j = jstart, jend
      DO i = istart, iend
!       zsprcor (i,j) = MAX( c0 , omyps(i,j) - c1 ) / romy
! (defining rhscal in this way and introducing it in the horizontal correlation
!  scale / function changes the weight function by less than 0.5 % , given one
!  single observation.)
        zsprcor (i,j) =           omyps(i,j)        / romy
        zsprcor (i,j) = (1 + zsprcor(i,j)) * EXP( - zsprcor(i,j) )
        zsprcor (i,j) = rhfpsdd + (c1 - rhfpsdd) *zsprcor(i,j)
! blending of the previous and current scale factor diminishes oscillations
! (wavelength: 2*tconbox) of the data density and the scale factor itself
        rhscal  (i,j) = 0.75_wp* zsprcor(i,j) + 0.25_wp* rhscal(i,j)
      ENDDO
    ENDDO
    nhisto = 0
    DO   j = jstart, jend
      DO i = istart, iend
        ihisto          = MIN( NINT( omyps(i,j) - 0.4999_wp ) + 1 , 50 )
        IF ((ihisto > 10) .AND. (ihisto <= 20))  ihisto = (ihisto - 1) /5  + 9
        IF ((ihisto > 20) .AND. (ihisto <= 50))  ihisto = (ihisto - 1) /10 + 11
        nhisto (ihisto) = nhisto(ihisto) + 1
      ENDDO
    ENDDO
!   IF (ntstep*dtdeh < 1.001_wp)  PRINT '("Hit:",10I4,I5,7I4)', nhisto
    IF (lwonl) WRITE( nupr,'("ps: Wtot, rhscal",2F9.4)') omyps(ionl,jonl)      &
                                                       , rhscal(ionl,jonl)
  ENDIF

!-------------------------------------------------------------------------------
!  Section 3: Computation of the "net" (weighted mean) observation increments,
!             "net" nudging weights, and analysis increments
!-------------------------------------------------------------------------------
 
  epsw = epsy / MAX( zgnups , epsy )

  IF (l2tls) THEN
    zdt = dt
  ELSE
    zdt = dt2
  ENDIF

  IF (kwtyp_ps == 1) THEN
    DO   j = jstart, jend
      DO i = istart, iend

        IF (omyps(i,j) > epsw) THEN
! 'smyps': sum over all ps-obs of my(ps) . my(ps) itself is 'gnudg' times a
!                                          relatively weighted weight.
          smyps       = zgnups * om2ps(i,j) / omyps(i,j)
          omyps (i,j) = smyps * faclbc(i,j,1)
          psanai(i,j) =   zdt *omyps(i,j) / (c1 + zdt *omyps(i,j))             &
                        * zpswi(i,j) / om2ps(i,j)
        ELSE
          psanai(i,j) = c0
        ENDIF
      ENDDO
    ENDDO

  ELSEIF (kwtyp_ps == 2) THEN
    DO   j = jstart, jend
      DO i = istart, iend

        IF (omyps(i,j) > epsw) THEN
! 'smyps': sum over all ps-obs of my(ps) . my(ps) itself is 'gnudg' times a
!                                          relatively weighted weight.
          smyps       = zgnups * (omyps(i,j) + om2ps(i,j))                     &
                               / (c1         + omyps(i,j))
          smyps       = smyps * faclbc(i,j,1)
          psanai(i,j) =   zdt *smyps / (c1 + zdt *smyps)                       &
                        * zpswi(i,j) / (omyps(i,j) + om2ps(i,j))

          omyps (i,j) = smyps * faclbc(i,j,1)
        ELSE
          psanai(i,j) = c0
        ENDIF
      ENDDO
    ENDDO
  ENDIF


  IF ((lwonl2) .AND. (     (lfirst)                                            &
                      .OR. (rmod( ntstep*dt,21600.0_wp ) < tconbox)))      &
    WRITE( nupr,'(''End of ps_ana_i: weight sqr'',F7.2                         &
                &,'', net weight'',F9.6,'', ana incr'',2F8.4)' )               &
           om2ps(ionl2,jonl2), omyps(ionl2,jonl2), psanai(ionl2,jonl2)         &
         ,faclbc(ionl2,jonl2,1)
  IF ((lwonl2) .AND. (ntstep*dtdeh < 2.001_wp))                            &
    WRITE( nupr,'(''ps-ana.incr. '',2I4, 4X, F9.3)' )                          &
           ionl2, jonl2, psanai(ionl2,jonl2)

  IF (ldump_ascii) THEN
    ! flush YUPRINT file
    IF ((lwonl2) .AND. (     (lfirst)                                          &
                        .OR. (rmod( ntstep*dt,21600.0_wp ) < tconbox))) THEN
      CLOSE (nupr)
      OPEN  (nupr,FILE=yuprint,FORM='FORMATTED',STATUS='OLD'                   &
                              ,POSITION='APPEND',IOSTAT=nstat)
      IF (nstat /= 0) PRINT '("OPENING OF FILE yuprint FAILED in ps_spreading")'
    ENDIF
  ENDIF


!-------------------------------------------------------------------------------
! End of module procedure ps_spreading
!-------------------------------------------------------------------------------

END SUBROUTINE ps_spreading


!-------------------------------------------------------------------------------
!+ Module procedure for spreading of single-level observation increments
!-------------------------------------------------------------------------------

SUBROUTINE surf_upair_spread ( lsurf , ivar , k , omymx )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure performs the weighting and spreading (extrapolation)
!   of (at most 2 surface-level or upper-air) single-level observation
!   increments from one observing station.
!
! Method:
!   Spreading of upper-air single-level data increments to the grid pts. with
!   the largest weighting along specified surfaces according to 'modespr' (see
!   below). Computation of the sum of 'relative nudging weights' by applying
!   the temporal and spatial correlations.
!   Option for non-isotropic lateral correlations (see below).
!   Humidity: Computation of specific humidity ('qd') increments from the
!             spreaded relative humidity ('RH') increments, so that the nudging
!             pulls towards observed 'RH' instead of observed 'qd'.
!   Mode of spreading:   'modespr' =
!     0: 'spreading along model levels' :
!        - vertical correlation as Gaussian of height differences ('dz') between
!          the obs. point and the grid pt. at level 'k' at the obs. location
!        - non-isotropy as Gaussian of pot. temp. differences on model levels 
!     1: 'spreading along horizontal levels' :
!        - vertical correlation as Gaussian of height differences ('dz') between
!          the obs. point and the target grid pt. at level 'k'
!        - non-isotropy as Gaussian of pot. temp. differences on horiz. levels
!        This is equivalent to mode '0' where the model levels are horizontal.
!     2: 'spreading along isentropic surfaces' :
!        - vertical correlation as Gaussian of potential temperature differences
!          between the obs. point and the target grid pt. at level 'k'
!        - non-isotropy as Gaussian of height differences on isentropic surfaces
!        Optionally, spreading may be along model levels for target grid points
!        above a specified model level. Vertical correlation is still a Gaussian
!        of potential temperature differences, and non-isotropic weights are not
!        possible at horizontal model levels in this case.
!     correlation functions: - in z     : EXP( -(g(z1-z2)/r*Tv1)**2 / dzscal) )
!                            - in theta : EXP( -(theta1-theta2)**2 / dthscal) )
!
! Written by        :  Christoph Schraff, DWD  (original version: 16.07.97)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE


! Subroutine arguments:
! --------------------

  LOGICAL                 , INTENT (IN)         ::       &
    lsurf               ! spreading of surface-level (not upper-air) obs. incr.

  INTEGER (KIND=iintegers), INTENT (IN)         ::       &
    ivar             ,& ! index of spreaded quantity: 1=(u,v); 3=T; 4=RH
    k                   ! index of current vertical model level

  REAL    (KIND=wp   ),     INTENT (INOUT)      ::       &
    omymx               ! max. spatial weight


! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    modespr          ,& ! mode of spreading, see above: 'Method'
    msi              ,& ! index of spreading parameter
    mni              ,& ! index of parameter defining non-isotropic weights
    meqi             ,& ! index of parameter defining equidistant points
    ityp             ,& ! index of weighted increment array
    mscalar          ,& ! if spreaded quantity is scalar : 0 , else (wind): 0
    iabl             ,& ! index of obs type for ABL weighting: 1=upair, 2=sfc
    ivrs             ,& ! index of spreaded quantity: 1=(u,v); 2=T; 3=RH
    ivay             ,& ! index of second component of spreaded vector quantity
!   ista             ,& ! index of observing station
    ixc    , ixw     ,& ! loop index over influenced grid points
    i      , j       ,& ! loop indices in horizontal direction
!   io     , jo      ,& ! local indices of location of observation
!   istaspr, jstaspr ,& ! lower left corner of domain containing area of infl.
!   iendspr, jendspr ,& ! upper right corner of domain containing area of infl
!   ijstot           ,& ! # grid pts. within domain containing area of influence
    kobtyp           ,& ! CMA observation type
    kcdtyp           ,& ! CMA observation code type
    mvb                 ! index of equidistant point

  REAL    (KIND=wp   )     ::  &
    rdsprs           ,& ! vertical correlation scale as by namelist
    rdsprsb, rdsprst ,& ! lower / upper vertical correlation scale
    rdspri           ,& ! vert. correl. scale for weighting the obs. incr. level
    rdsprni          ,& ! square of Gauss. non-isotropy radius,'zspr(mni)'-units
    rdsprsr          ,& ! reciproke of 'rdsprs'
    rdsprsbr,rdsprstr,& ! reciproke of 'rdsprsb', 'rdsprst'
    rdsprnir         ,& ! reciproke of 'rdsprni'
    zqts             ,& ! quotient of the tropospheric to the stratospheric mean
                        ! vertical gradient of the spreading parameter
    zcutni           ,& ! 'rdsprni' < 'zcutni' : correlation is non-isotropic
    zsprobs          ,& ! value of spreading param. (z, TH) at the obs.
    ztp              ,& ! value of spreading param. (z, TH) at (std atm) tropop.
    zdsobm           ,& ! diff. of spreading param. (z, TH) betw. obs. & grid p.
    zdsobt           ,& ! diff. of spreading param. (z, TH) betw. obs. & tropop.
    zsproni          ,& ! for non-isotr. w.: TH or z at spread. surf. & obs loc.
    zlopf            ,& ! quotient for interpolation
    zispr  , disprr     ! define vertic. equidist. pts, for non-isotropic weight

  REAL    (KIND=wp   )     ::  &
    fmultw           ,& ! 0. or 1., depending on multiple observation weighting
    fftho            ,& ! quotient pot. temperature / temperature at obs. pt.
    uincr1 , uincr2  ,& ! weighted obs. incr. of zonal wind (for 2 obs at 1 sta)
    vincr1 , vincr2  ,& ! wt. obs. incr. of meridional wind (for 2 obs at 1 sta)
    tincr1 , tincr2  ,& ! obs. incr. of pot. temperature  (for 2 obs at 1 sta.)
    rhincr1, rhincr2 ,& ! obs. incr. of relative humidity (for 2 obs at 1 sta.)
    rhincw1, rhincw2 ,& ! weighted obs. incr. of rel. hum.(for 2 obs at 1 sta.)
!   zewcor1, zewcor2 ,& ! targeted      vapour pressure   (for 2 obs at 1 sta.)
!   zqdcor1, zqdcor2 ,& ! targeted      specific humidity (for 2 obs at 1 sta.)
    uincrt , vincrt  ,& ! sum of the 2 partially weighted wind comp. obs. incr.
    uinct2 , vinct2  ,& ! sum of the 2 partially squared weighted wind obs. incr
    zuoem  , zvoem   ,& ! sum of the 2 spreaded obs. incr. of wind vector
    zuoem2 , zvoem2  ,& ! sum of the 2 squared spreaded obs. incr. of wind
    tincrt , tinct2  ,& ! sum of the 2 weighted (**2) obs. incr. of temperature
    rhincrt, rhinct2 ,& ! sum of the 2 weighted (**2) obs. incr of rel. humidity
    wtql1  , wtql2   ,& ! temporal nudging weight * quality factor * rel. coeff.
    omyvf            ,& ! local vertical correlation (if (lsigma))
    omyk1  , omyk2   ,& ! local part of the nudging weights
    omykf            ,& ! local part of the vertical nudging weight
    omykfij          ,& ! part of the vertical nudging weight
    omykt            ,& ! omyk1 + omyk2 : sum of the 2 local weights
    omykm            ,& ! omykt *omysc : sum of the 2 correlations (weight)
    omykh            ,& ! omykt *omysc *omysc
    omysc2           ,& ! omysc **2
!   omyk12 , omyk22  ,& ! omyk1 **2      ;   omyk2 **2
    omyk2t           ,& ! omyk1 **2 + omyk2 **2
    omytt            ,& ! weight used for spreading of potential temperature
    omykfty          ,& ! weight depending on observation type
    omysci              ! isotropic (part of) spatial correlations (weights)


  LOGICAL                  ::  &
    lupair           ,& ! spreading of upper-air (not surface-level) obs. incr.
    lsigma           ,& ! spreading along model levels
    ladjtp           ,& ! spreading of upper-air obs. incr. along isentropes
    lnoniso          ,& ! non-isotropic horizontal correlation
    lti2                ! linear temporal interpolation as time window


! Local (automatic) arrays:
! -------------------------
  REAL    (KIND=wp   )     ::  &
    omysc   (istart:iend,jstart:jend) ,& ! spatial correlations (weights)
                                         ! (without non-div. corr)
    zniseq  (mxispr)    ! equidistant profile for non-isotropic weights

!------------ End of header ----------------------------------------------------

 
!-------------------------------------------------------------------------------
! Begin Subroutine surf_upair_spread
!-------------------------------------------------------------------------------
  
!-------------------------------------------------------------------------------
!  Section 1: General preparations
!-------------------------------------------------------------------------------
 
  lupair  = (.NOT. lsurf)
  ivrs    = MAX( ivar - 1 , 1 )
  mscalar = MIN( ivar - 1 , 1 )

  IF (lsurf)  modespr = msprpsu
  IF (lupair) modespr = msprpar
  msi  = MAX( modespr , i1 )
  mni  = 3 - msi
  meqi = mni + 1

  IF (lsurf) THEN
    kobtyp = kobtysu(ista)
    kcdtyp = kcdtysu(ista)
  ELSE
    kobtyp = kobtyua(ista)
    kcdtyp = kcdtyua(ista)
  ENDIF
  ityp = isetyp0
  DO ixw = SUM( niwtyp ) , 1 , -1
    IF (     ((iwtyp(ixw) > 0) .AND. ( iwtyp(ixw) == kobtyp))                  &
        .OR. ((iwtyp(ixw) < 0) .AND. (-iwtyp(ixw) == kcdtyp)))                 &
      ityp = isetyp(ixw)
  ENDDO
  fmultw = REAL ( kwtyp(ityp) - 1, wp )
! IF (      (ntstep == 0) .AND. (MOD( ista, 10 ) == 0)  &
!     .AND. (jo <= jend) .AND. (jo >= jstart))  &
!   WRITE (0,*) 'ZZml2 ', ista, kobtyp, kcdtyp, ityp, kwtyp(ityp), lsurf


  lsigma =     ( modespr == 0)                                                 &
          .OR. ((modespr == 1) .AND. (k <= ktp))                               &
          .OR. ((modespr == 2) .AND. (k < ktth))
!         .OR. ((modespr == 2) .AND. (     ((.NOT. ltopth) .AND. (k < ktopth)) &
!                                     .OR. ((ltopth) .AND. (k <= ktp))))

! get the spreading parameter at the obs. point (and at the standard atmospheric
! tropopause if required), and the Gaussian radii of the vertical weight funct.
!   (IF (modespr == 2) THEN account for greater stability in the stratosphere)

  IF (lsurf) THEN
    zsprobs = zsposu (ista)
    rdsprs  = SQRT( vcorlsu(ivar) )
    IF (modespr <= 1) rdsprs = rdsprs * zrtdgsu(ista)
    ztp     = 100000.0_wp
  ELSEIF (lupair) THEN
    zsprobs = zsprob (ista)
    ztp     = zsprtp (ista)
    omykf   = c1
    IF (modespr <= 1) THEN
      rdsprs = SQRT( vcorls(ivar) ) * zrtdgua(ista)
      rdspri = SQRT( MAX( vcorls(ivar) , 0.04_wp ) )
    ELSEIF (lsigma) THEN
      rdsprs = SQRT( vcorls(ivar) )
      IF (zsprua(ista,k) > zsprtp(ista)) THEN
        omykf   = EXP( - (MAX( zsprtp(ista) -zsprobs, c0 ) /rdsprs) **2 )
        zsprobs = MAX( zsprobs , zsprtp(ista) )
        rdsprs  = qst * rdsprs
      ELSE
        omykf   = EXP( - (MAX( zsprobs -zsprtp(ista), c0 ) /rdsprs) **2 / qst )
        zsprobs = MIN( zsprobs , zsprtp(ista))
      ENDIF
      rdspri = SQRT( MAX( vcorls(ivar) , 33.0_wp ) )
    ELSE
      rdsprs = SQRT( vcorls(ivar) )
      IF (zsprobs > zsprtp(ista)) rdsprs = qst * rdsprs
      omykf  = EXP( - ((zsprobs - zsprtp(ista)) /rdsprs) **2 )
      rdspri = SQRT( MAX( vcorls(ivar) , 33.0_wp ) )
    ENDIF
    rdsprsb = rdsprs * fcorlua(ista,1,ivrs)
    rdsprst = rdsprs * fcorlua(ista,2,ivrs)
  ENDIF

!-------------------------------------------------------------------------------
!  Section 2: (Horizontally) local contributions to the weights:
!             temporal weights (implicitly includes threshold quality control),
!             quality factors, relative nudging coefficient, and (if spreading
!             along model levels) vertical weight
!-------------------------------------------------------------------------------

  omyvf  =  c1
  IF (lsurf) THEN
    wtql1  =  zwtsu (ista,1,ivrs)  *  zqualsu (ista,1,ivrs)                    &
            * gnudgsu(ivar) / MAX( gnudg(ivar) ,epsy )
    wtql2  =  zwtsu (ista,2,ivrs)  *  zqualsu (ista,2,ivrs)                    &
            * gnudgsu(ivar) / MAX( gnudg(ivar) ,epsy )
    IF (lsigma) omyvf = EXP( -((zsprobs -zsprsu(ista,k)) /rdsprs)**2 )
    omyvf = omyvf * EXP( -(zpblsu(ista,k) /vpblsu(ivar))**2 )
  ELSEIF (lupair) THEN
    IF (kcdtyp == nmodes) THEN
      wtql1  =  zwtua (ista,1,ivrs)  *  zqualua (ista,1,ivrs)                  &
              * gnudgms(ivar) / MAX( gnudg(ivar) ,epsy )
      wtql2  =  zwtua (ista,2,ivrs)  *  zqualua (ista,2,ivrs)                  &
              * gnudgms(ivar) / MAX( gnudg(ivar) ,epsy )
    ELSE
      wtql1  =  zwtua (ista,1,ivrs)  *  zqualua (ista,1,ivrs)                  &
              * gnudgar(ivar) / MAX( gnudg(ivar) ,epsy )
      wtql2  =  zwtua (ista,2,ivrs)  *  zqualua (ista,2,ivrs)                  &
              * gnudgar(ivar) / MAX( gnudg(ivar) ,epsy )
    ENDIF
    IF (lsigma)                                                                &
      omyvf = omykf * EXP( -( MIN( (zsprobs-zsprua(ista,k)) /rdsprst ,c0 )     &
                             +MAX( (zsprobs-zsprua(ista,k)) /rdsprsb ,c0 ))**2 )
    omyvf = omyvf * EXP( - (zvidua(ista) /rdspri) **2 )
  ENDIF
  omyk1  =  omyvf * wtql1
  omyk2  =  omyvf * wtql2
  IF (omyk1 <= epsy) omyk1 = c0
  IF (omyk2 <= epsy) omyk2 = c0

!-------------------------------------------------------------------------------
!  Section 3: Preparations for lateral spreading: get quantities used for the
!             vertical weighting and non-isotropic correction, and determine
!             whether linear temporal interpolation is applied.
!-------------------------------------------------------------------------------
 
  omymx = c0

  IF (modespr <= 1) THEN
    zispr  = xezispr
    disprr = c1 / dezispr
    zcutni = vcutnit
  ELSE
    zispr  = xthispr
    disprr = c1 / dthispr
    zcutni = vcutnip
  ENDIF
  zqts = c1 / qst

  rdsprsr  = c1 /      rdsprs
  IF (lsurf) THEN
    rdsprni = vcsnisu(ivar)
    lnoniso = (rdsprni <= zcutni) .AND. (lnissu)
    IF (modespr == 2) rdsprni = rdsprni * zrtdgsu(ista) **2
    iabl = 2
!   lti2 = (ltisu(ista,ivrs)) .AND. (omyk1 >= epsy) .AND. (omyk2 >= epsy)
    lti2 = (ltisu(ista,ivrs))
    rdsprstr = rdsprsr
    rdsprsbr = rdsprsr
  ELSEIF (lupair) THEN
    rdsprni = vcsni  (ivar)
    lnoniso = (rdsprni <= zcutni) .AND. (lnisua)
    IF (modespr == 2) rdsprni = rdsprni * zrtdgua(ista) **2
    IF ((modespr <= 1) .AND. (zsprua(ista,k) > ztp)) rdsprni = qst**2 * rdsprni
    iabl = 1
!   lti2 = (ltiua(ista,ivrs)) .AND. (omyk1 >= epsy) .AND. (omyk2 >= epsy)
    lti2 = (ltiua(ista,ivrs))
    rdsprsbr = c1 / MAX( rdsprsb , epsy )
    rdsprstr = c1 / MAX( rdsprst , epsy )
  ENDIF
  rdsprnir = c1 / MAX( rdsprni , epsy )

! If local weights are zero, further computations (spreading) will be skipped
  IF ((omyk1 <= epsy) .AND. (omyk2 <= epsy)) ivrs = 0

  IF ((lnoniso) .AND. (.NOT. lsigma) .AND. (ivrs > 0)) THEN
    DO mvb = 1 , mxispr
      IF (lsurf)  zniseq(mvb) = znissq(ista,mvb)
      IF (lupair) zniseq(mvb) = znisuq(ista,mvb)
    ENDDO
  ENDIF
  ladjtp = (lupair) .AND. (modespr == 2)

!-------------------------------------------------------------------------------
!  Section 4: Spatial correlations (weights)
!             (except for the non-divergent correction with the 2-dim.
!              horizontal wind vector correlation)
!-------------------------------------------------------------------------------

  IF ((lnoniso) .AND. (ivrs > 0)) THEN
    IF (lsigma) THEN
      zdsobm  = c0
      IF (lsurf ) zsproni = znissu(ista,k)
      IF (lupair) zsproni = znisua(ista,k)
      omykfij = c1
!CDIR NODEP,VOVERTAKE,VOB
      DO ixc = 1 , ncutof(1)
        i = icutof(ixc,1)
        j = jcutof(ixc,1)
        omysc (i,j) =  (c1 + mscalar *zdds(i,j,ivrs))               * omykfij  &
                         * EXP(- ( MIN( zdsobm * rdsprstr , c0 )               &
                                  +MAX( zdsobm * rdsprsbr , c0 ))**2           &
                               - (zsproni -zspr(i,j,mni))**2 * rdsprnir        &
                               - zdds(i,j,ivrs)  )                             &
                         * zablfc(i,j,ivrs,iabl)
        omymx = MAX( omysc(i,j) , omymx )
      ENDDO
    ELSE
      zdsobt = zsprobs - ztp
!CDIR NODEP,VOVERTAKE,VOB
      DO ixc = 1 , ncutof(1)
        i = icutof(ixc,1)
        j = jcutof(ixc,1)
        omykfij = c1
        zdsobm = zsprobs - zspr(i,j,msi)
        IF ((ladjtp) .AND. (zdsobt*(zspr(i,j,msi)-ztp) < epsy)) THEN
          omykfij = omykf
          zdsobm = ztp - zspr(i,j,msi)
          IF (zdsobt < c0) zdsobm = zdsobm * zqts
          IF (zdsobt > c0) zdsobm = zdsobm * qst
        ENDIF
        mvb = INT( (zspr(i,j,meqi) - zispr) * disprr ) + 1
        mvb = MIN( MAX( mvb , i1 ) , mxispr - 1 )
        zlopf   = (zsdnis(mvb,msi) - zspr(i,j,msi)) * zsdnid(mvb,msi)
        zsproni = (c1-zlopf) *zniseq(mvb) + zlopf*zniseq(mvb+1)
        omysc (i,j) =  (c1 + mscalar *zdds(i,j,ivrs))              * omykfij   &
                     * EXP(- ( MIN( zdsobm * rdsprstr , c0 )                   &
                              +MAX( zdsobm * rdsprsbr , c0 ))**2               &
                           - (zsproni -zspr(i,j,mni))**2 * rdsprnir            &
                           - zdds(i,j,ivrs)  )                                 &
                     * zablfc(i,j,ivrs,iabl)
        omymx = MAX( omysc(i,j) , omymx )
      ENDDO
    ENDIF
  ELSEIF ((.NOT. lsigma) .AND. (ladjtp) .AND. (ivrs > 0)) THEN
    zdsobt = zsprobs - ztp
!CDIR NODEP,VOVERTAKE,VOB
    DO ixc = 1 , ncutof(1)
      i = icutof(ixc,1)
      j = jcutof(ixc,1)
      omykfij = c1
      zdsobm  = zsprobs - zspr(i,j,msi)
      IF (zdsobt*(zspr(i,j,msi)-ztp) < epsy) THEN
        omykfij = omykf
        zdsobm = ztp - zspr(i,j,msi)
        IF (zdsobt < c0) zdsobm = zdsobm * zqts
        IF (zdsobt > c0) zdsobm = zdsobm * qst
      ENDIF
      omysc (i,j) =  (c1 + mscalar *zdds(i,j,ivrs))           * omykfij        &
                   * EXP(- ( MIN( zdsobm * rdsprstr , c0 )                     &
                            +MAX( zdsobm * rdsprsbr , c0 ))**2                 &
                         - zdds(i,j,ivrs)  )                                   &
                   * zablfc(i,j,ivrs,iabl)
      omymx = MAX( omysc(i,j) , omymx )
    ENDDO
  ELSEIF ((.NOT. lsigma) .AND. (lupair) .AND. (ivrs > 0)) THEN
!CDIR NODEP,VOVERTAKE,VOB
    DO ixc = 1 , ncutof(1)
      i = icutof(ixc,1)
      j = jcutof(ixc,1)
      zdsobm = zsprobs - zspr(i,j,msi)
      omysc (i,j) =  (c1 + mscalar *zdds(i,j,ivrs))                            &
                   * EXP(- ( MIN( zdsobm * rdsprstr , c0 )                     &
                            +MAX( zdsobm * rdsprsbr , c0 ))**2                 &
                         - zdds(i,j,ivrs)  )                                   &
                   * zablfc(i,j,ivrs,iabl)
      omymx = MAX( omysc(i,j) , omymx )
    ENDDO
  ELSEIF ((.NOT. lsigma) .AND. (ivrs > 0)) THEN
!CDIR NODEP,VOVERTAKE,VOB
    DO ixc = 1 , ncutof(1)
      i = icutof(ixc,1)
      j = jcutof(ixc,1)
      omysc (i,j) =  (c1 + mscalar *zdds(i,j,ivrs))                            &
                   * EXP(- ((zsprobs - zspr(i,j,msi)) *rdsprsr)**2             &
                         - zdds(i,j,ivrs)  )                                   &
                   * zablfc(i,j,ivrs,iabl)
      omymx = MAX( omysc(i,j) , omymx )
    ENDDO
  ELSEIF (ivrs == 1) THEN
!CDIR NODEP,VOVERTAKE,VOB
    DO ixc = 1 , ncutof(1)
      i = icutof(ixc,1)
      j = jcutof(ixc,1)
      omysc (i,j) =  EXP(- zdds(i,j,ivrs)  )                                   &
                   * zablfc(i,j,ivrs,iabl)
      omymx = MAX( omysc(i,j) , omymx )
    ENDDO
  ELSEIF (ivrs >= 2) THEN
!CDIR NODEP,VOVERTAKE,VOB
    DO ixc = 1 , ncutof(1)
      i = icutof(ixc,1)
      j = jcutof(ixc,1)
!     omysc (i,j) =  (c1 + mscalar *zdds(i,j,ivrs))
      omysc (i,j) =  (c1 +          zdds(i,j,ivrs))                            &
                   * EXP(- zdds(i,j,ivrs)  )                                   &
                   * zablfc(i,j,ivrs,iabl)
      omymx = MAX( omysc(i,j) , omymx )
    ENDDO
  ENDIF

! IF ((lwonl) .AND. (k == ke) .AND. (ivrs == 2) .AND. (ntstep >= nstop-1)) THEN
!   PRINT *,'zspr ',io, jo, ystid, zspr(ionl,jonl,msi), zsprobs                &
!                                , lcutof(ionl,jonl,1)
!   PRINT *,'zdds ',io, jo, ystid, zdds(ionl,jonl,ivrs), omysc(ionl,jonl)
! ENDIF

!-------------------------------------------------------------------------------
!  Section 5: Spreading of obs. incr. of generalized relative humidity
!-------------------------------------------------------------------------------

  IF (ivrs == 3) THEN
    rhincr1 = c0
    rhincr2 = c0
    IF (lsurf) THEN
      IF (omyk1 >= epsy) rhincr1 = xoisu(ista,1,ivar)
      IF (omyk2 >= epsy) rhincr2 = xoisu(ista,2,ivar)
    ELSE
      IF (omyk1 >= epsy) rhincr1 = xoiua(ista,1,ivar)
      IF (omyk2 >= epsy) rhincr2 = xoiua(ista,2,ivar)
    ENDIF
!   IF (lti2) THEN
!     rhincw1 = omyk1 * rhincr1
!     rhincw2 = omyk2 * rhincr2
!   ELSE
!     omyk12  = omyk1 **2
!     omyk22  = omyk2 **2
!     rhincw1 = omyk12 * rhincr1
!     rhincw2 = omyk22 * rhincr2
!     omyk2t  = omyk12 + omyk22
!   ENDIF
!   omykt   = omyk1 + omyk2
!   rhincrt = rhincw1 + rhincw2

    rhincw1 = omyk1 * rhincr1
    rhincw2 = omyk2 * rhincr2
    omykt   = omyk1 + omyk2
    rhincrt = rhincw1 + rhincw2
    IF (.NOT. lti2) THEN
      rhincw1 = omyk1 * rhincw1
      rhincw2 = omyk2 * rhincw2
      omyk2t  = omyk1 *omyk1 + omyk2 *omyk2
      rhinct2 = rhincw1 + rhincw2
    ENDIF

    IF (lti2) THEN
!CDIR NODEP,VOVERTAKE,VOB
      DO ixc = 1 , ncutof(1)
        i = icutof(ixc,1)
        j = jcutof(ixc,1)
!       zewcor1 = zewmod(i,j) + zewsat(i,j) *rhincr1
!       zewcor2 = zewmod(i,j) + zewsat(i,j) *rhincr2
!       zqdcor1 = rdv *zewcor1 / (zpk(i,j) - o_m_rdv *zewcor1 )
!       zqdcor2 = rdv *zewcor2 / (zpk(i,j) - o_m_rdv *zewcor2 )
!       rhincrt =  omyk1 *(zqdcor1 - zqdmod(i,j))                              &
!                + omyk2 *(zqdcor2 - zqdmod(i,j))
        omykm  =  omysc(i,j) * omykt
        omykh  =  omysc(i,j) *(omykm + fmultw)
        omy (i,j,ivrs,ityp) = omy (i,j,ivrs,ityp) + omykm
        om2 (i,j,ivrs,ityp) = om2 (i,j,ivrs,ityp) + omykm *omykm
        zwi (i,j,ivar,ityp) = zwi (i,j,ivar,ityp) + omykh *rhincrt
      ENDDO
    ELSE
!CDIR NODEP,VOVERTAKE,VOB
      DO ixc = 1 , ncutof(1)
        i = icutof(ixc,1)
        j = jcutof(ixc,1)
!       rhincrt =  omyk1 *omyk1 *(zqdcor1 - zqdmod(i,j))                       &
!                + omyk2 *omyk2 *(zqdcor2 - zqdmod(i,j))
        omysc2 =  omysc(i,j) * omysc(i,j)
        omy (i,j,ivrs,ityp) = omy (i,j,ivrs,ityp) + omysc(i,j) *omykt
        om2 (i,j,ivrs,ityp) = om2 (i,j,ivrs,ityp) + omysc2     *omyk2t
        zwi (i,j,ivar,ityp) = zwi (i,j,ivar,ityp) + omysc2     *rhinct2        &
                                                  + omysc(i,j) *rhincrt *fmultw
      ENDDO
    ENDIF

!-------------------------------------------------------------------------------
!  Section 6: Spreading of observation increments of potential temperature
!-------------------------------------------------------------------------------

  ELSEIF (ivrs == 2) THEN
    tincr1 = c0
    tincr2 = c0
    IF (lsurf) THEN
      fftho  = EXP( - rdocp * LOG( zpobsu(ista) ) )
      IF (omyk1 >= epsy) tincr1 = xoisu(ista,1,ivar) * fftho
      IF (omyk2 >= epsy) tincr2 = xoisu(ista,2,ivar) * fftho
    ELSE
      fftho  = EXP( - rdocp * LOG( zpobua(ista) ) )
      IF (omyk1 >= epsy) tincr1 = xoiua(ista,1,ivar) * fftho
      IF (omyk2 >= epsy) tincr2 = xoiua(ista,2,ivar) * fftho
    ENDIF
    tincr1 = omyk1 * tincr1
    tincr2 = omyk2 * tincr2
    omykt  = omyk1 + omyk2
    tincrt = tincr1 + tincr2
    IF (.NOT. lti2) THEN
      tincr1 = omyk1 * tincr1
      tincr2 = omyk2 * tincr2
      omyk2t = omyk1 *omyk1 + omyk2 *omyk2
      tinct2 = tincr1 + tincr2
    ENDIF

    IF (lti2) THEN
!CDIR NODEP,VOVERTAKE,VOB
      DO ixc = 1 , ncutof(1)
        i = icutof(ixc,1)
        j = jcutof(ixc,1)
        omykm  =  omysc(i,j) * omykt
        omytt  =  omysc(i,j) *(omykm + fmultw) * zeklop(i,j)
        omy (i,j,ivrs,ityp) = omy (i,j,ivrs,ityp) + omykm
        om2 (i,j,ivrs,ityp) = om2 (i,j,ivrs,ityp) + omykm *omykm
        zwi (i,j,ivar,ityp) = zwi (i,j,ivar,ityp) + omytt *tincrt
      ENDDO
    ELSE
!CDIR NODEP,VOVERTAKE,VOB
      DO ixc = 1 , ncutof(1)
        i = icutof(ixc,1)
        j = jcutof(ixc,1)
        omysc2 =  omysc(i,j) * omysc(i,j)
        omytt  =  omysc(i,j) * zeklop(i,j)
        omy (i,j,ivrs,ityp) = omy (i,j,ivrs,ityp) + omysc(i,j) *omykt
        om2 (i,j,ivrs,ityp) = om2 (i,j,ivrs,ityp) + omysc2 *omyk2t
        zwi (i,j,ivar,ityp) = zwi (i,j,ivar,ityp) + omytt  *tincrt *fmultw     &
                                                  + omysc2 *tinct2 *zeklop(i,j)
      ENDDO
    ENDIF

!-------------------------------------------------------------------------------
!  Section 7: Spreading of observation increments of the horizontal wind vector
!-------------------------------------------------------------------------------

  ELSEIF (ivrs == 1) THEN
    ivay = ivar + 1
    uincr1 = c0
    vincr1 = c0
    uincr2 = c0
    vincr2 = c0
    IF (lsurf) THEN
      IF (omyk1 >= epsy) uincr1 = xoisu(ista,1,ivar)
      IF (omyk1 >= epsy) vincr1 = xoisu(ista,1,ivay)
      IF (omyk2 >= epsy) uincr2 = xoisu(ista,2,ivar)
      IF (omyk2 >= epsy) vincr2 = xoisu(ista,2,ivay)
    ELSE
      IF (omyk1 >= epsy) uincr1 = xoiua(ista,1,ivar)
      IF (omyk1 >= epsy) vincr1 = xoiua(ista,1,ivay)
      IF (omyk2 >= epsy) uincr2 = xoiua(ista,2,ivar)
      IF (omyk2 >= epsy) vincr2 = xoiua(ista,2,ivay)
    ENDIF
    uincr1 = omyk1 * uincr1
    vincr1 = omyk1 * vincr1
    uincr2 = omyk2 * uincr2
    vincr2 = omyk2 * vincr2
    omykt  = omyk1 + omyk2
    uincrt = uincr1 + uincr2
    vincrt = vincr1 + vincr2
    IF (.NOT. lti2) THEN
      uincr1 = omyk1 * uincr1
      vincr1 = omyk1 * vincr1
      uincr2 = omyk2 * uincr2
      vincr2 = omyk2 * vincr2
      omyk2t = omyk1 *omyk1 + omyk2 *omyk2
      uinct2 = uincr1 + uincr2
      vinct2 = vincr1 + vincr2
    ENDIF

    IF (lti2) THEN
!CDIR NODEP,VOVERTAKE,VOB
      DO ixc = 1 , ncutof(1)
        i = icutof(ixc,1)
        j = jcutof(ixc,1)
! weighted sum of spreaded incr. using 2-dim horizontal wind correlations
        zuoem  =  zcoruu(i,j) *uincrt  +  zcoruv(i,j) *vincrt
        zvoem  =  zcorvu(i,j) *uincrt  +  zcorvv(i,j) *vincrt
        omykm  =  omysc(i,j) * omykt
        omykh  =  omysc(i,j) *(omykm + fmultw)
        omy (i,j,ivrs,ityp) = omy (i,j,ivrs,ityp) + omykm
        om2 (i,j,ivrs,ityp) = om2 (i,j,ivrs,ityp) + omykm *omykm
        zwi (i,j,ivar,ityp) = zwi (i,j,ivar,ityp) + omykh *zuoem
        zwi (i,j,ivay,ityp) = zwi (i,j,ivay,ityp) + omykh *zvoem
      ENDDO
    ELSEIF (kwtyp(ityp) == 1) THEN
! (this 'elseif' (kwtyp == 1) has been introduced for speed-up)
!CDIR NODEP,VOVERTAKE,VOB
      DO ixc = 1 , ncutof(1)
        i = icutof(ixc,1)
        j = jcutof(ixc,1)
! weighted sum of spreaded incr. using 2-dim horizontal wind correlations
        zuoem2 =  zcoruu(i,j) *uinct2  +  zcoruv(i,j) *vinct2
        zvoem2 =  zcorvu(i,j) *uinct2  +  zcorvv(i,j) *vinct2
        zuoem  =  zcoruu(i,j) *uincrt  +  zcoruv(i,j) *vincrt
        zvoem  =  zcorvu(i,j) *uincrt  +  zcorvv(i,j) *vincrt
        omysc2 =  omysc(i,j) * omysc(i,j)
        omy (i,j,ivrs,ityp) = omy (i,j,ivrs,ityp) + omysc(i,j) *omykt
        om2 (i,j,ivrs,ityp) = om2 (i,j,ivrs,ityp) + omysc2     *omyk2t
        zwi (i,j,ivar,ityp) = zwi (i,j,ivar,ityp) + omysc2     *zuoem2
        zwi (i,j,ivay,ityp) = zwi (i,j,ivay,ityp) + omysc2     *zvoem2
      ENDDO
    ELSE
!CDIR NODEP,VOVERTAKE,VOB
      DO ixc = 1 , ncutof(1)
        i = icutof(ixc,1)
        j = jcutof(ixc,1)
! weighted sum of spreaded incr. using 2-dim horizontal wind correlations
        zuoem2 =  zcoruu(i,j) *uinct2  +  zcoruv(i,j) *vinct2
        zvoem2 =  zcorvu(i,j) *uinct2  +  zcorvv(i,j) *vinct2
        zuoem  =  zcoruu(i,j) *uincrt  +  zcoruv(i,j) *vincrt
        zvoem  =  zcorvu(i,j) *uincrt  +  zcorvv(i,j) *vincrt
        omysc2 =  omysc(i,j) * omysc(i,j)
        omy (i,j,ivrs,ityp) = omy (i,j,ivrs,ityp) + omysc(i,j) *omykt
        om2 (i,j,ivrs,ityp) = om2 (i,j,ivrs,ityp) + omysc2     *omyk2t
        zwi (i,j,ivar,ityp) = zwi (i,j,ivar,ityp) + omysc2     *zuoem2         &
                                                  + omysc(i,j) *zuoem
!                                                 + omysc(i,j) *zuoem * fmultw
        zwi (i,j,ivay,ityp) = zwi (i,j,ivay,ityp) + omysc2     *zvoem2         &
                                                  + omysc(i,j) *zvoem
!                                                 + omysc(i,j) *zvoem * fmultw
      ENDDO
    ENDIF
! ENDIF


! for geostrophic pressure increments balancing wind increments at the lowest
! model level derived from surface-level wind observations, the spreading of
! wind observation increments is re-done here, however with isotropic weighting
! for the lowest model level only. This corresponds to spreading along sigma
! surfaces.
! (At k == ke follows that: zsprobs = zsprsu(ista,k) = zspr(i,j,msi) ,
!  and the vertical weighting consists only of:
!                             omyvf = EXP( -(zpblsu(ista,k) /vpblsu(ivar))**2 )
!                     where: zpblsu = MAX( 0, theta(ke) - theta(2m) )
!                     where theta is the potential temperature)
! An additional weighting depends on the observation code type.

    IF ((mpsgcor >= 1) .AND. (k == ke) .AND. (lsurf)) THEN
      omykfty = c1
      ! no correction for wind increments from in-situ data if (mpsgcor == 1)
      IF ((kobtyp /= nscatt) .AND. (mpsgcor == 1))  omykfty = c0
      IF (lti2) THEN
!CDIR NODEP,VOVERTAKE,VOB
        DO ixc = 1 , ncutof(1)
          i = icutof(ixc,1)
          j = jcutof(ixc,1)
!         omysci = (c1+mscalar*zdds(i,j,ivrs))* EXP( -zdds(i,j,ivrs) ) * omykfty
          omysci =                              EXP( -zdds(i,j,ivrs) ) * omykfty
! weighted sum of spreaded incr. using 2-dim horizontal wind correlations
          zuoem  =  zcoruu(i,j) *uincrt  +  zcoruv(i,j) *vincrt
          zvoem  =  zcorvu(i,j) *uincrt  +  zcorvv(i,j) *vincrt
          omykm  =  omysci * omykt
          omykh  =  omysci *(omykm + fmultw)
          omy (i,j,5,ityp) = omy (i,j,5,ityp) + omykm
          om2 (i,j,4,ityp) = om2 (i,j,4,ityp) + omykm *omykm
          zwi (i,j,5,ityp) = zwi (i,j,5,ityp) + omykh *zuoem
          zwi (i,j,6,ityp) = zwi (i,j,6,ityp) + omykh *zvoem
        ENDDO
      ELSE
!CDIR NODEP,VOVERTAKE,VOB
        DO ixc = 1 , ncutof(1)
          i = icutof(ixc,1)
          j = jcutof(ixc,1)
          omysci =                              EXP( -zdds(i,j,ivrs) ) * omykfty
! weighted sum of spreaded incr. using 2-dim horizontal wind correlations
          zuoem2 =  zcoruu(i,j) *uinct2  +  zcoruv(i,j) *vinct2
          zvoem2 =  zcorvu(i,j) *uinct2  +  zcorvv(i,j) *vinct2
          zuoem  =  zcoruu(i,j) *uincrt  +  zcoruv(i,j) *vincrt
          zvoem  =  zcorvu(i,j) *uincrt  +  zcorvv(i,j) *vincrt
          omysc2 =  omysci * omysci
          omy (i,j,5,ityp) = omy (i,j,5,ityp) + omysci *omykt
          om2 (i,j,4,ityp) = om2 (i,j,4,ityp) + omysc2 *omyk2t
          zwi (i,j,5,ityp) = zwi (i,j,5,ityp) + omysc2 *zuoem2                 &
                                              + omysci *zuoem * fmultw
          zwi (i,j,6,ityp) = zwi (i,j,6,ityp) + omysc2 *zvoem2                 &
                                              + omysci *zvoem * fmultw
        ENDDO
      ENDIF
    ENDIF
  ENDIF
!-------------------------------------------------------------------------------
!  Section 8: Finalize
!-------------------------------------------------------------------------------

! IF ((lwonl) .AND. (k == ke) .AND. (ivrs == 2) .AND. (ntstep >= nstop-1)) THEN
!   IF (lupair) PRINT *,'zwt  ',zwtua(ista,1,ivrs), zwtua(ista,2,ivrs)
!   IF (lsurf ) PRINT *,'zwt  ',zwtsu(ista,1,ivrs), zwtsu(ista,2,ivrs)
!   PRINT *,'zwilv',zwi(ionl,jonl,ivar), omy(ionl,jonl,ivrs)
! ENDIF


  IF (lsigma) omymx = omymx * omyvf

  IF ((ntstep <= 1) .AND. (lwonl) .AND. (k == ke-1) .AND. (lupair)) THEN
    IF (      (io  <= iend) .AND. (io  >= istart)                              &
        .AND. (jo  <= jend) .AND. (jo  >= jstart)) THEN
      IF ((lcutof(io,jo,1)) .AND. (ivrs > 0)) omyvf = omyvf * omysc(io,jo)
      IF (.NOT. lcutof(io,jo,1)) omyvf = c0
      zdsobm = zsprobs - zsprua(ista,k)
      IF (modespr <= 1) THEN
        WRITE( nupr,'(''airmy: '',A,2I4,I2,'',omyk1/2'',2F6.3,'' vCS-fac''     &
                    &,2F6.3,'' dz'',F7.0,'' omym'',F9.6)' )   ystid, io, jo    &
             , ivar, omyk1, omyk2, fcorlua(ista,1,ivrs), fcorlua(ista,2,ivrs)  &
             , zdsobm, omyvf
      ELSE
        WRITE( nupr,'(''airmy: '',A,2I4,I2,'',omyk1/2'',2F6.3,'' vCS-fac''     &
                    &,2F6.3,'' dth'',F6.1,'' omym'',F9.6)' )   ystid, io, jo   &
             , ivar, omyk1, omyk2, fcorlua(ista,1,ivrs), fcorlua(ista,2,ivrs)  &
             , zdsobm, omyvf
      ENDIF
    ENDIF
  ENDIF

  IF ((ntstep <= 1) .AND. (lwonl) .AND. (ystid(1:3) == '066')) THEN
    WRITE( nupr,'(''omyk: '',A ,2I2,5F6.3,F8.2,L2,F5.2,2F8.1,F4.1)' )          &
           ystid, ivar, ivrs, omyk1, omyk2, omyvf, wtql1, wtql2, rdsprs        &
         , lcutof(ionl2,jonl2,1), zdds(ionl2,jonl2,ivrs)                       &
         , zspr  (ionl2,jonl2,msi ), zsprobs, zablfc(ionl2,jonl2,ivrs,iabl)

! valgrind complained about too small or uninitialized values
!US    WRITE( nupr,'(''omyk: '',A ,2I2,5F6.3,F8.2,L2,     2F8.1,F4.1)' )          &
!US           ystid, ivar, ivrs, omyk1, omyk2, omyvf, wtql1, wtql2, rdsprs        &
!US         , lcutof(ionl2,jonl2,1)                                               & !, zdds(ionl2,jonl2,ivrs)
!US         , zspr  (ionl2,jonl2,msi ), zsprobs, zablfc(ionl2,jonl2,ivrs,iabl)
  ENDIF

! IF ((ntstep <= 2) .AND. (lwobs) .AND. (lupair)) THEN
!   IF (((io == ionl) .AND. (jo == jonl)) .OR. (ista == 1)) THEN
!     ijstot = (iendspr - istaspr + 1) * (jendspr - jstaspr + 1)
!     WRITE( nupr,'(''air-gp: '',A ,3I4,'' level'',I3'': # of infl. g.p.:''    &
!                  ,2I5)' ) ystid, ivar, io_tot, jo_tot, k, ijstot, ncutof(1)
!   ENDIF
! ENDIF

!-------------------------------------------------------------------------------
! End of module procedure surf_upair_spread
!-------------------------------------------------------------------------------

END SUBROUTINE surf_upair_spread


!-------------------------------------------------------------------------------
!+ Module procedure to organize spreading of screen-level obs increments
!-------------------------------------------------------------------------------

SUBROUTINE surf_org_spread ( k )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure organizes the spreading of observation increments
!   from surface-level data after the processor has received all the required
!   'local' information from all observations.
!
! Method:
!   Setting up the loop over the surface stations.
!   Exclusion of all reports / data, for which no grid pts. of the current node
!   and vertical and temporal model level exist in the 4-dim area of influence.
!   Area of influence, spatial correlations, and spreading by call of other
!   procedures. 
!
! Written by        :  Christoph Schraff, DWD  (original version: 30.07.97)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER (KIND=iintegers), INTENT (IN)         ::       &
    k                   ! index of current vertical model level

! Local parameters:
! ----------------

  REAL (KIND = wp)         , PARAMETER  :: &
    thrmymx = 0.01_wp    ! threshold: if the largest spatial weight at layer
                             ! k for a given surface obs is below this threshold
                             ! then the obs is not used for spread. to layer k-1

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    istasu           ,& ! index of observing station
!   ista             ,& ! index of observing station
    ivrs             ,& ! index of spreaded quantity: 1=(u,v); 2=T; 3=RH
    ivar             ,& ! index of spreaded quantity: 1=(u,v); 3=T; 4=RH
    ityp             ,& ! index of weighted increment array
    irange , jrange  ,& ! zonal , meridional 'radii' of the domain which
                        ! contains the area of influence
    izmin  , jzmin   ,& ! coord. of lowest grid pt. within horiz. area of infl.
    ivex             ,& ! 'ivar', or (compressed) coord. (izmin,jzmin), or
                        ! 1 (0) if node has (no) grid pts. in area of influence
    nactio           ,& ! action executed in called proc. 'cutoff_wind_correl'
    kk     , kvifa   ,& ! loop index, upper limit for influenced model levels
    kaob             ,& ! mean of lowest and highest influenced model levels
    istat  , nstat      ! status of memory (de-)allocation

  REAL    (KIND=wp   )     ::  &
    acthr            ,& ! actual forecast hour
    wtvmx            ,& ! max. temporal nudging weight of obs. with index 'ista'
    zrcutof          ,& ! max. horizontal cut-off radius [km]
    fcnodiv          ,& ! weight to the non-divergence correction (wind correl.)
    tnodivs          ,& ! temporal factor to weight to non-divergence correction
    z2cut            ,& ! 2 * zvcutsu (zvcutsu: height at vertical cut-off)
    xincrp , yincrp  ,& ! components of the sum of the 2 wind obs. incr.
    omymx               ! max. spatial weight

  LOGICAL                  ::  &
    lsurf            ,& ! spreading of surface-level (not upper-air) obs. incr.
    larea            ,& ! first guess for area of influence partly in sub-domain
    lnearonl         ,& ! obs. location is near (ionl2,jonl2)
    lwobs            ,& ! (TRUE if obs. station (io,jo) lies in local domain):
                        !  replaced by lwobs = lwonl
    lprsys , lprthr     ! for printing for control

! Local (automatic) arrays:
! -------------------------

  LOGICAL                  ::  &
    lsprsv  (nsutot,3)   ! temporal and horiz. domain of node influenced by obs.

  INTEGER (KIND=iintegers) ::  &
    ijspr   (4,nsutot)   ! domain containing area of influence

  REAL    (KIND=wp   )     ::  &
    zrcutov (nsutot,3)   ! horizontal cut-off radius [km]

  REAL    (KIND=wp   )     , ALLOCATABLE , SAVE ::  &
    zmymx   (:,:)       ! max. spatial weight

  LOGICAL                  , ALLOCATABLE , SAVE ::  &
    lsprssu (:)         ! node has part of area of influence of station 'ista'
!
!------------ End of header ----------------------------------------------------


 
!-------------------------------------------------------------------------------
! Begin Subroutine surf_org_spread
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1: Exclude all surface-level reports / data, for which no grid pts.
!             of the current node and vertical and temporal model level exist
!             within the 4-dim area of influence.
!             (The temporal check implicitly includes quality control.)
!-------------------------------------------------------------------------------

  acthr  = ntstep * dtdeh
  lwobs  = .FALSE.
  lprsys = (lfirst) .OR. (ntstep <= 2)                                         &
                    .OR.        (rmod(  ntstep   *dt,c3600 ) < tconbox)
  lprthr = (lfirst) .OR. (      (rmod( (ntstep+1)*dt,c3600 ) < tconbox)        &
                          .AND. (REAL(ntstep,wp)*dtdeh <=                      &
                                 REAL(NINT(REAL(MIN(nudgend,nstop),wp)*dtdeh),wp)-0.9_wp))

  IF ((nsutot == 0) .AND. ((lfirst) .OR. (rmod( ntstep*dt,c3600 ) < tconbox))  &
                    .AND. (lcroot) .AND. (k == ke))                            &
    PRINT '(''hour '',F3.0,'' : no surface-level data'')' , acthr


! allocate memory of fields used only in current procedure
! --------------------------------------------------------

  IF (k == ke) THEN
    IF (ALLOCATED( lsprssu )) THEN
      PRINT '("CAUTION in src_sing_spread: lsprssu is already allocated "      &
            &,"at time ",I6)', ntstep
      DEALLOCATE ( lsprssu  , STAT=istat )
    ENDIF
    ALLOCATE ( lsprssu (nsutot  ) , STAT=istat )
    ALLOCATE ( zmymx   (nsutot,3) , STAT=istat )

! check if this sub-domain may be influenced by the current report (station)
! (The check includes horizontal check, temporal check,
!  and check of the nudging coefficients.)
! --------------------------------------------------------------------------

    DO istasu = 1 , nsutot
      ista = isortsu(istasu)
      io  = iosu (ista)
      jo  = josu (ista)
      lsprssu (ista) = .FALSE.
      zrcutov (istasu,1) = -c1
      zrcutov (istasu,2) = -c1
      zrcutov (istasu,3) = -c1
      IF (      ((zwtsu(ista,1,1) > epsy) .OR. (zwtsu(ista,2,1) > epsy))       &
          .AND. (gnudgsu(1) > epsy))                                           &
        zrcutov (istasu,1) = zriflsu(ista,1) * cutofsu(1)
      IF (      ((zwtsu(ista,1,2) > epsy) .OR. (zwtsu(ista,2,2) > epsy))       &
          .AND. (gnudgsu(3) > epsy))                                           &
        zrcutov (istasu,2) = zriflsu(ista,2) * cutofsu(3)
      IF (      ((zwtsu(ista,1,3) > epsy) .OR. (zwtsu(ista,2,3) > epsy))       &
          .AND. (gnudgsu(4) > epsy))                                           &
        zrcutov (istasu,3) = zriflsu(ista,3) * cutofsu(4)
      zrcutof = MAX( zrcutov(istasu,1) , zrcutov(istasu,2) , zrcutov(istasu,3) )
      jrange  = INT( zrcutof * gppkmj )
      jstaspr = MAX( INT( jo - jrange ,iintegers) , jrs_tot )
      jendspr = MIN( INT( jo + jrange ,iintegers) , jre_tot )
      irange  = INT( zrcutof * MAX( gppkmi(jstaspr) , gppkmi(jendspr) ) )
!     irange  = INT( zrcutof * MAX( gppkmi(MIN( jstaspr,jend   ))              &
!                                 , gppkmi(MAX( jendspr,jstart )) ) )
      jstaspr = MAX( INT( jo - jrange ,iintegers) , jstart )
      jendspr = MIN( INT( jo + jrange ,iintegers) , jend )
      istaspr = MAX( INT( io - irange ,iintegers) , istart )
      iendspr = MIN( INT( io + irange ,iintegers) , iend )
!     lsprssu (ista) = (      (jstaspr <= jend) .AND. (jendspr >= jstart)      &
!                       .AND. (istaspr <= iend) .AND. (iendspr >= istart)      &
!                       .AND. (zrcutof > epsy) )
      larea             = (      (jstaspr <= jend) .AND. (jendspr >= jstart)   &
                           .AND. (istaspr <= iend) .AND. (iendspr >= istart)   &
                           .AND. (zrcutof > epsy) )
      lsprsv (istasu,1) = (larea) .AND. (zrcutov(istasu,1) > epsy)
      lsprsv (istasu,2) = (larea) .AND. (zrcutov(istasu,2) > epsy)
      lsprsv (istasu,3) = (larea) .AND. (zrcutov(istasu,3) > epsy)
    ENDDO

! check for which variables this sub-domain may
! be influenced by the current report (station)
! ---------------------------------------------

    nactio = 6
    DO ivrs = 1 , 3
      ivex = ivrs + MIN( ivrs-1 , 1 )
      DO istasu = 1 , nsutot
        IF (lsprsv(istasu,ivrs)) THEN
          ista = isortsu(istasu)
          io  = iosu (ista)
          jo  = josu (ista)
          jrange  = INT( zrcutov(istasu,ivrs) * gppkmj )
          jstaspr = MAX( INT( jo - jrange ,iintegers) , jrs_tot )
          jendspr = MIN( INT( jo + jrange ,iintegers) , jre_tot )
          irange  = INT( zrcutov(istasu,ivrs) *MAX( gppkmi(jstaspr)            &
                                                  , gppkmi(jendspr) ))
!         irange  = INT( zrcutov(istasu,ivrs) *MAX( gppkmi(MIN( jstaspr,jend)) &
!                                             , gppkmi(MAX( jendspr,jstart)) ) )
          ijspr (3,istasu) = MAX( INT( jo - jrange ,iintegers) , jstart )
          ijspr (4,istasu) = MIN( INT( jo + jrange ,iintegers) , jend )
          ijspr (1,istasu) = MAX( INT( io - irange ,iintegers) , istart )
          ijspr (2,istasu) = MIN( INT( io + irange ,iintegers) , iend )
          lsprsv (istasu,ivrs)  = (      (ijspr(3,istasu) <= jend)             &
                                   .AND. (ijspr(4,istasu) >= jstart)           &
                                   .AND. (ijspr(1,istasu) <= iend)             &
                                   .AND. (ijspr(2,istasu) >= istart) )
        ENDIF
      ENDDO
      DO istasu = 1 , nsutot
        ista = isortsu(istasu)
        IF (lsprsv(istasu,ivrs)) THEN
          io  = iosu (ista)
          jo  = josu (ista)
          istaspr = ijspr(1,istasu)
          iendspr = ijspr(2,istasu)
          jstaspr = ijspr(3,istasu)
          jendspr = ijspr(4,istasu)
          rtinfl  = zrcutov(istasu,ivrs)
!         nactio  = 6      ! set further above

          CALL cutoff_wind_correl ( nactio , k , c0 , ivex )
!         =======================

          lsprsv (istasu,ivrs) = (ivex > 0)
          IF (lsprsv(istasu,ivrs))  lsprssu (ista) = .TRUE.
        ENDIF
        IF (.NOT. lsprsv(istasu,ivrs))  kviflsu (ista,ivrs) = ke + 1
      ENDDO
    ENDDO

!-------------------------------------------------------------------------------
!  Section 2: Get the vertical range of influenced model levels
!             (if spreading is horizontal)
!-------------------------------------------------------------------------------

    DO istasu = 1 , nsutot

      ista = isortsu(istasu)
      io  = iosu (ista)
      jo  = josu (ista)

      IF (lsprssu(ista)) THEN
        ystid  = ystidsu (ista) (1:ilstidp)
        lwobs  = lwonl
        lnearonl = (ABS(io-ionl2) <= 3) .AND. (ABS(jo-jonl2) <= 3)

        DO ivrs = 1, 3
          ivar = ivrs + MIN( ivrs-1 , 1 )
!         IF ((lwonl) .AND. (ivrs == 2) .AND. (ntstep >= nstop-1))             &
!           PRINT *,'kvis ',io, jo, '  ', ystid, '  '                          &
!                          ,kviflsu(ista,ivrs),zvcutsu(ista,ivrs)
          IF ((msprpsu == 1) .AND. (kviflsu(ista,ivrs) < ke)) THEN
            zrcutof = zriflsu(ista,ivrs) * cutofsu(ivar)
            jrange  = INT( zrcutof * gppkmj )
            jstaspr = MAX( INT( jo - jrange ,iintegers) , jrs_tot )
            jendspr = MIN( INT( jo + jrange ,iintegers) , jre_tot )
            irange  = INT( zrcutof * MAX( gppkmi(jstaspr) , gppkmi(jendspr) ) )
!           irange  = INT( zrcutof * MAX( gppkmi(MIN( jstaspr,jend   ))        &
!                                       , gppkmi(MAX( jendspr,jstart )) ) )
            jstaspr = MAX( INT( jo - jrange ,iintegers) , jstart )
            jendspr = MIN( INT( jo + jrange ,iintegers) , jend )
            istaspr = MAX( INT( io - irange ,iintegers) , istart )
            iendspr = MIN( INT( io + irange ,iintegers) , iend )
            izmin = 0
            jzmin = 0
            IF (      (jstaspr <= jend) .AND. (jendspr >= jstart)              &
                .AND. (istaspr <= iend) .AND. (iendspr >= istart) ) THEN
              rtinfl = zrcutof
              nactio = 5
              ivex   = ivar

              CALL cutoff_wind_correl ( nactio , k , c0 , ivex )
!             =======================

              izmin = MOD( ivex , 10000 )
              jzmin = ivex / 10000
            ENDIF
            IF (izmin == 0) THEN
              kviflsu (ista,ivrs) = ke + 1
            ELSEIF (msprpsu >= 1) THEN
              z2cut = c2 * zvcutsu(ista,ivrs)
              kvifa = ABS( kviflsu(ista,ivrs) )
              cutsu: DO kk = kvifa , ke + 1
                IF (kk == ke+1)                                       EXIT cutsu
                IF ( hhl(izmin,jzmin,kk)                                       &
                    +hhl(izmin,jzmin,kk+1) <= z2cut)                  EXIT cutsu
              ENDDO cutsu
              kviflsu (ista,ivrs) = kk
            ENDIF
          ENDIF
!         IF ((lwonl) .AND. (ivrs == 2) .AND. (ntstep >= nstop-1))             &
!           PRINT *,'kvis ',io, jo, '  ', ystid, '  ', kviflsu(ista,ivrs)
          IF ((lfirst) .AND. (lwobs) .AND. (lnearonl))                         &
            WRITE( nupr,'(''sfc-kl: '',A ,I2,F7.0,'', height at obs/cutoff:''  &
                        &,2F7.1,'', k-range'',I3)' ) ystid, ivrs, zpobsu(ista) &
                 , zsposu(ista), zvcutsu(ista,ivrs), kviflsu(ista,ivrs)
          IF (k == ke) zmymx (ista,ivrs) = c1

        ENDDO
      ENDIF
    ENDDO

  ENDIF   !  (k == ke)


!-------------------------------------------------------------------------------
!  Section 3: Organize the spreading of observation increments
!-------------------------------------------------------------------------------

! start loop over surface stations
! ----------------------------------

loop_over_surface_stations:  DO istasu = 1 , nsutot

  ista = isortsu(istasu)
  io  = iosu (ista)
  jo  = josu (ista)

  IF (lsprssu(ista)) THEN

    io_tot = iosu_tot(ista)
    jo_tot = josu_tot(ista)
    ystid  = ystidsu (ista) (1:ilstidp)

!   lwobs  =       (io  <= iend) .AND. (io  >= istart)                         &
!            .AND. (jo  <= jend) .AND. (jo  >= jstart)
    lwobs  = lwonl
    lnearonl = (ABS(io-ionl2) <= 3) .AND. (ABS(jo-jonl2) <= 3)

! ! NO loop over time index
! -------------------------
!   loop_over_time_index:  DO itim = 1, 2

! start loop over variables
! -------------------------

    loop_over_variables:  DO ivrs = 1, 3

      ivar = ivrs + MIN( ivrs-1 , 1 )

! check if vertical and temporal level may be influenced by current report
! ------------------------------------------------------------------------

!     IF ((k >= kviflsu(ista,ivrs)) .AND. (zmymx(ista,ivrs) > thrmymx)) THEN
      IF (k >= kviflsu(ista,ivrs)) THEN

! get the scaled horizontal distance to the station location, the geometrical
! factors of the 2-dim horizontal wind correlations, and the area of influence
! ----------------------------------------------------------------------------

        wtvmx = MAX( zwtsu(ista,1,ivrs) , zwtsu(ista,2,ivrs) )
        fcnodiv = c0
        IF (ivar == 1) THEN
          tnodivs = tnondiv * MIN( wtuksua / wtukrsa , c1 )
          fcnodiv = MIN( c1 ,  (c1 + (tnodivs-c1) *(c1-wtvmx))                 &
                             * (cnondiv + fnondiv *fnodivc(1)) )
        ENDIF
        rtinfl  = MAX( epsy , (c1 + (rhtfsu(ivar) -c1) *(c1- wtvmx))           &
                             *rhiflsu(ivar) )
        zrcutof = cutofsu(ivar) * rtinfl
        jrange  = INT( zrcutof * gppkmj )
        jstaspr = MAX( INT( jo - jrange ,iintegers) , jrs_tot )
        jendspr = MIN( INT( jo + jrange ,iintegers) , jre_tot )
        irange  = INT( zrcutof * MAX( gppkmi(jstaspr) , gppkmi(jendspr) ) )
!       irange  = INT( zrcutof * MAX( gppkmi(MIN( jstaspr,jend   ))            &
!                                   , gppkmi(MAX( jendspr,jstart )) ) )
        jstaspr = MAX( INT( jo - jrange ,iintegers) , jstart )
        jendspr = MIN( INT( jo + jrange ,iintegers) , jend )
        istaspr = MAX( INT( io - irange ,iintegers) , istart )
        iendspr = MIN( INT( io + irange ,iintegers) , iend )
        nactio  = 3
        IF (     (jstaspr > jend) .OR. (jendspr < jstart)                      &
            .OR. (istaspr > iend) .OR. (iendspr < istart) )  nactio  =  0

! IF ((lwonl) .AND. (k == ke) .AND. (ivrs == 2) .AND. (ntstep >= nstop-1)) THEN
!   PRINT *,'ispr ',io, jo, '  ', ystid, '  ', iend, jend, ivrs, zrcutof       &
!                  ,irange, jrange, istaspr, jstaspr, iendspr, jendspr
! ENDIF


        CALL cutoff_wind_correl ( nactio , k , fcnodiv , ivar )
!       =======================

        omymx  = c0

        lsurf  = .TRUE.

! spreading of observation increments
! -----------------------------------

        IF ((k >= kviflsu(ista,ivrs)) .AND. (nactio > 0))                      &

          CALL surf_upair_spread ( lsurf , ivar , k , omymx )
!         ======================

! printout for control
! --------------------

        IF ((lprsys) .AND. (lwobs)) THEN
          kaob = (kviflsu(ista,ivrs) + ke) / 2
          ityp = isetyp0
          DO kk = SUM( niwtyp ) , 1 , -1
            IF (   ((iwtyp(kk) > 0) .AND. ( iwtyp(kk) == kobtysu(ista)))       &
               .OR.((iwtyp(kk) < 0) .AND. (-iwtyp(kk) == kcdtysu(ista))))      &
              ityp = isetyp(kk)
          ENDDO
          IF ((k == kaob) .AND. (lnearonl) .AND. (ivar == 4))                  &
            WRITE( nupr,'(''sfc-infl(km)'',F5.0, 1X,A , 2I4, F7.0              &
                        &,'', zqwi,omyq,2'',F9.6,2F8.4)' )                     &
                   zrcutof, ystid, io_tot, jo_tot, zpobsu(ista)                &
                 , zwi(ionl2,jonl2,4,ityp), omy(ionl2,jonl2,3,ityp)            &
                 , om2(ionl2,jonl2,3,ityp)
        ENDIF
        IF ((lprthr) .AND. (lwobs)) THEN
!         IF (      ((k <= kviflsu(ista,ivrs)) .OR. (omymx <= thrmymx))        &
          IF (      (k <= kviflsu(ista,ivrs))                                  &
              .AND. ((lfirst) .OR. (k < ke) .OR. (ivar == 4))) THEN
            IF (ivar == 1) THEN
              xincrp = xoisu(ista,1,ivar) + xoisu(ista,2,ivar)
              yincrp = xoisu(ista,1,2   ) + xoisu(ista,2,2   )
              WRITE( nupr,'(''thresh uv k'',2I3   ,2I4,1X, A , F6.0            &
                          &,'',mx omy:k+1,k:'',2F6.3,'',uv incr'',2F6.1)' )    &
                     k, kviflsu(ista,ivrs), io_tot,jo_tot, ystid, zzobsu(ista) &
                   , zmymx(ista,ivrs), omymx, xincrp, yincrp
            ELSEIF (ivar == 3) THEN
              WRITE( nupr,'(''thresh T  k'',2I3   ,2I4,1X, A , F6.0            &
                          &,'',mx omy:k+1,k:'',2F6.3,'', T incr'',2F6.1)' )    &
                     k, kviflsu(ista,ivrs), io_tot,jo_tot, ystid, zzobsu(ista) &
                   , zmymx(ista,ivrs), omymx, xoisu(ista,1,ivar)               &
                                            , xoisu(ista,2,ivar)
            ELSEIF (ivar == 4) THEN
              WRITE( nupr,'(''thresh RH k'',2I3   ,2I4,1X, A , F6.0            &
                          &,'',mx omy:k+1,k:'',2F6.3,'',RH incr'',2F6.3)' )    &
                     k, kviflsu(ista,ivrs), io_tot,jo_tot, ystid, zzobsu(ista) &
                   , zmymx(ista,ivrs), omymx, xoisu(ista,1,ivar)               &
                                            , xoisu(ista,2,ivar)
            ENDIF
          ENDIF
        ENDIF

        IF (k >= kviflsu(ista,ivrs)) zmymx (ista,ivrs) = omymx

      ENDIF

    ENDDO loop_over_variables

!   ENDDO loop_over_time_index

  ENDIF ! lsprssu(ista)

ENDDO loop_over_surface_stations


! deallocate memory
! -----------------

  IF (k == 1) THEN
    DEALLOCATE ( lsprssu , STAT=istat )
    DEALLOCATE ( zmymx   , STAT=istat )
  ENDIF

  IF (ldump_ascii) THEN
    ! flush YUPRINT file
    IF ((lwobs) .AND. ((lprthr) .OR. (lprsys))                                 &
                .AND. ((k == ke) .OR. (k == 1) .OR. (k >= ke-6))) THEN
      CLOSE (nupr)
      OPEN  (nupr,FILE=yuprint,FORM='FORMATTED',STATUS='OLD'                   &
                              ,POSITION='APPEND',IOSTAT=nstat)
      IF (nstat /= 0) PRINT '("OPENING OF FILE yuprint FAILED, surf_org_spread")'
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure surf_org_spread
!-------------------------------------------------------------------------------

END SUBROUTINE surf_org_spread


!-------------------------------------------------------------------------------
!+ Module procedure to organize spreading of upper-air single-level obs increm.
!-------------------------------------------------------------------------------

SUBROUTINE upair_org_spread ( k )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure organizes the spreading of observation increments
!   from upper-air single-level data after the processor has received all the
!   required 'local' information from all observations.
!
! Method:
!   Setting up the loop over the upper-air single-level stations.
!   Exclusion of all reports / data, for which no grid pts. of the current node
!   and vertical and temporal model level exist in the 4-dim area of influence.
!   Area of influence, spatial correlations, and spreading by call of other
!   procedures. 
!
! Written by        :  Christoph Schraff, DWD  (original version: 16.07.97)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER (KIND=iintegers), INTENT (IN)         ::       &
    k                   ! index of current vertical model level

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    istaua           ,& ! index of observing station
!   ista             ,& ! index of observing station
    ivrs             ,& ! index of spreaded quantity: 1=(u,v); 2=T; 3=RH
    ivar             ,& ! index of spreaded quantity: 1=(u,v); 3=T; 4=RH
    ityp             ,& ! index of weighted increment array
    kcdtyp           ,& ! CMA observation code type
    irange , jrange  ,& ! zonal , meridional 'radii' of the domain which
                        ! contains the area of influence
    izext  , jzext   ,& ! coord. of lowest or highest grid pt. within infl. area
    ivex             ,& ! 'ivar', or (compressed) coord. (izext,jzext), or
                        ! 1 (0) if node has (no) grid pts. in area of influence
    ilva   , ilvb    ,& ! indices of correlation scale table levels that are
                        ! adjacent to the specified model level
    nactio           ,& ! action executed in called proc. 'cutoff_wind_correl'
    nba    , kk      ,& ! loop indices
    kaob             ,& ! mean of lowest and highest influenced model levels
    istat  , nstat      ! status of memory (de-)allocation

  REAL    (KIND=wp   )     ::  &
    acthr            ,& ! actual forecast hour
    wtvmx            ,& ! max. temporal nudging weight of obs. with index 'ista'
    zrcutof          ,& ! max. horizontal cut-off radius [km]
    fvnodiv          ,& ! vertically dependent part of weight to non-div. corr.
    rvinfl           ,& ! vertically dependent part of horizontal correl scale
    tnodivs          ,& ! temporal factor to weight to non-divergence correction
    z2cut            ,& ! 2 * zvcutua (zvcutua: height at vertical cut-off)
    zlopf            ,& ! weight factor for vertical interpolations from the
                        ! correlation scale table
    xincrp , yincrp  ,& ! components of the sum of the 2 wind obs. incr.
    omymx               ! max. spatial weight

  LOGICAL                  ::  &
    lsurf            ,& ! spreading of surface-level (not upper-air) obs. incr.
    larea            ,& ! first guess for area of influence partly in sub-domain
    lobtyp           ,& ! unknown single-level upper-air observation type
    lnearonl         ,& ! obs. location is near (ionl2,jonl2)
    lwobs            ,& ! (TRUE if obs. station (io,jo) lies in local domain):
                        !  replaced by lwobs = lwonl
    lprsys , lprthr     ! for printing for control

! Local (automatic) arrays:
! -------------------------

  LOGICAL                  ::  &
    lsprsv  (nuatot,3)   ! temporal and horiz. domain of node influenced by obs.

  INTEGER (KIND=iintegers) ::  &
    ijspr   (4,nuatot),& ! domain containing area of influence
    ilv     (  nuatot)   ! nearest standard level above spreading point

  REAL    (KIND=wp   )     ::  &
    zrcutov (nuatot,3),& ! horizontal cut-off radius [km]
    fcnodiv (nuatot,3),& ! non-divergence correction factor for 2-D wind correl.
    rhtinfl (nuatot,3)   ! horizontal correlation scale

  LOGICAL                  , ALLOCATABLE , SAVE ::  &
    lsprsua (:)          ! node has part of area of influence of station 'ista'
!
!------------ End of header ----------------------------------------------------


 
!-------------------------------------------------------------------------------
! Begin Subroutine upair_org_spread
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1: Exclude all upper-air reports / data, for which no grid pts. of
!             the current node and vertical and temporal model level exist
!             within the 4-dim area of influence.
!             (Since upper limits for horizontal radii of influence are used,
!              these checks will find most, but not all data that have no in-
!              fluence. The temporal check implicitly includes quality control.)
!-------------------------------------------------------------------------------

  acthr  = ntstep * dtdeh
  lwobs  = .FALSE.   ! US initialize it somehow
  lprsys = (lfirst) .OR. (ntstep <= 2)                                         &
                    .OR.        (rmod(  ntstep   *dt,c3600 ) < tconbox)
  lprthr = (lfirst) .OR. (      (rmod( (ntstep+1)*dt,c3600 ) < tconbox)        &
                          .AND. (REAL(ntstep,wp)*dtdeh <=                      &
                                 REAL(NINT(REAL(MIN(nudgend,nstop),wp)*dtdeh),wp)-0.9_wp))

  IF ((nuatot == 0) .AND. ((lfirst) .OR. (rmod( ntstep*dt,c3600 ) < tconbox))  &
                    .AND. (lcroot) .AND. (k == ke))                            &
    PRINT '(''hour '',F3.0,'' : no upper-air single-level data'')' , acthr


! allocate memory of fields used only in current procedure
! --------------------------------------------------------

  IF (k == ke) THEN

    IF (ALLOCATED( lsprsua )) THEN
      PRINT '("CAUTION in src_sing_spread: lsprsua is already allocated "      &
            &,"at time ",I6)', ntstep
      DEALLOCATE ( lsprsua  , STAT=istat )
    ENDIF
    ALLOCATE ( lsprsua (nuatot) , STAT=istat )

    lobtyp = .FALSE.

! check if this sub-domain may be influenced by the current report (station)
! (The check includes horizontal check, temporal check,
!  and check of the nudging coefficients.)
! --------------------------------------------------------------------------

    DO istaua = 1 , nuatot
      ista = isortua(istaua)
      io  = ioua (ista)
      jo  = joua (ista)
      kcdtyp = kcdtyua(ista)
      lsprsua (ista) = .FALSE.
      zrcutov (istaua,1) = -c1
      zrcutov (istaua,2) = -c1
      zrcutov (istaua,3) = -c1
      IF (      ((zwtua(ista,1,1) > epsy) .OR. (zwtua(ista,2,1) > epsy))       &
          .AND. (     ((kcdtyp /= nmodes) .AND. (gnudgar(1) > epsy))           &
                 .OR. ((kcdtyp == nmodes) .AND. (gnudgms(1) > epsy))))         &
        zrcutov (istaua,1) = zriflua(ista,2,1) * cutofr(1)
      IF (      ((zwtua(ista,1,2) > epsy) .OR. (zwtua(ista,2,2) > epsy))       &
          .AND. (     ((kcdtyp /= nmodes) .AND. (gnudgar(3) > epsy))           &
                 .OR. ((kcdtyp == nmodes) .AND. (gnudgms(3) > epsy))))         &
        zrcutov (istaua,2) = zriflua(ista,2,2) * cutofr(3)
      IF (      ((zwtua(ista,1,3) > epsy) .OR. (zwtua(ista,2,3) > epsy))       &
          .AND. (     ((kcdtyp /= nmodes) .AND. (gnudgar(4) > epsy))           &
                 .OR. ((kcdtyp == nmodes) .AND. (gnudgms(4) > epsy))))         &
        zrcutov (istaua,3) = zriflua(ista,2,3) * cutofr(4)
      zrcutof = MAX( zrcutov(istaua,1) , zrcutov(istaua,2) , zrcutov(istaua,3) )
      jrange  = INT( zrcutof * gppkmj )
      jstaspr = MAX( INT( jo - jrange ,iintegers) , jrs_tot )
      jendspr = MIN( INT( jo + jrange ,iintegers) , jre_tot )
      irange  = INT( zrcutof * MAX( gppkmi(jstaspr) , gppkmi(jendspr) ) )
!     irange  = INT( zrcutof * MAX( gppkmi(MIN( jstaspr,jend   ))              &
!                                 , gppkmi(MAX( jendspr,jstart )) ) )
      jstaspr = MAX( INT( jo - jrange ,iintegers) , jstart )
      jendspr = MIN( INT( jo + jrange ,iintegers) , jend )
      istaspr = MAX( INT( io - irange ,iintegers) , istart )
      iendspr = MIN( INT( io + irange ,iintegers) , iend )
!     lsprsua (ista) = (      (jstaspr <= jend) .AND. (jendspr >= jstart)      &
!                       .AND. (istaspr <= iend) .AND. (iendspr >= istart)      &
!                       .AND. (zrcutof > epsy) )
      larea             = (      (jstaspr <= jend) .AND. (jendspr >= jstart)   &
                           .AND. (istaspr <= iend) .AND. (iendspr >= istart)   &
                           .AND. (zrcutof > epsy) )
      lsprsv (istaua,1) = (larea) .AND. (zrcutov(istaua,1) > epsy)
      lsprsv (istaua,2) = (larea) .AND. (zrcutov(istaua,2) > epsy)
      lsprsv (istaua,3) = (larea) .AND. (zrcutov(istaua,3) > epsy)
      IF (      (kobtyua(ista) /= nairep) .AND. (lcroot)                       &
          .AND. (kobtyua(ista) /= nsatob) .AND. (lprsys))  lobtyp = .TRUE.
    ENDDO
    IF (lobtyp)  PRINT *,'upair_local_info: obs type not AIREP or SATOB !!!'

! check for which variables this sub-domain may
! be influenced by the current report (station)
! ---------------------------------------------

    nactio = 6
    DO ivrs = 1 , 3
      ivex = ivrs + MIN( ivrs-1 , 1 )
      DO istaua = 1 , nuatot
!       IF (ivrs == 1)  lsprsua (ista) = .FALSE.
!       IF (lsprsv(istaua,ivrs))                                               &
!         lsprsv (istaua,ivrs) = (zrcutov(istaua,ivrs) > epsy)
        IF (lsprsv(istaua,ivrs)) THEN
          ista = isortua(istaua)
          io  = ioua (ista)
          jo  = joua (ista)
          jrange  = INT( zrcutov(istaua,ivrs) * gppkmj )
          jstaspr = MAX( INT( jo - jrange ,iintegers) , jrs_tot )
          jendspr = MIN( INT( jo + jrange ,iintegers) , jre_tot )
          irange  = INT( zrcutov(istaua,ivrs) *MAX( gppkmi(jstaspr)            &
                                                  , gppkmi(jendspr) ))
!         irange  = INT( zrcutov(istaua,ivrs) *MAX( gppkmi(MIN( jstaspr,jend)) &
!                                             , gppkmi(MAX( jendspr,jstart)) ) )
          ijspr (3,istaua) = MAX( INT( jo - jrange ,iintegers) , jstart )
          ijspr (4,istaua) = MIN( INT( jo + jrange ,iintegers) , jend )
          ijspr (1,istaua) = MAX( INT( io - irange ,iintegers) , istart )
          ijspr (2,istaua) = MIN( INT( io + irange ,iintegers) , iend )
          lsprsv (istaua,ivrs)  = (      (ijspr(3,istaua) <= jend)             &
                                   .AND. (ijspr(4,istaua) >= jstart)           &
                                   .AND. (ijspr(1,istaua) <= iend)             &
                                   .AND. (ijspr(2,istaua) >= istart) )
        ENDIF
      ENDDO
      DO istaua = 1 , nuatot
        IF (lsprsv(istaua,ivrs)) THEN
          ista = isortua(istaua)
          io  = ioua (ista)
          jo  = joua (ista)
          istaspr = ijspr(1,istaua)
          iendspr = ijspr(2,istaua)
          jstaspr = ijspr(3,istaua)
          jendspr = ijspr(4,istaua)
          rtinfl  = zrcutov(istaua,ivrs)
!         nactio  = 6      ! set further above

          CALL cutoff_wind_correl ( nactio , k , c0 , ivex )
!         =======================

          lsprsv (istaua,ivrs) = (ivex > 0)
          IF (lsprsv(istaua,ivrs))  lsprsua (ista) = .TRUE.
        ENDIF
      ENDDO
      DO istaua = 1 , nuatot
        IF (.NOT. lsprsv(istaua,ivrs)) THEN
          ista = isortua(istaua)
          kviflua (ista,1,ivrs) = 0
          kviflua (ista,2,ivrs) = ke + 1
        ENDIF
      ENDDO
    ENDDO

!-------------------------------------------------------------------------------
!  Section 2: Get the vertical range of influenced model levels
!             (if spreading is horizontal)
!-------------------------------------------------------------------------------

    DO istaua = 1 , nuatot

      ista = isortua(istaua)
      io  = ioua (ista)
      jo  = joua (ista)

      IF (lsprsua(ista)) THEN
        ystid  = ystidua (ista) (1:ilstidp)
        lwobs  = lwonl
        lnearonl = (ABS(io-ionl2) <= 5) .AND. (ABS(jo-jonl2) <= 5)

        DO ivrs = 1, 3
          ivar = ivrs + MIN( ivrs-1 , 1 )
!         IF ((lwonl) .AND. (ivrs == 2) .AND. (ntstep >= nstop-1))             &
!           PRINT *,'kvia ',io, jo, '  ', ystid, '  '                          &
!                          ,kviflua(ista,1,ivrs), kviflua(ista,2,ivrs)         &
!                          ,zvcutua(ista,1,ivrs), zvcutua(ista,2,ivrs)
          DO nba = 2 , 1 , -1
            IF (kviflua(ista,nba,ivrs) < 0) THEN
              zrcutof = zriflua(ista,nba,ivrs) * cutofr(ivar)
              jrange  = INT( zrcutof * gppkmj )
              jstaspr = MAX( INT( jo - jrange ,iintegers) , jrs_tot )
              jendspr = MIN( INT( jo + jrange ,iintegers) , jre_tot )
              irange  = INT( zrcutof * MAX( gppkmi(jstaspr), gppkmi(jendspr) ) )
!             irange  = INT( zrcutof * MAX( gppkmi(MIN( jstaspr,jend   ))      &
!                                         , gppkmi(MAX( jendspr,jstart )) ) )
              jstaspr = MAX( INT( jo - jrange ,iintegers) , jstart )
              jendspr = MIN( INT( jo + jrange ,iintegers) , jend )
              istaspr = MAX( INT( io - irange ,iintegers) , istart )
              iendspr = MIN( INT( io + irange ,iintegers) , iend )
              izext = 0
              jzext = 0
              IF (      (jstaspr <= jend) .AND. (jendspr >= jstart)            &
                  .AND. (istaspr <= iend) .AND. (iendspr >= istart) ) THEN
                rtinfl = zrcutof
                nactio = 3 + nba
                ivex   = ivar

                CALL cutoff_wind_correl ( nactio , k , c0 , ivex )
!               =======================

                izext = MOD( ivex , 10000 )
                jzext = ivex / 10000
              ENDIF
              IF ((izext == 0) .AND. (nba == 2)) THEN
                kviflua (ista,1,ivrs) = 0
                kviflua (ista,2,ivrs) = ke + 1
              ELSEIF (izext > 0) THEN
                z2cut = c2 * zvcutua(ista,nba,ivrs)
                cut: DO kk = 1 , ke + 1
                  kviflua (ista,nba,ivrs) = kk
                  IF (kk == ke+1)                                       EXIT cut
                  IF ( hhl(izext,jzext,kk)                                     &
                      +hhl(izext,jzext,kk+1) <= z2cut)                  EXIT cut
                ENDDO cut
                IF (nba == 1) kviflua(ista,nba,ivrs) = kviflua(ista,nba,ivrs) -1
!               IF ((lwonl) .AND. (ivrs == 2) .AND. (ntstep >= nstop-1))       &
!                 PRINT *,'hhext',izext,jzext, '  ', ystid, '  '               &
!                                ,hhl(izext,jzext,kviflua(ista,nba,ivrs)-1)    &
!                                ,hhl(izext,jzext,kviflua(ista,nba,ivrs))      &
!                                ,hhl(izext,jzext,kviflua(ista,nba,ivrs)+1)
!               IF ((lwonl) .AND. (ivrs == 2) .AND. (ntstep >= nstop-1))       &
!                 PRINT *,'hhonl',ionl ,jonl , '  ', ystid, '  '               &
!                                ,hhl(ionl ,jonl ,kviflua(ista,nba,ivrs)-1)    &
!                                ,hhl(ionl ,jonl ,kviflua(ista,nba,ivrs))      &
!                                ,hhl(ionl ,jonl ,kviflua(ista,nba,ivrs)+1)
              ELSE
                kviflua (ista,nba,ivrs) = ABS( kviflua(ista,nba,ivrs) ) - 1
              ENDIF
            ENDIF
          ENDDO
!         IF ((lwonl) .AND. (ivrs == 2) .AND. (ntstep >= nstop-1))             &
!           PRINT *,'kvia ',io, jo, '  ', ystid, '  '                          &
!                          ,kviflua(ista,1,ivrs), kviflua(ista,2,ivrs)
          IF ((lfirst) .AND. (lwobs) .AND. (lnearonl))                         &
            WRITE( nupr,'(''air-kl: '',A ,I2,F7.0,'', height at obs/cutoff:''  &
                        &,3F7.0,'', k-range'',2I3)') ystid, ivrs, zpobua(ista) &
                 , zzobua(ista), (zvcutua(ista,nba,ivrs), nba=1,2)             &
                               , (kviflua(ista,nba,ivrs), nba=1,2)
        ENDDO
      ENDIF
    ENDDO

  ENDIF   !  (k == ke)


!-------------------------------------------------------------------------------
!  Section 3: Organize the spreading of observation increments
!-------------------------------------------------------------------------------

! Get the horizontal correlation scales
! and the non-divergent correction factor to the 2-dim. wind correlations
! -----------------------------------------------------------------------

  DO istaua = 1 , nuatot
    ilv (istaua) = 0
    ista = isortua(istaua)
    IF (lsprsua(ista)) THEN
      IF (      (k <= MAX( kviflua(ista,1,1) , kviflua(ista,1,2)               &
                         , kviflua(ista,1,3) ))                                &
          .AND. (k >= MIN( kviflua(ista,2,1) , kviflua(ista,2,2)               &
                         , kviflua(ista,2,3) ))) THEN
        IF (zlopua(ista,k) <= tabcolp(ncolev)) THEN
          ilv (istaua) = ncolev
        ELSE
          ilv (istaua) = 2
          DO WHILE (zlopua(ista,k) <= tabcolp(ilv(istaua)))
            ilv (istaua) = ilv(istaua) + 1
          ENDDO
        ENDIF
      ENDIF
    ENDIF
  ENDDO

  DO istaua = 1 , nuatot
    ista = isortua(istaua)
    IF ((lsprsua(ista)) .AND. (ilv(istaua) >= 1)) THEN
      ilva = ilv(istaua)
      ilvb = ilva - 1
      zlopf = MIN( c1 , MAX( c0 , (tabcolp(ilvb) -zlopua(ista,k))              &
                                / (tabcolp(ilvb) -tabcolp(ilva)) ) )
      IF ((k <= kviflua(ista,1,1)) .AND. (k >= kviflua(ista,2,1))) THEN
        wtvmx               =  MAX( zwtua(ista,1,1) , zwtua(ista,2,1) )
        tnodivs             =  tnondiv * MIN( wtukara / wtukrsa , c1 )
        fvnodiv             = (c1-zlopf) *fnodivc(ilvb) + zlopf *fnodivc(ilva)
        fcnodiv (istaua,1)  =  MIN( c1 ,  (c1 + (tnodivs-c1) *(c1-wtvmx))      &
                                        * (cnondiv + fnondiv *fvnodiv) )
        rvinfl              = (c1-zlopf) *rhvsond(ilvb) + zlopf *rhvsond(ilva)
        rhtinfl (istaua,1)  =  MAX( epsy , (c1 + (rhtfac(1) -c1) *(c1- wtvmx)) &
                                          *(rhinfl(1) + rhvfac(1) *rvinfl))
      ENDIF
      IF ((k <= kviflua(ista,1,2)) .AND. (k >= kviflua(ista,2,2))) THEN
        wtvmx               =  MAX( zwtua(ista,1,2) , zwtua(ista,2,2) )
        fcnodiv (istaua,2)  =  c0
        rvinfl              = (c1-zlopf) *rhtsond(ilvb) + zlopf *rhtsond(ilva)
        rhtinfl (istaua,2)  =  MAX( epsy , (c1 + (rhtfac(3) -c1) *(c1- wtvmx)) &
                                          *(rhinfl(3) + rhvfac(3) *rvinfl))

      ENDIF
      IF ((k <= kviflua(ista,1,3)) .AND. (k >= kviflua(ista,2,3))) THEN
        wtvmx               =  MAX( zwtua(ista,1,3) , zwtua(ista,2,3) )
        fcnodiv (istaua,3)  =  c0
        rvinfl              = (c1-zlopf) *rhqsond(ilvb) + zlopf *rhqsond(ilva)
        rhtinfl (istaua,3)  =  MAX( epsy , (c1 + (rhtfac(4) -c1) *(c1- wtvmx)) &
                                          *(rhinfl(4) + rhvfac(4) *rvinfl))
      ENDIF
    ENDIF
  ENDDO

! start loop over upper-air stations
! ----------------------------------

loop_over_upper_air_stations:  DO istaua = 1 , nuatot

  ista = isortua(istaua)

  io  = ioua (ista)
  jo  = joua (ista)

  IF ((lsprsua(ista)) .AND. (ilv(istaua) >= 1)) THEN

    ystid  = ystidua (ista) (1:ilstidp)
    lwobs  = lwonl
    lnearonl = (ABS(io-ionl2) <= 5) .AND. (ABS(jo-jonl2) <= 5)

    io_tot = ioua_tot(ista)
    jo_tot = joua_tot(ista)

! ! NO loop over time index
! -------------------------
!   loop_over_time_index:  DO itim = 1, 2

! start loop over variables
! -------------------------

    loop_over_variables:  DO ivrs = 1, 3

      IF (ivrs == 1) ivar = 1
      IF (ivrs == 2) ivar = 3
      IF (ivrs == 3) ivar = 4

! check if vertical and temporal level may be influenced by current report
! ------------------------------------------------------------------------

      IF ((k <= kviflua(ista,1,ivrs)) .AND. (k >= kviflua(ista,2,ivrs))) THEN

! get the scaled horizontal distance to the station location, the geometrical
! factors of the 2-dim horizontal wind correlations, and the area of influence
! ----------------------------------------------------------------------------

        rtinfl  = rhtinfl(istaua,ivrs)
        zrcutof = cutofr(ivar) * rtinfl
        jrange  = INT( zrcutof * gppkmj )
        jstaspr = MAX( INT( jo - jrange ,iintegers) , jrs_tot )
        jendspr = MIN( INT( jo + jrange ,iintegers) , jre_tot )
        irange  = INT( zrcutof * MAX( gppkmi(jstaspr) , gppkmi(jendspr) ) )
!       irange  = INT( zrcutof * MAX( gppkmi(MIN( jstaspr,jend   ))            &
!                                   , gppkmi(MAX( jendspr,jstart )) ) )
        jstaspr = MAX( INT( jo - jrange ,iintegers) , jstart )
        jendspr = MIN( INT( jo + jrange ,iintegers) , jend )
        istaspr = MAX( INT( io - irange ,iintegers) , istart )
        iendspr = MIN( INT( io + irange ,iintegers) , iend )
        nactio  = 2
        IF (     (jstaspr > jend) .OR. (jendspr < jstart)                      &
            .OR. (istaspr > iend) .OR. (iendspr < istart) )  nactio  =  0

! IF ((lwonl) .AND. (k == ke) .AND. (ivrs == 2) .AND. (ntstep >= nstop-1)) THEN
!!  IF (ivrs == 2) PRINT *,'gppk ', jrs_tot,jre_tot, gppkmi
!   PRINT *,'ispr ',io, jo, '  ', ystid, '  ', iend, jend, ivrs, zrcutof       &
!                  ,irange, jrange, istaspr, jstaspr, iendspr, jendspr
! ENDIF

        CALL cutoff_wind_correl ( nactio , k , fcnodiv(istaua,ivrs) , ivar )
!       =======================

        omymx  = c0

        lsurf  = .FALSE.

        IF ((lprsys) .AND. (lwobs)) THEN
          kaob = (kviflua(ista,1,ivrs) + kviflua(ista,2,ivrs)) / 2
          ityp = isetyp0
          DO kk = SUM( niwtyp ) , 1 , -1
            IF (   ((iwtyp(kk) > 0) .AND. ( iwtyp(kk) == kobtyua(ista)))       &
               .OR.((iwtyp(kk) < 0) .AND. (-iwtyp(kk) == kcdtyua(ista))))      &
              ityp = isetyp(kk)
          ENDDO
          IF ((k == kaob) .AND. (lnearonl) .AND. (ivrs == 1))                  &
            WRITE( nupr,'(''air-infl(km)'',F5.0, 1X,A , 2I4, F7.0              &
                        &,'', zu/vwi,omyu'',2F8.3,F8.4)' )                     &
                   zrcutof, ystid, io_tot, jo_tot, zpobua(ista)                &
                 , zwi(ionl2,jonl2,1,ityp), zwi(ionl2,jonl2,2,ityp)            &
                 , omy(ionl2,jonl2,1,ityp)
        ENDIF

! spreading of observation increments
! -----------------------------------

        IF ((k >= kviflua(ista,2,ivrs)) .AND. (nactio > 0))                    &

          CALL surf_upair_spread ( lsurf , ivar , k , omymx )
!         ======================

! printout for control
! --------------------

        IF ((lprsys) .AND. (lwobs)) THEN
          kaob = (kviflua(ista,1,ivrs) + kviflua(ista,2,ivrs)) / 2
          IF ((k == kaob) .AND. (lnearonl) .AND. (ivrs == 1))                  &
            WRITE( nupr,'(''air-infl(km)'',F5.0, 1X,A , 2I4, F7.0              &
                        &,'', zu/vwi,omyu'',2F8.3,F8.4)' )                     &
                   zrcutof, ystid, io_tot, jo_tot, zpobua(ista)                &
                 , zwi(ionl2,jonl2,1,ityp), zwi(ionl2,jonl2,2,ityp)            &
                 , omy(ionl2,jonl2,1,ityp)
        ENDIF
        IF ((lprthr) .AND. (lwobs)) THEN
          IF (      (     (k <= kviflua(ista,2,ivrs))                          &
                     .OR. (k == kviflua(ista,1,ivrs)))                         &
              .AND. ((lfirst) .OR. (k < ke) .OR. (ivar == 1))) THEN
            IF (ivar == 1) THEN
              xincrp = xoiua(ista,1,ivar) + xoiua(ista,2,ivar)
              yincrp = xoiua(ista,1,2   ) + xoiua(ista,2,2   )
              WRITE( nupr,'(''thrair uv k'',3I3,2I4,1X, A , 2(1X,F6.0)         &
                          &,'',max omy:'',F6.3,'',uv incr'',2F6.1)' )      k   &
                   , kviflua(ista,1,ivrs), kviflua(ista,2,ivrs), io_tot,jo_tot &
                   , ystid, zsprob(ista), zsprua(ista,k), omymx, xincrp, yincrp
            ELSEIF (ivar == 3) THEN
              WRITE( nupr,'(''thrair T  k'',3I3,2I4,1X, A , 2(1X,F6.0)         &
                          &,'',max omy:'',F6.3,'',T  incr'',2F6.1)' )      k   &
                   , kviflua(ista,1,ivrs), kviflua(ista,2,ivrs), io_tot,jo_tot &
                   , ystid, zsprob(ista), zsprua(ista,k), omymx                &
                   , xoiua(ista,1,ivar), xoiua(ista,2,ivar)
            ENDIF
          ENDIF
        ENDIF

      ENDIF

    ENDDO loop_over_variables

!   ENDDO loop_over_time_index

  ENDIF ! lsprsua(ista)

ENDDO loop_over_upper_air_stations


! deallocate memory
! -----------------

  IF (k == 1) DEALLOCATE ( lsprsua , STAT=istat )

  IF (ldump_ascii) THEN
    ! flush YUPRINT file
    IF ((lwobs) .AND. ((lprthr) .OR. (lprsys))                                 &
                .AND. ((k == ke) .OR. (k == 1) .OR. (k >= ke-6))) THEN
      CLOSE (nupr)
      OPEN  (nupr,FILE=yuprint,FORM='FORMATTED',STATUS='OLD'                   &
                              ,POSITION='APPEND',IOSTAT=nstat)
      IF (nstat /= 0) PRINT '("OPENING OF FILE yuprint FAILED, upair_org_spread")'
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure upair_org_spread
!-------------------------------------------------------------------------------

END SUBROUTINE upair_org_spread


!===============================================================================

ELEMENTAL REAL (KIND=wp) FUNCTION rmod  ( ztti, ziv )
  !----------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  ztti, ziv
  !  uses parameter:                           epsy
  !---------------------------------------------------------------------------
  ! MOD function for positive REALS
  !---------------------------------------------------------------------------
  !
  rmod  =  ztti  -  ziv * REAL(INT( ztti/ziv + epsy ),wp) + epsy
  !
END FUNCTION rmod

!-------------------------------------------------------------------------------

END MODULE src_sing_spread
