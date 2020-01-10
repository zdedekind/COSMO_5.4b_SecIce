!+ Source module for spreading multi-level obs. increments for the nudging
!-------------------------------------------------------------------------------

MODULE src_mult_spread

!-------------------------------------------------------------------------------
!
! Description:
!   The module "src_mult_spread" performs the spreading of observational
!   information from all the multi-level data from the total model domain
!   to the model grid points on the local sub-domain.
!
! Method:
!   This module contains the following procedures:
!    - mult_find_levels   : find obs increment levels for vertical interpolation
!    - mult_org_spread    : organize spreading of multi-level obs increments
!    - mult_spread_mass   : spreading scalar multi-level obs increments
!    - mult_spread_modlev : spreading multi-level obs incr. along model levels
!    - mult_spread_wind   : spreading multi-level obs increments of horiz. wind
!   Driving routine is the module procedure "mult_org_spread".
!
!   This module also contains an elemental function, formerly statement funct.:
!    - rmod               : MOD function for positive reals
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
! 1.15       1998/11/02 Christoph Schraff
!  Global declaration of allocatable arrays moved to module 'data_nudge_spread'.
! 1.19       1998/12/11 Christoph Schraff
!  ANSI violations removed.
! 1.20       1999/01/07 Guenther Doms
!  Renaming of some global variables
! 1.31       1999/07/01 Christoph Schraff
!  Quantities related to MPI communicator 'icomm_world' replaced.
! 2.4        2001/01/29 Christoph Schraff
!  Reorganisation of loops to enhance efficiency on vector processors and to
!  render spreading loops more uniform.
! 2.5        2001/06/01 Christoph Schraff
!  Savety test at array allocations.
! 2.13       2002/01/18 Christoph Schraff
!  For 'hot' IBM compiler option: Data statements replaced by direct assignment.
! 2.14       2002/02/15 Ulrich Schaettler
!  Introduced compiler directives for vectorization on the NEC
!  (after Mauro Ballabio, CSCS)
! 3.3        2003/04/22 Maria Tomassini + Christoph Schraff
!  Extension for assimilation of GPS-derived IWV.
! 3.6        2003/12/11 Christoph Schraff
!  Horizontal correlation scale for GPS data given by namelist variable 'rhfgps'
! 3.12       2004/09/15 Christoph Schraff
!  Extension to (prepare to) include assimilation of satellite retrievals.
!  Minor bug correction related to the spreading in the lowest model layer.
! 3.18       2006/03/03 Christoph Schraff
!  New option for alternative weighting for multiple observations.
!  Option for separate weighting for different observation types. To allow for
!  reproducibility, use of same loops if 'ltiml', irrespective of 'lsprx'.
! V3_23        2007/03/30 Ulrich Schaettler
!  Introduced ldump_ascii to flush the ASCII files
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_22        2012/01/31 Christoph Schraff
!  - Optimisations for decreased cost on NEC-SX9 (e.g. by putting indirect
!     addressing in lower-dimensional loop in routine mult_find_levels;
!     or by vectorisation in sections 1, 3 of subroutine mult_org_spread).
!  - Some variables moved from 'data_nudge_all' to 'data_obs_lib_cosmo'.
!  - Subroutine arguments lists extended by vertical model level 'k'.
! V4_28        2013/07/12 Christoph Schraff
!  Statement function 'rmod' replaced by elemental function.
! V5_1         2014-11-28 Christoph Schraff, Oliver Fuhrer
!  Processing of Mode-S aircraft obs introduced (CS)
!  Replaced ireals by wp (working precision) (OF)
! V5_2         2015-05-21 Ulrich Schaettler
!  Initialized lnoniso with FALSE in SR mult_spread_mass and mult_spread_wind
!   because there could be an uninitialized use otherwise
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

    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction
    ke1,          & ! ke+1

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

    hhl               ! geometical height of half model levels        ( m )

! end of data_fields

!-------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------

    nstop,        & ! last time step of the forecast
    ntstep,       & ! actual time step

! 7. additional control variables
! -------------------------------

    ldump_ascii     ! for flushing (close and re-open) the ASCII files

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

    niwtyp       ,& ! 1,0,0,..: number of obs or code types which belong to
                    !           the sets of observing systems
    iwtyp        ,& ! 0,0,0,..: obs or code types belonging to a set of obs
                    !           system, specified successively for each set
    kwtyp        ,& ! 1   : function for weights W for multiple observations
    nudgend      ,& ! 0   : end of nudging period in timesteps
    tconbox      ,& ! 6*dt: timestep [s] for computing analysis increments
                    !       (i.e. time box of constant analysis increments)
    gnudg        ,& ! 6,12,6,6*10^-4: nudging coefficients for TEMP / PILOT data
    gnudgar      ,& ! 6, 0,6,0*10^-4: nudging coeffic. for AIRCRAFT data [1/s]
    gnudgms      ,& ! 6, 0,6,0*10^-4: nudging coeffic. for Mode-S aircraft [1/s]
    gnudgtv      ,& ! 0, 0,6,6*10^-4: nudging coeffic. for sat retrivals [1/s]
    gnudggp      ,& !       0 *10^-4: nudging coeffic. for GPS-derived IWV [1/s]
    msprpar      ,& ! 1   : switch specifying the surfaces along which observat.
                    !       increments of upper-air data are (primarily) spread 
    vcorls       ,& ! 2*.333,: square of the vertical correlation scale,
                    ! 2*.04    i.e. of the Gaussian vertical influence 'radius'
                    !          in log( pressure ) if (msprpar <= 1), or
                    !          in potential temperature if (msprpar == 2)
    wablua       ,& ! 4*  1. : factor to weights within the ABL
                    !          (atmospheric boundary layer, i.e. mixed layer)
    rhinfl       ,& ! 0.,70.,: constant part of the 'correlation scale of the
                    ! 0.,0.    autoregressive horiz. correlation function'
                    !          (='COSAC') [km]
    rhvfac       ,& ! 1., 0.,: multiplication factor to the vertically varying
                    ! 2* .83   part of the 'COSAC' (as def. in 'data_nudge_all')
    rhtfac       ,& ! 1.3,   : scaling factor of the total 'COSAC' for the
                    ! 1.43,    beginning and end of the nudging period for 1 obs
                    ! 1.3,     relative to the 'COSAC' at the obs. time as given
                    ! 1.3      by 'rhinfl', 'rhvfac')
    rhfrtv       ,& ! 1., 1.,: scaling (reduction) factor of the total 'COSAC'
                    ! 1., 0.5  for temperature / humidity satellite retrievals
    rhfgps       ,& ! 0.45   : scaling (reduction) factor of the total 'COSAC'
                    !          for humidity derived from GPS IWV
    cutofr       ,& ! 4* 3.5 : cut-off in 'COSAC' units of the horizontal
                    !          correlation function
    vcsni        ,& ! 4*2500.: square of Gaussian vertical influence 'radius'
                    !          in potential temperature (if msprpar <= 1) or
                    !          log( pressure ) (if msprpar == 2) on surfaces
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
    ionl         ,& ! 167    : / grid point coordinates
    jonl         ,& ! 103    : \ for standard output on nudging

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
    rmdi       ,& ! =-1.E31_wp : commonly used missing data indicator
    rmdich     ,& ! =-1.E30_wp : commonly used check value for miss data
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

    nairep     ,& ! AIREP reports (all aircraft reports)
    ntemp      ,& ! TEMP  reports
    npilot     ,& ! PILOT reports
    nsattv     ,& ! SATEM reports
    ngps       ,& ! GPS   reports
    nmodes        !   mode-s report

! end of data_obs_lib_cosmo

!-------------------------------------------------------------------------------

USE data_nudge_gather , ONLY :   &

! 1. Parameters and general variables
! -----------------------------------

    ilstidp      ,& ! char. length used for printing the station ID
                    ! Note: (ilstid >= ilstidg >= ilstidp) !!!
    ystid        ,& ! obs. station identity to be printed
    zalllow      ,& ! smallest allowed height within model domain
    zallhig      ,& ! largest  allowed height within model domain
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

! local information on multi-level reports
    oiml         ,& ! observation increment record OIR for multi-level reports
    zwtml        ,& ! temporal weights
    zspobml      ,& ! spreading parameter at the base / top of the profile
    fcorlml      ,& ! reduction to vertical correl. at the profile's base /top
    zvcutml      ,& ! vertical cut-off at the base / top of the profile
    zrtdgml      ,& ! convertor from height to pressure units for vertical
                    ! correlations at the base / top of the profile
    zpobml       ,& ! pressure (p) at the base / top of the profile
    zemkpml      ,& ! EXP( -R/cp *p ) at the base / top of the profile
    zriflml      ,& ! upper estimate to horizontal radius of influence
    zstpml       ,& ! spreading parameter at the tropopause & obs. location
    zsprml       ,& ! spreading parameter at model levels & obs. location
    zlopml       ,& ! log( pressure ) at model levels & obs. location
    zrtgkml      ,& ! convertor from height to pressure units for vertical
                    ! correlations at model levels & obs. location
    znisml       ,& ! parameter used for non-isotropic horizontal correlations
                    ! at model levels & obs. location
    znismq       ,& ! parameter used for non-isotropic horizontal correlations
                    ! at vertically equidistant points & obs. location
    ioml         ,& ! longitudinal index of station location on local domain
    joml         ,& ! latitudinal  index of station location on local domain
    ioml_tot     ,& ! longitudinal index of station location on global domain
    joml_tot     ,& ! latitudinal  index of station location on global domain
    kviflml      ,& ! vertical range of possibly influenced model levels
    kkoiml       ,& ! level index in the OIR (observation increment record)
    ksprml       ,& ! lowest level with spreading along model levels
    kobtyml      ,& ! CMA observation type
    kcdtyml      ,& ! CMA observation code type
    mszlev       ,& ! number of vertical levels with obs. incr. in the OIR
    mbotlv       ,& ! number of model levels above the surface without
                    ! observation increment
    mtoplv       ,& ! number of model levels at the top of the model without
                    ! observation increment
    ltiml        ,& ! .TRUE if temporal linear interpol. for multi-level rep.
    ystidml         ! station identity of multi-level station

USE data_nudge_gather , ONLY :   &

! 3. Observation increment record for multi-level reports 'oiml'
! --------------------------------------------------------------

    maxnoi       ,& ! length of observation increment record
    noiu         ,& ! zonal wind observation increment
    noiv         ,& ! meridional wind observation increment
    noit         ,& ! temperature observation increment
    noiqd        ,& ! specific humidity observation increment
    noirh        ,& ! relative humidity observation increment
    noiz         ,& ! height (of obs. increment level)
    noith        ,& ! potential temperature
    noilp        ,& ! log( pressure )
    noiuqc       ,& ! quality weight to wind obs. increment
    noitqc       ,& ! quality weight to temperature obs. increment
    noiqqc       ,& ! quality weight to humidity obs. increment
    noivcr       ,& ! normalization factor to vertical interpolation distances
    noiulr       ,& ! vertical interpolat. distance to obs. incr. level for wind
    noitlr       ,& ! vertical int'pol distance to obs. incr. level for temper.
    noiqlr       ,& ! vertical int'pol distance to obs. incr. level for humidity

! 5. Variables defining the size of the arrays containing the local information
! -----------------------------------------------------------------------------

    nmltot       ,& ! total number of active multi-level stations
    lnisua          ! non-isotrophic correlat. for upper-air single-lev. data

! end of data_nudge_gather

!-------------------------------------------------------------------------------

USE data_nudge_spread , ONLY :   &

! 3. Output of spreading procedures
! ---------------------------------

    omy          ,& ! sum of  spatial * temporal * quality 'spreading weights'
    om2          ,& ! sum of  squares of 'spreading weights'
    zwi          ,& ! sum of weighted (observation) increments (weights are
                    !        squares of 'spreading weights')

! 6. Further horizontal input fields for the spreading of obs. increments
! -----------------------------------------------------------------------

    zeklop       ,& ! exp( R/cp * log(p) )
    zthvg        ,& ! vertical gradient of potential temperature
    zspr         ,& ! spreading parameter , param. def. non-isotropic weights
    zdds         ,& ! scaled horizontal distance betw. obs. and target grid pt
    zcoruu       ,& ! zonal  wind - zonal  wind correlation  \  (without
    zcoruv       ,& ! zonal  wind - merid. wind correlation   \  EXP( -zdds )
    zcorvu       ,& ! merid. wind - zonal  wind correlation   /  -term )
    zcorvv       ,& ! merid. wind - merid. wind correlation  /
    zablfc       ,& ! reduced weighting inside/above ABL for upper-air/sfc obs

! 7. Further fields and variables used for or during the spreading
! ----------------------------------------------------------------

    klva         ,& ! indices of upper  __\  obs. increment levels used for
    klvb         ,& ! indices of lower    /  vertical interpolation to grid pt
    zoic1        ,& ! spreaded obs. incr. from first multi-level report
    zoiv1        ,& ! 2nd component of spreaded obs. incr. from first report
    omys1        ,& ! weight from first report
    zsdnis       ,& ! equidistant pts. used for non-isotropic horiz. correlat.
    zsdnid       ,& ! distance between two adjacent 'equidistant' pts.
    gppkmi       ,& ! convertor for zonal dist. from 'km' to rotated grid pts.
    gppkmj       ,& ! convertor for merid. dist. from 'km' to rotated grid pts
    rtinfl       ,& ! horizontal correlation scale
    lcutof       ,& ! .TRUE if grid pt. is within area of influence of report
    icutof       ,& ! x-coordinate of grid pt. within area of influence
    jcutof       ,& ! x-coordinate of grid pt. within area of influence
    ncutof       ,& ! number of grid points within area of influence

! 8. Indices and lists of indices
! -------------------------------

    isortml      ,& ! (sorted) list of stations with multi-level data
    io    ,jo    ,& ! local  indices of location of observation
    io_tot,jo_tot,& ! global indices of location of observation
    istaspr      ,& ! /  lower left corner of domain
    jstaspr      ,& ! \  containing area of influence
    iendspr      ,& ! /  upper right corner of domain
    jendspr      ,& ! \  containing area of influence
    jrs_tot      ,& ! /  index range for convertor
    jre_tot      ,& ! \  for zonal distances 'gppkmi'
    ista         ,& ! index of observing station
    itim            ! (time) index over obs. at one observing station

USE data_nudge_spread , ONLY :   &

! 9. Observation types in weighted increment arrays 'zwi'
! -------------------------------------------------------

    mxotyp       ,& ! number of observation types
    noiras       ,& ! radiosonde
    noisfc       ,& ! surface-level (SYNOP, SHIP, BUOY)
    noiair       ,& ! aircraft
    noisat       ,& ! satellite (SATOB, ATOVS, MSG)
    noigps          ! GPS

! end of data_nudge_spread

!-------------------------------------------------------------------------------

USE src_correl_cutoff,        ONLY :  &
      cutoff_wind_correl        ! areas of influence of observations


!-------------------------------------------------------------------------------
! Local declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE

! 1. Parameters
! -------------

  REAL (KIND = wp)         , PARAMETER    , PRIVATE :: &
    thresh1 = 16.39_wp ,& ! threshold for the vertical potential temperature
                              ! gradient above which increments are spread pure-
                              ! ly along isentropic surfaces (if msprpar == 2)
    thresh2 =  0.00_wp ,& ! threshold for the vertical potential temperature
                              ! gradient above which increments are spread along
                              ! horizontal surfaces, if (msprpar == 2)
    qs2     = qst * qst       ! qst **2

  LOGICAL                  , PARAMETER    , PRIVATE :: &
    lautorg = .TRUE.          ! .FALSE. ==> Cressman instead of autoregressive
                              !             horizontal correlation function


! 2. Variables
! ------------

  INTEGER (KIND=iintegers)                , PRIVATE :: &
    kobtyp         ,& ! observation type
    kcdtyp            ! CMA observation code type


!===============================================================================

!-------------------------------------------------------------------------------
! Public and Private Subroutines
!-------------------------------------------------------------------------------

CONTAINS

! INCLUDE "mult_find_levels.incf"
! INCLUDE "mult_org_spread.incf"
! INCLUDE "mult_spread_mass.incf"
! INCLUDE "mult_spread_modlev.incf"
! INCLUDE "mult_spread_wind.incf"

!-------------------------------------------------------------------------------
!+ Module procedure to find obs increment levels for vertical interpolation
!-------------------------------------------------------------------------------

SUBROUTINE mult_find_levels ( ivrc , k , modex )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure finds the level indices of the observation increments
!   used for vertical interpolation and subsequent lateral spreading
!   (extrapolation) to the target grid points along surfaces
!   (horizontal or isentropic) which differ from the model surfaces.
!
! Method:
!   First limit as far as possible the vertical range of indices used for the
!   search. Then find the required indices by comparing the values of the 
!   spreading parameter at the observation points successively in an outer
!   vertical loop with the values of the spreading parameter at the target
!   grid points from an inner loop over the horizontal area of influence.
!   For each target grid point, the required indices belong those (max. 2)
!   observation points, which contain increments of type 'ivrc', and which
!   have values of the spreading parameter which are closest (above resp. below)
!   to the value at the target grid point.
!   Spreading parameter is height (modex=1) or potential temperature (modex=2).
!
! Written by        :  Christoph Schraff, DWD  (original version: 03.09.97)
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
    ivrc             ,& ! index of spreaded quantity and counter of report:
                        ! 1=(u,v); 2=T; 3=RH for the 1st or the only report;
                        ! 4=(u,v); 5=T; 6=RH for the 2nd relevant report
    k                   ! index of current vertical model level


  INTEGER (KIND=iintegers), INTENT (INOUT)      ::       &
    modex               ! defines spreading parameter (see above 'method';
                        ! (modex=1): height ; (modex=2): potential temperature


! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
!   modespr          ,& ! mode of spreading, see above: 'Method'
    irep             ,& ! index of '?cutof': counter of report per station
    ivrs             ,& ! index of spreaded quantity: 1=(u,v); 2=T; 3=RH
    ivar             ,& ! index of spreaded quantity: 1=(u,v); 3=T; 4=RH
    isvv             ,& ! index of spreaded quantity depend. on spreading param.
!   nmloi            ,& ! index of report in the observation increment field
    noix             ,& ! index of spreaded quantity in obs. increment field
    nois             ,& ! index of (selected) spreading param. in obs incr field
    msbotx , mstopx  ,& ! index range used for the vertical search
    mbtopx           ,& ! auxiliary used to find that index range
    kkbot  , kktop   ,& ! index range used for the present report
    i      , j       ,& ! loop indices in horizontal direction
    klv              ,& ! loop indices in vertical direction
    ixc              ,& ! loop index over influenced grid points
!   io     , jo      ,& ! local indices of location of observation
!   istaspr, jstaspr ,& ! lower left corner of domain containing area of infl.
!   iendspr, jendspr ,& ! upper right corner of domain containing area of infl
    nprae            ,& ! number of printed grid pts.
    jpra   , jpre    ,& ! lowest , highest meridional index of printed grid pts.
    joc                 ! loop index for printout

  INTEGER (KIND=iintegers) ::  &
    klva_1d  (ie*je) ,& ! indices of upper __\  obs. incr. levels used for
    klvb_1d  (ie*je)    ! indices of lower   /  vert. interpol. to grid pt

  REAL    (KIND=wp   )     ::  &
    zspr_1d  (ie*je) ,& ! spreading parameter (auxilliary 1-dim array)
    zzmax            ,& ! height of highest __\ target grid pt. within area of
    zzmin            ,& !         / lowest    / influence on current model level
    zvgmin              ! min. vertical poten. temperature gradient in that area

  CHARACTER (LEN=15) yformat

! Local (automatic) arrays: None
! -------------------------

! Initialized variables:
! ----------------------

! LOGICAL                  , SAVE ::  &
!   lfrsloc = .TRUE.    ! .TRUE if local output is to be printed the first time

!------------ End of header ----------------------------------------------------

 
!-------------------------------------------------------------------------------
! Begin Subroutine mult_find_levels
!-------------------------------------------------------------------------------
  
!-------------------------------------------------------------------------------
!  Section 1: Preliminaries
!-------------------------------------------------------------------------------
 
! preset some values, and determine the height of the lowest and highest target
! grid pt. within the horizontal area of influence at the current model level
! -----------------------------------------------------------------------------

  ivrs   = MOD( ivrc-1 , 3 ) + 1
  irep   =     (ivrc-1)/ 3   + 1
! isvv   = 3*(modex-1) + ivrs
  isvv   =    modex
! nmloi  = knoiml(ista,itim)

  zzmax  = -400.0_wp
  zzmin  = 8849.0_wp
  zvgmin = thresh1 + c1
  DO i   = istaspr , iendspr
    DO j = jstaspr , jendspr
      klva (i,j,isvv) = -1
      klvb (i,j,isvv) = -1
    ENDDO
  ENDDO
  IF (modex == 2) THEN
!CDIR NODEP
    DO ixc = 1 , ncutof(irep)
      i = icutof(ixc,irep)
      j = jcutof(ixc,irep)
      zvgmin = MIN( zvgmin , zthvg(i,j) )
    ENDDO
  ELSE
!CDIR NODEP
    DO ixc = 1 , ncutof(irep)
      i = icutof(ixc,irep)
      j = jcutof(ixc,irep)
      zzmax = MAX( zzmax , zspr(i,j,modex) )
      zzmin = MIN( zzmin , zspr(i,j,modex) )
    ENDDO
  ENDIF

! determine the vertical range of indices used for the subsequent search
! ----------------------------------------------------------------------

  IF (ivrs == 1) noix = noiu
  IF (ivrs == 2) noix = noit
  IF (ivrs == 3) noix = noirh
  IF (modex == 2) THEN
    nois   = noith
    msbotx = kkoiml(ista,itim)
    mstopx = kkoiml(ista,itim) + mszlev(ista,itim) - 1
  ELSE
    nois   = noiz
    msbotx = -1
    mstopx = -2
    mbtopx = -2
    kkbot  = kkoiml(ista,itim)
    kktop  = kkoiml(ista,itim) + mszlev(ista,itim) - 1
    DO klv = kkbot , kktop
      IF (oiml(noix,klv) > rmdich) THEN
        IF ((oiml(nois,klv) <  zzmin)  .OR. (msbotx == -1)) msbotx = klv
        IF ((oiml(nois,klv) >  zzmax) .AND. (mstopx == -2)) mstopx = klv
        IF  (oiml(nois,klv) <= zzmax)                       mbtopx = klv
      ENDIF
    ENDDO
    mstopx = MAX( mstopx , mbtopx )
  ENDIF

!-------------------------------------------------------------------------------
!  Section 2: Find the indices of the required observation increments
!-------------------------------------------------------------------------------

!CDIR NODEP,VOVERTAKE,VOB
!CDIR ON_ADB(zspr_1d)
  DO ixc = 1 , ncutof(irep)
    i = icutof(ixc,irep)
    j = jcutof(ixc,irep)
    zspr_1d (ixc) = zspr(i,j,modex) 
    klva_1d (ixc) = -1
    klvb_1d (ixc) = -1
  END DO
  ! in the following 3-dim loop, indirect addressing 
  DO klv = mstopx , msbotx , -1
    IF (oiml(noix,klv) > rmdich) THEN
!CDIR NODEP,VOVERTAKE,VOB
!CDIR ON_ADB(zspr_1d)
!CDIR ON_ADB(klva_1d)
!CDIR ON_ADB(klvb_1d)
      DO ixc = 1 , ncutof(irep)
        IF (oiml(nois,klv) > zspr_1d(ixc)) THEN
          klva_1d (ixc) = klv
        ELSEIF (klvb_1d(ixc) == -1) THEN
          klvb_1d (ixc) = klv
        ENDIF
!     ! the above replaces indirect addressing below within this costly 3-d loop
!       i = icutof(ixc,irep)
!       j = jcutof(ixc,irep)
!       IF (oiml(nois,klv) > zspr(i,j,modex)) THEN
!         klva (i,j,isvv) = klv
!       ELSEIF (klvb(i,j,isvv) == -1) THEN
!         klvb (i,j,isvv) = klv
!       ENDIF
      ENDDO
    ENDIF
  ENDDO
!CDIR ON_ADB(zspr_1d)
!CDIR ON_ADB(klva_1d)
!CDIR ON_ADB(klvb_1d)
  DO ixc = 1 , ncutof(irep)
    i = icutof(ixc,irep)
    j = jcutof(ixc,irep)
    klva (i,j,isvv) = klva_1d (ixc)
    klvb (i,j,isvv) = klvb_1d (ixc)
  ENDDO

!! The following construction with vector registers is not faster
!!INTEGER (KIND=iintegers), PARAMETER :: n_vec_reg_size = 256  ! size of vec-reg
!!REAL    (KIND=wp   )     :: zspr_reg(n_vec_reg_size) ! vector registers
!!INTEGER (KIND=iintegers) :: klva_reg(n_vec_reg_size), klvb_reg(n_vec_reg_size)
!!!CDIR VREG(zspr_reg); !!!CDIR VREG(klva_reg); !!!CDIR VREG(klvb_reg)
!! IF (ncutof(irep) > c05* n_vec_reg_size) THEN
!!   DO ii = 1 , ncutof(irep) , n_vec_reg_size
!!!CDIR NODEP,SHORTLOOP
!!     DO ixc = ii , MIN(ii+n_vec_reg_size-1,ncutof(irep))     ! 1 , ncutof(irep)
!!       i = icutof(ixc,irep)
!!       j = jcutof(ixc,irep)
!!       zspr_reg (ixc+1-ii) = zspr(i,j,modex)
!!       klva_reg (ixc+1-ii) = -1
!!       klvb_reg (ixc+1-ii) = -1
!!     ENDDO
!!     DO klv = mstopx , msbotx , -1
!!       IF (oiml(noix,klv) > rmdich) THEN
!!!CDIR NODEP,SHORTLOOP
!!         DO ixc = ii , MIN(ii+n_vec_reg_size-1,ncutof(irep)) ! 1 , ncutof(irep)
!!           IF (oiml(nois,klv) > zspr_reg(ixc+1-ii)) THEN
!!             klva_reg (ixc+1-ii) = klv
!!           ELSEIF (klvb_reg(ixc+1-ii) == -1) THEN
!!             klvb_reg (ixc+1-ii) = klv
!!           ENDIF
!!         ENDDO
!!       ENDIF
!!     ENDDO
!!!CDIR NODEP,SHORTLOOP
!!     DO ixc = ii , MIN(ii+n_vec_reg_size-1,ncutof(irep))     ! 1 , ncutof(irep)
!!       i = icutof(ixc,irep)
!!       j = jcutof(ixc,irep)
!!       IF (klva_reg(ixc+1-ii) /= -1) klva(i,j,isvv) = klva_reg(ixc+1-ii)
!!       IF (klvb_reg(ixc+1-ii) /= -1) klvb(i,j,isvv) = klvb_reg(ixc+1-ii)
!!     ENDDO
!!   ENDDO
!! ENDIF

! printout for control
! --------------------

  IF ((ntstep <= 1) .AND. (lwonl)) THEN
!   io = imladm(ista,3)   !  io,jo known globally
!   jo = imladm(ista,4)
!   IF (((io == ionl) .AND. (jo == jonl) .AND. (k >= ke-6)) .OR. (lfrsloc)) THEN
    IF  ((io == ionl) .AND. (jo == jonl) .AND. (k >= ke-6)) THEN
      kkbot  = kkoiml(ista,itim)
      kktop  = kkoiml(ista,itim) + mszlev(ista,itim) - 1
      jpra = MAX( INT( jo  -5 ,iintegers) , jstaspr )
      jpre = MIN( INT( jpra+9 ,iintegers) , jendspr )
      jpra = MAX( INT( jpre-9 ,iintegers) , jstaspr )
      nprae = jo   - jpra + 1
      WRITE( nupr,'(''findlev: k='',I2,'', modespr='',I1,'', msbot/top='',2I3  &
                  &,'', zzmax/min'',2F7.0,'', jo'',I2,'', klvb/a'')' )         &
             k, modex, msbotx, mstopx, zzmax, zzmin, nprae
      nprae = jpre - jpra + 1
      ivar = ivrs + MIN( ivrs-1 , 1 )
      i    = MAX( INT( MIN( io , iendspr ) ,iintegers) , istaspr )
      WRITE( yformat, '(''('',I2,''I7,3I4     )'')' ) nprae
      WRITE( nupr, yformat )  (klvb (i,joc,isvv ), joc=jpra,jpre), ivar, i, jpra
      WRITE( nupr, yformat )  (klva (i,joc,isvv ), joc=jpra,jpre), ivar, io, jo
      WRITE( yformat, '(''('',I2,''F7.1,I4    )'')' ) nprae
      WRITE( nupr, yformat )  (zspr (i,joc,modex), joc=jpra,jpre), ivar
      nprae = MIN( kktop - kkbot + 1 , 10 )
      WRITE( yformat, '(''('',I2,''F7.1,I4    )'')' ) nprae
      WRITE( nupr, yformat )  (oiml (nois,klv), klv= kkbot,kkbot+nprae-1), ivar
      WRITE( nupr, yformat )  (oiml (noix,klv), klv= kkbot,kkbot+nprae-1), ivar
      IF (kktop >= kkbot+10) THEN
        nprae = MIN( kktop - kkbot - 9 , 10 )
        WRITE( yformat, '(''('',I2,''F7.1,I4    )'')' ) nprae
        WRITE( nupr,yformat ) (oiml (nois,klv), klv=kkbot+10,kkbot+nprae+9),ivar
        WRITE( nupr,yformat ) (oiml (noix,klv), klv=kkbot+10,kkbot+nprae+9),ivar
      ENDIF
!     lfrsloc = .FALSE.
    ENDIF
  ENDIF

! prepare horizontal spreading, if it is required
! although spreading along isentropes is selected
! -----------------------------------------------

  IF ((zvgmin <= thresh1) .OR. (k == ktth)) modex = 1

!-------------------------------------------------------------------------------
! End of module procedure mult_find_levels
!-------------------------------------------------------------------------------

END SUBROUTINE mult_find_levels


!-------------------------------------------------------------------------------
!+ Module procedure to organize spreading of multi-level obs increments
!-------------------------------------------------------------------------------

SUBROUTINE mult_org_spread ( k )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure organizes the spreading of observation increments
!   from multi-level reports after the processor has received all the required
!   'local' information from all observations.
!
! Method:
!   Setting up the loop over the stations with multi-level reports.
!   Exclusion of all reports / data, for which no grid pts. of the current node
!   and vertical and temporal model level exist in the 4-dim area of influence.
!   Area of influence, spatial correlations, and spreading by call of other
!   procedures. 
!
! Written by        :  Christoph Schraff, DWD  (original version: 04.09.97)
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
    modespr          ,& ! selected mode of spreading multi-level data
    modex            ,& ! mode of spreading used to find indices of obs. incr.
    istaml           ,& ! index of observing station
!   ista             ,& ! index of observing station
    ivrs             ,& ! index of spreaded quantity: 1=(u,v); 2=T; 3=RH
    ivar             ,& ! index of spreaded quantity: 1=(u,v); 3=T; 4=RH
    ivrc             ,& ! index of spreaded quantity combined with report count:
                        ! 1=(u,v); 2=T; 3=RH for the 1st or the only report;
                        ! 4=(u,v); 5=T; 6=RH for the 2nd relevant report
    ityp             ,& ! index of weighted increment array
    irange , jrange  ,& ! zonal , meridional 'radii' of the domain which
                        ! contains the area of influence
    i      , j       ,& ! loop indices in horizontal direction
    izext  , jzext   ,& ! coord. of lowest or highest grid pt. within infl. area
    ivex             ,& ! 'ivar', or (compressed) coord. (izext,jzext), or
                        ! 1 (0) if node has (no) grid pts. in area of influence
    ilva   , ilvb    ,& ! indices of correlation scale table levels that are
                        ! adjacent to the specified model level
    nactio           ,& ! action executed in called proc. 'cutoff_wind_correl'
    nba    , nbah    ,& ! loop indices for lower /upper cut-off and 1./2. report
    nbotop           ,& ! indicator for lower resp. upper cut-off
    kk               ,& ! vertical loop index
    kviflbot         ,& ! /  range (indices) of vertical levels
    kvifltop         ,& ! \  influenced by current report
    istat  , nstat      ! status of memory (de-)allocation

  REAL    (KIND=wp   )     ::  &
    acthr            ,& ! actual forecast hour
    zrwt             ,& ! temporal nudging weight used for correlation scales
    zrcutof          ,& ! max. horizontal cut-off radius [km]
    z2cut            ,& ! 2 * zvcutml (zvcutml: height at vertical cut-off)
    zlopf               ! weight factor for vertical interpolations from the
                        ! correlation scale table

  LOGICAL                  ::  &
    lsprx            ,& ! non-zero temporal and vertical weights in current rep.
    lraso            ,& ! spreading of obs. incr. from TEMP or PILOT (not AIREP)
    larea            ,& ! first guess for area of influence partly in sub-domain
    lwobs               ! (TRUE if obs. station (io,jo) lies in local domain):
                        !  replaced by lwobs = lwonl
!   lprsys , lprthr     ! for printing for control

! Local (automatic) arrays:
! -------------------------

  LOGICAL                  ::  &
    lsprsv  (nmltot,3)   ! temporal and horiz. domain of node influenced by obs.

  INTEGER (KIND=iintegers) ::  &
    ijspr   (4,nmltot),& ! domain containing area of influence
    ilv     (  nmltot)   ! nearest standard level above spreading point

  REAL    (KIND=wp   )     ::  &
    zrcutov (nmltot,3),& ! horizontal cut-off radius [km]
    fcnodiv (nmltot,3),& ! non-divergence correction factor for 2-D wind correl.
    rhtinfl (nmltot,3),& ! horizontal correlation scale
    rvinfl  (nmltot,3),& ! vertically dependent part of horizontal correl. scale
    fvnodiv (nmltot  )   ! vertically dependent part of weight to non-div. corr.

  LOGICAL                  , ALLOCATABLE , SAVE ::  &
    lsprsml (:)         ! node has part of area of influence of station 'ista'
!
!------------ End of header ----------------------------------------------------


!     THRESHOLD '1' FOR VERTICAL THETA(E) GRADIENT : + 2K /.122
!       EQUALS + 2K / LN(1.130) , WHICH IS + 2K / 1000M FOR TV = 280K
!      THRESH1 = + 2./.122
!      THRESH2 =   0./.122
!     FOR SPREADING ALONG ISOBARS TO BE PERFORMED IN *LEVSPRD* AND
!     *HSPREAD* : SET 'THRESH2' SO HIGH (20K / 10M) THAT VERTICAL
!       GRADIENT OF SPREADING PARAMETER WILL ALWAYS BE LESS
!     IF ((MSPRPAR.EQ.2).OR.(MSPRPAR.EQ.3)) THRESH2 = 2000./.122
!     THRESH1 = AMAX1( THRESH1 , THRESH2 + 1.E-8 )
!  -->   COMPARE 'CARR, MWR NOV. 93, P. 3109; BARWELL AND LORENC, 1985:'
!          CORRELATION = EXP( -3.*(LN(P1/P2))**2 ) , FOR WIND
!  -->   COMPARE 'RABIER, MCNALLY, ECMWF TECH. MEMO. 195, SEPT. 93:'
!          CORRELATION FOR T IS LESS, AND IS LESS CLOSER TO THE GROUND
!     ALODPS (IVAR) = .333
!  -->   COMPARE 'BENJAMIN, MWR JULY 89, P. 1597:'
!          CORRELATION = .32 + .68*(1.+.19*DTH) *EXP(-.19*DTH)
!          MY CHOICE AS FOR NOW:   EXP((-.19*DTH)**2)
!     DSPRDS (IVAR) = (1./.19)**2

 
!-------------------------------------------------------------------------------
! Begin Subroutine mult_org_spread
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1: Exclude all multi-level reports / data, for which no grid pts. of
!             the current node and vertical and temporal model level exist
!             within the 4-dim area of influence.
!             (Since upper limits for horizontal radii of influence are used,
!              these checks will find most, but not all data that have no in-
!              fluence. The temporal check implicitly includes quality control.)
!-------------------------------------------------------------------------------

  modespr = msprpar

  acthr  = ntstep * dtdeh
! lprsys = (lfirst) .OR. (ntstep <= 2)                                         &
!                   .OR.        (rmod(  ntstep   *dt,c3600 ) < tconbox)
! lprthr = (lfirst) .OR. (      (rmod( (ntstep+1)*dt,c3600 ) < tconbox)        &
!                         .AND. (ntstep*dtdeh <=                               &
!                                NINT(MIN(nudgend,nstop)*dtdeh)-0.9_wp))

  lwobs = lwonl
  IF ((nmltot == 0) .AND. ((lfirst) .OR. (rmod( ntstep*dt,c3600 ) < tconbox))  &
                    .AND. (lcroot) .AND. (k == ke))                            &
    PRINT '(''hour '',F3.0,'' : no multi-level data'')' , acthr


! allocate memory of fields used only in procedures called by current procedure
! -----------------------------------------------------------------------------

  ALLOCATE ( zoic1 (istart:iend,jstart:jend)   , STAT=istat )
  ALLOCATE ( zoiv1 (istart:iend,jstart:jend)   , STAT=istat )
  ALLOCATE ( omys1 (istart:iend,jstart:jend)   , STAT=istat )
  ALLOCATE ( klva  (istart:iend,jstart:jend,2) , STAT=istat )
  ALLOCATE ( klvb  (istart:iend,jstart:jend,2) , STAT=istat )

  DO j   = jstart , jend
    DO i = istart , iend
      zoic1 (i,j)   = c0
      zoiv1 (i,j)   = c0
      omys1 (i,j)   = c0
      klva  (i,j,1) = 0
      klva  (i,j,2) = 0
      klvb  (i,j,1) = 0
      klvb  (i,j,2) = 0
    ENDDO
  ENDDO

  IF (k == ke) THEN

    IF (ALLOCATED( lsprsml )) THEN
      PRINT '("CAUTION in src_mult_spread: lsprsml is already allocated "      &
            &,"at time ",I6)', ntstep
      DEALLOCATE ( lsprsml  , STAT=istat )
    ENDIF
    ALLOCATE ( lsprsml (nmltot) , STAT=istat )

! check if this sub-domain may be influenced by the current report (station)
! (The check includes horizontal check, temporal check,
!  and check of the nudging coefficients.)
! --------------------------------------------------------------------------

    DO ivrs = 1 , 3
      ivar = ivrs + MIN( ivrs-1 , 1 )
      DO istaml = 1 , nmltot
        ista = isortml(istaml)
        kobtyp = kobtyml(ista)
        kcdtyp = kcdtyml(ista)
        lraso  = (kobtyp == ntemp) .OR. (kobtyp == npilot)
        zrcutov (istaml,ivrs) = -c1
        IF (      (MAX( zwtml(ista,1,ivrs) , zwtml(ista,2,ivrs)                &
                      , zwtml(ista,3,ivrs) , zwtml(ista,4,ivrs) ) > epsy)      &
            .AND. (     ((lraso           ) .AND. (gnudg  (ivar) > epsy))      &
                   .OR. ((kobtyp == ngps)   .AND. (gnudggp       > epsy)       &
                                            .AND. (        ivar == 4   ))      &
                   .OR. (     (kobtyp == nairep)                               &
                        .AND. (kcdtyp /= nmodes) .AND. (gnudgar(ivar) > epsy)) &
                   .OR. (     (kobtyp == nairep)                               &
                        .AND. (kcdtyp == nmodes) .AND. (gnudgms(ivar) > epsy)) &
                   .OR. ((kobtyp == nsattv) .AND. (gnudgtv(ivar) > epsy))))    &
          zrcutov (istaml,ivrs) = zriflml(ista,2,ivrs) * cutofr(ivar)
      ENDDO
    ENDDO

    DO istaml = 1 , nmltot
      ista = isortml(istaml)
      io  = ioml (ista)
      jo  = joml (ista)
!     kobtyp = kobtyml(ista)
      lsprsml (ista) = .FALSE.
      zrcutof = MAX( zrcutov(istaml,1) , zrcutov(istaml,2) , zrcutov(istaml,3) )
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
      larea             = (      (jstaspr <= jend) .AND. (jendspr >= jstart)   &
                           .AND. (istaspr <= iend) .AND. (iendspr >= istart)   &
                           .AND. (zrcutof > epsy) )
      lsprsv (istaml,1) = (larea) .AND. (zrcutov(istaml,1) > epsy)
      lsprsv (istaml,2) = (larea) .AND. (zrcutov(istaml,2) > epsy)
      lsprsv (istaml,3) = (larea) .AND. (zrcutov(istaml,3) > epsy)
    ENDDO

! check for which variables this sub-domain may
! be influenced by the current report (station)
! ---------------------------------------------

    nactio = 6
    DO ivrs = 1 , 3
      ivex = ivrs + MIN( ivrs-1 , 1 ) 
      DO istaml = 1 , nmltot
        IF (lsprsv(istaml,ivrs)) THEN 
          ista = isortml(istaml)
          io  = ioml (ista)
          jo  = joml (ista)
          jrange  = INT( zrcutov(istaml,ivrs) * gppkmj )
          jstaspr = MAX( INT( jo - jrange ,iintegers) , jrs_tot )
          jendspr = MIN( INT( jo + jrange ,iintegers) , jre_tot ) 
          irange  = INT( zrcutov(istaml,ivrs) *MAX( gppkmi(jstaspr)            &
                                                  , gppkmi(jendspr) ))
!         irange  = INT( zrcutov(istaml,ivrs) *MAX( gppkmi(MIN( jstaspr,jend)) &
!                                             , gppkmi(MAX( jendspr,jstart)) ) )
          ijspr (3,istaml) = MAX( INT( jo - jrange ,iintegers) , jstart )
          ijspr (4,istaml) = MIN( INT( jo + jrange ,iintegers) , jend )
          ijspr (1,istaml) = MAX( INT( io - irange ,iintegers) , istart )
          ijspr (2,istaml) = MIN( INT( io + irange ,iintegers) , iend )
          lsprsv (istaml,ivrs)  = (      (ijspr(3,istaml) <= jend)             &
                                   .AND. (ijspr(4,istaml) >= jstart)           &
                                   .AND. (ijspr(1,istaml) <= iend)             &
                                   .AND. (ijspr(2,istaml) >= istart) )
        ENDIF
      ENDDO
      DO istaml = 1 , nmltot
        IF (lsprsv(istaml,ivrs)) THEN 
          ista = isortml(istaml)
          io  = ioml (ista)
          jo  = joml (ista)
          istaspr = ijspr(1,istaml)
          iendspr = ijspr(2,istaml)
          jstaspr = ijspr(3,istaml)
          jendspr = ijspr(4,istaml)
          rtinfl  = zrcutov(istaml,ivrs) 
!         nactio  = 6      ! set further above

          CALL cutoff_wind_correl ( nactio , k , c0 , ivex )
!         =======================

          lsprsv (istaml,ivrs) = (ivex > 0) 
          IF (lsprsv(istaml,ivrs))  lsprsml (ista) = .TRUE.
        ENDIF
      ENDDO
      DO istaml = 1 , nmltot
        IF (.NOT. lsprsv(istaml,ivrs)) THEN
          ista = isortml(istaml)
          kviflml (ista,1,ivrs) = 0
          kviflml (ista,2,ivrs) = ke + 1
          kviflml (ista,3,ivrs) = 0
          kviflml (ista,4,ivrs) = ke + 1
        ENDIF
      ENDDO
    ENDDO

!-------------------------------------------------------------------------------
!  Section 2: Get the vertical range of influenced model levels
!             (if spreading is horizontal)
!-------------------------------------------------------------------------------

    DO istaml = 1 , nmltot

      ista = isortml(istaml)
      io  = ioml (ista)
      jo  = joml (ista)

      DO ivrs = 1, 3
        IF (lsprsv(istaml,ivrs)) THEN 
          ivar = ivrs + MIN( ivrs-1 , 1 )
          DO nba = 4 , 1 , -1
            IF (kviflml(ista,nba,ivrs) < 0) THEN
              nbotop = 2 - MOD( nba , 2 )
              zrcutof = zriflml(ista,nbotop,ivrs) * cutofr(ivar)
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
                nactio = 3 + nbotop
                ivex   = ivar
  
                CALL cutoff_wind_correl ( nactio , k , c0 , ivex )
!               =======================
  
                izext = MOD( ivex , 10000 )
                jzext = ivex / 10000
              ENDIF
              IF ((izext == 0) .AND. (nbotop == 2)) THEN
                kviflml (ista,nba-1,ivrs) = 0
                kviflml (ista,nba  ,ivrs) = ke + 1
              ELSEIF (izext > 0) THEN
                z2cut = c2 * zvcutml(ista,nba,ivrs)
                cut: DO kk = 1 , ke + 1
                  kviflml (ista,nba,ivrs) = kk
                  IF (kk == ke+1)                                       EXIT cut
                  IF ( hhl(izext,jzext,kk)                                     &
                      +hhl(izext,jzext,kk+1) <= z2cut)                  EXIT cut
                ENDDO cut
                IF (nbotop == 1)                                               &
                  kviflml(ista,nba,ivrs) = kviflml(ista,nba,ivrs) - 1
              ELSE
                kviflml (ista,nba,ivrs) = ABS( kviflml(ista,nba,ivrs) ) - 1
              ENDIF
            ENDIF
          ENDDO
          IF ((lfirst) .AND. (lwobs) .AND. (kviflml(ista,1,ivrs) == 0)         &
                                     .AND. (kviflml(ista,3,ivrs) == 0))        &
            WRITE( nupr,'(''mul-kl: '',A ,I2, '', p,z: obs/v-cut:''            &
                        &,8F7.0,'': no k'')' ) ystidml(ista) (1:ilstidp), ivrs &
                 ,  zpobml (ista,1,ivrs)  , zpobml (ista,3,ivrs)               &
                 ,  zspobml(ista,2,ivrs)  , zspobml(ista,4,ivrs)               &
                 , (zvcutml(ista,nba,ivrs), nba=1,4)
!                , (kviflml(ista,nba,ivrs), nba=1,2)
        ENDIF
      ENDDO
    ENDDO

  ENDIF   !  (k == ke)


!-------------------------------------------------------------------------------
!  Section 3: Organize the spreading of observation increments
!-------------------------------------------------------------------------------


! Get the horizontal correlation scales
! and the non-divergent correction factor to the 2-dim. wind correlations
! -----------------------------------------------------------------------

  DO istaml = 1 , nmltot
    ilv (istaml) = 0
    ista = isortml(istaml)
    IF (lsprsml(ista)) THEN
      IF (      (k <= MAX( kviflml(ista,1,1) , kviflml(ista,3,1)               &
                         , kviflml(ista,1,2) , kviflml(ista,3,2)               &
                         , kviflml(ista,1,3) , kviflml(ista,3,3) ))            &
          .AND. (k >= MIN( kviflml(ista,2,1) , kviflml(ista,4,1)               &
                         , kviflml(ista,2,2) , kviflml(ista,4,2)               &
                         , kviflml(ista,2,3) , kviflml(ista,4,3) ))) THEN
        IF (zlopml(ista,k) <= tabcolp(ncolev)) THEN
          ilv (istaml) = ncolev
        ELSE
          ilv (istaml) = 2
          DO WHILE (zlopml(ista,k) <= tabcolp(ilv(istaml)))
            ilv (istaml) = ilv(istaml) + 1
          ENDDO
        ENDIF
      ENDIF
    ENDIF
  ENDDO

  DO istaml = 1 , nmltot
    ista = isortml(istaml)
    IF ((lsprsml(ista)) .AND. (ilv(istaml) >= 1)) THEN
      ilva = ilv(istaml)
      ilvb = ilva - 1
      zlopf = MIN( c1 , MAX( c0 , (tabcolp(ilvb) -zlopml(ista,k))              &
                                / (tabcolp(ilvb) -tabcolp(ilva)) ) )
      fvnodiv (istaml)   = (c1-zlopf) *fnodivc(ilvb) + zlopf *fnodivc(ilva)
      rvinfl  (istaml,1) = (c1-zlopf) *rhvsond(ilvb) + zlopf *rhvsond(ilva)
      rvinfl  (istaml,2) = (c1-zlopf) *rhtsond(ilvb) + zlopf *rhtsond(ilva)
      rvinfl  (istaml,3) = (c1-zlopf) *rhqsond(ilvb) + zlopf *rhqsond(ilva)
    ENDIF
  ENDDO

  DO ivrs = 1, 3
    ivar = ivrs + MIN( ivrs-1 , 1 )
    DO istaml = 1 , nmltot
      ista = isortml(istaml)
      kviflbot = MAX( kviflml(ista,1,ivrs) , kviflml(ista,3,ivrs) )
      kvifltop = MIN( kviflml(ista,2,ivrs) , kviflml(ista,4,ivrs) )
      IF ((lsprsml(ista)) .AND. (k <= kviflbot) .AND. (k >= kvifltop)) THEN
        zrwt = MAX( zwtml(ista,1,ivrs) , zwtml(ista,2,ivrs) )
!       IF (      (.NOT. ltiml(ista,ivrs))                                     &
!           .AND. (MIN( zwtml(ista,1,ivrs) , zwtml(ista,2,ivrs) ) > epsy))     &
!         zrwt = MIN( zwtml(ista,1,ivrs) , zwtml(ista,2,ivrs) )
        fcnodiv(istaml,ivrs)  = c0
        IF (ivrs == 1)                                                         &
          fcnodiv (istaml,ivrs) = MIN( c1, (c1 + (tnondiv-c1) *(c1-zrwt))      &
                                          *(cnondiv + fnondiv *fvnodiv(istaml)))
        rhtinfl (istaml,ivrs) = MAX( epsy,                                     &
                           (c1 + (rhtfac(ivar)-c1) *(c1-zrwt))                 &
                          *(rhinfl(ivar) + rhvfac(ivar)  *rvinfl(istaml,ivrs)) )
        IF (kobtyml(ista) == ngps)     rhtinfl (istaml,ivrs) =                 &
                                       rhtinfl (istaml,ivrs) * rhfgps
        IF (kobtyml(ista) == nsattv)   rhtinfl (istaml,ivrs) =                 &
                                       rhtinfl (istaml,ivrs) * rhfrtv(ivar)
      ENDIF
    ENDDO
  ENDDO


! start loop over multi-level reporting stations
! ----------------------------------------------

loop_over_multi_level_stations:  DO istaml = 1 , nmltot
 
  ista = isortml(istaml)

  io  = ioml (ista)
  jo  = joml (ista)

  kobtyp = kobtyml(ista)
  kcdtyp = kcdtyml(ista)
  lraso  = (kobtyp == ntemp) .OR. (kobtyp == npilot)

  ityp = isetyp0
  DO nba = SUM( niwtyp ) , 1 , -1
    IF (     ((iwtyp(nba) > 0) .AND. ( iwtyp(nba) == kobtyp))                  &
        .OR. ((iwtyp(nba) < 0) .AND. (-iwtyp(nba) == kcdtyp)))                 &
      ityp = isetyp(nba)
  ENDDO

  IF (lsprsml(ista)) THEN

    io_tot = ioml_tot(ista)
    jo_tot = joml_tot(ista)
    ystid  = ystidml (ista) (1:ilstidp)

!   lwobs  =       (io  <= iend) .AND. (io  >= istart)                         &
!            .AND. (jo  <= jend) .AND. (jo  >= jstart)
    lwobs  = lwonl

! flush YUPRINT file
!   IF (ldump_ascii) THEN
!     IF (      (lwonl) .AND. ((lfirst) .OR. (ntstep <= 1))                    &
!         .AND. (     ((io == ionl) .AND. (jo == jonl) .AND. (k >= ke-6))      &
!                .OR. (k == ke))) THEN
!       WRITE( nupr,'(''mul-pr: '',A )' )  ystid
!       CLOSE (nupr)
!       OPEN  (nupr,FILE=yuprint,FORM='FORMATTED',STATUS='OLD'                 &
!                               ,POSITION='APPEND',IOSTAT=nstat)
!       IF (nstat /= 0) PRINT '("OPENING OF FILE yuprint FAILS, mult_org_spread")'
!     ENDIF
!   ENDIF

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

      kviflbot = MAX( kviflml(ista,1,ivrs) , kviflml(ista,3,ivrs) )
      kvifltop = MIN( kviflml(ista,2,ivrs) , kviflml(ista,4,ivrs) )

      IF ((k <= kviflbot) .AND. (k >= kvifltop)) THEN

! get the scaled horizontal distance to the station location, the geometrical
! factors of the 2-dim horizontal wind correlations, and the area of influence
! ----------------------------------------------------------------------------

        rtinfl  = rhtinfl(istaml,ivrs)
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
        nactio  = 1
        IF (     (jstaspr > jend) .OR. (jendspr < jstart)                      &
            .OR. (istaspr > iend) .OR. (iendspr < istart) )  nactio  =  0

        CALL cutoff_wind_correl ( nactio , k , fcnodiv(istaml,ivrs) , ivar )
!       =======================

        IF (modespr == 2)                                                      &
! note: 'kviflml(ista,2or4,ivrs)' can be set to k+1 within *cutoff_wind_correl*
          kvifltop = MIN( kviflml(ista,2,ivrs) , kviflml(ista,4,ivrs) )

! spreading of observation increments along horizontal or isentropic surfaces
! ---------------------------------------------------------------------------

        IF ((k > ksprml(ista)) .AND. (k >= kvifltop) .AND. (nactio > 0)) THEN
          lsprx = .FALSE.

          loop_over_time_index:  DO itim = 1, 2

            ivrc  = ivrs
            IF ((itim == 2) .AND. (lsprx)) ivrc = ivrs + 3
!           lsprx =       (MAX( zwtml(ista,itim  ,ivrs)                        &
!                             , zwtml(ista,itim+2,ivrs) ) > epsy)              &
            lsprx =       (k <= kviflml(ista,2*itim-1,ivrs))                   &
                    .AND. (k >= kviflml(ista,2*itim  ,ivrs))

            IF (lsprx) THEN
              modex = modespr

              IF (modex == 2) CALL mult_find_levels ( ivrc , k , modex )
!                             =====================

! note: 'modex' can be set to 1 within *mult_find_levels*

              IF (modex == 1) CALL mult_find_levels ( ivrc , k , modex )
!                             =====================

            ENDIF

            IF (ivrs == 1) CALL mult_spread_wind ( lsprx , ityp , k )
!                          =====================

            IF (ivrs >= 2) CALL mult_spread_mass ( lsprx , ivar , ityp , k )
!                          =====================

          ENDDO loop_over_time_index

        ELSEIF ((k >= kvifltop) .AND. (nactio > 0)) THEN

! spreading of observation increments along model levels
! ------------------------------------------------------

          CALL mult_spread_modlev ( ivar , ityp , k )
!         =======================

        ENDIF

! printout for control
! --------------------

        IF ((lfirst) .AND. (lwonl) .AND. (k == ke-5)) THEN
          IF (ivrs == 1) WRITE(nupr,'(''zuwi,om?u: '',A,F7.4,4X,2F7.4)') ystid &
             ,zwi(ionl,jonl,1,ityp), omy(ionl,jonl,1,ityp),om2(ionl,jonl,1,ityp)
          IF (ivrs == 2) WRITE(nupr,'(''ztwi,om?t: '',A,F7.4,4X,2F7.4)') ystid &
             ,zwi(ionl,jonl,3,ityp), omy(ionl,jonl,2,ityp),om2(ionl,jonl,2,ityp)
          IF (ivrs == 3) WRITE(nupr,'(''zqwi,om?q: '',A,F11.8,  2F7.4)') ystid &
             ,zwi(ionl,jonl,4,ityp), omy(ionl,jonl,3,ityp),om2(ionl,jonl,3,ityp)
        ENDIF
        IF (lfirst) THEN
!         IF ((k == (kviflbot+kvifltop)/2) .AND. (ivrs == 1) .AND. (lwobs)) THEN
          IF ((     (k == kviflbot)                                            &
               .OR. (k <= kvifltop)) .AND. (ivrs == 1) .AND. (lwobs)) THEN
            nba  = 2
            nbah = 4
            IF (kviflml(ista,1,ivrs) < 1) nba  = 3
            IF (kviflml(ista,2,ivrs) < 1) nba  = 1
            IF (kviflml(ista,2,ivrs) < 1) nbah = 2
            IF (modespr <= 1) THEN
              WRITE( nupr,'(''mul-infl(km)'',F5.0, 1X,A ,2I4                   &
                          &,'', z: obs,v-cut:'',4F7.0 ,'', k-range:'',2I3)' )  &
                     zrcutof, ystid, io_tot, jo_tot                            &
                   , zspobml(ista,nba,ivrs)  , zspobml(ista,nbah,ivrs)         &
                   , zvcutml(ista,nba,ivrs)  , zvcutml(ista,nbah,ivrs)         &
                   , kviflbot, kvifltop
!                  , (mszlev(ista,nba), nba=1,2)
            ELSE
              WRITE( nupr,'(''mul-infl(km)'',F5.0, 1X,A ,2I4                   &
                          &,'', z: obs,v-cut:'',4F7.2 ,'', k-range:'',2I3)' )  &
                     zrcutof, ystid, io_tot, jo_tot                            &
                   , zspobml(ista,nba,ivrs)  , zspobml(ista,nbah,ivrs)         &
                   , zvcutml(ista,nba,ivrs)  , zvcutml(ista,nbah,ivrs)         &
                   , kviflbot, kvifltop
            ENDIF
          ENDIF
        ENDIF

      ENDIF

    ENDDO loop_over_variables

!   ENDDO loop_over_time_index

  ENDIF ! lsprsml(ista)

ENDDO loop_over_multi_level_stations


! deallocate memory
! -----------------

  DEALLOCATE ( zoic1 , STAT=istat )
  DEALLOCATE ( zoiv1 , STAT=istat )
  DEALLOCATE ( omys1 , STAT=istat )
  DEALLOCATE ( klva  , STAT=istat )
  DEALLOCATE ( klvb  , STAT=istat )

  IF (k == 1) DEALLOCATE ( lsprsml , STAT=istat )

! flush YUPRINT file
  IF (ldump_ascii) THEN
    IF ((lwonl) .AND. ((lfirst) .OR. (ntstep <= 1))                            &
                .AND. ((k == ke) .OR. (k == 1) .OR. (k >= ke-6))) THEN
      CLOSE (nupr)
      OPEN  (nupr,FILE=yuprint,FORM='FORMATTED',STATUS='OLD'                   &
                              ,POSITION='APPEND',IOSTAT=nstat)
      IF (nstat /= 0) PRINT '("OPENING OF FILE yuprint FAILED, mult_org_spread")'
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure mult_org_spread
!-------------------------------------------------------------------------------

END SUBROUTINE mult_org_spread


!-------------------------------------------------------------------------------
!+ Module procedure for spreading scalar multi-level obs increments
!-------------------------------------------------------------------------------

SUBROUTINE mult_spread_mass ( lsprx , ivar , ityp , k )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure performs the weighting and spreading (extrapolation)
!   of scalar observation increments (potential temperature, or specific
!   humidity) from (at most 2) multi-level reports from one observing station.
!   The lateral spreading takes place along surfaces (horizontal or isentropic)
!   which differ from the model surfaces.
!
! Method:
!   Spreading of multi-level data increments to the target grid pts. with the
!   largest weighting along specified surfaces according to 'modespr' (msprpar).
!   Computation of the sum of 'relative nudging weights' by applying the
!   temporal and spatial correlations and quality factors.
!   Option for non-isotropic lateral correlations (see below).
!   Humidity: Computation of specific humidity ('qd') increments from the
!             spreaded relative humidity ('RH') increments, so that the nudging
!             pulls towards observed 'RH' instead of observed 'qd'.
!   Mode of spreading:   'modespr' =
!     1: 'spreading along horizontal levels' :
!        - vertical correlation as Gaussian of height differences ('dz') between
!          the obs. point and the target grid pt. at level 'k'
!        - non-isotropy as Gaussian of pot. temp. differences on horiz. levels
!     2: 'spreading along isentropic surfaces' :
!        - vertical correlation as Gaussian of potential temperature differences
!          between the obs. point and the target grid pt. at level 'k'
!        - non-isotropy as Gaussian of height differences on isentropic surfaces
!        Optionally, spreading may be along model levels for target grid points
!        above a specified model level (see proc. *mult_spread_modlev*). At the
!        level below that level, or in the presence of weak thermal stability
!        at the target grid pt., spreading along isentropic surfaces is blended
!        with horizontal spreading for which in this case, however, vertical
!        correlations are still a Gaussian of potential temperature differences,
!        and non-isotropic weights are not possible.
!     correlation functions: - in z     : EXP( -(g(z1-z2)/r*Tv1)**2 / dzscal) )
!                            - in theta : EXP( -(theta1-theta2)**2 / dthscal) )
!     For target grid points within the vertical extent of the observation
!     increment profile, the weights are reduced in the following way:
!     1. The vertical distance to the nearest upper and lower observation is
!        assigned to be the sum of the 'distance' (as in '2.', see below) of the
!        observation increment level to its closest observation, and the (verti-
!        cal) distance of the obs. increment level to the target grid point.
!     2. An effective vertical distance is computed from the 2 (i.e. upper and
!        lower) distances in the same way as the resistance resulting from 2
!        parallel resistances:   d_eff = d_upper * d_lower / (d_upper + d_lower)
!     3. Weight function (like correlations):   EXP( -(d_eff)**2 / d?scal )
!
! Written by        :  Christoph Schraff, DWD  (original version: 10.09.97)
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
    lsprx               ! a report (i.e. obs. increments) exists with indices =
                        ! (ista,itim) and with non-zero temporal and vertical
                        ! weights for at least one datum in current report
                        ! ==> there are obs. increments to be spreaded

  INTEGER (KIND=iintegers), INTENT (IN)         ::       &
    ivar             ,& ! index of spreaded scalar quantity: 3=T; 4=RH
    ityp             ,& ! index of weighted increment array
    k                   ! index of current vertical model level


! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    modespr          ,& ! mode of spreading, see above: 'Method'
    msi    , msip    ,& ! index of (effective) spreading parameter
    mni              ,& ! index of parameter defining non-isotropic weights
    meqi             ,& ! index of parameter defining equidistant points
    ivrs             ,& ! index of spreaded scalar quantity: 2=T; 3=RH
    irep             ,& ! index of report per station for '?cutof'
    isvv             ,& ! index of spreaded quantity depend. on spreading param.
!   nmloi            ,& ! index of report in the observation increment field
    noix             ,& ! index of spreaded quantity in obs. increment field
    nois   , noic    ,& ! index of (selected) spreading param. in obs incr field
    noin             ,& ! index of non-isotropic weights para. in obs incr field
    noiqc            ,& ! index of quality factor in obs. increment field
    noivlr           ,& ! index of vertical weight given to the obs. incr. due
                        ! to the distance to the obs. at the vertical interpol.
    nba1   , nba2    ,& ! temporal and vertical indices
    ixc              ,& ! loop index over influenced grid points
    i      , j       ,& ! loop indices in horizontal direction
!   io     , jo      ,& ! local indices of location of observation
!   istaspr, jstaspr ,& ! lower left corner of domain containing area of infl.
!   iendspr, jendspr ,& ! upper right corner of domain containing area of infl
    nprae            ,& ! number of printed grid pts.
    jpra   , jpre    ,& ! lowest , highest meridional index of printed grid pts.
    mvb              ,& ! index of equidistant point
    joc                 ! loop index for printout

  REAL    (KIND=wp   )     ::  &
    fmultw           ,& ! 0. or 1., depending on multiple observation weighting
    rdsprs           ,& ! vertical correlation scale as by namelist
    rdsprsb, rdsprst ,& ! lower / upper vertical correlation scale
    rdsprsrb         ,& ! /  reciproke of square of lower /
    rdsprsrt         ,& ! \  upper vertical correlation scale
    rdsprir          ,& ! reciproke of vert. correl. scale used for interpolat.
    rdsprni          ,& ! square of Gauss. non-isotropy radius,'zspr(mni)'-units
    rdsprnir         ,& ! reciproke of 'rdsprni'
    zcutni           ,& ! 'rdsprni' < 'zcutni' : correlation is non-isotropic
    ztp              ,& ! value of spreading param. (z, TH) at (std atm) tropop.
    zsproni          ,& ! non-isotr. w.: TH or z at spreading surface & obs loc.
    zqts             ,& ! quotient of the tropospheric to the stratospheric mean
                        ! vertical gradient of the spreading parameter
    zfvip            ,& ! quotient for vertical interpolation
    zsvob            ,& ! spreading parameter at the obs. point
    zwvip            ,& ! -log of weight reduction due to vertical interpolation
    zwvip2           ,& ! vertical distance of obs. incr. level to target g.pt.
    zfvipni          ,& ! quotient for vertical interpol. from equidistant pts.
    omywth , omywp   ,& ! weighted nudging weight for a single obs. spreaded
                        ! along isentropes, resp. horizontally
    dthreshr         ,& ! diff. of stability thresholds defining the fraction of
                        ! spreading along isentropes (reciproke of item)
    sprmx            ,& ! max. allowed fraction of spreading along isentropes
    zdds2            ,& ! squared scaled horizontal distance to the obs. station
    omycorl          ,& ! local vertical correlation
    zqgnudg          ,& ! nudging coeff. relative to the coeff. for TEMP / PILOT
    qmyk1  , qmyk2   ,& ! fractions of temporal linear interpolation
    zispr  , disprr     ! define vertic. equidist. pts, for non-isotropic weight

  LOGICAL                  ::  &
    lnoniso             ! non-isotrophic horizontal correlation

  CHARACTER (LEN=15) yformat


  LOGICAL                  , SAVE ::  &
    lsprx1              ! non-zero temporal and vertical weights from 1st report


! Local (automatic) arrays:
! -------------------------
  REAL    (KIND=wp   )     ::  &
    zoics   (istart:iend,jstart:jend) ,& ! scalar obs. increment at target grid
                                         ! pt. from 1 report
    omysc   (istart:iend,jstart:jend) ,& ! nudging weight for a single obs.
    zvlop   (istart:iend,jstart:jend) ,& ! log(pressure) at the obs. incr. point
    zoizs   (istart:iend,jstart:jend) ,& ! scalar obs. increment at the target
                                         ! grid pt. using horizontal spreading
    omyzc   (istart:iend,jstart:jend) ,& ! nudging weight for a single obs.
                                         ! spread horizontally
    sprdv   (istart:iend,jstart:jend) ,& ! fraction of obs incr. to be spread
                                         ! along isentropes
    zsoem   (istart:iend,jstart:jend) ,& ! sum of the 2 spreaded obs. increments
    omykm   (istart:iend,jstart:jend) ,& ! sum of the 2 weights
    omy12   (istart:iend,jstart:jend) ,& ! square of the weight for first obs.
    omyc2   (istart:iend,jstart:jend)    ! square of the weight for second obs.

! Initialized variables:
! ----------------------

  LOGICAL                  , SAVE ::  &
    lfrsloc = .TRUE.    ! .TRUE if local output is to be printed the first time

!------------ End of header ----------------------------------------------------


 
!-------------------------------------------------------------------------------
! Begin Subroutine mult_spread_mass
!-------------------------------------------------------------------------------
  
!-------------------------------------------------------------------------------
!  Section 1: Preliminaries
!-------------------------------------------------------------------------------

  lnoniso = .FALSE.    ! initialize in any case for safety

  IF (itim == 1) lsprx1 = lsprx

  ivrs    = MAX( ivar - 1 , 1 )
  modespr = msprpar
! nmloi   = knoiml(ista,itim)

  irep    = 1
  IF ((itim == 2) .AND. (lsprx1)) irep = 2

  fmultw = REAL ( kwtyp(ityp) - 1, wp )

! lsigma =     ( modespr == 0)                                                 &
!         .OR. ((modespr == 1) .AND. (k <= ktp))                               &
!         .OR. ((modespr == 2) .AND. (k < ktth))
!!        .OR. ((modespr == 2) .AND. (     ((.NOT. ltopth) .AND. (k < ktopth)) &
!!                                    .OR. ((ltopth) .AND. (k <= ktp))))

  IF (lsprx) THEN

    nba1 = 2*itim - 1
    nba2 = 2*itim

! get the spreading parameter at standard atmospheric tropopause if required,
! and the Gaussian radii of the vertical weight function

    ztp     = zstpml (ista)
    rdsprs  = SQRT( vcorls(ivar) )
    IF (modespr <= 1) THEN
      rdsprsb = rdsprs * fcorlml(ista,nba1,ivrs) * zrtdgml(ista,nba1,ivrs)
      rdsprst = rdsprs * fcorlml(ista,nba2,ivrs) * zrtdgml(ista,nba2,ivrs)
      rdsprir = c1 / SQRT( MAX( vcorls(ivar) , 0.04_wp ) )
    ELSE
      rdsprsb = rdsprs * fcorlml(ista,nba1,ivrs)
      rdsprst = rdsprs * fcorlml(ista,nba2,ivrs)
      rdsprir = c1 / SQRT( MAX( vcorls(ivar) , 33.0_wp ) )
    ENDIF
    rdsprsrb = c1 / rdsprsb **2
    rdsprsrt = c1 / rdsprst **2

! (Horizontally) local contributions to the weights: temporal weights, (quality
! factors), relative nudging coeff., vertical weight (spread. along model lev.)

! wtql1  =  zwtua (ista,1,ivrs)  *  zqualua (ista,1,ivrs)                      &
!         * gnudgar(ivar) / MAX( gnudg(ivar) ,epsy )

! get quantities used for the vertical weighting and non-isotropic correction,
! and the relative nudging coefficient
 
    IF (modespr <= 1) THEN
      zispr  = xezispr
      disprr = c1 / dezispr
      zcutni = vcutnit
    ELSE
      zispr  = xthispr
      disprr = c1 / dthispr
      zcutni = vcutnip
    ENDIF

    rdsprni = vcsni  (ivar)
    lnoniso = (rdsprni <= zcutni) .AND. (lnisua)
    IF (modespr == 2) rdsprni = rdsprni * zrtgkml(ista,k) **2
    IF ((modespr <= 1) .AND. (zsprml(ista,k) > ztp)) rdsprni = qs2 * rdsprni
    rdsprnir = c1 / MAX( rdsprni , epsy )
    zqgnudg = c1
    IF (kobtyp == nairep) THEN
      IF (kcdtyp /= nmodes) zqgnudg = gnudgar(ivar) / MAX( gnudg(ivar) ,epsy )
      IF (kcdtyp == nmodes) zqgnudg = gnudgms(ivar) / MAX( gnudg(ivar) ,epsy )
    ENDIF
    IF (kobtyp == nsattv) zqgnudg = gnudgtv(ivar) / MAX( gnudg(ivar) ,epsy )
    IF (kobtyp == ngps  ) zqgnudg = gnudggp       / MAX( gnudg(ivar) ,epsy )
    IF ((kobtyp == ngps ) .AND. (ivar /= 4))  zqgnudg = c0
!   lti2 = (ltiml(ista,ivrs)) .AND. (omyk1 >= epsy) .AND. (omyk2 >= epsy)

! get indices for observation increment and other arrays

    msi  = MAX( modespr , i1 )
    mni  = 3 - msi
    meqi = mni + 1
!   isvv = 3*(msi-1) + ivrs
    isvv =    msi
    IF (ivrs == 2) THEN
      noix   = noit
      noiqc  = noitqc
      noivlr = noitlr
    ELSEIF (ivrs == 3) THEN
      noix   = noirh
      noiqc  = noiqqc
      noivlr = noiqlr
    ENDIF
    IF (modespr <= 1) THEN
      nois = noiz
      noin = noith
      zqts = c1
    ELSE
      nois = noith
      noin = noiz
      zqts = c1 / qst
    ENDIF

  ENDIF   ! lsprx

! initialize

  IF (ltiml(ista,ivrs)) THEN
    DO j   = jstaspr , jendspr
      DO i = istaspr , iendspr
        zoics (i,j) = c0
        omysc (i,j) = c0
      ENDDO
    ENDDO
  ENDIF

!-------------------------------------------------------------------------------
!  Section 2: Computation of the 'spreaded observation increment' at the target
!             grid pts. by vertical interpolation of the observation increments
!             at the obs. location and application of the horizontal correla-
!             tion. Computation of the nudging weights for single observations.
!-------------------------------------------------------------------------------

! case a: Horizontal correlation non-isotropic
! --------------------------------------------

  IF ((lnoniso) .AND. (lsprx)) THEN
!CDIR NODEP,VOVERTAKE,VOB
    DO ixc = 1 , ncutof(irep)
      i = icutof(ixc,irep)
      j = jcutof(ixc,irep)
! if increments are available for vertical interpolation
      IF ((klva(i,j,isvv) > -1) .AND. (klvb(i,j,isvv) > -1)) THEN
        zfvip       =   ( oiml (nois,klvb(i,j,isvv)) - zspr(i,j,msi) )         &
                      / ( oiml (nois,klvb(i,j,isvv))                           &
                         -oiml (nois,klva(i,j,isvv)) )
! vertical interpolation of observation increments
        zoics (i,j) =     zfvip     * oiml (noix ,klva(i,j,isvv))              &
                       + (c1-zfvip) * oiml (noix ,klvb(i,j,isvv))
! log(pressure) at obs. point for conveying temperature incr. to pot. temp. incr
        zvlop (i,j) =     zfvip     * oiml (noilp,klva(i,j,isvv))              &
                       + (c1-zfvip) * oiml (noilp,klvb(i,j,isvv))
! value of the parameter used for non-isotropic weights at the sounding station
! and at the height where the value of the spreading parameter equals that at
! the target grid point
        zsproni     =     zfvip     * oiml (noin ,klva(i,j,isvv))              &
                       + (c1-zfvip) * oiml (noin ,klvb(i,j,isvv))
! reduction of weight due to vertical interpolation within vertical profile
        zsvob       =  oiml (nois,klva(i,j,isvv))
        zwvip       =  (  MAX( MIN( zsvob,ztp ) - zspr(i,j,msi) , c0 )         &
                        + MAX( zsvob - MAX( zspr(i,j,msi),ztp ) , c0 )*zqts )  &
                      * oiml (noivcr,klva(i,j,isvv))                           &
                     +  oiml (noivlr,klva(i,j,isvv))
        zsvob       =  oiml (nois,klvb(i,j,isvv))
        zwvip2      =  (  MAX( MIN( zspr(i,j,msi),ztp ) - zsvob , c0 )         &
                        + MAX( zspr(i,j,msi) - MAX( zsvob,ztp ) , c0 )*zqts )  &
                      * oiml (noivcr,klvb(i,j,isvv))                           &
                     +  oiml (noivlr,klvb(i,j,isvv))
        zwvip       =  (zwvip * zwvip2 / MAX( zwvip+zwvip2,epsy ) *rdsprir) **2
! computation of the nudging weight
        omysc (i,j) =   (c1 + zdds(i,j,ivrs))                                  &
                      * EXP(- (zsproni -zspr(i,j,mni))**2 * rdsprnir           &
                            - zdds(i,j,ivrs)  - zwvip )                        &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * (  zfvip     * oiml (noiqc,klva(i,j,isvv))             &
                         +(c1-zfvip) * oiml (noiqc,klvb(i,j,isvv)) )
! no obs. increment available below 'zspr(i,j)'. Take the obs. increment at the
! lowest obs. point with data of type 'ivrs' and apply vertical correlation
      ELSEIF (klva(i,j,isvv) > -1) THEN
! obs. increment
        zoics (i,j) =  oiml (noix ,klva(i,j,isvv))
! log(pressure) at obs. point for conveying temperature incr. to pot. temp. incr
        zvlop (i,j) =  oiml (noilp,klva(i,j,isvv))
! value of the parameter used for non-isotropic weights at the sounding station
! and at the height where the value of the spreading parameter equals that at
! the target grid point (approximate computation by vertical interpolation from
! equidistant point (in terms of spreading parameter units)) 
        mvb         =  INT( (zspr(i,j,meqi) - zispr) * disprr ) + 1
        mvb         =  MIN( MAX( mvb , i1 ) , mxispr - 1 )
        zfvipni     =    (zsdnis(mvb,msi) - zspr(i,j,msi)) * zsdnid(mvb,msi)
        zsproni     =    (c1-zfvipni) *znismq(ista,mvb)                        &
                       +  zfvipni     *znismq(ista,mvb+1)
! computation of the nudging weight:
!   spreading parameter at the obs. point
        zsvob       =  oiml (nois,klva(i,j,isvv))
!   part of horizontal correlation
        omysc (i,j) =   (c1 + zdds(i,j,ivrs))                                  &
!   vertical correlation: pot. temperature differences are scaled by a larger
!   constant in the stratosphere than in the troposphere (larger mean stability)
                      * EXP(- (  MAX( MIN( zsvob,ztp ) - zspr(i,j,msi) , c0 )  &
                               + MAX( zsvob - MAX( zspr(i,j,msi),ztp ) , c0 )  &
                                 * zqts ) **2  *  rdsprsrb                     &
!   further weighting for non-isotropy
                            - (zsproni -zspr(i,j,mni))**2 * rdsprnir           &
!   factor of horizontal correlation
                            - zdds(i,j,ivrs)  )                                &
!   reduction of weight within atmospheric boundary layer; relat. nudging coeff.
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
!   quality factor
                      * oiml (noiqc,klva(i,j,isvv))
! no obs. increment available above 'zspr(i,j)'
      ELSEIF (klvb(i,j,isvv) > -1) THEN
        zoics (i,j) =  oiml (noix ,klvb(i,j,isvv))
        zvlop (i,j) =  oiml (noilp,klvb(i,j,isvv))
        mvb         =  INT( (zspr(i,j,meqi) - zispr) * disprr ) + 1
        mvb         =  MIN( MAX( mvb , i1 ) , mxispr - 1 )
        zfvipni     =    (zsdnis(mvb,msi) - zspr(i,j,msi)) * zsdnid(mvb,msi)
        zsproni     =    (c1-zfvipni) *znismq(ista,mvb)                        &
                       +  zfvipni     *znismq(ista,mvb+1)
        zsvob       =  oiml (nois,klvb(i,j,isvv))
        omysc (i,j) =   (c1 + zdds(i,j,ivrs))                                  &
                      * EXP(- (  MAX( MIN( zspr(i,j,msi),ztp ) - zsvob , c0 )  &
                               + MAX( zspr(i,j,msi) - MAX( zsvob,ztp ) , c0 )  &
                                 * zqts ) **2  *  rdsprsrt                     &
                            - (zsproni -zspr(i,j,mni))**2 * rdsprnir           &
                            - zdds(i,j,ivrs)  )                                &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * oiml(noiqc,klvb(i,j,isvv))
! no spreading along isentropes where vertical gradient of potential temperature
! is less than a threshold
      ELSE
        zoics (i,j) =  c0
        zvlop (i,j) =  c0
        omysc (i,j) =  c0
      ENDIF
    ENDDO
! conveying temperature increments to potential temperature increments
    IF (ivrs == 2) THEN
!CDIR NODEP,VOVERTAKE,VOB
      DO ixc = 1 , ncutof(irep)
        i = icutof(ixc,irep)
        j = jcutof(ixc,irep)
        zoics (i,j) = zoics(i,j) * zeklop(i,j) * EXP( -rdocp *zvlop(i,j) )
      ENDDO
    ENDIF
  ENDIF


! case b: Spreading along isentropes
! ----------------------------------

  IF ((.NOT. lnoniso) .AND. (modespr == 2) .AND. (lsprx)) THEN
!CDIR NODEP,VOVERTAKE,VOB
    DO ixc = 1 , ncutof(irep)
      i = icutof(ixc,irep)
      j = jcutof(ixc,irep)
      IF ((klva(i,j,isvv) > -1) .AND. (klvb(i,j,isvv) > -1)) THEN
        zfvip       =   ( oiml (nois,klvb(i,j,isvv)) - zspr(i,j,msi) )         &
                      / ( oiml (nois,klvb(i,j,isvv))                           &
                         -oiml (nois,klva(i,j,isvv)) )
! vertical interpolation of observation increments
        zoics (i,j) =     zfvip     * oiml (noix ,klva(i,j,isvv))              &
                       + (c1-zfvip) * oiml (noix ,klvb(i,j,isvv))
! log(pressure) at obs. point for conveying temperature incr. to pot. temp. incr
        zvlop (i,j) =     zfvip     * oiml (noilp,klva(i,j,isvv))              &
                       + (c1-zfvip) * oiml (noilp,klvb(i,j,isvv))
! reduction of weight due to vertical interpolation within vertical profile
        zsvob       =   oiml (nois,klva(i,j,isvv))
        zwvip       =  (  MAX( MIN( zsvob,ztp ) - zspr(i,j,msi) , c0 )         &
                        + MAX( zsvob - MAX( zspr(i,j,msi),ztp ) , c0 )*zqts )  &
                      * oiml (noivcr,klva(i,j,isvv))                           &
                     +  oiml (noivlr,klva(i,j,isvv))
        zsvob       =   oiml (nois,klvb(i,j,isvv))
        zwvip2      =  (  MAX( MIN( zspr(i,j,msi),ztp ) - zsvob , c0 )         &
                        + MAX( zspr(i,j,msi) - MAX( zsvob,ztp ) , c0 )*zqts )  &
                      * oiml (noivcr,klvb(i,j,isvv))                           &
                     +  oiml (noivlr,klvb(i,j,isvv))
        zwvip       =  (zwvip * zwvip2 / MAX( zwvip+zwvip2,epsy ) *rdsprir) **2
! computation of the nudging weight
        omysc (i,j) =   (c1 + zdds(i,j,ivrs))                                  &
                      * EXP(- zdds(i,j,ivrs)  - zwvip )                        &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * (  zfvip     * oiml (noiqc,klva(i,j,isvv))             &
                         +(c1-zfvip) * oiml (noiqc,klvb(i,j,isvv)) )
      ELSEIF (klva(i,j,isvv) > -1) THEN
        zoics (i,j) =   oiml (noix ,klva(i,j,isvv))
        zvlop (i,j) =   oiml (noilp,klva(i,j,isvv))
        zsvob       =   oiml (nois ,klva(i,j,isvv))
        omysc (i,j) =   (c1 + zdds(i,j,ivrs))                                  &
!   vertical correlation: pot. temperature differences are scaled by a larger
!   constant in the stratosphere than in the troposphere (larger mean stability)
                      * EXP(- (  MAX( MIN( zsvob,ztp ) - zspr(i,j,msi) , c0 )  &
                               + MAX( zsvob - MAX( zspr(i,j,msi),ztp ) , c0 )  &
                                 * zqts ) **2  *  rdsprsrb                     &
                            - zdds(i,j,ivrs)  )                                &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * oiml (noiqc,klva(i,j,isvv))
      ELSEIF (klvb(i,j,isvv) > -1) THEN
        zoics (i,j) =   oiml (noix ,klvb(i,j,isvv))
        zvlop (i,j) =   oiml (noilp,klvb(i,j,isvv))
        zsvob       =   oiml (nois ,klvb(i,j,isvv))
        omysc (i,j) =   (c1 + zdds(i,j,ivrs))                                  &
                      * EXP(- (  MAX( MIN( zspr(i,j,msi),ztp ) - zsvob , c0 )  &
                               + MAX( zspr(i,j,msi) - MAX( zsvob,ztp ) , c0 )  &
                                 * zqts ) **2  *  rdsprsrt                     &
                            - zdds(i,j,ivrs)  )                                &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * oiml (noiqc,klvb(i,j,isvv))
! no spreading along isentropes where vertical gradient of potential temperature
! is less than a threshold
      ELSE
        zoics (i,j) =  c0
        zvlop (i,j) =  c0
        omysc (i,j) =  c0
      ENDIF
    ENDDO
! conveying temperature increments to potential temperature increments
    IF (ivrs == 2) THEN
!CDIR NODEP,VOVERTAKE,VOB
      DO ixc = 1 , ncutof(irep)
        i = icutof(ixc,irep)
        j = jcutof(ixc,irep)
        zoics (i,j) = zoics(i,j) * zeklop(i,j) * EXP( -rdocp *zvlop(i,j) )
      ENDDO
    ENDIF
  ENDIF

! for cases a and b:
! If spreading along isentropes: spreading along horizontal surfaces must also
! be applied in areas of weak stability. However, the vertical correlations
! are still expressed as function of potential temperature differences, and
! non-isotropic weights (as function of height differences) are not possible.

  IF ((modespr == 2) .AND. (lsprx)) THEN
    dthreshr= c1 / MAX( thresh1 - thresh2 , epsy )
    sprmx   = c1
    IF (k == ktth) sprmx = c05
    msip = 3 - msi
    nois = noiz
    noic = noith
!   isvv = ivrs
    isvv = 1
!CDIR NODEP,VOVERTAKE,VOB
    DO ixc = 1 , ncutof(irep)
      i = icutof(ixc,irep)
      j = jcutof(ixc,irep)
! 'sprdv': fraction of obs. increment to be spread along isentropes;
! 1-sprdv: fraction of obs. increment to be spread horizontally
      sprdv (i,j) = MAX( c0 , MIN( sprmx, (zthvg(i,j)-thresh2) *dthreshr ) )
      IF (c1-sprdv(i,j) < epsy) THEN
        zoizs (i,j) =  c0
        zvlop (i,j) =  c0
        omyzc (i,j) =  c0
      ELSEIF ((klva(i,j,isvv) > -1) .AND. (klvb(i,j,isvv) > -1)) THEN
        zfvip       =   ( oiml (nois,klvb(i,j,isvv)) - zspr(i,j,msip) )        &
                      / ( oiml (nois,klvb(i,j,isvv))                           &
                         -oiml (nois,klva(i,j,isvv)) )
        zoizs (i,j) =     zfvip     * oiml (noix ,klva(i,j,isvv))              &
                       + (c1-zfvip) * oiml (noix ,klvb(i,j,isvv))
        zvlop (i,j) =     zfvip     * oiml (noilp,klva(i,j,isvv))              &
                       + (c1-zfvip) * oiml (noilp,klvb(i,j,isvv))
        zsvob       =   oiml (noic,klva(i,j,isvv))
        zwvip       =  (  MAX( MIN( zsvob,ztp ) - zspr(i,j,msi) , c0 )         &
                        + MAX( zsvob - MAX( zspr(i,j,msi),ztp ) , c0 )*zqts )  &
                      * oiml (noivcr,klva(i,j,isvv))                           &
                     +  oiml (noivlr,klva(i,j,isvv))
        zsvob       =   oiml (noic,klvb(i,j,isvv))
        zwvip2      =  (  MAX( MIN( zspr(i,j,msi),ztp ) - zsvob , c0 )         &
                        + MAX( zspr(i,j,msi) - MAX( zsvob,ztp ) , c0 )*zqts )  &
                      * oiml (noivcr,klvb(i,j,isvv))                           &
                     +  oiml (noivlr,klvb(i,j,isvv))
        zwvip       =  (zwvip * zwvip2 / MAX( zwvip+zwvip2,epsy ) *rdsprir) **2
        omyzc (i,j) =   (c1 + zdds(i,j,ivrs))                                  &
                      * EXP(- zdds(i,j,ivrs)  - zwvip )                        &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * (  zfvip     * oiml (noiqc,klva(i,j,isvv))             &
                         +(c1-zfvip) * oiml (noiqc,klvb(i,j,isvv)) )
      ELSEIF (klva(i,j,isvv) > -1) THEN
        zoizs (i,j) =   oiml (noix ,klva(i,j,isvv))
        zvlop (i,j) =   oiml (noilp,klva(i,j,isvv))
        zsvob       =   oiml (noic ,klva(i,j,isvv))
        omyzc (i,j) =   (c1 + zdds(i,j,ivrs))                                  &
                      * EXP(- (  MAX( MIN( zsvob,ztp ) - zspr(i,j,msi) , c0 )  &
                               + MAX( zsvob - MAX( zspr(i,j,msi),ztp ) , c0 )  &
                                 * zqts ) **2  *  rdsprsrb                     &
                            - zdds(i,j,ivrs)  )                                &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * oiml (noiqc,klva(i,j,isvv))
      ELSEIF (klvb(i,j,isvv) > -1) THEN
        zoizs (i,j) =   oiml (noix ,klvb(i,j,isvv))
        zvlop (i,j) =   oiml (noilp,klvb(i,j,isvv))
        zsvob       =   oiml (noic ,klvb(i,j,isvv))
        omyzc (i,j) =   (c1 + zdds(i,j,ivrs))                                  &
                      * EXP(- (  MAX( MIN( zspr(i,j,msi),ztp ) - zsvob , c0 )  &
                               + MAX( zspr(i,j,msi) - MAX( zsvob,ztp ) , c0 )  &
                                 * zqts ) **2  *  rdsprsrt                     &
                            - zdds(i,j,ivrs)  )                                &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * oiml (noiqc,klvb(i,j,isvv))
      ELSE
        zoizs (i,j) =  c0
        zvlop (i,j) =  c0
        omyzc (i,j) =  c0
      ENDIF
    ENDDO
!   isvv = 3*(msi-1) + ivrs
    isvv =    msi

! Merging of spreading along isentropes and of horizontal spreading

!CDIR NODEP,VOVERTAKE,VOB
    DO ixc = 1 , ncutof(irep)
      i = icutof(ixc,irep)
      j = jcutof(ixc,irep)
! conveying temperature increments to potential temperature increments
      IF (ivrs == 2)                                                           &
        zoizs (i,j)   = zoizs(i,j) * zeklop(i,j) * EXP( -rdocp *zvlop(i,j) )
      omywth        =     sprdv(i,j)  *omysc(i,j)
      omywp         = (c1-sprdv(i,j)) *omyzc(i,j)
      omysc (i,j)   =  omywth  +  omywp
      IF (omysc(i,j) > epsy)                                                   &
        zoics (i,j)   = (  omywth *zoics(i,j)                                  &
                         + omywp  *zoizs(i,j)) / omysc(i,j)
    ENDDO
  ENDIF


! case c: Horizontal spreading without non-isotropic corretion
!         (this case is treated separately for speed-up)
! ------------------------------------------------------------

! case c1: specific humdity (other than potential temperature)

  IF ((modespr <= 1) .AND. (.NOT. lnoniso) .AND. (lsprx) .AND. (ivrs /= 2)) THEN
!CDIR NODEP,VOVERTAKE,VOB
    DO ixc = 1 , ncutof(irep)
      i = icutof(ixc,irep)
      j = jcutof(ixc,irep)
      IF ((klva(i,j,isvv) > -1) .AND. (klvb(i,j,isvv) > -1)) THEN
        zfvip       =   ( oiml (nois  ,klvb(i,j,isvv)) - zspr(i,j,msi) )       &
                      / ( oiml (nois  ,klvb(i,j,isvv))                         &
                         -oiml (nois  ,klva(i,j,isvv)) )
        zoics (i,j) =     zfvip     * oiml (noix,klva(i,j,isvv))               &
                       + (c1-zfvip) * oiml (noix,klvb(i,j,isvv))
        zwvip       =    (oiml (nois  ,klva(i,j,isvv)) - zspr(i,j,msi))        &
                        * oiml (noivcr,klva(i,j,isvv))                         &
                       +  oiml (noivlr,klva(i,j,isvv))
        zwvip2      =  - (oiml (nois  ,klvb(i,j,isvv)) - zspr(i,j,msi))        &
                        * oiml (noivcr,klvb(i,j,isvv))                         &
                       +  oiml (noivlr,klvb(i,j,isvv))
        zwvip       =  (zwvip * zwvip2 / MAX( zwvip+zwvip2,epsy ) *rdsprir) **2
        omysc (i,j) =   (c1 + zdds(i,j,ivrs))                                  &
                      * EXP(- zdds(i,j,ivrs)  - zwvip )                        &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * (  zfvip     * oiml (noiqc,klva(i,j,isvv))             &
                         +(c1-zfvip) * oiml (noiqc,klvb(i,j,isvv)) )
      ELSEIF (klva(i,j,isvv) > -1) THEN
        zoics (i,j) =   oiml (noix ,klva(i,j,isvv))
        omysc (i,j) =   (c1 + zdds(i,j,ivrs))                                  &
                      * EXP(- ( oiml (nois,klva(i,j,isvv))                     &
                               -zspr(i,j,msi)) **2  *  rdsprsrb                &
                            - zdds(i,j,ivrs)  )                                &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * oiml (noiqc,klva(i,j,isvv))
      ELSEIF (klvb(i,j,isvv) > -1) THEN
        zoics (i,j) =   oiml (noix ,klvb(i,j,isvv))
        omysc (i,j) =   (c1 + zdds(i,j,ivrs))                                  &
                      * EXP(- ( oiml (nois,klvb(i,j,isvv))                     &
                               -zspr(i,j,msi)) **2  *  rdsprsrt                &
                            - zdds(i,j,ivrs)  )                                &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * oiml (noiqc,klvb(i,j,isvv))
      ELSE
        zoics (i,j) =  c0
        omysc (i,j) =  c0
      ENDIF
    ENDDO

! case c2: potential temperature

  ELSEIF ((modespr <= 1) .AND. (.NOT. lnoniso) .AND. (lsprx)) THEN
!CDIR NODEP,VOVERTAKE,VOB
    DO ixc = 1 , ncutof(irep)
      i = icutof(ixc,irep)
      j = jcutof(ixc,irep)
      IF ((klva(i,j,isvv) > -1) .AND. (klvb(i,j,isvv) > -1)) THEN
        zfvip       =   ( oiml (nois  ,klvb(i,j,isvv)) - zspr(i,j,msi) )       &
                      / ( oiml (nois  ,klvb(i,j,isvv))                         &
                         -oiml (nois  ,klva(i,j,isvv)) )
        zvlop (i,j) =      zfvip     * oiml (noilp,klva(i,j,isvv))             &
                        + (c1-zfvip) * oiml (noilp,klvb(i,j,isvv))
        zoics (i,j) =  (   zfvip     * oiml (noix ,klva(i,j,isvv))             &
                        + (c1-zfvip) * oiml (noix ,klvb(i,j,isvv)) )           &
                      * zeklop(i,j) * EXP( -rdocp *zvlop(i,j) )
        zwvip       =    (oiml (nois  ,klva(i,j,isvv)) - zspr(i,j,msi))        &
                        * oiml (noivcr,klva(i,j,isvv))                         &
                       +  oiml (noivlr,klva(i,j,isvv))
        zwvip2      =  - (oiml (nois  ,klvb(i,j,isvv)) - zspr(i,j,msi))        &
                        * oiml (noivcr,klvb(i,j,isvv))                         &
                       +  oiml (noivlr,klvb(i,j,isvv))
        zwvip       =  (zwvip * zwvip2 / MAX( zwvip+zwvip2,epsy ) *rdsprir) **2
        omysc (i,j) =   (c1 + zdds(i,j,ivrs))                                  &
                      * EXP(- zdds(i,j,ivrs)  - zwvip )                        &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * (  zfvip     * oiml (noiqc,klva(i,j,isvv))             &
                         +(c1-zfvip) * oiml (noiqc,klvb(i,j,isvv)) )
      ELSEIF (klva(i,j,isvv) > -1) THEN
        zoics (i,j) =   oiml (noix ,klva(i,j,isvv))                            &
                      * zeklop(i,j) * zemkpml(ista,nba1,ivrs)
        omysc (i,j) =   (c1 + zdds(i,j,ivrs))                                  &
                      * EXP(- ( oiml (nois,klva(i,j,isvv))                     &
                               -zspr(i,j,msi)) **2  *  rdsprsrb                &
                            - zdds(i,j,ivrs)  )                                &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * oiml (noiqc,klva(i,j,isvv))
      ELSEIF (klvb(i,j,isvv) > -1) THEN
        zoics (i,j) =   oiml (noix ,klvb(i,j,isvv))                            &
                      * zeklop(i,j) * zemkpml(ista,nba2,ivrs)
        omysc (i,j) =   (c1 + zdds(i,j,ivrs))                                  &
                      * EXP(- ( oiml (nois,klvb(i,j,isvv))                     &
                               -zspr(i,j,msi)) **2  *  rdsprsrt                &
                            - zdds(i,j,ivrs)  )                                &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * oiml (noiqc,klvb(i,j,isvv))
      ELSE
        zoics (i,j) =  c0
        omysc (i,j) =  c0
      ENDIF
    ENDDO
  ENDIF


! Correction for Cressman structure function instead of auto-
! regressive function (in this case, 'cnondiv' should be zero)
! ------------------------------------------------------------
  IF ((.NOT. lautorg) .AND. (lsprx)) THEN
!CDIR NODEP,VOVERTAKE,VOB
    DO ixc = 1 , ncutof(irep)
      i = icutof(ixc,irep)
      j = jcutof(ixc,irep)
      zdds2 = zdds(i,j,ivrs) **2
      omysc (i,j) = omysc(i,j) * MAX( c0, c1-zdds2 ) / (c1 + zdds2)            &
                               * EXP(zdds(i,j,ivrs)) / (c1 + zdds(i,j,ivrs))
    ENDDO
  ENDIF


!-------------------------------------------------------------------------------
!  Section 3: Storage of the spreaded increments and the weights and updating
!             the quantities needed for multiple observations
!-------------------------------------------------------------------------------

  IF (itim == 2) THEN
!   IF ((lsprx) .AND. (lsprx1) .AND. (ltiml(ista,ivrs))) THEN
    IF ((ltiml(ista,ivrs)) .AND. ((lsprx) .OR. (lsprx1))) THEN
! For temporal linear interpolation selected: if the nudging weight exclusive
! temporal weight for one obs. increment 'omykl' is larger than that 'omyks' for
! the other, then the temporal weighting for this first increment will be a mix-
! ture of temporal interpolation and of the sawtooth-shaped function valid for
! single observations. The fraction of temporal interpolation is 'omyks /omykl'.
! This allows to deal with temporal interpolation of a 'complete' and an
! incomplete sounding.
      DO j   = jstaspr , jendspr
        DO i = istaspr , iendspr
!         IF ((lcutof(i,j,1)) .AND. (lcutof(i,j,irep))) THEN
          IF ((lcutof(i,j,1)) .OR. (lcutof(i,j,irep))) THEN
            qmyk1    = MIN( omysc(i,j) / MAX( omys1(i,j),epsy ) , c1 )
            qmyk2    = MIN( omys1(i,j) / MAX( omysc(i,j),epsy ) , c1 )
            omys1 (i,j)  =  omys1(i,j) * (      qmyk1  *zwtml(ista,1,ivrs)     &
                                          + (c1-qmyk1) *zwtml(ista,3,ivrs))
            omysc (i,j)  =  omysc(i,j) * (      qmyk2  *zwtml(ista,2,ivrs)     &
                                          + (c1-qmyk2) *zwtml(ista,4,ivrs))
            omykm (i,j)  =  omys1(i,j) + omysc(i,j)
            zsoem (i,j)  = (omykm(i,j) + fmultw) * (  omys1(i,j) *zoic1(i,j)   &
                                                    + omysc(i,j) *zoics(i,j))
            omyc2 (i,j)  =  omykm(i,j) * omykm(i,j)
!         ELSEIF (lcutof(i,j,irep)) THEN
!           omykm (i,j)  =  omysc(i,j) * zwtml(ista,4,ivrs)
!           omyc2 (i,j)  =  omykm(i,j) * omykm(i,j)
!           zsoem (i,j)  = (omyc2(i,j) + fmultw *omykm(i,j)) * zoics(i,j)
!         ELSEIF (lcutof(i,j,1)) THEN
!           omykm (i,j)  =  omys1(i,j) * zwtml(ista,3,ivrs)
!           omyc2 (i,j)  =  omykm(i,j) * omykm(i,j)
!           zsoem (i,j)  = (omyc2(i,j) + fmultw *omykm(i,j)) * zoic1(i,j)
!         ENDIF
!         IF ((lcutof(i,j,1)) .OR. (lcutof(i,j,2))) THEN
            omy (i,j,ivrs,ityp) = omy(i,j,ivrs,ityp) + omykm(i,j)
            om2 (i,j,ivrs,ityp) = om2(i,j,ivrs,ityp) + omyc2(i,j)
            zwi (i,j,ivar,ityp) = zwi(i,j,ivar,ityp) + zsoem(i,j)
          ENDIF
        ENDDO
      ENDDO
    ELSE
      IF (lsprx) THEN
!CDIR NODEP,VOVERTAKE,VOB
        DO ixc = 1 , ncutof(irep)
          i = icutof(ixc,irep)
          j = jcutof(ixc,irep)
          omysc (i,j)    = omysc(i,j) * zwtml(ista,4,ivrs)
          omyc2 (i,j)    = omysc(i,j) * omysc(i,j)
          omy (i,j,ivrs,ityp) = omy(i,j,ivrs,ityp) +   omysc(i,j)
          om2 (i,j,ivrs,ityp) = om2(i,j,ivrs,ityp) +   omyc2(i,j)
          zwi (i,j,ivar,ityp) = zwi(i,j,ivar,ityp) + ( omysc(i,j) *fmultw      &
                                                      +omyc2(i,j) ) * zoics(i,j)
        ENDDO
      ENDIF
      IF (lsprx1) THEN
!CDIR NODEP,VOVERTAKE,VOB
        DO ixc = 1 , ncutof(1)
          i = icutof(ixc,1)
          j = jcutof(ixc,1)
          omys1 (i,j)    = omys1(i,j) * zwtml(ista,3,ivrs)
          omy12 (i,j)    = omys1(i,j) * omys1(i,j)
          omy (i,j,ivrs,ityp) = omy(i,j,ivrs,ityp) +   omys1(i,j)
          om2 (i,j,ivrs,ityp) = om2(i,j,ivrs,ityp) +   omy12(i,j)
          zwi (i,j,ivar,ityp) = zwi(i,j,ivar,ityp) + ( omys1(i,j) *fmultw      &
                                                      +omy12(i,j) ) * zoic1(i,j)
        ENDDO
      ENDIF
    ENDIF

  ENDIF

! 'zoic1', 'omys1' need not have an index for the variable type 'ivar', 'ivrs',
! as long as the inner loop in the calling procedure is over the time index
! 'itim' and the outer loop is over 'ivrs'

  IF ((itim == 1) .AND. (lsprx)) THEN
!CDIR NODEP,VOVERTAKE,VOB
    DO ixc = 1 , ncutof(1)
      i = icutof(ixc,1)
      j = jcutof(ixc,1)
      zoic1 (i,j) = zoics(i,j)
      omys1 (i,j) = omysc(i,j)
    ENDDO
  ELSEIF ((itim == 2) .AND. (lsprx1)) THEN
!CDIR NODEP,VOVERTAKE,VOB
    DO ixc = 1 , ncutof(1)
      i = icutof(ixc,1)
      j = jcutof(ixc,1)
      zoic1 (i,j) = c0
      omys1 (i,j) = c0
    ENDDO
  ENDIF

! printout for control   (if statements must equal those at the
! --------------------    initialization of fields 'omysc', 'zoics')

  IF ((ntstep <= 1) .AND. (lwonl) .AND. (lsprx)) THEN
!   io = imladm(ista,3)   ! io,jo known globally
!   jo = imladm(ista,4)
!   IF (((io == ionl) .AND. (jo == jonl) .AND. (k >= ke-6)) .OR. (lfrsloc)) THEN
    IF  ((io == ionl) .AND. (jo == jonl) .AND. (k >= ke-6)) THEN
      i    = MAX( INT( MIN( io , iendspr ) ,iintegers) , istaspr )
      j    = MAX( INT( MIN( jo , jendspr ) ,iintegers) , jstaspr )
      jpra = MAX( INT( jo  -5 ,iintegers) , jstaspr )
      jpre = MIN( INT( jpra+9 ,iintegers) , jendspr )
      jpra = MAX( INT( jpre-9 ,iintegers) , jstaspr )
      IF (klva(i,j,isvv) > -1) THEN
        IF (modespr == 2) nois = noith
        zsvob   = oiml (nois,klva(i,j,isvv))
        omycorl = EXP( -(zsvob - zspr(i,j,msi))**2 * rdsprsrb )
        WRITE( nupr,'(''spr_tqz: k='',I2,'', omysc, zoic''                     &
                    &,2F8.2, 2F7.1, 2F8.5)' )    k, rdsprsb, rdsprst           &
             , zsvob, zspr(i,j,msi), omycorl, omysc(i,j)
      ENDIF
      j     = jo   - jpra + 1
      nprae = jpre - jpra + 1
      WRITE( yformat, '(''('',I2,''(6X,L1),2I4)'')' ) nprae
      WRITE( nupr, yformat )  (lcutof(i,joc,irep), joc=jpra,jpre), ivar, j
      WRITE( yformat, '(''('',I2,''F7.4,I4    )'')' ) nprae
      WRITE( nupr, yformat )  (omysc (i,joc     ), joc=jpra,jpre), ivar
!     IF (ivrs == 3) WRITE( yformat, '(''('',I2,''F7.5,I4   )'')' ) nprae
      WRITE( nupr, yformat )  (zoics (i,joc     ), joc=jpra,jpre), ivar
      IF (modespr == 2) WRITE( nupr,yformat ) (sprdv(i,joc),joc=jpra,jpre),ivar
      WRITE( yformat, '(''('',I2,''F7.2,I4    )'')' ) nprae
      WRITE( nupr, yformat )  (zthvg (i,joc     ), joc=jpra,jpre), ivar
      WRITE( nupr, yformat )  (zdds  (i,joc,ivrs), joc=jpra,jpre), ivar
      WRITE( yformat, '(''('',I2,''F7.1,I4    )'')' ) nprae
      WRITE( nupr, yformat )  (zspr  (i,joc,msi ), joc=jpra,jpre), ivar
      lfrsloc = .FALSE.
    ENDIF
    IF (      (((io == ionl) .AND. (jo == jonl)) .OR. (MOD(ista ,4) == 0))     &
        .AND. (k == ke-5) .AND. (lsprx1))                                      &
      WRITE( nupr,'(''lti_mult '',A ,'', itim,ivar'',2I2,'', zoics/omys1 ''    &
                  &,F8.5,F7.4)' )                                              &
             ystid, itim, ivar, zoic1(ionl,jonl), omys1(ionl,jonl)
  ENDIF


!-------------------------------------------------------------------------------
! End of module procedure mult_spread_mass
!-------------------------------------------------------------------------------

END SUBROUTINE mult_spread_mass


!-------------------------------------------------------------------------------
!+ Module procedure for spreading multi-level obs incr. along model levels
!-------------------------------------------------------------------------------

SUBROUTINE mult_spread_modlev ( ivar , ityp , k )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure performs the weighting and spreading (extrapolation)
!   along model levels of observation increments from (at most 2) multi-level
!   reports from one observing station.
!
! Method:
!   Spreading of multi-level data increments to the target grid pts. with the
!   largest weighting along model levels.
!   Computation of the sum of 'relative nudging weights' by applying the
!   temporal and spatial correlations and quality factors.
!   Option for a non-isotropic correction to horizontal correlation (see below).
!   Humidity: Computation of specific humidity ('qd') increments from the
!             spreaded relative humidity ('RH') increments, so that the nudging
!             pulls towards observed 'RH' instead of observed 'qd'.
!   Vertical correlation as Gaussian of height differences ('dz') (modespr <= 1)
!   or potential temperature differences (modespr == 2) between the obs. point
!   and the grid pt. at level 'k' at the obs. location
!   Non-isotropy as Gaussian of potential temperature resp. height differences
!   on model levels
!     correlation functions: - in z     : EXP( -(g(z1-z2)/r*Tv1)**2 / dzscal) )
!                            - in theta : EXP( -(theta1-theta2)**2 / dthscal) )
!
! Written by        :  Christoph Schraff, DWD  (original version: 12.09.97)
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
    ivar             ,& ! index of spreaded quantity: 1=(u,v); 3=T; 4=RH
    ityp             ,& ! index of weighted increment array
    k                   ! index of current vertical model level


! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    modespr          ,& ! mode of spreading, see above: 'Method'
    mni              ,& ! index of parameter defining non-isotropic weights
    ivay             ,& ! index of 2nd comp. of spreaded vector quantity: 2=v
    ivrs             ,& ! index of spreaded vector quantity: 1=(u,v); 2=T; 3=RH
    noix   , noiy    ,& ! indices of spreaded quantity in obs. increment field
    nois             ,& ! index of (selected) spreading param. in obs incr field
    noiqc            ,& ! index of quality factor in obs. increment field
    noivlr           ,& ! index of vertical weight given to the obs. incr. due
                        ! to the distance to the obs. at the vertical interpol.
!   nmloi            ,& ! index of report in the observation increment field
!   nmloi1 , nmloi2  ,& ! dito
    kdat   , kdact   ,& ! level index of the required obs incr. in the obs incr.
    kdat1  , kdat2   ,& ! dito                                      \  profile
!   kd               ,& ! index of model level with the required obs. increment
    kkbot  , kktop   ,& ! index range used for the present report
    klv              ,& ! loop index in vertical direction
    i      , j       ,& ! loop indices in horizontal direction
    itima  , itime   ,& ! range of time loop over reports with non-zero weights
!   io     , jo      ,& ! local indices of location of observation
!   istaspr, jstaspr ,& ! lower left corner of domain containing area of infl.
!   iendspr, jendspr ,& ! upper right corner of domain containing area of infl
    ixc              ,& ! loop index over influenced grid points
    nprae            ,& ! number of printed grid pts.
    jpra   , jpre    ,& ! lowest , highest meridional index of printed grid pts.
    joc                 ! loop index for printout

  REAL    (KIND=wp   )     ::  &
    fmultw           ,& ! 0. or 1., depending on multiple observation weighting
    zqgnudg          ,& ! nudging coeff. relative to the coeff. for TEMP / PILOT
    ztp              ,& ! value of spreading param. (z, TH) at (std atm) tropop.
    zsvob            ,& ! spreading param. (z, TH) at the obs. incr. point
    zvdist           ,& ! vertical distance [spreading paramater units] between
                        ! obs. incr. point and the model level at the obs. loc.
    zvdspr           ,& ! vertical difference of spreading parameter between the
                        ! obs incr. & the current model level at the obs locat.
    rdsprs           ,& ! vertical correlation scale as by namelist
    rdsprsb, rdsprst ,& ! lower / upper vertical correlation scale
    rdspri           ,& ! vert. correl. scale for weighting the obs. incr. level
    rdsprni          ,& ! square of Gauss. non-isotropy radius,'zspr(mni)'-units
    rdsprnir         ,& ! reciproke of 'rdsprni'
    uincr1 , uincr2  ,& ! weighted obs. incr. of zonal wind (for 2 obs at 1 sta)
    vincr1 , vincr2  ,& ! wt. obs. incr. of meridional wind (for 2 obs at 1 sta)
    xincr1           ,& ! /  obs. incr. of scalar quanitiy (for 2 obs at 1 sta):
    xincr2           ,& ! \  pot. temperature or generalized relative humidity
!   zewcor1, zewcor2 ,& ! targeted      vapour pressure   (for 2 obs at 1 sta.)
!   zqdcor1, zqdcor2 ,& ! targeted      specific humidity (for 2 obs at 1 sta.)
    uincrt , vincrt  ,& ! sum of the 2 partially weighted wind comp. obs. incr.
    uinct2 , vinct2  ,& ! sum of the 2 partially squared weighted wind obs. incr
    xincrt , xinct2  ,& ! sum of the 2 partially weighted (**2) scalar obs. incr
    zdds2            ,& ! squared scaled horizontal distance to the obs. station
    zuoem  , zvoem   ,& ! sum of the 2 spreaded obs. incr. of wind vector
    zuoem2 , zvoem2  ,& ! sum of the 2 squared spreaded obs. incr. of wind
    omyk             ,& ! local part of the nudging weight
    omyk1  , omyk2   ,& ! local part of the nudging weights
    qmyk1            ,& ! /  local part of the nudging weights
    qmyk2            ,& ! \  taking into account partial temporal interpolation
    omyvc            ,& ! vertical correlation, i.e. vertical nudging weight
    omykt            ,& ! omyk1 + omyk2 : sum of the 2 local weights
    omykm            ,& ! omykt *omysc : sum of the 2 correlations (weight)
    omykh            ,& ! omykt *omysc *omysc
    omysc2           ,& ! omysc **2
    omyk2t           ,& ! omyk1 **2 + omyk2 **2
    omytt               ! weight used for spreading of potential temperature

  LOGICAL                  ::  &
    lnoniso          ,& ! non-isotrophic horizontal correlation
    lwabl               ! nudging weight reduced in the atmospheric bound. layer

  CHARACTER (LEN=15) yformat


! LOGICAL                  , SAVE ::  &
!   lsprx1              ! non-zero temporal and vertical weights from 1st report


! Local (automatic) arrays:
! -------------------------
  REAL    (KIND=wp   )     ::  &
    omysc (istart:iend,jstart:jend)   ! spatial nudging weight for a single obs.

! Initialized variables:
! ----------------------

  LOGICAL                  , SAVE ::  &
    lfrsloc = .TRUE.    ! .TRUE if local output is to be printed the first time

!------------ End of header ----------------------------------------------------


 
!-------------------------------------------------------------------------------
! Begin Subroutine mult_spread_modlev
!-------------------------------------------------------------------------------
  
!-------------------------------------------------------------------------------
!  Section 1: Preliminaries
!-------------------------------------------------------------------------------

  ivay    = ivar + 1
  ivrs    = MAX( ivar - 1 , 1 )
  modespr = msprpar

  zqgnudg = c1
  IF (kobtyp == nairep) THEN
    IF (kcdtyp /= nmodes) zqgnudg = gnudgar(ivar) / MAX( gnudg(ivar) ,epsy )
    IF (kcdtyp == nmodes) zqgnudg = gnudgms(ivar) / MAX( gnudg(ivar) ,epsy )
  ENDIF
  IF (kobtyp == nsattv) zqgnudg = gnudgtv(ivar) / MAX( gnudg(ivar) ,epsy )
  IF (kobtyp == ngps  ) zqgnudg = gnudggp       / MAX( gnudg(ivar) ,epsy )
  IF ((kobtyp == ngps ) .AND. (ivar /= 4))  zqgnudg = c0
  ztp     = zstpml (ista)

  fmultw = REAL ( kwtyp(ityp) - 1, wp )

! get indices for observation increment and other arrays

  mni  = 3 - MAX( modespr , i1 )
  IF (ivrs == 1) THEN
    noix   = noiu
    noiy   = noiv
    noiqc  = noiuqc
    noivlr = noiulr
  ELSEIF (ivrs == 2) THEN
    noix   = noit
    noiqc  = noitqc
    noivlr = noitlr
  ELSEIF (ivrs == 3) THEN
    noix   = noirh
    noiqc  = noiqqc
    noivlr = noiqlr
  ENDIF
  nois = noiz
  IF (modespr == 2) nois = noith

  omyk1 = c0
  omyk2 = c0

! checking if reports may influence current time and model level

  itima = 2
  itime = 1
! IF (      (MAX( zwtml(ista,1,ivrs) , zwtml(ista,3,ivrs) ) > epsy)            &
  IF (      (k <= kviflml(ista,1,ivrs))                                        &
      .AND. (k >= kviflml(ista,2,ivrs)))   itima = 1
! IF (      (MAX( zwtml(ista,2,ivrs) , zwtml(ista,4,ivrs) ) > epsy)            &
  IF (      (k <= kviflml(ista,3,ivrs))                                        &
      .AND. (k >= kviflml(ista,4,ivrs)))   itime = 2

  DO itim = itima , itime

!   nmloi = knoiml(ista,itim)


!-------------------------------------------------------------------------------
!  Section 2: Determination of the index of the appropriate level
!             in the obs. increment vertical profile
!-------------------------------------------------------------------------------

!   IF (      (k <= ke-mbotlv(ista,itim,ivrs))                                 &
!       .AND. (k >=  1+mtoplv(ista,itim,ivrs))) THEN
!     kdat = ke + 1 - k
!   ELSEIF (k > ke-mbotlv(ista,itim,ivrs)) THEN
!     kdat =  1 + mbotlv(ista,itim,ivrs)
!   ELSE
!     kdat = ke - mtoplv(ista,itim,ivrs)
!   ENDIF
! get 'kdat' if (msprpar > 0)
!   IF (levdiff) THEN
!     kd = ke + 1 - MAX( kdat,1 )
      kdat   = 0
      zvdist = zallhig
      kkbot  = kkoiml(ista,itim)
      kktop  = kkoiml(ista,itim) + mszlev(ista,itim) - 1
      DO klv = kkbot , kktop
!       IF (      (ABS( oiml(nois,klv) -zsprml(ista,kd) ) < epsy)              &
!           .AND. (oiml(noix,klv) > rmdich)) kdat = klv
        IF (      (ABS( oiml(nois,klv) -zsprml(ista,k ) ) < zvdist)            &
            .AND. (oiml(noix,klv) > rmdich)) THEN
          kdat   = klv
          zvdist = ABS( oiml(nois,klv) - zsprml(ista,k) )
        ENDIF
      ENDDO
      IF ((kdat == 0) .AND. (ntstep <= 1) .AND. (lwonl))                       &
        WRITE( nupr, '(''mult_spr_ml: NO DATA of var.'',I2,'' found at ''      &
                     &,''station '',A,'' for model level'',I3)' ) ivar, ystid, k
!   ENDIF

!-------------------------------------------------------------------------------
!  Section 3: (Horizontally) local contributions to the nudging weights:
!             temporal weights, vertical weight (correlation),
!             quality factors, relative nudging coefficient
!-------------------------------------------------------------------------------

    IF (kdat >= 1) THEN

! Vertical correlation: /=1 only if no obs. incr. exists on current model level
! -----------------------------------------------------------------------------

      zvdspr = oiml(nois,kdat) - zsprml(ista,k)
      omyvc  = c1
      IF (ABS( zvdspr ) > epsy) THEN
        rdsprs  = SQRT( vcorls(ivar) )
        IF (modespr <= 1) THEN
          rdsprsb = rdsprs * fcorlml(ista,2*itim-1,ivrs)                       &
                           * zrtdgml(ista,2*itim-1,ivrs)
          rdsprst = rdsprs * fcorlml(ista,2*itim  ,ivrs)                       &
                           * zrtdgml(ista,2*itim  ,ivrs)
          omyvc   = EXP( -( MIN( zvdspr /rdsprst , c0 )                        &
                           +MAX( zvdspr /rdsprsb , c0 ))**2 )
        ELSE
          rdsprsb = rdsprs * fcorlml(ista,2*itim-1,ivrs)
          rdsprst = rdsprs * fcorlml(ista,2*itim  ,ivrs)
          zsvob   = oiml(nois,kdat)
          IF (zvdspr < 0) THEN
            omyvc = EXP(- (  MAX( MIN( zsprml(ista,k), ztp ) - zsvob , c0 )    &
                           + MAX( zsprml(ista,k) - MAX( zsvob, ztp ) , c0 )    &
                             / qst )  /  rdsprst ) **2
          ELSE
            omyvc = EXP(- (  MAX( MIN( zsvob, ztp ) - zsprml(ista,k) , c0 )    &
                           + MAX( zsvob - MAX( zsprml(ista,k), ztp ) , c0 )    &
                             / qst )  /  rdsprsb ) **2
          ENDIF
        ENDIF
      ENDIF

! Vertical correlation within the profile
! ---------------------------------------

      IF (modespr <= 1) rdspri = SQRT( MAX( vcorls(ivar) , 0.04_wp ) )
      IF (modespr == 2) rdspri = SQRT( MAX( vcorls(ivar) , 33.0_wp ) )
      omyvc  = omyvc * EXP( - (oiml(noivlr,kdat) /rdspri) **2 )


! Compute and store the total local contribution to the weights
! -------------------------------------------------------------

      omyk  =  oiml(noiqc,kdat) * zqgnudg * omyvc
    ELSE
      omyk  =  c0
    ENDIF

    IF (itim == 1) THEN
      omyk1  = omyk
      kdat1  = kdat
!     nmloi1 = nmloi
    ELSE
      omyk2  = omyk
      kdat2  = kdat
!     nmloi2 = nmloi
    ENDIF

  ENDDO

  IF (.NOT. ltiml(ista,ivrs)) THEN
    omyk1 = omyk1 * zwtml(ista,3,ivrs)
    omyk2 = omyk2 * zwtml(ista,4,ivrs)
  ELSE

! For temporal linear interpolation selected: if the nudging weight exclusive
! temporal weight for one obs. increment 'omykl' is larger than that 'omyks' for
! the other, then the temporal weighting for this first increment will be a mix-
! ture of temporal interpolation and of the sawtooth-shaped function valid for
! single observations. The fraction of temporal interpolation is 'omyks /omykl'.
! This allows to deal with temporal interpolation of a 'complete' and an
! incomplete sounding.
    qmyk1 = MIN( omyk2 / MAX( omyk1,epsy ) , c1 )
    qmyk2 = MIN( omyk1 / MAX( omyk2,epsy ) , c1 )
    omyk1 = omyk1 * (qmyk1 *zwtml(ista,1,ivrs) + (c1-qmyk1) *zwtml(ista,3,ivrs))
    omyk2 = omyk2 * (qmyk2 *zwtml(ista,2,ivrs) + (c1-qmyk2) *zwtml(ista,4,ivrs))
  ENDIF

! If local weights are zero, further computations (spreading) will be skipped
! Note: 'kdat?' needs to be known subsequently only if (omyk? > epsy)

  IF (omyk1 <= epsy) omyk1 = c0
  IF (omyk2 <= epsy) omyk2 = c0
  IF ((omyk1 <= epsy) .AND. (omyk2 <= epsy)) ivrs = 0

  IF ((ntstep <= 1) .AND. (lwonl) .AND. (k == ke-1) .AND. (ivrs > 0))          &
    WRITE( nupr,'(''mulmy: '',A ,2I4,I2,'',omyk1/2'',2F6.3,'' v-corr'',F8.3,I4 &
                &,'' dz'',F7.0,'' omyvc'',F9.6)' ) ystid, io_tot, jo_tot, ivar &
         , omyk1, omyk2, rdsprst, kdat, zvdspr, omyvc

!-------------------------------------------------------------------------------
!  Section 4: Horizontal correlations (weights), and influence of boundary layer
!-------------------------------------------------------------------------------

  IF (ivrs > 0) THEN

! set parameters used for the vertical weighting and non-isotropic correction
 
    rdsprni = vcsni  (ivar)
    IF (modespr <= 1) lnoniso = (rdsprni <= vcutnit) .AND. (lnisua)
    IF (modespr == 2) lnoniso = (rdsprni <= vcutnip) .AND. (lnisua)
    IF (modespr == 2) rdsprni = rdsprni * zrtgkml(ista,k) **2
    IF ((modespr <= 1) .AND. (zsprml(ista,k) > ztp)) rdsprni = qs2 * rdsprni
    rdsprnir = c1 / MAX( rdsprni , epsy )
    lwabl = (wablua(ivar)-c1 > epsy)
!   lti2 = (ltiml(ista,ivrs)) .AND. (omyk1 >= epsy) .AND. (omyk2 >= epsy)

! 1-dim. autoregressive horizontal correlation

    IF ((ivrs == 1) .AND. (lautorg)) THEN
!CDIR NODEP,VOVERTAKE,VOB
      DO ixc = 1 , ncutof(1)
        i = icutof(ixc,1)
        j = jcutof(ixc,1)
        omysc (i,j) =                         EXP(- zdds(i,j,ivrs) )
      ENDDO
    ELSEIF (lautorg) THEN
!CDIR NODEP,VOVERTAKE,VOB
      DO ixc = 1 , ncutof(1)
        i = icutof(ixc,1)
        j = jcutof(ixc,1)
        omysc (i,j) = (c1 + zdds(i,j,ivrs)) * EXP(- zdds(i,j,ivrs) )
      ENDDO
    ELSE

! 1-dim. Cressman-type horizontal correlation

!CDIR NODEP,VOVERTAKE,VOB
      DO ixc = 1 , ncutof(1)
        i = icutof(ixc,1)
        j = jcutof(ixc,1)
        zdds2 = zdds(i,j,ivrs) **2
        omysc (i,j) = MAX( c0 , c1-zdds2 ) / (c1 + zdds2)
      ENDDO
    ENDIF

! Non-isotropic correction to horizontal correlation,
! and reduction of weights within atmospheric boundary layer

    IF ((lnoniso) .OR. (lwabl)) THEN
!CDIR NODEP,VOVERTAKE,VOB
      DO ixc = 1 , ncutof(1)
        i = icutof(ixc,1)
        j = jcutof(ixc,1)
        IF (lnoniso) omysc (i,j) = omysc(i,j) * EXP(- ( znisml(ista,k)         &
                                                       -zspr(i,j,mni))**2      &
                                                      * rdsprnir )
        IF (lwabl)   omysc (i,j) = omysc(i,j) * zablfc(i,j,ivrs,1)
      ENDDO
    ENDIF

  ENDIF


!-------------------------------------------------------------------------------
!  Section 5: Spreading of observation increments of the horizontal wind vector
!             (updating the quantities needed for multiple observations)
!-------------------------------------------------------------------------------

  IF (ivrs == 1) THEN
    ivay = ivar + 1
    IF ((omyk1 > epsy) .AND. (omyk2 > epsy)) THEN
! order of calculations must be identical to the 'else'block for reproducibility
      uincr1 = omyk1 * oiml(noix,kdat1)
      vincr1 = omyk1 * oiml(noiy,kdat1)
      uincr2 = omyk2 * oiml(noix,kdat2)
      vincr2 = omyk2 * oiml(noiy,kdat2)
      omykt  = omyk1 + omyk2
      uincrt = uincr1 + uincr2
      vincrt = vincr1 + vincr2
      IF (.NOT. ltiml(ista,ivrs)) THEN
        uincr1 = omyk1 * uincr1
        vincr1 = omyk1 * vincr1
        uincr2 = omyk2 * uincr2
        vincr2 = omyk2 * vincr2
        omyk2t = omyk1 *omyk1 + omyk2 *omyk2
        uinct2 = uincr1 + uincr2
        vinct2 = vincr1 + vincr2
      ENDIF
    ELSE
      IF (omyk1 > epsy) THEN
        omykt  = omyk1
        kdact  = kdat1
      ELSE
        omykt  = omyk2
        kdact  = kdat2
      ENDIF
      uincrt = omykt * oiml(noix,kdact)
      vincrt = omykt * oiml(noiy,kdact)
      IF (.NOT. ltiml(ista,ivrs)) THEN
        omyk2t = omykt * omykt
        uinct2 = omykt * uincrt
        vinct2 = omykt * vincrt
      ENDIF
    ENDIF

    IF (ltiml(ista,ivrs)) THEN
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
                                                  + omysc(i,j) *zuoem * fmultw
        zwi (i,j,ivay,ityp) = zwi (i,j,ivay,ityp) + omysc2     *zvoem2         &
                                                  + omysc(i,j) *zvoem * fmultw
      ENDDO
    ENDIF


!-------------------------------------------------------------------------------
!  Section 6: Spreading of observation increments of potential temperature
!                                      or of generalized relative humidity
!-------------------------------------------------------------------------------

  ELSEIF ((ivrs == 2) .OR. (ivrs == 3)) THEN
    IF ((omyk1 > epsy) .AND. (omyk2 > epsy)) THEN
! order of calculations must be identical to the 'else'block for reproducibility
      xincr1  =  oiml(noix,kdat1)
      xincr2  =  oiml(noix,kdat2)
      IF (ivrs == 2)  xincr1  =  xincr1 * EXP( -rdocp *oiml(noilp,kdat1) )
      IF (ivrs == 2)  xincr2  =  xincr2 * EXP( -rdocp *oiml(noilp,kdat2) )
      xincr1 = omyk1 * xincr1
      xincr2 = omyk2 * xincr2
      omykt   = omyk1 + omyk2
      xincrt  = xincr1 + xincr2
      IF (.NOT. ltiml(ista,ivrs)) THEN
        xincr1 = omyk1 * xincr1
        xincr2 = omyk2 * xincr2
        omyk2t = omyk1 *omyk1 + omyk2 *omyk2
        xinct2 = xincr1 + xincr2
      ENDIF
    ELSE
      IF (omyk1 > epsy) THEN
        omykt  = omyk1
        kdact  = kdat1
      ELSE
        omykt  = omyk2
        kdact  = kdat2
      ENDIF
      xincrt = oiml(noix,kdact)
      IF (ivrs == 2)  xincrt  =  xincrt * EXP( -rdocp *oiml(noilp,kdact) )
      xincrt = omykt * xincrt
      IF (.NOT. ltiml(ista,ivrs)) THEN
        omyk2t = omykt * omykt
        xinct2 = omykt * xincrt
      ENDIF
    ENDIF
  ENDIF

! potential temperature
! ---------------------

  IF (ivrs == 2) THEN
    IF (ltiml(ista,ivrs)) THEN
!CDIR NODEP,VOVERTAKE,VOB
      DO ixc = 1 , ncutof(1)
        i = icutof(ixc,1)
        j = jcutof(ixc,1)
        omykm  =  omysc(i,j) * omykt
        omytt  =  omysc(i,j) *(omykm + fmultw) * zeklop(i,j)
        omy (i,j,ivrs,ityp) = omy (i,j,ivrs,ityp) + omykm
        om2 (i,j,ivrs,ityp) = om2 (i,j,ivrs,ityp) + omykm *omykm
        zwi (i,j,ivar,ityp) = zwi (i,j,ivar,ityp) + omytt *xincrt
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
        zwi (i,j,ivar,ityp) = zwi (i,j,ivar,ityp) + omytt  *xincrt *fmultw     &
                                                  + omysc2 *xinct2 *zeklop(i,j)
      ENDDO
    ENDIF

! generalized relative humidity
! -----------------------------

  ELSEIF (ivrs == 3) THEN
    IF (ltiml(ista,ivrs)) THEN
!CDIR NODEP,VOVERTAKE,VOB
      DO ixc = 1 , ncutof(1)
        i = icutof(ixc,1)
        j = jcutof(ixc,1)
        omykm  =  omysc(i,j) * omykt
        omykh  =  omysc(i,j) *(omykm + fmultw)
        omy (i,j,ivrs,ityp) = omy (i,j,ivrs,ityp) + omykm
        om2 (i,j,ivrs,ityp) = om2 (i,j,ivrs,ityp) + omykm *omykm
        zwi (i,j,ivar,ityp) = zwi (i,j,ivar,ityp) + omykh *xincrt
      ENDDO
    ELSE
!CDIR NODEP,VOVERTAKE,VOB
      DO ixc = 1 , ncutof(1)
        i = icutof(ixc,1)
        j = jcutof(ixc,1)
        omysc2 =  omysc(i,j) * omysc(i,j)
        omy (i,j,ivrs,ityp) = omy (i,j,ivrs,ityp) + omysc(i,j) *omykt
        om2 (i,j,ivrs,ityp) = om2 (i,j,ivrs,ityp) + omysc2 *omyk2t
        zwi (i,j,ivar,ityp) = zwi (i,j,ivar,ityp) + omysc2     *xinct2         &
                                                  + omysc(i,j) *xincrt *fmultw
      ENDDO
    ENDIF
  ENDIF


! printout for control
! --------------------

  IF ((ntstep <= 1) .AND. (lwonl) .AND. (ivrs > 0)) THEN
!   io = imladm(ista,3)   ! io,jo knowm globally
!   jo = imladm(ista,4)
!   IF ((      (io == ionl) .AND. (jo == jonl)                                 &
!        .AND. ((k >= ke-4) .OR. (k == ktp))) .OR. (lfrsloc)) THEN
    IF ((io == ionl) .AND. (jo == jonl) .AND. ((k > ke-6) .OR. (k == ktp))) THEN
      i    = MAX( INT( MIN( io , iendspr ) ,iintegers) , istaspr )
      j    = MAX( INT( MIN( jo , jendspr ) ,iintegers) , jstaspr )
      jpra = MAX( INT( jo  -5 ,iintegers) , jstaspr )
      jpre = MIN( INT( jpra+9 ,iintegers) , jendspr )
      jpra = MAX( INT( jpre-9 ,iintegers) , jstaspr )
      IF (lcutof(i,j,1))                                                       &
        WRITE( nupr,'(''spr_uvml: k='',I2,1X,A , '', omysc, zoic'',F8.5)' )    &
               k, ystid, omysc(i,j)
      j     = jo   - jpra + 1
      nprae = jpre - jpra + 1
      WRITE( yformat, '(''('',I2,''(6X,L1),2I4)'')' ) nprae
      WRITE( nupr, yformat )  (lcutof(i,joc,1), joc=jpra,jpre), ivar, j
      WRITE( yformat, '(''('',I2,''F7.4,I4    )'')' ) nprae
      WRITE( nupr, yformat )  (omysc (i,joc          ), joc=jpra,jpre), ivar
      WRITE( nupr, yformat )  (omy   (i,joc,ivrs,ityp), joc=jpra,jpre), ivar
      WRITE( nupr, yformat )  (zwi   (i,joc,ivar,ityp), joc=jpra,jpre), ivar
      lfrsloc = .FALSE.
    ENDIF
  ENDIF


!-------------------------------------------------------------------------------
! End of module procedure mult_spread_modlev
!-------------------------------------------------------------------------------

END SUBROUTINE mult_spread_modlev


!-------------------------------------------------------------------------------
!+ Module procedure for spreading multi-level obs increments of horiz. wind
!-------------------------------------------------------------------------------

SUBROUTINE mult_spread_wind ( lsprx , ityp , k )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure performs the weighting and spreading (extrapolation)
!   of wind vector observation increments from (at most 2) multi-level reports
!   from one observing station along surfaces (horizontal or isentropic) which
!   differ from the model surfaces.
!
! Method:
!   Spreading of multi-level data increments to the target grid pts. with the
!   largest weighting along specified surfaces according to 'modespr' (msprpar).
!   Computation of the sum of 'relative nudging weights' by applying the
!   temporal and spatial correlations and quality factors.
!   Option for non-isotropic lateral correlations (see below).
!   Mode of spreading:   'modespr' =
!     1: 'spreading along horizontal levels' :
!        - vertical correlation as Gaussian of height differences ('dz') between
!          the obs. point and the target grid pt. at level 'k'
!        - non-isotropy as Gaussian of pot. temp. differences on horiz. levels
!     2: 'spreading along isentropic surfaces' :
!        - vertical correlation as Gaussian of potential temperature differences
!          between the obs. point and the target grid pt. at level 'k'
!        - non-isotropy as Gaussian of height differences on isentropic surfaces
!        Optionally, spreading may be along model levels for target grid points
!        above a specified model level (see proc. *mult_spread_modlev*). At the
!        level below that level, or in the presence of weak thermal stability
!        at the target grid pt., spreading along isentropic surfaces is blended
!        with horizontal spreading for which in this case, however, vertical
!        correlations are still a Gaussian of potential temperature differences,
!        and non-isotropic weights are not possible.
!     correlation functions: - in z     : EXP( -(g(z1-z2)/r*Tv1)**2 / dzscal) )
!                            - in theta : EXP( -(theta1-theta2)**2 / dthscal) )
!     For target grid points within the vertical extent of the observation
!     increment profile, the weights are reduced in the following way:
!     1. The vertical distance to the nearest upper and lower observation is
!        assigned to be the sum of the 'distance' (as in '2.', see below) of the
!        observation increment level to its closest observation, and the (verti-
!        cal) distance of the obs. increment level to the target grid point.
!     2. An effective vertical distance is computed from the 2 (i.e. upper and
!        lower) distances in the same way as the resistance resulting from 2
!        parallel resistances:   d_eff = d_upper * d_lower / (d_upper + d_lower)
!     3. Weight function (like correlations):   EXP( -(d_eff)**2 / d?scal )
!
! Written by        :  Christoph Schraff, DWD  (original version: 28.08.97)
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
    lsprx               ! a report (i.e. obs. increments) exists with indices =
                        ! (ista,itim) and with non-zero temporal and vertical
                        ! weights for at least one datum in current report
                        ! ==> there are obs. increments to be spreaded

  INTEGER (KIND=iintegers), INTENT (IN)         ::       &
    ityp             ,& ! index of weighted increment array
    k                   ! index of current vertical model level


! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    modespr          ,& ! mode of spreading, see above: 'Method'
    msi    , msip    ,& ! index of (effective) spreading parameter
    mni              ,& ! index of parameter defining non-isotropic weights
    meqi             ,& ! index of parameter defining equidistant points
    ivar   , ivay    ,& ! indices of spreaded vector quantity: 1=u, 2=v
    ivrs             ,& ! index of spreaded vector quantity: 1=(u,v)
    irep             ,& ! index of report per station for '?cutof'
    isvv             ,& ! index of spreaded quantity depend. on spreading param.
!   nmloi            ,& ! index of report in the observation increment field
    noix   , noiy    ,& ! indices of spreaded quantity in obs. increment field
    nois   , noic    ,& ! index of (selected) spreading param. in obs incr field
    noin             ,& ! index of non-isotropic weights para. in obs incr field
    noiqc            ,& ! index of quality factor in obs. increment field
    noivlr           ,& ! index of vertical weight given to the obs. incr. due
                        ! to the distance to the obs. at the vertical interpol.
    i      , j       ,& ! loop indices in horizontal direction
!   io     , jo      ,& ! local indices of location of observation
!   istaspr, jstaspr ,& ! lower left corner of domain containing area of infl.
!   iendspr, jendspr ,& ! upper right corner of domain containing area of infl
    ixc              ,& ! loop index over influenced grid points
    nprae            ,& ! number of printed grid pts.
    jpra   , jpre    ,& ! lowest , highest meridional index of printed grid pts.
    mvb              ,& ! index of equidistant point
    joc                 ! loop index for printout

  REAL    (KIND=wp   )     ::  &
    fmultw           ,& ! 0. or 1., depending on multiple observation weighting
    rdsprs           ,& ! vertical correlation scale as by namelist
    rdsprsb, rdsprst ,& ! lower / upper vertical correlation scale
    rdsprsrb         ,& ! /  reciproke of square of lower /
    rdsprsrt         ,& ! \  upper vertical correlation scale
    rdsprir          ,& ! reciproke of vert. correl. scale used for interpolat.
    rdsprni          ,& ! square of Gauss. non-isotropy radius,'zspr(mni)'-units
    rdsprnir         ,& ! reciproke of 'rdsprni'
    zcutni           ,& ! 'rdsprni' < 'zcutni' : correlation is non-isotropic
    ztp              ,& ! value of spreading param. (z, TH) at (std atm) tropop.
    zsproni          ,& ! non-isotr. w.: TH or z at spreading surface & obs loc.
    zqts             ,& ! quotient of the tropospheric to the stratospheric mean
                        ! vertical gradient of the spreading parameter
    zfvip            ,& ! quotient for vertical interpolation
    zoiu   , zoiv    ,& ! vertically interpolated obs. increments
    zsvob            ,& ! spreading parameter at the obs. point
    zwvip            ,& ! -log of weight reduction due to vertical interpolation
    zwvip2           ,& ! vertical distance of obs. incr. level to target g.pt.
    zfvipni          ,& ! quotient for vertical interpol. from equidistant pts.
    omywth , omywp   ,& ! weighted nudging weight for a single obs. spreaded
                        ! along isentropes, resp. horizontally
    dthreshr         ,& ! diff. of stability thresholds defining the fraction of
                        ! spreading along isentropes (reciproke of item)
    sprmx            ,& ! max. allowed fraction of spreading along isentropes
    zdds2            ,& ! squared scaled horizontal distance to the obs. station
    omycorl          ,& ! local vertical correlation
    zqgnudg          ,& ! nudging coeff. relative to the coeff. for TEMP / PILOT
    qmyk1  , qmyk2   ,& ! fractions of temporal linear interpolation
    zispr  , disprr     ! define vertic. equidist. pts, for non-isotropic weight

  LOGICAL                  ::  &
    lnoniso             ! non-isotrophic horizontal correlation

  CHARACTER (LEN=15) yformat


  LOGICAL                  , SAVE ::  &
    lsprx1              ! non-zero temporal and vertical weights from 1st report


! Local (automatic) arrays:
! -------------------------
  REAL    (KIND=wp   )     ::  &
    zoicu   (istart:iend,jstart:jend) ,& ! zonal  \  wind increment at target
    zoicv   (istart:iend,jstart:jend) ,& ! merid. /  grid pt. from 1 report
    omysc   (istart:iend,jstart:jend) ,& ! nudging weight for a single obs.
    zoizu   (istart:iend,jstart:jend) ,& ! zonal  \  wind incr. at target grid
    zoizv   (istart:iend,jstart:jend) ,& ! merid. /  pt. using horiz. spreading
    omyzc   (istart:iend,jstart:jend) ,& ! nudging weight for a single obs.
                                         ! spreaded horizontally
    sprdv   (istart:iend,jstart:jend) ,& ! fraction of obs incr. to be spread
                                         ! along isentropes
    zuoem   (istart:iend,jstart:jend) ,& ! sum of the 2 spreaded obs. increments
    zvoem   (istart:iend,jstart:jend) ,& ! sum of the 2 spreaded obs. increments
    omykm   (istart:iend,jstart:jend) ,& ! sum of the 2 weights
    omy12   (istart:iend,jstart:jend) ,& ! square of the weight for first obs.
    omyc2   (istart:iend,jstart:jend)    ! square of the weight for second obs.

! Initialized variables:
! ----------------------

  LOGICAL                  , SAVE ::  &
    lfrsloc = .TRUE.    ! .TRUE if local output is to be printed the first time

!------------ End of header ----------------------------------------------------


 
!-------------------------------------------------------------------------------
! Begin Subroutine mult_spread_wind
!-------------------------------------------------------------------------------
  
!-------------------------------------------------------------------------------
!  Section 1: Preliminaries
!-------------------------------------------------------------------------------

  lnoniso = .FALSE.    ! initialize in any case for safety

  ivar    = 1
 
  IF (itim == 1) lsprx1 = lsprx

  ivay    = ivar + 1
  ivrs    = MAX( ivar - 1 , 1 )
  modespr = msprpar
! nmloi   = knoiml(ista,itim)

  irep    = 1
  IF ((itim == 2) .AND. (lsprx1)) irep = 2

  fmultw = REAL ( kwtyp(ityp) - 1, wp )

! lsigma =     ( modespr == 0)                                                 &
!         .OR. ((modespr == 1) .AND. (k <= ktp))                               &
!         .OR. ((modespr == 2) .AND. (k < ktth))
!!        .OR. ((modespr == 2) .AND. (     ((.NOT. ltopth) .AND. (k < ktopth)) &
!!                                    .OR. ((ltopth) .AND. (k <= ktp))))

  IF (lsprx) THEN

! get the spreading parameter at standard atmospheric tropopause if required,
! and the Gaussian radii of the vertical weight function

    ztp     = zstpml (ista)
    rdsprs  = SQRT( vcorls(ivar) )
    IF (modespr <= 1) THEN
      rdsprsb = rdsprs *fcorlml(ista,2*itim-1,ivrs) *zrtdgml(ista,2*itim-1,ivrs)
      rdsprst = rdsprs *fcorlml(ista,2*itim  ,ivrs) *zrtdgml(ista,2*itim  ,ivrs)
      rdsprir = c1 / SQRT( MAX( vcorls(ivar) , 0.04_wp ) )
    ELSE
      rdsprsb = rdsprs *fcorlml(ista,2*itim-1,ivrs)
      rdsprst = rdsprs *fcorlml(ista,2*itim  ,ivrs)
      rdsprir = c1 / SQRT( MAX( vcorls(ivar) , 33.0_wp ) )
    ENDIF
    rdsprsrb = c1 / rdsprsb **2
    rdsprsrt = c1 / rdsprst **2

! (Horizontally) local contributions to the weights: temporal weights, (quality
! factors), relative nudging coeff., vertical weight (spread. along model lev.)

! wtql1  =  zwtua (ista,1,ivrs)  *  zqualua (ista,1,ivrs)                      &
!         * gnudgar(ivar) / MAX( gnudg(ivar) ,epsy )

! get quantities used for the vertical weighting and non-isotropic correction,
! and the relative nudging coefficient
 
    IF (modespr <= 1) THEN
      zispr  = xezispr
      disprr = c1 / dezispr
      zcutni = vcutnit
    ELSE
      zispr  = xthispr
      disprr = c1 / dthispr
      zcutni = vcutnip
    ENDIF

    rdsprni = vcsni  (ivar)
    lnoniso = (rdsprni <= zcutni) .AND. (lnisua)
    IF (modespr == 2) rdsprni = rdsprni * zrtgkml(ista,k) **2
    IF ((modespr <= 1) .AND. (zsprml(ista,k) > ztp)) rdsprni = qs2 * rdsprni
    rdsprnir = c1 / MAX( rdsprni , epsy )
    zqgnudg = c1
    IF (kobtyp == nairep) THEN
      IF (kcdtyp /= nmodes) zqgnudg = gnudgar(ivar) / MAX( gnudg(ivar) ,epsy )
      IF (kcdtyp == nmodes) zqgnudg = gnudgms(ivar) / MAX( gnudg(ivar) ,epsy )
    ENDIF
    IF (kobtyp == nsattv) zqgnudg = gnudgtv(ivar) / MAX( gnudg(ivar) ,epsy )
    IF (kobtyp == ngps  ) zqgnudg = gnudggp       / MAX( gnudg(ivar) ,epsy )
    IF ((kobtyp == ngps ) .AND. (ivar /= 4))  zqgnudg = c0
!   lti2 = (ltiml(ista,ivrs)) .AND. (omyk1 >= epsy) .AND. (omyk2 >= epsy)

! get indices for observation increment and other arrays

    msi    = MAX( modespr , i1 )
    mni    = 3 - msi
    meqi   = mni + 1
!   isvv   = 3*(msi-1) + ivrs
    isvv   =    msi
    noix   = noiu
    noiy   = noiv
    noiqc  = noiuqc
    noivlr = noiulr
    IF (modespr <= 1) THEN
      nois = noiz
      noin = noith
      zqts = c1
    ELSE
      nois = noith
      noin = noiz
      zqts = c1 / qst
    ENDIF

  ENDIF   ! lsprx

! initialize for printout

  IF (ltiml(ista,ivrs)) THEN
    DO j   = jstaspr , jendspr
      DO i = istaspr , iendspr
        zoicu (i,j) = c0
        zoicv (i,j) = c0
        omysc (i,j) = c0
      ENDDO
    ENDDO
  ENDIF

!-------------------------------------------------------------------------------
!  Section 2: Computation of the 'spreaded observation increment' at the target
!             grid pts. by vertical interpolation of the observation increments
!             at the obs. location and application of the horizontal correla-
!             tion. Computation of the nudging weights for single observations.
!-------------------------------------------------------------------------------

! case a: Horizontal correlation non-isotropic also with full divergence
! ----------------------------------------------------------------------

  IF ((lnoniso) .AND. (lsprx)) THEN
!CDIR NODEP,VOVERTAKE,VOB
    DO ixc = 1 , ncutof(irep)
      i = icutof(ixc,irep)
      j = jcutof(ixc,irep)
! if increments are available for vertical interpolation
      IF ((klva(i,j,isvv) > -1) .AND. (klvb(i,j,isvv) > -1)) THEN
        zfvip       =   ( oiml (nois,klvb(i,j,isvv)) - zspr(i,j,msi) )         &
                      / ( oiml (nois,klvb(i,j,isvv))                           &
                         -oiml (nois,klva(i,j,isvv)) )
! vertical interpolation of observation increments
        zoiu        =     zfvip     * oiml (noix,klva(i,j,isvv))               &
                       + (c1-zfvip) * oiml (noix,klvb(i,j,isvv))
        zoiv        =     zfvip     * oiml (noiy,klva(i,j,isvv))               &
                       + (c1-zfvip) * oiml (noiy,klvb(i,j,isvv))
! 2-dim horizontal wind correlations
        zoicu (i,j) =   zcoruu(i,j) *zoiu  +  zcoruv(i,j) *zoiv
        zoicv (i,j) =   zcorvu(i,j) *zoiu  +  zcorvv(i,j) *zoiv
! value of the parameter used for non-isotropic weights at the sounding station
! and at the height where the value of the spreading parameter equals that at
! the target grid point
        zsproni     =     zfvip     * oiml (noin,klva(i,j,isvv))               &
                       + (c1-zfvip) * oiml (noin,klvb(i,j,isvv))
! reduction of weight due to vertical interpolation within vertical profile
        zsvob       =   oiml (nois  ,klva(i,j,isvv))
        zwvip       =  (  MAX( MIN( zsvob,ztp ) - zspr(i,j,msi) , c0 )         &
                        + MAX( zsvob - MAX( zspr(i,j,msi),ztp ) , c0 )*zqts )  &
                      * oiml (noivcr,klva(i,j,isvv))                           &
                     +  oiml (noivlr,klva(i,j,isvv))
        zsvob       =   oiml (nois  ,klvb(i,j,isvv))
        zwvip2      =  (  MAX( MIN( zspr(i,j,msi),ztp ) - zsvob , c0 )         &
                        + MAX( zspr(i,j,msi) - MAX( zsvob,ztp ) , c0 )*zqts )  &
                      * oiml (noivcr,klvb(i,j,isvv))                           &
                     +  oiml (noivlr,klvb(i,j,isvv))
        zwvip       =  (zwvip * zwvip2 / MAX( zwvip+zwvip2,epsy ) *rdsprir) **2
! computation of the nudging weight
        omysc (i,j) =   EXP(- (zsproni -zspr(i,j,mni))**2 * rdsprnir           &
                            - zdds(i,j,ivrs)  - zwvip )                        &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * (  zfvip     * oiml (noiqc,klva(i,j,isvv))             &
                         +(c1-zfvip) * oiml (noiqc,klvb(i,j,isvv)) )
! no obs. increment available below 'zspr(i,j)'. Take the obs. increment at the
! lowest obs. point with data of type 'ivrs' and apply vertical correlation
      ELSEIF (klva(i,j,isvv) > -1) THEN
! 2-dim horizontal wind correlations
        zoicu (i,j) =    zcoruu(i,j) * oiml (noix,klva(i,j,isvv))              &
                       + zcoruv(i,j) * oiml (noiy,klva(i,j,isvv))
        zoicv (i,j) =    zcorvu(i,j) * oiml (noix,klva(i,j,isvv))              &
                       + zcorvv(i,j) * oiml (noiy,klva(i,j,isvv))
! value of the parameter used for non-isotropic weights at the sounding station
! and at the height where the value of the spreading parameter equals that at
! the target grid point (approximate computation by vertical interpolation from
! equidistant point (in terms of spreading parameter units)) 
        mvb         =  INT( (zspr(i,j,meqi) - zispr) * disprr ) + 1
        mvb         =  MIN( MAX( mvb , i1 ) , mxispr - 1 )
        zfvipni     =    (zsdnis(mvb,msi) - zspr(i,j,msi)) * zsdnid(mvb,msi)
        zsproni     =    (c1-zfvipni) *znismq(ista,mvb)                        &
                       +  zfvipni     *znismq(ista,mvb+1)
! computation of the nudging weight:
!   spreading parameter at the obs. point
        zsvob       =   oiml (nois,klva(i,j,isvv))
!   vertical correlation: pot. temperature differences are scaled by a larger
!   constant in the stratosphere than in the troposphere (larger mean stability)
        omysc (i,j) =   EXP(- (  MAX( MIN( zsvob,ztp ) - zspr(i,j,msi) , c0 )  &
                               + MAX( zsvob - MAX( zspr(i,j,msi),ztp ) , c0 )  &
                                 * zqts ) **2  *  rdsprsrb                     &
!   further weighting for non-isotropy
                            - (zsproni -zspr(i,j,mni))**2 * rdsprnir           &
!   factor of horizontal correlation
                            - zdds(i,j,ivrs)  )                                &
!   reduction of weight within atmospheric boundary layer, relat. nudging coeff.
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
!   quality factor
                      * oiml (noiqc,klva(i,j,isvv))
! no obs. increment available above 'zspr(i,j)'
      ELSEIF (klvb(i,j,isvv) > -1) THEN
        zoicu (i,j) =    zcoruu(i,j) * oiml (noix,klvb(i,j,isvv))              &
                       + zcoruv(i,j) * oiml (noiy,klvb(i,j,isvv))
        zoicv (i,j) =    zcorvu(i,j) * oiml (noix,klvb(i,j,isvv))              &
                       + zcorvv(i,j) * oiml (noiy,klvb(i,j,isvv))
        mvb         =  INT( (zspr(i,j,meqi) - zispr) * disprr ) + 1
        mvb         =  MIN( MAX( mvb , i1 ) , mxispr - 1 )
        zfvipni     =    (zsdnis(mvb,msi) - zspr(i,j,msi)) * zsdnid(mvb,msi)
        zsproni     =    (c1-zfvipni) *znismq(ista,mvb)                        &
                       +  zfvipni     *znismq(ista,mvb+1)
        zsvob       =   oiml (nois,klvb(i,j,isvv))
        omysc (i,j) =   EXP(- (  MAX( MIN( zspr(i,j,msi),ztp ) - zsvob , c0 )  &
                               + MAX( zspr(i,j,msi) - MAX( zsvob,ztp ) , c0 )  &
                                 * zqts ) **2  *  rdsprsrt                     &
                            - (zsproni -zspr(i,j,mni))**2 * rdsprnir           &
                            - zdds(i,j,ivrs)  )                                &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * oiml (noiqc,klvb(i,j,isvv))
! no spreading along isentropes where vertical gradient of potential temperature
! is less than a threshold
      ELSE
        zoicu (i,j) =  c0
        zoicv (i,j) =  c0
        omysc (i,j) =  c0
      ENDIF
    ENDDO
  ENDIF


! case b: Spreading along isentropes
! ----------------------------------

  IF ((.NOT. lnoniso) .AND. (modespr == 2) .AND. (lsprx)) THEN
!CDIR NODEP,VOVERTAKE,VOB
    DO ixc = 1 , ncutof(irep)
      i = icutof(ixc,irep)
      j = jcutof(ixc,irep)
      IF ((klva(i,j,isvv) > -1) .AND. (klvb(i,j,isvv) > -1)) THEN
        zfvip       =   ( oiml (nois,klvb(i,j,isvv)) - zspr(i,j,msi) )         &
                      / ( oiml (nois,klvb(i,j,isvv))                           &
                         -oiml (nois,klva(i,j,isvv)) )
! vertical interpolation of observation increments
        zoiu        =     zfvip     * oiml (noix,klva(i,j,isvv))               &
                       + (c1-zfvip) * oiml (noix,klvb(i,j,isvv))
        zoiv        =     zfvip     * oiml (noiy,klva(i,j,isvv))               &
                       + (c1-zfvip) * oiml (noiy,klvb(i,j,isvv))
! 2-dim horizontal wind correlations
        zoicu (i,j) =   zcoruu(i,j) *zoiu  +  zcoruv(i,j) *zoiv
        zoicv (i,j) =   zcorvu(i,j) *zoiu  +  zcorvv(i,j) *zoiv
! reduction of weight due to vertical interpolation within vertical profile
        zsvob       =   oiml (nois  ,klva(i,j,isvv))
        zwvip       =  (  MAX( MIN( zsvob,ztp ) - zspr(i,j,msi) , c0 )         &
                        + MAX( zsvob - MAX( zspr(i,j,msi),ztp ) , c0 )*zqts )  &
                      * oiml (noivcr,klva(i,j,isvv))                           &
                     +  oiml (noivlr,klva(i,j,isvv))
        zsvob       =   oiml (nois  ,klvb(i,j,isvv))
        zwvip2      =  (  MAX( MIN( zspr(i,j,msi),ztp ) - zsvob , c0 )         &
                        + MAX( zspr(i,j,msi) - MAX( zsvob,ztp ) , c0 )*zqts )  &
                      * oiml (noivcr,klvb(i,j,isvv))                           &
                     +  oiml (noivlr,klvb(i,j,isvv))
        zwvip       =  (zwvip * zwvip2 / MAX( zwvip+zwvip2,epsy ) *rdsprir) **2
! computation of the nudging weight
        omysc (i,j) =   EXP(- zdds(i,j,ivrs)  - zwvip )                        &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * (  zfvip     * oiml (noiqc,klva(i,j,isvv))             &
                         +(c1-zfvip) * oiml (noiqc,klvb(i,j,isvv)) )
      ELSEIF (klva(i,j,isvv) > -1) THEN
! 2-dim horizontal wind correlations
        zoicu (i,j) =    zcoruu(i,j) * oiml (noix,klva(i,j,isvv))              &
                       + zcoruv(i,j) * oiml (noiy,klva(i,j,isvv))
        zoicv (i,j) =    zcorvu(i,j) * oiml (noix,klva(i,j,isvv))              &
                       + zcorvv(i,j) * oiml (noiy,klva(i,j,isvv))
        zsvob       =   oiml (nois,klva(i,j,isvv))
!   vertical correlation: pot. temperature differences are scaled by a larger
!   constant in the stratosphere than in the troposphere (larger mean stability)
        omysc (i,j) =   EXP(- (  MAX( MIN( zsvob,ztp ) - zspr(i,j,msi) , c0 )  &
                               + MAX( zsvob - MAX( zspr(i,j,msi),ztp ) , c0 )  &
                                 * zqts ) **2  *  rdsprsrb                     &
                            - zdds(i,j,ivrs)  )                                &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * oiml (noiqc,klva(i,j,isvv))
      ELSEIF (klvb(i,j,isvv) > -1) THEN
        zoicu (i,j) =    zcoruu(i,j) * oiml (noix,klvb(i,j,isvv))              &
                       + zcoruv(i,j) * oiml (noiy,klvb(i,j,isvv))
        zoicv (i,j) =    zcorvu(i,j) * oiml (noix,klvb(i,j,isvv))              &
                       + zcorvv(i,j) * oiml (noiy,klvb(i,j,isvv))
        zsvob       =   oiml (nois,klvb(i,j,isvv))
        omysc (i,j) =   EXP(- (  MAX( MIN( zspr(i,j,msi),ztp ) - zsvob , c0 )  &
                               + MAX( zspr(i,j,msi) - MAX( zsvob,ztp ) , c0 )  &
                                 * zqts ) **2  *  rdsprsrt                     &
                            - zdds(i,j,ivrs)  )                                &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * oiml (noiqc,klvb(i,j,isvv))
! no spreading along isentropes where vertical gradient of potential temperature
! is less than a threshold
      ELSE
        zoicu (i,j) =  c0
        zoicv (i,j) =  c0
        omysc (i,j) =  c0
      ENDIF
    ENDDO
  ENDIF

! for cases a and b:
! If spreading along isentropes: spreading along horizontal surfaces must also
! be applied in areas of weak stability. However, the vertical correlations
! are still expressed as function of potential temperature differences, and
! non-isotropic weights (as function of height differences) are not possible.

  IF ((modespr == 2) .AND. (lsprx)) THEN
    dthreshr= c1 / MAX( thresh1 - thresh2 , epsy )
    sprmx   = c1
    IF (k == ktth) sprmx = c05
    msip = 3 - msi
    nois = noiz
    noic = noith
!   isvv = ivrs
    isvv = 1
!CDIR NODEP,VOVERTAKE,VOB
    DO ixc = 1 , ncutof(irep)
      i = icutof(ixc,irep)
      j = jcutof(ixc,irep)
! 'sprdv': fraction of obs. increment to be spread along isentropes;
! 1-sprdv: fraction of obs. increment to be spread horizontally
      sprdv (i,j) = MAX( c0 , MIN( sprmx, (zthvg(i,j)-thresh2) *dthreshr ) )
      IF (c1-sprdv(i,j) < epsy) THEN
        zoizu (i,j) =  c0
        zoizv (i,j) =  c0
        omyzc (i,j) =  c0
      ELSEIF ((klva(i,j,isvv) > -1) .AND. (klvb(i,j,isvv) > -1)) THEN
        zfvip       =   ( oiml (nois,klvb(i,j,isvv)) - zspr(i,j,msip) )        &
                      / ( oiml (nois,klvb(i,j,isvv))                           &
                         -oiml (nois,klva(i,j,isvv)) )
        zoiu        =     zfvip     * oiml (noix,klva(i,j,isvv))               &
                       + (c1-zfvip) * oiml (noix,klvb(i,j,isvv))
        zoiv        =     zfvip     * oiml (noiy,klva(i,j,isvv))               &
                       + (c1-zfvip) * oiml (noiy,klvb(i,j,isvv))
        zoizu (i,j) =   zcoruu(i,j) *zoiu  +  zcoruv(i,j) *zoiv
        zoizv (i,j) =   zcorvu(i,j) *zoiu  +  zcorvv(i,j) *zoiv
        zsvob       =   oiml (noic  ,klva(i,j,isvv))
        zwvip       =  (  MAX( MIN( zsvob,ztp ) - zspr(i,j,msi) , c0 )         &
                        + MAX( zsvob - MAX( zspr(i,j,msi),ztp ) , c0 )*zqts )  &
                      * oiml (noivcr,klva(i,j,isvv))                           &
                     +  oiml (noivlr,klva(i,j,isvv))
        zsvob       =   oiml (noic  ,klvb(i,j,isvv))
        zwvip2      =  (  MAX( MIN( zspr(i,j,msi),ztp ) - zsvob , c0 )         &
                        + MAX( zspr(i,j,msi) - MAX( zsvob,ztp ) , c0 )*zqts )  &
                      * oiml (noivcr,klvb(i,j,isvv))                           &
                     +  oiml (noivlr,klvb(i,j,isvv))
        zwvip       =  (zwvip * zwvip2 / MAX( zwvip+zwvip2,epsy ) *rdsprir) **2
        omyzc (i,j) =   EXP(- zdds(i,j,ivrs)  - zwvip )                        &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * (  zfvip     * oiml (noiqc,klva(i,j,isvv))             &
                         +(c1-zfvip) * oiml (noiqc,klvb(i,j,isvv)) )
      ELSEIF (klva(i,j,isvv) > -1) THEN
        zoizu (i,j) =    zcoruu(i,j) * oiml (noix,klva(i,j,isvv))              &
                       + zcoruv(i,j) * oiml (noiy,klva(i,j,isvv))
        zoizv (i,j) =    zcorvu(i,j) * oiml (noix,klva(i,j,isvv))              &
                       + zcorvv(i,j) * oiml (noiy,klva(i,j,isvv))
        zsvob       =   oiml (noic,klva(i,j,isvv))
        omyzc (i,j) =   EXP(- (  MAX( MIN( zsvob,ztp ) - zspr(i,j,msi) , c0 )  &
                               + MAX( zsvob - MAX( zspr(i,j,msi),ztp ) , c0 )  &
                                 * zqts ) **2  *  rdsprsrb                     &
                            - zdds(i,j,ivrs)  )                                &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * oiml (noiqc,klva(i,j,isvv))
      ELSEIF (klvb(i,j,isvv) > -1) THEN
        zoizu (i,j) =    zcoruu(i,j) * oiml (noix,klvb(i,j,isvv))              &
                       + zcoruv(i,j) * oiml (noiy,klvb(i,j,isvv))
        zoizv (i,j) =    zcorvu(i,j) * oiml (noix,klvb(i,j,isvv))              &
                       + zcorvv(i,j) * oiml (noiy,klvb(i,j,isvv))
        zsvob       =   oiml (noic,klvb(i,j,isvv))
        omyzc (i,j) =   EXP(- (  MAX( MIN( zspr(i,j,msi),ztp ) - zsvob , c0 )  &
                               + MAX( zspr(i,j,msi) - MAX( zsvob,ztp ) , c0 )  &
                                 * zqts ) **2  *  rdsprsrt                     &
                            - zdds(i,j,ivrs)  )                                &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * oiml (noiqc,klvb(i,j,isvv))
      ELSE
        zoizu (i,j) =  c0
        zoizv (i,j) =  c0
        omyzc (i,j) =  c0
      ENDIF
    ENDDO
!   isvv = 3*(msi-1) + ivrs
    isvv =    msi

! Merging of spreading along isentropes and of horizontal spreading

!CDIR NODEP,VOVERTAKE,VOB
    DO ixc = 1 , ncutof(irep)
      i = icutof(ixc,irep)
      j = jcutof(ixc,irep)
      omywth      =     sprdv(i,j)  *omysc(i,j)
      omywp       = (c1-sprdv(i,j)) *omyzc(i,j)
      omysc (i,j) =  omywth             + omywp
      zoicu (i,j) = (omywth *zoicu(i,j) + omywp *zoizu(i,j)) / omysc(i,j)
      zoicv (i,j) = (omywth *zoicv(i,j) + omywp *zoizv(i,j)) / omysc(i,j)
    ENDDO
  ENDIF


! case c: Horizontal spreading without non-isotropic corretion
!         (this case is treated separately for speed-up)
! ------------------------------------------------------------

  IF ((modespr <= 1) .AND. (.NOT. lnoniso) .AND. (lsprx)) THEN
!CDIR NODEP,VOVERTAKE,VOB
    DO ixc = 1 , ncutof(irep)
      i = icutof(ixc,irep)
      j = jcutof(ixc,irep)
      IF ((klva(i,j,isvv) > -1) .AND. (klvb(i,j,isvv) > -1)) THEN
        zfvip       =   ( oiml (nois,klvb(i,j,isvv)) - zspr(i,j,msi) )         &
                      / ( oiml (nois,klvb(i,j,isvv))                           &
                         -oiml (nois,klva(i,j,isvv)) )
        zoiu        =     zfvip     * oiml (noix,klva(i,j,isvv))               &
                       + (c1-zfvip) * oiml (noix,klvb(i,j,isvv))
        zoiv        =     zfvip     * oiml (noiy,klva(i,j,isvv))               &
                       + (c1-zfvip) * oiml (noiy,klvb(i,j,isvv))
        zoicu (i,j) =   zcoruu(i,j) *zoiu  +  zcoruv(i,j) *zoiv
        zoicv (i,j) =   zcorvu(i,j) *zoiu  +  zcorvv(i,j) *zoiv
        zwvip       =    (oiml (nois  ,klva(i,j,isvv)) - zspr(i,j,msi))        &
                        * oiml (noivcr,klva(i,j,isvv))                         &
                       +  oiml (noivlr,klva(i,j,isvv))
        zwvip2      =  - (oiml (nois  ,klvb(i,j,isvv)) - zspr(i,j,msi))        &
                        * oiml (noivcr,klvb(i,j,isvv))                         &
                       +  oiml (noivlr,klvb(i,j,isvv))
        zwvip       =  (zwvip * zwvip2 / MAX( zwvip+zwvip2,epsy ) *rdsprir) **2
        omysc (i,j) =   EXP(- zdds(i,j,ivrs)  - zwvip )                        &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * (  zfvip     * oiml (noiqc,klva(i,j,isvv))             &
                         +(c1-zfvip) * oiml (noiqc,klvb(i,j,isvv)) )
      ELSEIF (klva(i,j,isvv) > -1) THEN
        zoicu (i,j) =    zcoruu(i,j) * oiml (noix,klva(i,j,isvv))              &
                       + zcoruv(i,j) * oiml (noiy,klva(i,j,isvv))
        zoicv (i,j) =    zcorvu(i,j) * oiml (noix,klva(i,j,isvv))              &
                       + zcorvv(i,j) * oiml (noiy,klva(i,j,isvv))
        omysc (i,j) =   EXP(- ( oiml (nois,klva(i,j,isvv))                     &
                               -zspr(i,j,msi)) **2  *  rdsprsrb                &
                            - zdds(i,j,ivrs)  )                                &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * oiml (noiqc,klva(i,j,isvv))
      ELSEIF (klvb(i,j,isvv) > -1) THEN
        zoicu (i,j) =    zcoruu(i,j) * oiml (noix,klvb(i,j,isvv))              &
                       + zcoruv(i,j) * oiml (noiy,klvb(i,j,isvv))
        zoicv (i,j) =    zcorvu(i,j) * oiml (noix,klvb(i,j,isvv))              &
                       + zcorvv(i,j) * oiml (noiy,klvb(i,j,isvv))
        omysc (i,j) =   EXP(- ( oiml (nois,klvb(i,j,isvv))                     &
                               -zspr(i,j,msi)) **2  *  rdsprsrt                &
                            - zdds(i,j,ivrs)  )                                &
                      * zablfc(i,j,ivrs,1)  *  zqgnudg                         &
                      * oiml (noiqc,klvb(i,j,isvv))
      ELSE
        zoicu (i,j) =  c0
        zoicv (i,j) =  c0
        omysc (i,j) =  c0
      ENDIF
    ENDDO
  ENDIF


! Correction for Cressman structure function instead of auto-
! regressive function (in this case, 'cnondiv' should be zero)
! ------------------------------------------------------------
  IF ((.NOT. lautorg) .AND. (lsprx)) THEN
!CDIR NODEP,VOVERTAKE,VOB
    DO ixc = 1 , ncutof(irep)
      i = icutof(ixc,irep)
      j = jcutof(ixc,irep)
      zdds2 = zdds(i,j,ivrs) **2
      omysc (i,j) = omysc(i,j) * MAX( c0 , c1-zdds2 ) / (c1 + zdds2)           &
                               * EXP( + zdds(i,j,ivrs) )
    ENDDO
  ENDIF


!-------------------------------------------------------------------------------
!  Section 3: Storage of the spreaded increments and the weights and updating
!             the quantities needed for multiple observations
!-------------------------------------------------------------------------------

  IF (itim == 2) THEN
!   IF ((lsprx) .AND. (lsprx1) .AND. (ltiml(ista,ivrs))) THEN
    IF ((ltiml(ista,ivrs)) .AND. ((lsprx) .OR. (lsprx1))) THEN
! For temporal linear interpolation selected: if the nudging weight exclusive
! temporal weight for one obs. increment 'omykl' is larger than that 'omyks' for
! the other, then the temporal weighting for this first increment will be a mix-
! ture of temporal interpolation and of the sawtooth-shaped function valid for
! single observations. The fraction of temporal interpolation is 'omyks /omykl'.
! This allows to deal with temporal interpolation of a 'complete' and an
! incomplete sounding.
      DO j   = jstaspr , jendspr
        DO i = istaspr , iendspr
!         IF ((lcutof(i,j,1)) .AND. (lcutof(i,j,2))) THEN
          IF ((lcutof(i,j,1)) .OR. (lcutof(i,j,2))) THEN
            qmyk1    = MIN( omysc(i,j) / MAX( omys1(i,j),epsy ) , c1 )
            qmyk2    = MIN( omys1(i,j) / MAX( omysc(i,j),epsy ) , c1 )
            omys1 (i,j)  = omys1(i,j) * (      qmyk1  *zwtml(ista,1,ivrs)      &
                                         + (c1-qmyk1) *zwtml(ista,3,ivrs))
            omysc (i,j)  = omysc(i,j) * (      qmyk2  *zwtml(ista,2,ivrs)      &
                                         + (c1-qmyk2) *zwtml(ista,4,ivrs))
            omykm (i,j) =  omys1(i,j) + omysc(i,j)
            zuoem (i,j) = (omykm(i,j) + fmultw) * (  omys1(i,j) *zoic1(i,j)    &
                                                   + omysc(i,j) *zoicu(i,j))
            zvoem (i,j) = (omykm(i,j) + fmultw) * (  omys1(i,j) *zoiv1(i,j)    &
                                                   + omysc(i,j) *zoicv(i,j))
            omyc2 (i,j) =  omykm(i,j) * omykm(i,j)
!         ELSEIF (lcutof(i,j,2)) THEN
!           omykm (i,j) =  omysc(i,j) * zwtml(ista,4,ivrs)
!           omyc2 (i,j) =  omykm(i,j) * omykm(i,j)
!           zuoem (i,j) = (omyc2(i,j) + fmultw *omykm(i,j)) * zoicu(i,j)
!           zvoem (i,j) = (omyc2(i,j) + fmultw *omykm(i,j)) * zoicv(i,j)
!         ELSEIF (lcutof(i,j,1)) THEN
!           omykm (i,j) =  omys1(i,j) * zwtml(ista,3,ivrs)
!           omyc2 (i,j) =  omykm(i,j) * omykm(i,j)
!           zuoem (i,j) = (omyc2(i,j) + fmultw *omykm(i,j)) * zoic1(i,j)
!           zvoem (i,j) = (omyc2(i,j) + fmultw *omykm(i,j)) * zoiv1(i,j)
!         ENDIF
!         IF ((lcutof(i,j,1)) .OR. (lcutof(i,j,2))) THEN
            omy (i,j,ivrs,ityp) = omy(i,j,ivrs,ityp) + omykm(i,j)
            om2 (i,j,ivrs,ityp) = om2(i,j,ivrs,ityp) + omyc2(i,j)
            zwi (i,j,ivar,ityp) = zwi(i,j,ivar,ityp) + zuoem(i,j)
            zwi (i,j,ivay,ityp) = zwi(i,j,ivay,ityp) + zvoem(i,j)
          ENDIF
        ENDDO
      ENDDO
    ELSE
      IF (lsprx) THEN
!CDIR NODEP,VOVERTAKE,VOB
        DO ixc = 1 , ncutof(irep)
          i = icutof(ixc,irep)
          j = jcutof(ixc,irep)
          omysc (i,j)    = omysc(i,j) * zwtml(ista,4,ivrs)
          omyc2 (i,j)    = omysc(i,j) * omysc(i,j)
          omy (i,j,ivrs,ityp) = omy(i,j,ivrs,ityp) +   omysc(i,j)
          om2 (i,j,ivrs,ityp) = om2(i,j,ivrs,ityp) +   omyc2(i,j)
          zwi (i,j,ivar,ityp) = zwi(i,j,ivar,ityp) + ( omysc(i,j) *fmultw      &
                                                      +omyc2(i,j) ) * zoicu(i,j)
          zwi (i,j,ivay,ityp) = zwi(i,j,ivay,ityp) + ( omysc(i,j) *fmultw      &
                                                      +omyc2(i,j) ) * zoicv(i,j)
        ENDDO
      ENDIF
      IF (lsprx1) THEN
!CDIR NODEP,VOVERTAKE,VOB
        DO ixc = 1 , ncutof(1)
          i = icutof(ixc,1)
          j = jcutof(ixc,1)
          omys1 (i,j)    = omys1(i,j) * zwtml(ista,3,ivrs)
          omy12 (i,j)    = omys1(i,j) * omys1(i,j)
          omy (i,j,ivrs,ityp) = omy(i,j,ivrs,ityp) +   omys1(i,j)
          om2 (i,j,ivrs,ityp) = om2(i,j,ivrs,ityp) +   omy12(i,j)
          zwi (i,j,ivar,ityp) = zwi(i,j,ivar,ityp) + ( omys1(i,j) *fmultw      &
                                                      +omy12(i,j) ) * zoic1(i,j)
          zwi (i,j,ivay,ityp) = zwi(i,j,ivay,ityp) + ( omys1(i,j) *fmultw      &
                                                      +omy12(i,j) ) * zoiv1(i,j)
        ENDDO
      ENDIF
    ENDIF

  ENDIF

! 'zoic,v1', 'omys1' need not have an index for the variable type 'ivar', 'ivrs'
! as long as the inner loop in the calling procedure is over the time index
! 'itim' and the outer loop is over 'ivrs'

  IF ((itim == 1) .AND. (lsprx)) THEN
!CDIR NODEP,VOVERTAKE,VOB
    DO ixc = 1 , ncutof(1)
      i = icutof(ixc,1)
      j = jcutof(ixc,1)
      zoic1 (i,j) = zoicu(i,j)
      zoiv1 (i,j) = zoicv(i,j)
      omys1 (i,j) = omysc(i,j)
    ENDDO
  ELSEIF ((itim == 2) .AND. (lsprx1)) THEN
!CDIR NODEP,VOVERTAKE,VOB
    DO ixc = 1 , ncutof(1)
      i = icutof(ixc,1)
      j = jcutof(ixc,1)
      zoic1 (i,j) = c0
      zoiv1 (i,j) = c0
      omys1 (i,j) = c0
    ENDDO
  ENDIF

! printout for control   (if statements must equal those at the
! --------------------    initialization of fields 'omysc', 'zoics')

  IF ((ntstep <= 1) .AND. (lwonl) .AND. (lsprx)) THEN
!   io = imladm(ista,3)   ! io,jo knowm globally
!   jo = imladm(ista,4)
!   IF (((io == ionl) .AND. (jo == jonl) .AND. (k >= ke-4)) .OR. (lfrsloc)) THEN
    IF  ((io == ionl) .AND. (jo == jonl) .AND. (k >= ke-6)) THEN
      i    = MAX( INT( MIN( io , iendspr ) ,iintegers) , istaspr )
      j    = MAX( INT( MIN( jo , jendspr ) ,iintegers) , jstaspr )
      jpra = MAX( INT( jo  -5 ,iintegers) , jstaspr )
      jpre = MIN( INT( jpra+9 ,iintegers) , jendspr )
      jpra = MAX( INT( jpre-9 ,iintegers) , jstaspr )
      IF (klva(i,j,isvv) > -1) THEN
        IF (modespr == 2) nois = noith
        zsvob   = oiml (nois,klva(i,j,isvv))
        omycorl = EXP( -(zsvob - zspr(i,j,msi))**2 * rdsprsrb )
        WRITE( nupr,'(''spr_uvz: k='',I2,'', omysc, zoic''                     &
                    &,2F8.2, 2F7.1, 2F8.5)' )    k, rdsprsb, rdsprst           &
             , zsvob, zspr(i,j,msi), omycorl, omysc(i,j)
      ENDIF
      j     = jo   - jpra + 1
      nprae = jpre - jpra + 1
      WRITE( yformat, '(''('',I2,''(6X,L1),2I4)'')' ) nprae
      WRITE( nupr, yformat )  (lcutof(i,joc,irep), joc=jpra,jpre), ivar, j
      WRITE( yformat, '(''('',I2,''F7.4,I4    )'')' ) nprae
      WRITE( nupr, yformat )  (omysc (i,joc     ), joc=jpra,jpre), ivar
      WRITE( nupr, yformat )  (zoicu (i,joc     ), joc=jpra,jpre), ivar
      WRITE( nupr, yformat )  (zoicv (i,joc     ), joc=jpra,jpre), ivar
      IF (modespr == 2) WRITE( nupr,yformat ) (sprdv(i,joc),joc=jpra,jpre),ivar
      WRITE( yformat, '(''('',I2,''F7.2,I4    )'')' ) nprae
      WRITE( nupr, yformat )  (zthvg (i,joc     ), joc=jpra,jpre), ivar
      WRITE( nupr, yformat )  (zdds  (i,joc,ivrs), joc=jpra,jpre), ivar
      WRITE( yformat, '(''('',I2,''F7.1,I4    )'')' ) nprae
      WRITE( nupr, yformat )  (zspr  (i,joc,msi ), joc=jpra,jpre), ivar
      lfrsloc = .FALSE.
    ENDIF
    IF (      (((io == ionl) .AND. (jo == jonl)) .OR. (MOD(ista ,4) == 0))     &
        .AND. (k == ke-5) .AND. (lsprx1))                                      &
      WRITE( nupr,'(''lti_mult '',A ,'', itim'',I2,'', zoicu/v1,om1'',2F7.2    &
                  &,F7.4)' )                                                   &
             ystid, itim, zoic1(ionl,jonl), zoiv1(ionl,jonl), omys1(ionl,jonl)
  ENDIF


!-------------------------------------------------------------------------------
! End of module procedure mult_spread_wind
!-------------------------------------------------------------------------------

END SUBROUTINE mult_spread_wind


!===============================================================================

ELEMENTAL REAL (KIND=wp) FUNCTION rmod  ( ztti, ziv )
  !----------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  ztti, ziv
  !  uses parameter:                           epsy
  !---------------------------------------------------------------------------
  ! MOD function for positive REALS
  !---------------------------------------------------------------------------
  !
  rmod  =  ztti  -  ziv * REAL(INT( ztti/ziv + epsy ), wp) + epsy
  !
END FUNCTION rmod

!-------------------------------------------------------------------------------

END MODULE src_mult_spread
