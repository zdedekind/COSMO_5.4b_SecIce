!+ Source module for gathering the required local information for the nudging
!-------------------------------------------------------------------------------

MODULE src_gather_info

!-------------------------------------------------------------------------------
!
! Description:
!   The module "src_gather_info" performs all the communication between PE's
!   related to the computation of the local information (including the obser-
!   vation increments) and the spreading of the observational information to
!   model grid points. In particular, it gathers all the local information
!   from all observations required for the spreading, and computes from the
!   prognostic fields all fields, which are required for the spreading at any
!   model level except the lowest level.
!
! Method:
!   This module contains the following module procedures:
!    - gather_local_info  : collect all local information on obs from all PE's
!    - gather_spread_aux  : compute auxiliary quantities for spreading
!    - gather_nudge_aux   : compute auxiliary fields for nudging
!    - gather_varia       : auxiliary actions mostly related to communication
!   Driving routine is "organize_nudging" of module "src_obs_use_org.f90".
!
!   This module also contains an elemental function, formerly statement funct.:
!    - rmod               : MOD function for positive reals
!
!   Note: This module contains MPI-based communication between PE's.
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! History:
! Version    Date       Name
! -------    ----       ----
! 1.13       1998/10/22 Christoph Schraff
!  Initial release
! 1.15       1998/11/02 Christoph Schraff
!  Global declaration of allocatable arrays moved to individual procedures.
! 1.19       1998/12/11 Christoph Schraff
!  Introduction of a verification mode. Revised argument list to 'gather_all'.
!  ANSI violations removed.
! 1.20       1999/01/07 Guenther Doms
!  Renaming of some global variables
! 1.27       1999/03/29 Christoph Schraff
!  Modified control output format for the threshold quality control.
! 1.29       1999/05/11 Ulrich Schaettler
!  Adapted interfaces to utility-modules.
! 1.31       1999/07/01 Christoph Schraff
!  MPI calls related to fixed-length arrays replaced by parallel_utility calls;
!  quantities related to MPI communicator 'icomm_world' removed; input to
!  'global_values' adjusted.
! 1.34       1999/12/10 Ulrich Schaettler
!  Changed calls to timing routines
! 1.38       2000/04/06 Christoph Schraff
!  Quality control output for hydrostatic height; adaptation (output only) to
!  new surface layer parameterization; corrected call of printing OIR bounds.
! 1.39       2000/05/03 Ulrich Schaettler
!  Changed names for variables concerned to latitude or longitude and
!  adapted calls to get_timings.
! 1.42       2000/06/19 Christoph Schraff
!  Quality control output for multi-level checks.
! 2.4        2001/01/29 Christoph Schraff
!  Procedure 'gather_nudge_aux' added: determination of specified areas around
!  convectively precipitating grid points. Addition of 'nactio=0' actions for
!  'gather_local_info' and 'gather_spread_aux' to enable the spatial consistency
!  check for surface pressure observations.
! 2.5        2001/06/01 Christoph Schraff
!  Bug correction at array de-allocation. Savety check at array allocations.
!  Re-organisation of 'gather_varia'. Replacement of erroneous variable top_con
!  in order to make results independent from domain decomposition.
! 2.11       2001/09/28 Christoph Schraff
!  Modification for porting to IBM
! 2.12       2001/11/07 Christoph Schraff
!  Bug corr: 'isrtpqc' deallocated also if verification is done without nudging.
! 2.13       2002/01/18 Christoph Schraff
!  Bug corr: limit at computing 'ntimysu', for printing only.
! 2.19       2002/10/24 Christoph Schraff
!  Adaptions for geostrophic wind correction and for surface-level data nudging.
! 3.3        2003/04/22 Christoph Schraff
!  Quality control output for GPS IWV data.
! 3.5        2003/09/02 Ulrich Schaettler
!  Adapted interface for routine exchg_boundaries
! 3.6        2003/12/11 Christoph Schraff
!  Improved statistics on assimilated observation increments for ps, and
!  statistics on multi-level reports discriminated for observation types.
!  Additional print-out on the processing of GPS-derived IWV observations.
!  Ensuring the printing of error messages at calls of model_abort.
! 3.7        2004/02/18 Ulrich Schaettler
!  Renamed cphi by crlat (due to global renaming)
! 3.18       2006/03/03 Christoph Schraff
!  Adaptations including distribution of information for IWV spatial consistency
!  check. Adaptations for new multiple weighting option. New QC messages.
!  Caution messages on insufficient ODR size to specific file unit 'nucautn'.
!  Flushing of other output files.
! V3_23        2007/03/30 Ulrich Schaettler
!  Introduced ldump_ascii for flushing the ASCII files
! V4_1         2007/12/04 Christoph Schraff
!  Use of 'startlat' replaced by 'startlat_tot' to ensure reproducibility.
! V4_4         2008/07/16 Ulrich Schaettler
!  Eliminated timing variables which are unused
!  Enclosed barrier calls with ltime_barrier
!  Adapted interface of get_timings
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_17        2011/02/24 Ulrich Blahak
!  Adapted interface of exchg_boundaries; corrected kzdims(1:20) -> kzdims(1:24);
!  eliminated my_peri_neigh
! V4_22        2012/01/31 Christoph Schraff
!  - Increased efficiency by use of quicksort instead of bubblesort.
!  - Introduction of 'psaigeo', 'ztpgeo', 'kobtysu' (for geostrophic pressure
!    increments) and of 'kcdtyml', 'kcdtyua', 'kobtysu', 'kcdtysu' (e.g. for
!    for separate net increments for different set of observing systems).
!  - Modified use of variable 'zablfc'.
!  - Some variables moved from 'data_nudge_all' to 'data_obs_lib_cosmo'.
!  - Open and close of file yuverif moved to module 'src_obs_print_vof'.
!  - Keep file yuquctl closed except when is it used.
!  - Additional checks and writing related to insufficient ODR size ('nexcess').
!  - Bug fix for sorting co-located multi-level aircraft reports.
!  - Bug fix for printing number of GPS stations with active obs increments.
!  - Direct calls of MPI routines (mpi_allgatherv) replaced by call of
!    'gather_all'.
! V4_25        2012/09/28 Ulrich Schaettler
!  Introduced nexch_tag for MPI boundary exchange tag to replace ntstep
! V4_28        2013/07/12 Christoph Schraff, Ulrich Schaettler
!  Messages added for ps-obs rejected due to check against lateral BC fields.
!  Statement function 'rmod' replaced by elemental function. Improved comments.
!  Use parameters for vertical grid from module vgrid_refatm_utils (US)
! V4_29        2013/10/04 Ulrich Schaettler
!  For the COSMO-Model only use vcoord from vgrid_refatm_utils
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
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
    istartpar,    & ! start index for computations in the parallel program
    iendpar,      & ! end index for computations in the parallel program

!   meridional direction
    jstart,       & ! start index for the forecast of w, t, qd, qw and pp
    jend,         & ! end index for the forecast of w, t, qd, qw and pp
    jstartpar,    & ! start index for computations in the parallel program
    jendpar,      & ! end index for computations in the parallel program

! 4. constants for the horizontal rotated grid and related variables
! ------------------------------------------------------------------

    dlon,         & ! grid point distance in zonal direction (in degrees)
    dlat,         & ! grid point distance in meridional direction (in degrees)
    startlat_tot, & ! transformed latitude of the lower left grid point
                    ! of the total domain (in degrees, N>0)
    startlat,     & ! transformed latitude of the lower left grid point
                    ! of this subdomain (in degrees, N>0)
    degrad,       & ! factor for transforming degree to rad

! 4. variables for the time discretization and related variables
! --------------------------------------------------------------

    dt,           & ! long time-step
    dtdeh           ! dt / 3600 seconds

! end of data_modelconfig

!-------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 2. physical constants and related variables
! -------------------------------------------

    rdocp,        & ! r_d / cp_d
    r_earth         ! mean radius of the earth

! end of data_constants

!-------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------

    p0         ,    & ! reference pressure at main levels             ( Pa)
    hhl        ,    & ! geometical height of half model levels        ( m )

! 2. external parameter fields                                        (unit)
! ----------------------------

    rmy        ,    & ! Davis-parameter for boundary relaxation         --
    crlat      ,    & ! cosine of transformed latitude
    fr_land    ,    & ! fraction of land in a grid element

! 3. prognostic variables                                             (unit)
! -----------------------

    u          ,    & ! zonal wind speed                              ( m/s )
    v          ,    & ! meridional wind speed                         ( m/s )
    t          ,    & ! temperature                                   (  k  )
    pp         ,    & ! deviation from the reference pressure         ( pa  )

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------

    clc_con    ,    & ! cloud cover due to convection                   --
!   top_con    ,    & ! level index of convective cloud top             --
    qcvg_con   ,    & ! moisture convergence for Kuo-type closure     (    1/s)
!   mflx_con   ,    & ! cloud base massflux                           (kg/m2*s)
!   prr_con    ,    & ! precipitation rate of rain, convective        (kg/m2*s)
!   prs_con    ,    & ! precipitation rate of snow, convective        (kg/m2*s)
    prne_con          ! precipitation rate, no evaporat., convective  (kg/m2*s)

! end of data_fields

!-------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------

    nstop,        & ! last time step of the forecast
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
    nnew,         & ! corresponds to ntstep + 1

! 3. controlling the physics
! --------------------------

    lconv,        & ! forecast with convection
    nincconv,     & ! time step increment for running the convection scheme
    itype_tran,   & ! type of surface-atmosphere transfer

! 7. additional control variables
! -------------------------------

    ltime,        & ! detailled timings of the program are given
    lreproduce,   & ! the results are reproducible in parallel mode
    ldump_ascii,  & ! for flushing (close and re-open) the ASCII files
    lperi_x,      & ! if lartif_data=.TRUE.: periodic boundary conditions in x-dir.
                    ! or with Davies conditions (.FALSE.)
    lperi_y,      & ! if lartif_data=.TRUE.: periodic boundary conditions in y-dir.
                    ! or with Davies conditions (.FALSE.)
    l2dim           ! 2 dimensional runs

! end of data_runcontrol 

!-------------------------------------------------------------------------------

USE data_parallel,      ONLY :  &
    nprocx,          & ! number of processors in x-direction
    nprocy,          & ! number of processors in y-direction
    num_compute,   & ! number of compute PEs
    nboundlines,   & ! number of boundary lines of the domain for which
                     ! no forecast is computed = overlapping boundary
                     ! lines of the subdomains
    ldatatypes,      & ! if .TRUE.: use MPI-Datatypes for some communications
    ltime_barrier,   & ! if .TRUE.: use additional barriers for determining the
                       ! load-imbalance
    ncomm_type,      & ! type of communication
    my_cart_id,    & ! rank of this subdomain in the cartesian communicator
!   my_cart_pos,   & ! position of this subdomain in the cartesian grid
                     ! in x- and y-direction
    my_cart_neigh, & ! neighbors of this subdomain in the cartesian grid
!   my_peri_neigh, & ! periodic neighbors of this subdomain in the periodic grid
                     ! if it lies at the boundary of the total domain
    isubpos,       & ! positions of the subdomains in the total domain. Given
                     ! are the i- and the j-indices of the lower left and the
                     ! upper right grid point in the order
                     !                  i_ll, j_ll, i_ur, j_ur.
                     ! Only the interior of the domains are considered, not
                     ! the boundary lines.
    iexch_req,     & ! stores the sends requests for the neighbor-exchange
    icomm_cart,    & ! communicator for the virtual cartesian topology
    imp_reals,     & ! determines the correct REAL type used in the model
                     ! for MPI
    nexch_tag,     & ! tag to be used for MPI boundary exchange
                     !  (in calls to exchg_boundaries)
    sendbuf,       & ! sending buffer for boundary exchange:
                     ! 1-4 are used for sending, 5-8 are used for receiving
    isendbuflen,   & ! length of one column of sendbuf
    imp_integers     ! determines the correct INTEGER type used in the model
                     ! for MPI

! end of data_parallel

!-------------------------------------------------------------------------------

USE data_nudge_all , ONLY :   &

! 1. Parameters and related variables
! -----------------------------------

    lwonl        ,& ! .TRUE if grid pt. (ionl ,jonl ) lies in the local domain
    lwonl2       ,& ! .TRUE if grid pt. (ionl2,jonl2) lies in the local domain
    lcroot       ,& ! .TRUE if (my_cart_id  == 0) (for print on std. output)
    lfirst       ,& ! .TRUE if 'organize_nudging' is called the first time
    onl_cart_id  ,& ! 'my_cart_id'  of node with area containing (ionl,jonl)
    acthr        ,& ! actual forecast hour
    ltnudge      ,& ! nudging (not only verification) at current timestep

! 2. Namelist variables controlling the data assimilation
! -------------------------------------------------------

    lnudge       ,& ! .f. : on - off switch for nudging
    lverif       ,& ! .f. : on - off switch for verification
    nwtyp        ,& ! 1   : if > 1 then compute net obs. increments for 'nwtyp'
                    !       different sets of observing systems separately
    nudgend      ,& ! 0   : end of nudging period in timesteps
    nverend      ,& ! 0   : end of verification period in timesteps
    tconbox      ,& ! 6*dt: timestep [s] for computing analysis increments
                    !       (i.e. time box of constant analysis increments)
    luvgcor      ,& ! .t. : .t. ==> geostrophic wind correction applied
    mpsgcor      ,& ! 1   : mode to apply geostrophic pressure correction
    khumbal      ,& ! 100 : range around convectively precipitating grid pts, at
                    !       which specified (not relative) humidity is preserved
    ptpstop      ,& ! 400.: pressure [hPa] at upper boundary of the temperature
                    !       correction for 'surface' pressure nudging
    msprpar      ,& ! 1   : switch specifying the surfaces along which observat.
                    !       increments of upper-air data are (primarily) spread
    msprpsu      ,& ! 0   : switch specifying the surface along which surface-
                    !       level data increments are primarily spreaded
    wablua       ,& ! 4*  1. : factor to weights within the ABL
                    !          (atmospheric boundary layer, i.e. mixed layer)
    wablsu       ,& ! 4*  1. : factor to weights above the ABL
                    !          (atmospheric boundary layer, i.e. mixed layer)
    topobs       ,& !  849., : threshold [hPa]: above this level (p < topobs),
                    ! 1099.,   only obs. increments at model levels are used,
                    !  799.,699.  i.e. obs. incr. at obs. levels are not used
    botmod          ! 3*1099.: threshold [hPa]: below this level (p > botmod),
                    !    899.  only obs. increments at obs. levels are used,
                    !          i.e. obs. incr. at model levels are not computed
 
USE data_nudge_all , ONLY :   &

    maxmlo       ,& !  300   : max. number of multi-level reports in the ODR
    maxsgo       ,& ! 3000   : max. number of (surface-level and upper-air)
                    !          single-level reports in the ODR
    maxuso       ,& !  900   : max. number of upper-air single-level rep. in ODR
    maxgpo       ,& ! 3000   : max. number of GPS reports within total domain
    maxtvo       ,& !    1   : max. number of sat retrievals within total domain
    ionl         ,& !  167   : / grid point coordinates
    jonl         ,& !  103   : \ for standard output on nudging
    ionl2        ,& !  167   : / 2nd grid pt coordinates
    jonl2        ,& !  103   : \ for other standard output on nudging

! 5.0  Full model fields on Arakawa A grid
! ----------------------------------------
    a_u          ,& ! zonal wind speed            on Arakawa A grid ( m/s )
    a_v          ,& ! meridional wind speed       on Arakawa A grid ( m/s )
    a_p          ,& ! pressure (full value)       on main levels    ( Pa  )
    a_z          ,& ! geometrical height          of main levels    (  m  )             

! 5.2  Areal arrays
! -----------------
    zconbas      ,& ! height of base of convectively instable region
    zcontop         ! height of top  of convectively instable region

! end of data_nudge_all

!-------------------------------------------------------------------------------

USE data_obs_lib_cosmo, ONLY :  &

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

! 4. I/O device numbers for obs processing / nudging and file names
! -----------------------------------------------------------------

    ! file names
    yucautn    ,& ! caution messages if too many obs for ODR size
    yuquctl    ,& ! data rejected by threshold QC at obs time
    yustats    ,& ! statistics of processed reports
!   yuverif    ,& ! VOF (output): verification observation file (obs.
                  !      incl. quality control flag for verification)
    yuprint    ,& ! all the remaining information
    ! device numbers
    nucautn    ,& ! caution messages if too many obs for current ODR size
    nuqc       ,& ! data rejected by threshold quality control at obs time
    nustat     ,& ! statistics of processed reports
!   nuverif    ,& ! VOF (output): verification observation file (observations
                  !   incl. quality control flag for verification)
    nupr          ! all the remaining information

USE data_obs_lib_cosmo, ONLY :  &

! 5. CMA observation type and code type numbers
! ---------------------------------------------

    mxobtp     ,& ! number of observation types
    nairep     ,& ! AIRCRAFT reports
    ntemp      ,& ! TEMP  reports
    npilot     ,& ! PILOT reports
    ngps       ,& ! GPS   reports
    nsattv     ,& ! SATELLITE retrievals
    nscatt     ,& ! SCATT reports
    nldpcd     ,& !   pilot land   report
    nshpcd     ,& !   pilot ship   report
    nmopcd     ,& !   pilot mobile report
    nwp_eu     ,& !   wind profiler report (European)
    nra_eu     ,& !   sodar/rass report (European)
    npr_us     ,& !   wind profiler/rass report (USA)
    nravad     ,& !   radar vad wind profile report
    nsmsg1     ,& !   MSG_1 (METEOSAT-8) satellite retrieval
    nsmsg2     ,& !   MSG_2 (METEOSAT-9) satellite retrieval
    nnoa15     ,& !   NOAA15 satellite retrieval
    nnoa16     ,& !   NOAA16 satellite retrieval
    nnoa17     ,& !   NOAA17 satellite retrieval
    nnoa18        !   NOAA18 satellite retrieval

! end of data_obs_lib_cosmo

!-------------------------------------------------------------------------------

USE data_nudge_gather , ONLY :   &

! 1. Parameters and general variables
! -----------------------------------

    ilstidg      ,& ! char. length used for gathering the station ID
    ilstidp      ,& ! char. length used for printing the station ID
                    ! Note: (ilstid >= ilstidg >= ilstidp) !!!
    ystid        ,& ! obs. station identity to be printed
    zalllow      ,& ! smallest allowed height within model domain
    zallhig      ,& ! largest  allowed height within model domain
    p0r          ,& ! reference pressure for potential temperature
    mxispr       ,& !     number of equidistant vertical points
    xezispr      ,& !   / parameters used to
    dezispr      ,& !  /  define the vertically
    sezispr      ,& ! <   equidistant points
    xthispr      ,& !  \  (for non-isotropic
    dthispr         !   \ horizontal correlations)

USE data_nudge_gather , ONLY :   &

! 2. The required local information on the observations and their location
! ------------------------------------------------------------------------

! local information on multi-level reports
    maxmloi_tot  ,& ! length of arrays with local info on multi-level reports
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
    ystidsu         ! station identity of surface-level station

USE data_nudge_gather , ONLY :   &

! local information on 'surface' (i.e. on lowest model level) pressure data
    zoips        ,& ! observation increments
    omykps       ,& ! local part of nudging weight, i.e. exclud. horiz. weight
    r1ifps       ,& ! horizontal correlation scale
    zmassps      ,& ! total mass affected by the 'temperature correction'
    ztdpps       ,& ! 'temperature / pressure' at obs. point
    qcimps       ,& ! model pressure interpolated to obs. point
    qctmps       ,& ! time of all observations
    qcqfps       ,& ! weight factor related to obs. quality (representiveness)
    zoips_b      ,& ! increments: observations minus lateral boundary fields
    iops         ,& ! longitudinal index of station location on local domain
    jops         ,& ! latitudinal  index of station location on local domain
    iops_tot     ,& ! longitudinal index of station location on global domain
    jops_tot     ,& ! latitudinal  index of station location on global domain
    iqclps       ,& ! indicator for current obs. to undergo spatial check now
    iqcfps       ,& ! indicator for obs. to have passed latest threshold QC
    iqcnps       ,& ! index of local administrator of (local) ODR
    ltips        ,& ! .TRUE if temporal linear interpol. for surface pressure
    lmlps        ,& ! .TRUE if datum is derived from multi-level report
    ystidps         ! station identity of 'surface' pressure station

USE data_nudge_gather , ONLY :   &

! local information on integrated water vapour increments (for humidity QC)
    zoiciv       ,& ! observation increments
    zqcfiv       ,& ! local part of nudging weight, i.e. exclud. horiz. weight
    zmodiv       ,& ! model integrated water vapour IWV
    zsativ       ,& ! IWV of saturated model temperature profile
    zactiv       ,& ! observation time
    ioiv         ,& ! longitudinal index of station location on local domain
    joiv         ,& ! latitudinal  index of station location on local domain
    ioiv_tot     ,& ! longitudinal index of station location on global domain
    joiv_tot     ,& ! latitudinal  index of station location on global domain
    iqcliv       ,& ! indicator for current obs. to undergo spatial check now
    iqcfiv       ,& ! indicator for obs. to have passed latest threshold QC
    kiobiv       ,& ! index of local administrator of (local) ODR
    kioiiv       ,& ! index of local information array for multi-level data
    ktypiv       ,& ! observation type
    ystidiv      ,& ! station identity of integrated water vapour station

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

! 4. (Local) information gathered by 1 (or 2) nodes for printing for control
! --------------------------------------------------------------------------

    oyqc         ,& ! on data rejected by the threshold quality control
    myqc         ,& ! on data rejected by the threshold quality control
    yyqc         ,& ! on data rejected by the threshold quality control
    oysu         ,& ! on interpol. of surface-level data to lowest model level
    yysu            ! on interpol. of surface-level data to lowest model level

USE data_nudge_gather , ONLY :   &

! 5. Variables defining the size of the arrays containing the local information
! -----------------------------------------------------------------------------

    maxoil       ,& ! max.  number of vertical levels in the OIR 'oiml'
    maxpso       ,& ! max.  number of active 'surface' pressure data
    maxivq       ,& ! max.  number of IWV reports used for spatial check
    maxqcp       ,& ! max.  number of rejected data to be printed per timestep
    maxysu       ,& ! max.  number of extrapolated printed surface-level reports
    nmloit       ,& ! total number of active profiles of obs. incr. in the OIR
    nmltot       ,& ! total number of active multi-level stations
    nuatot       ,& ! total number of active upper-air single-level stations
    nsutot       ,& ! total number of active surface-level stations
    npstot       ,& ! total number of active surface-pressure stations
    nivtot       ,& ! total number of IWV reports used for spatial check
    ntotqc       ,& ! total number of rejected data to be printed per timestep
    ntotys       ,& ! total number of printed interpol. surface-level reports
    nuaex        ,& ! number of local upper-air obs not used due to ODR size 
    ktopsu       ,& ! uppermost model level possibly influenced by surface obs
    lnissu       ,& ! non-isotrophic correlations for surface-level data
    lnisua       ,& ! non-isotrophic correlat. for upper-air single-lev. data

! 6. Geometrics and variables used for spatial consistency checks
! ---------------------------------------------------------------

    pyyd         ,& ! latitudinal (merid.) distance on tangent cone projection
    pxxd2        ,& ! square of zonal distance on tangent cone projection
    pxsalpa      ,& ! factor used for distances on tangent cone projection
    isrtpqc      ,& ! (sorted) list of stations with 'surface' pressure data
    isrtvqc         ! (sorted) list of stations with IWV

! end of data_nudge_gather

!-------------------------------------------------------------------------------

USE data_nudge_spread , ONLY :   &

! 1. Analysis increment fields (to be kept in the long-term storage)
! ------------------------------------------------------------------

    uanai        ,& ! analysis increments of zonal wind
    vanai        ,& ! analysis increments of meridional wind
    tanai        ,& ! analysis increments of temperature
    qanai        ,& ! analysis increments of specific humidity
    psanai       ,& ! analysis increments of pressure at lowest model level
    psaigeo      ,& ! analysis increments of pressure at the lowest model level
                    !   from geostrophic balancing of wind analysis increments
    taips        ,& ! temperature analysis increments due to pressure nudging

! 2. Analysis increment fields
! ----------------------------

    zpai         ,& ! analysis incr. of pressure, condens./evapor. not included
    zroi         ,& ! analysis incr. of density , condens./evapor. not included
    ztwips       ,& ! temperature increments implied (physic'ly or statistic.)
                    ! by pressure analysis increments at the lowest model level
    ztpgeo       ,& ! temperature increments implied (physic'ly or statistic.)
                    ! by geostrophic (surface) pressure analysis increments

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
    c2alpa       ,& !  / further factors used to compute
    scalpa       ,& ! /  the 2-dim. horizontal wind correlations
    xcalpa       ,& ! \  and the horizontal distances
    xsalpa       ,& !  \ on the tangent cone projection

! 6. Further horizontal input fields for the spreading of obs. increments
! -----------------------------------------------------------------------

    zlop         ,& ! log( pressure )
    zeklop       ,& ! exp( R/cp * log(p) )
    zthvg        ,& ! vertical gradient of potential temperature
    zpk          ,& ! pressure
    zspr         ,& ! spreading parameter , param. def. non-isotropic weights
    zdds         ,& ! scaled horizontal distance betw. obs. and target grid pt
    zcoruu       ,& ! zonal  wind - zonal  wind correlation  \  (without
    zcoruv       ,& ! zonal  wind - merid. wind correlation   \  EXP( -zdds )
    zcorvu       ,& ! merid. wind - zonal  wind correlation   /  -term )
    zcorvv       ,& ! merid. wind - merid. wind correlation  /
    zablpo       ,& ! =1 if grid pt. is (fully) above the ABL, =0 otherwise
    zablfc       ,& ! reduced weighting inside/above ABL for upper-air/sfc obs
    rhscal          ! data density-dependent factor to ps horiz. correl. scale

USE data_nudge_spread , ONLY :   &

! 7. Further fields and variables used for or during the spreading
! ----------------------------------------------------------------

    zsdnis       ,& ! equidistant pts. used for non-isotropic horiz. correlat.
    zsdnid       ,& ! distance between two adjacent 'equidistant' pts.
    gppkmi       ,& ! convertor for zonal dist. from 'km' to rotated grid pts.
    gppkmj       ,& ! convertor for merid. dist. from 'km' to rotated grid pts
    lcutof       ,& ! .TRUE if grid pt. is within area of influence of report
    icutof       ,& ! x-coordinate of grid pt. within area of influence
    jcutof       ,& ! x-coordinate of grid pt. within area of influence
    icutmp       ,& ! x-coordinate of grid pt. within area of influence, tmp
    jcutmp       ,& ! x-coordinate of grid pt. within area of influence, tmp
    idimcut      ,& ! dimension of 'icutof', 'jcutof', etc.

! 8. Indices and lists of indices
! -------------------------------

    isortps      ,& ! (sorted) list of stations with 'surface' pressure data
    isortsu      ,& ! (sorted) list of stations with surface-level data
    isortua      ,& ! (sorted) list of stations with upper-air single-level data
    isortml      ,& ! (sorted) list of stations with multi-level data
    jrs_tot      ,& ! /  index range for convertor
    jre_tot      ,& ! \  for zonal distances 'gppkmi'
    ista         ,& ! index of observing station
!   kml250       ,& ! index of full model level corresponding to about 250 hPa

! 9. Observation types in weighted increment arrays 'zwi'
! -------------------------------------------------------

    mxotyp       ,& ! number of observation types

! 10. Other variables
! -------------------

    wvgeo           ! vertical weights to the geostrophic wind correction

! end of data_nudge_spread

!-------------------------------------------------------------------------------

USE environment,              ONLY :  &
      comm_barrier,  & ! explicit synchronization point
      model_abort,   & ! aborts the program in case of errors
      exchg_boundaries ! performs the boundary exchange betw. neighboring PE's

!-------------------------------------------------------------------------------

USE parallel_utilities,       ONLY :  &
      global_values, & ! computes global values by operating on local arrays
      gather_values, & ! gathers a set of values from all nod. to 1 or all nodes
      gather_all       ! gathers at all nodes a set of arrays from all nodes

!-------------------------------------------------------------------------------

USE vgrid_refatm_utils, ONLY :   &
    vcoord

!-------------------------------------------------------------------------------

USE time_utilities,           ONLY :  get_timings, i_communications_nud,     &
      i_barrier_waiting_nud, i_local_info, i_spread_ps_pre

!-------------------------------------------------------------------------------

USE utilities,                ONLY :  &
      sortrx           ! permutation vector for (quick-)sorting a given array

!-------------------------------------------------------------------------------
! Local declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE

! 1. Variables
! ------------

  INTEGER (KIND=iintegers)                , PRIVATE ::  &
    nmlall         ,& ! vertical profiles (stations)   /  number of
    nuaall         ,& ! upper-air single-level        /   active
    nsuall         ,& ! surface-level                /    observation
    npsall         ,& ! surface pressure             \    increments
    nqcall         ,& ! threshold quality control     \   from all
    nivall            ! IWV spatial consistency check  \  nodes

  REAL    (KIND = wp)                     , PRIVATE :: &
    zmaxrmy  (3)   ,& ! max. lateral boundary relaxation coefficients
    zrealdiff         ! elapsed time since latest measuring


!===============================================================================

!-------------------------------------------------------------------------------
! Public and Private Subroutines
!-------------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------------
!+ Module procedure to collect all local information on obs from all PE's
!-------------------------------------------------------------------------------

SUBROUTINE gather_local_info ( nactio , nexcess )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure gathers (collects) all the required, previously
!   computed local information on all observations (of all report types)
!   from each node, and distributes this information to all nodes.
!   Thus at the beginning, each node only has the local information from all
!   observations in his own (sub-)domain, and at the end, each node has the
!   local information from all observations in the whole of the model domain.
!   Optionally, the time needed for this communication is measured.
!
! Method:
!   For each type of reports separately (surface pressure, other surface-level
!   data, upper-air single-level data, multi-level data), the local information
!   is written into a buffer for each data type (integer, real, character).
!   Data of type logical are first converted to integer.
!   These buffers are the input for procedure 'gather_all', which distributes
!   the local contents of the buffers to all nodes (by calling MPI-routines).
!   Finally the contents of the buffers are written back to the data records
!   for further processing (spreading of observation increments).
!   For multi-level data, the observation increment record is distributed in
!   a separate step by direct calls of MPI-routines).
!
! Written by        :  Christoph Schraff, DWD  (original version: 06.07.97)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and 
!     Documenting Exchangeable Fortran 90 Code".
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!
! Modules used:       These are declared in the module declaration section
!
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER (KIND=iintegers), INTENT (IN)         ::       &
    nactio              ! determines, which quantities are to be computed (when)

  INTEGER (KIND=iintegers), INTENT (INOUT), OPTIONAL    ::       &
    nexcess (5)         ! number of reports (obs incr.) in excess of array size

! Local parameters: None
! ----------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    i1 = 1_iintegers    !

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    mxobs            ,& ! max. number of obs. of current obs. type
    mxiobty          ,& ! no. of integer elements to be gathered per obs & level
    mxrobty          ,& ! number of real elements to be gathered per obs & level
    mxyobty          ,& ! number of charac. elements to be gathered per obs.
    ichlen           ,& ! length of characters elements to be gathered
    ivrs             ,& ! index of observed quantity: 1=(u,v); 2=T; 3=RH
    itim             ,& ! (time) index over obs. at one observing station
    mxvrs            ,& ! number of observed quantities: 3 ((u,v), T, RH)
    mxrante          ,& ! number of (previous) real elements in buffer
    kk      , mv     ,& ! vertical loop indices
    ni2     , ni4    ,& ! loop indices (over time, upper/lower,...)
    mxi2    , mxi4   ,& ! upper limit of above loop indices
    mxlo             ,& ! number of data per information type and station
    noix             ,& ! variable index for the obs. increment record
    ilev             ,& ! vertical loop index for the obs. increment record
    nexceed          ,& ! number of reports in excess of array size
!   npsalin          ,& ! number of assimilated surface pressure increments
    izmplcode        ,& ! for MPI error code
!   iscount          ,& ! number of observations on local process
!   iscdim , ircdim  ,& ! dimensions of sending and receiving counters of obs.
!   irsdim           ,& ! number of real    elements sent to each processor
!   icc              ,& ! loop index
    nlevtot          ,& ! total number of observation increments levels in 'oiml
!   ncount (mxobtp+1),& ! counter for reports with obs increments (obs type)
    nzerr  , nstat   ,& ! indicator of status for (de-)allocation of fields
    ierr   , jerr       ! error stati for model_abort

  INTEGER (KIND=iintegers) ::  &
                        ! sum of reports from all nodes
                        !   (note: the other variables 'nmlall', 'nuaall',
                        !    'nsuall', 'npsall', 'nqcall', and 'nivall' are
                        !    defined globally within this module because they
                        !    are also used in 'gather_varia' for printout)
    noiall           ,& ! vertical profiles (reports)
    nysall              ! extrapolation of surface-level

! Local (automatic) arrays:
! -------------------------

! INTEGER (KIND=iintegers)   ::       &
!   ircount(num_compute), &! number of observations (or active obs. stations)
!   irrdim (num_compute), &! number of real elements received from each process
!   irdispl(num_compute)   ! displacement of incoming real data within 'rrvec'

  INTEGER (KIND=iintegers), ALLOCATABLE  ::       &
    ibufps   (:,:)    ,& ! for intermediate storage
    ibufsu   (:,:)    ,& ! for intermediate storage
    ibufua   (:,:)    ,& ! for intermediate storage
    ibufml   (:,:)    ,& ! for intermediate storage
    ibufqc   (:,:)       ! for intermediate storage

  REAL    (KIND=wp)       , ALLOCATABLE  ::       &
    rbufps   (:,:)    ,& ! for intermediate storage
    rbufsu   (:,:)    ,& ! for intermediate storage
    rbufua   (:,:)    ,& ! for intermediate storage
    rbufml   (:,:)    ,& ! for intermediate storage
    rbufqc   (:,:)    ,& ! for intermediate storage
    rbufoi   (:,:)    ,& ! for intermediate storage
    rbufps1    (:)       ! for intermediate storage

  CHARACTER (LEN=ilstidg) , ALLOCATABLE  ::       &
    ybufps   (:,:)    ,& ! for intermediate storage
    ybufsu   (:,:)    ,& ! for intermediate storage
    ybufua   (:,:)    ,& ! for intermediate storage
    ybufml   (:,:)       ! for intermediate storage

  CHARACTER (LEN=ilstidp) , ALLOCATABLE  ::       &
    ybufqc   (:,:)       ! for intermediate storage

  CHARACTER (LEN=75) yerrmsg , yerr
  CHARACTER (LEN=25) yroutine

!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin Subroutine gather_local_info
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Section 00: Synchronizing nodes
!-------------------------------------------------------------------------------

  yroutine  = 'gather_local_info'
  nzerr     = 0

  IF (ltime) THEN
    CALL get_timings (i_local_info, ntstep, dt, nzerr)
    IF (ltime_barrier) THEN
      CALL comm_barrier (icomm_cart, nzerr, yerrmsg)
      CALL get_timings (i_barrier_waiting_nud, ntstep, dt, nzerr)
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
!- Section 0: Gathering required info. for the spatial consistency checks
!-------------------------------------------------------------------------------
!- Section 0a: Write alerting messages for insufficient ODR size, if required
!-------------------------------------------------------------------------------

IF (nactio == 0) THEN

  IF (num_compute > 1) THEN

    CALL global_values ( nexcess, 5, 'SUM',imp_integers,icomm_cart, 0,yerr,ierr)
!   ==================
    IF (ierr /= 0)  CALL model_abort (my_cart_id, 11015, yerr, yroutine)
!                   ----------------
  ENDIF

  IF ((my_cart_id == 0) .AND. (MAXVAL(nexcess(1:5)) > 0)) THEN
    OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'            &
                               , POSITION='APPEND', IOSTAT=nstat)
    IF (nstat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
    IF (nstat /= 0) CALL model_abort (my_cart_id, 1409, yerr, yroutine)
    IF (nexcess(1) > 0) THEN
      WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," LOCAL MULTI-LEVEL OBS."  &
                     &," INCR. BEYOND maxmloi_tot ",I5)' )                     &
             ntstep, nexcess(1), maxmloi_tot
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxmlo OR maxtvo "  &
                     &,"BY AT LEAST",I6)' ) nexcess(1)
      nexcess (1) = nexcess(1) * 2
      WRITE( nucautn,'("   ==>  or INCREASE NAMELIST VARIABLE maxgpo BY AT L"  &
                     &,"EAST",I6)' ) nexcess(1)
    ENDIF
    IF (nexcess(2) > 0) THEN
      WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," LOCAL UPPER-AIR SINGLE-" &
                     &,"LEVEL OBS INCR. > maxuso ",I5)' )                      &
             ntstep, nexcess(2), maxuso
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxuso BY AT LEAST" &
                     &,I6)' ) nexcess(2)
    ENDIF
    IF (nexcess(3) > 0) THEN
      WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," LOCAL SURFACE-LEVEL OBS" &
                     &,". INCR. BEYOND maxsgo ",I5)' )                         &
             ntstep, nexcess(3), maxsgo
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxsgo BY AT LEAST" &
                     &,I6)' ) nexcess(3)
    ENDIF
    IF (nexcess(4) > 0) THEN
      WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," LOCAL SURFACE PRESSURE " &
                     &,"OBS. INCR. BEYOND maxpso ",I5)' )                      &
             ntstep, nexcess(4), maxpso
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxsgo BY AT LEAST" &
                     &,I6)' ) nexcess(4)
    ENDIF
    IF (nexcess(5) > 0) THEN
      WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," LOCAL IWV "              &
                     &,"OBS. INCR. BEYOND maxivq ",I5)' )                      &
             ntstep, nexcess(5), maxivq
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxgpo BY AT LEAST" &
                     &,I6)' ) nexcess(5)
    ENDIF
    CLOSE (nucautn)
  ENDIF

!-------------------------------------------------------------------------------
!- Section 0b: Gathering info for spatial consistency check of surface pressure 
!-------------------------------------------------------------------------------

! IF (npstot > 0) THEN

! allocate space for buffers
! --------------------------

    mxobs    = maxpso
    mxiobty  =  2 + 2 + 3
    mxrobty  = 11 + 4
    mxyobty  =  1
    ichlen   = ilstidg

    IF (nzerr == 0) ALLOCATE (ibufps (mxiobty, mxobs), STAT=nzerr)
    IF (nzerr == 0) ALLOCATE (rbufps (mxrobty, mxobs), STAT=nzerr)
    IF (nzerr == 0) ALLOCATE (ybufps (mxyobty, mxobs), STAT=nzerr)
    jerr = ABS( nzerr )
    CALL global_values ( jerr, 1,'MAX', imp_integers,icomm_cart, -1, yerr,ierr )
!   ==================
    IF ((jerr /= 0) .OR. (ierr /= 0)) THEN
!     yerrmsg = ' *** Allocation of sending buffers failed ***'
      WRITE( yerrmsg,'(" *** Allocation of xbufps failed ***",2I5)' ) jerr, ierr
      CALL model_abort (my_cart_id, 13072, yerrmsg, yroutine)
!     ================
    ENDIF

! write 'local' information to buffers
! ------------------------------------

    DO ista = 1 , npstot
      ibufps ( 1,ista) = iops(ista) + isubpos(my_cart_id,1) - nboundlines - 1
      ibufps ( 2,ista) = jops(ista) + isubpos(my_cart_id,2) - nboundlines - 1
      ibufps ( 3,ista) = 0
      ibufps ( 4,ista) = 0
      ibufps ( 5,ista) = iqclps (ista)
      ibufps ( 6,ista) = iqcfps (ista)
      ibufps ( 7,ista) = iqcnps (ista)
      rbufps ( 1,ista) = zoips  (ista,1)
      rbufps ( 2,ista) = zoips  (ista,2)
      rbufps ( 3,ista) = omykps (ista,1)
      rbufps ( 4,ista) = omykps (ista,2)
      rbufps ( 5,ista) = r1ifps (ista)
      rbufps ( 6,ista) = zmassps(ista)
      rbufps ( 7,ista) = ztdpps (ista)
      rbufps ( 8,ista) = qcimps (ista,1)
      rbufps ( 9,ista) = qcimps (ista,2)
      rbufps (10,ista) = qctmps (ista,1)
      rbufps (11,ista) = qctmps (ista,2)
      rbufps (12,ista) = qcqfps (ista,1)
      rbufps (13,ista) = qcqfps (ista,2)
      rbufps (14,ista) = zoips_b(ista,1)
      rbufps (15,ista) = zoips_b(ista,2)
      ybufps ( 1,ista) = ystidps(ista)
      IF (ltips(ista)) ibufps (3,ista) = 1
      IF (lmlps(ista)) ibufps (4,ista) = 1
    ENDDO

! gather information globally
! ---------------------------
 
    CALL gather_all ( mxiobty, mxrobty, mxyobty, mxobs , ichlen                &
                    , ibufps , rbufps , ybufps , npstot, npsall, yerrmsg, nzerr)
!   ===============

! read global information from buffers
! ------------------------------------

    DO ista = 1 , npstot
      iops_tot (ista) = ibufps(1,ista)
      jops_tot (ista) = ibufps(2,ista)
      iops (ista) = iops_tot(ista) - isubpos(my_cart_id,1) + nboundlines + 1
      jops (ista) = jops_tot(ista) - isubpos(my_cart_id,2) + nboundlines + 1

      ltips   (ista)   = (ibufps( 3,ista) == 1)
      lmlps   (ista)   = (ibufps( 4,ista) == 1)
      iqclps  (ista)   = ibufps ( 5,ista)
      iqcfps  (ista)   = ibufps ( 6,ista)
      iqcnps  (ista)   = ibufps ( 7,ista)
      zoips   (ista,1) = rbufps ( 1,ista)
      zoips   (ista,2) = rbufps ( 2,ista)
      omykps  (ista,1) = rbufps ( 3,ista)
      omykps  (ista,2) = rbufps ( 4,ista)
      r1ifps  (ista)   = rbufps ( 5,ista)
      zmassps (ista)   = rbufps ( 6,ista)
      ztdpps  (ista)   = rbufps ( 7,ista)
      qcimps  (ista,1) = rbufps ( 8,ista)
      qcimps  (ista,2) = rbufps ( 9,ista)
      qctmps  (ista,1) = rbufps (10,ista)
      qctmps  (ista,2) = rbufps (11,ista)
      qcqfps  (ista,1) = rbufps (12,ista)
      qcqfps  (ista,2) = rbufps (13,ista)
      zoips_b (ista,1) = rbufps (14,ista)
      zoips_b (ista,2) = rbufps (15,ista)
      ystidps (ista)   = ybufps (1,ista)
    ENDDO

! deallocate space for buffers
! ----------------------------

    DEALLOCATE (ibufps, STAT=nzerr)
    DEALLOCATE (rbufps, STAT=nzerr)
    DEALLOCATE (ybufps, STAT=nzerr)

! checking for array bounds
! -------------------------

    IF ((npsall > npstot) .AND. (lcroot)) THEN
      OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'          &
                                 , POSITION='APPEND', IOSTAT=nstat)
      IF (nstat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
      IF (nstat /= 0) CALL model_abort (my_cart_id, 1408, yerr, yroutine)
      WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," SURFACE PRESSURE OBS. "  &
                     &,"INCREMENTS, ARRAY SIZE ",I5)' ) ntstep, npsall, npstot
      nexceed  =  npsall - npstot
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxsgo BY AT LEAST" &
                     &,I6)' ) nexceed
      CLOSE (nucautn)
    ENDIF

! ENDIF

!-------------------------------------------------------------------------------
!- Section 0c: Gathering info for spatial consistency check of IWV
!-------------------------------------------------------------------------------

! IF (nivtot > 0) THEN

! allocate space for buffers
! --------------------------

    mxobs    = maxivq
    mxiobty  = 2 + 5
    mxrobty  = 6 + 2
    mxyobty  = 1
    ichlen   = ilstidg

    IF (nzerr == 0) ALLOCATE (ibufqc (mxiobty, mxobs), STAT=nzerr)
    IF (nzerr == 0) ALLOCATE (rbufqc (mxrobty, mxobs), STAT=nzerr)
    IF (nzerr == 0) ALLOCATE (ybufqc (mxyobty, mxobs), STAT=nzerr)
    jerr = ABS( nzerr )
    CALL global_values ( jerr, 1,'MAX', imp_integers,icomm_cart, -1, yerr,ierr )
!   ==================
    IF ((jerr /= 0) .OR. (ierr /= 0)) THEN
!     yerrmsg = ' *** Allocation of sending buffers failed ***'
      WRITE( yerrmsg,'(" *** Allocation of xbufiv failed ***",2I5)' ) jerr, ierr
      CALL model_abort (my_cart_id, 13072, yerrmsg, yroutine)
!     ================
    ENDIF

! write 'local' information to buffers
! ------------------------------------

    DO ista = 1 , nivtot
      ibufqc ( 1,ista) = ioiv(ista) + isubpos(my_cart_id,1) - nboundlines - 1
      ibufqc ( 2,ista) = joiv(ista) + isubpos(my_cart_id,2) - nboundlines - 1
      ibufqc ( 3,ista) = iqcliv (ista)
      ibufqc ( 4,ista) = iqcfiv (ista)
      ibufqc ( 5,ista) = kiobiv (ista)
      ibufqc ( 6,ista) = kioiiv (ista)
      ibufqc ( 7,ista) = ktypiv (ista)
      rbufqc ( 1,ista) = zoiciv (ista,1)
      rbufqc ( 2,ista) = zoiciv (ista,2)
      rbufqc ( 3,ista) = zqcfiv (ista,1)
      rbufqc ( 4,ista) = zqcfiv (ista,2)
      rbufqc ( 5,ista) = zmodiv (ista)
      rbufqc ( 6,ista) = zsativ (ista)
      rbufqc ( 7,ista) = zactiv (ista,1)
      rbufqc ( 8,ista) = zactiv (ista,2)
      ybufqc ( 1,ista) = ystidiv(ista)
    ENDDO

! gather information globally
! ---------------------------
 
    CALL gather_all ( mxiobty, mxrobty, mxyobty, mxobs , ichlen                &
                    , ibufqc , rbufqc , ybufqc , nivtot, nivall, yerrmsg, nzerr)
!   ===============

! read global information from buffers
! ------------------------------------

    DO ista = 1 , nivtot
      ioiv_tot (ista) = ibufqc(1,ista)
      joiv_tot (ista) = ibufqc(2,ista)
      ioiv (ista) = ioiv_tot(ista) - isubpos(my_cart_id,1) + nboundlines + 1
      joiv (ista) = joiv_tot(ista) - isubpos(my_cart_id,2) + nboundlines + 1
      iqcliv  (ista)   = ibufqc ( 3,ista)
      iqcfiv  (ista)   = ibufqc ( 4,ista)
      kiobiv  (ista)   = ibufqc ( 5,ista)
      kioiiv  (ista)   = ibufqc ( 6,ista)
      ktypiv  (ista)   = ibufqc ( 7,ista)
      zoiciv  (ista,1) = rbufqc ( 1,ista)
      zoiciv  (ista,2) = rbufqc ( 2,ista)
      zqcfiv  (ista,1) = rbufqc ( 3,ista)
      zqcfiv  (ista,2) = rbufqc ( 4,ista)
      zmodiv  (ista)   = rbufqc ( 5,ista)
      zsativ  (ista)   = rbufqc ( 6,ista)
      zactiv  (ista,1) = rbufqc ( 7,ista)
      zactiv  (ista,2) = rbufqc ( 8,ista)
      ystidiv (ista)   = ybufqc ( 1,ista)
    ENDDO

! deallocate space for buffers
! ----------------------------

    DEALLOCATE (ibufqc, STAT=nzerr)
    DEALLOCATE (rbufqc, STAT=nzerr)
    DEALLOCATE (ybufqc, STAT=nzerr)

! checking for array bounds
! -------------------------

    IF ((nivall > nivtot) .AND. (lcroot)) THEN
      OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'          &
                                 , POSITION='APPEND', IOSTAT=nstat)
      IF (nstat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
      IF (nstat /= 0) CALL model_abort (my_cart_id, 1408, yerr, yroutine)
      WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," IWV INCREMENTS FOR HUM"  &
                     &,"IDITY CHECK, ARRAY SIZE ",I5)' ) ntstep, nivall, nivtot
      nexceed  =  nivall - nivtot
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxmlo OR maxgpo "  &
                     &,"BY AT LEAST",I6)' ) nexceed
      CLOSE (nucautn)
    ENDIF

! ENDIF

ENDIF   !  if (nactio == 0)

!-------------------------------------------------------------------------------
!- Section 1: Gathering info for printing for control
!-------------------------------------------------------------------------------

IF (nactio == 1) THEN

! Gathering number of local upper-air obs not used due to ODR size
! ----------------------------------------------------------------

  CALL global_values ( nuaex, 1,'MAX', imp_integers, icomm_cart, 0, yerr, ierr)
! ==================
  IF (ierr /= 0)  CALL model_abort (my_cart_id, 13014, yerr, yroutine)
!                 ----------------

!-------------------------------------------------------------------------------
!- Section 1a: Gathering data rejected by the threshold quality control
!-------------------------------------------------------------------------------
 
! IF (ntotqc > 0) THEN

! allocate space for buffers
! --------------------------

    mxobs   = maxqcp
    mxiobty = 2
    mxrobty = 11
    mxyobty = 1
    ichlen  = ilstidp

    IF (nzerr == 0) ALLOCATE (ibufqc (mxiobty, mxobs), STAT=nzerr)
    IF (nzerr == 0) ALLOCATE (rbufqc (mxrobty, mxobs), STAT=nzerr)
    IF (nzerr == 0) ALLOCATE (ybufqc (mxyobty, mxobs), STAT=nzerr)
    jerr = ABS( nzerr )
    CALL global_values ( jerr, 1,'MAX', imp_integers,icomm_cart, -1, yerr,ierr )
!   ==================
    IF ((jerr /= 0) .OR. (ierr /= 0)) THEN
      WRITE( yerrmsg,'(" *** Allocation of xbufqc failed ***",2I5)' ) jerr, ierr
      CALL model_abort (my_cart_id, 13070, yerrmsg, yroutine)
!     ================
    ENDIF

! write 'local' information to buffers
! ------------------------------------

    DO   ista = 1 , ntotqc
      DO ivrs = 1 , 11
        rbufqc (ivrs,ista) = oyqc(ista,ivrs)
      ENDDO
      ibufqc (1,ista) = myqc(ista,1)
      ibufqc (2,ista) = myqc(ista,2)
      ybufqc (1,ista) = yyqc(ista  )
    ENDDO

! gather information globally
! ---------------------------
 
    CALL gather_all ( mxiobty, mxrobty, mxyobty, mxobs , ichlen                &
                    , ibufqc , rbufqc , ybufqc , ntotqc, nqcall, yerrmsg, nzerr)
!   ===============

! read global information from buffers
! ------------------------------------

    DO   ista = 1 , ntotqc
      DO ivrs = 1 , 11
        oyqc (ista,ivrs) = rbufqc (ivrs,ista)
      ENDDO
      myqc (ista,1) = ibufqc (1,ista)
      myqc (ista,2) = ibufqc (2,ista)
      yyqc (ista  ) = ybufqc (1,ista)
    ENDDO

! deallocate space for buffers
! ----------------------------

    DEALLOCATE (ibufqc, STAT=nzerr)
    DEALLOCATE (rbufqc, STAT=nzerr)
    DEALLOCATE (ybufqc, STAT=nzerr)

! checking for array bounds
! -------------------------

    IF ((nqcall > ntotqc) .AND. (lcroot)) THEN
      OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'          &
                                 , POSITION='APPEND', IOSTAT=nstat)
      IF (nstat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
      IF (nstat /= 0) CALL model_abort (my_cart_id, 1408, yerr, yroutine)
      WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," REPORTS REJECTED BY QU"  &
                     &,"ALITY CTRL., ARRAY SIZE ",I5)' ) ntstep, nqcall, ntotqc
      nexceed  = (nqcall - ntotqc) *20
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxsgo BY AT LEAST" &
                     &,I6)' ) nexceed
      nexceed  =  nexceed / maxoil + 1
      WRITE( nucautn,'("     OR INCREASE NAMELIST VARIABLE maxmlo BY AT LEAST" &
                     &,I6)' ) nexceed
      CLOSE (nucautn)
    ENDIF

! ENDIF

!-------------------------------------------------------------------------------
!- Section 1b: Gathering interpolated surface-level data
!-------------------------------------------------------------------------------
 
! IF (ntotys > 0) THEN
  IF (lfirst) THEN

! allocate space for buffers
! --------------------------

    mxobs   = maxysu
    mxiobty = 0
    mxrobty = 19
    mxyobty = 1
    ichlen  = ilstidp

    IF (nzerr == 0) ALLOCATE (ibufqc (MAX(mxiobty,i1), mxobs), STAT=nzerr)
    IF (nzerr == 0) ALLOCATE (rbufqc (mxrobty, mxobs), STAT=nzerr)
    IF (nzerr == 0) ALLOCATE (ybufqc (mxyobty, mxobs), STAT=nzerr)
    jerr = ABS( nzerr )
    CALL global_values ( jerr, 1,'MAX', imp_integers,icomm_cart, -1, yerr,ierr )
!   ==================
    IF ((jerr /= 0) .OR. (ierr /= 0)) THEN
      WRITE( yerrmsg,'(" *** Allocation of xbufqc failed ***",2I5)' ) jerr, ierr
      CALL model_abort (my_cart_id, 13071, yerrmsg, yroutine)
!     ================
    ENDIF

! write 'local' information to buffers
! ------------------------------------

    DO   ista = 1 , ntotys
      DO ivrs = 1 , 19
        rbufqc (ivrs,ista) = oysu(ista,ivrs)
      ENDDO
      ybufqc (1,ista) = yysu(ista)
    ENDDO

! gather information globally
! ---------------------------
 
    CALL gather_all ( mxiobty, mxrobty, mxyobty, mxobs , ichlen                &
                    , ibufqc , rbufqc , ybufqc , ntotys, nysall, yerrmsg, nzerr)
!   ===============

! read global information from buffers
! ------------------------------------

    DO   ista = 1 , ntotys
      DO ivrs = 1 , 19
        oysu (ista,ivrs) = rbufqc (ivrs,ista)
      ENDDO
      yysu (ista) = ybufqc (1,ista)
    ENDDO

! deallocate space for buffers
! ----------------------------

    DEALLOCATE (ibufqc, STAT=nzerr)
    DEALLOCATE (rbufqc, STAT=nzerr)
    DEALLOCATE (ybufqc, STAT=nzerr)

! ENDIF
  ENDIF

ENDIF   !  if (nactio == 1)

!-------------------------------------------------------------------------------

IF ((nactio == 2) .AND. (ltnudge)) THEN

!-------------------------------------------------------------------------------
!- Section 2: Gathering required info. on 'surface' pressure obs. increments:
!             modifications to local weights by spatial consistency check only
!-------------------------------------------------------------------------------

! gather local numbers of observation levels in 'oiml' from each node
! -------------------------------------------------------------------

  IF (npstot > 0) THEN

! allocate space for buffers
! --------------------------

    mxrobty = 2 * npstot

    IF (nzerr == 0) ALLOCATE (rbufps1 (mxrobty), STAT=nzerr)
    jerr = ABS( nzerr )
    CALL global_values ( jerr, 1,'MAX', imp_integers,icomm_cart, -1, yerr,ierr )
!   ==================
    IF ((jerr /= 0) .OR. (ierr /= 0)) THEN
      WRITE( yerrmsg,'(" *** Gathering of rbufps1 failed ***",2I5)' ) jerr, ierr
      CALL model_abort (my_cart_id, 13079, yerrmsg, yroutine)
!     ================
    ENDIF

! write ('local') information to buffers
! --------------------------------------

    DO ista = 1 , npstot
      rbufps1 (       ista) = omykps (ista,1)
      rbufps1 (npstot+ista) = omykps (ista,2)
    ENDDO

! gather information globally
! ---------------------------
 
    CALL global_values ( rbufps1, mxrobty, 'MAX', imp_reals, icomm_cart, -1    &
                       , yerrmsg, izmplcode )
!   ==================

    jerr = ABS( izmplcode )
    CALL global_values ( jerr, 1,'MAX', imp_integers,icomm_cart, -1, yerr, ierr)
!   ==================
    IF ((jerr /= 0) .OR. (ierr /= 0)) THEN
      WRITE( yerrmsg,'(" *** Gathering of rbufps1 failed ***",2I5)' ) jerr, ierr
      CALL model_abort (my_cart_id, 13081, yerrmsg, yroutine, izmplcode)
!     ================
    ENDIF

! read global information from buffers
! ------------------------------------

    DO ista = 1 , npstot
      omykps  (ista,1) = rbufps1 (       ista)
      omykps  (ista,2) = rbufps1 (npstot+ista)
    ENDDO

! deallocate space for buffers
! ----------------------------

    DEALLOCATE (rbufps1, STAT=nzerr)

  ENDIF

!-------------------------------------------------------------------------------
!- Section 3: Gathering required info. on other surface-level obs. increments
!-------------------------------------------------------------------------------
 
! IF (nsutot > 0) THEN

! determine globally uppermost level influenced by surface-level data
! -------------------------------------------------------------------

    CALL global_values ( ktopsu, 1, 'MIN', imp_integers, icomm_cart, -1        &
                       , yerrmsg, nzerr)
!   ==================

! allocate space for buffers
! --------------------------

    mxvrs = 3

    mxobs    = maxsgo
    mxiobty  = 4 + mxvrs + mxvrs
    mxrobty  = 6 + 8*mxvrs + 2*(ke - ktopsu + 1)
    IF (lnissu) mxrobty = mxrobty + ke - ktopsu + 1
    IF ((lnissu) .AND. (msprpsu >= 1)) mxrobty = mxrobty + mxispr
    mxyobty  = 1
    ichlen   = ilstidg

    IF (nzerr == 0) ALLOCATE (ibufsu (mxiobty, mxobs), STAT=nzerr)
    IF (nzerr == 0) ALLOCATE (rbufsu (mxrobty, mxobs), STAT=nzerr)
    IF (nzerr == 0) ALLOCATE (ybufsu (mxyobty, mxobs), STAT=nzerr)
    jerr = ABS( nzerr )
    CALL global_values ( jerr, 1,'MAX', imp_integers,icomm_cart, -1, yerr,ierr )
!   ==================
    IF ((jerr /= 0) .OR. (ierr /= 0)) THEN
      WRITE( yerrmsg,'(" *** Allocation of xbufsu failed ***",2I5)' ) jerr, ierr
      CALL model_abort (my_cart_id, 13073, yerrmsg, yroutine)
!     ================
    ENDIF

! write 'local' information to buffers
! ------------------------------------

    DO ista = 1 , nsutot
      ibufsu (1,ista) = iosu(ista) + isubpos(my_cart_id,1) - nboundlines - 1
      ibufsu (2,ista) = josu(ista) + isubpos(my_cart_id,2) - nboundlines - 1
      DO ivrs  = 1 , mxvrs
        ibufsu (ivrs+2      ,ista) = kviflsu(ista,  ivrs)
        ibufsu (ivrs+2+mxvrs,ista) = 0
        rbufsu (ivrs        ,ista) = xoisu  (ista,1,ivrs)
        rbufsu (ivrs+  mxvrs,ista) = xoisu  (ista,2,ivrs)
        rbufsu (ivrs+2*mxvrs,ista) = zwtsu  (ista,1,ivrs)
        rbufsu (ivrs+3*mxvrs,ista) = zwtsu  (ista,2,ivrs)
        rbufsu (ivrs+4*mxvrs,ista) = zvcutsu(ista,  ivrs)
        rbufsu (ivrs+5*mxvrs,ista) = zriflsu(ista,  ivrs)
        rbufsu (ivrs+6*mxvrs,ista) = zqualsu(ista,1,ivrs)
        rbufsu (ivrs+7*mxvrs,ista) = zqualsu(ista,2,ivrs)
        IF (ltisu(ista,ivrs)) ibufsu (ivrs+2+mxvrs,ista) = 1
      ENDDO
      rbufsu (8*mxvrs+1,ista) = xoisu  (ista,1,4)
      rbufsu (8*mxvrs+2,ista) = xoisu  (ista,2,4)
      rbufsu (8*mxvrs+3,ista) = zsposu (ista)
      rbufsu (8*mxvrs+4,ista) = zrtdgsu(ista)
      rbufsu (8*mxvrs+5,ista) = zpobsu (ista)
      rbufsu (8*mxvrs+6,ista) = zzobsu (ista)
      ibufsu (2*mxvrs+3,ista) = kobtysu(ista)
      ibufsu (2*mxvrs+4,ista) = kcdtysu(ista)
      mxrante = 8 *mxvrs + 6 - ktopsu + 1
      mv = ke - ktopsu + 1
      DO kk = ktopsu , ke
        rbufsu (mxrante+kk   ,ista) = zsprsu(ista,kk)
        rbufsu (mxrante+kk+mv,ista) = zpblsu(ista,kk)
        IF (lnissu) rbufsu (mxrante+kk+2*mv,ista) = znissu(ista,kk)
      ENDDO
      IF ((lnissu) .AND. (msprpsu >= 1)) THEN
        mxrante = mxrante + 3*mv + ktopsu - 1
        DO mv = 1 , mxispr
          rbufsu (mxrante+mv,ista) = znissq(ista,mv)
        ENDDO
      ENDIF
      ybufsu (1,ista) = ystidsu (ista)
    ENDDO

! gather information globally
! ---------------------------
 
    CALL gather_all ( mxiobty, mxrobty, mxyobty, mxobs , ichlen                &
                    , ibufsu , rbufsu , ybufsu , nsutot, nsuall, yerrmsg, nzerr)
!   ===============

! read global information from buffers
! ------------------------------------

    DO ista = 1 , nsutot
      iosu_tot (ista) = ibufsu(1,ista)
      josu_tot (ista) = ibufsu(2,ista)
      iosu (ista) = iosu_tot(ista) - isubpos(my_cart_id,1) + nboundlines + 1
      josu (ista) = josu_tot(ista) - isubpos(my_cart_id,2) + nboundlines + 1
      DO ivrs  = 1 , mxvrs
        kviflsu(ista,  ivrs) = ibufsu (ivrs+2      ,ista)
        xoisu  (ista,1,ivrs) = rbufsu (ivrs        ,ista)
        xoisu  (ista,2,ivrs) = rbufsu (ivrs+  mxvrs,ista)
        zwtsu  (ista,1,ivrs) = rbufsu (ivrs+2*mxvrs,ista)
        zwtsu  (ista,2,ivrs) = rbufsu (ivrs+3*mxvrs,ista)
        zvcutsu(ista,  ivrs) = rbufsu (ivrs+4*mxvrs,ista)
        zriflsu(ista,  ivrs) = rbufsu (ivrs+5*mxvrs,ista)
        zqualsu(ista,1,ivrs) = rbufsu (ivrs+6*mxvrs,ista)
        zqualsu(ista,2,ivrs) = rbufsu (ivrs+7*mxvrs,ista)
        ltisu  (ista,  ivrs) = (ibufsu(ivrs+2+mxvrs,ista) == 1)
      ENDDO
      xoisu  (ista,1,4) = rbufsu (8*mxvrs+1,ista)
      xoisu  (ista,2,4) = rbufsu (8*mxvrs+2,ista)
      zsposu (ista)     = rbufsu (8*mxvrs+3,ista)
      zrtdgsu(ista)     = rbufsu (8*mxvrs+4,ista)
      zpobsu (ista)     = rbufsu (8*mxvrs+5,ista)
      zzobsu (ista)     = rbufsu (8*mxvrs+6,ista)
      kobtysu(ista)     = ibufsu (2*mxvrs+3,ista)
      kcdtysu(ista)     = ibufsu (2*mxvrs+4,ista)
      mxrante = 8 *mxvrs + 6 - ktopsu + 1
      mv = ke - ktopsu + 1
      DO kk = ktopsu , ke
        zsprsu (ista,kk) = rbufsu (mxrante+kk   ,ista)
        zpblsu (ista,kk) = rbufsu (mxrante+kk+mv,ista)
        IF (lnissu) znissu (ista,kk) = rbufsu (mxrante+kk+2*mv,ista)
      ENDDO
      IF ((lnissu) .AND. (msprpsu >= 1)) THEN
        mxrante = mxrante + 3*mv + ktopsu - 1
        DO mv = 1 , mxispr
          znissq (ista,mv) = rbufsu (mxrante+mv,ista)
        ENDDO
      ENDIF
      ystidsu (ista) = ybufsu (1,ista)
    ENDDO

! deallocate space for buffers
! ----------------------------

    DEALLOCATE (ibufsu, STAT=nzerr)
    DEALLOCATE (rbufsu, STAT=nzerr)
    DEALLOCATE (ybufsu, STAT=nzerr)

! checking for array bounds
! -------------------------

    IF ((nsuall > nsutot) .AND. (lcroot)) THEN
      OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'          &
                                 , POSITION='APPEND', IOSTAT=nstat)
      IF (nstat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
      IF (nstat /= 0) CALL model_abort (my_cart_id, 1408, yerr, yroutine)
      WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," SURFACE-LEVEL OBS. INC"  &
                     &,"REMENTS, ARRAY SIZE ",I5)' ) ntstep, nsuall, nsutot
      nexceed  =  nsuall - nsutot
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxsgo BY AT LEAST" &
                     &,I6)' ) nexceed
      CLOSE (nucautn)
    ENDIF

! ENDIF


!-------------------------------------------------------------------------------
!- Section 4: Gathering required info. on upper-air single-level obs. increments
!-------------------------------------------------------------------------------

! IF (nuatot > 0) THEN

! allocate space for buffers
! --------------------------

    mxvrs = 3
    mxi2  = 2
    mxlo  = mxvrs * mxi2

!   mxobs    = maxsgo
    mxobs    = maxuso
    mxiobty  = 4 + mxlo + mxvrs
    mxrobty  = 8 + 6*mxlo + 2*ke
    IF (lnisua) mxrobty = mxrobty + ke
    IF ((lnisua) .AND. (msprpar >= 1)) mxrobty = mxrobty + mxispr
    mxyobty  = 1
    ichlen   = ilstidg

    IF (nzerr == 0) ALLOCATE (ibufua (mxiobty, mxobs), STAT=nzerr)
    IF (nzerr == 0) ALLOCATE (rbufua (mxrobty, mxobs), STAT=nzerr)
    IF (nzerr == 0) ALLOCATE (ybufua (mxyobty, mxobs), STAT=nzerr)
    jerr = ABS( nzerr )
    CALL global_values ( jerr, 1,'MAX', imp_integers,icomm_cart, -1, yerr,ierr )
!   ==================
    IF ((jerr /= 0) .OR. (ierr /= 0)) THEN
      WRITE( yerrmsg,'(" *** Allocation of xbufua failed ***",2I5)' ) jerr, ierr
      CALL model_abort (my_cart_id, 13074, yerrmsg, yroutine)
!     ================
    ENDIF

! write 'local' information to buffers
! ------------------------------------

    DO ista = 1 , nuatot
      ibufua (1,ista) = ioua(ista) + isubpos(my_cart_id,1) - nboundlines - 1
      ibufua (2,ista) = joua(ista) + isubpos(my_cart_id,2) - nboundlines - 1
      DO ivrs  = 1 , mxvrs
        DO ni2 = 1 , mxi2
          ibufua ((ivrs-1)*mxi2+ni2+     2,ista) = kviflua(ista,ni2,ivrs)
          rbufua ((ivrs-1)*mxi2+ni2       ,ista) = xoiua  (ista,ni2,ivrs)
          rbufua ((ivrs-1)*mxi2+ni2+  mxlo,ista) = zwtua  (ista,ni2,ivrs)
          rbufua ((ivrs-1)*mxi2+ni2+2*mxlo,ista) = fcorlua(ista,ni2,ivrs)
          rbufua ((ivrs-1)*mxi2+ni2+3*mxlo,ista) = zvcutua(ista,ni2,ivrs)
          rbufua ((ivrs-1)*mxi2+ni2+4*mxlo,ista) = zriflua(ista,ni2,ivrs)
          rbufua ((ivrs-1)*mxi2+ni2+5*mxlo,ista) = zqualua(ista,ni2,ivrs)
        ENDDO
      ENDDO
      rbufua (6*mxlo+1,ista) = xoiua  (ista,1,4)
      rbufua (6*mxlo+2,ista) = xoiua  (ista,2,4)
      rbufua (6*mxlo+3,ista) = zsprob (ista)
      rbufua (6*mxlo+4,ista) = zsprtp (ista)
      rbufua (6*mxlo+5,ista) = zrtdgua(ista)
      rbufua (6*mxlo+6,ista) = zpobua (ista)
      rbufua (6*mxlo+7,ista) = zzobua (ista)
      rbufua (6*mxlo+8,ista) = zvidua (ista)
      ibufua (  mxlo+3,ista) = kobtyua(ista)
      ibufua (  mxlo+4,ista) = kcdtyua(ista)
      DO ivrs  = 1 , mxvrs
        ibufua (ivrs+mxlo+4,ista) = 0
        IF (ltiua(ista,ivrs)) ibufua (ivrs+mxlo+4,ista) = 1
      ENDDO
      mxrante = 6 *mxlo + 8
      DO kk = 1 , ke
        rbufua (mxrante+kk   ,ista) = zsprua(ista,kk)
        rbufua (mxrante+kk+ke,ista) = zlopua(ista,kk)
        IF (lnisua) rbufua (mxrante+kk+2*ke,ista) = znisua(ista,kk)
      ENDDO
      IF ((lnisua) .AND. (msprpar >= 1)) THEN
        mxrante = mxrante + 3* ke
        DO mv = 1 , mxispr
          rbufua (mxrante+mv,ista) = znisuq(ista,mv)
        ENDDO
      ENDIF
      ybufua (1,ista) = ystidua (ista)
    ENDDO

! gather information globally
! ---------------------------
 
    CALL gather_all ( mxiobty, mxrobty, mxyobty, mxobs , ichlen                &
                    , ibufua , rbufua , ybufua , nuatot, nuaall, yerrmsg, nzerr)
!   ===============

! read global information from buffers
! ------------------------------------

    DO ista = 1 , nuatot
      ioua_tot (ista) = ibufua(1,ista)
      joua_tot (ista) = ibufua(2,ista)
      ioua (ista) = ioua_tot(ista) - isubpos(my_cart_id,1) + nboundlines + 1
      joua (ista) = joua_tot(ista) - isubpos(my_cart_id,2) + nboundlines + 1
      DO ivrs  = 1 , mxvrs
        DO ni2 = 1 , mxi2
          kviflua (ista,ni2,ivrs) = ibufua ((ivrs-1)*mxi2+ni2+     2,ista)
          xoiua   (ista,ni2,ivrs) = rbufua ((ivrs-1)*mxi2+ni2       ,ista)
          zwtua   (ista,ni2,ivrs) = rbufua ((ivrs-1)*mxi2+ni2+  mxlo,ista)
          fcorlua (ista,ni2,ivrs) = rbufua ((ivrs-1)*mxi2+ni2+2*mxlo,ista)
          zvcutua (ista,ni2,ivrs) = rbufua ((ivrs-1)*mxi2+ni2+3*mxlo,ista)
          zriflua (ista,ni2,ivrs) = rbufua ((ivrs-1)*mxi2+ni2+4*mxlo,ista)
          zqualua (ista,ni2,ivrs) = rbufua ((ivrs-1)*mxi2+ni2+5*mxlo,ista)
        ENDDO
      ENDDO
      xoiua  (ista,1,4) = rbufua (6*mxlo+1,ista)
      xoiua  (ista,2,4) = rbufua (6*mxlo+2,ista)
      zsprob (ista)     = rbufua (6*mxlo+3,ista)
      zsprtp (ista)     = rbufua (6*mxlo+4,ista)
      zrtdgua(ista)     = rbufua (6*mxlo+5,ista)
      zpobua (ista)     = rbufua (6*mxlo+6,ista)
      zzobua (ista)     = rbufua (6*mxlo+7,ista)
      zvidua (ista)     = rbufua (6*mxlo+8,ista)
      kobtyua(ista)     = ibufua (  mxlo+3,ista)
      kcdtyua(ista)     = ibufua (  mxlo+4,ista)
      DO ivrs  = 1 , mxvrs
        ltiua (ista,ivrs) = (ibufua(ivrs+mxlo+4,ista) == 1)
      ENDDO
      mxrante = 6 *mxlo + 8
      DO kk = 1 , ke
        zsprua (ista,kk) = rbufua (mxrante+kk   ,ista)
        zlopua (ista,kk) = rbufua (mxrante+kk+ke,ista)
        IF (lnisua) znisua (ista,kk) = rbufua (mxrante+kk+2*ke,ista)
      ENDDO
      IF ((lnisua) .AND. (msprpar >= 1)) THEN
        mxrante = mxrante + 3* ke
        DO mv = 1 , mxispr
          znisuq (ista,mv) = rbufua (mxrante+mv,ista)
        ENDDO
      ENDIF
      ystidua (ista) = ybufua (1,ista)
    ENDDO

! deallocate space for buffers
! ----------------------------

    DEALLOCATE (ibufua, STAT=nzerr)
    DEALLOCATE (rbufua, STAT=nzerr)
    DEALLOCATE (ybufua, STAT=nzerr)

! checking for array bounds
! -------------------------

    IF ((nuaall > nuatot) .AND. (lcroot)) THEN
      OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'          &
                                 , POSITION='APPEND', IOSTAT=nstat)
      IF (nstat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
      IF (nstat /= 0) CALL model_abort (my_cart_id, 1408, yerr, yroutine)
      WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," UPPER-AIR SINGLE-LEVEL"  &
                     &," OBS. INCR., ARRAY SIZE ",I5)' ) ntstep, nuaall, nuatot
      nexceed  =  nuaall - nuatot + nuaex
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxuso BY AT LEAST" &
                     &,I6)' ) nexceed
      CLOSE (nucautn)
    ENDIF

! ENDIF


!-------------------------------------------------------------------------------
!- Section 5: Gathering required info. on multi-level obs. increments
!-------------------------------------------------------------------------------
!- Section 5a: Gathering general local information on multi-level reports
!-------------------------------------------------------------------------------

! IF (nmltot > 0) THEN

! allocate space for buffers
! --------------------------

    mxvrs = 3
    mxi4  = 4
    mxlo  = mxvrs * mxi4

    mxobs    = maxmloi_tot
    mxiobty  =   mxlo + 4*mxvrs + 9 + mxvrs
    mxrobty  = 7*mxlo + 2*mxvrs + 1 + 3*ke
    IF (lnisua) mxrobty = mxrobty + ke
    IF ((lnisua) .AND. (msprpar >= 1)) mxrobty = mxrobty + mxispr
    mxyobty  = 1
    ichlen   = ilstidg

    IF (nzerr == 0) ALLOCATE (ibufml (mxiobty, mxobs), STAT=nzerr)
    IF (nzerr == 0) ALLOCATE (rbufml (mxrobty, mxobs), STAT=nzerr)
    IF (nzerr == 0) ALLOCATE (ybufml (mxyobty, mxobs), STAT=nzerr)
    jerr = ABS( nzerr )
    CALL global_values ( jerr, 1,'MAX', imp_integers,icomm_cart, -1, yerr,ierr )
!   ==================
    IF ((jerr /= 0) .OR. (ierr /= 0)) THEN
      WRITE( yerrmsg,'(" *** Allocation of xbufml failed ***",2I5)' ) jerr, ierr
      CALL model_abort (my_cart_id, 13075, yerrmsg, yroutine)
!     ================
    ENDIF

! write 'local' information to buffers
! ------------------------------------

    DO ista = 1 , nmltot
      ibufml (1,ista) = ioml(ista) + isubpos(my_cart_id,1) - nboundlines - 1
      ibufml (2,ista) = joml(ista) + isubpos(my_cart_id,2) - nboundlines - 1
      DO ivrs  = 1 , mxvrs
        DO ni4 = 1 , mxi4
          ibufml ((ivrs-1)*mxi4+ni4+     2,ista) = kviflml(ista,ni4,ivrs)
          rbufml ((ivrs-1)*mxi4+ni4       ,ista) = zspobml(ista,ni4,ivrs)
          rbufml ((ivrs-1)*mxi4+ni4+  mxlo,ista) = zwtml  (ista,ni4,ivrs)
          rbufml ((ivrs-1)*mxi4+ni4+2*mxlo,ista) = fcorlml(ista,ni4,ivrs)
          rbufml ((ivrs-1)*mxi4+ni4+3*mxlo,ista) = zvcutml(ista,ni4,ivrs)
          rbufml ((ivrs-1)*mxi4+ni4+4*mxlo,ista) = zrtdgml(ista,ni4,ivrs)
          rbufml ((ivrs-1)*mxi4+ni4+5*mxlo,ista) = zpobml (ista,ni4,ivrs)
          rbufml ((ivrs-1)*mxi4+ni4+6*mxlo,ista) = zemkpml(ista,ni4,ivrs)
        ENDDO
      ENDDO
      DO ivrs  = 1 , mxvrs
        rbufml (7*mxlo+          ivrs,ista) = zriflml(ista,1,ivrs)
        rbufml (7*mxlo+    mxvrs+ivrs,ista) = zriflml(ista,2,ivrs)
        ibufml (  mxlo+2+        ivrs,ista) = mbotlv (ista,1,ivrs)
        ibufml (  mxlo+2+  mxvrs+ivrs,ista) = mbotlv (ista,2,ivrs)
        ibufml (  mxlo+2+2*mxvrs+ivrs,ista) = mtoplv (ista,1,ivrs)
        ibufml (  mxlo+2+3*mxvrs+ivrs,ista) = mtoplv (ista,2,ivrs)
        ibufml (  mxlo+2+4*mxvrs+ivrs,ista) = 0
        IF (ltiml(ista,ivrs))  ibufml (mxlo+2+4*mxvrs+ivrs,ista) = 1
      ENDDO
      rbufml (7*mxlo+2*mxvrs+1,ista) = zstpml (ista)
      ibufml (  mxlo+5*mxvrs+3,ista) = mszlev (ista,1)
      ibufml (  mxlo+5*mxvrs+4,ista) = mszlev (ista,2)
      ibufml (  mxlo+5*mxvrs+5,ista) = ksprml (ista)
      ibufml (  mxlo+5*mxvrs+6,ista) = kkoiml (ista,1)
      ibufml (  mxlo+5*mxvrs+7,ista) = kkoiml (ista,2)
      ibufml (  mxlo+5*mxvrs+8,ista) = kobtyml(ista)
      ibufml (  mxlo+5*mxvrs+9,ista) = kcdtyml(ista)
      mxrante = 7 *mxlo + 2 *mxvrs + 1
      DO kk = 1 , ke
        rbufml (mxrante+kk     ,ista) = zsprml (ista,kk)
        rbufml (mxrante+kk+  ke,ista) = zlopml (ista,kk)
        rbufml (mxrante+kk+2*ke,ista) = zrtgkml(ista,kk)
        IF (lnisua) rbufml (mxrante+kk+3*ke,ista) = znisml(ista,kk)
      ENDDO
      IF ((lnisua) .AND. (msprpar >= 1)) THEN
        mxrante = mxrante + 4* ke
        DO mv = 1 , mxispr
          rbufml (mxrante+mv,ista) = znismq(ista,mv)
        ENDDO
      ENDIF
      ybufml (1,ista) = ystidml (ista)
    ENDDO

! gather information globally
! ---------------------------
 
    CALL gather_all ( mxiobty, mxrobty, mxyobty, mxobs , ichlen                &
                    , ibufml , rbufml , ybufml , nmltot, nmlall, yerrmsg, nzerr)
!   ===============

! read global information from buffers
! ------------------------------------

    DO ista = 1 , nmltot
      ioml_tot (ista) = ibufml(1,ista)
      joml_tot (ista) = ibufml(2,ista)
      ioml (ista) = ioml_tot(ista) - isubpos(my_cart_id,1) + nboundlines + 1
      joml (ista) = joml_tot(ista) - isubpos(my_cart_id,2) + nboundlines + 1
      DO ivrs  = 1 , mxvrs
        DO ni4 = 1 , mxi4
          kviflml (ista,ni4,ivrs) = ibufml ((ivrs-1)*mxi4+ni4+     2,ista)
          zspobml (ista,ni4,ivrs) = rbufml ((ivrs-1)*mxi4+ni4       ,ista)
          zwtml   (ista,ni4,ivrs) = rbufml ((ivrs-1)*mxi4+ni4+  mxlo,ista)
          fcorlml (ista,ni4,ivrs) = rbufml ((ivrs-1)*mxi4+ni4+2*mxlo,ista)
          zvcutml (ista,ni4,ivrs) = rbufml ((ivrs-1)*mxi4+ni4+3*mxlo,ista)
          zrtdgml (ista,ni4,ivrs) = rbufml ((ivrs-1)*mxi4+ni4+4*mxlo,ista)
          zpobml  (ista,ni4,ivrs) = rbufml ((ivrs-1)*mxi4+ni4+4*mxlo,ista)
          zemkpml (ista,ni4,ivrs) = rbufml ((ivrs-1)*mxi4+ni4+6*mxlo,ista)
        ENDDO
      ENDDO
      DO ivrs  = 1 , mxvrs
        zriflml(ista,1,ivrs) = rbufml (7*mxlo+          ivrs,ista)
        zriflml(ista,2,ivrs) = rbufml (7*mxlo+    mxvrs+ivrs,ista)
        mbotlv (ista,1,ivrs) = ibufml (  mxlo+2+        ivrs,ista)
        mbotlv (ista,2,ivrs) = ibufml (  mxlo+2+  mxvrs+ivrs,ista)
        mtoplv (ista,1,ivrs) = ibufml (  mxlo+2+2*mxvrs+ivrs,ista)
        mtoplv (ista,2,ivrs) = ibufml (  mxlo+2+3*mxvrs+ivrs,ista)
        ltiml  (ista  ,ivrs) = (ibufml(  mxlo+2+4*mxvrs+ivrs,ista) == 1)
      ENDDO
      zstpml (ista)   = rbufml (7*mxlo+2*mxvrs+1,ista)
      mszlev (ista,1) = ibufml (  mxlo+5*mxvrs+3,ista)
      mszlev (ista,2) = ibufml (  mxlo+5*mxvrs+4,ista)
      ksprml (ista)   = ibufml (  mxlo+5*mxvrs+5,ista)
      kkoiml (ista,1) = ibufml (  mxlo+5*mxvrs+6,ista)
      kkoiml (ista,2) = ibufml (  mxlo+5*mxvrs+7,ista)
      kobtyml(ista)   = ibufml (  mxlo+5*mxvrs+8,ista)
      kcdtyml(ista)   = ibufml (  mxlo+5*mxvrs+9,ista)
      mxrante = 7 *mxlo + 2 *mxvrs + 1
      DO kk = 1 , ke
        zsprml (ista,kk) = rbufml (mxrante+kk     ,ista)
        zlopml (ista,kk) = rbufml (mxrante+kk+  ke,ista)
        zrtgkml(ista,kk) = rbufml (mxrante+kk+2*ke,ista)
        IF (lnisua) znisml (ista,kk) = rbufml (mxrante+kk+3*ke,ista)
      ENDDO
      IF ((lnisua) .AND. (msprpar >= 1)) THEN
        mxrante = mxrante + 4* ke
        DO mv = 1 , mxispr
          znismq (ista,mv) = rbufml (mxrante+mv,ista)
        ENDDO
      ENDIF
      ystidml (ista) = ybufml (1,ista)
    ENDDO

! deallocate space for buffers
! ----------------------------

    DEALLOCATE (ibufml, STAT=nzerr)
    DEALLOCATE (rbufml, STAT=nzerr)
    DEALLOCATE (ybufml, STAT=nzerr)

! checking for array bounds
! -------------------------

    IF ((nmlall > nmltot) .AND. (lcroot)) THEN
      OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'          &
                                 , POSITION='APPEND', IOSTAT=nstat)
      IF (nstat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
      IF (nstat /= 0) CALL model_abort (my_cart_id, 1408, yerr, yroutine)
      WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," MULTI-LEVEL STATIONS O"  &
                     &,"F OBS. INCR., ARRAY SIZE ",I5)' ) ntstep, nmlall, nmltot
      nexceed  =  nmlall - nmltot
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxmlo, maxgpo OR"  &
                     &," maxtvo BY AT LEAST",I6)' ) nexceed
      CLOSE (nucautn)
    ENDIF

! ENDIF

!-------------------------------------------------------------------------------
!- Section 5b: Gathering the observation increment record 'oiml' (OIR)
!-------------------------------------------------------------------------------
 
! IF (nmloit > 0) THEN

! get global number of observation levels in 'oiml'
! -------------------------------------------------

    nlevtot = nmloit

    CALL global_values ( nlevtot, 1,'SUM',imp_integers,icomm_cart, -1,yerr,ierr)
!   ==================
    IF (ierr /= 0) THEN
      WRITE( yerrmsg,'(" *** Gathering on nlevtot failed ***",2I5)' ) jerr, ierr
      CALL model_abort (my_cart_id, 13075, yerrmsg, yroutine)
!     ================
    ENDIF
    IF (nlevtot > 0) THEN

! adjust size of OIR 'oiml' by re-allocation and allocate space for buffer
! ------------------------------------------------------------------------

      IF (nzerr == 0) ALLOCATE (rbufoi (maxnoi , MAX( nmloit,i1 )) , STAT=nzerr)
      jerr = ABS( nzerr )
      CALL global_values ( jerr, 1,'MAX',imp_integers,icomm_cart, -1,yerr,ierr )
!     ==================
      IF ((jerr /= 0) .OR. (ierr /= 0)) THEN
        WRITE( yerrmsg,'(" *** Allocation of xbufml failed ***",2I5)') jerr,ierr
        CALL model_abort (my_cart_id, 13075, yerrmsg, yroutine)
!       ================
      ENDIF

      DO ilev   = 1 , nmloit
        DO noix = 1 , maxnoi
          rbufoi (noix,ilev) = oiml(noix,ilev)
        ENDDO
      ENDDO
                    DEALLOCATE ( oiml                       , STAT=nzerr)
      IF (nzerr == 0) ALLOCATE ( oiml   (maxnoi , nlevtot)  , STAT=nzerr)
      IF (nzerr == 0) ALLOCATE ( ibufml (1      , nlevtot)  , STAT=nzerr)
      IF (nzerr == 0) ALLOCATE ( ybufml (1      , nlevtot)  , STAT=nzerr)

      jerr = ABS( nzerr )
      CALL global_values ( jerr, 1,'MAX',imp_integers,icomm_cart, -1,yerr,ierr )
!     ==================
      IF ((jerr /= 0) .OR. (ierr /= 0)) THEN
        WRITE( yerrmsg,'(" *** Allocation of oiml failed ***",2I5)' ) jerr, ierr
        CALL model_abort (my_cart_id, 13077, yerrmsg, yroutine)
!       ================
      ENDIF

      DO ilev   = 1 , nmloit
        DO noix = 1 , maxnoi
          oiml (noix,ilev) = rbufoi(noix,ilev)
        ENDDO
      ENDDO
                    DEALLOCATE ( rbufoi                     , STAT=nzerr)

! gather information globally
! ---------------------------

      CALL gather_all ( 0     , maxnoi, 0     , nlevtot, ilstidg               &
                      , ibufml, oiml  , ybufml, nmloit , noiall, yerrmsg, nzerr)
!     ===============

    ENDIF

! ENDIF

!-------------------------------------------------------------------------------
!- Section 5c: Making the level index pointer of 5a consistent with 'oiml' of 5b
!-------------------------------------------------------------------------------

    mxrante = 0
    DO ista = 1 , nmltot
      DO itim = 1 , 2
        IF (kkoiml(ista,itim) > 0) THEN
          kkoiml (ista,itim) = mxrante + 1
          mxrante            = mxrante + mszlev(ista,itim)
        ENDIF
      ENDDO
    ENDDO

!-------------------------------------------------------------------------------
!- Section 6: Checking for array bounds
!-------------------------------------------------------------------------------

  IF (lcroot) THEN
!   IF (noiall > nmloit)                                                       &
!     PRINT       '(''CAUTION !!!!! t='',I5,'':'',I5,'' VERTICAL PROFILES OF'' &
!                 &,'' OBS. INCR., ARRAY SIZE '',I5)' , ntstep, noiall, nmloit
    IF (nmlall > nmltot)                                                       &
      PRINT       '(''CAUTION !!!!! t='',I5,'':'',I5,'' MULTI-LEVEL STATIONS'' &
                  &,''OF OBS. INCR., ARRAY SIZE '',I5)' , ntstep, nmlall, nmltot
    IF (nuaall > nuatot)                                                       &
      PRINT       '(''CAUTION !!!!! t='',I5,'':'',I5,'' UPPER-AIR SINGLE-LEV'' &
                  &,''EL OBS. INCR., ARRAY SIZE '',I5)' , ntstep, nuaall, nuatot
    IF (nsuall > nsutot)                                                       &
      PRINT       '(''CAUTION !!!!! t='',I5,'':'',I5,'' SURFACE-LEVEL OBS. I'' &
                  &,''NCREMENTS, ARRAY SIZE '',I5)' , ntstep,  nsuall, nsutot
    IF (npsall > npstot)                                                       &
      PRINT       '(''CAUTION !!!!! t='',I5,'':'',I5,'' SURFACE PRESSURE OBS'' &
                  &,''. INCREMENTS, ARRAY SIZE '',I5)' , ntstep, npsall, npstot
    IF (nivall > nivtot)                                                       &
      PRINT       '(''CAUTION !!!!! t='',I5,'':'',I5,'' IWV INCREMENTS FOR H'' &
                  &,''UMIDITY CHECK, ARRAY SIZE '',I5)' , ntstep, nivall, nivtot
    IF (nqcall > ntotqc)                                                       &
      PRINT       '(''CAUTION !!!!! t='',I5,'':'',I5,'' REPORTS REJECTED BY '' &
                  &,''QUALITY CTRL., ARRAY SIZE '',I5)' , ntstep, nqcall, ntotqc

!   npsalin = 0
!   DO ista = 1 , npstot
!     IF (MAX( omykps(ista,1),omykps(ista,2) ) > epsy)  npsalin = npsalin + 1
!   ENDDO
!   IF (npsalin == 0) npsall = 0
!   ncount  = 0
!   DO ista = 1 , nmltot
!     ! this 'IF' condition is required here to exclude GPS obs, which are set
!     ! passive in the threshold QC but still used in the spatial consistency
!     ! check (SCC), or which are set passive in the SCC
!     IF (MAXVAL( zwtml(ista,:,:) ) > epsy) THEN
!       ncount (kobtyml(ista))  =  ncount(kobtyml(ista)) + 1
!       ncount (mxobtp+1)       =  ncount(mxobtp+1)      + 1
!     ENDIF
!   ENDDO
!   PRINT       '(I4,": # STATIONS WITH :",I4," multi-level:"                  &
!               &,I4," TEMP,",I4," PILOT,",I4," AMDAR,",I4," GPS")' ,          &
!          ntstep, ncount(mxobtp+1), ncount(ntemp) , ncount(npilot)            &
!                                  , ncount(nairep), ncount(ngps)
!   PRINT       '(4X,"   OBS INCREMENTS :",I4," upper-air sing-lv,"            &
!               &,I5," surface,",I5," surf. pressure")' ,                      &
!          nuaall, nsuall, npsalin
  ENDIF

  IF (lwonl) THEN
!   IF (noiall > nmloit)                                                       &
!     WRITE( nupr,'(''CAUTION !!!!! t='',I5,'':'',I5,'' VERTICAL PROFILES OF'' &
!                 &,'' OBS. INCR., ARRAY SIZE '',I5)' ) ntstep, noiall, nmloit
    IF (nmlall > nmltot)                                                       &
      WRITE( nupr,'(''CAUTION !!!!! t='',I5,'':'',I5,'' MULTI-LEVEL STATIONS'' &
                  &,''OF OBS. INCR., ARRAY SIZE '',I5)' ) ntstep, nmlall, nmltot
    IF (nuaall > nuatot)                                                       &
      WRITE( nupr,'(''CAUTION !!!!! t='',I5,'':'',I5,'' UPPER-AIR SINGLE-LEV'' &
                  &,''EL OBS. INCR., ARRAY SIZE '',I5)' ) ntstep, nuaall, nuatot
    IF (nsuall > nsutot)                                                       &
      WRITE( nupr,'(''CAUTION !!!!! t='',I5,'':'',I5,'' SURFACE-LEVEL OBS. I'' &
                  &,''NCREMENTS, ARRAY SIZE '',I5)' ) ntstep, nsuall, nsutot
    IF (npsall > npstot)                                                       &
      WRITE( nupr,'(''CAUTION !!!!! t='',I5,'':'',I5,'' SURFACE PRESSURE OBS'' &
                  &,''. INCREMENTS, ARRAY SIZE '',I5)' ) ntstep, npsall, npstot
    IF (nivall > nivtot)                                                       &
      WRITE( nupr,'(''CAUTION !!!!! t='',I5,'':'',I5,'' IWV INCREMENTS FOR H'' &
                  &,''UMIDITY CHECK, ARRAY SIZE '',I5)' ) ntstep, nivall, nivtot
    IF (nqcall > ntotqc)                                                       &
      WRITE( nupr,'(''CAUTION !!!!! t='',I5,'':'',I5,'' REPORTS REJECTED BY '' &
                  &,''QUALITY CTRL., ARRAY SIZE '',I5)' ) ntstep, nqcall, ntotqc

! flush YUPRINT file
    IF (ldump_ascii) THEN
      CLOSE (nupr)
      OPEN  (nupr,FILE=yuprint,FORM='FORMATTED',STATUS='OLD'                   &
                              ,POSITION='APPEND',IOSTAT=nstat)
      IF (nstat /= 0) yerrmsg = 'OPENING OF FILE yuprint FAILED'
      IF (nstat /= 0) CALL model_abort (my_cart_id, 1009, yerrmsg, yroutine)
    ENDIF
 
  ENDIF


!-------------------------------------------------------------------------------

ENDIF   !  if (nactio == 2) .and. (ltnudge)

!-------------------------------------------------------------------------------
!- Section 7: Measuring time for communication
!-------------------------------------------------------------------------------

  IF (ltime) CALL get_timings (i_communications_nud, ntstep, dt, nzerr)

!-------------------------------------------------------------------------------
!- End of the Subroutine gather_local_info
!-------------------------------------------------------------------------------

END SUBROUTINE gather_local_info


!-------------------------------------------------------------------------------
!+ Module procedure to compute auxiliary quantities for spreading
!-------------------------------------------------------------------------------

SUBROUTINE gather_spread_aux ( nactio , k )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure computes various auxiliary quantities (mostly fields)
!   for the spreading of observation increments, for the sake of computational
!   efficiency.
!
! Method:
!   Straightforward computations.
!
! Written by        :  Christoph Schraff, DWD  (original version: 14.10.97)
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

! Subroutine arguments:
! --------------------

  INTEGER (KIND=iintegers), INTENT (IN)         ::       &
    nactio        ! determines, which quantities are to be computed (when)

  INTEGER (KIND=iintegers), INTENT (IN)  , OPTIONAL    ::    &
    k             ! index of current vertical model level
                  !   (required for nactio == 2)

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    i      , j       ,& ! loop indices in horizontal direction
    ie_t2  , je_t2   ,& ! loop indices for horiz. differences on total domain
    mv     , kk      ,& ! loop indices
    km1              ,& ! index of model level used for vertical gradients
    joffset          ,& ! latitudinal offset of subdomain in grid point units
    ityp             ,& ! index of set of obs system in weighted increment array
    nzerr               ! status of memory (de-)allocation

  REAL    (KIND=wp   )     ::  &
    zrerdkm          ,& ! earth radius in [km]
    zrlats           ,& ! latitude in rotated model coord.
    zcrlats          ,& ! cos( zrlats)
    zsrlats          ,& ! sin( zrlats)
    zalpha           ,& ! half the angle subtended by the arc as given by the
                        ! longitudinal distance on the tangent cone projection
                        ! from the earth's center
    salpha , calpha  ,& ! sin( zalpha )   ,   cos( zalpha )
    zalphac          ,& ! scaling factor for computing 'zalpha'
    zxxdfac          ,& ! scaling factor for computing longit. distance on proj.
    xxd              ,& ! longitudinal (zonal) distance on tangent cone project.
    zhml             ,& ! height z
    zemhh            ,& ! exp( -z / scale )   , z = height
    zeklop0          ,& ! scaling factor for computing potential temperature
!   zfablmin         ,& ! min. degree at which grid pts are assigned to the ABL
!   zhabl               ! absolute height of the ABL (atmosph. bound. layer) top
    c6r                 ! 1 / 6

! LOGICAL                  ::  &
!   lwabl               ! nudging weights reduced in ABL (atmosph. bound. layer)


! LOGICAL                  , SAVE ::  &
!   lwablk              ! model level not entirely above ABL

  REAL    (KIND = wp)      , ALLOCATABLE , SAVE :: &
    zlopb    (:,:) ,& ! log( pressure )   at level k+1
    zlopa    (:,:) ,& ! log( pressure )   at level k-1
    zeklpa   (:,:) ,& ! exp( R/cp * log(p) )   at level k-1
    zth      (:,:) ,& ! potential temperature
    zthb     (:,:) ,& ! potential temperature   at level k+1
    ztha     (:,:)    ! potential temperature   at level k-1


! Local (automatic) arrays: None
! -------------------------
!
!------------ End of header ----------------------------------------------------


 
!-------------------------------------------------------------------------------
! Begin Subroutine gather_spread_aux
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1: Geometrical fields used to computed horizontal distances on the
!             tangent cone projection, used only for the spatial consistency
!             check of surface pressure observations
!-------------------------------------------------------------------------------

  IF (nactio == 0) THEN

    je_t2 = 2*je_tot - 1
    ie_t2 = 2*ie_tot - 1
    IF (ALLOCATED( pyyd )) THEN
      PRINT '("CAUTION in src_gather_info: pyyd is already allocated "         &
            &,"at time / PE",I6,I5)', ntstep, my_cart_id
      DEALLOCATE ( pyyd , pxxd2 , pxsalpa , STAT=nzerr )
    ENDIF
    ALLOCATE ( pyyd    (je_t2             )    , STAT = nzerr )
    ALLOCATE ( pxxd2   (ie_t2, jstart:jend)    , STAT = nzerr )
    ALLOCATE ( pxsalpa (ie_t2, jstart:jend)    , STAT = nzerr )

    zrerdkm = r_earth / 1000.0_wp
    c6r = c1 / 6.0_wp

    DO j = 1 , je_t2
      pyyd (j) = zrerdkm * (j-je_tot)*dlat *degrad
    ENDDO
    joffset = NINT( (startlat - startlat_tot) /dlat )
    DO j  = jstart , jend
      zrlats   = startlat_tot + (joffset+j-1)*dlat
      zsrlats  = SIN( zrlats * degrad )
      IF (ABS( zsrlats ) > epsy) THEN
        zalphac = c05 * zsrlats * dlon * degrad
        zxxdfac = c2 * zrerdkm * crlat(j,1) / zsrlats
        DO i = 1 , ie_t2
          zalpha        = (i-ie_tot) * zalphac
!   approxi. for sin(zalpha)
          salpha        = zalpha - zalpha *zalpha *zalpha *c6r
          xxd           = salpha * zxxdfac
          pxxd2   (i,j) = xxd    * xxd
          pxsalpa (i,j) = xxd    * salpha * c2
        ENDDO
      ELSE
        zxxdfac = zrerdkm * dlon * degrad
        DO i = 1 , ie_t2
          xxd           = (i-ie_tot) * zxxdfac
          pxxd2   (i,j) = xxd * xxd
          pxsalpa (i,j) = c0
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  
!-------------------------------------------------------------------------------
!  Section 1: Quantities to be computed once at every timestep after exchanging
!             the local information on the observations between the nodes
!    !!!!     Note: in sections 1b and 2, all grid point fields are computed, 
!    !!!!           which are required subsequently in the spreading procedures
!    !!!!           except for:  - hhl: used in several procedures at k == ke !
!                                - several fields in ps_spreading (k == ke !)
!-------------------------------------------------------------------------------

  IF (nactio == 1) THEN

! allocation of arrays used for the spreading
! -------------------------------------------

    je_t2 = 2*je_tot - 1
    ie_t2 = 2*ie_tot - 1

    jrs_tot = jstart - NINT( (startlat - startlat_tot) /dlat )
    jre_tot = jrs_tot + je_tot - (jstart-1) - (je-jend)

    idimcut = (iend - istart + 1) * (jend - jstart + 1)

    IF (ALLOCATED( ztwips )) THEN
      PRINT '("CAUTION in src_gather_info: ztwips is already allocated "       &
            &,"at time / PE",I6,I5)', ntstep, my_cart_id
      DEALLOCATE ( ztwips , ztpgeo , faclbc , omy , om2 , zwi , STAT=nzerr )
    ENDIF
    ALLOCATE ( ztwips (istart:iend,jstart:jend,ke) , STAT = nzerr )
    ALLOCATE ( ztpgeo (istart:iend,jstart:jend,ke) , STAT = nzerr )
    ALLOCATE ( faclbc (ie         ,je         ,3 ) , STAT = nzerr )
    ztwips (:,:,:)  = c0
    ztpgeo (:,:,:)  = c0
    IF (mpsgcor >= 1) THEN
      ALLOCATE ( omy    (istart:iend,jstart:jend,5,nwtyp)  , STAT = nzerr )
      ALLOCATE ( om2    (istart:iend,jstart:jend,4,nwtyp)  , STAT = nzerr )
      ALLOCATE ( zwi    (istart:iend,jstart:jend,6,nwtyp)  , STAT = nzerr )
    ELSE
      ALLOCATE ( omy    (istart:iend,jstart:jend,4,nwtyp)  , STAT = nzerr )
      ALLOCATE ( om2    (istart:iend,jstart:jend,3,nwtyp)  , STAT = nzerr )
      ALLOCATE ( zwi    (istart:iend,jstart:jend,4,nwtyp)  , STAT = nzerr )
    ENDIF

    IF (ALLOCATED( yyd )) THEN
      PRINT '("CAUTION in src_gather_info: yyd is already allocated "          &
            &,"at time / PE",I6,I5)', ntstep, my_cart_id
      DEALLOCATE ( yyd    , yyd2   , xxd2            , STAT=nzerr )
      DEALLOCATE ( c2alpa , scalpa , xcalpa , xsalpa , STAT=nzerr )
      DEALLOCATE ( zlop   , zeklop , zthvg  , zpk    , STAT=nzerr )
      DEALLOCATE ( zspr   , zdds   , zablpo , zablfc , STAT=nzerr )
      DEALLOCATE ( zcoruu , zcoruv , zcorvu , zcorvv , STAT=nzerr )
      DEALLOCATE ( lcutof , icutof , jcutof , icutmp , jcutmp , STAT=nzerr )
      DEALLOCATE ( gppkmi , zsdnis , zsdnid          , STAT=nzerr )
      DEALLOCATE ( zlopb  , zlopa  , zeklpa , zth, zthb, ztha , STAT=nzerr )
    ENDIF
    ALLOCATE ( yyd    (je_t2                  )    , STAT = nzerr )
    ALLOCATE ( yyd2   (je_t2                  )    , STAT = nzerr )
    ALLOCATE ( xxd2   (ie_t2      ,jstart:jend)    , STAT = nzerr )
    ALLOCATE ( c2alpa (ie_t2      ,jstart:jend)    , STAT = nzerr )
    ALLOCATE ( scalpa (ie_t2      ,jstart:jend)    , STAT = nzerr )
    ALLOCATE ( xcalpa (ie_t2      ,jstart:jend)    , STAT = nzerr )
    ALLOCATE ( xsalpa (ie_t2      ,jstart:jend)    , STAT = nzerr )

    ALLOCATE ( zlop   (istart:iend,jstart:jend)    , STAT = nzerr )
    ALLOCATE ( zeklop (istart:iend,jstart:jend)    , STAT = nzerr )
    ALLOCATE ( zthvg  (istart:iend,jstart:jend)    , STAT = nzerr )
    ALLOCATE ( zpk    (istart:iend,jstart:jend)    , STAT = nzerr )
!   ALLOCATE ( zqdmod (istart:iend,jstart:jend)    , STAT = nzerr )
!   ALLOCATE ( zewmod (istart:iend,jstart:jend)    , STAT = nzerr )
!   ALLOCATE ( zewsat (istart:iend,jstart:jend)    , STAT = nzerr )
    ALLOCATE ( zspr   (istart:iend,jstart:jend,3)  , STAT = nzerr )
    ALLOCATE ( zdds   (istart:iend,jstart:jend,3)  , STAT = nzerr )
    ALLOCATE ( zcoruu (istart:iend,jstart:jend)    , STAT = nzerr )
    ALLOCATE ( zcoruv (istart:iend,jstart:jend)    , STAT = nzerr )
    ALLOCATE ( zcorvu (istart:iend,jstart:jend)    , STAT = nzerr )
    ALLOCATE ( zcorvv (istart:iend,jstart:jend)    , STAT = nzerr )
    ALLOCATE ( zablpo (istart:iend,jstart:jend)    , STAT = nzerr )
    ALLOCATE ( zablfc (istart:iend,jstart:jend,3,2), STAT = nzerr )

    ALLOCATE ( lcutof (istart:iend,jstart:jend,2)  , STAT = nzerr )
    ALLOCATE ( icutof (idimcut                ,2)  , STAT = nzerr )
    ALLOCATE ( jcutof (idimcut                ,2)  , STAT = nzerr )
    ALLOCATE ( icutmp (idimcut                )    , STAT = nzerr )
    ALLOCATE ( jcutmp (idimcut                )    , STAT = nzerr )

    ALLOCATE ( gppkmi (            jrs_tot:jre_tot), STAT = nzerr )
    ALLOCATE ( zsdnis (            mxispr     ,2)  , STAT = nzerr )
    ALLOCATE ( zsdnid (            mxispr     ,2)  , STAT = nzerr )

    ALLOCATE ( zlopb  (istart:iend,jstart:jend)    , STAT = nzerr )
    ALLOCATE ( zlopa  (istart:iend,jstart:jend)    , STAT = nzerr )
    ALLOCATE ( zeklpa (istart:iend,jstart:jend)    , STAT = nzerr )
    ALLOCATE ( zth    (istart:iend,jstart:jend)    , STAT = nzerr )
    ALLOCATE ( zthb   (istart:iend,jstart:jend)    , STAT = nzerr )
    ALLOCATE ( ztha   (istart:iend,jstart:jend)    , STAT = nzerr )

! computation of fields used for the spreading
! --------------------------------------------

! for non-isotropic 1-dimensional horizontal correlations

    DO mv = 1 , mxispr
      zsdnis (mv,1) = -sezispr * LOG( xezispr + (mv-1) *dezispr )
      zsdnis (mv,2) =                 xthispr + (mv-1) *dthispr
    ENDDO
    DO mv = 1 , mxispr-1
      zsdnid (mv,1) = c1 / (zsdnis(mv,1) - zsdnis(mv+1,1))
      zsdnid (mv,2) = c1 / (zsdnis(mv,2) - zsdnis(mv+1,2))
    ENDDO

! conversion factors for horiz. distances from 'km' to rotated grid point units

    gppkmj = 360.0_wp / (40042.5_wp * dlat)
    joffset = NINT( (startlat - startlat_tot) /dlat )
    DO j = jrs_tot , jre_tot
      zcrlats     = COS( (startlat_tot + (joffset+j-1) *dlat) * degrad )
      IF (zcrlats > epsy) THEN
        gppkmi (j) = 360.0_wp / (40042.5_wp * zcrlats * dlon)
      ELSE
        gppkmi (j) = je_tot / dlon
      ENDIF
    ENDDO

! Auxiliary fields containing all computations with trigonometric functions
! for horizontal distances and correlations

    zrerdkm = r_earth / 1000.0_wp

    DO j = 1 , je_t2
      yyd   (j)   = zrerdkm * (j-je_tot)*dlat *degrad
      yyd2  (j)   = yyd(j) **2
    ENDDO

    joffset = NINT( (startlat - startlat_tot) /dlat )
    DO j  = jstart , jend
! use 'startlat_tot' instead of 'startlat' to ensure reproducibility
      zrlats   = startlat_tot + (joffset+j-1)*dlat
      zsrlats  = SIN( zrlats * degrad )
      IF (ABS( zsrlats ) > epsy) THEN
        zalphac = c05 * zsrlats * dlon * degrad
        zxxdfac = c2 * zrerdkm * crlat(j,1) / zsrlats
        DO i = 1 , ie_t2
          zalpha       = (i-ie_tot) * zalphac
          salpha       = SIN( zalpha )
          calpha       = COS( zalpha )
          xxd          = salpha * zxxdfac
          xxd2   (i,j) = xxd    **2
          c2alpa (i,j) = calpha **2
          scalpa (i,j) = salpha * calpha
          xcalpa (i,j) = xxd    * calpha
          xsalpa (i,j) = xxd    * salpha * c2
        ENDDO
      ELSE
        zxxdfac = zrerdkm * dlon * degrad
        DO i = 1 , ie_t2
          xxd          = (i-ie_tot) * zxxdfac
          xxd2   (i,j) = xxd **2
          c2alpa (i,j) = c1
          scalpa (i,j) = c0
          xcalpa (i,j) = xxd
          xsalpa (i,j) = c0
        ENDDO
      ENDIF
    ENDDO

! reduction factors for nudging weights (near the lateral boundaries)

    DO   j = jstartpar , jendpar
      DO i = istartpar , iendpar
        faclbc(i,j,1) = MAX( c1 - rmy(i,j,1) /zmaxrmy(1) , c0 )
        faclbc(i,j,2) = MAX( c1 - rmy(i,j,2) /zmaxrmy(2) , c0 )
        faclbc(i,j,3) = MAX( c1 - rmy(i,j,3) /zmaxrmy(3) , c0 )
      ENDDO
    ENDDO

  ENDIF
 
!-------------------------------------------------------------------------------
!  Section 2: Quantities to be computed once for every model level prior to the
!             spreading of observation increments to that level
!    !!!!     (Note: a child-node only needs to gather:
!                    - yyd, xxd and the other geometrical fields
!                    - all fields required to compute section 4
!-------------------------------------------------------------------------------

  IF (nactio == 2) THEN

! allocation of 'zpai' and 'zroi' (after subr. 'ps_temperatur_corr'), and
! setting zpai = zroi = 0 in outermost 'nboundlines' grid lines of total domain

    IF (k == ke) THEN
      IF (ALLOCATED( zpai )) THEN
        PRINT '("CAUTION in src_gather_info: zpai is already allocated "       &
              &,"at time / PE",I6,I5)', ntstep, my_cart_id
        DEALLOCATE ( zpai , zroi , STAT=nzerr )
      ENDIF
      ALLOCATE ( zpai   (ie,je,ke) , STAT = nzerr )
      ALLOCATE ( zroi   (ie,je,ke) , STAT = nzerr )
      IF (luvgcor) THEN
        DO    kk = 1 , ke
!         DO   j = jstartpar , jendpar
!           DO i = istartpar , iendpar
          DO   j = 1 , je
            DO i = 1 , ie
              zpai (i,j,kk) = c0
              zroi (i,j,kk) = c0
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDIF

! compute quantities which are used for the spreading of observation increments

    zeklop0 = EXP( rdocp * LOG( p0r ) )

    IF (k == ke) THEN
      DO   j = jstart, jend
        DO i = istart, iend
          zlopa  (i,j) = LOG( p0(i,j,k) + pp(i,j,k,nnew) )
          zeklpa (i,j) = EXP( rdocp *zlopa(i,j) )
          ztha   (i,j) = t(i,j,k,nnew) * zeklop0 / zeklpa(i,j)
          zlop   (i,j) = zlopa(i,j)
          zth    (i,j) = ztha (i,j)
!         zrhke  (i,j) =   fq2pv( qd(i,j,k,nnew) , zpk(i,j) )                  &
!                        / fpvsw( t (i,j,k,nnew) )
        ENDDO
      ENDDO
    ENDIF
    km1 = MAX( k-1 , 1 )
    DO   j = jstart, jend
      DO i = istart, iend
        zlopb  (i,j) = zlop  (i,j)
        zthb   (i,j) = zth   (i,j)
        zlop   (i,j) = zlopa (i,j)
        zth    (i,j) = ztha  (i,j)
        zeklop (i,j) = zeklpa(i,j)
        zlopa  (i,j) = LOG( p0(i,j,km1) + pp(i,j,km1,nnew) )
        zeklpa (i,j) = EXP( rdocp *zlopa(i,j) )
        ztha   (i,j) = t(i,j,km1,nnew) * zeklop0 / zeklpa(i,j)
        zthvg  (i,j) = (ztha(i,j) - zthb(i,j)) / (zlopb(i,j) - zlopa(i,j))
        zpk    (i,j) = p0(i,j,k) + pp(i,j,k,nnew)
        zhml         = c05 * (hhl(i,j,k) + hhl(i,j,k+1))
        zemhh        = EXP( - zhml / sezispr )
!       zqdmod (i,j) = qd(i,j,k,nnew)
!       zewmod (i,j) = fq2pv( qd(i,j,k,nnew) , zpk(i,j) )
!       zewsat (i,j) = fpvsw( t (i,j,k,nnew) )
        zspr   (i,j,1) = zhml
        zspr   (i,j,2) = zth (i,j)
        zspr   (i,j,3) = zemhh
        zablpo (i,j)   = c1
        zablfc (i,j,1,1) = wablua(1) - zablpo(i,j) *(wablua(1) - c1)
        zablfc (i,j,2,1) = wablua(3) - zablpo(i,j) *(wablua(3) - c1)
        zablfc (i,j,3,1) = wablua(4) - zablpo(i,j) *(wablua(4) - c1)
        zablfc (i,j,1,2) = c1        - zablpo(i,j) *(c1 - wablsu(1))
        zablfc (i,j,2,2) = c1        - zablpo(i,j) *(c1 - wablsu(3))
        zablfc (i,j,3,2) = c1        - zablpo(i,j) *(c1 - wablsu(4))
      ENDDO
    ENDDO

! compute the degree at which the grid points are assigned to the ABL
! (needed for reduction of weights within / above the ABL)

!   lwabl = ( wablua(1) *wablua(3) *wablua(4)                                  &
!            *wablsu(1)            *wablsu(4) < c1-epsy)
!   IF (k == ke) lwablk = .TRUE.
!   IF ((lwabl) .AND. (k == ke)) CALL ablheit ( zdabl )
!   IF ((lwabl) .AND. (lwablk)) THEN
!     zfablmin = c1
!     DO   j = jstart, jend
!       DO i = istart, iend
!         zhabl        = hhl(i,j,ke+1) + zdabl(i,j)
!!        zablpo (i,j) = MAX( c0 , MIN( c1 ,  (hhl(i,j,k) - zhabl)             &
!!                                          / (hhl(i,j,k) - hhl(i,j,k+1)) ) )
!! no nudging as soon as a fraction of the model layer is within the ABL
!         zablpo (i,j) = MAX( c0 , SIGN( c1 , hhl(i,j,k+1) - zhabl ) )
!         zfablmin     = MIN( zfablmin , zablpo(i,j) )
!         zablfc (i,j,1,1) = wablua(1) - zablpo(i,j) *(wablua(1) - c1)
!         zablfc (i,j,2,1) = wablua(3) - zablpo(i,j) *(wablua(3) - c1)
!         zablfc (i,j,3,1) = wablua(4) - zablpo(i,j) *(wablua(4) - c1)
!         zablfc (i,j,1,2) = c1        - zablpo(i,j) *(c1 - wablsu(1))
!         zablfc (i,j,2,2) = c1        - zablpo(i,j) *(c1 - wablsu(3))
!         zablfc (i,j,3,2) = c1        - zablpo(i,j) *(c1 - wablsu(4))
!       ENDDO
!     ENDDO
!     IF (zfablmin >= c1-epsy) lwablk = .FALSE.
!     IF (      (lwonl) .AND. (k > ke/2) .AND. (ntstep >= 484)                 &
!         .AND. (ntstep <= 484+INT( tconbox/dt-epsy ))) THEN
!       zhabl = hhl(ionl,jonl,ke+1) + zdabl(ionl,jonl)
!       WRITE( nupr,'(''zablpo='',F4.2,'' , zhabl='',F6.0,'' , hhl at k+1,k''  &
!                   &,2F7.0)' )                                                &
!              zablpo(ionl,jonl), zhabl, hhl(ionl,jonl,k+1), hhl(ionl,jonl,k)
!     ENDIF
!   ENDIF

! initialize the quantities which are updated during the spreading of obs. incr.

    DO kk = 1 , 4
      km1 = MIN( 3, kk )
      DO ityp = 1 , nwtyp
        DO   j = jstart, jend
          DO i = istart, iend
            zwi (i,j,kk ,ityp) = c0
            omy (i,j,kk ,ityp) = c0
            om2 (i,j,km1,ityp) = c0
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    IF (mpsgcor >= 1) THEN
      DO ityp = 1 , nwtyp
        DO   j = jstart, jend
          DO i = istart, iend
            zwi (i,j,5  ,ityp) = c0
            zwi (i,j,6  ,ityp) = c0
            omy (i,j,5  ,ityp) = c0
            om2 (i,j,4  ,ityp) = c0
          ENDDO
        ENDDO
      ENDDO
    ENDIF

  ENDIF


!-------------------------------------------------------------------------------
!  Section 3: De-allocation of arrays used for the spreading and nudging
!-------------------------------------------------------------------------------

  IF ((nactio == 3) .AND. (ltnudge)) THEN

! de-allocation of fields used for the spreading

    DEALLOCATE ( yyd    , STAT = nzerr )
    DEALLOCATE ( yyd2   , STAT = nzerr )
    DEALLOCATE ( xxd2   , STAT = nzerr )
    DEALLOCATE ( c2alpa , STAT = nzerr )
    DEALLOCATE ( scalpa , STAT = nzerr )
    DEALLOCATE ( xcalpa , STAT = nzerr )
    DEALLOCATE ( xsalpa , STAT = nzerr )

    DEALLOCATE ( zlop   , STAT = nzerr )
    DEALLOCATE ( zeklop , STAT = nzerr )
    DEALLOCATE ( zthvg  , STAT = nzerr )
    DEALLOCATE ( zpk    , STAT = nzerr )
!   DEALLOCATE ( zqdmod , STAT = nzerr )
!   DEALLOCATE ( zewmod , STAT = nzerr )
!   DEALLOCATE ( zewsat , STAT = nzerr )
    DEALLOCATE ( zspr   , STAT = nzerr )
    DEALLOCATE ( zdds   , STAT = nzerr )
    DEALLOCATE ( zcoruu , STAT = nzerr )
    DEALLOCATE ( zcoruv , STAT = nzerr )
    DEALLOCATE ( zcorvu , STAT = nzerr )
    DEALLOCATE ( zcorvv , STAT = nzerr )
    DEALLOCATE ( zablpo , STAT = nzerr )
    DEALLOCATE ( zablfc , STAT = nzerr )

    DEALLOCATE ( lcutof , STAT = nzerr )
    DEALLOCATE ( icutof , STAT = nzerr )
    DEALLOCATE ( jcutof , STAT = nzerr )
    DEALLOCATE ( icutmp , STAT = nzerr )
    DEALLOCATE ( jcutmp , STAT = nzerr )

    DEALLOCATE ( gppkmi , STAT = nzerr )
    DEALLOCATE ( zsdnis , STAT = nzerr )
    DEALLOCATE ( zsdnid , STAT = nzerr )

    DEALLOCATE ( zlopb  , STAT = nzerr )
    DEALLOCATE ( zlopa  , STAT = nzerr )
    DEALLOCATE ( zeklpa , STAT = nzerr )
    DEALLOCATE ( zth    , STAT = nzerr )
    DEALLOCATE ( zthb   , STAT = nzerr )
    DEALLOCATE ( ztha   , STAT = nzerr )

! de-allocation of fields used for the nudging

    DEALLOCATE ( zpai   , STAT = nzerr )
    DEALLOCATE ( zroi   , STAT = nzerr )

    DEALLOCATE ( ztwips , STAT = nzerr )
    DEALLOCATE ( ztpgeo , STAT = nzerr )
    DEALLOCATE ( omy    , STAT = nzerr )
    DEALLOCATE ( om2    , STAT = nzerr )
    DEALLOCATE ( zwi    , STAT = nzerr )
    DEALLOCATE ( faclbc , STAT = nzerr )

  ENDIF

! de-allocation of arrays containing the local info. used for the spreading

  IF (nactio == 3) THEN

    DEALLOCATE ( oiml    , STAT = nzerr )
    DEALLOCATE ( zwtml   , STAT = nzerr )
    DEALLOCATE ( zspobml , STAT = nzerr )
    DEALLOCATE ( fcorlml , STAT = nzerr )
    DEALLOCATE ( zvcutml , STAT = nzerr )
    DEALLOCATE ( zrtdgml , STAT = nzerr )
    DEALLOCATE ( zpobml  , STAT = nzerr )
    DEALLOCATE ( zemkpml , STAT = nzerr )
    DEALLOCATE ( kviflml , STAT = nzerr )
    DEALLOCATE ( zriflml , STAT = nzerr )
    DEALLOCATE ( mbotlv  , STAT = nzerr )
    DEALLOCATE ( mtoplv  , STAT = nzerr )
    DEALLOCATE ( ltiml   , STAT = nzerr )
    DEALLOCATE ( mszlev  , STAT = nzerr )
    DEALLOCATE ( kkoiml  , STAT = nzerr )
    DEALLOCATE ( ksprml  , STAT = nzerr )
    DEALLOCATE ( kobtyml , STAT = nzerr )
    DEALLOCATE ( kcdtyml , STAT = nzerr )
    DEALLOCATE ( zstpml  , STAT = nzerr )
    DEALLOCATE ( ystidml , STAT = nzerr )
    DEALLOCATE ( ioml    , STAT = nzerr )
    DEALLOCATE ( joml    , STAT = nzerr )
    DEALLOCATE ( ioml_tot, STAT = nzerr )
    DEALLOCATE ( joml_tot, STAT = nzerr )
    DEALLOCATE ( zsprml  , STAT = nzerr )
    DEALLOCATE ( zlopml  , STAT = nzerr )
    DEALLOCATE ( zrtgkml , STAT = nzerr )
    IF (lnisua) DEALLOCATE ( znisml  , STAT = nzerr )
    IF ((lnisua) .AND. (msprpar >= 1)) DEALLOCATE ( znismq  , STAT = nzerr )

    DEALLOCATE ( xoiua   , STAT = nzerr )
    DEALLOCATE ( zwtua   , STAT = nzerr )
    DEALLOCATE ( fcorlua , STAT = nzerr )
    DEALLOCATE ( zvcutua , STAT = nzerr )
    DEALLOCATE ( zriflua , STAT = nzerr )
    DEALLOCATE ( zqualua , STAT = nzerr )
    DEALLOCATE ( kviflua , STAT = nzerr )
    DEALLOCATE ( ltiua   , STAT = nzerr )
    DEALLOCATE ( zsprob  , STAT = nzerr )
    DEALLOCATE ( zsprtp  , STAT = nzerr )
    DEALLOCATE ( zrtdgua , STAT = nzerr )
    DEALLOCATE ( zpobua  , STAT = nzerr )
    DEALLOCATE ( zzobua  , STAT = nzerr )
    DEALLOCATE ( zvidua  , STAT = nzerr )
    DEALLOCATE ( kobtyua , STAT = nzerr )
    DEALLOCATE ( kcdtyua , STAT = nzerr )
    DEALLOCATE ( ystidua , STAT = nzerr )
    DEALLOCATE ( ioua    , STAT = nzerr )
    DEALLOCATE ( joua    , STAT = nzerr )
    DEALLOCATE ( ioua_tot, STAT = nzerr )
    DEALLOCATE ( joua_tot, STAT = nzerr )
    DEALLOCATE ( zsprua  , STAT = nzerr )
    DEALLOCATE ( zlopua  , STAT = nzerr )
    IF (lnisua) DEALLOCATE ( znisua  , STAT = nzerr )
    IF ((lnisua) .AND. (msprpar >= 1)) DEALLOCATE ( znisuq  , STAT = nzerr )

    DEALLOCATE ( xoisu   , STAT = nzerr )
    DEALLOCATE ( zwtsu   , STAT = nzerr )
    DEALLOCATE ( zvcutsu , STAT = nzerr )
    DEALLOCATE ( zriflsu , STAT = nzerr )
    DEALLOCATE ( zqualsu , STAT = nzerr )
    DEALLOCATE ( kviflsu , STAT = nzerr )
    DEALLOCATE ( kobtysu , STAT = nzerr )
    DEALLOCATE ( kcdtysu , STAT = nzerr )
    DEALLOCATE ( ltisu   , STAT = nzerr )
    DEALLOCATE ( zsposu  , STAT = nzerr )
    DEALLOCATE ( zrtdgsu , STAT = nzerr )
    DEALLOCATE ( zpobsu  , STAT = nzerr )
    DEALLOCATE ( zzobsu  , STAT = nzerr )
    DEALLOCATE ( ystidsu , STAT = nzerr )
    DEALLOCATE ( iosu    , STAT = nzerr )
    DEALLOCATE ( josu    , STAT = nzerr )
    DEALLOCATE ( iosu_tot, STAT = nzerr )
    DEALLOCATE ( josu_tot, STAT = nzerr )
    DEALLOCATE ( zsprsu  , STAT = nzerr )
    DEALLOCATE ( zpblsu  , STAT = nzerr )
    IF (lnissu) DEALLOCATE ( znissu  , STAT = nzerr )
    IF ((lnissu) .AND. (msprpsu >= 1)) DEALLOCATE ( znissq  , STAT = nzerr )

    DEALLOCATE ( zoips   , STAT = nzerr )
    DEALLOCATE ( omykps  , STAT = nzerr )
    DEALLOCATE ( r1ifps  , STAT = nzerr )
    DEALLOCATE ( zmassps , STAT = nzerr )
    DEALLOCATE ( ztdpps  , STAT = nzerr )
    DEALLOCATE ( qcimps  , STAT = nzerr )
    DEALLOCATE ( qctmps  , STAT = nzerr )
    DEALLOCATE ( qcqfps  , STAT = nzerr )
    DEALLOCATE ( zoips_b , STAT = nzerr )
    DEALLOCATE ( ltips   , STAT = nzerr )
    DEALLOCATE ( lmlps   , STAT = nzerr )
    DEALLOCATE ( ystidps , STAT = nzerr )
    DEALLOCATE ( iops    , STAT = nzerr )
    DEALLOCATE ( jops    , STAT = nzerr )
    DEALLOCATE ( iops_tot, STAT = nzerr )
    DEALLOCATE ( jops_tot, STAT = nzerr )
    DEALLOCATE ( iqclps  , STAT = nzerr )
    DEALLOCATE ( iqcfps  , STAT = nzerr )
    DEALLOCATE ( iqcnps  , STAT = nzerr )

    DEALLOCATE ( isrtpqc , STAT = nzerr )
    DEALLOCATE ( isrtvqc , STAT = nzerr )
    DEALLOCATE ( isortps , STAT = nzerr )
    DEALLOCATE ( isortsu , STAT = nzerr )
    DEALLOCATE ( isortua , STAT = nzerr )
    DEALLOCATE ( isortml , STAT = nzerr )

  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure gather_spread_aux
!-------------------------------------------------------------------------------

END SUBROUTINE gather_spread_aux


!-------------------------------------------------------------------------------
!+ Module procedure to compute auxiliary fields for nudging
!-------------------------------------------------------------------------------

SUBROUTINE gather_nudge_aux

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure computes auxiliary fields for the nudging, namely
!   the vertical extent of convective instability near any grid point used for
!   balancing the humidity when nudging temperature.
!
! Method:
!   Straightforward computations. Exchange of boundaries between PE's.
!
! Written by        :  Christoph Schraff, DWD  (original version: 14.08.00)
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

! Subroutine arguments: None
! --------------------

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    i      , j       ,& ! loop indices in horizontal direction
    kk, kzdims(24)   ,& ! vertical index
    nseekln          ,& ! # lines in each direction used to seek convective g.p.
    nbdlines         ,& ! number of boundary lines of extended domain
    naddlin          ,& ! shift of inner domain within the extended domain
    nadd             ,& ! length extension of local domain
    ieadd  , jeadd   ,& ! index range of total extended local domain
    jstartex,jendex  ,& ! zonal  index range of extended shifted inner domain
    iad    , jad     ,& ! loop indices over shifted inner domain
    istad  , iendad  ,& ! merid. index range of shifted inner domain
    jstad  , jendad  ,& ! zonal  index range of shifted inner domain
    ncon             ,& ! number of convective grid points
    ncircle          ,& ! radius of area of influence of a convective grid pt.
    noffset          ,& ! indix offset
    ic               ,& ! counter index
    ix               ,& ! indices used for control output
    nzerr               ! status of memory (de-)allocation

  INTEGER (KIND=iintegers)               ::       &
    mconv      (32)     ! counter of (non-)convective grid points

  REAL    (KIND=wp   )     , ALLOCATABLE ::       &
    zconbgp   (:,:)  ,& ! height of base \ of convectively instable region
    zcontgp   (:,:)  ,& ! height of top  / at 1 grid point
    zconbmi (:,:,:)  ,& ! height of base \ of maximum extent of convectively
    zcontmi (:,:,:)     ! height of top  / instable region at 3 grid points 

  INTEGER (KIND=iintegers) , ALLOCATABLE ::       &
    mconbas     (:)  ,& ! layer index of base of convectively instable region
    mcontop     (:)  ,& ! layer index of top  of convectively instable region
    iix         (:)  ,& ! merid. index
    jix         (:)  ,& ! zonal  index
    nads        (:)     ! array defining (shape of) area of influence

! For error handling
! ------------------
  INTEGER (KIND=iintegers) ::  izerror

  CHARACTER (LEN=80)       ::  yzerrmsg


! Local (automatic) arrays: None
! -------------------------
!
!------------ End of header ----------------------------------------------------


 
!-------------------------------------------------------------------------------
! Begin Subroutine gather_nudge_aux
!-------------------------------------------------------------------------------

  kzdims(:) = 0_iintegers

!-------------------------------------------------------------------------------
!  Section 1: Determination of vertical extent of convection at grid points
!             with convective precipitation (if evaporation is neglected)
!-------------------------------------------------------------------------------

! Determine grid pts. with convective precipitation (evaporation neglected)
! -------------------------------------------------------------------------

  IF ((lconv) .AND. ((ntstep < 2) .OR. (MOD(ntstep+1,nincconv) == 0))) THEN

!  'nseekln' must be < ke_tot*7 (because of 'isendbuflen')
!  'nseekln' must be < MIN( ie_tot/nprocx , je_tot/nprocy )
    nseekln  =  MIN( khumbal , 7*ke , ie_tot/nprocx , je_tot/nprocy )

    nbdlines =  MAX( nboundlines , nseekln )
    naddlin  =  MAX( nseekln - nboundlines , 0 )
    nadd     =  2 * naddlin
    ieadd    =  ie + nadd
    jeadd    =  je + nadd
    ncon     =  (iend - istart + 1) * (jend - jstart + 1)
    ALLOCATE ( nads    (nseekln) , STAT = nzerr )
    ALLOCATE ( mconbas (ncon)    , STAT = nzerr )
    ALLOCATE ( mcontop (ncon)    , STAT = nzerr )
    ALLOCATE ( iix     (ncon)    , STAT = nzerr )
    ALLOCATE ( jix     (ncon)    , STAT = nzerr )
    ALLOCATE ( zconbgp (ieadd,jeadd) , STAT = nzerr )
    ALLOCATE ( zcontgp (ieadd,jeadd) , STAT = nzerr )
    ALLOCATE ( zconbmi (ieadd,jeadd,nseekln+1) , STAT = nzerr )
    ALLOCATE ( zcontmi (ieadd,jeadd,nseekln+1) , STAT = nzerr )

    ic = 0
    DO   j = jstart, jend
      DO i = istart, iend
!       IF ((top_con(i,j) > c05) .AND. (mflx_con(i,j) > epsy)) THEN
!       IF ((top_con(i,j) > c05) .AND. (mflx_con(i,j) > 0.1_wp)) THEN
!       IF ((top_con(i,j) > c05) .AND. (prr_con(i,j)+prs_con(i,j) > epsy)) THEN
!       IF ((top_con(i,j) > c05) .AND. (prne_con(i,j) > epsy)) THEN
        IF (prne_con(i,j) > epsy) THEN
          ic  =  ic + 1
          mconbas (ic)  =  1
          mcontop (ic)  =  ke
          iix     (ic)  =  i
          jix     (ic)  =  j
!   deep (penetrative) or shallow convection: instability from surface
          IF (qcvg_con(i,j) > epsy)  mconbas (ic) = ke
        ENDIF
      ENDDO
    ENDDO
    ncon  =  ic

    DO   kk = ke-1 , 1 , -1
      DO ic = 1 , ncon
        IF (clc_con(iix(ic),jix(ic),kk) > epsy) THEN
!   required for all types of convection
          mcontop (ic)  =  kk
!   required for mid-level convection only;
!   the model layer below convective clouds is also included in the instability
          IF (mconbas(ic) == 1)  mconbas (ic)  =  kk + 1
        ENDIF
      ENDDO
    ENDDO

! Determine vertical extent of convection (with precip)
! -----------------------------------------------------

    DO   j = 1 , jeadd
      DO i = 1 , ieadd
        zconbgp (i,j) = 100000.0_wp
        zcontgp (i,j) =   -400.0_wp
      ENDDO
    ENDDO
    DO ic = 1 , ncon
      IF (mconbas(ic) > 1) THEN
        i = iix(ic) + naddlin
        j = jix(ic) + naddlin
        zconbgp (i,j)  =  c05* (  hhl(iix(ic),jix(ic),mconbas(ic)  )           &
                                + hhl(iix(ic),jix(ic),mconbas(ic)+1))
        zcontgp (i,j)  =  c05* (  hhl(iix(ic),jix(ic),mcontop(ic)  )           &
                                + hhl(iix(ic),jix(ic),mcontop(ic)+1))
      ENDIF
    ENDDO

! Exchange 'boundaries' of vertical extent fields between PE's
! ------------------------------------------------------------

    IF ((num_compute > 1) .AND. (nseekln > 0)) THEN
      jstartex = jstartpar + naddlin
      jendex   = jendpar   + naddlin
      IF (jstartpar == 1 )  jstartex = 1
      IF (jendpar   == je)  jendex   = jeadd

      IF (ltime) THEN
        CALL get_timings (i_spread_ps_pre, ntstep, dt, izerror)
        IF (ltime_barrier) THEN
          CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
          CALL get_timings (i_barrier_waiting_nud, ntstep, dt, izerror)
        ENDIF
      ENDIF

      kzdims(1:24)=(/1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)

      CALL exchg_boundaries ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart   &
                   , num_compute, ieadd, jeadd, kzdims, jstartex, jendex       &
                   , nseekln, nbdlines, my_cart_neigh, lperi_x, lperi_y, l2dim &
                   , 20000+nexch_tag, .FALSE., ncomm_type, izerror, yzerrmsg   &
                   , zconbgp(:,:), zcontgp(:,:) )
!     =====================

      IF (ltime) CALL get_timings (i_communications_nud,ntstep, dt, izerror)
    ENDIF

!-------------------------------------------------------------------------------
!  Section 2: Determination of specified areas around convective grid points
!-------------------------------------------------------------------------------

    istad   =  istart + naddlin
    jstad   =  jstart + naddlin
    iendad  =  iend   + naddlin
    jendad  =  jend   + naddlin

! specify (shape of) area of influence
    ncircle =  NINT( (nseekln + c05) * (nseekln + c05) )
    DO jad = 1 , nseekln
!   squares
!     nads(jad) = nseekln + 1
!   diamonds
!     nads(jad) = nseekln + 1 - jad
!   circles
      nads(jad) = 1
      DO iad = 1 , nseekln
        IF (jad*jad+iad*iad < ncircle)  nads(jad) = iad
      ENDDO
    ENDDO

! compute MIN / MAX over 1,3,5,...,2*nseekln+1 grid points in east-west dir.
    DO   j = jstad-nseekln, jendad+nseekln
      DO i = istad        , iendad
        zconbmi (i,j,1)  =  zconbgp(i,j)
        zcontmi (i,j,1)  =  zcontgp(i,j)
      ENDDO
    ENDDO
    DO   iad = 1            , nseekln
      DO   j = jstad-nseekln, jendad+nseekln
        DO i = istad        , iendad
          zconbmi (i,j,iad+1)  =  MIN( zconbgp(i-iad,j) , zconbmi(i,j,iad)     &
                                     , zconbgp(i+iad,j) )
          zcontmi (i,j,iad+1)  =  MAX( zcontgp(i-iad,j) , zcontmi(i,j,iad)     &
                                     , zcontgp(i+iad,j) )
        ENDDO
      ENDDO
    ENDDO
! compute MIN / MAX over 1,3,5,...,2*nseekln+1 previously computed
! appropriate 'elements' (i.e. east-westerly MIN / MAX) in north-south direction
    DO   i = istad, iendad
      DO j = jstad, jendad
        zconbgp (i,j)  =  zconbmi(i,j,nseekln+1)
        zcontgp (i,j)  =  zcontmi(i,j,nseekln+1)
      ENDDO
    ENDDO
    DO   jad = 1    , nseekln
      DO   i = istad, iendad
        DO j = jstad, jendad
          zconbgp (i,j)  =  MIN( zconbmi(i,j-jad,nads(jad)) , zconbgp(i,j)     &
                               , zconbmi(i,j+jad,nads(jad)) )
          zcontgp (i,j)  =  MAX( zcontmi(i,j-jad,nads(jad)) , zcontgp(i,j)     &
                               , zcontmi(i,j+jad,nads(jad)) )
        ENDDO
      ENDDO
    ENDDO
! shift to local domain
    DO   j = jstart, jend
      DO i = istart, iend
        zconbas (i,j)  =  zconbgp(i+naddlin,j+naddlin)
        zcontop (i,j)  =  zcontgp(i+naddlin,j+naddlin)
      ENDDO
    ENDDO

    DEALLOCATE ( nads    , STAT = nzerr )
    DEALLOCATE ( zconbgp , STAT = nzerr )
    DEALLOCATE ( zcontgp , STAT = nzerr )
    DEALLOCATE ( zconbmi , STAT = nzerr )
    DEALLOCATE ( zcontmi , STAT = nzerr )

! printout for control
! --------------------

    mconv = 0

! single precipitating convective grid pts.
    DO   j = jstart, jend
      DO i = istart, iend
        noffset = (1 - NINT( fr_land(i,j) )) * 7
        mconv(noffset+1) = mconv(noffset+1) + 1
      ENDDO
    ENDDO
    DO ic = 1 , ncon
      i = iix(ic)
      j = jix(ic)
      noffset = (1 - NINT( fr_land(i,j) )) * 7
      IF (mconbas(ic) > 1) THEN
        mconv(noffset+1) = mconv(noffset+1) - 1
        IF (qcvg_con(i,j) > epsy) THEN
          mconv(noffset+2) = mconv(noffset+2) + 1
        ELSEIF (mconbas(ic) > ke-9) THEN
          mconv(noffset+3) = mconv(noffset+3) + 1
        ELSE
          mconv(noffset+4) = mconv(noffset+4) + 1
        ENDIF
        IF (mcontop(ic) > ke-9) THEN
          mconv(noffset+5) = mconv(noffset+5) + 1
        ELSEIF (mcontop(ic) > 15) THEN
          mconv(noffset+6) = mconv(noffset+6) + 1
        ELSE
          mconv(noffset+7) = mconv(noffset+7) + 1
        ENDIF
      ENDIF
    ENDDO

! grid pts. in specified areas around precipitating convective grid pts.
    DO   j = jstart, jend
      DO i = istart, iend
        noffset = (1 - NINT( fr_land(i,j) )) * 7
        noffset = noffset + 14
        IF (zconbas(i,j) > hhl(i,j,2)) THEN
          mconv(noffset+1) = mconv(noffset+1) + 1
        ELSE
          IF (zconbas(i,j) < hhl(i,j,ke)) THEN
            mconv(noffset+2) = mconv(noffset+2) + 1
          ELSEIF (zconbas(i,j) < 1500.0_wp) THEN
            mconv(noffset+3) = mconv(noffset+3) + 1
          ELSE
            mconv(noffset+4) = mconv(noffset+4) + 1
          ENDIF
          IF (zcontop(i,j) < 1500.0_wp) THEN
            mconv(noffset+5) = mconv(noffset+5) + 1
          ELSEIF (zcontop(i,j) < 5000.0_wp) THEN
            mconv(noffset+6) = mconv(noffset+6) + 1
          ELSE
            mconv(noffset+7) = mconv(noffset+7) + 1
          ENDIF
          noffset = 28
          IF (zcontop(i,j)-zconbas(i,j) < 1500.0_wp) THEN
            mconv(noffset+1) = mconv(noffset+1) + 1
          ELSEIF (zcontop(i,j)-zconbas(i,j) < 3000.0_wp) THEN
            mconv(noffset+2) = mconv(noffset+2) + 1
          ELSEIF (zcontop(i,j)-zconbas(i,j) < 5000.0_wp) THEN
            mconv(noffset+3) = mconv(noffset+3) + 1
          ELSE
            mconv(noffset+4) = mconv(noffset+4) + 1
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    DEALLOCATE ( mconbas , STAT = nzerr )
    DEALLOCATE ( mcontop , STAT = nzerr )
    DEALLOCATE ( iix     , STAT = nzerr )
    DEALLOCATE ( jix     , STAT = nzerr )

    IF (num_compute > 1) THEN
      IF (ltime) THEN
        CALL get_timings (i_spread_ps_pre, ntstep, dt, izerror)
        IF (ltime_barrier) THEN
          CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
          CALL get_timings (i_barrier_waiting_nud, ntstep, dt, izerror)
        ENDIF
      ENDIF

      CALL global_values ( mconv, 32, 'SUM', imp_integers, icomm_cart          &
                         , onl_cart_id, yzerrmsg, izerror)
!     ==================

      IF (ltime) CALL get_timings (i_communications_nud, ntstep,dt, izerror)
    ENDIF

    IF (lwonl) THEN
      WRITE( nupr,'(''Step '',I4,'': # Convect. land: NO'',I6,'', BASE'',3I6   &
                  &,'', TOP'',3I6)' ) ntstep, (mconv(ix), ix=1,7)
      WRITE( nupr,'(''Step '',I4,'': # Convect. land: NO'',I6,'', BASE'',3I6   &
                  &,'', TOP'',4I6)' ) ntstep, (mconv(ix), ix=15,21), nseekln
      WRITE( nupr,'(''Step '',I4,'': # Convect. sea : NO'',I6,'', BASE'',3I6   &
                  &,'', TOP'',3I6)' ) ntstep, (mconv(ix), ix=8,14)
      WRITE( nupr,'(''Step '',I4,'': # Convect. sea : NO'',I6,'', BASE'',3I6   &
                  &,'', TOP'',3I6)' ) ntstep, (mconv(ix), ix=22,28)
      WRITE( nupr,'(''Step '',I4,'': # Convect. height: < 1500m / 3000m'',2I6  &
                  &,'', < / > 5000m'',2I6)' ) ntstep, (mconv(ix), ix=29,32)
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure gather_nudge_aux
!-------------------------------------------------------------------------------

END SUBROUTINE gather_nudge_aux


!-------------------------------------------------------------------------------
!+ Module procedure for auxiliary actions mostly related to communication
!-------------------------------------------------------------------------------

SUBROUTINE gather_varia ( nactio , lconai )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure performes various actions, most of which imply the
!   exchange of information between PE's (nodes), or are a consequence of
!   distributed memory architecture. Namely:
!   - checking and adjustment of namelist variables where this cannot be done
!     immediately after reading the namelist
!   - opening / closing of files for control output on the nudging
!   - initialization of other variables
!   - printing control output on the computation of the 'local information'
!     on the observations: - data rejected by threshold quality control
!                          - vertical profiles of observation increments
!                          - interpolated surface-level reports
!
! Method:
!   Straightforward. Use of procedure 'global_values' for the exchange of
!   information between nodes.
!
! Written by        :  Christoph Schraff, DWD  (original version: 09.01.98)
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

! Subroutine arguments:
! --------------------

  INTEGER (KIND=iintegers), INTENT (IN)         ::       &
    nactio          ! determines, which quantities are to be computed (when)

  LOGICAL                 , INTENT (IN) , OPTIONAL  ::   &
    lconai          ! TRUE if analysis increments are const. during 'tconbox'

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    keai             ,& ! number of vertical level in analysis increment fields
    i      , j       ,& ! loop indices in horizontal direction
    iyqc   , iysu    ,& ! loop indices for printing
    kk               ,& ! loop index
    myqcvar          ,& ! type of variable of rejected datum to be printed
    itim             ,& ! (time) index over obs. at one observing station
    ivrs             ,& ! index of spreaded quantity: 1=(u,v); 2=T; 3=RH
    ivar             ,& ! index of spreaded quantity: 1=(u,v); 3=T; 4=RH
    ntimysu          ,& ! time of interpolated surface-level report
    klvi             ,& ! level index in the OIR (observation increment record)
    itmp             ,& ! temporary station index in sorted list
    ist0   , ist1    ,& ! station indices in sorted station list
    nhend            ,& ! number of hours with nudging (dimension of 'noitot')
    nhour            ,& ! loop index
    nulast           ,& ! timestep to complete observation statistics
    ncount (mxobtp+1),& ! counter for reports with obs increments (obs type)
    ncount_pil (4)   ,& ! counter for 'PILOT' reports with obs increments
    ncount_stv (3)   ,& ! counter for satellite retrievals with obs increments
    npsalin          ,& ! number of assimilated surface pressure increments
    ninsitu          ,& ! number of assimilated in-situ reports with obs incr.
    nnscatm          ,& ! number of assimilated scatterometer obs increments
    nalltot          ,& ! total number of reports with obs increments
    maxneed          ,& ! additional array size needed
    nstat  , ierr    ,& ! status variable for opening of files
    nzerr               ! status of memory (de-)allocation

  REAL    (KIND=wp   )     ::  &
    zppr             ,& ! pressure at observation increment level
    zpkemin          ,& ! min. pressure on lowest model level
    zwts1t1, zwts1t2 ,& ! max. temporal weight at station 1: past, future
    zwts0t1, zwts0t2    ! max. temporal weight at station 0: past, future

  LOGICAL                  ::  &
    ltclose          ,& ! files are closed at current timestep
    lchange          ,& ! the order of at least one pair of stations is changed
    lexch               ! exchange current pair of stations


  CHARACTER (LEN=8)     :: yvar       ! variable type
  CHARACTER (LEN=1)     :: yscaldp    ! = 's' if station height diff. is scaled
  CHARACTER (LEN=20)    :: yroutine   ! name of this subroutine
  CHARACTER (LEN=75)    :: yerrmsg    ! error message
  CHARACTER (LEN=75)    :: yerr       ! error message


! Local (automatic) arrays:
! -------------------------

  REAL    (KIND=wp   )     ::  &
    zgather (6)         ! quantities to be collected from all nodes

  INTEGER (KIND=iintegers)  , ALLOCATABLE :: &
    isortys    (:)   ,& ! sorted list of sta. with single-level data for output
    icriterion (:)      ! criterion according to which the index list is sorted

  INTEGER (KIND = iintegers), ALLOCATABLE , SAVE :: &
    noitot   (:,:)    ! (hourly) max. number of stations with active obs. incr.

!
!------------ End of header ----------------------------------------------------


 
!-------------------------------------------------------------------------------
! Begin Subroutine gather_varia
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1a+b: Actions to be done at the beginning of the first nudging cycle
!-------------------------------------------------------------------------------
!  Section 1a: Actions within the local sub-domains
!-------------------------------------------------------------------------------

  IF ((nactio == 1) .AND. (lfirst)) THEN

! allocation of analysis increment fields for long-term storage
! -------------------------------------------------------------

    IF (lnudge) THEN
      keai = ke
      IF (.NOT. lconai) keai = 1

      IF (ALLOCATED( uanai )) THEN
        PRINT '("CAUTION in src_gather_info: uanai is already allocated "      &
              &,"at time / PE",I6,I5)', ntstep, my_cart_id
        DEALLOCATE ( uanai  , vanai  , tanai  , qanai  , STAT=nzerr )
        DEALLOCATE ( psanai , taips  , zconbas, zcontop, STAT=nzerr )
        DEALLOCATE ( rhscal , psaigeo                  , STAT=nzerr )
      ENDIF
      ALLOCATE ( uanai  (ie,je,keai) , STAT = nzerr )
      ALLOCATE ( vanai  (ie,je,keai) , STAT = nzerr )
      ALLOCATE ( tanai  (ie,je,keai) , STAT = nzerr )
      ALLOCATE ( qanai  (ie,je,keai) , STAT = nzerr )
      ALLOCATE ( psanai (ie,je     ) , STAT = nzerr )
      ALLOCATE ( psaigeo(ie,je     ) , STAT = nzerr )
      ALLOCATE ( taips  (ie,je,keai) , STAT = nzerr )
      ALLOCATE ( zconbas(ie,je     ) , STAT = nzerr )
      ALLOCATE ( zcontop(ie,je     ) , STAT = nzerr )
      ALLOCATE ( rhscal (ie,je     ) , STAT = nzerr )

      uanai   = c0
      vanai   = c0
      tanai   = c0
      qanai   = c0
      psanai  = c0
      psaigeo = c0
      taips   = c0
      zconbas = zallhig
      zcontop = zalllow
      rhscal  = c1
    ENDIF

! check and adjustment of namelist variables
! ------------------------------------------

    ionl = MAX( MIN( ionl+0 , ie_tot - nboundlines ) , nboundlines + 1 )
    jonl = MAX( MIN( jonl+0 , je_tot - nboundlines ) , nboundlines + 1 )

! convert global indices  (ionl*,jonl*)  to local indices
    IF (num_compute > 1) THEN
      ionl   =  ionl  - isubpos(my_cart_id,1) + nboundlines + 1
      jonl   =  jonl  - isubpos(my_cart_id,2) + nboundlines + 1
      ionl2  =  ionl2 - isubpos(my_cart_id,1) + nboundlines + 1
      jonl2  =  jonl2 - isubpos(my_cart_id,2) + nboundlines + 1
    ENDIF

! convert hPa to Pa
    DO ivar = 1 , 4
      topobs (ivar) = topobs(ivar) * 100.0_wp
      botmod (ivar) = botmod(ivar) * 100.0_wp
    ENDDO
    ptpstop = ptpstop * 100.0_wp

! check for ptpstop (NOTE that it needs values for pp,p0 !!!
!                    ==> check here rather than in 'input_nudging')
    zpkemin = 109900.0_wp
    DO   j = jstart , jend
      DO i = istart , iend
!       zpkemin = MIN( zpkemin , pp(i,j,ke,nnew) + p0(i,j,ke) )
        zpkemin = MIN( zpkemin , p0(i,j,ke) )
      ENDDO
    ENDDO

! prepare check for (ionl2,jonl2): must be in the same sub-domain as (ionl,jonl)
    lwonl  =        (ionl  <= iend) .AND. (ionl  >= istart)                    &
              .AND. (jonl  <= jend) .AND. (jonl  >= jstart)
    lwonl2 =        (ionl2 <= iend) .AND. (ionl2 >= istart)                    &
              .AND. (jonl2 <= jend) .AND. (jonl2 >= jstart)

! determination of other quantities
! ---------------------------------

    lcroot =  (my_cart_id  == 0)

! Indices for the (full) model levels corresponding to about 250 hPa
! Vertical weights (depend. on model layer only) to geostrophic wind correction
    IF (lnudge) THEN
!     kml250 = 0
      IF (ALLOCATED( wvgeo )) THEN
        PRINT '("CAUTION in src_gather_info: wvgeo is already allocated "      &
              &,"at time / PE",I6,I5)', ntstep, my_cart_id
        DEALLOCATE ( wvgeo  , STAT=nzerr )
      ENDIF
      ALLOCATE ( wvgeo (ke) , STAT = nzerr )
      DO kk = 1 , ke
        wvgeo (kk) = c1
        IF (c05*(vcoord%sigm_coord(kk)+vcoord%sigm_coord(kk+1)) > 0.95_wp)                         &
          wvgeo (kk) = MIN( (c1 - c05 *(vcoord%sigm_coord(kk) + vcoord%sigm_coord(kk+1))) * 20._wp &
                          , c1 )
!       IF (      (vcoord%sigm_coord(kk)  <= 0.25_wp)                                  &
!           .AND. (vcoord%sigm_coord(kk+1) > 0.25_wp)) kml250 = kk
      ENDDO
    ENDIF

! max. lateral boundary relaxation coeff., for reduction of nudging weights
    zmaxrmy (1) = c0
    zmaxrmy (2) = c0
    zmaxrmy (3) = c0
    DO   j = jstartpar , jendpar
      DO i = istartpar , iendpar
        zmaxrmy (1) = MAX( zmaxrmy(1) , rmy(i,j,1) )
        zmaxrmy (2) = MAX( zmaxrmy(2) , rmy(i,j,2) )
        zmaxrmy (3) = MAX( zmaxrmy(3) , rmy(i,j,3) )
      ENDDO
    ENDDO

! initialization for control output of array bound checking for obs. increments 
!   IF ( (num_compute > 1) .AND. (lcroot) .AND. (lnudge)) THEN
    IF ((lcroot) .AND. (lnudge)) THEN
      nhend = INT( REAL(MIN( nstop , nudgend ),wp) * dtdeh + epsy ) + 2
      IF (ALLOCATED( noitot )) THEN
        PRINT '("CAUTION in src_gather_info: noitot is already allocated "     &
              &,"at time / PE",I6,I5)', ntstep, my_cart_id
        DEALLOCATE ( noitot , STAT=nzerr )
      ENDIF
      ALLOCATE ( noitot (nhend,5) , STAT = nzerr )
      noitot (:,:) = 0
    ENDIF

! opening of files used in the nudging (and not only in the obs. processing)
! --------------------------------------------------------------------------

    yroutine = 'gather_varia'
    nstat    = 0
    ierr     = 0

    IF (lwonl) THEN
      OPEN (nupr,FILE=yuprint,FORM='FORMATTED',IOSTAT=nstat)
      nstat = ABS( nstat )
      IF (nstat == 0) REWIND nupr
    ENDIF

!   IF (lcroot) THEN
!     OPEN (nuqc,FILE=yuquctl,FORM='FORMATTED',IOSTAT=nstat)
!     IF (nstat /= 0) THEN
!       yerrmsg = 'OPENING OF FILE yuquctl FAILED'
!       CALL model_abort (my_cart_id, 7001, yerrmsg, yroutine)
!     ENDIF
!     REWIND nuqc
!   ENDIF

!   IF ((lcroot) .AND. (lverif)) THEN
!     OPEN (nuverif,FILE=yuverif,FORM='FORMATTED',IOSTAT=nstat)
!     IF (nstat /= 0) THEN
!       yerrmsg = 'OPENING OF FILE yuverif FAILED'
!       CALL model_abort (my_cart_id, 7001, yerrmsg, yroutine)
!     ENDIF
!     REWIND nuverif
!   ENDIF

!-------------------------------------------------------------------------------
!  Section 1b: Gathering of values to render the quantities globally valid
!-------------------------------------------------------------------------------

    zgather (1) =   zmaxrmy(1)
    zgather (2) =   zmaxrmy(2)
    zgather (3) =   zmaxrmy(3)
    zgather (4) = - zpkemin
    zgather (5) =   c0
    zgather (6) =   c0
    IF ((lwonl) .AND. (.NOT. lwonl2)) zgather (5) = c1
    IF (lwonl) zgather (6) = REAL(my_cart_id, wp)

    IF (num_compute > 1) THEN
      IF (ltime) THEN
        CALL get_timings (i_local_info, ntstep, dt, nzerr)
        IF (ltime_barrier) THEN
          CALL comm_barrier (icomm_cart, nzerr, yerrmsg)
          CALL get_timings (i_barrier_waiting_nud, ntstep, dt, nzerr)
        ENDIF
      ENDIF

      CALL global_values ( zgather, 6, 'MAX', imp_reals, icomm_cart, -1        &
                         , yerrmsg, nzerr)
!     ==================

      CALL global_values ( nstat, 1,'MAX',imp_integers,icomm_cart, -1,yerr,ierr)
!     ==================

      IF (ltime) CALL get_timings (i_communications_nud, ntstep, dt, nzerr)
    ENDIF

    IF ((nstat /= 0) .OR. (ierr /= 0)) THEN
      WRITE( yerrmsg,'("* OPENING OF FILE yuprint FAILED *",2I5)' ) nstat, ierr
      CALL model_abort (my_cart_id, 7001, yerrmsg, yroutine)
!     ================
    ENDIF

    zmaxrmy (1) = zgather(1)
    zmaxrmy (2) = zgather(2)
    zmaxrmy (3) = zgather(3)

! adjust namelist variables
    ptpstop = MIN( ptpstop , 70000.0_wp , -zgather(4) -10000.0_wp )
    IF (zgather(5) > epsy) THEN
      ionl2     = ionl
      jonl2     = jonl
      lwonl2    = lwonl
    ENDIF
    onl_cart_id = NINT(zgather(6))

  ENDIF

!-------------------------------------------------------------------------------
!  Section 1c: Actions done at each timestep
!-------------------------------------------------------------------------------

  IF (nactio == 1) THEN
    ALLOCATE ( a_u  (ie,je,ke) , STAT = nzerr )
    ALLOCATE ( a_v  (ie,je,ke) , STAT = nzerr )
    ALLOCATE ( a_p  (ie,je,ke) , STAT = nzerr )
    ALLOCATE ( a_z  (ie,je,ke) , STAT = nzerr )

    DO kk = 1 , ke
      DO j = 1, je
        DO i = 1, ie
          a_z  (i,j,kk)  =  c05 * (hhl(i,j,kk+1) + hhl(i,j,kk))
          a_p  (i,j,kk)  =  p0(i,j,kk) + pp(i,j,kk,nnew)
          ! for safety:
          a_u  (i,j,kk)  =  u(i,j,kk,nnew)
          a_v  (i,j,kk)  =  v(i,j,kk,nnew)
        ENDDO
      ENDDO
      DO j = jstart, jend
        DO i = istart, iend
          a_u  (i,j,kk)  =  c05 * (u(i,j,kk,nnew) + u(i-1,j,kk,nnew))
          a_v  (i,j,kk)  =  c05 * (v(i,j,kk,nnew) + v(i,j-1,kk,nnew))
        ENDDO
      ENDDO
    ENDDO
  ENDIF

!-------------------------------------------------------------------------------
!  Section 2: Making of a (sorted) list of stations with surface pressure obs.
!             resp. IWV data (for the spatial consistency checks)
!-------------------------------------------------------------------------------

  IF (nactio == 2) THEN

! allocate list and make unsorted lists of the stations
! -----------------------------------------------------

    IF (ALLOCATED( isrtpqc )) THEN
      PRINT '("CAUTION in src_gather_info: isrtpqc is already allocated "      &
            &,"at time / PE",I6,I5)', ntstep, my_cart_id
      DEALLOCATE ( isrtpqc , isrtvqc , STAT=nzerr )
    ENDIF
    ALLOCATE ( isrtpqc (maxpso) , STAT = nzerr )
    ALLOCATE ( isrtvqc (maxivq) , STAT = nzerr )
    DO ista = 1 , npstot
      isrtpqc (ista) = ista
    ENDDO
    DO ista = 1 , nivtot
      isrtvqc (ista) = ista
    ENDDO

! sort stations so that their order is independent from the domain decomposition
! ------------------------------------------------------------------------------
! (sorting with a hybrid quicksort - insertion sort bubblesort algorithm;
!  sorting is according to longitude then latitude (in grid point units);
!  note that the order of co-located reports is always independent from
!  the domain decomposition (except for multi-level aircraft reports).)
! -----------------------------------------------------------------------

! 'surface' pressure data
! -----------------------
    IF ((lreproduce) .AND. (npstot > 1)) THEN
      ALLOCATE ( icriterion (npstot) , STAT = nzerr )
      DO ista = 1 , npstot
        icriterion (ista) = iops_tot(ista) *(je_tot+1) + jops_tot(ista)
      ENDDO

      CALL sortrx ( npstot, icriterion, isrtpqc )
!     ===========

      DEALLOCATE ( icriterion , STAT = nzerr )
    ENDIF

! IWV data
! --------
    IF ((lreproduce) .AND. (nivtot > 1)) THEN
      ALLOCATE ( icriterion (nivtot) , STAT = nzerr )
      DO ista = 1 , nivtot
        icriterion (ista) = ioiv_tot(ista) *(je_tot+1) + joiv_tot(ista)
      ENDDO

      CALL sortrx ( nivtot, icriterion, isrtvqc )
!     ===========

      DEALLOCATE ( icriterion , STAT = nzerr )
    ENDIF

  ENDIF


!-------------------------------------------------------------------------------
!  Section 3: Control output on the computation of the local information
!             Part I: - the results of the threshold quality control
!                     - the single-level data (at first timestep)
!-------------------------------------------------------------------------------

  IF (nactio == 3) THEN

! De-allocatation of geometrical fields used for spatial consistency checks
! -------------------------------------------------------------------------

    DEALLOCATE ( pyyd    , STAT = nzerr )
    DEALLOCATE ( pxxd2   , STAT = nzerr )
    DEALLOCATE ( pxsalpa , STAT = nzerr )

! Output of threshold quality control
! -----------------------------------

    IF ((lcroot) .AND. ((lfirst) .OR. (ntotqc > 0))) THEN
      OPEN (nuqc   , FILE=yuquctl, FORM='FORMATTED', STATUS='UNKNOWN'          &
                                 , POSITION='APPEND', IOSTAT=nstat)
      IF (nstat /= 0) THEN
        yerrmsg = 'OPENING OF FILE yuquctl FAILED'
        CALL model_abort (my_cart_id, 7001, yerrmsg, yroutine)
      ENDIF
    ENDIF

    IF ((ntotqc == maxqcp) .AND. (lcroot)) THEN
      OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'          &
                                 , POSITION='APPEND', IOSTAT=nstat)
      IF (nstat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
      IF (nstat /= 0) CALL model_abort (my_cart_id, 1408, yerr, yroutine)
      WRITE( nucautn,'("CAUTION: PRINTING OF QUALITY CONTROL INCOMPLETE: "     &
                     &,"ntotqc = maxqcp =",I5)' )  ntotqc
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxsgo OR maxmlo")' )
      WRITE( nuqc   ,'("CAUTION: PRINTING OF QUALITY CONTROL INCOMPLETE: "     &
                     &,"ntotqc = maxqcp =",I5)' )  ntotqc
      CLOSE (nucautn)
    ENDIF

    IF ((lcroot) .AND. (lfirst)) THEN
      WRITE( nuqc,'('' '')' )
      WRITE( nuqc,'('' '')' )
      WRITE( nuqc,'(''List of Observations Rejected by the Threshold ''        &
                  &,''Quality Control at the Observation Time'')' )
      WRITE( nuqc,'(''===============================================''        &
                  &,''======================================='')' )
      WRITE( nuqc,'(88X,''Obs.'')' )
      WRITE( nuqc,'(''     Station ID Code Time Pressure Lat. Lon. ''          &
                  &,''Thresh. Var: Obs /Model  Var: Obs /Model  Diff Op.'')' )
      WRITE( nuqc,'('' '')' )
    ENDIF

    IF ((lcroot) .OR. (lwonl)) THEN
      DO iyqc = 1 , ntotqc
        myqcvar = ABS( myqc(iyqc,2) )
        IF (myqcvar ==  1) yvar = 'uv    : '
        IF (myqcvar ==  2) yvar = 'p-TEMP: '
        IF (myqcvar ==  3) yvar = 'T     : '
        IF (myqcvar ==  4) yvar = 'RH    : '
        IF (myqcvar ==  5) yvar = 'uv-10m: '
        IF (myqcvar ==  6) yvar = 'ps    : '
        IF (myqcvar ==  7) yvar = 'T -2m : '
        IF (myqcvar ==  8) yvar = 'RH-2m : '
        IF (myqcvar ==  9) yvar = 'z     : '
        IF (myqcvar == 10) yvar = 'dz    : '
        IF (myqcvar == 11) yvar = 'V-mult: '
        IF (myqcvar == 12) yvar = 'T-mult: '
        IF (myqcvar == 13) yvar = 'q-mult: '
        IF (myqcvar == 14) yvar = 'z-mult: '
        IF (myqcvar == 15) yvar = 'ps-scc: '
        IF (myqcvar == 16) yvar = 'IWV   : '
        IF (myqcvar == 17) yvar = 'IWV-sc: '
        IF (myqcvar == 18) yvar = 'ps-lbc: '
        IF (myqcvar == 19) yvar = 'p-lbc : '
        IF (myqcvar == 20) yvar = 'ps-sbc: '

! printing data rejected at their observation time to a separate file
        IF ((lcroot) .AND. (myqc(iyqc,1) > 0)) THEN
! horizontal wind, upper-air
          IF ((myqcvar == 1) .OR. ((myqcvar == 5) .AND. (itype_tran /= 1))) THEN
            WRITE( nuqc,'(2A , I4, F5.1, 1X,F6.1, 2F6.1, F5.1                  &
                        &,'' : u '',2F7.1, '' , v '',3F6.1)' )                 &
                   yvar, yyqc(iyqc), myqc(iyqc,1), (oyqc(iyqc,i), i=1,10)
! horizontal wind, 10 m
          ELSEIF (myqcvar == 5) THEN
            WRITE( nuqc,'(2A , I4, F5.1, 1X,F6.1, 2F6.1, F5.1                  &
                        &,'' : u '',2F7.1, '' , v '',3F6.1,F4.1)' )            &
                   yvar, yyqc(iyqc), myqc(iyqc,1), (oyqc(iyqc,i), i=1,11)
! 'surface' pressure
          ELSEIF (     (myqcvar == 2 ) .OR. (myqcvar == 6 )                    &
                  .OR. (myqcvar == 18) .OR. (myqcvar == 19)) THEN
            WRITE( nuqc,'(2A , I4, F5.1, 1X,3F6.1, F5.1, '' : ps'',2F7.1)' )   &
                   yvar, yyqc(iyqc), myqc(iyqc,1), (oyqc(iyqc,i), i=1,7)
! 'surface' pressure, spatial consistency check
          ELSEIF (     (myqcvar == 15) .OR. (myqcvar == 20)) THEN
            WRITE( nuqc,'(2A , I4, F5.1, 1X,3F6.1, F5.1, '' : ps'',2F7.1       &
                        &,'' , bias|w2:'',F6.1,F6.1)' )   &
                   yvar, yyqc(iyqc), myqc(iyqc,1), (oyqc(iyqc,i), i=1,9)
! temperature, upper-air
          ELSEIF (myqcvar == 3) THEN
            WRITE( nuqc,'(2A , I4, F5.1, 1X,3F6.1, F5.1, '' : T '',2F7.1)' )   &
                   yvar, yyqc(iyqc), myqc(iyqc,1), (oyqc(iyqc,i), i=1,7)
! temperature, 2 m
          ELSEIF (myqcvar == 7) THEN
            WRITE( nuqc,'(2A , I4, F5.1, 1X,3F6.1, F5.1, '' : T '',3F7.1)' )   &
                   yvar, yyqc(iyqc), myqc(iyqc,1), (oyqc(iyqc,i), i=1,8)
! relative humidity, upper-air
          ELSEIF ((myqcvar == 4) .AND. (myqc(iyqc,2) > 0)) THEN
            WRITE( nuqc,'(2A , I4, F5.1, 1X,3F6.1, F6.2, '': RH '',2F7.2)' )   &
                   yvar, yyqc(iyqc), myqc(iyqc,1), (oyqc(iyqc,i), i=1,7)
! relative humidity due to temperature, upper-air
          ELSEIF (myqcvar == 4) THEN
            WRITE( nuqc,'(2A , I4, F5.1, 1X,3F6.1, F6.2, '': RH '',2F7.2       &
                        &,'', T '',2F6.1)' )                                   &
                   yvar, yyqc(iyqc), myqc(iyqc,1), (oyqc(iyqc,i), i=1,9)
! relative humidity, 2 m
          ELSEIF ((myqcvar == 8) .AND. (myqc(iyqc,2) > 0)) THEN
            WRITE( nuqc,'(2A , I4, F5.1, 1X,3F6.1, F6.2, '': RH '',3F7.2)' )   &
                   yvar, yyqc(iyqc), myqc(iyqc,1), (oyqc(iyqc,i), i=1,8)
! relative humidity due to temperature, 2 m
          ELSEIF (myqcvar == 8) THEN
            WRITE( nuqc,'(2A , I4, F5.1, 1X,3F6.1, F6.2, '': RH '',2F7.2       &
                        &,'', T '',2F6.1, F7.2)' )                             &
                   yvar, yyqc(iyqc), myqc(iyqc,1), (oyqc(iyqc,i), i=1,10)
! hydrostatic height, upper-air
          ELSEIF (myqcvar == 9) THEN
            WRITE( nuqc,'(2A , I4, F5.1, 1X,3F6.1, F6.1                        &
                        &,'' (p-top:'',F6.1,'')    : z '',12X,F6.1)' )         &
                   yvar, yyqc(iyqc), myqc(iyqc,1), (oyqc(iyqc,i), i=1,7)
! hydrostatic thickness, upper-air
          ELSEIF (myqcvar == 10) THEN
            WRITE( nuqc,'(2A , I4, F5.1, 1X,3F6.1, F6.1                        &
                        &,''  p-top:'',F6.1,''  : dz'',F7.1                    &
                        &,'', dT-mean'',F6.1)' )                               &
                   yvar, yyqc(iyqc), myqc(iyqc,1), (oyqc(iyqc,i), i=1,8)
! multi-level check, upper-air (exclude 'z-mult' for wind profilers)
          ELSEIF (     ((myqcvar >= 11) .AND. (myqcvar <= 13))                 &
                  .OR. ((myqcvar == 14) .AND. (myqc(iyqc,1) < nwp_eu))) THEN
            WRITE( nuqc,'(2A , I4, F5.1, 1X,3F6.1, 6X ,''  p-top:'',F6.1)' )   &
                   yvar, yyqc(iyqc), myqc(iyqc,1), (oyqc(iyqc,i), i=1,5)
! integrated water vapour
          ELSEIF (myqcvar == 16) THEN
            WRITE( nuqc,'(2A , I4, F5.1, 1X,3F6.1, F6.2, '':IWV '',2F7.2)' )   &
                   yvar, yyqc(iyqc), myqc(iyqc,1), (oyqc(iyqc,i), i=1,7)
! IWV spatial consistency check
          ELSEIF (myqcvar == 17) THEN
            WRITE( nuqc,'(2A , I4, F5.1, 1X,3F6.1, F5.1, '' :IWV'',2F7.1       &
                        &,'' , bias|w2:'',F6.1,F6.1)' )   &
                   yvar, yyqc(iyqc), myqc(iyqc,1), (oyqc(iyqc,i), i=1,9)
          ENDIF
        ENDIF

! printing data which are rejected for the first time
        IF (lwonl) THEN
! horizontal wind
          IF ((myqcvar == 1) .OR. (myqcvar == 5)) THEN
            WRITE( nupr,'(A,''-QC: '',A ,I4,'' Obs/Mod/Thr'',F5.1,3F6.1,F5.1   &
                        &,'', Time O/M'',2F5.1  ,'', P '',F5.0)' )             &
                   yvar(1:6)   , yyqc(iyqc)  , myqc(iyqc,1)                    &
                 , oyqc(iyqc,6), oyqc(iyqc,8), oyqc(iyqc,7), oyqc(iyqc,9)      &
                 , oyqc(iyqc,5), oyqc(iyqc,1), acthr       , oyqc(iyqc,2)
! 'surface' pressure
          ELSEIF (     (myqcvar == 2 ) .OR. (myqcvar == 6 )                    &
                  .OR. (myqcvar == 18) .OR. (myqcvar == 19)) THEN
            WRITE( nupr,'(A ,''-QC: '',A ,I4,'' Obs/Mod/Thresh'',2F8.2,F6.2    &
                        &,'', Time Obs/Mod'',2F6.1)' )                         &
                   yvar(1:6)   , yyqc(iyqc)  , myqc(iyqc,1), oyqc(iyqc,6)      &
                 , oyqc(iyqc,7), oyqc(iyqc,5), oyqc(iyqc,1), acthr
! 'surface' pressure or IWV, spatial consistency check
          ELSEIF (     (myqcvar == 15) .OR. (myqcvar == 17)                    &
                  .OR. (myqcvar == 20)) THEN
            WRITE( nupr,'(A ,''-QC: '',A ,I4,'' Obs/Mod/Thr'',2F8.2,F6.2       &
                        &,'', T O/M'',2F5.1  ,'', bias|w'',F6.1,F5.1)' )       &
                   yvar(1:6)   , yyqc(iyqc)  , myqc(iyqc,1), oyqc(iyqc,6)      &
                 , oyqc(iyqc,7), oyqc(iyqc,5), oyqc(iyqc,1), acthr             &
                 , oyqc(iyqc,8), oyqc(iyqc,9)
! temperature
          ELSEIF ((myqcvar == 3) .OR. (myqcvar == 7)) THEN
            WRITE( nupr,'(A ,''-QC: '',A ,I4,'' Obs/Mod/Thresh'',3F6.1         &
                        &,'', Time Obs/Mod'',2F6.1  ,'', P '',F5.0)' )         &
                   yvar(1:6)   , yyqc(iyqc)  , myqc(iyqc,1), oyqc(iyqc,6)      &
                 , oyqc(iyqc,7), oyqc(iyqc,5), oyqc(iyqc,1), acthr, oyqc(iyqc,2)
! relative humidity
          ELSEIF ((myqcvar == 4) .OR. (myqcvar == 8)) THEN
            WRITE( nupr,'(A ,''-QC: '',A ,I4,'' Obs/Mod/Thresh'',3F6.2         &
                        &,'', Time Obs/Mod'',2F6.1  ,'', P '',F5.0)' )         &
                   yvar(1:6)   , yyqc(iyqc)  , myqc(iyqc,1), oyqc(iyqc,6)      &
                 , oyqc(iyqc,7), oyqc(iyqc,5), oyqc(iyqc,1), acthr, oyqc(iyqc,2)
! hydrostatic height or thickness
          ELSEIF ((myqcvar == 9) .OR. (myqcvar == 10)) THEN
            WRITE( nupr,'(A ,''-QC: '',A ,I4,'' Differ./Thresh'',2F6.1         &
                        &,'', Time Obs/Mod'',2F6.1,'', P '',F5.0,''-'',F5.0)') &
                   yvar(1:6)   , yyqc(iyqc)  , myqc(iyqc,1), oyqc(iyqc,7)      &
                 , oyqc(iyqc,5), oyqc(iyqc,1), acthr, oyqc(iyqc,2), oyqc(iyqc,6)
! multi-level check
          ELSEIF ((myqcvar >= 11) .AND. (myqcvar <= 14)) THEN
            WRITE( nupr,'(A ,''-QC: '',A ,I4, 27X                              &
                        &,'', Time Obs/Mod'',2F6.1,'', P '',F5.0,''-'',F5.0)') &
                   yvar(1:6)   , yyqc(iyqc)  , myqc(iyqc,1), oyqc(iyqc,1)      &
                 , acthr       , oyqc(iyqc,2), oyqc(iyqc,5)
! integrated water vapour
          ELSEIF (myqcvar == 16) THEN
            WRITE( nupr,'(A ,''-QC: '',A ,I4,'' Obs/Mod/Thresh'',3F6.2         &
                        &,'', Time Obs/Mod'',2F6.1  ,'', P '',F5.0)' )         &
                   yvar(1:6)   , yyqc(iyqc)  , myqc(iyqc,1), oyqc(iyqc,6)      &
                 , oyqc(iyqc,7), oyqc(iyqc,5), oyqc(iyqc,1), acthr, oyqc(iyqc,2)
          ENDIF
        ENDIF

      ENDDO
    ENDIF
    DEALLOCATE ( oyqc    , STAT = nzerr )
    DEALLOCATE ( myqc    , STAT = nzerr )
    DEALLOCATE ( yyqc    , STAT = nzerr )

! close YUQUCTL file
    IF ((lcroot) .AND. ((lfirst) .OR. (ntotqc > 0)))   CLOSE (nuqc)

!   IF (lwonl .AND. ldump_ascii) THEN
!     CLOSE (nupr)
!     OPEN  (nupr,FILE=yuprint,FORM='FORMATTED',STATUS='OLD'                   &
!                             ,POSITION='APPEND',IOSTAT=nstat)
!     IF (nstat /= 0) yerrmsg = 'OPENING OF FILE yuprint FAILED'
!     IF (nstat /= 0) CALL model_abort (my_cart_id, 1009, yerrmsg, yroutine)
!   ENDIF

! (quick-)sort messages of single-level data so that their order
! for output is independent from the domain decomposition
! --------------------------------------------------------------

    IF (ALLOCATED( isortys )) THEN
      PRINT '("CAUTION in src_gather_info: isortys is already allocated "      &
            &,"at time / PE",I6,I5)', ntstep, my_cart_id
      DEALLOCATE ( isortys , STAT=nzerr )
    ENDIF
    ALLOCATE ( isortys (maxysu) , STAT = nzerr )

    DO ista = 1 , ntotys
      isortys (ista) = ista
    ENDDO

    IF (      (lreproduce) .AND. (ntotys > 1)                                  &
        .AND. ((lfirst) .OR. (rmod( ntstep*dt,c3600,epsy ) < tconbox))) THEN
      ALLOCATE ( icriterion (ntotys) , STAT = nzerr )
      DO ista = 1 , ntotys
        icriterion (ista) =   NINT( oysu(ista,18) ) *(je_tot+1)                &
                            + NINT( oysu(ista,19) )
      ENDDO

      CALL sortrx ( ntotys, icriterion, isortys )
!     ===========

      DEALLOCATE ( icriterion , STAT = nzerr )
    ENDIF

! printing surface-level observations interpolated to the lowest model level
! and observation increments of single-level reports
! --------------------------------------------------------------------------

    IF ((lfirst) .OR. (rmod( ntstep*dt,c3600,epsy ) < tconbox)) THEN
      IF (lwonl) THEN
        WRITE( nupr,'('' '')' )
        WRITE( nupr,'(''Interpolation of surface-level observations to the ''  &
                    &,''lowest model level ke'')' )
        WRITE( nupr,'(''===================================================''  &
                    &,''====================='')' )
        WRITE( nupr,'(''(upper-air single-level reports are also listed her''  &
                    &,''e)'')' )
        WRITE( nupr,'(''(station height differences with the flag "s" are s''  &
                    &,''caled as for the extrapolation of pressure)'')' )
        WRITE( nupr,'(''       |  obs: observed values / ke: obs. interpola''  &
                    &,''ted to level ke / inc: obs. increments at ke |'')' )
        WRITE( nupr,'(''Sta.   | height  | surf. pressure |    10m - horizo''  &
                    &,''ntal wind    | 2m-temperature | 2m-rel humid.|t'')' )
        WRITE( nupr,'('' id    | diff. * | obs / ke / inc | 10m-obs /  at k''  &
                    &,''e  /   inc   | obs / ke / inc | obs/ ke /inc | '')' )
        DO itmp = 1 , ntotys
          iysu = isortys (itmp)
          IF ((oysu(iysu,3) > rmdich) .AND. (oysu(iysu,3) < c0)) THEN
            oysu (iysu,3)  = - oysu(iysu,3)
            yscaldp        = 's'
          ELSE
            yscaldp        = ' '
          ENDIF
          ntimysu = NINT( MAX( oysu(iysu,17) , -9.0_wp ) )
          WRITE( nupr,'(A,1X,A,F5.0, 2F7.1,F5.1, 6F5.1, 2F6.1,F5.1, 3F5.2,I2)')&
                 yysu(iysu), yscaldp, (oysu(iysu,i), i=1,16), ntimysu
        ENDDO
        IF (ldump_ascii) THEN
          ! flush YUPRINT file
          CLOSE (nupr)
          OPEN  (nupr,FILE=yuprint,FORM='FORMATTED',STATUS='OLD'               &
                                  ,POSITION='APPEND',IOSTAT=nstat)
          IF (nstat /= 0) yerrmsg = 'OPENING OF FILE yuprint FAILED'
          IF (nstat /= 0) CALL model_abort (my_cart_id, 1009, yerrmsg, yroutine)
        ENDIF
      ENDIF
      DEALLOCATE ( oysu    , STAT = nzerr )
      DEALLOCATE ( yysu    , STAT = nzerr )
    ENDIF
    DEALLOCATE ( isortys , STAT = nzerr )

    IF ((.NOT. ltnudge) .AND. (lwonl)) THEN
      WRITE( nupr,'(" Verification mode: observations are quality controled")' )
!     PRINT  *    , ' Verification mode: observations are quality controled'
    ENDIF

! printing number of local upper-air obs incr. not used due to ODR size
! ---------------------------------------------------------------------

    IF ((lcroot) .AND. (nuaex > 0)) THEN
      OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'          &
                                 , POSITION='APPEND', IOSTAT=nstat)
      IF (nstat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
      IF (nstat /= 0) CALL model_abort (my_cart_id, 1408, yerr, yroutine)
      WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," LOCAL UP-AIR SING-L."    &
                     &," OBS. INCR. BEYOND maxuso ",I5)' ) ntstep, nuaex, maxuso
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxuso BY AT LEAST" &
                     &,I6)' ) nuaex
      CLOSE (nucautn)
    ENDIF

  ENDIF


!-------------------------------------------------------------------------------
!  Section 4: Actions to be taken prior to the spreading of obs. information
!-------------------------------------------------------------------------------

  IF (nactio == 4) THEN

!-------------------------------------------------------------------------------
!  Section 4a: Making of (sorted) lists of stations for all types of data
!-------------------------------------------------------------------------------

! allocate lists of the stations
! ------------------------------

    IF (ALLOCATED( isortps )) THEN
      PRINT '("CAUTION in src_gather_info: isortps is already allocated "      &
            &,"at time / PE",I6,I5)', ntstep, my_cart_id
      DEALLOCATE ( isortps , isortsu , isortua , isortml , STAT=nzerr )
    ENDIF
    ALLOCATE ( isortps (maxpso) , STAT = nzerr )
    ALLOCATE ( isortsu (maxsgo) , STAT = nzerr )
    ALLOCATE ( isortua (maxuso) , STAT = nzerr )
    ALLOCATE ( isortml (maxmloi_tot) , STAT = nzerr )

! make unsorted lists of the stations
! -----------------------------------

    DO ista = 1 , npstot
      isortps (ista) = ista
    ENDDO
    DO ista = 1 , nsutot
      isortsu (ista) = ista
    ENDDO
    DO ista = 1 , nuatot
      isortua (ista) = ista
    ENDDO
    DO ista = 1 , nmltot
      isortml (ista) = ista
    ENDDO

! sort stations so that their order is independent from the domain decomposition
! ------------------------------------------------------------------------------
! (iterative sorting with bubblesort algorithm; sorting is according to horizon-
!  tal coordinates. Note that the order of co-located reports is always indepen-
!  dent from the domain decomposition (except for multi-level aircraft reports.)
! ------------------------------------------------------------------------------

! 'surface' pressure data
! -----------------------
    IF ((lreproduce) .AND. (npstot > 1)) THEN
      DO ista = 1 , npstot
        isortps (ista) = isrtpqc (ista)
      ENDDO
    ENDIF

! (other) surface-level data
! --------------------------
    IF ((lreproduce) .AND. (nsutot > 1)) THEN
      ALLOCATE ( icriterion (nsutot) , STAT = nzerr )
      DO ista = 1 , nsutot
        icriterion (ista) = iosu_tot(ista) *(je_tot+1) + josu_tot(ista)
      ENDDO

      CALL sortrx ( nsutot, icriterion, isortsu )
!     ===========

      DEALLOCATE ( icriterion , STAT = nzerr )
    ENDIF

! upper-air single-level data
! ---------------------------
    IF ((lreproduce) .AND. (nuatot > 1)) THEN
      ALLOCATE ( icriterion (nuatot) , STAT = nzerr )
      DO ista = 1 , nuatot
        icriterion (ista) = ioua_tot(ista) *(je_tot+1) + joua_tot(ista)
      ENDDO

      CALL sortrx ( nuatot, icriterion, isortua )
!     ===========

      DEALLOCATE ( icriterion , STAT = nzerr )
    ENDIF

! multi-level data
! ----------------
    IF ((lreproduce) .AND. (nmltot > 1)) THEN
      ALLOCATE ( icriterion (nmltot) , STAT = nzerr )
      DO ista = 1 , nmltot
        icriterion (ista) =   ioml_tot(ista) *2*(je_tot+1)                     &
                            + joml_tot(ista) *2                                &
                            + MAX( 1-ABS( kobtyml(ista)-nairep ) ,0 )
      ENDDO

      CALL sortrx ( nmltot, icriterion, isortml )
!     ===========

      DEALLOCATE ( icriterion , STAT = nzerr )

! sort co-located multi-level aircraft reports
      lchange = .TRUE.
      DO WHILE (lchange)
        lchange = .FALSE.
        DO ista = 1 , nmltot-1
          IF (      (ioml_tot(isortml(ista+1)) == ioml_tot(isortml(ista)))     &
              .AND. (joml_tot(isortml(ista+1)) == joml_tot(isortml(ista)))     &
              .AND. (kobtyml (isortml(ista  )) == nairep)                      &
              .AND. (kobtyml (isortml(ista+1)) == nairep)) THEN
            lexch = .FALSE.
            ist0  = isortml(ista)
            ist1  = isortml(ista+1)
            zwts1t1 = MAX( zwtml(ist1,1,1) ,zwtml(ist1,1,2), zwtml(ist1,1,3) )
            zwts0t1 = MAX( zwtml(ist0,1,1) ,zwtml(ist0,1,2), zwtml(ist0,1,3) )
            zwts1t2 = MAX( zwtml(ist1,2,1) ,zwtml(ist1,2,2), zwtml(ist1,2,3) )
            zwts0t2 = MAX( zwtml(ist0,2,1) ,zwtml(ist0,2,2), zwtml(ist0,2,3) )
            IF (      (ABS( zwts1t1 - zwts0t1 ) <= epsy)                       &
                .AND. (ABS( zwts1t2 - zwts0t2 ) <= epsy)) THEN
! the (one) report of both stations is done at the same time
              itim = 2
              IF (MAX( zwtml(ist0,1,1) ,zwtml(ist0,1,2) ) > epsy) itim = 1
              IF (mszlev(ist1,itim) == mszlev(ist0,itim)) THEN
! stations have the same number of vertical levels
! if loop 'check_oiml' is not quitted by the EXIT, the two stations contain
! identical reports (except possibly for the station identity and for 'fcorlml')
                check_oiml: DO ivrs = 1 , maxnoi
                  DO kk = 0 , mszlev(ist0,itim) - 1
                    IF (ABS( oiml(ivrs,kkoiml(ist0,itim)+kk)                   &
                            -oiml(ivrs,kkoiml(ist1,itim)+kk) ) >= epsy) THEN
                      IF (oiml(ivrs,kkoiml(ist1,itim)+kk) >=                   &
                          oiml(ivrs,kkoiml(ist0,itim)+kk) + epsy) lexch = .TRUE.
                      EXIT check_oiml
                    ENDIF
                  ENDDO
                ENDDO check_oiml
              ELSEIF (mszlev(ist1,itim) > mszlev(ist0,itim)) THEN
! exchange stations if the 2nd station contains more vertical levels
                lexch = .TRUE.
              ENDIF
            ELSEIF (MAX( zwts1t1,zwts1t2 ) >= MAX( zwts0t1,zwts0t2 )+epsy) THEN
! exchange stations if the 2nd contains a larger maximum temporal weight
              lexch = .TRUE.
            ELSEIF (   (MAX( zwts1t1,zwts1t2 ) >= MAX( zwts0t1,zwts0t2 )-epsy) &
                 .AND. (zwts1t2 > zwts0t2+epsy)) THEN
! exchange stations if the 2nd station has a future report with the same
! maximum temporal weight as the past report of the 1st station
              lexch = .TRUE.
            ENDIF
            IF (lexch) THEN
              itmp             = isortml (ista+1)
              isortml (ista+1) = isortml (ista)
              isortml (ista)   = itmp
              lchange = .TRUE.
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDIF

!-------------------------------------------------------------------------------
!  Section 4b: Control output on the computation of the local information
!              Part II: - the multi-level data (at first timestep)
!                       - the number of currently active reports or stations
!-------------------------------------------------------------------------------

! printing vertical profiles of observation increments
! ----------------------------------------------------

    IF ((lfirst) .AND. (ltnudge) .AND. (lwonl)) THEN
      WRITE( nupr,'('' '')' )
      WRITE( nupr,'(''Vertical profiles of observation increments from TE''    &
                  &,''MP, PILOT, AIRCRAFT'')' )
      WRITE( nupr,'(''===================================================''    &
                  &,''==================='')' )
      DO itmp = 1 , nmltot
        ista = isortml (itmp)
        DO itim = 1, 2
          IF ((kkoiml(ista,itim) > 0) .AND. (mszlev(ista,itim) > 0)) THEN
            ystid = ystidml(ista) (1:ilstidp)
            WRITE( nupr,'(''station '',A ,'', ('',I3,'','',I3,''),'',I3,       &
                         &'' levels ; mbotlv, mtoplv '',6I4)')                 &
                   ystid, ioml_tot(ista), joml_tot(ista), mszlev(ista,itim)    &
                 , (mbotlv(ista,itim,ivrs),ivrs=1,3)                           &
                 , (mtoplv(ista,itim,ivrs),ivrs=1,3)
            IF (kobtyml(ista) /= ngps) THEN
              WRITE( nupr,'(''  u-incr  v-incr  T-incr   Qd-incr RH-incr''     &
                          &,'' pressure height pot. T. '',A )' )   ystid
              DO kk = 1 , mszlev(ista,itim)
                klvi  = kk + kkoiml(ista,itim) - 1
                zppr = EXP( oiml(noilp,klvi) )
                WRITE( nupr,'(3F8.3, F10.5,F8.3, 2F8.0,F7.2)' )                &
                       oiml(noiu ,klvi), oiml(noiv ,klvi)                      &
                     , oiml(noit ,klvi), oiml(noiqd,klvi)                      &
                     , oiml(noirh,klvi), zppr                                  &
                     , oiml(noiz ,klvi), oiml(noith,klvi)
              ENDDO
            ELSE
              WRITE( nupr,'(''        quality weight    Qd-incr RH-incr''      &
                          &,'' pressure height pot. T. '',A )' )   ystid
              DO kk = 1 , mszlev(ista,itim)
                klvi  = kk + kkoiml(ista,itim) - 1
                zppr = EXP( oiml(noilp,klvi) )
                WRITE( nupr,'(8X,F8.5,8X, F10.5,F8.3, 2F8.0,F7.2)' )           &
                       oiml(noiqqc,klvi), oiml(noiqd,klvi)                     &
                     , oiml(noirh ,klvi), zppr                                 &
                     , oiml(noiz  ,klvi), oiml(noith,klvi)
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDIF

! printing IWV from GPS reports
! -----------------------------

    IF ((     (ntstep <= 1  ) .OR. (ntstep == 48) .OR. (ntstep == 138)         &
         .OR. (ntstep == 228)) .AND. (ltnudge) .AND. (lwonl)) THEN
      WRITE( nupr,'('' '')' )
      WRITE( nupr,'(''Integrated water vapour IWV values from processing ''    &
                  &,''of ground-based GPS'')' )
      WRITE( nupr,'(''===================================================''    &
                  &,''==================='')' )
      WRITE( nupr,'("timestep:",I5)' )  ntstep
      WRITE( nupr,'("    obs  #  satu-  obs IWV:   adjusted:    IWV      m"    &
                  &,"odel w.  satu- extra-            ")' )
      WRITE( nupr,'("   time it- rated  repor-  re-     ice,   from  model"    &
                  &,"  cloud  rated polat.            ")' )
      WRITE( nupr,'("STA  hr erat. lev    ted trieved   bias qv-obs    IWV"    &
                  &,"    ice  model quot  p-int  p-rep")' )
      DO itmp = 1 , nmltot
        ista = isortml (itmp)
        DO itim = 1, 2
          IF (      (kkoiml(ista,itim) > 0) .AND. (mszlev(ista,itim) > 0)      &
              .AND. (kobtyml(ista) == ngps) .AND. (ke >= 13)) THEN
            ystid = ystidml(ista) (1:ilstidp)
            ist0  = NINT( oiml(noiv,kkoiml(ista,itim)+1) )
            ist1  = NINT( oiml(noiv,kkoiml(ista,itim)+2) )
            kk    = kkoiml(ista,itim) + 3
            WRITE( nupr,'(A5, F4.1,2I4,1X, 7F7.2,F5.1,2F7.1)' )                &
                   ystid, oiml(noiv,kkoiml(ista,itim)), ist0, ist1             &
                        ,(oiml(noiv,klvi), klvi=kk,kk+9)
!             equiv.to: ,tobs, iter, isat, ziqgps, ziqint, ziqobs, ziqite,     &
!                        ziqmod, ziqmic, ziqsat, zfracb, pint, pgps
!           remove temporal misuse of 'noiv' for printout
            oiml(noiv,kkoiml(ista,itim):kkoiml(ista,itim)+12) = rmdi
          ENDIF
        ENDDO
      ENDDO
    ENDIF

! printing number of currently active reports or stations
! -------------------------------------------------------

!   IF (num_compute == 1) THEN
      ncount     = 0
!     ncounts    = 0
      ncount_pil = 0
      ncount_stv = 0
      DO ista = 1 , nmltot
        ! this 'IF' condition is required here to exclude GPS obs, which are set
        ! passive in the threshold QC but still used in the spatial consistency
        ! check (SCC), or which are set passive in the SCC
        IF (MAXVAL( zwtml(ista,:,:) ) > epsy) THEN
          ncount (kobtyml(ista))  =  ncount(kobtyml(ista)) + 1
          ncount (mxobtp+1)       =  ncount(mxobtp+1)      + 1
        ENDIF
      ENDDO
      IF (ncount(npilot) > 0) THEN
        DO ista = 1 , nmltot
          IF (      (kobtyml(ista) == npilot)                                   &
              .AND. (MAXVAL( zwtml(ista,:,:) ) > epsy)) THEN
            IF (     (kcdtyml(ista) == nldpcd) .OR. (kcdtyml(ista) == nshpcd)   &
                .OR. (kcdtyml(ista) == nmopcd)) THEN
              ncount_pil (1) = ncount_pil(1) + 1
            ELSEIF ((kcdtyml(ista) == nwp_eu) .OR.(kcdtyml(ista) == npr_us)) THEN
              ncount_pil (2) = ncount_pil(2) + 1
            ELSEIF  (kcdtyml(ista) == nra_eu) THEN
              ncount_pil (3) = ncount_pil(3) + 1
            ELSEIF  (kcdtyml(ista) == nravad) THEN
              ncount_pil (4) = ncount_pil(4) + 1
            ENDIF
          ENDIF
        ENDDO
      ENDIF
      IF (ncount(nsattv) > 0) THEN
        DO ista = 1 , nmltot
          IF (      (kobtyml(ista) == nsattv)                                   &
              .AND. (MAXVAL( zwtml(ista,:,:) ) > epsy)) THEN
            IF (    (kcdtyml(ista) == nnoa15) .OR.(kcdtyml(ista) == nnoa16)     &
               .OR. (kcdtyml(ista) == nnoa17) .OR.(kcdtyml(ista) == nnoa18)) THEN
              ncount_stv (1) = ncount_stv(1) + 1
            ELSEIF ((kcdtyml(ista) == nsmsg1) .OR.(kcdtyml(ista) == nsmsg2)) THEN
              ncount_stv (2) = ncount_stv(2) + 1
!           ELSEIF  (kcdtyml(ista) == nmetop)                                THEN
!             ncount_stv (3) = ncount_stv(3) + 1
            ENDIF
          ENDIF
!         IF     (     (ystidml(ista) (1:3) == 'msg')                           &
!                 .OR. (ystidml(ista) (1:3) == 'MSG')) THEN
!           ncounts (1) = ncounts(1) + 1
!         ELSEIF (     (ystidml(ista) (1:3) == 'noaa')                          &
!                 .OR. (ystidml(ista) (1:3) == 'NOAA')) THEN
!           ncounts (2) = ncounts(2) + 1
!         ELSEIF (     (ystidml(ista) (1:3) == 'metop')                         &
!                 .OR. (ystidml(ista) (1:3) == 'METOP')) THEN
!           ncounts (3) = ncounts(3) + 1
!         ENDIF
        ENDDO
      ENDIF
      npsalin = 0
      DO ista = 1 , npstot
        ! this 'IF' condition is required here to exclude PS obs, which are set
        ! passive in the threshold QC but still used in the spatial consistency
        ! check (SCC), or which are set passive in the SCC
        IF (MAX( omykps(ista,1),omykps(ista,2) ) > epsy)  npsalin = npsalin + 1
      ENDDO
      ninsitu = 0
      nnscatm = 0
      DO ista = 1 , nsutot
        IF (kobtysu(ista) == nscatt)  nnscatm = nnscatm + 1
        IF (kobtysu(ista) /= nscatt)  ninsitu = ninsitu + 1
      ENDDO
      nalltot = ncount(mxobtp+1) + nuatot + nsutot + npsalin
      IF (lcroot) THEN
        PRINT          '(4X,"STEP",I5,": NUMBER OF (SINGLE OR PAIRS OF) "      &
                       &,"REPORTS WITH OBS INCREMENTS:",I6)' ,  ntstep, nalltot
        IF (nalltot >= 1) THEN
          PRINT        '(13X,I6," multi-level:",I4," TEMP,",I3," PILOT,"       &
                       &,I4," WINDPROF,",I4," RADAR-VAD")' ,                   &
                     ncount(mxobtp+1), ncount(ntemp), ncount_pil(1)            &
                                     , ncount_pil(2), ncount_pil(4)
          PRINT        '(32X,I4," AIRCRAFT,",I2," RASS,",I5," GPS,"            &
                       &,I5," RETRIEVALS")' ,                                  &
                     ncount(nairep), ncount_pil(3), ncount(ngps), ncount(nsattv)
          IF (ncount(nsattv) > 0)                                              &
            PRINT      '(21X,"RETRIEVALS:",I4," ATOVS,",I5," SEVIRI,"          &
                       &,I4," IASI")' ,                                        &
                     ncount_stv(1), ncount_stv(2), ncount_stv(3)
          PRINT        '(13X,I6," sing-lev aircraft,",I5," in-situ surface,"   &
                       &,I5," surf. pressure")' ,                              &
                     nuatot, ninsitu, npsalin
          IF (nnscatm > 0)                                                     &
            PRINT      '(13X,I6," scatterometer")' , nnscatm
        ENDIF
!       PRINT       '(I4,": # STATIONS WITH:",I5," multi-level:"               &
!                   &,I4," TEMP,",I4," PILOT,",I4," AMDAR,",I4," GPS")' ,      &
!              ntstep, nmltot, ncount(ntemp), ncount(npilot), ncount(nairep)   &
!                            , ncount(ngps)
!       IF (ncount(nsattv) > 0)                                                &
!         PRINT       '(I4,"                  ",I5," satellite:"               &
!                     &,I5," MSG,",I4," NOAA,",I4," METOP")' ,                 &
!                ncount(nsattv), ncounts(1), ncounts(2), ncounts(3)
!       PRINT       '(4X,"   OBS INCREMENTS:",I5," upper-air sing-lv,"         &
!                   &,I5," surface,",I5," surf. pressure")' ,                  &
!              nuatot, ninsitu, npsalin
      ENDIF
      IF (lwonl ) THEN
        WRITE ( nupr , '(4X,"STEP",I5,": NUMBER OF (SINGLE OR PAIRS OF) "      &
                       &,"REPORTS WITH OBS INCREMENTS:",I6)' )  ntstep, nalltot
        IF (nalltot >= 1) THEN
          WRITE ( nupr,'(13X,I6," multi-level:",I4," TEMP,",I3," PILOT,"       &
                       &,I4," WINDPROF,",I4," RADAR-VAD")' )                   &
                     ncount(mxobtp+1), ncount(ntemp), ncount_pil(1)            &
                                     , ncount_pil(2), ncount_pil(4)
          WRITE ( nupr,'(32X,I4," AIRCRAFT,",I2," RASS,",I5," GPS,"            &
                       &,I5," RETRIEVALS")' )                                  &
                     ncount(nairep), ncount_pil(3), ncount(ngps), ncount(nsattv)
          IF (ncount(nsattv) > 0)                                              &
            WRITE(nupr,'(21X,"RETRIEVALS:",I4," ATOVS,",I5," SEVIRI,"          &
                       &,I4," IASI")' )                                        &
                     ncount_stv(1), ncount_stv(2), ncount_stv(3)
          WRITE ( nupr,'(13X,I6," sing-lev aircraft,",I5," in-situ surface,"   &
                       &,I5," surf. pressure")' )                              &
                     nuatot, ninsitu, npsalin
          IF (nnscatm > 0)                                                     &
            WRITE(nupr,'(13X,I6," scatterometer")' ) nnscatm
        ENDIF
!       WRITE( nupr,'("t=",I4,": # OBS INCREMENTS (STATIONS):",I5," multi-lev" &
!                   &,"el:",I4," TEMP,",I4," PILOT,",I4," AMDAR,",I4," GPS")') &
!              ntstep, nmlall, ncount(ntemp), ncount(npilot), ncount(nairep)   &
!                            , ncount(ngps)
!       IF (ncount(nsattv) > 0)                                                &
!         WRITE( nupr,'("t=",I4,": # OBS INCREMENTS (STATIONS):",I5," satelli" &
!                     &,"te:",I5," MSG,",I4," NOAA,",I4," METOP")' )           &
!                ntstep, ncount(nsattv), ncounts(1), ncounts(2), ncounts(3)
!       WRITE( nupr,'("t=",I4,": # OBS INCREMENTS (STATIONS):",I5," upper-air" &
!                   &," sing-lv,",I5," surface,",I5," surf. pressure")' )      &
!              ntstep, nuaall, nsuall, npsalin
      ENDIF
      ! flush YUPRINT file
      IF ((lwonl) .AND. (ldump_ascii)) THEN
        CLOSE (nupr)
        OPEN  (nupr,FILE=yuprint,FORM='FORMATTED',STATUS='OLD'                 &
                                ,POSITION='APPEND',IOSTAT=nstat)
        IF (nstat /= 0) yerrmsg = 'OPENING OF FILE yuprint FAILED'
        IF (nstat /= 0) CALL model_abort (my_cart_id, 1008, yerrmsg, yroutine)
      ENDIF
!   ENDIF

! storing hourly maximum number of active reports or stations
! -----------------------------------------------------------

!   IF ( (num_compute > 1)  .AND. (lcroot)                                     &
!        .AND. (lnudge) .AND. (ntstep <= MIN( nudgend,nstop ))) THEN
    IF ((lcroot) .AND. (lnudge) .AND. (ntstep <= MIN( nudgend,nstop ))) THEN
      nhour = INT( ntstep * dtdeh ) + 1
      IF (num_compute > 1) THEN
!       noitot (nhour,1) = MAX( noiall , noitot(nhour,1) )
        noitot (nhour,2) = MAX( nmlall , noitot(nhour,2) )
        noitot (nhour,3) = MAX( nuaall , noitot(nhour,3) )
        noitot (nhour,4) = MAX( nsuall , noitot(nhour,4) )
        noitot (nhour,5) = MAX( npsall , noitot(nhour,5) )
!       noitot (nhour,6) = MAX( nqcall , noitot(nhour,6) )
      ELSE
!       noitot (nhour,1) = MAX( nmloit , noitot(nhour,1) )
        noitot (nhour,2) = MAX( nmltot , noitot(nhour,2) )
        noitot (nhour,3) = MAX( nuatot , noitot(nhour,3) )
        noitot (nhour,4) = MAX( nsutot , noitot(nhour,4) )
        noitot (nhour,5) = MAX( npstot , noitot(nhour,5) )
!       noitot (nhour,6) = MAX( ntotqc , noitot(nhour,6) )
      ENDIF
    ENDIF

! flush YUPRINT file
    IF ((lwonl) .AND. (lfirst) .AND. (ldump_ascii)) THEN
      CLOSE (nupr)
      OPEN  (nupr,FILE=yuprint,FORM='FORMATTED',STATUS='OLD'                   &
                              ,POSITION='APPEND',IOSTAT=nstat)
      IF (nstat /= 0) yerrmsg = 'OPENING OF FILE yuprint FAILED'
      IF (nstat /= 0) CALL model_abort (my_cart_id, 1009, yerrmsg, yroutine)
    ENDIF

  ENDIF

!-------------------------------------------------------------------------------
!  Section 5: Control output for obs. increment array bounds, and cleanup
!-------------------------------------------------------------------------------

  IF (nactio == 5) THEN

    nulast  =  MIN( INT( MAX( nverend , nudgend ) ,iintegers) , nstop )
!   IF (       (num_compute > 1)  .AND. (lcroot)                               &
!        .AND. (lnudge) .AND. (ntstep == nulast)) THEN
    IF ((lcroot) .AND. (lnudge) .AND. (ntstep == nulast)) THEN
      yroutine = 'gather_varia'
      OPEN (nustat,FILE=yustats,FORM='FORMATTED',STATUS='OLD'                  &
                               ,POSITION='APPEND',IOSTAT=nstat)
      IF (nstat /= 0) THEN
        yerrmsg = 'OPENING OF FILE yustats FAILED'
        CALL model_abort (my_cart_id, 7005, yerrmsg, yroutine)
      ENDIF
      WRITE( nustat,'(''0'')' )
      WRITE( nustat,'(''1'')' )
      WRITE( nustat,'(''0 ---HOURLY MAX. TOTAL NUMBER OF STATIONS WITH ''      &
                    &,''ACTIVE OBSERVATION INCREMENTS'')' )
      WRITE( nustat,'(''+                 (only for MPP applications)'')' )
      WRITE( nustat,'(''+ ---                    multi-   |''                  &
                    &,''  upper-air   | surface- | surface '')' )
      WRITE( nustat,'(''+                        level    |''                  &
                    &,'' single-level |  level   | pressure'')' )
      nhend = INT( REAL(MIN( nstop , nudgend ),wp) * dtdeh + epsy ) + 1
      DO nhour = 1 , nhend
        WRITE( nustat,'(''+'',5X,I2,''. hour : '',3X,I10,I14,I13,I11)' )       &
               nhour, (noitot(nhour,iyqc), iyqc=2,5)
        DO kk = 2 , 5
          noitot (nhend+1,kk) = MAX( noitot(nhend+1,kk) , noitot(nhour,kk) )
        ENDDO
      ENDDO
      nhour = nhend + 1
      WRITE( nustat,'(''+     total max.  : '',  I10,I14,I13,I11)' )           &
                    (noitot(nhour,iyqc), iyqc=2,5)
      WRITE( nustat,'(''+     array bounds: '',  I10,I14,I13,I11)' )           &
             maxmloi_tot, maxuso, maxsgo, maxpso
      IF (     (noitot(nhour,2) > maxmloi_tot) .OR. (noitot(nhour,3) > maxuso) &
          .OR. (noitot(nhour,4) > maxsgo)  .OR. (noitot(nhour,5) > maxpso)) THEN
        WRITE( nustat,'(''1'')' )
        WRITE( nustat,'(''0     !!! CAUTION !!!!! CAUTION !!!!! ''             &
                      &,''CAUTION !!!!! CAUTION !!!!! CAUTION !!!!!'')' )
        WRITE( nustat,'(''+     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~''             &
                      &,''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'')' )
      ENDIF
      IF (noitot(nhour,2) > maxmloi_tot) THEN
        WRITE( nustat,'(''+     WARNING: array size for multi-level obs. ''    &
                      &,''increments is'', /                                   &
                      &,''      =======  too small to accommodate all ''       &
                      &,''obs. increments'')' )
        maxneed = maxmlo +    noitot(nhour,2) - maxmloi_tot
        WRITE( nustat,'(''        -->  Increase "MAXMLO" (namelist) from''     &
                      &,I4 ,'' to at least'',I4)' )                            &
               maxmlo, maxneed
        maxneed = maxgpo + 2*(noitot(nhour,2) - maxmloi_tot)
        WRITE( nustat,'(''          or increase "MAXGPO" (namelist) from''     &
                      &,I4 ,'' to at least'',I4)' )                            &
               maxgpo, maxneed
        maxneed = maxtvo +    noitot(nhour,2) - maxmloi_tot
        WRITE( nustat,'(''          or increase "MAXTVO" (namelist) from''     &
                      &,I4 ,'' to at least'',I4, /                             &
                      &,''             ================='')' )                 &
               maxtvo, maxneed
      ENDIF
      IF (noitot(nhour,3) > maxuso)                                            &
        WRITE( nustat,'(''+     WARNING: array size for upper-air single-''    &
                      &,''level obs. increments is'', /                        &
                      &,''      =======  too small to accommodate all ''       &
                      &,''obs. increments'', /                                 &
                      &,''        -->  Increase "MAXUSO" (namelist) from''     &
                      &,I5 ,'' to at least'',I5, /                             &
                      &,''             ================='')' )                 &
               maxuso, noitot(nhour,3)
      IF ((noitot(nhour,4) > maxsgo) .OR. (noitot(nhour,5) > maxpso)) THEN
        noitot (nhour,5) = MAX(      noitot(nhour,4)                           &
                              , INT( noitot(nhour,5)-maxmlo , iintegers ) )
        WRITE( nustat,'(''+     WARNING: array size for surface-level obs''    &
                      &,''. increments is'', /                                 &
                      &,''      =======  too small to accommodate all ''       &
                      &,''obs. increments'', /                                 &
                      &,''        -->  Increase "MAXSGO" (namelist) from''     &
                      &,I5 ,'' to at least'',I5, /                             &
                      &,''             ================='')' )                 &
               maxsgo, noitot(nhour,5)
      ENDIF
      DEALLOCATE ( noitot , STAT = nzerr )
      CLOSE( nustat )
    ENDIF

    DEALLOCATE ( a_u   , STAT = nzerr )
    DEALLOCATE ( a_v   , STAT = nzerr )
    DEALLOCATE ( a_p   , STAT = nzerr )
    DEALLOCATE ( a_z   , STAT = nzerr )

    IF ((lnudge) .AND. (ntstep == MIN( nudgend,nstop ))) THEN
      DEALLOCATE ( uanai  , STAT = nzerr )
      DEALLOCATE ( vanai  , STAT = nzerr )
      DEALLOCATE ( tanai  , STAT = nzerr )
      DEALLOCATE ( qanai  , STAT = nzerr )
      DEALLOCATE ( psanai , STAT = nzerr )
      DEALLOCATE ( psaigeo, STAT = nzerr )
      DEALLOCATE ( taips  , STAT = nzerr )
      DEALLOCATE ( zconbas, STAT = nzerr )
      DEALLOCATE ( zcontop, STAT = nzerr )
      DEALLOCATE ( rhscal , STAT = nzerr )
      DEALLOCATE ( wvgeo  , STAT = nzerr )
    ENDIF

    ltclose = (ntstep == MIN(INT( MAX( nverend , nudgend ),iintegers) , nstop ))

    IF ((lcroot) .AND. (ltclose)) CLOSE( nuqc )
    IF ((lwonl ) .AND. (ltclose)) CLOSE( nupr )

!   IF ((lcroot) .AND. (lverif) .AND. (ntstep == MIN( nverend,nstop ))) THEN
!     kk = 0
!     WRITE( nuverif,'(''  -9 xxxxx     '', I6,I5,3I6, I2,I4, I11,I11,I2,I3    &
!                    &, 2I5)' )   kk,kk,kk,kk,kk,kk,kk,kk,kk,kk,kk,kk,kk
!     CLOSE( nuverif )
!   ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure gather_varia
!-------------------------------------------------------------------------------

END SUBROUTINE gather_varia

!===============================================================================

ELEMENTAL REAL (KIND=wp) FUNCTION rmod  ( ztti, ziv, epsy )
  !----------------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  ztti, ziv, epsy
  !---------------------------------------------------------------------------
  ! MOD function for positive reals
  !---------------------------------------------------------------------------
  !
  rmod  =  ztti  -  ziv * REAL(INT( ztti/ziv + epsy ),wp) + epsy
  !
END FUNCTION rmod

!-------------------------------------------------------------------------------

END MODULE src_gather_info
