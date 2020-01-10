!+ Source module  "src_sfcana"
!-------------------------------------------------------------------------------

MODULE src_sfcana

!-------------------------------------------------------------------------------
!
! Description:
!   The module "sfcana" performs the analysis of near surface parameters,
!      such as screen-level temperature and humidity, 10m wind speed and
!      precipitation total.
!      The method used is
!   a) successive correction of the model fields with available observations
!      in case of 10 m wind speed, screen-level temperature and humitity analysis
!   b) simple distance weighted interpolation of observations in case of
!      precipitation analysis.
!   Driving routine is the module procedure "organize_sfcana"
!
! Current Code Owner: DWD, Michael Buchhold 
!  phone:  +49  69  8062 2726
!  fax:    +49  69  8062 3721
!  email:  michael.buchhold@.dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.22       1999/02/08 Michael Buchhold
!  Initial release
! 1.27       1999/03/29 Christoph Schraff
!  Revised ODR format and missing data indicator. All include decks included
!  in the module source deck. Allocatable arrays used in several module
!  procedures moved to 'data_obs_record'.
! 1.29       1999/05/11 Ulrich Schaettler
!  Adapted interfaces to utility-modules, prepared use of MPE_IO and 
!  corrected calls to intrinsics for use on non-Cray machines
! 1.31       1999/07/01 Christoph Schraff
!  MPI calls related to fixed-length arrays replaced by parallel_utility calls.
!  Bug corrections: memory deallocation for all arrays; double declaration of
!  'nobtot' removed.
! 1.34       1999/12/10 Ulrich Schaettler
!  Allocate and deallocate grib arrays in organize_sfcana
! 1.35       1999/12/16 Ulrich Schaettler
!  Bug correction: Allocation of grib fields after all possible RETURNs
! 1.37       2000/03/24 Michael Buchhold
!  Correction for printout
! 1.39       2000/05/03 Ulrich Schaettler
!  Changed names for variables concerned to latitude or longitude.
! 1.43       2000/07/18 Jean-Marie Bettems
!  Corrections for use on NEC SX5: call of 'copen', increase of 'maxobs'.
! 2.4        2001/01/29 Christoph Schraff
!  Reorganisation of loops for enhanced efficiency on vector processors.
! 2.7        2001/06/26 Michael Buchhold
!  Bug correction: Number of obs may not exceed array size. 
!  Data selection: Skip passive observations
! 2.8        2001/07/06 Ulrich Schaettler
!  Eliminated non-necessary variables from the USE-lists and dependency for 
!  combine_subarrays
! 2.11       2001/09/28 Ulrich Schaettler
!  Corrected a bug when gribbing the pole of the rotation grid (igds_out(21))
! 2.13       2002/01/18 Michael Buchhold
!  Replace subroutine for sorting observation array. Close analysis GRIB file.
!  Replace DATA statements by PARAMETER lists. Check precip obs time range.
! 2.14       2002/02/15 Ulrich Schaettler
!  Bug correction: Use of 'ibits' corrected.
!  Correction in the grib-coding of the upper right corner (igds_out(11))
! 2.17       2002/05/08 Ulrich Schaettler
!  Implemented a check for the range of longitude variables
! 2.19       2002/10/24 Michael Buchhold
!  Introduce analysis of 10m wind speed.
!  Slight modification of observation weight function and
!  observation and background error estimates
! 3.5        2003/09/02 Ulrich Schaettler
!  Adapted interface for routine exchg_boundaries
! 3.7        2004/02/18 Ulrich Schaettler
!  Renamed phi (rlat), rla (rlon);
!  Changed treatment of unit-numbers for ASCII-file handling
! 3.12       2004/09/15 Christoph Schraff
!  Bug correction on timing of closing files.
! 3.13       2004/12/03 Ulrich Schaettler
!  Put KIND-parameters for Grib library to data_parameters
! 3.15       2005/03/03 Christoph Schraff
!  Size of obs arrays made dependent on (old) namelist param. 'maxsgo' (LME).
!  Replaced FLOAT by REAL (Ulrich Schaettler)
! 3.16       2005/07/22 Michael Buchhold
!  Erroneous sorting of abs array fixed. Subroutine ia_comp modified.
!  Sort on station ID cancelled.
!  Size of obs array modified.
!  Number of sucsessive correction scans increased for t2m and rf2m
! 3.17       2005/12/12 Michael Buchhold
!  Changed unit-of-timerange for Grib records to 1 (instead of 0 before)
! 3.18       2006/03/03 Christoph Schraff / Ulrich Schaettler
!  Avoid model crash due to insufficient size of arrays 'indxi_', 'indxj_' ...
!  Use of lyear_360 for changed interface to get_utc_date; 
!  Adaptations to changed names in data_io: ds => ds_grib
!  Determine length of grib record with idims_out(19) (length in bytes)
! 3.19       2006/04/25 Ulrich Schaettler
!  Corrected specification of igds_out(9)
! 3.21       2006/12/04 Ulrich Schaettler
!  Put statement functions which are not used to a comment
! V4_4         2008/07/16 Ulrich Schaettler
!  Changed NL parameter lyear_360 to itype_calendar, to have several options
! V4_5         2008/09/10 Ulrich Schaettler
!  Adaptations for new reference atmosphere
! V4_8         2009/02/16 Guenther Zaengl
!  Modify grib encoding of new reference atmosphere for more flexibility
!  Bug fix for grib encoding of LM output
!  Add l_ke_in_gds to partly replace ldwd_grib_use
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_17        2011/02/24 Ulrich Blahak
!  Adapted interface of exchg_boundaries; corrected kzdims(1:20) -> kzdims(1:24)
! V4_19        2011/08/01 Ulrich Schaettler
!  Introduced conditional compilation for GRIBDWD
! V4_22        2012/01/31 Christoph Schraff 
!  Some variables moved from 'data_nudge_all' to 'data_obs_lib_cosmo'.
! V4_23        2012/05/10 Oliver Fuhrer
!  Removed obsolete Fortran features
! V4_24        2012/06/22 Hendrik Reich
!  Adapted length of strings for date variables
! V4_25        2012/09/28 Hendrik Reich, Ulrich Schaettler
!  Read also the minutes out of the date string
!  Introduced nexch_tag for MPI boundary exchange tag to replace ntstep (HJP)
!  Adapted control of when to do sfc-analysis to time steps that do not fit
!    into a full hour (US)
! V4_28        2013/07/12 Ulrich Schaettler, Christoph Schraff
!  Implemented grib_api for writing the surface analysis
!  Adapted interface to grib_api routines with special grib_api integer
!  Use subroutines make_grib_grid, make_grib_product from new module io_metadata
!   to set GRIB meta data (US)
!  Direct call of mpi routine replaced by call of 'gather_values'.
!  Statement functions replaced by elemental or intrinsic functions. (CS)
! V4_30        2013/11/08 Ulrich Schaettler
!  Renamed ipds to ipds_out to reflect usage for output
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
!  Remove explicit calls to open files code and use io_utilities.f90 instead
! V5_3         2015-10-09 Hans-Juergen Panitz
!  Put ifdef GRIBDWD around declaration of EXTERNAL irefts
!  Set pp_lansfc%lsfc_ana in SR sfcout depending on ktr (tri) input value
!   to distinguish between the different SP_10M fields
! V5_4         2016-03-10 Christoph Schraff
!  - Precipitation analysis according to 'raintp' instead of 12 hours fixed.
!  - Old index used for AOF read is replaced by 'nbswwe' (for NetCDF read).
!  - Variables related to the AOF interface removed.
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code  + DWD Extensions
!===============================================================================
!
! Declarations:
!
! Modules used:
!
!-------------------------------------------------------------------------------

USE data_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    iintegers, & ! KIND-type parameter for standard integer variables
    irealgrib, & ! KIND-type parameter for real variables in the grib library
    iwlength,  & ! length of an integer word in byte
    intgribf,  & ! KIND-type parameter for fortran files in the grib library
    intgribc,  & ! KIND-type parameter for C files in the grib library
    int_ga       ! integer precision for grib_api: length of message in bytes

!-------------------------------------------------------------------------------

USE data_parallel,      ONLY :  &
    num_compute,     & ! number of compute PEs
    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains
    ldatatypes,      & ! if .TRUE.: use MPI-Datatypes for some communications
    ltime_barrier,   & ! if .TRUE.: use additional barriers for determining the
                       ! load-imbalance
    ncomm_type,      & ! type of communication
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    isubpos,         & ! positions of the subdomains in the total domain. Given
                       ! are the i- and the j-indices of the lower left and the
                       ! upper right grid point in the order
                       !                  i_ll, j_ll, i_ur, j_ur.
                       ! Only the interior of the domains are considered, not
                       ! the boundary lines.
    icomm_cart,      & ! communicator for the virtual cartesian topology
    my_cart_neigh,   & ! neighbors of this subdomain in the cartesian grid
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    imp_integers,    & ! determines the correct INTEGER type used in the model
                       ! for MPI
    iexch_req,       & ! stores the sends requests for the neighbor-exchange
    nexch_tag,       & ! tag to be used for MPI boundary exchange
                       !  (in calls to exchg_boundaries)
    sendbuf,         & ! sending buffer for boundary exchange:
                       !   1-4 are used for sending, 5 is used for receiving
                       ! both buffers are allocated in organize_setup
    isendbuflen        ! length of one column of sendbuf

! end of data_parallel

!--------------------------------------------------------------------------

USE data_io,            ONLY :  &
    ldwd_grib_use,& ! use some DWD specific Grib settings
    l_ke_in_gds,  & ! explicit GDS entry for number of model levels
    npds,         & ! Dimension for product definition section (pds)
    ngds,         & ! Dimension for grid description section (gds)
    nbms,         & ! Dimension for bit map section (bms)
    nbds,         & ! Dimension for binary data section
    ndsup,        & ! Dimension for dsup
    ndims,        & ! Dimension for idims (contains all dime
    lfd,          & !
    lfa,          & ! Dimension for grib_api message in bytes
    lbm,          & !
    lds,          & !
    idwdednr,     & ! grib edition number for DWD library
    undefgrib,    & ! value for "undefined" in the grib routines
    ncenter,      & ! originating center
    nprocess_ini_in, & ! generating process identification for analysis data

!   Global arrays
    iblock,       & ! array for gribed data
    ymessage,     & ! array for grib-api message (in characters)
    idims_out,    & ! array for all dimensions
    ibmap,        & ! array for
    ipds_out,     & ! product definition section for output
    igds_out,     & ! grid description section
    ibms,         & ! bit map section
    ibds,         & ! binary data section
    dsup,         & ! Parameter for grib routines
    ds_grib,      & ! array for unpacked data
    igrib1_id,    & ! grib1 sample
    igrib2_id,    & ! grib1 sample
    var,          & !
    num_gribtabs, & ! number of GRIB tables used in LM variable table
    pp_nl,        & ! structure for gribout namelist
    pp_lansfc,    & ! output group for near surface analysis fields

    ydate_ini       ! start of the forecast
                    ! yyyymmddhhmmss (year, month, day, hour, min., sec.)

! end of data_io

!-----------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------

    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction
    ie_tot,       & ! number of grid points in zonal direction
    je_tot,       & ! number of grid points in meridional direction
    ie_max,       & ! Max. of ie on all processors
    je_max,       & ! Max. of je on all processors

! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from
!    the other ones because of the use of the staggered Arakawa-C-grid.
!
    jstartpar,    & ! start index for computations in the parallel program
    jendpar,      & ! end index for computations in the parallel program


! 4. constants for the horizontal rotated grid and related variables
! ------------------------------------------------------------------

    pollon,       & ! longitude of the rotated north pole (in degrees, E>0)
    pollat,       & ! latitude of the rotated north pole (in degrees, N>0)
    polgam,       & ! angle between the north poles of the systems
    dlon,         & ! grid point distance in zonal direction (in degrees)
    dlat,         & ! grid point distance in meridional direction (in degrees)
    startlat_tot, & ! transformed latitude of the lower left grid point
                    ! of the total domain (in degrees, N>0)
    startlon_tot, & ! transformed longitude of the lower left grid point
                    ! of the total domain (in degrees, E>0)
    endlon_tot,   & ! transformed longitude of the upper right grid point
                    ! of the total domain (in degrees, E>0)
    endlat_tot,   & ! transformed latitude of the upper right grid point
                    ! of the total domain (in degrees, N>0)
    startlon,     & ! transformed longitude of the lower left grid point
                    ! of this subdomain (in degrees, E>0)
    startlat,     & ! transformed latitude of the lower left grid point
                    ! of this subdomain (in degrees, N>0)
    degrad,       & ! factor for transforming degree to rad

! 5. variables for the time discretization and related variables
! --------------------------------------------------------------

    dt,           & ! long time-step
    dtdeh           ! dt / 3600 seconds

! end of data_modelconfig

!---------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 1. mathematical constants
! -------------------------

    pi,           & ! circle constant

! 2. physical constants and related variables
! -------------------------------------------

!   t0_melt,      & ! absolute zero for temperature
    r_d,          & ! gas constant for dry air
    g,            & ! acceleration due to gravity
    r_earth,      & ! mean radius of the earth

! 3. constants for parametrizations
! ---------------------------------

    b1,           & ! variables for computing the saturation steam pressure
    b2w,          & ! over water (w) and ice (e)
    b3,           & !               -- " --
    b4w             !               -- " --

! end of data_constants


!-------------------------------------------------------------------------------
USE data_fields     , ONLY :   &

! 1. fields for model output and diagnostics                          (unit )
! ------------------------------------------

    t_2m       ,    & ! temperature in 2m                             (  k  )
    qv_2m      ,    & ! specific water vapor content in 2m            (kg/kg)
    td_2m      ,    & ! dew-point in 2m                               (  k  )
    u_10m      ,    & ! zonal wind in 10m                             ( m/s )
    v_10m      ,    & ! meridional wind in 10m                        ( m/s )


! 2. external parameter fields                                        (unit)
! ----------------------------
    hsurf         , & ! height of surface topography
    fr_land       , & ! fraction of land in a grid element              --
    rlat          , & ! geographical latitude                         ( rad )
    rlon              ! geographical longitude                        ( rad )

! end of data_fields

!-------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

    itype_calendar,&! for specifying the calendar used
    nstart,       & ! first time step of the forecast
    nstop,        & ! last time step of the forecast
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
    nnew,         & ! corresponds to ntstep + 1
    lperi_x,      & ! if lartif_data=.TRUE.: periodic boundary conditions in x-dir.
                    ! or with Davies conditions (.FALSE.)
    lperi_y,      & ! if lartif_data=.TRUE.: periodic boundary conditions in y-dir.
                    ! or with Davies conditions (.FALSE.)
    l2dim,        & ! 2 dimensional runs
    idbg_level      ! to control the verbosity of debug output

! end of data_runcontrol

!-------------------------------------------------------------------------------

USE data_nudge_all , ONLY :   &

! 1. Parameters and related variables
! -----------------------------------

    lwonl        ,& ! .TRUE if grid pt. (ionl ,jonl ) lies in the local domain

! 2. Namelist variables controlling the data assimilation
! -------------------------------------------------------

    nudgend      ,& ! 0      : end of nudging period in timesteps
    nverend      ,& ! 0      : end of verification period in timesteps
    maxsgo       ,& ! 3000   : max. number of (surface-level and upper-air)
                    !                         single-level reports in the ODR
    ionl         ,& ! 167    : / grid point coordinates
    jonl         ,& ! 103    : \ for standard output on nudging
    lsurfa       ,& !  if surface fields to be analysed
    lt2m         ,& ! .f.    : .t. if 2m temperat. field is analysed
    lrh2m        ,& ! .f.    : .t. if 2m rel. hum. field is analysed
    lprecp       ,& ! .f.    : .t. if precipitation is analysed
    lff10m       ,& ! .f.    : .t. if 10m wind speed is analysed
    ht2a         ,& ! 999.   : time of 1. T2m-ana in hours since model start
    ht2i         ,& ! 999.   : time increment to next T2m analysis
    ht2next      ,& ! 999.   : next hour, when T2m-ana shall be done
    hh2a         ,& ! 999.   : time of 1. RH2m-ana in hours since model start
    hh2i         ,& ! 999.   : time increment to next RH2m analysis
    hh2next      ,& ! 999.   : next hour, when RH2m-ana shall be done
    hprc         ,& ! 999.   : time of prec-ana in hours since model start
    hprcnext     ,& ! 999.   : time of prec-ana in hours since model start
    hffa         ,& ! 999.   : time of 1. 10m wind speed analysis
    hffi         ,& ! 999.   : time increment to next wind speed analysis
    hffnext      ,& ! 999.   : next hour, when wind-ana shall be done
    nt2next      ,& !        : next time step, when T2m-ana shall be done
    nh2next      ,& !        : next time step, when RH2m-ana shall be done
    nffnext      ,& !        : next time step, when wind-ana shall be done
    nprcnext     ,& !        : next time step, when prec-ana shall be done
    raintp       ,& ! 12.    : length of precipitation period to be analysed
    ydir_lansfc  ,& ! './'   : directory where to write the 2-D analyses
    yform_lansfc ,& ! 'grb1' : format for the 2-D analyses files
    ldiasa       ,& ! .f.    : .t. for diagnostics of surface analysis

! 3.2 Device numbers
! ------------------

    nsfc

USE data_nudge_all , ONLY :   &

! 4.1 Surface analysis limits
! ---------------------------

    rmaxdp       ,& ! max. value for scaled normalised departure

! 4.2 Allocable arrays used for surface analysis
! ----------------------------------------------

    alat         ,& ! model latitudes
    alon         ,& ! model longitudes
    anal         ,& ! global analysis field
    firstg       ,& ! global first guess field
    rsurin       ,& ! increment field
    deprow       ,& ! observations or observed increments for each model row
    rpalto       ,& ! original station latitude
    rpalno       ,& ! original station longitude
    rpalat       ,& ! station latitude (assigned model g.p.)
    rpalon       ,& ! station longitude (assigned model g.p.)
    rcslon       ,& ! cosine of station longitude
    rsnlon       ,& ! sine of station longitude
    rcslat       ,& ! cosine of station latitude 
    rsnlat       ,& ! sine of station latitude 
    rpahgt       ,& ! station height   
    rpaobe       ,& ! observation error
    wa           ,& ! weight * observation increment
    wwa          ,& ! weights
    bufall       ,& ! buffer array containing observations from all PEs
    noatab       ,& ! index of the first obs in a model row
    indexv       ,& ! index describing the state of the analysis at a g.p.
    nnodep       ,& ! number of observations per row
    npagpt       ,& ! sea/land type of assigned grid point
    nuprow       ,& ! northern row of local data selection area
    ndnrow       ,& ! southern row of local data selection area
    yidall       ,& ! buffer of station ids of all observations from all PEs
    ypasid          ! station id

! end of data_nudge_all

!-------------------------------------------------------------------------------
    
USE data_obs_lib_cosmo , ONLY :   &

! 1. General parameters
! ---------------------

    rmdi       ,& ! =-1.E31_wp : commonly used missing data indicator
    rmdich     ,& ! =-1.E30_wp : commonly used check value for miss data
    epsy       ,& ! = 1.E-8_wp : commonly used very small value > 0 

! 4. I/O device numbers for obs processing / nudging and file names
! -----------------------------------------------------------------

    nupr       ,& ! unit number of file for all the remaining information

! 5. CMA observation type and code type numbers
! ---------------------------------------------

    nsynop     ,& ! SYNOP reports
    ndribu        ! DRIBU reports

! end of data_obs_lib_cosmo

!-------------------------------------------------------------------------------

USE data_obs_record , ONLY :   &

! 1.1 ODR header format
! ---------------------

    nhilon       ,& ! longitude of observing station
    nhjlat       ,& ! latitude  of observing station
    nhalt        ,& ! station altitude [m]
    nhtime       ,& ! time of observat. in forecast hours
    nhio         ,& ! (local) x-coord. of grid pt. assigned to obs
    nhjo         ,& ! (local) y-coord. of grid pt. assigned to obs
    nhitot       ,& ! global x-coord. of grid pt. assigned to obs
    nhjtot       ,& ! global y-coord. of grid pt. assigned to obs
    nhobtp       ,& ! observation type
    nhcode       ,& ! code type
    nhschr       ,& ! station characteristics                      (see 1.1.4)
    nhpass       ,& ! flag for report being set to 'passive'       (see 1.1.4)

! 1.2 Bit patterns for packed information in ODR header
! -----------------------------------------------------

    nvsebp       ,& ! bit pos. for report located at sea grid pt.     nvhsch

! 1.3 ODR body format
! -------------------

    nbsu         ,& ! u wind component [m/s]
    nbsv         ,& ! v wind component [m/s]
    nbst         ,& ! temperature [k]
    nbsrh        ,& ! relative humidity [/]
    nbsrr1       ,& ! precipitation amount over 1  hour                  [mm]
    nbsrr3       ,& ! precipitation amount over 3  hours                 [mm]
    nbsrr6       ,& ! precipitation amount over 6  hours                 [mm]
    nbsr12       ,& ! precipitation amount over 12 hours                 [mm]
    nbsr24       ,& ! precipitation amount over 24 hours                 [mm]
    nbswwe       ,& ! SYNOP: weather and ground group word  (below)

! 1.4 Bit patterns for packed information in ODR body
! ---------------------------------------------------

    nrtrbp       ,& ! bit position for code of precipitation
                    !     measurement duration       [Code table 4019, keys 0-7]
    nrtroc       ,& ! no. bits occupied by precip obs. duration code

! 1.5 Further quantities related to ODR
! -------------------------------------

    imdi         ,& ! missing data indicator for ODR integers (2^31-1)
    ntotsg       ,& ! tot. number of stored single-level reports
    ilstid       ,& ! character length of the station identity
    ilstidp      ,& ! char. length used for printing the station ID
                    ! Note: (ilstid >= ilstidg >= ilstidp) cf. data_nudge_gather

! 2. Observation data records (ODR)
! ---------------------------------

    osgbdy       ,& ! body of single-level ODR
    osghed       ,& ! header of single-level ODR
    mosghd       ,& ! header of single-level ODR
    mosgbd       ,& ! body of single-level ODR
    yosghd          ! header of single-level ODR

! end of data_obs_record

!---------------------------------------------------------------------------

USE utilities,                ONLY :  &
    get_utc_date                !

!-------------------------------------------------------------------------------

USE io_utilities, ONLY :   &
    open_file, close_file

!-------------------------------------------------------------------------------

USE environment,              ONLY :  &
    model_abort,              & ! aborts the program in case of errors
    exchg_boundaries,         & !
    comm_barrier                ! explicit synchronization point

!-------------------------------------------------------------------------------

USE parallel_utilities,       ONLY :  &
    gather_field,    & ! gathers the parts of a total field from all subdomains
    gather_all,      & ! gathers at all nodes a set of arrays from all nodes
    gather_values,   & ! gathers a set of values from all nod. to 1 or all nodes
    combine_subarrays  ! 

!-------------------------------------------------------------------------------

USE io_metadata,              ONLY :  &
    make_grib_init,  & ! setting constant grib meta data
    make_grib_grid,  & ! setting grib meta data for grid definition
    make_grib_product  ! setting grib meta data for product definition

!===============================================================================

#ifdef GRIBAPI
! grib_api interface
USE grib_api
#endif

!===============================================================================

IMPLICIT NONE
 
!===============================================================================

  CHARACTER (LEN=14)     ::                           &
    yandat1         ! analysis date in the form   yyyymmddhhmmss

   
  REAL (KIND=wp)           , PARAMETER  :: &
    tmxdif(4) = (/ 0.251_wp   ,& ! accepted time diff. between obs time
                   0.251_wp   ,& ! and analysis time  of the parameters:
                   0.01_wp    ,& ! t2m, rh2m, precipitation, ff10m
                   0.251_wp /),& ! 
    rmxhds    =    400._wp    ,& ! max. height difference (obs.-model)
                                     ! for precipitation observations
    rmxhdt    =    300._wp    ,& ! max. height difference (obs.-model)
                                     ! for T2m and RH2m observations
    rmxhdv    =    300._wp    ,& ! max. height difference (obs.-model)
                                     ! for ff10m observations
    rmprcs(3) = (/ 40000._wp  ,& ! precipitation scan radius in m.
                   70000._wp  ,& 
                  110000._wp/),& 
    rmt2ms(3) = (/200000._wp  ,& ! 2m temp. scan radius in m
                  100000._wp  ,&
                   50000._wp/),&
    rmr2ms(3) = (/200000._wp  ,& ! 2m hum. scan radius in m.
                  100000._wp  ,&
                   50000._wp/),&
    rmfffs(3) = (/200000._wp  ,& ! 10m wind  scan radius in m.
                  100000._wp  ,&
                   50000._wp/)

 

  INTEGER (KIND=iintegers) , PARAMETER  :: &
   nt2m    =  1  ,& ! 1; variable number for t2m
   nh2m    =  2  ,& ! 2; variable number for rh2m
   nprc    =  3  ,& ! 3; variable number for prc
   nff10m  =  4  ,& ! 4; variable number for ff10m
   nprcsc  =  3  ,& ! number of prec. successive correction scans
   nt2msc  =  3  ,& ! number of 2m temp. successive correction scans
   nr2msc  =  3  ,& ! number of 2m hum. successive correction scans
   nfffsc  =  3     ! number of 10m wind successive correction scans

   INTEGER (KIND=iintegers)                 , PRIVATE :: &
   isfco         ,& ! descriptor for surface analysis GRIB output
   nobtot           ! total number of observations gathered from all PEs


! 2. Variables
! ------------

!   Observation buffer format
!   -------------------------
   
  INTEGER (KIND=iintegers) , PARAMETER  :: &
   nboxnr =  1 ,& ! model grid point number; j_tot*10000+i_tot
   noblon =  2 ,& ! observation longitude
   noblat =  3 ,& ! observation latitude
   nstalt =  4 ,& ! station altitude
   nindsl =  5 ,& ! type of assigned grid point, land/sea indicator
   nt2mob =  6 ,& ! observed screen level temperature
   nt2oer =  7 ,& ! error of observed temperature
   nh2mob =  8 ,& ! observed screen level rel humidity
   nh2oer =  9 ,& ! error of observed rel. humidity
   nprcob = 10 ,& ! observed precipitation amount
   nprcer = 11 ,& ! error of observed precipitation
   nf10ob = 12 ,& ! observed 10m wind speed
   nf10er = 13 ,& ! error of observed wind speed
   nidpe  = 14 ,& ! identification of the PE the observation belongs to
   nobnum = 15 ,& ! ODR array index of the current observation
   lenobs = 15    ! length of a single report

  INTEGER (KIND=iintegers)                , PRIVATE :: &
   maxobs      ,& ! maximum number of observations to be selected for
                  ! the analysis of surface parameters. maxobs=f(maxsgo)
   maxdep         ! maximum number of increments within a model row
                  !   = MAX (NINT(maxobs / (0.4_wp*je_tot)), 20)
                  ! ( = maxobs / (je_tot **0.75) )

!===============================================================================

CONTAINS

!===============================================================================
!+ for controlling the analysis of any near surface parameter
!-------------------------------------------------------------------------------

SUBROUTINE organize_sfcana

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure in "src_sfcana" is the driving routine of the
!   analysis of screen-level temperature, humidity and precipitation
!
! Method:
!   1. Array space required for the surface analysis is allocated and
!      and the arrays are initialized.
!   2. Observations are extracted from ODR and exchanged between all PEs.     
!   3. The grid point analysis is called as a subroutine
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

! Subroutine arguments: none
! --------------------


! Local parameters: none
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    icomm            , & ! MPI communicator
    nanday           , & ! day
    i, j, js, iob    , & ! loop indices 
          jgp, jgpo  , &
    jsta             , &
    keys             , & ! number of sort keys
    igp_ob, jgp_ob   , & ! model grid point assigned to an abservation
    nobpro           , & ! number of observations at local process
    noball           , & ! total number of observations available
    noblast          , & ! last timestep, when observations are available
    iobtyp, icdtyp   , & ! observation/ code type
    ipass            , & ! active/passive status of an observation
    nbspr            , & ! ODR index for precipitation amount [mm]
    iobrtr           , & ! time range of precip observation
    istat                ! status variable


  REAL    (KIND=wp   )     ::  &
    hhan             , & ! hour
    acthr            , & ! actual forecast hour
    z2, ztgranul     , & ! to compute correct reference time
    ztdiff           , & ! time difference
    zhdif                ! difference between station altitude and model height
    

  LOGICAL                  ::  &
    lt2now           , & ! .true. if t2m analysis is to be done at the actual
                         !                                         time step
    lh2now           , & ! .true. if rh2m analysis is to be done 
    lprnow           , & ! .true. if precip. analysis is to be done 
    lffnow               ! .true. if wind speed analysis is to be done 

  LOGICAL, SAVE            ::  &
    lopen  = .TRUE.      ! .true. if surface analysis GRIB and diagnostic files
                         !        to be opened

  CHARACTER (LEN=11)     ::                           &
    sdate            , & ! date in the form dd.mm.yyyy
    stime                ! time in the form hh:min

  CHARACTER (LEN=28)     ::                           &
    yandat2    ! analysis date in the form   wd   dd.mm.yy  hh mm ss UTC

  CHARACTER (LEN=8)      :: date
  CHARACTER (LEN=10)     :: time
  CHARACTER (LEN=16)     :: yroutine
  CHARACTER (LEN=75)     :: yerrmsg


! Local (automatic) arrays:
! -------------------------

  REAL    (KIND=wp   )      , ALLOCATABLE :: &
    rbuf     (:,:)     ! buffer for 'gather_all'

  INTEGER (KIND = iintegers), ALLOCATABLE :: &
    indxo      (:),  & ! index array
    bufloc   (:,:),  & ! buffer array containing observations at local PE
    idat     (:,:)     ! intermediate obs array to be sorted
 
  CHARACTER (LEN = ilstid)  , ALLOCATABLE :: &
    yidloc   (:,:)     ! buffer of station ids of all observations from local PE

  CHARACTER (LEN = 256)                   :: &
    ydatname           ! directory + filename for lansfc-file
!
!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine organize_sfcana
!-------------------------------------------------------------------------------

!----------------------------------------------------------------------------
!  Section 1.0: Initialization
!----------------------------------------------------------------------------

  yroutine  = 'organize_sfcana'

  IF (.NOT.lsurfa)  RETURN

! Initialize the variables for doing the next analysis
! This is done at the very first model step
  IF (ntstep == nstart) THEN
    ht2next  = ht2a
    nt2next  = NINT (ht2next  * 3600.0_wp / dt)
    hh2next  = hh2a
    nh2next  = NINT (hh2next  * 3600.0_wp / dt)
    hffnext  = hffa
    nffnext  = NINT (hffnext  * 3600.0_wp / dt)
    hprcnext = hprc
    nprcnext = NINT (hprcnext * 3600.0_wp / dt)
  ENDIF

! actual hour of the forecast; analysis date and time (reference time)

  ! When using a flexible dt, the hour has eventually to be updated
  ! to the nearest output step (which should be a multiple of 0.25 h = 900.0 s
  ztgranul = 900.0_wp
  z2 = (REAL(ntstep, wp) * dt) / ztgranul
  IF (ABS(REAL(NINT(z2),wp) - z2) > 1E-5_wp) THEN
    ! determine date again with time step ztgranul and number of steps
    ! necessary (z2) to reach the same forecast time that wie have now.
    CALL get_utc_date(NINT(z2), ydate_ini, ztgranul, itype_calendar, yandat1, &
                      yandat2, nanday, hhan)
  ELSE
    CALL get_utc_date(ntstep, ydate_ini, dt, itype_calendar, yandat1,    &
                      yandat2, nanday, hhan)
  ENDIF

  acthr = ntstep*dtdeh

! Set switches for analyses of near surface parameters (t2m,rh2m,prec).
  lt2now = .FALSE.
  lh2now = .FALSE.
  lprnow = .FALSE.
  lffnow = .FALSE.
  IF ( lt2m)        THEN
    !IF (ABS(acthr-ht2a) < epsy .OR.  (acthr > ht2a .AND. MOD (acthr-ht2a,ht2i) < epsy))
    IF (ntstep == nt2next) THEN
      lt2now = .TRUE.
      ht2next = ht2next + ht2i
      nt2next = NINT (ht2next * 3600.0_wp / dt)
    ENDIF
  ENDIF
  IF (lrh2m)        THEN
    !IF (ABS(acthr-hh2a) < epsy .OR. (acthr > hh2a .AND. MOD (acthr-hh2a,hh2i) < epsy))
    IF (ntstep == nh2next) THEN
      lh2now = .TRUE.
      hh2next = hh2next + hh2i
      nh2next = NINT (hh2next * 3600.0_wp / dt)
    ENDIF
  ENDIF
  IF (lprecp)       THEN
    !IF (ABS(acthr-hprc) < epsy )
    IF (ntstep == nprcnext) THEN
      lprnow = .TRUE.
      ! no next is computed: but this was also not the case before changing the control (US)???
    ENDIF
  ENDIF
  IF (lff10m)       THEN
    !IF (ABS(acthr-hffa) < epsy .OR.  (acthr > hffa .AND. MOD (acthr-hffa,hffi) < epsy))
    IF (ntstep == nffnext) THEN
      lffnow = .TRUE.
      hffnext = hffnext + hffi
      nffnext = NINT (hffnext * 3600.0_wp / dt)
    ENDIF
  ENDIF
! last timestep, when observations are available
  noblast = MIN( INT( MAX( nverend , nudgend ) ,iintegers) , nstop )

! set index for precip entry according to 'raintp'
                                   nbspr = nbsr12
  IF (ABS(raintp- 1._wp) <= epsy)  nbspr = nbsrr1
  IF (ABS(raintp- 3._wp) <= epsy)  nbspr = nbsrr3
  IF (ABS(raintp- 6._wp) <= epsy)  nbspr = nbsrr6
  IF (ABS(raintp-24._wp) <= epsy)  nbspr = nbsr24

! invoke condition for doing any surface analysis at current timestep
IF ((lt2now .OR. lh2now .OR. lprnow .OR. lffnow) .AND. (ntstep <= noblast)) THEN

  IF (lopen) THEN

!   open GRIB output file for surface analysis fields
    IF (my_cart_id == 0) THEN
      ydatname = TRIM(pp_lansfc%ydir)//'lansfc'
      CALL open_file( isfco, ydatname, 'w  ', pp_lansfc%yform_write, 0, 0, 1, &
                      .FALSE., idbg_level, yerrmsg, istat )
      IF (istat /= 0) THEN
        yerrmsg = 'opening file lansfc failed ('//pp_lansfc%yform_write//')'
        CALL model_abort (my_cart_id, 2022, yerrmsg, yroutine)
      ENDIF
      PRINT *," Surface analysis GRIB output file <lansfc> opened"
    ENDIF

! open ASCII output file for surface analysis diagnostics
    IF (lwonl) THEN
      OPEN (nsfc ,FILE='YUSURF',FORM='FORMATTED',IOSTAT=istat)
      IF (istat /= 0) THEN
        yerrmsg = 'OPENING OF FILE YUSURF FAILED'
        CALL model_abort (my_cart_id, 2032, yerrmsg, yroutine)
      ENDIF
    ENDIF

    lopen = .FALSE. 
  ENDIF

  icomm   =  icomm_cart

! Maximum allowed time difference between obs time and analysis time
  ztdiff  = MAX (tmxdif(nt2m),tmxdif(nh2m),tmxdif(nprc),tmxdif(nff10m))

! maximum number of observations used (depends on 'ztdiff' resp. 'tmxdif')
  maxobs  = maxsgo * INT((1.0_wp + ztdiff))

! Maximum number of increments within a model row
  maxdep  = MAX ( NINT( maxobs / (0.4_wp*je_tot)) , 20)
! maxdep  = NINT( maxobs / SQRT(SQRT(REAL(je_tot*je_tot*je_tot, wp))) )

  ALLOCATE (bufloc (lenobs,maxobs), STAT=istat)
  ALLOCATE (yidloc (1     ,maxobs), STAT=istat)
  ALLOCATE (bufall (lenobs,maxobs), STAT=istat)
  ALLOCATE (yidall (       maxobs), STAT=istat)
  ALLOCATE (noatab (je_tot),        STAT=istat)
  ALLOCATE (alat   (je_tot),        STAT=istat)
  ALLOCATE (alon   (ie_tot),        STAT=istat)

  bufloc = imdi
  bufall = imdi

  DO j = 1,je_tot
    alat (j) = startlat_tot + (j-1)*dlat
  ENDDO
  DO i = 1,ie_tot
    alon (i) = startlon_tot + (i-1)*dlon

    ! The longitude values have to be limited to the range (-180.0,+180.0)
    IF (alon(i) > 180.0_wp) THEN
      alon(i) = alon(i) - 360.0_wp
    ENDIF
  ENDDO

! ! Set lfd, lds and lbm
! lds = ie_tot * je_tot
! lbm = 1875
! lfd = lds * 2 / iwlength + 2000   ! the "2" means 2 bytes per word

  ! Allocate GRIB arrays
  ALLOCATE (iblock(lfd),  ibmap(lbm),  STAT=istat)
  ALLOCATE (ds_grib(lds), dsup(ndsup), STAT=istat)
  ALLOCATE (ymessage(lfa),             STAT=istat)

!-------------------------------------------------------------------------------
! Section 2.0: Print run description
!-------------------------------------------------------------------------------

  IF (lwonl  .OR. (num_compute == 1)) THEN

    WRITE(nsfc,'(40X," SURFACE (ANY 2-D) ANALYSIS RUN",/,                    &
               & 41X," ------------------------------")')

    IF(lt2now)           THEN
       WRITE(nsfc,'("   -2M TEMPERATURE ANALYSIS ")')
    ENDIF
    IF(lh2now)            THEN
       WRITE(nsfc,'("   -2M RELATIVE HUMIDITY ANALYSIS ")')
    ENDIF
    IF(lprnow)            THEN
       WRITE(nsfc,'("   -PRECIPITATION ANALYSIS ")')
    ENDIF
    IF(lffnow)            THEN
       WRITE(nsfc,'("   -10M WIND SPEED ANALYSIS ")')
    ENDIF

    WRITE(nsfc,'(1x,/," ANALYSIS DATE: ",a,/,1x)') yandat2

!   Execution date and time
    CALL DATE_AND_TIME (date, time)
    sdate = date(7:8)//"."//date(5:6)//". "//date(1:4)
    stime = time(1:2)//":"//time(3:4)
    WRITE(nsfc,'(" Execution commenced on ",a," at ",a5,/,1x)') sdate, stime

!   Precipitation Analysis

    IF(lprecp .AND. ABS(acthr-hprc) < epsy) THEN
     WRITE(nsfc,'(1x,/," Total precipitation amount in the last ",           &
               & F3.0," hours")') raintp
     WRITE(nsfc,'(     " ------------------------------------------------")')
     WRITE(nsfc,'(" Precipitation analysis parameters:")')
     WRITE(nsfc,'(" Weights function:   W   = H(R) * V(H)")')
     WRITE(nsfc,'(" where           :  ",                                    &
               &  " H(R)= (RMAX**2-R**2)/(RMAX**2+R**2)")')
     WRITE(nsfc,'("                    ",                                    &
               &  " V(H)= 0.5*(COS(H/HMAX*PI)+0.5")')
     WRITE(nsfc,'("                    ",                                    &
               &  " HMAX= MAX(MOD_ORO , ",F7.0,") m")')  rmxhds
     WRITE(nsfc,'(" Analysis function:  A = SUM(W*D)/SUM(W)")')
     WRITE(nsfc,'(1x,/," Number of scans: ",I5)') nprcsc
     DO js = 1 , nprcsc
     WRITE(nsfc,'(" Scan",I2,", radius:",F9.0," m")') js,rmprcs(js)
     ENDDO
     WRITE(nsfc,'(" No smoothing")')
    ENDIF

!   2m Temperature Analysis

    IF(lt2m .AND. ABS(acthr-ht2a) < epsy) THEN
     WRITE(nsfc,'(1x,/," 2m temperature analysis sucessive correction",      &
               &  " parameters:")')
     WRITE(nsfc,'(      "--------------------------------------------",      &
               &  "------------")')
     WRITE(nsfc,'(" Weights function:   W   = H(R) * V(H)")')
     WRITE(nsfc,'(" where           :  ",                                    &
               &  " H(R)= (RMAX**2-R**2)/(RMAX**2+R**2)")')
     WRITE(nsfc,'("                    ",                                    &
               &  " V(H)= (HMAX**2-H**2)/(HMAX**2+H**2)")')
     WRITE(nsfc,'("                    ",                                    &
               &  " HMAX= MAX(MOD_ORO/2.7 , ",F7.0,") m")')  rmxhdt

     WRITE(nsfc,'(" Increment function: I = SUM(W*D)/SUM(W)")')
     WRITE(nsfc,'(1x,/," Number of scans:",I5)') nt2msc
     DO  js = 1 , nt2msc
     WRITE(nsfc,'(" Scan",I2,", radius:",F9.0," m")') js,rmt2ms(js+3-nt2msc)
     ENDDO
    ENDIF

!  2m Relative Humidity Analysis

    IF(lrh2m .AND. ABS(acthr-hh2a) < epsy) THEN
     WRITE(nsfc,'(1x,/," 2m Rel. Humid. Analysis Sucessive Correction",      &
               &  " parameters:")')
     WRITE(nsfc,'(      "--------------------------------------------",      &
               &  "------------")')
     WRITE(nsfc,'(" Weights function:   W   = H(R) * V(H)")')
     WRITE(nsfc,'(" where           :  ",                                    &
               &  " H(R)= (RMAX**2-R**2)/(RMAX**2+R**2)")')
     WRITE(nsfc,'("                    ",                                    &
               &  " V(H)= (HMAX**2-H**2)/(HMAX**2+H**2)")')
     WRITE(nsfc,'("                    ",                                    &
               &  " HMAX= MAX(MOD_ORO/2.7 , ",F7.0,") m")')  rmxhdt
     WRITE(nsfc,'(" Increment function: I = SUM(W*D)/SUM(W)")')
     WRITE(nsfc,'(1x,/," Number of scans:",I5)') nr2msc
     DO  js = 1 , nr2msc
     WRITE(nsfc,'(" Scan",I2,", radius:",F9.0," m")') js,rmr2ms(js+3-nr2msc)
     ENDDO
    ENDIF

!   10m Wind Speed Analysis

    IF(lff10m .AND. ABS(acthr-hffa) < epsy) THEN
     WRITE(nsfc,'(1x,/," 10m wind speed analysis sucessive correction",      &
               &  " parameters:")')
     WRITE(nsfc,'(      "--------------------------------------------",      &
               &  "------------")')
     WRITE(nsfc,'(" Weights function:   W   = H(R) * V(H)")')
     WRITE(nsfc,'(" where           :  ",                                    &
               &  " H(R)= (RMAX**2-R**2)/(RMAX**2+R**2)")')
     WRITE(nsfc,'("                    ",                                    &
               &  " V(H)= (HMAX**2-H**2)/(HMAX**2+H**2)")')
     WRITE(nsfc,'("                    ",                                    &
               &  " HMAX= MAX(MOD_ORO/3.0 , ",F7.0,") m")')  rmxhdv

     WRITE(nsfc,'(" Increment function: I = SUM(W*D)/SUM(W)")')
     WRITE(nsfc,'(1x,/," Number of scans:",I5)') nfffsc
     DO  js = 1 , nfffsc
     WRITE(nsfc,'(" Scan",I2,", radius:",F9.0," m")') js,rmfffs(js+3-nfffsc)
     ENDDO
    ENDIF
  ENDIF

!----------------------------------------------------------------------------
!  Section 3.0: Extraction of SYNOP observations from ODR
!----------------------------------------------------------------------------

  nobpro = 0
  Syn_rep: DO  iob = 1,ntotsg

  IF (nobpro >= maxobs) THEN
     PRINT  '(" CAUTION !!!!!  Subroutine ",a,/,                              &
             &" Array size for surface observations too small",               &
             &"   maxobs=",i5)', yroutine, maxobs
     EXIT Syn_rep
  ENDIF

  iobtyp = mosghd (iob,nhobtp)
  icdtyp = mosghd (iob,nhcode)
  ipass  = mosghd (iob,nhpass)

! check active/passive status and observation type
  IF (ipass == 0 .AND. (iobtyp == nsynop .OR. iobtyp == ndribu)) THEN

!   check observation time
    IF (ABS( osghed(iob,nhtime)-acthr ) <= ztdiff+epsy) THEN

      iobrtr = IBITS( mosgbd(iob,nbswwe), nrtrbp, nrtroc )
      IF ((iobrtr >= 5) .AND. (iobrtr <= 7)) THEN
        iobrtr = iobrtr - 4
      ELSE
        iobrtr = iobrtr * 6
      ENDIF

!     check if observation contains usable data
      IF (     (osgbdy(iob,nbst  ) /= rmdi)                                    &
          .OR. (osgbdy(iob,nbsrh ) /= rmdi)                                    &
          .OR. (osgbdy(iob,nbsu  ) /= rmdi .AND. osgbdy(iob,nbsv  ) /= rmdi)   &
          .OR.((osgbdy(iob,nbspr ) /= rmdi) .AND. (ABS( REAL (iobrtr, wp)      &
                                                       -raintp) <= epsy))) THEN

!       check height difference between station altitude and model orography
        igp_ob = mosghd (iob,nhio)
        jgp_ob = mosghd (iob,nhjo)
        zhdif  = osghed (iob,nhalt) - hsurf (igp_ob,jgp_ob)
        IF (zhdif > MAX (170.0_wp,osghed (iob,nhalt)/5.0_wp))   THEN
!         PRINT '(" height difference too big, observation rejected",/,        &
!           & " STAT_ID=",a," OBS TIME=",f5.1,"  ALTITUDE=",f5.0," Mod_Oro=",  &
!           & f6.0,"  LON=",f5.2,"  LAT=",f5.2)',                              &
!             yosghd(iob), osghed (iob,nhtime), osghed (iob,nhalt ),           &
!             hsurf (igp_ob,jgp_ob), osghed (iob,nhilon), osghed (iob,nhjlat)
          CYCLE Syn_rep
        ENDIF

        nobpro = nobpro + 1

!       station id
        yidloc (1,nobpro) = yosghd(iob)

!       model box in the total domain: jgp_tot*10000+igp_tot
        bufloc (nboxnr,nobpro) = mosghd (iob,nhjtot)*10000 + mosghd (iob,nhitot)

!       observation longitude and latitude
        bufloc (noblon,nobpro) = NINT (osghed (iob,nhilon)*100._wp)
        bufloc (noblat,nobpro) = NINT (osghed (iob,nhjlat)*100._wp)

!       station altitude
        bufloc (nstalt,nobpro) = NINT (osghed (iob,nhalt))

!       type of assigned gridpoint (land/sea indicator)
        IF (fr_land  (igp_ob,jgp_ob) >= 0.5_wp) THEN
          bufloc (nindsl,nobpro) = 1
        ELSE
          bufloc (nindsl,nobpro) = 0
        ENDIF

!       screen level temperature and observation error
!       adapt observed temperature to model orography height
        IF (osgbdy (iob, nbst  ) /= rmdi) THEN
          bufloc (nt2mob,nobpro) = NINT ( (osgbdy (iob, nbst  )+              &
                                   0.006_wp*zhdif)*100._wp)
          bufloc (nt2oer,nobpro) = NINT (0.8_wp *                             &
                                   (1._wp + 0.003_wp*ABS(zhdif))*100._wp)
        ELSE
          bufloc (nt2mob,nobpro) = imdi
          bufloc (nt2oer,nobpro) = imdi
        ENDIF

!       screen level humidity and observation error
        IF (osgbdy (iob, nbsrh ) /= rmdi) THEN
          bufloc (nh2mob,nobpro) = NINT (osgbdy (iob, nbsrh )*100._wp)
          bufloc (nh2oer,nobpro) = NINT (0.08_wp *                            &
                                   (1._wp + 0.003_wp*ABS(zhdif))*100._wp)
        ELSE
          bufloc (nh2mob,nobpro) = imdi
          bufloc (nh2oer,nobpro) = imdi
        ENDIF

!       precipitation amount and observation error
        IF ((osgbdy (iob,nbspr ) /= rmdi) .AND. (ABS( REAL (iobrtr, wp)       &
                                                     -raintp) <= epsy)) THEN
          bufloc (nprcob,nobpro) = NINT (osgbdy (iob, nbspr )*100._wp)
          bufloc (nprcer,nobpro) = NINT (1.0_wp * 100._wp) 
        ELSE
          bufloc (nprcob,nobpro) = imdi
          bufloc (nprcer,nobpro) = imdi
        ENDIF

!       10m wind speed and observation error
        IF (osgbdy (iob, nbsu  ) /= rmdi .AND. osgbdy (iob, nbsv  ) /= rmdi)   &
          THEN
          bufloc (nf10ob,nobpro) = NINT (                                      &
           SQRT (osgbdy (iob, nbsu)**2 + osgbdy (iob, nbsv)**2 )*100._wp)
          bufloc (nf10er,nobpro) = NINT (1.4_wp*                               &
                                   (1._wp + 0.003_wp*ABS(zhdif))*100._wp)
        ELSE
          bufloc (nf10ob,nobpro) = imdi
          bufloc (nf10er,nobpro) = imdi
        ENDIF

!       identification of the PE the observation belongs to
        bufloc (nidpe,nobpro) = my_cart_id
!       ODR array index of the current observation
        bufloc (nobnum,nobpro)   = iob
      ENDIF
    ENDIF
  ENDIF
  ENDDO Syn_rep
! CALL comm_barrier (icomm_cart, istat, yerrmsg)

!-------------------------------------------------------------------------------
!- Section 4: Gathering the observation buffer from & for all nodes
!-------------------------------------------------------------------------------

  nobtot  = nobpro

  IF (num_compute > 1)  THEN
    ALLOCATE (rbuf (1,maxobs), STAT=istat)
    rbuf (1,:) = 0._wp

    CALL gather_all ( lenobs,    0,      1, maxobs, ilstid                     &
                    , bufloc, rbuf, yidloc, nobtot, noball, yerrmsg, istat )
!   ===============

    IF (istat /= 0) CALL model_abort (my_cart_id, 13011, yerrmsg,yroutine,istat)
!                   ================
    IF (noball > nobtot) PRINT  '(" CAUTION !!!!!  Subroutine ",A , /          &
                                &," Array size for SYNOP observations too"     &
                                &," small: noball=",I5,", maxobs=",I5)' ,      &
                                yroutine, noball, maxobs       
    IF (lwonl) WRITE( nsfc,'(1X,/," Total number of surface observations"      &
                           &," extracted from ODR:  nobtot=",I5)')  nobtot
    DEALLOCATE (rbuf, STAT=istat)
  ENDIF

!-------------------------------------------------------------------------------
!   Section 5. Sort observation array on station id and model boxes
!------------------------------------------------------------------------------


  IF (nobtot > 0) THEN
!   allocate auxiliary array to be sorted
!   no. of sort keys: ilstid for station id, 1 for box number

!!  ALLOCATE (idat  (nobtot,ilstid+1), STAT=istat )
    ALLOCATE (idat  (nobtot,       1), STAT=istat )
    ALLOCATE (indxo (nobtot)         , STAT=istat )

!   fill sorting array with station id characters and boxnumber,
!   highest key last

!   23.06.2005: sort on station ID cancelled

!!  DO icc = 1,ilstid
!!  idat (1:nobtot,icc) = ICHAR( yidloc(1,1:nobtot) (icc:icc) )
!!  ENDDO
!!  idat (1:nobtot,ilstid+1) = bufloc(nboxnr,1:nobtot) 
    idat (1:nobtot,       1) = bufloc(nboxnr,1:nobtot)

!   sort on station id and boxnumber, most important key last

!!  keys=ilstid+1
    keys=       1
    call ia_orders(idat(1,1),nobtot,nobtot,keys,indxo)

    DO iob = 1,nobtot
      bufall (1:lenobs,iob) = bufloc(1:lenobs,indxo(iob))
      yidall (         iob) = yidloc(1       ,indxo(iob))
    ENDDO
    DEALLOCATE (indxo, STAT=istat)
    DEALLOCATE (idat , STAT=istat)
  ENDIF

! Deallocate the sending buffers
  DEALLOCATE (bufloc, STAT=istat)
  DEALLOCATE (yidloc, STAT=istat)

! Built up a table containing pointers to the first observation within
! each model row
! If a row contains no data the pointer to the next observation is filled in.
! After the last obs, the table entries are set to the total length of
! the observation array + 1.

  jgpo = 0
  jsta = 1
  DO iob = 1,nobtot
    jgp = bufall (nboxnr,iob)/10000
    IF (jgp /= jgpo) THEN
      DO j = jsta,jgp
        noatab (j) = iob
      ENDDO
      jsta = jgp + 1
      jgpo = jgp
    ENDIF
  ENDDO
  DO j = jsta , je_tot
    noatab (j) = nobtot + 1
  ENDDO


!-----------------------------------------------------------------------------
! Section 6: Surface analysis space allocation
!-----------------------------------------------------------------------------

  ALLOCATE (anal   (ie_tot, je_tot), STAT=istat)
  ALLOCATE (firstg (ie_tot, je_tot), STAT=istat)
  ALLOCATE (rsurin (ie, je),         STAT=istat)
  ALLOCATE (indexv (ie, je),         STAT=istat)
  ALLOCATE (wa     (maxobs),         STAT=istat)
  ALLOCATE (wwa    (maxobs),         STAT=istat)
  ALLOCATE (deprow (maxdep,je_tot),  STAT=istat)
  ALLOCATE (rpalat (maxdep,je_tot),  STAT=istat)
  ALLOCATE (rpalon (maxdep,je_tot),  STAT=istat)
  ALLOCATE (rcslon (maxdep,je_tot),  STAT=istat)
  ALLOCATE (rcslat (maxdep,je_tot),  STAT=istat)
  ALLOCATE (rsnlon (maxdep,je_tot),  STAT=istat)
  ALLOCATE (rsnlat (maxdep,je_tot),  STAT=istat)
  ALLOCATE (rpahgt (maxdep,je_tot),  STAT=istat)
  ALLOCATE (rpaobe (maxdep,je_tot),  STAT=istat)
  ALLOCATE (rpalto (maxdep,je_tot),  STAT=istat)
  ALLOCATE (rpalno (maxdep,je_tot),  STAT=istat)
  ALLOCATE (ypasid (maxdep,je_tot),  STAT=istat)
  ALLOCATE (npagpt (maxdep,je_tot),  STAT=istat)
  ALLOCATE (nuprow (je_tot),         STAT=istat)
  ALLOCATE (ndnrow (je_tot),         STAT=istat)
  ALLOCATE (nnodep (je_tot),         STAT=istat)

!-----------------------------------------------------------------------------
! Section 7: Start analysis of surface parameters
!-----------------------------------------------------------------------------

  IF (lt2now) CALL sfc_grpeva ('T_2M      ')
!             ===============

  IF (lh2now) CALL sfc_grpeva ('RELHUM_2M ')
!             ===============

  IF (lprnow) CALL sfc_grpeva ('TOT_PREC  ')
!             ===============

  IF (lffnow) CALL sfc_grpeva ('SP_10M    ')
!             ===============

!-----------------------------------------------------------------------------
! Section 8: Memory de-allocation
!-----------------------------------------------------------------------------

  DEALLOCATE (anal   , STAT=istat)
  DEALLOCATE (firstg , STAT=istat)
  DEALLOCATE (rsurin , STAT=istat)
  DEALLOCATE (indexv , STAT=istat)
  DEALLOCATE (wa     , STAT=istat)
  DEALLOCATE (wwa    , STAT=istat)
  DEALLOCATE (deprow , STAT=istat)
  DEALLOCATE (rpalat , STAT=istat)
  DEALLOCATE (rpalon , STAT=istat)
  DEALLOCATE (rcslon , STAT=istat)
  DEALLOCATE (rcslat , STAT=istat)
  DEALLOCATE (rsnlon , STAT=istat)
  DEALLOCATE (rsnlat , STAT=istat)
  DEALLOCATE (rpahgt , STAT=istat)
  DEALLOCATE (rpaobe , STAT=istat)
  DEALLOCATE (rpalto , STAT=istat)
  DEALLOCATE (rpalno , STAT=istat)
  DEALLOCATE (ypasid , STAT=istat)
  DEALLOCATE (npagpt , STAT=istat)
  DEALLOCATE (nuprow , STAT=istat)
  DEALLOCATE (ndnrow , STAT=istat)
  DEALLOCATE (nnodep , STAT=istat)

  DEALLOCATE (bufall , STAT=istat)
  DEALLOCATE (yidall , STAT=istat)
  DEALLOCATE (noatab , STAT=istat)
  DEALLOCATE (alat   , STAT=istat)
  DEALLOCATE (alon   , STAT=istat)

  ! Deallocate arrays for IO
  DEALLOCATE (iblock, ymessage, ibmap, ds_grib, dsup, STAT=istat)

ENDIF  ! doing any surface analysis at current timestep

  IF ((ntstep == noblast) .AND. (.NOT. lopen)) THEN

!   close GRIB output file for surface analysis fields
    IF (my_cart_id == 0) THEN
      CALL close_file( isfco, pp_lansfc%yform_write, 0, 0, 1, &
                       .FALSE., idbg_level, yerrmsg, istat )
      IF (istat /= 0) THEN
        yerrmsg = 'closing GRIB file lansfc failed'
        CALL model_abort (my_cart_id, 2023, yerrmsg, yroutine)
      ENDIF
      PRINT *," Surface analysis GRIB output file <lansfc> closed"
    ENDIF

! open ASCII output file for surface analysis diagnostics
    IF (lwonl) THEN
      CLOSE (nsfc ,IOSTAT=istat)
      IF (istat /= 0) THEN
        yerrmsg = 'CLOSING OF FILE YUSURF FAILED'
        CALL model_abort (my_cart_id, 2033, yerrmsg, yroutine)
      ENDIF
    ENDIF

  ENDIF

!-------------------------------------------------------------------------------
!  End of the Subroutine
!-------------------------------------------------------------------------------

END SUBROUTINE organize_sfcana

!===============================================================================
!===============================================================================
!+ for initializing the output control structure for sfcana
!-------------------------------------------------------------------------------

SUBROUTINE init_sfcana

!-------------------------------------------------------------------------------
!
! Description:
!  This subroutine initializes the special output group for the variables of
!  the near-surface analysis.  In case of grib_api output, the samples, which
!  are read in init_output, are cloned to the variable gribapi_id of this
!  output structure.
!
!-------------------------------------------------------------------------------

INTEGER(KIND=iintegers) :: n, n1, iz1, iz2, iz3, izstat

!-------------------------------------------------------------------------------

  izstat = 0_iintegers
  ALLOCATE(pp_lansfc, STAT = izstat)

  ! Initialize necessary components of pp_lansfc
  pp_lansfc%yvarml(:) = '          '    ! just a list of 4 variables
  pp_lansfc%yvarpl(:) = '          '    ! not used
  pp_lansfc%yvarzl(:)        = ''       ! not used
  pp_lansfc%yvarsl(:)        = ''       ! not used
  pp_lansfc%yvarc (:)        = ''       ! not used
  pp_lansfc%ilist_ml(:,:)    =  0       ! will be set later
  pp_lansfc%ilist_pl(:,:)    =  0       ! not used
  pp_lansfc%ilist_zl(:,:)    = -1       ! not used
  pp_lansfc%ilist_sl(:,:)    = -1       ! not used
  pp_lansfc%ilist_c (:,:)    = -1       ! not used
  pp_lansfc%nyvar_m          =  0       ! will be set later
  pp_lansfc%nyvar_p          =  0       ! not used
  pp_lansfc%nyvar_z          =  0       ! not used
  pp_lansfc%nyvar_s          =  0       ! not used
  pp_lansfc%nyvar_c          =  0       ! not used

  AllOCATE (pp_lansfc%ngrib(1), STAT=izstat)
  pp_lansfc%ngrib(:)         =  0       ! not used
  pp_lansfc%outsteps         =  0       ! not used
  pp_lansfc%nexthour         =  0       ! not used
  pp_lansfc%nextstep         =  0       ! not used
  pp_lansfc%lhour            = .TRUE.
  pp_lansfc%nprocess_ini_out = -999999
  pp_lansfc%nprocess_bd_out  = -999999
  pp_lansfc%nunit_of_time    = 1        ! must that be reset?
  pp_lansfc%slon             = startlon_tot
  pp_lansfc%slat             = startlat_tot
  pp_lansfc%elon             = endlon_tot
  pp_lansfc%elat             = endlat_tot
  pp_lansfc%i_out_start      = 1
  pp_lansfc%j_out_start      = 1
  pp_lansfc%i_out_end        = ie_tot
  pp_lansfc%j_out_end        = je_tot
  pp_lansfc%ie_out_tot       = ie_tot
  pp_lansfc%je_out_tot       = je_tot
  pp_lansfc%yform_write      = yform_lansfc
  pp_lansfc%ydir             = ydir_lansfc
  pp_lansfc%ysuffix          = ' '      ! not used
  pp_lansfc%ytunit           = ' '      ! not used
  pp_lansfc%ydomain          = 'f'      ! not used
  pp_lansfc%nrbit            = 16       ! do we need a namelist variable for that????
  pp_lansfc%plev(:)          = -1.0_wp  ! not used
  pp_lansfc%zlev(:)          = -1.0_wp  ! not used
  pp_lansfc%kepin            =    0     ! not used
  pp_lansfc%kezin            =    0     ! not used
  pp_lansfc%lcheck           = .TRUE.   ! not used
  pp_lansfc%lwrite_const     = .FALSE.  ! not used
  pp_lansfc%luvmasspoint     = .FALSE.  ! not used
  pp_lansfc%lanalysis        = .FALSE.
  pp_lansfc%lsfc_ana         = .TRUE.
  pp_lansfc%l_p_filter       = .FALSE.  ! not used
  pp_lansfc%l_z_filter       = .FALSE.  ! not used
  pp_lansfc%l_pmsl_filter    = .FALSE.  ! not used
  pp_lansfc%l_fi_filter      = .FALSE.  ! not used
  pp_lansfc%l_fi_pmsl_smooth = .FALSE.  ! not used

  NULLIFY (pp_lansfc%next)

  pp_lansfc%nyvar_m   =  4
  pp_lansfc%yvarml(1) = 'T_2M'
  pp_lansfc%yvarml(2) = 'RELHUM_2M'
  pp_lansfc%yvarml(3) = 'SP_10M'
  pp_lansfc%yvarml(4) = 'TOT_PREC'

  n1 = 0_iintegers
  loop_list_m: DO n = 1, pp_lansfc%nyvar_m
    DO iz3 = 1, num_gribtabs
      DO iz2 = 0, 255
        DO iz1 = 1,4
          IF (TRIM(var(iz1,iz2,iz3)%name) ==  TRIM(pp_lansfc%yvarml(n)) ) THEN
            n1 = n1+1
            pp_lansfc%ilist_ml(1,n1) = iz1
            pp_lansfc%ilist_ml(2,n1) = iz2
            pp_lansfc%ilist_ml(3,n1) = iz3

            CYCLE loop_list_m
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    ! If this point is reached, no variable with name pp_lansfc%yvarml(n) was found
    IF (my_cart_id == 0) THEN
      PRINT *, 'Variable ', TRIM(pp_lansfc%yvarml(n)), ' was not found in the GRIB table'
    ENDIF
  ENDDO loop_list_m

#ifdef GRIBAPI
  IF     (pp_lansfc%yform_write == 'api1') THEN
    ! clone this sample to igrib1_sample 
    CALL grib_clone(igrib1_id, pp_lansfc%igribapi_id, izstat)
    IF (izstat /= GRIB_SUCCESS) THEN
      PRINT *,   ' *** Error in grib_clone: from sample 1 ', izstat
    ENDIF
  ELSEIF (pp_lansfc%yform_write == 'api2') THEN
    CALL grib_clone(igrib2_id, pp_lansfc%igribapi_id, izstat)
    IF (izstat /= GRIB_SUCCESS) THEN
      PRINT *,   ' *** Error in grib_clone: from sample 2 ', izstat
    ENDIF
  ENDIF
#endif

  CALL make_grib_init(pp_lansfc)

END SUBROUTINE init_sfcana

!===============================================================================
!===============================================================================
!+ for performing the grid point analysis of any near-surface parameter
!-------------------------------------------------------------------------------

SUBROUTINE sfc_grpeva ( yvar )

!-------------------------------------------------------------------------------
!
! Description:
!
!   In case of screen level temperature or humidity a successive correction
!   scheme is adopted. Here the deviations of the observations from a guess  
!   field  are interpolated to each grid point using all data within a
!   given radius of influence.
!   The analysis method for precipitation amount is a simple distance-weighted
!   average of observations within a relative small area around the grid
!   points. If there are no observations the search radius is somewhat enlarged
!   and if there are still no observations no analyses is done at the g.p.
!
! Method:
!   The outer loop of the grid point analysis controlls the number
!   of passes with different data selection radii. The first step of
!   each cycle is the exchange of the background field between all PEs
!   in order to have global guess fields. The next loop runs
!   over model rows. For each model row the observations within a
!   latitude band are selected and the deviations from the guess field
!   are calculated. For each grid point the weighted mean of the increments
!   or absolute values are calculated. The increment field is added
!   to the background field. Then the next cycle is performed but now a
!   smaller radius of influence is used . After the last scan the increment
!   field is smoothed out. Finally some statistics on the analysis increments
!   are carried out and the analysed field is written out.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

! Subroutine arguments:
! --------------------
  CHARACTER (LEN=10), INTENT(IN)        :: &
   yvar                    ! analysis variable

! Local parameters: none
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
   izlocstat, irm      , & ! status variable
   implcode            , & ! status variable
   nscan               , & ! number of successive correction scans
   iobval              , & ! index of the observed variable 
   ioberr              , & ! index of the observation error
   imaxvi              , & ! index of the maximum departure
   jscan               , & ! scan number loop index
   nopass              , & ! number of passes over field to be smoothed
   nosmot              , & ! number of smoothings within a pass
   jpass, jsmot        , & ! loop indices
   jup1, jup2          , & ! loop indices
   jdw1, jdw2          , & ! loop indices
   jrow                , & ! loop indices
   jlat, jlon          , & ! loop indices
   i, j                , & ! loop indices
   jlb, jdep               ! loop indices
  INTEGER (KIND=iintegers) ::  &
   jlat_tot, jlon_tot  , & ! model lat/lon in the total domain
   jsta, jend          , & ! start and end of DO loop
   ista, iend          , &
   kzdims(24)          , & ! vertical dimensions of fields to be exchanged
   ino                 , & ! index of an obs. within the nudging obs array 
   iob, iobid          , & ! index of an obs. within the surface ana obs array
   igp, jgp            , & ! model g.p. assigned to obs.
   imbox               , & ! model box; jgp*10000+igp
   igptyp              , & ! sea/land type of assigned grid point
   idep                , & ! actual obs. number in a row
   islgp               , & ! sea/land identification
   inodat              , & ! observation counter
   isum                , & ! counter
   irowmx, irowmn      , & ! model row of maximum and minimum change
   ipoimx, ipoimn      , & ! model column of maximum and minimum change
   izerror, istat      , & !
   ntri                    ! time range indicator

  REAL    (KIND=wp   )     ::  &
   fpvsw                ,& ! statement function; Magnus formula
   zt                   ,& ! temperature in K
   zmxhdf              , & ! maximum height difference between obs and g.p.
   zhdfac              , & ! height dependent scaling factor for vertical
                           ! range of influence
   ztdiff              , & ! maximum time difference
   rdumpf              , & ! amplification factor
   zmscan              , & ! scan radius in m
   zdisdg              , & ! scan radius in earth radii
   zlat, zlon          , & ! model latitude, longitude
   zlonwb, zloneb      , & ! western/ eastern boundary of observation area
   zstalt              , & ! station altitude
   zobval              , & ! observed value
   zoberr              , & ! observation error
   zfgerr              , & ! first guess error
   zdep                , & ! observed departure from guess field
   zrsmrs              , & ! square of radius of influence
   zcoslt, zsinlt      , & ! cos and sin of model lat.
   zcosln, zsinln      , & ! cos and sin of model lon.
   zmodor                  ! model height
  REAL    (KIND=wp   )     ::  &
   zdisop              , & ! auxiliary variable
   zdistm              , & ! distance between obs and g.p. in m
   zdissq              , & ! square of distance
   zrormx              , & ! ABS(zdistm)/zmscan
   zvdmax              , & ! maximum vertical displacement
   zhdmhd              , & ! square of zvdmax
   zvdis, zavdis       , & ! actual vertical displacement
   zverd2              , & ! square of zvdis
   zvdovx              , & ! zavdis/zvdmax
   w, ww               , & ! sum of weights; sum of weighted increments
   znu                 , & ! amplification factor +-rdumpf
   zw1, zw2            , & ! weights for smoothing
   zinc                , & ! analysis increment
   zsum                , & ! sum of increments
   zlonmx, zlatmx      , & ! model lon and lat of maximum increment
   zlonmn, zlatmn      , & ! model lon and lat of minimum increment
   zaverg, zmaxch, zsd , & ! Average increment, standard deviation and max. inc.
   zmaxvl, zminvl          ! maximum and minimum increment


  LOGICAL                  ::  &
   lmyobs              , & ! .true. if obs. belongs to the local domain
   lseaob              , & ! .true. if observation platform on sea
   lrej                    ! .true. if observation is rejected

  CHARACTER (LEN=30) ::  ylab      ! lab for diagnostic print out
  CHARACTER (LEN= 8) ::  yunits    ! unit of printed values
  CHARACTER (LEN=16) ::  yvarib    ! type of variable
  CHARACTER (LEN=75) ::  yerrmsg   ! error message
  CHARACTER (LEN=25) ::  yroutine  ! name of subroutine


! Local (automatic) arrays:
! -------------------------

  INTEGER (KIND=iintegers), ALLOCATABLE       :: &
   indxi_(:), indxj_(:), indxk_(:)

  REAL (KIND=wp),     ALLOCATABLE       :: &
   zana (:,:)          , & ! local analysis field
   fgss (:,:)          , & ! model first guess field
   zarr (:,:,:)        , & ! array containing the subarrays from all PEs
   zarr1(:,:)          , & ! intermediate array to be gathered
   zauxar (:)              ! auxillary array

! variables for reorganisation of loops, to enhance efficiency on vector proces.
   LOGICAL                  :: lprecip_
   INTEGER (KIND=iintegers) :: i_, j_, ii_, jdep2_
   REAL (KIND=wp),     DIMENSION(10000) ::  zdistm_, zrormx_, wwa_, zmodor_    &
                                          , rpahgt_, zvdis_, zvdmax_, zwi_
!
!------------ End of header ----------------------------------------------------

! Statement functions
!--------------------

! Magnus formula for water 
  fpvsw (zt)       = b1 * EXP( b2w*(zt-b3)/(zt-b4w) ) 
   
!! Statement functions
!! ichbpt(invar,ibp,ibo) = IAND (invar, &
!!                     IOR (ishft (2**32-1,ibp+ibo), ishft(2**32-1,ibp-32) )
 
!-------------------------------------------------------------------------------
! Begin Subroutine sfc_grpeva
!-------------------------------------------------------------------------------

  kzdims(:) = 0_iintegers

! yvar is not a vectorizable type !
  lprecip_ = (yvar == 'TOT_PREC')

  yroutine = 'sfc_grpeva'
  irm      = 0
! Allocate  arrays
  ALLOCATE ( zana  (ie,je)                    , STAT=izlocstat )
  ALLOCATE ( fgss  (ie,je)                    , STAT=izlocstat )
  ALLOCATE ( zarr  (ie_max,je_max,num_compute), STAT=izlocstat )
  ALLOCATE ( zarr1 (ie_max,je_max)            , STAT=izlocstat )
  ALLOCATE ( zauxar(MAX(ie,je))               , STAT=izlocstat )
  IF (izlocstat /= 0) THEN
    yerrmsg = ' *** Allocation of local arrays failed ***'
    CALL model_abort (my_cart_id, 15000, yerrmsg, yroutine)
!   ================
  ENDIF

!------------------------------------------------------------------------------
!  Section 1.0: Initialization
!-------------------------------------------------------------------------------

  IF (yvar == 'T_2M')   THEN
     ylab   = '2M TEMPERATURE; UNIT=DEGREES'
     nscan  = nt2msc
     zhdfac = 2.7_wp
     zmxhdf = rmxhdt
     zfgerr = 1.2_wp
     iobval = nt2mob
     ioberr = nt2oer
     imaxvi = 1
     ztdiff = tmxdif(nt2m)
     yunits = 'DEGREE C'
     yvarib = 'T2M-OBSERVATION'
!    smoothing characteristics
     nopass      =   1
     rdumpf      =   0.5_wp
     nosmot      =   2
!    first guess field
     fgss        = t_2m

  ELSEIF (yvar == 'RELHUM_2M')   THEN
     ylab   = '2M REL. HUMIDITY; UNIT=%'
     nscan  = nr2msc
     zhdfac = 2.7_wp
     zmxhdf = rmxhdt
     zfgerr = 0.12_wp
     iobval = nh2mob
     ioberr = nh2oer
     imaxvi = 2
     ztdiff = tmxdif(nh2m)
     yunits = '% RELHUM'
     yvarib = 'RH2M-OBSERVATION'
!    smoothing characteristics
     nopass      =   1
     rdumpf      =   0.5_wp
     nosmot      =   2
!    first guess field
     DO j = 1,je
       DO i = 1,ie
         fgss(i,j) =   fpvsw( td_2m(i,j) ) / fpvsw( t_2m (i,j) )
!        fgss(i,j) =   fpvsw( td_2m(i,j), b1, b2w, b3, b4w )                   &
!                    / fpvsw( t_2m (i,j), b1, b2w, b3, b4w )
       ENDDO
     ENDDO

  ELSEIF (yvar == 'SP_10M')   THEN
     ylab   = '10M WIND SPEED; UNIT=M/S%'
     nscan  = nfffsc
     zhdfac = 3.0_wp
     zmxhdf = rmxhdv
     zfgerr = 2.1_wp
     iobval = nf10ob
     ioberr = nf10er
     imaxvi = 3
     ztdiff = tmxdif(nff10m)
     yunits = 'M/S'
     yvarib = 'W10M-OBSERVATION'
!    smoothing characteristics
     nopass      =   1
     rdumpf      =   0.5_wp
     nosmot      =   2
!    first guess field
     DO j = 1,je
     DO i = 1,ie
      fgss(i,j) = SQRT (u_10m(i,j)**2 + v_10m(i,j)**2 )
     ENDDO
     ENDDO

  ELSEIF (yvar == 'TOT_PREC')   THEN
     ylab   = 'PRECIPITATION AMOUNT; UNIT=mm'
     nscan  = nprcsc
     zmxhdf = rmxhds
     iobval = nprcob
     ioberr = nprcer
     imaxvi = 4
     ztdiff = tmxdif(nprc)
     yunits = 'mm WATER'
     yvarib = 'PRECIP. AMOUNT  '
!    smoothing characteristics
     nopass      =   0
     rdumpf      =   0.5_wp
     nosmot      =   2
!    first guess field
     anal(:,:)   = 0._wp
     zana(:,:)   = 0._wp
     fgss(:,:)   = 0._wp
  ENDIF

  IF (TRIM(yvar) == 'T_2M' .OR. TRIM(yvar) == 'RELHUM_2M' .OR. TRIM(yvar) == 'SP_10M') THEN
    kzdims(1:24)=(/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                  &
       (  1  ,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
        ie, je, kzdims, jstartpar, jendpar, 2, nboundlines, my_cart_neigh, &
        lperi_x, lperi_y, l2dim,                                           &
        20000+nexch_tag, ldatatypes, ncomm_type, izerror, yerrmsg,         &
        fgss )
    zana = fgss
  ENDIF
  indexv = 0
  rsurin = 0._wp

!------------------------------------------------------------------------------
!  Section 2.0:  Start up loop over successive correction cycles
!------------------------------------------------------------------------------

  SucCor: DO jscan = 1,nscan
! ==========================
! ++++++++++++++++++++++++++
  
 
! preset arrays initial values
 
  nuprow = 0
  ndnrow = 0
  nnodep = 0
  npagpt = 0
  rcslon = 0._wp
  rsnlon = 0._wp
  rcslat = 0._wp
  rsnlat = 0._wp
  rpalat = 0._wp
  rpalon = 0._wp
  rpalto = 0._wp
  rpalno = 0._wp
  rpahgt = 0._wp
  rpaobe = 0._wp
  deprow = 0._wp
  wa     = 0._wp
  wwa    = 0._wp
  ypasid = ' '
  IF (yvar /= 'TOT_PREC') THEN
    rsurin(:,:) = 0._wp
  ENDIF

!------------------------------------------------------------------------------
! Section 2.1  Compose total background fields
!-------------------------------------------------------------------------------

  IF (yvar /= 'TOT_PREC')                THEN
    IF (num_compute > 1)  THEN
!     Gather the data
      zarr1(1:ie,1:je) = zana (1:ie,1:je)

!     Gather the data

      CALL gather_values ( zarr1, zarr, ie_max, je_max, num_compute, imp_reals &
                         , -1, icomm_cart, yerrmsg, implcode )
!     ==================
      IF (implcode /= 0)                                                       &
        CALL model_abort (my_cart_id, 11151, yerrmsg, yroutine, implcode)

!     Combine the subarrays

      CALL combine_subarrays (zarr , anal)
!     ======================
    ELSE
      anal = zana
    ENDIF
  ENDIF  ! yvar /= 'TOT_PREC'

  IF (jscan == 1) THEN
     firstg = anal
  ENDIF

!  Write out 10m wind speed model forecast in grib code format
  IF (my_cart_id == 0 .AND. yvar == 'SP_10M' .AND. jscan == 1) THEN
!    time range indicator ntri=13 means nudging forecast
     ntri = 13
     CALL sfcout ( yvar, ntri)
!    ===========
  ENDIF


  IF (nobtot == 0)                        EXIT SucCor

!---------------------------------------------------------------------------
! Section 2.2: Find out latitude bands for each model row
!---------------------------------------------------------------------------
 
  IF (yvar == 'T_2M')  THEN
    zmscan    =  rmt2ms(jscan+3-nt2msc)
  ELSE IF(yvar == 'RELHUM_2M') THEN
    zmscan    =  rmr2ms(jscan+3-nr2msc)
  ELSE IF (yvar == 'TOT_PREC') THEN
    zmscan    = rmprcs(jscan)
  ELSE IF(yvar == 'SP_10M') THEN
    zmscan    =  rmfffs(jscan+3-nfffsc)
  ENDIF
  zdisdg      =  zmscan / (r_earth*degrad)
 
!  upper limit
 
   DO  jup1   =  je_tot,  1,  -1
     zlat     =  MIN (alat(je_tot),alat(jup1)+zdisdg)
     DO  jup2 =  je_tot-1,1,  -1
       IF (alat(jup2) <  zlat .AND. alat(jup2+1) >= zlat) THEN
         nuprow(jup1) = jup2
       EXIT
       ENDIF
     ENDDO
   ENDDO

!  lower limit
 
   DO  jdw1   =     1,   je_tot
   zlat       =  MAX(alat(1),alat(jdw1)-zdisdg)
     DO  jdw2 =     2,   je_tot
       IF (alat(jdw2) > zlat .AND. alat(jdw2-1) <= zlat) THEN
         ndnrow(jdw1) = jdw2 - 1 
         EXIT
       ENDIF
     ENDDO
   ENDDO

!  Find western and eastern boarder of observation area for this PE
   zlonwb     =  MAX (startlon_tot,startlon-zdisdg)
   zloneb     =  MIN (startlon_tot+(ie_tot-1)*dlon,                           &
                      startlon    +(ie    -1)*dlon + zdisdg)

!---------------------------------------------------------------------------
! Section 2.3:  Select observations and calculate departures  
!---------------------------------------------------------------------------

  jsta  = ndnrow (isubpos (my_cart_id,2)-nboundlines)
  jend  = nuprow (isubpos (my_cart_id,4)+nboundlines)
  ista  = MAX (INT ((zlonwb-startlon_tot)/dlon) + 1 , 1)
  iend  = MIN (INT ((zloneb-startlon_tot)/dlon) + 2 , ie_tot)


! Loop over data rows

  Data_rows: DO  jrow  =  jsta , jend    
! ++++++++++++++++++++++++++++++++++++++
   iob  =  noatab(jrow)
 
!  Loop over observations in a row

   All_obs: DO WHILE ((jrow <  je_tot .AND. iob <  noatab(jrow+1)) .OR.  &
                      (jrow == je_tot .AND. iob <= nobtot))
!  ==================

   lrej = .false.
!  assigned model grid point
   imbox  = bufall (nboxnr,iob)
   jgp  = imbox/10000
   igp  = imbox - jgp*10000

   IF (igp < ista .OR. igp > iend) THEN
      iob = iob + 1
      CYCLE All_obs
   ENDIF

!  find out if current obs belongs to this PE. Get ODR index and
!  type of observation 
   ino      = -1
   IF (bufall(nidpe,iob) == my_cart_id) THEN
!    index of an obs. within the nudging obs array
     ino    = bufall(nobnum,iob)
     lmyobs = .true.
     lseaob = (BTEST( mosghd(ino,nhschr), nvsebp ))
   ELSE
     lmyobs = .FALSE.
   ENDIF 

!  station identification
   iobid  = iob

!  sea/land type of assigned grid point
   igptyp = bufall(nindsl,iob)

!  station altitude
   zstalt = REAL  (bufall(nstalt,iob), wp)

!  observed value
   IF (bufall(iobval,iob) == imdi .OR.                                       &
       bufall(ioberr,iob) == imdi)       THEN
      iob = iob + 1
      CYCLE All_obs
   ELSE
     zobval = bufall(iobval,iob) *0.01_wp

!    observation error
     zoberr = bufall(ioberr,iob) *0.01_wp
   ENDIF

!  Check departure from guess field
   zdep   = zobval - anal(igp,jgp)

   IF (ABS(zdep) > rmaxdp(imaxvi))                THEN
     IF (jscan == nscan  .AND. lmyobs)                                      &
       PRINT '(a16," FAILED CHECK AGAINST GUESS FIELD",                     &
          & "  DEPARTURE=",e9.2, a8," ZOBVAL=",e9.2/,                       &
          & " STAT_ID=",a," HEIGHT=",f5.0,"  LON=",f5.2,                    &
          & " LAT=",f5.2,"  OBS TIME=",f5.1)' ,                             &
              yvarib, zdep,yunits, zobval, yidall(iobid), osghed(ino,nhalt),&
              osghed (ino,nhilon), osghed (ino,nhjlat), osghed(ino,nhtime)   
     lrej = .true. 
   ENDIF

   IF (.NOT.lrej) THEN
!  -------------------

!    Update number of data in a row
     IF (nnodep(jrow) >= maxdep)   THEN
       print  '(" C A U T I O N     Insufficient space for surface ",       &
                & "observations; Subroutine ",a)', yroutine
       print *,'maxdep, jrow, nnodep(jrow): ', maxdep,jrow,nnodep(jrow)
       DO i=1,nnodep(jrow)
       print *,'no, stid,lat,lon,depart:', i, ypasid(i,jrow), &
                & rpalto(i,jrow), rpalno(i,jrow), deprow(i,jrow)
       ENDDO
       EXIT  All_obs
     ENDIF
     nnodep(jrow)     = nnodep(jrow) + 1
     idep             = nnodep(jrow)
 
!    Station characteristics
     ypasid(idep,jrow) = yidall(iobid)
     npagpt(idep,jrow) = igptyp
     rpalat(idep,jrow) = startlat_tot + (jgp-1)*dlat
     rpalon(idep,jrow) = startlon_tot + (igp-1)*dlon

     ! The longitude values have to be limited to the range (-180.0,+180.0)
     IF (rpalon(idep,jrow) > 180.0_wp) THEN
       rpalon(idep,jrow) = rpalon(idep,jrow) - 360.0_wp
     ENDIF

     rpalto(idep,jrow) = bufall(noblat,iob) *0.01_wp
     rpalno(idep,jrow) = bufall(noblon,iob) *0.01_wp
     rpahgt(idep,jrow) = REAL  (bufall(nstalt,iob), wp)
     rpaobe(idep,jrow) = zoberr
     rcslon(idep,jrow) = COS (rpalon(idep,jrow)*degrad)
     rcslat(idep,jrow) = COS (rpalat(idep,jrow)*degrad)
     rsnlon(idep,jrow) = SIN (rpalon(idep,jrow)*degrad)
     rsnlat(idep,jrow) = SIN (rpalat(idep,jrow)*degrad)

!    departures
     deprow(idep,jrow)= zdep

     IF (jscan == nscan  .AND.lmyobs .AND. lseaob .AND. igptyp == 1) THEN
       PRINT '(" OBSERVATION TYPE AND TYPE OF ASSIGNED",                    &
       &" GRID POINT ARE DIFFERENT",/,                                      &
       &" STAT_ID=",A,"  GP TYPE=",I2,"  LSEAOB=",L2)' ,                    &
          yidall(iobid), igptyp, lseaob
     ENDIF
 
!    Print observation for diagnostic purposes
!!   IF (lwonl .AND. ldiasa)                    THEN
!!     WRITE(nsfc,'(" JROW=",I3," IDEP=",I3," STID=",a," IGP=",I3,          &
!!    &" JGP=",I3," LON=",F6.2," LAT=",F6.2," ALT=",F6.0," DEP=",F5.1)')    &
!!     jrow,idep,yidall(iobid),igp,jgp,rpalon(idep,jrow),rpalat(idep,jrow), &
!!     rpahgt(idep,jrow),zdep
!!   ENDIF
 
   ENDIF ! .NOT. lrej  
!  ------------------

!  increase pointer
   iob = iob + 1

  ENDDO All_obs
! =============
 
! Close row loop
 
 ENDDO Data_rows
!+++++++++++++++

!---------------------------------------------------------------------------
! Section 3.0  Successive correction analysis scan
!---------------------------------------------------------------------------
!
!  For all grid points calculate scan increment through a formula:
!
!                                 N
!                                 __
!                                \   W * D
!                                /__  I   I
!                                I=1
!        I(scan)=I(scan-1) +   ------------
!                                  N
!                                  __
!                                 \   W
!                                 /__  I
!                                 I=1
!
!  where I       =  increment (I(0)=0)
!        scan    =  scan number
!        W       =  weight
!        D       =  departure
!    weights can actually be calculated as:
!
!              W = H * V
!
!  where H = horizontal part of weight function [=H(dist)]
!        V = vertical part of weight function [=H(dz,...)]
!
!
! Scan constants
! square of radius
  zrsmrs    =    zmscan * zmscan
 
  i_ = 0
  DO jlb = 1 , je_tot
     i_ = i_ + nnodep(jlb)
  ENDDO
  IF (i_ > maxobs)                                                            &
     PRINT  '(" CAUTION !!!!!  Subroutine ",A,2X,A/,                          &
             &" Array size for surface obs too small",                        &
             &"   maxobs=",2I6)', yroutine, yvar, maxobs, i_
  ALLOCATE (indxi_ (i_ +1),  STAT=istat)
  ALLOCATE (indxj_ (i_ +1),  STAT=istat)
  ALLOCATE (indxk_ (i_ +1),  STAT=istat)
  indxi_ = 0
  indxj_ = 0
  indxk_ = 0

! Start up loop over model latitudes
! ----------------------------------
 
  Model_lat:   DO  jlat  =  1 , je
! ================================

!  row number in total domain
   jlat_tot = isubpos(my_cart_id,2) - nboundlines - 1 + jlat

!  cos and sin of model lat.
   zcoslt   =   COS(alat(jlat_tot)*degrad)
   zsinlt   =   SIN(alat(jlat_tot)*degrad)
 
!  Start up loop over model lon. points
!  ------------------------------------
 
   Model_lon:  DO  jlon  =  1 , ie
!  ===============================

!   colum number in total domain
    jlon_tot = isubpos(my_cart_id,1) - nboundlines - 1 + jlon

    zcosln=   COS(alon(jlon_tot)*degrad)
    zsinln=   SIN(alon(jlon_tot)*degrad)
    zmodor=   hsurf(jlon,jlat) 
    IF (fr_land (jlon,jlat) < 0.5_wp) THEN
       islgp = 0
    ELSE
       islgp = 1
    ENDIF
 
!   check index array in case of precipitation analysis
    IF (yvar == 'TOT_PREC' .AND. indexv(jlon,jlat) == 1)     CYCLE Model_lon
 
!--------------------------------------------------------------------------
!  Section  3.1: Observations' distance, weights and contributions
!                 to the current grid poinT
!--------------------------------------------------------------------------
 
!-------------------------------------------------------------------------
!  SECTION  3.1.1  Horizontal part
!-------------------------------------------------------------------------
!        Weighting functions used are of a form:
!
!                 2      2
!                rmax - r
!              W=--------        (analysis of t2m r2m ff10m)
!                 2      2
!                rmax + r
!
!         and  W=(0.5*(COS(r/rmax*pi)+1)/(1+r/rmax)   (precipitation)
!
!     where W       =  weight
!           rmax    =  max. radius distance
!           r       =  actual distance
!     The expression used to calculate the distance is of the form:
!
!              D=(COS(LON )*COS(LON )*COS(LAT )*COS(LAT ) +
!                        I         J         I         J
!                 SIN(LON )*COS(LAT )*SIN(LON )*COS(LAT ) +
!                        I         I         J         J
!                 SIN(LAT )*SIN(LAT ))                    *
!                        I         J
!                 REARTH
!
!     where D         =    distance
!           LON /LAT  = lon./lat. of anal. point
!              I    I
!           LON /LAT  = lon./lat. of obs. point
!              J    J
!           REARTH    = earth's radius
!
!----------------------------------------------------------------------
!   pre-set number of data to 0
    inodat       = 0

! determine indices of influencing observations
! (this allows 1-d loops over influencing observations subsequently)
    i_ = 0
    DO jlb = ndnrow(jlat_tot), nuprow(jlat_tot)
       DO jdep = 1, nnodep(jlb)
          i_ = i_ + 1
          indxi_(i_) = jdep
          indxj_(i_) = jlb
          indxk_(i_) = jdep+inodat
       ENDDO
       inodat = inodat + nnodep(jlb)
    ENDDO

!   Loop over influencing rows

    j_ = 0
!CDIR NODEP
    Rows_inf:  DO ii_ = 1, i_

!     Loop over data within row
!     Calculate horizontal distance in m to all
!     observations for given latitude row

        jdep   = indxi_(ii_)
        jlb    = indxj_(ii_)
        jdep2_ = indxk_(ii_)

        zdisop    = zcosln*rcslon(jdep,jlb)*zcoslt*rcslat(jdep,jlb) +      &
                    zsinln*zcoslt*rsnlon(jdep,jlb)*rcslat(jdep,jlb) +      &
                    zsinlt*rsnlat(jdep,jlb)
        zdistm    = MAX(ACOS(MIN(zdisop,1._wp)),1.e-5_wp) * r_earth
        zdissq    = zdistm * zdistm
        zrormx    = ABS(zdistm)/zmscan

!       Calculate weights
        IF (zdistm <= zmscan) THEN
          IF (lprecip_) THEN
            wwa (jdep2_) = .5_wp*(COS(zrormx*pi)+                 &
                                 1._wp)/(1._wp+0.3_wp*zrormx)
          ELSE
            wwa (jdep2_) = (zrsmrs-zdissq) /(zrsmrs+zdissq)
          ENDIF
        ELSE
            wwa (jdep2_) = 0._wp
        ENDIF

        IF (lwonl .AND. ldiasa .AND. (jlat == jonl) .AND. (jlon == ionl)       &
                  .AND. (wwa(jdep2_) > 0.001_wp))   THEN
           j_ = j_ + 1
           zdistm_(j_) = zdistm
           zrormx_(j_) = zrormx
           wwa_(j_)    = wwa(jdep2_) 
        ENDIF

!   Close outer loop over influencing rows

    ENDDO Rows_inf

!   separate loop for non-vectorizable 'write'
    DO ii_ = 1, j_
       WRITE( nsfc,'(" zdistm=",f7.0," zrormx=",f5.3," wwa=",f5.3)')     &
              zdistm_(ii_),zrormx_(ii_),wwa_(ii_)
    ENDDO
 
!-------------------------------------------------------------------------
!  SECTION  3.1.2: Vertical part as a function of height difference
!-------------------------------------------------------------------------

!        Weighting function used is of a form:
!
!                 2      2
!                zmax - z
!              W=--------       (analysis of t2m r2m, ff10m)
!                 2      2
!                zmax + z
!
!         and  W=(0.5*(COS(z/zmax*pi)+1)/(1+0.8*z/zmax)   (precipitation)
!
!     where W       =  weight
!           zmax    =  max height distance
!           z       =  actual height distance
! -----------------------------------------------------------------------

!   Pre-set
    inodat       = 0
 
!   Square of max. height difference
    zhdmhd       = zmxhdf * zmxhdf

!   Loop over influancing rows

    IF (yvar == 'TOT_PREC') THEN
       j_ = 0
       zvdmax   = MAX (zmodor,zmxhdf)

!   Loop over influancing rows
!   Loop over data within row
!CDIR NODEP
       Rows_inf_2_prec: DO ii_ = 1, i_

          jdep   = indxi_(ii_)
          jlb    = indxj_(ii_)
          jdep2_ = indxk_(ii_)

!         Calculate vertical distance in m to all
!         observations for given latitude row
          zvdis       =  rpahgt(jdep,jlb) - zmodor
          zavdis      =  ABS (zvdis )
          zvdovx      =  zavdis/zvdmax

!         Calculate weights
          IF (zavdis <= zvdmax) THEN
            wwa (jdep2_) = 0.5_wp*(COS(zvdovx*pi)+1._wp)/     &
                          (1._wp+0.8_wp*zvdovx)* wwa(jdep2_)
          ELSE
            wwa (jdep2_) = 0._wp
          ENDIF


          IF (lwonl .AND. ldiasa .AND. (jlat == jonl) .AND. (jlon == ionl)     &
                    .AND. (wwa(jdep2_) > 0.001_wp))   THEN
             j_ = j_ + 1
             zwi_(j_)    = 0.5_wp *COS(zavdis/zvdmax*pi) + 0.5_wp
             zmodor_(j_) = zmodor
             rpahgt_(j_) = rpahgt(jdep,jlb)
             zvdis_(j_)  = zvdis
             zvdmax_(j_) = zvdmax
             wwa_(j_)    = wwa(jdep2_)
          ENDIF
       ENDDO Rows_inf_2_prec

!      separate loop for non-vectorizable 'write'
       DO ii_ = 1, j_
          WRITE(nsfc,'(" zmodor=",f5.0," rpaght=",f5.0," zvdis=",f5.0,    &
      &      " zvdmax=",f5.0," zwi=",f5.3," wwa=",f5.3)')                 &
      &      zmodor_(ii_),rpahgt_(ii_),zvdis_(ii_),zvdmax_(ii_),          &
      &      zwi_(ii_),wwa_(ii_)
       ENDDO

    ELSE

       zvdmax  = MAX(zmodor/zhdfac,zmxhdf)
       zhdmhd  = zvdmax*zvdmax

!CDIR NODEP
       Rows_inf_2_noprec: DO ii_ = 1, i_

          jdep   = indxi_(ii_)
          jlb    = indxj_(ii_)
          jdep2_ = indxk_(ii_)

!         Calculate vertical distance in m to all
!         observations for given latitude row
          zvdis       =  rpahgt(jdep,jlb) - zmodor
          zavdis      =  ABS (zvdis )
          zverd2      =  zvdis * zvdis

!         Calculate weights
          IF (zavdis <= zvdmax) THEN
             wwa (jdep2_) = (zhdmhd-zverd2) /(zhdmhd+zverd2) *wwa(jdep2_)
          ELSE
             wwa (jdep2_) = 0._wp
          ENDIF

       ENDDO Rows_inf_2_noprec
    ENDIF

!---------------------------------------------------------------------------
!  Section 3.2:   Contributions of each datum, number of data and
!                 list of zero departure. Sum of contributions and weights
!---------------------------------------------------------------------------

   inodat    =   0
   w         =   0._wp
   ww        =   0._wp

   IF (.NOT. lprecip_) THEN
!CDIR NODEP
      DO ii_ = 1, i_
         jdep   = indxi_(ii_)
         jlb    = indxj_(ii_)
         jdep2_ = indxk_(ii_)
!        Use land obs for land gridpoints and water obs for water gp only
!        Take into account observation errors
         wwa(jdep2_) = REAL(MIN( ABS(1-islgp-npagpt(jdep,jlb)) , 1),wp) *         &
             wwa (jdep2_) * zfgerr**2/rpaobe(jdep,jlb)**2
         wa(jdep2_) = wwa(jdep2_) * deprow(jdep,jlb)
!        Sum up contributions and weights
         w  = w  + wa(jdep2_)
         ww = ww + wwa(jdep2_)
      ENDDO
   ELSE ! lprecip_
!CDIR NODEP
      DO ii_ = 1, i_
         jdep   = indxi_(ii_)
         jlb    = indxj_(ii_)
         jdep2_ = indxk_(ii_)
         wa(jdep2_) = wwa(jdep2_) * deprow(jdep,jlb)
!        Sum up contributions and weights
         w  = w  + wa(jdep2_)
         ww = ww + wwa(jdep2_)
      ENDDO
   ENDIF

!----------------------------------------------------------------------------
!  Section 3.3:   Evaluate increments and set index array
!----------------------------------------------------------------------------
 
   IF (ww   <= 0._wp)                                CYCLE Model_lon

   IF (yvar == 'TOT_PREC') THEN
     IF (ww >  0.20_wp)                    THEN
        rsurin(jlon,jlat) = w/ww
        indexv(jlon,jlat) = 1
     ELSE IF(jscan == nprcsc .AND. ww >  0.02_wp)      THEN
        rsurin(jlon,jlat) = w/ww
        indexv(jlon,jlat) = 1
     ELSE
        rsurin(jlon,jlat) = 0._wp
        CYCLE Model_lon
     ENDIF
     IF (rsurin(jlon,jlat) <  0.04_wp) rsurin(jlon,jlat) = 0._wp

   ELSEIF (yvar /= 'TOT_PREC') THEN
     rsurin(jlon,jlat) = w/(ww+1._wp)
     indexv(jlon,jlat) = 1
   ENDIF  ! yvar

!----------------------------------------------------------------------------
!  Section 3.4: Print out diagnostics
!----------------------------------------------------------------------------
  
   IF (ldiasa .AND. lwonl)                                       THEN
     IF(jlat == jonl  .AND. jlon == ionl )                 THEN
       zlat = rlat(jlon,jlat)/degrad
       zlon = rlon(jlon,jlat)/degrad
       WRITE(nsfc,'(1x,/," Diagnostic ",a14," analysis   scan ",i1,/,         &
         & "  grid point  ix=",i3,"  iy=",i3,"  lat=",f6.2,"  lon=",f7.2,     &
         & "  height=",f5.0)')                                                &
           ylab(1:14),jscan,jlon_tot,jlat_tot,zlat,zlon,zmodor

       WRITE(nsfc,'(" Observations influencing this grid point:",/,/,         &
         & " stat_id     lat      lon   height  increment  weight")')
       inodat = 0

!CDIR NODEP
       DO ii_ = 1, i_
          jdep   = indxi_(ii_)
          jlb    = indxj_(ii_)
          jdep2_ = indxk_(ii_)
          IF (wwa(jdep2_) >= 0.001_wp) THEN
             WRITE(nsfc,'(1X,a8,3x,f6.2,2x,f7.2,3x,f5.0,4x,f5.2,4x,f4.2)')     &
                   ypasid(jdep,jlb),rpalto(jdep,jlb),rpalno(jdep,jlb)          &
                  ,rpahgt(jdep,jlb),deprow(jdep,jlb),wwa(jdep2_)
          ENDIF
       ENDDO

       WRITE(nsfc,'(" Resulting ",a3," analysis increment: ",f6.2)')           &
             yvar, rsurin(jlon,jlat)
     ENDIF  ! selected g.p.
   ENDIF  ! ldiasa

!----------------------------------------------------------------------------
!  Section 3.5:  Add up increments to background field
!----------------------------------------------------------------------------
 
  IF (yvar /= 'TOT_PREC')        THEN
    zana(jlon,jlat) = anal(jlon_tot,jlat_tot)+rsurin(jlon,jlat)
  ENDIF
! ---------------------------------------------------------------
 
! Close longitude loop

  ENDDO Model_lon
! ===============
 
! CLose latitude loop
 
  ENDDO Model_lat
! ===============
 
  DEALLOCATE (indxi_ , STAT=istat)
  DEALLOCATE (indxj_ , STAT=istat)
  DEALLOCATE (indxk_ , STAT=istat)

! Close the scan loop
 
  ENDDO SucCor
! ++++++++++++
 
!-----------------------------------------------------------------------------
! Section 4:   Smooth final increments
!-----------------------------------------------------------------------------

  IF (yvar /= 'TOT_PREC') THEN
    rsurin(:,:) = zana(:,:) - fgss(:,:)
 
!   Loop over no. of passes and smoothings
    znu         =   rdumpf
    DO      jpass  =    1   ,    nopass
    DO      jsmot  =    1   ,    nosmot
!   -----------------------------------

      zw1         =   1.0_wp - znu
      zw2         =   0.5_wp * znu
 
!     Smoothing over lon. points
 
!       Loop over lat. rows
      DO    jlat   =    1   ,    je
 
!       Loop over lon. points
        DO    jlon   =    2   ,    ie-1
         zauxar(jlon)       = rsurin(jlon-1,jlat) + rsurin(jlon+1,jlat)
        ENDDO
        DO    jlon   =    2   ,    ie-1
         rsurin(jlon,jlat)  = zw1*rsurin(jlon,jlat) + zw2*zauxar(jlon)
        ENDDO
      ENDDO  ! jlat
 
!     Smoothing over lat. points
 
!       Loop over lon. points
      DO   jlon   =    1   ,    ie
 
!       Loop over lat. points
        DO   jlat   =    2   ,    je-1
         zauxar(jlat)       = rsurin(jlon,jlat-1) + rsurin(jlon,jlat+1)
        ENDDO
        DO   jlat   =    2   ,    je-1
         rsurin(jlon,jlat)  = zw1*rsurin(jlon,jlat) + zw2*zauxar(jlat)
        ENDDO
      ENDDO  ! jlon
 
!   Close the loops
      znu      =    -rdumpf
    ENDDO  ! jsmot
    ENDDO  ! jpass

  ENDIF  ! yvar /= TOT_PREC

!-----------------------------------------------------------------------------
! Section 6:  Final analysis
!-----------------------------------------------------------------------------

! Add up smoothed increments to background field
  zana(:,:) = fgss(:,:) + rsurin(:,:)

  IF (num_compute > 1) THEN
!   Store arrays to be gathered in buffer zarr1
    zarr1(1:ie,1:je) = zana(1:ie,1:je)

!   Gather the data

    CALL gather_values ( zarr1, zarr, ie_max, je_max, num_compute, imp_reals &
                       , -1, icomm_cart, yerrmsg, implcode )
!   ==================
    IF (implcode /= 0)                                                       &
      CALL model_abort (my_cart_id, 11152, yerrmsg, yroutine, implcode)

!   Combine the subarrays in the correct order on PE 0

    CALL combine_subarrays (zarr, anal)
!   ======================
  ELSE
    anal = zana
  ENDIF

! De-allocate local allocatable arrays
  DEALLOCATE ( zana   , STAT=izlocstat )
  DEALLOCATE ( fgss   , STAT=izlocstat )
  DEALLOCATE ( zarr   , STAT=izlocstat )
  DEALLOCATE ( zarr1  , STAT=izlocstat )
  DEALLOCATE ( zauxar , STAT=izlocstat )

!-----------------------------------------------------------------------------
! Section 6: Statistics
!-----------------------------------------------------------------------------

!IF (my_cart_id == 0) THEN
 IF (lwonl)           THEN

! Preset statistics variables
  isum     =  0
  zsum     =  0._wp
  zaverg   =  0._wp
  zsd      =  0._wp
  zmaxch   =  0._wp
  zmaxvl   =  -999999._wp
  irowmx   =  0
  ipoimx   =  0
  zminvl   =  +999999._wp
  irowmn   =  0
  ipoimn   =  0
 
! Loop over model latitudes
  DO   jlat  =  1  , je_tot

!  Row statistics calculation
  DO   jlon  = 1 ,  ie_tot
 
!  Number of updated points, sum of changes and standard deviation
    isum     =  isum + 1
    zinc     =  anal(jlon,jlat)-firstg(jlon,jlat)
    zsum     =  zsum  + zinc
    zsd      =  zsd   + zinc * zinc
 
!  max. positve change and lat/lon index

    IF(zinc >  zmaxvl)                 THEN
      zmaxvl   =  zinc
      irowmx   =  jlat
      ipoimx   =  jlon
    ENDIF
 
!  max. negative change and lat/lon index
 
    IF(zinc <  zminvl)                 THEN
      zminvl   =  zinc
      irowmn   =  jlat
      ipoimn   =  jlon
    ENDIF
  ENDDO  ! jlon
  ENDDO  ! jlat
   
!-----------------------------------------------------------------------
 
! Average change, standard deviation and max. change
 
  zaverg   =  zsum / isum
  zsd      =  sqrt((zsd-2*zsum*zaverg+zaverg*zaverg*isum)/isum)
  zmaxch   =  MAX (ABS(zmaxvl),ABS(zminvl))

! Print out statistics
 
  zlonmx     = alon  (ipoimx)
  zlatmx     = alat  (irowmx)
 
  zlonmn     = alon  (ipoimn)
  zlatmn     = alat  (irowmn)
 
  WRITE(nsfc,'(1x,/," *** STATISTICS ON ANALYSIS OF ",a)')  ylab   
  WRITE(nsfc,'("      SUM OF CHANGES       = " ,f11.4)') zsum
  WRITE(nsfc,'("      NO. OF ANAL. POINTS  = " ,i11  )') isum
  WRITE(nsfc,'("      AVERAGRE CHANGE      = " ,f11.4)') zaverg
  WRITE(nsfc,'("      STANDARD DEVIATION   = " ,f11.4)') zsd
  WRITE(nsfc,'("      MAX. CHANGE          = " ,f11.4)') zmaxch
  WRITE(nsfc,'("      MAX. POSITIVE CHANGE = " ,f11.4)') zmaxvl
  WRITE(nsfc,'("      AT MODEL LAT./LON.   = " ,2(f11.4,1x))')  zlatmx,zlonmx
  WRITE(nsfc,'("      MAX. NEGATIVE CHANGE = " ,f11.4)') zminvl
  WRITE(nsfc,'("      AT MODEL LAT./LON.   = " ,2(f11.4,1x))')  zlatmn,zlonmn
  ENDIF  ! lwonl==.t.

!-----------------------------------------------------------------------------
! Section 7. Write out analysis in grib code format
!-----------------------------------------------------------------------------

 IF (my_cart_id == 0) THEN
! time range indicator ntri=0 means uninitialzed analysis
  ntri = 0
  CALL sfcout ( yvar, ntri)
! ===========
 
 ENDIF  ! my_cart_id == 0        

!-------------------------------------------------------------------------------
!  End of the Subroutine
!-------------------------------------------------------------------------------

END SUBROUTINE sfc_grpeva

!-------------------------------------------------------------------------------

!+ Module procedure in "SFCANA" for grid pt. analysis of surface-level paramet.
!-------------------------------------------------------------------------------

SUBROUTINE sfcout (yvar, ktr)

!-------------------------------------------------------------------------------
! Description:
! sfcout packs surface analyses into grib format and writes them to disk.
!
! Method:
! This module procedure creates the product definition block and grid
! definition block according to the description of putpd1 and putgd1 (dwd) 
! and the WMO. Then the records are packed and written to disk.
!
!=======================================================================
!
! Declarations:
!
!=======================================================================

!=======================================================================
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
  CHARACTER (LEN=10),       INTENT(IN) :: &
   yvar                    ! analysis variable
  INTEGER (KIND=iintegers), INTENT(IN) :: &
  ktr                      ! time range indicator

!-----------------------------------------------------------------------
!
! Local parameters:

! GRIB output fields
!  name, element num, factor, bias, gribtab, lev.typ, levtop,  levbot

! (yen(i), ngbnr(i), gbfc(i), gbbs(i), ntab(i), nltp(i), nlvt(i), nlvb(i),   &
!                                                             i = 1,iofld) / &
!  't2m'   ,  11   ,   1.0   ,  0.0  ,    2   ,   105  ,   0  ,    2 ,       &
!  'r2m'   ,  52   ,   1.0   ,  0.0  ,    2   ,   105  ,   0  ,    2 ,       &
!  'fff'   ,  32   ,   1.0   ,  0.0  ,    2   ,   105  ,   0  ,   10 ,       &
!  'prc'   ,  61   ,   1.0   ,  0.0  ,    2   ,     1  ,   0  ,    0 /

! Local scalars:
#ifdef GRIBDWD
  INTEGER (KIND=intgribf), EXTERNAL    :: IREFTS
#endif

  INTEGER  (KIND=iintegers)            :: ivar, iv, i, j, igrl, ierrf, &
                                          i1, i2, i3, izerror, igrbid

  INTEGER (KIND=intgribf)              :: iz_ps=1

  CHARACTER (LEN=75)                   :: yerrmsg
  CHARACTER (LEN=25)                   :: yroutine

  INTEGER  (KIND=iintegers)            :: izntri

  INTEGER (KIND=int_ga)                :: irecord_lga      ! length of grib record

  REAL     (KIND=wp)                   :: zbias, zgribfactor

#ifdef GRIBAPI
INTEGER (KIND=kindOfSize) :: ibyte_size_out
#endif
!   
!
!- End of header
!
!===============================================================================
!-------------------------------------------------------------------------------
!  Section 1: Initializations
!-------------------------------------------------------------------------------

  yroutine    = 'sfcout'
  izerror     = 0_iintegers
  irecord_lga = 0_int_ga

! find table index of output variable
  ivar  = -1
  DO iv = 1, pp_lansfc%nyvar_m
    IF (TRIM(yvar) == TRIM(pp_lansfc%yvarml(iv)) ) THEN  
      ivar = iv
      i1 = pp_lansfc%ilist_ml(1,iv)
      i2 = pp_lansfc%ilist_ml(2,iv)
      i3 = pp_lansfc%ilist_ml(3,iv)
      PRINT *, ' found variable: ', TRIM(yvar), i1, i2, i3
      EXIT
    ENDIF
  ENDDO
  IF (ivar == -1) THEN
    WRITE(nsfc,'(A,A)') ' Undefined output variable:  ', TRIM(yvar)
    yerrmsg=' Undefined output variable'
    CALL model_abort (my_cart_id, 2022, yerrmsg, yroutine)
  ENDIF

!  moving arraydimensions into idims
!  the following settings depend on the grib size and could be different
  idims_out( 7) = lds

!  real dimensions
  idims_out(15) = ie_tot*je_tot
  idims_out(17) = ie_tot*je_tot

!-------------------------------------------------------------------------------
!  Section 2: Set the grid description section
!-------------------------------------------------------------------------------

  IF (pp_lansfc%yform_write == 'grb1') THEN
    ! set the namelist-dependent idims_out values
    idims_out( 7) = pp_lansfc%ie_out_tot * pp_lansfc%je_out_tot
    idims_out(15) = pp_lansfc%ie_out_tot * pp_lansfc%je_out_tot
    idims_out(17) = pp_lansfc%ie_out_tot * pp_lansfc%je_out_tot

    ! gridpoints, simple packing, floating point data
    ibds(2)   = 0

    ! nrbit, number of bits
    ibds(5)   = pp_lansfc%nrbit

    ! no bitmap
    ibms(3)   = -2
  ENDIF

#ifdef GRIBAPI
  IF (pp_lansfc%yform_write(1:3) == 'api') THEN
    CALL grib_clone(pp_lansfc%igribapi_id, igrbid, izerror)
    IF (izerror /= GRIB_SUCCESS) THEN
      PRINT *,   ' *** Error in grib_clone: from pp_lansfc sample ', izerror
    ENDIF
  ENDIF
#endif

  ! number of gridpoints and coordinates of corners
  CALL make_grib_grid (pp_lansfc, igrbid, i1, i2, i3, ' ', izerror)

  IF (izerror /= 0) THEN
    WRITE(nsfc,'(A)') ' Error in setting GRIB meta data: grid definition'
    WRITE( *  ,'(A)') ' Error in setting GRIB meta data: grid definition'
    yerrmsg=' Error in setting GRIB meta data: grid definition'
    CALL model_abort (my_cart_id, 2022, yerrmsg, yroutine)
  ENDIF

!------------------------------------------------------------------------
!  Section 3: Set the product description section
!------------------------------------------------------------------------

  ! set pp_lansfc%lanalysis depending on time range indicator
  IF (ktr == 13) THEN
    pp_lansfc%lanalysis = .TRUE.
    pp_lansfc%lsfc_ana  = .FALSE.
  ELSE
    pp_lansfc%lanalysis = .FALSE.
    pp_lansfc%lsfc_ana  = .TRUE.
    IF (TRIM(yvar) == 'TOT_PREC' .AND. (pp_lansfc%yform_write /= 'api2')) THEN
      ! save original timeRangeIndicator from setup_vartab (which is 4)
      ! and temporarily set it to 0, to get the usual GRIB(1) coding
      ! for this product
      izntri = var(i1,i2,i3)%ntri
      var(i1,i2,i3)%ntri = ktr
    ENDIF
  ENDIF

  CALL make_grib_product (pp_lansfc, igrbid, i1, i2, i3, 1, yandat1,    &
                          0, ' ', 0.0_wp, .FALSE., izerror)

  IF (izerror /= 0) THEN
    WRITE(nsfc,'(A)') ' Error in setting GRIB meta data: product definition'
    WRITE( *  ,'(A)') ' Error in setting GRIB meta data: product definition'
    yerrmsg=' Error in setting GRIB meta data'
    CALL model_abort (my_cart_id, 2022, yerrmsg, yroutine)
  ENDIF

  IF (ktr /= 13 .AND. TRIM(yvar) == 'TOT_PREC' .AND. (pp_lansfc%yform_write /= 'api2')) THEN
    ! reset original ntri again
      var(i1,i2,i3)%ntri = izntri
  ENDIF

!-------------------------------------------------------------------------
!  Section 4.  Scaling and packing the output field
!--------------------------------------------------------------------------

! Change the REAL format for the grib library and scale ds_grib

  zbias       = var(i1,i2,i3)%bias
  IF (TRIM(yvar) == 'RELHUM_2M') THEN
    ! the rest of the model already has the correct unit for relhum,
    ! therefore the scaling factor in setup_vartab is 1.0, but here we need
    zgribfactor = 100.0_wp
  ELSE
    zgribfactor = var(i1,i2,i3)%factor
  ENDIF

  DO j=1,je_tot
    DO i=1,ie_tot
      ds_grib((j-1)*ie_tot+i) =                                                   &
         REAL ( ((anal(i,j)+ zbias) * zgribfactor), irealgrib)
    ENDDO
  ENDDO

  IF     (pp_lansfc%yform_write == 'grb1') THEN

#ifdef GRIBDWD
    CALL grbex1(idwdednr, iz_ps, undefgrib, ndims, idims_out, ipds_out, igds_out, &
                ibms , ibds , ibmap, dsup , ds_grib  , iblock, ierrf)
    IF (ierrf /= 0) THEN
      yerrmsg = 'error in grbex1'
      CALL model_abort (my_cart_id, 2022, yerrmsg, yroutine)
    ENDIF

    ! length of GRIB record in bytes
    igrl = idims_out(19)
#endif

#ifdef GRIBAPI
  ELSEIF (pp_lansfc%yform_write(1:3) == 'api') THEN

    ! Set these values again explicit, because they could be overwritten
    ! by other data
    CALL grib_set (igrbid, 'bitsPerValue',   pp_lansfc%nrbit)

    CALL grib_set (igrbid, 'values',         ds_grib(:), ierrf)
    IF (ierrf /= GRIB_SUCCESS) THEN
      yerrmsg = 'error in lansfc for grib_set: values'
      CALL model_abort (my_cart_id, 2022, yerrmsg, yroutine)
    ENDIF

    ! length of GRIB record in bytes
    CALL grib_get_message_size(igrbid, ibyte_size_out, ierrf)

    IF (ibyte_size_out <= lfa) THEN
      CALL grib_get(igrbid, 'totalLength', irecord_lga)
      CALL grib_copy_message(igrbid,  ymessage)
    ELSE
      yerrmsg = 'error with message length: ymessage too small: '
      CALL model_abort (my_cart_id, 2022, yerrmsg, yroutine)
    ENDIF
#endif

  ENDIF

!--------------------------------------------------------------------------
!  Section 5:  Write GRIB field to disk
!--------------------------------------------------------------------------

  IF     (pp_lansfc%yform_write == 'grb1') THEN
#ifdef GRIBDWD
    CALL cuegex (isfco,iblock,igrl,ierrf)
    IF (ierrf /= 0) THEN
      yerrmsg = 'error in lansfc cuegex'
      CALL model_abort (my_cart_id, 2022, yerrmsg, yroutine)
    ENDIF
    PRINT *, TRIM(yvar), '- DWD-GRIB field written on file'
#endif
  ELSEIF (pp_lansfc%yform_write(1:3) == 'api') THEN
#ifdef GRIBAPI
    CALL grib_write_bytes (isfco, ymessage, irecord_lga, ierrf)
    IF (ierrf /= GRIB_SUCCESS) THEN
      yerrmsg = 'Error in lansfc grib_write_bytes'
      CALL model_abort (my_cart_id, 2022, yerrmsg, yroutine)
    ENDIF
#endif
  ENDIF

!=================================================================

END SUBROUTINE sfcout
!-------------------------------------------------------------------------------

!+ Module procedure in "SFCANA" for grid pt. analysis of surface-level paramet.
!-------------------------------------------------------------------------------

      SUBROUTINE ia_orders ( idat, nx, n, m, index)
 
!=======================================================================
!
! Description:
! Orders a multi-dimensional array
!
! Method:
! Orders a multi-dimensional array by recursive
! call of subroutine ia_qsort
!
!=======================================================================
!
! Current Code Owner: DWD, R. Hess
!    phone: +49-69-8062-2723, fax: +49-69-8062-2723
!    email: Reinhold.Hess@dwd.de
!
!=======================================================================
!
! Subroutine History:
! S-Version  Date       Name
! ---------- ---------- ----
! 1.1        2001/03/24 R. Hess
! 1.2        2001/08/22 G.Paul no type specification
!  Initial Release
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!
!=======================================================================
!
! Subroutine arguments
! --------------------
! Scalar arguments:

  INTEGER nx,n,m
! nx:      dimension of idat (:,m)
! n:       number of elements to be sorted
! m:       length of elements
!
! Array arguments :
! INTEGER(KIND=int_packed) :: idat (nx,m)
  INTEGER                  :: idat (nx,m)
  INTEGER              ::  index (nx)

! index:   integer array of length N containing pointers to the elements

! Local scalar

  INTEGER i
!-----------------------------------------------------------------------
!
!- End of header
!
!
  do i=1,N
    index(i) = i
  end do

  CALL ia_qsort(idat, nx, n, m, 1, n, index )
!
!
  END SUBROUTINE ia_orders
!
!=======================================================================
!
  RECURSIVE SUBROUTINE ia_qsort( idat, nx, n, m, l, r, index )
!
  INTEGER nx, n, m, l, r

!   l: left  boundary of sumdomain
!   r: right boundary of sumdomain

  INTEGER             :: idat (nx,m)
  INTEGER             :: index (nx)
!

!
  INTEGER i, j, v, t
!
  IF (l .lt. r) THEN

!   devide subdomain by central element
!   move it right boundary

     v              = index((l+r)/2)
     index((l+r)/2) = index(r)
     index(r)       = v

     i = l
     j = r - 1

     DO WHILE (.true.)

        DO WHILE (ia_comp(idat,index(i),v,m) .lt. 0)
           i = i + 1
        END DO

        DO WHILE ((ia_comp(idat,index(j),v,m) .ge. 0)  &
                   .and. (j .gt. l))
           j = j - 1
        END DO

        if (i .ge. j) EXIT

        t         = index(i)
        index(i)  = index(j)
        index(j) = t

     END DO

     index(r) = index(i)
     index(i) = v

     CALL ia_qsort(idat, nx, n, m, l,  i-1,  index )
     CALL ia_qsort(idat, nx, n, m, i+1,  r,  index )

  END IF

  END SUBROUTINE ia_qsort
!
!=======================================================================
!
  INTEGER FUNCTION ia_comp(ia,j,k,m)

!    compares two vectors
!    returns -1, if 1st vector is smaller
!    returns  1, if 1st vector is bigger
!    returns  0, if vectors are equal

  INTEGER j,k,m
  INTEGER                  :: ia(:,:)

  INTEGER i

  ia_comp = 0

!    compare highest priority first

  DO i = m, 1, -1
     IF (ia(j,i) .lt. ia(k,i)) THEN
        ia_comp = -1
        EXIT
     END IF
     IF (ia(j,i) .gt. ia(k,i)) THEN
        ia_comp = 1
        EXIT
     END IF

  END DO

  END FUNCTION ia_comp

!===============================================================================

! ELEMENTAL REAL FUNCTION fpvsw  ( zt, b1, b2w, b3, b4w )
  !------------------------------------------------------
! REAL    (KIND=wp)         , INTENT (IN)  ::  zt, b1, b2w, b3, b4w
  !----------------------------------------------------------------------------
  ! Magnus formula for water:  input  'zt'  : temperature
  !                            output 'fpvsw': saturation water vapour pressure
  ! Magnus formula for ice  :  if constants 'b2i', 'b4i' for ice are used for
  !                            'b2w', 'b4w' in the call
  !----------------------------------------------------------------------------
  !
! fpvsw  =  b1 * EXP( b2w*(zt-b3)/(zt-b4w) )
  !
! END FUNCTION fpvsw

!-------------------------------------------------------------------------------

END MODULE src_sfcana
!
!- End of module
