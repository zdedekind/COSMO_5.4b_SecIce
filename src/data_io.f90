!+ Data module for variables concerned with I/O
!------------------------------------------------------------------------------

MODULE data_io

!------------------------------------------------------------------------------
!
! Description:
!  This data module contains all data necessary for input and output of grib
!  files.
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Ulrich Schaettler
!  Initial release
! 1.7        1998/07/16 Guenther Doms
!  Increase of the number 'nvar' of output variables to 90.
! 1.8        1998/08/03 Ulrich Schaettler
!  Got parameters for Grib library from module data_parameters.f90
! 1.9        1998/09/16 Guenther Doms
!  Increase of the number 'nvar' of output variables to 91.
! 1.10       1998/09/29 Ulrich Schaettler
!  New control variable for writing analysis-files in structure pp_nl.
! 1.14       1998/10/26 Ulrich Schaettler
!  Split igds into igds_in and igds_out. Increase nvar to 93.
! 1.17       1998/11/17 Ulrich Schaettler
!  New control variable for writing constant fields in structure pp_nl.
! 1.20       1999/01/07 Guenther Doms
!  nvar has been increased to 96
! 1.24       1999/03/01 Guenther Doms
!  nvar has veen increased to 97
! 1.26       1999/03/25 Guenther Doms
!  nzmxc has been decreased from 14 to 13
! 1.30       1999/06/24 Matthias Raschendorfer
!  nlevels has been increased to 40.
!  nvar has been increased to 102 (TKE, convection closure)
! 1.32       1999/08/24 Guenther Doms
!  nvar has been increased to 105
! 1.33       1999/10/14 Reinhold Hess
!  nvar has been increased to 109
! 1.34       1999/12/10 Ulrich Schaettler
!  GRIB variables lfd, lbm and lds are not declared as PARAMETERs any more.
!  The Kind parameters for GRIB are now declared with defaults
! 1.39       2000/05/03 Ulrich Schaettler
!  Add variable n_num for specifying a certain nest during output
!  Increase nvar to 110
! 2.8        2001/07/06 Ulrich Schaettler
!  Added new structures for handling the I/O and got some variables from 
!  data_runcontrol.f90
! 2.11       2001/09/28 Ulrich Schaettler
!  Added control variables for using cloud ice as initial and boundary fields
! 2.18       2002/07/16 Ulrich Schaettler
!  Added new Namelist variables lbd_frame, npstrframe to run LM with boundary
!  data defined on frames (work by Lucio Torrisi, UGM Rome).
! 3.2        2003/02/07 Ulrich Schaettler
!  Added new Namelist variable ilevbotnoframe to run LM with frames and 
!  Rayleigh damping (work by Lucio Torrisi, UGM Rome).
! 3.6        2003/12/11 Ulrich Schaettler
!  Introduced new variable ngribout for indicating, how many GRIBOUT namelist
!  groups are present. Increased character length for directory names.
! 3.7        2004/02/18 Ulrich Schaettler
!  New variables for specifying additional GRIB tables and an output list for
!  synthetic satellite images;
!  increased number of possible pressure levels for output
! 3.12       2004/09/15 Christoph Schraff
!  New namelist variable 'lana_qr_qs'to read prognostic rain and snow.
! 3.13       2004/12/03 Ulrich Schaettler
!  Eliminated Kind-parameters for Grib-library and put it to data_parameters
! 3.15       2005/03/03 Ulrich Schaettler
!  Introduced new NL Parameter: nunit_of_time
! 3.16       2005/07/22 Ulrich Schaettler
!  Introduced new NL parameter: lana_qg; ysuffix
!  Added new entries in structure ar_des (idef_stat) and pp_nl (ysuffix)
! 3.17       2005/12/12 Ulrich Schaettler
!  New Namelist variable lana_rho_snow to read prognostic snow density
!  Enlarged lanfld, to select snow density from analyses or interpolated
!  GME data.
! 3.18       2006/03/03 Ulrich Schaettler (Thomas Prenosil, CLM Community)
!  Control variables for writing restart files
!  Introduction of NetCDF and related variables, dimensions and structures.
!  New switch l_fi_ps_smooth for smoothing fi, pmsl in mountainous terrain
!  New variable in structure pp_nl for output of subdomain-fields
!  Increase maximum number of input model variables and constant variables
!      to incorporate fields of the lake model FLake  (Dmitrii Mironov)
! 3.21       2006/12/04 Jochen Foerstner, Burkhardt Rockel
!  New Namelist variables llb_qr_qs, llb_qg to read qr-, qs- and qg-values
!  from lateral boundaries
!  Changes in idims_id_out definitions; Added "undefncdf"
! V3_23        2007/03/30 Ulrich Schaettler
!  Added variables nexthour, lhour in structure pp_nl
! V4_1         2007/12/04 Ulrich Schaettler
!  Editorial changes
! V4_8         2009/02/16 Ulrich Schaettler
!  Increase value of nlevels and introduce a new variable for maximum of
!  actual existing output levels (gives possibility to use many p- or z-levels)
!  New NL variable for end of (total) simulation ydate_end
!  Add l_ke_in_gds to partly replace ldwd_grib_use
! V4_9         2009/07/16 Ulrich Schaettler
!  Increase number of IO-variables to 400 (for COSMO_ART)
! V4_11        2009/11/30 Jan-Peter Schulz
!  Enlarge lanfld to select sea ice temperature and thickness from analyses or
!  interpolated GME data.
! V4_12        2010/05/11 Ulrich Schaettler
!  Increase nzmxc to 40 (number of variables in list of constant fields)
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_17        2011/02/24 Ulrich Blahak
!  Increased max. number of height level variables to 400
!  Introduced namelist parameter "itype_vertint" in GRIBOUT namelist group to 
!    specify the type of vertical interpolation to p- and z-levels 
! V4_18        2011/05/26 Ulrich Schaettler
!  Moved NL parameter yform_write to the group(s) /GRIBOUT/ to be able to 
!     specify the output format differently for every group.
!  Initialize ndims_id_in here and do not read it in src_input, because the
!     different dimensions are really hard coded here (Anne Roches)
!  Introduced additional NetCDF dimension for 3D external parameter for topographical
!     correction: sectors of the horizon (by Anne Roches)
!  Introduced additional NetCDF dimensions for output of synthetic satellite
!     images (Anne Roches, et al.)
! V4_23        2012/05/10 Ulrich Schaettler (for CLM), Ulrich Blahak, Lucio Torrisi
!  Increased ndims_id_out to 18 to include ke_snow
!  Changed name of l_fi_ps_smooth to l_fi_pmsl_smooth and added l_fi_filter and
!   l_pmsl_filter in order to be able to independently smooth FI and PMSL
!   with a digital FIR filter, as for all other fields with l_z_filter / l_p_filter. (UB)
!   Introduced a new Namelist switch lbdsst, to update only SST over sea with
!      boundary data
! V4_24        2012/06/22 Burkhardt Rockel, Hendrik Reich
!  Adapted definition of dim_ids to INT2LM:
!    changed ID for topo corrections from 11 to 15; added some more IDs
!  Increased character strings for dates to 14 digits to include minutes and seconds
!  Introduced (internal) variable lmmss to indicate whether the 14 digits or the
!    10 digits format is used
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer, Ulrich Blahak, Carlos Osuna
!  Implemented a new parameter max_gribrep for the maximum allowed repetition
!    of one GRIB parameter (in src_setup_vartab)
!  Implemented namelist switch itype_gather to use different methods for gathering
!    output fields (OF)
!  Added new namelist switch "loutput_q_densities" in the TYPE pp_nl
!   for the gribout namelist(s). If set to .true., hydrometeor contents
!   are output in units of kg/m**3 resp. 1/m**3 instead of kg/kg resp. 1/kg. (UB)
!  Add idimvert to ar_des structure for asynchronous netcdf I/O. This integer
!   determines the length of vertical dimension of a field (-1 if undefined like 
!   for a 2d field)  (CO)
!  Added an integer variable to pp_nl that assigns an index that identifies each
!   namelist gribout group (CO)
! V4_26        2012/12/06 Burkhardt Rockel (CLM)
!  Changes in NetCDF I/O: change some variables of dimension with lenth 1 to scalar
!  and reduce the number of dimension IDs
! V4_28        2013/07/12 Ulrich Schaettler
!  Introduced new organizational variables for Grib handling:
!    ylevltypes1, ylevltypes2, ysteptypes, itabletypes, idwdednr (old: iednr), igrbednr
!    variables for grib_api samples,
!  Introduced new entries in type pp_nl: 
!    igribapi_id: to store the grib-handler for the output list
!    lsfc_ana:    to indicate whether a surface analysis field will be written
!  Eliminated other entries which are not used any more: 
!    n_num, ysystem, ydbid, ydbtype, ydbpw
!  New namelist variable in /GRIBIN/ for HHL file: ydirhhl, ynamhhl
! V4_29        2013-10-02 Astrid Kerkweg, Ulrich Schaettler
!  Introduced parameter variable clen for string length of shortnames
!  Introduced namelist switch lan_w_so to choose W_SO from nudging or external 
!    analysis; increased nanfld to 15.
! V4_30        2013/11/08 Ulrich Schaettler
!  Split ipds into ipds_in and ipds_out, to separate for input and output
! V5_1         2014-11-28 Ulrich Schaettler, Ulrich Blahak, Oliver Fuhrer
!  Eliminated ydirhhl, ynamhhl again (no more needed)
!  New variable nrbit_phhl as packrate for hhl and p-fields (US)
!  New variable ngribednr to save input GRIB edition number (US)
!  Added new components "ydir_restart_in" and "ydir_restart_out" to structure pp_nl
!   to store the namelist parameter values for the directories where restart files 
!   are read and written (UB)
!  Add namelist variable dbz of TYPE(dbzcalc_params) to type pp_nl for
!   output GRIBOUT namelist. It contains many configuration parameters for
!   radar reflectivity calculations.
!  Replaced ireals by wp (working precision) (OF)
! V5_4         2016-03-10 Oliver Fuhrer, Pascal Spoerri, Xavier Lapillonne
!  Added a data structure to track boundary fields in case of GPU compilation
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
USE data_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    iintegers, & ! KIND-type parameter for standard integer variables
    irealgrib, & ! KIND-type parameter for real variables in the grib library
    intgribf     ! KIND-type parameter for fortran files in the grib library

#ifdef RADARFWO
USE data_radar, ONLY: dbzcalc_params  ! type to hold radar reflectivity meta data
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================

! 1.1. Global parameters, dimensions and arrays for the GRIB-library
! ------------------------------------------------------------------

#ifdef MESSY
INTEGER(KIND=iintegers), PARAMETER :: clen = 49
#else
INTEGER(KIND=iintegers), PARAMETER :: clen = 10
#endif

INTEGER (KIND=intgribf)             :: &
  idwdednr =  1,    & ! grib edition number for DWD library
  igrbednr            ! "working" grib edition number (to be set during run-time)

INTEGER (KIND=intgribf),  PARAMETER :: &
  npds   =    321_intgribf, & ! Dimension for product definition section (pds)
  ngds   =    626_intgribf, & ! Dimension for grid description section (gds)
  nbms   =      3_intgribf, & ! Dimension for bit map section (bms)
  nbds   =     11_intgribf, & ! Dimension for binary data section
  ndsup  =     73_intgribf, & ! Dimension for dsup
  ndims  =     20_intgribf    ! Dimension for idims (contains all dimensions)

! Grib handles for the grib_api sample data
INTEGER (KIND=intgribf)             :: &
  igrib1_id, igrib2_id

! The following dimensions are set during program execution
INTEGER (KIND=intgribf)             :: &
  lfd,                      & ! Dimension for iblock
  lfa,                      & ! Dimension for grib_api message in bytes
  lbm,                      & ! Dimension for bitmap: this value has to be at
                              ! least the same as lbmax in subroutine grbin1
  lds,                      & ! Dimension for unpacked data
  inrvert_in,               & ! number of vertical coordinate parameters of input data
  inrvert_out                 ! number of vertical coordinate parameters of output data

! Global Arrays:

INTEGER (KIND=intgribf)              :: &
  idims_in  (ndims),  & ! array for all dimensions (input)
  idims_out (ndims),  & ! array for all dimensions (output)
  ipds_in   (npds),   & ! product definition section for input
  ipds_out  (npds),   & ! product definition section for output
  igds_in   (ngds),   & ! grid description section for input
  igds_out  (ngds),   & ! grid description section for output
  ibms      (nbms),   & ! bit map section
  ibds      (nbds)      ! binary data section

! packed grib_field -> grib_api
CHARACTER(LEN=1), DIMENSION(:), ALLOCATABLE :: ymessage

! Arrays that are allocated during program execution
INTEGER (KIND=intgribf), ALLOCATABLE :: &
  iblock    (:),      & ! array for gribed data
  ibmap     (:)         ! array for bit map

REAL   (KIND=irealgrib), ALLOCATABLE :: &
  dsup      (:),      & ! array for special data
  ds_grib   (:)         ! array for unpacked data

REAL   (KIND=wp),        ALLOCATABLE :: &
  pv_in     (:),      & ! array for vertical coordinate parameters for input data
  pv_out    (:),      & ! array for vertical coordinate parameters for output data
  ds_real   (:)         ! array for unpacked data

! Arrays to convert GRIB1 level types and time range indicators to grib_api strings
! for GRIB1 (ylevltypes1) and GRIB2 (ylevltypes2), resp.
CHARACTER(LEN=30)  :: ylevltypes1(255), ylevltypes2(255)
CHARACTER(LEN=10)  :: ysteptypes(0:10)

! Array to convert GRIB2 scale factors to real numbers
REAL (KIND=wp)     :: rscalefac(0:9)

!==============================================================================

! 1.2. Global parameters, dimensions and arrays for the netcdf-library
! --------------------------------------------------------------------

INTEGER (KIND=iintegers), PARAMETER   :: &
  ndims_id_in =14, & ! number of dimensionID's for netCDF formatted input
  ndims_id_out=14    ! number of dimensionID's for netCDF formatted output


INTEGER (KIND=iintegers), ALLOCATABLE :: &
  idims_id_in(:)  ! for the IDs of the dimensions of netCDF formatted input
                  ! The different dimensions are:
                  !   1: ID for ie,                                 "rlon"
                  !   2: ID for je,                                 "rlat"
                  !   3: ID for ke,                                 "level"
                  !   4: ID for ke1,                                "level1"
                  !   5: ID for ntime,                              "time"
                  !   6: ID for nbnds,                              "nbds"
                  !   7: ID for ke_soil,                            "soil"
                  !   8: ID for ke_soil+1,                          "soil1"
                  !   9: ID for staggered ie,                       "srlon"
                  !  10: ID for staggered je                        "srlat"
                  !  11: ID for sections of topo ext. parameter     "nhori"
                  !  12: ID for all products of syn. sat. data      "nsynmsg"
                  !  13: ID for one group of products of MSG        "msgchan"
                  !  14: ID for ke_snow

INTEGER (KIND=iintegers)              :: &
  idims_id_out (ndims_id_out)
                  ! for the IDs of the dimensions of netCDF formatted output
                  ! The different dimensions are:
                  !   1: ID for ie,                                 "rlon"
                  !   2: ID for je,                                 "rlat"
                  !   3: ID for ke,                                 "level" ("pressure", "altitude")
                  !   4: ID for ke1,                                "level1"
                  !   5: ID for ntime,                              "time"
                  !   6: ID for nbnds,                              "bnds"
                  !   7: ID for ke_soil,                            "soil"
                  !   8: ID for ke_soil+1,                          "soil1"
                  !   9: ID for staggered ie,                       "srlon"
                  !  10: ID for staggered je,                       "srlat"
                  !  11: ID for sections of topo ext. parameter     "nhori"
                  !  12: ID for all products of syn. sat. data      "nsynmsg"
                  !  13: ID for one group of products of MSG        "msgchan"
                  !  14: ID for ke_snow

! NetCDF global attributes
CHARACTER (LEN=200)  :: &
  yncglob_institution,   & ! originating center name
  yncglob_title,         & ! title string for the output
  yncglob_source,        & ! program name and version
  yncglob_project_id,    & ! identification of the project of simulation
  yncglob_experiment_id, & ! identification of the experiment of simulation
  yncglob_contact,       & ! contact e.g. email address
  yncglob_references       ! URL, report etc.

INTEGER   (KIND=iintegers)          :: &
  ncglob_realization       ! number of the realization of the experiment

!==============================================================================

! 2. Parameters and dimensions for I/O
! ------------------------------------

INTEGER (KIND=iintegers), PARAMETER :: &
  ntrip = 30, &       ! maximum number of timing triples
  noutst_max = 1000

INTEGER (KIND=iintegers), PARAMETER :: &
  nzmxin = 400, & ! maximum number of input model variables
  nzmxml = 400, & ! maximum number of output model-level variables
  nzmxpl = 400, & ! maximum number of pressure-level variables
  nzmxzl = 400, & ! maximum number of height-level variables
  nzmxc  =  40, & ! maximum number of constant variables
  nzmxid = 400    ! maximum number of NetCDF variabe IDs
                  !  (maximum of all the above)

INTEGER (KIND=iintegers), PARAMETER :: &
  nlevels= 500    ! maximum number of pressure or height levels
     
INTEGER (KIND=iintegers)            :: &
  noutlevels      ! maximum actual existing number of output levels
     
INTEGER (KIND=iintegers), PARAMETER :: &
  max_gribtabs =  20, & ! maximum number of GRIB tables in LM variable table
  max_gribrep  =   4, & ! maximum number of allowed repetition of the same
                        ! GRIB parameter number - GRIB table number combination
  nvar         = 110    ! maximum number of variables in LM variable table

LOGICAL                             :: &
  lvar            ! indicates, whether the LM variable table is allocated

REAL    (KIND=irealgrib), PARAMETER :: &
  undefgrib =  -1E7_irealgrib ! value for "undefined" in the grib routines

REAL    (KIND=irealgrib)            :: &
! value for "undefined" in the netcdf routines
  undefncdf = -1.E20_irealgrib

REAL    (KIND=wp)                   :: &
  undef             ! the same as undefgrib but with other KIND-Parameter

! Global Scalars:

INTEGER (KIND=iintegers)            :: &
  nlocaldefnr,    & ! local definition number for GRIB2 local section (Namelist parameter)
  nactlocdefnr,   & ! to overwrite Namelist parameter with some center default
  ngribednr,      & ! to store GRIB edition number for input
  nprocess_ini_in,& ! process generating identification for initial (analysis) 
  nprocess_bd_in, & ! and for boundary (forecasts) data from input data
  ngribout,       & ! number of GRIBOUT namelist groups
  ncenter,        & ! originating center identification
  nsubcenter,     & ! originating sub-center identification
  num_gribtabs,   & ! number of GRIB tables used in LM variable table
  lst_gribtabs(max_gribtabs), & ! IDs of GRIB tables use
  itabletypes(255)  ! Array to convert GRIB1 table types to the index 
                    ! in structure lst_gribtabs

INTEGER   (KIND=iintegers)          :: &
  nr_griboutnl, & ! number of the gribout namelists
  noutst          ! number of output timesteps (to be calculated in calc_ngrib)

INTEGER (KIND=iintegers)            :: &
  itype_gather    ! Switch to determine gather method to use
                  !  = 1 use MPI_GATHER to gather each 2D field seperately
                  !  = 2 use MPI_ALL2ALLV to gather num_compute 2D fields at once

!==============================================================================

! 3. controlling the netcdf/grib-I/O
! ----------------------------------

! initial data
! ------------

  CHARACTER (LEN= 250) :: ydirini         ! catalog-name of the initial file 
  CHARACTER (LEN=  14) :: ydate_ini       ! start of the forecast
  CHARACTER (LEN=  14) :: ydate_end       ! end of the (total) forecast
  CHARACTER (LEN=clen) :: yvarini(nzmxin) ! list of initial fields


  INTEGER   (KIND=iintegers)       ::           &
    nsma_stat ! status for soil humidity analysis

  INTEGER   (KIND=iintegers) ,PARAMETER  ::     &
    nanfld=15       ! max. number of input fields to be checked for
                    ! time range indicator itri=0

  INTEGER   (KIND=iintegers)       ::           &
    nyvar_i,      & ! number of variables in the file with the initial data
    nversini,     & ! version number of initial data (in grib-code)
    ndabaini        ! minimal number of fields in initial data base

  LOGICAL                          ::           &
    lchkini,      & ! checking the initial data
    lanfld(nanfld)  ! contains switches for all fields to be checked for
                    ! time range indicator = 0

  CHARACTER (LEN=clen)             ::           &
    ynaman(nanfld)  ! name of fields to be checked for time range indicator

  DATA ynaman /   &
        'T_S       ', 'T_SNOW    ', 'T_CL      ',         &
        'W_SNOW    ', 'W_I       ', 'W_CL      ',         &
        'VIO3      ', 'HMO3      ', 'PLCOV     ',         &
        'LAI       ', 'ROOTDP    ', 'RHO_SNOW  ',         &
        'T_ICE     ', 'H_ICE     ', 'W_SO      ' /

! boundary data
! -------------

  CHARACTER (LEN= 250) :: ydirbd          ! catalog-name of the file
  CHARACTER (LEN=   1) :: ytunitbd        ! time unit for boundary data
  CHARACTER (LEN=  14) :: ydate_bd        ! start of the forecast from which
                                          ! the boundary fields are used
  CHARACTER (LEN=clen) :: yvarbd (nzmxin) ! list of boundary fields

  INTEGER   (KIND=iintegers)       ::           &
    nyvar_b,      & ! number of variables in the file with the boundary data
    nversbd,      & ! version number of boundary data (in grib-code)
    ndababd         ! minimal number of fields in boundary data base

  LOGICAL                          ::           &
    lchkbd,       & ! checking the boundary data
    lbdclim,      & ! boundary data in climate model     ! PIK  (D.Hauffe)
                    ! (in climate mode also some external parameters have
                    !  to be updated, which are held constant in forecast
                    !  mode; e.g. plant cover, root depth)
    lbdsst          ! T_S boundary data are provided only over sea
                    ! (SST is not maintained constant during the integration)

! reading and writing of ready files
! ----------------------------------

  CHARACTER (LEN=250)              ::           &
    ytrans_in,    & ! directory for reading ready-files
    ytrans_out      ! directory for writing ready-files

  INTEGER   (KIND=iintegers)       ::           &
    nincwait,     & ! if ready-file is not available wait nincwait seconds
                    ! until next attempt
    nmaxwait        ! if ready-file is not available after nmaxwait seconds,
                    ! abort the program

! reading and writing of restart files
! ------------------------------------

  INTEGER   (KIND=iintegers)       ::           &
    nhour_restart(3)! start-, stop-, increment of writing restart files
                    ! (in time steps)

  CHARACTER (LEN=250) :: ydir_restart_in    ! directory for reading restart file
  CHARACTER (LEN=250) :: ydir_restart_out   ! directory for writing restart file(s)
  CHARACTER (LEN=  1) :: ytunit_restart     ! unit of timescale

#ifdef RADARFWO
! reading and writing of Mie lookup tables:
! ------------------------------------

  CHARACTER (LEN=250) :: ydir_mielookup     ! Directory to store/read Mie lookuptables
#endif

! Variables for handling the NetCDF-/Gribfile I/O:
! ------------------------------------------------

  CHARACTER (LEN=  4) ::  &
    yform_read      ! format of the (read) files

  CHARACTER (LEN=  3) :: ymode_read     ! mode for opening the (read) Grib files
  CHARACTER (LEN=  3) :: ymode_write    ! mode for opening the (write) Grib files
  CHARACTER (LEN=  8) :: yuchkdat = 'YUCHKDAT'       ! checking the I/O data

  LOGICAL                          ::           &
    ldwd_grib_use,& ! use some DWD specific Grib settings
    l_ke_in_gds,  & ! explicit GDS entry for number of model levels
    l_ke_in_input,& ! indicates whether GRIB1 input data contains ke in meta data
    lbdana,       & ! boundary data are analysed data
    lana_qi,      & ! if .TRUE., take qi-values from analysis file
                    ! else, qi is set in the model
    llb_qi,       & ! if .TRUE., take qi_bd-values from lateral boundaries file
                    ! else, qi_bd is set in the model
    lana_qr_qs,   & ! if .TRUE., take qr- and qs-values from analysis file
                    ! else, qr and qs are set in the model
    llb_qr_qs,    & ! if .TRUE., take qr_bd- and qs_bd-values from lateral
                    ! bound. file else, qr_bd and qs_bd are set in the model
    lana_qg,      & ! if .TRUE., take qg-values from analysis file
                    ! else, qg is set in the model
    llb_qg,       & ! if .TRUE., take qg_bd-values from lateral boundaries file
                    ! else, qg_bd is set in the model
    lana_rho_snow,& ! if .TRUE., take rho_snow-values from analysis file
                    ! else, it is set in the model
    lbd_frame       ! if .TRUE., boundary data are on a frame

  INTEGER (KIND=iintegers)     ::    &
    npstrframe,   & ! width (number of points) of the strip around the b.d.frame
    ilevbotnoframe  ! bottom model level with b.d. defined on the whole grid
                    ! (model levels below ilevbotnoframe are defined on a frame)

  ! Unit Numbers for ASCII files related to I/O
  INTEGER (KIND=iintegers)     ::    &
    nuin,       & ! Unit number for Namelist INPUT files
    nuchkdat,   & ! Unit number for checking the I/O data
    ntrans_out    ! Unit Number for writing ready-Files during output

  LOGICAL                     ::    &
    lmmss         ! if .TRUE.  14 digits date format (YYYYMMDDHHMMSS)
                  ! if .FALSE. 10 digits date format (YYYYMMDDHH)

! 4. Defined data types for I/O
! -----------------------------

! Type for LM variable table
TYPE ar_des
  CHARACTER (LEN=clen)           :: name
  INTEGER(KIND=intgribf)         :: levtyp
  INTEGER(KIND=intgribf)         :: levtop
  INTEGER(KIND=intgribf)         :: levbot
  REAL(KIND=wp)                  :: factor
  REAL(KIND=wp)                  :: bias
  INTEGER(KIND=intgribf)         :: ntri
  INTEGER(KIND=intgribf)         :: rank
  INTEGER(KIND=intgribf)         :: idimvert
  REAL(KIND=wp), POINTER         :: p4   (:,:,:,:)
  REAL(KIND=wp), POINTER         :: p4_bd(:,:,:,:)
  REAL(KIND=wp), POINTER         :: p3   (:,:,:)
  REAL(KIND=wp), POINTER         :: p3_bd(:,:,:)
  REAL(KIND=wp), POINTER         :: p2   (:,:)
  INTEGER(KIND=iintegers)        :: idef_stat
! the following is for NetCDF
  CHARACTER (LEN=20)             :: units
  CHARACTER (LEN=80)             :: standard_name
  CHARACTER (LEN=80)             :: long_name
  CHARACTER (LEN=1)              :: lsm
END TYPE ar_des

TYPE (ar_des),     ALLOCATABLE   :: var(:,:,:)

! Type for the list descriptions
TYPE list_description
  CHARACTER (LEN=clen)           :: name                ! name of variable
  INTEGER (KIND=iintegers)       :: iloc1, iloc2, iloc3 ! location in vartab
  INTEGER (KIND=iintegers)       :: idimvert            ! vertical dimension
END TYPE

! Variables for describing the different lists
TYPE (list_description), ALLOCATABLE :: list_ini   (:), list_bd    (:), &
                                        list_res_o (:), list_res_n (:)

! packrate for HHL and P:
INTEGER (KIND=iintegers)   :: nrbit_phhl

! Type for output namelist
TYPE pp_nl
  INTEGER (KIND=iintegers)   :: nl_index           ! index of the current output group
  INTEGER (KIND=iintegers)   :: igribapi_id        ! id of grib record for grib_api
  CHARACTER (LEN=clen)       :: yvarml(nzmxml)     ! list of m-variables
  CHARACTER (LEN=clen)       :: yvarpl(nzmxpl)     ! variables of p-levels
  CHARACTER (LEN=clen)       :: yvarzl(nzmxzl)     ! variables of z-levels
  CHARACTER (LEN=clen)       :: yvarsl(nzmxzl)     ! variables of satellite-levels
  CHARACTER (LEN=clen)       :: yvarc (nzmxc )     ! constant variables
  INTEGER (KIND=iintegers)   :: ilist_ml(3,nzmxml) ! location in variable table
  INTEGER (KIND=iintegers)   :: ilist_pl(3,nzmxpl) ! location in variable table
  INTEGER (KIND=iintegers)   :: ilist_zl(3,nzmxzl) ! location in variable table
  INTEGER (KIND=iintegers)   :: ilist_sl(3,nzmxzl) ! location in TYPE sat_org_table
  INTEGER (KIND=iintegers)   :: ilist_c (3,nzmxc ) ! location in variable table
! INTEGER (KIND=iintegers)   :: n_num              ! current grid number
  INTEGER (KIND=iintegers)   :: nyvar_m            ! number of variables
  INTEGER (KIND=iintegers)   :: nyvar_p            ! number of variables
  INTEGER (KIND=iintegers)   :: nyvar_z            ! number of variables
  INTEGER (KIND=iintegers)   :: nyvar_s            ! number of variables
  INTEGER (KIND=iintegers)   :: nyvar_c            ! number of variables
  INTEGER (KIND=iintegers), POINTER  :: ngrib(:)   ! list of outputsteps
  INTEGER (KIND=iintegers)   :: outsteps           ! number of outputsteps
  INTEGER (KIND=iintegers)   :: nextstep           ! next output step
  INTEGER (KIND=iintegers)   :: nexthour           ! next output hour
  LOGICAL                    :: lhour              ! output steps are rounded
                                                   ! to the next multiple of 0.25

  ! process generating identification for analysis and for forecast
  ! data that are written
  INTEGER (KIND=iintegers)   :: nprocess_ini_out   !
  INTEGER (KIND=iintegers)   :: nprocess_bd_out    !

  ! indicator for unit-of-time (1hr, 15min, 30min,...)
  INTEGER (KIND=iintegers)   :: nunit_of_time      !

  ! the following variables are for specifying the output domain
  REAL    (KIND=wp)          :: slon               ! left longitude
  REAL    (KIND=wp)          :: slat               ! botton latitude
  REAL    (KIND=wp)          :: elon               ! right longitude
  REAL    (KIND=wp)          :: elat               ! top latitude
  INTEGER (KIND=iintegers)   :: i_out_start        ! start i-index
  INTEGER (KIND=iintegers)   :: j_out_start        ! start j-index
  INTEGER (KIND=iintegers)   :: i_out_end          ! end   i-index
  INTEGER (KIND=iintegers)   :: j_out_end          ! end   j-index
  INTEGER (KIND=iintegers)   :: ie_out_tot         ! # i grid points
  INTEGER (KIND=iintegers)   :: je_out_tot         ! # j grid points

  ! additional stuff
! CHARACTER (LEN=  4)        :: ysystem            ! file, db or ecfs
  CHARACTER (LEN=250)        :: ydir               ! directory of file or db
  CHARACTER (LEN=250)        :: ydir_restart_in    ! directory for reading restart files
  CHARACTER (LEN=250)        :: ydir_restart_out   ! directory for writing restart files
  CHARACTER (LEN= 25)        :: ysuffix            ! optional filename suffix
  CHARACTER (LEN=  1)        :: ytunit             ! unit of timescale
  CHARACTER (LEN=  1)        :: ydomain            ! sign for the domain
  CHARACTER (LEN=  4)        :: yform_write        ! format of the output files
  INTEGER (KIND=iintegers)   :: nrbit              ! packrate
! CHARACTER (LEN=  5)        :: ydbid              ! database-id
! CHARACTER (LEN=  5)        :: ydbtype            ! database-typ
! CHARACTER (LEN= 10)        :: ydbpw              ! database password
  REAL    (KIND=wp)          :: plev(nlevels)      ! pressure levels
  REAL    (KIND=wp)          :: zlev(nlevels)      ! z-levels
  INTEGER (KIND=iintegers)   :: kepin              ! number of pressurelevels
  INTEGER (KIND=iintegers)   :: kezin              ! number of z-levels
  LOGICAL                    :: lcheck             ! check output ?
  LOGICAL                    :: lwrite_const       ! write constant fields
  LOGICAL                    :: luvmasspoint       ! interpolate horizontal 
                                                   ! winds to mass grid points
  LOGICAL                    :: lanalysis          ! write an analysis file
                                                   ! (possible with nudging)
  LOGICAL                    :: lsfc_ana           ! write a surface analysis
  LOGICAL                    :: l_p_filter         ! switch for smoothing p-interpol. fields
  LOGICAL                    :: l_z_filter         ! switch for smoothing z-interpol. fields
  LOGICAL                    :: l_fi_filter        ! switch for independent smoothing of fi
  LOGICAL                    :: l_pmsl_filter      ! switch for independent smoothing of pmsl
  LOGICAL                    :: l_fi_pmsl_smooth   ! switch for additional smoothing of fi, pmsl
                                                   !  in mountaneous regions
  LOGICAL                    :: loutput_q_densities! switch for unit of QX, NCX:
                                                   !     kg/kg    resp. 1/kg   (.FALSE.) or
                                                   !     kg/m**3  resp. 1/m**3 (.TRUE.)
  INTEGER (KIND=iintegers)   :: itype_vertint      ! type of vertical interpolation  
  TYPE(pp_nl), POINTER       :: next               ! pointer to next nl

#ifdef RADARFWO
  TYPE(dbzcalc_params)       :: dbz
#endif

END TYPE pp_nl

TYPE (pp_nl), POINTER      :: root, pp_restart, pp_lansfc

! data structure to track boundary fields read from file
TYPE bd_field
  CHARACTER(LEN=10)        :: name
  INTEGER (KIND=iintegers) :: rank
  INTEGER (KIND=iintegers) :: iloc1
  INTEGER (KIND=iintegers) :: iloc2
  INTEGER (KIND=iintegers) :: iloc3
  INTEGER (KIND=iintegers) :: ntlev
  REAL (KIND=wp), pointer  :: p3(:,:,:)
  REAL (KIND=wp), pointer  :: p2(:,:)
end TYPE bd_field

INTEGER, PARAMETER       :: max_bd_fields = 100_iintegers
INTEGER (KIND=iintegers) :: num_bd_fields = 0_iintegers
TYPE(bd_field)           :: bd_list(max_bd_fields)

!==============================================================================

END MODULE data_io
