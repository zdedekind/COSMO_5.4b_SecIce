!+ Source Module for reading in Grib files
!------------------------------------------------------------------------------

MODULE src_input

!------------------------------------------------------------------------------
!
! Description:
!   This module contains subroutines necessary for reading the initial and the
!   boundary data for LM. It uses also routines from the module "io_utilities".
!
! Current Code Owner: DWD, Ulrich Schaettler
!    phone:  +49  69  8062 2739
!    fax:    +49  69  8062 3721
!    email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Ulrich Schaettler
!  Initial release
! 1.4        1998/05/22 Guenther Doms
!  Adaptions for two time-levels time integration scheme
! 1.5        1998/06/29 Guenther Doms
!  Use of the new fiels hhlr (half-level height refering to flat terrain)
! 1.7        1998/07/16 Guenther Doms
!  Removal of routines caldpst and calrssk from use in meteo_utilities.
! 1.8        1998/08/03 Ulrich Schaettler
!  Use grib parameters from module data_io.f90
! 1.9        1998/09/16 Guenther Doms
!  Use of parameters 'nincmxt' and 'nincmxu' (replacing 'nincmxn') from
!  data module 'data_runcontrol.f90'.
! 1.10       1998/09/29 Ulrich Schaettler
!  Adaptions of the Use lists.
! 1.11       1998/10/13 Christoph Schraff
!  Additional variables for selecting analysis fields.
! 1.14       1998/10/26 Ulrich Schaettler
!  Changed igds to igds_in.
! 1.17       1998/11/17 Ulrich Schaettler
!  Changes for reading ready files.
! 1.20       1999/01/07 Guenther Doms
!  Renaming of some global variables
! 1.29       1999/05/11 Ulrich Schaettler
!  Adapted interfaces to utility-modules and prepared use of MPE_IO
! 1.30       1999/06/24 Matthias Raschendofer
!  Use an additional parameter form module data_runcontrol: nstop
! 1.34       1999/12/10 Ulrich Schaettler
!  Renamed some boundary fields and use new Namelist variables
!  Include all module procedures in this file.
! 1.37       2000/03/24 Guenther Doms
!  Additional job-log printout for input fields with additional element number
! 1.39       2000/05/03 Ulrich Schaettler
!  Changed names for variables concerned to latitude or longitude.
!  Introduced possibility for parallel asynchronous and/or database IO.
! 1.41       2000/06/02 Ulrich Schaettler
!  Correction in Subroutine id_mix_bd.
! 2.8        2001/07/06 Ulrich Schaettler
!  Eliminated non-necessary variables from the USE-lists;
!  Adapted input to new organization of I/O
! 2.11       2001/09/28 Ulrich Schaettler
!  Corrected a bug for reading a ready file for initial data.
!  Corrected reading variables for multi-layer soil model.
! 2.14       2002/02/15 Ulrich Schaettler
!  Modifications to allow the use of the SLEVE coordinate
! 2.17       2002/05/08 Ulrich Schaettler
!  Modifications to perform I/O-communications in irealgrib-format
! 2.18       2002/07/16 Ulrich Schaettler
!  Eliminated variable rhde; use cf_snow instead; 
!  changed dimensions of igds_in ipds for calling subroutine check_lm_grid;
!  Changes for working with boundaries defined on frames (by Lucio Torrisi)
!  Use lmulti_layer to define boundary values for old soil model correctly
! 3.2        2003/02/07 Ulrich Schaettler
!  If Rayleigh damping and working with frames is activated, the full grid of 
!  the boundary fields is checked in the Rayleigh damping layer.
! 3.6        2003/12/11 Ulrich Schaettler
!  Adaptations for multi-layer soil model and additional checks for data input
! 3.7        2004/02/18 Ulrich Schaettler
!  Changed treatment of unit-numbers for ASCII-file handling
! 3.13       2004/12/03 Ulrich Schaettler
!  Put KIND-parameters for Grib library to data_parameters;
!  Adapted interfaces of gather_field and reference_atmosphere
!  Explicit formulation of relaxation for lateral boundaries (factor rmy)
!                                                 (Jochen Foerstner)
! 3.14       2005/01/25 Ulrich Schaettler
!  Initialization of t_so(0) for all time levels;
!  Initializations for Rayleigh-damping only for itype_spubc==1
!  Adjusted computation of llandmask: fr_land=0.5 is a land point now!!!
! 3.15       2005/03/03 Ulrich Schaettler
!  Calculate surface pressure again after reading first boundary data set
!  (to be consistent with atmospheric pressure after id_mix_bd)
!  Modifications to read boundary files also in 0.5/0.25 hourly increments
! 3.16       2005/07/22 Ulrich Schaettler, Helmut Frank
!  Adapted length of pathname for input-directory (U.S.)
!  Adaptation to read W_SO also from the soil moisture analysis
!  Read T_SO(1:ke_soil+1) from the Nudging, even if T_SO(0) is read from 
!  the SST analysis                                                 (H.F.)
! 3.18       2006/03/03 Ulrich Schaettler
!  Changes to read NetCDF data format and restart files (binary data format)
!  Changed treatment of ASCII files for introducing restart possibility
!  Changes to introduce new type of coding the vertical coordinate parameters
!    (introduced by MeteoSwiss for the SLEVE vertical coordinate)
! 3.19       2006/04/25 Ulrich Schaettler
!  Adaptations to read "old" NetCDF data and restart data
! 3.21       2006/12/04 Ulrich Schaettler, Burkhardt Rockel
!  Initialize error code (implcode) variable in scatter_data
!  Some adaptations for reading T_S_LAKE for the FLake-Model
!  Changed "tbnds" to "bnds" in read_nc_gdefs
!  Changed from soil layer boundaries to mid soil layer
!  Deleted Z0 -> GZ0 replacement as special case for netCDF
!  Included klv950, klv700 in interface to reference_atmosphere
! 3.22       2007/01/24 Jochen Foerstner
!  Corrections for Runge-Kutta Restart: variables have to be put to "nnew"
! V3_23        2007/03/30 Ulrich Schaettler
!  Adaptation of check_required for Runge-Kutta Restart
!  Introduced idbg_level
! V3_24        2007/04/26 Ulrich Schaettler
!  SR check_required has to work with global idbg_level, not the local izdebug
! V3_25        2007/05/21 Ulrich Schaettler
!  Corrections for reading data from the old 2-layer soil model
! V4_1         2007/12/04 Ulrich Schaettler
!  Call to SR check_record: changed argument for 3rd dimension of variables
!  Call to SR sleve_split_oro: introduced my_cart_id as argument
! V4_3         2008/02/25 Ulrich Schaettler
!  Corrected reading of boundary data with ytunit_bd='d'
! V4_4         2008/07/16 Ulrich Schaettler
!  More detailed timing for input
!  Adapted interface of get_timings
!  Changed NL parameter lyear_360 to itype_calendar, to have several options
!  De-grib nfltvc (for SLEVE coordinate) by REFSTF in any case
! V4_5         2008/09/10 Guenther Zaengl, Christoph Gebhardt
!  Adaptations for new reference atmosphere (G. Zaengl)
!  Modifications for lai, plcov, rootdp in case of EPS mode (C. Gebhardt)
! V4_8         2009/02/16 Guenther Zaengl, Oliver Fuhrer
!  Use p0hl (reference pressure at half levels) for full consistency with
!  new reference atmosphere implementation (Guenther)
!  Read vertical coordinate parameters only from 3d variables (levtyp 109/110)
!  Bug correction in check_required for Grib table number 205 (Oliver)
! V4_9         2009/07/16 Ulrich Schaettler, Burkhardt Rockel
!  Introduced new GRIB tables for COSMO_ART (241,242) and MCH (250)
!  Optimized loops in SR id_mix_bd
!  Compute tgcom for the whole domain (otherwise t_g is not set at the boundaries)
!  Bug fix for new reference atmosphere in netCDF input
! V4_11        2009/11/30 Jan-Peter Schulz, Ekaterina Machulskaya
!  Correction of computation for nlandpoints (fr_land must be >= 0.5) (JPS)
!  Adaptations for multi-layer snow model (EM)
! V4_12        2010/05/11 Jan-Peter Schulz, Ulrich Schaettler
!  Adapted interface parameter t_s for SR tgcom to account for seaice
!  Introduced SR check_vertcoord to check consistency of vertical coordinates
!    for every field (US)
!  Implemented a check for soiltyp - frland compatibility (Oli Fuhrer)
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_15        2010/11/19 Ulrich Schaettler, Lucio Torrisi
!  The field t_s_lake is removed again after adaptations in the SST Analysis
!  Increase lbm to 25000 for working with big frames (LT)
! V4_16        2010/12/07 Jan-Peter Schulz
!  Initialization of t_ice for all time levels (according to t_so)
! V4_17        2011/02/24 Ulrich Blahak
!  Replaced INT by NINT and format 2I2 by 2I2.2 in organize_input(), section 3.3
! V4_18        2011/05/26 Ulrich Schaettler
!  Adapted some printouts in the check of vertical coordinate parameters
!  Bug fix for ytunitbd/='d' (Burkhardt Rockel, et. al)
!  Extended interface to organize_input to allow missing values in list of input data
!    (Christoph Knote)
!  Adapted NetCDF I/O to deal with 3D external parameter field for sectors of
!   the horizon (for topographical corrections) and its attributes (Anne Roches)
! V4_19        2011/08/01 Ulrich Schaettler
!  Introduced conditional compilation for NetCDF and GRIBDWD
! V4_23        2012/05/10 Ulrich Schaettler, CLM, Lucio Torrisi
!  Adapted interfaces of check_input_grid, read_netcdf according to changes 
!    in io_utilities.f90
!  Adapted call to SR distribute_fields (added sender PE) (Uli Schaettler)
!  Add support for climatological year with 365 days
!    (replaced calls to SR difmin, difmin_360 by call to new SR diff_minutes)
!  In case of restart: change input directory of restart file from ydirini to ydir_restart
!  For 'binary' input allow T_SO(ksoil=0) to be read instead of T_S
!  Added field rho_snow_mult, RHO_SNOW_M
!  Correction in fill_realarray for level determination of multi-layer snow model
!  In case of lbdsst, T_S must not be defined on frames
! V4_24        2012/06/22 Ulrich Schaettler, Burkhardt Rockel,
!                         Michael Baldauf, Hendrik Reich
!  Reading restart files: Also check the possibility of irefatm=3 for restarts
!  Corrected an error message in SR fill_realarray
!  Adapted definition of dim_ids to INT2LM:
!      changed ID for topo corrections from 11 to 15 (by Burkhardt Rockel)
!  Adapted calls to SR reference_atmosphere_x  (Michael Baldauf)
!  Extension of date variables to minutes and seconds
!  ndiff_ini_bdref introduced to account for difference between boundary date
!      and reference date          (Hendrik Reich)
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer, Ulrich Blahak
!                         Florian Prill, CLM
!  Replaced qx-variables by using them from the tracer module
!  Bug fix for processing restart- and boundary files with ytunit='d': To get
!  the proper file name in the old 10-digit format, the minutes and seconds
!  of the date variable have to be set to '0000' (Uli Schaettler)
!  Implemented initialization of additional tracers (QHAIL, NCXXX) for
!   the 2-moment cloud microphysics scheme (Uli Blahak).
!  Reading and distribution of parameter bvref for irefatm=3 was missing for
!   binary restart files (Uli Blahak).
!  Bugfix: Calculation of qrs and rho was missing for initial data. (UB)
!  Adapted call to make_fn according to changes in io_utilities (Uli Schaettler)
!  Check soiltyp with NINT-function
!  Corrections in the 10/14 digit format, when working with NetCDF I/O (by CLM, HJP)
!  Make read of reference date in netCDF input more flexible; now years below 
!       1000 can also be read (BR)
!  Enhanced interface to organize_input due to possible prefetching of data (FP)
!  New Subroutine "create_file_name" (Section 1.3), which is now also used 
!    in organize_data.f90
! V4_26        2012/12/06 Ulrich Blahak, Burkhardt Rockel, Ulrich Schaettler
!                         Oliver Fuhrer
!  Changed "ytrans_in /= '    ' " to "LEN_TRIM(ytrans_in) > 0"
!  Changes in NetCDF I/O: change some variables of dimension with length 1 to scalar
!  and reduce the number of dimension IDs (Burkhardt Rockel)
!  Read time level ntke for TKE scheme from binary restart file for correct
!   TKE restart and write TKE-restart data into corresponding time level (US)
!  Introduced call to get_free_unit for restart files again (which was erroneously 
!   moved outside to organize_data, because of prefetching of grib files) (US)
!  Correct grib table numbers for multi-layer snow variables (BR)
!  Adapted I/O short names for multi-layer snow model (US)
!  Correction of error message (by Oli Fuhrer)
! V4_27        2013/03/19 Michael Baldauf, Astrid Kerkweg, Ulrich Schaettler
!  Introduced p0ref in argument list to SR reference_atmosphere_BV
!  MESSy interface introduced: switch off input for online coupling (AK)
!  some initialisations added (AK)
! V4_28        2013/07/12 Ulrich Schaettler
!  Implemented grib_api for reading GRIB(1/2) data
!  In particular: read GRIB2 HHL-file for COSMO-grid, if necessary
!  Use subroutines and variables for vertical grid and reference atmospheres 
!    from module vgrid_refatm_utils
!  Adapted interface to read_gribapi with special grib_api integer
! V4_29        2013-10-02 Ulrich Schaettler, Astrid Kerkweg
!  Corrections for NetCDF: 
!   - put vertical coordinate parameters to pv_in for later use in the model
!   - allocate zvc_params in SR read_nc_gdefs in all tasks
!  Corrections for GRIB (setting .not. l_ke_in_gds):
!   - compute the correct value of ivctype in SR get_vertcoord
!   - if GRIB1 data are read and GRIB2 data written, a HHL-file need not be read
!  Bug fix when running with more tasks than GRIB records in input file:
!    not all tasks called the global communication then
!  Introduced namelist switch lan_w_so to choose W_SO from nudging or external
!    analysis; increased nanfld to 15.
!  For the COSMO-Model only use vcoord and refatm from vgrid_refatm_utils
!  Unification of MESSy interfaces and COSMO Tracer structure (AK)
! V4_30        2013/11/08 Ulrich Schaettler
!  Renamed ipds to ipds_in to reflect usage for input data
!  Removed check of leveltype for Flake variables
!  Initialized pv_in also for pure GRIB2 I/O
!  Initialized refatm%refatm_id also for binary output during restarts
! V5_1         2014-11-28 Ulrich Schaettler, Helmut Frank, Ulrich Blahak, Oliver Fuhrer
!  Check iscalval1 instead of ilevbot to decide if T_SO at 0 is required.
!  In case of grib_api, check vertical coordinate parameters only for 
!    hybrid, hybridLayer
!  When checking reference atmosphere parameters for GRIB2, check izexed==2
!   (and not igrbednr==2)
!  Modifications to read HHL and P from laf-file and to read vcoord from 
!   typeOfSecondFixedSurface from HHL
!  Removed reading of extra HHL-file
!  Modifications to read a reference profile, which is at the lowest grid point
!   available above sea level
!  In case of restart: change input directory of restart file from ydir_restart to ydir_restart_in
!  Replaced ireals by wp (working precision) (OF)
!  Removed option lexpl_lbc=.FALSE.  (MB)
!  Renamed QN-variables with NC... (US)
!  Do not recalculate hhl_prof and vcflat for idealized restart runs (lartif_data=.true.). (UB)
! V5_2         2015-05-21 Ulrich Blahak, Ulrich Schaettler
!  Added computation of kflat for grib2 input after the calls to reference_atmosphere_xxx().
!  Replaced imp_ireals by imp_integers in call to global_values for exchanging igpsl,jgpsl (US)
!  Only use nzlist_xxx as indices, if these values are really set.
! V5_3         2015-10-09 Ulrich Schaettler
!  From now on, HHL is written to restart files and taken from there, if available 
!   (is necessary for GRIB2 based simulations, because HHL cannot be reconstructed then)
!  If HHL is not available, it is re-computed
!  Moved computation of ndiff_ini_bd to organize_data, to do it for all cases once
!   (including restarts)
!  Fix computation of vcflat when reading GRIB2 data (had to initialize izexvc)
!  Fix in SR create_file_name in case of lbdana=.TRUE. and a dt, which does not always
!   fit to a full hour (US)
! V5_4         2016-03-10 Oliver Fuhrer, Pascal Spoerri, Xavier Lapillonne, Ulrich Schaettler
!  Added tracking of boundary fields in case of GPU compilation
!  Switched off date checking for HHL fields (US)
! V5_4a        2016-05-10 Ulrich Schaettler
!  Bug fix for reading NetCDF Data without HHL: lzchecklist must not be set then
! V5_4b        2016-07-12 Ulrich Schaettler, Jochen Foerstner
!  Compute field isoiltyp (for ICON version of TERRA, but can be used everywhere
!   in the COSMO-Model) (US)
!  Compute field llakemask (for Tiedtke-Bechtold convection) (JF)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Modules used:
#ifdef GRIBAPI
USE grib_api
#endif

USE data_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    iintegers, & ! KIND-type parameter for standard integer variables
    irealgrib, & ! KIND-type parameter for real variables in the grib library
    iwlength,  & ! length of an integer word in byte
    intgribf,  & ! KIND-type parameter for fortran files in the grib library
    intgribc,  & ! KIND-type parameter for C files in the grib library
    int_ga       ! integer precision for grib_api: length of message in bytes

USE data_fields, ONLY :  &
    fr_land,      & ! fraction of land in a grid element 
    fr_lake,      & ! fraction of lakes in a grid element 
    hsurf,        & ! height of surface topography
    plcov,        & ! fraction of plant cover                         --
    lai,          & ! leaf area index of plants                       --
    rootdp,       & ! depth of the roots                            ( m  )
    hhl,          & ! geometrical height of model half levels
    llandmask,    & ! landpoint mask 
    llakemask,    & ! lake point mask 
    soiltyp,      & ! type of the soil (keys 0-9)
    isoiltyp,     & ! type of the soil (keys 0-9) (in integers)
    ps,           & ! surface pressure
    pp,           & ! deviation from the reference pressure 
    pp_bd,        & ! deviation from the reference pressure 
    qrs,          & ! precipitation water (water loading)
    rho,          & ! total density of moist air                   (kg/m3)
    rho0,         & ! reference density at the full model levels   (kg/m3)
    dp0,          & ! pressure thickness of model layers
    p0,           & ! base-state pressure of full model levels
    p0hl,         & ! base-state pressure of half model levels
    t0,           & ! base-state temperature of full model levels
    t0hl,         & ! base-state temperature of half model levels
    rmy,          & ! Davis-parameter for boundary relaxation 
    t,            & ! temperature
    t_s,          & ! temperature of the ground surface
    t_so,         & ! multi-layer soil temperature
    t_g,          & ! weighted surface temperature
    t_snow,       & ! temperature of the snow-surface
    t_snow_mult,  & ! temperature of the snow-surface
    dzh_snow_mult,& !
    w_snow_mult,  & !
    wliq_snow,    & !
    rho_snow_mult,& ! snow density !_br 23.01.12
    w_g1,         & ! water content of the upper soil layer
    w_g1_bd,      & ! boundary field for w_g1
    w_g2,         & ! water content of the medium soil layer
    w_g2_bd,      & ! boundary field for w_g2
    w_g3,         & ! water content of the lower soil layer
    w_g3_bd,      & ! boundary field for w_g3
    w_snow,       & ! water content of snow
    t_ice,        & ! temperature of ice/water surface              (  K  )
    h_ice           ! lake/sea ice thickness                        (  m  )

USE data_modelconfig, ONLY : &
    dlon,         & ! grid point distance in zonal direction (in degrees)
    dlat,         & ! grid point distance in meridional direction (in degrees)
    dt,           & ! long time-step
    istart,       & ! start index for computations
    jstart,       & ! start index for computations
    iend,         & ! end index for computation
    iendu,        & ! end index for the forecast of u
    jend,         & ! end index for computations
    jendv,        & ! end index for the forecast of v
    istartpar,    & ! start index for computations in the parallel program
    jstartpar,    & ! start index for computations in the parallel program
    iendpar,      & ! end index for computations in the parallel program
    jendpar,      & ! end index for computations in the parallel program
    ie,           & ! number of grid points in zonal direction
    ie_tot,       & ! number of grid points in zonal direction total
    je,           & ! number of grid points in meridional direction
    je_tot,       & ! number of grid points in meridional direction total
    ie_max,       & ! Max. of ie on all processors
    je_max,       & ! Max. of je on all processors
    ke,           & ! number of grid points in vertical direction
    ke_tot,       & ! number of grid points in vertical direction total
    nlandpoints,  & ! number of land points in the grid
    nlandpoints_tot,  & ! number of land points in the total domain
    ke1,          & ! KE+1
    ke_soil,      & ! number of layers in the multi-layer soil model
    ke_snow,      & ! number of layers in the multi-layer snow model
    czmls,        & ! depth of the soil main layers in meters
    czhls,        & ! depth of the soil half layers in meters
    msoilgrib,    & ! grib coded depth of main soil levels in centimeters
    pollon,       & ! longitude of the rotated north pole (in degrees, E>0)
    pollat,       & ! latitude of the rotated north pole (in degrees, N>0)
    polgam,       & ! angle between the north poles of the systems
    startlon_tot, & ! transformed longitude of the lower left grid point
                    ! of the total domain (in degrees, E>0)
    startlat_tot    ! transformed latitude of the lower left grid point
                    ! of the total domain (in degrees, N>0)

USE data_modelconfig, ONLY : &
    hhl_prof,     & ! a special hhl-profile
    klv950,       & ! k index of the LM-mainlevel, on 950 HPa
    klv850,       & ! k index of the LM-mainlevel, on 850 HPa
    klv800,       & ! k index of the LM-mainlevel, on 800 HPa
    klv700,       & ! k index of the LM-mainlevel, on 700 HPa
    klv500,       & ! k index of the LM-mainlevel, on 500 HPa
    klv400,       & ! k index of the LM-mainlevel, on 400 HPa
    klv300,       & ! k index of the LM-mainlevel, on 300 HPa
    idt_qv, idt_qc, idt_qs,  idt_qg,  idt_qr,    &
    idt_qi, idt_qh, idt_qnc, idt_qnr, idt_qni,   &
    idt_qns, idt_qng, idt_qnh

USE data_constants, ONLY : &
    p0ref,        & ! reference pressure for Exner-function (Pa)
    r_d,          & ! gas constant for dry air
    g,            & ! gravity acceleration
    rvd_m_o         ! r_v/r_d - 1

USE data_runcontrol, ONLY : &
    idbg_level,   & ! to control the verbosity of debug output
    ldebug_io ,   & ! if .TRUE., debug output for I/O
    lprintdeb_all,& ! .TRUE.:  all tasks print debug output
                    ! .FALSE.: only task 0 prints debug output
    ltime,        & ! detailed timings of the program are given
    ndiff_ini_bd, & ! difference between start date and date of boundary data
    llm,          & ! if .TRUE., running with lowered upper boundary
    l2tls,        & ! forecast with 2-TL integration scheme
    itype_fast_waves,& ! Type of fast waves solver for Runge-Kutta dynamics
    lmulti_layer, & ! run multi-layer soil model
    lmulti_snow,  & ! run multi-layer snow model
    lseaice,      & ! forecast with sea ice model
    llake,        & ! forecast with lake model
    itype_gscp,   & ! type of grid-scale precipitation physics
    lradtopo,     & ! if .TRUE., calculate topographic correction of radiation
    nhori,        & ! number of sectors for the horizont array by the topographic
                    ! correction of the radiation
    lartif_data,  & ! forecast with self-defined artificial data
    itype_calendar,&! for specifying the calendar used
    itype_spubc,  & ! type of Rayleigh damping in the upper levels
    nnew,         & ! corresponds to ntstep + 1
    nnow,         & ! corresponds to ntstep
    nold,         & ! corresponds to ntstep - 1
    ntke,         & ! time level for TKE
    nstart,       & ! first time step of the forecast
    nstop,        & ! last time step of the forecast
    ntstep,       & ! actual time step
    lspubc,       & ! with Rayleigh damping in the upper levels
    nuspecif,     & ! output unit for protocolling the task
    psm0,         & ! initial value for mean surface pressure ps
    dsem0,        & ! initial value for mean dry static energy
    msem0,        & ! initial value for mean moist static energy
    kem0,         & ! initial value for mean kinetic energy
    qcm0,         & ! initial value for mean cloudwater content
    yakdat1         ! actual date (ydate_ini+ntstep/dt) in the form
                    ! yyyymmddhhmmss (year, month, day, hour, min, sec)

USE data_runcontrol, ONLY : &
    leps,         & ! switch ensemble mode on/off
    fac_plcov,    & ! modification factor for PLCOV
    rmin_plcov,   & ! lower limit of PLCOV
    rmax_plcov,   & ! upper limit of PLCOV
    fac_rootdp,   & ! modification factor for ROOTDP
    rmin_rootdp,  & ! lower limit of ROOTDP
    rmax_rootdp,  & ! upper limit of ROOTDP
    fac_lai,      & ! modification factor for LAI
    rmin_lai,     & ! lower limit of LAI
    rmax_lai        ! upper limit of LAI

USE data_soil,      ONLY : &
    cdzw12,  & !  thickness of upper soil water layer in two-layer model
    cdzw13,  & !  thickness of upper soil water layer in three-layer model
    cdzw23,  & !  thickness of middle soil water layer in three-layer model
    cf_snow    !  parameter for the calculation of the fractional snow coverage

USE data_io,        ONLY : &
  ydirini,           & ! directory of initial data
  ydirbd,            & ! directory of boundary data
  ydir_restart_in,   & ! directory for reading restart file
  ytunitbd,          & ! unit of time for the boundary data
  ydate_ini,         & ! start of the forecast 
  ydate_bd,          & ! start of the forecast for boundary data
  lbdclim,           & ! boundary data in climate model
  lbdsst,            & ! T_S boundary data are used only over the sea
                       ! (SST is not maintained constant during the integration)
  l_ke_in_input,     & ! indicates whether GRIB1 input data contains ke in meta data
  lanfld,            & ! contains switches for all fields to be checked for
                       ! time range indicator = 0
  lbdana,            & ! boundary data are analysed data
  lchkbd,            & ! checking the boundary data
  lchkini,           & ! checking the initial data
  nanfld,            & ! max. number of input fields to be checked for
                       ! time range indicator itri=0
  ytrans_in,         & ! directory for reading ready-files
  nincwait,          & ! if ready-file is not available wait nincwait seconds
                       ! until next attempt
  nmaxwait,          & ! if ready-file is not available after nmaxwait seconds,
  ynaman,            & ! name of fields to be checked for time range indicator
  nsma_stat,         & ! status for soil moisture analysis
  ymode_read,        & ! mode for opening the (read) Grib files
  yform_read,        & ! format for reading files
  nuchkdat,          & ! checking the I/O data
  yuchkdat,          & ! checking the I/O data
  lmmss,             & ! 10/14 digits date format
  npds,              & ! Dimension for product definition section (pds)
  ngds,              & ! Dimension for grid description section (gds)
  nbms,              & ! Dimension for bit map section (bms)
  nbds,              & ! Dimension for binary data section
  ndsup,             & ! Dimension for dsup
  ndims,             & ! Dimension for idims (contains all dimensions)
  ndims_id_in,       & ! number of dimensionID's for netCDF formatted input
  idims_id_in,       & ! for IDs of the dimensions of netCDF formatted input
  lfd, lbm, lds,     & ! Array sizes for GRIB variables
  nvar,              & ! maximum number of variables in LM variable table
  inrvert_in,        & ! number of vertical coordinate parameters of input data
  pv_in,             & ! array for vertical coordinate parameters of input data
  max_bd_fields,     & ! maximum number of boundary fields which can be stored
  num_bd_fields,     & ! number of boundary fields read
  bd_list              ! array of structures containing bd fields

USE data_io,        ONLY : &
  idwdednr,          & ! grib edition number for DWD library
  igrbednr,          & ! "working" grib edition number (to be set during run-time)
  undefgrib,         & ! value for "undefined" in the grib routines
  undefncdf,         & ! value for "undefined" in the NetCDF routines
  undef,             & ! the same as undefgrib but with other KIND-Parameter
  ngribednr,         & ! to store GRIB edition number for input
  nprocess_ini_in,   & ! process gener. identification for initial (analysis)
  nprocess_bd_in,    & ! and for boundary (forecasts) data from input data
  lbd_frame,         & ! if .TRUE., boundary data are on a frame
  ilevbotnoframe,    & ! bottom model level with b.d. defined on the whole grid
                       ! (model levels below are defined on a frame)
  npstrframe,        & ! width (number of points) of the strip around
                       ! the b.d. frame
  ytunit_restart,    & ! unit of timescale

! Global arrays
  ylevltypes1,       & ! to convert GRIB1 level types to grib_api string typeOfLevel for GRIB1
  ylevltypes2,       & ! to convert GRIB1 level types to grib_api string typeOfLevel for GRIB2
  ysteptypes,        & ! to convert GRIB1 time range indicator to grib_api stinf stepType
  rscalefac,         & ! Array to convert GRIB2 scale factors to real numbers
  itabletypes,       & ! Array to convert GRIB1 table types to the index
                       ! in structure lst_gribtabs
  iblock,            & ! array for gribed data
  idims_in,          & ! array for all dimensions
  ibmap,             & ! array for
  ipds_in,           & ! product definition section for input
  igds_in,           & ! grid description section
  ibms,              & ! bit map section
  ibds,              & ! binary data section
  dsup,              & ! Parameter for grib routines
  ds_grib,           & ! array for unpacked data
  ds_real,           & ! array for unpacked data

! Global types
  list_description,  & ! structure for list descriptions
  var                  ! array for LM variable table

USE data_parallel,      ONLY :  &
    nproc,          & ! total number of processors: nprocx * nprocy
    lasync_io,      & ! if .TRUE.: the model runs with extra PEs for
                      ! asynchronous IO
    ltime_barrier,  & ! if .TRUE.: use additional barriers
    my_world_id,    & ! rank of this subdomain in the global communicator
    my_cart_id,     & ! rank of this subdomain in the cartesian communicator
    my_cart_neigh,  & ! neighbors
    isubpos,        & ! positions of the subdomains in the total domain. Given
                      ! are the i- and the j-indices of the lower left and the
                      ! upper right grid point in the order
                      !                  i_ll, j_ll, i_ur, j_ur.
                      ! Only the interior of the domains are considered, not
                      ! the boundary lines.
    num_compute,    & ! number of compute PEs
    nboundlines,    & ! number of boundary lines of the domain for which
                      ! no forecast is computed = overlapping boundary
                      ! lines of the subdomains
    icomm_cart,     & ! communicator that belongs to the cartesian grid 
    imp_reals,      & ! determines the correct REAL type used in the model
                      ! for MPI
    imp_grib,       & ! determines the REAL type for the GRIB library
    imp_integers,   & ! determines the correct INTEGER type used in the model
                      ! for MPI
    imp_character,  & ! determines the correct CHARACTER type used in the model
                      ! for MPI
    imp_logical       ! determines the correct LOGICAL   type used in the
                      ! model for MPI

!------------------------------------------------------------------------------

USE utilities,          ONLY: get_utc_date, sleve_split_oro, diff_minutes, &
                              dosleep

USE environment,        ONLY: model_abort, comm_barrier, get_free_unit,    &
                              release_unit

USE parallel_utilities, ONLY: distribute_values, scatter_values,           &
                              global_values, gather_field, ij_local,       &
                              distribute_field

USE meteo_utilities,    ONLY: tgcom, calrho, calps

USE io_utilities,       ONLY: open_file, read_grib, read_gribapi,          &
                              read_netcdf, read_restart, close_file,       &
                              make_fn, check_input_grid, check_record,     &
                              compute_grib_intbuffer_length

USE time_utilities,     ONLY: get_timings, i_read_data, i_distribute_data, &
                              i_computations_I, i_initializations,         &
                              i_meta_data_r

USE vgrid_refatm_utils, ONLY:                                              &
    reference_atmosphere, reference_atmosphere_2,  lanalyt_calc_t0p0,      &
    vcoord, refatm, rundefined, nfltvc, svc1, svc2,                        &
    reference_atmosphere_BVconst, k_index_of_pressure_levels,              &
    vcoord_type, refatm_defaults, vcoord_defaults, imax_refatmtype,        &
    imax_vcoordtype, uuid_in_string, uuid_2char, lhhl_hasbeenread,         &
    uuid_create, uuid_out, uuid_out_string
 
USE src_tracer,         ONLY: trcr_get, trcr_errorstr
USE data_tracer,        ONLY: T_ERR_NOTFOUND

!------------------------------------------------------------------------------

#ifdef NETCDF
USE netcdf,           ONLY :   &
  nf90_enotvar,            &
  nf90_enotatt,            &
  nf90_get_att,            &
  nf90_get_var,            &
  nf90_inq_dimid,          &
  nf90_inq_varid,          &
  nf90_inquire,            &
  nf90_inquire_dimension,  &
  nf90_noerr,              &
  nf90_strerror
#endif

#ifdef TWOMOM_SB
USE src_twomom_sb_interface, ONLY:  &
     set_qnc_from_qc_sb, set_qnr_from_qr_sb, set_qni_from_qi_sb, &
     set_qns_from_qs_sb, set_qng_from_qg_sb
#endif

#ifdef MESSY
USE messy_main_data_bi,       ONLY: L_IS_CLIENT
USE messy_main_channel_bi,    ONLY: messy_channel_read_restart
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================

! string variable to hold grid information
  CHARACTER (LEN=200)      :: grid_mapping

! indices of grid point with orography nearest to sea-level 
! (or sea level point, if it exists in the model domain)
! to get a reference hhl-profile closest to sea level
  INTEGER(KIND=iintegers)  :: igpsl, jgpsl

!==============================================================================

PUBLIC  organize_input, create_file_name

PRIVATE bd_from_id, id_mix_bd, get_vertcoord, fill_procarray,              &
        scatter_data, read_nc_vdefs

PUBLIC :: read_nc_gdefs      ! needed for Messy

!==============================================================================
! Module procedures
!==============================================================================

CONTAINS

!==============================================================================
!+ Driver routine for the input of Grib files
!------------------------------------------------------------------------------

SUBROUTINE organize_input (lopen_file, nufile, ydata, nlgw_lm, nlgw_input, &
                           listin, nvarin, yformat, yextension, lfirst,    &
                           ntlev, nbdstep, lallow_missing, lchecklist)

!------------------------------------------------------------------------------
!
! Description:
!  organize_input is the driver routine for reading grib files. Whether initial
!  or boundary data are read is determined by the parameter ydata (either
!  `initial' or `boundary'). The parameters lfirst, ntlev and nbdstep are 
!  optional but have to be present when reading boundary data.
!
!  The organization of this routine is such, that it works on parallel and
!  on sequential platforms. The strategy has been developed by Rainer Johanni
!  from SGI/CrayResearch.
!
! Method:
!  After creating the filename and opening the grib file, all records are read
!  in `read_loop' (endless loop over records until all records are read).
!  Every processor gets a total record for de-gribing (up to the maximum number
!  of records or processors). 
!  In `distribute_loop' (loop over processors), the records are distributed
!  according to the domain decomposition.
!  Depending on which data are read (`initial' or `boundary') additional
!  actions are taken.
!
! Input files:
!  Grib-files with initial or boundary data. The name of the files depends on 
!  the date of the forecast. The files are read record by record with the
!  routine read_grib.
!
!------------------------------------------------------------------------------

! Scalar arguments with intent(in):
  LOGICAL                 , INTENT(IN)  ::                      &
    lopen_file ! Flag. If true, was opened before (due to pre-fetching)
  
  INTEGER (KIND=iintegers), INTENT(INOUT) ::  &
    nufile

  CHARACTER (LEN=*)       , INTENT(IN)  ::                      &
    ydata          ! determines whether `initial' or `boundary' data

  INTEGER (KIND=iintegers), INTENT(IN)  ::                      &
    nvarin,      & ! number of variables in input list
    nlgw_lm,     & ! number of ground water levels in the LM 
    nlgw_input     ! number of ground water levels in the input file

  TYPE(list_description)  , INTENT(IN)  ::                      &
    listin(nvarin) ! List of fields for reading in

  CHARACTER (LEN= 4),       INTENT(IN)  ::    &
    yformat        ! determines the format that has to be read

  CHARACTER (LEN=1),        INTENT(IN)     ::    &
    yextension     ! indicates actual restart-file

  LOGICAL,                  INTENT(IN), OPTIONAL  ::            &
    lfirst         ! if .TRUE., first boundary data set

  INTEGER (KIND=iintegers), INTENT(IN), OPTIONAL  ::            &
    ntlev,       & ! time level for boundary data
    nbdstep        ! time step to create the file name

! CK 20101117 we need a way to read data without check if
! everything has been read.
  LOGICAL, INTENT(IN), OPTIONAL  ::                             &
    lallow_missing ! if set, it will omit the check if everything has been read.
  LOGICAL, INTENT(OUT), OPTIONAL ::                             &
    lchecklist(0:ke1, nvarin) ! if set, will return a copy of lzchecklist
                              ! to the calling routine

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER  (KIND=intgribf)   ::  &
! iniyy , ibdyy  ,& ! year     / of the initial
! inimm , ibdmm  ,& ! month   /  date and time
! inidd , ibddd  ,& ! day    <   of the runs providing
! inihh , ibdhh  ,& ! hour    \  the initial state
! inimin, ibdmin ,& ! minute   \ ctd.
! inisec, ibdsec, & ! second    \ ctd.
! imindif        ,& ! time difference [min] between these two runs
  ieef           ,& ! 
  ierrf             ! error status for grib routines

INTEGER  (KIND=iintegers)  ::  &
  ierr, ierrstat, niostat,     & ! status and error status variables
  j1, j2, js, nztlev, nzsteps, izlocsl, jzlocsl, itask!,          &
! ndiff_ini_bdref   ! corrects ndiff_ini_bd if ydate_ini not equal to
!                   ! reference time of boundary data 

INTEGER  (KIND=iintegers)  ::  &
  iz_rsize, izloc1, izloc2, izloc3, iz_rank, iz_dim3, iz_loclist, nzbytes, &
                                 ! characteristics of read records
  iz_lfd, iz_ngds                ! for compatibility

INTEGER  (KIND=iintegers)  ::  &
  iz_countl,                   & ! counter of required records
  izvar_count, izlev_count,    & ! counter for NetCDF variables and levels
  myzvar, myzlev, myzlevtot,   & ! organization indices returned by write_netcdf
  nzlist_p, nzlist_pp,         & ! to save the index of p, pp
  nzlist_hhl,                  & ! to save the index of hhl
  nzlist_t_so, nzlist_t_s,     & ! to save the index of t_s, t_so in the list
  nzlist_frsn,                 & ! to save the index of freshsnow
  nzlist_tsnow,                & ! to save the index of t_snow
  nzlist_tsnow_mult,           & ! to save the index of t_snow_mult
  nzlist_wtot, nzlist_wliq,    & ! to save the index of w_snow_mult and wliq_snow
  nzlist_dzh,                  & ! to save the index of dzh_snow_mult
  nzlist_rsnow_mult,           & ! to save the index of rho_snow_mult   !_br 23.01.12
  izvctype_read                  ! check vertical coordinate type in restarts

INTEGER  (KIND=iintegers)  ::  &
  nzday           ,& ! century day of the current date
  istat,iproc, nland, izerror,        &
  izdebug, i, j, ip, ie_p, je_p, it, jt, ij,        &  ! additional variables
  k, ksn, izexch, igribid, ireturn, igriblen,       &
  ilevtyp, ilevel, ilevtop, ilevbot, izpv,          &
  igenproc, izexed,                                 &
  iscalval, iscalfac, ij_min, izexvc

INTEGER (KIND=int_ga)    ::  &
  izmaxlen, iz_rsize_ga

INTEGER (KIND=iintegers)   :: &
  ivar_id(nvarin)    ! variable IDs of NetCDF

INTEGER (KIND=iintegers)   ::  &
  intbuf      (ngds) ! an integer buffer for sending

CHARACTER (LEN=100)        ::  &
  charbuf

REAL (KIND=irealgrib)      ::  &
  rr_min, zgribarray (ie_max,je_max,0:num_compute-1)

REAL (KIND=wp)             ::  &
  zrealarray (ie_max,je_max,0:num_compute-1)

REAL (KIND=wp)             ::  &
  zfac1, zfac2, zhour, zrbuf(ke1+10), &
  zvc_params(ke1), zvck    ! Additional variables

REAL (KIND=irealgrib)      ::  &
#ifdef GRIBDWD
  REFSTF,    & ! from DWD Grib library
#endif
  zundef       ! to pass undef-value to subroutines

LOGICAL                    ::  &
  lrequired,  & ! indicates whether a record from the grib file is required
  leof,       & ! indicates the end of file
  lgot_vcoord,& ! indicates that vertical coordinate parameters have been read
  lgot_vparam,& ! indicates that vertical coordinate parameters have been read
  lgot_vcuuid,& ! indicates that a uuid has already been read
  lgot_refatm,& ! indicates that reference atmosphere parameters have been read
  lzcomp_pp,  & ! if pp has to be computed when full pressure is read
  lzcomp_hhl, & ! to re-compute hhl in reference_atmosphere_xx
  lzexch_ednr,& ! to exchange grib edition number
  lzdate,     & ! to determine, whether date has to be checked
  linit,      & ! for the first cycle of read_loop
  lwait,      & ! flag to signal waiting for ready files
  lgexist,    & ! for inquiring grib files
  lwb3r         ! indicates whether three soil water levels are present

TYPE(list_description)  :: list_hhl(1)
INTEGER(KIND=iintegers) :: nlist_hhl = 1

LOGICAL                    ::  &
  lzchecklist (0:ke1, nvarin), & ! list for checking which variables have been read
  lzcheck_hhl (0:ke1,      1)    ! list for checking which variables have been read

CHARACTER (LEN=250)        ::  &
  yname      ! name of the grib file (is determined in make_fn)
CHARACTER (LEN=260)        ::  &
  yready     ! name of the ready file
CHARACTER (LEN=260)        ::  &
  yready2     !name of second ready file
CHARACTER (LEN=  3)        ::  &
  yzhead     ! header of grib-file name
CHARACTER (LEN= 25)        ::  &
  yroutine   ! name of this routine for error handling
CHARACTER (LEN= 28)        ::  &
  ydum       ! for date-checking
CHARACTER (LEN= 14)        ::  &
  ydatchk    ! for date-checking
CHARACTER (LEN= 14)        ::  &
  ydatchk2   ! for date-checking
CHARACTER (LEN=255)        ::  &
  yerrmsg    ! error message for error handling
CHARACTER (LEN= 30)        ::  &
  yzmyname, ytyofle     ! type of level from grib_api

! local fields for the SLEVE coordinate
REAL (KIND=wp)                          ::  &
  zhsurfs    (ie,je,2)  ! height of splitted topography parts

! local field needed for the computation of t_g
REAL (KIND=wp)                          ::  &
  zt_s(ie,je)           ! = t_s   on land and sea
                        ! = t_ice on sea ice (if present)

REAL (KIND=wp),     ALLOCATABLE         ::  &
  zhsurfs_tot(:,:,:), & ! height of splitted topography parts of full domain
  zhsurf_tot (:,:)      ! topography of full domain

REAL (KIND=wp),     POINTER :: &
  qv  (:,:,:)=> NULL(),    & ! QV at tlev=nnew
  qc  (:,:,:)=> NULL(),    & ! QC at tlev=nnew
  qr  (:,:,:)=> NULL(),    & ! QR at tlev=nnew
  qi  (:,:,:)=> NULL(),    & ! QI at tlev=nnew
  qs  (:,:,:)=> NULL(),    & ! QS at tlev=nnew
  qg  (:,:,:)=> NULL()       ! QG at tlev=nnew

#ifdef TWOMOM_SB

REAL (KIND=wp),     POINTER :: & 
  qh  (:,:,:)=> NULL(),    &
  qnc  (:,:,:)=> NULL(),    &
  qnr  (:,:,:)=> NULL(),    &
  qni  (:,:,:)=> NULL(),    &
  qns  (:,:,:)=> NULL(),    &
  qng  (:,:,:)=> NULL(),    &
  qnh  (:,:,:)=> NULL(),    &
  qnc_now  (:,:,:)=> NULL(),    &
  qnr_now  (:,:,:)=> NULL(),    &
  qni_now  (:,:,:)=> NULL(),    &
  qns_now  (:,:,:)=> NULL(),    &
  qng_now  (:,:,:)=> NULL(),    &
  qnh_now  (:,:,:)=> NULL(),    &
  qv_bd  (:,:,:,:)=> NULL(),    &
  qc_bd  (:,:,:,:)=> NULL(),    &
  qr_bd  (:,:,:,:)=> NULL(),    &
  qi_bd  (:,:,:,:)=> NULL(),    &
  qs_bd  (:,:,:,:)=> NULL(),    &
  qg_bd  (:,:,:,:)=> NULL(),    &
  qh_bd  (:,:,:,:)=> NULL(),    &
  qnc_bd  (:,:,:,:)=> NULL(),    &
  qnr_bd  (:,:,:,:)=> NULL(),    &
  qni_bd  (:,:,:,:)=> NULL(),    &
  qns_bd  (:,:,:,:)=> NULL(),    &
  qng_bd  (:,:,:,:)=> NULL(),    &
  qnh_bd  (:,:,:,:)=> NULL()

REAL (KIND=wp)             :: zqncmax, zqcmax

REAL (KIND=wp)             ::  &
     zrho(ie,je,ke)  ! for diagnosis of NCXXXX

LOGICAL                    ::  &
  lzcalc_qncloud, lzcalc_qnrain, lzcalc_qnsnow, lzcalc_qnice,  &
  lzcalc_qngraup, lzcalc_qnhail
#endif

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  ! Initialize, whether additional debug output shall be done
  IF (ldebug_io) THEN
    IF (lprintdeb_all) THEN
      izdebug = idbg_level
    ELSE
      IF (my_cart_id == 0) THEN
        izdebug = idbg_level
      ELSE
        izdebug = 0
      ENDIF
    ENDIF
  ELSE
    izdebug = 0
  ENDIF

  IF (izdebug >= 10) PRINT *, '  ENTERING:  organize_input  ', ydata

  yroutine    = 'organize_input'
  istat       =  0
  izerror     =  0
  ireturn     =  0
  izexed      =  0
  izexvc      = 1000000 ! a very high number for height in meters
  igpsl       = -1
  jgpsl       = -1  ! will be set only in the task that reads HSURF

  IF ((ydata == 'initial') .OR. (ydata == 'restart')) THEN
    lgot_vcoord = .FALSE.
    lgot_vparam = .FALSE.
    lgot_vcuuid = .FALSE.
    lgot_refatm = .FALSE.
    uuid_in_string='xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx'
  ELSE
    lgot_vcoord = .TRUE.
    lgot_vparam = .TRUE.
    lgot_vcuuid = (uuid_in_string /= 'xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx')
    lgot_refatm = .TRUE.
  ENDIF

  IF ( (yformat == 'grb1') .OR. (yformat == 'apix') .OR. (yformat == 'bina') ) THEN
    lzexch_ednr = .TRUE.
  ELSE
    lzexch_ednr = .FALSE.
  ENDIF
  IF ( (.NOT. lmulti_layer) .AND. (nlgw_lm == 3 .OR. nlgw_input == 3) ) THEN
    lwb3r = .TRUE.
  ELSE
    lwb3r = .FALSE.
  ENDIF

  ! Preparations to process P and HHL
  lzcomp_pp   = .FALSE.
  lzcomp_hhl  = .TRUE.
  nlist_hhl   = 1
  list_hhl%name  = 'HHL       '
  list_hhl%iloc1 =  2
  list_hhl%iloc2 =  8
  list_hhl%iloc3 =  1
  list_hhl%idimvert = ke+1

  lzcheck_hhl (:,:) = .FALSE.
  lzchecklist (:,:) = .FALSE.
  ! CK 20101117 we also need to "falsify" the output checklist if appropriate
  IF (present(lchecklist)) THEN
    lchecklist(:,:) = .FALSE.
  ENDIF 
  
  ! when reading emissions or other ART fields some checks fail, because not
  ! all fields are present in these files.
  ! Work around: set the indizes here to default -1,
  ! then check only for their existence if they have
  ! been found in the variables requested (listin).
  nzlist_t_s        = -1
  nzlist_t_so       = -1
  nzlist_frsn       = -1
  nzlist_tsnow      = -1
  nzlist_tsnow_mult = -1
  nzlist_wtot       = -1
  nzlist_wliq       = -1
  nzlist_dzh        = -1
  nzlist_rsnow_mult = -1
  nzlist_p          = -1
  nzlist_pp         = -1
  nzlist_hhl        = -1

  DO j2 = 1, nvarin
    IF (listin(j2)%name == 'PP        ') nzlist_pp         = j2
    IF (listin(j2)%name == 'P         ') nzlist_p          = j2
    IF (listin(j2)%name == 'HHL       ') nzlist_hhl        = j2

    IF (lmulti_layer) THEN
      IF (listin(j2)%name == 'T_S       ') nzlist_t_s        = j2
      IF (listin(j2)%name == 'FRESHSNW  ') nzlist_frsn       = j2
      IF (listin(j2)%name == 'T_SNOW    ') nzlist_tsnow      = j2
      IF (listin(j2)%name == 'T_SNOW_M  ') nzlist_tsnow_mult = j2
      IF (listin(j2)%name == 'H_SNOW_M  ') nzlist_dzh        = j2
      IF (listin(j2)%name == 'W_SNOW_M  ') nzlist_wtot       = j2
      IF (listin(j2)%name == 'WLIQ_SNOW ') nzlist_wliq       = j2
      IF (listin(j2)%name == 'RHO_SNOW_M') nzlist_rsnow_mult = j2

      IF (listin(j2)%name == 'T_SO      ') THEN
        nzlist_t_so = j2
        DO j1 = 0, listin(j2)%idimvert
          lzchecklist (j1,j2) = .FALSE.
        ENDDO
      ELSE
        lzchecklist(0,j2) = .TRUE.
        DO j1 = 1, listin(j2)%idimvert
          lzchecklist (j1,j2) = .FALSE.
        ENDDO
      ENDIF
    ELSE
      lzchecklist(0,j2) = .TRUE.
      DO j1 = 1, listin(j2)%idimvert
        lzchecklist (j1,j2) = .FALSE.
      ENDDO
    ENDIF
    DO j1 = listin(j2)%idimvert+1, ke1
      ! nothing has to be read on these levels
      lzchecklist (j1,j2) = .TRUE.
    ENDDO
    IF ( (listin(j2)%name == 'TKVH      ') .OR.                         &
         (listin(j2)%name == 'TKVM      ') ) THEN
      ! these variables are defined from level 2 onward
      lzchecklist (1,j2) = .TRUE.
    ENDIF
  ENDDO

#ifdef I2CINC 
  ! SKIP READ-IN-PROCEDURE IN CASE OF I2CINC FOR 'initial' and 'boundary'
  IF ((ydata /= 'initial' .AND. ydata /= 'boundary')  &
       .OR. (.NOT. L_IS_CLIENT)) THEN
#endif

  ! Set lfd, lds and lbm
  lbm = 25000                       ! set by Lucio for working with frames
  lds = ie_tot * je_tot
  nzbytes = 8                       ! to be on the safe side
  lfd = compute_grib_intbuffer_length(ie_tot, je_tot, nzbytes, iwlength)
    ! should be large enough also for higher packing rates

  ! Allocate GRIB arrays
  iz_ngds   = INT (ngds, iintegers)
  iz_lfd    = INT (lfd , iintegers)

  ! the value that has to be given to the subroutine check_record depends
  ! on the input format, because this determines which value is used for
  ! undefined points
  IF ((yformat == 'grb1') .OR. (yformat == 'apix') .OR. (yformat == 'bina')) THEN
    zundef      =       undefgrib
    undef       = REAL (undefgrib, wp)
  ELSEIF (yformat == 'ncdf') THEN
    zundef      =       undefncdf
    undef       = REAL (undefncdf, wp)
  ENDIF

  ALLOCATE (iblock(lfd), ibmap(lbm), STAT=istat)
  ALLOCATE (ds_grib(lds), ds_real(lds), dsup(ndsup), STAT=istat)

  ! dimensions for grib routines
  idims_in( 1) = npds
  idims_in( 2) = ngds
  idims_in( 3) = nbms
  idims_in( 4) = nbds
  idims_in( 5) = lbm
  idims_in( 6) = ndsup
  idims_in( 7) = lds
  idims_in( 8) = lfd
  idims_in(9:20) = 0

  IF (ydata == 'boundary') THEN
    ! check, whether optional parameters are present for boundary data
    IF (.NOT. (PRESENT(lfirst)) .OR. .NOT. (PRESENT(ntlev)) .OR.        &
        .NOT. (PRESENT(nbdstep)) ) THEN
      yerrmsg = 'Optional parameters not present for boundary data'
      CALL model_abort (my_cart_id, 2013, yerrmsg, yroutine)
    ENDIF

    ! if analysed data are used, the initial data are copied to the 
    ! boundary data
    IF((lfirst .OR. nstart == nstop) .AND. lbdana) THEN
       CALL bd_from_id (listin, nvarin)
       ! Deallocate arrays for IO
       DEALLOCATE (iblock, ibmap, ds_grib, ds_real, dsup)
       PRINT *, 'leaving organize_input because of lbdana = ', lbdana
       RETURN
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 1.2: Check ready files
!------------------------------------------------------------------------------

  ! by default we don't wait for ready files
  lwait = .FALSE.

  ! If required, generate names of ready files
  IF ( (LEN_TRIM(ytrans_in) > 0) .AND. (ydata /= 'restart') ) THEN

    yready = ''
    yready2 = ''

    IF (ydata == 'boundary') THEN
      IF (lfirst .AND. lbdana) THEN
        ! first boundary data are taken from initial data and
        ! nothing has to be read (and waited for)
        lwait = .FALSE.
      ELSE
        lwait = .TRUE.
      ENDIF
    ELSE
      lwait = .TRUE.
    ENDIF

    IF (ydata == 'initial') THEN
      ! Create the filename for initial files: LMA_date
      yzhead = 'LMA'
      nzsteps = 0
    ELSE
      ! Create the filename for boundary files: LMB_forecasttime
      yzhead = 'LMB'
      nzsteps = nbdstep+ndiff_ini_bd
    ENDIF

    ! Create file name with a 14-digit date
    CALL make_fn (yzhead, yakdat1, ydate_ini, 'f', yextension, nzsteps, dt,    &
             .TRUE., itype_calendar, ytrans_in, yready, .TRUE., izdebug, ierr)

    ! For the initial file there could also be a ready file with a 10-digit date
    IF (yzhead == 'LMA' .AND. yakdat1(11:14) == '0000') THEN
      CALL make_fn (yzhead, yakdat1, ydate_ini, 'f', yextension, nzsteps, dt, &
                    .TRUE., itype_calendar, ytrans_in, yready2, .FALSE., izdebug, ierr)
    ENDIF

  ENDIF

  IF (lwait) THEN
    IF (my_cart_id == 0 .AND. izdebug >= 10) PRINT *, '  CHECKING ready files '
    CALL wait_for_file(yready, yready2)
  END IF

#ifdef I2CINC
  ENDIF
  IF (yformat /= 'ncdf') THEN
    undef     = REAL(undefgrib, wp)
  ELSE
    undef     = REAL(undefncdf, wp)
  ENDIF
#endif

!------------------------------------------------------------------------------
! Section 1.3: Create file name and open YUCHKDAT
!------------------------------------------------------------------------------

#ifdef I2CINC
  ! DO NOT OPEN INITIAL OR BOUNDARY FILES !
  ! DATA IS MADE AVAILABLE BY MESSYMMD
  IF ((ydata /= 'initial' .AND. ydata /= 'boundary')  &
       .OR. (.NOT. L_IS_CLIENT))  THEN
#endif

  IF (izdebug >= 10) PRINT *, '  CREATE file names'
  CALL create_file_name (nbdstep, ydata, yextension, yformat, yname, ydatchk2)

#ifdef I2CINC
  ENDIF
#endif

  ! open file YUCHKDAT and print a headline
  IF ( ((lchkini) .OR. (lchkbd)) .AND. (my_cart_id == 0) ) THEN
    IF (izdebug >= 10) PRINT *, '  OPEN YUCHKDAT'
    OPEN(nuchkdat, FILE=yuchkdat, FORM=  'FORMATTED', STATUS='UNKNOWN',  &
                   POSITION='APPEND', IOSTAT=niostat)
    IF(niostat /= 0) THEN
      yerrmsg = ' ERROR    *** Error while opening file YUCHKDAT *** '
      ierrstat = 2005
      CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
    ENDIF
  
    ! print the headline for this grib file
    IF ( (lchkini) .AND. (ydata == 'initial') ) THEN
      WRITE (nuchkdat,'(A)') 'Check the data in initial file: '
      WRITE (nuchkdat,'(A,A)')                                             &
            '    File:   ',TRIM(yname)
      WRITE (nuchkdat,'(A,I5,A,I5,A,I5)')                                  &
            '    ie_tot =',ie_tot,'   je_tot =',je_tot,'   ke_tot =',ke_tot
      WRITE (nuchkdat,'(A)') '    '
      WRITE (nuchkdat,'(A,A)')                                             &
            '     var       ee    lev         min      ',                  &
            'imin   jmin          max      imax   jmax         mean  '
    ELSEIF ( (lchkbd) .AND. (ydata == 'boundary') ) THEN
      WRITE (nuchkdat,'(A,I7)')                                            &
                    'Check the data in boundary file for step: ', nbdstep
      WRITE (nuchkdat,'(A,A)')                                             &
            '    File:   ',TRIM(yname)
      WRITE (nuchkdat,'(A,I5,A,I5,A,I5)')                                  &
            '    ie_tot =',ie_tot,'   je_tot =',je_tot,'   ke_tot =',ke_tot
      WRITE (nuchkdat,'(A)') '    '
      WRITE (nuchkdat,'(A,A)')                                             &
            '     var       ee    lev         min      ',                  &
            'imin   jmin          max      imax   jmax         mean  '
    ELSEIF (ydata == 'restart') THEN
      WRITE (nuchkdat,'(A,A)')                                             &
                    'Check the data in restart  file for:   ', yextension
      WRITE (nuchkdat,'(A,A)')                                             &
            '    File:   ',TRIM(yname)
      WRITE (nuchkdat,'(A,I5,A,I5,A,I5)')                                  &
            '    ie_tot =',ie_tot,'   je_tot =',je_tot,'   ke_tot =',ke_tot
      WRITE (nuchkdat,'(A)') '    '
      WRITE (nuchkdat,'(A,A)')                                             &
            '     var       ee    lev         min      ',                  &
            'imin   jmin          max      imax   jmax         mean  '
    ENDIF
  ENDIF

#ifdef I2CINC
  ! DO NOT OPEN INITIAL OR BOUNDARY FILES !
  ! DATA IS MADE AVAILABLE BY MESSYMMD
  IF ((ydata /= 'initial' .AND. ydata /= 'boundary')  &
       .OR. (.NOT. L_IS_CLIENT))  THEN
#endif

#ifdef MESSY
  ! DO NOT OPEN COSMO RESTART FILES
  ! RESTART IS DONE VIA CHANNEL IN MESSY
  IF ( yextension =='o' .OR. yextension == 'n') THEN
     IF (yextension =='o')  CALL messy_channel_read_restart('COSMO_ORI')
  ELSE
#endif

  ! initialize level counter and logical flags
  iz_countl  = 0
  leof       = .FALSE.
  linit      = .TRUE.

!------------------------------------------------------------------------------
! Section 2: Opening the file
!------------------------------------------------------------------------------
 
  IF (izdebug >= 10) PRINT *, '  OPEN file '

  ! open file, if this has not already been done during prefetching
  IF (.NOT. lopen_file) THEN
    IF (izdebug >= 10) PRINT *, '  OPEN file '

    IF (yformat == 'bina') THEN
      ! get a free unit-number for Fortran OPEN for restart-files
      ! (this has not been done before, because there is no prefetching of restart files)
      CALL get_free_unit (nufile)
    ENDIF

    CALL open_file (nufile, yname, ymode_read, yformat, icomm_cart,    &
                    my_cart_id, num_compute, lasync_io, idbg_level,    &
                    yerrmsg, ierr)
    IF (ierr /= 0) THEN
      CALL model_abort (my_cart_id, 2014, yerrmsg, yroutine)
    ENDIF
  END IF

  IF ( (yformat == 'bina') .AND. (yextension == 'o')) THEN
    ! Get the first 2 (or 3) records containing the information about the
    ! vertical coordinate
    IF (my_cart_id == 0) THEN
      ! read the initial values for the meanvalues
      READ (nufile,IOSTAT=niostat) psm0, dsem0, msem0, kem0, qcm0, ntke

      ! read the vertical coordinate parameters
      ! for restart, still GRIB1 conventions are used
      READ (nufile,IOSTAT=niostat) izvctype_read, refatm%p0sl,   refatm%t0sl,   &
                                 refatm%dt0lp, vcoord%vcflat, zvc_params
      IF     ( (izvctype_read >   0) .AND. (izvctype_read <= 100) ) THEN
        refatm%irefatm = 1
        vcoord%ivctype = izvctype_read
      ELSEIF ( (izvctype_read > 100) .AND. (izvctype_read <= 200) ) THEN
        refatm%irefatm = 2
        vcoord%ivctype = izvctype_read - 100
        READ (nufile,IOSTAT=niostat) refatm%delta_t, refatm%h_scal
      ELSEIF ( (izvctype_read > 200) .AND. (izvctype_read <= 300) ) THEN
        refatm%irefatm = 3
        vcoord%ivctype = izvctype_read - 200
        READ (nufile,IOSTAT=niostat) refatm%bvref
      ENDIF
      vcoord%nlevels           = ke1
      IF     (vcoord%ivctype == 1) THEN
        vcoord%sigm_coord(1:ke1) = zvc_params(1:ke1)
      ELSEIF ( ANY( vcoord%ivctype == (/2,3,4/) ) ) THEN
        vcoord%vert_coord(1:ke1) = zvc_params(1:ke1)
      ENDIF
    ENDIF

    ! distribute vertical coordinate parameters (still in GRIB1 meta data format)
    IF (my_cart_id == 0) THEN
      intbuf(1) = vcoord%ivctype
      intbuf(2) = refatm%irefatm
      intbuf(3) = izvctype_read     ! for convenience: then we do not need to recompute
      intbuf(4) = ntke
      zrbuf (1) = refatm%p0sl
      zrbuf (2) = refatm%t0sl
      zrbuf (3) = refatm%dt0lp
      zrbuf (4) = vcoord%vcflat
      DO k = 1, ke1
        zrbuf (4+k) = zvc_params(k)
      ENDDO
      inrvert_in = 4+ke1
      IF (refatm%irefatm == 2) THEN
        zrbuf (4+ke1+1) = refatm%delta_t
        zrbuf (4+ke1+2) = refatm%h_scal
        inrvert_in = 4+ke1+2
      ENDIF
      IF (refatm%irefatm == 3) THEN
        zrbuf (4+ke1+1) = refatm%bvref
        inrvert_in = 4+ke1+1
      ENDIF
    ENDIF

    IF (num_compute > 1) THEN
      CALL distribute_values (intbuf, 4, 0, imp_integers, icomm_cart, istat)
      CALL distribute_values (zrbuf , ke1+10, 0, imp_reals, icomm_cart, istat)
    ENDIF

    ! Now fill the vcoord and refatm-type, but also pv_in, inrvert_in for later use
    ! do it in the new style of gds coding
    ! (is also done by my_cart_id == 0, to fill pv_in there
    vcoord%ivctype   = intbuf(1)
    refatm%irefatm   = intbuf(2)
    izvctype_read    = intbuf(3)
    ntke             = intbuf(4)
    refatm%p0sl      = zrbuf (1)
    refatm%t0sl      = zrbuf (2)
    refatm%dt0lp     = zrbuf (3)
    vcoord%vcflat    = zrbuf (4)

    pv_in( 1) = REAL (izvctype_read, wp)    ! has been read above
    pv_in( 2) = REAL (ke, wp)
    pv_in( 3) = refatm%p0sl
    pv_in( 4) = refatm%t0sl
    pv_in( 5) = refatm%dt0lp
    pv_in( 6) = vcoord%vcflat

    vcoord%nlevels           = ke1
    IF     (vcoord%ivctype == 1) THEN
      DO k = 1, ke1
        vcoord%sigm_coord(k) = zrbuf (4+k)
        pv_in(6+k)           = zrbuf (4+k)
      ENDDO
    ELSEIF ( ANY( vcoord%ivctype == (/2,3,4/) ) ) THEN
      DO k = 1, ke1
        vcoord%vert_coord(k) = zrbuf (4+k)
        pv_in(6+k)           = zrbuf (4+k)
      ENDDO
    ENDIF
    inrvert_in = 6+ke1
    IF (refatm%irefatm == 2) THEN
      refatm%delta_t = zrbuf (4+ke1+1)
      refatm%h_scal  = zrbuf (4+ke1+2)
      pv_in(6+ke1+1) = 0.0_wp
      pv_in(6+ke1+2) = 0.0_wp
      pv_in(6+ke1+3) = 0.0_wp
      pv_in(6+ke1+4) = refatm%delta_t
      pv_in(6+ke1+5) = refatm%h_scal
      inrvert_in = 6+ke1+5
    ENDIF
    IF (refatm%irefatm == 3) THEN
      refatm%bvref = zrbuf (4+ke1+1)
      pv_in(6+ke1+1) = 0.0_wp
      pv_in(6+ke1+2) = 0.0_wp
      pv_in(6+ke1+3) = 0.0_wp
      pv_in(6+ke1+4) = refatm%bvref
      inrvert_in = 6+ke1+4
    ENDIF

    lgot_vcoord   = .TRUE.
    l_ke_in_input = .TRUE.
  ENDIF    ! yformat = 'bina' and yextension = 'o'

#ifdef MESSY
! um_ak_20130925+
ENDIF
#endif

#ifdef I2CINC
ENDIF
! um_ak_20130925-
#endif

#ifdef NETCDF

#ifndef I2CINC
  IF (yformat == 'ncdf') THEN
#else
! um_ak_20130925+
  ! DO NOT OPEN INITIAL OR BOUNDARY FILES !
  ! DATA IS MADE AVAILABLE BY MESSYMMD
  IF ((ydata /= 'initial' .AND. ydata /= 'boundary')  &
       .OR. (.NOT. L_IS_CLIENT))  THEN
#endif

#ifdef MESSY
  ! DO NOT OPEN COSMO RESTART FILES
  ! RESTART IS DONE VIA CHANNEL IN MESSY
  IF ( yextension /='o' .OR. yextension /= 'n') THEN
#endif

    ! the number of dimensions ndims_id_in is now hardcoded in data_io,
    ! because there just are certain different dimension IDs

    ALLOCATE (idims_id_in(1:ndims_id_in))

    ! read the global attributes and definitions
    ! Again bug fix: use ydatchk2 defined above
    CALL read_nc_gdefs(nufile,ie_tot, je_tot, ke_tot, ke_soil, ydatchk2, &
                       startlat_tot, startlon_tot, dlon, dlat,           &
                       icomm_cart, my_cart_id, num_compute, yerrmsg, ierr)
    IF (ierr /= 0) THEN
      CALL model_abort (my_cart_id, 8102, yerrmsg, yroutine)
    ENDIF

    CALL read_nc_vdefs(nufile, listin, nvarin, ivar_id, pollon, pollat, &
                     icomm_cart, my_cart_id, num_compute, yerrmsg, ierr)
    IF (ierr /= 0) THEN
      CALL model_abort (my_cart_id, 8103, yerrmsg, yroutine)
    ENDIF

#ifdef MESSY
ENDIF
#endif

#ifdef I2CINC
ENDIF
#endif

    IF (.NOT. lgot_vcoord) THEN
      ! set the array pv_in (in style l_ke_in_input=.TRUE.) with the
      ! vertical coordinate parameters read in
      IF (refatm%irefatm == 1) THEN
        pv_in( 1) = REAL (vcoord%ivctype, wp)
      ELSE IF (refatm%irefatm == 2) THEN
        pv_in( 1) = REAL (vcoord%ivctype+100, wp)
      ENDIF
      pv_in( 2) = REAL (ke, wp)
      pv_in( 3) = refatm%p0sl
      pv_in( 4) = refatm%t0sl
      pv_in( 5) = refatm%dt0lp
      pv_in( 6) = vcoord%vcflat

      IF     (vcoord%ivctype == 1) THEN
        DO k = 1, ke1
          pv_in(6+k) = vcoord%sigm_coord(k)
        ENDDO
      ELSEIF ( ANY( vcoord%ivctype == (/2,3,4/) ) ) THEN
        DO k = 1, ke1
          pv_in(6+k) = vcoord%vert_coord(k)
        ENDDO
      ENDIF

      inrvert_in = 6+ke1
      IF (refatm%irefatm == 2) THEN
        pv_in(6+ke1+1)    = 0.0_wp
        pv_in(6+ke1+2)    = 0.0_wp
        pv_in(6+ke1+3)    = 0.0_wp
        pv_in(6+ke1+4)    = refatm%delta_t
        pv_in(6+ke1+5)    = refatm%h_scal
        inrvert_in = 6+ke1+5
      ENDIF
      IF (refatm%irefatm == 3) THEN
        pv_in(6+ke1+1)    = 0.0_wp
        pv_in(6+ke1+2)    = 0.0_wp
        pv_in(6+ke1+3)    = 0.0_wp
        pv_in(6+ke1+4)    = refatm%bvref
        inrvert_in = 6+ke1+4
      ENDIF

      lgot_vcoord   = .TRUE.
      l_ke_in_input = .TRUE.

    ENDIF
#ifndef I2CINC
  ENDIF
#endif
#endif

!------------------------------------------------------------------------------
! Section 3: (Endless) loop over all records in the file
!------------------------------------------------------------------------------

#ifdef I2CINC
! um_ak_20130925+
  ! DO NOT OPEN INITIAL OR BOUNDARY FILES !
  ! DATA IS MADE AVAILABLE BY MESSYMMD
  IF ((ydata /= 'initial' .AND. ydata /= 'boundary')  &
       .OR. (.NOT. L_IS_CLIENT))  THEN
#endif

#ifdef MESSY
  ! DO NOT OPEN COSMO RESTART FILES
  ! RESTART IS DONE VIA CHANNEL IN MESSY
  IF ( yextension /='o' .OR. yextension /= 'n') THEN
! um_ak_20130925-
#endif

  IF (izdebug >= 10)                &
    PRINT *, '     Start endless loop to read records: ', ydata, ';  ', yformat

  ! for NetCDF files this is nevertheless done in an ordered manner by going
  ! through all variables in listin. Some organizational variables have to
  ! be set for this:
  izvar_count = 1
  izlev_count = 0

  read_loop: DO

    IF (ntstep == nstart) THEN
      IF (ltime) CALL get_timings (i_initializations, ntstep, dt, izerror)
    ELSE
      IF (ltime) CALL get_timings (i_computations_I, ntstep, dt, izerror)
    ENDIF

  !----------------------------------------------------------------------------
  ! 3.1: Get a record
  !----------------------------------------------------------------------------

    ! some variables have to be initialized
    lrequired = .FALSE.
    igrbednr  = -1         ! to indicate that no record has been read

    ! Every PE gets one record from the file in an ordered manner
    ! (determined by the rank in the communicator icomm_rank). How this is 
    ! done exactly is determined in the routines read_*format*. These routines
    ! have to be called by every PE.
   
    SELECT CASE (yformat)

    CASE ('grb1')

      IF (izdebug >= 20) THEN
        PRINT *, '      Calling read_grib'
      ENDIF

      izmaxlen  = INT (iz_lfd*iwlength, int_ga)
      CALL read_grib    (nufile, izmaxlen, iz_lfd, icomm_cart, iblock,         &
                         iz_rsize_ga, num_compute, lasync_io, ltime_barrier,   &
                         yerrmsg, ierr)

      ! for GRIB1 there is no problem with too big message sizes
      iz_rsize = INT(iz_rsize_ga, iintegers)

      IF (idbg_level >= 20) THEN
        IF (iz_rsize == 0) THEN
          PRINT *, '       EOF reached'
        ELSE
          PRINT *, '       Read a dwdlib record with size ', iz_rsize
        ENDIF
      ENDIF

    CASE ('apix')

      IF (izdebug >= 20) THEN
        PRINT *, '      Calling read_gribapi'
      ENDIF

      izmaxlen  = INT (iz_lfd*iwlength, int_ga)
      CALL read_gribapi  (nufile, izmaxlen, iz_lfd, icomm_cart, iblock,         &
                          iz_rsize_ga, num_compute, lasync_io, ltime_barrier,   &
                          yerrmsg, ierr)

      IF (iz_rsize_ga == 0) THEN
        iz_rsize = 0_iintegers
        IF (idbg_level >= 20) THEN
          PRINT *, '       EOF readched'
        ENDIF
      ELSE
        iz_rsize = 1_iintegers ! just to indicate > 0
        IF (idbg_level >= 20) THEN
          PRINT *, '       Read a gribapi record with size ', iz_rsize_ga, my_cart_id, ierr
        ENDIF
      ENDIF

    CASE ('ncdf')

      CALL read_netcdf  (nufile, ie_tot, je_tot, izvar_count, izlev_count,  &
                         ivar_id, nvarin, idims_id_in, ndims_id_in,         &
                         icomm_cart, my_cart_id, num_compute,               &
                         0_iintegers, 0_iintegers, 0_iintegers, 0_iintegers,&
                         ds_grib, myzvar, myzlev, myzlevtot, iz_rsize,      &
                         imp_grib, lasync_io, yerrmsg, ierr)

    CASE ('bina')

      CALL read_restart (nufile, ie_tot, je_tot, ds_real, ipds_in, npds,    &
                         igds_in, ngds, icomm_cart, my_cart_id, num_compute,&
                         iz_rsize, imp_reals, imp_grib,                     &
                         lasync_io, yerrmsg, ierr)

      IF (izdebug >= 20) THEN
        IF (iz_rsize == 0) THEN
          PRINT *, '       EOF reached'
        ELSE
          PRINT *, '       Got record ', ipds_in(2), ipds_in(7), ipds_in(8)
        ENDIF
      ENDIF

    END SELECT

    IF (ierr /= 0) THEN
       CALL model_abort (my_cart_id, 2015, yerrmsg, yroutine)
    ENDIF

    IF (iz_rsize == 0) THEN
      ! this PE has got no more record because the end of file is reached
      leof = .TRUE.
    ENDIF
 
    IF (ntstep == nstart) THEN
      IF (ltime) CALL get_timings (i_initializations, ntstep, dt, izerror)
    ELSE
      IF (ltime) CALL get_timings (i_read_data, ntstep, dt, izerror)
    ENDIF

  !----------------------------------------------------------------------------
  ! 3.2: If this PE has got a record, de-grib it and put it to zgribarray
  !----------------------------------------------------------------------------

    IF (.NOT. leof) THEN 
      IF     (yformat == 'grb1') THEN
#ifdef GRIBDWD
        ! ds has to have the same precision as the REALs in the grib library
        ds_grib = 0.0_irealgrib

        ! de-grib the data with the routine grbin1 from the grib library
        CALL grbin1(idwdednr, undefgrib, ndims, idims_in, iblock, ibmap,       &
                    ipds_in, igds_in, ibms, ibds, dsup, ds_grib, ierrf)
        IF (ierrf /= 0) THEN
          yerrmsg = 'Error in grbin1'
          CALL model_abort (my_cart_id, 2016, yerrmsg, yroutine)
        ENDIF

        IF (idbg_level >= 20) THEN
          PRINT *, '       Got record ', iz_rsize, ipds_in(2), ipds_in(7), ipds_in(8)
        ENDIF

        ilevtyp  = ipds_in(8)
        igenproc = ipds_in(4)
        igrbednr = 1
        igribid  = 0
#endif

      ELSEIF (yformat == 'apix') THEN

#ifdef GRIBAPI
        ! Build the grib handle

        CALL grib_new_from_message (igribid, iblock, ireturn)
        IF (ireturn /= GRIB_SUCCESS) THEN
          PRINT *,   ' *** Error in grib_api grib_new_from_message  ', ireturn
          yerrmsg  = ' *** Error in grib_api grib_new_from_message'
          CALL model_abort (my_cart_id, 2016, yerrmsg, yroutine)
        ENDIF

        ! get edition number
        CALL grib_get (igribid, 'editionNumber',          igrbednr, ireturn)
        CALL grib_get (igribid, 'shortName',              yzmyname, ireturn)

        ! Get some keys
        CALL grib_get (igribid, 'generatingProcessIdentifier', igenproc  , ireturn)

        ! Get size of data and data itself
        ! Set missing_value before
        CALL grib_set (igribid, 'missingValue', undefgrib)
        CALL grib_get_size (igribid, 'values', igriblen, ireturn)
        IF (igriblen > lds) THEN
          PRINT *, ' *** ERROR: size of message is too big for allocated field: ', igriblen, lds
          yerrmsg = ' *** ERROR: Wrong size of field ***'
          CALL model_abort (my_cart_id, 2037, yerrmsg, yroutine)
        ENDIF

        CALL grib_get (igribid, 'values', ds_grib, ireturn)

        IF (ireturn /= GRIB_SUCCESS) THEN
          PRINT *,   ' *** Error in grib_api grib_get values', ireturn
          yerrmsg  = ' *** Error in grib_api grib_get values'
          CALL model_abort (my_cart_id, 2038, yerrmsg, yroutine)
        ENDIF
#endif
      ELSEIF (yformat == 'bina') THEN
        igrbednr = 1
        igribid  = 0
        igenproc = ipds_in(4)
        ilevtyp  = ipds_in(8)
      ENDIF

      IF (yformat /= 'ncdf') THEN
        ! get the database type from the boundary files. The result files 
        ! will be given the same database type
        IF ( linit .AND. (ydata == 'initial' .OR. ydata == 'restart') ) THEN
          nprocess_ini_in = igenproc
        ELSEIF ( linit .AND. (ydata == 'boundary') ) THEN
          nprocess_bd_in  = igenproc
        ENDIF
        linit = .FALSE.

        ! Determine whether this record is required or not: set lrequired;
        ! typeOfLevel (ytyofle), and the level information (ilevel, ilevtop, ilevbot)
        CALL check_required (listin, nvarin, ydata, yformat, igribid,        &
                             igrbednr, ytyofle, ilevel, ilevtop, ilevbot,    &
                             lrequired, izloc1, izloc2, izloc3, iz_loclist,  &
                             idbg_level, izerror)

        IF (izerror /= 0_iintegers) THEN
          WRITE (yerrmsg,'(A)') 'Errors in check_required'
          CALL model_abort (my_cart_id, izerror, yerrmsg, yroutine)
        ENDIF

        ! and put the data to the structure procarray.
        IF (lrequired) THEN

          ! Check the frame and adapt the b.d. at npstrframe definition, if it is necessary
          IF ( ydata == 'boundary' .AND. lbd_frame ) THEN
            CALL check_frame (ds_grib, izloc1, izloc2, izloc3, ytyofle, ilevtop)
          ENDIF

          IF ( yformat == 'grb1' .OR. yformat == 'apix') THEN

            CALL fill_procarray (izloc1, izloc2, izloc3, igribid, igrbednr,    &
                                 ytyofle, ilevel, ilevtop, ilevbot,            &
                                 iz_rank, iz_dim3, izdebug,                    &
                                 field_grib=ds_grib, proc_gribarray=zgribarray)

          ELSEIF (yformat == 'bina') THEN

            CALL fill_procarray (izloc1, izloc2, izloc3, igribid, igrbednr,    &
                                 ytyofle, ilevel, ilevtop, ilevbot,            &
                                 iz_rank, iz_dim3, izdebug,                    &
                                 field_real=ds_real, proc_realarray=zrealarray)

          ENDIF

          IF  ((iz_rank == -1) .OR. (iz_dim3 == -1)) THEN
            izerror = 5
            WRITE (yerrmsg,'(A,2I3)') 'iz_rank, iz_dim3 could not be determined:  ', iz_rank, iz_dim3
            CALL model_abort (my_cart_id, izerror, yerrmsg, yroutine)
          ENDIF
        ENDIF

      ELSE   ! (yformat == 'ncdf')

        ! set lrequired to .TRUE., because all records are needed
        lrequired = .TRUE.

        ! set izloc1-3 and other organizational variables
        izloc1     = listin(myzvar)%iloc1
        izloc2     = listin(myzvar)%iloc2
        izloc3     = listin(myzvar)%iloc3
        iz_loclist = myzvar
        iz_rank    = var(izloc1,izloc2,izloc3)%rank
        iz_dim3    = myzlev


        DO ip = 0,num_compute-1
          ie_p = isubpos(ip,3) - isubpos(ip,1) + 1 + 2*nboundlines
          je_p = isubpos(ip,4) - isubpos(ip,2) + 1 + 2*nboundlines
          DO j = 1, je_p
            DO i = 1, ie_p
              it = (isubpos(ip,1) - nboundlines) - 1 + i
              jt = (isubpos(ip,2) - nboundlines) - 1 + j
              ij = (jt-1) * ie_tot + it
              zgribarray(i,j,ip) = ds_grib(ij)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

    ENDIF    ! leof

    IF (.NOT. lrequired) THEN
      izloc1 = 1; izloc2 = 2; izloc3 = 1
    ENDIF

  !----------------------------------------------------------------------------
  ! 3.3: Get the vertical coordinate parameters from GRIB
  !----------------------------------------------------------------------------

    IF (lzexch_ednr) THEN

      ! First exchange information on grib edition number from the first PE which
      ! has an atmospheric field to all others
      IF (igrbednr == 1) THEN
        IF     ((yformat == 'grb1') .OR. (yformat == 'bina')) THEN
          IF ((ytyofle == 'hybrid') .OR. (ytyofle == 'hybridLayer')) THEN
            IF (.NOT. lgot_vcoord) THEN
              inrvert_in = igds_in(2)
            ENDIF
            izexch = my_cart_id
          ELSE
            izexch = nproc
          ENDIF
        ELSEIF (yformat == 'apix') THEN
#ifdef GRIBAPI
          CALL grib_get (igribid, 'PVPresent', izpv, ireturn)
          IF ((ytyofle == 'hybrid') .OR. (ytyofle == 'hybridLayer') .AND. (izpv == 1)) THEN
            IF (.NOT. lgot_vcoord) THEN
              CALL grib_get (igribid, 'NV',   inrvert_in, ireturn)
            ENDIF
            izexch = my_cart_id
          ELSE
            izexch = nproc
          ENDIF
#endif
        ENDIF

      ELSEIF (igrbednr == 2) THEN
        ! can only be for 'apix'
        IF ((ytyofle == 'generalVertical') .OR. (ytyofle == 'generalVerticalLayer')) THEN
          izexch = my_cart_id
        ELSE
          izexch = nproc
        ENDIF

      ELSE   ! if no record has been read
        ! this task has got no record to process
        izexch = nproc
      ENDIF

      ! Get the minimum of all values in izexch
      ! this call has to be done by ALL tasks
      IF (num_compute > 1) THEN
        CALL global_values (izexch, 1, 'MIN', imp_integers, icomm_cart,  &
                            -1, yerrmsg, istat)
      ENDIF

      IF (izexch < nproc) THEN
        ! Send grib edition number from task izexch to all others, to check if they
        ! are not using different editions
        IF (num_compute > 1) THEN
          IF (izexch == my_cart_id) THEN
            izexed = igrbednr
          ENDIF
          CALL distribute_values (izexed, 1, izexch, imp_integers, icomm_cart, ierr)
          lzexch_ednr = .FALSE.
        ELSE
          izexed = igrbednr
        ENDIF
      ENDIF

      IF ((idbg_level > 20) .AND. ldebug_io) THEN
        IF (lrequired) THEN
          PRINT *, ' Processor ', my_cart_id, ' has got field ',           &
             var(izloc1,izloc2,izloc3)%name, ytyofle, ' with exchg status: ', izexch
        ELSE
          PRINT *, ' Processor ', my_cart_id, ' has got no required field '
        ENDIF
      ENDIF

      IF ((ydata == 'boundary') .AND. (izexed == 2) .AND.         &
          (uuid_in_string == 'xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx')) THEN
        ! have to get the UUID from the boundary files
        lgot_vcoord = .FALSE.
      ENDIF

    ENDIF ! lzexch_ednr

    ! check, if all multi-level fields have same grib edition number
    IF ((ytyofle == 'generalVertical') .OR. (ytyofle == 'generalVerticalLayer') .OR. &
        (ytyofle == 'hybrid') .OR. (ytyofle == 'hybridLayer')) THEN
      IF ((igrbednr /= izexed) .AND. (igrbednr /= -1)) THEN
        print *, 'whatson:  ', my_cart_id, ytyofle, igrbednr, izexed
        yerrmsg = 'multi-level records with different edition numbers in GRIB file'
        CALL model_abort (my_cart_id, 2004, yerrmsg, yroutine)
      ENDIF
    ENDIF

    ! Now get the vertical coordinate information
    IF (.NOT. lgot_vcoord) THEN

      IF ( (yformat == 'grb1') .OR. (yformat == 'apix') ) THEN
        ! If this is nproc, no processor has got a multi level field
        IF (izexch < nproc) THEN
          IF (izexed == 1) THEN
            ! The processor with lowest id, that has got a multi level field
            ! distributes the vertical coordinate parameters to all others

            ! First distribute number of vertical coordinate parameters
            IF (num_compute > 1) THEN
              CALL distribute_values (inrvert_in, 1, izexch, imp_integers, icomm_cart, ierr)
            ENDIF

            ! Check, whether pv_in is big enough
            IF (inrvert_in > UBOUND(pv_in,1)) THEN
              yerrmsg = 'number of vertical coordinate parameters not consistent with ke:  '
              CALL model_abort (my_cart_id, 2004, yerrmsg, yroutine)
            ENDIF

            IF (izexch == my_cart_id) THEN
              ! get vertical coordinate parameters
              IF     (yformat == 'grb1') THEN
                DO k = 1, inrvert_in
#ifdef GRIBDWD
                  ! for some time the meta data for vertical coordinate parameters
                  ! was not set correct in COSMO and INT2LM
                  ! NOTE: the below can not be replaced by ABS(...) since the smallest
                  !       negative INTEGER may occurr which has no positive equivalent
                  IF ( (igds_in(25+k) > 1000) .OR. (igds_in(25+k) < -1000) ) THEN
                    pv_in(k) = REAL (REFSTF(igds_in(25+k)), wp)
                  ELSE
                    pv_in(k) = REAL (       igds_in(25+k) , wp)
                  ENDIF
#endif
                ENDDO
              ELSEIF (yformat == 'apix') THEN
#ifdef GRIBAPI
                CALL grib_get (igribid, 'pv',    pv_in  , ireturn)
#endif
              ENDIF
            ENDIF

            IF (num_compute > 1) THEN
              CALL distribute_values (pv_in, inrvert_in, izexch, imp_reals, icomm_cart, ierr)
            ENDIF

            CALL get_vertcoord (inrvert_in, pv_in)
            lgot_vcoord = .TRUE.

!US            ! this would only be necessary, if any GRIBOUT block writes GRIB2
!US            IF (my_cart_id == 0) THEN
!US              ! Create a new UUID for the HHL file identifier
!US              CALL uuid_create(uuid_out)
!US
!US              ! Distribute this uuid to all PEs
!US              DO i = 1, 16
!US                charbuf(i:i) = uuid_out(i)
!US              ENDDO
!US            ENDIF
!US
!US            IF (num_compute > 1) THEN
!US              CALL distribute_values (charbuf,  1, 0, imp_character, icomm_cart, ierr)
!US            ENDIF
!US
!US            IF (my_cart_id /= 0) THEN
!US              DO i = 1, 16
!US                uuid_out(i) = charbuf(i:i)
!US              ENDDO
!US            ENDIF
!US
!US            CALL uuid_2char (uuid_out, uuid_out_string)
!US            vcoord%vc_uuid(:) = uuid_out(:)
!US            IF (my_cart_id == 0) THEN
!US              PRINT *, '  Created new UUID for HHL-VGrid:  ', uuid_out_string
!US            ENDIF

          ELSEIF (izexed == 2) THEN
#ifdef GRIBAPI
            ! The processor with lowest id, that has got a multi level field
            ! distributes the meta data to all others

            IF (izexch == my_cart_id) THEN
              CALL grib_get (igribid, 'nlev',                vcoord%nlevels)
              CALL grib_get (igribid, 'numberOfVGridUsed',   vcoord%ivctype)
              CALL grib_get (igribid, 'uuidOfVGrid',         vcoord%vc_uuid)
              intbuf(1) = vcoord%nlevels
              intbuf(2) = vcoord%ivctype
              DO i = 1, 16
                charbuf(i:i) = vcoord%vc_uuid(i)
              ENDDO
            ENDIF

            ! First distribute number of vertical coordinate parameters
            IF (num_compute > 1) THEN
              CALL distribute_values (intbuf,   2, izexch, imp_integers, icomm_cart, ierr)
              CALL distribute_values (charbuf,  1, izexch, imp_character, icomm_cart, ierr)
            ENDIF

            IF (izexch /= my_cart_id) THEN
              vcoord%nlevels = intbuf(1)
              vcoord%ivctype = intbuf(2)
              DO i = 1, 16
                vcoord%vc_uuid(i) = charbuf(i:i)
              ENDDO
            ENDIF
            CALL uuid_2char (vcoord%vc_uuid, uuid_in_string)
            lgot_vcuuid  = .TRUE.
            lgot_vcoord  = .TRUE.    ! but we still have to get the vert_coord parameters!

            ! in the else case, no processor has got a field with leveltyp 150
            ! and the vertical coordinate parameters cannot be determined
#endif
          ENDIF  ! izexed = 1/2/-1
        ENDIF  ! izexch < nproc

      ENDIF ! yformat
    ENDIF   ! lgot_vcoord

  !----------------------------------------------------------------------------
  ! 3.4: Check the record
  !----------------------------------------------------------------------------

    IF ( (yformat == 'grb1') .OR. (yformat == 'apix') .OR. (yformat == 'bina') ) THEN
      IF (ydata == 'boundary') THEN
        CALL get_utc_date (nbdstep+ndiff_ini_bd, ydate_bd, dt, &
                           itype_calendar, ydatchk, ydum, nzday, zhour)
        ydatchk2 = ydatchk
      ELSE
        ! set the actual date to ydatchk2 (which is the initial date here)
        ydatchk2(1:14) = yakdat1(1:14)
      ENDIF

      SELECT CASE (TRIM(var(izloc1,izloc2,izloc3)%name))
      CASE ('HHL')
        lzdate=.FALSE.
      CASE DEFAULT
        lzdate=.TRUE.
      END SELECT

      CALL check_input_grid (igrbednr, igds_in, ngds, ipds_in, npds, igribid,  &
                     yformat, var(izloc1,izloc2,izloc3)%name, ydatchk2,        &
                     ie_tot, je_tot, ke_tot, startlat_tot, startlon_tot,       &
                     dlon, dlat, pollon, pollat, inrvert_in, pv_in, vcoord,    &
                     (.NOT.leof).AND.(lrequired), num_compute,                 &
                     icomm_cart, my_cart_id, lzdate, itype_calendar,           &
                     'COSMO', yerrmsg, ierr)

      IF (ierr /= 0) THEN
        yerrmsg = 'wrong grid description section or Namelist parameters'
        CALL model_abort (my_cart_id, 2004, yerrmsg, yroutine)
      ENDIF
    ENDIF

    ! print the maximum, minimum and meanvalues of each record
    IF ( ((ydata == 'initial') .AND. lchkini)  .OR.        &
         ((ydata == 'boundary') .AND. lchkbd) )  THEN
      ! just to make this call save:
      ieef = INT(izloc2, intgribf)
      CALL check_record (ds_grib, 1, ie_tot, 1, je_tot, 1, 1,              &
               1, ie_tot, 1, je_tot, 1, 1, zundef,                         &
               var(izloc1,izloc2,izloc3)%name, ieef,                       &
               iz_dim3,                                                    &
               (.NOT.leof) .AND. (lrequired), nuchkdat, num_compute,       &
               icomm_cart, my_cart_id, yerrmsg, ierr)
    ELSEIF (ydata == 'restart') THEN
      ds_grib(:) = REAL (ds_real(:), irealgrib)
      ! just to make this call save:
      IF (.NOT. lrequired) THEN
        izloc1 = 1; izloc2 = 1; izloc3 = 1
      ENDIF
      ieef = INT(izloc2, intgribf)
      CALL check_record (ds_grib, 1, ie_tot, 1, je_tot, 1, 1,              &
               1, ie_tot, 1, je_tot, 1, 1, zundef,                         &
               var(izloc1,izloc2,izloc3)%name, ieef,                       &
               iz_dim3,                                                    &
               (.NOT.leof) .AND. (lrequired), nuchkdat, num_compute,       &
               icomm_cart, my_cart_id, yerrmsg, ierr)
    ENDIF

    ! Now check, if it is a HHL record (to get the vertical coordinate of that level)
    IF (TRIM(var(izloc1,izloc2,izloc3)%name) == 'HHL') THEN
#ifdef GRIBAPI
      IF (igrbednr == 2) THEN
        CALL grib_get (igribid, 'scaleFactorOfSecondFixedSurface',  iscalfac,    ireturn)
        IF (iscalfac /= 2) THEN
          yerrmsg = 'wrong scale factor for HHL vertical coordinate'
          CALL model_abort (my_cart_id, 2004, yerrmsg, yroutine)
        ENDIF
        CALL grib_get (igribid, 'scaledValueOfSecondFixedSurface',  iscalval,    ireturn)
        vcoord%vert_coord(ilevtop) = REAL (iscalval, wp) * 0.01_wp
      ENDIF
#endif
    ELSEIF (TRIM(var(izloc1,izloc2,izloc3)%name) == 'HSURF') THEN
      ! or if it is hsurf (to get the minimal value above sea level for a reference profile)
      rr_min = 1.0E10_irealgrib
      DO ij = 1, ie_tot*je_tot
        IF (ds_grib(ij) <= rr_min .AND. ds_grib(ij) >= 0.0_irealgrib) THEN
          ij_min = ij
          rr_min = ds_grib(ij)
        ENDIF
      ENDDO
      ! convert ij to 2-D indices i,j
      igpsl = MOD(ij_min, ie_tot)
      IF (igpsl == 0) THEN
        igpsl = ie_tot
        jgpsl = ij_min / ie_tot
      ELSE
        jgpsl = (ij_min - igpsl) / ie_tot + 1
      ENDIF
      ! print *, 'got a min point:  ', ij_min, rr_min, igpsl, jgpsl, my_cart_id
    ENDIF

    IF (ntstep == nstart) THEN
      IF (ltime) CALL get_timings (i_initializations, ntstep, dt, izerror)
    ELSE
      IF (ltime) CALL get_timings (i_meta_data_r   , ntstep, dt, izerror)
    ENDIF

  !----------------------------------------------------------------------------
  ! 3.5: Distribute record to all PEs and put values to memory
  !----------------------------------------------------------------------------

    distribute_loop: DO iproc = 0, num_compute-1

      ! The following routine handles the distribution of one record to all
      ! others. If the record is not required, the distribute_loop is cycled
      ! (return status istat = -2), if all records are done, the read_loop
      ! is exited (return status istat = -1).
      ! The sender is the processor with rank=iproc, all others receive the
      ! corresponding data into the structure zlocalarray.
      IF (PRESENT(ntlev)) THEN
        nztlev = ntlev
      ELSE
        nztlev = 0
      ENDIF
      IF (yformat == 'bina') THEN
        CALL scatter_data (iproc, izloc1, izloc2, izloc3, iz_loclist,        &
                           iz_rank, iz_dim3, leof, lrequired,                &
                           lzchecklist, nvarin, ydata, yformat, nztlev,      &
                           yextension, istat, realarrays=zrealarray)
      ELSE
        CALL scatter_data (iproc, izloc1, izloc2, izloc3, iz_loclist,        &
                           iz_rank, iz_dim3, leof, lrequired,                &
                           lzchecklist, nvarin, ydata, yformat, nztlev,      &
                           yextension, istat, gribarrays=zgribarray)
      ENDIF

      IF (istat == -1) THEN
#ifdef GRIBAPI
        IF (yformat == 'apix') THEN
          ! Release the grib_api handle
          CALL grib_release(igribid)
        ENDIF
#endif
        EXIT  read_loop          ! all records are done
      ENDIF
      IF (istat == -2) CYCLE distribute_loop    ! record not required
      IF (istat > 0) THEN
        CALL model_abort (my_cart_id, 2017, 'problem in scatter_data', yroutine)
      ENDIF

      ! increase level counter
      iz_countl = iz_countl + 1

    ENDDO distribute_loop

    IF (ntstep == nstart) THEN
      IF (ltime) CALL get_timings (i_initializations, ntstep, dt, izerror)
    ELSE
      IF (ltime) CALL get_timings (i_distribute_data, ntstep, dt, izerror)
    ENDIF

  ENDDO read_loop

!------------------------------------------------------------------------------
! Section 4: Cleanup
!------------------------------------------------------------------------------
 
  ! Save the edition number to the global variable
  IF (ydata == 'initial') THEN
    IF     (yform_read == 'apix') THEN
      ngribednr = izexed
      IF ((idbg_level >= 2) .AND. (my_cart_id == 0)) THEN
        PRINT *, ' GRIB edition number used for atmospheric initial data: ', izexed
      ENDIF
    ELSEIF (yform_read == 'grb1') THEN
      ngribednr = 1
      IF ((idbg_level >= 2) .AND. (my_cart_id == 0)) THEN
        PRINT *, ' GRIB edition number used for atmospheric initial data: ', 1
      ENDIF
    ELSE
      ngribednr = -1
    ENDIF

    ! Just for information about reference atmosphere:
    IF ((idbg_level >= 2) .AND. (my_cart_id == 0)) THEN
      PRINT *, ' Working with reference atmosphere: ', refatm%irefatm
    ENDIF
  ELSEIF (ydata == 'restart') THEN
    ! set ngribednr to -1 for restart data
    ngribednr = -1
  ENDIF

  ! close grib file (has to be called by all PEs)
  CALL close_file (nufile, yformat, icomm_cart, my_cart_id, num_compute, &
                   lasync_io, idbg_level, yerrmsg, ierr)
  IF (ierr /= 0) THEN
     CALL model_abort (my_cart_id, 2017, yerrmsg, yroutine)
  ENDIF

#ifdef MESSY
ENDIF
#endif

#ifdef I2CINC
ENDIF
#endif

!------------------------------------------------------------------------------
! Section 5: Check if all data have been read
!------------------------------------------------------------------------------
 
  IF ( idbg_level > 5 .AND. ldebug_io ) THEN
    PRINT *, ' src_input: check completeness of input data', my_cart_id
  END IF

  IF (ydata == 'initial' .AND. lmulti_layer) THEN
    ! Check, if T_SO(0) or T_S has been read
    !CK This test assumes nzlist_t_so and nzlist_t_s to be positively defined
    !US           and this need not be true in COSMOART
    !CK workaround: Check if they are not default (-1)

    IF (nzlist_t_s >= 0 .OR. nzlist_t_so >= 0 .OR. nzlist_frsn >= 0) THEN
      ! if lmulti_layer, all these variables are included in the initial variable list
      ! and if one variable is set, we assume that all others are also set, because
      ! it really is the meteorological "initial"-list
      IF     ( (.NOT. lzchecklist(0, nzlist_t_so)) .AND.           &
               (      lzchecklist(1, nzlist_t_s )) ) THEN
        ! T_S has been read and is put to T_S0(0)
        t_so(:,:,0,:) = t_s(:,:,:)
        lzchecklist(0, nzlist_t_so) = .TRUE.
      ENDIF

      ! Now we can also set lzchecklist(1,nzlist_t_s) to .TRUE.,
      ! even if it has not been read
      IF (.NOT. lzchecklist(1,nzlist_t_s)) lzchecklist(1,nzlist_t_s) = .TRUE.

      ! Check, if freshsnow has been read; if not, the default setting is used
      IF (.NOT. lzchecklist(1, nzlist_frsn)) THEN
        IF (izdebug > 0) THEN
          PRINT *,' WARNING: ** FRESHSNW could not be read; ',     &
                                'Default setting is used **'
        ENDIF
        lzchecklist(1, nzlist_frsn) = .TRUE.
      ENDIF
    ENDIF
  ENDIF

  IF (ydata == 'initial' .AND. lmulti_snow) THEN

    ! Check, if T_SNOW_M has been read
    IF (nzlist_tsnow_mult >= 0 .OR. nzlist_tsnow >= 0) THEN

      IF     ( (.NOT. lzchecklist(1, nzlist_tsnow_mult)) .AND.   &
               (      lzchecklist(1, nzlist_tsnow )) ) THEN
        ! T_SNOW is put to T_SNOW_M
        DO ksn = 0, ke_snow
          t_snow_mult(:,:,ksn,:) = t_snow(:,:,:)
        END DO
        lzchecklist(:, nzlist_tsnow_mult) = .TRUE.
      ENDIF
    ENDIF

    IF (nzlist_dzh >= 0) THEN
      IF (.NOT. lzchecklist(1, nzlist_dzh)) THEN
        IF (izdebug > 0) THEN
          PRINT *,' WARNING: ** H_SNOW_M could not be read; ',     &
                                'will be initialized in TERRA **'
        ENDIF
        dzh_snow_mult(:,:,:,:) = 0._wp
        lzchecklist(:, nzlist_dzh) = .TRUE.
      ENDIF
    ENDIF

    IF (nzlist_wtot >= 0) THEN
      IF (.NOT. lzchecklist(1, nzlist_wtot)) THEN
        IF (izdebug > 0) THEN
          PRINT *,' WARNING: ** W_SNOW_M could not be read; ',     &
                                'will be initialized in TERRA **'
        ENDIF
        w_snow_mult(:,:,:,:) = 0._wp
        lzchecklist(:, nzlist_wtot) = .TRUE.
      ENDIF
    ENDIF

    IF (nzlist_wliq >= 0) THEN
      IF (.NOT. lzchecklist(1, nzlist_wliq)) THEN
        IF (izdebug > 0) THEN
          PRINT *,' WARNING: ** WLIQ_SNOW could not be read; ',     &
                                'Default setting is used **'
        ENDIF
        wliq_snow(:,:,:,:) = 0._wp
        lzchecklist(:, nzlist_wliq) = .TRUE.
      ENDIF
    ENDIF

    IF (nzlist_rsnow_mult >= 0) THEN
      IF (.NOT. lzchecklist(1, nzlist_rsnow_mult)) THEN
        IF (izdebug > 0) THEN
          PRINT *,' WARNING: ** RHO_SNOW_M could not be read; ',     &
                                'Default setting is used **'
        ENDIF
        rho_snow_mult(:,:,:,:) = 0._wp
        lzchecklist(:, nzlist_rsnow_mult) = .TRUE.
      ENDIF
    ENDIF
  ENDIF  ! linitial and lmulti_snow

  ! check for HHL and P/PP reading
  IF (nzlist_p >= 0 .OR. nzlist_pp >= 0 .OR. nzlist_hhl >= 0) THEN
    IF ((yformat == 'bina') .OR. (yformat == 'ncdf')) THEN
      ! in this case only PP necessary
      IF (nzlist_p   >= 0) lzchecklist(:,nzlist_p  ) = .TRUE.
      IF (nzlist_hhl >= 0) THEN
        IF (ANY(lzchecklist(:,nzlist_hhl) .EQV. .FALSE.)) THEN
          ! we are reading an old restart files without HHL. We assume that we
          ! can reproduce HHL correctly
          lzchecklist(:,nzlist_hhl) = .TRUE.
          lzcomp_hhl = .TRUE.
        ELSE
          ! all levels of HHL have been read and must not be recomputed
          lzcomp_hhl = .FALSE.
        ENDIF
      ENDIF
!     IF (nzlist_hhl >= 0) THEN
!       lzchecklist(:,nzlist_hhl) = .TRUE.
!       lzcomp_hhl = .FALSE.
!     ENDIF
    ELSE
      ! if one of these variables is set, all others are also included in the 
      ! initial variable list, hence the corresponding values are all set, because
      ! it really is the meteorological "initial"-list
      IF     (izexed    == 1) THEN
        ! atmospheric fields are taken from GRIB1, only PP necessary
        IF (nzlist_p   >= 0)    lzchecklist(:,nzlist_p  ) = .TRUE.
        IF (ydata == 'initial' .AND. nzlist_hhl >= 0) lzchecklist(:,nzlist_hhl) = .TRUE.
      ELSEIF (izexed    == 2) THEN
        ! the full pressure p has been stored to pp: this has to be corrected 
        ! after constructing the reference atmosphere
        lzcomp_pp  = .TRUE.

        !  calculate constant fields of the reference atmosphere
        lzcomp_hhl = .FALSE.
        IF (nzlist_pp >= 0) lzchecklist(:,nzlist_pp ) = .TRUE.
      ENDIF
    ENDIF
  ENDIF

#ifdef TWOMOM_SB
  IF (ydata == 'initial' .OR. TRIM(ydata) == 'boundary') THEN
    ! check, if NCxxx variables have been read
    ! if not, do not stop but calculate them afterwards

    DO j2 = 1, nvarin

      SELECT CASE (TRIM(listin(j2)%name))

      CASE ('NCCLOUD')
        IF ( .NOT.ALL(lzchecklist(1:listin(j2)%idimvert,j2)) ) THEN
          lzchecklist(1:listin(j2)%idimvert,j2) = .TRUE.
          lzcalc_qncloud = .TRUE.
        ELSE
          lzcalc_qncloud = .FALSE.
        ENDIF

      CASE ('NCRAIN' )
        IF ( .NOT.ALL(lzchecklist(1:listin(j2)%idimvert,j2)) ) THEN
          lzchecklist(1:listin(j2)%idimvert,j2) = .TRUE.
          lzcalc_qnrain  = .TRUE.
        ELSE
          lzcalc_qnrain  = .FALSE.
        ENDIF

      CASE ('NCSNOW' )
        IF ( .NOT.ALL(lzchecklist(1:listin(j2)%idimvert,j2)) ) THEN
          lzchecklist(1:listin(j2)%idimvert,j2) = .TRUE.
          lzcalc_qnsnow  = .TRUE.
        ELSE
          lzcalc_qnsnow  = .FALSE.
        ENDIF

      CASE ('NCICE'  )
        IF ( .NOT.ALL(lzchecklist(1:listin(j2)%idimvert,j2)) ) THEN
          lzchecklist(1:listin(j2)%idimvert,j2) = .TRUE.
          lzcalc_qnice   = .TRUE.
        ELSE
          lzcalc_qnice   = .FALSE.
        ENDIF

      CASE ('NCGRAUPEL')
        IF ( .NOT.ALL(lzchecklist(1:listin(j2)%idimvert,j2)) ) THEN
          lzchecklist(1:listin(j2)%idimvert,j2) = .TRUE.
          lzcalc_qngraup = .TRUE.
        ELSE
          lzcalc_qngraup = .FALSE.
        ENDIF

      CASE ('NCHAIL' )
        IF ( .NOT.ALL(lzchecklist(1:listin(j2)%idimvert,j2)) ) THEN
          lzchecklist(1:listin(j2)%idimvert,j2) = .TRUE.
          lzcalc_qnhail  = .TRUE.
        ELSE
          lzcalc_qnhail  = .FALSE.
        ENDIF

      END SELECT

    ENDDO
  ENDIF
#endif

  ! check, if all required modellevels are read
  ! CK 20101117 we need an option which does not check if everything has been read.
#ifndef MESSY
  IF (present(lallow_missing)) THEN
#else
  IF (.TRUE.) THEN
#endif
    IF (izdebug > 0) THEN
      PRINT *, '    Omitting check if all ', ydata, ' data has been read.'
    ENDIF
  ELSE
    ! this is the former way
    IF (ALL(lzchecklist(:,:))) THEN
      IF (izdebug > 0) THEN
        PRINT *, '    All variables and levels are read for ', ydata, ' data'
      ENDIF
    ELSE
      IF (my_cart_id == 0) THEN
        PRINT *, '  *** Not all variables / levels could be read for ',  &
                 ydata, ' data   ***'
        PRINT *, '  *** The following levels are missing:    *** '
        DO j2 = 1, nvarin
          IF (listin(j2)%name == 'T_SO      ' .OR. listin(j2)%name == 'T_SNOW_M') THEN
            js = 0
          ELSE
            js = 1
          ENDIF
          DO j1 = js, listin(j2)%idimvert
            IF (lzchecklist(j1,j2) .EQV. .FALSE.) THEN
              PRINT *, '           ', listin(j2)%name, ', level:  ',j1, '   ', &
                        lzchecklist(j1,j2)
            ENDIF
          ENDDO
        ENDDO
      ENDIF
      CALL comm_barrier (icomm_cart, ierr, yerrmsg)
      yerrmsg = 'Not all data available'
      CALL model_abort (my_cart_id, 2004, yerrmsg, yroutine)
    ENDIF
  ENDIF
  
  ! CK 20101117 fill the list to be returned containing the information what has been read
  IF (present(lchecklist)) THEN
    lchecklist(:,:) = lzchecklist(:,:)
  ENDIF

#ifndef MESSY
  IF (yformat == 'bina') THEN
    ! release the unit-number again
    CALL release_unit (nufile)
    IF (nzlist_t_so >= 0 .OR. nzlist_t_s >= 0) THEN
      IF     ( (       lzchecklist(0, nzlist_t_so)) .AND.           &
               ( .NOT. lzchecklist(1, nzlist_t_s )) ) THEN
        ! T_S0(0) has been read and is copied to T_S
         t_s(:,:,:) = t_so(:,:,0,:)
        lzchecklist(1, nzlist_t_s) = .TRUE.
      ENDIF
    ENDIF
  ENDIF
#endif

!------------------------------------------------------------------------------
! Section 6: Vertical grid and reference atmosphere settings
!------------------------------------------------------------------------------

  IF ( idbg_level > 5 .AND. ldebug_io ) THEN
    PRINT *, ' src_input: vertical grid and reference atmosphere', my_cart_id
  END IF

  IF ((ydata == 'initial') .OR. (ydata == 'restart' .AND. yextension == 'o')) THEN

    ! only for grib2 and not for idealized restart runs:
    IF (ngribednr == 2 .AND. (.NOT. lartif_data)) THEN

      IF ( idbg_level > 5 .AND. ldebug_io ) THEN
        PRINT *, ' src_input: gather vert_coord parameters for GRIB2  ', my_cart_id
      END IF

      ! gather all vert_coord parameters from the different tasks, where
      ! they are available
      CALL global_values (vcoord%vert_coord, ke+1, 'MAX', imp_reals, icomm_cart, -1, yerrmsg, ierr)

      ! determine a vcflat out of the given hhl values
      ! (this is only numerically exact for ivctype=2; for ivctype=3,4 it is only an approximation)
      gp_loop: DO j = 1, je
        DO i = 1, ie
          IF (hhl(i,j,ke1) > 10.0_wp) THEN
            ! solve the equation for vcflat: use a higher level, to avoid extinction
            ! in the denominator 
            !  hhl_k = vert_coord_k + (vcflat - vert_coord_k)/vcflat * hsurf
            ! vcflat = vert_coord_k / (1 - (hhl_ijk - vert_coord_k) / hsurf_ij)

            vcflat_k_loop: DO k = 1, ke
              zvck = vcoord%vert_coord(k) / (1.0_wp - (hhl(i,j,k) - vcoord%vert_coord(k)) / hhl(i,j,ke1))
!US           This condition is too sharp: because of the GRIB packing we loose
!US           precision for HHL and there could be bigger differences between
!US           vcoord%vert_coord(k) and hhl(:,:,k) even for flat surfaces
!US           IF (NINT (zvck - vcoord%vert_coord(k)) /= 0) THEN
              IF (ABS(NINT (zvck - vcoord%vert_coord(k))) > 1) THEN
                vcoord%vcflat = zvck
                izexvc = NINT (vcoord%vcflat, iintegers)
                EXIT vcflat_k_loop
                ! otherwise this level is above vcflat
              ENDIF
            ENDDO vcflat_k_loop
!           vcoord%vcflat = (-vcoord%vert_coord(ke1-5) * hhl(i,j,ke1)) /     &
!                        (hhl(i,j,ke1-5) - vcoord%vert_coord(ke1-5) - hhl(i,j,ke1))
!           izexvc = NINT (vcoord%vcflat, iintegers)
            EXIT gp_loop
          ENDIF
        ENDDO
      ENDDO gp_loop

      IF (num_compute > 1) THEN
        CALL global_values (izexvc, 1, 'MIN', imp_integers, icomm_cart,  &
                            -1, yerrmsg, ierr)
      ENDIF
      IF (izexvc == 1000000) THEN
        ! This means that we have a model domain with very flat orography
        !  (hsurf <= 10.0 meters everywhere), or all levels are flat,
        ! so we set vcflat to 10 Meters:
        vcoord%vcflat = 10.0_wp
      ELSE
        vcoord%vcflat = REAL(izexvc, wp)
      END IF
    ENDIF

    ! Additional calculations necessary for ivctype == 3
    IF ( (lzcomp_hhl) .AND. ((vcoord%ivctype == 3) .OR. (vcoord%ivctype == 4)) )THEN
      IF ( idbg_level > 5 .AND. ldebug_io ) THEN
        PRINT *, ' src_input: split topography hsurf  ', my_cart_id
      END IF

      ! split topography hsurf in a large-scale part zhsurfs(:,:,1) and
      ! a small-scale part zhsurfs(:,:,2).

      ! allocate and initialize fields for SLEVE
      ALLOCATE (zhsurf_tot  (ie_tot, je_tot   ),                       &
                zhsurfs_tot (ie_tot, je_tot, 2),    STAT=istat)
      zhsurfs(:,:,:) = 0.0_wp

      ! collect full topo from all PE's
      IF (num_compute == 1 ) THEN     ! we are running on one PE
        zhsurf_tot(:,:) = hsurf(:,:)
      ELSE                            ! we are running on multiple PE's
        CALL gather_field(hsurf     , ie    , je    ,                   &
                          zhsurf_tot, ie_tot, je_tot, 0, ierr)
      ENDIF

      ! split topo on PE 0
      IF (my_cart_id == 0) THEN
        CALL sleve_split_oro(zhsurf_tot, zhsurfs_tot, ie_tot, je_tot, nfltvc, &
                             0_iintegers, svc1, svc2, vcoord%vcflat,          &
                             nuspecif, my_cart_id, ierr, yerrmsg)
      ENDIF

      ! distribute splitted topo to all PE's
      IF (num_compute == 1) THEN      ! we are running on one PE
        zhsurfs(:,:,:) = zhsurfs_tot(:,:,:)
      ELSE                            ! we are running on multiple PE's
        CALL distribute_field (zhsurfs_tot(:,:,1), ie_tot, je_tot,          &
                               zhsurfs    (:,:,1), ie,     je,     0, ierr)
        CALL distribute_field (zhsurfs_tot(:,:,2), ie_tot, je_tot,          &
                               zhsurfs    (:,:,2), ie,     je,     0, ierr)
      ENDIF

      DEALLOCATE (zhsurf_tot, zhsurfs_tot)
    ENDIF

    IF ( idbg_level > 5 .AND. ldebug_io ) THEN
      PRINT *, ' src_input: compute reference atmosphere  ', refatm%irefatm, my_cart_id
    END IF

    SELECT CASE (refatm%irefatm)
    CASE (1)

      ! preliminary choice to be compatible with the current fast_waves solver
      IF ( itype_fast_waves == 2 ) THEN
        lanalyt_calc_t0p0 = .TRUE.   ! necessary only if irefatm=1
      ELSE
        lanalyt_calc_t0p0 = .FALSE.  ! necessary only if irefatm=1
      END IF

      CALL reference_atmosphere                                                &
       ( hhl, p0, p0hl, rho0, t0, t0hl, dp0, hsurf, zhsurfs, ie, je, ke,       &
         refatm, vcoord, svc1, svc2, r_d, g, lanalyt_calc_t0p0,                &
         lzcomp_hhl, yerrmsg, ierr)

    CASE (2)
      CALL reference_atmosphere_2                                              &
       ( hhl, p0, p0hl, rho0, t0, t0hl, dp0, hsurf, zhsurfs, ie, je, ke,       &
         refatm, vcoord, svc1, svc2, r_d, g, lzcomp_hhl, yerrmsg, ierr)
    CASE (3)
      CALL reference_atmosphere_BVconst                                        &
       ( hhl, p0, p0hl, rho0, t0, t0hl, dp0, hsurf, zhsurfs, ie, je, ke,       &
         refatm, vcoord, svc1, svc2, r_d, g, lzcomp_hhl, yerrmsg, ierr)

    CASE default
      yerrmsg = 'ERROR in reference atmosphere specification (wrong value for irefatm)'
      CALL model_abort (my_cart_id, 10007, yerrmsg, yroutine)
    END SELECT

    IF (ierr /= 0) THEN
      CALL model_abort (my_cart_id, 10000+ierr, yerrmsg, 'reference_atmosphere_x')
    ENDIF

    ! For grib2 input respectively lzcomp_hhl=.false., the vcoord%kflat could not 
    ! yet be correctly determined in the reference_atmosphere_xxx routines, because 
    ! it involves a global communication.
    ! We do this now here:

    IF ( .NOT. lzcomp_hhl .AND. (yformat /= 'bina') ) THEN
      ! Determine vcoord%vcflat from given hhl field:

      ! 1) set a missing value for initialization of the search:
      vcoord%kflat = -1

      ! 2) search the level from top to bottom, where for the first time the
      !    minimum and maximum of the local hhl (local PE) are different.
      !    The level above is then the last level with flat coord surfaces,
      !    and this is kflat for the local PE. Ideally, this search should result in the same
      !    kflat on all PEs having variations in their local orography.
      kloop1: DO k = 2, ke+1
        IF (MAXVAL(hhl(:,:,k)) /= MINVAL(hhl(:,:,k))) THEN
          vcoord%kflat = k-1
          EXIT kloop1
        ENDIF
      ENDDO kloop1

      ! 3) In the global grid, the global kflat must be the maximum over
      !    the local kflat's, because for PEs with only sea points (=flat orography), the
      !    above search did not find the correct kflat yet. Therefore, do an MPI exchange
      !    of the information.
      IF (num_compute > 1) THEN
        CALL global_values (vcoord%kflat, 1, 'MAX', imp_integers, icomm_cart,  &
                            -1, yerrmsg, ierr)
      ENDIF

      ! 4) If now the global kflat is still -1, this means that the whole domain
      !    only contains sea points and the coordinate surfaces are really flat everywhere.
      !    In this case, kflat has no real meaning but we set it to ke-1, so that
      !    subsequent computations involving kflat still work:
      IF (vcoord%kflat == -1) THEN
        vcoord%kflat = ke-1
      END IF

    END IF

    ! Exchange a special hhl-profile to all tasks
    ! check, which task has the grid points igpsl, jgpsl  (only one task has these values
    !  set, all others have -1)
    !  (if any:  in case of COSMO-ART it could be that no HSURF is present and no
    !   task has igpsl/jgpsl set)
    IF (num_compute > 1) THEN
      intbuf(1) = igpsl
      intbuf(2) = jgpsl
      CALL global_values (intbuf, 2, 'MAX', imp_integers, icomm_cart, -1, yerrmsg, ierr)
      IF ((igpsl > 0) .AND. (jgpsl > 0)) THEN
        CALL ij_local(igpsl, jgpsl, izlocsl, jzlocsl, itask, ierr)
      ELSE
        itask = -1
      ENDIF
    ELSE
      IF ((igpsl > 0) .AND. (jgpsl > 0)) THEN
        izlocsl = igpsl
        jzlocsl = jgpsl
        itask   = 0
      ELSE
        itask = -1
      ENDIF
    ENDIF
    IF (itask > -1 .AND. (.NOT. lartif_data)) THEN    ! else: nothing can be done
      ! also not idealized restart runs. hhl_prof is defined in src_artifdata.f90 then
      IF (my_cart_id == itask) THEN
        hhl_prof(0)     = vcoord%vcflat
        hhl_prof(1:ke1) = hhl(izlocsl,jzlocsl,1:ke1)
      ENDIF
      IF (num_compute > 1) THEN
        CALL distribute_values (hhl_prof, ke1+1, itask, imp_reals, icomm_cart, istat)
      ENDIF
    ENDIF

    CALL k_index_of_pressure_levels(  refatm%p0sl, vcoord%sigm_coord,   &
            ke, llm, klv950, klv850, klv800, klv700, klv500, klv400, klv300 )

    ! Now set inrvert_in, pv_in, in case of GRIB2 input
    IF (ngribednr == 2) THEN
      IF (refatm%irefatm == 1) THEN
        pv_in( 1) = REAL (vcoord%ivctype, wp)
      ELSE IF (refatm%irefatm == 2) THEN
        pv_in( 1) = REAL (vcoord%ivctype+100, wp)
      ENDIF
      pv_in( 2) = REAL (ke, wp)
      pv_in( 3) = refatm%p0sl
      pv_in( 4) = refatm%t0sl
      pv_in( 5) = refatm%dt0lp
      pv_in( 6) = vcoord%vcflat

      IF     (vcoord%ivctype == 1) THEN
        DO k = 1, ke1
          pv_in(6+k) = vcoord%sigm_coord(k)
        ENDDO
      ELSEIF ( ANY( vcoord%ivctype == (/2,3,4/) ) ) THEN
        DO k = 1, ke1
          pv_in(6+k) = vcoord%vert_coord(k)
        ENDDO
      ENDIF

      inrvert_in = 6+ke1
      IF (refatm%irefatm == 2) THEN
        pv_in(6+ke1+1)    = 0.0_wp
        pv_in(6+ke1+2)    = 0.0_wp
        pv_in(6+ke1+3)    = 0.0_wp
        pv_in(6+ke1+4)    = refatm%delta_t
        pv_in(6+ke1+5)    = refatm%h_scal
        inrvert_in = 6+ke1+5
      ENDIF
      IF (refatm%irefatm == 3) THEN
        pv_in(6+ke1+1)    = 0.0_wp
        pv_in(6+ke1+2)    = 0.0_wp
        pv_in(6+ke1+3)    = 0.0_wp
        pv_in(6+ke1+4)    = refatm%bvref
        inrvert_in = 6+ke1+4
      ENDIF
      l_ke_in_input = .TRUE.

    ENDIF
  ENDIF

  IF (lzcomp_pp) THEN
    ! Compute the pressure deviation, if only full pressure has been read:
    IF     (ydata == 'initial') THEN
      pp(:,:,:,nnow) = pp(:,:,:,nnow)  -  p0(:,:,:)
      pp(:,:,:,nnew) = pp(:,:,:,nnew)  -  p0(:,:,:)
    ELSEIF (ydata == 'boundary') THEN
      pp_bd(:,:,:,ntlev) = pp_bd(:,:,:,ntlev)  -  p0(:,:,:)
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 7: And now additional calculations
!------------------------------------------------------------------------------

  IF ( idbg_level > 5 .AND. ldebug_io ) THEN
    PRINT *, ' src_input: additional calculations: tracers ', my_cart_id
  END IF

 ! Retrieve the required microphysics tracers
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnew, ptr = qv)
  IF (izerror /= 0) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
#ifndef TWOMOM_SB
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nnew, ptr = qc)
  IF (izerror /= 0) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qr, ptr_tlev = nnew, ptr = qr)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qi, ptr_tlev = nnew, ptr = qi)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qs, ptr_tlev = nnew, ptr = qs)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qg, ptr_tlev = nnew, ptr = qg)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
#else
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nnew, ptr = qc, ptr_bd = qc_bd)
  IF (izerror /= 0) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qr, ptr_tlev = nnew, ptr = qr, ptr_bd = qr_bd)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qi, ptr_tlev = nnew, ptr = qi, ptr_bd = qi_bd)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qs, ptr_tlev = nnew, ptr = qs, ptr_bd = qs_bd)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qg, ptr_tlev = nnew, ptr = qg, ptr_bd = qg_bd)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qh, ptr_tlev = nnew, ptr = qh, ptr_bd = qh_bd)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qnc, ptr_tlev = nnew, ptr = qnc, ptr_bd = qnc_bd)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qnr, ptr_tlev = nnew, ptr = qnr, ptr_bd = qnr_bd)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qni, ptr_tlev = nnew, ptr = qni, ptr_bd = qni_bd)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qns, ptr_tlev = nnew, ptr = qns, ptr_bd = qns_bd)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qng, ptr_tlev = nnew, ptr = qng, ptr_bd = qng_bd)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qnh, ptr_tlev = nnew, ptr = qnh, ptr_bd = qnh_bd)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qnc, ptr_tlev = nnow, ptr = qnc_now)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qnr, ptr_tlev = nnow, ptr = qnr_now)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qni, ptr_tlev = nnow, ptr = qni_now)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qns, ptr_tlev = nnow, ptr = qns_now)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qng, ptr_tlev = nnow, ptr = qng_now)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qnh, ptr_tlev = nnow, ptr = qnh_now)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
#endif

  ! Additional calculations for initial data
  IF (ydata == 'initial' .OR.                                &
     (ydata == 'restart') .AND. (yextension == 'o') ) THEN

    IF ( idbg_level > 5 .AND. ldebug_io ) THEN
      PRINT *, ' src_input: additional calculations: initial/restart ', my_cart_id
    END IF

    !  create land/sea mask (corrected by JPS: fr_land must be >= 0.5)
    nlandpoints = COUNT(fr_land(:,:) >= 0.5_wp)

    ! determine nlandpoints_tot (sum over all nlandpoints)
    IF (num_compute > 1) THEN
      nland = COUNT(fr_land(istartpar:iendpar,jstartpar:jendpar) >= 0.5_wp)
      CALL global_values (nland, 1, 'SUM', imp_integers, icomm_cart, -1,   &
                          yerrmsg, istat)
      nlandpoints_tot = nland
    ELSE
      nlandpoints_tot = nlandpoints
    ENDIF

    WHERE (fr_land(:,:) >= 0.5_wp)
       llandmask(:,:) = .TRUE.
    ELSEWHERE
       llandmask(:,:) = .FALSE.
    END WHERE
    
    IF (llake) THEN
      WHERE (fr_lake(:,:) >= 0.5_wp)
        llakemask(:,:) = .TRUE.
      ELSEWHERE
        llakemask(:,:) = .FALSE.
      END WHERE
    END IF

    isoiltyp(:,:) = NINT(soiltyp(:,:))

    ! Check if land/sea mask is consistent with soiltyp.
    ! Note: if these are inconsistent, terra will simply crash. This has
    !       been repeatedly a problem working with different external parameter
    !       datasets
    i = 0
    i = COUNT(llandmask(:,:) .AND. (isoiltyp(:,:) > 8 .OR. isoiltyp(:,:) < 1))
    IF (i > 1) THEN
      DO j = 1, je
        DO i = 1, ie
          IF (fr_land(i,j) < 0.5_wp) THEN
             PRINT *, i,j, '*******', fr_land(i,j), llandmask(i,j), isoiltyp(i,j)
          ELSE
             PRINT *, i,j, '       ', fr_land(i,j), llandmask(i,j), isoiltyp(i,j)
          ENDIF
        ENDDO
      ENDDO
      yerrmsg = 'FR_LAND and SOILTYP in input data are not compatible'
      CALL model_abort (my_cart_id, 2044, yerrmsg, yroutine)
    ENDIF
  ENDIF

  IF (ydata == 'initial') THEN

    ! Compute initial qrs for timelevel nnew:
    qrs(:,:,:) = qr(:,:,:)
    IF (itype_gscp > 1) THEN
      ! snow is considered
      qrs(:,:,:) = qrs(:,:,:) + qs(:,:,:)
    ENDIF
    IF (itype_gscp > 2) THEN
      ! cloud ice is considered
      qrs(:,:,:) = qrs(:,:,:) + qi(:,:,:)
    ENDIF
    IF (itype_gscp > 3) THEN
      ! graupel is considered
      qrs(:,:,:) = qrs(:,:,:) + qg(:,:,:)
    ENDIF
#ifdef TWOMOM_SB
    IF (itype_gscp >= 2000) THEN
      ! also hail is considered
      qrs(:,:,:) = qrs(:,:,:) + qh(:,:,:)
    ENDIF
#endif

    ! For safety: Compute initial rho for timelevel nnew:
    !  (necessary for TWOMOM_SB, if NCXXX are not in IC and BC)
    IF ( idbg_level > 5 .AND. ldebug_io ) THEN
      WRITE (*,*) " src_input: calling calrho (for initial data)"
    END IF
    CALL calrho( t(:,:,:,nnew), pp(:,:,:,nnew), qv(:,:,:), &
         qc(:,:,:), qrs(:,:,:), p0, rho, ie, je, ke, r_d, &
         rvd_m_o )

  ENDIF

#ifdef TWOMOM_SB
  IF (ydata == 'initial' .OR. TRIM(ydata) == 'boundary') THEN
    ! In case of two-moment-microphysics: if NCX could not be read for analysis 
    ! and bc. data, calculate them from QX in a way which is consistent with the 
    ! parameterization assumptions drawn in itype_gscp = 4:

    ! First, calculate total density in cases when it is needed below 
    ! (NOTE: there is no hail in initial data):
    IF (itype_gscp >= 100) THEN
      IF (TRIM(ydata) == 'initial') THEN
        zrho(:,:,:) = rho(:,:,:)
      ELSEIF (TRIM(ydata) == 'boundary') THEN
        IF ( .NOT.(lfirst .AND. lbdana) ) THEN
          ! .. For simplicity, take zrho for boundary data from the "normal" model fields
          !    at the timelevel of its last computation. This is robust if frame-data are read and is
          !    accurate enough for diagnosing NCXXXX below:
          zrho(:,:,:) = rho(:,:,:)
        ENDIF
      ENDIF

      ! .. Lower limit of very small values of rho:
      zrho = MAX(zrho, 1e-10_wp)

    ENDIF

    DO j2 = 1, nvarin

    SELECT CASE (TRIM(listin(j2)%name))

    CASE ('NCCLOUD')

      IF (itype_gscp >= 100) THEN
        IF (lzcalc_qncloud) THEN
          ! .. NCCLOUD could not be read. In case of initial and boundary data, 
          !     set values consistent WITH the assumptions in Seifert/Beheng scheme:
          ! .. Nothing is done for restart data. The model will terminate soon anyways.
          IF (TRIM(ydata) == 'initial') THEN
            IF ((idbg_level > 20) .AND. ldebug_io) THEN
              WRITE (*,'(i4,a,/,9(i4,a,es12.5,a,es12.5,/))') &
                   my_cart_id, " src_input: calculating qnc with zrho(initial)", &
                   my_cart_id, " src_input: min t   ",MINVAL(t),   "  max t    ",MAXVAL(t), &
                   my_cart_id, " src_input: min qv  ",MINVAL(qv),  "  max qv   ",MAXVAL(qv), &
                   my_cart_id, " src_input: min qc  ",MINVAL(qc),  "  max qc   ",MAXVAL(qc), &
                   my_cart_id, " src_input: min qr  ",MINVAL(qr),  "  max qr   ",MAXVAL(qr), &
                   my_cart_id, " src_input: min qi  ",MINVAL(qi),  "  max qi   ",MAXVAL(qi), &
                   my_cart_id, " src_input: min qs  ",MINVAL(qs),  "  max qs   ",MAXVAL(qs), &
                   my_cart_id, " src_input: min qg  ",MINVAL(qg),  "  max qg   ",MAXVAL(qg), &
                   my_cart_id, " src_input: min qh  ",MINVAL(qh),  "  max qh   ",MAXVAL(qh), &
                   my_cart_id, " src_input: min zrho",MINVAL(zrho),"  max zrho ",MAXVAL(zrho)
            ENDIF
            IF ((idbg_level > 20) .AND. ldebug_io) WRITE (*,*) " src_input: calling set_qnc_from_qc_sb"
            qnc(:,:,:) = set_qnc_from_qc_sb(qc(:,:,:)*zrho) / zrho
            IF (.NOT. l2tls) THEN
              qnc_now(:,:,:) = qnc(:,:,:)
              zqncmax = MAXVAL(qnc(:,:,:))
            ELSE
              zqncmax = MAXVAL(qnc(:,:,:))
            END IF
          ELSE
            IF (.NOT.(lfirst .AND. lbdana)) THEN
              IF ((idbg_level > 20) .AND. ldebug_io) WRITE (*,*) " src_input: calculating qnc"
              qnc_bd(:,:,:,ntlev) = set_qnc_from_qc_sb(qc_bd(:,:,:,ntlev)*zrho) / zrho
            END IF
            zqncmax = MAXVAL(qnc_bd(:,:,:,ntlev))
          END IF
          lzchecklist(1:listin(j2)%idimvert,j2) = .TRUE.
          IF (num_compute > 1) THEN
            CALL global_values(zqncmax,1,'MAX',imp_reals,icomm_cart,0,yerrmsg,istat)
          ENDIF
          IF (my_cart_id == 0) THEN
            WRITE (*,*) '===== '//TRIM(listin(j2)%name)// &
                 ' is set with interpreted values for '// &
                 & TRIM(ydata)//' data. MAX = ',zqncmax
          ENDIF
        ENDIF
      ELSE
        WRITE (*,*) '===== '//TRIM(listin(j2)%name)// &
             ' is not needed and allowed for '//TRIM(ydata)//' data!'
        lzchecklist(1:listin(j2)%idimvert,j2) = .FALSE.
      END IF

    CASE ('NCRAIN')

      IF (itype_gscp >= 100) THEN
        IF (lzcalc_qnrain) THEN
          ! NCRAIN could not be read. Set values consistent with the assumptions in Seifert/Beheng scheme:
          IF (TRIM(ydata) == 'initial') THEN
            IF ((idbg_level > 20) .AND. ldebug_io) WRITE (*,*) " src_input: calculating qnr (initial)"
            qnr(:,:,:) = set_qnr_from_qr_sb(qr(:,:,:)*zrho) / zrho
            IF (.NOT. l2tls) THEN
              qnr_now(:,:,:) = qnr(:,:,:)
              zqncmax = MAXVAL(qnr(:,:,:))
            ELSE
              zqncmax = MAXVAL(qnr(:,:,:))
            END IF
          ELSE
            IF (.NOT.(lfirst .AND. lbdana)) THEN
              qnr_bd(:,:,:,ntlev) = set_qnr_from_qr_sb(qr_bd(:,:,:,ntlev)*zrho) / zrho
            END IF
            zqncmax = MAXVAL(qnr_bd(:,:,:,ntlev))
          END IF
          lzchecklist(1:listin(j2)%idimvert,j2) = .TRUE.
          IF (num_compute > 1) THEN
            CALL global_values(zqncmax,1,'MAX',imp_reals,icomm_cart,0,yerrmsg,istat)
          ENDIF
          IF (my_cart_id == 0) THEN
            WRITE (*,*) '===== '//TRIM(listin(j2)%name)// &
                 ' is set with interpreted values for '//&
                 & TRIM(ydata)//' data. MAX = ',zqncmax
          ENDIF
        ENDIF
      ELSE
        WRITE (*,*) '===== '//TRIM(listin(j2)%name)// &
             ' is not needed and allowed for '//TRIM(ydata)//' data!'
        lzchecklist(1:listin(j2)%idimvert,j2) = .FALSE.
      END IF

    CASE ('NCICE')

      IF (itype_gscp >= 100) THEN
        IF (lzcalc_qnice) THEN
          ! NCICE could not be read. Set values consistent with the assumptions in itype_gscp = 4:
          IF (TRIM(ydata) == 'initial') THEN
            IF ((idbg_level > 20) .AND. ldebug_io) WRITE (*,*) " src_input: calculating qni (initial)"
            qni(:,:,:) = set_qni_from_qi_sb(qi(:,:,:)*zrho) / zrho
            IF (.NOT. l2tls) THEN
              qni_now(:,:,:) = qni(:,:,:)
              zqncmax = MAXVAL(qni(:,:,:))
            ELSE
              zqncmax = MAXVAL(qni(:,:,:))
            END IF
          ELSE
            IF (.NOT.(lfirst .AND. lbdana)) THEN
              qni_bd(:,:,:,ntlev) = set_qni_from_qi_sb(qi_bd(:,:,:,ntlev)*zrho) / zrho
            END IF
            zqncmax = MAXVAL(qni_bd(:,:,:,ntlev))
          END IF
          lzchecklist(1:listin(j2)%idimvert,j2) = .TRUE.
          IF (num_compute > 1) THEN
            CALL global_values(zqncmax,1,'MAX',imp_reals,icomm_cart,0,yerrmsg,istat)
          ENDIF
          IF (my_cart_id == 0) THEN
            WRITE (*,*) '===== '//TRIM(listin(j2)%name)// &
                 ' is set with interpreted values for '//&
                 & TRIM(ydata)//' data. MAX = ',zqncmax
          ENDIF
        ENDIF
      ELSE
        WRITE (*,*) '===== '//TRIM(listin(j2)%name)// &
             ' is not needed and allowed for '//TRIM(ydata)//' data!'
        lzchecklist(1:listin(j2)%idimvert,j2) = .FALSE.
      END IF

    CASE ('NCSNOW')

      IF (itype_gscp >= 100) THEN
        IF (lzcalc_qnsnow) THEN
          ! NCSNOW could not be read. Set values consistent with the assumptions in itype_gscp = 4:
          IF (TRIM(ydata) == 'initial') THEN
            IF ((idbg_level > 20) .AND. ldebug_io) WRITE (*,*) " src_input: calculating qns (initial)"
            qns(:,:,:) = set_qns_from_qs_sb(qs(:,:,:)*zrho) / zrho
            IF (.NOT. l2tls) THEN
              qns_now(:,:,:) = qns(:,:,:)
              zqncmax = MAXVAL(qns(:,:,:))
            ELSE
              zqncmax = MAXVAL(qns(:,:,:))
            END IF
          ELSE
            IF (.NOT.(lfirst .AND. lbdana)) THEN
              qns_bd(:,:,:,ntlev) = set_qns_from_qs_sb(qs_bd(:,:,:,ntlev)*zrho) / zrho
            END IF
            zqncmax = MAXVAL(qns_bd(:,:,:,ntlev))
          END IF
          lzchecklist(1:listin(j2)%idimvert,j2) = .TRUE.
          IF (num_compute > 1) THEN
            CALL global_values(zqncmax,1,'MAX',imp_reals,icomm_cart,0,yerrmsg,istat)
          ENDIF
          IF (my_cart_id == 0) THEN
            WRITE (*,*) '===== '//TRIM(listin(j2)%name)//' is set with interpreted values for '//&
                   & TRIM(ydata)//' data. MAX = ',zqncmax
          ENDIF
        ENDIF
      ELSE
        WRITE (*,*) '===== '//TRIM(listin(j2)%name)// &
             ' is not needed and allowed for '//TRIM(ydata)//' data!'
        lzchecklist(1:listin(j2)%idimvert,j2) = .FALSE.
      END IF

    CASE ('NCGRAUPEL')

      IF (itype_gscp >= 100) THEN
        IF (lzcalc_qngraup) THEN
          ! NCGRAUPEL could not be read. Set values consistent with the assumptions in itype_gscp = 4:
          IF (TRIM(ydata) == 'initial') THEN
            IF ((idbg_level > 20) .AND. ldebug_io) WRITE (*,*) " src_input: calculating qng (initial)"
            qng(:,:,:) = set_qng_from_qg_sb(qg(:,:,:)*zrho) / zrho
            IF (.NOT. l2tls) THEN
              qng_now(:,:,:) = qng(:,:,:)
              zqncmax = MAXVAL(qng(:,:,:))
            ELSE
              zqncmax = MAXVAL(qng(:,:,:))
            END IF
          ELSE
            IF (.NOT.(lfirst .AND. lbdana)) THEN
              qng_bd(:,:,:,ntlev) = set_qng_from_qg_sb(qg_bd(:,:,:,ntlev)*zrho) / zrho
            END IF
            zqncmax = MAXVAL(qng_bd(:,:,:,ntlev))
          END IF
          lzchecklist(1:listin(j2)%idimvert,j2) = .TRUE.
          IF (num_compute > 1) THEN
            CALL global_values(zqncmax,1,'MAX',imp_reals,icomm_cart,0,yerrmsg,istat)
          ENDIF
          IF (my_cart_id == 0) THEN
            WRITE (*,*) '===== '//TRIM(listin(j2)%name)// &
                 ' is set with interpreted values for '//&
                 & TRIM(ydata)//' data. MAX = ',zqncmax
          ENDIF
        ENDIF
      ELSE
        WRITE (*,*) '===== '//TRIM(listin(j2)%name)// &
             ' is not needed and allowed for '//TRIM(ydata)//' data!'
        lzchecklist(1:listin(j2)%idimvert,j2) = .FALSE.
      END IF

    CASE ('NCHAIL')

       IF (itype_gscp >= 2000) THEN
        IF (lzcalc_qnhail) THEN
          ! NCHAIL could not be read. Set values to 0.0:
          IF (TRIM(ydata) == 'initial') THEN
            IF ((idbg_level > 20) .AND. ldebug_io) WRITE (*,*) " src_input: setting qnh to zero (initial)"
            qnh(:,:,:) = 0.0_wp
            IF (.NOT. l2tls) THEN
              qnh_now(:,:,:) = qnh(:,:,:)
              zqncmax = MAXVAL(qnh(:,:,:))
            ELSE
              zqncmax = MAXVAL(qnh(:,:,:))
            END IF
          ELSE
            IF (.NOT.(lfirst .AND. lbdana)) THEN
              qnh_bd(:,:,:,ntlev) = 0.0_wp
            END IF
            zqncmax = MAXVAL(qnh_bd(:,:,:,ntlev))
          END IF
          lzchecklist(1:listin(j2)%idimvert,j2) = .TRUE.
          IF (num_compute > 1) THEN
            CALL global_values(zqncmax,1,'MAX',imp_reals,icomm_cart,0,yerrmsg,istat)
          ENDIF
          IF (my_cart_id == 0) THEN
            WRITE (*,*) '===== '//TRIM(listin(j2)%name)// &
                 ' is set with interpreted values for '//&
                    &        TRIM(ydata)//' data. MAX = ',zqncmax
          ENDIF
        ENDIF
      ELSE
        WRITE (*,*) '===== '//TRIM(listin(j2)%name)// &
             ' is not needed and allowed for '//TRIM(ydata)//' data!'
        lzchecklist(1:listin(j2)%idimvert,j2) = .FALSE.
      END IF

    CASE default

      IF (my_cart_id == 0 .AND. .FALSE.) THEN
        WRITE (*,*) '+++++++++ '//TRIM(listin(j2)%name)// &
             ' does not need adjustment ++++++++++++++++++++++++++'
      ENDIF

    END SELECT

    ENDDO
  ENDIF
#endif

#ifdef I2CINC
  IF (.NOT. L_IS_CLIENT) THEN
#endif

  IF (yformat == 'grb1' .OR. yformat == 'apix' .OR. yformat == 'bina') THEN
    ! distribute the database type to all nodes, if iz_countl < num_compute
    IF (iz_countl < num_compute) THEN
      ! not all PEs have nprocess*
      IF (ydata == 'initial') THEN
        intbuf(1) = nprocess_ini_in
        CALL distribute_values (intbuf, 1, 0, imp_integers, icomm_cart, ierr)
        nprocess_ini_in = intbuf(1)
      ELSEIF (ydata == 'boundary') THEN
        intbuf(1) = nprocess_bd_in
        CALL distribute_values (intbuf, 1, 0, imp_integers, icomm_cart, ierr)
        nprocess_bd_in  = intbuf(1)
      ENDIF
    ENDIF
  ENDIF

#ifdef I2CINC
ENDIF
#endif

  IF (ydata == 'initial') THEN

    IF ( idbg_level > 5 .AND. ldebug_io ) THEN
      PRINT *, ' src_input: additional calculations: calps', my_cart_id
    END IF

    !  calculate surface pressure at initial time (timelevel 1 and 2)
    CALL calps ( ps(:,:,nnew), pp(:,:,ke,nnew), t(:,:,ke,nnew),  &
                 qv(:,:,ke  ), qc(:,:,ke     ), qrs(:,:,ke),     &
                 rho0(:,:,ke), p0(:,:,ke), dp0(:,:,ke), ie, je,  &
                 rvd_m_o, r_d, 1, ie, 1, je)
    
    ! calculate temperature tg at initial time (timelevel 1 and 2)
    IF (lmulti_layer) THEN
      ! In case of multi-layer soil model, also t_s has to be set
      ! because it is not read by initial data
      t_s(:,:,nnew) = t_so(:,:,0,nnew)
      t_s(:,:,nnow) = t_so(:,:,0,nnow)

      ! For t_so(:,:,0) all time levels have to be initialized to "save"
      ! the values over sea for the whole LM run (just for output)
      IF (.NOT. l2tls) THEN
        ! 3 time level Leapfrog scheme
        t_so(:,:,0,nold) = t_so(:,:,0,nnew)
      ELSE
        ! 2 time level Runge-Kutta scheme
        t_so(:,:,0,nnow) = t_so(:,:,0,nnew)
      ENDIF
    ENDIF

    IF (lseaice) THEN
      ! For t_ice(:,:,:) all time levels have to be initialized to "save"
      ! the values over land for the whole LM run (just for output)
      IF (.NOT. l2tls) THEN
        ! 3 time level Leapfrog scheme
        t_ice(:,:,nold) = t_ice(:,:,nnew)
      ELSE
        ! 2 time level Runge-Kutta scheme
        t_ice(:,:,nnow) = t_ice(:,:,nnew)
      ENDIF
    ENDIF

    ! Adapt interface parameter t_s for SR tgcom to account for seaice (JPS)
    DO j = 1,je
      DO i = 1,ie
        zt_s(i,j) = t_s(i,j,nnew)
        IF (lseaice) THEN
          IF (.NOT. llandmask(i,j) .AND. h_ice(i,j,nnew) > 0.0_wp) THEN
            zt_s(i,j) = t_ice(i,j,nnew)
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    IF(lmulti_snow) THEN
      CALL tgcom ( t_g(:,:,nnew), t_snow_mult(:,:,1,nnew),         &
                   zt_s(:,:)    , w_snow(:,:,nnew),                &
                   llandmask(:,:) , ie, je, cf_snow,               &
                   1, ie, 1, je)
    ELSE
      CALL tgcom ( t_g(:,:,nnew), t_snow(:,:,nnew),                &
                   zt_s(:,:)    , w_snow(:,:,nnew),                &
                   llandmask(:,:) , ie, je, cf_snow,               &
                   1, ie, 1, je)
    ENDIF

    ! copy the initial data (nnew) into timelevel nnow for leapfrog integration
    IF (.NOT.l2tls) THEN
      ps (:,:,nnow) = ps (:,:,nnew)
      t_g(:,:,nnow) = t_g(:,:,nnew)
    ENDIF

    IF (.NOT. lmulti_layer) THEN
      ! combine or split ground water levels, if necessary
      IF ( (nlgw_lm == 2) .AND. (nlgw_input == 3) ) THEN
        ! combine:  w_g1 = w_g1_input + w_g2_input
        !           w_g2 = w_g3_input
        w_g1(:,:,nnew) = w_g1(:,:,nnew) + w_g2(:,:,nnew)
        w_g2(:,:,nnew) = w_g3(:,:,nnew)
        IF(.NOT.l2tls) THEN
          w_g1(:,:,nnow) = w_g1(:,:,nnew)
          w_g2(:,:,nnow) = w_g2(:,:,nnew)
        ENDIF
      ELSEIF ( (nlgw_lm == 3) .AND. (nlgw_input == 2) ) THEN
        ! split:  w_g1, w_g2 = factor * w_g1_input
        !         w_g3 = wg2_input
        zfac1      = cdzw13 / cdzw12
        zfac2      = cdzw23 / cdzw12
        w_g3(:,:,nnew) = w_g2(:,:,nnew)
        w_g2(:,:,nnew) = w_g1(:,:,nnew) * zfac2
        w_g1(:,:,nnew) = w_g1(:,:,nnew) * zfac1
        IF(.NOT.l2tls) THEN
          w_g1(:,:,nnow) = w_g1(:,:,nnew)
          w_g2(:,:,nnow) = w_g2(:,:,nnew)
          w_g3(:,:,nnow) = w_g3(:,:,nnew)
        ENDIF
      ENDIF
    ENDIF
  ENDIF

  IF (leps) THEN
    IF ( (ydata == 'initial') .AND. (yformat /= 'bina') ) THEN
      ! scale lai, rootdp, plcov with an ensemble-member factor
      lai(:,:) = fac_lai * lai(:,:)
      WHERE (lai(:,:) > rmax_lai) lai(:,:) = rmax_lai
      WHERE (lai(:,:) < rmin_lai) lai(:,:) = rmin_lai

      plcov(:,:) = fac_plcov * plcov(:,:)
      WHERE (plcov(:,:) > rmax_plcov) plcov(:,:) = rmax_plcov
      WHERE (plcov(:,:) < rmin_plcov) plcov(:,:) = rmin_plcov

      rootdp(:,:) = fac_rootdp * rootdp(:,:)
      WHERE (rootdp(:,:) > rmax_rootdp) rootdp(:,:) = rmax_rootdp
      WHERE (rootdp(:,:) < rmin_rootdp) rootdp(:,:) = rmin_rootdp
    ENDIF
  ENDIF

  IF ( (ydata == 'restart') .AND. (yextension == 'o') ) THEN
    ! calculate temperature tg at initial time (timelevel 1 and 2)
    IF (lmulti_layer) THEN
      ! In case of multi-layer soil model, t_s must not be set, 
      ! because it is read by restart data

      ! For t_so(:,:,0) also the 3rd level has to be initialized to "save"
      ! the values over sea for the whole LM run (only for leapfrog scheme)
      IF (.NOT. l2tls) THEN
        t_so(:,:,0,nold) = t_so(:,:,0,nnow)
      ENDIF
    ENDIF
  ENDIF

  ! Additional calculations for boundary data
  IF (ydata == 'boundary') THEN

    IF ( idbg_level > 5 .AND. ldebug_io ) THEN
      PRINT *, ' src_input: additional calculations: boundary  ', my_cart_id
    END IF

    IF (lfirst .AND. nstart == 0) THEN
      ! If the boundary data is from forecast data, the initial data are
      ! mixed with the boundary data and the ground temperature is
      ! calculated again
      CALL id_mix_bd(listin, nvarin)
  
      CALL calps ( ps(:,:,nnew), pp(:,:,ke,nnew), t(:,:,ke,nnew),  &
                   qv(:,:,ke  ), qc(:,:,ke     ), qrs(:,:,ke),     &
                   rho0(:,:,ke), p0(:,:,ke), dp0(:,:,ke), ie, je,  &
                   rvd_m_o, r_d, 1, ie, 1, je)

    ! Adapt interface parameter t_s for SR tgcom to account for seaice (JPS)
      DO j = 1,je
        DO i = 1,ie
          zt_s(i,j) = t_s(i,j,nnew)
          IF (lseaice) THEN
            IF (.NOT. llandmask(i,j) .AND. h_ice(i,j,nnew) > 0.0_wp) THEN
              zt_s(i,j) = t_ice(i,j,nnew)
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      IF(lmulti_snow) THEN
        CALL tgcom ( t_g(:,:,nnew), t_snow_mult(:,:,1,nnew),       &
                     zt_s(:,:)    , w_snow(:,:,nnew),              &
                     llandmask(:,:) , ie, je, cf_snow,             &
                     1, ie, 1, je)
      ELSE
        CALL tgcom ( t_g(:,:,nnew), t_snow(:,:,nnew),              &
                     zt_s(:,:)    , w_snow(:,:,nnew),              &
                     llandmask(:,:) , ie, je, cf_snow,             &
                     1, ie, 1, je)
      ENDIF

      IF(.NOT.l2tls) THEN
        ps (:,:,nnow) = ps (:,:,nnew)
        t_g(:,:,nnow) = t_g(:,:,nnew)
      ENDIF
  
      IF (.NOT. lmulti_layer) THEN
        ! combine or split ground water levels, if necessary
        IF ( (nlgw_lm == 2) .AND. (nlgw_input == 3) ) THEN
          ! combine:  w_g1 = w_g1_input + w_g2_input
          !           w_g2 = w_g3_input
          WHERE (w_g1_bd(:,:,2) /= undef)
            w_g1_bd(:,:,1) = w_g1_bd(:,:,1) + w_g2_bd(:,:,1)
            w_g1_bd(:,:,2) = w_g1_bd(:,:,2) + w_g2_bd(:,:,2)
            w_g2_bd(:,:,1) = w_g3_bd(:,:,1)
            w_g2_bd(:,:,2) = w_g3_bd(:,:,2)
          ENDWHERE
        ELSEIF ( (nlgw_lm == 3) .AND. (nlgw_input == 2) ) THEN
          ! split:  w_g1, w_g2 = factor * w_g1_input
          !         w_g3 = w_g2_input
          zfac1       = cdzw13 / cdzw12
          zfac2       = cdzw23 / cdzw12
          WHERE (w_g1_bd(:,:,2) /= undef)
            w_g3_bd(:,:,1) = w_g2_bd(:,:,1)
            w_g3_bd(:,:,2) = w_g2_bd(:,:,2)
            w_g2_bd(:,:,1) = w_g1_bd(:,:,1) * zfac2
            w_g2_bd(:,:,2) = w_g1_bd(:,:,2) * zfac2
            w_g1_bd(:,:,1) = w_g1_bd(:,:,1) * zfac1
            w_g1_bd(:,:,2) = w_g1_bd(:,:,2) * zfac1
          ENDWHERE
        ENDIF
      ENDIF
    ELSE
      IF (.NOT. lmulti_layer) THEN
        ! combine or split ground water levels, if necessary
        IF ( (nlgw_lm == 2) .AND. (nlgw_input == 3) ) THEN
          ! combine:  w_g1 = w_g1_input + w_g2_input
          !           w_g2 = w_g3_input
          WHERE (w_g1_bd(:,:,2) /= undef)
            w_g1_bd(:,:,ntlev) = w_g1_bd(:,:,ntlev) + w_g2_bd(:,:,ntlev)
            w_g2_bd(:,:,ntlev) = w_g3_bd(:,:,ntlev)
          ENDWHERE
        ELSEIF ( (nlgw_lm == 3) .AND. (nlgw_input == 2) ) THEN
          ! split:  w_g1, w_g2 = factor * w_g1_input
          !         w_g3 = w_g2_input
          zfac1       = cdzw13 / cdzw12
          zfac2       = cdzw23 / cdzw12
          WHERE (w_g1_bd(:,:,2) /= undef)
            w_g3_bd(:,:,ntlev) = w_g2_bd(:,:,ntlev)
            w_g2_bd(:,:,ntlev) = w_g1_bd(:,:,ntlev) * zfac2
            w_g1_bd(:,:,ntlev) = w_g1_bd(:,:,ntlev) * zfac1
          ENDWHERE
        ENDIF
      ENDIF
    ENDIF
  ENDIF

  IF ( idbg_level > 5 .AND. ldebug_io ) THEN
    PRINT *, ' src_input: additional calculations: cleanup   ', my_cart_id
  END IF

  IF (my_cart_id == 0) THEN
    ! Write a blank line to YUCHKDAT
    IF ( lchkini .OR. lchkbd) THEN
      WRITE (nuchkdat,'(A)') '   '
      WRITE (nuchkdat,'(A)') '   '
    ENDIF
    CLOSE (nuchkdat, STATUS='KEEP')
  ENDIF

#ifndef I2CINC
  ! Deallocate arrays for IO
  DEALLOCATE (iblock, ibmap, ds_grib, ds_real, dsup)
  IF (yformat == 'ncdf') THEN
    DEALLOCATE (idims_id_in)
  ENDIF
#else
  IF (ALLOCATED(iblock))      DEALLOCATE(iblock)
  IF (ALLOCATED(ibmap))       DEALLOCATE(ibmap)
  IF (ALLOCATED(ds_grib))     DEALLOCATE(ds_grib)
  IF (ALLOCATED(ds_real))     DEALLOCATE(ds_real)
  IF (ALLOCATED(dsup))        DEALLOCATE(dsup)
  IF (ALLOCATED(idims_id_in)) DEALLOCATE(idims_id_in)
#endif

  IF (ntstep == nstart) THEN
    IF (ltime) CALL get_timings (i_initializations, ntstep, dt, izerror)
  ELSE
    IF (ltime) CALL get_timings (i_computations_I, ntstep, dt, izerror)
  ENDIF

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------
 
END SUBROUTINE organize_input

!==============================================================================
!==============================================================================
!+ Check whether a read (grib) field is required
!------------------------------------------------------------------------------

SUBROUTINE check_required (listin, nlistin, ydata, ylib, igribid, igrbednr,   &
                   ytyofle, ilevel, ilevtop, ilevbot,                         &
                   lrequired, iloc1, iloc2, iloc3, iloclist, idebug, ierrchk)

!------------------------------------------------------------------------------
!
! Description:
!  The following variables are set:
!   - lrequired: logical flag that indicates whether this record is required
!   - iloc1-3:   location of this variable in the LM variable table
!   - ytyofle:   typeOfLevel in grib_api notation
!   - ilevel, ilevbot, ilevtop: level information
!
! Method:
!  For every read record the location in the LM variable table is determined.
!  If no corresponding variable in the list listin is present, the record is 
!  not required. The flag lrequired is set to .FALSE. 
!  If the record is required, iloc1-3 are determined.
!  If inconsistencies are detected, a warning message is printed.
!
!------------------------------------------------------------------------------

! Subroutine input arguments
INTEGER (KIND=iintegers), INTENT(IN)   ::              &
  nlistin         ! number of variables that are read in

TYPE(list_description)  , INTENT(IN)  ::               &
  listin(nlistin) ! List of fields for reading in

CHARACTER (LEN= *)      , INTENT(IN)  ::               &
  ydata,        & ! initial, boundary or restart
  ylib            ! library used (format of the data)
                  !   'grb1','bina': use pds and gds
                  !   'apix'       : use grib_api keys

INTEGER (KIND=intgribf), INTENT(IN)   ::               &
  igribid,      & ! grib handle
  igrbednr,     & ! grib edition number
  idebug          ! if additional debug output shall be printed

! Subroutine output arguments
LOGICAL,                 INTENT(OUT)  ::               &
    lrequired     ! indicates whether this record is required or not

INTEGER (KIND=iintegers), INTENT(OUT) ::               &
    ilevel ,             & ! level
    ilevtop,             & ! top level
    ilevbot,             & ! bot level
    iloc1, iloc2, iloc3, & ! location of this variable in the LM variable table
    iloclist               ! location of this variable in the input list

CHARACTER (LEN=30),       INTENT(OUT) ::               &
  ytyofle                  ! typeOfLevel

INTEGER (KIND=iintegers), INTENT(OUT) ::               &
  ierrchk

! Local scalars:
INTEGER  (KIND=iintegers) :: iloc, n, izret,           &
  itabtyp,    & ! type of the grib table
  ilevtyp,    & ! integer type of level
  iee,        & ! grib number of the variable
  itri,       & ! time range indicator
  iaee,       & ! additional element number
  ityprda,    & ! typeOfProcessedData
  itygepr,    & ! typeOfGeneratingProcess
  iscalval1     ! scaledValueOfFirstFixedSurface

LOGICAL                               ::               &
  lzexternal

CHARACTER (LEN=30)                    ::               &
  yshortname                ! 

!- End of header
!==============================================================================
 
!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  ierrchk   = 0
  izret     = 0
  lrequired = .FALSE.
  iloc1     = 0

!------------------------------------------------------------------------------
! Section 2: Get necessary meta data
!------------------------------------------------------------------------------
  
  IF     (ylib == 'grb1' .OR. ylib == 'bina') THEN

    itabtyp = ipds_in( 2)
    iee     = ipds_in( 7)
    ilevtyp = ipds_in( 8)
    ilevtop = ipds_in( 9)
    ilevbot = ipds_in(10)
    IF     (ilevtyp == 109) THEN
      ilevel = ilevbot
    ELSEIF (ilevtyp == 110) THEN
      ilevel = ilevtop
    ELSE 
      ilevel = 1
    ENDIF
    itri    = ipds_in(19)
    iaee    = ipds_in(41)

    ! set string for leveltyp:
    ytyofle = ylevltypes1(ilevtyp)

    ! set iscalval1 = ilevbot, because this variable is checked for SMA
    iscalval1 = ilevbot

  ELSEIF (ylib == 'apix') THEN
#ifdef GRIBAPI
    IF     (igrbednr == 1) THEN
      CALL grib_get (igribid, 'table2Version',              itabtyp, izret) 
      CALL grib_get (igribid, 'indicatorOfParameter',       iee    , izret) 
      CALL grib_get (igribid, 'timeRangeIndicator',         itri   , izret)
      CALL grib_get (igribid, 'localElementNumber',         iaee   , izret)
      CALL grib_get (igribid, 'indicatorOfTypeOfLevel',     ilevtyp, izret) 

      ! set string for leveltyp:
      ytyofle = ylevltypes1(ilevtyp)

    ELSEIF (igrbednr == 2) THEN
      CALL grib_get (igribid, 'typeOfLevel',                ytyofle, izret) 
      CALL grib_get (igribid, 'localInformationNumber',     iaee   , izret) 
      CALL grib_get (igribid, 'typeOfProcessedData',        ityprda, izret) 
      CALL grib_get (igribid, 'typeOfGeneratingProcess',    itygepr, izret) 
      CALL grib_get (igribid, 'scaledValueOfFirstFixedSurface', iscalval1, izret)
    ENDIF

    ! and the independent keys:
    CALL grib_get (igribid, 'shortName',              yshortname, izret)
    CALL grib_get (igribid, 'level',                      ilevel, izret)
    CALL grib_get (igribid, 'topLevel',                  ilevtop, izret)
    CALL grib_get (igribid, 'bottomLevel',               ilevbot, izret)

    IF (igrbednr == 1) THEN
      ! set scaledValueOfFirstFixedSurface to ilevbot,
      ! because this variable is checked for SMA
      iscalval1 = ilevbot
    ENDIF
#endif
  ENDIF

  ! set error code for output variable
  ierrchk = izret

!------------------------------------------------------------------------------
! Section 3: Get position in variable table and filter out level types for GRIB1
!------------------------------------------------------------------------------
  
  IF     (igrbednr == 1) THEN
    ! determine the short name
    iloc2     = INT(iee, iintegers)
    iloc3     = itabletypes(itabtyp)
    IF (iloc3 == -1) THEN
      ! this table is not known
      IF (idebug > 0) THEN
        PRINT *, 'Warning from check_required:  **table typ unknown**', itabtyp
      ENDIF
      CALL model_abort (my_cart_id, 2004, 'table type unknown', 'check_required')
    ENDIF
  
    ! This search gives the location iloc1 and hence the name of the variable
    ! (name of variable = yshortname already known for grib_api, but we have to 
    !  filter out the level types which are not used, anyhow)

    SELECT CASE(ilevtyp)
    CASE(1,109,110)
      IF (ydata == 'restart') THEN
        ! from some fields, only tri=0 is used
        find_loop_1: DO iloc=1,4
          IF ((var(iloc,iloc2,iloc3)%levtyp == ilevtyp) .AND.     &
              (var(iloc,iloc2,iloc3)%ntri   == itri   ) ) THEN
            iloc1 = iloc
            EXIT find_loop_1
          ENDIF
        ENDDO find_loop_1
      ELSE
        find_loop_2: DO iloc=1,4
          IF  (var(iloc,iloc2,iloc3)%levtyp == ilevtyp) THEN
            iloc1 = iloc
            EXIT find_loop_2
          ENDIF
        ENDDO find_loop_2
      ENDIF
    CASE(211)
      IF ( (itabtyp==  2) .AND. ( (iee== 65) .OR. (iee== 66)) ) THEN
        iloc1 = 2               ! W_SNOW_M        H_SNOW_M
      END IF
      IF ( (itabtyp==201) .AND. ( (iee==133) .OR. (iee==203)) ) THEN
        iloc1 = 2               ! RHO_SNOW_M      T_SNOW_M
      END IF
      IF ( (itabtyp==201) .AND. (iee==137) ) THEN
        iloc1 = 1               ! WLIQ_SNOW
      END IF
    CASE(111,112)
      IF ( (itabtyp==201) .AND. ( (iee==197) .OR. (iee==198) .OR. (iee==199)) ) THEN
        ! these are the multi-dimensional variables for the multi-layer 
        ! model: for these variables, levbot and levtop cannot be compared
        ! But (at the moment) only the first location is used
        iloc1 = 1
      ELSE
        find_loop_3: DO iloc=1,4
          IF( (var(iloc,iloc2,iloc3)%levtyp == ilevtyp) .AND.                &
              (var(iloc,iloc2,iloc3)%levbot == ilevbot) .AND.                &
              (var(iloc,iloc2,iloc3)%levtop == ilevtop)       ) THEN
            IF ((lmulti_layer) .AND. (iee == 85) .AND. (ilevtop == 0)) THEN
              IF (var(iloc,iloc2,iloc3)%name == 'T_S') THEN
                iloc1 = iloc
                EXIT find_loop_3
              ENDIF
            ELSE
              iloc1 = iloc
              EXIT find_loop_3
            ENDIF
          ENDIF
        ENDDO find_loop_3
      ENDIF
    CASE(4,100,102,103)
      ! data on these leveltypes are not used.
      iloc1 = 0
    CASE(2,3,8,105)
      IF (ydata == 'restart') THEN
        ! data on these leveltypes are used for restart files
        find_loop_4: DO iloc=1,4
          IF ((var(iloc,iloc2,iloc3)%levtyp == ilevtyp) .AND.     &
              (var(iloc,iloc2,iloc3)%ntri   == itri   ) ) THEN
            iloc1 = iloc
            EXIT find_loop_4
          ENDIF
        ENDDO find_loop_4
      ELSE
        ! this record is not needed
        iloc1 = 0
      ENDIF
    CASE DEFAULT
      ! data on other than the above leveltypes are unknown in the model.
      iloc1 = 0
      IF (idebug > 0) THEN
        PRINT *, 'Warning from check_required:  **leveltyp unknown**',       &
                                                        ilevtyp, iee, itabtyp
      ENDIF
    END SELECT
  
    ! If iloc1 still is 0, no location was found
    IF (iloc1 == 0) THEN
      yshortname = 'unknown'
    ELSE
      yshortname = var(iloc1,iloc2,iloc3)%name
    ENDIF

  ENDIF    ! igrbednr == 1

!------------------------------------------------------------------------------
! Section 4: Check whether corresponding variable is in input list
!------------------------------------------------------------------------------
  
  DO n = 1, nlistin
    IF (TRIM(listin(n)%name) == TRIM(yshortname)) THEN
      lrequired = .TRUE.
      iloclist = n
      EXIT
    ENDIF
  ENDDO
  
!------------------------------------------------------------------------------
! Section 5: filter out level types which are not used for GRIB2
!------------------------------------------------------------------------------
  
  IF (lrequired) THEN
    IF     (igrbednr == 1) THEN
      ! For GRIB1, we must not read P, but only PP
      IF (TRIM(yshortname) == 'P') THEN
        lrequired = .FALSE.
      ENDIF
    ELSEIF (igrbednr == 2) THEN
      ! Set location in variable table from entry in input list
      iloc1 = listin(iloclist)%iloc1
      iloc2 = listin(iloclist)%iloc2
      iloc3 = listin(iloclist)%iloc3

      IF (yshortname(1:4) == 'W_SO') THEN
        ! mismatch between GRIB1 and GRIB2 leveltypes for soil humidity
        IF (TRIM(ytyofle) /= 'depthBelowLandLayer') THEN
          lrequired = .FALSE.
        ENDIF
      ELSE
        ! filter out unnecessary leveltypes
        IF (yshortname(MAX(1,LEN_TRIM(yshortname)-2):LEN_TRIM(yshortname)) /= '_LK') THEN
          IF (TRIM(ylevltypes2(var(iloc1,iloc2,iloc3)%levtyp)) /= TRIM(ytyofle)) THEN
            lrequired = .FALSE.
          ENDIF
        
          ! for the FLake variables we do not need to check the leveltypes, because 
          ! there are no FLake products with same parameterNumber but different leveltypes
          ! So if the shortname is recognized, it is the correct product
        ENDIF
      ENDIF

      ! For GRIB2, we must not read PP, but only P
      IF (TRIM(yshortname) == 'PP') THEN
        lrequired = .FALSE.
      ENDIF
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 6: Check additional values (like aee and tri) for special input fields
!------------------------------------------------------------------------------
  
  IF (lrequired) THEN
    IF (ydata=='initial') THEN

      IF (idebug > 20) THEN
        PRINT  *, 'Field is required:  ', var(iloc1,iloc2,iloc3)%name, &
                                          yshortname, ytyofle, ilevel, ilevtop, ilevbot
      ENDIF

      ! Check time range indicator (grib1) or typeOfGeneratingProcess (grib2)
      ! of some special fields. If it does not meet the required one the input 
      ! record is discarded

      ! depending on GRIB1/GRIB2 set lzexternal
      IF     (igrbednr == 1) THEN
        IF (itri == 13) THEN
          lzexternal = .FALSE.
        ELSE
          lzexternal = .TRUE.
        ENDIF
      ELSEIF (igrbednr == 2) THEN
        lzexternal = .TRUE.
        IF (ityprda == 0 .AND. itygepr == 202) THEN
          lzexternal = .FALSE.
        ENDIF
      ENDIF

      DO n = 1 , nanfld
        IF (ynaman(n) == var(iloc1,iloc2,iloc3)%name) THEN
  
          IF (lanfld(n)) THEN
            ! external analyses or GME-interpolated fields must be used
            ! (these fields do have tri == 0)
            IF (.NOT. lzexternal) THEN
              lrequired = .FALSE.
              ! not good for GRIB2: ilevbot will be 0 for all layers in the first meter
              ! IF (TRIM(ynaman(n)) == 'T_SO' .AND. ilevbot /= 0) lrequired =.TRUE.
              IF (TRIM(ynaman(n)) == 'T_SO' .AND. iscalval1/= 0) lrequired =.TRUE.
            ELSE
              ! for soil moisture we have to check, if it has to be a field from SMA
              ! for that check the additional element number

              IF ( (TRIM(yshortname) == 'W_G1') .OR.                    &
                   (TRIM(yshortname) == 'W_G2') .OR.                    &
                   (TRIM(yshortname) == 'W_G3') .OR.                    &
                   (TRIM(yshortname) == 'W_SO') ) THEN
                IF (iaee /= nsma_stat) THEN
                  lrequired = .FALSE.
                  IF (idebug > 0) THEN
                    PRINT *, 'Note: analysis field ',var(iloc1,iloc2,iloc3)%name, &
                             ' with add. element number ',iaee,' is discarded'
                  ENDIF
                ELSE
                  IF (idebug > 0) THEN
                    PRINT *, 'Note: analysis field ',var(iloc1,iloc2,iloc3)%name, &
                             ' with add. element number ',iaee,' is used'
                  ENDIF
                ENDIF
              ENDIF

            ENDIF
          ELSE
            ! no external analyses or GME-interpolated fields but fields from
            ! the LM Nudging must be used
            ! (these fields do have tri == 13)
            IF (lzexternal) THEN
              lrequired = .FALSE.
            ENDIF
          ENDIF

          ! Print message, which analysis fields are used / discarded
          IF (idebug > 0) THEN
            IF     (igrbednr == 1) THEN
              IF (lrequired) THEN
                PRINT *,  'Note: analysis field ', ynaman(n) (1:8),           &
                          ' with time range indicator', itri, ' is used'
              ELSE
                PRINT *,  'Note: analysis field ', ynaman(n) (1:8),           &
                          ' with time range indicator', itri, ' is discarded'
              ENDIF
            ELSEIF (igrbednr == 2) THEN
              IF (lrequired) THEN
                PRINT *,  'Note: analysis field ', ynaman(n) (1:8),           &
                          ' with type of generating process ', itygepr, ' is used'
              ELSE
                PRINT *,  'Note: analysis field ', ynaman(n) (1:8),           &
                          ' with type of generating process ', itygepr, ' is discarded'
              ENDIF
            ENDIF
          ENDIF

        ENDIF
      ENDDO

    ELSE

      ! ydata = boundary, restart
      IF (idebug > 20) THEN
        IF     (igrbednr == 1) THEN
          PRINT  *, 'Field is required:  ', yshortname, itabtyp, iee, ilevtyp
        ELSEIF (igrbednr == 2) THEN
          PRINT  *, 'Field is required:  ', yshortname, ytyofle
        ENDIF
      ENDIF

    ENDIF

  ELSE

    IF (idebug > 20) THEN
      IF     (igrbednr == 1) THEN
        PRINT  *, 'Field is NOT required:  ', yshortname, itabtyp, iee, ilevtyp
      ELSEIF (igrbednr == 2) THEN
        PRINT  *, 'Field is NOT required:  ', yshortname, ytyofle
      ENDIF
    ENDIF

  ENDIF

!------------------------------------------------------------------------------
! End of the subroutine 
!------------------------------------------------------------------------------

END SUBROUTINE check_required

!==============================================================================
!==============================================================================
!+ Put data from the record to special structure for distribution
!------------------------------------------------------------------------------

SUBROUTINE fill_procarray                                                  &
             (iloc1, iloc2, iloc3, igribid, iednr, ylevtyp, ilevel,        &
              ilevtop, ilevbot, irank, idim3, idebug,                      &
              field_grib, field_real, proc_gribarray, proc_realarray)

!------------------------------------------------------------------------------
!
! Description:
!  This routine puts the read data record into a special array for 
!  distributing to the other processors. The special array has the form
!  (ie_max, je_max, 0:num_compute-1), where ie_max and je_max are the maximum
!  possible dimensions of a subdomain. The part for every processor is
!  written into an extra level of this array.
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine input arguments
INTEGER (KIND=iintegers), INTENT(IN)   ::               &
    igribid,              & ! grib handler for grib_api records
    iednr,                & ! edition number of grib record to be processed
    iloc1, iloc2, iloc3,  & ! location in the variable table
    idebug                  ! if additional debug output shall be printed

INTEGER (KIND=intgribf) , INTENT(IN)   ::               &
    ilevel,     & ! level 
    ilevtop,    & ! top level
    ilevbot       ! bot level

CHARACTER (LEN=*)       , INTENT(IN)   ::               &
    ylevtyp       ! typeOfLevel

INTEGER (KIND=iintegers),  INTENT(OUT) ::               &
    irank,      & ! rank of this variable
    idim3         ! 

! optional arguments: only the fields for "grib" or for "real" are passed
REAL    (KIND=irealgrib), OPTIONAL, INTENT(IN)   ::               &
    field_grib(ie_tot, je_tot)   ! data of the record

REAL    (KIND=wp),        OPTIONAL, INTENT(IN)   ::               &
    field_real(ie_tot,je_tot)    ! data of the record

! record splitted for subdomains
REAL    (KIND=irealgrib), OPTIONAL, INTENT(OUT)  ::               &
    proc_gribarray(ie_max,je_max,0:num_compute-1)

REAL    (KIND=wp),        OPTIONAL, INTENT(OUT)  ::               &
    proc_realarray(ie_max,je_max,0:num_compute-1)

! Local scalars:
INTEGER  (KIND=iintegers) :: i, ie_p, je_p, n, iscalval1, iscalval2,   &
                             iscalfac1, iscalfac2, ireturn
REAL     (KIND=wp)        :: rdepth1, rdepth2

CHARACTER (LEN=25)        :: yroutine
CHARACTER (LEN=80)        :: yerrmsg

!- End of header
!==============================================================================
 
!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  yroutine='fill_procarray'
  yerrmsg = '           '
  ireturn = 0_iintegers
  irank   = -1
  idim3   = -1

!------------------------------------------------------------------------------
! Section 2: Put the record into the array "procarray"
!------------------------------------------------------------------------------

  ! depending on the rank of the variable the corresponding pointer is used
  SELECT CASE (var(iloc1,iloc2,iloc3)%rank)
  CASE(4)
    irank = 4
!US before:    SELECT CASE (var(iloc1,iloc2,iloc3)%levtyp)
    SELECT CASE (ylevtyp)
    CASE('hybrid')
      !US so wars in Grib1:  idim3 = ilevbot
      idim3 = ilevel
    CASE('hybridLayer')
      idim3 = ilevtop
    CASE('generalVertical')
      idim3 = ilevel
    CASE('generalVerticalLayer')
      idim3 = ilevtop
    CASE('depthBelowLand')
      ! Variables: T_SO:             from GRIB1 and GRIB2
      !            W_SO, W_SO_ICE:   from GRIB2

      IF     (iednr == 1) THEN

        ! In GRIB1, T_SO and W_SO are both coded with leveltype 111 (depthBelowLand);
        ! Normally, W_SO should be coded with leveltype 112 (depthBelowLandLayer),
        ! but the soil layers cannot properly be coded then
        ! The levels are coded in cm (integer values) in 'bottomlevel' (ipds_in(10))
        ! for these variables, ilevbot has to be compared to msoilgrib_in
        DO n=0,ke_soil+1
          IF (msoilgrib(n) == ilevbot) idim3 = n
        ENDDO

      ELSEIF (iednr == 2) THEN
        ! should only be T_SO and using grib_api

#ifdef GRIBAPI
        iscalval1   = -1
        iscalfac1   = -1

        ! in grib2, depth of soil layers are in meter
        ! to get the proper value in cm, scaledValueOfFirstFixedSurface is needed
        CALL grib_get (igribid, 'scaleFactorOfFirstFixedSurface',  iscalfac1,    ireturn)
        CALL grib_get (igribid, 'scaledValueOfFirstFixedSurface',  iscalval1,    ireturn)
        rdepth1 = iscalval1 * rscalefac(iscalfac1)

        IF (ABS(0.0_wp - rdepth1)  < 1.0E-5_wp ) THEN
          idim3 = 0
        ELSE
          ! Determine the level idim3
          DO i=1,ke_soil+1
            IF ( ABS(czmls(i) - rdepth1) < 1.0E-5_wp ) idim3 = i
          ENDDO
        ENDIF
#endif
      ENDIF   ! iednr

    CASE('depthBelowLandLayer')
      ! Variable:  W_SO, W_SO_ICE  for GRIB2 (and using grib_api)
#ifdef GRIBAPI
      iscalval1   = -1
      iscalval2   = -1
      iscalfac1   = -1
      iscalfac2   = -1

      ! in grib2, depth of soil layers are in meter
      ! to get the proper value in cm, scaledValueOfFirstFixedSurface is needed
      CALL grib_get (igribid, 'scaleFactorOfFirstFixedSurface',  iscalfac1,    ireturn)
      CALL grib_get (igribid, 'scaledValueOfFirstFixedSurface',  iscalval1,    ireturn)
      rdepth1 = iscalval1 * rscalefac(iscalfac1)

      ! to get the proper value in cm, scaledValueOfSecondFixedSurface is needed
      CALL grib_get (igribid, 'scaleFactorOfSecondFixedSurface',  iscalfac2,    ireturn)
      CALL grib_get (igribid, 'scaledValueOfSecondFixedSurface',  iscalval2,    ireturn)
      rdepth2 = iscalval2 * rscalefac(iscalfac2)

      ! Determine idim3
      DO i=1,ke_soil+1
        IF ( (ABS(czhls(i-1) - rdepth1) < 1.0E-5_wp ) .AND.  &
             (ABS(czhls(i  ) - rdepth2) < 1.0E-5_wp ) ) idim3 = i
      ENDDO
#endif
    CASE('snowLayer')
      ! these are the multidimensional variables for the multi layer snow model
      DO n=0,ke_snow
        IF (n == ilevtop) idim3 = n
      ENDDO
    CASE DEFAULT
      IF (idebug > 0) THEN
        PRINT *, 'Error in fill_procarray:  **wrong leveltyp**  ', ylevtyp
      ENDIF
    END SELECT
  CASE(3)
    irank = 3
    SELECT CASE (var(iloc1,iloc2,iloc3)%name)
    CASE ('PS        ','T_G       ','QV_S      ','W_SNOW    ',  &
          'T_S       ','T_M       ','W_G1      ','W_G2      ',  &
          'W_G3      ','W_I       ','T_SNOW    ','RHO_SNOW  ',  &
          'H_SNOW    ','H_ICE     ','T_ICE     ','T_MNW_LK  ',  &
          'T_WML_LK  ','T_BOT_LK  ','T_B1_LK   ','C_T_LK    ',  &
          'H_ML_LK   ','H_B1_LK   ')
      ! these are 2D variables with rank=3 because of time dependency
      idim3 = 1
    CASE DEFAULT
!US before:      SELECT CASE (var(iloc1,iloc2,iloc3)%levtyp)
      SELECT CASE (ylevtyp)
      CASE('hybrid')
        !US vorher   idim3 = ilevbot
        idim3 = ilevel
      CASE('hybridLayer')
        idim3 = ilevtop
        ! these are real 3D variables
      CASE('generalVertical')
        idim3 = ilevel
      CASE('generalVerticalLayer')
        idim3 = ilevtop
      END SELECT
    END SELECT
  CASE(2)
    irank = 2
    idim3 = 1
  CASE DEFAULT
    IF (idebug > 0) THEN
      PRINT *, 'Error in fill_procarray: **wrong rank**  ',               &
                var(iloc1,iloc2,iloc3)%rank
    ENDIF
  END SELECT

  IF     (PRESENT(field_grib) .AND. PRESENT(proc_gribarray)) THEN
    DO i=0,num_compute-1
      ie_p = isubpos(i,3) - isubpos(i,1) + 1 + 2*nboundlines
      je_p = isubpos(i,4) - isubpos(i,2) + 1 + 2*nboundlines
      proc_gribarray(1:ie_p,1:je_p,i)       =                             &
          field_grib(isubpos(i,1)-nboundlines : isubpos(i,3)+nboundlines, &
                     isubpos(i,2)-nboundlines : isubpos(i,4)+nboundlines)
    ENDDO
  ELSEIF (PRESENT(field_real) .AND. PRESENT(proc_realarray)) THEN
    DO i=0,num_compute-1
      ie_p = isubpos(i,3) - isubpos(i,1) + 1 + 2*nboundlines
      je_p = isubpos(i,4) - isubpos(i,2) + 1 + 2*nboundlines
      proc_realarray(1:ie_p,1:je_p,i)       =                             &
          field_real(isubpos(i,1)-nboundlines : isubpos(i,3)+nboundlines, &
                     isubpos(i,2)-nboundlines : isubpos(i,4)+nboundlines)
    ENDDO
  ENDIF

!------------------------------------------------------------------------------
! End of the subroutine 
!------------------------------------------------------------------------------

END SUBROUTINE fill_procarray

!==============================================================================
!==============================================================================

SUBROUTINE check_frame (field_grib, iloc1, iloc2, iloc3, ylevtyp, ilevtop)

!------------------------------------------------------------------------------
!
! Description:
!  This routine checks whether fields with a frame are written and sets 
!  the corresponding values to undefined in the variable field.
!
! Method:
!
!------------------------------------------------------------------------------

REAL    (KIND=irealgrib), INTENT(INOUT)::               &
    field_grib(ie_tot, je_tot)   ! data of the record

INTEGER (KIND=iintegers), INTENT(IN)   ::               &
    iloc1, iloc2, iloc3          ! location in the variable table

INTEGER (KIND=intgribf) , INTENT(IN)   ::               &
    ilevtop                      ! top level

CHARACTER (LEN=*)       , INTENT(IN)   ::               &
    ylevtyp                      ! typeOfLevel

!------------------------------------------------------------------------------

! Local Scalars
INTEGER(KIND=iintegers) :: i,j

CHARACTER (LEN=25)        :: yroutine
CHARACTER (LEN=80)        :: yerrmsg

!------------------------------------------------------------------------------

yroutine = 'check_frame'
yerrmsg  = '           '

  IF ( lspubc   .AND. itype_spubc == 1 .AND.                           &
       ( ((ylevtyp == 'hybridLayer') .AND. (ilevtop <= ilevbotnoframe  )) .OR.   &
         ((ylevtyp == 'hybrid'     ) .AND. (ilevtop <= ilevbotnoframe+1)) ) ) THEN
    ! The fields in the Rayleigh damping layer have to be defined on
    ! the full grid. For full level variables these are the fields
    ! above ilevbotnoframe and for half level variables these are the
    ! fields above ilevbotnoframe+1
    DO j = 1, je_tot
      DO i = 1, ie_tot
          IF ( field_grib(i,j) == undefgrib ) THEN
            yerrmsg='found undefined points '
            PRINT*, 'Error in check_frame:    ',yerrmsg
            CALL model_abort (my_cart_id, 2004, yerrmsg, yroutine)
          ENDIF
      ENDDO
    ENDDO
  ELSE
    ! all other fields can be defined with frames
    DO j = 1, je_tot
      DO i = 1, ie_tot
        IF (.NOT.(i > npstrframe .AND. i <= ie_tot-npstrframe .AND. &
                  j > npstrframe .AND. j <= je_tot-npstrframe)) THEN
          IF ( field_grib(i,j) == undefgrib ) THEN
            yerrmsg = 'wrong frame definition or npstrframe '
            CALL model_abort (my_cart_id, 2004, yerrmsg, yroutine)
          ENDIF
        ELSE
          IF (.NOT.(lbdsst .AND. var(iloc1,iloc2,iloc3)%name == 'T_S')) THEN
            field_grib(i,j) = undefgrib
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDIF

END SUBROUTINE check_frame

!==============================================================================
!US  !==============================================================================
!US  !+ Put data from the record to special structure for distribution
!US  !------------------------------------------------------------------------------
!US  
!US  SUBROUTINE fill_gribarray (field, iloc1, iloc2, iloc3, igribid, iednr,      &
!US                             ylevtyp, ilevel, ilevtop, ilevbot, ydata,        &
!US                             procarray, irank, idim3, idebug)
!US  
!US  !------------------------------------------------------------------------------
!US  !
!US  ! Description:
!US  !  This routine puts the read data record into a special array for 
!US  !  distributing to the other processors. The special array has the form
!US  !  (ie_max, je_max, 0:num_compute-1), where ie_max and je_max are the maximum
!US  !  possible dimensions of a subdomain. The part for every processor is
!US  !  written into an extra level of this array.
!US  !
!US  ! Method:
!US  !
!US  !------------------------------------------------------------------------------
!US  
!US  ! Subroutine input arguments
!US  REAL    (KIND=irealgrib), INTENT(INOUT)::               &
!US      field(ie_tot, je_tot)   ! data of the record
!US  
!US  INTEGER (KIND=iintegers), INTENT(IN)   ::               &
!US      igribid,              & ! grib handler for grib_api records
!US      iednr,                & ! edition number of grib record to be processed
!US      iloc1, iloc2, iloc3,  & ! location in the variable table
!US      idebug                  ! if additional debug output shall be printed
!US  
!US  INTEGER (KIND=intgribf) , INTENT(IN)   ::               &
!US      ilevel,     & ! level 
!US      ilevtop,    & ! top level
!US      ilevbot       ! bot level
!US  
!US  CHARACTER (LEN=*)       , INTENT(IN)   ::               &
!US      ylevtyp,    & ! typeOfLevel
!US      ydata         ! determines whether `initial' or `boundary' data
!US  
!US  ! Subroutine output arguments
!US  
!US  REAL    (KIND=irealgrib), INTENT(OUT)  ::               &
!US      procarray(ie_max,je_max,0:num_compute-1) ! record splitted for subdomains
!US  
!US  INTEGER (KIND=iintegers),  INTENT(OUT) ::               &
!US      irank,      & ! rank of this variable
!US      idim3         ! 
!US  
!US  ! Local scalars:
!US  INTEGER  (KIND=iintegers) :: i,j, ie_p, je_p, n, iscalval1, iscalval2,   &
!US                               iscalfac1, iscalfac2, ireturn
!US  REAL     (KIND=wp)        :: rfact, rdepth1, rdepth2
!US  
!US  CHARACTER (LEN=25)        :: yroutine
!US  CHARACTER (LEN=80)        :: yerrmsg
!US  
!US  !- End of header
!US  !==============================================================================
!US   
!US  !------------------------------------------------------------------------------
!US  ! Section 1: Initializations
!US  !------------------------------------------------------------------------------
!US  
!US    yroutine='fill_gribarray'
!US    ireturn = 0_iintegers
!US    irank   = -1
!US    idim3   = -1
!US  
!US  !------------------------------------------------------------------------------
!US  ! Section 2: Put the record into the array "procarray"
!US  !------------------------------------------------------------------------------
!US  
!US    ! depending on the rank of the variable the corresponding pointer is used
!US    SELECT CASE (var(iloc1,iloc2,iloc3)%rank)
!US    CASE(4)
!US      irank = 4
!US  !US before:    SELECT CASE (var(iloc1,iloc2,iloc3)%levtyp)
!US      SELECT CASE (ylevtyp)
!US      CASE('hybrid')
!US        !US so wars in Grib1:  idim3 = ilevbot
!US        idim3 = ilevel
!US      CASE('hybridLayer')
!US        idim3 = ilevtop
!US      CASE('generalVertical')
!US        idim3 = ilevel
!US      CASE('generalVerticalLayer')
!US        idim3 = ilevtop
!US      CASE('depthBelowLand')
!US        ! Variables: T_SO:             from GRIB1 and GRIB2
!US        !            W_SO, W_SO_ICE:   from GRIB2
!US  
!US        IF     (iednr == 1) THEN
!US  
!US          ! In GRIB1, T_SO and W_SO are both coded with leveltype 111 (depthBelowLand);
!US          ! Normally, W_SO should be coded with leveltype 112 (depthBelowLandLayer),
!US          ! but the soil layers cannot properly be coded then
!US          ! The levels are coded in cm (integer values) in 'bottomlevel' (ipds_in(10))
!US          ! for these variables, ilevbot has to be compared to msoilgrib_in
!US          DO n=0,ke_soil+1
!US            IF (msoilgrib(n) == ilevbot) idim3 = n
!US          ENDDO
!US  
!US        ELSEIF (iednr == 2) THEN
!US          ! should only be T_SO and using grib_api
!US  
!US  #ifdef GRIBAPI
!US          iscalval1   = -1
!US          iscalfac1   = -1
!US  
!US          ! in grib2, depth of soil layers are in meter
!US          ! to get the proper value in cm, scaledValueOfFirstFixedSurface is needed
!US          CALL grib_get (igribid, 'scaleFactorOfFirstFixedSurface',  iscalfac1,    ireturn)
!US          CALL grib_get (igribid, 'scaledValueOfFirstFixedSurface',  iscalval1,    ireturn)
!US          rdepth1 = iscalval1 * rscalefac(iscalfac1)
!US  
!US          ! Determine the level idim3
!US          DO i=0,ke_soil+1
!US            IF ( ABS(czmls(i) - rdepth1) < 1.0E-5_wp ) idim3 = i
!US          ENDDO
!US  #endif
!US        ENDIF   ! iednr
!US  
!US      CASE('depthBelowLandLayer')
!US        ! Variable:  W_SO, W_SO_ICE  for GRIB2 (and using grib_api)
!US  #ifdef GRIBAPI
!US        iscalval1   = -1
!US        iscalval2   = -1
!US        iscalfac1   = -1
!US        iscalfac2   = -1
!US  
!US        ! in grib2, depth of soil layers are in meter
!US        ! to get the proper value in cm, scaledValueOfFirstFixedSurface is needed
!US        CALL grib_get (igribid, 'scaleFactorOfFirstFixedSurface',  iscalfac1,    ireturn)
!US        CALL grib_get (igribid, 'scaledValueOfFirstFixedSurface',  iscalval1,    ireturn)
!US        rdepth1 = iscalval1 * rscalefac(iscalfac1)
!US  
!US        ! to get the proper value in cm, scaledValueOfSecondFixedSurface is needed
!US        CALL grib_get (igribid, 'scaleFactorOfSecondFixedSurface',  iscalfac2,    ireturn)
!US        CALL grib_get (igribid, 'scaledValueOfSecondFixedSurface',  iscalval2,    ireturn)
!US        rdepth2 = iscalval2 * rscalefac(iscalfac2)
!US  
!US        ! Determine idim3
!US        DO i=1,ke_soil+1
!US          IF ( (ABS(czhls(i-1) - rdepth1) < 1.0E-5_wp ) .AND.  &
!US               (ABS(czhls(i  ) - rdepth2) < 1.0E-5_wp ) ) idim3 = i
!US        ENDDO
!US  #endif
!US      CASE('snowLayer')
!US        ! these are the multidimensional variables for the multi layer snow model
!US        DO n=0,ke_snow
!US          IF (n == ilevtop) idim3 = n
!US        ENDDO
!US      CASE DEFAULT
!US        IF (idebug > 0) THEN
!US          PRINT *, 'Error in fill_gribarray:  **wrong leveltyp**  ', ylevtyp
!US        ENDIF
!US      END SELECT
!US    CASE(3)
!US      irank = 3
!US      SELECT CASE (var(iloc1,iloc2,iloc3)%name)
!US      CASE ('PS        ','T_G       ','QV_S      ','W_SNOW    ',  &
!US            'T_S       ','T_M       ','W_G1      ','W_G2      ',  &
!US            'W_G3      ','W_I       ','T_SNOW    ','RHO_SNOW  ',  &
!US            'H_SNOW    ','H_ICE     ','T_ICE     ','T_MNW_LK  ',  &
!US            'T_WML_LK  ','T_BOT_LK  ','T_B1_LK   ','C_T_LK    ',  &
!US            'H_ML_LK   ','H_B1_LK   ')
!US        ! these are 2D variables with rank=3 because of time dependency
!US        idim3 = 1
!US      CASE DEFAULT
!US  !US before:      SELECT CASE (var(iloc1,iloc2,iloc3)%levtyp)
!US        SELECT CASE (ylevtyp)
!US        CASE('hybrid')
!US          !US vorher   idim3 = ilevbot
!US          idim3 = ilevel
!US        CASE('hybridLayer')
!US          idim3 = ilevtop
!US          ! these are real 3D variables
!US        END SELECT
!US      END SELECT
!US    CASE(2)
!US      irank = 2
!US      idim3 = 1
!US    CASE DEFAULT
!US      IF (idebug > 0) THEN
!US        PRINT *, 'Error in fill_gribarray: **wrong rank**  ',              &
!US                  var(iloc1,iloc2,iloc3)%rank
!US      ENDIF
!US    END SELECT
!US  
!US    ! Check the frame and adapt the b.d. at npstrframe definition,
!US    ! if it is necessary
!US    IF ( ydata == 'boundary' .AND. lbd_frame ) THEN
!US      IF ( lspubc   .AND. itype_spubc == 1 .AND.                           &
!US           ( ((ylevtyp == 'hybridLayer') .AND. (ilevtop <= ilevbotnoframe  )) .OR.   &
!US             ((ylevtyp == 'hybrid'     ) .AND. (ilevtop <= ilevbotnoframe+1)) ) ) THEN
!US        ! The fields in the Rayleigh damping layer have to be defined on
!US        ! the full grid. For full level variables these are the fields
!US        ! above ilevbotnoframe and for half level variables these are the
!US        ! fields above ilevbotnoframe+1
!US        DO j = 1, je_tot
!US          DO i = 1, ie_tot
!US              IF ( field(i,j) == undefgrib ) THEN
!US                yerrmsg='found undefined points '
!US                PRINT*, 'Error in fill_gribarray: ',yerrmsg
!US                CALL model_abort (my_cart_id, 2004, yerrmsg, yroutine)
!US              ENDIF
!US          ENDDO
!US        ENDDO
!US      ELSE
!US        ! all other fields can be defined with frames
!US        DO j = 1, je_tot
!US          DO i = 1, ie_tot
!US            IF (.NOT.(i > npstrframe .AND. i <= ie_tot-npstrframe .AND. &
!US                      j > npstrframe .AND. j <= je_tot-npstrframe)) THEN
!US              IF ( field(i,j) == undefgrib ) THEN
!US                yerrmsg = 'wrong frame definition or npstrframe '
!US                CALL model_abort (my_cart_id, 2004, yerrmsg, yroutine)
!US              ENDIF
!US            ELSE
!US              IF (.NOT.(lbdsst .AND. var(iloc1,iloc2,iloc3)%name == 'T_S')) THEN
!US                field(i,j) = undefgrib
!US              ENDIF
!US            ENDIF
!US          ENDDO
!US        ENDDO
!US      ENDIF
!US    ENDIF
!US  
!US    DO i=0,num_compute-1
!US      ie_p = isubpos(i,3) - isubpos(i,1) + 1 + 2*nboundlines
!US      je_p = isubpos(i,4) - isubpos(i,2) + 1 + 2*nboundlines
!US      procarray(1:ie_p,1:je_p,i)       =                                     &
!US          field(isubpos(i,1)-nboundlines : isubpos(i,3)+nboundlines,         &
!US                isubpos(i,2)-nboundlines : isubpos(i,4)+nboundlines)
!US    ENDDO
!US  
!US  !------------------------------------------------------------------------------
!US  ! End of the subroutine 
!US  !------------------------------------------------------------------------------
!US  
!US  END SUBROUTINE fill_gribarray
!US  
!US  !==============================================================================
!US  !==============================================================================
!US  !+ Put data from the record to special structure for distribution
!US  !------------------------------------------------------------------------------
!US  
!US  SUBROUTINE fill_realarray (field, iloc1, iloc2, iloc3, ilevtyp, ilevtop,     &
!US                             ilevbot, ydata, procarray, irank, idim3, idebug)
!US  
!US  !------------------------------------------------------------------------------
!US  !
!US  ! Description:
!US  !  This routine puts the read data record into a special array for 
!US  !  distributing to the other processors. The special array has the form
!US  !  (ie_max, je_max, 0:num_compute-1), where ie_max and je_max are the maximum
!US  !  possible dimensions of a subdomain. The part for every processor is
!US  !  written into an extra level of this array.
!US  !
!US  ! Method:
!US  !
!US  !------------------------------------------------------------------------------
!US  
!US  ! Subroutine input arguments
!US  REAL    (KIND=wp),        INTENT(IN)   ::               &
!US      field(ie_tot,je_tot)    ! data of the record
!US  
!US  INTEGER (KIND=iintegers), INTENT(IN)   ::               &
!US      iloc1, iloc2, iloc3,  & ! location in the variable table
!US      idebug                  ! if additional debug output shall be printed
!US  
!US  INTEGER (KIND=intgribf) , INTENT(IN)   ::               &
!US      ilevtyp,    & ! level typ
!US      ilevtop,    & ! top level
!US      ilevbot       ! bot level
!US  
!US  CHARACTER (LEN=*)       , INTENT(IN)   ::               &
!US      ydata          ! determines whether `initial' or `boundary' data
!US  
!US  ! Subroutine output arguments
!US  
!US  REAL    (KIND=wp),        INTENT(OUT)  ::               &
!US      procarray(ie_max,je_max,0:num_compute-1) ! record splitted for subdomains
!US  
!US  INTEGER (KIND=iintegers),  INTENT(OUT) ::               &
!US      irank,      & ! rank of this variable
!US      idim3         ! 
!US  
!US  ! Local scalars:
!US  INTEGER  (KIND=iintegers) :: i,j, ie_p, je_p, n, it, jt, ip, ij
!US  CHARACTER (LEN=25)        :: yroutine
!US  CHARACTER (LEN=80)        :: yerrmsg
!US  
!US  !- End of header
!US  !==============================================================================
!US   
!US  !------------------------------------------------------------------------------
!US  ! Section 1: Initializations
!US  !------------------------------------------------------------------------------
!US  
!US    yroutine='fill_realarray'
!US    irank = -1
!US    idim3 = -1
!US  
!US  !------------------------------------------------------------------------------
!US  ! Section 2: Put the record into the array "procarray"
!US  !------------------------------------------------------------------------------
!US  
!US    ! depending on the rank of the variable the corresponding pointer is used
!US    SELECT CASE (var(iloc1,iloc2,iloc3)%rank)
!US    CASE(4)
!US      irank = 4
!US      SELECT CASE (var(iloc1,iloc2,iloc3)%levtyp)
!US      CASE(109)
!US        idim3 = ilevbot
!US      CASE(110)
!US        idim3 = ilevtop
!US      CASE(111)
!US        ! these are the multidimensional variables for the multi layer
!US        ! soil model
!US        DO n=0,ke_soil+1
!US          IF (msoilgrib(n) == ilevbot) idim3 = n
!US        ENDDO
!US      CASE(211)
!US        ! these are the multidimensional variables for the multi layer
!US        ! snow model
!US        DO n=0,ke_snow
!US          IF (n == ilevtop) idim3 = n
!US        ENDDO
!US      CASE DEFAULT
!US        IF (idebug > 0) THEN
!US          PRINT *, 'Error in fill_realarray:  **wrong leveltyp**  ',     &
!US                   var(iloc1,iloc2,iloc3)%levtyp
!US        ENDIF
!US      END SELECT
!US    CASE(3)
!US      irank = 3
!US      SELECT CASE (var(iloc1,iloc2,iloc3)%name)
!US      CASE ('PS        ','T_G       ','QV_S      ','W_SNOW    ',  &
!US            'T_S       ','T_M       ','W_G1      ','W_G2      ',  &
!US            'W_G3      ','W_I       ','T_SNOW    ','RHO_SNOW  ',  &
!US            'H_SNOW    ','H_ICE     ','T_ICE     ','T_MNW_LK  ',  &
!US            'T_WML_LK  ','T_BOT_LK  ','T_B1_LK   ','C_T_LK    ',  &
!US            'H_ML_LK   ','H_B1_LK   ')
!US        ! these are 2D variables with rank=3 because of time dependency
!US        idim3 = 1
!US      CASE DEFAULT
!US        SELECT CASE (var(iloc1,iloc2,iloc3)%levtyp)
!US        CASE(109)
!US          idim3 = ilevbot
!US        CASE(110)
!US          idim3 = ilevtop
!US          ! these are real 3D variables
!US        END SELECT
!US      END SELECT
!US    CASE(2)
!US      irank = 2
!US      idim3 = 1
!US    CASE DEFAULT
!US      IF (idebug > 0) THEN
!US        PRINT *, 'Error in fill_realarray: **wrong rank**  ',              &
!US                  var(iloc1,iloc2,iloc3)%rank
!US      ENDIF
!US    END SELECT
!US  
!US  !  DO ip=0,num_compute-1
!US  !    ie_p = isubpos(ip,3) - isubpos(ip,1) + 1 + 2*nboundlines
!US  !    je_p = isubpos(ip,4) - isubpos(ip,2) + 1 + 2*nboundlines
!US  !    DO j = 1, je_p
!US  !      DO i = 1, ie_p
!US  !        it = (isubpos(ip,1) - nboundlines) - 1 + i
!US  !        jt = (isubpos(ip,2) - nboundlines) - 1 + j
!US  !        ij = (jt-1) * je_tot + it
!US  !        procarray(i,j,ip) = field(ij)
!US  !      ENDDO
!US  !    ENDDO
!US  
!US    DO i=0,num_compute-1
!US      ie_p = isubpos(i,3) - isubpos(i,1) + 1 + 2*nboundlines
!US      je_p = isubpos(i,4) - isubpos(i,2) + 1 + 2*nboundlines
!US      procarray(1:ie_p,1:je_p,i)       =                                     &
!US          field(isubpos(i,1)-nboundlines : isubpos(i,3)+nboundlines,         &
!US                isubpos(i,2)-nboundlines : isubpos(i,4)+nboundlines)
!US    ENDDO
!US  
!US  !------------------------------------------------------------------------------
!US  ! End of the subroutine 
!US  !------------------------------------------------------------------------------
!US  
!US  END SUBROUTINE fill_realarray
!US  
!US  !==============================================================================
!==============================================================================
!+ get the vertical coordinate parameters from the grid description section
!------------------------------------------------------------------------------

SUBROUTINE get_vertcoord (inv, pv)

!------------------------------------------------------------------------------
!
! Description:
!   Read the vertical coordinate parameters (p0sl, t0sl, dt0lp, vcflat and
!   coordinate parameters) from the grid description section of a record.
!
! Method:
!   The values are de-gribed with the routine REFSTF from the grib library.
!
!------------------------------------------------------------------------------

! Subroutine arguments
INTEGER (KIND=iintegers), INTENT(IN) :: inv
REAL    (KIND=wp),        INTENT(IN) :: pv(inv)

! Local Variables
INTEGER (KIND=iintegers)    :: k, kec, ierrstat, idummy

CHARACTER (LEN=25) yroutine
CHARACTER (LEN=80) yerrmsg
!
!- End of header
!==============================================================================

  yroutine = 'get_vertcoord'

  idummy = NINT (pv( 1), iintegers)

  IF ((idummy >= 1) .AND. (idummy <= 200)) THEN

    ! This is the new grib GDS coding style introduced with LM 3.18

    ! new reference atmosphere has vertical coordinate types 101-103 in Grib
    ! (which correspond to types 1-3 in old reference atmosphere)
    IF (idummy <= 100) THEN
      vcoord%ivctype = idummy
      refatm%irefatm = 1
    ELSE
      ! if the reference atmosphere is determined, vertical coordinate type
      ! is set to values from 1 - 3
      vcoord%ivctype = idummy - 100
      refatm%irefatm = 2
    ENDIF

    ! Check the number of vertical levels:
    IF (NINT(pv(2), iintegers) /= ke) THEN
      WRITE (yerrmsg,'(A)')                                       &
           ' ERROR *** Number of vertical levels of input data is not correct *** '
      PRINT *, ' ERROR *** Number of vertical levels of input data is not correct *** ', &
                NINT(pv(2), iintegers), ke
      ierrstat = 2008
      CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
    ENDIF
    
    refatm%p0sl    = pv( 3)
    refatm%t0sl    = pv( 4)
    refatm%dt0lp   = pv( 5)
    vcoord%vcflat  = pv( 6)
    vcoord%nlevels = ke1

    IF     (vcoord%ivctype == 1) THEN
      DO k = 1,ke1
        vcoord%sigm_coord(k)  = pv(6+k)
      ENDDO
    ELSEIF ( ANY( vcoord%ivctype == (/2,3,4/) ) ) THEN
      DO k = 1,ke1
        vcoord%vert_coord(k)  = pv(6+k)
      ENDDO
    ENDIF

    IF ( ANY( vcoord%ivctype == (/3,4/) ) ) THEN
      ! read three more SLEVE parameters
      svc1   = pv(6 + ke1 + 1)
      svc2   = pv(6 + ke1 + 2)
      nfltvc = NINT (pv(6 + ke1 + 3), iintegers)
    ELSE
      ! just to have an initialization
      svc1   = 0.0_wp
      svc2   = 0.0_wp
      nfltvc = 0_iintegers
    ENDIF

    IF (refatm%irefatm == 2) THEN
      ! read additional parameters for new reference atmosphere
      refatm%delta_t = pv(6 + ke1 + 4)
      refatm%h_scal  = pv(6 + ke1 + 5)
    ENDIF

    l_ke_in_input = .TRUE.
  ELSE

    ! this is the old style of coding the vertical coordinate parameters
    refatm%p0sl   = pv( 1)
    refatm%t0sl   = pv( 2)
    refatm%dt0lp  = pv( 3)
    vcoord%vcflat = pv( 4)

    ! check, how many vertical levels are in input data
    ! check whether pv is descending or ascending
    IF     (pv(5) < pv(6)) THEN
      ! ascending (pressure based): count number of values until 1.0 is reached
      kec = 2
      checkp: DO k = 6, inv
        IF ( (pv(k) > pv(k-1)) .AND. (pv(k) <= 0.9999_wp) ) THEN
          kec = kec + 1
        ELSE
          EXIT checkp
        ENDIF
      ENDDO checkp
    ELSEIF (pv(5) > pv(6)) THEN
      ! descending (height based): count number of values until 0.0 is reached
      kec = 2
      checkh: DO k = 6, inv
        IF ( (pv(k) < pv(k-1)) .AND. (pv(k) >= 0.5_wp) ) THEN
          kec = kec + 1
        ELSE
          EXIT checkh
        ENDIF
      ENDDO checkh
    ENDIF

    IF (kec /= ke1) THEN
      WRITE (yerrmsg,'(A)')                                       &
           ' ERROR *** Number of vertical levels of input data is not correct *** '
      PRINT *, ' ERROR *** Number of vertical levels of input data is not correct *** ', &
                kec, ke
      ierrstat = 2008
      CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
    ELSE
      vcoord%nlevels = ke1
    ENDIF
    
    ! Get the type of the vertical coordinate and the reference atmosphere
    ! As this is not coded properly in the old coding style, we have to use
    ! a rather crude method here: see, how many vertical coordinate parameters
    ! are coded:

    SELECT CASE (inv - ke1 - 4)
    CASE (0, 1)
      ! case=1 can occur because between INT2LM and COSMO we wrote 
      ! ivctype as another meta data (but not between COSMO-Nudging and COSMO)
      refatm%irefatm = 1   ! for irefatm > 1, more coordinate parameters are coded
      IF     (pv(5) < pv(6)) THEN
        ! ascending (pressure based):
        vcoord%ivctype = 1
      ELSEIF (pv(5) > pv(6)) THEN
        ! descending (height based):
        vcoord%ivctype = 2
      ENDIF

    CASE (4)
      refatm%irefatm = 1   ! for irefatm > 1, more coordinate parameters are coded
      ! NOTE: here we allow only ivctype=3, since for SLEVE2 we impose coding using 
      !       l_ke_in_gds=.TRUE. (both in int2lm and cosmo), so this part of the
      !       code should never be called with ivctype=4
      vcoord%ivctype = 3   ! Sleve coordinate  

      ! read three more SLEVE parameters
      svc1   = pv(4 + ke1 + 2)
      svc2   = pv(4 + ke1 + 3)
      nfltvc = NINT (pv(4 + ke1 + 4), iintegers)

    CASE (5)
      refatm%irefatm = 3   ! for irefatm > 1, more coordinate parameters are coded
      vcoord%ivctype = NINT (pv(4 + ke1 + 1), iintegers)

    CASE (6)
      refatm%irefatm = 2   ! for irefatm > 1, more coordinate parameters are coded
      vcoord%ivctype = NINT (pv(4 + ke1 + 1), iintegers)
      IF (vcoord%ivctype > 100) THEN
        vcoord%ivctype = vcoord%ivctype - 100
      ENDIF

    CASE DEFAULT

      yerrmsg = ' ERROR *** ivctype and irefatm for input data not available ***'
      ierrstat = 2008
      CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)

    END SELECT

    IF     (vcoord%ivctype == 1) THEN
      DO k = 1,ke1
        vcoord%sigm_coord(k)  = pv(4+k)
      ENDDO
    ELSEIF ( ANY( vcoord%ivctype == (/2,3,4/) ) ) THEN
      DO k = 1,ke1
        vcoord%vert_coord(k)  = pv(4+k)
      ENDDO
    ENDIF
    IF ( ANY( vcoord%ivctype == (/3,4/) ) ) THEN
      ! read three more SLEVE parameters
      svc1   = pv(4 + ke1 + 2)
      svc2   = pv(4 + ke1 + 3)
      nfltvc = NINT (pv(4 + ke1 + 4), iintegers)
 
      ! Check for meaningful values of svc1, svc2 and nfltvc
      IF ((svc1 > vcoord%vert_coord(1)) .OR. (svc1 < 0.0_wp)) THEN
         yerrmsg = ' ERROR *** svc1 not in allowed range for ivctype = 3/4 ***'
         ierrstat = 2008
         CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
      ENDIF

      IF ((svc2 > vcoord%vert_coord(1)) .OR. (svc2 < 0.0_wp)) THEN
         yerrmsg = ' ERROR *** svc2 not in allowed range for ivctype = 3/4 ***'
         ierrstat = 2008
         CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
      ENDIF

      IF (nfltvc <= 0) THEN
         yerrmsg = ' ERROR *** nfltvc must be greater than or equal to zero ***'
         ierrstat = 2008
         CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
      ENDIF

    ELSE
       ! just to have an initialization
       svc1   = 0.0_wp
       svc2   = 0.0_wp
       nfltvc = 0_iintegers
    ENDIF

    IF (refatm%irefatm == 2) THEN
      ! read additional parameters for new reference atmosphere
      refatm%delta_t = pv(4 + ke1 + 5)
      refatm%h_scal  = pv(4 + ke1 + 6)
    ENDIF

    l_ke_in_input = .FALSE.
  ENDIF

END SUBROUTINE get_vertcoord

!==============================================================================
!==============================================================================
!+ Get first boundary data set from initial data
!------------------------------------------------------------------------------

SUBROUTINE bd_from_id(listin, nvarin)

!------------------------------------------------------------------------------
!
! Description:
!  If analysis data are used in the initial data set, these are copied
!  to the first boundary data set.
!
! Method:
!  For every boundary variable the corresponding initial variable is searched
!  and the values are copied to the boundary variable.
!
!------------------------------------------------------------------------------

! Subroutine arguments
INTEGER (KIND=iintegers), INTENT(IN)  ::                      &
  nvarin          ! number of variables in list
TYPE(list_description)  , INTENT(IN)  ::                      &
  listin(nvarin)  ! List of fields for reading in

! Local variables
INTEGER  (KIND=iintegers) :: ii, i1, i2, i3

!- End of header
!==============================================================================
 
!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 2: Loop over all variables in the input list 
!------------------------------------------------------------------------------

  variable_loop:  DO ii=1,nvarin

    ! Location of this variable in the variable table
    i1 = listin(ii)%iloc1
    i2 = listin(ii)%iloc2
    i3 = listin(ii)%iloc3

    SELECT CASE (var(i1,i2,i3)%rank)
    CASE (4)
      IF (ASSOCIATED(var(i1,i2,i3)%p4) .AND. ASSOCIATED(var(i1,i2,i3)%p4_bd)) THEN
        var(i1,i2,i3)%p4_bd(:,:,:,1) = var(i1,i2,i3)%p4(:,:,:,nnew)
      ENDIF
    CASE (3)
      IF (ASSOCIATED(var(i1,i2,i3)%p3) .AND. ASSOCIATED(var(i1,i2,i3)%p3_bd)) THEN
        var(i1,i2,i3)%p3_bd(:,:,1)   = var(i1,i2,i3)%p3(:,:,nnew)
      ENDIF
    CASE (2)
      IF (lbdclim) THEN
        SELECT CASE (var(i1,i2,i3)%name)
        CASE ('PLCOV     ','LAI       ', 'ROOTDP    ', 'VIO3      ', 'HMO3      ')
          IF (ASSOCIATED(var(i1,i2,i3)%p2) .AND. ASSOCIATED(var(i1,i2,i3)%p3_bd)) THEN
            var(i1,i2,i3)%p3_bd(:,:,1)   = var(i1,i2,i3)%p2(:,:)
          ENDIF
        CASE ('T_CL      ','W_CL      ')
          IF (lmulti_layer) THEN
            IF (ASSOCIATED(var(i1,i2,i3)%p2) .AND. ASSOCIATED(var(i1,i2,i3)%p3_bd)) THEN
              var(i1,i2,i3)%p3_bd(:,:,1)   = var(i1,i2,i3)%p2(:,:)
            ENDIF
          ENDIF
        END SELECT
      ENDIF
    END SELECT

  ENDDO  variable_loop

!------------------------------------------------------------------------------
! End of the subroutine 
!------------------------------------------------------------------------------

END SUBROUTINE bd_from_id

!==============================================================================
!+ Mixing initial and first boundary data set
!------------------------------------------------------------------------------

SUBROUTINE id_mix_bd(listin, nvarin)

!------------------------------------------------------------------------------
!
! Description:
!  If forecast data are used for the boundary data, the initial data set
!  is mixed with the first boundary data set.
!
! Method:
!  For every boundary variable the corresponding initial variable is searched
!  for in the table of the LM variables. The values of these both variables
!  are mixed according to the relaxation scheme.
!
!------------------------------------------------------------------------------
!
! Subroutine arguments
INTEGER (KIND=iintegers), INTENT(IN)  ::                      &
  nvarin          ! number of variables in list
TYPE(list_description)  , INTENT(IN)  ::                      &
  listin(nvarin)  ! List of fields for reading in

! Local variables
INTEGER (KIND=iintegers) :: ii, k, izrmy, i1, i2, i3
INTEGER (KIND=iintegers) :: iz_is, iz_ie, iz_js, iz_je
REAL (KIND=wp)           :: zdtrddt

! Local arrays:
REAL (KIND=wp)        :: zmy(ie,je)
!
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------
 
  zdtrddt = 300.0_wp / dt * MAX(dlon,dlat) / 0.5_wp

!------------------------------------------------------------------------------
! Section 2: Loop over all variables in the input list
!------------------------------------------------------------------------------
 
  variable_loop:  DO ii=1,nvarin

    ! Location of this variable in the variable table
    i1 = listin(ii)%iloc1
    i2 = listin(ii)%iloc2
    i3 = listin(ii)%iloc3

    ! check whether U or V are processed
    IF     (var(i1,i2,i3)%name(1:LEN_TRIM(var(i1,i2,i3)%name)) == 'U') THEN
       izrmy = 2
    ELSEIF (var(i1,i2,i3)%name(1:LEN_TRIM(var(i1,i2,i3)%name)) == 'V') THEN
       izrmy = 3
    ELSE
       izrmy = 1
    ENDIF

    ! define start- and end-indices for computation
    IF(my_cart_neigh(1) == -1) THEN
      iz_is = istart
    ELSE
      iz_is = 1
    ENDIF

    IF(my_cart_neigh(2) == -1) THEN
      IF (var(i1,i2,i3)%name(1:LEN_TRIM(var(i1,i2,i3)%name)) == 'V') THEN
        iz_je = jendv
      ELSE
        iz_je = jend
      ENDIF
    ELSE
      iz_je = je
    ENDIF

    IF(my_cart_neigh(3) == -1) THEN
      IF (var(i1,i2,i3)%name(1:LEN_TRIM(var(i1,i2,i3)%name)) == 'U') THEN
        iz_ie = iendu
      ELSE
        iz_ie = iend
      ENDIF
    ELSE
      iz_ie = ie
    ENDIF
      
    IF(my_cart_neigh(4) == -1) THEN
      iz_js = jstart
    ELSE
      iz_js = 1
    ENDIF

    ! compute zmy
    zmy(:,:) = 1.0_wp
    zmy(iz_is:iz_ie,iz_js:iz_je) = rmy(iz_is:iz_ie,iz_js:iz_je,izrmy)

  !----------------------------------------------------------------------------
  ! 2.3: Mix the fields
  !----------------------------------------------------------------------------
 
    SELECT CASE (var(i1,i2,i3)%rank)
    CASE (4)
      IF (ASSOCIATED(var(i1,i2,i3)%p4) .AND. ASSOCIATED(var(i1,i2,i3)%p4_bd)) THEN
        DO k= LBOUND(var(i1,i2,i3)%p4,3), UBOUND(var(i1,i2,i3)%p4,3)
          WHERE (var(i1,i2,i3)%p4_bd(1:ie,1:je,k,1) /= undef)
            var(i1,i2,i3)%p4(1:ie,1:je,k,nnew) =                            &
                     zmy(1:ie,1:je) * var(i1,i2,i3)%p4_bd(1:ie,1:je,k,1) +  &
           (1.0_wp - zmy(1:ie,1:je))* var(i1,i2,i3)%p4(1:ie,1:je,k,nnew)
          ENDWHERE
          IF (.NOT.l2tls) THEN
             var(i1,i2,i3)%p4(:,:,k,nnow) = var(i1,i2,i3)%p4(:,:,k,nnew)
          ENDIF
        ENDDO
      ENDIF
    CASE (3)
      IF (ASSOCIATED(var(i1,i2,i3)%p3) .AND. ASSOCIATED(var(i1,i2,i3)%p3_bd)) THEN
        WHERE (var(i1,i2,i3)%p3_bd(1:ie,1:je,1) /= undef)
          var(i1,i2,i3)%p3(1:ie,1:je,nnew) =                                &
                   zmy(1:ie,1:je) * var(i1,i2,i3)%p3_bd(1:ie,1:je,1) +      &
         (1.0_wp - zmy(1:ie,1:je))* var(i1,i2,i3)%p3(1:ie,1:je,nnew)
        ENDWHERE
        IF (.NOT.l2tls) THEN
          var(i1,i2,i3)%p3(:,:,nnow) = var(i1,i2,i3)%p3(:,:,nnew)
        ENDIF
      ENDIF
    END SELECT

  ENDDO  variable_loop
         
!------------------------------------------------------------------------------
! End of the subroutine 
!------------------------------------------------------------------------------

END SUBROUTINE id_mix_bd

!==============================================================================
!+ Subroutine for distributing the initial- and the boundary data
!------------------------------------------------------------------------------

SUBROUTINE scatter_data (index, iloc1, iloc2, iloc3, iloclist, irank, idim3,  &
                         leof, lrequired, lchecklist, nlistin, ydata,         &
                         yformat, ntlev, yextension, istat,                   &
                         gribarrays, realarrays)
!------------------------------------------------------------------------------
!
! Description:
!  This subroutine is called within the read-loop over all records (with the 
!  loop index "index"). In this loop up to num_compute (number of compute PE) 
!  records are read and distributed to the processors. Each processor gets a 
!  total record for decoding. After the decoding distribute_subarrays does the 
!  distribution (with scatter_values). It works for all numbers of processors.
!  The distributed subarrays are then copied to the corresponding variables.
!  Only one of the optional variables gribarrays and realarrays is passed,
!  to determine, which presicion is used
!
! Method:
!  Distribute the appropriate part of array gribarrays to the corresponding
!  processors (in parallel mode) or just copy it into subarray (in sequential
!  mode).
!
!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER (KIND=iintegers), INTENT(IN)    ::  &
  iloc1, iloc2, iloc3, iloclist,     &    ! characteristics of the variable
  irank, idim3,                      &    ! characteristics of the variable
  nlistin,                           &    ! number of variables in listin
  index,                             &    ! index of the distribution loop
  ntlev                                   ! layer for boundary variables

LOGICAL,                  INTENT(IN)    ::  &
  leof, lrequired                         ! end of file and required field

LOGICAL,                  INTENT(INOUT) ::  &
  lchecklist (0:ke1,nlistin)              ! for checking that everything is read

CHARACTER (LEN= *),       INTENT(IN)    ::  &
  ydata                                   ! initial or boundary data

CHARACTER (LEN= 4),       INTENT(IN)    ::    &
  yformat                                 ! format that has to be read

CHARACTER (LEN= 1),       INTENT(IN)    ::  &
  yextension                              ! initial or boundary data

! Scalar arguments with intent(out):
INTEGER (KIND=iintegers), INTENT(OUT)   ::  &
  istat                                   ! go on, cycle or exit the loop

! Array arguments with intent(in):
REAL    (KIND=irealgrib), OPTIONAL, INTENT(IN)    ::  &
  gribarrays (ie_max*je_max, num_compute) ! decomposed total field

REAL    (KIND=wp),        OPTIONAL, INTENT(IN)    ::  &
  realarrays (ie_max*je_max, num_compute) ! decomposed total field

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER (KIND=iintegers)   :: implcode, i, j, ij
INTEGER (KIND=iintegers)   :: iz_info(9), iz1, iz2, iz3
CHARACTER (LEN=25)         :: yroutine
CHARACTER (LEN=80)         :: yerrmsg

! Local arrays:
REAL    (KIND=irealgrib)   :: zsubarray_1d_grib(ie_max*je_max)
REAL    (KIND=wp)          :: zsubarray_1d_real(ie_max*je_max)
REAL    (KIND=wp)          :: zsubarray_2d(ie_max,je_max), zbias, zfactor
!
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine scatter_data
!------------------------------------------------------------------------------

  iz_info(:) = -999
  istat    = 0
  implcode = 0
  yroutine = 'scatter_data'

  ! check and distribute the next action
  IF (index == my_cart_id) THEN
    IF ( (.NOT.leof) .AND. (lrequired) ) THEN
      iz_info(1) = 0 ! Processor has data
      iz_info(2) = iloc1
      iz_info(3) = iloc2
      iz_info(4) = iloc3
      iz_info(5) = irank
      iz_info(6) = idim3
      iz_info(7) = iloclist
      IF (TRIM(var(iloc1,iloc2,iloc3)%name) == 'HSURF') THEN
        iz_info(8) = igpsl
        iz_info(9) = jgpsl
      ENDIF
    ELSE
      IF (leof) THEN
        iz_info(1) = -1 ! No data because EOF reached
      ELSE
        iz_info(1) = -2 ! No data because lunreq was set
      ENDIF
    ENDIF
  ENDIF

  IF (num_compute > 1) THEN
    CALL distribute_values (iz_info, 9, index, imp_integers, icomm_cart,   &
                            implcode)
    IF (implcode /= 0) THEN
      istat   = 1
      RETURN
    ENDIF
  ENDIF
  
  ! get some values out of iz_info
  istat = iz_info(1)
  iz1   = iz_info(2)
  iz2   = iz_info(3)
  iz3   = iz_info(4)
  IF (iz_info(8) /= -999) THEN
    igpsl = iz_info(8)
    jgpsl = iz_info(9)
  ENDIF
 
  ! Distribute the records
  IF (istat == 0) THEN
    IF (PRESENT(gribarrays)) THEN
      IF (num_compute > 1) THEN
        CALL scatter_values (gribarrays, zsubarray_1d_grib, ie_max*je_max, &
                             num_compute, imp_grib, index, icomm_cart,     &
                             yerrmsg, implcode)
      ELSE
        zsubarray_1d_grib(:) = gribarrays(:,1)
      ENDIF
    ELSEIF (PRESENT(realarrays)) THEN
      IF (num_compute > 1) THEN
        CALL scatter_values (realarrays, zsubarray_1d_real, ie_max*je_max, &
                             num_compute, imp_reals, index, icomm_cart,    &
                             yerrmsg, implcode)
      ELSE
        zsubarray_1d_real(:) = realarrays(:,1)
      ENDIF
    ENDIF

    IF (implcode /= 0) THEN
      istat = 2
      RETURN
    ENDIF
  
    zbias   = var(iz1,iz2,iz3)%bias
    zfactor = var(iz1,iz2,iz3)%factor

    IF (PRESENT(gribarrays)) THEN
      IF ( ydata == 'boundary' .AND. lbd_frame ) THEN
        ! This is done only in Grib Code
        DO j = 1, je
          DO i = 1, ie
            ij = (j-1) * ie_max + i
            IF (zsubarray_1d_grib(ij) /= undefgrib) THEN
              zsubarray_2d(i,j) = REAL (zsubarray_1d_grib(ij), wp) / zfactor - zbias
            ELSE
              zsubarray_2d(i,j) = undef
            ENDIF
          ENDDO
        ENDDO
      ELSEIF (ydata == 'initial' .OR. ydata == 'boundary') THEN
        ! Scale the field (only in Grib)
        ! transform it to the precision used (Grib and NetCDF)
        IF ((yformat == 'grb1') .OR. (yformat == 'apix')) THEN
          DO j = 1, je
            DO i = 1, ie
              ij = (j-1) * ie_max + i
              zsubarray_2d(i,j) = REAL (zsubarray_1d_grib(ij), wp) / zfactor - zbias
            ENDDO
          ENDDO
        ELSEIF (yformat == 'ncdf') THEN
          IF (TRIM(var(iz1,iz2,iz3)%name) == 'Z0') THEN
            DO j = 1, je
              DO i = 1, ie
                ij = (j-1) * ie_max + i
                zsubarray_2d(i,j) = REAL (zsubarray_1d_grib(ij), wp) / zfactor - zbias
              ENDDO
            ENDDO
          ELSE
            DO j = 1, je
              DO i = 1, ie
                ij = (j-1) * ie_max + i
                zsubarray_2d(i,j) = REAL (zsubarray_1d_grib(ij), wp)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDIF  ! ydata

    ELSEIF (PRESENT(realarrays)) THEN

      DO j = 1, je
        DO i = 1, ie
          ij = (j-1) * ie_max + i
          zsubarray_2d(i,j) = zsubarray_1d_real(ij)
        ENDDO
      ENDDO

    ENDIF ! PRESENT grib/real-arrays

    ! put data into the variables
    IF (ydata == 'initial') THEN
      IF(iz_info(5) == 4 ) THEN
        IF (ASSOCIATED (var(iz1,iz2,iz3)%p4)) THEN
          var(iz1,iz2,iz3)%p4(1:ie,1:je,iz_info(6),nnew) = zsubarray_2d(1:ie,1:je)
          IF(.NOT.l2tls) THEN
            var(iz1,iz2,iz3)%p4(1:ie,1:je,iz_info(6),nnow) = zsubarray_2d(1:ie,1:je)
          ENDIF
        ELSE
          PRINT *, '*** No memory allocated for ', var(iz1,iz2,iz3)%name
          yerrmsg = 'Pointer of variable table is not associated '
          CALL model_abort(my_cart_id, 2011, yerrmsg, yroutine)
        ENDIF
      ELSE IF(iz_info(5) == 3 ) THEN
        SELECT CASE (var(iz1,iz2,iz3)%name)
        CASE ('PS        ','T_G       ','QV_S      ','W_SNOW    ',  &
              'T_S       ','T_M       ','W_G1      ','W_G2      ',  &
              'W_G3      ','W_I       ','T_SNOW    ','RHO_SNOW  ',  &
              'H_SNOW    ','H_ICE     ','T_ICE     ','T_MNW_LK  ',  &
              'T_WML_LK  ','T_BOT_LK  ','T_B1_LK   ','C_T_LK    ',  &
              'H_ML_LK   ','H_B1_LK   ')
          IF (ASSOCIATED (var(iz1,iz2,iz3)%p3)) THEN
            var(iz1,iz2,iz3)%p3(1:ie,1:je,nnew) = zsubarray_2d(1:ie,1:je)
            IF(.NOT.l2tls) THEN
              var(iz1,iz2,iz3)%p3(1:ie,1:je,nnow) = zsubarray_2d(1:ie,1:je)
            ENDIF
          ELSE
            PRINT *, '*** No memory allocated for ', var(iz1,iz2,iz3)%name
            yerrmsg = 'Pointer of variable table is not associated '
            CALL model_abort(my_cart_id, 2011, yerrmsg, yroutine)
          ENDIF
        CASE DEFAULT
          IF (ASSOCIATED (var(iz1,iz2,iz3)%p3)) THEN
            var(iz1,iz2,iz3)%p3(1:ie,1:je,iz_info(6)) = zsubarray_2d(1:ie,1:je)
          ELSE
            PRINT *, '*** No memory allocated for ', var(iz1,iz2,iz3)%name
            yerrmsg = 'Pointer of variable table is not associated '
            CALL model_abort(my_cart_id, 2011, yerrmsg, yroutine)
          ENDIF
        END SELECT
      ELSE
        IF (ASSOCIATED (var(iz1,iz2,iz3)%p2)) THEN
          var(iz1,iz2,iz3)%p2(1:ie,1:je) = zsubarray_2d(1:ie,1:je)
        ELSE
          PRINT *, '*** No memory allocated for ', var(iz1,iz2,iz3)%name
          yerrmsg = 'Pointer of variable table is not associated '
          CALL model_abort(my_cart_id, 2011, yerrmsg, yroutine)
        ENDIF
      ENDIF
    ELSE IF (ydata == 'boundary') THEN
      IF(iz_info(5) == 4 ) THEN
        IF (ASSOCIATED (var(iz1,iz2,iz3)%p4_bd)) THEN
          var(iz1,iz2,iz3)%p4_bd(1:ie,1:je,iz_info(6),ntlev) = zsubarray_2d(1:ie,1:je)
        ELSE
          PRINT *, '*** No memory allocated for ', var(iz1,iz2,iz3)%name
          yerrmsg = 'Pointer of variable table is not associated '
          CALL model_abort(my_cart_id, 2011, yerrmsg, yroutine)
        ENDIF
      ELSEIF(iz_info(5) == 3 ) THEN
        IF (ASSOCIATED (var(iz1,iz2,iz3)%p3_bd)) THEN
          var(iz1,iz2,iz3)%p3_bd(1:ie,1:je,ntlev) = zsubarray_2d(1:ie,1:je)
        ELSE
          PRINT *, '*** No memory allocated for ', var(iz1,iz2,iz3)%name
          yerrmsg = 'Pointer of variable table is not associated '
          CALL model_abort(my_cart_id, 2011, yerrmsg, yroutine)
        ENDIF
      ELSEIF(iz_info(5) == 2 ) THEN
        ! these are the external parameters that have to be read for
        ! climate simulations
        IF (ASSOCIATED (var(iz1,iz2,iz3)%p3_bd)) THEN
          var(iz1,iz2,iz3)%p3_bd(1:ie,1:je,ntlev) = zsubarray_2d(1:ie,1:je)
        ELSE
          PRINT *, '*** No memory allocated for ', var(iz1,iz2,iz3)%name
          yerrmsg = 'Pointer of variable table is not associated '
          CALL model_abort(my_cart_id, 2011, yerrmsg, yroutine)
        ENDIF
      ELSE
        PRINT *, '*** Wrong value for rank of boundary variable ***'  
        PRINT *, '*** only rank=2/3/4 allowed, but rank = ', iz_info(5), '    ***'
        yerrmsg = 'Wrong value for rank of boundary variable'
        CALL model_abort(my_cart_id, 2011, yerrmsg, yroutine)
      ENDIF
    ELSEIF ( (ydata == 'restart') .AND. (yextension == 'o') ) THEN
      IF(iz_info(5) == 4 ) THEN
        IF (ASSOCIATED (var(iz1,iz2,iz3)%p4)) THEN
          IF (var(iz1,iz2,iz3)%name == 'TKE       ') THEN
            var(iz1,iz2,iz3)%p4(1:ie,1:je,iz_info(6),ntke) = zsubarray_2d(1:ie,1:je)
          ELSE
            IF (l2tls) THEN
              var(iz1,iz2,iz3)%p4(1:ie,1:je,iz_info(6),nnew) = zsubarray_2d(1:ie,1:je)
            ELSE
              var(iz1,iz2,iz3)%p4(1:ie,1:je,iz_info(6),nnow) = zsubarray_2d(1:ie,1:je)
            ENDIF
          ENDIF
        ELSE
          PRINT *, '*** No memory allocated for ', var(iz1,iz2,iz3)%name
          yerrmsg = 'Pointer of variable table is not associated '
          CALL model_abort(my_cart_id, 2011, yerrmsg, yroutine)
        ENDIF
      ELSE IF(iz_info(5) == 3 ) THEN
        SELECT CASE (var(iz1,iz2,iz3)%name)
        CASE ('PS        ','T_G       ','QV_S      ','W_SNOW    ',  &
              'T_S       ','T_M       ','W_G1      ','W_G2      ',  &
              'W_G3      ','W_I       ','T_SNOW    ','RHO_SNOW  ',  &
              'H_SNOW    ','H_ICE     ','T_ICE     ','T_MNW_LK  ',  &
              'T_WML_LK  ','T_BOT_LK  ','T_B1_LK   ','C_T_LK    ',  &
              'H_ML_LK   ','H_B1_LK   ')
          IF (ASSOCIATED (var(iz1,iz2,iz3)%p3)) THEN
            IF (l2tls) THEN
              var(iz1,iz2,iz3)%p3(1:ie,1:je,nnew) = zsubarray_2d(1:ie,1:je)
            ELSE
              var(iz1,iz2,iz3)%p3(1:ie,1:je,nnow) = zsubarray_2d(1:ie,1:je)
            ENDIF
          ELSE
            PRINT *, '*** No memory allocated for ', var(iz1,iz2,iz3)%name
            yerrmsg = 'Pointer of variable table is not associated '
            CALL model_abort(my_cart_id, 2011, yerrmsg, yroutine)
          ENDIF
        CASE DEFAULT
          IF (ASSOCIATED (var(iz1,iz2,iz3)%p3)) THEN
            var(iz1,iz2,iz3)%p3(1:ie,1:je,iz_info(6)) = zsubarray_2d(1:ie,1:je)
          ELSE
            PRINT *, '*** No memory allocated for ', var(iz1,iz2,iz3)%name
            yerrmsg = 'Pointer of variable table is not associated '
            CALL model_abort(my_cart_id, 2011, yerrmsg, yroutine)
          ENDIF
        END SELECT
      ELSE
        IF (ASSOCIATED (var(iz1,iz2,iz3)%p2)) THEN
          var(iz1,iz2,iz3)%p2(1:ie,1:je) = zsubarray_2d(1:ie,1:je)
        ELSE
          PRINT *, '*** No memory allocated for ', var(iz1,iz2,iz3)%name
          yerrmsg = 'Pointer of variable table is not associated '
          CALL model_abort(my_cart_id, 2011, yerrmsg, yroutine)
        ENDIF
      ENDIF
    ELSEIF ( (ydata == 'restart') .AND. (yextension == 'n') ) THEN
      IF(iz_info(5) == 4 ) THEN
        IF (ASSOCIATED (var(iz1,iz2,iz3)%p4)) THEN
          IF (var(iz1,iz2,iz3)%name == 'TKE       ') THEN
            IF (ntke == 1) THEN
              ! level 1 already set with o-restart file
              var(iz1,iz2,iz3)%p4(1:ie,1:je,iz_info(6),   3) = zsubarray_2d(1:ie,1:je)
              var(iz1,iz2,iz3)%p4(1:ie,1:je,iz_info(6),   2) = zsubarray_2d(1:ie,1:je)
            ELSE
              var(iz1,iz2,iz3)%p4(1:ie,1:je,iz_info(6),ntke-1) = zsubarray_2d(1:ie,1:je)
            ENDIF
          ELSE
            var(iz1,iz2,iz3)%p4(1:ie,1:je,iz_info(6),nnew) = zsubarray_2d(1:ie,1:je)
          ENDIF
        ELSE
          PRINT *, '*** No memory allocated for ', var(iz1,iz2,iz3)%name
          yerrmsg = 'Pointer of variable table is not associated '
          CALL model_abort(my_cart_id, 2011, yerrmsg, yroutine)
        ENDIF
      ELSE IF(iz_info(5) == 3 ) THEN
        SELECT CASE (var(iz1,iz2,iz3)%name)
        CASE ('PS        ','T_G       ','QV_S      ','W_SNOW    ',  &
              'T_S       ','T_M       ','W_G1      ','W_G2      ',  &
              'W_G3      ','W_I       ','T_SNOW    ','RHO_SNOW  ',  &
              'H_SNOW    ','H_ICE     ','T_ICE     ','T_MNW_LK  ',  &
              'T_WML_LK  ','T_BOT_LK  ','T_B1_LK   ','C_T_LK    ',  &
              'H_ML_LK   ','H_B1_LK   ')
          IF (ASSOCIATED (var(iz1,iz2,iz3)%p3)) THEN
            var(iz1,iz2,iz3)%p3(1:ie,1:je,nnew) = zsubarray_2d(1:ie,1:je)
          ELSE
            PRINT *, '*** No memory allocated for ', var(iz1,iz2,iz3)%name
            yerrmsg = 'Pointer of variable table is not associated '
            CALL model_abort(my_cart_id, 2011, yerrmsg, yroutine)
          ENDIF
        CASE DEFAULT
          IF (ASSOCIATED (var(iz1,iz2,iz3)%p3)) THEN
            var(iz1,iz2,iz3)%p3(1:ie,1:je,iz_info(6)) = zsubarray_2d(1:ie,1:je)
          ELSE
            PRINT *, '*** No memory allocated for ', var(iz1,iz2,iz3)%name
            yerrmsg = 'Pointer of variable table is not associated '
            CALL model_abort(my_cart_id, 2011, yerrmsg, yroutine)
          ENDIF
        END SELECT
      ELSE
        IF (ASSOCIATED (var(iz1,iz2,iz3)%p2)) THEN
          var(iz1,iz2,iz3)%p2(1:ie,1:je) = zsubarray_2d(1:ie,1:je)
        ELSE
          PRINT *, '*** No memory allocated for ', var(iz1,iz2,iz3)%name
          yerrmsg = 'Pointer of variable table is not associated '
          CALL model_abort(my_cart_id, 2011, yerrmsg, yroutine)
        ENDIF
      ENDIF
    ENDIF 

    ! Set the corresponding check variable to true:
    lchecklist (iz_info(6), iz_info(7)) = .TRUE.

#ifdef _OPENACC
    ! We need to save the list of fields which are read from boundary
    ! condition files in order to be able to save them to GPU
    ! setup entry in bd_list structure (only for 'boundary' mode)
    IF (trim(ydata) == 'boundary') THEN
      ! check if this field is already present
      DO i = 1, num_bd_fields
        IF (bd_list(i)%iloc1 == iz1 .AND. &
            bd_list(i)%iloc2 == iz2 .AND. &
            bd_list(i)%iloc3 == iz3 .AND. &
            bd_list(i)%rank  == iz_info(5) .AND. &
            bd_list(i)%ntlev == ntlev ) EXIT
      ENDDO
      IF (i > max_bd_fields) THEN
        yerrmsg = 'ERROR *** Increase max_bd_fields in data_io.f90!!! *** ERROR'
        CALL model_abort(my_cart_id, 5432, yerrmsg, yroutine)
      ENDIF
      IF (i > num_bd_fields) THEN
        num_bd_fields = i
        bd_list(i)%name  = TRIM(var(iz1,iz2,iz3)%name)
        bd_list(i)%rank  = iz_info(5)
        bd_list(i)%iloc1 = iz1
        bd_list(i)%iloc2 = iz2
        bd_list(i)%iloc3 = iz3
        bd_list(i)%ntlev = ntlev
        IF (iz_info(5) == 4 ) THEN
          ! rank 4 fields with timelevel
          bd_list(i)%p3 => var(iz1,iz2,iz3)%p4_bd(:,:,:,ntlev)
          bd_list(i)%p2 => NULL()
        ELSEIF (iz_info(5) == 3 ) THEN
          ! rank 3 field with timelevel
          bd_list(i)%p3 => NULL()
          bd_list(i)%p2 => var(iz1,iz2,iz3)%p3_bd(:,:,ntlev)
        ELSEIF (iz_info(5) == 2 ) THEN
          ! rank 2 field with timelevel
          bd_list(i)%p3 => NULL()
          bd_list(i)%p2 => var(iz1,iz2,iz3)%p3_bd(:,:,ntlev)
        ELSE
          yerrmsg = 'ERROR *** Unsupported field type encountered in src_input.f90!!! *** ERROR'
          CALL model_abort(my_cart_id, 5433, yerrmsg, yroutine)
        ENDIF
        !IF ( my_cart_id == 0 ) THEN
        !  WRITE(*,*) ' GPUINFO: Tracked BD-field input ' // &
        !    TRIM(bd_list(i)%name), bd_list(i)%rank, bd_list(i)%ntlev
        !ENDIF
      ENDIF
    ENDIF
#endif

  ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE scatter_data

!==============================================================================
!==============================================================================
!+ Module procedure in "io_utilities" for writing a NetCDF file
!------------------------------------------------------------------------------

SUBROUTINE read_nc_gdefs (ncid, iedim, jedim, kedim, kesoildim, ydate, &
                           startlat, startlon, dlon, dlat,  &
                           icomm, myid, npes, yerrmsg, istatus)
!
!------------------------------------------------------------------------------
!
! Description:
!   This routine reads global attributes from NetCDF formatted input file.
!
!------------------------------------------------------------------------------
!
!
! Subroutine / Function arguments

! Scalar arguments with intent(in):
  INTEGER (KIND=iintegers),   INTENT(IN) :: &
    ncid,              & ! NetCDF file ID
    icomm,             & ! MPI communicator
    myid,              & ! ID of this PE in communicator icomm
    npes,              & ! number of PEs
    iedim, jedim, kedim, kesoildim  ! dimensions of the input fields

  CHARACTER (LEN=*),        INTENT(IN)  ::  &
    ydate               ! actual date from Namelist parameter

  REAL    (KIND=wp),        INTENT(IN)  ::  &
    startlat, startlon, & ! coordinates of the lower left grid point
    dlon,  dlat           ! grid resolution in lambda and phi direction


! Scalar arguments with intent(out):
  INTEGER (KIND=iintegers), INTENT(OUT)     ::  &
    istatus          ! error index

  CHARACTER (LEN= *),        INTENT(OUT)   ::  &
    yerrmsg     ! string for error messages

!-----------------------------------------------------------------------
!
! Local scalars:

  INTEGER (KIND=iintegers), PARAMETER  :: ntime=1

  INTEGER (KIND=iintegers) ::  &
    i, j,                                    & ! loop variable
    ierror, izmplcode,                       & ! error status variable
    iedim_in, jedim_in, kedim_in, ke1dim_in, & ! input dimensions
    kesoildim_in, ke1soildim_in,             & ! input dimension of multi soil layers
    nhori_in,                                & ! input dimension of nhori
    ntime_in, ntbds_in,                      & !
    idb                                        ! index for finding a string

  INTEGER (KIND=iintegers)::   &
    jgridVarID,     & ! NetCDF ID for grid mapping
    jlonVarID,      & ! NetCDF ID for longitude
    jlatVarID,      & ! NetCDF ID for latitude
    jvcVarID ,      & ! NetCDF ID for the vertical component
    jsoilVarID,     & ! NetCDF ID for the multi soil layer component
    jsectVarID,     & ! NetCDF ID for the topo. correction field
    jtimeID           ! NetCDF ID for the time

  INTEGER (KIND=iintegers)   ::                                              &
    iyear_ref, imonth_ref, iday_ref, ihour_ref, imin_ref, isec_ref,          &
    iyear, imonth, iday, ihour, imin, isec, iref_act, ierrf, itimepassed
    
  CHARACTER (LEN=40) :: &
    ydate1           ! time unit and reference date from input

  CHARACTER (LEN=14) ::  &
    ydate_ref        ! reference date from input

  CHARACTER (LEN=20)  ::  &
    calendar         ! type of calendar

  REAL (KIND=irealgrib) :: &
    zpollat,       & ! latitude of North Pole in input data
    zpollon,       & ! longitude of North Pole in input data
    zpolgam          ! angle between the north poles of the systems

  CHARACTER (LEN=80) :: &
    grid_name        ! name of grid mapping

  REAL    (KIND=irealgrib)     ::  &
    zdlon_in, zdlat_in      ! grid resolution in latitude and longitude

  REAL    (KIND=irealgrib)    :: p0slgrib, t0slgrib, dt0lpgrib, vcflatgrib, &
                                 deltatgrib, hscalgrib, svc1grib, svc2grib !_br 25.11.08


!  Local arrays:
  
  REAL (KIND=wp)     :: &
    time(ntime)      ! forecast time

  REAL (KIND=wp),        ALLOCATABLE   :: &
    zvc_params(:)    ! vertical coordinate

  REAL (KIND=irealgrib), ALLOCATABLE   :: &
    vcoordgrib(:), & ! vertical coordinate
    longitude(:),  & ! rotated longitudes
    latitude(:),   & ! rotated latitudes
    zczml_soil(:)    ! height layers of multi layer soil model

!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

#ifdef NETCDF
  ! initialize error status variable
  yerrmsg       = '   '
  ierror        = 0
  istatus       = 0

! Processor 0 does the job
  IF (myid == 0) THEN
! Get the dimension ID's and the length of the dimensions
  istatus = nf90_inq_dimid (ncid, 'rlon', idims_id_in(1))
  IF (istatus /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(istatus))
    RETURN
  ENDIF
  istatus = nf90_inquire_dimension (ncid, idims_id_in(1), len=iedim_in)
  IF (istatus /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(istatus))
    RETURN
  ENDIF
  IF (iedim /= iedim_in) ierror = 1

  istatus = nf90_inq_dimid (ncid, 'rlat', idims_id_in(2))
  IF (istatus /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(istatus))
    RETURN
  ENDIF
  istatus = nf90_inquire_dimension (ncid, idims_id_in(2), len=jedim_in)
  IF (istatus /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(istatus))
    RETURN
  ENDIF
  IF (jedim /= jedim_in) ierror = 2

  istatus = nf90_inq_dimid (ncid, 'level', idims_id_in(3))
  IF (istatus /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(istatus))
    RETURN
  ENDIF
  istatus = nf90_inquire_dimension (ncid, idims_id_in(3), len=kedim_in)
  IF (istatus /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(istatus))
    RETURN
  ENDIF
  IF (kedim /= kedim_in) ierror = 3

  istatus = nf90_inq_dimid (ncid, 'level1', idims_id_in(4))
  IF (istatus /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(istatus))
    RETURN
  ENDIF
  istatus = nf90_inquire_dimension (ncid, idims_id_in(4), len=ke1dim_in)
  IF (istatus /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(istatus))
    RETURN
  ENDIF
  IF ((kedim+1) /= ke1dim_in) ierror = 3

  ALLOCATE (longitude(iedim_in), latitude(jedim_in), vcoordgrib(ke1dim_in), &
            zvc_params(ke1dim_in), STAT=istatus)
  IF (istatus /= 0) THEN    
    yerrmsg = 'Allocation error in read_nc_gdefs'
    RETURN
  ENDIF

  istatus = nf90_inq_dimid (ncid, 'time', idims_id_in(5))
  IF (istatus /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(istatus))
    RETURN
  ENDIF
  istatus = nf90_inquire_dimension (ncid, idims_id_in(5), len=ntime_in)
  IF (istatus /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(istatus))
    RETURN
  ENDIF
  IF (ntime_in /= 1) ierror = 4

  istatus = nf90_inq_dimid (ncid, 'bnds', idims_id_in(6))
  IF (istatus /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(istatus))
    RETURN
  ENDIF
  istatus = nf90_inquire_dimension (ncid, idims_id_in(6), len=ntbds_in)
  IF (istatus /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(istatus))
    RETURN
  ENDIF
  IF (ntbds_in /= 2) ierror = 4

  IF (lmulti_layer) THEN
  
 !   istatus = nf90_inq_dimid(ncid, "soil", idims_id_in(7))
 !   IF (istatus /= NF90_NOERR) THEN
 !     yerrmsg = TRIM(NF90_strerror(istatus))
 !     RETURN
 !   ENDIF
 !   istatus = nf90_inquire_dimension (ncid, idims_id_in(7), len=kesoildim_in)
 !   IF (istatus /= NF90_NOERR) THEN
 !     yerrmsg = TRIM(NF90_strerror(istatus))
 !     RETURN
 !   ENDIF
 !   IF (kesoildim /= kesoildim_in) ierror = 5

     istatus=nf90_inq_dimid(ncid, "soil1", idims_id_in(8))
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    istatus = nf90_inquire_dimension (ncid, idims_id_in(8), len=ke1soildim_in)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    kesoildim_in = ke1soildim_in - 1
    IF (kesoildim+1 /= ke1soildim_in) ierror = 5
    
    ! get the values of the soil layers
    istatus = nf90_inq_varid (ncid, 'soil1', jsoilVarID)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF

    ALLOCATE (zczml_soil(ke1soildim_in), STAT=istatus)
    IF (istatus /= 0) THEN
      yerrmsg = 'Allocation error in read_nc_gdefs'
      RETURN
    ENDIF

    istatus = nf90_get_var (ncid, jsoilVarID, zczml_soil)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    

    istatus = 0
    DO j = 1, ke1soildim_in
      IF (ABS(zczml_soil(j) - czmls(j)) > 0.001_wp) istatus = -1
    ENDDO
    IF (istatus /= 0) THEN
      PRINT *, 'soil layer     data file     namelist input'
      DO j = 1, ke1soildim_in
        WRITE(*,'(I8,2F14.3)') j, zczml_soil(j), czmls(j)
      ENDDO
      istatus = 30
      yerrmsg = 'Error in multi soil layer heights'
      RETURN
    ENDIF

    DEALLOCATE (zczml_soil)

 ENDIF

! Get dimension for 3D external parameter field for topographical corrections

 IF(lradtopo) THEN
    istatus = nf90_inq_dimid (ncid, "nhori", idims_id_in(11))
    IF (istatus /= NF90_NOERR) THEN
       yerrmsg = TRIM(NF90_strerror(istatus))
       RETURN
    ENDIF
    istatus = nf90_inquire_dimension (ncid, idims_id_in(11), len=nhori_in)
    IF (istatus /= NF90_NOERR) THEN
       yerrmsg = TRIM(NF90_strerror(istatus))
       RETURN
    ENDIF

    IF ((nhori) /= nhori_in) THEN
      ierror = 11
    ELSE
      istatus = nf90_inq_varid (ncid, "nhori", jsectVarID)
      IF (istatus /= NF90_NOERR) THEN
         yerrmsg = TRIM(NF90_strerror(istatus))
         RETURN
      ENDIF

    ENDIF
 ENDIF

! Get longitude and latitude data

  istatus = nf90_inq_varid (ncid, 'rlon', jlonVarID)
  IF (istatus /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(istatus))
    RETURN
  ENDIF
  istatus = nf90_get_var (ncid, jlonVarID, longitude)  
  IF (istatus /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(istatus))
    RETURN
  ENDIF

  istatus = nf90_inq_varid (ncid, 'rlat', jlatVarID)
  IF (istatus /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(istatus))
    RETURN
  ENDIF
  istatus = nf90_get_var (ncid, jlatVarID, latitude)  
  IF (istatus /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(istatus))
    RETURN
  ENDIF

  zdlat_in   = (latitude(jedim_in) - latitude(1))/ REAL(jedim_in-1,irealgrib)
  IF (longitude(iedim_in) - longitude(1) < 0.0_wp) THEN
    ! If the area is located around the 180-Meridian, longitude(iedim_in) - longitude(1)
    ! will be negative and 360 degrees have to be added
    zdlon_in   = (longitude(iedim_in) - longitude(1) + 360.0_irealgrib) / (iedim_in-1)
  ELSE
    zdlon_in   = (longitude(iedim_in) - longitude(1))/ (iedim_in-1)
  ENDIF

! get grid mapping values

  istatus = nf90_inq_varid(ncid, 'rotated_pole', jgridVarid)
  IF (istatus == NF90_NOERR) THEN   ! rotated latitude-longitude
    istatus = nf90_get_att(ncid, jgridVarid, 'grid_mapping_name', grid_name)
    IF (istatus /= NF90_NOERR) THEN
      PRINT *, 'Error in read_nc_gdefs / nf90_get_att'
      PRINT *, 'Attribute "grid_mapping_name"'
      yerrmsg = TRIM(NF90_strerror(istatus))
    ENDIF
    IF (grid_name(1:26) /= 'rotated_latitude_longitude') THEN
      PRINT *, 'Error in read_nc_gdefs'
      PRINT *, 'Invalid value for attribute "grid_mapping_name"'
      yerrmsg = TRIM(grid_name)
    ENDIF
    istatus = nf90_get_att(ncid, jgridVarid, 'grid_north_pole_latitude', zpollat)
    IF (istatus /= NF90_NOERR) THEN
      PRINT *, 'Error in read_nc_gdefs / nf90_get_att'
      PRINT *, 'Attribute "grid_north_pole_latitude"'
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    istatus = nf90_get_att(ncid, jgridVarid, 'grid_north_pole_longitude', zpollon)
    IF (istatus /= NF90_NOERR) THEN
      PRINT *, 'Error in read_nc_gdefs / nf90_get_att'
      PRINT *, 'Attribute "grid_north_pole_longitude"'
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    istatus = nf90_get_att(ncid, jgridVarid, 'north_pole_grid_longitude', zpolgam)
    IF (istatus == NF90_ENOTATT) THEN
      zpolgam = 0.0_irealgrib
    ELSE IF (istatus /= NF90_NOERR) THEN
      PRINT *, 'Error in read_nc_gdefs / nf90_get_att'
      PRINT *, 'Attribute "north_pole_grid_longitude"'
      PRINT *, TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
  ELSE IF (istatus == NF90_ENOTVAR) THEN  ! true latitude-longitude
    zpollon =  0.0_irealgrib
    zpollat = 90.0_irealgrib
    zpolgam =  0.0_irealgrib
  ELSE
     PRINT *, 'Error in read_nc_gdefs / nf90_inq_varid'
     yerrmsg = 'Variable "rotated_pole"'
     RETURN
  ENDIF

! compare the input values with the Namelist parameters

  IF (ABS(REAL(zpollat,wp) - pollat) > 1.0E-3_wp) ierror = 4
  IF (ABS(REAL(zpollon,wp) - pollon) > 1.0E-3_wp) ierror = 5
  IF (ABS(REAL(zpolgam,wp) - polgam) > 1.0E-3_wp) ierror = 5

  IF (ABS(REAL(latitude(1),wp)  - startlat) > 1.0E-3_wp) ierror = 6
  IF (ABS(REAL(longitude(1),wp)  - startlon) > 1.0E-3_wp) ierror = 7

  IF (ABS(REAL(zdlat_in,wp) - dlat) > 1.0E-3_wp) ierror = 8
  IF (ABS(REAL(zdlon_in,wp) - dlon) > 1.0E-3_wp) ierror = 9

! Get the vertical co-ordinate

  istatus = nf90_inq_varid (ncid, 'vcoord', jvcVarID)
  IF (istatus /= NF90_NOERR) THEN
    yerrmsg = 'vcoord '// TRIM(NF90_strerror(istatus))
    RETURN
  ENDIF
  istatus = nf90_get_var (ncid, jvcVarID, vcoordgrib)
  IF (istatus /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(istatus))
    RETURN
  ENDIF
  zvc_params = REAL(vcoordgrib, wp)
  istatus = nf90_get_att (ncid, jvcVarID, 'p0sl', p0slgrib)
  IF (istatus /= NF90_NOERR) THEN
    yerrmsg = 'p0sl '//TRIM(NF90_strerror(istatus))
    RETURN
  ENDIF
  refatm%p0sl = REAL (p0slgrib, wp)
  istatus = nf90_get_att (ncid, jvcVarID, 't0sl', t0slgrib)
  IF (istatus /= NF90_NOERR) THEN
    yerrmsg = 't0sl '//TRIM(NF90_strerror(istatus))
    RETURN
  ENDIF
  refatm%t0sl = REAL (t0slgrib, wp)
  istatus = nf90_get_att (ncid, jvcVarID, 'dt0lp', dt0lpgrib)
  IF (istatus /= NF90_NOERR) THEN
    yerrmsg = 'dt0lp '//TRIM(NF90_strerror(istatus))
    RETURN
  ENDIF
  refatm%dt0lp = REAL (dt0lpgrib, wp)
  istatus = nf90_get_att (ncid, jvcVarID, 'vcflat', vcflatgrib)
  IF (istatus /= NF90_NOERR) THEN
    yerrmsg = 'vcflat '//TRIM(NF90_strerror(istatus))
    RETURN
  ENDIF
  vcoord%vcflat = REAL (vcflatgrib, wp)

  istatus = nf90_get_att (ncid, jvcVarID, 'irefatm', refatm%irefatm)
  IF (istatus == NF90_NOERR) THEN
    istatus = nf90_get_att (ncid, jvcVarID, 'ivctype', vcoord%ivctype)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = 'ivctype '//TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    IF (refatm%irefatm == 2) THEN ! reference atmosphere type 2
      istatus = nf90_get_att (ncid, jvcVarID, 'delta_t', deltatgrib)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'delta_t '//TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
      refatm%delta_t = REAL(deltatgrib, wp)
      istatus = nf90_get_att (ncid, jvcVarID, 'h_scal', hscalgrib)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'h_scal '//TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
      refatm%h_scal = REAL(hscalgrib, wp)
    ENDIF
    IF ( ANY( vcoord%ivctype == (/3,4/) ) ) THEN ! SLEVE coordinates
     istatus = nf90_get_att (ncid, jvcVarID, 'svc1', svc1grib)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'svc1 '//TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
      svc1 = REAL(svc1grib, wp)
      istatus = nf90_get_att (ncid, jvcVarID, 'svc2', svc2grib)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'svc2 '//TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
      svc2 = REAL(svc2grib, wp)
      istatus = nf90_get_att (ncid, jvcVarID, 'nfltvc', nfltvc)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'nfltvc '//TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
    ENDIF
  ELSE ! old version (i.e. int2lm output without irefatm and ivctype)
    refatm%irefatm = 1
    vcoord%ivctype = 1
  ENDIF

  vcoord%nlevels           = ke1dim_in
  IF     (vcoord%ivctype == 1) THEN
    vcoord%sigm_coord(1:ke1dim_in) = zvc_params(1:ke1dim_in)
  ELSEIF ( ANY( vcoord%ivctype == (/2,3,4/) ) ) THEN
    vcoord%vert_coord(1:ke1dim_in) = zvc_params(1:ke1dim_in)
  ENDIF

! Get time data

  istatus = nf90_inq_varid (ncid, 'time', jtimeID)
  IF (istatus /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(istatus))
    RETURN
  ENDIF
  istatus = nf90_get_var (ncid, jtimeID, time)  
  IF (istatus /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(istatus))
    RETURN
  ENDIF
  istatus = nf90_get_att (ncid, jtimeID, 'calendar', calendar)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF
    IF (((calendar(1:8) == 'standard' .OR. calendar(1:19) == 'proleptic_gregorian') &
                                     .AND. (itype_calendar == 0)) .OR. &
         (calendar(1:7) == '360_day' .AND. (itype_calendar == 1)) .OR.              &
         (calendar(1:7) == '365_day' .AND. (itype_calendar == 2))) THEN
       ! callendar attribute is valid
    ELSE
!US    IF (.NOT.(((calendar(1:8) == 'standard'.OR.calendar(1:19) == 'proleptic_gregorian') &
!US        .AND. .NOT. lyear_360) .OR. (calendar(1:7) == '360_day'  .AND. lyear_360))) THEN
      IF (myid == 0) THEN
        PRINT *, 'calendar attribute = >',calendar,'<  '
        PRINT *, ' but wrong itype_calendar = ', itype_calendar
      ENDIF
      yerrmsg = ' ERROR *** Wrong calendar attribute ***'
      istatus = -1
      RETURN
    ENDIF

  istatus = nf90_get_att (ncid, jtimeID, 'units', ydate1)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF

    ! check the date for which the product is valid
      idb = index(ydate1,'since') + 5 
      DO i = 1,LEN_TRIM(ydate1)
        IF(ICHAR(ydate1(i:i)) < 48 .OR. ICHAR(ydate1(i:i)) > 57) ydate1(i:i)=' '
      ENDDO
        READ(ydate1(idb:),*,IOSTAT=istatus) iyear_ref,imonth_ref,iday_ref,ihour_ref,imin_ref,isec_ref
        IF (istatus /= 0) THEN
          yerrmsg = "ERROR in units attribute of time in the netCDF input"
          RETURN
        ENDIF
      IF (lmmss) THEN
        WRITE(ydate_ref,'(I4.4,5(I2.2))') iyear_ref,imonth_ref,iday_ref,ihour_ref,imin_ref,isec_ref
      ELSE
        WRITE(ydate_ref,'(I4.4,5(I2.2))') iyear_ref,imonth_ref,iday_ref,ihour_ref
      ENDIF

      ! Determine the corresponding integer values for the actual date "ydate"
      ! no minutes are given here
      IF (lmmss) THEN
        READ(ydate    ,'(I4,5I2)') iyear, imonth, iday, ihour, imin, isec
        imin_ref = 0 
        isec_ref = 0 
      ELSE
        READ(ydate    ,'(I4,3I2)') iyear, imonth, iday, ihour
        imin = 0
        isec = 0
      ENDIF
  
      ! Determine the difference between reference date and actual date
      ! in minutes
      CALL diff_minutes(iyear_ref, imonth_ref, iday_ref, ihour_ref, imin_ref, &
                        iyear,     imonth,     iday,     ihour,     imin,     &
                        itype_calendar, iref_act,  ierrf)

      ! add difference in seconds
      iref_act = iref_act * 60_iintegers + isec - isec_ref

      ! Now determine the time passed (since the reference date) in seconds
      ! US: itimepassed is iintegers, but time is real: have to convert: use NINT
      itimepassed = NINT (time(1), iintegers)

      ! The timepassed and the difference of the dates must now be equal:
      IF (iref_act /= itimepassed) THEN
        ierror = 10
      ENDIF

      IF (ierror  == 10) THEN
        WRITE(*,'(A,I4.4,5(A,I2.2))') ' The actual date           ', &
             iyear,'-',imonth,'-',iday,' ',ihour,':',imin,':',isec
        WRITE(*,'(A,I4.4,5(A,I2.2))') ' and the reference date    ', &
             iyear_ref,'-',imonth_ref,'-',iday_ref,' ',ihour_ref,':',imin_ref,':',isec_ref
        PRINT *, 'do not match!!'
        PRINT *, 'time difference is ', iref_act,' s, but should be ',time(1), ' s'
        istatus = ierror
        PRINT *, 'istatus = ',istatus
        yerrmsg = ' The actual date and the reference date do not match !!'
        RETURN
      ENDIF

!_br 25.11.08
  ! Check for the type and consistency of the vertical coordinate parameters
  IF (vcoord%ivctype == 1) THEN
    ! For this type the vertical coordinates should be ascending
    IF ( vcoord%sigm_coord(2) < vcoord%sigm_coord(1) ) THEN
       WRITE (yerrmsg,'(A,I5,2F10.5)')                                       &
            ' ERROR *** Vertical coordinates not ascending for type *** ',   &
            vcoord%ivctype, vcoord%sigm_coord(1), vcoord%sigm_coord(2)
       istatus = 2008
       RETURN
    ENDIF

  ELSEIF (vcoord%ivctype == 2) THEN

  ! For this type the vertical coordinates should be descending
    IF ( vcoord%vert_coord(2) > vcoord%vert_coord(1) ) THEN
       WRITE (yerrmsg,'(A,I5,2F10.5)')                                       &
            ' ERROR *** Vertical coordinates not descending for type *** ',  &
            vcoord%ivctype, vcoord%vert_coord(1), vcoord%vert_coord(2)
       istatus = 2008
       RETURN
    ENDIF

  ELSEIF ( ANY( vcoord%ivctype == (/3,4/) ) ) THEN

   ! Check for meaningful values of svc1, svc2 and nfltvc
    IF ((svc1 > vcoord%vert_coord(1)) .OR. (svc1 < 0.0_wp)) THEN
       yerrmsg = ' ERROR *** svc1 not in allowed range for ivctype = 3/4 ***'
       istatus = 2008
       RETURN
    ENDIF

    IF ((svc2 > vcoord%vert_coord(1)) .OR. (svc2 < 0.0_wp)) THEN
       yerrmsg = ' ERROR *** svc2 not in allowed range for ivctype = 3/4 ***'
       istatus = 2008
       RETURN
    ENDIF

    IF (nfltvc <= 0) THEN
       yerrmsg = ' ERROR *** nfltvc must be greater than or equal to zero ***'
       istatus = 2008
       RETURN
    ENDIF
  ELSE
    yerrmsg = ' ERROR *** Type ivctype of vertical coordinate not available***'
    istatus = 20
    RETURN
  ENDIF
!_br 25.11.08 end

!------------------------------------------------------------------------------
! Section 4: If necessary, print errors
!------------------------------------------------------------------------------

        IF (ierror /= 0 .and. ierror <= 11) THEN
          PRINT *, '         data file          namelist input'

          PRINT *, 'ie_tot                ',iedim_in     ,'       ',iedim
          PRINT *, 'je_tot                ',jedim_in     ,'       ',jedim
          PRINT *, 'ke_tot                ',kedim_in     ,'       ',kedim
          PRINT *, 'ke_tot + 1            ',ke1dim_in    ,'       ',kedim+1

          IF (lmulti_layer) THEN
            PRINT *, 'ke_soil               ',kesoildim_in ,'       ',kesoildim
            PRINT *, 'ke_soil + 1           ',ke1soildim_in,'       ',kesoildim+1
          ENDIF

          IF (lradtopo) THEN
             PRINT *, 'nhori                ',nhori_in     ,'       ',nhori
          ENDIF

          PRINT *, 'startlat_tot          ',latitude(1)  ,'       ',startlat
          PRINT *, 'startlon_tot          ',longitude(1) ,'       ',startlon

          PRINT *, 'dlat                  ',zdlat_in     ,'       ',dlat
          PRINT *, 'dlon                  ',zdlon_in     ,'       ',dlon

          PRINT *, 'pollat                ',zpollat,'       ', pollat
          PRINT *, 'pollon                ',zpollon,'       ', pollon
          PRINT *, 'polgam                ',zpolgam,'       ', polgam

          istatus = ierror
          RETURN
        ENDIF

  DEALLOCATE (latitude, longitude, vcoordgrib, STAT=istatus)
  IF (istatus /= 0) THEN
    yerrmsg = 'Deallocation error in read_nc_gdefs'
    istatus = 21
    RETURN
  ENDIF

  ELSE   ! myid /= 0

    ! the other PEs have to allocate zvc_params also
    ! (but here with ke1 only, ke1dim_in is not known to other PEs)
    ALLOCATE (zvc_params(ke1), STAT=istatus)

  ENDIF  ! myid == 0

! send information to the other processors
  IF (npes > 1) THEN
    CALL MPI_BCAST (vcoord%ivctype,  1, imp_integers, 0, icomm, izmplcode)
    CALL MPI_BCAST (refatm%irefatm,  1, imp_integers, 0, icomm, izmplcode) !_br 25.11.08
    CALL MPI_BCAST (refatm%p0sl,     1, imp_reals,    0, icomm, izmplcode)
    CALL MPI_BCAST (refatm%t0sl,     1, imp_reals,    0, icomm, izmplcode)
    CALL MPI_BCAST (refatm%dt0lp,    1, imp_reals,    0, icomm, izmplcode)
    CALL MPI_BCAST (vcoord%vcflat,   1, imp_reals,    0, icomm, izmplcode)
    IF (refatm%irefatm == 2) THEN
      CALL MPI_BCAST (refatm%delta_t,  1, imp_reals,    0, icomm, izmplcode)
      CALL MPI_BCAST (refatm%h_scal,   1, imp_reals,    0, icomm, izmplcode)
    ENDIF

    CALL MPI_BCAST (zvc_params, ke1, imp_reals,    0, icomm, izmplcode)
    vcoord%nlevels           = ke1
    IF     (vcoord%ivctype == 1) THEN
      vcoord%sigm_coord(1:ke1) = zvc_params(1:ke1)
    ELSEIF ( ANY( vcoord%ivctype == (/2,3,4/) ) ) THEN
      vcoord%vert_coord(1:ke1) = zvc_params(1:ke1)
    ENDIF

    IF ( ANY( vcoord%ivctype == (/3,4/) ) ) THEN  ! SLEVE coordinate
      CALL MPI_BCAST (svc1,     1, imp_reals,    0, icomm, izmplcode) 
      CALL MPI_BCAST (svc2,     1, imp_reals,    0, icomm, izmplcode)
      CALL MPI_BCAST (nfltvc,   1, imp_integers, 0, icomm, izmplcode)
    ENDIF
  ENDIF
#endif

!------------------------------------------------------------------------------
!- End of the Subroutine 
!------------------------------------------------------------------------------

END SUBROUTINE read_nc_gdefs

!==============================================================================
!+ Module procedure in "io_utilities" for writing a NetCDF file
!------------------------------------------------------------------------------

SUBROUTINE read_nc_vdefs (ncid, listin, nvarin, var_id, pollon, pollat, &
                           icomm, myid, npes, yerrmsg, istatus)
!
!------------------------------------------------------------------------------
!
! Description:
!   This routine reads attributes for each input variable 
!   from NetCDF formatted input file.
!
!------------------------------------------------------------------------------
!
!
! Subroutine / Function arguments

! Scalar arguments with intent(in):
  INTEGER (KIND=iintegers),   INTENT(IN) :: &
    ncid,           & ! NetCDF file ID
    nvarin,         & ! number of variables in input list
    icomm,          & ! MPI communicator
    myid,           & ! ID of this PE in communicator icomm
    npes              ! number of PEs

  REAL    (KIND=wp),        INTENT(IN)  ::  &
    pollat, pollon        ! coordinates of the rotated north pole

! Array arguments with intent(in):
  TYPE(list_description)  , INTENT(IN)  ::                      &
    listin(nvarin) ! List of fields for reading in

! Scalar arguments with intent(out):
  INTEGER (KIND=iintegers), INTENT(OUT)     ::  &
    istatus       ! error index

  CHARACTER (LEN= *),        INTENT(OUT)   ::  &
  yerrmsg     ! string for error messages

! Array arguments with intent(out):
  INTEGER (KIND=iintegers), INTENT(OUT) :: &
    var_id(nvarin)       ! NetCDF ID for each variable in the input list


!-----------------------------------------------------------------------
!
! Local scalars:

  INTEGER (KIND=iintegers) ::  &
    izmplcode           ! error status variable

  INTEGER (KIND=iintegers) ::  &
    j1                  ! loop index

  CHARACTER (LEN=11) ::  &
    coordinates         ! string holding names of the grid coordinates"
    
  CHARACTER (LEN=10)::  &
    ylocalname

  LOGICAL :: &
    lgrid_info_int2lm=.FALSE.

! Local arrays: 

!- End of header
!------------------------------------------------------------------------------

#ifdef NETCDF
! Processor 0 does the job
  IF (myid == 0) THEN
  
  var_loop: DO j1 = 1, nvarin

    ylocalname = listin(j1)%name

    ! get the variable ID
    istatus = nf90_inq_varid (ncid, TRIM(ylocalname), var_id(j1))
    IF (istatus /= NF90_NOERR) THEN
      ! this variable could not be found: but checking this is done
      ! later in src_input
      istatus = 0  ! otherwise the program stops if this is the last iteration
      var_id(j1) = -1
      CYCLE var_loop
    ENDIF

    ! get the grid mapping
    istatus = nf90_get_att (ncid, var_id(j1), 'grid_mapping', grid_mapping)
    IF (istatus /= NF90_NOERR) THEN
      yerrmsg = 'grid_mapping: '//TRIM(NF90_strerror(istatus))
      RETURN
    ENDIF

    IF (lgrid_info_int2lm) THEN
    ! check, if the grid is rotated lat/lon
    IF (grid_mapping(1:12) /= 'rotated_pole') THEN
      istatus = -1
      yerrmsg = 'No rotated grid: '//grid_mapping(1:26)
      RETURN
    ENDIF
    
    ! check, if U coordinate is staggered
    IF ( TRIM(ylocalname) == 'U') THEN
      istatus = nf90_get_att (ncid, var_id(j1), 'grid_info', coordinates)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'U coordinates: '//TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
    ENDIF   

    ! check, if V coordinate is staggered
    IF ( TRIM(ylocalname) == 'V') THEN
      istatus = nf90_get_att (ncid, var_id(j1), 'grid_info', coordinates)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'U coordinates: '//TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
    ENDIF   

    ELSE

    ! check, if the grid is rotated lat/lon
    IF (grid_mapping(1:12) /= 'rotated_pole') THEN
      istatus = -1
      yerrmsg = 'No rotated grid: '//grid_mapping(1:26)
      RETURN
    ENDIF
    
    IF ( TRIM(ylocalname) == 'U') THEN
      istatus = nf90_get_att (ncid, var_id(j1), 'coordinates', coordinates)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'U coordinates: '//TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
      IF (coordinates(1:11) /= 'slonu slatu') THEN
        IF (myid == 0) THEN
          PRINT *, 'coordinates = ',coordinates 
        ENDIF
        yerrmsg = ' ERROR *** Wrong coordinates for U-velocity ***'
        istatus = -1
        RETURN
      ENDIF
    ENDIF   

    ! check, if V coordinate is staggered
    IF ( TRIM(ylocalname) == 'V') THEN
      istatus = nf90_get_att (ncid, var_id(j1), 'coordinates', coordinates)
      IF (istatus /= NF90_NOERR) THEN
        yerrmsg = 'U coordinates: '//TRIM(NF90_strerror(istatus))
        RETURN
      ENDIF
      IF (coordinates(1:11) /= 'slonv slatv') THEN
        IF (myid == 0) THEN
          PRINT *, 'coordinates = ',coordinates 
        ENDIF
        yerrmsg = ' ERROR *** Wrong coordinates for V-velocity ***'
        istatus = -1
        RETURN
      ENDIF
    ENDIF   

    ENDIF
  
  ENDDO  var_loop
  
  ENDIF
    
! send information to the other processors
  IF (npes > 1) THEN
    CALL MPI_BCAST(var_id, nvarin, imp_integers, 0, icomm, izmplcode)
  ENDIF
#endif

!------------------------------------------------------------------------------
!- End of the Subroutine 
!------------------------------------------------------------------------------

END SUBROUTINE read_nc_vdefs

!==============================================================================
!==============================================================================
!+ Preprocessing for input file name generation
!------------------------------------------------------------------------------
  
SUBROUTINE create_file_name (nbdstep, ydata, yextension, yformat, yname, ydatchk2)
    
  ! Parameters
    
  INTEGER (KIND=iintegers), INTENT(IN), OPTIONAL  ::            &
    nbdstep        ! time step to create the file name

  CHARACTER (LEN=*)       , INTENT(IN)  ::                      &
    ydata          ! determines whether `initial' or `boundary' data

  CHARACTER (LEN=1),        INTENT(IN)     ::    &
    yextension     ! indicates actual restart-file

  CHARACTER (LEN= 4),       INTENT(IN)  ::    &
    yformat        ! determines the format that has to be read

  CHARACTER (LEN=250) ,     INTENT(OUT) ::  &
    yname      ! name of the grib file (is determined in make_fn)

  CHARACTER (LEN= 14) ,     INTENT(INOUT)      ::  &
    ydatchk2   ! for date-checking

  ! Local parameters

  CHARACTER (LEN= 14)        ::  &
    ydatchk    ! for date-checking

  INTEGER  (KIND=iintegers)  ::  &
    ierr ! status and error status variables

  INTEGER  (KIND=iintegers)  ::  &
    izdebug

  INTEGER  (KIND=iintegers)  ::  &
    nzday            ! century day of the current date

  REAL (KIND=wp)             ::  &
    zhour

  CHARACTER (LEN= 28)        ::  &
    ydum       ! for date-checking

  CHARACTER (LEN=  3)        ::  &
    yzhead     ! header of grib-file name

  LOGICAL                    ::  &
    lgexist

!------------------------------------------------------------------------------

  ! Initialize, whether additional debug output shall be done
  IF (ldebug_io) THEN
    IF (lprintdeb_all) THEN
      izdebug = idbg_level
    ELSE
      IF (my_cart_id == 0) THEN
        izdebug = idbg_level
      ELSE
        izdebug = 0
      ENDIF
    ENDIF
  ELSE
    izdebug = 0
  ENDIF

  IF (izdebug >= 10) PRINT *, '  CREATE file names'

  ! Create the file name: this routine has to be called with different
  ! arguments for initial or boundary data, resp.

  IF     (ydata == 'initial') THEN
    yzhead = 'laf'
    ! check for a file name where the date has 14 digits
    CALL make_fn (yzhead, ydate_ini, ydate_ini, ytunitbd, yextension, 0, dt, &
                .TRUE., itype_calendar, ydirini, yname, .TRUE., izdebug, ierr)

    ! In case of NetCDF IO, add ".nc" to the filename
    IF (yformat == 'ncdf') THEN
     yname = TRIM(yname)//'.nc'
    ENDIF

    lgexist = .FALSE.
    INQUIRE (FILE=yname, EXIST=lgexist)

    IF (.NOT. lgexist .AND. ydate_ini(11:14) == '0000') THEN
      ! try whether a 10 digit file exists
      CALL make_fn (yzhead, ydate_ini, ydate_ini, ytunitbd, yextension, 0, dt, &
                 .TRUE., itype_calendar, ydirini, yname, .FALSE., izdebug, ierr)

      ! In case of NetCDF IO, add ".nc" to the filename
      IF (yformat == 'ncdf') THEN
        yname = TRIM(yname)//'.nc'
      ENDIF

      ! if this file does not exist, SR open_file will fail
      ! therefore we do not need to stop here
      INQUIRE (FILE=yname, EXIST=lgexist)
   ENDIF

    ! from Version 4.24 on, ydatchk has 14 digits
    ydatchk = yakdat1

    ! Bug Fix by Burkhardt Rockel
    ! the model would crash if you check a NetCDF reference date with
    ! lbdana = .FALSE., ytunitbd != "d" and ydata = "boundary".
    ! (CLM Mantis BugId 25)
    ydatchk2(1:14) = ydatchk(1:14) ! from Version 4.24 on, ydatchk has 14 digits
    !ydatchk2(11:14) = '0000'

  ELSEIF (ydata == 'boundary') THEN
    IF (.NOT. lbdana) THEN
      yzhead = 'lbf'

      ! have to compute actual date for that case
      CALL get_utc_date (nbdstep+ndiff_ini_bd, ydate_ini, dt, itype_calendar,&
                           ydatchk, ydum, nzday, zhour )

      IF (ytunitbd == 'd') THEN

        ! Still bug fix by BR
        ydatchk2(1:14) = ydatchk(1:14)
        !! from 4.24 on, ydatchk has 14 digits: ydatchk2(11:14) = '0000'

        IF (izdebug > 10) THEN
          PRINT *, '     Creating file name with:  ', nbdstep, ndiff_ini_bd, &
                    dt, ' for actual date:  ', ydatchk
        ENDIF

        ! If the 10-digit format is used, filenames have to be determined to full
        ! hours, therefore set the last 4 digits of ydatchk, which might not
        ! be initialized
        IF (.NOT. lmmss) ydatchk(11:14) = '0000'
        CALL make_fn (yzhead, ydatchk, ydate_ini, ytunitbd, yextension,         &
                      nbdstep+ndiff_ini_bd, dt, .TRUE., itype_calendar, ydirbd, &
                      yname, .TRUE., izdebug, ierr)

        ! In case of NetCDF IO, add ".nc" to the filename
        IF (yformat == 'ncdf') THEN
         yname = TRIM(yname)//'.nc'
        ENDIF

        ! check whether the file with 14 digits in the date exists
        INQUIRE (FILE=yname, EXIST=lgexist)
        IF (.NOT. lgexist .AND. ydatchk(11:14) == '0000') THEN
          ! check whether the file with 10 digits in the date exists
          CALL make_fn (yzhead, ydatchk, ydate_ini, ytunitbd, yextension,       &
                        nbdstep+ndiff_ini_bd, dt, .TRUE., itype_calendar,       &
                        ydirbd, yname, .FALSE., izdebug, ierr)

          ! In case of NetCDF IO, add ".nc" to the filename
          IF (yformat == 'ncdf') THEN
           yname = TRIM(yname)//'.nc'
          ENDIF

          ! if this file does not exist, SR open_file will fail
          ! INQUIRE (FILE=yname, EXIST=lgexist)
        ENDIF

      ELSE

        ydatchk2(1:14)  = ydatchk(1:14)

        IF (izdebug > 10) THEN
          PRINT *, '     Creating file name with:  ', nbdstep, ndiff_ini_bd, dt
        ENDIF
        CALL make_fn (yzhead, ydate_bd, ydate_ini, ytunitbd, yextension,        &
                      nbdstep+ndiff_ini_bd, dt, .TRUE., itype_calendar,         &
                      ydirbd, yname, .FALSE., izdebug, ierr)

        ! In case of NetCDF IO, add ".nc" to the filename
        IF (yformat == 'ncdf') THEN
         yname = TRIM(yname)//'.nc'
        ENDIF
      ENDIF
    ELSE
      yzhead = 'laf'
      CALL get_utc_date (nbdstep+ndiff_ini_bd, ydate_ini, dt, itype_calendar,&
                         ydatchk, ydum, nzday, zhour )

      ! Still bug fix by BR
      ydatchk2(1:14) = ydatchk(1:14) ! from Version 4.24 on ydatchk has 14 digits
      !ydatchk2(11:14) = '0000'

      IF (izdebug > 10) THEN
        PRINT *, '     Creating file name with:  ', nbdstep, ndiff_ini_bd, dt
      ENDIF
      ! check if file with 14 digits in the date exists
      CALL make_fn (yzhead, ydatchk, ydate_ini, ytunitbd, yextension,           &
                    nbdstep+ndiff_ini_bd, dt, .TRUE., itype_calendar,           &
                    ydirbd, yname, .TRUE., izdebug, ierr)

      ! In case of NetCDF IO, add ".nc" to the filename
      IF (yformat == 'ncdf') THEN
        yname = TRIM(yname)//'.nc'
      ENDIF

      lgexist = .FALSE.
      INQUIRE (FILE=yname, EXIST=lgexist)

! this does not work with a dt, where the full hour is not reached, because
! then ydatchk(11:14) /= '0000'
!US   IF (.NOT. lgexist .AND. ydatchk(11:14) == '0000') THEN
      IF (.NOT. lgexist) THEN
         ! check whether the file with 10 digits in the date exists
         CALL make_fn (yzhead, ydatchk, ydate_ini, ytunitbd, yextension,        &
                       nbdstep+ndiff_ini_bd, dt, .TRUE., itype_calendar,        &
                       ydirbd, yname, .FALSE., izdebug, ierr)

         ! In case of NetCDF IO, add ".nc" to the filename
         IF (yformat == 'ncdf') THEN
           yname = TRIM(yname)//'.nc'
         ENDIF

         ! if this file does not exist, open_file will fail
         ! INQUIRE (FILE=yname, EXIST=lgexist)
      ENDIF
    ENDIF
  ELSEIF (ydata == 'restart') THEN
    yzhead = 'lrf'

    IF (ytunit_restart == 'd') THEN
      ! If the 10-digit format is used, filenames have to be determined to full
      ! hours, therefore set the last 4 digits of ydatchk, which might not
      ! be initialized
      IF (.NOT. lmmss) yakdat1(11:14) = '0000'
    ENDIF

    CALL make_fn (yzhead, yakdat1, ydate_ini, ytunit_restart, yextension, nstart, dt, &
                  .TRUE., itype_calendar, ydir_restart_in, yname, .TRUE., izdebug, ierr)

    IF (ytunit_restart == 'd') THEN
       ! check if file with 14 digits in the date exists
       lgexist = .FALSE.
       INQUIRE (FILE=yname, EXIST=lgexist)
       IF (.NOT. lgexist .AND. yakdat1(11:14) == '0000') THEN
         ! check whether the file with 10 digits in the date exists
         CALL make_fn (yzhead, yakdat1, ydate_ini, ytunit_restart, yextension, nstart, &
                       dt, .TRUE., itype_calendar, ydir_restart_in, yname, .FALSE., izdebug, ierr)
         ! if this file does not exist, open_file will fail
         ! INQUIRE (FILE=yname, EXIST=lgexist)
       ENDIF
    ENDIF
  ENDIF

END SUBROUTINE create_file_name

!==============================================================================
!==============================================================================

SUBROUTINE wait_for_file(yname, yname2)

  IMPLICIT NONE

  ! arguments
  CHARACTER (LEN=*), INTENT(IN) :: yname               ! name of the file
  CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: yname2    ! name of optional second filename

  ! local variables
  INTEGER :: nwait, ierr, i
  LOGICAL :: lexist, lyname2
  CHARACTER (LEN=200) :: yerrmsg
  CHARACTER (LEN= 25) :: yroutine

  yroutine = 'wait_for_file'
  ierr = 0_iintegers
  yerrmsg = ''
  lexist  = .FALSE.

  ! only check for file on PE 0
  IF (my_cart_id == 0) THEN

    ! check whether optional second filename should be handled
    lyname2 = PRESENT(yname2)
    IF (lyname2) lyname2 = LEN_TRIM(yname2) > 0

    nwait = 0
    DO WHILE ( (.NOT. lexist) .AND. (nwait < nmaxwait) )
      ! check for the first file
      INQUIRE (FILE=yname, EXIST=lexist)

      ! file not present, check for optional alternative file
      IF (.NOT. lexist .AND. lyname2) THEN
         INQUIRE (FILE=yname2, EXIST=lexist)
      ENDIF

      ! if not present, wait nincwait seconds and try again
      IF (.NOT. lexist) THEN
        IF (idbg_level > 1) THEN
          IF (lyname2) THEN
            PRINT *, 'file not available: ', TRIM(yname2)
          ELSE
            PRINT *, 'file not available: ', TRIM(yname)
          ENDIF
          PRINT *, '      sleep ', nincwait,' seconds'
        ENDIF
        i = dosleep(nincwait)
        nwait = nwait + nincwait
      ENDIF

    ENDDO

  ENDIF

  IF (num_compute > 1) THEN
    ! Synchronize the processes again
    CALL comm_barrier (icomm_cart, ierr, yerrmsg)

    ! Distribute lexist to all nodes
    CALL distribute_values (lexist, 1, 0, imp_logical, icomm_cart, ierr)
  ENDIF

  IF (.NOT. lexist) THEN
    ierr = 2006
    IF (lyname2) THEN
      yerrmsg  = ' *** ERROR:  ready-file 2 not available: '//TRIM(yname2)
    ELSE
      yerrmsg  = ' *** ERROR:  ready-file not available: '//TRIM(yname)
    ENDIF
    CALL model_abort (my_cart_id, ierr, yerrmsg, yroutine)
  ENDIF

END SUBROUTINE wait_for_file

!==============================================================================

END MODULE src_input
