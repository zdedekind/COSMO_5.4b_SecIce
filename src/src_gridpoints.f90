!+ Source module for gridpoint output
!------------------------------------------------------------------------------

MODULE src_gridpoints 

!------------------------------------------------------------------------------
!
! Description:
!   This module provides routines for the calculation of variable for    
!   gridpoint output
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  Ulrich.Schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Guenter Doms
!  Initial release
! 1.8        1998/08/03 Ulrich Schaettler
!  Eliminated dependency from data_io by introducing a parameterlist for 
!  print_gplong and print_gpshort.
! 1.10       1998/09/29 Ulrich Schaettler
!  Eliminate dependency on routine remark; changed formats for printing;
!  introduced extended station list.
! 1.12       1998/10/19 Ulrich Schaettler
!  Changed formats for long grid point output.
! 1.20       1999/01/07 Guenther Doms
!  Renaming of some global variables
! 1.29       1999/05/11 Ulrich Schaettler
!  Adapted interfaces to utility-modules and included the collection of
!  grid points
! 1.30       1999/06/24 Matthias Raschendorfer
!  New subroutine print_gpspec for writing the 2 new output files YUSTART and YUSPROF,
!  which are used for initialization and forcing the 1-d-model.
!  Use from module model_config: vcoord, vcflat, p0sl, t0sl, dt0lp.
!  Use from module data_constants: rvd_m_o.
!  Use from module data_fields: lai, h_can, d_pat, c_big, c_sml, r_air.
!  Use from module data_diagnostics: gp(hhl,lai,h_con,d_pat,rbig,rsml,dbig,dsml).
!  Chang of the additional parameters in CALL collect_gplong introduced by version 1.30
!                                    (strt,nust,mesd,numd)grpt.
!  Use from module data_runcontrol: lgpspec.
!  Use form module data_soil: cdzw(12,22,13,23,33).
!  Introduction of FUNCTIOn flin for linear interpolation.
!  Additional parameters in CALL collect_gplong:
!  (z)gp(hhl,lai,h_con,d_pat,rdrg,sai,cbig,csml), which are needed in print_gpspec.
!  (z)gpwu is used with onlyb one dimension.
! 1.32       1999/08/24 Guenther Doms
!  Some REAL variables changed to x_ireals.
! 1.33       1999/10/14 Matthias Raschendorfer
!  Change of the additional parameters in CALL collect_gplong introduced 
!  by version 1.30 to: (z)gp(hhl,lai,h_con,d_pat,rdrg,sai,cbig,csml).
! 1.34       1999/12/10 Ulrich Schaettler
!  Initialization of canopy or turbulence fields only if these are allocated.
! 1.39       2000/05/03 Ulrich Schaettler
!  Included declaration of variables from former module data_diagnostics
!  or module data_runcontrol, resp and changed some variable names.
! 2.8        2001/07/06 Ulrich Schaettler
!  Now use actual surface fluxes instead of summarized
! 3.5        2003/09/02 Ulrich Schaettler
!  Prepared output for cloud ice and new soil model and changed output format.
! 3.6        2003/12/11 Ulrich Schaettler
!  Initiated output for multi-layer soil model
! 3.7        2004/02/18 Ulrich Schaettler
!  Renamed alb (alb_rad), phi (rlat), rla (rlon)
! 3.13       2004/12/03 Reinhold Schrodin
!  Changed name of w_ice to w_so_ice
! 3.17       2005/12/12 Reinhold Schrodin
!  Additional variables: rho_snow, h_snow
! 3.18       2006/03/03 Ulrich Schaettler, Matthias Raschendorfer
!  Splitting the output: every grid point is written to an extra file.
!  Eliminated print_gpshort and collect_gpshort (not necessary any more).
!  New output files 'MODLEV' and 'MESDAT' for the single column model.
! 3.21       2006/12/04 Ulrich Schaettler, Matthias Raschendorfer
!  Removed blanks from the station names; adapted search for stations
!  Introduced timelevel in the argument list of subroutine gridpoints
!  Splitting the special SCLM-output as well,
!  thus eliminating 'print_gplong' and collect_gplong.
!  Information of files 'MODLEV' and 'MESDAT' are now in only one file
!  for each grid point with the same nomenclature as for long and short output.
! V3_23        2007/03/30 Ulrich Schaettler
!  Treatment of gpsoil also for the short meteographs
!  SCLM-output now for height based model levels (M. Raschendorfer)
! V4_3         2008/02/25 Ulrich Schaettler
!  Eliminated 'LM' from ASCII output
! V4_4         2008/07/16 Ulrich Schaettler
!  Changed NL parameter lyear_360 to itype_calendar, to have several options
! V4_8         2009/02/16 Hans-Juergen Panitz
!  Adapted an output format for longer numbers in climate mode
! V4_10        2009/09/11 Davide Cesari
!  Add characters after a backslash in comments; g95 interprets that as
!  continuation line
! V4_11        2009/11/30 Ekaterina Machulskaya
!  Adaptations for multi-layer snow model
! V4_12        2010/05/11 Ulrich Schaettler
!  Renamed t0 to t0_melt because of conflicting names
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_18        2011/05/26 Ulrich Schaettler
!  Changed the code owner
! V4_20        2011/08/31 Ulrich Schaettler
!  Bug fix for computing total pressure on highest half level (J. Schmidli)
! V4_23        2012/05/10 Ulrich Schaettler (Thorsten Reinhardt, CLM)
!  Additional variables in long  meteographs: SWDIR_S, SWDIFD_S, SWDIFU_S
!                                             SNOW_MELT
!  Changed some output formats to avoid "stars" during climate runs
! V4_24        2012/06/22 Hendrik Reich
!  Adapted the formats of the date variables
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer
!  Replaced qx-variables by using them from the tracer module
!  Check soiltyp with NINT-function
!  Adapted unit for height of convective clouds in the meteograph output
!  Increased length of station names to 30 characters
! V4_26        2012/12/06 Denis Blinov, Gdaly Rivin
!  Corrected dimensions for wind speed in meteographs output, which is
!    computed now in m/s and not in knots
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  Introduced MESSy interface
! V4_28        2013/07/12 Ulrich Schaettler
!  Removed un-used variables
!  Use subroutines and variables for vertical grid and reference atmospheres
!    from module vgrid_refatm_utils
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
!  For the COSMO-Model only use vcoord from vgrid_refatm_utils
! V5_1         2014-11-28 Ulrich Schaettler, Oliver Fuhrer
!  Replaced vcoord(k) > vcflat by k > kflat to remove vertical coordinate parameters
!  Replaced ireals by wp (working precision) (OF)
! V5_2         2015-05-21 Heike Vogel
!  Added allocation of gp-fields for COSMO-ART
! V5_4a        2016-05-10 Matthias Raschendorfer
!  Changing the variable name for liquid water content from 'QW_Air' into 'QL_Air'
!   for the special 'mesdat'-output (to be used by a later single-column run).
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Declarations:
!
! Modules used:

USE data_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 1. vertical coordinate parameters and related variables
! --------------------------------------------------------------------
    czmls,        &! depth of the main soil layers in meters

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------

    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction
    ke1,          & ! KE+1
    ke_soil,      & ! number of layers in the multi-layer soil model
    ie_tot,       & ! number of grid points in zonal direction
    je_tot,       & ! number of grid points in meridional direction

! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the 
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from 
!    the other ones because of the use of the staggered Arakawa-B-grid.
!    
!   zonal direction
    istart,       & ! start index for the forecast of w, t, qv, qc and pp
    iend,         & ! end index for the forecast of w, t, qv, qc and pp
    istartpar,    & ! start index for computations in the parallel program
    iendpar,      & ! end index for computations in the parallel program

!   meridional direction
    jstart,       & ! start index for the forecast of w, t, qv, qc and pp
    jend,         & ! end index for the forecast of w, t, qv, qc and pp
    jstartpar,    & ! start index for computations in the parallel program
    jendpar,      & ! end index for computations in the parallel program

! 4. constants for the horizontal rotated grid and related variables
! ------------------------------------------------------------------

    dlon,         & ! grid point distance in zonal direction (in degrees)
    dlat,         & ! grid point distance in meridional direction (in degrees)
    pollon,       & ! longitude of the rotated north pole (in degrees, E>0)
    pollat,       & ! latitude of the rotated north pole (in degrees, N>0)
    raddeg,       & ! factor for transforming rad to degree

! 5. variables for the time discretization and related variables
! --------------------------------------------------------------

    dt,           & ! long time-step
    dtdeh,        & ! dt / 3600 seconds

! 8. Organizational variables to handle the COSMO humidity tracers
! ----------------------------------------------------------------
      idt_qc, idt_qv, idt_qi

!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 2. physical constants and related variables
! -------------------------------------------

    t0_melt,      & ! melting temperature of ice   
    r_d,          & ! gas constant for dry air
    rdv,          & ! r_d / r_v
    o_m_rdv,      & ! 1 - r_d/r_v
    rvd_m_o,      & ! r_v/r_d - 1
    g,            & ! acceleration due to gravity
    day_len,      & ! mean length of the day

! 3. constants for parametrizations
! ---------------------------------

    b1,           & ! variables for computing the saturation vapour pressure
    b2w,          & ! over water (w) and ice (i)
    b3,           & !               -- " --
    b4w             !               -- " --

! end of data_constants

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------

    p0 ,            & ! reference pressure at main levels             ( pa  )
    hhl,            & ! height of model half levels                   ( m )

! 2. external parameter fields                                        (unit)
! ----------------------------

    hsurf,          & ! height of surface topography                  (m)
    fr_land,        & ! fraction of land in a grid element              -- 
    soiltyp,        & ! type of the soil (keys 0-9)                     --
    vio3,           & ! vertical integrated ozone contents            (pa O3)
    hmo3,           & ! ozone maximum                                 ( pa  )
    fc,             & ! coriolis-parameter                            ( 1/s )
    plcov,          & ! fraction of plant cover                         --
    lai,            & ! plant area index                                --
    sai,            & ! surface area index                              --
    rootdp,         & ! depth of the roots                            ( m  )
    gz0,            & ! roughness length * g                          (m2/s2)
    rlat,           & ! geographical latitude                         ( rad )
    rlon,           & ! geographical longitude                        ( rad )
    h_can,          & ! hight of the vertically resolved canopy       ( m )
    d_pat,          & ! horizontal pattern length scale               ( m )
    c_big,          & ! effective drag coefficient of canopy elements
                      ! larger than or equal to the tubulent length
                      ! scale                                         (1/m)
    c_sml,          & ! effective drag coefficient of canopy elements
                      ! smaller than the tubulent length scale        (1/m)
    r_air             ! log of air containing fraction of a gridbox
                      ! inside the canopy                             ( 1 )

USE data_fields     , ONLY :   &

! 3. prognostic variables                                             (unit)
! -----------------------

    u,              & ! zonal wind speed                              ( m/s )
    v,              & ! meridional wind speed                         ( m/s )
    w,              & ! vertical wind speed (defined on half levels)  ( m/s )
    t,              & ! temperature                                   (  k  )
    pp,             & ! deviation from the reference pressure         ( pa  )

! 5. fields for surface values and soil model variables               (unit )
! -----------------------------------------------------

    ps ,            & ! surface pressure                              ( pa  )
    t_snow,         & ! temperature of the snow-surface               (  k  )
    t_snow_mult,    & ! temperature of the snow-surface               (  k  )
    t_s,            & ! temperature of the ground surface             (  k  )
    t_g,            & ! weighted surface temperature                  (  k  )
    qv_s,           & ! specific water vapor content on the surface   (kg/kg)
    t_m,            & ! temperature between upper and medium 
                      ! soil layer                                    (  k  )
    t_cl,           & ! temperature between medium and lower
                      ! soil layer                                    (  k  )
    w_snow,         & ! water content of snow                         (m H2O)
    w_i,            & ! water content of interception water           (m H2O)
    w_g1,           & ! water content of the upper soil layer         (m H2O)
    w_g2,           & ! water content of the medium soil layer        (m H2O)
    w_g3,           & ! water content of the lower soil layer         (m H2O)
    w_cl,           & ! climatological water content                  (m H2O)
    t_so,           & ! multi-layer soil temperature                  (  k  )
    w_so,           & ! multi-layer soil moisture                     (m H2O)
    w_so_ice,       & ! multi-layer soil ice                          (m H2O)
    freshsnow,      & ! weighting function indicating 'freshness' of snow
    rho_snow ,      & ! prognostic snow density                       (kg/m3)
    rho_snow_mult,  & ! prognostic snow density                       (kg/m3)
    h_snow            ! snow height                                   (  m  )

USE data_fields     , ONLY :   &

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------

!   coefficients for turbulent diffusion in the atmosphere
!   (defined on half levels)
    tkvm,           & ! turbulent diffusion coefficient for momentum  (m2/s) 
    tkvh,           & ! turbulent diffusion coefficient for heat      (m2/s) 
                      ! and moisture
    tcm      ,      & ! transfer coefficient for momentum             ( -- ) 
    tch      ,      & ! transfer coefficient for heat and moisture    ( -- ) 

!   fields from the radiation scheme
    clc_sgs,        & ! subgrid-scale stratiform cloud cover            -- 
    alb_rad,        & ! albedo of the ground                            --
    swdir_s  ,      & ! direct comp. of solar radiative flux at surface ( W/m2)
    swdifd_s ,      & ! diffuse downward comp. of short wave rad. flux  ( W/m2)
    swdifu_s ,      & ! diffuse upward   comp. of short wave rad. flux  ( W/m2)
    sobs ,          & ! solar radiation at the ground                 ( w/m2)
    thbs ,          & ! thermal radiation at the ground               ( w/m2)
    pabs ,          & ! photosynthetic active radiation at the ground ( w/m2)
    sobt ,          & ! solar radiation at the upper boundary         ( w/m2)
                      ! of the atmosphere
    thbt ,          & ! thermal radiation at the upper boundary       ( w/m2)
                      ! of the atmosphere
    clch ,          & ! cloud cover with high clouds                    --   
    clcm ,          & ! cloud cover with medium clouds                  --   
    clcl ,          & ! cloud cover with low clouds                     --   
    clct ,          & ! total cloud cover                               --   

!   fields from the convection scheme 
    clc_con    ,    & ! cloud cover due to convection                   --
    prr_con    ,    & ! precipitation rate of rain, convective        (kg/m2*s)
    prs_con    ,    & ! precipitation rate of snow, convective        (kg/m2*s)
    bas_con    ,    & ! level index of convective cloud base            --
    top_con    ,    & ! level index of convective cloud top             --

!   fields from the grid-scale precipitation scheme 
    prr_gsp    ,    & ! precipitation rate of rain, grid-scale        (kg/m2*s)
    prs_gsp    ,    & ! precipitation rate of snow, grid-scale        (kg/m2*s

!   fields that are computed in the dynamics
    dpsdt      ,    & ! tendency of the surface pressure              ( pa/s)
    umfl_s     ,    & ! u-momentum flux (surface)                     ( N/m2)
    vmfl_s     ,    & ! v-momentum flux (surface)                     ( N/m2)
    shfl_s     ,    & ! average sensible heat flux (surface)          ( W/m2)
    lhfl_s            ! average latent heat flux (surface)            ( W/m2)

USE data_fields     , ONLY :   &

! 7. fields for model output and diagnostics                          (unit )
! ------------------------------------------

    t_2m       ,    & ! temperature in 2m                             (  K  )
    qv_2m      ,    & ! specific water vapor content in 2m            (kg/kg)
    td_2m      ,    & ! dew-point in 2m                               (  K  )
    u_10m      ,    & ! zonal wind in 10m                             ( m/s )
    v_10m      ,    & ! meridional wind in 10m                        ( m/s )
    tmin_2m    ,    & ! minimum temperature in 2m                     (  K  )
    tmax_2m    ,    & ! maximum temperature in 2m                     (  K  )
    vmax_10m   ,    & ! maximal windspeed in 10m                      ( m/s )
    rain_gsp   ,    & ! amount of rain from grid-scale precip. (sum)  (kg/m2)
    snow_gsp   ,    & ! amount of snow from grid-scale precip. (sum)  (kg/m2)
    rain_con   ,    & ! amount of rain from convective precip. (sum)  (kg/m2)
    snow_con   ,    & ! amount of snow from convective precip. (sum)  (kg/m2)
    runoff_s   ,    & ! surface water drainage; sum over forecast     (kg/m2)
    runoff_g   ,    & ! soil water drainage; sum over forecast        (kg/m2)
    snow_melt         ! amount of snow melt; sum over forecast        (kg/m2)

! end of data_fields
!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    nstart,       & ! first time step of the forecast
    ntstep,       & ! actual time step

! 3. controlling the physics
! --------------------------
    nlgw,         & ! number of prognostic soil water levels
    lmulti_layer, & ! run multi-layer soil model
    lmulti_snow , & ! run multi-layer snow model
    ltur,         & ! forecast with vertical diffusion
    lphys,        & ! forecast with physical parametrizations
    lgsp,         & ! forecast with grid scale precipitation
    l_cosmo_art,  & ! if .TRUE., run the COSMO_ART

! 9. Variables for Ascii file handling, time measuring, ...
! ---------------------------------------------------------
    itype_calendar  ! for specifying the calendar used

! end of data_runcontrol 

!------------------------------------------------------------------------------

USE data_parallel,      ONLY :  &
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    my_cart_pos,     & ! position of this subdomain in the cartesian grid
                       ! in x- and y-direction
    my_cart_neigh,   & ! neighbors of this subdomain in the cartesian grid
    isubpos,         & ! positions of the subdomains in the total domain. Given
                       ! are the i- and the j-indices of the lower left and the
                       ! upper right grid point in the order
                       !                  i_ll, j_ll, i_ur, j_ur.
                       ! Only the interior of the domains are considered, not
                       ! the boundary lines.
    icomm_cart,      & ! communicator for the virtual cartesian topology
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    nboundlines,     & ! number of boundary lines
    num_compute        ! number of compute PEs

!------------------------------------------------------------------------------

USE data_soil,          ONLY :  &
     cdzw12,    & !  thickness of upper soil water layer in two-layer model
     cdzw22,    & !  thickness of lower soil water layer in two-layer model
     cdzw13,    & !  thickness of upper soil water layer in three-layer model
     cdzw23,    & !  thickness of middle soil water layer in three-layer model
     cdzw33       !  thickness of lower soil water layer in three-layer model

!------------------------------------------------------------------------------

USE utilities,          ONLY :  &
      get_utc_date,        & ! actual date of the forecast in different forms
      uvrot2uv               !

!------------------------------------------------------------------------------

USE environment,        ONLY :  model_abort, get_free_unit

!------------------------------------------------------------------------------

USE vgrid_refatm_utils, ONLY :  vcoord

!------------------------------------------------------------------------------

USE src_tracer,         ONLY:   &
      trcr_get,                 &
#ifdef COSMOART
      trcr_get_block,           &
#endif
      trcr_errorstr

USE data_tracer,        ONLY:   T_ERR_NOTFOUND

!------------------------------------------------------------------------------

#ifdef COSMOART
USE data_cosmo_art,     ONLY:   lvolcano, ldust, lrad_dust, lonly_dust,       &
                                tau_dust_550nm, isp_aerotrans, trcr_idx_aero
USE art_species_data,   ONLY:   vsoila, vsoilc, vsoila0, vsoilc0
USE data_volcano,       ONLY:   isp_ash, ash_scale, trcr_idx_ash
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Variables for the module src_gridpoints

! `gp' stands for gridpoint; the other letters indicate the concerned variable.
! The first index is for referring to the selected gridpoints, the second
! index refers to the time-level.

! 1a) Arrays for the meteograms  (short gridpoint output)
!--------------------------------------------------------

  REAL   (KIND=wp),     ALLOCATABLE     ::             &
    gpps  (:,:) , & ! surface pressure
    gpu10m(:,:) , & ! zonal wind in 10m height
    gpv10m(:,:) , & ! meridional wind in 10m height
    gpu950(:,:) , & ! zonal wind in 950 hPa
    gpv950(:,:) , & ! meridional wind in 950 hPa
    gpu850(:,:) , & ! zonal wind in 850 hPa
    gpv850(:,:) , & ! meridional wind in 850 hPa
    gpu700(:,:) , & ! zonal wind in 700 hPa
    gpv700(:,:) , & ! meridional wind in 700 hPa
    gpu500(:,:) , & ! zonal wind in 500 hPa
    gpv500(:,:) , & ! meridional wind in 500 hPa
    gpt_g (:,:) , & ! ground temperature
    gpt2m (:,:) , & ! temperature in 2m height
    gptp  (:,:) , & ! temperature on lowest model level
    gpt850(:,:) , & ! temperature in 850 hPa
    gpt700(:,:) , & ! temperature in 700 hPa
    gpt500(:,:) , & ! temperature in 500 hPa
    gptd2m(:,:) , & ! dew-point in 2m height
    gpclcl(:,:) , & ! cloud cover     low
    gpclcm(:,:) , & ! cloud cover     medium
    gpclch(:,:) , & ! cloud cover     high
    gpclct(:,:) , & ! cloud cover     total
    gpfog (:,:) , & ! cloud cover of lowest model layer
    gphbas(:,:) , & ! base height of main convective cloud
    gphtop(:,:) , & ! top  height of main convective cloud
    gprrsn(:,:) , & ! rain; scale precipitation;       summation since start
    gprssn(:,:) , & ! snow; convective precipitation;  summation since start
    gprrkn(:,:) , & ! rain; scale precipitation;       summation since start
    gprskn(:,:) , & ! snow; convective precipitation;  summation since start
    gpwsno(:,:) , & !
    gpsmlt(:,:)     ! snow melt;                       summation since start

! 1b) Arrays for long gridpoint output
!-------------------------------------

! The third dimension, if present, refers to the vertical level.

  REAL   (KIND=wp),     ALLOCATABLE     ::             &
    gpu   (:,:,:),  & ! prognostic variable for the zonal wind
    gpv   (:,:,:),  & ! prognostic variable for the meridional wind
    gpw   (:,:,:),  & ! prognostic variable for the vertical velocity
    gpt   (:,:,:),  & ! prognostic variable for the temperature
    gpqv  (:,:,:),  & ! prognostic variable for the water vapor
    gpqc  (:,:,:),  & ! prognostic variable for the cloud water
    gpqi  (:,:,:),  & ! prognostic variable for the cloud ice
    gppp  (:,:,:),  & ! prognostic variable for the deviation of the
                      ! reference pressure
    gphhl (:,:)  ,  & ! height of half levels
    gphml (:,:)  ,  & ! height of main levels
    gptkvm(:,:,:),  & ! turbulent diffusion coefficient for momentum
    gptkvh(:,:,:),  & ! turbulent diffusion coefficient for heat and moisture
    gpclcs(:,:,:),  & ! cloud cover (subgrid scale stratiform)
    gpclcc(:,:,:),  & ! cloud cover (convective)

    gptsno   (:,:), & ! snow temperature
    gpt_s    (:,:), & ! surface temperature
    gpt_m    (:,:), & ! temperature of medium soil layer (old soil model)
    gpqvs    (:,:), & ! specific water vapor content at the surface
    gpw_i    (:,:), & ! water content of interception water
    gpwg1    (:,:), & ! water content of the upper  soil layer   \   old
    gpwg2    (:,:), & ! water content of the medium soil layer    >  soil
    gpwg3    (:,:), & ! water content of the lower  soil layer   /   model
    gptso  (:,:,:), & ! temperatures for multi-layer soil model
    gpwso  (:,:,:), & ! water contents for multi-layer soil model
    gpwice (:,:,:), & ! soil ice for multi-layer soil model
    gpfrsn   (:,:), & ! weighting function indicating 'freshness' of snow
    gprhsn   (:,:), & ! prognostic snow desity
    gph_sn   (:,:), & ! snow height (m)
    gptcm    (:,:), & ! transfer coefficient for momentum 
    gptch    (:,:), & ! transfer coefficient for heat and moisture
    gpgz0    (:,:), & ! surface roughness * g 
    gpdpdt   (:,:), & ! tendency of the surface pressure
    gpswdir_s(:,:), & ! direct comp. of solar radiative flux at surface ( W/m2)
    gpswdifd_s(:,:),& ! diffuse downward comp. of short wave rad. flux  ( W/m2)
    gpswdifu_s(:,:),& ! diffuse upward   comp. of short wave rad. flux  ( W/m2)
    gpsobs   (:,:), & ! solar radiation at the ground
    gpthbs   (:,:), & ! thermal radiation at the ground
    gppabs   (:,:), & ! photosynthetic active radiation at the ground
    gpsobt   (:,:)    ! solar radiation at upper boundary of the atmosphere

  REAL   (KIND=wp),     ALLOCATABLE     ::             &
    gpthbt   (:,:), & ! thermal radiation at upper boundary of the atmosphere
    gpalb    (:,:), & ! albedo of the ground
    gpprrs   (:,:), & ! precipitation rate of rain, grid-scale
    gpprss   (:,:), & ! precipitation rate of snow, grid-scale
    gpprrc   (:,:), & ! precipitation rate of snow, convective
    gpprsc   (:,:), & ! precipitation rate of rain, convective
    gptmn2   (:,:), & ! minimum temperature in 2m
    gptmx2   (:,:), & ! maximum temperature in 2m
    gpvb10   (:,:), & ! maximal windspeed in 10m
    gpabsf   (:,:), & ! surface water runoff; sum over forecast
    gpabgr   (:,:), & ! soil water runoff; sum over forecast
    gpshfl   (:,:), & ! sensible heat flux (surface)
    gplhfl   (:,:), & ! latent heat flux (surface)
    gpumfl   (:,:), & ! u-momentum flux (surface)
    gpvmfl   (:,:), & ! v-momentum flux (surface)
    gpp0     (:,:), & ! reference pressure
    gpplco   (:,:), & ! plant cover
    gplai    (:,:), & ! leaf area index
    gproot   (:,:), & ! root depth
    gpt_cl   (:,:), & ! temperature medium/lower soil layer (climatology) 
    gpw_cl   (:,:), & ! climatological water content   (old soil model)
    gpvio3   (:,:), & ! vertical integrated ozone contents
    gphmo3   (:,:), & ! ozone maximum
    gphsur   (:)  , & ! height of surface
    gpfrla   (:)  , & ! fraction of land
    gplati   (:)  , & ! latitude
    gplong   (:)  , & ! longitude
    gpsoil   (:)  , & ! soil-type
    gpfc     (:)  , & ! coriolis parameter
    gphcan   (:)  , & ! height of the vertically resolved canopy
    gpdpat   (:)  , & ! horizontal pattern length scale
    gprdrg   (:)  , & ! air containing fraction of a gridbox inside the canopy
    gpsai    (:)  , & ! surface area index
    gpcbig   (:)  , & ! effective drag coefficient of canopy elements
                      ! larger than or equal to the tubulent length
                      ! scale
    gpcsml   (:)      ! effective drag coefficient of canopy elements
                      ! smaller than the tubulent length scale

#ifdef COSMOART
  REAL   (KIND=wp), ALLOCATABLE     ::             &
    gpcash(:,:,:,:),   &
    gpashmc(:,:,:),    &
    gpcaero(:,:,:,:)
#endif

! Organizational variables
  INTEGER (KIND=iintegers), PARAMETER ::                        &
    nmaxgp=100         ! maximal number of grid points for grid point output

  TYPE list_of_stations
    INTEGER (KIND=iintegers)  :: igp, jgp
    REAL (KIND=wp)            :: rlatgp, rlongp
    CHARACTER (LEN=30)        :: ystation_name
  END TYPE list_of_stations

  TYPE (list_of_stations)          ::           &
    stationlist_all(nmaxgp),  & ! list of station specifications for the
                                ! total domain (only for checking purposes)
    stationlist_tot(nmaxgp),  & ! list of station specifications for the
                                ! total domain
    stationlist    (nmaxgp)     ! list of station specifications for the
                                ! subdomain

  INTEGER   (KIND=iintegers)       ::           &
    ntot_in(nmaxgp)    ! index of a special grid point in the station list

  INTEGER   (KIND=iintegers)       ::           &
    n0gp,         & ! time step of the first grid point calculation
    nincgp,       & ! time step increment of grid point calculations
    nnextgp,      & ! next step for grid point output
    ngp_tot,      & ! number of chosen grid points in the total domain
    ngp,          & ! number of chosen grid points in this subdomain
    nstepsgp        ! number of time steps in which grid point calculations
                    ! are performed

  REAL    (KIND=wp)                :: &
    h0gp,         & ! first step of grid point calculation in hours
    hincgp,       & ! hour increment of grid point calculations
    hnextgp         ! next step for grid point output in hours

  LOGICAL                          ::           &
    lgpshort,     & ! calculate and print a short form of the
                    ! grid point calculations   (1 line/step)
    lgplong,      & ! calculate and print a long form of the
                    ! grid point calculations   (1 page/step)
    lgpspec         ! calculate and print a special form of the
                    ! grid point output used by the 1-d-model

  INTEGER (KIND=iintegers)     ::    &
    nustgrpt,   & ! grid point output for 1-d-MODLEV-file
    numdgrpt      ! grid point output for 1-d-MESDAT-file

  CHARACTER (LEN= 8)  ::   &
    strtgrpt='YUMODLEV',   & ! grid point output for 1-d-MODLEV-file
    mesdgrpt='YUMESDAT'      ! grid point output for 1-d-MESDAT-file

  ! organization of grid point output, if every grid point is written
  ! to an extra file by each processor
  CHARACTER (LEN=32)  ::   &
    yugpname(nmaxgp)     ! file name for meteographs

  INTEGER (KIND=iintegers)     ::    &
    nugpunit(nmaxgp)     ! unit numbers for these files

  ! organization of stations
  INTEGER (KIND=iintegers) , PARAMETER    , PRIVATE :: &
    nstatgp = 131        ! number of listed sounding stations

  REAL    (KIND=wp)                       , PRIVATE :: &
    latstat (nstatgp) ,& ! longitude  \  of the listed
    lonstat (nstatgp)    ! latitude   /  sounding stations

  CHARACTER   (LEN=30)                    , PRIVATE :: &
    ystatgp (nstatgp)    ! name of the listed sounding stations

#ifdef COSMOART
  ! definition of number of ash diameter classes and output format
  CHARACTER(LEN = 80)  :: header_format
  CHARACTER(LEN = 150) :: header_line
  CHARACTER(LEN = 150) :: header_unit
#endif

!------------------------------------------------------------------------------

INCLUDE "mpif.h"

!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure in "diagnostics" for initializing variables
!------------------------------------------------------------------------------

SUBROUTINE init_gridpoints (ydate_ini)

!------------------------------------------------------------------------------
!
! Description:
!   This routine performs the initializations for the gridpoint computations.
!   Tasks are opening of files, allocation of space, computing initial values
!   and initializing the list of sounding stations.
!
! Method:
!   
!------------------------------------------------------------------------------

! Parameterlist
CHARACTER (LEN=14), INTENT(IN)     ::   &
     ydate_ini      ! start of the forecast  yyyymmddhh (year, month, day, hour)

! Local variables
REAL (KIND=wp)             :: hh1, zchkgp

INTEGER (KIND=iintegers)   ::       &
  i, k, l, n, izgp, jzgp,           & ! support variables
  nstat, istat, ilo, iup, jlo, jup, njata1, irlat, irlon, izlen

CHARACTER (LEN=14)  ::   ydat11
CHARACTER (LEN=28)  ::   ydat12
CHARACTER (LEN=30)  ::   ystation
CHARACTER (LEN=10)  ::   ysoiltyp(10)  ! soil type names

CHARACTER (LEN=20)  ::   yroutine
CHARACTER (LEN=80)  ::   yerrmsg 

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE init_gridpoints
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  Section 1: Initializations
!------------------------------------------------------------------------------

  yroutine = 'init_gridpoints'

  ysoiltyp(1:10) = (/'ice       ', 'rock      ', 'sand      ', 'sandy loam', &
                     'loam      ', 'clay loam ', 'clay      ', 'peat      ', &
                     'sea water ', 'sea ice   '/)

!------------------------------------------------------------------------------
!  Section 2: Decomposition of the grid points
!------------------------------------------------------------------------------

  ! Grid points can be chosen on the boundary lines of the total domain.
  ! Therefore determine ilo, iup, jlo, jup to look in which subdomain a
  ! certain grid point is located.

  IF (num_compute > 1) THEN

    IF (my_cart_neigh(1) == -1) THEN
      ilo = 1
    ELSE
      ilo = isubpos(my_cart_id,1)
    ENDIF

    IF (my_cart_neigh(2) == -1) THEN
      jup = je_tot
    ELSE
      jup = isubpos(my_cart_id,4)
    ENDIF

    IF (my_cart_neigh(3) == -1) THEN
      iup = ie_tot
    ELSE
      iup = isubpos(my_cart_id,3)
    ENDIF

    IF (my_cart_neigh(4) == -1) THEN
      jlo = 1
    ELSE
      jlo = isubpos(my_cart_id,2)
    ENDIF
 
    ngp = 0

    DO l=1,ngp_tot
      IF ( (ilo <= stationlist_tot(l)%igp) .AND. (stationlist_tot(l)%igp <= iup) ) THEN
        IF ( (jlo <= stationlist_tot(l)%jgp) .AND. (stationlist_tot(l)%jgp <= jup) ) THEN

          ! grid point l is located in this subdomain (my_cart_id)
          ngp = ngp + 1
          stationlist(ngp)%igp = nboundlines +                             &
                        (stationlist_tot(l)%igp - isubpos(my_cart_id,1) + 1)
          stationlist(ngp)%jgp = nboundlines +                             &
                        (stationlist_tot(l)%jgp - isubpos(my_cart_id,2) + 1)
          stationlist(ngp)%ystation_name = stationlist_tot(l)%ystation_name
          ntot_in(ngp) = l
        ENDIF
      ENDIF
    ENDDO

  ELSE    ! num_compute = 1

    ! set the *_tot variables to the local variables
    ngp                = ngp_tot
    stationlist(:)%igp = stationlist_tot(:)%igp
    stationlist(:)%jgp = stationlist_tot(:)%jgp
    stationlist(:)%ystation_name = stationlist_tot(:)%ystation_name
    DO l = 1, ngp_tot
      ntot_in(l)   = l
    ENDDO

  ENDIF

!------------------------------------------------------------------------------
!  Section 2: Initializations of sounding stations
!------------------------------------------------------------------------------

  ystatgp( :) = '                    '

  ! Northern Germany
  ystatgp( 1) = 'SCHLESWIG         '; latstat( 1) = 54.53_wp; lonstat( 1) =  9.55_wp
  ystatgp( 2) = 'GREIFSWALD        '; latstat( 2) = 54.10_wp; lonstat( 2) = 13.40_wp
  ystatgp( 3) = 'EMDEN             '; latstat( 3) = 53.38_wp; lonstat( 3) =  7.23_wp
  ystatgp( 4) = 'HAMBURG           '; latstat( 4) = 53.55_wp; lonstat( 4) =  9.97_wp
  ystatgp( 5) = 'WITTSTOCK         '; latstat( 5) = 53.20_wp; lonstat( 5) = 12.52_wp
  ystatgp( 6) = 'MEPPEN            '; latstat( 6) = 52.72_wp; lonstat( 6) =  7.32_wp
  ystatgp( 7) = 'BERGEN            '; latstat( 7) = 52.82_wp; lonstat( 7) =  9.93_wp
  ystatgp( 8) = 'HANNOVER          '; latstat( 8) = 52.47_wp; lonstat( 8) =  9.70_wp
  ystatgp( 9) = 'MAGDEBURG         '; latstat( 9) = 52.12_wp; lonstat( 9) = 11.58_wp
  ystatgp(10) = 'LINDENBERG        '; latstat(10) = 52.22_wp; lonstat(10) = 14.12_wp
  ystatgp(11) = 'MOENCHENGLADBACH  '; latstat(11) = 51.23_wp; lonstat(11) =  6.50_wp
  ystatgp(12) = 'ESSEN             '; latstat(12) = 51.40_wp; lonstat(12) =  6.97_wp
  ystatgp(13) = 'KASSEL            '; latstat(13) = 51.13_wp; lonstat(13) =  9.27_wp
  ystatgp(14) = 'LEIPZIG           '; latstat(14) = 51.32_wp; lonstat(14) = 12.42_wp
  ystatgp(15) = 'WAHNSDORF         '; latstat(15) = 51.12_wp; lonstat(15) = 13.68_wp
  ystatgp(16) = 'GOERLITZ          '; latstat(16) = 51.17_wp; lonstat(16) = 14.95_wp
  ystatgp(17) = 'MEININGEN         '; latstat(17) = 50.57_wp; lonstat(17) = 10.38_wp
  ystatgp(18) = 'ERFURT            '; latstat(18) = 50.98_wp; lonstat(18) = 10.97_wp
  ystatgp(19) = 'CHEMNITZ          '; latstat(19) = 50.80_wp; lonstat(19) = 12.87_wp
 
  ! Southern Germany
  ystatgp(20) = 'FRANKFURT         '; latstat(20) = 50.05_wp; lonstat(20) =  8.60_wp
  ystatgp(21) = 'OFFENBACH         '; latstat(21) = 50.12_wp; lonstat(21) =  8.73_wp
  ystatgp(22) = 'HOF               '; latstat(22) = 50.32_wp; lonstat(22) = 11.88_wp
  ystatgp(23) = 'IDAR-OBERSTEIN    '; latstat(23) = 49.70_wp; lonstat(23) =  7.33_wp
  ystatgp(24) = 'SAARBRUECKEN      '; latstat(24) = 49.22_wp; lonstat(24) =  7.12_wp
  ystatgp(25) = 'MANNHEIM          '; latstat(25) = 49.52_wp; lonstat(25) =  8.35_wp
  ystatgp(26) = 'NUERNBERG         '; latstat(26) = 49.50_wp; lonstat(26) = 11.08_wp
  ystatgp(27) = 'KUEMMERSBRUCK     '; latstat(27) = 49.43_wp; lonstat(27) = 11.90_wp
  ystatgp(28) = 'STUTTGART         '; latstat(28) = 48.83_wp; lonstat(28) =  9.20_wp
  ystatgp(29) = 'SIGMARINGEN       '; latstat(29) = 48.10_wp; lonstat(29) =  9.25_wp
  ystatgp(30) = 'ULM               '; latstat(30) = 48.40_wp; lonstat(30) =  9.98_wp
  ystatgp(31) = 'OBERPFAFFENHOFEN  '; latstat(31) = 48.08_wp; lonstat(31) = 11.28_wp
  ystatgp(32) = 'MUENCHEN          '; latstat(32) = 48.25_wp; lonstat(32) = 11.55_wp
  ystatgp(33) = 'DONAUESCHINGEN    '; latstat(33) = 47.97_wp; lonstat(33) =  8.52_wp
  ystatgp(34) = 'FRIEDRICHSHAFEN   '; latstat(34) = 47.67_wp; lonstat(34) =  9.52_wp
  ystatgp(35) = 'ALTENSTADT        '; latstat(35) = 47.83_wp; lonstat(35) = 10.87_wp
  ystatgp(36) = 'HOHENPEISSENBERG  '; latstat(36) = 47.80_wp; lonstat(36) = 11.02_wp
 
  ! Eastern Alpine countries
  ystatgp(37) = 'PAYERNE           '; latstat(37) = 46.82_wp; lonstat(37) =  6.95_wp
  ystatgp(38) = 'ZURICH-KLOTEN     '; latstat(38) = 47.48_wp; lonstat(38) =  8.53_wp
  ystatgp(39) = 'GENEVE-COINTRIN   '; latstat(39) = 46.25_wp; lonstat(39) =  6.13_wp
  ystatgp(40) = 'LINZ              '; latstat(40) = 48.23_wp; lonstat(40) = 14.18_wp
  ystatgp(41) = 'INNSBRUCK         '; latstat(41) = 47.27_wp; lonstat(41) = 11.35_wp
  ystatgp(42) = 'WIEN              '; latstat(42) = 48.25_wp; lonstat(42) = 16.37_wp
  ystatgp(43) = 'GRAZ              '; latstat(43) = 47.00_wp; lonstat(43) = 15.43_wp
  ystatgp(44) = 'LJUBLJANA         '; latstat(44) = 46.07_wp; lonstat(44) = 14.52_wp
  ystatgp(45) = 'ZAGREB            '; latstat(45) = 45.82_wp; lonstat(45) = 16.03_wp
  ystatgp(46) = 'MILANO            '; latstat(46) = 45.43_wp; lonstat(46) =  9.27_wp
  ystatgp(47) = 'UDINE             '; latstat(47) = 46.03_wp; lonstat(47) = 13.18_wp
  ystatgp(48) = 'BOLOGNA           '; latstat(48) = 44.65_wp; lonstat(48) = 11.62_wp
  ystatgp(49) = 'PRATICA           '; latstat(49) = 41.65_wp; lonstat(49) = 12.43_wp
  ystatgp(50) = 'BRINDISI          '; latstat(50) = 40.65_wp; lonstat(50) = 17.95_wp
  ystatgp(51) = 'CAGLIARI          '; latstat(51) = 39.25_wp; lonstat(51) =  9.07_wp
  ystatgp(52) = 'TRAPANI           '; latstat(52) = 37.92_wp; lonstat(52) = 12.50_wp
 
  ! France and Iberia
  ystatgp(53) = 'BREST             '; latstat(53) = 48.45_wp; lonstat(53) = -4.42_wp
  ystatgp(54) = 'TRAPPES           '; latstat(54) = 48.77_wp; lonstat(54) =  2.02_wp
  ystatgp(55) = 'NANCY             '; latstat(55) = 48.68_wp; lonstat(55) =  6.22_wp
  ystatgp(56) = 'BOURGES           '; latstat(56) = 47.07_wp; lonstat(56) =  2.37_wp
  ystatgp(57) = 'LYON              '; latstat(57) = 45.73_wp; lonstat(57) =  5.08_wp
  ystatgp(58) = 'BORDEAUX          '; latstat(58) = 44.83_wp; lonstat(58) = -0.70_wp
  ystatgp(59) = 'NIMES             '; latstat(59) = 43.87_wp; lonstat(59) =  4.40_wp
  ystatgp(60) = 'AJACCIO           '; latstat(60) = 41.92_wp; lonstat(60) =  8.80_wp
  ystatgp(61) = 'LA_CORUNA         '; latstat(61) = 43.37_wp; lonstat(61) = -8.42_wp
  ystatgp(62) = 'SANTANDER         '; latstat(62) = 43.47_wp; lonstat(62) = -3.82_wp
  ystatgp(63) = 'ZARAGOZA          '; latstat(63) = 41.67_wp; lonstat(63) = -1.02_wp
  ystatgp(64) = 'MADRID            '; latstat(64) = 40.50_wp; lonstat(64) = -3.58_wp
  ystatgp(65) = 'PALMA_DE_MALLORCA '; latstat(65) = 39.55_wp; lonstat(65) =  2.62_wp
  ystatgp(66) = 'LISBOA            '; latstat(66) = 38.77_wp; lonstat(66) = -9.13_wp
  ystatgp(67) = 'MURCIA            '; latstat(67) = 38.00_wp; lonstat(67) = -1.17_wp
  ystatgp(68) = 'GIBRALTAR         '; latstat(68) = 36.25_wp; lonstat(68) = -5.55_wp
 
  ! Great Britain
  ystatgp(69) = 'LERWICK           '; latstat(69) = 60.13_wp; lonstat(69) = -1.18_wp
  ystatgp(70) = 'STORNOWAY         '; latstat(70) = 58.22_wp; lonstat(70) = -6.32_wp
  ystatgp(71) = 'BOULMER           '; latstat(71) = 55.42_wp; lonstat(71) = -1.60_wp
  ystatgp(72) = 'ESKMEALS          '; latstat(72) = 54.32_wp; lonstat(72) = -3.40_wp
  ystatgp(73) = 'LONG_KESH         '; latstat(73) = 54.48_wp; lonstat(73) = -6.10_wp
  ystatgp(74) = 'AUGHTON           '; latstat(74) = 53.55_wp; lonstat(74) = -2.92_wp
  ystatgp(75) = 'NOTTINGHAM        '; latstat(75) = 53.00_wp; lonstat(75) = -1.25_wp
  ystatgp(76) = 'CAPEL_DEWI        '; latstat(76) = 52.42_wp; lonstat(76) = -4.00_wp
  ystatgp(77) = 'ABERPORTH         '; latstat(77) = 52.13_wp; lonstat(77) = -4.57_wp
  ystatgp(78) = 'HEMSBY            '; latstat(78) = 52.68_wp; lonstat(78) =  1.68_wp
  ystatgp(79) = 'VALENTIA          '; latstat(79) = 51.93_wp; lonstat(79) =-10.25_wp
  ystatgp(80) = 'LARKHILL          '; latstat(80) = 51.20_wp; lonstat(80) = -1.80_wp
  ystatgp(81) = 'CRAWLEY           '; latstat(81) = 51.08_wp; lonstat(81) =  0.22_wp
  ystatgp(82) = 'SHOEBURYNESS      '; latstat(82) = 51.55_wp; lonstat(82) =  0.83_wp
  ystatgp(83) = 'CAMBORNE          '; latstat(83) = 50.22_wp; lonstat(83) = -5.32_wp
  ystatgp(84) = 'HERSTMONCEUX      '; latstat(84) = 50.90_wp; lonstat(84) =  0.32_wp
 
  ! Western Scandinavia and Benelux
  ystatgp(85) = 'THORSHAVN         '; latstat(85) = 62.02_wp; lonstat(85) = -6.77_wp
  ystatgp(86) = 'BODO              '; latstat(86) = 67.25_wp; lonstat(86) = 14.40_wp
  ystatgp(87) = 'LULEA             '; latstat(87) = 65.55_wp; lonstat(87) = 22.13_wp
  ystatgp(88) = 'ORLAND            '; latstat(88) = 63.70_wp; lonstat(88) =  9.60_wp
  ystatgp(89) = 'SUNDSVALL         '; latstat(89) = 62.53_wp; lonstat(89) = 17.45_wp
  ystatgp(90) = 'OSLO              '; latstat(90) = 60.20_wp; lonstat(90) = 11.10_wp
  ystatgp(91) = 'STAVANGER         '; latstat(91) = 58.87_wp; lonstat(91) =  5.67_wp
  ystatgp(92) = 'GOTEBORG          '; latstat(92) = 57.67_wp; lonstat(92) = 12.30_wp
  ystatgp(93) = 'VISBY             '; latstat(93) = 57.65_wp; lonstat(93) = 18.65_wp
  ystatgp(94) = 'EKOFISK           '; latstat(94) = 56.53_wp; lonstat(94) =  3.22_wp
  ystatgp(95) = 'ALBORG            '; latstat(95) = 57.10_wp; lonstat(95) =  9.88_wp
  ystatgp(96) = 'KOPENHAGEN        '; latstat(96) = 55.77_wp; lonstat(96) = 12.53_wp
  ystatgp(97) = 'VLIELAND          '; latstat(97) = 53.25_wp; lonstat(97) =  4.92_wp
  ystatgp(98) = 'AMSTERDAM-SCHIPHOL'; latstat(98) = 52.30_wp; lonstat(98) =  4.77_wp
  ystatgp(99) = 'VALKENBURG        '; latstat(99) = 52.18_wp; lonstat(99) =  4.42_wp
  ystatgp(100)= 'DE_BILT           '; latstat(100)= 52.10_wp; lonstat(100)=  5.18_wp
  ystatgp(101)= 'UCCLE             '; latstat(101)= 50.80_wp; lonstat(101)=  4.35_wp
  ystatgp(102)= 'ST-HUBERT         '; latstat(102)= 50.03_wp; lonstat(102)=  5.40_wp
 
  ! East-central to southeastern Europe
  ystatgp(103)= 'LEBA              '; latstat(103)= 54.75_wp; lonstat(103)= 17.53_wp
  ystatgp(104)= 'POZNAN            '; latstat(104)= 52.42_wp; lonstat(104)= 16.83_wp
  ystatgp(105)= 'LEGIONOWO         '; latstat(105)= 52.40_wp; lonstat(105)= 20.97_wp
  ystatgp(106)= 'WROCLAW           '; latstat(106)= 51.13_wp; lonstat(106)= 16.98_wp
  ystatgp(107)= 'PRAHA             '; latstat(107)= 50.02_wp; lonstat(107)= 14.45_wp
  ystatgp(108)= 'BRNO              '; latstat(108)= 49.10_wp; lonstat(108)= 16.64_wp
  ystatgp(109)= 'POPRAD            '; latstat(109)= 49.05_wp; lonstat(109)= 20.53_wp
  ystatgp(110)= 'BUDAPEST          '; latstat(110)= 47.43_wp; lonstat(110)= 19.18_wp
  ystatgp(111)= 'SZEGED            '; latstat(111)= 46.25_wp; lonstat(111)= 20.10_wp
  ystatgp(112)= 'CLUJ              '; latstat(112)= 46.78_wp; lonstat(112)= 23.57_wp
  ystatgp(113)= 'SOFIA             '; latstat(113)= 42.65_wp; lonstat(113)= 23.38_wp
  ystatgp(114)= 'THESSALONIKI      '; latstat(114)= 40.52_wp; lonstat(114)= 22.97_wp
  ystatgp(115)= 'ATHENS            '; latstat(115)= 37.90_wp; lonstat(115)= 23.73_wp
 
  ! Northeastern Europe
  ystatgp(116)= 'SODANKYLA         '; latstat(116)= 67.37_wp; lonstat(116)= 26.65_wp
  ystatgp(117)= 'KANDALAKSA        '; latstat(117)= 67.15_wp; lonstat(117)= 32.35_wp
  ystatgp(118)= 'JYVASKYLA         '; latstat(118)= 62.40_wp; lonstat(118)= 25.68_wp
  ystatgp(119)= 'JOKIOINEN         '; latstat(119)= 60.82_wp; lonstat(119)= 23.50_wp
  ystatgp(120)= 'TALLIN            '; latstat(120)= 59.38_wp; lonstat(120)= 24.80_wp
  ystatgp(121)= 'ST.PETERBURG      '; latstat(121)= 59.95_wp; lonstat(121)= 30.70_wp
  ystatgp(122)= 'PSKOV             '; latstat(122)= 57.82_wp; lonstat(122)= 28.42_wp
  ystatgp(123)= 'RIGA              '; latstat(123)= 56.97_wp; lonstat(123)= 24.05_wp
  ystatgp(124)= 'VELIKIE_LUKI      '; latstat(124)= 56.35_wp; lonstat(124)= 30.62_wp
  ystatgp(125)= 'KALININGRAD       '; latstat(125)= 54.72_wp; lonstat(125)= 20.55_wp
  ystatgp(126)= 'KAUNAS            '; latstat(126)= 54.88_wp; lonstat(126)= 23.83_wp
  ystatgp(127)= 'MINSK             '; latstat(127)= 53.93_wp; lonstat(127)= 27.63_wp
  ystatgp(128)= 'BREST-LITOVSK     '; latstat(128)= 52.12_wp; lonstat(128)= 23.68_wp
  ystatgp(129)= 'LVIV              '; latstat(129)= 49.82_wp; lonstat(129)= 23.95_wp
  ystatgp(130)= 'UZHHOROD          '; latstat(130)= 48.63_wp; lonstat(130)= 22.27_wp
  ystatgp(131)= 'CHERNIVTSI        '; latstat(131)= 48.37_wp; lonstat(131)= 25.90_wp

!------------------------------------------------------------------------------
!  Section 3: Initializations for gridpoint calculations
!------------------------------------------------------------------------------

!----------------------------------------------------------------------------
!  Section 3.1: Allocation of space
!----------------------------------------------------------------------------

  IF (ngp > 0) THEN
    CALL alloc_gridpoints  (istat)
  ENDIF

!----------------------------------------------------------------------------
!  Section 3.2: Initialize values for constant fields
!----------------------------------------------------------------------------

  DO n = 1,ngp
    izgp=stationlist(n)%igp
    jzgp=stationlist(n)%jgp
    gphsur(n) = hsurf  (izgp,jzgp)
    gpfrla(n) = fr_land(izgp,jzgp)
    gpsoil(n) = soiltyp(izgp,jzgp)
    gplati(n) = rlat   (izgp,jzgp)
    gplong(n) = rlon   (izgp,jzgp)
    DO k = 1, ke
      gpp0(n,k) = p0(izgp,jzgp,k)
    ENDDO
  ENDDO

  IF (lgplong .OR. lgpspec) THEN
    DO n = 1,ngp
      izgp=stationlist(n)%igp
      jzgp=stationlist(n)%jgp

      gpfc  (n) = fc     (izgp,jzgp)
      gphcan(n) = h_can  (izgp,jzgp)
      gpdpat(n) = d_pat  (izgp,jzgp)
      gpsai(n)  = sai    (izgp,jzgp) 

      IF (lphys .AND. ltur) THEN
        gprdrg(n) = 1.0_wp-EXP(r_air(izgp,jzgp,ke1))
        gpcbig(n) = c_big(izgp,jzgp,ke1)
        gpcsml(n) = c_sml(izgp,jzgp,ke1)
      ENDIF

      DO k = 1, ke1
        gphhl (k,n) =  hhl(izgp,jzgp,k)
        IF (k < ke1)                                         &
        gphml (k,n) = (hhl(izgp,jzgp,k) + hhl(izgp,jzgp,k+1)) * 0.5_wp
      ENDDO

    ENDDO
  ENDIF

!------------------------------------------------------------------------------
!  Section 3.3: Definition (and opening) of the file(s)
!------------------------------------------------------------------------------

  IF (lgplong .OR. lgpshort .OR. lgpspec) THEN
    ! initialize file names and units
    yugpname(:) = '                      '
    nugpunit(:) = -1

    ! open a file for every grid point: loop over all grid points
    DO n = 1,ngp
      ! for further references
      l = ntot_in(n)

      ! test, whether a station name is given in the Namelist input
      ! and use this as file name
      IF (stationlist(n)%ystation_name(1:1) /= '-') THEN
        izlen = LEN_TRIM(stationlist(n)%ystation_name)
        yugpname(n)(1:izlen+2) = 'M_'//TRIM(stationlist(n)%ystation_name)
        ystation = TRIM(stationlist(n)%ystation_name)

      ! else, test, whether a station name is available
      ! then the station name will be used as file name
      ELSE
        zchkgp  =  MAX (dlat,dlon)
        station_loop: DO i = 1, nstatgp
          IF (      (ABS (latstat(i) - gplati(n)*raddeg) < zchkgp)         &
              .AND. (ABS (lonstat(i) - gplong(n)*raddeg) < zchkgp) ) THEN
            ystation = TRIM(ystatgp(i))
            yugpname(n)(1:LEN_TRIM(ystatgp(i))+2) = 'M_'//TRIM(ystatgp(i))
            EXIT station_loop
          ELSE
            ystation = '               '
          ENDIF
        ENDDO station_loop

        ! if no station name is available, longitude and latitude values
        ! are used instead (in 1/10 degree)
        IF (yugpname(n)(1:1) == ' ') THEN
          yugpname(n)(1:12) = 'M_0000_x000_'

          irlat = INT (ABS (gplati(n)*raddeg * 10.0_wp), iintegers)
          irlon = INT (ABS (gplong(n)*raddeg * 10.0_wp), iintegers)

          IF (gplati(n) < 0.0_wp) THEN
            yugpname(n)(12:12) = 'S'
          ELSE
            yugpname(n)(12:12) = 'N'
          ENDIF

          IF (gplong(n) < 0.0_wp) THEN
            yugpname(n)(7:7) = 'W'
          ELSE
            yugpname(n)(7:7) = 'E'
          ENDIF

          IF     (irlon <   10) THEN
            WRITE (yugpname(n)(6:6),'(I1)') irlon
          ELSEIF (irlon >=  10 .AND. irlon <  100) THEN
            WRITE (yugpname(n)(5:6),'(I2)') irlon
          ELSEIF (irlon >= 100 .AND. irlon < 1000) THEN
            WRITE (yugpname(n)(4:6),'(I3)') irlon
          ELSE
            WRITE (yugpname(n)(3:6),'(I4)') irlon
          ENDIF

          IF     (irlat <   10) THEN
            WRITE (yugpname(n)(11:11),'(I1)') irlat
          ELSEIF (irlat >=  10 .AND. irlat <  100) THEN
            WRITE (yugpname(n)(10:11),'(I2)') irlat
          ELSE
            WRITE (yugpname(n)( 9:11),'(I3)') irlat
          ENDIF
        ENDIF
      ENDIF

      ! get a unit number for every file
      CALL get_free_unit (nugpunit(n))

      ! write a headline in the corresponding file
      OPEN (nugpunit(n), FILE=TRIM(yugpname(n)), FORM='FORMATTED',        &
                          STATUS='UNKNOWN', POSITION='APPEND', IOSTAT=nstat)
      IF (nstat /= 0) THEN
        yerrmsg = 'opening of file '//TRIM(yugpname(n))//' failed'
        CALL model_abort (my_cart_id, 7004, yerrmsg, yroutine)
      ENDIF

      ! headline for a grid point
      CALL get_utc_date (0, ydate_ini, dt, itype_calendar, ydat11, ydat12, &
                         njata1, hh1)

      IF (lgplong) THEN
        WRITE (nugpunit(n), '(A,A,A,I5,A,I5)') 'GRID POINT:        ',      &
                     ystation, '    I: ',stationlist_tot(l)%igp,           &
                               '    J: ',stationlist_tot(l)%jgp
        WRITE (nugpunit(n), '(A,A)')    'Initial Date:      ',ydat12
        WRITE (nugpunit(n), '(A,F9.3)')                                    &
                            '    HSURF   ( m ): ', gphsur(n)
        WRITE (nugpunit(n), '(A,F9.3)')                                    &
                            '    FR_LAND ( % ): ', gpfrla(n)*100.0_wp
        WRITE (nugpunit(n), '(A,F9.3)')                                    &
                            '    LAT     (dgr): ', gplati(n)*raddeg
        WRITE (nugpunit(n), '(A,F9.3)')                                    &
                            '    LON     (dgr): ', gplong(n)*raddeg
        WRITE (nugpunit(n), '(A,F9.3)')                                    &
                            '    FC    (1E4/s): ', gpfc  (n)*1.0E4_wp
        WRITE (nugpunit(n), '(A,A,I6)')                                    &
            '    SOIL TYPE    : ', ysoiltyp(NINT(gpsoil(n))), NINT(gpsoil(n))
        WRITE (nugpunit(n), '(A)')      '    '
     
      ELSEIF (lgpshort) THEN

        ! Headline for the output
        WRITE (nugpunit(n), '(A,A)')                                          &
        '          *** Model: COSMO ***        Start of the forecast: ',ydat12
        WRITE (nugpunit(n), '(A)')                                            &
        '          Short meteograph of the COSMO-Model at selected grid points'
        WRITE (nugpunit(n), '(A)') '   '
      
        WRITE (nugpunit(n), '(A)')                                            &
        '     Meaning of the parameters:'
        WRITE (nugpunit(n), '(A)')                                            &
        '       HH:       forecast hour                               (hours)'
        WRITE (nugpunit(n), '(A)')                                            &
        '       PS:       mean sea level pressure                       (hpa)'
        WRITE (nugpunit(n), '(A)')                                            &
        '       DF10M:    wind direction and speed at 10m      (degree/ m/s )'
        WRITE (nugpunit(n), '(A)')                                            &
        '       DF500M:   wind direction and speed at 950 hPa  (degree/ m/s )'
        WRITE (nugpunit(n), '(A)')                                            &
        '       DF850:    wind direction and speed at 850 hPa  (degree/ m/s )'
        WRITE (nugpunit(n), '(A)')                                            &
        '       DF700     wind direction and speed at 700 hPa  (degree/ m/s )'
        WRITE (nugpunit(n), '(A)')                                            &
        '       DF500     wind direction and speed at 500 hPa  (degree/ m/s )'
        WRITE (nugpunit(n), '(A)') '  '
        WRITE (nugpunit(n), '(A)')                                            &
        '       TG:       ground temperature                              (C)'
        WRITE (nugpunit(n), '(A)')                                            &
        '       T2M:      temperature at 2m                               (C)'
        WRITE (nugpunit(n), '(A)')                                            &
        '       TD2M:     dew point temperature at 2m                     (C)'
        WRITE (nugpunit(n), '(A)')                                            &
        '       T30M:     temperature at ~ 30m (lowest level)             (C)'
        WRITE (nugpunit(n), '(A)')                                            &
        '       T850:     temperature at 850 hpa                          (C)'
        WRITE (nugpunit(n), '(A)')                                            &
        '       T700:     temperature at 700 hpa                          (C)'
        WRITE (nugpunit(n), '(A)')                                            &
        '       T500:     temperature at 500 hpa                          (C)'
        WRITE (nugpunit(n), '(A)') '  '
        WRITE (nugpunit(n), '(A)')                                            &
        '       HML:      cloud cover (high, medium, low)              (0..8)'
        WRITE (nugpunit(n), '(A)')                                            &
        '        = :      ground fog                                   (0..8)'
        WRITE (nugpunit(n), '(A)')                                            &
        '       HBAS:     base height of con. cloud above msl          (m*10)'
        WRITE (nugpunit(n), '(A)')                                            &
        '       HTOP      top  height of con. cloud above msl          (m*10)'
        WRITE (nugpunit(n), '(A)')                                            &
        '       RR:       rain amount (grid scale and convective)        (mm)'
        WRITE (nugpunit(n), '(A)')                                            &
        '       RS:       snow amount (grid scale and convective)        (mm)'
        WRITE (nugpunit(n), '(A)')                                            &
        '       WS:       water content of snow                           (m)'
        WRITE (nugpunit(n), '(A)') '  '
        WRITE (nugpunit(n), '(A)') '  '

        ! headline for the variables
        WRITE (nugpunit(n), '(A)') '  '
        WRITE (nugpunit(n), '(A)') '  '
        WRITE (nugpunit(n), '(A,A,A,I5,A,I5)') 'GRID POINT:        ',      &
                     ystation, '    I: ',stationlist_tot(l)%igp,           &
                               '    J: ',stationlist_tot(l)%jgp
        WRITE (nugpunit(n), '(A,A)')    'Initial Date:      ',ydat12
        WRITE (nugpunit(n), '(A,F9.3)') '    HSURF   ( m ): ', gphsur(n)
        WRITE (nugpunit(n), '(A,F9.3)') '    FR_LAND ( % ): ', gpfrla(n)*100.0_wp
        WRITE (nugpunit(n), '(A,F9.3)') '    LAT     (dgr): ', gplati(n)*raddeg
        WRITE (nugpunit(n), '(A,F9.3)') '    LON     (dgr): ', gplong(n)*raddeg
        WRITE (nugpunit(n), '(2A,I6)')  '    SOIL TYPE    : ',             &
                                    ysoiltyp(NINT(gpsoil(n))), NINT(gpsoil(n))
        WRITE (nugpunit(n), '(A)')      '    '

! if SNOW_MELT should be included in short meteographs
!       WRITE (nugpunit(n), '(A,A,A)')                                        &
!         '    HH   PMSL    DF10M  DF500M   DF850   DF700   DF500   TG   T2M',&
!         '   TD2M  T30M  T850  T700  T500 HML =  HBAS HTOP  RR    RS     WS',&
!         '  SNOW_MELT'
!       WRITE (nugpunit(n), '(A,A)')                                          &
!         '     h   hpa   dgr m/s dgr m/s  dgr m/s dgr m/s dgr m/s          ',&
!         '    degree centrigrade         octas     10*m        mm         m  mm'

        WRITE (nugpunit(n), '(A,A)')                                          &
          '    HH   PMSL    DF10M  DF500M   DF850   DF700   DF500   TG   T2M',&
          '   TD2M  T30M  T850  T700  T500 HML =  HBAS HTOP  RR    RS     WS'
        WRITE (nugpunit(n), '(A,A)')                                          &
          '     h   hpa   dgr m/s dgr m/s  dgr m/s dgr m/s dgr m/s          ',&
          '    degree centrigrade         octas     10*m        mm         m'
        WRITE (nugpunit(n), '(A)') '  '

      ELSEIF (lgpspec) THEN

        IF (ystation(1:1).EQ.' ') THEN
          WRITE(nugpunit(n),'("Location ",a)') yugpname(n)
        ELSE
          WRITE(nugpunit(n),'("Location ",a)') ystation
        END IF
        WRITE (nugpunit(n),*)

      ENDIF

      CLOSE (nugpunit(n), STATUS='KEEP')
    ENDDO

  ENDIF

!------------------------------------------------------------------------------
!  End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE init_gridpoints

!==============================================================================
!==============================================================================
!+ Module procedure to store values for selected gridpoints in special variable
!------------------------------------------------------------------------------

SUBROUTINE gridpoints (ydate_ini, ntl)

!------------------------------------------------------------------------------
!
! Description:
!   This routine stores values for selected grid points in special variables.
!   These variables are written to a special file at the end of the forecast.
!
! Method:
!
!------------------------------------------------------------------------------

! Parameterlist
CHARACTER (LEN=14),       INTENT(IN)    ::       &
  ydate_ini      ! start of the forecast  yyyymmddhh (year, month, day, hour)

INTEGER (KIND=iintegers), INTENT(IN)    ::       &
  ntl            ! time level for which to perform the grid point output
!------------------------------------------------------------------------------

! Local variables
REAL (KIND=wp)           :: zvb(ke), zrf(ke), rzpm(ke), rzph(ke1)
REAL (KIND=wp)           ::     &
  zphf, zupu, zupo, zvpu, zvpo, ztpu, ztpo, zpu, zpo, zpx,                  &
  zu10m, zu500, zu700, zu850, zu950, zv10m, zv500, zv700, zv850, zv950,     &
  zfgew, zfgqd, zt, zge, zp, zgew, zpr, zrr, zth2, hh2,                     &
  zpb, zpm, zhb, zhm, dzh, ztv, d_hor, dzw(4)

INTEGER (KIND=iintegers) :: izpm(ke), izph(ke1)
INTEGER (KIND=iintegers) ::     &
  mzg, izgp, jzgp, izgpl, jzgpu, nzpa, k, k500m, k850, k700, k500, kbas,    &
  ktop, nstat, i, izch  , izcm  , izcl, izd10m, izd500, izd700, izd850,     &
  izd950, izf10m, izf500, izf700, izf850, izf950, izfog , izpsr,            &
  njata2, kso, k1, l, nlen, i_hor, prec

CHARACTER (LEN=99)       ::  yline (ke1+32)
CHARACTER (LEN=14)       ::  ydat11, ydat21
CHARACTER (LEN=28)       ::  ydat12, ydat22

CHARACTER (LEN=20)       ::  yroutine
CHARACTER (LEN=255)      ::  yerrmsg

CHARACTER (LEN= 8)       :: form0
CHARACTER (LEN=10)       :: form1
CHARACTER (LEN=11)       :: form2, form3
CHARACTER (LEN=20)       :: dezi

LOGICAL   :: flat, lhydpres, lpreslev

REAL (KIND=wp),     POINTER :: &
  qv  (:,:,:)=> NULL(),        &   ! QV at ntl
  qc  (:,:,:)=> NULL(),        &   ! QC at ntl
  qi  (:,:,:)=> NULL()             ! QI at ntl

#ifdef COSMOART
INTEGER (KIND=iintegers) :: isp, loc_isp
REAL (KIND=wp),     POINTER :: &
  cash(:,:,:,:)=> NULL(),      &    ! cash at ntl
  caero(:,:,:,:)=> NULL()           ! caero at ntl
#endif

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

! Statementfunctions
zfgew (  zt  ) = b1 * EXP (b2w * (zt - b3)/(zt - b4w))
zfgqd (zge,zp) = rdv * zge / (zp - o_m_rdv * zge)

!------------------------------------------------------------------------------
!- Begin SUBROUTINE gridpoints
!------------------------------------------------------------------------------

  nstat = 0 
  yroutine = 'gridpoints'

  ! Calculate the gridpoint time index for the present time step
  IF (lgpspec) THEN
    IF (nincgp == 1) THEN
      nzpa = (ntstep + 1 - n0gp) / nincgp 
    ELSE
      nzpa = (ntstep + 1 - n0gp) / nincgp + 1
    ENDIF
  ELSEIF (lgplong .OR. lgpshort) THEN
    nzpa = 1
  ENDIF

  ! Get actual date
  zth2    = ntstep * dtdeh
  CALL get_utc_date (     0, ydate_ini, dt, itype_calendar, ydat11, ydat12, &
                     njata2, hh2)
  CALL get_utc_date (ntstep, ydate_ini, dt, itype_calendar, ydat21, ydat22, &
                     njata2, hh2)

  ! retrieve the required microphysics tracers (at timelevel ntl)
  CALL trcr_get(nstat, idt_qv, ptr_tlev = ntl, ptr = qv)
  IF (nstat /= 0) THEN
    yerrmsg = trcr_errorstr(nstat)
    CALL model_abort(my_cart_id, nstat, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(nstat, idt_qc, ptr_tlev = ntl, ptr = qc)
  IF (nstat /= 0) THEN
    yerrmsg = trcr_errorstr(nstat)
    CALL model_abort(my_cart_id, nstat, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(nstat, idt_qi, ptr_tlev = ntl, ptr = qi)
  IF (nstat /= 0 .AND. nstat /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(nstat)
    CALL model_abort(my_cart_id, nstat, yerrmsg, yroutine)
  ENDIF

#ifdef COSMOART
  IF(l_cosmo_art) THEN
    IF (lvolcano) THEN
      CALL trcr_get_block(nstat, idx_start=trcr_idx_ash(1), idx_end=trcr_idx_ash(isp_ash), &
                          ptr_tlev=ntl, ptr=cash)
      IF (nstat /= 0) THEN
        yerrmsg = trcr_errorstr(nstat)
        CALL model_abort(my_cart_id, nstat, yerrmsg, yroutine)
      ENDIF
    ELSE IF (ldust .AND. lonly_dust) THEN
      CALL trcr_get_block(nstat, idx_start=trcr_idx_aero(vsoila), idx_end=trcr_idx_aero(vsoilc0), &
                          ptr_tlev=ntl, ptr=caero)
      IF (nstat /= 0) THEN
        yerrmsg = trcr_errorstr(nstat)
        CALL model_abort(my_cart_id, nstat, yerrmsg, yroutine)
      ENDIF
    ENDIF
  ENDIF
#endif

!------------------------------------------------------------------------------
!- Section 1: Store variables for long or special gridpoint output  
!------------------------------------------------------------------------------

IF (lgplong .OR. lgpspec) THEN

  ! Loop over all chosen gridpoints
  DO  mzg = 1, ngp

    izgp   = stationlist(mzg)%igp
    jzgp   = stationlist(mzg)%jgp
    izgpl  = MAX(1,izgp-1)
    jzgpu  = MAX(1,jzgp-1)

    ! variables that vary in climate runs
    gpplco(mzg,nzpa) = plcov  (izgp,jzgp)
    gplai (mzg,nzpa) = lai    (izgp,jzgp)
    gproot(mzg,nzpa) = rootdp (izgp,jzgp)
    gpvio3(mzg,nzpa) = vio3   (izgp,jzgp)
    gphmo3(mzg,nzpa) = hmo3   (izgp,jzgp)

    ! Prognostic and diagnostic variables
    DO  k  = 1, ke
      gpu (k,mzg,nzpa)   = 0.5_wp*(u(izgp ,jzgp , k,ntl)+u(izgpl,jzgp ,k,ntl))
      gpv (k,mzg,nzpa)   = 0.5_wp*(v(izgp ,jzgp , k,ntl)+v(izgp ,jzgpu,k,ntl))
      gpw (k,mzg,nzpa)   = w    (izgp, jzgp, k, ntl)
      gpt (k,mzg,nzpa)   = t    (izgp, jzgp, k, ntl)
      gpqv(k,mzg,nzpa)   = qv   (izgp, jzgp, k     )
      gpqc(k,mzg,nzpa)   = qc   (izgp, jzgp, k     )
      IF ( ASSOCIATED(qi) ) THEN
        gpqi(k,mzg,nzpa) = qi   (izgp, jzgp, k     )
      ELSE
        gpqi(k,mzg,nzpa) = 0.0_wp
      ENDIF
      gppp(k,mzg,nzpa)   = pp   (izgp, jzgp, k, ntl)
      gpclcs(k,mzg,nzpa) = clc_sgs (izgp,jzgp,k)
      gpclcc(k,mzg,nzpa) = clc_con (izgp,jzgp,k)
#ifdef COSMOART
      IF (lgplong) THEN
        IF(l_cosmo_art) THEN
          IF (lvolcano) THEN
            gpashmc(k,mzg,nzpa) = 0.0_wp
            DO isp = 1, isp_ash
              gpcash(k,mzg,nzpa,isp) = cash(izgp,jzgp,k,isp)
              gpashmc(k,mzg,nzpa) = gpashmc(k,mzg,nzpa)                         &
                                  + ash_scale(isp) * gpcash(k,mzg,nzpa,isp)
            ENDDO
          ELSE IF (ldust .AND. lonly_dust) THEN
            DO isp = 1, isp_aerotrans
              gpcaero(k,mzg,nzpa,isp) = caero(izgp,jzgp,k,isp)
            ENDDO
            IF (lrad_dust) THEN
               gpcaero(k,mzg,nzpa,isp_aerotrans+1) = tau_dust_550nm(izgp,jzgp,k)
            END IF
          ENDIF
        ENDIF
      ENDIF
#endif

    ENDDO
    gpw (ke1,mzg,nzpa) = w    (izgp, jzgp, ke1, ntl)

    gpps   (mzg,nzpa)  = ps    (izgp, jzgp, ntl)
    gpt_g  (mzg,nzpa)  = t_g   (izgp, jzgp, ntl)
    gpqvs  (mzg,nzpa)  = qv_s  (izgp, jzgp, ntl)
    gpdpdt (mzg,nzpa)  = dpsdt (izgp, jzgp      )
    gpshfl (mzg,nzpa)  = shfl_s(izgp, jzgp      )
    gplhfl (mzg,nzpa)  = lhfl_s(izgp, jzgp      )
    gpumfl (mzg,nzpa)  = umfl_s(izgp, jzgp      )
    gpvmfl (mzg,nzpa)  = vmfl_s(izgp, jzgp      )

    ! variables from the turbulence
    DO k = 2, ke
      gptkvm (k,mzg,nzpa) = tkvm  (izgp,jzgp,k)
      gptkvh (k,mzg,nzpa) = tkvh  (izgp,jzgp,k)
    ENDDO

    gptcm (mzg,nzpa) = tcm(izgp,jzgp)
    gptch (mzg,nzpa) = tch(izgp,jzgp)
    gpgz0 (mzg,nzpa) = MAX( gz0 (izgp,jzgp) , 1.0E-10_wp )

    ! variables from grid scale precipitation
    gpprrs(mzg,nzpa) = prr_gsp(izgp,jzgp)
    gpprss(mzg,nzpa) = prs_gsp(izgp,jzgp)

    ! variables from the radiation scheme
    gpswdir_s (mzg,nzpa)  = swdir_s (izgp,jzgp)
    gpswdifd_s(mzg,nzpa)  = swdifd_s(izgp,jzgp)
    gpswdifu_s(mzg,nzpa)  = swdifu_s(izgp,jzgp)
    gpsobs    (mzg,nzpa) = sobs    (izgp,jzgp)
    gpthbs    (mzg,nzpa) = thbs    (izgp,jzgp)
    gppabs    (mzg,nzpa) = pabs    (izgp,jzgp)
    gpalb     (mzg,nzpa) = alb_rad (izgp,jzgp)
    gpsobt    (mzg,nzpa) = sobt    (izgp,jzgp)
    gpthbt    (mzg,nzpa) = thbt    (izgp,jzgp)
    gpclcl    (mzg,nzpa) = clcl    (izgp,jzgp)
    gpclcm    (mzg,nzpa) = clcm    (izgp,jzgp)
    gpclch    (mzg,nzpa) = clch    (izgp,jzgp)
    gpclct    (mzg,nzpa) = clct    (izgp,jzgp)

    ! variables from the convection scheme
    gpprrc(mzg,nzpa) = prr_con(izgp,jzgp)
    gpprsc(mzg,nzpa) = prs_con(izgp,jzgp)

    ! variables from the soil model terra
    IF(lmulti_snow) THEN
      gptsno(mzg,nzpa) = t_snow_mult(izgp,jzgp,1,ntl)
    ELSE
      gptsno(mzg,nzpa) = t_snow(izgp,jzgp,ntl)
    ENDIF
    gpt_s (mzg,nzpa) = t_s (izgp,jzgp,ntl)
    gpwsno(mzg,nzpa) = w_snow(izgp,jzgp,ntl)
    gpsmlt(mzg,nzpa) = snow_melt(izgp,jzgp)
    gpw_i (mzg,nzpa) = w_i (izgp,jzgp,ntl)
    IF (.NOT. lmulti_layer) THEN
      gpt_m (mzg,nzpa) = t_m (izgp,jzgp,ntl)
      gpt_cl(mzg,nzpa) = t_cl(izgp,jzgp)
      gpwg1 (mzg,nzpa) = w_g1(izgp,jzgp,ntl)
      gpwg2 (mzg,nzpa) = w_g2 (izgp,jzgp,ntl)
      IF ( nlgw == 3 )  gpwg3 (mzg,nzpa) = w_g3 (izgp,jzgp,ntl)
      gpw_cl(mzg,nzpa) = w_cl(izgp,jzgp)
    ELSE
      gptso (:,mzg,nzpa) = t_so (izgp,jzgp,:,ntl)
      gpwso (:,mzg,nzpa) = w_so (izgp,jzgp,:,ntl)
      gpwice(:,mzg,nzpa) = w_so_ice(izgp,jzgp,:,ntl)
      gpfrsn  (mzg,nzpa) = freshsnow(izgp,jzgp)
      IF(lmulti_snow) THEN
        gprhsn  (mzg,nzpa) = rho_snow_mult(izgp,jzgp,1,ntl)
      ELSE
        gprhsn  (mzg,nzpa) = rho_snow(izgp,jzgp,ntl)
      ENDIF
      gph_sn  (mzg,nzpa) = h_snow(izgp,jzgp,ntl)
    ENDIF
    gpabsf(mzg,nzpa) = runoff_s (izgp,jzgp)
    gpabgr(mzg,nzpa) = runoff_g (izgp,jzgp)

    ! variables from routine near_surface 
    gpu10m (mzg,nzpa) = u_10m (izgp,jzgp)
    gpv10m (mzg,nzpa) = v_10m (izgp,jzgp)
    gpt2m  (mzg,nzpa) = t_2m  (izgp,jzgp)
    gptd2m (mzg,nzpa) = td_2m (izgp,jzgp)
    gprrsn (mzg,nzpa) = rain_gsp(izgp,jzgp)
    gprssn (mzg,nzpa) = snow_gsp(izgp,jzgp)
    gprrkn (mzg,nzpa) = rain_con(izgp,jzgp)
    gprskn (mzg,nzpa) = snow_con(izgp,jzgp)
    gptmn2 (mzg,nzpa) = tmin_2m (izgp,jzgp)
    gptmx2 (mzg,nzpa) = tmax_2m (izgp,jzgp)
    gpvb10 (mzg,nzpa) = vmax_10m(izgp,jzgp)

  ENDDO
ENDIF

!------------------------------------------------------------------------------
! Section 2:  Write variables for long gridpoint output
!------------------------------------------------------------------------------

IF (lgplong) THEN
  ! Loop over all chosen gridpoints
  DO  mzg = 1, ngp
    ! open the file
    OPEN (nugpunit(mzg), FILE=TRIM(yugpname(mzg)), FORM='FORMATTED',      &
                         STATUS='UNKNOWN', POSITION='APPEND', IOSTAT=nstat)
    IF (nstat /= 0) THEN
      yerrmsg = 'opening of file '//TRIM(yugpname(mzg))//' failed'
      CALL model_abort (my_cart_id, 7004, yerrmsg, yroutine)
    ENDIF

    ! Headline for the first group of variables
    WRITE (nugpunit(mzg), '(A,A,A,F10.2,A,I10,A)')                           &
      '      Actual date:           ',ydat12,' + ',zth2,                     &
      '     (time step:  ',ntstep,')'
    WRITE (nugpunit(mzg),'(A,F6.1,A,F6.1)')                                  &
      '    PS(hpa): ',     gpps  (mzg,nzpa)*0.01_wp,                     &
      '    DPSDT(hpa/h): ',gpdpdt(mzg,nzpa)*36.0_wp
    WRITE(nugpunit(mzg), '(A)')  '  '
#ifdef COSMOART
    IF (l_cosmo_art) THEN
      IF (lvolcano) THEN
        header_line=''
        header_unit=''
        DO isp = 1, isp_ash
          WRITE(header_line,'(A,A,I1)') TRIM(header_line), '          ASH', isp
          WRITE(header_unit,'(A,A)') TRIM(header_unit), '          1/m3'
        END DO
        WRITE(nugpunit(mzg), '(A,A,A,A)')                                           &
          '  K   Pmain    T        QV       QC       QI REL_HUM     CLC   CLC_CON', &
          '      U      V    SPEED', TRIM(header_line), '    ASH_MC    HML'

        WRITE(nugpunit(mzg), '(A,A,A,A)')                                           &
          '       hPa   Grd C     g/kg         mg/kg          %      %       %   ', &
          '            m/s        ', TRIM(header_unit), '     ug/m3       m'
      ELSE IF (ldust .AND. lonly_dust) THEN
        header_line=''
        header_unit=''
        WRITE(header_line,'(A,A)') TRIM(header_line), '        VSOILA'
        WRITE(header_line,'(A,A)') TRIM(header_line), '        VSOILB'
        WRITE(header_line,'(A,A)') TRIM(header_line), '        VSOILC'
        WRITE(header_unit,'(A,A)') TRIM(header_unit), '         ug/m3'
        WRITE(header_unit,'(A,A)') TRIM(header_unit), '         ug/m3'
        WRITE(header_unit,'(A,A)') TRIM(header_unit), '         ug/m3'
        WRITE(header_line,'(A,A)') TRIM(header_line), '       VSOILA0'
        WRITE(header_line,'(A,A)') TRIM(header_line), '       VSOILB0'
        WRITE(header_line,'(A,A)') TRIM(header_line), '       VSOILC0'
        WRITE(header_unit,'(A,A)') TRIM(header_unit), '          1/m3'
        WRITE(header_unit,'(A,A)') TRIM(header_unit), '          1/m3'
        WRITE(header_unit,'(A,A)') TRIM(header_unit), '          1/m3'
        IF (lrad_dust) THEN
          WRITE(header_line,'(A,A)') TRIM(header_line), '  TAU_DUST'
          WRITE(header_unit,'(A,A)') TRIM(header_unit), '         1'
        END IF
        WRITE(nugpunit(mzg), '(A,A,A,A)')                                           &
          '  K   Pmain    T        QV       QC       QI REL_HUM     CLC   CLC_CON', &
          '      U      V    SPEED', TRIM(header_line), '    HML'

        WRITE(nugpunit(mzg), '(A,A,A,A)')                                           &
          '       hPa   Grd C     g/kg         mg/kg          %      %       %   ', &
          '            m/s        ', TRIM(header_unit), '       m'
      ELSE
        WRITE(nugpunit(mzg), '(A,A)')                                            &
          '  K   Pmain    T       QV      QC      QI REL_HUM     CLC   CLC_CON', &
          '      U      V    SPEED    HML'
        WRITE(nugpunit(mzg), '(A,A)')                                            &
          '       hPa   Grd C    g/kg       mg/kg          %      %       %   ', &
          '            m/s             m'
      END IF
    END IF
#else
    WRITE(nugpunit(mzg), '(A,A)')                                            &
      '  K   Pmain    T       QV      QC      QI REL_HUM     CLC   CLC_CON', &
      '      U      V    SPEED    HML'
    WRITE(nugpunit(mzg), '(A,A)')                                            &
      '       hPa   Grd C    g/kg       mg/kg          %      %       %   ', &
      '            m/s             m'
#endif
    WRITE(nugpunit(mzg), '(A)')  '  '

    ! Pressure at main and half levels (hpa)
    DO k = 1 , ke
      rzpm(k) = 0.01_wp * (gpp0(mzg,k) + gppp(k,mzg,nzpa))
      izpm(k) = NINT( rzpm(k) , iintegers )
    ENDDO
    ! Use total pressure computed in rzpm above:
    rzph(1)   = rzpm(1) - 0.5_wp*(rzpm(2)-rzpm(1))
    ! Using only p0 could lead to problems in some situations
    ! rzph(1)   = 0.01_wp *                                                &
    !                (gpp0(mzg,1) - 0.5_wp*(gpp0(mzg,2)-gpp0(mzg,1)))
    izph(1)   = NINT ( rzph(1) , iintegers )
    rzph(ke1) = 0.01_wp * gpps(mzg,nzpa)
    izph(ke1) = NINT ( rzph(ke1) , iintegers )
    DO k = 2 , ke
      rzph(k) = 0.5_wp * (rzpm(k)+rzpm(k-1))
      izph(k) = INT ( rzph(k) , iintegers )
    ENDDO

    ! Wind speed and relative humidity
    DO k = 1 , ke
      zvb(k)  = SQRT (gpu(k,mzg,nzpa)**2 + gpv(k,mzg,nzpa)**2)
      zgew    = zfgew (gpt(k,mzg,nzpa))
      zrf(k)  = ( gpqv(k,mzg,nzpa)+gpqc(k,mzg,nzpa))                           &
                / zfgqd(zgew, rzpm(k)*100.0_wp)
    ENDDO

    ! Output of the first group of variables
#ifdef COSMOART
    IF (l_cosmo_art) THEN
      IF (lvolcano) THEN
        WRITE(header_format, '(A,I1,A)') "(I3,F8.2,F7.2,3F9.3,6F8.2,", isp_ash, "F14.2,F10.2,F9.1)"
        DO k = 1, ke
          WRITE (nugpunit(mzg), TRIM(header_format))       &
               k, rzpm(k), gpt(k,mzg,nzpa)-t0_melt,                              &
               gpqv(k,mzg,nzpa)*1.0E3_wp, gpqc(k,mzg,nzpa)*1.0E6_wp,             &
               gpqi(k,mzg,nzpa)*1.0E6_wp, zrf(k)*100.0_wp,                       &
               gpclcs(k,mzg,nzpa)*100.0_wp, gpclcc(k,mzg,nzpa)*100.0_wp,         &
               gpu(k,mzg,nzpa), gpv(k,mzg,nzpa), zvb(k),                         &
               gpcash(k,mzg,nzpa,1:isp_ash),gpashmc(k,mzg,nzpa),                 &
               gphml(k,mzg)
        ENDDO
      ELSE IF (ldust .AND. lonly_dust) THEN
        IF (lrad_dust) THEN
          WRITE(header_format, '(A,I1,A)') "(I3,F8.2,F7.2,3F9.3,6F8.2,", isp_aerotrans, "F14.2,F10.6,F9.1)"
          loc_isp = isp_aerotrans + 1
        ELSE
          WRITE(header_format, '(A,I1,A)') "(I3,F8.2,F7.2,3F9.3,6F8.2,", isp_aerotrans, "F14.2,F9.1)"
          loc_isp = isp_aerotrans
        END IF
        DO k = 1, ke
          WRITE (nugpunit(mzg), TRIM(header_format))       &
               k, rzpm(k), gpt(k,mzg,nzpa)-t0_melt,                              &
               gpqv(k,mzg,nzpa)*1.0E3_wp, gpqc(k,mzg,nzpa)*1.0E6_wp,             &
               gpqi(k,mzg,nzpa)*1.0E6_wp, zrf(k)*100.0_wp,                       &
               gpclcs(k,mzg,nzpa)*100.0_wp, gpclcc(k,mzg,nzpa)*100.0_wp,         &
               gpu(k,mzg,nzpa), gpv(k,mzg,nzpa), zvb(k),                         &
               gpcaero(k,mzg,nzpa,1:loc_isp),                                    &
               gphml(k,mzg)
        ENDDO
      ELSE
        DO k = 1, ke
          WRITE (nugpunit(mzg), '(I3,F8.2,F7.2,3F8.3,6F8.2,F9.1)')               &
               k, rzpm(k), gpt(k,mzg,nzpa)-t0_melt,                              &
               gpqv(k,mzg,nzpa)*1.0E3_wp, gpqc(k,mzg,nzpa)*1.0E6_wp,             &
               gpqi(k,mzg,nzpa)*1.0E6_wp, zrf(k)*100.0_wp,                       &
               gpclcs(k,mzg,nzpa)*100.0_wp, gpclcc(k,mzg,nzpa)*100.0_wp,         &
               gpu(k,mzg,nzpa), gpv(k,mzg,nzpa), zvb(k), gphml(k,mzg)
        ENDDO
      END IF
    END IF
#else
    DO k = 1, ke
      WRITE (nugpunit(mzg), '(I3,F8.2,F7.2,3F8.3,6F8.2,F9.1)')                 &
             k, rzpm(k), gpt(k,mzg,nzpa)-t0_melt,                              &
             gpqv(k,mzg,nzpa)*1.0E3_wp, gpqc(k,mzg,nzpa)*1.0E6_wp,             &
             gpqi(k,mzg,nzpa)*1.0E6_wp, zrf(k)*100.0_wp,                       &
             gpclcs(k,mzg,nzpa)*100.0_wp, gpclcc(k,mzg,nzpa)*100.0_wp,         &
             gpu(k,mzg,nzpa), gpv(k,mzg,nzpa), zvb(k), gphml(k,mzg)
    ENDDO
#endif
    WRITE (nugpunit(mzg), '(A)')  '  '

    ! Headline for the second group of variables
    WRITE (nugpunit(mzg), '(A,A)')                                           &
      '  K   Phalf       W        TKVM      TKVH    HHL'
    WRITE (nugpunit(mzg), '(A,A)')                                           &
      '       hPa      cm/s           m**2/s         m'
    WRITE (nugpunit(mzg), '(A)')  '  '

    ! Output of the second group of variables
    DO k = 1, ke1
      IF ( (k == 1) .OR. (k == ke1) ) THEN
        WRITE (nugpunit(mzg), '(I3,F8.2,F10.3,A,F9.1)')                      &
         k, rzph(k), gpw(k,mzg,nzpa)*100.0_wp, '                    ',   &
         gphhl(k,mzg)
      ELSE
        WRITE (nugpunit(mzg), '(I3,F8.2,3F10.3,F9.1)')                       &
         k, rzph(k), gpw(k,mzg,nzpa)*100.0_wp, gptkvm(k,mzg,nzpa),       &
                     gptkvh(k,mzg,nzpa), gphhl(k,mzg)
      ENDIF
    ENDDO
    WRITE (nugpunit(mzg), '(A)')  '  '

    ! Prepare output of third group of variables
    ! write everything into yline first
    DO k = 1 , ke1+32
      yline(k) = ' '
    ENDDO

    k = 0
    WRITE (yline(k+ 1),'(A,F12.5)')                                          &
        '  Surface variables:        TCM     : ', zvb(ke)*gptcm (mzg,nzpa)
    WRITE (yline(k+ 2),'(A,F12.5)')                                          &
        '              ( m/s )       TCH     : ', zvb(ke)*gptch (mzg,nzpa)
    WRITE (yline(k+ 3),'(A,F12.5)')                                          &
        '              (  m  )       Z0      : ', gpgz0 (mzg,nzpa)/g
    WRITE (yline(k+ 4),'(A,F10.3)')                                          &
        '              ( w/m2)       SHFL    : ', gpshfl(mzg,nzpa)
    WRITE (yline(k+ 4)(56:99),'(A,F10.3)')                                   &
             '         ( N/m2)       UMFL    : ', gpumfl(mzg,nzpa)
    WRITE (yline(k+ 5),'(A,F10.3)')                                          &
        '                            LHFL    : ', gplhfl(mzg,nzpa)
    WRITE (yline(k+ 5)(56:99),'(A,F10.3)')                                   &
             '         ( N/m2)       VMFL    : ', gpvmfl(mzg,nzpa)
    WRITE (yline(k+ 6),'(A,F10.3)')                                          &
        '              ( g/kg)       QV_S    : ', gpqvs(mzg,nzpa)*1.0E3_wp


    k = 6
    IF (gpfrla(mzg) >= 0.5_wp) THEN
      k = k+1
      WRITE (yline(k),'(A,F10.3)')                                           &
        '              (kg/m2)       RUNOFF_S: ', gpabsf(mzg,nzpa)
      k = k+1
      WRITE (yline(k),'(A,F10.3)')                                           &
        '                            RUNOFF_G: ', gpabgr(mzg,nzpa)
      k = k+1
      WRITE (yline(k),'(A,F10.3)')                                           &
        '  Plants:                   LAI     : ', gplai (mzg,nzpa)
      WRITE (yline(k)(50:99),'(A,F10.3)')                                    &
        '     Ozone:                 VIO3    : ', gpvio3(mzg,nzpa)
      k = k+1
      WRITE (yline(k),'(A,F10.3)')                                           &
        '                            PLCOV   : ', gpplco(mzg,nzpa)
      WRITE (yline(k)(50:99),'(A,F10.3)')                                    &
        '                            HMO3    : ', gphmo3(mzg,nzpa)
      k = k+1
      WRITE (yline(k),'(A,F10.3)')                                           &
        '                            ROOTDP  : ', gproot(mzg,nzpa)

      IF (lmulti_layer) THEN
        ! landpoints with new soil variables

        k = k+2
        WRITE (yline(k),'(A,F10.3)')                                         &
        '  Soil temperatures         T_SNOW  : ', gptsno(mzg,nzpa)-t0_melt
        WRITE (yline(k)(50:99),'(A,F10.3)')                                  &
        '     Soil moistures/Snow    W_SNOW  : ', gpwsno(mzg,nzpa)*1.0E3_wp
        k = k+1
        WRITE (yline(k   ),'(A,F10.3)')                                      &
        '              (dgr C)       T_S     : ', gpt_s (mzg,nzpa)-t0_melt
        WRITE (yline(k)(50:99),'(A,F10.3)')                                  &
        '              (mm H2O)      W_I     : ', gpw_i (mzg,nzpa)*1.0E3_wp
        k = k+1
        WRITE (yline(k   ),'(A,F10.3)')                                      &
        '                            T_G     : ', gpt_g (mzg,nzpa)-t0_melt
        WRITE (yline(k)(50:99),'(A,F10.3)')                                  &
        '                            FRESHSNW: ', gpfrsn(mzg,nzpa)
        k = k+1
        WRITE (yline(k)(50:99),'(A,F10.3)')                                  &
        '              (kg/m3)       RHO_SNOW: ', gprhsn(mzg,nzpa)
        k = k+1
        WRITE (yline(k)(50:99),'(A,F10.3)')                                  &
        '              (m)           H_SNOW  : ', gph_sn(mzg,nzpa)
        k = k+1
        WRITE (yline(k)(50:99),'(A,F10.3)')                                  &
        '              (kg/m2)       SNOW_MELT:', gpsmlt(mzg,nzpa)

        DO kso = 0, ke_soil+1
          k = k+1
          WRITE (yline(k   ),'(A,I2,A,F10.3)')                               &
        '                            T_SO(',kso,'): ', gptso (kso,mzg,nzpa)-t0_melt
          IF (kso > 0) THEN
            WRITE (yline(k)(50:99),'(A,I2,A,F10.3,A,I2,A,F9.3)')             &
        '     W_SO(',kso,'): ', gpwso (kso,mzg,nzpa)*1.0E3_wp,             &
        ' W_SO_ICE(',kso,'): ', gpwice(kso,mzg,nzpa)*1.0E3_wp
          ENDIF
        ENDDO

      ELSE
        ! landpoints with old soil variables
        k = k+2
        WRITE (yline(k),'(A,F10.3)')                                         &
        '  Soil temperatures         T_SNOW  : ', gptsno(mzg,nzpa)-t0_melt
        WRITE (yline(k)(50:99),'(A,F10.3)')                                  &
        '       Soil moistures       W_SNOW  : ', gpwsno(mzg,nzpa)*1.0E3_wp
        k = k+1
        WRITE (yline(k),'(A,F10.3)')                                         &
        '              (dgr C)       T_S     : ', gpt_s (mzg,nzpa)-t0_melt
        WRITE (yline(k)(50:99),'(A,F10.3)')                                  &
        '              (mm H2O)      W_I     : ', gpw_i (mzg,nzpa)*1.0E3_wp
        k = k+1
        WRITE (yline(k),'(A,F10.3)')                                         &
        '                            T_G     : ', gpt_g (mzg,nzpa)-t0_melt
        WRITE (yline(k)(50:99),'(A,F10.3)')                                  &
        '                            W_G1    : ', gpwg1 (mzg,nzpa)*1.0E3_wp
        k = k+1
        WRITE (yline(k),'(A,F10.3)')                                         &
        '                            T_M     : ', gpt_m (mzg,nzpa)-t0_melt
        WRITE (yline(k)(50:99),'(A,F10.3)')                                  &
        '                            W_G2    : ', gpwg2 (mzg,nzpa)*1.0E3_wp
        IF (nlgw==2) THEN
          k = k+1
          WRITE (yline(k),'(A,F10.3)')                                       &
        '                            T_CL    : ', gpt_cl(mzg,nzpa)-t0_melt
          WRITE (yline(k)(50:99),'(A,F10.3)')                                &
        '                            W_CL    : ', gpw_cl(mzg,nzpa)*1.0E3_wp
        ELSE
          k = k+1
          WRITE (yline(k),'(A,F10.3)')                                       &
        '                            T_CL    : ', gpt_cl(mzg,nzpa)-t0_melt
          WRITE (yline(k)(50:99),'(A,F10.3)')                                &
        '                            W_G3    : ', gpwg3 (mzg,nzpa)*1.0E3_wp
          k = k+1
          WRITE (yline(k)(50:99),'(A,F10.3)')                                &
        '                            W_CL    : ', gpw_cl(mzg,nzpa)-t0_melt
        ENDIF
      ENDIF
    ELSE ! (fr_land < 0.5)
      k = k+2
      WRITE (yline(k),'(A,F12.5)')                                           &
        '  Soil temperatures         T_SNOW  : ', gptsno(mzg,nzpa)-t0_melt
      k = k+1
      WRITE (yline(k),'(A,F10.3)')                                           &
        '              (dgr C)       T_S     : ', gpt_s (mzg,nzpa)-t0_melt
      k = k+1
      WRITE (yline(k),'(A,F10.3)')                                           &
        '                            T_G     : ', gpt_g (mzg,nzpa)-t0_melt
    ENDIF

    k = k+2
    WRITE (yline(k),'(A,F10.3)')                                             &
        '  Temperatures              T2M     : ', gpt2m (mzg,nzpa)-t0_melt
    WRITE (yline(k)(50:99),'(A,F10.3)')                                      &
        '       Winds                U10M    : ', gpu10m(mzg,nzpa)
    k = k+1
    WRITE (yline(k),'(A,F10.3)')                                             &
        '              (dgr C)       TD2M    : ', gptd2m(mzg,nzpa)-t0_melt
    WRITE (yline(k)(50:99),'(A,F10.3)')                                      &
        '              ( m/s )       V10M    : ', gpv10m(mzg,nzpa)
    k = k+1
    WRITE (yline(k),'(A,F10.3)')                                             &
        '                            TMIN2M  : ', gptmn2(mzg,nzpa)-t0_melt
    WRITE (yline(k)(50:99),'(A,F10.3)')                                      &
        '                            VBMAX10M: ', gpvb10(mzg,nzpa)
    k = k+1
    WRITE (yline(k),'(A,F10.3)')                                             &
        '                            TMAX2M  : ', gptmx2(mzg,nzpa)-t0_melt

    k = k+2
    WRITE (yline(k),'(A,F10.3)')                                             &
        '  Solar radiation           SOBT    : ', gpsobt(mzg,nzpa)
    WRITE (yline(k)(50:99),'(A,F10.3)')                                      &
        '       Thermal radiation    THBT    : ', gpthbt(mzg,nzpa)
    k = k+1
    WRITE (yline(k),'(A,F10.3)')                                             &
        '              (w/m**2)      SOBS    : ', gpsobs(mzg,nzpa)
    WRITE (yline(k)(50:99),'(A,F10.3)')                                      &
        '              (w/m**2)      THBS    : ', gpthbs(mzg,nzpa)
    k = k+1
    WRITE (yline(k),'(A,F10.3)')                                           &
        '              (w/m**2)      SWDIR_S : ', gpswdir_s(mzg,nzpa)
    k = k+1
    WRITE (yline(k),'(A,F10.3)')                                           &
        '              (w/m**2)      SWDIFD_S: ', gpswdifd_s(mzg,nzpa)
    k = k+1
    WRITE (yline(k),'(A,F10.3)')                                           &
        '              (w/m**2)      SWDIFU_S: ', gpswdifu_s(mzg,nzpa)
    k = k+1
    WRITE (yline(k),'(A,F10.3)')                                             &
        '  Photosynt. active Rad.    PABS    : ', gppabs(mzg,nzpa)
    WRITE (yline(k)(50:99),'(A,F10.3)')                                      &
        '       Surface albedo (%)   ALB     : ', gpalb (mzg,nzpa)*100.0_wp

    k = k+2
    WRITE (yline(k),'(A,A)') '  Precipitation        rates     and    amount', &
        '                Cloud Cover'
    k = k+1
    WRITE (yline(k),'(A,A)') '                       (mm/d)            (mm) ', &
        '                    (%)'

    k = k+1
    WRITE (yline(k),'(A,F10.3,A,F10.3)')                           &
        '        RAIN_GSP: ', gpprrs(mzg,nzpa)*day_len, '       ', gprrsn(mzg,nzpa)
    WRITE (yline(k)(50:99),'(A,F10.3)')                                      &
        '             CLCH    : ', gpclch(mzg,nzpa)*100.0_wp
    k = k+1
    WRITE (yline(k),'(A,F10.3,A,F10.3)')                           &
        '        SNOW_GSP: ', gpprss(mzg,nzpa)*day_len, '       ', gprssn(mzg,nzpa)
    WRITE (yline(k)(50:99),'(A,F10.3)')                                      &
        '             CLCM    : ', gpclcm(mzg,nzpa)*100.0_wp
    k = k+1
    WRITE (yline(k),'(A,F10.3,A,F10.3)')                           &
        '        RAIN_CON: ', gpprrc(mzg,nzpa)*day_len, '       ', gprrkn(mzg,nzpa)
    WRITE (yline(k)(50:99),'(A,F10.3)')                                      &
        '             CLCL    : ', gpclcl(mzg,nzpa)*100.0_wp
    k = k+1
    WRITE (yline(k),'(A,F10.3,A,F10.3)')                           &
        '        SNOW_CON: ', gpprsc(mzg,nzpa)*day_len, '       ', gprskn(mzg,nzpa)
    WRITE (yline(k)(50:99),'(A,F10.3)')                                      &
        '             CLCT    : ', gpclct(mzg,nzpa)*100.0_wp

    zpr = gpprrs(mzg,nzpa) + gpprss(mzg,nzpa) + gpprrc(mzg,nzpa) + gpprsc(mzg,nzpa)
    zrr = gprrsn(mzg,nzpa) + gprssn(mzg,nzpa) + gprrkn(mzg,nzpa) + gprskn(mzg,nzpa)
    k = k+1
    WRITE (yline(k),'(A,F10.3,A,F10.3)')                           &
        '        TOTAL:    ', zpr*day_len, '       ', zrr

    k1 = k

    ! Output to the file
    DO k = 1 , k1
      WRITE (nugpunit(mzg), '(A,A)') ' ',yline(k)
    ENDDO
    WRITE (nugpunit(mzg), '(A)')  '  '
    WRITE (nugpunit(mzg), '(A)')  '  '

    CLOSE (nugpunit(mzg), STATUS='KEEP')
  ENDDO
ENDIF

!------------------------------------------------------------------------------
! Section 3:  Write variables for special SCLM gridpoint output
!------------------------------------------------------------------------------

IF (lgpspec) THEN

  form0='(f'
  form1='ss,es14.7)'
  form2=form1//')'
  form3='('//form1

  lpreslev=.FALSE. !use height based model levels
  lhydpres=.FALSE. !use non hydrostatic model pressure

  ! Loop over all chosen gridpoints
  DO  mzg = 1, ngp
    ! open the file
    OPEN (nugpunit(mzg), FILE=TRIM(yugpname(mzg)), FORM='FORMATTED',      &
                         STATUS='UNKNOWN', POSITION='APPEND', IOSTAT=nstat)
    IF (nstat /= 0) THEN
      yerrmsg = 'opening of file '//TRIM(yugpname(mzg))//' failed'
      CALL model_abort (my_cart_id, 7004, yerrmsg, yroutine)
    ENDIF

    !------------------------------------------------------------------------
    !Data that are constant in time:
    !------------------------------------------------------------------------

    IF (nzpa.EQ.1) THEN

      !Model level information
      !------------------------------------------------------------------------

      IF (lpreslev) THEN
         WRITE (nugpunit(mzg),'(a)') &
               'Reference pressure of atmospheric half levels [hPa] (including the surface):'
         WRITE (nugpunit(mzg),*)

         flat=.FALSE. !transition level not yet found
         zpb=gpps(mzg, 1); WRITE (nugpunit(mzg),form3) zpb*0.01_wp !surface pressure
         DO k=ke,1,-1
            IF (lhydpres) THEN !use hydrostatic pressure on boundary levels (half levels)
               ztv=gpt(k,mzg,1)*(1+rvd_m_o*gpqv(k,mzg,1)) !virtual temperature
               zpb=zpb*EXP((gphhl(k+1,mzg)-gphhl(k,mzg))*g/r_d/ztv) !upper bd. pressure
            ELSE !use non hydrostatic model pressure
               zpb=2.0_wp*(gpp0(mzg,k)+gppp(k,mzg, 1))-zpb !upper bd. pressure
            END IF
            WRITE (nugpunit(mzg),form3, ADVANCE='NO') zpb*0.01_wp

            IF ((k > vcoord%kflat) .OR. flat) THEN
               !Terrain following levels or tranition already found:
               WRITE (nugpunit(mzg),*) !next line
            ELSE !transition level
               flat=.TRUE. !transition level found
               WRITE (nugpunit(mzg),'(a)') ' *'
            END IF
         END DO
      ELSE
         WRITE (nugpunit(mzg),'(a)') &
               'Fixed height of atmospheric half levels [m] (above the surface):'
         WRITE (nugpunit(mzg),*)

         DO k=ke,1,-1
            WRITE (nugpunit(mzg),form3) gphhl(k,mzg)-gphsur(mzg)
         END DO
      END IF

      WRITE (nugpunit(mzg),*)

      WRITE (nugpunit(mzg),'(a)') &
            'Fixed height of soil half levels [m] (except the surface):'
      WRITE (nugpunit(mzg),*)

      zhb=0.0_wp !upper boundary height of the soil
      DO k=1, ke_soil+1
         zhb=2.0_wp*czmls(k)-zhb !boundary height of soil layer
         WRITE (nugpunit(mzg),form3) -zhb 
      END DO
      WRITE (nugpunit(mzg),*)

      !External parameters:
      !------------------------------------------------------------------------

      WRITE (nugpunit(mzg),'("G_Lat[Deg]     ",'//form1) gplati(mzg)*raddeg
      WRITE (nugpunit(mzg),'("G_Lon[Deg]     ",'//form1) gplong(mzg)*raddeg
      WRITE (nugpunit(mzg),'("HS_NN[m]       ",'//form1) gphsur(mzg) 
      WRITE (nugpunit(mzg),'("FR_Land[1]     ",'//form1) gpfrla(mzg) 
      WRITE (nugpunit(mzg),'("SO_Type[]      ",'//form1) gpsoil(mzg) 
      WRITE (nugpunit(mzg),*)

    END IF   

    !------------------------------------------------------------------------
    !Data that may change in time
    !------------------------------------------------------------------------

    i_hor=INT(zth2) !INTEGER value of the hour of the actual day
    d_hor=zth2-REAL(i_hor,wp) !fractional part of the hour
    prec=PRECISION(d_hor) !precision of numbers

    WRITE (form0(3:4),'(i2.2)') prec+1
    form0(5:5)='.'
    WRITE (form0(6:7),'(i2.2)') prec
    form0(8:8)=')'

    DO i=1,3+15*7
       WRITE (nugpunit(mzg),'(a)',ADVANCE='NO') '-'
    END DO
    WRITE (nugpunit(mzg),*)
        
    WRITE (dezi,form0) d_hor
    dezi=dezi(INDEX(dezi,'.'):)
    nlen=LEN_TRIM(dezi)
    DO WHILE (dezi(nlen-1:nlen).EQ.'00')
       nlen=nlen-1
    END DO

    WRITE (nugpunit(mzg),'("Date           ",a,a)') ydat21,dezi(1:nlen)
    WRITE (nugpunit(mzg),*)

    ! Surface variables:

    WRITE (nugpunit(mzg),'("PL_Cov[1]      ",'//form1) gpplco(mzg,nzpa)
    WRITE (nugpunit(mzg),'("LA_Ind[1]      ",'//form1) gplai (mzg,nzpa)
    WRITE (nugpunit(mzg),'("RO_Dept[m]     ",'//form1) gproot(mzg,nzpa)
    WRITE (nugpunit(mzg),'("O3_Int[Pa]     ",'//form1) gpvio3(mzg,nzpa)
    WRITE (nugpunit(mzg),'("O3_Max[Pa]     ",'//form1) gphmo3(mzg,nzpa)
    WRITE (nugpunit(mzg),'("RL_Mom[m]      ",'//form1) gpgz0(mzg,nzpa)/g

    IF (.NOT.lmulti_layer) THEN
       IF ( nlgw == 3 ) THEN !three hydrol. layers
          dzw(1)=cdzw13; dzw(2)=cdzw23; dzw(3)=cdzw33; dzw(4)=cdzw33
       ELSE
          dzw(1)=cdzw12; dzw(2)=cdzw22; dzw(3)=cdzw22
       END IF

       WRITE (nugpunit(mzg),'("TE_Clim[C]     ",'//form1) gpt_cl(mzg,nzpa)-t0_melt
       WRITE (nugpunit(mzg),'("WE_Clim[SVF]   ",'//form1) gpw_cl(mzg,nzpa)/dzw(nlgw+1)
    END IF   

    WRITE (nugpunit(mzg),'("PS_Mean[hPa]   ",'//form1) gpps(mzg,nzpa)*0.01_wp
    WRITE (nugpunit(mzg),'("TS_Erth[C]     ",'//form1) gpt_s(mzg,nzpa)-t0_melt
    WRITE (nugpunit(mzg),'("TS_Snow[C]     ",'//form1) gptsno(mzg,nzpa)-t0_melt
    WRITE (nugpunit(mzg),'("QVS_Mean[1]    ",'//form1) gpqvs(mzg,nzpa)
    WRITE (nugpunit(mzg),'("WS_Intc[m]     ",'//form1) gpw_i(mzg,nzpa)
    WRITE (nugpunit(mzg),'("WS_Snow[m]     ",'//form1) gpwsno(mzg,nzpa)

    IF (.NOT.lmulti_layer) THEN
       WRITE (nugpunit(mzg),'("TE_Mean[C]     ",'//form1) gpt_m(mzg,nzpa)-t0_melt
       WRITE (nugpunit(mzg),'("WE_1[SVF]      ",'//form1) gpwg1(mzg,nzpa)/dzw(1)
       WRITE (nugpunit(mzg),'("WE_2[SVF]      ",'//form1) gpwg2(mzg,nzpa)/dzw(2)
       IF ( nlgw == 3 ) THEN
          WRITE (nugpunit(mzg),'("WE_3[SVF]      ",'//form1) gpwg3(mzg,nzpa)/dzw(3)
       END IF
    END IF

  ! write(nugpunit(mzg),'("Regen           ",'//form1) gpprrs(mzg,1)+gpprrc(mzg,1)
  ! write(nugpunit(mzg),'("Schnee          ",'//form1) gpprss(mzg,1)+gpprsc(mzg,1)

  ! write(nugpunit(mzg),'("Rdrg            ",'//form1) gprdrg(mzg)
  ! write(nugpunit(mzg),'("Cbig            ",'//form1) gpcbig(mzg)
  ! write(nugpunit(mzg),'("Csml            ",'//form1) gpcsml(mzg)
  ! write(nugpunit(mzg),'("Hdrg            ",'//form1) gphcan(mzg)
  ! write(nugpunit(mzg),'("Dpat            ",'//form1) gpdpat(mzg)
  ! write(nugpunit(mzg),'("Ccan            ",'//form1) gpsai(mzg)

    WRITE (nugpunit(mzg),*)

    ! Soil profiles

    IF (lmulti_layer) THEN
       WRITE (nugpunit(mzg),'(a)') 'k   HM_Erth[mAG]   TO_Erth[C]     '// &
                                  'WC_Erth[SVF]   IC_Erth[SVF]'
       WRITE (nugpunit(mzg),*)

       zhb=0.0_wp !soil depth at the surface
       DO k=1,ke_soil+1
          dzh=zhb !upper boundary level height
          zhb=2.0_wp*czmls(k)-zhb !lower boundary level height
          dzh=zhb-dzh !depth of the soil layer
          WRITE (nugpunit(mzg),'(i3,4(x,'//form2) k, -czmls(k), &
                gptso(k,mzg,nzpa)-t0_melt, gpwso(k,mzg,nzpa)/dzh, gpwice(k,mzg,nzpa)/dzh
       END DO
       WRITE (nugpunit(mzg),*)
    END IF   

    ! Atmospheric profiles

    IF (lpreslev) THEN
       WRITE (nugpunit(mzg),'(a)',ADVANCE='NO') 'k   '//           &
            'PM_Air[hPa]    HM_Air[mAG]    VZ_Air[m/s]    VM_Air[m/s]    '// &
            'TO_Air[C]      QV_Air[1]      QL_Air[1]      '
    ELSE
       WRITE (nugpunit(mzg),'(a)',ADVANCE='NO') 'k   '//           &
            'HM_Air[mAG]    PM_Air[hPa]    VZ_Air[m/s]    VM_Air[m/s]    '// &
            'TO_Air[C]      QV_Air[1]      QL_Air[1]      '
    END IF

    IF ( ASSOCIATED(qi) ) THEN
       WRITE (nugpunit(mzg),'(a)') 'QI_Air[1]'
    ELSE
       WRITE (nugpunit(mzg),*) !next line
    END IF
    WRITE (nugpunit(mzg),*) 

    zpb=gpps(mzg,nzpa) !surface pressure
    DO k=ke,0,-1

       IF (k.EQ.0) THEN !additional measurement until the top of the model domain
          l=k+1 !use the mid level values of the upper most level

          zpm=zpb                      !use the upper most bd. press.
          zhm=gphhl(l,mzg)-gphsur(mzg) !use the upper most bd. height
       ELSE
          l=k !use mid level values of the current level

          IF (lhydpres) THEN !use hydrostatic pressure on mid levels (main levels)
             zpm=zpb !lower hydrostatic boundary pressure
             ztv=gpt(l,mzg,nzpa)*(1+rvd_m_o*gpqv(l,mzg,nzpa)) !virtual temperature
             zpb=zpm*EXP((gphhl(l+1,mzg)-gphhl(l,mzg))*g/r_d/ztv) !upper hydr. bd. pres.
             zpm=(zpm+zpb)*0.5_wp !hydrostatic mid level pressure
          ELSE !use non hydrostatic model pressure
             zpm=gpp0(mzg,l)+gppp(l,mzg,nzpa) !non hydr. mid level pressure
             zpb=2.0_wp*zpm-zpb !non hydr. upper boundary pressure
          END IF
          zhm=gphml(l,mzg)-gphsur(mzg)
       END IF 

       IF (lpreslev) THEN
          WRITE (nugpunit(mzg),'(i3,7(x,'//form2,ADVANCE='NO') &
                k, zpm*0.01_wp, zhm, &
                gpu(l,mzg,nzpa) , gpv(l,mzg,nzpa), &
                gpt(l,mzg,nzpa)-t0_melt, gpqv(l,mzg,nzpa), gpqc(l,mzg,nzpa)
       ELSE
          WRITE (nugpunit(mzg),'(i3,7(x,'//form2,ADVANCE='NO') &
                k, zhm, zpm*0.01_wp, &
                gpu(l,mzg,nzpa) , gpv(l,mzg,nzpa), &
                gpt(l,mzg,nzpa)-t0_melt, gpqv(l,mzg,nzpa), gpqc(l,mzg,nzpa)
       END IF

       IF ( ASSOCIATED(qi) ) THEN
          WRITE (nugpunit(mzg),'(x,'//form1) gpqi(l,mzg,nzpa)
       ELSE
          WRITE (nugpunit(mzg),*) !next line
       END IF

    END DO   
    WRITE (nugpunit(mzg),*)

  CLOSE (nugpunit(mzg), STATUS='KEEP')    
  END DO
   
END IF

!------------------------------------------------------------------------------
! Section 4:  Store variables for short gridpoint output
!------------------------------------------------------------------------------

IF (lgpshort) THEN

  ! Loop over all chosen gridpoints
  DO mzg = 1,ngp

    izgp = stationlist(mzg)%igp
    jzgp = stationlist(mzg)%jgp

    ! k-index of the level 500m above surface
    DO k = 1,ke
      IF ( 0.5_wp*(hhl(izgp,jzgp,k)+hhl(izgp,jzgp,k+1))               &
                                  -hsurf(izgp,jzgp)  >= 500.0_wp ) THEN
        k500m = k
      ENDIF
    ENDDO

    gpclcl(mzg,nzpa) = clcl(izgp,jzgp)
    gpclcm(mzg,nzpa) = clcm(izgp,jzgp)
    gpclch(mzg,nzpa) = clch(izgp,jzgp)
    gpfog (mzg,nzpa) = clc_sgs(izgp,jzgp,ke)
    gpwsno(mzg,nzpa) = w_snow(izgp,jzgp,ntl)
    gpsmlt(mzg,nzpa) = snow_melt(izgp,jzgp)

    ! variables from routine near_surface 
    gpps  (mzg,nzpa) = ps    (izgp,jzgp, ntl)
    gpu10m(mzg,nzpa) = u_10m (izgp,jzgp)
    gpv10m(mzg,nzpa) = v_10m (izgp,jzgp)
    gpu950(mzg,nzpa) = 0.5_wp*(u(izgp,jzgp,k500m,ntl) + u(izgp-1,jzgp,k500m,ntl))
    gpv950(mzg,nzpa) = 0.5_wp*(v(izgp,jzgp,k500m,ntl) + v(izgp,jzgp-1,k500m,ntl))
    gptp  (mzg,nzpa) = t     (izgp,jzgp,ke,ntl)
    gpt2m (mzg,nzpa) = t_2m  (izgp,jzgp)
    gpt_g (mzg,nzpa) = t_g   (izgp,jzgp, ntl)
    gptd2m(mzg,nzpa) = td_2m (izgp,jzgp)
    gprrsn(mzg,nzpa) = rain_gsp(izgp,jzgp)
    gprssn(mzg,nzpa) = snow_gsp(izgp,jzgp)
    gprrkn(mzg,nzpa) = rain_con(izgp,jzgp)
    gprskn(mzg,nzpa) = snow_con(izgp,jzgp)

    ! wind components and temperatures in 850, 700 and 500 hpa
    k850 = ke + 1
    k700 = ke + 1
    k500 = ke + 1
    DO k = ke, 1, -1
      zphf = p0(izgp,jzgp,k) + pp(izgp,jzgp,k,ntl)
      IF (zphf*0.01_wp > 850.0_wp) k850 = k
      IF (zphf*0.01_wp > 700.0_wp) k700 = k
      IF (zphf*0.01_wp > 500.0_wp) k500 = k
    ENDDO
    IF (k850 < ke+1) THEN
      zupu = 0.5_wp * (u(izgp,jzgp,k850,ntl) + u(izgp-1,jzgp,k850,ntl))
      zvpu = 0.5_wp * (v(izgp,jzgp,k850,ntl) + v(izgp,jzgp-1,k850,ntl))
      ztpu = t(izgp,jzgp, k850, ntl)
      zpu  = p0(izgp,jzgp,k850) + pp(izgp,jzgp,k850,ntl)
      zupo = 0.5_wp * (u(izgp,jzgp,k850-1,ntl) + u(izgp-1,jzgp,k850-1,ntl))
      zvpo = 0.5_wp * (v(izgp,jzgp,k850-1,ntl) + v(izgp,jzgp-1,k850-1,ntl))
      ztpo = t(izgp,jzgp, k850-1, ntl)
      zpo  = p0(izgp,jzgp,k850-1) + pp(izgp,jzgp,k850-1,ntl)
      zpx  = 85000.0_wp
      gpu850 (mzg,nzpa) = zupu + (zupu-zupo)/(zpu-zpo)*(zpx-zpu)
      gpv850 (mzg,nzpa) = zvpu + (zvpu-zvpo)/(zpu-zpo)*(zpx-zpu)
      gpt850 (mzg,nzpa) = ztpu + (ztpu-ztpo)/(zpu-zpo)*(zpx-zpu)
    ELSE IF (k850 >= ke+1) THEN
      gpu850 (mzg,nzpa) = gpu950 (mzg,nzpa)
      gpv850 (mzg,nzpa) = gpv950 (mzg,nzpa)
      gpt850 (mzg,nzpa) = gptp   (mzg,nzpa)
    ENDIF
    IF (k700 < ke+1) THEN
      zupu = 0.5_wp * (u(izgp,jzgp,k700,ntl) + u(izgp-1,jzgp,k700,ntl))
      zvpu = 0.5_wp * (v(izgp,jzgp,k700,ntl) + v(izgp,jzgp-1,k700,ntl))
      ztpu = t(izgp,jzgp, k700, ntl)
      zpu  = p0(izgp,jzgp,k700) + pp(izgp,jzgp,k700,ntl)
      zupo = 0.5_wp * (u(izgp,jzgp,k700-1,ntl) + u(izgp-1,jzgp,k700-1,ntl))
      zvpo = 0.5_wp * (v(izgp,jzgp,k700-1,ntl) + v(izgp,jzgp-1,k700-1,ntl))
      ztpo = t(izgp,jzgp, k700-1, ntl)
      zpo  = p0(izgp,jzgp,k700-1) + pp(izgp,jzgp,k700-1,ntl)
      zpx  = 70000.0_wp
      gpu700 (mzg,nzpa) = zupu + (zupu-zupo)/(zpu-zpo)*(zpx-zpu)
      gpv700 (mzg,nzpa) = zvpu + (zvpu-zvpo)/(zpu-zpo)*(zpx-zpu)
      gpt700 (mzg,nzpa) = ztpu + (ztpu-ztpo)/(zpu-zpo)*(zpx-zpu)
    ELSE IF (k700 >= ke+1) THEN
      gpu700 (mzg,nzpa) = gpu950 (mzg,nzpa)
      gpv700 (mzg,nzpa) = gpv950 (mzg,nzpa)
      gpt700 (mzg,nzpa) = gptp   (mzg,nzpa)
    ENDIF
    IF (k500 < ke+1) THEN
      zupu = 0.5_wp * (u(izgp,jzgp,k500,ntl) + u(izgp-1,jzgp,k500,ntl))
      zvpu = 0.5_wp * (v(izgp,jzgp,k500,ntl) + v(izgp,jzgp-1,k500,ntl))
      ztpu = t(izgp,jzgp, k500, ntl)
      zpu  = p0(izgp,jzgp,k500) + pp(izgp,jzgp,k500,ntl)
      zupo = 0.5_wp * (u(izgp,jzgp,k500-1,ntl) + u(izgp-1,jzgp,k500-1,ntl))
      zvpo = 0.5_wp * (v(izgp,jzgp,k500-1,ntl) + v(izgp,jzgp-1,k500-1,ntl))
      ztpo = t(izgp,jzgp, k500-1, ntl)
      zpo  = p0(izgp,jzgp,k500-1) + pp(izgp,jzgp,k500-1,ntl)
      zpx  = 50000.0_wp
      gpu500 (mzg,nzpa) = zupu + (zupu-zupo)/(zpu-zpo)*(zpx-zpu)
      gpv500 (mzg,nzpa) = zvpu + (zvpu-zvpo)/(zpu-zpo)*(zpx-zpu)
      gpt500 (mzg,nzpa) = ztpu + (ztpu-ztpo)/(zpu-zpo)*(zpx-zpu)
    ELSE IF (k500 >= ke+1) THEN
      gpu500 (mzg,nzpa) = gpu950 (mzg,nzpa)
      gpv500 (mzg,nzpa) = gpv950 (mzg,nzpa)
      gpt500 (mzg,nzpa) = gptp   (mzg,nzpa)
    ENDIF

    ! base and top height of convective cloud above msl
    kbas = NINT( bas_con(izgp,jzgp) )
    IF(kbas > 0) THEN
      gphbas(mzg,nzpa) = hhl(izgp,jzgp,kbas)
    ELSE
      gphbas(mzg,nzpa) = 0.0_wp
    ENDIF

    ktop = NINT( top_con(izgp,jzgp) )
    IF(ktop > 0) THEN
      gphtop(mzg,nzpa) = 0.5_wp*(hhl(izgp,jzgp,ktop)+hhl(izgp,jzgp,ktop+1))
    ELSE
      gphtop(mzg,nzpa) = 0.0_wp
    ENDIF

  ENDDO
ENDIF

!------------------------------------------------------------------------------
! Section 5:  Write variables for short gridpoint output
!------------------------------------------------------------------------------

IF (lgpshort) THEN
  ! Loop over all chosen gridpoints
  DO  mzg = 1, ngp
    ! open the file
    OPEN (nugpunit(mzg), FILE=TRIM(yugpname(mzg)), FORM='FORMATTED',      &
                         STATUS='UNKNOWN', POSITION='APPEND', IOSTAT=nstat)
    IF (nstat /= 0) THEN
      yerrmsg = 'opening of file '//TRIM(yugpname(mzg))//' failed'
      CALL model_abort (my_cart_id, 7004, yerrmsg, yroutine)
    ENDIF

    ! cloud cover in 1/8
    izch    = INT( (gpclch(mzg,nzpa) + 0.001_wp) * 800.0_wp ) / 100
    izcm    = INT( (gpclcm(mzg,nzpa) + 0.001_wp) * 800.0_wp ) / 100
    izcl    = INT( (gpclcl(mzg,nzpa) + 0.001_wp) * 800.0_wp ) / 100
    izfog   = INT( (gpfog (mzg,nzpa) + 0.001_wp) * 800.0_wp ) / 100
  
    ! surface pressure reduced to sealevel
    izpsr = NINT (gpps(mzg,nzpa) *                                        &
                   EXP (0.03418_wp * gphsur(mzg) / (gpt2m(mzg,nzpa))) )
  
    ! wind components in W-E and N-S components for the geographical grid
    CALL uvrot2uv (gpu10m(mzg,nzpa), gpv10m(mzg,nzpa), gplati(mzg)*raddeg,    &
                   gplong(mzg)*raddeg, pollat, pollon, zu10m, zv10m )
 
    CALL uvrot2uv (gpu950(mzg,nzpa), gpv950(mzg,nzpa), gplati(mzg)*raddeg,    &
                   gplong(mzg)*raddeg, pollat, pollon, zu950, zv950 )
 
    CALL uvrot2uv (gpu850(mzg,nzpa), gpv850(mzg,nzpa), gplati(mzg)*raddeg,    &
                   gplong(mzg)*raddeg, pollat, pollon, zu850, zv850 )
 
    CALL uvrot2uv (gpu700(mzg,nzpa), gpv700(mzg,nzpa), gplati(mzg)*raddeg,    &
                   gplong(mzg)*raddeg, pollat, pollon, zu700, zv700 )
 
    CALL uvrot2uv (gpu500(mzg,nzpa), gpv500(mzg,nzpa), gplati(mzg)*raddeg,    &
                   gplong(mzg)*raddeg, pollat, pollon, zu500, zv500 )
 
    ! wind components in direction and speed
    IF (zv10m == 0.0_wp) zv10m = 1.0E-10_wp
    IF (zv950 == 0.0_wp) zv950 = 1.0E-10_wp
    IF (zv850 == 0.0_wp) zv850 = 1.0E-10_wp
    IF (zv700 == 0.0_wp) zv700 = 1.0E-10_wp
    IF (zv500 == 0.0_wp) zv500 = 1.0E-10_wp
 
    izd10m = NINT( raddeg * ATAN2 (zu10m, zv10m) + 180.0_wp )
    IF(izd10m == 360) izd10m = 0
 
    izd950 = NINT( raddeg * ATAN2 (zu950, zv950) + 180.0_wp )
    IF(izd950 == 360) izd950 = 0
 
    izd850 = NINT( raddeg * ATAN2 (zu850, zv850) + 180.0_wp )
    IF(izd850 == 360) izd850 = 0
 
    izd700 = NINT( raddeg * ATAN2 (zu700, zv700) + 180.0_wp )
    IF(izd700 == 360) izd700 = 0
 
    izd500 = NINT( raddeg * ATAN2 (zu500, zv500) + 180.0_wp )
    IF(izd500 == 360) izd500 = 0
 
    ! wind speed in m/s
    izf10m = NINT ( SQRT (zu10m**2 + zv10m**2) )
    izf950 = NINT ( SQRT (zu950**2 + zv950**2) )
    izf850 = NINT ( SQRT (zu850**2 + zv850**2) )
    izf700 = NINT ( SQRT (zu700**2 + zv700**2) )
    izf500 = NINT ( SQRT (zu500**2 + zv500**2) )
  
    ! Output of the values for one grid point and one time step
    WRITE (nugpunit(mzg),'(F6.1, F8.2, 5(I4,A,I3), 7F6.1, 1X, 3I1, I2,  &
  &                       2I5,2F6.2,F7.3,F7.3)')                        &
      hh2, izpsr*0.01_wp, izd10m,'/', izf10m, izd950,'/', izf950,   &
      izd850,'/', izf850, izd700,'/', izf700, izd500,'/', izf500,       &
      gpt_g (mzg,nzpa) - t0_melt, gpt2m (mzg,nzpa) - t0_melt,           &
      gptd2m(mzg,nzpa) - t0_melt, gptp  (mzg,nzpa) - t0_melt,           &
      gpt850(mzg,nzpa) - t0_melt, gpt700(mzg,nzpa) - t0_melt,           &
      gpt500(mzg,nzpa) - t0_melt, izch, izcm, izcl, izfog,              &
      NINT (gphbas(mzg,nzpa)*0.1_wp),                               &
      NINT (gphtop(mzg,nzpa)*0.1_wp),                               &
      gprrsn(mzg,nzpa) + gprrkn(mzg,nzpa),                              &
      gprssn(mzg,nzpa) + gprskn(mzg,nzpa),                              &
      gpwsno(mzg,nzpa)
      !for writing snow_melt     gpwsno(mzg,nzpa),  gpsmlt(mzg,nzpa)

    CLOSE (nugpunit(mzg), STATUS='KEEP')
  ENDDO
ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE gridpoints

!==============================================================================
!==============================================================================
!+ Module Routine to allocate fields for grid point output 
!------------------------------------------------------------------------------

SUBROUTINE alloc_gridpoints   (istat)

!------------------------------------------------------------------------------
!
! Description:
!   This routine allocates space for all fields necessary for grid point output
!   and initializes them with 0.
!
! Method:
!   All ALLOCATABLE fields are allocated with the ALLOCATE statement.
!   The status of the allocation is checked. In case of an error, the
!   error variable is set:  istat = 1
!
!------------------------------------------------------------------------------

! Parameters
INTEGER (KIND=iintegers), INTENT(OUT)   ::       &
  istat              ! for local error-code

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE alloc_gridpoints
!------------------------------------------------------------------------------

  istat = 0

  ! Arrays for the short and the long form of grid point output
  ! -----------------------------------------------------------

  ALLOCATE ( gpps     (ngp, nstepsgp) , STAT=istat )  ; gpps   = 0.0_wp
  ALLOCATE ( gppp  (ke,ngp, nstepsgp) , STAT=istat )  ; gppp   = 0.0_wp
  ALLOCATE ( gpt2m    (ngp, nstepsgp) , STAT=istat )  ; gpt2m  = 0.0_wp
  ALLOCATE ( gptd2m   (ngp, nstepsgp) , STAT=istat )  ; gptd2m = 0.0_wp
  ALLOCATE ( gpu10m   (ngp, nstepsgp) , STAT=istat )  ; gpu10m = 0.0_wp
  ALLOCATE ( gpv10m   (ngp, nstepsgp) , STAT=istat )  ; gpv10m = 0.0_wp
  ALLOCATE ( gpt_g    (ngp, nstepsgp) , STAT=istat )  ; gpt_g  = 0.0_wp
  ALLOCATE ( gprrsn   (ngp, nstepsgp) , STAT=istat )  ; gprrsn = 0.0_wp
  ALLOCATE ( gprssn   (ngp, nstepsgp) , STAT=istat )  ; gprssn = 0.0_wp
  ALLOCATE ( gprrkn   (ngp, nstepsgp) , STAT=istat )  ; gprrkn = 0.0_wp
  ALLOCATE ( gprskn   (ngp, nstepsgp) , STAT=istat )  ; gprskn = 0.0_wp
  ALLOCATE ( gpwsno   (ngp, nstepsgp) , STAT=istat )  ; gpwsno = 0.0_wp
  ALLOCATE ( gpclcl   (ngp, nstepsgp) , STAT=istat )  ; gpclcl = 0.0_wp
  ALLOCATE ( gpclcm   (ngp, nstepsgp) , STAT=istat )  ; gpclcm = 0.0_wp
  ALLOCATE ( gpclch   (ngp, nstepsgp) , STAT=istat )  ; gpclch = 0.0_wp
  ALLOCATE ( gpclct   (ngp, nstepsgp) , STAT=istat )  ; gpclct = 0.0_wp
  ALLOCATE ( gpsmlt   (ngp, nstepsgp) , STAT=istat )  ; gpsmlt = 0.0_wp


  ALLOCATE ( gpp0     (ngp,ke)        , STAT=istat )  ; gpp0   = 0.0_wp
  ALLOCATE ( gphsur   (ngp)           , STAT=istat )  ; gphsur = 0.0_wp
  ALLOCATE ( gpfrla   (ngp)           , STAT=istat )  ; gpfrla = 0.0_wp
  ALLOCATE ( gpsoil   (ngp)           , STAT=istat )  ; gpsoil = 0.0_wp
  ALLOCATE ( gplati   (ngp)           , STAT=istat )  ; gplati = 0.0_wp
  ALLOCATE ( gplong   (ngp)           , STAT=istat )  ; gplong = 0.0_wp

  IF (lgpshort) THEN

    ! Arrays only for the short form of grid point output
    ! ---------------------------------------------------

    ALLOCATE ( gpu950   (ngp, nstepsgp) , STAT=istat )  ;  gpu950 = 0.0_wp
    ALLOCATE ( gpv950   (ngp, nstepsgp) , STAT=istat )  ;  gpv950 = 0.0_wp
    ALLOCATE ( gpu850   (ngp, nstepsgp) , STAT=istat )  ;  gpu850 = 0.0_wp
    ALLOCATE ( gpv850   (ngp, nstepsgp) , STAT=istat )  ;  gpv850 = 0.0_wp
    ALLOCATE ( gpu700   (ngp, nstepsgp) , STAT=istat )  ;  gpu700 = 0.0_wp
    ALLOCATE ( gpv700   (ngp, nstepsgp) , STAT=istat )  ;  gpv700 = 0.0_wp
    ALLOCATE ( gpu500   (ngp, nstepsgp) , STAT=istat )  ;  gpu500 = 0.0_wp
    ALLOCATE ( gpv500   (ngp, nstepsgp) , STAT=istat )  ;  gpv500 = 0.0_wp
    ALLOCATE ( gptp     (ngp, nstepsgp) , STAT=istat )  ;  gptp   = 0.0_wp
    ALLOCATE ( gpt850   (ngp, nstepsgp) , STAT=istat )  ;  gpt850 = 0.0_wp
    ALLOCATE ( gpt700   (ngp, nstepsgp) , STAT=istat )  ;  gpt700 = 0.0_wp
    ALLOCATE ( gpt500   (ngp, nstepsgp) , STAT=istat )  ;  gpt500 = 0.0_wp
    ALLOCATE ( gpfog    (ngp, nstepsgp) , STAT=istat )  ;  gpfog  = 0.0_wp
    ALLOCATE ( gphbas   (ngp, nstepsgp) , STAT=istat )  ;  gphbas = 0.0_wp
    ALLOCATE ( gphtop   (ngp, nstepsgp) , STAT=istat )  ;  gphtop = 0.0_wp

  ELSE IF (lgplong .OR. lgpspec) THEN

    ! Arrays only for the long form of grid point output
    ! --------------------------------------------------

    ALLOCATE ( gpu   (ke,  ngp, nstepsgp) , STAT=istat );  gpu    = 0.0_wp
    ALLOCATE ( gpv   (ke,  ngp, nstepsgp) , STAT=istat );  gpv    = 0.0_wp
    ALLOCATE ( gpw   (ke1, ngp, nstepsgp) , STAT=istat );  gpw    = 0.0_wp
    ALLOCATE ( gpt   (ke,  ngp, nstepsgp) , STAT=istat );  gpt    = 0.0_wp
    ALLOCATE ( gpqv  (ke,  ngp, nstepsgp) , STAT=istat );  gpqv   = 0.0_wp
    ALLOCATE ( gpqc  (ke,  ngp, nstepsgp) , STAT=istat );  gpqc   = 0.0_wp
    ALLOCATE ( gpqi  (ke,  ngp, nstepsgp) , STAT=istat );  gpqi   = 0.0_wp
    ALLOCATE ( gphhl (ke1, ngp)           , STAT=istat );  gphhl  = 0.0_wp
    ALLOCATE ( gphml (ke,  ngp)           , STAT=istat );  gphml  = 0.0_wp
    ALLOCATE ( gptkvm(ke,  ngp, nstepsgp) , STAT=istat );  gptkvm = 0.0_wp
    ALLOCATE ( gptkvh(ke,  ngp, nstepsgp) , STAT=istat );  gptkvh = 0.0_wp
    ALLOCATE ( gpclcs(ke,  ngp, nstepsgp) , STAT=istat );  gpclcs = 0.0_wp
    ALLOCATE ( gpclcc(ke,  ngp, nstepsgp) , STAT=istat );  gpclcc = 0.0_wp
    ALLOCATE ( gptsno     (ngp, nstepsgp) , STAT=istat );  gptsno = 0.0_wp
    ALLOCATE ( gpt_s      (ngp, nstepsgp) , STAT=istat );  gpt_s  = 0.0_wp
    ALLOCATE ( gpt_m      (ngp, nstepsgp) , STAT=istat );  gpt_m  = 0.0_wp
    ALLOCATE ( gpqvs      (ngp, nstepsgp) , STAT=istat );  gpqvs  = 0.0_wp
    ALLOCATE ( gpw_i      (ngp, nstepsgp) , STAT=istat );  gpw_i  = 0.0_wp
    ALLOCATE ( gpwg1      (ngp, nstepsgp) , STAT=istat );  gpwg1  = 0.0_wp
    ALLOCATE ( gpwg2      (ngp, nstepsgp) , STAT=istat );  gpwg2  = 0.0_wp
    IF (nlgw == 3) THEN
      ALLOCATE ( gpwg3    (ngp, nstepsgp) , STAT=istat );  gpwg3  = 0.0_wp
    ENDIF
    ALLOCATE ( gptso (0:ke_soil+1,ngp, nstepsgp),                  &
                                            STAT=istat );  gptso  = 0.0_wp
    ALLOCATE ( gpwso (  ke_soil+1,ngp, nstepsgp),                  &
                                            STAT=istat );  gpwso  = 0.0_wp
    ALLOCATE ( gpwice(  ke_soil+1,ngp, nstepsgp),                  &
                                            STAT=istat );  gpwice = 0.0_wp
    ALLOCATE ( gpfrsn     (ngp, nstepsgp) , STAT=istat );  gpfrsn = 0.0_wp
    ALLOCATE ( gprhsn     (ngp, nstepsgp) , STAT=istat );  gprhsn = 0.0_wp
    ALLOCATE ( gph_sn     (ngp, nstepsgp) , STAT=istat );  gph_sn = 0.0_wp
    ALLOCATE ( gptcm      (ngp, nstepsgp) , STAT=istat );  gptcm  = 0.0_wp
    ALLOCATE ( gptch      (ngp, nstepsgp) , STAT=istat );  gptch  = 0.0_wp
    ALLOCATE ( gpgz0      (ngp, nstepsgp) , STAT=istat );  gpgz0  = 0.0_wp
    ALLOCATE ( gpdpdt     (ngp, nstepsgp) , STAT=istat );  gpdpdt = 0.0_wp
    ALLOCATE ( gpsobs     (ngp, nstepsgp) , STAT=istat );  gpsobs = 0.0_wp
    ALLOCATE ( gpswdir_s  (ngp, nstepsgp) , STAT=istat );  gpswdir_s = 0.0_wp
    ALLOCATE ( gpswdifd_s (ngp, nstepsgp) , STAT=istat );  gpswdifd_s = 0.0_wp
    ALLOCATE ( gpswdifu_s (ngp, nstepsgp) , STAT=istat );  gpswdifu_s = 0.0_wp

    ALLOCATE ( gpthbs     (ngp, nstepsgp) , STAT=istat );  gpthbs = 0.0_wp
    ALLOCATE ( gppabs     (ngp, nstepsgp) , STAT=istat );  gppabs = 0.0_wp
    ALLOCATE ( gpsobt     (ngp, nstepsgp) , STAT=istat );  gpsobt = 0.0_wp
    ALLOCATE ( gpthbt     (ngp, nstepsgp) , STAT=istat );  gpthbt = 0.0_wp
    ALLOCATE ( gpalb      (ngp, nstepsgp) , STAT=istat );  gpalb  = 0.0_wp
    ALLOCATE ( gpprrs     (ngp, nstepsgp) , STAT=istat );  gpprrs = 0.0_wp
    ALLOCATE ( gpprss     (ngp, nstepsgp) , STAT=istat );  gpprss = 0.0_wp
    ALLOCATE ( gpprrc     (ngp, nstepsgp) , STAT=istat );  gpprrc = 0.0_wp
    ALLOCATE ( gpprsc     (ngp, nstepsgp) , STAT=istat );  gpprsc = 0.0_wp
    ALLOCATE ( gptmn2     (ngp, nstepsgp) , STAT=istat );  gptmn2 = 0.0_wp
    ALLOCATE ( gptmx2     (ngp, nstepsgp) , STAT=istat );  gptmx2 = 0.0_wp
    ALLOCATE ( gpvb10     (ngp, nstepsgp) , STAT=istat );  gpvb10 = 0.0_wp
    ALLOCATE ( gpabsf     (ngp, nstepsgp) , STAT=istat );  gpabsf = 0.0_wp
    ALLOCATE ( gpabgr     (ngp, nstepsgp) , STAT=istat );  gpabgr = 0.0_wp
    ALLOCATE ( gpshfl     (ngp, nstepsgp) , STAT=istat );  gpshfl = 0.0_wp
    ALLOCATE ( gplhfl     (ngp, nstepsgp) , STAT=istat );  gplhfl = 0.0_wp
    ALLOCATE ( gpumfl     (ngp, nstepsgp) , STAT=istat );  gpumfl = 0.0_wp
    ALLOCATE ( gpvmfl     (ngp, nstepsgp) , STAT=istat );  gpvmfl = 0.0_wp

    ALLOCATE ( gpplco     (ngp, nstepsgp) , STAT=istat );  gpplco = 0.0_wp
    ALLOCATE ( gplai      (ngp, nstepsgp) , STAT=istat );  gplai  = 0.0_wp
    ALLOCATE ( gproot     (ngp, nstepsgp) , STAT=istat );  gproot = 0.0_wp
    ALLOCATE ( gpt_cl     (ngp, nstepsgp) , STAT=istat );  gpt_cl = 0.0_wp
    ALLOCATE ( gpw_cl     (ngp, nstepsgp) , STAT=istat );  gpw_cl = 0.0_wp
    ALLOCATE ( gpvio3     (ngp, nstepsgp) , STAT=istat );  gpvio3 = 0.0_wp
    ALLOCATE ( gphmo3     (ngp, nstepsgp) , STAT=istat );  gphmo3 = 0.0_wp

    ALLOCATE ( gpfc       (ngp)           , STAT=istat );  gpfc   = 0.0_wp
    ALLOCATE ( gphcan     (ngp)           , STAT=istat );  gphcan= 0.0_wp
    ALLOCATE ( gpdpat     (ngp)           , STAT=istat );  gpdpat= 0.0_wp
    ALLOCATE ( gprdrg     (ngp)           , STAT=istat );  gprdrg= 0.0_wp
    ALLOCATE ( gpsai      (ngp)           , STAT=istat );  gpsai = 0.0_wp
    ALLOCATE ( gpcbig     (ngp)           , STAT=istat );  gpcbig= 0.0_wp
    ALLOCATE ( gpcsml     (ngp)           , STAT=istat );  gpcsml= 0.0_wp

#ifdef COSMOART
    ALLOCATE ( gpashmc (ke, ngp, nstepsgp) , STAT=istat ); gpashmc= 0.0_wp
    ALLOCATE ( gpcash  (ke, ngp, nstepsgp, isp_ash) , STAT=istat ); gpcash= 0.0_wp
    ALLOCATE ( gpcaero (ke, ngp, nstepsgp, isp_aerotrans+1), STAT=istat); gpcaero= 0.0_wp
#endif

  ENDIF

  IF (istat /= 0) THEN
    istat = 1
  ENDIF

!------------------------------------------------------------------------------
!  End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE alloc_gridpoints

!==============================================================================
!==============================================================================
!+ Module Routine to deallocates the allocated fields for gridpoint values
!------------------------------------------------------------------------------

SUBROUTINE dealloc_gridpoints  (istat)

!------------------------------------------------------------------------------
!
! Description:
!   This routine deallocates the allocated fields.
!
! Method:
!   DEALLOCATE statement
!
!------------------------------------------------------------------------------

! Parameters
INTEGER (KIND=iintegers), INTENT(OUT)   ::       &
  istat              ! for local error-code

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE dealloc_gridpoints  
!------------------------------------------------------------------------------

  istat = 0

  ! Arrays for the short and the long form of grid point output
  ! -----------------------------------------------------------

  DEALLOCATE ( gpps     , STAT=istat )
  DEALLOCATE ( gppp     , STAT=istat )
  DEALLOCATE ( gpt2m    , STAT=istat )
  DEALLOCATE ( gptd2m   , STAT=istat )
  DEALLOCATE ( gpu10m   , STAT=istat )
  DEALLOCATE ( gpv10m   , STAT=istat )
  DEALLOCATE ( gpt_g    , STAT=istat )
  DEALLOCATE ( gprrsn   , STAT=istat )
  DEALLOCATE ( gprssn   , STAT=istat )
  DEALLOCATE ( gprrkn   , STAT=istat )
  DEALLOCATE ( gprskn   , STAT=istat )
  DEALLOCATE ( gpwsno   , STAT=istat )
  DEALLOCATE ( gpclcl   , STAT=istat )
  DEALLOCATE ( gpclcm   , STAT=istat )
  DEALLOCATE ( gpclch   , STAT=istat )
  DEALLOCATE ( gpclct   , STAT=istat )
  DEALLOCATE ( gpsmlt   , STAT=istat )

  DEALLOCATE ( gpp0     , STAT=istat )
  DEALLOCATE ( gphsur   , STAT=istat )
  DEALLOCATE ( gpfrla   , STAT=istat )
  DEALLOCATE ( gpsoil   , STAT=istat )
  DEALLOCATE ( gplati   , STAT=istat )
  DEALLOCATE ( gplong   , STAT=istat )

  IF (lgpshort) THEN

    ! Arrays only for the short form of grid point output
    ! ---------------------------------------------------

    DEALLOCATE ( gpu950   , STAT=istat )
    DEALLOCATE ( gpv950   , STAT=istat )
    DEALLOCATE ( gpu850   , STAT=istat )
    DEALLOCATE ( gpv850   , STAT=istat )
    DEALLOCATE ( gpu700   , STAT=istat )
    DEALLOCATE ( gpv700   , STAT=istat )
    DEALLOCATE ( gpu500   , STAT=istat )
    DEALLOCATE ( gpv500   , STAT=istat )
    DEALLOCATE ( gptp     , STAT=istat )
    DEALLOCATE ( gpt850   , STAT=istat )
    DEALLOCATE ( gpt700   , STAT=istat )
    DEALLOCATE ( gpt500   , STAT=istat )
    DEALLOCATE ( gpfog    , STAT=istat )
    DEALLOCATE ( gphbas   , STAT=istat )
    DEALLOCATE ( gphtop   , STAT=istat )

 ELSE IF (lgplong .OR. lgpspec) THEN

    ! Arrays only for the long form of grid point output
    ! --------------------------------------------------

    DEALLOCATE ( gpu   , STAT=istat )
    DEALLOCATE ( gpv   , STAT=istat )
    DEALLOCATE ( gpw   , STAT=istat )
    DEALLOCATE ( gpt   , STAT=istat )
    DEALLOCATE ( gpqv  , STAT=istat )
    DEALLOCATE ( gpqc  , STAT=istat )
    DEALLOCATE ( gpqi  , STAT=istat )
    DEALLOCATE ( gphhl , STAT=istat )
    DEALLOCATE ( gphml , STAT=istat )
    DEALLOCATE ( gptkvm, STAT=istat )
    DEALLOCATE ( gptkvh, STAT=istat )
    DEALLOCATE ( gpclcs, STAT=istat )
    DEALLOCATE ( gptsno, STAT=istat )
    DEALLOCATE ( gpt_s , STAT=istat )
    DEALLOCATE ( gpt_m , STAT=istat )
    DEALLOCATE ( gpqvs , STAT=istat )
    DEALLOCATE ( gpw_i , STAT=istat )
    DEALLOCATE ( gpwg1 , STAT=istat )
    DEALLOCATE ( gpwg2 , STAT=istat )
    IF (nlgw == 3) THEN
      DEALLOCATE ( gpwg3 , STAT=istat )
    ENDIF
    DEALLOCATE ( gptso , STAT=istat )
    DEALLOCATE ( gpwso , STAT=istat )
    DEALLOCATE ( gpwice, STAT=istat )
    DEALLOCATE ( gpfrsn, STAT=istat )
    DEALLOCATE ( gprhsn, STAT=istat )
    DEALLOCATE ( gph_sn, STAT=istat )
    DEALLOCATE ( gptcm , STAT=istat )
    DEALLOCATE ( gptch , STAT=istat )
    DEALLOCATE ( gpgz0 , STAT=istat )
    DEALLOCATE ( gpdpdt, STAT=istat )
    DEALLOCATE ( gpsobs, STAT=istat )
    DEALLOCATE ( gpswdir_s, STAT=istat )
    DEALLOCATE ( gpswdifd_s, STAT=istat )
    DEALLOCATE ( gpswdifu_s, STAT=istat )

    DEALLOCATE ( gpthbs, STAT=istat )
    DEALLOCATE ( gppabs, STAT=istat )
    DEALLOCATE ( gpsobt, STAT=istat )
    DEALLOCATE ( gpthbt, STAT=istat )
    DEALLOCATE ( gpalb , STAT=istat )
    DEALLOCATE ( gpprrs, STAT=istat )
    DEALLOCATE ( gpprss, STAT=istat )
    DEALLOCATE ( gpprrc, STAT=istat )
    DEALLOCATE ( gpprsc, STAT=istat )
    DEALLOCATE ( gptmn2, STAT=istat )
    DEALLOCATE ( gptmx2, STAT=istat )
    DEALLOCATE ( gpvb10, STAT=istat )
    DEALLOCATE ( gpabsf, STAT=istat )
    DEALLOCATE ( gpabgr, STAT=istat )
    DEALLOCATE ( gpshfl, STAT=istat )
    DEALLOCATE ( gplhfl, STAT=istat )
    DEALLOCATE ( gpumfl, STAT=istat )
    DEALLOCATE ( gpvmfl, STAT=istat )

    DEALLOCATE ( gpplco, STAT=istat )
    DEALLOCATE ( gplai , STAT=istat )
    DEALLOCATE ( gproot, STAT=istat )
    DEALLOCATE ( gpfc  , STAT=istat )
    DEALLOCATE ( gpt_cl, STAT=istat )
    DEALLOCATE ( gpw_cl, STAT=istat )
    DEALLOCATE ( gpvio3, STAT=istat )
    DEALLOCATE ( gphmo3, STAT=istat )

    DEALLOCATE ( gphcan, STAT=istat )
    DEALLOCATE ( gpdpat, STAT=istat )
    DEALLOCATE ( gprdrg, STAT=istat )
    DEALLOCATE ( gpsai , STAT=istat )
    DEALLOCATE ( gpcbig, STAT=istat )
    DEALLOCATE ( gpcsml, STAT=istat )

#ifdef COSMOART
    IF (l_cosmo_art) THEN
      IF (lvolcano) THEN
        DEALLOCATE ( gpcash , STAT=istat )
        DEALLOCATE ( gpashmc, STAT=istat )
      ELSE IF (ldust .AND. lonly_dust) THEN
        DEALLOCATE ( gpcaero, STAT=istat )
      END IF
    END IF
#endif

  ENDIF

!------------------------------------------------------------------------------
!  End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE dealloc_gridpoints 

!==============================================================================
!==============================================================================
!+ Module function to make linear interpolations
!------------------------------------------------------------------------------

REAL (KIND=wp) FUNCTION flin(a,t,x,xs,ne)

!------------------------------------------------------------------------------
!
! Description:
!   flin returns the value of linear interpolation belonging to 
!   the argument a by the help of the fixed points (t,x).
!   The abszissa t must be monotonly increasing!
!
! Method:
!   linear interpolation
!
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!
! Modules used: none
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------

    INTEGER (KIND=iintegers)   ::  i,ne
    REAL    (KIND=wp)          ::  a,q,xs,t(1:ne),x(1:ne)

    x(ne)=xs
    i=1
1   IF (i.le.ne) THEN
       IF (t(i).lt.a) THEN
          i=i+1
          GOTO 1
       ELSEIF (i.eq.1) THEN
          flin=x(1) ! untere Extrapolation, falls a<t(1)
       ELSE
          q=(x(i)-x(i-1))/(t(i)-t(i-1))
          flin=x(i-1)+q*(a-t(i-1))
       ENDIF
    ELSE
       flin=x(ne) ! upper point of extrapolation, if a>t(ne)
    ENDIF

END FUNCTION flin

!===============================================================================

END MODULE src_gridpoints 
