!+ Source module for computing the leapfrog dynamical time stepping
!------------------------------------------------------------------------------

MODULE src_leapfrog

!------------------------------------------------------------------------------
!
! Description:
!   The module "src_leapfrog" performs one time step of the integration
!   of the spatially disretised thermodynamical equations. 
!   Driving routine is the model procedure "org_leapfrog", which
!   calls the rountines required and does some diagnostics.
!   The additional subroutine gauss solves a system of tridiagonal matrix
!   equations; it is called from slow_tendencies.
!   The additional subroutine satad does the saturation adjustment for
!   t, qv and qc; it is called from sardass.
!
! Current Code Owner: DWD, Michael Baldauf
!  phone:  +49  69  8062 2733
!  fax:    +49  69  8062 3721
!  email:  Michael.Baldauf@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.34       1999/12/10 Guenther Doms
!  Initial release
! 1.37       2000/03/24 Guenther Doms
!  Horizontal diffusion is now applied as a numerical filter at the
!  end of a timestep for timelevel nnew. This allows the timestep to be closer
!  to the CFL limit.
! 1.39       2000/05/03 Ulrich Schaettler
!  Changed some variable names and adapted call to timing routine get_timings.
!  Variables depended on dt are now calculated every time step (needed for
!  the nesting version).
! 2.2        2000/08/18 Guenther Doms
!  An error in the discretization of the v-velocity advection (kinetic energy
!  term) was corrected (only negligible impact).
! 2.6        2001/06/12 Guenther Doms
!  Correction in a data exchange related to the cloud-ice scheme.
! 2.8        2001/07/06 Ulrich Schaettler
!  Now use actual surface fluxes instead of summarized;
!  Corrected bug in time-measuring
! 2.9        2001/07/16 Guenther Doms
!  Introduction of the optional use of monotonic horizontal diffusion with
!  orographic limiting for the thermodynamic variables T, qv and qc.
!  The hole-filling scheme for negative qv/qc values has been changed to
!  a pure vertical redistribution scheme.
! 2.10       2001/07/24 Guenther Doms
!  The module routine slow_tendencies has been has been removed from
!  this module to become an external subroutine (slow_tendencies.f90).
! 2.17       2002/05/08 Ulrich Schaettler
!  Adaptations to the Kain-Fritsch convection scheme
! 2.18       2002/07/16 Ulrich Schaettler
!  Modifications due to changes in the 2 time level scheme.
! 3.2        2003/02/07 Ulrich Schaettler
!  Moved part from horizontal diffusion to external subroutine hori_diffusion.
! 3.5        2003/09/02 Ulrich Schaettler
!  Optimizations in the global communications
! 3.7        2004/02/18 Michael Baldauf
!  Added semi-Lagrange advection of precipitation
!  Determination of time level for qi depending on lprog_qi and lprogprec
!  Replaced local variable xkd by global Namelist parameter xkd
!     (by Jochen Foerstner)
!  Renamed cphi (crlat), acphir (acrlat)
! 3.12       2004/09/15 Ulrich Schaettler
!  Bug fix for variable zqiwt (has to be allocated with ke+1).
! 3.13       2004/12/03 Ulrich Schaettler
!  Modifications to run with latent heat nudging (Klaus Stephan, et.al.)
!  Introduction of graupel scheme; Modifications for 2D Version
!                                                (Thorsten Reinhardt)
! 3.16       2005/07/22 Ulrich Schaettler
!  Adapted call to procedure hori_diffusion to module procedure comp_hori_diff
! 3.18       2006/03/03 Ulrich Schaettler / Klaus Stephan
!  Changed treatment of ASCII files for introducing restart possibility
!  LHN namelist parameter moved to data_lheat_nudge to avoid to many dependencies
!  LHN bug fix: itype_gscp >= 3 to consider the graupel scheme as well
!  Introduction of a dynamical bottom boundary condition (after Gassmann)
! 3.21       2006/12/04 MeteoSwiss / Ulrich Schaettler
!  Modifications to use Bechtold convection scheme
! V3_23        2007/03/30 Ulrich Schaettler
!  Introduced idbg_level and debug output
! V4_4         2008/07/16 Ulrich Schaettler, Dmitrii Mironov
!  Convective tendencies of qc and qi (computed by Tiedtke scheme) are added 
!  to the total tendencies qctens and qitens (DM)
!  Adapted interface of get_timings
!  Moved SQRT out of loops while computing CFL (by. J.-O. Beismann)
!  Replaced lkainfri, lbechtol by itype_conv
! V4_8         2009/02/16 Guenther Zaengl
!  Use p0hl (reference pressure at half levels) for full consistency with
!  new reference atmosphere implementation
! V4_9         2009/07/16 Ulrich Schaettler, Heike Vogel
!  Implement dynamics for COSMO-ART and POLLEN in a beta-version
! V4_12        2010/05/11 Ulrich Schaettler, Oli Fuhrer
!  Removed t0(_melt)
!  Adapted call to SR comp_hori_diffusion
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_17        2011/02/24 Ulrich Blahak
!  Adapted interface of exchg_boundaries; corrected kzdims(1:20) -> kzdims(1:24);
!  eliminated my_peri_neigh; added possibility to output max. V_h and min./max. W
!  to stdout for diagnostic purposes, if (ldebug_dyn .AND. idbg_level > 3)
! V4_18        2011/05/26 Ulrich Schaettler
!   Introduced conditional compilation for Nudging
!   Add (Tiedtke) convective tendencies for Pollen and COSMO-ART (Christoph Knote)
! V4_23        2012/05/10 Ulrich Schaettler, Oliver Fuhrer
!  Removed switches lprogprec, ltrans_prec
!  Removed computations of total physical tendencies and moved them to 
!    organize_physics (Oli Fuhrer)
!  Bug Fix for computing MAX(w): this should only be done in the interior of the 
!    domain to give reproducible results
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer, Ulrich Blahak, Hans-Juergen Panitz
!  Replaced qx-variables by using them from the tracer module
!  UB: Implemented internal switch "l2mom_do_extra_satads" to be able
!   to switch on the extra saturation adjustments outside the microphysics parts.
!  Introduced nexch_tag for MPI boundary exchange tag to replace ntstep (HJP)
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  MESSy interface introduced
! V4_28        2013/07/12 KIT, Ulrich Schaettler
!  Changes to adapt COSMO-ART to new tracer module: all dependencies to 
!  COSMOART and POLLEN deleted, because this is now handled by the tracer module
!  Use subroutines and variables for vertical grid and reference atmospheres 
!    from module vgrid_refatm_utils
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler, Davide Cesari
!  Unification of MESSy interfaces and COSMO Tracer structure
!  For the COSMO-Model only use vcoord from vgrid_refatm_utils
!  Added missing argument kflat to fast-waves-call (DC)
! V4_30        2013/11/08 Ulrich Schaettler (Igor Andruska)
!  Added missing argument kflat also to semi-implicit-call
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
! V5_2         2015-05-21 Ulrich Schaettler
!  Reorganized numbering of MPI datatypes ID
! V5_3         2015-10-09 Ulrich Blahak, Michael Baldauf
!  Added storage of LH by satad on ttens_diab for output of pure diabatic tendencies (UB)
!  New Namelist switch l_euler_dynamics (MB)
! V5_4         2016-03-10 Oliver Fuhrer
!  Editorial changes
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

USE data_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------
    je_tot,       & ! number of grid points in meridional direction
    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ieje,         & ! ie * je
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
    istart,       & ! start index for the forecast of w, t, qv, qc and pp
    iend,         & ! end index for the forecast of w, t, qv, qc and pp
    istartu,      & ! start index for the forecast of u
    iendu,        & ! end index for the forecast of u
    istartv,      & ! start index for the forecast of v
    iendv,        & ! end index for the forecast of v
    istartpar,    & ! start index for computations in the parallel program
    iendpar,      & ! end index for computations in the parallel program

!   meridional direction
    jstart,       & ! start index for the forecast of w, t, qv, qc and pp
    jend,         & ! end index for the forecast of w, t, qv, qc and pp
    jstartu,      & ! start index for the forecast of u
    jendu,        & ! start index for the forecast of u
    jstartv,      & ! start index for the forecast of v
    jendv,        & ! end index for the forecast of v
    jstartpar,    & ! start index for computations in the parallel program
    jendpar,      & ! end index for computations in the parallel program

! 4. constants for the horizontal rotated grid and related variables
! ------------------------------------------------------------------

    degrad,       & ! factor for transforming degree to rad
    eddlon,       & ! 1 / dlon
    eddlat,       & ! 1 / dlat
    edadlat,      & ! 1 / (radius of the earth * dlat)

! 4. variables for the time discretization and related variables
! --------------------------------------------------------------

    dt,           & ! long time-step
    ed2dt,        & ! 1 / (2 * dt)
    dt2,          & ! 2 * dt
    epsass,       & ! eps for the Asselin-filter
    betasw,       & ! beta-variable for treatment of soundwaves
    vhmx_vol,     & ! maximum absolute horizontal wind in total model domain
    vhmx_cfl,     & ! maximum absolute horizontal wind velocity from CFL

! 8. Organizational variables to handle the COSMO humidity tracers
! ----------------------------------------------------------------
    idt_qv, idt_qc

! end of data_modelconfig
!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 1. mathematical constants
! -------------------------

    pi,           & ! circle constant

! 2. physical constants and related variables
! -------------------------------------------

    r_d,          & ! gas constant for dry air
    r_v,          & ! gas constant for water vapor
    rdv,          & ! r_d / r_v
    o_m_rdv,      & ! 1 - r_d/r_v
    rvd_m_o,      & ! r_v/r_d - 1
    cp_d,         & ! specific heat of dry air at constant pressure
    cpdr,         & ! 1 / cp_d
    rdocp,        & ! r_d / cp_d
    rcpv,         & ! wcp_d/wcp_v - 1
    rcpl,         & ! wcp_d/wcp_l - 1
    gamma,        & ! 1 / (1 - rdocp)   ( = cp_d/cv_d)
    lh_v,         & ! latent heat of vapourization
    lh_f,         & ! latent heat of fusion
    lh_s,         & ! latent heat of sublimation
    g,            & ! acceleration due to gravity
    gq,           & ! g*g                        
    gr,           & ! 1 / g                      
    r_earth,      & ! mean radius of the earth

! 3. constants for parametrizations
! ---------------------------------

    b1,           & ! variables for computing the saturation vapour pressure
    b2w,          & ! over water (w) and ice (i)
    b2i,          & !               -- " --
    b3,           & !               -- " --
    b4w,          & !               -- " --
    b4i,          & !               -- " --
    b234w,        & ! b2w * (b3 - b4w)
    aks2,         & ! variable for horizontal diffusion of second order
    aks4            ! variable for horizontal diffusion of fourth order

! end of data_constants

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------

    rho0       ,    & ! reference density at the full model levels    (kg/m3)
    dp0        ,    & ! reference pressure thickness of layers        ( Pa)
    p0         ,    & ! reference pressure at main levels             ( Pa)
    p0hl       ,    & ! reference pressure at half levels             ( Pa)
    hhl        ,    & ! geometical height of half model levels        ( m )

! 2. external parameter fields                                        (unit)
! ----------------------------

    fc         ,    & ! coriolis-parameter                            ( 1/s )
    crlat      ,    & ! cosine of transformed latitude
    acrlat     ,    & ! 1 / ( crlat * radius of the earth )           ( 1/m )

! 3. prognostic variables                                             (unit)
! -----------------------

    u          ,    & ! zonal wind speed                              ( m/s )
    v          ,    & ! meridional wind speed                         ( m/s )
    w          ,    & ! vertical wind speed (defined on half levels)  ( m/s )
    t          ,    & ! temperature                                   (  k  )
    pp         ,    & ! deviation from the reference pressure         ( pa  )

! 4. tendency fields for the prognostic variables                     (unit )
! -----------------------------------------------
!    timely deviation  by diabatic and adiabatic processes 
!    without sound-wave terms

    utens        ,  & ! u-tendency without sound-wave terms           ( m/s2)
    vtens        ,  & ! v-tendency without sound-wave terms           ( m/s2)
    wtens        ,  & ! w-tendency without sound-wave terms           ( m/s2       
                      ! (defined on half levels )
    ttens        ,  & ! t-tendency without sound-wave terms           ( K/s )
    pptens       ,  & ! pp-tendency without sound-wave terms          (Pa/s )
    ttens_diab   ,  & ! temperature tendency due to pure diabatic processes for output ( K/s )

! 5. fields for surface values and soil model variables               (unit )
! -----------------------------------------------------

    ps                ! surface pressure                              ( pa  )

USE data_fields     , ONLY :   &

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------

    rho       ,     & ! density of moist air

!   fields of the precipitation
    qrs        ,    & ! precipitation water (water loading)           (kg/kg)

!   fields that are computed in the dynamics
    dqvdt             ! threedimensional moisture convergence         (  1/s)

! end of data_fields

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    nstart,       & ! first time step of the forecast
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
    nold,         & ! corresponds to ntstep - 1
    nnow,         & ! corresponds to ntstep
    nnew,         & ! corresponds to ntstep + 1

! 3. controlling the physics
! --------------------------
    itype_gscp,   & ! type of grid-scale precipitation physics
    l_cosmo_art,  & ! if .TRUE., run the COSMO_ART
    l_pollen,     & ! of pollen

! 4. controlling the dynamics
! ---------------------------
    l_euler_dynamics, & ! on/off Euler-solver; however, Coriolis force and 
                        ! tracer advection may be performed
    lsemi_imp,    & ! if .TRUE.,  running with semi-implicit scheme,
                    ! if .FALSE., running with split-explicit scheme
    ldyn_bbc,     & ! dynamical bottom boundary condition

! 6. controlling the upper boundary condition
! -------------------------------------------
    lrubc ,       & ! with radiative upper boundary condition

! 7. additional control variables
! -------------------------------
    lcond,        & ! forecast with condensation/evaporation
    loutput_diab, & ! internal switch for output of temperature tendencies due 
                    ! to pure diabatic processes
    l2mom_satads, & ! in case of 2-moment scheme, do all the satads
                    ! (like for the 1-moment schemes), not just the
                    ! satad after the microphysics at the end of the timestep.
    lperi_x,        & ! if lartif_data=.TRUE.: periodic boundary conditions (.TRUE.) in x-dir.
                    !                 or with Davies conditions (.FALSE.)
    lperi_y,        & ! if lartif_data=.TRUE.: periodic boundary conditions (.TRUE.) in y-dir.
                    !                 or with Davies conditions (.FALSE.)
    l2dim,        & ! if lartif_data=.TRUE.: 2dimensional model version (.TRUE.) or
                    !                 full 3dimensional version (.FALSE.)
    lmetr,        & ! if lartif_data=.TRUE.: with metric terms (.TRUE.)
                    !                 or without metric terms (.FALSE.)
    ltime,        & ! detailled timings of the program are given
    lreproduce,   & ! the results are reproducible in parallel mode
    lhordiff,     & ! running with horizontal diffusion
    xkd,          & ! coefficient for divergence damping
    itype_hdiff,  & ! type of horizontal diffusion (=1: 4th order linear),
                    ! =2: 4th order linear monotonic with orographic limit)
    idbg_level,   & ! to control the verbosity of debug output
    ldebug_dyn,   & ! if .TRUE., debug output for dynamics
    lprintdeb_all,& ! .TRUE.:  all tasks print debug output
                    ! .FALSE.: only task 0 prints debug output

! 12. controlling verbosity of debug output and miscellaneous
! -----------------------------------------------------------
    linit_fields    ! to initialize also local variables with a default value
                    ! (some compilers seem to have problems, if fields are not
                    !  initialized, even if values are never used;
                    !  also, debugging could be easier then)

! end of data_runcontrol 

!------------------------------------------------------------------------------

USE data_parallel,      ONLY :  &
    ldatatypes,      & ! if .TRUE.: use MPI-Datatypes for some communications
    ltime_barrier,   & ! if .TRUE.: use additional barriers for determining the
                       ! load-imbalance
    num_compute,     & ! number of compute PEs
    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains
    ncomm_type,      & ! type of communication for boundary exchange
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    my_cart_pos,     & ! position of this subdomain in the cartesian grid
                       ! in x- and y-direction
    my_cart_neigh,   & ! neighbors of this subdomain in the cartesian grid
!    my_peri_neigh,   & ! periodic neighbors of this subdomain in the periodic 
                       ! grid; if it lies at the boundary of the total domain
    isubpos,         & ! positions of the subdomains in the total domain. Given
                       ! are the i- and the j-indices of the lower left and the
                       ! upper right grid point in the order
                       !                  i_ll, j_ll, i_ur, j_ur.
                       ! Only the interior of the domains are considered, not
                       ! the boundary lines.
    icomm_cart,      & ! communicator for the virtual cartesian topology
    iexch_req,       & ! stores the sends requests for the neighbor-exchange
                       ! that can be used by MPI_WAIT to identify the send
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    imp_integers,    & ! determines the correct INTEGER type used in the model
                       ! for MPI
    nexch_tag,       & ! tag to be used for MPI boundary exchange
                       !  (in calls to exchg_boundaries)
    sendbuf,         & ! sending buffer for boundary exchange:
                       ! 1-4 are used for sending, 5-8 are used for receiving
    isendbuflen        ! length of one column of sendbuf

!------------------------------------------------------------------------------

USE data_tracer,         ONLY :  T_ADV_ID, T_ADV_ON,                          &
                                 T_ERR_NOTFOUND

!------------------------------------------------------------------------------

USE environment,         ONLY :  comm_barrier, model_abort, exchg_boundaries
USE time_utilities,      ONLY :  get_timings, i_dyn_computations,           &
      i_communications_dyn, i_barrier_waiting_dyn, i_horizontal_advection, &
      i_complete_tendencies, i_slow_tendencies
USE meteo_utilities,     ONLY :  satad
USE parallel_utilities,  ONLY :  global_values
USE numeric_utilities,   ONLY :  hadv_pd, vadv_pd, hadv_cd2,               &
      backtraj_trilin_dt1_3tl, backtraj_trilin_dt2_3tl, interpol_sl_trilin
USE hori_diffusion,      ONLY :  comp_hori_diff
USE vgrid_refatm_utils,  ONLY :  vcoord

USE src_tracer,          ONLY :  trcr_get, trcr_meta_get, trcr_errorstr,      &
                                 trcr_get_ntrcr, trcr_get_index

#ifdef NUDGING
USE data_lheat_nudge,    ONLY :     &
    llhn,         & ! main switch for latent heat nudging (lhn)
    tt_lheat        ! profile of t-increments due to latent heating   ( K/s )
                    ! (stored for current and previous timestep)

USE src_lheating,        ONLY :  get_gs_lheating
#endif

!------------------------------------------------------------------------------
! Declarations
!------------------------------------------------------------------------------

IMPLICIT NONE

REAL (KIND = wp),     ALLOCATABLE :: & 
                    ! only used in Leapfrog integration scheme
  zwcon (:,:,:)     ! mass weighted contravariant vertical velocity
                    ! only used per time step

REAL (KIND = wp),     ALLOCATABLE :: &
                    ! only used in the semi-implicit integration scheme
  zqu (:,:,:),    & ! u right-hand side
  zqv (:,:,:),    & ! v right-hand side
  zqw (:,:,:),    & ! w right-hand side
  zqp (:,:,:),    & ! p right-hand side
  zqt (:,:,:),    & ! T right-hand side
  zqws(:,:,:),    & ! auxilary w right-hand side
  zqps(:,:,:)       ! wave equation right-hand side

REAL (KIND = wp),     ALLOCATABLE :: &
                    ! only used in the semi-implicit integration scheme
  za (:,:,:),     & ! lower diagonal of N matrix
  zb (:,:,:),     & ! main  diagonal of N matrix
  zc (:,:,:),     & ! upper diagonal of N matrix
  za1(:,:,:),     & ! lower diagonal of Dz1
  zb1(:,:,:),     & ! main  diagonal of Dz1
  za2(:,:,:),     & ! lower diagonal of Dz2
  zb2(:,:,:),     & ! main  diagonal of Dz2
  za3(:,:,:),     & ! lower diagonal of Dz2.Ni.Dz1 matrix
  zb3(:,:,:),     & ! main  diagonal of Dz2.Ni.Dz1 matrix
  zc3(:,:,:)        ! upper diagonal of Dz2.Ni.Dz1 matrix

REAL (KIND = wp)                  :: &
    dts,          & ! small time step for time splitting     
    epsray          ! Rayleigh damping coefficient for vhmx_vol > vhmx_cfl

INTEGER (KIND=iintegers)          :: &
    ismtstep,     & ! number of small time steps in leapfrog intervall 2*dt
    ismtstep_h,   & ! number of small time steps in time-intervall dt/2
    isp             ! COSMO-ART species index

!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure in "src_leapfrog" for organizing the time stepping  
!------------------------------------------------------------------------------

SUBROUTINE org_leapfrog

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "dynamics" is the driving routine of the
!   module, i.e. it acts as interface to the main program.
!
! Method:
!
!------------------------------------------------------------------------------

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    i,j,k,               & !  Loop indices
    kitpro,              & !  number of iterations in the saturation adjustment
    izstata,             & !  error status at allocation
    izstatd,             & !  error status at deallocation
    izerror, izdebug       !

  REAL    (KIND=wp)        ::  &
    zcs, zdtsmax, zvb, zdtsr

  CHARACTER (LEN=255)      ::  &
    yzerrmsg

  CHARACTER (LEN=25)       :: yzroutine

! Local (automatic) arrays:
! -------------------------
  REAL    (KIND=wp)        ::  &
    zuvwmax(0:ke)      ! maximum vertical velocity in a horizontal layer

! for fast_waves
  REAL    (KIND=wp)        ::  &
  zubdt_west (je,ke),& ! u-boundary tendency at west boundary (total domain)
  zubdt_east (je,ke),& ! u-boundary tendency at east boundary (total domain)
  zvbdt_south(ie,ke),& ! v-boundary tendency at south boundary (total domain)
  zvbdt_north(ie,ke)   ! v-boundary tendency at north boundary (total domain)

! from complete_tendencies
  REAL    (KIND=wp)        ::  &
    zgrhoc, zbqt  , zbqb , zbqa, vhmx_save

! Tracer pointers
REAL (KIND=wp),     POINTER :: &
  qv    (:,:,:) => NULL(),     &          ! QV at nnew
  qv_now(:,:,:) => NULL(),     &          ! QV at nnow
  qc_now(:,:,:) => NULL(),     &          ! QC at nnow
  qc    (:,:,:) => NULL()                 ! QC at nnew

! For debug output of MAX(U) and MIN/MAX(W) every timestep to stdout:
REAL    (KIND=wp)         ::   &
     vhmx_loc, wmx_loc, wmn_loc

! End of header
!==============================================================================
 
!------------------------------------------------------------------------------
! Begin Subroutine org_leapfrog
!------------------------------------------------------------------------------

izstata = 0_iintegers
izstatd = 0_iintegers
izerror = 0_iintegers

yzroutine = 'org_leapfrog'

! Initialize, whether debug output shall be done
IF (ldebug_dyn) THEN
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

! Retrieve the required microphysics tracers
CALL trcr_get(izerror, idt_qv, ptr_tlev = nnow, ptr = qv_now)
IF (izerror /= 0) THEN
  yzerrmsg = trcr_errorstr(izerror)
  CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
ENDIF
CALL trcr_get(izerror, idt_qv, ptr_tlev = nnew, ptr = qv)
IF (izerror /= 0) THEN
  yzerrmsg = trcr_errorstr(izerror)
  CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
ENDIF
CALL trcr_get(izerror, idt_qc, ptr_tlev = nnow, ptr = qc_now)
IF (izerror /= 0) THEN
  yzerrmsg = trcr_errorstr(izerror)
  CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
ENDIF
CALL trcr_get(izerror, idt_qc, ptr_tlev = nnew, ptr = qc)
IF (izerror /= 0) THEN
  yzerrmsg = trcr_errorstr(izerror)
  CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
ENDIF

!------------------------------------------------------------------------------
! Section 1: Some preparations at the initial time step(s):
!------------------------------------------------------------------------------

IF (izdebug > 10) THEN
  PRINT *, '        START org_leapfrog'
ENDIF

! Calculation of the small time step size dts and the number of small steps
! (ismtstep) to cover a 2*dt leapfrog time intervall.

! Constant factor for divergence-type damping
zcs       = SQRT( gamma*r_d*303.0_wp )
zdtsmax   = 1.0_wp/zcs/(1.0_wp+xkd)/edadlat/SQRT(2.0_wp)
zdtsmax   = MIN ( zdtsmax, 32.0_wp )
ismtstep  = INT (ABS(dt2)/zdtsmax, iintegers) + 1
ismtstep  = MAX (1,ismtstep)
dts       = dt2/ismtstep
zdtsr     = 1.0_wp/dts

IF (ntstep <= nstart+1) THEN
  IF ( (my_cart_id == 0) .AND. (idbg_level > 1) ) THEN
    WRITE (*,'(A,F9.2,A,F9.2,A,F9.2)') 'CS=',zcs,'DT2=',dt2,'DTSMAX=' ,zdtsmax
    WRITE (*,'(A,F9.2,A,I5,A,F9.2,A,F9.2)') 'SHORT TIME STEP ',dts,          &
               ' ismtstep=', ismtstep,' BETA = ', betasw, ' XKD = ',xkd
    WRITE (*,'(A,I5)') 'LEVELINDEX WHERE LEVELS BECOME FLAT, KFLAT = ', vcoord%kflat
  ENDIF
ENDIF

IF ( (ntstep == nstart) .AND. lrubc .AND. my_cart_id == 0) THEN
    PRINT *,' CAUTION!!: RADIATIVE UPPER BOUNDARY CONDITION NOT IMPLEMENTED'
    PRINT *,' CAUTION!!: MODEL RUN USES RIGID UPPER LID'
ENDIF

! Determine the min/.max. W in the model domain and the max |vh| and write to stdout
! for diagnostic purposes:
IF (ldebug_dyn .AND. idbg_level > 3) THEN
  vhmx_loc = SQRT(MAXVAL( u(:,:,:,nnow)**2 + v(:,:,:,nnow)**2 ))
  wmx_loc  = MAXVAL(w(:,:,:,nnow))
  wmn_loc  = MINVAL(w(:,:,:,nnow))
  IF (num_compute > 1) THEN
    CALL global_values (vhmx_loc, 1, 'MAX', imp_reals, icomm_cart, -1,    &
         yzerrmsg, izerror)
    CALL global_values (wmx_loc, 1, 'MAX', imp_reals, icomm_cart, -1,    &
         yzerrmsg, izerror)
    CALL global_values (wmn_loc, 1, 'MIN', imp_reals, icomm_cart, -1,    &
         yzerrmsg, izerror)
  ENDIF
  IF ( my_cart_id == 0 ) THEN
    WRITE (*,'(a,es12.5)') 'Maximum V_h : ', vhmx_loc
    WRITE (*,'(a,es12.5,a,es12.5)') 'Maximum W : ', wmx_loc, '     Minimum W : ', wmn_loc
  END IF
END IF

! Determine the maximum horizontal wind velocity at timelevel nnow
! and set the Rayleigh damping coefficients if vhmx_vol approaches the
! stability limit (vhmx_cfl).
 
! Determine the maxima of w(:,:,k,:) for reproducible results in the
! parallel program when computing the saturation adjustment at the end
! of the leapfrog scheme
DO k = 1,ke
  zuvwmax(k) = MAXVAL ( ABS(w(istartpar:iendpar,jstartpar:jendpar,k,nnow)) )
ENDDO

vhmx_vol = 0.0_wp
DO k = 1,ke
  DO j = jstartu,jendu
    DO i = istartv,iendv
      zvb      = u(i,j,k,nnow)**2 + v(i,j,k,nnow)**2
      vhmx_vol = MAX ( vhmx_vol, zvb)
    ENDDO
  ENDDO
ENDDO
vhmx_vol   = SQRT (vhmx_vol)
zuvwmax(0) = vhmx_vol
vhmx_save  = vhmx_vol

IF (num_compute > 1) THEN
  IF (ltime) THEN
    CALL get_timings (i_dyn_computations, ntstep, dt, izerror)
    IF (ltime_barrier) THEN
      CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
      CALL get_timings (i_barrier_waiting_dyn, ntstep, dt, izerror)
    ENDIF
  ENDIF

  IF (lreproduce) THEN
    CALL global_values (zuvwmax(0:ke), ke+1, 'MAX', imp_reals, icomm_cart, -1,&
                        yzerrmsg, izerror)
  ELSE
    CALL global_values (zuvwmax(0), 1, 'MAX', imp_reals, icomm_cart, -1,    &
                        yzerrmsg, izerror)
  ENDIF
  vhmx_vol = zuvwmax(0)

  IF (ltime) CALL get_timings (i_communications_dyn, ntstep, dt, izerror)
ENDIF

IF ( vhmx_vol > 0.95_wp*vhmx_cfl ) THEN
  epsray = 1.0_wp/(2000.0_wp*dt)
  epsray = epsray*(vhmx_vol - 0.95_wp*vhmx_cfl)/(0.05_wp*vhmx_cfl)

  IF ( (vhmx_vol > vhmx_cfl) .AND. (vhmx_vol == vhmx_save) ) THEN
    ! this PE has a violation: look for the indices
    DO k = 1,ke
      DO j = jstartu,jendu
        DO i = istartv,iendv
          zvb      = SQRT( u(i,j,k,nnow)**2 + v(i,j,k,nnow)**2 )
          IF (zvb == vhmx_save) THEN
            PRINT*,' *** WARNING *** vhmax exceeds CFL-value: ', &
             vhmx_vol, ':  i=',i,'; j=',j,'; PE=',my_cart_id
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDIF

ELSE
  epsray = 0.0_wp
ENDIF

!------------------------------------------------------------------------------
! Section 2: Compute the time tendencies due to horizontal diffusion
!            and add them to the tendency fields
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Section 3: Allocate and store boundary tendencies for u and v in case of
!            non-periodic boundary conditions.
!------------------------------------------------------------------------------

! These tendencies are required to compute the wind divergence along the
! boundary lines in the fast-wave solver.
! The boundary tendencies are premultiplied by the small time step

IF ( .NOT.lperi_x ) THEN
  DO k = 1, ke
    IF (my_cart_neigh(1) == -1) THEN
      DO  j = jstartu, jendu
        zubdt_west(j,k) = (u(istartu-1,j,k,nnew)-u(istartu-1,j,k,nold))*dts/dt2
      ENDDO
    ENDIF
    IF (my_cart_neigh(3) == -1) THEN
      DO  j = jstartu, jendu
        zubdt_east(j,k) = (u(iendu+1,j,k,nnew) - u(iendu+1,j,k,nold))*dts/dt2
      ENDDO
    ENDIF
  END DO
END IF
IF ( .NOT.lperi_y ) THEN
  DO k = 1, ke
    IF (my_cart_neigh(4) == -1) THEN
      DO i = istartv, iendv
       zvbdt_south(i,k) = (v(i,jstartv-1,k,nnew)-v(i,jstartv-1,k,nold))*dts/dt2
      ENDDO
    ENDIF
    IF (my_cart_neigh(2) == -1) THEN
      DO i = istartv, iendv
        zvbdt_north(i,k) = (v(i,jendv+1,k,nnew) - v(i,jendv+1,k,nold))*dts/dt2
      ENDDO
    ENDIF
  ENDDO
ENDIF

!------------------------------------------------------------------------------
! Section 4: Do forecast using the leapfrog / Klemp-Wilhelmson scheme
!------------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Section 4a: Compute the time tendencies due to horizontal advection
  !             and add them to the tendency fields
  !----------------------------------------------------------------------------

  ! The array zwcon is calculated in horizontal_advection and used also in
  ! slow_tendencies.f90. It is defined as a private array for the module and 
  ! is allocated here. Zero boundary conditions are set.

  ALLOCATE ( zwcon(ie,je,ke1), STAT=izstata )
  zwcon(:,:,:) = 0.0_wp
!  zwcon(:,:,ke1) = 0.0_wp

  IF (ltime) CALL get_timings (i_dyn_computations, ntstep, dt, izerror)

  IF (izdebug > 10) THEN
    PRINT *, '        horizontal_advection'
  ENDIF

  IF ( l_euler_dynamics ) THEN
    CALL horizontal_advection
  END IF

  IF (ltime) CALL get_timings (i_horizontal_advection, ntstep, dt, izerror)

  !----------------------------------------------------------------------------
  ! Section 4b: Compute the total time tendencies for the dynamic fields
  !             except the tendencies from sound and gravity wave terms
  !             and calculate the new values of the humidity fields
  !----------------------------------------------------------------------------

  IF (izdebug > 10) THEN
    PRINT *, '        total time tendencies'
  ENDIF

!+++++++++++++++++++ included from complete_tendencies +++++++++++++++++++++++

  ! Add Rayleigh friction to the velocity field.
  ! ---------------------------------------------------------------------

  ! Set the time level for Rayleigh friction and add damping terms

  ! only update tendencies if required
  IF (epsray > 0.0_wp) THEN
    DO  k = 1 , ke
      DO   j = jstartu, jendu
        DO i = istartu, iendu
          utens(i,j,k) = utens(i,j,k) - epsray*u(i,j,k,nold)
        ENDDO
      ENDDO
      DO   j = jstartv, jendv
        DO i = istartv, iendv
          vtens(i,j,k) = vtens(i,j,k) - epsray*v(i,j,k,nold)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! Completition of the vertical wind tendency due to bouyancy
  ! and water loading
  ! ----------------------------------------------------------------

  DO  k = 2, ke
    DO   j = jstart-1, jend+1
      DO i = istart-1, iend+1
        zgrhoc = (dp0(i,j,k-1)*rho0(i,j,k)+dp0(i,j,k)*rho0(i,j,k-1)) /     &
                 (dp0(i,j,k-1)*rho (i,j,k)+dp0(i,j,k)*rho (i,j,k-1)) *g
        zbqb  =  rvd_m_o*qv_now(i,j,k  ) - qc_now(i,j,k  ) - qrs(i,j,k  )
        zbqt  =  rvd_m_o*qv_now(i,j,k-1) - qc_now(i,j,k-1) - qrs(i,j,k-1)
        zbqa  = ( dp0(i,j,k)*zbqt + dp0(i,j,k-1)*zbqb )   &
                 / ( dp0(i,j,k) + dp0(i,j,k-1) )
        wtens(i,j,k) = wtens(i,j,k) + zgrhoc*zbqa
      ENDDO
    ENDDO
  ENDDO

  IF (ldyn_bbc) THEN
    ! wtens(:,:,ke1) is used to store the w-tendency (buoyancy eff.) at full
    ! level ke needed for the dynamical bottom boundary condition
    DO   j = jstart-1, jend+1
      DO i = istart-1, iend+1
        zgrhoc = rho0(i,j,ke)/rho (i,j,ke) *g
        zbqa  = rvd_m_o*qv_now(i,j,ke) - qc_now(i,j,ke) - qrs(i,j,ke  )
        wtens(i,j,ke1) = wtens(i,j,ke1) + zgrhoc*zbqa
      ENDDO
    ENDDO
  ENDIF

!+++++++++++++++++++ included from complete_tendencies +++++++++++++++++++++++

  IF (ltime) CALL get_timings (i_complete_tendencies, ntstep, dt, izerror)

  IF (izdebug > 10) THEN
    PRINT *, '        slow_tendencies'
  ENDIF

  CALL slow_tendencies (zwcon)

  IF (ltime) CALL get_timings (i_slow_tendencies, ntstep, dt, izerror)

  ! Deallocate the private array zwcon
  DEALLOCATE ( zwcon, STAT=izstatd )

  !----------------------------------------------------------------------------
  ! Section 4c: Compute the new values of the dynamic fields by integrating
  !             the sound and gravity wave terms with small time steps using
  !             the split-explicit scheme
  !----------------------------------------------------------------------------

  IF ( l_euler_dynamics ) THEN

    IF (lsemi_imp) THEN
      IF (izdebug > 10) THEN
        PRINT *, '        semi_implicit'
      ENDIF

      CALL semi_implicit (vcoord%kflat)
    ELSE
      IF (izdebug > 10) THEN
        PRINT *, '        fast_waves'
      ENDIF

      ! Pre-multiplication of the slow tendencies by the small timestep
      utens (:,:,:) = utens (:,:,:)*dts
      vtens (:,:,:) = vtens (:,:,:)*dts
      pptens(:,:,:) = pptens(:,:,:)*dts
      ttens (:,:,:) = ttens (:,:,:)*dts
      IF (ldyn_bbc) THEN
        wtens (:,:,2:ke1) = wtens (:,:,2:ke1)*dts
      ELSE
        wtens (:,:,2:ke ) = wtens (:,:,2:ke )*dts
      ENDIF

      CALL fast_waves (ismtstep, dts, xkd, utens, vtens, wtens, ttens, pptens, &
                       zubdt_west, zubdt_east, zvbdt_north, zvbdt_south,       &
                       ie, je, ke, vcoord%kflat )
    ENDIF

    ! Restore the pressure tendency for later use (calculation of omega
    ! in src_output)
    pptens(:,:,:) = pptens(:,:,:)*zdtsr

    IF (lhordiff) THEN
      IF (izdebug > 10) THEN
        PRINT *, '        horizontal_diffusion'
      ENDIF

      CALL comp_hori_diff (itype_hdiff)
    ENDIF

  ELSE

    IF ( my_cart_id == 0 ) THEN
      WRITE(*,*) "WARNING: l_euler_dynamics == .FALSE.!"
    END IF

    ! only Euler forward with physics tendencies:

    u (:,:,:,nnew) = u (:,:,:,nold) + dt2 * utens(:,:,:)
    v (:,:,:,nnew) = v (:,:,:,nold) + dt2 * vtens(:,:,:)
    w (:,:,:,nnew) = w (:,:,:,nold) + dt2 * wtens(:,:,:)
    T (:,:,:,nnew) = T (:,:,:,nold) + dt2 * ttens(:,:,:)
    pp(:,:,:,nnew) = pp(:,:,:,nold) + dt2 * pptens(:,:,:)

  ENDIF

!------------------------------------------------------------------------------
! Section 5: First saturation adjustment for t, qv, qc at time level nnew
!------------------------------------------------------------------------------

IF (lcond) THEN

  IF (itype_gscp < 100 .OR. l2mom_satads) THEN

    IF ( loutput_diab ) THEN
      ! prepare storing of the latent heat release by the following saturation adjustment:
      ttens_diab(:,:,:) = ttens_diab(:,:,:) - t(:,:,:,nnew) / dt
    ENDIF

    IF (izdebug > 10) THEN
      PRINT *, '        saturation adjustment'
    ENDIF

#ifdef NUDGING
    IF (llhn) THEN
      ! For calculation of latent heating rate T has to be stored before 
      ! saturation adjustment
      IF (itype_gscp >= 3) THEN
        ! in case of itype_gscp >= 3 (hydci or hydci_pp) get latent heat rate of 
        ! the last time step and store it on the next time step. Both Routines 
        ! are runnig at the end of the time step and therefore after the latent
        ! heat nudging.
        tt_lheat(:,:,:,nnew)=tt_lheat(:,:,:,nnow)
      ENDIF

      CALL get_gs_lheating ('add', 1, ke)
    ENDIF
#endif

    ! The array zwcon is allocated again here but just as intermediate storage
    ! for routine satad.
    ALLOCATE (zwcon(ie,je,9), STAT=izstata )

    DO  k = 1, ke
      DO j = jstartpar, jendpar
        DO i = istartpar, iendpar
          zwcon(i,j,1) =  p0(i,j,k) + pp(i,j,k,nnew)
        ENDDO
      ENDDO
      kitpro = 1
      IF (zuvwmax(k) >=  2.0_wp) kitpro = 2
      IF (zuvwmax(k) >= 10.0_wp) kitpro = 3

      CALL satad ( kitpro, t(:,:,k,nnew), qv(:,:,k),                       &
                   qc(:,:,k), t(:,:,k,nnow), zwcon(:,:,1),                 &
                   zwcon(:,:,2), zwcon(:,:,3), zwcon(:,:,4), zwcon(:,:,5), &
                   zwcon(:,:,6), zwcon(:,:,7), zwcon(:,:,8), zwcon(:,:,9), &
                   b1, b2w, b3, b4w, b234w, rdv, o_m_rdv,                  &
                   rvd_m_o, lh_v, cpdr, cp_d,                              &
                   ie, je, istartpar, iendpar, jstartpar, jendpar  )
    ENDDO

    DEALLOCATE (zwcon, STAT=izstatd)

#ifdef NUDGING
    ! calculate gridscale latent heating rate from the saturation adjustment
    IF (llhn) CALL get_gs_lheating ('inc', 1, ke)
#endif

    IF ( loutput_diab ) THEN
      ! complete the storing of the latent heat release by the above saturation adjustment:
      ttens_diab(:,:,:) = ttens_diab(:,:,:) + t(:,:,:,nnew) / dt
    ENDIF

  ENDIF

ENDIF

!------------------------------------------------------------------------------
! End of module procedure org_leapfrog  
!------------------------------------------------------------------------------

END SUBROUTINE org_leapfrog

!==============================================================================
!+ Module procedure in "src_leapfrog" for computing horizontal advection
!------------------------------------------------------------------------------

SUBROUTINE horizontal_advection

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "src_leapfrog" computes the time tendencies
!   of all prognostic variables due to horizontal advection. The advective
!   tendencies are then added to the tendency fields provided for each 
!   variable.
!
! Method:
!   Horizontal advection is written in advection form and discretizised
!   with second order centred differences on an Arakawa C-grid.
!   For u and v, the potential vorticity/kinetic energy form is used.
!
!------------------------------------------------------------------------------

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    i,  j,  k, istat, izerror, iznumsl, iznumpd2dt, iza, iztrcr

  REAL    (KIND=wp)        ::  &
    z2 , z3 ,            & !
    zdpq,                & !
    zuadv , zvadv ,      & !
    zgvq  , zguq  , zpvq,& !
    zzlam , zzphi ,      & !
    zws                    !
   
! Local (automatic) arrays:
! ------------
  REAL    (KIND=wp)        ::  &
    zgu   (ie,je,ke),    & !
    zgv   (ie,je,ke),    & !
    zdphr (ie,je,ke),    & !
    zdp0u (ie,je,ke),    & !
    zdp0v (ie,je,ke),    & !
    zke   (ie,je   ),    & !
    zpv   (ie,je   ),    & !
    zphl  (ie,je   ),    & !
    zvs   (ie,je,ke),    & !
    zrdx  (   je   ),    & !
    zrdy  (   je   ),    & !
    zfadsx(   je   ),    & !
    zfadsy(   je   ),    & !
    zfadux(   je   ),    & !
    zfadvx(   je   ),    & !
    zfpvoy(   je   ),    & !
    zfpvuy(   je   ),    & !
    zfkeo (   je   ),    & !
    zfkeu (   je   )       !

  REAL    (KIND=wp)     ::  &
    ut    (ie,je,ke),    & !
    vt    (ie,je,ke)       !

  INTEGER (KIND=iintegers) :: kzdims(24)

  INTEGER (KIND=iintegers) :: &
    izadv   (trcr_get_ntrcr()) , &
    izsp_adv_lf (trcr_get_ntrcr())

  INTEGER (KIND=iintegers), ALLOCATABLE :: btrj_idx (:,:,:,:)
  REAL    (KIND=wp)       , ALLOCATABLE :: btrj_wght(:,:,:,:)
  REAL    (KIND=wp)       , ALLOCATABLE :: zw_ke(:,:)
  REAL    (KIND=wp)       , ALLOCATABLE :: ztrwt(:,:,:)
  REAL    (KIND=wp)       , ALLOCATABLE :: zprecwt(:,:,:)

  CHARACTER (LEN=255) ::  yzerrmsg

  CHARACTER (LEN=25)  ::  yzroutine

! Tracer pointers:
!-----------------
  REAL (KIND=wp),     POINTER   :: &
    qv        (:,:,:) => NULL(),    & ! QV at nnow
    ztrcr_old (:,:,:) => NULL(),    & ! tracer variable at nold
    ztrcr_new (:,:,:) => NULL(),    & ! tracer variable at nnew
    ztrcr_now (:,:,:) => NULL(),    & ! tracer variable at nnow
    ztrcr_tens(:,:,:) => NULL()       ! tracer tendency variable

!
! End of header 
!==============================================================================
 
!------------------------------------------------------------------------------
! Begin Subroutine horizontal_advection
!------------------------------------------------------------------------------
 
  kzdims(:) = 0_iintegers

  yzroutine = 'horizontal_advection'

!------------------------------------------------------------------------------
! Section 0: Define some tracer properties and calculate some coefficients
!            accordingly
!------------------------------------------------------------------------------

  iznumsl    = 0_iintegers  ! number of tracers being advected using Semi-Lagr.
  iznumpd2dt = 0_iintegers  ! number of tracers being advected using pos. def.
                            ! scheme over 2 timelevels

  ! Retrieval of the metadata
  CALL trcr_meta_get(izerror, 'SP_ADV_LF', izsp_adv_lf)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_meta_get(izerror, T_ADV_ID, izadv)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  ! Definition of the number of tracers having special advection
  ! (pd over 2 timesteps, pd over 1 timestep, semi-Lagrange)
  DO iza = 1, trcr_get_ntrcr ()
    IF ( izsp_adv_lf (iza) == 2_iintegers ) THEN
      iznumpd2dt = iznumpd2dt + 1
    ELSEIF ( izsp_adv_lf (iza) == 3_iintegers ) THEN
      iznumsl = iznumsl + 1
    ENDIF
  ENDDO

  ! transformed velocity required for positive definite advection
  IF ( iznumpd2dt > 0_iintegers ) THEN
    DO  k = 1, ke
      DO j = 1, je
        DO i = 1, ie
          zvs(i,j,k) = v(i,j,k,nnow)*crlat(j,2)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! tmp variable required for pd over 2 timesteps
  IF ( iznumpd2dt > 0_iintegers ) THEN
    ALLOCATE ( ztrwt (ie,je,ke), STAT = izerror )
  ENDIF

!------------------------------------------------------------------------------
! Section 1:  Setup of metrical factors and interpolation coefficients
!------------------------------------------------------------------------------
   
  ! Metrical factors for horizontal discretization  

  IF (linit_fields) THEN
    zfadsx (:) = 0.0_wp
    zfadsy (:) = 0.0_wp
    zgu(:,:,:) = 0.0_wp
    zgv(:,:,:) = 0.0_wp
  ENDIF

  DO  j = jstart-1 , jend+1
    zfadsx(j) = 0.5_wp*acrlat(j,1)*eddlon
    zfadsy(j) = 0.5_wp*acrlat(j,1)*eddlat
    zfadux(j) = acrlat(j,1)*eddlon
    zfadvx(j) = acrlat(j,2)*eddlon
    zfpvoy(j) = acrlat(j,2)*eddlat*crlat(j+1,1)
    zfpvuy(j) = acrlat(j,2)*eddlat*crlat(j  ,1)
    zfkeo (j) = crlat(j  ,2)/crlat(j,1)
    zfkeu (j) = crlat(j-1,2)/crlat(j,1)
  ENDDO          
  zrdx(:) = acrlat(:,1)*eddlon
  zrdy(:) = acrlat(:,1)*eddlat

!------------------------------------------------------------------------------
! Section 2:  Calculation of horizontal advection on full levels using
!             centered differences O(2)
!------------------------------------------------------------------------------

loop_over_full_levels:  DO  k = 1, ke

  ! Calculation of mass weighted wind velocities (zgu, zgv), potential
  ! vorticity and kinetic energy 
  ! ----------------------------------------------------

  DO    j = jstart-2 , je-1
    DO  i = istart-2 , ie-1
      zdp0u(i,j,k) =  0.5_wp*( dp0(i,j,k) + dp0(i+1,j  ,k) )
      zdp0v(i,j,k) =  0.5_wp*( dp0(i,j,k) + dp0(i  ,j+1,k) )
      zgv(i,j,k)   =  zdp0v(i,j,k) * v(i,j,k,nnow)*crlat(j,2)
      zgu(i,j,k)   =  zdp0u(i,j,k) * u(i,j,k,nnow)                 
    ENDDO      
  ENDDO       
  DO  i = istart-1 , ie
    zdp0u(i,je,k) = dp0(i,je,k)
    zdp0v(i,je,k) = dp0(i,je,k)
  ENDDO       
  DO  j = jstart-1 , je
    zdp0u(ie,j,k) = dp0(ie,j,k)
    zdp0v(ie,j,k) = dp0(ie,j,k)
  ENDDO       

  IF (lmetr) THEN
    DO    j = jstart-1, jend+1
      DO  i = istart-1, iend+1
        z2           = zfadvx(j)*( v(i+1,j,k,nnow)-v(i,j,k,nnow) )
        z3           = zfpvoy(j)*u(i,j+1,k,nnow) - zfpvuy(j)*u(i,j,k,nnow)
        zdpq         = 0.5_wp*(zdp0u(i,j,k) + zdp0u(i,j+1,k))
        zpv(i,j)     = ( fc(i,j) + z2 - z3 ) / zdpq
        zke(i,j)     = 0.25_wp*(  u(i-1,j,k,nnow)**2 + u(i,j,k,nnow)**2  &
                              + zfkeo(j)*v(i,j  ,k,nnow)**2              &
                              + zfkeu(j)*v(i,j-1,k,nnow)**2 )
      ENDDO
    ENDDO
  ELSE
    DO    j = jstart-1, jend+1
      DO  i = istart-1, iend+1
        z2           = zfadvx(j)*( v(i+1,j,k,nnow)-v(i,j,k,nnow) )
        z3           = eddlat /r_earth * (u(i,j+1,k,nnow) - u(i,j,k,nnow) )
        zdpq         = 0.5_wp*(zdp0u(i,j,k) + zdp0u(i,j+1,k))
        zpv(i,j)     = ( fc(i,j) + z2 - z3 ) / zdpq
        zke(i,j)     = 0.25_wp*(  u(i-1,j,k,nnow)**2 + u(i,j,k,nnow)**2  &
                              + zfkeo(j)*v(i,j  ,k,nnow)**2              &
                              + zfkeu(j)*v(i,j-1,k,nnow)**2 )
      ENDDO
    ENDDO
  ENDIF

  ! Advective tendencies for the horizontal wind components (u, v)
  ! --------------------------------------------------------------

  DO    j = jstartu, jendu 
    DO  i = istartu, iendu
      zgvq  =  zgv(i,j,k) + zgv(i+1,j,k) + zgv(i,j-1,k) + zgv(i+1,j-1,k)
      zpvq  = (zpv(i,j)   + zpv(i,j-1)  ) / crlat(j,1)                    
      zuadv = - zfadux(j)*( zke(i+1,j)   - zke(i,j)   ) + 0.125_wp*zpvq*zgvq  
      utens(i,j,k) = utens(i,j,k) + zuadv
    ENDDO   
  ENDDO       
  DO    j = jstartv, jendv 
    DO  i = istartv, iendv
      zguq  =  zgu(i-1,j,k) + zgu(i,j,k) + zgu(i-1,j+1,k) + zgu(i,j+1,k)
      zpvq  =  zpv(i-1,j)   + zpv(i,j)      
      zvadv = - edadlat*( zke(i,j+1)   - zke(i,j)   ) - 0.125_wp*zpvq*zguq
      vtens(i,j,k) = vtens(i,j,k) + zvadv
    ENDDO   
  ENDDO       

  ! advective tendencies for prognostic scalars (t, pp, w)
  ! ------------------------------------------------------

  zdphr (:,:,k) = 1.0_wp/dp0(:,:,k)
  dqvdt (:,:,k) = 0.0_wp

! advective tendency for w at the full level ke
! needed to the dynamical bottom boundary condition
! ----------------------------------------------------

  IF (ldyn_bbc .AND. (k == ke)) THEN
    ALLOCATE (zw_ke(ie,je))
    DO   j = jstart-2, jend+2
      DO i = istart-2, iend+2
        zw_ke(i,j)= 0.5_wp * (w (i,j,k,nnow) + w (i,j,k+1,nnow))
      ENDDO
    ENDDO

    ! wtens(:,:,ke1) is used to store the w-tendency (buoyancy eff.) at full
    ! level ke needed for the dynamical bottom boundary condition
    wtens(:,:,ke1) = 0.0_wp

    CALL hadv_cd2 ( zw_ke(:,:), wtens(:,:,ke1), zgu(:,:,k), zgv(:,:,k),       &
                    zdphr(:,:,k), zfadsx, zfadsy, ie, je, istart-1, iend+1,   &
                    jstart-1, jend+1)
    DEALLOCATE (zw_ke)
  ENDIF

  CALL hadv_cd2 ( pp(:,:,k,nnow), pptens(:,:,k), zgu(:,:,k), zgv(:,:,k),      &
                  zdphr(:,:,k), zfadsx, zfadsy, ie, je, istart, iend, jstart, &
                  jend )
  CALL hadv_cd2 ( t (:,:,k,nnow), ttens (:,:,k), zgu(:,:,k), zgv(:,:,k),      &
                  zdphr(:,:,k), zfadsx, zfadsy, ie, je, istart, iend, jstart, &
                  jend )

ENDDO loop_over_full_levels

! advective  tendency for tracers for hadv_cd2
! --------------------------------------------

! loop over tracers
DO iztrcr = 1, trcr_get_ntrcr()

  ! do advection?
  IF ( izadv(iztrcr) == T_ADV_ON ) THEN

    ! advect using centered differences
    IF ( izsp_adv_lf(iztrcr) == 1_iintegers ) THEN

      ! retrieve pointer to tracer (at nnow)
      CALL trcr_get(izerror, iztrcr, ptr_tlev=nnow, ptr=ztrcr_now,            &
                    ptr_tens=ztrcr_tens)
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF

      ! compute horiz. advection
      DO k = 1, ke
        CALL  hadv_cd2 ( ztrcr_now(:,:,k),                                    &
                         ztrcr_tens(:,:,k), zgu(:,:,k), zgv(:,:,k),           &
                         zdphr(:,:,k),zfadsx, zfadsy, ie, je, istart,         &
                         iend, jstart, jend )
      ENDDO

    ENDIF

  ENDIF 

ENDDO

! treat special case of QV: save the advective tendency also in dqvdt
IF ( izadv(idt_qv) == T_ADV_ON ) THEN
  CALL trcr_get(izerror, idt_qv, ptr_tlev=nnow, ptr=qv)
  DO k = 1, ke
    CALL  hadv_cd2 ( qv(:,:,k),                                               &
                     dqvdt(:,:,k), zgu(:,:,k), zgv(:,:,k),                    &
                     zdphr(:,:,k),zfadsx, zfadsy, ie, je, istart,             &
                     iend, jstart, jend )
  ENDDO
ENDIF

!------------------------------------------------------------------------------
! Section 3:  Calculation of horizontal advection on half levels
!------------------------------------------------------------------------------

loop_over_half_levels:  DO  k = 2, ke

  ! First, interpolate the mass-weighted velocities to half levels
  DO j = jstart-1, jend+1
  DO i = istart-1, iend+1
    zgu (i,j,k) = u(i,j,k,nnow)*zdp0u(i,j,k-1) +  u(i,j,k-1,nnow)*zdp0u(i,j,k)
    zgv (i,j,k) = v(i,j,k,nnow)*zdp0v(i,j,k-1) +  v(i,j,k-1,nnow)*zdp0v(i,j,k)
    zgv (i,j,k) = zgv(i,j,k)  *crlat(j,2)
  ENDDO
  ENDDO
  zdphr(:,:,k) = 1.0_wp/(dp0(:,:,k) + dp0(:,:,k-1))

  ! 3.1 Advective tendencies for the vertical wind component (w)
  ! --------------------------------------------------------------

  CALL hadv_cd2 ( w(:,:,k,nnow), wtens(:,:,k), zgu(:,:,k), zgv(:,:,k),        &
                  zdphr(:,:,k), zfadsx, zfadsy, ie, je, istart, iend, jstart, &
                  jend )


  ! 3.2 Calculate the mass weighted contravariant vertical velocity zwcon
  ! --------------------------------------------------------------

  DO    j = jstart, jend+1
    DO  i = istart, iend+1
      zws  =  dp0(i,j,k-1)*rho0(i,j,k) + dp0(i,j,k)*rho0(i,j,k-1)
      zwcon(i,j,k) = - g*w(i,j,k,nnow)*zws*zdphr(i,j,k)
    ENDDO
  ENDDO

  IF ( k > vcoord%kflat ) THEN
    zphl(:,:) = p0hl(:,:,k)
    DO j = jstart-1, jend+1
      DO i = istart-1, iend+1
        zgu(i,j,k) = zgu(i,j,k)/((zdp0u(i,j,k)+zdp0u(i,j,k-1))**2)
        zgv(i,j,k) = zgv(i,j,k)/((zdp0v(i,j,k)+zdp0v(i,j,k-1))**2)
        zgu(i,j,k) = zgu(i,j,k)*(zphl(i+1,j) - zphl(i,j))
        zgv(i,j,k) = zgv(i,j,k)*(zphl(i,j+1) - zphl(i,j))
      ENDDO
    ENDDO
    DO    j = jstart, jend+1
      DO  i = istart, iend+1
        zzlam = zfadsx(j)*( zgu(i,j,k) + zgu(i-1,j,k) )
        zzphi = zfadsy(j)*( zgv(i,j,k) + zgv(i,j-1,k) )
        zwcon(i,j,k) = zwcon(i,j,k) - (zzlam + zzphi)*(dp0(i,j,k)+dp0(i,j,k-1))
      ENDDO     
    ENDDO       
  ENDIF

ENDDO loop_over_half_levels
 
! For periodic boundary conditions, the diagnostic variable zwcon doesn't have
! to be exchanged. For 2-D runs, we copy zwcon onto the y-boundaries.
IF ( l2dim ) THEN
  DO k = 2, ke
     DO j = 1,nboundlines
        zwcon (:,jstart-j,k) = zwcon (:,jend  +1-j,k)
        zwcon (:,jend  +j,k) = zwcon (:,jstart-1+j,k)
     ENDDO
  ENDDO
ENDIF

!------------------------------------------------------------------------------
! Section 4:  Calculation of horizontal advection on full levels using
!             an Euler time stepping (positive definite)
!------------------------------------------------------------------------------

! advective  tendency for tracers for hadv_pd
! --------------------------------------------

! loop over tracers
DO iztrcr = 1, trcr_get_ntrcr()

  ! do advection?
  IF ( izadv(iztrcr) == T_ADV_ON ) THEN

    ! Advect using positive definite scheme over 2 time steps
    IF ( izsp_adv_lf(iztrcr) == 2_iintegers ) THEN

      ! retrieve pointer to tracer (at nnew)
      CALL trcr_get(izerror, iztrcr, ptr_tlev=nnew, ptr=ztrcr_new)
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF

      ! retrieve pointer to tracer (at nold)
      CALL trcr_get(izerror, iztrcr, ptr_tlev=nold, ptr=ztrcr_old)
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF

      ! positive advection is done in two steps (first advancing from
      ! old->now, then advancing from now->new

      ! STEP 1 (old -> now)

      ! loop over vertical layers
      DO  k = 1, ke

        ! set boundary values (by copying in full field)
        ztrwt(:,:,k) = ztrcr_new(:,:,k)

        ! horizontal avection (result stored in ztrwt)
        CALL hadv_pd( u(:,:,k,nnow), zvs(:,:,k), ztrcr_old(:,:,k),         &
                      ztrwt(:,:,k) , zrdx(:), zrdy(:), dt,                  &
                      ie, je, istart, iend, jstart, jend )

      ENDDO

      ! vertical advection (result in-place ztrwt)
      CALL vadv_pd( zwcon(:,:,:), dp0(:,:,:) , ztrwt(:,:,:),                  &
                  dt, ie, je, ke, ke1, istart, iend, jstart, jend )

      ! perform halo-update of ztrwt before next advection step
      IF (num_compute > 1 ) THEN
        IF (ltime)                                                            &
          CALL get_timings (i_horizontal_advection, ntstep, dt, izerror)
        kzdims(1:20)=(/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                 &
         ( 4,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,ie,   &
          je, kzdims, jstartpar, jendpar, 2, nboundlines, my_cart_neigh,      &
          lperi_x, lperi_y, l2dim,                                            &
          20000+nexch_tag,ldatatypes, ncomm_type, izerror, yzerrmsg,          &
          ztrwt(:,:,:) )
        IF (ltime) CALL get_timings (i_communications_dyn, ntstep, dt, izerror)
      ENDIF

      ! STEP 2 (now -> new)

      ! loop over vertical layers
      DO k = 1, ke

        ! horizontal advection (result stored in ztrcr_new)
        CALL hadv_pd( u (:,:,k,nnow), zvs(:,:,k), ztrwt(:,:,k),             &
                    ztrcr_new(:,:,k), zrdx(:), zrdy(:), dt,                 &
                    ie, je, istart, iend, jstart, jend )

      ENDDO

      ! vertical advection (result in-place ztrcr_new)
      CALL vadv_pd( zwcon(:,:,:), dp0(:,:,:) , ztrcr_new(:,:,:),              &
                dt, ie, je, ke, ke1, istart, iend, jstart, jend )

    ENDIF ! izsp_adv_lf = 2
  ENDIF ! advection on?

ENDDO ! loop over tracers


IF ( ALLOCATED(ztrwt) ) THEN
  DEALLOCATE(ztrwt)
ENDIF

!------------------------------------------------------------------------------
! Section 5: Calculate transformed velocity components for SL advection
!------------------------------------------------------------------------------

IF (iznumsl > 0_iintegers) THEN
  DO k=1, ke
    DO j=1, je
      DO i=1, ie
         ut(i,j,k) = acrlat(j,1) * u(i,j,k,nnow)
         vt(i,j,k) = v(i,j,k,nnow) / r_earth
      ENDDO
    ENDDO
  ENDDO

  ALLOCATE (zprecwt(ie,je,ke+1), STAT=istat)
  zprecwt(:,:,:) = 0.0_wp

  DO k = 2, ke
    DO j = jstart, jend
      DO i = istart, iend

         ! vertical velocity only by air movement
         zprecwt(i,j,k) = zwcon(i,j,k) / ( p0(i,j,k) - p0(i,j,k-1) )

         ! with sedimentation velocity for rain
         ! zprecwt(i,j,k,2) = zwcon(i,j,k) / ( p0(i,j,k) - p0(i,j,k-1) ) &
         !             + g * rho0(i,j,k) / ( p0(i,j,k)-p0(i,j,k-1) ) * &
         !               vSediR(i,j,k)
      ENDDO
    ENDDO
  ENDDO

  IF ( num_compute > 1  ) THEN
    kzdims(1:24)=(/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                     &
       ( 4,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,&
        kzdims, jstartpar, jendpar, 2, nboundlines, my_cart_neigh,            &
        lperi_x, lperi_y, l2dim,                                              &
        20000+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,           &
        zprecwt(:,:,:) )
  ENDIF
ENDIF

ALLOCATE ( btrj_wght(ie, je, ke, 3), stat=istat )
ALLOCATE ( btrj_idx (ie, je, ke, 3), stat=istat )
btrj_wght(:,:,:,:) = 0.0_wp
btrj_idx (:,:,:,:) = 0_iintegers

CALL backtraj_trilin_dt2_3tl(ut,vt,zprecwt(:,:,:), eddlon, eddlat,            &
                             1.0_wp, dt, ie, je, ke, istart, iend,        &
                             jstart, jend, btrj_idx, btrj_wght )

DEALLOCATE( zprecwt )

!------------------------------------------------------------------------------
! Section 6:  Calculation of horizontal advection on full levels using
!             Semi-Lagrange advection
!------------------------------------------------------------------------------

  ! loop over tracers
  DO iztrcr = 1, trcr_get_ntrcr()

    ! advection on?
    IF ( izadv(iztrcr) == T_ADV_ON ) THEN

      ! advect with semi-Lagrange scheme?
      IF ( izsp_adv_lf(iztrcr) == 3_iintegers ) THEN

        ! get pointers to tracer (at nold and nnew)
        CALL trcr_get(izerror, iztrcr, ptr_tlev=nold, ptr=ztrcr_old)
        IF (izerror /= 0) THEN
          yzerrmsg = trcr_errorstr(izerror)
          CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
        ENDIF
        CALL trcr_get(izerror, iztrcr, ptr_tlev=nnew, ptr=ztrcr_new)
        IF (izerror /= 0) THEN
          yzerrmsg = trcr_errorstr(izerror)
          CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
        ENDIF

        ! semi-Lagrange advection (trilinear interpolation)
        CALL interpol_sl_trilin( ztrcr_old(:,:,:), ztrcr_new(:,:,:), btrj_idx,&
                                 btrj_wght, dt, ie, je, ke, istart, iend,     &
                                 jstart, jend )

      ENDIF ! izsp_adv_lf = 4?
    ENDIF ! advection on?

  ENDDO ! loop over tracers

  IF (iznumsl > 0_iintegers) THEN
    DEALLOCATE( btrj_idx, btrj_wght )
  ENDIF

!------------------------------------------------------------------------------
! End of module procedure "horizontal advection"
!------------------------------------------------------------------------------

END SUBROUTINE horizontal_advection

!==============================================================================

END MODULE src_leapfrog
