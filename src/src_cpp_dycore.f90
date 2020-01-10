!+ Source module for C++ dycore (rewrite of 2-timelevel Runge-Kutta dynamics)
!------------------------------------------------------------------------------

! uncomment the following define in order to debug the C++ Dycore
!#define CPP_DEBUG

MODULE src_cpp_dycore

#ifdef CPP_DYCORE

!------------------------------------------------------------------------------
!
! Description:
!
! Current Code Owner: C2SM, Pascal Spoerri
!  phone:  +41 58 460 9436
!  email:  pascal.spoerri@meteoswiss.ch
!
! History:
! Version      Date       Name
! ----------   ---------- ----
! V5_3         2015-10-09 Oliver Fuhrer
!  Initial release
! V5_4         2016-03-10 Oliver Fuhrer, Pascal Spoerri, Xavier Lapillonne
!  Added the radius of the earth to the constants that are given to the dycore.
!  The dycore will check if the constants match with its own constant. 
!  We also increased the dycore version for this change.
!  Added the cold pool diffusion threshold (thresh_cold_pool) to the dycore.
!  Added simple clipping for SL3 advection scheme.
!  Update of serialization statements
! V5_4a        2016-05-10 Pascal Spoerri
!  Replacement of t_s with t_g in the dycore wrapper interface. Adapted the 
!  dycore_verify_pointers interface accordingly. Increasing the 
!  dycoreInterfaceVersion to 11.
! V5_4b        2016-07-12 Pascal Spoerri
!  Use imode_turb from turb_data
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

USE kind_parameters , ONLY :   &
  wp           ! KIND-type parameters for real variables

USE data_parallel, ONLY :      &
  nprocx, nprocy,              &
  my_cart_neigh,               &
  my_cart_pos,                 &
  my_cart_id,                  &
  my_world_id,                 &
  num_compute,                 &
  icomm_compute,               &
  isubpos

USE data_io, ONLY :            &
    llb_qi,       & ! if .TRUE., take qi_bd-values from lateral boundaries file
                    ! else, qi_bd is set in the model
    llb_qr_qs,    & ! if .TRUE., take qr_bd- and qs_bd-values from lateral
                    ! bound. file else, qr_bd and qs_bd are set in the model
    llb_qg          ! if .TRUE., take qg_bd-values from lateral boundaries file
                    ! else, qg_bd is set in the model

USE environment, ONLY : model_abort

USE data_tracer, ONLY : T_MISSING,     T_NAME_ID,          &
    T_ADV_ID, T_DIFF_ID, T_TURB_ID, T_CONV_ID, T_LBC_ID,   &
    T_BBC_ID, T_RELAX_ID, T_DAMP_ID, T_CLP_ID
    !PS TODO Integrate SEDIM and SPPT into tracers
    !PS T_SEDIM_ID, T_SPPTPERT_ID

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

  ! 2. horizontal and vertical sizes of the fields and related variables
  ! --------------------------------------------------------------------
  ie,           & ! number of grid points in zonal direction
  je,           & ! number of grid points in meridional direction
  ke,           & ! number of grid points in vertical direction
  ke1,          & ! KE+1

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

  dlon,         & ! grid point distance in zonal direction (in degrees)
  dlat,         & ! grid point distance in meridional direction (in degrees)

  ! 5. variables for the time discretization and related variables
  ! --------------------------------------------------------------

  dt,           & ! long time-step

! 6. variables for the reference atmosphere
! -----------------------------------------

  ie_tot,       & ! number of grid points in zonal direction
  je_tot          ! number of grid points in meridional direction

! end of data_modelconfig

!------------------------------------------------------------------------------

USE vgrid_refatm_utils,       ONLY :  &
  vcoord

! end of vgrid_refatm_utils

!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

  ! 2. physical constants and related variables
  ! -------------------------------------------
  r_earth,      & ! mean radius of the earth (m)
  r_d,          & ! gas constant for dry air
  cp_d,         & ! specific heat of dry air at constant pressure
  g               ! acceleration due to gravity


! end of data_constants

!------------------------------------------------------------------------------
  
USE data_fields     , ONLY :   &

  ! 1. constant fields for the reference atmosphere                 (unit)
  ! -----------------------------------------------
  rho0       ,    & ! reference density at the full model levels    (kg/m3)
  p0         ,    & ! reference pressure at main levels             ( Pa  )
  dt0dz      ,    & ! temperature gradient of reference atmosphere  ( K/m )
  t0         ,    & ! reference temperature                         ( K   )
  p0hl,         & ! reference pressure at half levels           ( Pa  )
  t0hl,         & ! reference pressure at half levels           ( Pa  )
  hd_mask,      & ! 3D-domain mask for horizontal diffusion         --
                    ! vertical   turbulent diffusion coefficients
  tinc_lh,      & ! temperature increment due to latent heat      (  K  )
  tkvm,         & ! ... for momentum                              (m2/s)
  tkvh,         & ! ... for heat and moisture                     (m2/s)
  fc,           & ! coriolis-parameter                            ( 1/s )
!   turbulent coefficients at the surface
  tcm,          & ! transfer coefficient for momentum             ( -- )
  tch,          & ! transfer coefficient for heat and moisture    ( -- )
  u,            & ! zonal wind speed                              ( m/s )
  v,            & ! meridional wind speed                         ( m/s )
  w,            & ! vertical wind speed (defined on half levels)  ( m/s )
  t,            & ! temperature                                   (  k  )
  pp,           & ! deviation from the reference pressure         ( pa  )
  rho,          & ! total density of air                          (kg/m3)
  utens,        & ! u-tendency without sound-wave terms           ( m/s2)
  vtens,        & ! v-tendency without sound-wave terms           ( m/s2)
  wtens,        & ! w-tendency without sound-wave terms           ( m/s2)
  ttens,        & ! t-tendency without sound-wave terms           ( K/s )
  pptens,       & ! pp-tendency without sound-wave terms          (Pa/s )
  qrs,          & ! precipitation water content (water loading)   (kg/kg)
  dqvdt,        & ! threedimensional moisture convergence         ( 1/s )
  qvsflx,       & ! surface flux of water vapour                  (kg/m2s)
  umfl_s,       & ! u-momentum flux (surface)                     ( N/m2)
  vmfl_s,       & ! v-momentum flux (surface)                     ( N/m2)
  shfl_s,       & ! sensible heat flux (surface)                  ( W/m2)
  lhfl_s,       & ! latent heat flux (surface)                    ( W/m2)
  ps,           & ! surface pressure                           ( pa  )
  t_s,          & ! temperature of the ground surface (soil)   (  K  )
  t_g,          & ! weighted surface temperature               (  K  )
  a1t,          & !                                               ( -- )
  a2t,          & !   
  crlat,        & ! cosine of transformed latitude                  --
  acrlat,       & ! 1 / ( crlat * radius of the earth )           ( 1/m )
  tgrlat,       & ! tangens of transformed latitude             
  hhl,          & ! geometical height of half model levels        ( m )
  rdcoef,       & ! Rayleigh damping coefficient
  rmy,          & ! Davis-parameter for relaxation (mass, qv, qc)   --
  rmyq,         & ! Davis-parameter for relaxation (qr, qs, qg)     --
  u_bd,         & ! boundary field for u                        ( m/s )
  v_bd,         & ! boundary field for v                        ( m/s )
  t_bd,         & ! boundary field for t                        (  K  )
  pp_bd           ! boundary field for pp                       (  pa )

! end of data_fields

!------------------------------------------------------------------------------

USE src_stoch_physics,ONLY : &
  pertstoph       ! stochastic multiplier of physics tendencies   ( -- )

! end of src_stoch_physics

!------------------------------------------------------------------------------

USE grid_metrics_utilities, ONLY: &
  sqrtg_r_s,    & ! reciprocal square root of G at skalar points  ( 1/m )
  sqrtg_r_u,    & ! reciprocal square root of G at u points       ( 1/m )
  sqrtg_r_v,    & ! reciprocal square root of G at v points       ( 1/m )
  sqrtg_r_w       ! reciprocal square root of G at w points       ( 1/m )

! end of grid_metrics_utilities

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

  ! 1. start and end of the forecast
  ! --------------------------------
  nstart,        & ! first time step of the forecast
  ntstep,        & ! actual time step
  ! indices for permutation of three time levels
  nnow,          & ! corresponds to ntstep
  nnew,          & ! corresponds to ntstep + 1

  ! 3. controlling the physics
  ! --------------------------
  itype_gscp,    & ! type of microphys. parametrization

  ! 4. controlling the dynamics
  ! ---------------------------
  lcpp_dycore,  & ! whether to use C++ dycore on not
  ltadv_limiter,& ! use limiter for temperature advection (itheta_adv=2 only)
  ldyn_bbc,     & ! dynamical bottom boundary condition
  itype_spubc,  & ! type of Rayleigh damping in the upper levels

  ! 7. additional control variables
  ! -------------------------------
  y_scalar_advect,  & ! type of scalar advection scheme
                      ! "SL3_SC", "SL3_MF", "SL3_SFD", "Bott2", "Bott4"
                      ! "Bott2_Strang", "Bott4_Strang", "vanLeer", "PPM"
  iadv_order,       & ! order of the horizontal advection scheme for dynamics
  ltime,            & ! detailed timings of the program are given
  itheta_adv,       & ! =0: use T' (perturbation temperature) for advection
                      ! =1: use theta' (perturbation potential temperature)
                      ! =2: use theta (full potential temperature)
  ltur,             & ! forecast with vertical diffusion
  lcori,            & ! Coriolis force
  ldiabf_lh,        & ! include diabatic forcing due to latent heat in RK-scheme
  lhordiff,         & ! running with horizontal diffusion
  lperi_x,          & ! switch for periodicity in x direction
  lperi_y,          & ! switch for periodicity in y direction
  lsppt,            & ! switch for stochastic perturbation of parametrization tendencies (SPPT)
  itype_fast_waves, & ! to select a Runge-Kutta scheme
  divdamp_slope,    & ! exceed the theoretical slope stability criterion of
                      ! the divergence damping (only for itype_fast_waves=2)
  irunge_kutta,     & ! to select a Runge-Kutta scheme
  irk_order,        & ! order of runge-kutta scheme
  xkd,              & ! coefficient for divergence damping
  itype_hdiff,      & ! type of horizontal diffusion (=1: 4th order linear),
  hd_corr_u_bd,     & ! correction factor for horizontal diffusion flux of u,v,w in boundary zone
  hd_corr_t_bd,     & ! correction factor for horizontal diffusion flux of t in boundary zone
  hd_corr_trcr_bd,  & ! correction factor for horizontal diffusion fluc of qv,qc in boundary zone
  hd_corr_p_bd,     & ! correction factor for horizontal diffusion fluc of p in boundary zone
  hd_corr_u_in,     & ! correction factor for horizontal diffusion flux of u,v,w in domain
  hd_corr_t_in,     & ! correction factor for horizontal diffusion flux of t in domain
  hd_corr_trcr_in,  & ! correction factor for horizontal diffusion fluc of qv,qc in domain
  hd_corr_p_in,     & ! correction factor for horizontal diffusion fluc of p in domain
  hd_dhmax,         & ! maximum gridpoint height difference for applying
                      ! horizontal diffusion fluxes between them
  l_diff_smag,      & ! use Smagorinsky-diffusion for u and v
  l_diff_cold_pools,& ! use targeted diffusion for cold pools
  thresh_cold_pool, & ! threshold for targeted diffusion for cold pools
  rdheight,         & ! bottom height of Rayleigh damping layer
  nstop,            & ! last time step of the forecast period
  nnextbound,       & ! time step of the next boundary update
  nincbound,        & ! time step increment of boundary update
  nlastbound,       & ! time step of the last boundary update
  rlwidth,          & ! width of relaxation layer (if lexpl_lbc=.TRUE.)
  idbg_level,       & ! to control the verbosity of debug output
  lprintdeb_all,    & ! .TRUE.:  all tasks print debug output
                      ! .FALSE.: only task 0 prints debug output
  nbd1,             &
  nbd2,             &
  lcond
! end of data_runcontrol 


USE turb_data, ONLY: &
  imode_turb          ! mode of turbulent diffusion parametrization

#ifdef NUDGING
USE data_lheat_nudge, ONLY : &
  llhn         ,& ! on/off switch for latent heat nudging (lhn)
  llhnverif    ,& ! on/off switch for verification against radar
  tt_lheat        ! profile of t-increments due to latent heating   ( K/s )
                  ! (stored for current and previous timestep) 
#endif

! end of data_lheat_nudge

USE src_artifdata , ONLY :   &
  lnosurffluxes_m  ! switch to turn on free-slip lower BC by setting tcm to 0.0

! end of src_artifdata

#ifdef POLLEN
USE data_pollen, ONLY :  &
  trcr_meta_sedimvel
#endif

!------------------------------------------------------------------------------

USE parallel_utilities,       ONLY : i_global, j_global

USE time_utilities,           ONLY :  get_timings, i_dyn_computations, &
  i_cppdycore_copyin,                                         &
  i_cppdycore_step,                                           &
  i_cppdycore_copyout,                                        &
  i_cppdycore_swap

USE utilities,           ONLY:  to_upper

USE iso_c_binding , ONLY :  &
  C_PTR, C_NULL_PTR, C_ASSOCIATED, C_NULL_CHAR, C_LOC, &
  C_INT

USE src_tracer          ,     ONLY :  &
  trcr_get, trcr_meta_get, trcr_errorstr, trcr_get_ntrcr


USE iso_c_binding        ,    ONLY : C_FLOAT, C_DOUBLE

!==============================================================================

IMPLICIT NONE

! default private
PRIVATE

#ifdef SINGLEPRECISION
# define cwp C_FLOAT
#else
# define cwp C_DOUBLE
#endif

INTERFACE set_dycore_param
  MODULE PROCEDURE                     &
    set_dycore_param_real,             &
    set_dycore_param_int,              &
    set_dycore_param_logical,          &
    set_dycore_param_string
END INTERFACE

INTERFACE dycore_set_2dfield_ij
  MODULE PROCEDURE                     &
    dycore_set_2dfield_ij_withconnectorname, &
    dycore_set_2dfield_ij_simple
END INTERFACE

INTERFACE dycore_tracer_define_metainfo
  MODULE PROCEDURE                     &
    dycore_tracer_define_metainfo_i
END INTERFACE

INTERFACE dycore_tracer_set_metainfo
  MODULE PROCEDURE                     &
    dycore_tracer_set_metainfo_i
END INTERFACE

!==============================================================================

!------------------------------------------------------------------------------
! Declarations
!------------------------------------------------------------------------------

INTEGER                  :: &
  izdebug      !  MPI task specific debug level (usually muted on all PEs except PE 0)

INTEGER                  :: &
  icomm_gclcart      !  3D cartesian communicator for GCL
 
INTEGER                  :: nnextbound_local = -1

!==============================================================================

! Declare public entities
PUBLIC :: &
  cpp_dycore_init,              &
  cpp_dycore_compute,           &
  cpp_dycore_swap_pointers,     &
  cpp_dycore_cleanup

CONTAINS

!==============================================================================
!+ Module procedure in "src_cpp_dycore" for initializing the C++ dycore
!------------------------------------------------------------------------------

SUBROUTINE cpp_dycore_init()

!------------------------------------------------------------------------------
!
! Description: This routine initializes the C++ Dycore which connects
!        fortran COSMO with the C++ Dycore. 
!
! Method: It calls C++ Dycore functions that connect data fields in fortran
!        with data fields resident in a repository in the cpp dycore. 
!        After setup, the C++ Dycore is initialized
!------------------------------------------------------------------------------

IMPLICIT NONE

! external functions
!-------------------------- 
 INTERFACE

    ! Dycore init
    SUBROUTINE dycorewrapper_init_setup( ie, je, ke, ie_tot, je_tot, icomm_gclcart, lperi_x, lperi_y, &
         pe_at_wb, pe_at_eb, pe_at_sb, pe_at_nb, num_compute, isubpos) &
         BIND(c, name='dycorewrapper_init_setup')
      USE, INTRINSIC     :: iso_c_binding
      INTEGER(C_INT), value     :: ie, je, ke, ie_tot, je_tot, icomm_gclcart, num_compute
      LOGICAL(C_BOOL), value    :: lperi_x, lperi_y, pe_at_wb, pe_at_eb, pe_at_sb, pe_at_nb
      INTEGER(C_INT), DIMENSION(*)   :: isubpos 
    END SUBROUTINE dycorewrapper_init_setup

    ! Dycore freeze configuration
    SUBROUTINE dycorewrapper_param_freeze()  &
         BIND(c, name='dycorewrapper_param_freeze')
    END SUBROUTINE dycorewrapper_param_freeze

    ! Dycore finalize setup
    SUBROUTINE dycorewrapper_finalize_setup()  &
         BIND(c, name='dycorewrapper_finalize_setup') 
    END SUBROUTINE dycorewrapper_finalize_setup

    ! Check decomposition
    SUBROUTINE dycorewrapper_domaindecomposition_check( pos1, pos2, pos3 ) &
         BIND(c, name='dycorewrapper_domaindecomposition_check')
      USE, INTRINSIC         :: iso_c_binding
      INTEGER(C_INT), value  :: pos1, pos2, pos3
    END SUBROUTINE dycorewrapper_domaindecomposition_check
    
    ! Dycore finalize setup
    SUBROUTINE dycorewrapper_copy_initial_bdfields()  &
         BIND(c, name='dycorewrapper_copy_initial_bdfields')
    END SUBROUTINE dycorewrapper_copy_initial_bdfields

    ! Dycore setup tracers
    SUBROUTINE dycorewrapper_setup_tracers()  &
         BIND(c, name='dycorewrapper_setup_tracers')
    END SUBROUTINE dycorewrapper_setup_tracers
    
 END INTERFACE

! Local scalars:
! -------------

INTEGER ::  &
  i, j, k, izerror

CHARACTER (LEN=255) :: yzerrmsg
CHARACTER (LEN=25) :: yzroutine

LOGICAL(KIND=1) ::  l_lperi_x, l_lperi_y

INTEGER :: idumm, irank

INTEGER :: cartdims(3)  ! # procs for each dimension of cart comm

LOGICAL :: lreorder             ! defines if rank can be reordered in cartesian communicator

LOGICAL :: cartperiods(3)       ! defines if cart comm is periodic

INTEGER :: nznumdims  !  number of dimensions of cart comm

INTEGER :: izmplcode            ! error code

LOGICAL(KIND=1) :: pe_at_wb, pe_at_eb, pe_at_nb, pe_at_sb ! determines whether the PE is
                                                  ! at the W, E, N or S global boundary

INTEGER :: iScalarAdvectionScheme

LOGICAL :: lStrangSplitting

REAL (KIND=wp), TARGET :: dummy_2d(1,1), dummy_3d(1,1,1)

! Tracers metainformation IDs
INTEGER :: adv_id
INTEGER :: diff_id
INTEGER :: turb_id
INTEGER :: conv_id
INTEGER :: lbc_id
INTEGER :: bbc_id
INTEGER :: relax_id
INTEGER :: damp_id
INTEGER :: clp_id
INTEGER :: sedim_id

CHARACTER (LEN=255)      :: tracername
INTEGER :: tracerid
INTEGER :: iztrcr

REAL (KIND=wp), POINTER :: &
  trcr_new(:,:,:)      => NULL(),  & ! tracer data field pointer (tlev=nnew)
  trcr_now(:,:,:)      => NULL(),  & ! tracer data field pointer (tlev=nnow)
  trcr_s(:,:,:)        => NULL(),  & ! tracer surface field pointer
  trcr_sedimvel(:,:,:) => NULL(),  & ! tracer sedimentation velocity pointer
  trcr_s_new(:,:)      => NULL(),  & ! tracer surface field pointer (tlev=nnew)
  trcr_s_now(:,:)      => NULL(),  & ! tracer surface field pointer (tlev=nnow)
  trcr_tens(:,:,:)     => NULL(),  & ! tracer tendency field pointer
  trcr_bd(:,:,:,:)     => NULL(),  & ! tracer boundary field pointer
  trcr_bd1(:,:,:)      => NULL(),  & ! tracer boundary field pointer at nbd1
  trcr_bd2(:,:,:)      => NULL()     ! tracer boundary field pointer at nbd2

LOGICAL(KIND=1) :: lValidBD, lValidSurf, lValidSedimvel

! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin module procedure cpp_dycore_init
!------------------------------------------------------------------------------

yzroutine = 'cpp_dycore_init'
yzerrmsg  = ''
izerror   = 0

! Initialize, whether debug output shall be done
IF (lprintdeb_all) THEN
  izdebug = idbg_level
ELSE
  IF (my_cart_id == 0) THEN
    izdebug = idbg_level
  ELSE
    izdebug = 0
  ENDIF
ENDIF

! save value of nnextbound in a local copy
nnextbound_local = nnextbound

! declare GPU pointers
!$acc data create(dummy_2d,dummy_3d)

!------------------------------------------------------------------------------
! Section 1: initialize the C++ Dycore
!------------------------------------------------------------------------------ 

  ! check that Fortran and dycore code are compatible
  CALL dycore_check()

  ! creating a new 3D cartesian communicator (degenerated copy of icomm_compute)
  cartdims(1) = nprocx
  cartdims(2) = nprocy
  cartdims(3) = 1
  nznumdims = 3

  cartperiods(1) = lperi_x
  cartperiods(2) = lperi_y
  cartperiods(3) = .TRUE.
  lreorder = .FALSE.

  CALL MPI_CART_CREATE( icomm_compute, nznumdims, cartdims, cartperiods , lreorder, &
      icomm_gclcart, izmplcode);

  CALL MPI_COMM_RANK( icomm_gclcart, irank, izmplcode );

  pe_at_wb = .false.
  pe_at_eb = .false.
  pe_at_nb = .false.
  pe_at_sb = .false.

  ! west
  IF ( my_cart_neigh(1) == -1) THEN 
    pe_at_wb = .true.
  ENDIF
  ! east
  IF ( my_cart_neigh(3) == -1) THEN
    pe_at_eb = .true.
  ENDIF
  ! south
  IF ( my_cart_neigh(4) == -1) THEN
    pe_at_sb = .true.
  ENDIF
  ! north
  IF ( my_cart_neigh(2) == -1) THEN
    pe_at_nb = .true.
  ENDIF

  ! echo PE and domain setup to stdout
  !DO k = 0, num_compute-1
  !  IF (irank == k) THEN
  !    WRITE(*,'(a,i4)')  '*** GRID SETUP OF RANK ', my_cart_id
  !    WRITE(*,'(a15,3(i8))') 'pos', my_cart_pos(1), my_cart_pos(2), my_cart_pos(3)
  !    WRITE(*,'(a15,4(a8))') 'dir', 'west', 'east', 'south', 'north'
  !    WRITE(*,'(a15,4(l8))') 'border', pe_at_wb, pe_at_eb, pe_at_sb, pe_at_nb
  !    WRITE(*,'(a15,4(i8))') 'ij-range', i_global(1), i_global(ie), j_global(1), j_global(je)
  !    WRITE(*,'(a15,4(i8))') 'neigh', my_cart_neigh(1), my_cart_neigh(3), my_cart_neigh(4), my_cart_neigh(2)
  !  END IF
  !  CALL MPI_BARRIER(icomm_compute, izmplcode)
  !END DO

  ! format convertion from LOGICAL to LOGICAL*1
  l_lperi_x = lperi_x
  l_lperi_y = lperi_y

  ! initialize the C++ Dycore
  CALL dycorewrapper_init_setup( ie, je, ke, ie_tot, je_tot, icomm_gclcart, l_lperi_x, l_lperi_y, &
          pe_at_wb, pe_at_eb, pe_at_sb, pe_at_nb, num_compute, isubpos)

!------------------------------------------------------------------------------
! Section 2: pass configuration parameters to the Dycore
!  
!------------------------------------------------------------------------------ 

  CALL set_dycore_param("dt", dt)
  CALL set_dycore_param("dlon", dlon)
  CALL set_dycore_param("dlat", dlat)
  CALL set_dycore_param("xkd", xkd)
  CALL set_dycore_param("hd_corr_u_bd", hd_corr_u_bd)
  CALL set_dycore_param("hd_corr_u_in", hd_corr_u_in)
  CALL set_dycore_param("hd_corr_trcr_bd", hd_corr_trcr_bd)
  CALL set_dycore_param("hd_corr_trcr_in", hd_corr_trcr_in)
  CALL set_dycore_param("hd_corr_t_bd", hd_corr_t_bd)
  CALL set_dycore_param("hd_corr_t_in", hd_corr_t_in)
  CALL set_dycore_param("hd_corr_p_bd", hd_corr_p_bd)
  CALL set_dycore_param("hd_corr_p_in", hd_corr_p_in)
  CALL set_dycore_param("hd_dhmax", hd_dhmax)
  CALL set_dycore_param("rdheight", rdheight)
  CALL set_dycore_param("rlwidth", rlwidth)
  CALL set_dycore_param("divdamp_slope", divdamp_slope)

  CALL set_dycore_param("itype_fast_waves", itype_fast_waves)
  CALL set_dycore_param("itype_hdiff", itype_hdiff)
  CALL set_dycore_param("irk_order", irk_order)
  CALL set_dycore_param("irunge_kutta", irunge_kutta)
  CALL set_dycore_param("nincbound", nincbound)
  CALL set_dycore_param("kflat", vcoord%kflat)
  CALL set_dycore_param("iadv_order", iadv_order)
  CALL set_dycore_param("imode_turb", imode_turb)
  CALL set_dycore_param("itype_gscp", itype_gscp)

  CALL set_dycore_param("lcpp_dycore", lcpp_dycore)
  CALL set_dycore_param("ltadv_limiter", ltadv_limiter)
  CALL set_dycore_param("ltur", ltur)
  CALL set_dycore_param("lnosurffluxes_m", lnosurffluxes_m)
  CALL set_dycore_param("l_diff_cold_pools", l_diff_cold_pools)
  CALL set_dycore_param("thresh_cold_pool", thresh_cold_pool)
  CALL set_dycore_param("lcond", lcond)
  CALL set_dycore_param("lcori", lcori)
  CALL set_dycore_param("lhordiff", lhordiff)
  CALL set_dycore_param("ldiabf_lh", ldiabf_lh)
  CALL set_dycore_param("ldyn_bbc", ldyn_bbc)
  CALL set_dycore_param("lsppt", lsppt)
#ifdef NUDGING
  CALL set_dycore_param("llhn", llhn)
  CALL set_dycore_param("llhnverif", llhnverif)
#else
  CALL set_dycore_param("llhn", .FALSE.)
  CALL set_dycore_param("llhnverif", .FALSE.)
#endif
  CALL set_dycore_param("y_scalar_advect", y_scalar_advect)
  CALL set_dycore_param("l_diff_smag", l_diff_smag)
  CALL set_dycore_param("itype_spubc", itype_spubc)
  CALL set_dycore_param("r_earth_cosmo", r_earth)

  ! TODO: this should be enforced to be consistent with l_calc_lhs_at_1st_RKstep in fast_waves_sc.f90
  CALL set_dycore_param("lLHSAt1stRKStage", .TRUE.)

  CALL dycorewrapper_param_freeze()

!------------------------------------------------------------------------------
! Section 3: setup the C++ Dycore with all the fields needed by the
!    dycore. Missing fields will trigger assertions in the dycore. 
!    Field names passed to dycore should pass naming convention adopted
!    by the dycore (see DycoreWrapperRepository::Init);
!    otherwise an assertion is triggered. 
!------------------------------------------------------------------------------ 

  ! feed the C++ Dycore with set of fields and corresponding pointers
  CALL dycore_set_3dfield   ("hhl", hhl)
  CALL dycore_set_3dfield   ("rho0", rho0)
  CALL dycore_set_3dfield   ("p0", p0)
  CALL dycore_set_3dfield   ("p0hl", p0hl)
  CALL dycore_set_3dfield   ("t0hl", t0hl)
  CALL dycore_set_3dfield   ("dt0dz", dt0dz)
  CALL dycore_set_3dfield   ("t0", t0)
  CALL dycore_set_3dfield   ("sqrtgrs", sqrtg_r_s)
  CALL dycore_set_3dfield   ("sqrtgru", sqrtg_r_u)
  CALL dycore_set_3dfield   ("sqrtgrv", sqrtg_r_v)
  CALL dycore_set_3dfield   ("sqrtgrw", sqrtg_r_w)
  CALL dycore_set_3dfield   ("hdmask", hd_mask)
  CALL dycore_set_1dfield_k ("hhlr", vcoord%vert_coord(1:ke1)) ! AAAAA TODO This does not work if je == ke
  CALL dycore_set_1dfield_k ("rdcoef", rdcoef)
  CALL dycore_set_2dfield_ij("rmyScalar", rmy(:,:,1))
  CALL dycore_set_2dfield_ij("rmyU", rmy(:,:,2))
  CALL dycore_set_2dfield_ij("rmyV", rmy(:,:,3))
  CALL dycore_set_2dfield_ij("rmyQrQsQg", rmyq)

  CALL dycore_set_3dfield   ("u_nnew", u(:,:,:,nnew))
  CALL dycore_set_3dfield   ("u_nnow", u(:,:,:,nnow))
  CALL dycore_set_3dfield   ("v_nnew", v(:,:,:,nnew))
  CALL dycore_set_3dfield   ("v_nnow", v(:,:,:,nnow))
  CALL dycore_set_3dfield   ("w_nnew", w(:,:,:,nnew))
  CALL dycore_set_3dfield   ("w_nnow", w(:,:,:,nnow))
  CALL dycore_set_3dfield   ("t_nnew", t(:,:,:,nnew))
  CALL dycore_set_3dfield   ("t_nnow", t(:,:,:,nnow))
  CALL dycore_set_3dfield   ("pp_nnew", pp(:,:,:,nnew))
  CALL dycore_set_3dfield   ("pp_nnow", pp(:,:,:,nnow))
  CALL dycore_set_3dfield   ("rho", rho)
  CALL dycore_set_3dfield   ("u_tens", utens)
  CALL dycore_set_3dfield   ("v_tens", vtens)
  CALL dycore_set_3dfield   ("w_tens", wtens)
  CALL dycore_set_3dfield   ("tp_tens", ttens)
  CALL dycore_set_3dfield   ("pp_tens", pptens)
  CALL dycore_set_3dfield   ("qrs", qrs)
  CALL dycore_set_3dfield   ("dqvdt", dqvdt)
  CALL dycore_set_3dfield   ("tkvm", tkvm)
  CALL dycore_set_3dfield   ("tkvh", tkvh)
  IF (ldiabf_lh) THEN
    CALL dycore_set_3dfield   ("tinc_lh", tinc_lh)
  ENDIF
  IF (lsppt) THEN
    CALL dycore_set_3dfield   ("pertstoph", pertstoph)
  ENDIF
#ifdef NUDGING
  IF (llhn .OR. llhnverif) THEN
    CALL dycore_set_3dfield   ("tt_lheat_nnew", tt_lheat(:,:,:,nnew))
    CALL dycore_set_3dfield   ("tt_lheat_nnow", tt_lheat(:,:,:,nnow))
  END IF
#endif
  CALL dycore_set_2dfield_ij("fc", fc)
  CALL dycore_set_2dfield_ij("tcm", tcm)
  CALL dycore_set_2dfield_ij("tch", tch)
  CALL dycore_set_2dfield_ij("p_s_nnow", "p_s_nnow", ps(:,:,nnow))
  CALL dycore_set_2dfield_ij("t_g_nnow", "t_g_nnow", t_g(:,:,nnow))
  CALL dycore_set_2dfield_ij("p_s_nnew", "p_s_nnew", ps(:,:,nnew))
  CALL dycore_set_2dfield_ij("t_g_nnew", "t_g_nnew", t_g(:,:,nnew))
  CALL dycore_set_2dfield_ij("qvsflx", "qvsflx", qvsflx(:,:))
  CALL dycore_set_2dfield_ij("lhfl_s", "lhfl_s", lhfl_s(:,:))
  CALL dycore_set_2dfield_ij("shfl_s", "shfl_s", shfl_s(:,:))
  CALL dycore_set_2dfield_ij("umfl_s", "umfl_s", umfl_s(:,:))
  CALL dycore_set_2dfield_ij("vmfl_s", "vmfl_s", vmfl_s(:,:))
  CALL dycore_set_1dfield_k("a1t", a1t)
  CALL dycore_set_1dfield_k("a2t", a2t)
  CALL dycore_set_1dfield_j("crlat0", crlat(:,1))
  CALL dycore_set_1dfield_j("crlat1", crlat(:,2))
  CALL dycore_set_1dfield_j("acrlat0", acrlat(:,1))
  CALL dycore_set_1dfield_j("acrlat1", acrlat(:,2))
  CALL dycore_set_1dfield_j("tgrlat0", tgrlat(:,1))
  CALL dycore_set_1dfield_j("tgrlat1", tgrlat(:,2))
  CALL dycore_set_3dfield("u_nnew_bd", u(:,:,:,nnew))
  CALL dycore_set_3dfield("v_nnew_bd", v(:,:,:,nnew))
  CALL dycore_set_3dfield("t_nnew_bd", t(:,:,:,nnew))
  CALL dycore_set_3dfield("pp_nnew_bd", pp(:,:,:,nnew))
  CALL dycore_set_3dfield("u_nnow_bd", u(:,:,:,nnow))
  CALL dycore_set_3dfield("v_nnow_bd", v(:,:,:,nnow))
  CALL dycore_set_3dfield("t_nnow_bd", t(:,:,:,nnow))
  CALL dycore_set_3dfield("pp_nnow_bd", pp(:,:,:,nnow))

  CALL dycore_set_3dfield("u_bd1", u_bd(:,:,:,nbd1))
  CALL dycore_set_3dfield("u_bd2", u_bd(:,:,:,nbd2))
  CALL dycore_set_3dfield("v_bd1", v_bd(:,:,:,nbd1))
  CALL dycore_set_3dfield("v_bd2", v_bd(:,:,:,nbd2))
  CALL dycore_set_3dfield("t_bd1", t_bd(:,:,:,nbd1))
  CALL dycore_set_3dfield("t_bd2", t_bd(:,:,:,nbd2))
  CALL dycore_set_3dfield("pp_bd1", pp_bd(:,:,:,nbd1))
  CALL dycore_set_3dfield("pp_bd2", pp_bd(:,:,:,nbd2))

!------------------------------------------------------------------------------
! Section 3a: setup the tracers. The dycore supports now tracers with
!             metainformation and there is not difference between moisture
!             and regular tracers anymore.
!------------------------------------------------------------------------------ 

  ! Register the required metainformation
  CALL dycore_tracer_define_metainfo("ADV",   0, adv_id)
  CALL dycore_tracer_define_metainfo("DIFF",  0, diff_id)
  CALL dycore_tracer_define_metainfo("TURB",  0, turb_id)
  CALL dycore_tracer_define_metainfo("CONV",  0, conv_id)
  CALL dycore_tracer_define_metainfo("LBC",   0, lbc_id)
  CALL dycore_tracer_define_metainfo("BBC",   0, bbc_id)
  CALL dycore_tracer_define_metainfo("RELAX", 0, relax_id)
  CALL dycore_tracer_define_metainfo("DAMP",  0, damp_id)
  CALL dycore_tracer_define_metainfo("CLP",   0, clp_id)
!  CALL dycore_tracer_define_metainfo("SEDIM", 0, sedim_id)
!  CALL dycore_tracer_define_metainfo("SPPT_PERT", 0, spptpert_id)

  ! Register tracers and set metainfo
  DO iztrcr = 1, trcr_get_ntrcr()

    ! Register tracer
    CALL trcr_meta_get( izerror, iztrcr, T_NAME_ID, tracername )
    IF ( izerror /= 0 ) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort (my_world_id, 9031 , yzerrmsg, yzroutine )
    ENDIF
    CALL dycore_tracer_add(tracername, tracerid)

    ! Check that C++ ID is Fortran ID - 1
    IF (tracerid /= iztrcr-1) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort (my_world_id, 9035 , yzerrmsg, yzroutine )
    ENDIF

    ! Set metainfo
    CALL dycore_tracer_set_metainfo(tracerid, adv_id, iztrcr, T_ADV_ID)
    CALL dycore_tracer_set_metainfo(tracerid, diff_id, iztrcr, T_DIFF_ID)
    CALL dycore_tracer_set_metainfo(tracerid, turb_id, iztrcr, T_TURB_ID)
    CALL dycore_tracer_set_metainfo(tracerid, conv_id, iztrcr, T_CONV_ID)
    CALL dycore_tracer_set_metainfo(tracerid, lbc_id, iztrcr, T_LBC_ID)
    CALL dycore_tracer_set_metainfo(tracerid, bbc_id, iztrcr, T_BBC_ID)
    CALL dycore_tracer_set_metainfo(tracerid, relax_id, iztrcr, T_RELAX_ID)
    CALL dycore_tracer_set_metainfo(tracerid, damp_id, iztrcr, T_DAMP_ID)
    CALL dycore_tracer_set_metainfo(tracerid, clp_id, iztrcr, T_CLP_ID)
!    CALL dycore_tracer_set_metainfo(tracerid, spptpert_id, iztrcr, T_SPPTPERT_ID)
!    CALL dycore_tracer_set_metainfo(tracerid, sedim_id, iztrcr, T_SEDIM_ID)

  END DO

  ! freeze metainformation and allocate fields
  CALL dycorewrapper_setup_tracers()

  ! connect tracer data fields
  DO iztrcr = 1, trcr_get_ntrcr()

    ! Get name
    CALL trcr_meta_get( izerror, iztrcr, T_NAME_ID, tracername )
    IF ( izerror /= 0 ) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort (my_world_id, 9031 , yzerrmsg, yzroutine )
    ENDIF

    ! Get nnew, tens and nnow fields
    CALL trcr_get( izerror, iztrcr, ptr_tlev=nnew, ptr=trcr_new, ptr_tens=trcr_tens)
    CALL trcr_get( izerror, iztrcr, ptr_tlev=nnow, ptr=trcr_now)
    IF ( (.NOT. ASSOCIATED(trcr_new)) .OR. (.NOT. ASSOCIATED(trcr_now)) ) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort (my_world_id, 9037 , yzerrmsg, yzroutine )
    ENDIF

    ! Try to get boundary field
    CALL trcr_get( izerror, iztrcr, ptr_bd_notlev=trcr_bd1,ptr_tlev=nbd1)
    CALL trcr_get( izerror, iztrcr, ptr_bd_notlev=trcr_bd2,ptr_tlev=nbd2)
    IF (izerror == 0) THEN
      lValidBD = .TRUE.      
    ELSE
      lValidBD = .FALSE.
      trcr_bd1 => dummy_3d
      trcr_bd2 => dummy_3d
    END IF

    ! Try to get surface field
    CALL trcr_meta_get(izerror, iztrcr, "SURF_FIELD", trcr_s)
    IF ( ASSOCIATED(trcr_s) ) THEN
      lValidSurf = .TRUE.
      trcr_s_new => trcr_s(:,:,nnew)
      trcr_s_now => trcr_s(:,:,nnow)
    ELSE
      lValidSurf = .FALSE.
      trcr_s_new => dummy_2d
      trcr_s_now => dummy_2d
    END IF

    lValidSedimvel = .FALSE.
#ifdef POLLEN
    ! Try to get sedim velocity field
    CALL trcr_meta_get(izerror, iztrcr, trcr_meta_sedimvel, trcr_sedimvel)
    IF ( ASSOCIATED(trcr_sedimvel) ) THEN
      lValidSedimvel = .TRUE.
    ELSE
      ! FUO HACK: this is currently required since it is not legal to pass a NULL pointer to
      !           the dycore_tracer_connect subroutine.
      trcr_sedimvel => trcr_new
    END IF
#endif

    CALL trcr_meta_get( izerror, iztrcr, T_NAME_ID, tracername )

    ! Call connect routine
    CALL dycore_tracer_connect(tracername, iztrcr-1, trcr_new, trcr_now, trcr_tens, &
                               trcr_s_new, trcr_s_now, trcr_bd1, trcr_bd2, trcr_sedimvel, &
                               lValidSurf, lValidBD, lValidSedimvel)

  END DO

!------------------------------------------------------------------------------
! Section 4: finalize the C++ Dycore setup
!    Here the C++ Dycore is initialized.
!------------------------------------------------------------------------------ 

!!$  CALL check_dycore_setup(800)

  CALL dycorewrapper_finalize_setup()

#ifdef CPP_DEBUG
  IF (num_compute > 1) THEN 
    idumm = 0
    ! If symbol to next function is not found, compile Dycore in DEBUG mode
    CALL dycorewrapper_domaindecomposition_check( my_cart_pos(1), my_cart_pos(2), idumm )
  ENDIF
#endif

  ! at first time step we copy nbd1 & nbd2 into C++ dycore
  CALL dycorewrapper_copy_initial_bdfields()
 
!$acc end data

  IF (ltime) CALL get_timings (i_dyn_computations, ntstep, dt, izerror)

!----------------------------------------------------------------------------
! End of module procedure cpp_dycore_init
!----------------------------------------------------------------------------

END SUBROUTINE cpp_dycore_init


SUBROUTINE dycore_set_3dfield( fieldname, field )

  ! arguments
  CHARACTER (LEN=*), INTENT(IN) :: fieldname  ! name of field
  REAL (KIND=wp), INTENT(IN), TARGET :: field(:,:,:)

  ! external functions
  INTERFACE
     ! Pass pointer and dimensions of ijk fields to dycore
     SUBROUTINE dycorewrapper_set_ijkfield(field, isize, jsize, ksize, fieldname_size, fieldname_wrapper) &
          BIND(c, name='dycorewrapper_set_ijkfield')
       USE, INTRINSIC :: iso_c_binding
       INTEGER(C_INT), value    :: isize,jsize,ksize
       TYPE(C_PTR), value       :: field
       INTEGER(C_INT), value    :: fieldname_size
       CHARACTER(KIND=C_CHAR), DIMENSION(*)  :: fieldname_wrapper
     END SUBROUTINE dycorewrapper_set_ijkfield
  END INTERFACE
  
  ! locals
  CHARACTER (LEN=16) :: fieldname_wrapper        !XL_TODO changed from 12 to 16, this should be done for other set method, and a check should be added if LEN_TRIM(fieldname) > LEN(fieldname_wrapper)
  INTEGER                :: fieldname_size
  REAL(KIND=wp), POINTER :: padd  !XL_HACK : required for gnu < 4.9

  fieldname_wrapper = '                '
  fieldname_wrapper = fieldname
  fieldname_size = LEN_TRIM(fieldname_wrapper)

  !$acc data present(field)
  !$acc host_data use_device(field)

  IF (izdebug > 20) THEN
    WRITE(*,*) 'INFO-DBG: dycore_set_3dfield for '//TRIM(fieldname)
  ENDIF

  padd=>field(1,1,1) !XL_HACK : required for gnu < 4.9
  CALL dycorewrapper_set_ijkfield( C_LOC(padd), SIZE(field,1), SIZE(field,2), SIZE(field,3), fieldname_size, fieldname_wrapper ) !XL_HACK
  !XL_HACK CALL dycorewrapper_set_ijkfield( C_LOC(field(1,1,1)), SIZE(field,1), SIZE(field,2), SIZE(field,3), fieldname_size, fieldname_wrapper )
  !$acc end host_data
  !$acc end data

END SUBROUTINE dycore_set_3dfield


SUBROUTINE dycore_set_2dfield_ij_withconnectorname( connectorName, fieldname, field )
  
  ! arguments
  CHARACTER (LEN=*), INTENT(IN) :: connectorName ! name of connector
  CHARACTER (LEN=*), INTENT(IN) :: fieldname  ! name of field
  REAL (KIND=wp), INTENT(IN), TARGET :: field(:,:)

  ! external functions
  INTERFACE
     ! Pass pointer and dimensions of ij fields to dycore
     SUBROUTINE dycorewrapper_set_ijfield(field, isize, jsize,  fieldname_size, connectorname_size, fieldname_wrapper, connectorname_wrapper ) &
          BIND(c, name='dycorewrapper_set_ijfield')
       USE, INTRINSIC :: iso_c_binding
       INTEGER(C_INT), value      :: isize, jsize
       TYPE(C_PTR), value         :: field
       INTEGER(C_INT), value      :: fieldname_size, connectorname_size
       CHARACTER(KIND=C_CHAR), DIMENSION(*)  :: fieldname_wrapper, connectorname_wrapper
     END SUBROUTINE dycorewrapper_set_ijfield
  END INTERFACE

  ! locals
  CHARACTER (LEN=12) :: fieldname_wrapper
  INTEGER            :: fieldname_size

  CHARACTER (LEN=12) :: connectorname_wrapper
  INTEGER            :: connectorname_size
  REAL(KIND=wp), POINTER :: padd  !XL_HACK : required for gnu < 4.9

  fieldname_wrapper = '            '
  fieldname_wrapper = fieldname
  fieldname_size = LEN_TRIM(fieldname_wrapper)

  connectorname_wrapper = '            '
  connectorname_wrapper = connectorName
  connectorname_size = LEN_TRIM(connectorname_wrapper)

  !$acc data present(field)
  !$acc host_data use_device(field)

  IF (izdebug > 20) THEN
    WRITE(*,*) 'INFO-DBG: dycore_set_2dfield_ij_withconnectorname for '//  &
               TRIM(fieldname)//' and connector '//TRIM(connectorName)
  ENDIF

  padd=>field(1,1) !XL_HACK : required for gnu < 4.9
  CALL dycorewrapper_set_ijfield( C_LOC(padd), SIZE(field,1), SIZE(field,2), fieldname_size, connectorname_size, fieldname_wrapper, connectorname_wrapper ) !XL_HACK
  !XL_HACK CALL dycorewrapper_set_ijfield( C_LOC(field(1,1)), SIZE(field,1), SIZE(field,2), fieldname_size, connectorname_size, fieldname_wrapper, connectorname_wrapper )
  !$acc end host_data
  !$acc end data

END SUBROUTINE dycore_set_2dfield_ij_withconnectorname


SUBROUTINE dycore_set_2dfield_ij_simple( fieldname, field )

  CHARACTER (LEN=*), INTENT(IN) :: fieldname  ! name of field
  REAL (KIND=wp), INTENT(IN) :: field(:,:)

  IF (izdebug > 20) THEN
    WRITE(*,*) 'INFO-DBG: dycore_set_2dfield_ij_simple for '//TRIM(fieldname)
  ENDIF

  CALL dycore_set_2dfield_ij_withconnectorname( fieldname, fieldname, field )

END SUBROUTINE dycore_set_2dfield_ij_simple


SUBROUTINE dycore_set_1dfield_j( fieldname, field )
  
  ! arguments
  CHARACTER (LEN=*), INTENT(IN) :: fieldname  ! name of field
  REAL (KIND=wp), INTENT(IN), TARGET :: field(:)

    ! external functions
  INTERFACE

     ! Pass pointer and dimension j fields to dycore
     SUBROUTINE dycorewrapper_set_jfield(field, jsize, fieldname_size, fieldname_wrapper) &
          BIND(c, name='dycorewrapper_set_jfield')
       USE, INTRINSIC :: iso_c_binding
       INTEGER(C_INT), value      :: jsize
       TYPE(C_PTR), value         :: field
       INTEGER(C_INT), value      :: fieldname_size
       CHARACTER(KIND=C_CHAR), DIMENSION(*)  :: fieldname_wrapper
     END SUBROUTINE dycorewrapper_set_jfield

  END INTERFACE
  
  ! locals
  CHARACTER (LEN=12) :: fieldname_wrapper
  INTEGER                :: fieldname_size
  REAL(KIND=wp), POINTER :: padd  !XL_HACK : required for gnu < 4.9

  fieldname_wrapper = '            '
  fieldname_wrapper = fieldname
  fieldname_size = LEN_TRIM(fieldname_wrapper)

  !$acc data present(field)
  !$acc host_data use_device(field)

  IF (izdebug > 20) THEN
    WRITE(*,*) 'INFO-DBG: dycore_set_1dfield_j for '//TRIM(fieldname)
  ENDIF

  padd=>field(1) !XL_HACK : required for gnu < 4.9
  CALL dycorewrapper_set_jfield( C_LOC(padd), SIZE(field,1), fieldname_size, fieldname_wrapper ) !XL_HACK
  !XL_HACK CALL dycorewrapper_set_jfield( C_LOC(field(1)), SIZE(field,1), fieldname_size, fieldname_wrapper )

  !$acc end host_data
  !$acc end data

END SUBROUTINE dycore_set_1dfield_j


SUBROUTINE dycore_set_1dfield_k( fieldname, field )
  
  ! arguments
  CHARACTER (LEN=*), INTENT(IN) :: fieldname  ! name of field
  REAL (KIND=wp), INTENT(IN), TARGET :: field(:)

    ! external functions
  INTERFACE

     ! Pass pointer and dimension k fields to dycore
     SUBROUTINE dycorewrapper_set_kfield(field, ksize, fieldname_size, fieldname_wrapper) &
          BIND(c, name='dycorewrapper_set_kfield')
       USE, INTRINSIC :: iso_c_binding
       INTEGER(C_INT), value      :: ksize
       TYPE(C_PTR), value         :: field
       INTEGER(C_INT), value      :: fieldname_size
       CHARACTER(KIND=C_CHAR), DIMENSION(*)  :: fieldname_wrapper
     END SUBROUTINE dycorewrapper_set_kfield

  END INTERFACE
  
  ! locals
  CHARACTER (LEN=12) :: fieldname_wrapper
  INTEGER                :: fieldname_size
  REAL(KIND=wp), POINTER :: padd  !XL_HACK : required for gnu < 4.9

  fieldname_wrapper = '            '
  fieldname_wrapper = fieldname
  fieldname_size = LEN_TRIM(fieldname_wrapper)

  !$acc data present(field)
  !$acc host_data use_device(field)

  IF (izdebug > 20) THEN
    WRITE(*,*) 'INFO-DBG: dycore_set_1dfield_k for '//TRIM(fieldname)
  ENDIF


  padd=>field(1) !XL_HACK : required for gnu < 4.9 
  CALL dycorewrapper_set_kfield( C_LOC(padd), SIZE(field,1), fieldname_size, fieldname_wrapper ) !XL_HACK
  !XL_HACK CALL dycorewrapper_set_kfield( C_LOC(field(1)), SIZE(field,1), fieldname_size, fieldname_wrapper )

  !$acc end host_data
  !$acc end data

END SUBROUTINE dycore_set_1dfield_k


SUBROUTINE dycore_tracer_add( tracername, tracerid )

  ! arguments
  CHARACTER (LEN=*)       , INTENT(IN)  :: tracername  ! name of field
  INTEGER                 , INTENT(OUT) :: tracerid

  ! external functions
  INTERFACE
     ! register tracer variable
     SUBROUTINE dycorewrapper_add_tracervariable(tracername, namelength, tracerid) &
          BIND(c, name='dycorewrapper_add_tracervariable')
       USE, INTRINSIC :: iso_c_binding
       CHARACTER(KIND=C_CHAR), DIMENSION(*)  :: tracername
       INTEGER(C_INT), value                 :: namelength
       INTEGER(C_INT)                        :: tracerid
     END SUBROUTINE dycorewrapper_add_tracervariable
  END INTERFACE
  
  ! locals
  INTEGER                  :: namelength
  INTEGER (KIND=C_INT)     :: tracerid_

  IF (izdebug > 20) THEN
    WRITE(*,*) 'INFO-DBG: dycore_tracer_add for '//TRIM(tracername)
  ENDIF

  namelength = LEN_TRIM(tracername)
  CALL dycorewrapper_add_tracervariable( tracername, namelength, tracerid_ )
  tracerid = tracerid_

END SUBROUTINE dycore_tracer_add

SUBROUTINE dycore_tracer_define_metainfo_i( metainfoname, defaultvalue, metainfoid )

  ! arguments
  CHARACTER (LEN=*)       , INTENT(IN)  :: metainfoname
  INTEGER                 , INTENT(IN)  :: defaultvalue
  INTEGER                 , INTENT(OUT) :: metainfoid

  ! external functions
  INTERFACE
     ! define metainfo with integer value
     SUBROUTINE dycorewrapper_define_metainfo_i(metainfoname, namelength, defaultvalue, metainfoid) &
          BIND(c, name='dycorewrapper_define_metainfo_i')
       USE, INTRINSIC :: iso_c_binding
       CHARACTER(KIND=C_CHAR), DIMENSION(*)  :: metainfoname
       INTEGER(C_INT), value                 :: namelength
       INTEGER(C_INT), value                 :: defaultvalue
       INTEGER(C_INT)                        :: metainfoid
     END SUBROUTINE dycorewrapper_define_metainfo_i
  END INTERFACE
  
  ! locals
  INTEGER                  :: namelength
  INTEGER (KIND=C_INT)     :: metainfoid_

  IF (izdebug > 20) THEN
    WRITE(*,*) 'INFO-DBG: dycore_tracer_define_metainfo_i for '//TRIM(metainfoname), defaultvalue
  ENDIF

  namelength = LEN_TRIM(metainfoname)
  CALL dycorewrapper_define_metainfo_i( metainfoname, namelength, defaultvalue, metainfoid_ )
  metainfoid = metainfoid_

END SUBROUTINE dycore_tracer_define_metainfo_i


SUBROUTINE dycore_tracer_set_metainfo_i( traceriddycore, metainfoiddycore, &
                                          tracerid, metainfoid )

  ! arguments
  INTEGER , INTENT(IN) :: traceriddycore
  INTEGER , INTENT(IN) :: metainfoiddycore
  INTEGER , INTENT(IN) :: tracerid
  INTEGER , INTENT(IN) :: metainfoid

  ! external functions
  INTERFACE
    ! set metainfo with integer value
    SUBROUTINE dycorewrapper_set_metainfo_i( tracerid, metainfoid, value ) &
         BIND(c, name='dycorewrapper_set_metainfo_i')
      USE, INTRINSIC :: iso_c_binding
      INTEGER(C_INT), value :: tracerid
      INTEGER(C_INT), value :: metainfoid
      INTEGER(C_INT), value :: value
    END SUBROUTINE dycorewrapper_set_metainfo_i
  END INTERFACE

  ! local variables
  INTEGER             :: value
  INTEGER             :: izerror
  CHARACTER (LEN=255) :: yzerrmsg

  IF (izdebug > 20) THEN
    WRITE(*,*) 'INFO-DBG: dycore_tracer_set_metainfo_i for ', tracerid, metainfoid
  ENDIF

  CALL trcr_meta_get(izerror, tracerid, metainfoid, value)
  IF ( izerror /= 0 ) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort (my_world_id, 9091 , yzerrmsg, 'dycore_tracer_set_metainfo')
  ENDIF

  CALL dycorewrapper_set_metainfo_i( traceriddycore, metainfoiddycore, value )

END SUBROUTINE dycore_tracer_set_metainfo_i


SUBROUTINE dycore_tracer_connect(tracername, traceriddycore, nnew_ptr, nnow_ptr, tens_ptr,   &
                                 surf_new_ptr, surf_now_ptr, bd1_ptr, bd2_ptr, sedimvel_ptr, &
                                 valid_surf, valid_bd, valid_sedimvel )

  IMPLICIT NONE

  ! arguments
  CHARACTER (LEN=*)       , INTENT(IN)  :: tracername  ! name of field
  INTEGER                 , INTENT(IN) :: traceriddycore
  REAL (KIND=wp), POINTER ::    &
    nnew_ptr(:,:,:),   & ! tracer data field pointer (tlev=nnew)
    nnow_ptr(:,:,:),   & ! tracer data field pointer (tlev=nnow)
    tens_ptr(:,:,:),   & ! tracer tendency field pointer
    surf_new_ptr(:,:), & ! tracer surface field pointer at nnew
    surf_now_ptr(:,:), & ! tracer surface field pointer at nnow
    bd1_ptr(:,:,:),    & ! tracer boundary field pointer at nbd1
    bd2_ptr(:,:,:),    & ! tracer boundary field pointer at nbd2
    sedimvel_ptr(:,:,:)  ! tracer sedimentation velocity field pointer
  LOGICAL(KIND=1), INTENT(IN) :: valid_surf, valid_bd, valid_sedimvel

  CHARACTER(LEN=32) :: hexstr
  REAL(KIND=wp), POINTER :: padd_nnewField, padd_nnowField, padd_tensField, &      !XL_HACK : required for gnu < 4.9
       padd_surfNnewField, padd_surfNnowField, padd_bd1Field, padd_bd2Field,&      !XL_HACK : required for gnu < 4.9
       padd_sedimvelField                                                          !XL_HACK : required for gnu < 4.9

  ! external functions
  INTERFACE
    ! internal C interface
    SUBROUTINE dycorewrapper_connect_tracer( tracerid, isize, jsize, ksize, &
                                     nnewField, nnowField, tensField, &
                                     surfNnewField, surfNnowField, bd1Field, bd2Field, svField, &
                                     useSurf, useBd, useSedimvel ) &
         BIND(c, name='dycorewrapper_connect_tracer')
      USE, INTRINSIC :: iso_c_binding
      INTEGER(C_INT), value :: tracerid
      INTEGER(C_INT), value :: isize, jsize, ksize
      TYPE(C_PTR), value    :: nnewField, &
                               nnowField, &
                               tensField, &
                               surfNnewField, &
                               surfNnowField, &
                               bd1Field, &
                               bd2Field, &
                               svField
      LOGICAL(KIND=C_BOOL), value :: useSurf, useBd, useSedimvel

    END SUBROUTINE dycorewrapper_connect_tracer
  END INTERFACE

  !$acc data &
  !$acc present(nnew_ptr, nnow_ptr, tens_ptr, bd1_ptr, bd2_ptr) &
  !$acc present(surf_new_ptr, surf_now_ptr, sedimvel_ptr)

  !$acc host_data &
  !$acc use_device(nnew_ptr, nnow_ptr, tens_ptr, bd1_ptr, bd2_ptr) &
  !$acc use_device(surf_new_ptr, surf_now_ptr, sedimvel_ptr)

  IF (izdebug > 20) THEN
    WRITE(*,*) 'INFO-DBG: dycore_tracer_connect for '//TRIM(tracername), traceriddycore
    WRITE(hexstr,'(Z20)') loc(nnew_ptr(1,1,1))
    WRITE(*,*) '          nnew address: 0x' // TRIM(ADJUSTL(hexstr))
    WRITE(hexstr,'(Z20)') loc(nnow_ptr(1,1,1))
    WRITE(*,*) '          nnow address: 0x' // TRIM(ADJUSTL(hexstr))
    WRITE(hexstr,'(Z20)') loc(tens_ptr(1,1,1))
    WRITE(*,*) '          tens address: 0x' // TRIM(ADJUSTL(hexstr))
    IF (valid_surf) THEN
      WRITE(hexstr,'(Z20)') loc(surf_new_ptr(1,1))
      WRITE(*,*) '     surf nnew address: 0x' // TRIM(ADJUSTL(hexstr))
      WRITE(hexstr,'(Z20)') loc(surf_now_ptr(1,1))
      WRITE(*,*) '     surf nnow address: 0x' // TRIM(ADJUSTL(hexstr))
    ENDIF
    IF (valid_bd) THEN
      WRITE(hexstr,'(Z20)') loc(bd1_ptr(1,1,1))
      WRITE(*,*) '           bd1 address: 0x' // TRIM(ADJUSTL(hexstr))
      WRITE(hexstr,'(Z20)') loc(bd2_ptr(1,1,1))
      WRITE(*,*) '           bd2 address: 0x' // TRIM(ADJUSTL(hexstr))
    ENDIF
    IF (valid_sedimvel) THEN
      WRITE(hexstr,'(Z20)') loc(sedimvel_ptr(1,1,1))
      WRITE(*,*) '           sedimvel address: 0x' // TRIM(ADJUSTL(hexstr))
    ENDIF
  ENDIF

!XL_HACK>: Required for gnu < 4.9.0
  padd_nnewField => nnew_ptr(1,1,1)
  padd_nnowField => nnow_ptr(1,1,1)
  padd_tensField => tens_ptr(1,1,1)
  padd_surfNnewField => surf_new_ptr(1,1)
  padd_surfNnowField => surf_now_ptr(1,1)
  padd_bd1Field => bd1_ptr(1,1,1)
  padd_bd2Field => bd2_ptr(1,1,1)
  IF (valid_sedimvel) THEN
    padd_sedimvelField => sedimvel_ptr(1,1,1)
  ELSE
    NULLIFY(padd_sedimvelField)
  ENDIF
  CALL dycorewrapper_connect_tracer( traceriddycore, &
                                      size(nnew_ptr,1), &
                                      size(nnew_ptr,2), &
                                      size(nnew_ptr,3), &
                                      C_LOC(padd_nnewField)     , &
                                      C_LOC(padd_nnowField)     , &
                                      C_LOC(padd_tensField)     , &
                                      C_LOC(padd_surfNnewField) , &
                                      C_LOC(padd_surfNnowField) , &
                                      C_LOC(padd_bd1Field)      , &
                                      C_LOC(padd_bd2Field)      , &
                                      C_LOC(padd_sedimvelField) , &
                                      valid_surf, valid_bd, valid_sedimvel )
!XL_HACK>>
!XL_HACK  CALL dycorewrapper_connect_tracer( traceriddycore, &
!XL_HACK                                      size(nnew_ptr,1), &
!XL_HACK                                      size(nnew_ptr,2), &
!XL_HACK                                      size(nnew_ptr,3), &
!XL_HACK                                      C_LOC(nnew_ptr(1,1,1))     , &
!XL_HACK                                      C_LOC(nnow_ptr(1,1,1))     , &
!XL_HACK                                      C_LOC(tens_ptr(1,1,1))     , &
!XL_HACK                                      C_LOC(surf_new_ptr(1,1))   , &
!XL_HACK                                      C_LOC(surf_now_ptr(1,1))   , &
!XL_HACK                                      C_LOC(bd1_ptr(1,1,1))      , &
!XL_HACK                                      C_LOC(bd2_ptr(1,1,1))      , &
!XL_HACK                                      C_LOC(sedimvel_ptr(1,1,1)) , &
!XL_HACK                                      valid_surf, valid_bd, valid_sedimvel )
!XL_HACK<
  !$acc end host_data
  !$acc end data

END SUBROUTINE dycore_tracer_connect

SUBROUTINE dycore_tracer_verify( tracerid )

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: tracerid

  ! external functions
  INTERFACE
    ! internal C interface
    SUBROUTINE dycorewrapper_verify_tracer_pointers( tracername, namelength, &
                   tracerFieldNnew, tracerFieldNnow, tracerFieldTens, &
                   hasSurfaceField, tracerFieldSurfNnew, tracerFieldSurfNnow, &
                   hasBoundaryField, tracerFieldBD1, tracerFieldBD2 ) &
         BIND(c, name='dycorewrapper_verify_tracer_pointers')
      USE, INTRINSIC :: iso_c_binding
      CHARACTER(KIND=C_CHAR)      :: tracername
      INTEGER(KIND=C_INT), value  :: namelength
      LOGICAL(KIND=C_BOOL), value :: hasSurfaceField, hasBoundaryField
      TYPE(C_PTR), value          :: tracerFieldNnew, &
                                     tracerFieldNnow, &
                                     tracerFieldTens, &
                                     tracerFieldSurfNnew, &
                                     tracerFieldSurfNnow, &
                                     tracerFieldBD1, &
                                     tracerFieldBD2
    END SUBROUTINE dycorewrapper_verify_tracer_pointers
  END INTERFACE

  ! Internal variables
  CHARACTER(LEN=255) :: tracername
  REAL (KIND=wp), POINTER :: &
    trcr_new(:,:,:)  => NULL(),  & ! tracer data field pointer (tlev=nnew)
    trcr_now(:,:,:)  => NULL(),  & ! tracer data field pointer (tlev=nnow)
    trcr_s(:,:,:)    => NULL(),  & ! tracer surface field pointer
    trcr_s_new(:,:)  => NULL(),  & ! tracer surface field pointer (tlev=nnew)
    trcr_s_now(:,:)  => NULL(),  & ! tracer surface field pointer (tlev=nnow)
    trcr_tens(:,:,:) => NULL(),  & ! tracer tendency field pointer
    trcr_bd(:,:,:,:) => NULL(),  & ! tracer boundary field pointer
    trcr_bd1(:,:,:)  => NULL(),  & ! tracer boundary field pointer at nbd1
    trcr_bd2(:,:,:)  => NULL()     ! tracer boundary field pointer at nbd2
  REAL(KIND=wp), POINTER :: padd_nnewField, padd_nnowField, padd_tensField, &   !HACK : required for gnu < 4.9
       padd_surfNnewField, padd_surfNnowField, padd_bd1Field, padd_bd2Field     !HACK : required for gnu < 4.9
  LOGICAL(KIND=1) :: lValidSurf, lValidBD
  INTEGER :: izerror
  CHARACTER (LEN=255) :: yzerrmsg
  CHARACTER (LEN=25) :: yzroutine
  REAL (KIND=wp), TARGET :: dummy_2d(1,1), dummy_3d(1,1,1)

  yzroutine = "dycore_tracer_verify"

  ! Get tracer name
  CALL trcr_meta_get( izerror, tracerid, T_NAME_ID, tracername )
  IF ( izerror /= 0 ) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort (my_world_id, 9031 , yzerrmsg, yzroutine )
  ENDIF

  ! Get nnew, tens and nnow fields
  CALL trcr_get( izerror, tracerid, ptr_tlev=nnew, ptr=trcr_new, ptr_tens=trcr_tens)
  CALL trcr_get( izerror, tracerid, ptr_tlev=nnow, ptr=trcr_now)
  IF ( (.NOT. ASSOCIATED(trcr_new)) .OR. (.NOT. ASSOCIATED(trcr_now)) ) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort (my_world_id, 9037 , yzerrmsg, yzroutine )
  ENDIF

  ! Try to get surface field 
  CALL trcr_meta_get(izerror, tracerid, "SURF_FIELD", trcr_s)
  IF ( ASSOCIATED(trcr_s) ) THEN
    lValidSurf = .TRUE.
    trcr_s_new => trcr_s(:,:,nnew)
    trcr_s_now => trcr_s(:,:,nnow)
  ELSE
    lValidSurf = .FALSE.
    trcr_s_new => dummy_2d
    trcr_s_now => dummy_2d
  END IF

  ! Try to get boundary field
  CALL trcr_get( izerror, tracerid, ptr_bd_notlev=trcr_bd1,ptr_tlev=nbd1)
  CALL trcr_get( izerror, tracerid, ptr_bd_notlev=trcr_bd2,ptr_tlev=nbd2)
  IF (izerror == 0) THEN
    lValidBD = .TRUE.
  ELSE
    lValidBD = .FALSE.
    trcr_bd1 => dummy_3d
    trcr_bd2 => dummy_3d
  END IF

  !$acc data &
  !$acc present(trcr_new, trcr_now, trcr_tens, trcr_bd1, trcr_bd2) &
  !$acc present(trcr_s_new, trcr_s_now)

  !$acc host_data &
  !$acc use_device(trcr_new, trcr_now, trcr_tens, trcr_bd1, trcr_bd2) &
  !$acc use_device(trcr_s_new, trcr_s_now)

  !HACK: Required for gnu < 4.9.0
  padd_nnewField     => trcr_new(1,1,1)
  padd_nnowField     => trcr_now(1,1,1)
  padd_tensField     => trcr_tens(1,1,1)
  padd_surfNnewField => trcr_s_new(1,1)
  padd_surfNnowField => trcr_s_now(1,1)
  padd_bd1Field      => trcr_bd1(1,1,1)
  padd_bd2Field      => trcr_bd2(1,1,1)

  ! Call C function
  CALL dycorewrapper_verify_tracer_pointers(tracername, LEN(tracername), &
        C_LOC(padd_nnewField), C_LOC(padd_nnowField), C_LOC(padd_tensField), &
        lValidSurf, C_LOC(padd_surfNnewField), C_LOC(padd_surfNnewField), &
        lValidBD, C_LOC(padd_bd1Field), C_LOC(padd_bd2Field))

  !$acc end host_data
  !$acc end data
  
END SUBROUTINE dycore_tracer_verify

SUBROUTINE set_dycore_param_real( paramname, realvalue)

  ! arguments
  CHARACTER (LEN=*), INTENT(IN) :: paramname
  REAL (KIND=wp), INTENT(IN) :: realvalue

  ! external functions
  INTERFACE
     ! Pass real parameter to dycore
     SUBROUTINE dycorewrapper_set_realparam(realvalue, paramname_length, paramname) &
          BIND(c, name='dycorewrapper_set_realparam')
       USE, INTRINSIC :: iso_c_binding
       REAL(KIND=cwp), value :: realvalue
       INTEGER(C_INT), value      :: paramname_length
       CHARACTER(KIND=C_CHAR), DIMENSION(*)  :: paramname
     END SUBROUTINE dycorewrapper_set_realparam
  END INTERFACE
  
  ! locals
  INTEGER :: paramname_length

  IF (izdebug > 20) THEN
    WRITE(*,*) 'INFO-DBG: set_dycore_param_real for '//TRIM(paramname), realvalue
  ENDIF

  paramname_length = LEN_TRIM(paramname)
  CALL dycorewrapper_set_realparam(realvalue, paramname_length, paramname)

END SUBROUTINE set_dycore_param_real


SUBROUTINE set_dycore_param_int( paramname, intvalue)
  ! arguments
  CHARACTER (LEN=*), INTENT(IN) :: paramname
  INTEGER          , INTENT(IN) :: intvalue

  ! external functions
  INTERFACE
     ! Pass integer parameter to dycore
     SUBROUTINE dycorewrapper_set_intparam(intvalue, paramname_length, paramname) &
          BIND(c, name='dycorewrapper_set_intparam')
       USE, INTRINSIC :: iso_c_binding
       INTEGER(C_INT), value      :: intvalue, paramname_length
       CHARACTER(KIND=C_CHAR), DIMENSION(*)  :: paramname
     END SUBROUTINE dycorewrapper_set_intparam
  END INTERFACE
  
  ! locals
  INTEGER :: paramname_length

  IF (izdebug > 20) THEN
    WRITE(*,*) 'INFO-DBG: set_dycore_param_int for '//TRIM(paramname), intvalue
  ENDIF

  paramname_length = LEN_TRIM(paramname)
  CALL dycorewrapper_set_intparam(intvalue, paramname_length, paramname)

END SUBROUTINE set_dycore_param_int


SUBROUTINE set_dycore_param_logical( paramname, logicalvalue)
  ! arguments
  CHARACTER (LEN=*), INTENT(IN) :: paramname
  LOGICAL, INTENT(IN) :: logicalvalue

  ! external functions
  INTERFACE
     ! Pass logical parameter to dycore
     SUBROUTINE dycorewrapper_set_boolparam(logicalvalue, paramname_length, paramname) &
          BIND(c, name='dycorewrapper_set_boolparam')
       USE, INTRINSIC :: iso_c_binding
       LOGICAL(KIND=C_BOOL), value :: logicalvalue
       INTEGER(C_INT), value       :: paramname_length
       CHARACTER(KIND=C_CHAR), DIMENSION(*)  :: paramname
     END SUBROUTINE dycorewrapper_set_boolparam
  END INTERFACE
  
  ! locals
  INTEGER         :: paramname_length
  LOGICAL(KIND=1) :: l_logicalvalue

  IF (izdebug > 20) THEN
    WRITE(*,*) 'INFO-DBG: set_dycore_param_logical for '//TRIM(paramname), logicalvalue
  ENDIF

  paramname_length = LEN_TRIM(paramname)
  ! format convertion from LOGICAL to LOGICAL*1
  l_logicalvalue=logicalvalue
  CALL dycorewrapper_set_boolparam(l_logicalvalue, paramname_length, paramname)

END SUBROUTINE set_dycore_param_logical

SUBROUTINE set_dycore_param_string(paramname, stringvalue)
  ! arguments
  CHARACTER (LEN=*), INTENT(IN) :: paramname
  CHARACTER (LEN=*), INTENT(IN) :: stringvalue

  ! external functions
  INTERFACE
     ! Pass string parameter to dycore
     SUBROUTINE dycorewrapper_set_strparam(stringvalue_length, stringvalue, paramname_length, paramname) &
          BIND(c, name='dycorewrapper_set_strparam')
       USE, INTRINSIC :: iso_c_binding
       INTEGER(C_INT), value       :: stringvalue_length
       CHARACTER(KIND=C_CHAR), DIMENSION(*)  :: stringvalue
       INTEGER(C_INT), value       :: paramname_length
       CHARACTER(KIND=C_CHAR), DIMENSION(*)  :: paramname
     END SUBROUTINE dycorewrapper_set_strparam
  END INTERFACE
  
  ! locals
  INTEGER :: stringvalue_length
  INTEGER :: paramname_length

  stringvalue_length = LEN_TRIM(stringvalue)
  paramname_length = LEN_TRIM(paramname)
  CALL dycorewrapper_set_strparam(stringvalue_length, stringvalue, paramname_length, paramname)

END SUBROUTINE set_dycore_param_string


!==============================================================================
!+ Module procedure in "src_cpp_dycore" for organizing the time stepping  
!------------------------------------------------------------------------------

SUBROUTINE cpp_dycore_compute

!------------------------------------------------------------------------------
!
! Description: This subroutine executes C++ dycore functions
!      to run a C++ Dycore time step
!
! Method: 1. Copy input prognostic fields from Fortran to C++ Dycore
!         2. Run a time step of the C++ Dycore
!         3. Copy output prognostic fields from C++ Dycore to Fortran fields. 
!         4. Swap the pointers kept in data structures in the C++ Dycore, to
!            keep C++ Dycore fields synchronized with the Fortran fields 
!            (which will be swapped later one)
!         In debug mode, call a C++ Dycore function that checks for same
!         pointers in Fortran fields and those of the C++ Dycore.
!
!------------------------------------------------------------------------------

IMPLICIT NONE

! external functions
!-------------------------- 
 INTERFACE

    ! Reset meters
    SUBROUTINE dycorewrapper_reset_meters()  &
         BIND(c, name='dycorewrapper_reset_meters')
    END SUBROUTINE dycorewrapper_reset_meters

    ! CopyIn PrognosticFields
    SUBROUTINE dycorewrapper_copyin_prognosticfields()  &
         BIND(c, name='dycorewrapper_copyin_prognosticfields')
    END SUBROUTINE dycorewrapper_copyin_prognosticfields

    ! Do a dycore step
    SUBROUTINE dycorewrapper_do_dycorestep(ntstep, nlastbound) &
         BIND(c, name='dycorewrapper_do_dycorestep')
      USE, INTRINSIC :: iso_c_binding
      INTEGER(C_INT), value   :: ntstep, nlastbound
    END SUBROUTINE dycorewrapper_do_dycorestep

    ! CopyOut PrognosticFields
    SUBROUTINE dycorewrapper_copyout_prognosticfields()  &
         BIND(c, name='dycorewrapper_copyout_prognosticfields')
    END SUBROUTINE dycorewrapper_copyout_prognosticfields
    
 END INTERFACE

  ! Local scalars:
  ! -------------

  INTEGER :: k, izerror, iztrcr

  LOGICAL :: hasSurfField, hasTensField

  ! End of header
  !==============================================================================

  !------------------------------------------------------------------------------
  ! Begin module procedure cpp_dycore_compute
  !------------------------------------------------------------------------------

  IF (ltime) CALL get_timings (i_dyn_computations, ntstep, dt, izerror)

  IF ( ntstep == nstart + 1) THEN
    CALL dycorewrapper_reset_meters()
  ENDIF

  ! *********** If there is new bd data, swap pointers and copy nbd2 into the C++ dycore
  IF ( (ntstep+1 > nnextbound_local) .AND. (ntstep < nstop) ) THEN
    IF (izdebug > 8) THEN
      WRITE(*,*) 'INFO: swapping and copying new boundary fields at nbd2 =',nbd2
    ENDIF
    CALL dycore_swap_bdpointers()
    nnextbound_local = nnextbound
  ENDIF

  ! *********** copy in data for C++ dycore
  CALL dycorewrapper_copyin_prognosticfields()
  IF (ltime) CALL get_timings (i_cppdycore_copyin, ntstep, dt, izerror)

  ! *********** CALL timestepping of C++ dycore
  CALL dycorewrapper_do_dycorestep(ntstep, nlastbound)
  IF (ltime) CALL get_timings (i_cppdycore_step, ntstep, dt, izerror)

  ! *********** copy out data for C++ dycore
  CALL dycorewrapper_copyout_prognosticfields()
  IF (ltime) CALL get_timings (i_cppdycore_copyout, ntstep, dt, izerror)

#ifdef CPP_DEBUG
#ifndef _OPENACC

  IF (my_cart_id==0) PRINT *,'      C++ DYCORE FORTRAN: Verifying pointers'

  !$acc data &
  !$acc present(hhl, rho0, p0, p0hl, t0hl, dt0dz, t0,  sqrtg_r_s,  sqrtg_r_u,  sqrtg_r_v,  sqrtg_r_w)     &
  !$acc present(hd_mask, u, v,  w, t,  pp,  rho,  utens,  vtens,  wtens,  ttens,  pptens )                &
  !$acc present(qrs,  tkvm,  tkvh,  fc,  tcm,  tch,  ps, t_g, crlat)                                      &
  !$acc present(dqvdt , lhfl_s, shfl_s, umfl_s, vmhl_s)                                                   &
  !$acc present( a1t, a2t, acrlat,  tgrlat, qv, qvtens, qv_s, qc, qctens )                                &
  !$acc present( u_bd, v_bd, t_bd, pp_bd, qv_bd, qc_bd)                                                   &
  !$acc present_or_create(qi, qitens, qr)                                                                 &
  !$acc present_or_create(qs, qg)                                                                         &
  !$acc present_or_create(qi_bd, qr_bd, qs_bd, qg_bd)

  !$acc host_data &
  !$acc use_device(hhl, rho0, p0, p0hl, t0hl, dt0dz, t0,  sqrtg_r_s,  sqrtg_r_u,  sqrtg_r_v,  sqrtg_r_w)  &
  !$acc use_device(hd_mask, u, v,  w, t,  pp,  rho,  utens,  vtens,  wtens,  ttens,  pptens )             &
  !$acc use_device(qrs, tkvm,  tkvh,  fc,  tcm,  tch,  ps, t_g,  a1t,  a2t,  crlat)                       &
  !$acc use_device(dqvdt , lhfl_s, shfl_s, umfl_s, vmhl_s)                                                &
  !$acc use_device( acrlat,  tgrlat, qv, qvtens, qv_s, qc, qctens )                                       &
  !$acc use_device(u_bd, v_bd, t_bd, pp_bd, qv_bd, qc_bd)                                                 &
  !$acc use_device(qi, qitens, qr)                                                                        &
  !$acc use_device(qs, qg)                                                                                &
  !$acc use_device(qi_bd, qr_bd, qs_bd, qg_bd)

  CALL dycore_verify_pointers(hhl, rho0, p0, p0hl, t0hl, dt0dz, t0,                               &
        sqrtg_r_s, sqrtg_r_u, sqrtg_r_v, sqrtg_r_w, hd_mask,                                      &
        u(:,:,:,nnew), u(:,:,:,nnow), v(:,:,:,nnew), v(:,:,:,nnow), w(:,:,:,nnew), w(:,:,:,nnow), &
        t(:,:,:,nnew), t(:,:,:,nnow), pp(:,:,:,nnew), pp(:,:,:,nnow), rho,                        &
        utens, vtens, wtens, ttens, pptens(:,:,:), qrs(:,:,:),                                    &
        tkvm, tkvh, fc, tcm, tch,                                                                 &
        ps(:,:,nnew), t_g(:,:,nnew), ps(:,:,nnow), t_g(:,:,nnow),                                 &
        a1t, a2t, crlat(:,1), crlat(:,2), acrlat(:,1), acrlat(:,2),                               &
        trglat(:,1), trglat(:,2), rmy(:,:,1), rmy(:,:,2), rmy(:,:,3), rmyq,                       &
        u_bd(:,:,:,nbd1), u_bd(:,:,:,nbd2), v_bd(:,:,:,nbd1), v_bd(:,:,:,nbd2),                   &
        t_bd(:,:,:,nbd1), t_bd(:,:,:,nbd2), pp_bd(:,:,:,nbd1), pp_bd(:,:,:,nbd2),                 &
       dqvdt(:,:,:), qvsflx(:,:), lhfl_s(:,:), shfl_s(:,:), umfl_s(:,:), vmfl_s(:,:))

  !$acc end host_data
  !$acc end data
  
  ! Check tracers
  DO iztrcr = 1, trcr_get_ntrcr()
    CALL dycore_tracer_verify(iztrcr)
  END DO

#endif
#endif

  !----------------------------------------------------------------------------
  ! End of module procedure cpp_dycore_compute
  !----------------------------------------------------------------------------

END SUBROUTINE cpp_dycore_compute

!==============================================================================
!+ Module procedure in "src_cpp_dycore" for swapping pointers of all prognostic
! fields in the C++ dycore
!------------------------------------------------------------------------------

SUBROUTINE cpp_dycore_swap_pointers
  IMPLICIT NONE

  INTERFACE
     SUBROUTINE dycorewrapper_swap_fields()  &
          BIND(c, name='dycorewrapper_swap_fields')
     END SUBROUTINE dycorewrapper_swap_fields
  END INTERFACE

  INTEGER ::  izerror

  IF (izdebug > 20) THEN
    WRITE(*,*) 'INFO-DBG: cpp_dycore_swap_pointers'
  ENDIF

  CALL dycorewrapper_swap_fields()
  IF (ltime) CALL get_timings (i_cppdycore_swap, ntstep, dt, izerror)

END SUBROUTINE cpp_dycore_swap_pointers

!==============================================================================
!+ Module procedure in "src_cpp_dycore" for swapping pointers of boundary
! fields in the C++ dycore
!------------------------------------------------------------------------------

SUBROUTINE dycore_swap_bdpointers
  IMPLICIT NONE

  INTERFACE
     SUBROUTINE dycorewrapper_swapandcopy_bdfields()  &
          BIND(c, name='dycorewrapper_swapandcopy_bdfields')
     END SUBROUTINE dycorewrapper_swapandcopy_bdfields
  END INTERFACE

  IF (izdebug > 20) THEN
    WRITE(*,*) 'INFO-DBG: cpp_dycore_swap_bdpointers'
  ENDIF

  CALL dycorewrapper_swapandcopy_bdfields()

END SUBROUTINE dycore_swap_bdpointers

!==============================================================================
!+ Module procedure in "src_cpp_dycore" for checking compatibility between
!  C++ dycore and Fortran code
!------------------------------------------------------------------------------

SUBROUTINE dycore_check
  
  ! This routine verify the consistency between the dycore and the Fortran code
  ! The tests are performed in the dycore and retrun a check flag. If the flag value 
  ! is greater than 0 an error message is printed and an MPI_ABORD is called.
  ! Note that the format of integers passed here is fixed to int32

  USE, INTRINSIC :: iso_c_binding

  ! external functions
  INTERFACE

     SUBROUTINE dycorewrapper_check(interfaceVersion, intLength, floatLength, &
          logicTestTrue, logicTestFalse, kflat, checkFlag, errorstr, errorstrLength) &
! Different function name for cpu and gpu
! This ensures an error at linking time if one attemps to link a CPU Fortran 
! with GPU dycore and vice versa
#ifndef _OPENACC
       BIND(c, name='dycorewrapper_check_cpu')
#else 
       BIND(c, name='dycorewrapper_check_gpu')
#endif
       USE, INTRINSIC :: iso_c_binding
       ! in
       INTEGER(KIND=C_INT32_T),  value :: interfaceVersion, intLength, floatLength
       LOGICAL(KIND=C_BOOL), value :: logicTestTrue, logicTestFalse
       INTEGER(KIND=C_INT32_T),  value :: kflat
       ! out 
       INTEGER(KIND=C_INT32_T) :: checkFlag
       CHARACTER(KIND=C_CHAR), DIMENSION(256)  :: errorstr
       INTEGER(KIND=C_INT32_T) :: errorstrLength
     END SUBROUTINE dycorewrapper_check

  END INTERFACE
  
  ! locals

  ! Version number of the interface
  ! The value has to be kept manually synchronized between this file and the c++
  ! header file DycoreWrapper.h. It should be incremented for every changes that
  ! make the dycore and a previous version of the Fortran code incompatible.
  INTEGER(KIND=C_INT32_T), PARAMETER :: dycoreInterfaceVersion = 11

  INTEGER(KIND=C_INT32_T)  ::  intLength, floatLength
  LOGICAL(KIND=C_BOOL) :: logicTestTrue, logicTestFalse
  INTEGER(KIND=C_INT32_T)  :: ikflat
  INTEGER(KIND=C_INT32_T) :: checkFlag, errorstrLength
  CHARACTER(256) :: errorstr
  
! FUO TODO: sizeof is a GNU extension!
  intLength = sizeof(1)
  floatLength = sizeof(1.0_wp)
  
  ! use to check that Fortran and c logical are compatible
  logicTestTrue = .True.
  logicTestFalse = .False.

  ! check kflat
  ikflat=vcoord%kflat-1 !substract 1 because kflat is 0 based in the c++ dycore

  IF (izdebug > 20) THEN
    WRITE(*,*) 'INFO-DBG: dycorewrapper_check'
  ENDIF

  CALL dycorewrapper_check(dycoreInterfaceVersion, intLength, floatLength, &
          logicTestTrue, logicTestFalse, ikflat, checkFlag, errorstr, errorstrLength)
  
  ! WARNING
  IF (checkFlag < 0) THEN
    IF (my_cart_id == 0) THEN
      WRITE(*,*) '===== WARNING from C++ Dycore ===='
      WRITE(*,*) errorstr(1:errorstrLength)
      WRITE(*,*) '=================================='
    ENDIF
  END IF

  ! ERROR
  IF (checkFlag > 0) THEN
     CALL model_abort (my_world_id, 9062 , errorstr(1:errorstrLength), 'dycore_check')
  END IF

END SUBROUTINE dycore_check

!==============================================================================
!+ Module procedure in "src_cpp_dycore" for cleaning up
!------------------------------------------------------------------------------

SUBROUTINE cpp_dycore_cleanup

!------------------------------------------------------------------------------
!
! Description: This subroutine executes C++ dycore wrapper functions
!      to cleanup the C++ Dycore
!
! Method: 1. Output timers
!         2. Cleanup
!
!------------------------------------------------------------------------------

IMPLICIT NONE

 INTERFACE

    ! Print meters
    SUBROUTINE dycorewrapper_print_meters()  &
         BIND(c, name='dycorewrapper_print_meters')
    END SUBROUTINE dycorewrapper_print_meters
    
    ! Finalize dycore
    SUBROUTINE dycorewrapper_finalize()  &
         BIND(c, name='dycorewrapper_finalize')
    END SUBROUTINE dycorewrapper_finalize

 END INTERFACE

! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin module procedure cpp_dycore_cleanup
!------------------------------------------------------------------------------

  IF (izdebug > 20) THEN
    WRITE(*,*) 'INFO-DBG: cpp_dycore_cleanup'
  ENDIF

  IF (ltime) THEN
    CALL dycorewrapper_print_meters()
  ENDIF
  CALL dycorewrapper_finalize()

!----------------------------------------------------------------------------
! End of module procedure cpp_dycore_cleanup
!----------------------------------------------------------------------------

END SUBROUTINE cpp_dycore_cleanup

#endif

END MODULE src_cpp_dycore

