!+ Source module for the interface between organize_physics and turbulence
!----------------------------------------------------------------------------------

MODULE turb_interface

!-------------------------------------------------------------------------------
!
! Description:
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  Ulrich.Schaettler@dwd.de
!
! History:
! Version      Date       Name
! ----------   ---------- ----
! V5_4a        2016-05-10 Ulrich Schaettler
!  Initial release
! V5_4b        2016-07-12 Ulrich Schaettler, Ulrich Blahak
!  Comment out the use of tracer, which does not work with Cray CONTIGUOUS at
!  the moment. (US)
!  Implement turb_finalize for deallocation work arrays when running with OPENACC. (US)
!  Added SR gen_H0noise() and hnoise_b for idealized surface fluxes for the 
!    new blocked transfer scheme code. (UB)
!
! Code Description:
! Language: Fortran 90 
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Modules used:

USE kind_parameters, ONLY :   &
    wp           ! KIND-type parameter for real variables

!-------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------

    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction
    ke1,          & ! ke + 1
    kcm,          & ! index of the lowest model layer higher than the canopy

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
    istartu,      & ! start index for the forecast of u
    iendu,        & ! end index for the forecast of u
    istartv,      & ! start index for the forecast of v
    iendv,        & ! end index for the forecast of v
    istartpar,    & ! start index for computations in the parallel program
    iendpar,      & ! end index for computations in the parallel program

!   meridional direction
    jstart,       & ! start index for the forecast of w, t, qd, qw and pp
    jend,         & ! end index for the forecast of w, t, qd, qw and pp
    jstartu,      & ! start index for the forecast of u
    jendu,        & ! start index for the forecast of u
    jstartv,      & ! start index for the forecast of v
    jendv,        & ! end index for the forecast of v
    jstartpar,    & ! start index for computations in the parallel program
    jendpar,      & ! end index for computations in the parallel program

! 4. constants for the horizontal rotated grid and related variables
! ------------------------------------------------------------------

!   dlam,         & ! grid point distance in zonal direction (in degrees)
!   dphi,         & ! grid point distance in meridional direction (in degrees)
    eddlon,       & ! 1 / dlon
    eddlat,       & ! 1 / dlat
    edadlat,      & ! 1 / (radius of the earth * dlat)

! 5. variables for the time discretization and related variables
! --------------------------------------------------------------

    dt,           & ! long time-step
    dt2,          & ! 2 * dt
    dtdeh,        & ! dt / 3600 seconds
    lalloc_h_ice, & !

! 8. Organizational variables to handle the COSMO humidity tracers
! ----------------------------------------------------------------
    idt_qv, idt_qc, idt_qi

! end of data_modelconfig

!-------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 2. physical constants and related variables
! -------------------------------------------
    ! not yet necessary
    lh_v,         & ! evaporation heat
    lh_s,         & ! sublimation heat
    t0_melt,      & ! absolute freezing temperature
    r_earth         ! mean radius of the earth

! end of data_constants

!-------------------------------------------------------------------------------

USE data_fields     , ONLY :   &
    tke,            & ! SQRT(2*TKE); TKE='turbul. kin. energy'        ( m/s )
!   gz0,            & ! roughness length * g of the vertically not
!                     ! resolved canopy                               (m2/s2)
!   fr_land,        & ! land portion of a grid point area             ( 1 )
!   plcov,          & ! fraction of plant cover                       ( 1 )
!   depth_lk ,      & ! lake depth                                    ( m )
!   sai,            & ! surface area index                            ( 1 )
!   lai,            & ! leaf area index                               ( 1 )
!   tai,            & ! transpiration area index                      ( 1 )
!   eai,            & ! (evaporative) earth area index                ( 1 )
!   qvsflx     ,    & ! vater vapor   flux density                    ( Kg/m2s)
    h_ice      ,    & ! ice thickness                                 (  m  )
    hhl        ,    & !
    u_m, v_m   ,    & !
    u, v, w    ,    & !
    acrlat     ,    & !
    tkred_sfc  ,    & ! reduction factor for minimum diffusion coefficients near the surface
    l_hori     ,    & ! horizontal grid spacing (location dependent horizontal scale)
    hdiv       ,    & ! horizontal divergence
    hdef2      ,    & ! horizontal deformation square
    dwdx       ,    & ! horizontal gradients of vertical wind
    dwdy       ,    & ! horizontal gradients of vertical wind
    h0noise           ! random noise on idealized surface fluxes

!-------------------------------------------------------------------------------

USE data_block_fields  !, ONLY :   &

!-------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. Time indices                  
! ---------------                  
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
    nold,         & ! corresponds to ntstep - 1
    nnow,         & ! corresponds to ntstep
    nnew,         & ! corresponds to ntstep + 1
                    ! indices for permutation of the tke time levels
    ntke,         & ! time step index of current TKE values
                    ! (corresponds to nold for "leap-frog" or nnow for 2 time levels)
    ninctura,     & ! time step increment for running the vertical diffusion
    l3dturb,      & ! 3D-turbulence: CALL explicit_horizontal_diffusion (RK)
    lartif_data,  & ! forecast with self-defined artificial data
    l2tls,        & ! time integration by two timelevel RK-scheme (.TRUE.)
                    ! or by default three-time level KW-scheme (.FALSE.)
    nproma,       & ! block size for physical parameterizations

! 3. controlling the physics
! --------------------------
    itype_tran,   & ! type of surface-atmosphere transfer
    itype_turb,   & ! type of turbulent diffusion parametrization
    itype_synd,   & ! type of diagnosis of near surface synop. station values
    itype_vdif,   & ! type of vertical diffusion calculation
                    ! < 0: vertical diffusion calculated by the specific COSMO-routines in 'slow_tendencies'
                    !  -2: apply former (non-blocked) COSMO-version of Raschendorfer-scheme
                    !      (only temporary and meaningful, if "itype_turb=3")
                    ! >=0: vertical diffusion calculated by the common COSMO/ICON-routine contained in
                    !      in MODULE 'trubulence_turbdiff' (not yet completely implemented, only if "itype_turb=3")
                    !   0: vertical diffusion within the CALL of the turbulence model (ICON-mode)
                    ! > 0: vertical diffusion CALLed at the end of physics-section and considering explicit
                    !      tendencies in the semi-implicit solution.
    imode_mdif,   & ! mode of momentum diffusion (in case of "itype_vdif>0")
                    !   1: extra call for diffusion of each wind-component applied to the staggered wind positions
                    !   0: diffusion of momentum together with that of scalars on horiz. mass-positions
                    !      and never using any explicit wind-tendencies in implicit diffusion equation
                    !  -1: as " 0", but using the expl. wind-tend. already present at mass-positions (if requested)
                   !!! -2: as "-1", but using even the staggered wind-tend. (after horiz. interp. and if requested)
                   !!!     Notice that this option is not implemented (yet), as it seems really dispensable!!!
    lprintdeb_all, idbg_level, ldebug_tur

USE turb_data, ONLY : &

    modvar,       & ! structure type containing model variable data
    nmvar,        & ! = nscal+nvel
    itype_sher,   & ! type of shear production for TKE
    imode_pat_len,& ! mode of determining the length scale of surface patterns (related to 'pat_len')
    ltkeshs,      & ! consider separ. horiz. shear production of TKE
    loutshs,      & ! consider separ. horiz. shear production of TKE for output
    turb_wkarr_dealloc

! end of data_runcontrol 

!-------------------------------------------------------------------------------

#ifndef SCLM
USE src_artifdata,      ONLY:  tkvhfix, tkhhfix, tkvmfix, tkhmfix, &
                               lnosurffluxes_m, lnosurffluxes_h, lsensiflux_fix, &
                               set_idealized_surffluxes_block, gen_H0noise
#else
USE data_1d_global, ONLY : lsclm, im,jm
#endif

!-------------------------------------------------------------------------------

USE data_parallel,      ONLY :  &
    num_compute,     & ! number of compute PEs
    my_cart_neigh,   & ! neighbors of this subdomain in the cartesian grid
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    ldatatypes,      & ! if .TRUE.: use MPI-Datatypes for some communications
    ltime_barrier,   & ! if .TRUE.: use additional barriers for determining the
                       ! load-imbalance
    ncomm_type,      & ! type of communication
    icomm_cart,      & ! communicator for the virtual cartesian topology
    iexch_req,       & ! stores the sends requests for the neighbor-exchange
                       ! that can be used by MPI_WAIT to identify the send 
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains
    nexch_tag,       & ! tag to be used for MPI boundary exchange
                       !  (in calls to exchg_boundaries)
    sendbuf,         & ! sending buffer for boundary exchange:
                       ! 1-4 are used for sending, 5-8 are used for receiving
    isendbuflen        ! length of one column of sendbuf

!-------------------------------------------------------------------------------

! for blocked version of physics
USE turb_transfer,  ONLY : turbtran
USE turb_diffusion, ONLY : turbdiff
USE turb_utilities, ONLY : turb_param

!-------------------------------------------------------------------------------

USE src_tracer,       ONLY : trcr_get, trcr_errorstr, trcr_meta_get, &
                             trcr_get_ntrcr !number of all tracers

USE data_tracer,      ONLY : T_TURB_ID, T_TURB_1D, T_BBC_ID, &
                             T_BBC_ZEROFLUX, T_BBC_ZEROVAL, T_BBC_SURF_VAL

!-------------------------------------------------------------------------------

USE src_block_fields, ONLY  :  register_copy, CopylistStruct,               &
                               copyToBlockF, copyFromBlockF,                &
                               init_copy_list  !, mind_ilon, mind_jlat

!-------------------------------------------------------------------------------

USE environment,      ONLY :   model_abort, exchg_boundaries

!-------------------------------------------------------------------------------
! Declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE

PUBLIC :: turb_init, turb_prepare, turb_organize

!-------------------------------------------------------------------------------

! local variables for turb_prepare
REAL (KIND=wp), ALLOCATABLE :: lays(:,:,:), vgrd(:,:,:)

REAL (KIND=wp), PARAMETER :: &
  z1  =1.0_wp, &
  z2  =2.0_wp, &
  z1d2=0.5_wp, &
  z1d4=0.25_wp

INTEGER              :: &
  ntimtke,              & ! number of TKE time levels
  ntr_totl,             & ! total number of tracers available
  ntr_diff                ! number of tracers for vertical diffusion

INTEGER, ALLOCATABLE :: &
  itr_diff(:), & ! index array for tracers dedicated for vertical diffusion
  izturb  (:), & ! index for turb mix for all tracers
  izbbc   (:)    ! bottom BC for all tracers

TYPE (modvar), ALLOCATABLE :: ptr(:) !structure for passive tracers being diffused

TYPE (CopylistStruct) :: turCopyList

!===============================================================================

CONTAINS

!===============================================================================

!+ turbulence_organize calls the different schemes depending on options set
!  it works on the blocked version of the turbulence

SUBROUTINE turb_organize (ib, ipend, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!  This routine organizes the calls to the different turbulence schemes and
!  acts as interface between the COSMO-Model (routine organize_physics) and
!  the different turbulence packages. It uses data from the corresponding
!  COSMO-Modules and passes these on the parameterizations.
!
! Method:
!
!  Implementing the CALL of SUB 'organize_turbdiff' as it is already used in ICON,
!   which is (according to the ICON rules) CALLed by a long argmument list,
!   uses only a single horizontal array-index and doesn't contain any
!   "horizontal operations", which now are contained in SUB 'horiz_wind_calcul'.
!  SUB 'organize_turbdiff' and all related routines are contained in the single
!   MODULE 'src_turbdiff' for easy exchange between COSMO and ICON.
!  SUB 'organize_turbdiff' contains the following actions: 
!   surface-atmospher transfer, atmospheric turbulence and vertical diffusion,
!   and it can be CALLed for any combination of these three actions.
!  Particularly it contains new routines for semi-implicit vertical diffusion, which
!   are employed in ICON just with the CALL of SUB 'turbdiff' in case of "itype_vdif=0".
!  Otherwise, the vertical diffusion is employed later, after all physical tendencies
!   are collected. For the time being, this remains within the 'slow_tendencies'-routines.
!  In the case "itype_vdif>0" 'apply_turbdiff' is called in 'slow_tendencies' only for
!   the calculation of vertical diffusion, where 'itype_vdif' controls the mode of considering 
!   the already updated (explicit) tendencies in the diffusion calculation.
!  In the case "itype_vdif<=-1", the previous routines for vertical diffusion 
!   in 'slow_tendencies' are used.
!  In the (only temporary) case "itype_vdif=-2", the former SUB 'organize_turbulence' is used
!   instead of 'turbulence_organize' and the former (non-blocked) COSMO-version of turbulence
!   is applied!
!  In case of "itype_vdif=0" no explicit tendencies are considered.
!  By use of SUB 'organize_turbdiff' also arbitrary passive tracers can be included into the
!   calculation of vertical diffusion by use of the data-structure 'ptr' of TYPE 'modvar',
!   which allows to include all necessary information particularly of defining the lower boundary.
!  Compact calls of 'apply_canopy_init' and 'apply_turbdiff', avoiding the explicit
!    editing of long argument lists
!   being more comfortabe when called several times in the code.
!  As SUB 'organize_turbdiff' (and also 'init_canopy') may be called for different actions,
!   the preparing compact CALLs 'apply_turbdiff (and 'apply_canopy_init') can be used,
!   in order to avoid editing and evaluating the long argument lists several times.
!  The meaning of 'imode_turbdif' and 'imode_turbtran' has changed and refers to the mode
!   of solving the TKE-equation now:
!   0: diagnostic; 1: prognostic (previous version); 2: prognostic (new positive def. version)
!  The lower boundary condition for vertical diffusion is now controlled by 'lsflcnd': 
!   .TRUE.: explicit flux condition; .FALSE.: concentration condition 
!  Preliminarly the facility for copying horizontal arrays to the 1D block-data structure
!   using a dynamical estimated block length is not used. Rather 1D partial arrays of the
!   zonal dimension are used as block-vectors.
!  The application of the case 'lartif_data' has been slightly adapted.
!
!------------------------------------------------------------------------------

! Formal Parameters:
INTEGER, INTENT(IN)  :: ib                 ! current block index
INTEGER, INTENT(IN)  :: ipend              ! length of current block

INTEGER, INTENT(OUT)           ::                      &
  ierror       ! error status

CHARACTER (LEN= *),       INTENT(OUT)           ::                      &
  yerrmsg      ! error message

!------------------------------------------------------------------------------

! Local Variables:
INTEGER            :: &
  ndiff,              & ! number of 1-st order variables used in the turbulence model
  izdebug,            & ! for debug printout
  izerrstat,          & ! for error checking
  nvar,               & ! indicating time level for physics
  kzdims(24),         & ! vertical dimensions for fields to be exchanged
  iv, i,j,k,n, iztrcr   ! loop indices

INTEGER, SAVE      :: &
  nprv,               & ! previous time step index for tke
  nvor,               & ! and what is passed to the routines
  iini                  ! type of initialization for turbulence scheme

LOGICAL            :: &
  lstfnct               ! calculation of a new stability function in sub. turbdiff

REAL (KIND=wp)     :: &
  zdt,                & ! time step to be used (depends on Leapfrog / Runge-Kutta)
  zfs_ice               !frozen fraction of the surface

!US REAL (KIND=wp), TARGET :: &
!US   zut_turb(nproma,ke),     & ! turbulenct u-tendency
!US   zvt_turb(nproma,ke)        ! turbulenct v-tendency

REAL (KIND=wp), POINTER  :: &
! qv     (:,:,:) => NULL(), & ! QV at nvar
! qv_tens(:,:,:) => NULL(), & ! QV tendency
! qc     (:,:,:) => NULL(), & ! QC at nvar
! qc_tens(:,:,:) => NULL(), & ! QC tendency
! Note: qi is not necessary, unless ice clouds are considered in the turbulence model!
  v2d    (:,:)   => NULL(), & ! dummy for a 2D variable array
  v3d    (:,:,:) => NULL()    ! dummy for a 3D variable array

REAL (KIND=wp), POINTER  :: &
  ! pointer to velocity components and their tendencies
  u_pt(:,:), utens_pt(:,:), v_pt(:,:), vtens_pt(:,:)

REAL (KIND=wp), POINTER  :: &
  ! pointer to surface variables
  tvm_pt(:), t_g_pt(:), qv_s_pt(:), ps_pt(:) 

REAL (KIND=wp), POINTER  :: &
  ! pointer to atmospheric variables
  tkvm_pt(:,:), rho_pt(:,:), hhl_pt(:,:), dp0_pt(:,:) 

REAL (KIND=wp), POINTER  :: &
  h_ice_p(:,:)

CHARACTER(LEN=250) :: &
  yerrormsg             ! for error strings

CHARACTER(LEN=250) :: &
  yroutine              ! name of this routine

!-End of Header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  ierror       = 0
  izerrstat    = 0
  yroutine     = 'turb_interface'
  yerrmsg(:)   = ' '
  yerrormsg(:) = ' '

  IF (ldebug_tur) THEN
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

  ! get the correct timelevel for dynamic variables
  IF (l2tls) THEN
    nvar = nnow
    zdt  = dt
  ELSE
    nvar = nold
    zdt  = dt2
  ENDIF

  ! check, whether atmospheric turbulence / stability function has to be 
  ! computed
  IF ( (ntstep < 10) .OR. ( MOD(ntstep+1,ninctura) == 0 ) ) THEN
    ! calculate turbulence or at least stability functions
    lstfnct = .TRUE.
  ELSE
    ! don't calculate stability functions or even turbulence
    lstfnct = .FALSE.
  ENDIF

!------------------------------------------------------------------------------
! Surface transfer and turbulent diffusion schemes:
!------------------------------------------------------------------------------

  IF (izdebug > 5 .AND. ib==1) THEN
    PRINT *, '  new TURBULENCE with itype_tran, itype_turb, lstfnct = ', itype_turb, itype_tran, lstfnct, ib
  ENDIF
  
  IF (itype_tran.EQ.2 .AND. itype_turb.EQ.3) THEN
    ! Default scheme for surface-to-atmosphere transfer and atmospheric turbulence by Raschendorfer 
    ! based on a generalized prognostic TKE equation and turbulent budget equations for quasi 
    ! conservative variables:
   
    ! determine initialization for tke: only for first block
    IF (ib == 1) THEN
      IF (ntke == 0) THEN ! the very fist time step
        ntke=ntimtke
        iini=2 ! initialization within this time step
      ELSE ! any other time step
        iini=0 ! no initialization
      ENDIF

      ! independent 2-timelevel stepping for turbulence model
      nprv = ntke
      ntke = MOD(ntke, ntimtke)+1
    ENDIF

    ! get the tracers: QV, QC + tendencies
    !US: we do not check for other tracer that should be vertically diffused here:

! in the blocked version we have the qv_b, qvtens_b, qc_b, qctens_b available

    ! here: itnd=0    don't consider explic. tendenc. in implic. vert. diff.-calc.
    ! No consideration of explicit horizontal wind tendencies for vertical diffusion
    ! set full pressure zpres
    DO k=1, ke
      DO iv = 1, ipend
        ut_turb_b(iv,k) = 0.0_wp
        vt_turb_b(iv,k) = 0.0_wp
        ptot_b   (iv,k) = p0_b(iv,k) + pp_b(iv,k)
      END DO
    END DO

    ! loop over all tracers to be diffused
    ! How to do this in the block structure?
    ! Up to now we do not have a tracer structure in blocked form
    ! only the microphysics tracers are there as allocatable fields
    ! AND: right now only cloud ice QI is needed

    DO n=1, ntr_diff
      iztrcr=itr_diff(n) !specific tracer index

      ! as long as we do not have a specific tracer structure in blocked version,
      ! only cloud ice can be treated here (QR, QS, QG are not vertically diffused)
      IF (iztrcr == idt_qi) THEN

        ptr(n)%av => qi_b(:,:)
        ptr(n)%at => qitens_b(:,:)
        IF (izbbc(iztrcr) == T_BBC_SURF_VAL) THEN
          ! should not be the case for QI!!!
          !USptr(n)%sv => v2d(:,j)
          ptr(n)%sv => NULL()
          ptr(n)%fc=.FALSE.   !so far only concentration surface variables active for tracers
        ELSE
          ptr(n)%sv => NULL() !use zero surface value
          IF (izbbc(iztrcr) == T_BBC_ZEROVAL) THEN
            ptr(n)%fc=.FALSE. !surface value is a concentration
          ELSEIF (izbbc(iztrcr) == T_BBC_ZEROFLUX) THEN
            ptr(n)%fc=.TRUE.  !surface value is a flux density
          ELSE
            WRITE (yerrormsg,'(I3)') iztrcr
            yerrormsg = 'ERROR: NOT VALID BBC-OPTION FOR TRACER NUMBER '//TRIM(yerrormsg)
            CALL model_abort(my_cart_id, izerrstat, yerrormsg, yroutine)
          ENDIF
        ENDIF
      ELSE
!US this case is not possible at the moment!!!
!   we need a modification in the (block) tracer module for that
        ! retrieve the required passive tracer for vertical diffusion at specified time level:
!       CALL trcr_get(izerrstat, iztrcr, ptr_tlev = nvar, ptr = v3d)
!       IF (izerrstat /= 0) THEN
!         yerrormsg = trcr_errorstr(izerrstat)
!         CALL model_abort(my_cart_id, izerrstat, yerrormsg, yroutine)
!       ENDIF
!       ptr(n)%av => v3d(:,j,:)

!       CALL trcr_get(izerrstat, iztrcr, ptr_tens = v3d)
!       IF (izerrstat /= 0) THEN
!         yerrormsg = trcr_errorstr(izerrstat)
!         CALL model_abort(my_cart_id, izerrstat, yerrormsg, yroutine)
!       ENDIF
!       ptr(n)%at => v3d(:,j,:)

!       IF (izbbc(iztrcr) == T_BBC_SURF_VAL) THEN
!         CALL trcr_meta_get(izerrstat, iztrcr, 'SURF_FIELD', v2d)
!         IF (izerrstat /= 0) THEN
!           yerrormsg = trcr_errorstr(izerrstat)
!           CALL model_abort(my_cart_id, izerrstat, yerrormsg, yroutine)
!         ENDIF
!         ptr(n)%sv => v2d(:,j)
!         ptr(n)%fc=.FALSE.   !so far only concentration surface variables active for tracers
!       ELSE
!         ptr(n)%sv => NULL() !use zero surface value
!         IF (izbbc(iztrcr) == T_BBC_ZEROVAL) THEN
!           ptr(n)%fc=.FALSE. !surface value is a concentration
!         ELSEIF (izbbc(iztrcr) == T_BBC_ZEROFLUX) THEN
!           ptr(n)%fc=.TRUE.  !surface value is a flux density
!         ELSE
!           WRITE (yerrormsg,'(I3)') iztrcr
!           yerrormsg = 'ERROR: NOT VALID BBC-OPTION FOR TRACER NUMBER '//TRIM(yerrormsg)
!           CALL model_abort(my_cart_id, izerrstat, yerrormsg, yroutine)
!         ENDIF
!       END IF
      ENDIF ! iztrcr==idt_qi
    END DO

#ifdef SCLM
    IF (j.EQ.jm) THEN !this data-line contains the central point of SCLM
      lsclm=.TRUE. !apply SCLM-Modifications
      !Note: This needs to be adapted, if the 1D data-lines are no longer
      !      j-vectors. Then also the value of 'im' may vary and it has
      !      to be reset after the data-line loop again.
    ELSE !this data-line does not contain the central point of SCLM
      lsclm=.FALSE. !don't apply SCLM-modifications
    END IF
#endif

! ptr is allocated in turb_init and available here
  ndiff=nmvar+ntr_diff !number of 1-st order variables used in the turbulence model
                       !note that cloud ice is treated like a passive trace here

!US now pass nvor instead of nprv, to have the same value in turbtran and turbdiff
 nvor=nprv !Eingangsbelegung von 'nvor' (wird bei Iterationen auf 'ntur' gesetzt)

    CALL turbtran (                                                  &
       iini=iini, lstfnct   = lstfnct,                               &
                  lnsfdia   = (itype_synd.EQ.2),                     &
                  ltkeinp   = .FALSE.,                               &
                  lgz0inp   = .FALSE.,                               &
                  lsrflux   = (itype_vdif==0),                       &
                  dt_tke    = dt,                                    &
                  nvor      = nvor,                                  &
                  ntur      = ntke,                                  &
                  ntim      = ntimtke,                               &
                  nvec      = nproma,                                &
                  ke        = ke,                                    &
                  ke1       = ke1,                                   &
                  kcm       = kcm,                                   &
                  iblock    = ib,                                    &
                  ivstart   = 1,                                     &
                  ivend     = ipend,                                 &
                  l_hori    = l_hori_b(:),                           &
                  hhl       = hhl_b(:,:),                            &
                  fr_land   = fr_land_b(:),                          &
                  depth_lk  = depth_lk_b(:),                         &
                  h_ice     = h_ice_b(:),                            &
                  gz0       = gz0_b(:),                              &
                  sai       = sai_b(:),                              &
                  t_g       = t_g_b(:),                              &
                  qv_s      = qv_s_b(:),                             &
                  ps        = ps_b(:),                               &
                  u         = u_m_b(:,:),                            &
                  v         = v_m_b(:,:),                            &
                  w         = w_b(:,:),                              &
                  t         = t_b(:,:),                              &
                  qv        = qv_b(:,:),                             &
                  qc        = qc_b(:,:),                             &
                  prs       = ptot_b(:,:),                           &
!                 epr       = not available
                  tcm       = tcm_b(:),                              &
                  tch       = tch_b(:),                              &
                  tvm       = tvm_b(:),                              &
                  tvh       = tvh_b(:),                              &
                  tfm       = tfm_b(:),                              &
                  tfh       = tfh_b(:),                              &
                  tfv       = tfv_b(:),                              &
                  tkr       = tkr_b(:),                              &
                  tke       = tke_b(:,:,:),                          &
                  tkvm      = tkvm_b(:,:),                           &
                  tkvh      = tkvh_b(:,:),                           &
                  rcld      = rcld_b(:,:),                           &
                  hdef2     = hdef2_b(:,:),                          &
!                 hdiv      = not available
                  dwdx      = dwdx_b(:,:),                           &
                  dwdy      = dwdy_b(:,:),                           &
                  edr       = edr_b(:,:),                            &
                  tketens   = tketens_b(:,:),                        &
                  t_2m      = t_2m_b(:),                             &
                  qv_2m     = qv_2m_b(:),                            &
                  td_2m     = td_2m_b(:),                            &
                  rh_2m     = rh_2m_b(:),                            &
                  u_10m     = u_10m_b(:),                            &
                  v_10m     = v_10m_b(:),                            &
!Achtung: not yet activated
                  shfl_s    = shfl_s_b(:),                           &
!                 lhfl_s    = not available
                  qvfl_s    = qvsflx_b(:),                           &
!                 umfl_s    = not available
!                 vmfl_s    = not available
                  ierrstat=izerrstat, yerrormsg=yerrormsg, yroutine=yroutine)

#ifndef SCLM
  ! special option for idealized test cases
  IF (lartif_data) THEN

    ! Specify tch, t_g and qv_s in a way that predefined surface
    ! sensible and latent heat fluxes result, based on
    ! some NAMELIST parameters:
    IF (lsensiflux_fix) THEN
      IF (itype_vdif < 0) THEN
        CALL set_idealized_surffluxes_block (1, ipend, ke, t_b, ptot_b, qv_b, u_m_b, v_m_b, ps_b, h0noise_b, &
             t_g_b, t_g_new_b, t_s_b, t_s_new_b, qv_s_b, qv_s_new_b, tch_b)
      ELSE
        PRINT *,' ERROR    *** lsensiflux_fix=.TRUE. only possible if itype_vdif < 0 *** '
        ierror = 2798
        RETURN
      END IF
    END IF

    ! Set surface diffusion coefficients and therefore surface fluxes 
    ! of momentum, heat and moisture to 0:
    IF (lnosurffluxes_m) THEN
      IF (itype_vdif < 0) THEN
        tcm_b(:) = 0.0_wp
      ELSE
        PRINT *,' ERROR    *** lnosurffluxes_m=.TRUE. only possible if itype_vdif < 0 *** '
        ierror = 2796
        RETURN
      END IF
    END IF
    IF (lnosurffluxes_h) THEN
      IF (itype_vdif < 0) THEN
        tch_b(:) = 0.0_wp
      ELSE
        PRINT *,' ERROR    *** lnosurffluxes_h=.TRUE. only possible if itype_vdif < 0 *** '
        ierror = 2797
        RETURN
      ENDIF
    END IF

  END IF
#endif


    CALL turbdiff (                                                  &
       iini=iini, lturatm   = (itype_turb==3),                       &
                  lstfnct   = lstfnct,                               &
                  ltkeinp   = .FALSE.,                               &
                  itnd      = 0,                                     &
                  lum_dif   = (itype_vdif==0),                       &
                  lvm_dif   = (itype_vdif==0),                       &
                  lscadif   = (itype_vdif==0),                       &
                  lsfluse   = .FALSE.,                               &
                  lqvcrst   = .FALSE.,                               &
                  dt_var    = zdt,                                   &
                  dt_tke    = dt,                                    &
                  nvor      = nvor,                                  &
                  ntur      = ntke,                                  &
                  ntim      = ntimtke,                               &
                  nvec      = nproma,                                &
                  ke        = ke,                                    &
                  ke1       = ke1,                                   &
                  kcm       = kcm,                                   &
                  iblock    = ib,                                    &
                  ivstart   = 1,                                     &
                  ivend     = ipend,                                 &
                  l_hori    = l_hori_b(:),                           &
                  hhl       = hhl_b(:,:),                            &
                  fr_land   = fr_land_b(:),                          &
                  l_pat     = l_pat_b(:),                            &
                  dp0       = dp0_b(:,:),                            &
!                 trop_mask = not available                          &
                  gz0       = gz0_b(:),                              &
!                 c_big     = not available                          &
!                 c_sml     = not available                          &
!                 r_air     = not available                          &
                  t_g       = t_g_b(:),                              &
                  qv_s      = qv_s_b(:),                             &
                  ps        = ps_b(:),                               &
                  u         = u_m_b(:,:),                            &
                  v         = v_m_b(:,:),                            &
                  w         = w_b(:,:),                              &
                  t         = t_b(:,:),                              &
                  qv        = qv_b(:,:),                             &
                  qc        = qc_b(:,:),                             &
                  prs       = ptot_b(:,:),                           &
                  rho       = rho_b(:,:),                            &
!                 epr       = not available                          &
                  ptr       = ptr(1:ntr_diff),                       &
                  ndtr      = ntr_diff,                              &
                  ndiff     = ndiff,                                 &
                  tvm       = tvm_b(:),                              &
                  tvh       = tvh_b(:),                              &
                  tfm       = tfm_b(:),                              &
                  tfh       = tfh_b(:),                              &
                  tfv       = tfv_b(:),                              &
                  tkr       = tkr_b(:),                              &
                  tkred_sfc = tkred_sfc_b(:),                        &
                  tke       = tke_b(:,:,:),                          &
                  tkvm      = tkvm_b(:,:),                           &
                  tkvh      = tkvh_b(:,:),                           &
                  rcld      = rcld_b(:,:),                           &
                  tkhm      = tkhm_b(:,:),                           &
                  tkhh      = tkhh_b(:,:),                           &
                  hdef2     = hdef2_b(:,:),                          &
                  hdiv      = hdiv_b(:,:),                           &
                  dwdx      = dwdx_b(:,:),                           &
                  dwdy      = dwdy_b(:,:),                           &
                  edr       = edr_b(:,:),                            &
                  tket_sso  = tket_sso_b(:,:),                       &
                  tket_conv = tket_conv_b(:,:),                      &
                  tket_hshr = tket_hshr_b(:,:),                      &
                  u_tens    = ut_turb_b(:,:),                        &
                  v_tens    = vt_turb_b(:,:),                        &
                  t_tens    = ttens_b(:,:),                          &
                  qv_tens   = qvtens_b(:,:),                         &
                  qc_tens   = qctens_b(:,:),                         &
                  tketens   = tketens_b(:,:),                        &
                  tketadv   = tket_adv_b(:,:),                       &
                  qv_conv   = dqvdt_b(:,:),                          &
                  ut_sso    = ut_sso_b(:,:),                         &
                  vt_sso    = vt_sso_b(:,:),                         &
!Achtung: not yet activated
                  shfl_s    = shfl_s_b(:),                           &
!                 lhfl_s    = not available                          &
                  qvfl_s    = qvsflx_b(:),                           &
!                 umfl_s    = not available                          &
!                 vmfl_s    = not available                          &

                  ierrstat=izerrstat, yerrormsg=yerrormsg, yroutine=yroutine)

!US----------------------------------------------------------------------------------------------

    IF (izerrstat.NE.0) THEN
      CALL model_abort (my_cart_id, izerrstat, yerrormsg, yroutine)
    END IF

!Achtung: not yet activated
!  why is it in the code then????
    DO k = 1, ke
    DO iv = 1, ipend
!     ! Set dqvdt_b to 0.0 at the boundary of the domain
!     i = mind_ilon(iv,ib)
!     j = mind_jlat(iv,ib)
!     IF (i <= nboundlines        .OR.  j <= nboundlines .OR.                 &
!         i >= ie-nboundlines+1   .OR.  j >= je-nboundlines+1)   THEN
!       dqvdt_b (iv,k) = 0.0_wp
!     ENDIF

      ! Estimation of the frozen fraction of the surface:
      IF (t_snow_b(iv) < t0_melt) THEN
        zfs_ice=1.0_wp
      ELSEIF (t_snow_b(iv) /= t_s_b(iv)) THEN
        zfs_ice=MAX( 0.0_wp, MIN( 1.0_wp, &  !only for security
              (t_g_b(iv)-t_s_b(iv))/(t_snow_b(iv)-t_s_b(iv)) ) )
      ELSE
        zfs_ice=0.0_wp
      END IF
      lhfl_s_b(iv) = (zfs_ice*lh_s + (1.0_wp-zfs_ice)*lh_v) * qvsflx_b(iv)
     !lhfl_s_b(i,j)=qvsflx_b(i,j)*lh_v
    END DO
    END DO

#ifdef SCLM
    lsclm=.TRUE. !reset lsclm
#endif

     !Note:
     !If, in case of active surface tiles, the tile loop for surface transfer would not be
     ! carried out within SUB 'organize_turbdiff', SUB 'apply_turbdiff' can be CALLed only for 
     ! surface transfer above each tile ("ltursrf=T", "lturatm=F", "lum_dif=F", "lvm_dif=F", "lscadif=F").
     ! In this case, the surface transfer would be excluded here.
     !If semi-implicit vertical diffusion shall be executed considering all the explicit tendencies
     ! of other parameterizations, SUB 'apply_turbdiff' can be CALLed again later only for vertical
     ! diffusion ("ltursrf=F", "lturatm=F", "lum_dif=T", "lvm_dif=T", "lscadif=T").
     ! In this case, the vertical diffusion would be excluded here.

  ELSE

    PRINT *,' ERROR *** New blocked turbulence code only possible if both itype_tran == 2 and itype_turb == 3 *** '
    ierror = 3797
    RETURN

  ENDIF       ! IF (itype_tran.EQ.2 .OR. itype_turb.EQ.3) THEN
  

#ifndef SCLM
!MR: The concept of forcing the model by predefined surface fluxes or diffusion coefficients
!     needs to be reconsidered.

  ! special option for idealized test cases
  IF (lartif_data) THEN
    IF (itype_turb.EQ.100) THEN
      ! Set constant diffusion coefficients from namelist parameters:
      WRITE (*,'(a,3(es8.1,","),es8.1,a)')   &
             '!!!!!!! ATTENTION: tkvh, tkhh, tkvm, tkhm = const = ', &
             tkvhfix,tkhhfix,tkvmfix,tkhmfix,' !!!!!!!!!!!!!!'
      tkvh_b(:,:) = tkvhfix
      tkvm_b(:,:) = tkvmfix
      IF (l3dturb) THEN
         tkhh_b(:,:) = tkhhfix
         tkhm_b(:,:) = tkhmfix
      END IF
    END IF
  ELSEIF (itype_turb.EQ.100) THEN
    yerrmsg  = ' ERROR    *** itype_turb = 100 (const. diff. coeff.) only '//&
               'if lartif_data=.TRUE. *** '
    ierror = 2798
    RETURN
  ENDIF
#endif

END SUBROUTINE turb_organize    

!===============================================================================
!===============================================================================

SUBROUTINE turb_init

!-------------------------------------------------------------------------------
!
! Description:
!  turb_init executes tasks for the turbulence scheme, necessary at the 
!  beginning of the model
!  - determine ntimtke: number of time levels for tke (depending on NL switches)
!  - set ntr_totl:      total number of tracers
!  - set ntr_diff:      number of tracers for vertical diffusion
!  - allocate structures to keep these tracers
!
!-------------------------------------------------------------------------------

INTEGER            :: &
  izdebug,            & ! for debug printout
  iztrcr,             & ! loop index
  izerrstat,          & ! for error checking
  i,j                   ! loop indices

REAL (KIND=wp)     :: &
  rzsethori             ! to initialize l_hori

CHARACTER(LEN=250) :: &
  yerrormsg             ! for error strings

CHARACTER(LEN=250) :: &
  yroutine              ! name of this routine

! End of Header
!-------------------------------------------------------------------------------

  izerrstat     = 0
  yerrormsg(:)  = ' '
  yroutine      = 'turb_init'

  ! Get number of time levels for tke
  ntimtke = UBOUND(tke,4)

  CALL turb_param

  ! Set tkred_sfc to 1.0 (is not used at the moment)
  tkred_sfc(:,:) = 1.0_wp

  ! Get total number of tracers
  ntr_totl = trcr_get_ntrcr()

  ! Allocate structures to keep the tracers for vertical diffusion
  ALLOCATE (itr_diff(ntr_totl))   ! all tracers to be diffused
  ALLOCATE (izbbc   (ntr_totl))   ! all tracers with bottom boundary condition
  ALLOCATE (izturb  (ntr_totl))   ! all tracers for turbulent mixing
  ALLOCATE (ptr     (ntr_totl))   ! is of TYPE modvar

  ! Check tracer list for the ones that should be vertically diffused
  CALL trcr_meta_get(izerrstat, T_TURB_ID, izturb)
  IF (izerrstat /= 0) THEN
    yerrormsg = trcr_errorstr(izerrstat)
    CALL model_abort(my_cart_id, izerrstat, yerrormsg, yroutine)
  ENDIF
  CALL trcr_meta_get(izerrstat, T_BBC_ID, izbbc)
  IF (izerrstat /= 0) THEN
    yerrormsg = trcr_errorstr(izerrstat)
    CALL model_abort(my_cart_id, izerrstat, yerrormsg, yroutine)
  ENDIF

  ntr_diff=0

  trcr_loop: DO iztrcr = 1, ntr_totl
!print *, 'checking tracers:  ', iztrcr, izturb(iztrcr), T_TURB_1D
    ! retrieve the required microphysics tracers (at specified timelevel):
    IF ((iztrcr == idt_qv) .OR. (iztrcr == idt_qc) ) THEN 
      ! qv and qc are not passive in the turbulence model and will be pointed to 
      ! in the turbulence scheme
      CYCLE trcr_loop
    ELSEIF ( izturb(iztrcr) == T_TURB_1D) THEN
      ! This tracer is valid for vertical diffusion
      ntr_diff=ntr_diff+1
      itr_diff(ntr_diff)=iztrcr
!print *, 'getting the next tracer:  ', ntr_diff, itr_diff(ntr_diff)
    END IF

    ptr(iztrcr)%av => NULL()
    ptr(iztrcr)%at => NULL()
    ptr(iztrcr)%sv => NULL()
    ptr(iztrcr)%fc =  .FALSE.

  ENDDO trcr_loop

  ! Initialize the field l_hori
  rzsethori = r_earth/SQRT(eddlat*eddlon)
  DO j = 1, je
    DO i = 1, ie
      ! Note: a location dependent horizontal scale may be used in the future
      ! l_hori(i,j)=z1/SQRT(edadlat*eddlon*acrlat(j,1))

      l_hori(i,j) = rzsethori
     END DO
   END DO

!-------------------------------------------------------------------------------

END SUBROUTINE turb_init

!===============================================================================
!===============================================================================

SUBROUTINE turb_prepare (ntlv)

!-------------------------------------------------------------------------------
!
! Description:
!  turb_prepare precomputes fields, for which neighbourhood relations are 
!  necessary (gradients). The fields are:
!    hdiv:       horizontal divergence
!    hdef2:      horizontal deformation square
!    dwdx, dwdy: horizontal gradients of vertical wind
!
!  These fields are allocated in a special routine: turb_wkarr_prep_alloc
!
! Method:
!
!-------------------------------------------------------------------------------

INTEGER, INTENT(IN) :: ntlv    ! time level used

REAL (KIND=wp) :: edh, zs11, zs12, zs13, zs21, zs22, zs23, zs33

INTEGER        :: i,j,k

! End of Header
!-------------------------------------------------------------------------------

  IF ((itype_turb==3) .AND. ((itype_sher /= 0) .OR. ltkeshs .OR. loutshs)) THEN
    ! 3D sherar correction requested

    ! Horizontal divergence and horizontal deformation square:
    DO k = 2, ke

      DO j = 2, je-1
        DO i = 2, ie-1
          !Model layer slope (times respective horizontal grid spacing):
          lays(i,j,1) = ( hhl(i+1,j,k) - hhl(i-1,j,k) ) * 0.5_wp
          lays(i,j,2) = ( hhl(i,j+1,k) - hhl(i,j-1,k) ) * 0.5_wp

          !Vertical wind gradient:
          edh=2.0_wp/(hhl(i,j,k-1)-hhl(i,j,k+1))
          vgrd(i,j,1) = ( u_m(i,j,k-1) - u_m(i,j,k) ) *edh !d3v1
          vgrd(i,j,2) = ( v_m(i,j,k-1) - v_m(i,j,k) ) *edh !d2v2
        END DO
      END DO

      DO j = 2, je-1
        DO i = 2, ie-1
          !Aligned horizontal wind gradients:
          zs11 = ( u(i,j,k  ,ntlv) - u(i-1,j,k  ,ntlv)   +       &
                   u(i,j,k-1,ntlv) - u(i-1,j,k-1,ntlv) ) * z1d2
          zs11 = ( zs11 - lays(i,j,1)*vgrd(i,j,1) ) * eddlon*acrlat(j,1) !d1v1

          zs22 = ( v(i,j,k  ,ntlv) - v(i,j-1,k  ,ntlv)   +       &
                   v(i,j,k-1,ntlv) - v(i,j-1,k-1,ntlv) ) * z1d2
          zs22 = ( zs22 - lays(i,j,2)*vgrd(i,j,2) ) * edadlat            !d2v2

         !Mixed horizontal wind gradients :
          zs12 = ( v_m(i+1,j,k  ) - v_m(i-1,j,k  ) + &
                   v_m(i+1,j,k-1) - v_m(i-1,j,k-1) ) * z1d4
          zs12 = ( zs12 - lays(i,j,1)*vgrd(i,j,2) ) * eddlon*acrlat(j,1) !d1v2

          zs21 = ( u_m(i,j+1,k  ) - u_m(i,j-1,k  ) + &
                   u_m(i,j+1,k-1) - u_m(i,j-1,k-1) ) * z1d4
          zs21 = ( zs21 - lays(i,j,2)*vgrd(i,j,1) ) * edadlat            !d2v1 

          !Horizontal divergence:
          hdiv(i,j,k) = zs11 + zs22

          !Horizontal deformation square:
          hdef2(i,j,k) = ( zs12 + zs21 )**2 + ( zs11 - zs22 )**2
        END DO
      END DO

      DO i = 1, ie
        hdiv   (i, 1,k) = hdiv  (i,   2,k)
        hdiv   (i,je,k) = hdiv  (i,je-1,k)
        hdef2  (i, 1,k) = hdef2 (i,   2,k)
        hdef2  (i,je,k) = hdef2 (i,je-1,k)
      ENDDO

      DO j = 1, je
        hdiv   (1, j,k) = hdiv  (2,   j,k)
        hdiv   (ie,j,k) = hdiv  (ie-1,j,k)
        hdef2  (1, j,k) = hdef2 (2,   j,k)
        hdef2  (ie,j,k) = hdef2 (ie-1,j,k)
      ENDDO


      IF (itype_sher.EQ.2) THEN !shear correction by vertical wind requested

        ! Horizontal gradients of vertical wind-speed:
        DO j = 2, je-1
          DO i = 1, ie-1
            ! Horizontal differences of vertical wind:
            zs13 = ( w(i+1,j,k,ntlv) - w(i-1,j,k,ntlv) ) * 0.5_wp
            zs23 = ( w(i,j+1,k,ntlv) - w(i,j-1,k,ntlv) ) * 0.5_wp

            ! Aligned vertical wind gradient
            zs33 = ( w(i,j,k-1,ntlv) - w(i,j,k+1,ntlv) ) / ( hhl(i,j,k-1) - hhl(i,j,k+1) )

            ! Horizontal gradients of vertical wind
            dwdx(i,j,k) = ( zs13 - lays(i,j,1)*zs33 ) * acrlat(j,1) * eddlon !d1v3
            dwdy(i,j,k) = ( zs23 - lays(i,j,2)*zs33 ) * edadlat              !d2v3
          END DO
        END DO

        DO i = 1, ie
          dwdx   (i, 1,k) = dwdx  (i,   2,k)
          dwdx   (i,je,k) = dwdx  (i,je-1,k)
          dwdy   (i, 1,k) = dwdy  (i,   2,k)
          dwdy   (i,je,k) = dwdy  (i,je-1,k)
        ENDDO

        DO j = 1, je
          dwdx   (1, j,k) = dwdx  (2,   j,k)
          dwdx   (ie,j,k) = dwdx  (ie-1,j,k)
          dwdy   (1, j,k) = dwdy  (2,   j,k)
          dwdy   (ie,j,k) = dwdy  (ie-1,j,k)
        ENDDO

      END IF

    END DO  ! k-loop from 2, ke-1

  END IF    ! lturatm .AND. (itype_sher.NE.0 .OR. ltkeshs .OR. loutshs

  ! k=1
  ! Calculation of needed effective expressions from horizontal wind gradients
  ! at the surface level and excluding the most outer horizontal boundary lines:

  IF ((itype_tran==2) .AND. itype_sher.EQ.2) THEN
    ! surface-atmospher transfer active and shear correction by vertical wind requested:

    DO j = 2, je-1      !???? by Matthias:  jstart, jend
      DO i = 2, ie-1    !???? by Matthias:  istart, iend
        ! Surface slopes:
        lays(i,j,1) = z1d2 * ( hhl(i+1,j,ke1) - hhl(i-1,j,ke1) ) * eddlon*acrlat(j,1)
        lays(i,j,2) = z1d2 * ( hhl(i,j+1,ke1) - hhl(i,j+1,ke1) ) * edadlat

        ! effective vertical difference of vertical wind speed at lowest full level:
        zs33 = z1d2 * ( w(i,j,ke,ntlv) - w(i,j,ke1,ntlv) ) != z1d2*(wm(ke) + wm(ke+1)) - wm(ke1) 

        ! effective horizontal differences of vertical wind speed:
        dwdx(i,j,ke1) = lays(i,j,1)*zs33
        dwdy(i,j,ke1) = lays(i,j,2)*zs33

        ! effective squared difference of deformation:
        hdef2(i,j,ke1) = ( lays(i,j,2)*u_m(i,j,ke) - lays(i,j,1)*v_m(i,j,ke) )**2
      END DO
    END DO

    ! Note: 
    ! The effective differences correspond to gradients in relation to the preliminary
    !  layer depth of 1m.
    ! In 'turbtran' these effektive differences will be divided by the effective
    !  depth of the surface layer in order to get the desired gradient values.

    DO i = 1, ie
      dwdx   (i, 1,ke1) = dwdx  (i,   2,ke1)
      dwdx   (i,je,ke1) = dwdx  (i,je-1,ke1)
      dwdy   (i, 1,ke1) = dwdy  (i,   2,ke1)
      dwdy   (i,je,ke1) = dwdy  (i,je-1,ke1)
      hdef2  (i, 1,ke1) = hdef2 (i,   2,ke1)
      hdef2  (i,je,ke1) = hdef2 (i,je-1,ke1)
    ENDDO

    DO j = 1, je
      dwdx   (1, j,ke1) = dwdx  (2,   j,ke1)
      dwdx   (ie,j,ke1) = dwdx  (ie-1,j,ke1)
      dwdy   (1, j,ke1) = dwdy  (2,   j,ke1)
      dwdy   (ie,j,ke1) = dwdy  (ie-1,j,ke1)
      hdef2  (1, j,ke1) = hdef2 (2,   j,ke1)
      hdef2  (ie,j,ke1) = hdef2 (ie-1,j,ke1)
    ENDDO

  END IF

!-------------------------------------------------------------------------------

  ! Compute random noise for idealized surface fluxes in ijk structure:
  !  Will be copied to block on h0noise_b(:).
  IF (lartif_data .AND. lsensiflux_fix) THEN
    CALL gen_H0noise (h0noise)
  END IF

!-------------------------------------------------------------------------------

END SUBROUTINE turb_prepare

!===============================================================================
!===============================================================================

SUBROUTINE turb_init_copy

!-------------------------------------------------------------------------------
!
! Description:
!   Register all required copies to/from block for turbulence schemes 
!
!-------------------------------------------------------------------------------

  ! init copy list
  CALL init_copy_list(turCopyList)

  ! Register required copy to block
  ! Variables with in intent IN
  CALL register_copy (l_hori_b      , turCopyList, copyToBlockF)
  CALL register_copy (dp0_b         , turCopyList, copyToBlockF)
  CALL register_copy (hhl_b         , turCopyList, copyToBlockF)
  CALL register_copy (fr_land_b     , turCopyList, copyToBlockF)
  CALL register_copy (depth_lk_b    , turCopyList, copyToBlockF)
  CALL register_copy (gz0_b         , turCopyList, copyToBlockF)
  CALL register_copy (sai_b         , turCopyList, copyToBlockF)
  CALL register_copy (l_pat_b       , turCopyList, copyToBlockF)
  CALL register_copy (t_snow_b      , turCopyList, copyToBlockF)
  CALL register_copy (t_s_b         , turCopyList, copyToBlockF)
  CALL register_copy (t_g_b         , turCopyList, copyToBlockF)
  CALL register_copy (qv_s_b        , turCopyList, copyToBlockF)
  CALL register_copy (ps_b          , turCopyList, copyToBlockF)
  CALL register_copy (u_m_b         , turCopyList, copyToBlockF)
  CALL register_copy (v_m_b         , turCopyList, copyToBlockF)
  CALL register_copy (w_b           , turCopyList, copyToBlockF)
  CALL register_copy (t_b           , turCopyList, copyToBlockF)
  CALL register_copy (p0_b          , turCopyList, copyToBlockF)
  CALL register_copy (pp_b          , turCopyList, copyToBlockF)
  CALL register_copy (qv_b          , turCopyList, copyToBlockF)
  CALL register_copy (qc_b          , turCopyList, copyToBlockF)
  IF (itype_vdif >= 0) THEN
  CALL register_copy (qi_b          , turCopyList, copyToBlockF)
  ENDIF
  CALL register_copy (rho_b         , turCopyList, copyToBlockF)
  CALL register_copy (tcm_b         , turCopyList, copyToBlockF)
  CALL register_copy (tch_b         , turCopyList, copyToBlockF)
  CALL register_copy (tvm_b         , turCopyList, copyToBlockF)
  CALL register_copy (tvh_b         , turCopyList, copyToBlockF)
  CALL register_copy (tfm_b         , turCopyList, copyToBlockF)
  CALL register_copy (tfh_b         , turCopyList, copyToBlockF)
  CALL register_copy (tfv_b         , turCopyList, copyToBlockF)
  CALL register_copy (tkr_b         , turCopyList, copyToBlockF)
  CALL register_copy (tkred_sfc_b   , turCopyList, copyToBlockF)
  CALL register_copy (tkvm_b        , turCopyList, copyToBlockF)
  CALL register_copy (tkvh_b        , turCopyList, copyToBlockF)
  CALL register_copy (rcld_b        , turCopyList, copyToBlockF)
  CALL register_copy (hdef2_b       , turCopyList, copyToBlockF)
  CALL register_copy (hdiv_b        , turCopyList, copyToBlockF)
  CALL register_copy (dwdx_b        , turCopyList, copyToBlockF)
  CALL register_copy (dwdy_b        , turCopyList, copyToBlockF)
  CALL register_copy (tket_conv_b   , turCopyList, copyToBlockF)
  CALL register_copy (utens_b       , turCopyList, copyToBlockF)
  CALL register_copy (vtens_b       , turCopyList, copyToBlockF)
  CALL register_copy (ttens_b       , turCopyList, copyToBlockF)
  CALL register_copy (qvtens_b      , turCopyList, copyToBlockF)
  CALL register_copy (qctens_b      , turCopyList, copyToBlockF)
  CALL register_copy (qitens_b      , turCopyList, copyToBlockF)
  CALL register_copy (tketens_b     , turCopyList, copyToBlockF)
  CALL register_copy (tket_adv_b    , turCopyList, copyToBlockF)
  CALL register_copy (dqvdt_b       , turCopyList, copyToBlockF)
  CALL register_copy (ut_turb_b     , turCopyList, copyToBlockF)
  CALL register_copy (vt_turb_b     , turCopyList, copyToBlockF)
  CALL register_copy (ut_sso_b      , turCopyList, copyToBlockF)
  CALL register_copy (vt_sso_b      , turCopyList, copyToBlockF)
  CALL register_copy (ut_conv_b     , turCopyList, copyToBlockF)
  CALL register_copy (vt_conv_b     , turCopyList, copyToBlockF)
  CALL register_copy (t_2m_b        , turCopyList, copyToBlockF)
  CALL register_copy (qv_2m_b       , turCopyList, copyToBlockF)
  CALL register_copy (td_2m_b       , turCopyList, copyToBlockF)
  CALL register_copy (rh_2m_b       , turCopyList, copyToBlockF)
  CALL register_copy (u_10m_b       , turCopyList, copyToBlockF)
  CALL register_copy (v_10m_b       , turCopyList, copyToBlockF)
  CALL register_copy (shfl_s_b      , turCopyList, copyToBlockF)
  CALL register_copy (lhfl_s_b      , turCopyList, copyToBlockF)
  CALL register_copy (qvsflx_b      , turCopyList, copyToBlockF)

  IF (ALLOCATED(h_ice)) CALL register_copy(h_ice_b    ,turCopyList,copyToBlockF)

  IF (lartif_data .AND. lsensiflux_fix) THEN
    CALL register_copy (h0noise_b     , turCopyList, copyToBlockF)
  ENDIF

!XL_TODO tke not supported  CALL register_copy(tke_b    ,turCopyList,copyToBlockF)


  ! Variables with intent OUT
  CALL register_copy (gz0_b         , turCopyList, copyFromBlockF)
  CALL register_copy (tcm_b         , turCopyList, copyFromBlockF)
  CALL register_copy (tch_b         , turCopyList, copyFromBlockF)
  CALL register_copy (tvm_b         , turCopyList, copyFromBlockF)
  CALL register_copy (tvh_b         , turCopyList, copyFromBlockF)
  CALL register_copy (tfm_b         , turCopyList, copyFromBlockF)
  CALL register_copy (tfh_b         , turCopyList, copyFromBlockF)
  CALL register_copy (tfv_b         , turCopyList, copyFromBlockF)
  CALL register_copy (tkr_b         , turCopyList, copyFromBlockF)
  CALL register_copy (tkvm_b        , turCopyList, copyFromBlockF)
  CALL register_copy (tkvh_b        , turCopyList, copyFromBlockF)
  CALL register_copy (rcld_b        , turCopyList, copyFromBlockF)
  CALL register_copy (tkhm_b        , turCopyList, copyFromBlockF)
  CALL register_copy (tkhh_b        , turCopyList, copyFromBlockF)
  CALL register_copy (edr_b         , turCopyList, copyFromBlockF)
  CALL register_copy (tket_sso_b    , turCopyList, copyFromBlockF)
  CALL register_copy (tket_hshr_b   , turCopyList, copyFromBlockF)
  CALL register_copy (ut_turb_b     , turCopyList, copyFromBlockF)
  CALL register_copy (vt_turb_b     , turCopyList, copyFromBlockF)
  CALL register_copy (ttens_b       , turCopyList, copyFromBlockF)
  CALL register_copy (qvtens_b      , turCopyList, copyFromBlockF)
  CALL register_copy (qctens_b      , turCopyList, copyFromBlockF)
  IF (itype_vdif >= 0) THEN
  CALL register_copy (qitens_b      , turCopyList, copyFromBlockF)
  ENDIF
  CALL register_copy (tketens_b     , turCopyList, copyFromBlockF)
  CALL register_copy (tket_adv_b    , turCopyList, copyFromBlockF)
  CALL register_copy (dqvdt_b       , turCopyList, copyFromBlockF)
  CALL register_copy (t_2m_b        , turCopyList, copyFromBlockF)
  CALL register_copy (qv_2m_b       , turCopyList, copyFromBlockF)
  CALL register_copy (td_2m_b       , turCopyList, copyFromBlockF)
  CALL register_copy (rh_2m_b       , turCopyList, copyFromBlockF)
  CALL register_copy (u_10m_b       , turCopyList, copyFromBlockF)
  CALL register_copy (v_10m_b       , turCopyList, copyFromBlockF)
  CALL register_copy (shfl_s_b      , turCopyList, copyFromBlockF)
  CALL register_copy (lhfl_s_b      , turCopyList, copyFromBlockF)
  CALL register_copy (qvsflx_b      , turCopyList, copyFromBlockF)
  IF (lartif_data .AND. lsensiflux_fix) THEN
  CALL register_copy (t_s_b         , turCopyList, copyFromBlockF)
  CALL register_copy (t_s_new_b     , turCopyList, copyFromBlockF)
  CALL register_copy (t_g_b         , turCopyList, copyFromBlockF)
  CALL register_copy (t_g_new_b     , turCopyList, copyFromBlockF)
  CALL register_copy (qv_s_b        , turCopyList, copyFromBlockF)
  CALL register_copy (qv_s_new_b    , turCopyList, copyFromBlockF)
  ENDIF

!-------------------------------------------------------------------------------

END SUBROUTINE turb_init_copy

!===============================================================================
!===============================================================================

SUBROUTINE turb_prepare_wkarr_alloc

!-------------------------------------------------------------------------------
!
! Description:
!  allocates necessary data for turb_prepare 
!
!US: But because of copyin/copyout registration, we need these variables 
!    allocated at the very beginning
! Method:
!
!-------------------------------------------------------------------------------

  INTEGER :: izstat=0,           &
             izerrstat=0
  CHARACTER(LEN=80) :: yerrormsg

!-------------------------------------------------------------------------------

  yerrormsg(:) = ' '

! ! fields needed in turbulence scheme
! ALLOCATE (hdiv(ie,je,ke1),   hdef2(ie,je,ke1),                           &
!           dwdx(ie,je,ke1),   dwdy (ie,je,ke1),    STAT=izstat)
! IF (izstat /= 0) THEN
!   PRINT *, '*** Error allocating variables for turb_prepare ***'
!   izerrstat = 2801
!   yerrormsg = '*** Error allocating variables for turb_prepare ***'
!   CALL model_abort(my_cart_id, izerrstat, yerrormsg, 'turb_prepare_wkarr_alloc')
! ENDIF

  ! fields needed only for the preparation
  ALLOCATE (lays(ie,je,ke1),   vgrd(ie,je,ke1),     STAT=izstat)
  IF (izstat /= 0) THEN
    PRINT *, '*** Error allocating variables for turb_prepare ***'
    izerrstat = 2801
    yerrormsg = '*** Error allocating variables for turb_prepare ***'
    CALL model_abort(my_cart_id, izerrstat, yerrormsg, 'turb_prepare_wkarr_alloc')
  ENDIF

!-------------------------------------------------------------------------------

END SUBROUTINE turb_prepare_wkarr_alloc

!===============================================================================
!===============================================================================

SUBROUTINE turb_prepare_wkarr_dealloc

!-------------------------------------------------------------------------------
!
! Description:
!  allocates necessary data for turb_prepare 
!
! Method:
!
!-------------------------------------------------------------------------------

  INTEGER :: izstat=0,           &
             izerrstat=0
  CHARACTER(LEN=80) :: yerrormsg

!-------------------------------------------------------------------------------

  yerrormsg(:) = ' '

! DEALLOCATE (hdiv, hdef2, dwdx, dwdy, lays, vgrd,    STAT=izstat)
  DEALLOCATE (lays, vgrd,                             STAT=izstat)
  IF (izstat /= 0) THEN
    PRINT *, '*** Error deallocating variables for turb_prepare ***'
    izerrstat = 2801
    yerrormsg = '*** Error allocating variables for turb_prepare ***'
    CALL model_abort(my_cart_id, izerrstat, yerrormsg, 'turb_prepare_wkarr_dealloc')
  ENDIF

!-------------------------------------------------------------------------------

END SUBROUTINE turb_prepare_wkarr_dealloc

!===============================================================================
!==============================================================================
!+ To clean up at the end of the program
!------------------------------------------------------------------------------    

SUBROUTINE turb_finalize

!------------------------------------------------------------------------------
!
! Description:
!
! Method:
!
!------------------------------------------------------------------------------

  INTEGER :: izerror
  CHARACTER(LEN=80) :: yerrmsg

!----------- End of header ----------------------------------------------------

  izerror = 0

  ! Deallocate the working arrays

#ifdef _OPENACC
  CALL turb_prepare_wkarr_dealloc

  CALL turb_wkarr_dealloc(izerror)
  IF (izerror /= 0) THEN
      yerrmsg = 'Error in turb_wkarr_dealloc'
      CALL model_abort(my_cart_id, 666, yerrmsg, 'turb_finalize')
  ENDIF
#endif

!------------------------------------------------------------------------------
! End of the subroutine 
!------------------------------------------------------------------------------

END SUBROUTINE turb_finalize

!==============================================================================
!===============================================================================

END MODULE turb_interface
