!+ Module for boundary relaxation, Asselin filtering and saturation adjustment
!------------------------------------------------------------------------------

MODULE src_relaxation

!------------------------------------------------------------------------------
!
! Description:
!   The module "relaxation" performs the lateral boundary relaxation,
!   the Asselin filtering and the saturation adjustment at the end of 
!   a large time step.
!
! Current Code Owner: DWD, Michael Baldauf
!  phone:  +49  69  8062 2733
!  fax:    +49  69  8062 3721
!  email:  Michael.Baldauf@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Guenther Doms     
!  Initial release
! 1.2        1998/03/30 Ulrich Schaettler
!  Introduction of new routine init_relaxation.
! 1.4        1998/05/22 Guenther Doms
!  Modifications for optional use of 2-TL integration scheme
! 1.5        1998/06/29 Guenther Doms
!  New formulation of the upper boundary Rayleigh damping layer.
!  See Section 5.3 in the Scientific Documentation.
! 1.7        1998/07/16 Guenther Doms
!  Change in the calling list of routine 'calps'.
! 1.11       1998/10/13 Guenther Doms
!  Code optimization as suggested by U.Schaettler.
! 1.12       1998/10/19 Ulrich Schaettler
!  First saturation adjustment moved to organize_dynamics.
! 1.20       1999/01/07 Guenther Doms
!  Renaming of some global variables
! 1.23       1999/02/08 Guenther Doms
!  Additional Rayleigh damping on qv and qc to avoid
!  inconsitencies with nuging of temperature in the strotosphere.
! 1.24       1999/03/01 Guenther Doms
!  Relaxation of cloud ice included (optionally). Boundary values for cloud ice
!  are interpreted from qc and the specific humidity is recalculated.
! 1.29       1999/05/11 Ulrich Schaettler
!  Eliminated unnecessary dependencies
! 1.34       1999/12/10 Ulrich Schaettler
!  Eliminated timing routines
! 1.39       2000/05/03 Ulrich Schaettler
!  Changed some variable names and calculate variables depended on dt 
!  every time step (needed for the interactive nesting version).
! 2.11       2001/09/28 Ulrich Schaettler
!  Changed treatment of cloud ice, if boundary fields are provided.
! 2.18       2002/07/16 Ulrich Schaettler
!  If working with frames, boundary update is done only for defined points
!  (based on the work by Lucio Torrisi, UGM Rome).
! 3.2        2003/02/07 Ulrich Schaettler
!  If working with frames and Rayleigh damping, check whether enough levels
!  are defined on the full grid (based on the work by Lucio Torrisi, UGM Rome).
! 3.5        2003/09/02 Ulrich Schaettler
!  Eliminated part for computing the surface pressure. This is now done in
!  near_surface to avoid the communication for this variable
! 3.7        2004/02/18 Ulrich Schaettler
!  Introduced lprog_qi for treatment of cloud-ice
! 3.13       2004/12/03 Ulrich Schaettler
!  Modifications to run with latent heat nudging (Klaus Stephan, et.al.)
!  Explicit formulation of lateral boundary relaxation (Jochen Foerstner)
!  Preparations for another option of Rayleigh damping (Lucio Torrisi)
! 3.14       2005/01/25 Lucio Torrisi
!  Introduction of new type of Rayleigh damping: it uses filtered LM fields
!  instead of boundary fields
! 3.15       2005/03/03 Ulrich Schaettler
!  Replaced FLOAT by REAL
! 3.16       2005/07/22 Ulrich Schaettler
!  Correction for lhn: compute time level nold only in 3tl scheme
! 3.18       2006/03/03 Klaus Stephan
!  LHN namelist parameter moved to data_lheat_nudge to avoid to many dependencies
!  Corrected setting of hfw_m_nb (for itype_spubc = 2)
!  Introduced temperature increment due to latent heat (Jochen Foerstner)
! 3.21       2006/12/04 Lucio Torrisi
!  Corrections for using the new Rayleigh Damping option itype_spubc=2
!  New Namelist variable nfi_spubc2 to control number of smoothing applications
! V3_23        2007/03/30 Ulrich Schaettler, Lucio Torrisi, Jochen Foerstner
!  Skip boundary relaxation and Rayleigh damping for the last time step
!  (if last time step is just for a full hour, the next boundary values are missing)
!  Added Rayleigh damping also for qr, qs, qg (for itype_spubc=2 at the moment)
!  Added new treatment of lateral boundaries for qr, qs, qg
!  Introduced new namelist variable itype_lbcqx
! V4_4         2008/07/16 Ulrich Schaettler
!  Eliminated timing variables which are unused
! V4_5         2008/09/10 Guenther Zaengl
!  Add possibility to combine radiative lateral boundaries with weak Davies relaxation
! V4_9         2009/07/16 Ulrich Schaettler, Heike Vogel, Christian Bollmann
!  Introduced treatment of variables for COSMO-ART
!  Implemented 3D versions of operators
!  Renamed switch itype_lbcqx to more meaningful itype_outflow_qrsg
!  Inserted call to SR collapse, to change loop indices
!  Eliminate reduction of qv in case llb_qi = .TRUE.
! V4_11        2009/11/30 Oli Fuhrer
!  Implementation of modified Asselin Filter after Williams 2009
! V4_12        2010/05/11 Ulrich Schaettler, Oli Fuhrer
!  Eliminated statement functions fpvsw, fpvsi, fqvs
!  Eliminated lhdiff_mask; compute hd_mask in any case
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_17        2011/02/24 Ulrich Blahak
!  Adapted interface of exchg_boundaries; corrected kzdims(1:20) -> kzdims(1:24);
!  eliminated my_peri_neigh
! V4_18        2011/05/26 Ulrich Schaettler
!  Introduced conditional compilation for Nudging
!  Separation of time interpolation for COSMO-ART / Pollen fields from the rest
!   (Christoph Knote, for COSMO-ART)
! V4_23        2012/05/10 Ulrich Schaettler
!  Added _ireals to literal constants
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer, Ulrich Blahak, Hans-Juergen Panitz
!  Replaced qx-variables by using them from the tracer module
!  UB: Implemented internal switch "l2mom_satads" to be able
!   to switch on the extra saturation adjustments outside the microphysics parts.
!  Introduced nexch_tag for MPI boundary exchange tag to replace ntstep (HJP)
! V4_26        2012/12/06 Anne Roches
!  Changes and technical adaptations to the tracer handling
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  MESSy interface introduced: tracer relaxation
! V4_28        2013/07/12 KIT, Ulrich Schaettler
!  Changes to adapt COSMO-ART to new tracer module: some dependencies to
!  COSMOART and POLLEN deleted, because this is now handled by the tracer module
!  Use subroutines and variables for vertical grid and reference atmospheres
!    from module vgrid_refatm_utils
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
!  For the COSMO-Model only use vcoord from vgrid_refatm_utils
! V5_1         2014-11-28 Ulrich Schaettler, Ulrich Blahak, Oliver Fuhrer
!                         Michael Baldauf, Anne Roches
!  Compute the Rayleigh damping coefficients only for the layer above kflat
!   and use the reference hhl-profile (US)
!  QNICE in a below "exception" should only be touched in case
!   of running the 2-moment scheme. This has been fixed. (UB)
!  Replaced ireals by wp (working precision) (OF)
!  Removed option lexpl_lbc=.FALSE. (MB)
!  Removed hacks for the tracer module (AR)
!  Bugfix: for QNICE, use idt_qni instead of the tracer index for qc
! V5_2         2015-05-21 Ulrich Blahak
!  Replaced vcoord%kflat by ke in 2 loops for the computation of rdcoef. (UB)
! V5_3         2015-10-09 Oliver Fuhrer
!  Splitted SR sardass in two subroutines relaxation and timefilter
!  Implemented serialization comments (OF)
! V5_3a        2015-11-24 Klaus Stephan
!  Moved update of tt_lheat to src_lheat_nudge
! V5_4b        2016-07-12 Pascal Spoerri, Oliver Fuhrer
!  Restricted saturation adjustment to the inner domain
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

! 1. vertical coordinate parameters and related variables
! -------------------------------------------------------

    hhl_prof,     & ! a special hhl-profile

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------
    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ieje,         & ! ie * je
    ke,           & ! number of grid points in vertical direction
    ke1,          & ! ke+1
    ke_rd,        & ! lowest level with Rayleigh-damping

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

! 4. variables for the time discretization and related variables
! --------------------------------------------------------------
    dt,           & ! long time-step
    dt2,          & ! 2 * dt
    epsass,       & ! eps for the Asselin-filter
    alphaass,     & ! weight for Williams 2009 modification (0.5 < alpha <= 1)

! 8. Organizational variables to handle the COSMO humidity tracers
! ----------------------------------------------------------------
    idt_qv,  idt_qc,  idt_qi

#ifdef TWOMOM_SB
USE data_modelconfig, ONLY :   &
     idt_qni
#endif

! end of data_modelconfig

!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 1. mathematical constants
! -------------------------
    pi,           & ! circle constant

! 2. physical constants and related variables
! -------------------------------------------
    r_d,          & ! gas constant for dry air
    rdv,          & ! r_d / r_v
    o_m_rdv,      & ! 1 - r_d/r_v
    rvd_m_o,      & ! r_v/r_d - 1
    cp_d,         & ! specific heat of dry air at constant pressure
    cpdr,         & ! 1 / cp_d
    lh_v,         & ! latent heat of vapourization
    g,            & ! acceleration due to gravity

! 3. constants for parametrizations
! ---------------------------------
    b1,           & ! variables for computing the saturation vapour pressure
    b2w,          & ! over water (w) and ice (i)
    b2i,          & ! over water (w) and ice (i)
    b3,           & !               -- " --
    b4w,          & !               -- " --
    b4i,          & !               -- " --
    b234w           ! b2w * (b3 - b4w)

! end of data_constants

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------
    rho0       ,    & ! reference density at the full model levels    (kg/m3)
    dp0        ,    & ! reference pressure thickness of layers        ( Pa)
    p0         ,    & ! reference pressure at main levels             ( Pa)

! 2. external parameter fields                                        (unit)
! ----------------------------
    rmy        ,    & ! Davis-parameter for relaxation (mass, qv, qc)   --
    rmyq       ,    & ! Davis-parameter for relaxation (qr, qs, qg)     --
    hd_mask    ,    & ! 3D-domain mask for horizontal diffusion         --
    least_lbdz ,    & ! mask for eastern  lateral boundary zone
    lwest_lbdz ,    & ! mask for western  lateral boundary zone
    lnorth_lbdz,    & ! mask for northern lateral boundary zone
    lsouth_lbdz,    & ! mask for southern lateral boundary zone

! 3. prognostic variables                                             (unit)
! -----------------------
    u          ,    & ! zonal wind speed                              ( m/s )
    v          ,    & ! meridional wind speed                         ( m/s )
    w          ,    & ! vertical wind speed (defined on half levels)  ( m/s )
    t          ,    & ! temperature                                   (  k  )
    pp         ,    & ! deviation from the reference pressure         ( pa  )

! 4. tendency fields for the prognostic variables                     (unit )
! -----------------------------------------------
    utens      ,    & ! u-tendency without sound-wave terms           ( m/s2)
    vtens      ,    & ! v-tendency without sound-wave terms           ( m/s2)

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------

    rdcoef     ,    & ! Rayleigh damping coefficients
    tinc_lh    ,    & ! temperature increment due to latent heat      (  K  )

! 8. fields for the boundary values                                   (unit )
! ---------------------------------
    u_bd          , & ! boundary field for u                          ( m/s )
    v_bd          , & ! boundary field for v                          ( m/s )
    w_bd          , & ! boundary field for w                          ( m/s )
    t_bd          , & ! boundary field for t                          (  k  )
    pp_bd             ! boundary field for pp                         (  pa )

! end of data_fields

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    nstart,       & ! first time step of the forecast
    nstop,        & ! last time step of the forecast
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
    nold,         & ! corresponds to ntstep - 1
    nnow,         & ! corresponds to ntstep
    nnew,         & ! corresponds to ntstep + 1

! 2. boundary definition and update
! ---------------------------------
    nlastbound,   & ! time step of the last boundary update
    nincbound,    & ! time step increment of boundary update
    nbd1,         & ! indices for permutation of the
    nbd2,         & ! two boundary time levels
    lperi_x,      & ! if lartif_data=.TRUE.: periodic boundary conditions in x-dir.
                    ! or with Davies conditions (.FALSE.)
    lperi_y,      & ! if lartif_data=.TRUE.: periodic boundary conditions in y-dir.
                    ! or with Davies conditions (.FALSE.)
    l2dim,        & ! 2 dimensional runs

! 3. controlling the physics
! --------------------------
    lphys,        & ! forecast with physical parametrizations
    lgsp,         & ! forecast with grid-scale precipitation
    itype_gscp,   & ! type of grid-scale precipitaiton physics

! 4. controlling the upper boundary condition
! -------------------------------------------
    lrubc ,       & ! with radiative upper boundary condition
    lspubc ,      & ! with Rayleigh damping in the upper levels
    itype_outflow_qrsg,  & ! type of relaxation treatment for qr, qs, qg
    itype_spubc , & ! type of Rayleigh damping in the upper levels
    nfi_spubc2,   & ! Number of applications of smoother for the determination
                    !  of the large scale field used in the Rayleigh damping
                    !  with itype_spubc=2
    rdheight,     & ! bottom height of Rayleigh damping layer
    nrdtau,       & ! number of time steps in Rayleigh damping time scale

! 5. additional control variables
! -------------------------------
    l2tls,        & ! forecast with 2-TL integration scheme
    l_cosmo_art,  & ! if .TRUE., run the COSMO_ART
    l_pollen,     & ! of pollen
    lcond,        & ! forecast with condensation/evaporation
    ldiabf_satad, & ! include diabatic forcing due to saturation adjustment
    lreproduce,   & ! the results are reproducible in parallel mode
    lw_freeslip,  & ! if .TRUE.: with free slip lateral boundary condition and
                    ! if .FALSE. specified lateral boundary values for w
    lradlbc,      & ! if lartif_data=.TRUE.: radiative lateral boundary conditions (.TRUE.)
                    !                 or with Davies conditions (.FALSE.)
    relax_fac,    & ! reduction factor for strength of lateral boundary relaxation
                    ! (relevant for radiative lateral boundary conditions)

! 7. additional control variables
! -------------------------------
    l2mom_satads, & ! in case of 2-moment scheme, do all the satads
                    ! (like for the 1-moment schemes), not just the
                    ! satad after the microphysics at the end of the timestep.

! 12. controlling verbosity of debug output
! -----------------------------------------
    idbg_level,   & ! to control the verbosity of debug output
    ldebug_dyn,   & ! if .TRUE., debug output for dynamics
    lprintdeb_all   ! .TRUE.:  all tasks print debug output
                    ! .FALSE.: only task 0 prints debug output

! end of data_runcontrol 

!------------------------------------------------------------------------------

USE data_parallel,      ONLY :  &
    num_compute,     & ! number of compute PEs
    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains
    ncomm_type,      & ! type of communication
    my_cart_neigh,   & ! neighbors of this subdomain in the cartesian grid
    icomm_cart,      & ! communicator for the virtual cartesian topology
                       ! that can be used by MPI_WAIT to identify the send
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    nexch_tag,       & ! tag to be used for MPI boundary exchange
                       !  (in calls to exchg_boundaries)
    sendbuf,         & ! sending buffer for boundary exchange:
                       ! 1-4 are used for sending, 5-8 are used for receiving
    isendbuflen,     & ! length of one column of sendbuf
    my_cart_id         ! rank of this subdomain in the cartesian communicator

!------------------------------------------------------------------------------

USE data_io,            ONLY :  &
    lbd_frame,  & ! if .TRUE., boundary data are on a frame
    ilevbotnoframe, & ! bottom model level with b.d. defined on the whole grid
                      ! (model levels below ilevbotnoframe are defined on a frame)
    undef         ! value for "undefined" points

!------------------------------------------------------------------------------

USE data_tracer,        ONLY :                                                &
    T_LBC_ID      ,  T_LBC_ZERO  , T_LBC_CST     , T_LBC_FILE  ,              &
    T_LBC_ZEROGRAD,  T_LBC_USER  ,                                            &
    T_RELAX_ID    ,  T_RELAX_FULL, T_RELAX_INFLOW,                            &
    T_DAMP_ID     ,  T_DAMP_ON   ,                                            &
    T_ADV_ID      ,  T_ADV_ON    ,                                            &
    T_MISSING     ,  T_ERR_NOTFOUND

!------------------------------------------------------------------------------

#ifdef COSMOART
USE data_cosmo_art,      ONLY:  &
    lgas,         & ! of gases
    laero,        & ! of aerosols
    lgasbd,       & ! boundary values of gases
    laerobd,      & ! boundary values of aerosols
    isp_gas,      & ! number of gas phase species
    isp_gastrans, & ! number of transported gas phase species
    isp_aero,     & ! number of aerosol variables
    isp_aerotrans,& ! number of transported aerosol variables
    ! CK 20101204 interpolation in time for boundary fields now separate for ART boundaries
    nlastbound_art, &
    nincbound_art,  &
    nbd1_art,       & ! indices for permutation of the
    nbd2_art,       & ! two boundary time levels
    trcr_idx_gas,   &
    trcr_idx_aero
#endif

!------------------------------------------------------------------------------

USE environment,              ONLY :  comm_barrier, exchg_boundaries,         &
                                      collapse, model_abort
USE meteo_utilities,          ONLY :  satad, calps
USE utilities,                ONLY :  smooth9
USE vgrid_refatm_utils,       ONLY :  vcoord

!------------------------------------------------------------------------------

USE src_tracer,               ONLY :  trcr_get, trcr_get_ntrcr,               &
                                      trcr_meta_get, trcr_errorstr,           &
                                      trcr_get_index, trcr_get_block

!------------------------------------------------------------------------------

#ifdef NUDGING
USE data_lheat_nudge,    ONLY:  &
    llhn,       & ! main switch for LHN
    llhnverif,  & ! main switch for LHN
    tt_lheat      ! profile of t-increments due to latent heating   ( K/s )
                  ! (stored for current and previous timestep)

!------------------------------------------------------------------------------

USE src_lheating,               ONLY :  &
      get_gs_lheating           ! storage of grid scale latent heating for lhn
#endif

#ifdef TWOMOM_SB
USE src_twomom_sb_interface,    ONLY :   &
    set_qni_from_qi_sb_scalar
#endif

!==============================================================================
! Declarations
!==============================================================================

IMPLICIT NONE

!==============================================================================
! Module Procedures
!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure in "relaxation" for initializing the relaxation
!------------------------------------------------------------------------------

SUBROUTINE init_relaxation (ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module initializes values necessary for computing
!   the relaxation, the Asselin-Filter and the saturation adjustment.
!
! Method:
!
!------------------------------------------------------------------------------

! Parameterlist:
!---------------
INTEGER (KIND=iintegers), INTENT(INOUT) ::  &
  ierror

CHARACTER (LEN=*)       , INTENT(INOUT) ::  &
  yerrmsg      ! error message

! Local scalars:
! -------------
  REAL    (KIND=wp   )     ::  &
    zlev
  INTEGER (KIND=iintegers) ::  &
    k, istat, ilevbotrd,   &
    hdm_ubcext, ke_hdm

!----------- End of header ----------------------------------------------------
 
!------------------------------------------------------------------------------
! Begin Subroutine init_relaxation     
!------------------------------------------------------------------------------
  
  ! Allocation of the field rdcoeff (Rayleigh damping coefficient)
  ALLOCATE ( rdcoef(ke)       , STAT=istat )

  ! Parameter for the Rayleigh damping for the simulation of a radiative
  ! upper boundary condition (lspubc = .TRUE.)
  rdcoef(:) = 0.0_wp

  IF (lspubc .EQV. .TRUE.) THEN
    ! Calculation of the Rayleigh damping coefficients as function of
    ! height (see Section 5.4 of the Documentation(Part I)). The reference 
    ! height (refering to z=0) is used, i.e. it is assumed that the
    ! bottom height of the damping layer ('rdheight') is located above the
    ! model layer height which separates terrain-following and flat levels.
    ! Then all hhl-values in these layers are constant and refer to the 
    ! reference height. Only the half level below kflat has to be treated 
    ! different: this is taken from a reference hhl-profile chosen to start
    ! as low as possible (if a sea point is available, it is just the profile
    ! above sea level)

    DO k = 1, ke
      ! all relevant layers have to be above the model layer which separates
      ! terrain-following and flat levels.
!US   zlev = 0.5_wp*(vcoord%vert_coord(k)+vcoord%vert_coord(k+1))
      zlev = 0.5_wp*(hhl_prof(k) + hhl_prof(k+1))
      IF( zlev >= rdheight ) THEN
        rdcoef(k) = 1.0_wp - COS( pi*(zlev-rdheight)/(hhl_prof(1)-rdheight) )
        rdcoef(k) = rdcoef(k) / ( 2.0_wp*nrdtau*dt )
        ke_rd = k
      ENDIF
    ENDDO

    ! The coefficients for Rayleigh damping are premultiplied by the time
    ! step for integration (2.0*dt for Leapfrog, dt for 2-TL scheme)
    ! This is dependent on the timestep and therefore on the grid used.
    ! The multiplication with dt therefore is done in routine relaxation
    ! during the time stepping
    IF (.NOT. l2tls) THEN
      rdcoef(:) = 2.0_wp * rdcoef(:)
    ENDIF

    ! Print the damping coefficients to standard output
    IF ( (my_cart_id == 0) .AND. (idbg_level > 1) ) THEN
      PRINT*, ' '
      PRINT*, ' Damping coefficients in Rayleigh damping layer'
      PRINT*, ' level   height    damping coefficient'
      DO k = 1, ke
        zlev = 0.5_wp*(hhl_prof(k) + hhl_prof(k+1))
        IF( zlev >= rdheight ) THEN
          ilevbotrd = k   ! lowest level with Rayleigh damping
          WRITE(*,'(I6,F12.2,F14.6)') k, zlev, dt * rdcoef(k) 
        ENDIF
      ENDDO 
      PRINT*, ' '
      PRINT*, '    Lowest Level with Rayleigh-damping: ', ke_rd
      PRINT*, ' '
      IF ( lbd_frame .AND. itype_spubc == 1 ) THEN
        IF (ilevbotnoframe < ilevbotrd ) THEN
          PRINT *,' ERROR    *** Not enough levels defined on full grid   ***'
          PRINT *,'          ***          for Rayleigh damping !!         ***'
          PRINT *,'          *** Got only ilevbotnoframe = ', ilevbotnoframe
          PRINT *,'          *** But need ilevbotrd      = ', ilevbotrd
          ierror = 1
          yerrmsg = 'Too few levels defined on full grid for Rayleigh damping'
        ENDIF
      ENDIF
    ENDIF

    ! set up mask for horizontal diffusion
    ! 2nd step: HD at points in levels according to the upper relaxation zone
    ! number of additional vertical levels where HD should be applied
    ! (extension of upper boundary cond. zone)
    hdm_ubcext = 4

    ! set remaining vertical levels of hd_mask
    ke_hdm = MIN( ke_rd+hdm_ubcext, ke )
    DO k = 1, ke_hdm
      hd_mask(:,:,k) = 1.0_wp
    END DO

  ENDIF

!------------------------------------------------------------------------------
! End of module procedure init_relaxation
!------------------------------------------------------------------------------

END SUBROUTINE init_relaxation

!==============================================================================
!+ Module procedure in "relaxation" for saturation adjustment and relaxation
!------------------------------------------------------------------------------

SUBROUTINE relaxation (dt_alter)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine computes the impact of condensation and evaporation of 
!   cloud water on the prognostic variables t, qv and qc at time level nnew.
!   All prognostic variable except w are relaxed towards their boundary
!   values at time level nnew.
!   Optionally, Rayleigh damping in the upper model layer is included.

! Method:
!   The saturation adjustment technique is applied for condensation and  
!   evaporation. The lateral boundary relaxation scheme follows Davies.
!   Rayleigh damping is according to Clark.
!
!------------------------------------------------------------------------------

! Parameterlist:
! --------------
  REAL (KIND=wp),     INTENT(IN), OPTIONAL         ::         &
    dt_alter      ! needed for the digital filtering

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    i, j, k,           & !
#ifdef COSMOART
    isp,               & !  
#endif
    iztrcr,            & !
    izcount,           & !
    izr,               & !
#ifdef TWOMOM_SB
    izqc, izqi,        & !
#endif
    izrelaxqi,         & !
    izlbcqi,           & !
    iztrcr_damp,       & !
    izerror,           & !  error status
    kitbnd,            & !  Number of iterations in the saturation adjustment
                         !  after boundary updating at time level nnew
!   kitasf,            & !  Number of iterations in the saturation adjustment
!                        !  after Asselin filtering at time level nnow
    hfwidth, istata, n,& !
    num_var,           & !
    izdebug              !  for additional debug output

  REAL    (KIND=wp   )     ::  &
    zdtird,            & !  Coefficient for linear interpolation of boundary
                         !  values
    zemzdr,            & !  Coefficient for linear interpolation of boundary
                         !  values
!   zdd,               & ! Asselin filter displacement of midpoint
    zwbe, zdtrdcoef, zdt, zrelfac_u, zrelfac_p, zrelfac_s
   
! Local (automatic) arrays:
! -------------------------
  REAL    (KIND=wp   )     ::  &
    zmys    (ie,je),   & !  Relaxation coefficient at the scalar gridpoint
    zmyp    (ie,je),   & !  Extra coefficient for pressure
    zmyq    (ie,je,ke),& !  Relaxation coefficient for qr, qs and qg
    zmyu    (ie,je),   & !  Relaxation coefficient at the u-point
    zmyv    (ie,je),   & !  Relaxation coefficient at the v-point
    ztaur   (ke   ),   & !  Coefficient for Rayleigh damping
    zubd    (ie,je,ke),& !  boundary value of u
    zvbd    (ie,je,ke),& !  boundary value of v
    ztbd    (ie,je,ke),& !  boundary value of t
    zppbd   (ie,je,ke),& !  boundary value of pp
    zu_west (je   ),   & !  value of u at western  boundary
    zu_east (je   ),   & !  value of u at eastern  boundary
    zv_south(ie   ),   & !  value of v at southern boundary
    zv_north(ie   ),   & !  value of v at northern boundary
    zphfe   (ie,je),   & !
!   ztstart (ie,je),   & !
    zr1     (ie,je),   & !
    zr2     (ie,je),   & !
    zr3     (ie,je),   & !
    zr4     (ie,je),   & !
    zr5     (ie,je),   & !
    zr6     (ie,je),   & !
    zr7     (ie,je),   & !
    zr8     (ie,je)      !

  INTEGER  (KIND=iintegers) :: &
    izlbc  (trcr_get_ntrcr()), & ! type of lateral BC for tracers
    izrelax(trcr_get_ntrcr()), & ! type of relaxation for tracers
    izdamp (trcr_get_ntrcr()), & ! type of Rayleigh damping for tracers
    izadv  (trcr_get_ntrcr())    ! type of advection for tracers
!   izsp_adv_lf(trcr_get_ntrcr())


! Local (allocatable) arrays:
! ---------------------------
  REAL (KIND=wp),     ALLOCATABLE ::  &
    field_tmp  (:,:,:),   & !
#ifdef COSMOART
    zcgasbd    (:,:,:,:), & !  boundary value of cgas
    zcaerobd   (:,:,:,:), & !  boundary value of caero
#endif
    ztrbd      (:,:,:,:)    !  interpolated boundary field for tracers

  CHARACTER (LEN=255)      ::  &
    yzerrmsg

  CHARACTER (LEN=25)       :: yzroutine

  INTEGER (KIND=iintegers) ::  &
    kzdims  (24)          !

! Tracer pointers:
!-----------------
  REAL (KIND=wp),     POINTER :: &
    ztrcr_new  (:,:,:)  => NULL(),   & ! tracer variable at nnew
!   ztrcr_now  (:,:,:)  => NULL(),   & ! tracer variable at nnow
!   ztrcr_old  (:,:,:)  => NULL(),   & ! tracer variable at nold
    ztrcr_bd   (:,:,:,:)=> NULL(),   & ! tracer boundary variable
    qv_new     (:,:,:)  => NULL(),   & ! QV at tlev=nnew
!   qv_now     (:,:,:)  => NULL(),   & ! QV at tlev=nnow
    qc_new     (:,:,:)  => NULL()      ! QC at tlev=nnew
!   qc_now     (:,:,:)  => NULL()      ! QC at tlev=nnow

#ifdef COSMOART
! CK 20101204 ART bd update frequency does not need to be the same as meteo.
! therefore, the weights calculated for meteo might be wrong.

! Variables for COSMO-ART:
! ---------------------------
  REAL    (KIND=wp   )     ::  &
    zdtird_art,        & !  Coefficient for linear interpolation of boundary
                         !  values wrt ART boundary update frequency
    zemzdr_art           !  Coefficient for linear interpolation of boundary
                         !  values wrt ART boundary update frequency

  REAL (KIND=wp),     POINTER :: &
    cgas_new      (:,:,:,:)   => NULL(),   & ! cgas at tlev=nnew
    cgas_now      (:,:,:,:)   => NULL(),   & ! cgas at tlev=nnow
    cgas_old      (:,:,:,:)   => NULL(),   & ! cgas at tlev=nold
    cgas_bd       (:,:,:,:,:) => NULL(),   & ! cgas at boundary
    caero_new     (:,:,:,:)   => NULL(),   & ! caero at tlev=nnew
    caero_now     (:,:,:,:)   => NULL(),   & ! caero at tlev=nnow
    caero_old     (:,:,:,:)   => NULL(),   & ! caero at tlev=nold
    caero_bd      (:,:,:,:,:) => NULL()      ! caero at boundary
#endif

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine relaxation
!------------------------------------------------------------------------------
  
izerror = 0_iintegers
yzroutine = 'relaxation'

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

#ifdef COSMOART
IF (l_cosmo_art) THEN
  IF (lgas) THEN
    ALLOCATE (zcgasbd(ie,je,ke,isp_gas),     STAT = izerror)
  ENDIF
  IF (laero) THEN
    ALLOCATE (zcaerobd(ie,je,ke,isp_aero),   STAT = izerror)
  ENDIF
ENDIF
#endif

!------------------------------------------------------------------------------
!  Section 1: Initial setting of local variables
!------------------------------------------------------------------------------
 
  IF (izdebug > 10) THEN
    PRINT *, '        START relaxation'
  ENDIF

! Calculation of the coefficients for the linear interpolation of
! boundary values:

  zdtird = REAL(ntstep + 1 - nlastbound, wp) / REAL(nincbound, wp)
  ! avoid weights > 1 for the last timestep    ! does this change the results of last time step???
  zdtird = MIN( 1.0_wp, zdtird )
  zemzdr = 1.0_wp - zdtird
 
#ifdef COSMOART
! CK 20101204 ART boundary conditions might have different weights, as update
! frequency of meteorology is not necessarily equal to chemistry update freq.
  IF (l_cosmo_art) THEN
    IF (lgasbd .OR. laerobd) THEN
      zdtird_art = REAL(ntstep + 1 - nlastbound_art, wp) / REAL(nincbound_art, wp)
      zemzdr_art = 1.0_wp - zdtird_art

      IF (izdebug > 10) THEN
        WRITE(*,"(A,F4.2,A,I2,A,F4.2,A,I2)") "         ART boundary weighting: ",  &
                                      zemzdr_art," bd_tlev ",nbd1_art,", ",        &
                                      zdtird_art," bd_tlev ",nbd2_art
      ENDIF
    ENDIF
  ENDIF
#endif

! Number of interpolations in saturation adjustment   

  kitbnd = 1
 
! Reduction coefficients for lateral boundary relaxation if lateral radiation condition is used

  IF (lradlbc) THEN
    zrelfac_u = MIN(1.0_wp,2.0_wp*relax_fac)
    zrelfac_s = MIN(1.0_wp,relax_fac)
    zrelfac_p = 0.0_wp
  ELSE
    zrelfac_u = 1.0_wp
    zrelfac_s = 1.0_wp
    zrelfac_p = 1.0_wp
  ENDIF

! Coefficients for implicit calculation of lateral boundary relaxation scheme

  DO j = jstartpar, jendpar
    DO i = istartpar, iendpar
      zmys(i,j) = zrelfac_s*rmy(i,j,1)
      zmyp(i,j) = zrelfac_p*rmy(i,j,1)
      zmyu(i,j) = zrelfac_u*rmy(i,j,2)
      zmyv(i,j) = zrelfac_u*rmy(i,j,3)
    ENDDO
  ENDDO

  IF (ntstep == 0) THEN
    zdt = 2.0_wp * dt
  ELSE
    zdt = dt
  ENDIF
  IF (PRESENT(dt_alter)) THEN
    zdt = dt_alter
  ENDIF
  DO k = 1,ke
    zdtrdcoef = zdt * rdcoef(k)
    ztaur(k)  = zdtrdcoef / (1.0_wp + zdtrdcoef)
  ENDDO

  ! retrieve the required metadata
  CALL trcr_meta_get(izerror, T_LBC_ID, izlbc)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_meta_get(izerror, T_RELAX_ID, izrelax)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_meta_get(izerror, T_DAMP_ID, izdamp)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_meta_get(izerror, T_ADV_ID, izadv)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

!------------------------------------------------------------------------------
!  Section 2:  Lateral boundary relaxation and Rayleigh damping for mass
!              variables
!------------------------------------------------------------------------------

IF (izdebug > 10) THEN
  PRINT *, '        Lateral boundary relaxation with itype_outflow_qrsg = ', &
                                                           itype_outflow_qrsg
  PRINT *, '        Rayleigh damping            with itype_spubc = ',        &
                                                           itype_spubc
ENDIF

! Subsection 2.1: lateral boundary relaxation (u, v, w, t, pp)
!-------------------------------------------------------------

! Do boundary relaxation and Rayleigh damping not for the last time step
! If no more boundary values are present, this could lead to bad values
! US in MCH version, they do also the last step and limit relaxation weight to 1
IF (ntstep < nstop) THEN

  IF ( .NOT. lw_freeslip ) THEN
    IF (lbd_frame) THEN
      DO k = 1, ke
        DO j = jstartpar, jendpar
          DO i = istartpar, iendpar
            IF (w_bd(i,j,k+1,nbd2) /= undef) THEN
              zwbe = zemzdr*w_bd(i,j,k+1,nbd1) + zdtird*w_bd(i,j,k+1,nbd2)
              w(i,j,k+1,nnew) = w(i,j,k+1,nnew)                             &
                                  -zmys(i,j) * (w(i,j,k+1,nnew)-zwbe)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ELSE
      DO k = 1, ke
        DO j = jstartpar, jendpar
          DO i = istartpar, iendpar
            zwbe      = zemzdr*w_bd(i,j,k+1,nbd1) + zdtird*w_bd(i,j,k+1,nbd2)
            w(i,j,k+1,nnew) = w(i,j,k+1,nnew)                               &
                                -zmys(i,j) * (w(i,j,k+1,nnew)-zwbe)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  IF (lbd_frame) THEN
    DO k = 1, ke
      DO j = jstartpar, jendpar
        DO i = istartpar, iendpar
          IF (u_bd(i,j,k,nbd2) /= undef) THEN
            zubd(i,j,k)    = zemzdr*u_bd(i,j,k,nbd1) + zdtird*u_bd(i,j,k,nbd2)
            u(i,j,k,nnew)  = u(i,j,k,nnew)                                  &
                               - zmyu(i,j) * (u(i,j,k,nnew) - zubd(i,j,k))
            zvbd(i,j,k)    = zemzdr*v_bd(i,j,k,nbd1) + zdtird*v_bd(i,j,k,nbd2)
            v(i,j,k,nnew)  = v(i,j,k,nnew)                                  &
                               - zmyv(i,j) * (v(i,j,k,nnew) - zvbd(i,j,k))
          ELSE
            zubd(i,j,k)    = u(i,j,k,nnew)
            zvbd(i,j,k)    = v(i,j,k,nnew)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ELSE
    DO k = 1, ke
      DO j = jstartpar, jendpar
        DO i = istartpar, iendpar
          zubd(i,j,k)    = zemzdr*u_bd(i,j,k,nbd1) + zdtird*u_bd(i,j,k,nbd2)
          u(i,j,k,nnew)  = u(i,j,k,nnew)                                    &
                             - zmyu(i,j) * (u(i,j,k,nnew) - zubd(i,j,k))
          zvbd(i,j,k)    = zemzdr*v_bd(i,j,k,nbd1) + zdtird*v_bd(i,j,k,nbd2)
          v(i,j,k,nnew)  = v(i,j,k,nnew)                                    &
                             - zmyv(i,j) * (v(i,j,k,nnew) - zvbd(i,j,k))
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF (lbd_frame) THEN
    DO k = 1, ke
      DO j = jstartpar, jendpar
        DO i = istartpar, iendpar
          IF (t_bd(i,j,k,nbd2) /= undef) THEN
            zppbd(i,j,k) = zemzdr*pp_bd(i,j,k,nbd1) + zdtird*pp_bd(i,j,k,nbd2)
            ztbd(i,j,k)  = zemzdr* t_bd(i,j,k,nbd1) + zdtird* t_bd(i,j,k,nbd2)
          ELSE
            zppbd(i,j,k) = pp(i,j,k,nnew)
            ztbd(i,j,k)  = t(i,j,k,nnew)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ELSE
    DO k = 1, ke
      DO j = jstartpar, jendpar
        DO i = istartpar, iendpar
          zppbd(i,j,k)   = zemzdr*pp_bd(i,j,k,nbd1) + zdtird*pp_bd(i,j,k,nbd2)
          ztbd(i,j,k)    = zemzdr* t_bd(i,j,k,nbd1) + zdtird* t_bd(i,j,k,nbd2)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  DO k = 1, ke
    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar
        pp(i,j,k,nnew) = pp(i,j,k,nnew)                                     &
                            - zmyp(i,j) * (pp(i,j,k,nnew)-zppbd(i,j,k))
        t (i,j,k,nnew) = t (i,j,k,nnew)                                     &
                            - zmys(i,j) * ( t(i,j,k,nnew)- ztbd(i,j,k))
      ENDDO
    ENDDO
  ENDDO

! Subsection 2.2: Rayleigh damping (itype_spubc=1) (u, v, w, t, pp)
!------------------------------------------------------------------

  IF (lspubc .AND. (itype_spubc == 1)) THEN
    ! for k > ke_rd it is multiplication with 0.0
    DO k = 1, ke_rd
      DO j = jstartpar, jendpar
        DO i = istartpar, iendpar
          u (i,j,k,nnew)   = u (i,j,k,nnew)                                 &
                                - ztaur(k ) * ( u(i,j,k,nnew)- zubd(i,j,k))
          v (i,j,k,nnew)   = v (i,j,k,nnew)                                 &
                                - ztaur(k ) * ( v(i,j,k,nnew)- zvbd(i,j,k))
          w (i,j,k+1,nnew) = w (i,j,k+1,nnew)                               &
                                - ztaur(k ) *   w(i,j,k+1,nnew)
          pp(i,j,k,nnew)   = pp(i,j,k,nnew)                                 &
                                - ztaur(k ) * (pp(i,j,k,nnew)-zppbd(i,j,k))
          t (i,j,k,nnew)   = t (i,j,k,nnew)                                 &
                                - ztaur(k ) * (t (i,j,k,nnew)-ztbd (i,j,k))
        ENDDO
      ENDDO
    ENDDO
  ENDIF

!------------------------------------------------------------------------------
!  Section 3: Lateral boundary relaxation and Rayleigh damping for tracers
!------------------------------------------------------------------------------

  ! Subsection 3.0: Preparations
  !------------------------------

  ! allocation of required fields
  ALLOCATE ( ztrbd(ie, je, ke, trcr_get_ntrcr() ), STAT=izerror )
  ztrbd = 0.0_wp

  ! determine if the special coefficients for relaxation 
  ! at inflow only are needed and compute them if required
  izcount = 0_iintegers
  DO izr = 1, trcr_get_ntrcr()
    IF (izrelax (izr) == T_RELAX_INFLOW) THEN
      izcount = izcount + 1
    ENDIF
  ENDDO 

  IF ( izcount > 0_iintegers ) THEN
    ! change zmyq according to inflow / outflow boundaries
    ! no relaxation of qr, qs and qg is done at outflow boundary points
    DO k = 1, ke
      zu_west(1:je)  = 0.0_wp
      zv_south(1:ie) = 0.0_wp
      zu_east(1:je)  = 0.0_wp
      zv_north(1:ie) = 0.0_wp

      ! western boundary
      IF (my_cart_neigh(1) == -1) THEN
        DO j = jstartpar, jendpar
          zu_west(j)  = u(istart,j,k,nnew)
        ENDDO
      ENDIF
      ! southern boundary
      IF (my_cart_neigh(4) == -1) THEN
        DO i = istartpar, iendpar
          zv_south(i) = v(i,jstart,k,nnew)
        ENDDO
      ENDIF
      ! eastern boundary
      IF (my_cart_neigh(3) == -1) THEN
        DO j = jstartpar, jendpar
          zu_east(j)  = u(iend,  j,k,nnew)
        ENDDO
      ENDIF
      ! northern boundary
      IF (my_cart_neigh(2) == -1) THEN
        DO i = istartpar, iendpar
          zv_north(i) = v(i,jend,  k,nnew)
        ENDDO
      ENDIF

      ! set up zmyq
      DO j = jstartpar, jendpar
        DO i = istartpar, iendpar
          ! point outside of lateral boundary zone
          IF ( .NOT.( lwest_lbdz(i,j) .OR. lsouth_lbdz(i,j) .OR.  &
                      least_lbdz(i,j) .OR. lnorth_lbdz(i,j) ) ) THEN
            zmyq(i,j,k) = 0.0_wp
          ! southwest corner
          ELSE IF ( lsouth_lbdz(i,j) .AND. lwest_lbdz(i,j) .AND.  &
            zv_south(i) < 0.0_wp .AND. zu_west(j) < 0.0_wp ) THEN
            zmyq(i,j,k) = 0.0_wp
          ! southeast corner
          ELSE IF ( lsouth_lbdz(i,j) .AND. least_lbdz(i,j) .AND.  &
            zv_south(i) < 0.0_wp .AND. zu_east(j) > 0.0_wp ) THEN
            zmyq(i,j,k) = 0.0_wp
          ! northwest corner
          ELSE IF ( lnorth_lbdz(i,j) .AND. lwest_lbdz(i,j) .AND.  &
            zv_north(i) > 0.0_wp .AND. zu_west(j) < 0.0_wp ) THEN
            zmyq(i,j,k) = 0.0_wp
          ! northeast corner
          ELSE IF ( lnorth_lbdz(i,j) .AND. least_lbdz(i,j) .AND.  &
            zv_north(i) > 0.0_wp .AND. zu_east(j) > 0.0_wp ) THEN
            zmyq(i,j,k) = 0.0_wp
          ! western boundary
          ELSE IF ( lwest_lbdz(i,j)  .AND. zu_west(j)  < 0.0_wp ) THEN
            zmyq(i,j,k) = 0.0_wp
          ! southern boundary
          ELSE IF ( lsouth_lbdz(i,j) .AND. zv_south(i) < 0.0_wp ) THEN
            zmyq(i,j,k) = 0.0_wp
          ! eastern boundary
          ELSE IF ( least_lbdz(i,j)  .AND. zu_east(j)  > 0.0_wp ) THEN
            zmyq(i,j,k) = 0.0_wp
          ! northern boundary
          ELSE IF ( lnorth_lbdz(i,j) .AND. zv_north(i) > 0.0_wp ) THEN
            zmyq(i,j,k) = 0.0_wp
          ELSE
            ! otherwise set zmyq to value of zmys for variables at scalar points,
            ! i.e. do relaxation
            zmyq(i,j,k) = rmyq(i,j)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! Subsection 3.1: Definition of the boundary field
  !-------------------------------------------------
 

! fuo: COSMO-ART has two options:
!      1) if lgasbd is .true., then they do exactly what we do below
!      2) if lgasbd is .false., they relax against a constant field in ztrcr_bd(i,j,k,1)
!      At one point we should probably implement also 2)

  ! loop over tracers
  DO iztrcr = 1, trcr_get_ntrcr()

    ! we need to compute the time interolated value against which a
    ! field is relaxed if either lateral relaxation is switched on or
    ! if Rayleigh damping in the vertical is applied
    IF ( ( izrelax(iztrcr) == T_RELAX_FULL .OR.                               &
           izrelax(iztrcr) == T_RELAX_INFLOW ) .OR.                               &
         ( izdamp (iztrcr) == T_DAMP_ON      .AND.                               &
           lspubc    .AND.  itype_spubc == 1 ) ) THEN
        
      CALL trcr_get( izerror, iztrcr, ptr_tlev=nnew,                        &
                     ptr=ztrcr_new, ptr_bd=ztrcr_bd )
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr( izerror )
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF

      IF ( izlbc(iztrcr) == T_LBC_ZERO ) THEN
        ztrbd(:,:,:,iztrcr) = 0.0_wp
      ELSEIF ( ANY(izlbc(iztrcr) == (/T_LBC_FILE, T_LBC_USER/) ) ) THEN
        ! loop over levels and latitudes
        DO k = 1, ke
          DO j = jstartpar, jendpar
            ! lin. combination of the old and new BC
            IF ( lbd_frame ) THEN
              DO i = istartpar, iendpar
                IF ( ztrcr_bd(i,j,k,nbd2) /= undef ) THEN
                  ztrbd(i,j,k,iztrcr) = zemzdr * ztrcr_bd(i,j,k,nbd1) +     &
                                        zdtird * ztrcr_bd(i,j,k,nbd2)
                ENDIF
              ENDDO
            ELSE
              DO i = istartpar, iendpar
                ztrbd(i,j,k,iztrcr) = zemzdr * ztrcr_bd(i,j,k,nbd1) +       &
                                      zdtird * ztrcr_bd(i,j,k,nbd2)
              ENDDO
            ENDIF
          ENDDO ! j
        ENDDO ! k
      ELSEIF ( ANY(izlbc(iztrcr) == (/ T_LBC_CST, T_LBC_ZEROGRAD /) ) ) THEN
        !nothing to do
      ELSE
        yzerrmsg = "Boundary relaxation illegal for this type of lateral BC"
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF
    ENDIF
  ENDDO ! iztrcr


  ! Subsection 3.2: Boundary relaxation
  !-----------------------------------

! ROA REMARK: this is not a "reproducible feature" in the sense
!             that it is not handled by an option in the tracer
!             structure.

  ! Treat the QI exception
  IF ( idt_qi > 0 ) THEN
    CALL trcr_meta_get(izerror, idt_qi, T_LBC_ID, izlbcqi)
    IF (izerror /= 0) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
    IF ( izlbcqi /= T_LBC_FILE ) THEN
      CALL trcr_meta_get(izerror, idt_qi, T_RELAX_ID, izrelaxqi)
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF
      IF ( izrelaxqi == T_RELAX_FULL ) THEN
        DO k = 1, ke
          DO j = jstartpar, jendpar
            DO i = istartpar, iendpar
              IF (t_bd(i,j,k,nbd2) /= undef) THEN
                ztrbd(i,j,k,idt_qi) = 0.0_wp
                IF ( ztbd(i,j,k) < 248.15_wp ) THEN
                  ztrbd(i,j,k,idt_qi) = ztrbd(i,j,k,idt_qc)
                  ztrbd(i,j,k,idt_qc) = 0.0_wp
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO

#ifdef TWOMOM_SB
        ! In this case, boundary values for cloud ice content qi have been
        ! interpreted from cloud water content above. Set boundary values
        ! for qni in a way that values are consistent with a fixed mean volume particle
        ! diameter given in the below SR set_qni_from_qi_sb_scalar().
        IF (itype_gscp >= 100) THEN
          DO k = 1, ke
            DO j = jstartpar, jendpar
              DO i = istartpar, iendpar
                ztrbd(i,j,k,idt_qni) = set_qni_from_qi_sb_scalar(ztrbd(i,j,k,idt_qi))
              ENDDO
            ENDDO
          ENDDO
        ENDIF
#endif
      ENDIF
    ENDIF
  ENDIF

  ! loop over tracers
  DO iztrcr = 1, trcr_get_ntrcr() 
    ! get tracer
    IF (  izrelax(iztrcr) == T_RELAX_FULL     .OR.                            &
          izrelax(iztrcr) == T_RELAX_INFLOW  ) THEN
      CALL trcr_get( izerror, iztrcr, ptr_tlev=nnew,                          &
                     ptr=ztrcr_new, ptr_bd=ztrcr_bd )
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF

      IF ( ANY(izlbc(iztrcr) == (/ T_LBC_ZERO, T_LBC_FILE, T_LBC_USER /) ) ) THEN
        ! loop over levels and latitudes
        IF ( izrelax(iztrcr) == T_RELAX_FULL ) THEN
          DO k = 1, ke
            DO j = jstartpar, jendpar
              ! boundary relaxation
              DO i = istartpar, iendpar
                IF ( ztrcr_bd(i,j,k,nbd2) /= undef ) THEN
                  ztrcr_new(i,j,k) = ztrcr_new(i,j,k) - zmys(i,j) *           &
                                    (ztrcr_new(i,j,k) -ztrbd(i,j,k,iztrcr))
                ENDIF
              ENDDO
            ENDDO ! j
          ENDDO ! k
        ELSEIF ( izrelax(iztrcr) == T_RELAX_INFLOW ) THEN
          DO k = 1, ke
            DO j = jstartpar, jendpar
              ! boundary relaxation
              DO i = istartpar, iendpar
                IF ( ztrcr_bd(i,j,k,nbd2) /= undef ) THEN
                  ztrcr_new(i,j,k) = ztrcr_new(i,j,k) - zmyq(i,j,k) *         &
                                    (ztrcr_new(i,j,k) -ztrbd(i,j,k,iztrcr))
                ENDIF
              ENDDO
            ENDDO ! j
          ENDDO ! k
        ENDIF
      ELSE
        ! nothing to do for zerogradient or constant
      ENDIF
    ENDIF
  ENDDO ! iztrcr

  ! Subsection 3.3: Rayleigh damping
  !---------------------------------

  ! loop over tracers
  DO iztrcr = 1, trcr_get_ntrcr() 

  ! check if damping is on for this tracer
    IF ( izdamp(iztrcr) == T_DAMP_ON  .AND. lspubc .AND. itype_spubc ==1 ) THEN

      ! get pointer to tracer (at nnew)
      CALL trcr_get( izerror, iztrcr, ptr_tlev=nnew,                          &
                     ptr=ztrcr_new, ptr_bd=ztrcr_bd )
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF

      ! loop over Rayleigh damping layers
      DO k = 1, ke_rd
        DO j = jstartpar, jendpar
          DO i = istartpar, iendpar
            ! Rayleigh damping
            ztrcr_new(i,j,k) = ztrcr_new(i,j,k)                             &
                           - ztaur(k) * (ztrcr_new(i,j,k)-ztrbd(i,j,k,iztrcr))
          ENDDO
        ENDDO
      ENDDO
    ENDIF

  ENDDO ! loop over tracers


  ! Subsection 3.4: Cleaning
  !---------------------------------

  DEALLOCATE (   ztrbd   )

!------------------------------------------------------------------------------
!  Section 4: Lateral boundary relaxation and Rayleigh damping for COSMO-ART
!------------------------------------------------------------------------------

! Subsection 4.1 : lateral boundary relaxation (ART)
!----------------------------------------------------

#ifdef COSMOART
  IF (l_cosmo_art) THEN
    IF (lgas) THEN

      CALL trcr_get_block(izerror, idx_start=trcr_idx_gas(1), idx_end=trcr_idx_gas(isp_gas), &
                 ptr_tlev = nnew, ptr = cgas_new)
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF
      CALL trcr_get_block(izerror, idx_start=trcr_idx_gas(1), idx_end=trcr_idx_gas(isp_gas), &
                 ptr_tlev = nnow, ptr = cgas_now, ptr_bd = cgas_bd)
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF

      IF (lgasbd) THEN
        DO isp = 1, isp_gastrans
          DO k = 1, ke
            DO j = jstartpar, jendpar
              DO i = istartpar, iendpar
                zcgasbd(i,j,k,isp) = zemzdr_art*cgas_bd(i,j,k,isp,nbd1_art) +   &
                                     zdtird_art*cgas_bd(i,j,k,isp,nbd2_art)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO isp = 1, isp_gastrans
          DO k = 1, ke
            DO j = jstartpar, jendpar
              DO i = istartpar, iendpar
                zcgasbd(i,j,k,isp) = cgas_bd(i,j,k,isp,1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      DO isp = 1, isp_gastrans
        DO k = 1, ke
          DO j = jstartpar, jendpar
            DO i = istartpar, iendpar
              cgas_new(i,j,k,isp) = cgas_new(i,j,k,isp)             &
                - zmys(i,j)*(cgas_new(i,j,k,isp)-zcgasbd(i,j,k,isp))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF (laero) then

      CALL trcr_get_block(izerror, idx_start=trcr_idx_aero(1), idx_end=trcr_idx_aero(isp_aero), &
                 ptr_tlev = nnew, ptr = caero_new)
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF
      CALL trcr_get_block(izerror, idx_start=trcr_idx_aero(1), idx_end=trcr_idx_aero(isp_aero), &
                 ptr_tlev = nnow, ptr = caero_now, ptr_bd = caero_bd)
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF

      IF (laerobd) THEN
        DO isp = 1, isp_aerotrans
          DO k = 1, ke
            DO j = jstartpar, jendpar
              DO i = istartpar, iendpar
                zcaerobd(i,j,k,isp) = zemzdr_art*caero_bd(i,j,k,isp,nbd1_art) +     &
                                      zdtird_art*caero_bd(i,j,k,isp,nbd2_art)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO isp = 1, isp_aerotrans
          DO k = 1, ke
            DO j = jstartpar, jendpar
              DO i = istartpar, iendpar
                zcaerobd(i,j,k,isp) = caero_bd(i,j,k,isp,1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      DO isp = 1, isp_aerotrans
        DO k = 1, ke
          DO j = jstartpar, jendpar
            DO i = istartpar, iendpar
              caero_new(i,j,k,isp) = caero_new(i,j,k,isp)           &
                - zmys(i,j)*(caero_new(i,j,k,isp)-zcaerobd(i,j,k,isp))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDIF
#endif

! Subsection 4.2 : Rayleigh damping (ART)
!----------------------------------------

#ifdef COSMOART
  IF (l_cosmo_art) THEN
  
    IF (lgas) THEN
      DO isp = 1,isp_gastrans
        DO k = 1, ke
          DO j = jstartpar, jendpar
            DO i = istartpar, iendpar
              cgas_new(i,j,k,isp) = cgas_new(i,j,k,isp)         &
                             - ztaur(k)*(cgas_new(i,j,k,isp)-zcgasbd(i,j,k,isp))
            ENDDO
          ENDDO 
        ENDDO
      ENDDO
    ENDIF
    
    IF (laero) THEN
      DO isp = 1,isp_aerotrans
        DO k = 1, ke
          DO j = jstartpar, jendpar
            DO i = istartpar, iendpar
              caero_new(i,j,k,isp) = caero_new(i,j,k,isp)         &
                - ztaur(k)*(caero_new(i,j,k,isp)-zcaerobd(i,j,k,isp))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDIF
#endif

  !$ser verbatim IF (ntstep>0) THEN
  !$ser savepoint RelaxationUnittest.Apply-out LargeTimeStep=ntstep
  !$ser data u=u(:,:,:,nnew) u_tens=utens
  !$ser data v=v(:,:,:,nnew) v_tens=vtens
  !$ser data w=w(:,:,:,nnew)
  !$ser data t=t(:,:,:,nnew)
  !$ser data pp=pp(:,:,:,nnew)
  !$ser tracer %all@nnew
  !$ser verbatim ENDIF ! ntstep>0

!------------------------------------------------------------------------------
!  Section 5: Rayleigh damping for itype_spubc==2 (restore to filtered fields)
!------------------------------------------------------------------------------

  IF (lspubc .AND. (itype_spubc == 2)) THEN 

    ! R.D. with itype_spubc=2  tends to restore the smoothed LM prognostic 
    ! fields. It needs updated gridline halo for the computation of
    ! the large scale component of the LM prognostic fields.

    ! The procedure is the same for the following variables:
    ! u, v, t, pp, qv, qc, qi, qr, qs, qg. w is treated seperately.
    ! The first ke_rd levels from all these fields are stored to field_tmp.
    ! Then all levels of field_tmp are smoothed and exchanged simultaneously.
    num_var = 4

    ! set the dimensions for the boundary exchange
    kzdims(1)    = num_var*ke_rd
    kzdims(2:24) = 0_iintegers

    ! Number of applications of smoother for the determination of the 
    ! large scale field used in the Rayleigh damping with itype_spubc=2

    ALLOCATE( field_tmp(ie,je,num_var*ke_rd), STAT = istata )

    ! Store the first ke_rd levels of all necessary variables to field_tmp
    DO k = 1, ke_rd
      field_tmp(:,:, k          ) =  u(:,:,k,nnew)
      field_tmp(:,:, k +   ke_rd) =  v(:,:,k,nnew)
      field_tmp(:,:, k + 2*ke_rd) =  t(:,:,k,nnew)
      field_tmp(:,:, k + 3*ke_rd) = pp(:,:,k,nnew)
    ENDDO

    ! Do the filtering now on all levels of field_tmp
    DO n = 1, nfi_spubc2
      hfwidth = 1
      CALL exchg_boundaries                                              &
        ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
         ie, je, kzdims, jstartpar, jendpar, hfwidth, nboundlines,       &
         my_cart_neigh, lperi_x, lperi_y, l2dim,                         &
         2000+nexch_tag+1, .FALSE., ncomm_type, izerror,                 &
         yzerrmsg, field_tmp(:,:,:) )
      CALL smooth9( field_tmp(:,:,:), istartpar, iendpar,                &
                    jstartpar, jendpar, ie, je, num_var*ke_rd)
    ENDDO

    ! And do the damping
    DO  k = 1, ke_rd
      DO j = jstartpar, jendpar
        DO i = istartpar, iendpar
           u(i,j,k,nnew) =  u(i,j,k,nnew) -                          &
                      ztaur(k)*( u(i,j,k,nnew)-field_tmp(i,j, k          ))
           v(i,j,k,nnew) =  v(i,j,k,nnew) -                          &
                      ztaur(k)*( v(i,j,k,nnew)-field_tmp(i,j, k +   ke_rd))
           t(i,j,k,nnew) =  t(i,j,k,nnew) -                          &
                      ztaur(k)*( t(i,j,k,nnew)-field_tmp(i,j, k + 2*ke_rd))
          pp(i,j,k,nnew) = pp(i,j,k,nnew) -                          &
                      ztaur(k)*(pp(i,j,k,nnew)-field_tmp(i,j, k + 3*ke_rd))

          ! w is treated extra 
          w(i,j,k+1,nnew)= w (i,j,k+1,nnew) - ztaur(k )*  w (i,j,k+1,nnew)

        ENDDO
      ENDDO
    ENDDO

    ! deallocate the temporary field again
    DEALLOCATE (field_tmp)

    ! Rayleigh damping for itype_spubc==2 (restore to filtered fields)
    ! for the tracers
    
    ! R.D. with itype_spubc=2  tends to restore the smoothed LM prognostic
    ! fields. It needs updated gridline halo for the computation of
    ! the large scale component of the LM prognostic fields.
    ! The first ke_rd levels from all tracers are stored to field_tmp.
    ! Then all levels of field_tmp are smoothed and exchanged simultaneously.

    num_var     = 0_iintegers

    ! number of variables (tracers) undergoing Rayleigh damping
    DO iztrcr = 1, trcr_get_ntrcr()
      IF ( izdamp(iztrcr) == T_DAMP_ON ) THEN
        num_var = num_var + 1
      ENDIF
    ENDDO

    ! set the dimensions for the boundary exchange
    kzdims(1)    = num_var*ke_rd
    kzdims(2:24) = 0_iintegers

    ! Number of applications of smoother for the determination of the
    ! large scale field used in the Rayleigh damping with itype_spubc=2

    ALLOCATE( field_tmp(ie,je,num_var*ke_rd), STAT = istata )

    ! loop over tracers
    iztrcr_damp = 0_iintegers
    DO iztrcr = 1, trcr_get_ntrcr()
      IF ( izdamp(iztrcr) == T_DAMP_ON ) THEN

        ! get pointer to tracer (at nnew)
        CALL trcr_get( izerror, iztrcr, ptr_tlev=nnew, ptr=ztrcr_new )
        IF ( izerror /= 0) THEN
          yzerrmsg = trcr_errorstr(izerror)
          CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
        ENDIF

        ! store the first ke_rd levels of all tracers to field_tmp
        iztrcr_damp = iztrcr_damp + 1
        DO k = 1, ke_rd
          field_tmp ( :,:,k+(iztrcr_damp-1)*ke_rd ) =  ztrcr_new (:,:,k)
        ENDDO

      ENDIF
    ENDDO

    ! Do the filtering now on all levels of field_tmp
    DO n = 1, nfi_spubc2
      hfwidth = 1
      CALL exchg_boundaries                                                  &
        ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,       &
          ie, je, kzdims, jstartpar, jendpar, hfwidth, nboundlines,          &
          my_cart_neigh, lperi_x, lperi_y, l2dim,                            &
          2000+nexch_tag+1, .FALSE., ncomm_type, izerror,                    &
          yzerrmsg, field_tmp(:,:,:) )
      CALL smooth9( field_tmp(:,:,:), istartpar, iendpar,                    &
                               jstartpar, jendpar, ie, je, num_var*ke_rd)
    ENDDO

    ! loop over tracers
    iztrcr_damp = 0_iintegers
    DO iztrcr = 1, trcr_get_ntrcr()
      IF ( izdamp(iztrcr) == T_DAMP_ON ) THEN

        ! get pointer to tracer (at nnew)
        CALL trcr_get( izerror, iztrcr, ptr_tlev=nnew, ptr=ztrcr_new )
        IF ( izerror /= 0) THEN
          yzerrmsg = trcr_errorstr(izerror)
          CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
        ENDIF

        ! do Rayleigh damping (against smoothed field)
        iztrcr_damp = iztrcr_damp + 1
        DO  k = 1, ke_rd
          DO j = jstartpar, jendpar
            DO i = istartpar, iendpar
              ztrcr_new(i,j,k) = ztrcr_new(i,j,k) - ztaur(k)                  &
                               * (ztrcr_new(i,j,k)                            &
                               - field_tmp(i,j,k+(iztrcr_damp-1)*ke_rd))
            ENDDO
          ENDDO
        ENDDO

      ENDIF
    ENDDO

    ! deallocate the temporary field again
    DEALLOCATE (field_tmp)

  ENDIF

!------------------------------------------------------------------------------
!  Section 6: Changes in the saturation conditions
!------------------------------------------------------------------------------

  IF (izdebug > 10) THEN
    PRINT *, '        Changes in the saturation conditions'
  ENDIF

  IF ( ldiabf_satad ) THEN
    ! initialize temperature increment due to latent heat
    DO  k = 1, ke
      tinc_lh(:,:,k) = tinc_lh(:,:,k) - t(:,:,k,nnew)
    ENDDO
  ENDIF

#ifdef NUDGING
  IF (llhn .OR. llhnverif) THEN
    ! For calculation of latent heating rate T has to be stored before
    ! saturation adjustment. 
    CALL get_gs_lheating ('add',1,ke)       ! add T to tt_lheat
  ENDIF
#endif

  IF (lcond) THEN
   IF (itype_gscp < 100 .OR. l2mom_satads) THEN
    CALL trcr_get( izerror, idt_qv, ptr_tlev=nnew, ptr=qv_new )
    IF (izerror /= 0) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get( izerror, idt_qc, ptr_tlev=nnew, ptr=qc_new )
    IF (izerror /= 0) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF

!$ser verbatim ! cannot serialize first step as there is some initialization
!$ser verbatim IF (ntstep>0) THEN
!$ser savepoint SaturationAdjustmentUnittest.Apply-in LargeTimeStep=ntstep
!$ser data t=t(:,:,:,nnew) t_nnow=t(:,:,:,nnow) pp=pp(:,:,:,nnew)
!$ser tracer QV@nnew QC@nnew
!$ser verbatim ENDIF

    DO k = 1, ke
      DO j = jstartpar, jendpar
        DO i = istartpar, iendpar
          zphfe(i,j) = p0(i,j,k) + pp(i,j,k,nnew)
        ENDDO
      ENDDO
      CALL satad( kitbnd, t(:,:,k,nnew), qv_new(:,:,k),           &
                  qc_new(:,:,k), t(:,:,k,nnow), zphfe,            &
                  zr1, zr2, zr3, zr4, zr5, zr6, zr7, zr8,         &
                  b1, b2w, b3, b4w, b234w, rdv, o_m_rdv,          &
                  rvd_m_o, lh_v, cpdr, cp_d,                      &
                  ie, je, istart, iend, jstart, jend )
    ENDDO

!$ser verbatim ! cannot serialize first step as there is some initialization
!$ser verbatim IF (ntstep>0) THEN
!$ser savepoint SaturationAdjustmentUnittest.Apply-out LargeTimeStep=ntstep
!$ser data t=t(:,:,:,nnew)
!$ser tracer QV@nnew QC@nnew
!$ser verbatim ENDIF

   ENDIF
  ENDIF

  IF ( ldiabf_satad ) THEN
    ! compute temperature increment due to latent heat
    DO k = 1, ke
      tinc_lh(:,:,k) = tinc_lh(:,:,k) + t(:,:,k,nnew)
    ENDDO
  ENDIF

#ifdef NUDGING
  IF (llhn .OR. llhnverif) CALL get_gs_lheating ('inc',1,ke) ! inc T to tt_lheat
#endif

ENDIF  ! ntstep < nstop

#ifdef COSMOART
IF (l_cosmo_art) THEN
  IF (lgas) THEN
    DEALLOCATE (zcgasbd)
  ENDIF
  IF (laero) THEN
    DEALLOCATE (zcaerobd)
  ENDIF
ENDIF
#endif

!------------------------------------------------------------------------------
! End of module procedure relaxation
!------------------------------------------------------------------------------

END SUBROUTINE relaxation

!==============================================================================
!==============================================================================
!+ Module procedure in "relaxation" for Asselin time-filtering
!------------------------------------------------------------------------------

SUBROUTINE timefilter()

!------------------------------------------------------------------------------
!
! Description:
!   The Asselin time filter is applied on all prognostic variables
!   at time level nnow.
!
! Method:
!   Time-filtering is applied according to Asselin. The value of alphaass
!   can be used to increase the order of the time-filter according to
!   Williams (2009, MWR)
!
!------------------------------------------------------------------------------

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    i, j, k,           & !  
    iztrcr,            & !
    izerror,           & !  error status
    kitasf,            & !  Number of iterations in the saturation adjustment
                         !  after Asselin filtering at time level nnow
    izdebug              !  for additional debug output

  REAL    (KIND=wp   )     ::  &
    zdd                  ! Asselin filter displacement of midpoint

! Local (automatic) arrays:
! -------------------------
  REAL    (KIND=wp   )     ::  &
    zphfe   (ie,je),   & !
    ztstart (ie,je),   & !
    zr1     (ie,je),   & !
    zr2     (ie,je),   & !
    zr3     (ie,je),   & !
    zr4     (ie,je),   & !
    zr5     (ie,je),   & !
    zr6     (ie,je),   & !
    zr7     (ie,je),   & !
    zr8     (ie,je)      !

  INTEGER  (KIND=iintegers) :: &
    izadv  (trcr_get_ntrcr()), & ! type of advection for tracers
    izsp_adv_lf(trcr_get_ntrcr())

! Local (allocatable) arrays:
! ---------------------------
  CHARACTER (LEN=80)       ::  &
    yzerrmsg

  CHARACTER (LEN=25)       :: yzroutine

! Tracer pointers:
!-----------------
  REAL (KIND=wp),     POINTER :: &
    ztrcr_new  (:,:,:)  => NULL(),   & ! tracer variable at nnew
    ztrcr_now  (:,:,:)  => NULL(),   & ! tracer variable at nnow
    ztrcr_old  (:,:,:)  => NULL(),   & ! tracer variable at nold
    qv_now     (:,:,:)  => NULL(),   & ! QV at tlev=nnow
    qc_now     (:,:,:)  => NULL()      ! QC at tlev=nnow

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine timefilter             
!------------------------------------------------------------------------------

izerror = 0_iintegers
yzroutine = 'timefilter'

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

!------------------------------------------------------------------------------
!  Section 1: Initial setting of local variables
!------------------------------------------------------------------------------

IF ( .NOT. l2tls ) THEN

  IF (izdebug > 10) THEN
    PRINT *, '        Start Asselin time filtering'
  ENDIF

! Number of interpolations in saturation adjustment   

  kitasf = 1

  CALL trcr_meta_get(izerror, T_ADV_ID, izadv)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_meta_get(izerror, 'SP_ADV_LF', izsp_adv_lf)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

!------------------------------------------------------------------------------
!  Section 7: Asselin time filtering in case of Leapfrog
!             Changes in the saturation conditions are taken into account
!------------------------------------------------------------------------------

  DO  k = 1, ke
    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar
        zdd = epsass*( u(i,j,k  ,nnew) - 2.0_wp*u(i,j,k  ,nnow) + u(i,j,k  ,nold) )
        u(i,j,k,nnow)   =  u(i,j,k  ,nnow) + alphaass * zdd
        u(i,j,k,nnew)   =  u(i,j,k  ,nnew) + (alphaass - 1.0_wp) * zdd
        zdd = epsass*( v(i,j,k  ,nnew) - 2.0_wp*v(i,j,k  ,nnow) + v(i,j,k  ,nold) )
        v(i,j,k,nnow)   =  v(i,j,k  ,nnow) + alphaass * zdd
        v(i,j,k,nnew)   =  v(i,j,k  ,nnew) + (alphaass - 1.0_wp) * zdd
        zdd = epsass*( w(i,j,k+1,nnew) - 2.0_wp*w(i,j,k+1,nnow) + w(i,j,k+1,nold) )
        w(i,j,k+1,nnow) =  w(i,j,k+1,nnow) + alphaass * zdd
        w(i,j,k+1,nnew) =  w(i,j,k+1,nnew) + (alphaass - 1.0_wp) * zdd
      ENDDO
    ENDDO

    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar
        zdd = epsass*( pp(i,j,k,nnew) - 2.0_wp*pp(i,j,k,nnow) + pp(i,j,k,nold) )
        pp(i,j,k,nnow) =  pp(i,j,k,nnow) + alphaass * zdd
        pp(i,j,k,nnew) =  pp(i,j,k,nnew) + (alphaass - 1.0_wp) * zdd
      ENDDO
    ENDDO

    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar
        zdd = epsass*(  t(i,j,k,nnew) - 2.0_wp* t(i,j,k,nnow) +  t(i,j,k,nold) )
        t(i,j,k,nnow)  =   t(i,j,k,nnow) + alphaass * zdd
        t(i,j,k,nnew)  =   t(i,j,k,nnew) + (alphaass - 1.0_wp) * zdd
      ENDDO
    ENDDO
  ENDDO

  ! loop over tracers
  DO iztrcr = 1, trcr_get_ntrcr()

    ! check if Asselin filtering is required
    IF ( izsp_adv_lf(iztrcr) == 1_iintegers .AND. izadv(iztrcr) == T_ADV_ON ) THEN

      ! get pointer to tracers (at nnew, nnow and nold)
      CALL trcr_get(izerror, iztrcr, ptr_tlev=nnew, ptr=ztrcr_new)
      IF ( izerror /= 0_iintegers )THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF
      CALL trcr_get(izerror, iztrcr, ptr_tlev=nnow, ptr=ztrcr_now)
      IF ( izerror /= 0_iintegers )THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF
      CALL trcr_get(izerror, iztrcr, ptr_tlev=nold, ptr=ztrcr_old)
      IF ( izerror /= 0_iintegers )THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF

      ! do Asselin filtering
      DO  k = 1, ke
        DO j = jstartpar, jendpar
          DO i = istartpar, iendpar
            zdd = epsass*( ztrcr_new(i,j,k) - 2.0_wp*ztrcr_now(i,j,k)      &
                         + ztrcr_old(i,j,k) )
            ztrcr_now(i,j,k) =  ztrcr_now(i,j,k) + alphaass * zdd
            ztrcr_new(i,j,k) =  ztrcr_new(i,j,k) +                            &
                               (alphaass - 1.0_wp) * zdd
          ENDDO
        ENDDO
      ENDDO

    ENDIF
  ENDDO

  ! saturation adjustment
  IF (lcond) THEN

    ! get pointer to qv and qc (at nnow)
    CALL trcr_get(izerror, idt_qv, ptr_tlev=nnow, ptr=qv_now)
    IF ( izerror /= 0_iintegers )THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get(izerror, idt_qc, ptr_tlev=nnow, ptr=qc_now)
    IF ( izerror /= 0_iintegers )THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF

    ! do saturation adjustment
    DO  k = 1, ke
      ztstart(:,:) = t(:,:,k,nnow)
      zphfe(:,:) = p0(:,:,k) + pp(:,:,k,nnow)
      CALL satad( kitasf, t(:,:,k,nnow), qv_now(:,:,k),           &
                  qc_now(:,:,k), ztstart, zphfe,                  &
                  zr1, zr2, zr3, zr4, zr5, zr6, zr7, zr8,         &
                  b1, b2w, b3, b4w, b234w, rdv, o_m_rdv,          &
                  rvd_m_o, lh_v, cpdr, cp_d,                      &
                  ie, je, istartpar, iendpar, jstartpar, jendpar  )
    ENDDO

  ENDIF

ENDIF

!------------------------------------------------------------------------------
! End of module procedure timefilter
!------------------------------------------------------------------------------

END SUBROUTINE timefilter

!==============================================================================

END MODULE src_relaxation
