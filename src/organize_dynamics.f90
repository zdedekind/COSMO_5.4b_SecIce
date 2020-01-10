!+ External procedure for organizing the dynamical time stepping
!------------------------------------------------------------------------------

SUBROUTINE organize_dynamics (yaction, ierror, yerrmsg, dt_alter, linit)

!------------------------------------------------------------------------------
!
! Description:
!   This procedure is the driving routine for calling the different
!   dynamical modules. It is called first at the beginning of the program just
!   to initialize the NAMELIST and several other variables. Later it is 
!   called during the time-stepping and in the initialization.
!
! Method:
!
! Current Code Owner: DWD, Michael Baldauf
!  phone:  +49  69  8062 2733
!  fax:    +49  69  8062 3721
!  email:  Michael.Baldauf@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.34       1999/12/10 Ulrich Schaettler
!  Initial release
! 1.39       2000/05/03 Ulrich Schaettler
!  Included subroutine for Namelist input and some technical changes.
! 2.9        2001/07/16 Guenther Doms
!  Introduction of new Namelist input parameters in input group INPUT_DYN
!  to control horizontal diffusion (itype_hdiff, hd_corr_t, hd_corr_q and
!  hd_dhmax).
! 2.13       2002/01/18 Ulrich Schaettler
!  Set lmulti_layer to .FALSE., if the physics are not computed
! 2.17       2002/05/08 Ulrich Schaettler
!  Bug correction in the computing of vhmx_cfl
! 2.18       2002/07/16 Ulrich Schaettler
!  Use a1t, a2t from data_fields
! 2.19       2002/10/24 Ulrich Schaettler
!  Adaptation for the 2 time level scheme.
! 3.2        2003/02/07 Ulrich Schaettler
!  Introduced a parameter list for subroutine init_relaxation.
! 3.5        2003/09/02 Ulrich Schaettler
!  Set parameter nincconv=1 if, lphys=.FALSE.
! 3.6        2003/12/11 Ulrich Schaettler
!  Modifications for checking the IOSTAT-value when reading the NAMELIST
! 3.7        2004/02/18 Ulrich Schaettler
!  New Namelist variables for selecting new Runge-Kutta time stepping
! 3.13       2004/12/03 Ulrich Schaettler
!  New Namelist variable itype_spubc for choosing kind of Rayleigh damping
!                                                           (Lucio Torrisi)
!  New Namelist variables lexpl_lbc, rlwidth for explicit formulation of
!  lateral boundary relaxation
!  Initialized new fields for sqrtg_r_* (Jochen Foerstner)
! 3.14       2005/01/25 Jochen Foerstner
!  Possibility of iadv_order = 1/2
! 3.15       2005/03/03 Ulrich Schaettler
!  Replaced FLOAT by REAL
! 3.16       2005/07/22 Ulrich Schaettler
!  New Namelist parameters for the Runge-Kutta scheme;
!  Changed default values for some other parameters
! 3.18       2006/03/03 Ulrich Schaettler
!  Changed treatment of ASCII files for introducing restart possibility
!  Add Namelist variables and action "specnudge" for spectral nudging
!  Add Namelist variables ldiabf_lh, ldyn_bbc
! 3.19       2006/04/25 Ulrich Schaettler
!  Corrected distribution of logbuf
! 3.21       2006/12/04 Ulrich Schaettler
!  Eliminate settings of physical Namelist Input for "dry" runs
!  Read namelist variables lcori_deep, ladv_deep for deep atmosphere (R. Petrik)
!  Read namelist variables nfi_spubc2 for itype_spubc=2 (Lucio Torrisi)
! V3_23        2007/03/30 Michael Baldauf
!  Initialization of the reciprocal square root of G at skalar points,
!  in case of 3TL scheme (necessary for computing integrals)
!  New NL variable itype_lbcqx for lateral boundary treatment of qr,qs,qg
!     (by Jochen Foerstner)
! V4_5         2008/09/10 Ronny Petrik, Michael Baldauf
!  Adaptations to src_integrals
!  Add possibility to combine radiative lateral boundary condition with weak
!  Davies relaxation (Guenther Zaengl)
! V4_8         2009/02/16 Ulrich Schaettler
!  More variables are initialized in init_dynamics
!  Comment print-statement for maximum stable wind velocity
!  Transferred Namelist variables lcori, lmetr, lradlbc from setup
! V4_9         2009/07/16 Guenther Zaengl
!  Add option for potential temperature advection
!  Renamed switch itype_lbcqx to more meaningful itype_outflow_qrsg
!  New NL switch itype_lbc_qrsg for lateral boundary treatment of qr, qs, qg
!  Introduced new mask arrays for t-, q- and u-fields
! V4_11        2009/11/30 Guenther Zaengl
!  Read namelist switches alphaass, ltadv_limiter, itype_w_lbc
! V4_12        2010/05/11 Ulrich Schaettler, Michael Baldauf, Oli Fuhrer
!  Renamed itype_lbc_w itype_bbc_w because of "bottom boundary condition"
!  Calculate metric coefficients, if needed
!  Input of new NL variables for horizontal diffusion; eliminated lhdiff_mask
!  Renamed hd_mask_dcoeff to hd_mask_dcoeff_p
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_15        2010/11/19 Michael Baldauf
!  Replace NAMELIST-Var. lva_impl_dyn  by  y_vert_adv_dyn
! V4_17        2011/02/24 Ulrich Blahak
!  Check the value of rlwidth only if it is needed
! V4_18        2011/05/26 Michael Baldauf
!  Introduced new NL switch y_scalar_advect (which replaces lsl_adv_qx and yef_adv_qx)
! V4_21        2011/12/06 Michael Baldauf
!  Introduced new NL switch l_diff_Smag for Smagorinsky diffusion
! V4_23        2012/05/10 Michael Baldauf, Ulrich Schaettler, Burkhardt Rockel
!  Allow new variants of the Bott-Advection schemes in Namelist Input:
!    BOTT2_STRANG_B, BOTT4_STRANG_B: Strang Splitting only at the bottom (lowest 5 levels)
!    BOTT2_XYZYX, BOTT4_XYZYX: modified sequence compared to the current Strang splitting
!  Call new SR init_grid_metrics, which collects some other initialization routines
!  Removed src_2timelevel.f90 and related option for irunge_kutta (==0)
!  New Namelist parameter nincsn introduced in Namelist-Group "DYNCTL"
!    It defines a time increment for spectral nudging. (BR)
! V4_24        2012/06/22 Michael Baldauf
!  Input of Namelist switch itype_fast_waves
!  Optional call of SR 'init_fast_waves_sc'
! V4_25        2012/09/28 Michael Baldauf
!  Adapted format of YUSPECIF-output
! V4_26        2012/12/06 Anne Roches
!  Replacement of hd_corr_q_XXX by hd_corr_trcr_XXX in order to be consistent
!  also with the naming of other switches (e.g. ltrcr_trilin, lef_adv_trcr_notpd)
! V4_27        2013/03/19 Michael Baldauf, Ulrich Schaettler
!  Read new namelist parameter divdamp_slope for new fast-waves solver
!  Modified call to SR init_fast_waves_sc to init_fast_waves_sc_1
!  Introduced new action "cleanup" for deallocation of fields
!  Modified default values of some tuning constants to reflect settings of COSMO-EU (US)
! V4_28        2013/07/12 Ulrich Blahak
!  Conditional namelist setting of lw_freeslip (Uli Blahak).
!    (no change of results; prevents model breaks  in idealised tests)
! V4_29        2013/10/04 Davide Cesari, Ulrich Schaettler
!  Call to finalize_runge_kutta only for l2tls
!     (otherwise it could crash during Leapfrog-runs)
! V5_1         2014-11-28 Michael Baldauf, Ulrich Blahak, Oliver Fuhrer
!  Introduced new NL variables lhor_pgrad_Mahrer, l_3D_div_damping (MB)
!  Removed Namelist-Parameter lexpl_lbc; replaced crltau by crltau_inv (MB)
!  Changed the format of some YUSPECIF entries for the CLM namelist tool. (UB)
!  Implemented F2003 IOMSG-mechanism for better namelist error messages. (UB)
!  Call to calps, to get consistent values for the surface pressure in ps(:,:,nnew) (US)
!  Replaced ireals by wp (working precision) (OF)
! V5_2         2015-05-21 Ulrich Schaettler
!  Deallocation of memory at the end of the program run
! V5_3         2015-10-09 Ulrich Schaettler (Work from MCH), Michael Baldauf
!  Splitted SR sardass in two subroutines relaxation and timefilter (US)
!  New namelist switch l_euler_dynamics (MB)
!  Added possibility for no scalar advection (y_scalar_advect="NONE") (MB)
!  Removed calculation of surface pressure (nnew) after relaxation, this is
!   now done in lmorg (OF)
!  But put it after calculation of Asselin Filter in leapfrog (US)
! V5_3a        2015-11-24 Ulrich Blahak
!  Corrections for ASCII output in YUSPECIF
! V5_4         2016-03-10 Oliver Fuhrer, Pascal Spoerri, Xavier Lapillonne
!  Configured the default for the cold pool diffusion threshold thresh_cold_pool_d 
!  and added the thresh_cold_pool variable to the distribution buffer.
!  Added simple clipping for SL3 advection scheme.
! V5_4a        2016-05-10 Michael Baldauf, Pascal Spoerri
!  Changed default value for divdamp_slope to 1.0 (because of bug fix in 
!   fast_waves_sc) (MB)
!  Integrated a namelist switch (l_satad_dyn_iter) to allow for a dynamic number 
!   of saturation adjustments. (PS) 
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

USE data_parameters,     ONLY:  wp, iintegers
USE data_constants,      ONLY:  aks2, aks4, pi, r_earth, rvd_m_o, r_d, g, cp_d
USE data_parallel,       ONLY:  my_cart_id, nproc, imp_reals,                 &
                                imp_integers, imp_logical, imp_character,     &
                                nboundlines, my_world_id, icomm_world,        &
                                realbuf, intbuf, logbuf, charbuf
USE data_modelconfig,    ONLY:  ie_tot, je_tot, ke, ke1, dt, dt2, epsass,     &
                                istart, iend, jstart, jend,                   &
                                istartpar, iendpar, jstartpar, jendpar,       &
                                alphaass, betasw, betagw, beta2sw, beta2gw,   &
                                startlat_tot, dlon, dlat, eddlon, eddlat,     &
                                degrad, vhmx_cfl, klv850, ie, je,             &
                                dx_min, dy_min, dz_min, idt_qv, idt_qc
USE data_runcontrol,     ONLY:  lsemi_imp, l2tls, lcpp_dycore, lrubc, lspubc, &
                                rdheight, lcond, nrdtau, ikrylow_si, maxit_si,&
                                eps_si, iprint_si, nusolver, yusolver,        &
                                nuspecif, lw_freeslip, itype_spubc,           &
                                itype_lbc_qrsg, itype_outflow_qrsg, lhordiff, &
                                itype_hdiff, crltau_inv, ntstep, xkd, nstart, &
                                hd_corr_u_bd, hd_corr_t_bd, hd_corr_trcr_bd,  &
                                hd_corr_p_bd, hd_corr_u_in, hd_corr_t_in,     &
                                hd_corr_trcr_in, hd_corr_p_in, hd_dhmax,      &
                                irk_order, iadv_order, irunge_kutta, rlwidth, &
                                intcr_max, y_vert_adv_dyn, nbl_exchg, ldiabf_lh,&
                                y_scalar_advect, ieva_order, ldyn_bbc,        &
                                yvarsn, nzmxsn, nyvar_sn, isc_sn, jsc_sn,     &
                                nincsn, pp_sn, alpha_sn, lspecnudge,          &
                                lcori_deep, ladv_deep, nfi_spubc2, idbg_level,&
                                lprintdeb_all, relax_fac, ldiabf_satad,       &
                                l2dim, lcori, lmetr, lradlbc, lartif_data,    &
                                itheta_adv, itype_bbc_w, ltadv_limiter,       &
                                l_dzeta_d_needed, l_diff_Smag,                &
                                itype_fast_waves, lperi_x, lperi_y,           &
                                divdamp_slope, lhor_pgrad_Mahrer,             &
                                l_3D_div_damping, nnew, l_euler_dynamics,     &
                                l_diff_cold_pools, thresh_cold_pool,          &
                                l_satad_dyn_iter

USE data_fields,         ONLY:  a1t, a2t, hhl, hd_mask, hd_mask_dcoeff_p,     &
                                hd_mask_dcoeff_t, hd_mask_dcoeff_q,           &
                                hd_mask_dcoeff_u, ofa_hdx, ofa_hdy, ps, pp,   &
                                t, t0, rho0, p0, dp0, qrs, rdcoef, dt0dz

USE grid_metrics_utilities, ONLY:  init_grid_metrics, final_grid_metrics
USE parallel_utilities,  ONLY:  distribute_values
USE meteo_utilities,     ONLY:  calps
USE environment,         ONLY:  get_free_unit, model_abort
USE utilities,           ONLY:  to_upper
USE vgrid_refatm_utils,  ONLY:  dealloc_refatm_structure, dealloc_vcoord_structure, &
                                refatm

USE src_leapfrog,        ONLY:  org_leapfrog
USE src_relaxation,      ONLY:  init_relaxation, relaxation, timefilter
USE src_runge_kutta,     ONLY:  org_runge_kutta, finalize_runge_kutta
USE src_spectral_nudging,ONLY:  spectral_nudging, init_spectral_nudging
USE src_advection_rk,    ONLY:  advection_semi_lagrange_init
USE src_tracer,          ONLY:  trcr_get, trcr_errorstr
USE fast_waves_sc,       ONLY:  init_fast_waves_sc_1

#ifdef CPP_DYCORE
USE src_cpp_dycore,      ONLY:  cpp_dycore_compute, cpp_dycore_cleanup
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!                                                           
! Parameter list:
CHARACTER (LEN= *),       INTENT(IN)            ::                      &
  yaction      ! action to be performed

INTEGER (KIND=iintegers), INTENT(OUT)           ::                      &
  ierror       ! error status 

CHARACTER (LEN=  *),      INTENT(OUT)           ::                      &
  yerrmsg      ! error message

REAL (KIND=wp),     INTENT(IN)                  ::                      &
  dt_alter     ! alternative dt: needed for the digital filtering

LOGICAL, INTENT(IN)                             ::                      &
  linit        ! because the PRESENT-function does not work on the T3E

!------------------------------------------------------------------------------
!
! Local variables: 
INTEGER (KIND=iintegers)   ::   izerrstat, nuin, izdebug, i,j,k
CHARACTER (LEN= 9)         ::   yinput
CHARACTER (LEN=17)         ::   yzroutine
CHARACTER (LEN=255)        ::   yzerrmsg
REAL    (KIND=wp  )        ::   zaks2mx, zaks4mx, zdt, zdlon, &
                                zinv_hd, zhd_coeff, z_delta

! Tracer pointers:
REAL (KIND=wp),     POINTER ::                        &
  qv_new     (:,:,:)  => NULL(), & ! QV at nnew
  qc_new     (:,:,:)  => NULL()    ! QC at nnew

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------
!- Begin Subroutine organize_dynamics
!------------------------------------------------------------------------------

  yzroutine = 'organize_dynamics'
  izerrstat = 0_iintegers
  ierror    = 0_iintegers
  yerrmsg   = '    '

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

!------------------------------------------------------------------------------
! Section 1: Input of the NAMELIST Input
!------------------------------------------------------------------------------

  IF (yaction == 'input') THEN

    ! Open NAMELIST-INPUT file
    IF (my_world_id == 0) THEN
      IF (idbg_level > 0) THEN
        PRINT *,'    INPUT OF THE NAMELISTS FOR DYNAMICS'
      ENDIF
      yinput   = 'INPUT_DYN'
      nuin     =  1
      OPEN(nuin   , FILE=yinput  , FORM=  'FORMATTED', STATUS='UNKNOWN',  &
           IOSTAT=izerrstat)
      IF(izerrstat /= 0) THEN
        yerrmsg  = ' ERROR    *** Error while opening file INPUT_DYN *** '
        ierror   = 2
        RETURN
      ENDIF
    ENDIF

    ! read the NAMELIST group dynctl:
    CALL input_dynctl (nuspecif, nuin, izerrstat)

    IF (my_world_id == 0) THEN
      ! Close file for input of the NAMELISTS
      CLOSE (nuin    , STATUS='KEEP')
    ENDIF

    IF (izerrstat < 0) THEN
      yerrmsg = 'ERROR *** while reading NAMELIST Group /DYNCTL/ ***'
      ierror  = 3
    ELSEIF (izerrstat > 0) THEN
      yerrmsg = 'ERROR *** Wrong values occured in NAMELIST INPUT_DYN ***'
      ierror  = 4
    ENDIF

!------------------------------------------------------------------------------
! Section 2: Initialization of the packages
!------------------------------------------------------------------------------

  ELSEIF (yaction == 'init') THEN

    IF (izdebug > 1) THEN
      PRINT *, '    INITIALIZATIONS for DYNAMICS and RELAXATION'
    ENDIF
    ! additional initializations for the dynamics
    CALL init_dynamics   (izerrstat, yzerrmsg)
    CALL init_relaxation (izerrstat, yzerrmsg)
    IF ( itype_fast_waves == 2 ) THEN
      CALL init_fast_waves_sc_1
    END IF
    IF (izerrstat /= 0_iintegers) THEN
      ierror  = 5
      yerrmsg = yzerrmsg
    ENDIF

    IF (lspecnudge) THEN
      IF (izdebug > 1) THEN
        PRINT *, '    INITIALIZATIONS for SPECTRAL NUDGING'
      ENDIF
      CALL init_spectral_nudging
    ENDIF

    IF ( (y_scalar_advect=="SL3_SC") .OR.  &
         (y_scalar_advect=="SL3_MF") .OR.  &
         (y_scalar_advect=="SL3_SFD")   ) THEN
      CALL advection_semi_lagrange_init
    END IF

!------------------------------------------------------------------------------
! Section 3: Dynamical time stepping
!------------------------------------------------------------------------------

  ELSEIF (yaction == 'compute') THEN

  !----------------------------------------------------------------------------
  ! Section 3.1: Initialization of grid- (time-step-) dependent values
  !----------------------------------------------------------------------------


    IF (lcpp_dycore) THEN

      ! FUO TODO: insert a routine that prints the CFL number here

#ifdef CPP_DYCORE
      IF (izdebug > 5) THEN
        PRINT *, '      C++ DYCORE FORTRAN: STEP'
      END IF
      CALL cpp_dycore_compute()
#endif

    ELSE

      ! In the horizontal diffusion, the amplitude of waves with length 2*dx are
      ! damped to 1/E within 2*dt; the 2nd value in the MIN-function is the
      ! critical value for stability
      ! take care for the proper dt:
      IF ( (ntstep == 0) .AND. (.NOT. l2tls) ) THEN
        zdt = 2.0_wp * dt
      ELSE
        zdt = dt
      ENDIF
      IF (linit) THEN
        zdt = ABS (dt_alter)
      ENDIF

      zaks2mx = 1.0_wp / ( 16.0_wp * zdt)
      zaks4mx = 1.0_wp / (128.0_wp * zdt)
      aks2    = MIN ( 1.0_wp / (2.0_wp*zdt*pi**2) , zaks2mx )
      aks4    = MIN ( 1.0_wp / (2.0_wp*zdt*pi**4) , zaks4mx )

      ! Maximum absolute horizontal wind velocity vhmx_cfl which may not
      ! be exceeded for a given large timstep dt.
      zdlon  = r_earth * dlon * degrad
      dx_min = MIN (zdlon * COS( startlat_tot                 *degrad),       &
                    zdlon * COS((startlat_tot+(je_tot-1)*dlat)*degrad) )

      dy_min = r_earth * dlat * degrad

      dz_min = 1.0E30_wp
      DO k=1, ke
        DO j=jstart, jend
          DO i=istart, iend
            z_delta = hhl(i,j,k)-hhl(i,j,k+1)
            dz_min = MIN( dz_min, z_delta )
          END DO
        END DO
      END DO

      vhmx_cfl = MIN( dx_min, dy_min) / (SQRT(2.0_wp)*zdt)   ! not correct for 'Runge-Kutta-'advection

      IF (dt_alter < 0) THEN
        ! this is the backward integration during filtering: set aksx negative
        aks2 = -aks2
        aks4 = -aks4
      ENDIF

      ! Set the hd_mask_coeff matrix
      IF (ntstep == nstart) THEN
        ! for ntstep == 0 and also for restarts

        ! Definition of orographic limiter function
        zinv_hd = 1.0_wp / hd_dhmax
        DO  k = 1, ke
          DO j = 1, je-1
            DO i = 1, ie-1
               ofa_hdx(i,j,k) = MAX ( 0.0_wp, 1.0_wp - &
                                ( (hhl(i+1,j,k+1)-hhl(i,j,k+1)) * zinv_hd  )**2 )
               ofa_hdy(i,j,k) = MAX ( 0.0_wp, 1.0_wp - &
                                ( (hhl(i,j+1,k+1)-hhl(i,j,k+1)) * zinv_hd  )**2 )
            ENDDO
          ENDDO
        ENDDO

        IF (l2tls) THEN
          zhd_coeff = aks4 * dt
        ELSE
          zhd_coeff = aks4 * dt2
        ENDIF
        hd_mask_dcoeff_p(:,:,:) = zhd_coeff *                                  &
                                  ( hd_corr_p_bd *  hd_mask(:,:,:)             &
                                 + (hd_corr_p_in * (1.0_wp - hd_mask(:,:,:)) ))
        hd_mask_dcoeff_t(:,:,:) = zhd_coeff *                                  &
                                  ( hd_corr_t_bd *  hd_mask(:,:,:)             &
                                 + (hd_corr_t_in * (1.0_wp - hd_mask(:,:,:)) ))
        hd_mask_dcoeff_q(:,:,:) = zhd_coeff *                                  &
                                  ( hd_corr_trcr_bd *  hd_mask(:,:,:)          &
                                 + (hd_corr_trcr_in * (1.0_wp - hd_mask(:,:,:)) ))
        hd_mask_dcoeff_u(:,:,:) = zhd_coeff *                                  &
                                  ( hd_corr_u_bd *  hd_mask(:,:,:)             &
                                 + (hd_corr_u_in * (1.0_wp - hd_mask(:,:,:)) ))

      ELSEIF (ntstep == 1 .AND. (.NOT. l2tls)) THEN
        ! only for Leapfrog, because the first step is with dt/2
        hd_mask_dcoeff_p(:,:,:) = hd_mask_dcoeff_p(:,:,:) * 2.0_wp
        hd_mask_dcoeff_t(:,:,:) = hd_mask_dcoeff_t(:,:,:) * 2.0_wp
        hd_mask_dcoeff_q(:,:,:) = hd_mask_dcoeff_q(:,:,:) * 2.0_wp
        hd_mask_dcoeff_u(:,:,:) = hd_mask_dcoeff_u(:,:,:) * 2.0_wp
      ENDIF

      !----------------------------------------------------------------------------
      ! Section 3.2: Call the chosen dynamical package during time-stepping
      !----------------------------------------------------------------------------

      IF (l2tls) THEN
        ! dynamical time stepping with two timelevel scheme
        IF (izdebug > 5) THEN
          PRINT *, '      RUNGE-KUTTA scheme with irunge_kutta = ', irunge_kutta
        ENDIF
        SELECT CASE(irunge_kutta)
        CASE(0)
          ierror   = 22
          yerrmsg  = ' ERROR: This scheme is no more available: irunge_kutta = '
          WRITE (yerrmsg(58:60), '(I3)') irunge_kutta
        CASE(1,2)
          CALL org_runge_kutta
        END SELECT
      ELSE
        IF (izdebug > 5) THEN
          PRINT *, '      LEAPFROG scheme'
        ENDIF
        ! dynamical time stepping with leapfrog scheme
        CALL org_leapfrog
      ENDIF

    ENDIF  ! lcpp_dycore

  !----------------------------------------------------------------------------
  ! Section 3.3: Compute surface pressure for time level nnew
  !----------------------------------------------------------------------------

    ! retrieve the required microphysics tracers
    CALL trcr_get(ierror, idt_qv, ptr_tlev = nnew, ptr = qv_new)
    IF (ierror /= 0) THEN
      yzerrmsg = trcr_errorstr(ierror)
      CALL model_abort(my_cart_id, ierror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get(ierror, idt_qc, ptr_tlev = nnew, ptr = qc_new)
    IF (ierror /= 0) THEN
      yzerrmsg = trcr_errorstr(ierror)
      CALL model_abort(my_cart_id, ierror, yzerrmsg, yzroutine)
    ENDIF

    ! Calculate the surface pressure ps for the new time level nnew
    CALL calps ( ps(:,:   ,nnew), pp(:,:,ke,nnew), t(:,:,ke,nnew),     &
                 qv_new(:,:,ke ), qc_new(:,:,ke ), qrs(:,:,ke)   ,     &
                 rho0(:,:,ke), p0(:,:,ke), dp0(:,:,ke),                &
                 ie, je, rvd_m_o, r_d,                                 &
                 istartpar, iendpar, jstartpar, jendpar )

!------------------------------------------------------------------------------
! Section 4.1: Boundary relaxation
!------------------------------------------------------------------------------

  ELSEIF (yaction == 'relaxation') THEN

    IF (izdebug > 5) THEN
      PRINT *, '      RELAXATION'
    ENDIF

    IF (linit) THEN
      CALL relaxation (ABS(dt_alter))
    ELSE
      CALL relaxation ()
    ENDIF

!   ! Compute surface pressure for time level nnew
!   ! retrieve the required microphysics tracers
!   CALL trcr_get(ierror, idt_qv, ptr_tlev = nnew, ptr = qv_new)
!   IF (ierror /= 0) THEN
!     yzerrmsg = trcr_errorstr(ierror)
!     CALL model_abort(my_cart_id, ierror, yzerrmsg, yzroutine)
!   ENDIF
!   CALL trcr_get(ierror, idt_qc, ptr_tlev = nnew, ptr = qc_new)
!   IF (ierror /= 0) THEN
!     yzerrmsg = trcr_errorstr(ierror)
!     CALL model_abort(my_cart_id, ierror, yzerrmsg, yzroutine)
!   ENDIF

!   ! Calculate the surface pressure ps for the new time level nnew
!   CALL calps ( ps(:,:   ,nnew), pp(:,:,ke,nnew), t(:,:,ke,nnew),     &
!                qv_new(:,:,ke ), qc_new(:,:,ke ), qrs(:,:,ke)   ,     &
!                rho0(:,:,ke), p0(:,:,ke), dp0(:,:,ke),                &
!                ie, je, rvd_m_o, r_d,                                 &
!                istartpar, iendpar, jstartpar, jendpar )

!------------------------------------------------------------------------------
! Section 4.2: Asselin time-filter (only used for LF core)
!------------------------------------------------------------------------------

  ELSEIF (yaction == 'timefilter') THEN

    IF (.NOT. l2tls) THEN
      IF (izdebug > 5) THEN
        PRINT *, '      ASSELIN TIME-FILTER'
      ENDIF

      CALL timefilter()

      ! Compute surface pressure for time level nnew
      ! retrieve the required microphysics tracers
      CALL trcr_get(ierror, idt_qv, ptr_tlev = nnew, ptr = qv_new)
      IF (ierror /= 0) THEN
        yzerrmsg = trcr_errorstr(ierror)
        CALL model_abort(my_cart_id, ierror, yzerrmsg, yzroutine)
      ENDIF
      CALL trcr_get(ierror, idt_qc, ptr_tlev = nnew, ptr = qc_new)
      IF (ierror /= 0) THEN
        yzerrmsg = trcr_errorstr(ierror)
        CALL model_abort(my_cart_id, ierror, yzerrmsg, yzroutine)
      ENDIF

      ! Calculate the surface pressure ps for the new time level nnew
      CALL calps ( ps(:,:   ,nnew), pp(:,:,ke,nnew), t(:,:,ke,nnew),     &
                   qv_new(:,:,ke ), qc_new(:,:,ke ), qrs(:,:,ke)   ,     &
                   rho0(:,:,ke), p0(:,:,ke), dp0(:,:,ke),                &
                   ie, je, rvd_m_o, r_d,                                 &
                   istartpar, iendpar, jstartpar, jendpar )
    ENDIF

!------------------------------------------------------------------------------
! Section 5: Spectral nudging of boundary data
!------------------------------------------------------------------------------

  ELSEIF (yaction == 'specnudge') THEN

    IF (izdebug > 5) THEN
      PRINT *, '      SPECTRAL NUDGING'
    ENDIF

    CALL spectral_nudging

!------------------------------------------------------------------------------
! Section 6: Cleanup at the end
!------------------------------------------------------------------------------

  ELSEIF (yaction == 'cleanup') THEN

    IF (izdebug > 5) THEN
      PRINT *, '      CLEANUP of DYNAMICS'
    ENDIF

    IF (l2tls .AND. (.NOT. lcpp_dycore)) THEN
      ! is not initialized, if lcpp_dycore
      CALL finalize_runge_kutta
    ENDIF

    ! the following calls are necessary also for lcpp_dycore, because initialization
    ! is also done then (but necessary?)

    ! finalization of grid_metrics_utilities
    CALL final_grid_metrics (izerrstat)
    IF (izerrstat /= 0) THEN
      yzerrmsg = 'finalization of grid_metrics failed'
      CALL model_abort(my_cart_id, izerrstat, yzerrmsg, yzroutine)
    ENDIF

    ! Deallocate some variables
    DEALLOCATE (a1t, a2t, STAT=izerrstat)
    IF (izerrstat /= 0) THEN
      yzerrmsg = 'deallocation of a1t,a2t failed'
      CALL model_abort(my_cart_id, izerrstat, yzerrmsg, yzroutine)
    ENDIF

    ! rdcoef was allocated in init_relaxation
    DEALLOCATE (rdcoef,   STAT=izerrstat)
    IF (izerrstat /= 0) THEN
      yzerrmsg = 'deallocation of rdcoef failed'
      CALL model_abort(my_cart_id, izerrstat, yzerrmsg, yzroutine)
    ENDIF

    ! deallocation of vcoord-values
    CALL dealloc_vcoord_structure (izerrstat)
    IF (izerrstat /= 0) THEN
      yzerrmsg = 'deallocation of vcoord-structure failed'
      CALL model_abort(my_cart_id, izerrstat, yzerrmsg, yzroutine)
    ENDIF

    ! deallocation of refatm-values
    CALL dealloc_refatm_structure (izerrstat)
    IF (izerrstat /= 0) THEN
      yzerrmsg = 'deallocation of refatm-structure failed'
      CALL model_abort(my_cart_id, izerrstat, yzerrmsg, yzroutine)
    ENDIF

#ifdef CPP_DYCORE
    ! NOTE: the C++ dycore will always be initialized if -DCPP_DYCORE is set, thus
    !       we also aways have to clean it up
    IF (izdebug > 5) THEN
      PRINT *, '      HP2C DYCORE FORTRAN: CLEANUP'
    END IF
    CALL cpp_dycore_cleanup()
#endif

!------------------------------------------------------------------------------
! Section 7: All other actions are wrong
!------------------------------------------------------------------------------

  ELSE

    ierror  = 1
    yerrmsg = 'ERROR *** No valid action for the dynamics ***'

  ENDIF

!------------------------------------------------------------------------------
! Internal procedures
!------------------------------------------------------------------------------

CONTAINS

!==============================================================================
!+ Internal procedure in "organize_dynamics" for the input of NAMELIST dynctl
!------------------------------------------------------------------------------

SUBROUTINE input_dynctl (nuspecif, nuin, ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine organizes the input of the NAMELIST-group dynctl. 
!   The group dynctl contains variables for the dynamics: the long 
!   time-step dt, parameters for the Asselin-Filter and for the treatment
!   of the sound waves and logical variables for the upper boundary condition.
!
! Method:
!   All variables are initialized with default values and then read in from
!   the file INPUT. The input values are checked for errors and for
!   consistency. If wrong input values are detected the program prints
!   an error message. The program is not stopped in this routine but an
!   error code is returned to the calling routine that aborts the program after
!   reading in all other namelists.
!   In parallel mode, the variables are distributed to all nodes with the
!   environment-routine distribute_values.   
!   Both, default and input values are written to the file YUSPECIF
!   (specification of the run).
!
!------------------------------------------------------------------------------

! Parameters
  INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
    nuspecif,     & ! Unit number for protocolling the task
    nuin            ! Unit number for Namelist INPUT file

  INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
    ierrstat        ! error status variable

! Variables for default values
  REAL (KIND=wp)             ::       &
    epsass_d,      & ! eps for the Asselin-filter
    alphaass_d,    & ! weight for Williams 2009 modification (0.5 < alpha <= 1)
    betasw_d,      & ! beta-variable for treatment of soundwaves    in w
    betagw_d,      & ! beta-variable for treatment of gravity-waves in w
    beta2sw_d,     & ! beta-variable for treatment of soundwaves    in p*, T*
    beta2gw_d,     & ! beta-variable for treatment of gravity-waves in p*, T*
    rdheight_d,    & ! bottom-height of Rayleigh damping layer
    crltau_inv_d,  & ! time factor for relaxation time tau_r = crltau * dt
    rlwidth_d,     & ! width of relaxation layer
    relax_fac_d,   & ! reduction factor for lateral relaxation zone
    xkd_d,         & ! coefficient for divergence damping
    divdamp_slope_d, & ! exceed the theoretical slope stability criterion of
                       ! the divergence damping (only for itype_fast_waves=2)
    eps_si_d,      & ! epsilon for the elliptic solver
    hd_corr_u_bd_d,& ! correction factor for horizontal diffusion flux of u,v,w in boundary zone
    hd_corr_t_bd_d,& ! correction factor for horizontal diffusion flux of t in boundary zone
    hd_corr_trcr_bd_d,& ! correction factor for horizontal diffusion flux of tracers in boundary zone
    hd_corr_p_bd_d,& ! correction factor for horizontal diffusion flux of p in boundary zone
    hd_corr_u_in_d,& ! correction factor for horizontal diffusion flux of u,v,w in domain
    hd_corr_t_in_d,& ! correction factor for horizontal diffusion flux of t in domain
    hd_corr_trcr_in_d,& ! correction factor for horizontal diffusion flux of tracers in domain
    hd_corr_p_in_d,& ! correction factor for horizontal diffusion flux of p in domain
    hd_dhmax_d,    & ! maximum gridpoint height difference for applying
                     ! horizontal diffusion fluxes between them
    pp_sn_d,       & ! lowest pressure level for spectral nudging
    alpha_sn_d,    & ! amplification factor for spectral nudging (<=1.)
    thresh_cold_pool_d ! temperature threshold for targeted cold pools diffusion

  INTEGER (KIND=iintegers)   ::        &
    nrdtau_d,      & ! number of time steps in Rayleigh damping time scale
    ikrylow_si_d,  & ! dimension of Krylow space for elliptic solver
    maxit_si_d,    & ! maximum number of iteration for elliptic solver
    iprint_si_d,   & ! to control whether statistics of the solver are printed
    irunge_kutta_d,& ! =1: use RK scheme,
                     ! =2: use TVD-RK scheme
    irk_order_d,   & ! order of the Runge-Kutta scheme
    iadv_order_d,  & ! order of the horizontal advection scheme for dynamics
    ieva_order_d,  & ! order of the explicit vertical advection scheme
    itheta_adv_d,  & ! =0: use T' (perturbation temperature) for advection
                     ! =1: use theta' (perturbation potential temperature)
                     ! =2: use theta (full potential temperature)
    itype_fast_waves_d, & ! Type of fast waves solver (1=old, 2=new)
    itype_bbc_w_d, & ! Bottom boundary condition for vertical wind (RK only)
                     ! =0/1: RK-like method following iadv_order
                     ! =2/3: differencing following iadv_order without RK stepping
                     ! =4/5: Fourth-order centered differences
                     ! 0/2/4: include linear extrapolation of horizontal wind to sfc
                     ! 1/3/5: no extrapolation of horizontal wind to sfc
                     ! or new nomenclature in the form "1ed"
    intcr_max_d,   & ! max. integer courant number in cr-indepentent advection
    itype_outflow_qrsg_d, & ! type of relaxation treatment for qr, qs, qg
                     ! =1: (default) same treatment as all variables
                     ! =2: no relaxation for qr,qs,qg at outflow boundary points
    itype_lbc_qrsg_d, & ! type of lateral boundary treatment for qr, qs, qg,
                     ! =1: (default) zero-gradient condition
                     ! =2: set qr,qs,qg to 0.0 at the boundaries
                     ! =3: no presetting at the boundaries
                     !     (must not be chosen for Leapfrog applications)
    itype_spubc_d, & ! type of Rayleigh damping (=1: external bc fields, 
                     !                           =2: lm smoothed fields)
    nfi_spubc2_d,  & ! Number of applications of smoother for the determination
                     !  of the large scale field used in the Rayleigh damping
                     !  with itype_spubc=2
    itype_hdiff_d, & ! type of horizontal diffusion (=1: 4th order linear),
                     ! =2: 4th order linear monotonic with orographic limiter)
    nzmxsn_d,      & ! number of variables in default list for spectral nudging
    isc_sn_d,      & ! spectral nudging in i-direction
    jsc_sn_d,      & ! spectral nudging in j-direction
    nincsn_d         ! time increment for spectral nudging

  LOGICAL                    ::       &
    l_euler_dynamics_d,  & ! on/off Euler-solver; however, Coriolis force and tracer 
                           ! advection may be performed
    lsemi_imp_d,         & ! running with semi-implicit or split-explicit scheme,
    l2tls_d,             & ! forecast using a 2-TL integration scheme
    lcpp_dycore_d,       & ! use C++ dycore?
    ltadv_limiter_d,     & ! use limiter for temperature advection (itheta_adv=2 only)
    lhordiff_d,          & ! running with horizontal diffusion
    lcond_d,             & ! forecast with condensation/evaporation
    ldiabf_lh_d,         & ! include diabatic forcing due to latent heat in RK-scheme
    lspubc_d,            & ! with Rayleigh damping in the upper levels
    lrubc_d,             & ! radiative upper boundary condition
                           ! (not implemented yet)
    lw_freeslip_d,       & ! free slip lateral boundary condition for w
    ldyn_bbc_d,          & ! dynamical bottom boundary condition
    lspecnudge_d,        & ! spectral nudging of boundary data
    lcori_deep_d,        & ! take cos(phi) coriolis terms into account
    ladv_deep_d,         & ! use all metric advective terms
    lcori_d,             & ! to run with or without Coriolis force
    lmetr_d,             & ! =.TRUE.:  with metric terms
                           ! =.FALSE.: or without metric terms
    lradlbc_d,           & ! =.TRUE.:  radiative lateral boundary conditions
                           ! =.FALSE.: or with Davies conditions
    l_diff_Smag_d,       & ! =.TRUE.: use Smagorinsky-diffusion for u and v
    lhor_pgrad_Mahrer_d, & ! if =.TRUE., horizontal p-gradients are calculated 
                           ! analogous to Mahrer (1984) (only for itype_fast_waves==2)
    l_3D_div_damping_d,  & ! if =.TRUE., the fully 3D (=isotropic) divergence damping
                           ! is used (only for itype_fast_waves==2)
    l_diff_cold_pools_d, & ! =.TRUE.: use targeted diffusion for cold pools
    l_satad_dyn_iter_d     ! =.TRUE.: Dynamic number of iterations for the saturation adjustment
                           ! in the runge kutta dynamical core. 

  CHARACTER (LEN=10)          ::           &
    y_vert_adv_dyn_d ! choice between several explicit and implicit vertical advection 
                     ! schemes for the dynamic variables

  CHARACTER (LEN=20)          ::           &
    y_scalar_advect_d

  CHARACTER (LEN= 10)         ::           &
    yvarsn_d(nzmxsn), & ! list of the variables for spectral nudging
    tmpyvar, tmpyvar_d

! other local variables
  REAL (KIND=wp)              ::           &
    zdx, zdy, zldomain_x, zldomain_y

  INTEGER (KIND=iintegers)    ::           &
    ierr, iz_err, invar, i,                &
    zistart_tot, zjstart_tot, ziend_tot, zjend_tot

  CHARACTER(LEN=250)         :: iomsg_str

! Define the namelist group
  NAMELIST /dynctl/ epsass, betasw, lspubc, lrubc, rdheight, nrdtau,   &
                    lsemi_imp, l2tls, lhordiff, lcond, lw_freeslip,    &
                    ikrylow_si, xkd, eps_si, maxit_si, iprint_si,      &
                    irunge_kutta, irk_order, iadv_order, itype_hdiff,  &
                    hd_corr_u_bd, hd_corr_t_bd, hd_corr_trcr_bd,       &
                    hd_corr_p_bd, hd_corr_u_in, hd_corr_t_in,          &
                    hd_corr_trcr_in, hd_corr_p_in, hd_dhmax,           &
                    itype_outflow_qrsg, itype_lbc_qrsg, itype_spubc,   &
                    rlwidth, lcpp_dycore,                              &
                    intcr_max, y_vert_adv_dyn, y_scalar_advect,        &
                    ieva_order, ldiabf_lh, betagw, beta2sw, beta2gw,   &
                    ldyn_bbc, lspecnudge, yvarsn, nincsn, isc_sn,      &
                    jsc_sn, pp_sn, alpha_sn, crltau_inv, lcori_deep,   &
                    ladv_deep, nfi_spubc2, relax_fac, lcori, lmetr,    &
                    lradlbc, itheta_adv, itype_bbc_w, alphaass,        &
                    ltadv_limiter, l_diff_Smag, itype_fast_waves,      &
                    divdamp_slope, lhor_pgrad_Mahrer, l_3D_div_damping,&
                    l_euler_dynamics, l_diff_cold_pools,               &
                    thresh_cold_pool, l_satad_dyn_iter


!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE input_dynctl
!------------------------------------------------------------------------------

ierrstat = 0_iintegers
iz_err   = 0_iintegers

IF (my_world_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 1: Initialize the default variables
!------------------------------------------------------------------------------

  epsass_d             = 0.15_wp
  alphaass_d           = 1.0_wp
  betasw_d             = 0.4_wp
  betagw_d             = 0.4_wp
  beta2sw_d            = 0.4_wp
  beta2gw_d            = 0.4_wp
  lspubc_d             = .TRUE.
  lrubc_d              = .FALSE.
  lradlbc_d            = .FALSE.
  lcond_d              = .TRUE.
  ldiabf_lh_d          = .TRUE.
  lhordiff_d           = .TRUE.
  lsemi_imp_d          = .FALSE.
  l_euler_dynamics_d   = .TRUE.
  l2tls_d              = .TRUE.
  lcpp_dycore_d        = .FALSE.
  y_scalar_advect_d    = "BOTT2_STRANG"
  ltadv_limiter_d      = .FALSE.
  y_vert_adv_dyn_d     = "impl2"
  lw_freeslip_d        = .TRUE.
  ldyn_bbc_d           = .FALSE.
  lcori_d              = .TRUE.
  lmetr_d              = .TRUE.
  lcori_deep_d         = .FALSE.
  ladv_deep_d          = .FALSE.
  rdheight_d           = 11000.0_wp
  crltau_inv_d         = 1.0_wp
  rlwidth_d            = 85000.0_wp
  nrdtau_d             = 5
  xkd_d                = 0.1_wp
  eps_si_d             = 1.0E-8_wp
  ikrylow_si_d         = 20
  maxit_si_d           = 200
  iprint_si_d          = 0
  irunge_kutta_d       = 1
  irk_order_d          = 3
  iadv_order_d         = 3
  ieva_order_d         = 3
  itheta_adv_d         = 0
  itype_bbc_w_d        = 114
  intcr_max_d          = 1
  itype_outflow_qrsg_d = 1
  itype_lbc_qrsg_d     = 1
  itype_spubc_d        = 1
  nfi_spubc2_d         = 10
  itype_hdiff_d        = 2
  hd_corr_u_bd_d       = 0.25_wp
  hd_corr_u_in_d       = 0.25_wp
  hd_corr_t_bd_d       = 0.00_wp
  hd_corr_t_in_d       = 0.00_wp
  hd_corr_trcr_bd_d    = 0.0_wp
  hd_corr_trcr_in_d    = 0.0_wp
  hd_corr_p_bd_d       = 0.00_wp
  hd_corr_p_in_d       = 0.00_wp
  hd_dhmax_d           = 250.0_wp
  relax_fac_d          = 0.01_wp
  l_diff_Smag_d        = .TRUE.
  lhor_pgrad_Mahrer_d  = .FALSE.
  l_3D_div_damping_d   = .FALSE.
  itype_fast_waves_d   = 2
  divdamp_slope_d      = 1.0_wp
  l_diff_cold_pools_d  = .TRUE.
  thresh_cold_pool_d   = 10.0_wp
  l_satad_dyn_iter_d    = .TRUE.

  ! variables for the nudging
  lspecnudge_d         = .FALSE.
  yvarsn_d (:)         = '          '
  yvarsn_d ( 1)        = 'U         '; yvarsn_d ( 2)  = 'V         '
  nzmxsn_d             = 2
  isc_sn_d             = 2
  jsc_sn_d             = 2
  nincsn_d             = 1
  pp_sn_d              = 850.0_wp
  alpha_sn_d           = 0.05_wp

!------------------------------------------------------------------------------
!- Section 2: Initialize variables with default
!------------------------------------------------------------------------------

  epsass               = epsass_d
  alphaass             = alphaass_d
  betasw               = betasw_d
  betagw               = betagw_d
  beta2sw              = beta2sw_d
  beta2gw              = beta2gw_d
  lspubc               = lspubc_d
  lrubc                = lrubc_d
  lradlbc              = lradlbc_d
  lcond                = lcond_d
  ldiabf_lh            = ldiabf_lh_d
  lhordiff             = lhordiff_d
  lsemi_imp            = lsemi_imp_d
  l_euler_dynamics     = l_euler_dynamics_d
  l2tls                = l2tls_d
  lcpp_dycore          = lcpp_dycore_d
  y_scalar_advect      = y_scalar_advect_d
  ltadv_limiter        = ltadv_limiter_d
  y_vert_adv_dyn       = y_vert_adv_dyn_d
  lw_freeslip          = lw_freeslip_d
  ldyn_bbc             = ldyn_bbc_d
  lcori                = lcori_d
  lmetr                = lmetr_d
  lcori_deep           = lcori_deep_d
  ladv_deep            = ladv_deep_d
  rdheight             = rdheight_d
  crltau_inv           = crltau_inv_d
  rlwidth              = rlwidth_d
  relax_fac            = relax_fac_d
  nrdtau               = nrdtau_d
  xkd                  = xkd_d
  eps_si               = eps_si_d
  ikrylow_si           = ikrylow_si_d
  maxit_si             = maxit_si_d
  iprint_si            = iprint_si_d
  irunge_kutta         = irunge_kutta_d
  irk_order            = irk_order_d
  iadv_order           = iadv_order_d
  ieva_order           = ieva_order_d
  itheta_adv           = itheta_adv_d
  itype_bbc_w          = itype_bbc_w_d
  intcr_max            = intcr_max_d
  itype_outflow_qrsg   = itype_outflow_qrsg_d
  itype_lbc_qrsg       = itype_lbc_qrsg_d
  itype_spubc          = itype_spubc_d
  nfi_spubc2           = nfi_spubc2_d
  itype_hdiff          = itype_hdiff_d
  hd_corr_u_bd         = hd_corr_u_bd_d
  hd_corr_u_in         = hd_corr_u_in_d
  hd_corr_t_bd         = hd_corr_t_bd_d
  hd_corr_t_in         = hd_corr_t_in_d
  hd_corr_trcr_bd      = hd_corr_trcr_bd_d
  hd_corr_trcr_in      = hd_corr_trcr_in_d
  hd_corr_p_bd         = hd_corr_p_bd_d
  hd_corr_p_in         = hd_corr_p_in_d
  hd_dhmax             = hd_dhmax_d
  l_diff_Smag          = l_diff_Smag_d
  lhor_pgrad_Mahrer    = lhor_pgrad_Mahrer_d
  l_3D_div_damping     = l_3D_div_damping_d
  itype_fast_waves     = itype_fast_waves_d
  divdamp_slope        = divdamp_slope_d
  l_diff_cold_pools    = l_diff_cold_pools_d
  thresh_cold_pool     = thresh_cold_pool_d
  l_satad_dyn_iter     = l_satad_dyn_iter_d

! Variables for spectral nudging
  yvarsn (:)           = ''
  nyvar_sn             = 0
  lspecnudge           = lspecnudge_d
  isc_sn               = isc_sn_d
  jsc_sn               = jsc_sn_d
  nincsn               = nincsn_d
  pp_sn                = pp_sn_d
  alpha_sn             = alpha_sn_d

!------------------------------------------------------------------------------
!- Section 3: Input of the namelist values
!------------------------------------------------------------------------------

  iomsg_str(:) = ' '
  READ (nuin, dynctl, IOSTAT=iz_err, IOMSG=iomsg_str)

  IF (iz_err /= 0) WRITE (*,'(A,A)') 'Namelist-ERROR DYNCTL: ', TRIM(iomsg_str)
ENDIF

IF (nproc > 1) THEN
  ! distribute error status to all processors
  CALL distribute_values  (iz_err, 1, 0, imp_integers,  icomm_world, ierr)
ENDIF

IF (iz_err /= 0) THEN
  ierrstat = -1
  RETURN
ENDIF

IF (my_world_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 4: Check values for errors and consistency
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 4.1: Checks, if C++ DyCore is used
!------------------------------------------------------------------------------

#ifdef GCL_COMM
#ifndef CPP_DYCORE
    PRINT *, ' ERROR  *** GCL requires C++ Dycore and -DCPP_DYCORE at compile time *** '
    ierrstat = 4242
#endif
#endif

  IF ( lcpp_dycore ) THEN

#ifndef CPP_DYCORE
    PRINT *, ' ERROR  *** C++ Dycore requires -DCPP_DYCORE at compile time *** '
    ierrstat = 4241
    lcpp_dycore = .false.
#endif

! FUO_HACK
! Here are some hard-coded check which assure that the user does not
! choose some namelist settings that are not supported by the C++
! Dycore. These should be updated as new functionality is added or new
! restrictions are discovered

! FUO TODO: This list should be updated, especially in the view of the new
!           FW solver support. Also, it is not clear how to handle compile time switches.

    ! check dynamics parameters
    IF (.NOT.l2tls) THEN
      PRINT *, ' ERROR  *** C++ dycore requires l2tls=.true. *** '
      ierrstat = 4250
    ENDIF
    IF (irunge_kutta /= 1) THEN
      PRINT *, ' ERROR  *** C++ dycore requires irunge_kutta=1 *** '
      ierrstat = 4251
    ENDIF
    IF (irk_order /= 3) THEN
      PRINT *, ' ERROR  *** C++ dycore requires irk_order=3 *** '
      ierrstat = 4252
    ENDIF
    IF (iadv_order /= 3 .AND. iadv_order /= 5) THEN
      PRINT *, ' ERROR  *** C++ dycore requires iadv_order=3 or iadv_order=5 *** '
      ierrstat = 4253
    ENDIF
    IF (itheta_adv /= 0) THEN
      PRINT *, ' ERROR  *** C++ dycore requires itheta_adv=0 *** '
      ierrstat = 4254
    ENDIF
    IF (itype_fast_waves == 1 .AND. itype_bbc_w /= 2) THEN
      PRINT *, ' ERROR  *** C++ dycore requires itype_bbc_w=2 for the RK fast waves solver *** '
      ierrstat = 4256
    ENDIF
    IF (itype_fast_waves == 2 .AND. (itype_bbc_w /= 114 .AND. itype_bbc_w /= 124)) THEN
      PRINT *, ' ERROR  *** C++ dycore requires itype_bbc_w=114 or itype_bbc_w=124 for the SC fast waves solver *** '
      ierrstat = 4256
    ENDIF
    IF (lspecnudge) THEN
      PRINT *, ' ERROR  *** C++ dycore requires lspecnudge=.false. *** '
      ierrstat = 4256
    ENDIF
    CALL to_upper( y_scalar_advect )
    IF (TRIM(y_scalar_advect) /= "BOTT2" .AND. TRIM(y_scalar_advect) /= "BOTT2_STRANG") THEN
      PRINT *, ' ERROR  *** C++ dycore requires y_scalar_advect="BOTT2"/"BOTT2_STRANG" *** '
      ierrstat = 4258
    ENDIF
    IF (TRIM(y_vert_adv_dyn) /= "impl2") THEN
      PRINT *, ' ERROR  *** C++ dycore requires y_vert_adv_dyn="impl2" *** '
      ierrstat = 4259
    ENDIF
    IF( itype_spubc /= 1 .AND. itype_spubc /= 3 ) THEN
      PRINT *, ' ERROR  *** C++ dycore requires itype_spubc = 1 or 3 *** '
      ierrstat = 4296
    ENDIF
    IF (lradlbc) THEN
      PRINT *, ' ERROR  *** C++ dycore requires lradlbc=.false. *** '
      ierrstat = 4263
    ENDIF
    IF (.NOT. lw_freeslip) THEN
      PRINT *, ' ERROR  *** C++ dycore requires lw_freeslip=.true. *** '
      ierrstat = 4266
    ENDIF
    IF (lsemi_imp) THEN
      PRINT *, ' ERROR  *** C++ dycore requires lsemi_imp=.false. *** '
      ierrstat = 4267
    ENDIF
    IF (itype_fast_waves == 1 .AND. .NOT. ldyn_bbc) THEN
      PRINT *, ' ERROR  *** C++ dycore requires ldyn_bbc=.true. for the RK fast waves solver *** '
      ierrstat = 4269
    ENDIF
    IF (itype_fast_waves == 2 .AND. ldyn_bbc) THEN
      PRINT *, ' ERROR  *** C++ dycore requires ldyn_bbc=.false. for the SC fast waves solver *** '
      ierrstat = 4269
    ENDIF
    IF (.NOT. lmetr) THEN
      PRINT *, ' ERROR  *** C++ dycore requires lmetr=.true. *** '
      ierrstat = 4271
    ENDIF
    IF (lcori_deep) THEN
      PRINT *, ' ERROR  *** C++ dycore requires lcori_deep=.false. *** '
      ierrstat = 4272
    ENDIF
    IF (ladv_deep) THEN
      PRINT *, ' ERROR  *** C++ dycore requires ladv_deep=.false. *** '
      ierrstat = 4273
    ENDIF
    IF (intcr_max > 1) THEN
      PRINT *, ' ERROR  *** C++ dycore requires intcr_max<=1 *** '
      ierrstat = 4275
    ENDIF
    IF (itype_lbc_qrsg /= 1) THEN
      PRINT *, ' ERROR  *** C++ dycore requires itype_lbc_qrsg=1 *** '
      ierrstat = 4277
    ENDIF
    IF (betasw /= 0.4_wp) THEN
      PRINT *, ' ERROR  *** C++ dycore requires betasw=0.4 *** '
      ierrstat = 4289
    ENDIF
    IF (betagw /= 0.4_wp) THEN
      PRINT *, ' ERROR  *** C++ dycore requires betagw=0.4 *** '
      ierrstat = 4290
    ENDIF
    IF (beta2sw /= 0.4_wp) THEN
      PRINT *, ' ERROR  *** C++ dycore requires beta2sw=0.4 *** '
      ierrstat = 4291
    ENDIF
    IF (beta2gw /= 0.4_wp) THEN
      PRINT *, ' ERROR  *** C++ dycore requires beta2gw=0.4 *** '
      ierrstat = 4292
    ENDIF
    IF (l_satad_dyn_iter) THEN
      PRINT *, ' ERROR *** C++ dycore requires l_satad_dyn_iter=FALSE ***'
      ierrstat = 4293
    ENDIF

  ENDIF ! lcpp_dycore=.TRUE.

!------------------------------------------------------------------------------
!- Section 4.1: Usual checks for Fortran DyCore
!------------------------------------------------------------------------------

  IF ( (divdamp_slope < 0.1_wp) .OR. (divdamp_slope > 3.0_wp) ) THEN
    PRINT *, ' ERROR  *** Wrong value for divdamp_slope : ', divdamp_slope,' *** '
    PRINT *, '        ***   must be >= 0.1 and <= 3 !!!                    *** '
    ierrstat = 1002
  ENDIF
  IF ( (itype_outflow_qrsg < 1) .OR. (itype_outflow_qrsg > 2) ) THEN
    PRINT *, ' ERROR  *** Wrong value for itype_outflow_qrsg: ',itype_outflow_qrsg,' *** '
    PRINT *, '        ***   must be >= 1 and <= 2 !!!                    *** '
    ierrstat = 1002
  ENDIF

  IF ( (itype_lbc_qrsg < 1) .OR. (itype_lbc_qrsg > 3) ) THEN
    PRINT *, ' ERROR  *** Wrong value for itype_lbc_qrsg: ',itype_lbc_qrsg,' *** '
    PRINT *, '        ***   must be >= 1 and <= 3 !!!                    *** '
    ierrstat = 1002
  ENDIF

  IF ( (itype_lbc_qrsg == 3) .AND. (.NOT. l2tls) ) THEN
    PRINT *, ' ERROR  *** Wrong setting for itype_lbc_qrsg: ',itype_lbc_qrsg,' *** '
    PRINT *, '        *** For Leapfrog application itype_lbc_qrsg = 3 is not allowed  *** '
    ierrstat = 1002
  ENDIF

  IF ( (itype_spubc < 1) .OR. (itype_spubc > 3) ) THEN
    PRINT *, ' ERROR  *** Wrong value for itype_spubc: ',itype_spubc,' *** '
    PRINT *, '        ***   must be >= 1 and <= 3 !!!                  *** '
    ierrstat = 1002
  ENDIF

  IF ( (itype_hdiff < 1) .OR. (itype_hdiff > 2) ) THEN
    PRINT *, ' ERROR  *** Wrong value for itype_hdiff: ',itype_hdiff,' *** '
    PRINT *, '        ***   must be >= 1 and <= 2 !!!                  *** '
    ierrstat = 1002
  ENDIF

  IF ( nfi_spubc2 < 1) THEN
    PRINT *, ' ERROR  *** Wrong value for nfi_spubc2: ',nfi_spubc2,' *** '
    PRINT *, '        ***   must be >= 1  !!!                        *** '
    ierrstat = 1002
  ENDIF

  IF ( crltau_inv > 2.0_wp ) THEN
    PRINT *, ' ERROR  *** Wrong value for crltau_inv: ', crltau_inv,' *** '
    PRINT *, '        ***   must be <= 2.0 !!!                        *** '
    ierrstat = 1002
  ENDIF

  IF ( alphaass <= 0.5_wp .OR. alphaass > 1.0_wp ) THEN
    PRINT *, ' ERROR  *** Wrong value for alphaass: ', alphaass,' *** '
    PRINT *, '        ***   must be > 0.5 and <= 1.0 !!!           *** '
    ierrstat = 1002
  ENDIF

  ! plausibility of lateral boundary treatment:

  ! periodic BC's have highest priority:
  IF ( lperi_x .OR. lperi_y ) THEN
    IF ( lradlbc .AND. (my_cart_id == 0) ) THEN
      WRITE(*,*) "WARNING: lradlbc is set to .FALSE.!"
    END IF
    lradlbc = .FALSE.

    IF (lw_freeslip .AND. (my_cart_id == 0) ) THEN
      WRITE(*,*) "WARNING: lw_freeslip is set to .FALSE.!"
    END IF
    lw_freeslip = .FALSE.
  END IF

  ! Check for periodic boundary conditions, metric terms and 2D-version
  IF (lartif_data .EQV. .FALSE.) THEN
    IF (lradlbc .EQV. .TRUE.) THEN
      PRINT *,' ERROR    *** lradlbc = .TRUE. only if lartif_data = .TRUE. *** '
      ierrstat = 1002
    ENDIF
    IF (lmetr .EQV. .FALSE.) THEN
      PRINT *,' ERROR    *** lmetr = .FALSE. only if lartif_data = .TRUE. *** '
      ierrstat = 1002
    ENDIF
  ENDIF

  ! calculate (estimate of) domain size
  zistart_tot =  1     + nboundlines
  zjstart_tot =  1     + nboundlines
  ziend_tot   = ie_tot - nboundlines
  zjend_tot   = je_tot - nboundlines
  zdx = r_earth * dlon * degrad
  zdx = MIN( zdx * COS( (startlat_tot+(zjstart_tot-1)*dlat)*degrad ),       &
             zdx * COS( (startlat_tot+(zjend_tot  -1)*dlat)*degrad ) )
  zdy = r_earth * dlat * degrad
  zldomain_x = ( ziend_tot - zistart_tot + 1 ) * zdx
  zldomain_y = ( zjend_tot - zjstart_tot + 1 ) * zdy
  IF (idbg_level > 1) THEN
    PRINT *, '  DOMAIN SIZE (approx.) in m: L_x = ', zldomain_x
  ENDIF
  IF ( .NOT.l2dim ) THEN
    IF (idbg_level > 1) THEN
      PRINT *, '                              L_y = ', zldomain_y
    ENDIF
    ! zldomain_x = MIN( zldomain_x, zldomain_y )
  ENDIF

  ! check setting of rlwidth
  IF ( .NOT.lperi_x ) THEN
    IF ( rlwidth > 0.25_wp * zldomain_x ) THEN
      PRINT *, ' ERROR  *** Choosen value of rlwidth = ', rlwidth,' is to big *** '
      PRINT *, '        ***   it must be <= 0.25 * L_x !!!        *** '
      ierrstat = 1002
    ENDIF
  END IF
  IF ( .NOT.(lperi_y .OR. l2dim) ) THEN
    IF ( rlwidth > 0.25_wp * zldomain_y ) THEN
      PRINT *, ' ERROR  *** Choosen value of rlwidth = ', rlwidth,' is to big *** '
      PRINT *, '        ***   it must be <= 0.25 * L_y !!!        *** '
      ierrstat = 1002
    ENDIF
  END IF

  IF ( (irunge_kutta < 1) .OR. (irunge_kutta > 2) ) THEN
    PRINT *, ' ERROR  *** Wrong value for irunge_kutta: ',irunge_kutta,' *** '
    PRINT *, '        ***   must be >= 1 and <= 2 !!!                    *** '
    ierrstat = 1002
  ENDIF

  ! check only if new Runge-Kutta dynamical core is used
  IF ( l2tls ) THEN

    SELECT CASE( itheta_adv )
    CASE( 0 )
      ! reference setting... nothing to say.
    CASE( 1, 2 )
      PRINT *, ' *** Use potential temperature for advection *** '
    CASE default
      PRINT *, ' ERROR  *** Wrong value for itheta_adv: ',itheta_adv  ,' *** '
      PRINT *, '        ***   must be >= 0 and <= 2 !!!   *** '
      ierrstat = 1002
    END SELECT

    IF (ltadv_limiter .AND. itheta_adv == 2) THEN
      PRINT *, 'Use limiter for potential temperature for advection'
    ELSE IF (ltadv_limiter) THEN
      PRINT *, 'WARNING: The temperature advection limiter works only'
      PRINT *, 'for itheta_adv = 2 and therefore will not be used here!'
      ltadv_limiter = .FALSE.
    ENDIF

    IF ( itype_bbc_w < 100 ) THEN
      SELECT CASE( itype_bbc_w )
      CASE( 0, 1 )
        ! reference setting... nothing to say.
      CASE( 2, 3, 4, 5 )
        PRINT *, ' *** Use differencing following iadv_order without RK stepping *** '
        PRINT *, ' ***             for bottom boundary condition for w           *** '
      CASE default
        PRINT *, ' ERROR  *** Wrong value for itype_bbc_w: ',itype_bbc_w ,' *** '
        PRINT *, '        ***   must be between 0 and 5!!!   *** '
        ierrstat = 1002
      END SELECT
    ELSE
      ! obviously the alternative nomenclature for itype_bbc_w is used;
      ! range check is done in fast_waves_sc.f90
    ENDIF

    SELECT CASE( irk_order )
    CASE( 3 )
      ! good choice... nothing to say.
    CASE( 1, 2 )
      PRINT *, ' WARNING *** Recommended value for irk_order is 3 *** '
      PRINT *, '         ***   the value: ', irk_order, ' is used !!! *** '
    CASE default
      PRINT *, ' ERROR  *** Wrong value for irk_order: ',irk_order  ,' *** '
      PRINT *, '        ***   must be >= 1 and <= 3 !!!   *** '
      ierrstat = 1002
    END SELECT

    SELECT CASE( iadv_order )
    CASE( 3, 5 )
      ! good choice... nothing to say.
    CASE( 1, 2, 4, 6 )
      PRINT *, ' WARNING *** Recommended value for iadv_order is 5 or 3 *** '
      PRINT *, '         ***   the value: ', iadv_order, ' is used !!! *** '
    CASE default
      PRINT *, ' ERROR  *** Wrong value for iadv_order: ', iadv_order, ' *** '
      PRINT *, '        ***   must be >= 1 and <= 6 !!!   *** '
      ierrstat = 1002
    END SELECT

    IF ( y_vert_adv_dyn == "expl" ) THEN
      SELECT CASE( ieva_order )
      CASE( 3, 5 )
        ! good choice... nothing to say.
      CASE( 1, 2, 4, 6 )
        PRINT *, ' WARNING *** Recommended value for ieva_order is 5 or 3 *** '
        PRINT *, '         ***   the value: ', ieva_order, ' is used !!! *** '
      CASE default
        PRINT *, ' ERROR  *** Wrong value for ieva_order: ', ieva_order, ' *** '
        PRINT *, '        ***   must be >= 1 and <= 6 !!!   *** '
        ierrstat = 1002
      END SELECT
    ENDIF

    CALL to_upper( y_scalar_advect )

    SELECT CASE( TRIM(y_scalar_advect) )
    CASE( "SL3_MF" )
      ! good choice... nothing to say.
    CASE( "SL3_SC", "SL3_SFD")
      PRINT *, ' WARNING *** Instead of the SL3_SC or SL3_SFD-scheme, it is  *** '
      PRINT *, '         *** recommended to use y_scalar_advect="SL3_MF"     *** '
    CASE( "BOTT2", "BOTT2_STRANG", "BOTT4", "BOTT4_STRANG" )
      ! good choice... nothing to say.
    CASE( "BOTT2_XYZYX", "BOTT2_STRANG_B", "BOTT4_XYZYX", "BOTT4_STRANG_B" )
      ! new options for testing
    CASE( "VANLEER", "VANLEER_STRANG" )
      PRINT *, ' WARNING *** Instead of the vanLeer-scheme, it is recommended *** '
      PRINT *, '         *** to use y_scalar_advect="Bott2" or "Bott4"        *** '
    CASE( "PPM", "PPM_STRANG" )
      PRINT *, ' WARNING *** Instead of the PPM-scheme, it is recommended   *** '
      PRINT *, '         *** to use y_scalar_advect="Bott2" or "Bott4"      *** '
    CASE( "NONE" )
      PRINT *, ' WARNING *** tracer advection is completely switched off ***'
    CASE default
      PRINT *, ' ERROR  *** Wrong value for y_scalar_advect: "', y_scalar_advect,'" *** '
      ierrstat = 1002
    END SELECT

    SELECT CASE( intcr_max )
    CASE( 1 )
      ! good choice... nothing to say.
    CASE( 0 )
      PRINT *, ' WARNING *** Recommended value for intcr_max is 1 *** '
      PRINT *, '         ***   the value: ', intcr_max, ' is used !!! *** '
    CASE( 2, 3, 4 )
      PRINT *, ' WARNING *** Recommended value for intcr_max is 1         *** '
      PRINT *, '         ***   INT(Cr)-values>1 are possible for the      *** '
      PRINT *, '         ***     advection of qx but not for the dynamics *** '
      PRINT *, '         ***   the value: ', intcr_max, ' is used nevertheless !!! *** '
    CASE default
      PRINT *, ' ERROR  *** Wrong value for intcr_max: ',intcr_max  ,' *** '
      PRINT *, '        ***   must be >= 0 and <= 4 !!!                *** '
      ierrstat = 1002
    END SELECT

  ELSE

    ldiabf_lh = .FALSE.

  ENDIF

  ! Variables for the spectral nudging
  IF (lspecnudge) THEN

    ! determination of the spectral nudging variables
    invar = COUNT(yvarsn(:) /= '')
    IF( invar == 0) THEN
      ! nothing has been read and the default list is used
      yvarsn(1:nzmxsn_d)  = yvarsn_d(1:nzmxsn_d)
      yvarsn(nzmxsn_d+1:) = ''
      nyvar_sn            = nzmxsn_d
    ELSEIF ( (invar >= 1) .AND. (invar <= nzmxsn) ) THEN
      nyvar_sn         = invar
    ELSEIF (invar > nzmxsn) THEN
      PRINT *,         &
        ' ERROR  *** number of variables for spectral nudging is too big ***'
      ierrstat = 1002
    ELSE
      ! something went wrong: cannot determine number of variables for 
      ! spectral nudging
      PRINT *,         &
       ' ERROR  *** cannot determine number of variables for ',          &
                                                 'spectral nudging ***'
      ierrstat = 1002
    ENDIF

    ! this list has to be checked against the list of variables for
    ! boundary data, which has not been read up to now.
    ! This check is done in init_spectral_nudging

    IF (alpha_sn > 1.0_wp) THEN
      PRINT *,' ERROR    *** alpha_sn = ',alpha_sn,' > 1.0 *** '
      ierrstat = 1002
    ELSEIF (alpha_sn < 0._wp) THEN
      PRINT *,' ERROR    *** alpha_sn = ',alpha_sn,' < 0.0 *** '
      ierrstat = 1002
    ENDIF

  ENDIF

  IF ( itype_fast_waves /= 2 ) THEN
    IF ( lhor_pgrad_Mahrer ) THEN
      PRINT *, ' WARNING *** lhor_pgrad_Mahrer==.TRUE. has no effect if itype_fast_waves/=2  ***'
    END IF
    IF ( l_3D_div_damping ) THEN
      PRINT *, ' WARNING *** l_3D_div_damping==.TRUE. has no effect if itype_fast_waves/=2  ***'
    END IF
  END IF

ENDIF

!------------------------------------------------------------------------------
!- Section 5: Distribute variables to all nodes
!------------------------------------------------------------------------------

IF (nproc > 1) THEN

  IF (my_world_id == 0) THEN
    realbuf ( 1) = epsass
    realbuf ( 2) = betasw
    realbuf ( 3) = rdheight
    realbuf ( 4) = xkd
    realbuf ( 5) = eps_si
    realbuf ( 6) = hd_corr_u_bd
    realbuf ( 7) = hd_corr_t_bd
    realbuf ( 8) = hd_corr_trcr_bd
    realbuf ( 9) = hd_corr_p_bd
    realbuf (10) = hd_corr_u_in
    realbuf (11) = hd_corr_t_in
    realbuf (12) = hd_corr_trcr_in
    realbuf (13) = hd_corr_p_in
    realbuf (14) = hd_dhmax
    realbuf (15) = crltau_inv
    realbuf (16) = rlwidth
    realbuf (17) = pp_sn
    realbuf (18) = alpha_sn
    realbuf (19) = betagw
    realbuf (20) = beta2sw
    realbuf (21) = beta2gw
    realbuf (22) = relax_fac
    realbuf (23) = alphaass
    realbuf (24) = divdamp_slope
    realbuf (25) = thresh_cold_pool

    intbuf  ( 1) = nrdtau
    intbuf  ( 2) = maxit_si
    intbuf  ( 3) = ikrylow_si
    intbuf  ( 4) = iprint_si
    intbuf  ( 5) = irunge_kutta
    intbuf  ( 6) = irk_order
    intbuf  ( 7) = iadv_order
    intbuf  ( 8) = ieva_order
    intbuf  ( 9) = intcr_max
    intbuf  (10) = itype_hdiff
    intbuf  (11) = itype_outflow_qrsg
    intbuf  (12) = itype_lbc_qrsg
    intbuf  (13) = itype_spubc
    intbuf  (14) = isc_sn
    intbuf  (15) = jsc_sn
    intbuf  (16) = nyvar_sn
    intbuf  (17) = nfi_spubc2
    intbuf  (18) = itheta_adv
    intbuf  (19) = itype_bbc_w
    intbuf  (20) = nincsn
    intbuf  (21) = itype_fast_waves

    logbuf  ( 1) = lrubc
    logbuf  ( 2) = lspubc
    logbuf  ( 3) = lcond
    logbuf  ( 4) = lw_freeslip
    logbuf  ( 5) = l2tls
    logbuf  ( 6) = lradlbc
    logbuf  ( 7) = ltadv_limiter
    logbuf  ( 8) = lsemi_imp
    logbuf  ( 9) = lhordiff
    logbuf  (10) = ldiabf_lh
    logbuf  (11) = ldyn_bbc
    logbuf  (12) = lspecnudge
    logbuf  (13) = lcori_deep
    logbuf  (14) = ladv_deep
    logbuf  (15) = lcori
    logbuf  (16) = lmetr
    logbuf  (17) = l_diff_Smag
    logbuf  (18) = lhor_pgrad_Mahrer
    logbuf  (19) = l_3D_div_damping
    logbuf  (20) = lcpp_dycore
    logbuf  (21) = l_euler_dynamics
    logbuf  (22) = l_diff_cold_pools
    logbuf  (23) = l_satad_dyn_iter

    charbuf ( 1) = y_scalar_advect
    charbuf ( 2) = y_vert_adv_dyn
    DO i=1,nzmxsn
      charbuf(2+i) = yvarsn (i)
    ENDDO
  ENDIF

  CALL distribute_values (realbuf, 25, 0, imp_reals,     icomm_world, ierr)
  CALL distribute_values (intbuf , 21, 0, imp_integers,  icomm_world, ierr)
  CALL distribute_values (logbuf , 23, 0, imp_logical,   icomm_world, ierr)
  CALL distribute_values (charbuf,  2+nzmxsn, 0, imp_character,    &
                                                         icomm_world, ierr)

  IF (my_world_id /= 0) THEN
    epsass             = realbuf ( 1)
    betasw             = realbuf ( 2)
    rdheight           = realbuf ( 3)
    xkd                = realbuf ( 4)
    eps_si             = realbuf ( 5)
    hd_corr_u_bd       = realbuf ( 6)
    hd_corr_t_bd       = realbuf ( 7)
    hd_corr_trcr_bd    = realbuf ( 8)
    hd_corr_p_bd       = realbuf ( 9)
    hd_corr_u_in       = realbuf (10)
    hd_corr_t_in       = realbuf (11)
    hd_corr_trcr_in    = realbuf (12)
    hd_corr_p_in       = realbuf (13)
    hd_dhmax           = realbuf (14)
    crltau_inv         = realbuf (15)
    rlwidth            = realbuf (16)
    pp_sn              = realbuf (17)
    alpha_sn           = realbuf (18)
    betagw             = realbuf (19)
    beta2sw            = realbuf (20)
    beta2gw            = realbuf (21)
    relax_fac          = realbuf (22)
    alphaass           = realbuf (23)
    divdamp_slope      = realbuf (24)
    thresh_cold_pool   = realbuf (25)

    nrdtau             = intbuf  ( 1)
    maxit_si           = intbuf  ( 2)
    ikrylow_si         = intbuf  ( 3)
    iprint_si          = intbuf  ( 4)
    irunge_kutta       = intbuf  ( 5)
    irk_order          = intbuf  ( 6)
    iadv_order         = intbuf  ( 7)
    ieva_order         = intbuf  ( 8)
    intcr_max          = intbuf  ( 9)
    itype_hdiff        = intbuf  (10)
    itype_outflow_qrsg = intbuf  (11)
    itype_lbc_qrsg     = intbuf  (12)
    itype_spubc        = intbuf  (13)
    isc_sn             = intbuf  (14)
    jsc_sn             = intbuf  (15)
    nyvar_sn           = intbuf  (16)
    nfi_spubc2         = intbuf  (17)
    itheta_adv         = intbuf  (18)
    itype_bbc_w        = intbuf  (19)
    nincsn             = intbuf  (20)
    itype_fast_waves   = intbuf  (21)

    lrubc              = logbuf  ( 1)
    lspubc             = logbuf  ( 2)
    lcond              = logbuf  ( 3)
    lw_freeslip        = logbuf  ( 4)
    l2tls              = logbuf  ( 5)
    lradlbc            = logbuf  ( 6)
    ltadv_limiter      = logbuf  ( 7)
    lsemi_imp          = logbuf  ( 8)
    lhordiff           = logbuf  ( 9)
    ldiabf_lh          = logbuf  (10)
    ldyn_bbc           = logbuf  (11)
    lspecnudge         = logbuf  (12)
    lcori_deep         = logbuf  (13)
    ladv_deep          = logbuf  (14)
    lcori              = logbuf  (15)
    lmetr              = logbuf  (16)
    l_diff_Smag        = logbuf  (17)
    lhor_pgrad_Mahrer  = logbuf  (18)
    l_3D_div_damping   = logbuf  (19)
    lcpp_dycore        = logbuf  (20)
    l_euler_dynamics   = logbuf  (21)
    l_diff_cold_pools  = logbuf  (22)
    l_satad_dyn_iter   = logbuf  (23)

    y_scalar_advect    = TRIM(charbuf ( 1)(1:20))
    y_vert_adv_dyn     = TRIM(charbuf ( 2)(1:10))
    DO i=1,nyvar_sn
      yvarsn(i)        = TRIM(charbuf(2+i)(1:10))
    ENDDO
  ENDIF

ENDIF

!------------------------------------------------------------------------------
!- Section 6: Output of the namelist variables and their default values
!------------------------------------------------------------------------------

IF (my_world_id == 0) THEN

  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(A23)') '0     NAMELIST:  dynctl'
  WRITE (nuspecif, '(A23)') '      -----------------'
  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(T7,A,T33,A,T52,A,T70,A)') 'Variable', 'Actual Value',   &
                                               'Default Value', 'Format'

  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                                           'epsass ',epsass ,epsass_d ,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                                     'alphaass ',alphaass ,alphaass_d ,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                                           'betasw ',betasw ,betasw_d ,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                                           'betagw ',betagw ,betagw_d ,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                                        'beta2sw ',beta2sw ,beta2sw_d ,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                                        'beta2gw ',beta2gw ,beta2gw_d ,' R '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                           'lrubc  ',lrubc  ,lrubc_d  ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                           'lspubc ',lspubc ,lspubc_d ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                       'lradlbc',lradlbc, lradlbc_d   ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                          'lcond', lcond, lcond_d     ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                              'ldiabf_lh', ldiabf_lh, ldiabf_lh_d     ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                         'lw_freeslip', lw_freeslip, lw_freeslip_d    ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
             'l_euler_dynamics', l_euler_dynamics, l_euler_dynamics_d ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                          'l2tls', l2tls, l2tls_d     ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                        'lcpp_dycore', lcpp_dycore, lcpp_dycore_d     ,' L '
  WRITE (nuspecif, '(T8,A,T33,A12  ,T52,A12  ,T71,A3)')                      &
       'y_vert_adv_dyn', TRIM(y_vert_adv_dyn), TRIM(y_vert_adv_dyn_d) ,'C*8'
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                'lsemi_imp', lsemi_imp, lsemi_imp_d   ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                'lhordiff',  lhordiff , lhordiff_d    ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                  'ldyn_bbc ',ldyn_bbc, ldyn_bbc_d    ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                             'lcori',lcori, lcori_d   ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                             'lmetr',lmetr, lmetr_d   ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                               'lcori_deep', lcori_deep, lcori_deep_d ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                               'ladv_deep' , ladv_deep , ladv_deep_d  ,' L '
  WRITE (nuspecif, '(T8,A,T33,A14  ,T52,A14  ,T71,A3)')                      &
    'y_scalar_advect', TRIM(y_scalar_advect), TRIM(y_scalar_advect_d)  ,'C*8'
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                                      'rdheight', rdheight, rdheight_d,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                              'crltau_inv', crltau_inv, crltau_inv_d  ,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                                         'rlwidth', rlwidth, rlwidth_d,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                                   'relax_fac', relax_fac, relax_fac_d,' R '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                              'nrdtau',nrdtau,nrdtau_d,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                     'maxit_si ',maxit_si ,maxit_si_d ,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                               'ikrylow_si', ikrylow_si, ikrylow_si_d ,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                               'iprint_si ', iprint_si , iprint_si_d  ,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                         'irunge_kutta ', irunge_kutta, irunge_kutta_d,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                'irk_order ', irk_order, irk_order_d  ,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                              'iadv_order', iadv_order, iadv_order_d  ,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                              'ieva_order', ieva_order, ieva_order_d  ,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                             'itheta_adv ', itheta_adv, itheta_adv_d  ,' I '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                    'ltadv_limiter ', ltadv_limiter, ltadv_limiter_d  ,' L '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                          'itype_bbc_w ', itype_bbc_w, itype_bbc_w_d  ,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                 'intcr_max', intcr_max, intcr_max_d  ,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
       'itype_outflow_qrsg', itype_outflow_qrsg, itype_outflow_qrsg_d ,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                   'itype_lbc_qrsg', itype_lbc_qrsg, itype_lbc_qrsg_d ,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                            'itype_spubc', itype_spubc, itype_spubc_d ,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                               'nfi_spubc2', nfi_spubc2, nfi_spubc2_d ,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                           'itype_hdiff', itype_hdiff, itype_hdiff_d  ,' I '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                           'l_diff_Smag',l_diff_Smag, l_diff_Smag_d   ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
         'l_diff_cold_pools',l_diff_cold_pools, l_diff_cold_pools_d   ,' L '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
             'thresh_cold_pool', thresh_cold_pool, thresh_cold_pool_d, ' R '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
         'l_satad_dyn_iter',l_satad_dyn_iter, l_satad_dyn_iter_d      ,' L '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                            'hd_corr_u_bd',hd_corr_u_bd,hd_corr_u_bd_d,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                            'hd_corr_t_bd',hd_corr_t_bd,hd_corr_t_bd_d,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                   'hd_corr_trcr_bd',hd_corr_trcr_bd,hd_corr_trcr_bd_d,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                            'hd_corr_p_bd',hd_corr_p_bd,hd_corr_p_bd_d,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                            'hd_corr_u_in',hd_corr_u_in,hd_corr_u_in_d,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                            'hd_corr_t_in',hd_corr_t_in,hd_corr_t_in_d,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                   'hd_corr_trcr_in',hd_corr_trcr_in,hd_corr_trcr_in_d,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                            'hd_corr_p_in',hd_corr_p_in,hd_corr_p_in_d,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                                        'hd_dhmax',hd_dhmax,hd_dhmax_d,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                                              'xkd    ', xkd, xkd_d   ,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                       'divdamp_slope', divdamp_slope, divdamp_slope_d,' R '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
            'lhor_pgrad_Mahrer',lhor_pgrad_Mahrer, lhor_pgrad_Mahrer_d,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
               'l_3D_div_damping',l_3D_div_damping, l_3D_div_damping_d,' L '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                                         'eps_si', eps_si, eps_si_d   ,' R '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
              'itype_fast_waves', itype_fast_waves, itype_fast_waves_d,' I '
  WRITE (nuspecif, '(A2)')  '  '

  ! variables for the spectral nudging
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A)')                       &
                  'lspecnudge    ',lspecnudge    , lspecnudge_d      , ' L '

  WRITE(nuspecif, '(A)')  'Variables for spectral nudging'
  DO i=1,MAX(nyvar_sn,nzmxsn_d)
    tmpyvar = yvarsn(i)
    IF (LEN_TRIM(yvarsn(i)) == 0) THEN
      tmpyvar = '-'
    END IF
    tmpyvar_d = yvarsn_d(i)
    IF (LEN_TRIM(yvarsn_d(i)) == 0) THEN
      tmpyvar_d = '-'
    END IF
    WRITE (nuspecif, '(T8,A,I3.3,A,T33,A10,T52,A10,T71,A)') 'yvarsn(',i,')',&
                                                tmpyvar, tmpyvar_d, 'C*10'
  ENDDO

  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                              'isc_sn      ',isc_sn   ,isc_sn_d ,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                              'jsc_sn      ',jsc_sn   ,jsc_sn_d ,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                              'nincsn      ',nincsn   ,nincsn_d ,' I '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                              'pp_sn       ',pp_sn    ,pp_sn_d  ,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                              'alpha_sn    ',alpha_sn ,alpha_sn_d  ,' R '

  WRITE (nuspecif, '(A2)')  '  '

ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE input_dynctl

!==============================================================================
!+ Internal procedure in "organize_dynamics" for initializing the dynamics
!------------------------------------------------------------------------------

SUBROUTINE init_dynamics (ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "organize_dynamics" is the initializing 
!   routine of the dynamical modules.
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
! --------------
  INTEGER (KIND=iintegers) ::  &
    i, j, k,             & !  loop indices
    istat, niostat         !  error status variable

  REAL    (KIND=wp       ) ::  &
    zaiw_surf, zaiw_back, zdlon, z_delta

  REAL    (KIND=wp   )     ::  &
    rovcp, rp00, govcp

! End of header
!==============================================================================
 
!------------------------------------------------------------------------------
! Begin Subroutine init_dynamics   
!------------------------------------------------------------------------------

  ! Initialize some variables
  ldiabf_satad = .FALSE.
  nbl_exchg    = 2_iintegers

  ! Set constants for potential temperature advection
  rovcp = r_d/cp_d
  rp00  = 1.E-5_wp
  govcp = g/cp_d

  IF ( l2tls ) THEN

    ! Set number of boundlines to exchange when using 2TL scheme
    SELECT CASE(iadv_order)
    CASE( 5, 6 )
      nbl_exchg = 3_iintegers
    END SELECT

    SELECT CASE( TRIM(y_scalar_advect) )
    CASE( "SL3_SC", "SL3_MF", "SL3_SFD" )
      nbl_exchg = MAX( nbl_exchg, intcr_max+2_iintegers )
    CASE( "VANLEER", "VANLEER_STRANG" )
      nbl_exchg = MAX( nbl_exchg, intcr_max+2_iintegers )
    CASE( "PPM", "PPM_STRANG" )
      nbl_exchg = MAX( nbl_exchg, intcr_max+3_iintegers )
    CASE( "BOTT2", "BOTT2_STRANG", "BOTT2_XYZYX", "BOTT2_STRANG_B" )
      nbl_exchg = MAX( nbl_exchg, intcr_max+2_iintegers )
    CASE( "BOTT4", "BOTT4_STRANG", "BOTT4_XYZYX", "BOTT4_STRANG_B" )
      nbl_exchg = MAX( nbl_exchg, intcr_max+3_iintegers )
    END SELECT

    IF ( nbl_exchg > nboundlines ) THEN
      IF (my_cart_id ==0) THEN
        PRINT *,' ERROR ***   VALUE OF NBOUNDLINES IS TOO SMALL   *** '
        PRINT *,'       *** ( CR-INDEP. ADVECTION OF QX and TKE ) *** '
        PRINT *,'       *** NECESSARY VALUE OF NBOUNDLINES: ', nbl_exchg,'     *** '
        PRINT *,'       *** FOR  RUNGE-KUTTA ADVECTION OPERATORS. *** '
        CALL model_abort (my_cart_id, 1001, 'wrong nboundlines', 'init_dynamics')
      END IF
    END IF

  ENDIF

  ! Allocation of the fields a1t, a2t
  ALLOCATE ( a1t    (ke1)      , STAT=istat )
  ALLOCATE ( a2t    (ke1)      , STAT=istat )

  ! Vertical varying implicitness of vertical diffusion (def. at half levels)
  zaiw_back = 0.75_wp  ! background value for implicit weight above 850 hpa
  zaiw_surf = 1.20_wp  ! value for implicit weight at the surface
  DO k = 1,klv850
    a1t (k) = zaiw_back  
  ENDDO
  DO k = klv850+1, ke
    a1t (k) = zaiw_back + (zaiw_surf-zaiw_back) * (k-klv850)  /           &
                                                REAL(ke-klv850, wp)
  ENDDO
  a1t (ke1) = zaiw_surf
  a2t(:)    = 1.0_wp - a1t(:)

  ! Maximum absolute horizontal wind velocity vhmx_cfl which may not
  ! be exceeded for a given large timstep dt.
  zdlon  = r_earth * dlon * degrad
  dx_min = MIN (zdlon * COS( startlat_tot                 *degrad),       &
                zdlon * COS((startlat_tot+(je_tot-1)*dlat)*degrad) )

  dy_min = r_earth * dlat * degrad

  dz_min = 1.0E30_wp
  DO k=1, ke
    DO j=jstart, jend
      DO i=istart, iend
        z_delta = hhl(i,j,k)-hhl(i,j,k+1)
        dz_min = MIN( dz_min, z_delta )
      END DO
    END DO
  END DO

  vhmx_cfl = dx_min / (SQRT(2.0_wp)*dt)

  IF (my_cart_id ==0) THEN
    PRINT *, ' Minimum horizontal grid spacing (dx, dy): ', dx_min, dy_min
    PRINT *, ' Minimum vertical   grid spacing (dz)    : ', dz_min
    PRINT *, ' Value of big time step used    : ', dt
!US PRINT *, ' Maximum stable horizontal wind velocity (CFL): ', vhmx_cfl
  ENDIF

  ! Open file for solver statistics in the semi-implicit scheme, if required
  IF (lsemi_imp) THEN
    IF (iprint_si > 0) THEN
      CALL get_free_unit (nusolver)
    ENDIF
    IF ( (iprint_si > 0) .AND. (my_cart_id == 0) ) THEN
      ! Open file YUSOLVER
      OPEN(nusolver, FILE=yusolver, FORM=  'FORMATTED', STATUS='UNKNOWN',     &
           IOSTAT=niostat)
      IF(niostat /= 0) THEN
        yerrmsg = ' ERROR    *** Error while opening file YUSOLVER *** '
        ierror = 4
        RETURN
      ENDIF

      ! Write semi-implicit characteristics
      WRITE(nusolver,'(A)')     ' Some characteristics for the elliptic solver'
      WRITE(nusolver,'(A,I5)')  '  Maximum number of iterations:  ', maxit_si
      WRITE(nusolver,'(A,I5)')  '  Dimension of the Krylow space: ', ikrylow_si
      WRITE(nusolver,'(A,F20.17)')  '  Precision: eps =               ', eps_si
      CLOSE (nusolver, STATUS='KEEP')
    ELSEIF ( (iprint_si < 0) .AND. (my_cart_id == 0) ) THEN
      ! Print some characteristics to standard out
      PRINT *, ' Some characteristics for the elliptic solver'
      PRINT *, '  Maximum number of iterations:  ', maxit_si
      PRINT *, '  Dimension of the Krylow space: ', ikrylow_si 
      PRINT *, '  Precision: eps =               ', eps_si
    ENDIF 
  ENDIF


  CALL init_grid_metrics (istat)
  IF (istat /= 0) THEN
    IF (my_cart_id == 0) THEN
      PRINT *, ' Initialization of Grid Metrics failed'
    ENDIF
    ierror = istat
  ENDIF

  ! FUO_NOTE: the following section has been moved out of src_runge_kutta.f90
  !           where it has been executed for IF ntstep==nstart

  ! compute reference profiles for RK-dynamical core
  IF ( l2tls .AND. ANY(irunge_kutta == (/1,2/)) ) THEN
    IF ( itheta_adv == 0 ) THEN
      ! calculate the gradient of t0

      IF (refatm%irefatm == 1) THEN
        DO  k = 1, ke
          dt0dz(:,:,k) = - refatm%dt0lp * g * rho0(:,:,k) / p0(:,:,k)
        ENDDO
      ELSE IF ( refatm%irefatm == 2 ) THEN
        DO  k = 1, ke
          dt0dz(:,:,k) =  -refatm%delta_t/refatm%h_scal*EXP(-0.5_wp*(hhl(:,:,k)+ &
                           hhl(:,:,k+1))/refatm%h_scal)
        ENDDO
      ELSE IF ( refatm%irefatm == 3 ) THEN
        DO  k = 1, ke
          dt0dz(:,:,k) = t0(:,:,k)*refatm%bvref**2/g - govcp
        ENDDO
      END IF

    ELSE IF ( itheta_adv >= 1 ) THEN
      ! when potential temperature is used, t0 and dt0dz carry
      ! theta_0 and dtheta_0/dz

      IF ( refatm%irefatm == 1 ) THEN
        DO  k = 1, ke
          t0(:,:,k)    = p0(:,:,k) / ( r_d * rho0(:,:,k) )* &
            EXP(-rovcp*LOG(p0(:,:,k)*rp00))
          dt0dz(:,:,k) = (govcp - refatm%dt0lp * g * rho0(:,:,k) / p0(:,:,k))* &
            EXP(-rovcp*LOG(p0(:,:,k)*rp00))
        ENDDO
      ELSE IF ( refatm%irefatm == 2 ) THEN
        DO  k = 1, ke
          t0(:,:,k)    = p0(:,:,k) / ( r_d * rho0(:,:,k) )* &
            EXP(-rovcp*LOG(p0(:,:,k)*rp00))
          dt0dz(:,:,k) =  (govcp - refatm%delta_t/refatm%h_scal*EXP(-0.5_wp*(hhl(:,:,k)+ &
            hhl(:,:,k+1))/refatm%h_scal)) * EXP(-rovcp*LOG(p0(:,:,k)*rp00))
        ENDDO
      ELSE IF ( refatm%irefatm == 3 ) THEN
        DO  k = 1, ke
          t0(:,:,k)    = p0(:,:,k) / ( r_d * rho0(:,:,k) )* &
            EXP(-rovcp*LOG(p0(:,:,k)*rp00))
          dt0dz(:,:,k) = t0(:,:,k)*refatm%bvref**2/g
        ENDDO
      ENDIF

    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! End of module procedure init_dynamics
!------------------------------------------------------------------------------

END SUBROUTINE init_dynamics

!==============================================================================

!------------------------------------------------------------------------------
! End of module procedure organize_dynamics
!------------------------------------------------------------------------------

END SUBROUTINE organize_dynamics
