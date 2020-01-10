!+ External procedure for computing values near the surface
!------------------------------------------------------------------------------

!option! -pvctl ifopt
SUBROUTINE near_surface ( nx )

!------------------------------------------------------------------------------
!
! Description:
!   This routine calculates values near the surface.
!   - winds in 10 m
!   - temperature, dew-point temperature and specific water vapor content in 2m
!   - accumulation of precipitation rates
!   - minimal and maximal temperature in 2m
!   - maximal expected squall
!   - mean values over forecast for solar and thermal heating and radiation
!   - in case of llm, some variables are set to the lowest level or to some
!     unrealistic values
!
! Method:
!   See Comments in the Sections.
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8236 1493
!  email:  uschaettler@dwd.d400.de
!
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.29       1999/05/11 Ulrich Schaettler
!  Initial release
! 1.30       1999/06/24 Matthias Raschendofer
!  Alternative call of the diagnosis of 2m- and 10m-values: call synop_diag.
! 1.32       1999/08/24 Guenther Doms
!  Reset of top_con and bas_con at every hour; 
!  direct inclusion of former routines 'low_winds_temp' and 'synop_diag'.
! 1.33       1999/10/14 Matthias Raschendorfer
!  The former routine 'synop_diag' is now included in routine 'turbtran' in
!  a more sophisticated way. Here, the previous code is left only for the
!  case of not running the turbulence routines.
! 2.8        2001/07/06 Ulrich Schaettler
!  Added summarization of surface fluxes
! 2.18       2002/07/16 Ulrich Schaettler
!  Changed treatment of dew-point temperature (by Lucio Torrisi, UGM Rome)
! 2.19       2002/10/24 Ulrich Schaettler
!  Moved re-initialization of vmax_10m to lmorg (initialize_loop)
! 3.5        2003/09/02 Ulrich Schaettler
!  Compute surface pressure here (instead of in the relaxation);
!  Thus the communication for ps can be avoided
! 3.13       2004/12/03 Thorsten Reinhardt
!  Changes for the new graupel scheme
! 3.15       2005/03/03 Jan-Peter Schulz
!  Adapt calculation of vmax_10m to the GME 40km/40L formulation
! 3.18       2006/03/03 Matthias Raschendorfer, Ulrich Schaettler
!  Introduction of rh_2m; Additional fields for Climate-LM Version.
!  Near surace levels also possible between the lowest two model levels
!  in the case of 'itype_synd==1'.
! 3.21       2006/12/04 Burkhardt Rockel, Ulrich Schaettler
!  Renamed sunshhrs, sodwdir to dursun, sodwddm
!  Do not calculate summations and meanvalues for ntstep==0
! V3_23        2007/03/30 Matthias Raschendorfer
!  Moved 'akt' to MOULE data_turbulence.
!  Accumulation of fields for topographic radiation correction (Matteo Buzzi)
! V3_24        2007/04/26 Ulrich Schaettler
!  Eliminated nincmxt and introduced control as for other increments
! V4_1         2007/12/04 Jan-Peter Schulz
!  Re-tuning of gusts by using wind speed at 10 m instead of 30 m
! V4_4         2008/07/16 Jan-Peter Schulz
!  Accumulation of fields for sub-grid scale orography scheme
! V4_8         2009/02/16 Ulrich Schaettler
!  Compute additional averaged values for radiation
!  Introduced several options for wind gusts
! V4_9         2009/07/16 Ulrich Schaettler, Christian Bollmann
!  Inserted a NEC compiler option directive
! V4_11        2009/11/30 Lucio Torrisi
!  Computation of averages for additional fields
! V4_12        2010/05/11 Oliver Fuhrer
!  Additional computations for sunshine duration
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_17        2011/02/24 Ulrich Blahak
!  Fixed index bounds (istartpar/jstartpar instead of istart/jstart) in computation
!  of zum1, zum2, ...
! V4_18        2011/05/26 Ulrich Schaettler
!  After adjusting all loop boundaries, the additional loops for setting some
!  fields at the boundary zone can be eliminated
!  Bug fix for itype_diag_gusts=2: use values in 10m instead of 30m (Oli Fuhrer)
! V4_21        2011/12/06 Jan-Peter Schulz
!  Introduce gust option itype_diag_gusts = 4. Here the gust factor weakly
!  depends on the mean wind speed at 10 m.
! V4_23        2012/05/10 Burkhardt Rockel (CLM)
!  Introduction of new diagnostic variable for maximum wind speed in 10m height
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer, Ulrich Blahak
!  Replaced qx-variables by using them from the tracer module
!  For the 2-moment microphysics scheme, added hail_gsp and prh_gsp (UB).
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  Bug Fix: summation of fluxes in multi-layer soil must only be done,
!    if lmulti_layer is TRUE
!  Introduced MESSy interface
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
! V5_4a        2016-05-10 Matthias Raschendorfer
!  Use akt from new module turb_data instead of data_turbulence
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

USE data_parameters,    ONLY:  wp, iintegers

USE data_fields,        ONLY:                                        &
       u, v, t, pp, p0, hsurf, hhl, gz0, tcm, tch, u_10m,            &
       v_10m, t_2m, td_2m, rh_2m, tmin_2m, tmax_2m, qv_2m, vmax_10m, &
       ps, t_g, qv_s, rain_gsp, rain_con, snow_gsp, snow_con,        &
       grau_gsp, thbs, pabs,                                         &
       apab_s, asob_s, asob_t, athb_s, athb_t, sobs, sobt, thbt,     &
       prg_gsp, prs_gsp, prs_con, prr_gsp, prr_con,tfh,tfm, top_con, &
       bas_con,                                                      &
       umfl_s, vmfl_s, shfl_s, lhfl_s, aumfl_s, avmfl_s, ashfl_s,    &
       alhfl_s, rho0, qrs, dpsdt, dp0, fr_land,                      &
       t_2m_av, td_2m_av, u_10m_av, v_10m_av, sodwddm, dursun,       &
       swdir_s, swdifd_s, swdifu_s, lwd_s, lwu_s, aswdir_s,          &
       aswdifd_s, aswdifu_s, alwd_s, alwu_s,                         &
       lhfl_bs,lhfl_pl,alhfl_bs,alhfl_pl,                            &
       ustr_sso, vstr_sso, vdis_sso, austr_sso, avstr_sso, avdis_sso,&
       vgust_con, vgust_dyn, sod_t, asod_t, tke, vabsmx_10m,         &
       dursun_r, dursun_m, sun_el, sun_azi, horizon,                 &
       hail_gsp, prh_gsp

USE data_modelconfig,   ONLY:                                        &
       ie, je, ke, ke1, istart, jstart, istartpar, iendpar,          &
       jstartpar, jendpar, dt, nehddt, dt2, ke_soil, idt_qv, idt_qc

USE data_runcontrol,    ONLY:                                        &
       lphys, ltur, llm, nold, nnow, nnew, ntstep, itype_synd, l2tls,&
       itype_gscp, hlastmxt, hnextmxt, hincmxt, nlastmxt, nnextmxt,  &
       lsso, itype_diag_gusts, ntke, nhori, lradtopo, lmulti_layer

USE data_parallel,      ONLY:                                        &
       my_cart_neigh, my_cart_id

USE data_constants,     ONLY:                                        &
       g, cp_d, cpdr, b1, b2w, b3, b4w, r_d, rdv, rvd_m_o,           &
       o_m_rdv, rdocp, lh_v, repsilon

USE turb_data,          ONLY:   akt

USE data_io,            ONLY:                                        &
       lbdclim

USE meteo_utilities,    ONLY :  calps

USE environment,        ONLY :  model_abort

USE src_tracer,         ONLY :  trcr_get, trcr_errorstr

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER, INTENT (IN)       :: nx
 
!------------------------------------------------------------------------------
!
! Local scalars:
REAL (KIND=wp)             ::    &
  zvmin, zumold, zvmold, zumold_t, zvmold_t, zrealdiff,                    &
  z10, zu10, zv10, zu30, zv30, z30, zumk, zvmk, zpil, zpiu, zthvl, zthvu,  &
  zthvd, zthv1, zbuon, zpeb, zitke, z_p, z_h, ztke, ztke_top, ztke_bot,    &
  zv_bra

REAL (KIND=wp)             ::    &
  z10g, z2g, zpke, zcm, zch, zh, zhs, zsqcm, zchdcm, zh1, zh2, zlnz1, zdz, &
  zh05m, zh2m, zp2m, zew2m, ze2m, zewke, zgqd2m, zgqdke, zqd2m, zfrac,     &
  zqd2mmin

REAL (KIND=wp)             ::    &
  zk, zk1, ztk1, zqk1, zpk1, zuk, zvk, zuk1, zvk1, zfk1,                   &
  zf2m, zf10m, z0m, zdzg, zval, zfac, zdtr, esat, zt, sec, ha_sun

INTEGER (KIND=iintegers)   ::     &
  zi_bra(ie,je)                     ! for the computation of Brasseur; values 0 or 1

INTEGER (KIND=iintegers)   ::    &
  i, j, k, ii, jj, ntl, kso, izerror

CHARACTER (LEN=255)        ::  yzerrmsg
CHARACTER (LEN=25)         ::  yzroutine

! Local arrays:
REAL (KIND=wp)             ::     &
  zhll (ie,je),                   & ! height of lowest model layer
  zhlt (ie,je),                   & ! height of the second lowest model layer
  zvpb (ie,je),                   & ! wind speed in the lowest model layer
  zvp10(ie,je),                   & ! wind speed in 10 m above ground
  zvp30(ie,je),                   & ! wind speed in 30 m above ground
  zum1 (ie,je),                   & ! u-wind on the lowest level at the mass position
  zvm1 (ie,je),                   & ! v-wind on the lowest level at the mass position
  zum2 (ie,je),                   & ! u-wind on level above the lowest at a mass pos.
  zvm2 (ie,je),                   & ! v-wind on level above the lowest at a mass pos.
  zhori(ie,je,nhori)

REAL (KIND=wp),     POINTER ::    &
  qv_new  (:,:,:) => NULL(),      & ! QV at nnew
  qv_now  (:,:,:) => NULL(),      & ! QV at nnow 
  qv_nx   (:,:,:) => NULL(),      & ! QV at nx
  qc_new  (:,:,:) => NULL(),      & ! QC at nnew 
  qc_now  (:,:,:) => NULL()         ! QC at nnow

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE near_surface
!------------------------------------------------------------------------------

! Statement functions
!--------------------

! Magnus formula for water
  esat (zt)       = b1 * EXP( b2w*(zt-b3)/(zt-b4w) )

  yzroutine = 'near_surface'
  izerror   = 0_iintegers
  yzerrmsg  = ''

!------------------------------------------------------------------------------
!  Section 0: Surface Pressure computation (from Relaxation)
!------------------------------------------------------------------------------

  ! retrieve the required microphysics tracers
  CALL trcr_get(izerror, idt_qv, ptr_tlev=nnew, ptr=qv_new)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qv, ptr_tlev=nnow, ptr=qv_now)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qv, ptr_tlev=nx, ptr=qv_nx)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev=nnew, ptr=qc_new)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev=nnow, ptr=qc_now)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  CALL calps ( ps(:,:,nnew), pp(:,:,ke,nnew), t(:,:,ke,nnew),     &
               qv_new(:,:,ke), qc_new(:,:,ke), qrs(:,:,ke),       &
               rho0(:,:,ke), p0(:,:,ke), dp0(:,:,ke),             &
               ie, je, rvd_m_o, r_d, 1, ie, 1, je)

  IF ( .NOT. l2tls ) THEN
    CALL calps ( ps(:,:,nnow), pp(:,:,ke,nnow), t(:,:,ke,nnow),     &
                 qv_now(:,:,ke), qc_now(:,:,ke), qrs(:,:,ke),       &
                 rho0(:,:,ke), p0(:,:,ke), dp0(:,:,ke),             &
                 ie, je, rvd_m_o, r_d, 1, ie, 1, je)
  ENDIF

  IF (l2tls) THEN
    zdtr = 1.0_wp/dt
    ntl  = nnow
  ELSE
    zdtr = 1.0_wp/dt2
    ntl  = nold
  ENDIF

  DO j = jstartpar, jendpar
    DO i = istartpar, iendpar
      dpsdt(i,j) = ( ps(i,j,nnew) - ps(i,j,ntl)) * zdtr
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
!  Section 1: Initializations
!------------------------------------------------------------------------------

  ! Heights of lowest two model layers
  zhll(:,:) = 0.5_wp*( hhl(:,:,ke)   + hhl(:,:,ke1) )
  zhlt(:,:) = 0.5_wp*( hhl(:,:,ke-1) + hhl(:,:,ke)  )

  ! minimal wind velocity
  zvmin = 0.01_wp

  ! Velocity of the wind in the Prandtl layer (k=ke) and in the layer
  ! above (k=ke-1) for time step nnow.
  ! Interpolate wind speed in 10 m above ground (for gust computation).
  ! And interpolation of u and v
  DO j = jstartpar,jendpar
    DO i = istartpar,iendpar
      ii = MAX(i-1,1)
      jj = MAX(j-1,1)

      zumold    = 0.5_wp * (u(i,j,ke,nnow) + u(ii,j,ke,nnow))
      zvmold    = 0.5_wp * (v(i,j,ke,nnow) + v(i,jj,ke,nnow))
      zvpb(i,j) = MAX ( SQRT (zumold**2 + zvmold**2) , zvmin )

      zumold_t  = 0.5_wp * (u(i,j,ke-1,nnow) + u(ii,j,ke-1,nnow))
      zvmold_t  = 0.5_wp * (v(i,j,ke-1,nnow) + v(i,jj,ke-1,nnow))

      IF     (itype_diag_gusts == 1 .OR. itype_diag_gusts == 4) THEN
        ! new computation of gusts
        z10  = hsurf(i,j) + 10.0_wp
        zu10 = ( zumold * (zhlt(i,j)-z10)                          &
                +zumold_t*(z10-zhll(i,j)) )/( zhlt(i,j)-zhll(i,j) )
        zv10 = ( zvmold * (zhlt(i,j)-z10)                          &
                +zvmold_t*(z10-zhll(i,j)) )/( zhlt(i,j)-zhll(i,j) )
        zvp10(i,j) = MAX ( SQRT (zu10**2 + zv10**2) , zvmin )
      ELSEIF (itype_diag_gusts == 2) THEN
        ! recompute wind at 30m level as in old model version
        ! note: we extrapolate from ke and ke-1 levels, which are
        !       currently roughly at 10m and 34.2m
        z30  = hsurf(i,j) + 30.0_wp
        zu30 = ( zumold * (zhlt(i,j)-z30)                          &
                +zumold_t*(z30-zhll(i,j)) )/( zhlt(i,j)-zhll(i,j) )
        zv30 = ( zvmold * (zhlt(i,j)-z30)                          &
                +zvmold_t*(z30-zhll(i,j)) )/( zhlt(i,j)-zhll(i,j) )
        zvp30(i,j) = MAX ( SQRT (zu30**2 + zv30**2) , zvmin )
      ENDIF

      zum1(i,j) = 0.5_wp * ( u(i,j,ke,nx) + u(ii,j,ke,nx) )
      zvm1(i,j) = 0.5_wp * ( v(i,j,ke,nx) + v(i,jj,ke,nx) )
      zum2(i,j) = 0.5_wp * ( u(i,j,ke-1,nx) + u(ii,j,ke-1,nx) )
      zvm2(i,j) = 0.5_wp * ( v(i,j,ke-1,nx) + v(i,jj,ke-1,nx) )

    ENDDO
  ENDDO

!------------------------------------------------------------------------------
!  Section 2: Reinitialize tmin_2m, tmax_2m and top_con, bas_con
!------------------------------------------------------------------------------

  IF (ntstep-1 == nnextmxt) THEN
    tmin_2m (:,:) = 500.0_wp
    tmax_2m (:,:) =   0.0_wp

    ! Determine next step for re-initializing
    hlastmxt = hnextmxt
    hnextmxt = hlastmxt + hincmxt
    nlastmxt = NINT (hlastmxt * 3600.0_wp / dt)
    nnextmxt = NINT (hnextmxt * 3600.0_wp / dt)
  ENDIF

  IF ( MOD ( ntstep-1, nehddt  ) == 0 ) THEN
    top_con (:,:) = 0.0_wp
    bas_con (:,:) = 0.0_wp
  ENDIF

!------------------------------------------------------------------------------
!  Section 3a: Compute wind, temperature and humidity at screen levels 
!              (10m and 2m) using a default scheme (itype_synopd = 1) based
!              on similarity theory
!------------------------------------------------------------------------------

  IF (.NOT.llm .AND. itype_synd == 1 ) THEN

    z10g = 10.0_wp * g
    z2g  =  2.0_wp * g

    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar

        ! some initializations
        ! check tcm and tch
        IF ( (lphys .EQV. .FALSE.) .OR. (ltur .EQV. .FALSE.) ) THEN
          zcm       = 1.0E-4_wp
          zch       = 1.0E-4_wp
          gz0 (i,j) = MAX ( gz0(i,j) , 1.0E-3_wp )
        ELSE
          zcm       = MAX ( tcm(i,j) , 5.0E-4_wp )
          zch       = MAX ( tch(i,j) , 4.0E-5_wp )
        ENDIF
        ! other variables
        zpke   = p0(i,j,ke) + pp(i,j,ke,nx)
        zh     = cp_d* t(i,j,ke,nx) + g*zhll(i,j)
        zhs    = cp_d* t_g(i,j,nx) + g*hsurf(i,j)
        zsqcm  = SQRT (zcm)
        zchdcm = zch / ( akt * zsqcm )
        zh1    = g*( zhll(i,j) - hsurf(i,j) )
        zh2    = g*( zhlt(i,j) - hsurf(i,j) )

        IF (zhs <= zh) THEN
          ! stable case
          zlnz1      = LOG ( (z10g + gz0(i,j)) / gz0(i,j) )  &
                     - z10g / zh1 * LOG ( (zh1+gz0(i,j)) / gz0(i,j) )
          zlnz1      =(z10g/zh1 + zsqcm/akt * zlnz1)
          zh05m      = zhs + 0.25_wp * (zh - zhs)
          zh2m       = zh05m + 1.5_wp * g * (zh - zh05m) / (zh1 - 0.5_wp*g)
        ELSE
          ! unstable case
          zsqcm     = MAX ( zsqcm  , 0.01_wp )
          zchdcm    = MAX ( zchdcm , 0.01_wp )

          zdz       = MAX(zh1-z10g, 0.0_wp)
          zlnz1     = 1.0_wp - zsqcm/akt * LOG ( 1.0_wp + (EXP ( akt/zsqcm ) - 1.0_wp) &
                    *gz0(i,j) * zdz / (zh1*(z10g + gz0(i,j))))
          zdz       = MAX(zh1 - z2g, 0.0_wp)
          zh2m      = zhs + (zh - zhs) * (1.0_wp - zchdcm *                            &
                      LOG ( 1.0_wp + (EXP ( 1.0_wp/zchdcm) - 1.0_wp) * (gz0(i,j) *     &
                            zdz / (zh1 * (z2g + gz0(i,j))) ) ))
        ENDIF

        ! wind in 10 m

        IF (z10g.LT.zh1) THEN !10m level below lowest model half level
           u_10m(i,j) = zum1(i,j) * zlnz1
           v_10m(i,j) = zvm1(i,j) * zlnz1
        ELSE !linear interpolation between the lowest two model levels
           u_10m(i,j) = (zum1(i,j)*(zh2 - z10g) + zum2(i,j)*(z10g - zh1))/ &
                        (zh2 - zh1)
           v_10m(i,j) = (zvm1(i,j)*(zh2 - z10g) + zvm2(i,j)*(z10g - zh1))/ &
                        (zh2 - zh1)
        ENDIF

        ! temperature in 2 m

        IF (z2g.LT.zh1) THEN !2m level below lowest model half level
           t_2m(i,j)  = ( zh2m - z2g - g*hsurf(i,j)) * cpdr
        ELSE !linear interpolation between the lowest two model levels
           t_2m(i,j)  = (t(i,j,ke,nx)*(zh2 - z2g) + t(i,j,ke-1,nx)*(z2g - zh1))/ &
                        (zh2 - zh1)
        ENDIF

        ! dew point and relative humidity in 2 m
        zp2m        = ps(i,j,nx) * (1.0_wp - z2g / &
                      (r_d * t_2m(i,j) * (1.0_wp + rvd_m_o* qv_nx(i,j,ke))))
        zew2m       = esat(t_2m(i,j))
        zewke       = esat(t(i,j,ke,nx))
        zgqd2m      = rdv * zew2m / (zp2m - o_m_rdv* zew2m)
        zgqdke      = rdv * zewke / (zpke - o_m_rdv* zewke)
        ! instead of using 1.0E-10 in the MAX, zqd2mmin is now used
        zqd2mmin    = 0.05_wp * zgqd2m
        zqd2m       = MAX ( zqd2mmin, qv_nx(i,j,ke) * zgqd2m/zgqdke )
        qv_2m(i,j)  = zqd2m

        ze2m        = zp2m*zqd2m / (rdv + o_m_rdv*zqd2m)
        rh_2m(i,j)  = 100.0_wp * MIN ( ze2m / zew2m , 1.0_wp )
        zfrac       = LOG ( ze2m / b1 )
        td_2m(i,j)  = (b2w*b3 - b4w*zfrac) / (b2w - zfrac)
        td_2m(i,j)  = MIN ( td_2m(i,j) , t_2m(i,j) )

      ENDDO
    ENDDO
  ENDIF

!------------------------------------------------------------------------------
!  Section 3b: Compute wind, temperature and humidity at screen levels
!              (10m and 2m) using a new scheme (itype_synd = 2) from
!              M. Raschendorfer at DWD which is also based on similarity 
!              theory. Geometrical height is counted with respect to the
!              lower boundary of the atmosphere, exluding the roughness
!              layer.
!              The very high vertical resulution case, where the lowest model
!              layer is below the 2m- or 10m-screen level, is also considered.
!------------------------------------------------------------------------------

  IF (.NOT.llm .AND. itype_synd == 2  &
               .AND. ( (lphys .EQV. .FALSE.) .OR. (ltur .EQV. .FALSE.) )) THEN

    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar

         IF (fr_land(i,j).LT.0.5_wp) THEN
            gz0(i,j) = MAX ( gz0(i,j) , 1.0E-3_wp )
         END IF
         z0m = gz0(i,j)/g

         ! Temperature and humidity at 2m screen-level

         k   = ke
         zk  = 0.5_wp * ( hhl(i,j,ke) - hsurf(i,j) )
         zk1 = 0.0_wp
         DO WHILE ( zk < 2.0_wp .AND. k > 1 )
           k   = k - 1
           zk1 = zk
           zk  = 0.5_wp*( hhl(i,j,k)+hhl(i,j,k+1) ) - hsurf(i,j)
         END DO
         zdzg = ( 2.0_wp - zk1 )*g
         IF ( k == ke ) THEN
           zfac = tfh (i,j)
           ztk1 = t_g (i,j,nx)
           zqk1 = qv_s(i,j,nx)
           zpk1 = ps  (i,j,nx)
         ELSE
           zfac = 1.0_wp
           ztk1 = t (i,j,k+1,nx)
           zqk1 = qv_nx(i,j,k+1)
           zpk1 = p0(i,j,k+1) + pp(i,j,k+1,nx)
         END IF
         zfk1 = LOG( (  zk + z0m    )/(zk1+z0m) )
         zf2m = LOG( (2.0_wp+z0m)/(zk1+z0m) ) / zfk1
         t_2m (i,j) = t (i,j,k,nx) + ( zf2m - 1.0_wp ) &
                                   * ( t (i,j,k,nx) - ztk1 )*zfac
         qv_2m(i,j) = qv_nx(i,j,k) + ( zf2m - 1.0_wp ) &
                                   * ( qv_nx(i,j,k) - zqk1 )*zfac
         zval = zpk1 * EXP( -zdzg / ( r_d*t_2m(i,j) &
                                   *(1.0_wp+rvd_m_o*qv_nx(i,j,k)) ) ) !P_2m
         zval = zval/(rdv/MAX(qv_2m(i,j),repsilon)-rdv+1.0_wp)        !e_2m
         rh_2m(i,j) = 100.0_wp * MIN ( zval/esat(t_2m(i,j)), 1.0_wp)
         zfac = LOG( zval/b1 )
         td_2m(i,j) = MIN ( ( b2w*b3 - b4w*zfac ) / ( b2w - zfac ) , t_2m(i,j) )

         ! Wind at 10m screen-level 

         DO WHILE ( zk < 10.0_wp .AND. k > 1 )
           k   = k - 1
           zk1 = zk
           zk  = 0.5_wp*( hhl(i,j,k)+hhl(i,j,k+1) ) - hsurf(i,j)  
         END DO
         ii  = MAX( i-1, 1 )
         jj  = MAX( j-1, 1 )
         zuk = 0.5_wp * ( u(i,j,k,nx) + u(ii,j,k,nx) )
         zvk = 0.5_wp * ( v(i,j,k,nx) + v(i,jj,k,nx) )
         IF ( k == ke ) THEN
           zfac = tfm(i,j)
           zuk1 = 0.0_wp
           zvk1 = 0.0_wp
         ELSE
           zfac = 1.0_wp
           zuk1 = 0.5_wp * ( u(i,j,k+1,nx) + u(ii,j,k+1,nx) )
           zvk1 = 0.5_wp * ( v(i,j,k+1,nx) + v(i,jj,k+1,nx) )
         END IF
         zfk1       = LOG( (zk+z0m) / (zk1+z0m) )
         zf10m      = LOG( (10.0_wp+z0m) / (zk1+z0m) ) / zfk1
         u_10m(i,j) = zuk + ( zf10m - 1.0_wp )*( zuk - zuk1 )*zfac
         v_10m(i,j) = zvk + ( zf10m - 1.0_wp )*( zvk - zvk1 )*zfac

       END DO
    END DO
  END IF

!------------------------------------------------------------------------------
!  Section 3c: Compute wind, temperature and humidity at screen levels
!              llm-runs
!------------------------------------------------------------------------------

  IF ( llm ) THEN
    ! in case of an llm_run,
    ! set u_10m, v_10m, vmax_10m to unrealistic values, td_2m, rh_2m to 0.0 and
    ! t_2m, qv_2m to the lowest level values of t and qv
    u_10m     (:,:) = 99.9_wp
    v_10m     (:,:) = 99.9_wp
    vmax_10m  (:,:) = 99.9_wp
    vabsmx_10m(:,:) = 99.9_wp
    vgust_dyn (:,:) = 99.9_wp
    vgust_con (:,:) = 99.9_wp
    t_2m      (:,:) = t (:,:,ke,nx)
    qv_2m     (:,:) = qv_nx(:,:,ke)
    td_2m     (:,:) = 0.0_wp
    rh_2m     (:,:) = 0.0_wp
  ELSE
    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar
        ! maximum wind speed in 10m without gust
        vabsmx_10m(i,j) = MAX (vabsmx_10m(i,j) ,                         &
                               SQRT(u_10m(i,j)*u_10m(i,j)+v_10m(i,j)*v_10m(i,j)) )
      ENDDO
    ENDDO
  ENDIF

!------------------------------------------------------------------------------
!  Section 4: 2m temperature, minima, maxima
!------------------------------------------------------------------------------

  DO j = jstartpar, jendpar
    DO i = istartpar, iendpar
      ! minimal and maximal temperature in 2m
      tmin_2m(i,j) = MIN ( tmin_2m(i,j) , t_2m(i,j) )
      tmax_2m(i,j) = MAX ( tmax_2m(i,j) , t_2m(i,j) )
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
!  Section 5: maximal expected squall
!------------------------------------------------------------------------------

  ! The maximum wind speed at 10 m includes an empirical formula for
  ! the gustiness based on the turbulent kinetic energy in the Prandtl
  ! layer
  !   zvmax_10m = SQRT ( u(i,j,ke,t)**2 + v(i,j,ke,t)**2 )*    &
  !                    (1._wp + 3._wp*2.4_wp*zsqcm)
  !
  ! For higher model levels, a speed in 30 m is estimated and:
  ! zvmax_10m = speed_30m + 3*2.4*ustar = speed_30m + 3*2.4*zsqcm*speed_ke
  ! This was used before and is now calculated for option itype_diag_gusts = 2
  ! 
  ! A re-tuning of the gusts in the COSMO-EU model domain (with a lowest 
  ! model level at about 10 m) resulted in the replacement of zvp30 
  ! (wind speed at 30 m) and zvpb (wind speed in the lowest model layer) 
  ! by zvp10 (wind speed at 10 m). This is computed for itype_diag_gusts = 1
  !
  ! itype_diag_gusts = 3 computes gusts after the method of Brasseur
  !
  ! itype_diag_gusts = 4 computes gusts like itype_diag_gusts = 1, but the
  ! gust factor (3*2.4) becomes now 3*2.4 + 0.09*zvp10, which makes it weakly
  ! dependent on the mean wind speed at 10 m.

  IF     (itype_diag_gusts == 1) THEN

    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar
        IF ( (.NOT. lphys) .OR. (.NOT. ltur) ) THEN
          zcm = 1.0E-4_wp
        ELSE
          zcm = MAX ( tcm(i,j) , 5.0E-4_wp )
        ENDIF
        zsqcm = SQRT (zcm)
        vgust_dyn(i,j) = MAX ( vgust_dyn(i,j) ,                             &
                zvp10(i,j) + 3.0_wp * 2.4_wp * zsqcm * zvp10(i,j) )
      ENDDO
    ENDDO

  ELSEIF (itype_diag_gusts == 2) THEN

    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar
        IF ( (.NOT. lphys) .OR. (.NOT. ltur) ) THEN
          zcm = 1.0E-4_wp
        ELSE
          zcm = MAX ( tcm(i,j) , 5.0E-4_wp )
        ENDIF
        zsqcm = SQRT (zcm)
        vgust_dyn(i,j) = MAX ( vgust_dyn(i,j) ,                             &
                zvp30(i,j) + 3.0_wp * 2.4_wp * zsqcm * zvpb(i,j) )
      ENDDO
    ENDDO

  ELSEIF (itype_diag_gusts == 3) THEN
    ! Compute wind gusts after Brasseur: but in this form it will not vectorize!!!

    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar
  
        zv_bra = SQRT( zum1(i,j)**2 + zvm1(i,j)**2 )
  
        ! tke in the code is sqrt(2*turb. kin. En.)
  
        ztke_bot  = 0.5_wp * tke(i,j,ke+1,ntke)**2
        ztke_top  = 0.01_wp * ztke_bot
        zpeb  = 0.0_wp
        zitke = 0.0_wp
        z_p   = 0.0_wp
  
        ! init loop exit condition
        zi_bra(:,:) = 0
  
        DO k=ke,2,-1      !! Loop over levels
  
          IF (zi_bra(i,j) == 0 ) THEN
  
            ! compute virtual potential temperature zthv (lower and upper)
            zpil   = (1.0E-5_wp * (p0(i,j,k  ) + pp(i,j,k  ,nx)))**rdocp
            zpiu   = (1.0E-5_wp * (p0(i,j,k-1) + pp(i,j,k-1,nx)))**rdocp
            zthvl  = t(i,j,k  ,nx)*(1.0_wp+rvd_m_o*qv_nx(i,j,k  ))/zpil
            zthvu  = t(i,j,k-1,nx)*(1.0_wp+rvd_m_o*qv_nx(i,j,k-1))/zpiu
  
            ! compute vertical integral of potential energy of buoyancy (zpeb)
            z_h    = hhl(i,j,k) - hhl(i,j,k+1)
            !zthvd  = abs(zthvu - zthvl)
            zthvd  = zthvu - zthvl
            zthv1  = 0.5_wp * (zthvl + zthvu)
            zbuon  = g*zthvd/zthv1
            zpeb   = zpeb  + zbuon * z_h
  
            ! compute vertical integral of tke (zitke)
            z_p    = z_p + z_h
            !ztke   = (0.5_wp * tke(i,j,k,ntke)**2+0.5_wp * tke(i,j,k-1,ntke)**2)*0.5_wp
            ztke   = 0.5_wp * tke(i,j,k,ntke)**2
            zitke  = zitke + ztke * z_h
  
            ! compute wind
            zum1(i,j) = 0.5_wp * ( u(i,j,k,nx) + u(MAX(1,i-1),j,k,nx) )
            zvm1(i,j) = 0.5_wp * ( v(i,j,k,nx) + v(i,MAX(1,j-1),k,nx) )
            zum2(i,j) = 0.5_wp * ( u(i,j,k-1,nx) + u(MAX(1,i-1),j,k-1,nx) )
            zvm2(i,j) = 0.5_wp * ( v(i,j,k-1,nx) + v(i,MAX(1,j-1),k-1,nx) )
            zumk   = 0.5_wp * (zum1(i,j) + zum2(i,j))
            zvmk   = 0.5_wp * (zvm1(i,j) + zvm2(i,j))
            zv_bra = SQRT( zumk**2 + zvmk**2 )
  
            ! loop exit condition
            IF (ztke <= ztke_top) zi_bra(i,j) = 1   ! Loop exit
            IF (zitke/z_p <= zpeb ) zi_bra(i,j) = 1 ! Loop exit
  
          ENDIF
  
        ENDDO  !! Loop over levels
  
        vgust_dyn(i,j) = MAX (vgust_dyn(i,j), zv_bra)
  
      ENDDO
    ENDDO

  ELSEIF (itype_diag_gusts == 4) THEN

    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar
        IF ( (.NOT. lphys) .OR. (.NOT. ltur) ) THEN
          zcm = 1.0E-4_wp
        ELSE
          zcm = MAX ( tcm(i,j) , 5.0E-4_wp )
        ENDIF
        zsqcm = SQRT (zcm)
        vgust_dyn(i,j) = MAX ( vgust_dyn(i,j) ,                             &
                zvp10(i,j)                                                  &
              + (3.0_wp * 2.4_wp + 0.09_wp * zvp10(i,j))        & ! gust factor
              * zsqcm * zvp10(i,j) )
      ENDDO
    ENDDO

  ENDIF

  ! compute combination of dynamical and convective gust 
  ! (old meaning of vmax_10m)
  vmax_10m(:,:) = MAX (vgust_dyn(:,:), vgust_con(:,:))

!------------------------------------------------------------------------------
!  Section 4b: summations, means
!------------------------------------------------------------------------------

  ! for the leapfrog scheme, these summations must not be done in the first
  ! step, because it is just an intermediate step. These values are calculated
  ! again in ntstep == 1.

  IF (l2tls .OR. (ntstep > 0)) THEN

    IF (lradtopo) THEN
      zhori(:,:,:) = horizon(:,:,:)
    ELSE
      zhori(:,:,:) = 0.0_wp
    ENDIF

    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar

        ! accumulation of precipitation rates
        rain_gsp(i,j) = rain_gsp(i,j) + dt * prr_gsp(i,j)
        snow_gsp(i,j) = snow_gsp(i,j) + dt * prs_gsp(i,j)
        IF (itype_gscp >= 4) THEN
          grau_gsp(i,j) = grau_gsp(i,j) + dt * prg_gsp(i,j)
        ENDIF
        ! Possible only when running the 2-moment microphysics
        IF (itype_gscp >= 2000) THEN
          hail_gsp(i,j) = hail_gsp(i,j) + dt * prh_gsp(i,j)
        ENDIF

        rain_con(i,j) = rain_con(i,j) + dt * prr_con(i,j)
        snow_con(i,j) = snow_con(i,j) + dt * prs_con(i,j)

        IF (lbdclim) THEN
          ! time mean values of  t_2m, td_2m, u_10m, v_10m
          t_2m_av (i,j)    = t_2m_av (i,j) + t_2m (i,j)
          td_2m_av(i,j)    = td_2m_av(i,j) + td_2m(i,j)
          u_10m_av(i,j)    = u_10m_av(i,j) + u_10m(i,j)
          v_10m_av(i,j)    = v_10m_av(i,j) + v_10m(i,j)
        ENDIF

        ! mean values of solar and thermal radiative heating
        asob_t(i,j) = asob_t(i,j) + sobt(i,j)
        athb_t(i,j) = athb_t(i,j) + thbt(i,j)
        asob_s(i,j) = asob_s(i,j) + sobs(i,j)
        athb_s(i,j) = athb_s(i,j) + thbs(i,j)
        apab_s(i,j) = apab_s(i,j) + pabs(i,j)
        asod_t(i,j) = asod_t(i,j) + sod_t(i,j)

        ! from Burkhardt Rockel, GKSS
        IF (sodwddm(i,j) > 120.0_wp) THEN
          dursun(i,j) = dursun(i,j) + dt
        ENDIF

        ! maximum possible sunshine duration buz
        sec=(360.0_wp/nhori)
        ii = INT(sun_azi(i,j)/sec)
        IF (ii >= nhori) THEN
           ii = nhori - 1
        ENDIF
        k  = MOD(ii+1,24)
        sec=(360.0_wp/nhori)
        ! interpolation of horizon in the direction of sun azimuth
        ! ( 1.0E-7_wp=securi due to zsmu0 > 1.0E-9_wp)
        ha_sun = 1.0E-7_wp + (zhori(i,j,k+1)*(sun_azi(i,j)-sec*ii)+ &
                 zhori(i,j,ii+1)*(sec*(ii+1)-sun_azi(i,j)))/sec
        IF (sun_el(i,j) > ha_sun) THEN
           dursun_m(i,j) = dursun_m(i,j) + dt
        ENDIF

        !relative sunshine duration buz
        IF (dursun_m(i,j) > 0.0_wp) THEN
            dursun_r(i,j) = 100.0_wp*dursun(i,j)/dursun_m(i,j)
        ENDIF

        ! accumulation of surface fluxes
        aumfl_s(i,j) = aumfl_s(i,j) + umfl_s(i,j)
        avmfl_s(i,j) = avmfl_s(i,j) + vmfl_s(i,j)
        ashfl_s(i,j) = ashfl_s(i,j) + shfl_s(i,j)
        alhfl_s(i,j) = alhfl_s(i,j) + lhfl_s(i,j)
        alhfl_bs(i,j)= alhfl_bs(i,j) + lhfl_bs(i,j)

      ENDDO
    ENDDO

    ! Accumulation of fluxes in the soil
    IF (lmulti_layer) THEN
      DO kso=1,ke_soil
        DO j = jstartpar, jendpar
          DO i = istartpar, iendpar
            alhfl_pl(i,j,kso) = alhfl_pl(i,j,kso) + lhfl_pl(i,j,kso)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF (lsso) THEN
      DO j = jstartpar, jendpar
        DO i = istartpar, iendpar

          ! accumulation of fluxes due to sub-grid scale orography
          austr_sso(i,j) = austr_sso(i,j) + ustr_sso(i,j)
          avstr_sso(i,j) = avstr_sso(i,j) + vstr_sso(i,j)
          avdis_sso(i,j) = avdis_sso(i,j) + vdis_sso(i,j)

        ENDDO
      ENDDO
    ENDIF

    ! accumulation of single solar (short wave) and thermal (long wave)
    ! surface components
    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar
        aswdir_s(i,j)  = aswdir_s(i,j)  + swdir_s(i,j)
        aswdifd_s(i,j) = aswdifd_s(i,j) + swdifd_s(i,j)
        aswdifu_s(i,j) = aswdifu_s(i,j) + swdifu_s(i,j)
        alwd_s(i,j)    = alwd_s(i,j)    + lwd_s(i,j)
        alwu_s(i,j)    = alwu_s(i,j)    + lwu_s(i,j)
      ENDDO
    ENDDO
  ENDIF

!------------------------------------------------------------------------------
!- End Subroutine near_surface
!------------------------------------------------------------------------------

END SUBROUTINE near_surface
