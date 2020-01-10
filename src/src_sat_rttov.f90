!+ Source Module for generation of synthetic satellite images using RTTOV9
!------------------------------------------------------------------------------

MODULE src_sat_rttov

!------------------------------------------------------------------------------
!
! Description:
!   1- Setup up profiles and call rttov_ifc_fill
!   2- Run rttov9.3 with interpolation and scattering parameterization
!
!
! Method:
!   see comments in program
!       
! Current Code Owner: DWD, Robin Faulwetter
!  phone:  +49  69  8062 2746
!  fax:    +49  69  8062 3721
!  email:  Robin.Faulwetter@dwd.de
!       
! History:
! Version    Date       Name
! ---------- ---------- ----
! V4_18        2011/05/26 Robin Faulwetter
!  Initial Release
! V4_24        2012/06/22 Robin Faulwetter
!  Added a NOMOVE Directive for the SX-compiler
! V4_25        2012/09/28 Ulrich Schaettler
!  Check soiltyp with NINT-function
! V4_26        2012/12/06 Andreas Messer
!  Modifications for using also RTTOV10
! V4_27        2013/03/19 Ulrich Schaettler
!  Technical changes in the use of src_obs_rad
! V4_28        2013/07/12 Ulrich Schaettler
!  Renamed nlist_chan to nchan_list
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
! V5_2         2015-05-21 Robin Faulwetter
!  Improved calculation of input profiles for RTTOV
! V5_3         2015-10-09 Robin Faulwetter
!  Corrected interpolation of input for RTTOV (linterp=T)
!
! Code Description:
!   Language:          Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:
!
!==============================================================================
#if defined(RTTOV9) || defined(RTTOV10)

USE data_constants  , ONLY :  &
    r_earth, b1, b2w, b3, b4w, rdv, o_m_rdv

USE data_fields,      ONLY :  &
    clw_con, p0, rho, clc_sgs, clc_con,         &
    t_2m, qv_2m, u_10m, v_10m,                  &
    fr_land, hsurf, soiltyp, lseamask,          &
    rlat, rlon, sun_el, sun_azi

USE data_modelconfig, ONLY :  &
    raddeg

USE data_parallel,    ONLY :  &
    my_cart_id,               &
    my_world_id

USE data_parameters , ONLY :  &
    wp,        &! KIND-type parameters for real variables
    iintegers   ! kind-type parameter for "normal" integer variables

USE data_runcontrol,  ONLY :  &
    lseaice, ldebug_dia, lprintdeb_all, idbg_level, ntstep

USE data_satellites , ONLY :  &
    sat_compute, num_sensors, numchans, rcnw,    &
    extrp_type, iceshape, iwc2effdiam,           &
    extrp_const, extrp_lin, extrp_clim,          &
    p_top, t_top, q_top, lcon_clw, zenmax9,      &
    zenmax10

#if defined(RTTOV9)
USE mo_rttov_ifc,     ONLY :                     &
    rttov_fill_input, rttov_direct_ifc,          &
    NO_ERROR, rttov_ifc_errMsg

#elif defined(RTTOV10)
USE mo_rttov_ifc,     ONLY :                     &
    rttov_fill_input, rttov_direct_ifc,          &
    NO_ERROR, rttov_ifc_errMsg,                  &
    rttov_ifc_version, OUT_ASB, OUT_CSB, OUT_ASR,&
    OUT_CSR
#endif

USE src_obs_rad,      ONLY:   &
    prepare_rttov_input,      &
    get_n_add

!==============================================================================

IMPLICIT NONE

!==============================================================================

CONTAINS

!==============================================================================
! NEW VERSION USING RTTOV9.3 WITH RTTOV_IFC
!==============================================================================

SUBROUTINE org_sat_rttov (t, qv, qc, qi, qs, qg, pp, ps, t_g, h_ice,       &
                          idim, jdim, ke, idims, idime, jdims, jdime,      &
                          nsensors, synme7, synmsg, nusyn, yerror, ierror)

!------------------------------------------------------------------------------

INTEGER (KIND=iintegers), INTENT(IN) ::           &
  idim, jdim, ke,             & ! dimensions of fields
  idims, idime, jdims, jdime, & ! start and end-indices for the loops
  nsensors,                   & ! Max. No. of sensors to be used
  nusyn

INTEGER (KIND=iintegers), INTENT(OUT) ::          &
  ierror                        ! error status variable

CHARACTER (LEN= *),       INTENT(OUT) ::          &
  yerror                        ! error message

REAL  (KIND=wp),          INTENT(IN) ::           &
  t      (idim,jdim,ke),      & !
  qv     (idim,jdim,ke),      & !
  qc     (idim,jdim,ke),      & !
  qi     (idim,jdim,ke),      & !
  qs     (idim,jdim,ke),      & !
  pp     (idim,jdim,ke),      & !
  ps     (idim,jdim),         & !
  t_g    (idim,jdim)            ! prognostic variables


REAL  (KIND=wp),          INTENT(IN) ::           &
  ! variables, that are possibly not available
  h_ice  (:,:),               & !
  qg     (:,:,:)

REAL  (KIND=wp),     INTENT(OUT), TARGET ::       &
  synme7 (idim,jdim, 8),      & !
  synmsg (idim,jdim,32)         ! synthetic satellite images

! Local variables:
!------------------

INTEGER (KIND=iintegers)   ::        &
  izdebug                              ! for debug messages

! Atmopheric profiles
REAL  (KIND=wp),     ALLOCATABLE ::  &
  pres    (:,:),              & !
  humidity(:,:),              & !
  temp    (:,:),              & !
  clc  (:,:,:),               & !
  cloud(:,:,:)                  !

! Other FOV (field of view) related variables
REAL  (KIND=wp)      ::       &
  psurf(idims:idime),         & !
  t2m  (idims:idime),         & !
  hum2m(idims:idime),         & !
  u10m (idims:idime),         & !
  v10m (idims:idime),         & !
  t_s  (idims:idime),         & !
  s_hgt(idims:idime),         & !
  sat_a(idims:idime),         & !
  sat_z(idims:idime)            !

INTEGER (KIND=iintegers)  ::  &
  stype(idims:idime),         & !
  wtype(idims:idime),         & !
  idg(idims:idime),           & !
  ish(idims:idime)

! Other local variables
REAL  (KIND=wp),     POINTER ::      &
  syn(:,:,:)

REAL  (KIND=wp)      :: lon

INTEGER (KIND=iintegers)   ::        &
  ncalc   (nsensors),             &
  iprof(nsensors, (idime - idims + 1)*maxval(numchans(:))), &
  ichan(nsensors, (idime - idims + 1)*maxval(numchans(:)))  

REAL  (KIND=wp)      ::              &
  t_b      (maxval(numchans(:)), (idime - idims + 1)), &
  t_b_clear(maxval(numchans(:)), (idime - idims + 1)), &
  rad      (maxval(numchans(:)), (idime - idims + 1)), &
  rad_clear(maxval(numchans(:)), (idime - idims + 1)), &
  emiss    (maxval(numchans(:)), (idime - idims + 1))


INTEGER (KIND=iintegers) ::          &
  i, j, k, l, n, nzp, isens,         &
  n_profs, im, jm, istatus,          &
  rad_output

INTEGER (KIND=iintegers), SAVE :: n_lev = 0
INTEGER (KIND=iintegers), SAVE :: n_add = 0

INTEGER (KIND=iintegers) ::          &
  iprint = 0_iintegers

LOGICAL ::                           &
  l_h_ice, l_qg

REAL (KIND=wp)     ::                &
  alpha_e, r_atm, r_sat

REAL (KIND=wp) :: T_b_m = 0._wp
REAL (KIND=wp) :: fac   = 0._wp

CHARACTER (LEN=80) :: yzerrmsg

! Explanantion of MSG SEVIRI channels
! number lambda(micrometers)  Useful for                                  Peak
! 1      IR 3.9               Surface, clouds, wind fields                surface
! 2      WV 6.2               Water vapor, high level clouds,             300 hPa
! 3      WV 7.3               Water vapor                                 500 hPa
! 4      IR 8.7               Surface, clouds                             surface
! 5      IR 9.7               Ozone                                       surface + 50hPa
! 6      IR 10.8              Surface, clouds, wind fields                surface
! 7      IR 12.0              Surface, clouds, wind fields                surface
! 8      IR 13.4              Cirrus cloud height                         850 hPa

!-----End of header------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initialisations
!------------------------------------------------------------------------------

  ! Initialize, whether debug output shall be done
  IF (ldebug_dia) THEN
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

  !----------------------------------------------------------------------------
  ! Section 1.1: Initialize some variables
  !----------------------------------------------------------------------------

  IF (izdebug > 9) THEN
    PRINT *, '      RTTOV 9 initializations'
  ENDIF

  ! Check avalability of variables
  l_h_ice = ALL(SHAPE(h_ice) > 0)
  l_qg    = ALL(SHAPE(qg) > 0)

  ! distance of satellite to middle of the earth
  r_sat       = 35880E3_wp + r_earth

  ierror     = 0_iintegers
  yerror     = '         '

  !----------------------------------------------------------------------------
  ! Section 1.2: Set up instrument, channel and profile numbers for rttov
  !----------------------------------------------------------------------------

  n_profs = (idime - idims + 1)
  DO isens = 1, num_sensors
    ncalc(isens)    = n_profs * numchans(isens)
    DO j = 1, numchans(isens)
      iprof(isens, j:ncalc(isens):numchans(isens)) = (/ (k, k=1,n_profs) /)
      ichan(isens, j:ncalc(isens):numchans(isens)) = &
           (/ (sat_compute(isens)%nchan_list(j), k=1,n_profs) /)
    ENDDO
  ENDDO

  !----------------------------------------------------------------------------
  ! Section 1.3: Set up instrument, channel and profile numbers for rttov
  !----------------------------------------------------------------------------
  IF (n_lev <= 0) THEN
    call get_n_add(n_add)
    n_lev = ke + 1 + n_add  ! ke levels + 1 top level + additional levels
                            ! for optional extrapolation above model top
    IF (my_world_id == 0) THEN
      WRITE (nusyn,*) '  '
      WRITE (nusyn,'(A,T44,I4)') 'Extrapolation type         :', extrp_type
      WRITE (nusyn,'(A,T44,I4)') 'Number of additional levels:', n_add
    END IF
  END IF
  ALLOCATE(pres    (n_lev, idims:idime))
  ALLOCATE(humidity(n_lev, idims:idime))
  ALLOCATE(temp    (n_lev, idims:idime))
  ALLOCATE(clc  (6, n_lev-1, idims:idime)) !RF: NEW
  ALLOCATE(cloud(6, n_lev-1, idims:idime)) !RF: NEW

!------------------------------------------------------------------------------
! Section 2: Compute the synthetic brightness temperatures
!------------------------------------------------------------------------------

  IF (izdebug > 9) THEN
    PRINT *, '      RTTOV 9 j-loop to compute brightness temperatures'
  ENDIF

!CDIR NOMOVE
  jloop:  DO j = jdims, jdime

    !--------------------------------------------------------------------------
    ! Section 2.1: Generate RTTOV fields and boundary data from current state
    !--------------------------------------------------------------------------

    CALL prepare_rttov_input((/(i, i=idims,idime)/),             &
                             (/(j, i=idims,idime)/),             &
                             n_lev,                              &
                             n_add,                              &
                             pp,ps,t,t_g,qv,qc,qi,qs,qg,h_ice,   &
                             pres,temp,humidity,                 &
                             t2m,hum2m,psurf,s_hgt,u10m,v10m,t_s,&
                             stype,wtype,                        &
                             cloud    = cloud,                   &
                             cfrac    = clc,                     &
                             idg      = idg,                     &
                             ish      = ish,                     &
                             ierrstat = istatus)

    IF (istatus /= 0) THEN
      WRITE(yerror,'(" ERROR *** while computing synthetic satellite &
           &images: prepare_rttov_input ",I3)') istatus
      ierror = 9005
      RETURN
    ENDIF

    !--------------------------------------------------------------------------
    ! Section 2.2: Supply the arrays to the RTTOV library
    !--------------------------------------------------------------------------

    istatus = rttov_fill_input(                          &
          press      = pres    (1:n_lev, idims:idime),   &
          temp       = temp    (1:n_lev, idims:idime),   &
          humi       = humidity(1:n_lev, idims:idime),   &
          t2m        = t2m  (idims:idime),               &
          q2m        = hum2m(idims:idime),               &
          psurf      = psurf(idims:idime),               &
          hsurf      = s_hgt(idims:idime),               &
          u10m       = u10m (idims:idime),               &
          v10m       = v10m (idims:idime),               &
          stemp      = t_s  (idims:idime),               &
          stype      = stype(idims:idime),               &
          latgroundp = rlat(idims:idime,j) * raddeg,     &
          satzenith  = (/(0.0_wp, i=idims,idime)/),      &
          sunZenith  = sun_el(idims:idime,j),            &
          cloud      = cloud(1:6, 1:n_lev-1, idims:idime), & !RF: NEW
          cfrac      = clc(1:6,   1:n_lev-1, idims:idime), & !RF: NEW
          idg        = idg(idims:idime),                 &
          ish        = ish(idims:idime),                 &
          watertype  = wtype(idims:idime),               &
#ifdef RTTOV10
          addsolar   = .false.,                          &
          addrefrac  = .true.,                           &
          rttov9_compat = .true.,                        &
#else
          addsolar   = (/(.false., i=idims,idime)/),     &
          addrefrac  = (/(.true., i=idims,idime)/),      &
#endif
          addinterp  = .true.,                           &
          ivect      = 1)

    IF (istatus /= NO_ERROR) THEN
      WRITE(*,'("*** RTTOV ERROR (rttov_fill_input, proc ",I4,") ",A)') &
           my_cart_id, TRIM(rttov_ifc_errMsg(istatus))
      WRITE(yerror,'(" ERROR *** while computing synthetic satellite &
           &images: rttov_fill_input ",I3)') istatus
      ierror = 9005
      RETURN
    ENDIF

    sensor_loop: DO isens = 1, nsensors

      !------------------------------------------------------------------------
      ! Section 2.3: Set/compute some sensor dependent quantities
      !------------------------------------------------------------------------

!CDIR nodep
!CDIR vector             
      DO i = idims, idime  
        ! Since the emissitivy is intent(inout) in RTTOV, we have to 
        ! reinitialize it
        DO k = 1,  numchans(isens)
          emiss(k, i-idims+1) = sat_compute(isens)%emissivity(k)
        ENDDO
        lon = rlon(i,j) - sat_compute(isens)%longitude
        ! Calculate the satellite zenith angle 
        alpha_e  = ACOS(COS(rlat(i,j)) * COS(lon))
        r_atm    = SQRT(r_sat**2 +r_earth**2 -2*r_sat*r_earth*COS(alpha_e))
        sat_z(i) = ASIN(SIN(alpha_e)*r_sat/r_atm) * raddeg
#ifdef RTTOV10
        IF (rttov_ifc_version < 10) THEN
#else
        IF (.TRUE.) THEN
#endif
          sat_z(i) = MIN(ABS(sat_z(i)), zenmax9)
        ELSE
          sat_z(i) = MIN(ABS(sat_z(i)), zenmax10)
        ENDIF

        ! Calculate the satellite azimuth angle
        sat_a(i) = raddeg * (1.0_wp + ATAN2(TAN(lon), SIN(rlat(i,j))))
      ENDDO

      rad_output = 0
      IF (sat_compute(isens)%lcloud_tem) rad_output = ibset(rad_output, OUT_ASB)
      IF (sat_compute(isens)%lclear_tem) rad_output = ibset(rad_output, OUT_CSB)
      IF (sat_compute(isens)%lcloud_rad) rad_output = ibset(rad_output, OUT_ASR)
      IF (sat_compute(isens)%lclear_rad) rad_output = ibset(rad_output, OUT_CSR)

      !------------------------------------------------------------------------
      ! Section 2.4: Compute the brightness temperatures
      !------------------------------------------------------------------------

      istatus = rttov_direct_ifc(                                 &
             isens,                                               &
             iprof(isens,1:ncalc(isens)),                         &
             ichan(isens,1:ncalc(isens)),                         &
             emiss(1:numchans(isens), 1:n_profs),                 &
             satAzim   = sat_a(idims:idime),                      &
             satZenith = sat_z(idims:idime),                      &
             T_b       = T_b(1:numchans(isens), 1:n_profs),       &
             T_b_clear = T_b_clear(1:numchans(isens), 1:n_profs), &
             rad       = rad      (1:numchans(isens), 1:n_profs), &
             radClear  = rad_clear(1:numchans(isens), 1:n_profs), &
             rad_out_flg = rad_output                             , &
             iprint    = iprint)

      IF (istatus /= NO_ERROR) THEN
        CALL print_profile(idims,j)
        WRITE(*,*) '*** RTTOV ERROR (rttov_direct_ifc, proc ', &
             my_cart_id,') ',TRIM(rttov_ifc_errMsg(istatus))
        WRITE(yerror,'(" ERROR *** while computing synthetic satellite &
             &images: rttov_fill_input ",I3)') istatus
        ierror = 9006
        RETURN
      ENDIF

      IF ((sat_compute(isens)%ysatellite == 'METEOSAT').and.&
          (sat_compute(isens)%nsat_id == 7)) THEN
        syn => synme7
      ELSEIF ((sat_compute(isens)%ysatellite == 'MSG').and.&
              (sat_compute(isens)%nsat_id == 2)) THEN
        syn => synmsg
      ELSE
        PRINT *,'*** ERROR in src_sat_tbs: satellite '//&
                TRIM(sat_compute(isens)%ysatellite), &
                sat_compute(isens)%nsat_id,' not implemented.'
        CYCLE
      ENDIF

      IF (sat_compute(isens)%lcloud_tem) THEN
        DO k = 1, numchans(isens)
          syn(idims:idime,j, (k-1)*4+1) = T_b(k, 1:n_profs) 
          IF (ANY(T_b(k, 1:n_profs) > 1.E30_wp)) PRINT*,my_cart_id,'invalid T_b',T_b(k, 1:n_profs)
        ENDDO
      ENDIF
      IF (sat_compute(isens)%lclear_tem) THEN
        DO k = 1, numchans(isens)
          syn(idims:idime,j, (k-1)*4+2) = T_b_clear(k, 1:n_profs) 
          IF (ANY(T_b_clear(k, 1:n_profs) > 1.E30_wp)) PRINT*,my_cart_id,'invalid T_b_clear',T_b_clear(k, 1:n_profs)
        ENDDO
      ENDIF
      IF (sat_compute(isens)%lcloud_rad) THEN
        DO k = 1, numchans(isens)
          syn(idims:idime,j, (k-1)*4+3) = Rad(k, 1:n_profs) 
          IF (ANY(Rad(k, 1:n_profs) > 1.E30_wp)) PRINT*,my_cart_id,'invalid Rad',Rad(k, 1:n_profs)
        ENDDO
      ENDIF
      IF (sat_compute(isens)%lclear_rad) THEN
        DO k = 1, numchans(isens)
          syn(idims:idime,j, (k-1)*4+4) = Rad_clear(k, 1:n_profs) 
          IF (ANY(Rad_clear(k, 1:n_profs) > 1.E30_wp)) PRINT*,my_cart_id,'invalid Rad_clear',Rad_clear(k, 1:n_profs)
        ENDDO
      ENDIF

    END DO sensor_loop
  END DO jloop

  IF (izdebug > 9) THEN
    PRINT *, '      RTTOV 9 end of j-loop to compute brightness temperatures'
  ENDIF

!==============================================================================

CONTAINS

!==============================================================================

!------------------------------------------------------------------------------
!+ Internal procedure in "org_sat_tbs" for preventing supersaturation
!------------------------------------------------------------------------------

ELEMENTAL REAL(KIND=wp) FUNCTION q_sat(p, t)
  REAL(KIND=wp), intent(in) :: p, t

  REAL(KIND=wp) :: p_sat

  p_sat = B1 * exp(B2w * (t - B3) / (t - B4w))
  q_sat = Rdv * p_sat / max((p - O_m_rdv * p_sat), 1.0_wp)

END FUNCTION q_sat

!==============================================================================

!------------------------------------------------------------------------------
!+ Internal procedure in "org_sat_tbs" for debugging purposes
!------------------------------------------------------------------------------

SUBROUTINE print_profile(i,j)

  INTEGER (KIND=iintegers), INTENT(IN) :: i,j
  INTEGER  (KIND=iintegers):: k,l

  DO k = 1, size(pres,1)
    PRINT *, my_cart_id, 'pres', k, pres(k,i)
  ENDDO
  DO k = 1, size(temp,1)
    PRINT *, my_cart_id, 'temp', k, temp(k,i)
  ENDDO
  DO k = 1, size(humidity,1)
    PRINT *, my_cart_id, 'humidity', k, humidity(k,i)
  ENDDO
  PRINT *, my_cart_id, 't2m ',t2m(i)
  PRINT *, my_cart_id, 'hum2m ',hum2m(i)
  PRINT *, my_cart_id, 'psurf ',psurf(i)
  PRINT *, my_cart_id, 's_hgt ',s_hgt(i)
  PRINT *, my_cart_id, 'u10m ',u10m(i)
  PRINT *, my_cart_id, 'v10m ',v10m(i)
  PRINT *, my_cart_id, 't_s ',t_s(i)
  PRINT *, my_cart_id, 'stype ',stype(i)
  PRINT *, my_cart_id, 'lat ',rlat(i,j) * raddeg
  PRINT *, my_cart_id, 'sat_z ',sat_z(i)
  PRINT *, my_cart_id, 'sun_z ',sun_el(i,j)
  PRINT *, my_cart_id, 'sat_a ',sat_a(i)
  PRINT *, my_cart_id, 'cloud'
  DO k = 1, size(cloud,2)
     WRITE (*,'(6(1x,E13.6))') cloud(:, k, i)
  ENDDO
  PRINT*,'cloud cover'
  DO k = 1, size(cloud,2)
     WRITE (*,'(6(1x,E13.6))') clc(:,k , i)
  ENDDO

  ! output for tstrad.exe (part of RTTOV)
  OPEN(10,file='input.txt')  !For SEVIRI MSG-2
  WRITE(10,*) 12
  WRITE(10,*) 2
  WRITE(10,*) 21
  WRITE(10,*) 0
  WRITE(10,*) 1
  WRITE(10,*) stype(i)
  WRITE(10,*) 1
  WRITE(10,*) 0
  WRITE(10,*) 1
  WRITE(10,*) 4
  WRITE(10,*) 2
  WRITE(10,*) sat_z(i)
  WRITE(10,*) sat_a(i)
  WRITE(10,*) 0
  WRITE(10,*) rlat(i,j) * raddeg
  WRITE(10,*) s_hgt(i)
  WRITE(10,*) 0
  CLOSE(10)

  OPEN(10, file='plevs.dat')
  WRITE(10,*) 1,SIZE(pres,1)
  DO k=1, SIZE(pres,1)
    WRITE(10,*) pres(k,i)
  ENDDO
  CLOSE(10)

  OPEN(10, file='prof.dat')
  DO k=1, SIZE(pres,1)
    WRITE(10,*) temp(k,i)
  ENDDO
  DO k=1, SIZE(pres,1)
    WRITE(10,*) humidity(k,i)
  ENDDO
  DO k=1, SIZE(pres,1)
    WRITE(10,*) 0.0_wp
  ENDDO
  DO k=1, SIZE(pres,1)
    WRITE(10,*) 0.0_wp
  ENDDO
  WRITE(10,'(6(" ",D13.6))') t2m(i), hum2m(i), psurf(i), u10m(i), v10m(i), &
         100000.0_wp
  WRITE(10,'(7(" ",D13.6))') t_s(i), 3.0_wp, 5.0_wp, 15.0_wp, 0.1_wp, 0.30_wp
  WRITE(10,*) 500.0_wp, 0.0_wp
  CLOSE(10)

  OPEN(10, file='prof_cloud.dat')
  DO k=1, SIZE(cloud,2)
    WRITE(10,'(6(" ",D13.6))') (cloud(l,k,i), l=1,6)
  ENDDO
  close(10)
  open(10, file='prof_cfrac.dat')
  DO k=1, SIZE(clc,2)
    WRITE(10,'(6(" ",D13.6))') (clc(l,k,i), l=1,6)
  ENDDO
  CLOSE(10)

END SUBROUTINE print_profile

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE org_sat_rttov

!==============================================================================
#endif

END MODULE src_sat_rttov
