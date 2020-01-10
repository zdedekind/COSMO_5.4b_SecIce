!+ Source Module for generation of synthetic satellite images using RTTOV7
!------------------------------------------------------------------------------

MODULE src_sat_tbs

!------------------------------------------------------------------------------
!
! Description:
!   1- interpolate the T, qv, ... profiles to the 43-level RTTOV grid
!   2- run rttovcld direct model
!
!    This software was developed within the context of 
!    the EUMETSAT Satellite Application Facility on 
!    Numerical Weather Prediction (NWP SAF), under the 
!    Cooperation Agreement dated 25 November 1998, between 
!    EUMETSAT and the Met Office, UK, by one or more partners 
!    within the NWP SAF. The partners in the NWP SAF are 
!    the Met Office, ECMWF, KNMI and MeteoFrance.
! 
!    Copyright 2002, EUMETSAT, All Rights Reserved.
!    
! Method:
!   see comments in program
!       
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!       
! History:
! Version    Date       Name
! ---------- ---------- ----
! 3.7        2004/02/18 Christian Keil (DLR) / Ulrich Schaettler
!  Initial Release
! 3.10       2004/06/23 Ulrich Schaettler
!  Changed treatment of error values from RTTOV-library
! 3.21       2006/12/04 Ulrich Schaettler
!  Bug Correction for satellite zenith angle  (F. Ament, C. Keil)
! V3_23        2007/03/30 Ulrich Schaettler, Lucio Torrisi
!  Limit the value for the zenith angle to -89.5 ... 89.5
! V4_3         2008/02/25 Ulrich Schaettler
!  Fit linear extrapolation of temperature above modellevels into soft limits
!  Avoid too much debug printout
!  Introduced a limitation for the surface pressure.
! V4_5         2008/09/10 Ulrich Schaettler
!  Restructured org_sat_tbs not to work on a single, but on a vector of profiles
!  For that, a modified RTTOV library is needed!
! V4_8         2009/02/16 Guenther Zaengl
!  Use p0hl (reference pressure at half levels) for full consistency with
!  new reference atmosphere implementation
! V4_9         2009/07/16 Ulrich Schaettler
!  Adaptations to (additional) vectorizations in RTTOV7
!  Included dimensional and computational indices to parameter list of org_sat_tbs
!   (which gives the option to call org_sat_tbs for one long vector)
!  Moved call to RTTVI to organize_satellites
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_18        2011/05/26 Ulrich Schaettler
!  Bug fix for computing the extrapolation below the model surface pressure;
!  the variable index must be set to ke+1
! V4_19        2011/08/01 Ulrich Schaettler
!  Introduced conditional compilation for RTTOV7
! V4_28        2013/07/12 Ulrich Schaettler
!  Renamed nlist_chan to nchan_list
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
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
!------------------------------------------------------------------------------

USE data_constants  , ONLY :  &
    r_earth

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :  &
    raddeg

!------------------------------------------------------------------------------

USE data_parallel   , ONLY :  &
    my_cart_id, imp_integers, icomm_cart

!------------------------------------------------------------------------------

USE data_parameters , ONLY :  &
    wp,        &! KIND-type parameters for real variables
    dp,        &! KIND-type parameters for real variables (double precision)
    iintegers   ! kind-type parameter for "normal" integer variables

!------------------------------------------------------------------------------

USE data_satellites , ONLY :  &
    sat_compute, num_sensors, numchans, maxknchpf, jpch, rcnv, rcnw,        &
    utmn, uqmn, ppres_d, o3_ref, const_aheat, r_sat

!------------------------------------------------------------------------------

USE parallel_utilities, ONLY :  global_values

!==============================================================================

IMPLICIT NONE

!==============================================================================

#ifdef RTTOV7
INCLUDE "RTTOV7_SPLINE_vec.interface"
INCLUDE "RTTOV7_SUPSAT_vec.interface"
INCLUDE "RTTOV7_RTTOVCLD.interface"
#endif

!==============================================================================

CONTAINS

!==============================================================================        

SUBROUTINE org_sat_tbs (t, qv, clw, qi, pp, p0, p0hl, clc, ps, t_g, t_2m,   &
                        qv_2m, u10m, v10m, fr_land, rlat, rlon,             &
                        idim, jdim, ke, idims, idime, jdims, jdime,         &
                        nlevels, nprof, nsensors , npav, nsav, nssv,        &
                        synme7, synmsg, yerror, ierror)

!------------------------------------------------------------------------------

INTEGER (KIND=iintegers), INTENT(IN) ::           &
  idim, jdim, ke, & !
  idims, idime, jdims, jdime, & ! LM org. variables
  nlevels,  & ! No. of pressure levels used in the rttov-library
  nprof,    & ! No. of profiles to be computed per rttov-CALL
  nsensors, & ! Max. No. of sensors to be used
  npav,     & ! No. of profile variables       (from Mod_CParam)
  nsav,     & ! No. of surface air variables   (from Mod_CParam)
  nssv        ! No. of skin variables          (from Mod_CParam)

INTEGER (KIND=iintegers), INTENT(OUT) ::           &
  ierror      ! error status variable

CHARACTER (LEN= *),       INTENT(OUT) ::           &
  yerror       ! error message

REAL  (KIND=wp),     INTENT(IN) ::           &
  t      (idim,jdim,ke),     &
  qv     (idim,jdim,ke),     &
  clw    (idim,jdim,ke),     &
  qi     (idim,jdim,ke),     &
  pp     (idim,jdim,ke),     &
  p0     (idim,jdim,ke),     &
  p0hl   (idim,jdim,ke+1),   &
  clc    (idim,jdim,ke),     &
  ps     (idim,jdim),        &
  t_g    (idim,jdim),        &
  t_2m   (idim,jdim),        &
  qv_2m  (idim,jdim),        &
  u10m   (idim,jdim),        &
  v10m   (idim,jdim),        &
  fr_land(idim,jdim),        &
  rlat   (idim,jdim),        &
  rlon   (idim,jdim)

REAL  (KIND=wp),     INTENT(OUT) ::          &
  synme7 (idim,jdim, 8),     &
  synmsg (idim,jdim,32)

! Local parameters:
!------------------

REAL  (KIND=wp)      ::              &
  presfin    (idim,nlevels+ke),      &
  wvmmfin    (idim,nlevels+ke),      &
  wvvmfin    (idim,nlevels+ke),      &
  tfin       (idim,nlevels+ke),      &
  psurfmod   (idim),                 &

  as         (idim,5,nlevels+ke),    &

  aryint     (idim,nlevels),         &
  textr      (idim,nlevels),         &
  tint       (idim,nlevels),         &
  wvvmint    (idim,nlevels),         &

  wvmm       (idim,ke+1),            &
  wvs        (idim),                 &
  rh         (idim),                 &
  gradt      (idim),                 &
  gradq      (idim),                 &

  pav_d      (nprof,nlevels,npav),   &
  psav_d     (nprof,nsav),           &
  pssv_d     (nprof,nssv),           &
  pangl_d    (nprof),                &
  pangs_d    (nprof),                &
  pcvm_d     (nprof,ke,4),           &
  papf_d     (idim,ke),              &
  papf_dd    (nprof,ke),             &
  paph_d     (idim,ke+1),            &
  paph_dd    (nprof,ke+1)

INTEGER (KIND=iintegers)   ::        &
  nlchan      (maxknchpf)          , &
  kprof       (maxknchpf)          , &
  kchan       (maxknchpf)

REAL  (KIND=wp)            ::        &
  emis        (maxknchpf)          , &
  pemis_d     (maxknchpf)          , &
  tau_d       (maxknchpf,nlevels)  , &
  tausfc_d    (maxknchpf)          , &
  ptb_d       (maxknchpf)          , &
  ptbcld_d    (maxknchpf)          , &
  prad_d      (maxknchpf)          , &
  pradcld_d   (maxknchpf)          , &
  radovm_d    (maxknchpf,2*ke+2)   , &
  pcldemis_d  (maxknchpf,ke)       , &
  pait_d      (maxknchpf,ke+1)     , &
  pais_d      (maxknchpf,ke+1)

INTEGER (KIND=iintegers) ::          &
  ksurf      (nprof ),               & ! Surface type index
  ifail      (nprof,nsensors),       & ! error status from rttovcld
  ierrorcode (idim,jdim,nsensors),       & ! for gathering the error status per grid point
  itwarnflag (idim,jdim),                & ! for warning if temperature is adjusted
  iqwarnflag (idim,jdim)                   ! for warning if humidity is adjusted

REAL  (KIND=wp)     ::           &
  tempsup(idim*jdim*nlevels), wsup(idim*jdim*nlevels),psup(idim*jdim*nlevels)

! Local scalars:
!---------------

REAL  (KIND=wp)     ::           &
  t2m, q2m, psurf, wvmmextr

! Interface to RTTOVCLD
INTEGER (KIND=iintegers) ::           &

  ! Input
!! that is nprof:    knpf,      & ! Number of profiles
!! that is nlevels:  klenpf,    & ! Length of input  profile vectors
  knchpf       ! Number of output radiances
               !    (= channels used * profiles)
!! that is ke:  klevm,     & ! Number of model (native) levels


INTEGER (KIND=iintegers) ::                  &
  levbot(idim),                              &
  levtop(idim),                              &
  index (idim)

INTEGER (KIND=iintegers) ::                  &
  i, j, k, kl, n, isens, izmaxerror, izmaxtwarn, izmaxqwarn,  &
  iatm, ichan, izerr, igatherbuf(3), klu, itest

REAL (KIND=wp)     :: alpha_e, r_atm

CHARACTER (LEN=80) :: yzerrmsg

!-----End of header------------------------------------------------------------

#ifdef RTTOV7
!------------------------------------------------------------------------------
! Section 1: Initialisations
!------------------------------------------------------------------------------

  ! distance of satellite to middle of the earth
  r_sat       = 35880E3_wp + r_earth

  ierror     = 0_iintegers
  yerror     = '         '
  ifail(:,:) = 0_iintegers
  ierrorcode(:,:,:) = 0_iintegers
  itwarnflag(:,:)   = 0_iintegers
  iqwarnflag(:,:)   = 0_iintegers

  pangs_d(:) = 0.0_wp
  pangl_d(:) = 0.0_wp

  izerr      = 0_iintegers
       
  nlchan(1:numchans(1)) = sat_compute(1)%nchan_list(1:numchans(1))
  nlchan(1:numchans(2)) = sat_compute(2)%nchan_list(1:numchans(2))
  emis  (1:numchans(1)) = sat_compute(1)%emissivity(1:numchans(1))
  emis  (1:numchans(2)) = sat_compute(2)%emissivity(1:numchans(2))
        
  !----------------------------------------------------------------------------
  ! Section 2: Compute the profiles
  !----------------------------------------------------------------------------

  jloop:  DO j = jdims, jdime    ! order of the loop affects the order of output

    ! Some explanations to the meaning of variables:
    !  nlevels: (=43) number of pressure levels in RTTOV (called SAT levels later on)
    !  levtop:        highest SAT level of ppres_d below model top (papf_d(1))
    !  levbot:        lowest SAT level of ppres_d above model surface pressure (psurfmod)
    !  index:         number of active model levels, which have been eventually enlarged
    !                 by an extrapolation above highest level and / or extrapolation
    !                 below surface pressure

    !--------------------------------------------------------------------------
    ! Section 2.1: Preset local variables (on ke modellevels + surface values)
    !--------------------------------------------------------------------------

    DO k = 1, ke
      DO i = idims, idime  ! of the Sat-BTs: longitudes the innermost

        papf_d (i,k) = (p0 (i,j,k) + pp(i,j,k)) / 100.0_wp ! pressure in hPa

        ! Specific humidity(kg/kg) > Mass mixing ratio (g/g)
        wvmm   (i,k) = (qv(i,j,k) / (1.0_wp - qv(i,j,k)))
        wvmmfin(i,k) =  wvmm  (i,k)
        presfin(i,k) =  papf_d(i,k)
        tfin   (i,k) =  t     (i,j,k)
      ENDDO
    ENDDO

    DO k = 2, ke
      DO i = idims, idime  ! of the Sat-BTs: longitudes the innermost
        paph_d(i,k) = 0.5_wp * (papf_d(i,k-1)+papf_d(i,k))
      ENDDO
    ENDDO

    ! and all surface (and toa) fields
    DO i = idims, idime  ! of the Sat-BTs: longitudes the innermost
      ! specific ozone (kg/kg) set to climat. background value   -----CK

      psurf = ps     (i,j) / 100.0_wp     ! surface pressure (hPa)
      IF (psurf >= 1100.0_wp) THEN
        ! happened once at the Dead Sea
        ! an inconsistency to the upper pressure levels could occur; 
        ! but this should do no harm
        psurf = 1099.9_wp
      ENDIF
      psurfmod(i) = psurf    ! modified surface pressure

      t2m   = t_2m   (i,j)     ! 2-meter temperature (K)
      q2m   = qv_2m  (i,j)     ! 2-meter specific humidity (kg/kg)

      paph_d  (i,   1) = (p0hl(i,j,1) + 0.5_wp * pp(i,j,1)) / 100.0_wp
      paph_d  (i,ke+1) = psurf
      presfin (i,ke+1) = psurf
      wvmm    (i,ke+1) = (q2m / (1.0_wp - q2m))
      wvmmfin (i,ke+1) = wvmm(i,ke+1)
      tfin    (i,ke+1) = t2m
    ENDDO

  !----------------------------------------------------------------------------
  ! Section 2.2: Extrapolation below surface pressure
  !----------------------------------------------------------------------------

    levbot(:) = 0_iintegers
    index (:) = ke + 1     ! US (May 2011): this was ke before, but the surface
                           ! values have to be considered as an extra level!!!

    DO kl = nlevels, 1, -1
      DO i = idims, idime  ! of the Sat-BTs: longitudes the innermost
        IF (ppres_d(kl,1) > psurfmod(i)) THEN
          levbot(i) = kl
        ENDIF
      ENDDO
    ENDDO
        
    DO kl = 1, nlevels
      DO i = idims, idime  ! of the Sat-BTs: longitudes the innermost
        IF ( (levbot(i) /= 0_iintegers) .AND. (kl >= levbot(i)) ) THEN

          textr(i,kl) = t_2m(i,j) * EXP (const_aheat * LOG (ppres_d(kl,1) / psurfmod(i)))
      
          index(i)  =  nlevels - levbot(i) + 1     ! number of SAT levels below model surface
          index(i)  =  index(i) + ke + 1           ! number of "extended model levels"

          tfin   (i,ke+2-levbot(i)+kl) = textr  (i,kl)
          presfin(i,ke+2-levbot(i)+kl) = ppres_d(  kl,1)

        ENDIF
      ENDDO
    ENDDO

    CALL RTTOV7_supsat_vec (t_2m(:,j), wvs, psurfmod, idim, idims, idime)

    ! Extrapolate water vapour below surface pressure
    ! => constant relative humidity
        
    n = 0
    DO kl = 1, nlevels
      DO i = idims, idime  ! of the Sat-BTs: longitudes the innermost
        IF (levbot(i) /= 0_iintegers) THEN
          rh(i) = wvmm(i,ke+1) * 1000.0_wp / wvs(i)
          IF (kl >= levbot(i)) THEN
            n = n+1
            ! gather values for RTTOV7_supsat
            tempsup(n) = textr  (i,kl)
            psup(n)    = ppres_d(  kl,1)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
        
    CALL RTTOV7_supsat_vec (tempsup, wsup, psup, n, 1, n)

    n = 0
    DO kl = 1, nlevels
      DO i = idims, idime  ! of the Sat-BTs: longitudes the innermost
        IF ( (levbot(i) /= 0_iintegers) .AND. (kl >= levbot(i)) ) THEN
          n = n+1
          wvmmextr = rh(i) * wsup(n) / 1000.0_wp
          wvmmfin(i,ke+2-levbot(i)+kl) = wvmmextr
        ENDIF
      ENDDO
    ENDDO

    !----------------------------------------------------------------------------
    ! Section 6: Extrapolation above highest declared level
    !----------------------------------------------------------------------------
        
    levtop(:) = 1_iintegers

    DO kl = 1, nlevels + ke  ! +2
      DO i = idims, idime  ! of the Sat-BTs: longitudes the innermost
        IF (kl <= index(i)) THEN
          wvvmfin(i,kl)=wvmmfin(i,kl) * rcnw
        ENDIF
      ENDDO
    ENDDO
        
    ! => linear extrapolation
    !US only possible from nlevels on  DO kl = nlevels + ke+2, 1, -1
    DO kl = nlevels, 1, -1
      DO i = idims, idime  ! of the Sat-BTs: longitudes the innermost
        IF (kl <= index(i)) THEN
          IF (ppres_d(kl,1) >= papf_d(i,1) ) levtop(i) = kl
        ENDIF
      ENDDO
    ENDDO
        
    DO kl = nlevels + ke, 1, -1
      ! move values down in the structures, to free "levtop" levels above
      DO i = idims, idime  ! of the Sat-BTs: longitudes the innermost
        IF ( (levtop(i) /= 1_iintegers) .AND. (kl <= index(i)) ) THEN
          tfin   (i,kl+levtop(i)-1) = tfin   (i,kl)
          wvvmfin(i,kl+levtop(i)-1) = wvvmfin(i,kl)
          presfin(i,kl+levtop(i)-1) = presfin(i,kl)
        ENDIF
      ENDDO
    ENDDO

    DO i = idims, idime  ! of the Sat-BTs: longitudes the innermost
      IF (levtop(i) /= 1) THEN
        index(i) = index(i) + levtop(i) -1   ! increase number of "active" levels

        ! Take the gradient from the 2 highest modellevels
        ! This gradient can sometimes be rather big for the temperature.
        ! Linear extrapolation to higher pressure levels can then be
        ! problematic!
        gradt(i) = (tfin(i,1) - tfin(i,2)) / (presfin(i,1)-presfin(i,2))
        gradq(i) = (wvvmfin(i,1) - wvvmfin(i,2)) / (presfin(i,1)-presfin(i,2))
      ENDIF
    ENDDO

    ! Extrapolation above highest model level
    DO kl = 1, nlevels + ke
      ! set values for the "levtop" new values above
      DO i = idims, idime  ! of the Sat-BTs: longitudes the innermost
        IF ( (levtop(i) /= 1_iintegers) .AND. (kl <= levtop(i)-1) ) THEN
          klu = MIN (kl, levtop(i)-1)

          presfin(i,kl) = ppres_d(kl,1)
          tfin   (i,kl) = tfin(i,levtop(i)) + gradt(i) * (presfin(i,kl)-presfin(i,levtop(i)))

          IF (tfin(i,kl) < (utmn(klu,1) + 5.0_wp)) THEN
            ! add a security value to the final temperature and set warnflag
            tfin(i,kl) = utmn(klu,1) + 5.0_wp
            itwarnflag(i,j) = 1
          ENDIF

          wvvmfin(i,kl) = wvvmfin(i,levtop(i)) + gradq(i) * (presfin(i,kl)-presfin(i,levtop(i)))
        ENDIF
      ENDDO
    ENDDO
        
    !--------------------------------------------------------------------------
    ! Section 7: Interpolation to given pressure grid
    !--------------------------------------------------------------------------
        
    CALL RTTOV7_spline_vec (presfin(:,:), tfin(:,:), ppres_d(:,1),        &
                        aryint(:,:), as(:,:,:), index(:),                 &
                        nlevels+ke, nlevels, idim, idims, idime, 1)
    tint(idims:idime,:)      = aryint(idims:idime,:)

    CALL RTTOV7_spline_vec (presfin(:,:), wvvmfin(:,:), ppres_d(:,1),     &
                        aryint(:,:), as(:,:,:), index(:),                 &
                        nlevels+ke, nlevels, idim, idims, idime, 1)
    wvvmint(idims:idime,:) = aryint(idims:idime,:) / rcnw

    DO kl = 1, nlevels
      DO i = idims, idime  ! of the Sat-BTs: longitudes the innermost
        ! Spline interpolation may exceptionnally give negative values
        ! for water vapour profiles   => correction
        IF (wvvmint(i,kl) <= uqmn(kl,1)) THEN
          ! add a security value to the final humidity and set warnflag
          wvvmint(i,kl) = uqmn(kl,1) + 1.0E-7_wp
          iqwarnflag(i,j) = 1
        ENDIF
      ENDDO
    ENDDO
        
    ! Exploitation of the profile
    ! At this stage of the IATM loop, the following arrays are 
    ! available on the ppres_d grid :
    !        - ppres_d   in Pa             
    !        - tint      in K               
    !        - wvvmint   in ppmv           
    !        - o3vmint   in ppmv          
    ! In addition, the following scalars are still available : 
    !        - rlsm                  (fr_land)
    !        - st        in K        (t_g)
    !        - psurf     in hPa 
    !        - t2m       in K     
    !        - q2m       in kg/kg
    
    !----------------------------------------------------------------------------
    ! Section 8:  Fill in RTTOV input arrays
    !----------------------------------------------------------------------------
        
    DO kl = 1, nlevels
      iatm = 0
      DO i = idims, idime  ! of the Sat-BTs: longitudes the innermost
        iatm = iatm + 1
        pav_d(iatm,kl,1) = tint   (i,kl)
        pav_d(iatm,kl,2) = wvvmint(i,kl)
        pav_d(iatm,kl,3) =  o3_ref(  kl)
        pav_d(iatm,kl,4) = 0.0_wp
      ENDDO
    ENDDO

    iatm = 0
    DO i = idims, idime  ! of the Sat-BTs: longitudes the innermost
      iatm = iatm + 1

      psav_d(iatm,1) = t_2m (i,j)
      psav_d(iatm,2) = qv_2m(i,j)
      psav_d(iatm,3) = psurfmod(i)
      psav_d(iatm,4) = u10m(i,j)
      psav_d(iatm,5) = v10m(i,j)
      pssv_d(iatm,1) = t_g(i,j)    ! st
      pssv_d(iatm,2) = 3.0_wp
      pssv_d(iatm,3) = 5.0_wp
      pssv_d(iatm,4) = 15.0_wp
      pssv_d(iatm,5) = 0.1_wp
      pssv_d(iatm,6) = 0.3_wp
        
      pcvm_d(iatm,:,1) =   t(i,j,:)
      pcvm_d(iatm,:,2) = clc(i,j,:)
      pcvm_d(iatm,:,3) = clw(i,j,:)
      pcvm_d(iatm,:,4) =  qi(i,j,:)

      papf_dd(iatm,:)  = papf_d(i,:)
      paph_dd(iatm,:)  = paph_d(i,:)

      ! 0=land    1=sea    2=snow/ice
      ksurf(iatm) = INT(1.0_wp - fr_land(i,j))   ! RWS mod to sea
        
      ! Calculate the zenith angle pangl_d
      alpha_e = ACOS (COS(rlat(i,j)) * COS(rlon(i,j)))
      r_atm   = SQRT (r_sat**2 +r_earth**2 -2*r_sat*r_earth*COS(alpha_e))
      pangl_d(iatm) = ASIN (SIN(alpha_e)*r_sat/r_atm)*raddeg

      ! Limit the zenith angle to values between -89.5 ... 89.5
      ! (otherwise there will be problems at the poles; and with these values
      !  the results at the poles should not be too wrong)
      IF (pangl_d(iatm) < -89.5_wp) pangl_d(iatm) = -89.5_wp
      IF (pangl_d(iatm) >  89.5_wp) pangl_d(iatm) =  89.5_wp
    ENDDO

    !------------------------------------------------------------------------
    ! Section 9:  Call to RTTOVCLD
    !------------------------------------------------------------------------
        
    sensor_loop: DO isens = 1, num_sensors

      IF ( (j == jdims) .AND. (my_cart_id == 0) ) THEN
        PRINT *, '   Computing: ', sat_compute(isens)%ysatellite,      &
                                   sat_compute(isens)%nsat_id,         &
                                   sat_compute(isens)%ysensor
      ENDIF

      ! Initial values for the surface emissivities (have to be reinitialized for
      ! every grid point, because pemis_d is INOUT in rttovcld
      DO ichan = 1, numchans(isens)
        iatm = 0
        DO i = idims, idime  ! of the Sat-BTs: longitudes the innermost
          iatm = iatm+1
          pemis_d(numchans(isens)*(iatm-1)+ichan) = emis(ichan)

          kprof(numchans(isens)*(iatm-1)+1:numchans(isens)*iatm) = iatm
          kchan(numchans(isens)*(iatm-1)+1:numchans(isens)*iatm) = nlchan(1:numchans(isens))
        ENDDO
      ENDDO

      knchpf = nprof * numchans(isens)

      CALL RTTOV7_rttovcld                                              &
        (nprof,nlevels,ke, ppres_d(:,isens),pangl_d,                    &
         pangs_d,ksurf,isens,knchpf,kchan,kprof,                        &
         pav_d,psav_d,pssv_d,pcvm_d, papf_dd,paph_dd,                   &
         pemis_d,ifail,prad_d,ptb_d,pradcld_d,ptbcld_d,tau_d,tausfc_d,  &
         radovm_d,pcldemis_d,pait_d,pais_d)
        
      iatm = 0
      DO i = idims, idime  ! of the Sat-BTs: longitudes the innermost
        iatm = iatm + 1

        ! Set errorcode
        IF (ifail(iatm,isens) /= 0_iintegers) THEN
          ierrorcode(i,j,isens) = ifail(iatm,isens)
        ENDIF

        ! main output:
        ! ------------
        ! pradcld_d (knchpf) = cloud-affected radiances
        ! prad_d    (knchpf) = clear-sky radiances
        ! ptbcld_d  (knchpf) = cloud-affected Tbs 
        ! ptb_d     (knchpf) = clear-sky Tbs
        ! ichan    : WV
        ! ichan+1  : IR

        SELECT CASE (sat_compute(isens)%ysatellite)
        CASE('METEOSAT')
          IF (sat_compute(isens)%lcloud_tem) THEN
            synme7(i,j, 1) = ptbcld_d (2*iatm - 1)     ! Tb cloudy     WV Channel 6.4
          ENDIF
          IF (sat_compute(isens)%lclear_tem) THEN
            synme7(i,j, 2) = ptb_d    (2*iatm - 1)     ! Tb clearsky   WV Channel
          ENDIF
          IF (sat_compute(isens)%lcloud_rad) THEN
            synme7(i,j, 3) = pradcld_d(2*iatm - 1)     ! Rad cloudy    WV Channel
          ENDIF
          IF (sat_compute(isens)%lclear_rad) THEN
            synme7(i,j, 4) = prad_d   (2*iatm - 1)     ! Rad clearsky  WV Channel
          ENDIF

          IF (sat_compute(isens)%lcloud_tem) THEN
            synme7(i,j, 5) = ptbcld_d (2*iatm)     ! Tb cloudy     IR Channel  11.5
          ENDIF
          IF (sat_compute(isens)%lclear_tem) THEN
            synme7(i,j, 6) = ptb_d    (2*iatm)     ! Tb clearsky   IR Channel
          ENDIF
          IF (sat_compute(isens)%lcloud_rad) THEN
            synme7(i,j, 7) = pradcld_d(2*iatm)     ! Rad cloudy    IR Channel
          ENDIF
          IF (sat_compute(isens)%lclear_rad) THEN
            synme7(i,j, 8) = prad_d   (2*iatm)     ! Rad clearsky  IR Channel
          ENDIF

        CASE('MSG     ')
          synmsg(i,j, 1) = ptbcld_d  (8*iatm-7)     ! Tb  cloudy    msg( 3) WV1.6
          synmsg(i,j, 2) = ptb_d     (8*iatm-7)     ! Tb  clear sky msg( 3) WV1.6
          synmsg(i,j, 3) = pradcld_d (8*iatm-7)     ! Rad cloudy    msg( 3) WV1.6
          synmsg(i,j, 4) = prad_d    (8*iatm-7)     ! Rad clear sky msg( 3) WV1.6

          synmsg(i,j, 5) = ptbcld_d  (8*iatm-6)     ! Tb  cloudy    msg( 4) WV6.2
          synmsg(i,j, 6) = ptb_d     (8*iatm-6)     ! Tb  clear sky msg( 4) WV6.2
          synmsg(i,j, 7) = pradcld_d (8*iatm-6)     ! Rad cloudy    msg( 4) WV1.6
          synmsg(i,j, 8) = prad_d    (8*iatm-6)     ! Rad clear sky msg( 4) WV1.6

          synmsg(i,j, 9) = ptbcld_d  (8*iatm-5)     ! Tb  cloudy    msg( 5) WV7.3
          synmsg(i,j,10) = ptb_d     (8*iatm-5)     ! Tb  clear sky msg( 5) WV7.3
          synmsg(i,j,11) = pradcld_d (8*iatm-5)     ! Rad cloudy    msg( 5) WV1.6
          synmsg(i,j,12) = prad_d    (8*iatm-5)     ! Rad clear sky msg( 5) WV1.6

          synmsg(i,j,13) = ptbcld_d  (8*iatm-4)     ! Tb  cloudy    msg( 6) IR3.8
          synmsg(i,j,14) = ptb_d     (8*iatm-4)     ! Tb  clear sky msg( 6) IR3.8
          synmsg(i,j,15) = pradcld_d (8*iatm-4)     ! Rad cloudy    msg( 6) WV1.6
          synmsg(i,j,16) = prad_d    (8*iatm-4)     ! Rad clear sky msg( 6) WV1.6

          synmsg(i,j,17) = ptbcld_d  (8*iatm-3)     ! Tb  cloudy    msg( 7) IR8.7
          synmsg(i,j,18) = ptb_d     (8*iatm-3)     ! Tb  clear sky msg( 7) IR8.7
          synmsg(i,j,19) = pradcld_d (8*iatm-3)     ! Rad cloudy    msg( 7) WV1.6
          synmsg(i,j,20) = prad_d    (8*iatm-3)     ! Rad clear sky msg( 7) WV1.6

          synmsg(i,j,21) = ptbcld_d  (8*iatm-2)     ! Tb  cloudy    msg( 8) IR10.8
          synmsg(i,j,22) = ptb_d     (8*iatm-2)     ! Tb  clear sky msg( 8) IR10.8
          synmsg(i,j,23) = pradcld_d (8*iatm-2)     ! Rad cloudy    msg( 8) WV1.6
          synmsg(i,j,24) = prad_d    (8*iatm-2)     ! Rad clear sky msg( 8) WV1.6

          synmsg(i,j,25) = ptbcld_d  (8*iatm-1)     ! Tb  cloudy    msg( 9) IR12.0
          synmsg(i,j,26) = ptb_d     (8*iatm-1)     ! Tb  clear sky msg( 9) IR12.0
          synmsg(i,j,27) = pradcld_d (8*iatm-1)     ! Rad cloudy    msg( 9) WV1.6
          synmsg(i,j,28) = prad_d    (8*iatm-1)     ! Rad clear sky msg( 9) WV1.6

          synmsg(i,j,29) = ptbcld_d  (8*iatm  )     ! Tb  cloudy    msg(10) IR9.7
          synmsg(i,j,30) = ptb_d     (8*iatm  )     ! Tb  clear sky msg(10) IR9.7
          synmsg(i,j,31) = pradcld_d (8*iatm  )     ! Rad cloudy    msg(10) WV1.6
          synmsg(i,j,32) = prad_d    (8*iatm  )     ! Rad clear sky msg(10) WV1.6
        END SELECT

      ENDDO

    ENDDO sensor_loop

  ENDDO  jloop

!----------------------------------------------------------------------------
! Section 10: Cleanup of RTTOV
!----------------------------------------------------------------------------

! Check error and warning flags
igatherbuf(1) = MAXVAL(itwarnflag(:,:))
igatherbuf(2) = MAXVAL(iqwarnflag(:,:))
igatherbuf(3) = MAXVAL(ierrorcode(:,:,:))

CALL global_values (igatherbuf, 3, 'MAX', imp_integers, icomm_cart,    &
                    0_iintegers, yzerrmsg, izerr)

izmaxtwarn = igatherbuf(1)
izmaxqwarn = igatherbuf(2)
izmaxerror = igatherbuf(3)

! Gather these values on PE 0

IF (my_cart_id == 0) THEN
  IF (izmaxtwarn > 0_iintegers) THEN
    PRINT *, '   REMARK:  Some temperature profiles adjusted to soft limits',  &
                          ' during Extrapolation'
  ENDIF

  IF (izmaxqwarn > 0_iintegers) THEN
    PRINT *, '   REMARK:  Some humidity    profiles adjusted to soft limits'
  ENDIF

  IF ( (izmaxerror > 0_iintegers) .AND. (izmaxerror <= 19_iintegers) ) THEN
    PRINT *, '   WARNING: Some soft limits have been violated'
  ENDIF

  IF (izmaxerror > 19_iintegers) THEN
    PRINT *, '   ERROR:   Some hard limits have been violated'
    ierror = 1_iintegers
  ENDIF
ENDIF

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------
#endif

END SUBROUTINE org_sat_tbs

!==============================================================================

END MODULE src_sat_tbs
