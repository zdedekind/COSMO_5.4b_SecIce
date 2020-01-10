!+ Source module for turbulence utility routines
!==============================================================================

MODULE  turbulence_utilities

!==============================================================================
!
! Description:
! 
!   Routines (module procedures) currently contained:
!  
! History:
! Version      Date       Name
! ------------ ---------- ----
! V4_20        2011/08/31 Ulrich Schaettler
!  Initial release
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
! V5_4a        2016-05-10 Matthias Raschendorfer
!  Use new module turb_data instead of data_turbulence
!  Define variables from former data_turbulence, which are only used in the
!   old turbulence scheme
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
!------------------------------------------------------------------------------

USE data_parameters , ONLY :   &
  wp,        & ! KIND-type parameters for real variables
  iintegers    ! KIND-type parameter for standard integer variables

USE turb_data       , ONLY :   &
  pat_len, c_sea, c_soil, c_lnd, e_surf, q_crit,      &
  clc_diag, d_mom, a_mom, d_heat, tur_len, a_heat

!------------------------------------------------------------------------------

IMPLICIT NONE

!------------------------------------------------------------------------------

! Definitions from former data_turbulence
REAL (KIND=wp),     TARGET :: &
     c_tke,tet_g,c_g,rim, &
     d_0,d_1,d_2,d_3,d_4,d_5,d_6, &
     a_3,a_5,a_6,b_1,b_2, &
     l_scal,                      &
     l_hori                         ! horizontal grid spacing

!------------------------------------------------------------------------------

CONTAINS

!==============================================================================
!==============================================================================
!+ Module procedure canopy_source in "src_turbdiff" for calculation                   
!+ of scalar source terms inside the model canopy        

SUBROUTINE canopy_source

!------------------------------------------------------------------------------
!
! Description:
!
! Method:
!
!------------------------------------------------------------------------------

END SUBROUTINE canopy_source

!==============================================================================
!==============================================================================
!+ Module procedure diag_level in "src_turbdiff" for computing the upper level index
!+ used for near surface diganostics

SUBROUTINE diag_level (i_st, i_en, j_st, j_en, ke1, zdia_2d, k_2d,      &
                       hk_2d, hk1_2d, hhl)

!------------------------------------------------------------------------------

   INTEGER (KIND=iintegers), INTENT(IN) :: &
!
      i_st, i_en, j_st, j_en,   & ! start end end indices of horizontal domain
      ke1                         ! vertical dimension of hhl

   REAL (KIND=wp),     INTENT(IN) :: &
!
      zdia_2d(:,:)  !diagnostic height

   INTEGER (KIND=iintegers), INTENT(INOUT) :: &
!
      k_2d(:,:)     !index field of the upper level index
                    !to be used of near surface diagnostics

   REAL (KIND=wp),     INTENT(INOUT) :: &
!
      hk_2d(:,:), & ! mid level height above ground belonging to 'k_2d'
     hk1_2d(:,:)    ! mid level height above ground of the previous layer (below)

   REAL (KIND=wp),     INTENT(IN)    :: &
      hhl(:,:,:)    ! height of half levels

   INTEGER (KIND=iintegers) :: i, j

   LOGICAL :: lcheck

!------------------------------------------------------------------------------

   lcheck=.TRUE. !check whether a diagnostic level is above the current layer

   DO WHILE (lcheck) !loop while previous layer had to be checked
      lcheck=.FALSE. !check next layer ony, if diagnostic level is at least once
                     !above the current layer
      DO j=j_st,j_en
      DO i=i_st,i_en
         IF (hk_2d(i,j)<zdia_2d(i,j) .AND. k_2d(i,j)>1) THEN !diagnostic level is above current layer
            lcheck=lcheck .OR. .TRUE. !for this point or any previous one the
                                      !diagnostic level is above the current layer
            k_2d(i,j)=k_2d(i,j)-1
            hk1_2d(i,j)=hk_2d(i,j)
            hk_2d(i,j)=(hhl(i,j,k_2d(i,j))+hhl(i,j,k_2d(i,j)+1))*0.5_wp-hhl(i,j,ke1)
          END IF
       END DO   
       END DO   

   END DO   

END SUBROUTINE diag_level

!==============================================================================
!==============================================================================
!+ Module procedure init_volume_canopy for initialization
!+ of special external parameters describing the surface canopy needed for the
!+ description of surface-to-atmosphere transfer and within canopy diffusion:

SUBROUTINE init_volume_canopy (ie, je, ke, ke1, kcm,                   &
                               istartpar, iendpar, jstartpar, jendpar, &
                               fr_land, d_pat, c_big, c_sml, r_air)

!------------------------------------------------------------------------------
!
! Description:
!
!   In the module 'init_volume_canopy' additional external parametr fields, 
!   which are used in the turbulence scheme 'turbdiff' (especially parameters 
!   for the physical description of the roughness canopy) are covered by the 
!   appropriate values by reading the refering parameter files and/or by making 
!   some diagnostic calculations using the known parameters. 
!   The 3-d fields of the canopy parameters are dynamically allocated using the maximum 
!   canopy hight within the model domain.
!
! Method:
!
!   For the present there exists no concept of generating those additional data. 
!   Thus the model will run either with an artificial canopy arcitecture or 
!   (for simplicity) without any vertically resolved canopy. But even in the 
!   latter case at least the allocation of the 3-d canopy data fields must be 
!   done, because the canopy concept is incorporated in the tubulent diffusion 
!   scheme 'turbdiff'.
!
!   NOTE: The height of the canopy layer and the canopy variables have
!         to be allocated before this routine is called!
!
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Declarations
!------------------------------------------------------------------------------

!Formal Parameters:

INTEGER (KIND=iintegers), INTENT(IN) :: &

! Horizontal and vertical sizes of the fields and related variables:
! --------------------------------------------------------------------

    ie,  & ! number of grid points in zonal      direction 
    je,  & ! number of grid points in meridional direction
    ke,  & ! number of main model levels (start index is 1)
    ke1, & ! number of half model levels (start index is 1) 
    kcm, & ! index of the lowest model layer higher than the canopy
    istartpar, iendpar, & ! zonal      start and end index including the model boundary lines
    jstartpar, jendpar    ! meridional start and end index including the model boundary lines

! External parameter fields:
! --------------------------

REAL (KIND=wp),           INTENT(IN)    :: &
    fr_land(ie,je)          !

REAL (KIND=wp),           INTENT(INOUT) :: &

    d_pat(ie,je),         & ! horizontal pattern length scale
    c_big(ie,je,kcm:ke1), & ! effective drag coefficient of canopy elements
                            ! larger than or equal to the turbulent length scale (1/m)
    c_sml(ie,je,kcm:ke1), & ! effective drag coefficient of canopy elements
                            ! smaller than the turbulent length scale            (1/m)
    r_air(ie,je,kcm:ke1)    ! log of air containing fraction of a gridbox inside
                            ! the canopy                                          (1)

! ----------------
! Local variables:
! ----------------

  INTEGER (KIND=iintegers) ::  &
    i,  j         !  loop indices

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine init_volume_canopy
!------------------------------------------------------------------------------

! Input of the external canopy-parameters:

! At this stage there is no concept of generating those Parameters (c_big, c_sml, r_air).
! They may be derived as functions of rbig, dbig, rsml, dsml, which either come from
! the primary external parameter files or may be derived from other primary external
! parameters like canopy-hight and -type:

! Provisional values for the canopy parameters:
  ! Uppermost canopy level has been determined before 
  ! and canopy fields are allocated yet
  c_big(:,:,:) = 0.0_wp ! isotr. drag-coeff. of big canopy-elem.
  c_sml(:,:,:) = 0.0_wp !  ,,       ,,       ,, small     ,,
  r_air(:,:,:) = 0.0_wp !log(1-rdrg) !log of the volume-fraction being not covered

! Provisional values for pattern lenth array:
! this should be a 2D external parameter field!
  DO j=jstartpar, jendpar
    DO i=istartpar, iendpar
      IF (fr_land(i,j) < 0.5_wp) THEN
        d_pat(i,j) = 0.0_wp
      ELSE
        d_pat(i,j) = pat_len 
        END IF
    END DO
  END DO

END SUBROUTINE init_volume_canopy

!==============================================================================
!==============================================================================
!+ Module procedure init_surface_canopy in for initialization
!+ of special external parameters describing the effective values of surface 
!+ areas indices.

SUBROUTINE init_surface_canopy (ie, je, itran_scheme,                      &
                                istartpar, iendpar, jstartpar, jendpar,    &
                                fr_land, plcov, lai, sai, tai, eai)

!------------------------------------------------------------------------------
!
! Description:
!
! Method:
!
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Declarations
!------------------------------------------------------------------------------

!Formal Parameters:

INTEGER (KIND=iintegers), INTENT(IN) :: &

! Horizontal and vertical sizes of the fields and related variables:
! --------------------------------------------------------------------

    ie,  & ! number of grid points in zonal      direction 
    je,  & ! number of grid points in meridional direction
    istartpar, iendpar, & ! zonal      start and end index including the model boundary lines
    jstartpar, jendpar, & ! meridional start and end index including the model boundary lines
    itran_scheme

! External parameter fields:
! --------------------------

REAL (KIND=wp),           INTENT(IN) :: &

    fr_land(ie,je)  ! land portion of a grid point area             ( 1 )

REAL (KIND=wp),           INTENT(INOUT) :: &

    plcov(ie,je), & ! fraction of plant cover                       ( 1 )
    lai  (ie,je), & ! leaf area index                               ( 1 )
    sai  (ie,je), & ! surface area index                            ( 1 )
    tai  (ie,je), & ! transpiration area index                      ( 1 )
    eai  (ie,je)    ! (evaporative) earth area index                ( 1 )

! ----------------
! Local variables:
! ----------------

  INTEGER (KIND=iintegers) ::  &
    i,  j         !  loop indices

  REAL    (KIND=wp)        ::  fakt

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine init_surface_canopy
!------------------------------------------------------------------------------

! Effective values of the surface area indices:

  DO j = jstartpar, jendpar
    DO i = istartpar, iendpar
      IF (fr_land(i,j) < 0.5_wp) THEN
        sai(i,j) = c_sea
      ELSE
        tai(i,j) = MAX( 1.0E-6_wp, lai(i,j) )
      END IF
    END DO
  END DO

  IF (itran_scheme.EQ.1) THEN
    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar
        IF (fr_land(i,j) >= 0.5_wp) THEN
          sai(i,j) = tai(i,j)
          eai(i,j) = (1.0_wp-plcov(i,j))*sai(i,j)
          tai(i,j) = plcov(i,j)*tai(i,j)
        ENDIF
      ENDDO
    ENDDO
  ELSE
    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar
        IF (fr_land(i,j) >= 0.5_wp) THEN
          tai(i,j) = plcov(i,j)*tai(i,j)  ! transpiration area index
          eai(i,j) = c_soil               ! evaporation area index
          sai(i,j) = c_lnd+tai(i,j)       ! surface area index

          ! effective area indeces by multiplication with the reduction factor fakt:
          fakt     = EXP (e_surf*LOG( sai(i,j)) )/sai(i,j)
          sai(i,j) = fakt * sai(i,j)
          eai(i,j) = fakt * eai(i,j)
          tai(i,j) = fakt * tai(i,j)
        ENDIF
      ENDDO
    ENDDO
  ENDIF

END SUBROUTINE init_surface_canopy

!==============================================================================
!==============================================================================
!+ Module procedure stab_funct in "src_turbdiff" for computing the stability
!+ function for vertical diffusion

SUBROUTINE stab_funct (sm, sh, fm2, fh2, frc, tvs, tls, dd, i_st,i_en, j_st,j_en)

   REAL (KIND=wp)            , INTENT(INOUT) :: &
        sm(:,:),  & !stablility function for momentum             [1] 
        sh(:,:)     !stablility function for scalars (heat)       [1] 

   REAL (KIND=wp)            , INTENT(IN)    :: &
        dd (:,:,0:)   ! local derived turbulence parameter

   REAL (KIND=wp)            , INTENT(IN) :: &
        fm2(:,:),   & !squared forcing frequency for momentum       [1/s2]
        fh2(:,:),   & !squared forcing frequency for scalars (heat) [1/s2]
        frc(:,:),   & !forcing function for TKE                     [m/s2]
        tvs(:,:),   & !turbulent velocity scale SQRT(2*TKE)         [m/s]
        tls(:,:)      !turbulent length scale

   INTEGER (KIND=iintegers), INTENT(IN) :: &
        i_st, i_en, & !start- and end index of horizontal i-loop
        j_st, j_en    !start- and end index of horizontal j-loop

   INTEGER (KIND=iintegers) :: i,j !loop indices

   REAL (KIND=wp)              :: &
        d1, d2, d3, d4 ,d5 ,d6

   REAL (KIND=wp)     :: &
        a11, a12, a22, a21, &
        a3, a5, a6, be1, be2, &
        gama, fakt, val1, val2, &
        gm, gh

!------------------------------------------------------------------------------

   DO j=j_st, j_en
   DO i=i_st, i_en

      d1=dd(i,j,1); d2=dd(i,j,2); d3=dd(i,j,3)
      d4=dd(i,j,4); d5=dd(i,j,5); d6=dd(i,j,6)

!US      gama=tls(i,j)*frc(i,j)/tvs(i,j)**2 !entspr. 1/d_mom im Gleichgewicht
!US                                         !und ausserh. des Bestandes
      gama = tls(i,j) * frc(i,j) / (tvs(i,j) * tvs(i,j))

!US      fakt=(tls(i,j)/tvs(i,j))**2 
      fakt = tls(i,j)/tvs(i,j)
      fakt = fakt*fakt

!     Folgnde Fallunterscheidung muss gemacht werden,
!     um positiv definite Loesungen fuer die Stabilitaets-
!     funktionen zu ermoeglichen:

      IF (fh2(i,j).GE.0.0_wp) THEN ! stab. Schichtung
!        Allgemeinste der hier verwendeten Loesungen:

         be1=1.0_wp
         be2=be1-c_g
       
         gh=fh2(i,j)*fakt
         gm=fm2(i,j)*fakt
         
         a11=d1+(d5-d4)*gh
         a12=d4*gm
         a21=(d6-d4)*gh
         a22=d2+d3*gh+d4*gm

         fakt=a11*a22-a12*a21
         sh(i,j)=(be1*a22-be2*a12)/fakt
         sm(i,j)=(be2*a11-be1*a21)/fakt
      ELSE ! labile Schichtung
!        Weiter eingeschraenkte Loesungen, bei denen Gm u.
!        Gh unter Einfuehrung der Ri-Zahl eliminiert werden.
!        Dabei wird gama aus der zuvor geloesten TKE-Gleich.
!        genommen. Wegen frc=tls*(sm*fm2-sh*fh2)
!        ist dies aber von den Vorgaengerwerten von sh u. sm
!        aghaengig:

!        Um physikalisch unsinnige Loesungen bez. Singulari-
!        taeten zu vermeiden, muessen be1 u. be2 > 0 sein:

         be1=1.0_wp-d4*gama
         be2=be1-c_g

!        Weil tvs schon vorher entsprechend nach unten beschraenkt wurde,
!        ist eine explizite Beschraenkung hier unnoetig!

         a3=d3*gama/d2
         a5=d5*gama/d1
         a6=d6*gama/d2
         be1=be1/d1
         be2=be2/d2

         val1=(fm2(i,j)*be2+(a5-a3+be1)*fh2(i,j))/(2.0_wp*be1)
!US      val2=val1+sqrt(val1**2-(a6+be2)*fh2(i,j)*fm2(i,j)/be1)
         val2=val1+sqrt(val1*val1-(a6+be2)*fh2(i,j)*fm2(i,j)/be1)
         fakt=fh2(i,j)/(val2-fh2(i,j))
         sh(i,j)=be1-a5*fakt
         sm(i,j)=sh(i,j)*(be2-a6*fakt)/(be1-(a5-a3)*fakt)
      END IF   

   END DO
   END DO

END SUBROUTINE stab_funct

!==============================================================================
!==============================================================================

SUBROUTINE turb_cloud (                      &
!
   ie, je, ke, ke1,            kcs,    kce,  &
   istart, iend, jstart, jend, kstart, kend, &
!
   prs, ps, rcld, t, qv, qc,                 &
!
   clc, clwc,                                &
   itype_wcld, lhocp, uc1, uc2, ucl,         &
   rdv, o_m_rdv, rvd_m_o, b1, b2w, b3, b4w, b234w)
   
!------------------------------------------------------------------------------
!
! Description:
!
!     This routine calculates the area fraction of a grid box covered
!     by stratiform (non-convective) clouds.
!     If subgrid-scale condensation is required, an additional
!     saturation adjustment is done.
!
! Method:
!
!     itype_wcld = 1 :
!     The fractional cloud cover clc is determined empirically from
!     relative humidity. Also, an in-cloud water content of sugrid-scale
!     clouds is determined as a fraction of the saturation specific
!     humidity. Both quantities define the grid-volume mean cloud water
!     content.
!     itype_wcld=2:
!     A Gaussion distribution is assumed for the saturation deficit
!     dq = qt - qs where qt = qv + ql is the total water content and
!     qs is the saturation specific humidity. Using the standard deviation
!     rcld of this distribution (on input) and the conservative grid-scale
!     quantities qt and tl (liquid water temperature), a corrected liquid
!     water content is determined which contains alse the contributions from
!     subgrid-scale clouds. A corresponding cloudiness is also calculated.
!
!------------------------------------------------------------------------------

! Subroutine arguments
!----------------------

! Scalar arguments with intent(in):

INTEGER (KIND=iintegers), INTENT(IN) :: &  ! indices used for allocation of arrays
  ie, je,            & ! number of grib points of horizontal dimensions
  ke, ke1,           & ! number of full and half level dimension
  kcs, kce,          & ! lower and upper level index bound for output variabls
  istart, iend, & ! zonal      start and end index
  jstart, jend, & ! meridional start and end index
  kstart, kend, & ! vertical   start and end index
  itype_wcld

! Array arguments with intent(in):

REAL (KIND=wp),     INTENT(IN) :: & !
  prs (ie,je,ke ),    &    ! base state pressure (")
  rcld(ie,je,ke1),    &    ! standard deviation of saturation deficit
  ps  (ie,je),        &    ! surface pressure
  t   (ie,je,ke ),    &    ! temperature (main levels)
  qv  (ie,je,ke )          ! water vapour (")

REAL (KIND=wp),     OPTIONAL, INTENT(IN) :: & !
  qc  (ie,je,ke )          ! cloud water  (")

! Array arguments with intent(out):

REAL (KIND=wp),     INTENT(OUT) :: &
  clc (ie,je,kcs:kce),  & ! stratiform subgrid-scale cloud cover
  clwc(ie,je,kcs:kce)     ! liquid water content of ""

REAL (KIND=wp),     INTENT(IN) :: & !
  uc1, uc2, ucl, lhocp,                            &
  rdv, o_m_rdv, rvd_m_o, b1, b2w, b3, b4w, b234w

! Local variables and constants
! -----------------------------

INTEGER (KIND=iintegers) :: &
  i,j,k ! loop indices

REAL (KIND=wp),     PARAMETER :: &
  zsig_max = 1.0E-3_wp,  & ! max. standard deviation of saturation deficit
  zclwfak  = 0.005_wp,   & ! fraction of saturation specific humidity
  zuc      = 0.95_wp       ! constant for critical relative humidity

REAL (KIND=wp)     :: &
  qt(ie,je), tl(ie,je), qs(ie,je), gam(ie,je) , &
  pres, ql, dq, q, sig, uc, & !
  zsigma, zclc1, zq_max, zx   !

REAL (KIND=wp)           :: &
  zpsat_w, zqvap, zdqsdt,       & !statement functions and
  zpvap, zqsat, ztemp, zpres      !their formal arguments

!------------ End of header ---------------------------------------------------

! Definition of statement functions:

! saturation vapour pressure over water (zpsat_w) and over ice (zpsat_i):
  zpsat_w(ztemp) = b1 * exp( b2w*(ztemp-b3)/(ztemp-b4w) )
! zpsat_i(ztemp) = b1 * exp( b2i*(ztemp-b3)/(ztemp-b4i) )

! specific humidity:
  zqvap(zpvap,zpres) = rdv * zpvap / ( zpres - o_m_rdv*zpvap )

! Derivation of zqsat with respect to temperature:
  zdqsdt(ztemp,zqsat) = b234w * ( 1.0_wp + rvd_m_o*zqsat ) * zqsat &
                             / (ztemp-b4w)**2

! Begin Subroutine turb_cloud
! ---------------------------

  zq_max = q_crit*(1.0_wp/clc_diag - 1.0_wp)

  DO k = kstart, kend

     IF (PRESENT(qc)) THEN
        DO j = jstart, jend
        DO i = istart, iend
           qt(i,j) = qc(i,j,k) +       qv(i,j,k) ! total water content
           tl(i,j) =  t(i,j,k) - lhocp*qc(i,j,k) ! liquid water temperature
        END DO 
        END DO 
     ELSE !qv and t already contain the conservation variables
        DO j = jstart, jend
        DO i = istart, iend
           qt(i,j) = qv(i,j,k)
           tl(i,j) =  t(i,j,k)
        END DO
        END DO
     END IF

     DO j = jstart, jend
     DO i = istart, iend
        qs(i,j) = zqvap( zpsat_w( tl(i,j) ), prs(i,j,k) )        ! saturation mixing ratio
       gam(i,j) = 1.0_wp / ( 1.0_wp + lhocp*zdqsdt( tl(i,j), qs(i,j) ) ) ! slope factor
     END DO
     END DO

     DO j = jstart, jend
     DO i = istart, iend

        pres = prs(i,j,k)        ! pressure
        dq   = qt(i,j) - qs(i,j) ! saturation deficit

        IF ( itype_wcld .EQ. 1 ) THEN

        ! Calculation of cloud cover and cloud water content
        ! using an empirical relative humidity criterion

          zsigma = pres / ps(i,j)

          ! critical relative humidity
          uc = zuc - uc1 * zsigma * ( 1.0_wp - zsigma )  &
                         * ( 1.0_wp + uc2*(zsigma-0.5_wp) )

          ! cloud cover
!US       clc(i,j,k) = MAX( 0.0_wp, &
!US                         MIN( 1.0_wp, clc_diag * ((qt(i,j)/qs(i,j)-uc)/(ucl-uc))**2 ) )
          zx = (qt(i,j)/qs(i,j) - uc) / (ucl-uc)
          clc(i,j,k) = MAX( 0.0_wp, MIN( 1.0_wp, clc_diag *  zx*zx ) )

          ! in-cloud water content
          ql = qs(i,j) * zclwfak

          ! grid-volume water content
          IF ( dq > 0.0_wp ) THEN
            zclc1 = clc_diag * ( (1.0_wp-uc)/(ucl-uc) )**2
            ql    = ql + (gam(i,j)*dq-ql)*(clc(i,j,k)-zclc1)/(1.0_wp-zclc1)
          END IF
          clwc(i,j,k) = clc(i,j,k) * ql

        ELSEIF ( itype_wcld .EQ. 2 ) THEN

        ! Statistical calculation of cloud cover and cloud water content
        ! using the standard deviation of the saturation deficit

          sig = MIN ( zsig_max, rcld(i,j,k) )

          ! in case of sig=0, the method is similar to grid-scale
          ! saturation adjustment. Otherwise, a fractional cloud cover
          ! is diagnosed.
          IF ( sig <= 0.0_wp ) THEN
            clc(i,j,k)  = ABS ( (SIGN(1.0_wp,dq)+1.0_wp)*0.5_wp )
            clwc(i,j,k) = clc(i,j,k) * gam(i,j) * dq
          ELSE
            q = dq/sig
            clc(i,j,k) = MIN ( 1.0_wp, MAX ( 0.0_wp, clc_diag * ( 1.0_wp + q/q_crit) ) )
            IF ( q <= - q_crit ) THEN
              clwc(i,j,k) = 0.0_wp
            ELSEIF ( q >= zq_max ) THEN
              clwc(i,j,k) = gam(i,j) * dq
            ELSE
              clwc(i,j,k) = gam(i,j) * sig * ( q + q_crit ) &
                                           * ( q + zq_max ) / ( 2.0_wp*( q_crit + zq_max) )
            ENDIF
          ENDIF

        ENDIF

     ENDDO
     ENDDO

  ENDDO

END SUBROUTINE turb_cloud

!==============================================================================
!==============================================================================
!+ Module procedure turb_param in "src_turbdiff" for computing some deduced parameters
!+ for turbulent transfer and - diffusion

SUBROUTINE turb_param (istart, iend, jstart, jend, grav, cp_d, dd)

   INTEGER (KIND=iintegers), INTENT(IN)    ::    &
     istart, iend, jstart, jend     ! start- and end-indices

   REAL (KIND=wp)          , INTENT(IN)    ::    &
     grav, cp_d

   REAL (KIND=wp)          , INTENT(INOUT) ::    &
     dd (:,:,0:)    ! local derived turbulence parameter


   INTEGER (KIND=iintegers) :: i,j !loop indices

!     Belegung abgeleiteter Konstanten:

!US      c_tke=d_mom**(1.0_wp/3.0_wp)
      c_tke = EXP (1.0_wp/3.0_wp * LOG(d_mom))
      c_g=1.0_wp-1.0_wp/(a_mom*c_tke)-6.0_wp*a_mom/d_mom !=3*0.08

      d_0=d_mom

      d_1=1.0_wp/a_heat
      d_2=1.0_wp/a_mom
      d_3=9.0_wp*a_heat
      d_4=6.0_wp*a_mom
      d_5=3.0_wp*(d_heat+d_4)
      d_6=d_3+3.0_wp*d_4

      rim=1.0_wp/(1.0_wp+(d_mom-d_4)/d_5)

      a_3=d_3/(d_2*d_mom)
      a_5=d_5/(d_1*d_mom)
      a_6=d_6/(d_2*d_mom)

      b_1=(1.0_wp-d_4/d_mom)/d_1
      b_2=(1.0_wp-d_4/d_mom-c_g)/d_2

      tet_g=grav/cp_d   ! g should be grav

!US   l_hori=r_earth/sqrt(eddlon*eddlat)    ! this is set before
      l_scal=MIN( l_hori, tur_len )
!print *,"in turb_param l_hori=",l_hori," l_scal=",l_scal

      DO j=jstart,jend
      DO i=istart,iend
         dd(i,j,0)=d_mom

         dd(i,j,1)=d_1
         dd(i,j,2)=d_2
         dd(i,j,3)=d_3
         dd(i,j,4)=d_4
         dd(i,j,5)=d_5
         dd(i,j,6)=d_6

         dd(i,j,7)=rim
      END DO
      END DO

END SUBROUTINE turb_param

!==============================================================================

END MODULE turbulence_utilities
