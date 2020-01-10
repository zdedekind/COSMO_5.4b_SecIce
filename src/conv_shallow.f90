!==============================================================================
!+ Source module  "conv_shallow"
!------------------------------------------------------------------------------

MODULE conv_shallow

!------------------------------------------------------------------------------
!
! Description:
!   The module "conv_shallow" performs calculations related to the
!   parameterization of subgrid-scale moist (shallow) convection for use in
!   high-resolution runs on the meso-gamma scale (convection permitting).
!   The present approach is based on a simple mass-flux scheme with moisture 
!   convergence closure.
!
!   All global fields of the model that are used by the convection routines
!   are passed via argument lists in blocked data format, while scalars are
!   imported by USE statements below. The interface of the convection
!   routines and the model is provided by the organizational routine
!   "conv_organize" in the module conv_interface.f90.
!
!   The blocked version of the convection schemes together with the interface
!   routine has been introduced in version 5.04b. The "history" before this
!   version refers to the module src_conv_tiedtke in ijk-data format.
!
! Current Code Owner: DWD, Dmitrii Mironov
!  phone:  +49  69  8062 2705
!  fax:    +49  69  8062 3721
!  email:  Dmitrii.Mironov@dwd.de
!
! History from former version src_conv_shallow:
! ---------------------------------------------
! 3.16       2005/07/22 Jochen Foerstner
!  Initial release
! 3.18       2006/03/03 Jochen Foerstner
!  Corrections for writing instantaneous values (Jochen Foerstner)
! V4_5         2008/09/10 Ulrich Schaettler
!  Moved declaration of entr_sc (before: entrscv) to new module data_convection
! V4_8         2009/02/16 Guenther Zaengl
!  Use p0hl (reference pressure at half levels) for full consistency with
!  new reference atmosphere implementation
! V4_12        2010/05/11 Ulrich Schaettler
!  Removed t0(_melt), qvsflx
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_18        2011/05/26 Martin Koehler
!  Replaced hard coded value by Namelist variable: thick_sc
! V4_20        2011/08/31 Matthias Raschendorfer
!  Introducing calculation of convective buoyant TKE production 'tket_conv'
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer
!  Replaced qx-variables by using them from the tracer module
! V4_26        2012/12/06 Matthias Raschendorfer
!  Allowing only non negative buoyant production terms of TKE
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  Introduced MESSy interface
! V4_28        2013/07/12 Ulrich Schaettler
!  Eliminated reference to sigmr (not used)
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
! V5_1         2014-11-28 Matthias Raschendorfer, Oliver Fuhrer, Michael Baldauf
!  Calculation of 'tket_conv' corrected.
!  Replaced ireals by wp (working precision) (OF)
!  Removed option lexpl_lbc=.FALSE.  (MB)
! V5_3         2015-10-09 Ulrich Blahak
!  Added computation of ttdiab_conv, the pure diabatic tendency due to convection
! VXXX         2016-05-02 Katherine Osterreid, Xavier Lapillonne
!  Adapt code to block structure. Add OpenACC statements
!  Local automatic fortran array have been replaced with allocatables
!  New allocate, deallocate routines have been introduced
!  Add y_conv_closure parameter for new closure of shallow convection
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V5_4b        2016-07-12 Katherine Osterried, Xavier Lapillonne, 
!                         Ulrich Schaettler
!  Initial release for blocked version of convection in COSMO:
!  Adapt code to block structure. Add OpenACC statements
!  Local automatic fortran array have been replaced with allocatables for OPENACC
!  Add y_conv_closure parameter for new closure of shallow convection
! @VERSION@    @DATE@     Ulrich Schaettler
!  Updated the acc present list
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

USE kind_parameters, ONLY :   &
    wp              ! KIND-type parameter for real variables

!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

    r_d,          & ! gas constant for dry air
    rdv,          & ! r_d / r_v
    rvd_m_o,      & ! r_v/r_d - 1
    cp_d,         & ! specific heat of dry air at constant pressure
    lh_v,         & ! latent heat of vapourization
    lh_s,         & ! latent heat of sublimation
    g,            & ! acceleration due to gravity
    b1,           & ! variables for computing the saturation vapour pressure
    b2w,          & ! over water (w) and ice (i)
    b2i,          & !               -- " --
    b3,           & !               -- " --
    b4w,          & !               -- " --
    b4i             !               -- " --

!------------------------------------------------------------------------------

USE conv_data, ONLY:                                                        &
    ctmelt               , & ! tripel point
    cc2                  , & ! definition of utitility constants
    c5hlccp              , & ! for saturation humidity
    c5hlscp              , & !
    chlcdcp              , & !
    chlsdcp              , & !
    lmfmid , entrmid, cmfcmax, cmfcmin, cmfctop

!==============================================================================


#ifdef _OPENACC
USE conv_data, ONLY:                                                        &
    mlab   , mtype  , ztu    , zqu    , zlu    , zlude  , zmfu   , zmfus  , &
    zmfuq  , zmful  , zqsen  , ztenh  , zqenh  , zqsenh , zrhgdz , z1dp   , &
    zrho   , zqold  , zentr  , zmfub  , zdqpbl , zdmfen , zdmfde , loflag
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure to compute shallow convection
!------------------------------------------------------------------------------

SUBROUTINE cu_shallow (                                          &
           pt     , pqv     ,                                    &
           pfif   , pfih    , ppf     , pph   ,                  &
           pt_g   , pshfl_s , plhfl_s , pdqvdt,                  &
           y_conv_closure   , thick_sc, entr_sc,                 &
           nvec   , kdim    , isc     , iec   ,                  &
           pcumflx, pcuconv , pco_clw , pmfu  ,                  &
           pdt_con,pdtdiab_con, pdqv_con, pdtke_con,             &
           mbas_con, mtop_con,  locum)

!------------------------------------------------------------------------------
!
! Description:
!
!   The module procedure cu_shallow organizes the massflux (shallow) cumulus
!   parameterization scheme.
!   This routine computes the physical tendencies of the prognostic variables
!   t and qv due to convective processes. These are calculated from
!     - convective fluxes due to updrafts
!     - condensation of cloud water in the updrafts
!     - evaporation of cloud water in the environment
!
!   Input for the scheme are the grid scale values of T, qv, qc, w, p, fi
!   and the moisture tendency due to turbulent mixing and 3-d advection
!
!   Output of the scheme are
!     - tendencies of T and qv
!     - cloud base and cloud top model level indices and hights
!
!
! Method:
!
!   The parameterisation is based on a mass flux representation of
!   convective processes and proceeds along the following steps:
!
!     - definition of constants and parameters
!     - specification of half level values and initialization of
!       updraft values
!     - determination of cloud base (if existing)
!       and specification of cloud base mass flux from PBL moisture
!       budget (lcape = .false.), or from convective available
!       energy (lcape = .true. )
!     - cloud ascent calculations (no downdrafts are considered)
!     - final adjustments to convective fluxes
!     - calculation of the tendencies of T and qv
!
! Externals
!
!   cu_asc  : cloud ascent for entraining plume
!   cu_cond : condensation/evaporation calculations
!
! Switches
!
!   lmfmid = .T.  MIDLEVEL CONVECTION IS SWITCHED ON
!
! Model parameters
!
!   Have been defined within the module
!
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

! Input data
! ----------
  INTEGER, INTENT (IN) ::  &
     nvec ,       & ! array dimension in zonal direction
     kdim ,       & ! array dimension in vertical direction
     isc  ,       & ! start index for first  array computation
     iec            ! end   index for first  array computation

  CHARACTER (LEN=10),       INTENT(IN) ::   &
    y_conv_closure  ! type of shallow convection closure

  REAL    (KIND=wp   ),     INTENT (IN) ::  &
     pt      (nvec,kdim),     & ! temperature at full levels
     pqv     (nvec,kdim),     & ! specific humidiy at full levels
     pfif    (nvec,kdim  ),   & ! geopotential at full levels
     pfih    (nvec,kdim+1),   & ! geopotential at half levels
     ppf     (nvec,kdim  ),   & ! pressure at full levels
     pph     (nvec,kdim+1),   & ! pressure at half levels
     pt_g    (nvec),          & ! weighted surface temperature
     pshfl_s (nvec),          & ! surface sensible heat flux in cu_base massflux
     plhfl_s (nvec),          & ! surface latent heat flux in cu_base massflux
     pdqvdt  (nvec,kdim)        ! moisture tencency

  REAL    (KIND=wp   ),     INTENT (IN) ::  &
     thick_sc,                & ! limit for convective clouds to be "shallow" (in Pa)
     entr_sc                    ! mean entrainment rate for shallow convection

! Output data
! -----------
  REAL    (KIND=wp   ),     INTENT (OUT) ::  &
     pdt_con     (nvec, kdim), & ! convective tendency of temperature
     pdtdiab_con (nvec, kdim), & ! conv T-tend. due to LH Exchanges
     pdqv_con    (nvec, kdim), & ! convective tendency of specific humidity
     pdtke_con   (nvec, kdim), & ! convective boyant TKE prpoduction on half levels
     pco_clw     (nvec, kdim), & ! convective cloud liquid water
     pmfu        (nvec, kdim), & ! convective mass flux, now has vertical dependence
     pcumflx     (nvec),       & ! cu_base massflux
     pcuconv     (nvec)          ! moisture convergence used for cu_base massflux
     ! These output-arrays have to be set to zero on call of this routine !

  INTEGER, INTENT (OUT)    ::  &
     mbas_con    (nvec)      , & ! cloud base level index
     mtop_con    (nvec)          ! cloud top  level index

  LOGICAL                 , INTENT (OUT) ::  &
     locum(nvec)        ! indicatior for convection at gridpoints

! Local scalars and automatic arrays (also to be used in lower level routines):
! ----------------------------------
  INTEGER ::  &
    ip, k,                 & ! loop indices over spatial dimensions
    mkb                      !


  REAL    (KIND=wp   )     ::  &
    zbuo, zc3, zc4, zzs, zqumqe, zdqmin, zpbmpt, zzp, zl,  & !
    zdmfmax, zdprho, zzentr, zqeen, zseen, zscde, zqude,   & !
    zmfusk, zmfuqk, zmfulk, zzdmf, invmasterlen, w_star,   &
    pa_frac  (nvec)          ! fractional area of core

  LOGICAL                  ::  &
    llo1, lfixedfrac         !

  REAL    (KIND=wp   ) ::  &
    zexn, zthlp, zthle, zthvp, zthve, zqsatenv, zqsatenvplus, zcor, &
    zqsatenvmin, zdqsat, zgamma, zalpha, zbeta, zchic, zqsatenvcb

#ifndef _OPENACC
  ! Declarations of working  arrays for the shallow convection
  INTEGER               :: &
    mlab   (nvec,kdim),    & !
    mtype  (nvec)            !

  REAL (KIND=wp)        :: &
    ztu    (nvec,kdim),    & !
    zqu    (nvec,kdim),    & !
    zlu    (nvec,kdim),    & !
    zlude  (nvec,kdim),    & !
    zmfu   (nvec,kdim),    & !
    zmfus  (nvec,kdim),    & !
    zmfuq  (nvec,kdim),    & !
    zmful  (nvec,kdim),    & !
    zqsen  (nvec,kdim),    & ! saturation specific humidiy at full levels
    ztenh  (nvec,kdim),    & !
    zqenh  (nvec,kdim),    & !
    zqsenh (nvec,kdim),    & !
    zrhgdz (nvec,kdim),    & ! rho*g*dz
    z1dp   (nvec)     ,    & !
    zrho   (nvec)     ,    & !
    zqold  (nvec)     ,    & !
    zentr  (nvec)     ,    & !
    zmfub  (nvec)     ,    & !
    zdqpbl (nvec)     ,    & !
    zdmfen (nvec)     ,    & !
    zdmfde (nvec)            !

  LOGICAL               :: &
    loflag (nvec)            !
#endif

!------------ End of header ---------------------------------------------------

  !$acc data                                                                  &
  !$acc present (pt, pqv, pfif, pfih, ppf, pph, pdqvdt, pt_g, pshfl_s)        &
  !$acc present (pdt_con,pdtdiab_con, pdqv_con, pdtke_con, plhfl_s)           &
  !$acc present (pcumflx, pcuconv, pco_clw, pmfu)                             &
  !$acc present (mbas_con, mtop_con, locum)                                   &
  !$acc present (mlab, mtype, ztu, zqu, zlu, zlude, zmfu, zmfus, zmfuq)       &
  !$acc present (zmful, zqsen, ztenh, zqenh, zqsenh, zrhgdz, z1dp, zrho)      &
  !$acc present (zqold, zentr, zmfub, zdqpbl, zdmfen, zdmfde, loflag)

!------------------------------------------------------------------------------
! Begin Subroutine cu_shallow
!------------------------------------------------------------------------------

  ! the highest level of these variables is not set below!
  pdt_con  (:,1)   = 0.0_wp
  pdqv_con (:,1)   = 0.0_wp
  pdtke_con(:,1)   = 0.0_wp
  pdtdiab_con(:,1) = 0.0_wp

!******************************************************************************
!*                                                                            *
!* Section 1: Initialisation of some variables and computation of additional  *
!*            half level properties ( former subroutine cu_ini )              *
!*                                                                            *
!******************************************************************************

  ! Initialization of arrays needed in the parameterisation of
  ! convection; in particular interpolation of large-scale
  ! (environmental) fields to model half levels and determination
  ! of the level of maximum (upward) vertical velocity

  ! Choose if you want the GRANT version of the modified shallow convection SB14
  lfixedfrac = .FALSE.

  ! 1. Specify saturation specific humidity, large scale parameter at half
  !    levels, adjust temperature if statically unstable, find level of maximum
  !    vertical velocity.
  !    Half level geopotential set to central value between full layers !
  !    First guess for half level temperature from adiabatic interpolation
  !    of adjacent full layer values (maximum of two values is taken)
  !    First guess for half level saturation humidity is saturation
  !    humidity at upper full level

  !$acc parallel
  DO k = 1, kdim
    !$acc loop gang vector
    DO ip = isc, iec
      IF( pt(ip,k) - ctmelt > 0.0_wp) THEN
        zc3    = b2w
        zc4    = b4w
      ELSE
        zc3    = b2i
        zc4    = b4i
      END IF
      zqsen(ip,k) = cc2*EXP( zc3*(pt(ip,k)-B3)/(pt(ip,k)-zc4) )/ppf(ip,k)
      zqsen(ip,k) = zqsen(ip,k) / ( 1.0_wp - rvd_m_o*zqsen(ip,k) )
    ENDDO
  ENDDO
  !$acc end parallel


  !$acc parallel
  DO k = 2, kdim
    !$acc loop gang vector
    DO ip = isc, iec
      ztenh(ip,k) = ( MAX (cp_d*pt(ip,k-1)+pfif(ip,k-1),cp_d*pt(ip,k)+pfif(ip,k) ) &
                    - pfih(ip,k) )/cp_d
      zqsenh(ip,k)= zqsen(ip,k-1)
      z1dp  (ip)  = 1.0_wp/pph(ip,k)
      loflag(ip)  = .TRUE.

      ! adjust temperature and humidity at half level to value for moist
      ! adiabatic interpolation from neighbour levels (ip.e. consider
      ! condensation/evaporation process)
      CALL cu_cond ( ztenh(ip,k), zqsenh(ip,k), z1dp(ip),        &
                     loflag(ip) , .TRUE.  , .TRUE. )

      ! interpolation of specific humidity to half levels, following
      ! a moist adiabate from the upper full level and avoiding supersaturation
      zqenh(ip,k) = MIN( pqv(ip,k-1), zqsen(ip,k-1) ) + zqsenh(ip,k) - zqsen(ip,k-1)
      zqenh(ip,k) = MAX( zqenh(ip,k), 0.0_wp )

    ENDDO
  END DO    ! vertical loop
  !$acc end parallel

  !$acc parallel
  !$acc loop gang vector
  DO ip = isc, iec
    ztenh (ip,1   ) = pt  (ip,1)
    zqenh (ip,1   ) = pqv (ip,1)
  END DO
  !$acc end parallel

  ! avoid unstable structure in half level temperature profile
  !$acc parallel
  !$acc loop seq
  DO k = kdim-1, 2, -1
    !$acc loop gang vector
    DO ip = isc, iec
      zzs = MAX( cp_d*ztenh(ip,k)+pfih(ip,k), cp_d*ztenh(ip,k+1)+pfih(ip,k+1) )
      ztenh(ip,k) = ( zzs - pfih(ip,k) )/cp_d  ! ztenh may become incompatible
    END DO                                    ! with zqenh
  END DO
  !$acc end parallel

  ! Initialize updraft and downdraft values
  !$acc parallel
  !$acc loop seq
  DO k = 1, kdim
    !$acc loop gang vector
    DO ip = isc, iec
      ztu    (ip,k) = ztenh (ip,k)
      zqu    (ip,k) = zqenh (ip,k)
      zlu    (ip,k) = 0.0_wp
      zmfu   (ip,k) = 0.0_wp
      zmfus  (ip,k) = 0.0_wp
      zmfuq  (ip,k) = 0.0_wp
      zlude  (ip,k) = 0.0_wp
      mlab   (ip,k) = 0
      zrhgdz (ip,k) = ppf(ip,k)*( pfih(ip,k) - pfih(ip,k+1) )/(r_d*pt(ip,k))
    END DO
  END DO
  !$acc end parallel

!******************************************************************************
!*                                                                            *
!* Section 2: Cloud base calculations                                         *
!*                                                                            *
!******************************************************************************

  ! a) Determination of cloud-base values
  !     Method
  !       humidity, pressure and geopotential at half levels this
  !       routine computes corresponding cloud base values and the
  !       following flag indices:
  !             mlab=1 as indicator of sub-cloud levels
  !             mlab=3 as indicator of unstable in-cloud levels
  !     - to define a cloud base, surface air is lifted dry-adiaba-
  !       tically up to cloud base (non-entraining plume, i.e. for
  !       constant massflux)
  !     - temperature and humidity are adjusted inside of the cloud
  !       to take into account condensation effects

  ! Initialize lifting level variables
  !$acc parallel
  !$acc loop gang vector
  DO ip = isc, iec
    mlab (ip,kdim ) = 1        ! lowest layer is below cloud base
    mbas_con(ip)    = kdim-1   ! cloud base index
    mtop_con(ip)    = kdim-1   ! cloud base index
    locum(ip)       =.false.
  END DO
  !$acc end parallel

  !$acc parallel
  !$acc loop seq
  DO k = kdim-1, 2, -1     ! Vertical loop over layers
    !$acc loop gang vector
    DO ip = isc, iec
      z1dp (ip) = 1.0_wp/pph(ip,k)
      IF ( mlab(ip,k+1).EQ.1 ) THEN
        loflag(ip) = .TRUE.   ! ascent continues
      ELSE
        loflag(ip) = .FALSE.  ! ascent complete
      ENDIF

      IF( loflag(ip) ) THEN        ! cloud base not found yet
        zqu(ip,k) = zqu(ip,k+1) ! retain parcel's humdity
                      ! parcel temperature after ascent to next layer
        ztu(ip,k) = ( cp_d*ztu(ip,k+1) + pfih(ip,k+1) - pfih(ip,k) )/cp_d

        ! difference between parcel (virtual) temperature and environmental
        ! (virtual) temperature determines buoancy of parcel (0.5K added)
        zbuo =   ztu  (ip,k) * ( 1.0_wp + rvd_m_o*zqu  (ip,k) )  &
               - ztenh(ip,k) * ( 1.0_wp + rvd_m_o*zqenh(ip,k) ) + 0.5_wp
        IF(zbuo > 0.0_wp) mlab(ip,k) = 1 ! sub-cloud indicator for positve buoyancy
        zqold(ip) = zqu(ip,k) ! store parcel humidity in local variable
      END IF

      !Check for condensation and adjust parcel temperature and humidity
      !if condensation/sublimation occurs

      CALL cu_cond ( ztu(ip,k), zqu(ip,k), z1dp(ip)  ,                 &
                     loflag(ip) , .TRUE.  , .FALSE.  )

!     If ascent calculations are still active and parcel humidity
!     changed due to condensation:
!
      IF( loflag(ip) .AND. zqu(ip,k).NE.zqold(ip) ) THEN
        mlab(ip,k) = 2 ! Index for lifting condensation level
        zlu(ip,k) = zlu(ip,k) + zqold(ip) - zqu(ip,k)
        zbuo =   ztu  (ip,k)*( 1.0_wp + rvd_m_o*zqu  (ip,k) )  &
               - ztenh(ip,k)*( 1.0_wp + rvd_m_o*zqenh(ip,k) ) + 0.5_wp
        IF(zbuo > 0.0_wp) THEN ! Test for buoyancy
          mbas_con(ip) = k        ! define cloud base index
          locum(ip)    = .TRUE.   ! indicate existence of unstable cloud base
          mlab(ip,k) = 3          ! index for unstable LCL
        END IF
      END IF
    END DO
  END DO
  !$acc end parallel


  ! b) total moisture convergence and decision on type of convection
  !$acc parallel
  !$acc loop gang vector
  DO ip = isc, iec
    zdqpbl(ip) = 0.0_wp
  END DO
  !$acc end parallel

  !$acc parallel
  !$acc loop seq
  DO k = 2, kdim
    !$acc loop gang vector
    DO ip = isc, iec
      IF ( k >= mbas_con(ip) )  THEN
        zdqpbl(ip) = zdqpbl(ip) + pdqvdt(ip,k)*zrhgdz(ip,k)
      ENDIF
    ENDDO
  ENDDO
  !$acc end parallel

  !$acc parallel
  !$acc loop gang vector
  DO ip = isc, iec
    mtype(ip) = 2     ! only shallow convection allowed here
  ENDDO
  !$acc end parallel

!******************************************************************************
!*                                                                            *
!* Section 3: Moisture supply in boundary layer and preliminary cloud base    *
!*            mass flux (excluding downdraft effects)                         *
!*                                                                            *
!******************************************************************************

! Moisture convergence closure
SELECT CASE (y_conv_closure)

CASE ("standard")
  !$acc parallel
  !$acc loop gang vector
  DO ip = isc, iec
    mkb = mbas_con(ip)
    zqumqe = zqu(ip,mkb) + zlu(ip,mkb) - zqenh(ip,mkb)
    zdqmin = MAX( 0.01_wp*zqenh(ip,mkb), 1.E-10_wp )
    llo1 = zdqpbl(ip) > (1.E-3_wp*g*zqumqe) & ! minimum positive moist. conv.
          .AND. zqumqe > zdqmin             & ! parcel humidity exceeds env.hum.
          .AND. locum(ip)                     ! convective grid point
    IF ( llo1 ) THEN
      ! SB14 OLD TIEDKE MOISTURE CONVERGENCE FORMULATION
      zmfub(ip) = zdqpbl(ip) / (g*MAX(zqumqe, zdqmin))
    ELSE
      zmfub(ip) = 0.0_wp
      locum(ip) = .FALSE.            ! GP switched off
    ENDIF
    zmfub(ip) = MIN( zmfub(ip), cmfcmax)   !  massflux upper limit
    pcumflx(ip)= zmfub(ip)
    IF (zdqpbl(ip).gt.0.0001_wp) THEN
      pcuconv(ip) = zdqpbl(ip)
    ELSE
      pcuconv(ip) = 0.0_wp
    END IF
    zentr(ip) = entr_sc
  ENDDO
  !$acc end parallel

CASE ("Boeing")
  !$acc parallel
  !$acc loop gang vector
  DO ip = isc, iec
    mkb = mbas_con(ip)
    zqumqe = zqu(ip,mkb) + zlu(ip,mkb) - zqenh(ip,mkb)
    zdqmin = MAX( 0.01_wp*zqenh(ip,mkb), 1.E-10_wp )
    ! avoid 'convective drizzle' by using a minimum moisture convergence
    ! llo1   = zdqpbl(i) > 0.0         & ! positive moisture convergence
    ! triggering criterium parcel humidity exceeds env.hum.
    llo1 = (zqumqe > zdqmin)
    pa_frac(ip) = 0.0_wp
    IF ( llo1 ) THEN
      ! SB14 CALCULATE w_star FROM SURFACE FLUXES
      w_star = MAX(0.0_wp,-pfih(ip,mkb)*(r_d*               &
                          pshfl_s(ip)/(cp_d*pph(ip,kdim+1))+&
                          rvd_m_o*r_d*pt_g(ip)*plhfl_s(ip)/ &
                          (lh_v*pph(ip,kdim+1)))            )**0.3333333_wp
      IF (lfixedfrac) THEN
        ! SB14 CLOSURE ACCORDING TO GRANT,
        ! PUTTING A SCALING FACTOR AT CLOUD BASE OF 0.03
        zmfub(ip) = 0.03_wp*w_star
      ELSE
        ! SB14 TIEDKE MOISTURE CONVERGENCE FORMULATION
        ! SB14 MORE COMPLICATED CLOSURE,
        ! WHICH LIFTS A SURFACE PARCEL TO ESTIMATE THE CLOUD FRACTION
        zqsatenv = cc2*EXP( b2w*(ztenh (ip,mkb)-b3)/(ztenh (ip,mkb)-b4w) )* &
                       (1.0_wp/pph(ip,mkb))
        zqsatenv = MIN( 0.5_wp, zqsatenv)
        zcor  = 1.0_wp / ( 1.0_wp - rvd_m_o*zqsatenv)
        zqsatenv = zqsatenv*zcor
        ! SB14 FOLLOWING VAN STRATUM ET AL 2014, JAS
        ! (FILL OUT THEIR EQN 7 FOR THE STD DEVIATION
        ! INTO CUIJPERS AND BECHTOLD 1995)
        ! PRECORRECTION FACTOR 0.30 FROM NEGGERS ET AL 2009
        pa_frac(ip)=0.15_wp+0.108_wp*                                       &
          ATAN(-1.55_wp*0.50_wp*(zqsatenv-zqenh(ip,mkb))/MAX(zqumqe, zdqmin))
        ! SB14 DIFFERENCE THOUGH, IS THAT THE UPDRAFT FRACTION IS CAPPED AT 15%
        ! SEEMS REASONABLE TO PREVENT LARGE UPDRAFT FRACTIONS IN STRATOCUMULUS
        pa_frac(ip)=MIN(0.15_wp,MAX(0.0_wp,pa_frac(ip)))
        zmfub(ip) = 0.84_wp*w_star*pa_frac(ip)
        IF(pa_frac(ip)<1.0e-6_wp) THEN
           pa_frac(ip) = 0.0_wp
           locum(ip) = .FALSE.            ! GP switched off
        ENDIF
      ENDIF ! lfixedfrac
      pcumflx(ip)= zmfub(ip)
    END IF
    IF (zdqpbl(ip).gt.0.0001_wp) THEN
      pcuconv(ip) = zdqpbl(ip)
    ELSE
      pcuconv(ip) = 0.0_wp
    END IF
    zentr(ip) = entr_sc
  ENDDO
  !$acc end parallel

END SELECT


!******************************************************************************
!*                                                                            *
!* Section 5: Cloud ascent for entraining plume                               *
!*                                                                            *
!******************************************************************************

! a) cloud ascent calcualation for an entraining plume to obtain the
!    vertical in-cloud profiles of T, q and l
!    Surface air is lifted dry-adiabatically to cloud base from there, a moist
!    ascent is calculated for an entraining/detraining plume;
!    cases where shallow convection does not occur, a check for possible
!    existence of mid-level (shallow) convection is made

  ! Security parameter
  zdmfmax = cmfcmax /( kdim/4 )

  !$acc parallel
  !$acc loop gang vector
  DO ip = isc, iec
    IF( .NOT.locum(ip) ) mtype(ip) = 0   ! type of convection
  ENDDO
  !$acc end parallel

  !$acc parallel
  DO k = 1, kdim
    !$acc loop gang vector
    DO ip = isc, iec
      pmfu  (ip,k) = 0.0_wp
      zmfus (ip,k) = 0.0_wp
      zmfuq (ip,k) = 0.0_wp
      zlude  (ip,k) = 0.0_wp
    ENDDO
  ENDDO
  !$acc end parallel

  !$acc parallel
  !$acc loop gang vector
  DO ip = isc, iec
    IF (.NOT.locum(ip) ) THEN
      DO k = 1, kdim
        mlab(ip,k) = 0
      ENDDO
      mbas_con(ip) = kdim -1
      mtop_con(ip) = kdim -1
      zmfub(ip)     = 0.0_wp
    ENDIF
  ENDDO
  !$acc end parallel

! Initialize values at lifting level
  !$acc parallel
  !$acc loop gang vector
  DO ip = isc, iec
    IF ( locum(ip) ) THEN
      mkb = mbas_con(ip)
      pmfu (ip,mkb) = zmfub(ip)
      zmfus(ip,mkb) = zmfub(ip)*( cp_d*ztu(ip,mkb) + pfih(ip,mkb) ) ! SB14: static energy
      zmfuq(ip,mkb) = zmfub(ip)*zqu(ip,mkb)                         ! SB14: total water
      zmful(ip,mkb) = zmfub(ip)*zlu(ip,mkb)                         ! SB14: liquid energy
      mtop_con(ip)  = mbas_con(ip)
    ENDIF
  ENDDO
  !$acc end parallel

! Perform ascent:  a) dry adiabatic lifting and
!                  b) allowing for condensation
!                  c) check buoancy and set flags
!                     (mlab=1 :sub-cloud layer, mlab=2 :cloud layer)

SELECT CASE (y_conv_closure)

CASE ("standard")
  !$acc parallel
  !$acc loop seq
  DO k = kdim-1, 2, -1
    !$acc loop gang vector
    DO ip = isc, iec
      IF( mlab(ip,k+1) == 0 ) mlab(ip,k) = 0
      IF( mlab(ip,k+1) == 3 ) THEN
        loflag(ip)= .TRUE.   ! active (unstable) grid point  in layer below
      ELSE
        loflag(ip)= .FALSE.  ! inactive grid point
      ENDIF

    ! calculation of entrainement/detrainement rates
      zdmfen(ip) = 0.0_wp
      zdmfde(ip) = 0.0_wp
      zrho  (ip) = pph(ip,k+1)/(r_d*ztenh(ip,k+1))
      z1dp  (ip) = 1.0_wp/pph(ip,k)

      IF( loflag(ip) ) THEN
        zdprho = ( pfih(ip,k) - pfih(ip,k+1) ) / g
        zzentr = zentr(ip)*pmfu(ip,k+1)*zdprho
        IF( k < mbas_con(ip) ) THEN
            zdmfde(ip) = zzentr
            zdmfen(ip) = zzentr
        ENDIF
      ENDIF

    ! adiabatic ascent for entraining/detraining plume

      IF(loflag(ip)) THEN
        zdmfen(ip) = MIN( zdmfen(ip), zdmfmax )
        zdmfde(ip) = MIN( zdmfde(ip), 0.75_wp*pmfu(ip,k+1) )
        pmfu(ip,k)= pmfu(ip,k+1) + zdmfen(ip) - zdmfde(ip)
        zqeen = zqenh(ip,k+1)*zdmfen(ip)
        zseen = ( cp_d*ztenh (ip,k+1) + pfih(ip,k+1) )*zdmfen(ip)
        zscde = ( cp_d*ztu(ip,k+1) + pfih(ip,k+1) )*zdmfde(ip)
        zqude = zqu(ip,k+1)*zdmfde(ip)
        zlude(ip,k) = zlu(ip,k+1)*zdmfde(ip)
        zmfusk = zmfus(ip,k+1) + zseen - zscde
        zmfuqk = zmfuq(ip,k+1) + zqeen - zqude
        zmfulk = zmful(ip,k+1) - zlude(ip,k)
        zlu  (ip,k) =  zmfulk*( 1.0_wp/MAX(cmfcmin,pmfu(ip,k)) )
        zqu(ip,k) =  zmfuqk*( 1.0_wp/MAX(cmfcmin,pmfu(ip,k)) )
        ztu (ip,k) = (zmfusk*( 1.0_wp/MAX(cmfcmin,pmfu(ip,k)) ) -      &
                                  pfih(ip,k))/cp_d
        ztu(ip,k)  = MAX( 100.0_wp, ztu(ip,k) )
        ztu(ip,k)  = MIN( 400.0_wp, ztu(ip,k) )
        zqold(ip)  = zqu(ip,k)
      ENDIF

      ! corrections for moist ascent by adjusting T,q and l
      ! calulation of condensation and corresponding adjustment of T and q

      CALL cu_cond ( ztu(ip,k), zqu(ip,k), z1dp(ip)  ,              &
                     loflag(ip) , .TRUE.  , .FALSE. )

      IF( loflag(ip) .AND.  zqu(ip,k).NE.zqold(ip) ) THEN
        zlu (ip,k) = zlu(ip,k) + zqold(ip) - zqu(ip,k)
        zbuo      =  ztu(ip,k)*(1.0_wp+rvd_m_o*zqu(ip,k))  &
                    - ztenh(ip,k)*(1.0_wp+rvd_m_o*zqenh(ip,k))
        IF(mlab(ip,k+1)==1) zbuo = zbuo + 0.5_wp
        IF( zbuo > 0.0_wp .AND. pmfu(ip,k) >= 0.1_wp*zmfub(ip) ) THEN
          mtop_con(ip) = k
          mlab(ip,k)   = 3
        ELSE
          mlab(ip,k)  = 0
          pmfu(ip,k) = 0.0_wp
        ENDIF
      ENDIF

      IF( loflag(ip) ) THEN
        zmful(ip,k) = zlu(ip,k)*pmfu(ip,k)
        zmfus(ip,k) = ( cp_d*ztu(ip,k) + pfih(ip,k) )*pmfu(ip,k)
        zmfuq(ip,k) = zqu(ip,k)*pmfu(ip,k)
      ENDIF
    ENDDO
  ENDDO         ! vertical loop
  !$acc end parallel

CASE ("Boeing")
  !$acc parallel
  !$acc loop seq
  DO k = kdim-1, 2, -1
    !$acc loop gang vector
    DO ip = isc, iec
      IF( mlab(ip,k+1) == 0 ) mlab(ip,k) = 0
      IF( mlab(ip,k+1) == 3 ) THEN
        loflag(ip)= .TRUE.   ! active (unstable) grid point  in layer below
      ELSE
        loflag(ip)= .FALSE.  ! inactive grid point
      ENDIF

    ! calculation of entrainement/detrainement rates
      zdmfen(ip) = 0.0_wp
      zdmfde(ip) = 0.0_wp
      zrho  (ip) = pph(ip,k+1)/(r_d*ztenh(ip,k+1))
      z1dp  (ip) = 1.0_wp/pph(ip,k)

      IF( loflag(ip) ) THEN
        zdprho = ( pfih(ip,k) - pfih(ip,k+1) ) / g
        zzentr = zentr(ip)*pmfu(ip,k+1)*zdprho
        IF( k < mbas_con(ip) ) THEN
          ! calculate chi_critical at level below (k+1)
          zexn = (pph(ip,k+1)/1.0e5_wp)**(r_d/cp_d)
          zthlp = (ztu (ip,k+1)-(lh_v/cp_d)*zlu(ip,k+1))/zexn
          zthle = (ztenh (ip,k+1))/zexn
          zthvp = ztu(ip,k+1)*(1.0_wp+rvd_m_o*zqu(ip,k+1)-zlu(ip,k+1))/zexn
          zthve = ztenh(ip,k+1)*(1.0_wp+rvd_m_o*zqenh(ip,k+1))/zexn
          zqsatenvplus = cc2 * EXP(b2w*(ztenh(ip,k+1)+0.001_wp-b3)/ &
                                  (ztenh(ip,k+1)+0.001_wp-b4w))   * &
                               (1.0_wp/pph(ip,k+1))
          zqsatenvplus = MIN( 0.5_wp, zqsatenvplus )
          zcor  = 1.0_wp / ( 1.0_wp - rvd_m_o*zqsatenvplus )
          zqsatenvplus = zqsatenvplus*zcor
          zqsatenvmin = cc2 * EXP(b2w*(ztenh(ip,k+1)-0.001_wp-b3)/ &
                                 (ztenh(ip,k+1)-0.001_wp-b4w))   * &
                              (1.0_wp/pph(ip,k+1))
          zqsatenvmin = MIN( 0.5_wp, zqsatenvmin )
          zcor  = 1.0_wp / ( 1.0_wp - rvd_m_o*zqsatenvmin )
          zqsatenvmin = zqsatenvmin*zcor
          zqsatenv = cc2 * EXP(b2w*(ztenh(ip,k+1)-b3)/ &
                              (ztenh(ip,k+1)-b4w))   * &
                           (1.0_wp/pph(ip,k+1))
          zqsatenv = MIN( 0.5_wp, zqsatenv)
          zcor  = 1.0_wp / ( 1.0_wp - rvd_m_o*zqsatenv)
          zqsatenv = zqsatenv*zcor
          zdqsat = (zqsatenvplus-zqsatenvmin)/0.002_wp
          zgamma = (lh_v/cp_d)*zdqsat
          zalpha = zexn*zthlp/(lh_v/cp_d)
          zbeta = (1.0_wp+(1.0_wp+rvd_m_o)*zgamma*zalpha)/(1.0_wp+zgamma)
          zchic= (zthvp-zthve)/                                          &
                 (zbeta*(zthlp-zthle) + (zbeta-zalpha)*                  &
                 (lh_v/cp_d)*(zqu(ip,k+1)+zlu(ip,k+1)-zqenh(ip,k+1))/zexn)
          ! SB 14 hard coded entrainment: scales as 1/z above 100 m above
          ! surface for standard entrainment
          ! entr_sc = 0.00030_wp, hence zzentr*3333.3_wp=1
          ! g/0.5_wp*(pfih(i,k)+pfih(i,k+1)-(pfih(i,ke)+pfih(i,ke))) is the
          ! height above the surface calculated from the geopotential difference
          ! Inverse master length scale
          invmasterlen = g/MAX(0.5_wp*(pfih(ip,k)+pfih(ip,k+1)-  &
                         (pfih(ip,kdim)+pfih(ip,kdim))),1000.0_wp)
          zdmfen(ip) = (zzentr*3333.3_wp)*invmasterlen
          !SB 14: old formulation
          ! zdmfen(i) = zzentr*
          !      (1.3_wp-zqenh(i,k+1)/(zqsatenv))*((zqsatenv)/(zqsatenvcb))**3
          IF( (k < (mbas_con(ip)-1)) .or. ((zthvp-zthve)>0.0_wp)) THEN
            ! SB 14 Following Boing 2012, with prefactors increased to account
            ! for shallow convection. The prefactor here follows the relation
            ! between z_lcl and cloud area in Stirling and Stratton (2012,
            ! actually they have congestus)
            ! Estimate of typical cloud radius
            invmasterlen=invmasterlen/sqrt(0.09_wp/3.141592_wp)
            zdmfde(ip)  = zdmfen(ip)+pmfu(ip,k+1)*zdprho*               &
              MAX(invmasterlen-3.0*invmasterlen*MAX(zchic,0.0_wp),0.0_wp)
          ELSE
            zdmfde(ip)=0.0_wp
          ENDIF
        ENDIF
      ENDIF

    ! adiabatic ascent for entraining/detraining plume

      IF(loflag(ip)) THEN
        zdmfen(ip) = MIN( zdmfen(ip), zdmfmax )
        zdmfde(ip) = MIN( zdmfde(ip), 0.75_wp*pmfu(ip,k+1) )
        pmfu(ip,k)= pmfu(ip,k+1) + zdmfen(ip) - zdmfde(ip)
        zqeen = zqenh(ip,k+1)*zdmfen(ip)
        zseen = ( cp_d*ztenh (ip,k+1) + pfih(ip,k+1) )*zdmfen(ip)
        zscde = ( cp_d*ztu(ip,k+1) + pfih(ip,k+1) )*zdmfde(ip)
        zqude = zqu(ip,k+1)*zdmfde(ip)
        zlude(ip,k) = zlu(ip,k+1)*zdmfde(ip)
        zmfusk = zmfus(ip,k+1) + zseen - zscde
        zmfuqk = zmfuq(ip,k+1) + zqeen - zqude
        zmfulk = zmful(ip,k+1) - zlude(ip,k)
        zlu  (ip,k) =  zmfulk*( 1.0_wp/MAX(cmfcmin,pmfu(ip,k)) )
        zqu(ip,k) =  zmfuqk*( 1.0_wp/MAX(cmfcmin,pmfu(ip,k)) )
        ztu (ip,k) = (zmfusk*( 1.0_wp/MAX(cmfcmin,pmfu(ip,k)) ) -      &
                                  pfih(ip,k))/cp_d
        ztu(ip,k)  = MAX( 100.0_wp, ztu(ip,k) )
        ztu(ip,k)  = MIN( 400.0_wp, ztu(ip,k) )
        zqold(ip)  = zqu(ip,k)
      ENDIF

      ! corrections for moist ascent by adjusting T,q and l
      ! calulation of condensation and corresponding adjustment of T and q

      CALL cu_cond ( ztu(ip,k), zqu(ip,k), z1dp(ip)  ,              &
                     loflag(ip) , .TRUE.  , .FALSE. )

      IF( loflag(ip) .AND.  zqu(ip,k).NE.zqold(ip) ) THEN
        zlu (ip,k) = zlu(ip,k) + zqold(ip) - zqu(ip,k)
        zbuo      =  ztu(ip,k)*(1.0_wp+rvd_m_o*zqu(ip,k))  &
                    - ztenh(ip,k)*(1.0_wp+rvd_m_o*zqenh(ip,k))
        IF(mlab(ip,k+1)==1) zbuo = zbuo + 0.5_wp
        IF( zbuo > 0.0_wp .AND. pmfu(ip,k) >= 0.1_wp*zmfub(ip) ) THEN
          mtop_con(ip) = k
          mlab(ip,k)   = 3
        ELSE
          mlab(ip,k)  = 0
          pmfu(ip,k) = 0.0_wp
        ENDIF
      ENDIF

      IF( loflag(ip) ) THEN
        zmful(ip,k) = zlu(ip,k)*pmfu(ip,k)
        zmfus(ip,k) = ( cp_d*ztu(ip,k) + pfih(ip,k) )*pmfu(ip,k)
        zmfuq(ip,k) = zqu(ip,k)*pmfu(ip,k)
      ENDIF
    ENDDO
  ENDDO         ! vertical loop
  !$acc end parallel

END SELECT


    ! convective fluxes above non-buoancy level
    !
    !   (NOTE: CLOUD VARIABLES LIKE T,Q and L ARE NOT
    !          AFFECTED BY DETRAINMENT and ARE ALREADY KNOWN
    !          FROM PREVIOUS CALCULATIONS ABOVE)

  !$acc parallel
  !$acc loop gang vector
  DO ip = isc, iec
    IF( mtop_con(ip) == kdim-1 ) locum(ip) = .FALSE.
    mbas_con(ip) = MAX( mbas_con(ip), mtop_con(ip) )
  ENDDO
  !$acc end parallel

  !$acc parallel
  !$acc loop gang vector
  DO ip = isc, iec
    IF( locum(ip) ) THEN
      k          = mtop_con(ip) - 1
      zzdmf      = cmfctop
      zdmfde(ip)  = ( 1.0_wp - zzdmf )*pmfu(ip,k+1)
      zlude(ip,k) = zdmfde(ip)*zlu(ip,k+1)
      pmfu(ip,k) = pmfu(ip,k+1) - zdmfde(ip)
      zmfus(ip,k)  = ( cp_d*ztu(ip,k) + pfih(ip,k) )*pmfu(ip,k)
      zmfuq(ip,k)  = zqu(ip,k)*pmfu(ip,k)
      zmful(ip,k)  = zlu(ip,k)*pmfu(ip,k)
      zlude(ip,k-1) = zmful(ip,k)
    ENDIF
  ENDDO
  !$acc end parallel


  ! Additional checks
  ! Check cloud depth and delete deep convection
  !$acc parallel
  !$acc loop gang vector
  DO ip = isc, iec
    IF (locum(ip)) THEN
      ! cloud thickness (in Pa)
      zpbmpt = pph(ip,mbas_con(ip)) - pph(ip,mtop_con(ip))
      IF (zpbmpt > thick_sc) THEN
        locum(ip) = .FALSE.
        mbas_con(ip) = 0
        mtop_con(ip) = 0
      ENDIF
    ENDIF
  ENDDO
  !$acc end parallel

  ! Store in-cloud water content
  !$acc parallel
  DO k = 1, kdim
    !$acc loop gang vector
    DO ip = isc, iec
      IF (locum(ip)) THEN
        pco_clw(ip,k) = MAX ( 0.0_wp, MIN(100.0_wp, zlu(ip,k)) )
      ENDIF
    ENDDO
  ENDDO
  !$acc end parallel


!******************************************************************************
!*                                                                            *
!* Section 7: Determine final convective fluxes (former subroutine cu_flux)   *
!*            and conside evaporation of rain in sub-cloud layers             *
!*                                                                            *
!******************************************************************************

  !$acc parallel
  !$acc loop gang vector
  DO ip = isc, iec
    IF( .NOT.locum(ip))  mtype(ip) = 0
  ENDDO
  !$acc end parallel

  !$acc parallel
  DO k = 2, kdim
    !$acc loop gang vector
    DO ip = isc, iec
      IF( locum(ip) .AND. k >= mtop_con(ip)-1 ) THEN
         zmfus(ip,k) = zmfus(ip,k) - pmfu (ip,k)*(cp_d*ztenh(ip,k)+pfih(ip,k))
         zmfuq(ip,k) = zmfuq(ip,k) - pmfu (ip,k)*zqenh(ip,k)
      ELSE              ! above cloud and non-convective points
         pmfu  (ip,k  ) = 0.0_wp
         zmfus (ip,k  ) = 0.0_wp
         zmfuq (ip,k  ) = 0.0_wp
         zmful (ip,k  ) = 0.0_wp
         zlude (ip,k-1) = 0.0_wp
      ENDIF
    ENDDO
  ENDDO
  !$acc end parallel

  !$acc parallel
  !$acc loop seq
  DO k = 2, kdim
    !$acc loop gang vector
    DO ip = isc, iec
      IF( locum(ip) .AND. k > mbas_con(ip) ) THEN
        mkb = mbas_con(ip)
        zzp = ( pph(ip,kdim+1) - pph(ip,k) )/( pph(ip,kdim+1) - pph(ip,mkb) )
        IF ( mtype(ip).EQ.3 ) zzp=zzp**2
        zmfus(ip,k) = zmfus(ip,mkb) * zzp
        zmfuq(ip,k) = zmfuq(ip,mkb) * zzp
        zmful(ip,k) = zmful(ip,mkb) * zzp
      ENDIF
    ENDDO
  ENDDO
  !$acc end parallel

!******************************************************************************
!*                                                                            *
!* Section 8: Compute the final tendencies for grid scale variables T and qv  *
!*            (former subroutine cu_dtdq)                                     *
!*                                                                            *
!******************************************************************************

  !$acc parallel
  DO k = 2, kdim - 1  ! above lowest model layer
    !$acc loop gang vector
    DO ip = isc, iec
      IF(locum(ip)) THEN ! for convective grid points only
        llo1 = (pt(ip,k)-ctmelt) > 0.0_wp  !
        IF ( llo1 ) THEN
          zl = lh_v
        ELSE
          zl = lh_s
        ENDIF
        pdt_con(ip,k) =    g/zrhgdz(ip,k)/cp_d  &
                      * (  zmfus (ip,k+1) - zmfus (ip,k) &
                         - zl*( zmful(ip,k+1)-zmful(ip,k) ) )

        pdtdiab_con(ip,k) =    g/zrhgdz(ip,k)/cp_d  &
                      * ( - zl*( zmful(ip,k+1)-zmful(ip,k) ) )

        pdqv_con(ip,k)=   g/zrhgdz(ip,k)                   &
                       *(  zmfuq(ip,k+1) - zmfuq(ip,k)   &
                         + zmful(ip,k+1) - zmful(ip,k)   )
      END IF
    END DO
  END DO
  !$acc end parallel

  !$acc parallel
  !$acc loop gang vector
  DO ip = isc, iec
    IF(locum(ip)) THEN      ! for convective grid points only
      llo1 = (pt(ip,kdim)-ctmelt) > 0.0_wp
      IF ( llo1 ) THEN
         zl = lh_v
      ELSE
         zl = lh_s
      ENDIF
      pdt_con(ip,kdim) = - g /zrhgdz(ip,kdim)/cp_d * ( zmfus(ip,kdim) -zl*zmful(ip,kdim) )
      pdtdiab_con(ip,k) = - g /zrhgdz(ip,k)/cp_d * ( -zl*zmful(ip,k) )
      pdqv_con(ip,kdim)= - g /zrhgdz(ip,kdim) * ( zmfuq(ip,kdim)  + zmful(ip,kdim) )
    END IF
  END DO
  !$acc end parallel

  !$acc parallel
  DO k = 2, kdim ! for all model half levels (except the lower model boundary)
    !$acc loop gang vector
    DO ip = isc, iec
      IF (locum(ip)) THEN ! for convective grid points only
        ! RA: Computation of the convective buoyant TKE production:
        pdtke_con(ip,k) = MAX( 0.0_wp, g*r_d/pph(ip,k) *    &
                         ( (1.0_wp+rvd_m_o*zqenh(ip,k))*zmfus(ip,k)/cp_d + &
                           rvd_m_o*ztenh(ip,k)*zmfuq(ip,k) ) )
        ! RA
      END IF
    END DO
  END DO
  !$acc end parallel

  !$acc end data

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------
END SUBROUTINE cu_shallow

!==============================================================================
!+ Module procedure in "Convection"
!------------------------------------------------------------------------------

SUBROUTINE cu_cond (  pt , pqv , p1op , plflag , pcflag , peflag )

!-----------------------------------------------------------------------------!
! Description:
!
!   The module procedure cu_cond does a saturation adjustment for
!   temperature and specific humidity.
!
!   Method:    Thermodynamic adjustment by instantaneous condensation
!              at constant pressure using a double iteration method.
!              Release of latent heat of condendation and of deposition
!              is considered depending on temperature.
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------
! Subroutine arguments:
! --------------------
! Input data
! ----------

  REAL     (KIND=wp   ),     INTENT (IN) ::  &
     p1op   ! reciprocal of pressure, 1.0/p

  LOGICAL                 ,  INTENT (IN) ::  &
     plflag,& ! switch for points where adjustment shall be made
     pcflag,       & ! condensation only (.TRUE)
     peflag          ! evaporation only  (.TRUE)

! Input/Output data
! -----------
  REAL     (KIND=wp   ),     INTENT (INOUT) ::  &
     pt    , & ! temperature on input, adjusted on output
     pqv       ! specific humidity on input, adjusted on ouput

! Local scalars and automatic arrays:
! ----------------------------------
  INTEGER ::  &
    ip               ! loop indix

  REAL    (KIND=wp   )     ::  &
    zcond            ! condensation amount

  REAL    (KIND=wp   )     ::  &
    zc3, zc4, zc5, zhldcp, zqsat, zcor, zcond1, zfacc, zface  ! local storage

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine cu_cond
!------------------------------------------------------------------------------
  zfacc = 0.0_wp
  zface = 0.0_wp
  IF (pcflag) zfacc = 1.0_wp
  IF (peflag) zface = 1.0_wp

  zcond = 0.0_wp     ! Initialize condensation variable
  IF(plflag) THEN        ! only, if ascent still continues
    IF( pt - ctmelt > 0.0_wp) THEN     ! condensation
      zc3    = b2w
      zc4    = b4w
      zc5    = c5hlccp
      zhldcp = chlcdcp
    ELSE                                       ! deposition
      zc3    = b2i
      zc4    = b4i
      zc5    = c5hlscp
      zhldcp = chlsdcp
    END IF
    zqsat = cc2*EXP( zc3*(pt-b3)/(pt-zc4) )*p1op
    zqsat = MIN( 0.5_wp, zqsat )
    zcor  = 1.0_wp / ( 1.0_wp - rvd_m_o*zqsat )
    zqsat = zqsat*zcor
    zcond1   = (pqv-zqsat)/(1.0_wp+zc5*zqsat*zcor/(pt-zc4)**2)
    zcond  = zfacc*MAX( zcond1, 0.0_wp )  + &
               zface*MIN( zcond1, 0.0_wp )
    pt    = pt + zhldcp*zcond
    pqv   = pqv - zcond
  END IF

  !Second iteration
  IF( plflag .AND. zcond.NE.0.0_wp) THEN  !saturation adjustment
    IF( pt - ctmelt > 0.0_wp ) THEN
      zc3    = b2w
      zc4    = b4w
      zc5    = c5hlccp
      zhldcp = chlcdcp
    ELSE
      zc3    = b2i
      zc4    = b4i
      zc5    = c5hlscp
      zhldcp = chlsdcp
    END IF
    zqsat = cc2*EXP( zc3*(pt-b3)/(pt-zc4) )*p1op
    zqsat = MIN( 0.5_wp, zqsat )
    zcor  = 1.0_wp / ( 1.0_wp - rvd_m_o*zqsat )
    zqsat = zqsat * zcor
    zcond1=(pqv-zqsat)/(1.0_wp+zc5*zqsat*zcor/(pt-zc4)**2)
    pt   = pt + zhldcp*zcond1
    pqv  = pqv - zcond1
  END IF

!-----------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------
END SUBROUTINE cu_cond

!==============================================================================

!------------------------------------------------------------------------------
! End of module conv_shallow
!------------------------------------------------------------------------------

END MODULE conv_shallow
