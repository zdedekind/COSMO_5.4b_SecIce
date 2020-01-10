!+ Dummy interfaces for RTTOV library routines
!------------------------------------------------------------------------------
!
! Description:
!   This file provides dummy interfaces for the calls to the RTTOV-library.
!   This library, which was written by UKMO et.al., is licensed by the 
!   ESA NWP-SAF. If the library is not available, these dummy interfaces 
!   have to be compiled and linked with the LM in order to avoid "unresolved
!   symbols" error while linking the LM.
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 3.7        2004/02/18 Ulrich Schaettler
!  Initial release
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!==============================================================================

SUBROUTINE rttov7_mpi_settings (nproc_app, myproc_app, ncomm_app,      &
                                imp_real_app, ierr)

!------------------------------------------------------------------------------

INTEGER, PARAMETER       ::                                       &
   ireals    = SELECTED_REAL_KIND (12,200),                       &
   iintegers = KIND  (1)

INTEGER(KIND=iintegers), INTENT(IN)  ::       &
   nproc_app,    & ! number of processors in communicator of application
   myproc_app,   & ! id of this processor in communicator of application
   ncomm_app,    & ! communicator of application
   imp_real_app    ! Real-datatype used in application

INTEGER(KIND=iintegers), INTENT(OUT) ::       &
   ierr            ! error-status variable

!------------------------------------------------------------------------------

ierr = -99999

END SUBROUTINE rttov7_mpi_settings

!==============================================================================
!==============================================================================

SUBROUTINE RTTOV7_RTTVI(                                   &
     kerr, kppf, kpnsat, kplev, kpch, kpchus,              &
     kpnav, kpnsav, kpnssv, kpncv,                         &
     nrttovid, platform, satellite, instrument , numchans, &
     preslev, otmin, otmax, oqmin, oqmax, oozmin, oozmax,  &
     ivch, niu1)

!------------------------------------------------------------------------------

INTEGER, PARAMETER       ::                                       &
   ireals    = SELECTED_REAL_KIND (12,200),                       &
   iintegers = KIND  (1)

! Subroutine arguments
! scalar arguments with intent(in):
INTEGER(KIND=iintegers), INTENT(IN) ::       &
   nrttovid     ,  & ! number of RTTOV ids  requested
   platform(*)  ,  & ! number of platform. id's
   satellite(*) ,  & ! number of satellite. id's
   instrument(*),  & ! number of instrument. id's
   kpnsat       ,  & ! max no. of satellites
   kplev        ,  & ! no of rt levels
   kpch              ! max no. of channels

! Array  arguments with intent(inout):
INTEGER(KIND=iintegers), INTENT(INOUT) ::    &
   numchans(*)  ,  & ! Number of channels initialised
   niu1(*)           ! optional unit number to read
                                       ! rt_coef... files.
! Scalar arguments with intent(out):
INTEGER(KIND=iintegers), INTENT(OUT) ::      &
   kerr  ,         & ! error flag, returns kerr /= 0 if error
   kppf  ,         & ! max no. profiles processed in parallel
   kpchus,         & ! max no. of channels used
   kpnav ,         & ! max no of profile variables
   kpnsav,         & ! max no of surface variables
   kpnssv,         & ! max no of skin variables
   kpncv             ! max no of cloud variables

! Array  arguments with intent(out):
REAL (KIND=ireals), INTENT(OUT) ::           &
   preslev(kplev), & ! 43 pressure levels  (Pa)
   otmin  (kplev), & ! min temp array (K)
   otmax  (kplev), & ! max temp array (K)
   oqmin  (kplev), & ! min q array   (kg/kg)
   oqmax  (kplev), & ! max q array   (kg/kg)
   oozmin (kplev), & ! min ozone array   (kg/kg)
   oozmax (kplev)    ! max ozone array   (kg/kg)

INTEGER(KIND=iintegers), INTENT(INOUT) ::    &
   ivch(kpch,kpnsat)  ! array of valid channel numbers

!------------------------------------------------------------------------------

kerr = -99999

END SUBROUTINE RTTOV7_RTTVI

!==============================================================================
!==============================================================================

SUBROUTINE RTTOV7_RTTOVCLD &
     (knpf, klenpf, klevm,                                            &
      ppres, pangl, pangs, ksurf, ksat, knchpf,                       &
      kchan, kprof, pav, psav, pssv, pcvm, pap, paph,                 &
      pemis, ifail, prad, ptb, pradcld, ptbcld, tau, tausfc,          &
      pradovm, pcldemis, pait, pais)

!------------------------------------------------------------------------------

INTEGER, PARAMETER       ::                                       &
  ireals    = SELECTED_REAL_KIND (12,200),                       &
  iintegers = KIND  (1)

INTEGER(KIND=iintegers), INTENT(IN) ::     &
  knpf         , &  ! Number of profiles
  klenpf       , &  ! Length of input  profile vectors
  ksat         , &  ! Satellite index (see rttvi)
  knchpf       , &  ! Number of output radiances (= channels used * profiles)
  klevm        , &  ! Number of model (native) levels
  kchan(knchpf), &  ! Channel indices
  kprof(knchpf), &  ! Profiles indices
  ksurf(knpf)       ! Surface type index

REAL(KIND=ireals)      , INTENT(IN) ::     &
  ppres(klenpf)        , & ! Pressure levels (hpa) of atmospheric profile vectors
  pangl(knpf)          , & ! Satellite local zenith angle (deg)
  pangs(knpf)          , & ! Solar zenith angle at surface (deg)
  pav  (klenpf,4,knpf) , & ! Atmosp. profile variables
  psav (5,knpf)        , & ! Surface air variables
  pssv (6,knpf)        , & ! Surface skin variables
  pap  (knpf,klevm)    , & ! Full-level model pressures (hPa) of
                           !   atmospheric profile vectors
  paph (knpf,klevm+1)  , & ! Half-level model pressures (hPa) of
                           !   atmospheric profile vectors
  pcvm (knpf,klevm,4)      ! Temperature and cloud variables on klevm layers

REAL(KIND=ireals)      , INTENT(INOUT) ::  &
  pemis(knchpf)         !  surface emissivities


INTEGER(KIND=iintegers), INTENT(OUT)   ::  &
  ifail(knpf,9)         !  return flag

REAL(KIND=ireals)      , INTENT(OUT)   ::  &
  prad    (knchpf)          , & ! clear-sky radiances (mw/cm-1/ster/sq.m)
  ptb     (knchpf)          , & ! clear-sky brightness temperatures (K)
  pradcld (knchpf)          , & ! cloud-affected radiance
  ptbcld  (knchpf)          , & ! cloud-affected brightness temperature
  tau     (knchpf,klenpf)   , & ! clear-sky transmittance from each
                                ! standard pressure level
  tausfc  (knchpf)          , & ! clear-sky transmittance from surface
  pradovm (knchpf,2*klevm+2), & ! RT quantities for cloud computation 
  pcldemis(knchpf,klevm)    , & ! cloud emissivity
  pait    (knchpf,klevm+1)  , & ! toa weights of the cloud layers
  pais    (knchpf,klevm+1)      ! surface weights of the cloud layers

!------------------------------------------------------------------------------

ifail(:,:) = -99999

END SUBROUTINE RTTOV7_RTTOVCLD

!==============================================================================
!==============================================================================

SUBROUTINE RTTOV7_supsat(temp,wv,pres)

!------------------------------------------------------------------------------

INTEGER, PARAMETER       ::                                      &
  ireals    = SELECTED_REAL_KIND (12,200),                       &
  iintegers = KIND  (1)

REAL(KIND=ireals), INTENT(IN)  :: temp, pres
REAL(KIND=ireals), INTENT(OUT) :: wv

!------------------------------------------------------------------------------

wv = -999999.0_ireals

END SUBROUTINE RTTOV7_supsat

!==============================================================================
!==============================================================================

SUBROUTINE RTTOV7_spline(xi,yi,xo,yo,as,ni,no,ii)

!------------------------------------------------------------------------------

INTEGER, PARAMETER       ::                                      &
  ireals    = SELECTED_REAL_KIND (12,200),                       &
  iintegers = KIND  (1)

INTEGER(KIND=iintegers), INTENT(IN)  :: ni, no, ii

REAL(KIND=ireals),       INTENT(IN)  ::     &
   xi(ni), yi(ni), xo(no)

REAL(KIND=ireals),       INTENT(OUT) ::     &
   yo(no), as(5,ni)


yo(:)   = -999999.0_ireals
as(:,:) = -999999.0_ireals

END SUBROUTINE RTTOV7_spline

!==============================================================================
