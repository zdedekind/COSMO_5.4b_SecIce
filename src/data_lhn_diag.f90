!+ Data module for diagnostic fields for Latent Heat Nudging
!------------------------------------------------------------------------------

MODULE data_lhn_diag

!------------------------------------------------------------------------------
!
! Description:
!  This module declares all diagnostic fields for the diagnosis of
!   the LHN approach that have to reside in
!
!  All fields are declared as allocatable arrays. They are allocated in the
!  setup of the model and deallocated in the cleanup at the end of the
!  program.
!
! Current Code Owner: DWD, Klaus Stephan
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 3.21       2006/12/04 Klaus Stephan
!  Initial Release
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  MESSy interface introduced
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Eliminated MESSy interface, because MESSy cannot be used with data assimilation
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
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
USE data_parameters , ONLY :   &
   wp,       & ! KIND-type parameters for real variables
   iintegers   ! kind-type parameter for "normal" integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Global (i.e. public) Declarations:

  REAL  (KIND=wp), TARGET, ALLOCATABLE ::           &
    tt_lheat_o(:,:,:),& ! GRIB output variable for tt_lheat
    tinc_lhn_o(:,:,:),& ! GRIB output variable for tinc_lhn
    ttm_cv_o(:,:,:)     ! GRIB output variable for ttm_cv

  REAL  (KIND=wp), TARGET, ALLOCATABLE ::           &
    t_lh_mod1(:,:,:),& ! model latent heating rate due to
                       ! rain processes (src_gscp.f90)                ( K/s )
    qv_lh_mod1(:,:,:),&! rate of specific humidity change due to lh   (kg/(kg*s))
    qc_lh_mod1(:,:,:),&! rate of cloud water change due to lh   (kg/(kg*s))

    t_lh_mod2(:,:,:),& ! model latent heating rate due to
                       ! condesation/evaporation of cloud water (satad)  ( K/s )
    qv_lh_mod2(:,:,:),&! rate of specific humidity change due to lh   (kg/(kg*s))
    qc_lh_mod2(:,:,:),&! rate of cloud water change due to lh   (kg/(kg*s))

    t_lh_mod3(:,:,:),& ! model latent heating rate due to
                       ! corrections in src_nudging.f90 (satad)       ( K/s )
    qv_lh_mod3(:,:,:),&! rate of specific humidity change due to lh   (kg/(kg*s))
    qc_lh_mod3(:,:,:),&! rate of cloud water change due to lh   (kg/(kg*s))

    t_lh_mod4(:,:,:),& ! model latent heating rate due to
                       ! corrections in src_relaxation.f90 (satad)       ( K/s )
    qv_lh_mod4(:,:,:),&! rate of specific humidity change due to lh   (kg/(kg*s))
    qc_lh_mod4(:,:,:)  ! rate of cloud water change due to lh   (kg/(kg*s))

  REAL  (KIND=wp), TARGET, ALLOCATABLE ::           &
    dt_gscp(:,:,:)    ,& ! diagnostic variable for the rate of change in temperature
    dt_turb(:,:,:)    ,& ! diagnostic variable for the rate of change in temperature
    dt_rad(:,:,:)     ,& ! diagnostic variable for the rate of change in temperature
    dt_horad(:,:,:)   ,& ! diagnostic variable for the rate of change in temperature
    dt_slow(:,:,:)    ,& ! diagnostic variable for the rate of change in temperature
    dt_fast(:,:,:)    ,& ! diagnostic variable for the rate of change in temperature
    dt_hordif(:,:,:)  ,& ! diagnostic variable for the rate of change in temperature
    dt_verdif(:,:,:)  ,& ! diagnostic variable for the rate of change in temperature
    dt_assi(:,:,:)    ,& ! diagnostic variable for the rate of change in temperature
    dt_rel(:,:,:)     ,& ! diagnostic variable for the rate of change in temperature
    dt_rel2(:,:,:)    ,& ! diagnostic variable for the rate of change in temperature
    dqv_gscp(:,:,:)   ,& ! diagnostic variable for the rate of change in qv
    dqv_turb(:,:,:)   ,& ! diagnostic variable for the rate of change in qv
    dqv_horad(:,:,:)  ,& ! diagnostic variable for the rate of change in qv
    dqv_slow(:,:,:)   ,& ! diagnostic variable for the rate of change in qv
    dqv_hordif(:,:,:) ,& ! diagnostic variable for the rate of change in qv
    dqv_verdif(:,:,:) ,& ! diagnostic variable for the rate of change in qv
    dqv_assi(:,:,:)   ,& ! diagnostic variable for the rate of change in qv
    dqv_lhn(:,:,:)    ,& ! diagnostic variable for the rate of change in qv
    dqv_rel(:,:,:)    ,& ! diagnostic variable for the rate of change in qv
    dqv_rel2(:,:,:)   ,& ! diagnostic variable for the rate of change in qv
    dqc_gscp(:,:,:)   ,& ! diagnostic variable for the rate of change in qc
    dqc_turb(:,:,:)   ,& ! diagnostic variable for the rate of change in qc
    dqc_horad(:,:,:)  ,& ! diagnostic variable for the rate of change in qc
    dqc_slow(:,:,:)   ,& ! diagnostic variable for the rate of change in qc
    dqc_hordif(:,:,:) ,& ! diagnostic variable for the rate of change in qc
    dqc_verdif(:,:,:) ,& ! diagnostic variable for the rate of change in qc
    dqc_assi(:,:,:)   ,& ! diagnostic variable for the rate of change in qc
    dqc_lhn(:,:,:)    ,& ! diagnostic variable for the rate of change in qc
    dqc_rel(:,:,:)    ,& ! diagnostic variable for the rate of change in qc
    dqc_rel2(:,:,:)   ,& ! diagnostic variable for the rate of change in qc
    dw_horad(:,:,:)   ,&
    dw_hordif(:,:,:)   ,&
    dw_bouy(:,:,:)    ,&
    dw_slow(:,:,:)    ,&
    dw_fast(:,:,:)    ,&
    dw_rel(:,:,:)

  LOGICAL :: lhn_diag_tens

!==============================================================================

END MODULE data_lhn_diag
