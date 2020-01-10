!+ Interface module for organizing the microphysics in blocked data structure
!------------------------------------------------------------------------------

MODULE gscp_interface

!------------------------------------------------------------------------------
!
! Description:
!
! This module contains calls to initialize the "copy to block" and organizes
! the calls to the block version of the microphysics schemes.
! 
!
! Current Code Owner: MeteoSwiss, Xavier Lapillonne
!  phone: +41 58 460 9237
!  fax  : +41 58 460 9278
!  email:  xavier.lapillonne@meteoswiss.ch
!
! History:
! Version      Date       Name
! ----------   ---------- ----
! V5_1         2014-11-28 Xavier Lapillonne
!  Initial release
! V5_2         2015-05-21 Ulrich Schaettler
!  Introduced ifdef NUDGING directive to be able to compile without Nudging
! V5_3         2015-10-09 Ulrich Blahak
!  Also put tt_lheat_new_b and qrsflux_b to ifdef NUDGING
!  Introduced possibility of calling microphysics at the beginning of all
!   parameterizations
! V5_4         2016-03-10 Xavier Lapillonne
!  Deactivate OpenACC statements for the moment being
!  Updated code owner information
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Declarations:
!
! Modules used:
!
USE kind_parameters,  ONLY : wp

USE data_modelconfig, ONLY : dt, dt2, ke, idt_qg, idt_qr, idt_qs, idt_qi

USE data_runcontrol,  ONLY : itype_gscp, idbg_level, ldebug_gsp, lprintdeb_all, &
                             l2tls, nproma, nlastproma, nblock, ntstep, ltime,  &
                             lgsp_first

USE data_constants,   ONLY : r_d, rvd_m_o

USE data_parallel,    ONLY : my_cart_id

!------------------------------------------------------------------------------

USE data_fields,      ONLY : tinc_lh, prg_gsp
USE src_stoch_physics,ONLY : pertstoph

#ifdef NUDGING
USE data_lheat_nudge, ONLY : tt_lheat,qrsflux
USE data_block_fields,ONLY : tt_lheat_nx_b, qrsflux_b
#endif


USE data_block_fields,ONLY : t_nx_b, ptot_b, pp_nx_b, p0_b, rho0_b, dp0_b,  &
                             hhl_b, dz_b, qv_nx_b, qc_nx_b, qi_nx_b,        &
                             qr_nx_b, qs_nx_b, qg_nx_b, qrs_b, prr_gsp_b,   &
                             prs_gsp_b, prg_gsp_b, tinc_lh_b,               &
                             rho_b, pertstoph_b, ps_b

!------------------------------------------------------------------------------

USE gscp_data,        ONLY : gscp_set_coefficients
USE gscp_kessler,     ONLY : kessler
USE gscp_hydor,       ONLY : hydor
USE gscp_cloudice,    ONLY : cloudice
USE gscp_graupel,     ONLY : graupel

USE src_block_fields, ONLY :  register_copy, CopylistStruct,               &
                              copyToBlockF, copyFromBlockF,                &
                              mind_ilon, mind_jlat, init_copy_list

USE meteo_utilities,  ONLY :  calrho_block, calps_block

!                             mind_ilon, mind_jlat, request_copy,          &
!                             copy_to_block, copy_from_block,              &
!                             finalize_copy, init_copy_list
!USE time_utilities,   ONLY : get_timings, i_copyblocks, i_precipitation

!==============================================================================

IMPLICIT NONE

!------------------------------------------------------------------------------

PRIVATE

!==============================================================================

!------------------------------------------------------------------------------
! Global (i.e. public) Declarations:
PUBLIC :: gscp_init, gscp_organize, gscp_init_copy, gscpCopyList

!==============================================================================
! Module variables
TYPE(CopylistStruct) :: gscpCopyList

!==============================================================================
! Module procedures in "gscp_interface" 
!==============================================================================

CONTAINS 

!==============================================================================
!+ Module procedure "gscp_init" in "gscp_interface" 
!------------------------------------------------------------------------------    

SUBROUTINE gscp_init

!------------------------------------------------------------------------------
!
! Description:
!   initialize coefficients for microphysics
!
!------------------------------------------------------------------------------

! Locals
! ------

INTEGER :: izdebug      ! local debug level

!----------- End of header ----------------------------------------------------

  ! init debugging message level
  IF (lprintdeb_all) THEN
    izdebug = idbg_level
  ELSE
    IF (my_cart_id == 0) THEN
      izdebug = idbg_level
    ELSE
      izdebug = 0
    ENDIF
  ENDIF

  ! initializations for microphysics
  SELECT CASE(itype_gscp)
  CASE (1, 2, 3, 4)
     CALL gscp_set_coefficients (izdebug)
  END SELECT

END SUBROUTINE gscp_init

!==============================================================================
!==============================================================================
!+ Module procedure "gscp_init_copy" in "gscp_interface" 
!------------------------------------------------------------------------------    

SUBROUTINE gscp_init_copy

!------------------------------------------------------------------------------
!
! Description:
!   Register all required copies to/from block for gscp schemes 
!
!------------------------------------------------------------------------------

  ! init copy list
  CALL init_copy_list(gscpCopyList)

  ! Register required copy to block
  ! Variables with in intent IN
  CALL register_copy(hhl_b    ,gscpCopyList,copyToBlockF)
  CALL register_copy(t_nx_b   ,gscpCopyList,copyToBlockF)
  CALL register_copy(p0_b     ,gscpCopyList,copyToBlockF)
  CALL register_copy(pp_nx_b  ,gscpCopyList,copyToBlockF)
  CALL register_copy(rho_b    ,gscpCopyList,copyToBlockF)
  CALL register_copy(qv_nx_b  ,gscpCopyList,copyToBlockF)
  CALL register_copy(qc_nx_b  ,gscpCopyList,copyToBlockF)
  IF (idt_qi>=0)            CALL register_copy(qi_nx_b    ,gscpCopyList,copyToBlockF)
  IF (idt_qr>=0)            CALL register_copy(qr_nx_b    ,gscpCopyList,copyToBlockF)
  IF (idt_qs>=0)            CALL register_copy(qs_nx_b    ,gscpCopyList,copyToBlockF)
  IF (idt_qg>=0)            CALL register_copy(qg_nx_b    ,gscpCopyList,copyToBlockF)
  IF (ALLOCATED(tinc_lh))   CALL register_copy(tinc_lh_b  ,gscpCopyList,copyToBlockF)
  IF (ALLOCATED(pertstoph)) CALL register_copy(pertstoph_b,gscpCopyList,copyToBlockF)
#ifdef NUDGING
  IF (ALLOCATED(tt_lheat))  CALL register_copy(tt_lheat_nx_b,gscpCopyList,copyToBlockF)
  IF (ALLOCATED(qrsflux))   CALL register_copy(qrsflux_b    ,gscpCopyList,copyToBlockF)
#endif
  !fields required to update rho and ps after microphysics (when called at the beginning
  !of the time loop)
  CALL register_copy(ps_b     ,gscpCopyList,copyToBlockF)
  CALL register_copy(rho_b    ,gscpCopyList,copyToBlockF)
  CALL register_copy(rho0_b   ,gscpCopyList,copyToBlockF)
  CALL register_copy(dp0_b    ,gscpCopyList,copyToBlockF)


  !Variables with intent OUT
  CALL register_copy(t_nx_b  ,gscpCopyList,copyFromBlockF)
  CALL register_copy(qv_nx_b ,gscpCopyList,copyFromBlockF)
  CALL register_copy(qc_nx_b ,gscpCopyList,copyFromBlockF)
  IF (idt_qi>=0)           CALL register_copy(qi_nx_b  ,gscpCopyList,copyFromBlockF)
  IF (idt_qr>=0)           CALL register_copy(qr_nx_b  ,gscpCopyList,copyFromBlockF)
  IF (idt_qs>=0)           CALL register_copy(qs_nx_b  ,gscpCopyList,copyFromBlockF)
  IF (idt_qg>=0)           CALL register_copy(qg_nx_b  ,gscpCopyList,copyFromBlockF)
  IF (ALLOCATED(tinc_lh))  CALL register_copy(tinc_lh_b,gscpCopyList,copyFromBlockF)
#ifdef NUDGING
  IF (ALLOCATED(tt_lheat)) CALL register_copy(tt_lheat_nx_b, gscpCopyList,copyFromBlockF)
  IF (ALLOCATED(qrsflux))  CALL register_copy(qrsflux_b    , gscpCopyList,copyFromBlockF)
#endif  
  CALL register_copy(prr_gsp_b ,gscpCopyList,copyFromBlockF)
  IF (idt_qs>=0)           CALL register_copy(prs_gsp_b ,gscpCopyList,copyFromBlockF)
  IF (ALLOCATED(prg_gsp))  CALL register_copy(prg_gsp_b ,gscpCopyList,copyFromBlockF)
  CALL register_copy(qrs_b     ,gscpCopyList,copyFromBlockF)
  CALL register_copy(ps_b      ,gscpCopyList,copyFromBlockF)
  CALL register_copy(rho_b     ,gscpCopyList,copyFromBlockF)

END SUBROUTINE gscp_init_copy

!==============================================================================
!==============================================================================
!+ Module procedure to select the scheme for grid-scale clouds and precipitation
!------------------------------------------------------------------------------    

SUBROUTINE gscp_organize (ib, ipend, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   This module selects the parameterization scheme to calculate the effects
!   of grid-scale clouds and precipitation. Depending on the namelist input
!   parameter "itype_gscp", the following routines may be called, which calculate
!   the rates of change of temperature, cloud water, water vapor and:
!
!     (1) "kessler" : rain
!     (2) "hydor"   : rain, snow
!     (3) "cloudice": rain, snow, cloud ice
!     (4) "graupel" : rain, snow, cloud ice and graupel
!
!------------------------------------------------------------------------------

! Subroutines arguments
! --------------------

  INTEGER, INTENT(IN)  :: ib                 ! current block index
  INTEGER, INTENT(IN)  :: ipend              ! length of current block
  INTEGER, INTENT(OUT) :: ierror
  CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg

! Locals
! ------

  INTEGER :: izdebug      !local debug level
  INTEGER :: ip, k
  REAL(KIND=wp) :: zdt

!----------- End of header ----------------------------------------------------

!------------------------------------------------------------------------------
! Begin subroutine gscp_organize
!------------------------------------------------------------------------------
  
  !NOacc data                                          &
  !NOacc present(ptot_b,pp_nx_b,p0_b,dz_b,hhl_b,qrs_b) &
  !NOacc present(qr_nx_b,qs_nx_b,qi_nx_b,qg_nx_b)

  ierror     = 0
  yerrmsg(:) = ' '

  ! init debugging message level
  IF (lprintdeb_all) THEN
    izdebug = idbg_level
  ELSE
    IF (my_cart_id == 0) THEN
      izdebug = idbg_level
    ELSE
      izdebug = 0
    ENDIF
  ENDIF

  ! adapt the time step
  IF ( l2tls ) THEN
     zdt   = dt
  ELSE
     zdt   = dt2
  ENDIF

  !compute auxiliary fields used in microphysics
  !NOacc parallel
  DO k = 1, ke
  !NOacc loop gang vector
    DO ip=1,ipend
      ptot_b(ip,k) = pp_nx_b (ip,k) + p0_b(ip,k)
      dz_b(ip,k)   = hhl_b(ip,k) - hhl_b(ip,k+1)
    END DO
  END DO
  !NOacc end parallel

  SELECT CASE (itype_gscp)

  CASE(1)  !  (warm rain scheme: kessler)

    CALL kessler  (                      &
      nvec   = nproma               ,    & !> in:  actual array size
      ke     = ke                   ,    & !< in:  actual array size
      ivstart=  1                   ,    & !< in:  start index of calculation
      ivend  = ipend                ,    & !< in:  end index of calculation
      kstart =  1                   ,    & !< in:  vertical start index
      zdt    = zdt                  ,    & !< in:  timestep
      dz     = dz_b                 ,    & !< in:  vertical layer thickness
      t      = t_nx_b               ,    & !< in:  temp,tracer,...
      p      = ptot_b               ,    & !< in:  full level pres
      rho    = rho_b                ,    & !< in:  density
      qv     = qv_nx_b              ,    & !< in:  spec. humidity
      qc     = qc_nx_b              ,    & !< in:  cloud water
      qr     = qr_nx_b              ,    & !< in:  rain water
      prr_gsp= prr_gsp_b            ,    & !< out: precipitation rate of rain
      tinc_lh= tinc_lh_b            ,    & !  t-increment due to latent heat 
!     pstoph = pertstoph_b          ,    & !> in:  stochastic multiplier of physics tendencies
#ifdef NUDGING
      tt_lheat= tt_lheat_nx_b       ,    & !  t-increments due to latent heating (nud) 
!     qrsflux = qrsflux_b           ,    & !  total precipitation flux
#endif
      idbg=izdebug                  ,    &
      l_cv=.TRUE. )

    ! update qrs
    !NOacc parallel
    DO k=1,ke
    !NOacc loop gang vector
      DO ip=1,ipend
        qrs_b(ip,k) = qr_nx_b(ip,k)
      END DO
    END DO
    !NOacc end parallel

  CASE(2)  !  (warm rain scheme with snow: hydor)

    CALL hydor    (                      &
      nvec   = nproma               ,    & !> in:  actual array size
      ke     = ke                   ,    & !< in:  actual array size
      ivstart=  1                   ,    & !< in:  start index of calculation
      ivend  = ipend                ,    & !< in:  end index of calculation
      kstart =  1                   ,    & !< in:  vertical start index
      zdt    = zdt                  ,    & !< in:  timestep
      dz     = dz_b                 ,    & !< in:  vertical layer thickness
      t      = t_nx_b               ,    & !< in:  temp,tracer,...
      p      = ptot_b               ,    & !< in:  full level pres
      rho    = rho_b                ,    & !< in:  density
      qv     = qv_nx_b              ,    & !< in:  spec. humidity
      qc     = qc_nx_b              ,    & !< in:  cloud water
      qr     = qr_nx_b              ,    & !< in:  rain water
      qs     = qs_nx_b              ,    & !< in:  snow
      prr_gsp= prr_gsp_b            ,    & !< out: precipitation rate of rain
      prs_gsp= prs_gsp_b            ,    & !< out: precipitation rate of snow
      tinc_lh= tinc_lh_b            ,    & !  t-increment due to latent heat 
!     pstoph = pertstoph_b          ,    & !> in:  stochastic multiplier of physics tendencies
#ifdef NUDGING
      tt_lheat= tt_lheat_nx_b       ,    & !  t-increments due to latent heating (nud) 
!     qrsflux = qrsflux_b           ,    & !  total precipitation flux
#endif
      idbg=izdebug                  ,    &
      l_cv=.TRUE. )

    ! update qrs
    !NOacc parallel
    DO k=1,ke
    !NOacc loop gang vector
      DO ip=1,ipend
        qrs_b(ip,k) = qr_nx_b(ip,k) + qs_nx_b(ip,k)
      END DO
    END DO
    !NOacc end parallel

  CASE(3)  !  (2-cat ice: cloud ice, snow)

    CALL cloudice (                      &
      nvec   = nproma               ,    & !> in:  actual array size
      ke     = ke                   ,    & !< in:  actual array size
      ivstart=  1                   ,    & !< in:  start index of calculation
      ivend  = ipend                ,    & !< in:  end index of calculation
      kstart =  1                   ,    & !< in:  vertical start index
      zdt    = zdt                  ,    & !< in:  timestep
      dz     = dz_b                 ,    & !< in:  vertical layer thickness
      t      = t_nx_b               ,    & !< in:  temp,tracer,...
      p      = ptot_b               ,    & !< in:  full level pres
      rho    = rho_b                ,    & !< in:  density
      qv     = qv_nx_b              ,    & !< in:  spec. humidity
      qc     = qc_nx_b              ,    & !< in:  cloud water
      qi     = qi_nx_b              ,    & !< in:  cloud ice
      qr     = qr_nx_b              ,    & !< in:  rain water
      qs     = qs_nx_b              ,    & !< in:  snow
      prr_gsp= prr_gsp_b            ,    & !< out: precipitation rate of rain
      prs_gsp= prs_gsp_b            ,    & !< out: precipitation rate of snow
      tinc_lh= tinc_lh_b            ,    & !  t-increment due to latent heat 
      pstoph = pertstoph_b          ,    & !> in:  stochastic multiplier of physics tendencies
#ifdef NUDGING
      tt_lheat= tt_lheat_nx_b       ,    & !  t-increments due to latent heating (nud) 
      qrsflux = qrsflux_b           ,    & !  total precipitation flux
#endif
      idbg=izdebug                  ,    &
      l_cv=.TRUE. )

    ! update qrs
    !NOacc parallel
    DO k=1,ke
    !NOacc loop gang vector
      DO ip=1,ipend
        qrs_b(ip,k) = qr_nx_b(ip,k)+qs_nx_b(ip,k)+qi_nx_b(ip,k)
      END DO
    END DO
    !NOacc end parallel

  CASE(4)  

    CALL graupel  (                      &
      nvec   = nproma               ,    & !> in:  actual array size
      ke     = ke                   ,    & !< in:  actual array size
      ivstart=  1                   ,    & !< in:  start index of calculation
      ivend  = ipend                ,    & !< in:  end index of calculation
      kstart =  1                   ,    & !< in:  vertical start index
      zdt    = zdt                  ,    & !< in:  timestep
      dz     = dz_b                 ,    & !< in:  vertical layer thickness
      t      = t_nx_b               ,    & !< in:  temp,tracer,...
      p      = ptot_b               ,    & !< in:  full level pres
      rho    = rho_b                ,    & !< in:  density
      qv     = qv_nx_b              ,    & !< in:  spec. humidity
      qc     = qc_nx_b              ,    & !< in:  cloud water
      qi     = qi_nx_b              ,    & !< in:  cloud ice
      qr     = qr_nx_b              ,    & !< in:  rain water
      qs     = qs_nx_b              ,    & !< in:  snow
      qg     = qg_nx_b              ,    & !< in:  graupel
      prr_gsp= prr_gsp_b            ,    & !< out: precipitation rate of rain
      prs_gsp= prs_gsp_b            ,    & !< out: precipitation rate of snow
      prg_gsp= prg_gsp_b            ,    & !< out: precipitation rate of snow
      tinc_lh= tinc_lh_b            ,    & !  t-increment due to latent heat 
      pstoph = pertstoph_b          ,    & !> in:  stochastic multiplier of physics tendencies
#ifdef NUDGING
      tt_lheat= tt_lheat_nx_b      ,    & !  t-increments due to latent heating (nud) 
      qrsflux = qrsflux_b          ,    & !  total precipitation flux
#endif
      idbg=izdebug                 ,    &
      l_cv=.TRUE. )  
          
    ! update qrs
    !NOacc parallel
    DO k=1,ke
    !NOacc loop gang vector
      DO ip=1,ipend
        qrs_b(ip,k) = qr_nx_b(ip,k)+qs_nx_b(ip,k)+qi_nx_b(ip,k)+qg_nx_b(ip,k)
      END DO
    END DO
    !NOacc end parallel

  CASE DEFAULT
    ! This shoud not happen
    PRINT*, 'No valid itype_gscp:  ', itype_gscp

  END SELECT

  IF (lgsp_first) THEN
    !Need to update rho and ps when calling GSP before the other physics
    CALL calrho_block ( t_nx_b, pp_nx_b, qv_nx_b, qc_nx_b,  &
                  qrs_b, p0_b, rho_b, nproma, ke, r_d, rvd_m_o, lacc=.TRUE.)

    !ps at nnow
    CALL calps_block ( ps_b, pp_nx_b(:,ke), t_nx_b(:,ke),          &
                 qv_nx_b(:,ke ), qc_nx_b(:,ke ), qrs_b(:,ke),      &
                 rho0_b(:,ke), p0_b(:,ke), dp0_b(:,ke),            &
                 nproma, rvd_m_o, r_d,                             &
                 1, nproma, lacc=.TRUE.)
  ENDIF

  !NOacc end data

!------------------------------------------------------------------------------
! End of module procedure gscp_organize
!------------------------------------------------------------------------------

END SUBROUTINE gscp_organize
       
!==============================================================================
!==============================================================================

END MODULE gscp_interface
