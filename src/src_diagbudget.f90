!+ Source module for budget diagnostics computations in LM subdomains
!------------------------------------------------------------------------------
 
MODULE src_diagbudget
 
!------------------------------------------------------------------------------
!
! Description:
!   This module provides routines for budget diagnostic computations
!   in subdomains of the LM.
!
! Current Code Owner: DWD, Christina Koepken
!    phone: +49  69  8062 2757
!    fax:   +49  69  8236 1493
!    email: christina.koepken@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.33       1999/10/14 Reinhold Hess
!  Initial release
! 1.34       1999/12/10 Ulrich Schaettler
!  Changed calls to timing routines
! 1.39       2000/05/03 Ulrich Schaettler
!  Changed some variable names.
! 2.8        2001/07/06 Ulrich Schaettler
!  Eliminated non-necessary variables from the USE-lists
! 3.7        2004/02/18 Ulrich Schaettler
!  Replaced type_gscp by global variable lprog_qi
!  Renamed idiv_hum (tdiv_hum), cphi (crlat), acphir (acrlat)
! V4_10        2009/09/11 Christian Bollmann
!  Added compiler directive to use option _on_adb for NEC
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer
!  Replaced qx-variables by using them from the tracer module
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  Introduced MESSy interface
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
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
    wp,         & ! kind-type parameter for "normal" integer variables
    iintegers     ! KIND-type parameters for real variables
!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------

    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction

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

!   meridional direction
    jstart,       & ! start index for the forecast of w, t, qv, qc and pp
    jend,         & ! end index for the forecast of w, t, qv, qc and pp

! 4. constants for the horizontal rotated grid and related variables
! ------------------------------------------------------------------
    eddlon,       & ! 1 / dlon
    eddlat,       & ! 1 / dlat

! 5. variables for the time discretization and related variables
! --------------------------------------------------------------
    dt,           & ! long time-step

! 8. Organizational variables to handle the COSMO humidity tracers
! ----------------------------------------------------------------
    idt_qv, idt_qc, idt_qi

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------
    hhl,            & ! height of model half levels                   ( m )

! 2. external parameter fields                                        (unit)
! ----------------------------
    crlat      ,    & ! cosine of transformed latitude                  --
    acrlat     ,    & ! 1 / ( crlat * radius of the earth )           ( 1/m )

! 3. prognostic variables                                             (unit)
! -----------------------
    u,              & ! zonal wind speed                              ( m/s )
    v,              & ! meridional wind speed                         ( m/s )
    w,              & ! vertical wind speed (defined on half levels)  ( m/s )
    t,              & ! temperature                                   (  k  )
    pp                ! deviation from the reference pressure         ( pa  )

USE data_fields     , ONLY :   &

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------
    rho  ,          & ! density of moist air

! 7. fields for model output and diagnostics                          (unit )
! ------------------------------------------
    qvsflx     ,    & ! surface moisture flux                         (kg/m2*s)
    tdiv_hum   ,    & ! vertical sum for  divergence of humidity      (kg/m2)
    aevap_s           ! average evaporation flux (surface)            (kg/m2)

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    nold,         & ! corresponds to ntstep - 1
    nnow,         & ! corresponds to ntstep
    nnew,         & ! corresponds to ntstep + 1

! 3. controlling the physics
! --------------------------
    lphys,        & ! forecast with physical parametrizations
    lgsp            ! forecast with grid scale precipitation

!------------------------------------------------------------------------------

USE data_parallel, ONLY : my_cart_id
    
!------------------------------------------------------------------------------

USE environment  , ONLY : model_abort

!------------------------------------------------------------------------------

USE src_tracer   , ONLY : trcr_get, trcr_errorstr
USE data_tracer  , ONLY : T_ERR_NOTFOUND

!==============================================================================

IMPLICIT NONE

!==============================================================================
! Module procedures
!==============================================================================

CONTAINS

!==============================================================================
!+ Organization routine
!------------------------------------------------------------------------------

!option! -pvctl _on_adb
SUBROUTINE organize_diagbudget

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure is the driving routine for the diagnostic budget 
!   computations. It is called from the main program during the time stepping.
!   The vertically integrated divergence of the water content is computed in
!     subroutine 'horizontal divergence' and then added to 'tdiv_hum'.
!   Additionally evaporation and integrated water content and integrated
!     water vapout content are computed.
!
! Method:
!   Mainly call of subroutines
!
!------------------------------------------------------------------------------
!
! Local scalars:
  INTEGER (KIND=iintegers)   ::      izerror
  CHARACTER (LEN=255)        ::      yzerrmsg

! Local (automatic) arrays:
  REAL    (KIND=wp   )     ::  &
    psi  (ie,je,ke),   & !  variable to compute divergence
    div_h(ie,je   )      !  vertical integral horizontal divergence 

  REAL (KIND=wp),     POINTER :: &
    qv  (:,:,:)=> NULL(),       & ! QV at nnow
    qc  (:,:,:)=> NULL(),       & ! QC at nnow
    qi  (:,:,:)=> NULL()          ! QI at nnow

  CHARACTER(LEN=25)          :: yzroutine = 'organize_diagbudget'

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------
!- Begin Subroutine organize_diagbudget
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Horizontal Divergences of Selected Values
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1.1: Water Budget
!------------------------------------------------------------------------------

  ! retrieve the required microphysics tracers
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnow, ptr = qv)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nnow, ptr = qc)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qi, ptr_tlev = nnow, ptr = qi)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

! rho is computed for nnow in 'initialize_loop.incf', therefore this routine
! is called before the Asselin-Time filtering in the time-stepping

  IF ( ASSOCIATED(qi) ) THEN
    psi = rho * (qv(:,:,:) + qc(:,:,:) + qi(:,:,:))
  ELSE
    psi = rho * (qv(:,:,:) + qc(:,:,:))
  END IF

  CALL horizontal_divergence(psi,div_h)

  tdiv_hum(istart:iend,jstart:jend) =                 &
  tdiv_hum(istart:iend,jstart:jend) + dt * div_h(istart:iend,jstart:jend)

  ! Integrate the surface moisture flux in time
  aevap_s(istart:iend,jstart:jend) =                 &
    aevap_s(istart:iend,jstart:jend) + dt * qvsflx(istart:iend,jstart:jend)

!------------------------------------------------------------------------------
! End of module procedure organize_diagbudget
!------------------------------------------------------------------------------

END SUBROUTINE organize_diagbudget

!==============================================================================

!------------------------------------------------------------------------------

!option! -pvctl _on_adb
SUBROUTINE horizontal_divergence ( psi, div_h )

!------------------------------------------------------------------------------
!
! Description:
!  This subroutine computes the vertically integrated divergence of the
!   variable 'psi' and returns the result in 'div_h'.
!   Numerical conservativity of 'div_h' is provided.
!
! Method:
!  Multiply 'psi' with layer thickness, multiplication with 'u' and 'v'.
!  Then sum up all layers and compute divergence of the sum.
!
!------------------------------------------------------------------------------

! Parameterlist:
! -------------

REAL (KIND=wp),     INTENT (IN)          ::    &
  psi  (ie,je,ke)      !  variable to compute divergence

REAL (KIND=wp),     INTENT (OUT)         ::    &
  div_h(ie,je   )      !  vertical integral horizontal divergence of psi

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    i,  j,  k              !  Loop indices in longitudinal, latitudinal and
                           !  vertical direction and error status variable

! Local (automatic) arrays:
! ------------
  REAL    (KIND=wp   )     ::  &
    psih  (ie,je,ke),    & !
    psi_u (ie,je,ke),    & !
    psi_v (ie,je,ke),    & !
    psi_uh(ie,je   ),    & !
    psi_vh(ie,je   ),    & !
    zfadsx(   je   )       !

!
! End of header 
!==============================================================================
 
!------------------------------------------------------------------------------
! Begin Subroutine horizontal_divergence
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1:  Setup of metrical factors and interpolation coefficients
!------------------------------------------------------------------------------
   
  ! Metrical factors for horizontal discretization  

  DO  j = jstart-1 , jend+1
    zfadsx(j) = acrlat(j,1)*eddlon
  ENDDO          

!------------------------------------------------------------------------------
! Section 2:  Calculation of horizontal divergence on full levels
!------------------------------------------------------------------------------

  DO k = 1, ke
    psih(:,:,k) = psi(:,:,k) * (hhl(:,:,k) - hhl(:,:,k+1))
  ENDDO 

  DO      k = 1,         ke
    DO    j = jstart-1 , jend
      DO  i = istart-1 , iend
        psi_u(i,j,k) = 0.5_wp * (psih(i,j,k)+psih(i+1,j,k)) * u(i,j,k,nnow)
        psi_v(i,j,k) = 0.5_wp * (psih(i,j,k)+psih(i,j+1,k)) * v(i,j,k,nnow)    &
                                                            * crlat(j,2)
      ENDDO
    ENDDO
  ENDDO

  psi_uh(istart-1:iend,jstart-1:jend) =        &
            SUM(psi_u(istart-1:iend,jstart-1:jend,:),3)
  psi_vh(istart-1:iend,jstart-1:jend) =        &
            SUM(psi_v(istart-1:iend,jstart-1:jend,:),3)

  DO    j = jstart , jend
    DO  i = istart , iend
      div_h(i,j) = zfadsx(j) * (psi_uh(i,j)-psi_uh(i-1,j)) +           &
                   eddlat * (acrlat(j,2)*psi_vh(i,j)-acrlat(j-1,2)     &
                                                        * psi_vh(i,j-1))
    END DO
  END DO
  
!------------------------------------------------------------------------------
! End of module procedure "horizontal divergence"
!------------------------------------------------------------------------------
END SUBROUTINE horizontal_divergence

!==============================================================================

END MODULE src_diagbudget
