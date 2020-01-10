!+ Module for latent heat nudging 
!-------------------------------------------------------------------------------

MODULE src_lheating

!-------------------------------------------------------------------------------
!
! Description:
! The module contains only one subroutine which collect all parts of the latent
! heating of the modell.
! The module is therefore called from following modules:
!  grid scale latent heat release:
!  - src_gscp: get the latent heat release due to cloud processes
!  - src_leapfrog | src_runge_kutta : latent heat release due to SATAD
!  - src_nudging : latent heat release due to SATAD
!  - src_relaxation : latent heat release due to SATAD
!  convective scale if the convection scheme is use so far
!  - src_conv_tiedke
!
! Current Code Owner: DWD, Christina Koepken
!    phone: +49  69  8062 2757
!    fax:   +49  69  8236 1493
!    email: christina.koepken@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 3.13       2004/12/03 Christina Koepken
!  Initial code
! 3.18       2006/03/03 Klaus Stephan
!  LHN namelist parameter moved to data_lheat_nudge to avoid to many dependencies
! V4_4         2008/07/16 Jens-Olaf Beismann
!  Introduced NEC compiler directives
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_23        2012/05/10 Ulrich Schaettler
!  Removed src_2timelevel and modified comments accordingly
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!===============================================================================

! Modules used:

USE data_parameters , ONLY:  wp, iintegers
USE data_runcontrol,  ONLY:  nnew
USE data_fields,      ONLY:  t
USE data_parallel,    ONLY:  my_cart_id
USE data_lheat_nudge, ONLY:  tt_lheat, ktop_lhn, kbot_lhn

USE environment,      ONLY:  model_abort

!===============================================================================

IMPLICIT NONE

!===============================================================================

CONTAINS

!===============================================================================
!+ Calculation and storage of gridscale latent heating rate (for lhn)
!-------------------------------------------------------------------------------

SUBROUTINE get_gs_lheating (cmode,kup,klow,tinc)

!-------------------------------------------------------------------------------
!
! Description:
! This subroutine collects all parts of model latent heating calculated at
! different locations where diabatic processes take place.
!
! Method:
! There are two different ways for calculating the dT from latent heating:
!
! Incremental updating by adding a temperature difference
! -------------------------------------------------------
! If the diabatic temperature contribution is directly added to the 
! temperature field  the routine is called before the process with cmode='add' 
! and after the process with cmode='inc' resulting in adding to the model 
! latent heating the temperature increment via  
! temperature difference after minus before the process.
!
!
! Direct updating of the latent heating with a temperature increment
! ------------------------------------------------------------------
! (added by Daniel Leuenberger)
!
! If the diabatic temperature contribution is indirectly (via ttens) added
! to the temperature field, the contribution of latent heating by the process
! cannot be added via temperature differences.
! Hence the model latent heating is updated by a direct temperature increment
! provided as optional argument (tinc) and setting cmode='dir'
!
!-------------------------------------------------------------------------------

! Subroutine / Function arguments
!-------------------------------------------------------------------------------
! Scalar arguments with intent(in):
  CHARACTER (LEN=3), INTENT(IN)        ::  cmode
  INTEGER (KIND=iintegers), INTENT(IN) ::       &
    kup,klow         ! upper and lower level index
  REAL (KIND=wp),     OPTIONAL, INTENT(IN) ::       &
       tinc(:,:,:)   ! direct temperature increment

! Local parameters, scalars, arrays :
!-------------------------------------------------------------------------------
! Local scalars:
!---------------

  CHARACTER (LEN=20)    ::  yroutine    ! name of the subroutine
  CHARACTER (LEN=80)    ::  yerrmsg     ! text message for model abort


!- End of header
!===============================================================================
!-------------------------------------------------------------------------------
! Begin of subroutine
!-------------------------------------------------------------------------------
  yroutine = 'get_gs_lheating'
!-------------------------------------------------------------------------------
! Section 1 : Storage of temperature field of the current timestep
!-------------------------------------------------------------------------------
!         1a: Store temperature field for later calculation of heating increment
!-------------------------------------------------------------------------------

  IF (cmode == "new" ) THEN
     IF (PRESENT(tinc)) THEN
!CDIR OUTERUNROLL=16
        tt_lheat(:,:,kup:klow,nnew) = tinc(:,:,kup:klow)
     ELSE
        yerrmsg =' ERROR  Optional argument tinc must be present for cmode = "new"'
        CALL model_abort (my_cart_id,6001,yerrmsg,yroutine)
     ENDIF

!-------------------------------------------------------------------------------
!         1b: Store t field - added to previously calculated increments
!-------------------------------------------------------------------------------
  ELSEIF(cmode == "add" ) THEN

!CDIR OUTERUNROLL=16
     tt_lheat(:,:,kup:klow,nnew) = tt_lheat(:,:,kup:klow,nnew)    &
                                        - t(:,:,kup:klow,nnew)

!-------------------------------------------------------------------------------
! Section 2 : Calculation of the temperature increment
!-------------------------------------------------------------------------------
  ELSEIF (cmode == "inc") THEN

!CDIR OUTERUNROLL=16
     tt_lheat(:,:,kup:klow,nnew) = tt_lheat(:,:,kup:klow,nnew)    &
                                        + t(:,:,kup:klow,nnew)
!-------------------------------------------------------------------------------
! Section 3 : Direct addition of temperature increment
!-------------------------------------------------------------------------------
  ELSEIF (cmode == "dir") THEN
     IF (PRESENT(tinc)) THEN
!CDIR OUTERUNROLL=16
        tt_lheat(:,:,kup:klow,nnew) = tt_lheat(:,:,kup:klow,nnew)    &
                                        + tinc(:,:,kup:klow)
     ELSE
        yerrmsg =' ERROR  Optional argument tinc must be present for cmode = "dir"'
        CALL model_abort (my_cart_id,6001,yerrmsg,yroutine)
     ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! End of subroutine
!-------------------------------------------------------------------------------

END SUBROUTINE get_gs_lheating

!===============================================================================

! End of module
!===============================================================================

END MODULE src_lheating
