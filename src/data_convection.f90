!+ Data module for variables of the convection parameterizations
!------------------------------------------------------------------------------

MODULE data_convection

!------------------------------------------------------------------------------
!
! Description:
!  This module contains variables that are used in the various convection
!  parameterizations.
!
! Current Code Owner: DWD, Dmitrii Mironov
!  phone:  +49  69  8062 2705
!  fax:    +49  69  8062 3721
!  email:  Dmitrii.Mironov@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V4_5         2008/09/10 Ulrich Schaettler
!  Initial release
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_18        2011/05/26 Ulrich Schaettler
!  New variable thick_sc (by Martin Koehler)
!  Changed the code owner
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!
! Declarations:
!
! Modules used:

USE data_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================


! Variables for shallow convection (also used in Tiedtke)
! ----------------------

  REAL (KIND=wp)     ::           &
    ! mean entrainment rate for shallow convection
    entr_sc  = 0.00030_wp

  REAL (KIND=wp)     ::           &
    ! limit for convective clouds to be "shallow" (in Pa)
    ! Shallow convection parameterization becomes active only if cloud
    ! thickness from cloud base to cloud top exceeds a threshold.  To evaluate
    ! this condition a parcel is launched.  This threshold is typically set to
    ! values between 200hPa and 300hPa with a COSMO DE default of 250hPa.
    thick_sc = 25000.0_wp

    ! COSMO-DE default (by Guenther Doms) : 250 hPa
    ! IFS      default (by Peter Bechtold): 200 hPa
    ! reasonable values:  between 100 hPa and 450 hPa

!==============================================================================

END MODULE data_convection
