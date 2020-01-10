!
!+ 3DVAR/COSMO module to define kind parameters.
!
! $Id: mo_kind.f90,v 5.0 2013-11-08 10:08:04 for0adm Exp $
!
!==============================================================================
!> Define kind parameters
!>
MODULE mo_kind
!
!-------------------------------------------------------------------------------
! Description:
!   Define kind parameters.
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V4_22        2012/01/31 Christoph Schraff
!  Initial release, based on 3DVAR version V1_10.
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
! Authors:
! Andreas Rhodin  MPIfM/DWD  1999  initial revision
!                            2007  changed for NEC SX
!------------------------------------------------------------------------------
   implicit none
   public

   !---------------------
   ! Real kind parameters
   !---------------------
   !
   !>  8 byte real kind parameter
   !
   integer, parameter :: dp = selected_real_kind(13)
   !
   !>  4 byte real kind parameter
   !
   integer, parameter :: sp = selected_real_kind(6)
   !
   !>  working precision used in the program
   !
   integer, parameter :: wp = dp

   !------------------------
   ! Integer kind parameters
   !------------------------
   !
   !>  8 byte integer kind parameter
   !
   integer, parameter :: i8 = selected_int_kind(14)  ! 8 byte integer
   !
   !>  4 byte integer kind parameter
   !
   integer, parameter :: i4 = selected_int_kind (9)  ! 4 byte integer
   !
   !>  2 byte integer kind parameter
   !
   integer, parameter :: i2 = selected_int_kind (4)  ! 2 byte integer
   !
   !>  1 byte integer kind parameter
   !
   integer, parameter :: i1 = selected_int_kind (2)  ! 1 byte integer

   !---------------------------------------------------------
   !> kind parameter for integer sufficient to hold a pointer
   !---------------------------------------------------------
#ifdef CP4
   integer, parameter :: cp = i4
#else
   integer, parameter :: cp = i8
#endif

end module mo_kind

!==============================================================================
