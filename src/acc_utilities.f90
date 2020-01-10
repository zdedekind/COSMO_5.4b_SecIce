!+ Module acc_utilities providing utility functions for OpenACC execution
!------------------------------------------------------------------------------

MODULE acc_utilities

!------------------------------------------------------------------------------
!
! Description:
!
! This module provides a set of utility routines which can be used for
! easier OpenACC execution.
! Some tools to synchronize host and cpu memory are provided.
!
! Method:
!
! OpenACC directives are used to provide functionality
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
! V5_4         2016-03-10 Xavier Lapillonne, Oliver Fuhrer
!  Initial release
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

USE kind_parameters,    ONLY: wp

USE data_runcontrol,    ONLY: idbg_level

USE data_io,            ONLY: num_bd_fields, bd_list

USE data_parallel,      ONLY: my_cart_id

!------------------------------------------------------------------------------

IMPLICIT NONE

!------------------------------------------------------------------------------

PRIVATE

!==============================================================================

!------------------------------------------------------------------------------
! Global (i.e. public) Declarations:

PUBLIC :: acc_update_device, acc_update_bd_fields, acc_copyin_1d

! Interfaces for overloaded subroutines

INTERFACE acc_update_device
   MODULE PROCEDURE acc_update_device_3d, acc_update_device_2d
END INTERFACE

!==============================================================================
! Module Procedures in src_tracer
!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure "acc_update_bd_fields" in "acc_utilites" to copy BD fields
!------------------------------------------------------------------------------

SUBROUTINE acc_update_bd_fields( ierr, yerrmsg )

!------------------------------------------------------------------------------
!
! Description:
!
! This subroutine copies the boundary fields read from disc to GPU
! using the data gathered in the src_input subroutine.
!
! Method:
!
! Use OpenACC !$update device directives
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

 INTEGER           , INTENT(OUT) :: ierr      ! error status
 CHARACTER (LEN=80), INTENT(OUT) :: yerrmsg   ! error message

! Local parameters:
! ----------------

 INTEGER  :: i

!------------------------------------------------------------------------------
! Begin Subroutine acc_update_bd_fields
!------------------------------------------------------------------------------

  ierr = 0
  yerrmsg = ''

#ifdef _OPENACC
  DO i = 1, num_bd_fields
    IF ( ( my_cart_id == 0 ) .AND. ( idbg_level > 8 ) ) THEN
      WRITE(*,*) ' GPUINFO: DEVICE UPDATE BD-FIELD: ' // &
        TRIM(bd_list(i)%name) // ' ', bd_list(i)%rank, bd_list(i)%ntlev
    ENDIF
    IF ( ASSOCIATED( bd_list(i)%p2 ) ) THEN
       CALL acc_update_device( bd_list(i)%p2 )
    ELSEIF ( ASSOCIATED( bd_list(i)%p3 ) ) THEN
       CALL acc_update_device( bd_list(i)%p3 )
    ELSE
       ierr = 9042
       yerrmsg  = ' ERROR *** Neither bd_list(i)%p2 nor bd_list(i)%p3 is associated  ***'
       RETURN
    END IF
  END DO
#endif

!------------------------------------------------------------------------------
! End of module procedure acc_update_bd_fields
!------------------------------------------------------------------------------

END SUBROUTINE acc_update_bd_fields

!==============================================================================
  
!==============================================================================
!+ Module procedure "acc_update_device" in "acc_utilities"
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!
! Description:
!
! The subroutines below are used for updating pointers in structures since this
! is currently not supported by the OpenACC standard.
!
! Method:
!
! Using OpenACC update directives
!------------------------------------------------------------------------------

SUBROUTINE acc_update_device_3d( fld )
  
  REAL (KIND=wp) :: fld(:,:,:)

#ifdef _OPENACC
  IF ( ( my_cart_id == 0 ) .AND. ( idbg_level > 20 ) ) THEN
    WRITE(*,*) ' GPUINFO: CPU-GPU copy, one 3d field'
  ENDIF
#endif

  !$acc data present(fld)
  !$acc update device(fld)
  !$acc end data

END SUBROUTINE acc_update_device_3d

SUBROUTINE acc_update_device_2d( fld )
  
  REAL (KIND=wp) :: fld(:,:)

#ifdef _OPENACC
  IF ( ( my_cart_id == 0 ) .AND. ( idbg_level > 20 ) ) THEN
    WRITE(*,*) ' GPUINFO: CPU-GPU copy, one 2d field'
  ENDIF
#endif

  !$acc data present(fld)
  !$acc update device(fld)
  !$acc end data

END SUBROUTINE acc_update_device_2d

!==============================================================================
  
!==============================================================================
!+ Module procedure "acc_create_device" in "acc_utilities"
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!
! Description:
!
! The subroutines below are used for allocating arrays in structures since this
! is currently not supported by the OpenACC standard.
!
! Method:
!
! Using OpenACC update directives
!------------------------------------------------------------------------------

SUBROUTINE acc_copyin_1d( fld )
  
  REAL (KIND=wp) :: fld(:)

  !$acc enter data copyin(fld)

END SUBROUTINE acc_copyin_1d

SUBROUTINE acc_delete_1d( fld )
  
  REAL (KIND=wp) :: fld(:)

  !$acc exit data delete(fld)

END SUBROUTINE acc_delete_1d

!==============================================================================

END MODULE acc_utilities
