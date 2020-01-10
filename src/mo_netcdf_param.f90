!
!+ 3DVAR/COSMO module to interface the NetCDF include file 'netcdf.inc'
!
! $Id: mo_netcdf_param.f90,v 5.0 2013-11-08 10:08:04 for0adm Exp $
!-------------------------------------------------------------------------------
!
MODULE mo_netcdf_param
!
!-------------------------------------------------------------------------------
! Description:
!   This module provides acess to the parameters defined in the
!   file 'netcdf.inc' of the NetCDF package
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
!------------------------------------------------------------------------------

public

include 'netcdf.inc'

end module mo_netcdf_param
