!+ Source Module providing utility routines for grib I/O
!==============================================================================

MODULE io_utilities

!==============================================================================
!
! Description:
!  This module provides all routines that have to deal with input and output
!  of grib files. Some of these routines use the message-passing library MPI. 
!  If possible the routines are written plug-compatible (with the exception
!  of setup_io).
!
!  Routines (module procedures) currently contained:
!
!    - make_fn
!       creates file names for binary data and for ready files
!
!    - open_file
!       opens a file with routines from special libraries (Grib, NetCDF, ...)
!
!    - compute_grib_intbuffer_length
!       estimates the length of the buffer necessary to hold a Grib record
!
!    - close_file
!       closes a file with routines from special libraries (Grib, NetCDF, ...)
!
!    - read_grib
!       reads records from a grib file using dwdlib and distributes them to the PEs
!
!    - read_gribapi
!       reads records from a grib file using grib_api and distributes them to the PEs
!
!    - read_netcdf
!       reads records from a NetCDF file and distributes them to the PEs
!
!    - read_restart
!       reads records from a restart file and distributes them to the PEs
!
!    - write_grib
!       collects gribed records from the processors and writes them to disk.
!
!    - write_gribapi
!       collects gribed records from the processors and writes them to disk.
!
!    - write_netcdf
!       collects records from the processors and writes them to netcdf files.
!
!    - write_restart
!       collects records from the processors and writes them to a binary
!       restart file
!
!    - write_ready_final
!       creates a ready file; i.e. a small file indicating completion of write step
!
!    - check_record
!       prints MIN, MAX and meanvalues of all records to the file YUCHKDAT
!
!    - check_input_grid
!       checks the NAMELIST parameters against the meta data from the file
!       (unified now for all possible input models, except GME)
!
!    - check_gme_grid 
!       checks the NAMELIST parameters against the GME values from the input.
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Ulrich Schaettler
!  Initial release
! 1.4        1998/05/22 Guenther Doms
!  Adaptions for the two time-level integration scheme
! 1.7        1998/07/16 Guenther Doms
!  Removal of global array 'rrssk'.
! 1.8        1998/08/03 Ulrich Schaettler
!  Use grib parameters from module data_io.f90
! 1.10       1998/09/29 Ulrich Schaettler
!  Updated Use lists from the modules
! 1.17       1998/11/17 Ulrich Schaettler
!  Correct some ANSI violations.
! 1.20       1999/01/07 Guenhter Doms
!  Renaming of some global variables.
! 1.29       1999/05/11 Ulrich Schaettler
!  Adaptations to use this module also in GME2LM
! 1.32       1999/08/24 Guenther Doms
!  Correct a REAL-declaration
! 1.34       1999/12/10 Ulrich Schaettler
!  Changed routine check_record to use it also in GME2LM properly.
! 1.39       2000/05/03 Ulrich Schaettler
!  Adaptations to include asynchronous IO.
! 2.8        2001/07/06 Ulrich Schaettler
!  Checked and corrected the interfaces to mpe_io;
!  Introduced ymode_read and ymode_write;
!  In check_lm_grid and check_gme_grid now the actual date of the forecast is
!  checked instead of the initial date (or reference date).
! 2.9        2001/07/16 Ulrich Schaettler
!  Correction for reading Grib data.
! 2.12       2001/11/07 Ulrich Schaettler
!  Corrected the check for the pole of the rotated grid; adaptation to GME2LM
! 2.14       2002/02/15 Ulrich Schaettler
!  Modifications to SR check_lm_grid to allow the use of the SLEVE coordinate
! 2.15       2002/03/20 Ulrich Schaettler
!  Bug correction in check_lm_grid
! 2.17       2002/05/08 Ulrich Schaettler
!  Adaptations to perform communications for I/O in irealgrib-format;
!  in SR check_lm_grid, only bit 7 of igds(9) is tested for 0
! 3.6        2003/12/11 Ulrich Schaettler
!  Additional checks for "time range indicator" and "unit of time range"
!  in the Subroutines check_lm_grid and check_gme_grid
! 3.7        2004/02/18 Ulrich Schaettler
!  Adaptations in routine make_fn for checking the extension of file names
! 3.13       2004/12/03 Ulrich Schaettler
!  Adaptations to new version of INT2LM
! 3.15       2005/03/03 Ulrich Schaettler
!  Replaced FLOAT by REAL; Implemented unit_of_time in check-routines
! 3.16       2005/07/22 Ulrich Schaettler
!  Eliminated time check for T_CL and W_CL in case of LM2LM
! 3.18       2006/03/03 Ulrich Schaettler
!  Modifications to implement NetCDF and Restart I/O
!    (renamed and changed: open_grib, close_grib -> open_file, close_file
!     new: read_netcdf, read_restart, write_netcdf, write_restart)
!  New ytunit 'd' in SR make_fn (for CLM: file name with date as in laf-files)
!  New SR difmin_360: computes difference for 2 dates using a year with 360 days
!  Added parameter lyear_360 in SRs check_*_grid to choose proper SR difmin
!  Changes in SR check_lm_grid for new type of coding the vertical coordinates
!    (introduced by MeteoSwiss for SLEVE vertical coordinate)
! 3.19       2006/04/25 Ulrich Schaettler
!  In read_netcdf: skip variables in the initial list, that are not present
! 3.21       2006/12/04 Ulrich Schaettler
!  In case of ymode_write=append, restart files are opened as NEW
!  Adapt check_record to use of undefined values in netcdf
!  Bug corrections in read_restart, write_restart
! V3_23        2007/03/30 Ulrich Schaettler
!  Introduced idbg_level for verbosity of output
!  Changed SR make_fn to adapt the file name to quarterly hours, if necessary
! V3_26        2007/05/29 Ulrich Schaettler
!  Adaptations in SR make_fn for writing analysis files with flexible dt
! V4_1         2007/12/04 Ulrich Schaettler
!  Some editorial changes (print outs for debugging)
! V4_4         2008/07/16 Ulrich Schaettler, Burkhardt Rockel
!  Changed NL parameter lyear_360 to itype_calendar, to have several options
!  Allow 2D fields in netcdf input
! V4_5         2008/09/10 Guenther Zaengl
!  Adaptations for new reference atmosphere
! V4_8         2009/02/16 Guenther Zaengl
!  Change grib encoding for new reference atmosphere to facilitate future
!  extensions
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_17        2011/02/24 Ulrich Blahak
!  In subroutine make_fn, introduced parameter ztgranul to make the
!  allowed time distances of output files more flexible in case of lhour=.true.
!  Before, the allowed time distance was set fixed to a quarter hour.
! V4_18        2011/05/26 Ulrich Schaettler
!  Adapted SR read_netcdf to identify the 3D external parameter field for
!  topographical corrections (Anne Roches)
! V4_19        2011/08/01 Ulrich Schaettler
!  Introduced conditional compilation for NetCDF and GRIBDWD
! V4_20        2011/08/31 Ulrich Schaettler
!  Safer communication between IO and compute PEs for little endian problem with MPI_BYTE
!  (send 4 additional bytes)
! V4_23        2012/05/10 Ulrich Schaettler, CLM
!  Unifications with INT2LM Version 1.19
!   => usage of grib_api; new SR read_gribapi
!   => new (generic) SR check_input_grid, replacing check_[lm,ec,um]_grid
!   => changes to check_gme_grid because of grib_api usage
!   => modified interface to read_restart (has to be adapted in COSMO-Model calls)
!  nf90_clobber replaced by NF90_64BIT_OFFSET for large netCDF file support
! V4_24        2012/06/22 Ulrich Schaettler, Burkhardt Rockel, Hendrik Reich
!  Reading restart files: Adaptations are necessary in check_input_grid to deal
!   with binary restart files
! Adaptations for reading NetCDF files: Unification of IDs with INT2LM
!   Changed dim_id for topo corrections from 11 to 15 (by Burkhardt)
! Adapted date variables to longer character strings (Hendrik Reich)
! V4_25        2012/09/28 Ulrich Schaettler, Florian Prill
!  Modified make_fn to create the correct date string, if this has to be adapted
!   because of a non-fitting dt
! Florian Prill:
!  Removed database support
!  Adjusted calls to mpe_io2 interfaces
!  Added 2 new procedures: write_ready_final, compute_grib_intbuffer_length
! V4_26        2012/12/06 Ulrich Blahak
!  Bugfix: added missing keyword "lmmss" to write_ready_final, so that
!   it can be passed on to "make_fn" and ready files can be written in the
!   old 10-digit date format instead of hardwired 14 digits.
!  In make_fn, added "/" to the end of the pathname if necessary.
! V4_27        2013/03/19 Ulrich Schaettler
!  Adaptations to read ivctype from old style of gds-coding, so that it can
!  also be read as a real vertical coordinate parameter (as all other variables)
! V4_28        2013/07/12 Ulrich Schaettler
!  Unification with INT2LM Version (implemented api-blocks)
!  Modified check_input_grid for a more consistent checking of vertical
!   coordinate parameters for GRIB1
!  Adapted interfaces to grib_api routines with integer kind int_ga (which
!   could be 4- or 8-byte integer)
! V4_29        2013/10/04 Ulrich Schaettler
!  Rename vcoord_in to vcoord in the COSMO-Model
!  Before decoding vertical coordinate parameters with REFSTF, check the value
!  of the corresponding gds-entry, because some were not really coded
! V4_30        2013/11/08 Ulrich Schaettler
!  Corrected usage of uninitialized variable yzhead in SR write_ready_final
!  Modified check for pv_in in SR check_input_grid: only check vertical
!    coordinate parameters
! V5_1         2014-11-28 Ulrich Schaettler, Oliver Fuhrer, Jochen Foerstner
!  Enlarged check value for comparing pv-values to 5.0E-3
!  SR check_gme_grid: initialize ylevtyp
!  Added new argument nbytes to function compute_grib_intbuffer_length
!   (to indicate bytes per value for grib packing)
!  Technical adaptations to INT2LM Version
!  Replaced ireals by wp (working precision) (OF)
!  Changed format of reals for YUCHKDAT, because volcanic ash needs larger values (JF)
! V5_2         2015-05-21 Ulrich Schaettler
!  Changed limit for checking grib meta values of rotated pole to 1E-4, because
!   there were rounding problems with 1E-5
!  Unification with INT2LM version: modified make_fn, check_record, write_ready_final
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Modules used:
USE data_parameters , ONLY :   &
  wp,        & ! KIND-type parameters for real variables
  iintegers, & ! KIND-type parameter for standard integer variables
  irealgrib, & ! KIND-type parameter for the real variables in the grib library
  intgribf,  & ! KIND-type parameter for the fortran files of the grib library
  intgribc,  & ! KIND-type parameter for the c files of the grib library
  iwlength,  & ! length of a integer word of the grib library in byte
  int_ga       ! integer precision for grib_api: length of message in bytes

!==============================================================================

#ifdef GRIBAPI
! grib_api interface
USE grib_api
#endif

! declare routines used from the parallel IO interface MPE_IO
USE mpe_io2,            ONLY :   &
    mpe_io_open, mpe_io_read, mpe_io_write, mpe_io_close, mpe_io_complete,   &
    mpe_io_ready

USE utilities,          ONLY :   diff_minutes, get_utc_date
USE vgrid_refatm_utils, ONLY :   vcoord_type, uuid_match,                    &
                                 uuid_string_length, uuid_2char

!==============================================================================

#ifdef NETCDF
! declare netcdf variables and functions
USE netcdf,           ONLY :   &
  NF90_clobber,            &
  NF90_64BIT_OFFSET,       &
  NF90_close ,             &
  NF90_create,             &
  NF90_noerr,              &
  NF90_nofill,             &
  NF90_nowrite,            &
  NF90_open,               &
  NF90_inquire_variable,   &
  NF90_inquire_dimension,  &
  NF90_put_var,            &
  NF90_get_var,            &
  NF90_set_fill,           &
  NF90_strerror,           &
  NF90_write
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Include statements
  INCLUDE "mpif.h"

!==============================================================================

CONTAINS

!==============================================================================
!+ Creates the file name for grib-files or for ready-files
!------------------------------------------------------------------------------

SUBROUTINE make_fn (yhead, ydate, yinidate, ytunit, yexten, nstep, dt, lhour, &
                    icalendar, ydir, ypathname, lmmss, idbglv, ierror, icon_prefix)

!------------------------------------------------------------------------------
!
! Description:
!   *make_fn* creates the file names for the grib files or for ready-files
!   depending on yhead and ytunit.
!     Analysis:                 yhead//ydate_fn//yexten
!     Forecast and Boundary:    yhead//ytunit//yforrange//yexten
!
!   File name for forecasts ('f') or restart files ('r')
!   The form of the forecast range depends on "ytunit":
!     For ytunit = 't':  yforrange is the number of timesteps
!     For ytunit = 'f':  yforrange is given in the form "ddhhmmss" where
!                        dd: days, hh: hours, mm: minutes, ss: seconds
!     For ytunit = 'c':  yforrange is given in the form "yyydddhh" where
!                        yyy: years, ddd: days, hh: hours
!     For ytunit = 'd':  yforrange is given as ydate (for Climate LM Version)
!
! Method:
!   Necessary adaptations, if the actual time step used for I/O does not fit
!   to the chosen time for I/O:
!   Depending on the time step dt, it can happen, that the actual time step
!   does not have the correct time, but is the "nearest one" possible for the
!   chosen I/O step. Then the forecast range is adapted to the chosen I/O step
!   by using a fixed time granularity ztgranul (which is 900.0 seconds for the
!   moment being). Also the date is adapted to the chosen I/O step.
!
!------------------------------------------------------------------------------

!     Input
CHARACTER (LEN= 3),       INTENT(IN)  ::    &
  yhead        ! 3-Character header of the file, e.g.
               !   'gaf': analysis of GME for the full domain
               !   'laf': analysis of LM for the full domain
               !   'lbf': boundary for LM on the full domain
               !   'GME' or 'LMA': for ready files

CHARACTER (LEN=14),       INTENT(IN)  ::    &
  ydate,     & ! Initial date of the forecast in the form
               !   yyyymmddhh (yyyy: year, mm: month, dd: day, hh: time,
               !               mm: minute, ss: second)
  yinidate     ! Initial date, if the actual date has to be adapted

CHARACTER (LEN= 1),       INTENT(IN)  ::    &
  ytunit       ! 1-Character time unit for forecast range, e.g.
               !   't': timesteps, 'f': forecast mode, 
               !   'c': climate mode

CHARACTER (LEN= 1),       INTENT(IN)  ::    &
  yexten       ! 1-Character extension of file name, e.g.
               !   'c': constant files
               !   'z': data on z-levels
               !   'p': data on pressure levels
               !   's': data on "satellite" levels

CHARACTER (LEN= *),       INTENT(IN)  ::    &
  ydir         ! directory of the file

INTEGER (KIND=iintegers), INTENT(IN)  ::    &
  nstep,     & ! Current timestep number
  icalendar, & ! type of calendar
  idbglv       ! debug level for verbosity of output

REAL (KIND=wp),           INTENT(IN)  ::    &
  dt           ! Timestep of the model (in seconds)

LOGICAL,                  INTENT(IN)  ::    &
  lhour,     & ! if .TRUE., the filename has to be constructed as multiples
               ! of quarterly hours (only for ytunit = 'f' or 'c'
  lmmss        ! if .TRUE. explicit dates have 14 digits, 10 otherwise

CHARACTER (LEN=*),       INTENT(IN), OPTIONAL  ::    &
  icon_prefix  ! header of the file in case of ICON data. This replaces
               ! yhead in the filename if yhead(1:1) = 'i'

!------------------------------------------------------------------------------

! Output
CHARACTER (LEN= *),       INTENT(OUT) ::    &
  ypathname    ! Name of the file with full path name for access

INTEGER (KIND=iintegers), INTENT(OUT) ::    &
  ierror       ! Error flag; set to 0 if no error occured

!------------------------------------------------------------------------------

! Local variables
REAL (KIND=wp)                        ::    &
  zforrange, zforrange1, zforrange2, & ! Forecast range in seconds
  z2, zacthour

INTEGER (KIND=iintegers)              ::    &
  mfor_s,    & ! Forecast range (seconds)
  mfor_m,    & ! Forecast range (minutes)
  mfor_h,    & ! Forecast range (hours)
  mfor_d,    & ! Forecast range (days)
  mfor_y       ! Forecast range (years)

INTEGER (KIND=iintegers)              ::    &
  izh, nzjulianday

CHARACTER (LEN=14)                    ::    &
  yforrange, & ! Forecast range as character
  yzlocdate    ! if ydate has to be modified

CHARACTER (LEN=20)                    ::    &
  yfname       ! name of the file

CHARACTER (LEN=14)    ::           &
  yzakdat1    ! actual date (ydate_ini+ntstep/dt) in the form
              ! yyyymmddhhmmss (year, month, day, hour, minutes, seconds)
CHARACTER (LEN=28)    ::           &
  yzakdat2    ! actual date (ydate_ini+ntstep/dt) in the form
              ! wd dd.mm.yyyy  hh mm ss UTC  (weekday, ...)

REAL(KIND=wp),     PARAMETER          ::    &
  ztgranul = 900.0_wp    ! in s, 15 min for now.
       ! for lhour=.true.: 
       ! this parameter gives the granularity of the I/O time in seconds;
       ! which is rounded to the nearest multiple of ztgranul
       ! NOTE: ztgranul has to fit exactly into 3600 s, i.e., ztgranul might be,
       !       for example, something like 60.0, 300.0, 600.0, 900.0 or 1800.0

! End of header
!------------------------------------------------------------------------------

  ierror = 0
  yforrange = '              '
  zforrange2 = 0.0_wp

  IF (lmmss) THEN
    ! date and file names include minutes and seconds
    yzlocdate(1:14) = ydate (1:14)
  ELSE
    ! date and file names do not include minutes and seconds
    yzlocdate(1:14) = ydate (1:10)//'    '
  ENDIF

  zforrange1 = REAL(nstep, wp)*dt        ! forecast time in seconds
  zforrange  = zforrange1
  z2         = zforrange1 / ztgranul         ! forecast time in multiples of ztgranul

  IF (lhour) THEN
    IF (ABS(REAL(NINT(z2),wp) - z2) > 1.0E-5_wp) THEN
      ! zforrange has to be adapted to a multiple of ztgranul
      zforrange2 = REAL(NINT(z2),wp) * ztgranul
      zforrange  = zforrange2

      ! then also the date has to be adapted!!!
      CALL get_utc_date(NINT(z2), yinidate, ztgranul, icalendar,       &
                        yzlocdate, yzakdat2, nzjulianday, zacthour)

      IF (.NOT. lmmss) THEN
        ! date and file names do not include minutes and seconds
        yzlocdate(1:14) = yzlocdate (1:10)//'    '
      ENDIF

    ENDIF
  ENDIF

  IF (idbglv > 50) THEN
    PRINT *, '     in make_fn:  ', REAL(nstep, wp)*dt, z2, zforrange1, zforrange2, &
                                    zforrange, '  ', TRIM(yzlocdate)
  ENDIF

  ! Compute forecast range as character string depending on ytunit
  SELECT CASE (ytunit)
  CASE ('t')       ! forecast range in timesteps
    WRITE (yforrange,'(I8.8)') nstep

  CASE ('f')       ! forecast range in forecast-mode

    mfor_d    = INT ( zforrange/86400.0_wp, iintegers)
    IF (mfor_d > 99_iintegers) THEN
      ierror = 5
    ENDIF
    mfor_h    = INT ((zforrange - REAL(mfor_d, wp)*86400.0_wp) /    &
                                                   3600.0_wp, iintegers)
    mfor_m    = INT ((zforrange - REAL(mfor_d, wp)*86400.0_wp       &
                                - REAL(mfor_h, wp)* 3600.0_wp) /    &
                                                     60.0_wp, iintegers)
    mfor_s    = NINT( zforrange - REAL(mfor_d, wp)*86400.0_wp       &
                                - REAL(mfor_h, wp)* 3600.0_wp       &
                                - REAL(mfor_m, wp)*   60.0_wp,      &
                                                                  iintegers)
    WRITE (yforrange,'(4(I2.2))') mfor_d, mfor_h, mfor_m, mfor_s

  CASE ('c')

    IF     (icalendar == 0) THEN
      mfor_y    = INT ( zforrange / (365.0_wp*86400.0_wp), iintegers)
      mfor_d    = INT ((zforrange -                                           &
                          REAL (mfor_y,wp)*365.0_wp*86400.0_wp)               &
                                                  / 86400.0_wp, iintegers)
      mfor_h    = NINT((zforrange -                                           &
                          REAL (mfor_y, wp)*365.0_wp*86400.0_wp               &
                        - REAL (mfor_d, wp)*86400.0_wp) /                     &
                                                     3600.0_wp, iintegers)
    ELSEIF (icalendar == 1) THEN
      mfor_y    = INT ( zforrange / (360.0_wp*86400.0_wp), iintegers)
      mfor_d    = INT ((zforrange -                                           &
                          REAL (mfor_y,wp)*360.0_wp*86400.0_wp)               &
                                                  / 86400.0_wp, iintegers)
      mfor_h    = NINT((zforrange -                                           &
                          REAL (mfor_y, wp)*360.0_wp*86400.0_wp               &
                        - REAL (mfor_d, wp)*86400.0_wp) /                     &
                                                     3600.0_wp, iintegers)
    ELSEIF (icalendar == 2) THEN
      mfor_y    = INT ( zforrange / (365.0_wp*86400.0_wp), iintegers)
      mfor_d    = INT ((zforrange -                                           &
                          REAL (mfor_y,wp)*365.0_wp*86400.0_wp)               &
                                                  / 86400.0_wp, iintegers)
      mfor_h    = NINT((zforrange -                                           &
                          REAL (mfor_y, wp)*365.0_wp*86400.0_wp               &
                        - REAL (mfor_d, wp)*86400.0_wp) /                     &
                                                     3600.0_wp, iintegers)
    ENDIF
    IF (mfor_y > 999_iintegers) THEN
      ierror = 6
    ENDIF
    WRITE (yforrange,'(I3.3, I3.3, I2.2)') mfor_y, mfor_d, mfor_h

  CASE ('d')
    IF (lmmss) THEN
      ! date and file names include minutes and seconds
      yforrange(1:14) = yzlocdate (1:14)
    ELSE
      ! date and file names do not include minutes and seconds
      yforrange(1:14) = yzlocdate (1:10)//'    '
    ENDIF

  CASE DEFAULT
    ! Unknown file type or time unit, error exit
    IF (idbglv > 0) THEN
      PRINT *,'  Error in *make_fn*, unknown time unit,',' ytunit: ', ytunit
    ENDIF
    ierror = 1
  END SELECT

  ! put file name together: try to make a complete distinction of cases:
  SELECT CASE (yhead(1:1))

  CASE ('G','L','E','J','N','C')

    ! ready files for GME or COSMO
    SELECT CASE (yhead(3:3))

    CASE ('B','E','F','_')
      ! file name for forecast mode
      yfname = yhead//'_'//TRIM(yforrange)

    CASE ('A')
      ! file name for analysis mode
      yfname = yhead//'_'//TRIM(yzlocdate)

    CASE DEFAULT
      PRINT *,'  Error in *make_fn*, unknown third character to indicate &
                        &forecast or analysis mode:  ', yhead(1:1)
      ierror = 2
    END SELECT

  CASE ('I')
    ! ready files for ICON
    SELECT CASE (yhead(3:3))

    CASE ('L')
      ! file name for forecast mode:  IGLO_
      yfname = 'IGLO'//'_'//TRIM(yforrange)

    CASE ('A')
      ! file name for analysis mode:  IGA_
      yfname = yhead//'_'//TRIM(yzlocdate)

    CASE DEFAULT
      PRINT *,'  Error in *make_fn*, unknown third character to indicate &
                        &forecast or analysis mode:  ', yhead(1:1)
      ierror = 2
    END SELECT

  CASE ('g','l','e','c','j','n')
    ! model data files
    SELECT CASE (yhead(2:2))

    CASE ('i','b','f','r')
      ! file name for forecast mode
      IF (yexten == ' ') THEN
        yfname = yhead//ytunit//TRIM(yforrange)
      ELSE
        yfname = yhead//ytunit//TRIM(yforrange)//yexten
      ENDIF

    CASE ('a')
      ! file name for analysis mode
      IF (yexten == ' ') THEN
        yfname = yhead//TRIM(yzlocdate)
      ELSE
        yfname = yhead//TRIM(yzlocdate)//yexten
      ENDIF

    CASE DEFAULT
      PRINT *,'  Error in *make_fn*, unknown second character to indicate &
                        &forecast or analysis mode:  ', yhead(1:1)
      ierror = 2
    END SELECT

  CASE ('i')
    ! model data files for ICON: here the filename prefix has been
    !  specified via namelist and is overtaken from the optional input string
    !  "icon_prefix":

    IF (.NOT. PRESENT(icon_prefix)) THEN
      yfname = 'noprefix'//TRIM(yforrange)//TRIM(yexten)
      PRINT *,'  Error in *make_fn*: icon_prefix not specified! '
      ierror = 2
    ELSE
      yfname = TRIM(icon_prefix)//TRIM(yforrange)//TRIM(yexten)
    END IF

!!$ UB: this was the old method without "icon_prefix". We leave it in the code because
!!$     it is not clear how ICON analyses will be named, so we might revive it again later.

!!$    ! model data files for ICON: have 2 digits to indicate model, 
!!$    !    so we have to add another digit (only 'ig' possible at the moment)
!!$    SELECT CASE (yhead(2:2))
!!$
!!$    CASE ('i','b','f','r')
!!$      ! file name for forecast mode
!!$      IF (yexten == ' ') THEN
!!$        yfname = 'ig'//yhead(2:3)//ytunit//TRIM(yforrange)
!!$      ELSE
!!$        yfname = 'ig'//yhead(2:3)//ytunit//TRIM(yforrange)//yexten
!!$      ENDIF
!!$
!!$    CASE ('a')
!!$      ! file name for analysis mode
!!$      IF (yexten == ' ') THEN
!!$        yfname = 'ig'//yhead(2:3)//TRIM(yzlocdate)
!!$      ELSE
!!$        yfname = 'ig'//yhead(2:3)//TRIM(yzlocdate)//yexten
!!$      ENDIF
!!$
!!$    CASE DEFAULT
!!$      PRINT *,'  Error in *make_fn*, unknown second character to indicate &
!!$                        &forecast or analysis mode:  ', yhead(2:2)
!!$      ierror = 2
!!$    END SELECT

  CASE DEFAULT
    PRINT *,'  Error in *make_fn*, unknown first character to indicate model:  ', yhead(1:1)
    ierror = 2
  END SELECT

  ! Add the directory
  IF (LEN_TRIM(ydir) > 0) THEN
    ! Add "/" if necessary:
    IF (ydir(LEN_TRIM(ydir):LEN_TRIM(ydir)) /= '/') THEN
      ypathname = ydir(1:LEN_TRIM(ydir))//'/'//yfname
    ELSE
      ypathname = ydir(1:LEN_TRIM(ydir))//yfname
    END IF
  ELSE
    ypathname = yfname
  END IF

END SUBROUTINE make_fn

!==============================================================================
!==============================================================================
!+ Module procedure in "io_utilities" for opening a file
!------------------------------------------------------------------------------

SUBROUTINE open_file (nudat, datname, ymode, yformat, icomm, my_id, npes,   &
                      lasync, idbglv, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   Open a file for reading or writing. The calling routine has to specify
!   the mode of the file (read, write, append, additional (Grib) information),
!   the format of the file (Grib1, NetCDF, binary) and the filename.
!   The routine returns the file specifier (nudat). Only in case of binary
!   restart files, the file specifier is already determined before.
!   This is a collective routine, i.e. it must be called by all compute
!   processors.
!
! Method:
!   Calls to special libraries.
!
!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):

! Parameter list 
CHARACTER (LEN= *)       , INTENT(IN)    ::  &
  datname     ! filename to open

CHARACTER (LEN= 3)       , INTENT(IN)    ::  &
  ymode       ! access mode: read, write, append: 'r  ', 'w  ', ...

CHARACTER (LEN= 4)       , INTENT(IN)    ::  &
  yformat     ! format of the file (grb1, ncdf, bina)

INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
  icomm,    & ! MPI communicator
  my_id,    & ! id in communicator icomm
  npes,     & ! number of PEs
  idbglv      ! for verbosity of debug output

! Scalar arguments with intent(out):
INTEGER  (KIND=iintegers), INTENT(INOUT) ::  &
  nudat    ! internal file descriptor

LOGICAL                  , INTENT(IN)    ::  &
  lasync   ! indicates whether asynchronous IO is used or not

CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
  yerrmsg     ! string for error messages

INTEGER  (KIND=iintegers), INTENT(OUT)   ::  &
  ierror      ! error status

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER  (KIND=iintegers)           :: &
  joldMode,                            & ! for netcdf routines
  implcode                               ! Error variable for MPI

INTEGER  (KIND=intgribc)            :: &
  nudatc,                              & ! file descriptor for C-routine
  ierrc                                  ! error code for C-routine

LOGICAL                             :: &
  stop_dummy                             ! for calling mpe_io_open

CHARACTER (LEN= 3)                  ::  &
  yzmode      ! modified access mode for mpe_io_open

!- End of header
!------------------------------------------------------------------------------

ierror   = 0
yerrmsg  = '         '
ierrc    = 0_intgribc
nudatc   = 0_intgribc

IF (my_id == 0 .AND. idbglv > 1)                                           &
   PRINT *,'OPEN: ',yformat,'-file: ', datname(1:LEN_TRIM(datname))

SELECT CASE (yformat)

! At the moment, asynchronous IO is possible only for grib1-files
! I/O for all other formats is handled by processor 0

#ifdef GRIBDWD
CASE ('grb1')

  IF (lasync .OR. (npes > 1) ) THEN
    ! Open file using the parallel IO interface
    yzmode = ymode
    yzmode(2:2) = 'g'
    CALL mpe_io_open (nudat, datname(1:LEN_TRIM(datname)), yzmode,         &
                      stop_dummy, ierror)
    IF (ierror /= 0) THEN
      yerrmsg = 'Error opening '//TRIM(datname)
      ! this error message has to be broadcasted to all other PEs
      ierror  = 3
    ENDIF
  ELSE
    ! Open file using griblib calls
    IF (my_id == 0) THEN
      CALL copen(nudatc,datname(1:LEN_TRIM(datname)),ymode,ierrc)
      nudat = nudatc
      IF (ierrc /= 0) THEN
        yerrmsg = 'Error opening '//datname
        ! this error message has to be broadcasted to all other PEs
        ierror  = 3
      ENDIF
    ENDIF
  ENDIF
#endif

#ifdef GRIBAPI
CASE ('apix','api1','api2')

  IF (lasync .OR. (npes > 1) ) THEN
    ! Open file using the parallel IO interface
    yzmode = ymode
    IF     (yformat == 'api1') THEN
      yzmode(2:2) = '1'
    ELSEIF (yformat == 'api2') THEN
      yzmode(2:2) = '2'
    ELSEIF (yformat == 'apix') THEN
      yzmode(2:2) = 'x'
    ENDIF
    CALL mpe_io_open (nudat, datname(1:LEN_TRIM(datname)), yzmode,         &
                      stop_dummy, ierror)
    IF (ierror /= 0) THEN
      yerrmsg = 'Error opening '//datname
      ! this error message has to be broadcasted to all other PEs
      ierror  = 3
    ENDIF
  ELSE
    ! Open file using grib-api calls
    IF (my_id == 0) THEN

      IF (ymode(1:1) == 'r') THEN

        CALL grib_open_file(nudat, TRIM(datname), 'r', ierror)
        IF (ierror /= GRIB_SUCCESS) THEN
          yerrmsg = 'Error opening '//datname
          ! this error message has to be broadcasted to all other PEs
        ENDIF

      ELSEIF (ymode(1:1) == 'w') THEN

        CALL grib_open_file(nudat, TRIM(datname), 'w', ierror)
        IF (ierror /= GRIB_SUCCESS) THEN
          yerrmsg = 'Error opening '//datname
          ! this error message has to be broadcasted to all other PEs
        ENDIF

      ENDIF

    ENDIF

  ENDIF
#endif

#ifdef NETCDF
CASE ('ncdf')

  IF (my_id == 0) THEN
    IF (ymode(1:1) == 'r') THEN
  
      ierror = nf90_open(TRIM(datname), nf90_nowrite, nudat)
      IF (ierror /= nf90_noerr) THEN
        PRINT *, TRIM(NF90_strerror(ierror))
      ENDIF
  
    ELSEIF (ymode(1:1) == 'w') THEN
  
      ierror = nf90_create(TRIM(datname), IOR(nf90_clobber,nf90_64bit_offset), nudat)
      ierror = nf90_set_fill(nudat, NF90_NOFILL, joldMode)
      IF (ierror /= nf90_noerr) THEN
        PRINT *, TRIM(NF90_strerror(ierror))
      ENDIF
  
    ELSEIF (ymode(1:1) == 'a') THEN
  
      ierror = nf90_open(TRIM(datname), nf90_write, nudat)
      ierror = nf90_set_fill(nudat, NF90_NOFILL, joldMode)
      IF (ierror /= nf90_noerr) THEN
        PRINT *, TRIM(NF90_strerror(ierror))
      ENDIF
  
    ELSE
      PRINT *, ' *** ERROR:  Wrong mode for netcdf-format: ', ymode(1:1), &
               ' *** '
      ierror = -1
    ENDIF
  ENDIF
#endif

CASE ('bina')
 
  IF (my_id == 0) THEN
    IF (ymode(1:1) == 'r') THEN

      OPEN (nudat, FILE=TRIM(datname), FORM='UNFORMATTED', STATUS='OLD', &
                   ACTION='READ', IOSTAT=ierror)
    ELSEIF ((ymode(1:1) == 'w') .OR. (ymode(1:1) == 'a')) THEN

      OPEN (nudat, FILE=TRIM(datname), FORM='UNFORMATTED', STATUS='UNKNOWN', &
                   ACTION='WRITE', IOSTAT=ierror)

    ELSE
      ierror = -1
    ENDIF
  ENDIF

CASE DEFAULT

  IF (my_id == 0) THEN
    ierror = -10
    PRINT *, ' *** ERROR: File format ', yformat, ' is not known to LM! ***'
  ENDIF

END SELECT


! The error message has to be broadcasted to all other PEs
IF (npes > 1) THEN
  CALL MPI_BCAST (ierror, 1, MPI_INTEGER, 0, icomm, implcode)
  IF (ierror /= 0) RETURN
ENDIF

END SUBROUTINE open_file

!==============================================================================
!==============================================================================
!+ Module procedure in "io_utilities" for closing a file
!------------------------------------------------------------------------------

SUBROUTINE close_file (nudat, yformat, icomm, my_id, npes, lasync,    &
                       idbglv, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   Close a file. The calling routine has to specify the format of the file.
!   Depending on this format, the necessary closing routines are called.
!   At the moment, asynchronous IO is implemented only for the Grib1 format.
!   This is a collective routine, i.e. it must be called by all processors
!
! Method:
!   Calls to special libraries.
!
!------------------------------------------------------------------------------
!
! Scalar arguments with intent(in):
INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
  nudat, & ! internal file descriptor
  icomm, & ! MPI communicator
  my_id, & ! id in communicator icomm
  npes,  & ! number of PEs
  idbglv   ! for verbosity of output

CHARACTER (LEN= 4)       , INTENT(IN)    ::  &
  yformat     ! format of the file (grib1, netcdf, binary)

LOGICAL                  , INTENT(IN)    ::  &
  lasync   ! indicates whether asynchronous IO is used or not

! Scalar arguments with intent(out):
CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
  yerrmsg     ! string for error messages

INTEGER  (KIND=iintegers), INTENT(OUT)   ::  &
  ierror      ! error status

!------------------------------------------------------------------------------

! Local scalars:
INTEGER  (KIND=iintegers)           :: &
  implcode                               ! Error variable for MPI

INTEGER  (KIND=intgribc)            ::  &
  nudatc,                               & ! file descriptor for C-routine
  ierrc                                   ! error code for C-routine

! End of Header
!------------------------------------------------------------------------------

ierror   = 0
implcode = 0
yerrmsg  = '         '

IF (my_id == 0 .AND. idbglv > 1) print *,'CLOSING ', yformat, ' FILE'

SELECT CASE (yformat)

! At the moment, asynchronous IO is possible only for grib1-files
! I/O for all other formats is handled by processor 0

#ifdef GRIBDWD
CASE ('grb1')

  IF (lasync .OR. (npes > 1) ) THEN
    ! Close file using the parallel IO interface
    CALL mpe_io_close (nudat, implcode)
  ELSE
    ! Close file using griblib calls
    IF (my_id == 0) THEN
      nudatc = nudat
      CALL cclose(nudatc,'exi',ierrc)
      implcode = INT (ierrc, iintegers)
    ENDIF
  ENDIF
#endif

#ifdef GRIBAPI
CASE ('apix','api1','api2')

  IF (lasync .OR. (npes > 1) ) THEN
    ! Close file using the parallel IO interface
    CALL mpe_io_close(nudat, implcode)
  ELSE
    ! Close file using grib-api calls
    IF (my_id == 0) THEN
      CALL grib_close_file(nudat)
    ENDIF
  ENDIF
#endif

#ifdef NETCDF
CASE ('ncdf')

  IF (my_id == 0) THEN
    implcode = nf90_close(nudat)

    IF (implcode /= NF90_NOERR) THEN
      PRINT *, TRIM(NF90_strerror(implcode))
    ENDIF
  ENDIF
#endif

CASE ('bina')

  IF (my_id == 0) THEN
    CLOSE (nudat, IOSTAT=implcode)
  ENDIF

CASE DEFAULT

END SELECT

  IF (implcode /= 0) THEN
    yerrmsg = 'Error closing '//yformat//' file'
    ierror  = 3
    RETURN
  ENDIF

END SUBROUTINE close_file

!==============================================================================
!==============================================================================
!+ Module procedure in "io_utilities" for reading a grib record from a file
!------------------------------------------------------------------------------

SUBROUTINE read_grib (nudat, maxlen, ilfd, icomm, idata, ilen, npes,         &
                      lasync, ltime_barrier, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   Processor 0 reads all records and distributes them to the other processors. 
!   As many records are read from the grib file as there are processors 
!   present.
!
! Some Remarks:  (from R. Johanni, SGI/Cray):
!   "maxlen" and "ilen" are lengthes in BYTES !!!!
!   "data" should also be considered as a BYTE-array since all the I/O Routines
!   are dealing with BYTE-numbers only.
!   Also the GRIB1 standard explicitly only speaks about bytes and
!   never about integer words - GRIB1 files therefore should always
!   be considered as byte-streams!
!   "data" is declared as (KIND=intgribf) here, since more obvious data types
!   (as CHARACTER or INTEGER*1 would be) lead to numerous problems with
!   different compilers (esp. CRAY).
!
!------------------------------------------------------------------------------
!
! Scalar arguments with intent(in):
INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
  nudat,     & ! internal file descriptor
  ilfd,      & ! dimension for grib routines
  icomm,     & ! MPI communicator
  npes         ! number of PE

INTEGER  (KIND=int_ga),    INTENT(IN)    ::  &
  maxlen       ! max. length of data array

! Scalar arguments with intent(out):
INTEGER  (KIND=int_ga),    INTENT(OUT)   ::  &
  ilen         ! length of data actually read

! Array arguments with intent(out):
INTEGER  (KIND=intgribf ), TARGET, INTENT(OUT)   ::  &
  idata(:)     ! array to be filled

LOGICAL                  , INTENT(IN)    ::  &
  lasync,    & ! indicates whether asynchronous IO is used or not
  ltime_barrier! whether barrier call should be done

CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
  yerrmsg      ! string for error messages

INTEGER  (KIND=iintegers), INTENT(OUT)   ::  &
  ierror       ! error status

!------------------------------------------------------------------------------

! Local scalars
INTEGER  (KIND=iintegers)                ::  &
  ilr,           & ! record length read
  iw,            & ! for reading controlword
  implcode         ! Error variable for MPI

INTEGER  (KIND=intgribf )                ::  &
  ibuf(ilfd)

INTEGER  (KIND=intgribc)                 ::  &
  nudatc, maxlenc, ilenc, ierrc     ! corresponding variables for C-routines

!- End of header
!------------------------------------------------------------------------------

#ifdef GRIBDWD
ierrc    = 0_intgribc
ierror   = 0
yerrmsg  = '         '

! Safety first
IF ( maxlen > iwlength*ilfd ) THEN
  yerrmsg = 'lfd too small'
  ierror  = 1
  RETURN
ENDIF

! This is in all cases a collective Routine, therefore:
IF (npes > 1) THEN
  IF (ltime_barrier) THEN
    CALL MPI_BARRIER(icomm,implcode)
    IF ( implcode /= 0 ) THEN
      yerrmsg = 'MPI_BARRIER failed'
      ierror  = 2
      RETURN
    ENDIF
  ENDIF
ENDIF

IF (lasync .OR. (npes > 1) ) THEN
  CALL mpe_io_read(nudat, idata, ilen, ierror)
ELSE
  ! Read using DWD library calls
  nudatc  = INT (nudat, intgribc)
  maxlenc = INT (maxlen, intgribc)

  ! Read record
  CALL cuegin (nudatc, maxlenc, idata, ilenc, ierrc)
  ilen = INT (ilenc, int_ga)
  IF ( ierrc /= 0 ) THEN
    yerrmsg = 'Error in cuegin'
    ierror  = 5
    RETURN
  ENDIF
ENDIF

IF (npes > 1) THEN
  IF (ltime_barrier) THEN
    CALL MPI_BARRIER(icomm,implcode)
    IF ( implcode /= 0 ) THEN
      yerrmsg = 'Error in MPI_BARRIER'
      ierror  = 10
      RETURN
    ENDIF
  ENDIF
ENDIF
#endif

END SUBROUTINE read_grib

!==============================================================================
!==============================================================================
!+ Module procedure in "io_utilities" for reading a grib record from a file
!------------------------------------------------------------------------------

SUBROUTINE read_gribapi (nudat, maxlen, ilfd, icomm, idata, ilen, npes,     &
                         lasync, ltime_barrier, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   Processor 0 reads all records and distributes them to the other processors.
!   As many records are read from the grib file as there are processors
!   present.
!
! Some Remarks:  (from R. Johanni, SGI/Cray):
!   "maxlen" and "ilen" are lengthes in BYTES !!!!
!   "data" should also be considered as a BYTE-array since all the I/O Routines
!   are dealing with BYTE-numbers only.
!   Also the GRIB1 standard explicitly only speaks about bytes and
!   never about integer words - GRIB1 files therefore should always
!   be considered as byte-streams!
!   "data" is declared as (KIND=intgribf) here, since more obvious data types
!   (as CHARACTER or INTEGER*1 would be) lead to numerous problems with
!   different compilers (esp. CRAY).
!
!------------------------------------------------------------------------------
!
! Scalar arguments with intent(in):
INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
  nudat,     & ! internal file descriptor
  ilfd,      & ! dimension for grib routines
  icomm,     & ! MPI communicator
  npes         ! number of PE

INTEGER  (KIND=int_ga),    INTENT(IN)    ::  &
  maxlen       ! max. length of data array

! Scalar arguments with intent(out):
INTEGER  (KIND=int_ga),    INTENT(OUT)   ::  &
  ilen         ! length of data actually read

! Array arguments with intent(out):
INTEGER  (KIND=intgribf ), INTENT(OUT)   ::  &
  idata(ilfd)  ! array to be filled

LOGICAL                  , INTENT(IN)    ::  &
  lasync,    & ! indicates whether asynchronous IO is used or not
  ltime_barrier! whether barrier call should be done

CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
  yerrmsg      ! string for error messages

INTEGER  (KIND=iintegers), INTENT(OUT)   ::  &
  ierror       ! error status

!------------------------------------------------------------------------------

! Local scalars
INTEGER  (KIND=iintegers)                ::  &
  iw,            & ! for reading controlword
  implcode         ! Error variable for MPI

! Integers for grib_api interface should all be the same: here intgribf
INTEGER  (KIND=intgribf )                ::  &
  ierrf, nudatf, ibuf(3*maxlen)

INTEGER  (KIND=int_ga)                   ::  &
  ilen_ga, ilr

!------------------------------------------------------------------------------

#ifdef GRIBAPI
ierror   = 0
yerrmsg  = '         '
nudatf   = INT(nudat, intgribf)

! Safety first
IF ( maxlen > iwlength*ilfd ) THEN
  yerrmsg = 'lfd too small'
  ierror  = 1
  RETURN
ENDIF

! This is in all cases a collective Routine, therefore:
IF (npes > 1) THEN
  IF (ltime_barrier) THEN
    CALL MPI_BARRIER(icomm,implcode)
    IF ( implcode /= 0 ) THEN
      yerrmsg = 'MPI_BARRIER failed'
      ierror  = 2
      RETURN
    ENDIF
  ENDIF
ENDIF

IF (lasync .OR. (npes > 1) ) THEN
  CALL mpe_io_read(nudat, idata, ilen, ierror)
ELSE
  ! Read using grib-api calls
  ilen_ga = maxlen   ! ilen is intent(inout)
  CALL grib_read_from_file (nudatf, idata, ilen_ga, ierrf)
  ilen = INT (ilen_ga, int_ga)
  IF ( ierrf /= GRIB_SUCCESS ) THEN
    IF (ierrf /= GRIB_END_OF_FILE) THEN
      CALL grib_check(ierrf, 'read_from_file', '')
      yerrmsg = 'Error in grib_read_from_file'
      ierror  = 5
      RETURN
    ELSE
      ilen  = 0_int_ga
    ENDIF
  ENDIF
ENDIF

IF (npes > 1) THEN
  IF (ltime_barrier) THEN
    CALL MPI_BARRIER(icomm,implcode)
    IF ( implcode /= 0 ) THEN
      yerrmsg = 'Error in MPI_BARRIER'
      ierror  = 10
      RETURN
    ENDIF
  ENDIF
ENDIF

#else
ilen    = -1
ierror  = 11
yerrmsg = 'Model not compiled for grib_api'
#endif

END SUBROUTINE read_gribapi

!==============================================================================
!==============================================================================
!+ Module procedure in "io_utilities" for reading a grib record from a file
!------------------------------------------------------------------------------

SUBROUTINE read_netcdf (ncid, idim, jdim, ivar_count, ilev_count, ivar_id,   &
                        nvarin, idims, ndims, icomm, my_id, npes,            &
                        west_add_in, east_add_in, south_add_in, north_add_in,&
                        record, myvar, mylev, mylevtot, irsize, idattyp,     &
                        lasync, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   Processor 0 reads all records and distributes them to the other processors. 
!   As many records are read from the file as there are processors present.
!
!   In case of NetCDF IO, the reading can be done in an ordered way, because
!   it is known which variables are in the file and at which location. This
!   has been investigated by PE 0 in the subroutine read_nc_vdefs.
!   So the variables are read in the order in which they appear in the list
!   of the input variables (listin). ivar_count gives the index of the 
!   variable that has to be processed.
!
!------------------------------------------------------------------------------
!
! Scalar arguments with intent(in):
INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
  ncid,            & ! internal file descriptor
  idim, jdim,      & ! dimensions of the 2D record
  nvarin,          & ! number of elements in the input variable list
  ndims,           & ! dimension of idims
  icomm,           & ! MPI communicator
  my_id,           & ! ID of this PE in communicator icomm
  npes,            & ! number of PEs
  idattyp,         & ! type of data for MPI calls
  idims(ndims),    & ! for the dimension IDs
  ivar_id(nvarin), & ! array with the variable IDs
  west_add_in,     & !
  east_add_in,     & !
  south_add_in,    & !
  north_add_in

! Scalar arguments with intent(inout):
INTEGER  (KIND=iintegers), INTENT(INOUT) ::  &
  ivar_count,  & ! actual index in the variable list
  ilev_count     ! actual level of a 3D variable

! Array arguments with intent(out):
INTEGER  (KIND=intgribf ), INTENT(OUT)   ::  &
  myvar,       & ! index of the variable this PE has to process (in listin)
  mylev,       & ! level of the 3D variable this PE has to process
  mylevtot,    & ! total levels of the variable this PE has to process
  irsize         ! length of data actually read

REAL     (KIND=irealgrib), INTENT(OUT)   ::  &
  record(idim*jdim) ! data to be read

LOGICAL                  , INTENT(IN)    ::  &
  lasync         ! indicates whether asynchronous IO is used or not

CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
  yerrmsg        ! string for error messages

INTEGER  (KIND=iintegers), INTENT(OUT)   ::  &
  ierror         ! error status

!------------------------------------------------------------------------------

! Local scalars
INTEGER  (KIND=iintegers)                ::  &
  ipr,           & ! loop counter
  iedims, jedims,& ! dimensions read from the NetCDF file
  ilen,          & ! length of a record
  ilenbuf,       & ! length of the buffer to receive data
  ndimids,       & ! number of dimension id's
  dimsids(4),    & ! for NetCDF dimension id's
  klev,          & ! number of vertical levels of a variable
  implcode         ! Error variable for MPI

INTEGER  (KIND=iintegers)                ::  &
  istatus(MPI_STATUS_SIZE)                     ! for status in MPI calls

INTEGER (KIND=iintegers)                ::  &
  i, j, ij   ! loop variables

REAL (KIND=irealgrib)                    ::  &
  rbuf(idim*jdim+3), reof, &
  rbuf2d(idim-east_add_in-west_add_in, jdim-south_add_in-north_add_in)

!- End of header
!------------------------------------------------------------------------------

#ifdef NETCDF
ierror   = 0
yerrmsg  = '         '
reof     = -999999.0_irealgrib
ilen     = (idim-east_add_in-west_add_in) * (jdim-south_add_in-north_add_in)
ilenbuf  = ilen + 3

IF(my_id == 0) THEN

  ! Processor 0 reads the data and sends it to all others
  ! Every processor gets a 2D record (as long as there are enough records)

  processor_loop: DO ipr = 0, npes-1

  IF (ivar_count <= nvarin) THEN
    increase_loop: DO WHILE (ivar_id(ivar_count) == -1)
      ivar_count = ivar_count + 1
      IF (ivar_count > nvarin) THEN
        ! we are at the end
        EXIT increase_loop
      ENDIF
    ENDDO increase_loop
  ENDIF

  IF (ivar_count <= nvarin) THEN

    ! there is still something to read
    ! get number of dimension id's for variable ivar_count
    !  ndimids = 2: means 3 dimensions; a 2D space array
    !  ndimids = 3: means 3 dimensions; a 2D space array + time
    !  ndimids = 4: means 4 dimensions; a 3D space array + time
    ierror = NF90_inquire_variable (ncid, ivar_id(ivar_count), ndims=ndimids)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      PRINT *, '*** ERROR in read_netcdf:  ', ierror, TRIM(yerrmsg)
      RETURN
    ENDIF

    dimsids(:) = 0_iintegers

    ! get the values of the dimension ids: iedims, jedims, klev
    ierror = NF90_inquire_variable (ncid, ivar_id(ivar_count),     &
                                                 dimids = dimsids(1:ndimids))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      PRINT *, '*** ERROR in read_netcdf:  ', ierror, TRIM(yerrmsg)
      RETURN
    ENDIF

    ierror = NF90_inquire_dimension (ncid, dimsids(1), len=iedims)
    ierror = NF90_inquire_dimension (ncid, dimsids(2), len=jedims)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      PRINT *, '*** ERROR in read_netcdf:  ', ierror, TRIM(yerrmsg)
      RETURN
    ENDIF
    IF ( (iedims /= idim-west_add_in-east_add_in) .OR. &
         (jedims /= jdim-south_add_in-north_add_in) ) THEN
      ierror  = 10
      yerrmsg = '  got wrong dimensions from NetCDF file ***'
      PRINT *, '*** ERROR in read_netcdf:  ', ierror, TRIM(yerrmsg)
      PRINT *, '      iedim = ', iedims, ';    jedim = ', jedims
      RETURN
    ENDIF

    ! read the data, depending on its dimension in space
    SELECT CASE (ndimids)

    CASE (2) ! this is a 2D field
      klev = 1
      ierror = NF90_get_var (ncid, ivar_id(ivar_count), rbuf2d,     &
                             start=(/1,1/), count=(/iedims, jedims/))
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        PRINT *, '*** ERROR in read_netcdf:  ', ierror, TRIM(yerrmsg)
        RETURN
      ENDIF

    CASE (3) ! this is a 3D field
      ierror = NF90_inquire_dimension (ncid, dimsids(3), len=klev)
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        PRINT *, '*** ERROR in read_netcdf:  ', ierror, TRIM(yerrmsg)
        RETURN
      ENDIF
      ilev_count = ilev_count + 1
      ierror = NF90_get_var (ncid, ivar_id(ivar_count), rbuf2d,     &
                        start=(/1,1,ilev_count/), count=(/iedims, jedims, 1/))
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        PRINT *, '*** ERROR in read_netcdf:  ', ierror, TRIM(yerrmsg)
        RETURN
      ENDIF

    CASE (4) ! this is a 4D field
      ierror = NF90_inquire_dimension (ncid, dimsids(3), len=klev)
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        PRINT *, '*** ERROR in read_netcdf:  ', ierror, TRIM(yerrmsg)
        RETURN
      ENDIF
      ilev_count = ilev_count + 1
      ierror = NF90_get_var (ncid, ivar_id(ivar_count), rbuf2d,     &
                    start=(/1,1,ilev_count,1/), count=(/iedims,jedims,1,1/))
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        PRINT *, '*** ERROR in read_netcdf:  ', ierror, TRIM(yerrmsg)
        RETURN
      ENDIF

    END SELECT

    irsize = ilen

    DO j = 1, jedims
      DO i = 1, iedims
        ij = (j-1)*iedims + i
        rbuf(ij) = rbuf2d(i,j)
      ENDDO
    ENDDO

    ! Distribute the record to other processors (or keep it for yourself)
    IF (ipr == 0) THEN
      record(1:ilen) = rbuf(1:ilen)
      myvar    = ivar_count
      IF (ilev_count == 0) THEN
        ! this is a 2D array
        mylev    = 1
      ELSE
        ! this is a 3D array
        mylev    = ilev_count
      ENDIF
      mylevtot = klev
    ELSE
      rbuf(ilen + 1) = REAL(ivar_count, irealgrib)
      IF (ilev_count == 0) THEN
        ! this is a 2D array
        rbuf(ilen + 2) = 1.0_irealgrib
      ELSE
        ! this is a 3D array
        rbuf(ilen + 2) = REAL(ilev_count, irealgrib)
      ENDIF
      rbuf(ilen + 3) = REAL(klev      , irealgrib)
      CALL MPI_SEND(rbuf, ilen+3, idattyp, ipr, 1002, icomm, implcode)
      IF ( implcode /= 0 ) THEN
        yerrmsg = 'Error in MPI_SEND record'
        ierror  = 7
        RETURN
      ENDIF
    ENDIF

    ! reset ilev_count and increase ivar_count, if necessary
    IF (ilev_count == klev) THEN
      ! reset ilev_count
      ilev_count = 0
    ENDIF
    IF (ilev_count == 0) THEN
      ! either a 2D field has been read or all levels of a 3D field
      ivar_count = ivar_count + 1
    ENDIF

  ELSE  !  (ivar_count > nvarin)

    ! the end of the file is reached: send reof to the rest of the processors
    IF (ipr == 0) THEN
      irsize = 0
    ELSE
      CALL MPI_SEND (reof, 1, idattyp, ipr, 1002, icomm, implcode)
      IF ( implcode /= 0 ) THEN
        yerrmsg = 'Error in MPI_SEND record'
        ierror  = 7
        RETURN
      ENDIF
    ENDIF

  ENDIF

  ENDDO processor_loop

ELSE   ! (my_id /= 0)

  ! All other Processors just receive the data
  CALL MPI_RECV(rbuf, ilenbuf, idattyp, 0, 1002, icomm, istatus, implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'Error in MPI_RECV record'
    ierror  = 9
    RETURN
  ENDIF

  IF (rbuf(1) /= reof) THEN
    record(1:ilen) = rbuf(1:ilen)
    myvar          = NINT(rbuf(ilen + 1),intgribf)
    mylev          = NINT(rbuf(ilen + 2),intgribf)
    mylevtot       = NINT(rbuf(ilen + 3),intgribf)
    irsize   = ilen
  ELSE
    irsize   = 0
  ENDIF

ENDIF
#endif

END SUBROUTINE read_netcdf

!==============================================================================
!==============================================================================
!+ Module procedure in "io_utilities" for reading a grib record from a file
!------------------------------------------------------------------------------

SUBROUTINE read_restart (nudat, idim, jdim, record, ipds, npds, igds, ngds, &
                         icomm, my_id, npes, irsize, irealtyp, igribtyp,    &
                         lasync, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!
!------------------------------------------------------------------------------

INTEGER (KIND=iintegers),  INTENT(IN)     :: &
  nudat,       & ! unit number
  idim, jdim,  & ! dimensions of the 2D record
  icomm,       & ! MPI communicator
  my_id,       & ! ID of this PE in communicator icomm
  npes,        & ! number of PEs
  irealtyp,    & ! type of real data for MPI calls
  igribtyp       ! type of real data for MPI calls

INTEGER (KIND=intgribf),   INTENT(IN)     :: &
  npds, ngds     ! dimension of pds and gds

INTEGER (KIND=intgribf),   INTENT(OUT)    :: &
  ipds(npds),  & ! pds: product definition section
  igds(ngds)     ! pds: product definition section

INTEGER (KIND=iintegers),  INTENT(OUT)    :: &
  irsize         ! length of data actually read

REAL     (KIND=wp),        INTENT(OUT)   ::  &
  record(idim*jdim) ! data to be read

LOGICAL                  , INTENT(IN)    ::  &
  lasync         ! indicates whether asynchronous IO is used or not

CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
  yerrmsg        ! string for error messages

INTEGER  (KIND=iintegers), INTENT(OUT)   ::  &
  ierror         ! error status

!------------------------------------------------------------------------------

! Local scalars
INTEGER  (KIND=iintegers)                ::  &
  ipr,           & ! loop counter
  ilen, ilr,     & ! length of a record
  ilenbuf,       & ! length of the buffer to receive data
  izerr,         & ! local error variable
  implcode         ! Error variable for MPI

INTEGER  (KIND=iintegers)                ::  &
  istatus(MPI_STATUS_SIZE)                     ! for status in MPI calls

REAL (KIND=wp)                           ::  &
  rbuf(idim*jdim), reof

INTEGER(KIND=intgribf)                   ::  &
  ipdsbuf(npds), igdsbuf(ngds)

LOGICAL                                  ::  &
  lzeof            ! for end-of-file

!- End of header
!------------------------------------------------------------------------------

ierror   = 0
yerrmsg  = '         '
reof     = -999999.0_wp
izerr    = 0_iintegers
implcode = 0_iintegers
ilen     = idim * jdim
ilenbuf  = idim * jdim

IF(my_id == 0) THEN

  ! Processor 0 reads the data and sends it to all others
  ! Every processor gets a 2D record (as long as there are enough records)

  lzeof = .FALSE.

  processor_loop: DO ipr = 0, npes-1

    IF (.NOT. lzeof) THEN
      READ (nudat, IOSTAT=izerr) ipdsbuf, igdsbuf, rbuf

      IF (izerr < 0) THEN
        ! this is the end-of-file condition
        lzeof  = .TRUE.
        ilr    = 0
        IF (ipr > 0) THEN
          ! send reof to processor ipr
          CALL MPI_SEND (reof, 1, irealtyp, ipr, 1002, icomm, implcode)
          IF ( implcode /= 0 ) THEN
            yerrmsg = 'Error in MPI_SEND record'
            ierror  = 2
            RETURN
          ENDIF
        ELSE
          ! set irsize for processor 0 to eof:
          irsize = ilr
        ENDIF
      ELSEIF (izerr > 0) THEN
        ! an error occured
        ierror  = 1
        yerrmsg = 'error while reading restart file'
        RETURN
      ELSE  
        ! no error and no end-of-file occured

        ! Distribute the record to other processors (or keep it for yourself)
        IF (ipr == 0) THEN
          record(1:ilen) = rbuf(1:ilen)
          ipds(1:npds)   = ipdsbuf(1:npds)
          igds(1:ngds)   = igdsbuf(1:ngds)
          irsize         = ilen
        ELSE
          CALL MPI_SEND(rbuf,    ilen, irealtyp, ipr, 1002, icomm, implcode)
          CALL MPI_SEND(ipdsbuf, npds, igribtyp, ipr, 1003, icomm, implcode)
          CALL MPI_SEND(igdsbuf, ngds, igribtyp, ipr, 1004, icomm, implcode)
          IF ( implcode /= 0 ) THEN
            yerrmsg = 'Error in MPI_SEND record'
            ierror  = 7
            RETURN
          ENDIF
        ENDIF
      ENDIF

    ELSE

      ! send reof to all other processors
      IF (ipr /= 0) THEN
        CALL MPI_SEND (reof, 1, irealtyp, ipr, 1002, icomm, implcode)
        IF ( implcode /= 0 ) THEN
          yerrmsg = 'Error in MPI_SEND record'
          ierror  = 3
          RETURN
        ENDIF
      ENDIF

    ENDIF

  ENDDO processor_loop

ELSE   ! (my_id /= 0)

  ! All other Processors just receive the data
  CALL MPI_RECV(record, ilenbuf, irealtyp, 0, 1002, icomm, istatus, implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'Error in MPI_RECV record'
    ierror  = 9
    RETURN
  ENDIF

  IF (record(1) /= reof) THEN
    irsize = ilen
    ! also receive the pds and gds
    CALL MPI_RECV(ipds, npds, igribtyp, 0, 1003, icomm, istatus, implcode)
    CALL MPI_RECV(igds, ngds, igribtyp, 0, 1004, icomm, istatus, implcode)
  ELSE
    irsize = 0
  ENDIF

ENDIF

END SUBROUTINE read_restart

!==============================================================================
!==============================================================================
!+ Module procedure in "io_utilities" for writing a grib file
!------------------------------------------------------------------------------

SUBROUTINE write_grib (nudat, idata, ilen, ilfd, icomm, npes,            &
                       lflush, lasync, ltime_barrier, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   Processor 0 gathers all fields from the other processors and writes them
!   to disk.
!
! Some Remarks:  (from R. Johanni, SGI/Cray):
!   "ilen" is in BYTES !!!!
!   "data" should also be considered as a BYTE-array since all the I/O Routines
!   are dealing with BYTE-numbers only.
!   Also the GRIB1 standard explicitly only speaks about bytes and
!   never about integer words - GRIB1 files therefore should allways
!   be considered as byte-streams!
!   "data" is declared as (KIND=intgribf) here, since more obvious data types
!   (as CHARACTER or INTEGER*1 would be) lead to numerous problems with
!   different compilers (esp. CRAY).
!
!------------------------------------------------------------------------------
!
! Subroutine arguments
! Scalar arguments with intent(in):
INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
  nudat,     & ! internal file descriptor
  ilfd,      & ! dimension for grib routines
  icomm,     & ! MPI communicator
  npes         ! number of compute PE

INTEGER  (KIND=int_ga),    INTENT(INOUT) ::  &
  ilen         ! length of data actually read

! Array arguments with intent(in):
INTEGER  (KIND=intgribf ), TARGET, INTENT(INOUT) ::  &
  idata(:)     ! array to be written

LOGICAL                  , INTENT(IN)    ::  &
  lasync,    & ! indicates whether asynchronous IO is used or not
  lflush,    & ! at the end of the writing
  ltime_barrier!

! Scalar arguments with intent(out):
CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
  yerrmsg      ! string for error messages

INTEGER  (KIND=iintegers), INTENT(OUT)   ::  &
  ierror       ! error status

!------------------------------------------------------------------------------

! Local scalars
INTEGER  (KIND=iintegers)                ::  &
  ilr,           & ! record length to write
  implcode         ! Error variable for MPI

INTEGER  (KIND=iintegers)                ::  &
  istatus(MPI_STATUS_SIZE)                     ! for status in MPI calls

INTEGER  (KIND=intgribf )                ::  &
  ibuf(ilfd)

INTEGER  (KIND=intgribc)                 ::  &
  nudatc, ilenc, ierrc       ! corresponding variables for C-routines

!- End of header
!------------------------------------------------------------------------------

#ifdef GRIBDWD
ierror   = 0
yerrmsg  = '         '

! Safety first
IF ( ilen > iwlength*ilfd ) THEN
  yerrmsg = 'lfd too small'
  ierror  = 1
  RETURN
ENDIF

! This is in all cases a collective Routine, therefore:
IF (npes > 1) THEN
  IF (ltime_barrier) THEN
    CALL MPI_BARRIER(icomm,implcode)
    IF ( implcode /= 0 ) THEN
      yerrmsg = 'MPI_BARRIER failed'
      ierror  = 2
      RETURN
    ENDIF
  ENDIF
ENDIF

nudatc = nudat

IF (lasync .OR. (npes > 1) ) THEN
  CALL mpe_io_write(idata, ilen, ierror)
  IF( lflush ) THEN
    CALL mpe_io_complete(ierror)
  ENDIF
ELSE
  ! this is a pure sequential program without any MPI
  ierrc = 0 ! should be set in cuegex, but is not !!!!!!
  ilenc = INT (ilen, intgribc)
  CALL cuegex (nudatc, idata, ilenc, ierrc)
  IF ( ierrc /= 0 ) THEN
    yerrmsg = 'Error in cuegex'
    ierror  = 5
    RETURN
  ENDIF
ENDIF

IF (npes > 1) THEN
  IF (ltime_barrier) THEN
    CALL MPI_Barrier(icomm,implcode)
    IF ( implcode /= 0 ) THEN
      yerrmsg = 'Error in MPI_BARRIER'
      ierror  = 11
      RETURN
    ENDIF
  ENDIF
ENDIF
#endif

END SUBROUTINE write_grib

!==============================================================================
!==============================================================================
!+ Module procedure in "io_utilities" for writing a grib file
!------------------------------------------------------------------------------

SUBROUTINE write_gribapi (nudat, ydata, ilen, ilfd, npes, icomm,        &
                          lflush, lasync, ltime_barrier, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   Processor 0 gathers all fields from the other processors and writes them
!   to disk using grib_write from grib_api.
!
! Some Remarks:  (from R. Johanni, SGI/Cray):
!   "ilen" is in BYTES !!!!
!   "data" is a 1-character array, which are just bytes
!
!------------------------------------------------------------------------------
!
! Subroutine arguments
! Scalar arguments with intent(in):
INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
  nudat,     & ! internal file descriptor
  ilfd,      & ! dimension for grib routines
  npes,      & ! number of PE
  icomm        ! MPI communicator

INTEGER  (KIND=int_ga),    INTENT(INOUT) ::  &
  ilen         ! length of data actually to write

! Array arguments with intent(in):
CHARACTER(LEN=1),                  INTENT(INOUT) ::  &
  ydata(ilfd)  ! array to be written (needs explicit length ilfd, otherwise
               ! compiler complains in connection with grib_write_bytes

LOGICAL                  , INTENT(IN)    ::  &
  lasync,    & ! indicates whether asynchronous IO is used or not
  lflush,    & ! at the end of the writing
  ltime_barrier!

! Scalar arguments with intent(out):
CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
  yerrmsg      ! string for error messages

INTEGER  (KIND=iintegers), INTENT(OUT)   ::  &
  ierror       ! error status

!------------------------------------------------------------------------------

! Local scalars
INTEGER  (KIND=iintegers)                ::  &
  isize_buffer,  & ! length of integer buffer
  implcode         ! Error variable for MPI

INTEGER  (KIND=intgribf ), ALLOCATABLE, TARGET  ::  &
  idata(:)     ! array to be written

INTEGER  (KIND=intgribc)                 ::  &
  nudatc, ierrc       ! corresponding variables for C-routines

!- End of header
!------------------------------------------------------------------------------

#ifdef GRIBAPI
ierror   = 0
yerrmsg  = '         '

! Safety first
IF ( ilen > ilfd ) THEN
  yerrmsg = 'lfd too small'
  ierror  = 1
  RETURN
ENDIF

! This is in all cases a collective Routine, therefore:
IF (npes > 1) THEN
  IF (ltime_barrier) THEN
    CALL MPI_BARRIER(icomm,implcode)
    IF ( implcode /= 0 ) THEN
      yerrmsg = 'MPI_BARRIER failed'
      ierror  = 2
      RETURN
    ENDIF
  ENDIF
ENDIF

nudatc = nudat

IF (lasync .OR. (npes > 1) ) THEN

! Somehow interpret ydata as idata!!!
  isize_buffer  = NINT (REAL(ilen/4, wp), iintegers) + 1
  ALLOCATE (idata(isize_buffer), STAT = ierror)
  IF (ilen > 0) THEN
    idata = TRANSFER (ydata, idata, isize_buffer)
    ! has to be re-transfered in mpe_io2
  ELSE
    idata(1) = 0
  ENDIF

  CALL mpe_io_write(idata, ilen, ierror)
  IF( lflush ) THEN
    CALL mpe_io_complete(ierror)
  ENDIF
ELSE
  IF (ilen > 0) THEN
    CALL grib_write_bytes (nudatc, ydata, ilen, ierrc)
    IF (ierrc /= GRIB_SUCCESS) THEN
      ierror = 1
      yerrmsg = 'Error in grib_write_bytes'
      RETURN
    ENDIF
  ENDIF
ENDIF

IF (npes > 1) THEN
  IF (ltime_barrier) THEN
    CALL MPI_Barrier(icomm,implcode)
    IF ( implcode /= 0 ) THEN
      yerrmsg = 'Error in MPI_BARRIER'
      ierror  = 11
      RETURN
    ENDIF
  ENDIF
ENDIF
#endif

END SUBROUTINE write_gribapi

!==============================================================================
!==============================================================================
!+ Module procedure in "io_utilities" for writing a netcdf file
!------------------------------------------------------------------------------

SUBROUTINE write_netcdf (nudat, data, i_tot, j_tot, irec_len, ncorg, icomm,  &
                         my_id, npes, idattyp, lasync, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   Processor 0 gathers all fields from the other processors and writes them
!   to netcdf files.
!
!------------------------------------------------------------------------------
!
! Subroutine arguments
! Scalar arguments with intent(in):
INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
  nudat,     & ! internal file descriptor
  i_tot,     & ! i-dimension for the array data
  j_tot,     & ! j-dimension for the array data
  irec_len,  & ! length of the record (0 if no record available)
  icomm,     & ! MPI communicator
  my_id,     & ! Processor ID in the communicator icomm
  npes,      & ! number of total PEs
  idattyp      ! data type for MPI interface

INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
  ncorg(3,0:npes-1) 
               ! organizational variables for netcdf. Only Processor 0 has 
               ! correct information, which is necessary for writing the data
               !   ncorg(1,:): actual number of level
               !   ncorg(2,:): maximum number of levels for this variable
               !   ncorg(3,:): netcdf ID of this variable

! Array arguments with intent(in):
REAL (KIND=irealgrib),     INTENT(IN)    ::  &
  data(i_tot*j_tot)   
               ! array to be written

LOGICAL                  , INTENT(IN)    ::  &
  lasync       ! indicates whether asynchronous IO is used or not
               ! (not used at the moment)

! Scalar arguments with intent(out):
CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
  yerrmsg      ! string for error messages

INTEGER  (KIND=iintegers), INTENT(OUT)   ::  &
  ierror       ! error status

!------------------------------------------------------------------------------

! Local scalars
INTEGER  (KIND=iintegers)                ::  &
  ipr,           & ! loop counter
  ilen,          & ! i_tot*j_tot
  implcode,      & ! Error variable for MPI
  izerr            ! another error variable

INTEGER  (KIND=iintegers)                ::  &
  istatus(MPI_STATUS_SIZE)                     ! for status in MPI calls

REAL(KIND=irealgrib)                     ::  &
  rbuf(i_tot*j_tot), & ! buffer for receiving the data
  reof                 ! for sending end-of-file message

!- End of header
!------------------------------------------------------------------------------

#ifdef NETCDF
ilen     = i_tot * j_tot
izerr    = 0_iintegers
ierror   = 0
reof     = -999999.0_irealgrib
yerrmsg  = '         '

IF(my_id == 0) THEN

  ! Processor 0 receives the data and writes it out
  DO ipr = 0, npes-1

    IF (ipr == 0) THEN

     IF (irec_len > 0) THEN
       ! Processor 0 has the data already and can write it to disk
       IF     (ncorg(2,ipr) == 1) THEN
         ! this is a 2D variable
         izerr = nf90_put_var (nudat, ncorg(3,ipr), data,              &
                    start=(/ 1, 1, 1 /), count=(/ i_tot, j_tot, 1 /))
         IF (izerr /= nf90_noerr) THEN
           yerrmsg = 'Error writing netcdf 2D variable'
           ierror  = 1
           RETURN
         ENDIF
       ELSEIF (ncorg(2,ipr) >  1) THEN
         ! this is a 2D-slice of a 3D variable
         izerr = nf90_put_var (nudat, ncorg(3,ipr), data,              &
            start=(/ 1, 1, ncorg(1,ipr),1 /), count=(/ i_tot, j_tot, 1, 1 /))
         IF (izerr /= nf90_noerr) THEN
           yerrmsg = 'Error writing netcdf 3D variable'
           ierror  = 2
           RETURN
         ENDIF
       ENDIF
     ENDIF

    ELSE

      ! Processor 0 gets the data from the other processors
      CALL MPI_RECV(rbuf, ilen, idattyp, ipr, 1002, icomm, istatus, implcode)
      IF ( implcode /= 0 ) THEN
        yerrmsg = 'Error in MPI_RECV record'
        ierror  = 3
        RETURN
      ENDIF

      IF (rbuf(1) /= reof) THEN
        ! processor ipr really sent a field that has to be written
        ! (otherwise rbuf(1) = reof: then just an end message is sent)
        IF     (ncorg(2,ipr) == 1) THEN
          ! this is a 2D variable
          izerr = nf90_put_var (nudat, ncorg(3,ipr), rbuf,             &
                     start=(/ 1, 1, 1 /), count=(/ i_tot, j_tot, 1 /))
          IF (izerr /= nf90_noerr) THEN
            yerrmsg = 'Error writing netcdf 2D variable'
            ierror  = 4
            RETURN
          ENDIF
        ELSEIF (ncorg(2,ipr) >  1) THEN
          ! this is a 2D-slice of a 3D variable
          izerr = nf90_put_var (nudat, ncorg(3,ipr), rbuf,             &
             start=(/ 1, 1, ncorg(1,ipr),1 /), count=(/ i_tot, j_tot, 1, 1 /))
          IF (izerr /= nf90_noerr) THEN
            yerrmsg = 'Error writing netcdf 3D variable'
            ierror  = 5
            RETURN
          ENDIF
        ENDIF
      ENDIF

    ENDIF

  ENDDO

ELSE

  ! All other Processors just send the data
  IF (irec_len > 0) THEN
    CALL MPI_SEND(data, ilen, idattyp, 0, 1002, icomm, implcode)
    IF ( implcode /= 0 ) THEN
      yerrmsg = 'Error in MPI_SEND record'
      ierror  = 6
      RETURN
    ENDIF
  ELSE
    CALL MPI_SEND(reof, 1, idattyp, 0, 1002, icomm, implcode)
    IF ( implcode /= 0 ) THEN
      yerrmsg = 'Error in MPI_SEND end message'
      ierror  = 7
      RETURN
    ENDIF
  ENDIF

ENDIF
#endif

END SUBROUTINE write_netcdf

!==============================================================================
!==============================================================================
!+ Module procedure in "io_utilities" for writing a binary restart file
!------------------------------------------------------------------------------

SUBROUTINE write_restart (nudat, data, i_tot, j_tot, irec_len, ipds_arr,     &
                          npds_dim, igds_arr, ngds_dim, icomm, my_id, npes,  &
                          idattyp, lasync, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   Processor 0 gathers all fields from the other processors and writes them
!   to a binary restart file
!
!------------------------------------------------------------------------------
!
! Subroutine arguments
! Scalar arguments with intent(in):
INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
  nudat,     & ! internal file descriptor
  i_tot,     & ! i-dimension for the array data
  j_tot,     & ! j-dimension for the array data
  irec_len,  & ! length of the record (0 if no record available)
  icomm,     & ! MPI communicator
  my_id,     & ! Processor ID in the communicator icomm
  npes,      & ! number of total PEs
  idattyp      ! data type for MPI interface

INTEGER  (KIND=intgribf),  INTENT(IN)    ::  &
  npds_dim,  & ! dimension for the product definition section
  ngds_dim     ! dimension for the grid definition section

INTEGER  (KIND=intgribf),  INTENT(IN)    ::  &
  ipds_arr(npds_dim), & ! product definition section
  igds_arr(ngds_dim)    ! grid definition section

! Array arguments with intent(in):
REAL (KIND=wp),            INTENT(IN)    ::  &
  data(i_tot*j_tot)   
               ! array to be written

LOGICAL                  , INTENT(IN)    ::  &
  lasync       ! indicates whether asynchronous IO is used or not
               ! (not used at the moment)

! Scalar arguments with intent(out):
CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
  yerrmsg      ! string for error messages

INTEGER  (KIND=iintegers), INTENT(OUT)   ::  &
  ierror       ! error status

!------------------------------------------------------------------------------

! Local scalars
INTEGER  (KIND=iintegers)                ::  &
  ipr,           & ! loop counter
  ilen,          & ! i_tot*j_tot
  implcode,      & ! Error variable for MPI
  izerr            ! another error variable

INTEGER  (KIND=iintegers)                ::  &
  istatus(MPI_STATUS_SIZE)                     ! for status in MPI calls

REAL(KIND=wp)                            ::  &
  rbuf(i_tot*j_tot),  & ! buffer for receiving the data
  reof                  ! to indicate end-of-file

INTEGER  (KIND=intgribf)                 ::  &
  ibuf_pds(npds_dim), & ! buffer for product definition section
  ibuf_gds(ngds_dim)    ! buffer for grid definition section

!- End of header
!------------------------------------------------------------------------------

ilen     = i_tot * j_tot
izerr    = 0_iintegers
ierror   = 0
reof     = HUGE (1.0_wp)
yerrmsg  = '         '

IF(my_id == 0) THEN

  ! Processor 0 receives the data and writes it out
  DO ipr = 0, npes-1

    IF (ipr == 0) THEN

      IF (irec_len > 0) THEN
        ! Processor 0 has the data already and can write it to disk
        WRITE (nudat,IOSTAT=izerr) ipds_arr, igds_arr, data

        IF (izerr /= 0_iintegers) THEN
          yerrmsg = 'Error writing restart variable'
          ierror  = 1
          RETURN
        ENDIF
      ENDIF

    ELSE

      ! Processor 0 gets the data from the other processors
      CALL MPI_RECV (ibuf_pds, npds_dim, MPI_INTEGER, ipr, 1003, icomm,    &
                     istatus, implcode)
      IF ( implcode /= 0 ) THEN
        yerrmsg = 'Error in MPI_RECV pds'
        ierror  = 2
        RETURN
      ENDIF

      CALL MPI_RECV (ibuf_gds, ngds_dim, MPI_INTEGER, ipr, 1004, icomm,    &
                     istatus, implcode)
      IF ( implcode /= 0 ) THEN
        yerrmsg = 'Error in MPI_RECV gds'
        ierror  = 3
        RETURN
      ENDIF

      CALL MPI_RECV (rbuf,    ilen,      idattyp,     ipr, 1005, icomm,    &
                     istatus, implcode)
      IF ( implcode /= 0 ) THEN
        yerrmsg = 'Error in MPI_RECV record'
        ierror  = 4
        RETURN
      ENDIF

      IF (rbuf(1) /= reof) THEN
        ! processor ipr really sent a field that has to be written
        WRITE (nudat, IOSTAT=izerr) ibuf_pds, ibuf_gds, rbuf

        IF (izerr /= 0_iintegers) THEN
          yerrmsg = 'Error writing restart variable'
          ierror  = 5
          RETURN
        ENDIF

      ENDIF

    ENDIF

  ENDDO

ELSE

  ! All other Processors just send the data
  CALL MPI_SEND (ipds_arr, npds_dim, MPI_INTEGER, 0, 1003, icomm, implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'Error in MPI_SEND pds'
    ierror  = 6
    RETURN
  ENDIF

  CALL MPI_SEND (igds_arr, ngds_dim, MPI_INTEGER, 0, 1004, icomm, implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'Error in MPI_SEND gds'
    ierror  = 7
    RETURN
  ENDIF

  IF (irec_len > 0) THEN
    CALL MPI_SEND (data, ilen, idattyp, 0, 1005, icomm, implcode)
    IF ( implcode /= 0 ) THEN
      yerrmsg = 'Error in MPI_SEND record'
      ierror  = 8
      RETURN
    ENDIF
  ELSE
    CALL MPI_SEND (reof, 1, idattyp, 0, 1005, icomm, implcode)
    IF ( implcode /= 0 ) THEN
      yerrmsg = 'Error in MPI_SEND end message'
      ierror  = 9
      RETURN
    ENDIF
  ENDIF

ENDIF


END SUBROUTINE write_restart

!==============================================================================
!==============================================================================
!+ checks the input and the output data
!------------------------------------------------------------------------------

SUBROUTINE check_record                                                   &
           (field, idim1s, idim1e, idim2s, idim2e, idim3s, idim3e,        &
                   icom1s, icom1e, icom2s, icom2e, icom3s, icom3e, undef, &
            yname, iee, ilevel, lwork, nout, npe, icomm, myid, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!  Check_record calculates the minimum, the maximum and the meanvalue of a 
!  two-dimensional total field and prints them to a file with unit number
!  nout. Three dimensions have to be given for the two-dimensional field,
!  because this routine is also used for GME-records, where the third dimension
!  represents the diamond. In case of LM, the parameters of the third 
!  dimensions are 1. 
!  
!  This routine has to be called by all compute PEs. Every processor writes 
!  its results to a line of characters. These lines are gathered by compute
!  PE 0, which prints them.
!
! Method:
!  Fortran 90 intrinsic functions are used to determine the minimum, maximum
!  and the meanvalue of ufeld. The data from all processors are gathered with
!  MPI_GATHER to processor 0. For the T3E implementation, the characters have
!  to be transformed to integers (at the moment).
!
! Output files:
!  Ascii file with unit number nout
!
!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):

INTEGER (KIND=iintegers), INTENT(IN) ::  &
  idim1s, idim1e, idim2s, idim2e, idim3s, idim3e, & 
            ! start- and end-indices for the dimensions of field
  icom1s, icom1e, icom2s, icom2e, icom3s, icom3e
            ! start- and end-indices for the computations

INTEGER (KIND=intgribf),  INTENT(IN) ::  &
  iee,             & ! grib number                        "
  ilevel             ! level of this field

INTEGER (KIND=iintegers), INTENT(IN) ::  &
  nout               ! unit number of the output file

LOGICAL,                  INTENT(IN) ::  &
  lwork              ! is .TRUE. if there is work to do

INTEGER (KIND=iintegers), INTENT(IN) ::  &
  npe,             & ! number of pes in communicator icomm
  icomm,           & ! communicator of compute PEs
  myid               ! ID of this PE in communicator icomm

CHARACTER (LEN= *),       INTENT(IN) ::  &
  yname              ! name of the field in variable table

! Array arguments with intent(in):
REAL   (KIND=irealgrib),  INTENT(IN) ::  &
  field (idim1s:idim1e, idim2s:idim2e, idim3s:idim3e),  &
  undef              ! for points that are not defined by the bitmap

! Variables for output
CHARACTER (LEN= *),      INTENT(OUT) ::  &
  yerrmsg            ! for error message

INTEGER (KIND=iintegers),INTENT(OUT) ::  &
  ierror             ! error status

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER (KIND=iintegers)   ::   &
  j1min, j1max, j2min, j2max, jdmin, jdmax, idiv
                                  ! location of min and max in the 2D field

INTEGER (KIND=iintegers)   ::   &
  i, n, j1, j2, jd, implcode      ! Additional variables

REAL    (KIND=irealgrib)   ::   &
  zmin, zmax, zsum, zdiv,       & ! values of min, max and meanvalue
  undef_eps                       ! slightly changed undef-value to cope with
                                  ! numerical rounding
! Local arrays:
CHARACTER (LEN=200)        ::   &
  yline,                        & ! line to be written for one processor
  ylines(npe)                     ! array of lines to collect all lines from 
                                  ! the other processors

INTEGER                    ::   & ! really the standard integers!!
  iline (200),                  & ! the same for the integer values
  ilines(200,npe)   
!
!- End of header
!------------------------------------------------------------------------------

  ierror  = 0
  yerrmsg = '   '
  undef_eps = undef + 100.0_irealgrib

  j1min   = 1
  j2min   = 1
  j1max   = 1
  j2max   = 1

  IF (lwork) THEN
    ! if this processor has a required record then
    ! determine maximum and minimum with their locations
    zsum = 0.0_irealgrib
    zmin =   ABS (undef)
    zmax = - ABS (undef)
    idiv = 0_iintegers

    DO jd = icom3s, icom3e
      DO j2 = icom2s, icom2e
        DO j1 = icom1s, icom1e
          IF (field(j1,j2,jd) > undef_eps) THEN
            idiv = idiv + 1
            zsum = zsum + field(j1,j2,jd)
            IF (field(j1,j2,jd) > zmax) THEN
              zmax  = field(j1,j2,jd)
              j1max = j1
              j2max = j2
              jdmax = jd
            ENDIF
            IF (field(j1,j2,jd) < zmin) THEN
              zmin  = field(j1,j2,jd)
              j1min = j1
              j2min = j2
              jdmin = jd
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    IF (idiv /= 0_iintegers) THEN
      zdiv = REAL (idiv, irealgrib)
      zsum = zsum / zdiv
    ELSE
      ! this means that the whole field is undefined
      zsum = undef
    ENDIF

    ! Adjust the indices, if idim* is not equal to icom*, so that the
    ! indices are related to icom* for the first and the second dimension
    IF (icom1s > idim1s) THEN
      j1max = j1max - (icom1s-idim1s)
      j1min = j1min - (icom1s-idim1s)
    ENDIF
    IF (icom2s > idim2s) THEN
      j2max = j2max - (icom2s-idim2s)
      j2min = j2min - (icom2s-idim2s)
    ENDIF
   
    ! print the values to yline
    IF (icom3e == 1) THEN
      IF (icom2e == 1) THEN
        ! this must be an ICON field
        IF (zsum == undef) THEN
          WRITE(yline,'(A,A,I4,I4,A)')  '  ', yname, iee, ilevel, &
               '                       WHOLE FIELD UNDEFINED    '
        ELSE
          WRITE(yline,'(A,A,I4,I4,F18.6,I10,A,F18.6,I10,A,F18.7)')           &
                      '  ', yname, iee, ilevel, zmin,        j1min, '     ', &
                      zmax,        j1max, '     ', zsum
        ENDIF
      ELSE
        ! likely it is an LM field
        IF (zsum == undef) THEN
          WRITE(yline,'(A,A,I4,I4,A)')  '  ', yname, iee, ilevel, &
               '                       WHOLE FIELD UNDEFINED    '
        ELSE
          WRITE(yline,'(A,A,I4,I4,F18.6,2I5,A,F18.6,2I5,A,F18.7)')           &
                      '  ', yname, iee, ilevel, zmin, j1min, j2min, '     ', &
                      zmax, j1max, j2max, '     ', zsum
        ENDIF
      ENDIF
    ELSE
      ! this really is an GME-field
      IF (zsum == undef) THEN
        WRITE(yline,'(A,A,I4,I4,A)')  '  ', yname, iee, ilevel, &
             '                       WHOLE FIELD UNDEFINED    '
      ELSE
        WRITE(yline,'(A,A,I4,I4,F18.6,3I5,F18.6,3I5,F18.7)')               &
                    '  ', yname, iee, ilevel, zmin, j1min, j2min, jdmin,   &
                    zmax, j1max, j2max, jdmax, zsum
      ENDIF
    ENDIF
  ELSE
    ! print just one blank to yline
    yline = ' '
  ENDIF
  
  ! collect all yline to ylines in processor 0
  IF (npe > 1) THEN
    ! for T3E's sake, convert the characters to standard integers
    DO i=1,200
       iline(i) = ICHAR(yline(i:i))
    ENDDO
    CALL MPI_Gather(iline, 200, MPI_INTEGER, ilines, 200, MPI_INTEGER,    &
                    0, icomm, implcode)
  
  ! CALL MPI_Gather(yline, 200, MPI_CHARACTER, ylines, 200, MPI_CHARACTER,&
  !                 0, icomm, implcode)
    IF ( implcode /= 0 ) THEN
      ierror  = implcode
      yerrmsg = 'Error in MPI_GATHER'
      RETURN
    ENDIF
  ELSE
    ylines(1) = yline
  ENDIF
  
  ! processor 0 prints all ylines
  IF( myid == 0) THEN
    IF(npe > 1) THEN
      DO n=1,npe
        DO i=1,200
          ylines(n)(i:i) = CHAR(ilines(i,n))
        ENDDO
      ENDDO
    ENDIF
    DO i=1,npe
      IF(ylines(i) /= ' ') THEN
        WRITE(nout,'(A)') ylines(i)(1:LEN_TRIM(ylines(i)))
      ENDIF
    ENDDO
  ENDIF
  
END SUBROUTINE check_record

!==============================================================================
!==============================================================================
!+ Check LM grid description section against the Namelist parameters
!------------------------------------------------------------------------------

SUBROUTINE check_input_grid (igrbed, igds, idimgds, ipds, idimpds, igrbhand,  &
                 ylib, yname, ydate, iedim, jedim, kedim, startlat, startlon, &
                 dlon, dlat, pollon, pollat, inv_in, pvz_in, vc_type,         &
                 lwork, npe, icomm, myid, lcheckdate, itype_calendar, ymodel, &
                 yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!  This routine checks the input grid against the Namelist parameters.
!  It is a collective routine and has to be called by all compute PEs.
!  For the input grid the following 10 values are checked:
!    ie_tot, je_tot, ke_tot, startlat, startlon, pollat, pollon, dlon, dlat
!    and the actual date of the forecast.
!
! Method:
!  If the values in the grid description section and the Namelist parameters
!  are not equal, an error status variable is set (/= 0). The error status
!  variables from all processors are gathered and scattered to all. If an
!  error occurs, the processor with lowest rank prints its (wrong) values.
!
!------------------------------------------------------------------------------

! Subroutine arguments
INTEGER (KIND=intgribf), INTENT(IN)   ::  &
  igrbed,           & ! grib edition number (from messages read)
  igrbhand,         & ! grib handler (only for grib_api)
  idimgds,          & ! dimension of the grid description section
  idimpds             ! dimension of the product definition section

INTEGER (KIND=intgribf), INTENT(IN)   ::  &
  igds(idimgds),    & ! grid description section
  ipds(idimpds)       ! product definition section

CHARACTER (LEN=*),        INTENT(IN)  ::  &
  ylib,             & ! used grib library (grb1 or apix)
  ymodel,           & ! name of the input model
  yname,            & ! short name of the product
  ydate               ! actual date from Namelist parameter

INTEGER (KIND=iintegers), INTENT(IN)  ::  &
  iedim, jedim, kedim, & ! dimensions of the fields
  inv_in                 ! number of vertical coordiante parameters

REAL    (KIND=wp),        INTENT(IN)  ::  &
  startlat, startlon, & ! coordinates of the lower left grid point
  dlon,  dlat,        & ! grid resolution in lambda and phi direction
  pollat, pollon        ! coordinates of the rotated north pole

REAL    (KIND=wp),        INTENT(IN)  ::  &
  pvz_in(inv_in)        ! vertical coordinate parameters already read in

TYPE(vcoord_type),        INTENT (IN) ::  &
  vc_type               ! vertical coordinate parameters

LOGICAL,                  INTENT(IN) ::   &
  lwork               ! is .TRUE. if there is work to do

INTEGER (KIND=iintegers), INTENT(IN)  ::  &
  itype_calendar,   & ! for specifying the calendar used
  npe,              & ! number of PEs in communicator icomm
  icomm,            & ! communicator for compute PEs
  myid                ! ID of this PE in communicator icomm

LOGICAL, INTENT(IN)       ::  &
  lcheckdate          ! indicates constant data for which no date is checked

! arguments with intent(out):
CHARACTER (LEN=*),        INTENT(OUT) ::  &
  yerrmsg             ! error message

INTEGER (KIND=iintegers), INTENT(OUT) ::  &
  ierror              ! error status variable

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER (KIND=iintegers)  ::  &
  i, izmplcode, ireturn, ic, ilenvt, ilasopo, ilosopo

REAL    (KIND=wp)         ::  &
  zstartlat_g, zstartlon_g,   & ! coordinates of the lower left grid point
  zendlat_g,   zendlon_g,     & ! coordinates of the lower left grid point
  zpollat_grib, zpollon_grib, & ! coordinates of the rotated north pole
  zdlon_grib, zdlat_grib,     & ! grid resolution in latitude and longitude
  zsouthp_lat, zsouthp_lon      ! south pole coordinates from NAMELIST input

REAL    (KIND=wp),        ALLOCATABLE     ::  zpv(:)

REAL    (KIND=irealgrib)  :: refstf, rlasopo, rlosopo, rla1, rlo1, rla2, rlo2

CHARACTER  (LEN=14)       ::  &
  ydate_ref     ! reference date from the grib

CHARACTER  (LEN=12)       ::  &
  ydate_val     ! reference date from the grib

! Local arrays:
INTEGER                   ::  & ! really the standard integer
  isenderror (0:npe-1),       & !
  irecverror (0:npe-1)          !

INTEGER (KIND=intgribf)   ::                                               &
  iyear_ref, imonth_ref, iday_ref, ihour_ref, imin_ref,                    &
  iyear, imonth, iday, ihour, imin, isec, iref_act, ierrf, ifactor,        &
  itimepassed, ie_grib, je_grib, ke_grib, ilevtyp, inv, iyear_val,         &
  imonth_val, iday_val, ihour_val, imin_val, ijincr, idatreptyp,           &
  idifftime, k, iznlev, iznrvg

CHARACTER (LEN=  1)       :: yzuuid(16)
CHARACTER (LEN= 30)       :: ylevtyp, ydatreptyp, yrefdate, yreftime,      &
                             yvaldate, yvaltime, yvaltime2
CHARACTER(LEN=uuid_string_length) :: yzuuidin_string, yzuuid_string

LOGICAL                   ::  &
  lzclimate     ! to indicate whether T_CL or W_CL are processed

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  ! initialize error status variable
  yerrmsg       = '   '
  ierror        = 0
  isenderror(:) = 0
  irecverror(:) = 0
  ierrf         = 0_intgribf

  IF (lwork) THEN

    IF     ((ylib == 'grb1') .OR. (ylib == 'bina')) THEN

      ! Check whether the grid is rotated_ll or regular_ll
      idatreptyp = igds(4)

      ! But these are now Grib1 codings!!!!
      IF     (idatreptyp == 10) THEN
        ! this is a rotated lat-lon grid
        ydatreptyp = 'rotated_ll'
      ELSEIF (idatreptyp ==  0) THEN
        ! this is a regular lat-lon grid
        ydatreptyp = 'regular_ll'
      ELSE
        ierror  = 21
        yerrmsg = 'Wrong type of grid in Grib'
        ydatreptyp = 'wrong_grid'
      ENDIF

#ifdef GRIBAPI
    ELSEIF (ylib == 'apix') THEN
      CALL grib_get (igrbhand, 'dataRepresentationType', idatreptyp,    ireturn)
      CALL grib_get (igrbhand, 'typeOfGrid',             ydatreptyp,      ireturn)
#endif
    ENDIF

!------------------------------------------------------------------------------
! Section 2.1: Compare the values for the grid dimensions
!------------------------------------------------------------------------------

    ! get the meta data for the grid:
    IF ((ylib == 'grb1') .OR. (ylib == 'bina')) THEN
      inv     = igds(2)
      ie_grib = igds(5)
      je_grib = igds(6)
      ke_grib =   -1
      ilevtyp = ipds(8)

      ylevtyp = '                '
      IF     (ilevtyp == 109) THEN
        ylevtyp = 'hybrid'
      ELSEIF (ilevtyp == 110) THEN
        ylevtyp = 'hybridLayer'
      ENDIF

#ifdef GRIBDWD
      ! the vertical coordinate parameters have already been read
      ! in case of restart files and must not be read again
      IF (((ilevtyp == 109) .OR. (ilevtyp == 110)) .AND. (ylib /= 'bina') ) THEN
        IF (ymodel == 'COSMO') THEN
          ! get vertical coordinate parameters
          ALLOCATE(zpv(inv), STAT=ierror)
          DO k = 1, inv
            ! NOTE: the below can not be replaced by ABS(...) since the smallest
            !       negative INTEGER may occurr which has no positive equivalent
            IF ( (igds(25+k) > 1000) .OR. (igds(25+k) < -1000) ) THEN
              zpv(k) = REAL (REFSTF(igds(25+k)), wp)
            ELSE
              zpv(k) = REAL (       igds(25+k) , wp)
            ENDIF
          ENDDO
        ENDIF
      ENDIF
#endif

    ELSEIF (ylib == 'apix') THEN

#ifdef GRIBAPI
      CALL grib_get (igrbhand, 'NV',                     inv    , ireturn)
      CALL grib_get (igrbhand, 'Ni',                     ie_grib, ireturn)
      CALL grib_get (igrbhand, 'Nj',                     je_grib, ireturn)
      CALL grib_get (igrbhand, 'indicatorOfTypeOfLevel', ilevtyp, ireturn)
      CALL grib_get (igrbhand, 'typeOfLevel'           , ylevtyp, ireturn)
      ke_grib =   -1

      IF ( (ymodel(1:2) == 'UM') .AND. (TRIM(yname) == 'V') ) THEN
        ! horizontal staggered V has one grid point less in west-east
        je_grib = je_grib+1
      ENDIF 

      IF ((ylevtyp == 'hybrid') .OR. (ylevtyp == 'hybridLayer') .AND. (inv > 0)) THEN
        IF (ymodel == 'COSMO') THEN
          ALLOCATE(zpv(inv), STAT=ierror)
          CALL grib_get (igrbhand, 'pv',                    zpv     , ireturn)
        ENDIF
      ELSEIF ((ylevtyp == 'generalVertical') .OR. (ylevtyp == 'generalVerticalLayer')) THEN
        CALL grib_get (igrbhand, 'nlev',                iznlev)
        CALL grib_get (igrbhand, 'numberOfVGridUsed',   iznrvg)
        CALL grib_get (igrbhand, 'uuidOfVGrid',         yzuuid)
      ENDIF
#endif
    ENDIF

    IF ( ((ylevtyp == 'hybrid') .OR. (ylevtyp == 'hybridLayer')) .AND. (ylib /= 'bina') ) THEN
      SELECT CASE (ymodel)
      CASE ('COSMO')
        ! Check that zpv are identical with already read pvz_in

        ! But it could be that initial and boundary data have different coding style
        IF     ( (zpv(1) < 500.0_wp) .AND. (pvz_in(1) < 500.0_wp) ) THEN
          ! both are of the same new coding style
          IF ( (inv_in /= inv) .AND. (inv_in+1 /= inv) ) THEN
            ke_grib = -1
            ierror = 23
          ELSE
            ! just check the vertical coordinate parameters
            DO k = 1, kedim+1
              IF (ABS(zpv(6+k)-pvz_in(6+k)) > 5.0E-3_wp) THEN
                ierror = 24
              ENDIF
            ENDDO
          ENDIF
        ELSEIF ( (zpv(1) > 500.0_wp) .AND. (pvz_in(1) > 500.0_wp) ) THEN
          ! both are of the same old coding style
          IF ( (inv_in /= inv) .AND. (inv_in+1 /= inv) ) THEN
            ke_grib = -1
            ierror = 23
          ELSE
            ! just check the vertical coordinate parameters
            DO k = 1, kedim+1
              IF (ABS(zpv(4+k)-pvz_in(4+k)) > 5.0E-3_wp) THEN
                ierror = 24
              ENDIF
            ENDDO
          ENDIF
        ELSEIF ( (zpv(1) < 500.0_wp) .AND. (pvz_in(1) > 500.0_wp) ) THEN
          ! just check the vertical coordinate parameters
          DO k = 1, kedim+1
            IF (ABS(zpv(6+k)-pvz_in(4+k)) > 5.0E-3_wp) THEN
              ierror = 24
            ENDIF
          ENDDO
        ELSEIF ( (zpv(1) > 500.0_wp) .AND. (pvz_in(1) < 500.0_wp) ) THEN
          ! just check the vertical coordinate parameters
          DO k = 1, kedim+1
            IF (ABS(zpv(4+k)-pvz_in(6+k)) > 5.0E-3_wp) THEN
              ierror = 24
            ENDIF
          ENDDO
        ENDIF

      CASE ('GME','IFS', 'UMR','UMG','GSM','GFS')
        IF (inv /= 0) THEN
          ke_grib = INT (0.5_wp * inv) - 1
        ELSE
          ke_grib = -1
        ENDIF
      CASE ('HIRLM')
        IF (inv /= 0) THEN
          ke_grib = inv / 2
        ELSE
          ke_grib = -1
        ENDIF
      END SELECT
    ELSE
      ke_grib = -1
    ENDIF

    ! check the field dimensions
    IF (iedim /= ie_grib) ierror = 1
    IF (jedim /= je_grib) ierror = 2

    IF ((ylevtyp == 'hybrid') .OR. (ylevtyp == 'hybridLayer')) THEN
      IF (ke_grib /= -1) THEN
        IF (kedim /= ke_grib ) ierror = 3
      ENDIF
    ELSEIF ( ((ylevtyp == 'generalVertical') .OR. (ylevtyp == 'generalVerticalLayer')) ) THEN
      IF (iznlev /= vc_type%nlevels) ierror = 3
      IF (iznrvg /= vc_type%ivctype) ierror = 3
      IF (.NOT. uuid_match(vc_type%vc_uuid, yzuuid)) ierror = 3
      CALL uuid_2char(yzuuid,               yzuuid_string)
      CALL uuid_2char(vc_type%vc_uuid,      yzuuidin_string)
    ENDIF

!------------------------------------------------------------------------------
! Section 2.2: Compare the values for the rotated pole
!------------------------------------------------------------------------------

    IF (ydatreptyp(1:10) == 'rotated_ll') THEN
      ! Instead of recomputing the given grib coordinates of the rotated south
      ! pole to the north pole, the namelist coordinates of the north pole are
      ! computed to the south pole in the same way as it is done in GME2LM.
      IF     ((ylib == 'grb1') .OR. (ylib == 'bina')) THEN
        zpollat_grib =   REAL(igds(20), wp)*0.001_wp
        zpollon_grib =   REAL(igds(21), wp)*0.001_wp

      ELSEIF (ylib == 'apix') THEN
#ifdef GRIBAPI
        IF (ymodel == 'UMR') THEN
          ! just for UM2LM, but why???
          CALL grib_get (igrbhand, 'latitudeOfSouthernPole',  ilasopo, ireturn)
          rlasopo = REAL(ilasopo, irealgrib) * 1E-6_irealgrib
          CALL grib_get (igrbhand, 'longitudeOfSouthernPole', ilosopo, ireturn)
          rlosopo = REAL(ilosopo, irealgrib) * 1E-6_irealgrib
        ELSE
          CALL grib_get (igrbhand, 'latitudeOfSouthernPoleInDegrees',  rlasopo, ireturn)
          CALL grib_get (igrbhand, 'longitudeOfSouthernPoleInDegrees', rlosopo, ireturn)
        ENDIF

        zpollat_grib =   REAL(rlasopo, wp)
        zpollon_grib =   REAL(rlosopo, wp)
        ! Make zpollon_grib conforming to COSMO-norm
        IF (zpollon_grib > 180.0_wp) THEN
          zpollon_grib = zpollon_grib - 360.0_wp
        ENDIF
#endif
      ENDIF

      zsouthp_lat  =   - pollat
      zsouthp_lon  =   pollon + 180.0_wp
      IF (zsouthp_lon > 180.0_wp) THEN
        zsouthp_lon = zsouthp_lon - 360.0_wp
      ENDIF

      ! change to 1E-4, because there really were rounding problems with 1E-5
      IF (ABS(zpollat_grib - zsouthp_lat) > 1.0E-4_wp) ierror = 6
      IF (ABS(zpollon_grib - zsouthp_lon) > 1.0E-4_wp) ierror = 7
    ENDIF

!------------------------------------------------------------------------------
! Section 2.3: Compare the values for the grid domain
!------------------------------------------------------------------------------

    ! Get the domain values from the meta data
    IF     ((ylib == 'grb1') .OR. (ylib == 'bina')) THEN
      SELECT CASE (ymodel)
      CASE ('IFS','GSM')
        ! models with different scanning mode
        zstartlat_g = REAL(igds(10), wp)*0.001_wp
        zendlat_g   = REAL(igds( 7), wp)*0.001_wp
      CASE DEFAULT
        zstartlat_g = REAL(igds( 7), wp)*0.001_wp
        zendlat_g   = REAL(igds(10), wp)*0.001_wp
      END SELECT
      zstartlon_g = REAL(igds( 8), wp)*0.001_wp
      zendlon_g   = REAL(igds(11), wp)*0.001_wp

      ! for the increments check the resolution and component flag
      IF (IBITS(igds(9),7,1) /= 0) THEN
        zdlat_grib   = REAL(igds(13), wp)*0.001_wp
        zdlon_grib   = REAL(igds(12), wp)*0.001_wp
      ELSE
        ! ABS takes care of the scanning mode ??
        ! zdlat_grib   = ABS(zendlat_g - zstartlat_g) / (jedim - 1)
        zdlat_grib   = (zendlat_g - zstartlat_g) / (jedim - 1)
        IF ((zendlon_g-zstartlon_g) < 0.0_wp) THEN
          ! If the area is located around the 180-Meridian, (zendlon_g-zstartlon_g)
          ! will be negative and 360 degrees have to be added
          zdlon_grib   = (zendlon_g-zstartlon_g + 360.0_wp) / (iedim-1)
        ELSE
          zdlon_grib   = (zendlon_g - zstartlon_g) / (iedim - 1)
        ENDIF
      ENDIF

    ELSEIF (ylib == 'apix') THEN
#ifdef GRIBAPI
      CALL grib_get (igrbhand, 'latitudeOfFirstGridPointInDegrees',  rla1,  ireturn)
      CALL grib_get (igrbhand, 'longitudeOfFirstGridPointInDegrees', rlo1,  ireturn)
      CALL grib_get (igrbhand, 'latitudeOfLastGridPointInDegrees',   rla2,  ireturn)
      CALL grib_get (igrbhand, 'longitudeOfLastGridPointInDegrees',  rlo2,  ireturn)

      SELECT CASE (ymodel)
      CASE ('IFS','GSM')
        ! models with different scanning mode
        zstartlat_g = REAL(rla2,wp)
        zendlat_g   = REAL(rla1,wp)
      CASE DEFAULT
        zstartlat_g = REAL(rla1,wp)
        zendlat_g   = REAL(rla2,wp)
      END SELECT

      ! Careful: for grib2, longitudes are between    0...360 degrees,
      !          for grib1, longitudes are between -180...180
      ! the COSMO convention still is -180...180
      IF     (igrbed == 1) THEN
        zstartlon_g = REAL(rlo1, wp)
        zendlon_g   = REAL(rlo2, wp)
      ELSEIF (igrbed == 2) THEN
        zstartlon_g = REAL(rlo1, wp)
        IF (rlo1 >= 180.0_wp) THEN
          zstartlon_g = REAL(rlo1, wp) - 360.0_wp
        ENDIF

        zendlon_g   = REAL(rlo2, wp)
        IF (rlo2 >  180.0_wp) THEN
          zendlon_g   = REAL(rlo2, wp) - 360.0_wp
        ENDIF
      ENDIF

      ! for the increments check the resolution and component flag
      CALL grib_get (igrbhand, 'ijDirectionIncrementGiven', ijincr,  ireturn)
      IF (ijincr == 1) THEN
        CALL grib_get (igrbhand, 'iDirectionIncrementInDegrees', zdlon_grib, ireturn)
        CALL grib_get (igrbhand, 'jDirectionIncrementInDegrees', zdlat_grib, ireturn)
      ELSE
        zdlat_grib   = (zendlat_g - zstartlat_g) / (jedim - 1)
        IF ((zendlon_g-zstartlon_g) < 0.0_wp) THEN
          ! If the area is located around the 180-Meridian, (zendlon_g-zstartlon_g)
          ! will be negative and 360 degrees have to be added
          zdlon_grib   = (zendlon_g-zstartlon_g + 360.0_wp) / (iedim-1)
        ELSE
          zdlon_grib   = (zendlon_g - zstartlon_g) / (iedim - 1)
        ENDIF
      ENDIF
#endif
    ENDIF

    ! Check the values against the Namelist variables and allow for a staggered grid
    ! (these usually are U and V from the COSMO-Model and the Unified-Model;
    !  check only for dlat/2, dlon/2 for all U and V variables)
    IF ((TRIM(yname) == 'U') .OR. (TRIM(yname) == 'AUMFL_S')) THEN
      IF (ABS(zstartlat_g  - startlat) > 1.0E-3_wp) ierror = 4
      IF (ABS(zstartlon_g  - startlon) > zdlon_grib + 1.0E-3_wp) ierror = 5
    ELSEIF ((TRIM(yname) == 'V').OR. (TRIM(yname) == 'AVMFL_S')) THEN
      IF (ABS(zstartlat_g  - startlat) > zdlat_grib + 1.0E-3_wp) ierror = 4
      IF (ABS(zstartlon_g  - startlon) > 1.0E-3_wp) ierror = 5
    ELSE
      IF (ABS(zstartlat_g  - startlat) > 1.0E-3_wp) ierror = 4
      IF (ABS(zstartlon_g  - startlon) > 1.0E-3_wp) ierror = 5
    ENDIF

    IF (ABS(zdlat_grib - dlat) > 1.0E-3_wp) ierror = 8
    IF (ABS(zdlon_grib - dlon) > 1.0E-3_wp) ierror = 9

!------------------------------------------------------------------------------
! Section 2.4: Check the date
!------------------------------------------------------------------------------

    IF (lcheckdate) THEN
      ! check the date
      ! Before, only the reference date was checked here.
      ! But what really is needed is to check the time for which the product
      ! is valid. This has to be calculated using ipds (16-19): unit of time
      ! range, time stamps and time range indicator

      ! Only for the fields T_CL and W_CL (climatological soil values) it is
      ! possible that older fields are taken (if LM initial data are computed
      ! from an older GME-file, e.g. from gfff00120000 from the run 12 hours
      ! before). T_CL and W_CL would then be taken from the corresponding
      ! giff00000000-file, and thus would be 12 hours older).

      ! Check whether T_CL or W_CL are processed
      IF (TRIM(yname) == 'T_CL' .OR. TRIM(yname) == 'W_CL') THEN
        lzclimate = .TRUE.
      ELSE
        lzclimate = .FALSE.
      ENDIF

      ! Determine the reference date from grib values
      IF ((ylib == 'grb1') .OR. (ylib == 'bina')) THEN

        IF (ipds(11) == 100) THEN
          WRITE(ydate_ref,'(7I2.2)')                          &
                           ipds(22)  ,       0 , ipds(12), ipds(13), ipds(14), &
                           ipds(15)  ,       0
          iyear_ref = ipds(22) * 100
        ELSE
          WRITE(ydate_ref,'(7I2.2)')                          &
                           ipds(22)-1, ipds(11), ipds(12), ipds(13), ipds(14), &
                           ipds(15)  ,       0
          iyear_ref = (ipds(22)-1) * 100 + ipds(11)
        ENDIF
        imonth_ref = ipds(12)
        iday_ref   = ipds(13)
        ihour_ref  = ipds(14)
        imin_ref   = ipds(15)

        ! Determine the corresponding integer values for the actual date "ydate"
        ! no minutes are given here
        READ(ydate    ,'(I4,5I2)') iyear, imonth, iday, ihour, imin, isec

        ! Determine the difference between reference date and actual date
        ! in minutes
        CALL diff_minutes (iyear_ref, imonth_ref, iday_ref, ihour_ref, imin_ref, &
                           iyear,     imonth,     iday,     ihour,     imin,     &
                           itype_calendar, iref_act,  ierrf)

        ! Determine the factor for calculating the time passed (since the
        ! reference date) in minutes. For that purpose the unit of time range
        ! (ipds(16)) is needed.
        SELECT CASE (ipds(16))
        CASE (0)
          ! unit is minute
          ifactor = 1
        CASE (1)
          ! unit is hour
          ifactor = 60
        CASE (2)
          ! unit is day
          ifactor = 1440
        CASE (10)
          ! unit is 3 hours
          ifactor = 180
        CASE (11)
          ! unit is 6 hours
          ifactor = 360
        CASE (12)
          ! unit is 12 hours
          ifactor = 720
        CASE (13)
          ! unit is 15 minutes (DWD)
          ifactor = 15
        CASE (14)
          ! unit is 30 minutes (DWD)
          ifactor = 30
        CASE DEFAULT
          ! this is not implemented or even not defined
          ierror  = 11
        END SELECT

        IF (ierror == 0) THEN
          ! Now determine the time passed (since the reference date) in minutes.
          ! For that the time range indicator (ipds(19)) is needed to select the
          ! correct grib octet (ipds(17) or ipds(18)).
          SELECT CASE (ipds(19))
          CASE (0,1,13)
            itimepassed = ifactor * ipds(17)
          CASE (10)
            itimepassed = ifactor * ipds(18)
          CASE (2,3,4,5)
            ! this is a product which is valid for a time period and cannot be
            ! used as input field, but for restart fields
            !US:  ierror = 12
            itimepassed = ifactor * ipds(18)
          CASE DEFAULT
            ! these time range indicators are not used in the LM-Package
            ierror = 13
          END SELECT
        ENDIF

        IF (ierror == 0) THEN
          ! The timepassed and the difference of the dates are sometimes equal,
          ! but with hincbound=0.5/0.25, there can be differences of a timestep
          ! So allow for a difference of 2 minutes
          idifftime = ABS(iref_act - itimepassed)
          IF (idifftime > 2_iintegers) THEN
            IF (.NOT. lzclimate) THEN
              ierror = 10
            ELSE
              PRINT *, 'Using older climate values (T_CL, W_CL) from ', ydate_ref
            ENDIF
          ENDIF
        ENDIF

      ELSEIF (ylib == 'apix') THEN
#ifdef GRIBAPI
        CALL grib_get (igrbhand, 'dataDate',     yrefdate, ireturn)    ! edition independent
        CALL grib_get (igrbhand, 'dataTime',     yreftime, ireturn)    ! edition independent
        CALL grib_get (igrbhand, 'validityDate', yvaldate, ireturn)    ! edition independent
        CALL grib_get (igrbhand, 'validityTime', yvaltime, ireturn)    ! edition independent
        !US adapt ydate_ref to 14 digits: dataTime should have 4 digits:  ddmm
        !                                 so set ss to 0
        ydate_ref(1:14) = yrefdate(1:8)//yreftime(1:4)//'00'

        ! Leading 0 are missing in yvaltime!! It should be a 4 character string
        yvaltime2 = '0000'
        ilenvt = LEN_TRIM(yvaltime)
        DO ic = 1, ilenvt
          yvaltime2(5-ic:5-ic) = yvaltime(ilenvt+1-ic:ilenvt+1-ic)
        ENDDO
        ydate_val(1:12) = yvaldate(1:8)//yvaltime2(1:4)


        ! Determine the integer values for the actual date "ydate"
        ! no minutes are given here
        READ(ydate    ,'(I4,5I2)') iyear, imonth, iday, ihour, imin, isec

        ! and for the validity date
        READ(ydate_val,'(I4,4I2)') iyear_val, imonth_val, iday_val, ihour_val, imin_val

        ! Determine the difference between these two dates in minutes
        ! because of a dt that does not fit into a full hour, there might be differences
        CALL diff_minutes (iyear_val, imonth_val, iday_val, ihour_val, imin_val, &
                           iyear,     imonth,     iday,     ihour,     imin,     &
                           itype_calendar, iref_act,  ierrf)

        ! The timepassed and the difference of the dates are sometimes equal,
        ! but with hincbound=0.5/0.25, there can be differences of a timestep
        ! So allow for a difference of 2 minutes
        idifftime = ABS(iref_act)
        IF (idifftime > 2_iintegers) THEN
          IF (.NOT. lzclimate) THEN
            ierror = 10
            PRINT *, 'Actual   date:  ', iyear, imonth, iday, ihour, imin
            PRINT *, 'Validity date:  ', iyear_val, imonth_val, iday_val, ihour_val, imin_val
            PRINT *, 'Difference   :  ', idifftime, '  minutes!'
          ELSE
            PRINT *, 'Using older climate values (T_CL, W_CL) from ', ydate_ref
          ENDIF
        ENDIF
#endif
      ENDIF  ! grb1 or apix

    ENDIF ! lcheckdate

  ENDIF   ! lwork

!------------------------------------------------------------------------------
! Section 3: Exchange error status variable
!------------------------------------------------------------------------------

  isenderror (myid) = ierror

  IF (npe > 1) THEN
    CALL MPI_ALLGATHER (isenderror(myid), 1, MPI_INTEGER,              &
                        irecverror,       1, MPI_INTEGER, icomm, izmplcode)
  ELSE
    irecverror(:) = isenderror(:)
  ENDIF

!------------------------------------------------------------------------------
! Section 4: If necessary, print errors
!------------------------------------------------------------------------------

  DO i = 0, npe-1
    IF (irecverror(i)/= 0) THEN
      ierror = irecverror(i)

      ! only one processor prints the error message
      IF (myid == i) THEN
        IF (ierror < 10) THEN
          PRINT *, '         data file          namelist input', ierror

          PRINT *, 'ie_tot                ',ie_grib     ,'       ',iedim
          PRINT *, 'je_tot                ',je_grib     ,'       ',jedim
        IF ((ylevtyp == 'hybrid') .OR. (ylevtyp == 'hybridLayer')) THEN
          IF (ke_grib /= -1) THEN
          PRINT *, 'ke_tot                ',ke_grib     ,'       ',kedim
          ENDIF
        ELSEIF ( ((ylevtyp == 'generalVertical') .OR. (ylevtyp == 'generalVerticalLayer')) ) THEN
          PRINT *, 'nlevels               ',vc_type%nlevels  ,'       ',iznlev
          PRINT *, 'Number of Grid Used   ',vc_type%ivctype  ,'       ',iznrvg
          PRINT *, 'UUID                  ',yzuuidin_string  ,'       ',yzuuid_string
        ENDIF

          PRINT *, 'startlat_tot          ',zstartlat_g ,'       ',startlat
          PRINT *, 'startlon_tot          ',zstartlon_g ,'       ',startlon

          IF (ydatreptyp(1:10) == 'rotated_ll') THEN
          PRINT *, 'rot. south pole (lat) ',zpollat_grib,'       ',zsouthp_lat
          PRINT *, 'rot. south pole (lon) ',zpollon_grib,'       ',zsouthp_lon
          ENDIF

          PRINT *, 'dlat                  ',zdlat_grib  ,'       ',dlat
          PRINT *, 'dlon                  ',zdlon_grib  ,'       ',dlon
        ELSEIF (ierror == 10) THEN
          PRINT *, ' The actual date           ', ydate
          PRINT *, ' and the reference date    ', ydate_ref, ' + ', &
                     idifftime, 'minutes    do not match!!'
        ELSEIF (ierror == 11) THEN
          PRINT *, ' Wrong (or not implemented) unit of time (uot)', ipds(16)
          PRINT *, '    Implemented values are uot =  0 (minute)     '
          PRINT *, '                           uot =  1 (hour)       '
          PRINT *, '                           uot =  2 (day)        '
          PRINT *, '                           uot = 10 ( 3 hours)   '
          PRINT *, '                           uot = 11 ( 6 hours)   '
          PRINT *, '                           uot = 12 (12 hours)   '
          PRINT *, '                           uot = 13 (15 minutes) '
          PRINT *, '                           uot = 14 (30 minutes) '
        ELSEIF (ierror == 12) THEN
          PRINT *, ' Wrong time range indicator (tri)', ipds(19)
          PRINT *, '    The values tri = 2,3,4,5 mean that the product is '
          PRINT *, '    valid for a  time period and cannot  be used as a '
          PRINT *, '    input field'
        ELSEIF (ierror == 13) THEN
          PRINT *, ' This time range indicator (tri)', ipds(19)
          PRINT *, ' is not used in the LM Package'
        ELSEIF (ierror == 21) THEN
          PRINT *, ' Wrong type of grid: ', ydatreptyp
        ELSEIF (ierror == 23) THEN
          PRINT *, ' Number of vertical coordinate parameters not consistent in:  ', &
                     myid, inv_in, inv
        ELSEIF (ierror == 24) THEN
          PRINT *, ' vertical coordinate parameters not consistent in:  ', myid
          PRINT *, ' pv from first data set:  ', pvz_in
          PRINT *, ' pv from later data set:  ', zpv
        ELSE
          PRINT *, ' ??? This is an unknown error ??? ', ierror, i, myid, npe
        ENDIF
      ENDIF
      EXIT
    ENDIF
  ENDDO

  IF (ALLOCATED(zpv)) DEALLOCATE(zpv)

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE check_input_grid

!==============================================================================
!==============================================================================
!+ Check GME grid description section against the Namelist parameters
!------------------------------------------------------------------------------

SUBROUTINE check_gme_grid (igds, idimgds, ipds, idimpds, igrbhand,           &
                 ylib, yname, ydate, i3edim, nidim, ni2dim, ni3dim, ndiam,   &
                 lwork, npe, icomm, myid, lcheckdate, itype_calendar,        &
                 yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!  This routine checks the GME-grid against the Namelist parameters. 
!  It is a collective routine and has to be called by all compute PEs. 
!  For GME the following 5 values are checked:
!    i3e_gme, ni_gme, ni2, ni3, nd
!
! Method:
!  If the values in the grid description section and the Namelist parameters
!  are not equal, an error status variable is set (/= 0). The error status
!  variables from all processors are gathered and scattered to all. If an 
!  error occurs, the processor with lowest rank prints its (wrong) values.
!
!------------------------------------------------------------------------------

! Subroutine arguments
INTEGER (KIND=intgribf), INTENT(IN)   ::  &
  igrbhand,         & ! grib handler (only for grib_api)
  idimgds,          & ! dimension of the grid description section
  idimpds             ! dimension of the product definition section

INTEGER (KIND=intgribf), INTENT(IN)   ::  &
  igds(idimgds),    & ! grid description section
  ipds(idimpds)       ! product definition section

CHARACTER (LEN=*),        INTENT(IN)  ::  &
  ylib,             & ! used grib library (grb1 or apix)
  yname,            & ! short name of the product
  ydate               ! actual date from Namelist parameter

INTEGER (KIND=iintegers), INTENT(IN)  ::  &
  i3edim, nidim, ni2dim, ni3dim, ndiam
 
LOGICAL,                  INTENT(IN) ::   &
  lwork               ! is .TRUE. if there is work to do

INTEGER (KIND=iintegers), INTENT(IN)  ::  &
  itype_calendar,   & ! for specifying the calendar used
  npe,              & ! number of PEs in communicator icomm
  icomm,            & ! communicator for compute PEs
  myid                ! ID of this PE in communicator icomm
 
LOGICAL, INTENT(IN)       ::  &
  lcheckdate          ! indicates constant data for which no date is checked

! arguments with intent(out):
CHARACTER (LEN=*),        INTENT(OUT) ::  &
  yerrmsg             ! error message

INTEGER (KIND=iintegers), INTENT(OUT) ::  &
  ierror              ! error status variable

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER (KIND=iintegers)  ::  &
  i3egrib, i, izmplcode, nigrib, ni2grib, ni3grib, ndiamgrib, inv, ireturn, &
  ic, ilenvt

LOGICAL                   ::  &
  lzclimate     ! to indicate whether T_CL or W_CL are processed

CHARACTER  (LEN=14)       ::  &
  ydate_ref     ! reference date from the grib

CHARACTER  (LEN=12)       ::  &
  ydate_val     ! reference date from the grib   (!does not contain seconds!!!)

CHARACTER (LEN= 30)       :: yrefdate, yvaldate, ylevtyp
CHARACTER (LEN= 10)       :: yreftime, yvaltime
CHARACTER (LEN=  4)       :: yvaltime2

! Local arrays: 
INTEGER                   ::  & ! really the standart integer
  isenderror (0:npe-1),       & !
  irecverror (0:npe-1)          !

INTEGER (KIND=intgribf)   ::                                               &
  iyear_ref, imonth_ref, iday_ref, ihour_ref, imin_ref,                    &
  iyear, imonth, iday, ihour, imin, isec, iref_act, ierrf, ifactor,        &
  itimepassed, iyear_val, imonth_val, iday_val, ihour_val, imin_val,       &
  iednr

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  ! initialize error status variable
  yerrmsg       = '   '
  ierror        = 0
  isenderror(:) = 0
  irecverror(:) = 0
  ierrf         = 0_intgribf

!------------------------------------------------------------------------------
! Section 2.1: Compare the values for the grid dimensions
!------------------------------------------------------------------------------

  IF (lwork) THEN
    ! check the field dimensions
    IF (ylib == 'grb1') THEN
      iednr     = 1
      IF     (ipds(8) == 109) THEN
        ylevtyp = 'hybrid'
      ELSEIF (ipds(8) == 110) THEN
        ylevtyp = 'hybridLayer'
      ELSE
        ylevtyp = 'others'
      ENDIF

      inv       = igds(2)
      nigrib    = igds(8)
      ni2grib   = igds(5)
      ni3grib   = igds(6)
      ndiamgrib = igds(7)

#ifdef GRIBAPI
    ELSEIF (ylib == 'apix') THEN
      CALL grib_get (igrbhand, 'editionNumber',     iednr, ireturn)
      CALL grib_get (igrbhand, 'typeOfLevel'  ,   ylevtyp, ireturn)

      CALL grib_get (igrbhand, 'Ni', nigrib   , ireturn)
      CALL grib_get (igrbhand, 'n2', ni2grib  , ireturn)
      CALL grib_get (igrbhand, 'n3', ni3grib  , ireturn)
      CALL grib_get (igrbhand, 'nd', ndiamgrib, ireturn)
      CALL grib_get (igrbhand, 'NV', inv      , ireturn)
#endif
    ENDIF

    ! check vertical levels and vertical coordinates only from
    ! levtyp = 'hybrid' or 'hybridlayer'

    IF ( (ylevtyp == 'hybrid') .OR. (ylevtyp == 'hybridLayer') ) THEN
      i3egrib = inv/2 - 1
      IF (i3egrib   /= i3edim ) ierror = 1
    ENDIF

    IF (nigrib    /= nidim  ) ierror = 2
    IF (ni2grib   /= ni2dim ) ierror = 3
    IF (ni3grib   /= ni3dim ) ierror = 4
    IF (ndiamgrib /= ndiam  ) ierror = 5

!------------------------------------------------------------------------------
! Section 2.2: Check the date
!------------------------------------------------------------------------------

    IF (lcheckdate) THEN

      ! check the date
      ! Before, only the reference date was checked here.
      ! But what really is needed is to check the time for which the product
      ! is valid. This has to be calculated using ipds (16-19): unit of time
      ! range, time stamps and time range indicator
      ! Only for the fields T_CL and W_CL (climatological soil values) it is
      ! possible that older fields are taken (if LM initial data are computed
      ! from an older GME-file, e.g. from gfff00120000 from the run 12 hours
      ! before). T_CL and W_CL would then be taken from the corresponding
      ! giff00000000-file, and thus would be 12 hours older).

      ! Check whether T_CL or W_CL are processed
      IF (TRIM(yname) == 'T_CL' .OR. TRIM(yname) == 'W_CL') THEN
        lzclimate = .TRUE.
      ELSE
        lzclimate = .FALSE.
      ENDIF

      ! Determine the reference date from grib values
      IF (ylib == 'grb1') THEN
#ifdef GRIBDWD
        IF (ipds(11) == 100) THEN
          WRITE(ydate_ref,'(7I2.2)')                          &
                           ipds(22)  ,       0 , ipds(12), ipds(13), ipds(14), &
                           ipds(15)  ,       0
          iyear_ref = ipds(22) * 100
        ELSE
          WRITE(ydate_ref,'(7I2.2)')                          &
                           ipds(22)-1, ipds(11), ipds(12), ipds(13), ipds(14), &
                           ipds(15)  ,       0
          iyear_ref = (ipds(22)-1) * 100 + ipds(11)
        ENDIF
        imonth_ref = ipds(12)
        iday_ref   = ipds(13)
        ihour_ref  = ipds(14)
        imin_ref   = ipds(15)

        ! Determine the corresponding integer values for the actual date "ydate"
        ! no minutes are given here
        READ(ydate    ,'(I4,5I2)') iyear, imonth, iday, ihour, imin, isec

        ! Determine the difference between reference date and actual date
        ! in minutes
        CALL diff_minutes (iyear_ref, imonth_ref, iday_ref, ihour_ref, imin_ref, &
                           iyear,     imonth,     iday,     ihour,     imin,     &
                           itype_calendar, iref_act,  ierrf)

        ! Determine the factor for calculating the time passed (since the
        ! reference date) in minutes. For that purpose the unit of time range
        ! (ipds(16)) is needed.
        SELECT CASE (ipds(16))
        CASE (0)
          ! unit is minute
          ifactor = 1
        CASE (1)
          ! unit is hour
          ifactor = 60
        CASE (2)
          ! unit is day
          ifactor = 1440
        CASE (10)
          ! unit is 3 hours
          ifactor = 180
        CASE (11)
          ! unit is 6 hours
          ifactor = 360
        CASE (12)
          ! unit is 12 hours
          ifactor = 720
        CASE (13)
          ! unit is 15 minutes (DWD)
          ifactor = 15
        CASE (14)
          ! unit is 30 minutes (DWD)
          ifactor = 30
        CASE DEFAULT
          ! this is not implemented or even not defined
          ierror  = 11
        END SELECT

        IF (ierror == 0) THEN
          ! Now determine the time passed (since the reference date) in minutes.
          ! For that the time range indicator (ipds(19)) is needed to select the
          ! correct grib octet (ipds(17) or ipds(18)).
          SELECT CASE (ipds(19))
          CASE (0,1,13)
            itimepassed = ifactor * ipds(17)
          CASE (10)
            itimepassed = ifactor * ipds(18)
          CASE (2,3,4,5)
            ! this is a product which is valid for a time period and cannot be
            ! used as input field
            ierror = 12
          CASE DEFAULT
            ! these time range indicators are not used in the LM-Package
            ierror = 13
          END SELECT
        ENDIF

        IF (ierror == 0) THEN
          ! The timepassed and the difference of the dates are sometimes equal,
          ! but with hincbound=0.5/0.25, there can be differences of a timestep
          ! So allow for a difference of 2 minutes
          IF ( ABS(iref_act - itimepassed) > 2_iintegers) THEN
            IF (.NOT. lzclimate) THEN
              ierror = 10
            ELSE
              PRINT *, 'Using older climate values (T_CL, W_CL) from ', ydate_ref
            ENDIF
          ENDIF
        ENDIF
#endif

      ELSEIF (ylib == 'apix') THEN
#ifdef GRIBAPI
        CALL grib_get (igrbhand, 'dataDate',     yrefdate, ireturn)    ! edition independent
        CALL grib_get (igrbhand, 'dataTime',     yreftime, ireturn)    ! edition independent
        CALL grib_get (igrbhand, 'validityDate', yvaldate, ireturn)    ! edition independent
        CALL grib_get (igrbhand, 'validityTime', yvaltime, ireturn)    ! edition independent
        !US adapt ydate_ref to 14 digits: dataTime should have 4 digits:  ddmm
        !                                 so set ss to 0
        ydate_ref(1:14) = yrefdate(1:8)//yreftime(1:4)//'00'

        ! Leading 0 are missing in yvaltime!! It should be a 4 character string
        yvaltime2 = '0000'
        ilenvt = LEN_TRIM(yvaltime)
        DO ic = 1, ilenvt
          yvaltime2(5-ic:5-ic) = yvaltime(ilenvt+1-ic:ilenvt+1-ic)
        ENDDO
        ydate_val(1:12) = yvaldate(1:8)//yvaltime2(1:4)


        ! Determine the integer values for the actual date "ydate"
        ! no minutes are given here
        READ(ydate    ,'(I4,5I2)') iyear, imonth, iday, ihour, imin, isec

        READ(ydate_val,'(I4,4I2)') iyear_val, imonth_val, iday_val, ihour_val, imin_val

        IF ( (ABS(iyear_val -iyear ) > 0)     .OR.    &
             (ABS(imonth_val-imonth) > 0)     .OR.    &
             (ABS(iday_val  -iday  ) > 0)     .OR.    &
             (ABS(ihour_val -ihour ) > 0)     .OR.    &
             (ABS(imin_val  -imin  ) > 0) ) THEN
           IF (.NOT. lzclimate) THEN
             ierror = 10
             PRINT *, 'Actual   date:  ', iyear, imonth, iday, ihour, imin
             PRINT *, 'Validity date:  ', iyear_val, imonth_val, iday_val, ihour_val, imin_val
           ELSE
             PRINT *, 'Using older climate values (T_CL, W_CL) from ', ydate_ref
           ENDIF
        ENDIF
#endif
      ENDIF  ! grb1 or apix

    ENDIF ! lcheckdate

  ENDIF   ! lwork

!------------------------------------------------------------------------------
! Section 3: Exchange error status variable
!------------------------------------------------------------------------------

  isenderror (myid) = ierror  

  IF (npe > 1) THEN
    CALL MPI_ALLGATHER (isenderror(myid), 1, MPI_INTEGER,              &
                        irecverror,       1, MPI_INTEGER, icomm, izmplcode)
  ELSE 
    irecverror(:) = isenderror(:)
  ENDIF

!------------------------------------------------------------------------------
! Section 4: If necessary, print errors
!------------------------------------------------------------------------------

  DO i = 0, npe-1
    IF (irecverror(i)/= 0) THEN
      ierror  = irecverror(i)
      yerrmsg = 'Wrong grid specification or date for GME data'

      ! only one processor prints the error message
      IF (myid == i) THEN
        IF (ierror < 6) THEN
          PRINT *, '         data file          namelist input'

          PRINT *, 'ni_gme     ',nigrib   ,'       ',nidim
          PRINT *, 'ni2        ',ni2grib  ,'       ',ni2dim
          PRINT *, 'ni3        ',ni3grib  ,'       ',ni3dim
          PRINT *, 'nd         ',ndiamgrib,'       ',ndiam
          PRINT *, 'i3e_gme    ',i3egrib  ,'       ',i3edim
        ELSEIF (ierror == 10) THEN
          PRINT *, ' The actual date           ', ydate
          PRINT *, ' and the reference date    ', ydate_ref, ' + ', &
                     itimepassed, 'minutes    do not match!!'
        ELSEIF (ierror == 11) THEN
          PRINT *, ' Wrong (or not implemented) unit of time (uot)', ipds(16)
          PRINT *, '    Implemented values are uot =  0 (minute)     '
          PRINT *, '                           uot =  1 (hour)       '
          PRINT *, '                           uot =  2 (day)        '
          PRINT *, '                           uot = 10 ( 3 hours)   '
          PRINT *, '                           uot = 11 ( 6 hours)   '
          PRINT *, '                           uot = 12 (12 hours)   '
          PRINT *, '                           uot = 13 (15 minutes) '
          PRINT *, '                           uot = 14 (30 minutes) '
        ELSEIF (ierror == 12) THEN
          PRINT *, ' Wrong time range indicator (tri)', ipds(19)
          PRINT *, '    The values tri = 2,3,4,5 mean that the product is '
          PRINT *, '    valid for a  time period and cannot  be used as a '
          PRINT *, '    input field'
        ELSEIF (ierror == 13) THEN
          PRINT *, ' This time range indicator (tri)', ipds(19)
          PRINT *, ' is not used in the LM Package'
        ELSE
          PRINT *, ' ??? This is an unknown error ??? '
        ENDIF
      ENDIF
      EXIT
    ENDIF
  ENDDO

!------------------------------------------------------------------------------
!- End of the Subroutine 
!------------------------------------------------------------------------------

END SUBROUTINE check_gme_grid

!==============================================================================
!==============================================================================
!+ Creates a "ready file", i.e. a small file indicating some event.
!------------------------------------------------------------------------------

SUBROUTINE write_ready_final(nstep, dt, yzhead_ready, lasync, npes, my_comm_id, &
                             ytrans_out, ntrans_out, ydate, ydate_ini,          &
                             yextension, itype_cal, lmmss, izdebug, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   Ready files are used within DWD's NWP suite to handle inter-dependencies
!   of programs.
!   Subroutine "write_ready_final" is called by all compute PEs.
!   It sets a flag, s.t. the last IO process creates a ready file.
!   The subroutine is called at the end of an output stage.
!
!------------------------------------------------------------------------------

  INTEGER (KIND=iintegers), INTENT(IN) :: nstep      ! actual time-step
  REAL  (KIND=wp)         , INTENT(IN) :: dt         ! long time-step
  CHARACTER (LEN=3)       , INTENT(IN) :: yzhead_ready ! head of the ready file
  LOGICAL                 , INTENT(IN) ::  &
    lasync       ! indicates whether asynchronous IO is used or not
  INTEGER (KIND=iintegers), INTENT(IN) ::  &
    ntrans_out,                            &         ! Unit Number for writing ready-Files during output
    npes,                                  &         ! number of compute PE
    my_comm_id,                            &         ! ID in communicator
    itype_cal                                        ! calendar used

  CHARACTER (LEN=250)     , INTENT(IN) ::  &
    ytrans_out      ! directory for writing ready-files

  CHARACTER (LEN=14),       INTENT(IN)  ::    &
    ydate,       & ! Actual date of the forecast in the form
                   !   yyyymmddhhmmss (yyyy: year, mm: month, dd: day, hh: time)
    ydate_ini      ! Initial date

  CHARACTER (LEN=1)     , INTENT(IN) ::  yextension

  LOGICAL, INTENT(IN) :: lmmss  ! .true.=new 14-digit date format, .false.=old 10-digit format

  INTEGER (KIND=iintegers), INTENT(IN) ::  &
    izdebug
  INTEGER                 , INTENT(OUT):: ierror     ! error code

  ! Local variables:
  CHARACTER (LEN=260)        ::   &
    yzname             ! full path- and file-name of ready files

  INTEGER (KIND=iintegers)   ::   &
    izerrstat,       & ! error status variable
    ifn_len

  ! nothing to do?
  IF (LEN_TRIM(ytrans_out) == 0) RETURN

  IF (izdebug > 10) THEN
    WRITE (*,*) "proc ", my_comm_id, ": Enter write_ready_final for file name head: ", yzhead_ready
  END IF

  ! Create the filename LMF_forecasttime
  CALL make_fn (yzhead_ready, ydate, ydate_ini, 'f', yextension, nstep, dt,     &
                .TRUE., itype_cal, ytrans_out, yzname, lmmss, izdebug, izerrstat)

  ifn_len = LEN_TRIM(yzname) ! file name length

  ! Write ready-file, if required
  ISPARALLEL: IF (lasync .OR. (npes > 1) ) THEN
    ! Open file using the parallel IO interface
    CALL mpe_io_ready(yzname(1:ifn_len), ifn_len, ntrans_out, ierror)
  ELSE
    ! Write the file
    OPEN  (ntrans_out, FILE=yzname, FORM='FORMATTED')
    WRITE (ntrans_out, '(A)') 'ready'
    CLOSE (ntrans_out)
  ENDIF ISPARALLEL

END SUBROUTINE write_ready_final

!==============================================================================
!==============================================================================
!+ calculate length of integer buffer for grib fields
!------------------------------------------------------------------------------

FUNCTION compute_grib_intbuffer_length(ie_tot, je_tot, nbytes, iwlength) RESULT(lfd)

  IMPLICIT NONE

  ! Parameters

  INTEGER, INTENT(IN) :: ie_tot, je_tot  ! horizontal grid dimensions
  INTEGER, INTENT(IN) :: nbytes          ! bytes per value for grib packing
  INTEGER, INTENT(IN) :: iwlength        ! length of integers used in the griblib in byte

  ! Local variables:

  INTEGER :: lfd

  lfd = ie_tot * je_tot * nbytes / iwlength + 5000   ! the "5000" just is safety

END FUNCTION compute_grib_intbuffer_length

!==============================================================================

END MODULE io_utilities
