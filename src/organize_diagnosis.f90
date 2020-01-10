!+ External procedure for organizing the diagnostic computations
!------------------------------------------------------------------------------

SUBROUTINE organize_diagnosis (yaction, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure is the driving routine for the diagnostic 
!   computations. It is called from the main program first to initialize the
!   diagnostic packages and then during the time stepping
!   to determine which diagnostic computations have to be performed.
!   The diagnostic computations are
!     - diagnosis:    mean values and maxima of certain subdomains
!     - differences:  differences of  certain values between the LM-forecast
!                     and the given boundary fields
!     - meanvalues:   mean values and maxima of the total domain
!     - grid points:  variables at selected grid points are stored for output
!
! Method:
!   Determine whether a certain action has to be performed in this time step
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.29       1999/05/11 Ulrich Schaettler
!  Initial release
! 1.34       1999/12/10 Ulrich Schaettler
!  Included Call for printing values and call to src_diagbudget
! 1.39       2000/05/03 Ulrich Schaettler
!  Included subroutine for Namelist input and some technical changes.
! 2.12       2001/11/07 Ulrich Schaettler
!  Changed my_cart_id to my_world_id for file treatment of INPUT_DIA
! 3.5        2003/09/02 Ulrich Schaettler
!  Eliminated diagnostics "differences" and "diagnosis"
! 3.6        2003/12/11 Ulrich Schaettler
!  Adaptations for multi-layer soil model;
!  Modifications for checking the IOSTAT-value when reading the NAMELISTs
! 3.7        2004/02/18 Ulrich Schaettler
!  Changed treatment of unit-numbers for ASCII-file handling
! 3.15       2005/03/03 Ulrich Schaettler
!  Changes to open and close file YUPRMEAN in every output step
! 3.18       2006/03/03 Ulrich Schaettler
!  Changes to split the meteograph output in a single file for every gridpoint
!    (file names contain either the name of a nearby station or geographical
!     coordinates in 0.1 degrees).
!  'lmulti_layer .AND. lgpspec' is possible now.
! 3.19       2006/04/25 Ulrich Schaettler
!  Corrected a bug in computing grid point indices for given lat-lon values
! 3.21       2006/12/04 Burkhardt Rockel, Ulrich Schaettler
!  Introduced polgam
!  Call Subroutine gridpoints with a timelevel (nnow) in the argument list
!  Splitted output for PRMEAN to PRMASS and PRHUMI
! V3_23        2007/03/30 Michael Baldauf, Ulrich Schaettler
!  New Namelist variables for calculation of integrals
!  Bug fix when computing nearest grid point for meteograph stations
!  Changes in calling grid point calculation to allow a flexible dt
! V3_24        2007/04/26 Ulrich Schaettler
!  Removed identical printouts
! V3_26        2007/05/29 Ulrich Schaettler
!  Bug correction when using nincgp (hincgp must be 0.0 then)
! V4_3         2008/02/25 Ulrich Schaettler
!  Correction for computation of nstepsgp (must not be done with integer arithmetic)
! V4_4         2008/07/16 Ulrich Schaettler
!  Read 2 additional NL Parameters itype_diag_t2m, itype_diag_gusts
! V4_8         2009/02/16 Ulrich Schaettler
!  Added options for itype_diag_t2m and itype_diag_gusts
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_21        2011/12/06 Jan-Peter Schulz
!  Added option 4 for itype_diag_gusts
! V4_23        2012/05/10 Oliver Fuhrer
!  Additional ASCII Output necessary for testsuite purposes
!  Add a new NL switch lyuprdbg to /DIACTL/ to activate this output
! V4_25        2012/09/28 Ulrich Schaettler
!  Increased length of station names to 30 characters
! V5_1         2014-11-28 Ulrich Blahak, Oliver Fuhrer
!  Changed the format of some YUSPECIF entries for the CLM namelist tool.
!  Implemented F2003 IOMSG-mechanism for better namelist error messages. (UB)
!  Replaced ireals by wp (working precision) (OF)
! V5_2         2015-05-21 Ulrich Schaettler
!  Include ntke in argument list for mean_testsuite, to select proper TKE time level
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

USE data_parameters,  ONLY:  wp, iintegers

USE data_runcontrol,  ONLY:                                                  &
       hstart, nstart, nstop, ntstep, nuspecif, nnow, lmulti_layer,          &
       idbg_level, ldebug_dia, lprintdeb_all, l2tls, itype_diag_t2m,         &
       itype_diag_gusts, ntke

USE data_modelconfig, ONLY:  ie_tot, je_tot, ke_tot, dt, startlat_tot,       &
                             startlon_tot, dlat, dlon, pollat, pollon,       &
                             polgam

USE data_parallel,    ONLY:                                                  &
       my_cart_id, num_compute, imp_integers, icomm_cart, my_world_id,       &
       icomm_world, imp_logical, nproc, intbuf, logbuf, realbuf, charbuf,    &
       imp_character, imp_reals, nboundlines

USE data_io,          ONLY: ydate_ini

USE parallel_utilities, ONLY:  distribute_values
USE environment,        ONLY:  get_free_unit
USE utilities,          ONLY:  rla2rlarot, phi2phirot

USE src_gridpoints,     ONLY:                                                &
       lgplong, lgpshort, lgpspec,nmaxgp, ngp, ngp_tot, n0gp, nincgp,        &
       h0gp, hincgp, nnextgp, hnextgp, nstepsgp,                             &
       stationlist_all, stationlist_tot, stationlist, list_of_stations
USE src_gridpoints,     ONLY:  init_gridpoints,  gridpoints,                 &
                               dealloc_gridpoints, nustgrpt, numdgrpt

USE src_meanvalues,     ONLY:  n0meanval, nincmeanval, nextmeanval,          &
                               nuprmass, nuprhumi, ltestsuite, nuprtest
USE src_meanvalues,     ONLY:  init_meanvalues,  meanvalues, mean_testsuite

USE src_diagbudget,     ONLY:  organize_diagbudget

USE src_integrals,      ONLY:  imin_integ, imax_integ, jmin_integ,           &
                               jmax_integ, kmin_integ, kmax_integ,           &
                               l_integrals
USE src_integrals,      ONLY:  organize_integrals

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Local organizational variables
  LOGICAL, SAVE                    ::           &
    lcomgp         ! perform grid point calculations

!==============================================================================

! Parameterlist
CHARACTER (LEN= *),       INTENT(IN)      ::  yaction
CHARACTER (LEN= *),       INTENT(OUT)     ::  yerrmsg
INTEGER (KIND=iintegers), INTENT(OUT)     ::  ierror

! Local scalars:
INTEGER (KIND=iintegers)                  ::                        &
  izerrstat, nuin, izdebug

REAL (KIND=wp)                            ::  zdt

CHARACTER (LEN= 9)                        ::  yinput

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------
!- Begin Subroutine organize_diagnosis
!------------------------------------------------------------------------------

ierror  = 0_iintegers
yerrmsg = '   '

! Initialize, whether debug output shall be done
IF (ldebug_dia) THEN
  IF (lprintdeb_all) THEN
    izdebug = idbg_level
  ELSE
    IF (my_cart_id == 0) THEN
      izdebug = idbg_level
    ELSE
      izdebug = 0
    ENDIF
  ENDIF
ELSE
  izdebug = 0
ENDIF

!------------------------------------------------------------------------------
! Section 1: Input of the Namelists
!------------------------------------------------------------------------------

  IF (yaction == 'input') THEN

    ! Open NAMELIST-INPUT file
    IF (my_world_id == 0) THEN
      IF (idbg_level > 0) THEN
        PRINT *,'    INPUT OF THE NAMELISTS FOR DIAGNOSTICS'
      ENDIF
      yinput   = 'INPUT_DIA'
      nuin     =  1
      OPEN(nuin   , FILE=yinput  , FORM=  'FORMATTED', STATUS='UNKNOWN',  &
           IOSTAT=izerrstat)
      IF(izerrstat /= 0) THEN
        yerrmsg  = ' ERROR    *** Error while opening file INPUT_DIA *** '
        ierror   = 2
        RETURN
      ENDIF
    ENDIF

    ! read in the NAMELIST-group
    CALL input_diactl (nuspecif, nuin, izerrstat)

    ! Close INPUT file
    IF (my_world_id == 0) THEN
      ! Close file for input of the NAMELISTS
      CLOSE (nuin    , STATUS='KEEP')
    ENDIF

    IF (izerrstat < 0) THEN
      yerrmsg = 'ERROR *** while reading NAMELIST Group /DIACTL/ ***'
      ierror  = 3
    ELSEIF (izerrstat > 0) THEN
      yerrmsg = 'ERROR *** Wrong values occured in NAMELIST INPUT_DIA ***'
      ierror  = 4
    ENDIF

!------------------------------------------------------------------------------
! Section 2: Initializations
!------------------------------------------------------------------------------

  ELSEIF (yaction == 'init') THEN

    IF (my_cart_id == 0) THEN
      PRINT *, '    INITIALIZE CONTROL VARIABLES AND MEAN VALUES'
    ENDIF
    IF (lcomgp) THEN
      ! get unit numbers for Ascii grid point output
      CALL get_free_unit (nustgrpt)
      CALL get_free_unit (numdgrpt)

      CALL init_gridpoints (ydate_ini)
    ENDIF

    ! get unit numbers for mean value output
    CALL get_free_unit (nuprmass)
    CALL get_free_unit (nuprhumi)
    IF (ltestsuite) CALL get_free_unit (nuprtest)
    CALL init_meanvalues (ydate_ini)

    IF ( l_integrals ) THEN
      CALL organize_integrals('init')
    END IF

!------------------------------------------------------------------------------
! Section 3: Compute the values during time stepping
!------------------------------------------------------------------------------

  ELSEIF (yaction == 'compute') THEN

  !----------------------------------------------------------------------------
  ! Section 3.1: mean values and maxima of the total domain
  !----------------------------------------------------------------------------
  
    ! check whether it has to be performed in this step
    IF ( (ntstep == nextmeanval) .OR. (ntstep == nstop) ) THEN
      IF (izdebug > 10) THEN
        PRINT *, 'calling mean values calculation in step ', ntstep
      ENDIF

      CALL meanvalues

      nextmeanval = nextmeanval + nincmeanval
    ENDIF
  
    ! check whether testsuite output has to be performed in this step
    ! note: Leapfrog should not output nnow at last timestep, as it will still
    !       be modified by time-filter and will not be reproducible as compared
    !       to a longer run
    IF (ltestsuite) THEN
      IF ( ( ANY(ntstep == (/ 0, 1, 2, 3 /)) .OR. (MOD(ntstep,10) == 0) ) &
           .AND. (l2tls .OR. ntstep /= nstop) ) THEN
        IF (izdebug > 10) THEN
          PRINT *, 'calling testsuite values calculation in step ', ntstep
        ENDIF
        CALL mean_testsuite(nnow, ntke)
      ENDIF
    ENDIF

  !----------------------------------------------------------------------------
  ! Section 3.2: grid point calculations
  !----------------------------------------------------------------------------
  
    IF (lcomgp) THEN
      IF ( (ntstep == nnextgp) .OR. (ntstep == nstop) ) THEN

        IF (izdebug > 10) THEN
          PRINT *, 'calling grid point calculation in step ', ntstep
        ENDIF

        CALL gridpoints (ydate_ini, nnow)

        ! determine next time step for grid point output
        IF (hincgp > 0.0_wp) THEN
          hnextgp = hnextgp + hincgp
          IF (ntstep == 0 .AND. (.NOT. l2tls)) THEN
            zdt = 2 * dt
          ELSE
            zdt =     dt
          ENDIF
          nnextgp = NINT (3600.0_wp * hnextgp / zdt)
        ELSE
          nnextgp = nnextgp + nincgp
        ENDIF
      ENDIF
    ENDIF

  !----------------------------------------------------------------------------
  ! Section 3.3: mean values and maxima of the total domain
  !----------------------------------------------------------------------------

    IF ( l_integrals ) THEN
      CALL organize_integrals('compute')
    END IF

!------------------------------------------------------------------------------
! Section 4: Calculation of budget diagnostics
!------------------------------------------------------------------------------

  ELSEIF (yaction == 'diagbudget') THEN

    ! This routine has to be called before the boundary relaxation in order
    ! to use the same values for u and v as in initialize_loop to compute
    ! the density

    CALL organize_diagbudget

!------------------------------------------------------------------------------
! Section 5: Deallocation of special fields for grid point output
!------------------------------------------------------------------------------

  ELSEIF (yaction == 'dealloc') THEN

    IF (lcomgp .AND. (ngp > 0) ) THEN
      CALL dealloc_gridpoints  (izerrstat)
    ENDIF

    IF ( l_integrals ) THEN
      CALL organize_integrals('dealloc')
    END IF

  ELSE

    ierror  = 1
    yerrmsg = 'ERROR *** No valid action for the diagnostics ***'

  ENDIF

!------------------------------------------------------------------------------
! Internal procedures
!------------------------------------------------------------------------------

CONTAINS

!==============================================================================
!+ Internal procedure in "organize_diagnosis" for the input of NAMELIST diactl
!------------------------------------------------------------------------------

SUBROUTINE input_diactl (nuspecif, nuin, ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine organizes the input of the NAMELIST-group diactl. 
!   The group diactl contains variables for controlling the diagnostic
!   computations:
!     - grid point calculations
!     - diagnosis:   area and volume means of selected areas
!     - differences: differences between the LM-forecast and the boundary
!                    values (full fields without spectral decomposition)
!     - mean values: calculates control variables
!   Note that the NAMELIST-groups lmgrid, runctl and dynctl have to be read 
!   before diactl!
!
! Method:
!   All variables are initialized with default values and then read in from
!   the file INPUT. The input values are checked for errors and for
!   consistency. If wrong input values are detected the program prints
!   an error message. The program is not stopped in this routine but an
!   error code is returned to the calling routine that aborts the program after
!   reading in all other namelists.
!   In parallel mode, the variables are distributed to all nodes with the
!   environment-routine distribute_values.
!   Both, default and input values are written to the file YUSPECIF
!   (specification of the run).
!
!------------------------------------------------------------------------------

! Parameter list
  INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
    nuspecif,     & ! Unit number for protocolling the task
    nuin            ! Unit number for Namelist INPUT file

  INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
    ierrstat        ! error status variable

!------------------------------------------------------------------------------

! Local variables

! Variables for default values

! grid point calculations
  REAL (KIND=wp)             ::       &
    h0gp_d,       & ! hour of the first grid point calculation (default)
    hincgp_d        ! hour increment of grid point calculations (default)

  REAL (KIND=wp)             ::       &
    rhuge,        & ! for default settings
    zlons, zlats, & ! 
    zlatgp_rot, zlongp_rot

  REAL (KIND=wp)             ::       &
    zrlon(ie_tot,je_tot), & ! for determination of igp_tot/jgp_tot
    zrlat(ie_tot,je_tot)    ! out of rlongp_tot/rlatgp_tot

  INTEGER   (KIND=iintegers)       ::           &
    n0gp_d,           & ! time step of the first grid point calculation
    nincgp_d,         & ! time step increment of grid point calculations
    itype_diag_t2m_d, & ! type of T_2M diagnostics
    itype_diag_gusts_d  ! type of gusts diagnostics

  TYPE (list_of_stations)    ::       &
    stationlist_tot_d (nmaxgp)

  LOGICAL                    ::       &
    lgpshort_d,   & ! calculate and print a short form of the
                    ! grid point calculations   (1 line/step)
    lgplong_d,    & ! calculate and print a long form of the
                    ! grid point calculations   (1 page/step)
    lgpspec_d       ! calculate and print a special form of the
                    ! grid point output used by the 1-d-model

  LOGICAL                    ::       &
    lzfound, lzignore

  CHARACTER (LEN=1) :: ychar

! calculation of integrals
  LOGICAL :: l_integrals_d

! coordinates of the corners of the cuboid for integration:
  INTEGER   (KIND=iintegers)       ::           &
    imin_integ_d, imax_integ_d, &
    jmin_integ_d, jmax_integ_d, &
    kmin_integ_d, kmax_integ_d

! mean values
  INTEGER   (KIND=iintegers)       ::           &
    n0meanval_d,  & ! time step of the first output
    nincmeanval_d   ! time step increment of outputs

! print out of values for the testsuite
  LOGICAL                          ::           &
    ltestsuite_d

! Other variables
  INTEGER (KIND=iintegers)   ::       &
    n, noffset, ierr, iz_err, i, j, ngp_tot_d

  CHARACTER(LEN=250)         :: iomsg_str

! Define the namelist group
  NAMELIST /diactl/ h0gp, hincgp, n0gp, nincgp, lgpshort, lgplong, lgpspec, &
                    n0meanval, nincmeanval, stationlist_tot,                &
                    imin_integ, imax_integ, jmin_integ, jmax_integ,         &
                    kmin_integ, kmax_integ, l_integrals, itype_diag_t2m,    &
                    itype_diag_gusts, ltestsuite

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE input_diactl
!------------------------------------------------------------------------------

ierrstat = 0_iintegers
iz_err   = 0_iintegers
rhuge    = HUGE (1.0_wp)

IF (my_world_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 1: Initialize the default variables
!------------------------------------------------------------------------------

! grid point calculations
  h0gp_d     = 0.0_wp
  hincgp_d   = 0.0_wp
  n0gp_d     = 0
  nincgp_d   = HUGE (1_iintegers)
  lgpshort_d = .FALSE.
  lgplong_d  = .FALSE.
  lgpspec_d  = .FALSE.

  ngp_tot_d  = 5
  ngp_tot    = ngp_tot_d
  stationlist_tot_d(1)%igp    =  5
  stationlist_tot_d(1)%jgp    =  5
  stationlist_tot_d(2)%igp    = ie_tot - 2
  stationlist_tot_d(2)%jgp    =  5
  stationlist_tot_d(3)%igp    = ie_tot / 2
  stationlist_tot_d(3)%jgp    = je_tot / 2
  stationlist_tot_d(4)%igp    =  5
  stationlist_tot_d(4)%jgp    = je_tot - 5
  stationlist_tot_d(5)%igp    = ie_tot - 5
  stationlist_tot_d(5)%jgp    = je_tot - 5

  stationlist_tot_d(:)%rlongp = rhuge
  stationlist_tot_d(:)%rlatgp = rhuge

  stationlist_tot_d(:)%ystation_name = '--------------------'

! mean values
  n0meanval_d   = 0
  nincmeanval_d = 10

! calculation of integrals
  ! default: integration area is the whole LM domain
  ! (important: NAMELISTs  'LMGRID' and 'RUNCTL' have to be known here!)
  l_integrals_d = .FALSE.
  imin_integ_d =      1 + nboundlines
  imax_integ_d = ie_tot - nboundlines
  jmin_integ_d =      1 + nboundlines
  jmax_integ_d = je_tot - nboundlines
  kmin_integ_d =      1
  kmax_integ_d = ke_tot

! different types of diagnostics
  itype_diag_t2m_d   = 1
  itype_diag_gusts_d = 1

! testsuite
  ltestsuite_d       = .FALSE.

!------------------------------------------------------------------------------
!- Section 2: Initialize default variables and error status variable
!------------------------------------------------------------------------------

! grid point calculations
  h0gp     = h0gp_d
  hincgp   = hincgp_d
  n0gp     = n0gp_d
  nincgp   = nincgp_d
  lgpshort = lgpshort_d
  lgplong  = lgplong_d
  lgpspec  = lgpspec_d

  stationlist_tot(:)%igp    = 0
  stationlist_tot(:)%jgp    = 0
  stationlist_tot(:)%rlongp = rhuge
  stationlist_tot(:)%rlatgp = rhuge
  stationlist_tot(:)%ystation_name = '--------------------'

! mean values
  n0meanval   = n0meanval_d
  nincmeanval = nincmeanval_d

! calculation of integrals
  l_integrals = l_integrals_d
  imin_integ  = imin_integ_d
  imax_integ  = imax_integ_d
  jmin_integ  = jmin_integ_d
  jmax_integ  = jmax_integ_d
  kmin_integ  = kmin_integ_d
  kmax_integ  = kmax_integ_d

! different types of diagnostics
  itype_diag_t2m     = itype_diag_t2m_d
  itype_diag_gusts   = itype_diag_gusts_d

! testsuite
  ltestsuite = ltestsuite_d

!------------------------------------------------------------------------------
!- Section 3: Input of the namelist values
!------------------------------------------------------------------------------

  iomsg_str(:) = ' '
  READ (nuin, diactl, IOSTAT=iz_err, IOMSG=iomsg_str)

  IF (iz_err /= 0) WRITE (*,'(A,A)') 'Namelist-ERROR DIACTL: ', TRIM(iomsg_str)
ENDIF

IF (nproc > 1) THEN
  ! distribute error status to all processors
  CALL distribute_values  (iz_err, 1, 0, imp_integers,  icomm_world, ierr)
ENDIF

IF (iz_err /= 0) THEN
  ierrstat = -1
  RETURN
ENDIF

IF (my_world_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 4: Check values for errors and consistency
!             Check grid points and domains for diagnosis
!------------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  !- Section 4.1: Grid point calculations
  !----------------------------------------------------------------------------

  ! Only one kind of grid point output (long or short) is permitted;
  ! if both are true, only lgpshort will be done
  IF( (lgpshort .EQV. .TRUE.) .AND. (lgplong .EQV. .TRUE.) ) THEN
    PRINT *,' WARNING  *** Only the short grid point output will be done ***'
    lgplong = .FALSE.
  ENDIF

  ! Check the input for the grid points:
  ! if i- and j-values have been read (/=0), then these values are taken
  ! if these values have not been read, the rlatgp-,rlongp values are
  ! checked 

  ! First compute rotated latitude and longitude values of all grid points 
  ! in the total domain
  DO j = 1, je_tot
    zlats = startlat_tot + (j-1) * dlat

    DO i = 1, ie_tot
      zlons  = startlon_tot + (i-1) * dlon
      IF (zlons  > 180.0_wp) THEN
        zlons  = zlons  - 360.0_wp
      ENDIF

      zrlat(i,j) = zlats
      zrlon(i,j) = zlons
    ENDDO
  ENDDO

  ngp_tot = 0

  gridpoint_loop: DO n = 1, nmaxgp
    IF ((stationlist_tot(n)%igp /= 0) .AND. (stationlist_tot(n)%jgp /= 0)) THEN
      ! check for correctness
      lzignore = .FALSE.
      IF ( stationlist_tot(n)%igp > ie_tot-1  ) THEN
        IF (izdebug > 1) THEN
          PRINT *,' WARNING *** stationlist_tot(',n,')%igp = ',        &
                 stationlist_tot(n)%igp, ' > ie_tot-1   = ', ie_tot-1 ,' *** '
        ENDIF
        lzignore = .TRUE.
      ENDIF
      IF ( stationlist_tot(n)%igp < 2 ) THEN
        IF (izdebug > 1) THEN
          PRINT *,' WARNING *** stationlist_tot(',n,')%igp = ',        &
                 stationlist_tot(n)%igp, ' < 2 *** '
        ENDIF
        lzignore = .TRUE.
      ENDIF
      IF ( stationlist_tot(n)%jgp > je_tot-1  ) THEN
        IF (izdebug > 1) THEN
          PRINT *,' WARNING *** stationlist_tot(',n,')%jgp = ',        &
                 stationlist_tot(n)%jgp, ' > je_tot-1   = ', je_tot-1 ,' *** '
        ENDIF
        lzignore = .TRUE.
      ENDIF
      IF ( stationlist_tot(n)%jgp < 2 ) THEN
        IF (izdebug > 1) THEN
          PRINT *,' WARNING *** stationlist_tot(',n,')%jgp = ',        &
                 stationlist_tot(n)%jgp, ' < 2 *** '
        ENDIF
        lzignore = .TRUE.
      ENDIF
      IF (lzignore) THEN
        IF (izdebug > 1) THEN
          PRINT *, ' WARNING *** grid point ', n, &
               stationlist_tot(n)%ystation_name, ' is ignored ***'
        ENDIF
        CYCLE gridpoint_loop
      ENDIF

      ! determination of  rlatgp- and rlongp-values for this grid point
      ! is done in the subdomain after the decomposition; but reset the
      ! values before the decomposition
      stationlist_tot(n)%rlatgp = rhuge
      stationlist_tot(n)%rlongp = rhuge

    ! at least one of the grid point indices was 0, so the rlatgp and
    ! rlongp values are checked
    ELSEIF ((stationlist_tot(n)%rlatgp /= rhuge) .AND.               &
            (stationlist_tot(n)%rlongp /= rhuge))    THEN
      ! try to find the next grid point and get its indices
      lzfound = .FALSE.
      DO j = 1, je_tot-1
        DO i = 1, ie_tot-1
          ! compute the rotated coordinates for the grid points
          zlatgp_rot = phi2phirot (stationlist_tot(n)%rlatgp,     &
                                   stationlist_tot(n)%rlongp, pollat, pollon)
          zlongp_rot = rla2rlarot (stationlist_tot(n)%rlatgp,     &
                                   stationlist_tot(n)%rlongp, pollat, pollon, polgam)

          IF (  (zrlon(i  ,j  ) <= zlongp_rot   ) .AND.                  &
                                  (zlongp_rot    < zrlon(i+1,j  )) .AND. &
                (zrlon(i  ,j+1) <= zlongp_rot   ) .AND.                  &
                                  (zlongp_rot    < zrlon(i+1,j+1)) .AND. &
                (zrlat(i  ,j  ) <= zlatgp_rot   ) .AND.                  &
                                  (zlatgp_rot    < zrlat(i  ,j+1)) .AND. &
                (zrlat(i+1,j  ) <= zlatgp_rot   ) .AND.                  &
                                  (zlatgp_rot    < zrlat(i+1,j+1)) ) THEN
            ! now we are in the correct grid cell and have to search the
            ! nearest grid point
            lzfound = .TRUE.

            IF ( (zlongp_rot    - zrlon(i,j)) <                          &
                              (zrlon(i+1,j) - zlongp_rot   ) ) THEN
              stationlist_tot(n)%igp = i
            ELSE
              stationlist_tot(n)%igp = i+1
            ENDIF

            IF ( (zlatgp_rot    - zrlat(i,j)) <                          &
                              (zrlat(i,j+1) - zlatgp_rot   ) ) THEN
              stationlist_tot(n)%jgp = j
            ELSE
              stationlist_tot(n)%jgp = j+1
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      IF (.NOT. lzfound) THEN
        IF (izdebug > 1) THEN
          WRITE (*,'(A,I4,A,2F8.3,A)') ' WARNING *** Grid point ',n, ': (',  &
                   stationlist_tot(n)%rlatgp, stationlist_tot(n)%rlongp,    &
                   ') is outside the total domain  ***'
          WRITE (*,'(A,I4,A,A)') ' WARNING *** grid point ', n, &
               stationlist_tot(n)%ystation_name, ' is ignored ***'
        ENDIF
        CYCLE gridpoint_loop
      ENDIF

    ELSEIF ((stationlist_tot(n)%igp    == 0)     .AND.       &
            (stationlist_tot(n)%jgp    == 0)     .AND.       &
            (stationlist_tot(n)%rlatgp == rhuge) .AND.       &
            (stationlist_tot(n)%rlongp == rhuge))    THEN
      ! no more grid point specification has been given
      EXIT gridpoint_loop
    ELSE
      ! otherwise, no correct specifications have been given
      IF (izdebug > 1) THEN
        PRINT *, ' WARNING *** No correct specifications are given',   &
                      ' for grid point ', n, ' ***'
        PRINT *, ' WARNING *** grid point ', n, &
               stationlist_tot(n)%ystation_name, ' is ignored ***'
      ENDIF
      CYCLE gridpoint_loop
    ENDIF

    ! If this point is reached, this grid point is taken
    ! store all settings in stationlist_all
    ngp_tot = ngp_tot + 1
    stationlist_all(ngp_tot)%igp           = stationlist_tot(n)%igp
    stationlist_all(ngp_tot)%jgp           = stationlist_tot(n)%jgp
    stationlist_all(ngp_tot)%rlongp        = stationlist_tot(n)%rlongp
    stationlist_all(ngp_tot)%rlatgp        = stationlist_tot(n)%rlatgp
    stationlist_all(ngp_tot)%ystation_name = stationlist_tot(n)%ystation_name
  ENDDO gridpoint_loop

  IF (ngp_tot == 0) THEN
    ! no specification of grid points has been given. 
    ! Take the default values
    stationlist_tot(1:5) = stationlist_tot_d(1:5)
    ngp_tot = 5
  ELSE
    ! replace all settings in stationlist_tot
    DO n = 1, ngp_tot
      stationlist_tot(n)%igp           = stationlist_all(n)%igp
      stationlist_tot(n)%jgp           = stationlist_all(n)%jgp
      stationlist_tot(n)%rlongp        = stationlist_all(n)%rlongp
      stationlist_tot(n)%rlatgp        = stationlist_all(n)%rlatgp
      stationlist_tot(n)%ystation_name = stationlist_all(n)%ystation_name
    ENDDO
  ENDIF

  ! Check the stationnames: there must be no blanks in it
  DO n = 1, ngp_tot
    ychar = stationlist_tot(n)%ystation_name(1:1)
    IF ( (ychar == ' ') .OR. (ychar == '-') .OR. (ychar == '')) THEN
      ! no station name has been given
      stationlist_tot(n)%ystation_name = '--------------------'
    ELSE
      ! check for blanks in station names and convert them to '_'
      DO j = 1, LEN_TRIM(stationlist_tot(n)%ystation_name)
        IF (stationlist_tot(n)%ystation_name(j:j) == ' ') THEN
          stationlist_tot(n)%ystation_name(j:j) = '_'
        ENDIF
      ENDDO
    ENDIF
  ENDDO

  ! Check whether the values for start and end of grid point calculations are
  ! given in hours and calculate the values in time steps

  IF ( hincgp /= hincgp_d ) THEN
    ! Check that only full hours, 0.5 or 0.25 are used
    IF ( ABS(REAL(NINT(hincgp), wp) - hincgp) > 1.0E-5_wp) THEN
      ! then it is not a full hour, only allow 0.25 and 0.5
      IF ( (hincgp /= 0.50_wp) .AND. (hincgp /= 0.25_wp) ) THEN
        PRINT *, 'ERROR: *** This is not a valid hincgp: ', hincgp, ' ***'
        PRINT *, '       *** only values = n.0 / 0.5 / 0.25 are allowed   ***'
        ierrstat = 1002
      ENDIF
    ENDIF
  ENDIF

  IF ( h0gp /= h0gp_d ) THEN
    ! Check that only full hours, 0.5 or 0.25 are used
    IF ( ABS(REAL(NINT(h0gp), wp) - h0gp) > 1.0E-5_wp) THEN
      ! then it is not a full hour, only allow 0.25 and 0.5
      IF ( (h0gp /= 0.50_wp) .AND. (h0gp /= 0.25_wp) ) THEN
        PRINT *, 'ERROR: *** This is not a valid h0gp: ', h0gp, ' ***'
        PRINT *, '       *** only values = n.0 / 0.5 / 0.25 are allowed   ***'
        ierrstat = 1002
      ENDIF
    ENDIF
  ENDIF
    
  ! Set hnextgp, nnextgp and adapt the values in case of restarts
  IF (h0gp >= hstart) THEN
    hnextgp = h0gp
  ELSE
    hnextgp = hstart   ! could be the case for restarts
  ENDIF
  nnextgp = NINT ( hnextgp   * 3600.0_wp/dt)

  ! Check, whether nincgp is set by namelist: then nnextgp has to be 
  ! re-evaluated
  IF ( nincgp /= nincgp_d ) THEN
    ! adapt n0gp in case of restarts
    IF (n0gp >= nstart) THEN
      nnextgp = n0gp
    ELSE
      ! determine the next step for grid point output with nincgp
      nnextgp = n0gp
      findn0: DO WHILE (nnextgp < nstart)
        nnextgp = nnextgp + nincgp
      ENDDO findn0
    ENDIF
  ELSE
    ! nincgp was not set, take hincgp as default
    nincgp  = NINT (hincgp * 3600.0_wp / dt)

    ! and also re-evaluate n0gp corresponding to value of h0gp
    n0gp = NINT ( h0gp   * 3600.0_wp/dt)
  ENDIF

  IF (n0gp < 0) THEN
    PRINT *,' WARNING  *** n0gp is set to 0 *** '
    n0gp = 0
  ENDIF

  IF (nincgp < 1) THEN
    PRINT *,' WARNING  *** nincgp is set to 1 hour *** '
    nincgp = NINT ( 3600.0_wp/dt)
  ENDIF
 
  IF (lgpspec) THEN
    nstepsgp = NINT( REAL((nstop - n0gp),wp) / REAL(nincgp,wp) ) + 1
    IF ( nstepsgp*ngp_tot > 5000 ) THEN
      PRINT *,' WARNING  *** Much memory is used for grid point output *** '
    ENDIF
  ELSEIF (lgplong .OR. lgpshort) THEN
    nstepsgp = 1    ! every computed step is written to file at once
  ELSE
    nstepsgp = 0
  ENDIF

  IF ( nstop < n0gp) THEN
    PRINT *,' WARNING  ***  n0gp > nstop:  no output of grid points  *** '
    nstepsgp = 0
    lgpshort = .FALSE.
    lgplong  = .FALSE.
    lgpspec  = .FALSE.
  ENDIF

  !----------------------------------------------------------------------------
  !- Section 4.2: Mean values
  !----------------------------------------------------------------------------

  IF (n0meanval < 0) THEN
    PRINT *,' WARNING  *** n0meanval is set to 0 *** '
    n0meanval = 0
  ENDIF

  IF (nincmeanval < 1) THEN
    PRINT *,' WARNING  *** nincmeanval is set to 10 steps *** '
    nincmeanval = 10
  ENDIF

  IF (n0meanval >= nstart) THEN
    nextmeanval = n0meanval
  ELSE
    nextmeanval = n0meanval
    findnext: DO WHILE (nextmeanval < nstart)
      nextmeanval = nextmeanval + nincmeanval
    ENDDO findnext
  ENDIF

  !----------------------------------------------------------------------------
  !- Section 4.3: different types of diagnostics
  !----------------------------------------------------------------------------

  IF ( (itype_diag_t2m < 1) .OR. (itype_diag_t2m > 2) ) THEN
    PRINT *, 'ERROR: *** Not a valid value for itype_diag_t2m: ', itype_diag_t2m
    PRINT *, '       *** only values 1 <= itype_diag_t2m <= 2 are allowed ***'
    ierrstat = 1002
  ENDIF

  IF ( (itype_diag_gusts < 1) .OR. (itype_diag_gusts > 4) ) THEN
    PRINT *, 'ERROR: *** Not a valid value for itype_diag_gusts: ', itype_diag_gusts
    PRINT *, '       *** only values 1 <= itype_diag_gusts <= 4 are allowed ***'
    ierrstat = 1002
  ENDIF
ENDIF

!------------------------------------------------------------------------------
!- Section 5: Distribute variables to all nodes
!------------------------------------------------------------------------------

IF (nproc > 1) THEN

  IF (my_world_id == 0) THEN
    intbuf ( 1) = n0meanval
    intbuf ( 2) = nincmeanval
    intbuf ( 3) = nextmeanval
    intbuf ( 4) = n0gp
    intbuf ( 5) = nincgp
    intbuf ( 6) = nnextgp
    intbuf ( 7) = nstepsgp
    intbuf ( 8) = ngp_tot
    intbuf ( 9) = imin_integ
    intbuf (10) = imax_integ
    intbuf (11) = jmin_integ
    intbuf (12) = jmax_integ
    intbuf (13) = kmin_integ
    intbuf (14) = kmax_integ
    intbuf (15) = itype_diag_t2m
    intbuf (16) = itype_diag_gusts

    realbuf( 1) = h0gp
    realbuf( 2) = hincgp
    realbuf( 3) = hnextgp

    ! indices of selected grid points
    DO n = 1,ngp_tot
      intbuf (20 + n)           = stationlist_tot(n)%igp
      intbuf (20 + ngp_tot + n) = stationlist_tot(n)%jgp
      realbuf(20 + n)           = stationlist_tot(n)%rlatgp
      realbuf(20 + ngp_tot + n) = stationlist_tot(n)%rlongp
      charbuf(n)(1:30)          = stationlist_tot(n)%ystation_name(1:30)
    ENDDO

    logbuf ( 1) = lgpshort
    logbuf ( 2) = lgplong
    logbuf ( 3) = lgpspec
    logbuf ( 4) = l_integrals
    logbuf ( 5) = ltestsuite
  ENDIF
  noffset = 20 + 2 * nmaxgp

  CALL distribute_values (intbuf,  noffset,  0, imp_integers, icomm_world, ierr)
  CALL distribute_values (realbuf, noffset,  0, imp_reals,    icomm_world, ierr)
  CALL distribute_values (charbuf, nmaxgp,   0, imp_character,icomm_world, ierr)
  CALL distribute_values (logbuf,  5,        0, imp_logical,  icomm_world, ierr)

  IF (my_world_id /= 0) THEN
    n0meanval     = intbuf ( 1)
    nincmeanval   = intbuf ( 2)
    nextmeanval   = intbuf ( 3)
    n0gp          = intbuf ( 4)
    nincgp        = intbuf ( 5)
    nnextgp       = intbuf ( 6)
    nstepsgp      = intbuf ( 7)
    ngp_tot       = intbuf ( 8)
    imin_integ    = intbuf ( 9)
    imax_integ    = intbuf (10)
    jmin_integ    = intbuf (11)
    jmax_integ    = intbuf (12)
    kmin_integ    = intbuf (13)
    kmax_integ    = intbuf (14)
    itype_diag_t2m= intbuf (15)
    itype_diag_gusts= intbuf (16)

    h0gp          = realbuf( 1)
    hincgp        = realbuf( 2)
    hnextgp       = realbuf( 3)

    ! indices of selected grid points
    DO n = 1,ngp_tot
      stationlist_tot(n)%igp = intbuf (20 + n)
      stationlist_tot(n)%jgp = intbuf (20 + ngp_tot + n)
      stationlist_tot(n)%rlatgp = realbuf (n)
      stationlist_tot(n)%rlongp = realbuf (ngp_tot + n)
      stationlist_tot(n)%ystation_name(1:30) = charbuf (n)(1:30)
    ENDDO

    lgpshort      = logbuf ( 1)
    lgplong       = logbuf ( 2)
    lgpspec       = logbuf ( 3)
    l_integrals   = logbuf ( 4)
    ltestsuite    = logbuf ( 5)
  ENDIF

ENDIF

! Set some additional control variables
lcomgp = lgpshort .OR. lgplong .OR. lgpspec

!------------------------------------------------------------------------------
!- Section 6: Output of the namelist variables and their default values
!------------------------------------------------------------------------------

IF (my_world_id == 0) THEN

  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(A23)') '0     NAMELIST:  diactl'
  WRITE (nuspecif, '(A23)') '      -----------------'
  WRITE (nuspecif, '(A2)')  '  '

  ! grid point calculations
  WRITE (nuspecif, '(A30)') '      Grid point calculations:'
  WRITE (nuspecif, '(T7,A,T33,A,T52,A,T70,A)') 'Variable', 'Actual Value',   &
                                               'Default Value', 'Format'

  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                     'lgpshort ',lgpshort , lgpshort_d,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                     'lgplong  ',lgplong  , lgplong_d ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                     'lgpspec  ',lgpspec  , lgpspec_d ,' L '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                                'n0gp  ',n0gp  ,n0gp_d,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                          'nincgp  ',nincgp  ,nincgp_d,' I '
  WRITE (nuspecif, '(T8,A,T33,F12.2,T52,F12.2,T71,A3)')                      &
                                                'h0gp  ',h0gp  ,h0gp_d,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.2,T52,F12.2,T71,A3)')                      &
                                          'hincgp  ',hincgp  ,hincgp_d,' R '
  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(T8,A)') 'Selected grid points:'

  WRITE (nuspecif, '(A)') '       Nr.    i-index   j-index    lon    lat    name'
  DO n = 1, ngp_tot
    IF (stationlist_tot(n)%rlatgp == rhuge) THEN
      WRITE (nuspecif,'(T10,A,I3.3,A,2X,2(I10,2X),2(F10.3,2X),A)')                           &
             'stationlist_tot(',n,')  ', stationlist_tot(n)%igp, stationlist_tot(n)%jgp,      &
             -999.999, -999.999,                                               &
             stationlist_tot(n)%ystation_name
    ELSE
      WRITE (nuspecif,'(T10,A,I3.3,A,2X,2(I10,2X),2(F10.3,2X),A)')                           &
             'stationlist_tot(',n,')  ', stationlist_tot(n)%igp,    stationlist_tot(n)%jgp,   &
                stationlist_tot(n)%rlatgp, stationlist_tot(n)%rlongp,        &
                stationlist_tot(n)%ystation_name
    ENDIF
  ENDDO
  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(T8,A)') 'Default selected grid points:'

  WRITE (nuspecif, '(A)') '       Nr.    i-index   j-index    lon    lat    name'
  DO n = 1, ngp_tot_d
    IF (stationlist_tot_d(n)%rlatgp == rhuge) THEN
      WRITE (nuspecif,'(T10,A,I3.3,A,2X,2(I10,2X),2(F10.3,2X),A)')                           &
             'stationlist_tot_d(',n,')', stationlist_tot(n)%igp, stationlist_tot(n)%jgp,      &
             -999.999, -999.999,                                               &
             stationlist_tot(n)%ystation_name
    ELSE
      WRITE (nuspecif,'(T10,A,I3.3,A,2X,2(I10,2X),2(F10.3,2X),A)')                           &
             'stationlist_tot_d(',n,')', stationlist_tot(n)%igp,    stationlist_tot(n)%jgp,   &
                stationlist_tot(n)%rlatgp, stationlist_tot(n)%rlongp,        &
                stationlist_tot(n)%ystation_name
    ENDIF
  ENDDO
  WRITE (nuspecif, '(A2)')  '  '

  ! mean values
  WRITE (nuspecif, '(A18)') '      Mean values:'
  WRITE (nuspecif, '(T7,A,T33,A,T52,A,T70,A)') 'Variable', 'Actual Value',   &
                                               'Default Value', 'Format'

  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                   'n0meanval',n0meanval  ,n0meanval_d,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                             'nincmeanval',nincmeanval  ,nincmeanval_d,' I '
  WRITE (nuspecif, '(A2)')  '  '

  ! calculation of integrals
  WRITE (nuspecif, '(T7,A)') 'Volume or area integrals:'
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                       &
                           'l_integrals ',l_integrals ,l_integrals_d, ' L '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                       &
                           'imin_integ  ',imin_integ  ,imin_integ_d , ' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                       &
                           'imax_integ  ',imax_integ  ,imax_integ_d , ' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                       &
                           'jmin_integ  ',jmin_integ  ,jmin_integ_d , ' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                       &
                           'jmax_integ  ',jmax_integ  ,jmax_integ_d , ' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                       &
                           'kmin_integ  ',kmin_integ  ,kmin_integ_d , ' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                       &
                           'kmax_integ  ',kmax_integ  ,kmax_integ_d , ' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                       &
                  'itype_diag_t2m', itype_diag_t2m, itype_diag_t2m_d, ' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                       &
            'itype_diag_gusts', itype_diag_gusts, itype_diag_gusts_d, ' I '

  WRITE (nuspecif, '(A2)')  '  '

  ! printing of testsuite
  WRITE (nuspecif, '(A24)') '      Testsuite Output: '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                       &
                            'ltestsuite', ltestsuite,  ltestsuite_d, ' L '
ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE input_diactl

!==============================================================================

!------------------------------------------------------------------------------
! End of external procedure organize_diagnosis                              
!------------------------------------------------------------------------------

END SUBROUTINE organize_diagnosis
