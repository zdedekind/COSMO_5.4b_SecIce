!+ External procedure organize_eps for organizing the ensemble prediction mode
!------------------------------------------------------------------------------

SUBROUTINE organize_eps(yaction, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!
!   This procedure is the driving routine for EPS mode set-up
!
!   Internal Routines/Functions currently contained:
!
!    - subroutine input_eps
!      This subroutine reads, checks and distributes the NAMELIST input
!      associated with the EPS
!
! Method:
!   Determine whether a certain action has to be performed in this time step
!
! Current Code Owner: DWD/FE15, Susanne Theis
!  phone:  +49  069 8062 2741
!  fax:    +49  69  8062 3721
!  email:  susanne.theis@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 3.21       2006/12/04 Susanne Theis
!  Initial release
! V4_5         2008/09/10 Christoph Gebhardt
!  Add namelist parameters for modifying values of lai, plcov, rootdp
! V4_8         2009/02/16 Ulrich Schaettler
!  Corrections in the settings of default values
! V4_10        2009/09/11 Ulrich Schaettler
!  Corrections in the settings of rmin_plcov
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V5_1         2014-11-28 Ulrich Blahak, Oliver Fuhrer, Lucio Torrisi, Christoph Schraff
!  Changed the format of some YUSPECIF entries for the CLM namelist tool.
!  Implemented F2003 IOMSG-mechanism for better namelist error messages. (UB)
!  Replaced ireals by wp (working precision) (OF)
!  Namelist parameters added controlling the stochastic perturbation of physics
!    tendencies scheme (SPPT)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

USE data_parameters,    ONLY :  &
    iintegers,        & ! KIND-type parameter for integer variables
    wp                  ! KIND-type parameters for real variables

USE parallel_utilities, ONLY :  &
    distribute_values     ! subroutine for MPI value distribution

USE data_parallel,      ONLY :  &
    my_world_id,   & ! rank of this subdomain in the global communicator
    icomm_world,   & ! communicator that belongs to igroup_world, i.e.
                     ! = MPI_COMM_WORLD
    nproc,         & ! total number of processors
    intbuf,        & ! buffers for distributing the Namelist
    realbuf,       & !
    logbuf,        & !
    imp_reals,     & !correct type for MPI
    imp_integers,  & !
    imp_logical

USE data_runcontrol,    ONLY :  &
    ldfi,          & ! 
    l2tls,         & ! forecast with 2-TL integration scheme
    lsppt,         & ! switch: if .true., perturb physics tendencies stochastic.
    iepsmem,       & ! ID of the member in the ensemble (ID >= 0) 
    iepstyp,       & ! ensemble type ID
    iepstot,       & ! total number of EPS members
    fac_plcov,     & ! modification factor for PLCOV
    rmin_plcov,    & ! lower limit of PLCOV
    rmax_plcov,    & ! upper limit of PLCOV
    fac_rootdp,    & ! modification factor for ROOTDP
    rmin_rootdp,   & ! lower limit of ROOTDP
    rmax_rootdp,   & ! upper limit of ROOTDP
    fac_lai,       & ! modification factor for LAI
    rmin_lai,      & ! lower limit of LAI
    rmax_lai,      & ! upper limit of LAI
    lhorint_rn,    & ! random numbers (defined on a rn horiz. coarse grid)
                     !   horizontally interpolated on model grid (otherwise
                     !   model grid points contained in the same rn coarse grid
                     !   point have the same random number value)
    ltimeint_rn,   & ! 3-D fields of random numbers are interpolated linearly
                     !   in time, otherwise remain constant in time until the
                     !   validity time of the next 3-D field of random numbers
    lgauss_rn,     & ! use a gaussian distribution of random numbers
                     !   (otherwise a uniform distribution is used)
    imode_rn,      & ! = 0: use only 1 stream of rn for all rn time steps
                     ! = 1: use a new stream of rn for every rn time step,
                     !      this enables temporal corrations in DA cycles
    npatmax,       & ! max. number of rn patterns with different scale/correl.
    npattern_rn,   & ! number of rn patterns with different scale/correlation
    hinc_rn,       & ! time increment (in [hrs]) for drawing a new field of
                     !   random numbers (used if (imode_rn > 0))
    nseed_rn,      & ! external part of seed for random number generation
    nseed_rn2,     & ! external seed to generate rn at ntstep>0 (if imode_rn>0)
    dlat_rn,       & ! random number coarse grid point distance in
                     !   meridional direction (in degrees)
    dlon_rn,       & ! random number coarse grid point distance in
                     !   zonal direction (in degrees)
    stdv_rn,       & ! standard deviation of the gaussian distribution of 
                     !   random numbers
    range_rn,      & ! max magnitude of random numbers 
    rvtaper_rn,    & ! externally specified function for vertical tapering of
                     !  the random number
    itype_vtaper_rn,&! type of tapering near surface and in stratosphere
    itype_qxpert_rn,&! defines which hydrometeor tendencies are perturbed
    itype_qxlim_rn,& ! type of reduction/removal of the perturbation 
                     !   in case of negative (qv, qc, qi) or 
                     !   supersaturated (qv) values
    ie_rn,         & ! number of horiz. coarse grid points in zonal direction
    je_rn,         & ! number of horiz. coarse grid points in meridional dir.
                     !   where random numbers are defined
    nuspecif         ! file unit number of YUSPECIF

USE data_modelconfig, ONLY :   &
    dlon,          & ! grid point distance in zonal direction (in degrees)
    dlat,          & ! grid point distance in meridional direction (in degrees)
    ke_tot,        & ! number of grid points in vertical direction
    dt               ! long time-step

USE vgrid_refatm_utils, ONLY:  &
    khmax            ! max. number of vertical coordinate parameters

USE src_stoch_physics,  ONLY :  &
    init_stoch_phys ,& ! initializes rn generator, produces first rn field
    compute_stoch_phys ! computes a 3-D random field for SPPT for current time

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Parameter list:
CHARACTER (LEN= *),       INTENT(IN)            ::                      &
  yaction      ! action to be performed

INTEGER (KIND=iintegers), INTENT(OUT)           ::                      &
  ierror       ! error status

CHARACTER (LEN= *),       INTENT(OUT)           ::                      &
  yerrmsg      ! error message


! local variables:

  INTEGER (KIND=iintegers) ::                       &
    nuin,                     &! unit-number of INPUT-File
    izerrstat                  ! error control


  CHARACTER (LEN= 9)       ::                       &
    yinput                     ! Namelist INPUT file

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine organize_eps
!------------------------------------------------------------------------------

ierror  = 0_iintegers
yerrmsg = '   '

!------------------------------------------------------------------------------
! Section 1: Input of the Namelist
!------------------------------------------------------------------------------

IF (yaction == 'input') THEN

  ! Open NAMELIST-INPUT file
  IF (my_world_id == 0) THEN
    PRINT *,'    INPUT OF THE NAMELISTS FOR ENSEMBLE'
    yinput   = 'INPUT_EPS'
    nuin     =  1
    OPEN(nuin   , FILE=yinput  , FORM=  'FORMATTED', STATUS='UNKNOWN',  &
         IOSTAT=izerrstat)
    IF(izerrstat /= 0) THEN
      yerrmsg  = ' ERROR    *** Error while opening file INPUT_EPS *** '
      ierror   = 2001
      RETURN
    ENDIF
  ENDIF

  ! Read the NAMELIST-group
  CALL input_epsctl (nuspecif, nuin, izerrstat)

  IF (my_world_id == 0) THEN
    ! Close file for input of the NAMELISTS
    CLOSE (nuin    , STATUS='KEEP')
  ENDIF

  IF (izerrstat < 0) THEN
    yerrmsg = 'ERROR *** while reading NAMELIST Group /EPSCTL/ ***'
    ierror  = 2003
  ELSEIF (izerrstat > 0) THEN
    yerrmsg = 'ERROR *** Wrong values occured in NAMELIST INPUT_EPS ***'
    ierror  = 2004
  ENDIF

ELSEIF (yaction == 'init') THEN

  IF (lsppt) THEN
        CALL init_stoch_phys
  ENDIF

ELSEIF (yaction == 'compute') THEN

  IF (lsppt) THEN
      CALL compute_stoch_phys
  ENDIF

ENDIF

!------------------------------------------------------------------------------
! Internal Procedures
!------------------------------------------------------------------------------

CONTAINS

!==============================================================================
!+ Internal procedure in "organize_diagnosis" for the input of NAMELIST diactl
!------------------------------------------------------------------------------

SUBROUTINE input_epsctl (nuspecif, nuin, ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine checks and distributes the NAMELIST input
!   associated with the EPS
!
! Method:
!   copied from organize_physics.f90 (input_phyctl)
!   and modified according to the current needs
!
!   All variables are initialized with default values and then read in from
!   the file INPUT. The input values are checked for errors and for
!   consistency. If wrong input values are detected the program prints
!   an error message. The program is not stopped in this routine but an
!   error code is returned to the calling routine that aborts the program
!   after reading in all other namelists.
!   In parallel mode, the variables are distributed to all nodes with the
!   environment-routine distribute_values.
!   Both, default and input values are written to the file YUSPECIF
!   (specification of the run).
!
!------------------------------------------------------------------------------

! Parameter list:
  INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
    nuspecif,     & ! Unit number for protocolling the task
    nuin            ! Unit number for Namelist INPUT file

  INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
    ierrstat        ! error status variable

! local variables:

! local namelist variables
  INTEGER (KIND=iintegers) ::                       &
    ninc_rn (npatmax)   ! timestep increment for drawing a new field of
                        !   random numbers (used if (imode_rn == 0))
! variables for default values
  INTEGER (KIND=iintegers) ::                       &
    iepsmem_d,        & ! Ensemble member-ID
    iepstot_d,        & ! total number of ensemble members
    iepstyp_d,        & ! Ensemble typ-ID
    imode_rn_d,       & ! = 0: use only 1 stream of rn for all rn time steps
                        ! = 1: use a new stream of rn for every rn time step,
                        !      this enables temporal corrations in DA cycles
    npattern_rn_d,    & ! number of rn patterns with different scale/correlation
    ninc_rn_d,        & ! timestep increment for drawing a new field of 
                        !   random numbers (used if (imode_rn == 0))
    nseed_rn_d,       & ! external part of seed for random number generation
    nseed_rn2_d,      & ! external seed to generate rn at ntstep>0 (imode_rn>0)
    itype_vtaper_rn_d,& ! type of tapering near surface and in stratosphere
    itype_qxpert_rn_d,& ! defines which hydrometeor tendencies are perturbed
    itype_qxlim_rn_d, & ! type of reduction/removal of the perturbation 
                        !   in case of negative (qv, qc, qi) or 
                        !   supersaturated (qv) values
    ie_rn_d,          & ! number of horiz. coarse grid points in zonal dir.
    je_rn_d             ! number of horiz. coarse grid points in meridional
                        !   direction where random numbers are defined

  REAL    (KIND=wp)        ::                       &
    fac_plcov_d,   & ! modification factor for PLCOV
    rmin_plcov_d,  & ! lower limit of PLCOV
    rmax_plcov_d,  & ! upper limit of PLCOV
    fac_rootdp_d,  & ! modification factor for ROOTDP
    rmin_rootdp_d, & ! lower limit of ROOTDP
    rmax_rootdp_d, & ! upper limit of ROOTDP
    fac_lai_d,     & ! modification factor for LAI
    rmin_lai_d,    & ! lower limit of LAI
    rmax_lai_d,    & ! upper limit of LAI
    dlat_rn_d,     & ! random number coarse grid point distance in
                     !   meridional direction (in degrees)
    dlon_rn_d,     & ! random number coarse grid point distance in
                     !   zonal direction (in degrees)
    stdv_rn_d,     & ! standard deviation of the gaussian distribution of 
                     !   random numbers
    range_rn_d,    & ! max magnitude of random numbers 
    vtaper_rn_d,   & ! specified function for vertical tapering of the random
                     !   numbers (only in stratosphere)
    hinc_rn_d        ! time increment (in [hrs]) for drawing a new field of
                     !   random numbers (used if (imode_rn > 0))

  LOGICAL ::       &
    lhorint_rn_d,  & ! random numbers (defined on a rn hor. coarse grid)
                     ! horizontally interpolated on model grid (otherwise
                     ! model grid points contained in the same rn coarse grid
                     ! point have the same random number value)
    ltimeint_rn_d, & ! random numbers (defined on a rn hor. coarse grid) are 
                     ! interpolated in time
    lgauss_rn_d      ! use a gaussian distribution of random numbers
                     ! (otherwise a uniform distribution is used)

! miscellaneous variables
  INTEGER (KIND=iintegers) ::                       &
    istat, &
    i,     &
    ierr             ! error variable for distribute_values
                     ! in this form useless (no input/output)

  REAL (KIND=wp)     ::                             &
    range_tot     ,& ! max magnitude of the random numbers at any grid point
    vtaper_rn (khmax)

  CHARACTER (LEN= 8)       ::                       &
    yinput           ! Namelist INPUT file

  CHARACTER(LEN=250)         :: iomsg_str

! Define the namelist group
  NAMELIST /epsctl/ iepsmem, iepstot, iepstyp, fac_lai, rmin_lai, rmax_lai,  &
                    fac_plcov, rmin_plcov, rmax_plcov,                       &
                    fac_rootdp, rmin_rootdp, rmax_rootdp,                    &
                    lhorint_rn, ltimeint_rn, itype_qxlim_rn, itype_qxpert_rn,&
                    npattern_rn, nseed_rn, nseed_rn2, imode_rn,              &
                    ninc_rn, hinc_rn, lgauss_rn, stdv_rn, range_rn,          &
                    itype_vtaper_rn, vtaper_rn, dlat_rn, dlon_rn

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE input_epsctl
!------------------------------------------------------------------------------

ierrstat = 0_iintegers
ALLOCATE (rvtaper_rn(ke_tot), STAT=istat)

IF (my_world_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 1: Initialize the default variables
!------------------------------------------------------------------------------

  iepsmem_d         = -1_iintegers
  iepstot_d         = -1_iintegers
  iepstyp_d         = -1_iintegers
  itype_vtaper_rn_d =  1_iintegers
  itype_qxpert_rn_d =  0_iintegers
  itype_qxlim_rn_d  =  0_iintegers
  ninc_rn_d         = -1_iintegers
  nseed_rn_d        =  0_iintegers
  nseed_rn2_d       =  0_iintegers
  npattern_rn_d     =  1_iintegers
  imode_rn_d        =  1_iintegers

  fac_plcov_d    = 1.0_wp
  rmin_plcov_d   = 0.0_wp
  rmax_plcov_d   = 1.0_wp
  fac_rootdp_d   = 1.0_wp
  rmin_rootdp_d  = 0.0_wp
  rmax_rootdp_d  = 2.0_wp
  fac_lai_d      = 1.0_wp
  rmin_lai_d     = 0.0_wp
  rmax_lai_d     = 8.0_wp
  hinc_rn_d      = 6.0_wp
  stdv_rn_d      = 0.5_wp
  range_rn_d     = 2.0_wp * stdv_rn_d
  dlat_rn_d      = 2.5_wp
  dlon_rn_d      = 2.5_wp
  vtaper_rn_d    = 1.0_wp

  ltimeint_rn_d  = .TRUE.
  lhorint_rn_d   = .TRUE.
  lgauss_rn_d    = .TRUE.

!------------------------------------------------------------------------------
!- Section 2: Initialize variables with defaults
!------------------------------------------------------------------------------

  iepsmem         = iepsmem_d
  iepstot         = iepstot_d
  iepstyp         = iepstyp_d
  itype_vtaper_rn = itype_vtaper_rn_d
  itype_qxpert_rn = itype_qxpert_rn_d
  itype_qxlim_rn  = itype_qxlim_rn_d
  ninc_rn         = ninc_rn_d
  nseed_rn(:)     = nseed_rn_d
  nseed_rn2(:)    = nseed_rn2_d
  npattern_rn     = npattern_rn_d
  imode_rn        = imode_rn_d

  fac_plcov      = fac_plcov_d
  rmin_plcov     = rmin_plcov_d
  rmax_plcov     = rmax_plcov_d
  fac_rootdp     = fac_rootdp_d
  rmin_rootdp    = rmin_rootdp_d
  rmax_rootdp    = rmax_rootdp_d
  fac_lai        = fac_lai_d
  rmin_lai       = rmin_lai_d
  rmax_lai       = rmax_lai_d
  hinc_rn        = hinc_rn_d
  dlat_rn        = dlat_rn_d
  dlon_rn        = dlon_rn_d
  stdv_rn        = stdv_rn_d
  range_rn       = range_rn_d
  vtaper_rn(:)   = vtaper_rn_d

  lgauss_rn      = lgauss_rn_d
  lhorint_rn     = lhorint_rn_d
  ltimeint_rn    = ltimeint_rn_d

!------------------------------------------------------------------------------
!- Section 3: Input of the namelist values
!------------------------------------------------------------------------------

  iomsg_str(:) = ' '
  READ (nuin, epsctl, IOSTAT=ierr, IOMSG=iomsg_str)

  IF (ierr /= 0) WRITE (*,'(A,A)') 'Namelist-ERROR EPSCTL: ', TRIM(iomsg_str)
ENDIF

IF (nproc > 1) THEN
  ! distribute error status to all processors
  CALL distribute_values  (ierr, 1, 0, imp_integers,  icomm_world, ierr)
ENDIF

IF (ierr /= 0) THEN
  ierrstat = -1
  RETURN
ENDIF

IF (my_world_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 4: Check values for errors and consistency
!------------------------------------------------------------------------------

  rvtaper_rn(:)   = vtaper_rn(1:ke_tot)

! Check whether values are OK for EPS parameters
  IF ( iepsmem < 0 .OR. iepstyp < 0 .OR. iepstot < 0 .OR. &
       iepsmem > iepstot ) THEN
    PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
    PRINT *,' *** iepsmem, iepstot, iepstyp *** ',iepsmem,iepstot,iepstyp
    ierrstat = 1002
    RETURN
  ENDIF

  IF ( fac_plcov < 0.0_wp .OR. rmin_plcov > rmax_plcov ) THEN
    PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
    PRINT *,' *** fac_plcov, rmin_plcov, rmax_plcov *** ',  &
            fac_plcov, rmin_plcov, rmax_plcov
    ierrstat = 1002
    RETURN
  ENDIF

  IF ( rmax_plcov > 1.0_wp ) THEN
    PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
    PRINT *,' *** rmax_plcov >1.0 is not allowed *** ', rmax_plcov
    ierrstat = 1002
    RETURN
  ENDIF

  IF ( fac_rootdp < 0.0_wp .OR. rmin_rootdp > rmax_rootdp ) THEN
    PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
    PRINT *,' *** fac_rootdp, rmin_rootdp, rmax_rootdp *** ',  &
            fac_rootdp, rmin_rootdp, rmax_rootdp
    ierrstat = 1002
    RETURN
  ENDIF

  IF ( fac_lai < 0.0_wp .OR. rmin_lai > rmax_lai ) THEN
    PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
    PRINT *,' *** fac_lai, rmin_lai, rmax_lai *** ',  &
            fac_lai, rmin_lai, rmax_lai
    ierrstat = 1002
    RETURN
  ENDIF

! Check whether values are OK for SPPT parameters
  IF (iepsmem < 0 .AND. lsppt) THEN
    PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
    PRINT *,' *** iepsmem is needed, if lsppt = T *** ',  &
       iepsmem
    ierrstat = 1002
    RETURN
  ENDIF
  IF (ldfi .AND. lsppt) THEN
    PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
    PRINT *,' *** with stoch phys -> ldfi = F *** '
    ierrstat = 1002
    RETURN
  ENDIF
  IF (.NOT.l2tls .AND. lsppt) THEN
    PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
    PRINT *,' *** with stoch phys -> l2tls = T *** '
    ierrstat = 1002
    RETURN
  ENDIF

  IF (.NOT. lsppt) THEN
    npattern_rn = 0
  ELSE
    IF ((itype_qxpert_rn < 0 .OR. itype_qxpert_rn > 2)) THEN
      PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
      PRINT *,' *** specify 0<=itype_qxpert_rn<=2  *** '
      ierrstat = 1002
      RETURN
    ENDIF
    IF ((itype_qxlim_rn < 0 .OR. itype_qxlim_rn > 1)) THEN
      PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
      PRINT *,' *** specify 0<=itype_qxlim_rn<=1  *** '
      ierrstat = 1002
      RETURN
    ENDIF
    IF ((itype_vtaper_rn < 0 .OR. itype_vtaper_rn > 3)) THEN
      PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
      PRINT *,' *** specify 0<=itype_vtaper_rn<=3  *** '
      ierrstat = 1002
      RETURN
    ENDIF
    IF (itype_vtaper_rn==3 .AND. (     ANY(rvtaper_rn(:) < 0.0_wp)      &
                                  .OR. ANY(rvtaper_rn(:) > 1.0_wp))) THEN
      PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
      PRINT *,' *** with itype_vtaper_rn==3 specify 0<=vtaper_rn<=1  *** '
      ierrstat = 1002
      RETURN
    ENDIF
    IF ((npattern_rn > npatmax) .OR. (npattern_rn < 1)) THEN
      PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
      PRINT *,' *** npattern_rn =', npattern_rn,' > npatmax = ', npatmax         &
             ,' or < 1 not allowed *** '
      ierrstat = 1002
      RETURN
    ENDIF
    imode_rn = MAX( imode_rn, 0 )
  ENDIF

  !   note: this loop is automatically skipped if .not.lsppt
  DO i=1,npattern_rn
    !   the namelist value of 'hinc_rn' is preferred over that of 'ninc_rn'
    IF (hinc_rn(i) <= 0.0001_wp)  hinc_rn(i)  =  dt* ninc_rn(i) /3600.0_wp
    !   'hinc_rn' is rounded to the nearest minute to allow for exact time steps
    !   of some fractions of hours which do not have a finite decimal expansion
    !   such as 1/3 hrs (20 min), 1/12 hrs (5 min), etc.
    hinc_rn(i)  =  NINT( hinc_rn(i) *60.0_wp) /60.0_wp
    !   check positiveness of ninc_rn, hinc_rn (both exist and correspoond)
    IF (hinc_rn(i) <= 0.0001_wp) THEN
      PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
      PRINT *,' *** n, ninc_rn, hinc_rn *** ',i, ninc_rn(i), hinc_rn(i)
      ierrstat = 1002
      RETURN
    ENDIF
    IF (imode_rn >  0) THEN
      !   it is required that hinc_rn is chosen such that 24 hours is
      !     divisible by hinc_rn (where hinc_rn need not be an integer);
      !     and it is assumed that hinc_rn is a multiple of 1 second
      !   (otherwise, the random number fields, that are produced in
      !    subsequent runs for the same time, differ from each other,
      !    and this destroys temporal corralations of the random fields)
      IF (MOD( 86400, NINT( hinc_rn(i) *3600.0_wp) ) > 0) THEN
        PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
        PRINT *,' *** 24 hours not divisible by hinc_rn(',i,') = ',hinc_rn(i)  &
               ,' *** '
        ierrstat = 1002
        RETURN
      ENDIF
    ENDIF
    IF ((nseed_rn(i) < 0) .OR. ((nseed_rn2(i) < 0) .AND. (imode_rn >  0))) THEN
      PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
      IF (nseed_rn(i) < 0)  PRINT *,' *** nseed_rn(', i,') = ',nseed_rn(i)     &
                                   ,' must be positive *** '
      IF (nseed_rn2(i) < 0) PRINT *,' *** nseed_rn2(',i,') = ',nseed_rn2(i)    &
                                  ,' must be positive *** '
      ierrstat = 1002
      RETURN
    ENDIF
    IF ((dlat_rn(i) < dlat) .OR. (dlon_rn(i) < dlon)) THEN
      PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
      PRINT *,' *** , dlat_rn(',i,')/dlon_rn(',i,') must be >= dlat/dlon *** ' &
             , dlat_rn(i), dlon_rn(i), dlat, dlon
      ierrstat = 1002
      RETURN
    ENDIF
    IF ((range_rn(i) < 0.0_wp) .OR. (stdv_rn(i) < 0.0_wp)) THEN
      PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
      PRINT *,' *** , range_rn(',i,') and stdv_rn(',i,') must be > 0 *** '     &
             , range_rn(i), stdv_rn(i)
      ierrstat = 1002
      RETURN
    ENDIF
  ENDDO

  !   range_tot must be <= 1; otherwise the stochastic perturbations can change
  !   the sign of physics tendencies, which may easily cause model instabilities
  range_tot=0.0_wp
  DO i=1,npattern_rn
    range_tot = range_tot + (range_rn(i)**2)
  ENDDO
  range_tot = SQRT(range_tot)
  IF (range_tot > 1.0_wp) THEN
    PRINT *,' ERROR    *** Bad values in NAMELIST epsctl *** '
    PRINT *,' *** range_tot = SQRT( SUM( range_rn(n)^2 ) ) must be <= 1 *** '
    ierrstat = 1002
    RETURN
  ENDIF

ENDIF

!------------------------------------------------------------------------------
!- Section 5: Distribute variables to all nodes
!------------------------------------------------------------------------------

IF (nproc > 1) THEN

  IF (my_world_id == 0) THEN
    intbuf  ( 1) = iepsmem
    intbuf  ( 2) = iepstot
    intbuf  ( 3) = iepstyp
    intbuf  ( 4) = itype_vtaper_rn
    intbuf  ( 5) = itype_qxlim_rn
    intbuf  ( 6) = itype_qxpert_rn
    intbuf  ( 7) = npattern_rn
    intbuf  ( 8) = imode_rn

    realbuf  ( 1) = fac_plcov
    realbuf  ( 2) = rmin_plcov
    realbuf  ( 3) = rmax_plcov
    realbuf  ( 4) = fac_rootdp
    realbuf  ( 5) = rmin_rootdp
    realbuf  ( 6) = rmax_rootdp
    realbuf  ( 7) = fac_lai
    realbuf  ( 8) = rmin_lai
    realbuf  ( 9) = rmax_lai

    logbuf   ( 1) = lgauss_rn
    logbuf   ( 2) = lhorint_rn
    logbuf   ( 3) = ltimeint_rn
  ENDIF

  CALL distribute_values (intbuf,  8, 0, imp_integers, icomm_world, ierr)
  CALL distribute_values (realbuf, 9, 0, imp_reals,    icomm_world, ierr)
  CALL distribute_values (logbuf,  3, 0, imp_logical,  icomm_world, ierr)
  IF (my_world_id /= 0)  npattern_rn = intbuf  ( 7)
  CALL distribute_values (nseed_rn, npattern_rn,0,imp_integers,icomm_world,ierr)
  CALL distribute_values (nseed_rn2,npattern_rn,0,imp_integers,icomm_world,ierr)
  CALL distribute_values (hinc_rn,  npattern_rn,0, imp_reals, icomm_world, ierr)
  CALL distribute_values (dlat_rn,  npattern_rn,0, imp_reals, icomm_world, ierr)
  CALL distribute_values (dlon_rn,  npattern_rn,0, imp_reals, icomm_world, ierr)
  CALL distribute_values (stdv_rn,  npattern_rn,0, imp_reals, icomm_world, ierr)
  CALL distribute_values (range_rn, npattern_rn,0, imp_reals, icomm_world, ierr)
  CALL distribute_values (rvtaper_rn,    ke_tot,0, imp_reals, icomm_world, ierr)

  IF (my_world_id /= 0) THEN
    iepsmem  = intbuf  ( 1)
    iepstot  = intbuf  ( 2)
    iepstyp  = intbuf  ( 3)
    itype_vtaper_rn = intbuf  ( 4)
    itype_qxlim_rn  = intbuf  ( 5)
    itype_qxpert_rn = intbuf  ( 6)
!   npattern_rn     = intbuf  ( 7)
    imode_rn        = intbuf  ( 8)

    fac_plcov   = realbuf ( 1)
    rmin_plcov  = realbuf ( 2)
    rmax_plcov  = realbuf ( 3)
    fac_rootdp  = realbuf ( 4)
    rmin_rootdp = realbuf ( 5)
    rmax_rootdp = realbuf ( 6)
    fac_lai     = realbuf ( 7)
    rmin_lai    = realbuf ( 8)
    rmax_lai    = realbuf ( 9)

    lgauss_rn   = logbuf   ( 1)
    lhorint_rn  = logbuf   ( 2)
    ltimeint_rn = logbuf   ( 3)
  ENDIF

ENDIF

!------------------------------------------------------------------------------
!- Section 6: Output of the namelist variables and their default values
!------------------------------------------------------------------------------

IF (my_world_id == 0) THEN

  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(A23)') '0     NAMELIST:  epsctl'
  WRITE (nuspecif, '(A23)') '      -----------------'
  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(T7,A,T33,A,T51,A,T58,A)') 'Variable', 'Actual Value', &
                                               'Default Value', 'Format'
  WRITE (nuspecif, '(T8,A,T33,I12,T52,I12  ,T71,A3)')                      &
                               'iepsmem',iepsmem, iepsmem_d  ,' I '
  WRITE (nuspecif, '(T8,A,T33,I12,T52,I12  ,T71,A3)')                      &
                               'iepstot',iepstot, iepstot_d  ,' I '
  WRITE (nuspecif, '(T8,A,T33,I12,T52,I12  ,T71,A3)')                      &
                               'iepstyp',iepstyp, iepstyp_d  ,' I '
  WRITE (nuspecif, '(T8,A,T33,F12.2,T52,F12.2  ,T71,A3)')                  &
                               'fac_plcov',fac_plcov, fac_plcov_d  ,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.2,T52,F12.2  ,T71,A3)')                  &
                            'rmin_plcov',rmin_plcov, rmin_plcov_d  ,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.2,T52,F12.2  ,T71,A3)')                  &
                            'rmax_plcov',rmax_plcov, rmax_plcov_d  ,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.3,T52,F12.3  ,T71,A3)')                  &
                            'fac_rootdp',fac_rootdp, fac_rootdp_d  ,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.3,T52,F12.3  ,T71,A3)')                  &
                         'rmin_rootdp',rmin_rootdp, rmin_rootdp_d  ,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.3,T52,F12.3  ,T71,A3)')                  &
                         'rmax_rootdp',rmax_rootdp, rmax_rootdp_d  ,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.2,T52,F12.2  ,T71,A3)')                  &
                                     'fac_lai',fac_lai, fac_lai_d  ,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.2,T52,F12.2  ,T71,A3)')                  &
                                  'rmin_lai',rmin_lai, rmin_lai_d  ,' R '
  WRITE (nuspecif, '(T8,A,T33,F12.2,T52,F12.2  ,T71,A3)')                  &
                                  'rmax_lai',rmax_lai, rmax_lai_d  ,' R '

  WRITE (nuspecif, '(T8,A,T33,I12,T52,I12  ,T71,A3)')                        &
                'itype_vtaper_rn',itype_vtaper_rn, itype_vtaper_rn_d  ,' I '
  WRITE (nuspecif, '(T8,A,T33,I12,T52,I12  ,T71,A3)')                        &
                'itype_qxpert_rn',itype_qxpert_rn, itype_qxpert_rn_d  ,' I '
  WRITE (nuspecif, '(T8,A,T33,I12,T52,I12  ,T71,A3)')                        &
                'itype_qxlim_rn', itype_qxlim_rn,  itype_qxlim_rn_d   ,' I '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                           'lgauss_rn'  , lgauss_rn  , lgauss_rn_d    ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                           'ltimeint_rn', ltimeint_rn, ltimeint_rn_d  ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                           'lhorint_rn' , lhorint_rn , lhorint_rn_d   ,' L '
  WRITE (nuspecif, '(T8,A,T33,I12,T52,I12  ,T71,A3)')                        &
                           'imode_rn   ', imode_rn   , imode_rn_d     ,' I '
  WRITE (nuspecif, '(T8,A,T33,I12,T52,I12  ,T71,A3)')                        &
                           'npattern_rn', npattern_rn, npattern_rn_d  ,' I '
  DO i=1,npattern_rn
    WRITE (nuspecif, '(T8,A,I3,A,T33,I12,T52,I12  ,T71,A3)')                 &
                           'nseed_rn(',i,')', nseed_rn(i), nseed_rn_d  ,' I '
    IF (imode_rn > 0)                                                        &
      WRITE (nuspecif, '(T8,A,I3,A,T33,I12,T52,I12  ,T71,A3)')               &
                           'nseed_rn2(',i,')',nseed_rn2(i),nseed_rn2_d ,' I '
    WRITE (nuspecif, '(T8,A,I3,A,T33,F12.2,T52,F12.2  ,T71,A3)')             &
                           'hinc_rn(',i,')' , hinc_rn(i) , hinc_rn_d   ,' R '
    WRITE (nuspecif, '(T8,A,I3,A,T33,F12.2,T52,F12.2  ,T71,A3)')             &
                           'dlat_rn(',i,')' , dlat_rn(i) , dlat_rn_d   ,' R '
    WRITE (nuspecif, '(T8,A,I3,A,T33,F12.2,T52,F12.2  ,T71,A3)')             &
                           'dlon_rn(',i,')' , dlon_rn(i) , dlon_rn_d   ,' R '
    WRITE (nuspecif, '(T8,A,I3,A,T33,F12.2,T52,F12.2  ,T71,A3)')             &
                           'stdv_rn(',i,')' , stdv_rn(i) , stdv_rn_d   ,' R '
    WRITE (nuspecif, '(T8,A,I3,A,T33,F12.2,T52,F12.2  ,T71,A3)')             &
                           'range_rn(',i,')', range_rn(i), range_rn_d  ,' R '
  ENDDO
  DO i=1,ke_tot
    WRITE (nuspecif, '(T8,A,I3,A,T33,F12.2,T52,F12.2  ,T71,A3)')             &
                        'vtaper_rn(',i,')',rvtaper_rn(i),vtaper_rn_d  ,' R '
  ENDDO
  WRITE (nuspecif, '(A2)')  '  '

ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE input_epsctl

!==============================================================================

!------------------------------------------------------------------------------
! End of external procedure organize_eps
!------------------------------------------------------------------------------

END SUBROUTINE organize_eps
