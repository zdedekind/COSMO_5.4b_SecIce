!+ Source Module providing utility routines for setting I/O meta data
!==============================================================================

MODULE io_metadata

!==============================================================================
!
! Description:
!  This module provides routines to set meta data for I/O
!
!  Routines (module procedures) currently contained:
!
!    - make_grib_init:       initialization of constant grib meta data
!    - make_grib_grid:       setting additional product dependent grid meta data
!    - make_grib_product:    setting product meta data
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version      Date       Name
! ----------   ---------- ----
! V4_28        2013/07/12 Ulrich Schaettler
!  Initial release
! V4_29        2013-10-02 Ulrich Schaettler, Ulrich Blahak
!  For the COSMO-Model only use vcoord and refatm from vgrid_refatm_utils
!  SR make_grib_product: move block to set productDefinitionTemplateNumber
!  downwards to properly code the ensembles in GRIB2
!  GRIB2: Adaptations to write endStep correct for quarterly output
!    and also set stepUnits for nunit_of_time=0 (Minutes in GRIB2)
!  Correct transformation from old to new style of coding GRIB1 gds (UB)
! V4_30        2013/11/08 Ulrich Schaettler
!  Renamed ipds to ipds_out to reflect usage for output
! V5_1         2014-11-28 Ulrich Schaettler, Ulrich Blahak, Oliver Fuhrer
!  Set GRIB2 typeOfGeneratingProcess for surface analysis products to 0
!    instead of 202 (Nugding)
!  SR make_grib_init: set ipds_out variables in any case, which is necessary, 
!   because there could be several GRIBOUT blocks and original settings could
!   be overwritten
!  GRIB2: set shortname and productDefinitionTemplateNumber for generalVertical
!   as last keys to get a proper coding the keys for this level type
!   (same problem as in INT2LM)
!  Do not set leveltypes for the FLake variables. These are only set implicitely
!   by setting the shortname
!  Do not write vertical coordinates to GRIB1 metadata, when using grib_api
!  Do not write refatm_id, vcoord_id to localInformationNumber of P, PP, HHL, resp.
!  Write height of levels above mean sea level to typeOfSecondFixedSurface for HHL
!  Added final coding of other sections to make_grib_product 
!     (was in src_output before)
!  Use nrbit_phhl for grib(2) packing of P and HHL
!  Initializing of pv_out must not be done for artificial data (UB)
!  Replaced ireals by wp (working precision) (OF)
! V5_2         2015-05-21 Burkhardt Rockel
!  Bug fix for writing variables for multi-layer snow model
! V5_3         2015-10-09 Doerte Liermann
!   grib_api meta data: modified
!    - setting of generating process identifier
!    - setting of local information number
!    - setting of Local Use Section
!      remove local section before changing the center and set it again afterwards
!    - use key "backgroundProcess" instead of backgroundGeneratingProcessIdentifier
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Modules used:

USE data_io,         ONLY:    &
  ydate_ini,         & ! start of the forecast 
  ncenter,           & ! originating center identification
  nsubcenter,        & ! originating sub-center identification
  lst_gribtabs,      & ! IDs of GRIB tables use
  nlocaldefnr,       & ! local definition number for GRIB local section
  nactlocdefnr,      & ! to overwrite Namelist parameter with some center default
  nprocess_ini_in,   & ! process gener. identification for initial (analysis)
  nprocess_bd_in,    & ! and for boundary (forecasts) data from input data
  nsma_stat,         & ! status for soil moisture analysis
  nrbit_phhl,        & ! packing Rate for P and HHL in GRIB2
  idims_out,         & ! array for all dimensions
  ipds_out,          & ! product definition section for output
  igds_out,          & ! grid description section
  ibms,              & ! bit map section
  ibds,              & ! binary data section
  inrvert_in,        & ! number of vertical coordinate parameters of input data
  inrvert_out,       & ! number of vertical coordinate parameters of output data
  l_ke_in_gds,       & ! explicit GDS entry for number of model levels
  l_ke_in_input,     & ! indicates whether GRIB1 input data contains ke in meta data
  pv_in,             & ! vertical coordinate parameters for input data
  pv_out,            & ! vertical coordinate parameters for output data
  ds_grib,           & ! array for unpacked data
  ds_real,           & ! array for unpacked data
  ldwd_grib_use,     & ! use some DWD specific Grib settings
  lbdclim,           & ! boundary data in climate model     ! PIK  (D.Hauffe)
  ylevltypes1,       & ! to convert GRIB1 level types to grib_api string typeOfLevel for GRIB1
  ylevltypes2,       & ! to convert GRIB1 level types to grib_api string typeOfLevel for GRIB2
  ysteptypes,        & ! to convert GRIB1 time range indicator to grib_api stinf stepType
  pp_nl,             & ! structure for gribout namelist
  var                  ! array for LM variable table

!==============================================================================

USE data_modelconfig,ONLY : &
    czmls,        & ! depth of the main soil layers in m
    czhls,        & ! depth of the half soil layers in m
    msoilgrib,    & ! grib coded depth of main soil levels in centimeters
    dlon,         & ! grid point distance in zonal direction (in degrees)
    dlat,         & ! grid point distance in meridional direction (in degrees)
    dt,           & ! long time-step
    ke,           & ! number of grid points in vertical direction
    pollon,       & ! longitude of the rotated north pole (in degrees, E>0)
    pollat,       & ! latitude of the rotated north pole (in degrees, N>0)
    polgam,       & ! angle between the north poles of the systems
    startlon_tot, & ! transformed longitude of the lower left grid point
                    ! of the total domain (in degrees, E>0)
    startlat_tot    ! transformed latitude of the lower left grid point
                    ! of the total domain (in degrees, N>0)

!==============================================================================

USE data_parallel,      ONLY :  &
    my_cart_id        ! rank of this subdomain in the cartesian communicator

!==============================================================================

USE data_parameters, ONLY :   &
  wp,        & ! KIND-type parameters for real variables
  iintegers, & ! KIND-type parameter for standard integer variables
  irealgrib, & ! KIND-type parameter for the real variables in the grib library
  intgribf,  & ! KIND-type parameter for the fortran files of the grib library
  intgribc,  & ! KIND-type parameter for the c files of the grib library
  iwlength     ! length of a integer word of the grib library in byte

!==============================================================================

USE data_runcontrol, ONLY : &
    nlastmxu,     & ! last step when vbmax was "nullified"
    nlastmxt,     & ! last step when tmin, tmax were "nullified"
    ntstep,       & ! actual time step
    nvers,        & ! version number of experiment for documentation
    lroutine,     & ! if .TRUE., run an operational forecast
    luseobs,      & ! switch for using observational data for nudging/nudgecast
    leps,         & ! switch ensemble mode on/off
    iepsmem,      & ! ID of ensemble member (EPS)
    iepstot,      & ! total number ensemble members (EPS)
    iepstyp,      & ! ID of ensemble generation type (EPS)
    ldfi,         & ! whether digital filtering or not
    lartif_data     ! forecast with self-defined artificial data

!==============================================================================

#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
USE data_satellites,          ONLY :  &
    sat_compute, num_sensors, nmsgchan
#endif

!==============================================================================

USE environment,              ONLY :  &
    model_abort 

USE vgrid_refatm_utils,       ONLY :  &
    refatm, vcoord, nfltvc, svc1, svc2

!==============================================================================

#ifdef GRIBAPI
! grib_api interface
USE grib_api
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================

CONTAINS

!==============================================================================
! set array pv_out of vertical coordinate parameters for output
!------------------------------------------------------------------------------

SUBROUTINE set_vcoord_refatm_out

!------------------------------------------------------------------------------
!
! Description:
!
!------------------------------------------------------------------------------

INTEGER (KIND=iintegers)  :: izerrstat, k
REAL (KIND=wp),     ALLOCATABLE   :: zpvloc(:)

!------------------------------------------------------------------------------

izerrstat = 0

!------------------------------------------------------------------------------
! Section 1: set pv_out
!------------------------------------------------------------------------------

IF (ALLOCATED(pv_in) .AND. (.NOT. lartif_data)) THEN

  IF (l_ke_in_input .EQV. l_ke_in_gds) THEN
    inrvert_out = inrvert_in

    ALLOCATE (pv_out(inrvert_out), STAT=izerrstat)

    pv_out(1:inrvert_out) = pv_in(1:inrvert_out)

  ELSEIF ((l_ke_in_gds) .AND. (.NOT. l_ke_in_input)) THEN

    ! change pv_in from old coding style to pv_out in new coding style
    IF (my_cart_id == 0) THEN
      PRINT *, '    Mapping old gds coding style to new one'
    ENDIF

    ! use a local pv, which is a bit bigger than pv_in (5 should be enough)
    ALLOCATE (zpvloc(inrvert_in+5), STAT=izerrstat)

    ! Reset pv_out
    ! The first two entries in output pv_out are: ivctype and ke, which are not 
    ! present in pv_in  
    IF (refatm%irefatm == 1) THEN
      zpvloc( 1) = REAL (vcoord%ivctype, wp)
    ELSE IF (refatm%irefatm == 2) THEN
      zpvloc( 1) = REAL (vcoord%ivctype+100, wp)
    ENDIF
    zpvloc( 2) = REAL (ke            , wp)

    ! the next values are for reference atmosphere and vcoord
    DO k = 3, 6+ke+1
      zpvloc(k) = pv_in(k-2)
    ENDDO

    inrvert_out = 6 + ke+1

    IF ( ANY( vcoord%ivctype == (/3,4/) ) ) THEN
      ! the first value after the vertical coordinate parameters in the old
      ! coding style is ivctype: in the new style this is already in pv(1)
      zpvloc (6+ke+1+1) = pv_in (4+ke+1+2)
      zpvloc (6+ke+1+2) = pv_in (4+ke+1+3)
      zpvloc (6+ke+1+3) = pv_in (4+ke+1+4)
      inrvert_out = 6 + ke+1 + 5
    ENDIF

    IF (refatm%irefatm == 2) THEN
      IF ( ALL( vcoord%ivctype /= (/3,4/) ) ) THEN
        zpvloc (6+ke+1+1) = 0.0_wp
        zpvloc (6+ke+1+2) = 0.0_wp
        zpvloc (6+ke+1+3) = 0.0_wp
      ENDIF
      zpvloc (6+ke+1+4) = pv_in (4+ke+1+5)
      zpvloc (6+ke+1+5) = pv_in (4+ke+1+6)
      inrvert_out = 6 + ke+1 + 5
    ENDIF

    IF (refatm%irefatm == 3) THEN
      IF ( ALL( vcoord%ivctype /= (/3,4/) ) ) THEN
        zpvloc (6+ke+1+1) = 0.0_wp
        zpvloc (6+ke+1+2) = 0.0_wp
        zpvloc (6+ke+1+3) = 0.0_wp
      ENDIF
      zpvloc (6+ke+1+4) = pv_in (4+ke+1+5)
      inrvert_out = 6 + ke+1 + 4
    ENDIF

    ! allocate pv_out and deallocate zpvloc again
    ALLOCATE (pv_out(inrvert_out), STAT=izerrstat)
    pv_out(1:inrvert_out) = zpvloc(1:inrvert_out)

    DEALLOCATE (zpvloc)

  ELSEIF ((.NOT. l_ke_in_gds) .AND. (l_ke_in_input)) THEN
    ! change pv_in from new coding style to pv_out in old coding style

    ! to use a local pv with the same size of pv_in is ok
    ALLOCATE (zpvloc(inrvert_in), STAT=izerrstat)

    DO k = 1, 4+ke+1
      zpvloc(k) = pv_in(k+2)
    ENDDO

    inrvert_out = 4 + ke+1

    IF ( ANY( vcoord%ivctype == (/3,4/) ) ) THEN
      ! the first value after the vertical coordinate parameters in the old
      ! coding style is ivctype: in the new style this is already in pv(1)
      IF (refatm%irefatm == 1) THEN
        zpvloc (4+ke+1+1) = REAL (vcoord%ivctype, wp)
      ELSE IF (refatm%irefatm == 2) THEN
        zpvloc (4+ke+1+1) = REAL (vcoord%ivctype+100, wp)
      ENDIF
      zpvloc (4+ke+1+2) = pv_in (6+ke+1+1)
      zpvloc (4+ke+1+3) = pv_in (6+ke+1+2)
      zpvloc (4+ke+1+4) = pv_in (6+ke+1+3)
      inrvert_out = 4 + ke+1 + 4
    ENDIF

    IF (refatm%irefatm == 2) THEN
      IF ( ALL( vcoord%ivctype /= (/3,4/) ) ) THEN
        IF (refatm%irefatm == 1) THEN
          zpvloc (4+ke+1+1) = REAL (vcoord%ivctype, wp)
        ELSE IF (refatm%irefatm == 2) THEN
          zpvloc (4+ke+1+1) = REAL (vcoord%ivctype+100, wp)
        ENDIF
        zpvloc (4+ke+1+2) = 0.0_wp
        zpvloc (4+ke+1+3) = 0.0_wp
        zpvloc (4+ke+1+4) = 0.0_wp
      ENDIF
      zpvloc (4+ke+1+5) = pv_in (6+ke+1+4)
      zpvloc (4+ke+1+6) = pv_in (6+ke+1+5)
      inrvert_out = 4 + ke+1 + 6
    ENDIF

    IF (refatm%irefatm == 3) THEN
      IF ( ALL( vcoord%ivctype /= (/3,4/) ) ) THEN
        IF (refatm%irefatm == 1) THEN
          zpvloc (4+ke+1+1) = REAL (vcoord%ivctype, wp)
        ELSE IF (refatm%irefatm == 2) THEN
          zpvloc (4+ke+1+1) = REAL (vcoord%ivctype+100, wp)
        ENDIF
        zpvloc (4+ke+1+2) = 0.0_wp
        zpvloc (4+ke+1+3) = 0.0_wp
        zpvloc (4+ke+1+4) = 0.0_wp
      ENDIF
      zpvloc (4+ke+1+5) = pv_in (6+ke+1+4)
      inrvert_out = 4 + ke+1 + 5
    ENDIF

    ! allocate pv_out and deallocate zpvloc again
    ALLOCATE (pv_out(inrvert_out), STAT=izerrstat)
    pv_out(1:inrvert_out) = zpvloc(1:inrvert_out)

    DEALLOCATE (zpvloc)

  ENDIF

ELSE   
  ! if pv_in is not allocated, we are running with artificial data and pv has to be set

  ! number of the vertical coordinate parameters
  IF (.NOT. l_ke_in_gds) THEN
    ! old style of coding the vertical coordinate parameters
    IF ( ANY( vcoord%ivctype == (/3,4/) ) ) THEN   ! SLEVE coordinates
      inrvert_out = 8 + ke+1
    ELSE                     ! NON SLEVE coordinates
      inrvert_out = 5 + ke+1
    ENDIF
    IF (refatm%irefatm == 2) inrvert_out = 10 + ke+1
  ELSE
    ! new style of coding the vertical coordinate parameters
    IF ( ANY( vcoord%ivctype == (/3,4/) ) ) THEN   ! SLEVE coordinates
      inrvert_out = 9 + ke+1
    ELSE                     ! NON SLEVE coordinates
      inrvert_out = 6 + ke+1
    ENDIF
    IF (refatm%irefatm == 2) inrvert_out = 11 + ke+1
  ENDIF

  ALLOCATE (pv_out(inrvert_out),   STAT = izerrstat)

  IF (.NOT. l_ke_in_gds) THEN
    ! old style of coding: first 4 values for the reference atmosphere,
    ! then the vertical coordinate parameters and last eventually
    ! additional parameters for the SLEVE coordinate
    pv_out( 1) = refatm%p0sl
    pv_out( 2) = refatm%t0sl
    pv_out( 3) = refatm%dt0lp
    pv_out( 4) = vcoord%vcflat

    IF     (vcoord%ivctype == 1) THEN
      DO k = 1, ke+1
        pv_out (4 + k) = vcoord%sigm_coord(k)
      ENDDO
    ELSEIF ( ANY( vcoord%ivctype == (/2,3,4/) ) ) THEN
      DO k = 1, ke+1
        pv_out (4 + k) = vcoord%vert_coord(k)
      ENDDO
    ENDIF

    ! But this must be coded as integer later on because of 
    ! backwards compatibility
    pv_out (4+ke+1+1) = REAL (vcoord%ivctype, wp)

    IF ( ANY( vcoord%ivctype == (/3,4/) ) ) THEN
      pv_out (4+ke+1+2) = svc1
      pv_out (4+ke+1+3) = svc2
      pv_out (4+ke+1+4) = REAL(nfltvc,wp)
    ENDIF

    IF (refatm%irefatm == 2) THEN ! Write parameters for new reference atmosphere
      pv_out (4+ke+1+1) = REAL(vcoord%ivctype+100, wp)
      IF ( ALL( vcoord%ivctype /= (/3,4/) ) ) THEN
        pv_out (4+ke+1+2) = 0.0_wp
        pv_out (4+ke+1+3) = 0.0_wp
        pv_out (4+ke+1+4) = 0.0_wp
      ENDIF
      pv_out (4+ke+1+5) = refatm%delta_t
      pv_out (4+ke+1+6) = refatm%h_scal
    ENDIF
  ELSE
    ! new style of coding: the first 2 values are the vertical coordinate
    ! type and the number of main levels used. Then 4 values for the
    ! reference atmosphere, the vertical coordinate parameters and
    ! last eventually additional parameters for the SLEVE coordinate
    IF (refatm%irefatm == 1) THEN
      pv_out ( 1) = REAL (vcoord%ivctype, wp)
    ELSE IF (refatm%irefatm == 2) THEN
      pv_out ( 1) = REAL (vcoord%ivctype+100, wp)
    ENDIF

    pv_out ( 2) = REAL (ke,  wp)
    pv_out ( 3) = refatm%p0sl
    pv_out ( 4) = refatm%t0sl
    pv_out ( 5) = refatm%dt0lp
    pv_out ( 6) = vcoord%vcflat

    IF     (vcoord%ivctype == 1) THEN
      DO k = 1, ke+1
        pv_out (6 + k) = vcoord%sigm_coord(k)
      ENDDO
    ELSEIF ( ANY( vcoord%ivctype == (/2,3,4/) ) ) THEN
      DO k = 1, ke+1
        pv_out (6 + k) = vcoord%vert_coord(k)
      ENDDO
    ENDIF

    IF ( ANY( vcoord%ivctype == (/3,4/) ) ) THEN
      ! SLEVE coordinate: 3 new parameters are written to GDS
      pv_out (6+ke+1+1) =      svc1
      pv_out (6+ke+1+2) =      svc2
      pv_out (6+ke+1+3) = REAL(nfltvc, wp)
    ENDIF

    IF (refatm%irefatm == 2) THEN ! Write parameters for new reference atmosphere
      IF ( ALL( vcoord%ivctype /= (/3,4/) ) ) THEN
        pv_out (6+ke+1+1) = 0.0_wp
        pv_out (6+ke+1+2) = 0.0_wp
        pv_out (6+ke+1+3) = 0.0_wp
      ENDIF
      pv_out (6+ke+1+4) = refatm%delta_t
      pv_out (6+ke+1+5) = refatm%h_scal
    ENDIF
  ENDIF

ENDIF

END SUBROUTINE set_vcoord_refatm_out

!==============================================================================
!==============================================================================
! initialization of constant grib meta data

SUBROUTINE make_grib_init (ptr_to_out)

!------------------------------------------------------------------------------
!
! Description:
!   This routine sets the constant meta data when using grib_api. 
!   For the grid definition, only few data have to be set, because most data
!   are taken from the samples.
!   For GRIB2, also the constant meta data for the product definition are set.
!
!------------------------------------------------------------------------------

! Subroutine arguments
TYPE(pp_nl),              INTENT(INOUT)  ::    &
  ! INOUT: because the generating process identifier might be reset here
  ptr_to_out           ! pointer to the namelist group

!------------------------------------------------------------------------------

! Local variables
INTEGER (KIND=iintegers)                 ::    &
  izgrbid, izgeneprocid, izbackprocid, izpollat, izpollon, izmodnvers,    &
  nzdate, nztime, nzsecond, nzstatus

REAL (KIND=wp)                           ::    &
  pollat_sp, pollon_sp

CHARACTER  (LEN=80)          :: yzerrmsg
CHARACTER  (LEN=25)          :: yzroutine

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initialization; Set local use information
!------------------------------------------------------------------------------

yzroutine   = 'make_grib_init'
yzerrmsg(:) = ' '
izgrbid     = ptr_to_out%igribapi_id

! Initialisation of ipds_out
ipds_out(:)  = -999999

! information about generating process: this is local use for every center
! Values can be read for every group GRIBOUT by 
!  - nprocess_ini_out: for specifying generating process for Nudging runs
!       (for DWD: if not read, value is taken from initial data)
!  - nprocess_bd_out:  for specifying generating process for forecast runs
!       (for DWD: if not read, value is taken from boundary data)

SELECT CASE (ncenter)
CASE (78)   ! DWD

  ! generating process identification
  IF (ptr_to_out%lanalysis .OR. ptr_to_out%lsfc_ana) THEN
    IF (ptr_to_out%nprocess_ini_out == -999999) THEN
      ptr_to_out%nprocess_ini_out =  nprocess_ini_in
    ENDIF
    ! has to be set in any case
    izgeneprocid =  ptr_to_out%nprocess_ini_out
  ELSE
    IF (ptr_to_out%nprocess_bd_out == -999999) THEN
      ptr_to_out%nprocess_bd_out =  nprocess_bd_in
    ENDIF
    ! has to be set in any case
    izgeneprocid =  ptr_to_out%nprocess_bd_out
  ENDIF

  ! background generating process identification (for GRIB2)
  izbackprocid = INT (nvers/16384_iintegers)

CASE DEFAULT

  IF (ptr_to_out%nprocess_ini_out == -999999) THEN
      izgeneprocid = nprocess_ini_in
  ELSE
      izgeneprocid = ptr_to_out%nprocess_ini_out
  ENDIF
  izbackprocid = 255

END SELECT

!------------------------------------------------------------------------------
! Section 2: GRIB1 Grid Description Section for DWD GRIB1 library
!------------------------------------------------------------------------------

! set these values in any case, because there is only one pds_out array for
! all GRIBOUT blocks
ipds_out(3)  = ncenter            ! centre identification
ipds_out(5)  = 255_intgribf       ! number of used grid
ipds_out(6)  = 128_intgribf       ! additional flags

SELECT CASE (ptr_to_out%yform_write)

!------------------------------------------------------------------------------
! Section 3: GRIB1 Grid Description Section for grib_api
!------------------------------------------------------------------------------

CASE ('api1')

#ifdef GRIBAPI
  ! here the constant gds- and pds-values are set for the sample of this outblock

  ! data representation type 
  CALL grib_set (izgrbid,'dataRepresentationType',           10)  ! igds_out( 4)

  ! calculation of the left bottom and right upper corner in millidegrees:
  ! igds_out(7, 8, 10, 11)
  ! This depends on the variable (U- and V-related variables are shifted in
  ! the Arakawa C-grid) and is determined during the output step.

  ! number of gridpoints 
  CALL grib_set (izgrbid,'Ni',            ptr_to_out%ie_out_tot)  ! igds_out( 5)
  CALL grib_set (izgrbid,'Nj',            ptr_to_out%je_out_tot)  ! igds_out( 6)
  CALL grib_set (izgrbid,'resolutionAndComponentFlags',       8)  ! igds_out( 9)

  ! increments: are set to 0 explicitly, because of igds_out(9) = 8
  CALL grib_set (izgrbid,'Di',                                0)  ! igds_out(12)
  CALL grib_set (izgrbid,'Dj',                                0)  ! igds_out(13)

  CALL grib_set (izgrbid,'scanningMode',                     64)  ! igds_out(14)
  ! igds_out(15:19) = 0

  ! coordinates of the pole: in the grib code the southern pole has to be
  ! specified: for the latitude this is the negative value, for the
  ! longitude it is the value + 180.0 and then limiting the value again
  ! to -180.0 ... 180.0
  pollon_sp = pollon + 180.0_wp
  IF (pollon_sp > 180.0_wp) THEN
    pollon_sp = pollon_sp - 360.0_wp
  ENDIF
  izpollat = NINT(-pollat      * 1000.0_wp)
  izpollon = NINT( pollon_sp   * 1000.0_wp)

  CALL grib_set (izgrbid, 'latitudeOfSouthernPole',    izpollat)  ! igds_out(20)
  CALL grib_set (izgrbid, 'longitudeOfSouthernPole',   izpollon)  ! igds_out(21)
  CALL grib_set (izgrbid, 'angleOfRotation',             polgam)  ! igds_out(22)

  ! vertical coordinate parameters
  ! these should be set per record; but because of historic reasons, we set it
  ! for all GRIB1 products
  ! (changed in Feb 2014: grib_api only writes vertical coordinates for
  ! leveltype 109/110)
! CALL grib_set (izgrbid,'PVPresent',     1)
! CALL grib_set (izgrbid,'NV',  inrvert_out)
! CALL grib_set (izgrbid,'pv',       pv_out)


  ! Product definition section
  CALL grib_set (izgrbid,'centre',                                     ncenter)
  CALL grib_set (izgrbid,'generatingProcessIdentifier',           izgeneprocid)

!------------------------------------------------------------------------------
! Section 4: Constant GRIB2 meta data
!------------------------------------------------------------------------------

CASE ('api2')

  ! Indicator Section
  ! -----------------

  ! local section has to be deleted first before changing the centre 
  !    -> it has to be defined afterwards with new localDefinitionNumber
  ! (to be able to use sample files with unknown local section, we do this for all
  !  centers, also for DWD)
  CALL grib_set (izgrbid, 'grib2LocalSectionPresent',        0)   ! delete local section

  CALL grib_set (izgrbid, 'centre',                    ncenter)   ! originating centre
  CALL grib_set (izgrbid, 'subCentre',              nsubcenter)   ! originating subcentre

  ! Identification Section
  ! ----------------------

  IF (ptr_to_out%lanalysis .OR. ptr_to_out%lsfc_ana) THEN
    CALL grib_set (izgrbid, 'significanceOfReferenceTime',     0)   ! Analysis
  ELSE
    CALL grib_set (izgrbid, 'significanceOfReferenceTime',     1)   ! start of forecast
  ENDIF

  ! Reference time of data
  ! In case of nudging runs, this is overwritten later on
  READ(ydate_ini( 1: 8),'(I8)') nzdate
  READ(ydate_ini( 9:12),'(I4)') nztime
  READ(ydate_ini(13:14),'(I2)') nzsecond
  CALL grib_set (izgrbid,'dataDate',                     nzdate)   ! yyyymmdd
  CALL grib_set (izgrbid,'dataTime',                     nztime)   ! hhmm
  CALL grib_set (izgrbid,'second',                     nzsecond)   ! ss

  ! Production Status: depends on the center and is computed in the local use section

  ! Type of processed data
  IF (leps) THEN
    ! this is set to "control and perturbed forecast products"
    CALL grib_set (izgrbid,'typeOfProcessedData',             5)
  ELSE
    IF (ptr_to_out%lanalysis .OR. ptr_to_out%lsfc_ana) THEN
      CALL grib_set (izgrbid,'typeOfProcessedData',             0)  ! Analysis products
    ELSE
      CALL grib_set (izgrbid,'typeOfProcessedData',             1)  ! Forecast products
    ENDIF
  ENDIF

  ! Local Use Section
  ! -----------------

  ! For every COSMO application we require a local use section
  CALL grib_set (izgrbid, 'grib2LocalSectionPresent',1)

  ! Set the local definition number: every center has to set something!
  SELECT CASE (nlocaldefnr)
  CASE (-1)
    ! local definition number not set per namelist: use DWD defaults:
    IF (leps) THEN
      IF (iepstyp >= 0) THEN
        CALL grib_set (izgrbid, 'localDefinitionNumber',         253)  ! ensemble BC
        CALL grib_set (izgrbid, 'localTypeOfEnsembleForecast', iepstyp)
        nactlocdefnr = 253
      ENDIF
    ELSE
      CALL grib_set (izgrbid,'localDefinitionNumber',            254)  ! deterministic system
      nactlocdefnr = 254
    ENDIF
  CASE (250,252,253,254)
    CALL grib_set (izgrbid,'localDefinitionNumber', nlocaldefnr)
    nactlocdefnr = nlocaldefnr
  CASE DEFAULT
    PRINT *, ' *** ERROR in localDefinitionNumber *** ', nlocaldefnr
    PRINT *, ' ***       This is not a valid definition number! *** '
    yzerrmsg = 'Error in settings for GRIB2 local use section'
    CALL model_abort (my_cart_id, 2022, yzerrmsg, yzroutine)
  END SELECT

  SELECT CASE (ncenter)
  CASE (78)
    ! Production Status: check nvers:
    izmodnvers = MODULO(nvers, 16384)

    SELECT CASE (izmodnvers)
    CASE (  1: 50)
      nzstatus = 0    ! Operational
    CASE ( 51: 99)
      nzstatus = 1    ! Operational Test
    CASE (100:16384)
      nzstatus = 2    ! Research
    END SELECT

    CALL grib_set (izgrbid, 'productionStatusOfProcessedData', nzstatus)

    ! local number of experiment: set "pure" experiment number from nvers
    CALL grib_set (izgrbid, 'localNumberOfExperiment',      izmodnvers)  ! GRIB1: ipds_out(47)

  CASE DEFAULT

    ! Production Status
    IF (lroutine) THEN
      CALL grib_set (izgrbid,'productionStatusOfProcessedData', 0)  ! Operational
    ELSE
      CALL grib_set (izgrbid,'productionStatusOfProcessedData', 2)  ! Research
    ENDIF

    ! local number of experiment: set "pure" experiment number from nvers
    CALL grib_set (izgrbid, 'localNumberOfExperiment',      nvers)  ! GRIB1: ipds_out(47)

  END SELECT

  ! local information number: is already set in shortName (if present)
  ! CALL grib_set (izgrbid, 'localInformationNumber',                0)  !  ipds_out(41)


  ! Grid Definition Section
  ! -----------------------

  ! Most keys are already set in the sample grib message

  ! Number of points Ni, Nj
  CALL grib_set (izgrbid,'Ni',          ptr_to_out%ie_out_tot)
  CALL grib_set (izgrbid,'Nj',          ptr_to_out%je_out_tot)

  ! start and end points have to be set per product because of the staggered grid

  ! resolution and component flags
  ! the data in the sample are set to:
  !   Bit 3 = 1: i direction increments given
  !   Bit 4 = 1: j direction increments given
  !   Bit 5 = 1: Resolved u- and v-components of vector quantities relative to the defined grid
  !     value of octet 55 therefore is: 56
  CALL grib_set (izgrbid,'ijDirectionIncrementGiven',       1)
  CALL grib_set (izgrbid,'uvRelativeToGrid',                1)

! US: if we set the increments here, the grib_clone later could round the value 
!     and lead to uncorrect values: set it with the first and last lat-lon-values
! CALL grib_set (izgrbid,'iDirectionIncrementInDegrees', dlon)
! CALL grib_set (izgrbid,'jDirectionIncrementInDegrees', dlat)

  ! scanning mode is set in the sample

  ! specifications of rotated pole
  ! convert COSMO north pole to rotated south pole
  pollat_sp = -pollat
  ! pollon is specified in the range -180.0...+180.0
  ! but GRIB2 gives values in the range 0.0...+360.0
  pollon_sp =  pollon + 180.0_wp
  ! above statement converts to south pole and automatically to range 0.0...+360.0
  ! so no additional correction as in GRIB1 is necessary

  CALL grib_set (izgrbid, 'latitudeOfSouthernPoleInDegrees',  pollat_sp)
  CALL grib_set (izgrbid, 'longitudeOfSouthernPoleInDegrees', pollon_sp)
  CALL grib_set (izgrbid, 'angleOfRotationInDegrees',         polgam)

  ! Product Definition (constant meta data
  ! ------------------

  CALL grib_set (izgrbid,'generatingProcessIdentifier',           izgeneprocid)
  CALL grib_set (izgrbid,'backgroundProcess', izbackprocid)
! CALL grib_set (izgrbid,'backgroundGeneratingProcessIdentifier', izbackprocid)

!  ! this is also local use
!  SELECT CASE (ncenter)
!  CASE (78)
!
!    ! Type of generating Process (there are some DWD defined table entries)
!    IF (leps) THEN
!      CALL grib_set (izgrbid,'typeOfGeneratingProcess',           4) ! Ensemble Forecast
!    ELSE
!      IF (ptr_to_out%lanalysis) THEN
!        CALL grib_set (izgrbid,'typeOfGeneratingProcess',       202) ! Nudging Analysis
!      ELSEIF (ptr_to_out%lsfc_ana) THEN
!        CALL grib_set (izgrbid,'typeOfGeneratingProcess',         0) ! External Analysis
!      ELSE
!        IF (ldfi) THEN
!          CALL grib_set (izgrbid,'typeOfGeneratingProcess',       1) ! Initialization
!        ELSE
!          IF (luseobs) THEN
!            CALL grib_set (izgrbid,'typeOfGeneratingProcess',   203) ! Nudgecast
!          ELSE
!            CALL grib_set (izgrbid,'typeOfGeneratingProcess',     2) ! Forecast
!          ENDIF
!        ENDIF
!      ENDIF
!    ENDIF
!
!  CASE DEFAULT
!
!    IF (leps) THEN
!      CALL grib_set (izgrbid,'typeOfGeneratingProcess',           4) ! Ensemble Forecast
!    ELSE
!      IF (ptr_to_out%lanalysis .OR. ptr_to_out%lsfc_ana) THEN
!        CALL grib_set (izgrbid,'typeOfGeneratingProcess',         0) ! Analysis
!      ELSE
!        IF (ldfi) THEN
!          CALL grib_set (izgrbid,'typeOfGeneratingProcess',       1) ! Initialization
!        ELSE
!          CALL grib_set (izgrbid,'typeOfGeneratingProcess',       2) ! Forecast
!        ENDIF
!      ENDIF
!    ENDIF
!
!  END SELECT
#endif

END SELECT

END SUBROUTINE make_grib_init

!==============================================================================
!==============================================================================
!+ Module procedure in io_metadata to create the product definition block
!------------------------------------------------------------------------------

SUBROUTINE make_grib_grid (ptr_to_out, igrbid, n1, n2, n3, yextension, idebug)

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure sets the still missing parts for the grib grid 
!   definition. It distinguishes between GRIB1 and GRIB2 and whether DWD libgrib
!   or grib_api is used.
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
! Arguments with intent(in):
TYPE(pp_nl),              INTENT(IN)     ::    &
  ptr_to_out        ! pointer to the namelist group

INTEGER  (KIND=iintegers) , INTENT(IN) :: n1, n2, n3, idebug, igrbid
CHARACTER(LEN=1),           INTENT(IN) :: yextension

! Local scalars:
INTEGER  (KIND=iintegers)              ::      &
  idlon, jdlat, istartlon, jstartlat, iendlon, jendlat, izreturn

REAL     (KIND=wp)                     ::      &
  zstartlon_tot, zendlon_tot, zstartlat_tot, zendlat_tot

CHARACTER(LEN=10)                      ::      &
  yzname

! Local scalars
INTEGER  (KIND=iintegers)              ::      &
  izlen

!- End of header
!==============================================================================

  ! Initialization:
  izreturn = 0
  yzname   = var(n1,n2,n3)%name

!------------------------------------------------------------------------------
! Section 1: Characteristics of (staggered) grid
!------------------------------------------------------------------------------

  IF (idebug > 15) THEN
    WRITE (*,'(A)') '  io_metadata, make_grib_grid: grid description'
  ENDIF

  ! calculation of the left bottom and upper right corner
  ! these depend also on the namelist group
  zstartlon_tot = startlon_tot + REAL (ptr_to_out%i_out_start-1,wp)*dlon
  zstartlat_tot = startlat_tot + REAL (ptr_to_out%j_out_start-1,wp)*dlat

  ! Careful: In case of p- or z-levels, u and v are interpolated to 
  !          to the mass grid point. The gds has not to be "staggered"
  !          in this case
  izlen  = LEN_TRIM(yzname)
  IF( ( (yzname(1:izlen) == 'U')            .OR.        &
        (yzname(1:izlen) == 'AUMFL_S') )    .AND.       &
      (.NOT. ptr_to_out%luvmasspoint)         .AND.     &
      (yextension /= 'p') .AND. (yextension /= 'z') ) THEN
    zstartlon_tot = zstartlon_tot + 0.5_wp * dlon
  ENDIF
  zendlon_tot = zstartlon_tot + (ptr_to_out%ie_out_tot-1)*dlon

  IF( ( (yzname(1:izlen) == 'V')            .OR.        &
        (yzname(1:izlen) == 'AVMFL_S') )    .AND.       &
      (.NOT. ptr_to_out%luvmasspoint)         .AND.     &
      (yextension /= 'p') .AND. (yextension /= 'z') ) THEN
    zstartlat_tot = zstartlat_tot + 0.5_wp * dlat
  ENDIF
  zendlat_tot = zstartlat_tot + (ptr_to_out%je_out_tot-1)*dlat

  IF (ptr_to_out%yform_write == 'api2') THEN
    ! For GRIB2 all longitude values have to be limited 
    ! to the range (0.0,360.0)
    IF (zstartlon_tot < 0.0_wp) THEN
      zstartlon_tot = zstartlon_tot + 360.0_wp
    ENDIF
    IF (zendlon_tot < 0.0_wp) THEN
      zendlon_tot = zendlon_tot + 360.0_wp
    ENDIF
  ELSE
    ! For GRIB1 all longitude values have to be limited 
    ! to the range (-180.0,+180.0)
    IF (zstartlon_tot > 180.0_wp) THEN
      zstartlon_tot = zstartlon_tot - 360.0_wp
    ENDIF
    IF (zendlon_tot > 180.0_wp) THEN
      zendlon_tot = zendlon_tot - 360.0_wp
    ENDIF
  ENDIF

  ! this is for GRIB1 (DWDLIB and grib_api)
  jstartlat = NINT(zstartlat_tot * 1000.0_wp)
  istartlon = NINT(zstartlon_tot * 1000.0_wp)
  jendlat   = NINT(zendlat_tot   * 1000.0_wp)
  iendlon   = NINT(zendlon_tot   * 1000.0_wp)

  SELECT CASE (ptr_to_out%yform_write)

  CASE ('grb1','bina')
 
    ! complete grid definition section
    ! --------------------------------

    ! number of gridpoints (depend on namelist group)
    igds_out( 5) = ptr_to_out%ie_out_tot
    igds_out( 6) = ptr_to_out%je_out_tot

    igds_out( 7) = jstartlat
    igds_out( 8) = istartlon
    igds_out(10) = jendlat
    igds_out(11) = iendlon


#ifdef GRIBAPI
  CASE ('api1')

    ! ie, je already set with constant meta data

    CALL grib_set (igrbid,'La1',           jstartlat)  ! igds_out( 7)
    CALL grib_set (igrbid,'Lo1',           istartlon)  ! igds_out( 8)
    CALL grib_set (igrbid,'La2',             jendlat)  ! igds_out(10)
    CALL grib_set (igrbid,'Lo2',             iendlon)  ! igds_out(11)

    ! vertical coordinate parameters are set per record now for writing GRIB1 with
    ! grib_api (only levtyp == 109/110 will have vertical coordinates)
    IF ( (var(n1,n2,n3)%levtyp == 109) .OR. (var(n1,n2,n3)%levtyp == 110) ) THEN
      CALL grib_set (igrbid,'PVPresent',     1)
      CALL grib_set (igrbid,'NV',  inrvert_out)
      CALL grib_set (igrbid,'pv',       pv_out)
    ELSE
      CALL grib_set (igrbid,'PVPresent',     0)
      CALL grib_set (igrbid,'NV',            0)
    ENDIF

  CASE ('api2')

    ! ie, je already set with constant meta data, but not increments, because 
    ! of precision loss during cloning

    ! As long as we do not know about ongoings in grib_api, we set the original
    ! keys with integer values here.
    ! Multiplying with 1E6_wp should not lead to rounding errors, as long as 
    ! wp (and the related variables dlon, dlat, etc.) are really double precision
    idlon     = NINT(dlon          * 1.0E6_wp, iintegers)
    jdlat     = NINT(dlat          * 1.0E6_wp, iintegers)
    jstartlat = NINT(zstartlat_tot * 1.0E6_wp, iintegers)
    istartlon = NINT(zstartlon_tot * 1.0E6_wp, iintegers)
    jendlat   = NINT(zendlat_tot   * 1.0E6_wp, iintegers)
    iendlon   = NINT(zendlon_tot   * 1.0E6_wp, iintegers)

    CALL grib_set (igrbid, 'iDirectionIncrement',            idlon )
    CALL grib_set (igrbid, 'jDirectionIncrement',            jdlat )
    CALL grib_set (igrbid, 'latitudeOfFirstGridPoint',   jstartlat )
    CALL grib_set (igrbid, 'longitudeOfFirstGridPoint',  istartlon )
    CALL grib_set (igrbid, 'latitudeOfLastGridPoint',      jendlat )
    CALL grib_set (igrbid, 'longitudeOfLastGridPoint',     iendlon )

!US CALL grib_set (igrbid, 'latitudeOfFirstGridPointInDegrees',  zstartlat_tot)
!US CALL grib_set (igrbid, 'longitudeOfFirstGridPointInDegrees', zstartlon_tot)
!US CALL grib_set (igrbid, 'latitudeOfLastGridPointInDegrees',     zendlat_tot)
!US CALL grib_set (igrbid, 'longitudeOfLastGridPointInDegrees',    zendlon_tot)
#endif
  END SELECT

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE make_grib_grid

!==============================================================================
!==============================================================================
!+ Module procedure in io_metadata to create the product definition block
!------------------------------------------------------------------------------

SUBROUTINE make_grib_product (ptr_to_out, igrbid, n1, n2, n3, nlevel, ydatact,       &
                           nlastout, yextension, slev, lrestart, idebug)

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure sets the still missing parts for the grib product
!   sections. It distinguishes between GRIB1 and GRIB2 and whether DWD libgrib
!   or grib_api is used.
!
!    - product definition:
!         => coding of product
!         => vertical information
!         => time and date
!         => statistical processing
!         => local information
!         => final coding of other sections (bds, bms)
!    
!------------------------------------------------------------------------------

! Subroutine / Function arguments
! Arguments with intent(in):
TYPE(pp_nl),              INTENT(IN)     ::    &
  ptr_to_out        ! pointer to the namelist group

INTEGER  (KIND=iintegers) , INTENT(IN) :: n1, n2, n3, nlevel, nlastout, igrbid, idebug
CHARACTER (LEN=14),         INTENT(IN) :: ydatact
REAL     (KIND=wp)        , INTENT(IN) :: slev 
CHARACTER (LEN=1),          INTENT(IN) :: yextension 
LOGICAL                   , INTENT(IN) :: lrestart

! Local scalars:
INTEGER  (KIND=iintegers)              ::                                    &
  ndgb, mmgb, jjgb, nhgb, jcgb, nmingb, ijj, imm, idd, ihh, imi, iss, izlen, &
  istartstep, iendstep, iprdeftnr, izednr, mlev, ndate, ntime, nsecond,      &
  izlevtyp, izlevtop, izlevbot, nztri, izreturn, ilin

CHARACTER (LEN= 14)    :: ydatref
CHARACTER (LEN=  8)    :: ydate
CHARACTER (LEN= 10)    :: ytime, yzsteptyp
CHARACTER (LEN= 25)    :: yzname
CHARACTER (LEN= 30)    :: yzlevtyp

!- End of header
!==============================================================================

  ! Initialization:
  izreturn = 0
  yzname   = '                         '

  ! set edition number
  SELECT CASE (ptr_to_out%yform_write)
  CASE ('grb1','api1','bina')
    izednr = 1
  CASE ('api2')
    izednr = 2
  END SELECT

  yzname = TRIM(var(n1,n2,n3)%name)
  izlen    = LEN_TRIM(TRIM(yzname))

!------------------------------------------------------------------------------
! Section 1: Basic definitions: stepType, productDefinitionTemplateNumber
!------------------------------------------------------------------------------

  IF (idebug > 15) THEN
    WRITE (*,'(A)') '  io_metadata, make_grib: product coding'
  ENDIF

  ! set steptype
  yzsteptyp = ysteptypes(var(n1,n2,n3)%ntri)

  IF (yzsteptyp == 'diff') THEN
    IF (TRIM(var(n1,n2,n3)%name) == 'TMIN_2M') THEN
      yzsteptyp = 'min'
    ELSE
      yzsteptyp = 'max'
    ENDIF
  ENDIF

  ! put values to grib meta data
  SELECT CASE (ptr_to_out%yform_write)
  CASE ('grb1','bina')

    ipds_out(2)  = lst_gribtabs(n3)   ! GRIB1 table number
    IF (ptr_to_out%lanalysis .OR. ptr_to_out%lsfc_ana) THEN
      ipds_out(4)  = ptr_to_out%nprocess_ini_out ! generating process identification
    ELSE
      ipds_out(4)  = ptr_to_out%nprocess_bd_out  ! generating process identification
    ENDIF
    ipds_out(7)  = n2                 ! element number / identification of parameter

#ifdef GRIBAPI
  CASE ('api1')

    ! number of used grid and additional flags are set in the sample
    CALL grib_set (igrbid,'table2Version',            lst_gribtabs(n3))
    CALL grib_set (igrbid,'indicatorOfParameter',                   n2)

  CASE ('api2')

    ! set the typeOfGeneratingProcess: this is not constant, because lanalysis / lsfc_ana
    ! are modified in the sfcana. But this is also local use.
    SELECT CASE (ncenter)
    CASE (78)

      ! Type of generating Process (there are some DWD defined table entries)
      IF (leps) THEN
        CALL grib_set (igrbid,'typeOfGeneratingProcess',           4) ! Ensemble Forecast
      ELSE
        IF (ptr_to_out%lanalysis) THEN
          CALL grib_set (igrbid,'typeOfGeneratingProcess',       202) ! Nudging Analysis
        ELSEIF (ptr_to_out%lsfc_ana) THEN
          CALL grib_set (igrbid,'typeOfGeneratingProcess',         0) ! External Analysis
        ELSE
          IF (ldfi) THEN
            CALL grib_set (igrbid,'typeOfGeneratingProcess',       1) ! Initialization
          ELSE
            IF (luseobs) THEN
              CALL grib_set (igrbid,'typeOfGeneratingProcess',   203) ! Nudgecast
            ELSE
              CALL grib_set (igrbid,'typeOfGeneratingProcess',     2) ! Forecast
            ENDIF
          ENDIF
        ENDIF
      ENDIF

    CASE DEFAULT

      IF (leps) THEN
        CALL grib_set (igrbid,'typeOfGeneratingProcess',           4) ! Ensemble Forecast
      ELSE
        IF (ptr_to_out%lanalysis .OR. ptr_to_out%lsfc_ana) THEN
          CALL grib_set (igrbid,'typeOfGeneratingProcess',         0) ! Analysis
        ELSE
          IF (ldfi) THEN
            CALL grib_set (igrbid,'typeOfGeneratingProcess',       1) ! Initialization
          ELSE
            CALL grib_set (igrbid,'typeOfGeneratingProcess',       2) ! Forecast
          ENDIF
        ENDIF
      ENDIF

    END SELECT

    ! We have to be very careful, in which order the different keys are set:
    ! to get the correct product definition template, this key has to be set
    ! before the shortname. It depends on the stepType, which is (still)
    ! determined by Grib1 timeRangeIndicator. 
    ! And even before, the stepType must be set
    ! To make it more complicated, the synthetic satellite images and the aerosol
    ! fields have different templates, which can be determined as follows:
    !   - the synthetic satellite images: by the GRIB1 leveltype 222
    !   - the aerosol values:             by the first four characters: 'AER_'

    ! First set the stepType
    CALL grib_set (igrbid,'stepType',                 yzsteptyp) ! 

    ! determine the correct productDefinitionTemplateNumber now
    SELECT CASE (yzsteptyp)
    CASE ('instant')

      SELECT CASE (var(n1,n2,n3)%name(1:3))
      ! different product definition templates have to be set for special variables
      ! (synthetic satellite images; aerosol variables (chemical constituents))
      ! Up to now these variables can be identified by their name

      CASE ('SYN')
        IF (leps) THEN
          iprdeftnr = 33     ! individual ensemble forecast
        ELSE
          iprdeftnr = 32     ! Standard products
        ENDIF
      CASE ('AER')
        ! atmospheric chemical constituents
        IF (leps) THEN
          iprdeftnr = 41     ! individual ensemble forecast
        ELSE
          iprdeftnr = 40     ! Standard products
        ENDIF
      CASE DEFAULT
        ! for all other variables
        IF (leps) THEN
          iprdeftnr =  1     ! individual ensemble forecast
        ELSE
          iprdeftnr =  0     ! Standard products
        ENDIF
      END SELECT

    CASE ('diff','min','max','avg','accum')

      ! statistically processed data
      IF (leps) THEN
        iprdeftnr =   11     ! individual ensemble forecast
      ELSE
        iprdeftnr =    8     ! Standard products
      ENDIF

    END SELECT

! here we had the productDefinitionTemplateNumber before
! but then we could not encode the new vertical coordinate in ensemble mode
! because of: NLEV: key/value not found error in grib_api
#endif
  END SELECT

!------------------------------------------------------------------------------
! Section 2: Level type information and vertical coordinate parameters
!------------------------------------------------------------------------------

  IF (idebug > 15) THEN
    WRITE (*,'(A)') '  io_metadata, make_grib: level type and vertical coordinate'
  ENDIF

  ! determine yzlevtyp and izlevtyp
  SELECT CASE (yextension)
  CASE ('s')
#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
    ! satellite channels on output
    ! leveltyp and value of level
    izlevtyp = 222_intgribf
    yzlevtyp = 'syntheticSatellites'
#endif
  CASE ('p')
    ! pressure level on output
    izlevtyp = 100_intgribf
    yzlevtyp = 'isobaricInhPa'
  CASE ('z')
    ! z- levels on output
    izlevtyp = 103_intgribf
    yzlevtyp = 'heightAboveSea'
  CASE DEFAULT
    ! model-levels on output
    izlevtyp = var(n1,n2,n3)%levtyp
    SELECT CASE (ptr_to_out%yform_write)
    CASE ('grb1','api1','bina')
      yzlevtyp = ylevltypes1(var(n1,n2,n3)%levtyp)
    CASE ('api2')
      yzlevtyp = ylevltypes2(var(n1,n2,n3)%levtyp)
    END SELECT
  END SELECT

  ! yzlevtyp has been set above, only W_SOxxx and FLake variables have to be adapted
  IF (izednr == 2) THEN
    IF (yzname(1:4) == 'W_SO') THEN
      yzlevtyp = 'depthBelowLandLayer'
    ENDIF

    ! FLake variables are treated different: their level information is only set
    ! by the shortName
    IF (izlen > 2) THEN
      IF (TRIM(yzname(izlen-2:izlen)) == '_LK') THEN
        yzlevtyp = 'FLAKE'
      ENDIF
    ENDIF
  ENDIF

  ! now set upper and lower boundary of the level
  SELECT CASE(izlevtyp)
  CASE (1,2,3,4,8,102,200)
    izlevtop = 0_intgribf
    izlevbot = 0_intgribf
  CASE (100)
    izlevtop = 0_intgribf
    izlevbot = NINT(slev*0.01_wp)
  CASE (103)
    izlevtop = 0_intgribf
    izlevbot = NINT(slev)
  CASE (105)
    izlevtop = 0_intgribf
    izlevbot = var(n1,n2,n3)%levbot
  CASE (109)
    izlevtop = 0_intgribf
    izlevbot = nlevel
  CASE (110,211)
    ! 211: also for snow levels in multi-layer snow model
    izlevtop = INT (nlevel  , intgribf)
    izlevbot = INT (nlevel+1, intgribf)
  CASE (111)
    IF ( (var(n1,n2,n3)%name == 'ALHFL_PL  ') .OR. &
         (var(n1,n2,n3)%name == 'W_SO      ') .OR. &
         (var(n1,n2,n3)%name == 'T_SO      ') .OR. &
         (var(n1,n2,n3)%name == 'W_SO_ICE  ') ) THEN
      ! multi-layer soil model
      izlevtop = 0_intgribf        ! levtop
      izlevbot = msoilgrib(nlevel) ! depth of main soil level in cm
    ELSE
      ! old 2 layer soil model
      izlevtop = 0_intgribf
      izlevbot = var(n1,n2,n3)%levbot
    ENDIF
  CASE (112)
    izlevtop = var(n1,n2,n3)%levtop
    izlevbot = var(n1,n2,n3)%levbot
#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
  CASE (222)
    izlevtop =   0_intgribf
    izlevbot = sat_compute(NINT(slev))%ngrib_chan(nlevel)
#endif
  CASE DEFAULT
    PRINT *, ' ERROR *** Wrong value for GRIB1 level type:    ', izlevtyp
    CALL model_abort (my_cart_id, 2015, 'wrong level type', 'make_grib')
  END SELECT 

  SELECT CASE (ptr_to_out%yform_write)

  CASE ('grb1','bina')

    ipds_out( 8) = izlevtyp
    ipds_out( 9) = izlevtop
    ipds_out(10) = izlevbot

#ifdef GRIBAPI
  CASE ('api1')

    CALL grib_set (igrbid,'indicatorOfTypeOfLevel',      izlevtyp)
    CALL grib_set (igrbid,'topLevel',                    izlevtop)
    CALL grib_set (igrbid,'bottomLevel',                 izlevbot)

  CASE ('api2')

    ! Now set shortName and productDefinitionTemplateNumber for non-generalVertical
!   IF (TRIM(yzlevtyp(1:15)) /= 'generalVertical' .AND.              &
!       TRIM(yzlevtyp)       /= 'syntheticSatellites') THEN
    IF (TRIM(yzlevtyp(1:15)) /= 'generalVertical') THEN
      CALL grib_set (igrbid, 'productDefinitionTemplateNumber', iprdeftnr)
!     CALL grib_set (igrbid, 'shortName',        TRIM(yzname))
    ELSE
      IF (TRIM(yzname) == 'HHL') THEN
        ! set shortname now, because otherwise no level information for hhl are written???
        CALL grib_set (igrbid, 'shortName',        TRIM(yzname))
      ENDIF
    ENDIF

    IF (yzlevtyp(1:15) == 'generalVertical') THEN
      CALL grib_set (igrbid, 'genVertHeightCoords',        1)
    ENDIF

    IF (TRIM(yzlevtyp) /= 'syntheticSatellites') THEN
      CALL grib_set (igrbid, 'typeOfLevel',      TRIM(yzlevtyp))
    ENDIF

    SELECT CASE (yzlevtyp)

    CASE ('surface')

      CALL grib_set (igrbid,'topLevel',                    izlevtop)
      CALL grib_set (igrbid,'bottomLevel',                 izlevbot)
      CALL grib_set_missing (igrbid, 'scaleFactorOfFirstFixedSurface')
      CALL grib_set_missing (igrbid, 'scaledValueOfFirstFixedSurface')

      IF     (TRIM(var(n1,n2,n3)%name) == 'HSURF') THEN
        ! some additional things have to be set for HSURF
        CALL grib_set (igrbid,'typeOfSecondFixedSurface',          101)
        CALL grib_set (igrbid,'scaleFactorOfSecondFixedSurface',     0)
        CALL grib_set (igrbid,'scaledValueOfSecondFixedSurface',     0)
!     ELSEIF (TRIM(var(n1,n2,n3)%name) == 'QV_S') THEN
!       CALL grib_set (igrbid, 'scaleFactorOfFirstFixedSurface',     0)
!       CALL grib_set_missing (igrbid, 'scaleFactorOfSecondFixedSurface')
!       CALL grib_set_missing (igrbid, 'scaledValueOfSecondFixedSurface')
      ELSE
        CALL grib_set_missing (igrbid, 'scaleFactorOfSecondFixedSurface')
        CALL grib_set_missing (igrbid, 'scaledValueOfSecondFixedSurface')
      ENDIF

      CALL grib_set (igrbid,'PVPresent', 0)
      CALL grib_set (igrbid,'NV',        0)

    CASE ('cloudBase','cloudTop','isothermZero')

      CALL grib_set (igrbid,'PVPresent', 0)
      CALL grib_set (igrbid,'NV',        0)

      CALL grib_set_missing (igrbid, 'scaleFactorOfFirstFixedSurface')
      CALL grib_set_missing (igrbid, 'scaledValueOfFirstFixedSurface')
      CALL grib_set (igrbid,'scaleFactorOfSecondFixedSurface',     0)
      CALL grib_set (igrbid,'scaledValueOfSecondFixedSurface',     0)

    CASE ('nominalTop')

      CALL grib_set (igrbid,'PVPresent', 0)
      CALL grib_set (igrbid,'NV',        0)

      CALL grib_set_missing (igrbid, 'scaleFactorOfFirstFixedSurface')
      CALL grib_set_missing (igrbid, 'scaledValueOfFirstFixedSurface')
      CALL grib_set_missing (igrbid, 'scaleFactorOfSecondFixedSurface')
      CALL grib_set_missing (igrbid, 'scaledValueOfSecondFixedSurface')

    CASE ('isobaricInhPa')          ! pressure levels

      CALL grib_set (igrbid,'level',               slev*0.01_wp)

      CALL grib_set (igrbid,'PVPresent', 0)
      CALL grib_set (igrbid,'NV',        0)

    CASE ('isobaricLayer')          ! pressure levels

      CALL grib_set (igrbid,'PVPresent', 0)
      CALL grib_set (igrbid,'NV',        0)

    CASE ('meanSea')                ! PMSL

      CALL grib_set (igrbid,'PVPresent', 0)
      CALL grib_set (igrbid,'NV',        0)

      CALL grib_set (igrbid,'scaleFactorOfFirstFixedSurface',     0)
      CALL grib_set (igrbid,'scaledValueOfFirstFixedSurface',     0)
      CALL grib_set_missing (igrbid, 'scaleFactorOfSecondFixedSurface')
      CALL grib_set_missing (igrbid, 'scaledValueOfSecondFixedSurface')

    CASE ('heightAboveSea')         ! z-levels

      CALL grib_set (igrbid,'level',                           slev)

      CALL grib_set (igrbid,'PVPresent', 0)
      CALL grib_set (igrbid,'NV',        0)

    CASE ('heightAboveGround')      ! T_2M, etc.

!!!      CALL grib_set (igrbid,'level',                     )

      CALL grib_set (igrbid,'PVPresent', 0)
      CALL grib_set (igrbid,'NV',        0)

      CALL grib_set_missing (igrbid, 'scaleFactorOfSecondFixedSurface')
      CALL grib_set_missing (igrbid, 'scaledValueOfSecondFixedSurface')

    CASE ('hybrid')

      CALL grib_set (igrbid,'typeOfFirstFixedSurface',           105)
      CALL grib_set (igrbid,'scaledValueOfFirstFixedSurface', nlevel)

      IF     (TRIM(var(n1,n2,n3)%name) == 'HHL') THEN
        ! some additional things have to be set for HHL
        CALL grib_set (igrbid,'typeOfSecondFixedSurface',          101)
      ELSE
        !US CALL grib_set (igrbid,'typeOfSecondFixedSurface',      255)
        CALL grib_set (igrbid,'typeOfSecondFixedSurface',            1)   ! such it is in converted data????
        CALL grib_set_missing (igrbid,'scaleFactorOfSecondFixedSurface')
        CALL grib_set_missing (igrbid,'scaledValueOfSecondFixedSurface')
      ENDIF

      ! vertical coordinate parameters
      ! vertical coordinate parameters
      CALL grib_set (igrbid,'PVPresent',           1)
      CALL grib_set (igrbid,'NV',        inrvert_out)
      CALL grib_set (igrbid,'pv',             pv_out)

    CASE ('hybridLayer')

      CALL grib_set (igrbid, 'topLevel',                   nlevel)
      CALL grib_set (igrbid, 'bottomLevel',                nlevel + 1)

      ! vertical coordinate parameters
      CALL grib_set (igrbid,'PVPresent',           1)
      CALL grib_set (igrbid,'NV',        inrvert_out)
      CALL grib_set (igrbid,'pv',             pv_out)

    CASE ('generalVertical')

      CALL grib_set (igrbid, 'typeOfFirstFixedSurface',           150)
      CALL grib_set (igrbid, 'scaleFactorOfFirstFixedSurface',      0)
      CALL grib_set (igrbid, 'scaledValueOfFirstFixedSurface', nlevel)
!     CALL grib_set (igrbid, 'level',                          nlevel)

      CALL grib_set (igrbid, 'NV',                                  6)

      ! Only for HHL
      IF (TRIM(var(n1,n2,n3)%name) == 'HHL') THEN
        CALL grib_set (igrbid, 'typeOfSecondFixedSurface',          101)

        mlev = NINT (vcoord%vert_coord(nlevel) *  100.0_wp, iintegers)   ! only 2 digits after decimal point
        CALL grib_set (igrbid, 'scaleFactorOfSecondFixedSurface',     2)  ! values in cm (more or less)
        CALL grib_set (igrbid, 'scaledValueOfSecondFixedSurface',  mlev)

!US     CALL grib_set_missing (igrbid, 'scaleFactorOfSecondFixedSurface')
!US     CALL grib_set_missing (igrbid, 'scaledValueOfSecondFixedSurface')
      ENDIF

      CALL grib_set (igrbid, 'nlev',                   vcoord%nlevels)
      CALL grib_set (igrbid, 'numberOfVGridUsed',      vcoord%ivctype)
      CALL grib_set (igrbid, 'uuidOfVGrid',            vcoord%vc_uuid)

    CASE ('generalVerticalLayer')

      CALL grib_set (igrbid, 'typeOfFirstFixedSurface',           150)
      CALL grib_set (igrbid, 'topLevel',                       nlevel)

      CALL grib_set (igrbid,'typeOfSecondFixedSurface',           150)
      CALL grib_set (igrbid, 'bottomLevel',                nlevel + 1)

      CALL grib_set (igrbid,'NV',                                   6)
!     CALL grib_set (igrbid,'numberOfVerticalGridDescriptors',      6)

      CALL grib_set (igrbid, 'nlev',                   vcoord%nlevels)
      CALL grib_set (igrbid, 'numberOfVGridUsed',      vcoord%ivctype)
      CALL grib_set (igrbid, 'uuidOfVGrid',            vcoord%vc_uuid)

    CASE ('depthBelowLand')

      SELECT CASE (yzname)

      CASE ('T_S','T_M','T_CL')

        CALL grib_set (igrbid, 'level',                0)

        CALL grib_set (igrbid, 'scaleFactorOfFirstFixedSurface',          0)
        CALL grib_set (igrbid, 'scaledValueOfFirstFixedSurface',   var(n1,n2,n3)%levbot)
        CALL grib_set_missing (igrbid, 'scaleFactorOfSecondFixedSurface')
        CALL grib_set_missing (igrbid, 'scaledValueOfSecondFixedSurface')

      CASE DEFAULT                                 !for T_SO (0...8)

        IF (nlevel == 0) THEN
          mlev = 0_iintegers
          CALL grib_set (igrbid, 'scaleFactorOfFirstFixedSurface',     0)  ! T_SO(0)
        ELSEIF (nlevel == 1) THEN
          mlev = NINT (czmls(nlevel) * 1000.0_wp, iintegers)
          CALL grib_set (igrbid, 'scaleFactorOfFirstFixedSurface',     3)  ! 0.005 m
        ELSE
          mlev = NINT (czmls(nlevel) *  100.0_wp, iintegers)
          CALL grib_set (igrbid, 'scaleFactorOfFirstFixedSurface',     2)  ! values in cm
        ENDIF

        CALL grib_set (igrbid, 'scaledValueOfFirstFixedSurface',    mlev)

        CALL grib_set_missing (igrbid, 'scaleFactorOfSecondFixedSurface')
        CALL grib_set_missing (igrbid, 'scaledValueOfSecondFixedSurface')

      END SELECT

      CALL grib_set (igrbid,'PVPresent', 0)
      CALL grib_set (igrbid,'NV',        0)

    CASE ('depthBelowLandLayer')

      IF ( (var(n1,n2,n3)%levbot == 0) .AND. (var(n1,n2,n3)%levtop == 0) .AND. &
                                                (nlevel /= 0)) THEN ! W_SO/W_SO_ICE

        mlev = NINT (czhls(nlevel-1) * 100.0_wp, iintegers)
        CALL grib_set (igrbid, 'scaleFactorOfFirstFixedSurface',       2)
        CALL grib_set (igrbid, 'scaledValueOfFirstFixedSurface',    mlev)

        mlev = NINT (czhls(nlevel  ) * 100.0_wp, iintegers)
        CALL grib_set (igrbid, 'scaleFactorOfSecondFixedSurface',      2)
        CALL grib_set (igrbid, 'scaledValueOfSecondFixedSurface',   mlev)

      ELSE         ! W_Gx, RUNOFF_x

        CALL grib_set (igrbid, 'scaleFactorOfFirstFixedSurface',       2)
        CALL grib_set (igrbid, 'scaledValueOfFirstFixedSurface',    var(n1,n2,n3)%levtop)
        CALL grib_set (igrbid, 'scaleFactorOfSecondFixedSurface',      2)
        CALL grib_set (igrbid, 'scaledValueOfSecondFixedSurface',   var(n1,n2,n3)%levbot)

      ENDIF

      CALL grib_set (igrbid,'PVPresent', 0)
      CALL grib_set (igrbid,'NV',        0)

    CASE ('entireAtmosphere')

      CALL grib_set (igrbid,'PVPresent', 0)
      CALL grib_set (igrbid,'NV',        0)

      CALL grib_set_missing (igrbid, 'scaleFactorOfFirstFixedSurface')
      CALL grib_set_missing (igrbid, 'scaledValueOfFirstFixedSurface')
      CALL grib_set_missing (igrbid, 'scaleFactorOfSecondFixedSurface')
      CALL grib_set_missing (igrbid, 'scaledValueOfSecondFixedSurface')

    CASE ('snowLayer')

      CALL grib_set (igrbid,'PVPresent', 0)
      CALL grib_set (igrbid,'NV',        0)

    CASE ('syntheticSatellites')
#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
      ! for synthetic satellite images the productDefinitionTemplateNumber
      ! has to be defined here, otherwise the CentralWaveNumber is not defined
!     CALL grib_set (igrbid,'productDefinitionTemplateNumber', 32) ! Standard products

      ! set vertical coordinate parameters off
      CALL grib_set (igrbid,'PVPresent', 0)
      CALL grib_set (igrbid,'NV',        0)

      ! Set even longer shortName
      ! Add label of BT/RAD, CL/CS
      SELECT CASE (sat_compute(NINT(slev))%ngrib_aees(nlevel))
      CASE (1)
        yzname = TRIM(yzname)//'_BT_CL_'
      CASE (2)
        yzname = TRIM(yzname)//'_BT_CS_'
      CASE (3)
        yzname = TRIM(yzname)//'_RAD_CL_'
      CASE (4)
        yzname = TRIM(yzname)//'_RAD_CS_'
      END SELECT

      ! Add name of channel
      yzname = TRIM(yzname)//sat_compute(NINT(slev))%ychan_name(sat_compute(NINT(slev))%ngrib_chan(nlevel))

      ! set scale factor (is not set by shortName.def
      CALL grib_set (igrbid, 'scaleFactorOfCentralWaveNumber',      0)
!     CALL grib_set (igrbid, 'shortName',        TRIM(yzname))
#endif

    CASE ('FLAKE')

      ! nothing has to be done

    CASE DEFAULT

      CALL grib_set (igrbid,'PVPresent', 0)
      CALL grib_set (igrbid,'NV',        0)

      CALL grib_set_missing (igrbid, 'scaleFactorOfFirstFixedSurface')
      CALL grib_set_missing (igrbid, 'scaledValueOfFirstFixedSurface')
      CALL grib_set_missing (igrbid, 'scaleFactorOfSecondFixedSurface')
      CALL grib_set_missing (igrbid, 'scaledValueOfSecondFixedSurface')

    END SELECT

    ! And now set the product definition template number for generalVertical
    ! this seems to work also for the ensemble mode
    IF (yzlevtyp(1:15) == 'generalVertical') THEN
      CALL grib_set (igrbid, 'productDefinitionTemplateNumber', iprdeftnr)
!     CALL grib_set (igrbid, 'shortName',        TRIM(yzname))
    ENDIF

    ! Set ensemble information
    IF (leps) THEN
      CALL grib_set (igrbid,'typeOfEnsembleForecast',             192)
      CALL grib_set (igrbid,'perturbationNumber',          iepsmem) ! Ensemble Forecast
      CALL grib_set (igrbid,'numberOfForecastsInEnsemble', iepstot) ! Ensemble Forecast

      ! Only for HHL: set these values again, because grib_api re-computes them
      !               when switching the template
      IF (TRIM(var(n1,n2,n3)%name) == 'HHL') THEN
        CALL grib_set (igrbid, 'typeOfSecondFixedSurface',          101)

        mlev = NINT (vcoord%vert_coord(nlevel) *  100.0_wp, iintegers)   ! only 2 digits after decimal point
        CALL grib_set (igrbid, 'scaleFactorOfSecondFixedSurface',     2)  ! values in cm (more or less)
        CALL grib_set (igrbid, 'scaledValueOfSecondFixedSurface',  mlev)
      ENDIF
    ENDIF

! should be ok with above settings
!   IF (TRIM(yzname) /= 'HHL') THEN
      ! for ensemble mode, it has to be set again
      CALL grib_set (igrbid,'shortName',        TRIM(yzname))
!   ENDIF
#endif

  END SELECT

!------------------------------------------------------------------------------
! Section 3: Time and Date
!     GRIB1: ipds_out(11,12,13,14,15,16,17,18,19,22)
!     grib_api: dataDate, dataTime, indicatorOfUnitOfTimeRange,
!               startStep,endStep,timeRangeIndicator
!------------------------------------------------------------------------------

  IF (idebug > 15) THEN
    WRITE (*,'(A)') '  io_metadata, make_grib: time and date'
  ENDIF

  ! Reference time of data
  ydatref(1:14) = ydate_ini(1:14)

  IF (ptr_to_out%lanalysis .OR. ptr_to_out%lsfc_ana) THEN
    IF ( (var(n1,n2,n3)%ntri <= 1) .OR. (ncenter == 78) ) THEN
      ydatref(1:14) = ydatact(1:14)
    ENDIF    
  ENDIF

  READ(ydatref,'(6I2)') jcgb, jjgb, mmgb, ndgb, nhgb, nmingb
  READ(ydatref( 1: 8),'(I8)') ndate
  READ(ydatref( 9:12),'(I4)') ntime
  READ(ydatref(13:14),'(I2)') nsecond

  IF(jjgb == 0) THEN
    ipds_out(11) = 100_intgribf
    ipds_out(22) = INT (jcgb, intgribf)
  ELSE  
    ipds_out(11) = INT (jjgb, intgribf)
    ipds_out(22) = INT (jcgb, intgribf) + 1
  ENDIF
  ipds_out(12) = INT (mmgb,   intgribf)
  ipds_out(13) = INT (ndgb,   intgribf)
  ipds_out(14) = INT (nhgb,   intgribf)
  ipds_out(15) = INT (nmingb, intgribf)
   
  ! Indicator of unit of time range
  ipds_out(16) = INT (ptr_to_out%nunit_of_time, intgribf)

  ! Time Range Indicator
  ! determine a local time range indicator
  nztri = var(n1,n2,n3)%ntri

  IF( ntstep == 0 ) THEN
    ! In step 0, all variables but the near-surface analysis fileds
    ! are written with either tri=0 or with tri=1 (if ldfi is used)
    IF (ptr_to_out%lsfc_ana) THEN
!US the ncenter==78 is local use, but this way it is also used
!   by other centers, therefore comment that out for the moment being
!   (also few lines below)
!US   IF ( (ptr_to_out%lanalysis) .AND. (ncenter == 78) ) THEN
      IF ( (ptr_to_out%lanalysis)  ) THEN
        nztri = 13_intgribf   ! for all products
      ELSE
        nztri =  0_intgribf
      ENDIF
    ELSE
      IF (ldfi) THEN
        nztri = 1_intgribf   ! this means initialized analyses. 
      ELSE
        nztri = 0_intgribf   ! this means uninitialized analyses. 
      ENDIF
    ENDIF
  ELSE
!US IF ( (ptr_to_out%lanalysis) .AND. (ncenter == 78) ) THEN
    IF ( (ptr_to_out%lanalysis) ) THEN
      nztri = 13_intgribf   ! for all products
    ENDIF
  ENDIF
   
  ! determine startStep and endStep
  ! Time p1 and p2:  
  !  p1:       actual time for which a forecast product is valid
  !            last deletion of min, max-values (tri=2)
  !       or:  begin of averaging period        (tri=3)
  !  p2:       actual time for which a forecast product is valid
  !       or:  end of averaging period          (tri=3)

  SELECT CASE(nztri)
  CASE(0,1)
    IF (ptr_to_out%lanalysis .OR. ptr_to_out%lsfc_ana) THEN
      istartstep = 0_intgribf
      iendstep = 0_intgribf
    ELSE
      istartstep = inrsteps (ntstep,   dt, ptr_to_out%nunit_of_time)
      iendstep = 0_intgribf
    ENDIF
  CASE(2)
    IF(TRIM(var(n1,n2,n3)%name) == 'VMAX_10M'  .OR. &
       TRIM(var(n1,n2,n3)%name) == 'VGUST_DYN' .OR. &
       TRIM(var(n1,n2,n3)%name) == 'VGUST_CON' .OR. &
       TRIM(var(n1,n2,n3)%name) == 'VABSMX_10M') THEN
      istartstep = inrsteps (nlastmxu, dt, ptr_to_out%nunit_of_time)
    ELSE   ! tmin_2m, tmax_2m
      istartstep = inrsteps (nlastmxt, dt, ptr_to_out%nunit_of_time)
    ENDIF
    iendstep = inrsteps (ntstep,   dt, ptr_to_out%nunit_of_time)
  CASE(3)
    IF (lbdclim) THEN
      istartstep = inrsteps (nlastout, dt, ptr_to_out%nunit_of_time)
    ELSE
      istartstep = 0_intgribf
    ENDIF
    iendstep = inrsteps (ntstep,   dt, ptr_to_out%nunit_of_time)
  CASE(4)
    IF (lbdclim) THEN
      istartstep = inrsteps (nlastout, dt, ptr_to_out%nunit_of_time)
    ELSE
      istartstep = 0_intgribf
    ENDIF
    iendstep = inrsteps (ntstep,   dt, ptr_to_out%nunit_of_time)
  CASE(13)
    ! assimilation runs at DWD
    istartstep = 0_intgribf
    iendstep = 0_intgribf

    IF (var(n1,n2,n3)%ntri >= 2) THEN
      iendstep = inrsteps (ntstep,   dt, ptr_to_out%nunit_of_time)
    ENDIF

    IF (var(n1,n2,n3)%name(1:8) == 'VMAX_10M'  .OR. &
        var(n1,n2,n3)%name(1:9) == 'VGUST_DYN' .OR. &
        var(n1,n2,n3)%name(1:9) == 'VGUST_CON' .OR. &
       TRIM(var(n1,n2,n3)%name) == 'VABSMX_10M') THEN
      iendstep = iendstep- inrsteps (nlastmxu, dt, ptr_to_out%nunit_of_time)
    ENDIF
    IF ( (var(n1,n2,n3)%name(1:7) == 'TMIN_2M')   .OR.     &
         (var(n1,n2,n3)%name(1:7) == 'TMAX_2M') ) THEN
      iendstep = iendstep- inrsteps (nlastmxt, dt, ptr_to_out%nunit_of_time)
    ENDIF
  END SELECT

  ! Check, if ipds_out(17,18) are between 0 and 254
  ! (but for restart-files in climate mode this will be violated, 
  !  therefore omit the warnings; the same might be true for very
  !  high resolution artificial cases with a high output frequency)
  IF (.NOT. lbdclim .AND. .NOT. lartif_data .AND.                      &
      (ptr_to_out%yform_write == 'grb1' .OR. ptr_to_out%yform_write == 'api1')) THEN
    IF (.NOT. lrestart) THEN
      IF ( (istartstep < 0_intgribf) .OR. (istartstep > 254_intgribf) ) THEN
        PRINT *, ' WARNING: *** istartstep = ', istartstep,      &
                 ' is outside the valid range *** ', n1, n2, n3, ptr_to_out%nunit_of_time, nztri
      ENDIF
      IF ( (iendstep < 0_intgribf) .OR. (iendstep > 254_intgribf) ) THEN
        PRINT *, ' WARNING: *** iendstep = ', iendstep,      &
                 ' is outside the valid range *** ', n1, n2, n3, ptr_to_out%nunit_of_time, nztri
      ENDIF
    ENDIF
  ENDIF

  ipds_out(17) = istartstep
  ipds_out(18) = iendstep
  ipds_out(19) = nztri

#ifdef GRIBAPI
  SELECT CASE (ptr_to_out%yform_write)

  CASE ('api1','api2')

    CALL grib_set (igrbid,'dataDate',                     ndate)     ! yyyymmdd
    CALL grib_set (igrbid,'dataTime',                     ntime)     ! hhmm
    CALL grib_set (igrbid,'second',                     nsecond)     ! ss

    ! must be set before startstep, endstep!
    IF     (ptr_to_out%yform_write == 'api1') THEN
      CALL grib_set (igrbid,'timeRangeIndicator',           nztri)     !
    ENDIF

    SELECT CASE (yzsteptyp)

    CASE ('instant')

!print *, 'setting steps:     ', yzname, '  ', yzsteptyp, '  ', yzlevtyp, '  ', &
!         ptr_to_out%nunit_of_time, istartstep, iendstep
      CALL grib_set (igrbid,'startStep',             istartstep)     !
      IF (ptr_to_out%lsfc_ana .AND. TRIM(yzname) == 'TOT_PREC') THEN
        CALL grib_set (igrbid,'endStep',             iendstep)     !
      ENDIF

    CASE ('min','max','avg','accum')

!print *, 'setting statproc:  ', yzname, '  ', yzsteptyp, '  ', yzlevtyp, '  ', &
!         ptr_to_out%nunit_of_time, istartstep, iendstep
      CALL grib_set (igrbid,'indicatorOfUnitOfTimeRange', ptr_to_out%nunit_of_time)
      CALL grib_set (igrbid,'stepUnits',                  ptr_to_out%nunit_of_time)
      CALL grib_set (igrbid,'startStep',             istartstep)     !
      CALL grib_set (igrbid,'endStep',                 iendstep)     !

      IF (ptr_to_out%yform_write == 'api2') THEN
        CALL grib_set (igrbid,'typeOfTimeIncrement',                    2)                  !
        CALL grib_set (igrbid,'indicatorOfUnitForTimeRange',  ptr_to_out%nunit_of_time)     !
        CALL grib_set (igrbid,'numberOfMissingInStatisticalProcess',    0)
      ENDIF

    END SELECT

    ! must be called after setting startStep, endStep
    ! otherwise it is set to default 1 again
    CALL grib_set (igrbid,'indicatorOfUnitOfTimeRange', ptr_to_out%nunit_of_time)

  END SELECT
#endif

!------------------------------------------------------------------------------
! Section 4: local use information
!     GRIB1: ipds_out(37,41,42-46,47,50-52)
!     grib_api2: local use section
!------------------------------------------------------------------------------

  IF (idebug > 15) THEN
    WRITE (*,'(A)') '  io_metadata, make_grib: local information'
  ENDIF

  ipds_out(20) = 0_intgribf
  ipds_out(21) = 0_intgribf
! ipds_out(23) = 0_intgribf
  ipds_out(24) = 0_intgribf
   
  ! local use area
  IF (leps) THEN
    idims_out(11) =  54_intgribf
    ipds_out( 1)      =  66_intgribf
    ipds_out(37)      = 253_intgribf
    ipds_out(48:49)   =   0_intgribf
    ipds_out(53:54)   =   0_intgribf
  ELSE
    idims_out(11) =  47_intgribf
    ipds_out(37)      = 254_intgribf
  ENDIF
  ipds_out(38:40) = 0_intgribf

  IF ( yextension == 's') THEN
#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
    ipds_out(41) = sat_compute(NINT(slev))%ngrib_aees(nlevel)
#endif
  ELSE
    SELECT CASE(nsma_stat)
      CASE(21,23,27)
        ipds_out(41)    = INT (nsma_stat, intgribf)
      CASE DEFAULT
        SELECT CASE (ptr_to_out%yform_write)
        CASE ('api1','api2')
#ifdef GRIBAPI
          !DL localInformationNumber may be already defined via 
          !   setting of shortName, DO NOT override
          CALL grib_get (igrbid, 'localInformationNumber', ilin)
          ipds_out(41)    = INT (ilin, intgribf)
#endif
        CASE DEFAULT
          ! necessary for example for writing restarts
          ipds_out(41)    = 0_intgribf
        END SELECT
    END SELECT
  ENDIF

  ! write actual date and time to local use area of pds
  CALL DATE_AND_TIME(ydate,ytime)
  READ(ydate,'(I4,2I2)') ijj, imm, idd
  READ(ytime,'(3I2,4X)') ihh, imi, iss
  ! According to DWD standard, ipds_out(42) contains the offset to the year 1900
  ipds_out(42)    = INT (ijj-1900, intgribf)
  ipds_out(43)    = INT (imm, intgribf)
  ipds_out(44)    = INT (idd, intgribf)
  ipds_out(45)    = INT (ihh, intgribf)
  ipds_out(46)    = INT (imi, intgribf)

  ipds_out(47)    = INT (nvers, intgribf)

  IF (leps) THEN
    ipds_out(50)   = INT (iepstyp, intgribf)
    ipds_out(51)   = INT (iepstot, intgribf)
    ipds_out(52)   = INT (iepsmem, intgribf)
  ENDIF

#ifdef GRIBAPI
  SELECT CASE (ptr_to_out%yform_write)

  CASE ('api1')

    SELECT CASE (ncenter)
    CASE (78)  ! DWD
      CALL grib_set (igrbid,'localDefinitionNumber',      ipds_out(37))
      CALL grib_set (igrbid,'localInformationNumber',     ipds_out(41))
 
      CALL grib_set (igrbid,'localDecodeDateYear',        ipds_out(42))
      CALL grib_set (igrbid,'localDecodeDateMonth',       ipds_out(43))
      CALL grib_set (igrbid,'localDecodeDateDay',         ipds_out(44))
      CALL grib_set (igrbid,'localDecodeDateHour',        ipds_out(45))
      CALL grib_set (igrbid,'localDecodeDateMinute',      ipds_out(46))
      CALL grib_set (igrbid,'localVersionNumber',         ipds_out(47))

      IF (leps) THEN
        CALL grib_set (igrbid,'etyp',                     ipds_out(50))
        CALL grib_set (igrbid,'etot',                     ipds_out(51))
        CALL grib_set (igrbid,'enum',                     ipds_out(52))
      ENDIF

    CASE DEFAULT

      CALL grib_set (igrbid,'localInformationNumber',     ipds_out(41))
      CALL grib_set (igrbid,'localVersionNumber',         ipds_out(47))

    END SELECT

  CASE ('api2')

!DL SELECT CASE (ncenter)
!DL CASE (78)  ! DWD
    SELECT CASE (nactlocdefnr)
    CASE (252,253,254)

      CALL grib_set (igrbid,'localCreationDateYear',            ijj)
      CALL grib_set (igrbid,'localCreationDateMonth',           imm)
      CALL grib_set (igrbid,'localCreationDateDay',             idd)
      CALL grib_set (igrbid,'localCreationDateHour',            ihh)
      CALL grib_set (igrbid,'localCreationDateMinute',          imi)
      CALL grib_set (igrbid,'localCreationDateSecond',          iss)

      ! the value of ipds_out(41) is set above depending on the fields 
      ! which really need it, so here it can be just set to the GRIB2 key
      CALL grib_set (igrbid,'localInformationNumber',     ipds_out(41))

    CASE DEFAULT

      CALL grib_set (igrbid,'localInformationNumber',     ipds_out(41))

    END SELECT


  END SELECT
#endif

!------------------------------------------------------------------------------
! Section 5: Final coding of other sections: bms, bds
!------------------------------------------------------------------------------

  SELECT CASE (ptr_to_out%yform_write)

  CASE ('grb1')

    ! set the namelist-dependent idims_out values
    idims_out( 7) = ptr_to_out%ie_out_tot * ptr_to_out%je_out_tot
    idims_out(15) = ptr_to_out%ie_out_tot * ptr_to_out%je_out_tot
    idims_out(17) = ptr_to_out%ie_out_tot * ptr_to_out%je_out_tot

    ! gridpoints, simple packing, floating point data
    ibds(2)   = 0

    ! nrbit, number of bits
    ibds(5)   = ptr_to_out%nrbit

    ! no bitmap
    ibms(3)   = -2

  CASE ('api1')

#ifdef GRIBAPI
    ! Set these values again explicit, because they could be overwritten
    ! by input data (reduce size of ibmap)
    ! ibds (2)     = 0
    CALL grib_set (igrbid, 'bitmapPresent',      0)
    CALL grib_set (igrbid, 'bitsPerValue',   ptr_to_out%nrbit)
#endif

  CASE ('api2')

#ifdef GRIBAPI
    ! Set these values again explicit, because they could be overwritten
    ! by input data (reduce size of ibmap)
    ! ibds (2)     = 0
    CALL grib_set (igrbid, 'bitmapPresent',      0)

    IF ((TRIM(var(n1,n2,n3)%name) == 'HHL') .OR. (TRIM(var(n1,n2,n3)%name) == 'P')) THEN
      CALL grib_set (igrbid, 'bitsPerValue',         nrbit_phhl)
    ELSE
      CALL grib_set (igrbid, 'bitsPerValue',   ptr_to_out%nrbit)
    ENDIF
#endif

  END SELECT

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE make_grib_product

!==============================================================================
!==============================================================================
!+ Module procedure in io_metadata to create the grid definition block
!------------------------------------------------------------------------------

SUBROUTINE makegds

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure creates the grid definition block according to the
!   description of putgd1 (dwd) and the WMO. This is done for GRIB1 and the
!   resulting field igds_out can be used for DWD GRIB library and for grib_api
!
! Method:
!
!==============================================================================

! Subroutine arguments

! Local scalars:
INTEGER  (KIND=iintegers)          :: k
REAL     (KIND=wp)                 :: pollon_grib

#ifdef GRIBDWD
INTEGER (KIND=intgribf), EXTERNAL  :: IREFTS
#endif

CHARACTER(LEN=100) :: yzerrmsg

!- End of header
!==============================================================================

igds_out(:) = -999999

! The following numbers can be specified for DWD GRIB library, where the
! vertical coordinate parameters can be put in all GRIB messages. With grib_api
! these parameters are only put in messages with level type 109/110 and have to
! be specified accordingly.

! length of igds_out in bytes
igds_out(1) = 42 + (inrvert_out) * 4

IF ( ke > UBOUND( igds_out, DIM=1 ) - 37 ) THEN
  ! Some people like to test with several hundred levels. Then the dimensions
  ! of the DWD grib library are not big enough and have to be adapted
  ! Grib output with DWD Grib library is then not possible any more, but
  ! NetCDF is
  yzerrmsg = "ke > ngds-37: Choose a bigger value for ngds in data_io!"
  CALL model_abort (my_cart_id, 1234, yzerrmsg, "makegds")
END IF

igds_out(2) = inrvert_out

! location of the list of vertical coordinate parameters in bytes
igds_out(3) = 43

! data representation type
igds_out(4) = 10

! calculation of the left bottom corner and
! left bottom and right upper corner in millidegrees:
! This depends on the variable (U- and V-related variables are shifted in
! the Arakawa C-grid) and is determined during the output step.

igds_out( 9) = 8       ! this was 0 before; but due to Grib 1, Code Table 7
                       ! bit no. 5 must be set to 1, to indicate that u,v
                       ! components are relative to defined (rotated) grid

! increments
igds_out(12) = 0
igds_out(13) = 0


igds_out(14) = 64
igds_out(15:19) = 0

! coordinates of the pole: in the grib code the southern pole has to be
! specified: for the latitude this is the negative value, for the
! longitude it is the value + 180.0 and then limiting the value again
! to -180.0 ... 180.0
pollon_grib = pollon + 180.0_wp
IF (pollon_grib > 180.0_wp) THEN
  pollon_grib = pollon_grib - 360.0_wp
ENDIF

igds_out(20) = NINT(-pollat      * 1000.0_wp)
igds_out(21) = NINT( pollon_grib * 1000.0_wp)

#ifdef GRIBDWD
! in case of restarts, this value is hopefully not needed
! at least it is not checked
igds_out(22) = IREFTS(REAL(polgam, irealgrib))

! vertical coordinate parameters

! In case of restarts and not using DWD Grib-library, these values
! will not be written to the Restart Metadata, but in the first records
!  of the restart files.
! And they are also not checked in check_input_grid

DO k = 1, inrvert_out
  igds_out(25 + k) = IREFTS( REAL(pv_out(k), irealgrib) )
ENDDO
#endif

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE makegds

!==============================================================================
!==============================================================================

INTEGER (KIND=intgribf) FUNCTION inrsteps (nstep, dt, nuot)

!------------------------------------------------------------------------------
!
! Description:
!   This function computes the number of steps in units of nuot, necessary
!   to reach the forecast time given by nstep*dt. These are the steps 
!   that have to be specified for GRIB meta data
!
!------------------------------------------------------------------------------

! Function arguments
INTEGER (KIND=iintegers), INTENT(IN) ::   &
  nstep,     & ! time steps done during the simulation
  nuot         ! unit_of_time (GRIB indicatorOfUnitOfTimeRange)

REAL (KIND=wp),           INTENT(IN) ::   &
  dt           ! model time step

!------------------------------------------------------------------------------

REAL (KIND=wp)                       ::   hm

!------------------------------------------------------------------------------

  SELECT CASE (nuot)
  CASE ( 0);       hm  =    60.0_wp    !  1 minute
  CASE ( 1);       hm  =  3600.0_wp    !  1 hour
  CASE ( 2);       hm  = 86400.0_wp    !  1 day
  CASE (10);       hm  = 10800.0_wp    !  3 hours
  CASE (11);       hm  = 21600.0_wp    !  6 hours
  CASE (12);       hm  = 43200.0_wp    ! 12 hours
  CASE (13);       hm  =   900.0_wp    ! 15 minutes
  CASE (14);       hm  =  1800.0_wp    ! 30 minutes
  CASE (15);       hm  =   600.0_wp    ! 10 minutes
  CASE DEFAULT
    inrsteps = -1
    RETURN
  END SELECT

  inrsteps = NINT (nstep * dt / hm, intgribf)

END FUNCTION inrsteps

!==============================================================================
!==============================================================================

END MODULE io_metadata
