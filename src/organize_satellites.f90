!+ External procedure for organizing the computations of satellite images
!  and observation processsing
!------------------------------------------------------------------------------

SUBROUTINE organize_satellites (yaction, ierror, yerror)

!------------------------------------------------------------------------------
!
! Description:
!   This procedure is the driving routine for calling the computation of 
!   satellite radiances and brightness temperatures for synthetic satellite
!   images and data assimilations. It is called first at the beginning of the
!   program just to read the Namelist Input. During the time stepping it is
!   called at each step when output for the satellite images is required.
!
! Method:
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 3.7        2004/02/18 Ulrich Schaettler
!  Initial release
! 3.8        2004/03/23 Ulrich Schaettler
!  Changed communicator icomm_cart to icomm_world for the Namelist distribution
! 3.10       2004/06/23 Ulrich Schaettler
!  Added status of synthetic satellite images to transfer to output routines
! 3.13       2004/12/03 Christian Keil, Thorsten Reinhardt
!  Corrected computation of total cloud water content (syclw)
! 3.14       2005/01/25 Thorsten Reinhardt, Christian Keil
!  Use of qi+qs for computation of LMSynSat
! 3.16       2005/07/22 Davide Cesari
!  Use of qs only, if it is allocated
! 3.19       2006/04/25 Ulrich Schaettler
!  If an error occurs during the computations, the model will abort
! 3.21       2006/12/04 Ulrich Schaettler
!  Bug correction for computing the satellite zenith angle
!  Use timelevel nnow for the computations (is done also now in the output)
! V4_3         2008/02/25 Ulrich Schaettler
!  Correction for computing nsat_next in case of restart runs
! V4_5         2008/09/10 Ulrich Schaettler
!  Changed call to org_sat_tbs to call RTTOV library for a vector of profiles
! V4_8         2009/02/16 Guenther Zaengl
!  Add p0hl (reference pressure on half levels) for full consistency with new
!  implementation of reference atmosphere
! V4_9         2009/07/16 Ulrich Schaettler
!  Included call to RTTVI to init-routine; allocate additional memory for that
!  Set reference ozone profile
!  Adapt computational indices to vector version, if istartpar==iendpar
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_18        2011/05/26 Robin Faulwetter, Ulrich Schaettler
!  Implemented usage of RTTOV9
! V4_20        2011/08/31 Oliver Fuhrer
!  Allocation of nsat_steps in any case, because it is accessed in the 
!  init- and compute-phase of organize_satellites
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer
!  Replaced qx-variables by using them from the tracer module
! V4_26        2012/12/06 Andreas Messer
!  Modifications for using also RTTOV10 and for satellite observation processing
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  MESSy interface introduced (AK)
!  h_ice_dummy, qg_dummy  consistently nullified and deallocated (AK)
!  Introduced conditional compilation for Nudging, to separate nudging parts 
!  (lobsrad) from SYNSAT parts (US)
! V4_28        2013/07/12 Ulrich Schaettler
!  Extensions to provide GRIB2 shortnames and additional meta data
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure: put allocation
!   and deallocation of memory for synsat variables to global memory
! V5_1         2014-11-28 Ulrich Blahak, Oliver Fuhrer
!  Changed the format of some YUSPECIF entries for the CLM namelist tool.
!  Implemented F2003 IOMSG-mechanism for better namelist error messages. (UB)
!  Replaced ireals by wp (working precision) (OF)
! V5_2         2015-05-21 Robin Faulwetter, Ulrich Schaettler
!  3 new NL variables in /SATCTL/: lread_ct, yclouddir, linterp
!  Deallocate additional fields at the end of the program (US)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

USE data_parameters,    ONLY:  wp, iintegers

!------------------------------------------------------------------------------

USE data_fields,      ONLY :              &
    pp, p0, p0hl, ps,                     &
    t, t_g, t_2m, qv_2m, u_10m, v_10m,    &
    fr_land, clw_con, clc_con, clc_sgs,   &
    h_ice, synme7, synmsg, rlat, rlon

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :  &
    ie, je, ke, istartpar, iendpar, jstartpar, jendpar, degrad, raddeg, &
    lalloc_h_ice, idt_qv, idt_qc, idt_qi, idt_qg, idt_qs, dt

!------------------------------------------------------------------------------

USE data_runcontrol,  ONLY :              &
    luse_rttov, itype_gscp,               &
    nstart, ntstep, nnow, nuspecif,       &
    isynsat_stat,                         &
    ldebug_dia, lprintdeb_all, idbg_level,&
    itype_calendar

!------------------------------------------------------------------------------

USE data_parallel,    ONLY :  &
    nproc,           & ! total number of processors: nprocx * nprocy
    num_compute,     & ! number of compute PEs
    my_world_id,     & ! rank of this subdomain in the global communicator
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    icomm_cart,      & ! communicator for the virtual cartesian topology
    icomm_world,     & ! communicator that belongs to igroup_world, i.e.
                       ! = MPI_COMM_WORLD
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    imp_integers,    & ! determines the correct INTEGER type used in the
                       ! model for MPI
    imp_character,   & ! determines the correct CHARACTER type used in the
                       ! model for MPI
    imp_logical,     & ! determines the correct LOGICAL   type used in the
                       ! model for MPI
    nboundlines,     & !
    intbuf, realbuf, logbuf, charbuf ! Buffer for distributing the namelists

!------------------------------------------------------------------------------

USE data_io,          ONLY :  &
    root, pp_nl, noutst, ydate_ini

!------------------------------------------------------------------------------

USE data_satellites , ONLY :  &
    sat_input_type, sat_compute, yrttov_satell_table, yrttov_sensor_table,   &
    num_sensors, numchans, maxknchpf, kiu1,                                  &
    nsat_steps, nsat_next, jpch, jpnsat, jplev, jppf, jpnav, jpnsav, jpnssv, &
    utmx, utmn, uqmx, uqmn, uomx, uomn, ppres_d, o3_ref, ivch,               &
    const_aheat, lcon_clw, lsynsat, lobsrad,                                 &
    instruments, channels, mchans, n_chans, addclouds,                       &
    extrp_type, itype_rttov, iceshape, iwc2effdiam, lread_ct, yclouddir,     &
    linterp

!------------------------------------------------------------------------------

USE parallel_utilities, ONLY :  &
    distribute_values, global_values

!------------------------------------------------------------------------------

USE src_sat_tbs,      ONLY :  org_sat_tbs

!------------------------------------------------------------------------------

USE src_obs_rad,      ONLY :  &
    nsensors, sensors

#ifdef NUDGING
USE src_obs_rad,      ONLY :  &
    input_obs_radctl, input_obs_satpp, calc_obs_satpp, cloud_read_prep
#endif

!------------------------------------------------------------------------------

USE environment,      ONLY :  model_abort, get_free_unit

!------------------------------------------------------------------------------

USE utilities,        ONLY :  get_utc_date

!------------------------------------------------------------------------------

USE src_tracer,       ONLY :  trcr_get, trcr_errorstr
USE data_tracer,      ONLY :  T_ERR_NOTFOUND

!------------------------------------------------------------------------------

#if   defined RTTOV9
USE src_sat_rttov,    ONLY :  org_sat_rttov
USE mo_rttov_ifc,     ONLY :         &
    rttov_init, rttov_cleanup, NO_ERROR, rttov_ifc_errMsg
#elif defined RTTOV10
USE src_sat_rttov,    ONLY :  org_sat_rttov
USE mo_rttov_ifc,     ONLY :         &
    rttov_init, rttov_cleanup, NO_ERROR, rttov_ifc_errMsg, rttov_ifc_version
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================

#ifdef RTTOV7
INCLUDE "RTTOV7_RTTVI_new.interface"
#endif

!==============================================================================
!
! Parameter list:
! ---------------

CHARACTER (LEN= *),       INTENT(IN)            ::                      &
  yaction      ! action to be performed

INTEGER (KIND=iintegers), INTENT(OUT)           ::                      &
  ierror       ! error status

CHARACTER (LEN= *),       INTENT(OUT)           ::                      &
  yerror       ! error message

! Local parameters: 
#if defined RTTOV9
  INTEGER, PARAMETER :: rttov_ifc_version = 9 ! For backwards compat with old rttov ifc
#endif  

! Local variables: 
! ----------------

CHARACTER(LEN=*), PARAMETER :: &
  yusyn = 'YUSYNSAT'                ! synsat diagnostic file

! General variables
INTEGER (KIND=iintegers)   ::        &
  izerrstat, izstatus,               & ! Error status
  izdebug,                           & ! for debug messages
  i, j, k,                           & ! Loop variables
  izdim, jzdim,                      & ! Horizontal dimensions
  izdimstart, izdimend,              & ! Zonal start/end for processor
  jzdimstart, jzdimend                 ! Meridional start/end for processor

! Namelist INPUT file
INTEGER (KIND=iintegers)   :: nuin
CHARACTER (LEN=16)         :: yinput
CHARACTER (LEN=80)         :: yzerror
CHARACTER (LEN=25)         :: yzroutine

! Fields for RTTOV7
REAL (KIND=wp)             ::   &
  syclc(ie,je,ke), syclw(ie,je,ke), zfield(ie,je,ke)

! Pointers to fields, that are possibly not available
REAL(KIND=wp),     POINTER ::        &
  h_ice_dummy(:,:,:)  => NULL(),     &
  qg_dummy   (:,:,:)  => NULL()

! Interfaces to RTTVI / RTTOVCLD
INTEGER (KIND=iintegers), SAVE ::           &
  nprof,     & ! no. of profiles processed in parallel
  kppf,      & ! max no. profiles processed in parallel
  kpnsat,    & ! max no. of satellites
  kplev,     & ! no of rt levels
  kpch,      & ! max no. of channels
  kpchus,    & ! max no. of channels used
  kpnav,     & ! max no of profile variables
  kpnsav,    & ! max no of surface variables
  kpnssv,    & ! max no of skin variables
  kpncv        ! max no of cloud variables

! Tracer pointers
REAL (KIND=wp),     POINTER :: &
  qv  (:,:,:) => NULL(),    &   ! QV at nnow
  qc  (:,:,:) => NULL(),    &   ! QC at nnow
  qg  (:,:,:) => NULL(),    &   ! QG at nnow
  qi  (:,:,:) => NULL(),    &   ! QI at nnow
  qs  (:,:,:) => NULL()         ! QS at nnow


!for reading NWCSAF SEVIRI cloud type and cloud top height
CHARACTER (LEN=14)         :: yzactdate
CHARACTER (LEN=28)         :: yzdum  ! actual date in the form   wd   dd.mm.yy  hh UTC

REAL (KIND=wp)             ::     &
  rzdum,                          &       ! dummy variable
  ctype(ie,je),                   &       ! cloud type
  cth(ie,je)

INTEGER (KIND=iintegers)   ::  izdum                                ! dummy variable

INTEGER (KIND=iintegers), SAVE :: nusyn

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------
!- Begin Subroutine organize_satellites
!------------------------------------------------------------------------------

! Parameterlist
ierror    = 0_iintegers
yerror    = '      '
yzroutine = 'organize_satellites'

! Local variables
izerrstat = 0_iintegers
izstatus  = 0_iintegers
kiu1(:)   = 0_iintegers
yzerror   = '      '

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

IF (yaction == 'input') THEN

  !----------------------------------------------------------------------------
  ! Section 1: Input of the Namelist
  !----------------------------------------------------------------------------

  ! Open NAMELIST-INPUT file
  IF (my_world_id == 0) THEN
    PRINT *,'    INPUT OF THE NAMELISTS FOR RTTOV SYNSAT'
    yinput   = 'INPUT_SAT'
    nuin     =  1

    OPEN(nuin   , FILE=yinput  , FORM=  'FORMATTED', STATUS='UNKNOWN',  &
         IOSTAT=izerrstat)
    IF(izerrstat /= 0) THEN
      yerror = ' ERROR    *** Error while opening file INPUT_SAT *** '
      ierror = 9001
      RETURN
    ENDIF
  ENDIF

  CALL input_satctl  (nuspecif, nuin, izerrstat)

  IF (my_world_id == 0) THEN
    ! Close file for input of the NAMELISTS
    CLOSE (nuin    , STATUS='KEEP')
  ENDIF

  IF (izerrstat < 0) THEN
    yerror = 'ERROR *** while reading NAMELIST Group /SATCTL/ ***'
    ierror = 9002
  ELSEIF (izerrstat > 0) THEN
    yerror = ' ERROR *** Wrong values occured in NAMELIST INPUT_SAT ***'
    ierror = 9003
  ENDIF

#ifdef NUDGING
  IF (lobsrad) THEN
    IF (my_world_id == 0) THEN
      PRINT *,'    INPUT OF THE NAMELISTS FOR RTTOV RAD OBS'
      yinput   = 'INPUT_OBS_RAD'
      nuin     =  1

      OPEN(nuin   , FILE=yinput  , FORM=  'FORMATTED', STATUS='UNKNOWN',  &
           IOSTAT=izerrstat)
      IF(izerrstat /= 0) THEN
        yerror = ' ERROR    *** Error while opening file INPUT_OBS_RAD *** '
        ierror = 9004
        RETURN
      ENDIF
    ENDIF

    CALL input_obs_radctl ( nuspecif, nuin, izerrstat)

    IF (my_world_id == 0) THEN
      ! Close file for input of the NAMELISTS
      CLOSE (nuin    , STATUS='KEEP')
    ENDIF

    IF (izerrstat < 0) THEN
      yerror = 'ERROR *** while reading NAMELIST group /TOVS_OBS/ ***'
      ierror = 9005
    ELSEIF (izerrstat > 0) THEN
      yerror = ' ERROR *** Wrong values occured in NAMELIST INPUT_OBS_RAD ***'
      ierror = 9006
    ENDIF
  ENDIF
#endif

ELSEIF (yaction == 'input-satpp') THEN

  !----------------------------------------------------------------------------
  ! Section 2: Reading of satpp files
  !----------------------------------------------------------------------------

#ifdef NUDGING
  IF (lobsrad) THEN
    CALL input_obs_satpp(izerrstat)
    
    IF (izerrstat /= 0) THEN
      PRINT *, ' ERROR *** while reading satpp files *** '
      ierror = 9010
      RETURN
    ENDIF

    ! Reading of NWCSAF SEVIRI cloud type
    IF (lread_ct) THEN
      CALL get_utc_date ( ntstep, ydate_ini, dt, itype_calendar, yzactdate       &
                        , yzdum, izdum,rzdum)
!     =================

      CALL cloud_read_prep(yzactdate,yaction,ctype,cth,izerrstat)
      !cloud classification: 
        ! 0     non-processed
        ! 1     cloud-free land
        ! 2     cloud free sea
        ! 3     land contaminated by snow
        ! 4     sea contaminated by snow/ice   
        ! 5,6   very low clouds
        ! 7,8   low clouds
        ! 9,10  medium clouds
        ! 11,12 high opaque clouds
        ! 13,14 very high opaque clouds
        ! 15    high semitransparent thin clouds
        ! 16    high semitransparent meanly thick clouds
        ! 17    high semitransparent thick clouds
        ! 18    high semitransparent above low or medium clouds
        ! 19    fractional clouds (sub-pixel water clouds)
        ! 20    undefined   
      IF (izerrstat /= 0) THEN
        PRINT *, ' ERROR *** while reading cloud type GRIB file *** '
        ierror = 9100
        RETURN
      ENDIF
    ENDIF

  ENDIF  !lobsrad

#endif

ELSEIF (yaction == 'init') THEN

  !----------------------------------------------------------------------------
  ! Section 3: Initialization of the RTTOV library
  !----------------------------------------------------------------------------

  ! set status to 0 initially; 
  ! isynsat_stat is checked by the output routine whether output can be done
  ! (isynsat_stat = 0), or whether something has happened and no images are
  ! available for output (isynsat_stat /= 0)
  isynsat_stat = 0

  ! initialize the counter for the output steps
  ! (take into account restart runs)
  IF (nstart == 0) THEN
    nsat_next = 1
  ELSE
    nsat_next = 1
    DO WHILE ( (nsat_steps(nsat_next) <  nstart) .AND. (nsat_steps(nsat_next) >= 0) )
      nsat_next = nsat_next + 1
    ENDDO
  ENDIF

  ! Set number of channels for all satellites
  numchans  (1) = sat_compute(1)%num_chan
  numchans  (2) = sat_compute(2)%num_chan

  IF (itype_rttov == 7) THEN

#ifdef RTTOV7
    ! initialize the MPI-variables in the library
    CALL rttov7_mpi_settings (nproc, my_cart_id, icomm_cart, imp_reals, izerrstat)
    IF (izerrstat /= 0_iintegers) THEN
      yerror = ' ERROR *** while initializing MPI values for RTTOV7 ***'
      ierror = 9011
    ENDIF

    ! Allocate and initialize the minima and maxima of the profiles for
    ! temperature, humidity and ozone

    ! These variables could be different for different instruments.
    ! In fact, for Meteosat and MSG, which are only considered in the COSMO-Model
    ! there are only slight differences in the upper atmosphere, which do
    ! influence the results not significantly
    ! Therefore, no effort is made, to distinguish this
    ALLOCATE (utmx(jplev,num_sensors), utmn(jplev,num_sensors),                &
              uqmx(jplev,num_sensors), uqmn(jplev,num_sensors),                &
              uomx(jplev,num_sensors), uomn(jplev,num_sensors),                &
              ppres_d(jplev,num_sensors), o3_ref(jplev),                       &
              ivch(jpch,jpnsat),  STAT=izerrstat)

    ! this is only checked here, if RTTOV7 would have different numbers,
    ! it will hopefully complain, otherwise the namelist input is used
    ! for computations
    ivch(:,:)     = 0_iintegers

    CALL RTTOV7_rttvi_new (izerrstat, kppf, kpnsat, kplev, kpch, kpchus, kpnav,    &
         kpnsav, kpnssv, kpncv, num_sensors, sat_compute(:)%nsatell_table_id,  &
         sat_compute(:)%nsat_id, sat_compute(:)%nsensor_table_id, numchans,    &
         ppres_d, utmn, utmx, uqmn, uqmx, uomn, uomx, ivch, kiu1)

    IF (izerrstat /= 0_iintegers) THEN
      yerror = ' ERROR *** initializing RTTOV7_rttvi ***'
      ierror = 10_iintegers
      RETURN
    ENDIF

    ppres_d(:,:) = ppres_d(:,:) / 100.0_wp

    DEALLOCATE (ivch)    ! this is not really needed

    ! Set the reference ozone profile
    o3_ref(1:jplev) =                                                         &
          (/0.969339E-05_wp, 0.100043E-04_wp, 0.101194E-04_wp,    &
            0.101751E-04_wp, 0.102181E-04_wp, 0.102463E-04_wp,    &
            0.102127E-04_wp, 0.102454E-04_wp, 0.101497E-04_wp,    &
            0.935746E-05_wp, 0.809846E-05_wp, 0.672055E-05_wp,    &
            0.519183E-05_wp, 0.372400E-05_wp, 0.258175E-05_wp,    &
            0.172204E-05_wp, 0.119041E-05_wp, 0.845131E-06_wp,    &
            0.649756E-06_wp, 0.526723E-06_wp, 0.412770E-06_wp,    &
            0.303262E-06_wp, 0.210869E-06_wp, 0.156281E-06_wp,    &
            0.123483E-06_wp, 0.107569E-06_wp, 0.100696E-06_wp,    &
            0.959905E-07_wp, 0.916880E-07_wp, 0.891324E-07_wp,    &
            0.846609E-07_wp, 0.811542E-07_wp, 0.778113E-07_wp,    &
            0.757306E-07_wp, 0.711752E-07_wp, 0.663465E-07_wp,    &
            0.616040E-07_wp, 0.567638E-07_wp, 0.521222E-07_wp,    &
            0.479273E-07_wp, 0.443969E-07_wp, 0.419678E-07_wp,    &
            0.409984E-07_wp/)

    ! Set some constant numbers
    const_aheat = 287.0_wp / 1005.0_wp

    IF (izdebug > 5) THEN
      PRINT *, '     INITIALIZED RTTOV 7'
    ENDIF
#endif

#if defined(RTTOV9) || defined(RTTOV10)
  ELSEIF (itype_rttov == rttov_ifc_version) THEN
    ALLOCATE(instruments(3,num_sensors+nsensors), &
             addclouds(num_sensors+nsensors))

    ! instruments 1 to (1/2) are for synsat
    instruments(1,1:num_sensors) = sat_compute(1:num_sensors)%nsatell_table_id
    instruments(2,1:num_sensors) = sat_compute(1:num_sensors)%nsat_id
    instruments(3,1:num_sensors) = sat_compute(1:num_sensors)%nsensor_table_id
    addclouds(1:num_sensors)      = .True.
    mchans = max(numchans(1), numchans(2))

    ! instruments (2/3) to ... are for radiance feedbk
    DO i = 1,nsensors
      instruments(:,num_sensors + i) = sensors(i)%id
      addclouds(num_sensors + i)     = sensors(i)%addcloud
      mchans = MAX(mchans, SIZE(sensors(i)%channels))
    ENDDO  

    ALLOCATE(channels(mchans, num_sensors + nsensors))
    ALLOCATE(n_chans(num_sensors + nsensors))

    DO i = 1, num_sensors
      n_chans(i)                = sat_compute(i)%num_chan
      channels(1:n_chans(i), i) = sat_compute(i)%nchan_list(1:n_chans(i))
    ENDDO

    DO i = 1,nsensors
      n_chans(num_sensors+i) = SIZE(sensors(i)%channels)
      channels(1:n_chans(num_sensors+i),num_sensors +i) = sensors(i)%channels(:)
    ENDDO

    izstatus = rttov_init(  &
         instruments      , &
         channels         , &
         n_chans          , &
         my_cart_id       , &
         nproc            , &
         0                , &
         icomm_cart       , &
         appRegLim=.TRUE. , &
#ifdef RTTOV9
         readCloud=.TRUE.)
#else
         readCloud=addclouds)
#endif

    IF (izstatus /= NO_ERROR) THEN
      IF (my_cart_id == 0) THEN
        PRINT *, ' ERROR *** while initializing RTTOV', rttov_ifc_version, ' *** '
        PRINT *, ' RTTOV ERROR *** ', TRIM(rttov_ifc_errMsg(izstatus))
      ENDIF
      yerror = ' ERROR *** while initializing RTTOV   *** '
      WRITE (yerror(36:37),'(I2)') rttov_ifc_version
      ierror = 9005
    ELSE
      IF (izdebug > 5) THEN
        PRINT *, '     INITIALIZED RTTOV ', rttov_ifc_version
      ENDIF
    ENDIF

#endif

  ENDIF ! itype_rttov

ELSEIF (yaction == 'compute') THEN

  !----------------------------------------------------------------------------
  ! Section 4: Computations
  !----------------------------------------------------------------------------

  IF (lsynsat .AND. (nsat_steps(nsat_next) == ntstep) ) THEN
    ! set status to 0 again
    ! isynsat_stat is checked by the output routine whether output can be done
    ! (isynsat_stat = 0), or whether something has happened and no images are
    ! available for output (isynsat_stat /= 0)
    isynsat_stat = 0

    izdim      = ie
    jzdim      = je
    izdimstart = istartpar
    izdimend   = iendpar
    jzdimstart = jstartpar
    jzdimend   = jendpar

    ! retrieve the required microphysics tracers
    CALL trcr_get(ierror, idt_qv, ptr_tlev = nnow, ptr=qv)
    IF (ierror /= 0) THEN
      yerror = trcr_errorstr(ierror)
      CALL model_abort(my_cart_id, ierror, yerror, yzroutine)
    ENDIF
    CALL trcr_get(ierror, idt_qc, ptr_tlev=nnow, ptr=qc)
    IF (ierror /= 0) THEN
      yerror = trcr_errorstr(ierror)
      CALL model_abort(my_cart_id, ierror, yerror, yzroutine)
    ENDIF
    CALL trcr_get(ierror, idt_qg, ptr_tlev=nnow, ptr=qg)
    IF (ierror /= 0 .AND. ierror /= T_ERR_NOTFOUND) THEN
      yerror = trcr_errorstr(ierror)
      CALL model_abort(my_cart_id, ierror, yerror, yzroutine)
    ENDIF
    CALL trcr_get(ierror, idt_qi, ptr_tlev=nnow, ptr=qi)
    IF (ierror /= 0 .AND. ierror /= T_ERR_NOTFOUND) THEN
      yerror = trcr_errorstr(ierror)
      CALL model_abort(my_cart_id, ierror, yerror, yzroutine)
    ENDIF
    CALL trcr_get(ierror, idt_qs, ptr_tlev=nnow, ptr=qs)
    IF (ierror /= 0 .AND. ierror /= T_ERR_NOTFOUND) THEN 
      yerror = trcr_errorstr(ierror)
      CALL model_abort(my_cart_id, ierror, yerror, yzroutine)
    ENDIF  

    ! no. of profiles and output radiances
    nprof     = iendpar-istartpar+1
    maxknchpf = nprof * MAX (numchans(1), numchans(2))

    IF (itype_rttov == 7) THEN

#ifdef RTTOV7
      IF (izdebug > 5) THEN
        PRINT *, '      RTTOV 7 PREPARATION'
      ENDIF

      ! Check, whether number of profiles is not too large:
      IF (nprof > kppf) THEN
        PRINT *, ' ERROR *** number of profiles is too large   ***'
        PRINT *, '       *** ', nprof, ' > ', kppf, '  ***'
        ierror = 9020
        RETURN
      ENDIF

      syclc(:,:,:)=clc_sgs(:,:,:)+clc_con(:,:,:)*(1.0_wp - clc_sgs(:,:,:))
      syclc(:,:,:)=MAX( 0.0_wp, MIN( 1.0_wp,syclc(:,:,:)))

      IF ( ASSOCIATED(qi) ) THEN
        IF ( ASSOCIATED(qs) ) THEN
          DO  k = 1, ke
            DO i=istartpar,iendpar
              DO j=jstartpar,jendpar
                IF (qs(i,j,k) > 1.0E-7_wp) syclc(i,j,k)=1.0_wp
              ENDDO
            ENDDO
          ENDDO
          zfield(:,:,:) = qi(:,:,:)+qs(:,:,:)
        ELSE
          zfield(:,:,:) = qi(:,:,:)
        ENDIF
      ENDIF

      IF (lcon_clw) THEN
        syclw(:,:,:)=qc(:,:,:)+clw_con(:,:,:)
      ELSE
        syclw(:,:,:)=qc(:,:,:)
      ENDIF

      ! compute the synthetic pictures
      CALL org_sat_tbs(t(:,:,:,nnow), qv(:,:,:), syclw(:,:,:),               &
                      zfield(:,:,:), pp(:,:,:,nnow), p0(:,:,:), p0hl(:,:,:), &
                      syclc(:,:,:), ps(:,:,nnow), t_g(:,:,nnow), t_2m(:,:),  &
                      qv_2m(:,:), u_10m(:,:), v_10m(:,:), fr_land(:,:),      &
                      rlat(:,:), rlon(:,:), izdim, jzdim, ke,                &
                      izdimstart, izdimend, jzdimstart, jzdimend,            &
                      jplev, nprof, jpnsat, jpnav, jpnsav, jpnssv,           &
                      synme7(:,:,:), synmsg(:,:,:), yzerror, izerrstat)

      IF (izerrstat /= 0_iintegers) THEN
        yerror = ' ERROR *** while computing synthetic satellite images ***'
        ierror = 9021
      ENDIF
#endif

#if defined(RTTOV9) || defined(RTTOV10)
    ELSEIF (itype_rttov == rttov_ifc_version) THEN

      IF (izdebug > 5) THEN
        PRINT *, '      RTTOV 9 PREPARATION'
      ENDIF

      IF (lalloc_h_ice) THEN
         h_ice_dummy => h_ice(:,:,:)
      ELSE
         ALLOCATE(h_ice_dummy(0,0,nnow))
      ENDIF
      IF (ASSOCIATED(qg)) THEN
         qg_dummy => qg(:,:,:)
      ELSE
         ALLOCATE(qg_dummy(0,0,0))
      ENDIF

      ! compute the synthetic pictures
      CALL org_sat_rttov(t(:,:,:,nnow), qv(:,:,:), qc(:,:,:),                &
           qi(:,:,:), qs(:,:,:), qg_dummy(:,:,:),                            &
           pp(:,:,:,nnow), ps(:,:,nnow), t_g(:,:,nnow),                      &
           h_ice_dummy(:,:,nnow), izdim, jzdim, ke, izdimstart, izdimend,    &
           jzdimstart, jzdimend, num_sensors, synme7(:,:,:), synmsg(:,:,:),  &
           nusyn, yzerror, izstatus)
      IF (izstatus /= 0_iintegers) THEN
         yerror = ' ERROR *** while computing synthetic satellite images ***'
         ierror = 9021
      ENDIF

      IF (.NOT. lalloc_h_ice)    DEALLOCATE(h_ice_dummy)
      IF (.NOT.ASSOCIATED  (qg)) DEALLOCATE(qg_dummy)
#endif

    ENDIF ! itype_rttov

  ENDIF

#ifdef NUDGING
  IF (lobsrad) THEN
    CALL calc_obs_satpp(izerrstat)

    IF (izerrstat /= 0) THEN
      PRINT *, ' ERROR *** while calculating obs *** '
      ierror = 9010
      RETURN 
    ENDIF
  ENDIF
#endif

ELSEIF (yaction == 'dealloc') THEN

  !----------------------------------------------------------------------------
  ! Section 5: Deallocation at the end of each time step
  !----------------------------------------------------------------------------

  IF (lsynsat .AND. (nsat_steps(nsat_next) == ntstep) ) THEN
    nsat_next = nsat_next+1
  ENDIF

ELSEIF (yaction == 'cleanup') THEN

  !----------------------------------------------------------------------------
  ! Section 6: Cleanup memory at the end of the program
  !----------------------------------------------------------------------------

#ifdef NUDGING
  ! Reading of NWCSAF SEVIRI cloud type
  IF (lread_ct) THEN
    CALL get_utc_date ( ntstep, ydate_ini, dt, itype_calendar, yzactdate       &
                      , yzdum, izdum,rzdum)
!   =================
    CALL cloud_read_prep(yzactdate,yaction,ctype,cth,izerrstat)
    IF (izerrstat /= 0) THEN
      PRINT *, ' ERROR *** while reading cloud type GRIB file *** '
      ierror = 9100
      RETURN
    ENDIF
  ENDIF
#endif

#ifdef RTTOV7
  ! ivch has already been deallocated
  DEALLOCATE (utmx, utmn, uqmx, uqmn, uomx, uomn, ppres_d, o3_ref,        STAT=izerrstat)
  IF (izerrstat /= 0) THEN
    PRINT *, ' ERROR *** while deallocating memory for RTTOV7 *** '
    ierror = 9102
    RETURN
  ENDIF
#endif

#if defined(RTTOV9) || defined(RTTOV10)
  IF (itype_rttov == rttov_ifc_version) THEN
    izstatus = rttov_cleanup()
    IF (izstatus /= NO_ERROR) THEN
      write(*,'("*** RTTOV ERROR (proc ",I2,") ",A)') my_world_id,&
           trim(rttov_ifc_errMsg(izstatus))
      write(yerror,'(" ERROR *** while computing synthetic satellite images: &
           &rttov_cleanup ",I3)') izstatus
      ierror = 9008
    ENDIF
  ENDIF ! itype_rttov

  ! deallocate some memory from organize_satellites
  DEALLOCATE (instruments, addclouds, channels, n_chans, STAT=izerrstat)
  IF (izerrstat /= 0) THEN
    PRINT *, ' ERROR *** while deallocating memory for RTTOV', rttov_ifc_version,' *** '
    ierror = 9102
    RETURN
  ENDIF
#endif

  DEALLOCATE (nsat_steps, STAT=izerrstat)
  IF (izerrstat /= 0) THEN
    PRINT *, ' ERROR *** while deallocating memory for RTTOV: nsat_steps *** '
    ierror = 9102
    RETURN
  ENDIF

!------------------------------------------------------------------------------
! Section 6: All other actions are wrong
!------------------------------------------------------------------------------

ELSE

  ierror = 1
  yerror = 'ERROR *** Action "'//TRIM(yaction)//&
          &'" not valid action for satellite computations ***'

ENDIF

!------------------------------------------------------------------------------
! Internal Procedures
!------------------------------------------------------------------------------

CONTAINS

!==============================================================================
!+ Internal procedure in "organize_data" for the input of NAMELIST satctl
!------------------------------------------------------------------------------

SUBROUTINE input_satctl (nuspecif, nuin, ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine organizes the input of the NAMELIST-group satctl.
!   The group contains variables for the computation of satellite radiances
!   and brightness temperatures
!
! Method:
!
!------------------------------------------------------------------------------

! Parameters
  INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
    nuspecif,     & ! Unit number for protocolling the task
    nuin            ! Unit number for Namelist INPUT file

  INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
    ierrstat        ! error status variable

! Variables for default values
  TYPE (sat_input_type)      ::       &
    sat_input_01, sat_input_02    ! for specifying up to 2 different sensors

  INTEGER (KIND=iintegers)   ::       &
    num_sensors_d,        & ! number of sensors used during the computations
                            ! should be less than MIN(jpnsat, 2)
                            ! NOTE: we are assuming that jpnsat <= 2!!
    itype_rttov_d,        & ! Version of the RTTOV library to be used
    extrp_type_d,         & !
    iceshape_d,           & !
    iwc2effdiam_d,        & !
    nchan_input_01(jpch), & ! input channel list for first   sensor
    nchan_input_02(jpch)    ! input channel list for second sensor

  REAL    (KIND=wp)          ::       &
    emiss_input_01(jpch), & ! for reading input emissivities
    emiss_input_02(jpch)    ! for reading input emissivities

  REAL    (KIND=wp)        ::       &
    sat_long_01,     & !
    sat_long_02,     & !
    sat_long_01_d,   & !
    sat_long_02_d

  LOGICAL                    ::       &
    ldouble, lchange, lcon_clw_d, lsynsat_d, lobsrad_d, lread_ct_d, linterp_d

  CHARACTER (LEN=  12) :: yname
  CHARACTER (LEN= 100) :: yclouddir_d

  INTEGER (KIND=iintegers)   :: ierr, ic, ic_int, ic_real, isens, ipic,     &
                                nc, nc2, i, tmp

  ! Pointer for the linked list of namelist GRIBOUTs
  TYPE(pp_nl) , POINTER      ::  now, next
  INTEGER (KIND=iintegers), ALLOCATABLE   :: nstep_sat(:)

  CHARACTER(LEN=250)         :: iomsg_str

! Define the namelist group
  NAMELIST /satctl/ num_sensors, sat_input_01, sat_input_02,            &
                    nchan_input_01, nchan_input_02,                     &
                    emiss_input_01, emiss_input_02,                     &
                    sat_long_01,    sat_long_02,                        &
                    itype_rttov, lsynsat, lobsrad,                      &
                    lcon_clw, extrp_type, iceshape, iwc2effdiam,        &
                    lread_ct, yclouddir, linterp

!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE input_satctl
!-------------------------------------------------------------------------------

! Initialize the RTTOV-tables
  yrttov_satell_table( 1) = 'NOAA    '
  yrttov_satell_table( 2) = 'DMSP    '
  yrttov_satell_table( 3) = 'METEOSAT'
  yrttov_satell_table( 4) = 'GOES    '
  yrttov_satell_table( 5) = 'GMS     '
  yrttov_satell_table( 6) = 'FY-2    '
  yrttov_satell_table( 7) = 'TRMM    '
  yrttov_satell_table( 8) = 'ERS     '
  yrttov_satell_table( 9) = 'EOS     '
  yrttov_satell_table(10) = 'xxxxxxxx'
  yrttov_satell_table(11) = 'ENVISAT '
  yrttov_satell_table(12) = 'MSG     '
  yrttov_satell_table(13) = 'FY-1    '

  yrttov_sensor_table( 0) = 'HIRS        '
  yrttov_sensor_table( 1) = 'MSU         '
  yrttov_sensor_table( 2) = 'SSU         '
  yrttov_sensor_table( 3) = 'AMSU-A      '
  yrttov_sensor_table( 4) = 'AMSU-B      '
  yrttov_sensor_table( 5) = 'AVHRR       '
  yrttov_sensor_table( 6) = 'SSMI        '
  yrttov_sensor_table( 7) = 'VTPR1       '
  yrttov_sensor_table( 8) = 'VTPR2       '
  yrttov_sensor_table( 9) = 'TMI         '
  yrttov_sensor_table(10) = 'SSMIS       '
  yrttov_sensor_table(11) = 'AIRS        '
  yrttov_sensor_table(12) = 'xxxxxxxxxxxx'
  yrttov_sensor_table(13) = 'MODIS       '
  yrttov_sensor_table(14) = 'ATSR        '
  yrttov_sensor_table(15) = 'xxxxxxxxxxxx'
  yrttov_sensor_table(16) = 'xxxxxxxxxxxx'
  yrttov_sensor_table(17) = 'xxxxxxxxxxxx'
  yrttov_sensor_table(18) = 'xxxxxxxxxxxx'
  yrttov_sensor_table(19) = 'xxxxxxxxxxxx'
  yrttov_sensor_table(20) = 'MVIRI       '
  yrttov_sensor_table(21) = 'SEVIRI      '
  yrttov_sensor_table(22) = 'GOES-IMAGER '
  yrttov_sensor_table(23) = 'GOES-SOUNDER'
  yrttov_sensor_table(24) = 'GMS-IMAGER  '
  yrttov_sensor_table(25) = 'FY2-VISSR   '
  yrttov_sensor_table(26) = 'FY2-MVISR   '

IF (my_world_id == 0) THEN

!-------------------------------------------------------------------------------
!- Section 1: Initialize the default variables
!-------------------------------------------------------------------------------

  num_sensors_d     = 0_iintegers
  lcon_clw_d        = .FALSE.
  itype_rttov_d     = itype_rttov
  extrp_type_d      = extrp_type
  iceshape_d        = iceshape
  iwc2effdiam_d     = iwc2effdiam
  linterp_d         = .TRUE.
  lsynsat_d         = .TRUE.
  lobsrad_d         = .FALSE.
  lread_ct_d        = .FALSE.
  yclouddir_d       = ''
  sat_long_01_d     = -999.0_wp
  sat_long_02_d     = -999.0_wp

!-------------------------------------------------------------------------------
!- Section 2: Initialize variables with defaults
!-------------------------------------------------------------------------------

  num_sensors       = num_sensors_d
  lcon_clw          = lcon_clw_d
  sat_input_01      = sat_input_type('yyyyyyyy', 0, 'yyyyyyyyyyyy', 0,       &
                                      .FALSE., .FALSE., .FALSE., .FALSE.)
  sat_input_02      = sat_input_01
  nchan_input_01(:) = 0_iintegers
  nchan_input_02(:) = 0_iintegers
  emiss_input_01(:) = 0.0_wp
  emiss_input_02(:) = 0.0_wp
  linterp           = linterp_d
  lsynsat           = lsynsat_d
  lobsrad           = lobsrad_d
  lread_ct          = lread_ct_d
  yclouddir         = yclouddir_d
  sat_long_01       = sat_long_01_d
  sat_long_02       = sat_long_02_d

!-------------------------------------------------------------------------------
!- Section 3: Input of the namelist values
!-------------------------------------------------------------------------------

  iomsg_str(:) = ' '
  READ (nuin, satctl, IOSTAT=ierr, IOMSG=iomsg_str)

  IF (ierr /= 0) WRITE (*,'(A,A)') 'Namelist-ERROR SATCTL: ', TRIM(iomsg_str)
ENDIF

IF (nproc > 1) THEN
  ! distribute error status to all processors
  CALL distribute_values  (ierr, 1, 0, imp_integers,  icomm_world, ierr)
ENDIF

IF (ierr /= 0) THEN
  PRINT *, 'Reading of SATCTL failed  ', ierr
  ierrstat = -1
  RETURN
ENDIF

IF (my_world_id == 0) THEN

!-------------------------------------------------------------------------------
!- Section 4: Check values for errors and consistency
!-------------------------------------------------------------------------------

  SELECT CASE (itype_rttov)
  CASE (7)

#ifndef RTTOV7
    PRINT *, ' ERROR    *** itype_rttov = 7, but RTTOV7 is not defined for compilation *** '
    ierrstat = 7002
#endif

  CASE (9)

#ifndef RTTOV9
    PRINT *, ' ERROR    *** itype_rttov = 9, but RTTOV9 is not defined for compilation *** '
    ierrstat = 7002
#endif

  CASE (10)

#ifndef RTTOV10
    PRINT *, ' ERROR    *** itype_rttov = 10, but RTTOV10 is not defined for compilation *** '
    ierrstat = 7002
#endif

  CASE DEFAULT
    PRINT *, ' ERROR    *** Wrong Version of RTTOV model *** '
    PRINT *, '          *** itype_rttov = ', itype_rttov
    PRINT *, '          *** Only itype_rttov = 7 / 9 / 10 are valid, if the proper -DRTTOVx is set ***'
    ierrstat = 7002
  END SELECT

#ifndef RTTOV10
  IF (lobsrad) THEN
    PRINT *, ' ERROR    *** lobsrad must not be TRUE, because model not compiled for RTTOV10 ***'
    ierrstat = 7002
  ENDIF
#endif

#ifndef NUDGING
  IF (lobsrad) THEN
    PRINT *, ' ERROR    *** lobsrad must not be TRUE, because model not compiled for NUDGING ***'
    ierrstat = 7002
  ENDIF
#endif

  IF (num_sensors > jpnsat) THEN
    PRINT *, ' ERROR    *** Number of sensors is bigger than allowed *** '
    PRINT *, '          *** num_sensors = ', num_sensors
    PRINT *, '          *** Max. number of sensors = ', jpnsat
    ierrstat = 7002
  ENDIF

  IF (num_sensors <= 0) THEN
    PRINT *, ' ERROR    *** Number of sensors is not set *** '
  ENDIF

  IF (itype_gscp < 3) THEN
    PRINT *, ' ERROR    *** SYNSAT with itype_gscp < 3 not supported *** '
    PRINT *, '          *** Change itype_gscp or do not set luse_rttov *** '
    ierrstat = 7002
  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!- Section 5: Distribute variables to all nodes
!-------------------------------------------------------------------------------

IF (nproc > 1) THEN

  IF (my_world_id == 0) THEN
    intbuf ( 1) = num_sensors
    intbuf ( 2) = sat_input_01%nsat_id
    intbuf ( 3) = sat_input_01%num_chan
    intbuf ( 4) = sat_input_02%nsat_id
    intbuf ( 5) = sat_input_02%num_chan
    intbuf ( 6) = extrp_type
    intbuf ( 7) = itype_rttov
    intbuf ( 8) = iwc2effdiam
    intbuf ( 9) = iceshape

    ic = 9
    intbuf (ic+1:ic+jpch) = nchan_input_01(1:jpch)
    ic = ic + jpch
    intbuf (ic+1:ic+jpch) = nchan_input_02(1:jpch)
    ic = ic + jpch

    realbuf( 1) = sat_long_01
    realbuf( 2) = sat_long_02

    ic = 2
    realbuf(ic+1:ic+jpch) = emiss_input_01(1:jpch)
    ic = ic + jpch
    realbuf(ic+1:ic+jpch) = emiss_input_02(1:jpch)
    ic = ic + jpch

    logbuf ( 1) = lcon_clw
    logbuf ( 2) = lsynsat
    logbuf ( 3) = lobsrad
    logbuf ( 4) = sat_input_01%lclear_rad
    logbuf ( 5) = sat_input_01%lcloud_rad
    logbuf ( 6) = sat_input_01%lclear_tem
    logbuf ( 7) = sat_input_01%lcloud_tem
    logbuf ( 8) = sat_input_02%lclear_rad
    logbuf ( 9) = sat_input_02%lcloud_rad
    logbuf (10) = sat_input_02%lclear_tem
    logbuf (11) = sat_input_02%lcloud_tem
    logbuf (12) = lread_ct
    logbuf (13) = linterp

    charbuf( 1) = sat_input_01%ysatellite
    charbuf( 2) = sat_input_01%ysensor
    charbuf( 3) = sat_input_02%ysatellite
    charbuf( 4) = sat_input_02%ysensor
    charbuf( 5) = yclouddir
  ENDIF
  ic_int  = 9 + 2*jpch
  ic_real = 2 + 2*jpch

  CALL distribute_values (intbuf,  ic_int, 0, imp_integers,  icomm_world, ierr)
  CALL distribute_values (realbuf,ic_real, 0, imp_reals,     icomm_world, ierr)
  CALL distribute_values (logbuf,      13, 0, imp_logical,   icomm_world, ierr)
  CALL distribute_values (charbuf,      5, 0, imp_character, icomm_world, ierr)

  IF (my_world_id /= 0) THEN
    num_sensors                = intbuf ( 1)
    sat_input_01%nsat_id       = intbuf ( 2)
    sat_input_01%num_chan      = intbuf ( 3)
    sat_input_02%nsat_id       = intbuf ( 4)
    sat_input_02%num_chan      = intbuf ( 5)
    extrp_type                 = intbuf ( 6)
    itype_rttov                = intbuf ( 7)
    iwc2effdiam                = intbuf ( 8)
    iceshape                   = intbuf ( 9)

    ic = 9
    nchan_input_01(1:jpch) = intbuf (ic+1:ic+jpch)
    ic = ic + jpch
    nchan_input_02(1:jpch) = intbuf (ic+1:ic+jpch)

    sat_long_01                = realbuf( 1)
    sat_long_02                = realbuf( 2)

    ic = 2
    emiss_input_01(1:jpch) = realbuf(ic+1:ic+jpch)
    ic = ic + jpch
    emiss_input_02(1:jpch) = realbuf(ic+1:ic+jpch)

    lcon_clw                 = logbuf ( 1)
    lsynsat                  = logbuf ( 2)
    lobsrad                  = logbuf ( 3)
    sat_input_01%lclear_rad  = logbuf ( 4)
    sat_input_01%lcloud_rad  = logbuf ( 5)
    sat_input_01%lclear_tem  = logbuf ( 6)
    sat_input_01%lcloud_tem  = logbuf ( 7)
    sat_input_02%lclear_rad  = logbuf ( 8)
    sat_input_02%lcloud_rad  = logbuf ( 9)
    sat_input_02%lclear_tem  = logbuf (10)
    sat_input_02%lcloud_tem  = logbuf (11)
    lread_ct                 = logbuf (12)
    linterp                  = logbuf (13)

    sat_input_01%ysatellite  = charbuf( 1)(1:8)
    sat_input_01%ysensor     = charbuf( 2)(1:12)
    sat_input_02%ysatellite  = charbuf( 3)(1:8)
    sat_input_02%ysensor     = charbuf( 4)(1:12)
    yclouddir                = charbuf( 5)
  ENDIF

ENDIF

IF (my_world_id == 0) THEN
  call get_free_unit(nusyn)
  OPEN(nusyn, FILE=yusyn, FORM='FORMATTED', STATUS='UNKNOWN', &
       POSITION='APPEND', IOSTAT=ierr)
  IF (ierr /= 0) THEN
    ierrstat = 99101
    PRINT *, ' ERROR   *** OPENING OF FILE yusyn FAILED. ***'
    RETURN
  ENDIF
END IF

! Now fill the sat_org_type
IF (num_sensors > 0 .AND. (SIZE(sat_compute) >= 1)) THEN
  sat_compute(1)%ysatellite         = sat_input_01%ysatellite
  sat_compute(1)%nsat_id            = sat_input_01%nsat_id
  sat_compute(1)%wmo_satid          = wmo_satid(sat_compute(1)%ysatellite,sat_compute(1)%nsat_id)
  IF ((sat_long_01 >= -180.0_wp) .AND. (sat_long_01 <= 360.0_wp)) THEN
    sat_compute(1)%longitude        = sat_long_01 * degrad
  ELSE
    sat_compute(1)%longitude        = sat_longitude(sat_compute(1)%wmo_satid)
  ENDIF
  sat_compute(1)%ysensor            = sat_input_01%ysensor
  sat_compute(1)%num_chan           = sat_input_01%num_chan
  sat_compute(1)%nchan_list(1:jpch) = nchan_input_01(1:jpch)
  DO nc = 1, sat_input_01%num_chan
    sat_compute(1)%ychan_name(nc)   = set_channel_name(sat_compute(1)%ysensor, sat_compute(1)%nchan_list(nc))
  ENDDO
  sat_compute(1)%emissivity(1:jpch) = emiss_input_01(1:jpch)
  sat_compute(1)%lclear_rad         = sat_input_01%lclear_rad
  sat_compute(1)%lcloud_rad         = sat_input_01%lcloud_rad
  sat_compute(1)%lclear_tem         = sat_input_01%lclear_tem
  sat_compute(1)%lcloud_tem         = sat_input_01%lcloud_tem
  IF (my_world_id == 0) THEN
    WRITE (nusyn,*) '  '
    WRITE(nusyn,'("Sat 1              = ",A,1x,I3,1x,I3)') &
         TRIM(sat_input_01%ysatellite),sat_input_01%nsat_id,sat_compute(1)%wmo_satid
    WRITE(nusyn,'("Sat 1 longitude    = ",F7.3)') sat_compute(1)%longitude
    WRITE(nusyn,'("Sat 1 sensor       = ",A)') trim(sat_compute(1)%ysensor)
    WRITE(nusyn,'("Sat 1 n_chan       = ",I4)') sat_compute(1)%num_chan
    WRITE(nusyn,'("Sat 1 channels     = ")')
    DO nc = 1, sat_input_01%num_chan
      WRITE(nusyn,'("                     ",I3,1x,A)') sat_compute(1)%nchan_list(nc),sat_compute(1)%ychan_name(nc)
    ENDDO
  ENDIF

ENDIF

IF (num_sensors > 1 .AND. (SIZE(sat_compute) >= 2)) THEN
  sat_compute(2)%ysatellite         = sat_input_02%ysatellite
  sat_compute(2)%nsat_id            = sat_input_02%nsat_id
  sat_compute(2)%wmo_satid          = wmo_satid(sat_compute(2)%ysatellite,sat_compute(2)%nsat_id)
  IF ((sat_long_02 >= -180.0_wp) .AND. (sat_long_02 <= 360.0_wp)) THEN
    sat_compute(2)%longitude        = sat_long_02 * degrad
  ELSE
    sat_compute(2)%longitude        = sat_longitude(sat_compute(2)%wmo_satid)
  ENDIF
  sat_compute(2)%ysensor            = sat_input_02%ysensor
  sat_compute(2)%num_chan           = sat_input_02%num_chan
  sat_compute(2)%nchan_list(1:jpch) = nchan_input_02(1:jpch)
  DO nc = 1, sat_input_02%num_chan
    sat_compute(2)%ychan_name(nc)   = set_channel_name(sat_compute(2)%ysensor, sat_compute(2)%nchan_list(nc))
  ENDDO
  sat_compute(2)%emissivity(1:jpch) = emiss_input_02(1:jpch)
  sat_compute(2)%lclear_rad         = sat_input_02%lclear_rad
  sat_compute(2)%lcloud_rad         = sat_input_02%lcloud_rad
  sat_compute(2)%lclear_tem         = sat_input_02%lclear_tem
  sat_compute(2)%lcloud_tem         = sat_input_02%lcloud_tem
  IF (my_world_id == 0) THEN
    WRITE (nusyn,*) '  '
    WRITE(nusyn,'("Sat 2              = ",A,1x,I3,1x,I3)') &
         TRIM(sat_input_01%ysatellite),sat_input_01%nsat_id,sat_compute(2)%wmo_satid
    WRITE(nusyn,'("Sat 2 longitude    = ",F7.3)') sat_compute(2)%longitude
    WRITE(nusyn,'("Sat 2 sensor       = ",A)') trim(sat_compute(2)%ysensor)
    WRITE(nusyn,'("Sat 2 n_chan       = ",I4)') sat_compute(2)%num_chan
    WRITE(nusyn,'("Sat 2 channels     = ")')
    DO nc = 1, sat_input_01%num_chan
      WRITE(nusyn,'("                     ",I3,1x,A)') sat_compute(2)%nchan_list(nc),sat_compute(2)%ychan_name(nc)
    ENDDO
  ENDIF
ENDIF

DO isens = 1, num_sensors
 
  ! Find the entries in the rttov satellite table and sensor table:
  ! satellites
  yname = sat_compute(isens)%ysatellite(1:LEN_TRIM(sat_compute(isens)%ysatellite))
  DO ic = 1, 13
    IF (yrttov_satell_table(ic) == yname) THEN
      sat_compute(isens)%nsatell_table_id = ic
      EXIT
    ENDIF
  ENDDO
 
  ! sensors
  yname = sat_compute(isens)%ysensor(1:LEN_TRIM(sat_compute(isens)%ysensor))
  DO ic = 1, 26
    IF (yrttov_sensor_table(ic) == yname) THEN
      sat_compute(isens)%nsensor_table_id = ic
      EXIT
    ENDIF
  ENDDO

  ipic = 0
  DO ic = 1, sat_compute(isens)%num_chan
    ! Specify additional element numbers and channels for the grib output
    IF (sat_compute(isens)%lcloud_tem) THEN
      ipic = ipic + 1
      sat_compute(isens)%ngrib_chan(ipic) = sat_compute(isens)%nchan_list(ic)
      sat_compute(isens)%ngrib_aees(ipic) = 1
    ENDIF
  
    IF (sat_compute(isens)%lclear_tem) THEN
      ipic = ipic + 1
      sat_compute(isens)%ngrib_chan(ipic) = sat_compute(isens)%nchan_list(ic)
      sat_compute(isens)%ngrib_aees(ipic) = 2
    ENDIF

    IF (sat_compute(isens)%lcloud_rad) THEN
      ipic = ipic + 1
      sat_compute(isens)%ngrib_chan(ipic) = sat_compute(isens)%nchan_list(ic)
      sat_compute(isens)%ngrib_aees(ipic) = 3
    ENDIF

    IF (sat_compute(isens)%lclear_rad) THEN
      ipic = ipic + 1
      sat_compute(isens)%ngrib_chan(ipic) = sat_compute(isens)%nchan_list(ic)
      sat_compute(isens)%ngrib_aees(ipic) = 4
    ENDIF
  ENDDO
  sat_compute(isens)%ndim3_field = ipic

ENDDO

! Determine the time steps for which satellite computations must be done
nc = 0
ALLOCATE (nstep_sat(10*noutst))

IF ( ASSOCIATED(root) .AND. (root%outsteps > 0) ) THEN
  now => root
  IF (now%yvarsl(1) /= '') THEN
    nc  = now%outsteps
    nstep_sat(1:nc) = now%ngrib(1:nc)
  ELSE
    nc  = 0
  ENDIF

  DO WHILE (ASSOCIATED(now%next))
    now => now%next
    IF (now%yvarsl(1) /= '') THEN
      nstep_sat(nc+1:nc+now%outsteps) = now%ngrib(1:now%outsteps)
      nc2 = nc + now%outsteps

      ! bubblesort
      DO
        lchange = .FALSE.
        DO i=1,nc2-1
          IF(nstep_sat(i) > nstep_sat(i+1)) THEN
            tmp        = nstep_sat(i+1)
            nstep_sat(i+1) = nstep_sat(i)
            nstep_sat(i)   = tmp
            lchange = .TRUE.
          ENDIF
        ENDDO
        IF (.NOT. lchange) EXIT
      ENDDO

      ! eliminate identical entries
      DO
        ldouble = .FALSE.
        DO i=1,nc2-1
          IF(nstep_sat(i) == nstep_sat(i+1)) THEN
            nstep_sat(i:nc2-1) = nstep_sat(i+1:nc2)
            ldouble = .TRUE.
            EXIT
          ENDIF
        ENDDO
        IF (.NOT. ldouble) EXIT
        nc2 = nc2-1
      ENDDO

      nc = nc2
    ENDIF

  ENDDO
ENDIF

! If no satellite output is required, nc=0
! nsat_steps should be allocated also in this case, because it is accessed
! during initialization and compute-phase.
IF (nc > 0) THEN
  ALLOCATE (nsat_steps(nc+1), STAT=nc2)
  nsat_steps(1:nc) = nstep_sat(1:nc)
  nsat_steps(nc+1) = -1   ! this is the end condition
ELSE
  ALLOCATE (nsat_steps(1), STAT=nc2)
  nsat_steps(1) = -1      ! this is the end condition
ENDIF

DEALLOCATE(nstep_sat)

!-------------------------------------------------------------------------------
!- Section 6: Output of the namelist variables and their default values
!-------------------------------------------------------------------------------

IF (my_world_id == 0) THEN

  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(A)') '0     NAMELIST:  satctl'
  WRITE (nuspecif, '(A)') '      -----------------'
  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(T7,A,T33,A,T51,A,T70,A)') 'Variable', 'Actual Value',   &
                                               'Default Value', 'Format'

  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                              'num_sensors',num_sensors,num_sensors_d, ' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                            'itype_rttov', itype_rttov, itype_rttov_d, ' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                 'extrp_type',extrp_type,extrp_type_d, ' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                              'iwc2effdiam',iwc2effdiam,iwc2effdiam_d, ' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                       'iceshape',iceshape,iceshape_d, ' I '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                            'lcon_clw  ',  lcon_clw  , lcon_clw_d  ,   ' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                             'lsynsat   ',  lsynsat   , lsynsat_d  ,   ' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                             'lobsrad   ',  lobsrad   , lobsrad_d  ,   ' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                             'lread_ct  ',  lread_ct , lread_ct_d ,    ' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                             'linterp   ',  linterp   , linterp_d  ,   ' L '
  WRITE (nuspecif, '(A2)')  '  '

  DO ic = 1, num_sensors
    WRITE (nuspecif, '(T8,A,i2.2,A,T33,A8  ,T52,A8  ,T71,A3)')                      &
         'sat_input_',ic,'%ysatellite',TRIM(sat_input_01%ysatellite), 'yyyyyyyy',' C8 '
    WRITE (nuspecif, '(T8,A,i2.2,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
         'sat_input_',ic,'%nsat_id',sat_input_01%nsat_id, 0, ' I '
    WRITE (nuspecif, '(T8,A,i2.2,A,T33,A12  ,T52,A12  ,T71,A3)')                      &
         'sat_input_',ic,'%ysensor',TRIM(sat_input_01%ysensor), 'yyyyyyyyyyyy',' C12 '
    WRITE (nuspecif, '(T8,A,i2.2,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
         'sat_input_',ic,'%num_chan',sat_input_01%num_chan, 0, ' I '
    WRITE (nuspecif, '(T8,A,i2.2,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
         'sat_input_',ic,'%lclear_rad',sat_input_01%lclear_rad, .FALSE., ' I '
    WRITE (nuspecif, '(T8,A,i2.2,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
         'sat_input_',ic,'%lcloud_rad',sat_input_01%lcloud_rad, .FALSE., ' I '
    WRITE (nuspecif, '(T8,A,i2.2,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
         'sat_input_',ic,'%lclear_tem',sat_input_01%lclear_tem, .FALSE., ' I '
    WRITE (nuspecif, '(T8,A,i2.2,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
         'sat_input_',ic,'%lcloud_tem',sat_input_01%lcloud_tem, .FALSE., ' I '
    WRITE (nuspecif, '(A2)')  '  '
  END DO
  WRITE (nuspecif, '(A2)')  '  '

  WRITE (nuspecif, '(A)') '0 Info on selected Satellites:'
  WRITE (nuspecif, '(A2)')  '  '
  DO ic = 1, num_sensors
    WRITE (nuspecif, '("0",T10,A,T46,A)')                                       &
         'Name of satellite:',  TRIM(sat_compute(ic)%ysatellite)
    WRITE (nuspecif, '("0",T10,A,T46,I8)')                                       &
         'Satellite ID:',  sat_compute(ic)%nsat_id
    WRITE (nuspecif, '("0",T10,A,T46,A)')                                       &
         'Name of sensor:',  TRIM(sat_compute(ic)%ysensor)
    WRITE (nuspecif, '("0",T10,A,T46,I8)')                                       &
          'Entry in RTTOV satellite table:    ', sat_compute(ic)%nsatell_table_id
    WRITE (nuspecif, '("0",T10,A,T46,I8)')                                       &
          'Entry in RTTOV sensor table:       ', sat_compute(ic)%nsensor_table_id
    WRITE (nuspecif, '("0",T10,A,T46,I8)')                                       &
          'Number of channels used:           ', sat_compute(ic)%num_chan
    WRITE (nuspecif, '("0",T10,A,T44,F8.3,A)')                                   &
                     'Satellite longitude:    ', sat_compute(ic)%longitude*raddeg, 'deg E'
    WRITE (nuspecif, '("0",T10,A,T46,L8)')                                       &
          'Clear Sky Radiance selected:       ', sat_compute(ic)%lclear_rad
    WRITE (nuspecif, '("0",T10,A,T46,L8)')                                       &
          'Cloudy Sky Radiance selected:      ', sat_compute(ic)%lcloud_rad
    WRITE (nuspecif, '("0",T10,A,T46,L8)')                                       &
          'Clear Sky Temperature selected:    ', sat_compute(ic)%lclear_tem
    WRITE (nuspecif, '("0",T10,A,T46,L8)')                                       &
          'Cloudy Sky Temperature selected:   ', sat_compute(ic)%lcloud_tem
    WRITE (nuspecif, '("0",T10,A,T30,10I4)')                                     &
          'Selected Channels:',                                              &
          sat_compute(ic)%nchan_list(1:sat_compute(ic)%num_chan)
    WRITE (nuspecif, '("0",T10,A,T30,10F4.1)')                                   &
          'Chosen Emissivities:',                                            &
          sat_compute(ic)%emissivity(1:sat_compute(ic)%num_chan)
    WRITE (nuspecif, '(A2)')  '  '
  ENDDO
  WRITE (nuspecif, '(A)') '0     END OF NAMELIST  satctl'
  WRITE (nuspecif, '(A2)')  '  '

ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE input_satctl

!==============================================================================
!==============================================================================
!+ Internal procedure for defining the satellite ID
!------------------------------------------------------------------------------

INTEGER FUNCTION wmo_satid(yname, id)
  CHARACTER(LEN=*),        INTENT(IN)           :: yname
  INTEGER(KIND=iintegers), INTENT(IN), OPTIONAL :: id

  INTEGER(KIND=iintegers) :: idh
  CHARACTER(len=80)       :: hname

  IF (.not.present(id)) THEN
    read(yname,*) hname, idh
  ELSE
    hname = yname
    idh   = id
  ENDIF

  SELECT CASE(trim(hname))
  CASE('METEOSAT','meteosat')
    SELECT CASE(idh)
    CASE(1)
      wmo_satid = 58
    CASE(2)
      wmo_satid = 59
    CASE(3)
      wmo_satid = 50
    CASE(4)
      wmo_satid = 51
    CASE(5)
      wmo_satid = 52
    CASE(6)
      wmo_satid = 53
    CASE(7)
      wmo_satid = 54
    CASE(8)
      wmo_satid = 55
    CASE(9)
      wmo_satid = 56
    CASE(10)
      wmo_satid = 57
    CASE(11)
      wmo_satid = 70
    END SELECT
  CASE('MSG','msg')
    SELECT CASE(idh)
    CASE(1)
      wmo_satid = 55
    CASE(2)
      wmo_satid = 56
    CASE(3)
      wmo_satid = 57
    CASE(4)
      wmo_satid = 70
    END SELECT
  CASE('ENVISAT','envisat')
    wmo_satid = 60
  END SELECT

END FUNCTION wmo_satid

!==============================================================================
!==============================================================================
!+ Internal procedure for defining the channel name
!------------------------------------------------------------------------------

CHARACTER (LEN=10) FUNCTION set_channel_name (yname, nchan)

  CHARACTER(LEN=*),        INTENT(IN) :: yname
  INTEGER(KIND=iintegers), INTENT(IN) :: nchan

  CHARACTER(LEN=10)  :: ychan


  ychan = '          '

  SELECT CASE(TRIM(yname))
  CASE('MVIRI','mviri')

    SELECT CASE(nchan)
    CASE(1)
      ychan = 'WV6.4     '
    CASE(2)
      ychan = 'IR11.5    '
    END SELECT

  CASE('SEVIRI','seviri')

    SELECT CASE(nchan)
    CASE(1)
      ychan = 'IR3.9     '
    CASE(2)
      ychan = 'WV6.2     '
    CASE(3)
      ychan = 'WV7.3     '
    CASE(4)
      ychan = 'IR8.7     '
    CASE(5)
      ychan = 'IR9.7     '
    CASE(6)
      ychan = 'IR10.8    '
    CASE(7)
      ychan = 'IR12.1    '
    CASE(8)
      ychan = 'IR13.4    '
    END SELECT

  END SELECT

  set_channel_name = ychan

END FUNCTION set_channel_name

!==============================================================================
!==============================================================================
!+ Internal procedure for defining the satellite longitude
!------------------------------------------------------------------------------

REAL(KIND=wp)     FUNCTION sat_longitude(wmo_satid)

  INTEGER(KIND=iintegers), INTENT(IN) :: wmo_satid

  SELECT CASE(wmo_satid)
  CASE(54)
    sat_longitude = 57.5_wp * degrad
  CASE(55)
    sat_longitude = 9.5_wp * degrad
  CASE(56,57)
    sat_longitude = 0.0_wp * degrad
  CASE default
    ierror = 1
    yerror = 'ERROR *** Satellite position unknown ***'
  END SELECT

END FUNCTION sat_longitude

!=============================================================================
!=============================================================================

!------------------------------------------------------------------------------
! End of module procedure organize_satellites
!------------------------------------------------------------------------------

END SUBROUTINE organize_satellites
