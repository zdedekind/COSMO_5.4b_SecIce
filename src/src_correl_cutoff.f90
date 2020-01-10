!+ Source module for horiz. correlations and areas of influence for the nudging
!-------------------------------------------------------------------------------

MODULE src_correl_cutoff

!-------------------------------------------------------------------------------
!
! Description:
!   The module "CORREL_CUTOFF" computes the scaled horizontal distances,
!   geometrical factors of the 2-dimensional horizontal wind correlations,
!   and the areas of influence of observations depending on the cutoffs.
!   This module contains only the following procedure:
!   - cutoff_wind_correl
!   This is called by procedures from the modules "src_mult_spread" and
!   "src_sing_spread", which organize the spreading of the observational
!   information.
!   Note: This module does not contain any communication between PE's.
!   08.12.97
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2734
!  fax:    +49  69  8236 1493
!  email:  schraff@fe1.dwd.d400.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.13       1998/10/22 Christoph Schraff
!  Initial release
! 1.19       1998/12/11 Christoph Schraff
!  Logical 'IF' expressions containing arrays specifically for multi-level resp.
!  single-level data are used only after checking for 'nactio' (array bounds !).
! 2.4        2001/01/29 Christoph Schraff
!  Reorganisation of loops to enhance efficiency on vector processors.
! 3.6        2003/12/11 Christoph Schraff
!  Printout for SATOB data. Bug correction for print-out.
! 3.12       2004/09/15 Christoph Schraff
!  Another small bug fix for print-out to prevent array bound violations.
! 3.18       2006/03/03 Christoph Schraff
!  Bug correction at determination of 'lsigma' for model layer 'ktp'.
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_22        2012/01/31 Christoph Schraff
!  - Optimisations for decreased cost on NEC-SX9 (particularly in Part B of
!    sections 1 and 2: loop splitting, or 1-D loops where beneficial).
!  - Some variables moved from 'data_nudge_all' to 'data_obs_lib_cosmo'.
!  - Subroutine arguments list extended by vertical model level 'k'.
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!===============================================================================
!
! Declarations:
!
! Modules used:
!
!-------------------------------------------------------------------------------

USE data_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!-------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------

!   ie,           & ! number of grid points in zonal direction
!   je,           & ! number of grid points in meridional direction
!   ke,           & ! number of grid points in vertical direction
    ie_tot,       & ! number of grid points in zonal direction
    je_tot,       & ! number of grid points in meridional direction

! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------

    istart,       & ! start index for the forecast of w, t, qd, qw and pp
    iend,         & ! end index for the forecast of w, t, qd, qw and pp

!   meridional direction
    jstart,       & ! start index for the forecast of w, t, qd, qw and pp
    jend            ! end index for the forecast of w, t, qd, qw and pp

! end of data_modelconfig

!-------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------

    nstop,        & ! last time step of the forecast
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels

! 5. additional control variables
! -------------------------------

    lreproduce      ! the results are reproducible in parallel mode

! end of data_runcontrol 

!-------------------------------------------------------------------------------

USE data_nudge_all , ONLY :   &

! 1. parameters and related variables
! -----------------------------------

    lwonl        ,& ! .TRUE if grid pt. (ionl ,jonl ) lies in the local domain
    lwonl2       ,& ! .TRUE if grid pt. (ionl2,jonl2) lies in the local domain
    lfirst       ,& ! .TRUE if 'organize_nudging' is called the first time

! 2. namelist variables controlling the data assimilation
! -------------------------------------------------------

    tconbox      ,& ! 6*dt: timestep [s] for computing analysis increments
                    !       (i.e. time box of constant analysis increments)
    msprpar      ,& ! 1   : switch specifying the surfaces along which observat.
                    !       increments of upper-air data are (primarily) spread 
    msprpsu      ,& ! 0   : switch specifying the surface along which surface-  
                    !       level data increments are primarily spreaded
    cutofr       ,& ! 4* 3.5 : cut-off in 'COSAC' units of the horizontal
                    !          correlation function
    cutofsu      ,& ! 2., 3.5: cut-off in 'COSAC' units of the horizontal
                    ! 2., 2.   correlation function
    ionl         ,& ! 167    : / grid point coordinates
    jonl         ,& ! 103    : \ for standard output on nudging
    ionl2        ,& ! 167    : / 2nd grid pt coordinates
    jonl2           ! 103    : \ for other standard output on nudging

! end of data_nudge_all

!-------------------------------------------------------------------------------

USE data_obs_lib_cosmo , ONLY :   &

! 1. General parameters
! ---------------------

    c0         ,& ! standard real constant 0.0
    c1         ,& ! standard real constant 1.0
    c2         ,& ! standard real constant 2.0
    c05        ,& ! standard real constant 0.5
    c3600      ,& ! standard real constant 3600.0
    rmdi       ,& ! =-1.E31_wp : commonly used missing data indicator
    rmdich     ,& ! =-1.E30_wp : commonly used check value for miss data
    epsy       ,& ! = 1.E-8_wp : commonly used very small value > 0
    i1         ,& ! standard integer constant 1

! 4. I/O device numbers for obs processing / nudging and file names
! -----------------------------------------------------------------

    nupr       ,& ! unit number of file for all the remaining information

! 5. CMA observation type and code type numbers
! ---------------------------------------------

    nairep     ,& ! AIREP reports (all aircraft reports)
    nsatob        ! SATOB reports

! end of data_obs_lib_cosmo

!-------------------------------------------------------------------------------

USE data_nudge_gather , ONLY :   &

! 1. Parameters and general variables
! -----------------------------------

    ystid        ,& ! obs. station identity to be printed
    zalllow      ,& ! smallest allowed height within model domain
    zallhig      ,& ! largest  allowed height within model domain
    ktp          ,& ! lowermost purely horizontal model main level
    ktth         ,& ! top model level with spreading along isentropic surfaces

! 2. The required local information on the observations and their location
! ------------------------------------------------------------------------

! local information on multi-level reports
    zvcutml      ,& ! vertical cut-off at the base / top of the profile
    kviflml      ,& ! vertical range of possibly influenced model levels
    zspobml      ,& ! spreading parameter at the base / top of the profile
    zsprml       ,& ! spreading parameter at model levels & obs. location
    zriflml      ,& ! upper estimate to horizontal radius of influence
    kobtyml      ,& ! observation type
    ltiml        ,& ! .TRUE if temporal linear interpol. for multi-level rep.
! local information on upper-air single-level reports
    zvcutua      ,& ! vertical cut-off below / above the obs.
    kviflua      ,& ! vertical range of possibly influenced model levels
    zsprob       ,& ! spreading parameter at the obs. point
    zsprua       ,& ! spreading parameter at model levels & obs. location
    zriflua      ,& ! upper estimate to horizontal radius of influence
! local information on surface-level reports
    zvcutsu      ,& ! vertical cut-off above the obs.
    kviflsu      ,& ! vertical range of possibly influenced model levels
    zsposu       ,& ! spreading parameter at the obs. point
    zsprsu          ! spreading parameter at model levels & obs. location

! end of data_nudge_gather

!-------------------------------------------------------------------------------

USE data_nudge_spread , ONLY :   &

! 5. Geometrical fields used for horizontal distances and wind correlations
! -------------------------------------------------------------------------

    yyd          ,& ! latitudinal (merid.) distance on tangent cone projection
    yyd2         ,& ! yyd **2
    xxd2         ,& ! square of zonal distance on tangent cone projection
    c2alpa       ,& !  / further factors used to compute
    scalpa       ,& ! /  the 2-dim. horizontal wind correlations
    xcalpa       ,& ! \  and the horizontal distances
    xsalpa       ,& !  \ on the tangent cone projection

! 6. Further horizontal input fields for the spreading of obs. increments
! -----------------------------------------------------------------------

    zspr         ,& ! spreading parameter , param. def. non-isotropic weights
    zdds         ,& ! scaled horizontal distance betw. obs. and target grid pt
    zcoruu       ,& ! zonal  wind - zonal  wind correlation  \  (without
    zcoruv       ,& ! zonal  wind - merid. wind correlation   \  EXP( -zdds )
    zcorvu       ,& ! merid. wind - zonal  wind correlation   /  -term )
    zcorvv       ,& ! merid. wind - merid. wind correlation  /

! 7. Further fields and variables used for or during the spreading
! ----------------------------------------------------------------

    gppkmi       ,& ! convertor for zonal dist. from 'km' to rotated grid pts.
    gppkmj       ,& ! convertor for merid. dist. from 'km' to rotated grid pts
    rtinfl       ,& ! horizontal correlation scale
    lcutof       ,& ! .TRUE if grid pt. is within area of influence of report
    icutof       ,& ! x-coordinate of grid pt. within area of influence
    jcutof       ,& ! x-coordinate of grid pt. within area of influence
    icutmp       ,& ! x-coordinate of grid pt. within area of influence, tmp
    jcutmp       ,& ! x-coordinate of grid pt. within area of influence, tmp
    ncutof       ,& ! number of grid points within area of influence
    idimcut      ,& ! dimension of 'icutof', 'jcutof', etc.

! 8. Indices
! ----------

    io    ,jo    ,& ! local  indices of location of observation
    io_tot,jo_tot,& ! global indices of location of observation
    istaspr      ,& ! /  lower left corner of domain
    jstaspr      ,& ! \  containing area of influence
    iendspr      ,& ! /  upper right corner of domain
    jendspr      ,& ! \  containing area of influence
    jrs_tot      ,& ! /  index range for convertor
    jre_tot      ,& ! \  for zonal distances 'gppkmi'
    ista            ! index of observing station

! end of data_nudge_spread

!-------------------------------------------------------------------------------

IMPLICIT NONE

!===============================================================================

!-------------------------------------------------------------------------------
! Public and Private Subroutines
!-------------------------------------------------------------------------------

CONTAINS

!+ Module procedure cutoff_wind_correl in "CORREL_CUTOFF" for computing
!  horizontal correlations and areas of influence of observations
!-------------------------------------------------------------------------------

SUBROUTINE cutoff_wind_correl ( nactio , k , fcnodiv , ivarex )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "CORREL_CUTOFF" determines the scaled hori-
!   zontal distance to the observation location, the geometrical factors of the
!   2-dim horizontal wind correlations, and the areas of influence.
!   Optionally, it only determines the location of the highest or lowest grid
!   point within the horizontal area of influence, or checks if grid points
!   exist within this area at all.
!
! Method:
!   Horizontal distances and the 2-dimenstional horizontal wind correlations
!   are computed by using the tangent cone projection, following the ideas of
!   Bell et al. (1996) for the UKMO Analysis Correction Scheme.
!   The area of influence is determined by using the above horizontal distances,
!   by computing the vertical distances to the top / base of the obs. reports,
!   and by subsequent checking against specified cut-off radii.
!   (For lateral spreading of observational information along model levels, the
!    vertical check, which needs to be done only at the obs. location in this
!    case, has been done before.)
!   (For lateral spreading along isentropes, different vertical cut-off radii
!    are used in the troposphere and the stratosphere.)
!
!
! Current Code Owner:  Christoph Schraff
!
! S-story:
! Version   Date       Comment
! -------   ----       -------
! 1.0       24.07.97   Original code.    Christoph Schraff
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE


! Subroutine arguments:
! --------------------

  INTEGER (KIND=iintegers), INTENT (IN)         ::       &
    k     ,& ! index of current vertical model level
    nactio   ! = 1  : for multi-level reports at the same location
             ! = 2  : for upper-air single-level data
             ! = 3  : for surface-level data
             ! =4,5 : for determining only the location of the highest (4) or
             !        lowest (5) grid pt. within the horiz. area of influence
             ! = 6  : for checking if grid pts. exist in horiz. area of infl.
             ! = 0  : no action

  REAL    (KIND=wp   ),     INTENT (IN)         ::       &
    fcnodiv  ! weight to non-divergence correction in the 2-dim wind correlation

  INTEGER (KIND=iintegers), INTENT (INOUT)      ::       &
    ivarex   ! as input : index of variable: 1:u,v; 2:T; 3:RH
             ! as output: (nactio = 6   : if grid pts. exist in area of infl.
             !                            for any (0,1) variable: 1, else: 0
             !            (nactio = 4,5): location (compressed coord.) of
             !                            lowest / highest grid pt. within
             !                            area of influence

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    i      , j       ,& ! loop indices in horizontal direction
!   io     , jo      ,& ! local indices of location of observation
!   istaspr, jstaspr ,& ! lower left corner of domain containing area of infl.
!   iendspr, jendspr ,& ! upper right corner of domain containing area of infl.
    ilenspr, jlenspr ,& ! size of 2-dimensional domain containing area of infl.
    idix   , jdix    ,& ! indices depending on the horiz. distance to the obs.
    izext  , jzext   ,& ! coord. of lowest or highest grid pt. within infl. area
!   ista             ,& ! index of observing station
    ivar             ,& ! index of spreaded quantity: 1=(u,v); 3=T; 4=RH
    ivrs             ,& ! index of spreaded quantity: 1=(u,v); 2=T; 3=RH
    msi              ,& ! index of spreading parameter
    mvifl1 , mvifl2  ,& ! indicators if report may influence current model level
    irange , jrange  ,& ! zonal , meridional 'radii' of the domain which
                        ! contains maximum possible the area of influence
    istamx , jstamx  ,& ! lower left corner of domain containing this area
    iendmx , jendmx  ,& ! upper right corner of domain containing this area
    ici    , ici2    ,& ! index counter over influenced grid points
    ixc    , ii      ,& ! loop index over influenced grid points
    j_shft           ,& ! number of j-lines before current j-line
    kk               ,& ! loop indices in vertical   direction
    ncutmp           ,& ! number of grid points within area of influence, tmp
    mijtot           ,& ! # grid pts. within area of influence
    ijstot           ,& ! # grid pts. within domain containing area of influence
    nsign               ! sign depending if check against lower or upper cut-off

  REAL    (KIND=wp   )     ::  &
    zrtinflr         ,& ! reciproke of horizontal correlation scale
    zrtinfl2         ,& ! square    of horizontal correlation scale
    zrcutof          ,& ! horizontal cut-off radius [km]
    zrcuto2          ,& ! square of horizontal cut-off radius [km^2]
    zdist2           ,& ! square of horizontal distance to obs [km^2]
    zcuta            ,& ! upper vertical cut-off (units of spreading parameter)
    zcutb            ,& ! lower vertical cut-off (units of spreading parameter)
    zcuta1 , zcutb1  ,& ! upper , lower cut-off for 1st multi-level report
    zcuta2 , zcutb2  ,& ! upper , lower cut-off for 2nd multi-level report
    rnondiv          ,& ! auxiliary quantity used for 2-dim wind correlations
    zmargin          ,& ! pot. temperature margin accounting for instabilities,
                        ! when excluding influence to model levels further above
    zsob             ,& ! value of spreading parameter at top obs. of type ivar
    zhext            ,& ! extreme height of orography within area of influence
    zorosig             ! scaled orography

  LOGICAL                  ::  &
    lsigma           ,& ! spreading along model levels
    ladj                ! range of influenced model levels to be adjusted

  CHARACTER (LEN= 3) yobsys

  LOGICAL                  ::  &
    l_aux (iendspr*jendspr)  ! auxilliary array: 'lcutof' as 1-dim array

! Local (automatic) arrays: None
! -------------------------
!
!------------ End of header ----------------------------------------------------


 
!-------------------------------------------------------------------------------
! Begin Subroutine cutoff_wind_correl
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1: Upper-air (multi-level & single-level) data
!-------------------------------------------------------------------------------
 
  IF ((nactio == 1) .OR. (nactio == 2)) THEN

! Part A: get horizontal cutoff radius, & vertical cutoffs (in z or theta units)
! ------------------------------------------------------------------------------

    ivar = ivarex
    ivrs = MAX( ivar-1 , 1 )
    msi  = MAX( msprpar , 1_iintegers )
! lsigma = .true. if the MAIN level 'k' is a sigma level
    lsigma  = (msprpar == 0)  .OR.  ((msprpar == 1) .AND. (k <= ktp-1))        &
                              .OR.  ((msprpar == 2) .AND. (k < ktth-1))
    zrtinflr= c1 / rtinfl
    zrtinfl2= rtinfl * rtinfl
    zrcutof = rtinfl * cutofr(ivar)
    zrcuto2 = zrcutof **2
    mvifl1  = 0
    mvifl2  = 0
    zcuta1  = c0
    zcuta2  = c0
    IF (lsigma) THEN
      zcutb =   rmdi
      zcuta = - rmdi
    ELSEIF (nactio == 1) THEN
! 'mvifl?' == 1 if report may influence current model level, else 'mvifl?' == 0
      mvifl1 = MIN( 1 ,  MAX(  kviflml(ista,1,ivrs) - k + 1 , 0 )              &
                       * MAX( -kviflml(ista,2,ivrs) + k + 1 , 0 ) )
      mvifl2 = MIN( 1 ,  MAX(  kviflml(ista,3,ivrs) - k + 1 , 0 )              &
                       * MAX( -kviflml(ista,4,ivrs) + k + 1 , 0 ) )
! If both reports may influence current model level, then 'zcutb' is the minimum
! of the 2 lower cut-offs, else 'zcutb' = lower cut-off of the relevant report.
      zcutb = MIN( mvifl1 *(zvcutml(ista,1,ivrs) - zallhig)                    &
                 , mvifl2 *(zvcutml(ista,3,ivrs) - zallhig) ) + zallhig
      zcuta = MAX( mvifl1 *(zvcutml(ista,2,ivrs) - zalllow)                    &
                 , mvifl2 *(zvcutml(ista,4,ivrs) - zalllow) ) + zalllow
    ELSE
      zcutb = zvcutua(ista,1,ivrs)
      zcuta = zvcutua(ista,2,ivrs)
    ENDIF
! Preparation for parts C & D
    IF (mvifl1*mvifl2 == 1) THEN
      zcuta1 = zvcutml(ista,2,ivrs)
      zcuta2 = zvcutml(ista,4,ivrs)
    ENDIF

! Part B: get area of influence, horizontal distances, geometrical factors
!         to the 2-dim wind correlation.
! ------------------------------------------------------------------------

!   get area of influence and horizontal distance
!   ---------------------------------------------

    ! indirect addressing is moved from the following loop
    ! to a separate one because this runs faster on NEC-SX9
    DO j   = jstaspr , jendspr
!CDIR NODEP
!CDIR ON_ADB(lcutof)
      DO i = istaspr , iendspr
        lcutof (i,j,1) = .FALSE.
        lcutof (i,j,2) = .FALSE.
        IF ((zspr(i,j,msi) < zcuta) .AND. (zspr(i,j,msi) > zcutb) ) THEN
          idix = i - io + ie_tot
          jdix = j - jo + je_tot
          zdist2 = xxd2(idix,j) + yyd2(jdix) - yyd(jdix) *xsalpa(idix,j)
          IF (zdist2 < zrcuto2) THEN
            lcutof (i,j,1) = .TRUE.
            zdds (i,j,ivrs) =   SQRT( zdist2 ) * zrtinflr
          ENDIF
        ENDIF
      ENDDO
    ENDDO

!   determine indirect addressing
!   -----------------------------

    ! map from 2-d array 'lcutof' to 1-d array 'l_aux'
    ixc = 0
    DO j   = jstaspr , jendspr
!CDIR ON_ADB(lcutof)
!CDIR ON_ADB(l_aux)
      DO i = istaspr , iendspr
        ixc = ixc + 1
        l_aux (ixc) = lcutof(i,j,1)
      ENDDO
    ENDDO

    ! use 1-d array 'l_aux' in 1-d loop for indirect addressing
    ilenspr = iendspr - istaspr + 1
    jlenspr = jendspr - jstaspr + 1
    ici = 0
!CDIR ON_ADB(l_aux)
    DO ixc = 1 , jlenspr *ilenspr
      j = (ixc-1) / ilenspr + jstaspr
      i =  ixc-1            + istaspr - (j-jstaspr)*ilenspr
      IF (l_aux(ixc)) THEN
        ici = ici + 1
        icutof (ici,1) = i
        jcutof (ici,1) = j
      ENDIF
    ENDDO
! the above loop replaces the following one because it runs faster on NEC-SX9
!   DO j   = jstaspr , jendspr
!     DO i = istaspr , iendspr
!       IF (lcutof(i,j,1)) THEN
!         ici = ici + 1
!         icutof (ici,1) = i
!         jcutof (ici,1) = j
!       ENDIF
!     ENDDO
!   ENDDO
    ncutof (1) = ici

!   get geometrical factors to the 2-dim wind correlation
!   -----------------------------------------------------

    IF (ivrs == 1) THEN

!!    DO ixc = 1 , jlenspr *ilenspr
!!      j = (ixc-1) / ilenspr + jstaspr
!!      i = ixc - (j-jstaspr)*ilenspr + istaspr - 1
!! (on NEC: the above 1-d set up of this loop takes much longer)
      DO j   = jstaspr , jendspr
!CDIR NODEP
        DO i = istaspr , iendspr
          IF (lcutof(i,j,1)) THEN
            idix = i - io + ie_tot
            jdix = j - jo + je_tot
            zdist2      = zdds(i,j,ivrs) *zdds(i,j,ivrs) *zrtinfl2
            rnondiv     = fcnodiv * zdds(i,j,ivrs) / MAX( zdist2 , epsy )
            zcoruu(i,j) =  c2 * c2alpa(idix,j) - c1                          &
                         - rnondiv *( (zdist2 +yyd2(jdix)) *c2alpa(idix,j)   &
                                     - zdist2 )
            zcorvv(i,j) =  c2 * c2alpa(idix,j) - c1                          &
                         - rnondiv *( (zdist2 -yyd2(jdix)) *c2alpa(idix,j) )
            zcoruv(i,j) =  c2 * scalpa(idix,j)                               &
                         - rnondiv *( (zdist2 +yyd2(jdix)) *scalpa(idix,j)   &
                                     -      yyd(jdix)      *xcalpa(idix,j) )
            zcorvu(i,j) =- c2 * scalpa(idix,j)                               &
                         + rnondiv *( (zdist2 -yyd2(jdix)) *scalpa(idix,j)   &
                                     +      yyd(jdix)      *xcalpa(idix,j) )
          ENDIF
        ENDDO
      ENDDO

    ENDIF


! Part C: adjust range of influenced model levels, if (a) spreading is along
!         isentropes, (b) the potential temperature on present level exceeds
!         the cut-off within the maximum possible horizontal area of influence,
!         (c) reproducibility using different numbers of nodes is not required
! ----------------------------------------------------------------------------

    ladj  =  (msprpar == 2) .AND. (.NOT. lreproduce)
    IF ((ladj) .AND. (nactio == 1)) ladj =      (kviflml(ista,2,ivrs) == ktth) &
                                           .OR. (kviflml(ista,4,ivrs) == ktth)
    IF ((ladj) .AND. (nactio == 2)) ladj =      (kviflua(ista,2,ivrs) == ktth)
    IF  (ladj) THEN
      IF (mvifl1*mvifl2 == 1) zcuta = MIN( zcuta1 , zcuta2 )
      IF (nactio == 1) zmargin = MAX( c0,         zsprml(ista,k)               &
                                         -MINVAL( zsprml(ista,k:ktth) ))
      IF (nactio == 2) zmargin = MAX( c0,         zsprua(ista,k)               &
                                         -MINVAL( zsprua(ista,k:ktth) ))
! 4K margin: instabilities (in pot. temperature) of up to 4K are safely treated
      zmargin = zmargin + 4.0_wp
      IF (nactio == 1) zrcutof = zriflml(ista,2,ivrs) * cutofr(ivar)
      IF (nactio == 2) zrcutof = zriflua(ista,2,ivrs) * cutofr(ivar)
      jrange  = INT( zrcutof * gppkmj )
      jstamx  = MAX( INT( jo - jrange, iintegers ) , jrs_tot )
      jendmx  = MIN( INT( jo + jrange, iintegers ) , jre_tot )
      irange  = INT( zrcutof * MAX( gppkmi(jstamx) , gppkmi(jendmx) ) )
      jstamx  = MAX( INT( jo - jrange, iintegers ) , jstart  )
      jendmx  = MIN( INT( jo + jrange, iintegers ) , jend    )
      istamx  = MAX( INT( io - irange, iintegers ) , istart  )
      iendmx  = MIN( INT( io + irange, iintegers ) , iend    )
      ici = 0
      loopj: DO   j = jstamx , jendmx
        loopi: DO i = istamx , iendmx
          IF (zspr(i,j,msi)-zmargin < zcuta) ici = ici + 1
        ENDDO loopi
        IF ( ici > 0 ) EXIT loopj
      ENDDO loopj
      IF (mvifl1*mvifl2 == 1) zcuta = MAX( zcuta1 , zcuta2 )
      IF ((ici == 0) .AND. (nactio == 2)) kviflua(ista,2,ivrs) = k + 1
      IF ((ici == 0) .AND. (nactio == 1)) THEN
        IF (      ((mvifl2 == 0) .OR. (zcuta1 < zcuta2+epsy))                  &
            .AND. (kviflml(ista,2,ivrs) == ktth)) kviflml(ista,2,ivrs) = k + 1
        IF (      ((mvifl1 == 0) .OR. (zcuta2 < zcuta1+epsy))                  &
            .AND. (kviflml(ista,4,ivrs) == ktth)) kviflml(ista,4,ivrs) = k + 1
      ENDIF
    ENDIF

! Part D: If there are 2 multi-level reports from the current obs. station
!         influencing the current model level (and time), then the above
!         quantities have been computed for the area which is influenced
!         by the one OR the other report. Subsequently, the area of
!         influence is computed for each report separately.
! ------------------------------------------------------------------------

    IF (mvifl1*mvifl2 == 1) THEN
! 'mvifl?' == 1 if report may influence current model level, else 'mvifl?' == 0
! (Note: the lower cut-off is left unaltered after the first check in part A)
      mvifl1 = MIN( 1 , MAX( -kviflml(ista,2,ivrs) + k + 1 , 0 ) )
      mvifl2 = MIN( 1 , MAX( -kviflml(ista,4,ivrs) + k + 1 , 0 ) )
      zcutb1 = zvcutml(ista,1,ivrs)
      zcuta1 = zvcutml(ista,2,ivrs)
      zcutb2 = zvcutml(ista,3,ivrs)
      zcuta2 = zvcutml(ista,4,ivrs)
      IF (mvifl1*mvifl2 == 1) THEN
        DO j   = jstaspr , jendspr
!CDIR NODEP
          DO i = istaspr , iendspr
            lcutof (i,j,1) = .FALSE.
            lcutof (i,j,2) = .FALSE.
          ENDDO
        ENDDO
        ncutmp = ncutof(1)
        icutmp (1:ncutmp) = icutof(1:ncutmp,1)
        jcutmp (1:ncutmp) = jcutof(1:ncutmp,1)
        ici = 0
        ici2 = 0
!CDIR NODEP
        DO ixc = 1, ncutmp
          i = icutmp(ixc)
          j = jcutmp(ixc)
          IF ((zspr(i,j,msi) < zcuta1) .AND. (zspr(i,j,msi) > zcutb1)) THEN
            ici = ici + 1
            icutof (ici,1) = i
            jcutof (ici,1) = j
! set 'lcutof' (for linear temporal interpolation and printout)
            lcutof (i,j,1) = .TRUE.
          ENDIF
          IF ((zspr(i,j,msi) < zcuta2) .AND. (zspr(i,j,msi) > zcutb2)) THEN
            ici2 = ici2 + 1
            icutof (ici2,2) = i
            jcutof (ici2,2) = j
            lcutof (i,j ,2) = .TRUE.
          ENDIF
        ENDDO
        ncutof (1) = ici
        ncutof (2) = ici2

! If one (but only one) of the 2 reports have been cancelled in part C, and if
! the lower cut-off of the remaining report (which must have the higher upper
! cut-off) is above the lower cut-off of the cancelled report, then the area
! of influence has to be recomputed
      ELSEIF (      ((zcutb1-zcutb2)*(zcuta1-zcuta2) > epsy)                   &
              .AND. (mvifl1 + mvifl2 == 1)) THEN
        DO j   = jstaspr , jendspr
!CDIR NODEP
          DO i = istaspr , iendspr
            lcutof (i,j,1) = .FALSE.
          ENDDO
        ENDDO
!   if the 2nd report remains, re-set 'zcuta1', 'zcutb1'
        IF (mvifl2 == 1) THEN
          zcutb1 = zvcutml(ista,3,ivrs)
          zcuta1 = zvcutml(ista,4,ivrs)
        ENDIF
        ncutmp = ncutof(1)
        icutmp (1:ncutmp) = icutof(1:ncutmp,1)
        jcutmp (1:ncutmp) = jcutof(1:ncutmp,1)
        ici = 0
!CDIR NODEP
        DO ixc = 1, ncutmp
          i = icutmp(ixc)
          j = jcutmp(ixc)
          IF ((zspr(i,j,msi) < zcuta1) .AND. (zspr(i,j,msi) > zcutb1)) THEN
            ici = ici + 1
            icutof (ici,1) = i
            jcutof (ici,1) = j
            lcutof (i,j,1) = .TRUE.
          ENDIF
        ENDDO
        ncutof (1) = ici
      ENDIF
    ENDIF

! Part E: printout for control
! ----------------------------

    IF ((ntstep <= 2) .AND. (lwonl)) THEN
      IF (     ((io == ionl ) .AND. (jo == jonl ))                             &
          .OR. ((io == ionl2) .AND. (jo == jonl2)) .OR. (ista == 1)) THEN
        ijstot = (iendspr - istaspr + 1) * (jendspr - jstaspr + 1)
        yobsys = 'air'
        IF (nactio == 2) zsob = zsprob (ista)
        IF (nactio == 1) THEN
          zsob  = MAX( mvifl1 * zspobml(ista,2,ivrs)                           &
                     , mvifl2 * zspobml(ista,4,ivrs) )
          IF (lsigma)  zsob = MAX( zspobml(ista,2,ivrs) , zspobml(ista,4,ivrs) )
          IF (kobtyml(ista) /= nairep)  yobsys = 'ras'
        ELSEIF (nactio == 2) THEN
          zsob  = zsprob (ista)
          IF (kobtyml(ista) == nsatob)  yobsys = 'amv'
        ENDIF
        zrcutof = rtinfl * cutofr(ivar)
        zhext   = c0
        IF (      (io >= istart) .AND. (io <= iend)                            &
            .AND. (jo >= jstart) .AND. (jo <= jend))   zhext = zspr(io,jo,msi)
        IF (msi == 1) THEN
          WRITE( nupr,'(A3,''-cor '',A ,I2,I3,2I4,F6.0,'',spr-par: ob,cutof''  &
                      &,3F7.0,'',at k'',F7.0,2I5)' )    yobsys, ystid          &
               , ivar, k, io_tot, jo_tot, zrcutof, zsob                        &
               , zcutb, zcuta, zhext, ijstot, ncutof(1)
        ELSE
          WRITE( nupr,'(A3,''-cor '',A ,I2,I3,2I4,F6.0,'',spr-par: ob,cutof''  &
                      &,3F7.2,'',at k'',F7.2,2I5)' )    yobsys, ystid          &
               , ivar, k, io_tot, jo_tot, zrcutof, zsob                        &
               , zcutb, zcuta, zhext, ijstot, ncutof(1)
        ENDIF
      ENDIF
    ENDIF


!-------------------------------------------------------------------------------
!  Section 2: Surface-level data
!-------------------------------------------------------------------------------
 
  ELSEIF (nactio == 3) THEN

! Part A: get horizontal cutoff radius, & vertical cutoffs (in z or theta units)
! ------------------------------------------------------------------------------

    ivar = ivarex
    ivrs = MAX( ivar-1 , 1 )
    msi  = MAX( msprpsu , 1_iintegers )
    lsigma =  (msprpsu == 0)  .OR.  ((msprpsu == 1) .AND. (k <= ktp-1))        &
                              .OR.  ((msprpsu == 2) .AND. (k < ktth-1))
    zrtinflr= c1 / rtinfl
    zrtinfl2= rtinfl * rtinfl
    zrcutof = rtinfl * cutofsu(ivar)
    zrcuto2 = zrcutof **2
    zcuta   = zvcutsu(ista,ivrs)
    IF (lsigma) zcuta = - rmdi

! Part B: get area of influence, horizontal distances, geometrical factors
!         to the 2-dim wind correlation.
! ------------------------------------------------------------------------

    IF ((msprpsu <= 1) .OR. (lreproduce)) THEN

      ! indirect addressing is moved from the following loop
      ! to a separate one because this runs faster on NEC-SX9
      DO j   = jstaspr , jendspr
!CDIR ON_ADB(lcutof)
        DO i = istaspr , iendspr
          lcutof (i,j,1) = .FALSE.
          IF (zspr(i,j,msi) < zcuta) THEN
            idix = i - io + ie_tot
            jdix = j - jo + je_tot
            zdist2 = xxd2(idix,j) + yyd2(jdix) - yyd(jdix) *xsalpa(idix,j)
            IF (zdist2 < zrcuto2) THEN
              lcutof (i,j,1) = .TRUE.
              zdds (i,j,ivrs) =  SQRT( zdist2 ) * zrtinflr
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      ! indirect addressing as in section 1
      ixc=0
      DO j   = jstaspr , jendspr
!CDIR ON_ADB(lcutof)
!CDIR ON_ADB(l_aux)
        DO i = istaspr , iendspr
          ixc = ixc + 1
          l_aux (ixc) = lcutof(i,j,1)
        ENDDO
      ENDDO
      ilenspr = iendspr - istaspr + 1
      jlenspr = jendspr - jstaspr + 1
      ici = 0
!CDIR ON_ADB(l_aux)
      DO ixc = 1 , jlenspr *ilenspr
        j = (ixc-1) / ilenspr + jstaspr
        i =  ixc-1            + istaspr - (j-jstaspr)*ilenspr
        IF (l_aux(ixc)) THEN
          ici = ici + 1
          icutof (ici,1) = i
          jcutof (ici,1) = j
        ENDIF
      ENDDO
      ncutof (1) = ici

      ! get geometrical factors to the 2-dim wind correlation
      IF (ivrs == 1) THEN
        DO j   = jstaspr , jendspr
          DO i = istaspr , iendspr
            IF (lcutof(i,j,1)) THEN
              idix = i - io + ie_tot
              jdix = j - jo + je_tot
              zdist2      =  zdds(i,j,ivrs) *zdds(i,j,ivrs) *zrtinfl2
              rnondiv     =  fcnodiv * zdds(i,j,ivrs) / MAX( zdist2 , epsy )
              zcoruu(i,j) =  c2 * c2alpa(idix,j) - c1                          &
                           - rnondiv *( (zdist2 +yyd2(jdix)) *c2alpa(idix,j)   &
                                       - zdist2 )
              zcorvv(i,j) =  c2 * c2alpa(idix,j) - c1                          &
                           - rnondiv *( (zdist2 -yyd2(jdix)) *c2alpa(idix,j) )
              zcoruv(i,j) =  c2 * scalpa(idix,j)                               &
                           - rnondiv *( (zdist2 +yyd2(jdix)) *scalpa(idix,j)   &
                                       -      yyd(jdix)      *xcalpa(idix,j) )
              zcorvu(i,j) =- c2 * scalpa(idix,j)                               &
                           + rnondiv *( (zdist2 -yyd2(jdix)) *scalpa(idix,j)   &
                                       +      yyd(jdix)      *xcalpa(idix,j) )
            ENDIF
          ENDDO
        ENDDO
      ENDIF
! if spreading parameter on the present model level exceeds the cut-off within
! the horizontal area of influence, the obs. is not assumed to influence any
! levels further above (if horizontal correlation scale is constant with height)
! ==> adjust vertical levels of influence
      IF ((ncutof(1) == 0) .AND. (msprpsu == 1))   kviflsu(ista,ivrs) = k + 1

    ELSE
! if (msprpsu == 2) and (.NOT. lreproduce) :
! different treatment of adjustment of vertical levels of influence 
      zmargin = MAX( c0,         zsprsu(ista,k)                                &
                        -MINVAL( zsprsu(ista,k:kviflsu(ista,ivrs)) ))
! 4K margin: instabilities (in pot. temperature) of up to 4K are safely treated
      zmargin = zmargin + 4.0_wp
      ici = 0
      ici2 = 0
      DO j   = jstaspr , jendspr
        DO i = istaspr , iendspr
          lcutof (i,j,1) = .FALSE.
          IF (zspr(i,j,msi)-zmargin < zcuta) THEN
            idix = i - io + ie_tot
            jdix = j - jo + je_tot
            zdist2 = xxd2(idix,j) + yyd2(jdix) - yyd(jdix) *xsalpa(idix,j)
            IF ((zdist2 < zrcuto2) .AND. (zspr(i,j,msi) < zcuta)) THEN
              ici = ici + 1
              icutof (ici,1) = i
              jcutof (ici,1) = j
              lcutof (i,j,1) = .TRUE.
              zdds (i,j,ivrs) =  SQRT( zdist2 ) * zrtinflr
              IF (ivrs == 1) THEN
                rnondiv     = fcnodiv * zdds(i,j,ivrs) / MAX( zdist2 , epsy )
                zcoruu(i,j) = c2 * c2alpa(idix,j) - c1                         &
                            - rnondiv *( (zdist2 + yyd2(jdix)) *c2alpa(idix,j) &
                                         - zdist2 )
                zcorvv(i,j) = c2 * c2alpa(idix,j) - c1                         &
                            - rnondiv *( (zdist2 - yyd2(jdix)) *c2alpa(idix,j) )
                zcoruv(i,j) = c2 * scalpa(idix,j)                              &
                            - rnondiv *( (zdist2 + yyd2(jdix)) *scalpa(idix,j) &
                                        -      yyd(jdix)       *xcalpa(idix,j) )
                zcorvu(i,j) =-c2 * scalpa(idix,j)                              &
                            + rnondiv *( (zdist2 - yyd2(jdix)) *scalpa(idix,j) &
                                        +      yyd(jdix)       *xcalpa(idix,j) )
              ENDIF
            ELSEIF (zdist2 < zrcuto2) THEN
              ici2 = ici2 + 1
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      ncutof (1) = ici
      IF (ici2 == 0)   kviflsu(ista,ivrs) = k + 1
    ENDIF

! Part C: printout for control
! ----------------------------

    IF ((ntstep <= 1) .AND. (lwonl2)) THEN
      IF (((io == ionl2) .AND. (jo == jonl2)) .OR. (ista == 1)) THEN
        ijstot = (iendspr - istaspr + 1) * (jendspr - jstaspr + 1)
        zhext   = c0
        IF (      (io >= istart) .AND. (io <= iend)                            &
            .AND. (jo >= jstart) .AND. (jo <= jend))   zhext = zspr(io,jo,msi)
        IF (msprpsu <= 1) THEN
          WRITE( nupr,'(''sfc-cor '',A ,I2,I3,2I4,F6.0,'',spr-par: ob,cutof''  &
                      &,2F7.0,'',at k'',F7.0,2I5)' )    ystid                  &
               , ivar, k, io_tot, jo_tot, zrcutof, zsposu(ista)                &
               , zcuta, zhext, ijstot, ncutof(1)
        ELSE
          WRITE( nupr,'(''sfc-cor '',A ,I2,I3,2I4,F6.0,'',spr-par: ob,cutof''  &
                      &,2F7.2,'',at k'',F7.2,2I5)' )    ystid                  &
               , ivar, k, io_tot, jo_tot, zrcutof, zsposu(ista)                &
               , zcuta, zhext, ijstot, ncutof(1)
        ENDIF
      ENDIF
    ENDIF


!-------------------------------------------------------------------------------
!  Section 3: Only: highest or lowest grid point within horiz. area of influence
!                   (called only if (k == ke)  --> hhl replaced by zhml=zspr(1)
!-------------------------------------------------------------------------------

  ELSEIF ((nactio == 4) .OR. (nactio == 5)) THEN

! nactio = 4 : highest grid pt. for finding lowest  influenced model level
! nactio = 5 : lowest  grid pt. for finding highest influenced model level
    ivar  = ivarex
    izext = 0
    jzext = 0
    nsign = 9 - 2 *nactio
    IF (nactio == 4) zhext =  zalllow
    IF (nactio == 5) zhext = -zallhig
! here, 'rtinfl' is the cut-off radius rather than the correlation scale
    zrcuto2 = rtinfl **2

    ! use of 1-D loop or of 2 separate loops for nactio == 4 resp. 5
    ! instead of the following 2-D loop is not faster on NEC-SX9
    DO j   = jstaspr , jendspr
      DO i = istaspr , iendspr
        idix = i - io + ie_tot
        jdix = j - jo + je_tot
        zdist2 = xxd2(idix,j) + yyd2(jdix) - yyd(jdix) *xsalpa(idix,j)
        IF (zdist2 < zrcuto2) THEN
          zorosig = nsign * zspr(i,j,1)
          IF (zorosig > zhext) THEN
            zhext = zorosig
            izext = i
            jzext = j
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    ivarex  =  10000* jzext  +  izext


!-------------------------------------------------------------------------------
!  Section 4: Only: check if grid pts. exist within horizontal area of influence
!-------------------------------------------------------------------------------

  ELSEIF (nactio == 6) THEN

    ivar    = ivarex
! here, 'rtinfl' is the cut-off radius rather than the correlation scale
    zrcuto2 = rtinfl **2
    zdist2  = zrcuto2 + c1

    loopjc: DO j   = jstaspr , jendspr
      ici = 0
      loopic: DO i = istaspr , iendspr
        idix = i - io + ie_tot
        jdix = j - jo + je_tot
        zdist2 = xxd2(idix,j) + yyd2(jdix) - yyd(jdix) *xsalpa(idix,j)
        IF (zdist2 < zrcuto2)  ici = ici + 1
      ENDDO loopic
      IF ( ici > 0 )  EXIT loopjc
    ENDDO loopjc
    ivarex = MIN( ici , 1 )

  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure cutoff_wind_correl
!-------------------------------------------------------------------------------

END SUBROUTINE cutoff_wind_correl


!-------------------------------------------------------------------------------

END MODULE src_correl_cutoff
