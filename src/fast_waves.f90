!+ External procedure for computing fast wave terms      
!------------------------------------------------------------------------------

SUBROUTINE fast_waves (nsmsteps, dts, xkd,                                  &
                       suten, svten, swten, stten, sppten,                  &
                       zubdt_west, zubdt_east, zvbdt_north, zvbdt_south,    &
                       ie, je, ke, kzflat )

!------------------------------------------------------------------------------
!
! Description:
!   This external procedure (for the use in src_leapfrog)
!   calculates the new values of the prognostic variables u,v,w,pp and T 
!   at time level n+1 (nnew).  
!
! Method:
!   The Leapfrog/split explicit time integration scheme is used: Within  
!   a Leapfrog intervall 2*dt the slow tendencies of the prognostic variables
!   due to advection and diffusion (which have been calculated before and
!   are stored in the corresponding tendency fieleds) are held constant
!   wheras the fast terms, which describe sound and gravity wave propagation,
!   are stepped forward by using a small timestep dts.
!   Horizontal wave propagation is treated explicitly, for vertical wave
!   propagation an implicit scheme is used.
!   In case of the two time-level integrations scheme, this routind is called
!   twice: first for the time interval dt/2 to obtain intermediate values of
!   the variables, and then for the interval dt to calculate the final values
!   for time level nnew.
!
! Current Code Owner: DWD, Michael Baldauf
!  phone:  +49  69  8062 2733
!  fax:    +49  69  8062 3721
!  email:  Michael.Baldauf@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.34       1999/12/10 Guenther Doms
!  Initial release
! 1.39       2000/05/03 Ulrich Schaettler
!  Adapted Calls to timing-routines and changed variable names related to
!  latitude and longitude. Introduced switch lw_freeslip for boundary
!  treatment of w.
! 2.8        2001/07/06 Ulrich Schaettler
!  Eliminated non-necessary variables from the USE-lists
! 2.18       2002/07/16 Ulrich Schaettler
!  Adaptations to use fast_waves also in the new 2 time level scheme.
!  (free slip lower boundary condition)
! 3.5        2003/09/02 Ulrich Schaettler
!  Adaptation of interface for exchg_boundaries 
! 3.7        2004/02/18 Ulrich Schaettler
!  Renamed cphi (crlat), acphir (acrlat)
! 3.18       2006/03/03 Ulrich Schaettler
!  Eliminated variables from the USE lists that were not necessary
!  Implemented use of a dynamical bottom boundary condition
!  (formulation from Almut Gassmann)
! V4_4         2008/07/16 Ulrich Schaettler
!  Adapted interface of get_timings
!  Introduced NEC Compiler Directives and re-ordered some loops for vector
! V4_9         2009/07/16 Ulrich Schaettler, Uwe Boehm (CLM, PIK)
!  Adapted end indices for variable xrhsy of dynamical bottom boundary condition
! V4_10        2009/09/11 Andy Dobler (CLM, Uni Fra)
!  Another adaptation of these indices
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_17        2011/02/24 Ulrich Blahak
!  Adapted interface of exchg_boundaries; corrected kzdims(1:20) -> kzdims(1:24);
!  eliminated my_peri_neigh
! V4_18        2011/05/26 Ulrich Schaettler
!  Changed the code owner
! V4_23        2012/05/10 Ulrich Schaettler
!  Removed src_2timelevel and modified comments accordingly
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer, Hans-Juergen Panitz
!  Replaced qx-variables by using them from the tracer module
!  Introduced nexch_tag for MPI boundary exchange tag to replace ntstep (HJP)
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  MESSy interface introduced
! V4_28        2013/07/12 Ulrich Schaettler
!  Pass kflat as parameter to avoid dependency on vgrid_refatm_utils
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
! V5_2         2015-05-21 Ulrich Schaettler
!  Reorganized numbering of MPI datatypes ID
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

USE data_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------
    ke1,          & ! ke+1

! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from
!    the other ones because of the use of the staggered Arakawa-B-grid.
!
!   zonal direction
    istart,       & ! start index for the forecast of w, t, qv, qc and pp
    iend,         & ! end index for the forecast of w, t, qv, qc and pp
    istartu,      & ! start index for the forecast of u
    iendu,        & ! end index for the forecast of u
    istartv,      & ! start index for the forecast of v
    iendv,        & ! end index for the forecast of v
    istartpar,    & ! start index for computations in the parallel program
    iendpar,      & ! end index for computations in the parallel program

!   meridional direction
    jstart,       & ! start index for the forecast of w, t, qv, qc and pp
    jend,         & ! end index for the forecast of w, t, qv, qc and pp
    jstartu,      & ! start index for the forecast of u
    jendu,        & ! start index for the forecast of u
    jstartv,      & ! start index for the forecast of v
    jendv,        & ! end index for the forecast of v
    jstartpar,    & ! start index for computations in the parallel program
    jendpar,      & ! end index for computations in the parallel program

! 4. constants for the horizontal rotated grid and related variables
! ------------------------------------------------------------------

    eddlon,       & ! 1 / dlon
    eddlat,       & ! 1 / dlat
    edadlat,      & ! 1 / (radius of the earth * dlat)

! 5. variables for the time discretization and related variables
! --------------------------------------------------------------

    dt,           & ! long time-step
    betasw,       & ! beta-variable for treatment of soundwaves

! 8. Organizational variables to handle the COSMO humidity tracers
! ----------------------------------------------------------------
    idt_qv, idt_qc

! end of data_modelconfig

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------
    rho0       ,    & ! reference density at the full model levels    (kg/m3)
    dp0        ,    & ! reference pressure thickness of layers        ( Pa)
    p0         ,    & ! reference pressure at main levels             ( Pa)
    hhl        ,    & ! geometical height of half model levels        ( m )

! 2. external parameter fields                                        (unit)
! ----------------------------
    crlat      ,    & ! cosine of transformed latitude
    acrlat     ,    & ! 1 / ( crlat * radius of the earth )           ( 1/m )

! 3. prognostic variables                                             (unit)
! -----------------------
    u          ,    & ! zonal wind speed                              ( m/s )
    v          ,    & ! meridional wind speed                         ( m/s )
    w          ,    & ! vertical wind speed (defined on half levels)  ( m/s )
    t          ,    & ! temperature                                   (  k  )
    pp         ,    & ! deviation from the reference pressure         ( pa  )

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------
    rho               ! density of moist air

! end of data_fields

!------------------------------------------------------------------------------

USE data_parallel,      ONLY :  &
    num_compute,     & ! number of compute PEs
    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains
    ldatatypes,      & ! if .TRUE.: use MPI-Datatypes for some communications
    ltime_barrier,   & ! if .TRUE.: use additional barriers for determining the
                       ! load-imbalance
    ncomm_type,      & ! type of communication
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    my_cart_neigh,   & ! neighbors of this subdomain in the cartesian grid
    icomm_cart,      & ! communicator for the virtual cartesian topology
    iexch_req,       & ! stores the sends requests for the neighbor-exchange
                       ! that can be used by MPI_WAIT to identify the send
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    nexch_tag,       & ! tag to be used for MPI boundary exchange
                       !  (in calls to exchg_boundaries)
    sendbuf,         & ! sending buffer for boundary exchange:
                       ! 1-4 are used for sending, 5-8 are used for receiving
    isendbuflen        ! length of one column of sendbuf

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
    nold,         & ! corresponds to ntstep - 1
    nnow,         & ! corresponds to ntstep
    nnew,         & ! corresponds to ntstep + 1

! 7. additional control variables
! -------------------------------
    ldyn_bbc,     & ! if .TRUE., dynamical bottom boundary condition
    llm,          & ! if .TRUE., running with a lowered upper boundary
    l2tls,        & ! forecast with 2-TL integration scheme
    lperi_x,      & ! if lgen=.TRUE.: periodic boundary conditions (.TRUE.) in x-dir.
                    !                 or with Davies conditions (.FALSE.)
    lperi_y,      & ! if lgen=.TRUE.: periodic boundary conditions (.TRUE.) in y-dir.
                    !                 or with Davies conditions (.FALSE.)
    l2dim,        & ! if lgen=.TRUE.: 2dimensional model version (.TRUE) or
                    !                 full 3dimensional version (.FALSE)
    ltime,        & ! detailled timings of the program are given
    lw_freeslip     ! if .TRUE.: with free slip lateral boundary condition and
                    ! if .FALSE. specified lateral boundary values for w

!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 2. physical constants and related variables
! -------------------------------------------
    r_d,          & ! gas constant for dry air
    cpdr,         & ! 1 / cp_d
    gamma,        & ! 1 / (1 - rdocp)   ( = cp_d/cv_d)
    g               ! acceleration due to gravity

!------------------------------------------------------------------------------

USE src_tracer,       ONLY: trcr_get, trcr_errorstr

!------------------------------------------------------------------------------

USE environment     , ONLY :  exchg_boundaries, comm_barrier, model_abort
USE time_utilities  , ONLY :  get_timings, i_fast_waves, i_fast_waves_comm,  &
                              i_fast_waves_barrier, i_barrier_waiting_dyn

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Subroutine arguments: 
! --------------------
  INTEGER (KIND=iintegers), INTENT (IN) ::  &
  nsmsteps,             & ! number of small time steps for integration
    ie, je, ke,         & ! dimension indices of input arrays
    kzflat                ! level-index where half-levels bcome flat

  REAL    (KIND=wp   ),     INTENT (IN)    ::  &
    dts,                & ! small time step
    xkd                   ! constant factor for divergence-type damping

  REAL    (KIND=wp   ),     INTENT (INOUT) ::  &
  suten (ie,je,ke),   & ! slow tendency of u-velocity (*dts)
  svten (ie,je,ke),   & ! slow tendency of v-velocity (*dts)
  swten (ie,je,ke+1), & ! slow tendency of w-velocity (*dts)
  stten (ie,je,ke),   & ! slow tendency of temperature (*dts)
  sppten(ie,je,ke),   & ! slow tendency of pressure perturbation (*dts)
  zubdt_west (je,ke), & ! u-boundary tendency at west boundary (total domain)
  zubdt_east (je,ke), & ! u-boundary tendency at east boundary (total domain)
  zvbdt_south(ie,ke), & ! v-boundary tendency at south boundary (total domain)
  zvbdt_north(ie,ke)    ! v-boundary tendency at north boundary (total domain)

! Local scalars:
! -----------------------------
  INTEGER (KIND=iintegers) ::  &
    i,  j,  k,           & !  Loop indices
    kzdims(24),          & !  Vertical dimensions for exchg_datatypes
    nx  ,                & !  timelevel index for start values
    izts,                & !  small timestep loop index
    kp1, km1,            & !  k+1, k-1
    ilowu, jlowv,        & !  start- and end-indices for subdomains with
    iup,   jup             !  or without neighbours

  REAL    (KIND=wp   )     ::  &
    zbp, zbm, zbp2,      & !
    zgbm, zgbp,          & !
    zgdtsh, zrofac,      & !
    zfy, zfyt,           & !
    zpdx, zpdy,          & !
    zdppz, zfact,        & !
    zdpq1, zdpq2, zdpq3, & !
    zppqx, zppqy,        & !
    zpgradx, zpgrady,    & !
    zgu, zgo,            & !
    zaw, zcw,            & !
    zwdiv, zpten,        & !
    zrhoqx, zrhoqy,      & !
    zjpvn, zjpvs,        & !
    zjlur, zjlul,        & !
    zdenom,              & !
    zsim2, zsim1, zsi,   & !
    zsip1, zsip2,        & !
    zvelo_rdx, zdfdx       !

   
! Local (automatic) arrays:
! -----------------------------
  REAL    (KIND=wp   )     ::  &
    za      (ie,je,2:ke),      & !
    zb      (ie,je,2:ke),      & !
    zc      (ie,je,2:ke),      & !
    zd      (ie,je,ke1) ,      & !
    zcp1    (ie,je,ke)  ,      & !
    zcpr    (ie,je,ke)  ,      & !
    zcw1    (ie,je,ke)  ,      & !
    zgfakt  (ie,je,ke)  ,      & !
    zppqz   (ie,je)     ,      & !
    zpps    (ie,ke)     ,      & !
    ztb     (ie,ke)     ,      & !
    zpnj    (ie,je,ke)  ,      & !
    zfx     (je   )     ,      & !
    zfxt    (je   )     ,      & !
    zfydn   (je   )     ,      & !
    zfyd    (je   )     ,      & !
    zfyds   (je   )     ,      & !
    zpi     (ie,je,ke)  ,      & !
    zp0u    (ie,je,ke)  ,      & !
    zp0v    (ie,je,ke)  ,      & !
    zrhoqx_i(ie,je,ke)  ,      & !
    zrhoqy_i(ie,je,ke)  ,      & !
    zwa     (ie,je)              !

  REAL (KIND=wp),     POINTER :: &
    qv  (:,:,:)=> NULL(),      & ! QV at tlev=nnew
    qc  (:,:,:)=> NULL()         ! QC at tlev=nnew

 ! Variables for bottom BC
 !------------------------
   REAL (KIND=wp)      :: &
     zmit1, zmit2,        &
     x1 (ie,je),          &
     xrhsy(ie,je),        &
     xrhsx(ie,je),        &
     xp0y(ie,je),         &
     xp0x(ie,je)

! For error handling
! ------------------
  INTEGER (KIND=iintegers) ::  izerror

  CHARACTER (LEN=255)      ::  yzerrmsg
  CHARACTER (LEN=25)       ::  yzroutine

! Local parameters:
! ----------------
  REAL (KIND=wp),     PARAMETER :: &
  c1 = 1.0_wp/12.0_wp

!
! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Statement Functions
!------------------------------------------------------------------------------

! 3rd order velo*df/dx operator (statement function)
zdfdx(zsim2, zsim1, zsi, zsip1, zsip2, zvelo_rdx) =  &
       c1 * ( zvelo_rdx*(zsim2-8.0_wp*(zsim1-zsip1)      -zsip2) &
        +abs(zvelo_rdx)*(zsim2-4.0_wp*(zsim1+zsip1)+6.0_wp*zsi+zsip2) )

!------------------------------------------------------------------------------
! Begin Subroutine fast_waves    
!------------------------------------------------------------------------------

yzroutine = 'fast_waves'
izerror  = 0_iintegers
yzerrmsg = '   '
kzdims(:)= 0_iintegers

! At the outer boundaries the RK2 differences degenerate
zwa(:,:) = hhl(:,:,ke1)

!------------------------------------------------------------------------------
!  Section 1: Setup of parameters for time integration
!------------------------------------------------------------------------------

! zpi has to be set to 0 (at least at the boundary)
  zpi(:,:,:) = 0.0_wp 

! betasw is the Ikawa beta parameter (=0: centered, =1 backward, =-1 forward)

zbp  = (1.0_wp+betasw)*0.5_wp
zbm  = (1.0_wp-betasw)*0.5_wp
zbp2 = zbp**2
 
IF (my_cart_neigh(1) == -1) THEN
  ilowu = istartu
ELSE
  ilowu = istartu-1
ENDIF

IF (my_cart_neigh(2) == -1) THEN
  jup = jend
ELSE
  jup = jend+1
ENDIF

IF (my_cart_neigh(3) == -1) THEN
  iup = iend
ELSE
  iup = iend+1
ENDIF

IF (my_cart_neigh(4) == -1) THEN
  jlowv = jstartv
ELSE
  jlowv = jstartv-1
ENDIF

!------------------------------------------------------------------------------
!  Section 2: Some preparations for time integration
!------------------------------------------------------------------------------

! Exchange one boundary line of suten and svten
IF (ltime) THEN
  CALL get_timings (i_fast_waves, ntstep, dt, izerror)
  IF (ltime_barrier) THEN
    CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
    CALL get_timings (i_barrier_waiting_dyn, ntstep, dt, izerror)
  ENDIF
ENDIF

kzdims(1:24)=(/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
CALL exchg_boundaries                                                    &
     (15  ,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
      ie, je, kzdims, jstartpar, jendpar, 2, nboundlines,                &
      my_cart_neigh, lperi_x, lperi_y, l2dim,                            &
      10000+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,        &
      suten(:,:,:), svten(:,:,:))

IF (ltime) CALL get_timings (i_fast_waves_comm, ntstep, dt, izerror)

! The start values of the variables for the integration (nx=nold for the
! leapfrog integration and nx=nnow for the 2-TL scheme) are stored on the
! n+1 (nnew) timelevel (including the boundary lines of u and v).
! The corresponding array is then used as working array for timestepping
IF (l2tls) THEN
  nx = nnow
ELSE
  nx = nold
ENDIF

DO  k = 1, ke
  DO    j = jstartu  , jendu  
    DO  i = istartu-1, iendu+1
        u (i,j,k,nnew)  = u (i,j,k,nx)
    ENDDO     
  ENDDO
  DO    j = jstartv-1, jendv+1
    DO  i = istartv  , iendv  
        v (i,j,k,nnew)  = v (i,j,k,nx)
    ENDDO     
  ENDDO
  DO    j = jstart, jend
    DO  i = istart, iend
        w (i,j,k,nnew)  = w (i,j,k,nx)
        t (i,j,k,nnew)  = t (i,j,k,nx)
    ENDDO
  ENDDO       
  DO    j = jlowv, jup
    DO  i = ilowu, iup
        pp(i,j,k,nnew)  = pp(i,j,k,nx)
    ENDDO
  ENDDO       
ENDDO
w(istart:iend,jstart:jend,ke1,nnew ) = w(istart:iend,jstart:jend,ke1,nx)
 
! Factors for the setup of the tridiagonal gaussian matrix and some
! local arrays for horizontal and vertical discretization
 
DO  k = 2, ke
  DO    j = jstart , jend
    DO  i = istart , iend
      zrofac      =  ( rho0(i,j,k)*dp0(i,j,k-1) + rho0(i,j,k-1)*dp0(i,j,k) ) &
                    /( rho (i,j,k)*dp0(i,j,k-1 )+ rho (i,j,k-1)*dp0(i,j,k) )
      zcw1(i,j,k) = 2.0_wp*g*dts*zrofac/(dp0(i,j,k-1)+dp0(i,j,k))
    ENDDO
  ENDDO       
ENDDO 

zfy     = edadlat
zfyt    = edadlat*dts
zgdtsh  = g*dts/2.0_wp
zgbp    = zgdtsh*zbp
zgbm    = zgdtsh*zbm

DO j = jstart - 1, jend + 1
   zfx(j)   = eddlon*acrlat(j,1)
   zfxt(j)  = eddlon*acrlat(j,1)*dts
   zfydn(j) = eddlat*acrlat(j,1)*crlat(j  ,2)
   zfyds(j) = eddlat*acrlat(j,1)*crlat(j-1,2)
   zfyd(j)  = eddlat*acrlat(j,1)*crlat(j  ,1)
ENDDO          

! Get needed tracers out of the tracer structure
CALL trcr_get (izerror, idt_qv, ptr_tlev = nnew, ptr = qv)
IF (izerror /= 0) THEN
  yzerrmsg = trcr_errorstr(izerror)
  CALL model_abort (my_cart_id, izerror, yzerrmsg, yzroutine)
ENDIF
CALL trcr_get (izerror, idt_qc, ptr_tlev = nnew, ptr = qc)
IF (izerror /= 0) THEN
  yzerrmsg = trcr_errorstr(izerror)
  CALL model_abort (my_cart_id, izerror, yzerrmsg, yzroutine)
ENDIF

!------------------------------------------------------------------------------
! RJ: Some more preparations for optimization
!------------------------------------------------------------------------------

! Avoid division by zrhoqx and zrhoqy in the time integration
DO  k = 1, ke
  DO  j = jstart-1, jend+1
    DO  i = istart-1, iend
      zrhoqx          = 0.5_wp*( rho(i,j,k) + rho(i+1,j,k) )
      zrhoqx_i(i,j,k) = zfxt(j)/zrhoqx
    ENDDO
  ENDDO
  DO  j = jstart-1, jend
    DO  i = istart-1, iend+1
      zrhoqy          = 0.5_wp*( rho(i,j,k) + rho(i,j+1,k) )
      zrhoqy_i(i,j,k) = zfyt/zrhoqy
    ENDDO
  ENDDO
ENDDO

DO  k = 1, ke
  kp1 = MIN(k+1, ke)
  km1 = MAX(k-1, 1 )
  DO  j = jstart-1, jend
    DO  i = istart-1, iend
      zdpq1   =   p0(i  ,j  ,kp1) - p0(i  ,j  ,km1)
      zdpq2   =   p0(i+1,j  ,kp1) - p0(i+1,j,  km1)
      zdpq3   =   p0(i  ,j+1,kp1) - p0(i  ,j+1,km1)
      zp0u(i,j,k) =  (p0(i+1,j,k) - p0(i  ,j,k)) / (zdpq2 + zdpq1)
      zp0v(i,j,k) =  (p0(i,j+1,k) - p0(i  ,j,k)) / (zdpq3 + zdpq1)
    ENDDO
  ENDDO
ENDDO


! Setup of some constant arrays
DO    k = 1, ke
  zpnj  (:,:,k) = gamma*dts*( p0(:,:,k) + pp(:,:,k,nnow) )
  zcpr  (:,:,k) = cpdr/rho(:,:,k)
  DO j = jstart-1, jend+1
    DO i = istart-1, iend+1
      zgfakt(i,j,k) = 0.5_wp/(r_d*rho0(i,j,k)*t(i,j,k,nnow))
      zcp1  (i,j,k) = zpnj(i,j,k)*rho0(i,j,k)*g/dp0(i,j,k)
    ENDDO
  ENDDO

  ! Setup of the lower (za), the main (zb) and the upper (zc) diagonals
  ! of the Gauss matrix at the first small timestep
  IF ( k > 1 ) THEN
  DO  j = jstart, jend
    DO  i = istart, iend
      zgu       = 1.0_wp - zgfakt(i,j,k  )*dp0(i,j,k-1)
      zgo       = 1.0_wp + zgfakt(i,j,k-1)*dp0(i,j,k)
      zaw       = zcw1(i,j,k)*zgo*zbp2
      zcw       = zcw1(i,j,k)*zgu*zbp2
      za(i,j,k) =    - zaw*( zcp1(i,j,k-1)-rho0(i,j,k-1)*zgdtsh )
      zb(i,j,k) = 1.0_wp + zcw*( zcp1(i,j,k  )-rho0(i,j,k  )*zgdtsh )    &
                        + zaw*( zcp1(i,j,k-1)+rho0(i,j,k-1)*zgdtsh )
      zc(i,j,k) =       - zcw*( zcp1(i,j,k  )+rho0(i,j,k  )*zgdtsh )
    ENDDO
  ENDDO
  ENDIF
ENDDO

  ! LU-Decomposition of the matrix

  DO  j = jstart, jend
  DO  i = istart, iend
    zdenom    =  zb(i,j,ke)
    zb(i,j,ke) = 1.0_wp/zdenom
    za(i,j,ke) = -za(i,j,ke)*zb(i,j,ke)
  ENDDO
  ENDDO
  DO  k = ke-1, 2, -1
    DO  j = jstart, jend
    DO  i = istart, iend
      zdenom    =  zc(i,j,k)*za(i,j,k+1) + zb(i,j,k)
      zb(i,j,k) = 1.0_wp/zdenom
      za(i,j,k) = -za(i,j,k)*zb(i,j,k)
    ENDDO
    ENDDO
  ENDDO

DO  k = 1, ke
  zpi(istart-1:iend+1,jstart-1:jend+1,k) =                              &
                pp(istart-1:iend+1,jstart-1:jend+1,k,nnew)
ENDDO         

!------------------------------------------------------------------------------
!     ***************************************************************
!     *   begin of time integration using small time steps          *
!     ***************************************************************
!------------------------------------------------------------------------------

loop_over_small_timesteps:  DO  izts = 1, nsmsteps


  !----------------------------------------------------------------------------
  !  Section 3: Integration of the wind components u and v
  !----------------------------------------------------------------------------

  ! Periodic boundary conditions on pressure within the time loop !

  ! exchange pi:
  IF (lperi_x .OR. lperi_y .OR. l2dim) THEN
    IF (ltime) CALL get_timings (i_fast_waves, ntstep, dt, izerror)
    kzdims(1:24)=(/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                  &
         ( 0,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
          ie, je, kzdims, jstartpar, jendpar, 2, nboundlines,              &
          my_cart_neigh, lperi_x, lperi_y, l2dim,                          &
          11000+nexch_tag, .FALSE., ncomm_type, izerror, yzerrmsg,         &
          zpi(:,:,:) )
    IF (ltime) CALL get_timings (i_fast_waves_comm, ntstep, dt, izerror)
  END IF

  loop_over_levels_3:  DO  k = 1, ke

    kp1   = MIN ( ke, k+1 )
    km1   = MAX ( 1 , k-1 )

    IF ( k < kzflat ) THEN

      DO  j = jstartu, jendu
        DO  i = ilowu, iendu   
          zpgradx =  (zpi(i+1,j,k) - zpi(i,j,k)) * zrhoqx_i(i,j,k)
          u(i,j,k,nnew) = u(i,j,k,nnew) - zpgradx + suten(i,j,k)
        ENDDO 
      ENDDO
      DO  j = jlowv, jendv
        DO  i = istartv, iendv   
          zpgrady =  (zpi(i,j+1,k) - zpi(i,j,k)) * zrhoqy_i(i,j,k)
          v(i,j,k,nnew) = v(i,j,k,nnew) - zpgrady + svten(i,j,k)
        ENDDO 
      ENDDO   

    ELSEIF ( (.NOT. ldyn_bbc) .OR.  (ldyn_bbc .AND. k < ke) ) THEN

      zppqz(istart-1:iend+1,jstart-1:jend+1) =                  &
                    zpi(istart-1:iend+1,jstart-1:jend+1,kp1)    &
                   -zpi(istart-1:iend+1,jstart-1:jend+1,km1)

      DO  j = jstartu, jendu
        DO  i = ilowu, iendu   
          zppqx   =  zppqz(i+1,j) + zppqz(i,j)
          zdppz   =  zppqx * zp0u(i,j,k)
          zpdx    =  zpi(i+1,j,k) - zpi(i,j,k)
          zpgradx =  zpdx - zdppz
          u(i,j,k,nnew) = u(i,j,k,nnew) - zpgradx*zrhoqx_i(i,j,k) &
                                        + suten(i,j,k)
        ENDDO 
      ENDDO 
!CDIR OUTERUNROLL=16
      DO  j = jlowv, jendv
        DO  i = istartv, iendv   
          zppqy   =  zppqz(i,j+1) + zppqz(i,j)
          zdppz   =  zppqy * zp0v(i,j,k)
          zpdy    =  zpi(i,j+1,k) - zpi(i,j,k)
          zpgrady =  zpdy - zdppz
          v(i,j,k,nnew) = v(i,j,k,nnew) - zpgrady*zrhoqy_i(i,j,k)   &
                                        + svten(i,j,k)
        ENDDO 
      ENDDO   

    ELSE

      ! Dynamical Bottom boundary condition (formulation of Almut Gassmann)
      !
      !  swten(:,:,ke1) is used to store the w-tendency (hor. adv. and 
      !       buoyancy eff.) at the full level ke 
      !       (computed in src_leapfrog.f90)
      DO  j = jlowv, jendv+1
        DO  i = istartv-1, iendv
          xp0x(i,j) = (p0(i+1,j,ke) - p0(i,j,ke)) * zfx(j)
        ENDDO
      ENDDO
      DO  j = jstartu-1, jendu
        DO  i = ilowu, iendu+1
          xp0y(i,j) = (p0(i,j+1,ke) - p0(i,j,ke)) * zfy
        ENDDO
      ENDDO
      DO j = jstart-1, jend+1
        DO i = istart-1, iend+1
          x1(i,j) = 1.0_wp                                          &
                  - (p0(i,j,ke)+zpi(i,j,ke))*2.0_wp*zgfakt(i,j,ke)  &
                  + swten(i,j,ke1)*rho(i,j,ke) / rho0(i,j,ke)/g/dts
        ENDDO
      ENDDO

      DO  j = jlowv, jendv+1
        DO  i = istartv-1, iendv
          xrhsx(i,j) = - zrhoqx_i(i,j,ke)/dts * (          &
                         ( zpi(i+1,j,ke)-zpi(i,j,ke) )     &
                       + ( x1(i+1,j)+x1(i,j) )*0.5_wp  &
                       * ( p0(i+1,j,ke) - p0(i,j,ke) ) )   &
                       + suten(i,j,ke) / dts
        ENDDO
      ENDDO

!US_CLM      DO  j = jlowv, jendv+1
!US_CLM        DO  i = istartv-1, iendv
!AD_CLM      DO  j = jlowv, jendv
!AD_CLM        DO  i = istartv-1, iendv+1
      DO  j = jstartu-1, jendu
        DO  i = ilowu, iendu+1
          xrhsy(i,j) = - zrhoqy_i(i,j,ke)/dts * (          &
                         ( zpi(i,j+1,ke)-zpi(i,j,ke) )     &
                       + ( x1(i,j+1)+x1(i,j) )*0.5_wp  &
                       * ( p0(i,j+1,ke) - p0(i,j,ke) ) )   &
                       + svten(i,j,ke) / dts
        ENDDO
      ENDDO

      DO  j = jstartu, jendu
        DO  i = ilowu, iendu
          zmit1          = 0.5_wp * ( xp0y(i+1,j)   + xp0y(i,j)          &
                                        + xp0y(i+1,j-1) + xp0y(i,j-1) )      &
                           / ( g*( rho0(i+1,j,ke)+rho0(i,j,ke) ) )
          zmit2          = xp0x(i,j) * 2.0_wp                            &
                           / ( g*( rho0(i+1,j,ke)+rho0(i,j,ke) ) )
          u(i,j,ke,nnew) = u(i,j,ke,nnew)    &
                         + dts / ( 1.0_wp+zmit2**2+zmit1**2 ) * (        &
                           - zmit2 * zmit1 *                                 &
                           0.25_wp * ( xrhsy(i+1,j)   + xrhsy(i,j)       &
                                         + xrhsy(i+1,j-1) + xrhsy(i,j-1) )   &
                           + ( 1.0_wp+zmit1**2 ) * xrhsx(i,j) )
        ENDDO
      ENDDO

      DO  j = jlowv, jendv
        DO  i = istartv, iendv
          zmit1          = 0.5_wp * ( xp0x(i,j+1)   + xp0x(i,j)          &
                                        + xp0x(i-1,j+1) + xp0x(i-1,j) )      &
                           / ( g*( rho0(i,j+1,ke)+rho0(i,j,ke) ) )
          zmit2          = xp0y(i,j) * 2.0_wp                            &
                           / ( g*( rho0(i,j+1,ke)+rho0(i,j,ke) ) )
          v(i,j,ke,nnew) = v(i,j,ke,nnew)    &
                         + dts / ( 1.0_wp+zmit2**2+zmit1**2 ) * (        &
                         - zmit1 * zmit2 *                                   &
                         0.25_wp * ( xrhsx(i,j)   + xrhsx(i-1,j)         &
                                       + xrhsx(i,j+1) + xrhsx(i-1,j+1))      &
                         + ( 1.0_wp+zmit1**2 ) * xrhsy(i,j) )
        ENDDO
      ENDDO

    ENDIF

    ! Set boundary values for u and v from the boundary tendencies
    ! to compute the wind divergence at the scalar gridpoints along
    ! the boundary lines
      IF (.NOT.lperi_x) THEN
      IF (my_cart_neigh(1) == -1) THEN
        DO  j = jstartu, jendu
          u(istartu-1,j,k,nnew) = u(istartu-1,j,k,nnew) + zubdt_west(j,k)
        ENDDO
      ENDIF
      IF (my_cart_neigh(3) == -1) THEN
        DO  j = jstartu, jendu
          u(iendu+1,j,k,nnew) = u(iendu+1,j,k,nnew) + zubdt_east(j,k)
        ENDDO
      ENDIF
      END IF
      IF (.NOT.lperi_y) THEN
      IF (my_cart_neigh(4) == -1) THEN
        DO i = istartv, iendv
          v(i,jstartv-1,k,nnew) = v(i,jstartv-1,k,nnew) + zvbdt_south(i,k)
        ENDDO
      ENDIF
      IF (my_cart_neigh(2) == -1) THEN
        DO i = istartv, iendv
          v(i,jendv+1,k,nnew) = v(i,jendv+1,k,nnew) + zvbdt_north(i,k)
        ENDDO
      ENDIF
    ENDIF

  ENDDO loop_over_levels_3

  ! Exchange u and v for periodic boundary conditions
  IF ( lperi_x .OR. lperi_y .OR. l2dim) THEN
    IF (ltime) THEN
      CALL get_timings (i_fast_waves, ntstep, dt, izerror)
      IF (ltime_barrier) THEN
        CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
        CALL get_timings (i_fast_waves_barrier, ntstep, dt, izerror)
      ENDIF
    ENDIF
    
    kzdims(1:24)=(/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                     &
         (0      , sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
          ie, je, kzdims, jstartpar, jendpar, 2, nboundlines,                 &
          my_cart_neigh, lperi_x, lperi_y, l2dim,                             &
          12000+nexch_tag, .FALSE., ncomm_type, izerror, yzerrmsg,            &
          u(:,:,:,nnew), v(:,:,:,nnew))
    IF (ltime) CALL get_timings (i_fast_waves_comm, ntstep, dt, izerror)
    
  ENDIF

  !----------------------------------------------------------------------------
  !  Section 4: Integration of vertical velocity, temperature and
  !             pertubation pressure using a vertical implicit scheme
  !----------------------------------------------------------------------------

  DO k = 1, ke
    kp1 = MIN( ke, k+1 )
    km1 = MAX( 1 , k-1 )
    zfact = 1.0_wp
    IF(k==1 .OR. k==ke) zfact=0.5_wp
    IF ( k < kzflat ) THEN
!CDIR OUTERUNROLL=16
      DO j = jstart, jend
      DO i = istart, iend
        zpi(i,j,k) = zfx(j)*( u(i,j,k,nnew) - u(i-1,j,k,nnew) ) &
                    + zfydn(j)*v(i,j,k,nnew) - zfyds(j)*v(i,j-1,k,nnew)

        ! Caution: the same block occurs below in the next loop again
        ! (this has been made for better optimization)
        ! Calculation of the explicit part of the 3-d divergence
        zwdiv = zcp1(i,j,k)*( w(i,j,k+1,nnew) - w(i,j,k,nnew) )*zbm

        ! Store the 3-d divergence on zpi
        zpi(i,j,k) = zwdiv - zpnj(i,j,k) * zpi(i,j,k)
      ENDDO
      ENDDO
    ELSE
!CDIR OUTERUNROLL=16
      DO j = jstart, jend
      DO i = istart, iend
        zpi(i,j,k) =                                                         &
         zfx(j)*( u(i,j,k,nnew) - u(i-1,j,k,nnew)                            &
          - zfact * (u(i-1,j,kp1,nnew) - u(i-1,j,km1,nnew))* zp0u(i-1,j,k)   &
          - zfact * (u(i  ,j,kp1,nnew) - u(i  ,j,km1,nnew))* zp0u(i  ,j,k) ) &
          + zfydn(j)*(  v(i,j  ,k,nnew)                                      &
             - zfact * (v(i,j,kp1,nnew) - v(i,j,km1,nnew)) * zp0v(i,j,k) )   &
          + zfyds(j)*(- v(i,j-1,k,nnew)                                      &
             - zfact * (v(i,j-1,kp1,nnew) - v(i,j-1,km1,nnew)) * zp0v(i,j-1,k) )

          ! Calculation of the explicit part of the 3-d divergence
          zwdiv = zcp1(i,j,k)*( w(i,j,k+1,nnew) - w(i,j,k,nnew) )*zbm

          ! Store the 3-d divergence on zpi
          zpi(i,j,k) = zwdiv - zpnj(i,j,k) * zpi(i,j,k)
      ENDDO
      ENDDO
    ENDIF
  ENDDO      

  loop_from_south_to_north:  DO  j = jstart, jend
  !----------------------------------------------
 
    DO   k = 1, ke
      DO  i = istart, iend
        ! Explicit part of total pressure tendency
        zpten      = zpi(i,j,k) + sppten(i,j,k)                &
                   + rho0(i,j,k)*( w(i,j,k+1,nnew) + w(i,j,k,nnew) )*zgbm

        ! Explicit part of pp(n+1)
        zpps (i,k)        = pp(i,j,k,nnew) + zpten*zbp 
        pp   (i,j,k,nnew) = pp(i,j,k,nnew) + zpten
      ENDDO   
      DO  i = istart, iend
        ztb   (i,k) = ( t(i,j,k,nnew)*r_d*rho0(i,j,k)-p0(i,j,k) )*zgfakt(i,j,k)
      ENDDO   

      ! Calculate the righthand side of the matrix eqation
      IF ( k > 1 ) THEN
        DO  i = istart, iend
          zgu       = 1.0_wp - zgfakt(i,j,k  )*dp0(i,j,k-1)
          zgo       = 1.0_wp + zgfakt(i,j,k-1)*dp0(i,j,k  )
          zd(i,j,k) = w(i,j,k,nnew) + swten(i,j,k) &
                         + zcw1(i,j,k)*( zgu*zpps(i,k) - zgo*zpps(i,k-1) &
                         + ztb(i,k)*dp0(i,j,k-1) + ztb(i,k-1)*dp0(i,j,k) )
        ENDDO
      ENDIF
    ENDDO
  ENDDO loop_from_south_to_north

  ! Set the upper (w=0) and lower (w=u*dz/dx) boundary conditions on w
  IF (llm) THEN
    DO  j = jstart, jend
      DO  i = istart, iend
        w (i,j,  1,nnew) = 0.0_wp
        w (i,j,ke1,nnew) = 0.0_wp
        zd(i,j,ke1)      = w(i,j,ke1,nnew)
      ENDDO
    ENDDO
  ELSE
    IF( .NOT. l2tls ) THEN
      DO  j = jstart, jend
        DO  i = istart, iend
          zjpvn = ( hhl(i  ,j+1,ke1) - hhl(i  ,j  ,ke1) )*v(i,j  ,ke,nnew)
          zjpvs = ( hhl(i  ,j  ,ke1) - hhl(i  ,j-1,ke1) )*v(i,j-1,ke,nnew)
          zjlur = ( hhl(i+1,j  ,ke1) - hhl(i  ,j  ,ke1) )*u(i  ,j,ke,nnew)
          zjlul = ( hhl(i  ,j  ,ke1) - hhl(i-1,j  ,ke1) )*u(i-1,j,ke,nnew)
          w(i,j,ke1,nnew) =   0.5_wp*( zfydn(j)*zjpvn + zfyds(j)*zjpvs  &
                                  + zfx(j)*( zjlur + zjlul ) )
          zd(i,j,ke1)     = w(i,j,ke1,nnew)
          w(i,j,  1,nnew) = 0.0_wp
        ENDDO
      ENDDO
    ELSE
      DO  j = jstart, jend
        DO  i = istart, iend
          zjpvn = 0.5_wp*(v(i,j,ke,nnew)+v(i,j-1,ke,nnew))*zfyd(j)
          zjlur = 0.5_wp*(u(i,j,ke,nnew)+u(i-1,j,ke,nnew))*zfx(j)
          zwa(i,j) = (zdfdx(hhl(i-2,j,ke1),hhl(i-1,j,ke1),hhl(i,j,ke1), &
                      hhl(i+1,j,ke1),hhl(i+2,j,ke1),zjlur) + &
                      zdfdx(hhl(i,j-2,ke1),hhl(i,j-1,ke1),hhl(i,j,ke1), &
                      hhl(i,j+1,ke1),hhl(i,j+2,ke1),zjpvn))    &
                      *dt*0.5_wp + hhl(i,j,ke1)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  DO  j = jstart, jend

    ! Up-down algorithm to solve the tridiagonal matrix equation for w

!CDIR OUTERUNROLL=8
    DO  k = ke, 2, -1
      DO  i = istart, iend
        zd(i,j,k) = ( zd(i,j,k) - zd(i,j,k+1)*zc(i,j,k) )*zb(i,j,k)
      ENDDO
    ENDDO

!CDIR OUTERUNROLL=8
    DO  k = 1, ke-1
      DO  i = istart, iend
        w(i,j,k+1,nnew) = za(i,j,k+1)*w(i,j,k,nnew) + zd(i,j,k+1)
      ENDDO
    ENDDO
 
  ENDDO

  ! Calculation of pp(n+1) and t(n+1)
  ! Completition of the divergence damping term zpi

  DO k = 1, ke
    DO  j = jstart, jend
      DO  i = istart, iend
        zwdiv           =  zcp1(i,j,k)*(w(i,j,k+1,nnew)-w(i,j,k,nnew))*zbp  
        pp (i,j,k,nnew) =  pp(i,j,k,nnew) + zwdiv  &
                         + rho0(i,j,k)*(w(i,j,k+1,nnew)+w(i,j,k,nnew))*zgbp
        zpi(i,j,k     ) =  zpi(i,j,k) + zwdiv
      ENDDO
    ENDDO      
    DO  j = jstart, jend
      DO  i = istart, iend
        t  (i,j,k,nnew) =  t(i,j,k,nnew) + stten(i,j,k)                     & 
                         + zpi(i,j,k)*zcpr(i,j,k)
      ENDDO
    ENDDO      
  ENDDO

  IF (izts < nsmsteps) THEN
    ! Addition of 3-d divergence to perturbation pressure for divergence damping
    DO  k = 1, ke
     zpi(istart:iend,jstart:jend,k) =                                 &
                       pp(istart:iend,jstart:jend,k,nnew)             &
                         + xkd*zpi(istart:iend,jstart:jend,k)
    ENDDO         
  
    ! Exchange one boundary line of zpi
    IF (ltime) THEN
      CALL get_timings (i_fast_waves, ntstep, dt, izerror)
      IF (ltime_barrier) THEN
        CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
        CALL get_timings (i_fast_waves_barrier, ntstep, dt, izerror)
      ENDIF
    ENDIF
    kzdims(1:24)=(/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                 &
         (19,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
          ie, je, kzdims, jstartpar, jendpar, 1, nboundlines,             &
          my_cart_neigh, lperi_x, lperi_y, l2dim,                         &
          15000+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,     &
          zpi(:,:,:) )
    IF (ltime) CALL get_timings (i_fast_waves_comm, ntstep, dt, izerror)

  ENDIF

  !----------------------------------------------------------------------------
  !     ***************************************************************
  !     *   end of time integration using small time steps            *
  !     ***************************************************************
  !----------------------------------------------------------------------------

ENDDO loop_over_small_timesteps

!------------------------------------------------------------------------------
!  Section 5a: Set free-slip lateral boundary conditions on w for
!              non-periodic boundary conditions
!              (applies also for single-processor and/or 2D case)
!------------------------------------------------------------------------------

IF (lw_freeslip .AND. .NOT.lperi_x) THEN
  DO  k = 1, ke1
    
    IF (my_cart_neigh(1) == -1) THEN
      DO  j = jstart, jend
        DO i = 1, nboundlines
          w(i,j,k,nnew) = w(istart,j,k,nnew)
        ENDDO
      ENDDO
    ENDIF
    
    IF (my_cart_neigh(3) == -1) THEN
      DO  j = jstart, jend
        DO i = ie-nboundlines+1, ie
          w(i,j,k,nnew) = w(iend  ,j,k,nnew)
        ENDDO
      ENDDO
    ENDIF
    
  END DO
  
ENDIF

IF (lw_freeslip .AND. .NOT.lperi_y) THEN

  DO  k = 1, ke1

    IF (my_cart_neigh(4) == -1) THEN
      DO  j = 1, nboundlines
        w(:,j,k,nnew) = w(:,jstart,k,nnew)
      ENDDO
    ENDIF
    
    IF (my_cart_neigh(2) == -1) THEN
      DO  j = je-nboundlines+1, je
        w(:,j,k,nnew) = w(:,jend  ,k,nnew)
      ENDDO
    ENDIF
    
  ENDDO
ENDIF

!------------------------------------------------------------------------------
!  Section 5b: Set periodic lateral boundary conditions on w and all
!              all other variables if required
!------------------------------------------------------------------------------

! Boundary exchange:
IF ( lperi_x .OR. lperi_y .OR. l2dim) THEN
  IF (ltime) THEN
    CALL get_timings (i_fast_waves, ntstep, dt, izerror)
    IF (ltime_barrier) THEN
      CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
      CALL get_timings (i_fast_waves_barrier, ntstep, dt, izerror)
    ENDIF
  ENDIF
  
! fuo: is the exchange of qv and qc at nnew really needed here?
  kzdims(1:24)=(/ke1,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
  CALL exchg_boundaries                                                       &
       ( 0      ,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
        ie, je, kzdims, jstartpar, jendpar, 2, nboundlines,                   &
        my_cart_neigh, lperi_x, lperi_y, l2dim,                               &
        20000+nexch_tag, .FALSE.   , ncomm_type, izerror, yzerrmsg,           &
        w (:,:,:,nnew), t (:,:,:,nnew), qv(:,:,:), qc(:,:,:),                 &
        pp(:,:,:,nnew) )
  IF (ltime) CALL get_timings (i_fast_waves_comm, ntstep, dt, izerror)
  
ENDIF

IF (ltime) CALL get_timings (i_fast_waves, ntstep, dt, izerror)

!------------------------------------------------------------------------------
!  End of the external procedure fast_waves
!------------------------------------------------------------------------------

END SUBROUTINE fast_waves
