!+ Source module for computing transfer coefficients
!------------------------------------------------------------------------------

 MODULE turbulence_tran

!------------------------------------------------------------------------------
!
! Description:
!
!  The module turbulence_tran calculates the tendencies for turbulent
!  vertical transport of momentum and heat and the coefficients
!  for turbulent transfer as well.
!
! The turbulence model (with some Prandtl-layer approximations is used 
! for the calculation of turbulent transfer between atmosphere and the
! lower boundary too.
!
! Current Code Owner: DWD, Matthias Raschendorfer
!  phone:  +49  69  8062 2708
!  fax:    +49  69  8062 3721
!  email:  matthias.raschendorfer@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.30       1999/06/24 Matthias Raschendorfer
!  Initial release
! 1.33       1999/10/14 Matthias Raschendorfer
!  Improved calculation of the vertical gradients at the z0-level, using a
!  simlified integration of the inverse of the turbulent diff. coeff. from
!  this level to the lowest model level.
!  Introduction of an optional correction of the vertical gradients controled
!  by the LOOGECAL parameter lprfcor.
!  Calculation of an seperate roughness parameter z0h for scalars as a
!  function of z0 and the effective surface area index sai.
!  Reformulation of the laminar resistance.
!  Reformulation of the Charnock-Formula using the additional parameter alpha1.
!  Using the redefined parameter c_g (including the factor 3).
!  Inclusion of the diagnosis of 2m- and 10m- values (former routine
!  'synop_diag').  Attention: this quite fundamental reformulation of the
!  surface scheme has not been implementet in turbdiff.incf yet.
! 1.34       1999/12/10 Matthias Raschendorfer
!  Reformulation of the consideration of a roughness canopy in the part 'synop_diag'.
! 1.37       2000/03/24 Matthias Raschendorfer
! Introduction of an canopy resistance for momentum 
! 1.39       2000/05/03 Ulrich Schaettler
!  Global variable names lam_m and lam_h changed to rlam_m and rlam_h.
! 2.2        2000/08/18 Matthias Raschendorfer
!  Introduction of the molecular diffusion coefficients as minimal values 
!  for the turbulent ones. Calculation of the SQRT of the wind energy in 
!  the modified Charnock formula at the lowest model level
!  using the TKE at that level and not (as before) at the lower boundary.
!  Now sai is defined by (sai_before + 1).
! 2.3        2000/11/15 Guenther Doms
!  Some local variable have been redifined to allow for reproducible
!  results on MPP platforms.
! 2.12       2001/11/07 Matthias Raschendorfer
!  Limitation of 'gama' in the calulation of the Prandtl-layer resistance
!  in order to avoid unrealistic deformation of the Prandtl-layer profiles
!  which were the reason of some "jumps" in the 2m-temperature during periods
!  of stabilisation at the surface.
!  Introduction of a minimal roughness length for the initialization of 'gz0'
!  over water surfaces.
! 2.15       2002/03/20 Matthias Raschendorfer
!  Modified interpolation of the 10m wind vector with the help of the roughness length
!  of a typical SYNOP-station (z0m_d).
!  Multiplication of z0m with the laminar correction factor. 
! 2.16       2002/03/28 Matthias Raschendorfer
!  Modification of the laminar limit without touching the value of z0m, as the actual
!  formulation became numerical unstable during a test assimilation run.
! 2.17       2002/05/08 Ulrich Schaettler
!  Optimizations for vectorization: splitted the big loop in Section 4
! 2.18       2002/07/16 Ulrich Schaettler
!  Added a test for variables dz_s0_h, fac_h_2d before a division
!  (because of problems on the NEC)
! 2.19       2002/10/24 Ulrich Schaettler
!  Deleted 2 lines of code that wrongly remained after update 2.17
!  Adaptations to 2-timelevel scheme (use of ntlev)
! 3.6        2003/12/11 Ulrich Schaettler
!  Optimizations for vectorization: modification of 2 DO WHILE loops
!  Modification of an IF...ELSEIF... Statement for vectorization in Section 4f
! 3.7        2004/02/18 Matthias Raschendorfer
!  Introduction of the parameter rat_sea
! 3.14       2005/01/25 Jochen Foerstner
!  Introduced new types of turbulent diffusion parameterizations
!  Replaced SIGN function
!  Adjusted upper boundary of DO-WHILE loops (by NEC)
! 3.16       2005/07/22 Matthias Raschendorfer
!  Some adaptations to the modifications in 'turbdiff.incf'
! 3.18       2006/03/03 Matthias Raschendorfer
!  Introduction of rh_2m and limitation of td_2m
! 3.19       2006/04/25 Jochen Foerstner / Matthias Raschendorfer
!  Application of the lower limit for the roughness length over sea
!  not only during the initialization of 'gz0'
! 3.21       2006/12/04 Dmitrii Mironov
!  Changes to use the FLake model, Ulrich Schaettler
!  Changed interface to subroutine cloud_diag from meteo_utilities
! V3_23        2007/03/30 Matthias Raschendorfer
!  Renaming some variables for better understanding.
!  Introducing the reduction factor for evaporation tfv.
!  Introduction of some output variables for SCLM.
!  Change of some parameter names.
!  Eliminated stab_funct.incf and turb_param.incf
!  (Substitution of file 'stab_funct.incf' by a SUBROUTINE)
!  New file 'statement_functs.incf' containing statement functions
!  Removing the final tfh- and tfm- calculation with respect to the lam. layer,
!  which was the wrong definition for its use in 'turbdiff.icnf'
! V3_25        2007/05/21 Ulrich Schaettler
!  Moved an IF-clause outside DO-Loops in line 1193ff
! V4_3         2008/02/25 Matthias Raschendorfer
!  Changing interpolation onto diagnostic levels,
!  in particular without an exponential canopy profile,
!  but with a dianostic Prandtl layer interpolation even for scalars,
!  using an adopted canopy layer resistance.
!  Calculation of the 3D diagnostic field 'edr' for the eddy dissipotion rate.
!  Changing the treatment of parameter field dd(:,:,:) to achiev better vectorisation
!  in SUBROUTINE 'stab_funct'.
! V4_8         2009/02/16 Ulrich Schaettler
!  Introduced itype_diag_t2m and "old" 2m temperature as an option
! V4_10        2009/09/11 Matthias Raschendorfer, Jan-Peter Schulz
!  Removing the horizontal loops for SCLM treatment. 
!  Modifications for seaice model: eliminated l_ls_ice, introduced lseaice
!   (Jan-Peter Schulz)
! V4_12        2010/05/11 Ulrich Schaettler
!  Renamed t0 to t0_melt because of conflicting names
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_20        2011/08/31 Ulrich Schaettler
!  Restructured SR turbdiff now as an extra module with argument list
!  for future unified COSMO-ICON physics
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Allocation of tkvm, tkvh for 1:ke1 in the vertical
! V5_1         2014-11-28 Matthias Raschendorfer, Oliver Fuhrer
!  Introduction of LOGICAL lscm for SC-runs and employing SC-data.
!  Replaced ireals by wp (working precision) (OF)
! V5_4a        2016-05-10 Matthias Raschendorfer
!  Use new module turb_data instead of data_turbulence
!  Moved several variables from data_runcontrol to turb_data
!
! Code Description:
! Language: Fortran 90 
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================


! Modules used:
!------------------------------------------------------------------------------

USE data_constants,  ONLY:   &
    t0_melt, grav=>g, con_m, con_h, b1, b2w, b3, b4w, cp_d, rvd_m_o, r_d,  &
    lhocp, rdv, lh_v, b234w, o_m_rdv, rdocp, uc1, uc2, ucl,                &
    b2i, b4i

USE data_flake,      ONLY:   &
    h_ice_min_flk

USE data_parameters, ONLY:   &
    wp,           & ! KIND-type parameter for real variables
    iintegers       ! KIND-type parameter for standard integer variables

USE data_runcontrol, ONLY:   &
    lseaice, llake, itype_synd, itype_diag_t2m, lscm

USE turb_data,       ONLY:   &
    lprfcor, icldm_tran, imode_tran, itype_wcld,                 &
    akt, alpha0, rlam_mom, rlam_heat, rat_lam, rat_sea, rat_can, &
    z0m_dia, z0_ice, zt_ice, len_min, vel_min, d_m=>d_mom,       &
    clc_diag, q_crit, epsi, tkesmot, it_end

USE meteo_utilities,       ONLY:   cloud_diag
USE turbulence_utilities,  ONLY:   stab_funct, turb_param,    &
                                   turb_cloud, diag_level,    &
                                   tet_g, rim, l_scal, c_g,   &
                                   a_3, a_5, a_6, b_1, b_2, d_4

!SCLM---------------------------------------------------------------------------
#ifdef SCLM
USE data_1d_global, ONLY : &

    i_cal, im, jm, &

    UUA, VVA, UWA, VWA, WWA, UST, TTA, TWA, SHF, LHF, &
    TKE_SCLM=>TKE, BOYPR, SHRPR, DISSI, TRANP
#endif
!SCLM---------------------------------------------------------------------------

!==============================================================================

IMPLICIT NONE

!==============================================================================

REAL (KIND=wp),     PARAMETER :: &
    z0=0.0_wp,&
    z1=1.0_wp,&
    z2=2.0_wp,&
    z3=3.0_wp,&
    z4=4.0_wp,&
    z5=5.0_wp,&
    z6=6.0_wp,&
    z7=7.0_wp,&
    z8=8.0_wp,&
    z9=9.0_wp,&
    z10=10.0_wp

REAL (KIND=wp)     :: &
    z1d2=z1/z2,&
    z1d3=z1/z3,&
    z2d3=z2/z3,&
    z3d2=z3/z2

#ifdef SCLM
REAL (KIND=wp)     :: &
    teta, qvap, rhos
#endif

INTEGER (KIND=iintegers) :: &
    istat=0

LOGICAL :: lerror=.FALSE.

!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure turbtran in "src_turbdiff" for computing the coefficients
!+ for turbulent transfer
!------------------------------------------------------------------------------

SUBROUTINE turbtran (dt_tke, lstfnct, iini, nprv, ntur, ntim, nvor,           &
                     ie, je, ke, ke1, kcm, vst,                               &
                     istartpar, iendpar, jstartpar, jendpar,                  &
                     hhl, fr_land, depth_lk, sai,                             &
                     u, v, t, qv, qc, prs, ps, qv_s, t_g, h_ice,              &
                     gz0, tcm, tch, tfm, tfh, tfv,                            &
                     tke, tkvm, tkvh, rcld,                                   &
                     t_2m, qv_2m, td_2m, rh_2m, u_10m, v_10m,                 &
                     edr, shfl_s, lhfl_s,                                     &
                     yerrormsg, ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!
!     Es werden die Transferkoeffizienten fuer den Austausch von Impuls,
!     sowie fuehlbarer und latenter Waerme bestimmt und die Modellwerte
!     fuer die bodennahen Messwerte (in 2m und 10m) berechnet.
!     
! Method:
!
!     Hierzu wird der gesamte Bereich von den festen Oberflachen am
!     Unterrand des Modells bis hin zur untersten Hauptflaeche in
!     die drei Teilbereiche:
!
!     - laminare Grenzschicht (L-Schicht)
!     - turbulente Bestandesschicht (B-Schicht)
!     - turbulnte Prandtl-Schicht (P-Schicht)
!
!     aufgeteilt. Fuer jeden dieser Teilbereiche wird (getrennt nach
!     skalaren Eigenschaften und Impuls) ein zugehoeriger Transport-
!     widerstand berechnet, der gleich einer effektiven Widerstands-
!     laenge ( dz_(sg, g0, 0a)_(h,m) ) dividiert durch den Diffusions-
!     koeffizienten am Unterrand der P-Schicht (Niveau '0') ist.
!     Die Konzentrationen am Unterrand der B-Schicht, also im
!     Abstand der L-Schicht-Dicke entlang der festen Oberflaechen,
!     haben den Index 'g' (ground) und die Oberflaechenkonzentrationen
!     den Index 's'. Groessen fuer den Imoulst haben den Index 'm' (momentum)
!     und solche fuer skalare Eigenschaften 'h' (heat).
!     Der Widerstand der P-Schicht vom Nieveau '0' bis zum Niveau 'a' 
!     (atmospheric) der untersten Hauptflaeche wird durch vertikale 
!     Integration der Modellgleichungen in P-Schicht-Approximation 
!     (vertikal konstante Flussdichten, turbulente Laengenskala lin. Funkt.
!      von der Hoehe) gewonnen.
!     Dabei wird das atmosphaerische Turbulenzschema aus der Subroutine 
!     'turbdiff' benutzt, so dass alo keine empirischen Profilfunktionen 
!     benutzt werden. Zur Vereinfachung der Integration  wird das Produkt 
!     aus turbulenter Geschwindigkeitsskala 'q' und der Stabilitaetsfunktion
!     s(h,m), alo die stabilitaetsabhaengige turb. Geschwindigkeitsskala 
!     'v' innerhalb der P-Schicht als linear angesehen.
!     Die turb. Laengenskala im Niveau '0' wird mit der Rauhigkeitslaenge 'z0'
!     (multipliziert mit der v.Kaman-Konstanten) gleichgesetzt. Formal werden
!     dann fuer das Nieveau '0' Vertikalgradienten und auch Diffusions-
!     koeffizienten abgeleitet.
!     Unter der Annahme, dass 'v' innerhalb der B-Schicht konstant bleibt,
!     ergibt sich die laminare Widerstandslaenge dz_sg als prop. zu 'z0' 
!     und die Widerstandsstrecke durch die B-Schicht als prop. zu 
!     'z0*ln(delta/z0)', wobei 'delta' die Dicke der L-Schicht ist, die der 
!     Abstand von einer ebenen Wand sein soll in dem der turbulente 
!     Diffusionskoeffizient fuer Impuls gleich dem molekularen ist.
!     Ferner wird angenommen, dass die Widerstaende durch die L- und 
!     B-Schicht prop. zur effektiven Quellflaech der Bestandeselemente 
!     zuzueglich der Grundflaeche des Erdbodens sind. Die Bestandesoberflaechen
!     werden durch den Wert 'sai' (surface area index) ausgedrueckt und setzt
!     sich aus dem Flaechenindex der transpirierenden Oberflaechen 'lai' 
!     (leaf area index) und dem fuer die nicht transpirierenden Flaechen 
!     zusammen. Im Falle nicht benetzter Oberlfaechen hat die latente Waerme 
!     i.a. eine kleinere Quellflaeche als die fuehlbare Waerme, so dass die 
!     Wiederstaende fuer beide Groessen unterschieden werden muessten.
!     Um dies zu vermeiden, wird nur der Widerstand fuer die fuehlbare Waerme
!     berechnet. Dafuer wird aber bei der Berechnung der effektiven 
!     Oberflaechenkonzentration 'qv_s' der spez. Feuchtigkeit in Subroutine 
!     'terra1' dieser Effekt beruecksichtigt.
!     Beim vertikalen Impulstransport ist aber noch die zusaetzliche 
!     Impulssenke innerhalb der B-Schicht durch die Wirkung der Formreibungs-
!     kraft zu beruecksichtigen, was durch einen zusaetzlichen Flaechenindex 
!     'dai' (drag area index) bewerkstelligt wird.
!
!     Die Vertikalprofile aller Eigenschaften innerhalb der P-Schicht ergeben
!     sich aus dem vertikal integrierten Turbulenzmodell in P-Schicht-
!     Approximation zu logarithmischen Funktionen, welche durch die 
!     thermische Schichtung modifiziert sind. Wie bereits erwaehnt, ist die 
!     Stabilitaetsfunktion nur noch von Konstanten des atmosphaerischen
!     Turbulenzmodells abhaengig. Das Transferschema ist somit auch automatisch
!     konsistent zum oben anschliessenden Turbulenzmodell formuliert.
!
!     Die Profilfunktionen innerhalb der B-Schicht ergeben sich aus der 
!     Annahme eines Gleichgewichtes zwischen vertikalen Flussdichtedivergenzen
!     und Quellstaerken durch die laminaren Grenzschichten der Rauhigkeits-
!     elemente bei vertikal konstanten Bestandeseigenschaften zu exponentiellen
!     Funktionen. Durch die Bedingung eines glatten Ueberganges zwischen beiden
!     Profiltypen im Niveau '0' und der Bedingung, dass im Abstand einer 
!     effektiven Bestandesdicke 'Hb' unterhalb des Nieveaus '0' die Bestandes-
!     profile in die Konzentration am Unterrand der B-Schicht (Niveau mit 
!     Index 'g') uebergehen, ist das gesamte Transferschema geschlossen und 
!     es kann auch der "drag area index" 'dai', sowie die Bestandeshoehe 
!     'Hb' selbst eliminiert werden.
!
!     Zur Charakterisierung des Oberflaechentransfers werden dann nur die 
!     externen Parameter 'z0', 'sai', 'lai' und je ein globaler Parameter 
!     fuer den laminaren Grenzschichtwiderstand des skalaren - und des 
!     Impulstransportes benoetigt '(lam_(h,m)'. Hieraus koennte auch eine 
!     aequivalente Rauhigkeitslaenge fuer Skalare 'z0h' berechnet werden. 
!     Die Oberfalaechenkonzentrationen (Niveau mit Index 's') fuer die skalaren
!     Groessen werden im Modul 'terra' berechnet. Fuer den Impuls gilt die 
!     Haftbedingung. Im Grundniveau des atmosphaerischen Modells 
!     (Niveau '0') verschwindet also der Wind i.a. nicht; dies ist erst 
!     entlang der festen Oberfalechen der Fall. Die bodennahen synoptischen 
!     Niveaus werden nun vom Niveau 'z=-Hb', also von der effektiven Bestandes-
!     grundflaeche (Umsatzniveau) aus gezaehlt. Ist z.B. 'Hb>2m', werden 
!     die 2m-Werte entlang der exponentiellen Bestandesprofile ausgeweret. 
!     Ist 'Hb<2m', wird das logarithmische Profil in der Hoehe '2m-Hb' entlang
!     dinnerhalb der P-Schicht ausgewertet.
!     Die resultierenden Transferkoeffizienten 'tc(h,m)' sind die Kehrwerte
!     des Gesamtwiderstandes von den festen Oberflaechen (Neviau 's') bis 
!     zur untersten Modellhauptflaeche (Niveau 'a').
!     Die turbulenten Diffusionskoeffizienten 'tkv(h,m)' fuer den vertikalen
!     Index 'ke1', beziehen sich aber auf den Unterrand des atmosphaerischen 
!     Modells (Niveau '0').
!     Mit Hilfe der Felder 'tf(mh)' werden noch Reduktionsfaktoren der 
!     Transferkoeffizienten durch die Wirkung der L-Schicht uebertragen. 
!     Diese koennen im Modul 'terra' benutzt werden, um ev. das effektive 
!     'qv_s' so zu bestimmen, als gaebe es fuer fuehlbare und latente Waerme
!     unterschiedliche Parameter fuer den laminaren Transportwiderstand.
!     Zu beachten ist, dass im Falle eines vertikal vom atmosphaerischen Modell
!     aufgeloesten 'Makrobestandes' (z.B. Bebauung, Wald oder subskalige 
!     Orographie) das Transferschema genauso wie im Falle eines nicht 
!     aufgeloesten Bestandes angewendet wird. Allerdings beziehen sich die 
!     den Bestand des Transferschemas charakterisierenden externen Parameter
!     dann auf den nicht vertikal aufgeloesten verbleibenden 'Mikrobestand',
!     der ev. allein durch niedrigen Bewuchs gebildet wird.
!     Im Transferschema eingearbeitet ist auch dei iterative Bestimmmung der 
!     Rauhigkeitslaenge der Meeresoberflaeche gemaess einer modifizierten 
!     Charnock-Formel, bei der die Wellenerzeugung bei verschwindenden 
!     mittleren Wind mit hilfe der zur TKE ausgedrueckt wird.
!
!==============================================================================

! Declarations:

!------------------------------------------------------------------------------

! Arguments with intent(in):

  REAL (KIND=wp),           INTENT (IN) :: &
    dt_tke     ! time step for TKE stepping

  INTEGER (KIND=iintegers), INTENT(IN) :: &
    iini,    & ! type of initialization (0: no, 1: separate before the time loop
               !                              , 2: within the first time step)
    nprv,    & ! previous    time step index of tke
    ntur,    & ! current new time step index of tke
    ntim       ! number of tke time levels

  INTEGER (KIND=iintegers), INTENT(INOUT) :: &
    nvor       ! actual time step index of TKE (also used in turbdiff)

  INTEGER (KIND=iintegers), INTENT (IN) :: &
    ie, je, ke, ke1, kcm,                  & ! dimensions of the array
    istartpar, iendpar, jstartpar, jendpar   ! start- and end-indices for computation

  INTEGER (KIND=iintegers), INTENT (IN) :: &
    vst             ! velocity component staggering

  LOGICAL,                  INTENT (IN) :: &
    lstfnct

!------------------------------------------------------------------------------

REAL (KIND=wp),             INTENT(IN) :: &
    hhl     (ie,je,ke1), & ! height of model half levels                   ( m )
    fr_land (ie,je),     & ! land portion of a grid point area             ( 1 )
    depth_lk(ie,je),     & ! lake depth                                    ( m )
    sai     (ie,je)        ! surface area index                            ( 1 )

!------------------------------------------------------------------------------

  REAL (KIND=wp),           INTENT(IN)    :: &

    ! Fields for surface values and soil/canopy model variables:
    ! ----------------------------------------------------------

    h_ice(ie,je),        & ! ice thickness                                 (  m  )
    ps   (ie,je),        & ! surface pressure                              ( pa  )
    qv_s (ie,je),        & ! specific water vapor content on the surface   (kg/kg)
    t_g  (ie,je)           ! specific water vapor content on the surface   (kg/kg)

!------------------------------------------------------------------------------

  REAL (KIND=wp),           INTENT(INOUT) :: &

    ! Atmospheric model variables:
    ! ----------------------------

!US the latest news from ICON: not only the microphysics, but also the turbulence
!   scheme (the fast physics) should update the prognostic variables, not the
!   tendencies

    u  (ie,je,ke),       & ! zonal wind speed                              ( m/s )
    v  (ie,je,ke),       & ! meridional wind speed                         ( m/s )
    t  (ie,je,ke),       & ! temperature                                   (  k  )
    qv (ie,je,ke),       & ! specific water vapor content                  (kg/kg)
    qc (ie,je,ke),       & ! specific cloud water content                  (kg/kg)
    prs(ie,je,ke)          ! full pressure                                 ( pa  )

!------------------------------------------------------------------------------

  REAL (KIND=wp),           INTENT(INOUT) :: &

    ! Atmospheric variables of the turbulence model:
    ! ----------------------------------------------

    tke(ie,je,ke1,ntim), & ! SQRT(2*TKE); TKE='turbul. kin. energy'        ( m/s )
                           ! (defined on half levels)
    tkvm(ie,je,ke1),     & ! turbulent diffusion coefficients for momentum (m/s2 )
    tkvh(ie,je,ke1),     & ! turbulent diffusion coefficients for heat     (m/s2 )
                           ! and moisture
    rcld(ie,je,ke1)        ! standard deviation of the saturation deficit
                           ! (as input and output)
                           ! fractional cloud cover (in turbdiff)            --

!------------------------------------------------------------------------------

  REAL (KIND=wp),           INTENT(INOUT), OPTIONAL :: &

    ! Tendency fields for the prognostic variables:
    ! ---------------------------------------------

     edr(ie,je,ke1)        ! eddy dissipation rate of TKE (EDR)            (m2/s3)


!------------------------------------------------------------------------------

  REAL (KIND=wp),           INTENT(INOUT) :: &

    ! Diagnostic surface variable of the turbulence model:
    ! ----------------------------------------------------

    gz0(ie,je),          & ! roughness length * g of the vertically not
                           ! resolved canopy                               (m2/s2)
    tcm(ie,je),          & ! turbulent transfer coefficients for momentum    --
    tch(ie,je),          & ! turbulent transfer coefficients for heat        --
    tfm(ie,je),          & ! factor of laminar transfer of momentum          --
    tfh(ie,je),          & ! factor of laminar transfer of scalars           --
    tfv(ie,je)             ! laminar reduction factor for evaporation        --

!------------------------------------------------------------------------------

  REAL (KIND=wp),           INTENT(OUT)   :: &

    ! Diagnostic near surface variables:
    ! ----------------------------------

    t_2m (ie,je),        & ! temperature in 2m                             (  K  )
    qv_2m(ie,je),        & ! specific water vapor content in 2m            (kg/kg)
    td_2m(ie,je),        & ! dew-point in 2m                               (  K  )
    rh_2m(ie,je),        & ! relative humidity in 2m                       (  %  )
    u_10m(ie,je),        & ! zonal wind in 10m                             ( m/s )
    v_10m(ie,je)           ! meridional wind in 10m                        ( m/s )


!------------------------------------------------------------------------------

  REAL (KIND=wp),           INTENT(OUT), OPTIONAL   :: &

     shfl_s(ie,je),      & ! sensible heat flux at the surface             (W/m2) (positive upward)
     lhfl_s(ie,je)         ! latent   heat flux at the surface             (W/m2) (positive upward)

!------------------------------------------------------------------------------

  INTEGER (KIND=iintegers), INTENT(OUT)   :: &
    ierrstat     ! error status  

  CHARACTER (LEN= *),       INTENT(OUT)   :: &
    yerrormsg    ! error message 

!------------------------------------------------------------------------------
! US: from now on old turbtran

! Local parameters:

      INTEGER (KIND=iintegers), PARAMETER :: &

!             Indexgrenzen:

              nred=4,      & !Anzahl der gegenueber vertikalen
                             !feuchtadiabatischen Verrueckungen
                             !invarianten Modellvariablen 

!             Zeiger fuer die Variablen :

              tet_l=1,     & !feucht-potentielle Temperatur
              h2o_g=2,     & !Gesamtwasseergehalt
              u_m=3,       & !zonale Geschw.komp. im Massenzentrum
              v_m=4          !meridionale  ,,      ,,     ,, 

!     Beachte:tet_l,h2o_g muessen in [1,ninv]; u_m,v_m in [ninv+1,nred],
!             liegen!

! Local scalars:

      LOGICAL lini           ! initialization step at model start

      INTEGER (KIND=iintegers) :: &

              i_st,i_en,    & !Start- und Endindices in zonale Richtung
              j_st,j_en,    & !Start- und Endindices in merid. Richtung
              i,j,k,        & !Diskretis.index fuer lamda, phi, sigma
                  ii,jj,    & !Indices fuer diverse Schleifen
              ntim1,        & !Index fuer die gueltige Zeitstufe der TKE
              it_start,     & !Startindex der Iterationen
              it_durch        !Durchgangsindex der Iterationen

      REAL (KIND=wp)     :: &

!          Hilfsvariablen:

           wert, val1, val2, & ! Platzhalter fuer beliebige Zwischenergebnisse
           fakt,             & !  ,,         ,,     ,,      Faktoren

!     Platzh. fuer therm. und mech. Antrieb der Turbulenz in (1/s)**2 (fh2,fm2):

           fh2,fm2, &

!     Platzh. fuer horiz. Geschw.-Komponenten und Geschw.betrag:

           vel1,vel2,velh, &

!     Platzh. fuer SQRT(2*TKE)-Werte:

           q0,q1,q2,q3, &

!     Platzh. fuer Hoehendifferenzen und Laengaenskalen:

           dh,l_turb,lh,lm,z0d,z_surf,len1,len2, &
           dz_sg_m, dz_sg_h, dz_g0_m, dz_g0_h, &
           dz_s0_m, dz_sa_m, &
           h_2m, h_10m, a_top, a_atm, a_2m, a_10m, &

!     sonstiges:

           a_tet,a_vap,xf,wf, &
           uk1,vk1,uk,vk, &
           rat_m,rat_h,fac_m,fac_h, fr_sd_h

      REAL (KIND=wp)     :: &
        pk1, dzg, x1, x2, alf, rc_store, ex_store, tet_store, &
        virt_2d  (ie,je)


! Local arrays:

      LOGICAL            :: lo_ice(ie,je) ! logical sea ice indicator

!     Lokale Hilfsfelder:

!US for a better vectorization (hopefully)
      REAL (KIND=wp)     ::  &
        l_tur_z0 (ie,je),    &

        h_top_2d (ie,je),    &
        h_atm_2d (ie,je),    &

        z0m_2d   (ie,je),    &
        z0d_2d   (ie,je),    &
        z2m_2d   (ie,je),    &
        z10m_2d  (ie,je),    &

        hk_2d    (ie,je),    &
        hk1_2d   (ie,je),    &
        h_can_2d (ie,je),    &

        rat_m_2d (ie,je),    &
        rat_h_2d (ie,je),    &
        fac_h_2d (ie,je),    &
        fm2_2d   (ie,je),    &
        fh2_2d   (ie,je),    &
        frc_2d   (ie,je),    &

        vel1_2d  (ie,je),    &
        vel2_2d  (ie,je),    &
        vel_2d   (ie,je),    &
        ta_2d    (ie,je),    &
        qda_2d   (ie,je),    &
        ts_2d    (ie,je),    &
        qds_2d   (ie,je),    &

        dz_0a_m  (ie,je),    &
        dz_0a_h  (ie,je),    &
        dz_sa_h  (ie,je),    &
        dz_s0_h  (ie,je)

INTEGER (KIND=iintegers) ::  &
        k_2d     (ie,je)

      REAL (KIND=wp)     :: &

!     Vertikale Gradienten verschiedener thermodyn. Variablen:

           grad(ie,je,nred), &


!     Zwischenspeicher fuer Wolkenwasser und Bedeckungsgrad:

           clc(ie,je,ke:ke),clcw(ie,je,ke:ke)

REAL (KIND=wp),     TARGET :: &
          dd(ie,je,0:7)    ! local derived turbulence parameter

!     Optimisation variables

INTEGER (KIND=iintegers) :: idone
LOGICAL notdone(iendpar,jendpar)


INCLUDE 'statement_functs.incf'

!- End of header --------------------------------------------------------------
!------------------------------------------------------------------------------

! 1)  Vorbereitungen:

      ierrstat=0
      yerrormsg = ''

 IF (iini.GT.0) THEN !an initialization run
    lini=.TRUE.
    IF (iini.EQ.1) THEN !separate initialization before the time loop
       it_start=1 !only 'it_end' iterations for initialization
                  !and an additional one at the first time loop
    ELSE !initialization within the first time step
       it_start=0 !"it_end+1" iterations for initializatzion
    END IF
 ELSE !not an initialization run
    lini=.FALSE.
    it_start=it_end !only one iteration
 END IF

!     Unterste Hauptflaeche halbiert den Abstand zur
!     untersten Nebenflaeche:

      xf=z2 
      wf=xf-z1
      xf=z1/xf

!     Festlegung der horiz. Schleifengrenzen: 

      i_st=istartpar 
      i_en=iendpar

      j_st=jstartpar
      j_en=jendpar

!     Berechnung abgeleiteter Parameter:

      CALL turb_param  (istartpar, iendpar, jstartpar, jendpar, grav, cp_d, dd)

! 2)  Initialisierung der z0-Werte ueber Meer
!     und der laminaren Transferfaktoren: 

      ! Set the logical mask lo_ice to distinguish between ice covered
      ! and open water sea or lake grid points.

      DO j=j_st,j_en
        DO i=i_st,i_en

          IF (fr_land(i,j).LT.z1d2) THEN
            ! Water point.
            IF (.NOT. lseaice) THEN
              ! Sea ice model is not used.
              ! Ice surface if SST is less than the salt water freezing temperature.
              lo_ice(i,j) = t_g(i,j) < t0_melt + zt_ice
            ELSE
              ! Sea ice model is used.
              ! Ice surface if ice is present.
              lo_ice(i,j) = h_ice(i,j) > 0.0_wp
            END IF
            IF (llake) THEN
              ! Lake model is used.
              ! Ice surface if this is a lake point AND ice is present.
              IF ((depth_lk(i,j) > 0.0_wp) .AND. (h_ice(i,j) >= h_Ice_min_flk)) &
              lo_ice(i,j) = .TRUE.
            END IF
          END IF

        END DO
      END DO

      IF (lini) THEN

         DO j=j_st,j_en
         DO i=i_st,i_en

            IF (fr_land(i,j).LT.z1d2) THEN

!              Ueber Wasserpunkten:

              ! Use ice surface roughness or open-water surface roughness
              ! according to lo_ice

               IF ( lo_ice(i,j) ) THEN

!                 Bei Eisdecke: 

                  gz0(i,j)=grav*z0_ice
               ELSE

!                 Bei von Schubspannung abhaengiger Wellenhoehe:

!                 Einfachste Schaetzung der Schubspannung als Impusls-
!                 flussdichte durch die Nebenflaeche ke mit Hilfe
!                 einer diagnostischen TKE ohne Beruecksichtigung von
!                 Feuchte-Effekten und mit neuchtralen Stabilitaets-
!                 funktion und Anwendung der Charnockflormel:

                  ii=MAX( i-vst, 1 )
                  jj=MAX( j-vst, 1 )

                  l_turb=hhl(i,j,ke)-hhl(i,j,ke1)
                  l_turb=akt*MAX( len_min, l_turb/(z1+l_turb/l_scal) )
!test
!                 l_turb=akt*MIN( l_scal, hhl(i,j,ke)-hhl(i,j,ke1) )
!test
                  dh=hhl(i,j,ke-1)-hhl(i,j,ke1)

                  vel1=(u(i,j,ke-1)+u(ii,j,ke-1))
                  vel2=(u(i,j,ke)+u(ii,j,ke))
                  grad(i,j,u_m)=(vel1-vel2)/dh

                  vel1=(v(i,j,ke-1)+v(i,jj,ke-1))
                  vel2=(v(i,j,ke)+v(i,jj,ke))
                  grad(i,j,v_m)=(vel1-vel2)/dh

                  grad(i,j,tet_l)=z2*(t(i,j,ke-1)-t(i,j,ke))/dh &
                                 +tet_g

                  fm2=grad(i,j,u_m)**2+grad(i,j,v_m)**2
                  fh2=grav*grad(i,j,tet_l)/t(i,j,ke)

                  ! Vereinfachte Loesung mit Rf=Ri:
                  IF (fh2.GE.(z1-rim)*fm2) THEN
                     ! Die krit. Ri-Zahl wird ueberschritten und lm, sowie lh
                     ! werden durch lm bei der krit. Ri-Zahl angenaehert:
                     fakt=z1/rim-z1
                     lm=l_turb*(b_2-(a_6+a_3)*fakt)
                     lh=lm
                  ELSE
                     fakt=fh2/(fm2-fh2)
                     lm=l_turb*(b_2-(a_6+a_3)*fakt)
                     lh=l_turb*(b_1-a_5*fakt)
                  END IF

                  val1=lm*fm2; val2=lh*fh2
                  wert=MAX( val1-val2, rim*val1 )

                  q0=SQRT(d_m*l_turb*wert)

                  tke(i,j,ke,nvor)=q0

                  wert=lm*q0*SQRT(fm2)
                  gz0(i,j)=MAX( grav*len_min, alpha0*wert )
!                 gz0(i,j)=MAX( grav*len_min, &
!                               alpha0*wert+alpha1*grav*con_m/SQRT(wert) )
               END IF
            END IF

            tkvm(i,j,ke)=con_m
            tkvh(i,j,ke)=con_h
            tkvm(i,j,ke1)=con_m
            tkvh(i,j,ke1)=con_h

         END DO
         END DO
      END IF


! 3)  Diagnose des Fluessigwassergehaltes und des Bedeckungsgrades
!     in der unterstersten Modellschicht:

      CALL cloud_diag(clc,clcw,                                        &
           i_st,i_en,j_st,j_en,ke,ke,                                  &
           1,ie,1,je,ke,ke,                                            &
           ie,je,ke,ke1,                                               &
           rdv,o_m_rdv,rvd_m_o,lhocp,t0_melt,                          &
           b1,b2w,b3,b4w,b234w,b2i,b4i,                                &
           uc1,uc2,ucl, clc_diag, q_crit,                              &
           t(:,:,:),qv(:,:,:),qc(:,:,:),prs(:,:,:), rcld,ps(:,:),      &
           itype_wcld)

! 4)  Berechnung der Transferkoeffizienten:

      DO it_durch=it_start, it_end !Iterationen

      DO j=j_st,j_en
         DO i=i_st,i_en

!           Berechnung der benoetigten Laengenskalen:

            ! Dicke der Modell-Prandtl-Schicht
            h_top_2d(i,j) = hhl(i,j,ke)-hhl(i,j,ke1)
            h_atm_2d(i,j) = h_top_2d(i,j)*xf

            ! Hoehe des 2m- und 10m-Niveaus
            h_2m  = z2
            h_10m = z10

            ! Rauhigkeitslaenge
            z0m_2d(i,j) = gz0(i,j)/grav
            z_surf      = z0m_2d(i,j)/sai(i,j)

            ! turbulente Laengenskala
            l_tur_z0(i,j) = akt*z0m_2d(i,j)

            ! turbulente Distanz auf der untersten Hauptflaeche
            a_atm = h_atm_2d(i,j)+z0m_2d(i,j)

            ! turbulente Distanz auf der untersten Nebenflaeche
            a_top = h_top_2d(i,j)+z0m_2d(i,j)

!           Laminare Korrektur der Diffusionskoeffizienten:
            tkvm(i,j,ke1)=MAX( con_m, tkvm(i,j,ke1) )
            tkvh(i,j,ke1)=MAX( con_h, tkvh(i,j,ke1) )

            fakt=z1+(z1-REAL(NINT(fr_land(i,j)),wp))*(rat_sea-z1)

            rat_m=tkvm(i,j,ke1)/con_m
            rat_h=tkvh(i,j,ke1)/con_h

!           Berechnung der effektiven Widerstandslaenge der L-Schicht:
            dz_sg_m=rlam_mom*z_surf
            dz_sg_h=fakt*rlam_heat*z_surf*(rat_h/rat_m)

!           Berechnung weiterer effektiver Widerstandslaengen fuer Skalare:

!           Bestandesschicht ohne lam. Grenzschicht:
            dz_g0_h=z_surf*LOG(rat_m)

!           Bestandesschicht inclusive lam. Grenzschicht:
            dz_s0_h(i,j)=dz_sg_h+dz_g0_h

!           Berechnung der effektiven Bestandeshoehe:
            IF (dz_sg_h.eq.z0) THEN
               h_can_2d(i,j)=rat_can*z0m_2d(i,j)
            ELSE
               h_can_2d(i,j)=rat_can*dz_s0_h(i,j)*LOG(dz_s0_h(i,j)/dz_sg_h)
            END IF

!           Berechnung weiterer effektiver Widerstandslaengen fuer Impuls:

!           Widerstandslaengen inclusive lam. Grenzschicht:
            wert=z1d2*dz_sg_m
            dz_s0_m=wert+SQRT(wert**2+h_can_2d(i,j)*dz_sg_m)

!           Widerstandslaengen ohne lam. Grenzschicht:
            dz_g0_m=dz_s0_m-dz_sg_m

!           der turb. Prandtl-Schicht:

            rat_m=(tkvm(i,j,ke)*z0m_2d(i,j))/(tkvm(i,j,ke1)*a_top)
            rat_m_2d(i,j)=MIN( z2, MAX( z1d2, rat_m ) )

            fac_m=(rat_m_2d(i,j)-z1)*z0m_2d(i,j)/h_top_2d(i,j)
            IF (fac_m.EQ.z1) THEN
               dz_0a_m(i,j)=z0m_2d(i,j)*h_atm_2d(i,j)/a_atm
            ELSE
               dz_0a_m(i,j)=z0m_2d(i,j)*LOG(a_atm/(z0m_2d(i,j)+fac_m*h_atm_2d(i,j)))/(z1-fac_m)
            END IF

            rat_h=(tkvh(i,j,ke)*z0m_2d(i,j))/(tkvh(i,j,ke1)*a_top)
            rat_h_2d(i,j)=MIN( z2, MAX( z1d2, rat_h ) )

            fac_h_2d(i,j)=(rat_h_2d(i,j)-z1)*z0m_2d(i,j)/h_top_2d(i,j)
            IF (fac_h_2d(i,j).EQ.z1) THEN
               dz_0a_h(i,j)=z0m_2d(i,j)*h_atm_2d(i,j)/a_atm
            ELSE
               dz_0a_h(i,j)=z0m_2d(i,j)*LOG(a_atm/(z0m_2d(i,j)+fac_h_2d(i,j)*h_atm_2d(i,j))) &
                                      /(z1-fac_h_2d(i,j))
            END IF

!           von den Oberflaechen bis zum Oberrand der Prandtl-Schicht
!           (unterste Modell-Hauptflaeche):
            dz_sa_m      = dz_s0_m      + dz_0a_m(i,j)
            dz_sa_h(i,j) = dz_s0_h(i,j) + dz_0a_h(i,j)

!           Reduktionsfaktoren fuer die Bestandesschicht
!           incl. lam. Grenzschicht:

            tfm(i,j)=dz_0a_m(i,j)/dz_sa_m
            tfh(i,j)=dz_0a_h(i,j)/dz_sa_h(i,j)

!           Reduktionsfaktor fuer die Verdunstung aufgrund eines um den
!           Faktor 'rat_lam' gegenueber fuehlbarer Waerme vergroesserten
!           laminaren Transpostwiderstandes:

            tfv(i,j)=z1/(z1+(rat_lam-z1)*dz_sg_h/dz_sa_h(i,j))
         END DO
      END DO

      DO j=j_st,j_en
         DO i=i_st,i_en
!           Berechnung der Windgeschwindigkeiten an den Massepunkten:

            ii=MAX( i-vst, 1 )
            jj=MAX( j-vst, 1 )

            vel1_2d(i,j)=(u(i,j,ke)+u(ii,j,ke))*z1d2
            vel2_2d(i,j)=(v(i,j,ke)+v(i,jj,ke))*z1d2

!           Berechnung der korrigierten Prandtl-Schicht-Werte:

            IF (lprfcor) THEN
               len1=z2*h_top_2d(i,j)
               len2=(h_top_2d(i,j)-h_atm_2d(i,j))**2 &
                 /((hhl(i,j,ke-1)+hhl(i,j,ke))*z1d2-hhl(i,j,ke1)-h_atm_2d(i,j))
               lm=len1-tfm(i,j)*h_atm_2d(i,j)-len2
               lh=len1-tfh(i,j)*h_atm_2d(i,j)-len2

               velh=(u(i,j,ke-1)+u(ii,j,ke-1))*z1d2
               vel1_2d(i,j)=(len1*vel1_2d(i,j)-len2*velh)/lm
               velh=(v(i,j,ke-1)+v(i,jj,ke-1))*z1d2
               vel2_2d(i,j)=(len1*vel2_2d(i,j)-len2*velh)/lm

               ta_2d(i,j)=(len1*t(i,j,ke)-h_atm_2d(i,j)*tfh(i,j)*t_g(i,j) &
                                    -len2*t(i,j,ke-1))/lh
               qda_2d(i,j)=(len1*qv(i,j,ke)-h_atm_2d(i,j)*tfh(i,j)*qv_s(i,j) &
                                      -len2*qv(i,j,ke-1))/lh

               vel_2d(i,j)=MAX(vel_min, SQRT(vel1_2d(i,j)**2+vel2_2d(i,j)**2) )
            ELSE
               ta_2d(i,j)=t(i,j,ke)
               qda_2d(i,j)=qv(i,j,ke)

               vel_2d(i,j)=MAX( vel_min, SQRT(vel1_2d(i,j)**2+vel2_2d(i,j)**2) )
            END IF

         END DO
      END DO

      DO j=j_st,j_en
         DO i=i_st,i_en
!           Berechnung der unteren Randwerte der Prandtl-Schicht:
             ts_2d(i,j)= ta_2d(i,j)*(z1-tfh(i,j))+ t_g(i,j)*tfh(i,j)
            qds_2d(i,j)=qda_2d(i,j)*(z1-tfh(i,j))+qv_s(i,j)*tfh(i,j)
         END DO
      END DO

!US      IF (it_durch.LE.it_end) THEN
!US         !Strange to say this seemingly needless condition makes the NEC-compiler
!US         !to inline the CALLs of thermodynamic FUNCTIONs and finally to do the
!US         !vectorization of the follwing loop!

      DO j=j_st,j_en
         DO i=i_st,i_en

!           Berechnung von thermodynamischen Groessen:

            IF (icldm_tran.LE.0) THEN
              rc_store=z0
            ELSEIF (icldm_tran.EQ.1) THEN
              IF ( qc(i,j,ke) > z0 ) THEN
                rc_store = z1
              ELSE
                rc_store = z0
              ENDIF
            ELSE
              rc_store=clc(i,j,ke)
            ENDIF

            ex_store=zexner(ps(i,j)) !Exner-Faktor
            virt_2d(i,j)=z1/(z1+rvd_m_o*qds_2d(i,j)-clcw(i,j,ke)) !rezip. virtueller Faktor
            tet_store=ts_2d(i,j)/ex_store !pot. Temp

!           Temp.tendenz der Saettigungsfeuchte:
            alf=zdqsdt(ts_2d(i,j),zqvap(zpsat_w(ts_2d(i,j)),ps(i,j)))

!           Berechn. der thermodyn. Faktoren a,b,c:
            x1=z1/(z1+alf*lhocp)
            x2=alf*ex_store

!           Berechn. des thermodyn. Faktors d:
            fakt=lhocp/ts_2d(i,j)-(z1+rvd_m_o)*virt_2d(i,j)

!           Berechn. der thermodyn. Faktoren A_tet und A_q:
            a_tet=grav*(z1/tet_store-rc_store*x1*x2*fakt)
            a_vap=grav*(rvd_m_o*virt_2d(i,j)+rc_store*x1*fakt)

!           Berechnung der benoetigten vertikalen Gradienten und
!           abgeleiteter Groessen:

            grad(i,j,u_m)=tfm(i,j)*vel1_2d(i,j)/dz_0a_m(i,j)
            grad(i,j,v_m)=tfm(i,j)*vel2_2d(i,j)/dz_0a_m(i,j)

!           Beachte: Fuer die Windkompnenten wurden die unteren
!           Randwerte der Prandtl-Schicht nicht berechnet, daher
!           muss der Faktor tfm hier bei grad(u_m) und grad(v_m)
!           auftauchen!

            grad(i,j,tet_l)=((ta_2d(i,j)-ts_2d(i,j))/dz_0a_h(i,j) &
                       +tet_g*(z1+clcw(i,j,ke)*lhocp*virt_2d(i,j)/ts_2d(i,j)))/ex_store
            grad(i,j,h2o_g)=(qda_2d(i,j)-qds_2d(i,j))/dz_0a_h(i,j)

            fm2_2d(i,j)=grad(i,j,u_m)**2+grad(i,j,v_m)**2
            fh2_2d(i,j)=a_tet*grad(i,j,tet_l)+a_vap*grad(i,j,h2o_g)
         END DO
      END DO

      DO j=j_st,j_en
         DO i=i_st,i_en
!           Berechnung der atmosphaerischen Forcing-Funktion:

            IF (it_durch.EQ.0) THEN

!              Startinitialisierung:
               IF (fh2_2d(i,j).GE.(z1-rim)*fm2_2d(i,j)) THEN
                  ! Die krit. Ri-Zahl wird ueberschritten und lm, sowie lh
                  ! werden durch lm bei der krit. Ri-Zahl angenaehert:
                  fakt=z1/rim-z1
                  lm=l_tur_z0(i,j)*(b_2-(a_6+a_3)*fakt)
                  lh=lm
               ELSE
                  fakt=fh2_2d(i,j)/(fm2_2d(i,j)-fh2_2d(i,j))
                  lm=l_tur_z0(i,j)*(b_2-(a_6+a_3)*fakt)
                  lh=l_tur_z0(i,j)*(b_1-a_5*fakt)
               END IF

               val1=lm*fm2_2d(i,j)
               val1=MAX(rim*val1,val1-lh*fh2_2d(i,j))

               q0=MAX(SQRT(d_m*l_tur_z0(i,j)*val1),vel_min)

               ntim1=nvor
            ELSE

!              mit Hilfe der vorhergehenden Werte:
               q0=tke(i,j,ke1,nvor)

               lm=tkvm(i,j,ke1)/q0
               lh=tkvh(i,j,ke1)/q0

               val1=lm*fm2_2d(i,j)
               val1=MAX(rim*val1,val1 -lh*fh2_2d(i,j))

               ntim1=ntur
            END IF

            val1=MIN((z1-c_g)*q0**2/(d_4*l_tur_z0(i,j)),val1)

! 4f)       Bestimmung des neuen SQRT(2*TKE)-Wertes:

            IF (imode_tran.LE.1) THEN
               q3=SQRT(d_m*l_tur_z0(i,j)*val1)
            ELSE
               q1=d_m*l_tur_z0(i,j)/dt_tke
               q2=q0+val1*dt_tke
               q3=q1*(sqrt(z1+z4*q2/q1)-z1)*z1d2
            END IF

            q1=SQRT(l_tur_z0(i,j)*val1*d_4/(z1-c_g))
            q3=MAX(vel_min,q1,q3*(z1-tkesmot)+q0*tkesmot)

!           Was die diversen Beschraenkungen durch MAX- und MIN-
!           Funktionen betrifft, siehe bei den entsprechenden Stellen
!           in 'turbdiff.incf'.

!           val1  wird noch zusaetzlich nach oben beschraenkt, damit
!           dieser Wert im Zusammenhang mit q0 keine Singularitaet
!           in der Stabilitaetsfunktion verursachen wuerde. Dies traegt zur
!           numerischen Stabilitaet der iterativ ueber den Zeitindex erfolgenden
!           Bestimmung der Vertikalgradienten bei.

            tke(i,j,ke1,ntim1)=q3
            edr(i,j,ke1)=q3**3/(d_m*l_tur_z0(i,j))

            frc_2d(i,j)=val1 

         END DO
      END DO

! 4g) Bestimmung der neuen stabilitaetsabh. Laengenskalen:

         CALL stab_funct(sm=tkvm(:,:,ke1), sh=tkvh(:,:,ke1), fm2=fm2_2d, fh2=fh2_2d, &
                         frc=frc_2d, tvs=tke(:,:,ke1,ntim1), tls=l_tur_z0, dd=dd,     &
                         i_st=i_st,i_en=i_en, j_st=j_st,j_en=j_en)

      DO j=j_st,j_en
         DO i=i_st,i_en

! 4h)       Bestimmung der durch Wirkung der L-Schicht
!           korrigierten Diffusionskoeffizienten
!           und der zugehoerigen Transferkoeffizienten:

!           unkorrigierte Diffusionskoeffizienten:

            fakt=l_tur_z0(i,j)*tke(i,j,ke1,ntim1)
            tkvm(i,j,ke1)=fakt*tkvm(i,j,ke1)
            tkvh(i,j,ke1)=fakt*tkvh(i,j,ke1)

!           Belegung der Felder fuer die Transferkoeffizienten:

            tcm(i,j)=tkvm(i,j,ke1)*tfm(i,j)/(dz_0a_m(i,j)*vel_2d(i,j))
            tch(i,j)=tkvh(i,j,ke1)*tfh(i,j)/(dz_0a_h(i,j)*vel_2d(i,j))

! 4i)       Diagnose von gz0 (fuer den naechsten Zeitschritt)
!           ueber Wasserflaechen mit der (angepassten) Charnock-Formel
!           und Einschraenkung von z0m_dia ueber Land:

            IF (fr_land(i,j).LT.z1d2) THEN

              ! Use ice surface roughness or open-water surface roughness
              ! according to lo_ice
               IF ( lo_ice(i,j) ) THEN
                  ! Ice-covered grid box
                  gz0(i,j)=grav*z0_ice
               ELSE
                  velh=(tke(i,j,ke1,ntim1)+tke(i,j,ke,nvor))*z1d2
                  wert=tcm(i,j)*vel_2d(i,j)*SQRT(vel_2d(i,j)**2+velh**2)
                  gz0(i,j)=MAX(grav*len_min,alpha0*wert)
!                 gz0(i,j)=MAX(grav*len_min,alpha0*wert+grav*alpha1*con_m/SQRT(wert))
               END IF
               !Ueber See gibt es keinen synoptischen Garten
               z0d_2d(i,j)=z0m_2d(i,j)
            ELSE
               !Die Rauhigkeitslaenge einer SYNOP Station soll immer
               !kleiner als 10m bleiben:
               z0d_2d(i,j)=MIN(z10,z0m_dia)
            END IF
         END DO
      END DO

!----------------------------------------

! 5)  Diagnose der meteorologischen Groessen im 2m- und 10m-Niveau:

      DO j=j_st,j_en
         DO i=i_st,i_en

            IF     (itype_diag_t2m == 1) THEN
              z2m_2d (i,j) = h_2m
              z10m_2d(i,j) = h_10m
            ELSEIF (itype_diag_t2m == 2) THEN
              z2m_2d (i,j) = z2-h_can_2d(i,j) !2m ueber dem Bodenniveau des Bestandes
              z10m_2d(i,j) = z10-z0m_2d(i,j) !Hoehe in der turbulente Distanz 10m betraegt
            ENDIF

            !Erste Belegung zweier benachbarter Modellniveaus:            

            hk_2d(i,j)=h_atm_2d(i,j)
            hk1_2d(i,j)=0.0_wp
            k_2d(i,j)=ke

         ENDDO
      ENDDO

!     Diagnose der 2m-Groessen:

!US   The following loop has been modified by NEC for a better vectorization
      idone = 0
      DO WHILE (idone < (i_en-i_st+1)*(j_en-j_st+1))
         idone = 0 ! points that are done, are counted again in every do while
                   ! loop. So idone must be tilted, not to forget some points
         DO j=j_st,j_en
            DO i=i_st,i_en
               notdone(i,j) = (hk_2d(i,j)<z2m_2d(i,j).AND.k_2d(i,j)>1)
               IF (notdone(i,j)) THEN
                  k_2d(i,j)=k_2d(i,j)-1
                  hk1_2d(i,j)=hk_2d(i,j)
                  hk_2d(i,j)=(hhl(i,j,k_2d(i,j))+hhl(i,j,k_2d(i,j)+1))*z1d2-hhl(i,j,ke1)
               ELSE
                  idone = idone + 1
               ENDIF
            ENDDO
         ENDDO
      END DO

      DO j=j_st,j_en
         DO i=i_st,i_en

            IF (k_2d(i,j).EQ.ke) THEN

!              2m-Niveau unterhalb der untersten Modell-Hauptflaeche
!              in der lokalen Prandtl-Schicht mit Rauhigkeitslaenge z0d:

               pk1=ps(i,j)

             IF     (itype_diag_t2m == 1) THEN
               z0d=z0d_2d(i,j)
               a_atm=h_atm_2d(i,j)+z0d
               a_2m=h_2m+z0d

!              Dimensionsloser Widerstand des Rauhigkeitsbestandes der SYNOP-Wiese,
!              wobei die turbulente Geschwindigkeitsskala und der Oberflaechenindex
!              gegenueber dem mittleren Rauhigkeitsbestand gleich bleiben:

               fr_sd_h=MAX(z0, dz_s0_h(i,j)/z0m_2d(i,j)+LOG(z0d/z0m_2d(i,j)))

               fac_h=(rat_h_2d(i,j)-z1)*z0d/h_top_2d(i,j)
               IF (fac_h.EQ.z1) THEN
                  fakt=fr_sd_h+h_atm_2d(i,j)/a_atm
                  wert=fr_sd_h+h_2m/a_2m
               ELSE
                  fakt=fr_sd_h+LOG(a_atm/(z0d+fac_h*h_atm_2d(i,j)))/(z1-fac_h)
                  wert=fr_sd_h+LOG(a_2m/(z0d+fac_h*h_2m))/(z1-fac_h)
               END IF

               fakt=wert/fakt

               t_2m(i,j) =  t_g(i,j) + (ta_2d(i,j)-t_g(i,j))*fakt &
                          + tet_g*(h_atm_2d(i,j)*fakt-h_2m)

               qv_2m(i,j)= qv_s(i,j) + (qda_2d(i,j)-qv_s(i,j))*fakt

             ELSEIF (itype_diag_t2m == 2) THEN
               IF (z2m_2d(i,j).LT.z0) THEN

!                 2m-Niveau liegt innerhalb der Bestandesschicht
!                 mit exponentiellen Vertikalprofilen:

                  fakt=dz_s0_h(i,j)/dz_sa_h(i,j)

                  ! In order to avoid overflow of the exp() below, we need
                  ! to make sure the exponent does not exceed ~80 or ~300 for
                  ! single precision or double precision, respectively. Since
                  ! z2m_2d is on the order of 2, we need to limit dz_s0_h
                  ! to values larger than ~0.02 or ~0.002, repectively,
                  ! in order to fulfill this condition.
                  IF (dz_s0_h(i,j) >= 0.02_wp) THEN
                    fakt=MIN(z1,MAX(z0,fakt*EXP(z2m_2d(i,j)/dz_s0_h(i,j))))
                  ELSE
! mr: Der Exponent x geht gegen -inf., und 0<=fakt<0.02, so dass fakt*EXP(x) gegen Null geht
                    fakt=z0
                  ENDIF

               ELSE

!                 2m-Niveau liegt innerhalb der Modell_Prandtl-Schicht
!                 mit logarithmischen Vertikalprofilen:

                  IF (ABS(z1-fac_h_2d(i,j)) < epsi ) THEN
                     wert=z0m_2d(i,j)*z2m_2d(i,j)/(z2m_2d(i,j)+z0m_2d(i,j))
                  ELSE
                     wert=z0m_2d(i,j)*LOG((z2m_2d(i,j)+z0m_2d(i,j))/            &
                         (z0m_2d(i,j)+fac_h_2d(i,j)*z2m_2d(i,j)))/(z1-fac_h_2d(i,j))
                  END IF
                  fakt=(dz_s0_h(i,j)+wert)/dz_sa_h(i,j)

               END IF

               t_2m(i,j) =  t_g(i,j) + (ta_2d(i,j)  -  t_g(i,j)) * fakt
               qv_2m(i,j)= qv_s(i,j) + (qda_2d(i,j) - qv_s(i,j)) * fakt

             ENDIF  ! itype_diag_t2m

            ELSE

!              2m-Niveau liegt oberhalb der untersten Hauptflaeche:

               k = k_2d(i,j)
               pk1=prs(i,j,k+1)

             IF     (itype_diag_t2m == 1) THEN
               wert=(h_2m+z0d_2d(i,j))/(hk1_2d(i,j)+z0d_2d(i,j))
               fakt=(hk_2d(i,j)+z0d_2d(i,j))/(hk1_2d(i,j)+z0d_2d(i,j))

               wert=LOG(wert)/LOG(fakt)

               t_2m(i,j) =  t(i,j,k+1) + wert*( t(i,j,k)- t(i,j,k+1)) &
                          + tet_g*( (hk_2d(i,j)-hk1_2d(i,j))*wert - (h_2m-hk1_2d(i,j)) )
               qv_2m(i,j) = qv(i,j,k+1) + wert*(qv(i,j,k)-qv(i,j,k+1))

             ELSEIF (itype_diag_t2m == 2) THEN
               wert=(z2m_2d(i,j)+z0m_2d(i,j))/(hk1_2d(i,j)+z0m_2d(i,j))
               fakt=(hk_2d(i,j)+z0m_2d(i,j))/(hk1_2d(i,j)+z0m_2d(i,j))

               wert=LOG(wert)/LOG(fakt)

                t_2m(i,j) =  t(i,j,k+1) + wert * ( t(i,j,k) -  t(i,j,k+1))
               qv_2m(i,j) = qv(i,j,k+1) + wert * (qv(i,j,k) - qv(i,j,k+1))
             ENDIF
            END IF

!           Berechnung der zugehoerigen Taupunktstemperatur:

            dzg=(z2m_2d(i,j)-hk1_2d(i,j))*grav
            wert=pk1 &
                *EXP(-dzg/(r_d*t_2m(i,j)*(z1+rvd_m_o*qv_2m(i,j))))  !P_2m
            wert=wert/(rdv/MAX(qv_2m(i,j),epsi)-rdv+z1)             !e_2m

            rh_2m(i,j)=100.0_wp*MIN(wert/zpsat_w(t_2m(i,j)),z1)

            fakt=LOG(wert/b1)

            td_2m(i,j)=MIN((b2w*b3-b4w*fakt)/(b2w-fakt),t_2m(i,j))

         ENDDO
      ENDDO

!     Diagnose der 10m-Groessen:

!US   The following loop has been modified by NEC for a better vectorization
      idone = 0
      DO WHILE (idone < (i_en-i_st+1)*(j_en-j_st+1))
         idone = 0 ! points that are done, are counted again in every do while
                   ! loop. So idone must be tilted, not to forget some points
         DO j=j_st,j_en
            DO i=i_st,i_en
               notdone(i,j) = (hk_2d(i,j)<z10m_2d(i,j).AND.k_2d(i,j)>1)
               IF (notdone(i,j)) THEN
                  k_2d(i,j)=k_2d(i,j)-1
                  hk1_2d(i,j)=hk_2d(i,j)
                  hk_2d(i,j)=(hhl(i,j,k_2d(i,j))+hhl(i,j,k_2d(i,j)+1))*z1d2-hhl(i,j,ke1)
               ELSE
                  idone = idone + 1
               ENDIF
            ENDDO
         ENDDO
      END DO

      DO j=j_st,j_en
         DO i=i_st,i_en

            ii=MAX(i-1,1)
            jj=MAX(j-1,1)

            IF (k_2d(i,j).EQ.ke) THEN

!              10m-Niveau unterhalb der untersten Modell-Hauptflaeche
!              in der lokalen Prandtl-Schicht mit Rauhigkeitslaenge z0d:

               z0d=z0d_2d(i,j)
               a_atm=h_atm_2d(i,j)+z0d
               a_10m=h_10m+z0d

               fac_m=(rat_m_2d(i,j)-z1)*z0d/h_top_2d(i,j)
               IF (fac_m.EQ.z1) THEN
                  fakt=h_atm_2d(i,j)/a_atm
                  wert=h_10m/a_10m
               ELSE
                  fakt=LOG(a_atm/(z0d+fac_m*h_atm_2d(i,j)))
                  wert=LOG(a_10m/(z0d+fac_m*h_10m))
               END IF

               fakt=wert/fakt

               u_10m(i,j)=vel1_2d(i,j)*fakt
               v_10m(i,j)=vel2_2d(i,j)*fakt
                 
            ELSE

               k = k_2d(i,j)

!              10m-Niveau oberhalb der untersten Modell-Hauptflaeche:

!US this has been forgotten before
               uk =(u(i,j,k)+u(ii,j,k))*z1d2
               vk =(v(i,j,k)+v(i,jj,k))*z1d2
!US

               uk1=(u(i,j,k+1)+u(ii,j,k+1))*z1d2
               vk1=(v(i,j,k+1)+v(i,jj,k+1))*z1d2

               wert=(h_10m+z0d_2d(i,j))/(hk1_2d(i,j)+z0d_2d(i,j))
               fakt=(hk_2d(i,j)+z0d_2d(i,j))/(hk1_2d(i,j)+z0d_2d(i,j))

               wert=LOG(wert)/LOG(fakt)

               u_10m(i,j)=uk1+wert*(uk-uk1)
               v_10m(i,j)=vk1+wert*(vk-vk1)
            END IF

!-----------------------------------------------------------------------
!           nur im alten 1-d-Modell:

!           Berechnung der der Flussdichten fuer Impuls,
!           fuehlbare und latente Energie und der Ri-Zahl:

!           flimpuls(i,j,ke1)=rho_2d(i,j)*tkvm(i,j,ke1)*SQRT(fm2_2d(i,j))
!           flfbwae(i,j,ke1)=cp_d*wert*tkvh(i,j,ke1)*grad(i,j,tet_l)
!           fllawae(i,j,ke1)=lh_v*wert*tkvh(i,j,ke1)*grad(i,j,h2o_g)
 
!           IF (fh2_2d(i,j).EQ.z0) THEN
!              rizahl(i,j,ke1)=z0
!           ELSE
!              rizahl(i,j,ke1)=fh2_2d(i,j)/MAX(fm2_2d(i,j),epsi*abs(fh2_2d(i,j)))
!           END IF

!           Beachte: flimpuls,flfbwae,fllawae und die Ri-Zahl
!                    werden nur fuer Diagnosezwecke berechnet.
!-----------------------------------------------------------------------

         END DO
      END DO

!SCLM-------------------------------------------------------------------
#ifdef SCLM
      IF (lscm) THEN

         TKE_SCLM%mod(ke1)%val=tke(im,jm,ke1,ntur); TKE_SCLM%mod(ke1)%vst=i_cal

         !Achtung: TKE%mod(ke1)%val zeigt z.Z. noch auf die alte TKE-Zeitstufe.
         !Somit wird also der alte tke-Wert mit dem neuen ueberschrieben,
         !was aber ohne Bedeutung ist, weil ab jetzt der alte tke-Wert nicht
         !mehr benoetigt wird. Beachte, dass die Modelleinheit hier [m/s] ist.

         velh=tkvm(im,jm,ke1)*SQRT(fm2_2d(im,jm))
         wert=tkvh(im,jm,ke1)*fh2_2d(im,jm)

         teta=-tkvh(im,jm,ke1)*grad(im,jm,tet_l) !teta cov. with w
         qvap=-tkvh(im,jm,ke1)*grad(im,jm,h2o_g) !qvap cov. with w
         rhos=ps(im,jm)*virt_2d(im,jm)/(r_d*ts_2d(im,jm))

         UST%mod(ke1)%val=SQRT(velh) ; UST%mod(ke1)%vst=i_cal
         TWA%mod(ke1)%val=teta       ; TWA%mod(ke1)%vst=i_cal

         SHF%mod(0)%val=cp_d*rhos*teta ; SHF%mod(0)%vst=i_cal
         LHF%mod(0)%val=lh_v*rhos*qvap ; LHF%mod(0)%vst=i_cal
       ! LMO%mod(0)%val=velh**3/wert ; LMO%mod(0)%vst=i_cal
      END IF
#endif
!SCLM-------------------------------------------------------------------

  ENDDO

END SUBROUTINE turbtran

!==============================================================================

END MODULE turbulence_tran
