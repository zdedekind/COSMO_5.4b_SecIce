!>
!! Source module for computing diffusion coefficients
!! and implicit vertical diffusion:
!!
!! @par Description of *turb_diffusion*:
!!   This  module calculates the tendencies for turbulent
!!   vertical transport of momentum and heat and the coefficients
!!   for turbulent diffusion as well.
!!
!!   The clousure is made on lever 2.5 (Mellor/Yamada) using a prognostic
!!   TKE-equation and includes the formulation of a flow through a porous
!!   medium (roughness layer)
!!
!!   The turbulence model (with some Prandtl-layer approximations is used
!!   for the calculation of turbulent transfer between atmosphere and the
!!   lower boundary too.
!!
!! The module contains the public subroutines :
!!
!!   turbdiff
!!
!! called from the turbulence interface routine of the model.
!!
!! Current Code Owner: DWD, Matthias Raschendorfer
!!  phone:  +49  69  8062 2708
!!  fax:    +49  69  8062 3721
!!  email:  matthias.raschendorfer@dwd.de
!
! History:
! Version      Date       Name
! ----------   ---------- ----
! V5_4a        2016-05-10 Matthias Raschendorfer
!  Initial Release, based on the turbdiff-module of the non-blocked version
!    Blocked and further developed version of the turbulence model
!    Adopted ICON developments as switchable (hardcoded) options
!    Partly new (more consistent) interpretation of already existing selectors
!    Gradually controlling numerical restrictions
!    Adopted 3D-options from COSMO
!    Correcting the horizontal limit of the turbulent length scale
!    Correction of some bugs of the latest ICON version related to optional 
!      lower concentration condition.
!    Using new turbulent transfer velocities (tvm, tvh, tkm), allowing an 
!      easier formulation of transfer-resistances
!
!! @par Copyright and License
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of this software is hereby granted free of charge for an unlimited
!! time, provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
!-------------------------------------------------------------------------------

MODULE turb_diffusion

!-------------------------------------------------------------------------------

! Modules used:
!-------------------------------------------------------------------------------
! Parameter for precision
!-------------------------------------------------------------------------------

#ifdef __COSMO__
USE kind_parameters, ONLY :   &
#elif defined(__ICON__)
USE mo_kind,         ONLY :   &
#endif
    wp              ! KIND-type parameter for real variables

!-------------------------------------------------------------------------------
! Mathematical and physical constants
!-------------------------------------------------------------------------------

#ifdef __COSMO__
USE data_constants, ONLY : &

! Physical constants and related variables:
! -------------------------------------------

    r_d,          & ! gas constant for dry air
    rdv,          & ! r_d / r_v
    rvd_m_o,      & ! r_v/r_d - 1
    cp_d,         & ! specific heat for dry air
    lh_v,         & ! evaporation heat
    lhocp,        & ! lh_v / cp_d
    con_m,        & ! kinematic vsicosity of dry air (m2/s)
    con_h,        & ! scalar conductivity of dry air (m2/s)
    t0_melt,      & ! absolute zero for temperature (K)
    grav => g,    & ! acceleration due to gravity
!
! Parameters for auxilary parametrizations:
! ------------------------------------------

    b1,           & ! variables for computing the saturation steam pressure
    b2w,          & ! over water (w) and ice (i)
    b3,           & !               -- " --
    b4w             !               -- " --
#endif

#ifdef __ICON__
USE mo_physical_constants, ONLY : &
!
! Physical constants and related variables:
! -------------------------------------------
!
    r_d      => rd,       & ! gas constant for dry air
    rvd_m_o  => vtmpc1,   & ! r_v/r_d - 1
    cp_d     => cpd,      & ! specific heat for dry air
    lh_v     => alv,      & ! evaporation heat
    lhocp    => alvdcp,   & ! lh_v / cp_d
    t0_melt  => tmelt,    & ! absolute zero for temperature (K)
    b3       => tmelt,    & !          -- " --

    rdv, con_m, con_h, grav


USE mo_convect_tables, ONLY : &
!
! Parameters for auxilary parametrizations:
! ------------------------------------------
!
    b1       => c1es,     & ! variables for computing the saturation steam pressure
    b2w      => c3les,    & ! over water (w) and ice (e)
    b4w      => c4les       !               -- " --
#endif

!-------------------------------------------------------------------------------
! From Flake model
!-------------------------------------------------------------------------------

#ifdef __COSMO__
USE data_flake,    ONLY: &
#elif defined(__ICON__)
USE mo_data_flake, ONLY: &
#endif
    h_Ice_min_flk      ! Minimum ice thickness [m]

!-------------------------------------------------------------------------------
! Turbulence data (should be the same in ICON and COSMO)
!-------------------------------------------------------------------------------

USE turb_data, ONLY : &

! Numerical constants and parameters:
! -----------------------------------

    tkhmin,       & ! minimal diffusion coefficients for heat
    tkmmin,       & ! minimal diffusion coefficients for momentum
    tkhmin_strat, & ! additional minimal diffusion coefficients for heat for stratosphere
    tkmmin_strat, & ! additional minimal diffusion coefficients for momentum for stratosphere
    ditsmot,      & ! smoothing factor for direct time-step iteration
    tndsmot,      & ! vertical smoothing factor for diffusion tendencies
    frcsmot,      & ! vertical smoothing factor for TKE forcing
    epsi,         & ! relative limit of accuracy for comparison of numbers
    it_end,       & ! number of initialization iterations (>=0)

! Parameters describing physical properties of the lower boundary 
! of the atmosphere:
!---------------------------------------------------------------

    rlam_mom,     & ! scaling factor of the laminar boudary layer for momentum
    rlam_heat,    & ! scaling factor of the laminar boudary layer for heat
    rat_can,      & ! factor for the canopy height
    rat_sea,      & ! ratio of laminar scaling factors for heat over sea and land
    rat_lam,      & ! ratio of laminar scaling factors for vapour and heat
    z0m_dia,      & ! roughness length of a typical synoptic station
    alpha0,       & ! lower bound for Charnock-parameter
    alpha1,       & ! parameter scaling the molek. roughness of water waves
    c_lnd,        & ! surface area index of the land exept the leaves
    c_soil,       & ! surface area index of the (evaporative) soil
    c_sea,        & ! surface area index of the waves over sea
    e_surf,       & ! exponent to get the effective surface area
    zt_ice,       & ! freezing temperature of sea ice
    z0_ice,       & ! roughness length of sea ice
    tur_len,      & ! maximal turbulent length scale [m]
    pat_len,      & ! lenth scale of subscale patterns over land [m]
    len_min,      & ! minimal turbulent length scale [m]
    vel_min,      & ! minimal velocity scale [m/s]
    akt,          & ! von Karman-constant
!
!?  a_h=>a_heat,  & ! factor for turbulent heat transport
!?  a_m=>a_mom,   & ! factor for turbulent momentum transport
    d_h=>d_heat,  & ! factor for turbulent heat dissipation
    d_m=>d_mom,   & ! factor for turbulent momentum dissipation
!
    c_diff,       & ! factor for turbulent diffusion of TKE
    a_hshr,       & ! factor for horizontal shear production of TKE
    c_scld,       & ! factor for liquid water flux density in sub grid scale clouds

    ! derived quantities from turb_param
    tet_g, rim, b_m, b_h, sm_0, sh_0,   &
    d_1, d_2, d_3, d_4, d_5, d_6,       &
    a_3, a_5 ,a_6,                      &
    tur_rcpv, tur_rcpl,                 &

    ! used derived types
    modvar, turvar, varprf, & !

! Switches controlling the turbulence model, turbulent transfer and diffusion:
! ----------------------------------------------------------------------------

    ltkesso,      & ! consider SSO-wake turbulence production of TKE
    ltkecon,      & ! consider convective buoyancy production of TKE
    ltkeshs,      & ! consider separ. horiz. shear production of TKE 
    loutshs,      & ! consider separ. horiz. shear production of TKE for output
    lnonloc,      & ! nonlocal calculation of vertical gradients used for turb. diff.
    lprfcor,      & ! using the profile values of the lowest main level instead of
    ltmpcor,      & ! consideration of thermal TKE-sources in the enthalpy budget
    lexpcor,      & ! explicit corrections of the implicit calculated turbul. diff.

!   for semi-implicit vertical diffusion:
    lsflcnd,      & ! lower flux condition
    ldynimp,      & ! dynamical calculation of implicit weights
    lprecnd,      & ! preconditioning of tridiagonal matrix
    lfreeslip       ! free-slip lower boundary condition (enforeced zero-flux condition for
                    ! for all diffused variables, only for idealized test cases)

USE turb_data, ONLY : &

! Selectors controlling the turbulence model, turbulent transfer and diffusion:
! ----------------------------------------------------------------------------
!
    imode_tran,   & ! mode of TKE-equation in transfer scheme                    (compare 'imode_turb')
    imode_turb,   & ! mode of TKE-equation in turbulence scheme
                    !  0: diagnostic equation
                    !  1: prognostic equation (default)
                    !  2: prognostic equation (implicitly positive definit)
    icldm_tran,   & ! mode of water cloud representation in transfer parametr.   (compare 'icldm_tran)
    icldm_turb,   & ! mode of water cloud representation in turbulence parametr.
                    ! -1: ignoring cloud water completely (pure dry scheme)
                    !  0: no clouds considered (all cloud water is evaporated)
                    !  1: only grid scale condensation possible
                    !  2: also sub grid (turbulent) condensation considered
    itype_sher,   & ! type of shear production for TKE
                    ! 0: only vertical shear of horizontal wind
                    ! 1: previous plus horizontal shear correction
                    ! 2: previous plus shear from vertical velocity
    ilow_def_cond,& !type of the default condition at the lower boundary
                    ! 1: zero surface flux density
                    ! 2: zero surface value
    imode_calcirc,& ! mode of treating the circulation term (related to 'pat_len', imode_pat_len')
                    ! 1: explicit calculation of the flux convergence
                    ! 2: quasi implicit treatment by calculation of effective TKE-gradients
    imode_shshear,& ! mode of calculat. the separated horizontal shear mode related to 'ltkeshs', 'a_hshr')
                    ! 1: with a constant lenght scale 
                    ! 2: with a Ri-dependent length sclale correction
    imode_tkvmini,& ! mode of calculating the minimal turbulent diff. coeffecients
                    ! 1: with a constant value
                    ! 2: with a stability dependent correction
    imode_rat_sea,& ! mode of scaling the laminar resistance for heat over sea (related to 'rat_sea')
                    ! 1: constant ratio compared to land surface
                    ! 2: with a correction for a strongly overheated SST
    imode_vel_min,& ! mode of calculating the minimal turbulent velocity scale (in the surface layer only)
                    ! 1: with a constant value
                    ! 2: with a stability dependent correction
    imode_charpar,& ! mode of estimating the Charnock-Parameter
                    ! 1: use a constant value 
                    ! 2: use a wind-dependent value with a constant lower bound
    imode_syndiag,& ! mode of diagnostics at the synoptic near surface levels (related to 'itype_diag_t2m')
                    ! 1: direct interpolation of temperature and specific humidity
                    ! 2: interpol. of conserved quantities and subsequent statistical saturation adjustm.,
                    !    allowing particularly for the diagnostic of cloud water at the 2m-level (fog)
    imode_trancnf,& ! mode of configuring the transfer-scheme 
                    ! 1: old version: start. with lamin. diffus.; with a lamin. correct. for profile-funct.;
                    !    interpol. T_s rather then Tet_l onto zero-level; calcul. only approx. Tet_l-grads.;
                    !    using an upper bound for TKE-forcing; without transmit. skin-layer depth to turbul.
                    ! 2: 1-st ConSAT: start. with estim. Ustar, without a laminar correct. for prof.-funct.;
                    !    interpol. Tet_l onto zero-level; calcul. Tet_l-gradients directly; 
                    !    without an upper bound for TKE-forcing; with transmit. skin-layer depth to turbul.
                    ! 3: 2-nd ConSAT: as "2", but with a hyperbolic interpol. of profile function
                    !    for stable stratification
                    ! 4: 3-rd ConSAT: as "3", but without using an upper interpolation node
    imode_tkediff,& ! mode of implicit TKE-Diffusion (related to 'c_diff')
                    ! 1: in terms of q=SQRT(2*TKE))
                    ! 2; in terms of TKE=0.5*TKE**2
    imode_adshear,& ! mode of considering additional shear by scale interaction
                    ! 1: not considered for stability functions
                    ! 2:     considered for stability functions
    imode_tkemini,& ! mode of fixing a lower limit of q=2TKE**2
                    ! 1: by using 'vel_min' only
                    ! 2: by adapting to minimal diffusion coefficients
    imode_lamdiff   ! mode of considering laminar diffusion at surface layer
                    ! 1: only when calculating the profile functions
                    ! 2: surface-layer diffusion coeff. always at least at laminar value

USE turb_data, ONLY : &
    mom     ,     & ! index for a momentum variable
    sca     ,     & ! index for a scalar   variable
    ! integers: indices for special layers in some variables (a, dicke)
    nvel    ,     & ! Geschwindigkeitskomponenten
    naux    ,     & !number of auxilary variables
    nmvar   ,     & !
    ntyp    ,     & ! Anzahl von Variablentypen (mom) und (sca)
    u_m     ,     & ! zonale Geschw.komp. im Massenzentrum
    v_m     ,     & ! meridionale  ,,      ,,     ,,
    tet_l   ,     & ! feucht-potentielle Temperatur
    tet     ,     & ! potentielle Temperatur
    tem     ,     & ! Temperatur
    h2o_g   ,     & ! Gesamtwasseergehalt
    vap     ,     & ! Wasserdampfmischungsverh.
    liq     ,     & ! Fluessigwasser  ,,
    ! and now the fields
    l_scal  ,     & !
    fc_min  ,     &
    tmps    ,     & ! surface temperature
    vaps    ,     & ! surface specific humidity
    prss    ,     & ! surface pressure
    eprs    ,     & ! surface Exner-factor
    liqs    ,     & ! liquid water content at the surface
    grad    ,     & ! vertical gradient
    vari    ,     & ! reduced set of variables in the first part of 'turbdiff'
                    ! and later their effective vertical gradients
    dicke   ,     & ! (effektive) Dicke einer Modellschicht
                    ! bzw. sonstiges Hilfsfeld, bei k=ke1 steht
                    ! eine effekt. Dicke der Modell-Prandtlsch.

    len_scale,    & ! turbulent length-scale
    hor_scale,    & ! effective hoprizontal length-scale used for sep. horiz. shear calc.
    rhon    ,     & ! Luftdichte auf Nebenflaechen (einschl. surface)

    shv     ,     & ! velocity scale of the separated horiz. shear mode
    frh     ,     & ! thermal forcing (1/s2) or thermal acceleration (m/s2)
    frm     ,     & ! mechan. forcing (1/s2) or mechan. accelaration (m/s2)
    ftm     ,     & ! mechan. forcing (1/s2) by pure turbulent shear 

    a       ,     & ! a() enthaelt zu beginn die Werte der 5 thermodynamischen
                    ! Koeffizienten dQs/dT,ex_fakt,cp_fakt,g_tet,g_vap
                    ! auf den Positionen 1 bis 5 des letzten Index.
                    !
                    ! Im Falle der impliziten Berechnung von Flussdivergenzen
                    ! enthaelt a() hierfuer benoetigte Hilfsgroessen
     hlp    ,     & !
     can    ,     & !
! Hilfsfelder fuer eine und zwei Variablenschichten:
     lay    ,     & !
     lays   ,     & !
     vars   ,     & !

     dzsm   ,     & !
     dzsh   ,     & ! effektive Dicke der Prandtl-Schicht
     src    ,     & ! arbitrary source term

! Seicher-Target fuer obige Pointer:
     rho_tar    , & !
     exner_tar  , & !
     diss_tar   , & !
     shflu_s_tar, & !
     qvflu_s_tar, & !
     cbig_tar   , & !
     csml_tar   , & !
     rair_tar        

!-------------------------------------------------------------------------------
! Control parameters for the run
!-------------------------------------------------------------------------------

! ICON data have to be declared for these variables, which is done later on
!-------------------------------------------------------------------------------

#ifdef __COSMO__
USE data_fields, ONLY :   &
    impl_weight => a1t    ! implicit weights for tridiagonal solver

USE data_runcontrol, ONLY:   &

! Switches controlling other physical parameterizations:

    lscm,            & ! if Single Column Model is used (default: FALSE)
    l3dturb,         & ! a model run with 3D-(turbulent)-diffusion
    lsso,            & ! SSO-Scheme is active
    lconv              ! confection scheme is active
#endif

#ifdef __ICON__
USE turb_data,         ONLY:   &
    lscm,           & ! if Single Column Model is used (default: FALSE)
    impl_weight,    & ! vertical column with implicit weights for tridiagonal solver
    l3dturb,        & !
    lsso,           & ! SSO-Scheme is active
    lconv             ! confection scheme is active
#endif

USE turb_utilities,          ONLY:   &
!US turb_param,                      &
    adjust_satur_equil,              &
    solve_turb_budgets,              &
    vert_grad_diff,                  &
    prep_impl_vert_diff,             &
    calc_impl_vert_diff,             &
    vert_smooth,                     &
    bound_level_interp,              &
    zexner
    
!-------------------------------------------------------------------------------
#ifdef SCLM
USE data_1d_global, ONLY : &
!
    lsclm, lsurflu, latmflu, i_cal, i_upd, i_mod, im, &
!
    UUA, VVA, UWA, VWA, WWA, UST, TWA, QWA, TTA, TQA, QQA, SHF, LHF, &
    TKE_SCLM=>TKE, BOYPR, SHRPR, DISSI, TRANP
#endif
!SCLM---------------------------------------------------------------------------

!===============================================================================

IMPLICIT NONE

PUBLIC  :: turbdiff

!===============================================================================

REAL (KIND=wp), PARAMETER :: &
    z0 = 0.0_wp,    &
    z1 = 1.0_wp,    &
    z2 = 2.0_wp,    &
    z3 = 3.0_wp,    &
    z4 = 4.0_wp,    &
    z5 = 5.0_wp,    &
    z6 = 6.0_wp,    &
    z7 = 7.0_wp,    &
    z8 = 8.0_wp,    &
    z9 = 9.0_wp,    &
    z10=10.0_wp,    &

    z1d2=z1/z2     ,&
    z1d3=z1/z3     ,&
    z2d3=z2/z3     ,&
    z3d2=z3/z2

REAL (KIND=wp) :: xx

INTEGER :: &
    istat=0, ilocstat=0

LOGICAL :: &
    lerror=.FALSE.

!===============================================================================

CONTAINS

!===============================================================================


SUBROUTINE turbdiff ( &
!
          iini, lturatm,          lstfnct,          ltkeinp,          &
          itnd, lum_dif, lvm_dif, lscadif,          lsfluse, lqvcrst, &
!
          dt_var,dt_tke, nvor, ntur, ntim, &
!
          nvec, ke, ke1, kcm, iblock, ivstart, ivend, &
!
          l_hori, hhl, fr_land, dp0, trop_mask, &
!
          gz0, l_pat, c_big, c_sml, r_air, &
!
          t_g, qv_s, ps, &
          u, v, w, t, qv, qc, prs, rho, epr, &
!
          ptr, ndtr, ndiff, &
!
          tvm, tvh, tfm, tfh, tfv, tkr, tkred_sfc, &
          tke, tkvm, tkvh, rcld, tkhm, tkhh, &
          hdef2, hdiv, dwdx, dwdy,           &
!
          edr, tket_sso, tket_conv, tket_hshr, &
          u_tens, v_tens, t_tens, &
          qv_tens, qc_tens, &
          tketens, tketadv, &
          qv_conv, ut_sso, vt_sso, &
!
          shfl_s, lhfl_s, qvfl_s, umfl_s, vmfl_s, &
!
          ierrstat, yerrormsg, yroutine)

!-------------------------------------------------------------------------------
!
! 
! All tendency parameters are OPTIONAL (except 'tketens' in case of "lturatm=T". If they are missing,
!  calculated tendencies of SUB 'turbdiff' are automatically added to the related prognostic variables.
! It is also possible to use only one time level for TKE using "ntim=1" and thus "nprv=1=ntur".
!
! Description:
!
!     Es werden die Diffusionskoeffizienten berechnet und ggf. Anteile
!     der zeitlichen Tendenzen der turbulenten Diffusion bestimmt
!     und zu den Tendenzfeldern hinzuaddiert.
!     Optional wird eine explizite oder (teil-)implizite Berechnung der
!     Diffusionstendenzen oder aber nur eine Berechnung der Diffusions-
!     koeffizienten durchgefuehrt. Im letzten Fall wird dann ein
!     implizit zu berechnender Anteil der Diffusionstendenzen an
!     anderer Stelle (slow_tendencies) bestimmt.
!     Allerdings koennen dann zusaetzliche explizite Korrekturtendenzen
!     hier in tubdiff bestimmt werden.
!
! Method:
!
!     Die Berechnung basiert auf einer Schliessung 2-ter Ordnung auf
!     dem level 2.5 (nach Mellor/Yamada). Demnach wird also eine
!     prognostische Gleichung fuer die TKE geloest.
!     Ausser der TKE-Advektion, die zusammen mit den Advektionstendenzen
!     der anderen prognostischen Variablen an anderer Stelle berechnet
!     wird, geschieht die gesamte TKE-Prognose in diesem Unterprogramm.

!     Die Formulierung des Schemas erfolgt mit thermodynamischen
!     Variablen, die bei feuchtadiabatischen vertikalen Verrueckungen
!     erhalten bleiben (pot. Fluessigw.temp. und Gesamtwassergehalt),
!     so dass der Kondesationseffekt auf subskalige Vertikalbewegungen
!     beruecksichtigt wird.
!     Die turbulenten Flussdichten der Erhaltungsgroessen werden in
!     solche der Modellvariablen konvertiert, so dass die thermodyn.
!     Kopplung der Flussdichten richtig erhalten bleibt.

!     Angeschlossen ist auch ein optionales statistisches Wolkenschema
!     (nach Sommeria und Deardorff), sub turb_cloud, welches auch
!     subskalige Bewoelkung mit Hilfe der ueber das Feld rcld ausge-
!     gebenen Standardabweichung des Saettigungsdefizites berechnet.

!     Das Turbulenzschema wurde so verallgemeinert, dass es auch bei
!     einer vertikal aufgeloesten Bestandesschicht gueltig ist, indem
!     idealisierend von der Durchstroemung eines poroesen Mediums
!     ausgegangen wird. Die Bilanzgleichungen 1-ter und 2-ter Ordnung
!     enthalten dann zusaetzliche Terme, welche die Wechselwirkungen
!     mit der Bestandes-Matrix beschreiben. Dies wirkt sich zum einen
!     auf die Stabilitaetsfunktionen und zum anderen vor allem auf die
!     TKE-Gleichung aus, welche einen auf den Formwiderstand der
!     Bestandeselemente zurueckzufuehrenden zusaetzlichen Quellterm
!     (Nachlaufturbulenz) enthaelt. Ausserdem werden die turbulenten
!     Flussdichtedivergenzen noch um einen Zusatzterm, welcher der
!     Reduktion des lufterfuellten Volumens im Gitterelement Rechnung
!     traegt, erweitert. Der Effekt des Formwiderstandes in der
!     Impulsgleichung ist ebenfalls beruecksichtigt. Die zusaetzlichen
!     Tendenzterme, die auf die Flussdichten zwischen Bestandes-Matrix
!     und umgebender Luft zurueckzufuehren sind (Bestandesquellen),
!     muessen noch in einem separaten Bestandesmodell parametrisiert
!     werden und sind nicht Gegenstand des Turbulenzmodells.
!     Schliesslich wird auch der Effekt der Transformation von Turbulenz
!     auf der dominierenden Skala in kleinskalige dissipative Turbulenz
!     durch Wirbelbrechen an Koerpern mit
!          Laengenskalen der Abmessungen << turbulente Laengenskala
!     beruecksichtigt; was sich durch eine (von der Laengenskala und
!     Volumendichte jener sehr kleinen Bestandeselemente aubhaengige)
!     Modifikation der Modellkonstanten ausdruecken laesst.

!     Es wird auch versucht den Effekt thermisch induzierter
!     Zirkulationen auf die TKE-Produktion zu beruecksichtigen
!     Hierdurch wird (vor allem) der Austauch in der naechtlichen
!     Grenzschicht erhoeht, was der Tendenz des alten Schemas,
!     in Bodennaehe zu kalte und nicht schnell genug anwachsende
!     Inversionen zu produzieren, entgegenwirkt.

!     Optional kann die Berechnung der vertikalen Gradienten durch eine
!     nicht-lokale Variante erfolgen. Hierbei werden die Gradienten
!     mit Profilen gebildet, die mit einem ueber die stabilitaets-
!     abhaengige Laengenskala gebildeten gleitenden Mittel behandelt
!     wurden.

!     Die Bildung der Anteile der Diffusionstendenzen, die numerisch
!     durch die Multiplikation mit einer Tridiagonalmatrix ausdrueckbar
!     sind, kann (neben der expliziten Variante) auch implizit erfolgen
!     (im Falle der Berechnung von nich-lokalen Gradienten und fuer die
!     TKE allerdings nur explizit).
!     Bei expliziter Rechnung ist, um auch bei Zeitschritten von
!     mehreren Minuten numerisch stabil zu bleiben, eine Limitierung
!     der Groesse der Diffusionskoeffezienten und eine im vertikalen Integral
!     quellenfreie numerische Glaettung der Vertikalprofile der Diffusions-
!     tendenzen erforderlich, sowie eine teilimplizite Behandlung der
!     Diffusionstendenzen in der untersten Modellschicht, erforderlich.

!     Die unteren Randwerte der turbulenten Flussdichten werden ueber
!     die Transferkoeffizienten zwischen Erdboden und unterster
!     Modellschicht (tcm und tch) bestimmt.
!     Optional koennen die Transferkoeffizienten auch mit diesem
!     Unterprogramm bestimmt werden, indem das Turbulenzmodell auch auf
!     das Niveau d+z0 angewandt wird, wobei vertikale Gradienten in
!     diesem Niveau mit Hilfe der Prandtl-Schicht-Hypothese berechnet
!     werden.
!     In diesem Zusammenhang wird auch die Wirkung der laminaren
!     Grenzschicht behandelt.
!
!     Turbulente Horizontaldiffusion (um ein 3-d-Schema zu erhalten)
!     ist noch nicht enthalten, kann aber integriert werden.
!     Uebergabevariablen:
!
!-------------------------------------------------------------------------------

! Declarations
!-------------------------------------------------------------------------------

!Formal Parameters:
!-------------------------------------------------------------------------------

! 0. Parameters controlling the call of 'organize_turbdiff':

LOGICAL, INTENT(IN) :: &

   lstfnct,      & !calculation of stability function required

   lturatm,      & !running turbulence model between atmosph. layers (updating diffusion coefficients)
   lum_dif,      & !running vertical gradient diffusion of horizontal u-momenum
   lvm_dif,      & !running vertical gradient diffusion of horizontal v-momenum
   lscadif,      & !running vertical gradient diffusion of scalar properties

   lsfluse,      & !use explicit heat flux densities at the suface
   lqvcrst,      & !qv-flux-divergence reset requested (only if 'qv_conv' is present)

   ltkeinp         !TKE present as input (at level k=ke1 for current time level 'ntur')

REAL (KIND=wp), INTENT(IN) :: &

   dt_var,       & !time step for ordinary prognostic variables
   dt_tke          !time step for the 2-nd order porgnostic variable 'tke'

INTEGER,        INTENT(IN) :: &

   iini,         & !type of initialization (0: no, 1: separate before the time loop
                   !                             , 2: within the first time step)
   itnd,         & !type of tendency cons. (0: no, 1: in implicit vertical diffusion equation
                   !                               2: by adding to current profile before vertical diffusion
                   !                               3: by using corrected virtual vertical profiles
   ntur,         & !current new time step index of tke
   ntim            !number of tke time levels

INTEGER,        INTENT(IN) :: &

! Horizontal and vertical sizes of the fields and related variables:
! --------------------------------------------------------------------

    nvec,         & ! number of grid points in zonal      direction
    ke,           & ! index of the lowest main model level
    ke1             ! index of the lowest model half level (=ke+1)


INTEGER,        INTENT(INOUT) :: &
    nvor,         & !previous    time step index of tke
    kcm             ! level index of the upper canopy bound

INTEGER,        INTENT(IN) :: &
    iblock

INTEGER,        INTENT(IN) :: &

! Start- and end-indices for the computations in the horizontal layers:
! -----------------------------------------------------------------------

    ivstart,      & ! start index in the nproma vector
    ivend           ! end index in the nproma vector

! Constants related to the earth, the coordinate system
! and the reference atmosphere:
! --------------------------------------------------------------------------

REAL (KIND=wp), DIMENSION(:,:), INTENT(IN) :: &
!
    hhl             ! height of model half levels                   ( m )

REAL (KIND=wp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: &
!
    dp0             ! pressure thickness of layer                   (pa )

REAL (KIND=wp), DIMENSION(:), INTENT(IN) :: &
!
    fr_land,      & ! land portion of a grid point area             ( 1 )
    l_pat  ,      & ! effektive Laengenskala der therm. Inhomogenitaeten 
                    ! der Erdbodenoberflaeche
    l_hori          ! horizontal grid spacing (m)
  
REAL (KIND=wp), DIMENSION(:,kcm:), TARGET, OPTIONAL, INTENT(IN) :: &
!
    c_big,        & ! effective drag coefficient of canopy elements
!                   ! larger than or equal to the turbulent length scale (1/m)
    c_sml           ! effective drag coefficient of canopy elements
                    ! smaller than the turbulent length scale            (1/m)

REAL (KIND=wp), DIMENSION(:,kcm-1:), TARGET, OPTIONAL, INTENT(IN) :: &
    r_air           ! log of air containing fraction of a gridbox inside
!                   ! the canopy                                          (1)

! Fields for surface values and soil/canopy model variables:
! ------------------------------------------------------------

REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(IN) :: &
!
    ps,           & ! surface pressure                              ( pa  )
    qv_s,         & ! specific water vapor content on the surface   (kg/kg)
    t_g             ! weighted surface temperature                  (  k  )

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
!
! Atmospheric model variables:
! ---------------------------------
!
     u,           & ! zonal wind speed       (at mass positions)    ( m/s )
     v,           & ! meridional wind speed  (at mass positions)    ( m/s )
     t,           & ! temperature                                   (  k  )
     qv,          & ! specific water vapor content                  (kg/kg)
     qc             ! specific cloud water content                  (kg/kg)

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(IN) :: &
!
     prs            ! atmospheric pressure                          ( pa  )

REAL (KIND=wp), DIMENSION(:,:), TARGET, OPTIONAL, INTENT(IN) :: &
!
     rho,         & ! total density of air                          (kg/m3)
     epr            ! exner pressure                                 (1)

REAL (KIND=wp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: &
!
     w              ! vertical wind speed (defined on half levels)  ( m/s )

TYPE (modvar),            OPTIONAL :: ptr(:) ! passive tracers
INTEGER,                  OPTIONAL :: ndtr   ! number of tracers to be diffused
INTEGER                            :: ndiff  ! number of 1-st order variables

REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(INOUT) :: &
!
! Diagnostic surface variable of the turbulence model:
! -----------------------------------------------------
!
     gz0,          & ! roughness length * g of the vertically not
                     ! resolved canopy                               (m2/s2)
!Achtung: Der g-Faktor ist ueberfluessig!

!    turbulent (transfer) velocity scales at the surface
     tvm,          & ! for momentum                                  ( m/s)
     tvh,          & ! for heat and moisture                         ( m/s)

     !Notice that 'tcm' and 'tch' are dispensable. The common use of the related
     !vecolities  'tvm' and 'tvh' makes live much easier!!               

!    turbulent transfer factors for laminar- and roughness-layer transfer
     tfm,          & ! of momentum                                     --
     tfh,          & ! of scalars                                      --
     tfv             ! of water vapor compared to heat                 --

REAL (KIND=wp), DIMENSION(:), TARGET, OPTIONAL, INTENT(INOUT) :: &

!    reference surface diffusion coefficient
!    (only if 'ltursrf' and "imode_trancnf.GE.2"):
     tkr             ! l*Ustar                                       (m2/s)

REAL (KIND=wp), DIMENSION(:), TARGET, OPTIONAL, INTENT(IN) :: &
     tkred_sfc       ! reduction factor for minimum diffusion coefficients near the surface

! Atmospheric variables of the turbulence model:
! ------------------------------------------------

REAL (KIND=wp), DIMENSION(nvec,ke1,ntim), TARGET, INTENT(INOUT) :: &

     tke             ! q:=SQRT(2*TKE); TKE='turbul. kin. energy'     ( m/s )
                     ! (defined on half levels)

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &

     tkvm,         & ! turbulent diffusion coefficient for momentum  (m2/s )
     tkvh            ! turbulent diffusion coefficient for heat      (m2/s )
                     ! (and other scalars)

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &

     rcld            ! standard deviation of the local oversaturation
                     ! (as input and output)
                     ! fractional cloud cover (in turbdiff)            --

REAL (KIND=wp), DIMENSION(:,:), OPTIONAL, TARGET, INTENT(IN) :: &

     hdef2,        & ! horizontal deformation square at half levels  ( 1/s2 )
     hdiv,         & ! horizontal divergence                   ,,    ( 1/s )

     dwdx,         & ! zonal      derivative of vertical wind  ,,    ( 1/s )
     dwdy            ! meridional derivative of vertical wind  ,,    ( 1/s )

REAL (KIND=wp), DIMENSION(:,:), OPTIONAL, TARGET, INTENT(INOUT) :: &
!
! Tendency fields for the prognostic variables:
! -----------------------------------------------
!
     u_tens,       & ! u-tendency                                    ( m/s2)
     v_tens,       & ! v-tendency                                    ( m/s2)
     t_tens,       & ! t-tendency                                    ( K/s )
     qv_tens,      & ! qv-tendency                                   ( 1/s )
     qc_tens         ! qc-tendency                                   ( 1/s )

REAL (KIND=wp), DIMENSION(:,:), TARGET, OPTIONAL, INTENT(INOUT) :: &
!
     tketens,      & ! diffusion tendency of q=SQRT(2*TKE)           ( m/s2)
     tketadv         ! advection tendency of q=SQRT(2*TKE)           ( m/s2)

REAL (KIND=wp), DIMENSION(:), TARGET, OPTIONAL, INTENT(IN) :: &
!
     trop_mask      ! mask-factor (1: within tropics; 0: within extra-tropics)
                    ! used for vertical smoothing of TKE forcing terms
REAL (KIND=wp), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: &
!
     qv_conv         ! qv-flux-convergence                            ( 1/s )

REAL (KIND=wp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: &
!
     ut_sso,       & ! u-tendency due to the SSO-Scheme              ( 1/s )
     vt_sso          ! v-tendency due to the SSO-Scheme              ( 1/s )

REAL (KIND=wp), DIMENSION(:,:), OPTIONAL, TARGET, INTENT(OUT) :: &
!
     edr,          & ! eddy dissipation rate of TKE (EDR)            (m2/s3)
     tket_sso,     & ! TKE-tendency due to SSO wake production       (m2/s3)
     tket_hshr,    & ! TKE-tendency due to separ. horiz. shear       (m2/s3)
!
     tkhm,         & ! horizontal diffusion coefficient for momentum ( m2/s )
     tkhh            ! horizontal diffusion coefficient for scalars  ( m2/s )

REAL (KIND=wp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: &
!
     tket_conv       ! TKE-tendency due to convective buoyancy       (m2/s3)

REAL (KIND=wp), DIMENSION(:), OPTIONAL, TARGET, INTENT(INOUT) :: &
!
     shfl_s,       & ! sensible heat flux at the surface             (W/m2)    (positive downward)
     lhfl_s,       & ! latent   heat flux at the surface             (W/m2)    (positive downward)
     qvfl_s,       & ! water vapor   flux at the surface             (kg/m2/s) (positive downward)
     umfl_s,       & ! u-momentum flux at the surface                (N/m2)    (positive downward)
     vmfl_s          ! v-momentum flux at the surface                (N/m2)    (positive downward)

INTEGER, INTENT(INOUT) :: ierrstat

CHARACTER (LEN=*), INTENT(INOUT) :: yroutine
CHARACTER (LEN=*), INTENT(INOUT) :: yerrormsg

!-------------------------------------------------------------------------------
!Local Parameters:
!-------------------------------------------------------------------------------

INTEGER :: &
    i, k,       & !horizontaler und vertikaler Laufindex
!US    nvor,       & !laufende Zeittstufe des bisherigen TKE-Feldes
    it_start,   & !Startindex der Iterationen
    kem,        & !ke oder ke1
    ntrac         !number of included passive tracers

LOGICAL ::  &
    lini,      & !initialization required
    ldovardif, & !berechne (teil-)implizite Vert.diff von Mod.var 1-ter Ordnung
    ldogrdcor, & !mache Gradientkorrektur bei Berechnung der vertikalen Diffusion
    lssintact    !trenne Skalen-Interaktionsterme vom mech. Forcing ab

REAL (KIND=wp) :: &
    fr_tke              ! z1/dt_tke

#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
REAL (KIND=wp), POINTER, CONTIGUOUS :: &
#else
REAL (KIND=wp), POINTER :: &
#endif
!
!    pointer for density, exner factor and eddy dissipation rate:
     rhoh(:,:), exner(:,:), ediss(:,:), &
!
!    pointer for densities of sh-flux and qv-flux:
     shflu_s(:), qvflu_s(:)

!-------------------------------------------------------------------------------

! Declarations:

!     Lokale logical Variablen:

      LOGICAL ldotkedif, & !berechne (teil-)implizite Vert.diff von TKE
              lcircterm, & !Zirkulationsterm muss berechnet werden
              lcircdiff, & !Zirkulationsterm wird zusammen mit TKE-Diffusion bestimmt
              linisetup, & !initiales setup bei impliziter Vertikaldiffusion
              lnewvtype, & !neuer Variablentyp muss vorbereitet werden
              lsflucond    !untere Flussrandbedingung

!     Lokale Integer-Hilfsvariablen:

!     Lokale Integer-Hilfsvariablen:

      INTEGER :: &
              ii,n,m,kk,    & !Indices fuer diverse Schleifen
              ku,k1,k2,     & !Schicht-Indices
              kgc,          & !oberer Schichtindex des Bereiches mit Gradientkorrektur
              it_durch,     & !Durchgangsindex der Iterationen
              nprim,        & !Erster  Index der diffundierenden Variablen
              nlast,        & !Letzter Index der diffundierenden Variablen
              ncorr,        & !Startindex der Variablen mit Gradientkorrektur
              igrdcon,      & !Index fuer Modus der Gradientberuecksichtigung
              itndcon,      & !Index fuer Modus der  Tendenzberuecksichtigung
              ivtype          !Index fuer Variablentyp

      INTEGER :: &
!
!             Eingrenzende Hoehenvieaus bei der Berechnung der
!             gemittelten Profile bei nicht-lokaler Gradient-Berechnung:
!
              lev(2)

!     Lokale real Variablen:

      REAL (KIND=wp) :: &
!
!          Hilfsvariablen:
!
           wert, val1, val2, & ! Platzhalter fuer beliebige Zwischenergebnisse
           fakt,             & !  ,,         ,,     ,,     Faktoren
!
!     Platzh. fuer thermodynamische Hilfsgreossen
!
           virt,      & !z1+(Rv/Rd-1)*qv-qc
           flw_h2o_g, & !                 rc/(1+(lh_v/cp_d)*d_qsat/d_T)
           flw_tet_l    !exner*d_qsat/d_T*rc/(1+(lh_v/cp_d)*d_qsat/d_T)

      REAL (KIND=wp) :: &
!
           fr_var       !1/dt_var

      REAL (KIND=wp) :: &
!
!     Platzh. fuer therm. und mech. Antrieb der Turbulenz in (1/s)**2
!     (fh2,fm2):
!
           fh2,fm2, &
!
!     Platzh. fuer horiz. Geschw.-Komponenten und bel. Geschw.:
!
           vel1,vel2,velo, &
!
!     Platzh. fuer die Hoehe ueber Grund, Hoehendifferenzen, obere und
!     untere Hoehe, turbulente Laengaenskalen, Kohaerenzlaenge,
!     Dicke der laminaren Grenzschicht,sowie eine Laenge allgemein:
!
           h,hu,l_turb,lh,lm,kohae_len,len, &
!
           edh, & ! Kehrwert von Schichtdicken
!
!     Zwischenspeicher fuer
!
           thermik, & !(negative) Auftriebsproduktion an TKE
           phasdif, & !Temperaturtendenz durch Phasendiffusion
!
!
!Tuning<
!          x1,x2,x3,
           x4,xri(nvec,ke)
!>Tuning

! Local arrays:

INTEGER :: &
    ivtp(ndiff) ! index of variable type

! Time increment and inverse time increment of ordinary prognostic variables:
REAL (KIND=wp), TARGET :: &
     tinc(ndiff)      ,& !
     tinv(ndiff)      ,& !

     hig(2)              ! obere und untere Referenzhoehe bei der Bildung
                         ! nicht-lokaler Gradienten

#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      REAL (KIND=wp), POINTER, CONTIGUOUS :: &
#else
      REAL (KIND=wp), POINTER :: &
#endif
!
!          Pointer fuer Tendenzfelder:
!
           utens(:,:), vtens(:,:), &
           ttens(:,:), qvtens(:,:), qctens(:,:), &
!
!          Pointer fuer Felder der Bestandes-Parameter:
!
           cbig(:,:), csml(:,:), rair(:,:)

! Note:
! The following buffers wouldn't be necessary, if the related pointers above
! were allowed to be allocated at run time:

#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      REAL (KIND=wp), DIMENSION(:,:), POINTER, CONTIGUOUS :: &
#else
      REAL (KIND=wp), DIMENSION(:,:), POINTER :: &
#endif
!
            cur_prof, upd_prof, sav_prof, &
!++++
            expl_mom, impl_mom, invs_mom, &
            diff_dep, eff_flux
!++++

      LOGICAL ::  &
!
           ltend(ndiff), &  !calculation of tendencies required
           lsfli(ndiff), &  !surface value input is a flux density instead of a concentration
           leff_flux,    &  !calculation of effective flux density required
           lcalc_frcsmot    !local control variable if smoothing of TKE forcing terms needs to be calculated

      TYPE (modvar) :: dvar(ndiff)  !model variables to be diffused

      TYPE (turvar) :: vtyp(ntyp)   !variable types (momentum and scalars)

      TYPE (varprf) :: pvar(naux+1) !vertical variable profiles

!---- End of header ------------------------------------------------------------

integer :: i2d, j2d

!===============================================================================
! these are still settings from organize_turbdiff

!All variables and their tendencies are defined at horizontal mass positions.

! for this call: itnd=0

 ldogrdcor=(lexpcor .AND. lturatm)             !gradient correction has to be done
 ldovardif=(lum_dif .OR. lvm_dif .OR. lscadif) !some variable has to be diffused

 lssintact=((ltkesso.OR.ltkeshs.OR.ltkecon) .AND. imode_adshear.EQ.1)

 IF (PRESENT(ptr)) THEN !passive tracers are present
    !number of tracers
    IF (PRESENT(ndtr)) THEN
       ntrac = ndtr
    ELSE
       ntrac = UBOUND(ptr,1)
    END IF
 ELSE
    ntrac=0
 END IF

 fr_tke=z1/dt_tke

 kem=ke

!DIR$ IVDEP
 DO i=ivstart, ivend
!Achtung: Korrektur durch Faktor 1/2 (wirkt bei sehr kleinen horiz. Gitterzellen
    l_scal(i)=MIN( z1d2*l_hori(i), tur_len )
!__________________________________________________________________________
!test: frm ohne fc_min-Beschraenkung: Bewirkt Unterschiede!
!     fc_min(i)=(vel_min/MAX( l_hori(i), tur_len ))**2
      fc_min(i)=z0
!__________________________________________________________________________
 END DO

 IF (iini.GT.0) THEN !an initialization run
    lini=.TRUE.
    IF (iini.EQ.1) THEN !separate initialization before the time loop
       it_start=1 !only 'it_end' iterations for initialization
                  !and an additional one at the first time loop
    ELSE !initialization within the first time step
       it_start=0 !"it_end+1" iterations for initializatzion
    END IF

!   Bestimmung der initialen Werte fuer die laminaren Reduktions-
!   koeffizienten und die Standardabw. des Saettigungsdef.:
!   Initializing some special variables:

#ifdef __ICON__
    DO k=1, ke1
!DIR$ IVDEP
       DO i=ivstart, ivend
          rcld(i,k)=z0 !no standard-deviat. of local over-saturation
       END DO
    END DO
!DIR$ IVDEP
    DO i=ivstart, ivend
       tfh(i)=z1 !no roughness- and laminar-layer-resistance for scalars
       tfm(i)=z1 !no roughness- and laminar-layer-resistance for momentum
    END DO
    !Notice that the above initalizations will stay, if another turbulence-
    !or transfer-scheme is used!
#endif
 ELSE !not an initialization run
    lini=.FALSE.
    it_start=it_end !only one iteration
!test:
!it_start=1
!test
 END IF

 !Note:
 !A call with "iini=2" (at the first time step) provides the same result as a
 !  call with "iini=1" (before the   time loop) followed by a second
 !  call with "iini=0" (at the first time step).

!US now pass nvor
!US nvor=nprv !Eingangsbelegung von 'nvor' (wird bei Iterationen auf 'ntur' gesetzt)

 !Note:
 !It is also possible to use only one time level for TKE ("ntim=1" and thus "nprv=1=ntur").

 IF (PRESENT(rho)) THEN
    rhoh => rho
 ELSE
    rhoh => rho_tar
  ! ALLOCATE ( rhoh(ie,ke),   STAT=ilocstat ); istat = istat + ilocstat
 END IF

 IF (PRESENT(epr)) THEN
    exner => epr
 ELSE
    exner => exner_tar
  ! ALLOCATE ( exner(ie,ke),  STAT=ilocstat ); istat = istat + ilocstat
 END IF

 IF (PRESENT(edr)) THEN
    ediss => edr
 ELSE
    ediss => diss_tar
  ! ALLOCATE ( ediss(ie,ke1),  STAT=ilocstat ); istat = istat + ilocstat 
 END IF
 IF (PRESENT(shfl_s)) THEN
    shflu_s => shfl_s
 ELSE
    shflu_s => shflu_s_tar
  ! ALLOCATE ( shflu_s(ie),    STAT=ilocstat ); istat = istat + ilocstat 
 END IF
 IF (PRESENT(qvfl_s)) THEN
    qvflu_s => qvfl_s
 ELSE
    qvflu_s => qvflu_s_tar
  ! ALLOCATE ( qvflu_s(ie),    STAT=ilocstat ); istat = istat + ilocstat 
 END IF

!Achtung: Korrektur: Konsistente Behandlung der Null-Fluss-Randbedingung
!als moeglicher 'default' (etwa fuer qc)
 IF (ilow_def_cond.EQ.2) THEN !zero surface value of liquid water
!DIR$ IVDEP
    DO i=ivstart, ivend
       liqs(i,ke1)=z0
    END DO
 ELSE !constant liquid water within the transfer-layer
!DIR$ IVDEP
    DO i=ivstart, ivend
       liqs(i,ke1)=qc(i,ke)
    END DO
 END IF

!-------------------------------------------------------------------------------

!US IF (lturatm .OR. ldovardif) THEN
!US    IF (.NOT.lerror) CALL turbdiff
!US END IF

!===============================================================================
! here turbdiff really starts

!     nur in LM-Umgebung:

      istat=0; ilocstat=0; ierrstat=0
      yerrormsg = ''; yroutine='turbdiff'; lerror=.FALSE.

!     Fuer die Turb.par. benutzter Variablensatz auf Hauptflaechen:
!     Bei k=ke1 stehen die unteren Randwerte der Prandtlschicht
!     (skin-layer-Werte)

!     Der letzte Index bezeichnet die physik. Bedeutung der Variablen
!     und hat einen der Werte u_m,v_m,tet_l,h2o_g,liq;
!                        bzw. tem,vap,liq.
!     Der letzte Index bezeichnet die physik. Bedeutung der Variablen
!     Bei k=ke1 stehen die unteren Randwerte der Prandtlschicht
!     (skin-layer-Werte)
!     Am Ende des Programms wird vari mit den (Co-Varianzen) der
!     Geschwindigkeitskomponenten ueberschrieben, die fuer die
!     Berechung der Horizontaldiffusion benoetigt werden.
!     vari() enthaelt spaeter auch die nmvar (nichtlokalen) vertikalen
!     Gradienten und auch die durch Wirkung der subskaligen Kondensation
!     veraenderten (effektiven) vertikalen Gradienten.
!     Zum Schluss enthaelt vari() fuer die turbulente Horizontaldiff.
!     benoetigte Komponenten des turbulenten Spannungstensors.

!print *,"in turbdiff c_diff=",c_diff," tkhmin=",tkhmin

!########################################################################

      fr_var=z1/dt_var

      IF (PRESENT(c_big)) THEN
         cbig => c_big
      ELSE
         cbig => cbig_tar
       ! ALLOCATE ( cbig(ie,kcm:ke1), STAT=ilocstat ); istat = istat + ilocstat
      END IF
      IF (PRESENT(c_sml)) THEN
         csml => c_sml
      ELSE
         csml => csml_tar
       ! ALLOCATE ( csml(ie,kcm:ke1), STAT=ilocstat ); istat = istat + ilocstat
      END IF
      IF (PRESENT(r_air)) THEN
         rair => r_air
      ELSE
         rair => rair_tar
       ! ALLOCATE ( rair(ie,kcm-1:ke1), STAT=ilocstat ); istat = istat + ilocstat
      END IF

      IF (istat /= 0) THEN
         ierrstat = 1004
         yerrormsg= &
         'ERROR *** Allocation of space for meteofields failed ***'
         lerror=.TRUE.; RETURN
      ENDIF

      ltend(u_m)=PRESENT(u_tens)
      IF (ltend(u_m)) THEN !calculation of tendencies required
         utens => u_tens !'utens' points to the tendency
      ELSE                 !update of ordinary prognostic variables required
         utens => u      !'utens' points to the prognostic variables
      END IF
      ltend(v_m)=PRESENT(v_tens)
      IF (ltend(v_m)) THEN
         vtens => v_tens
      ELSE
         vtens => v
      END IF
      ltend(tem)=PRESENT(t_tens)
      IF (ltend(tem)) THEN
         ttens => t_tens
      ELSE
         ttens => t
      END IF
      ltend(vap)=PRESENT(qv_tens)
      IF (ltend(vap)) THEN
         qvtens => qv_tens
      ELSE
         qvtens => qv
      END IF
      ltend(liq)=PRESENT(qc_tens)
      IF (ltend(liq)) THEN
         qctens => qc_tens
      ELSE
         qctens => qc
      END IF

      ! check if vertical smoothing of TKE forcing terms is needed
      IF (frcsmot > z0) THEN
        IF (.NOT. PRESENT(trop_mask)) THEN
          lcalc_frcsmot = .TRUE.
        ELSE IF (ANY(trop_mask(ivstart:ivend) > z0)) THEN
          lcalc_frcsmot = .TRUE.
        ELSE
          lcalc_frcsmot = .FALSE.
        ENDIF
      ELSE
        lcalc_frcsmot = .FALSE.
      ENDIF

      lsfli(:)=.FALSE. !surface values are concentrations by default

      dvar(u_m)%av  => u  ; dvar(u_m)%at => utens  ; dvar(u_m)%sv => NULL()
      dvar(v_m)%av  => v  ; dvar(v_m)%at => vtens  ; dvar(v_m)%sv => NULL()

!Note: Use                                           var(u_m)%sv => u(:,ke)
!      and                                           var(v_m)%sv => v(:,ke)
!      in order to force a "free-slip condition"!

      dvar(tem)%av  => t  ; dvar(tem)%at => ttens  ; dvar(tem)%sv => t_g
      dvar(vap)%av  => qv ; dvar(vap)%at => qvtens ; dvar(vap)%sv => qv_s
      dvar(liq)%av  => qc ; dvar(liq)%at => qctens ; dvar(liq)%sv => NULL()

!SCLM --------------------------------------------------------------------------------
#ifdef SCLM
      IF (lsclm) THEN
         IF (SHF%mod(0)%vst.GT.i_cal .AND. SHF%mod(0)%ist.EQ.i_mod) THEN
            !measured SHF has to be used for forcing:
            lsfli(tem)=.TRUE.
         END IF
         IF (LHF%mod(0)%vst.GT.i_cal .AND. LHF%mod(0)%ist.EQ.i_mod) THEN
            !measured LHF has to be used for forcing:
            lsfli(vap)=.TRUE.
         END IF
      END IF
      !Note: the measured SHF and LHF have to be present by shfl_s and qvfl_s!
#endif
!SCLM --------------------------------------------------------------------------------

      IF (lsfluse) THEN !use explicit heat flux densities at the surface
         lsfli(tem)=.TRUE.; lsfli(vap)=.TRUE.
      END IF

      IF ((lsfli(tem) .AND. .NOT.PRESENT(shfl_s)) .OR. &
          (lsfli(vap) .AND. .NOT.PRESENT(qvfl_s))) THEN
         ierrstat = 1004; lerror=.TRUE.
         yerrormsg='ERROR *** forcing with not present surface heat flux densities  ***'
         RETURN
      ENDIF

      IF (lsfli(tem)) dvar(tem)%sv => shfl_s
      IF (lsfli(vap)) dvar(vap)%sv => qvfl_s

      IF (PRESENT(ptr) .AND. ntrac .GE. 1) THEN !passive tracers are present
         DO m=1, ntrac
            n=liq+m
            dvar(n)%av => ptr(m)%av
            ltend(n)=ASSOCIATED(ptr(m)%at)
            IF (ltend(n)) THEN
               dvar(n)%at => ptr(m)%at
            ELSE
               dvar(n)%at => ptr(m)%av
            END IF
            IF (ASSOCIATED(ptr(m)%sv)) THEN
               dvar(n)%sv => ptr(m)%sv; lsfli(n)=ptr(m)%fc
            ELSE
               dvar(n)%sv => NULL()   ; lsfli(n)=.FALSE.
            END IF
         END DO
      END IF

      vtyp(mom)%tkv => tkvm ; vtyp(mom)%tsv => tvm
      vtyp(sca)%tkv => tkvh ; vtyp(sca)%tsv => tvh

      !Note:
      !If a tendency field of an ordinary prognostic variable is not present,
      !the related time step increment due to turbulent diffusion will be
      !added to the prognostic variable directly.

      DO n=1,ndiff
         IF (ltend(n)) THEN !calculation of tendencies required
            tinc(n)=z1        !no time increment multiplication for tendencies
            tinv(n)=fr_var    !division by time increment for variable increments
         ELSE               !update of prognostic variables required
            tinc(n)=dt_var    !time increment multiplication for tendencies
            tinv(n)=z1        !no division by time increment for variable increments
         END IF
         IF (n.LE.nvel) THEN
            ivtp(n)=mom
         ELSE
            ivtp(n)=sca
         END IF
!print *,"n=",n," tinc=",tinc(n)," tinv=",tinv(n)
      END DO

      IF (l3dturb .AND..NOT. (PRESENT(tkhm) .AND. PRESENT(tkhh))) THEN
         ierrstat = 1004; lerror=.TRUE.
         yerrormsg='ERROR *** 3D-diffusion with not present horiz. diff.coeffs. ***'
      END IF

!########################################################################

!-----------------------
   IF (lturatm) THEN
!-----------------------

! 0)  Berechnung der Erhaltungsvariablen (auf 'vari') samt des Bedeckungsgrades
!     und thermodynamischer Hilfgroessen, sowie der turbulenten Laengenskalen:

!-----------------------
!Achtung: Bei T-Gleichung in cv-Form
!         gesamte Thermodynamik ueberpruefen auch gegenueber satad

      lcircterm=(pat_len.GT.z0)
      ldotkedif=(c_diff .GT.z0)
      lcircdiff=(lcircterm .AND. imode_calcirc.EQ.2)

!     Berechnung abgeleiteter Parameter:

!US   CALL turb_param

!     Thermodynamische Hilfsvariablen auf Hauptflaechen:

      CALL adjust_satur_equil ( khi=1, ktp=1, &
!
           i_st=ivstart, i_en=ivend, k_st=1, k_en=ke, i1dim=nvec,   &
!
           lcalrho=.NOT.PRESENT(rho), lcalepr=.NOT.PRESENT(epr),    &
           lcaltdv=.TRUE., lpotinp=.FALSE., ladjout=.FALSE.,        &
!
           icldmod=icldm_turb,                                      &
!
           zrcpv=tur_rcpv, zrcpl=tur_rcpl,                          &
!
           prs=prs, t=t,     qv=qv,    qc=qc,                       &
!
           psf=ps,                                                  &
!
           exner=exner, rcld=rcld, dens=rhoh,                       &
!
           r_cpd=a(:,:,2), qst_t=a(:,:,3), g_tet=a(:,:,4),          &
                                           g_h2o=a(:,:,5),          &
!
           tet_l=vari(:,:,tet_l), q_h2o=vari(:,:,h2o_g),            &
                                  q_liq=vari(:,:,liq) )

!     Thermodynamische Hilfsvariablen auf Unterrand der Prandtl-Schicht:

!DIR$ IVDEP
      DO i=ivstart, ivend
         prss(i,ke1)=ps(i)
         tmps(i,ke1)=t_g(i)
         vaps(i,ke1)=qv_s(i)
      END DO

      CALL adjust_satur_equil ( khi=ke1, ktp=ke, &
!
           i_st=ivstart, i_en=ivend, k_st=ke1, k_en=ke1, i1dim=nvec,&
!
           lcalrho=.TRUE., lcalepr=.TRUE.,                          &
           lcaltdv=.TRUE., lpotinp=.FALSE., ladjout=.FALSE.,        &
!
           icldmod=icldm_turb,                                      &
!
           zrcpv=tur_rcpv, zrcpl=tur_rcpl,                          &
!
!Achtung: Korrektur: Konsistente Behandlung der unteren Null-Fluss-Randbedingung fuer qc
!         und an COSMO-Version angepasste Interpretation von "icldmod=-1":
           prs=prss, t=tmps, qv=vaps, qc=liqs,                      &
!
           fip=tfh,                                                 &
!
           exner=a(:,ke1:ke1,1), rcld=rcld(:,ke1:ke1),              &
           dens=rhon(:,ke1:ke1),                                    &
!
           r_cpd=a(:,ke1:ke1,2), qst_t=a(:,ke1:ke1,3),              &
           g_tet=a(:,ke1:ke1,4), g_h2o=a(:,ke1:ke1,5),              &
!
           tet_l=vari(:,ke:ke1,tet_l), q_h2o=vari(:,ke:ke1,h2o_g),  &
                                       q_liq=vari(:,ke:ke1,liq) )

!     Beachte:
!     'vari(:,ke1,tet_l)' und 'vari(:,ke1,h2o_g) sind jetzt die Erhaltungsvariablen
!      am Unterrand der Prandtl-Schicht (zero-level). Die Werte im Niveau 'ke' wurden dabei
!      zur Interpolation benutzt.
!     'a(:,ke1,1)' enthaelt den Exner-Faktor im zero-level. Das Feld 'a(:,:,1) wird im Folgenden
!      mit dem Exner-Faktor auf NF belegt.

!     Kommentar: 
!     Der 2-te Aufruf von 'adjust_satur_equil' stellt die unteren Randwerte der thermodyn. Variablen
!      zur Verfuegung. Dies koennte in den 1-ten Aufruf integriert werden, wenn alle thermodyn.
!      Modell-Variablen bis "k=ke1" allociert waeren. Dies wuere Vieles vereinfachen!

      IF (imode_trancnf.EQ.1) THEN !old version of zero-level-values requested
         !Transformation of Tet_l at zero-level into the value following from the old
         !treatment of interpolation in terms of T_l (rather than Tet_l):
         DO i=ivstart, ivend
            vari(i,ke1,tet_l) = vari(i,ke1,tet_l)  &
                            + ( (exner(i,ke)-a(i,ke1,1))*vari(i,ke,tet_l)*(z1-tfh(i)) ) &
                              / a(i,ke1,1)
         END DO
      END IF

!     Berechnung der horizontalen Windgeschwindigkeiten
!     im Massenzentrum der Gitterbox:

      DO k=1,ke
!DIR$ IVDEP
         DO i=ivstart, ivend
            vari(i,k,u_m)=u(i,k)
            vari(i,k,v_m)=v(i,k)
         END DO
      END DO
!DIR$ IVDEP
      DO i=ivstart, ivend
         vari(i,ke1,u_m)=vari(i,ke,u_m)*(z1-tfm(i))
         vari(i,ke1,v_m)=vari(i,ke,v_m)*(z1-tfm(i))
      END DO

!     Berechnung der Schichtdicken und der Dichte auf Nebenflaechen:

      DO k=1,ke
!DIR$ IVDEP
         DO i=ivstart, ivend
            dicke(i,k)=hhl(i,k)-hhl(i,k+1)
         END DO
      END DO

!     Interpolation der thermodyn. Hilfsgroessen im Feld a(),
!     der Wolkendichte und der Luftdichte auf Nebenflaechen:

      CALL bound_level_interp( ivstart, ivend, 2,ke, &
!___________________________________________________________________________
!test: mass weighted interpolation
!                              nvars=1, pvar=(/varprf(rhon,rhoh)/), depth=dicke)
!Achtung: Macht minimale Unterschiede
                               nvars=1, pvar=(/varprf(rhon,rhoh)/), depth=dp0)
!___________________________________________________________________________

      pvar(1)%bl => rcld     ; pvar(1)%ml => rcld  !NF-Werte wieder nach 'rcld'
      pvar(2)%bl => a(:,:,1) ; pvar(2)%ml => exner !NF-Werte nach 'a(:,:,1)', weil
                                                   !'exner' auf HF noch benoetigt wird.
      DO n=2,naux
         pvar(1+n)%bl => a(:,:,n) ; pvar(1+n)%ml => a(:,:,n)
      END DO
      !Note: 
      !Internal order of level looping in 'bound_level_interp' allows to store the 
      !'bl'-values (output) at the same place as the 'ml'-values (input).
      
      CALL bound_level_interp( ivstart, ivend, 2,ke, &
                               nvars=naux+1, pvar=pvar, &
                               depth=dp0, rpdep=hlp)

!print *,"nach interpol"

      !Spezifische effektive Dicke der Prandtlschicht:

!DIR$ IVDEP
      DO i=ivstart, ivend

!Achtung: < zum Vergleich mit alter Variante:
! velo=MAX( vel_min, SQRT(vari(i,ke,u_m)**2+vari(i,ke,v_m)**2) )
! tvm(i)=tcm(i)*velo
! tvh(i)=tch(i)*velo
!> zum Vergleich: in ICON werden zwischen turbtran und turbdiff 'u' und 'v' incrementiert!
!Achtung: Modifikation tcm -> tvm; tch -> tvh: macht Unterschiede
         dzsm(i)=tkvm(i,ke1)*tfm(i)/tvm(i)
         dzsh(i)=tkvh(i,ke1)*tfh(i)/tvh(i)
      END DO

      !Beachte: Auch wenn 'turbtran' nicht als Transfer-Schema benutzt wird, muessen 'tkvm', tkvh'
      !         und 'tke' und (falls es PRESENT ist) auch 'edr' fuer "k=ke1" belegt sein!

!     Berechnung der turbulenten Laengenscalen:

!DIR$ IVDEP
      DO i=ivstart, ivend
         len_scale(i,ke1)=gz0(i)/grav
      END DO
      DO k=ke,kcm,-1 !Innerhalb des Bestandesmodells
!DIR$ IVDEP
         DO i=ivstart, ivend
            IF (cbig(i,k).GT.z0) THEN
!              Die turbulente Laengenskala wird durch die Laengenskala 
!              der lufterfuellten Zwischenraeume limitiert:

               l_turb=z1/(cbig(i,k)*sqrt(z1/exp(rair(i,k))-z1))
               len_scale(i,k)=MIN( dicke(i,k)+len_scale(i,k+1), l_turb )
            ELSE
               len_scale(i,k)=dicke(i,k)+len_scale(i,k+1)
            END IF
         END DO
      END DO
      DO k=kcm-1,1,-1
!DIR$ IVDEP
         DO i=ivstart, ivend
            len_scale(i,k)=dicke(i,k)+len_scale(i,k+1)
         END DO
      END DO

!     Uebergang von der maximalen turbulenten Laengenskala zur
!     effektiven turbulenten Laengenskala:

      DO k=ke1,1,-1
!DIR$ IVDEP
         DO i=ivstart, ivend
            lay(i)=l_scal(i)
         END DO
!DIR$ IVDEP
         DO i=ivstart, ivend
            len_scale(i,k)=akt*MAX( len_min, &
                                  ! len_scale(i,k)/(z1+len_scale(i,k)/lay(i)) )
                                    lay(i)*len_scale(i,k)/(lay(i)+len_scale(i,k)) )
         END DO
      END DO

!     Initialisierung der Felder fuer tke,tkvh,tkvm:

      IF (lini) THEN  !nur beim allerersten Durchgang

!        Erste Schaetzwerte aus vereinfachtem TKE-Gleichgewicht:

         DO k=2,kem
            DO i=ivstart, ivend

!              der Einfachheit halber nur lokale Berechnung der
!              vertikalen Gradienten:

               len=len_scale(i,k)
               edh=z2/(hhl(i,k+1)-hhl(i,k-1))

               grad(u_m  )=(vari(i,k,u_m  )-vari(i,k-1,u_m  ))*edh
               grad(v_m  )=(vari(i,k,v_m  )-vari(i,k-1,v_m  ))*edh
               grad(tet_l)=(vari(i,k,tet_l)-vari(i,k-1,tet_l))*edh
               grad(h2o_g)=(vari(i,k,h2o_g)-vari(i,k-1,h2o_g))*edh

               fh2=a(i,k,4)*grad(tet_l)+a(i,k,5)*grad(h2o_g)
               fm2=MAX( grad(u_m)**2+grad(v_m)**2, fc_min(i) )

               ! Vereinfachte Loesung mit Rf=Ri:
               IF (fh2.GE.(z1-rim)*fm2) THEN
                  ! Die krit. Ri-Zahl wird ueberschritten und lm, sowie lh
                  ! werden durch lm bei der krit. Ri-Zahl angenaehert:
                  fakt=z1/rim-z1
                  lm=len*(sm_0-(a_6+a_3)*fakt)
                  lh=lm
               ELSE
                  fakt=fh2/(fm2-fh2)
                  lm=len*(sm_0-(a_6+a_3)*fakt)
                  lh=len*(sh_0-a_5*fakt)
               END IF

               val1=lm*fm2
               val2=lh*fh2
               hlp(i,k)=MAX( val1-val2, rim*val1 )

               IF (ltkeinp) THEN
!Achtung: Bug: Auch in ICON falsch !
                  tke(i,k,1)=tke(i,k,ntur)
                ! tke(i,ke,1)=tke(i,ke,ntur)
               ELSE
                  tke(i,k,1)=MAX( vel_min, SQRT(d_m*len*hlp(i,k)) ) !Initialwert fuer SQRT(2TKE)
               END IF

               val1=MAX ( con_m, tkmmin ); tkvm(i,k)=lm*tke(i,k,1)
               val2=MAX ( con_h, tkhmin ); tkvh(i,k)=lh*tke(i,k,1)
               IF (imode_tkemini.EQ.2) THEN
                  tke(i,k,1)=tke(i,k,1)*MAX( z1, val1/tkvm(i,k), & !adapted tke
                                                 val2/tkvh(i,k) )
               END IF
!Achtung: Bislang fehlte die Beschraenkung: Macht Unterschiede geg. ICON
               tkvm(i,k)=MAX( val1, tkvm(i,k) ) !corrected tkv
               tkvh(i,k)=MAX( val2, tkvh(i,k) )

!              Am Anfang konnte noch keine Diffusion von q=SQRT(2*TKE) berechnet werden:

               tketens(i,k)=z0

            END DO

         END DO

!DIR$ IVDEP
         DO i=ivstart, ivend
            tke(i,1,1)=tke(i,2,1)
         END DO

         DO n=2,ntim
            DO k=1,kem
!DIR$ IVDEP
               DO i=ivstart, ivend
                  tke(i,k,n)=tke(i,k,1)
               END DO
            END DO
         END DO

      END IF

! 1)  Berechnung der benoetigten vertikalen Gradienten und
!     Abspeichern auf 'vari':

!     Am unteren Modellrand:

!DIR$ IVDEP
      DO i=ivstart, ivend
         lays(i,mom)=z1/dzsm(i)
         lays(i,sca)=z1/dzsh(i)
      END DO
      DO n=1,nmvar
!DIR$ IVDEP
         DO i=ivstart, ivend
            vari(i,ke1,n)=(vari(i,ke,n)-vari(i,ke1,n))*lays(i,ivtp(n))
         END DO
      END DO

!     An den darueberliegenden Nebenflaechen:

      IF (lnonloc) THEN

!        Berechnung nicht-lokaler Gradienten:

         DO n=1,nmvar

!           Berechnung vertikalen Integralfunktionen in hlp():

            DO i=ivstart, ivend
               hlp(i,ke1)=z0
            END DO

            DO k=ke,2,-1
!DIR$ IVDEP
               DO i=ivstart, ivend
                  hlp(i,k)=hlp(i,k+1)+vari(i,k,n)*(hhl(i,k)-hhl(i,k+1))
               END DO
            END DO

            k1=1
            k2=2
!DIR$ IVDEP
            DO i=ivstart, ivend
               lays(i,k1)=hhl(i,1)-hhl(i,2)
            END DO

            DO k=2,ke

!              Berechnung der nicht-lokalen Gradienten als mittlere
!              Differenzenquotienten ueber die stabilitaetsabhaengige
!              Laengenskala in tkvh() bzw tkvm():

!              Speichern der stab.abh. Laengenskala unter lay():

               IF (n.LE.nvel) THEN
!DIR$ IVDEP
                  DO i=ivstart, ivend
                     lay(i)=tkvm(i,k)/tke(i,k,nvor)
                  END DO
               ELSE
!DIR$ IVDEP
                  DO i=ivstart, ivend
                     lay(i)=tkvh(i,k)/tke(i,k,nvor)
                  END DO
               END IF

!              Bestimmung der nicht-lokalen Gradienten und
!              Zwischenspeichern derselben auf dem Feld dicke():

!DIR$ IVDEP
               DO i=ivstart, ivend

                  lays(i,k2)=hhl(i,k)-hhl(i,k+1)

                  IF (lay(i).LE. &
                     z1d2*MIN( lays(i,k1), lays(i,k2)) ) THEN

!                    Die vertikalen Diffusionswege schneiden weder
!                    eine untere noch eine obere Hauptflaeche. Es
!                    koennen somit die lokalen Gradienten genommen
!                    werden. Bei sehr kleinen Diffusionslaengen, muss
!                    aus num. Gruenden sogar die lokale Berechnung
!                    gewaehlt werden:

                     dicke(i,k)=z2*(vari(i,k-1,n)-vari(i,k,n)) &
                                    /(hhl(i,k-1)-hhl(i,k+1))
                  ELSE

!                    Berechn. der benoetigten Referenzhoehen und -level:

                     h=hhl(i,k)
                     hu=MAX( h-lay(i), hhl(i,ke1) )
                     hig(1)=hu+lay(i)
                     hig(2)=h+lay(i)

                     kk=k
111                  IF (hhl(i,kk).GT.hu) THEN
                        kk=kk+1
                        GOTO 111
                     END IF
                     ku=kk
                     DO ii=1,2
112                     IF (kk.GT.1) THEN
                           IF (hhl(i,kk).LE.hig(ii)) THEN
                              kk=kk-1
                              GOTO 112
                           END IF
                        END IF
                        lev(ii)=kk+1
                     END DO

!                    Berechnung der gemittelten Differenzenquotienten
!                    als Ausdruck fuer die nicht-lokalen Gradienten:

                     wert=hlp(i,ku)-hlp(i,k) &
                         +hlp(i,lev(2))-hlp(i,lev(1)) &
                         +vari(i,ku-1,n)*(hu-hhl(i,ku)) &
                         -vari(i,lev(1)-1,n)*(hig(1)-hhl(i,lev(1)))&
                         +vari(i,lev(2)-1,n)*(hig(2)-hhl(i,lev(2)))

                     dicke(i,k)=wert/(lay(i)*(h-hu))
                  END IF
               END DO

               kk=k1
               k1=k2
               k2=kk

            END DO

!           Sichern der nicht-lokalen Gradienten im Feld vari():

            DO k=2,ke
!DIR$ IVDEP
               DO i=ivstart, ivend
                  vari(i,k,n)=dicke(i,k)
               END DO
            END DO

         END DO

!        Belegung von dicke() mit den Schichtdicken*rho/dt_tke
!        bzgl. Nebenflaechen:

         DO k=2,ke
!DIR$ IVDEP
            DO i=ivstart, ivend
               dicke(i,k)=rhon(i,k)*z1d2*(hhl(i,k-1)-hhl(i,k+1))*fr_tke
            END DO
         END DO

      ELSE    ! lnonloc

!        Berechnung lokaler Gradienten:

         DO k=ke,2,-1
!DIR$ IVDEP
            DO i=ivstart, ivend
               len=(hhl(i,k-1)-hhl(i,k+1))*z1d2
               hlp(i,k)=z1/len
               dicke(i,k)=rhon(i,k)*len*fr_tke
            END DO
         END DO

         DO n=1,nmvar
            DO k=ke,2,-1
!DIR$ IVDEP
               DO i=ivstart, ivend
                  vari(i,k,n)=(vari(i,k-1,n)-vari(i,k,n))*hlp(i,k)
               END DO
            END DO
         END DO

      END IF    ! lnonloc
!print *,"nach gradient"

! 2)  Berechnung der verallgemeinerten Antriebsfunktionen einschliesslich der
!     Korrekturen innerhalb der Rauhigkeitsschicht (samt der Windtendenz durch Formreibung)
!     und der Scherung durch nicht-turbulente subskalige Stroemungen:

!-------------------------
!     Thermal forcing:
!-------------------------

!Achtung:
      DO k=2,ke1 !'frh'(ke1) wird fuer Zirkulationsterm und Temperaturkorrektur benoetigt
!DIR$ IVDEP
         DO i=ivstart, ivend
            frh(i,k)=a(i,k,4)*vari(i,k,tet_l) &
                    +a(i,k,5)*vari(i,k,h2o_g)
         END DO
      END DO
      !Note: a(:,:,4) and a(:,:,5) are free now.

!-------------------------
!     Total mechanical forcing:  
!-------------------------

     !hdef2 = (d1v2+d2v1)^2 + (d1v1-d2v2)^2 !horizontal deformation square     (at half levels)
     !hdiv  = (d1v1+d2v2)                   !horizontal wind-divergence            ,,

     !dwdx !zonal      derivation of vertical wind                                 ,,
     !dwdy !meridional derivation of vertical wind                                 ,,

     !vel_div=hdiv+dzdw=0 !Incomressibility

     !itype_sher = 0 : only single column vertical shear
     !             1 : previous and additional 3D horiz. shear correction
     !             2 : previous and additional 3D vertc. shear correction

     !ltkeshs: consider separated non-turbulent horizontal shear mode for TKE forcing
     !loutshs: consider separated non-turbulent horizontal shear mode for output

!     Mechanical forcing by vertical shear:

      IF (itype_sher.EQ.2 .AND. (PRESENT(dwdx) .AND. PRESENT(dwdy) .AND. PRESENT(hdiv))) THEN
         !Include 3D-shear correction by the vertical wind (employing incomressibility):
         DO k=2,kem
!DIR$ IVDEP
            DO i=ivstart, ivend
               frm(i,k)=MAX( (vari(i,k,u_m)+dwdx(i,k))**2+(vari(i,k,v_m)+dwdy(i,k))**2 &
                                                         +z3*hdiv(i,k)**2, fc_min(i) )
            END DO
         END DO
      ELSE   
         !Load pure single column shear:
         DO k=2,kem 
!DIR$ IVDEP
            DO i=ivstart, ivend
               frm(i,k)=MAX( vari(i,k,u_m)**2+vari(i,k,v_m)**2, fc_min(i))
            END DO
         END DO
      END IF    

!     Mechanical forcing by horizontal shear:

      IF (PRESENT(hdef2)) THEN 

         IF (itype_sher.GE.1) THEN
            !Apply horizontal 3D-shear correction:
            DO k=2,kem 
!DIR$ IVDEP
               DO i=ivstart, ivend
                  frm(i,k)=frm(i,k)+hdef2(i,k) !extended shear
               END DO
            END DO
         END IF   

      END IF

      IF (lssintact) THEN !save pure turbulent shear
         DO k=2,kem
!DIR$ IVDEP
            DO i=ivstart, ivend
               ftm(i,k)=frm(i,k)
            END DO
         END DO
      END IF

!For_Tuning>
!     Preparation for Richardson-number-dependent factor used for correcting 
!     the minimum diffusion coefficient and the horizontal shear production term:

      DO k=2,ke
!DIR$ IVDEP
        DO i=ivstart, ivend
!Achtung: Mit Hilfe von 'frm' und 'frh' auszudruecken: <
!x1 = z1d2*(hhl(i,k-1)-hhl(i,k+1))
!x2 = MAX( 1.e-6_wp, ((u(i,k-1)-u(i,k))**2+(v(i,k-1)-v(i,k))**2)/x1**2 )          ! |du/dz|**2
!x3 = MAX( 1.e-5_wp, grav/(z1d2*(t(i,k-1)+t(i,k)))*((t(i,k-1)-t(i,k))/x1+tet_g) ) !       N**2

!xri(i,k) = EXP(z2d3*LOG(x2/x3))  ! 1/Ri**(2/3)

           xri(i,k)=EXP( z2d3*LOG( MAX( 1.e-6_wp, frm(i,k) ) / & !1/Ri**(2/3)
                                   MAX( 1.e-5_wp, frh(i,k) ) ) )
!>
        END DO
      END DO
!>For_Tuning

      IF (PRESENT(hdef2)) THEN

         !Additional impact by separated horizontal shear:
         IF ((ltkeshs .OR. (loutshs .AND. PRESENT(tket_hshr))) .AND. PRESENT(hdiv)) THEN 
            !Include separated horizontal shear mode:

            fakt=z1/(z2*sm_0)**2; wert=a_hshr*akt*z1d2

!DIR$ IVDEP
            DO i=ivstart, ivend
               lay(i)=wert*l_hori(i) !uncorrected effective horizontal length scale
            END DO

            IF (imode_shshear.EQ.2) THEN
!>Tuning
              DO k=2,kem
!DIR$ IVDEP
                DO i=ivstart, ivend
                  ! Factor for variable 3D horizontal-vertical length scale proportional to 1/SQRT(Ri),
                  ! decreasing to zero in the lowest kilometer above ground

                  x4 = MIN( 1._wp, 1.0e-3_wp*(hhl(i,k)-hhl(i,ke1)) ) ! low-level reduction factor
                  hor_scale(i,k) = lay(i)*MIN( 5.0_wp, MAX( 0.01_wp, x4*xri(i,k) ) )
                END DO
              END DO
!>Tuning: This kind of correction can be substituded by a less ad-hoc approach.
            ELSE
              DO k=2,kem
!DIR$ IVDEP
                DO i=ivstart, ivend
                  hor_scale(i,k) = lay(i)
                END DO
              END DO
            ENDIF

            !strain velocity (shv) of the separated horizontal shear mode:
            IF (imode_shshear.EQ.0) THEN !former variant based on 3D-shear and incompressibility
               DO k=2,kem
!DIR$ IVDEP
                  DO i=ivstart, ivend
                     shv(i,k)=hor_scale(i,k)*SQRT(hdef2(i,k)+hdiv(i,k)**2) !not equal to trace of 2D-strain tensor
                  END DO
               END DO
            ELSE !new variant in accordance with the trace constraint for the separated horizontal strain tensor
               DO k=2,kem
!DIR$ IVDEP
                  DO i=ivstart, ivend
                     wert=fakt*hdiv(i,k)
                     shv(i,k)=hor_scale(i,k)*(SQRT(wert**2+hdef2(i,k))-wert) !equal to trace of 2D-strain tensor
                  END DO
               END DO
            END IF

            DO k=2,kem 
!DIR$ IVDEP
               DO i=ivstart, ivend
                   hlp(i,k)=(shv(i,k))**3/hor_scale(i,k) !additional TKE-source by related shear
               END DO
            END DO
            IF (loutshs .AND. PRESENT(tket_hshr)) THEN
               !Load output variable for the TKE-source by separated horiz. shear:
               DO k=2,kem 
!DIR$ IVDEP
                  DO i=ivstart, ivend
                     tket_hshr(i,k)=hlp(i,k)
                  END DO
               END DO
            END IF
            IF (ltkeshs) THEN 
               !Consider separated horizontal shear mode in mechanical forcing:
               DO k=2,kem 
!DIR$ IVDEP
                  DO i=ivstart, ivend
                     frm(i,k)=frm(i,k)+hlp(i,k)/tkvm(i,k) !extended shear
                  END DO
               END DO
               IF (l3dturb) THEN
                  ! Load related horizontal diffusion coefficients:
                  fakt=sh_0/sm_0
                  DO k=2,kem
!DIR$ IVDEP
                     DO i=ivstart, ivend
                        tkhm(i,k)=hor_scale(i,k)*shv(i,k) !for momentum related to the sep. shear mode
                        tkhh(i,k)=fakt*tkhm(i,k)          !for scalars    ,,       ,,            ,,
                     END DO
                  END DO
               END IF
            END IF

         END IF   
 
      END IF

      !Addition verallgemeinerter Scherterme durch nicht-turbulente subskalige Stroemung:

      DO k=2,kem

         IF (.NOT.lini) THEN !nicht bei der Initialisierung

            IF (lsso .AND. PRESENT(ut_sso) .AND. PRESENT(vt_sso)) THEN
               !SSO-Schema ist aktiv und SSO-Tendenzen des Windes sind vorhanden:

!              Berechnung der TKE-Tendenz durch Nachlaufwirbelproduktion aus SSO-Tendenzen:

!DIR$ IVDEP
               DO i=ivstart, ivend
!Achtung: horizontale Pos. beachten!
                 vel1=-(ut_sso(i,k)  *u(i,k  )+vt_sso(i,k)  *v(i,k  ))*dp0(i,k-1)
                 vel2=-(ut_sso(i,k-1)*u(i,k-1)+vt_sso(i,k-1)*v(i,k-1))*dp0(i,k)

                 src(i)=MAX( z0, (vel1+vel2)/(dp0(i,k-1)+dp0(i,k)) )

!                Beachte:
!                Die SSO-Tendenzen beziehen sich tatsaechlich auf Massenpunkte, sie werden
!                erst spaeter in SUB 'organize_physics' auf die Zwischenpositionen interpoliert!
!                Obwohl vel1 und vel2 immer positiv sein muessten, wird zur Sicherheit MAX benutzt!
               END DO

               IF (PRESENT(tket_sso)) THEN
!DIR$ IVDEP
                  DO i=ivstart, ivend
                     tket_sso(i,k)=src(i)
                  END DO
               END IF

!              Addition des Scherterms durch Nachlaufwirbel:

               IF (ltkesso) THEN !Nachlaufwirbeltendenzen sollen beruecksichtigt werden
!DIR$ IVDEP
                  DO i=ivstart, ivend
                     frm(i,k)=frm(i,k) + src(i)/tkvm(i,k)
                  END DO
               END IF

            END IF

            IF (lconv .AND. ltkecon .AND. PRESENT(tket_conv)) THEN
               !Konvektionsschema ist aktiv, soll mit Turbulenz interagieren und conv. TKE-Tend. ist vorhanden:

!              Addition des Scherterms durch konvektive Zirkulation:

!DIR$ IVDEP
               DO i=ivstart, ivend
                  frm(i,k)=frm(i,k) + MAX( z0, tket_conv(i,k)/tkvm(i,k) )
               END DO

!              Beachte:  Obwohl tket_conv immer positiv sein muesste, wird zur Sicherheit MAX benutzt!
            END IF
         END IF

      END DO

!     Berechnung von Korrekturtermen innerhalb der Rauhigkeitsschicht
!     (ausser Volumenterme, die zur Diffusion gehoeren):

      DO k=kcm,kem !von oben nach unten durch Rauhiggkeitsschicht

!Achtung: Neue Behandlung der Rauhigkeitsschicht einfuehren

!        Vertikalwind und Formreibungskoeff. auf Hauptflaechen:

!DIR$ IVDEP
         DO i=ivstart, ivend
            lay(i)=z1d2*(w(i,k)+w(i,k+1))
            src(i)=z1d2*(cbig(i,k)  +csml(i,k) &
                        +cbig(i,k+1)+csml(i,k+1))
         END DO

!        Windbetrag auf der angrenzenden Hauptflaeche:

         IF (k.EQ.kcm) THEN
!DIR$ IVDEP
            DO i=ivstart, ivend
               velo=z1d2*(w(i,k-1)+w(i,k))
               lays(i,1)=SQRT(u(i,k-1)**2+v(i,k-1)**2+velo**2)
             END DO
         END IF

!        Reduzierte Konstanten und implizite Horizontalwind-Tendenzen durch Formreibung:

!DIR$ IVDEP
         DO i=ivstart, ivend

!           Windbetrag auf der aktuellen Hauptflaeche:
            lays(i,2)=SQRT(u(i,k)**2+v(i,k)**2+lay(i)**2)

!           Aufaddieren der Windtendenzen durch Fromreibung:
            wert=src(i)*lays(i,2)

            utens(i,k)=utens(i,k)-tinc(u_m)*wert*u(i,k)/(z1+dt_var*wert)
            vtens(i,k)=vtens(i,k)-tinc(v_m)*wert*v(i,k)/(z1+dt_var*wert)

!           Windbetrag auf Nebenflaechen:
            can(i,k)=(lays(i,2)*dp0(i,k-1)+lays(i,1)*dp0(i,k))/(dp0(i,k-1)+dp0(i,k))

!           Windbetrag unter Wirkung der Formreibung:
            velo=can(i,k)/(z1+can(i,k)*(cbig(i,k)+csml(i,k))*dt_var)

!           Addition des Scherterms durch Nachlaufwirbel an grossen Rauhigkeitselementen:
            frm(i,k)=frm(i,k)+cbig(i,k)*velo**3/tkvm(i,k)

!           Frequenz der kleinskaligen Rauhigkeitselemente:
            can(i,k)=csml(i,k)*can(i,k)

         END DO

      END DO  !k=kcm,kem !von oben nach unten druch Rauhigkeitsschicht

      !Optionale vertikale Glaettung des mechanischen Antriebs:

      IF (lcalc_frcsmot) THEN
         CALL vert_smooth ( &
              i_st=ivstart, i_en=ivend, k_tp=1, k_sf=ke1, &
              disc_mom=dicke, cur_tend=frm, vertsmot=frcsmot, smotfac=trop_mask )
      END IF

!     Optionale vertikale Glaettung des thermischen Antriebs:

      IF (lcalc_frcsmot) THEN
         CALL vert_smooth ( &
              i_st=ivstart, i_en=ivend, k_tp=1, k_sf=ke1, &
              disc_mom=dicke, cur_tend=frh, vertsmot=frcsmot, smotfac=trop_mask )
      END IF

!     Belegung von tkvh und tkvm mit den stabilitaetsabhaengigen Laengenmassen:

      DO k=2,kem
!DIR$ IVDEP
         DO i=ivstart, ivend
            tkvh(i,k)=tkvh(i,k)/tke(i,k,nvor)
            tkvm(i,k)=tkvm(i,k)/tke(i,k,nvor)
         END DO
      END DO

!return

! 3)  Loesung der turbulenten Bilanzgleichungen (Bestimmung von TKE und der Stabilitaetsfunktionen)
!     und Berechnung der turbulenten Diffusionskoeffizienten:

      DO it_durch=it_start, it_end

        !Die Schleife wird nur bei der Initialisierung (d.h. beim ersten Aufruf) wiederholt,
        !um TKE-Gleichgewicht anzunaehern. Die resultierenden TKE-Werte der Zeitstufe 'ntur'
        !gehoeren in diesem Fall dann zur Zeitstufe der uebrigen prognostischen Variablen
        !('nold' bei "leap-frog" oder 'nnow' bei 2-Zeitebenen).
        !Fuer die folgenden Aufrufe wird die Schleife nur einmal durchlaufen und liefert TKE-Werte
        !die gegenueber den Vorgaengerwerten um einen Zeitschritt weiter in der Zukunft liegen,
        !also wieder zur Zeitstufe der uebrigen prognostischen Variablen gehoeren.

         CALL solve_turb_budgets (khi=1, it_s=it_durch,                                 &
                                  i_st=ivstart, i_en=ivend, k_st=2, k_en=kem,       &
                                  kcm=kcm, ntur=ntur, nvor=nvor,                        &
                                  lssintact=lssintact, lupfrclim=.FALSE.,               &
                                  lpresedr=PRESENT(edr), lstfnct=lstfnct, ltkeinp=ltkeinp, &
                                  imode_stke=imode_turb, imode_vel_min=1,               &
                                  dt_tke=dt_tke, fr_tke=fr_tke, b_h=b_h, b_m=b_m,       &
                                  rim=rim,                                              &
                                  d_1=d_1, d_2=d_2, d_3=d_3, d_4=d_4, d_5=d_5, d_6=d_6, &
                                  tke=tke, ediss=ediss,                                 &
                                  fm2=frm, fh2=frh, ft2=ftm, lsm=tkvm, lsh=tkvh,        &
#ifdef SCLM
                                  grd=vari,                                             &
#endif
                                  fcd=can, tls=len_scale, tvt=tketens, avt=tketadv)

         IF (it_durch.LT.it_end .AND. .NOT.ltkeinp) THEN
            nvor=ntur !benutze nun aktuelle TKE-Werte als Vorgaengerwerte
         END IF

      END DO

!     Kein TKE-Gradient am Oberrand:

!DIR$ IVDEP
      DO i=ivstart, ivend
         tke(i,1,ntur)=tke(i,2,ntur)
      END DO

      IF (iini.EQ.1) THEN !only for separate initialization before the time loop

         DO k=2, kem
!DIR$ IVDEP
            DO i=ivstart, ivend
!Achtung: Bislang fehtlte die laminare Beschraenkung
               val1=MAX ( con_m, tkmmin ); tkvh(i,k)=tkvh(i,k)*tke(i,k,ntur)
               val2=MAX ( con_h, tkhmin ); tkvm(i,k)=tkvm(i,k)*tke(i,k,ntur)
               IF (imode_tkemini.EQ.2) THEN
                  tke(i,k,ntur)=tke(i,k,ntur)*MAX( z1, val1/tkvm(i,k), & !adapted tke
                                                       val2/tkvh(i,k) )
               END IF
               tkvm(i,k)=MAX( val1, tkvm(i,k) ) !corrected tkv
               tkvh(i,k)=MAX( val2, tkvh(i,k) )
            END DO
         END DO

       ! IF (.NOT.PRESENT(rhoh))  DEALLOCATE ( rhoh , STAT=ilocstat )
       ! IF (.NOT.PRESENT(epr))   DEALLOCATE ( exner, STAT=ilocstat )
       ! IF (.NOT.PRESENT(c_big)) DEALLOCATE ( cbig,  STAT=ilocstat )
       ! IF (.NOT.PRESENT(c_sml)) DEALLOCATE ( csml,  STAT=ilocstat )
       ! IF (.NOT.PRESENT(r_air)) DEALLOCATE ( rair,  STAT=ilocstat )

         RETURN !finish this subroutine

      END IF

!  4) Berechnung der effektiven turbulenten Vertikalgradienten,
!     Standardabweichnung des Saettigungsdefizites und der Diffusionskoeffizienten:

      IF (ltmpcor) THEN
      !  Berechnung des vert. Temp.grad. fuer den Phasendiffusionsterm:
         DO k=2, ke1
!DIR$ IVDEP
            DO i=ivstart, ivend
               frm(i,k)=a(i,k,1)*vari(i,k,tet_l)-tet_g  !vertical temperature gradient
            END DO
         END DO
         IF (icldm_turb.NE.-1) THEN !water phase changes are possible
            DO k=2, ke1
!DIR$ IVDEP
               DO i=ivstart, ivend
                  frm(i,k)=frm(i,k)+lhocp*vari(i,k,liq) !liquid water correction
               END DO
            END DO
         END IF

      !  Dies geschieht schon hier, weil im naechsten Schritt das Feld vari()
      !  durch die effiktiven Gradienten ueberschrieben wird.
      END IF

      DO k=2, ke1
!DIR$ IVDEP
         DO i=ivstart, ivend
            a(i,k,5)=a(i,k,1)*a(i,k,3)                           !exner*d_qsat/d_T
            a(i,k,4)=c_scld*rcld(i,k)/(z1+rcld(i,k)*(c_scld-z1)) !effective cloud cover
            rcld(i,k)=SQRT(len_scale(i,k)*tkvh(i,k)*d_h)* & !standard deviation
                           ABS(a(i,k,5)*vari(i,k,tet_l)   & !of local
                              -vari(i,k,h2o_g))             !oversaturation
         END DO   
      END DO   

      IF (icldm_turb.NE.-1) THEN !consideration of water phase changes
         DO k=2, ke1
!DIR$ IVDEP
            DO i=ivstart, ivend
             ! Effective vertical gradient of liquid water content:
   
               flw_h2o_g=a(i,k,4)/(z1+lhocp*a(i,k,3))  !weight of h2o_g-flux
               flw_tet_l=-flw_h2o_g*a(i,k,5)           !weight of tet_l-flux
   
               vari(i,k,liq)=  flw_tet_l*vari(i,k,tet_l) &                    ! eff_grad(liq)
                            +  flw_h2o_g*vari(i,k,h2o_g)   
   
             ! Effective vertical gradient of water vapor content and pot. temper.:
   
               vari(i,k,vap)= vari(i,k,h2o_g)-vari(i,k,liq)                   ! eff_grad(vap)
               vari(i,k,tet)=(vari(i,k,tet_l)+vari(i,k,liq)*lhocp/a(i,k,1)) & ! eff_grad(tet) 
                             *a(i,k,2)                                        !*(Cp/Cpd)

             ! Note: 
             ! -flux(h2o_g)/(rho*K)=grad(h2o_g)=eff_grad(vap)+eff_grad(liq)    
             ! -flux(tet_l)/(rho*K)=grad(tet_l)=eff_grad(tet)-eff_grad(liq)*lh_v/(cp_d*exner)
             ! Treating "lh_v/(cp_d*exner)" like a constnat, besides numerical effects, 
             !  vertical gradient diffusion of non conserved variables without a moist gradient correction
             !  does not change the resulting diffusion of conserved variables. Thus the redistribution of
             !  water constituents can (or should) be left to a final (sub grid scale) saturation adjustment.
            END DO
         END DO   
      ELSE !no water phase change possible
        !'vari(:,:,liq)' bleibt unveraendert, genauso wie auch
        !'vari(:,:,tet)'='vari(:,:,tet_l)' und 'vari(:,:,vap)'='vari(:,:,h2o_g)'.
         DO k=2, ke1
!DIR$ IVDEP
            DO i=ivstart, ivend
               vari(i,k,tet)=vari(i,k,tet_l)*a(i,k,2) ! grad(tet)*(Cp/Cpd)
            END DO
         END DO
      END IF 

      !Beachte:
      !Die eff_grad ergeben multipliziert mit 'rho*tkvh' die vertikalen Flussdichten
      ! von 'tet', 'vap' und 'liq' unter Beruecksichtigung der turbulenten Phasenuebergaenge.
      !'vari(:,:,tet)' ist der in der Teta-Gleichung benoetigte effective Gradient 
      ! und entspricht der 'tet'-Flussdichte*(Cp/Cpd).
      !Die Indices 'tet', 'tem' und 'tet_l' haben den gleichen Wert
      ! so wie auch 'vap' und 'h2o_g' den gleichen Wert haben.
      !Die Unterscheidungen sollen nur den jeweiligen Inhalt verdeutlichen.
      !Im Falle "icldm_turb.EQ.EQ.-1" verschwindet der eff. Bedeckungsgrad "a(:,:,4)=0",
      ! troltzdem steht auch dann die Standardabweichung der Uebersaettigung in 'rcld'.

!     Beschraenkung der Diffusionskoeffizienten:

      val1=tkmmin; val2=tkhmin !default minimum diffusion values

      DO k=2, ke
!DIR$ IVDEP
         DO i=ivstart, ivend
!           Berechn. der Diffusionskoeffizienten:

            IF (imode_tkvmini.EQ.2) THEN
!>Tuning
               ! Factor for variable minimum diffusion coefficient proportional to 1/SQRT(Ri);
               ! the namelist parameters tkhmin/tkmmin specify the value for Ri=1:
               fakt=MIN( z1, tkred_sfc(i)*(0.25_wp+7.5e-3_wp*(hhl(i,k)-hhl(i,ke1))) ) !low-level red.-fact.
!US old        fakt=MIN( z1,              (0.25_wp+7.5e-3_wp*(hhl(i,k)-hhl(i,ke1))) ) !low-level red.-fact.
               fakt=MIN( 2.5_wp, MAX( 0.01_wp, fakt*xri(i,k) ) )

               val1=tkmmin*fakt; val2=tkhmin*fakt

               IF (tkhmin_strat.GT.z0 .OR. tkmmin_strat.GT.z0) THEN
                  ! Enhanced diffusion in the stratosphere - very important for the data assimilation cycle,
                  ! but can also be used in forecasting mode because there is no detectable detrimental
                  ! impact on gravity waves:
                  fakt = MIN( z1, 2.e-4_wp*MAX( z0, hhl(i,k) - 25000._wp ) ) !lin. incr. betw. 25 and 30 km
                  fakt = fakt*MIN( 7.5_wp, MAX( 0.125_wp, xri(i,k) ) )

                  val1=MAX( val1, tkmmin_strat*fakt ) ; val2=MAX( val2, tkhmin_strat*fakt )
               END IF
!>Tuning: This kind of correction can be substituded by a less ad-hoc approach.
            END IF

!Achtung: Beschraenkung mit lam. diff.coef. fehlte bislang auch in ICON; macht ev. Unterschiede
            val1=MAX ( con_m, val1 ); tkvm(i,k)=tkvm(i,k)*tke(i,k,ntur)
            val2=MAX ( con_h, val2 ); tkvh(i,k)=tkvh(i,k)*tke(i,k,ntur)
            IF (imode_tkemini.EQ.2) THEN
               tke(i,k,ntur)=tke(i,k,ntur)*MAX( z1, val1/tkvm(i,k), & !adapted tke
                                                    val2/tkvh(i,k) )
            END IF
            tkvm(i,k)=MAX( val1, tkvm(i,k) ) !corrected tkv's
            tkvh(i,k)=MAX( val2, tkvh(i,k) )

!test: ohne tk?min (wird aber bei Diffusion benutzt!)
!test: ohne tk?min und ohen lam. Diffkof. (wird aber bei Diffusion benutzt!)

!           tkvh und tkvm enthalten jetzt nicht mehr Diffusions-
!           laengen, sondern wieder Diffusionskoeffizienten in m^2/s!

         END DO
      END DO

      IF (l3dturb) THEN
         !Consider horizontal diffusion coefficients:
         IF (PRESENT(hdef2) .AND. PRESENT(hdiv) .AND. ltkeshs) THEN
            !Add isotropic turbulent part to that part due to the sep. horiz. shear mode:
            DO k=2,kem
!DIR$ IVDEP
               DO i=ivstart, ivend
                  tkhm(i,k)=tkhm(i,k)+tkvm(i,k)
                  tkhh(i,k)=tkhh(i,k)+tkvh(i,k)
               END DO
            END DO
         ELSE !no treatment of sep. horiz. shear mode has taken place
            !Load only the isotropic turbulent part:
            DO k=2,kem
!DIR$ IVDEP
               DO i=ivstart, ivend
                  tkhm(i,k)=tkvm(i,k); tkhh(i,k)=tkvh(i,k)
               END DO
            END DO
         END IF
      END IF

!  5) Untere Randwerte der Turbulenzgroessen:

      !Beachte: Auch wenn 'turbtran' nicht als Transfer-Schema benutzt wird,
      !         muessen 'tkvm', tkvh' und 'tke' und (falls es PRESENT ist)
      !         auch 'edr' fuer "k=ke1" belegt sein!

      IF (lsfli(tem)) THEN !use explicit shfl_s
!DIR$ IVDEP
         DO i=ivstart, ivend
            vari(i,ke1,tet)=shfl_s(i)/(cp_d*rhon(i,ke1)*tkvh(i,ke1)*a(i,ke1,1))
         END DO
         !Note:'vari(tet)' belongs to potential temperature,
         !     and 'a(1)' contains Exner factor on half levels.
      END IF
      IF (lsfli(vap)) THEN !use explicit qvfl_s
!DIR$ IVDEP
         DO i=ivstart, ivend
            vari(i,ke1,vap)=qvfl_s(i)/(rhon(i,ke1)*tkvh(i,ke1))
         END DO
      END IF
      !Note: "tkvh(ke1) > 0.0" is required!

! 6)  Berechnung der zu TKE-Quellen gehoerigen Temperaturtendenzen
!     (ausser der Divergenz des Zirkulationsflusses):

      IF (ltmpcor) THEN
         IF (.NOT.PRESENT(edr)) THEN
!DIR$ IVDEP
            DO i=ivstart, ivend
               ediss(i,ke1)=tke(i,ke1,ntur)**3/(d_m*len_scale(i,ke1))
            END DO
         END IF

         DO k=2, ke1
!DIR$ IVDEP
            DO i=ivstart, ivend
               thermik=tkvh(i,k)*frh(i,k)
               phasdif=tkvh(i,k)*frm(i,k) &
                      *(tur_rcpv*vari(i,k,vap)+tur_rcpl*vari(i,k,liq))

               tketens(i,k)=len_scale(i,k)/a(i,k,2) &
                       *((ediss(i,k)+thermik)/cp_d+phasdif)
            END DO

!           Beachte:
!           In 'frm' steht an dieser Stelle der Temp.-Gradient!
!           Wegen der nachfolgenden Interpolation auf Hauptflaechen,
!           wird 'tketens' mit der turbulenten Laengenskala multipliziert.
         END DO

!DIR$ IVDEP
         DO i=ivstart, ivend
            ttens(i,1)=ttens(i,1)+tinc(tem)*tketens(i,2) &
                         /(len_scale(i,1)+len_scale(i,2))

         END DO
         DO k=2,ke
!DIR$ IVDEP
            DO i=ivstart, ivend
               ttens(i,k)=ttens(i,k)+tinc(tem)*(tketens(i,k)+tketens(i,k+1)) &
                            /(len_scale(i,k)+len_scale(i,k+1))
            END DO
         END DO
      END IF

! 7)  Bestimmung des Zirkulationstermes als zusaetzliche TKE-Flussdichte:
!Achtung: Zirkulationsterm revidieren:

      IF (lcircterm) THEN !Der Zirkulationsterm muss berechnet werden

         DO k=2,ke1

            IF (k.LT.ke1) THEN
!              Interpolation des Druckes auf die naechst hoehere
!              Nebenflaeche:

!DIR$ IVDEP
               DO i=ivstart, ivend
                  lay(i)=(prs(i,k)*dp0(i,k-1)+prs(i,k-1)*dp0(i,k)) &
                          /(dp0(i,k)+dp0(i,k-1))
               END DO
            ELSE
!DIR$ IVDEP
               DO i=ivstart, ivend
                  lay(i)=ps(i)
               END DO
            END IF

!DIR$ IVDEP
            DO i=ivstart, ivend

!Achtung: Variation durch Wolkenbedeckung
               fakt=z1-z2*ABS(a(i,k,4)-z1d2)
               len=MAX( l_pat(i), SQRT(fakt*len_scale(i,k)*l_hori(i)) )

!              Berechnung der lokalen Kohaerenzlaenge fuer die Parameterisierung
!              des Zirkulationstermes (kohae_len=l_pat*exnr*grad(tet_v)*r/g):

               fakt=frh(i,k)*lay(i)/(rhon(i,k)*grav**2)
               kohae_len=len*SIGN(z1,fakt)*MIN( ABS(fakt), z1 )

!              Belegung von 'frh' mit der TKE-Flussdichte des Zirkulationstermes
!              (Zirkulationsfluss oder Drucktransport):

               frh(i,k)=tkvh(i,k)*kohae_len*frh(i,k)

!              Die Divergenz dieser Flussdichte ist gleichzeitig
!              eine weitere Quelle fuer thermische Energie.

            END DO

!           Addition des Teta-Gradienten, welcher zur Teta-Flussdichte
!           durch den Zirkulationsterm gehoert:

            IF (ltmpcor) THEN
!DIR$ IVDEP
               DO i=ivstart, ivend
                  vari(i,k,tet)=vari(i,k,tet) &
                                -frh(i,k)/(tkvh(i,k)*a(i,k,1)*a(i,k,2)*cp_d)
               END DO
            END IF

            !Die kinetische Energie der Zirkulation wird der thermischen Energie
            !entzogen!

         END DO

      END IF   

! 8) Berechnung der Diffusionstendenz von q=SQRT(2*TKE) einschliesslich der
!    q-Tendenz durch den Zirkulationsterm:

!     Vorbereitung zur Bestimmung der zugehoerigen Incremente von TKE=(q**2)/2:

      cur_prof => hlp
      upd_prof => a(:,:,1)
      sav_prof => a(:,:,5)


      IF (ldotkedif .OR. lcircdiff) THEN

         expl_mom => a(:,:,2)
         impl_mom => a(:,:,3)
         invs_mom => a(:,:,4)

         ! Diffusions-Koeffizienten auf NF:
         DO k=2, ke1
            IF (imode_tkediff.EQ.2) THEN !Diffusion in terms of TKE
!DIR$ IVDEP
               DO i=ivstart, ivend
!___________________________________________________________________________
!test: TKE-Diffusion mit Stab.fnkt. fuer Skalare:
! sav_prof(i,k)=c_diff*tkvh(i,k)
                  sav_prof(i,k)=c_diff*len_scale(i,k)*tke(i,k,ntur)
!___________________________________________________________________________
               END DO
            ELSE !Diffusion in terms of q=SQRT(2*TKE)
!DIR$ IVDEP
               DO i=ivstart, ivend
                  sav_prof(i,k)=c_diff*len_scale(i,k)*tke(i,k,ntur)**2
               END DO
            END IF
         END DO

         DO k=3, ke1
!DIR$ IVDEP
            DO i=ivstart, ivend
               expl_mom(i,k)=rhoh(i,k-1)*z1d2*(sav_prof(i,k-1)+sav_prof(i,k)) &
                                             /(hhl(i,k-1)-hhl(i,k))
            END DO
            !Beachte: 
            !'expl_mom' bezieht sich auf HF, also die Fluss-Niveaus fuer die TKE (bzw. q-)-Diffusion.
            !Wegen der spaeteren Nutzung der SUBs 'prep_impl_vert_diff' und 'calc_impl_vert_diff'
            ! muss ein Fluss-Niveau (hier HF) ueber dem Variabl.-Niveau (hier NF) mit gleichem Index liegen.
!Achtung: In der COSMO-Version werden 'rhon', 'len_scale' und 'tke' einzeln auf HF interpoliert,
!         was geringe Unterschiede verursacht!
         END DO

      END IF

      IF (ldotkedif .OR. lcircterm) THEN
         IF (imode_tkediff.EQ.2) THEN !Diffusion in terms of TKE
            DO k=2, ke1
!DIR$ IVDEP
               DO i=ivstart, ivend
                  sav_prof(i,k)=z1d2*tke(i,k,ntur)**2 !TKE
               END DO
            END DO
         ELSE !Diffusion in terms of q=SQRT(2*TKE)
            DO k=2, ke
!DIR$ IVDEP
               DO i=ivstart, ivend
                  sav_prof(i,k)=tke(i,k,ntur)         !q=SQRT(2*TKE)
                  dicke(i,k)=dicke(i,k)*tke(i,k,ntur) !related effective discretization momentum
               END DO
            END DO
            !Das Feld 'dicke' wird bei "k=ke1" nicht benoetigt und war zuvor dort auch nicht belegt!
            DO i=ivstart, ivend
               sav_prof(i,ke1)=tke(i,ke1,ntur) !q at the surface layer 
            END DO
         END IF
         !Beachte: 
         !Das Feld 'tke' enthaelt nicht TKE sondern q!
      END IF

!     Aufnahme des Zirkulationstermes:

      IF (lcircterm) THEN !Der Zirkulationsterm muss berechnet werden

!        Interpolation der Zirulationsflussdichte auf Hauptflaechen:

         !Beachte:
         !'frm' ist frei. 
         !'frh' enthaelt bislang eine TKE-Flussdichte (bis auf den Faktor 'rhon'),
         ! deren Vertikalprofil im wesentlichen durch d_z(tet_v)**2 bestimmt ist,
         ! was zumindest in der Prandtl-Schicht prop. zu 1/len_scale ist.
         !Die Interpolation von "rhon*frh" auf Hauptflaechen erfolgt daher mit
         ! 'rhon*frh*len_scale'. Anderenfalls ist mit grossen Interpolationsfehlern zu rechnen.

         DO k=2,ke1
!DIR$ IVDEP
            DO i=ivstart, ivend
!Achtung: In COSMO-Version ist "imode_tkediff=1".
!         Es wird in dieser Version aber 'frh' an dieser Stelle durch 'tke' dividiert,
!         was geringe Unterschiede verursacht:"
               frh(i,k)=rhon(i,k)*frh(i,k)*len_scale(i,k) !skalierte Flussdichte auf NF
!frh(i,k)=rhon(i,k)*frh(i,k)/tke(i,k,ntur)*len_scale(i,k) !skalierte Flussdichte auf NF

!if (k.lt.ke1) then
!   frh(i,k)=rhon(i,k)*grav*tkvh(i,k)*t(i,k)/t_g(i)/tke(i,k,ntur)*len_scale(i,k)
!else
!   frh(i,k)=rhon(i,k)*grav*tkvh(i,k)/tke(i,k,ntur)*len_scale(i,k)
!endif
            END DO
         END DO

!        Korrektur der TKE-Profile durch die Zirkulations-Tendenz:

         IF (imode_calcirc.EQ.1) THEN !expliziten Berechnung der Zirkulationstendenz

            cur_prof => sav_prof

            DO k=3,ke1
!DIR$ IVDEP
               DO i=ivstart, ivend
                  frm(i,k)=(frh(i,k)+frh(i,k-1))/(len_scale(i,k)+len_scale(i,k-1)) !interpolierte Flussdichte auf HF
               END DO
            END DO

!Achtung: In COSMO-Version ist "imode_tkediff=1". 
!Da 'frh' bereits durch 'tke' dividiert wurde, muss fuer die exakte COSMO-Version der 'tke'-Faktor in 'dicke'
!wieder beseitigt werden:
            k=2
!DIR$ IVDEP
            DO i=ivstart, ivend
               upd_prof(i,k)=sav_prof(i,k)+frm(i,k+1)/dicke(i,k)
!upd_prof(i,k)=sav_prof(i,k)+frm(i,k+1)*tke(i,k,ntur)/dicke(i,k)

!upd_prof(i,k)=sav_prof(i,k)
!tketens(i,k)=frm(i,k+1)*tke(i,k,ntur)/dicke(i,k)

            END DO
            DO k=3,ke
!DIR$ IVDEP
               DO i=ivstart, ivend
                  upd_prof(i,k)=sav_prof(i,k)-(frm(i,k)-frm(i,k+1))/dicke(i,k)
!upd_prof(i,k)=sav_prof(i,k)-(frm(i,k)-frm(i,k+1))*tke(i,k,ntur)/dicke(i,k)

!Achtung: In der COSMO-Version wird die explizite q-Tendenz durch den Zirkulationsterm
!erst nach der TKE-Diffusion hinzuaddiert (verursacht geringe Unterschiede):
!upd_prof(i,k)=sav_prof(i,k)
!tketens(i,k)=-(frm(i,k)-frm(i,k+1))*tke(i,k,ntur)/dicke(i,k)

               END DO
            END DO
            !Beachte:
            !'dicke' enthaelt bereits den Faktor "1/dt_tke" (sowie den Faktor 'tke' bei "imode_tkediff=1").
            !'upd_prof' enthaelt das um die Zirkulations-Tendenz aufdatierte TKE-Profil
            ! (oder q-Profile bei "imode_tkediff=1")

            !Zuschlag durch Volumenterm aus der Divergenzbildung:
            DO k=ke,kcm,-1 !innerhalb der Rauhigkeitsschicht
!DIR$ IVDEP
               DO i=ivstart, ivend
!Achtung: Korrektur
                ! upd_prof(i,k)=upd_prof(i,k)+dt_tke*frh(i,k)*z1d2*(rair(i,k-1)-rair(i,k+1)) &
                  upd_prof(i,k)=upd_prof(i,k)+frh(i,k)*z1d2*(rair(i,k-1)-rair(i,k+1)) &
                                                           /(len_scale(i,k)*dicke(i,k))
                  !'frh' enthaelt die Zirkultions-Flussdichte der TKE auf NF (skaliert mit 'len_scale').
               END DO
            END DO

            !Bereucksichtige Zirkulations-Tendenz:
             itndcon=1 !indem 'upd_prof' auf rechter Seite der impliz. Diff.-Gl. benutzt wird.
!Achtung: Um die COSMO-Version exakt nachzubilden, darf die Zirkulationstendenz nicht bei der
!impliziten Diffusions-Gleichung benutzt werden:
!itndcon=0
            !Fuer die expliziten Diff.-Anteile wird aber 'cur_prof' benutzt.

         ELSE !quasi implizite Berechnung der Zirkulatinstendenz (entspricht "lcircdiff=T")
!DIR$ IVDEP
            DO i=ivstart, ivend
               cur_prof(i,2)=sav_prof(i,2)
            END DO
            DO k=3,ke1
!DIR$ IVDEP
               DO i=ivstart, ivend
                  wert=(frh(i,k)+frh(i,k-1))/((len_scale(i,k)+len_scale(i,k-1))*expl_mom(i,k))
                  cur_prof(i,k)=(cur_prof(i,k-1)-sav_prof(i,k-1)+wert)+sav_prof(i,k)
               END DO
            END DO
            !Beachte:
            !'cur_prof' enthaelt ein virtuelles TKE-Profil (oder q-Profile bei "imode_tkediff=1"), 
            ! dessen Diffusions-Tendenz die Zirkulations-Tendenz einschliesst.

            !Bereucksichtige Zirkulations-Tendenz:
            itndcon=0 !indem 'cur_prof' auf rechter Seite der impliz. Diff.-Gl. benutzt wird.
            !Fuer die expliziten Diff.-Anteile wird ebenfalls 'cur_prof' benutzt.
         END IF   

      ELSEIF (ldotkedif) THEN
 
         cur_prof => sav_prof

         itndcon=0 !'cur_prof' wird auf rechter Seite der impliz. Diff.-Gl.
                   ! und fuer explizite Diff.-Anteile benutzt.
      END IF

!     Aufdatieren des TKE-Profils durch die (erweiterte) Diffusions-Tendenz 

      IF (ldotkedif .OR. lcircdiff) THEN

         !'frm', 'frh' und 'len_scale' sind frei.
         !In den Diffusionroutinen wird vorausgesetzt, dass ein Flussniveau mit gleichem
         !Vertikalindex wie ein Konzentrationsniveau gerade ueber letzterem liegt.
         !Die bisherige Hauptflaechenindizierung musste daher fuer die Uebergabefelder
         !der Routinen 'prep_impl_vert_diff' und 'calc_impl_vert_diff' um "1" verschoben werden.

         CALL prep_impl_vert_diff( lsflucond=.FALSE., ldynimpwt=ldynimp, lprecondi=lprecnd, &
                                   i_st=ivstart, i_en=ivend, k_tp=1, k_sf=ke1,  &
!Achtung: q_Diff
              disc_mom=dicke, expl_mom=expl_mom, impl_mom=impl_mom, invs_mom=invs_mom, &
!disc_mom=sav_prof*dicke, expl_mom=expl_mom, impl_mom=impl_mom, invs_mom=invs_mom, &
              invs_fac=frh, scal_fac=frm, impl_weight=impl_weight )

!        Berechnung der vertikalen Diffusionstendenzen von TKE=z1d2*q**2:

         eff_flux => len_scale

         CALL calc_impl_vert_diff ( lsflucond=.FALSE.,lprecondi=lprecnd, &
                                    leff_flux=(kcm.LE.ke), itndcon=-itndcon, &
                                    i_st=ivstart, i_en=ivend,k_tp=1, k_sf=ke1, &
!Achtung: q_Diff
              disc_mom=dicke, expl_mom=expl_mom, impl_mom=impl_mom, invs_mom=invs_mom, &
!disc_mom=sav_prof*dicke, expl_mom=expl_mom, impl_mom=impl_mom, invs_mom=invs_mom, &
              invs_fac=frh, scal_fac=frm, cur_prof=cur_prof, upd_prof=upd_prof, eff_flux=eff_flux )

         !Beachte:
         !Bei "imode_diff=1" erfolgt q-Diffusion, so dass 'disc_mom' und 'expl_mom' den zusaetzlichen 
         ! Faktor 'sav_prof'='tke'=q enthalten!
         !'upd_prof' enthaelt jetzt die mit der Diffusionstendenz aufdatierten (modifizierte) TKE-Werte.
         !Weil "itndcon<=0", bleiben auf 'cur_prof' die Eingangsprofile erhalten.
         !'eff_flux' enthaelt die effektiven Flussdichten (positiv abwaerts) der (semi-)impliziten
         ! Vertikaldiffusion.

         IF (lcircdiff) THEN !es wurden virtuelle Effektiv-Profile benutzt
            DO k=2,ke
!DIR$ IVDEP
               DO i=ivstart, ivend
                  upd_prof(i,k)=sav_prof(i,k)+upd_prof(i,k)-cur_prof(i,k) !aufdatierte echte Profile
               END DO
            END DO
         END IF
!++++

!        Zuschlag durch Volumenterm innerhalb der Rauhigkeitsschicht:

         DO k=ke,kcm,-1 !innerhalb der Rauhigkeitsschicht
!DIR$ IVDEP
            DO i=ivstart, ivend
               wert=(eff_flux(i,k)*dp0(i,k)+eff_flux(i,k+1)*dp0(i,k-1))/(dp0(i,k)+dp0(i,k-1))
                   !effektive TKE-Flussdichte interpoliert auf die k-te Nebenflaeche,
                   !wobei sich 'eff_flux(:,k)' auf die (k-1)-te Hauptflaeche bezieht!
!Achtung: Korrektur
             ! upd_prof(i,k)=upd_prof(i,k)+dt_tke*wert*z1d2*(rair(i,k-1)-rair(i,k+1))/dicke(i,k)
               upd_prof(i,k)=upd_prof(i,k)+wert*z1d2*(rair(i,k-1)-rair(i,k+1))/dicke(i,k)
            END DO
         END DO
         !'upd_prof' enthaelt das mit dem Volunmenterm beaufschlagte aufdatierte Profil.

      END IF     

! 9)  Speichern der zugehoerigen q-Tendenzen:

      IF (ldotkedif .OR. lcircterm) THEN   

         IF (imode_tkediff.EQ.2) THEN !Diffusion in terms of TKE
            !'upd_prof' ist ein TKE-Profil:
            DO k=2,ke 
!DIR$ IVDEP
               DO i=ivstart, ivend
!___________________________________________________________________________
!test:
                ! tketens(i,k)=( MAX( upd_prof(i,k), z0 ) - sav_prof(i,k) )*fr_tke/tke(i,k,ntur)
                  tketens(i,k)=( SQRT( 2*MAX( upd_prof(i,k), z0 ) ) - tke(i,k,ntur) )*fr_tke
!___________________________________________________________________________
               END DO 
            END DO
         ELSE !Diffusion in terms of q=SQRT(2*TKE)
            !'upd_prof' ist ein q-Profil:
            DO k=2,ke
!DIR$ IVDEP
               DO i=ivstart, ivend
!Achtung:
                 tketens(i,k)=( MAX( upd_prof(i,k), z0 ) - tke(i,k,ntur) )*fr_tke

!Achtung: Bei der COSMO-Version wird erst hier die explizite Zirkulations-Tendenz hinzuaddiert:
!tketens(i,k)=MAX( -tke(i,k,ntur), tketens(i,k)+( upd_prof(i,k) - tke(i,k,ntur) ) )*fr_tke
               END DO
            END DO
         END IF
         !'tketens' enthaelt jetzt immer eine q-Tendenz!

         !Am Unterrand gibt es keine q-Tendenz durch Diffusionsterme:
!DIR$ IVDEP
         DO i=ivstart, ivend
            tketens(i,ke1)=z0
         END DO
                
!        Optionale vertikale Glaettung der erweiterten Diffusionstendenz von q=SQRT(2*TKE):

         IF (tndsmot.GT.z0) THEN
            CALL vert_smooth ( &
                 i_st=ivstart, i_en=ivend, k_tp=1, k_sf=ke1, &
                 disc_mom=dicke, cur_tend=tketens, vertsmot=tndsmot )
         END IF

      ELSE !keine q-Tendenzen, weder durch TKE-Diffusion noch durch den Zirkulationsterm

!        Zuruecksetzen der q-Tendenzen:

         DO k=2,ke1
!DIR$ IVDEP
            DO i=ivstart, ivend
               tketens(i,k)=z0
            END DO
         END DO

      END IF


! 10) Interpolationen auf Hauptflaechen fuer die Standardabweichnung
!     des Saettigungsdefizites:

!DIR$ IVDEP
      DO i=ivstart, ivend
         rcld(i,1)=rcld(i,2)
      END DO
      DO k=2,kem-1
!DIR$ IVDEP
         DO i=ivstart, ivend
            rcld(i,k)=(rcld(i,k)+rcld(i,k+1))*z1d2
         END DO
      END DO
!     Fuer die unterste Hauptflaeche (k=ke) wird bei kem=ke
!     der Wert auf der entspr. Nebenflaeche beibehalten.

!-----------------------
   END IF !lturatm   
!-----------------------

!########################################################################

!Achtung: ".NOT.lini" ergaenzt
!--------------------------------------------------
      IF ((ldovardif .OR. ldogrdcor) .AND. .NOT.lini) THEN !Vertikaldiffusion wird hier berechnet
!--------------------------------------------------
         !Note: 
         !If ".NOT.ldovardif .AND. ldogrdcor", only a correction of pure vertical gradient diffusion
         ! due to sub grid scale condensation or non-local gradients is performed.
      
!        Berechnung der Luftdichte und des Exner-Faktors am Unterrand:

!DIR$ IVDEP
         DO i=ivstart, ivend
            virt=z1+rvd_m_o*qv_s(i) !virtueller Faktor
            rhon(i,ke1)=ps(i)/(r_d*virt*t_g(i))
            eprs(i,ke1)=zexner(ps(i))
         END DO
         !Note:
         !In the turbulence model 'rhon(:,ke1)' belongs to the lower boundary of the
         !Prandtl-layer, rather than to the surface level.

!        Berechnung von Hilfsgroessen:

         IF (.NOT.lturatm) THEN !turbulence model was not running

            IF (.NOT.PRESENT(rho)) THEN
               DO k=1,ke
!DIR$ IVDEP
                  DO i=ivstart, ivend
                     virt=z1+rvd_m_o*qv(i,k)-qc(i,k) !virtueller Faktor
                     rhoh(i,k)=prs(i,k)/(r_d*virt*t(i,k))
                  END DO
               END DO
            END IF

            CALL bound_level_interp( ivstart, ivend, 2,ke, &
!___________________________________________________________________________
!test: mass weighted interpolation
!                              nvars=1, pvar=(/varprf(rhon,rhoh)/), depth=hhl, auxil=hlp)
                               nvars=1, pvar=(/varprf(rhon,rhoh)/), depth=dp0)
!___________________________________________________________________________

            IF (lscadif .AND. .NOT.PRESENT(epr)) THEN
               DO k=1, ke
!DIR$ IVDEP
                  DO i=ivstart, ivend
                     exner(i,k)=zexner(prs(i,k))
                  END DO
               END DO
            END IF

         END IF

!        Setzen von Steuerparametern:

         IF (ldogrdcor) THEN
            IF (lnonloc) THEN
               ncorr=1
            ELSE
               ncorr=nvel+1
            END IF
         ELSE
            ncorr=ndiff+1
         END IF

         ivtype=0

!-----------------------------------------------------------------
!        Berechnung der Vertikaldiffusion von Modellvariablen auf Hauptflaechen:
!-----------------------------------------------------------------

!        DO n=nprim, nlast !loop over all variables to be diffused
         DO n=1, ndiff !loop over all variables to be diffused potentially
         IF ( (lum_dif .AND. n.EQ.u_m)   .OR. &                   !u_m-diffusion or
              (lvm_dif .AND. n.EQ.v_m)   .OR. &                   !v_m-diffusion or
              (lscadif .AND. n.GT.nvel)  .OR. &                   !sca-diffusion or
            (ldogrdcor .AND. n.GE.ncorr .AND. n.LE.nmvar) ) THEN !gradien correction

            m=MIN(n,nmvar)

            IF (ivtype.EQ.0) THEN
               linisetup=.TRUE.; lnewvtype=.TRUE.
            ELSE
               linisetup=.FALSE.;lnewvtype=.FALSE.
            END IF

            IF (n.LE.nvel) THEN !a wind component
               IF (ivtype.EQ.sca) lnewvtype=.TRUE.

!Achtung:
               lsflucond=.FALSE. !no lower flux condition for momentum!
!lsflucond=lsflcnd !no lower flux condition for momentum!

               ivtype=mom
            ELSE
               IF (ivtype.EQ.mom) lnewvtype=.TRUE.

               lsflucond=lsflcnd !use chosen type of lower boundary condition

               ivtype=sca
            END IF

            IF (n.LT.ncorr .OR. n.GT.nmvar) THEN 
               igrdcon=0 !keine Gradientkorrektur der Profile
            ELSEIF ( (.NOT.lscadif .AND. ivtype.EQ.sca) .OR. &
                     (.NOT.lum_dif .AND.      n.EQ.u_m) .OR. &
                     (.NOT.lvm_dif .AND.      n.EQ.v_m) ) THEN
               igrdcon=1 !nur Profile aus Gradientkorrektur
            ELSE
               igrdcon=2 !korrigierte Profile aus effektiven Gradienten
            END IF

            kgc=2 !uppermost level of gradient correction

            IF (.NOT.ltend(n)) THEN !tendency array not present
               itndcon=0 !no explicit tendency consideration
            ELSE
               itndcon=itnd !use chosen mode of tendency consideration
            END IF
!test: never expl tendency consideration
!itndcon=0
!test

            IF (lsfli(n) .AND. (.NOT.lturatm .OR. (n.NE.tem .AND. n.NE.vap))) THEN
               !Load effective surface layer gradients due to given flux values:

!DIR$ IVDEP
               DO i=ivstart, ivend
                  vari(i,ke1,m)=dvar(n)%sv(i)/(rhon(i,ke1)*vtyp(ivtype)%tkv(i,ke1))
               END DO
               IF (n.EQ.tem) THEN !flux density is that of sensible heat
!DIR$ IVDEP
                  DO i=ivstart, ivend
                     vari(i,ke1,m)=vari(i,ke1,m)/(cp_d*eprs(i,ke1))
                  END DO
               END IF
               !Note:
               !In this case not the current surface concentration but the current flux-density
               ! at the surface is used in 'vert_grad_diff'!
               !Hoewever, 'vari' contains vertical gradients at this place, and for ".NOT.lsflucond"
               ! a related surface concentration is recalculated in 'vert_grad_diff'.
               !'tkv(ke1)' needs to be >0, which is always the case, if calculated by 'turbtran'!
               !Thus "tkvh(ke1)=0.0" should not be forced, if "lsfli=.TRUE"!
               !For tracers it is "m=nmvar"!
               !In case of "lturatm=T" 'vari(ke1)' has already been loaded using shlf_s or qvfl_s!
            END IF

!           Belegung der Eingangsprofile und -Tendenzen:

            cur_prof => hlp

            DO k=1,ke
!DIR$ IVDEP
               DO i=ivstart, ivend
                  cur_prof(i,k)=dvar(n)%av(i,k)
               END DO
            END DO

            IF (ASSOCIATED(dvar(n)%sv)) THEN !surface variable is present
!DIR$ IVDEP
               DO i=ivstart, ivend
                  cur_prof(i,ke1)=dvar(n)%sv(i)
               END DO
            ELSEIF (n.LE.nvel .OR. ilow_def_cond.EQ.2) THEN
               !No-slip-condition for momentum or zero-concentr.-condition as a default:
!DIR$ IVDEP
               DO i=ivstart, ivend
                  cur_prof(i,ke1)=z0
               END DO
            ELSE !enforce a zero flux condition as a default
!DIR$ IVDEP
               DO i=ivstart, ivend
                  cur_prof(i,ke1)=cur_prof(i,ke)
               END DO
            END IF
            IF (itndcon.GT.0) THEN !explicit tendencies have to be considered
               DO k=1,ke
!DIR$ IVDEP
                  DO i=ivstart, ivend
                     dicke(i,k)=dvar(n)%at(i,k)
                  END DO
               END DO
            END IF

            IF (n.EQ.tem) THEN !temperature needs to be transformed
               DO k=1,ke
!DIR$ IVDEP
                  DO i=ivstart, ivend
                     cur_prof(i,k)=cur_prof(i,k)/exner(i,k) !potential temperature
                  END DO
               END DO      
!DIR$ IVDEP
               DO i=ivstart, ivend
                  cur_prof(i,ke1)=cur_prof(i,ke1)/eprs(i,ke1)
               END DO
               IF (itndcon.GT.0) THEN !explicit tendencies to be considered
                  DO k=1,ke
!DIR$ IVDEP
                     DO i=ivstart, ivend
                        dicke(i,k)=dicke(i,k)/exner(i,k)
                     END DO
                  END DO
               END IF    
            END IF

            IF (.NOT.(lsfluse .AND. lsflcnd)) THEN ! calculation of effective flux density required
               IF ( ( n.EQ.tem .AND. PRESENT(shfl_s) ) .OR. ( n.EQ.vap .AND. PRESENT(qvfl_s) ) ) THEN
                  leff_flux = .TRUE.
               ELSE
                  leff_flux = .FALSE.
               END IF
            ELSE
               leff_flux = .FALSE.
            END IF

!           Berechnung der vertikalen Diffusionstendenzen:

            CALL vert_grad_diff( kcm, kgc,                            &
!
                 i_st=ivstart, i_en=ivend, k_tp=0, k_sf=ke1,          &
!
                 dt_var=dt_var, ivtype=ivtype, igrdcon=igrdcon, itndcon=itndcon, &
!
                 linisetup=linisetup, lnewvtype=lnewvtype,            &
                 lsflucond=lsflucond, lsfgrduse=lsfli(n),             &
                 ldynimpwt=ldynimp  , lprecondi=lprecnd,              &
                 leff_flux=leff_flux,                                 &
!
                 rho=rhoh, rho_n=rhon, hhl=hhl, r_air=rair,           &
!
                 tkv=vtyp(ivtype)%tkv, tsv=vtyp(ivtype)%tsv,          &
!
                 impl_weight=impl_weight,                             &
!
                 disc_mom=a(:,:,1), expl_mom=a(:,:,2),                &
                 impl_mom=a(:,:,3), invs_mom=a(:,:,4),                &
                 diff_dep=a(:,:,5), diff_mom=len_scale,               &
                 invs_fac=frh, scal_fac=frm,                          &
!
                 dif_tend=dicke, cur_prof=cur_prof, eff_flux=vari(:,:,m) )

            !Beachte:
            !'frh', 'frm' und 'len_scale' sind genauso wie 'a(:,:,1:5)' Hilfsspeicher in 'vert_grad_diff'.
            !Weil Fluesse ab "n>=liq=nmvar" nicht mehr benoetigt werden, bleibt 'vari' nur bis 
            ! 'nmvar' dimensioniert und vari(nmvar) wird auch fuer "n>nmvar" benutzt.

!           Sichern der Tendenzen:

            IF (n.EQ.tem) THEN
               DO k=1,ke
!DIR$ IVDEP
                  DO i=ivstart, ivend
                     dvar(n)%at(i,k)=dvar(n)%at(i,k)+exner(i,k)*dicke(i,k)*tinc(n)
                  END DO
               END DO
            ELSE
               DO k=1,ke
!DIR$ IVDEP
                  DO i=ivstart, ivend
                     dvar(n)%at(i,k)=dvar(n)%at(i,k)+dicke(i,k)*tinc(n)
                  END DO
               END DO
            END IF

            IF (n.EQ.vap .AND. PRESENT(qv_conv)) THEN
               !qv-flux-convergence (always a tendency) needs to be adapted:
               DO k=1,ke
                  IF (lqvcrst) THEN 
                     !by initializing 'qv_conv' with vertical qv-diffusion:
!DIR$ IVDEP
                     DO i=ivstart, ivend
                        qv_conv(i,k)=dicke(i,k)
                     END DO 
                  ELSE !by adding vertical qv-diffusion to 'qv_conv':
!DIR$ IVDEP
                     DO i=ivstart, ivend
                        qv_conv(i,k)=qv_conv(i,k)+dicke(i,k)
                     END DO
                  END IF
               END DO
            END IF       
                    
         END IF !diffusion calculation requested
         END DO !1, ndiff    

!-----------------------------------------------------------------

!Achtung:
!Ist cp-Fluss tatsaechlich der thermische Erdbodenantrieb?
!Was gilt im Falle der T-Gleichung in cv-Form?

!        Berechnung der effektiven Oberflaechenflussdichten:

!Achtung: "lscadif" ergaenzt
         IF (.NOT.(lsfluse .AND. lsflcnd) .AND. lscadif) THEN 
            !effektive Oberfl.flussdichten wurden neu bestimmt

            IF (PRESENT(shfl_s) .OR. lscm) THEN
!DIR$ IVDEP
               DO i=ivstart, ivend
                  shflu_s(i)=eprs(i,ke1)*cp_d*vari(i,ke1,tet)
               END DO
            END IF
            IF (PRESENT(qvfl_s) .OR. lscm) THEN
!DIR$ IVDEP
               DO i=ivstart, ivend
                  qvflu_s(i)=vari(i,ke1,vap)
               END DO
            END IF

!---------------------------------------------------------------------------------------
#ifdef SCLM
            IF (lsclm .AND. latmflu) THEN
               !Berechnung der Enthalpieflussdichten:

               SHF%mod(0)%val=shflu_s(im)     ; SHF%mod(0)%vst=i_cal
               LHF%mod(0)%val=qvflu_s(im)*lh_v; LHF%mod(0)%vst=i_cal

               !Note:
               !IF ".NOT.latmflu", SHF and LHF either are loaded by the fluxes used for
               ! the soil budget (lertflu) or they have been loaded above by the explicit 
               ! SHF and LHF at the surface (lsurflu).
               !SHF and LHF are positive downward and they may have been corrected with
               ! vertical integrated correction tendencies.
               !Thus they always refer to the used flux densities, which are only then equal
               ! to the explicit surface flux density, if a lower flux condition is used "lsflcnd=.TRUE.".
            END IF
#endif
!SCLM-----------------------------------------------------------------------------------

            !Bem: shflu_s und qvflu_s, sowie SHF und LHF sind positiv abwaerts!

         END IF

         IF (lum_dif .AND. PRESENT(umfl_s)) THEN
!DIR$ IVDEP
            DO i=ivstart, ivend
               umfl_s(i)=vari(i,ke1,u_m)
            END DO
         END IF
         IF (lvm_dif .AND. PRESENT(vmfl_s)) THEN
!DIR$ IVDEP
            DO i=ivstart, ivend
               vmfl_s(i)=vari(i,ke1,v_m)
            END DO
         END IF

!--------------------------------------------------
      END IF !Vertikaldiffusion wird hier berechnet
!--------------------------------------------------


! 11) Berechnung der Tendenzen infolge horizontaler turb. Diffusion
!     (und aufsummieren auf die Tendenzfelder):

! 16) Deallocierung lokaler dynamischer Felder:

    ! IF (.NOT.PRESENT(epr))   DEALLOCATE ( exner, STAT=ilocstat )
    ! IF (.NOT.PRESENT(c_big)) DEALLOCATE ( cbig,  STAT=ilocstat )
    ! IF (.NOT.PRESENT(c_sml)) DEALLOCATE ( csml,  STAT=ilocstat )
    ! IF (.NOT.PRESENT(r_air)) DEALLOCATE ( rair,  STAT=ilocstat )

END SUBROUTINE turbdiff

!==============================================================================

END MODULE turb_diffusion
