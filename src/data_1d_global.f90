!+ Data module for single column model
!------------------------------------------------------------------------------

MODULE data_1d_global

!------------------------------------------------------------------------------
!
! Description:
!
!   Definition of all global variables used for the single column mode,
!   especially the fundamental datastructure 'var_dat' and the model variables
!   that all have the TYPE 'var_dat'. If the model variable is present in the
!   linked general model, there is always a component of the structure pointing
!   to the respective REAL field There is a structure component for the error
!   variance of the variable as well. Further there are components for 
!   measurements and time series generated form the measurements, which can be
!   used to force the model. There are a couple of other components for special
!   variable properties. If there are other variables belonging to diagnostic
!   levels of the respective variable, there is also a component pointing to 
!   these variables related to diagnostic levels.
!
! Current Code Owner: DWD, Matthias Raschendorfer
!  phone:  +49  69  8062 2708
!  fax:    +49  69  8062 3721
!  email:  matthias.raschendorfer@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! V3_23        2007/03/30 Matthias Raschendorfer
!  Initial release
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_23        2012/05/10 Oliver Fuhrer
!  Removed obsolete Fortran features
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

USE data_parameters , ONLY :   &
    wp,         & ! kind-type parameter for "normal" integer variables
    iintegers     ! KIND-type parameters for real variables

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Configuration parameters:

INTEGER   (KIND=iintegers), PARAMETER :: &

 n_tlev=2,    & !always 2 timelevels are used for the single column run

 n_uni=10,    & !maximal number of additional units for a dimension

 textlen=120, & !maximal lenght of a text string
 wordlen=40,  & !maximal length of a word string

 if_loc=20,   & !number of the local opend file

! Field positions of the variable indices:

  p_sta=1,   & !variable status 
   i_nop=0,   & !not present
   i_est=1,   & !estimated value (e.g. from a previous time step)
   i_cal=2,   & !actualised value by a model calculation
   i_app=3,   & !approximate measurement by conversions or vertical interpolation
   i_dir=4,   & !directy measurement or assimilation including conversions of them
   i_def=5,   & !predefined value (like soil depth) that never must be changed 
  p_inp=2,   & !input index
   i_noi=0,   & !no input
   i_mes=1,   & !only measurement input
   i_int=2,   & !input also for variabel conversions during interpolation in time
   i_mod=3,   & !input also for model forcing
  p_out=3,   & !output index
  p_n_mes=4, & !number of measure points
  p_l_mes=5, & !measurement number just after model time
  p_l_cit=6, & !index of the latest conversion iteration

 n_ind=p_l_cit  !number of variable indices 

  !Note: ('i_nop', 'i_cal', 'i_dir') and ('i_noi', 'i_mes', 'i_int', 'i_mod')
  !      must have increasing values!
  !      'i_nop' has to be '0' and 'i_est' has to be '1'!

INTEGER   (KIND=iintegers), PARAMETER :: &

! Level positions:

  i_nol=0,    & !no level position

  i_mid=1,    & !mid layer position
  i_bnd=2,    & !layer boundary position

 n_ppos=i_bnd,   & !number of profile level positions

  i_low=0,    & !lowest model level
  i_sur=1,    & !surface level
  i_2m =2,    & ! 2m level
  i_10m=3,    & !10m level

 n_hpos=i_10m,   & !number of diagnostic atmosph. height positions

  i_850hPa=4, & !850hPa level

 n_dpos=i_850hPa   !number of diagnostic atmospheric positions

!Note: The diagn. atmo. height positions have to increase with increasing height.
!      The diagn. atmo. press. positions have to increase with decreasing press...

INTEGER   (KIND=iintegers), PARAMETER :: &

! Indicators of the effect relation between model and variable:

  i_noe=0, & !variables with (so far) no direct effect from the model as they are
             !no direct model variables. They may be assimilated anywhere
  i_lev=1, & !level variable that has not to be assimilated
             !(like soil layer depth or atm. level pressure)
  i_exe=2, & !variable with explicite effect at the beginning of the time loop
             !(has to be assimilated at the beginning of the time loop,
             ! but is not used for the initial model state, thus will be init. by the model)
  i_ini=3, & !variable with initial effect (used for the initial model state
             !and thus has an explicite effect as well)
  i_ext=4, & !external parameter (the same as initial effect,
             !but will not be updated by the model itself) 
  i_fix=5, & !fix parameter (the same as external parameter,
             !but with fixed values during the model run even in reality)
  i_tra=6, & !variable that is calculated within the transfer scheme
             !and has to be assimilated after it
  i_tur=7, & !variable that is calculated within the diffusion scheme
             !and has to be assimilated after it
  i_soi=8, & !variable that is calculated within the soil model  
             !and has to be assimilated after it
  i_phy=9, & !variable that is updated elsewhere within the physical routines
             !and has to be updated after the model physics
  i_dia=10,& !variable that is calculated within other diagnostic calculations
             !and has to be assimilated after the model diagnostics
  i_ten=11,& !tendency field, that has to be set to zero at the
             !beginning of the time loop. It has to be assimilated
             !before the general time step integration.
  i_cor=12,& !correction field, which has to be differentiated and then 
             !integrated for each time step

 n_meff=i_cor     !number of model effect indicators

INTEGER   (KIND=iintegers), PARAMETER :: &

! Indicators for variable dimensions:

  i_Leng=1 ,& !length (constant in time)
  i_Heig=2 ,& !height (of atm. model layers)
  i_SDep=3 ,& !soil depth (constant in time)
  i_VMag=4 ,& !magnitude of horizontal velocity
  i_VCmp=5 ,& !compount of the horizontal velocity vector
  i_Temp=6 ,& !temperature
  i_TDew=7 ,& !dew point temperature
  i_SpHu=8, & !specific humidity of the air
  i_SpWt=9 ,& !specific cloud water
  i_SpIc=10,& !specific cloud ice
  i_ReHu=11,& !relative humidity
  i_AtPr=12,& !atm. pressure
  i_O3Pr=13,& !O3 -  prat. pressure
  i_O3Lv=14,& !pressure level of the O3-maximum
  i_SoWt=15,& !soil water
  i_SoIc=16,& !soil ice
  i_Frac=17,& !fraction
  i_LAId=18,& !leave area index
  i_GAng=19,& !geographic angles
  i_WiDr=20,& !wind direction
  i_Numb=21,& !number
  i_Freq=22,& !frequency
  i_Dens=23,& !density
  i_IWtL=24,& !interception water level
  i_SWtl=25,& !snow water level
  i_CDif=26,& !diffusion coefficient
  i_LRog=27,& !roughness length
  i_LRot=28,& !root depth
  i_LAlt=29,& !Altitude of the surface above NN
  i_EFlu=30,& !energy flux density
  i_VVar=31,& !velocity variance
  i_TVar=32,& !velocity variance
  i_TpVl=33,& !temperat. veloc. covariance

 n_dim=i_TpVl  !number of physical variable dimensions

INTEGER   (KIND=iintegers), PARAMETER :: &

! Indicators of the conversion groups:

  i_wind=1,& !wind group
  i_ahum=2,& !atmospheric humidity group
  i_dhum=3,& !dry soil surface humidity group
  i_snow=4,& !snow group
  i_intc=5,& !interception water group
  i_csrf=6,& !surface concentration group
  i_vtem=7,& !virtual temperature group
  i_ptem=8,& !potential temperature group
  i_glaw=9,& !gas law group

 n_group=i_glaw, & !number of conversion groups

! Indices for the type of measurement storage:

  p_mtim=1,  & !time of measurement
  p_mval=2,  & !value of measurement
  p_koef=3,  & !interpolation coefficient

 n_mtyp=p_koef     !'n_mtyp' stored information for a measurement point

INTEGER   (KIND=iintegers), PARAMETER :: &

!Profile classes:

  i_sd=-1,&  !index for the soil depth profile
  i_sf=0, &  !index for the surface
  i_ap=1, &  !index for the measured atm. press. profile 
  i_ah=2, &  !index for the measured atm. height profile

 m_spr=IABS(i_sd), & !number of all   soil profile classes
 n_spr=m_spr,      & !number of model soil profile classes

 m_apr=i_ap, & !number of model atm. profile classes
 n_apr=i_ah, & !number of all   atm. profile classes

!Note: i_ap is the only model atm. profile class,
!      thus it has to have the index '1'!

!Variable types:

  i_l=0, & !level variables

  i_p=1, & !prognostic direct model variables (with dimens. for horiz. pos. and time)
  i_d=2, & !diagnostic direct model variables (with dimens. for horiz. pos. and not for time)
  i_g=3, & !global     direct model variables (with dimens. neither for horiz. pos. nor for time)

  i_s=4, & !special measurement variables
           !(no model variables, with dimens. neither for horiz. pos. nor for time)
  i_c=5, & !correction integral of prognostic variables

 m_typ=i_g, & !number of model variable types
 r_typ=i_s, & !number of all variable types
 a_typ=i_c    !number of all variable types including the correction integrals

!Note: "i_l=0" is obligatory!

INTEGER   (KIND=iintegers), PARAMETER :: &

!Index for the unit groups:

  u_num=0, & !for the number of used units
  u_int=1, & !for internal units (SI-units, proper for vert. interpolation)
  u_mod=2, & !for special model units (like 'int'-units but with some exeptions)
  u_out=3, & !for output units

 n_ung=u_out, & !number of unit groups

!Incex for the adjutment periods:

   i_loc=0, & !index for the local model in space
   i_day=1, & !index for the dayly period
   i_syn=2, & !index for the synoptical period
   i_ann=3, & !index for the annual period
   i_kli=4, & !index for the klimatological period

 n_cyc=i_kli, & !number of adjustment cycles (day, synoptic period, year)

!Index for the combination mode:

   c_exp=1, & !explicite combination of model values and measurements
   c_imp=2, & !implicite combination of model values and measurements

!Index for water or ice:

   i_wat=1, & !water
   i_ice=2, & !ice

!Index for code positions and purposes (aims) of code segments:

   p_meas_init=1, & !at initialisation of measurement input
   p_meas_prof=2, & !at input of measured vertical profiles
   p_dset_init=3, & !at initialisation of a input data set
   p_loop_init=4, & !at initialisation of the time loop
   p_step_init=5, & !at initialisation of a time step
   p_after_phys=6,& !after the physic package
   p_during_inp=7,& !during input (only for aim of 'variable_conversion')
   p_during_out=8,& !during output  ,,   ,,  ,, ,,          ,,

   n_code_posi=p_during_out

CHARACTER (LEN=textlen), ALLOCATABLE :: &

   posi_name(:) !names of the code positions

! Numbers:

REAL  (KIND=wp),     PARAMETER :: &

  z0=0.0_wp,   &
  z1=1.0_wp,   &
  z2=2.0_wp,   &
  z3=3.0_wp,   &
  z4=4.0_wp,   &
  z5=5.0_wp,   &
  z10=10.0_wp, &
  z1d2=z1/z2,  &
  z1d3=z1/z3,  &
  z1d4=z1/z4,  &
  z1d5=z1/z5,  &
  zp1=0.1_wp, zp2=0.2_wp, zp3=0.3_wp, zp4=0.4_wp, zp5=0.5_wp, &
  zp6=0.6_wp, zp7=0.7_wp, zp8=0.8_wp, zp9=0.9_wp, &

  kilo=1000.0_wp, &
  hect=100.0_wp,  &
  mili=0.001_wp,  &
  cent=0.01_wp

! Tuning parameter:

REAL  (KIND=wp)     :: &

  bld_hig=1000.0_wp,  & !blending height [m]
  geo_hig=5000.0_wp,  & !height above ground of the geostrophic wind in [m]

  ef_cmb=2.0_wp,      & !error factor in the combination formula
  ef_mod=0.1_wp,      & !error factor for the model calculations
  af_mes=0.01_wp,     & !amplitude factor for the registration error of measurements
  af_loc=1.0_wp         !amplitude factor for uncertainties in location

INTEGER   (KIND=iintegers) :: &

  cnv_it=3    !maximal number of conversion iterations

!Note:
!Most of these tuning parameters may be overwritten by the input of the 'para' NAMELIST.  

! Physical constants and parameters:

REAL  (KIND=wp),     POINTER :: & !taken form the used model

  T0,      & !zero point temperature [C]
  T0_sea,  & !freezing temperature of ocean water [C]
  Cp_d,    & !specific heat capacity of dry air [J/(Kg*K)]
  R_d,     & !gas constant for dry air [J/(Kg*K)]
  R_v,     & !gas constant for water vapor [J/(Kg*K)]
  RdoCp,   & !R_d/Cp_d
  Rdv,     & !R_d/R_v
  o_m_Rdv, & !1.0-R_d/R_v
  Rvd_m_o, & !R_v/R_d-1.0

  p1_e,    & !p_sat(T1_e) [Pa]
  T1_e,    & !reference temp. for evap. [C]

  g,       & !gravity

  FCap(:), & !field capacity of the soil types
  PVol(:)    !pore volume of the soil types

 REAL (KIND=wp)     :: &

! Special constants for evaporation of water and ice:

  Lh_e(2), & !latent heat of vapourization [J/Kg]

  T2_e(2), & !scaling temp. for evap. over water and ice [C]
  fc_e(2), & !factor for evap. over water and ice [1]
  dT_e(2), & !temp.const. in d(q_sat)/dT fc_e*(T1_e-T2_e) 

! Mathematical constants:

  pi         !circle number

! Special parameters:

REAL (KIND=wp),     PARAMETER :: &

  ef_sin=0.5_wp,           & !error factor for the sinus wave

  day_hor=24.0_wp,         & !hours of a day 
  hor_sec=3600.0_wp,       & !seconds of an hour
  day_sec=day_hor*hor_sec, & !seconds of a day
  syn_day=5.0_wp,          & !days of a synoptical period
  yar_day=365.0_wp,        & !days of a year

  sur_prs=100000.0_wp,     & !reference surface pressure in [Pa]

  abl_prs_frc=  0.85_wp,   & !pressure at the ABL-height  / surface pressur
  trp_prs_frc=  0.25_wp,   & !   ,,   ,,  ,,  tropopause  /   ,,      ,,
  stp_prs_frc= 0.001_wp,   & !   ,,   ,,  ,,  stratopause /   ,,      ,,
  mep_prs_frc=1.0E-5_wp,   & !   ,,   ,,  ,,  mesopause   /   ,,      ,,

  rog_hig=0.1_wp,          & !height of the roughn. layer in [m] a.G.
  abl_hig=1000.0_wp,       & !height of the ABL           in [m] a.G.
  trp_hig=10000.0_wp,      & !height of the tropopause    in [m] a.G.
  stp_hig=50000.0_wp,      & !height of the stratopause   in [m] a.G.
  mep_hig=80000.0_wp,      & !height of the mesopause     in [m] a.G.
 
  abl_tmp=  0.0_wp,        & !mean temperature at the ABL-height  in [C]
  trp_tmp=-50.0_wp,        & !mean temperature at the tropopause  in [C]
  stp_tmp=  0.0_wp,        & !mean temperature at the stratopause in [C]
  mep_tmp=-80.0_wp           !mean temperature at the stratopause in [C]

REAL (KIND=wp)     :: &

  !Adjustment periods and heights:

  dt_cyc(i_loc:n_cyc)=(/z0,      &       !zero period
                        z1,      &       !dayly period in days
                        syn_day, &       !synoptic period in days
                        yar_day, &       !annual period in days
                        3*yar_day/), &   !klimat. period in days

  hi_cyc(i_loc:n_cyc)=(/z0,      &       !zero height
                        abl_hig, &       !adjustm. height for the dayly period
                        trp_hig, &       !adjustm. height for the synptic period
                        stp_hig, &       !adjustm. height for the annual period 
                        mep_hig/)        !adjustm. height for the klimat. period

!Note:
!The adjustment periods and heights for the locale space mode will be defined later!
!The units used for all constants and parameters not necessarily are 'int'-units,
!but they belong to the (fixed) chose of 'int'-units!

! global variables:

LOGICAL :: &

  lsclm=.FALSE., & !SCLM is not yet running 

  lhlev=.FALSE., & !SCLM refers to height levels

  lscreen,  & !output on screen active
  ldetail,  & !detailed control output
  lwindpr,  & !geostrophic wind forcing active
  lsatadj,  & !saturation adjustment for condensation and evaporation activ
  lradtnd,  & !radiative tendencies in the atmosphere are calculated
  lpotloc,  & !using a local potential temperature with respect to the surface pressure
  lisobar,  & !calculating isobaric tendencies even for height based model levels
  lctl,     & !information in the file 'control' are plotted for the actual time step
  lplo,     & !model outpuit is plotted for the actual time step
  lmlsoil,  & !multi layer soil model active

  lprog,    & !a pure model run
  lforc,    & !a forced model run
  lassi,    & !an assimilation run
  lmeas       !a measurement run

CHARACTER (LEN=wordlen) :: &

  fil_suff, & !suffix for output files
  place       !name of the location

CHARACTER (LEN=4) :: &
  run_mod     !modus of the model run

INTEGER   (KIND=iintegers), POINTER :: &

  n_old, n_now, n_new, n_tke, & !time step indices
  n_tstep                       !nummber of time steps

INTEGER   (KIND=iintegers) :: &

  n_start,  & !time steps belonging to the delay of the model start with
              !respect to the previous full hour
  end_step, & !last model time step'
  plo_step, & !step interval of recording the results, if 'plo_hour=0'
  ctl_step, & !step interval of recording additional information, if 'ctl_hour=0'

  k_geo,    & !level number of the geostrophic wind
  n_hyd,    & !number of the hydrological soil levels (lmulti_layer=F)
  n_efr,    & !number of thermal efr levels (lmulti_layer=F)

  if_num,   & !highest actual used file number
  if_cnt,   & !filenumber of the control file

  n_hori,   & !number of grid points in a horizontal direction
              !depends on the configuration of the used 3d-model (LM)
  im, jm,   & !horizontal indices relevant for the 1d-model
  ke_erd,   & !number of the soil levels
  ke_atm,   & !number of the atmospheric levels

  m_var,    & !maximal number of variable names for a given
              !profile class and variable type

  n_lpos,   & !number of all dignostic level positions

  u_act       !unit index for the actual used units in the respective programm section

REAL  (KIND=wp),     TARGET :: &

  dtsec,    & !timestep for the physics [sec]
  dthor,    & !timestep for the physics [hor]
  nbegin,   & !startstep of recording the results after model start
  mod_tim,  & !actual model run time in [hr] after 'SDat'
  mes_tim,  & !actual measurem. time in [hr] after 'SDat'

  res_hori, & !horizontal resolution [DEG]
  PTop,     & !reference pressure at the top terrain following levle [Pa]

  soil_diff   !depth of the temperature wave into the soil in [m]
              !devided by the SQRT of the time period in [s]

REAL (KIND=wp),     ALLOCATABLE, TARGET :: &

  zhyd(:),  & !main level height of the special hydrological soil layers [m] (lbodwgl=F)
  dzhyd(:), & !depth of the special hydrological soil layers [m] (lbodwgl=F)

  zefr(:), &  !boundary level depth for the two efr soil layers [m] (lbodwgl=F)

  Ugeo(:), & !geostrophic wind in longitudinal direction
  Vgeo(:)    !geostrophic wind in latitudinal direction

INTEGER (KIND=iintegers), ALLOCATABLE :: &

  n_nam(:,:), & !number of variable names for a given profiles and type

  i_hyd(:),   & !indices of the diagnostic hydrological levels
  i_efr(:)      !indices of the diagnostic efr levels

INTEGER (KIND=iintegers), TARGET :: &

  unid(u_num:n_ung,n_dim)=-1 !unit index for all dimensions and unit groups

! Data structures for the model variables:

TYPE tot_val !total value
     SEQUENCE
     REAL (KIND=wp)     :: v !variable value
     REAL (KIND=wp)     :: e !variable error (variance)
END TYPE tot_val

TYPE tot_pnt !total pointer
     SEQUENCE
     REAL (KIND=wp),     POINTER :: v !variable value
     REAL (KIND=wp),     POINTER :: e !variable error (variance)
END TYPE tot_pnt

TYPE dat_time !date-time group
     CHARACTER (LEN=8) :: day
     REAL (KIND=wp)    :: hor
END TYPE dat_time

TYPE (dat_time) :: &

     SDat, & !start time
     IDat, & !initial time (belongs to the previous full hour of SDat)
     ADat    !actual time

TYPE mod_dat !variable data on model levels

     REAL (KIND=wp),           POINTER :: lev      !level value
     REAL (KIND=wp),           POINTER :: val      !variable value
     REAL (KIND=wp),           POINTER :: err      !variable error variance (belonging to internal units)
     INTEGER (KIND=iintegers), POINTER :: ind(:)   !index vector with various information
     INTEGER (KIND=iintegers), POINTER :: vst      !index for the model  variable status
     INTEGER (KIND=iintegers), POINTER :: bst      !index for the buffer variable status
     LOGICAL                           :: lnk      !true if a diagn. level linked to an other model variable
     TYPE (tot_pnt)                    :: tot      !total value, that is the tuple of value and error variance
     TYPE (tot_pnt)                    :: buf      !variable buffer for the total value

END TYPE mod_dat

TYPE mes_dat !measurement data of model variables

     INTEGER (KIND=iintegers), POINTER :: lev(:)   !level index profile for measurements
     REAL (KIND=wp),           POINTER :: val(:)   !measurement values on measurement levles 
     REAL (KIND=wp),           POINTER :: err(:)   !measurement error on measurement levles 
     INTEGER (KIND=iintegers), POINTER :: flg(:)   !read flag profile for measurements

END TYPE mes_dat

TYPE ser_dat !time serie of variable data
     REAL (KIND=wp),           POINTER :: tim(:)   !serie of the time values
     REAL (KIND=wp),           POINTER :: val(:)   !serie of the variable values
     REAL (KIND=wp),           POINTER :: err(:)   !serie of the variable errors
     REAL (KIND=wp),           POINTER :: kof(:)   !serie of the interpolation koefficients

END TYPE ser_dat

TYPE var_prop !variable property
     
     CHARACTER (LEN=wordlen)  :: na          !unit string
     REAL    (KIND=wp)        :: am          !maximal amplitude of the time serie
     REAL    (KIND=wp)        :: af(0:n_cyc) !amlitude factor for the adjustment cycle
                                             !(including the spacial mode)
END TYPE var_prop

TYPE uni_conv !parameter for unit conversions

     CHARACTER (LEN=wordlen)  :: nam       !unit string 
     REAL      (KIND=wp)      :: add       !additive conversion parameter into interp. units
     REAL      (KIND=wp)      :: mul       !multipl. conversion parameter into interp. units
     REAL      (KIND=wp)      :: pot       !potentl. conversion parameter into interp. units

END TYPE uni_conv

TYPE var_dat !all data of a model variable

     !Variable identification:

     CHARACTER (LEN=wordlen),  POINTER :: nam      !name(s) of the variable
     INTEGER (KIND=iintegers), POINTER :: dim(:)   !indicator for the physical dimensions
     INTEGER (KIND=iintegers), POINTER :: pos      !indicator for the level position 
     INTEGER (KIND=iintegers), POINTER :: low      !index for the lowest level position (closest to the surface)
     INTEGER (KIND=iintegers), POINTER :: hig      !index for the highest level position (farest from the surface)
     INTEGER (KIND=iintegers), POINTER :: stp      !index step to the next farer level from the sruface
     INTEGER (KIND=iintegers), POINTER :: dia(:)   !index for the (additional) diagnostic level
     INTEGER (KIND=iintegers), POINTER :: top      !index for the (additional) top level
     INTEGER (KIND=iintegers), POINTER :: act      !actual time step index
     INTEGER (KIND=iintegers), POINTER :: eff      !index for the effect on the model

     INTEGER (KIND=iintegers), POINTER :: cls      !variable class index
     INTEGER (KIND=iintegers), POINTER :: typ      !variable type index
     INTEGER (KIND=iintegers), POINTER :: ind      !variable indication index

     !Model variable value fields:

     REAL (KIND=wp),     POINTER :: lv0            !level value of a surface variable
     REAL (KIND=wp),     POINTER :: lv1(:)         !profile of the vertical level coordinate

     REAL (KIND=wp),     POINTER :: v0d            !0-d value
     REAL (KIND=wp),     POINTER :: v1d(:)         !1-d value field
     REAL (KIND=wp),     POINTER :: v2d(:,:)       !2-d value field
     REAL (KIND=wp),     POINTER :: v3d(:,:,:)     !3-d value field 
     REAL (KIND=wp),     POINTER :: v4d(:,:,:,:)   !4-d value field

     REAL (KIND=wp),     POINTER :: amp(:,:)       !amplitude of the time serie dep. on adjustment cycle and level

     !Variable value structers:

     TYPE (sing_var),  POINTER :: dvr(:)           !diagnostic variables

     TYPE (mod_dat), POINTER :: mod(:)             !variable data on model levels
     TYPE (ser_dat), POINTER :: ser(:)             !time serie of variable data on model levels
     TYPE (mes_dat), POINTER :: mes                !measurement data of model variables

END TYPE var_dat

TYPE sing_var !single variable

     TYPE (var_dat), POINTER :: var !model variable

END TYPE sing_var

TYPE var_list !variabel list

     TYPE (var_dat), POINTER :: var(:) !vector of variables

END TYPE var_list

TYPE conv_var !conversion variable

     TYPE (var_dat),           POINTER :: var    !model variable
     INTEGER (KIND=iintegers), POINTER :: nit(:) !number of conversion iterations for each level

END  TYPE conv_var
     

! Surface variables:
!-------------------------------------

INTEGER   (KIND=iintegers), PARAMETER :: &

 n_psrf=11 !maximal number of prognostic surface variables

TYPE (var_dat), POINTER :: &

 PSM, & !mean surface pressure
 TSE, & !surface temperature of the earth (sonw- and ice free part of the land)
 TSI, & !surface temperature of ice and snow
 TSM, & !mean temperature at the surface
 TEM, & !mean temperature of the soil
 QVSM,& !mean spec. humid. at the surface
 WSI, & !surface water of the interception storage
 WSS, & !surface water of the snow storage
 WE1, & !water content of the 1-st special hydr. soil layer
 WE2, & !water content of the 2-nd special hydr. soil layer
 WE3    !water content of the 3-st special hydr. soil layer

INTEGER   (KIND=iintegers), PARAMETER :: &

 n_dsrf=25 !maximal number of diagnostic surface variables

TYPE (var_dat), POINTER :: &

 GLA, & !geographic longitude
 GLO, & !geographic longitude
 HSN, & !height of the lower boundary above NN
 FRL, & !land fraction
 PLC, & !plant cover
 LAI, & !leave area index
 SOT, & !leave area index
 ROD, & !root depth
 RLM, & !roughness length for momentum
 O3I, & !vertical integrated ozone content
 O3M, & !pressure hight of ozone maximum
 TEC, & !climatological soil temperature
 WEC, & !climatological soil water content
 TA2M, & !synoptical temperature       2m above the surface
 QVA2M,& !synoptical specivic humidity 2m above the surface
 TDA2M,& !synoptical dew point temper. 2m above the surface
 VZ10M,& !synoptical zonal      velocity 10m above the surface
 VM10M,& !synoptical meridional velocity 10m above the surface
 RADSOS,&!solar radiation flux density at the surface
 RADTHS,&!ther. radiation flux density at the surface
 RADSOT,&!solar radiation flux density at the model top
 RADTHT,&!ther. radiation flux density at the model top
 TMIN2M,&!minimum of TA2M
 TMAX2M,&!maximum of TA2M
 VMAX10M !maximum of wind magnitude at 10m above the surface

INTEGER   (KIND=iintegers), PARAMETER :: &

 n_csrf=n_psrf !maximal number of correction integrals for the progn. surf. var.

TYPE (sing_var) :: &

 srf_cor(n_csrf) !correction integrals for the progn. surf. var.

INTEGER   (KIND=iintegers), PARAMETER :: &

 n_gsrf=0 !maximal number of global surface measurement variables

INTEGER   (KIND=iintegers), PARAMETER :: &

 n_ssrf=16 !maximal number of special surface measurement variables

TYPE (var_dat), POINTER :: &

 INC,   & !interception water cover of the snow free part of the land
 SNC,   & !snow cover of the land
 SIC,   & !sea ice cover
 TSS,   & !surface temperature of the sea part
 TSL,   & !surface temperature of the land part
 QVSL,  & !spec. humid. over the land part
 QVSD,  & !spec. humid. over the dry soil
 TDSD,  & !dew point temperature at the dry soil surface
 RHSD,  & !dew point temperature at the dry soil surface
 TDSM,  & !mean dew point temperature at the surface
 DSM,   & !mean atmosperic density at the surface
 P2M,   & !pressure  2m above the surface
 P10M,  & !pressure 10m above the surface
 WST10M,& !windstr. 10m above the surface
 SHF,   & !sensible heat flux
 LHF,   & !latent heat flux
 LMO      !Monin-Obuchov stability length

! Soil variable profiles:
!-------------------------------------

INTEGER   (KIND=iintegers), PARAMETER :: &

 n_perd=3 !maximal number of prognostic soil variables

TYPE (var_dat), POINTER :: &

 TOE, & !ordinary temperature prof. of the earth (of the multi layer model)
 WCE, & !soil water content (of the multi layer model)
 ICE    !soil ice profile (of the multi layer model)

INTEGER   (KIND=iintegers), PARAMETER :: &

 n_derd=0 !maximal number of diagnostic soil variables

INTEGER   (KIND=iintegers), PARAMETER :: &

 n_cerd=n_perd !maximal number of correction integrals for the progn. soil. var.

TYPE (sing_var) :: &

 erd_cor(n_cerd) !correction integrals for the progn. soil. var.

INTEGER   (KIND=iintegers), PARAMETER :: &

 n_gerd=3 !maximal number of global soil variables

TYPE (var_dat), POINTER :: &

 HME, & !height of the soil mid levels
 HBE, & !height of the soil boundary levels
 HDE    !height difference of the soil layer boundaries

INTEGER   (KIND=iintegers), PARAMETER :: &

 n_serd=0 !maximal number of special soil variables

! Atmospheric variable profiles:
!-------------------------------------

INTEGER   (KIND=iintegers), PARAMETER :: &

 n_patm=7  !maximal number of prognostic atmospheric variables

TYPE (var_dat), POINTER :: &

 VZA, & !zonal velocity of the air
 VMA, & !meridional velocity of the air
 TOA, & !ordinary temperature of the air
 QVA, & !specific humidity of the air
 QWA, & !specific cloud water content of the air
 QIA, & !specific cloud ice content of the air
 TKE    !turbulent kinetic energy

INTEGER   (KIND=iintegers), PARAMETER :: &

 n_datm=11+n_patm !maximal number of diagnostic atmospheric variables

TYPE (var_dat), POINTER :: &

 HMA, & !height of the atmospheric mid - or measurement levels
 HBA, & !height of the atmospheric boundary levels
 PMA, & !pressure of the atmospheric mid - or measurement levels 
 PBA, & !pressure of the atmospheric boundary levels 
 PDA, & !pressure difference within a atmothpheric layer
 DMA, & !density of the atmosphere
 TDM, & !turbulent diffusion coefficient for momentum
 TDH, & !turbulent diffusion coefficient for scalars (heat)
 CLC, & !cloud cover
 SRHA,& !solar radiative heating of the atmosphere
 TRHA   !thermal radiative heating of the atmosphere

TYPE (sing_var) :: &

 atm_cnv(n_patm)    !convective tendencies for the progn. atmo. var.

INTEGER   (KIND=iintegers), PARAMETER :: &

 n_catm=n_patm  !maximal number of correction integrals for the progn. atmo. var.

TYPE (sing_var) :: &

 atm_cor(n_catm) !correction integrals for the progn. atmo. var.

INTEGER   (KIND=iintegers), PARAMETER :: &

 n_gatm=0    !maximal number of global atmospheric measurement variables

INTEGER   (KIND=iintegers), PARAMETER :: &

 n_satm=18   !maximal number of special atmospheric measurement variables

TYPE (var_dat), POINTER :: &

  WST, & !wind strength
  WDR, & !wind direction
  TDA, & !dew point temperature of the air
  TVA, & !virtual   temperature of the air
  TPA, & !potential temperature of the air
  RHA, & !relative humidity of the air
  UUA, & !zonal      velocity variance (u'u')
  VVA, & !meridional velocity variance (v'v')
  WWA, & !vertical   velocity variance (w'w')
  UWA, & !zonal-vertical      velocity covar (u'w')
  VWA, & !meridional-vertical velocity covar (v'w')
  UST, & !friction velocity (u_star)
  TWA, & !covariance of liq. wat. pot. temp and vert. vel. (tet_l'w')
  TTA, & !variance of liq. wat. pot. temp (tet_l'tet_l')
  BOYPR, & !buoyant production of TKE
  SHRPR, & !shear production of TKE
  DISSI, & !dissipation of TKE
  TRANP    !total transport of TKE

TYPE (sing_var), POINTER ::  &

 Lvar(:,:), & !level variabale for a given profile class and level position

 WErd(:),&    !water content of the special hydr. soil layers
 TErd(:)      !temperature at the efr soil levels

TYPE conv_group !structure for a conversion group
     INTEGER (KIND=iintegers)  :: num !number of group members
     INTEGER (KIND=iintegers)  :: red !reduced number of group members (without pressure variables)
     TYPE (conv_var),  POINTER :: mem(:) !vector for the group members
END TYPE conv_group

TYPE (conv_group) :: group(n_group) !conversion groups

TYPE (var_list) :: dom(0:a_typ, -n_spr:n_apr) !model domain structure dependend on 
                                              !variable class and type

TYPE (var_prop) :: prop(n_dim) !variable properties for all dimensions

TYPE (uni_conv) :: &

  conv(0:n_uni)=uni_conv('', z0, z0, z0) !unit conversions for all units of a dimension

END MODULE data_1d_global
