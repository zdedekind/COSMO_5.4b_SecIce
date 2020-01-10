!
!+ 3DVAR/COSMO feedback file tables and table entry values
!
! $Id: mo_fdbk_tables.f90,v 4.29 2013-10-04 08:51:39 for0adm Exp $
!-------------------------------------------------------------------------------
!
MODULE mo_fdbk_tables
!
!-------------------------------------------------------------------------------
! Description:
!   3DVAR/COSMO feedback file tables and table entry values,
!   specifying the valid values of the variables in the feedback file.
!   This module is used commonly by the COSMO model and 3DVAR program packages !
!
! Current Code Owners:
!    For DWD 3DVAR:                        For COSMO:
!    DWD, Andreas Rhodin                   DWD, Christoph Schraff
!    phone: +49 69 8062 2722               phone: +49 69 8062 2725
!    fax:   +49 69 8062 3721               fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de          email: christoph.schraff@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V4_22        2012/01/31 Christoph Schraff
!  Initial release, based on 3DVAR version V1_10, with further updates.
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer
!  Splitted a declaration with too many continuation lines
! V4_26        2012/12/06 Andreas Messer
!  Modifications for using also RTTOV10 and update with 3DVAR-V1_20:
!  3DVAR-V1_13 (this 3DVAR version was included in COSMO 4.22), except:
!               new codetype OC_WP_JP = 134 (Japanese wind profiler);
!               new variable 'n_vn' (number of variable numbers).
!  3DVAR-V1_15: New flag "FL_OPERATOR" (observation operator not applicable).
!  3DVAR-V1_17: New entries: ST_DISMISS, VN_REFL, VN_RADIANCE, VN_RREFL.
!  3DVAR-V1_19: New veri run type VT_LIN_ANA: linear operator on analysis (Y^a).
!  3DVAR-V1_20: New flags FL_FG_LBC, VN_CTH, VN_TRH.
! V4_28        2013/07/12 Christoph Schraff / Andreas Rhodin
!  'OC_WP_JP' made public, and updates to correspond with 3DVAR-V1_23:
!  3DVAR-V1_22: New variable number VN_VGUST (vertical gust);
!               mnemonic fixed for VN_CTH (CTH).
!  3DVAR-V1_23: New variable 'ct_nwc' (Cloud Type according to NWC-SAF) and
!               related type 'ct_nwc_entries' and parameters 'CT_*', 'OC_CLPRD'.
! V5_1         2014-11-28 Christoph Schraff
!  Update to 3DVAR-V1_31:
!  3DVAR-V1_27: New entry VE_BIASCOR for variational bias correction.
!  3DVAR-V1_29: Constants defined for wind/solar power observation operator:
!               OT_POWER, OC_PWIND, OC_PWSOL, VN_PWIND, VN_PWSOL.
!  3DVAR-V1_31: OC_ACARS added; OC_WP_JP made public.
! V5_3         2015-10-09 Christoph Schraff
!  - 3DVAR-V1_35: obs code types 'OC_MODES', 'OC_GPSRO', 'OC_GPSGB' and variable
!                 'VN_ELEV' added.
!  - Variables 'VN_RAD_GL', 'VN_RAD_DF', and 'VN_RAD_LW' added.
!  - Observation operator flags 'OF_MISSING', 'OF_RAD_CLEAR_SKY' etc added.
!
! CAUTION: This module is used commonly by the 3DVAR and COSMO main programs.!!!
!!!        Therefore, anybody wanting to introduce a modification to this    !!!
!!!        module in the context of either of these programs must consult    !!!
!!!        the 'current code owner' of this module for the other program,    !!!
!!!        in order to allow for checking that the modification will comply  !!!
!!!        with both program packages. This must be done before the          !!!
!!!        modification is put into the Version Control System (VCS).        !!!
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  original source to be used in 3DVAR/LETKF/COSMO
!------------------------------------------------------------------------------
!-------------
! Modules used
!-------------
use mo_t_table, only :      t_table,   &! derived type definition for tables
                       e => t_entry,   &!   and table entries
                            init_table  ! initialise a table
implicit none

!================
! public entities
!================
private
!-------
! Tables
!-------
public :: status            ! report and observation status flags
public :: flags             ! report and observation check  flags
public :: obstype           ! observation types
public :: codetype          ! code types
public :: varno             ! variables
public :: runtype           ! run type
public :: runclass          ! run class
public :: satsens           ! satellite sensors
public :: rsondtype         ! radiosonde type
public :: trackteqn         ! tracking technique
public :: meas_equip        ! measureing equipment used
public :: radiation_corr    ! solar & infrared radiation correction
public :: surftype          ! surface type
public :: surf_char         ! model surface characteristics
public :: ct_nwc            ! cloud type according to NWC SAF
public :: flg_1dvar         ! 1dvar processing flags
public :: flg_cld           ! 1dvar cloud flags
public :: level_sig         ! vertical level significance
public :: phase             ! aircraft phase, radiances fov, GPSRO PCD
public :: rollangle         ! aircraft roll angle
public :: retrtype          ! AMV retrieval type
public :: ensmem            ! generalised ensemble member flag
public :: oflag             ! opservation operator processing flag

!------------
! subroutines
!------------
public :: init_fdbk_tables  ! initialise the table (fill in table entries)
public :: clean_fdbk_tables ! deallocate the table

!----------
! run class
!----------
public :: RC_HAUPT, RC_VOR, RC_ASS, RC_TEST

!-----------------------------------------------
! constants: observation or report status values
!-----------------------------------------------
public :: ST_ACCEPTED, ST_ACTIVE, ST_MERGED, ST_PASSIVE, ST_REJECTED, &
          ST_PAS_REJ, ST_OBS_ONLY, ST_DISMISS

!----------------------------------
! observation or report check flags
!----------------------------------
public :: FL_OBSTYPE, FL_BLACKLIST, FL_SUSP_LOCT, FL_TIME, FL_AREA,      &
          FL_HEIGHT, FL_SURF, FL_CLOUD, FL_PRACTICE, FL_DATASET,         &
          FL_REDUNDANT, FL_FLIGHTTRACK, FL_MERGE, FL_THIN, FL_RULE,      &
          FL_OBS_ERR, FL_GROSS, FL_NO_BIASCOR, FL_FG, FL_FG_LBC, FL_NONE,&
          FL_NO_OBS, FL_OPERATOR

!------------------
! observation types
!------------------
public :: OT_SYNOP, OT_AIREP, OT_SATOB, OT_DRIBU, OT_TEMP,  OT_PILOT, &
          OT_SATEM, OT_PAOB,  OT_SCATT, OT_RAD,   OT_GPSRO, OT_GPSGB, &
          OT_RADAR, OT_POWER, n_ot

!-----------------------
! observation code types
!-----------------------
public :: OC_SRSCD, OC_ATSCD, OC_AHSCD, OC_ATSHS, OC_AIRCD, OC_CODAR, &
          OC_AMDAR, OC_STBCD, OC_DRBCD, OC_TESAC, OC_LDTCD, OC_SHTCD, &
          OC_TDROP, OC_LDPCD, OC_SHPCD, OC_ATOVS, OC_WP_EU, OC_RA_EU, &
          OC_PR_US, OC_RAVAD, OC_SEVIR, OC_ASCAT, OC_QSCAT, OC_TMPMB, &
          OC_PLTMB, OC_AIRS,  OC_IASI,  OC_CLPRD, OC_PWIND, OC_PWSOL, &
          OC_WP_JP, OC_ACARS, OC_MODES, OC_GPSRO, OC_GPSGB

!-----------------
! variable numbers
!-----------------
public :: VN_U, VN_V, VN_Z, VN_DZ, VN_RH, VN_RH2M, VN_T, VN_TD, VN_T2M,   &
          VN_TD2M, VN_TS, VN_PTEND, VN_W1, VN_WW, VN_VV, VN_CH, VN_CM,    &
          VN_CL, VN_NH, VN_N_L, VN_C, VN_NS, VN_SDEPTH, VN_E, VN_TRTR,    &
          VN_RR, VN_JJ, VN_GCLG, VN_N, VN_SFALL, VN_PS, VN_DD, VN_FF,     &
          VN_RAWBT, VN_U10M, VN_V10M, VN_Q, VN_VT, VN_HEIGHT, VN_BENDANG, &
          VN_IMPPAR, VN_REFR, VN_ZPD, VN_ZWD, VN_SPD, VN_GUST, VN_P,      &
          VN_TMIN, VN_PRED, VN_N_M, VN_N_H, VN_W, VN_TURB, VN_NFXME,      &
          VN_ICLG, VN_PWC, VN_NUM, VN_RADVEL, VN_FLEV, VN_ELEV, VN_REFL,  &
          VN_RADIANCE, VN_RREFL, VN_CTH, VN_TRH, VN_VGUST,                &
          VN_PWIND, VN_PWSOL, VN_RAD_GL, VN_RAD_DF, VN_RAD_LW, n_vn

!--------------------------------
! satellite instruments (sensors)
!--------------------------------
public :: SI_HIRS, SI_AMSU_A, SI_AMSU_B, SI_AIRS, SI_IASI, SI_SEVIRI

!----------------
! radiosonde type
!----------------
public :: RS_GRAW, RS_BASORA, RS_RU_A_MRZ, RS_RU_MET1, RS_80, RS_VIZ_M2,  &
          RS_DC_MODEM, RS_RU_A_BAR, RS_90_DIG12, RS_RU_ARMA, RS_92_DIG12, &
          RS_92_DIG3, RS_92_AUTO, RS_RU_V_MRZ, RS_RU_V_BAR, RS_SA_DAT4G,  &
          RS_MISS

!-------------------
! tracking technique
!-------------------
public :: TT_NOWIND, TT_AUX_OPTIC, TT_AUX_RANGE, TT_LORANC, TT_SATNAV, &
          TT_NOTSPEC, TT_NORMAL, TT_MISS, TT_AIR_PHASE

!----------------------------
! type of measuring equipment
!----------------------------
public :: TME_PRESS, TME_OPTTHEO, TME_RADTHEO, TME_RADAR, TME_VLFOMEGA, &
          TME_LORANC, TME_WINDPROF, TME_SATNAV, TME_RASS, TME_SODAR, TME_MISS

!--------------------------------------
! solar & infrared radiation correction
!--------------------------------------
public :: RC_NO, RC_CS_CI, RC_CS_IN, RC_CS, RC_SO_IN_AUTO, RC_SO_AUTO, &
          RC_SO_IN_CNTRY, RC_SO_CNTRY, RC_MISS

!---------------------------------------
! ct_nwc cloud type according to NWC SAF
!---------------------------------------
public :: CT_NOPR,        CT_LAND_FREE,   CT_SEA_FREE,     CT_LAND_SNOW,  &
          CT_SEA_ICE,     CT_CUM_VLOW,    CT_STRAT_VLOW,   CT_LOW_CUM,    &
          CT_LOW_STRAT,   CT_MED_CUM,     CT_MED_STRAT,    CT_HI_OP_CUM,  &
          CT_HI_OP_STRAT, CT_VHI_OP_CUM,  CT_VHI_OP_STRAT, CT_HI_ST_THIN, &
          CT_HI_ST_MEAN,  CT_HI_ST_THICK, CT_HI_ST_ABOVE,  CT_FRAC,       &
          CT_UNDEF

!-----------------------------------------------------
! surftype, flg_1dvar, flg_cld (consistent with 1dvar)
!-----------------------------------------------------
public :: ST_SEA, ST_ICE, ST_LAND, ST_HIGHLAND, ST_MISMATCH
public :: D1_DATA, D1_MIN, D1_SUR, D1_CLD
public :: CL_CLEAR, CL_IR_CLOUDY, CL_MW_CLEAR, CL_MW_CLOUDY

!-----------------------------------------------------
! surf_char (model surface characteristics)
!-----------------------------------------------------
public :: MS_LAND, MS_SEA, MS_ICE, MS_NO_ICE, MS_SNOW, MS_NO_SNOW

!-------------------
! level significance
!-------------------
public :: LS_SURFACE, LS_STANDARD, LS_TROPO, LS_MAX, LS_SIGN

!------------------------------------------------
! phase of aircraft flight and roll angle quality
!------------------------------------------------
public :: PH_UNS, PH_LVR, PH_LVW, PH_ASC, PH_DES, PH_MIS
public :: RA_GOOD, RA_BAD, RA_MIS

!------------------------------------------------------------------
! retrieval type: for AMV satellite derived wind computation method
!------------------------------------------------------------------
public :: RT_IR,  RT_IR1,  RT_IR2,  RT_IR3,  &
          RT_VIS, RT_VIS1, RT_VIS2, RT_VIS3, &
          RT_WV,  RT_WV1,  RT_WV2,  RT_WV3

!-----------------------------------------------------
! specification of verification data (runtype, ensmem)
!-----------------------------------------------------
public :: VT_FORECAST, VT_FIRSTGUESS, VT_PREL_ANA, VT_ANALYSIS, VT_INIT_ANA, &
          VT_LIN_ANA
public :: VE_ENS_MEAN, VE_DETERM, VE_ENS_SPREAD, VE_BG_ERROR, VE_TALAGRAND,  &
          VE_VQC_WEIGHT, VE_MEMBER, VE_ENS_MEAN_OBS, VE_BIASCOR

!--------------------------
! observation operator flag
!--------------------------
public :: OF_MISSING, OF_RAD_CLEAR_SKY, OF_RAD_CLOUDY, OF_BT_CLEAR_SKY

!==============================================================================
! Module variables
!==============================================================================
save
logical :: init = .false.   ! flag if module was already initialised

!==============================================================================
! Run Class
!==============================================================================
!----------
! constants
!----------
integer ,parameter :: RC_HAUPT = 0 ,&
                      RC_VOR   = 1 ,&
                      RC_ASS   = 2 ,&
                      RC_TEST  = 3
!--------------
! table entries
!--------------
type(e) ,target :: runclass_entries (4) = (/ &
e(RC_HAUPT ,'HAUPT','','main forecast cycle'),&
e(RC_VOR   ,'VOR  ','','pre-run            '),&
e(RC_ASS   ,'ASS  ','','assimilation cycle '),&
e(RC_TEST  ,'TEST ','','test (offline)     ')/)

!------
! table
!------
type(t_table) ,pointer :: runclass

!==============================================================================
! Observation or Report Status Flags
!==============================================================================
!----------
! constants
!----------
integer ,parameter :: ST_ACCEPTED =  0 ,&
                      ST_ACTIVE   =  1 ,&
                      ST_MERGED   =  3 ,&
                      ST_PASSIVE  =  5 ,&
                      ST_REJECTED =  7 ,&
                      ST_PAS_REJ  =  9 ,&
                      ST_OBS_ONLY = 11 ,&
                      ST_DISMISS  = 13
!--------------
! table entries
!--------------
type (e) ,target :: status_entries (8) = (/                                     &
e(ST_ACCEPTED ,'ACCEPTED','','active and VQC accepted (used in 3D-Var only)  '),&
e(ST_ACTIVE   ,'ACTIVE  ','','used in the assimilation                       '),&
e(ST_MERGED   ,'MERGED  ','','not used, merged into multilevel report        '),&
e(ST_PASSIVE  ,'PASSIVE ','','not used, only monitored                       '),&
e(ST_REJECTED ,'REJECTED','','not used due to suspicious quality             '),&
e(ST_PAS_REJ  ,'PAS_REJ ','','passive and rejected                           '),&
e(ST_OBS_ONLY ,'OBS_ONLY','','observation only, no model equivalent available'),&
e(ST_DISMISS  ,'DISMISS' ,'','dismiss observation, should not appear in file ')/)

!------
! table
!------
type(t_table) ,pointer :: status

!==============================================================================
! Observation or Report Check Flags
!==============================================================================
!----------
! constants
!----------
integer ,parameter :: &
      FL_OBSTYPE     =  0 ,&! passive report type (at obs. location)
      FL_BLACKLIST   =  1 ,&! blacklist (or not whitelist)
      FL_SUSP_LOCT   =  2 ,&! suspicious location or date/time
      FL_TIME        =  3 ,&! time     not in valid range
      FL_AREA        =  4 ,&! location not in valid area
      FL_HEIGHT      =  5 ,&! location not in valid height range
      FL_SURF        =  6 ,&! incorrect surface (land,ice,etc)
      FL_CLOUD       =  7 ,&! cloud check
      FL_PRACTICE    =  8 ,&! bad reporting practice / insufficient data
      FL_DATASET     =  9 ,&! dataset qality flags
      FL_REDUNDANT   = 10 ,&! redundant report
      FL_FLIGHTTRACK = 11 ,&! flight track error flag
      FL_MERGE       = 12 ,&! merged reports (e.g. TEMP ABCD)
      FL_THIN        = 13 ,&! thinning
      FL_RULE        = 14 ,&! complex rule
      FL_OBS_ERR     = 15 ,&! observation error too large
      FL_GROSS       = 16 ,&! gross error flag
      FL_NO_BIASCOR  = 17 ,&! no bias correction available
      FL_FG          = 18 ,&! observation - first guess check
      FL_NO_OBS      = 19 ,&! no observations in report
      FL_OPERATOR    = 20 ,&! observation operator not applicable
      FL_FG_LBC      = 21, &! obs - lateral boundary condition check
      FL_NONE        = 32   ! no flag set
!--------------
! table entries
!--------------
type (e) ,target :: flag_entries (23) = (/                                    &
e(FL_SUSP_LOCT  ,'SUSP_LOCT  ','','suspicious location or date/time     ~ 1'),&
e(FL_TIME       ,'TIME       ','','time     not in valid range          ~ 2'),&
e(FL_AREA       ,'AREA       ','','location not in valid area           ~ 3'),&
e(FL_PRACTICE   ,'PRACTICE   ','','bad reporting practice/insuff. data  ~ 4'),&
e(FL_DATASET    ,'DATASET    ','','dataset qality flags                 ~ 5'),&
e(FL_BLACKLIST  ,'BLACKLIST  ','','blacklist (or not whitelist)         ~ 6'),&
e(FL_HEIGHT     ,'HEIGHT     ','','location not in valid height range   ~ 7'),&
e(FL_SURF       ,'SURF       ','','incorrect surface (land,ice,etc)     ~ 8'),&
e(FL_CLOUD      ,'CLOUD      ','','cloud check                          ~ 9'),&
e(FL_GROSS      ,'GROSS      ','','gross error flag                     ~10'),&
e(FL_OBSTYPE    ,'OBSTYPE    ','','passive report type (at obs.location)~11'),&
e(FL_REDUNDANT  ,'REDUNDANT  ','','redundant report                     ~12'),&
e(FL_FLIGHTTRACK,'FLIGHTTRACK','','flight track error flag              ~13'),&
e(FL_MERGE      ,'MERGE      ','','merged reports (e.g. TEMP ABCD)      ~14'),&
e(FL_THIN       ,'THIN       ','','thinning                             ~15'),&
e(FL_RULE       ,'RULE       ','','complex rule                         ~16'),&
e(FL_NO_BIASCOR ,'NO_BIASCOR ','','no bias correction available         ~17'),&
e(FL_OBS_ERR    ,'OBS_ERR    ','','observation error too large          ~18'),&
e(FL_NO_OBS     ,'NO_OBS     ','','no observations in report            ~19'),&
e(FL_FG         ,'FG         ','','observation - first guess check      ~20'),&
e(FL_FG_LBC     ,'FG_LB      ','','obs- lateral boundary condition check~21'),&
e(FL_OPERATOR   ,'OPERATOR   ','','observation operator not applicable  ~22'),&
e(FL_NONE       ,'NONE       ','','no flag set                             ')/)

!------
! table
!------
type(t_table) ,pointer :: flags

!==============================================================================
! Observation Types
!==============================================================================
!----------
! constants
!----------
integer ,parameter :: OT_SYNOP =  1 ,&
                      OT_AIREP =  2 ,&
                      OT_SATOB =  3 ,&
                      OT_DRIBU =  4 ,&
                      OT_TEMP  =  5 ,&
                      OT_PILOT =  6 ,&
                      OT_SATEM =  7 ,&
                      OT_PAOB  =  8 ,&
                      OT_SCATT =  9 ,&
                      OT_RAD   = 10 ,&
                      OT_GPSRO = 11 ,&
                      OT_GPSGB = 12 ,&
                      OT_RADAR = 13 ,&
                      OT_POWER = 14
integer ,parameter :: n_ot     = 14   ! number of observation types

!--------------
! table entries
!--------------
type (e) ,target :: obstype_entries (n_ot) = (/                           &
e(OT_SYNOP ,'SYNOP','','SYNOP report,               ~(ECMWF convention)'),&
e(OT_AIREP ,'AIREP','','AIREP report,               ~(ECMWF convention)'),&
e(OT_SATOB ,'SATOB','','SATOB report (AMV),         ~(ECMWF convention)'),&
e(OT_DRIBU ,'DRIBU','','DRIBU report,               ~(ECMWF convention)'),&
e(OT_TEMP  ,'TEMP ','','TEMP  report,               ~(ECMWF convention)'),&
e(OT_PILOT ,'PILOT','','PILOT report,               ~(ECMWF convention)'),&
e(OT_SATEM ,'SATEM','','SATEM report,               ~(ECMWF convention)'),&
e(OT_PAOB  ,'PAOB ','','PAOB  report.               ~(ECMWF convention)'),&
e(OT_SCATT ,'SCATT','','Scatterometer report,       ~(ECMWF convention)'),&
e(OT_RAD   ,'RAD'  ,'','Radiances                   ~(ECMWF convention)'),&
e(OT_GPSRO ,'GPSRO','','GPS Radio occultations,       ~(DWD convention)'),&
e(OT_GPSGB ,'GPSGB','','GPS ground based observations ~(DWD convention)'),&
e(OT_RADAR ,'RADAR','','RADAR (volume data)           ~(DWD convention)'),&
e(OT_POWER ,'POWER','','POWER (win, solar) data       ~(DWD convention)')/)

!------
! table
!------
type(t_table) ,pointer :: obstype

!==============================================================================
! Observation Code Types
!==============================================================================
!----------
! constants
!----------
integer ,parameter :: OC_SRSCD =  11 ,&! synop surface report
                      OC_ATSCD =  14 ,&! automatic synop surface report
                  !!!             16     french RADOME
                      OC_AHSCD =  21 ,&! ship synop report
                  !!! OC_ABSCD =  22 ,&! ship synop abbreviated report
                  !!! OC_SHRED =  23 ,&! shred reportKp synop report
                      OC_ATSHS =  24 ,&! automatic ship synop report
                      OC_GPS   = 110 ,&! GPS
                      OC_METAR = 140 ,&! METAR

                      OC_AIRCD = 141 ,&! AIREP report
                      OC_CODAR =  41 ,&! CODAR report
                  !!! OC_COLBA = 241 ,&! COLBA report
                      OC_AMDAR = 144 ,&! AMDAR report
                      OC_ACARS = 145 ,&! ACARS report
                      OC_MODES = 146 ,&! MODE-S report

                      OC_CLPRD =  87 ,&! Cloud (height) product
                      OC_STBCD =  88 ,&! SATOB report
                  !!!             89 ,&! High Resolution VIS wind
                      OC_AMV   =  90 ,&! AMV
                  !!! OC_SST   = 188 ,&! SST report

                      OC_DRBCD = 165 ,&! DRIBU report
                  !!! OC_BATHY =  63 ,&! BATHY report
                  !!!            160     ERS as DRIBU
                      OC_TESAC =  64 ,&! scatterometer

                      OC_LDTCD =  35 ,&! land TEMP report
                      OC_SHTCD =  36 ,&! ship TEMP report
                      OC_TMPMB =  37 ,&! TEMP  mobile
                      OC_PLTMB =  38 ,&! PILOT mobile
                      OC_TDROP = 135 ,&! temp-drop report
                  !!! OC_ROCOB =  39 ,&! ROCOB rep
                  !!! OC_ROCSH =  40 ,&! ROCOB ship report

                      OC_LDPCD =  32 ,&! land pilot report
                      OC_SHPCD =  33 ,&! ship pilot report

                      OC_WP_EU = 132 ,&! European  wind profiler
                      OC_RA_EU = 133 ,&! European sodar/rass report
                      OC_WP_JP = 134 ,&! Japanese  wind profiler
                      OC_PR_US = 136 ,&! wind/profiler/rass report (USA)
                      OC_RAVAD = 137 ,&! radar VAD wind profile report
                      OC_PWIND = 150 ,&! wind  power data
                      OC_PWSOL = 151 ,&! solar power data
                  !!!             34     American wind profiler

                  !!!      8 122 210     Scatterometer

                      OC_ATOVS = 210 ,&! ATOVS satellite data (radiances)
                  !!! OC_STMCD =  86 ,&! satem report
                  !!! OC_STOVS = 186 ,&! high resolution ATOVS satellite data
                  !!! OC_SMSG1 =  71 ,&! MSG_1 satellite retrieval
                  !!! OC_NOA15 = 206 ,&! NOAA15 satellite retrieval
                  !!! OC_NOA16 = 207 ,&! NOAA16 satellite retrieval
                  !!! OC_NOA17 = 208 ,&! NOAA17 satellite retrieval
                  !!! OC_NOA18 = 209 ,&! NOAA18 satellite retrieval
                  !!! OC_GPGFZ =  96 ,&! GPS report processed by GFZ
                  !!! OC_1DVAR = 999 ,&! satellite retrieval (1dvar)
                      OC_SEVIR = 218 ,&! SEVIRI
                      OC_ASCAT = 123 ,&! ASCAT scatterometer
                      OC_QSCAT = 122 ,&! QSCAT scatterometer
                      OC_AIRS  = 216 ,&! AIRS
                      OC_IASI  = 217 ,&! IASI
                      OC_GPSRO = 250 ,&! GPS Radio Occultation
                      OC_GPSGB = 251   ! GPS ground based observations
!--------------
! table entries
!--------------
type (e) ,target :: codetype_entries (38) = (/                        &
e(OC_SRSCD ,'SRSCD','','synop surface report                       '),&
e(OC_ATSCD ,'ATSCD','','automatic synop surface report             '),&
e(OC_AHSCD ,'AHSCD','','ship synop report                          '),&
e(OC_ATSHS ,'ATSHS','','automatic ship synop report                '),&
e(OC_METAR ,'METAR','','METAR                                      '),&
e(OC_GPS   ,'GPS'  ,'','GPS                                        '),&
e(OC_AIRCD ,'AIRCD','','airep report                               '),&
e(OC_CODAR ,'CODAR','','codar report                               '),&
e(OC_AMDAR ,'AMDAR','','amdar report                               '),&
e(OC_MODES ,'MODES','','mode-s report                              '),&
e(OC_ACARS ,'ACARS','','acars report                               '),&
e(OC_CLPRD ,'CLPRD','','cloud (height) product                     '),&
e(OC_STBCD ,'STBCD','','satob report                               '),&
e(OC_AMV   ,'AMV'  ,'','AMV                                        '),&
e(OC_DRBCD ,'DRBCD','','dribu report                               '),&
e(OC_TESAC ,'TESAC','','scatterometer                              '),&
e(OC_LDTCD ,'LDTCD','','land temp report                           '),&
e(OC_SHTCD ,'SHTCD','','temp ship report                           '),&
e(OC_TDROP ,'TDROP','','temp-drop report                           '),&
e(OC_TMPMB ,'TMPMB','','temp mobile                                '),&
e(OC_LDPCD ,'LDPCD','','land pilot report                          '),&
e(OC_SHPCD ,'SHPCD','','ship pilot report                          '),&
e(OC_PLTMB ,'PLTMB','','pilot mobile                               '),&
e(OC_ATOVS ,'ATOVS','','ATOVS satellite data (1dvar)               '),&
e(OC_WP_EU ,'WP_EU','','European  wind profiler                    '),&
e(OC_RA_EU ,'RA_EU','','European sodar/rass report                 '),&
e(OC_WP_JP ,'WP_JP','','Japanese  wind profiler                    '),&
e(OC_PR_US ,'PR_US','','wind/profiler/rass report (USA)            '),&
e(OC_RAVAD ,'RAVAD','','radar VAD wind profile report              '),&
e(OC_PWIND ,'PWIND','','wind power data                            '),&
e(OC_PWSOL ,'PWSOL','','solar power data                           '),&
e(OC_SEVIR ,'SEVIR','','SEVIRI                                     '),&
e(OC_ASCAT ,'ASCAT','','ASCAT scatterometer                        '),&
e(OC_QSCAT ,'QSCAT','','QSCAT scatterometer                        '),&
e(OC_AIRS  ,'AIRS' ,'','AIRS                                       '),&
e(OC_IASI  ,'IASI' ,'','IASI                                       '),&
e(OC_GPSRO ,'GPSRO','','GPS Radio Occultation                      '),&
e(OC_GPSGB ,'GPSGB','','GPS ground based observations              ')/)


!------
! table
!------
type(t_table) ,pointer :: codetype

!==============================================================================
! Variable Numbers
!==============================================================================
!----------
! constants
!----------
integer ,parameter :: VN_U       =   3 ,&! u-component of wind
                      VN_V       =   4 ,&! v-component of wind
                      VN_Z       =   1 ,&! geopotential
                      VN_DZ      =  57 ,&! thickness
                      VN_PWC     =   9 ,&! precipitable water content kg/m**2
                      VN_TRH     =  28 ,&! transformed relative humidity
                      VN_RH      =  29 ,&! relative humidity
                      VN_RH2M    =  58 ,&! 2 metre relative humidity
                      VN_T       =   2 ,&! upper air temperature
                      VN_TD      =  59 ,&! upper air dew point
                      VN_T2M     =  39 ,&! 2 metre temperature

                      VN_TD2M    =  40 ,&! 2 metre dew point
                      VN_TS      =  11 ,&! surface temperature
                      VN_PTEND   =  30 ,&! pressure tendency
                      VN_W1      =  60 ,&! past weather
                      VN_WW      =  61 ,&! present weather
                      VN_VV      =  62 ,&! visibility
                      VN_CH      =  63 ,&! type of high clouds
                      VN_CM      =  64 ,&! type of middle clouds
                      VN_CL      =  65 ,&! type of low clouds
                      VN_NH      =  66   ! cloud base height

integer ,parameter :: VN_N_L     =  67 ,&! low cloud amount
                  !!! VN_HSHS    =  68 ,&! additional cloud group height
                      VN_C       =  69 ,&! additional cloud group type
                      VN_NS      =  70 ,&! additional cloud group amount
                      VN_SDEPTH  =  71 ,&! snow depth
                      VN_E       =  72 ,&! state of ground
                  !!! VN_TGTG    =  73 ,&! ground temperature
                  !!! VN_SPSP1   =  74 ,&! special phenomena 1

                  !!! VN_SPSP2   =  75 ,&! special phenomena 2
                  !!! VN_RS      =  76 ,&! ice code type

                  !!! VN_ESES    =  77 ,&! ice thickness              (1751)
                  !!! VN_IS      =  78 ,&! ice                        (1751)
                      VN_TRTR    =  79 ,&! time period of information (h)
                      VN_RR      =  80 ,&! precipitation amount       (kg/m^2)
                      VN_JJ      =  81   ! maximum temperature        (K)
                  !!! VN_VS      =  82 ,&! ship speed                 (m/s)
                  !!! VN_DS      =  83 ,&! ship direction             (degree)
                  !!! VN_HWHW    =  84 ,&! wave height                (m)
                  !!! VN_PWPW    =  85 ,&! wave period                (s)
                  !!! VN_DWDW    =  86   ! wave direction             (degree)

integer ,parameter :: VN_GCLG    =  87 ,&! general cloud group         Table
                  !!! VN_RHLC    =  88 ,&! relative humidity from low clouds
                  !!! VN_RHMC    =  89 ,&! relative humidity from middle clouds
                  !!! VN_RHHC    =  90 ,&! relative humidity from high clouds
                      VN_N       =  91 ,&! total cloud amount         (20011)
                      VN_SFALL   =  92 ,&! 6h snow fall               (m)
                      VN_PS      = 110 ,&! surface pressure           (Pa)
                      VN_DD      = 111 ,&! wind direction             (degree)
                      VN_FF      = 112 ,&! wind speed                 (m/s)

                      VN_REFL    = 118 ,&! reflectivity               (0..1)
                      VN_RAWBT   = 119 ,&! brightness temperature     (K)
                      VN_RADIANCE= 120 ,&! radiance                   (W/sr/m^3)

                  !!! VN_SATCL   = 121 ,&! cloud amount from satellite(%)
                  !!! VN_SCATSS  = 122 ,&! backscatter                (dB)
                  !!! VN_DU      =   5 ,&! wind shear u-component     (1/s)
                  !!! VN_DV      =   6 ,&! wind shear v-component     (1/s)
                      VN_U10M    =  41 ,&! 10m u-component of wind    (m/s)
                      VN_V10M    =  42 ,&! 10m v-component of wind    (m/s)
                  !!! VN_RHLAY   =  19 ,&! layer relative humidity
                  !!! VN_AUXIL   = 200 ,&! auxiliary variable

                  !!! VN_CLLQW   = 123 ,&! cloud liquid water         (kg/kg)
                  !!! VN_SCATDD  = 124 ,&! ambiguous v-component      (m/s)
                  !!! VN_SCATFF  = 125 ,&! ambiguous u-component      (m/s)
                      VN_Q       =   7 ,&! specific humidity          (kg/kg)
                  !!! VN_SCATWD  = 126 ,&! ambiguous wind direction   (degree)
                  !!! VN_SCARWS  = 127 ,&! ambiguous wind speed       (m/s)
                      VN_W       =   8 ,&! vertical speed             (m/s)
                      VN_VT      =  56 ,&! virtual temperature        (K)
                  !!! VN_O3LAY   = 130 ,&! ozone layer density        (kg/m^2)
                      VN_CTH     = 155 ,&! cloud top height           (m)
                      VN_HEIGHT  = 156 ,&! height                     (m)
                      VN_FLEV    = 157   ! nominal flight level       (m)
                  !!! VN_XXXX    = 215   ! SSM/I multilevel variable
                  !!! --------- new --------------------
                  !!!             206     ozone                       (Dopson)
                  !!!             160     past weather                 numeric
                  !!!             130     pressure tendency charact.   numeric
                  !!!              12     sea water temperature        numeric
integer ,parameter :: VN_TURB    = 244 ,&! degree of turbulence         WMO table 011031
                      VN_NFXME   = 249 ,&! max wind speed (10min mean)(m/s)
                      VN_RREFL   = 192 ,&! radar reflectivity         (Db)
                      VN_RADVEL  = 193 ,&! radial velocity            (m/s)
                      VN_PDELAY  = 128 ,&! atmospheric path delay     (m)
                      VN_BENDANG = 162 ,&! radio occ. bending angle   (Rad)
                      VN_ICLG    =  95 ,&! individual cloud layer group Table
                      VN_N_M     =  93 ,&! middle cloud amount         WMO table 020011
                      VN_N_H     =  94 ,&! high   cloud amount         WMO table 020011
                  !!! ---------------------------------------------
                      VN_ELEV    = 158 ,&! elevation                  (degree)
                      VN_PWIND   = 230 ,&! wind  power data
                      VN_PWSOL   = 231 ,&! solar power data
                      VN_RAD_GL  = 237 ,&! 1-h global solar radiation (J/m2) 
                      VN_RAD_DF  = 238 ,&! 1-h diffuse solar radiat.  (J/m2) 
                      VN_RAD_LW  = 239 ,&! 1-h long-wave radiation    (J/m2) 
                      VN_IMPPAR  = 252 ,&! impact parameter           (m)
                      VN_REFR    = 248 ,&! refractivity
                      VN_ZPD     = 245 ,&! zenith path delay
                      VN_ZWD     = 246 ,&! zenith wet delay
                      VN_SPD     = 247 ,&! slant path delay
                      VN_VGUST   = 240 ,&! vertical gust (aircrafts)  (m/s)
                      VN_GUST    = 242 ,&! wind gust                  (m/s)
                      VN_P       = 251 ,&! pressure                   (Pa)
                      VN_TMIN    = 243 ,&! minimum temperature        (K)
                      VN_PRED    = 241 ,&! reduced pressure           (Pa)
                      VN_NUM     =   0   ! ordinal (channel) number   (  )
integer, parameter :: n_vn       =  72   ! number of variable numbers

!--------------
! table entries
!--------------
type (e) ,target :: varno_entries (n_vn) = (/                               &
e(VN_NUM     ,'NUM'     ,''          ,'ordinal (channel) number         ~L '),&
e(VN_U       ,'U'       ,'m/s'       ,'u-component of wind                 '),&
e(VN_V       ,'V'       ,'m/s'       ,'v-component of wind                 '),&
e(VN_W       ,'W'       ,'m/s'       ,'vertical velocity                   '),&
e(VN_Z       ,'Z'       ,'(m/s)**2'  ,'geopotential               ~(also L)'),&
e(VN_DZ      ,'DZ'      ,'(m/s)**2'  ,'thickness                           '),&
e(VN_PWC     ,'PWC'     ,'kg/m**2'   ,'precipitable water content          '),&
e(VN_TRH     ,'TRH'     ,'0..1'      ,'transformed relative humidity       '),&
e(VN_RH      ,'RH'      ,'0..1'      ,'relative humidity                   '),&
e(VN_RH2M    ,'RH2M'    ,'0..1'      ,'2 metre relative humidity           '),&
e(VN_T       ,'T'       ,'K'         ,'upper air temperature               '),&
e(VN_TD      ,'TD'      ,'K'         ,'upper air dew point                 '),&
e(VN_T2M     ,'T2M'     ,'K'         ,'2 metre temperature                 '),&
e(VN_TD2M    ,'TD2M'    ,'K'         ,'2 metre dew point                   '),&
e(VN_TS      ,'TS'      ,'K'         ,'surface temperature                 '),&
e(VN_PTEND   ,'PTEND'   ,'Pa/3h'     ,'pressure tendency                   '),&
e(VN_W1      ,'W1'      ,'WMO 020004','past weather                        '),&
e(VN_WW      ,'WW'      ,'WMO 020003','present weather                     '),&
e(VN_VV      ,'VV'      ,'m'         ,'visibility                          '),&
e(VN_CH      ,'CH'      ,'WMO 020012','type of high clouds                 '),&
e(VN_CM      ,'CM'      ,'WMO 020012','type of middle clouds               '),&
e(VN_CL      ,'CL'      ,'WMO 020012','type of low clouds                  '),&
e(VN_NH      ,'NH'      ,'m'         ,'cloud base height                   '),&
e(VN_N_L     ,'N_L'     ,'WMO 020011','low cloud amount                    '),&
e(VN_N_M     ,'N_M'     ,'WMO 020011','medium cloud amount                 '),&
e(VN_N_H     ,'N_H'     ,'WMO 020011','high cloud amount                   '),&
e(VN_C       ,'C'       ,'WMO  500'  ,'additional cloud group type      ~L '),&
e(VN_NS      ,'NS'      ,'WMO 2700'  ,'additional cloud group amount       '),&
e(VN_SDEPTH  ,'SDEPTH'  ,'m'         ,'snow depth                          '),&
e(VN_E       ,'E'       ,'WMO 020062','state of ground                     '),&
e(VN_TRTR    ,'TRTR'    ,'h'         ,'time period of information       ~L '),&
e(VN_RR      ,'RR'      ,'kg/m**2'   ,'precipitation amount                '),&
e(VN_JJ      ,'JJ'      ,'K'         ,'maximum temperature                 '),&
e(VN_GCLG    ,'GCLG'    ,'Table 6'   ,'general cloud group                 '),&
e(VN_N       ,'N'       ,'WMO 020011','total cloud amount                  '),&
e(VN_SFALL   ,'SFALL'   ,'m'         ,'6h snow fall                        '),&
e(VN_PS      ,'PS'      ,'Pa'        ,'surface (station) pressure          '),&
e(VN_DD      ,'DD'      ,'degree'    ,'wind direction                      '),&
e(VN_FF      ,'FF'      ,'m/s'       ,'wind speed                          '),&
e(VN_REFL    ,'REFL'    ,'0..1'      ,'reflectivity                        '),&
e(VN_RAWBT   ,'RAWBT'   ,'K'         ,'brightness temperature              '),&
e(VN_RADIANCE,'RADIANCE','W/sr/m**3' ,'radiance                            '),&
e(VN_U10M    ,'U10M'    ,'m/s'       ,'10m u-component of wind             '),&
e(VN_V10M    ,'V10M'    ,'m/s'       ,'10m v-component of wind             '),&
e(VN_Q       ,'Q'       ,'kg/kg'     ,'specific humidity                   '),&
e(VN_VT      ,'VT'      ,'K'         ,'virtual temperature                 '),&
e(VN_CTH     ,'CTH'     ,'m'         ,'cloud top height                    '),&
e(VN_HEIGHT  ,'HEIGHT'  ,'m'         ,'height                           ~L '),&
e(VN_FLEV    ,'FLEV'    ,'m'         ,'nominal flight level             ~L '),&
e(VN_ELEV    ,'ELEV'    ,'degree'    ,'elevation                        ~L '),&
e(VN_PWIND   ,'PWIND'   ,'W'         ,'wind power data                     '),&
e(VN_PWSOL   ,'PWSOL'   ,'W'         ,'solar power data                    '),&
e(VN_RREFL   ,'RREFL'   ,'Db'        ,'radar reflectivity                  '),&
e(VN_RADVEL  ,'RADVEL'  ,'m/s'       ,'radial velocity                     '),&
e(VN_PDELAY  ,'PDELAY'  ,'m'         ,'atmospheric path delay              '),&
e(VN_BENDANG ,'BENDANG' ,'rad'       ,'bending angle                       '),&
e(VN_IMPPAR  ,'IMPPAR'  ,'m'         ,'impact parameter                 ~L '),&
e(VN_REFR    ,'REFR'    ,''          ,'refractivity                        '),&
e(VN_ZPD     ,'ZPD'     ,''          ,'zenith path delay                   '),&
e(VN_ZWD     ,'ZWD'     ,''          ,'zenith wet delay                    '),&
e(VN_SPD     ,'SPD'     ,''          ,'slant path delay                    '),&
e(VN_VGUST   ,'VGUST'   ,'m/s'       ,'vertical gust (aircrafts)           '),&
e(VN_GUST    ,'GUST'    ,'m/s'       ,'wind gust                           '),&
e(VN_P       ,'P'       ,'Pa'        ,'pressure                         ~L '),&
e(VN_TMIN    ,'TMIN'    ,'K'         ,'minimum temperature                 '),&
e(VN_RAD_GL  ,'RAD_GL'  ,'J/m**2'    ,'global solar radiation              '),&
e(VN_RAD_DF  ,'RAD_DF'  ,'J/m**2'    ,'diffuse solar radiation             '),&
e(VN_RAD_LW  ,'RAD_LW'  ,'J/m**2'    ,'long-wave (downward) radiation      '),&
e(VN_PRED    ,'PRED'    ,'Pa'        ,'reduced pressure                    '),&
e(VN_TURB    ,'TURB'    ,'WMO 011031','degree of turbulence                '),&
e(VN_NFXME   ,'NFXME'   ,'m/s'       ,'max wind speed (10min mean)         '),&
e(VN_ICLG    ,'ICLG'    ,'Table 7'   ,'individual cloud layer group        ')/)
!------
! table
!------
type(t_table) ,pointer :: varno

!==============================================================================
! General cloud group
!==============================================================================
!----------
! constants
!----------
integer ,parameter :: &!
!       GC_SGBP =  0 ,&! bit position for vert.signif.      (WMO 008002)
        GC_CLBP =  0 ,&! bit position for cloud amount      (WMO 020011)
        GC_LCBP =  4 ,&! bit position for low cloud type    (WMO 020012)
        GC_MCBP = 10 ,&! bit position for middle cloud type (WMO 020012)
        GC_HCBP = 16 ,&! bit position for high cloud type   (WMO 020012)
!       GC_SGOC =  6 ,&! no.bits used for vert.signif.      (WMO 008002)
        GC_CLOC =  4 ,&! no.bits used for cloud amount      (WMO 020011)
        GC_LCOC =  6 ,&! no.bits used for low cloud type    (WMO 020012)
        GC_MCOC =  6 ,&! no.bits used for middle cloud type (WMO 020012)
        GC_HCOC =  6   ! no.bits used for high cloud type   (WMO 020012)

!--------------
! table entries
!--------------
type (e) ,target :: gen_cg_entries (8) = (/&
!e(GC_SGBP,'SGBP'   ,'WMO 008002' ,'bit position for vert.signif.'     ),&
 e(GC_CLBP,'CLBP'   ,'WMO 020011' ,'bit position for cloud amount'     ),&
 e(GC_LCBP,'LCBP'   ,'WMO 020012' ,'bit position for low cloud type'   ),&
 e(GC_MCBP,'MCBP'   ,'WMO 020012' ,'bit position for middle cloud type'),&
 e(GC_HCBP,'HCBP'   ,'WMO 020012' ,'bit position for high type'        ),&
!e(GC_SGOC,'SGOC'   ,'WMO 008002' ,'no.bits used for vert.signif.'     ),&
 e(GC_CLOC,'CLOC'   ,'WMO 020011' ,'no.bits used for cloud amount'     ),&
 e(GC_LCOC,'LCOC'   ,'WMO 020012' ,'no.bits used for low cloud type'   ),&
 e(GC_MCOC,'MCOC'   ,'WMO 020012' ,'no.bits used for middle cloud type'),&
 e(GC_HCOC,'HCOC'   ,'WMO 020012' ,'no.bits used for high type'        )/)

!------
! table
!------
type(t_table) ,pointer :: gen_cg

!==============================================================================
! Individual cloud group
!==============================================================================
!----------
! constants
!----------
integer ,parameter :: &!
!       IC_SGBP =  0 ,&! bit position for vert.signif. (WMO 008002)
        IC_CLBP =  0 ,&! bit position for cloud amount (WMO 020011)
        IC_CTBP =  4 ,&! bit position for cloud type   (WMO 020012)
        IC_BSBP = 10 ,&! bit position for cloud base height     (m)
!       IC_SGOC =  6 ,&! no.bits used for vert.signif. (WMO 008002)
        IC_CLOC =  4 ,&! no.bits used for cloud amount (WMO 020011)
        IC_CTOC =  6 ,&! no.bits used for cloud type   (WMO 020012)
        IC_BSOC = 14   ! no.bits used for cloud base height     (m)
!--------------
! table entries
!--------------
type (e) ,target :: ind_cg_entries (6) = (/&
!e(IC_SGBP,'SGBP'   ,'WMO 008002' ,'bit position for vert.signif.'     ),&
 e(IC_CLBP,'CLBP'   ,'WMO 020011' ,'bit position for cloud amount'     ),&
 e(IC_CTBP,'CTBP'   ,'WMO 020012' ,'bit position for cloud type'       ),&
 e(IC_BSBP,'BSBP'   ,'m'          ,'bit position for cloud base height'),&
!e(IC_SGOC,'SGOC'   ,'WMO 008002' ,'no.bits used for vert.signif.'     ),&
 e(IC_CLOC,'CLOC'   ,'WMO 020011' ,'no.bits used for cloud amount'     ),&
 e(IC_CTOC,'CTOC'   ,'WMO 020012' ,'no.bits used for cloud type'       ),&
 e(IC_BSOC,'BSOC'   ,'m'          ,'no.bits used for cloud base height')/)

!------
! table
!------
type(t_table) ,pointer :: ind_cg

!==============================================================================
! Satellite Instruments (Sensors)
!==============================================================================
!----------
! constants
!----------
integer ,parameter :: SI_HIRS      =  0 ,&! RTTOV8 channels numbers : 1-  19
                  !!! SI_MSU       =  1 ,&!                           1-   4
                  !!! SI_SSU       =  2 ,&!                           1-   3
                      SI_AMSU_A    =  3 ,&!                           1-  15
                      SI_AMSU_B    =  4 ,&!                           1-   5
                  !!! SI_AVHRR     =  5 ,&!                           1-   3
                  !!! SI_SSMI      =  6 ,&!                           1-   4
                  !!! SI_VTPR1     =  7 ,&!                           1-   8
                  !!! SI_VTPR2     =  8 ,&!                           1-   8
                  !!! SI_TMI       =  9 ,&!                           1-   9
                  !!! SI_SSMIS     = 10 ,&!                           1-  21
                      SI_AIRS      = 11 ,&!                           1-2378
                  !!! SI_HSB       = 12 ,&!                           1-   4
                  !!! SI_MODIS     = 13 ,&!                           1-  17
                  !!! SI_ATSR      = 14 ,&!                           1-   3
                  !!! SI_MHS       = 15 ,&!                           1-   5
                      SI_IASI      = 16 ,&!                           1-8461
                  !!! SI_AMSR      = 17 ,&!                           1-   7
                  !!! SI_MVIRI     = 20 ,&!                           1-   2
                      SI_SEVIRI    = 21   !                           1-   8
                  !!! SI_GOES_IM   = 22 ,&! GOES Imager               1-   4
                  !!! SI_GOES_SND  = 23 ,&! GOES Sounder              1-  18
                  !!! SI_GMS_IM    = 24 ,&! GMS/MTSAT imager          1-   4
                  !!! SI_FY2_VISSR = 25 ,&!                           1-   2
                  !!! SI_FY1_VISSR = 26 ,&!                           1-   3
                  !!! SI_CRIS      = 27 ,&!
                  !!! SI_CMSS      = 28 ,&!
                  !!! SI_VIIRS     = 29 ,&!
                  !!! SI_WINDSAT   = 30   !                           1-   5
!--------------
! table entries
!--------------
type (e) ,target :: satsens_entries (6) = (/&
e(SI_HIRS   ,'HIRS'   ,'' ,''),&
e(SI_AMSU_A ,'AMSU_A' ,'' ,''),&
e(SI_AMSU_B ,'AMSU_B' ,'' ,''),&
e(SI_AIRS   ,'AIRS'   ,'' ,''),&
e(SI_IASI   ,'IASI'   ,'' ,''),&
e(SI_SEVIRI ,'SEVIRI' ,'' ,'')/)

!------
! table
!------
type(t_table) ,pointer :: satsens


!==============================================================================
! radiosonde type
!==============================================================================
integer ,parameter :: RS_GRAW     =  17 ,&! Graw   (D)
                      RS_BASORA   =  26 ,&! Basora (CH)
                      RS_RU_A_MRZ =  27 ,&! AVK-MRZ (Russia)
                      RS_RU_MET1  =  28 ,&! Meteorit Marz2-1 (Russia)
                      RS_80       =  37 ,&! Vaisala RS 80
                      RS_VIZ_M2   =  49 ,&! VIZ MARK II (USA)
                      RS_DC_MODEM =  57 ,&! M2K2-DC Modem (France)
                      RS_RU_A_BAR =  58 ,&! AVK-BAR (Russia)
                      RS_90_DIG12 =  71 ,&! Vaisala RS 90 Digicora I,II or Marvin
                      RS_RU_ARMA  =  75 ,&! AVK-MRZ-ARMA (Russia)
                      RS_92_DIG12 =  79 ,&! Vaisala RS 92 Digicora I,II or Marvin
                      RS_92_DIG3  =  80 ,&! Vaisala RS 92 Digicora III
                      RS_92_AUTO  =  81 ,&! Vaisala RS 92 Autosonde
                      RS_RU_V_MRZ =  88 ,&! MARL-A or Vektor-M-MRZ (Russia)
                      RS_RU_V_BAR =  89 ,&! MARL-A or Vektor-M-BAR (Russia)
                      RS_SA_DAT4G =  99 ,&! BAT-4G (South Africa)
                      RS_MISS     = 255   ! missing value

!--------------
! table entries
!--------------
type (e) ,target :: rsondtype_entries (17) = (/                          &
e(RS_GRAW    ,'GRAW'       ,'' ,'Graw   (D)')                           ,&
e(RS_BASORA  ,'BASORA'     ,'' ,'Basora (CH)')                          ,&
e(RS_RU_A_MRZ,'RS_RU_A_MRZ','' ,'AVK-MRZ (Russia)')                     ,&
e(RS_RU_MET1 ,'RS_RU_MET1' ,'' ,'Meteorit Marz2-1 (Russia)')            ,&
e(RS_80      ,'RS_80'      ,'' ,'Vaisala RS 80')                        ,&
e(RS_VIZ_M2  ,'RS_VIZ_M2'  ,'' ,'VIZ MARK II (USA)')                    ,&
e(RS_DC_MODEM,'RS_DC_MODEM','' ,'M2K2-DC Modem (France)')               ,&
e(RS_RU_A_BAR,'RS_RU_A_BAR','' ,'AVK-BAR (Russia)')                     ,&
e(RS_90_DIG12,'RS_90_DIG12','' ,'Vaisala RS 90 Digicora I,II or Marvin'),&
e(RS_RU_ARMA ,'RS_RU_ARMA' ,'' ,'AVK-MRZ-ARMA (Russia)')                ,&
e(RS_92_DIG12,'RS_92_DIG12','' ,'Vaisala RS 92 Digicora I,II or Marvin'),&
e(RS_92_DIG3 ,'RS_92_DIG3' ,'' ,'Vaisala RS 92 Digicora III')           ,&
e(RS_92_AUTO ,'RS_92_AUTO' ,'' ,'Vaisala RS 92 Autosonde')              ,&
e(RS_RU_V_MRZ,'RS_RU_V_MRZ','' ,'MARL-A or Vektor-M-MRZ (Russia)')      ,&
e(RS_RU_V_BAR,'RS_RU_V_BAR','' ,'MARL-A or Vektor-M-BAR (Russia)')      ,&
e(RS_SA_DAT4G,'RS_SA_DAT4G','' ,'BAT-4G (South Africa)')                ,&
e(RS_MISS    ,'MISS'       ,'' ,'missing value')/)

!------
! table
!------
type(t_table) ,pointer :: rsondtype

!==============================================================================
! tracking technique
!==============================================================================
integer ,parameter :: TT_NOWIND    =   0 ,&! no windfinding
                      TT_AUX_OPTIC =   2 ,&! automatic with aux. optical direction finding
                      TT_AUX_RANGE =   3 ,&! automatic with auxiliary ranging
                      TT_LORANC    =   6 ,&! automatic cross chain Loran-C
                      TT_SATNAV    =   8 ,&! automatic satellite navigation
                      TT_NOTSPEC   =  19 ,&! tracking technique not specified
                      TT_NORMAL    =  70 ,&! all systems in normal operation
                      TT_MISS      = 127 ,&! missing value
                      TT_AIR_PHASE =  18   ! aircraft obs: flight phase
                                           !   determined in data assimilation
!--------------
! table entries
!--------------
type (e) ,target :: trackteqn_entries (9) = (/                                  &
e(TT_NOWIND   ,'NOWIND'   ,'' ,'no windfinding')                               ,&
e(TT_AUX_OPTIC,'AUX_OPTIC','' ,'automatic with aux. optical direction finding'),&
e(TT_AUX_RANGE,'AUX_RANGE','' ,'automatic with auxiliary ranging')             ,&
e(TT_LORANC   ,'LORANC'   ,'' ,'automatic cross chain Loran-C')                ,&
e(TT_SATNAV   ,'SATNAV'   ,'' ,'automatic satellite navigation')               ,&
e(TT_NOTSPEC  ,'NOTSPEC'  ,'' ,'tracking technique not specified')             ,&
e(TT_NORMAL   ,'NORMAL'   ,'' ,'all systems in normal operation')              ,&
e(TT_MISS     ,'MISS'     ,'' ,'missing value')                                ,&
e(TT_AIR_PHASE,'AIR_PHASE','' ,'aircraft obs: flight phase from data assimil.')/)

!------
! table
!------
type(t_table) ,pointer :: trackteqn

!==============================================================================
! type of measuring equipment used
!==============================================================================
integer ,parameter :: &
  TME_PRESS    = 0 ,&! pressure instrument associated with wind measuring equipment
  TME_OPTTHEO  = 1 ,&! optical theodolit
  TME_RADTHEO  = 2 ,&! radio theodolite
  TME_RADAR    = 3 ,&! radar
  TME_VLFOMEGA = 4 ,&! VLF-Omega
  TME_LORANC   = 5 ,&! Loran-C
  TME_WINDPROF = 6 ,&! wind profiler
  TME_SATNAV   = 7 ,&! satellite navigation
  TME_RASS     = 8 ,&! radio acoustic sounding system (RASS)
  TME_SODAR    = 9 ,&! SODAR
  TME_MISS     = 15  ! missing

!--------------
! table entries
!--------------
type (e) ,target :: meas_equip_entries (11) = (/     &
e(TME_PRESS   ,'PRESS'   ,'','pressure instrument associated with wind measuring equipment'),&
e(TME_OPTTHEO ,'OPTTHEO' ,'','optical theodolit'   ),&
e(TME_RADTHEO ,'RADTHEO' ,'','radio theodolite'    ),&
e(TME_RADAR   ,'RADAR'   ,'','radar'               ),&
e(TME_VLFOMEGA,'VLFOMEGA','','VLF-Omega'           ),&
e(TME_WINDPROF,'WINDPROF','','wind profiler'       ),&
e(TME_LORANC  ,'LORANC'  ,'','Loran-C'             ),&
e(TME_SATNAV  ,'SATNAV'  ,'','satellite navigation'),&
e(TME_RASS    ,'RASS'    ,'','radio acoustic sounding system (RASS)'),&
e(TME_SODAR   ,'SODAR'   ,'','SODAR'               ),&
e(TME_MISS    ,'MISS'    ,'','missing'             )/)

!------
! table
!------
type(t_table) ,pointer :: meas_equip

!==============================================================================
! solar & infrared radiation correction (NSR, 002013)
!==============================================================================
integer ,parameter ::  &
  RC_NO          =  0 ,&! no correction
  RC_CS_CI       =  1 ,&! CIMO solar + CIMO infrared corrected
  RC_CS_IN       =  2 ,&! CIMO solar + infrared corrected
  RC_CS          =  3 ,&! CIMO solar corrected only
  RC_SO_IN_AUTO  =  4 ,&! solar + infrared corr., automatic. by rsond. system
  RC_SO_AUTO     =  5 ,&! solar corrected automatically by radiosonde system
  RC_SO_IN_CNTRY =  6 ,&! solar + infrared corr. as specified by country
  RC_SO_CNTRY    =  7 ,&! solar corrected by country
  RC_MISS        = 15   ! missing

!--------------
! table entries
!--------------
type (e) ,target :: radiation_corr_entries (9) = (/&
e(RC_NO         ,'NO'         ,'','no correction'),&
e(RC_CS_CI      ,'CS_CI'      ,'','CIMO solar + CIMO infrared corrected'),&
e(RC_CS_IN      ,'CS_IN'      ,'','CIMO solar + infrared corrected'),&
e(RC_CS         ,'CS'         ,'','CIMO solar corrected only'),&
e(RC_SO_IN_AUTO ,'SO_IN_AUTO' ,'','solar + infrared corr., automatic. by rsond. system'),&
e(RC_SO_AUTO    ,'SO_AUTO'    ,'','solar corrected automatically by radiosonde system'),&
e(RC_SO_IN_CNTRY,'SO_IN_CNTRY','','solar + infrared corr. as specified by country'),&
e(RC_SO_CNTRY   ,'SO_CNTRY'   ,'','solar corrected by country'),&
e(RC_MISS       ,'MISS'       ,'','missing')/)

!------
! table
!------
type(t_table) ,pointer :: radiation_corr

!==============================================================================
! surf_char: model surface characteristics (bit pattern)
!==============================================================================

integer ,parameter :: MS_LAND    = 0  ! ( 1)
integer ,parameter :: MS_SEA     = 1  ! ( 2)
integer ,parameter :: MS_ICE     = 2  ! ( 4)
integer ,parameter :: MS_NO_ICE  = 3  ! ( 8)
integer ,parameter :: MS_SNOW    = 4  ! (16)
integer ,parameter :: MS_NO_SNOW = 5  ! (32)

type (e) ,target :: surf_char_entries (6) = (/               &
e(MS_LAND   ,'LAND'   ,'','set if some fraction is covered by land'       ),&
e(MS_SEA    ,'SEA'    ,'','set if some fraction is covered by sea'        ),&
e(MS_ICE    ,'ICE'    ,'','set if some fraction is covered by sea-ice'    ),&
e(MS_NO_ICE ,'NO_ICE' ,'','set if some fraction is not covered by sea-ice'),&
e(MS_SNOW   ,'SNOW'   ,'','set if some fraction is covered by snow'       ),&
e(MS_NO_SNOW,'NO_SNOW','','set if some fraction is not covered by snow'   )/)

type(t_table) ,pointer :: surf_char

!=======================================
! ct_nwc cloud type according to NWC SAF
!=======================================

integer ,parameter :: CT_NOPR         =  0 ! non-processed containing no data or corrupted data
integer ,parameter :: CT_LAND_FREE    =  1 ! cloud free land
integer ,parameter :: CT_SEA_FREE     =  2 ! cloud free sea
integer ,parameter :: CT_LAND_SNOW    =  3 ! land contaminated by snow
integer ,parameter :: CT_SEA_ICE      =  4 ! sea contaminated by snow/ice
integer ,parameter :: CT_CUM_VLOW     =  5 ! very low and cumuliform clouds
integer ,parameter :: CT_STRAT_VLOW   =  6 ! very low and stratiform clouds
integer ,parameter :: CT_LOW_CUM      =  7 ! low and cumuliform clouds
integer ,parameter :: CT_LOW_STRAT    =  8 ! low and stratiform clouds
integer ,parameter :: CT_MED_CUM      =  9 ! medium and cumuliform clouds
integer ,parameter :: CT_MED_STRAT    = 10 ! medium and stratiform clouds
integer ,parameter :: CT_HI_OP_CUM    = 11 ! high opaque and cumuliform clouds
integer ,parameter :: CT_HI_OP_STRAT  = 12 ! high opaque and stratiform clouds
integer ,parameter :: CT_VHI_OP_CUM   = 13 ! very high opaque and cumuliform clouds
integer ,parameter :: CT_VHI_OP_STRAT = 14 ! very high opaque and stratiform clouds
integer ,parameter :: CT_HI_ST_THIN   = 15 ! high semitransparent thin clouds
integer ,parameter :: CT_HI_ST_MEAN   = 16 ! high semitransparent meanly thick clouds
integer ,parameter :: CT_HI_ST_THICK  = 17 ! high semitransparent thick clouds
integer ,parameter :: CT_HI_ST_ABOVE  = 18 ! high semitransparent above low or medium clouds
integer ,parameter :: CT_FRAC         = 19 ! fractional clouds (sub-pixel water clouds)
integer ,parameter :: CT_UNDEF        = 20 ! undefined (undefined by CMa)

type (e) ,target :: ct_nwc_entries (21) = (/                                              &
e(CT_NOPR        ,'NOPR'        ,'','non-processed containing no data or corrupted data'),&
e(CT_LAND_FREE   ,'LAND_FREE'   ,'','cloud free land'                                   ),&
e(CT_SEA_FREE    ,'SEA_FREE'    ,'','cloud free sea'                                    ),&
e(CT_LAND_SNOW   ,'LAND_SNOW'   ,'','land contaminated by snow'                         ),&
e(CT_SEA_ICE     ,'SEA_ICE'     ,'','sea contaminated by snow/ice'                      ),&
e(CT_CUM_VLOW    ,'CUM_VLOW'    ,'','very low and cumuliform clouds'                    ),&
e(CT_STRAT_VLOW  ,'STRAT_VLOW'  ,'','very low and stratiform clouds'                    ),&
e(CT_LOW_CUM     ,'LOW_CUM'     ,'','low and cumuliform clouds'                         ),&
e(CT_LOW_STRAT   ,'LOW_STRAT'   ,'','low and stratiform clouds'                         ),&
e(CT_MED_CUM     ,'MED_CUM'     ,'','medium and cumuliform clouds'                      ),&
e(CT_MED_STRAT   ,'MED_STRAT'   ,'','medium and stratiform clouds'                      ),&
e(CT_HI_OP_CUM   ,'HI_OP_CUM'   ,'','high opaque and cumuliform clouds'                 ),&
e(CT_HI_OP_STRAT ,'HI_OP_STRAT' ,'','high opaque and stratiform clouds'                 ),&
e(CT_VHI_OP_CUM  ,'VHI_OP_CUM'  ,'','very high opaque and cumuliform clouds'            ),&
e(CT_VHI_OP_STRAT,'VHI_OP_STRAT','','very high opaque and stratiform clouds'            ),&
e(CT_HI_ST_THIN  ,'HI_ST_THIN'  ,'','high semitransparent thin clouds'                  ),&
e(CT_HI_ST_MEAN  ,'HI_ST_MEAN'  ,'','high semitransparent meanly thick clouds'          ),&
e(CT_HI_ST_THICK ,'HI_ST_THICK' ,'','high semitransparent thick clouds'                 ),&
e(CT_HI_ST_ABOVE ,'HI_ST_ABOVE' ,'','high semitransparent above low or medium clouds'   ),&
e(CT_FRAC        ,'FRAC'        ,'','fractional clouds (sub-pixel water clouds)'        ),&
e(CT_UNDEF       ,'UNDEF'       ,'','undefined (undefined by CMa)'                      )/)

type(t_table) ,pointer :: ct_nwc

!==============================================================================
! surftype, flg_1dvar, flg_cld (consistent with 1dvar)
!==============================================================================

integer ,parameter :: ST_SEA      = 0
integer ,parameter :: ST_ICE      = 1
integer ,parameter :: ST_LAND     = 2
integer ,parameter :: ST_HIGHLAND = 3
integer ,parameter :: ST_MISMATCH = 4

type (e) ,target :: surftype_entries (5) = (/&
e(ST_SEA      ,'SEA'      ,'' ,''),&
e(ST_ICE      ,'ICE'      ,'' ,''),&
e(ST_LAND     ,'LAND'     ,'' ,''),&
e(ST_HIGHLAND ,'HIGHLAND' ,'' ,''),&
e(ST_MISMATCH ,'MISMATCH' ,'' ,'')/)

type(t_table) ,pointer :: surftype

!------------------------------------------------------------------------------

integer ,parameter :: D1_DATA = 0  !
integer ,parameter :: D1_MIN  = 1  ! 1dvar minimisation failed
integer ,parameter :: D1_SUR  = 2  ! wrong surface type
integer ,parameter :: D1_CLD  = 3  ! cloud flag

type (e) ,target :: flg_1dvar_entries (4) = (/      &
e(D1_DATA ,'DATA' ,'' ,'                         '),&
e(D1_MIN  ,'MIN'  ,'' ,'1dvar minimisation failed'),&
e(D1_SUR  ,'SUR'  ,'' ,'wrong surface type       '),&
e(D1_CLD  ,'CLD'  ,'' ,'cloud flag               ')/)

type(t_table) ,pointer :: flg_1dvar

!------------------------------------------------------------------------------

integer ,parameter :: CL_CLEAR     = 0
integer ,parameter :: CL_IR_CLOUDY = 1
integer ,parameter :: CL_MW_CLEAR  = 2
integer ,parameter :: CL_MW_CLOUDY = 3

type (e) ,target :: flg_cld_entries (4) = (/&
e(CL_CLEAR     ,'CLEAR'     ,'' ,''),&
e(CL_IR_CLOUDY ,'IR_CLOUDY' ,'' ,''),&
e(CL_MW_CLEAR  ,'MW_CLEAR'  ,'' ,''),&
e(CL_MW_CLOUDY ,'MW_CLOUDY' ,'' ,'')/)

type(t_table) ,pointer :: flg_cld

!==============================================================================
! level significance
!==============================================================================
!----------
! constants
!----------
integer ,parameter :: LS_SURFACE       = 0 ! surface
integer ,parameter :: LS_STANDARD      = 1 ! standard level
integer ,parameter :: LS_TROPO         = 2 ! tropopause level
integer ,parameter :: LS_MAX           = 3 ! maximum wind level
integer ,parameter :: LS_SIGN          = 4 ! significant level

!--------------
! table entries
!--------------
type (e) ,target :: level_sig_entries (5) = (/                       &
e(LS_SURFACE       ,'SURFACE'       ,'' ,'surface                  '),&
e(LS_STANDARD      ,'STANDARD'      ,'' ,'standard level           '),&
e(LS_TROPO         ,'TROPO'         ,'' ,'tropopause level         '),&
e(LS_MAX           ,'MAX'           ,'' ,'maximum wind level       '),&
e(LS_SIGN          ,'SIGN'          ,'' ,'significant level        ')/)

!------
! table
!------
type(t_table) ,pointer :: level_sig

!==============================================================================
! phase of aircraft flight and roll angle quality
!==============================================================================
!----------
! constants
!----------
!                               0,1 ! Reserved
integer ,parameter :: PH_UNS  = 2   ! Unsteady
integer ,parameter :: PH_LVR  = 3   ! Level flight, routine observation
integer ,parameter :: PH_LVW  = 4   ! Level flight, highest wind encountered
integer ,parameter :: PH_ASC  = 5   ! Ascending
integer ,parameter :: PH_DES  = 6   ! Descending
integer ,parameter :: PH_MIS  = 7   ! Missing value

integer ,parameter :: RA_GOOD = 0   ! Good
integer ,parameter :: RA_BAD  = 1   ! Bad
!                               2   ! Reserved
integer ,parameter :: RA_MIS  = 3   ! Missing value

!--------------
! table entries
!--------------
type (e) ,target :: phase_entries (6) = (/                     &
e(PH_UNS  ,'UNS' ,'','Unsteady                              '),&
e(PH_LVR  ,'LVR' ,'','Level flight, routine observation     '),&
e(PH_LVW  ,'LVW' ,'','Level flight, highest wind encountere '),&
e(PH_ASC  ,'ASC' ,'','Ascending                             '),&
e(PH_DES  ,'DES' ,'','Descending                            '),&
e(PH_MIS  ,'MIS' ,'','Missing                               ')/)

type (e) ,target :: rollangle_entries (3) = (/&
e(RA_GOOD ,'GOOD','','Good                 '),&
e(RA_BAD  ,'BAD' ,'','Bad                  '),&
e(RA_MIS  ,'MIS' ,'','Missing value        ')/)

!-------
! tables
!-------
type(t_table) ,pointer :: phase
type(t_table) ,pointer :: rollangle

!==============================================================================
! retrtype: for AMV satellite derived wind computation method
!==============================================================================
!----------
! constants
!----------
!Wind derived from cloud motion observed in the ..
integer ,parameter :: RT_IR   =   1 ! infrared channel
integer ,parameter :: RT_VIS  =   2 ! visible channel
integer ,parameter :: RT_WV   =   3 ! water vapour channel
!                                 4 ! combination of spectral chans
!                                 5 ! water vapour ch. in clear air
!                                 6 ! ozon channel
!                                 7 ! water vap.ch., cloud or clear not spec.
integer ,parameter :: RT_IR1  = 101 !  8.7 um Meteos8-9, 10.7 um GOES 10-12
integer ,parameter :: RT_IR2  = 201 !  9.7 um Meteos8-9,  3.9 um GOES 10-12
integer ,parameter :: RT_IR3  = 301 ! 10.8 um Meteos8-9
integer ,parameter :: RT_VIS1 = 102 !  0.6 um Meteos8-9,  0.65um GOES 10-12
integer ,parameter :: RT_VIS2 = 302 !  0.75um Meteos8-9
integer ,parameter :: RT_VIS3 = 202 !  0.8 um Meteos8-9
integer ,parameter :: RT_WV1  = 103 !  6.2 um Meteos8-9,  7.4 um GOES 10-12
integer ,parameter :: RT_WV2  = 203 !  7.3 um Meteos8-9,  7.0 um GOES 10-12
integer ,parameter :: RT_WV3  = 303 !                     6.8 um GOES 10-12

!--------------
! table entries
!--------------
type (e) ,target :: retrtype_entries (12) = (/                    &
e(RT_IR  ,'IR     ','','infrared channel                          '),&
e(RT_VIS ,'VIS    ','','visible channel                           '),&
e(RT_WV  ,'WV     ','','water vapour channel                      '),&
e(RT_IR1 ,'IR1    ','',' 8.7 um Meteosat 8-9, ~ 10.7 um GOES 10-12'),&
e(RT_IR2 ,'IR2    ','',' 9.7 um Meteosat 8-9, ~  3.9 um GOES 10-12'),&
e(RT_IR3 ,'IR3    ','','10.8 um Meteosat 8-9                      '),&
e(RT_VIS1,'VIS1   ','',' 0.6 um Meteosat 8-9, ~  0.65um GOES 10-12'),&
e(RT_VIS3,'VIS3   ','',' 0.8 um Meteosat 8-9                      '),&
e(RT_VIS2,'VIS2-HR','',' 0.75um Meteosat 8-9                      '),&
e(RT_WV1 ,'WV1    ','',' 6.2 um Meteosat 8-9,                     '),&
e(RT_WV2 ,'WV2    ','',' 7.3 um Meteosat 8-9,                     '),&
e(RT_WV3 ,'WV3    ','','                ~ 6.8/6.5 um GOES 10,11/12')/)

!-------
! tables
!-------
type(t_table) ,pointer :: retrtype

!==============================================================================
! specification of verification data (runtype, ensmem)
!==============================================================================

integer ,parameter :: VT_FORECAST    = 0 ! forecast
integer ,parameter :: VT_FIRSTGUESS  = 1 ! first guess
integer ,parameter :: VT_PREL_ANA    = 2 ! preliminary analysis
integer ,parameter :: VT_ANALYSIS    = 3 ! analysis
integer ,parameter :: VT_INIT_ANA    = 4 ! initialised analysis
integer ,parameter :: VT_LIN_ANA     = 5 ! linear operator on analysis

type (e) ,target :: runtype_entries (6) = (/                                 &
e(VT_FORECAST   ,'FORECAST   ','','forecast                                 '),&
e(VT_FIRSTGUESS ,'FIRSTGUESS ','','first guess                              '),&
e(VT_PREL_ANA   ,'PREL_ANA   ','','preliminary analysis in observation space'),&
e(VT_ANALYSIS   ,'ANALYSIS   ','','analysis                                 '),&
e(VT_INIT_ANA   ,'INIT_ANA   ','','initialised analysis                     '),&
e(VT_LIN_ANA    ,'LIN_ANA    ','','linear operator on analysis (Y_a)        ')/)

type(t_table) ,pointer :: runtype

!------------------------------------------------------------------------------
integer ,parameter :: VE_ENS_MEAN     =  0 ! ensemble mean
integer ,parameter :: VE_DETERM       = -1 ! derterministic model run
integer ,parameter :: VE_ENS_SPREAD   = -2 ! ensemble spread
integer ,parameter :: VE_BG_ERROR     = -3 ! 3dvar background error
integer ,parameter :: VE_TALAGRAND    = -4 ! Talagrand index
integer ,parameter :: VE_VQC_WEIGHT   = -5 ! variational quality control weight
integer ,parameter :: VE_MEMBER       = -6 ! generic value for ensemble member
integer ,parameter :: VE_ENS_MEAN_OBS = -7 ! ensemble mean in observation space
integer ,parameter :: VE_BIASCOR      = -8 ! bias correction applied

type (e) ,target :: ensmem_entries (9) = (/                               &
e(VE_ENS_MEAN    ,'ENS_MEAN    ','','ensemble mean                     '),&
e(VE_DETERM      ,'DETERM      ','','derterministic model run          '),&
e(VE_ENS_SPREAD  ,'ENS_SPREAD  ','','ensemble spread                   '),&
e(VE_BG_ERROR    ,'BG_ERROR    ','','3dvar background error            '),&
e(VE_TALAGRAND   ,'TALAGRAND   ','','Talagrand index                   '),&
e(VE_VQC_WEIGHT  ,'VQC_WEIGHT  ','','variational quality control weight'),&
e(VE_MEMBER      ,'MEMBER      ','','generic value for ensemble member '),&
e(VE_ENS_MEAN_OBS,'ENS_MEAN_OBS','','ensemble mean in observation space'),&
e(VE_BIASCOR     ,'BIASCOR     ','','bias correction                   ')/)
type(t_table) ,pointer :: ensmem


!==============================================================================
! observation operator flag (observation type dependent)
!==============================================================================
integer ,parameter :: OF_MISSING       = 0 ! missing (not present in old files)
integer ,parameter :: OF_RAD_CLEAR_SKY = 1 ! clear sky radiances 
integer ,parameter :: OF_RAD_CLOUDY    = 2 ! cloudy radiances
integer ,parameter :: OF_BT_CLEAR_SKY  = 3 ! clear sky radiances 
type (e) ,target :: oflag_entries (4) = (/                   &
e(OF_MISSING      ,'MISSING      ','','missing value      '),&
e(OF_RAD_CLEAR_SKY,'RAD_CLEAR_SKY','','clear sky radiances'),&
e(OF_RAD_CLOUDY   ,'RAD_CLOUDY   ','','cloudy radiances   '),&
e(OF_BT_CLEAR_SKY ,'BT_CLEAR_SKY ','','clear sky br.temp. ')/)
type(t_table) ,pointer :: oflag

!==============================================================================
contains
!==============================================================================
  subroutine clean_fdbk_tables
    if (.not.init) return
    init = .false.
    deallocate (status)
    deallocate (flags)
    deallocate (obstype)
    deallocate (codetype)
    deallocate (varno)
    deallocate (satsens)
    deallocate (rsondtype)
    deallocate (trackteqn)
    deallocate (meas_equip)
    deallocate (radiation_corr)
    deallocate (surftype)
    deallocate (flg_1dvar)
    deallocate (surf_char)
    deallocate (flg_cld)
    deallocate (level_sig)
    deallocate (phase)
    deallocate (rollangle)
    deallocate (retrtype)
    deallocate (runtype)
    deallocate (runclass)
    deallocate (ensmem)
    deallocate (ind_cg)
    deallocate (gen_cg)
  end subroutine clean_fdbk_tables
!------------------------------------------------------------------------------
  subroutine init_fdbk_tables (latex)
  logical, intent(in), optional :: latex
  !-------------------------------------------------------------------
  ! join table entries and table meta data (name and caption).
  !   optionally (if 'latex' is given and true)
  !   write LaTeX file with tables for inclusion in the documentation.
  !-------------------------------------------------------------------

    if (init) return
    init = .true.

    call init_table (status                                ,&
                     status_entries                        ,&
                     'status'                              ,&
                     'Observation or report status values' ,&
                     .false.                               ,&
                     latex)

    call init_table (flags                         ,&
                     flag_entries                  ,&
                     'flags'                       ,&
                     'Report or Observation Flags' ,&
                     .true.                        ,&
                     latex)

    call init_table (obstype             ,&
                     obstype_entries     ,&
                     'obstype'           ,&
                     'Observation Types' ,&
                     .false.             ,&
                     latex)

    call init_table (codetype                 ,&
                     codetype_entries         ,&
                     'codetype'               ,&
                     'Observation Code Types' ,&
                     .false.                  ,&
                     latex)

    call init_table (varno                                    ,&
                     varno_entries                            ,&
                     'varno'                                  ,&
                     'Observation and level variable numbers' ,&
                     .false.                                  ,&
                     latex)

    call init_table (satsens                      ,&
                     satsens_entries              ,&
                     'sat-sensor'                 ,&
                     'RTTOV 8 satellite sensors'  ,&
                     .false.                      ,&
                     latex)

    call init_table (rsondtype                                          ,&
                     rsondtype_entries                                  ,&
                     'rsondtype'                                        ,&
'radiosonde type (NRARA), for other values see WMO common code table C2',&
                     .false.                                            ,&
                     latex)

    call init_table (trackteqn                                             ,&
                     trackteqn_entries                                     ,&
                     'trackteqn'                                           ,&
'tracking technique (NSASA), for other values see WMO common code table C7',&
                     .false.                                               ,&
                     latex)

    call init_table (meas_equip                                      ,&
                     meas_equip_entries                              ,&
                     'measequip'                                     ,&
                     'type of measuring equipment used (NA4, 002003)',&
                      .false.                                        ,&
                      latex)

    call init_table (radiation_corr                                         ,&
                     radiation_corr_entries                                 ,&
                     'radiationcorr'                                        ,&
                     'solar and infrared radiation correction (NSR, 002013)',&
                     .false.                                                ,&
                      latex)

    call init_table (surftype                               ,&
                     surftype_entries                       ,&
                     'surftype'                             ,&
                     'surface types consistent with 1d-Var' ,&
                     .true.                                 ,&
                     latex)

    call init_table (flg_1dvar               ,&
                     flg_1dvar_entries       ,&
                     'flg_1dvar'             ,&
                     '1dvar processing flag' ,&
                     .true.                  ,&
                     latex)

    call init_table (surf_char                      ,&
                     surf_char_entries              ,&
                     'surf_char'                    ,&
                     'model surface characteristics',&
                     .true.                         ,&
                     latex)

    call init_table (flg_cld               ,&
                     flg_cld_entries       ,&
                     'flg_cld'             ,&
                     'cloud flag'          ,&
                     .true.                ,&
                     latex)

    call init_table (level_sig                           ,&
                     level_sig_entries                   ,&
                     'level_sig'                         ,&
                     'level significance for TEMP/PILOT.',&
                     .true.                              ,&
                     latex)

    call init_table (phase                                     ,&
                     phase_entries                             ,&
                     'phase'                                   ,&
                     'aircraft phase, radiances fov, gpsro pcd',&
                     .false.                                   ,&
                     latex)

    call init_table (rollangle                        ,&
                     rollangle_entries                ,&
                     'rollangle'                      ,&
                     'aircraft roll angle'            ,&
                     .false.                          ,&
                     latex)

    call init_table (retrtype                                    ,&
                     retrtype_entries                            ,&
                     'retrtype'                                  ,&
                     'satellite derived wind computation method' ,&
                     .false.                                     ,&
                     latex)

    call init_table (runtype                    ,&
                     runtype_entries            ,&
                     'runtype'                  ,&
                     'type of verification run' ,&
                     .false.                    ,&
                     latex)

    call init_table (runclass                    ,&
                     runclass_entries            ,&
                     'runclass'                  ,&
                     'class of verification run' ,&
                     .false.                     ,&
                     latex)

    call init_table (ensmem                                   ,&
                     ensmem_entries                           ,&
                     'ensmem'                                 ,&
                     'specification of the verification data' ,&
                     .false.                                  ,&
                     latex)

    call init_table (ind_cg                                 ,&
                     ind_cg_entries                         ,&
                     'ind_cg'                               ,&
                     'individual cloud group bit positions ',&
                     .false.                                ,&
                     latex)

    call init_table (gen_cg                             ,&
                     gen_cg_entries                     ,&
                     'gen_cg'                           ,&
                     'general cloud group bit positions',&
                     .false.                            ,&
                     latex)

  end subroutine init_fdbk_tables
!==============================================================================
end module mo_fdbk_tables
