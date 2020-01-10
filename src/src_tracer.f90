!+ Methods for general treatment of tracers
!------------------------------------------------------------------------------
#ifdef MESSY
#define STATIC_FIELDS
#define IJNK_STORAGE
#endif

MODULE src_tracer

!------------------------------------------------------------------------------
!
! Description:
!
!   This modules performs several actions related to the tracers ("passive actions" only).
!   The following routines are available:
!
!     NAME                 DESCRIPTION                                ERROR CODE RANGE
!
!     trcr_init            initialize the tracer module               4000-4049
!     trcr_new             define a new tracer                        4050-4099
!     trcr_alloc           allocate the memory for the tracers        4100-4149
!     trcr_setup_vartab    define the setup_vartab for tracers        4150-4199
!     trcr_print           print the tracer list                      4200-4249
!     trcr_get             retrieve tracer data                       4250-4299
!     trcr_get_block       retrieve block of tracer data              4250-4289
!     trcr_swap            swap timelevels for tracers                4290-4299
!     trcr_cleanup         cleanup the tracer module                  4300-4349
!     trcr_calc            calculate nnew as nnow + dt * tend for MESSy
!     trcr_get_index       retrieve tracer index                      4350-4359
!     trcr_get_ntrcr       get number of tracers                          -
!     trcr_meta_define     define tracer information                     5000
!     trcr_meta_get        retrieve tracer information                   5000
!     trcr_meta_set        store tracer information                      5000
!     trcr_errorstr        get string for error number                    -
!
! Current Code Owner: MeteoSwiss, Oliver Fuhrer
!  phone:  +41 58 460 9359
!  fax:    +41 58 460 9278
!  email:  oliver.fuhrer@meteoswiss.ch
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer
!  Initial release
! V4_26        2012/12/06 Anne Roches
!  Changes and technical adaptations to the tracer handling
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  Adaptations for using this tracer module also with MESSy
! V4_28        2013/07/12 Oliver Fuhrer
!  Correction of a tracer-index
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
! V5_1         2014-11-28 Oliver Fuhrer, Anne Roches, Xavier Lapillonne
!  Replaced ireals by wp (working precision) (OF)
!  Removed hacks for the tracers (AR)
!  New pointer output option added: trcr_get_block (XL)
! V5_2         2015-05-21 Ulrich Schaettler
!  Reorganized numbering of MPI datatypes ID
! V5_3a        2015-11-24 Ulrich Schaettler
!  Included new routines trcr_acc_update_device, trcr_acc_update_host for GPUs
! V5_4         2016-03-10 Oliver Fuhrer
!  Updated code owner information
! V5_4a        2016-05-10 Pascal Spoerri, Xavier Lapillonne
!  Integrated optional in trcr_get_byindex argument to obtain tracer boundary without time level
!  (merged from pompa).
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:
!
USE data_parameters, ONLY :   &
  wp              , & ! KIND-type parameter for real variables
  sp              , & ! KIND-type parameter for real variables (single precision)
  dp              , & ! KIND-type parameter for real variables (double precision)
  iintegers           ! KIND-type parameter for standard integer variables

USE data_parallel,   ONLY :   &
  my_cart_id          ! rank of this subdomain in the cartesian communicator

USE data_runcontrol, ONLY :   &
  l2tls ,           & ! time integration by two timelevel RK-scheme (.TRUE.)
  itype_conv,       & ! type of convection parameterization
  lconv,            & ! forecast with convection
  ltur,             & ! forecast with turbulent mixing
  l3dturb,          & ! 3D turbulent mixing
  lhordiff,         & ! running with horizontal diffusion
  hd_corr_trcr_bd,  & ! correction factor for horizontal diffusion flux of tracers in boundary zone
  hd_corr_trcr_in,  & ! correction factor for horizontal diffusion flux of tracers in domain
  l2dim,            & ! lartif_data=.TRUE.:  2dimensional model version
  nnow,             & ! corresponds to ntstep
  nnew,             & ! corresponds to ntstep + 1
  nold,             & ! corresponds to ntstep - 1
  lspubc ,          & ! with Rayleigh damping in the upper levels
  itype_spubc,      & ! type of Rayleigh damping in the upper levels
  idbg_level,       & ! to control the verbosity of debug output
  lprintdeb_all       ! .TRUE.:  all tasks print debug output

USE data_modelconfig, ONLY :  &
  ie,               & ! number of grid points in zonal direction
  je,               & ! number of grid points in meridional direction
  ke,               & ! number of grid points in vertical direction
  dt                  ! long time-step

USE data_io,         ONLY :   &
  ar_des,           & ! structure for LM variable table
  var ,             & ! array for LM variable table
  max_gribrep,      & ! maximum number of allowed repetition of combination
                      ! grib number - grib table number
  num_gribtabs        ! number of GRIB tables used in LM variable table

USE data_tracer

!------------------------------------------------------------------------------

USE data_tracer_metadata,   ONLY :   &
  I_MAX_STRLEN        ! maximum string length for storage

USE src_tracer_metadata,    ONLY :   &
  metadata_init,    & ! inititalize metadata storage
  metadata_modify,  & ! modify metadata storage
  metadata_finish,  & ! cleanup metadata storage
  metadata_print,   & ! print metadata storage
  metadata_define,  & ! define metadata entry
  metadata_set,     & ! store metadata entry
  metadata_get,     & ! retrieve metadata entry
  metadata_delete,  & ! remove metadata entry
  metadata_error      ! error message

#ifdef MESSY
USE data_fields,        ONLY: qv_s
USE data_modelconfig,   ONLY: idt_qv

! MESSy/BMIL
USE messy_main_blather_bi,    ONLY: info_bi, error_bi
USE messy_main_tracer_mem_bi, ONLY: GPTRSTR, xt_tl, xtte, xt_bd &
                                  , xt, xtm1, xtf, ntrac_gp, ti_gp
USE messy_main_tracer_bi,     ONLY: tracer_halt
! MESSy/SMCL
USE messy_main_tracer,        ONLY: new_tracer, set_tracer, get_tracer     &
                                  , get_tracer_set_info, ON, OFF           &
                                  , full2base_sub, print_tracer_set        &
                                  , AMOUNTFRACTION, AIR, TR_NEXIST         &
                                  , I_GRIBTAB, I_GRIBPARAM, I_MMD_INIT     &
                                  , I_ADVECT,  I_CONVECT, I_VDIFF, I_HDIFF &
                                  , I_INITIAL, I_RELAX,   I_DAMP,  I_LBC   &
                                  , MAX_CASK_I, MAX_CASK_R, MAX_CASK_S
USE messy_main_constants_mem, ONLY: STRLEN_ULONG, STRLEN_LONG, STRLEN_MEDIUM
#endif

!------------------------------------------------------------------------------

IMPLICIT NONE

!------------------------------------------------------------------------------

PRIVATE

! Local variables

INTEGER(KIND=iintegers)                ::       &
  izerrind = -999_iintegers,        &  ! index of the tracer for which an error occurs
  izdebug                              ! for local debug output

INTEGER(KIND=iintegers) :: &
  i_status = T_STAT_START              ! status of tracer module (the tracer module 
                                       ! is a state machine and this integer
                                       ! defines the current state (see all
                                       ! T_STAT_XXX constants in data_tracer.f90
                                       ! for a list of all possible states)
                                       
INTEGER(KIND=iintegers) :: &
  n_trcr    = 0_iintegers  ,        &  ! number of tracers
  n_bd      = 0_iintegers  ,        &  ! number of tracers having boundary data
  n_tlev    = 0_iintegers              ! number of timelevels saved

CHARACTER(LEN=10)       :: yzerrind=''! string corresponding to izerrind

#ifdef MESSY
TYPE t_meta_relate_name
   CHARACTER(LEN=STRLEN_LONG)   :: cosmo_str = ''
   ! index of META INFORMATION cosmo
   INTEGER(KIND=iintegers)      :: cosmo_idx = 0
   ! index of META INFORMATION MESSy
   INTEGER(KIND=iintegers)      :: messy_idx = 0
END type t_meta_relate_name
INTEGER, PARAMETER :: NMAX_META = 9
TYPE(t_meta_relate_name), DIMENSION(NMAX_META), SAVE :: meta_list
#endif

#ifdef STATIC_FIELDS
INTEGER(KIND=iintegers), SAVE :: t_nnew, t_nnow, t_nold
#endif

! Declare public entities

PUBLIC :: &
  trcr_init,                        &
  trcr_new,                         &
  trcr_alloc,                       &
  trcr_setup_vartab,                & 
  trcr_print,                       &
  trcr_get,                         &
  trcr_get_block,                   &
  trcr_get_index,                   &
  trcr_swap,                        &
  trcr_cleanup,                     &
  trcr_calc,                        &
  trcr_get_ntrcr,                   &
  trcr_meta_define,                 &
  trcr_meta_get,                    &
  trcr_meta_set,                    &
  trcr_acc_update_device,           &
  trcr_acc_update_host,             &
  trcr_errorstr

! Interfaces for overloaded subroutines

INTERFACE trcr_get
  MODULE PROCEDURE          &
    trcr_get_byname,        &
    trcr_get_byindex
END INTERFACE

INTERFACE trcr_meta_define
  MODULE PROCEDURE trcr_meta_def_real4,    trcr_meta_def_real8,               &
                   trcr_meta_def_integer,  trcr_meta_def_string,              &
                   trcr_meta_def_logical,  trcr_meta_def_pointer2,            &
                   trcr_meta_def_pointer3, trcr_meta_def_pointer4
END INTERFACE

INTERFACE trcr_meta_set
  MODULE PROCEDURE trcr_meta_set_r4_byindind,   trcr_meta_set_r8_byindind,    &
                 trcr_meta_set_int_byindind,    trcr_meta_set_str_byindind,   &
                 trcr_meta_set_log_byindind,    trcr_meta_set_ptr2_byindind,  &
                 trcr_meta_set_ptr3_byindind,   trcr_meta_set_ptr4_byindind,  &
                 trcr_meta_set_r4_byindname,    trcr_meta_set_r8_byindname,   &
                 trcr_meta_set_int_byindname,   trcr_meta_set_str_byindname,  &
                 trcr_meta_set_log_byindname,   trcr_meta_set_ptr2_byindname, &
                 trcr_meta_set_ptr3_byindname,  trcr_meta_set_ptr4_byindname, &
                 trcr_meta_set_r4_bynameind,    trcr_meta_set_r8_bynameind,   &
                 trcr_meta_set_int_bynameind,   trcr_meta_set_str_bynameind,  &
                 trcr_meta_set_log_bynameind,   trcr_meta_set_ptr2_bynameind, &
                 trcr_meta_set_ptr3_bynameind,  trcr_meta_set_ptr4_bynameind, &
                 trcr_meta_set_r4_bynamename,   trcr_meta_set_r8_bynamename,  &
                 trcr_meta_set_int_bynamename,  trcr_meta_set_str_bynamename, &
                 trcr_meta_set_log_bynamename,  trcr_meta_set_ptr2_bynamename,&
                 trcr_meta_set_ptr3_bynamename, trcr_meta_set_ptr4_bynamename,&
                 trcr_meta_set_r4v_byindex,     trcr_meta_set_r8v_byindex,    &
                 trcr_meta_set_intv_byindex,    trcr_meta_set_strv_byindex,   &
                 trcr_meta_set_logv_byindex,                                  &
                 trcr_meta_set_r4v_byname,      trcr_meta_set_r8v_byname,     &
                 trcr_meta_set_intv_byname,     trcr_meta_set_strv_byname,    &
                 trcr_meta_set_logv_byname
END INTERFACE

INTERFACE trcr_meta_get
  MODULE PROCEDURE trcr_meta_get_r4_byindind,   trcr_meta_get_r8_byindind,    &
                 trcr_meta_get_int_byindind,    trcr_meta_get_str_byindind,   &
                 trcr_meta_get_log_byindind,    trcr_meta_get_ptr2_byindind,  &
                 trcr_meta_get_ptr3_byindind,   trcr_meta_get_ptr4_byindind,  &
                 trcr_meta_get_r4_byindname,    trcr_meta_get_r8_byindname,   &
                 trcr_meta_get_int_byindname,   trcr_meta_get_str_byindname,  &
                 trcr_meta_get_log_byindname,   trcr_meta_get_ptr2_byindname, &
                 trcr_meta_get_ptr3_byindname,  trcr_meta_get_ptr4_byindname, &
                 trcr_meta_get_r4_bynameind,    trcr_meta_get_r8_bynameind,   &
                 trcr_meta_get_int_bynameind,   trcr_meta_get_str_bynameind,  &
                 trcr_meta_get_log_bynameind,   trcr_meta_get_ptr2_bynameind, &
                 trcr_meta_get_ptr3_bynameind,  trcr_meta_get_ptr4_bynameind, &
                 trcr_meta_get_r4_bynamename,   trcr_meta_get_r8_bynamename,  &
                 trcr_meta_get_int_bynamename,  trcr_meta_get_str_bynamename, &
                 trcr_meta_get_log_bynamename,  trcr_meta_get_ptr2_bynamename,&
                 trcr_meta_get_ptr3_bynamename, trcr_meta_get_ptr4_bynamename,&
                 trcr_meta_get_r4v_byindex,     trcr_meta_get_r8v_byindex,    &
                 trcr_meta_get_intv_byindex,    trcr_meta_get_strv_byindex,   &
                 trcr_meta_get_logv_byindex,                                  &
                 trcr_meta_get_r4v_byname,      trcr_meta_get_r8v_byname,     &
                 trcr_meta_get_intv_byname,     trcr_meta_get_strv_byname,    &
                 trcr_meta_get_logv_byname
END INTERFACE

!==============================================================================
! Module Procedures in src_tracer
!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure "trcr_init" in "src_tracer" to initialize a tracer struct.
!------------------------------------------------------------------------------

SUBROUTINE trcr_init ( ierr )
!------------------------------------------------------------------------------
!
! Description:
!
! This routine initializes a tracer structure.
!
!
! Method:
!
! All elements of the derived datatype trcr_type are initialized with a value
! or a parameter value.
!------------------------------------------------------------------------------

 IMPLICIT NONE

! Subroutine arguments:
! --------------------

!error status
 INTEGER                 ,INTENT(OUT)         :: ierr 

! Local parameters:
! ----------------

 INTEGER(KIND=iintegers)   :: i, izerr

!------------------------------------------------------------------------------
! Begin Subroutine trcr_init
!------------------------------------------------------------------------------
  
  ! initialize error status
  ierr = 0_iintegers

  ! check workflow
  ierr = 4000_iintegers
  IF ( i_status /= T_STAT_START ) RETURN

  ! init debugging message level
  IF (lprintdeb_all) THEN
    izdebug = idbg_level
  ELSE
    IF (my_cart_id == 0_iintegers) THEN
      izdebug = idbg_level
    ELSE
      izdebug = 0_iintegers
    ENDIF
  ENDIF

  ! index for the variables holding the tracer data
  trcr(:)%idx_trcr = T_MISSING
  trcr(:)%idx_bd   = T_MISSING
  trcr(:)%igribrep = 1_iintegers

#ifdef STATIC_FIELDS
#ifndef MESSY
  ! setup fixed time indices
  t_nnew = nnew
  t_nnow = nnow
  IF (.NOT. l2tls) THEN
    t_nold = nold
  ENDIF
#else
  t_nnew = 1
  IF (l2tls) THEN
     t_nnow = 3
  ELSE
     t_nnow = 6
     t_nold = 3
  ENDIF
#endif
#endif

#ifndef MESSY
  ! init tracer meta-data storage
  ierr = 5000_iintegers
  CALL metadata_init(mtrcr, imaxbuf=n_trcr_max, imaxkey=n_mkey_max, &
         ibuflen=n_mbuf_max)
  IF (metadata_error() /= 0_iintegers) RETURN

  ! metadata
  ierr = 5000_iintegers
  CALL trcr_meta_define( izerr, T_NAME_KEY, T_NAME_DEFAULT, &
    iidx=T_NAME_ID, lprotect=.TRUE. )
  IF (izerr /= 0_iintegers) RETURN
  CALL trcr_meta_define( izerr, T_PARENT_KEY, T_PARENT_DEFAULT, &
    iidx=T_PARENT_ID, lprotect=.TRUE. )
  IF (izerr /= 0_iintegers) RETURN
  CALL trcr_meta_define( izerr, T_UNITS_KEY, T_UNITS_DEFAULT, &
    iidx=T_UNITS_ID, lprotect=.TRUE. )
  IF (izerr /= 0_iintegers) RETURN
  CALL trcr_meta_define( izerr, T_NCSTDNAME_KEY, T_NCSTDNAME_DEFAULT, &
    iidx=T_NCSTDNAME_ID, lprotect=.TRUE. )
  IF (izerr /= 0_iintegers) RETURN
  CALL trcr_meta_define( izerr, T_NCLONGNAME_KEY, T_NCLONGNAME_DEFAULT, &
    iidx=T_NCLONGNAME_ID, lprotect=.TRUE. )
  IF (izerr /= 0_iintegers) RETURN
  CALL trcr_meta_define( izerr, T_GRBPARAM_KEY, T_GRBPARAM_DEFAULT, &
    iidx=T_GRBPARAM_ID, lprotect=.TRUE. )
  IF (izerr /= 0_iintegers) RETURN
  CALL trcr_meta_define( izerr, T_GRBTABLE_KEY, T_GRBTABLE_DEFAULT, &
    iidx=T_GRBTABLE_ID, lprotect=.TRUE. )
  IF (izerr /= 0_iintegers) RETURN

  ! control switches
  ierr = 5000_iintegers
  CALL trcr_meta_define( izerr, T_ADV_KEY, T_ADV_DEFAULT, &
    iidx=T_ADV_ID, lprotect=.TRUE. )
  IF (izerr /= 0_iintegers) RETURN
  CALL trcr_meta_define( izerr, T_DIFF_KEY, T_DIFF_DEFAULT, &
    iidx=T_DIFF_ID, lprotect=.TRUE. )
  IF (izerr /= 0_iintegers) RETURN
  CALL trcr_meta_define( izerr, T_TURB_KEY, T_TURB_DEFAULT, &
    iidx=T_TURB_ID, lprotect=.TRUE. )
  IF (izerr /= 0_iintegers) RETURN
  CALL trcr_meta_define( izerr, T_CONV_KEY, T_CONV_DEFAULT, &
    iidx=T_CONV_ID, lprotect=.TRUE. )
  IF (izerr /= 0_iintegers) RETURN
  CALL trcr_meta_define( izerr, T_INI_KEY, T_INI_DEFAULT, &
    iidx=T_INI_ID, lprotect=.TRUE. )
  IF (izerr /= 0_iintegers) RETURN
  CALL trcr_meta_define( izerr, T_LBC_KEY, T_LBC_DEFAULT, &
    iidx=T_LBC_ID, lprotect=.TRUE. )
  IF (izerr /= 0_iintegers) RETURN
  CALL trcr_meta_define( izerr, T_BBC_KEY, T_BBC_DEFAULT, &
    iidx=T_BBC_ID, lprotect=.TRUE. )
  IF (izerr /= 0_iintegers) RETURN
  CALL trcr_meta_define( izerr, T_RELAX_KEY, T_RELAX_DEFAULT, &
    iidx=T_RELAX_ID, lprotect=.TRUE. )
  IF (izerr /= 0_iintegers) RETURN
  CALL trcr_meta_define( izerr, T_DAMP_KEY, T_DAMP_DEFAULT, &
    iidx=T_DAMP_ID, lprotect=.TRUE. )
  IF (izerr /= 0_iintegers) RETURN
  CALL trcr_meta_define( izerr, T_CLP_KEY, T_CLP_DEFAULT, &
    iidx=T_CLP_ID, lprotect=.TRUE. )
  IF (izerr /= 0_iintegers) RETURN
#else
  ! BUILD META ASSOCIATION LIST:

  meta_list(1)%cosmo_str = 'SP_ADV_LF'
  meta_list(1)%cosmo_idx = -1
  meta_list(1)%messy_idx = I_ADVECT

  meta_list(2)%cosmo_str = T_ADV_KEY
  meta_list(2)%cosmo_idx = T_ADV_ID
  meta_list(2)%messy_idx = I_ADVECT

  meta_list(3)%cosmo_str = T_DIFF_KEY
  meta_list(3)%cosmo_idx = T_DIFF_ID
  meta_list(3)%messy_idx = I_HDIFF

  meta_list(4)%cosmo_str = T_TURB_KEY
  meta_list(4)%cosmo_idx = T_TURB_ID
  meta_list(4)%messy_idx = I_VDIFF

  meta_list(5)%cosmo_str = T_CONV_KEY
  meta_list(5)%cosmo_idx = T_CONV_ID
  meta_list(5)%messy_idx = I_CONVECT

  meta_list(6)%cosmo_str = T_INI_KEY
  meta_list(6)%cosmo_idx = T_INI_ID
  meta_list(6)%messy_idx = I_INITIAL

  meta_list(7)%cosmo_str = T_LBC_KEY
  meta_list(7)%cosmo_idx = T_LBC_ID
  meta_list(7)%messy_idx = I_LBC

  meta_list(8)%cosmo_str = T_RELAX_KEY
  meta_list(8)%cosmo_idx = T_RELAX_ID
  meta_list(8)%messy_idx = I_RELAX

  meta_list(9)%cosmo_str = T_DAMP_KEY
  meta_list(9)%cosmo_idx = T_DAMP_ID
  meta_list(9)%messy_idx = I_DAMP

#endif

  ! update tracer module status
  i_status = T_STAT_DEFINE

  ! no errors
  ierr = 0_iintegers

!------------------------------------------------------------------------------
! End of module procedure trcr_init
!------------------------------------------------------------------------------

END SUBROUTINE trcr_init

!==============================================================================


!==============================================================================
!+ Module procedure "trcr_new" in "src_tracer" to define a new tracer.
!------------------------------------------------------------------------------

SUBROUTINE trcr_new ( ierr, yshort_name, igribparam, igribtable, yparent,     &
                      yunits, ystandard_name, ylong_name,                     &
                      itype_adv, itype_diff, itype_turbmix, itype_passconv,   &
                      itype_ini, itype_lbc, itype_bbc, itype_relax,           &
                      itype_damp, itype_clip, idx_trcr )
!------------------------------------------------------------------------------
!
! Description:
!
! This routine defines a new tracer (and its associated meta-data).
! This routine has to be called each time that a new tracer is required
! anywhere in the code. It should be called from the "init" section of
! each module which uses tracers and _before_ memory allocation. Here,
! we only prepare the necessary data structures.
! In addition, some basic consistency checks are done (e.g. name
! already used, grib table/parameter already occupied, etc.)
!
!
! Method:
!
!
!------------------------------------------------------------------------------

 IMPLICIT NONE

! Subroutine arguments:
! --------------------
 INTEGER                ,INTENT(OUT)         :: ierr           !error status        
 CHARACTER(LEN=*)       ,INTENT(IN)          :: yshort_name    !name of tracer    
 INTEGER(KIND=iintegers),INTENT(IN)          :: igribparam     !GRIB parameter number (= iloc2)
 INTEGER(KIND=iintegers),INTENT(IN)          :: igribtable     !GRIB table number     (= iloc3)
 CHARACTER(LEN=*)       ,INTENT(IN)          :: yparent        !name of parent module (=calling routine)                          
 CHARACTER(LEN=*)       ,INTENT(IN)          :: yunits         !units of tracer (needed for NetCDF)                 
 CHARACTER(LEN=*)       ,INTENT(IN) ,OPTIONAL:: ystandard_name !standard name of tracer (needed for NetCDF)   
 CHARACTER(LEN=*)       ,INTENT(IN) ,OPTIONAL:: ylong_name     !long name of tracer (needed for NetCDF)
 INTEGER(KIND=iintegers),INTENT(IN) ,OPTIONAL:: itype_adv      !advection specification
 INTEGER(KIND=iintegers),INTENT(IN) ,OPTIONAL:: itype_diff     !artificial hyperdiffusion specification
 INTEGER(KIND=iintegers),INTENT(IN) ,OPTIONAL:: itype_turbmix  !turbulent mixing specification  
 INTEGER(KIND=iintegers),INTENT(IN) ,OPTIONAL:: itype_passconv !convective passvie transport specification  
 INTEGER(KIND=iintegers),INTENT(IN) ,OPTIONAL:: itype_ini      !type of initial conditions (IC)                                                
 INTEGER(KIND=iintegers),INTENT(IN) ,OPTIONAL:: itype_lbc      !type of lat. boundary conditions (BC)              
 INTEGER(KIND=iintegers),INTENT(IN) ,OPTIONAL:: itype_bbc      !type of bottom boundary conditions (BC)                                    
 INTEGER(KIND=iintegers),INTENT(IN) ,OPTIONAL:: itype_relax    !type of relaxation
 INTEGER(KIND=iintegers),INTENT(IN) ,OPTIONAL:: itype_damp     !type of Rayleigh damping
 INTEGER(KIND=iintegers),INTENT(IN) ,OPTIONAL:: itype_clip     !type of clipping           
 INTEGER(KIND=iintegers),INTENT(OUT),OPTIONAL:: idx_trcr       !index for tracer data                                                     

!
! Local parameters:
! ----------------

 INTEGER(KIND=iintegers) :: &          !temporary indices for new 
   i_trcr, &                            ! tracer,
   i_bd                                 ! boundary,
                                                 
 INTEGER(KIND=iintegers) :: &
   i,      &  !index for looping inside tracer struct
   itmp, itmp1, itmp2, itmp3

 LOGICAL :: &
   ltmp

 CHARACTER(LEN=ilen_sn) :: &
   ytmp(n_trcr+1)
                                                 
#ifdef MESSY
  ! LOCAL Definitions MESSY
  CHARACTER(LEN=*), PARAMETER :: substr = 'trcr_new'
  INTEGER                     :: status
  CHARACTER(STRLEN_ULONG)     :: standard_name = ''

  INTEGER                     :: idt
#endif

!----------- End of header ----------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine trcr_new
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  Section 1: Initialization and preparation
!------------------------------------------------------------------------------

  ! initialize error status
  ierr     = 0_iintegers
  yzerrind = yshort_name

  ! check workflow
  ierr = 4050_iintegers
  IF ( i_status /= T_STAT_DEFINE ) RETURN

#ifndef MESSY
  ! increase number of tracers
  ierr = 4055_iintegers
  IF (n_trcr+1_iintegers > n_trcr_max) RETURN
  i_trcr = n_trcr + 1_iintegers

!------------------------------------------------------------------------------
!  Section 2: First checks
!------------------------------------------------------------------------------

  ! check length of strings (mandatory and optional)
  ierr = 4060_iintegers
  IF (LEN_TRIM(yshort_name) > ilen_sn) RETURN           !tracer name
  ierr = 4062_iintegers
  IF (LEN_TRIM(yparent) > ilen_ln) RETURN               !parent name
  ierr = 4064_iintegers
  IF (LEN_TRIM(yunits) > ilen_un) RETURN                !unit name
  IF (PRESENT(ystandard_name)) THEN                  
    ierr = 4066_iintegers
    IF (LEN_TRIM(ystandard_name) > ilen_ln) RETURN      !standard name
  ENDIF
  IF (PRESENT(ylong_name)) THEN
    ierr = 4067_iintegers
    IF (LEN_TRIM(ylong_name) > ilen_ln) RETURN          !long name
  ENDIF
 
!------------------------------------------------------------------------------
!  Section 3: Value assignments
!------------------------------------------------------------------------------

  ! add metadata buffer for this tracer
  CALL metadata_modify(mtrcr, inumbuf=i_trcr)

  ! metadata (mandatory)
  ierr = 5000_iintegers
  CALL trcr_meta_set( itmp, i_trcr, T_NAME_ID, TRIM(yshort_name) )
  IF (itmp /= 0_iintegers) RETURN
  CALL trcr_meta_set( itmp, i_trcr, T_PARENT_ID, TRIM(yparent) )
  IF (itmp /= 0_iintegers) RETURN
  CALL trcr_meta_set( itmp, i_trcr, T_UNITS_ID, TRIM(yunits) )
  IF (itmp /= 0_iintegers) RETURN
  CALL trcr_meta_set( itmp, i_trcr, T_GRBPARAM_ID, igribparam )
  IF (itmp /= 0_iintegers) RETURN
  CALL trcr_meta_set( itmp, i_trcr, T_GRBTABLE_ID, igribtable )
  IF (itmp /= 0_iintegers) RETURN

  ! metadata (optional)
  ierr = 5000_iintegers
  IF (PRESENT(ystandard_name)) THEN
    CALL trcr_meta_set( itmp, i_trcr, T_NCSTDNAME_ID, TRIM(ystandard_name) )
    IF (itmp /= 0_iintegers) RETURN
  ENDIF
  IF (PRESENT(ylong_name)) THEN
    CALL trcr_meta_set( itmp, i_trcr, T_NCLONGNAME_ID, TRIM(ylong_name) )
    IF (itmp /= 0_iintegers) RETURN
  ENDIF

  ! control switches (mandatory)

  ! control switches (optional)
  ierr = 5000_iintegers
  IF (PRESENT(itype_adv)) THEN
    CALL trcr_meta_set( itmp, i_trcr, T_ADV_ID, itype_adv )
    IF (itmp /= 0_iintegers) RETURN
  ENDIF
  IF (PRESENT(itype_diff)) THEN
    CALL trcr_meta_set( itmp, i_trcr, T_DIFF_ID, itype_diff )
    IF (itmp /= 0_iintegers) RETURN
  ENDIF
  IF (PRESENT(itype_turbmix)) THEN
    CALL trcr_meta_set( itmp, i_trcr, T_TURB_ID, itype_turbmix )
    IF (itmp /= 0_iintegers) RETURN
  ENDIF
  IF (PRESENT(itype_passconv)) THEN
    CALL trcr_meta_set( itmp, i_trcr, T_CONV_ID, itype_passconv )
    IF (itmp /= 0_iintegers) RETURN
  ENDIF
  IF (PRESENT(itype_ini)) THEN
    CALL trcr_meta_set( itmp, i_trcr, T_INI_ID, itype_ini )
    IF (itmp /= 0_iintegers) RETURN
  ENDIF
  IF (PRESENT(itype_lbc)) THEN
    CALL trcr_meta_set( itmp, i_trcr, T_LBC_ID, itype_lbc )
    IF (itmp /= 0_iintegers) RETURN
  ENDIF
  IF (PRESENT(itype_bbc)) THEN
    CALL trcr_meta_set( itmp, i_trcr, T_BBC_ID, itype_bbc )
    IF (itmp /= 0_iintegers) RETURN
  ENDIF
  IF (PRESENT(itype_relax)) THEN
    CALL trcr_meta_set( itmp, i_trcr, T_RELAX_ID, itype_relax )
    IF (itmp /= 0_iintegers) RETURN
  ENDIF
  IF (PRESENT(itype_damp)) THEN
    CALL trcr_meta_set( itmp, i_trcr, T_DAMP_ID, itype_damp )
    IF (itmp /= 0_iintegers) RETURN
  ENDIF
  IF (PRESENT(itype_clip)) THEN
    CALL trcr_meta_set( itmp, i_trcr, T_CLP_ID, itype_clip )
    IF (itmp /= 0_iintegers) RETURN
  ENDIF

  ! it is not possible to have convective transport if lconv is set to 
  ! FALSE or if itype_conv is not equal to 0.
  ierr = 5000_iintegers
  CALL trcr_meta_get( itmp, i_trcr, T_CONV_ID, itmp1 )
  IF (itmp /= 0_iintegers) RETURN
  IF ( itmp1 == T_CONV_ON .AND. (.NOT. (lconv .AND. itype_conv==0)) ) THEN
    IF (izdebug > 1_iintegers) THEN
      PRINT *,  ''
      PRINT *,  '  *** WARNING: Passive transport by convection specified but lconv=F or itype_conv/=0 ***'
      PRINT *,  '  (Your specification: itype_passconv is not respected for '//TRIM(yshort_name)//')'
      PRINT *,  ''
    ENDIF
    ierr = 5000_iintegers
    CALL trcr_meta_set( itmp, i_trcr, T_CONV_ID, T_CONV_OFF )
    IF (itmp /= 0_iintegers) RETURN
  ENDIF

  ! it is not possible to have horizontal diffusion if lhordiff is set to

  ierr = 5000_iintegers
  CALL trcr_meta_get( itmp, i_trcr, T_DIFF_ID, itmp1 )
  IF (itmp /= 0_iintegers) RETURN
  IF ( itmp1 == T_DIFF_ON .AND. ((.NOT. lhordiff ) .OR. &

       (hd_corr_trcr_bd == 0 .AND. hd_corr_trcr_in == 0)) ) THEN

    IF (izdebug > 1_iintegers) THEN
      PRINT *,  ''
      PRINT *,  '  *** WARNING: Horizontal diffusion specified but lhordiff=F or hd_corr_trcr_xx=0***'
      PRINT *,  '  (Your specification: itype_diff is not respected for '//TRIM(yshort_name)//')'
      PRINT *,  ''
    ENDIF
    ierr = 5000_iintegers
    CALL trcr_meta_set( itmp, i_trcr, T_DIFF_ID, T_DIFF_OFF )
    IF (itmp /= 0_iintegers) RETURN
  ENDIF

  ! it is not possible to have boundary relaxation or Rayleigh damping
  ! in case of zerograd boundary conditions or constant boundary conditions
  ierr = 5000_iintegers
  CALL trcr_meta_get( itmp, i_trcr, T_LBC_ID, itmp1 )
  IF (itmp /= 0_iintegers) RETURN
  ierr = 5000_iintegers
  CALL trcr_meta_get( itmp, i_trcr, T_RELAX_ID, itmp2 )
  IF (itmp /= 0_iintegers) RETURN
  ierr = 5000_iintegers
  CALL trcr_meta_get( itmp, i_trcr, T_DAMP_ID, itmp3 )
  IF (itmp /= 0_iintegers) RETURN

  IF ( (itmp1 == T_LBC_ZEROGRAD .OR. itmp1==T_LBC_CST) .AND.                  &
       (itmp2 == T_RELAX_FULL .OR. itmp2 == T_RELAX_INFLOW                    &
     .OR. (itmp3 == T_DAMP_ON .AND. itype_spubc == 1)) ) THEN
    IF (izdebug > 1_iintegers) THEN
      PRINT *,  ''
      PRINT *,  '  *** WARNING: Relaxation and damping not meaningful for this BD ***'
      PRINT *,  '  (Your specification: itype_relax or/and itype_damp are not respected for '//TRIM(yshort_name)//')'
      PRINT *,  ''
    ENDIF
    ierr = 5000_iintegers
    CALL trcr_meta_set( itmp, i_trcr, T_RELAX_ID, T_RELAX_OFF )
    IF (itmp /= 0_iintegers) RETURN
    ierr = 5000_iintegers
    CALL trcr_meta_set( itmp, i_trcr, T_DAMP_ID, T_DAMP_OFF )
    IF (itmp /= 0_iintegers) RETURN
  ENDIF

  ! it is not possible to have 3D turbulent mixing if l3dturb is set to 
  ! FALSE 
  ierr = 5000_iintegers
  CALL trcr_meta_get( itmp, i_trcr, T_TURB_ID, itmp1 )
  IF (itmp /= 0_iintegers) RETURN
  IF ( itmp1 == T_TURB_3D .AND. .NOT. l3dturb ) THEN
    IF (izdebug > 1_iintegers) THEN
      PRINT *,  ''
      PRINT *,  '  *** WARNING: 3D turbulent mixing specified but l3dturb=F ***'
      PRINT *,  '  (Your specification: itype_turbmix is not respected for '//TRIM(yshort_name)//')'
      PRINT *,  ''
    ENDIF
    ierr = 5000_iintegers
    CALL trcr_meta_set( itmp, i_trcr, T_TURB_ID, T_TURB_1D )
    IF (itmp /= 0_iintegers) RETURN
  ENDIF

  ! it is not possible to have turbulent mixing if ltur is set to 
  ! FALSE 
  ierr = 5000_iintegers
  CALL trcr_meta_get( itmp, i_trcr, T_TURB_ID, itmp1 )
  IF (itmp /= 0_iintegers) RETURN
  IF ( (itmp1 == T_TURB_1D .OR. itmp1 == T_TURB_3D) .AND. .NOT. ltur ) THEN
    IF (izdebug > 1_iintegers) THEN
      PRINT *,  ''
      PRINT *,  '  *** WARNING: turbulent mixing specified but ltur=F ***'
      PRINT *,  '  (Your specification: itype_turbmix is not respected for '//TRIM(yshort_name)//')'
      PRINT *,  ''
    ENDIF
    ierr = 5000_iintegers
    CALL trcr_meta_set( itmp, i_trcr, T_TURB_ID, T_TURB_OFF )
    IF (itmp /= 0_iintegers) RETURN
  ENDIF

  ! it is not possible to have 3D turbulent mixing in case of Leapfrog
  ierr = 5000_iintegers
  CALL trcr_meta_get( itmp, i_trcr, T_TURB_ID, itmp1 )
  IF (itmp /= 0_iintegers) RETURN
  IF ( itmp1 == T_TURB_3D .AND. .NOT. l2tls ) THEN
    IF (izdebug > 1_iintegers) THEN
      PRINT *,  ''
      PRINT *,  '  *** WARNING: 3D turbulent mixing chosen in combination with Leapfrog (not possible!) ***'
      PRINT *,  '  (Your specification: itype_turbmix is not respected for '//TRIM(yshort_name)//')'
      PRINT *,  ''
    ENDIF
    ierr = 5000_iintegers
    CALL trcr_meta_set( itmp, i_trcr, T_TURB_ID, T_TURB_1D )
    IF (itmp /= 0_iintegers) RETURN
  ENDIF

!------------------------------------------------------------------------------
!  Section 4: Indices determination
!------------------------------------------------------------------------------

  ! determine index of tracer field (in trcr_data) for current tracer
  trcr(i_trcr)%idx_trcr = i_trcr

  ! determine index of boundary field (in trcr_data_bd) for current tracer
  ierr = 5000_iintegers
  CALL trcr_meta_get( itmp, i_trcr, T_LBC_ID, itmp1 )
  IF (itmp /= 0_iintegers) RETURN
  IF ( ANY( itmp1 == (/ T_LBC_ZERO, T_LBC_FILE, T_LBC_CST, T_LBC_ZEROGRAD, T_LBC_USER /) ) ) THEN
    i_bd =  n_bd + 1_iintegers
    trcr(i_trcr)%idx_bd = i_bd
  ENDIF

  ! check for name uniqueness
  ierr = 5000_iintegers
  CALL trcr_meta_get( itmp, T_NAME_ID, ytmp )
  IF (itmp /= 0_iintegers) RETURN
  ierr = 4070_iintegers
  DO i = 1, n_trcr
    IF (TRIM(ytmp(i)) == TRIM(yshort_name)) RETURN
  ENDDO

  ! check validity range of integers (mandatory and optional)
  ierr = 5000_iintegers
  CALL trcr_meta_get( itmp, i_trcr, T_GRBPARAM_ID, itmp1 )
  IF (itmp /= 0_iintegers) RETURN
  ierr = 4080_iintegers
  IF ( itmp1 < 0_iintegers .OR. &     !GRIB param. num.
       itmp1 > 255_iintegers ) RETURN 
  ierr = 5000_iintegers
  CALL trcr_meta_get( itmp, i_trcr, T_ADV_ID, itmp1 )
  IF (itmp /= 0_iintegers) RETURN
  ierr = 4081_iintegers
  IF ( itmp1 < 0_iintegers .OR.  &     !advection type
       itmp1 > T_ADV_MAX ) RETURN
  ierr = 5000_iintegers
  CALL trcr_meta_get( itmp, i_trcr, T_DIFF_ID, itmp1 )
  IF (itmp /= 0_iintegers) RETURN
  ierr = 4082_iintegers
  IF ( itmp1 < 0_iintegers .OR. &     !horiz. diff. type
       itmp1 > T_DIFF_MAX ) RETURN
  ierr = 5000_iintegers
  CALL trcr_meta_get( itmp, i_trcr, T_TURB_ID, itmp1 )
  IF (itmp /= 0_iintegers) RETURN
  ierr = 4083_iintegers
  IF ( itmp1 < 0_iintegers .OR. &  !turb. mixing type
       itmp1 > T_TURB_MAX ) RETURN
  ierr = 5000_iintegers
  CALL trcr_meta_get( itmp, i_trcr, T_CONV_ID, itmp1 )
  IF (itmp /= 0_iintegers) RETURN
  ierr = 4084_iintegers
  IF ( itmp1 < 0_iintegers .OR. &     !convective passive transport
       itmp1 > T_CONV_MAX ) RETURN
  ierr = 5000_iintegers
  CALL trcr_meta_get( itmp, i_trcr, T_INI_ID, itmp1 )
  IF (itmp /= 0_iintegers) RETURN
  ierr = 4085_iintegers
  IF ( itmp1 < 0_iintegers .OR.  &     !IC type
       itmp1 > T_INI_MAX ) RETURN
  ierr = 5000_iintegers
  CALL trcr_meta_get( itmp, i_trcr, T_LBC_ID, itmp1 )
  IF (itmp /= 0_iintegers) RETURN
  ierr = 4086_iintegers
  IF ( itmp1 < 0_iintegers .OR.   &    !lat. BC type
       itmp1 > T_LBC_MAX ) RETURN
  ierr = 5000_iintegers
  CALL trcr_meta_get( itmp, i_trcr, T_BBC_ID, itmp1 )
  IF (itmp /= 0_iintegers) RETURN
  ierr = 4087_iintegers
  IF ( itmp1 < 0_iintegers .OR.   &    !bottom BC type
       itmp1 > T_BBC_MAX ) RETURN
  ierr = 5000_iintegers
  CALL trcr_meta_get( itmp, i_trcr, T_RELAX_ID, itmp1 )
  IF (itmp /= 0_iintegers) RETURN
  ierr = 4088_iintegers
  IF ( itmp1 < 0_iintegers .OR.   &    !relaxation type
       itmp1 > T_RELAX_MAX ) RETURN
  ierr = 5000_iintegers
  CALL trcr_meta_get( itmp, i_trcr, T_DAMP_ID, itmp1 )
  IF (itmp /= 0_iintegers) RETURN
  ierr = 4089_iintegers
  IF ( itmp1 < 0_iintegers .OR.   &    !damping type
       itmp1 > T_DAMP_MAX ) RETURN
  ierr = 5000_iintegers
  CALL trcr_meta_get( itmp, i_trcr, T_CLP_ID, itmp1 )
  IF (itmp /= 0_iintegers) RETURN
  ierr = 4090_iintegers
  IF ( itmp1 < 0_iintegers .OR. &     !clipping type
       itmp1 > T_CLP_MAX ) RETURN

!------------------------------------------------------------------------------
!  Section 5: Tracer registration and subroutine outputs
!------------------------------------------------------------------------------  

  ! all checks passed, we are safe to accept this tracer
  n_trcr = i_trcr
  n_bd   = i_bd

  ! return indices if requested
  IF (PRESENT(idx_trcr)) THEN
    idx_trcr = trcr(i_trcr)%idx_trcr
  END IF

  !******************************************************************

#else

  !******************************************************************
  ! MESSy Tracer Interface: define new tracer
  !******************************************************************

  IF (PRESENT(ystandard_name)) THEN
     standard_name = TRIM(ystandard_name)
  ELSE IF (PRESENT(ylong_name)) THEN
     standard_name = TRIM(ylong_name)
  ELSE
     standard_name = TRIM(yshort_name)
  ENDIF


  CALL new_tracer(status, GPTRSTR, TRIM(yshort_name), 'COSMO', &
       longname = TRIM(standard_name), &
       quantity = AMOUNTFRACTION,      &
       unit     = yunits,              &
       medium   = AIR,                 &
       idx      = idt                  )
  CALL tracer_halt(substr, status)
  IF (PRESENT(idx_trcr)) idx_trcr = idt

  CALL set_tracer(status, GPTRSTR, idt, I_GRIBTAB,   igribtable)
  CALL tracer_halt(substr, status)
  CALL set_tracer(status, GPTRSTR, idt, I_GRIBPARAM, igribparam)
  CALL tracer_halt(substr, status)
  IF (PRESENT(itype_adv)) THEN
     CALL set_tracer(status, GPTRSTR, idt, I_ADVECT,  itype_adv)
     CALL tracer_halt(substr, status)
  ENDIF
  IF (PRESENT(itype_turbmix)) THEN
     CALL set_tracer(status, GPTRSTR, idt, I_VDIFF,   itype_turbmix)
     CALL tracer_halt(substr, status)
  ENDIF
  IF (PRESENT(itype_diff)) THEN
     CALL set_tracer(status, GPTRSTR, idt, I_HDIFF,   itype_diff)
     CALL tracer_halt(substr, status)
  ENDIF
  IF (PRESENT(itype_passconv)) THEN
     CALL set_tracer(status, GPTRSTR, idt, I_CONVECT, itype_passconv)
     CALL tracer_halt(substr, status)
  ENDIF
  IF (PRESENT(itype_ini)) THEN
     IF (itype_ini == T_INI_FILE) THEN
        ! misuse MMD INIT SWITCH here for offline file initialisation
        CALL set_tracer(status, GPTRSTR, idt, I_MMD_INIT, ON)
        CALL tracer_halt(substr, status)
     ENDIF
     CALL set_tracer(status, GPTRSTR, idt, I_INITIAL, itype_ini)
     CALL tracer_halt(substr, status)
  ENDIF
  IF (PRESENT(itype_lbc)) THEN
     CALL set_tracer(status, GPTRSTR, idt, I_LBC, itype_lbc)
     CALL tracer_halt(substr, status)
  ENDIF
  IF (PRESENT(itype_bbc)) THEN
     CALL set_tracer(status, GPTRSTR, idt, I_RELAX, itype_relax)
     CALL tracer_halt(substr, status)
  ENDIF
  IF (PRESENT(itype_relax)) THEN
     CALL set_tracer(status, GPTRSTR, idt, I_RELAX, itype_relax)
     CALL tracer_halt(substr, status)
  ENDIF
  IF (PRESENT(itype_damp)) THEN
     CALL set_tracer(status, GPTRSTR, idt, I_DAMP,  itype_damp)
     CALL tracer_halt(substr, status)
  ENDIF
  IF (PRESENT(itype_clip)) THEN
     write (*,*) 'NO TRACE CLIPPING FOR MESSy TRACER ',TRIM(yshort_name)&
          , ' OMITTED !!'
  END IF
  ierr = status

  ! all checks passed, we are safe to accept this tracer
  CALL get_tracer_set_info(ierr, GPTRSTR,ntrac=n_trcr)
  CALL tracer_halt(substr,ierr)
  n_bd   = n_trcr
#endif

  ! no errors
  ierr = 0_iintegers

!------------------------------------------------------------------------------
! End of module procedure trcr_new
!------------------------------------------------------------------------------

END SUBROUTINE trcr_new

!==============================================================================


!==============================================================================
!+ Module procedure "trcr_alloc" in "src_tracer" to allocate memory.
!------------------------------------------------------------------------------

SUBROUTINE trcr_alloc ( ierr )
!------------------------------------------------------------------------------
!
! Description:
!
! This subroutine allocates the memory for the tracer module. After memory
! allocation, the pointers of the tracer structure need to be updated to
! point to the correct location. Also, the corresponding pointers in the
! variable table (src_setup_vartab.f90) needs to be updated.
!
! Method:
!
! Dynamic memory allocation using the Fortran ALLOCATE with dimensions
! specified by (ie,je,ke,n_trcr,nt).
! Allocates tracer structures on accelerator using OpenACC
!------------------------------------------------------------------------------

 IMPLICIT NONE

! Subroutine arguments:
! --------------------

 INTEGER(KIND=iintegers), INTENT(OUT) :: ierr    ! error status


! Local variables
 INTEGER(KIND=iintegers) ::  &
  izl   ,  & ! for local error-code
  i          ! index for looping over tracers

!------------------------------------------------------------------------------
! Begin Subroutine trcr_alloc
!------------------------------------------------------------------------------
  ! initialize error status
  ierr = 0_iintegers

  ! check workflow
  ierr = 4100_iintegers
  IF ( i_status /= T_STAT_DEFINE ) RETURN


  ! define number of time levels
  IF ( l2tls ) THEN
#ifndef MESSY
    n_tlev = 2_iintegers
  ELSE
    n_tlev = 3_iintegers
#else
    n_tlev = 3_iintegers
  ELSE
    n_tlev = 6_iintegers
#endif
  ENDIF

  ! allocate memory
#ifdef IJNK_STORAGE
#ifdef MESSY
  ! provide tracer meta data and memory earlier
  CALL messy_tracer_meta
  n_trcr = ntrac_gp
  
  DO i=1,n_trcr
     trcr(i)%idx_trcr =i
     trcr(i)%idx_bd   =i
  ENDDO
  trcr_data      => xt_tl
  trcr_data_bd   => xt_bd
  trcr_data_tens => xtte
  ierr = 0_iintegers
#else
  ierr = 0_iintegers
  ALLOCATE ( trcr_data(ie, je, n_trcr, ke, n_tlev), STAT=izl )
  ierr = ierr + izl
  ALLOCATE ( trcr_data_bd(ie, je, n_bd, ke, 2),     STAT=izl )
  ierr = ierr + izl
  ALLOCATE ( trcr_data_tens(ie, je, n_trcr, ke),    STAT=izl )
  ierr = ierr + izl
#endif
#else
  ierr = 0_iintegers
  ALLOCATE ( trcr_data(ie, je, ke, n_trcr, n_tlev), STAT=izl ) 
  ierr = ierr + izl
  ALLOCATE ( trcr_data_bd(ie, je, ke, n_bd, 2),     STAT=izl ) 
  ierr = ierr + izl
  ALLOCATE ( trcr_data_tens(ie, je, ke, n_trcr),    STAT=izl ) 
  ierr = ierr + izl
#endif

  ! debug message
  IF (izdebug > 10_iintegers) THEN
    PRINT *, '   ALLOCATED tracer fields:  ', ierr
  ENDIF

  ! check for error
  IF (ierr /= 0_iintegers) THEN
    ierr = 4105_iintegers
    RETURN
  ENDIF

  ! initialize arrays
  trcr_data      = 0.0_wp
  trcr_data_bd   = 0.0_wp
  trcr_data_tens = 0.0_wp

  ! update tracer module status
  i_status = T_STAT_ALLOC

  ! no errors
  ierr = 0_iintegers

  ! allocate tracer structures on accelerator
  !$acc enter data copyin(trcr_data,trcr_data_bd,trcr_data_tens)

!------------------------------------------------------------------------------
! End of module procedure trcr_alloc
!------------------------------------------------------------------------------

END SUBROUTINE trcr_alloc

!==============================================================================
!==============================================================================
!+ Module procedure "trcr_setup_vartab" in "src_tracer" to add tracers to the
!                   variable var.
!------------------------------------------------------------------------------

SUBROUTINE trcr_setup_vartab ( ierr )
!------------------------------------------------------------------------------
!
! Description:
!
! This routine enters the tracer variables in the variable table as it is done 
! for other model variables in src_setup_vartab.f90.
!
! Method:
!
! Identical as src_setup_vartab.f90 but with additional checks to ensure
! unicity.
!------------------------------------------------------------------------------

 IMPLICIT NONE

! Subroutine arguments:
! --------------------

 INTEGER(KIND=iintegers), INTENT(OUT) :: ierr    ! error status

! Local variables

 INTEGER(KIND=iintegers) :: i_trcr      ! index to loop inside tracer struct.
 INTEGER(KIND=iintegers) :: itmp, i1,i2,i3    ! variable table inidices
 INTEGER(KIND=iintegers) :: izbbc
 CHARACTER(LEN=ilen_sn)  :: yname
 CHARACTER(LEN=ilen_un)  :: yunits
 CHARACTER(LEN=ilen_ln)  :: ystandard_name, ylong_name

 ! The following numbers according to grib table definition in
 ! organize_data
 INTEGER(KIND=iintegers), PARAMETER :: messy_gribtab_start = 7  ! 230
 INTEGER(KIND=iintegers), PARAMETER :: messy_gribtab_max   = 15 ! 238

!------------------------------------------------------------------------------
! Begin Subroutine trcr_setup_vartab
!------------------------------------------------------------------------------
  
  ! initialize error status
  ierr = 0_iintegers

  ! check workflow
  ierr = 4150_iintegers
  IF ( i_status /= T_STAT_ALLOC ) RETURN

  ! loop over all registered tracers
  DO i_trcr = 1, n_trcr

    ! Check if GRIB table number allowed
    ierr = 5000_iintegers
    CALL trcr_meta_get( itmp, i_trcr, T_GRBTABLE_ID, i3 )
    IF (itmp /= 0_iintegers) RETURN
    ierr = 4155_iintegers
    IF ( (i3 < 1_iintegers &
#ifdef MESSY
         .AND. i3 /= -1_iintegers & ! DEFAULT FOR NOT SET GRIBTAB numbers
#endif
         ) .OR. i3 > num_gribtabs ) RETURN

    ! check if name already exists in variable table
    ierr = 5000_iintegers
    CALL trcr_meta_get( itmp, i_trcr, T_NAME_ID, yname )
    IF (itmp /= 0_iintegers) RETURN

    ierr = 4157_iintegers
    DO i3 = 1, num_gribtabs
      DO i2 = 0, 255
        DO i1 = 1, max_gribrep
          IF ( TRIM(var(i1,i2,i3)%name) == TRIM(yname) ) RETURN
        ENDDO
      ENDDO
    ENDDO

    ! get grib table and parameter number
    ierr = 5000_iintegers
    CALL trcr_meta_get( itmp, i_trcr, T_GRBPARAM_ID, i2 )
    IF (itmp /= 0_iintegers) RETURN
    CALL trcr_meta_get( itmp, i_trcr, T_GRBTABLE_ID, i3 )
    IF (itmp /= 0_iintegers) RETURN

#ifdef MESSY
    IF (i2 == -1 .AND. i3 == -1) THEN
       tab_loop: DO i3 = messy_gribtab_start, messy_gribtab_max
          DO i2 = 1, 254
              i1 = 1
              IF (LEN_TRIM(var(i1,i2,i3)%name) == 0) THEN
                 CALL set_tracer(ierr, GPTRSTR, i_trcr, I_GRIBTAB,   i3)
                 CALL tracer_halt('trcr_setup_vartab', ierr)
                 CALL set_tracer(ierr, GPTRSTR, i_trcr, I_GRIBPARAM, i2)
                 CALL tracer_halt('trcr_setup_vartab', ierr)

                 EXIT tab_loop
              END IF
          ENDDO
       ENDDO tab_loop
       ! to many tracers requested, end of dedicated grib-tables
       ierr = 4159_iintegers
       IF ( i3 > messy_gribtab_max) RETURN
    ELSE
#endif
    ! find empty position in variable table
    DO i1 = 1, max_gribrep
      IF ( LEN_TRIM(var(i1,i2,i3)%name) == 0 ) EXIT
    ENDDO
    ierr = 4159_iintegers
    IF ( i1 > max_gribrep ) RETURN
#ifdef MESSY
    ENDIF
#endif

    ! get other metadata
    ierr = 5000_iintegers
    CALL trcr_meta_get( itmp, i_trcr, T_UNITS_ID, yunits )
    IF (itmp /= 0_iintegers) RETURN
    CALL trcr_meta_get( itmp, i_trcr, T_NCSTDNAME_ID, ystandard_name )
    IF (itmp /= 0_iintegers) RETURN
    CALL trcr_meta_get( itmp, i_trcr, T_NCLONGNAME_ID, ylong_name )
    IF (itmp /= 0_iintegers) RETURN

    ! setup variable table entry
    ! fields: name levtyp levtop levbot factor bias ntri rank p4 p4_bd p3 p3_bd p2
    !         idef_stat units standard_name long_name land/sea
    var(i1,i2,i3) = ar_des(                                                                &
      yname, 110, 0, 0, 1.0_wp, 0.0_wp, 0, 4, ke,                                          &
#ifdef IJNK_STORAGE
      trcr_data(:,:,trcr(i_trcr)%idx_trcr,:,:), trcr_data_bd(:,:,trcr(i_trcr)%idx_bd,:,:), &
#else
      trcr_data(:,:,:,trcr(i_trcr)%idx_trcr,:), trcr_data_bd(:,:,:,trcr(i_trcr)%idx_bd,:), &
#endif
      NULL(), NULL(), NULL(), 1,                                                           &
      yunits, ystandard_name, ylong_name, ' '                                              &
     )

  ENDDO

  ! no errors
  ierr = 0_iintegers

!------------------------------------------------------------------------------
! End of module procedure trcr_setup_vartab
!------------------------------------------------------------------------------

END SUBROUTINE trcr_setup_vartab

!==============================================================================

!==============================================================================
!+ Module procedure "trcr_print" in "src_tracer" to print the list of tracers.
!------------------------------------------------------------------------------

SUBROUTINE trcr_print ( ierr ) 
!------------------------------------------------------------------------------
!
! Description:
!
! Print the list of the tracers currently handled by the tracer module into
! the standard output.
!
! Method:
!
!
!------------------------------------------------------------------------------

 IMPLICIT NONE

! Subroutine arguments:
! --------------------

 INTEGER(KIND=iintegers), INTENT(OUT) :: ierr    ! error status


! Local parameters:
! ----------------

 INTEGER(KIND=iintegers)   :: i

!----------- End of header ----------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine trcr_print
!------------------------------------------------------------------------------
  ! initialize error status
  ierr = 0_iintegers

  ! check workflow
  ierr = 4200_iintegers
  IF ( i_status /= T_STAT_ALLOC ) RETURN

#ifndef MESSY
  ! write tracer information on PE 0 in the standard output
  IF (my_cart_id == 0_iintegers) THEN
    IF (n_trcr > 0_iintegers) THEN

      WRITE (*, '(A2)')  '  '
      WRITE (*, '(A50)')                                                      &
        'The tracers currently handled in the model are:'
      WRITE (*, '(A50)')                                                      &
        ' ============================================'
      WRITE (*, '(A2)')  '  '
      
      CALL metadata_print(mtrcr)

      WRITE (*, '(A2)')  '  '

    ELSE

      WRITE (*, '(A2)')  '  '
      WRITE (*, '(A70)')                                                      &
        '    Currently there are no tracers handled in the model.'
      WRITE (*, '(A2)')  '  '
    ENDIF

  ENDIF
#else
  IF (my_cart_id == 0_iintegers) THEN
     CALL print_tracer_set
  ENDIF
#endif

  ! no errors
  ierr = 0_iintegers

!------------------------------------------------------------------------------
! End of module procedure trcr_print
!------------------------------------------------------------------------------

END SUBROUTINE trcr_print

!==============================================================================

!==============================================================================
!+ Subroutine "trcr_get_index" in "src_tracer" to get back a tracer index.
!------------------------------------------------------------------------------

SUBROUTINE trcr_get_index ( ierr, yname, idx )
!------------------------------------------------------------------------------
!
! Description:
!
! This routine gets back a tracer index from a given tracer name
!
!
! Method:
!
! Identification by name in the data structure.
!------------------------------------------------------------------------------

 IMPLICIT NONE

! Subroutine arguments:
! --------------------

 INTEGER(KIND=iintegers),   INTENT(OUT) :: ierr  ! error status
 CHARACTER(LEN=*)       ,   INTENT(IN)  :: yname ! name of tracer
 INTEGER(KIND=iintegers),   INTENT(OUT) :: idx   ! index of tracer

! Local parameters:
! ----------------

 INTEGER(KIND=iintegers) :: i                    ! looping index
 INTEGER(KIND=iintegers) :: itmp
 CHARACTER(LEN=ilen_sn)  :: ytmp(n_trcr)         ! storage for tracer names

!----------- End of header ----------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine trcr_get_index
!------------------------------------------------------------------------------

  ! initialize error status
  ierr     = 0_iintegers

  ! check workflow
  ierr = 4350_iintegers
  IF ( i_status /= T_STAT_DEFINE .AND. i_status /= T_STAT_ALLOC ) RETURN

  ! find matching tracer (by name)
  ierr = 5000_iintegers
  CALL trcr_meta_get( itmp, T_NAME_ID, ytmp )
  IF (itmp /= 0_iintegers) RETURN

  ! search for tracer
  DO i = 1, n_trcr
    IF ( ytmp(i) == yname ) EXIT
  ENDDO

  ! get index or return error if no tracer has been found
  izerrind = i

  IF ( i <= n_trcr ) then
    idx = i
    ierr = 0_iintegers
  ELSE
    idx  = T_MISSING
    ierr = T_ERR_NOTFOUND
  ENDIF

!------------------------------------------------------------------------------
! End of module procedure trcr_get_index
!------------------------------------------------------------------------------
  
END SUBROUTINE trcr_get_index

!==============================================================================

!==============================================================================
!+ Subroutine "trcr_get_byindex" in "src_tracer" to get back a tracer.
!------------------------------------------------------------------------------

SUBROUTINE trcr_get_byindex ( ierr, idx_trcr, ptr, ptr_bd, ptr_bd_notlev, ptr_tens, ptr_tlev )
!------------------------------------------------------------------------------
!
! Description:
!
! This routine gets back a given tracer from the structure trcr.
!
!
! Method:
!
! Identification by index in the data structure.
!------------------------------------------------------------------------------

 IMPLICIT NONE

! Subroutine arguments:
! --------------------

 INTEGER(KIND=iintegers), INTENT(OUT)          :: ierr            ! error status
 INTEGER(KIND=iintegers), INTENT(IN)           :: idx_trcr        ! index for tracer data
 REAL(KIND=wp)          , POINTER,    OPTIONAL :: ptr(:,:,:)      ! pointer to tracer
 REAL(KIND=wp)          , POINTER,    OPTIONAL :: ptr_bd(:,:,:,:) ! pointer to tracer boundary 
 REAL(KIND=wp)          , POINTER,    OPTIONAL :: ptr_bd_notlev(:,:,:) ! pointer to tracer boundary without time level
 REAL(KIND=wp)          , POINTER,    OPTIONAL :: ptr_tens(:,:,:) ! pointer to tracer tendency 
 INTEGER(KIND=iintegers), INTENT(IN), OPTIONAL :: ptr_tlev        ! time level for pointers

! Local parameters:
! ----------------
 LOGICAL :: lfound             ! flag if tracer has been found

 INTEGER(KIND=iintegers) :: &
   i,                     &    ! looping index
   iztlev                      ! time level

!----------- End of header ----------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine trcr_get_byindex
!------------------------------------------------------------------------------

  ! initialize error status 
  ierr     = 0_iintegers
  izerrind = idx_trcr

  ! initialize pointer (for safetey)
  IF (PRESENT(ptr)) NULLIFY(ptr)
  IF (PRESENT(ptr_bd)) NULLIFY(ptr_bd)
  IF (PRESENT(ptr_bd_notlev)) NULLIFY(ptr_bd_notlev)
  IF (PRESENT(ptr_tens)) NULLIFY(ptr_tens)
  
  ! check workflow
  ierr = 4250_iintegers
  IF ( i_status /= T_STAT_DEFINE .AND. i_status /= T_STAT_ALLOC ) RETURN
  IF ( ( PRESENT(ptr) .OR. PRESENT(ptr_bd) .OR. PRESENT(ptr_tens))            &
        .AND. i_status /= T_STAT_ALLOC ) RETURN

  ! init
  lfound = .FALSE.

  ! get time level (only needed if ptr requested)
  IF (PRESENT(ptr_tlev)) THEN
    ! time level is only relevant for ptr  and not _bd or _tens
    ierr = 4253_iintegers
    IF (.NOT. ( PRESENT(ptr) .OR. PRESENT(ptr_bd_notlev) ) ) RETURN
    ! set requested timelevel
#ifdef STATIC_FIELDS
    iztlev = -1_iintegers
    IF (ptr_tlev == nnew) THEN
      iztlev = t_nnew
    ELSEIF (ptr_tlev == nnow) THEN
      iztlev = t_nnow
    ELSEIF (.NOT. l2tls .AND. ptr_tlev == nold) THEN
      iztlev = t_nold
    ENDIF
#else
    iztlev = ptr_tlev
#endif
  ELSE
    ! set default time level
#ifdef STATIC_FIELDS
#ifndef MESSY
    iztlev = 1
#else
    iztlev = t_nnow
#endif
#else
    iztlev = nnow
#endif
  ENDIF

  ! check requested timelevel
  ierr = 4255_iintegers
  IF (iztlev < 0_iintegers .OR. iztlev > n_tlev) RETURN

#ifndef MESSY
  ! find matching tracer (by index)
  DO i = 1, n_trcr
     IF ( trcr(i)%idx_trcr == idx_trcr ) THEN
#else
        IF (idx_trcr < 0 .OR. idx_trcr > n_trcr) THEN
           ierr = T_ERR_NOTFOUND
           RETURN
        END IF
        i=idx_trcr
#endif
        
        ! we found a matching tracer
        lfound = .TRUE.
        
        ! provide pointer (if requested)
        IF ( PRESENT(ptr) ) THEN
          ierr = 4262_iintegers
          IF (trcr(i)%idx_trcr == T_MISSING) RETURN
#ifdef IJNK_STORAGE
          ptr => trcr_data(:,:,trcr(i)%idx_trcr,:,iztlev)
#else
          ptr => trcr_data(:,:,:,trcr(i)%idx_trcr,iztlev)
#endif
        ENDIF
      
        ! provide boundary pointer (if requested)
        IF ( PRESENT(ptr_bd) ) THEN
          ierr = 4265_iintegers
          IF (trcr(i)%idx_bd == T_MISSING) RETURN
#ifdef IJNK_STORAGE
          ptr_bd => trcr_data_bd(:,:,trcr(i)%idx_bd,:,:)
#else
          ptr_bd => trcr_data_bd(:,:,:,trcr(i)%idx_bd,:)
#endif
        ENDIF

        ! provide boundary pointer without time level (if requested)
        IF ( PRESENT(ptr_bd_notlev) ) THEN
          ierr = 4265_iintegers
          IF (trcr(i)%idx_bd == T_MISSING) RETURN
#ifdef IJNK_STORAGE
          ptr_bd_notlev => trcr_data_bd(:,:,trcr(i)%idx_bd,:,iztlev)
#else
          ptr_bd_notlev => trcr_data_bd(:,:,:,trcr(i)%idx_bd,iztlev)
#endif
        ENDIF

        ! provide tendency pointer (if requested)
        IF (PRESENT(ptr_tens)) THEN
          ierr = 4270_iintegers
          IF (trcr(i)%idx_trcr == T_MISSING) RETURN
#ifdef IJNK_STORAGE
          ptr_tens => trcr_data_tens(:,:,trcr(i)%idx_trcr,:)
#else
          ptr_tens => trcr_data_tens(:,:,:,trcr(i)%idx_trcr)
#endif
        ENDIF

        
        ! tracer indexes are unique, so we can safely exit
#ifndef MESSY
        EXIT
     ENDIF
  ENDDO
#endif
  ! return error if no tracer has been found
  ierr = T_ERR_NOTFOUND
  IF (.NOT. lfound) RETURN
  
  ! no errors
  ierr = 0_iintegers
  
!------------------------------------------------------------------------------
! End of module procedure trcr_get_byindex
!------------------------------------------------------------------------------
  
END SUBROUTINE trcr_get_byindex

!==============================================================================

!==============================================================================
!+ Subroutine "trcr_get_byname" in "src_tracer" to get back a tracer.
!------------------------------------------------------------------------------

SUBROUTINE trcr_get_byname ( ierr, yname, ptr, ptr_bd, ptr_tens, ptr_tlev)
!------------------------------------------------------------------------------
!
! Description:
!
! This routine gets back a given tracer from the structure trcr.
!
! Method:
!
!Identification by index in the data structure.
!------------------------------------------------------------------------------

 IMPLICIT NONE

! Subroutine arguments:
! --------------------

 INTEGER(KIND=iintegers), INTENT(OUT)           :: ierr            ! error status
 CHARACTER(LEN=*)       , INTENT(IN)            :: yname           ! name of tracer
 REAL(KIND=wp)          , POINTER    , OPTIONAL :: ptr(:,:,:)      ! pointer to tracer
 REAL(KIND=wp)          , POINTER    , OPTIONAL :: ptr_bd(:,:,:,:) ! pointer to tracer boundary 
 REAL(KIND=wp)          , POINTER    , OPTIONAL :: ptr_tens(:,:,:) ! pointer to tracer tendency 
 INTEGER(KIND=iintegers), INTENT(IN) , OPTIONAL :: ptr_tlev        ! time level for pointer


! Local parameters:
! ----------------
 INTEGER(KIND=iintegers) :: &
   i_trcr                         ! tracer index
 CHARACTER(LEN=*), PARAMETER :: ysubstr = 'trcr_get_byname'

!----------- End of header ----------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine trcr_get_byname
!------------------------------------------------------------------------------

  ! initialize error status
  ierr = 0_iintegers

  ! initialize pointer (for safetey)
  IF (PRESENT(ptr)) NULLIFY(ptr)
  IF (PRESENT(ptr_bd)) NULLIFY(ptr_bd)
  IF (PRESENT(ptr_tens)) NULLIFY(ptr_tens)

  ! get tracer index
  CALL trcr_get_index( ierr, yname, i_trcr )
  IF ( ierr /= 0_iintegers ) RETURN

  ! if found, call trcr_get (by index)
  CALL trcr_get_byindex( ierr, i_trcr, ptr=ptr, ptr_bd=ptr_bd, ptr_tens= ptr_tens, ptr_tlev=ptr_tlev)
  IF ( ierr /= 0_iintegers ) RETURN

END SUBROUTINE trcr_get_byname

!==============================================================================
!==============================================================================
!+ Subroutine "trcr_get_block" in "src_tracers" to get back a block of tracers.
!------------------------------------------------------------------------------

SUBROUTINE trcr_get_block( ierr, idx_start, idx_end, ptr, ptr_bd, ptr_tens, ptr_wtlev, ptr_tlev )
!------------------------------------------------------------------------------
!
! Description:
!
! This routine gets back a given tracer from the structure trcr.
!
!
! Method:
!
! Identification by index in the data structure.
!------------------------------------------------------------------------------

 IMPLICIT NONE

! Subroutine arguments:
! --------------------

 INTEGER(KIND=iintegers), INTENT(OUT)           :: ierr                 ! error status
 INTEGER(KIND=iintegers), INTENT(IN)            :: idx_start            ! start index for tracer data block
 INTEGER(KIND=iintegers), INTENT(IN)            :: idx_end              ! end index for tracer data block
 REAL(KIND=wp)          , POINTER    , OPTIONAL :: ptr(:,:,:,:)         ! pointer to tracer
 REAL(KIND=wp)          , POINTER    , OPTIONAL :: ptr_bd(:,:,:,:,:)    ! pointer to tracer boundary 
 REAL(KIND=wp)          , POINTER    , OPTIONAL :: ptr_tens(:,:,:,:)    ! pointer to tracer tendency 
 REAL(KIND=wp)          , POINTER    , OPTIONAL :: ptr_wtlev(:,:,:,:,:) ! pointer to tracer with time levels
 INTEGER(KIND=iintegers), INTENT(IN) , OPTIONAL :: ptr_tlev             ! time level for pointers

! Local parameters:
! ----------------
 INTEGER(KIND=iintegers) :: &
   i,i_start,i_end,       &    ! looping index
   iztlev                      ! time level

!----------- End of header ----------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine trcr_get_block
!------------------------------------------------------------------------------
#ifndef MESSY

#ifdef IJNK_STORAGE
  WRITE(0,*) 'ERROR: NON IJK-STORAGE FORMAT IS NOT SUPPORTED BY TRCR_GET_BLOCK'
  STOP
#else

  ! initialize error status 
  ierr     = 0_iintegers
  izerrind = idx_start

  ! initialize pointer (for safetey)
  IF (PRESENT(ptr)) NULLIFY(ptr)
  IF (PRESENT(ptr_bd)) NULLIFY(ptr_bd)
  IF (PRESENT(ptr_tens)) NULLIFY(ptr_tens)
  IF (PRESENT(ptr_wtlev)) NULLIFY(ptr_wtlev)

  ! check workflow
  ierr = 4250_iintegers
  IF ( i_status /= T_STAT_DEFINE .AND. i_status /= T_STAT_ALLOC ) RETURN
  IF ( ( PRESENT(ptr) .OR. PRESENT(ptr_bd) .OR. PRESENT(ptr_tens))            &
        .AND. i_status /= T_STAT_ALLOC ) RETURN

  ! get time level (only needed if ptr requested)
  IF (PRESENT(ptr_tlev)) THEN
    ! time level is only relevant for ptr and not _bd or _tens
    ierr = 4253_iintegers
    IF (.NOT. PRESENT(ptr)) RETURN
    ! set requested timelevel
    iztlev = ptr_tlev
  ELSE
    ! set default time level
    iztlev = nnow
  ENDIF
  ierr = 4255_iintegers
  IF (iztlev < 0_iintegers .OR. iztlev > n_tlev) RETURN

  ! find matching start tracer (by index)
  i_start = T_MISSING
  DO i = 1, n_trcr
     IF ( trcr(i)%idx_trcr == idx_start ) THEN
        i_start = i
     ENDIF
     IF ( trcr(i)%idx_trcr == idx_end ) THEN
        i_end = i
     ENDIF
  ENDDO

  ! return error if no tracer has been found
  ierr = T_ERR_NOTFOUND
  IF (i_start == T_MISSING .OR. i_end == T_MISSING) RETURN

  ! check that we have memory allocated and that it is contiguous
  ! if ok, then provide pointer
  IF ( PRESENT(ptr) ) THEN
    DO i = i_start, i_end
      ierr = 4262_iintegers
      IF (trcr(i)%idx_trcr == T_MISSING) RETURN
      ierr = 4282_iintegers
      IF (i > i_start) THEN
        IF ( PRESENT(ptr) .AND. (trcr(i)%idx_trcr /= trcr(i-1)%idx_trcr + 1) ) RETURN
      END IF
    END DO
    ptr => trcr_data(:,:,:,trcr(i_start)%idx_trcr:trcr(i_end)%idx_trcr,iztlev)
  END IF
  IF ( PRESENT(ptr_wtlev) ) THEN
    DO i = i_start, i_end
      ierr = 4262_iintegers
      IF (trcr(i)%idx_trcr == T_MISSING) RETURN
      ierr = 4282_iintegers
      IF (i > i_start) THEN
        IF ( PRESENT(ptr_wtlev) .AND. (trcr(i)%idx_trcr /= trcr(i-1)%idx_trcr + 1) ) RETURN
      END IF
    END DO
    ptr_wtlev => trcr_data(:,:,:,trcr(i_start)%idx_trcr:trcr(i_end)%idx_trcr,:)
  END IF
  IF ( PRESENT(ptr_bd) ) THEN
    DO i = i_start, i_end
      ierr = 4265_iintegers
      IF (trcr(i)%idx_bd == T_MISSING) RETURN
      ierr = 4285_iintegers
      IF (i > i_start) THEN
        IF ( PRESENT(ptr_bd) .AND. (trcr(i)%idx_bd /= trcr(i-1)%idx_bd + 1) ) RETURN
      END IF
    END DO
    ptr_bd => trcr_data_bd(:,:,:,trcr(i_start)%idx_bd:trcr(i_end)%idx_bd,:)
  END IF
  IF (PRESENT(ptr_tens)) THEN
    DO i = i_start, i_end
      ierr = 4270_iintegers
      IF (trcr(i)%idx_trcr == T_MISSING) RETURN
      ierr = 4290_iintegers
      IF (i > i_start) THEN
        IF ( PRESENT(ptr_tens) .AND. (trcr(i)%idx_trcr /= trcr(i-1)%idx_trcr + 1) ) RETURN
      END IF
    END DO
    ptr_tens => trcr_data_tens(:,:,:,trcr(i_start)%idx_trcr:trcr(i_end)%idx_trcr)
  END IF

  ! no errors
  ierr = 0_iintegers
#endif

#else
  CALL error_bi('subroutine trcr_get_block not implemented for MESSy', ' ')
#endif

!------------------------------------------------------------------------------
! End of module procedure trcr_get_block
!------------------------------------------------------------------------------
  
END SUBROUTINE trcr_get_block

!==============================================================================
!==============================================================================
!+ Subroutine "trcr_check_index" in "src_tracer" to check a tracer index.
!------------------------------------------------------------------------------

SUBROUTINE trcr_check_index ( ierr, idx )
!------------------------------------------------------------------------------
!
! Description:
!
! This routine checks the validity of a tracer index
!
!
! Method:
!
!------------------------------------------------------------------------------

 IMPLICIT NONE

! Subroutine arguments:
! --------------------

 INTEGER(KIND=iintegers),   INTENT(OUT)   :: ierr  ! error status
 INTEGER(KIND=iintegers),   INTENT(IN)    :: idx   ! index of tracer

! Local parameters:
! ----------------


!----------- End of header ----------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine trcr_check_index
!------------------------------------------------------------------------------

  ! initialize error status
  ierr     = 0_iintegers
  izerrind = idx

  ! check workflow
  ierr = 4350_iintegers
  IF ( i_status /= T_STAT_DEFINE .AND. i_status /= T_STAT_ALLOC ) RETURN

#ifndef MESSY
  ! check index and return error if no tracer has been found
  IF ( i_status == T_STAT_DEFINE ) THEN
    ! if we are in the middle of defining tracers, this subroutine might also
    ! be called when defining a new tracer but before n_trcr has been updated,
    ! so we relax the check (TODO: this might go awfully wrong)
    IF ((idx >= 1) .AND. (idx <= n_trcr+1)) THEN
      ierr = 0_iintegers
    ELSE
      ierr = T_ERR_NOTFOUND
    ENDIF
  ELSE
#endif
    IF ((idx >= 1) .AND. (idx <= n_trcr)) THEN
      ierr = 0_iintegers
    ELSE
      ierr = T_ERR_NOTFOUND
    ENDIF
#ifndef MESSY
  ENDIF
#endif
!------------------------------------------------------------------------------
! End of module procedure trcr_check_index
!------------------------------------------------------------------------------
  
END SUBROUTINE trcr_check_index

!==============================================================================
!==============================================================================
!+ Module procedure "trcr_swap" in "src_tracer" to allocate memory.
!------------------------------------------------------------------------------

SUBROUTINE trcr_swap ( ierr )
!------------------------------------------------------------------------------
!
! Description:
!
! This subroutine swaps the timelevels for the tracers in case of static
! memory allocation for individual timelevels
!
! Method:
!
! Copy the fields from nnew -> nnow (and nnow -> nold)
!------------------------------------------------------------------------------

 IMPLICIT NONE

! Subroutine arguments:
! --------------------

 INTEGER(KIND=iintegers), INTENT(OUT) :: ierr    ! error status

!------------------------------------------------------------------------------
! Begin Subroutine trcr_swap
!------------------------------------------------------------------------------

#ifdef STATIC_FIELDS

  ! initialize error status
  ierr = 0_iintegers

  ! check workflow
  ierr = 4290_iintegers
  IF ( i_status /= T_STAT_ALLOC ) RETURN

  ! swap timelevels
  IF ( .NOT. l2tls ) THEN
    trcr_data(:,:,:,:,t_nold) = trcr_data(:,:,:,:,t_nnow)
  ENDIF
  trcr_data(:,:,:,:,t_nnow) = trcr_data(:,:,:,:,t_nnew)

  ! debug message
  IF (izdebug > 10_iintegers) THEN
    PRINT *, '      Timlevel swapping for tracers'
  ENDIF

  ! no errors
  ierr = 0_iintegers

#endif

!------------------------------------------------------------------------------
! End of module procedure trcr_swap
!------------------------------------------------------------------------------

END SUBROUTINE trcr_swap

!==============================================================================
!==============================================================================
!+ Module procedure "trcr_cleanup" in "src_tracer" to clean memory.
!------------------------------------------------------------------------------

SUBROUTINE trcr_cleanup ( ierr )
!------------------------------------------------------------------------------
!
! Description:
!
! This subroutine de-allocates the memory for the tracer module.
!
!
! Method:
!
! Dynamic memory de-allocation using the Fortran DEALLOCATE
!------------------------------------------------------------------------------

 IMPLICIT NONE

! Subroutine arguments:
! --------------------

 INTEGER(KIND=iintegers)  , INTENT(OUT)           :: ierr          ! error status

#ifndef MESSY
! Local parameters:
! ----------------

 INTEGER(KIND=iintegers)   ::     &
  i           ,  &                        ! looping index
  izl                                     ! for local error-code
!------------------------------------------------------------------------------
! Begin Subroutine trcr_cleanup
!------------------------------------------------------------------------------

  ! initialize error status
  ierr = 0_iintegers

  ! check workflow
  ierr = 4300_iintegers
  IF ( i_status /= T_STAT_ALLOC ) RETURN

  ! deallocate memory
  ierr = 0_iintegers
  DEALLOCATE (trcr_data       , STAT=izl) ;  ierr = ierr + izl
  DEALLOCATE (trcr_data_bd    , STAT=izl) ;  ierr = ierr + izl
  DEALLOCATE (trcr_data_tens  , STAT=izl) ;  ierr = ierr + izl

  ! init tracer meta-data storage
  CALL metadata_finish(mtrcr)

  ! debug message
  IF (izdebug > 10_iintegers) THEN
    PRINT *, '   DEALLOCATED tracer fields:  ', ierr
  ENDIF

  ! check for error
  IF (ierr /= 0_iintegers) THEN
    ierr = 4305_iintegers
    RETURN
  ENDIF

  ! update tracer module status
  i_status = T_STAT_FINISH

  ! no errors
  ierr = 0_iintegers
#else
  ierr = 0_iintegers
  CALL info_bi('subroutine trcr_cleanup not required for MESSy', '' )
#endif

!------------------------------------------------------------------------------
! End of module procedure trcr_cleanup
!------------------------------------------------------------------------------

END SUBROUTINE trcr_cleanup
!==============================================================================
!==============================================================================
!+ Module procedure "trcr_cleanup" in "src_tracer" to clean memory.
!------------------------------------------------------------------------------

SUBROUTINE trcr_calc ( flag, itrcr)
!------------------------------------------------------------------------------
!
! Description:
!
! This subroutine de-allocates the memory for the tracer module.
!
!
! Method:
!
! Dynamic memory de-allocation using the Fortran DEALLOCATE
!------------------------------------------------------------------------------

 IMPLICIT NONE

! Subroutine arguments:
! --------------------

 INTEGER(KIND=iintegers),           INTENT(IN)   :: flag  !  direction
 INTEGER(KIND=iintegers), OPTIONAL, INTENT(IN)   :: itrcr ! tracer_index

#ifdef MESSY
 INTEGER(KIND=iintegers) :: i

 SELECT CASE(flag)

 CASE(-1)

    IF (PRESENT(itrcr)) THEN
       xtte(:,:,itrcr,:) = xtte(:,:,itrcr,:)  &
            + ( xt (:,:,itrcr,:) - (xtm1(:,:,itrcr,:)  &
            + dt * xtte(:,:,itrcr,:))) / dt
       xt(:,:,itrcr,:)   = xtm1(:,:,itrcr,:)
    ELSE
       xtte(:,:,:,:) = xtte(:,:,:,:)  &
            + ( xt (:,:,:,:) - (xtm1(:,:,:,:) + dt * xtte(:,:,:,:))) / dt
       xt(:,:,:,:)   = xtm1(:,:,:,:)
    END IF

 CASE(1)
    ! CALC nnew from nnow +  dt * xtte
    IF (PRESENT(itrcr)) THEN
       CALL exchg_tens(itrcr)
       xt(:,:,itrcr,:) = xtm1(:,:,itrcr,:) + dt * xtte(:,:,itrcr,:)
    ELSE
       DO i=1, n_trcr
          CALL exchg_tens(i)
       END DO
       xt(:,:,:,:) = xtm1(:,:,:,:) + dt * xtte(:,:,:,:)
    ENDIF
 CASE DEFAULT

    CALL error_bi('trcr_calc', 'unknown flag')

 END SELECT

CONTAINS

  SUBROUTINE exchg_tens(itrcr)

    ! performs the boundary exchange between neighboring processors
    USE environment,      ONLY : exchg_boundaries
    USE data_parallel,    ONLY : num_compute, nboundlines, ldatatypes,    &
                                 ncomm_type, my_cart_neigh, icomm_cart,   &
                                 imp_reals, sendbuf, isendbuflen
    USE data_modelconfig, ONLY : ie, je, jstartpar, jendpar
    USE data_runcontrol , ONLY : lperi_x, lperi_y, l2dim, ntstep

    IMPLICIT NONE

    INTEGER(KIND=iintegers), INTENT(IN) :: itrcr

    ! LOCAL
    INTEGER(KIND=iintegers) :: kzdims(24) = 0
    INTEGER(KIND=iintegers) :: izerror
    CHARACTER(LEN=80)       :: yzerrmsg

    !  vertical dimensions for exchg_boundaries
!    kzdims(1:24) =  (/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    kzdims(1) = ke

    izerror = 0
    yzerrmsg = ' '

    CALL exchg_boundaries                                              &
         ( 0     , sendbuf, isendbuflen, imp_reals, icomm_cart,        &
         num_compute,ie, je,                                           &
         kzdims, jstartpar, jendpar, 2, nboundlines, my_cart_neigh,    &
         lperi_x, lperi_y, l2dim,                                      &
         20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,         &
         xtte(:,:,itrcr,:))

  END SUBROUTINE exchg_tens

#endif

END SUBROUTINE trcr_calc

!==============================================================================
!==============================================================================
!+ Module function "trcr_get_ntrcr" in "src_tracer" to get the
!                    number of tracers.
!------------------------------------------------------------------------------

PURE FUNCTION trcr_get_ntrcr()
!------------------------------------------------------------------------------
!
! Description:
!
! This function returns the number of tracers currently handled in the model.
!
!
! Method:
!
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

 INTEGER(KIND=iintegers) :: trcr_get_ntrcr ! number of tracers

!------------------------------------------------------------------------------
! Begin Subroutine trcr_get_ntrcr
!------------------------------------------------------------------------------

  ! return number of tracers
  trcr_get_ntrcr = n_trcr

!------------------------------------------------------------------------------
! End of module function trcr_get_ntrcr
!------------------------------------------------------------------------------

END FUNCTION trcr_get_ntrcr
!==============================================================================

!==============================================================================
!+ Module procedures to define metadata item
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!
! Description:
!   The following procedures overload the procedure trcr_meta_define which is
!   called upon a request to define a metadata item. For each data type a 
!   separate procedure  has to be created. The procedures all have exactly
!   the same structure with the exception of type of the ydefault
!   argument as well as some local type-dependent parameters. The procedures rely
!   heavily on the metadata_define function which do most of the work.
!
! Subroutines:
!
!   meta_def_integer
!   meta_def_real4
!   meta_def_real8
!   meta_def_string
!   meta_def_logical
!   meta_def_pointer2
!   meta_def_pointer3
!   meta_def_pointer4
!
! Error numbers: 1000-1099
!
! Method:
!   Check define request
!   Define metadata entry
!
!==============================================================================

SUBROUTINE trcr_meta_def_integer(ierr,yname,ydefault,iidx,lprotect)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: yname
  INTEGER                , INTENT(IN)  :: ydefault
  INTEGER, OPTIONAL      , INTENT(OUT) :: iidx
  LOGICAL, OPTIONAL      , INTENT(IN)  :: lprotect
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_define(mtrcr, yname, ydefault, iidx, lprotect)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
!CALL info_bi('no definition of metadata '//TRIM(yname)//' possible in MESSy','')
#endif
END SUBROUTINE trcr_meta_def_integer

SUBROUTINE trcr_meta_def_real4(ierr,yname,ydefault,iidx,lprotect)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: yname
  REAL(KIND=sp)          , INTENT(IN)  :: ydefault
  INTEGER, OPTIONAL      , INTENT(OUT) :: iidx
  LOGICAL, OPTIONAL      , INTENT(IN)  :: lprotect
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_define(mtrcr, yname, ydefault, iidx, lprotect)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
!CALL info_bi('no definition of metadata '//TRIM(yname)//' possible in MESSy','')
#endif
END SUBROUTINE trcr_meta_def_real4

SUBROUTINE trcr_meta_def_real8(ierr,yname,ydefault,iidx,lprotect)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: yname
  REAL(KIND=dp)          , INTENT(IN)  :: ydefault
  INTEGER, OPTIONAL      , INTENT(OUT) :: iidx
  LOGICAL, OPTIONAL      , INTENT(IN)  :: lprotect
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_define(mtrcr, yname, ydefault, iidx, lprotect)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
!CALL info_bi('no definition of metadata '//TRIM(yname)//' possible in MESSy','')
#endif
END SUBROUTINE trcr_meta_def_real8

SUBROUTINE trcr_meta_def_string(ierr,yname,ydefault,iidx,lprotect)
  IMPLICIT NONE
  INTEGER(KIND=iintegers)    , INTENT(OUT) :: ierr
  CHARACTER(LEN=*)           , INTENT(IN)  :: yname
  CHARACTER(LEN=*)           , INTENT(IN)  :: ydefault
  CHARACTER(LEN=I_MAX_STRLEN)              :: yzstr
  INTEGER, OPTIONAL          , INTENT(OUT) :: iidx
  LOGICAL, OPTIONAL          , INTENT(IN)  :: lprotect
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_define(mtrcr, yname, ydefault, iidx, lprotect)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
!CALL info_bi('no definition of metadata '//TRIM(yname)//'p ossible in MESSy','')
#endif
END SUBROUTINE trcr_meta_def_string

SUBROUTINE trcr_meta_def_logical(ierr,yname,ydefault,iidx,lprotect)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: yname
  LOGICAL                , INTENT(IN)  :: ydefault
  INTEGER, OPTIONAL      , INTENT(OUT) :: iidx
  LOGICAL, OPTIONAL      , INTENT(IN)  :: lprotect
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_define(mtrcr, yname, ydefault, iidx, lprotect)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
!CALL info_bi('no definition of metadata '//TRIM(yname)//' possible in MESSy','')
#endif
END SUBROUTINE trcr_meta_def_logical

SUBROUTINE trcr_meta_def_pointer2(ierr,yname,ydefault,iidx,lprotect)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: yname
  REAL(KIND=wp) ,POINTER , INTENT(IN)  :: ydefault(:,:)
  INTEGER, OPTIONAL      , INTENT(OUT) :: iidx
  LOGICAL, OPTIONAL      , INTENT(IN)  :: lprotect
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_define(mtrcr, yname, ydefault, iidx, lprotect)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
!CALL info_bi('no definition of metadata '//TRIM(yname)//' possible in MESSy','')
#endif
END SUBROUTINE trcr_meta_def_pointer2

SUBROUTINE trcr_meta_def_pointer3(ierr,yname,ydefault,iidx,lprotect)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: yname
  REAL(KIND=wp), POINTER , INTENT(IN)  :: ydefault(:,:,:)
  INTEGER, OPTIONAL      , INTENT(OUT) :: iidx
  LOGICAL, OPTIONAL      , INTENT(IN)  :: lprotect
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_define(mtrcr, yname, ydefault, iidx, lprotect)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
!CALL info_bi('no definition of metadata '//TRIM(yname)//' possible in MESSy','')
#endif
END SUBROUTINE trcr_meta_def_pointer3

SUBROUTINE trcr_meta_def_pointer4(ierr,yname,ydefault,iidx,lprotect)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: yname
  REAL(KIND=wp), POINTER , INTENT(IN)  :: ydefault(:,:,:,:)
  INTEGER, OPTIONAL      , INTENT(OUT) :: iidx
  LOGICAL, OPTIONAL      , INTENT(IN)  :: lprotect
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_define(mtrcr, yname, ydefault, iidx, lprotect)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
!CALL info_bi('no definition of metadata '//TRIM(yname)//' possible in MESSy','')
#endif
END SUBROUTINE trcr_meta_def_pointer4

!==============================================================================

!==============================================================================
!+ Module procedures to get data from metadata storage
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!
! Description:
!   The following procedures overload the procedure trcr_meta_get which is
!   called upon a request to retrieve a data item from the metadata storage. For
!   each data type as well as for scalar/vector data a separate procedure
!   has to be created. The procedures all have exactly the same structure with
!   the exception of dimension and type of the ydata argument as well as
!   some local type-dependent parameters. The procedures rely heavily on the
!   trcr_meta_get_check function which do most of the work.
!
! Subroutines:
!
!  trcr_meta_get_int_byindind      trcr_meta_get_r4_byindind
!  trcr_meta_get_r8_byindind       trcr_meta_get_str_byindind
!  trcr_meta_get_log_byindind      trcr_meta_get_ptr2_byindind
!  trcr_meta_get_ptr3_byindind     trcr_meta_get_ptr4_byindind
!  trcr_meta_get_int_byindname     trcr_meta_get_r4_byindname
!  trcr_meta_get_r8_byindname      trcr_meta_get_str_byindname
!  trcr_meta_get_log_byindname     trcr_meta_get_ptr2_byindname
!  trcr_meta_get_ptr3_byindname    trcr_meta_get_ptr4_byindname
!  trcr_meta_get_int_bynameind     trcr_meta_get_r4_bynameind
!  trcr_meta_get_r8_bynameind      trcr_meta_get_str_bynameind
!  trcr_meta_get_log_bynameind     trcr_meta_get_ptr2_bynameind
!  trcr_meta_get_ptr3_bynameind    trcr_meta_get_ptr4_bynameind
!  trcr_meta_get_int_bynamename    trcr_meta_get_r4_bynamename
!  trcr_meta_get_r8_bynamename     trcr_meta_get_str_bynamename
!  trcr_meta_get_log_bynamename    trcr_meta_get_ptr2_bynamename
!  trcr_meta_get_ptr3_bynamename    trcr_meta_get_ptr4_bynamename
!  trcr_meta_get_intv_byindex      trcr_meta_get_r4v_byindex
!  trcr_meta_get_r8v_byindex       trcr_meta_get_strv_byindex
!  trcr_meta_get_logv_byindex
!  trcr_meta_get_intv_byname       trcr_meta_get_r4v_byname
!  trcr_meta_get_r8v_byname        trcr_meta_get_strv_byname
!  trcr_meta_get_logv_byname
!
! Error numbers: 800-899
!
! Method:
!
!==============================================================================

SUBROUTINE trcr_meta_get_int_byindind(ierr,itrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  INTEGER                , INTENT(IN)  :: imeta
  INTEGER                , INTENT(OUT) :: ydata
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    IF (imeta < 0 .OR. imeta > MAX_CASK_I) THEN
       CALL error_bi(' meta data index not defined in MESSy' &
            , 'trcr_meta_get_int_byindind' )
    ENDIF
    IF (imeta == T_ADV_ID) THEN
       IF (ti_gp(itrcr)%tp%meta%cask_i(imeta) >= 1) THEN
          ydata = T_ADV_ON
       ELSE
          ydata = T_ADV_OFF
       ENDIF
    ELSE IF (imeta == T_DAMP_ID) THEN
       IF (ti_gp(itrcr)%tp%meta%cask_i(imeta) >= 1) THEN
          ydata = T_DAMP_ON
       ELSE
          ydata = T_DAMP_OFF
       ENDIF
    ELSE
       ydata = ti_gp(itrcr)%tp%meta%cask_i(imeta)
    ENDIF
#endif
  ENDIF
END SUBROUTINE trcr_meta_get_int_byindind

SUBROUTINE trcr_meta_get_r4_byindind(ierr,itrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  INTEGER                , INTENT(IN)  :: imeta
  REAL(KIND=sp)          , INTENT(OUT) :: ydata
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    IF (imeta < 0 .OR. imeta > MAX_CASK_R) THEN
       CALL error_bi(' meta data index not defined in MESSy' &
            , 'trcr_meta_get_r4_byindind' )
    ENDIF
    ydata = REAL(ti_gp(itrcr)%tp%meta%cask_r(imeta),wp)
#endif
  ENDIF
END SUBROUTINE trcr_meta_get_r4_byindind

SUBROUTINE trcr_meta_get_r8_byindind(ierr,itrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  INTEGER                , INTENT(IN)  :: imeta
  REAL(KIND=dp)          , INTENT(OUT) :: ydata
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    IF (imeta < 0 .OR. imeta > MAX_CASK_R) THEN
       CALL error_bi(' meta data index not defined in MESSy' &
            , 'trcr_meta_get_r8_byindind' )
    ENDIF
    ydata = ti_gp(itrcr)%tp%meta%cask_r(imeta)
#endif
  ENDIF
END SUBROUTINE trcr_meta_get_r8_byindind

SUBROUTINE trcr_meta_get_str_byindind(ierr,itrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  INTEGER                , INTENT(IN)  :: imeta
  CHARACTER(LEN=*)       , INTENT(OUT) :: ydata
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    IF (imeta < 0 .OR. imeta > MAX_CASK_S) THEN
       SELECT CASE(imeta)
       CASE(-1, -4, -5) ! name
          ydata = TRIM(ti_gp(itrcr)%tp%ident%fullname)
       CASE (-2)
          ydata = TRIM(ti_gp(itrcr)%tp%ident%submodel)
       CASE (-3)
          ydata = TRIM(ti_gp(itrcr)%tp%ident%unit)
       CASE DEFAULT
          !CALL info_bi(&
          !'Meta data index unknown in MESSy retruning empty string' &
          !           , 'trcr_meta_get_str_byindind' )
          ydata(:) = ' '
       END SELECT
    ELSE
       ydata = ti_gp(itrcr)%tp%meta%cask_s(imeta)
    ENDIF
#endif
  ENDIF
END SUBROUTINE trcr_meta_get_str_byindind

SUBROUTINE trcr_meta_get_log_byindind(ierr,itrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  INTEGER                , INTENT(IN)  :: imeta
  LOGICAL                , INTENT(OUT) :: ydata
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    !CALL info_bi('logical meta data not available in MESSY', '')
    ierr = 5000_iintegers
#endif
  ENDIF
END SUBROUTINE trcr_meta_get_log_byindind

SUBROUTINE trcr_meta_get_ptr2_byindind(ierr,itrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  INTEGER                , INTENT(IN)  :: imeta
  REAL(KIND=wp), POINTER , INTENT(OUT) :: ydata(:,:)
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    !CALL info_bi('Due to different concepts: No 2D-meta data avaiable in MESSy'&
    !         , '' )
    ierr = 5000_iintegers
#endif
  ENDIF
END SUBROUTINE trcr_meta_get_ptr2_byindind

SUBROUTINE trcr_meta_get_ptr3_byindind(ierr,itrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  INTEGER                , INTENT(IN)  :: imeta
  REAL(KIND=wp), POINTER , INTENT(OUT) :: ydata(:,:,:)
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    !CALL info_bi('Due to different concepts: No 3D-meta data avaiable in MESSy'&
    !         , '' )
    ierr = 5000_iintegers
#endif
  ENDIF
END SUBROUTINE trcr_meta_get_ptr3_byindind

SUBROUTINE trcr_meta_get_ptr4_byindind(ierr,itrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  INTEGER                , INTENT(IN)  :: imeta
  REAL(KIND=wp), POINTER , INTENT(OUT) :: ydata(:,:,:,:)
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    !CALL info_bi('Due to different concepts: No 4D-meta data avaiable in MESSy'&
    !         , '' )
    ierr = 5000_iintegers
#endif
  ENDIF
END SUBROUTINE trcr_meta_get_ptr4_byindind

SUBROUTINE trcr_meta_get_int_byindname(ierr,itrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  INTEGER                , INTENT(OUT) :: ydata
  INTEGER                              :: imeta
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    imeta = func_imeta(ymeta)
    IF (imeta > 0 .AND. imeta <= MAX_CASK_I) THEN
       IF (imeta == T_ADV_ID) THEN
          IF (ti_gp(itrcr)%tp%meta%cask_i(imeta) >= 1) THEN
             IF (TRIM(ymeta) == 'SP_ADV_LF') THEN
                ydata = ti_gp(itrcr)%tp%meta%cask_i(imeta)
             ELSE
                ydata = T_ADV_ON
             END IF
          ELSE
             ydata = T_ADV_OFF
          ENDIF
       ELSE IF (imeta == T_DAMP_ID) THEN
          IF (ti_gp(itrcr)%tp%meta%cask_i(imeta) >= 1) THEN
             ydata = T_DAMP_ON
          ELSE
             ydata = T_DAMP_OFF
          ENDIF
       ELSE
          ydata = ti_gp(itrcr)%tp%meta%cask_i(imeta)
       END IF
    ELSE
       ! CALL info_bi('Meta data name '//TRIM(ymeta)//' unknown in MESSy' &
       !            , 'trcr_meta_get_int_byindname' )
       !       ierr = 5000_iintegers
    ENDIF
#endif
  ENDIF
END SUBROUTINE trcr_meta_get_int_byindname

SUBROUTINE trcr_meta_get_r4_byindname(ierr,itrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  REAL(KIND=sp)          , INTENT(OUT) :: ydata
  INTEGER                              :: imeta
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    imeta = func_imeta(ymeta)
    IF (imeta > 0 .AND. imeta <= MAX_CASK_R) THEN
       ydata = REAL(ti_gp(itrcr)%tp%meta%cask_r(imeta),wp)
    ELSE
       !CALL info_bi('Meta data name '//TRIM(ymeta)//' unknown in MESSy' &
       !         , 'trcr_meta_get_r4_byindname' )
       !ierr = 5000_iintegers
    ENDIF
#endif
  ENDIF
END SUBROUTINE trcr_meta_get_r4_byindname

SUBROUTINE trcr_meta_get_r8_byindname(ierr,itrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  REAL(KIND=dp)          , INTENT(OUT) :: ydata
  INTEGER                              :: imeta
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    imeta = func_imeta(ymeta)
    IF (imeta > 0 .AND. imeta <= MAX_CASK_R) THEN
       ydata = ti_gp(itrcr)%tp%meta%cask_r(imeta)
    ELSE
       !CALL info_bi('Meta data name  '//TRIM(ymeta)//' unknown in MESSy' &
       !           , 'trcr_meta_get_r8_byindname' )
       !ierr = 5000_iintegers
    ENDIF
#endif
  ENDIF
END SUBROUTINE trcr_meta_get_r8_byindname

SUBROUTINE trcr_meta_get_str_byindname(ierr,itrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  CHARACTER(LEN=*)       , INTENT(OUT) :: ydata
  INTEGER                              :: imeta
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    imeta = func_imeta(ymeta)
    IF (imeta > 0 .AND. imeta <= MAX_CASK_S) THEN
       ydata = ti_gp(itrcr)%tp%meta%cask_s(imeta)
    ELSE
       SELECT CASE(imeta)
       CASE(-1, -4, -5) ! name
          ydata = TRIM(ti_gp(itrcr)%tp%ident%fullname)
       CASE (-2)
          ydata = TRIM(ti_gp(itrcr)%tp%ident%submodel)
       CASE (-3)
          ydata = TRIM(ti_gp(itrcr)%tp%ident%unit)
       CASE DEFAULT
          !  CALL info_bi('Meta data index unknown '//&
          ! //TRIM(ymeta)//' in MESSy retruning empty string' &
          !                , 'trcr_meta_get_str_byindname' )
          ydata(:) = ' '
       END SELECT
    ENDIF
#endif
  ENDIF
END SUBROUTINE trcr_meta_get_str_byindname

SUBROUTINE trcr_meta_get_log_byindname(ierr,itrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  LOGICAL                , INTENT(OUT) :: ydata
  INTEGER                              :: imeta
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    imeta = func_imeta(ymeta)
    IF (imeta > 0 .AND. imeta <= MAX_CASK_I) THEN
       IF (ti_gp(itrcr)%tp%meta%cask_i(imeta) == OFF) THEN
          ydata = .FALSE.
       ELSE
          ydata = .TRUE.
       ENDIF
    ELSE
       ! CALL info_bi('Meta data name '//TRIM(ymeta)//' unknown in MESSy' &
       !            , 'trcr_meta_get_log_byindname' )
       !ierr = 5000_iintegers
    ENDIF
#endif
  ENDIF
END SUBROUTINE trcr_meta_get_log_byindname

SUBROUTINE trcr_meta_get_ptr2_byindname(ierr,itrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  REAL(KIND=wp), POINTER , INTENT(OUT) :: ydata(:,:)
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    !CALL info_bi('Due to different concepts: No 2D-meta data avaiable in MESSy'&
    !         , '' )
    ierr = 5000_iintegers
#endif
  ENDIF
END SUBROUTINE trcr_meta_get_ptr2_byindname

SUBROUTINE trcr_meta_get_ptr3_byindname(ierr,itrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  REAL(KIND=wp), POINTER , INTENT(OUT) :: ydata(:,:,:)
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    !CALL info_bi('Due to different concepts: No 3D-meta data avaiable in MESSy'&
    !         , '' )
    ierr = 5000_iintegers
    IF (itrcr == idt_qv .AND. ymeta=='SURF_FIELD') THEN
       ierr = 0_iintegers
       ydata => qv_s
    END IF
#endif
  ENDIF
END SUBROUTINE trcr_meta_get_ptr3_byindname

SUBROUTINE trcr_meta_get_ptr4_byindname(ierr,itrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  REAL(KIND=wp), POINTER , INTENT(OUT) :: ydata(:,:,:,:)
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    !CALL info_bi('Due to different concepts: No 4D-meta data avaiable in MESSy'&
    !         , '' )
    ierr = 5000_iintegers
#endif
  ENDIF
END SUBROUTINE trcr_meta_get_ptr4_byindname

SUBROUTINE trcr_meta_get_int_bynameind(ierr,ytrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  INTEGER(KIND=iintegers), INTENT(IN)  :: imeta
  INTEGER                , INTENT(OUT) :: ydata
  INTEGER                              :: itrcr
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    IF (imeta < 0 .OR. imeta > MAX_CASK_I) THEN
       CALL error_bi(' meta data index not defined in MESSy' &
            , 'trcr_meta_get_int_bynameind' )
    ENDIF
    ydata = ti_gp(itrcr)%tp%meta%cask_i(imeta)
#endif
  END IF
END SUBROUTINE trcr_meta_get_int_bynameind

SUBROUTINE trcr_meta_get_r4_bynameind(ierr,ytrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  INTEGER(KIND=iintegers), INTENT(IN)  :: imeta
  REAL(KIND=sp)          , INTENT(OUT) :: ydata
  INTEGER                              :: itrcr
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    IF (imeta < 0 .OR. imeta > MAX_CASK_R) THEN
       CALL error_bi(' meta data index not defined in MESSy' &
            , 'trcr_meta_get_r4_bynameind' )
    ENDIF
    ydata = REAL(ti_gp(itrcr)%tp%meta%cask_r(imeta),wp)
#endif
  END IF
END SUBROUTINE trcr_meta_get_r4_bynameind

SUBROUTINE trcr_meta_get_r8_bynameind(ierr,ytrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  INTEGER(KIND=iintegers), INTENT(IN)  :: imeta
  REAL(KIND=dp)          , INTENT(OUT) :: ydata
  INTEGER                              :: itrcr
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    IF (imeta < 0 .OR. imeta > MAX_CASK_R) THEN
       CALL error_bi(' meta data index not defined in MESSy' &
            , 'trcr_meta_get_r8_bynameind' )
    ENDIF
    ydata = ti_gp(itrcr)%tp%meta%cask_r(imeta)
#endif
  END IF
END SUBROUTINE trcr_meta_get_r8_bynameind

SUBROUTINE trcr_meta_get_str_bynameind(ierr,ytrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  INTEGER(KIND=iintegers), INTENT(IN)  :: imeta
  CHARACTER(LEN=*)       , INTENT(OUT) :: ydata
  INTEGER                              :: itrcr
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    IF (imeta < 0 .OR. imeta > MAX_CASK_S) THEN
       SELECT CASE(imeta)
       CASE(-1, -4, -5) ! name
           ydata = TRIM(ti_gp(itrcr)%tp%ident%fullname)
        CASE (-2)
           ydata = TRIM(ti_gp(itrcr)%tp%ident%submodel)
        CASE (-3)
           ydata = TRIM(ti_gp(itrcr)%tp%ident%unit)
        CASE DEFAULT
           ! CALL info_bi( &
           ! 'Meta data index unknown in MESSy retruning empty string' &
           !    , 'trcr_meta_get_str_bynameind' )
           ydata(:) = ' '
        END SELECT
     ELSE
        ydata = ti_gp(itrcr)%tp%meta%cask_s(imeta)
     ENDIF
#endif
  END IF
END SUBROUTINE trcr_meta_get_str_bynameind

SUBROUTINE trcr_meta_get_log_bynameind(ierr,ytrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  INTEGER(KIND=iintegers), INTENT(IN)  :: imeta
  LOGICAL                , INTENT(OUT) :: ydata
  INTEGER                              :: itrcr
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    IF (imeta < 0 .OR. imeta > MAX_CASK_I) THEN
       CALL error_bi(' meta data index not defined in MESSy' &
            , 'trcr_meta_get_log_bynameind' )
    ENDIF
    ! CALL info_bi('logical meta data not available in MESSY', '')
    ierr = 5000_iintegers
#endif
  END IF
END SUBROUTINE trcr_meta_get_log_bynameind

SUBROUTINE trcr_meta_get_ptr2_bynameind(ierr,ytrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  INTEGER(KIND=iintegers), INTENT(IN)  :: imeta
  REAL(KIND=wp), POINTER , INTENT(OUT) :: ydata(:,:)
  INTEGER :: itrcr
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    !CALL info_bi('Due to different concepts: No 2D-meta data avaiable in MESSy'&
    !         , '' )
    ierr = 5000_iintegers
#endif
  ENDIF
END SUBROUTINE trcr_meta_get_ptr2_bynameind

SUBROUTINE trcr_meta_get_ptr3_bynameind(ierr,ytrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  INTEGER(KIND=iintegers), INTENT(IN)  :: imeta
  REAL(KIND=wp), POINTER , INTENT(OUT) :: ydata(:,:,:)
  INTEGER :: itrcr
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    !CALL info_bi('Due to different concepts: No 3D-meta data avaiable in MESSy'&
    !         , '' )
    ierr = 5000_iintegers
#endif
  ENDIF
END SUBROUTINE trcr_meta_get_ptr3_bynameind

SUBROUTINE trcr_meta_get_ptr4_bynameind(ierr,ytrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  INTEGER(KIND=iintegers), INTENT(IN)  :: imeta
  REAL(KIND=wp), POINTER , INTENT(OUT) :: ydata(:,:,:,:)
  INTEGER :: itrcr
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    !CALL info_bi('Due to different concepts: No 4D-meta data avaiable in MESSy'&
    !           , '' )
    ierr = 5000_iintegers
#endif
  ENDIF
END SUBROUTINE trcr_meta_get_ptr4_bynameind

SUBROUTINE trcr_meta_get_int_bynamename(ierr,ytrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  INTEGER                , INTENT(OUT) :: ydata
  INTEGER                              :: itrcr
  INTEGER                              :: imeta
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    imeta = func_imeta(ymeta)
    IF (imeta > 0 .AND. imeta <= MAX_CASK_I) THEN
       IF (imeta == T_ADV_ID) THEN
          IF (ti_gp(itrcr)%tp%meta%cask_i(imeta) >= 1) THEN
             IF (TRIM(ymeta) == 'SP_ADV_LF') THEN
                ydata = ti_gp(itrcr)%tp%meta%cask_i(imeta)
             ELSE
                ydata = T_ADV_ON
             END IF
          ELSE
             ydata = T_ADV_OFF
          ENDIF
       ELSE IF (imeta == T_DAMP_ID) THEN
          IF (ti_gp(itrcr)%tp%meta%cask_i(imeta) >= 1) THEN
             ydata = T_DAMP_ON
          ELSE
             ydata = T_DAMP_OFF
          ENDIF
       ELSE
          ydata = ti_gp(itrcr)%tp%meta%cask_i(imeta)
       END IF
    ELSE
       ! CALL info_bi('Meta data name '//TRIM(ymeta)//' unknown in MESSy' &
       !            , 'trcr_meta_get_int_bynamename' )
       ! ierr = 5000_iintegers
    ENDIF
#endif
  END IF
END SUBROUTINE trcr_meta_get_int_bynamename

SUBROUTINE trcr_meta_get_r4_bynamename(ierr,ytrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  REAL(KIND=sp)          , INTENT(OUT) :: ydata
  INTEGER                              :: itrcr
  INTEGER                              :: imeta
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    imeta = func_imeta(ymeta)
    IF (imeta > 0 .AND. imeta <= MAX_CASK_R) THEN
       ydata = REAL(ti_gp(itrcr)%tp%meta%cask_i(imeta),wp)
    ELSE
       !CALL info_bi('Meta data name '//TRIM(ymeta)//' unknown in MESSy' &
       !           , 'trcr_meta_get_r4_bynamename' )
       !ierr = 5000_iintegers
    ENDIF
#endif
  END IF
END SUBROUTINE trcr_meta_get_r4_bynamename

SUBROUTINE trcr_meta_get_r8_bynamename(ierr,ytrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  REAL(KIND=dp)          , INTENT(OUT) :: ydata
  INTEGER                              :: itrcr
  INTEGER                              :: imeta
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    imeta = func_imeta(ymeta)
    IF (imeta > 0 .AND. imeta <= MAX_CASK_R) THEN
       ydata = ti_gp(itrcr)%tp%meta%cask_r(imeta)
    ELSE
       ! CALL info_bi('Meta data name '//TRIM(ymeta)//' unknown in MESSy' &
       !            , 'trcr_meta_get_r8_bynamename' )
       !ierr = 5000_iintegers
    ENDIF
#endif
  END IF
END SUBROUTINE trcr_meta_get_r8_bynamename

SUBROUTINE trcr_meta_get_str_bynamename(ierr,ytrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  CHARACTER(LEN=*)       , INTENT(OUT) :: ydata
  INTEGER                              :: itrcr
  INTEGER                              :: imeta
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    imeta = func_imeta(ymeta)
    IF (imeta > 0 .AND. imeta <= MAX_CASK_S) THEN
       ydata = TRIM(ti_gp(itrcr)%tp%meta%cask_s(imeta))
    ELSE
       SELECT CASE(imeta)
       CASE(-1, -4, -5) ! name
           ydata = TRIM(ti_gp(itrcr)%tp%ident%fullname)
        CASE (-2)
           ydata = TRIM(ti_gp(itrcr)%tp%ident%submodel)
        CASE (-3)
           ydata = TRIM(ti_gp(itrcr)%tp%ident%unit)
        CASE DEFAULT
           ! CALL info_bi( &
           ! 'Meta data name '//TRIM(ymeta)//' unknown in MESSy retruning empty string' &
           !                , 'trcr_meta_get_str_bynamename' )
           ydata(:) = ' '
        END SELECT
     ENDIF
#endif
  END IF
END SUBROUTINE trcr_meta_get_str_bynamename

SUBROUTINE trcr_meta_get_log_bynamename(ierr,ytrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  LOGICAL                , INTENT(OUT) :: ydata
  INTEGER                              :: itrcr
  INTEGER                              :: imeta
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    imeta = func_imeta(ymeta)
    IF (imeta > 0 .AND. imeta <= MAX_CASK_I) THEN
       IF (ti_gp(itrcr)%tp%meta%cask_i(imeta) == OFF) THEN
          ydata = .FALSE.
       ELSE
          ydata = .TRUE.
       ENDIF
    ELSE
       ! CALL info_bi('Meta data name '//TRIM(ymeta)//' unknown in MESSy' &
       !            , 'trcr_meta_get_log_byindname' )
       !ierr = 5000_iintegers
    ENDIF
#endif
  END IF
END SUBROUTINE trcr_meta_get_log_bynamename

SUBROUTINE trcr_meta_get_ptr2_bynamename(ierr,ytrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  REAL(KIND=wp), POINTER , INTENT(OUT) :: ydata(:,:)
  INTEGER :: itrcr
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    !CALL info_bi('Due to different concepts: No 2D-meta data avaiable in MESSy'&
    !           , '' )
    ierr = 5000_iintegers
#endif
  ENDIF
END SUBROUTINE trcr_meta_get_ptr2_bynamename

SUBROUTINE trcr_meta_get_ptr3_bynamename(ierr,ytrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  REAL(KIND=wp), POINTER , INTENT(OUT) :: ydata(:,:,:)
  INTEGER :: itrcr
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    !CALL info_bi('Due to different concepts: No 3D-meta data avaiable in MESSy'&
    !           , '' )
    ierr = 5000_iintegers
#endif
  ENDIF
END SUBROUTINE trcr_meta_get_ptr3_bynamename

SUBROUTINE trcr_meta_get_ptr4_bynamename(ierr,ytrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  REAL(KIND=wp), POINTER , INTENT(OUT) :: ydata(:,:,:,:)
  INTEGER :: itrcr
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_get(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    !CALL info_bi('Due to different concepts: No 4D-meta data avaiable in MESSy'&
    !         , '' )
    ierr = 5000_iintegers
#endif
  ENDIF
END SUBROUTINE trcr_meta_get_ptr4_bynamename

SUBROUTINE trcr_meta_get_intv_byindex(ierr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER, INTENT(IN) :: imeta
  INTEGER, INTENT(OUT) :: ydata(:)
  ! LOCAL
  INTEGER :: i
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_get(mtrcr,imeta,ydata)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers 
#else
  IF (imeta < 0 .OR. imeta > MAX_CASK_I) THEN
     ! special case for qv_s surface flux
     !     CALL info_bi(' meta data index not defined in MESSy' &
     !          , 'trcr_meta_get_intv_byindex' )
     !     CALL info_bi('returning OFF', 'trcr_meta_get_intv_byindex'  )
     ydata(:) = 0
     IF (imeta==T_BBC_ID) THEN
        DO i=1,n_trcr
           IF (TRIM(ti_gp(i)%tp%ident%fullname) == 'QV') THEN
              ydata(i) = T_BBC_SURF_VAL
              EXIT
           ENDIF
        END DO
     ENDIF
  ELSE IF (imeta == T_ADV_ID) THEN
     DO i=1,n_trcr
        IF (ti_gp(i)%tp%meta%cask_i(imeta) >= 1) THEN
           ydata(i) = T_ADV_ON
        ELSE
           ydata(i) = T_ADV_OFF
        ENDIF
     END DO
  ELSE IF (imeta == T_DAMP_ID) THEN
     DO i=1,n_trcr
        IF (ti_gp(i)%tp%meta%cask_i(imeta) >= 1) THEN
           ydata(i) = T_DAMP_ON
        ELSE
           ydata(i) = T_DAMP_OFF
        ENDIF
     END DO
  ELSE
     DO i=1,n_trcr
        ydata(i) = ti_gp(i)%tp%meta%cask_i(imeta)
     END DO
  ENDIF
#endif
END SUBROUTINE trcr_meta_get_intv_byindex

SUBROUTINE trcr_meta_get_r4v_byindex(ierr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: imeta
  REAL(KIND=sp)          , INTENT(OUT) :: ydata(:)
  INTEGER                              :: i
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_get(mtrcr,imeta,ydata)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers 
#else
  IF (imeta < 0 .OR. imeta > MAX_CASK_R) THEN
     !     CALL info_bi(' meta data index not defined in MESSy' &
     !          , 'trcr_meta_get_r4v_byindex' )
     !     CALL info_bi('returning zero', 'trcr_meta_get_r4v_byindex'  )
     ydata(:) = 0.0_wp
  ELSE
     DO i=1,n_trcr
        ydata(i) = REAL(ti_gp(i)%tp%meta%cask_r(imeta),wp)
     END DO
  ENDIF
#endif
END SUBROUTINE trcr_meta_get_r4v_byindex

SUBROUTINE trcr_meta_get_r8v_byindex(ierr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: imeta
  REAL(KIND=dp)          , INTENT(OUT) :: ydata(:)
  INTEGER                              :: i
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_get(mtrcr,imeta,ydata)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
  IF (imeta < 0 .OR. imeta > MAX_CASK_R) THEN
     !     CALL info_bi(' meta data index not defined in MESSy' &
     !          , 'trcr_meta_get_r8v_byindex' )
     !     CALL info_bi('returning zero', 'trcr_meta_get_r8v_byindex'  )
     ydata(:) = 0.0_wp
  ELSE
     DO i=1,n_trcr
        ydata(i) = ti_gp(i)%tp%meta%cask_r(imeta)
     END DO
  ENDIF
#endif
END SUBROUTINE trcr_meta_get_r8v_byindex

SUBROUTINE trcr_meta_get_strv_byindex(ierr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: imeta
  CHARACTER(LEN=*)       , INTENT(OUT) :: ydata(:)
  INTEGER                              :: i
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_get(mtrcr,imeta,ydata)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
  IF (imeta < 0 .OR. imeta > MAX_CASK_S) THEN
     SELECT CASE(imeta)
     CASE(-1, -4, -5) ! name
        DO i=1,n_trcr
           ydata(i) = TRIM(ti_gp(i)%tp%ident%fullname)
        END DO
     CASE (-2)
        DO i=1,n_trcr
           ydata(i) = TRIM(ti_gp(i)%tp%ident%submodel)
        END DO
     CASE (-3)
        DO i=1,n_trcr
           ydata(i) = TRIM(ti_gp(i)%tp%ident%unit)
        END DO
     CASE DEFAULT
        !  CALL info_bi(' meta data index not defined in MESSy' &
        !             , 'trcr_meta_get_strv_byindex' )
        !  CALL info_bi('returning empty string', 'trcr_meta_get_strv_byindex'  )
        ydata(:) = ' '
     END SELECT
  ELSE
     DO i=1,n_trcr
        ydata(i) = TRIM(ti_gp(i)%tp%meta%cask_s(imeta))
     END DO
  ENDIF
#endif
END SUBROUTINE trcr_meta_get_strv_byindex

SUBROUTINE trcr_meta_get_logv_byindex(ierr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: imeta
  LOGICAL                , INTENT(OUT) :: ydata(:)
  INTEGER                              :: i
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_get(mtrcr,imeta,ydata)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
  IF (imeta < 0 .OR. imeta > MAX_CASK_I) THEN
     !     CALL info_bi(' meta data index not defined in MESSy' &
     !                , 'trcr_meta_get_logv_byindex' )
     !     CALL info_bi('returning .FALSE.', 'trcr_meta_get_logv_byindex'  )
     ydata(:) = .FALSE.
  ELSE
     DO i=1,n_trcr
        IF (ti_gp(i)%tp%meta%cask_i(imeta) == OFF) THEN
           ydata(i) = .FALSE.
        ELSE
           ydata(i) = .TRUE.
        END IF
     END DO
  ENDIF
#endif
END SUBROUTINE trcr_meta_get_logv_byindex

SUBROUTINE trcr_meta_get_intv_byname(ierr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  INTEGER                , INTENT(OUT) :: ydata(:)
  INTEGER                              :: imeta
  INTEGER                              :: i
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_get(mtrcr,ymeta,ydata)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
  imeta = func_imeta(ymeta)
  IF (imeta > 0 .AND. imeta <= MAX_CASK_I) THEN
     IF (imeta == I_ADVECT) THEN
        IF (TRIM(ymeta) == 'SP_ADV_LF') THEN
           DO i=1,n_trcr
              ydata(i) = ti_gp(i)%tp%meta%cask_i(imeta)
           END DO
        ELSE
           DO i=1,n_trcr
              IF (ti_gp(i)%tp%meta%cask_i(imeta) >= 1) THEN
                 ydata = T_ADV_ON
              ELSE
                 ydata = T_ADV_OFF
              END IF
           END DO
        ENDIF
     ELSE IF (imeta == I_DAMP) THEN
         DO i=1,n_trcr
            IF (ti_gp(i)%tp%meta%cask_i(imeta) >= 1) THEN
               ydata(i) = T_DAMP_ON
            ELSE
               ydata(i) = T_DAMP_OFF
            ENDIF

         END DO
     ELSE
        DO i=1,n_trcr
           ydata(i) = ti_gp(i)%tp%meta%cask_i(imeta)
        END DO
     ENDIF
  ELSE
     ! CALL info_bi('Meta data name '//TRIM(ymeta)//' unknown in MESSy' &
     !            , 'trcr_meta_get_intv_byname' )
     !     CALL info_bi('returning OFF', 'trcr_meta_get_intv_byname'  )
     ydata(:) = 0
  ENDIF
#endif
END SUBROUTINE trcr_meta_get_intv_byname

SUBROUTINE trcr_meta_get_r4v_byname(ierr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  REAL(KIND=sp)          , INTENT(OUT) :: ydata(:)
  INTEGER                              :: imeta
  INTEGER                              :: i
#ifndef MESSY
  ierr = 0_iintegers
  CALL metadata_get(mtrcr,ymeta,ydata)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
  imeta = func_imeta(ymeta)
  IF (imeta > 0 .AND. imeta <= MAX_CASK_R) THEN
     DO i = 1, n_trcr
        ydata(i) = REAL(ti_gp(i)%tp%meta%cask_r(imeta),wp)
     END DO
  ELSE
     !     CALL info_bi('Meta data name '//TRIM(ymeta)//' unknown in MESSy' &
     !                , 'trcr_meta_get_r4v_byname' )
     !     CALL info_bi('returning zero', 'trcr_meta_get_r4v_byname'  )
     ydata(:) = 0.0_wp

  ENDIF
#endif
END SUBROUTINE trcr_meta_get_r4v_byname

SUBROUTINE trcr_meta_get_r8v_byname(ierr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  REAL(KIND=dp)          , INTENT(OUT) :: ydata(:)
  INTEGER                              :: imeta
  INTEGER                              :: i
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_get(mtrcr,ymeta,ydata)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
  imeta = func_imeta(ymeta)
  IF (imeta > 0 .AND. imeta <= MAX_CASK_R) THEN
     DO i = 1, n_trcr
        ydata(i) = ti_gp(i)%tp%meta%cask_r(imeta)
     END DO
  ELSE
     !     CALL info_bi('Meta data name '//TRIM(ymeta)//' unknown in MESSy' &
     !                , 'trcr_meta_get_r8v_byname' )
     !     CALL info_bi('returning zero', 'trcr_meta_get_r8v_byname'  )
     ydata(:) = 0.0_wp
  ENDIF
#endif
END SUBROUTINE trcr_meta_get_r8v_byname

SUBROUTINE trcr_meta_get_strv_byname(ierr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  CHARACTER(LEN=*)       , INTENT(OUT) :: ydata(:)
  INTEGER                              :: imeta
  INTEGER                              :: i
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_get(mtrcr,ymeta,ydata)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
  imeta = func_imeta(ymeta)
  IF (imeta > 0 .AND. imeta <= MAX_CASK_S) THEN
     DO i = 1, n_trcr
        ydata(i) = TRIM(ti_gp(i)%tp%meta%cask_s(imeta))
     END DO
  ELSE
     SELECT CASE(imeta)
     CASE(-1, -4, -5) ! name
        DO i=1,n_trcr
           ydata(i) = TRIM(ti_gp(i)%tp%ident%fullname)
        END DO
     CASE (-2)
        DO i=1,n_trcr
           ydata(i) = TRIM(ti_gp(i)%tp%ident%submodel)
        END DO
     CASE (-3)
        DO i=1,n_trcr
           ydata(i) = TRIM(ti_gp(i)%tp%ident%unit)
        END DO
     CASE DEFAULT
        ! CALL info_bi('Meta data name '//TRIM(ymeta)//' unknown in MESSy' &
        !            , 'trcr_meta_get_strv_byname' )
        ! CALL info_bi('returning empty string', 'trcr_meta_get_strv_byname'  )
        ydata(:) = ' '
     END SELECT
  ENDIF
#endif
END SUBROUTINE trcr_meta_get_strv_byname

SUBROUTINE trcr_meta_get_logv_byname(ierr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*), INTENT(IN) :: ymeta
  LOGICAL, INTENT(OUT) :: ydata(:)
  INTEGER :: imeta
  INTEGER :: i
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_get(mtrcr,ymeta,ydata)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
  imeta = func_imeta(ymeta)
  IF (imeta > 0 .AND. imeta <= MAX_CASK_I) THEN
     DO i = 1, n_trcr
        IF (ti_gp(i)%tp%meta%cask_i(imeta) == OFF) THEN
           ydata(i) = .FALSE.
        ELSE
           ydata(i) = .TRUE.
        ENDIF
     END DO
  ELSE
     !     CALL info_bi('Meta data name '//TRIM(ymeta)//' unknown in MESSy' &
     !          , 'trcr_meta_get_logv_byname' )
     !     CALL info_bi('returning .FALSE.', 'trcr_meta_get_logv_byname'  )
     ydata(:) = .FALSE.
  ENDIF
#endif
END SUBROUTINE trcr_meta_get_logv_byname

!==============================================================================

!==============================================================================
!+ Module procedures to set data into metadata storage
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!
! Description:
!   The following procedures overload the procedure trcr_meta_set which is
!   called upon a request to store a data item in the metadata storage. For
!   each data type as well as for scalar/vector data a separate procedure
!   has to be created. The procedures all have exactly the same structure with
!   the exception of dimension and type of the ydata argument as well as
!   some local type-dependent parameters. The procedures rely heavily on the
!   trcr_meta_set_check and trcr_meta_set_store functions which do most of the work.
!
! Subroutines:
!
!  trcr_meta_set_int_byindind     trcr_meta_set_r4_byindind
!  trcr_meta_set_r8_byindind      trcr_meta_set_str_byindind
!  trcr_meta_set_log_byindind     trcr_meta_set_ptr2_byindind
!  trcr_meta_set_ptr3_byindind    trcr_meta_set_ptr4_byindind
!  trcr_meta_set_int_byindname    trcr_meta_set_r4_byindname
!  trcr_meta_set_r8_byindname     trcr_meta_set_str_byindname
!  trcr_meta_set_log_byindname    trcr_meta_set_ptr2_byindname
!  trcr_meta_set_ptr3_byindname   trcr_meta_set_ptr4_byindname
!  trcr_meta_set_int_bynameind    trcr_meta_set_r4_bynameind
!  trcr_meta_set_r8_bynameind     trcr_meta_set_str_bynameind
!  trcr_meta_set_log_bynameind    trcr_meta_set_ptr2_bynameind
!  trcr_meta_set_ptr3_bynameind   trcr_meta_set_ptr4_bynameind
!  trcr_meta_set_int_bynamename   trcr_meta_set_r4_bynamename
!  trcr_meta_set_r8_bynamename    trcr_meta_set_str_bynamename
!  trcr_meta_set_log_bynamename   trcr_meta_set_ptr2_bynamename
!  trcr_meta_set_ptr3_bynamename  trcr_meta_set_ptr4_bynamename
!  trcr_meta_set_intv_byindex     trcr_meta_set_r4v_byindex
!  trcr_meta_set_r8v_byindex      trcr_meta_set_strv_byindex
!  trcr_meta_set_logv_byindex
!  trcr_meta_set_intv_byname      trcr_meta_set_r4v_byname
!  trcr_meta_set_r8v_byname       trcr_meta_set_strv_byname
!  trcr_meta_set_logv_byname
!
! Error numbers: 1200-1299
!
! Method:
!
!==============================================================================

SUBROUTINE trcr_meta_set_int_byindind(ierr,itrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER, INTENT(IN) :: itrcr
  INTEGER, INTENT(IN) :: imeta
  INTEGER, INTENT(IN) :: ydata
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    IF (imeta < 0 .OR. imeta > MAX_CASK_I) THEN
       CALL error_bi(' meta data index not defined in MESSy' &
            , 'trcr_meta_set_int_byindind' )
    ENDIF
    CALL set_tracer(ierr, GPTRSTR, itrcr, imeta, ydata)
    CALL tracer_halt('trcr_meta_set_int_byindind', ierr)
#endif
  ENDIF
END SUBROUTINE trcr_meta_set_int_byindind

SUBROUTINE trcr_meta_set_r4_byindind(ierr,itrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  INTEGER                , INTENT(IN)  :: imeta
  REAL(KIND=sp)          , INTENT(IN)  :: ydata
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    IF (imeta < 0 .OR. imeta > MAX_CASK_R) THEN
       CALL error_bi(' meta data index not defined in MESSy' &
            , 'trcr_meta_set_r4_byindind' )
    ENDIF
    CALL set_tracer(ierr, GPTRSTR, itrcr, imeta, REAL(ydata,wp))
    CALL tracer_halt('trcr_meta_set_r4_byindind', ierr)
#endif
  ENDIF
END SUBROUTINE trcr_meta_set_r4_byindind

SUBROUTINE trcr_meta_set_r8_byindind(ierr,itrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  INTEGER                , INTENT(IN)  :: imeta
  REAL(KIND=dp)          , INTENT(IN)  :: ydata
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    IF (imeta < 0 .OR. imeta > MAX_CASK_R) THEN
       CALL error_bi(' meta data index not defined in MESSy' &
            , 'trcr_meta_set_r8_byindind' )
    ENDIF
    CALL set_tracer(ierr, GPTRSTR, itrcr, imeta, ydata)
    CALL tracer_halt('trcr_meta_set_r8_byindind', ierr)
#endif
  ENDIF
END SUBROUTINE trcr_meta_set_r8_byindind

SUBROUTINE trcr_meta_set_str_byindind(ierr,itrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER, INTENT(IN) :: itrcr
  INTEGER, INTENT(IN) :: imeta
  CHARACTER(LEN=*), INTENT(IN) :: ydata
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    IF (imeta < 0 .OR. imeta > MAX_CASK_S) THEN
       CALL error_bi(' meta data index not defined in MESSy' &
            , 'trcr_meta_set_str_byindind' )
    ENDIF
    CALL set_tracer(ierr, GPTRSTR, itrcr, imeta, ydata)
    CALL tracer_halt('trcr_meta_set_str_byindind', ierr)
#endif
  ENDIF
END SUBROUTINE trcr_meta_set_str_byindind

SUBROUTINE trcr_meta_set_log_byindind(ierr,itrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER, INTENT(IN) :: itrcr
  INTEGER, INTENT(IN) :: imeta
  LOGICAL, INTENT(IN) :: ydata
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    IF (imeta < 0 .OR. imeta > MAX_CASK_I) THEN
       CALL error_bi(' meta data index not defined in MESSy' &
            , 'trcr_meta_set_log_byindind' )
    ENDIF
    IF (ydata) THEN
       CALL set_tracer(ierr, GPTRSTR, itrcr, imeta, ON)
       CALL tracer_halt('trcr_meta_set_log_byindind', ierr)
    ELSE
       CALL set_tracer(ierr, GPTRSTR, itrcr, imeta, OFF)
       CALL tracer_halt('trcr_meta_set_log_byindind', ierr)
    END IF
#endif
  ENDIF
END SUBROUTINE trcr_meta_set_log_byindind

SUBROUTINE trcr_meta_set_ptr2_byindind(ierr,itrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER, INTENT(IN) :: itrcr
  INTEGER, INTENT(IN) :: imeta
  REAL(KIND=wp),     POINTER, INTENT(IN) :: ydata(:,:)
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    ! CALL info_bi('Due to different concepts: No 2D-meta data in MESSy', '' )
    ierr = 5000_iintegers
#endif
  ENDIF
END SUBROUTINE trcr_meta_set_ptr2_byindind

SUBROUTINE trcr_meta_set_ptr3_byindind(ierr,itrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER, INTENT(IN) :: itrcr
  INTEGER, INTENT(IN) :: imeta
  REAL(KIND=wp),     POINTER, INTENT(IN) :: ydata(:,:,:)
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    ! CALL info_bi('Due to different concepts: No 3D-meta data in MESSy','' )
    ierr = 5000_iintegers
#endif
  ENDIF
END SUBROUTINE trcr_meta_set_ptr3_byindind

SUBROUTINE trcr_meta_set_ptr4_byindind(ierr,itrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER, INTENT(IN) :: itrcr
  INTEGER, INTENT(IN) :: imeta
  REAL(KIND=wp),     POINTER, INTENT(IN) :: ydata(:,:,:,:)
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    !CALL info_bi('Due to different concepts: No 4D-meta data in MESSy', '' )
    ierr = 5000_iintegers
#endif
  ENDIF
END SUBROUTINE trcr_meta_set_ptr4_byindind

SUBROUTINE trcr_meta_set_int_byindname(ierr,itrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER, INTENT(IN) :: itrcr
  CHARACTER(LEN=*), INTENT(IN) :: ymeta
  INTEGER, INTENT(IN) :: ydata
  INTEGER             :: imeta
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    imeta = func_imeta(ymeta)
    IF (imeta > 0 .AND. imeta <= MAX_CASK_I) THEN
       CALL set_tracer(ierr, GPTRSTR, itrcr, imeta, ydata)
       CALL tracer_halt('trcr_meta_set_int_byindname', ierr)
    ELSE
       ! CALL info_bi('Meta data name '//TRIM(ymeta)//' unknown in MESSy' &
       !            , 'trcr_meta_set_int_byindname' )
       !ierr = 5000_iintegers
    ENDIF
#endif
  ENDIF
END SUBROUTINE trcr_meta_set_int_byindname

SUBROUTINE trcr_meta_set_r4_byindname(ierr,itrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  REAL(KIND=sp)          , INTENT(IN)  :: ydata
  INTEGER                              :: imeta
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    imeta = func_imeta(ymeta)
    IF (imeta > 0 .AND. imeta <= MAX_CASK_R) THEN
       CALL set_tracer(ierr, GPTRSTR, itrcr, imeta, REAL(ydata,wp))
       CALL tracer_halt('trcr_meta_set_r4_byindname', ierr)
    ELSE
       ! CALL info_bi('Meta data name '//TRIM(ymeta)//' unknown in MESSy' &
       !            , 'trcr_meta_set_r4_byindname' )
       !ierr = 5000_iintegers
    ENDIF
#endif
  ENDIF
END SUBROUTINE trcr_meta_set_r4_byindname

SUBROUTINE trcr_meta_set_r8_byindname(ierr,itrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  REAL(KIND=dp)          , INTENT(IN)  :: ydata
  INTEGER                              :: imeta
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    imeta = func_imeta(ymeta)
    IF (imeta > 0 .AND. imeta <= MAX_CASK_R) THEN
       CALL set_tracer(ierr, GPTRSTR, itrcr, imeta, ydata)
       CALL tracer_halt('trcr_meta_set_r8_byindname', ierr)
    ELSE
       ! CALL info_bi('Meta data name '//TRIM(ymeta)//' unknown in MESSy' &
       !            , 'trcr_meta_set_r8_byindname' )
       !ierr = 5000_iintegers
    ENDIF
#endif
  ENDIF
END SUBROUTINE trcr_meta_set_r8_byindname

SUBROUTINE trcr_meta_set_str_byindname(ierr,itrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  CHARACTER(LEN=*)       , INTENT(IN)  :: ydata
  INTEGER                              :: imeta
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    imeta = func_imeta(ymeta)
    IF (imeta > 0 .AND. imeta <= MAX_CASK_S) THEN
       CALL set_tracer(ierr, GPTRSTR, itrcr, imeta, ydata)
       CALL tracer_halt('trcr_meta_set_str_byindname', ierr)
    ELSE
       ! CALL info_bi('Meta data name '//TRIM(ymeta)//' unknown in MESSy' &
       !            , 'trcr_meta_set_str_byindname' )
       !ierr = 5000_iintegers
    ENDIF
#endif
  ENDIF
END SUBROUTINE trcr_meta_set_str_byindname

SUBROUTINE trcr_meta_set_log_byindname(ierr,itrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  LOGICAL                , INTENT(IN)  :: ydata
  INTEGER                              :: imeta
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    imeta = func_imeta(ymeta)
    IF (imeta > 0 .AND. imeta <= MAX_CASK_S) THEN
       IF (ydata) THEN
          CALL set_tracer(ierr, GPTRSTR, itrcr, imeta, ON)
          CALL tracer_halt('trcr_meta_set_int_bynameind', ierr)
       ELSE
          CALL set_tracer(ierr, GPTRSTR, itrcr, imeta, OFF)
          CALL tracer_halt('trcr_meta_set_log_byindname', ierr)
       ENDIF
    ELSE
       ! CALL info_bi('Meta data name '//TRIM(ymeta)//' unknown in MESSy' &
       !            , 'trcr_meta_set_str_byindname' )
       !ierr = 5000_iintegers
    ENDIF
#endif
  ENDIF
END SUBROUTINE trcr_meta_set_log_byindname

SUBROUTINE trcr_meta_set_ptr2_byindname(ierr,itrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  REAL(KIND=wp), POINTER , INTENT(IN)  :: ydata(:,:)
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    !CALL info_bi('Due to different concepts: No 4D-meta data in MESSy', '' )
    ierr = 5000_iintegers
#endif
  ENDIF
END SUBROUTINE trcr_meta_set_ptr2_byindname

SUBROUTINE trcr_meta_set_ptr3_byindname(ierr,itrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  REAL(KIND=wp), POINTER , INTENT(IN)  :: ydata(:,:,:)
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    ! CALL info_bi('Due to different concepts: No 4D-meta data in MESSy' , '' )
    ierr = 5000_iintegers
#endif
  ENDIF
END SUBROUTINE trcr_meta_set_ptr3_byindname

SUBROUTINE trcr_meta_set_ptr4_byindname(ierr,itrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: itrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  REAL(KIND=wp), POINTER , INTENT(IN)  :: ydata(:,:,:,:)
  ierr = 0_iintegers
  CALL trcr_check_index(ierr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    ! CALL info_bi('Due to different concepts: No 4D-meta data in MESSy', '' )
    ierr = 5000_iintegers
#endif
  ENDIF
END SUBROUTINE trcr_meta_set_ptr4_byindname

SUBROUTINE trcr_meta_set_int_bynameind(ierr,ytrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  INTEGER                , INTENT(IN)  :: imeta
  INTEGER                , INTENT(IN)  :: ydata
  INTEGER                              :: itrcr
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    IF (imeta < 0 .OR. imeta > MAX_CASK_I) THEN
       CALL error_bi(' meta data index not defined in MESSy' &
            , 'trcr_meta_set_int_bynameind' )
    ENDIF
    CALL set_tracer(ierr, GPTRSTR, itrcr, imeta, ydata)
    CALL tracer_halt('trcr_meta_set_int_bynameind', ierr)
#endif
  END IF
END SUBROUTINE trcr_meta_set_int_bynameind

SUBROUTINE trcr_meta_set_r4_bynameind(ierr,ytrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  INTEGER                , INTENT(IN)  :: imeta
  REAL(KIND=sp)          , INTENT(IN)  :: ydata
  INTEGER                              :: itrcr
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    IF (imeta < 0 .OR. imeta > MAX_CASK_R) THEN
       CALL error_bi(' meta data index not defined in MESSy' &
            , 'trcr_meta_set_r4_bynameind' )
    ENDIF
    CALL set_tracer(ierr, GPTRSTR, itrcr, imeta, REAL(ydata,wp))
    CALL tracer_halt('trcr_meta_set_r4_bynameind', ierr)
#endif
  END IF
END SUBROUTINE trcr_meta_set_r4_bynameind

SUBROUTINE trcr_meta_set_r8_bynameind(ierr,ytrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  INTEGER                , INTENT(IN)  :: imeta
  REAL(KIND=dp)          , INTENT(IN)  :: ydata
  INTEGER                              :: itrcr
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    IF (imeta < 0 .OR. imeta > MAX_CASK_R) THEN
       CALL error_bi(' meta data index not defined in MESSy' &
            , 'trcr_meta_set_r8_bynameind' )
    ENDIF
    CALL set_tracer(ierr, GPTRSTR, itrcr, imeta, ydata)
    CALL tracer_halt('trcr_meta_set_r8_bynameind', ierr)
#endif
  END IF
END SUBROUTINE trcr_meta_set_r8_bynameind

SUBROUTINE trcr_meta_set_str_bynameind(ierr,ytrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  INTEGER                , INTENT(IN)  :: imeta
  CHARACTER(LEN=*)       , INTENT(IN)  :: ydata
  INTEGER                              :: itrcr
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    IF (imeta < 0 .OR. imeta > MAX_CASK_S) THEN
       CALL error_bi(' meta data index not defined in MESSy' &
            , 'trcr_meta_set_str_bynameind' )
    ENDIF
    CALL set_tracer(ierr, GPTRSTR, itrcr, imeta, ydata)
    CALL tracer_halt('trcr_meta_set_str_bynameind', ierr)
#endif
  END IF
END SUBROUTINE trcr_meta_set_str_bynameind

SUBROUTINE trcr_meta_set_log_bynameind(ierr,ytrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  INTEGER                , INTENT(IN)  :: imeta
  LOGICAL                , INTENT(IN)  :: ydata
  INTEGER                              :: itrcr
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    IF (imeta < 0 .OR. imeta > MAX_CASK_I) THEN
       CALL error_bi(' meta data index not defined in MESSy' &
            , 'trcr_meta_set_log_bynameind' )
    ENDIF
    IF (ydata) THEN
       CALL set_tracer(ierr, GPTRSTR, itrcr, imeta, ON)
       CALL tracer_halt('trcr_meta_set_log_bynameind', ierr)
    ELSE
       CALL set_tracer(ierr, GPTRSTR, itrcr, imeta, OFF)
       CALL tracer_halt('trcr_meta_set_log_bynameind', ierr)
    ENDIF
#endif
  END IF
END SUBROUTINE trcr_meta_set_log_bynameind

SUBROUTINE trcr_meta_set_ptr2_bynameind(ierr,ytrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  INTEGER                , INTENT(IN)  :: imeta
  REAL(KIND=wp), POINTER , INTENT(IN)  :: ydata(:,:)
  INTEGER                              :: itrcr
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    ! CALL info_bi('Due to different concepts: No 2D-meta data in MESSy', '' )
    ierr = 5000_iintegers
#endif
  END IF
END SUBROUTINE trcr_meta_set_ptr2_bynameind

SUBROUTINE trcr_meta_set_ptr3_bynameind(ierr,ytrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  INTEGER                , INTENT(IN)  :: imeta
  REAL(KIND=wp), POINTER , INTENT(IN)  :: ydata(:,:,:)
  INTEGER :: itrcr
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    ! CALL info_bi('Due to different concepts: No 3D-meta data in MESSy' , '' )
    ierr = 5000_iintegers
#endif
  END IF
END SUBROUTINE trcr_meta_set_ptr3_bynameind

SUBROUTINE trcr_meta_set_ptr4_bynameind(ierr,ytrcr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  INTEGER                , INTENT(IN)  :: imeta
  REAL(KIND=wp), POINTER , INTENT(IN)  :: ydata(:,:,:,:)
  INTEGER :: itrcr
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,imeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    ! CALL info_bi('Due to different concepts: No 4D-meta data in MESSy', '' )
    ierr = 5000_iintegers
#endif
  END IF
END SUBROUTINE trcr_meta_set_ptr4_bynameind

SUBROUTINE trcr_meta_set_int_bynamename(ierr,ytrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  INTEGER                , INTENT(IN)  :: ydata
  INTEGER                              :: itrcr
  INTEGER                              :: imeta
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    imeta = func_imeta(ymeta)
    IF (imeta > 0 .AND. imeta <= MAX_CASK_I) THEN
       CALL set_tracer(ierr, GPTRSTR, itrcr, imeta, ydata)
       CALL tracer_halt('trcr_meta_set_int_bynamename' , ierr)
    !ELSE
       ! CALL info_bi('Meta data name '//TRIM(ymeta)//'  unknown in MESSy' &
       !            , 'trcr_meta_set_int_bynamename' )
    ENDIF
#endif
  END IF
END SUBROUTINE trcr_meta_set_int_bynamename

SUBROUTINE trcr_meta_set_r4_bynamename(ierr,ytrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  REAL(KIND=sp)          , INTENT(IN)  :: ydata
  INTEGER                              :: itrcr
  INTEGER                              :: imeta
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    imeta = func_imeta(ymeta)
    IF (imeta > 0 .AND. imeta <= MAX_CASK_R) THEN
       CALL set_tracer(ierr, GPTRSTR, itrcr, imeta, REAL(ydata,wp))
       CALL tracer_halt('trcr_meta_set_r4_bynamename', ierr)
    ! ELSE
       !  CALL info_bi('Meta data name '//TRIM(ymeta)//'  unknown in MESSy' &
       !             , 'trcr_meta_set_r4_bynamename' )
    ENDIF
#endif
  END IF
END SUBROUTINE trcr_meta_set_r4_bynamename

SUBROUTINE trcr_meta_set_r8_bynamename(ierr,ytrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  REAL(KIND=dp)          , INTENT(IN)  :: ydata
  INTEGER                              :: itrcr
  INTEGER                              :: imeta
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    imeta = func_imeta(ymeta)
    IF (imeta > 0 .AND. imeta <= MAX_CASK_R) THEN
       CALL set_tracer(ierr, GPTRSTR, itrcr, imeta, ydata)
       CALL tracer_halt('trcr_meta_set_r8_bynamename', ierr)
    ! ELSE
       ! CALL info_bi('Meta data name '//TRIM(ymeta)//'  unknown in MESSy' &
       !            , 'trcr_meta_set_r8_bynamename' )
    ENDIF
#endif
  END IF
END SUBROUTINE trcr_meta_set_r8_bynamename

SUBROUTINE trcr_meta_set_str_bynamename(ierr,ytrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  CHARACTER(LEN=*)       , INTENT(IN)  :: ydata
  INTEGER                              :: itrcr
  INTEGER                              :: imeta
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    imeta = func_imeta(ymeta)
    IF (imeta > 0 .AND. imeta <= MAX_CASK_S) THEN
       CALL set_tracer(ierr, GPTRSTR, itrcr, imeta, ydata)
       CALL tracer_halt('trcr_meta_set_str_bynamename', ierr)
    ! ELSE
       ! CALL info_bi('Meta data name '//TRIM(ymeta)//'  unknown in MESSy' &
       !            , 'trcr_meta_set_str_bynamename' )
    ENDIF
#endif
  END IF
END SUBROUTINE trcr_meta_set_str_bynamename

SUBROUTINE trcr_meta_set_log_bynamename(ierr,ytrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  LOGICAL                , INTENT(IN)  :: ydata
  INTEGER                              :: itrcr
  INTEGER                              :: imeta
  ierr = 0_iintegers
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
#ifndef MESSY
    CALL metadata_set(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
    imeta = func_imeta(ymeta)
    IF (imeta > 0 .AND. imeta <= MAX_CASK_S) THEN
       IF (ydata) THEN
          CALL set_tracer(ierr, GPTRSTR, itrcr, imeta, ON)
          CALL tracer_halt('trcr_meta_set_log_bynamename', ierr)
       ELSE
          CALL set_tracer(ierr, GPTRSTR, itrcr, imeta, OFF)
          CALL tracer_halt('trcr_meta_set_log_bynamename', ierr)
       ENDIf
     ! ELSE
       !CALL info_bi('Meta data name '//TRIM(ymeta)//' unknown in MESSy: IGNORE'&
       !           , 'trcr_meta_set_str_bynamename' )
    ENDIF
#endif
  END IF
END SUBROUTINE trcr_meta_set_log_bynamename

SUBROUTINE trcr_meta_set_ptr2_bynamename(ierr,ytrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  REAL(KIND=wp), POINTER , INTENT(IN)  :: ydata(:,:)
  INTEGER :: itrcr
  ierr = 0_iintegers
#ifndef MESSY
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
    CALL metadata_set(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
  END IF
#else
  CALL error_bi('this subroutine should not be called in case of MESSy'  &
       , 'trcr_meta_set_ptr2_bynamename' )
#endif
END SUBROUTINE trcr_meta_set_ptr2_bynamename

SUBROUTINE trcr_meta_set_ptr3_bynamename(ierr,ytrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  REAL(KIND=wp), POINTER , INTENT(IN)  :: ydata(:,:,:)
  INTEGER :: itrcr
  ierr = 0_iintegers
#ifndef MESSY
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
    CALL metadata_set(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
  END IF
#else
  CALL error_bi('this subroutine should not be called in case of MESSy'  &
       , 'trcr_meta_set_ptr3_bynamename' )
#endif
END SUBROUTINE trcr_meta_set_ptr3_bynamename

SUBROUTINE trcr_meta_set_ptr4_bynamename(ierr,ytrcr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ytrcr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  REAL(KIND=wp), POINTER , INTENT(IN)  :: ydata(:,:,:,:)
  INTEGER :: itrcr
  ierr = 0_iintegers
#ifndef MESSY
  CALL trcr_get_index(ierr, ytrcr, itrcr)
  IF (ierr == 0_iintegers) THEN
    CALL metadata_set(mtrcr,itrcr,ymeta,ydata)
    IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
  END IF
#else
  CALL error_bi('this subroutine should not be called in case of MESSy'  &
       , 'trcr_meta_set_ptr4_bynamename' )
#endif
END SUBROUTINE trcr_meta_set_ptr4_bynamename

SUBROUTINE trcr_meta_set_intv_byindex(ierr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: imeta
  INTEGER                , INTENT(IN)  :: ydata(:)
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_set(mtrcr,imeta,ydata)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
  CALL error_bi('this subroutine should not be called in case of MESSy'  &
       , 'trcr_meta_set_intv_byindex' )
#endif
END SUBROUTINE trcr_meta_set_intv_byindex

SUBROUTINE trcr_meta_set_r4v_byindex(ierr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: imeta
  REAL(KIND=sp)          , INTENT(IN)  :: ydata(:)
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_set(mtrcr,imeta,ydata)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
  CALL error_bi('this subroutine should not be called in case of MESSy'  &
       , 'trcr_meta_set_r4v_byindex' )
#endif
END SUBROUTINE trcr_meta_set_r4v_byindex

SUBROUTINE trcr_meta_set_r8v_byindex(ierr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: imeta
  REAL(KIND=dp)          , INTENT(IN)  :: ydata(:)
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_set(mtrcr,imeta,ydata)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
  CALL error_bi('this subroutine should not be called in case of MESSy'  &
       , 'trcr_meta_set_r8v_byindex' )
#endif
END SUBROUTINE trcr_meta_set_r8v_byindex

SUBROUTINE trcr_meta_set_strv_byindex(ierr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER                , INTENT(IN)  :: imeta
  CHARACTER(LEN=*)       , INTENT(IN)  :: ydata(:)
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_set(mtrcr,imeta,ydata)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
  CALL error_bi('this subroutine should not be called in case of MESSy'  &
       , 'trcr_meta_set_strv_byindex' )
#endif
END SUBROUTINE trcr_meta_set_strv_byindex

SUBROUTINE trcr_meta_set_logv_byindex(ierr,imeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  INTEGER, INTENT(IN) :: imeta
  LOGICAL, INTENT(IN) :: ydata(:)
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_set(mtrcr,imeta,ydata)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
  CALL error_bi('this subroutine should not be called in case of MESSy'  &
       , 'trcr_meta_set_logv_byindex' )
#endif
END SUBROUTINE trcr_meta_set_logv_byindex

SUBROUTINE trcr_meta_set_intv_byname(ierr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  INTEGER                , INTENT(IN)  :: ydata(:)
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_set(mtrcr,ymeta,ydata)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
  CALL error_bi('this subroutine should not be called in case of MESSy'  &
       , 'trcr_meta_set_intv_byname' )
#endif
END SUBROUTINE trcr_meta_set_intv_byname

SUBROUTINE trcr_meta_set_r4v_byname(ierr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  REAL(KIND=sp)          , INTENT(IN)  :: ydata(:)
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_set(mtrcr,ymeta,ydata)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
  CALL error_bi('this subroutine should not be called in case of MESSy'  &
       , 'trcr_meta_set_r4v_byname' )
#endif
END SUBROUTINE trcr_meta_set_r4v_byname

SUBROUTINE trcr_meta_set_r8v_byname(ierr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  REAL(KIND=dp)          , INTENT(IN)  :: ydata(:)
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_set(mtrcr,ymeta,ydata)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
  CALL error_bi('this subroutine should not be called in case of MESSy'  &
       , 'trcr_meta_set_r8v_byname' )
#endif
END SUBROUTINE trcr_meta_set_r8v_byname

SUBROUTINE trcr_meta_set_strv_byname(ierr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  CHARACTER(LEN=*)       , INTENT(IN)  :: ydata(:)
  ierr = 0_iintegers
#ifndef MESSY
  CALL metadata_set(mtrcr,ymeta,ydata)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
  CALL error_bi('this subroutine should not be called in case of MESSy'  &
       , 'trcr_meta_set_strv_byname' )
#endif
END SUBROUTINE trcr_meta_set_strv_byname

SUBROUTINE trcr_meta_set_logv_byname(ierr,ymeta,ydata)
  IMPLICIT NONE
  INTEGER(KIND=iintegers), INTENT(OUT) :: ierr
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  LOGICAL                , INTENT(IN)  :: ydata(:)
#ifndef MESSY
  ierr = 0_iintegers
  CALL metadata_set(mtrcr,ymeta,ydata)
  IF (metadata_error() /= 0_iintegers) ierr = 5000_iintegers
#else
  CALL error_bi('this subroutine should not be called in case of MESSy'  &
       , 'trcr_meta_set_logv_byname' )
#endif
END SUBROUTINE trcr_meta_set_logv_byname

!==============================================================================
!==============================================================================
#ifdef MESSY
SUBROUTINE list_meta(ymeta, imetac, imetam)
  
  CHARACTER(LEN=*)       , INTENT(IN)  :: ymeta
  INTEGER(KIND=iintegers), INTENT(OUT) :: imetac
  INTEGER(KIND=iintegers), INTENT(OUT) :: imetam
  ! LOCAL
  INTEGER :: i
  
  DO i = 1, NMAX_META
      IF (TRIM(ymeta) == TRIM(meta_list(i)%cosmo_str)) THEN
         imetac = meta_list(i)%cosmo_idx
         imetam = meta_list(i)%messy_idx
         RETURN
      END IF
  END DO
  
  imetac = -99
  imetam = -99
  
  RETURN
END SUBROUTINE list_meta
INTEGER FUNCTION func_imeta(ymeta)

  CHARACTER(LEN=*) :: ymeta
  ! LOCAL
  INTEGER :: i
  DO i = 1, NMAX_META
      IF (TRIM(ymeta) == TRIM(meta_list(i)%cosmo_str)) THEN
         func_imeta = meta_list(i)%messy_idx
         RETURN
      END IF
  END DO

  func_imeta = -99

  RETURN
END FUNCTION func_imeta
#endif

!==============================================================================
!==============================================================================
!+ Module procedure "trcr_acc_update_device" in "src_tracer"
!------------------------------------------------------------------------------

SUBROUTINE trcr_acc_update_device
!------------------------------------------------------------------------------
!
! Description:
!
! This subroutine copies the host values of tracer structures to the device
!
!
! Method:
!
! Using OpenACC update directives
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine trcr_acc_update_device
!------------------------------------------------------------------------------

 !$acc update device(trcr_data,trcr_data_bd,trcr_data_tens)

END SUBROUTINE trcr_acc_update_device

!==============================================================================
!+ Module procedure "trcr_acc_update_device" in "src_tracer"
!------------------------------------------------------------------------------

SUBROUTINE trcr_acc_update_host
!------------------------------------------------------------------------------
!
! Description:
!
! This subroutine copies the host values of tracer structures to the host
!
!
! Method:
!
! Using OpenACC update directives
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine trcr_acc_update_host
!------------------------------------------------------------------------------

 !$acc update host(trcr_data,trcr_data_bd,trcr_data_tens)

END SUBROUTINE trcr_acc_update_host

!==============================================================================
!==============================================================================
!+ Module procedure "trcr_errorstr" in "src_tracer" to return an error message
!------------------------------------------------------------------------------

FUNCTION trcr_errorstr ( ierr )
!------------------------------------------------------------------------------
!
! Description:
!
! This function returns an error message in a string which corresponds to the
! error number supplied as an argument.
!
!
! Method:
!
!
!------------------------------------------------------------------------------

 IMPLICIT NONE

! Subroutine arguments:
! --------------------
 INTEGER(KIND=iintegers), INTENT(IN) :: ierr  ! error status
 CHARACTER(LEN=ilen_en)              :: trcr_errorstr ! number of tracers

! Local parameters:
! ----------------

 INTEGER :: itmp
 CHARACTER(LEN=256) :: zstr
 CHARACTER(LEN=16) :: zstr2

!------------------------------------------------------------------------------
! Begin Subroutine trcr_errorstr
!------------------------------------------------------------------------------

  SELECT CASE(ierr)

  ! no error
  CASE(0)
    zstr = 'NO ERROR'

  ! tracer not found error
  CASE(T_ERR_NOTFOUND)
    zstr = 'NO MATCHING TRACER FOUND'

  ! trcr_init errors         (4000-4049)
  CASE(4000)
    zstr = 'TRCR_INIT CALLED OUT OF SEQUENCE'

  ! trcr_new errors          (4050-4099)
  CASE(4050)
    zstr = 'TRCR_NEW CALLED OUT OF SEQUENCE'
  CASE(4055)
    zstr = 'TOTAL NUMBER OF TRACERS REACHED. INCREASE ntrcr_max IN data_tracer.f90'
  CASE(4060)
    zstr = 'TRACER NAME STRING TOO LONG'
  CASE(4062)
    zstr = 'PARENT NAME STRING TOO LONG' 
  CASE(4064)
    zstr = 'UNITS STRING TOO LONG' 
  CASE(4066)
    zstr = 'STANDARD NAME STRING TOO LONG' 
  CASE(4067)
    zstr = 'LONG NAME STRING TOO LONG'
  CASE(4068)
    zstr = 'LEVEL TYPE STRING TOO LONG' 
  CASE(4070)
    zstr = 'TRACER NAME NOT UNIQUE'
  CASE(4080)
    zstr = 'SPECIFIED GRIB NUMBER NOT IN RANGE 0..255'  
  CASE(4081)
    zstr = 'UNKNOWN ADVECTION METHOD (itype_adv) CHOSEN'  
  CASE(4082)
    zstr = 'UNKNOWN DIFFUSION METHOD (itype_diff) CHOSEN'
  CASE(4083)
    zstr = 'UNKNOWN TURBULENT MIXING METHOD (itype_turbmix) CHOSEN'
  CASE(4084)
    zstr = 'UNKNOWN CONVECTION PASSIVE TRANSPORT METHOD (itype_passconv) CHOSEN'
  CASE(4085)
    zstr = 'UNKNOWN INITIALIZATION METHOD (itype_ini) CHOSEN'
  CASE(4086)
    zstr = 'UNKNOWN LATERAL BOUNDARY CONDITIONS METHOD (itype_lbc) CHOSEN'
  CASE(4087)
    zstr = 'UNKNOWN BOTTOM BOUNDARY CONDITIONS METHOD (itype_bbc) CHOSEN'
  CASE(4088)
    zstr = 'UNKNOWN RELAXATION METHOD (itype_relax) CHOSEN'  
  CASE(4089)
    zstr = 'UNKNOWN DAMPING METHOD (itype_damp) CHOSEN'  
  CASE(4090)
    zstr = 'UNKNOWN CLIPPING METHOD (itype_clip) CHOSEN'    

  ! trcr_alloc errors        (4100-4149)
  CASE(4100)
    zstr = 'TRCR_ALLOC CALLED OUT OF SEQUENCE'
  CASE(4105)
    zstr = 'MEMORY ALLOCATION FAILED'
  
  ! trcr_setup_vartab errors (4150-4199)
  CASE(4150)
    zstr = 'TRCR_SETUP_VARTAB CALLED OUT OF SEQUENCE'
  CASE(4155)
    zstr = 'SPECIFIED GRIB TABLE NUMBER NOT IN RANGE 1..num_gribtabs'
  CASE(4157)
    zstr = 'NAME CONFLICT BETWEEN TRACER AND VARIABLE TABLE'
  CASE(4159)
    zstr = 'RAN OUT OF LOCATIONS IN VARIABLE TABLE. INCREASE max_gribrep IN data_io.f90'
  
  ! trcr_print errors        (4200-4249)
  CASE(4200)
    zstr = 'TRCR_PRINT CALLED OUT OF SEQUENCE'  

  ! trcr_get errors          (4250-4289)
  CASE(4250)
    zstr = 'TRCR_GET CALLED OUT OF SEQUENCE'
  CASE(4253)
    zstr = 'TIMELEVEL SPECIFIED BUT NO TIME DEPENDENT POINTER REQUIRED'
  CASE(4255)
    zstr = 'ILLEGAL TIMELEVEL SPECIFIED TO TRCR_GET'
  CASE(4262)
    zstr = 'DATA POINTER REQUESTED FOR TRACER WHICH DOES NOT HAVE ONE'
  CASE(4265)
    zstr = 'BOUNDARY DATA POINTER REQUESTED FOR TRACER WHICH DOES NOT HAVE ONE'
  CASE(4270)
    zstr = 'TENDENCY DATA POINTER REQUESTED FOR TRACER WHICH DOES NOT HAVE ONE'

  ! trcr_swap errors      (4290-4299)
  CASE(4290)
    zstr = 'TRCR_SWAP CALLED OUT OF SEQUENCE'

  ! trcr_cleanup errors      (4300-4349)
  CASE(4300)
    zstr = 'TRCR_CLEANUP CALLED OUT OF SEQUENCE'
  CASE(4305)
    zstr = 'MEMORY DEALLOCATION FAILED'

  ! trcr_get_index errors      (4350-4359)
  CASE(4350)
    zstr = 'TRCR_GET_IDX CALLED OUT OF SEQUENCE'

  ! metadata errors
  CASE(5000)
    itmp = metadata_error(zstr)
! FUO: check this, should probably contain a newline
    WRITE(zstr2,'(i4.4)') itmp
    zstr = 'ERROR METADATA ' // TRIM(zstr2) // ': ' // TRIM(zstr)

  END SELECT

  ! add index or name of the tracer to the error message if meaningful
  IF ((ierr == T_ERR_NOTFOUND)) THEN
    WRITE(yzerrind, '(i3)') izerrind
    zstr = TRIM(zstr) //' (TRCR= '//TRIM(yzerrind)//')'
  ENDIF
  IF ((ierr >= 4050_iintegers .AND. ierr <= 4099_iintegers)) THEN
    zstr = TRIM(zstr) //' (TRCR= '//TRIM(yzerrind)//')'
  ENDIF 
  IF ((ierr >= 4250_iintegers .AND. ierr <= 4299_iintegers)) THEN
    WRITE(yzerrind, '(i3)') izerrind
    zstr = TRIM(zstr) //' (TRCR= '//TRIM(yzerrind)//')'
  ENDIF 

  ! construct error string:
  WRITE(zstr2,'(i4.4)') ierr
  trcr_errorstr = 'ERROR TRCR ' // TRIM(zstr2) // ': ' // TRIM(zstr)

!------------------------------------------------------------------------------
! End of module function trcr_errorstr
!------------------------------------------------------------------------------

END FUNCTION trcr_errorstr
!==============================================================================

!==============================================================================

END MODULE src_tracer
