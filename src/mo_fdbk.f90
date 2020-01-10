!
!+ 3DVAR/COSMO feedback file specific derived type 't_fdbk' and routines
!
! $Id: mo_fdbk.f90,v 5.0.1.1 2014-02-21 09:36:32 for0adm Exp $
!-------------------------------------------------------------------------------
!
MODULE mo_fdbk
!
!-------------------------------------------------------------------------------
! Description:
!   This module defines the derived type 't_fdbk'.  In addition to
!   the information on dimensions and variables (within component 'nc'
!   of derived type 't_netcdf_file') the type holds feedback file
!   specific information (number of data items actually written, global
!   attributes). Subroutine 'setup_fdbk' sets up the information
!   within this data type, 'create_fdbk' writes a NetCDF file based
!   on the meta-data stored in a variable of this type, 'close_fdbk'
!   closes the NetCDF-file and 'cleanup_fdbk' deallocates pointer
!   components of 't_fdbk'. All the routines are based on those
!   defined in 'mo_t_netcdf_file' but extends them by the feedback file
!   specific parts.
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
! V4_28        2013/07/12 Christoph Schraff / Andreas Rhodin
!  Update to 3DVAR-V1_26:
!  3DVAR-V1_14: New routines write_veri, write_fdbk_var; add error handling.
!  3DVAR-V1_15: Increment feedback file version number to 1.01;
!               make error processing COSMO conform (use 'model_abort') from
!               module 'environment.f90');
!               fix subr. write_veri: pass real not real(wp) to nf_put_vara_real;
!               new specific subr. write_fdbk_char: write char. string variable.
!  3DVAR-V1_16: Make 'phase' a mandatory parameter for SCATTerometer (=fov).
!  3DVAR-V1_17: New optional parameters: add_history (..,write),
!               create_fdbk (..,create);
!               provide single + double precision version of subr. 'write_veri'.
!  3DVAR-V1_18: Extended comments and documentation.
!  3DVAR-V1_19: Add variables for RADAR operator;
!               revised IASI diagnostics in feedback file.
!  3DVAR-V1_20: 'spec_index' as NF_INT (not NF_SHORT); fix minor memory leaks.
!  3DVAR-V1_22: RADAR observation type specific changes.  Write l2c as float.
!               Make 'plevel' a mandatory entry.  Define "varno_back" as SHORT.
!               Add observation specific parameters (obs_par, t_datum).
!  3DVAR-V1_23: New variable 'ct_nwc' (Cloud Type according to NWC-SAF). Write
!               localisation length scales (v_loc/h_loc) to LETKF feedback files
!  3DVAR-V1_26: E.g.: Adjust explanation of veri_domain_size for ICON.
! V5_1         2014-11-28 Andreas Rhodin, Ulrich Schaettler, Christoph Schraff
!  - Added KIND parameters for standard real variables
!  - Replaced mo_kind by kind_parameters (US)
!  Update to 3DVAR-V1_31 (CS):
!  3DVAR-V1_27: dbkz: change SHORT to INT; surftype: changed from BYTE to SHORT.
!  3DVAR-V1_28: optional argument n to get_veri to read limited number of obs
!               from fof-file.
!  3DVAR-V1_31: new variable ch_nwc: Cloud Height according to NWC-SAF;
!               adaptions for ensemble entries.
! V5_3         2015-10-09 Christoph Schraff
!  Update to 3DVAR-V1_42:
!  3DVAR-V1_35: changes for GPSGB (SPD: body var. azimuth, level type elevation)
!  3DVAR-V1_37: change index_x from short to int (for ICON)
!  3DVAR-V1_42: define variables dlat, dlon for TEMP;
!               declare 'z_station' as NF_INT; add 'veri_operator_flag'.
! V5_4a        2016-05-10 Ulrich Blahak
!  Added the RADAR optional flag to 'dlon' and 'dlat' to add_var (nc, ...)
!   to make dlat and dlon available for radar file bodies (UB)
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
! Andreas Rhodin  DWD  2008  original code
!                      2008  minor modifications
!------------------------------------------------------------------------------

!-------------
! Modules used
!-------------
use kind_parameters,  only: sp, dp, wp            ! real kind parameters
use environment,      only: model_abort           ! abort in case of error
use mo_t_netcdf_file, only: t_netcdf_file          ,&
                            add_dim                ,&
                            add_var                ,&
                            create_netcdf_file     ,&
                            open_netcdf_file_read  ,&
                            open_netcdf_file_write ,&
                            close_netcdf_file      ,&
                            destruct_netcdf_file
use mo_netcdf_param,  only: NF_INT, NF_SHORT, NF_BYTE, NF_CHAR, NF_FLOAT,     &
                            NF_DOUBLE, NF_FILL_SHORT,                         &
                            nf_put_att_text, nf_put_att_int, nf_put_att_real, &
                            nf_get_att_text, nf_get_att_int, nf_get_att_real, &
                            NF_GLOBAL, nf_enddef, nf_redef, nf_put_vara_text, &
                            nf_put_vara_int, nf_put_vara_real, nf_inq_attlen, &
                            nf_inq_dimid, nf_inq_dimlen, nf_inq_varid,        &
                            nf_get_var_text, nf_get_var_int, nf_get_var_real, &
                            nf_get_vara_real, NF_NOERR, nf_strerror,          &
                            nf_put_var_int, nf_put_var_real, nf_put_var_text, &
                            NF_EINDEFINE
use mo_fdbk_tables,   only: init_fdbk_tables ,&
                            status           ,&
                            flags            ,&
                            obstype          ,&
                            codetype         ,&
                            varno            ,&
                            runtype          ,&
                            runclass         ,&
                            satsens          ,&
                            surftype         ,&
                            surf_char        ,&
                            ct_nwc           ,&
                            flg_1dvar        ,&
                            flg_cld          ,&
                            level_sig        ,&
                            retrtype         ,&
                            ensmem           ,&
                            CT_NOPR          ,&
                            VE_MEMBER
implicit none

!----------------
! Public entities
!----------------

  !----------------------------------------
  ! set up the meta data of a feedback file
  ! in the generic data type t_netcdf_file
  !----------------------------------------

private
public :: t_fdbk                  ! feedback file data type
public :: t_fdbk_meta             ! component of t_fdbk
public :: setup_fdbk              ! set up t_netcdf_file for feedback file
public :: create_fdbk             ! create new feedback file
public :: open_fdbk_read          ! open feedback file for read access
public :: open_fdbk_write         ! open feedback file for write access
public :: add_history             ! add history entry to global attributes
public :: add_verification        ! add verification data to file
public :: write_global_attributes ! write global attributes
public :: read_meta               ! read  attributes + verification meta data
public :: print_fdbk              ! print attributes + verification meta data
public :: close_fdbk              ! close a feedback file
public :: cleanup_fdbk            ! deallocate components of t_fdbk
public :: get_veri_index          ! return indices of verification runs
public :: get_veri                ! read verification entry
public :: get_varid               ! get NetCDF varid
public :: get_fillvalue           ! get NetCDF fillvalue
public :: write_fdbk_var          ! write variable to feedback file
public :: write_veri              ! write verification entry

interface write_fdbk_var
  module procedure write_fdbk_int
  module procedure write_fdbk_real
  module procedure write_fdbk_char
end interface write_fdbk_var

interface write_veri
  module procedure write_veri_sp
  module procedure write_veri_dp
end interface write_veri

interface get_veri
  module procedure get_veri
  module procedure get_veri_n
end interface get_veri

!---------------------------------------------
! len of global attributes (character strings)
!---------------------------------------------
integer, parameter ::  TLEN = 28 ! len of title       attribute
integer, parameter ::  HLEN = 80 ! len of history     attribute
integer, parameter ::  ILEN = 24 ! len of institution attribute
integer, parameter ::  SLEN = 16 ! len of source      attribute

!-----------------------------------------------
! len of meta data variables (character strings)
!-----------------------------------------------
integer, parameter ::  MLEN = 10 ! len of 'veri_model'
integer, parameter :: IDLEN = 12 ! len of 'veri_initial_date'
integer, parameter ::  DLEN = 64 ! len of 'veri_description'

character(len=5), parameter :: FEEDBACK_VERSION = ' 1.01'

type t_fdbk_meta
  character(len=MLEN) :: model          = ' '    ! model used for verification
  integer             :: run_type       = -1     ! type of model run
  integer             :: run_class      = -1     ! class of model run
  character(len=IDLEN):: initial_date   = ' '    ! start of verification period
  integer             :: forecast_time  =  0     ! forecast time at verification_ref_time
  real(sp)            :: resolution (2) =  0._sp ! model resolution (x,y)
  integer             :: domain_size(3) =  0     ! domain size      (x,y,z)
  character(len=DLEN) :: description    = ' '    ! detailed description
  integer             :: ens_member     =  0     ! ensemble member number / special values
  integer             :: exp_id         = -1     ! experiment Id
  integer             :: operator_flag  =  0     ! observation operator flag
end type t_fdbk_meta

type t_fdbk                                  ! specific feedback file data
  !--------------------------------
  ! global feedback file attributes
  !--------------------------------
  integer              :: n_hdr         = 0 ! number of header records    used
  integer              :: n_body        = 0 ! number of body   records    used
  integer              :: n_radar       = 0 ! number of radar  records    used
  integer              :: n_veri        = 0 ! number of verification runs used
  character(len=TLEN)  :: title         ='' ! (CF) title
  character(len=5)     :: version       ='' ! (CF) version.subversion
  character(len=ILEN)  :: institution   ='' ! (CF) institution
  character(len=HLEN)  &                    ! (CF) history
              ,pointer :: history(:) =>NULL()
  character(len=SLEN)  :: source        ='' ! (CF) source
  integer              :: refdate       = 0 ! reference date 'yyyymmdd'
  integer              :: reftime       = 0 ! reference time 'hhmm'
  integer              :: start         = 0 ! verification start (minutes)
  integer              :: end           = 0 ! verification end   (minutes)
  real(sp)             :: resolution (2)= 0 ! resolution : lat,lon  (degree)
  integer              :: domain     (3)= 0 ! domain size: x,y,z
  real(sp)             :: pole       (2)= 0 ! pole coordinates : lat,lon
  real(sp)             :: lower_left (2)= 0 ! lower left  corner of domain
  real(sp)             :: upper_right(2)= 0 ! upper right corner of domain
  !------------------------------
  ! generic NetCDF data structure
  !------------------------------
  type (t_netcdf_file) :: nc                ! generic NetCDF data structure
  !---------------------------
  ! verification run meta data
  !---------------------------
  type(t_fdbk_meta) ,pointer :: veri(:) =>NULL() ! verification run meta data
end type t_fdbk

contains
!------------------------------------------------------------------------------
  subroutine setup_fdbk (nc, latex)
  type (t_netcdf_file) ,intent(out)          :: nc    ! NetCDF data type
  logical              ,intent(in) ,optional :: latex ! write LaTEX tables
  !----------------------------------------
  ! set up the meta data of a feedback file
  ! in the generic data type t_netcdf_file
  !----------------------------------------
    integer ,parameter :: NDIMS      =  9
    integer ,parameter :: NHEAD      = 51
    integer ,parameter :: NBODY      = 22
    integer ,parameter :: NRADAR     =  6
    integer ,parameter :: NMETA      = 11
    integer ,parameter :: NVERY      =  1
    integer ,parameter :: NVARS      = NHEAD + NBODY + NVERY + NRADAR + NMETA

    integer :: p_h, p_b, p_v, p_2, p_3, p_m, p_i, p_d, p_r

    !------------------
    ! set up the tables
    !------------------
    call init_fdbk_tables (latex)
    !----------------------------------
    ! allocate dimensions and variables
    !----------------------------------
    allocate (nc% dims(NDIMS))
    allocate (nc% vars(NVARS))
    !----------------------
    ! define the dimensions
    !----------------------
    call add_dim (nc, 'd_hdr',      0, p_h)
    call add_dim (nc, 'd_body',     0, p_b)
    call add_dim (nc, 'd_radar',    0, p_r)
    call add_dim (nc, 'd_veri',     0, p_v,  unlimited=.true.)
    call add_dim (nc, 'd_2',        2, p_2)
    call add_dim (nc, 'd_3',        3, p_3)
    call add_dim (nc, 'char10',  MLEN, p_m)
    call add_dim (nc, 'char12', IDLEN, p_i)
    call add_dim (nc, 'char64',  DLEN, p_d)

    !===================================
    ! define the variables of the header
    !===================================
    call add_var (nc, 'i_body'       ,NF_INT,   (/    p_h/),'index of 1st entry in report body')
    call add_var (nc, 'l_body'       ,NF_SHORT, (/    p_h/),'number of entries in report body')
    call add_var (nc, 'n_level'      ,NF_SHORT, (/    p_h/),'number of levels in report')
    call add_var (nc, 'i_spec'       ,NF_INT,   (/    p_h/),'index of 1st entry in radar table',&
                                                      opt = 'RADAR'                             )
    call add_var (nc, 'l_spec'       ,NF_SHORT, (/    p_h/),'number of entries in radar table',&
                                                      opt = 'RADAR'                            )
    call add_var (nc, 'data_category',NF_SHORT, (/    p_h/),'BUFR4 data category')
    call add_var (nc, 'sub_category' ,NF_SHORT, (/    p_h/),'BUFR4 data sub-category')
    call add_var (nc, 'center'       ,NF_SHORT, (/    p_h/),'station processing center')
    call add_var (nc, 'sub_center'   ,NF_SHORT, (/    p_h/),'station processing sub-center')
    call add_var (nc, 'obstype'      ,NF_BYTE,  (/    p_h/),'observation type',&
                                                     table = obstype)
    call add_var (nc, 'codetype'     ,NF_SHORT, (/    p_h/),'code type',&
                                                     table = codetype)
    call add_var (nc, 'ident'        ,NF_INT,   (/    p_h/),'station or satellite id as integer')
    call add_var (nc, 'statid'       ,NF_CHAR,  (/p_m,p_h/),'station id as character string')
    call add_var (nc, 'lat'          ,NF_FLOAT, (/    p_h/),'latitude',&
                                                    units = 'degree')
    call add_var (nc, 'lon'          ,NF_FLOAT, (/    p_h/),'longitude',&
                                                    units = 'degree')
    call add_var (nc, 'time'         ,NF_SHORT, (/    p_h/),'observation minus reference time',&
                                                    units = 'min')
    call add_var (nc, 'time_nomi'    ,NF_SHORT, (/    p_h/),'nominal (synoptic) minus reference time',&
                                                    units = 'min')
    call add_var (nc, 'time_dbase'   ,NF_SHORT, (/    p_h/),'data base minus reference time',&
                                                    units = 'min')
    call add_var (nc, 'z_station'    ,NF_INT,   (/    p_h/),'station height',&
                                                    units = 'm', invalid=NF_FILL_SHORT)
    call add_var (nc, 'z_modsurf'    ,NF_SHORT, (/    p_h/),'model surface height',&
                                                    units = 'm')
    call add_var (nc, 'r_state'      ,NF_BYTE,  (/    p_h/),'status of the report',&
                                                     table = status)
    call add_var (nc, 'r_flags'      ,NF_INT,   (/    p_h/),'report quality check flags',&
                                                     table = flags)
    call add_var (nc, 'r_check'      ,NF_BYTE,  (/    p_h/),'check which raised the report status flag value',&
                                                     table = flags)
    call add_var (nc, 'sta_corr'     ,NF_BYTE,  (/    p_h/),'station correction indicator')
    call add_var (nc, 'index_x'      ,NF_INT,   (/    p_h/),'index x of model grid point assigned to report')
    call add_var (nc, 'index_y'      ,NF_SHORT, (/    p_h/),'index y of model grid point assigned to report')
    call add_var (nc, 'mdlsfc'       ,NF_BYTE,  (/    p_h/),'model surface characteristics',&
                                                     table = surf_char)
    call add_var (nc, 'instype'      ,NF_SHORT, (/    p_h/),'station type or satellite instrument type',&
                                                     table = satsens                                    )
    call add_var (nc, 'sun_zenit'    ,NF_FLOAT, (/    p_h/),'sun zenith angle')

    !--------------------------
    ! observation type specific
    !--------------------------
    call add_var (nc, 'retrtype' ,NF_SHORT, (/    p_h/),'retrieval type',&
                                                  table = retrtype,         &
                                                    opt = 'SATOB'           )
    call add_var (nc, 'phase'    ,NF_SHORT,  (/    p_h/),'aircraft phase',&
                                                    opt = 'AIREP RAD GPSRO SCATT')
    call add_var (nc, 'tracking' ,NF_BYTE,  (/    p_h/),'tracking technique',&
                                                    opt = 'TEMP PILOT AIREP')
    call add_var (nc, 'meas_type',NF_BYTE,  (/    p_h/),'type of measuring equipment used',&
                                                    opt = 'TEMP PILOT')
    call add_var (nc, 'rad_corr' ,NF_BYTE,  (/    p_h/),'solar and infrared radiation correction',&
                                                    opt = 'TEMP')
    call add_var (nc, 'ct_nwc',   NF_BYTE,  (/    p_h/),'Cloud Type according to NWC SAF',&
                                                  table = ct_nwc,                         &
                                                    opt = 'RAD',                          &
                                                invalid = CT_NOPR                         )
    call add_var (nc, 'ch_nwc',NF_FLOAT, (/    p_h/),'Cloud Height according to NWC SAF', &
                                                  units = 'm',                            &
                                                    opt = 'RAD'                           )
    call add_var (nc, 'flg_1dvar',NF_BYTE,  (/    p_h/),'1dvar processing flags',&
                                                  table = flg_1dvar,                &
                                                    opt = 'RAD'                   )
    call add_var (nc, 'flg_cld'  ,NF_BYTE,  (/    p_h/),'1dvar cloud flags',&
                                                  table = flg_cld,             &
                                                    opt = 'RAD'              )
    call add_var (nc, 'surftype' ,NF_SHORT, (/    p_h/),'surface type bit pattern',&
                                                  table = surftype,                &
                                                    opt = 'RAD'                    )
    call add_var (nc, 'sat_zenit',NF_FLOAT, (/    p_h/),'satellite zenith angle',  &
                                                    opt = 'RAD'                    )
    call add_var (nc, 'varno_back',NF_SHORT,(/    p_h/),'type of observed quantity',&
                                                  table = varno,                    &
                                                    opt = 'RADAR'                   )
    call add_var (nc, 'vnyquist',NF_FLOAT,(/      p_h/),'nyquist frequency',&
                                                  units = 'm/s',            &
                                                    opt = 'RADAR'           )
    call add_var (nc, 'spec_r_flags',  NF_INT,(/  p_h/),'observation type specific flags',&
                                                    opt = 'RADAR'                         )
    !---------------
    ! 3DVAR specific
    !---------------
    call add_var (nc, 'obs_id'   ,NF_INT,   (/    p_h/),'unique observation id',               opt = '3DVAR')
    call add_var (nc, 'source'   ,NF_BYTE,  (/    p_h/),'input file number',                   opt = '3DVAR')
    call add_var (nc, 'record'   ,NF_INT,   (/    p_h/),'record number in the input file',     opt = '3DVAR')
    call add_var (nc, 'subset'   ,NF_SHORT, (/    p_h/),'subset number in the record',         opt = '3DVAR')
    call add_var (nc, 'dbkz'     ,NF_INT,   (/    p_h/),'DWD data base id',                    opt = '3DVAR')
    call add_var (nc, 'index_d'  ,NF_BYTE,  (/    p_h/),'model grid diamond index assigned to report',      &
                                                                                               opt = '3DVAR')
    !=================================
    ! define the variables of the body
    !=================================
    call add_var (nc, 'varno'    ,NF_SHORT, (/    p_b/),'type of the observed quantity',&
                                                  table = varno)
    call add_var (nc, 'obs'      ,NF_FLOAT, (/    p_b/),'bias corrected observation')
    call add_var (nc, 'bcor'     ,NF_FLOAT, (/    p_b/),'bias correction, corrected minus observed')
    call add_var (nc, 'level'    ,NF_FLOAT, (/    p_b/),'level of observation')
    call add_var (nc, 'level_typ',NF_SHORT, (/    p_b/),'type of level information',&
                                                  table = varno                        )
    call add_var (nc, 'level_sig',NF_SHORT, (/    p_b/),'level significance',&
                                                  table = level_sig             )
    call add_var (nc, 'state'    ,NF_BYTE,  (/    p_b/),'status of the observation',&
                                                  table = status                       )
    call add_var (nc, 'flags'    ,NF_INT,   (/    p_b/),'observation quality check flags',&
                                                  table = flags                              )
    call add_var (nc, 'check'    ,NF_BYTE,  (/    p_b/),'check which raised the observation status flag value',&
                                                  table = flags                                                   )
    call add_var (nc, 'e_o'      ,NF_FLOAT, (/    p_b/),'observational error')
    call add_var (nc, 'qual'     ,NF_SHORT, (/    p_b/),'observation confidence from data provider')
    !--------------------------
    ! observation type specific
    !--------------------------
!   call add_var (nc, 'pcc'      ,NF_BYTE,  (/    p_b/),'observation confidence from data provider',&
!                                                units = '%',                                       &
!                                                  opt = 'SATOB AIREP'                              )
!   call add_var (nc, 'snr'      ,NF_SHORT, (/    p_b/),'signal to noise ratio (wind profiler)'    ,&
!                                                units = 'db',                                      &
!                                                  opt = 'PILOT'                                    )
    call add_var (nc, 'plevel'   ,NF_FLOAT, (/    p_b/),'nominal pressure level'                   ,&
                                                 units = 'Pa'                                       )
    call add_var (nc, 'dlat'     ,NF_FLOAT, (/    p_b/),'latitude (drifting)',                      &
                                                    units = 'degree', opt = 'GPSRO TEMP RADAR'      )
    call add_var (nc, 'dlon'     ,NF_FLOAT, (/    p_b/),'longitude (drifting)',                     &
                                                    units = 'degree', opt = 'GPSRO TEMP RADAR'      )
    call add_var (nc, 'accuracy' ,NF_FLOAT, (/    p_b/),'accuracy from data provider'              ,&
                                                   opt = 'GPSRO PILOT') ! wind profiler so far
!   call add_var (nc,'spec_index',NF_SHORT, (/    p_b/),'observation index in case of sparse coverage',&
    call add_var (nc,'spec_index',NF_INT,   (/    p_b/),'observation index in case of sparse coverage',&
                                                   opt = 'RADAR')
    call add_var (nc, 'pl_width' ,NF_SHORT, (/   p_b/),'width of sensitivity function'             ,&
                                                   opt = 'IASI'                                     )
    call add_var (nc, 'l2c'      ,NF_FLOAT, (/   p_b/),'level to channel assignment'               ,&
                                                   opt = 'IASI'                                     )
    call add_var (nc, 'azimuth'  ,NF_FLOAT, (/   p_b/),'azimuth'                                   ,&
                                                   units = 'degree', opt = 'GPSGB GPSRO'            )
    call add_var (nc, 'obs_par_1',NF_FLOAT, (/   p_b/),'observation specific parameter 1'          ,&
                                                   opt = 'RAD'                                      )
!   call add_var (nc, 'obs_par_2',NF_FLOAT, (/   p_b/),'observation specific parameter 2'          ,&
!                                                  opt = 'GPSGB'                                    )

!   !---------------
!   ! 3DVAR specific
!   !---------------
!   call add_var (nc, 'w_qc'     ,NF_FLOAT, (/    p_b/),'Variational Quality Control weight', opt = '3DVAR')

    !---------------
    ! LETKF specific
    !---------------
    call add_var (nc, 'v_loc'     ,NF_FLOAT, (/    p_b/),'vertical localisation radius',   opt = 'LETKF')
    call add_var (nc, 'h_loc'     ,NF_FLOAT, (/    p_b/),'horizontal localisation radius', opt = 'LETKF')

    !--------------------------
    ! RADAR type specific table
    !--------------------------
    call add_var (nc,'radar_azimuth'    ,NF_FLOAT, (/ p_r/),'azimuth of ray'    ,units = 'degree' ,opt = 'RADAR')
    call add_var (nc,'radar_elevation'  ,NF_FLOAT, (/ p_r/),'elevation of ray'  ,units = 'degree' ,opt = 'RADAR')
    call add_var (nc,'radar_nrange'     ,NF_SHORT, (/ p_r/),'number of bins per ray'              ,opt = 'RADAR')
    call add_var (nc,'radar_range_start',NF_FLOAT, (/ p_r/),'distance of 1.bin' ,units = 'm'      ,opt = 'RADAR')
    call add_var (nc,'radar_drange'     ,NF_FLOAT, (/ p_r/),'distance of bins'  ,units = 'm'      ,opt = 'RADAR')
    call add_var (nc,'radar_nbody'      ,NF_SHORT, (/ p_r/),'number of corresponding body entries',opt = 'RADAR')

    !=============================
    ! define the verification data
    !=============================
    call add_var (nc, 'veri_data',NF_FLOAT, (/p_b,p_v/),'modelled quantity (as indicated by veri_ens_member)')

    !==================================
    ! define the verification meta data
    !==================================
    call add_var(nc,'veri_model'        ,NF_CHAR, (/p_m,p_v/),'model used for verification, e.g. COSMO, GME ...')
    call add_var(nc,'veri_run_type'     ,NF_BYTE,     (/p_v/),'type of model run',&
                                                                  table = runtype)
    call add_var(nc,'veri_run_class'    ,NF_BYTE,     (/p_v/),'class of model run',&
                                                                  table = runclass)
    call add_var(nc,'veri_initial_date' ,NF_CHAR, (/p_i,p_v/),'start of verification period',&
                                                                  units = 'yyyymmddhhmm')
    call add_var(nc,'veri_forecast_time',NF_INT,      (/p_v/),'forecast time at verification_ref_time',&
                                                                  units = 'hhmm')
    call add_var(nc,'veri_resolution'   ,NF_FLOAT,(/p_2,p_v/),'model resolution, seperately for the x- and y-dimension',&
                                                                  units = 'degree')
    call add_var(nc,'veri_domain_size'  ,NF_INT,  (/p_3,p_v/),'domain size: nx,ny,nz for COSMO; ni,ni,nz for GME/ICON')
    call add_var(nc,'veri_description'  ,NF_CHAR, (/p_d,p_v/),'more detailed description than veri_run_type')
    call add_var(nc,'veri_ens_member'   ,NF_INT,      (/p_v/),'ensemble member number, special meaning for non-positive values',&
                                                                  table = ensmem)
    call add_var(nc,'veri_exp_id'       ,NF_INT,      (/p_v/),'experiment Id of the verification run')
    call add_var(nc,'veri_operator_flag',NF_INT,      (/p_v/),'observation operator flags')
  end subroutine setup_fdbk
!------------------------------------------------------------------------------
  subroutine create_fdbk (fb, path, model, version, institution, n_hdr,    &
                          n_body, refdate, reftime, start, end, resolution,&
                          domain, comment, time, runtime,                  &
                          pole, lower_left, upper_right, n_radar,          &
                          opt, create                                      )
  !-------------------------------------------------------------
  ! create new feedback file
  !   task 1: set up derived type according to actual parameters
  !   task 2: create file and write global attributes
  !-------------------------------------------------------------
  type(t_fdbk)               ,intent(inout) :: fb         ! feedback file data type
  character(len=*)           ,intent(in) :: path          ! pathname
  character(len=*)           ,intent(in) :: model         ! model string
  character(len=*)           ,intent(in) :: version       ! model version
  character(len=*)           ,intent(in) :: institution   ! institution string
  integer                    ,intent(in) :: n_hdr         ! allocated size of header
  integer                    ,intent(in) :: n_body        ! allocated size of body
  integer                    ,intent(in) :: refdate       ! reference time yyyymmdd
  integer                    ,intent(in) :: reftime       ! reference time hhmm
  integer                    ,intent(in) :: start         ! verification start (minutes)
  integer                    ,intent(in) :: end           ! verification stop  (minutes)
  real(sp)                   ,intent(in) :: resolution (2)! model resolution    (degree)
  integer                    ,intent(in) :: domain     (3)! domain size          (x,y,z)
  character(len=*)           ,intent(in) :: comment       ! comment  (for history)
  character(len=*)           ,intent(in) :: time          ! time     (for history)
  character(len=*) ,optional ,intent(in) :: runtime       ! run time (for history)
  real(sp)         ,optional ,intent(in) :: pole       (2)! location of pole    (degree)
  real(sp)         ,optional ,intent(in) :: lower_left (2)! lower left  (lat,lon degree)
  real(sp)         ,optional ,intent(in) :: upper_right(2)! upper right (lat,lon degree)
  integer          ,optional ,intent(in) :: n_radar       ! allocated radar specific table
  character(len=*) ,optional ,intent(in) :: opt           ! flag optional variables
  logical          ,optional ,intent(in) :: create        ! if .false. postpone task 2

    character(len=10) :: model10
    logical           :: lcreate
    integer           :: nradar
    model10 = model
    nradar  = 1;      if (present(n_radar)) nradar  = n_radar
    lcreate = .true.; if (present(create))  lcreate = create
    !---------------------------------------------------------
    ! create the NetCDF file (define dimensions and variables)
    !---------------------------------------------------------
    fb% nc% path    = path
    where(fb%nc% dims% name=='d_hdr')   fb% nc% dims% len = n_hdr
    where(fb%nc% dims% name=='d_body')  fb% nc% dims% len = n_body
    where(fb%nc% dims% name=='d_radar') fb% nc% dims% len = nradar

    if (lcreate) then
      call create_netcdf_file (fb%nc, opt=opt)
      if (fb%nc% error /= 0) return
    endif

    !-----------------------------------------
    ! set up global attributes in derived type
    !-----------------------------------------
    fb%     n_hdr   = 0                    ! size of header actually used
    fb%     n_body  = 0                    ! size of body   actually used
    fb%     n_radar = 0                    ! size of radar table     used
    fb%     n_veri  = 0                    ! number of verification runs
    fb%     title   = model10//' Verification Data'
    fb% institution = institution
    fb%     source  = model10//' '//version
    fb%     version = FEEDBACK_VERSION
    fb%     refdate = refdate
    fb%     reftime = reftime
    fb%     start   = start
    fb%     end     = end
    fb% resolution  = resolution
    fb%     domain  = domain
    if (present(pole))        fb% pole        = pole
    if (present(lower_left))  fb% lower_left  = lower_left
    if (present(upper_right)) fb% upper_right = upper_right
    call add_history (fb, model, time, comment, runtime, write=create)

    !--------------------------------
    ! write attributes to NetCDF file
    !--------------------------------
    if (lcreate) call write_global_attributes (fb)

  end subroutine create_fdbk
!------------------------------------------------------------------------------
  subroutine write_global_attributes (fb)
  type (t_fdbk), intent(inout) :: fb

    integer :: n
    integer :: ir
    n  = size     (fb% history)
    ir = nf_redef (fb% nc% ncid)
    if (ir /= NF_NOERR .and. ir /= NF_EINDEFINE) then
      write(0,*) 'FAILURE setting "redefine" mode in write_global_attributes'
      call model_abort (-1,-1,'failure setting "redefine" mode',&
                        'write_global_attributes')
    endif

    !------------------------
    ! write global attributes
    !------------------------
    fb%nc% error = nf_put_att_text (fb%nc% ncid, NF_GLOBAL, 'title', TLEN, &
                                    fb% title)

    fb%nc% error = nf_put_att_text (fb%nc% ncid, NF_GLOBAL, 'institution', &
                                    len_trim(fb% institution), fb% institution)

    fb%nc% error = nf_put_att_text (fb%nc% ncid, NF_GLOBAL, 'source',      &
                                    len_trim(fb% source),      fb% source)

    fb%nc% error = nf_put_att_text (fb%nc% ncid, NF_GLOBAL, 'history',n*HLEN, fb% history)

    fb%nc% error = nf_put_att_text (fb%nc% ncid, NF_GLOBAL, 'file_version_number',5, fb% version)

    fb%nc% error = nf_put_att_int  (fb%nc% ncid, NF_GLOBAL, 'n_hdr', NF_INT, 1, fb% n_hdr)

    fb%nc% error = nf_put_att_int  (fb%nc% ncid, NF_GLOBAL, 'n_body', NF_INT, 1, fb% n_body)

    fb%nc% error = nf_put_att_int  (fb%nc% ncid, NF_GLOBAL, 'n_radar', NF_INT, 1, fb% n_radar)

    fb%nc% error = nf_put_att_int  (fb%nc% ncid, NF_GLOBAL, 'verification_ref_date',NF_INT, 1, fb% refdate)

    fb%nc% error = nf_put_att_int  (fb%nc% ncid, NF_GLOBAL, 'verification_ref_time',NF_INT, 1, fb% reftime)

    fb%nc% error = nf_put_att_int  (fb%nc% ncid, NF_GLOBAL, 'verification_start',NF_INT, 1, fb% start)

    fb%nc% error = nf_put_att_int  (fb%nc% ncid, NF_GLOBAL, 'verification_end',NF_INT, 1, fb% end)

    fb%nc% error = nf_put_att_real (fb%nc% ncid, NF_GLOBAL, 'resolution',NF_FLOAT, 2, fb% resolution)

    fb%nc% error = nf_put_att_int  (fb%nc% ncid, NF_GLOBAL, 'domain_size',NF_INT, 3, fb% domain)

    if (any (fb% pole        /= 0._sp .or. &
             fb% lower_left  /= 0._sp .or. &
             fb% upper_right /= 0._sp)     ) then

      fb%nc% error = nf_put_att_real (fb%nc% ncid, NF_GLOBAL, 'pole_lat_lon',NF_FLOAT, 2, fb% pole)

      fb%nc% error = nf_put_att_real (fb%nc% ncid, NF_GLOBAL, 'lower_left_lat_lon',NF_FLOAT, 2, fb% lower_left)

      fb%nc% error = nf_put_att_real (fb%nc% ncid, NF_GLOBAL, 'upper_right_lat_lon',NF_FLOAT, 2, fb% upper_right)

    endif

  end subroutine write_global_attributes
!------------------------------------------------------------------------------
  subroutine read_meta (fb)
  type (t_fdbk), intent(inout) :: fb

    integer                :: n, varid, i, d_veri
    real(sp)  ,allocatable :: resolution  (:,:)
    integer   ,allocatable :: domain_size (:,:)

    fb%nc% error = nf_get_att_text (fb%nc% ncid, NF_GLOBAL, 'title', fb% title)
    if (fb%nc% error /= NF_NOERR) return

    fb%nc% error = nf_get_att_text (fb%nc% ncid, NF_GLOBAL, 'institution', fb% institution)
    if (fb%nc% error /= NF_NOERR) return

    fb%nc% error = nf_get_att_text (fb%nc% ncid, NF_GLOBAL, 'source', fb% source)
    if (fb%nc% error /= NF_NOERR) return

    fb%nc% error = nf_inq_attlen   (fb%nc% ncid, NF_GLOBAL, 'history', n)
    if (fb%nc% error /= NF_NOERR) return

    n = (n+(HLEN-1))/HLEN
    if(associated(fb% history)) deallocate (fb% history)
    allocate (fb% history (n))
    fb%nc% error = nf_get_att_text (fb%nc% ncid, NF_GLOBAL, 'history', fb% history)
    if (fb%nc% error /= NF_NOERR) return

    fb%nc% error = nf_get_att_text (fb%nc% ncid, NF_GLOBAL, 'file_version_number', fb% version)
    if (fb%nc% error /= NF_NOERR) return

    fb%nc% error = nf_get_att_int  (fb%nc% ncid, NF_GLOBAL, 'n_hdr', fb% n_hdr)
    if (fb%nc% error /= NF_NOERR) return

    fb%nc% error = nf_get_att_int  (fb%nc% ncid, NF_GLOBAL, 'n_body', fb% n_body)
    if (fb%nc% error /= NF_NOERR) return

    fb%nc% error = nf_get_att_int  (fb%nc% ncid, NF_GLOBAL, 'n_radar', fb% n_radar)
    if (fb%nc% error /= NF_NOERR) fb% n_radar = 0

    fb%nc% error = nf_get_att_int  (fb%nc% ncid, NF_GLOBAL, 'verification_ref_date', fb% refdate)
    if (fb%nc% error /= NF_NOERR) return

    fb%nc% error = nf_get_att_int  (fb%nc% ncid, NF_GLOBAL, 'verification_ref_time', fb% reftime)
    if (fb%nc% error /= NF_NOERR) return

    fb%nc% error = nf_get_att_int  (fb%nc% ncid, NF_GLOBAL, 'verification_start', fb% start)
    if (fb%nc% error /= NF_NOERR) return

    fb%nc% error = nf_get_att_int  (fb%nc% ncid, NF_GLOBAL, 'verification_end', fb% end)
    if (fb%nc% error /= NF_NOERR) return

    fb%nc% error = nf_get_att_real (fb%nc% ncid, NF_GLOBAL, 'resolution', fb% resolution)
    if (fb%nc% error /= NF_NOERR) return

    fb%nc% error = nf_get_att_int  (fb%nc% ncid, NF_GLOBAL, 'domain_size', fb% domain)
    if (fb%nc% error /= NF_NOERR) return

    fb%nc% error = nf_get_att_real (fb%nc% ncid, NF_GLOBAL, 'pole_lat_lon', fb% pole)
    if (fb%nc% error /= NF_NOERR) fb% pole = 0._sp

    fb%nc% error = nf_get_att_real (fb%nc% ncid, NF_GLOBAL, 'lower_left_lat_lon', fb% lower_left)
    if (fb%nc% error /= NF_NOERR) fb% lower_left = 0._sp

    fb%nc% error = nf_get_att_real (fb%nc% ncid, NF_GLOBAL, 'upper_right_lat_lon', fb% upper_right)
    if (fb%nc% error /= NF_NOERR) fb% upper_right = 0._sp

    !----------------------------
    ! read verification meta data
    !----------------------------
    fb%nc% error = nf_inq_dimid (fb%nc% ncid, 'd_veri',d_veri)
    if (fb%nc% error /= NF_NOERR) return
    fb%nc% error = nf_inq_dimlen (fb%nc% ncid, d_veri, fb% n_veri)
    if (fb%nc% error /= NF_NOERR) return

    if (associated (fb% veri)) deallocate (fb% veri)
    allocate (fb% veri        (   fb% n_veri))
    allocate (    resolution  (2, fb% n_veri))
    allocate (    domain_size (3, fb% n_veri))

    fb%nc% error = nf_inq_varid    (fb%nc% ncid, 'veri_model', varid)
    fb%nc% error = nf_get_var_text (fb%nc% ncid, varid, fb% veri% model)
    call blank_nul_chars (fb% veri% model)

    fb%nc% error = nf_inq_varid    (fb%nc% ncid, 'veri_run_type', varid)
    fb%nc% error = nf_get_var_int  (fb%nc% ncid, varid, fb% veri% run_type)

    fb%nc% error = nf_inq_varid    (fb%nc% ncid, 'veri_run_class', varid)
    fb%nc% error = nf_get_var_int  (fb%nc% ncid, varid, fb% veri% run_class)

    fb%nc% error = nf_inq_varid    (fb%nc% ncid, 'veri_initial_date', varid)
    fb%nc% error = nf_get_var_text (fb%nc% ncid, varid, fb% veri% initial_date)

    fb%nc% error = nf_inq_varid    (fb%nc% ncid, 'veri_forecast_time', varid)
    fb%nc% error = nf_get_var_int  (fb%nc% ncid, varid, fb% veri% forecast_time)

    fb%nc% error = nf_inq_varid    (fb%nc% ncid, 'veri_description', varid)
    fb%nc% error = nf_get_var_text (fb%nc% ncid, varid, fb% veri% description)

    fb%nc% error = nf_inq_varid    (fb%nc% ncid, 'veri_ens_member', varid)
    fb%nc% error = nf_get_var_int  (fb%nc% ncid, varid, fb% veri% ens_member)

    fb%nc% error = nf_inq_varid    (fb%nc% ncid, 'veri_exp_id', varid)
    fb%nc% error = nf_get_var_int  (fb%nc% ncid, varid, fb% veri% exp_id)

    fb% veri% operator_flag = 0
    fb%nc% error = nf_inq_varid    (fb%nc% ncid, 'veri_operator_flag', varid)
    if (fb%nc% error == NF_NOERR)                                              & 
    fb%nc% error = nf_get_var_int  (fb%nc% ncid, varid, fb% veri% operator_flag)

    resolution   = 0._sp
    fb%nc% error = nf_inq_varid    (fb%nc% ncid, 'veri_resolution', varid)
    fb%nc% error = nf_get_var_real (fb%nc% ncid, varid, resolution)

    domain_size  = 0
    fb%nc% error = nf_inq_varid    (fb%nc% ncid, 'veri_domain_size', varid)
    fb%nc% error = nf_get_var_int  (fb%nc% ncid, varid, domain_size)

    do i = 1, fb% n_veri
      fb% veri(i)% resolution  = resolution  (:,i)
      fb% veri(i)% domain_size = domain_size (:,i)
    end do

  end subroutine read_meta
!------------------------------------------------------------------------------
  subroutine print_fdbk (fdbk)
  type (t_fdbk)     ,intent(in) :: fdbk

    integer :: i
    integer :: hl1

    write(6,'(a)')           repeat('_',79)
    write(6,'()')
    write(6,'(a)')         '                      Feedback File Inventory'
    write(6,'(a,a)')       '  version          = ',trim(fdbk% version)
    write(6,'(a,a )')      '  path             = ',trim(fdbk% nc% path)
    write(6,'()')
    write(6,'(a,a)')       '  title            = ',trim(fdbk% title)
    write(6,'(a,a)')       '  institution      = ',trim(fdbk% institution)
    write(6,'(a,a )')      '  source           = ',trim(fdbk% source)
    if (associated(fdbk% history)) then
    do i = 1, size(fdbk% history)
    if (i==1) then
    hl1 = len(fdbk% history(1)) - 1
    write(6,'(a,a)')       '  history          = ',trim(fdbk% history(i)(1:hl1))
    else
    write(6,'(a,a)')       '                     ',trim(fdbk% history(i)(1:hl1))
    endif
    end do
    end if
    write(6,'()')
    write(6,'(a,i8)')       '  refdate          = ',     fdbk% refdate
    write(6,'(a,i8.4)')     '  reftime          = ',     fdbk% reftime
    write(6,'(a,i8)')       '  start (min)      = ',     fdbk% start
    write(6,'(a,i8)')       '  end   (min)      = ',     fdbk% end
    write(6,'()')
    write(6,'(a,2f8.3)')    '  resolution       = ',     fdbk% resolution
    write(6,'(a,3i8)')      '  domain           = ',     fdbk% domain
    write(6,'(a,2f8.3)')    '  pole             = ',     fdbk% pole
    write(6,'(a,2f8.3)')    '  lower_left       = ',     fdbk% lower_left
    write(6,'(a,2f8.3)')    '  upper_right      = ',     fdbk% upper_right
    write(6,'()')
    write(6,'(a,i8)')       '  n_hdr            = ',     fdbk% n_hdr
    write(6,'(a,i8)')       '  n_body           = ',     fdbk% n_body
    write(6,'(a,i8)')       '  n_radar          = ',     fdbk% n_radar
    write(6,'(a,i8)')       '  n_veri           = ',     fdbk% n_veri
    if (associated (fdbk% veri)) then
     do i=1,fdbk% n_veri
      write(6,'()')
      write(6,'(i3,a,a)')   i, ': model         = ',trim(fdbk% veri(i)% model)
      write(6,'(3x,a,2f8.3)')  '  resolution    = ',     fdbk% veri(i)% resolution
      write(6,'(3x,a,3i8)')    '  domain_size   = ',     fdbk% veri(i)% domain_size
      write(6,'(3x,a,a)')      '  initial_date  = ',     fdbk% veri(i)% initial_date
      write(6,'(3x,a,i12.4)')  '  forecast_time = ',     fdbk% veri(i)% forecast_time
      write(6,'(3x,a,i8)')     '  run_type      = ',     fdbk% veri(i)% run_type
      write(6,'(3x,a,i8)')     '  run_class     = ',     fdbk% veri(i)% run_class
      write(6,'(3x,a,i8)')     '  ens_member    = ',     fdbk% veri(i)% ens_member
      write(6,'(3x,a,a)')      '  description   = ',     fdbk% veri(i)% description
      write(6,'(3x,a,i8)')     '  exp_id        = ',     fdbk% veri(i)% exp_id
      write(6,'(3x,a,i8)')     '  operator_flag = ',     fdbk% veri(i)% operator_flag
     end do
    endif
    write(6,'()')
    write(6,'(a)')           repeat('_',79)

  end subroutine print_fdbk
!------------------------------------------------------------------------------
  subroutine add_history (fdbk, model, starttime, runtype, runtime, write)
  !---------------------------------------
  ! add history entry to global attributes
  !  1) store in derived type
  !  2) actually write to file
  !---------------------------------------
  type (t_fdbk)     ,intent(inout)        :: fdbk
  character(len=*)  ,intent(in)           :: model
  character(len=*)  ,intent(in)           :: starttime
  character(len=*)  ,intent(in)           :: runtype
  character(len=*)  ,intent(in) ,optional :: runtime   ! (yyyymmddhhmm)
  logical           ,intent(in) ,optional :: write     ! if .false. postpone 2)

    character(len=8)             :: yyyymmdd
    character(len=10)            :: hhmmss_mmm
    character(len=12)            :: yyyymmddhhmm
    integer                      :: n
    integer                      :: ir
    character(len=HLEN) ,pointer :: tmp (:)
    logical                      :: lwrite

    !----------------------------
    ! process optional parameters
    !----------------------------
    lwrite = .true.; if (present(write)) lwrite = write
    if (present (runtime)) then
      yyyymmddhhmm = runtime
    else
      call date_and_time (yyyymmdd, hhmmss_mmm)
      yyyymmddhhmm = yyyymmdd//hhmmss_mmm
    endif

    !---------------------------
    ! (re)allocate history array
    !---------------------------
    if (.not.associated(fdbk% history)) then
      n = 1
      allocate (fdbk% history (1))
    else
      n = size (fdbk% history) + 1
      tmp => fdbk% history
      allocate (fdbk% history (n))
      fdbk% history (1:n-1) = tmp
      deallocate (tmp)
    endif

    !--------------------
    ! add line to history
    !--------------------
    if(n>1) fdbk% history(n-1) (HLEN:HLEN) = achar (10) ! newline
    write (fdbk% history(n), '(a12,1x,a8,1x,a12,1x,a)') &
      yyyymmddhhmm, model, starttime, runtype

    !--------------
    ! write history
    !--------------
    if (lwrite) then
      ir = nf_redef (fdbk% nc% ncid)
      if (ir /= NF_NOERR .and. ir /= NF_EINDEFINE) then
        write(0,*) 'FAILURE setting "redefine" mode in add_history'
        write(0,*) trim(nf_strerror (ir))
        call model_abort (-1,-1,'failure setting "redefine" mode','add_history')
      endif
      fdbk%nc% error = nf_put_att_text (fdbk%nc% ncid, NF_GLOBAL,      &
                                        'history',n*HLEN, fdbk% history)
    endif

  end subroutine add_history
!------------------------------------------------------------------------------
  subroutine add_verification (fdbk, model, run_type, run_class, initial_date,&
                               fc_time, resolution, domain_size, description, &
                               ens_member, exp_id, id_veri, operator_flag,    &
                               replace, ierr                                  )
  !-------------------------------------------------------------------
  ! Add an verification entry to the feedback file.  All meta data
  ! entries (variables veri_*) are written.  Only the actual data
  ! (veri_data) has to be written in a subsequent step.  If 'replace'
  ! is present and .gt. 0, the entry with the respective entry is
  ! replaced. Otherwise the entry is appended. The index actually used
  ! is returned in 'replace'. By passing 'ierr' a program abort due
  ! to an error condition may be captured. In this case 'ierr' returns
  ! a value /= zero.
  !-------------------------------------------------------------------
  type (t_fdbk)     ,intent(inout) :: fdbk
  character(len=*)  ,intent(in)    :: model          ! verification model
  integer           ,intent(in)    :: run_type       ! run type
  integer           ,intent(in)    :: run_class      ! run class
  character(len=12) ,intent(in)    :: initial_date   ! yyyymmddhhmm
  integer           ,intent(in)    :: fc_time        ! hhhmm
  real(sp)          ,intent(in)    :: resolution (2) ! (degree) for x,y
  integer           ,intent(in)    :: domain_size(3) ! nx,ny,nz or ni,ni,nz
  character(len=*)  ,intent(in)    :: description    !
  integer           ,intent(in)    :: ens_member     ! ensemble member number
  integer           ,intent(in)    :: exp_id         ! experiment id
  integer           ,intent(out)   :: id_veri        ! var.id. of 'veri_data'
  integer ,optional ,intent(in)    :: operator_flag  ! observ. operator flag
  integer ,optional ,intent(inout) :: replace        ! replace entry
  integer ,optional ,intent(out)   :: ierr           ! error return parameter

    integer           :: vid
    character(len=10) :: c10
    character(len=64) :: c64
    character(len=24) :: name
    integer           :: n
    integer           :: op_flag

    n                        = 0
    if (present (replace)) n = replace
    if (n <= 0)            n = fdbk% n_veri + 1
    if (present (replace)) replace = n
    fdbk% n_veri = max (fdbk% n_veri, n)

    fdbk%nc% error = nf_enddef (fdbk% nc% ncid)

    name = 'veri_model'
    vid  = get_varid (fdbk, name)
    c10  = model
    fdbk%nc% error = nf_put_vara_text (fdbk% nc% ncid, vid, (/1,n/), (/10,1/), c10)
    if (fdbk%nc% error /= NF_NOERR) goto 999

    name = 'veri_run_type'
    vid  = get_varid (fdbk, name)
    fdbk%nc% error = nf_put_vara_int  (fdbk% nc% ncid, vid, n, 1, run_type)
    if (fdbk%nc% error /= NF_NOERR) goto 999

    name = 'veri_run_class'
    vid  = get_varid (fdbk, name)
    fdbk%nc% error = nf_put_vara_int  (fdbk% nc% ncid, vid, n, 1, run_class)
    if (fdbk%nc% error /= NF_NOERR) goto 999

    name = 'veri_initial_date'
    vid  = get_varid (fdbk, name)
    fdbk%nc% error = nf_put_vara_text (fdbk% nc% ncid, vid, (/1,n/), (/12,1/), initial_date)
    if (fdbk%nc% error /= NF_NOERR) goto 999

    name = 'veri_forecast_time'
    vid  = get_varid (fdbk, name)
    fdbk%nc% error = nf_put_vara_int  (fdbk% nc% ncid, vid, n, 1, fc_time)
    if (fdbk%nc% error /= NF_NOERR) goto 999

    name = 'veri_resolution'
    vid  = get_varid (fdbk, name)
    fdbk%nc% error = nf_put_vara_real (fdbk% nc% ncid, vid, (/1,n/), (/2,1/), resolution)
    if (fdbk%nc% error /= NF_NOERR) goto 999

    name = 'veri_domain_size'
    vid  = get_varid (fdbk, name)
    fdbk%nc% error = nf_put_vara_int  (fdbk% nc% ncid, vid, (/1,n/), (/3,1/), domain_size)
    if (fdbk%nc% error /= NF_NOERR) goto 999

    name = 'veri_description'
    vid  = get_varid (fdbk, name)
    c64  = description
    fdbk%nc% error = nf_put_vara_text (fdbk% nc% ncid, vid, (/1,n/), (/64,1/), c64)
    if (fdbk%nc% error /= NF_NOERR) goto 999

    name = 'veri_ens_member'
    vid  = get_varid (fdbk, name)
    fdbk%nc% error = nf_put_vara_int  (fdbk% nc% ncid, vid, n, 1, ens_member)
    if (fdbk%nc% error /= NF_NOERR) goto 999

    name = 'veri_exp_id'
    vid = get_varid (fdbk, name)
    fdbk%nc% error = nf_put_vara_int  (fdbk% nc% ncid, vid, n, 1, exp_id)
    if (fdbk%nc% error /= NF_NOERR) goto 999

    op_flag = 0; if (present (operator_flag)) op_flag = operator_flag 
    name    = 'veri_operator_flag'
    vid     =  get_varid (fdbk, name)
    fdbk%nc% error = nf_put_vara_int  (fdbk% nc% ncid, vid, n, 1, op_flag)
    if (fdbk%nc% error /= NF_NOERR) goto 999

    id_veri = get_varid (fdbk, 'veri_data')

    !-------
    ! return
    !-------
    if (present(ierr)) ierr = 0
    return


999 continue
    !---------------
    ! error handling
    !---------------
    write(6,*) 'add_verification: ERROR writing ',name
    write(6,*) trim(nf_strerror(fdbk%nc% error))
    if (present(ierr)) then
      ierr = fdbk%nc% error
    else
      call model_abort (-1,-1,'ERROR writing '//trim(name)//         &
                              ': '//trim(nf_strerror(fdbk%nc%error)),&
                              'add_verification'                     )
    endif

  end subroutine add_verification
!------------------------------------------------------------------------------
  function get_varid (fb, name) result (varid)
  !-----------------
  ! get NetCDF varid
  !-----------------
  type(t_fdbk)     ,intent(in) :: fb    ! feedback file data type
  character(len=*) ,intent(in) :: name  ! variable name
  integer                      :: varid ! NetCDF varid returned

    integer :: status
    status = nf_inq_varid (fb% nc% ncid, name, varid)
    if (status /= NF_NOERR) then
      write(6,*) 'get_varid: variable not present:',name
      varid = -999
    endif
  end function get_varid
!------------------------------------------------------------------------------
  function get_fillvalue (fb, name) result (fillvalue)
  !-------------------------------------------
  ! return fillvalue used for a given variable
  !-------------------------------------------
  type(t_fdbk)     ,intent(in) :: fb        ! feedback file data type
  character(len=*) ,intent(in) :: name      ! variable name
  real(sp)                     :: fillvalue ! fillvalue

    integer :: status
    integer :: varid

    status = nf_inq_varid (fb% nc% ncid, name, varid)
    if (status /= NF_NOERR) then
      write(6,*) 'get_varid: variable not present:',name
      fillvalue = -999._sp
    endif
    status = nf_get_att_real (fb% nc% ncid, varid, '_FillValue', fillvalue)

  end function get_fillvalue
!------------------------------------------------------------------------------
  subroutine open_fdbk_write (fb, path)
  !--------------------------------------
  ! open a feedback file for write access
  !--------------------------------------
  type(t_fdbk)     ,intent(inout) :: fb    ! feedback file data type
  character(len=*) ,intent(in)    :: path  ! pathname

    fb% nc% path    = path

    call open_netcdf_file_write (fb% nc)

  end subroutine open_fdbk_write
!------------------------------------------------------------------------------
  subroutine open_fdbk_read (fb, path)
  !-------------------------------------
  ! open a feedback file for read access
  !-------------------------------------
  type(t_fdbk)     ,intent(inout) :: fb    ! feedback file data type
  character(len=*) ,intent(in)    :: path  ! pathname

    fb% nc% path    = path
    call open_netcdf_file_read (fb% nc)

  end subroutine open_fdbk_read
!------------------------------------------------------------------------------
  subroutine close_fdbk (fb)
  type (t_fdbk) ,intent(inout) :: fb

    call close_netcdf_file (fb% nc)

  end subroutine close_fdbk
!------------------------------------------------------------------------------
  subroutine cleanup_fdbk (fb)
  type (t_fdbk) ,intent(inout) :: fb
    type (t_fdbk) :: empty
    call destruct_netcdf_file                (fb% nc)
    if (associated (fb% veri))    deallocate (fb% veri)
    if (associated (fb% history)) deallocate (fb% history)
    fb = empty
  end subroutine cleanup_fdbk
!------------------------------------------------------------------------------
  subroutine get_veri_index (idx, nidx, fb, model, run_type, run_class,     &
                             initial_date, forecast_time, ens_member, exp_id,&
                             description)
  !--------------------------------------------------------------
  ! Returns the indices of verification runs in array IDX(:) .
  ! Indices are filtered as specified by the optional parameters.
  ! The number of returned indices is given in NIDX.
  ! If the number of indices requested is larger than size(idx)
  !    NIDX is returned with a negative sign.
  !--------------------------------------------------------------
  integer          ,intent(out)          :: idx (:)  ! indices returned
  integer          ,intent(out)          :: nidx     ! number of indices
  type (t_fdbk)    ,intent(in)           :: fb
  character(len=*) ,intent(in) ,optional :: model
  integer          ,intent(in) ,optional :: run_type
  integer          ,intent(in) ,optional :: run_class
  character(len=*) ,intent(in) ,optional :: initial_date
  integer          ,intent(in) ,optional :: forecast_time
  integer          ,intent(in) ,optional :: ens_member
  integer          ,intent(in) ,optional :: exp_id
  character(len=*) ,intent(in) ,optional :: description
    integer :: n, i
    integer :: ix (fb% n_veri)
    n = 0
    do i = 1, min (fb% n_veri, size (fb% veri))
      if (present (model)) then
        if (model /= fb% veri(i)% model) cycle
      endif
      if (present (run_type)) then
        if (run_type /= fb% veri(i)% run_type) cycle
      endif
      if (present (run_class)) then
        if (run_class /= fb% veri(i)% run_class) cycle
      endif
      if (present (initial_date)) then
        if (initial_date /= fb% veri(i)% initial_date) cycle
      endif
      if (present (forecast_time)) then
        if (forecast_time /= fb% veri(i)% forecast_time) cycle
      endif
      if (present (ens_member)) then
        if (ens_member == VE_MEMBER) then
          if (         0 >= fb% veri(i)% ens_member) cycle
        else
          if (ens_member /= fb% veri(i)% ens_member) cycle
        endif
      endif
      if (present (exp_id)) then
        if (exp_id /= fb% veri(i)% exp_id) cycle
      endif
      if (present (description)) then
        if (trim(description) /= trim(fb% veri(i)% description)) cycle
      endif
      n = n + 1
      ix(n) = i
    end do
    idx = 0
    if (n <= size (idx)) then
      nidx     = n
      idx(1:n) = ix(1:n)
    else
      nidx     = -n
    endif
  end subroutine get_veri_index
!------------------------------------------------------------------------------
  function get_veri (fb, index) result (veri)
  type (t_fdbk)    ,intent(in) :: fb               ! file meta data
  integer          ,intent(in) :: index            ! verification index
  real(sp)                     :: veri(fb% n_body) ! output
  !-------------------------------------------
  ! read verification entry from feedback file
  !-------------------------------------------
    integer :: status, varid
    status = nf_inq_varid (fb% nc% ncid, 'veri_data', varid)
    if (status/=NF_NOERR) then
      write(6,*) 'ERROR: cannot read "veri_data"'
      call model_abort (-1,-1,'cannot read "veri_data"','get_veri')
    endif
    status = nf_get_vara_real(fb% nc% ncid, varid,               &
                              (/1,index/), (/fb% n_body,1/), veri)
    if (status/=NF_NOERR) then
      write(6,*) 'ERROR: cannot read "veri_data"'
      call model_abort (-1,-1,'cannot read "veri_data"','get_veri')
    endif
  end function get_veri
!------------------------------------------------------------------------------
  function get_veri_n (fb, index, n) result (veri)
  type (t_fdbk)    ,intent(in) :: fb               ! file meta data
  integer          ,intent(in) :: index            ! verification index
  integer                      :: n                ! # of entries to read
  real(sp)                     :: veri(n)          ! output
  !-------------------------------------------
  ! read verification entry from feedback file
  !-------------------------------------------
    integer :: status, varid
    status = nf_inq_varid (fb% nc% ncid, 'veri_data', varid)
    if (status/=NF_NOERR) then
      write(6,*) 'ERROR: cannot read "veri_data"'
      call model_abort (-1,-1,'cannot read "veri_data"','get_veri')
    endif
    status = nf_get_vara_real (fb% nc% ncid, varid,           &
                              (/1,index/), (/n,1/), veri)
    if (status/=NF_NOERR) then
      write(6,*) 'ERROR: cannot read "veri_data"'
      call model_abort (-1,-1,'cannot read "veri_data"','get_veri')
    endif
  end function get_veri_n
!------------------------------------------------------------------------------
  subroutine write_veri_dp (fb, index, veri)
  type (t_fdbk)    ,intent(in) :: fb
  integer          ,intent(in) :: index
  real(dp)                     :: veri(:)
  !-------------------------
  ! write verification entry
  !-------------------------
    integer  :: status, varid
    real(sp) :: tmp (size(veri))
    status = nf_inq_varid (fb% nc% ncid, 'veri_data', varid)
    if (status/=NF_NOERR) then
      write(6,*) 'ERROR: cannot write "veri_data"'
      call model_abort (-1,-1,'cannot write "veri_data"','write_veri')
    endif
    tmp = REAL(veri,sp)
    status = nf_put_vara_real(fb% nc% ncid, varid,               &
                              (/1,index/), (/fb% n_body,1/), tmp)
    if (status/=NF_NOERR) then
      write(6,*) 'ERROR: cannot write "veri_data"'
      call model_abort (-1,-1,'cannot write "veri_data"','write_veri')
    endif
  end subroutine write_veri_dp
!------------------------------------------------------------------------------
  subroutine write_veri_sp (fb, index, veri)
  type (t_fdbk)    ,intent(in) :: fb
  integer          ,intent(in) :: index
  real(sp)                     :: veri(:)
  !-------------------------
  ! write verification entry
  !-------------------------
    integer  :: status, varid
    real(sp) :: tmp (size(veri))
    status = nf_inq_varid (fb% nc% ncid, 'veri_data', varid)
    if (status/=NF_NOERR) then
      write(6,*) 'ERROR: cannot write "veri_data"'
      call model_abort (-1,-1,'cannot write "veri_data"','write_veri')
    endif
    tmp = veri
    status = nf_put_vara_real(fb% nc% ncid, varid,              &
                              (/1,index/), (/fb% n_body,1/), tmp)
    if (status/=NF_NOERR) then
      write(6,*) 'ERROR: cannot write "veri_data"'
      call model_abort (-1,-1,'cannot write "veri_data"','write_veri')
    endif
  end subroutine write_veri_sp
!------------------------------------------------------------------------------
  subroutine write_fdbk_int (fb, iv, var, fill)
  !----------------------------------------
  ! write integer variable to feedback file
  !----------------------------------------
  type(t_fdbk)      ,intent(in) :: fb     ! feedback file meta data
  integer           ,intent(in) :: iv     ! variable index
  integer           ,intent(in) :: var(:) ! variable to write
  integer ,optional ,intent(in) :: fill   ! fillvalue to replace

    integer :: tmp (size(var))
    integer :: status

    tmp = var
    if (present (fill)) then
      where (tmp == fill) tmp = fb% nc% vars(iv)% invalid
    endif

    status = nf_put_var_int (fb% nc% ncid, fb% nc% vars(iv)% varid, tmp)
    if (status/=NF_NOERR) then
      write(6,*) 'ERROR: cannot write',fb% nc% vars(iv)% name
      call model_abort (-1,-1,'cannot write '//fb% nc% vars(iv)% name,&
                              'write_fdbk_int'                         )
    endif

  end subroutine write_fdbk_int
!------------------------------------------------------------------------------
  subroutine write_fdbk_real (fb, iv, var, fill)
  !-------------------------------------
  ! write real variable to feedback file
  !-------------------------------------
  type(t_fdbk)       ,intent(in) :: fb     ! feedback file meta data
  integer            ,intent(in) :: iv     ! variable index
  real(wp)           ,intent(in) :: var(:) ! variable to write
  real(wp) ,optional ,intent(in) :: fill   ! fillvalue to replace

    real(sp) :: tmp (size(var))
    integer  :: status

    !-----------------------------------------
    ! copy to default real, replace fillvalues
    !-----------------------------------------
    if (present (fill)) then
      where (var == fill)
        tmp = REAL(fb% nc% vars(iv)% rinvalid,sp)
      elsewhere
        tmp = REAL(var,sp)
      endwhere
    else
      tmp = REAL(var,sp)
    endif
    !---------------
    ! write variable
    !---------------
    status = nf_put_var_real(fb% nc% ncid, fb% nc% vars(iv)% varid, tmp)
    if (status/=NF_NOERR) then
      write(6,*) 'ERROR: cannot write "data"'
      call model_abort (-1,-1,'cannot write "data"','write_fdbk_real')
    endif

  end subroutine write_fdbk_real
!------------------------------------------------------------------------------
  subroutine write_fdbk_char (fb, iv, var)
  !--------------------------------------------------
  ! write character string  variable to feedback file
  !--------------------------------------------------
  type(t_fdbk)       ,intent(in) :: fb     ! feedback file meta data
  integer            ,intent(in) :: iv     ! variable index
  character(len=*)   ,intent(in) :: var(:) ! variable to write

    integer :: status

    status = nf_put_var_text(fb% nc% ncid, fb% nc% vars(iv)% varid, var)

    if (status/=NF_NOERR) then
      write(6,*) 'ERROR: cannot write "data"'
      call model_abort (-1,-1,'cannot write "data"','write_fdbk_char')
    endif

  end subroutine write_fdbk_char
!------------------------------------------------------------------------------
  elemental subroutine blank_nul_chars (s)
    character(len=*), intent(inout) :: s

    integer              :: i
    character, parameter :: NUL = achar (0)

    do i=1, len (s)
       if (s(i:i) == NUL) s(i:i) = ' '
    end do
  end subroutine blank_nul_chars
!------------------------------------------------------------------------------
end module mo_fdbk
