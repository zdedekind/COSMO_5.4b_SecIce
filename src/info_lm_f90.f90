!+ Source module for storing and printing library information
!
MODULE info_lm_f90
!
! Description:
!   This module stores some library information and information of the rules
!   how the binaries have been built.
!   Additionally it provides a subroutine for printing some of this information.
!
! Routines (module procedures) currently contained:
!   - info_print:
!     Print some or all of the information to stdout.
!   - info_readnl:
!     Read some information from the named NAMELIST file.
!   - info_define:
!     Define some information once via the parameter list.
!   - info_print_recursive:
!     Whenever INFO_RECURSIVE is enabled, envelope calls to print the
!     information of subsequent loaded libraries.
!
! Current Code Owner: DWD, Michael Gertz
!  phone: +49 69-8062-2735
!  fax  : +49 69-8062-3721
!  email: Michael.Gertz@dwd.de
!
! History:
! Version       Date       Name
! ------------- ---------- ----
!  V4_13        2010/05/11 Michael Gertz
!  Initial release
! V4_15        2010/11/19 Ulrich Schaettler
!  Replaced LEN_TRIM construct by TRIM
!
! Code Description:
! Language: Fortran 90.
!==============================================================================
!
! Declarations:
!
IMPLICIT NONE
PRIVATE
PUBLIC :: info_readnl, info_define, info_print, info_getvalue
!
!==============================================================================
!
! Global (i.e. public) Declarations:
!
! Global Parameters:
!
! Unfortunately following information is not statically available.
! Besides most of it may change on every checkout or not available at all.
! Therefore the definitions within the repository have to be defined
! empty. They have to be (mostly manually) filled in after the checkout
! of the code or directly before the start of the compilation / linking
! of the binary.
!
! Following declarations have to be defined just after checkout from the version control system:
CHARACTER (LEN=*), PARAMETER :: INFO_LibraryName = 'cosmo5.4b.1'
CHARACTER (LEN=*), PARAMETER :: INFO_RevisionTag = ''
CHARACTER (LEN=*), PARAMETER :: INFO_CheckinDate = '2020-01-13 14:53:01 +0100'
CHARACTER (LEN=*), PARAMETER :: INFO_RevisionNumber = '80016f5 @ (HEAD -> master)'
CHARACTER (LEN=*), PARAMETER :: INFO_CheckoutDate = '(missing)'
CHARACTER (LEN=*), PARAMETER :: INFO_ProductionDate = '(missing)'
! Following declarations have to be defined just before compiling:
CHARACTER (LEN=*), PARAMETER :: INFO_CodeIsModified = 'true'
CHARACTER (LEN=*), PARAMETER :: INFO_CompilerCall = 'ftn -D__PGI_FORTRAN_ -Mpreprocess -DNUDGING -DGRIBDWD -DGRIBAPI ' // &
    '-DNETCDF -DTWOMOM_SB -DFOR_LM -Kieee -Mfree -Mdclchk -DHAS_IOMSG -D__COSMO__ -O3 -Mvect=sse -Mlre ' // &
    '-Mvect=noassoc -Mipa=fast,inline -Mnofma'
CHARACTER (LEN=*), PARAMETER :: INFO_CompilerVersion = 'pgf90 19.7-0 LLVM 64-bit target on x86-64 Linux ' // &
    '-tp haswell-64'
CHARACTER (LEN=*), PARAMETER :: INFO_DefinedMacros = '-I. ' // &
    '-I/users/dedekind/code/COSMO5.4b.1_2D_SecIce_Rates/src ' // &
    '-I/opt/cray/pe/mpt/7.7.10/gni/mpich-pgi/19.1/include -I/project/s799/geirund/software/lib ' // &
    '-D__MPICH2 -I/project/c14/install/daint/libgrib_api/v1.20.0.3/pgi/include'
CHARACTER (LEN=*), PARAMETER :: INFO_UndefinedMacros = '(missing)'
CHARACTER (LEN=*), PARAMETER :: INFO_DebugOptions = '(missing)'
CHARACTER (LEN=*), PARAMETER :: INFO_LinkOptions = 'ftn -D__PGI_FORTRAN_ = -L/project/s799/geirund/software/lib ' // &
    '-L/project/c14/install/daint/libgrib1 -lgrib1_pgi ' // &
    '-L/project/c14/install/daint/libgrib_api/v1.20.0.3/pgi/lib -lgrib_api_f90 -lgrib_api ' // &
    '-L/project/c14/install/daint/libjasper/lib -ljasper'
CHARACTER (LEN=*), PARAMETER :: INFO_CompiledBy = 'dedekind'
CHARACTER (LEN=*), PARAMETER :: INFO_CompileTime = 'Mon Jan 13 17:25:00 CET 2020'
CHARACTER (LEN=*), PARAMETER :: INFO_CompileMachine = 'daint'
CHARACTER (LEN=*), PARAMETER :: INFO_GribApiVersion  = ''
!
! Global Variables:
!
! Following information has to be evaluated during runtime. Therefore it will
! be available after the first call of info_print(), info_readnl() or info_define().
CHARACTER (LEN=16) :: INFO_StartTime = ''
!
! Currently there is no way to fill in the information with std-routines.
! Therefore they have to be filled by the user on calling info_define() or info_readnl().
CHARACTER (LEN=80) :: INFO_BinaryName = 'cosmo'
CHARACTER (LEN=80) :: INFO_RunMachine = ''
CHARACTER (LEN=80) :: INFO_Nodes = ''
CHARACTER (LEN=80) :: INFO_Domain = ''
!
! Define local default options on demand. This will be done by the user on
! calling info_readnl().
CHARACTER (LEN=80) :: INFO_Options = ''
!
!- End of module header

CONTAINS
!+
SUBROUTINE info_readnl ( NamelistFile )
!
! Description:
!   Read some information from the namelist input file.
!
!   Currently it is not possible with FORTRAN95 to get the information
!   of the full path of binary name like the $0 in C. Additionally
!   we cannot determine on which host(s) the binary is running and the
!   domain of the data spread through the nodes.
!   Therefore this information has to be defined manually. On using info_readnl()
!   this information may be defined within the segment /info_defaults/ which
!   has to reside within the named namelist of your choice. Missing information
!   will be ignored silently.
!   Currently following information may be defined within /info_defaults/:
!   INFO_Options ..: List of print options
!   INFO_BinaryName: Name (best: full path) of the binary
!   INFO_RunMachine: The machine (OS) where the program is running
!   INFO_Nodes ....: Description of the nodes the binary is running
!   INFO_Domain ...: The domain the binary is calculating
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
CHARACTER (LEN=*), INTENT(IN) :: NamelistFile  ! Full path name of the namelist file
!
! Local variables:
CHARACTER (LEN=8)  :: Date             ! Current date
CHARACTER (LEN=10) :: Time             ! Current time
LOGICAL            :: FileExist
INTEGER            :: IOStat
!
NAMELIST /info_defaults/ INFO_Options, INFO_BinaryName, INFO_RunMachine, INFO_Nodes, INFO_Domain
!
!- End of header
!
! Initialize the unknown runtime parameters as far as possible:
!
if ( LEN_TRIM(INFO_StartTime) == 0 ) THEN
  CALL DATE_AND_TIME(Date, Time)
  INFO_StartTime = Date(1:4) // '-' // Date(5:6) // '-' // Date(7:8) // &
                   ' ' // Time(1:2) // ':' // Time(3:4) // ':' // Time(5:6)
END IF
!
! Reading in the variables from the namelist file whenever available:
!
INQUIRE (FILE=NamelistFile, EXIST=FileExist, IOSTAT=IOStat)
IF ( FileExist ) THEN
  OPEN (UNIT=35, FILE=NamelistFile, STATUS='OLD', ACTION='READ', IOSTAT=IOStat)
  IF ( IOStat == 0 ) THEN
    READ (UNIT=35, NML=info_defaults, IOSTAT=IOStat)
    CLOSE (UNIT=35, IOSTAT=IOStat)
  END IF
END IF
!
RETURN
END SUBROUTINE info_readnl
!+
SUBROUTINE info_define ( BinaryName, RunMachine, Nodes, Domain )
!
! Description:
!   Define some information for later usage.
!
!   Currently it is not possible with FORTRAN95 to get the information
!   of the full path of binary name like the $0 in C. Additionally
!   we cannot determine on which host(s) the binary is running and the
!   domain of the data spread through the nodes.
!   Maybe the user is able with proprietary routines to fetch this
!   information and store ist into the common variables on using this routine.
!   Currently following information may be defined:
!   INFO_BinaryName: Name (best: full path) of the binary
!   INFO_RunMachine: The machine (OS) where the program is running
!   INFO_Nodes ....: Description of the nodes the binary is running
!   INFO_Domain ...: The domain the binary is calculating
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
CHARACTER (LEN=*), OPTIONAL, INTENT(IN) :: BinaryName ! Full path name of the binary
CHARACTER (LEN=*), OPTIONAL, INTENT(IN) :: RunMachine ! Machine (OS) where the binary is running
CHARACTER (LEN=*), OPTIONAL, INTENT(IN) :: Nodes      ! List of nodes where the binary is running
CHARACTER (LEN=*), OPTIONAL, INTENT(IN) :: Domain     ! Data domain
!
! Local variables:
CHARACTER (LEN=8)  :: Date             ! Current date
CHARACTER (LEN=10) :: Time             ! Current time
!
!- End of header
!
! Initialize the unknown runtime parameters as far as possible:
!
if ( LEN_TRIM(INFO_StartTime) == 0 ) THEN
  CALL DATE_AND_TIME(Date, Time)
  INFO_StartTime = Date(1:4) // '-' // Date(5:6) // '-' // Date(7:8) // &
                   ' ' // Time(1:2) // ':' // Time(3:4)
END IF
!
! Initialize local variables:
!
IF ( (PRESENT(BinaryName) ) .AND. (LEN_TRIM(INFO_BinaryName) == 0) ) THEN
  INFO_BinaryName = BinaryName(1:MIN(LEN_TRIM(BinaryName),80))
END IF
IF ( (PRESENT(RunMachine) ) .AND. (LEN_TRIM(INFO_RunMachine) == 0) ) THEN
  INFO_RunMachine = RunMachine(1:MIN(LEN_TRIM(RunMachine),80))
END IF
IF ( (PRESENT(Nodes) ) .AND. (LEN_TRIM(INFO_Nodes) == 0) ) THEN
  INFO_Nodes = Nodes(1:MIN(LEN_TRIM(Nodes),80))
END IF
IF ( (PRESENT(Domain) ) .AND. (LEN_TRIM(INFO_Domain) == 0) ) THEN
  INFO_Domain = Domain(1:MIN(LEN_TRIM(Domain),80))
END IF
!
RETURN
END SUBROUTINE info_define
!+
PURE FUNCTION info_getvalue ( Options )
!
! Description:
!   Get one value of information.
!
!   The character defines the identical value which will be printed by info_print().
!
!   Not every information may be available at any time on any machine. Whenever
!   no information is available it has to default to an empty string.
!
!   Currently following option characters are available:
!    '?' Special: Get a list of all Options available
!    'c' Get the time the code had been compiled
!    'd' Get the debug options used for compilation
!    'i' Get the time the code had been put into the version system
!    'm' Get whether the code has been marked locally modfied
!    'n' Get the library name
!    'o' Get the machine (OS) where the binary has been compiled
!    'p' Get the time the code had been put into production
!    'r' Get the revision number from where the code had been extracted
!    't' Get the revision tag from where the code had been extracted
!    'x' Get the time the code had been extracted from the version system
!    'A' Get the GRIB_API version number
!    'B' Get the full path of the binary name
!    'C' Get the user (login name) who did the compiling
!    'D' Get all the macros defined (-D...) during compilation
!    'L' Get the compiler options used for linking
!    'N' Get the compiler name used for compilation
!    'O' Get the machine (OS) where the binary is running
!    'R' Get the information of other used libraries, too
!    'S' Get the runtime date/time of the strart of the binary
!    'U' Get all the macros undefined (-U...) during compilation
!    'V' Get the compiler version used for compilation
!    'W' Get the nodes where the binary is running
!    'X' Get the data decomposition
!==============================================================================
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
CHARACTER (LEN=*), OPTIONAL, INTENT(IN) :: Options  ! Declare what to print
!
! Local variables:
CHARACTER (LEN=8)  :: Date             ! Current date
CHARACTER (LEN=10) :: Time             ! Current time
CHARACTER (LEN=1)  :: PrintOption
CHARACTER (LEN=80) :: info_getvalue
!
!- End of header
!
! Get the value:
!
IF ( LEN_TRIM(Options) > 0 ) THEN
  PrintOption = Options(1:1)
ELSE
  PrintOption = '?'
END IF
!
SELECT CASE ( PrintOption )
CASE ( 'c' )   ! Compiled timestamp
  info_getvalue = TRIM(INFO_CompileTime)
CASE ( 'd' )   ! Debug options
  info_getvalue = TRIM(INFO_DebugOptions)
CASE ( 'i' )   ! Checkin timestamp
  info_getvalue = TRIM(INFO_CheckinDate)
CASE ( 'm' )   ! Modified flag
  info_getvalue = TRIM(INFO_CodeIsModified)
CASE ( 'n' )   ! Library name
  info_getvalue = TRIM(INFO_LibraryName)
CASE ( 'o' )   ! machine (OS) of compilation
  info_getvalue = TRIM(INFO_CompileMachine)
CASE ( 'p' )   ! Production timestamp
  info_getvalue = TRIM(INFO_ProductionDate)
CASE ( 'r' )   ! Revision number
  info_getvalue = TRIM(INFO_Revisionnumber)
CASE ( 't' )   ! Version tag
  info_getvalue = TRIM(INFO_RevisionTag)
CASE ( 'x' )   ! Checkout timestamp
  info_getvalue = TRIM(INFO_CheckoutDate)
CASE ( 'A' )   ! GRIB_API version number
  info_getvalue = TRIM(INFO_GribApiVersion)
CASE ( 'B' )   ! Binary name (full path)
  info_getvalue = TRIM(INFO_BinaryName)
CASE ( 'C' )   ! Compiling user
  info_getvalue = TRIM(INFO_CompiledBy)
CASE ( 'D' )   ! Defined macros
  info_getvalue = TRIM(INFO_DefinedMacros)
CASE ( 'L' )   ! Linker call
  info_getvalue = TRIM(INFO_LinkOptions)
CASE ( 'N' )   ! Compiler name
  info_getvalue = TRIM(INFO_CompilerCall)
CASE ( 'O' )   ! machine (OS) during runtime
  info_getvalue = TRIM(INFO_RunMachine)
CASE ( 'S' )   ! Start time
  info_getvalue = TRIM(INFO_StartTime)
CASE ( 'U' )   ! Undefined macros
  info_getvalue = TRIM(INFO_UndefinedMacros)
CASE ( 'V' )   ! Compiler version
  info_getvalue = TRIM(INFO_CompilerVersion)
CASE ( 'W' )   ! Where running
  info_getvalue = TRIM(INFO_Nodes)
CASE ( 'X' )   ! Data decomposition
  info_getvalue = TRIM(INFO_Domain)
CASE ( '?' )   ! List of values
  info_getvalue = 'cdimnoprtxABCDLNOSUVWX'
CASE DEFAULT
  info_getvalue = 'Invalid option: ' // Options
END SELECT
!
END FUNCTION info_getvalue
!+
SUBROUTINE info_print ( Options )
!
! Description:
!   Print some of the information to stdout.
!
!   The information will be printet line by line. On naming some options within
!   the Options list you may decide which information should be printed at which
!   line. Every single character within the list defines on line of information.
!   The Options string will be processed from left to right. Empty lines may be
!   created with the blank character.
!
!   Not every information may be available at any time on any machine. Whenever
!   no information is available it has to default to an empty string.
!
!   Currently following option characters are available:
!    '!' Special: Suppress all output (do NOT print any information)
!    '?' Special: Print all information (no additional option allowed)
!    ' ' Print a empty line
!    'c' Print the time the code had been compiled
!    'd' Print the debug options used for compilation
!    'i' Print the time the code had been put into the version system
!    'm' Print whether the code has been marked locally modfied
!    'n' Print the library name
!    'o' Get the machine (OS) where the binary has been compiled
!    'p' Print the time the code had been put into production
!    'r' Print the revision number from where the code had been extracted
!    't' Print the revision tag from where the code had been extracted
!    'x' Print the time the code had been extracted from the version system
!    'A' Print the GRIB_API version number
!    'B' Print the full path of the binary name
!    'C' Print the user (login name) who did the compiling
!    'D' Print all the macros defined (-D...) during compilation
!    'L' Print the compiler options used for linking
!    'N' Print the compiler name used for compilation
!    'O' Get the machine (OS) where the binary is running
!    'R' Print the information of other used libraries, too
!    'S' Print the runtime date/time of the strart of the binary
!    'U' Print all the macros undefined (-U...) during compilation
!    'V' Print the compiler version used for compilation
!    'W' Print the nodes where the binary is running
!    'X' Print the data decomposition
!==============================================================================
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
CHARACTER (LEN=*), OPTIONAL, INTENT(IN) :: Options  ! Declare what to print
!
! Local parameters:
CHARACTER (LEN=*), PARAMETER :: AllOptions      = 'B ntripxmcCo NVDUdLA SOWXR'
CHARACTER (LEN=*), PARAMETER :: EnhancedOptions = 'B ntixmcCo DUdA SOWXR'
CHARACTER (LEN=*), PARAMETER :: DefaultOptions  = 'B ntimcCA SWXR'
!
! Local variables:
INTEGER            :: NumberOfOptions  ! Count of information to print
INTEGER            :: Count            ! Local counter
INTEGER            :: UseDefault       ! Default flag
INTEGER            :: Error            ! Number of errors
INTEGER            :: DoRecursive      ! Run recursively through subsequent libraries
CHARACTER (LEN=8)  :: Date             ! Current date
CHARACTER (LEN=10) :: Time             ! Current time
CHARACTER (LEN=1)  :: PrintOption
!
!- End of header
!
! Initialize local variables:
!
Error           = 0
DoRecursive     = 0
!
! Initialize the unknown runtime parameters as far as possible:
!
if ( LEN_TRIM(INFO_StartTime) == 0 ) THEN
  CALL DATE_AND_TIME(Date, Time)
  INFO_StartTime = Date(1:4) // '-' // Date(5:6) // '-' // Date(7:8) // &
                   ' ' // Time(1:2) // ':' // Time(3:4)
END IF
!
! Decide which Options should be used:
!
IF ( PRESENT(Options) ) THEN
  IF ( LEN_TRIM(Options) > 0 ) THEN
    IF ( Options(1:1) == '?' ) THEN ! All options requested
      UseDefault      = 2
      NumberOfOptions = LEN_TRIM(AllOptions)
    ELSE IF ( Options(1:1) == '!' ) THEN ! Do nothing
      RETURN
    ELSE
      UseDefault      = 0
      NumberOfOptions = LEN_TRIM(Options)
    END IF
  END IF
ELSE
  IF ( LEN_TRIM(INFO_Options) > 0 ) THEN ! User pre-defined options required
    IF ( INFO_Options(1:1) == '?' ) THEN ! All options requested
      UseDefault      = 2
      NumberOfOptions = LEN_TRIM(AllOptions)
    ELSE IF ( INFO_Options(1:1) == '!' ) THEN ! Do nothing
      RETURN
    ELSE
      UseDefault      = 1
      NumberOfOptions = LEN_TRIM(INFO_Options)
    END IF
  ELSE IF ( INDEX(INFO_CodeIsModified, 'true') > 0 ) THEN ! Enhanced options requested
    UseDefault      = 3
    NumberOfOptions = LEN_TRIM(EnhancedOptions)
  ELSE
    UseDefault      = 4
    NumberOfOptions = LEN_TRIM(DefaultOptions)
  END IF
END IF
!
! Print the header line:
!
PRINT *, ''
PRINT *, '==== Code information used to build this binary ===='
!
! Loop through the options:
!
DO Count = 1, NumberOfOptions
  IF ( UseDefault == 0 ) THEN            ! Use named options
    PrintOption = Options(Count:Count)
  ELSE IF ( UseDefault == 1 ) THEN       ! Use predefined options
    PrintOption = INFO_Options(Count:Count)
  ELSE IF ( UseDefault == 2 ) THEN       ! Use all options
    PrintOption = AllOptions(Count:Count)
  ELSE IF ( UseDefault == 3 ) THEN       ! Use optimized options
    PrintOption = EnhancedOptions(Count:Count)
  ELSE                                   ! Use default options
    PrintOption = DefaultOptions(Count:Count)
  END IF
!
  SELECT CASE ( PrintOption )
  CASE ( ' ' )   ! Empty line
    PRINT *, ''
  CASE ( 'c' )   ! Compiled timestamp
    PRINT *, 'Compile-Date ......: ', TRIM(INFO_CompileTime)
  CASE ( 'd' )   ! Debug options
    PRINT *, 'Debug options .....: ', TRIM(INFO_DebugOptions)
  CASE ( 'i' )   ! Checkin timestamp
    PRINT *, 'Checkin-Date ......: ', TRIM(INFO_CheckinDate)
  CASE ( 'm' )   ! Modified flag
    PRINT *, 'Code is modified ..: ', TRIM(INFO_CodeIsModified)
  CASE ( 'n' )   ! Library name
    PRINT *, 'Library name ......: ', TRIM(INFO_LibraryName)
  CASE ( 'o' )   ! machine (OS) of compilation
    PRINT *, 'Compiled on .......: ', TRIM(INFO_CompileMachine)
  CASE ( 'p' )   ! Production timestamp
    PRINT *, 'Put into production: ', TRIM(INFO_ProductionDate)
  CASE ( 'r' )   ! Revision number
    PRINT *, 'Revision number ...: ', TRIM(INFO_Revisionnumber)
  CASE ( 't' )   ! Version tag
    PRINT *, 'Tag name ..........: ', TRIM(INFO_RevisionTag)
  CASE ( 'x' )   ! Checkout timestamp
    PRINT *, 'Checkout-Date .....: ', TRIM(INFO_CheckoutDate)
  CASE ( 'A' )   ! GRIB_API version number
    PRINT *, 'GRIB_API version ..: ', TRIM(INFO_GribApiVersion)
  CASE ( 'B' )   ! Binary name (full path)
    PRINT *, 'Binary name .......: ', TRIM(INFO_BinaryName)
  CASE ( 'C' )   ! Compiling user
    PRINT *, 'Compiled by .......: ', TRIM(INFO_CompiledBy)
  CASE ( 'D' )   ! Defined macros
    PRINT *, 'Macros defined ....: ', TRIM(INFO_DefinedMacros)
  CASE ( 'L' )   ! Linker call
    PRINT *, 'Linker options ....: ', TRIM(INFO_LinkOptions)
  CASE ( 'N' )   ! Compiler name
    PRINT *, 'Compiler call .....: ', TRIM(INFO_CompilerCall)
  CASE ( 'O' )   ! machine (OS) during runtime
    PRINT *, 'Running on machine : ', TRIM(INFO_RunMachine)
  CASE ( 'R' )   ! Info of other used libraries
    DoRecursive = 1     ! Activate info_call of other libraries
  CASE ( 'S' )   ! Start time
    PRINT *, 'Current start time : ', TRIM(INFO_StartTime)
  CASE ( 'U' )   ! Undefined macros
    PRINT *, 'Macros undefined ..: ', TRIM(INFO_UndefinedMacros)
  CASE ( 'V' )   ! Compiler version
    PRINT *, 'Compiler version ..: ', TRIM(INFO_CompilerVersion)
  CASE ( 'W' )   ! Where running
    PRINT *, 'Running on nodes ..: ', TRIM(INFO_Nodes)
  CASE ( 'X' )   ! Data decomposition
    PRINT *, 'Data decomposition : ', TRIM(INFO_Domain)
  CASE DEFAULT
    Error = Error + 1    ! Set error flag
  END SELECT
END DO
!
! Print the trailer line:
!
PRINT *, '==== End of code information ===='
PRINT *, ''
!
! Whenever the error flag has been set, print help
!
IF ( Error > 0 ) THEN
  PRINT *, 'Warning: Invalid option within option string: ', Options
  PRINT *, Error, 'defined options are illegal.'
  PRINT *, ''
  PRINT *, 'Following print options are available:'
  PRINT *, '!: Special: Suppress all output (do NOT print any information)'
  PRINT *, '?: Special: Print all information (no additional option allowed)'
  PRINT *, ' : Print a empty line'
  PRINT *, 'c: Print the time the code had been compiled'
  PRINT *, 'd: Print the debug options used for compilation'
  PRINT *, 'i: Print the time when the code had been put into the version system'
  PRINT *, 'm: Print whether the code has been marked modified'
  PRINT *, 'n: Print the library name'
  PRINT *, 'o: Print the machine (OS) where the binary has been compiled'
  PRINT *, 'p: Print the time the code had been put into production'
  PRINT *, 'r: Print the revision number'
  PRINT *, 't: Print the revision tag'
  PRINT *, 'x: Print the extraction time'
  PRINT *, 'A: Print the GRIB_API version number'
  PRINT *, 'B: Print the binary name (full path)'
  PRINT *, 'C: Print the user (login name) who did the compiling'
  PRINT *, 'D: Print all the macros defined (-D...) during compilation'
  PRINT *, 'L: Print the compiler options used for linking'
  PRINT *, 'N: Print the compiler name used for compilation'
  PRINT *, 'O: Print the machine (OS) where the binary is running'
  PRINT *, 'R: Print the information of other used libraries, too'
  PRINT *, 'S: Print the start date/time of the binary'
  PRINT *, 'U: Print all the macros undefined (-U...) during compilation'
  PRINT *, 'V: Print the compiler version used for compilation'
  PRINT *, 'W: Print the nodes where the program is running'
  PRINT *, 'X: Print the data decomposition'
  PRINT *, ''
END IF
!
! Check whether we shall do recursion to other libraries:
!
IF ( DoRecursive > 0 ) THEN
  CALL info_print_recursive ( Options, UseDefault )
END IF
!
RETURN
END SUBROUTINE info_print
!+
! Providing envelope subroutines for additional used libraries:
!
SUBROUTINE info_print_recursive ( Options, UseDefault )
!
! Description:
!   Print some of the information of the library io to stdout.
!
!==============================================================================
! Declarations:
!
#ifdef INFO_RECURSIVE
! Load the info_print functions of the info_* modules of all subequent libraries:
! Define them with an alias to get unique names!
!
!Example USE info_io,       ONLY: info_print_io       => info_print
!Example USE info_math_dwd, ONLY: info_print_math_dwd => info_print
#endif
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER, OPTIONAL, INTENT(IN) :: UseDefault
CHARACTER (LEN=*), INTENT(IN) :: Options  ! Declare what to print
!
#ifdef INFO_RECURSIVE
! Local parameters:
! These parameters will be slightly different to the above ones, because
! we assume to use a production library without any changes.
! Therefore some of the variables are very useless for us.
CHARACTER (LEN=*), PARAMETER :: AllOptions      = 'ntricCo NVDUd'
CHARACTER (LEN=*), PARAMETER :: EnhancedOptions = 'ntic DUd'
CHARACTER (LEN=*), PARAMETER :: DefaultOptions  = 'ntic'
!
!- End of header
!
! Print the information on calling the print_info routine of every subsequent library:
!
IF ( PRESENT(UseDefault) ) THEN
  IF ( UseDefault == 0 ) THEN            ! Use named options
!Example     CALL info_print_io ( Options )
!Example     CALL info_print_math_dwd ( Options )
  ELSE IF ( UseDefault == 1 ) THEN       ! Use all options
!Example     CALL info_print_io ( AllOptions )
!Example     CALL info_print_math_dwd ( AllOptions )
  ELSE IF ( UseDefault == 2 ) THEN       ! Use optimized options
!Example     CALL info_print_io ( EnhancedOptions )
!Example     CALL info_print_math_dwd ( EnhancedOptions )
  ELSE                                   ! Use default options
!Example     CALL info_print_io ( DefaultOptions )
!Example     CALL info_print_math_dwd ( DefaultOptions )
  END IF
ELSE
!Example   CALL info_print_io ( Options )
!Example   CALL info_print_math_dwd ( Options )
END IF
!
#endif
RETURN
END SUBROUTINE info_print_recursive
!
END MODULE info_lm_f90
